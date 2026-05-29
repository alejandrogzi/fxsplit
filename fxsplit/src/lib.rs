// Copyright (c) 2026 Alejandro Gonzalez-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

#![allow(clippy::doc_overindented_list_items)]

use std::{
    collections::{HashMap, HashSet},
    fs::{create_dir_all, File},
    io::{self, BufWriter, Read, Write},
    ops::Range,
    os::unix::fs::symlink,
    path::{Path, PathBuf},
    sync::Arc,
};

use anyhow::{anyhow, Context, Result};
use bzip2::write::BzEncoder;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression as GzCompression;
use memchr::{memchr, memchr_iter};
use memmap2::Mmap;
use rayon::prelude::*;
use twobit::{nucleotide::Nucleotide, TwoBitFile};

pub mod cli;
use cli::Args;

const FASTA_HEADER: u8 = b'>';
const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];
const TWOBIT_MAGIC: [u8; 4] = [0x1A, 0x41, 0x27, 0x43];
const TWOBIT_MAGIC_REVERSED: [u8; 4] = [0x43, 0x27, 0x41, 0x1A];
const HEADER_SPLIT_WARN_THRESHOLD: usize = 100_000;
const TWOBIT_SIGNATURE: u32 = 0x1A41_2743;
/// Line width used when wrapping sequences emitted as FASTA from 2bit input.
const FASTA_LINE_WIDTH: usize = 60;

/// Detects the input format based on file content.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum InputFormat {
    /// Plain FASTA format (.fa, .fasta)
    Fasta,
    /// Gzipped FASTA format (.fa.gz, .fasta.gz)
    FastaGz,
    /// Plain FASTQ format (.fq, .fastq)
    Fastq,
    /// Gzipped FASTQ format (.fq.gz, .fastq.gz)
    FastqGz,
    /// UCSC 2bit binary format (.2bit)
    TwoBit,
}

/// Represents resolved input sources for splitting operations.
///
/// Distinguishes between file-based inputs (via mmap) and
/// memory-based inputs (from stdin or when passthrough optimization is used).
enum ResolvedInput {
    /// Plain FASTA file path
    FastaPath(PathBuf),
    /// Gzipped FASTA file path
    FastaGzPath(PathBuf),
    /// Plain FASTQ file path
    FastqPath(PathBuf),
    /// Gzipped FASTQ file path
    FastqGzPath(PathBuf),
    /// 2bit file path
    TwoBitPath(PathBuf),
    /// FASTA data in memory
    FastaMemory {
        /// The FASTA data bytes
        data: SharedBytes<Vec<u8>>,
    },
    /// FASTQ data in memory
    FastqMemory {
        /// The FASTQ data bytes
        data: SharedBytes<Vec<u8>>,
    },
    /// 2bit data in memory
    TwoBitMemory {
        /// The 2bit data bytes
        data: SharedBytes<Vec<u8>>,
    },
}

/// Output sequence format selected via `--output-format`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub enum OutputFormat {
    /// Plain or compressed FASTA
    Fasta,
    /// Plain or compressed FASTQ (only valid when the input is FASTQ)
    Fastq,
    /// UCSC 2bit binary (never compressed)
    #[value(name = "2bit")]
    TwoBit,
}

impl OutputFormat {
    /// Human-readable label used in diagnostics.
    fn label(self) -> &'static str {
        match self {
            Self::Fasta => "FASTA",
            Self::Fastq => "FASTQ",
            Self::TwoBit => "2bit",
        }
    }
}

/// Output compression codec selected via `--compression`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub enum Compression {
    /// gzip (DEFLATE) compression
    Gzip,
    /// bzip2 compression
    Bzip2,
    /// Zstandard compression
    Zstd,
}

impl Compression {
    /// Returns the filename suffix (including the leading dot) for this codec.
    fn extension_suffix(self) -> &'static str {
        match self {
            Self::Gzip => ".gz",
            Self::Bzip2 => ".bz2",
            Self::Zstd => ".zst",
        }
    }
}

/// Distinguishes the source byte layout passed to `write_chunk`.
#[derive(Debug, Clone, Copy)]
enum SourceKind {
    /// FASTA record bytes
    Fasta,
    /// FASTQ record bytes
    Fastq,
}

/// Defines how records should be distributed across output files.
#[derive(Debug, Clone, Copy)]
pub enum SplitMode {
    /// Split by fixed number of records per file
    ChunkSize(usize),
    /// Split into a fixed number of files
    NumFiles(usize),
    /// Create one output file per FASTA/2bit header
    FileHeader,
}

/// Represents a range of record indices for chunking.
///
/// Used to specify which records belong to each output chunk.
#[derive(Debug, Clone)]
pub struct ChunkRegion {
    /// Start index (inclusive) of the chunk
    pub start: usize,
    /// End index (exclusive) of the chunk
    pub end: usize,
}

/// Thread-safe wrapper for shared byte data using Arc.
struct SharedBytes<B> {
    inner: Arc<B>,
}

impl<B> SharedBytes<B> {
    /// Creates a new SharedBytes from an Arc.
    ///
    /// # Arguments
    /// * `inner` - The Arc-wrapped data to share
    fn from_arc(inner: Arc<B>) -> Self {
        Self { inner }
    }
}

impl<B> Clone for SharedBytes<B> {
    fn clone(&self) -> Self {
        Self {
            inner: Arc::clone(&self.inner),
        }
    }
}

impl<B: AsRef<[u8]>> SharedBytes<B> {
    fn as_slice(&self) -> &[u8] {
        self.inner.as_ref().as_ref()
    }
}

impl<B: AsRef<[u8]>> AsRef<[u8]> for SharedBytes<B> {
    fn as_ref(&self) -> &[u8] {
        self.as_slice()
    }
}

/// Writer for sequence files supporting plain and compressed output.
enum SequenceWriter {
    /// Plain text writer
    Plain(BufWriter<File>),
    /// Gzipped writer
    Gzip(GzEncoder<BufWriter<File>>),
    /// bzip2 writer
    Bzip2(BzEncoder<BufWriter<File>>),
    /// Zstandard writer
    Zstd(zstd::Encoder<'static, BufWriter<File>>),
}

impl SequenceWriter {
    /// Creates a new SequenceWriter for the given path and compression codec.
    ///
    /// # Arguments
    /// * `path` - Output file path
    /// * `compression` - Compression codec, or `None` for plain output
    ///
    /// # Example
    /// ```rust, ignore
    /// let writer = SequenceWriter::create(Path::new("output.fa"), None)?;
    /// ```
    fn create(path: &Path, compression: Option<Compression>) -> Result<Self> {
        let buffer = BufWriter::new(File::create(path)?);
        Ok(match compression {
            None => Self::Plain(buffer),
            Some(Compression::Gzip) => Self::Gzip(GzEncoder::new(buffer, GzCompression::fast())),
            Some(Compression::Bzip2) => {
                Self::Bzip2(BzEncoder::new(buffer, bzip2::Compression::fast()))
            }
            Some(Compression::Zstd) => Self::Zstd(zstd::Encoder::new(buffer, 3)?),
        })
    }

    /// Flushes and finalizes the writer.
    ///
    /// For compressed writers, this ensures all data is flushed and the
    /// codec trailer is written.
    fn finish(self) -> Result<()> {
        match self {
            Self::Plain(mut writer) => {
                writer.flush()?;
                Ok(())
            }
            Self::Gzip(mut encoder) => {
                encoder.flush()?;
                let mut writer = encoder.finish()?;
                writer.flush()?;
                Ok(())
            }
            Self::Bzip2(mut encoder) => {
                encoder.flush()?;
                let mut writer = encoder.finish()?;
                writer.flush()?;
                Ok(())
            }
            Self::Zstd(encoder) => {
                let mut writer = encoder.finish()?;
                writer.flush()?;
                Ok(())
            }
        }
    }
}

impl Write for SequenceWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            Self::Plain(writer) => writer.write(buf),
            Self::Gzip(writer) => writer.write(buf),
            Self::Bzip2(writer) => writer.write(buf),
            Self::Zstd(writer) => writer.write(buf),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            Self::Plain(writer) => writer.flush(),
            Self::Gzip(writer) => writer.flush(),
            Self::Bzip2(writer) => writer.flush(),
            Self::Zstd(writer) => writer.flush(),
        }
    }
}

/// Represents a sequence from a 2bit file.
#[derive(Debug, Clone)]
struct TwoBitSequence {
    /// Sequence name/chromosome identifier
    name: String,
    /// Sequence length in bases
    length: usize,
    /// Sequence bases as bytes (A=0, C=1, G=2, T=3)
    bases: Vec<u8>,
    /// Hard-masked regions (N)
    hard_blocks: Vec<Range<usize>>,
    /// Soft-masked regions (lowercase)
    soft_blocks: Vec<Range<usize>>,
}

/// Main entry point for splitting operations.
///
/// Automatically detects input format and routes to appropriate splitter.
///
/// # Arguments
/// * `args` - Parsed command-line arguments
///
/// # Example
/// ```rust, ignore
/// let args = Args::parse();
/// run(&args)?;
/// ```
pub fn run(args: &Args) -> Result<()> {
    match resolve_input(args)? {
        ResolvedInput::FastaPath(path) => split_fasta_path(args, &path),
        ResolvedInput::FastaGzPath(path) => split_fasta_gz_path(args, &path),
        ResolvedInput::FastqPath(path) => split_fastq_path(args, &path),
        ResolvedInput::FastqGzPath(path) => split_fastq_gz_path(args, &path),
        ResolvedInput::TwoBitPath(path) => split_2bit_path(args, &path),
        ResolvedInput::FastaMemory { data } => split_fasta_bytes(args, data, None),
        ResolvedInput::FastqMemory { data } => split_fastq_bytes(args, data),
        ResolvedInput::TwoBitMemory { data } => split_2bit_bytes(args, data),
    }
}

/// Library entry point for isolated splitting without CLI parsing.
///
/// Useful for programmatic use cases where Args is not needed.
///
/// # Arguments
/// * `args` - Command-line arguments as strings (excludes program name)
///
/// # Example
/// ```rust, ignore
/// fxsplit::lib_iso_split(vec!["-f".into(), "input.fa".into(), "-c".into(), "100".into()])?;
/// ```
pub fn lib_iso_split(args: Vec<String>) -> Result<()> {
    let args = cli::Args::from(args);
    run(&args)
}

/// Splits plain FASTA files only.
///
/// # Arguments
/// * `args` - Parsed command-line arguments
///
/// # Errors
/// Returns error if input is not plain FASTA
pub fn split_fa(args: &Args) -> Result<()> {
    match resolve_input(args)? {
        ResolvedInput::FastaPath(path) => split_fasta_path(args, &path),
        ResolvedInput::FastaMemory { data } => split_fasta_bytes(args, data, None),
        _ => anyhow::bail!("ERROR: input is not plain FASTA"),
    }
}

/// Splits gzipped FASTA files only.
///
/// # Arguments
/// * `args` - Parsed command-line arguments
///
/// # Errors
/// Returns error if input is not gzipped FASTA
pub fn split_fa_gz(args: &Args) -> Result<()> {
    match resolve_input(args)? {
        ResolvedInput::FastaGzPath(path) => split_fasta_gz_path(args, &path),
        _ => anyhow::bail!("ERROR: input is not gzipped FASTA"),
    }
}

/// Splits FASTQ files (plain or gzipped).
///
/// # Arguments
/// * `args` - Parsed command-line arguments
///
/// # Errors
/// Returns error if input is not FASTQ format
pub fn split_fq(args: &Args) -> Result<()> {
    match resolve_input(args)? {
        ResolvedInput::FastqPath(path) => split_fastq_path(args, &path),
        ResolvedInput::FastqGzPath(path) => split_fastq_gz_path(args, &path),
        ResolvedInput::FastqMemory { data } => split_fastq_bytes(args, data),
        _ => anyhow::bail!("ERROR: input is not FASTQ"),
    }
}

/// Splits 2bit files.
///
/// # Arguments
/// * `args` - Parsed command-line arguments
///
/// # Errors
/// Returns error if input is not 2bit format
pub fn split_2bit(args: &Args) -> Result<()> {
    match resolve_input(args)? {
        ResolvedInput::TwoBitPath(path) => split_2bit_path(args, &path),
        ResolvedInput::TwoBitMemory { data } => split_2bit_bytes(args, data),
        _ => anyhow::bail!("ERROR: input is not 2bit"),
    }
}

/// Resolves input source from args (file path or stdin).
fn resolve_input(args: &Args) -> Result<ResolvedInput> {
    match args.file.as_deref() {
        Some(path) if !is_stdin_path(path) => resolve_file_input(path),
        _ => resolve_stdin_input(),
    }
}

/// Resolves a file path to a ResolvedInput variant.
fn resolve_file_input(path: &Path) -> Result<ResolvedInput> {
    Ok(match detect_file_format(path)? {
        InputFormat::Fasta => ResolvedInput::FastaPath(path.to_path_buf()),
        InputFormat::FastaGz => ResolvedInput::FastaGzPath(path.to_path_buf()),
        InputFormat::Fastq => ResolvedInput::FastqPath(path.to_path_buf()),
        InputFormat::FastqGz => ResolvedInput::FastqGzPath(path.to_path_buf()),
        InputFormat::TwoBit => ResolvedInput::TwoBitPath(path.to_path_buf()),
    })
}

/// Reads stdin and resolves its format.
fn resolve_stdin_input() -> Result<ResolvedInput> {
    log::info!("INFO: reading input from stdin");
    let mut data = Vec::new();
    io::stdin().read_to_end(&mut data)?;

    if data.is_empty() {
        anyhow::bail!("ERROR: stdin is empty");
    }

    if is_twobit_bytes(&data) {
        return Ok(ResolvedInput::TwoBitMemory {
            data: SharedBytes::from_arc(Arc::new(data)),
        });
    }

    if is_gzip_bytes(&data) {
        let decompressed = decompress_gzip(&data)?;
        return Ok(match detect_text_format(&decompressed)? {
            InputFormat::Fasta => ResolvedInput::FastaMemory {
                data: SharedBytes::from_arc(Arc::new(decompressed)),
            },
            InputFormat::Fastq => ResolvedInput::FastqMemory {
                data: SharedBytes::from_arc(Arc::new(decompressed)),
            },
            _ => anyhow::bail!("ERROR: unsupported gzipped stdin format"),
        });
    }

    Ok(match detect_text_format(&data)? {
        InputFormat::Fasta => ResolvedInput::FastaMemory {
            data: SharedBytes::from_arc(Arc::new(data)),
        },
        InputFormat::Fastq => ResolvedInput::FastqMemory {
            data: SharedBytes::from_arc(Arc::new(data)),
        },
        _ => anyhow::bail!("ERROR: unsupported stdin format"),
    })
}

/// Detects input format from file path and magic bytes.
fn detect_file_format(path: &Path) -> Result<InputFormat> {
    if path_ends_with(path, ".2bit") {
        return Ok(InputFormat::TwoBit);
    }
    if path_ends_with(path, ".fa.gz") || path_ends_with(path, ".fasta.gz") {
        return Ok(InputFormat::FastaGz);
    }
    if path_ends_with(path, ".fq.gz") || path_ends_with(path, ".fastq.gz") {
        return Ok(InputFormat::FastqGz);
    }
    if path_ends_with(path, ".fa") || path_ends_with(path, ".fasta") {
        return Ok(InputFormat::Fasta);
    }
    if path_ends_with(path, ".fq") || path_ends_with(path, ".fastq") {
        return Ok(InputFormat::Fastq);
    }

    let mut file = File::open(path)?;
    let mut probe = [0_u8; 4096];
    let count = file.read(&mut probe)?;
    let probe = &probe[..count];

    if is_twobit_bytes(probe) {
        return Ok(InputFormat::TwoBit);
    }

    if is_gzip_bytes(probe) {
        anyhow::bail!(
            "ERROR: gzipped input {} requires a FASTA/FASTQ extension",
            path.display()
        );
    }

    detect_text_format(probe)
}

/// Detects text format (FASTA or FASTQ) from first non-whitespace byte.
fn detect_text_format(data: &[u8]) -> Result<InputFormat> {
    let first = data
        .iter()
        .copied()
        .find(|byte| !byte.is_ascii_whitespace())
        .ok_or_else(|| anyhow!("ERROR: unable to detect input format from empty content"))?;

    match first {
        b'>' => Ok(InputFormat::Fasta),
        b'@' => Ok(InputFormat::Fastq),
        _ => anyhow::bail!("ERROR: unsupported input format"),
    }
}

/// Checks if path represents stdin (dash).
fn is_stdin_path(path: &Path) -> bool {
    path == Path::new("-")
}

/// Checks if data starts with gzip magic bytes.
fn is_gzip_bytes(data: &[u8]) -> bool {
    data.starts_with(&GZIP_MAGIC)
}

/// Checks if data starts with 2bit magic bytes.
fn is_twobit_bytes(data: &[u8]) -> bool {
    data.starts_with(&TWOBIT_MAGIC) || data.starts_with(&TWOBIT_MAGIC_REVERSED)
}

/// Decompresses gzip data in memory.
fn decompress_gzip(data: &[u8]) -> Result<Vec<u8>> {
    let mut decoder = MultiGzDecoder::new(data);
    let mut decompressed = Vec::new();
    decoder.read_to_end(&mut decompressed)?;
    Ok(decompressed)
}

/// Checks if path filename ends with suffix.
fn path_ends_with(path: &Path, suffix: &str) -> bool {
    path.file_name()
        .and_then(|file_name| file_name.to_str())
        .map(|file_name| file_name.ends_with(suffix))
        .unwrap_or(false)
}

/// Returns the filename suffix (e.g. `.gz`) for the given compression codec.
fn compression_suffix(compression: Option<Compression>) -> &'static str {
    compression.map(Compression::extension_suffix).unwrap_or("")
}

/// Computes the output file extension from the output format and compression.
///
/// 2bit never carries a compression suffix.
fn output_extension(format: OutputFormat, compression: Option<Compression>) -> String {
    match format {
        OutputFormat::Fasta => format!("fa{}", compression_suffix(compression)),
        OutputFormat::Fastq => format!("fastq{}", compression_suffix(compression)),
        OutputFormat::TwoBit => "2bit".to_string(),
    }
}

/// Computes the FASTA-input output extension, preserving the longer `fasta`
/// spelling when the identity FASTA path used a `.fasta`/`.fasta.gz` source.
fn fasta_input_output_extension(
    format: OutputFormat,
    compression: Option<Compression>,
    source_path: Option<&Path>,
) -> String {
    if format == OutputFormat::Fasta
        && source_path
            .is_some_and(|path| path_ends_with(path, ".fasta") || path_ends_with(path, ".fasta.gz"))
    {
        format!("fasta{}", compression_suffix(compression))
    } else {
        output_extension(format, compression)
    }
}

/// Resolves and validates the output format for a given natural input format.
///
/// `--output-format` defaults to the natural format. FASTQ output is only
/// permitted when the input is FASTQ (no quality scores exist otherwise).
fn resolve_output_format(natural: OutputFormat, args: &Args) -> Result<OutputFormat> {
    let requested = args.output_format.unwrap_or(natural);
    if requested == OutputFormat::Fastq && natural != OutputFormat::Fastq {
        anyhow::bail!(
            "ERROR: cannot produce FASTQ from {} input (no quality scores)",
            natural.label()
        );
    }
    Ok(requested)
}

/// Determines the effective compression for an output, forcing 2bit outputs
/// to be uncompressed (and warning once if a codec was requested for them).
fn effective_compression(
    format: OutputFormat,
    compression: Option<Compression>,
) -> Option<Compression> {
    if format == OutputFormat::TwoBit {
        if compression.is_some() {
            log::warn!("WARNING: 2bit outputs are never compressed; ignoring --compression");
        }
        None
    } else {
        compression
    }
}

/// Splits plain FASTA file at given path.
fn split_fasta_path(args: &Args, path: &Path) -> Result<()> {
    log::info!("INFO: running in FASTA mode with args: {:?}", args);
    let data = SharedBytes::from_arc(Arc::new(mmap_file(path)?));
    split_fasta_bytes(args, data, Some(path))
}

/// Splits gzipped FASTA file at given path.
fn split_fasta_gz_path(args: &Args, path: &Path) -> Result<()> {
    log::info!("INFO: running in FASTA.GZ mode with args: {:?}", args);
    let data = SharedBytes::from_arc(Arc::new(decompress_gzip_file(path)?));
    split_fasta_bytes(args, data, Some(path))
}

/// Splits plain FASTQ file at given path.
fn split_fastq_path(args: &Args, path: &Path) -> Result<()> {
    log::info!("INFO: running in FASTQ mode with args: {:?}", args);
    let data = SharedBytes::from_arc(Arc::new(mmap_file(path)?));
    split_fastq_bytes(args, data)
}

/// Splits gzipped FASTQ file at given path.
fn split_fastq_gz_path(args: &Args, path: &Path) -> Result<()> {
    log::info!("INFO: running in FASTQ.GZ mode with args: {:?}", args);
    let data = SharedBytes::from_arc(Arc::new(decompress_gzip_file(path)?));
    split_fastq_bytes(args, data)
}

/// Splits 2bit file at given path.
fn split_2bit_path(args: &Args, path: &Path) -> Result<()> {
    log::info!("INFO: running in 2BIT mode with args: {:?}", args);
    let data = SharedBytes::from_arc(Arc::new(mmap_file(path)?));
    split_2bit_bytes(args, data)
}

/// Splits FASTA data into chunks using parallel processing.
///
/// # Arguments
/// * `args` - Command-line arguments
/// * `data` - FASTA data to split
/// * `source_path` - Original file path for passthrough optimization
fn split_fasta_bytes<B>(args: &Args, data: SharedBytes<B>, source_path: Option<&Path>) -> Result<()>
where
    B: AsRef<[u8]> + Send + Sync + 'static,
{
    let mode = args.mode()?;
    let out_fmt = resolve_output_format(OutputFormat::Fasta, args)?;
    let compression = effective_compression(out_fmt, args.compression);
    let extension = fasta_input_output_extension(out_fmt, compression, source_path);
    let pool = build_thread_pool(args.threads)?;
    let headers = find_fasta_headers(data.as_slice());

    if headers.is_empty() {
        anyhow::bail!("ERROR: No FASTA records found");
    }

    log::info!(
        "INFO: Found {} FASTA records in {}",
        headers.len(),
        input_label(args)
    );

    create_dir_all(&args.outdir)?;
    let suffix = args.suffix.clone().unwrap_or_default();

    if out_fmt == OutputFormat::Fasta
        && !args.no_mask
        && passthrough_compression_matches(source_path, compression)
        && maybe_passthrough_original_file(args, source_path, mode, headers.len(), &extension)?
    {
        return Ok(());
    }

    if matches!(mode, SplitMode::FileHeader) {
        warn_if_many_header_outputs(headers.len(), &input_label(args));
        let ids = headers
            .iter()
            .map(|&start| extract_fasta_id(data.as_slice(), start))
            .collect::<Vec<_>>();
        let filenames = make_unique_filenames_from_ids(&ids, &suffix, &extension);

        pool.install(|| {
            headers
                .par_iter()
                .enumerate()
                .try_for_each(|(index, &start)| -> Result<()> {
                    let end = *headers.get(index + 1).unwrap_or(&data.as_slice().len());
                    let output = PathBuf::from(&args.outdir).join(&filenames[index]);
                    write_chunk(
                        &output,
                        &data.as_slice()[start..end],
                        SourceKind::Fasta,
                        out_fmt,
                        compression,
                        args.no_mask,
                    )
                })
        })?;

        return Ok(());
    }

    let ranges = record_ranges_for_mode(headers.len(), mode)?;
    let prefix = output_prefix(mode);

    pool.install(|| {
        ranges
            .into_par_iter()
            .enumerate()
            .try_for_each(|(index, range)| -> Result<()> {
                let start = *headers.get(range.start).unwrap_or(&data.as_slice().len());
                let end = *headers.get(range.end).unwrap_or(&data.as_slice().len());
                let output = PathBuf::from(&args.outdir)
                    .join(numbered_output_name(prefix, index, &suffix, &extension));
                write_chunk(
                    &output,
                    &data.as_slice()[start..end],
                    SourceKind::Fasta,
                    out_fmt,
                    compression,
                    args.no_mask,
                )
            })
    })?;

    Ok(())
}

/// Splits FASTQ data into chunks using parallel processing.
///
/// # Arguments
/// * `args` - Command-line arguments
/// * `data` - FASTQ data to split
fn split_fastq_bytes<B>(args: &Args, data: SharedBytes<B>) -> Result<()>
where
    B: AsRef<[u8]> + Send + Sync + 'static,
{
    let mode = args.mode()?;
    if matches!(mode, SplitMode::FileHeader) {
        anyhow::bail!(
            "ERROR: --headers is only supported for FASTA/FASTA.GZ/2BIT files, not FASTQ/FASTQ.GZ"
        );
    }

    let out_fmt = resolve_output_format(OutputFormat::Fastq, args)?;
    let compression = effective_compression(out_fmt, args.compression);
    let extension = output_extension(out_fmt, compression);
    let pool = build_thread_pool(args.threads)?;
    let suffix = args.suffix.clone().unwrap_or_default();
    let record_starts = find_fastq_record_starts(data.as_slice())?;

    log::info!(
        "INFO: Found {} FASTQ records in {}",
        record_starts.len(),
        input_label(args)
    );

    create_dir_all(&args.outdir)?;
    let ranges = record_ranges_for_mode(record_starts.len(), mode)?;
    let prefix = output_prefix(mode);

    pool.install(|| {
        ranges
            .into_par_iter()
            .enumerate()
            .try_for_each(|(index, range)| -> Result<()> {
                let start = *record_starts
                    .get(range.start)
                    .unwrap_or(&data.as_slice().len());
                let end = *record_starts
                    .get(range.end)
                    .unwrap_or(&data.as_slice().len());
                let output = PathBuf::from(&args.outdir)
                    .join(numbered_output_name(prefix, index, &suffix, &extension));
                write_chunk(
                    &output,
                    &data.as_slice()[start..end],
                    SourceKind::Fastq,
                    out_fmt,
                    compression,
                    args.no_mask,
                )
            })
    })?;

    Ok(())
}

/// Splits 2bit data into chunks using parallel processing.
///
/// # Arguments
/// * `args` - Command-line arguments
/// * `data` - 2bit data to split
fn split_2bit_bytes<B>(args: &Args, data: SharedBytes<B>) -> Result<()>
where
    B: AsRef<[u8]> + Send + Sync + 'static,
{
    let mode = args.mode()?;
    let out_fmt = resolve_output_format(OutputFormat::TwoBit, args)?;
    let compression = effective_compression(out_fmt, args.compression);
    let extension = output_extension(out_fmt, compression);
    let pool = build_thread_pool(args.threads)?;
    let suffix = args.suffix.clone().unwrap_or_default();

    let reader = TwoBitFile::from_buf(data.clone()).map_err(|err| anyhow!(err.to_string()))?;
    let sequence_info = reader.sequence_info();
    if sequence_info.is_empty() {
        anyhow::bail!("ERROR: No 2bit records found");
    }

    let names = sequence_info
        .iter()
        .map(|info| info.chr.clone())
        .collect::<Vec<_>>();
    let lengths = sequence_info
        .iter()
        .map(|info| info.length)
        .collect::<Vec<_>>();

    log::info!(
        "INFO: Found {} 2BIT records in {}",
        names.len(),
        input_label(args)
    );

    create_dir_all(&args.outdir)?;

    if matches!(mode, SplitMode::FileHeader) {
        warn_if_many_header_outputs(names.len(), &input_label(args));
        let filenames = make_unique_filenames_from_ids(&names, &suffix, &extension);

        pool.install(|| {
            names
                .par_iter()
                .enumerate()
                .try_for_each(|(index, _)| -> Result<()> {
                    let output = PathBuf::from(&args.outdir).join(&filenames[index]);
                    write_twobit_group(
                        &output,
                        data.clone(),
                        &names,
                        &lengths,
                        ChunkRegion {
                            start: index,
                            end: index + 1,
                        },
                        out_fmt,
                        compression,
                        !args.no_mask,
                    )
                })
        })?;

        return Ok(());
    }

    let ranges = record_ranges_for_mode(names.len(), mode)?;
    let prefix = output_prefix(mode);

    pool.install(|| {
        ranges
            .into_par_iter()
            .enumerate()
            .try_for_each(|(index, range)| -> Result<()> {
                let output = PathBuf::from(&args.outdir)
                    .join(numbered_output_name(prefix, index, &suffix, &extension));
                write_twobit_group(
                    &output,
                    data.clone(),
                    &names,
                    &lengths,
                    range,
                    out_fmt,
                    compression,
                    !args.no_mask,
                )
            })
    })?;

    Ok(())
}

/// Writes a group of 2bit sequences to an output file.
///
/// # Arguments
/// * `output` - Output file path
/// * `data` - 2bit data source
/// * `names` - Sequence names
/// * `lengths` - Sequence lengths
/// * `range` - Which sequences to write
/// * `out_fmt` - Output format (2bit identity or FASTA conversion)
/// * `compression` - Compression codec for FASTA output (ignored for 2bit)
/// * `preserve_softmask` - Whether to preserve soft-masking
#[allow(clippy::too_many_arguments)]
fn write_twobit_group<B>(
    output: &Path,
    data: SharedBytes<B>,
    names: &[String],
    lengths: &[usize],
    range: ChunkRegion,
    out_fmt: OutputFormat,
    compression: Option<Compression>,
    preserve_softmask: bool,
) -> Result<()>
where
    B: AsRef<[u8]> + Send + Sync + 'static,
{
    let reader = TwoBitFile::from_buf(data).map_err(|err| anyhow!(err.to_string()))?;
    let mut reader = reader.enable_softmask(preserve_softmask);

    let mut records = Vec::with_capacity(range.end.saturating_sub(range.start));
    for index in range.start..range.end {
        if index >= names.len() {
            break;
        }

        let name = names[index].clone();
        let bases = reader
            .read_sequence(&name, ..)
            .map_err(|err| anyhow!(err.to_string()))?
            .into_bytes();
        let hard_blocks = reader
            .hard_masked_blocks(&name, ..)
            .map_err(|err| anyhow!(err.to_string()))?;
        let soft_blocks = if preserve_softmask {
            reader
                .soft_masked_blocks(&name, ..)
                .map_err(|err| anyhow!(err.to_string()))?
        } else {
            Vec::new()
        };

        records.push(TwoBitSequence {
            name,
            length: lengths[index],
            bases,
            hard_blocks,
            soft_blocks,
        });
    }

    match out_fmt {
        OutputFormat::TwoBit => {
            let file = File::create(output)?;
            let mut writer = BufWriter::new(file);
            write_twobit_file(&mut writer, &records)?;
            writer.flush()?;
        }
        OutputFormat::Fasta => {
            let mut writer = SequenceWriter::create(output, compression)?;
            let mut buffer = Vec::new();
            for record in &records {
                twobit_record_to_fasta(record, FASTA_LINE_WIDTH, &mut buffer);
            }
            writer.write_all(&buffer)?;
            writer.finish()?;
        }
        OutputFormat::Fastq => {
            // Unreachable: resolve_output_format rejects FASTQ for 2bit input.
            anyhow::bail!("ERROR: cannot produce FASTQ from 2bit input (no quality scores)");
        }
    }

    Ok(())
}

/// Finds all FASTA header positions (lines starting with '>').
fn find_fasta_headers(data: &[u8]) -> Vec<usize> {
    memchr_iter(FASTA_HEADER, data)
        .filter(|&index| index == 0 || data.get(index - 1) == Some(&b'\n'))
        .collect()
}

/// Extracts the sequence ID from a FASTA header line.
fn extract_fasta_id(data: &[u8], header_start: usize) -> String {
    let start = header_start.saturating_add(1);
    if start >= data.len() {
        return "record".to_string();
    }

    let line_end = data[start..]
        .iter()
        .position(|&byte| byte == b'\n')
        .map(|index| start + index)
        .unwrap_or(data.len());

    let mut line = &data[start..line_end];
    if line.last() == Some(&b'\r') {
        line = &line[..line.len().saturating_sub(1)];
    }

    let token_end = line
        .iter()
        .position(|byte| byte.is_ascii_whitespace())
        .unwrap_or(line.len());

    String::from_utf8_lossy(&line[..token_end]).to_string()
}

/// Sanitizes a record ID for use as a filename.
///
/// Replaces non-alphanumeric characters (except _ - .) with underscores
/// and collapses consecutive underscores.
fn sanitize_record_id(id: &str) -> String {
    let mut out = String::with_capacity(id.len());
    let mut previous_underscore = false;

    for ch in id.chars() {
        let mapped = if ch.is_ascii_alphanumeric() || matches!(ch, '_' | '-' | '.') {
            ch
        } else {
            '_'
        };

        if mapped == '_' {
            if previous_underscore {
                continue;
            }
            previous_underscore = true;
        } else {
            previous_underscore = false;
        }

        out.push(mapped);
    }

    let trimmed = out.trim_matches('_');
    if trimmed.is_empty() {
        "record".to_string()
    } else {
        trimmed.to_string()
    }
}

/// Generates unique filenames from record IDs.
///
/// Handles duplicates by appending numeric suffixes.
fn make_unique_filenames_from_ids(ids: &[String], suffix: &str, extension: &str) -> Vec<String> {
    let mut used = HashSet::with_capacity(ids.len());
    let mut output = Vec::with_capacity(ids.len());

    for id in ids {
        let sanitized = sanitize_record_id(id);
        let base = if suffix.is_empty() {
            sanitized
        } else {
            format!("{sanitized}_{suffix}")
        };

        let mut candidate = base.clone();
        let mut duplicate_index = 1usize;
        while used.contains(&candidate) {
            candidate = format!("{base}_{duplicate_index}");
            duplicate_index += 1;
        }

        used.insert(candidate.clone());
        output.push(format!("{candidate}.{extension}"));
    }

    output
}

/// Warns if --headers mode would create too many output files.
fn warn_if_many_header_outputs(num_records: usize, input: &str) {
    if num_records > HEADER_SPLIT_WARN_THRESHOLD {
        log::warn!(
            "WARNING: input {} contains {} records; --headers will create {} output files",
            input,
            num_records,
            num_records
        );
    }
}

/// Builds a Rayon thread pool for parallel processing.
fn build_thread_pool(threads: usize) -> Result<rayon::ThreadPool> {
    if threads == 0 {
        anyhow::bail!("ERROR: --threads must be greater than 0");
    }

    Ok(rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()?)
}

/// Calculates chunk regions based on split mode.
///
/// # Arguments
/// * `num_records` - Total number of records
/// * `mode` - The split mode to use
fn record_ranges_for_mode(num_records: usize, mode: SplitMode) -> Result<Vec<ChunkRegion>> {
    match mode {
        SplitMode::ChunkSize(records_per_file) => {
            if records_per_file == 0 {
                anyhow::bail!("ERROR: --chunks must be greater than 0");
            }

            let chunks = num_records.div_ceil(records_per_file);
            Ok((0..chunks)
                .map(|index| ChunkRegion {
                    start: index * records_per_file,
                    end: ((index + 1) * records_per_file).min(num_records),
                })
                .collect())
        }
        SplitMode::NumFiles(files) => {
            if files == 0 {
                anyhow::bail!("ERROR: --files must be greater than 0");
            }

            let base = num_records / files;
            let remainder = num_records % files;
            let mut start = 0usize;
            let mut ranges = Vec::with_capacity(files);

            for index in 0..files {
                let width = if index < remainder { base + 1 } else { base };
                let end = (start + width).min(num_records);
                ranges.push(ChunkRegion { start, end });
                start = end;
            }

            Ok(ranges)
        }
        SplitMode::FileHeader => Ok((0..num_records)
            .map(|index| ChunkRegion {
                start: index,
                end: index + 1,
            })
            .collect()),
    }
}

/// Finds FASTQ record start positions (lines starting with '@').
///
/// Validates FASTQ format (4 lines per record).
fn find_fastq_record_starts(data: &[u8]) -> Result<Vec<usize>> {
    let mut starts = Vec::new();
    let mut pos = 0usize;
    let mut line_count = 0usize;

    while pos < data.len() {
        let line_end = line_end(data, pos);
        let mut line = &data[pos..line_end];
        if line.last() == Some(&b'\r') {
            line = &line[..line.len().saturating_sub(1)];
        }

        match line_count % 4 {
            0 => {
                if !line.starts_with(b"@") {
                    anyhow::bail!(
                        "ERROR: malformed FASTQ input: record header does not start with '@'"
                    );
                }
                starts.push(pos);
            }
            2 => {
                if !line.starts_with(b"+") {
                    anyhow::bail!(
                        "ERROR: malformed FASTQ input: plus line does not start with '+'"
                    );
                }
            }
            _ => {}
        }

        line_count += 1;
        pos = next_line_start(data, line_end);
    }

    if starts.is_empty() {
        anyhow::bail!("ERROR: No FASTQ records found");
    }

    if !line_count.is_multiple_of(4) {
        anyhow::bail!(
            "ERROR: malformed FASTQ input: line count {} is not divisible by 4",
            line_count
        );
    }

    Ok(starts)
}

/// Writes one output chunk from a FASTA/FASTQ source slice, converting to the
/// requested output format and applying compression.
///
/// The identity paths (FASTA→FASTA, FASTQ→FASTQ) copy/transform the slice
/// directly; conversions build the target representation in memory first.
fn write_chunk(
    output: &Path,
    chunk: &[u8],
    src: SourceKind,
    out_fmt: OutputFormat,
    compression: Option<Compression>,
    no_mask: bool,
) -> Result<()> {
    match out_fmt {
        OutputFormat::TwoBit => {
            let records = match src {
                SourceKind::Fasta => fasta_slice_to_twobit_records(chunk, no_mask)?,
                SourceKind::Fastq => fastq_slice_to_twobit_records(chunk, no_mask)?,
            };
            let file = File::create(output)?;
            let mut writer = BufWriter::new(file);
            write_twobit_file(&mut writer, &records)?;
            writer.flush()?;
            Ok(())
        }
        OutputFormat::Fasta => {
            let mut writer = SequenceWriter::create(output, compression)?;
            match src {
                SourceKind::Fasta => write_fasta_slice(&mut writer, chunk, no_mask)?,
                SourceKind::Fastq => {
                    let mut buffer = Vec::with_capacity(chunk.len());
                    fastq_slice_to_fasta(chunk, no_mask, &mut buffer);
                    writer.write_all(&buffer)?;
                }
            }
            writer.finish()
        }
        OutputFormat::Fastq => {
            // Only valid for FASTQ sources (validated in resolve_output_format).
            let mut writer = SequenceWriter::create(output, compression)?;
            write_fastq_slice(&mut writer, chunk, no_mask)?;
            writer.finish()
        }
    }
}

/// Returns true when the original input file can be passed through unchanged:
/// the requested compression must match the source file's on-disk compression.
fn passthrough_compression_matches(
    source_path: Option<&Path>,
    compression: Option<Compression>,
) -> bool {
    let Some(path) = source_path else {
        return false;
    };
    let source_is_gzip = path_ends_with(path, ".gz");
    match compression {
        None => !source_is_gzip,
        Some(Compression::Gzip) => source_is_gzip,
        _ => false,
    }
}

/// Converts a FASTQ source slice to FASTA, dropping quality scores.
fn fastq_slice_to_fasta(data: &[u8], uppercase: bool, out: &mut Vec<u8>) {
    let mut pos = 0usize;
    let mut line_index = 0usize;

    while pos < data.len() {
        let end = line_end(data, pos);
        let next = next_line_start(data, end);
        let line = strip_cr(&data[pos..end]);

        match line_index % 4 {
            0 => {
                out.push(FASTA_HEADER);
                let rest = if line.first() == Some(&b'@') {
                    &line[1..]
                } else {
                    line
                };
                out.extend_from_slice(rest);
                out.push(b'\n');
            }
            1 => {
                if uppercase {
                    out.extend(line.iter().map(|byte| byte.to_ascii_uppercase()));
                } else {
                    out.extend_from_slice(line);
                }
                out.push(b'\n');
            }
            _ => {}
        }

        line_index += 1;
        pos = next;
    }
}

/// Converts a FASTA source slice into 2bit sequence records.
fn fasta_slice_to_twobit_records(data: &[u8], no_mask: bool) -> Result<Vec<TwoBitSequence>> {
    let headers = find_fasta_headers(data);
    if headers.is_empty() {
        anyhow::bail!("ERROR: No FASTA records found");
    }

    let mut records = Vec::with_capacity(headers.len());
    let mut seen = HashSet::with_capacity(headers.len());

    for (index, &start) in headers.iter().enumerate() {
        let end = *headers.get(index + 1).unwrap_or(&data.len());
        let name = extract_fasta_id(data, start);
        ensure_unique_twobit_name(&mut seen, &name)?;
        let sequence = concat_fasta_sequence(&data[start..end]);
        records.push(build_twobit_sequence(name, &sequence, no_mask)?);
    }

    Ok(records)
}

/// Converts a FASTQ source slice into 2bit sequence records.
fn fastq_slice_to_twobit_records(data: &[u8], no_mask: bool) -> Result<Vec<TwoBitSequence>> {
    let mut records = Vec::new();
    let mut seen = HashSet::new();
    let mut pos = 0usize;
    let mut line_index = 0usize;
    let mut current_name: Option<String> = None;

    while pos < data.len() {
        let end = line_end(data, pos);
        let next = next_line_start(data, end);
        let line = strip_cr(&data[pos..end]);

        match line_index % 4 {
            0 => {
                let header = if line.first() == Some(&b'@') {
                    &line[1..]
                } else {
                    line
                };
                let token_end = header
                    .iter()
                    .position(|byte| byte.is_ascii_whitespace())
                    .unwrap_or(header.len());
                current_name = Some(String::from_utf8_lossy(&header[..token_end]).to_string());
            }
            1 => {
                let name = current_name.take().unwrap_or_else(|| "record".to_string());
                ensure_unique_twobit_name(&mut seen, &name)?;
                records.push(build_twobit_sequence(name, line, no_mask)?);
            }
            _ => {}
        }

        line_index += 1;
        pos = next;
    }

    if records.is_empty() {
        anyhow::bail!("ERROR: No FASTQ records found");
    }

    Ok(records)
}

/// Writes a single 2bit sequence as a FASTA record, wrapping at `line_width`.
///
/// The record's `bases` already reflect masking decisions (hard-mask `N` and,
/// when softmask is preserved, lowercase soft-mask), so no extra masking is done.
fn twobit_record_to_fasta(record: &TwoBitSequence, line_width: usize, out: &mut Vec<u8>) {
    out.push(FASTA_HEADER);
    out.extend_from_slice(record.name.as_bytes());
    out.push(b'\n');

    for line in record.bases.chunks(line_width.max(1)) {
        out.extend_from_slice(line);
        out.push(b'\n');
    }
}

/// Builds a 2bit sequence record from a raw (concatenated) sequence.
///
/// Detects hard-mask (`N`/`n`) and soft-mask (lowercase) runs. Fails fast on
/// any base not representable in 2bit (anything other than A/C/G/T/U/N).
fn build_twobit_sequence(name: String, sequence: &[u8], no_mask: bool) -> Result<TwoBitSequence> {
    let bases: Vec<u8> = if no_mask {
        sequence
            .iter()
            .map(|byte| byte.to_ascii_uppercase())
            .collect()
    } else {
        sequence.to_vec()
    };

    for (position, &base) in bases.iter().enumerate() {
        if !is_twobit_base(base) {
            anyhow::bail!(
                "ERROR: cannot convert {} to 2bit: code '{}' at position {} is not representable in 2bit (only A/C/G/T/N)",
                name,
                base as char,
                position
            );
        }
    }

    let hard_blocks = find_runs(&bases, |byte| byte == b'N' || byte == b'n');
    let soft_blocks = if no_mask {
        Vec::new()
    } else {
        find_runs(&bases, |byte| byte.is_ascii_lowercase())
    };

    Ok(TwoBitSequence {
        name,
        length: bases.len(),
        bases,
        hard_blocks,
        soft_blocks,
    })
}

/// Returns true if a byte is a nucleotide representable in the 2bit format.
fn is_twobit_base(base: u8) -> bool {
    matches!(
        base,
        b'A' | b'C' | b'G' | b'T' | b'U' | b'N' | b'a' | b'c' | b'g' | b't' | b'u' | b'n'
    )
}

/// Finds maximal runs of bytes satisfying `predicate`, as half-open ranges.
fn find_runs(data: &[u8], predicate: impl Fn(u8) -> bool) -> Vec<Range<usize>> {
    let mut runs = Vec::new();
    let mut start: Option<usize> = None;

    for (index, &byte) in data.iter().enumerate() {
        if predicate(byte) {
            start.get_or_insert(index);
        } else if let Some(begin) = start.take() {
            runs.push(begin..index);
        }
    }

    if let Some(begin) = start {
        runs.push(begin..data.len());
    }

    runs
}

/// Concatenates the sequence body of a single FASTA record (skipping the
/// header line and stripping line breaks), handling multi-line sequences.
fn concat_fasta_sequence(record: &[u8]) -> Vec<u8> {
    let Some(newline) = memchr(b'\n', record) else {
        return Vec::new();
    };

    let body = &record[newline + 1..];
    let mut sequence = Vec::with_capacity(body.len());
    for &byte in body {
        if byte != b'\n' && byte != b'\r' {
            sequence.push(byte);
        }
    }
    sequence
}

/// Errors if a sequence name has already been used within a single 2bit output.
fn ensure_unique_twobit_name(seen: &mut HashSet<String>, name: &str) -> Result<()> {
    if !seen.insert(name.to_string()) {
        anyhow::bail!(
            "ERROR: duplicate sequence name '{}' in a single 2bit output; \
             use --headers or a smaller --chunks to keep names unique per file",
            name
        );
    }
    Ok(())
}

/// Strips a single trailing carriage return from a line slice.
fn strip_cr(line: &[u8]) -> &[u8] {
    if line.last() == Some(&b'\r') {
        &line[..line.len() - 1]
    } else {
        line
    }
}

/// Writes FASTA data with optional uppercase conversion.
fn write_fasta_slice(writer: &mut dyn Write, data: &[u8], uppercase: bool) -> Result<()> {
    if !uppercase {
        writer.write_all(data)?;
        return Ok(());
    }

    let mut transformed = Vec::with_capacity(data.len());
    let mut at_line_start = true;
    let mut in_header = false;

    for &byte in data {
        if at_line_start {
            in_header = byte == b'>';
        }

        transformed.push(if in_header {
            byte
        } else {
            byte.to_ascii_uppercase()
        });

        if byte == b'\n' {
            at_line_start = true;
            in_header = false;
        } else {
            at_line_start = false;
        }
    }

    writer.write_all(&transformed)?;
    Ok(())
}

/// Writes FASTQ data with optional uppercase conversion for sequences only.
fn write_fastq_slice(writer: &mut dyn Write, data: &[u8], uppercase: bool) -> Result<()> {
    if !uppercase {
        writer.write_all(data)?;
        return Ok(());
    }

    let mut transformed = Vec::with_capacity(data.len());
    let mut pos = 0usize;
    let mut line_index = 0usize;

    while pos < data.len() {
        let line_end = line_end(data, pos);
        let next = next_line_start(data, line_end);
        let line = &data[pos..next];

        if line_index % 4 == 1 {
            transformed.extend(line.iter().map(|byte| byte.to_ascii_uppercase()));
        } else {
            transformed.extend_from_slice(line);
        }

        line_index += 1;
        pos = next;
    }

    writer.write_all(&transformed)?;
    Ok(())
}

/// Optimizes output by symlinking/copying original file when records fit in one chunk.
fn maybe_passthrough_original_file(
    args: &Args,
    source_path: Option<&Path>,
    mode: SplitMode,
    num_records: usize,
    extension: &str,
) -> Result<bool> {
    let Some(source_path) = source_path else {
        return Ok(false);
    };

    let SplitMode::ChunkSize(records_per_file) = mode else {
        return Ok(false);
    };

    if num_records > records_per_file {
        return Ok(false);
    }

    let suffix = args.suffix.clone().unwrap_or_default();
    let output =
        PathBuf::from(&args.outdir).join(numbered_output_name("chunk", 0, &suffix, extension));

    log::warn!(
        "Only {} records found, <= --chunks {}, reusing original file...",
        num_records,
        records_per_file
    );

    if args.copy {
        std::fs::copy(source_path, &output)?;
    } else {
        symlink(source_path, &output)?;
    }

    Ok(true)
}

/// Returns the output prefix based on split mode.
fn output_prefix(mode: SplitMode) -> &'static str {
    match mode {
        SplitMode::NumFiles(_) => "part",
        _ => "chunk",
    }
}

/// Generates a numbered output filename.
fn numbered_output_name(prefix: &str, index: usize, suffix: &str, extension: &str) -> String {
    format!("tmp_{prefix}_{index:03}_{suffix}.{extension}")
}

/// Returns a human-readable label for the input source.
fn input_label(args: &Args) -> String {
    args.file
        .as_deref()
        .filter(|path| !is_stdin_path(path))
        .map(|path| path.display().to_string())
        .unwrap_or_else(|| "stdin".to_string())
}

/// Finds the end position of a line (next newline or end of data).
fn line_end(data: &[u8], start: usize) -> usize {
    memchr(b'\n', &data[start..])
        .map(|index| start + index)
        .unwrap_or(data.len())
}

/// Finds the start of the next line after the given line end.
fn next_line_start(data: &[u8], line_end: usize) -> usize {
    if line_end < data.len() {
        line_end + 1
    } else {
        data.len()
    }
}

/// Memory-maps a file for efficient reading.
fn mmap_file(path: &Path) -> Result<Mmap> {
    let file = File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    Ok(unsafe { Mmap::map(&file)? })
}

/// Decompresses a gzipped file into memory.
fn decompress_gzip_file(path: &Path) -> Result<Vec<u8>> {
    let file = File::open(path)?;
    let mut decoder = MultiGzDecoder::new(file);
    let mut decompressed = Vec::new();
    decoder.read_to_end(&mut decompressed)?;
    Ok(decompressed)
}

/// Writes a complete 2bit file with header and sequence records.
fn write_twobit_file(writer: &mut dyn Write, records: &[TwoBitSequence]) -> Result<()> {
    let version = 0_u32;
    let reserved = 0_u32;
    let sequence_count = usize_to_u32(records.len())?;

    for field in [TWOBIT_SIGNATURE, version, sequence_count, reserved] {
        writer.write_all(&field.to_ne_bytes())?;
    }

    let mut index_size = 0_u32;
    let mut record_sizes = 0_u32;
    let mut relative_offsets = HashMap::<&str, u32>::with_capacity(records.len());

    for record in records {
        if !record.name.is_ascii() {
            anyhow::bail!("ERROR: 2bit sequence name is not ASCII: {}", record.name);
        }
        if record.bases.len() != record.length {
            anyhow::bail!(
                "ERROR: 2bit sequence length mismatch for {}: expected {}, saw {}",
                record.name,
                record.length,
                record.bases.len()
            );
        }

        index_size = index_size
            .checked_add(usize_to_u32(
                1 + record.name.len() + std::mem::size_of::<u32>(),
            )?)
            .ok_or_else(|| anyhow!("ERROR: 2bit index size overflow"))?;

        relative_offsets.insert(&record.name, record_sizes);

        let packed_size = usize_to_u32(record.length.div_ceil(4))?;
        let blocks_size = usize_to_u32(
            4 + 4
                + (record.hard_blocks.len() * 4)
                + (record.hard_blocks.len() * 4)
                + 4
                + (record.soft_blocks.len() * 4)
                + (record.soft_blocks.len() * 4)
                + 4,
        )?;
        record_sizes = record_sizes
            .checked_add(blocks_size)
            .and_then(|size| size.checked_add(packed_size))
            .ok_or_else(|| anyhow!("ERROR: 2bit record size overflow"))?;
    }

    let sequence_records_start = 16_u32
        .checked_add(index_size)
        .ok_or_else(|| anyhow!("ERROR: 2bit header size overflow"))?;

    for record in records {
        let name_size = u8::try_from(record.name.len()).map_err(|_| {
            anyhow!(
                "ERROR: 2bit sequence name is too long: {}",
                record.name.len()
            )
        })?;
        writer.write_all(&[name_size])?;
        writer.write_all(record.name.as_bytes())?;
        let offset = sequence_records_start
            .checked_add(*relative_offsets.get(record.name.as_str()).unwrap_or(&0))
            .ok_or_else(|| anyhow!("ERROR: 2bit record offset overflow"))?;
        writer.write_all(&offset.to_ne_bytes())?;
    }

    for record in records {
        writer.write_all(&usize_to_u32(record.length)?.to_ne_bytes())?;
        write_twobit_blocks(writer, &record.hard_blocks)?;
        write_twobit_blocks(writer, &record.soft_blocks)?;
        writer.write_all(&reserved.to_ne_bytes())?;
        write_twobit_bases(writer, &record.bases)?;
    }

    Ok(())
}

/// Writes 2bit block data (start positions and lengths).
fn write_twobit_blocks(writer: &mut dyn Write, blocks: &[Range<usize>]) -> Result<()> {
    writer.write_all(&usize_to_u32(blocks.len())?.to_ne_bytes())?;

    let mut lengths = Vec::with_capacity(blocks.len() * std::mem::size_of::<u32>());
    for block in blocks {
        let start = usize_to_u32(block.start)?;
        let width = block
            .end
            .checked_sub(block.start)
            .ok_or_else(|| anyhow!("ERROR: invalid 2bit block with end before start"))?;
        writer.write_all(&start.to_ne_bytes())?;
        lengths.extend_from_slice(&usize_to_u32(width)?.to_ne_bytes());
    }

    writer.write_all(&lengths)?;
    Ok(())
}

/// Writes 2bit packed bases (4 bases per byte).
fn write_twobit_bases(writer: &mut dyn Write, bases: &[u8]) -> Result<()> {
    for chunk in bases.chunks(4) {
        let mut byte = 0_u8;
        for (index, base) in chunk.iter().enumerate() {
            let nucleotide = twobit_nucleotide(*base)?;
            byte |= nucleotide.bits() << (6 - (index * 2));
        }
        writer.write_all(&[byte])?;
    }

    Ok(())
}

/// Converts usize to u32 with error handling.
fn usize_to_u32(value: usize) -> Result<u32> {
    u32::try_from(value).map_err(|_| anyhow!("ERROR: value {value} does not fit in a 2bit field"))
}

/// Converts ASCII sequence byte to 2bit nucleotide, accepting lowercase n from overlapping masks.
fn twobit_nucleotide(base: u8) -> Result<Nucleotide> {
    if base == b'n' {
        Ok(Nucleotide::N)
    } else {
        base.try_into()
            .map_err(|err: twobit::Error| anyhow!(err.to_string()))
    }
}
