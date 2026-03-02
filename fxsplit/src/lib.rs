#![allow(clippy::while_let_on_iterator)]
#![allow(clippy::doc_overindented_list_items)]

//! Core module for splitting a .fa/.fq file into chunks
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This module contains the main function for splitting .fa/.fq files
//! based on custom requirements in parallel.
//!
//! In short, the module accepts any type of .fa or .fq file
//! and process the reads or sequences inside them in parallel
//! when is possible. Compressed files are also accepted. The
//! user has the ability to specify is the splitting process should
//! be done based on specific chunk sizes or number of files, and
//! the amount of parallelization that should be used in the process.

use std::{
    collections::HashSet,
    fs::{create_dir_all, File},
    io::{BufRead, BufReader, BufWriter, Read, Write},
    os::unix::fs::symlink,
    path::{Path, PathBuf},
    sync::Arc,
};

use anyhow::Result;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use memchr::memchr_iter;
use memmap2::Mmap;
use rayon::prelude::*;

const FA_NEEDLE: u8 = b'>';
const HEADER_SPLIT_WARN_THRESHOLD: usize = 100_000;

pub mod cli;
use cli::Args;

/// Dispatches file processing based on its suffix.
///
/// This macro inspects the file name of the given `Path` and attempts to match
/// its suffix against a predefined set of patterns. For each matched suffix,
/// it executes the corresponding action. If no suffix matches, it logs an error
/// and bails out, indicating an unrecognized file format.
///
/// This provides a convenient way to route different file types to specific
/// processing functions, often based on common bioinformatics file extensions
/// like `.fa.gz` or `.fastq`.
///
/// # Arguments
///
/// * `$file` - An expression that evaluates to a reference to a `Path` or `PathBuf`,
///             representing the input file.
/// * `$suffixes_and_actions` - A block containing `suffix => $action` pairs,
///                             where `suffix` is a literal string to match
///                             against the end of the file name, and `$action`
///                             is a block of code to execute if the suffix matches.
///
/// # Panics
///
/// This macro will `unwrap_or_default()` on `file_name().to_str()`, which
/// means it expects valid UTF-8 file names for matching. Non-UTF-8 file names
/// will result in an empty string for `f`, which might not match any suffix.
///
/// # Errors
///
/// Returns an `anyhow::Error` if no suffix matches the input file's name.
///
/// # Example
///
/// ```rust, ignore
/// use std::path::PathBuf;
/// use anyhow::Result;
///
/// // Assume these functions exist for the example
/// fn process_fasta(path: &PathBuf) -> Result<()> { /* ... */ Ok(()) }
/// fn process_fastq(path: &PathBuf) -> Result<()> { /* ... */ Ok(()) }
///
/// let my_file = PathBuf::from("data/sequences.fa.gz");
///
/// dispatch!(&my_file, {
///     "fa.gz" => process_fasta(&my_file)?,
///     "fq.gz" => process_fastq(&my_file)?,
///     "fa" => process_fasta(&my_file)?,
/// });
/// ```
#[macro_export]
macro_rules! dispatch {
    ($file:expr, { $($suffix:literal => $action:expr),* $(,)?}) => {{
        let f = $file.file_name().and_then(|f| f.to_str()).unwrap_or_default();
        let mut matched = false;
        $(
            if f.ends_with($suffix) {
                matched = true;
                $action
            }
        )*
        if !matched {
            dbg!(f);
            anyhow::bail!("ERROR: unrecognized file format: {}", $file.display());
        }
    }};
}

/// Checks if a file path ends with a specific suffix.
///
/// This helper function checks if the file name portion of a path ends with the given suffix.
/// It handles the case where the file name might not be valid UTF-8.
///
/// # Arguments
///
/// * `path` - The path to check
/// * `suffix` - The suffix to look for (e.g., "fa.gz", "fasta")
///
/// # Returns
///
/// * `true` if the file name ends with the suffix, `false` otherwise
///
/// # Example
///
/// ```rust, ignore
/// use std::path::Path;
///
/// let path = Path::new("data/sequences.fa.gz");
/// assert!(path_ends_with(&path, "fa.gz"));
/// ```
fn path_ends_with(path: &Path, suffix: &str) -> bool {
    path.file_name()
        .and_then(|f| f.to_str())
        .map(|name| name.ends_with(suffix))
        .unwrap_or(false)
}

/// Determines the output file extension for non-gzipped FASTA files.
///
/// # Arguments
///
/// * `path` - The input file path
///
/// # Returns
///
/// * `"fasta"` if the path ends with "fasta", otherwise `"fa"`
///
/// # Example
///
/// ```rust, ignore
/// use std::path::Path;
///
/// let ext = fasta_output_extension(&Path::new("sequences.fasta"));
/// assert_eq!(ext, "fasta");
/// ```
fn fasta_output_extension(path: &Path) -> &'static str {
    if path_ends_with(path, "fasta") {
        "fasta"
    } else {
        "fa"
    }
}

/// Determines the output file extension for gzipped FASTA files.
///
/// # Arguments
///
/// * `path` - The input file path
///
/// # Returns
///
/// * `"fasta.gz"` if the path ends with "fasta.gz", otherwise `"fa.gz"`
///
/// # Example
///
/// ```rust, ignore
/// use std::path::Path;
///
/// let ext = fasta_gz_output_extension(&Path::new("sequences.fasta.gz"));
/// assert_eq!(ext, "fasta.gz");
/// ```
fn fasta_gz_output_extension(path: &Path) -> &'static str {
    if path_ends_with(path, "fasta.gz") {
        "fasta.gz"
    } else {
        "fa.gz"
    }
}

/// Checks if the input file is gzipped based on its extension.
///
/// # Arguments
///
/// * `path` - The input file path
///
/// # Returns
///
/// * `true` if the path ends with ".gz", `false` otherwise
///
/// # Example
///
/// ```rust, ignore
/// use std::path::Path;
///
/// assert!(fastq_is_gz(&Path::new("sequences.fq.gz")));
/// assert!(!fastq_is_gz(&Path::new("sequences.fq")));
/// ```
fn fastq_is_gz(path: &Path) -> bool {
    path_ends_with(path, ".gz")
}

/// Finds all FASTA header positions in the given data.
///
/// Scans the byte array for the FASTA header marker (`>`) and returns
/// the byte positions where headers start. Filters to ensure headers
/// are at the start of the file or after a newline.
///
/// # Arguments
///
/// * `data` - The raw file content as bytes
///
/// # Returns
///
/// * A `Vec<usize>` containing byte positions of all FASTA headers
///
/// # Example
///
/// ```rust, ignore
/// let data = b">header1\nACGT\n>header2\nTGCA";
/// let headers = find_fasta_headers(data);
/// assert_eq!(headers, vec![0, 12]);
/// ```
fn find_fasta_headers(data: &[u8]) -> Vec<usize> {
    memchr_iter(FA_NEEDLE, data)
        .filter(|&i| i == 0 || data.get(i - 1) == Some(&b'\n'))
        .collect()
}

/// Extracts the FASTA record ID from a header line.
///
/// Parses the header starting at the given position and extracts the
/// first token (typically the sequence ID before any whitespace).
///
/// # Arguments
///
/// * `data` - The raw file content as bytes
/// * `header_start` - The byte position where the header starts
///
/// # Returns
///
/// * The extracted ID as a `String`, or "record" if extraction fails
///
/// # Example
///
/// ```rust, ignore
/// let data = b">seq1 description here\nACGT";
/// let id = extract_fasta_id(data, 0);
/// assert_eq!(id, "seq1");
/// ```
fn extract_fasta_id(data: &[u8], header_start: usize) -> String {
    let start = header_start.saturating_add(1);
    if start >= data.len() {
        return "record".to_string();
    }

    let line_end = data[start..]
        .iter()
        .position(|&b| b == b'\n')
        .map(|i| start + i)
        .unwrap_or(data.len());

    let mut line = &data[start..line_end];
    if line.last() == Some(&b'\r') {
        line = &line[..line.len().saturating_sub(1)];
    }

    let token_end = line
        .iter()
        .position(|b| b.is_ascii_whitespace())
        .unwrap_or(line.len());

    let token = &line[..token_end];
    String::from_utf8_lossy(token).to_string()
}

/// Sanitizes a FASTA ID for use as a filename.
///
/// Replaces non-alphanumeric characters (except underscore, hyphen, dot)
/// with underscores. Collapses multiple consecutive underscores into one
/// and trims leading/trailing underscores.
///
/// # Arguments
///
/// * `id` - The raw FASTA ID to sanitize
///
/// # Returns
///
/// * A sanitized string safe for use as a filename
///
/// # Example
///
/// ```rust, ignore
/// let sanitized = sanitize_fasta_id("seq_1 description");
/// assert_eq!(sanitized, "seq_1_description");
/// ```
fn sanitize_fasta_id(id: &str) -> String {
    let mut out = String::with_capacity(id.len());
    let mut prev_is_underscore = false;

    for ch in id.chars() {
        let mapped = if ch.is_ascii_alphanumeric() || matches!(ch, '_' | '-' | '.') {
            ch
        } else {
            '_'
        };

        if mapped == '_' {
            if prev_is_underscore {
                continue;
            }
            prev_is_underscore = true;
        } else {
            prev_is_underscore = false;
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

/// Generates unique filenames for each FASTA header in the input.
///
/// Creates filenames based on sanitized FASTA IDs, handling duplicates
/// by appending numeric suffixes. Each filename includes the specified
/// suffix and file extension.
///
/// # Arguments
///
/// * `data` - The raw file content as bytes
/// * `headers` - Vector of byte positions for each header
/// * `suffix` - Optional suffix to append to filenames
/// * `extension` - File extension to use (e.g., "fa", "fasta.gz")
///
/// # Returns
///
/// * A `Vec<String>` of unique filenames, one per header
///
/// # Example
///
/// ```rust, ignore
/// let data = b">seq1\nACGT\n>seq2\nTGCA";
/// let headers = vec![0, 12];
/// let names = make_unique_header_filenames(data, &headers, "part", "fa");
/// assert_eq!(names, vec!["seq1_part.fa", "seq2_part.fa"]);
/// ```
fn make_unique_header_filenames(
    data: &[u8],
    headers: &[usize],
    suffix: &str,
    extension: &str,
) -> Vec<String> {
    let mut used = HashSet::with_capacity(headers.len());
    let mut names = Vec::with_capacity(headers.len());

    for &header_start in headers {
        let raw_id = extract_fasta_id(data, header_start);
        let sanitized = sanitize_fasta_id(&raw_id);

        let base_name = if suffix.is_empty() {
            sanitized
        } else {
            format!("{sanitized}_{suffix}")
        };

        let mut candidate = base_name.clone();
        let mut duplicate_idx = 1usize;
        while used.contains(&candidate) {
            candidate = format!("{base_name}_{duplicate_idx}");
            duplicate_idx += 1;
        }

        used.insert(candidate.clone());
        names.push(format!("{candidate}.{extension}"));
    }

    names
}

/// Logs a warning if the number of output files would be excessive.
///
/// Warns the user when splitting by header would create more than
/// 100,000 output files, as this may cause filesystem issues.
///
/// # Arguments
///
/// * `num_headers` - The number of headers (and thus output files)
/// * `input` - The input file path for context in the warning
///
/// # Example
///
/// ```rust, ignore
/// use std::path::Path;
///
/// warn_if_many_header_outputs(150000, &Path::new("input.fa"));
/// // Logs warning if logging is enabled
/// ```
fn warn_if_many_header_outputs(num_headers: usize, input: &Path) {
    if num_headers > HEADER_SPLIT_WARN_THRESHOLD {
        log::warn!(
            "WARNING: input {} contains {} FASTA headers; --headers will create {} output files",
            input.display(),
            num_headers,
            num_headers
        );
    }
}

/// Calculates byte regions for each chunk based on the split mode.
///
/// # Arguments
///
/// * `headers` - Vector of byte positions for each FASTA header
/// * `data_len` - Total length of the data in bytes
/// * `mode` - The splitting mode (ChunkSize, NumFiles, or FileHeader)
///
/// # Returns
///
/// * A `Result<Vec<ChunkRegion>>` containing start/end byte positions
///
/// # Errors
///
/// Returns an error if chunk size or number of files is 0
///
/// # Example
///
/// ```rust, ignore
/// use crate::{SplitMode, ChunkRegion};
///
/// let headers = vec![0, 50, 100, 150];
/// let regions = chunk_regions_for_mode(&headers, 200, SplitMode::ChunkSize(2)).unwrap();
/// assert_eq!(regions.len(), 2);
/// ```
fn chunk_regions_for_mode(
    headers: &[usize],
    data_len: usize,
    mode: SplitMode,
) -> Result<Vec<ChunkRegion>> {
    match mode {
        SplitMode::ChunkSize(records_per_file) => {
            if records_per_file == 0 {
                anyhow::bail!("ERROR: --chunks must be greater than 0");
            }

            let nchunks = headers.len().div_ceil(records_per_file);

            Ok((0..nchunks)
                .map(|i| {
                    let start = *headers.get(i * records_per_file).unwrap_or(&data_len);
                    let end = *headers.get((i + 1) * records_per_file).unwrap_or(&data_len);
                    ChunkRegion { start, end }
                })
                .collect())
        }
        SplitMode::NumFiles(files) => {
            if files == 0 {
                anyhow::bail!("ERROR: --files must be greater than 0");
            }

            // let records_per_file = (headers.len() + files - 1) / files;
            let records_per_file = headers.len().div_ceil(files);

            Ok((0..files)
                .map(|i| {
                    let start_idx = i * records_per_file;
                    let end_idx = ((i + 1) * records_per_file).min(headers.len());
                    let start = *headers.get(start_idx).unwrap_or(&data_len);
                    let end = *headers.get(end_idx).unwrap_or(&data_len);
                    ChunkRegion { start, end }
                })
                .collect())
        }
        SplitMode::FileHeader => Ok((0..headers.len())
            .map(|i| {
                let start = headers[i];
                let end = *headers.get(i + 1).unwrap_or(&data_len);
                ChunkRegion { start, end }
            })
            .collect()),
    }
}

/// Splits large sequencing files (FASTA, FASTQ, gzipped versions) into smaller chunks or files.
///
/// This is the main entry point for the `iso-split` functionality. It parses
/// command-line arguments, determines the input file type based on its extension,
/// and dispatches to the appropriate splitting function (`split_fa`, `split_fa_gz`, `split_fq`).
///
/// The splitting logic depends on the `SplitMode` specified in the `Args`
/// (either by a fixed number of records per chunk/file or by a target number of output files).
///
/// # Arguments
///
/// * `args` - A `Vec<String>` representing the command-line arguments passed to the program.
///
/// # Returns
///
/// * `Result<()>` - An `Ok(())` on successful completion, or an `anyhow::Error` if
///                 argument parsing fails, the file format is unrecognized, or
///                 any splitting operation encounters an error.
///
/// # Errors
///
/// * Returns an error if argument parsing (`cli::Args::from`) fails.
/// * Returns an error if the input file's suffix is not recognized by the `dispatch!` macro.
/// * Returns any error propagated from the called splitting functions (`split_fa`, etc.).
///
/// # Example
///
/// ```rust, ignore
/// fn main() -> Result<()> {
///     lib_iso_split(vec!["--file".to_string(), "input.fasta.gz".to_string(), "--chunk-size".to_string(), "1000".to_string()])
///     Ok(())
/// }
/// ```
pub fn lib_iso_split(args: Vec<String>) -> Result<()> {
    let args = cli::Args::from(args);

    dispatch!(&args.file, {
        "fa.gz" => split_fa_gz(&args)?,
        "fasta.gz" =>  split_fa_gz(&args)?,
        "fq.gz" =>  split_fq(&args)?,
        "fastq.gz" =>  split_fq(&args)?,
        "fq" =>  split_fq(&args)?,
        "fastq" =>  split_fq(&args)?,
        "fa" =>  split_fa(&args)?,
        "fasta" =>  split_fa(&args)?,
    });

    Ok(())
}

/// Splits a non-gzipped FASTA file into multiple smaller FASTA files.
///
/// This function reads a FASTA file, identifies the start positions of all records
/// using `memchr_iter` to find `FA_NEEDLE` (typically `>`). It then divides the
/// file's content (memory-mapped for efficiency) into chunks based on the
/// specified `SplitMode` (either `ChunkSize` or `NumFiles`). Each chunk is then
/// written to a new output file within the designated output directory, utilizing
/// a Rayon thread pool for parallel processing.
///
/// # Arguments
///
/// * `args` - A reference to an `Args` struct containing the input file path,
///            output directory, number of threads, splitting mode, and an optional suffix.
///
/// # Returns
///
/// * `Result<()>` - An `Ok(())` on successful completion, or an `anyhow::Error` if
///                 any operation (file opening, memory mapping, directory creation,
///                 writing to files, or thread pool building) fails.
///
/// # Errors
///
/// * Returns an error if the input FASTA file cannot be opened or memory-mapped.
/// * Returns an error if no FASTA records are found in the input file.
/// * Returns an error if the output directory cannot be created.
/// * Returns an error if `SplitMode::NumFiles` is 0.
/// * Returns any `std::io::Error` during file writing.
///
/// # Parallelism
///
/// This function uses `rayon` for parallel processing of chunks, improving
/// performance for large files. The number of threads is configured via `args.threads`.
///
/// # Example
///
/// ```rust, ignore
/// use anyhow::Result;
/// use std::path::PathBuf;
/// // Assuming Args and SplitMode are defined as in lib_iso_split example
///
/// fn main() -> Result<()> {
///     let args = cli::Args {
///         file: PathBuf::from("input.fa"),
///         outdir: PathBuf::from("fa_chunks"),
///         threads: 4,
///         suffix: Some("part".to_string()),
///         mode_chunk_size: Some(100), // Split into chunks of 100 records
///         mode_num_files: None,
///     };
///     // split_fa(&args)?;
///     println!("Successfully split FASTA file.");
///     Ok(())
/// }
/// ```
pub fn split_fa(args: &Args) -> Result<()> {
    log::info!("INFO: running in FASTA mode with args: {:?}", &args);

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()?;

    let file = File::open(&args.file)?;
    let mmap = unsafe { Mmap::map(&file)? };
    let data = Arc::new(mmap);

    let header = find_fasta_headers(data.as_ref());
    if header.is_empty() {
        anyhow::bail!("ERROR: No FASTA records found");
    }

    let mode = args.mode()?;
    let suffix = args.suffix.clone().unwrap_or_default();
    let extension = fasta_output_extension(&args.file);
    create_dir_all(&args.outdir)?;

    if matches!(mode, SplitMode::FileHeader) {
        warn_if_many_header_outputs(header.len(), &args.file);
        let filenames = make_unique_header_filenames(data.as_ref(), &header, &suffix, extension);

        pool.install(|| {
            (0..header.len())
                .into_par_iter()
                .try_for_each(|i| -> Result<()> {
                    let start = header[i];
                    let end = *header.get(i + 1).unwrap_or(&data.len());
                    let output_file = PathBuf::from(&args.outdir).join(&filenames[i]);
                    let mut writer = BufWriter::new(File::create(output_file)?);
                    writer.write_all(&data[start..end])?;
                    writer.flush()?;
                    Ok(())
                })
        })?;

        return Ok(());
    }

    let chunks = chunk_regions_for_mode(&header, data.len(), mode)?;
    let prefix = match mode {
        SplitMode::NumFiles(_) => "part",
        _ => "chunk",
    };

    pool.install(|| {
        chunks
            .into_par_iter()
            .enumerate()
            .try_for_each(|(i, chunk)| -> Result<()> {
                let output_file = PathBuf::from(&args.outdir)
                    .join(format!("tmp_{prefix}_{:03}_{suffix}.{extension}", i));
                let mut writer = BufWriter::new(File::create(output_file)?);
                writer.write_all(&data[chunk.start..chunk.end])?;
                writer.flush()?;
                Ok(())
            })
    })?;

    Ok(())
}

/// Splits a gzipped FASTA file (`.fa.gz` or `.fasta.gz`) into multiple smaller gzipped FASTA files.
///
/// This function first decompresses the entire gzipped input file into memory.
/// It then identifies the start positions of all FASTA records using `memchr_iter`.
/// The decompressed data is divided into chunks based on the `SplitMode` (either
/// `ChunkSize` for records per file or `NumFiles` for a target number of files).
/// Each chunk is then individually compressed and written to a new gzipped output file
/// within the specified output directory, leveraging a global Rayon thread pool for parallel execution.
///
/// A special case is handled: if the total number of records is less than or equal to
/// the `records_per_file` when in `ChunkSize` mode, it creates a symlink to the original file
/// instead of splitting, to avoid unnecessary work.
///
/// # Arguments
///
/// * `args` - A reference to an `Args` struct containing the input file path,
///            output directory, number of threads, splitting mode, and an optional suffix.
///
/// # Returns
///
/// * `Result<()>` - An `Ok(())` on successful completion, or an `anyhow::Error` if
///                 any operation (file opening, decompression, directory creation,
///                 writing to files, or thread pool building) fails.
///
/// # Errors
///
/// * Returns an error if the input gzipped FASTA file cannot be opened or decompressed.
/// * Returns an error if no FASTA records are found in the decompressed data.
/// * Returns an error if the output directory cannot be created.
/// * Returns an error if `SplitMode::NumFiles` is 0.
/// * Returns any `std::io::Error` during file writing or compression.
///
/// # Parallelism
///
/// This function uses `rayon` for parallel processing of chunks, improving
/// performance for large files. It builds a global thread pool with `args.threads`.
///
/// # Example
///
/// ```rust, ignore
/// use anyhow::Result;
/// use std::path::PathBuf;
/// // Assuming Args and SplitMode are defined as in lib_iso_split example
///
/// fn main() -> Result<()> {
///     let args = cli::Args {
///         file: PathBuf::from("input.fa.gz"),
///         outdir: PathBuf::from("fa_gz_chunks"),
///         threads: 4,
///         suffix: Some("part".to_string()),
///         mode_num_files: Some(10), // Split into 10 output files
///         mode_chunk_size: None,
///     };
///     // split_fa_gz(&args)?;
///     println!("Successfully split gzipped FASTA file.");
///     Ok(())
/// }
/// ```
pub fn split_fa_gz(args: &Args) -> Result<()> {
    log::info!("INFO: running in FASTA.GZ mode with args: {:?}", &args);

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()?;

    let mut gz = MultiGzDecoder::new(File::open(&args.file)?);
    let mut decompressed = Vec::new();
    gz.read_to_end(&mut decompressed)?;
    let data = Arc::new(decompressed);

    let header = find_fasta_headers(data.as_ref());
    log::info!(
        "INFO: Found {} records in {}",
        header.len(),
        &args.file.display()
    );

    if header.is_empty() {
        log::error!("ERROR: No FASTA records found in decompressed data");
        anyhow::bail!("ERROR: No FASTA records found in decompressed data");
    }

    let mode = args.mode()?;
    let suffix = args.suffix.clone().unwrap_or_default();
    let extension = fasta_gz_output_extension(&args.file);
    create_dir_all(&args.outdir)?;

    if let SplitMode::ChunkSize(records_per_file) = mode {
        if records_per_file == 0 {
            anyhow::bail!("ERROR: --chunks must be greater than 0");
        }

        if header.len() <= records_per_file {
            let outpath =
                PathBuf::from(&args.outdir).join(format!("tmp_chunk_000_{suffix}.{extension}"));
            std::fs::create_dir_all(&args.outdir)?;
            log::warn!(
                "Only {} records found, <= --chunks {}, creating symlink to original file...",
                header.len(),
                records_per_file
            );
            symlink(&args.file, &outpath)?;
            return Ok(());
        }
    }

    if matches!(mode, SplitMode::FileHeader) {
        warn_if_many_header_outputs(header.len(), &args.file);
        let filenames = make_unique_header_filenames(data.as_ref(), &header, &suffix, extension);

        pool.install(|| {
            (0..header.len())
                .into_par_iter()
                .try_for_each(|i| -> Result<()> {
                    let start = header[i];
                    let end = *header.get(i + 1).unwrap_or(&data.len());
                    let output_file = PathBuf::from(&args.outdir).join(&filenames[i]);
                    let file = File::create(output_file)?;
                    let writer = BufWriter::new(file);
                    let mut encoder =
                        flate2::write::GzEncoder::new(writer, flate2::Compression::fast());
                    encoder.write_all(&data[start..end])?;
                    encoder.finish()?;
                    Ok(())
                })
        })?;

        return Ok(());
    }

    let chunks = chunk_regions_for_mode(&header, data.len(), mode)?;
    let prefix = match mode {
        SplitMode::NumFiles(_) => "part",
        _ => "chunk",
    };

    pool.install(|| {
        chunks
            .into_par_iter()
            .enumerate()
            .try_for_each(|(i, chunk)| -> Result<()> {
                let output_file = PathBuf::from(&args.outdir)
                    .join(format!("tmp_{prefix}_{:03}_{suffix}.{extension}", i));
                let file = File::create(output_file)?;
                let writer = BufWriter::new(file);
                let mut encoder =
                    flate2::write::GzEncoder::new(writer, flate2::Compression::fast());
                encoder.write_all(&data[chunk.start..chunk.end])?;
                encoder.finish()?;
                Ok(())
            })
    })?;

    Ok(())
}

/// Split mode
///
/// Specifies how to split the input file:
/// - `ChunkSize(n)`: Split into files with at most n records each
/// - `NumFiles(n)`: Split into exactly n output files
/// - `FileHeader`: Create one output file per FASTA record
///
/// # Arguments
///
/// * `ChunkSize(usize)` - Number of records per output file
/// * `NumFiles(usize)` - Number of output files to create
/// * `FileHeader` - Split by individual FASTA headers
///
/// # Example
///
/// ```rust, ignore
/// use fxsplit::SplitMode;
///
/// let mode = SplitMode::ChunkSize(1000);
/// let mode = SplitMode::NumFiles(10);
/// let mode = SplitMode::FileHeader;
/// ```
#[derive(Debug, Clone, Copy)]
pub enum SplitMode {
    ChunkSize(usize), // N records per output file
    NumFiles(usize),  // K output files
    FileHeader,       // split by file header
}

// public structs
/// Region [iso-split]
///
/// This enum is used to store region boundaries
/// of a fq/fa file
///
/// # Example
///
/// ```rust, ignore
/// use iso::Region;
///
/// let region = Region{start: 1, end: 43};
/// assert_eq!(43, region.end);
/// ```
#[derive(Debug, Clone)]
pub struct ChunkRegion {
    pub start: usize,
    pub end: usize,
}

/// Splits a gzipped FASTQ file (`.fq.gz` or `.fastq.gz`) into multiple smaller gzipped FASTQ files.
///
/// This function acts as a dispatcher for FASTQ splitting, determining the
/// records-per-file based on the `SplitMode` (`ChunkSize` or `NumFiles`).
/// If splitting by `NumFiles`, it first counts all records in the input file
/// to calculate the appropriate `records_per_file` for even distribution.
/// It then delegates the actual splitting and writing to `split_by_chunk_size`.
///
/// # Arguments
///
/// * `args` - A reference to an `Args` struct containing the input file path,
///            output directory, number of threads (though not directly used here,
///            passed to underlying functions), splitting mode, and an optional suffix.
///
/// # Returns
///
/// * `anyhow::Result<()>` - An `Ok(())` on successful completion, or an `anyhow::Error` if
///                         the splitting mode is invalid, record counting fails, or
///                         `split_by_chunk_size` encounters an error.
///
/// # Errors
///
/// * Returns an error if `args.mode()` returns an error (e.g., no split mode specified).
/// * Returns an error if `count_fastq_records` fails.
/// * Propagates any error from `split_by_chunk_size`.
///
/// # Example
///
/// ```rust, ignore
/// use anyhow::Result;
/// use std::path::PathBuf;
/// // Assuming Args and SplitMode are defined as in lib_iso_split example
///
/// fn main() -> Result<()> {
///     let args = cli::Args {
///         file: PathBuf::from("input.fastq.gz"),
///         outdir: PathBuf::from("fq_chunks"),
///         threads: 4, // Not directly used in split_fq, but passed to helpers
///         suffix: Some("batch".to_string()),
///         mode_chunk_size: None,
///         mode_num_files: Some(5), // Split into 5 output files
///     };
///     // split_fq(&args)?;
///     println!("Successfully split FASTQ file.");
///     Ok(())
/// }
/// ```
pub fn split_fq(args: &Args) -> anyhow::Result<()> {
    log::info!("INFO: running in FASTQ mode with args: {:?}", &args);

    let mode = args.mode()?;
    let suffix = args.suffix.clone().unwrap_or_default();

    match mode {
        SplitMode::ChunkSize(records_per_file) => {
            if records_per_file == 0 {
                anyhow::bail!("ERROR: --chunks must be greater than 0");
            }
            log::info!("INFO: splitting by chunk size!");
            split_fastq_by_chunk_size(&args.file, &args.outdir, records_per_file, &suffix)
        }
        SplitMode::NumFiles(num_files) => {
            if num_files == 0 {
                anyhow::bail!("ERROR: --files must be greater than 0");
            }
            log::info!("INFO: splitting by number of files!");
            let total_records = count_fastq_records(&args.file)?;
            split_fastq_by_num_files(&args.file, &args.outdir, num_files, total_records, &suffix)
        }
        SplitMode::FileHeader => anyhow::bail!(
            "ERROR: --headers is only supported for FASTA/FASTA.GZ files, not FASTQ/FASTQ.GZ"
        ),
    }
}

/// Represents a single FASTQ record.
///
/// # Fields
///
/// * `header` - The FASTQ header line (starts with @)
/// * `seq` - The nucleotide sequence
/// * `plus` - The plus line (usually just "+")
/// * `qual` - Quality scores string
///
/// # Example
///
/// ```rust, ignore
/// let record = FastqRecord {
///     header: "@seq1".to_string(),
///     seq: "ACGT".to_string(),
///     plus: "+".to_string(),
///     qual: "IIII".to_string(),
/// };
/// ```
struct FastqRecord {
    header: String,
    seq: String,
    plus: String,
    qual: String,
}

/// Writer type for FASTQ output files.
///
/// Handles both plain and gzipped FASTQ file output.
///
/// # Variants
///
/// * `Plain` - Plain text FASTQ writer
/// * `Gzip` - Gzipped FASTQ writer
enum FastqWriter {
    Plain(BufWriter<File>),
    Gzip(BufWriter<GzEncoder<File>>),
}

impl FastqWriter {
    /// Creates a new FASTQ writer for the specified file index.
    ///
    /// # Arguments
    ///
    /// * `out_dir` - Output directory path
    /// * `index` - File index number (used in filename)
    /// * `suffix` - Optional suffix for the filename
    /// * `gzip_output` - Whether to use gzip compression
    ///
    /// # Returns
    ///
    /// * `anyhow::Result<Self>` - The created writer
    ///
    /// # Example
    ///
    /// ```rust, ignore
    /// use std::path::Path;
    ///
    /// let writer = FastqWriter::create(Path::new("output"), 0, "part", true).unwrap();
    /// ```
    fn create(
        out_dir: &Path,
        index: usize,
        suffix: &str,
        gzip_output: bool,
    ) -> anyhow::Result<Self> {
        let extension = if gzip_output { "fastq.gz" } else { "fastq" };
        let path = out_dir.join(format!("tmp_chunk_{index:03}_{suffix}.{extension}"));
        let file = File::create(path)?;

        if gzip_output {
            let encoder = GzEncoder::new(file, Compression::fast());
            Ok(FastqWriter::Gzip(BufWriter::new(encoder)))
        } else {
            Ok(FastqWriter::Plain(BufWriter::new(file)))
        }
    }

    /// Writes a single FASTQ record to the output.
    ///
    /// # Arguments
    ///
    /// * `record` - The FastqRecord to write
    ///
    /// # Returns
    ///
    /// * `anyhow::Result<()>` - Success or error
    ///
    /// # Example
    ///
    /// ```rust, ignore
    /// use crate::FastqRecord;
    ///
    /// let record = FastqRecord {
    ///     header: "@seq1".to_string(),
    ///     seq: "ACGT".to_string(),
    ///     plus: "+".to_string(),
    ///     qual: "IIII".to_string(),
    /// };
    /// // writer.write_record(&record).unwrap();
    /// ```
    fn write_record(&mut self, record: &FastqRecord) -> anyhow::Result<()> {
        match self {
            FastqWriter::Plain(writer) => {
                writer.write_all(record.header.as_bytes())?;
                writer.write_all(b"\n")?;
                writer.write_all(record.seq.as_bytes())?;
                writer.write_all(b"\n")?;
                writer.write_all(record.plus.as_bytes())?;
                writer.write_all(b"\n")?;
                writer.write_all(record.qual.as_bytes())?;
                writer.write_all(b"\n")?;
                Ok(())
            }
            FastqWriter::Gzip(writer) => {
                writer.write_all(record.header.as_bytes())?;
                writer.write_all(b"\n")?;
                writer.write_all(record.seq.as_bytes())?;
                writer.write_all(b"\n")?;
                writer.write_all(record.plus.as_bytes())?;
                writer.write_all(b"\n")?;
                writer.write_all(record.qual.as_bytes())?;
                writer.write_all(b"\n")?;
                Ok(())
            }
        }
    }

    /// Finalizes and closes the writer, flushing any buffers.
    ///
    /// For gzip writers, also finishes the gzip stream.
    ///
    /// # Returns
    ///
    /// * `anyhow::Result<()>` - Success or error
    ///
    /// # Example
    ///
    /// ```rust, ignore
    /// // writer.finalize().unwrap();
    /// ```
    fn finalize(self) -> anyhow::Result<()> {
        match self {
            FastqWriter::Plain(mut writer) => {
                writer.flush()?;
                Ok(())
            }
            FastqWriter::Gzip(mut writer) => {
                writer.flush()?;
                let encoder = writer
                    .into_inner()
                    .map_err(|e| anyhow::anyhow!("ERROR: failed to flush gzip writer: {}", e))?;
                let _ = encoder.finish()?;
                Ok(())
            }
        }
    }
}

/// Opens a FASTQ file for reading, handling gzip decompression.
///
/// # Arguments
///
/// * `input` - Path to the input FASTQ file
///
/// # Returns
///
/// * `anyhow::Result<Box<dyn BufRead>>` - A buffered reader
///
/// # Example
///
/// ```rust, ignore
/// use std::path::Path;
///
/// let reader = open_fastq_reader(Path::new("input.fq.gz")).unwrap();
/// ```
fn open_fastq_reader(input: &Path) -> anyhow::Result<Box<dyn BufRead>> {
    let file = File::open(input)?;
    if fastq_is_gz(input) {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

/// Reads the next FASTQ record from an iterator of lines.
///
/// # Arguments
///
/// * `lines` - Iterator of lines from the FASTQ file
///
/// # Returns
///
/// * `anyhow::Result<Option<FastqRecord>>` - The record or None at EOF
///
/// # Errors
///
/// Returns an error if the FASTQ format is malformed
///
/// # Example
///
/// ```rust, ignore
/// use std::io::Cursor;
///
/// let data = "@seq1\nACGT\n+\nIIII\n@seq2\nTGCA\n+\nJJJJ";
/// let mut lines = Cursor::new(data).lines();
/// let record = next_fastq_record(&mut lines).unwrap();
/// ```
fn next_fastq_record<I>(lines: &mut I) -> anyhow::Result<Option<FastqRecord>>
where
    I: Iterator<Item = std::io::Result<String>>,
{
    let header = match lines.next() {
        Some(h) => h?,
        None => return Ok(None),
    };

    let seq = lines
        .next()
        .ok_or_else(|| anyhow::anyhow!("ERROR: unexpected EOF in FASTQ record"))??;
    let plus = lines
        .next()
        .ok_or_else(|| anyhow::anyhow!("ERROR: unexpected EOF in FASTQ record"))??;
    let qual = lines
        .next()
        .ok_or_else(|| anyhow::anyhow!("ERROR: unexpected EOF in FASTQ record"))??;

    Ok(Some(FastqRecord {
        header,
        seq,
        plus,
        qual,
    }))
}

/// Splits a FASTQ file into chunks of a fixed number of records.
///
/// # Arguments
///
/// * `input` - Path to the input FASTQ file
/// * `out_dir` - Output directory path
/// * `records_per_file` - Maximum records per output file
/// * `suffix` - Optional suffix for output filenames
///
/// # Returns
///
/// * `anyhow::Result<()>` - Success or error
///
/// # Example
///
/// ```rust, ignore
/// use std::path::Path;
///
/// split_fastq_by_chunk_size(
///     Path::new("input.fq.gz"),
///     Path::new("output"),
///     1000,
///     "part"
/// ).unwrap();
/// ```
fn split_fastq_by_chunk_size(
    input: &Path,
    out_dir: &Path,
    records_per_file: usize,
    suffix: &str,
) -> anyhow::Result<()> {
    create_dir_all(out_dir)?;

    let gzip_output = fastq_is_gz(input);
    let mut lines = open_fastq_reader(input)?.lines();
    let mut writer: Option<FastqWriter> = None;
    let mut file_index = 0usize;
    let mut record_index = 0usize;

    while let Some(record) = next_fastq_record(&mut lines)? {
        if writer.is_none() {
            writer = Some(FastqWriter::create(
                out_dir,
                file_index,
                suffix,
                gzip_output,
            )?);
        } else if record_index == records_per_file {
            if let Some(current_writer) = writer.take() {
                current_writer.finalize()?;
            }
            file_index += 1;
            record_index = 0;
            writer = Some(FastqWriter::create(
                out_dir,
                file_index,
                suffix,
                gzip_output,
            )?);
        }

        if let Some(current_writer) = writer.as_mut() {
            current_writer.write_record(&record)?;
        }
        record_index += 1;
    }

    if let Some(current_writer) = writer {
        current_writer.finalize()?;
    }

    Ok(())
}

/// Splits a FASTQ file into a specified number of output files.
///
/// Distributes records as evenly as possible across the output files.
///
/// # Arguments
///
/// * `input` - Path to the input FASTQ file
/// * `out_dir` - Output directory path
/// * `num_files` - Number of output files to create
/// * `total_records` - Total number of records in the input
/// * `suffix` - Optional suffix for output filenames
///
/// # Returns
///
/// * `anyhow::Result<()>` - Success or error
///
/// # Example
///
/// ```rust, ignore
/// use std::path::Path;
///
/// split_fastq_by_num_files(
///     Path::new("input.fq.gz"),
///     Path::new("output"),
///     10,
///     5000,
///     "part"
/// ).unwrap();
/// ```
fn split_fastq_by_num_files(
    input: &Path,
    out_dir: &Path,
    num_files: usize,
    total_records: usize,
    suffix: &str,
) -> anyhow::Result<()> {
    create_dir_all(out_dir)?;

    let gzip_output = fastq_is_gz(input);
    let mut lines = open_fastq_reader(input)?.lines();

    let base_records = total_records / num_files;
    let remainder = total_records % num_files;

    for file_idx in 0..num_files {
        let records_for_file = if file_idx < remainder {
            base_records + 1
        } else {
            base_records
        };

        let mut writer = FastqWriter::create(out_dir, file_idx, suffix, gzip_output)?;
        for _ in 0..records_for_file {
            let record = next_fastq_record(&mut lines)?
                .ok_or_else(|| anyhow::anyhow!("ERROR: unexpected EOF while splitting FASTQ"))?;
            writer.write_record(&record)?;
        }
        writer.finalize()?;
    }

    if next_fastq_record(&mut lines)?.is_some() {
        anyhow::bail!("ERROR: FASTQ reader still has remaining records after splitting");
    }

    Ok(())
}

/// Counts the total number of records in a FASTQ file.
///
/// FASTQ files have 4 lines per record, so this counts lines and divides by 4.
///
/// # Arguments
///
/// * `input` - Path to the input FASTQ file
///
/// # Returns
///
/// * `anyhow::Result<usize>` - The total number of records
///
/// # Errors
///
/// Returns an error if the file is malformed (line count not divisible by 4)
///
/// # Example
///
/// ```rust, ignore
/// use std::path::Path;
///
/// let count = count_fastq_records(Path::new("input.fq.gz")).unwrap();
/// ```
fn count_fastq_records(input: &Path) -> anyhow::Result<usize> {
    let mut lines = open_fastq_reader(input)?.lines();
    let mut line_count = 0usize;

    while let Some(line) = lines.next() {
        line?;
        line_count += 1;
    }

    if !line_count.is_multiple_of(4) {
        anyhow::bail!(
            "ERROR: malformed FASTQ input: line count {} is not divisible by 4",
            line_count
        );
    }

    Ok(line_count / 4)
}
