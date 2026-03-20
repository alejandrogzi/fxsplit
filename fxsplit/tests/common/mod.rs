#![allow(dead_code)]

use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::fs::{self, File};
use std::io::{BufWriter, Read, Write};
use std::ops::Range;
use std::path::{Path, PathBuf};
use tempfile::TempDir;
use twobit::{
    convert::{fasta::FastaReader, to_2bit},
    TwoBitFile,
};

/// Creates a temporary directory for test files.
pub fn tempdir() -> TempDir {
    tempfile::tempdir().expect("failed to create tempdir")
}

/// Writes text content to a file.
pub fn write_text(path: &Path, content: &str) {
    fs::write(path, content).expect("failed to write text file");
}

/// Writes gzipped content to a file.
pub fn write_gzip(path: &Path, content: &str) {
    let file = File::create(path).expect("failed to create gzip output");
    let mut encoder = GzEncoder::new(file, Compression::default());
    encoder
        .write_all(content.as_bytes())
        .expect("failed to write gzip content");
    encoder.finish().expect("failed to finish gzip encoder");
}

/// Writes FASTA content as a 2bit file.
pub fn write_2bit(path: &Path, fasta_content: &str) {
    let reader = FastaReader::mem_open(fasta_content.as_bytes().to_vec())
        .expect("failed to build fasta reader");
    let file = File::create(path).expect("failed to create 2bit file");
    let mut writer = BufWriter::new(file);
    to_2bit(&mut writer, &reader).expect("failed to convert fasta to 2bit");
    writer.flush().expect("failed to flush 2bit writer");
}

pub fn write_2bit_with_masks(
    path: &Path,
    name: &str,
    bases: &str,
    hard_blocks: &[Range<usize>],
    soft_blocks: &[Range<usize>],
) {
    let mut writer = BufWriter::new(File::create(path).expect("failed to create 2bit file"));
    let signature = 0x1A41_2743_u32;
    let version = 0_u32;
    let sequence_count = 1_u32;
    let reserved = 0_u32;
    let index_size = 1 + name.len() + 4;
    let sequence_offset = 16 + index_size as u32;

    for field in [signature, version, sequence_count, reserved] {
        writer
            .write_all(&field.to_ne_bytes())
            .expect("failed to write 2bit header");
    }

    writer
        .write_all(&[u8::try_from(name.len()).expect("name too long")])
        .expect("failed to write 2bit name len");
    writer
        .write_all(name.as_bytes())
        .expect("failed to write 2bit name");
    writer
        .write_all(&sequence_offset.to_ne_bytes())
        .expect("failed to write 2bit offset");

    writer
        .write_all(&(bases.len() as u32).to_ne_bytes())
        .expect("failed to write 2bit length");
    write_blocks(&mut writer, hard_blocks);
    write_blocks(&mut writer, soft_blocks);
    writer
        .write_all(&reserved.to_ne_bytes())
        .expect("failed to write 2bit reserved");
    write_packed_bases(&mut writer, bases.as_bytes());
    writer.flush().expect("failed to flush 2bit writer");
}

/// Runs fxsplit and expects success.
pub fn run_split(args: Vec<String>) {
    fxsplit::lib_iso_split(args).expect("split command failed");
}

/// Runs fxsplit and expects an error, returning the error message.
pub fn run_split_expect_err(args: Vec<String>) -> String {
    let err = fxsplit::lib_iso_split(args).expect_err("split command unexpectedly succeeded");
    err.to_string()
}

/// Lists all files in a directory, sorted alphabetically.
pub fn list_files(dir: &Path) -> Vec<String> {
    let mut out: Vec<String> = fs::read_dir(dir)
        .expect("failed to read output dir")
        .map(|e| {
            e.expect("bad dir entry")
                .file_name()
                .to_string_lossy()
                .to_string()
        })
        .collect();
    out.sort();
    out
}

/// Reads text content from a file (optionally gzipped).
pub fn read_text(path: &Path, gz: bool) -> String {
    read_all(path, gz)
}

/// Counts FASTA records (header lines starting with '>').
pub fn count_fasta_records(path: &Path, gz: bool) -> usize {
    read_all(path, gz)
        .lines()
        .filter(|line| line.starts_with('>'))
        .count()
}

/// Counts FASTQ records (4 lines per record).
pub fn count_fastq_records(path: &Path, gz: bool) -> usize {
    read_all(path, gz).lines().count() / 4
}

/// Counts sequences in a 2bit file.
pub fn count_2bit_records(path: &Path) -> usize {
    let tb = TwoBitFile::open(path).expect("failed to open 2bit file");
    tb.chrom_names().len()
}

/// Reads all sequences from a 2bit file.
pub fn read_2bit_sequences(path: &Path, preserve_mask: bool) -> Vec<(String, String)> {
    let tb = TwoBitFile::open(path).expect("failed to open 2bit file");
    let mut tb = tb.enable_softmask(preserve_mask);
    let names = tb.chrom_names();
    names
        .into_iter()
        .map(|name| {
            let seq = tb
                .read_sequence(&name, ..)
                .expect("failed to read 2bit sequence");
            (name, seq)
        })
        .collect()
}

/// Converts a path to a string.
pub fn path_str(path: &Path) -> String {
    path.to_string_lossy().to_string()
}

/// Joins a directory path with a filename.
pub fn join(dir: &Path, name: &str) -> PathBuf {
    dir.join(name)
}

/// Reads all content from a file (optionally gzipped).
fn read_all(path: &Path, gz: bool) -> String {
    if gz {
        let file = File::open(path).expect("failed to open gz file");
        let mut decoder = MultiGzDecoder::new(file);
        let mut out = String::new();
        decoder
            .read_to_string(&mut out)
            .expect("failed to read gz contents");
        out
    } else {
        fs::read_to_string(path).expect("failed to read file")
    }
}

fn write_blocks(writer: &mut BufWriter<File>, blocks: &[Range<usize>]) {
    writer
        .write_all(&(blocks.len() as u32).to_ne_bytes())
        .expect("failed to write block count");

    let mut lengths = Vec::with_capacity(blocks.len() * 4);
    for block in blocks {
        writer
            .write_all(&(block.start as u32).to_ne_bytes())
            .expect("failed to write block start");
        lengths.extend_from_slice(&((block.end - block.start) as u32).to_ne_bytes());
    }

    writer
        .write_all(&lengths)
        .expect("failed to write block lengths");
}

fn write_packed_bases(writer: &mut BufWriter<File>, bases: &[u8]) {
    for chunk in bases.chunks(4) {
        let mut byte = 0_u8;
        for (index, base) in chunk.iter().enumerate() {
            let bits = match *base {
                b'T' | b't' | b'N' | b'n' => 0,
                b'C' | b'c' => 1,
                b'A' | b'a' => 2,
                b'G' | b'g' => 3,
                other => panic!("unsupported test nucleotide: {}", other),
            };
            byte |= bits << (6 - (index * 2));
        }
        writer
            .write_all(&[byte])
            .expect("failed to write packed 2bit bases");
    }
}
