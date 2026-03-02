#![allow(dead_code)]

use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::fs::{self, File};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use tempfile::TempDir;

pub fn tempdir() -> TempDir {
    tempfile::tempdir().expect("failed to create tempdir")
}

pub fn write_text(path: &Path, content: &str) {
    fs::write(path, content).expect("failed to write text file");
}

pub fn write_gzip(path: &Path, content: &str) {
    let file = File::create(path).expect("failed to create gzip output");
    let mut encoder = GzEncoder::new(file, Compression::default());
    encoder
        .write_all(content.as_bytes())
        .expect("failed to write gzip content");
    encoder.finish().expect("failed to finish gzip encoder");
}

pub fn run_split(args: Vec<String>) {
    fxsplit::lib_iso_split(args).expect("split command failed");
}

pub fn run_split_expect_err(args: Vec<String>) -> String {
    let err = fxsplit::lib_iso_split(args).expect_err("split command unexpectedly succeeded");
    err.to_string()
}

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

pub fn count_fasta_records(path: &Path, gz: bool) -> usize {
    read_all(path, gz)
        .lines()
        .filter(|line| line.starts_with('>'))
        .count()
}

pub fn count_fastq_records(path: &Path, gz: bool) -> usize {
    read_all(path, gz).lines().count() / 4
}

pub fn path_str(path: &Path) -> String {
    path.to_string_lossy().to_string()
}

pub fn join(dir: &Path, name: &str) -> PathBuf {
    dir.join(name)
}

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
