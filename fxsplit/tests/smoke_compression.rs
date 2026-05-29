mod common;

use common::*;

const SAMPLE: &str = ">r1\nACGT\n>r2\nTTTT\n>r3\nGGGG\n";

fn split_with_compression(codec: &str) -> (tempfile::TempDir, std::path::PathBuf, Vec<String>) {
    let t = tempdir();
    let input = join(t.path(), "input.fa");
    let outdir = join(t.path(), "out");
    write_text(&input, SAMPLE);

    run_split(vec![
        "--file".to_string(),
        path_str(&input),
        "--chunks".to_string(),
        "1".to_string(),
        "--threads".to_string(),
        "2".to_string(),
        "--outdir".to_string(),
        path_str(&outdir),
        "--compression".to_string(),
        codec.to_string(),
    ]);

    let files = list_files(&outdir);
    (t, outdir, files)
}

#[test]
fn default_is_uncompressed() {
    let t = tempdir();
    let input = join(t.path(), "input.fa");
    let outdir = join(t.path(), "out");
    write_text(&input, SAMPLE);

    run_split(vec![
        "--file".to_string(),
        path_str(&input),
        "--chunks".to_string(),
        "1".to_string(),
        "--outdir".to_string(),
        path_str(&outdir),
    ]);

    let files = list_files(&outdir);
    assert_eq!(files.len(), 3);
    for file in &files {
        assert!(file.ends_with(".fa"), "expected plain .fa, got {file}");
    }
}

#[test]
fn gzip_compression_produces_decodable_output() {
    let (_t, outdir, files) = split_with_compression("gzip");
    assert_eq!(files.len(), 3);
    for file in &files {
        assert!(file.ends_with(".fa.gz"), "got {file}");
        assert_eq!(count_fasta_records(&join(&outdir, file), true), 1);
    }
}

#[test]
fn bzip2_compression_produces_decodable_output() {
    let (_t, outdir, files) = split_with_compression("bzip2");
    assert_eq!(files.len(), 3);
    for file in &files {
        assert!(file.ends_with(".fa.bz2"), "got {file}");
        let content = read_bzip2(&join(&outdir, file));
        assert_eq!(content.lines().filter(|l| l.starts_with('>')).count(), 1);
    }
}

#[test]
fn zstd_compression_produces_decodable_output() {
    let (_t, outdir, files) = split_with_compression("zstd");
    assert_eq!(files.len(), 3);
    for file in &files {
        assert!(file.ends_with(".fa.zst"), "got {file}");
        let content = read_zstd(&join(&outdir, file));
        assert_eq!(content.lines().filter(|l| l.starts_with('>')).count(), 1);
    }
}

#[test]
fn twobit_output_is_never_compressed() {
    let t = tempdir();
    let input = join(t.path(), "input.fa");
    let outdir = join(t.path(), "out");
    write_text(&input, SAMPLE);

    // Request zstd, but 2bit output must remain uncompressed.
    run_split(vec![
        "--file".to_string(),
        path_str(&input),
        "--chunks".to_string(),
        "3".to_string(),
        "--outdir".to_string(),
        path_str(&outdir),
        "--output-format".to_string(),
        "2bit".to_string(),
        "--compression".to_string(),
        "zstd".to_string(),
    ]);

    let files = list_files(&outdir);
    assert_eq!(files.len(), 1);
    assert!(files[0].ends_with(".2bit"), "got {}", files[0]);
    assert!(!files[0].ends_with(".zst"));
    // Still a valid, readable 2bit file.
    assert_eq!(count_2bit_records(&join(&outdir, &files[0])), 3);
}
