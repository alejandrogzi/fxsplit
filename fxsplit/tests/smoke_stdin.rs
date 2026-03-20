mod common;

use common::*;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::io::Write;
use std::process::{Command, Stdio};

fn run_binary(args: &[&str], stdin_bytes: &[u8]) {
    let mut child = Command::new(env!("CARGO_BIN_EXE_fxsplit"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .spawn()
        .expect("failed to spawn fxsplit");

    child
        .stdin
        .as_mut()
        .expect("missing stdin handle")
        .write_all(stdin_bytes)
        .expect("failed to write stdin");

    let status = child.wait().expect("failed to wait on fxsplit");
    assert!(status.success());
}

#[test]
fn plain_fasta_can_be_read_from_stdin() {
    let t = tempdir();
    let outdir = join(t.path(), "out_chunks");

    run_binary(
        &[
            "--chunks",
            "2",
            "--threads",
            "2",
            "--outdir",
            &path_str(&outdir),
        ],
        b">r1\nAAAA\n>r2\nCCCC\n>r3\nGGGG\n",
    );

    let files = list_files(&outdir);
    assert_eq!(files.len(), 2);

    let counts: Vec<usize> = files
        .iter()
        .map(|name| count_fasta_records(&join(&outdir, name), false))
        .collect();
    assert_eq!(counts, vec![2, 1]);
}

#[test]
fn gz_fastq_can_be_read_from_stdin() {
    let t = tempdir();
    let outdir = join(t.path(), "out_chunks");

    let mut gz = GzEncoder::new(Vec::new(), Compression::default());
    gz.write_all(b"@r1\naaaa\n+\n!!!!\n@r2\ncccc\n+\n####\n")
        .expect("failed to write gzip payload");
    let input = gz.finish().expect("failed to finish gzip payload");

    run_binary(
        &[
            "--chunks",
            "1",
            "--threads",
            "2",
            "--outdir",
            &path_str(&outdir),
        ],
        &input,
    );

    let files = list_files(&outdir);
    assert_eq!(files.len(), 2);

    let counts: Vec<usize> = files
        .iter()
        .map(|name| count_fastq_records(&join(&outdir, name), true))
        .collect();
    assert_eq!(counts, vec![1, 1]);
}
