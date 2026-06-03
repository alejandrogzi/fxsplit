mod common;

use common::*;

/// 25-base sequence used across the FASTA window tests.
const SEQ25: &str = "ACGTACGTACGTACGTACGTACGTA";

#[test]
fn fasta_basepairs_windows_one_file_per_window() {
    let t = tempdir();
    let input = join(t.path(), "input.fa");
    let outdir = join(t.path(), "out");

    write_text(&input, &format!(">chr1\n{SEQ25}\n"));

    run_split(vec![
        "--file".to_string(),
        path_str(&input),
        "--basepairs".to_string(),
        "10".to_string(),
        "--threads".to_string(),
        "2".to_string(),
        "--outdir".to_string(),
        path_str(&outdir),
    ]);

    let files = list_files(&outdir);
    assert_eq!(
        files,
        vec![
            "chr1_1-10.fa".to_string(),
            "chr1_11-20.fa".to_string(),
            "chr1_21-25.fa".to_string(),
        ]
    );

    // Headers keep the colon; coordinates are 1-based inclusive; tail is short.
    assert_eq!(
        read_text(&join(&outdir, "chr1_1-10.fa"), false),
        ">chr1:1-10\nACGTACGTAC\n"
    );
    assert_eq!(
        read_text(&join(&outdir, "chr1_11-20.fa"), false),
        ">chr1:11-20\nGTACGTACGT\n"
    );
    assert_eq!(
        read_text(&join(&outdir, "chr1_21-25.fa"), false),
        ">chr1:21-25\nACGTA\n"
    );

    for file in &files {
        assert_eq!(count_fasta_records(&join(&outdir, file), false), 1);
    }
}

#[test]
fn fasta_basepairs_window_ge_length_keeps_bare_name() {
    let t = tempdir();
    let input = join(t.path(), "input.fa");
    let outdir = join(t.path(), "out");

    write_text(&input, &format!(">chr1\n{SEQ25}\n"));

    run_split(vec![
        "--file".to_string(),
        path_str(&input),
        "--basepairs".to_string(),
        "1000".to_string(),
        "--outdir".to_string(),
        path_str(&outdir),
    ]);

    let files = list_files(&outdir);
    assert_eq!(files, vec!["chr1.fa".to_string()]);
    // No coordinate suffix when the window covers the whole record.
    assert_eq!(
        read_text(&join(&outdir, "chr1.fa"), false),
        format!(">chr1\n{SEQ25}\n")
    );
}

#[test]
fn fasta_basepairs_multiline_multirecord() {
    let t = tempdir();
    let input = join(t.path(), "input.fa");
    let outdir = join(t.path(), "out");

    // chrA spans two body lines (20 bp, soft-masked); chrB is 10 bp.
    write_text(&input, ">chrA\nacgtACGTAC\nGTACGTacgt\n>chrB\nNNNNNNNNNN\n");

    run_split(vec![
        "--file".to_string(),
        path_str(&input),
        "--basepairs".to_string(),
        "8".to_string(),
        "--threads".to_string(),
        "3".to_string(),
        "--outdir".to_string(),
        path_str(&outdir),
    ]);

    let files = list_files(&outdir);
    assert_eq!(
        files,
        vec![
            "chrA_1-8.fa".to_string(),
            "chrA_17-20.fa".to_string(),
            "chrA_9-16.fa".to_string(),
            "chrB_1-8.fa".to_string(),
            "chrB_9-10.fa".to_string(),
        ]
    );

    // Sequence is concatenated across body lines and original case is preserved.
    assert_eq!(
        read_text(&join(&outdir, "chrA_1-8.fa"), false),
        ">chrA:1-8\nacgtACGT\n"
    );
    assert_eq!(
        read_text(&join(&outdir, "chrA_9-16.fa"), false),
        ">chrA:9-16\nACGTACGT\n"
    );
    assert_eq!(
        read_text(&join(&outdir, "chrA_17-20.fa"), false),
        ">chrA:17-20\nacgt\n"
    );
    assert_eq!(
        read_text(&join(&outdir, "chrB_9-10.fa"), false),
        ">chrB:9-10\nNN\n"
    );
}

#[test]
fn fasta_basepairs_no_mask_uppercases() {
    let t = tempdir();
    let input = join(t.path(), "input.fa");
    let outdir = join(t.path(), "out");

    write_text(&input, ">chr1\nacgtACGTac\n");

    run_split(vec![
        "--file".to_string(),
        path_str(&input),
        "--basepairs".to_string(),
        "5".to_string(),
        "--no-mask".to_string(),
        "--outdir".to_string(),
        path_str(&outdir),
    ]);

    assert_eq!(
        read_text(&join(&outdir, "chr1_1-5.fa"), false),
        ">chr1:1-5\nACGTA\n"
    );
    assert_eq!(
        read_text(&join(&outdir, "chr1_6-10.fa"), false),
        ">chr1:6-10\nCGTAC\n"
    );
}

#[test]
fn fasta_basepairs_gzip_compression() {
    let t = tempdir();
    let input = join(t.path(), "input.fa");
    let outdir = join(t.path(), "out");

    write_text(&input, &format!(">chr1\n{SEQ25}\n"));

    run_split(vec![
        "--file".to_string(),
        path_str(&input),
        "--basepairs".to_string(),
        "10".to_string(),
        "--compression".to_string(),
        "gzip".to_string(),
        "--outdir".to_string(),
        path_str(&outdir),
    ]);

    let files = list_files(&outdir);
    assert_eq!(files.len(), 3);
    assert!(
        files.iter().all(|f| f.ends_with(".fa.gz")),
        "unexpected files: {files:?}"
    );
    assert_eq!(
        read_text(&join(&outdir, "chr1_1-10.fa.gz"), true),
        ">chr1:1-10\nACGTACGTAC\n"
    );
}

#[test]
fn fasta_basepairs_to_2bit_output() {
    let t = tempdir();
    let input = join(t.path(), "input.fa");
    let outdir = join(t.path(), "out");

    write_text(&input, &format!(">chr1\n{SEQ25}\n"));

    run_split(vec![
        "--file".to_string(),
        path_str(&input),
        "--basepairs".to_string(),
        "10".to_string(),
        "--output-format".to_string(),
        "2bit".to_string(),
        "--outdir".to_string(),
        path_str(&outdir),
    ]);

    let files = list_files(&outdir);
    assert_eq!(
        files,
        vec![
            "chr1_1-10.2bit".to_string(),
            "chr1_11-20.2bit".to_string(),
            "chr1_21-25.2bit".to_string(),
        ]
    );

    let first = read_2bit_sequences(&join(&outdir, "chr1_1-10.2bit"), true);
    assert_eq!(
        first,
        vec![("chr1:1-10".to_string(), "ACGTACGTAC".to_string())]
    );
}

#[test]
fn twobit_basepairs_windows() {
    let t = tempdir();
    let input = join(t.path(), "input.2bit");
    let outdir = join(t.path(), "out");

    // 24 bp: a multiple of 4. The `twobit` crate's `to_2bit` (used by the
    // `write_2bit` helper) corrupts the final base when a sequence length is
    // 1 mod 4, so the input fixture uses a length that round-trips cleanly.
    let seq24 = "ACGTACGTACGTACGTACGTACGT";
    write_2bit(&input, &format!(">chr1\n{seq24}\n"));

    run_split(vec![
        "--file".to_string(),
        path_str(&input),
        "--basepairs".to_string(),
        "10".to_string(),
        "--threads".to_string(),
        "2".to_string(),
        "--outdir".to_string(),
        path_str(&outdir),
    ]);

    let files = list_files(&outdir);
    assert_eq!(
        files,
        vec![
            "chr1_1-10.2bit".to_string(),
            "chr1_11-20.2bit".to_string(),
            "chr1_21-24.2bit".to_string(),
        ]
    );

    for file in &files {
        assert_eq!(count_2bit_records(&join(&outdir, file)), 1);
    }

    // First window (partial 2bit byte) and the short tail round-trip exactly.
    let first = read_2bit_sequences(&join(&outdir, "chr1_1-10.2bit"), true);
    assert_eq!(
        first,
        vec![("chr1:1-10".to_string(), "ACGTACGTAC".to_string())]
    );
    let tail = read_2bit_sequences(&join(&outdir, "chr1_21-24.2bit"), true);
    assert_eq!(tail, vec![("chr1:21-24".to_string(), "ACGT".to_string())]);
}

#[test]
fn fastq_basepairs_is_rejected() {
    let t = tempdir();
    let input = join(t.path(), "input.fastq");
    let outdir = join(t.path(), "out");

    write_text(&input, "@r1\nACGT\n+\n!!!!\n");

    let err = run_split_expect_err(vec![
        "--file".to_string(),
        path_str(&input),
        "--basepairs".to_string(),
        "2".to_string(),
        "--outdir".to_string(),
        path_str(&outdir),
    ]);

    assert!(err.contains("basepairs"), "unexpected error: {err}");
}

#[test]
fn basepairs_conflicts_with_other_modes() {
    // clap conflicts call process::exit, so drive the built binary as a subprocess.
    let t = tempdir();
    let input = join(t.path(), "input.fa");
    write_text(&input, &format!(">chr1\n{SEQ25}\n"));

    for other in [["--chunks", "2"], ["--files", "2"]] {
        let output = std::process::Command::new(env!("CARGO_BIN_EXE_fxsplit"))
            .arg("--file")
            .arg(path_str(&input))
            .arg("--basepairs")
            .arg("10")
            .args(other)
            .output()
            .expect("failed to spawn fxsplit");

        assert!(
            !output.status.success(),
            "--basepairs {other:?} should conflict but exited 0"
        );
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(
            stderr.contains("cannot be used with"),
            "unexpected stderr for {other:?}: {stderr}"
        );
    }
}
