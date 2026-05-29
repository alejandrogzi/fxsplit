mod common;

use common::*;

/// Builds a `run_split` argument vector with a trailing flag list.
fn args(input: &std::path::Path, outdir: &std::path::Path, extra: &[&str]) -> Vec<String> {
    let mut out = vec![
        "--file".to_string(),
        path_str(input),
        "--threads".to_string(),
        "2".to_string(),
        "--outdir".to_string(),
        path_str(outdir),
    ];
    out.extend(extra.iter().map(|s| s.to_string()));
    out
}

#[test]
fn fastq_input_to_fasta_output() {
    let t = tempdir();
    let input = join(t.path(), "input.fastq");
    let outdir = join(t.path(), "out");
    write_text(
        &input,
        "@r1\nAAAA\n+\n!!!!\n@r2\nCCCC\n+\n####\n@r3\nGGGG\n+\n$$$$\n",
    );

    run_split(args(
        &input,
        &outdir,
        &["--chunks", "2", "--output-format", "fasta"],
    ));

    let files = list_files(&outdir);
    assert_eq!(files.len(), 2);
    for file in &files {
        assert!(file.ends_with(".fa"), "unexpected extension: {file}");
    }
    let total: usize = files
        .iter()
        .map(|name| count_fasta_records(&join(&outdir, name), false))
        .sum();
    assert_eq!(total, 3);

    // The first record should be a faithful FASTA translation (no quality).
    let first = read_text(&join(&outdir, &files[0]), false);
    assert!(first.contains(">r1\nAAAA\n"));
    assert!(!first.contains('+'));
}

#[test]
fn fastq_input_to_2bit_output() {
    let t = tempdir();
    let input = join(t.path(), "input.fastq");
    let outdir = join(t.path(), "out");
    write_text(&input, "@r1\nACGT\n+\n!!!!\n@r2\nTTTT\n+\n####\n");

    run_split(args(
        &input,
        &outdir,
        &["--chunks", "1", "--output-format", "2bit"],
    ));

    let files = list_files(&outdir);
    assert_eq!(files.len(), 2);
    for file in &files {
        assert!(file.ends_with(".2bit"), "unexpected extension: {file}");
        assert_eq!(count_2bit_records(&join(&outdir, file)), 1);
    }
}

#[test]
fn fasta_input_to_2bit_output_roundtrips_sequence() {
    let t = tempdir();
    let input = join(t.path(), "input.fa");
    let outdir = join(t.path(), "out");
    write_text(&input, ">r1\nACGTACGT\n>r2\nTTTTGGGG\n");

    run_split(args(
        &input,
        &outdir,
        &["--chunks", "2", "--output-format", "2bit"],
    ));

    let files = list_files(&outdir);
    assert_eq!(files.len(), 1);
    let path = join(&outdir, &files[0]);
    assert!(files[0].ends_with(".2bit"));
    assert_eq!(count_2bit_records(&path), 2);

    let seqs = read_2bit_sequences(&path, false);
    assert_eq!(seqs[0], ("r1".to_string(), "ACGTACGT".to_string()));
    assert_eq!(seqs[1], ("r2".to_string(), "TTTTGGGG".to_string()));
}

#[test]
fn multiline_fasta_to_2bit_preserves_length() {
    let t = tempdir();
    let input = join(t.path(), "input.fa");
    let outdir = join(t.path(), "out");
    // 12 bases spread across three lines.
    write_text(&input, ">r1\nACGT\nACGT\nACGT\n");

    run_split(args(
        &input,
        &outdir,
        &["--chunks", "1", "--output-format", "2bit"],
    ));

    let files = list_files(&outdir);
    assert_eq!(files.len(), 1);
    let seqs = read_2bit_sequences(&join(&outdir, &files[0]), false);
    assert_eq!(seqs[0].1, "ACGTACGTACGT");
}

#[test]
fn twobit_input_to_fasta_output() {
    let t = tempdir();
    let input = join(t.path(), "input.2bit");
    let outdir = join(t.path(), "out");
    write_2bit(&input, ">chr1\nACGTACGT\n>chr2\nTTTTGGGG\n");

    run_split(args(
        &input,
        &outdir,
        &["--chunks", "2", "--output-format", "fasta"],
    ));

    let files = list_files(&outdir);
    assert_eq!(files.len(), 1);
    assert!(files[0].ends_with(".fa"));
    let content = read_text(&join(&outdir, &files[0]), false);
    assert!(content.contains(">chr1\nACGTACGT\n"));
    assert!(content.contains(">chr2\nTTTTGGGG\n"));
}

#[test]
fn twobit_input_to_fasta_preserves_softmask() {
    let t = tempdir();
    let input = join(t.path(), "input.2bit");
    let outdir = join(t.path(), "out");
    // Soft-mask bases 0..4 (lowercase), hard-mask 6..8 (N).
    write_2bit_with_masks(&input, "chr1", "ACGTACGT", &[6..8], &[0..4]);

    run_split(args(
        &input,
        &outdir,
        &["--chunks", "1", "--output-format", "fasta"],
    ));

    let files = list_files(&outdir);
    let content = read_text(&join(&outdir, &files[0]), false);
    // Lowercase soft-mask retained, hard-mask shown as N.
    assert!(content.contains(">chr1\nacgtACNN\n"), "got: {content:?}");
}

#[test]
fn twobit_input_to_fasta_no_mask_uppercases_softmask() {
    let t = tempdir();
    let input = join(t.path(), "input.2bit");
    let outdir = join(t.path(), "out");
    write_2bit_with_masks(&input, "chr1", "ACGTACGT", &[6..8], &[0..4]);

    run_split(args(
        &input,
        &outdir,
        &["--chunks", "1", "--output-format", "fasta", "--no-mask"],
    ));

    let files = list_files(&outdir);
    let content = read_text(&join(&outdir, &files[0]), false);
    // Soft-mask uppercased, hard-mask N preserved.
    assert!(content.contains(">chr1\nACGTACNN\n"), "got: {content:?}");
}

#[test]
fn fasta_input_to_2bit_no_mask_drops_softmask() {
    let t = tempdir();
    let input = join(t.path(), "input.fa");
    let outdir = join(t.path(), "out");
    write_text(&input, ">r1\nacgtACGT\n");

    run_split(args(
        &input,
        &outdir,
        &["--chunks", "1", "--output-format", "2bit", "--no-mask"],
    ));

    let files = list_files(&outdir);
    // With --no-mask, softmask is dropped: reading with mask still yields uppercase.
    let seqs = read_2bit_sequences(&join(&outdir, &files[0]), true);
    assert_eq!(seqs[0].1, "ACGTACGT");
}

#[test]
fn fasta_input_rejects_fastq_output() {
    let t = tempdir();
    let input = join(t.path(), "input.fa");
    let outdir = join(t.path(), "out");
    write_text(&input, ">r1\nACGT\n");

    let err = run_split_expect_err(args(
        &input,
        &outdir,
        &["--chunks", "1", "--output-format", "fastq"],
    ));
    assert!(
        err.contains("cannot produce FASTQ from FASTA"),
        "got: {err}"
    );
}

#[test]
fn twobit_input_rejects_fastq_output() {
    let t = tempdir();
    let input = join(t.path(), "input.2bit");
    let outdir = join(t.path(), "out");
    write_2bit(&input, ">chr1\nACGT\n");

    let err = run_split_expect_err(args(
        &input,
        &outdir,
        &["--chunks", "1", "--output-format", "fastq"],
    ));
    assert!(err.contains("cannot produce FASTQ from 2bit"), "got: {err}");
}

#[test]
fn fasta_to_2bit_rejects_iupac_ambiguity_codes() {
    let t = tempdir();
    let input = join(t.path(), "input.fa");
    let outdir = join(t.path(), "out");
    write_text(&input, ">r1\nACGTR\n");

    let err = run_split_expect_err(args(
        &input,
        &outdir,
        &["--chunks", "1", "--output-format", "2bit"],
    ));
    assert!(err.contains("not representable in 2bit"), "got: {err}");
}
