mod common;

use common::*;

fn sample_fastq() -> &'static str {
    "@r1\nAAAA\n+\n!!!!\n@r2\nCCCC\n+\n####\n@r3\nGGGG\n+\n$$$$\n@r4\nTTTT\n+\n%%%%\n@r5\nNNNN\n+\n&&&&\n"
}

#[test]
fn fastq_plain_chunks_mean_records_per_file() {
    let t = tempdir();
    let input = join(t.path(), "input.fastq");
    let outdir = join(t.path(), "out_chunks");

    write_text(&input, sample_fastq());

    run_split(vec![
        "--file".to_string(),
        path_str(&input),
        "--chunks".to_string(),
        "2".to_string(),
        "--threads".to_string(),
        "2".to_string(),
        "--outdir".to_string(),
        path_str(&outdir),
    ]);

    let files = list_files(&outdir);
    assert_eq!(files.len(), 3);

    let counts: Vec<usize> = files
        .iter()
        .map(|name| count_fastq_records(&join(&outdir, name), false))
        .collect();
    assert_eq!(counts, vec![2, 2, 1]);
}

#[test]
fn fastq_plain_files_mean_number_of_output_files() {
    let t = tempdir();
    let input = join(t.path(), "input.fastq");
    let outdir = join(t.path(), "out_files");

    write_text(&input, sample_fastq());

    run_split(vec![
        "--file".to_string(),
        path_str(&input),
        "--files".to_string(),
        "3".to_string(),
        "--threads".to_string(),
        "2".to_string(),
        "--outdir".to_string(),
        path_str(&outdir),
    ]);

    let files = list_files(&outdir);
    assert_eq!(files.len(), 3);

    let counts: Vec<usize> = files
        .iter()
        .map(|name| count_fastq_records(&join(&outdir, name), false))
        .collect();
    assert_eq!(counts, vec![2, 2, 1]);
}

#[test]
fn fastq_plain_rejects_headers_mode() {
    let t = tempdir();
    let input = join(t.path(), "input.fastq");
    let outdir = join(t.path(), "out_headers");

    write_text(&input, sample_fastq());

    let err = run_split_expect_err(vec![
        "--file".to_string(),
        path_str(&input),
        "--headers".to_string(),
        "--threads".to_string(),
        "2".to_string(),
        "--outdir".to_string(),
        path_str(&outdir),
    ]);

    assert!(err.contains("only supported for FASTA/FASTA.GZ/2BIT"));
}
