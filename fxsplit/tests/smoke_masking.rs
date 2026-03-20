mod common;

use common::*;

#[test]
fn fasta_no_mask_uppercases_sequence_content() {
    let t = tempdir();
    let input = join(t.path(), "input.fasta");
    let outdir = join(t.path(), "out_chunks");

    write_text(&input, ">r1\naCgtnn\n");

    run_split(vec![
        "--file".to_string(),
        path_str(&input),
        "--chunks".to_string(),
        "1".to_string(),
        "--threads".to_string(),
        "2".to_string(),
        "--no-mask".to_string(),
        "--outdir".to_string(),
        path_str(&outdir),
    ]);

    let files = list_files(&outdir);
    assert_eq!(files.len(), 1);
    let output = read_text(&join(&outdir, &files[0]), false);
    assert_eq!(output, ">r1\nACGTNN\n");
}

#[test]
fn fastq_no_mask_uppercases_sequence_only() {
    let t = tempdir();
    let input = join(t.path(), "input.fastq");
    let outdir = join(t.path(), "out_chunks");

    write_text(&input, "@r1\naaaa\n+\nzzzz\n");

    run_split(vec![
        "--file".to_string(),
        path_str(&input),
        "--chunks".to_string(),
        "1".to_string(),
        "--threads".to_string(),
        "2".to_string(),
        "--no-mask".to_string(),
        "--outdir".to_string(),
        path_str(&outdir),
    ]);

    let files = list_files(&outdir);
    assert_eq!(files.len(), 1);
    let output = read_text(&join(&outdir, &files[0]), false);
    assert_eq!(output, "@r1\nAAAA\n+\nzzzz\n");
}
