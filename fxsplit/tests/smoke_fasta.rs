mod common;

use common::*;

#[test]
fn fasta_chunks_mean_records_per_file() {
    let t = tempdir();
    let input = join(t.path(), "input.fasta");
    let outdir = join(t.path(), "out_chunks");

    write_text(
        &input,
        ">r1\nAAAA\n>r2\nCCCC\n>r3\nGGGG\n>r4\nTTTT\n>r5\nNNNN\n",
    );

    run_split(vec![
        "--file".to_string(),
        path_str(&input),
        "--chunks".to_string(),
        "2".to_string(),
        "--outdir".to_string(),
        path_str(&outdir),
    ]);

    let files = list_files(&outdir);
    assert_eq!(files.len(), 3);

    let counts: Vec<usize> = files
        .iter()
        .map(|name| count_fasta_records(&join(&outdir, name), false))
        .collect();
    assert_eq!(counts, vec![2, 2, 1]);
}

#[test]
fn fasta_files_mean_number_of_output_files() {
    let t = tempdir();
    let input = join(t.path(), "input.fasta");
    let outdir = join(t.path(), "out_files");

    write_text(
        &input,
        ">r1\nAAAA\n>r2\nCCCC\n>r3\nGGGG\n>r4\nTTTT\n>r5\nNNNN\n",
    );

    run_split(vec![
        "--file".to_string(),
        path_str(&input),
        "--files".to_string(),
        "3".to_string(),
        "--outdir".to_string(),
        path_str(&outdir),
    ]);

    let files = list_files(&outdir);
    assert_eq!(files.len(), 3);

    let counts: Vec<usize> = files
        .iter()
        .map(|name| count_fasta_records(&join(&outdir, name), false))
        .collect();
    assert_eq!(counts, vec![2, 2, 1]);
}

#[test]
fn fasta_headers_use_sanitized_ids() {
    let t = tempdir();
    let input = join(t.path(), "input.fasta");
    let outdir = join(t.path(), "out_headers");

    write_text(&input, ">seq one\nAAAA\n>seq/two\nCCCC\n>seq/two\nGGGG\n");

    run_split(vec![
        "--file".to_string(),
        path_str(&input),
        "--headers".to_string(),
        "--outdir".to_string(),
        path_str(&outdir),
    ]);

    let files = list_files(&outdir);
    assert_eq!(
        files,
        vec![
            "seq.fasta".to_string(),
            "seq_two.fasta".to_string(),
            "seq_two_1.fasta".to_string()
        ]
    );

    for file in files {
        assert_eq!(count_fasta_records(&join(&outdir, &file), false), 1);
    }
}
