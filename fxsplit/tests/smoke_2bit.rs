mod common;

use common::*;

fn sample_fasta() -> &'static str {
    ">r1\nAAAA\n>r2\ncccc\n>r3\nGGGG\n>r4\nTTTT\n>r5\nNNNN\n"
}

#[test]
fn twobit_chunks_mean_records_per_file() {
    let t = tempdir();
    let input = join(t.path(), "input.2bit");
    let outdir = join(t.path(), "out_chunks");

    write_2bit(&input, sample_fasta());

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
        .map(|name| count_2bit_records(&join(&outdir, name)))
        .collect();
    assert_eq!(counts, vec![2, 2, 1]);

    let first = read_2bit_sequences(&join(&outdir, &files[0]), true);
    assert_eq!(first[1], ("r2".to_string(), "cccc".to_string()));
}

#[test]
fn twobit_files_mean_number_of_output_files() {
    let t = tempdir();
    let input = join(t.path(), "input.2bit");
    let outdir = join(t.path(), "out_files");

    write_2bit(&input, sample_fasta());

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
        .map(|name| count_2bit_records(&join(&outdir, name)))
        .collect();
    assert_eq!(counts, vec![2, 2, 1]);
}

#[test]
fn twobit_headers_use_sanitized_ids() {
    let t = tempdir();
    let input = join(t.path(), "input.2bit");
    let outdir = join(t.path(), "out_headers");

    write_2bit(&input, ">seq one\nAAAA\n>seq/two\nCCCC\n>seq/two\nGGGG\n");

    run_split(vec![
        "--file".to_string(),
        path_str(&input),
        "--headers".to_string(),
        "--threads".to_string(),
        "2".to_string(),
        "--outdir".to_string(),
        path_str(&outdir),
    ]);

    let files = list_files(&outdir);
    assert_eq!(
        files,
        vec![
            "seq_one.2bit".to_string(),
            "seq_two.2bit".to_string(),
            "seq_two_1.2bit".to_string()
        ]
    );

    for file in files {
        assert_eq!(count_2bit_records(&join(&outdir, &file)), 1);
    }
}

#[test]
fn twobit_no_mask_drops_softmask() {
    let t = tempdir();
    let input = join(t.path(), "input.2bit");
    let outdir = join(t.path(), "out_chunks");

    write_2bit(&input, ">r1\nacgtACGT\n");

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

    let records = read_2bit_sequences(&join(&outdir, &files[0]), true);
    assert_eq!(records, vec![("r1".to_string(), "ACGTACGT".to_string())]);
}

#[test]
fn twobit_preserves_lowercase_n_from_overlapping_masks() {
    let t = tempdir();
    let input = join(t.path(), "input.2bit");
    let outdir = join(t.path(), "out_headers");

    write_2bit_with_masks(&input, "chr1", "ACnnGT", &[2..4], &[2..4]);

    run_split(vec![
        "--file".to_string(),
        path_str(&input),
        "--headers".to_string(),
        "--threads".to_string(),
        "2".to_string(),
        "--outdir".to_string(),
        path_str(&outdir),
    ]);

    let files = list_files(&outdir);
    assert_eq!(files, vec!["chr1.2bit".to_string()]);

    let records = read_2bit_sequences(&join(&outdir, &files[0]), true);
    assert_eq!(records, vec![("chr1".to_string(), "ACnnGT".to_string())]);
}
