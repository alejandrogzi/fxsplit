# fxsplit

`fxsplit` splits FASTX files (plain or gzipped) into smaller chunks/files/headers.

## Split modes

- `--chunks N`: write `N` records per output file.
- `--files K`: write exactly `K` output files with records as evenly distributed as possible.
- `--headers`: FASTA-only mode, one output file per FASTA header using sanitized FASTA IDs as filenames.

`--chunks`, `--files`, and `--headers` are mutually exclusive.

## Install

```bash
cargo install --path .
```

## Usage

```bash
fxsplit --help
fxsplit --file input.fasta --chunks 1000 --outdir chunks
fxsplit --file input.fastq.gz --files 8 --outdir parts
fxsplit --file input.fasta.gz --headers --outdir by_header
```

## Docker

Build:

```bash
docker build -t fxsplit:local .
```

Run:

```bash
docker run --rm fxsplit:local fxsplit --help
docker run --rm -v "$PWD:/data" fxsplit:local fxsplit --file /data/input.fasta --chunks 1000 --outdir /data/out
```
