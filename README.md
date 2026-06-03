<p align="center">
  <p align="center">
    <img width=200 align="center" src="./assets/fxsplit.png" >
  </p>

  <span>
    <h1 align="center">
        fxsplit
    </h1>
  </span>

  <p align="center">
    <a href="https://img.shields.io/badge/version-0.0.5-green" target="_blank">
      <img alt="Version Badge" src="https://img.shields.io/badge/version-0.0.5-green">
    </a>
    <a href="https://crates.io/crates/fxsplit" target="_blank">
      <img alt="Crates.io Version" src="https://img.shields.io/crates/v/fxsplit">
    </a>
    <a href="https://github.com/alejandrogzi/fxsplit" target="_blank">
      <img alt="GitHub License" src="https://img.shields.io/github/license/alejandrogzi/fxsplit?color=blue">
    </a>
    <a href="https://crates.io/crates/fxsplit" target="_blank">
      <img alt="Crates.io Total Downloads" src="https://img.shields.io/crates/d/fxsplit">
    </a>
  </p>


  <p align="center">

  </p>

  <p align="center">
    <samp>
        <span>split FASTX (and 2bit) files (plain or gzipped) into smaller chunks/files/headers</span>
        <br>
        <br>
        <a href="https://docs.rs/fxsplit/0.0.5/fxsplit/">docs</a> .
        <a href="https://github.com/alejandrogzi/fxsplit?tab=readme-ov-file#Usage">usage</a> .
        <a href="https://github.com/alejandrogzi/fxsplit?tab=readme-ov-file#Installation">install</a> .
        <a href="https://github.com/alejandrogzi/fxsplit?tab=readme-ov-file#Docker">docker</a>
    </samp>
  </p>
  
</p>

## Split modes

- `--chunks N`: write `N` records per output file.
- `--files K`: write exactly `K` output files with records as evenly distributed as possible.
- `--headers`: FASTA/2bit-only mode, one output file per header using sanitized FASTA IDs as filenames.
- `--basepairs N`: FASTA/2bit-only mode, split each record into fixed-length windows of `N` bases, one output file per window. Windows are named `{header}:{start}-{end}` using 1-based inclusive coordinates (e.g. `chr1:1-100`, `chr1:101-200`, with a shorter trailing window). When `N` is greater than or equal to a record's length, that record is emitted whole under its bare name.

`--chunks`, `--files`, `--headers`, and `--basepairs` are mutually exclusive.

## Install

```bash
cargo install fxsplit
```

## Usage

```bash
fxsplit --help
fxsplit --file input.fasta --chunks 1000 --outdir chunks
fxsplit --file input.fastq.gz --files 8 --outdir parts
fxsplit --file input.fasta.gz --headers --outdir by_header
fxsplit --file genome.fa --basepairs 1000000 --outdir windows
```

## Docker

build:

```bash
docker build -t fxsplit:local .
```

run:

```bash
docker run --rm fxsplit:local fxsplit --help
docker run --rm -v "$PWD:/data" fxsplit:local fxsplit --file /data/input.fasta --chunks 1000 --outdir /data/out
```

pull:

```bash
docker pull alejandrogzi/fxsplit:latest
```

## Nextflow

borrow fxsplit module from [fxsplit/main.nf](https://github.com/alejandrogzi/fxsplit/blob/main/assets/nextflow/fxsplit/main.nf) and use it in your pipeline.

## Conda

```bash
conda install -c bioconda fxsplit
```


