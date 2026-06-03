// Copyright (c) 2026 Alejandro Gonzalez-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

use crate::{Compression, OutputFormat, SplitMode};
use clap::Parser;
use std::path::PathBuf;

/// Command-line arguments for fxsplit.
#[derive(Debug, Parser)]
#[clap(
    name = "fxsplit",
    version = env!("CARGO_PKG_VERSION"),
    author = "Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>",
    about = "split FASTX/2BIT in N chunks/files/headers"
)]
pub struct Args {
    #[arg(
        short = 'f',
        long = "file",
        required = false,
        value_name = "PATH",
        help = ".fa/.fasta/.fq/.fastq/.2bit file to split; reads stdin if omitted or set to -"
    )]
    pub file: Option<PathBuf>,

    #[arg(
        short = 'c',
        long = "chunks",
        required = false,
        value_name = "CHUNKS",
        conflicts_with_all(["files", "headers"]),
        help = "Number of records per output file"
    )]
    pub chunks: Option<usize>,

    #[arg(
        short = 'F',
        long = "files",
        required = false,
        value_name = "FILES",
        conflicts_with_all(["chunks", "headers"]),
        help = "Number of files to split the input in"
    )]
    pub files: Option<usize>,

    #[arg(
        short = 'H',
        long = "headers",
        required = false,
        conflicts_with_all(["chunks", "files"]),
        help = "Split FASTA/2BIT by header (one output file per record)"
    )]
    pub headers: bool,

    #[arg(
        short = 'b',
        long = "basepairs",
        required = false,
        value_name = "BASEPAIRS",
        conflicts_with_all(["chunks", "files", "headers"]),
        help = "Split each FASTA/2BIT record into fixed-length base-pair windows; \
                one output file per window named {header}:{start}-{end} (1-based)"
    )]
    pub basepairs: Option<usize>,

    #[arg(
        short = 't',
        long = "threads",
        help = "Number of Rayon worker threads",
        value_name = "THREADS",
        default_value_t = num_cpus::get()
    )]
    pub threads: usize,

    #[arg(
        short = 'o',
        long = "outdir",
        required = false,
        value_name = "PATH",
        help = "Output directory path",
        default_value("chunks")
    )]
    pub outdir: PathBuf,

    #[arg(
        short = 's',
        long = "suffix",
        required = false,
        value_name = "VALUE",
        help = "Suffix to append at the end of the chunk file"
    )]
    pub suffix: Option<String>,

    #[arg(
        short = 'C',
        long = "copy",
        required = false,
        help = "Copy the input file instead of creating a symlink when a passthrough optimization is possible",
        action = clap::ArgAction::SetTrue
    )]
    pub copy: bool,

    #[arg(
        short = 'N',
        long = "no-mask",
        required = false,
        help = "Convert masked/lowercase sequence content to uppercase in the outputs",
        action = clap::ArgAction::SetTrue
    )]
    pub no_mask: bool,

    #[arg(
        short = 'Z',
        long = "output-format",
        required = false,
        value_name = "FORMAT",
        value_enum,
        help = "Output format: fasta, fastq, or 2bit (default: same as input). fastq output requires fastq input"
    )]
    pub output_format: Option<OutputFormat>,

    #[arg(
        short = 'G',
        long = "compression",
        required = false,
        value_name = "CODEC",
        value_enum,
        help = "Compress all fasta/fastq outputs with the given codec: gzip, bzip2, or zstd (default: uncompressed). 2bit outputs are never compressed"
    )]
    pub compression: Option<Compression>,
}

impl Args {
    /// Parses arguments from a vector of strings.
    ///
    /// # Arguments
    /// * `args` - Command-line arguments (excludes program name)
    pub fn from(args: Vec<String>) -> Self {
        let mut full_args = vec![env!("CARGO_PKG_NAME").to_string()];
        full_args.extend(args);
        Self::parse_from(full_args)
    }

    /// Returns the split mode based on provided arguments.
    ///
    /// # Errors
    /// Returns error if multiple or no split modes are specified
    pub fn mode(&self) -> anyhow::Result<SplitMode> {
        match (self.chunks, self.files, self.headers, self.basepairs) {
            (Some(n), None, false, None) => Ok(SplitMode::ChunkSize(n)),
            (None, Some(n), false, None) => Ok(SplitMode::NumFiles(n)),
            (None, None, true, None) => Ok(SplitMode::FileHeader),
            (None, None, false, Some(n)) => Ok(SplitMode::BasePairs(n)),
            _ => Err(anyhow::anyhow!(
                "You must provide exactly one split mode: --chunks, --files, --headers, or --basepairs"
            )),
        }
    }
}
