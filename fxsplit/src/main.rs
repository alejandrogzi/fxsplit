//! Core module for splitting a .fa/.fq file into chunks
//! Alejandro Gonzales-Irribarren, 2026
//!
//! This module contains the main function for splitting .fa/.fq files
//! based on custom requirements in parallel.
//!
//! In short, the module accepts any type of .fa or .fq file
//! and process the reads or sequences inside them in parallel
//! when is possible. Compressed files are also accepted. The
//! user has the ability to specify is the splitting process should
//! be done based on specific chunk sizes or number of files, and
//! the amount of parallelization that should be used in the process.

use anyhow::Result;
use clap::{self, Parser};
use fxsplit::cli::Args;
use log::{info, Level};
use simple_logger::init_with_level;

use fxsplit::*;

/// Parses command-line arguments and dispatches to the appropriate
/// splitting function based on input file type.
///
/// # Arguments
///
/// Takes command-line arguments from the shell:
///
/// * `-f, --file` - Input FASTA/FASTQ file (required)
/// * `-c, --chunks` - Records per output file
/// * `-F, --files` - Number of output files
/// * `-H, --headers` - Split by FASTA header
/// * `-t, --threads` - Number of threads (default: CPU count)
/// * `-o, --outdir` - Output directory (default: "chunks")
/// * `-s, --suffix` - Optional suffix for output files
///
/// # Returns
///
/// * `Result<()>` - Success or error
///
/// # Example
///
/// ```rust, ignore
/// // Run from command line:
/// // fxsplit -f input.fa.gz -c 1000 -o output -t 4
/// ```
fn main() -> Result<()> {
    let start = std::time::Instant::now();
    init_with_level(Level::Info).unwrap();

    let args: Args = Args::parse();

    dispatch!(&args.file, {
        "fa.gz" => split_fa_gz(&args)?,
        "fasta.gz" =>  split_fa_gz(&args)?,
        "fq.gz" =>  split_fq(&args)?,
        "fastq.gz" =>  split_fq(&args)?,
        "fq" =>  split_fq(&args)?,
        "fastq" =>  split_fq(&args)?,
        "fa" =>  split_fa(&args)?,
        "fasta" =>  split_fa(&args)?,
    });

    let elapsed = start.elapsed();
    info!("Elapsed time: {:?}", elapsed);

    Ok(())
}
