//! Core module for splitting a .fa/.fq file into chunks
//! Alejandro Gonzales-Irribarren, 2025
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

use crate::SplitMode;
use clap::Parser;
use std::path::PathBuf;

/// Command-line arguments for fxsplit.
///
/// # Fields
///
/// * `file` - Input FASTA/FASTQ file path (required)
/// * `chunks` - Number of records per output file
/// * `files` - Number of output files to create
/// * `headers` - Split by individual FASTA headers
/// * `threads` - Number of threads for parallel processing
/// * `outdir` - Output directory path (default: "chunks")
/// * `suffix` - Optional suffix for output filenames
///
/// # Example
///
/// ```rust, ignore
/// use std::path::PathBuf;
///
/// let args = Args::from(vec![
///     "-f".to_string(), "input.fasta.gz".to_string(),
///     "-c".to_string(), "1000".to_string(),
///     "-o".to_string(), "output".to_string(),
/// ]);
/// ```
#[derive(Debug, Parser)]
#[clap(
    name = "fxsplit",
    version = env!("CARGO_PKG_VERSION"),
    author = "Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>",
    about = "split FASTX in N chunks/files/headers"
)]
pub struct Args {
    #[arg(
        short = 'f',
        long = "file",
        required = true,
        value_name = "PATH",
        help = ".fa/.fasta/.fq/.fastq file to split (optionally .gz)"
    )]
    pub file: PathBuf,

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
        help = "Split FASTA by header (one output file per FASTA record)"
    )]
    pub headers: bool,

    #[arg(
        short = 't',
        long = "threads",
        help = "Number of threads",
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
        help = "Copy the input file instead of creating a symlink",
        action = clap::ArgAction::SetTrue
    )]
    pub copy: bool,
}

impl Args {
    /// Parses command-line arguments from a vector of strings.
    ///
    /// Prepends the program name automatically.
    ///
    /// # Arguments
    ///
    /// * `args` - Vector of command-line argument strings
    ///
    /// # Returns
    ///
    /// * `Self` - The parsed Args struct
    ///
    /// # Example
    ///
    /// ```rust, ignore
    /// let args = Args::from(vec![
    ///     "-f".to_string(), "input.fa".to_string(),
    ///     "-c".to_string(), "500".to_string(),
    /// ]);
    /// ```
    pub fn from(args: Vec<String>) -> Self {
        let mut full_args = vec![env!("CARGO_PKG_NAME").to_string()];
        full_args.extend(args);

        Args::parse_from(full_args)
    }

    /// Determines the split mode from the provided arguments.
    ///
    /// # Returns
    ///
    /// * `anyhow::Result<SplitMode>` - The determined split mode
    ///
    /// # Errors
    ///
    /// Returns an error if no split mode is specified or multiple are provided
    ///
    /// # Example
    ///
    /// ```rust, ignore
    /// let args = Args::from(vec!["-f".to_string(), "input.fa".to_string(), "-c".to_string(), "100".to_string()]);
    /// let mode = args.mode().unwrap();
    /// ```
    pub fn mode(&self) -> anyhow::Result<SplitMode> {
        match (self.chunks, self.files, self.headers) {
            (Some(n), None, false) => Ok(SplitMode::ChunkSize(n)),
            (None, Some(n), false) => Ok(SplitMode::NumFiles(n)),
            (None, None, true) => Ok(SplitMode::FileHeader),
            _ => Err(anyhow::anyhow!(
                "You must provide exactly one split mode: --chunks, --files, or --headers"
            )),
        }
    }
}
