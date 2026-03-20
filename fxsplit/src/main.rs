use clap::Parser;
use fxsplit::{cli::Args, run};
use log::{error, info, Level};
use simple_logger::init_with_level;

/// Main entry point for the fxsplit CLI.
fn main() {
    let start = std::time::Instant::now();
    init_with_level(Level::Info).unwrap();

    let args = Args::parse();
    if let Err(err) = run(&args) {
        error!("{err}");
        for cause in err.chain().skip(1) {
            error!("caused by: {cause}");
        }
        std::process::exit(1);
    }

    info!("Elapsed time: {:?}", start.elapsed());
}
