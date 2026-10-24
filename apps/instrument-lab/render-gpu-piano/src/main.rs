//! `render-gpu-piano`: the grand_piano command line, with contact-free stretches stepped on a GPU.
//!
//! It compiles the fs-couple `grand_piano` example's own sources (nothing is copied), installs an
//! accelerator (`accel.rs`), and runs the unchanged command line: same flags, same output. See README.md.
#[allow(dead_code, unused_imports, unused_variables, unused_mut, clippy::all)]
#[path = "../../../../crates/fs-couple/examples/grand_piano/main.rs"]
mod grand_piano;
/// The grand_piano sources name their bank module from the crate root.
pub(crate) use grand_piano::linear;

mod accel;
#[cfg(all(feature = "gpu", not(target_os = "macos")))]
mod cuda;
mod emulate;
#[cfg(all(feature = "gpu", target_os = "macos"))]
mod metal;
#[cfg(feature = "gpu")]
mod pack;
mod parallel;
mod plan;

fn main() {
    grand_piano::audio::register_accelerator(accel::factory);
    let args: Vec<String> = std::env::args().skip(1).collect();
    // `--precompile sm_XY`: fill the CUDA kernel cache (NVRTC only, no GPU needed), then exit
    if args.first().map(String::as_str) == Some("--precompile") {
        #[cfg(all(feature = "gpu", not(target_os = "macos")))]
        let done = cuda::precompile(args.get(1).map_or("sm_89", String::as_str));
        #[cfg(not(all(feature = "gpu", not(target_os = "macos"))))]
        let done: Result<(), String> = Err("this build has no CUDA backend".into());
        if let Err(error) = done {
            eprintln!("render-gpu-piano: {error}");
            std::process::exit(1);
        }
        return;
    }
    let requested = std::env::var("FS_PIANO_GPU").unwrap_or_else(|_| "auto".into());
    if matches!(requested.as_str(), "cuda" | "metal")
        && args.iter().any(|arg| {
            matches!(arg.as_str(), "--diagnostic-volume" | "--bridge-trace-csv")
        })
    {
        eprintln!("render-gpu-piano: required {requested} backend cannot render CPU-only diagnostic volume or bridge tracing");
        std::process::exit(2);
    }
    if let Err(error) = grand_piano::cli_main(&args) {
        eprintln!("render-gpu-piano: {error}");
        std::process::exit(1);
    }
}
