//! Runnable normalized 2-D elasticity/topology optimization marquee.
//!
//! This binary deliberately exposes the real fs-topols + canonical CutFEM path.
//! It is a numerical software demonstration, not a 3-D, experiment-validated,
//! KKT-certified, or physical-design claim.

use fs_topols::{
    Cantilever, GridSdf, GuardedSettings, GuardedStop, OptimizeSettings,
    optimize_compliance_guarded,
};
use std::error::Error;
use std::fs::{File, OpenOptions};
use std::io::{BufWriter, Write};
use std::path::Path;
use std::process::ExitCode;

#[path = "elasticity/projected.rs"]
mod projected;

fn writer(path: &Path) -> std::io::Result<BufWriter<File>> {
    OpenOptions::new().write(true).create_new(true).open(path).map(BufWriter::new)
}

fn json_string(value: &str) -> String {
    let mut out = String::with_capacity(value.len() + 2);
    out.push('"');
    for ch in value.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            ch if ch.is_control() => {
                use std::fmt::Write as _;
                let _ = write!(out, "\\u{:04x}", ch as u32);
            }
            ch => out.push(ch),
        }
    }
    out.push('"');
    out
}

fn parse<T: std::str::FromStr>(args: &[String], index: usize, default: T) -> Result<T, Box<dyn Error>>
where
    T::Err: Error + 'static,
{
    match args.get(index) {
        Some(value) => Ok(value.parse()?),
        None => Ok(default),
    }
}

fn run() -> Result<u8, Box<dyn Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.first().is_some_and(|arg| arg == "--projected") {
        return projected::run(&args[1..]);
    }
    if args.is_empty() || args.len() > 5 {
        return Err("usage: fs-marquee-elasticity OUTPUT_DIR [LEVEL=4] [ITERATIONS=12] [VOLFRAC=0.45] [MAX_CANDIDATES=5]".into());
    }
    let output_dir = Path::new(&args[0]);
    if output_dir.try_exists()? {
        return Err("output directory already exists; refusing to overwrite it".into());
    }

    let level: u32 = parse(&args, 1, 4)?;
    let iterations: usize = parse(&args, 2, 12)?;
    let volfrac: f64 = parse(&args, 3, 0.45)?;
    let max_candidates: usize = parse(&args, 4, 5)?;
    if !(2..=7).contains(&level) {
        return Err("LEVEL must lie in [2, 7] for this bounded executable".into());
    }
    if iterations == 0 || iterations > 200 {
        return Err("ITERATIONS must lie in [1, 200]".into());
    }
    if !(volfrac.is_finite() && volfrac > 0.0 && volfrac <= 1.0) {
        return Err("VOLFRAC must lie in (0, 1]".into());
    }
    if !(1..=16).contains(&max_candidates) {
        return Err("MAX_CANDIDATES must lie in [1, 16]".into());
    }

    let n = 1usize << level;
    let mut phi = GridSdf::from_fn(n, &|_, y| (y - 0.5).abs() - 0.42);
    let settings = OptimizeSettings {
        level,
        volfrac,
        iterations,
        move_cells: 0.35,
        nucleation_period: 4,
        hole_radius_cells: 1.5,
        ..OptimizeSettings::default()
    };
    let guarded = GuardedSettings {
        max_candidates,
        contraction: 0.5,
        volume_tolerance: 0.01,
        min_relative_improvement: 0.0,
    };
    let report = optimize_compliance_guarded(
        &mut phi,
        Cantilever { load: 1.0, band: 0.125 },
        settings,
        guarded,
    )?;

    std::fs::create_dir(output_dir)?;
    let mut candidates = writer(&output_dir.join("candidates.jsonl"))?;
    for candidate in &report.candidates {
        let state = match candidate.final_state {
            Some(state) => format!(
                "{{\"compliance\":{:.17e},\"volume\":{:.17e},\"snapshot\":\"{:#018x}\",\"trajectory_compliance_delta\":{:.17e},\"trajectory_volume_delta\":{:.17e}}}",
                state.compliance,
                state.volume,
                state.snapshot,
                state.compliance_delta,
                state.volume_delta,
            ),
            None => "null".to_string(),
        };
        let refusal = candidate.refusal.as_deref().map_or_else(|| "null".to_string(), json_string);
        writeln!(
            candidates,
            "{{\"index\":{},\"move_cells\":{:.17e},\"volume_feasible\":{},\"improvement_gate\":{},\"final_state\":{},\"refusal\":{}}}",
            candidate.index,
            candidate.move_cells,
            candidate.volume_feasible,
            candidate.improvement_gate,
            state,
            refusal,
        )?;
    }
    candidates.flush()?;

    let (status, exit) = match report.stop {
        GuardedStop::Accepted => ("accepted", 0),
        GuardedStop::NoFeasibleCandidate => ("no_feasible_candidate", 10),
        GuardedStop::NoImprovingCandidate => ("no_improving_candidate", 11),
        GuardedStop::AllCandidatesRefused => ("all_candidates_refused", 12),
    };
    let accepted = report.accepted.map_or_else(
        || "null".to_string(),
        |state| format!(
            "{{\"compliance\":{:.17e},\"volume\":{:.17e},\"snapshot\":\"{:#018x}\",\"trajectory_compliance_delta\":{:.17e},\"trajectory_volume_delta\":{:.17e}}}",
            state.compliance,
            state.volume,
            state.snapshot,
            state.compliance_delta,
            state.volume_delta,
        ),
    );
    let summary = format!(
        "{{\"model\":\"normalized_unit_square_plane_strain_cantilever\",\"authority\":\"estimated\",\"status\":\"{status}\",\"level\":{level},\"iterations\":{iterations},\"volume_limit\":{volfrac:.17e},\"candidate_budget\":{max_candidates},\"baseline_compliance\":{:.17e},\"baseline_volume\":{:.17e},\"accepted_move_cells\":{},\"accepted\":{},\"claims\":{{\"physical_validation\":false,\"kkt_convergence\":false,\"global_optimum\":false,\"three_dimensional\":false}}}}",
        report.baseline.compliance,
        report.baseline.volume,
        report.accepted_move_cells.map_or_else(|| "null".to_string(), |v| format!("{v:.17e}")),
        accepted,
    );
    let mut summary_file = writer(&output_dir.join("summary.json"))?;
    writeln!(summary_file, "{summary}")?;
    summary_file.flush()?;

    if report.accepted.is_some() {
        let mut field = writer(&output_dir.join("level-set.csv"))?;
        writeln!(field, "x_normalized,y_normalized,phi_normalized")?;
        for j in 0..=n {
            for i in 0..=n {
                let [x, y] = phi.pos(i, j);
                writeln!(field, "{x:.17e},{y:.17e},{:.17e}", phi.node(i, j))?;
            }
        }
        field.flush()?;
        if let Some(trajectory) = report.trajectory.as_ref() {
            let mut trace = writer(&output_dir.join("trajectory.jsonl"))?;
            for row in &trajectory.rows {
                writeln!(trace, "{row}")?;
            }
            trace.flush()?;
        }
    }

    println!("{summary}");
    Ok(exit)
}

fn main() -> ExitCode {
    match run() {
        Ok(code) => ExitCode::from(code),
        Err(error) => {
            eprintln!("fs-marquee-elasticity: {error}");
            ExitCode::FAILURE
        }
    }
}
