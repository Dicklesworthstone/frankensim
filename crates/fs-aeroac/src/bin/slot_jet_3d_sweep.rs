//! `slot_jet_3d_sweep` — the recorded heavy-run driver for bead
//! frankensim-music-v8-root-3ez8g.10.1: one geometry, a ladder of
//! central-moment second-order rates at fixed `u_jet` (the clean Re
//! actuator), each rung settled/recorded/classified independently
//! with a typed per-rung refusal, all receipts to ONE fail-closed
//! JSONL file plus a terminal regime-map line.
//!
//! Receipts are measurements, not claims: a broadband rung and a
//! fully tonal ladder are both retained outcomes. The FFT bin
//! disclosure rides on every rung line.

use std::fs;
use std::path::Path;

use fs_aeroac::slot_jet_3d::{SlotJet3dConfig, classify_rung, run_slot_jet_3d};
use fs_lbm::d3q19::CollisionModel3;

struct Args {
    nx: usize,
    ny: usize,
    nz: usize,
    slot_half: f64,
    u_jet: f64,
    higher_order_rate: f64,
    second_order_rates: Vec<f64>,
    seed_amplitude: f64,
    edge_distance: usize,
    plate_length: usize,
    fringe_width: usize,
    fringe_sigma: f64,
    steps_settle: usize,
    steps_record: usize,
    out: String,
}

fn parse_args() -> Result<Args, String> {
    let mut a = Args {
        nx: 128,
        ny: 40,
        nz: 32,
        slot_half: 2.5,
        u_jet: 0.04,
        higher_order_rate: 1.9,
        second_order_rates: vec![1.9, 1.95, 1.98],
        seed_amplitude: 0.05,
        edge_distance: 24,
        plate_length: 12,
        fringe_width: 16,
        fringe_sigma: 0.5,
        steps_settle: 20000,
        steps_record: 16384,
        out: "target/slot-jet-3d-sweep".to_string(),
    };
    let mut it = std::env::args().skip(1);
    while let Some(k) = it.next() {
        let v = it.next().ok_or_else(|| format!("missing value for {k}"))?;
        match k.as_str() {
            "--nx" => a.nx = v.parse().map_err(|_| k)?,
            "--ny" => a.ny = v.parse().map_err(|_| k)?,
            "--nz" => a.nz = v.parse().map_err(|_| k)?,
            "--slot-half" => a.slot_half = v.parse().map_err(|_| k)?,
            "--u-jet" => a.u_jet = v.parse().map_err(|_| k)?,
            "--r3" => a.higher_order_rate = v.parse().map_err(|_| k)?,
            "--r2" => {
                a.second_order_rates = v
                    .split(',')
                    .map(str::parse)
                    .collect::<Result<_, _>>()
                    .map_err(|_| k.clone())?;
            }
            "--seed" => a.seed_amplitude = v.parse().map_err(|_| k)?,
            "--fringe-width" => a.fringe_width = v.parse().map_err(|_| k)?,
            "--fringe-sigma" => a.fringe_sigma = v.parse().map_err(|_| k)?,
            "--settle" => a.steps_settle = v.parse().map_err(|_| k)?,
            "--record" => a.steps_record = v.parse().map_err(|_| k)?,
            "--edge-distance" => a.edge_distance = v.parse().map_err(|_| k)?,
            "--plate-length" => a.plate_length = v.parse().map_err(|_| k)?,
            "--out" => a.out = v,
            other => return Err(format!("unknown flag {other}")),
        }
    }
    Ok(a)
}

/// The pinned per-rung geometry: one nozzle column, the authored
/// edge/plate layout, the shared fringe.
fn rung_config(args: &Args, second_order_rate: f64) -> SlotJet3dConfig {
    SlotJet3dConfig {
        nx: args.nx,
        ny: args.ny,
        nz: args.nz,
        slot_half: args.slot_half,
        u_jet: args.u_jet,
        collision: CollisionModel3::CentralMoment {
            second_order_rate,
            higher_order_rate: args.higher_order_rate,
        },
        nozzle_thickness: 1,
        edge_distance: args.edge_distance,
        plate_length: args.plate_length,
        fringe_width: args.fringe_width,
        fringe_sigma: args.fringe_sigma,
        seed_amplitude: args.seed_amplitude,
        steps_settle: args.steps_settle,
        steps_record: args.steps_record,
    }
}

/// Execute one ladder rung and render its receipt line (or typed
/// refusal line). Returns `(line, classified, broadband_hit, refusal)`.
fn execute_rung(cfg: &SlotJet3dConfig, r2: f64) -> (String, usize, usize, usize) {
    match run_slot_jet_3d(cfg) {
        Ok(run) => match classify_rung(&run, cfg) {
            Ok(rung) => {
                let broadband = usize::from(!rung.tonal);
                (rung.to_jsonl(), 1, broadband, 0)
            }
            Err(e) => (
                format!(
                    "{{\"schema\":\"fs-aeroac.slot-jet-3d.refusal/v1\",\
\"reynolds_actuator\":{r2},\"stage\":\"classify\",\"error\":\"{e}\"}}"
                ),
                0,
                0,
                1,
            ),
        },
        Err(e) => (
            format!(
                "{{\"schema\":\"fs-aeroac.slot-jet-3d.refusal/v1\",\
\"reynolds_actuator\":{r2},\"stage\":\"run\",\"error\":\"{e}\"}}"
            ),
            0,
            0,
            1,
        ),
    }
}

fn main() {
    let args = match parse_args() {
        Ok(a) => a,
        Err(e) => {
            eprintln!("slot_jet_3d_sweep argument refusal: {e}");
            eprintln!(
                "usage: slot_jet_3d_sweep [--nx N] [--ny N] [--nz N] [--slot-half F] \
[--u-jet F] [--r2 R1,R2,...] [--r3 F] [--seed F] [--fringe-width N] [--fringe-sigma F] \
[--edge-distance N] [--plate-length N] [--settle N] [--record N] [--out DIR]"
            );
            std::process::exit(2);
        }
    };
    let out_path = Path::new(&args.out).join("run.jsonl");
    if out_path.exists() {
        eprintln!(
            "refusing to overwrite existing receipt file {} (fail-closed rerun)",
            out_path.display()
        );
        std::process::exit(3);
    }
    if let Err(e) = fs::create_dir_all(&args.out) {
        eprintln!("cannot create output dir {}: {e}", args.out);
        std::process::exit(3);
    }

    let header = format!(
        "{{\"schema\":\"fs-aeroac.slot-jet-3d.sweep/v1\",\"nx\":{},\"ny\":{},\"nz\":{},\
\"slot_half\":{},\"u_jet\":{},\"higher_order_rate\":{},\"second_order_rates\":{:?},\
\"seed_amplitude\":{},\"fringe_width\":{},\"fringe_sigma\":{},\"steps_settle\":{},\
\"steps_record\":{}}}",
        args.nx,
        args.ny,
        args.nz,
        args.slot_half,
        args.u_jet,
        args.higher_order_rate,
        args.second_order_rates,
        args.seed_amplitude,
        args.fringe_width,
        args.fringe_sigma,
        args.steps_settle,
        args.steps_record
    );

    let mut lines = vec![header];
    let mut refusals = 0usize;
    let mut broadband_rungs = 0usize;
    let mut classified_rungs = 0usize;
    for &r2 in &args.second_order_rates {
        let cfg = rung_config(&args, r2);
        let (line, classified, broadband, refusal) = execute_rung(&cfg, r2);
        classified_rungs += classified;
        broadband_rungs += broadband;
        refusals += refusal;
        println!("{line}");
        lines.push(line);
    }

    // Terminal regime-map line. Both outcomes are honest wins; the
    // verdict names what was measured, nothing more.
    let verdict = if broadband_rungs > 0 {
        "broadband-rungs-present"
    } else if classified_rungs > 0 {
        "tonal-across-ladder"
    } else {
        "all-rungs-refused"
    };
    let terminal = format!(
        "{{\"schema\":\"fs-aeroac.slot-jet-3d.terminal/v1\",\"rungs\":{},\"classified\":{},\
\"refusals\":{},\"broadband_rungs\":{},\"verdict\":\"{verdict}\",\
\"no_claim\":\"no experimental or video-backed flue-noise claim; lattice measurements only\"}}",
        args.second_order_rates.len(),
        classified_rungs,
        refusals,
        broadband_rungs
    );
    println!("{terminal}");
    lines.push(terminal);

    if let Err(e) = fs::write(&out_path, lines.join("\n") + "\n") {
        eprintln!("cannot write {}: {e}", out_path.display());
        std::process::exit(3);
    }
}
