//! Schelleng-style playable-region sweep for the bowed-string fixture
//! (bead frankensim-music-v8-root-3ez8g.7.5).
//!
//! Sweeps (normal force, bow speed) on fixed bow-station grids, classifies
//! every run with the SAME emergent detector the gates use, and emits one
//! deterministic JSONL receipt row per run plus per-column boundary summary
//! rows. The measured boundary slopes are THIS RIG'S data, compared
//! qualitatively to the literature shape; no Schelleng constant is
//! transcribed as truth.
//!
//! Usage:
//! ```text
//! bowed_sweep --out PATH [--steps N] [--quick]
//! ```
//!
//! Determinism: ONE-HOST (see `fs_couple::bowed_string` module docs).

use std::io::Write;
use std::path::PathBuf;

use fs_couple::bowed_string::{
    BowGesture, BowedRunConfig, BowedStringCard, FrictionIsland, Termination, classify,
    gate_metrics, run_bowed,
};
use fs_couple::stribeck_friction::StribeckFriction;

const SCHEMA: &str = "frankensim.bowed-sweep.v1";

fn card() -> BowedStringCard {
    BowedStringCard {
        length_m: 0.65,
        tension_n: 60.0,
        linear_density_kg_m: 6.0e-4,
        mode_count: 16,
        zetas: (0..16).map(|k| 1.0e-3 + 1.5e-4 * k as f64).collect(),
        sample_rate_hz: 48_000,
    }
}

/// JSON string for an f64 that may be NaN/inf without violating JSONL
/// consumers: non-finite maps to null by convention of this schema.
fn jnum(value: f64) -> String {
    if value.is_finite() {
        format!("{value}")
    } else {
        "null".to_string()
    }
}

struct Args {
    out: PathBuf,
    steps: usize,
}

fn parse_args() -> Result<Args, String> {
    let mut out = None;
    let mut steps = 14_400_usize;
    let mut quick = false;
    let mut iter = std::env::args().skip(1);
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--out" => out = Some(iter.next().ok_or("--out needs a path")?),
            "--steps" => {
                let raw = iter.next().ok_or("--steps needs a number")?;
                steps = raw.parse().map_err(|_| format!("bad --steps {raw}"))?;
            }
            "--quick" => quick = true,
            other => return Err(format!("unknown argument {other}")),
        }
    }
    if quick {
        steps = steps.min(7_200);
    }
    Ok(Args {
        out: PathBuf::from(out.ok_or("--out PATH is required")?),
        steps,
    })
}

fn main() -> std::process::ExitCode {
    let args = match parse_args() {
        Ok(a) => a,
        Err(e) => {
            eprintln!("bowed_sweep: refusal: {e}");
            return std::process::ExitCode::from(2);
        }
    };
    if let Some(parent) = args.out.parent() {
        if let Err(e) = std::fs::create_dir_all(parent) {
            eprintln!("bowed_sweep: cannot create output dir: {e}");
            return std::process::ExitCode::from(2);
        }
    }
    let sink_file = match std::fs::File::create(&args.out) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("bowed_sweep: cannot create {}: {e}", args.out.display());
            return std::process::ExitCode::from(2);
        }
    };
    let mut sink = std::io::BufWriter::new(sink_file);

    let card = card();
    let rosin = StribeckFriction {
        mu_static: 0.8,
        mu_dynamic: 0.4,
        stiction_m_s: 0.04,
    };

    // Coarse-but-honest grid; the artifact is the regime MAP.
    let stations = [0.08_f64, 0.11];
    let forces: Vec<f64> = if args.steps <= 7_200 {
        (1..=13).map(|i| 0.3 + 0.3 * i as f64).collect()
    } else {
        (1..=26).map(|i| 0.3 + 0.15 * i as f64).collect()
    };
    let speeds: Vec<f64> = (1..=9).map(|i| 0.05 * i as f64).collect();

    let mut rows_emitted = 0_u64;
    // (station index, speed in centi-m/s) -> playable forces in that column.
    let mut playable_by_column: std::collections::BTreeMap<(usize, u64), Vec<f64>> =
        std::collections::BTreeMap::new();

    for (station_index, &station) in stations.iter().enumerate() {
        for &v_bow in &speeds {
            for &force in &forces {
                let gesture = match BowGesture::admit(v_bow, force, station) {
                    Ok(g) => g,
                    Err(_) => continue,
                };
                let cfg = BowedRunConfig {
                    card: card.clone(),
                    island: FrictionIsland::Stribeck(rosin),
                    gesture,
                    steps: args.steps,
                    subsamples: 4,
                    termination: Termination::Rigid,
                    listener_m: 1.0,
                };
                let (class, peak_hz, slip, intervals, energy, note) = match run_bowed(&cfg) {
                    Ok(log) => {
                        let m = gate_metrics(&log, &card);
                        (
                            classify(&m),
                            m.peak_hz,
                            m.slip_frac,
                            m.intervals_per_period,
                            log.peak_total_energy_j,
                            String::new(),
                        )
                    }
                    Err(e) => (
                        "refused",
                        f64::NAN,
                        f64::NAN,
                        f64::NAN,
                        f64::NAN,
                        format!("{e}"),
                    ),
                };
                if class == "playable" {
                    playable_by_column
                        .entry((station_index, (v_bow * 100.0).round() as u64))
                        .or_default()
                        .push(force);
                }
                let note_field = if note.is_empty() {
                    String::new()
                } else {
                    format!(r#","note":"{}""#, note.replace('\\', "\\\\").replace('"', "'"))
                };
                let row = format!(
                    concat!(
                        r#"{{"schema":"{SCHEMA}","kind":"run","#,
                        r#""station_fraction":{station},"v_bow_m_s":{v_bow},"normal_force_n":{force},"#,
                        r#""class":"{class}","peak_hz":{peak},"fundamental_hz":{f1},"#,
                        r#""slip_fraction":{slip},"intervals_per_period":{intervals},"#,
                        r#""peak_energy_j":{energy},"steps":{steps}{note}}}"#
                    ),
                    SCHEMA = SCHEMA,
                    station = station,
                    v_bow = v_bow,
                    force = force,
                    class = class,
                    peak = jnum(peak_hz),
                    f1 = jnum(card.fundamental_hz()),
                    slip = jnum(slip),
                    intervals = jnum(intervals),
                    energy = jnum(energy),
                    steps = args.steps,
                    note = note_field,
                );
                if writeln!(sink, "{row}").is_err() {
                    eprintln!("bowed_sweep: sink write failed");
                    return std::process::ExitCode::from(2);
                }
                rows_emitted += 1;
            }
        }
    }

    // Boundary summary: minimum and maximum PLAYABLE force per column.
    let mut boundary_rows = 0_usize;
    for ((station_index, speed_centi), mut played) in playable_by_column {
        played.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let f_min = played.first().copied().unwrap_or(f64::NAN);
        let f_max = played.last().copied().unwrap_or(f64::NAN);
        let row = format!(
            r#"{{"schema":"{SCHEMA}","kind":"boundary","station_fraction":{},"v_bow_m_s":{},"min_playable_force_n":{},"max_playable_force_n":{}}}"#,
            stations[station_index],
            speed_centi as f64 / 100.0,
            jnum(f_min),
            jnum(f_max),
        );
        let _ = writeln!(sink, "{row}");
        boundary_rows += 1;
    }

    eprintln!(
        "bowed_sweep: {} run rows + {} boundary rows -> {}",
        rows_emitted,
        boundary_rows,
        args.out.display()
    );
    std::process::ExitCode::SUCCESS
}
