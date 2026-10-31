//! Bow-force sweep for one bowed-string performance card.
//!
//! Reads the string, damping and friction records of a
//! `frankensim-bowed-performance-v1` file (the body, ambient, listener and
//! schedule records are ignored: this observes the string and the bow with a rigid
//! bridge). For a fixed relative bow position and bow speed it runs the library's
//! own constant-gesture bowed run on a geometric grid of normal forces and
//! classifies every run with the SAME emergent detector the gates use. One JSONL
//! row per force reports the measured slip fraction beside beta, which ideal
//! Helmholtz motion predicts for it, the flyback intervals per period and the class;
//! a final row gives the playable band and Schelleng's maximum force
//! `2 Z0 v / (beta (mu_s - mu_d))` for the card's own impedance and friction drop.
//! The thresholds inside `classify` are the rig's declared tolerances, not
//! literature constants. By default the bow starts at full speed (a step onset);
//! with `--accel A` the bow speed instead ramps from rest at `A` m/s^2 and then
//! holds, refreshed on a 1 kHz control tick as the scheduled renderer does, so the
//! sweep reaches the regime a render with that onset reaches. Neither mode maps
//! hysteresis: every run starts from rest.
//!
//! With `--trace PATH` and `--points 1`, the run's bow-point string velocity and
//! string-minus-bow relative velocity are also written to PATH, one pair of
//! little-endian f64 per audio sample (bow-point first), so the onset can be read
//! period by period instead of only through the tail metrics.
//!
//! ```text
//! bowed_card_sweep CARD.performance --beta 0.0909 --speed 0.15 \
//!     [--force-min 0.02] [--force-max 3] [--points 21] [--seconds 1.0] [--accel 1.0] [--trace PATH]
//! ```

use std::process::ExitCode;

use fs_couple::bowed_string::runtime::BowedStringState;
use fs_couple::bowed_string::{
    BowGesture, BowedRunConfig, BowedRunLog, BowedStringCard, FrictionIsland, Termination,
    classify, gate_metrics, run_bowed,
};
use fs_couple::bowed_string::runtime::schedule::file::parse_friction_record;

const SCHEMA: &str = "frankensim.bowed-card-sweep.v1";
const MAX_CARD_BYTES: usize = 1024 * 1024;
const MAX_POINTS: usize = 200;
const MAX_STEPS: usize = 2_880_000;
/// Control tick of the bow gesture [Hz], the rate the example schedules declare.
const CONTROL_RATE_HZ: f64 = 1000.0;

struct Args {
    card: String,
    beta: f64,
    speed_m_s: f64,
    force_min_n: f64,
    force_max_n: f64,
    points: usize,
    seconds: f64,
    /// Constant bow acceleration from rest [m/s^2]; zero is a step onset.
    accel_m_s2: f64,
    /// Where to write the single run's velocity trace, if anywhere.
    trace: Option<String>,
}

struct Card {
    string: BowedStringCard,
    friction: fs_tribo::FrictionLaw,
    mu_static: f64,
    mu_dynamic: f64,
    subsamples: usize,
}

fn real(word: &str, what: &str) -> Result<f64, String> {
    let value: f64 = word
        .parse()
        .map_err(|_| format!("{what}: {word} is not a number"))?;
    if value.is_finite() {
        Ok(value)
    } else {
        Err(format!("{what} must be finite"))
    }
}

fn words<'a>(line: &'a str, key: &str) -> Result<Vec<&'a str>, String> {
    let mut parts = line.split_ascii_whitespace();
    if parts.next() != Some(key) {
        return Err(format!("expected a {key} record"));
    }
    Ok(parts.collect())
}

fn read_card(text: &str) -> Result<Card, String> {
    let (header, _) = text
        .split_once("\nschedule\n")
        .ok_or("expected an LF-delimited schedule record")?;
    let lines: Vec<&str> = header.lines().collect();
    if lines.len() < 6 || lines[0] != "frankensim-bowed-performance-v1" {
        return Err("expected the frankensim-bowed-performance-v1 header records".to_string());
    }
    let audio = words(lines[1], "audio")?;
    let string = words(lines[2], "string")?;
    let damping = words(lines[3], "damping")?;
    let friction = parse_friction_record(lines[4]).map_err(str::to_string)?;
    let sub = words(lines[5], "subsamples")?;
    if audio.len() != 3 || string.len() != 6 || sub.len() != 1 {
        return Err("a record has the wrong number of fields".to_string());
    }
    let sample_rate_hz: u32 = audio[0]
        .parse()
        .map_err(|_| "the audio rate is not an integer".to_string())?;
    let mode_count: usize = string[5]
        .parse()
        .map_err(|_| "the mode count is not an integer".to_string())?;
    if mode_count == 0 || mode_count > 512 || damping.len() != mode_count {
        return Err("the damping record must hold one ratio per retained mode".to_string());
    }
    let zetas = damping
        .iter()
        .map(|word| real(word, "damping"))
        .collect::<Result<Vec<f64>, String>>()?;
    let string_card = BowedStringCard {
        length_m: real(string[0], "length")?,
        tension_n: real(string[1], "tension")?,
        linear_density_kg_m: real(string[2], "linear density")?,
        bending_stiffness_n_m2: real(string[3], "bending stiffness")?,
        viscous_bending_n_m2_s: real(string[4], "viscous bending")?,
        mode_count,
        zetas,
        sample_rate_hz,
    };
    string_card
        .validate()
        .map_err(|error| format!("the string card is refused: {error}"))?;
    let (mu_static, mu_dynamic) = match friction {
        fs_tribo::FrictionLaw::Stribeck { static_mu, kinetic_mu, .. }
        | fs_tribo::FrictionLaw::Hyperbolic { static_mu, kinetic_mu, .. } => (static_mu, kinetic_mu),
        _ => return Err("unsupported shared bow friction law".to_string()),
    };
    let subsamples: usize = sub[0]
        .parse()
        .map_err(|_| "subsamples is not an integer".to_string())?;
    if subsamples == 0 || subsamples > 256 {
        return Err("subsamples must be in 1..=256".to_string());
    }
    Ok(Card {
        string: string_card,
        friction,
        mu_static,
        mu_dynamic,
        subsamples,
    })
}

fn parse_args() -> Result<Args, String> {
    let mut card = None;
    let mut beta = None;
    let mut speed = None;
    let mut force_min = 0.02;
    let mut force_max = 3.0;
    let mut points = 21_usize;
    let mut seconds = 1.0;
    let mut accel = 0.0;
    let mut trace = None;
    let mut iter = std::env::args().skip(1);
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--trace" => trace = Some(iter.next().ok_or("--trace needs a path")?),
            "--beta" | "--speed" | "--force-min" | "--force-max" | "--points" | "--seconds"
            | "--accel" => {
                let raw = iter.next().ok_or_else(|| format!("{arg} needs a value"))?;
                match arg.as_str() {
                    "--beta" => beta = Some(real(&raw, "--beta")?),
                    "--speed" => speed = Some(real(&raw, "--speed")?),
                    "--force-min" => force_min = real(&raw, "--force-min")?,
                    "--force-max" => force_max = real(&raw, "--force-max")?,
                    "--seconds" => seconds = real(&raw, "--seconds")?,
                    "--accel" => accel = real(&raw, "--accel")?,
                    _ => {
                        points = raw
                            .parse()
                            .map_err(|_| format!("--points needs an integer, got {raw}"))?;
                    }
                }
            }
            other if !other.starts_with("--") && card.is_none() => card = Some(other.to_string()),
            other => return Err(format!("unknown argument {other}")),
        }
    }
    let beta = beta.ok_or("--beta is required")?;
    let speed_m_s = speed.ok_or("--speed is required")?;
    if !(beta > 0.0 && beta < 1.0) {
        return Err("--beta must lie strictly between 0 and 1".to_string());
    }
    if !(speed_m_s > 0.0) {
        return Err("--speed must be positive".to_string());
    }
    if !(force_min > 0.0 && force_max >= force_min) {
        return Err("--force-min must be positive and not above --force-max".to_string());
    }
    if points == 0 || points > MAX_POINTS {
        return Err(format!("--points must be in 1..={MAX_POINTS}"));
    }
    if !(seconds > 0.0) {
        return Err("--seconds must be positive".to_string());
    }
    if !(accel >= 0.0) {
        return Err("--accel must be zero (a step onset) or positive".to_string());
    }
    if trace.is_some() && points != 1 {
        return Err("--trace records one run: it needs --points 1".to_string());
    }
    Ok(Args {
        card: card.ok_or("a CARD.performance path is required")?,
        beta,
        speed_m_s,
        force_min_n: force_min,
        force_max_n: force_max,
        points,
        seconds,
        accel_m_s2: accel,
        trace,
    })
}

/// The log's bow-point and relative velocities, interleaved per sample as little-endian f64.
fn trace_bytes(log: &BowedRunLog) -> Vec<u8> {
    log.bow_point_velocity_m_s
        .iter()
        .zip(&log.relative_velocity_m_s)
        .flat_map(|(bow, rel)| bow.to_le_bytes().into_iter().chain(rel.to_le_bytes()))
        .collect()
}

/// Bow speed at control tick `tick`: constant acceleration from rest, then held.
fn ramp_speed_m_s(tick: usize, accel_m_s2: f64, speed_m_s: f64) -> f64 {
    (accel_m_s2 * tick as f64 / CONTROL_RATE_HZ).min(speed_m_s)
}

/// The library's persistent runtime with the bow speed ramped at constant
/// acceleration from rest and refreshed once per control tick, the way the
/// scheduled renderer applies a gesture; force and position are constant.
/// Only the string channels the detector reads are logged.
fn run_ramped(
    config: &BowedRunConfig,
    speed_m_s: f64,
    accel_m_s2: f64,
) -> Result<BowedRunLog, String> {
    let mut state = BowedStringState::new(config, 1).map_err(|error| format!("{error}"))?;
    let samples_per_tick =
        ((f64::from(config.card.sample_rate_hz) / CONTROL_RATE_HZ).round() as usize).max(1);
    let force_n = config.gesture.normal_force_n;
    let station = config.gesture.station_fraction;
    let mut log = BowedRunLog {
        bow_point_velocity_m_s: Vec::with_capacity(config.steps),
        relative_velocity_m_s: Vec::with_capacity(config.steps),
        bridge_force_n: Vec::with_capacity(config.steps),
        body_volume_velocity_m3_s: Vec::new(),
        radiated_pressure_pa: Vec::new(),
        final_total_energy_j: 0.0,
        peak_total_energy_j: 0.0,
    };
    for i in 0..config.steps {
        if i % samples_per_tick == 0 {
            let speed = ramp_speed_m_s(i / samples_per_tick, accel_m_s2, speed_m_s);
            state
                .set_bow(speed, force_n, station)
                .map_err(|error| format!("{error}"))?;
        }
        let sample = state.step().map_err(|error| format!("{error}"))?;
        log.bow_point_velocity_m_s
            .push(sample.bow_point_velocity_m_s);
        log.relative_velocity_m_s.push(sample.relative_velocity_m_s);
        log.bridge_force_n.push(sample.bridge_force_n);
    }
    log.final_total_energy_j = state.total_modal_energy_j();
    log.peak_total_energy_j = log.final_total_energy_j;
    Ok(log)
}

/// JSON number; a non-finite value is `null`, as in the other sweep schema.
fn jnum(value: f64) -> String {
    if value.is_finite() {
        format!("{value}")
    } else {
        "null".to_string()
    }
}

fn run() -> Result<(), String> {
    let args = parse_args()?;
    let bytes =
        std::fs::read(&args.card).map_err(|error| format!("cannot read {}: {error}", args.card))?;
    if bytes.len() > MAX_CARD_BYTES {
        return Err("the card exceeds 1 MiB".to_string());
    }
    let text = std::str::from_utf8(&bytes).map_err(|_| "the card is not UTF-8".to_string())?;
    let card = read_card(text)?;
    let steps = (args.seconds * f64::from(card.string.sample_rate_hz)).round() as usize;
    if steps == 0 || steps > MAX_STEPS {
        return Err(format!("--seconds must give 1..={MAX_STEPS} samples"));
    }
    let f1 = card.string.fundamental_hz();
    let z0 = (card.string.tension_n * card.string.linear_density_kg_m).sqrt();
    let drop = card.mu_static - card.mu_dynamic;
    let schelleng_fmax_n = if drop > 0.0 {
        2.0 * z0 * args.speed_m_s / (args.beta * drop)
    } else {
        f64::NAN
    };
    let mut playable: Vec<f64> = Vec::new();
    for i in 0..args.points {
        let t = if args.points > 1 {
            i as f64 / (args.points - 1) as f64
        } else {
            0.0
        };
        let force_n = args.force_min_n * (args.force_max_n / args.force_min_n).powf(t);
        // A ramped onset starts at rest; a step onset starts at the held speed.
        let start_speed = if args.accel_m_s2 > 0.0 {
            0.0
        } else {
            args.speed_m_s
        };
        let gesture = BowGesture::admit(start_speed, force_n, args.beta)
            .map_err(|error| format!("gesture refused at {force_n} N: {error:?}"))?;
        let config = BowedRunConfig {
            card: card.string.clone(),
            island: FrictionIsland::Constitutive(card.friction.clone()),
            gesture,
            steps,
            subsamples: card.subsamples,
            termination: Termination::Rigid,
            listener_m: 1.0,
        };
        let log = if args.accel_m_s2 > 0.0 {
            run_ramped(&config, args.speed_m_s, args.accel_m_s2)
                .map_err(|error| format!("run refused at {force_n} N: {error}"))?
        } else {
            run_bowed(&config).map_err(|error| format!("run refused at {force_n} N: {error}"))?
        };
        if let Some(path) = &args.trace {
            std::fs::write(path, trace_bytes(&log))
                .map_err(|error| format!("cannot write the trace {path}: {error}"))?;
        }
        let metrics = gate_metrics(&log, &config.card);
        let class = classify(&metrics);
        if class == "playable" {
            playable.push(force_n);
        }
        println!(
            "{{\"schema\":\"{SCHEMA}\",\"card\":{:?},\"f1_hz\":{},\"beta\":{},\"speed_m_s\":{},\"accel_m_s2\":{},\"force_n\":{},\"slip_fraction\":{},\"flybacks_per_period\":{},\"peak_hz\":{},\"peak_to_semitone\":{},\"class\":\"{class}\"}}",
            args.card,
            jnum(f1),
            jnum(args.beta),
            jnum(args.speed_m_s),
            jnum(args.accel_m_s2),
            jnum(force_n),
            jnum(metrics.slip_frac),
            jnum(metrics.intervals_per_period),
            jnum(metrics.peak_hz),
            jnum(metrics.peak_to_semitone_ratio),
        );
    }
    println!(
        "{{\"schema\":\"{SCHEMA}\",\"summary\":true,\"card\":{:?},\"beta\":{},\"speed_m_s\":{},\"accel_m_s2\":{},\"steps\":{steps},\"playable_force_min_n\":{},\"playable_force_max_n\":{},\"schelleng_fmax_n\":{}}}",
        args.card,
        jnum(args.beta),
        jnum(args.speed_m_s),
        jnum(args.accel_m_s2),
        jnum(playable.first().copied().unwrap_or(f64::NAN)),
        jnum(playable.last().copied().unwrap_or(f64::NAN)),
        jnum(schelleng_fmax_n),
    );
    Ok(())
}

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(message) => {
            eprintln!("bowed_card_sweep: {message}");
            ExitCode::from(2)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const CARD: &str = "frankensim-bowed-performance-v1\naudio 48000 4801 1\nstring 0.65 60 0.0006 0 0 3\ndamping 0.001 0.001 0.001\nstribeck 0.8 0.4 0.04\nsubsamples 16\nbody 0.003 0.15 280 0.02\nambient 293.15 101325 0\nlistener_m 1\ncompile_limits 100000 1024\nschedule\nrest";

    #[test]
    fn reads_the_string_and_friction_records() {
        let card = read_card(CARD).expect("the example card must read");
        assert_eq!(card.string.mode_count, 3);
        assert_eq!(card.string.sample_rate_hz, 48_000);
        assert_eq!(card.subsamples, 16);
        assert!((card.string.length_m - 0.65).abs() < 1.0e-15);
        assert!((card.mu_static - 0.8).abs() < 1.0e-15);
        assert!((card.mu_dynamic - 0.4).abs() < 1.0e-15);
    }

    #[test]
    fn sweep_retains_the_performers_shared_hyperbolic_law() {
        let text = CARD.replace("stribeck 0.8 0.4 0.04", "friction hyperbolic 0.8 0.4 0.2");
        let card = read_card(&text).unwrap();
        assert_eq!(card.friction, fs_tribo::FrictionLaw::Hyperbolic {
            static_mu: 0.8, kinetic_mu: 0.4, characteristic_speed: 0.2,
        });
        assert_eq!((card.mu_static, card.mu_dynamic), (0.8, 0.4));
    }

    #[test]
    fn refuses_a_misordered_record_and_a_short_damping_row() {
        assert!(read_card(&CARD.replace("stribeck", "friction")).is_err());
        assert!(
            read_card(&CARD.replace("damping 0.001 0.001 0.001", "damping 0.001 0.001")).is_err()
        );
        assert!(read_card(&CARD.replace("\nschedule\nrest", "")).is_err());
    }

    #[test]
    fn the_onset_ramps_at_constant_acceleration_then_holds() {
        // 1 m/s^2 on a 1 kHz tick reaches 0.15 m/s at tick 150 and holds.
        assert_eq!(ramp_speed_m_s(0, 1.0, 0.15), 0.0);
        assert!((ramp_speed_m_s(75, 1.0, 0.15) - 0.075).abs() < 1.0e-15);
        assert!((ramp_speed_m_s(150, 1.0, 0.15) - 0.15).abs() < 1.0e-15);
        assert_eq!(ramp_speed_m_s(10_000, 1.0, 0.15), 0.15);
    }

    #[test]
    fn the_trace_interleaves_bow_point_and_relative_velocity_per_sample() {
        let log = BowedRunLog {
            bow_point_velocity_m_s: vec![0.25, -1.5],
            relative_velocity_m_s: vec![0.0, -1.65],
            bridge_force_n: vec![0.0, 0.0],
            body_volume_velocity_m3_s: Vec::new(),
            radiated_pressure_pa: Vec::new(),
            final_total_energy_j: 0.0,
            peak_total_energy_j: 0.0,
        };
        let bytes = trace_bytes(&log);
        assert_eq!(bytes.len(), 4 * 8);
        let words: Vec<f64> = bytes
            .chunks_exact(8)
            .map(|b| f64::from_le_bytes(b.try_into().unwrap()))
            .collect();
        assert_eq!(words, [0.25, 0.0, -1.5, -1.65]);
    }

    #[test]
    fn non_finite_numbers_print_as_null() {
        assert_eq!(jnum(f64::NAN), "null");
        assert_eq!(jnum(f64::INFINITY), "null");
        assert_eq!(jnum(0.5), "0.5");
    }
}
