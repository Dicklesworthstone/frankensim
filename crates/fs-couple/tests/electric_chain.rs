//! Electric chain gates + the split-vs-full-DAE bake-off (music bead
//! `frankensim-music-v8-root-3ez8g.9.5`).
//!
//! The composed electric claim: a plucked MODAL STRING (the one string
//! owner) drives a Faraday PICKUP (pure reader, .9.2) through a tone
//! RC into a TRIODE gain stage (device island against the Kirchhoff
//! DAE, .9.3), whose AC plate swing drives a THIELE-SMALL sealed-box
//! driver (.9.4). Every stage is a named, receipted object; the chain
//! test asserts the COMPOSITION laws:
//!
//! - ec-001: the chain executes end-to-end with per-stage logging
//!   (levels, island iterations, plate headroom) and the
//!   NO-STRING-REDUPLICATION witness: a bare string stepped alongside
//!   stays BITWISE identical — the chain reads, never owns.
//! - ec-002: bridge-vs-neck voicing SURVIVES the chain: at linear
//!   drive, per-harmonic output ratios equal the pickup window-gain
//!   ratios (the .9.2 oracle), so the comb physics is not erased by
//!   amp or box.
//! - ec-003: the harmonic ladder through the WHOLE chain: a
//!   single-mode string (pure-sine velocity, so any harmonic content
//!   is device-made) shows H2 rising with drive and H2 > H3 at the
//!   speaker output — the .9.3 constitutive signature, composed.
//! - ec-004: the box composes transparently: the chain-measured
//!   driver transfer (excursion line / terminal-voltage line) matches
//!   the SAME driver measured alone, at a near-resonance and an
//!   off-resonance frequency.
//! - ec-005/006: the split-vs-full-DAE bake-off on a tone-RC + gain
//!   stage fixture — split Gauss-Seidel island (budget image) vs a
//!   per-sample monolithic Newton solve of the same constitutive
//!   equations to 1e-14 (authority image). Pairwise receipt, EXPECTED
//!   KeepBoth; the committed receipt's verdict is re-asserted.
//! - ec-007/008: the listening artifact (clean vs driven render of
//!   the same pluck; far-field pressure at 1 m from the driver
//!   model) with sidecar + receipt chained by digest.
//! - ec-009: the electric gate-summary artifact enumerates every
//!   `electric-string` registry row.
//!
//! Budget rows (D25) wait on the budget lane covering this chain;
//! live-default stays `no` — disclosed, not fabricated.

use std::collections::BTreeMap;

use fs_blake3::hash_domain;
use fs_couple::bakeoff::{BakeoffOutcome, BakeoffReceipt, ContenderResult};
use fs_couple::pickup::{Pickup, PickupPose};
use fs_couple::speaker::{TsCard, TsDriver};
use fs_math::det;
use fs_phs::device::{TriodeCard, TriodeStage};
use fs_psycho::receipt::{ListeningReceipt, ListeningVerdict};

const RATE: f64 = 48_000.0;
/// 400 samples per cycle EXACTLY: the 16000-sample projection
/// windows then hold an integer cycle count of every harmonic, so
/// cross-harmonic projection leakage is zero by discrete
/// orthogonality. MEASURED at 110 Hz (non-integer cycles): the
/// fundamental leaked ~1.7% into every harmonic line —
/// drive-independent — inflating the weak neck mode-3 line 27% and
/// flooring the ladder's H2 at ~0.7%.
const F0: f64 = 120.0;
/// The bake-off pluck stays at the receipt's minted frequency.
const BAKEOFF_PLUCK_HZ: f64 = 110.0;
const N_MODES: usize = 6;
const SUPPLY_V: f64 = 300.0;
const RL_OHM: f64 = 100.0e3;
const BIAS_V: f64 = -1.5;
const BOX_M3: f64 = 0.010;

fn repo_root() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(std::path::Path::parent)
        .expect("workspace root")
        .to_path_buf()
}

/// Exact-ZOH modal string rotors — the ONE string owner in these
/// fixtures (the .9.2 fixture, verbatim physics).
#[derive(Clone, PartialEq)]
struct ModalString {
    rot: Vec<(f64, f64)>,
    state: Vec<(f64, f64)>,
    omega: Vec<f64>,
}

impl ModalString {
    fn pluck(n_modes: usize, f0: f64, pluck_frac: f64) -> ModalString {
        let dt = 1.0 / RATE;
        let mut rot = Vec::new();
        let mut state = Vec::new();
        let mut omega = Vec::new();
        for k in 1..=n_modes {
            let w = core::f64::consts::TAU * f0 * k as f64;
            let zeta = 1.0e-4_f64;
            let wd = w * (1.0 - zeta * zeta).sqrt();
            let decay = det::exp(-zeta * w * dt);
            rot.push((decay * det::cos(wd * dt), decay * det::sin(wd * dt)));
            let q0 = det::sin(k as f64 * core::f64::consts::PI * pluck_frac) / (k * k) as f64;
            state.push((q0, 0.0));
            omega.push(w);
        }
        ModalString { rot, state, omega }
    }

    fn step(&mut self) -> Vec<f64> {
        for ((c, s), st) in self.rot.iter().zip(self.state.iter_mut()) {
            let (re, im) = *st;
            *st = (c * re - s * im, s * re + c * im);
        }
        self.state
            .iter()
            .zip(&self.omega)
            .map(|(&(_, im), &w)| -w * im)
            .collect()
    }
}

/// Single-frequency projection magnitude over `signal`.
fn line_mag(signal: &[f64], freq: f64) -> f64 {
    let omega = core::f64::consts::TAU * freq / RATE;
    let (mut re, mut im) = (0.0f64, 0.0f64);
    for (n, &v) in signal.iter().enumerate() {
        re += v * det::cos(omega * n as f64);
        im -= v * det::sin(omega * n as f64);
    }
    (re * re + im * im).sqrt()
}

fn rms(signal: &[f64]) -> f64 {
    (signal.iter().map(|v| v * v).sum::<f64>() / signal.len().max(1) as f64).sqrt()
}

/// One-pole tone RC (exact ZOH; both bake-off contenders and the
/// chain share this law bit-for-bit).
#[derive(Clone, Copy)]
struct ToneRc {
    a: f64,
    y: f64,
}

impl ToneRc {
    fn new(cutoff_hz: f64) -> ToneRc {
        ToneRc {
            a: det::exp(-core::f64::consts::TAU * cutoff_hz / RATE),
            y: 0.0,
        }
    }

    fn step(&mut self, x: f64) -> f64 {
        self.y = self.a * self.y + (1.0 - self.a) * x;
        self.y
    }

    /// |H| at `freq` for the exact-ZOH pole (analytic, for oracles).
    fn mag_at(&self, freq: f64) -> f64 {
        let w = core::f64::consts::TAU * freq / RATE;
        let (c, s) = (det::cos(w), det::sin(w));
        let re = 1.0 - self.a * c;
        let im = self.a * s;
        (1.0 - self.a) / (re * re + im * im).sqrt()
    }
}

/// Per-sample chain record for the logging clause.
struct ChainSample {
    emf_v: f64,
    grid_v: f64,
    plate_ac_v: f64,
    spk_u_v: f64,
    excursion_m: f64,
    pressure_1m_pa: f64,
}

/// The composed chain. The string is passed IN per step (one owner,
/// outside this struct) — the chain holds no string coordinates.
struct ElectricChain {
    pickup: Pickup,
    input_gain: f64,
    tone: ToneRc,
    stage: TriodeStage,
    out_atten: f64,
    driver: TsDriver,
}

impl ElectricChain {
    fn new(pose: PickupPose, n_modes: usize, input_gain: f64, out_atten: f64) -> ElectricChain {
        ElectricChain {
            pickup: Pickup::bind(pose, n_modes).expect("pickup"),
            input_gain,
            tone: ToneRc::new(1_200.0),
            stage: TriodeStage::new(TriodeCard::koren_12ax7(), SUPPLY_V, RL_OHM, BIAS_V)
                .expect("stage bias admission"),
            out_atten,
            driver: TsDriver::new(TsCard::datasheet_class_6p5(), Some(BOX_M3), None)
                .expect("driver"),
        }
    }

    fn step(&mut self, modal_velocities: &[f64], dt: f64) -> ChainSample {
        let emf_v = self.pickup.emf_v(modal_velocities).expect("emf");
        let grid_v = self.tone.step(self.input_gain * emf_v);
        let plate_v = self.stage.step(grid_v, dt).expect("stage step");
        let plate_ac_v = plate_v - self.stage.bias_plate_v();
        let spk_u_v = self.out_atten * plate_ac_v;
        let (_, pressure_1m_pa) = self.driver.step(spk_u_v, dt).expect("driver step");
        ChainSample {
            emf_v,
            grid_v,
            plate_ac_v,
            spk_u_v,
            excursion_m: self.driver.excursion_m(),
            pressure_1m_pa,
        }
    }
}

/// Calibrate the input gain so the tone-filtered grid drive has the
/// requested RMS for THIS string trajectory (the volume knob, set by
/// measurement — no unit folklore).
fn calibrate_input_gain(
    pose: PickupPose,
    string: &ModalString,
    n_modes: usize,
    target_grid_rms: f64,
) -> f64 {
    let pickup = Pickup::bind(pose, n_modes).expect("pickup");
    let mut probe = string.clone();
    let mut tone = ToneRc::new(1_200.0);
    let mut samples = Vec::with_capacity(4_800);
    for _ in 0..4_800 {
        let v = probe.step();
        samples.push(tone.step(pickup.emf_v(&v).expect("emf")));
    }
    target_grid_rms / rms(&samples).max(1.0e-300)
}

fn bridge_pose() -> PickupPose {
    PickupPose {
        station_frac: 0.87,
        height_m: 3.0e-3,
        aperture_frac: 0.08,
    }
}

fn neck_pose() -> PickupPose {
    PickupPose {
        station_frac: 0.30,
        height_m: 3.0e-3,
        aperture_frac: 0.08,
    }
}

#[test]
fn ec_001_chain_composes_with_one_string_owner() {
    let dt = 1.0 / RATE;
    let string0 = ModalString::pluck(N_MODES, F0, 0.22);
    let gain = calibrate_input_gain(bridge_pose(), &string0, N_MODES, 0.10);
    let mut chain = ElectricChain::new(bridge_pose(), N_MODES, gain, 0.05);
    let mut string = string0.clone();
    let mut bare = string0.clone();
    let n = 9_600usize;
    let mut emf = Vec::with_capacity(n);
    let mut grid = Vec::with_capacity(n);
    let mut plate = Vec::with_capacity(n);
    let mut spk_u = Vec::with_capacity(n);
    let mut worst_excursion = 0.0f64;
    let mut plate_min = f64::INFINITY;
    let mut plate_max = f64::NEG_INFINITY;
    let mut worst_island_residual = 0.0f64;
    for _ in 0..n {
        let v = string.step();
        bare.step();
        let s = chain.step(&v, dt);
        emf.push(s.emf_v);
        grid.push(s.grid_v);
        plate.push(s.plate_ac_v);
        spk_u.push(s.spk_u_v);
        worst_excursion = worst_excursion.max(s.excursion_m.abs());
        let plate_abs = s.plate_ac_v + chain.stage.bias_plate_v();
        plate_min = plate_min.min(plate_abs);
        plate_max = plate_max.max(plate_abs);
        worst_island_residual = worst_island_residual.max(chain.stage.telemetry().worst_residual_a);
        assert!(s.pressure_1m_pa.is_finite());
    }
    // THE ONE-OWNER WITNESS: the chain read the string; it never
    // stepped, scaled, or copied it. Bitwise.
    assert!(string == bare, "the chain mutated the string owner's state");
    let headroom_v = (400.0 - plate_max).min(plate_min);
    assert!(headroom_v > 50.0, "plate headroom {headroom_v}");
    assert!(worst_excursion < 4.0e-3, "excursion {worst_excursion}");
    assert!(rms(&emf) > 0.0 && rms(&plate) > 0.0);
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"ec-001-compose\",\"emf_rms_v\":{:.3e},\
         \"grid_rms_v\":{:.3},\"plate_ac_rms_v\":{:.3},\"spk_u_rms_v\":{:.3},\
         \"excursion_peak_m\":{:.3e},\"plate_headroom_v\":{:.1},\
         \"island_iters_per_sample\":{},\"island_worst_residual_a\":{:.3e}}}",
        rms(&emf),
        rms(&grid),
        rms(&plate),
        rms(&spk_u),
        worst_excursion,
        headroom_v,
        chain.stage.telemetry().iterations,
        worst_island_residual
    );
}

#[test]
fn ec_002_voicing_survives_the_chain() {
    // Linear drive (grid RMS 0.02 V, gain ~61 -> plate swing ~1.7 V
    // peak): per-harmonic output ratios bridge/neck must equal the
    // pickup gain ratios — amp and box are common factors that cancel
    // line-by-line on the SAME string trajectory.
    let dt = 1.0 / RATE;
    let string0 = ModalString::pluck(N_MODES, F0, 0.22);
    let gain = calibrate_input_gain(bridge_pose(), &string0, N_MODES, 0.005);
    let n = 24_000usize;
    let mut lines = Vec::new();
    let mut gains_by_pose = Vec::new();
    for pose in [bridge_pose(), neck_pose()] {
        let mut chain = ElectricChain::new(pose, N_MODES, gain, 0.05);
        gains_by_pose.push(chain.pickup.mode_gains().to_vec());
        let mut string = string0.clone();
        let mut out = Vec::with_capacity(n);
        for _ in 0..n {
            let v = string.step();
            out.push(chain.step(&v, dt).excursion_m);
        }
        let tail = &out[n / 3..];
        lines.push(
            (1..=N_MODES)
                .map(|k| line_mag(tail, F0 * k as f64))
                .collect::<Vec<_>>(),
        );
    }
    let (bridge_lines, neck_lines) = (&lines[0], &lines[1]);
    let centroid = |ls: &[f64]| -> f64 {
        let num: f64 = ls
            .iter()
            .enumerate()
            .map(|(i, l)| F0 * (i + 1) as f64 * l)
            .sum();
        num / ls.iter().sum::<f64>()
    };
    let (cb, cn) = (centroid(bridge_lines), centroid(neck_lines));
    assert!(cb > cn, "bridge centroid {cb:.0} <= neck {cn:.0}");
    // The ratio law holds where the line is STRING-carried: the gain
    // floor excludes lines parked on a station node (where nothing
    // string-borne remains to compare), and the 0.005 V RMS drive
    // keeps constitutive H2/H3 far below every kept line. The
    // integer-cycle window (see F0) is what makes the per-line
    // comparison exact — at 110 Hz the leakage of the fundamental
    // into the weak neck mode-3 line read as a drive-independent 27%
    // "violation" of a law that in fact holds.
    let floor = |gains: &[f64]| 0.10 * gains.iter().fold(0.0f64, |m, g| m.max(g.abs()));
    let (fb, fn_) = (floor(&gains_by_pose[0]), floor(&gains_by_pose[1]));
    let mut worst_rel = 0.0f64;
    let mut kept = 0usize;
    for k in 0..N_MODES {
        if gains_by_pose[0][k].abs() < fb || gains_by_pose[1][k].abs() < fn_ {
            println!(
                "{{\"suite\":\"fs-couple\",\"case\":\"ec-002-voicing\",\"skipped_line\":{k},\
                 \"reason\":\"near a station node; not string-carried\"}}"
            );
            continue;
        }
        kept += 1;
        let predicted = (gains_by_pose[0][k] / gains_by_pose[1][k]).abs();
        let measured = bridge_lines[k] / neck_lines[k].max(1.0e-300);
        let rel = (measured - predicted).abs() / predicted.max(1.0e-300);
        worst_rel = worst_rel.max(rel);
        assert!(
            rel < 0.06,
            "line {k}: chain ratio {measured:.4} vs pickup-gain ratio {predicted:.4}"
        );
    }
    assert!(kept >= 4, "too few string-carried lines kept: {kept}");
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"ec-002-voicing\",\"bridge_centroid_hz\":{cb:.0},\
         \"neck_centroid_hz\":{cn:.0},\"worst_line_ratio_rel\":{worst_rel:.4}}}"
    );
}

#[test]
fn ec_003_harmonic_ladder_through_the_chain() {
    // SINGLE-mode string: the velocity into the pickup is a pure
    // decaying sine, so every 2f0/3f0 line at the speaker output is
    // DEVICE-made (the .9.3 constitutive curvature), not string
    // content.
    let dt = 1.0 / RATE;
    let string0 = ModalString::pluck(1, F0, 0.22);
    let n = 24_000usize;
    let mut ladders = Vec::new();
    for target_grid_rms in [0.03f64, 0.40] {
        let gain = calibrate_input_gain(bridge_pose(), &string0, 1, target_grid_rms);
        let mut chain = ElectricChain::new(bridge_pose(), 1, gain, 0.02);
        let mut string = string0.clone();
        let mut out = Vec::with_capacity(n);
        for _ in 0..n {
            let v = string.step();
            out.push(chain.step(&v, dt).excursion_m);
        }
        let tail = &out[n / 3..];
        let l1 = line_mag(tail, F0);
        let h2 = line_mag(tail, 2.0 * F0) / l1.max(1.0e-300);
        let h3 = line_mag(tail, 3.0 * F0) / l1.max(1.0e-300);
        ladders.push((target_grid_rms, h2, h3));
    }
    let (lo, hi) = (ladders[0], ladders[1]);
    assert!(
        hi.1 > 4.0 * lo.1,
        "H2 must rise with drive: {:.4} -> {:.4}",
        lo.1,
        hi.1
    );
    assert!(
        hi.1 > hi.2,
        "even-first asymmetry: H2 {:.4} H3 {:.4}",
        hi.1,
        hi.2
    );
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"ec-003-ladder\",\"h2_low\":{:.5},\"h2_high\":{:.5},\
         \"h3_low\":{:.5},\"h3_high\":{:.5}}}",
        lo.1, hi.1, lo.2, hi.2
    );
}

#[test]
fn ec_004_box_transfer_composes_transparently() {
    // The chain-measured driver transfer (excursion line over
    // terminal-voltage line) must equal the driver measured ALONE at
    // the same frequency — composition adds no hidden coupling. Near
    // resonance and off resonance.
    let dt = 1.0 / RATE;
    let n = 48_000usize;
    let mut worst_rel = 0.0f64;
    let mut measured = Vec::new();
    for f in [65.0f64, 260.0] {
        // Direct: the driver alone under a 1 V sine.
        let mut driver =
            TsDriver::new(TsCard::datasheet_class_6p5(), Some(BOX_M3), None).expect("driver");
        let mut x_direct = Vec::with_capacity(n);
        for k in 0..n {
            let u = det::sin(core::f64::consts::TAU * f * k as f64 / RATE);
            driver.step(u, dt).expect("driver step");
            x_direct.push(driver.excursion_m());
        }
        let tail = n / 3..n;
        let u_direct: Vec<f64> = (0..n)
            .map(|k| det::sin(core::f64::consts::TAU * f * k as f64 / RATE))
            .collect();
        let direct = line_mag(&x_direct[tail.clone()], f) / line_mag(&u_direct[tail.clone()], f);
        // Composed: single-mode string AT this frequency through the
        // whole chain; transfer from the same line projections.
        let string0 = ModalString::pluck(1, f, 0.22);
        let gain = calibrate_input_gain(bridge_pose(), &string0, 1, 0.02);
        let mut chain = ElectricChain::new(bridge_pose(), 1, gain, 0.05);
        let mut string = string0.clone();
        let mut spk_u = Vec::with_capacity(n);
        let mut x_chain = Vec::with_capacity(n);
        for _ in 0..n {
            let v = string.step();
            let s = chain.step(&v, dt);
            spk_u.push(s.spk_u_v);
            x_chain.push(s.excursion_m);
        }
        let composed = line_mag(&x_chain[tail.clone()], f) / line_mag(&spk_u[tail], f);
        let rel = (composed - direct).abs() / direct.max(1.0e-300);
        worst_rel = worst_rel.max(rel);
        measured.push((f, direct, composed));
        assert!(
            rel < 0.05,
            "at {f} Hz: composed transfer {composed:.4e} vs direct {direct:.4e}"
        );
    }
    // The resonance region is visible: more excursion per volt near
    // resonance than well above it.
    assert!(measured[0].1 > 3.0 * measured[1].1);
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"ec-004-box\",\"xfer_65hz_m_per_v\":{:.4e},\
         \"xfer_260hz_m_per_v\":{:.4e},\"worst_composed_vs_direct_rel\":{:.4}}}",
        measured[0].1, measured[1].1, worst_rel
    );
}

// ---------------------------------------------------------------------
// The split-vs-full-DAE bake-off (tone RC + gain stage fixture).

fn bakeoff_receipt_path() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/receipts/electric-split-vs-full-dae.bakeoff")
}

/// The authority image: per-sample monolithic Newton on the SAME
/// constitutive equations, `(supply − vp)/RL = Ip(vp, vg)`, solved to
/// 1e-14 relative — no split, no sweep budget, no descriptor
/// bookkeeping.
struct FullDaeNewton {
    card: TriodeCard,
    v_p: f64,
    bias_plate_v: f64,
    iterations: usize,
}

impl FullDaeNewton {
    fn new() -> FullDaeNewton {
        let card = TriodeCard::koren_12ax7();
        let mut v_p = SUPPLY_V * 0.6;
        for _ in 0..80 {
            let (ip, dip_dvp, _) = card.plate_current(v_p, BIAS_V);
            let f = (SUPPLY_V - v_p) / RL_OHM - ip;
            if f.abs() < 1.0e-16 {
                break;
            }
            v_p = (v_p - f / (-1.0 / RL_OHM - dip_dvp)).clamp(1.0, SUPPLY_V);
        }
        FullDaeNewton {
            card,
            v_p,
            bias_plate_v: v_p,
            iterations: 0,
        }
    }

    fn step(&mut self, v_in: f64) -> f64 {
        let v_g = BIAS_V + v_in;
        for _ in 0..60 {
            self.iterations += 1;
            let (ip, dip_dvp, _) = self.card.plate_current(self.v_p.max(0.0), v_g);
            let f = (SUPPLY_V - self.v_p) / RL_OHM - ip;
            if f.abs() < 1.0e-14 * (SUPPLY_V / RL_OHM) {
                break;
            }
            self.v_p = (self.v_p - f / (-1.0 / RL_OHM - dip_dvp)).clamp(0.0, SUPPLY_V);
        }
        self.v_p - self.bias_plate_v
    }
}

/// Run one contender over a grid-drive signal; returns the AC output
/// trajectory and total solver iterations.
fn run_contender(split: bool, drive: &[f64]) -> (Vec<f64>, usize) {
    let dt = 1.0 / RATE;
    let mut tone = ToneRc::new(1_200.0);
    let mut out = Vec::with_capacity(drive.len());
    if split {
        let mut stage =
            TriodeStage::new(TriodeCard::koren_12ax7(), SUPPLY_V, RL_OHM, BIAS_V).expect("stage");
        let mut iters = 0usize;
        for &x in drive {
            let g = tone.step(x);
            let v = stage.step(g, dt).expect("split step");
            iters += stage.telemetry().iterations;
            out.push(v - stage.bias_plate_v());
        }
        (out, iters)
    } else {
        let mut full = FullDaeNewton::new();
        for &x in drive {
            let g = tone.step(x);
            out.push(full.step(g));
        }
        let iters = full.iterations;
        (out, iters)
    }
}

fn bakeoff_drive_sine(amp: f64, freq: f64, n: usize) -> Vec<f64> {
    (0..n)
        .map(|k| amp * det::sin(core::f64::consts::TAU * freq * k as f64 / RATE))
        .collect()
}

fn bakeoff_drive_pluck(amp: f64, n: usize) -> Vec<f64> {
    (0..n)
        .map(|k| {
            let t = k as f64 / RATE;
            amp * det::exp(-t / 0.15) * det::sin(core::f64::consts::TAU * BAKEOFF_PLUCK_HZ * t)
        })
        .collect()
}

fn measure_qois(split: bool) -> (BTreeMap<String, f64>, usize, Vec<f64>) {
    let n = 24_000usize;
    // Small-signal gain at 220 Hz.
    let sine = bakeoff_drive_sine(0.02, 220.0, n);
    let (out, mut iters) = run_contender(split, &sine);
    let tail = n / 3..n;
    let gain = line_mag(&out[tail.clone()], 220.0) / line_mag(&sine[tail.clone()], 220.0);
    // Harmonic ladder at drive 0.35.
    let big = bakeoff_drive_sine(0.35, 220.0, n);
    let (out2, i2) = run_contender(split, &big);
    iters += i2;
    let h2 = line_mag(&out2[tail.clone()], 440.0) / line_mag(&out2[tail], 220.0);
    // Pluck transient trajectory (kept for the cross-image L2).
    let pluck = bakeoff_drive_pluck(0.4, n);
    let (out3, i3) = run_contender(split, &pluck);
    iters += i3;
    let mut qois = BTreeMap::new();
    qois.insert("gain-small-signal-220hz".to_string(), gain);
    qois.insert("h2-rel-drive-0p35".to_string(), h2);
    (qois, iters, out3)
}

#[test]
#[ignore = "minting run: measures both images and writes the bake-off receipt"]
#[allow(clippy::too_many_lines)] // one coherent minting run
fn ec_005_mint_split_vs_full_dae_receipt() {
    let (mut split_qois, split_iters, split_pluck) = measure_qois(true);
    let (mut full_qois, full_iters, full_pluck) = measure_qois(false);
    // Cross-image pluck-transient distance (relative L2; the authority
    // is the reference trajectory, so its own value is 0 by
    // construction — disclosed in the rationale).
    let num: f64 = split_pluck
        .iter()
        .zip(&full_pluck)
        .map(|(a, b)| (a - b) * (a - b))
        .sum();
    let den: f64 = full_pluck.iter().map(|v| v * v).sum();
    let l2 = (num / den.max(1.0e-300)).sqrt();
    split_qois.insert("pluck-l2-vs-authority".to_string(), l2);
    full_qois.insert("pluck-l2-vs-authority".to_string(), 0.0);
    // Reference: analytic small-signal gain from the card's own
    // partials at the bias point, times the tone-pole magnitude;
    // authority-measured for the nonlinear QoIs (disclosed).
    let stage = TriodeStage::new(TriodeCard::koren_12ax7(), SUPPLY_V, RL_OHM, BIAS_V).expect("s");
    let (_, gp, gm) = TriodeCard::koren_12ax7().plate_current(stage.bias_plate_v(), BIAS_V);
    let analytic_gain = ToneRc::new(1_200.0).mag_at(220.0) * gm * RL_OHM / (1.0 + gp * RL_OHM);
    let mut reference = BTreeMap::new();
    reference.insert("gain-small-signal-220hz".to_string(), analytic_gain);
    reference.insert(
        "h2-rel-drive-0p35".to_string(),
        full_qois["h2-rel-drive-0p35"],
    );
    reference.insert("pluck-l2-vs-authority".to_string(), 0.0);
    let mut cards = Vec::new();
    for v in [SUPPLY_V, RL_OHM, BIAS_V, 1_200.0, 0.02, 0.35, 0.4, RATE] {
        cards.extend_from_slice(&v.to_bits().to_le_bytes());
    }
    let receipt = BakeoffReceipt {
        filling: "electric-string".to_string(),
        fixture: "crates/fs-couple/tests/electric_chain.rs tone-RC(1200 Hz) + Koren 12AX7 \
                  gain stage (300 V, 100k, -1.5 V); sine 220 Hz at 0.02/0.35 V + 110 Hz \
                  pluck burst"
            .to_string(),
        shared_cards: hash_domain("org.frankensim.fs-couple.electric-bakeoff-cards.v1", &cards),
        reference,
        contenders: [
            ContenderResult {
                image: "split-circuit-image".to_string(),
                owner_crates: vec!["fs-phs".to_string(), "fs-couple".to_string()],
                measured: split_qois.clone(),
                states: 8,
                steps: 72_000,
                solver_iterations: split_iters,
                failure_modes: vec![
                    "fixed 3-sweep Gauss-Seidel budget leaves a bounded island residual \
                     (telemetry-disclosed); QoI drift vs the authority is the measured \
                     price of the per-sample cost bound"
                        .to_string(),
                ],
            },
            ContenderResult {
                image: "full-dae-authority".to_string(),
                owner_crates: vec!["fs-phs".to_string()],
                measured: full_qois.clone(),
                states: 1,
                steps: 72_000,
                solver_iterations: full_iters,
                failure_modes: vec![
                    "unbounded per-sample Newton count (data-dependent); no real-time \
                     budget claim"
                        .to_string(),
                ],
            },
        ],
        outcome: BakeoffOutcome::KeepBoth {
            scope_a: "split-circuit-image keeps the BUDGET lane: bounded per-sample cost \
                      (fixed sweeps), disclosed residual, QoI drift within the receipt's \
                      measured envelope"
                .to_string(),
            scope_b: "full-dae-authority keeps the AUTHORITY/verification lane: per-sample \
                      1e-14 solve of the same constitutive equations; the reference other \
                      images are measured against"
                .to_string(),
        },
        rationale: format!(
            "measured: split gain {:.3} vs authority {:.3} vs analytic {:.3}; split H2 \
             {:.5} vs authority {:.5}; split pluck L2 vs authority {:.3e}; iterations \
             {} (split, fixed) vs {} (authority, data-dependent). The authority's own \
             pluck-L2 is 0 BY CONSTRUCTION (it is the reference trajectory). KeepBoth: \
             the split image buys a hard per-sample budget with a measured, small QoI \
             drift; the monolith is the authority (D14 factorization, both lanes named)",
            split_qois["gain-small-signal-220hz"],
            full_qois["gain-small-signal-220hz"],
            analytic_gain,
            split_qois["h2-rel-drive-0p35"],
            full_qois["h2-rel-drive-0p35"],
            split_qois["pluck-l2-vs-authority"],
            split_iters,
            full_iters,
        ),
        listening_receipts: Vec::new(),
    };
    receipt.validate().expect("receipt validates");
    std::fs::write(
        bakeoff_receipt_path(),
        receipt.to_canonical_bytes().expect("encode"),
    )
    .expect("write receipt");
    println!(
        "minted {} hash {}",
        bakeoff_receipt_path().display(),
        receipt.content_hash().expect("hash").to_hex()
    );
}

#[test]
fn ec_006_committed_bakeoff_receipt_holds_its_verdict() {
    let bytes = std::fs::read(bakeoff_receipt_path())
        .expect("tests/receipts/electric-split-vs-full-dae.bakeoff (mint test)");
    let receipt = BakeoffReceipt::from_canonical_bytes(&bytes).expect("decode");
    receipt.validate().expect("validate");
    let split = &receipt.contenders[0].measured;
    let full = &receipt.contenders[1].measured;
    let analytic = receipt.reference["gain-small-signal-220hz"];
    // The verdict's load-bearing facts, re-assertable from the receipt:
    // both images sit on the analytic small-signal gain; the split
    // image's nonlinear drift vs the authority stays inside the
    // authored envelope that justifies KeepBoth.
    for (name, qois) in [("split", split), ("full", full)] {
        let rel = (qois["gain-small-signal-220hz"] - analytic).abs() / analytic;
        assert!(
            rel < 0.02,
            "{name} gain off the analytic oracle by {rel:.4}"
        );
    }
    let h2_rel = (split["h2-rel-drive-0p35"] - full["h2-rel-drive-0p35"]).abs()
        / full["h2-rel-drive-0p35"].max(1.0e-300);
    assert!(h2_rel < 0.10, "split H2 drift vs authority {h2_rel:.4}");
    assert!(
        split["pluck-l2-vs-authority"] < 0.05,
        "split pluck L2 {:.3e}",
        split["pluck-l2-vs-authority"]
    );
    assert!(
        matches!(&receipt.outcome, BakeoffOutcome::KeepBoth { .. }),
        "outcome drifted: {:?}",
        receipt.outcome
    );
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"ec-006-bakeoff\",\"verdict\":\"pass\",\
         \"split_gain\":{:.3},\"full_gain\":{:.3},\"analytic_gain\":{:.3},\
         \"split_h2\":{:.5},\"full_h2\":{:.5},\"pluck_l2\":{:.3e},\"hash\":\"{}\"}}",
        split["gain-small-signal-220hz"],
        full["gain-small-signal-220hz"],
        analytic,
        split["h2-rel-drive-0p35"],
        full["h2-rel-drive-0p35"],
        split["pluck-l2-vs-authority"],
        receipt.content_hash().expect("hash").to_hex()
    );
}

// ---------------------------------------------------------------------
// The listening artifact: the same pluck, clean then driven.

fn render_chain_pressure(target_grid_rms: f64, seconds: f64) -> Vec<f64> {
    let dt = 1.0 / RATE;
    let string0 = ModalString::pluck(N_MODES, F0, 0.22);
    let gain = calibrate_input_gain(bridge_pose(), &string0, N_MODES, target_grid_rms);
    let mut chain = ElectricChain::new(bridge_pose(), N_MODES, gain, 0.05);
    let mut string = string0;
    let n = (seconds * RATE) as usize;
    let mut out = Vec::with_capacity(n);
    for _ in 0..n {
        let v = string.step();
        out.push(chain.step(&v, dt).pressure_1m_pa);
    }
    out
}

#[test]
#[ignore = "minting run: renders data/listening/electric-clean-vs-driven.{wav,provenance.json} + receipt"]
fn ec_007_mint_electric_listening_artifact() {
    let root = repo_root();
    let clean = render_chain_pressure(0.03, 1.6);
    let driven = render_chain_pressure(0.40, 1.6);
    let mut signal = clean;
    signal.extend(std::iter::repeat_n(0.0, (0.1 * RATE) as usize));
    signal.extend(driven);
    let peak = signal.iter().fold(0.0f64, |m, v| m.max(v.abs()));
    let rms_pa = rms(&signal);
    let full_scale_pa = peak * 1.25;
    let (wav, clipped) =
        fs_couple::pcm_wav::encode_pcm16_wav(&signal, 48_000, full_scale_pa).expect("wav");
    assert_eq!(clipped, 0, "never clip a listening artifact");
    let hash = fs_blake3::hash_domain("org.frankensim.music-render.wav.v1", &wav);
    std::fs::write(
        root.join("data/listening/electric-clean-vs-driven.wav"),
        &wav,
    )
    .expect("wav");
    let provenance = format!(
        "{{\"schema\":\"frankensim-music-render-provenance-v1\",\"fixture\":\"electric-clean-\
         vs-driven (same 120 Hz pluck; grid RMS 0.03 V then 0.40 V through pickup -> tone RC \
         -> 12AX7 stage -> sealed-box TS driver; signal = model far-field pressure at 1 m)\",\
         \"sample_rate_hz\":48000,\"samples\":{},\"block\":480,\
         \"full_scale_pa\":{full_scale_pa:e},\"clipped_samples\":0,\"peak_pa\":{peak:e},\
         \"rms_pa\":{rms_pa:e},\"wav_blake3\":\"{}\",\"encoder\":\"fs_couple::pcm_wav (mono \
         PCM16, never peak-normalized)\"}}\n",
        signal.len(),
        hash.to_hex()
    );
    std::fs::write(
        root.join("data/listening/electric-clean-vs-driven.provenance.json"),
        provenance,
    )
    .expect("sidecar write");
    let lat = fs_psycho::log_attack_time(&signal, 48_000.0, 480).expect("attack time");
    let receipt = ListeningReceipt {
        listener: "pending".to_string(),
        session: "2026-08-17".to_string(),
        artifact_hex: hash.to_hex(),
        artifact_ref: "data/listening/electric-clean-vs-driven.provenance.json".to_string(),
        question: "does the second (driven) pluck read as tube drive — thicker, compressed, \
                   even-harmonic — rather than fuzz or clipping artifacts?"
            .to_string(),
        verdict: ListeningVerdict::Unadjudicated,
        observations: "same pluck twice: clean (grid 0.03 V RMS) then driven (0.40 V RMS); \
                       H2-led ladder measured in ec-003; awaiting the owner's ear"
            .to_string(),
        metrics: fs_psycho::receipt::AttachedMetrics {
            loudness_sone: None,
            sharpness_acum: None,
            log_attack_time: Some(lat),
            spl_db: None,
        },
    };
    std::fs::write(
        root.join("data/listening/electric-clean-vs-driven.listening-receipt"),
        receipt.to_canonical_bytes().expect("encode"),
    )
    .expect("receipt write");
    println!(
        "minted electric-clean-vs-driven artifact, wav_blake3 {}",
        hash.to_hex()
    );
}

#[test]
fn ec_008_committed_listening_chain_holds() {
    let root = repo_root();
    let receipt_bytes =
        std::fs::read(root.join("data/listening/electric-clean-vs-driven.listening-receipt"))
            .expect("committed listening receipt (mint test)");
    let receipt = ListeningReceipt::from_canonical_bytes(&receipt_bytes).expect("decode");
    let sidecar = std::fs::read_to_string(
        root.join("data/listening/electric-clean-vs-driven.provenance.json"),
    )
    .expect("sidecar");
    let wav = std::fs::read(root.join("data/listening/electric-clean-vs-driven.wav")).expect("wav");
    let hash = fs_blake3::hash_domain("org.frankensim.music-render.wav.v1", &wav);
    assert_eq!(
        receipt.artifact_hex,
        hash.to_hex(),
        "receipt/wav digest split"
    );
    assert!(sidecar.contains(&hash.to_hex()), "sidecar/wav digest split");
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"ec-008-listening-chain\",\"artifact\":\"{}\",\
         \"verdict\":\"chain-intact\"}}",
        hash.to_hex()
    );
}

#[test]
fn ec_009_gate_summary_enumerates_electric_rows() {
    // The committed electric gate-summary artifact must name EVERY
    // electric-string row in the registry — an omitted row is a
    // silently-missing gate.
    let root = repo_root();
    let registry = std::fs::read_to_string(root.join("instrument-claims.json")).expect("registry");
    let summary = std::fs::read_to_string(root.join("data/claims/electric-gate-summary.tsv"))
        .expect("committed electric gate summary (regenerate when electric rows change)");
    assert!(summary.starts_with("# frankensim-electric-gate-summary-v1"));
    let mut images = Vec::new();
    let mut cursor = 0usize;
    while let Some(hit) = registry[cursor..].find("\"filling\": \"electric-string\"") {
        let at = cursor + hit;
        let tail = &registry[at..];
        let image_tag = "\"image\": \"";
        let img_at = tail.find(image_tag).expect("image field") + image_tag.len();
        let img_end = tail[img_at..].find('"').expect("image end");
        images.push(tail[img_at..img_at + img_end].to_string());
        cursor = at + image_tag.len();
    }
    assert!(
        images.len() >= 5,
        "expected the full electric menu, found {images:?}"
    );
    for image in &images {
        assert!(
            summary.contains(image.as_str()),
            "gate summary omits electric image {image}"
        );
    }
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"ec-009-gate-summary\",\"verdict\":\"pass\",\
         \"rows\":{}}}",
        images.len()
    );
}
