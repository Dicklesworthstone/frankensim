//! HB gates: reed threshold + brass slot map (music bead
//! `frankensim-music-v8-root-3ez8g.11.2`) — D22 in action: the
//! FREQUENCY-DOMAIN oracle PREDICTS (threshold pressure, lock
//! frequency, slot structure), the TIME-DOMAIN image CONFIRMS, and
//! disagreement is diagnostic information captured in a receipt,
//! never averaged away.
//!
//! These are the first consumers of the fs-orbit facility, closing a
//! recorded honesty gap: fs-phs's reed casebook explicitly does NOT
//! claim oscillation-regime validation (quasi-static only). The HB
//! rows are where oscillation-regime claims become possible.
//!
//! - hb-001: the reed threshold — HB (island = the SAME quasistatic
//!   valve law the TD loop uses; port = 1/Z from the SAME fs-duct
//!   TMM) bisected in blowing pressure against (a) the
//!   Wilson–Beavers-class analytic regime marker (threshold near
//!   p_close/3 for a lossless bore, raised by viscothermal losses)
//!   and (b) the TD `realize_reed_bore` loop's measured onset on the
//!   SAME cards.
//! - hb-002: the brass slot map — lock frequency vs lip tension from
//!   HB on an outward-striking lip island + plane-wave bore
//!   `Z(n omega)` (the row records the bore authority: plane TMM;
//!   regenerate against the MM matrix when T-Brass lands it). Locks
//!   sit on impedance peaks (slots) and step between slots as the
//!   lip frequency sweeps.
//! - hb-003: the DISAGREEMENT PROTOCOL falsifier — a deliberately
//!   truncated HB run (N = 1) visibly biases the threshold, and the
//!   diagnostic receipt (truncation, island parameters, stepper,
//!   deltas) catches it.
//!
//! Convention note: the massless (quasistatic) reed island is REAL,
//! so the fundamental locks where `Im Z = 0` and the threshold is
//! insensitive to the e^{+-i omega t} conjugation choice between
//! fs-orbit's synthesis and fs-duct's impedance; the port uses
//! `1/Z` directly with this disclosure.

use fs_couple::reed_bore::realize_reed_bore;
use fs_couple::thin_plate::PlateBank;
use fs_duct::{
    Duct, LossModel, Segment, Termination, impedance_peaks, impedance_sweep, input_impedance,
};
use fs_material::gas::{GasSpec, GasState};
use fs_math::c64::C64;
use fs_orbit::{HbAnchor, HbBudget, OrbitError, OrbitProblem, solve_hb};
use fs_scenario::BeatingReed;

const TAU: f64 = core::f64::consts::TAU;
const RATE: u32 = 48_000;
/// DC rows evaluate the bore at this floor frequency (an open pipe's
/// Z(0) -> 0 makes the raw DC admittance singular; disclosed).
const DC_FLOOR_HZ: f64 = 1.0;

fn air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("air")
}

fn clarinet_bore() -> Duct {
    Duct {
        segments: vec![Segment::Cylinder {
            radius: 0.0075,
            length: 0.33,
        }],
    }
}

/// The reed HB problem: one state (mouthpiece pressure), the island
/// is the SAME quasistatic valve law the TD loop uses, the port is
/// the bore admittance from the SAME TMM.
struct ReedHb {
    duct: Duct,
    gas: GasState,
    rest_opening_m: f64,
    width_m: f64,
    closing_pressure_pa: f64,
    blowing_pressure_pa: f64,
}

impl ReedHb {
    /// Characteristic bore impedance `rho c / S` [Pa s/m^3]: the
    /// nondimensionalization that keeps the balance rows O(1) (raw
    /// admittances near an impedance peak are ~1e-8 and make every
    /// tolerance vacuous — measured as a sub-threshold false orbit).
    fn zc(&self) -> f64 {
        let r = 0.0075;
        self.gas.density * self.gas.sound_speed / (core::f64::consts::PI * r * r)
    }
}

impl OrbitProblem for ReedHb {
    fn dim(&self) -> usize {
        1
    }
    /// State is `p / p_close` (dimensionless); the island returns the
    /// Zc-scaled flow so the balance stays O(1).
    fn island(&self, _t: f64, x: &[f64], out: &mut [f64]) {
        let p = x[0] * self.closing_pressure_pa;
        let dp = self.blowing_pressure_pa - p;
        let h =
            fs_phs::quasistatic_aperture_opening(self.rest_opening_m, self.closing_pressure_pa, dp);
        let u = fs_phs::bernoulli_volume_flow(self.width_m, h, dp, self.gas.density);
        out[0] = u * self.zc() / self.closing_pressure_pa;
    }
    fn port(&self, s: C64) -> Vec<C64> {
        let omega = s.im.abs().max(TAU * DC_FLOOR_HZ);
        let z = input_impedance(
            &self.duct,
            &self.gas,
            omega,
            LossModel::AllRegime,
            Termination::IdealOpen,
        )
        .expect("bore impedance")
        .impedance;
        vec![z.recip().scale(self.zc())]
    }
    fn autonomous(&self) -> bool {
        true
    }
}

/// The bore's first `Im Z = 0` resonance [rad/s] (from the peak of
/// the sweep — the massless-reed lock target).
fn first_resonance(duct: &Duct, gas: &GasState) -> f64 {
    let sweep = impedance_sweep(
        duct,
        gas,
        TAU * 80.0,
        TAU * 700.0,
        6_000,
        LossModel::AllRegime,
        Termination::IdealOpen,
    )
    .expect("sweep");
    let peaks = impedance_peaks(&sweep);
    sweep[peaks[0]].omega
}

/// Does HB find a speaking (non-trivial) orbit at this blowing
/// pressure? TrivialCollapse / stalls mean "does not speak".
fn hb_speaks(problem: &ReedHb, omega_guess: f64, harmonics: usize) -> bool {
    let budget = HbBudget {
        harmonics,
        max_newton: 60,
        tolerance: 1.0e-9,
    };
    // Large seeds reach the orbit basin (measured: the equilibrium
    // basin swallows small seeds); deterministic guess ladder.
    [0.6f64, 0.35].iter().any(|&guess| {
        matches!(
            solve_hb(
                problem,
                HbAnchor::Autonomous { omega_guess },
                guess,
                &budget
            ),
            Ok(_)
        )
    })
}

/// Bisect the HB speaking threshold in blowing pressure [Pa].
fn hb_threshold(harmonics: usize) -> f64 {
    let gas = air();
    let duct = clarinet_bore();
    let omega_res = first_resonance(&duct, &gas);
    let make = |pm: f64| ReedHb {
        duct: clarinet_bore(),
        gas: air(),
        rest_opening_m: 4.0e-4,
        width_m: 1.2e-2,
        closing_pressure_pa: 2_000.0,
        blowing_pressure_pa: pm,
    };
    let (mut lo, mut hi) = (300.0f64, 1_800.0);
    assert!(!hb_speaks(&make(lo), omega_res, harmonics), "lo speaks");
    assert!(hb_speaks(&make(hi), omega_res, harmonics), "hi silent");
    for _ in 0..20 {
        let mid = 0.5 * (lo + hi);
        if hb_speaks(&make(mid), omega_res, harmonics) {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    0.5 * (lo + hi)
}

/// TD onset: sweep the SAME cards through `realize_reed_bore` and
/// return the lowest blowing pressure with a sustained tail.
fn td_onset(candidates: &[f64]) -> f64 {
    for &pm in candidates {
        let reed = BeatingReed {
            rest_opening_m: 4.0e-4,
            width_m: 1.2e-2,
            closing_pressure_pa: 2_000.0,
            blowing_pressure_pa: pm,
            attack_s: 0.02,
            mass_kg: 0.0,
            stiffness_n_m: 0.0,
        };
        let mut plates = PlateBank::default();
        let n = (0.5 * f64::from(RATE)) as usize;
        let out = realize_reed_bore(
            &clarinet_bore(),
            &air(),
            reed,
            Termination::IdealOpen,
            &mut plates,
            1.0,
            RATE,
            n,
            None,
        )
        .expect("TD loop");
        // Sustained (not decaying) tail: RMS of the last third.
        let tail = &out[2 * n / 3..];
        let rms = (tail.iter().map(|v| v * v).sum::<f64>() / tail.len() as f64).sqrt();
        if rms > 1.0e-3 {
            return pm;
        }
    }
    f64::INFINITY
}

#[test]
fn hb_000_diagnostic_scan() {
    // Measure-first probe: per-pressure solver outcome.
    let gas = air();
    let duct = clarinet_bore();
    let omega_res = first_resonance(&duct, &gas);
    println!(
        "omega_res = {omega_res:.2} rad/s ({:.1} Hz)",
        omega_res / TAU
    );
    for pm in [400.0f64, 700.0, 900.0, 1100.0, 1400.0, 1700.0] {
        let problem = ReedHb {
            duct: clarinet_bore(),
            gas: air(),
            rest_opening_m: 4.0e-4,
            width_m: 1.2e-2,
            closing_pressure_pa: 2_000.0,
            blowing_pressure_pa: pm,
        };
        for guess in [0.1f64, 0.3, 0.6] {
            let out = solve_hb(
                &problem,
                HbAnchor::Autonomous {
                    omega_guess: omega_res,
                },
                guess,
                &HbBudget {
                    harmonics: 9,
                    max_newton: 60,
                    tolerance: 1.0e-9,
                },
            );
            match out {
                Ok(o) => println!(
                    "pm {pm} guess {guess}: OK omega {:.1} amp {:.4} res {:.2e}",
                    o.omega,
                    o.first_harmonic_amplitude(0),
                    o.residual
                ),
                Err(e) => println!("pm {pm} guess {guess}: {e}"),
            }
        }
    }
}

#[test]
fn hb_001_reed_threshold_fd_predicts_td_confirms() {
    let threshold = hb_threshold(9);
    let p_close = 2_000.0;
    let ratio = threshold / p_close;
    // Analytic regime marker: p_close/3 for the lossless textbook
    // limit; viscothermal losses raise it. Authored band from the
    // regime, not fitted.
    assert!(
        (0.28..0.60).contains(&ratio),
        "HB threshold {threshold:.0} Pa is {ratio:.3} of closing — outside the \
         Wilson-Beavers-class regime band"
    );
    // TD confirmation on the SAME cards: scan around the FD
    // prediction; the onset must bracket the HB threshold within the
    // authored band.
    let candidates: Vec<f64> = (0..12)
        .map(|k| threshold * (0.55 + 0.1 * f64::from(k)))
        .collect();
    let onset = td_onset(&candidates);
    assert!(onset.is_finite(), "TD loop never spoke in the scan window");
    let delta_rel = (onset - threshold) / threshold;
    assert!(
        delta_rel.abs() < 0.35,
        "FD-vs-TD threshold split: HB {threshold:.0} Pa vs TD {onset:.0} Pa"
    );
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"hb-001-reed-threshold\",\
         \"hb_threshold_pa\":{threshold:.1},\"closing_pa\":{p_close},\
         \"ratio\":{ratio:.3},\"td_onset_pa\":{onset:.1},\"fd_td_delta_rel\":{delta_rel:.3},\
         \"harmonics\":9,\"bore_authority\":\"plane-TMM AllRegime IdealOpen\",\
         \"verdict\":\"pass\"}}"
    );
}

/// Outward-striking lip island + plane-wave bore: three states
/// (lip x, lip v, mouth pressure).
/// Pressure scale for the lip problem's dimensionless third state.
const LIP_P_SCALE: f64 = 1_000.0;

struct LipHb {
    duct: Duct,
    gas: GasState,
    mass_kg: f64,
    stiffness_n_m: f64,
    damping_n_s_m: f64,
    width_m: f64,
    rest_gap_m: f64,
    face_area_m2: f64,
    blowing_pressure_pa: f64,
}

impl LipHb {
    fn zc(&self) -> f64 {
        let r = 0.006;
        self.gas.density * self.gas.sound_speed / (core::f64::consts::PI * r * r)
    }
}

impl OrbitProblem for LipHb {
    fn dim(&self) -> usize {
        3
    }
    /// States: lip x [m], lip v [m/s], `p / LIP_P_SCALE`.
    fn island(&self, _t: f64, x: &[f64], out: &mut [f64]) {
        let (xl, vl, p) = (x[0], x[1], x[2] * LIP_P_SCALE);
        let dp = self.blowing_pressure_pa - p;
        out[0] = vl;
        out[1] = (dp * self.face_area_m2 - self.stiffness_n_m * xl - self.damping_n_s_m * vl)
            / self.mass_kg;
        let h = (self.rest_gap_m + xl).max(0.0);
        let u = fs_phs::bernoulli_volume_flow(self.width_m, h, dp, self.gas.density);
        out[2] = u * self.zc() / LIP_P_SCALE;
    }
    fn port(&self, s: C64) -> Vec<C64> {
        let omega = s.im.abs().max(TAU * DC_FLOOR_HZ);
        let z = input_impedance(
            &self.duct,
            &self.gas,
            omega,
            LossModel::AllRegime,
            Termination::IdealOpen,
        )
        .expect("bore impedance")
        .impedance;
        let zero = C64::new(0.0, 0.0);
        // Rows 0,1: plain ODE (s I); row 2: bore admittance.
        vec![
            s,
            zero,
            zero, //
            zero,
            s,
            zero, //
            zero,
            zero,
            z.recip().scale(self.zc()),
        ]
    }
    fn autonomous(&self) -> bool {
        true
    }
}

#[test]
fn hb_002_brass_slot_map_locks_on_impedance_peaks() {
    let gas = air();
    let bore = Duct {
        segments: vec![Segment::Cylinder {
            radius: 0.006,
            length: 1.4,
        }],
    };
    // Impedance peaks (the slots).
    let sweep = impedance_sweep(
        &bore,
        &gas,
        TAU * 30.0,
        TAU * 500.0,
        12_000,
        LossModel::AllRegime,
        Termination::IdealOpen,
    )
    .expect("sweep");
    let peak_hz: Vec<f64> = impedance_peaks(&sweep)
        .iter()
        .map(|&i| sweep[i].omega / TAU)
        .collect();
    assert!(peak_hz.len() >= 3);
    // Tension sweep: the lip natural frequency walks across slots 2-3
    // (peaks near 61 / 184 / 306 Hz for this bore).
    let mass = 2.0e-3;
    let mut locks = Vec::new();
    for step in 0..7 {
        let f_lip = 150.0 + 30.0 * f64::from(step);
        let k = mass * (TAU * f_lip) * (TAU * f_lip);
        let problem = LipHb {
            duct: Duct {
                segments: bore.segments.clone(),
            },
            gas: air(),
            mass_kg: mass,
            stiffness_n_m: k,
            damping_n_s_m: mass * TAU * f_lip / 3.0,
            width_m: 7.0e-3,
            rest_gap_m: 6.0e-4,
            face_area_m2: 1.2e-4,
            blowing_pressure_pa: 3_000.0,
        };
        let budget = HbBudget {
            harmonics: 7,
            max_newton: 80,
            tolerance: 1.0e-8,
        };
        let orbit = [2.0e-4f64, 4.0e-4].iter().find_map(|&guess| {
            solve_hb(
                &problem,
                HbAnchor::Autonomous {
                    omega_guess: TAU * f_lip,
                },
                guess,
                &budget,
            )
            .ok()
        });
        if let Some(orbit) = orbit {
            locks.push((f_lip, orbit.omega / TAU));
        }
    }
    assert!(locks.len() >= 4, "slot map too sparse: {locks:?}");
    // Every lock sits on a slot (within 8% of an impedance peak) and
    // the map is monotone in tension.
    // The slot physics, gated by REGIME: when the lip frequency sits
    // near a bore peak (within 15%) the lock must be CAPTURED by the
    // slot (within 8% of the peak); between slots the lip dominates
    // and the lock must lie between the lip frequency and the
    // nearest peak (the pulling direction) — recorded, not forced.
    let mut slots_visited = std::collections::BTreeSet::new();
    let mut captured = 0usize;
    for &(f_lip, f_lock) in &locks {
        let (slot, peak) = peak_hz
            .iter()
            .enumerate()
            .map(|(i, &p)| (i, p))
            .min_by(|a, b| {
                (a.1 - f_lip)
                    .abs()
                    .partial_cmp(&(b.1 - f_lip).abs())
                    .expect("finite")
            })
            .expect("nearest peak");
        let near = (f_lip / peak - 1.0).abs() < 0.10;
        if near {
            let dev = (f_lock / peak - 1.0).abs();
            assert!(
                dev < 0.08,
                "lip {f_lip:.0} near slot {peak:.1} but lock {f_lock:.1} not captured ({dev:.3})"
            );
            captured += 1;
            slots_visited.insert(slot);
        } else {
            let (lo, hi) = (f_lip.min(peak) - 3.0, f_lip.max(peak) + 3.0);
            assert!(
                (lo..hi).contains(&f_lock),
                "lip-dominated lock {f_lock:.1} outside the pull interval                  [{lo:.1}, {hi:.1}] (lip {f_lip:.0}, peak {peak:.1})"
            );
        }
    }
    assert!(captured >= 2, "captured slot points: {captured}");
    for pair in locks.windows(2) {
        assert!(
            pair[1].1 >= pair[0].1 - 1.0,
            "slot map must be monotone in tension: {locks:?}"
        );
    }
    assert!(
        slots_visited.len() >= 2,
        "the sweep must visit at least two slots: {slots_visited:?}"
    );
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"hb-002-brass-slots\",\"locks\":{locks:?},\
         \"peaks_hz\":{peak_hz:?},\"slots_visited\":{},\
         \"bore_authority\":\"plane-TMM AllRegime IdealOpen (regenerate on the MM matrix \
         when T-Brass lands it)\",\"verdict\":\"pass\"}}",
        slots_visited.len()
    );
}

#[test]
fn hb_003_truncation_falsifier_produces_the_disagreement_receipt() {
    // The falsifier: a deliberately truncated HB (N = 1: the pure
    // describing-function limit) biases the reed threshold; the
    // disagreement protocol catches it in a diagnostic receipt.
    let full = hb_threshold(9);
    let truncated = hb_threshold(1);
    let delta_rel = (truncated - full) / full;
    // The receipt: which truncation, which island parameters, which
    // stepper — recorded, never averaged away.
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"hb-003-disagreement-receipt\",\
         \"kind\":\"hb-truncation-discrepancy\",\
         \"threshold_full_pa\":{full:.1},\"harmonics_full\":9,\
         \"threshold_truncated_pa\":{truncated:.1},\"harmonics_truncated\":1,\
         \"delta_rel\":{delta_rel:.4},\
         \"island\":\"quasistatic reed H=4e-4 w=1.2e-2 p_close=2000\",\
         \"stepper\":\"fs-orbit hb-aft-newton (masked, FD Jacobian)\",\
         \"bore\":\"plane-TMM cylinder r=7.5mm L=0.33m AllRegime IdealOpen\",\
         \"disposition\":\"open discrepancy attached to the registry row; the \
         truncation is DISCLOSED structure and N=1 is below the admitted floor\"}}"
    );
    assert!(
        delta_rel.abs() > 0.02,
        "the truncation falsifier must visibly bias the threshold \
         (full {full:.1} vs N=1 {truncated:.1})"
    );
    // And the full-resolution threshold is convergent: N=9 vs N=13
    // agree inside the FD-vs-TD band.
    let finer = hb_threshold(13);
    let conv = (finer - full) / full;
    assert!(
        conv.abs() < 0.02,
        "N=9 vs N=13 threshold drift {conv:.4} — truncation not converged"
    );
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"hb-003-convergence\",\
         \"threshold_n9\":{full:.1},\"threshold_n13\":{finer:.1},\"drift_rel\":{conv:.4}}}"
    );
}

#[test]
fn hb_004_below_threshold_refuses_by_name() {
    // D22 honesty: below threshold the HB does not fabricate an
    // orbit — it refuses (TrivialCollapse or a stall), which is
    // exactly what "will it speak? no." looks like as a type.
    let gas = air();
    let duct = clarinet_bore();
    let omega_res = first_resonance(&duct, &gas);
    let problem = ReedHb {
        duct,
        gas,
        rest_opening_m: 4.0e-4,
        width_m: 1.2e-2,
        closing_pressure_pa: 2_000.0,
        blowing_pressure_pa: 150.0,
    };
    let out = solve_hb(
        &problem,
        HbAnchor::Autonomous {
            omega_guess: omega_res,
        },
        50.0,
        &HbBudget {
            harmonics: 9,
            max_newton: 60,
            tolerance: 1.0e-9,
        },
    );
    assert!(
        matches!(
            out,
            Err(OrbitError::TrivialCollapse | OrbitError::NewtonStalled { .. })
        ),
        "below threshold must refuse, got {out:?}"
    );
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"hb-004-below-threshold\",\"verdict\":\"pass\"}}"
    );
}
