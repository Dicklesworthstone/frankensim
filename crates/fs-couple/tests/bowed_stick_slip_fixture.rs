//! M4.3-W2 (IRP pexfp.152): the ideal-string stick-slip fixture and the
//! time-step ladder of the set-valued bow contact.
//!
//! **String.** The exact-ZOH modal string of `bowed_string` with no bending.
//! Its supports are rigid, and a resistive bridge enters as modal damping. A
//! dashpot `R` reflects velocity waves with `r = (Z0 - R)/(Z0 + R)`, so every
//! mode decays at `sigma = (c/2L) ln((R + Z0)/(R - Z0))` and
//! `zeta_k = sigma/omega_k`. Spatially this is uniform viscous damping: the
//! decay rates are the dashpot's, but the bridge itself does not move.
//!
//! **Start.** Schelleng's bounds are bounds on sustaining Helmholtz motion. Each
//! run therefore starts from ideal Helmholtz motion at mid-stick, restored
//! through the modal model's transactional seam, rather than from rest (from
//! rest the string mostly settles into multi-slip motion). In mass-normalized
//! modal coordinates that state is `q_k = 0` and
//! `v_k = sqrt(mu L/2) (Y omega_1/k) (-1)^(k+1)`, with `Y = 2 v_b/(beta pi omega_1)`.
//! It is the Fourier series of the duty-beta bow-point velocity.
//!
//! **Bounds.** The analytic bow-force bounds are Schelleng's (1973), as restated
//! on Euphonics section 9.3.1:
//! - Eq. (4): `F_max = 2 Z0 v_b/(beta (mu_s - mu_d))`.
//! - Eq. (9): `F_min = Z0^2 v_b/(2 R beta^2 (mu_s - mu_d))`.
//!
//! The bead's mistranscribed `F_max` is printed beside them, per the verifier's
//! ruling. Eq. (4) takes the slip speed as `v_b/beta`. Keeping the exact
//! `v_b (1 - beta)/beta` is printed as a diagnosis, not gated.
//!
//! **Observables.** All come from the contact receipts.

use fs_couple::bowed_string::runtime::{BowContactReceipt, BowedStringState};
use fs_couple::bowed_string::{
    BowGesture, BowedRunConfig, BowedStringCard, FrictionIsland, Termination,
};
use fs_couple::modal_acoustic_time::ModalAcousticState;
use fs_couple::stribeck_friction::{ContactBranch, StribeckFriction};
use fs_math::det;
use std::sync::OnceLock;

// String fixture values follow the cello D string of the Galluzzo bowing rig
// (Euphonics 9.6.3 Table 2: 146.8 Hz, 111 N, 2.7 g/m). They are used here as
// fixture values for an ideal string.
const F0_HZ: f64 = 146.8;
const TENSION_N: f64 = 111.0;
const DENSITY_KG_M: f64 = 2.7e-3;
/// Bridge resistance over string impedance; Schelleng's analysis needs
/// `R >> Z0`, and beta >= Z0/(4R) for a Helmholtz region at beta = 1/20.
const R_OVER_Z0: f64 = 20.0;
const BOW_SPEED_M_S: f64 = 0.1;
/// Static and asymptotic kinetic limits of Smith & Woodhouse (2000) Eq. (2).
const MU_S: f64 = 1.2;
const MU_D: f64 = 0.35;
/// Stribeck decay width, a fixture value well below every Helmholtz slip
/// speed `v_b (1/beta - 1)`, so the slip coefficient is mu_d.
const DECAY_M_S: f64 = 0.05;
const BETAS: [f64; 5] = [1.0 / 20.0, 1.0 / 14.0, 1.0 / 10.0, 1.0 / 7.0, 1.0 / 5.0];
/// Contact-rate ladder. The retained band scales with the rate (64 modes at
/// 48 kHz), so `omega_K h` and hence the discrete bow-point mobility stay
/// fixed while both refine towards the ideal string.
const RUNGS: [(u32, usize); 3] = [(48_000, 64), (96_000, 128), (192_000, 256)];
/// Periods run from the Helmholtz start, and the classified tail.
const RUN_PERIODS: f64 = 60.0;
const TAIL_PERIODS: usize = 20;

fn z0() -> f64 {
    det::sqrt(TENSION_N * DENSITY_KG_M)
}
fn wave_speed() -> f64 {
    det::sqrt(TENSION_N / DENSITY_KG_M)
}
fn length_m() -> f64 {
    wave_speed() / (2.0 * F0_HZ)
}
fn resistance() -> f64 {
    R_OVER_Z0 * z0()
}
fn schelleng_fmax(beta: f64) -> f64 {
    2.0 * z0() * BOW_SPEED_M_S / (beta * (MU_S - MU_D))
}
fn bead_text_fmax(beta: f64) -> f64 {
    2.0 * z0() * z0() * BOW_SPEED_M_S / (resistance() * beta * (MU_S - MU_D))
}
fn schelleng_fmin(beta: f64) -> f64 {
    z0() * z0() * BOW_SPEED_M_S / (2.0 * resistance() * beta * beta * (MU_S - MU_D))
}

fn card(sample_rate_hz: u32, modes: usize) -> BowedStringCard {
    let (c, l) = (wave_speed(), length_m());
    let sigma = c / (2.0 * l) * det::ln((R_OVER_Z0 + 1.0) / (R_OVER_Z0 - 1.0));
    BowedStringCard {
        length_m: l,
        tension_n: TENSION_N,
        linear_density_kg_m: DENSITY_KG_M,
        bending_stiffness_n_m2: 0.0,
        viscous_bending_n_m2_s: 0.0,
        mode_count: modes,
        zetas: (1..=modes)
            .map(|k| {
                #[allow(clippy::cast_precision_loss)]
                let omega = k as f64 * core::f64::consts::PI * c / l;
                sigma / omega
            })
            .collect(),
        sample_rate_hz,
    }
}

/// Ideal Helmholtz motion at mid-stick, in mass-normalized modal coordinates.
fn helmholtz_start(modes: usize, beta: f64) -> Vec<ModalAcousticState> {
    let omega_1 = core::f64::consts::PI * wave_speed() / length_m();
    let amplitude = 2.0 * BOW_SPEED_M_S / (beta * core::f64::consts::PI * omega_1);
    let scale = det::sqrt(DENSITY_KG_M * length_m() / 2.0);
    (1..=modes)
        .map(|k| {
            #[allow(clippy::cast_precision_loss)]
            let kk = k as f64;
            let sign = if k % 2 == 1 { 1.0 } else { -1.0 };
            ModalAcousticState {
                displacement_m_sqrt_kg: 0.0,
                velocity_m_sqrt_kg_per_s: scale * amplitude * omega_1 / kk * sign,
            }
        })
        .collect()
}

fn voice(
    sample_rate_hz: u32,
    modes: usize,
    substeps: usize,
    beta: f64,
    force: f64,
) -> BowedStringState {
    let config = BowedRunConfig {
        card: card(sample_rate_hz, modes),
        island: FrictionIsland::Stribeck(
            StribeckFriction::try_new(MU_S, MU_D, DECAY_M_S, DECAY_M_S).unwrap(),
        ),
        gesture: BowGesture::admit(BOW_SPEED_M_S, force, beta).unwrap(),
        steps: 0,
        subsamples: substeps,
        termination: Termination::Rigid,
        listener_m: 1.0,
    };
    let mut state = BowedStringState::new(&config, 1).unwrap();
    state
        .restore_string_states(&helmholtz_start(modes, beta))
        .unwrap();
    state
}

/// Every substep receipt of `periods` nominal periods.
fn receipts(
    state: &mut BowedStringState,
    sample_rate_hz: u32,
    periods: f64,
) -> Vec<BowContactReceipt> {
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let samples = (f64::from(sample_rate_hz) * periods / F0_HZ).round() as usize;
    let mut out = Vec::new();
    for _ in 0..samples {
        state.step().unwrap();
        out.extend_from_slice(state.last_contact_receipts());
    }
    out
}

/// Indices of the first sliding step after a stuck one.
fn onsets(receipts: &[BowContactReceipt]) -> Vec<usize> {
    receipts
        .windows(2)
        .enumerate()
        .filter(|(_, w)| {
            w[0].branch == ContactBranch::Stuck
                && matches!(
                    w[1].branch,
                    ContactBranch::SlidingLow { .. } | ContactBranch::SlidingHigh { .. }
                )
        })
        .map(|(i, _)| i + 1)
        .collect()
}

/// Helmholtz motion in the classified tail: one slip onset per nominal
/// period, every interval within 5 % of the period. Returns the verdict and
/// the tail's sliding fraction.
fn classify(receipts: &[BowContactReceipt], step_s: f64) -> (bool, f64) {
    #[allow(clippy::cast_precision_loss)]
    let window = (TAIL_PERIODS as f64 / F0_HZ / step_s).round();
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let tail = &receipts[receipts.len() - window as usize..];
    let found = onsets(tail);
    let periodic = found.len() + 1 >= TAIL_PERIODS
        && found.len() <= TAIL_PERIODS
        && found.windows(2).all(|w| {
            #[allow(clippy::cast_precision_loss)]
            let interval = (w[1] - w[0]) as f64 * step_s;
            (interval * F0_HZ - 1.0).abs() <= 0.05
        });
    let sliding = tail
        .iter()
        .filter(|r| r.branch != ContactBranch::Stuck)
        .count();
    #[allow(clippy::cast_precision_loss)]
    let fraction = sliding as f64 / tail.len() as f64;
    (periodic, fraction)
}

fn helmholtz_at(beta: f64, force: f64, (rate, modes): (u32, usize)) -> (bool, f64) {
    let mut state = voice(rate, modes, 1, beta, force);
    let trace = receipts(&mut state, rate, RUN_PERIODS);
    classify(&trace, 1.0 / f64::from(rate))
}

/// Geometric bisection between a Helmholtz force and a non-Helmholtz one.
fn boundary(beta: f64, mut inside: f64, mut outside: f64, rung: (u32, usize)) -> f64 {
    for _ in 0..10 {
        let mid = det::sqrt(inside * outside);
        if helmholtz_at(beta, mid, rung).0 {
            inside = mid;
        } else {
            outside = mid;
        }
    }
    det::sqrt(inside * outside)
}

struct Bounds {
    beta: f64,
    fmin: f64,
    fmax: f64,
    center_fraction: f64,
    center_helmholtz: bool,
}

/// Measured bounds per beta on the finest rung, computed once per binary.
/// The bracket starts at 1.5 times Schelleng's geometric mean, which must be
/// Helmholtz, and spans [F_min/4, 4 F_max].
fn measured_bounds() -> &'static [Bounds] {
    static BOUNDS: OnceLock<Vec<Bounds>> = OnceLock::new();
    BOUNDS.get_or_init(|| {
        let rung = RUNGS[2];
        std::thread::scope(|scope| {
            let handles: Vec<_> = BETAS
                .iter()
                .map(|&beta| {
                    scope.spawn(move || {
                        let start = 1.5 * det::sqrt(schelleng_fmin(beta) * schelleng_fmax(beta));
                        assert!(
                            helmholtz_at(beta, start, rung).0,
                            "start force {start} is not Helmholtz at beta {beta}"
                        );
                        let fmax = boundary(beta, start, 4.0 * schelleng_fmax(beta), rung);
                        let fmin = boundary(beta, start, schelleng_fmin(beta) / 4.0, rung);
                        let (center_helmholtz, center_fraction) =
                            helmholtz_at(beta, det::sqrt(fmin * fmax), rung);
                        Bounds {
                            beta,
                            fmin,
                            fmax,
                            center_fraction,
                            center_helmholtz,
                        }
                    })
                })
                .collect();
            handles.into_iter().map(|h| h.join().unwrap()).collect()
        })
    })
}

#[test]
fn g1_helmholtz_slip_fraction_equals_beta() {
    // Finest rung, at the geometric centre of the measured Helmholtz band.
    let (rate, modes) = RUNGS[2];
    let per_period = f64::from(rate) / F0_HZ;
    let mut failures = Vec::new();
    for b in measured_bounds() {
        let tolerance = 0.05 * b.beta + 2.0 / per_period;
        let error = b.center_fraction - b.beta;
        println!(
            "{{\"test\":\"helmholtz_slip_fraction\",\"rung_hz\":{rate},\"modes\":{modes},\"beta\":{},\"force_n\":{},\"slip_fraction\":{},\"error\":{error},\"tolerance\":{tolerance},\"helmholtz\":{}}}",
            b.beta,
            det::sqrt(b.fmin * b.fmax),
            b.center_fraction,
            b.center_helmholtz
        );
        if !(b.center_helmholtz && error.abs() <= tolerance) {
            failures.push(b.beta);
        }
    }
    assert!(
        failures.is_empty(),
        "slip fraction or Helmholtz failed at beta {failures:?}"
    );
}

fn slope(xs: &[f64], ys: &[f64]) -> f64 {
    #[allow(clippy::cast_precision_loss)]
    let n = xs.len() as f64;
    let (mx, my) = (xs.iter().sum::<f64>() / n, ys.iter().sum::<f64>() / n);
    let sxy: f64 = xs.iter().zip(ys).map(|(x, y)| (x - mx) * (y - my)).sum();
    let sxx: f64 = xs.iter().map(|x| (x - mx) * (x - mx)).sum();
    sxy / sxx
}

/// F_max against eq. (4) and its exponent, gated as written in DONE-WHEN (3).
/// The exponent gate is an open owner question: the exact construction's own
/// exponent over these beta is printed beside it. F_min is printed only. This
/// stand-in's equal-decay damping is not a resistive termination, so the
/// DONE-WHEN's F_min is gated in `bowed_resistive_termination.rs` (verifier
/// ruling, mail 5685).
#[test]
fn g1_schelleng_fmax_exponent_and_bounds() {
    let (rate, modes) = RUNGS[2];
    let bounds = measured_bounds();
    let ln_beta: Vec<f64> = bounds.iter().map(|b| det::ln(b.beta)).collect();
    let fmax_exponent = slope(
        &ln_beta,
        &bounds.iter().map(|b| det::ln(b.fmax)).collect::<Vec<_>>(),
    );
    let fmin_exponent = slope(
        &ln_beta,
        &bounds.iter().map(|b| det::ln(b.fmin)).collect::<Vec<_>>(),
    );
    let exact_exponent = slope(
        &ln_beta,
        &bounds
            .iter()
            .map(|b| det::ln((1.0 - b.beta) / b.beta))
            .collect::<Vec<_>>(),
    );
    let mut out_of_band = Vec::new();
    for b in bounds {
        let (rmin, rmax) = (
            b.fmin / schelleng_fmin(b.beta),
            b.fmax / schelleng_fmax(b.beta),
        );
        println!(
            "{{\"test\":\"schelleng\",\"rung_hz\":{rate},\"modes\":{modes},\"beta\":{},\"z0\":{},\"r\":{},\"delta_mu\":{},\"fmin_n\":{},\"fmin_eq9_n\":{},\"fmin_ratio_not_gated\":{rmin},\"fmax_n\":{},\"fmax_eq4_n\":{},\"fmax_ratio\":{rmax},\"fmax_ratio_exact_slip\":{},\"fmax_bead_text_n\":{},\"fmax_bead_text_ratio\":{},\"band\":1.5}}",
            b.beta,
            z0(),
            resistance(),
            MU_S - MU_D,
            b.fmin,
            schelleng_fmin(b.beta),
            b.fmax,
            schelleng_fmax(b.beta),
            rmax / (1.0 - b.beta),
            bead_text_fmax(b.beta),
            b.fmax / bead_text_fmax(b.beta)
        );
        if !(1.0 / 1.5..=1.5).contains(&rmax) {
            out_of_band.push((b.beta, rmax));
        }
    }
    println!(
        "{{\"test\":\"schelleng\",\"fmax_exponent\":{fmax_exponent},\"fmax_band\":[-1.1,-0.9],\"exact_slip_speed_exponent\":{exact_exponent},\"fmin_exponent_not_gated\":{fmin_exponent}}}"
    );
    assert!(
        (-1.1..=-0.9).contains(&fmax_exponent),
        "F_max exponent {fmax_exponent}"
    );
    assert!(
        out_of_band.is_empty(),
        "F_max outside 1.5x of Schelleng eq. 4: {out_of_band:?}"
    );
}

/// Slip-onset times [s], root index taken at each onset, and the root
/// counts seen, from the Helmholtz start at 0.6 F_max (inside the band).
fn ladder_rung(
    sample_rate_hz: u32,
    modes: usize,
    substeps: usize,
    beta: f64,
) -> (Vec<f64>, Vec<u8>, Vec<u8>) {
    let mut state = voice(
        sample_rate_hz,
        modes,
        substeps,
        beta,
        0.6 * schelleng_fmax(beta),
    );
    let trace = receipts(&mut state, sample_rate_hz, 30.0);
    #[allow(clippy::cast_precision_loss)]
    let h = 1.0 / (f64::from(sample_rate_hz) * substeps as f64);
    let found = onsets(&trace);
    #[allow(clippy::cast_precision_loss)]
    let times = found.iter().map(|&i| i as f64 * h).collect();
    let taken = found.iter().map(|&i| trace[i].selected).collect();
    let mut counts: Vec<u8> = trace.iter().map(|r| r.roots).collect();
    counts.sort_unstable();
    counts.dedup();
    (times, taken, counts)
}

#[test]
fn g1_dt_ladder_root_sequence_converges() {
    // beta = 1/10 from the Helmholtz start; the first 20 slip onsets and the
    // root each slip takes are compared across rungs.
    let beta = 0.1;
    let compared = 20;
    let rungs: Vec<_> = RUNGS
        .iter()
        .map(|&(rate, modes)| ladder_rung(rate, modes, 1, beta))
        .collect();
    for (rung, (times, taken, counts)) in RUNGS.iter().zip(&rungs) {
        println!(
            "{{\"test\":\"dt_ladder\",\"rung_hz\":{},\"modes\":{},\"onsets\":{},\"first_onsets_s\":{:?},\"taken_roots\":{:?},\"root_counts_seen\":{counts:?}}}",
            rung.0,
            rung.1,
            times.len(),
            &times[..times.len().min(compared)],
            &taken[..taken.len().min(compared)]
        );
        assert!(
            times.len() >= compared,
            "{} onsets at {} Hz",
            times.len(),
            rung.0
        );
    }
    let rms = |a: &[f64], b: &[f64]| -> f64 {
        #[allow(clippy::cast_precision_loss)]
        let n = compared as f64;
        det::sqrt(
            a.iter()
                .zip(b)
                .take(compared)
                .map(|(x, y)| (x - y) * (x - y))
                .sum::<f64>()
                / n,
        )
    };
    let coarse = rms(&rungs[0].0, &rungs[1].0);
    let fine = rms(&rungs[1].0, &rungs[2].0);
    let order = det::ln(coarse / fine) / det::ln(2.0);
    // A fixed band refined only by substeps: the discrete mobility falls with h.
    let fixed: Vec<_> = [1_usize, 2, 4]
        .iter()
        .map(|&s| ladder_rung(48_000, 64, s, beta))
        .collect();
    println!(
        "{{\"test\":\"dt_ladder\",\"onset_rms_48_96_s\":{coarse:e},\"onset_rms_96_192_s\":{fine:e},\"observed_order\":{order},\"fixed_band_48k_x1_2_4_root_counts\":{:?},\"fixed_band_onset_rms\":[{:e},{:e}]}}",
        fixed.iter().map(|r| r.2.clone()).collect::<Vec<_>>(),
        rms(&fixed[0].0, &fixed[1].0),
        rms(&fixed[1].0, &fixed[2].0)
    );
    assert!(
        fine < coarse,
        "onset times converge: {coarse:e} then {fine:e}"
    );
    assert_eq!(rungs[0].2, rungs[1].2, "root multiplicity 48 vs 96 kHz");
    assert_eq!(rungs[1].2, rungs[2].2, "root multiplicity 96 vs 192 kHz");
    assert_eq!(
        rungs[1].1[..compared],
        rungs[2].1[..compared],
        "selected-root sequence 96 vs 192 kHz"
    );
}
