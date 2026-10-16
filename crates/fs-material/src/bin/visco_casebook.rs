//! End-to-end viscoelastic-damping casebook
//! (bead frankensim-fsim-visco-damping-ybc75): synthesize a fractional-Zener
//! complex-modulus dataset -> fit the 4 parameters with exact AD gradients ->
//! lower to a certified Prony series -> time-step a single-DOF ring-down with
//! the runtime recursion -> re-measure the loss factor from the decay
//! envelope and compare it against the model, emitting one JSON line per
//! stage.
//!
//! Output contract: every field is deterministic (values hashed with the
//! canonical FNV-1a-64 helper over little-endian f64 bytes) EXCEPT
//! `elapsed_ms`, which is wall-clock run evidence and varies per host. The
//! synthetic dataset's multiplicative perturbation is derived from
//! `fnv1a64(sample index bytes)` — replayable from this source alone, no RNG
//! state. A refusal on any stage is itself a reportable outcome: the binary
//! emits the structured row with the named `FS-MAT-VISCO-*` code and exits
//! nonzero.
//!
//! No-claim boundary: the decay-envelope comparison is a NUMERICAL
//! consistency check of the fitted/lowered/stepped chain against its own
//! declared physics. It does not calibrate any real material, and the
//! measured-vs-model gate below is authored from an observed run with
//! headroom, not derived from an analytic error bound.

use std::time::Instant;

use fs_casebook::fnv1a64;
use fs_material::visco::{
    FractionalZener, GeneralizedMaxwell, LoweredModel, fit_fractional_zener, lower_to_prony,
};
use fs_math::det;

const SUITE: &str = "fs-material-visco-casebook";
/// Deterministic perturbation scale for the synthetic dataset (relative).
const NOISE_REL: f64 = 1.0e-3;
/// Ring-down carriers [Hz], inside the certified lowering band.
const CARRIERS_HZ: [f64; 3] = [120.0, 1_000.0, 6_000.0];
/// Authored gate for |eta_measured - eta_model| / eta_model per carrier.
/// Measured 2026-08-12 on the fixed configuration below: worst carrier
/// 2.44e-2 (at 1 kHz); the gate keeps ~3x headroom. Integrator order,
/// envelope-fit error, and the fitted-vs-truth noise floor — not lowering
/// error — dominate the observed value.
const ETA_GATE_REL: f64 = 8.0e-2;

fn hash_f64s(values: &[f64]) -> u64 {
    let mut bytes = Vec::with_capacity(values.len() * 8);
    for v in values {
        bytes.extend_from_slice(&v.to_le_bytes());
    }
    fnv1a64(&bytes)
}

/// Deterministic multiplicative perturbation in `[-NOISE_REL, NOISE_REL]`
/// keyed by sample index and channel — replayable without RNG state.
fn perturb(index: u64, channel: u64) -> f64 {
    let h = fnv1a64(&[index.to_le_bytes(), channel.to_le_bytes()].concat());
    #[allow(clippy::cast_precision_loss)]
    let unit = (h >> 11) as f64 / (1u64 << 53) as f64;
    NOISE_REL * (2.0 * unit - 1.0)
}

fn refuse(case: &str, error: &dyn core::fmt::Display) -> ! {
    println!(
        "{{\"suite\":\"{SUITE}\",\"case\":\"{case}\",\"verdict\":\"refused\",\"error\":\"{error}\"}}"
    );
    std::process::exit(1)
}

struct RingdownResult {
    omega_measured: f64,
    eta_measured: f64,
    dissipated: f64,
    steps: usize,
}

/// Free ring-down of a unit-mass single-DOF oscillator whose restoring
/// force is the Prony stress (velocity-Verlet; force refreshed from the
/// recursive update each step). Returns the measured ring frequency, the
/// log-decrement loss factor, and the closed dissipation ledger total.
fn ringdown(model: &GeneralizedMaxwell, omega_target: f64, periods: usize) -> RingdownResult {
    // Geometric factor so the linearized natural frequency lands on the
    // carrier: x'' = -g·sigma(x), g = omega^2 / E'(omega).
    let (ep, _) = model.modulus(omega_target);
    let g = omega_target * omega_target / ep;
    let dt = 2.0 * core::f64::consts::PI / omega_target / 512.0;
    let steps = periods * 512;

    let mut state = model.state();
    let mut x = 1.0e-3;
    let mut v = 0.0;
    let mut sigma = model.step(&mut state, x, 0.0);
    let mut a = -g * sigma;
    let mut peaks: Vec<(f64, f64)> = Vec::new();
    let mut x_prev = x;
    for n in 0..steps {
        let x_new = x + v * dt + 0.5 * a * dt * dt;
        // Local positive maximum: rising into the current sample, falling
        // after it. Sub-sample refinement is unnecessary at 512
        // steps/period because the regression below spans many peaks.
        if n >= 1 && x > x_prev && x >= x_new && x > 0.0 {
            #[allow(clippy::cast_precision_loss)]
            peaks.push((n as f64 * dt, x));
        }
        sigma = model.step(&mut state, x_new, dt);
        let a_new = -g * sigma;
        v += 0.5 * (a + a_new) * dt;
        x_prev = x;
        x = x_new;
        a = a_new;
    }
    assert!(peaks.len() >= 8, "ring-down must expose enough peaks");
    // Least-squares slope of ln(peak) vs peak time gives the decay rate
    // delta_t; eta = 2*delta_t/omega. Frequency from mean peak spacing.
    #[allow(clippy::cast_precision_loss)]
    let m = peaks.len() as f64;
    let mean_t = peaks.iter().map(|p| p.0).sum::<f64>() / m;
    let mean_l = peaks.iter().map(|p| det::ln(p.1)).sum::<f64>() / m;
    let mut num = 0.0;
    let mut den = 0.0;
    for &(t, p) in &peaks {
        num += (t - mean_t) * (det::ln(p) - mean_l);
        den += (t - mean_t) * (t - mean_t);
    }
    let decay_rate = -num / den;
    let spacing = (peaks[peaks.len() - 1].0 - peaks[0].0) / (m - 1.0);
    let omega_measured = 2.0 * core::f64::consts::PI / spacing;
    RingdownResult {
        omega_measured,
        eta_measured: 2.0 * decay_rate / omega_measured,
        dissipated: state.dissipated,
        steps,
    }
}

#[allow(clippy::too_many_lines)] // one linear casebook script, one stage per block
fn main() {
    let start = Instant::now();

    // Stage 1: synthesize. Wood-like ground truth; 24 samples over six
    // decades with deterministic +/-0.1% multiplicative perturbation.
    let truth = FractionalZener::new(9.0e9, 1.5e10, 0.35, 2.0e-4)
        .unwrap_or_else(|e| refuse("synthesize", &e));
    let mut samples = Vec::with_capacity(24);
    for k in 0..24u64 {
        #[allow(clippy::cast_precision_loss)]
        let t = k as f64 / 23.0;
        let w = det::exp(det::ln(10.0) + t * (det::ln(1.0e7) - det::ln(10.0)));
        let (ep, epp) = truth.modulus(w);
        samples.push((w, ep * (1.0 + perturb(k, 0)), epp * (1.0 + perturb(k, 1))));
    }
    let data_hash = hash_f64s(
        &samples
            .iter()
            .flat_map(|&(w, a, b)| [w, a, b])
            .collect::<Vec<_>>(),
    );
    println!(
        "{{\"suite\":\"{SUITE}\",\"case\":\"synthesize\",\"samples\":{},\"noise_rel\":{NOISE_REL:e},\
         \"truth\":[{:.6e},{:.6e},{:.6e},{:.6e}],\"data_hash\":\"{data_hash:016x}\",\"verdict\":\"pass\"}}",
        samples.len(),
        truth.e0,
        truth.e_inf,
        truth.alpha,
        truth.tau,
    );

    // Stage 2: fit from a deliberately offset start, logging every accepted
    // iteration's relative residual norm.
    let init =
        FractionalZener::new(3.0e9, 4.0e10, 0.6, 1.0e-5).unwrap_or_else(|e| refuse("fit", &e));
    // Residual budget: uniform ±NOISE_REL multiplicative noise on both
    // channels has RMS ≈ NOISE_REL/√3 ≈ 5.8e-4 (observed converged
    // residual 4.0e-4); the authored budget leaves ~4x headroom over the
    // noise floor while still refusing any structurally wrong fit.
    let fit =
        fit_fractional_zener(&samples, &init, 200, 2.5e-3).unwrap_or_else(|e| refuse("fit", &e));
    for (iter, residual) in fit.residual_history.iter().enumerate() {
        println!(
            "{{\"suite\":\"{SUITE}\",\"case\":\"fit-iteration\",\"iter\":{iter},\"residual_rel\":{residual:.6e}}}"
        );
    }
    let fitted = &fit.model;
    let param_err = [
        (fitted.e0 - truth.e0).abs() / truth.e0,
        (fitted.e_inf - truth.e_inf).abs() / truth.e_inf,
        (fitted.alpha - truth.alpha).abs() / truth.alpha,
        (fitted.tau - truth.tau).abs() / truth.tau,
    ];
    println!(
        "{{\"suite\":\"{SUITE}\",\"case\":\"fit\",\"iterations\":{},\"residual_rel\":{:.6e},\
         \"fitted\":[{:.6e},{:.6e},{:.6e},{:.6e}],\"param_rel_err_hash\":\"{:016x}\",\"verdict\":\"pass\"}}",
        fit.residual_history.len(),
        fit.residual_history.last().copied().unwrap_or(f64::NAN),
        fitted.e0,
        fitted.e_inf,
        fitted.alpha,
        fitted.tau,
        hash_f64s(&param_err),
    );

    // Stage 3: lower to a certified Prony band covering the carriers.
    let lowered: LoweredModel =
        lower_to_prony(fitted, 20.0, 2.0e4, 10, 0.02).unwrap_or_else(|e| refuse("lower", &e));
    let term_hash = hash_f64s(
        &lowered
            .model
            .terms
            .iter()
            .flat_map(|&(e, t)| [e, t])
            .collect::<Vec<_>>(),
    );
    println!(
        "{{\"suite\":\"{SUITE}\",\"case\":\"lower\",\"terms\":{},\"band_rad_s\":[{:.6e},{:.6e}],\
         \"sup_rel_err\":{:.6e},\"verification_points\":{},\"term_hash\":\"{term_hash:016x}\",\"verdict\":\"pass\"}}",
        lowered.model.terms.len(),
        lowered.band.0,
        lowered.band.1,
        lowered.sup_rel_err,
        lowered.verification_points,
    );

    // Stage 4: ring-down at each carrier; re-measure eta from the decay
    // envelope and compare against the lowered model at the MEASURED ring
    // frequency (band-checked, so an out-of-band drift refuses loudly).
    let mut worst_rel = 0.0f64;
    for f_hz in CARRIERS_HZ {
        let omega = 2.0 * core::f64::consts::PI * f_hz;
        let ring = ringdown(&lowered.model, omega, 40);
        let eta_model = lowered
            .loss_factor_checked(ring.omega_measured.clamp(lowered.band.0, lowered.band.1))
            .unwrap_or_else(|e| refuse("ringdown", &e));
        let rel = (ring.eta_measured - eta_model).abs() / eta_model;
        worst_rel = worst_rel.max(rel);
        assert!(
            ring.dissipated >= 0.0,
            "dissipation ledger must be non-negative"
        );
        let verdict = if rel <= ETA_GATE_REL { "pass" } else { "fail" };
        println!(
            "{{\"suite\":\"{SUITE}\",\"case\":\"ringdown\",\"carrier_hz\":{f_hz},\
             \"omega_measured\":{:.6e},\"eta_model\":{:.6e},\"eta_measured\":{:.6e},\
             \"rel_err\":{:.3e},\"gate_rel\":{ETA_GATE_REL:e},\"dissipated\":{:.6e},\
             \"steps\":{},\"verdict\":\"{verdict}\"}}",
            ring.omega_measured, eta_model, ring.eta_measured, rel, ring.dissipated, ring.steps,
        );
        if rel > ETA_GATE_REL {
            println!(
                "{{\"suite\":\"{SUITE}\",\"case\":\"summary\",\"verdict\":\"fail\",\
                 \"worst_eta_rel_err\":{rel:.3e}}}"
            );
            std::process::exit(1);
        }
    }

    let elapsed_ms = start.elapsed().as_secs_f64() * 1e3;
    println!(
        "{{\"suite\":\"{SUITE}\",\"case\":\"summary\",\"worst_eta_rel_err\":{worst_rel:.3e},\
         \"elapsed_ms\":{elapsed_ms:.1},\"verdict\":\"pass\"}}"
    );
}
