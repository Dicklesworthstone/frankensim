//! Compose a [`BernoulliAperture`] onto a characteristic line.
//!
//! The time port is the TMM reflectance FIR
//! ([`crate::driving_point`]). Reed lay is an [`fs_dcontact`] obstacle
//! whose Hunt–Crossley `χ` is matched to the existing viscous damper
//! (not a private impact law). Observer absorption is applied once
//! by realize (ISO 9613-1). Vocal folds, lip reeds, and relief
//! valves are the same objects.

use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::BernoulliAperture;
use crate::thin_plate::PlateBank;
use crate::unilateral_contact::slit_contact_force;
use fs_dcontact::Obstacle;
use fs_duct::{Duct, Termination};
use fs_material::gas::GasState;
use fs_math::det;
use fs_scenario::BeatingReed;

/// Realize mouthpiece pressure from a beating reed on a TMM bore.
///
/// Since the block render API landed (music bead 3ez8g.2.1) this is a
/// thin wrapper over [`crate::render::ReedBoreVoice`]: the voice owns the
/// per-sample loop; this one-shot path builds it, renders `n` samples in
/// one block, and swaps the caller's plate bank in and out so refusal
/// semantics (bank untouched on error) and results stay byte-identical to
/// the pre-API code. Two paths, one loop body — they cannot drift.
///
/// # Errors
/// Domain, TMM, or solver refusals.
#[allow(clippy::too_many_arguments)] // one coherent realization record
pub fn realize_reed_bore(
    physics: &Duct,
    gas: &GasState,
    reed: BeatingReed,
    termination: Termination,
    plates: &mut PlateBank,
    listener_m: f64,
    sample_rate_hz: u32,
    n: usize,
    wall: Option<&fs_phs::WallPin>,
) -> Result<Vec<f64>, AcousticRealizeError> {
    let mut voice = crate::render::ReedBoreVoice::new(
        physics,
        gas,
        reed,
        termination,
        PlateBank::default(),
        listener_m,
        sample_rate_hz,
        n,
        wall,
    )?;
    // Move the caller's bank in only after every fallible admission step
    // has succeeded, and back out unconditionally, so a refusal leaves it
    // untouched and success returns the post-render plate state exactly
    // as the pre-API loop did.
    core::mem::swap(voice.plate_bank_mut(), plates);
    let mut p_bore_hist = vec![0.0; n];
    let result = voice.step_block(&mut p_bore_hist);
    core::mem::swap(voice.plate_bank_mut(), plates);
    result?;
    Ok(p_bore_hist)
}

pub(crate) fn blowing_envelope(reed: BeatingReed, t: f64) -> f64 {
    if reed.attack_s <= 0.0 {
        return reed.blowing_pressure_pa;
    }
    if t >= reed.attack_s {
        return reed.blowing_pressure_pa;
    }
    let x = t / reed.attack_s;
    reed.blowing_pressure_pa * 0.5 * (1.0 - det::cos(core::f64::consts::PI * x))
}

pub(crate) fn reed_structural(reed: BeatingReed) -> (f64, f64) {
    let face = reed.width_m * 0.025;
    let k = if reed.stiffness_n_m > 0.0 {
        reed.stiffness_n_m
    } else {
        reed.closing_pressure_pa * face / reed.rest_opening_m
    };
    let r_damp = 2.0 * 0.35 * det::sqrt((k * reed.mass_kg).max(0.0));
    (k, r_damp)
}

pub(crate) fn aperture_of(reed: BeatingReed) -> BernoulliAperture {
    BernoulliAperture {
        rest_opening_m: reed.rest_opening_m,
        width_m: reed.width_m,
        closing_pressure_pa: reed.closing_pressure_pa,
    }
}

/// Aperture-junction solver mode (bead frankensim-2s4i5).
///
/// `Strict` is the certification default: the retained deterministic
/// bisection with its 21-point grid-argmin fallback. `FastNewton` is a
/// DECLARED FAST MODE (the fs-rand ziggurat precedent): island Newton
/// on the analytic Jacobian with a monotone-descent guard that hands
/// every cornered sample to the strict path. It is not bitwise-equal
/// to the strict path and must never become the default until it
/// earns the same proof.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ReedSolverMode {
    /// Deterministic bisection — the certification default.
    #[default]
    Strict,
    /// Declared fast mode: guarded analytic-Jacobian Newton.
    FastNewton,
}

/// Per-voice counters backing the fast-mode fallback-hit-rate receipt
/// (bead frankensim-2s4i5).
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct FastSolveStats {
    /// Samples resolved by Newton without fallback.
    pub newton_samples: u64,
    /// Samples handed to the strict path (cornered or non-monotone).
    pub fallback_samples: u64,
}

impl FastSolveStats {
    /// Fraction of solver-work samples that needed the strict path.
    /// The closed-branch early-out is shared structure, not solver
    /// work, and counts in neither numerator nor denominator.
    #[must_use]
    pub fn fallback_rate(&self) -> f64 {
        let total = self.newton_samples + self.fallback_samples;
        if total == 0 {
            0.0
        } else {
            self.fallback_samples as f64 / total as f64
        }
    }
}

/// Flow dead zone: `fs_phs::bernoulli_volume_flow` returns exactly 0
/// for `|dp| < 1e-12`; there the residual slope is the wave term alone.
const FLOW_DEAD_ZONE_PA: f64 = 1.0e-12;
/// Corner guard: below this |dp| the Bernoulli sqrt kink makes the
/// analytic slope unreliable, so the sample goes to the strict path.
/// Far below any physically meaningful reed pressure (uPa scale).
const NEWTON_CORNER_PA: f64 = 1.0e-6;
/// Newton cap before handing the sample to the strict path.
const NEWTON_MAX_ITERS: usize = 12;
/// Line-search halvings allowed before the sample is declared
/// cornered (4 halvings cover a 16x overshoot).
const NEWTON_BACKTRACK_MAX: u32 = 4;
/// Step-size convergence: stop when a full Newton step moves `p_plus`
/// by less than this relative amount. Converging on the shared
/// residual acceptance (`1e-8·(1+|p_m|)` in flow units) instead would
/// let the fast root sit `tol/|f'|` — order 100 Pa — from the strict
/// root, which is audible. Step-sized convergence bounds the deviation
/// at uPa scale for ~2 extra cheap iterations.
const NEWTON_STEP_TOL: f64 = 1.0e-9;

/// Analytic `d f / d p_plus` of [`reed_flow_mismatch`] at `p_plus`.
///
/// With `dp = p_m − ((1+r0)·p_plus + p_minus_hist)`:
/// `f = U(dp) + u_body − ((1−r0)·p_plus − p_minus_hist)/zc`, so
/// `f' = −(1+r0)·U'(dp) − (1−r0)/zc`. The Bernoulli law is smooth in
/// the aperture interior (`U' = w·(o'·g + o·g')`) but has an unbounded
/// slope as `|dp| → 0`; [`NEWTON_CORNER_PA`] refuses those samples.
#[allow(clippy::too_many_arguments)] // mirrors reed_flow_mismatch
fn reed_flow_jacobian(
    reed: BeatingReed,
    rho: f64,
    zc: f64,
    r0: f64,
    p_minus_hist: f64,
    p_m: f64,
    p_plus: f64,
) -> Option<f64> {
    let dp = p_m - ((1.0 + r0) * p_plus + p_minus_hist);
    if !dp.is_finite() || dp.abs() < NEWTON_CORNER_PA {
        return None;
    }
    let open_raw = 1.0 - dp / reed.closing_pressure_pa;
    if !open_raw.is_finite() {
        return None;
    }
    let opening = reed.rest_opening_m * open_raw.clamp(0.0, 1.0);
    // Interior derivative of the clamp; saturated ends are flat.
    let d_opening = if open_raw > 0.0 && open_raw < 1.0 {
        -reed.rest_opening_m / reed.closing_pressure_pa
    } else {
        0.0
    };
    let du_dp_wave_term = -(1.0 - r0) / zc;
    if dp.abs() < FLOW_DEAD_ZONE_PA {
        return Some(du_dp_wave_term);
    }
    let g = dp.signum() * fs_math::det::sqrt(2.0 * dp.abs() / rho);
    let g_prime = 1.0 / (rho * fs_math::det::sqrt(2.0 * dp.abs() / rho));
    let du_dp_flow = reed.width_m * (d_opening * g + opening * g_prime);
    let jacobian = -(1.0 + r0) * du_dp_flow + du_dp_wave_term;
    if jacobian.is_finite() && jacobian.abs() > f64::MIN_POSITIVE {
        Some(jacobian)
    } else {
        None
    }
}

/// Declared fast mode: guarded island Newton on the analytic Jacobian.
///
/// Convergence tolerance matches the strict path's acceptance test
/// exactly. Any cornered sample — sqrt kink, non-finite state,
/// vanishing slope, non-monotone residual, or iteration cap — is
/// handed untouched to the strict bisection and counted in
/// `stats.fallback_samples`. The strict path remains the
/// certification default; this mode exists so a measured budget row
/// can price the upgrade.
#[allow(clippy::too_many_arguments)] // one coherent junction record
pub(crate) fn solve_reed_wave_fast(
    reed: BeatingReed,
    rho: f64,
    zc: f64,
    r0: f64,
    p_minus_hist: f64,
    p_m: f64,
    guess: f64,
    u_body: f64,
    stats: &mut FastSolveStats,
) -> Result<f64, AcousticRealizeError> {
    // Shared closed-branch early-out: identical bytes in both modes
    // (no aperture flow at all), so it is not a solver-work sample.
    let denom = (1.0 - r0).clamp(-0.999, 0.999);
    let closed_plus = p_minus_hist / denom;
    let p_bore_closed = closed_plus + (p_minus_hist + r0 * closed_plus);
    if p_m - p_bore_closed >= reed.closing_pressure_pa {
        return Ok(closed_plus);
    }
    let mut p = guess;
    let mut f = reed_flow_mismatch(reed, rho, zc, r0, p_minus_hist, p_m, p, u_body);
    if !f.is_finite() {
        stats.fallback_samples += 1;
        return solve_reed_wave_strict(reed, rho, zc, r0, p_minus_hist, p_m, guess, u_body);
    }
    for _ in 0..NEWTON_MAX_ITERS {
        let Some(j) = reed_flow_jacobian(reed, rho, zc, r0, p_minus_hist, p_m, p) else {
            stats.fallback_samples += 1;
            return solve_reed_wave_strict(reed, rho, zc, r0, p_minus_hist, p_m, guess, u_body);
        };
        // Damped Newton: the raw step may overshoot the root from a
        // cold seed (measured: |f_new| can exceed |f| by <20% with the
        // true root lying BETWEEN the iterates). Halve the step until
        // the residual magnitude strictly decreases; a cornered sample
        // that never improves hands itself to the strict path.
        let full_step = -f / j;
        if !full_step.is_finite() {
            stats.fallback_samples += 1;
            return solve_reed_wave_strict(reed, rho, zc, r0, p_minus_hist, p_m, guess, u_body);
        }
        let mut applied = full_step;
        let mut p_new = p + applied;
        let mut f_new = reed_flow_mismatch(reed, rho, zc, r0, p_minus_hist, p_m, p_new, u_body);
        let mut halvings = 0_u32;
        while (!f_new.is_finite() || f_new.abs() >= f.abs()) && halvings < NEWTON_BACKTRACK_MAX {
            applied *= 0.5;
            p_new = p + applied;
            f_new = reed_flow_mismatch(reed, rho, zc, r0, p_minus_hist, p_m, p_new, u_body);
            halvings += 1;
        }
        if !f_new.is_finite() || f_new.abs() >= f.abs() {
            stats.fallback_samples += 1;
            return solve_reed_wave_strict(reed, rho, zc, r0, p_minus_hist, p_m, guess, u_body);
        }
        // Converged on an UNDOUBLED full step; a backtracked step means
        // "still descending", never "done".
        let converged = halvings == 0 && full_step.abs() <= NEWTON_STEP_TOL * (1.0 + p_new.abs());
        p = p_new;
        f = f_new;
        if converged {
            stats.newton_samples += 1;
            return Ok(p);
        }
    }
    stats.fallback_samples += 1;
    solve_reed_wave_strict(reed, rho, zc, r0, p_minus_hist, p_m, guess, u_body)
}

#[allow(clippy::too_many_arguments)] // one coherent reed-step record
pub(crate) fn step_massive_reed(
    reed: BeatingReed,
    rho: f64,
    zc: f64,
    p_minus: f64,
    p_m: f64,
    y: f64,
    v: f64,
    dt: f64,
    u_body: f64,
    lay: Option<&Obstacle>,
    mode: ReedSolverMode,
    stats: &mut FastSolveStats,
) -> Result<(f64, f64, f64), AcousticRealizeError> {
    let face = reed.width_m * 0.025;
    let (k, r_damp) = reed_structural(reed);
    let p_plus = match mode {
        ReedSolverMode::Strict => {
            solve_reed_wave(reed, rho, zc, 0.0, p_minus, p_m, 2.0 * p_minus, u_body)?
        }
        ReedSolverMode::FastNewton => solve_reed_wave_fast(
            reed,
            rho,
            zc,
            0.0,
            p_minus,
            p_m,
            2.0 * p_minus,
            u_body,
            stats,
        )?,
    };
    let p_bore = p_plus + p_minus;
    let dp = p_m - p_bore;
    let mut acc = (-k * (y - reed.rest_opening_m) - r_damp * v - face * dp) / reed.mass_kg;
    let mut y1 = y + dt * v;
    let mut v1 = v + dt * acc;
    if let Some(obstacle) = lay {
        let contact = slit_contact_force(obstacle, y1)
            .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
        let hc = obstacle.dissipative_modal_forces(1, &[y1, v], &[v])[0];
        acc += (contact + hc) / reed.mass_kg;
        v1 = v + dt * acc;
        y1 = y + dt * v1;
    }
    if !y1.is_finite() || !v1.is_finite() {
        return Err(AcousticRealizeError::Reed {
            what: "massive reed left the finite set",
        });
    }
    Ok((p_plus, y1, v1))
}

#[allow(clippy::too_many_arguments)] // one coherent junction record
#[allow(clippy::unnecessary_wraps)] // uniform Result surface across the reed solvers
pub(crate) fn solve_reed_wave_strict(
    reed: BeatingReed,
    rho: f64,
    zc: f64,
    r0: f64,
    p_minus_hist: f64,
    p_m: f64,
    guess: f64,
    u_body: f64,
) -> Result<f64, AcousticRealizeError> {
    let denom = (1.0 - r0).clamp(-0.999, 0.999);
    let closed_plus = p_minus_hist / denom;
    let p_bore_closed = closed_plus + (p_minus_hist + r0 * closed_plus);
    if p_m - p_bore_closed >= reed.closing_pressure_pa {
        return Ok(closed_plus);
    }
    let span = (2.0 * reed.closing_pressure_pa)
        .max(2.0 * p_m.abs())
        .max(1.0);
    let mut lo = guess - span;
    let mut hi = guess + span;
    let mut f_lo = reed_flow_mismatch(reed, rho, zc, r0, p_minus_hist, p_m, lo, u_body);
    let mut f_hi = reed_flow_mismatch(reed, rho, zc, r0, p_minus_hist, p_m, hi, u_body);
    if f_lo * f_hi > 0.0 {
        lo = -span;
        hi = span;
        f_lo = reed_flow_mismatch(reed, rho, zc, r0, p_minus_hist, p_m, lo, u_body);
        f_hi = reed_flow_mismatch(reed, rho, zc, r0, p_minus_hist, p_m, hi, u_body);
    }
    if f_lo * f_hi > 0.0 {
        let mut best = guess;
        let mut best_a =
            reed_flow_mismatch(reed, rho, zc, r0, p_minus_hist, p_m, guess, u_body).abs();
        for k in 0..21 {
            let x = -span + (2.0 * span) * f64::from(k) / 20.0;
            let a = reed_flow_mismatch(reed, rho, zc, r0, p_minus_hist, p_m, x, u_body).abs();
            if a < best_a {
                best_a = a;
                best = x;
            }
        }
        return Ok(best);
    }
    let mut mid = f64::midpoint(lo, hi);
    for _ in 0..48 {
        mid = f64::midpoint(lo, hi);
        let f_mid = reed_flow_mismatch(reed, rho, zc, r0, p_minus_hist, p_m, mid, u_body);
        if f_mid.abs() < 1.0e-8 * (1.0 + p_m.abs()) {
            return Ok(mid);
        }
        if f_lo * f_mid <= 0.0 {
            hi = mid;
            f_hi = f_mid;
        } else {
            lo = mid;
            f_lo = f_mid;
        }
    }
    let _ = f_hi;
    Ok(mid)
}

/// Certification default: the retained deterministic bisection path,
/// unchanged bit-for-bit from the pre-fast-mode code. All callers that
/// do not explicitly opt into [`ReedSolverMode::FastNewton`] land here.
#[allow(clippy::too_many_arguments)] // one coherent junction record
#[allow(clippy::unnecessary_wraps)] // uniform Result surface across the reed solvers
pub(crate) fn solve_reed_wave(
    reed: BeatingReed,
    rho: f64,
    zc: f64,
    r0: f64,
    p_minus_hist: f64,
    p_m: f64,
    guess: f64,
    u_body: f64,
) -> Result<f64, AcousticRealizeError> {
    solve_reed_wave_strict(reed, rho, zc, r0, p_minus_hist, p_m, guess, u_body)
}

#[allow(clippy::too_many_arguments)] // one coherent junction record
fn reed_flow_mismatch(
    reed: BeatingReed,
    rho: f64,
    zc: f64,
    r0: f64,
    p_minus_hist: f64,
    p_m: f64,
    p_plus: f64,
    u_body: f64,
) -> f64 {
    let p_minus = p_minus_hist + r0 * p_plus;
    let p_bore = p_plus + p_minus;
    let dp = p_m - p_bore;
    let flow = aperture_of(reed).volume_flow(dp, rho);
    let u_wave = (p_plus - p_minus) / zc;
    flow + u_body - u_wave
}

#[cfg(test)]
mod fast_mode_tests {
    use super::*;

    use fs_material::gas::{GasSpec, GasState};
    use fs_scenario::BeatingReed;

    const BORE_RADIUS_M: f64 = 0.0022;

    fn air20_at(pressure_pa: f64) -> GasState {
        GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, pressure_pa).expect("air")
    }

    fn air20() -> GasState {
        air20_at(101_325.0)
    }

    fn reed() -> BeatingReed {
        BeatingReed {
            rest_opening_m: 4.0e-4,
            width_m: 0.013,
            closing_pressure_pa: 6_000.0,
            blowing_pressure_pa: 2_800.0,
            attack_s: 0.008,
            mass_kg: 0.0,
            stiffness_n_m: 0.0,
        }
    }

    fn fixture_impedance(gas: &GasState) -> f64 {
        let area = core::f64::consts::PI * BORE_RADIUS_M * BORE_RADIUS_M;
        gas.characteristic_impedance / area
    }

    #[test]
    fn fixture_impedance_is_derived_from_declared_gas_state() {
        let sea_level = air20();
        let reduced_pressure = air20_at(90_000.0);
        let area = core::f64::consts::PI * BORE_RADIUS_M * BORE_RADIUS_M;
        for gas in [sea_level, reduced_pressure] {
            assert_eq!(
                fixture_impedance(&gas).to_bits(),
                (gas.characteristic_impedance / area).to_bits()
            );
        }
        assert_ne!(
            fixture_impedance(&sea_level).to_bits(),
            fixture_impedance(&reduced_pressure).to_bits(),
            "fixture impedance must retain the declared pressure dependence"
        );
    }

    #[test]
    fn jacobian_matches_central_finite_difference_in_smooth_regions() {
        let reed = reed();
        let gas = air20();
        let rho = gas.density;
        let zc = fixture_impedance(&gas);
        // Interior points: away from the dp ≈ 0 kink, both clamp ends,
        // and the flow dead zone.
        let cases = [
            (0.0_f64, 500.0_f64, 1_000.0_f64),
            (0.0, -300.0, 2_000.0),
            (-0.4, 150.0, 3_000.0),
            (0.9, -80.0, -500.0),
        ];
        for &(r0, h, pm) in &cases {
            for &p in &[-400.0, -50.0, 0.0, 250.0, 900.0] {
                let j = reed_flow_jacobian(reed, rho, zc, r0, h, pm, p)
                    .unwrap_or_else(|| panic!("corner refused at r0={r0} h={h} pm={pm} p={p}"));
                let eps = 1.0e-3 * (1.0 + p.abs());
                let f_hi = reed_flow_mismatch(reed, rho, zc, r0, h, pm, p + eps, 0.0);
                let f_lo = reed_flow_mismatch(reed, rho, zc, r0, h, pm, p - eps, 0.0);
                let fd = (f_hi - f_lo) / (2.0 * eps);
                let rel = ((j - fd) / fd.abs().max(1.0e-30)).abs();
                assert!(
                    rel < 1.0e-5,
                    "J mismatch at r0={r0} h={h} pm={pm} p={p}: analytic {j} vs FD {fd}"
                );
            }
        }
    }

    #[test]
    fn corner_regions_are_refused_not_guessed() {
        let reed = reed();
        let gas = air20();
        let rho = gas.density;
        let zc = fixture_impedance(&gas);
        // Drive dp → 0: with r0 = 0, dp = pm − p − h, so the
        // equal-pressure point is p = pm − h; the sqrt-kink guard must
        // refuse its neighborhood.
        let p_equal = reed.blowing_pressure_pa - 100.0; // h := 100
        for p in [p_equal - 1.0e-9, p_equal, p_equal + 1.0e-9] {
            assert!(
                reed_flow_jacobian(reed, rho, zc, 0.0, 100.0, reed.blowing_pressure_pa, p)
                    .is_none(),
                "kink neighborhood must refuse at p={p}"
            );
        }
    }

    #[test]
    fn fast_newton_meets_the_strict_residual_contract_within_conditioning_band() {
        const FLOW_ACCEPT_TOL: f64 = 1.0e-8;
        const POSITION_HEADROOM: f64 = 4.0;
        let reed = reed();
        let gas = air20();
        let rho = gas.density;
        let zc = fixture_impedance(&gas);
        let mut stats = FastSolveStats::default();
        // A sweep across open/interior/closing phases with varied
        // history and reflection. TWO gates, both derived from the
        // strict path's own acceptance test:
        // 1. RESIDUAL CONTRACT: every Newton-resolved sample satisfies
        //    |f| < 1e-8*(1+p_m) — byte-for-byte the same acceptance
        //    the strict bisection uses.
        // 2. CONDITIONING-SCALED POSITION BAND: the residual is nearly
        //    flat here (|J| ~ 1e-7 because zc ~ 2.7e7 makes the wave
        //    term tiny), so ANY solver's root position carries an
        //    uncertainty of flow_tol*(1+p_m)/|J| — tens of Pa at some
        //    operating points. The positional assertion uses four
        //    times that per-sample bound; a fixed microbar band would
        //    encode false precision about an ill-conditioned root.
        let mut max_dev = 0.0_f64;
        let mut max_allowed = 0.0_f64;
        for k in 0..64 {
            let pm = 200.0 + 120.0 * f64::from(k);
            let h = 40.0 * f64::from(k % 7) - 120.0;
            let guess = 0.25 * pm + 30.0 * f64::from(k % 5) + 7.0;
            let strict = solve_reed_wave_strict(reed, rho, zc, 0.0, h, pm, guess, 0.0)
                .expect("strict solves");
            let fast = solve_reed_wave_fast(reed, rho, zc, 0.0, h, pm, guess, 0.0, &mut stats)
                .expect("fast solves");
            let dev = (fast - strict).abs();
            let f_fast = reed_flow_mismatch(reed, rho, zc, 0.0, h, pm, fast, 0.0);
            assert!(
                f_fast.abs() < FLOW_ACCEPT_TOL * (1.0 + pm.abs()),
                "Newton root violates the strict residual contract at \
                 k={k}: |f|={:e} vs tol {:e}",
                f_fast.abs(),
                FLOW_ACCEPT_TOL * (1.0 + pm.abs())
            );
            let j = reed_flow_jacobian(reed, rho, zc, 0.0, h, pm, fast)
                .map_or(0.0, f64::abs)
                .max(reed_flow_jacobian(reed, rho, zc, 0.0, h, pm, strict).map_or(0.0, f64::abs));
            let allowed = POSITION_HEADROOM * (FLOW_ACCEPT_TOL * (1.0 + pm.abs())) / j.max(1.0e-12);
            max_dev = max_dev.max(dev);
            max_allowed = max_allowed.max(allowed);
            assert!(
                dev <= allowed,
                "root deviation {dev:.3} Pa exceeds the conditioning band \
                 {allowed:.3} Pa at k={k} (pm={pm}, h={h}, |J|={j:e})"
            );
        }
        println!(
            "receipt: aperture-newton battery max|p_fast-p_strict| = {max_dev:e} Pa \
(max conditioning-scaled allowance {max_allowed:e} Pa)"
        );
        // MEASURED (this battery, damped-Newton solver): ~43% of
        // cold-seeded samples defer to bisection - a flat residual
        // (|J|~1e-7) makes far-off guesses genuinely hard, and the
        // deferred samples return byte-identical strict output, so
        // the cost is performance-only. The production-shaped
        // fallback receipt lives in tests/reed_newton_fusion.rs
        // (warm-seeded render loop). Bound set from that measurement,
        // not from optimism.
        assert!(
            stats.fallback_rate() <= 0.5,
            "unexpected fallback rate on smooth battery: {stats:?}"
        );
    }

    #[test]
    fn cornered_guess_hands_to_strict_and_still_matches() {
        let reed = reed();
        let gas = air20();
        let rho = gas.density;
        let zc = fixture_impedance(&gas);
        let mut stats = FastSolveStats::default();
        // A guess landing EXACTLY on the kink (dp = 0) must hand the
        // sample to the strict path and return the strict root.
        let pm = 2_800.0;
        let h = 100.0;
        let corner_guess = pm - h; // dp := pm − 2·guess − 0·… = 0 here
        let strict = solve_reed_wave_strict(reed, rho, zc, 0.0, h, pm, corner_guess, 0.0)
            .expect("strict solves");
        let fast = solve_reed_wave_fast(reed, rho, zc, 0.0, h, pm, corner_guess, 0.0, &mut stats)
            .expect("cornered fast defers to strict");
        assert_eq!(stats.fallback_samples, 1);
        assert_eq!(stats.newton_samples, 0);
        assert_eq!(fast.to_bits(), strict.to_bits()); // byte-identical handoff, never weakened
    }
}
