//! Bounded, measured refinement of retained harmonic order.

use super::{
    AftWorkspace, HbAnchor, HbBudget, HbOrbit, HbScales, OrbitError, OrbitProblem, Packing, TAU,
    solve_hb, solve_hb_scaled_seeded, solve_hb_seeded,
};
use fs_math::{c64::C64, det};

/// Physical accuracy and work limits for [solve_hb_adaptive].
#[derive(Debug, Clone, PartialEq)]
pub struct HbRefinementBudget {
    /// Largest solved order; must exceed the starting HbBudget harmonic count.
    pub max_harmonics: usize,
    /// Absolute RMS equation-defect tolerance per equation, in that row's units.
    /// No state norm or DC preload rescales these tolerances.
    pub residual_tolerances: Vec<f64>,
    /// Absolute RMS change between successive solved waveforms per state,
    /// in that coordinate's units. Includes newly added harmonics.
    pub state_tolerances: Vec<f64>,
    /// Relative change in the solved angular frequency.
    pub frequency_tolerance: f64,
}

/// Measured finite-band evidence at one solved harmonic order.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HbRefinementCheck {
    /// Retained order of the solved candidate.
    pub harmonics: usize,
    /// Largest assessed equation harmonic, always 2 * harmonics + 1.
    pub checked_harmonics: usize,
    /// Samples in the independent finer assessment quadrature.
    pub aft_samples: usize,
    /// Largest per-equation total RMS defect across the full checked band /
    /// caller tolerance. Retained and omitted bands share this one budget.
    pub total_residual_ratio: f64,
    /// Largest per-equation retained RMS defect / caller tolerance.
    pub retained_residual_ratio: f64,
    /// Largest per-equation omitted-band RMS defect / caller tolerance.
    pub omitted_residual_ratio: f64,
    /// Largest per-equation RMS change between assessment grids / tolerance.
    pub sampling_defect_ratio: f64,
    /// Largest per-state RMS waveform change / caller tolerance.
    /// Absent until two distinct retained orders have actually been solved.
    pub state_change_ratio: Option<f64>,
    /// Relative frequency change between the same two solved orders.
    pub relative_frequency_change: Option<f64>,
}

/// An orbit with measured order-refinement evidence, not a tail-error bound.
#[derive(Debug, Clone)]
pub struct AdaptiveHbOrbit {
    /// Last independently solved orbit, with its original HB diagnostics.
    pub orbit: HbOrbit,
    /// Every assessed solved order, including the initial candidate.
    pub refinement_trace: Vec<HbRefinementCheck>,
}

fn admit(
    d: usize,
    anchor: HbAnchor,
    hb: &HbBudget,
    budget: &HbRefinementBudget,
) -> Result<(), OrbitError> {
    let positive = |x: f64| x.is_finite() && x > 0.0;
    let omega = match anchor {
        HbAnchor::Forced { omega } => omega,
        HbAnchor::Autonomous { omega_guess } | HbAnchor::Backbone { omega_guess, .. } => {
            omega_guess
        }
    };
    if d == 0
        || hb.harmonics == 0
        || budget.max_harmonics <= hb.harmonics
        || hb.max_newton == 0
        || !positive(hb.tolerance)
        || !positive(omega)
        || !positive(budget.frequency_tolerance)
        || budget.residual_tolerances.len() != d
        || budget.state_tolerances.len() != d
        || budget
            .residual_tolerances
            .iter()
            .chain(&budget.state_tolerances)
            .any(|&x| !positive(x))
        || matches!(anchor, HbAnchor::Backbone { amplitude, .. } if !positive(amplitude))
    {
        return Err(OrbitError::BadParameter {
            what: "adaptive HB needs distinct positive orders, finite positive tolerances per row/state, and a valid anchor",
        });
    }
    let checked = budget
        .max_harmonics
        .checked_mul(2)
        .and_then(|n| n.checked_add(1));
    let extent = checked
        .and_then(|n| n.checked_mul(2))
        .and_then(|n| n.checked_add(1))
        .and_then(|n| n.checked_mul(d));
    if extent.and_then(|n| n.checked_mul(n)).is_none()
        || checked
            .and_then(|n| n.checked_mul(4))
            .and_then(|n| n.checked_add(4))
            .and_then(usize::checked_next_power_of_two)
            .is_none()
    {
        return Err(OrbitError::BadParameter {
            what: "adaptive HB harmonic extent overflows",
        });
    }
    Ok(())
}

/// Refine retained order, warm-starting each larger solve from the last orbit.
///
/// Orders grow as min(2N+1, max_harmonics). Each candidate's original equation
/// is independently assessed through 2N+1; no port is evaluated above N,
/// because the omitted state coefficients are exactly zero. Acceptance requires
/// resolved quadrature, a combined retained/omitted equation defect within each
/// row's single RMS budget, and waveform
/// and frequency agreement between two actual solves. Backbone shape comparison
/// aligns its free phase; forced and autonomous warm-start gauges are preserved.
///
/// Finite-band agreement is measured evidence only. It cannot discover arbitrary
/// forcing outside the checked band, establish a forward-error bound near an
/// ill-conditioned response, or certify an infinite Fourier tail. Owners of
/// explicit forcing must include its support in the starting order.
///
/// # Errors
/// Invalid budgets, a fixed-order solver refusal, unresolved assessment sampling,
/// or [OrbitError::HarmonicsUnresolved] with the last measured check. A retained
/// residual failure refuses immediately: increasing order is not a repair for a
/// Newton tolerance too loose for the supplied physical equation tolerances.
pub fn solve_hb_adaptive<P: OrbitProblem>(
    problem: &P,
    anchor: HbAnchor,
    guess_amplitude: f64,
    hb: &HbBudget,
    budget: &HbRefinementBudget,
) -> Result<AdaptiveHbOrbit, OrbitError> {
    admit(problem.dim(), anchor, hb, budget)?;
    if !guess_amplitude.is_finite() {
        return Err(OrbitError::BadParameter {
            what: "adaptive HB seed must be finite",
        });
    }
    refine(
        problem,
        anchor,
        solve_hb(problem, anchor, guess_amplitude, hb)?,
        hb,
        budget,
        None,
    )
}

/// [solve_hb_adaptive] with the fixed-order solver's full canonical packed seed.
///
/// # Errors
/// The same refusals as [solve_hb_adaptive], plus malformed/nonfinite seeds.
pub fn solve_hb_adaptive_seeded<P: OrbitProblem>(
    problem: &P,
    anchor: HbAnchor,
    seed: Vec<f64>,
    hb: &HbBudget,
    budget: &HbRefinementBudget,
) -> Result<AdaptiveHbOrbit, OrbitError> {
    admit(problem.dim(), anchor, hb, budget)?;
    let coefficients = problem.dim() * (2 * hb.harmonics + 1);
    let unknown_frequency = !matches!(anchor, HbAnchor::Forced { .. });
    if seed.len() != coefficients + usize::from(unknown_frequency)
        || seed.iter().any(|x| !x.is_finite())
        || (unknown_frequency && seed[coefficients] <= 0.0)
    {
        return Err(OrbitError::BadParameter {
            what: "adaptive HB needs a finite full seed and positive unknown frequency",
        });
    }
    refine(
        problem,
        anchor,
        solve_hb_seeded(problem, anchor, seed, hb)?,
        hb,
        budget,
        None,
    )
}

/// [solve_hb_adaptive_seeded] with explicit physical Newton conditioning at
/// every solved order. Seeds, output coefficients and refinement tolerances
/// retain their original units; [`HbScales`] controls only the fixed-order
/// solver's numerical chart and absolute normalized balance/sampling gate.
///
/// Warm starts extend the same physical waveform, with no change of scaling
/// when modes are added. The independent retained/omitted-band assessment and
/// actual inter-order waveform comparison remain mandatory. A loose Newton
/// gate still cannot bypass the caller's physical refinement tolerances.
///
/// # Errors
/// The refusals of [solve_hb_adaptive_seeded] and [solve_hb_scaled_seeded].
pub fn solve_hb_adaptive_scaled_seeded<P: OrbitProblem>(
    problem: &P,
    anchor: HbAnchor,
    seed: Vec<f64>,
    hb: &HbBudget,
    budget: &HbRefinementBudget,
    scales: &HbScales,
) -> Result<AdaptiveHbOrbit, OrbitError> {
    admit(problem.dim(), anchor, hb, budget)?;
    let initial = solve_hb_scaled_seeded(problem, anchor, seed, hb, scales)?.orbit;
    refine(problem, anchor, initial, hb, budget, Some(scales))
}

// Parseval RMS for a real Fourier series, normalized BEFORE squaring so tiny
// physical units do not disappear by underflow. DC has weight one; AC weight two.
fn row_ratio(values: &[f64], pack: &Packing, first: usize, tolerances: &[f64]) -> f64 {
    let mut worst = 0.0_f64;
    for (k, &tol) in tolerances.iter().enumerate() {
        let mut sum = if first == 0 {
            (values[pack.re_index(0, k)] / tol).powi(2)
        } else {
            0.0
        };
        for n in first.max(1)..=pack.n_harm {
            sum += 2.0
                * ((values[pack.re_index(n, k)] / tol).powi(2)
                    + (values[pack.im_index(n, k)] / tol).powi(2));
        }
        if !sum.is_finite() {
            return f64::INFINITY;
        }
        worst = worst.max(det::sqrt(sum));
    }
    worst
}

// Assess the same original island. A missing state mode contributes -G_n;
// calling an unspecified P(in*omega) for that exact zero would add no evidence.
fn projected_residual<P: OrbitProblem>(
    problem: &P,
    orbit: &HbOrbit,
    pack: &Packing,
    aft: &mut AftWorkspace,
) -> Vec<f64> {
    let retained = orbit.coeffs[0].len() - 1;
    for j in 0..aft.samples {
        for k in 0..pack.d {
            let mut x = orbit.coeffs[k][0].re;
            for n in 1..=retained {
                let c = orbit.coeffs[k][n];
                let i = j * pack.n_harm + n - 1;
                x += 2.0 * (c.re * aft.cosine[i] - c.im * aft.sine[i]);
            }
            aft.xs[j * pack.d + k] = x;
        }
        aft.island.fill(0.0);
        problem.island(
            TAU / orbit.omega * j as f64 / aft.samples as f64,
            &aft.xs[j * pack.d..(j + 1) * pack.d],
            &mut aft.island,
        );
        aft.gs[j * pack.d..(j + 1) * pack.d].copy_from_slice(&aft.island);
    }
    let mut residual = vec![0.0; pack.coeff_len()];
    for n in 0..=pack.n_harm {
        let port = (n <= retained).then(|| problem.port(C64::new(0.0, n as f64 * orbit.omega)));
        for k in 0..pack.d {
            let mut g = C64::new(0.0, 0.0);
            for j in 0..aft.samples {
                let basis = if n == 0 {
                    C64::new(1.0, 0.0)
                } else {
                    let i = j * pack.n_harm + n - 1;
                    C64::new(aft.cosine[i], -aft.sine[i])
                };
                g = g + basis.scale(aft.gs[j * pack.d + k] / aft.samples as f64);
            }
            let mut r = C64::new(0.0, 0.0) - g;
            if let Some(p) = &port {
                for kk in 0..pack.d {
                    r = r + p[k * pack.d + kk] * orbit.coeffs[kk][n];
                }
            }
            residual[pack.re_index(n, k)] = r.re;
            if n > 0 {
                residual[pack.im_index(n, k)] = r.im;
            }
        }
    }
    residual
}

fn assess<P: OrbitProblem>(
    problem: &P,
    orbit: &HbOrbit,
    budget: &HbRefinementBudget,
) -> Result<HbRefinementCheck, OrbitError> {
    let n = orbit.coeffs[0].len() - 1;
    let pack = Packing {
        d: orbit.coeffs.len(),
        n_harm: 2 * n + 1,
        omega_slot: None,
    };
    let mut aft = AftWorkspace::new(&pack)?;
    if aft.samples < orbit.aft_samples {
        aft = AftWorkspace::with_samples(&pack, orbit.aft_samples, 0)?;
    }
    let mut coarse = projected_residual(problem, orbit, &pack, &mut aft);
    loop {
        let mut finer = aft.finer(&pack)?;
        let fine = projected_residual(problem, orbit, &pack, &mut finer);
        let difference: Vec<_> = fine.iter().zip(&coarse).map(|(a, b)| a - b).collect();
        let sampling = row_ratio(&difference, &pack, 0, &budget.residual_tolerances);
        if sampling <= 1.0 {
            let total = row_ratio(&fine, &pack, 0, &budget.residual_tolerances);
            let omitted = row_ratio(&fine, &pack, n + 1, &budget.residual_tolerances);
            let mut retained = fine;
            for mode in n + 1..=pack.n_harm {
                for k in 0..pack.d {
                    retained[pack.re_index(mode, k)] = 0.0;
                    retained[pack.im_index(mode, k)] = 0.0;
                }
            }
            return Ok(HbRefinementCheck {
                harmonics: n,
                checked_harmonics: pack.n_harm,
                aft_samples: finer.samples,
                total_residual_ratio: total,
                retained_residual_ratio: row_ratio(
                    &retained,
                    &pack,
                    0,
                    &budget.residual_tolerances,
                ),
                omitted_residual_ratio: omitted,
                sampling_defect_ratio: sampling,
                state_change_ratio: None,
                relative_frequency_change: None,
            });
        }
        if finer.refinement >= super::MAX_AFT_REFINEMENTS {
            return Err(OrbitError::SamplingUnresolved {
                samples: finer.samples,
                residual: super::norm(&fine),
                defect: super::norm(&difference),
            });
        }
        aft = finer;
        coarse = fine;
    }
}

fn changes(
    old: &HbOrbit,
    new: &HbOrbit,
    anchor: HbAnchor,
    budget: &HbRefinementBudget,
) -> (f64, f64) {
    let mut phase = C64::new(1.0, 0.0);
    if matches!(anchor, HbAnchor::Backbone { .. }) {
        let a = old.coeffs[0][1];
        let b = new.coeffs[0][1];
        phase = (a / C64::new(a.abs(), 0.0)) * (b.conj() / C64::new(b.abs(), 0.0));
    }
    let pack = Packing {
        d: new.coeffs.len(),
        n_harm: new.coeffs[0].len() - 1,
        omega_slot: None,
    };
    let mut difference = vec![0.0; pack.coeff_len()];
    for k in 0..pack.d {
        difference[pack.re_index(0, k)] = new.coeffs[k][0].re - old.coeffs[k][0].re;
        let mut rotation = C64::new(1.0, 0.0);
        for n in 1..=pack.n_harm {
            rotation = rotation * phase;
            let previous = old.coeffs[k].get(n).copied().unwrap_or(C64::new(0.0, 0.0));
            let delta = new.coeffs[k][n] * rotation - previous;
            difference[pack.re_index(n, k)] = delta.re;
            difference[pack.im_index(n, k)] = delta.im;
        }
    }
    (
        row_ratio(&difference, &pack, 0, &budget.state_tolerances),
        (new.omega - old.omega).abs() / new.omega.abs().max(old.omega.abs()),
    )
}

fn refine<P: OrbitProblem>(
    problem: &P,
    anchor: HbAnchor,
    mut orbit: HbOrbit,
    hb: &HbBudget,
    budget: &HbRefinementBudget,
    scales: Option<&HbScales>,
) -> Result<AdaptiveHbOrbit, OrbitError> {
    let mut previous: Option<HbOrbit> = None;
    let mut trace = Vec::new();
    loop {
        if !orbit.omega.is_finite() || orbit.omega <= 0.0 {
            return Err(OrbitError::BadParameter {
                what: "adaptive HB solved frequency must be positive and finite",
            });
        }
        let mut check = assess(problem, &orbit, budget)?;
        if let Some(old) = &previous {
            let (state, frequency) = changes(old, &orbit, anchor, budget);
            check.state_change_ratio = Some(state);
            check.relative_frequency_change = Some(frequency);
        }
        trace.push(check);
        if check.retained_residual_ratio > 1.0 {
            return Err(OrbitError::HarmonicsUnresolved { check });
        }
        if check.total_residual_ratio <= 1.0
            && check.state_change_ratio.is_some_and(|ratio| ratio <= 1.0)
            && check
                .relative_frequency_change
                .is_some_and(|change| change <= budget.frequency_tolerance)
        {
            return Ok(AdaptiveHbOrbit {
                orbit,
                refinement_trace: trace,
            });
        }
        if check.harmonics == budget.max_harmonics {
            return Err(OrbitError::HarmonicsUnresolved { check });
        }
        let next = (2 * check.harmonics + 1).min(budget.max_harmonics);
        let pack = Packing {
            d: problem.dim(),
            n_harm: next,
            omega_slot: (!matches!(anchor, HbAnchor::Forced { .. }))
                .then_some(problem.dim() * (2 * next + 1)),
        };
        let mut seed = vec![0.0; pack.len()];
        for k in 0..pack.d {
            seed[pack.re_index(0, k)] = orbit.coeffs[k][0].re;
            for n in 1..=check.harmonics {
                seed[pack.re_index(n, k)] = orbit.coeffs[k][n].re;
                seed[pack.im_index(n, k)] = orbit.coeffs[k][n].im;
            }
        }
        if let Some(slot) = pack.omega_slot {
            seed[slot] = orbit.omega;
        }
        let next_budget = HbBudget {
            harmonics: next,
            ..*hb
        };
        let next_orbit = if let Some(scales) = scales {
            solve_hb_scaled_seeded(problem, anchor, seed, &next_budget, scales)?.orbit
        } else {
            solve_hb_seeded(problem, anchor, seed, &next_budget)?
        };
        previous = Some(orbit);
        orbit = next_orbit;
    }
}
