//! Harmonic refinement with a fixed physical amplitude and full-waveform phase.

use super::{
    BackboneContinuationBudget, BackboneContinuationError, BackboneContinuationPoint, Chart,
    Checked, Constraint, Evaluation, Result, correct_with_phase, invalid, poll,
};
use crate::{
    HbAnchor, HbBudget, HbOrbit, HbRefinementBudget, HbRefinementCheck, HbScales,
    MAX_AFT_REFINEMENTS, OrbitError, OrbitProblem, adaptive,
};
use fs_math::c64::C64;

/// Physical-chart corrector controls shared across every solved harmonic order.
#[derive(Clone, Debug)]
pub struct BackboneRefinementControls {
    /// Newton iterations per order, in 1..=1024.
    pub max_newton: usize,
    /// Absolute norm limit for normalized original balances, phase and relative
    /// amplitude together. A separate physical residual audit is also mandatory.
    pub tolerance: f64,
    /// Minimum non-DC waveform norm in the declared normalized state chart.
    pub minimum_oscillation: f64,
    /// Complete real Fourier unknowns plus frequency at the largest order.
    /// Also limited to 2048 by the original dense conservative corrector.
    pub maximum_unknowns: usize,
    /// Conservative numerical payload, including correction, independent
    /// assessment, both inter-order waveforms and fixed phase. Excludes the
    /// caller-owned source, initial seed and allocator metadata.
    pub maximum_payload_bytes: usize,
    /// Complete original AFT residual grids across ALL orders, including
    /// corrector probes, rejected trials and independent harmonic assessments.
    pub maximum_residual_evaluations: usize,
    /// Maximum estimated condition of the column-equilibrated QR R factor.
    /// This numerical screen is neither a certified rank nor forward-error bound.
    pub maximum_condition_estimate: f64,
}

/// An original physical backbone point resolved at successive harmonic orders.
#[derive(Clone, Debug)]
pub struct AdaptiveBackboneOrbit {
    /// The final rank-admitted conservative corrector result. Its orbit retains
    /// every original state and physical coefficient; its step is zero.
    pub point: BackboneContinuationPoint,
    /// Independent original-equation and actual inter-order measurements.
    pub refinement_trace: Vec<HbRefinementCheck>,
    pub residual_evaluations: usize,
    pub tangent_preparations: usize,
}

impl BackboneRefinementControls {
    fn corrector(&self, harmonics: usize, scales: &HbScales) -> BackboneContinuationBudget {
        BackboneContinuationBudget {
            harmonics,
            steps: 0,
            max_newton: self.max_newton,
            max_attempts_per_step: 1,
            step_sizes: [1.0; 3],
            tolerance: self.tolerance,
            minimum_oscillation: self.minimum_oscillation,
            maximum_unknowns: self.maximum_unknowns,
            maximum_payload_bytes: self.maximum_payload_bytes,
            maximum_residual_evaluations: self.maximum_residual_evaluations,
            maximum_condition_estimate: self.maximum_condition_estimate,
            scales: scales.clone(),
        }
    }
}

// Exact coefficient prolongation: never resample, clip or synthesize a solved
// high-order waveform. Newly admitted harmonics start at zero and are corrected.
fn prolong(orbit: &HbOrbit, harmonics: usize) -> HbOrbit {
    let mut coeffs = orbit.coeffs.clone();
    for row in &mut coeffs {
        row.resize(harmonics + 1, C64::ZERO);
    }
    // Supplied convergence metadata is neither trusted nor copied. In
    // particular an arbitrary caller trace cannot consume unbudgeted memory.
    HbOrbit {
        coeffs,
        omega: orbit.omega,
        residual: 0.0,
        residual_trace: Vec::new(),
        aft_samples: 0,
        sampling_defect: 0.0,
    }
}

// Bound every possible assessment grid, including a corrector that already
// needed its finest grid. Exact-length FFT plan sizes need not be monotone, so
// inspect the finite grid ladder instead of estimating it at one endpoint.
fn assessment_payload(d: usize, first: usize, maximum: usize) -> Result<usize> {
    let mut peak = 0_u128;
    let mut h = first;
    loop {
        let checked = 2 * h + 1;
        let nc = d as u128 * (2 * checked + 1) as u128;
        let mut corrector_samples = (4 * h + 4).next_power_of_two();
        for _ in 0..=MAX_AFT_REFINEMENTS {
            let mut samples = (4 * checked + 4).next_power_of_two().max(corrector_samples);
            for _ in 0..=MAX_AFT_REFINEMENTS {
                let fft = crate::transform::payload(samples, checked, d)? as u128;
                // Two phase/sample/transform workspaces, full/omitted residual
                // arrays, linear action and a complete original port matrix.
                peak = peak.max(
                    2 * fft
                        + 8 * (4 * samples as u128 * (d as u128 + 1)
                            + 8 * nc
                            + 16 * d as u128
                            + 2 * (d as u128).pow(2)),
                );
                samples = 2 * samples + 1;
            }
            corrector_samples = 2 * corrector_samples + 1;
        }
        if h == maximum {
            break;
        }
        h = (2 * h + 1).min(maximum);
    }
    // Persistent actual waveforms, reference, scaled chart and trace entries.
    peak += 128 * (d as u128 * (2 * maximum + 1) as u128 + 1);
    usize::try_from(peak)
        .ok()
        .filter(|&v| v <= isize::MAX as usize)
        .ok_or_else(|| invalid("adaptive backbone assessment payload overflows"))
}

/// Resolve a conservative periodic point with measured harmonic convergence.
///
/// The seed supplies every original coefficient and unknown frequency; its
/// convergence metadata is ignored. Its COMPLETE waveform fixes one Parseval
/// phase hyperplane, zero-extended unchanged at every larger order. The chosen
/// original coordinate pins `2*abs(X_1)` to the same positive physical amplitude.
/// No state permutation or phase rotation between solved orders is introduced.
///
/// Each solve reuses the conservative continuation owner's original-row
/// rectangular Householder QR corrector. Rank/condition admission is required
/// even for residual-zero seeds. It appends phase and amplitude without dropping
/// a physical balance, adding jitter, or removing dissipation. An exact declared
/// state tangent is used without fallback; otherwise bounded centered probes
/// are retained. Frequency follows that owner's centered original-port probes.
/// A stalled corrector independently checks a finer grid: a measured sampling
/// defect can trigger quadrature refinement at the SAME physical iterate, within
/// the existing Newton/grid budgets, before increasing retained harmonic order.
///
/// Orders grow as `min(2H+1,max_harmonics)`. New state harmonics start at zero;
/// the actual nonlinear equations solve for them. The shared adaptive assessor
/// independently projects through `2H+1` AND measures the complete sampled
/// original time-domain equation defect on distinct grids. Acceptance requires
/// physical retained/omitted/full residual tolerances, sampling agreement, and
/// waveform/frequency agreement between two actual solved orders in the SAME
/// fixed phase. Large DC states or unrelated amplitudes cannot relax these
/// absolute per-equation/state tolerances. A retained-row failure refuses rather
/// than treating spectral growth as a substitute for sufficient Newton accuracy.
///
/// Payload at every possible order is admitted before source callbacks. The
/// shared residual-grid allowance spans all orders and assessments. Cancellation
/// is polled around original AFT/tangent/QR passes and between assessment grids;
/// each bounded complete numerical pass remains atomic. Refusal publishes no
/// partial result. Measured finite resolution is not an infinite-tail bound,
/// continuous-time existence proof, primitive-period test or stability claim.
#[allow(clippy::too_many_arguments)]
pub fn solve_backbone_adaptive<P: OrbitProblem>(
    problem: &P,
    coordinate: usize,
    amplitude: f64,
    seed: &HbOrbit,
    controls: &BackboneRefinementControls,
    refinement: &HbRefinementBudget,
    scales: &HbScales,
    mut cancelled: impl FnMut() -> bool,
) -> Result<AdaptiveBackboneOrbit> {
    poll(&mut cancelled)?;
    let d = problem.dim();
    let first = seed
        .coeffs
        .first()
        .map_or(0, |row| row.len().saturating_sub(1));
    let anchor = HbAnchor::Backbone {
        amplitude,
        omega_guess: seed.omega,
    };
    let hb = HbBudget {
        harmonics: first,
        max_newton: controls.max_newton,
        tolerance: controls.tolerance,
        deflate_equilibrium: false,
    };
    adaptive::admit(d, anchor, &hb, refinement)?;
    let maximum_unknowns = d
        .checked_mul(2 * refinement.max_harmonics + 1)
        .and_then(|v| v.checked_add(1))
        .ok_or_else(|| invalid("unknown extent overflow"))?;
    if !problem.autonomous()
        || coordinate >= d
        || maximum_unknowns > controls.maximum_unknowns
        || maximum_unknowns > 2048
    {
        return Err(invalid(
            "complete autonomous backbone and admitted maximum harmonic extent required",
        ));
    }
    let callback = problem.island_state_jacobian();
    let extra = assessment_payload(d, first, refinement.max_harmonics)?;
    let available = controls
        .maximum_payload_bytes
        .checked_sub(extra)
        .ok_or_else(|| invalid("adaptive backbone assessment exceeds payload allowance"))?;
    let mut budget = controls.corrector(first, scales);
    budget.maximum_payload_bytes = available;
    // Validate the original seed BEFORE prolongation can conceal malformed rows.
    let initial_chart = Chart::new(
        d,
        coordinate,
        &[(amplitude, seed)],
        callback.is_some(),
        &budget,
    )?;
    let initial = initial_chart.seed(seed)?;
    let phase = initial_chart.phase(problem, &initial, &budget)?;
    let largest = prolong(seed, refinement.max_harmonics);
    let mut maximum_budget = budget.clone();
    maximum_budget.harmonics = refinement.max_harmonics;
    Chart::new(
        d,
        coordinate,
        &[(amplitude, &largest)],
        callback.is_some(),
        &maximum_budget,
    )?;
    drop(largest);
    drop(initial_chart);
    drop(initial);

    let mut candidate = prolong(seed, first);
    let mut previous: Option<HbOrbit> = None;
    let mut trace = Vec::new();
    let mut evaluations = 0;
    let mut tangent_preparations = 0;
    loop {
        poll(&mut cancelled)?;
        let h = candidate.coeffs[0].len() - 1;
        budget.harmonics = h;
        let chart = Chart::new(
            d,
            coordinate,
            &[(amplitude, &candidate)],
            callback.is_some(),
            &budget,
        )?;
        let predictor = chart.seed(&candidate)?;
        let mut eval = Evaluation {
            problem,
            chart: &chart,
            cancelled: &mut cancelled,
            evaluations,
            limit: controls.maximum_residual_evaluations,
            physical: vec![0.0; chart.pack.len()],
            callback,
            tangent: None,
            tangent_preparations,
        };
        let (_, point) = correct_with_phase(
            &mut eval,
            &predictor,
            Constraint::Amplitude(amplitude),
            &budget,
            0.0,
            &phase,
            true,
        )?;
        evaluations = eval.evaluations;
        tangent_preparations = eval.tangent_preparations;
        drop(eval);
        let mut check = adaptive::assess_with(
            &Checked {
                source: problem,
                dimension: d,
            },
            &point.orbit,
            refinement,
            |_samples| {
                poll(&mut cancelled)?;
                if !problem.autonomous() || problem.dim() != d {
                    return Err(invalid(
                        "source must preserve its complete autonomous state chart",
                    ));
                }
                if evaluations >= controls.maximum_residual_evaluations {
                    return Err(BackboneContinuationError::EvaluationLimit { evaluations });
                }
                evaluations += 1;
                Ok(())
            },
        )?;
        poll(&mut cancelled)?;
        if let Some(old) = &previous {
            // Autonomous here selects NO phase alignment: both actual solves
            // already satisfy the one immutable complete-waveform phase plane.
            let (state, frequency) = adaptive::changes(
                old,
                &point.orbit,
                HbAnchor::Autonomous {
                    omega_guess: seed.omega,
                },
                refinement,
            );
            check.state_change_ratio = Some(state);
            check.relative_frequency_change = Some(frequency);
        }
        trace.push(check);
        if check.retained_residual_ratio > 1.0 {
            return Err(OrbitError::HarmonicsUnresolved { check }.into());
        }
        if check.total_residual_ratio <= 1.0
            && check.omitted_residual_ratio <= 1.0
            && check.state_change_ratio.is_some_and(|v| v <= 1.0)
            && check
                .relative_frequency_change
                .is_some_and(|v| v <= refinement.frequency_tolerance)
        {
            return Ok(AdaptiveBackboneOrbit {
                point,
                refinement_trace: trace,
                residual_evaluations: evaluations,
                tangent_preparations,
            });
        }
        if h == refinement.max_harmonics {
            return Err(OrbitError::HarmonicsUnresolved { check }.into());
        }
        candidate = prolong(&point.orbit, (2 * h + 1).min(refinement.max_harmonics));
        previous = Some(point.orbit);
    }
}
