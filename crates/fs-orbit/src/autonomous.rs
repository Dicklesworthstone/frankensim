//! Autonomous periodic branches with independent frequency and parameter.
//!
//! The original AFT balance and dense LU remain the numerical owners. The two
//! extra rows are a full-waveform integral phase condition and a scaled
//! pseudo-arclength plane. No physical equation, damping or memory is removed.

use super::{AftWorkspace, ContinuableProblem, HbOrbit, OrbitError, Packing,
    MAX_AFT_REFINEMENTS, hb_residual, pack_orbit, HbAnchor};
use fs_la::factor::lu;
use fs_math::det;

/// Physical units used only to condition the continuation equations and metric.
#[derive(Debug, Clone)]
pub struct AutonomousScales {
    /// Positive characteristic magnitudes in each original state coordinate.
    pub states: Vec<f64>,
    /// Positive characteristic magnitudes in each original balance row's units.
    /// These need not have state/time units for a general frequency-domain port.
    pub equations: Vec<f64>,
    pub angular_frequency: f64,
    pub parameter: f64,
}

/// Bounded cold continuation; all tolerances below refer to the declared scales.
#[derive(Debug, Clone)]
pub struct AutonomousContinuationBudget {
    /// Fixed retained order, including all explicit source and memory states.
    pub harmonics: usize,
    /// Number of additional accepted points after the two corrected seeds (at most 65,536).
    pub steps: usize,
    /// Corrector iteration cap in 1..=1024.
    pub max_newton: usize,
    /// Maximum attempts per additional point in 1..=64, including halving retries.
    pub max_attempts_per_step: usize,
    /// Initial/minimum/maximum positive dimensionless pseudo-arclength steps.
    pub step_sizes: [f64; 3],
    /// Absolute norm limit on scaled HB, phase and continuation rows together.
    /// DC magnitude does not relax this limit. Sampling has the same limit.
    pub tolerance: f64,
    /// Minimum scaled RMS non-DC waveform norm; equilibrium is not a cycle.
    pub minimum_oscillation: f64,
    /// Complete coefficients + frequency + parameter, also hard-capped at 2048.
    pub maximum_unknowns: usize,
    /// Conservative solver/AFT/retained-result payload; excludes caller models
    /// and allocator metadata. Checked before source evaluation or allocation.
    pub maximum_payload_bytes: usize,
    /// Total original AFT residual calls, including Jacobians and finer grids.
    pub maximum_residual_evaluations: usize,
    pub scales: AutonomousScales,
}

/// One accepted point, retaining physical coefficients and measured diagnostics.
#[derive(Debug, Clone)]
pub struct AutonomousContinuationPoint {
    pub parameter: f64,
    /// Original-state Fourier waveform and physical frequency. The native orbit
    /// residual/trace/defect fields contain unscaled HB norms, not mixed units.
    pub orbit: HbOrbit,
    /// Zero for either fixed-parameter seed; otherwise the actual accepted step.
    pub step: f64,
    pub scaled_residual: f64,
    pub scaled_sampling_defect: f64,
    pub phase_residual: f64,
    pub arclength_residual: f64,
}

#[derive(Debug, Clone)]
pub struct AutonomousBranch {
    pub points: Vec<AutonomousContinuationPoint>,
    pub residual_evaluations: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub enum AutonomousContinuationError {
    Invalid(&'static str),
    Cancelled,
    EvaluationLimit { evaluations: usize },
    Numerical(OrbitError),
}
impl core::fmt::Display for AutonomousContinuationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "autonomous continuation refused: {self:?}")
    }
}
impl std::error::Error for AutonomousContinuationError {}
impl From<OrbitError> for AutonomousContinuationError {
    fn from(error: OrbitError) -> Self { Self::Numerical(error) }
}
type Result<T> = core::result::Result<T, AutonomousContinuationError>;
fn invalid(what: &'static str) -> AutonomousContinuationError {
    AutonomousContinuationError::Invalid(what)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() { Err(AutonomousContinuationError::Cancelled) } else { Ok(()) }
}
fn length(values: &[f64]) -> f64 {
    values.iter().fold(0.0, |a, &b| det::hypot(a, b))
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    let (mut sum, mut correction) = (0.0, 0.0);
    for (&a, &b) in a.iter().zip(b) {
        let term = a * b - correction;
        let next = sum + term;
        correction = (next - sum) - term;
        sum = next;
    }
    sum
}

struct Chart {
    pack: Packing,
    /// physical packed coefficient = normalized coefficient * state_unit
    state_units: Vec<f64>,
    row_units: Vec<f64>,
}
impl Chart {
    fn new(d: usize, seeds: [(f64, &HbOrbit); 2], b: &AutonomousContinuationBudget) -> Result<Self> {
        let positive = |v: f64| v.is_finite() && v > 0.0;
        let nc = b.harmonics.checked_mul(2).and_then(|n| n.checked_add(1))
            .and_then(|n| n.checked_mul(d)).ok_or_else(|| invalid("harmonic extent overflow"))?;
        let nz = nc.checked_add(2).ok_or_else(|| invalid("unknown extent overflow"))?;
        if d == 0 || b.harmonics == 0 || nz > 2048 || nz > b.maximum_unknowns
            || b.steps > 65_536 || !(1..=1024).contains(&b.max_newton)
            || !(1..=64).contains(&b.max_attempts_per_step)
            || !positive(b.tolerance) || b.tolerance >= 1.0
            || !positive(b.minimum_oscillation)
            || b.step_sizes.iter().any(|&v| !positive(v))
            || b.step_sizes[1] > b.step_sizes[0] || b.step_sizes[0] > b.step_sizes[2]
            || !positive(b.scales.angular_frequency) || !positive(b.scales.parameter)
            || b.scales.states.len() != d || b.scales.equations.len() != d
            || b.scales.states.iter().chain(&b.scales.equations).any(|&v| !positive(v))
            || b.maximum_residual_evaluations == 0
            || !seeds[0].0.is_finite() || !seeds[1].0.is_finite() || seeds[0].0 == seeds[1].0
        { return Err(invalid("complete finite seeds, scales and continuation budgets required")); }
        for (_, seed) in seeds {
            if !positive(seed.omega) || !(super::TAU / seed.omega).is_finite()
                || seed.coeffs.len() != d || seed.coeffs.iter().any(|row|
                    row.len() != b.harmonics + 1 || row[0].im != 0.0
                    || row.iter().any(|v| !v.re.is_finite() || !v.im.is_finite()))
            { return Err(invalid("every original-state seed coefficient must match the declared order")); }
        }
        let samples = (4 * b.harmonics + 4).next_power_of_two() as u128;
        let finest = (samples + 1) * (1_u128 << MAX_AFT_REFINEMENTS) - 1;
        let count = b.steps as u128 + 2;
        let size = nz as u128;
        let bytes = 8 * (3 * size * size + 64 * size
            + 4 * finest * (d + b.harmonics) as u128
            + count * (2 * d as u128 * (b.harmonics + 1) as u128
                + 3 * d as u128 + b.max_newton as u128 + 32));
        if bytes > b.maximum_payload_bytes as u128 {
            return Err(invalid("complete continuation payload exceeds its allowance"));
        }
        let pack = Packing { d, n_harm: b.harmonics, omega_slot: Some(nc) };
        let mut state_units = vec![0.0; nz];
        let mut row_units = vec![0.0; nc];
        for k in 0..d {
            state_units[k] = b.scales.states[k];
            row_units[k] = b.scales.equations[k];
            for n in 1..=b.harmonics {
                for slot in [pack.re_index(n, k), pack.im_index(n, k)] {
                    state_units[slot] = b.scales.states[k] / core::f64::consts::SQRT_2;
                    row_units[slot] = b.scales.equations[k] / core::f64::consts::SQRT_2;
                }
            }
        }
        state_units[nc] = b.scales.angular_frequency;
        state_units[nc + 1] = b.scales.parameter;
        if state_units.iter().chain(&row_units).any(|&v| !positive(v)) {
            return Err(invalid("physical normalization is not representable"));
        }
        Ok(Self { pack, state_units, row_units })
    }
    fn seed(&self, parameter: f64, orbit: &HbOrbit) -> Result<Vec<f64>> {
        let mut physical = pack_orbit(orbit, &HbAnchor::Autonomous { omega_guess: orbit.omega });
        physical.push(parameter);
        for (v, unit) in physical.iter_mut().zip(&self.state_units) {
            let old = *v;
            *v /= unit;
            if !v.is_finite() || (old != 0.0 && *v == 0.0) {
                return Err(invalid("normalized seed is not representable"));
            }
        }
        Ok(physical)
    }
    fn phase(&self, reference: &[f64], floor: f64) -> Result<Vec<f64>> {
        let nc = self.pack.coeff_len();
        if length(&reference[self.pack.d..nc]) <= floor {
            return Err(OrbitError::TrivialCollapse.into());
        }
        let mut phase = vec![0.0; nc + 2];
        for n in 1..=self.pack.n_harm {
            for k in 0..self.pack.d {
                let (re, im) = (self.pack.re_index(n, k), self.pack.im_index(n, k));
                phase[re] = -(n as f64) * reference[im];
                phase[im] = n as f64 * reference[re];
            }
        }
        let norm = length(&phase);
        if !norm.is_finite() || norm == 0.0 { return Err(invalid("unresolved phase direction")); }
        for v in &mut phase { *v /= norm; }
        Ok(phase)
    }
    fn scaled(&self, raw: &[f64]) -> Vec<f64> {
        raw.iter().zip(&self.row_units).map(|(v, unit)| v / unit).collect()
    }
}

struct Evaluation<'a, P, F> {
    problem: &'a mut P,
    chart: &'a Chart,
    cancelled: &'a mut F,
    evaluations: usize,
    limit: usize,
    physical: Vec<f64>,
}
impl<P: ContinuableProblem, F: FnMut() -> bool> Evaluation<'_, P, F> {
    fn eval(&mut self, z: &[f64], aft: &mut AftWorkspace, out: &mut [f64]) -> Result<()> {
        poll(self.cancelled)?;
        if self.evaluations >= self.limit {
            return Err(AutonomousContinuationError::EvaluationLimit { evaluations: self.evaluations });
        }
        let nc = self.chart.pack.coeff_len();
        for ((target, &v), &unit) in self.physical.iter_mut().zip(z).zip(&self.chart.state_units) {
            *target = v * unit;
            if !target.is_finite() || (v != 0.0 && *target == 0.0) {
                return Err(invalid("physical trial coordinates are not representable"));
            }
        }
        if self.physical[nc] <= 0.0 || !(super::TAU / self.physical[nc]).is_finite() {
            return Err(invalid("physical trial period must be positive and finite"));
        }
        self.problem.set_parameter(self.physical[nc + 1]);
        if self.problem.dim() != self.chart.pack.d || !self.problem.autonomous() {
            return Err(invalid("continuation must preserve the complete autonomous state chart"));
        }
        self.evaluations += 1;
        hb_residual(self.problem, &self.chart.pack, &self.physical[..nc + 1],
            f64::NAN, aft, out);
        poll(self.cancelled)?;
        if out.iter().any(|v| !v.is_finite()) {
            return Err(invalid("original harmonic balance produced nonfinite values"));
        }
        Ok(())
    }
}

fn augmented(chart: &Chart, raw: &[f64], z: &[f64], reference: &[f64],
    phase: &[f64], predictor: &[f64], plane: &[f64]) -> Vec<f64>
{
    let mut r = chart.scaled(raw);
    r.push(dot(&z.iter().zip(reference).map(|(a, b)| a - b).collect::<Vec<_>>(), phase));
    r.push(dot(&z.iter().zip(predictor).map(|(a, b)| a - b).collect::<Vec<_>>(), plane));
    r
}

#[allow(clippy::too_many_arguments)]
fn correct<P: ContinuableProblem, F: FnMut() -> bool>(
    eval: &mut Evaluation<'_, P, F>, reference: &[f64], predictor: &[f64], plane: &[f64],
    b: &AutonomousContinuationBudget, step: f64,
) -> Result<(Vec<f64>, AutonomousContinuationPoint)> {
    let chart = eval.chart;
    let nc = chart.pack.coeff_len();
    let nz = nc + 2;
    let phase = chart.phase(reference, b.minimum_oscillation)?;
    let mut z = predictor.to_vec();
    let mut aft = AftWorkspace::new(&chart.pack)?;
    let mut raw = vec![0.0; nc];
    let mut plus = vec![0.0; nc];
    let mut minus = vec![0.0; nc];
    let mut jac = vec![0.0; nz * nz];
    let mut trace = Vec::with_capacity(b.max_newton);
    for _ in 0..b.max_newton {
        eval.eval(&z, &mut aft, &mut raw)?;
        let r = augmented(chart, &raw, &z, reference, &phase, predictor, plane);
        let rn = length(&r);
        trace.push(length(&raw));
        if !rn.is_finite() { return Err(invalid("scaled corrector residual overflow")); }
        if rn <= b.tolerance {
            let mut finer = aft.finer(&chart.pack)?;
            eval.eval(&z, &mut finer, &mut plus)?;
            let checked = augmented(chart, &plus, &z, reference, &phase, predictor, plane);
            let difference: Vec<_> = plus.iter().zip(&raw).map(|(a, b)| a - b).collect();
            let defect = length(&chart.scaled(&difference));
            let checked_norm = length(&checked);
            if checked_norm <= b.tolerance && defect <= b.tolerance {
                let _ = chart.phase(&z, b.minimum_oscillation)?;
                let (coeffs, omega) = chart.pack.unpack(&eval.physical[..nc + 1]);
                *trace.last_mut().expect("current corrector") = length(&plus);
                return Ok((z, AutonomousContinuationPoint {
                    parameter: eval.physical[nc + 1], step,
                    orbit: HbOrbit { omega, coeffs, residual: length(&plus),
                        residual_trace: trace, aft_samples: finer.samples,
                        sampling_defect: length(&difference) },
                    scaled_residual: checked_norm, scaled_sampling_defect: defect,
                    phase_residual: checked[nc], arclength_residual: checked[nc + 1],
                }));
            }
            if finer.refinement >= MAX_AFT_REFINEMENTS {
                return Err(OrbitError::SamplingUnresolved { samples: finer.samples,
                    residual: length(&plus), defect: length(&difference) }.into());
            }
            aft = finer;
            continue;
        }
        for col in 0..nz {
            poll(eval.cancelled)?;
            let old = z[col];
            let h = if col == nc { (1e-6 * old.abs().max(1.0)).min(0.25 * old) }
                else { 1e-6 * old.abs().max(1.0) };
            if !h.is_finite() || h <= 0.0 || old + h == old || old - h == old {
                return Err(invalid("finite-difference column is not representable"));
            }
            z[col] = old + h;
            eval.eval(&z, &mut aft, &mut plus)?;
            z[col] = old - h;
            eval.eval(&z, &mut aft, &mut minus)?;
            z[col] = old;
            for row in 0..nc {
                jac[row * nz + col] = ((plus[row] - minus[row]) / (2.0 * h)) / chart.row_units[row];
            }
            jac[nc * nz + col] = phase[col];
            jac[(nc + 1) * nz + col] = plane[col];
        }
        if jac.iter().any(|v| !v.is_finite()) { return Err(invalid("corrector Jacobian overflow")); }
        let factor = lu(&jac, nz).map_err(|_| OrbitError::SingularJacobian)?;
        let mut delta = r;
        factor.solve(&mut delta);
        poll(eval.cancelled)?;
        let mut accepted = None;
        for backtrack in 0..12_i32 {
            let alpha = 2.0_f64.powi(-backtrack);
            let trial: Vec<_> = z.iter().zip(&delta).map(|(v, dv)| v - alpha * dv).collect();
            if trial[nc] <= 0.0 || trial.iter().any(|v| !v.is_finite()) { continue; }
            eval.eval(&trial, &mut aft, &mut plus)?;
            let trial_r = augmented(chart, &plus, &trial, reference, &phase, predictor, plane);
            if length(&trial_r) < rn { accepted = Some(trial); break; }
        }
        match accepted {
            Some(trial) => z = trial,
            None => return Err(OrbitError::NewtonStalled {
                residual: *trace.last().expect("iteration"), iterations: trace.len(), trace }.into()),
        }
    }
    Err(OrbitError::NewtonStalled { residual: *trace.last().unwrap_or(&f64::INFINITY),
        iterations: b.max_newton, trace }.into())
}

fn tangent(chart: &Chart, old: &[f64], now: &[f64], floor: f64) -> Result<Vec<f64>> {
    let phase = chart.phase(now, floor)?;
    let mut t: Vec<_> = now.iter().zip(old).map(|(a, b)| a - b).collect();
    let unprojected = length(&t);
    let phase_drift = dot(&t, &phase);
    for (v, p) in t.iter_mut().zip(&phase) { *v -= phase_drift * p; }
    let norm = length(&t);
    if !norm.is_finite() || norm <= 64.0 * f64::EPSILON * unprojected {
        return Err(invalid("seeds or accepted points do not resolve a branch tangent"));
    }
    for v in &mut t { *v /= norm; }
    Ok(t)
}

/// Continue smooth self-excited cycles through folds in a physical parameter.
///
/// Both nearby full-state seeds are corrected at their authored parameters;
/// convergence metadata supplied in the seeds is never trusted. Their order
/// sets the branch direction. Each later corrector solves coefficients,
/// frequency AND parameter jointly. The normalized integral phase condition
/// uses every non-DC state harmonic and is refreshed at each accepted point.
/// Parseval weights make the metric independent of exact-zero harmonic padding.
///
/// Success returns exactly `steps + 2` points, not a parameter endpoint (which
/// can be crossed several times at folds). On any refusal, the mutable family's
/// parameter is restored to the last accepted value, or the first seed value
/// before the first acceptance. No partial branch is published. Source callbacks
/// must be deterministic and may not hide evolving state in parameter setters.
/// Cancellation is polled between residuals, Jacobian columns and steps; native
/// AFT evaluation and a dense LU remain bounded noninterruptible units.
///
/// Fixed harmonic truncation and finite-difference Newton remain disclosed.
/// This neither classifies stability nor certifies folds, branch identity,
/// fundamental period, omitted Fourier tails, hybrid cycles or conservative
/// families with more than one neutral direction. Source-band limits belong to
/// the original model. No numerical damping, equilibrium deflation or Jacobian
/// jitter is introduced to force a branch through a physical/numerical refusal.
pub fn continue_autonomous_branch<P: ContinuableProblem>(
    problem: &mut P, seeds: [(f64, &HbOrbit); 2], budget: &AutonomousContinuationBudget,
    mut cancelled: impl FnMut() -> bool,
) -> Result<AutonomousBranch> {
    poll(&mut cancelled)?;
    let chart = Chart::new(problem.dim(), seeds, budget)?;
    if !problem.autonomous() { return Err(invalid("an autonomous problem is required")); }
    let seed0 = chart.seed(seeds[0].0, seeds[0].1)?;
    let seed1 = chart.seed(seeds[1].0, seeds[1].1)?;
    let nz = chart.state_units.len();
    let mut eval = Evaluation { problem, chart: &chart, cancelled: &mut cancelled,
        evaluations: 0, limit: budget.maximum_residual_evaluations, physical: vec![0.0; nz] };
    let mut last_parameter = seeds[0].0;
    let result = (|| {
        let mut fixed = vec![0.0; nz];
        fixed[nz - 1] = 1.0;
        let (a, first) = correct(&mut eval, &seed0, &seed0, &fixed, budget, 0.0)?;
        last_parameter = first.parameter;
        let (mut z, second) = correct(&mut eval, &a, &seed1, &fixed, budget, 0.0)?;
        last_parameter = second.parameter;
        let mut direction = tangent(&chart, &a, &z, budget.minimum_oscillation)?;
        let mut points = Vec::with_capacity(budget.steps + 2);
        points.push(first);
        points.push(second);
        let mut step = budget.step_sizes[0];
        for _ in 0..budget.steps {
            poll(eval.cancelled)?;
            let mut accepted = None;
            for _ in 0..budget.max_attempts_per_step {
                let predictor: Vec<_> = z.iter().zip(&direction).map(|(a, t)| a + step * t).collect();
                match correct(&mut eval, &z, &predictor, &direction, budget, step) {
                    Ok(point) => { accepted = Some(point); break; }
                    Err(AutonomousContinuationError::Numerical(OrbitError::NewtonStalled { .. }
                        | OrbitError::SingularJacobian | OrbitError::TrivialCollapse)) => {
                        if step <= budget.step_sizes[1] { break; }
                        step = (0.5 * step).max(budget.step_sizes[1]);
                    }
                    Err(error) => return Err(error),
                }
            }
            let Some((next, point)) = accepted else {
                return Err(OrbitError::ContinuationExhausted { steps: points.len() }.into());
            };
            let mut next_direction = tangent(&chart, &z, &next, budget.minimum_oscillation)?;
            if dot(&next_direction, &direction) < 0.0 {
                for v in &mut next_direction { *v = -*v; }
            }
            direction = next_direction;
            z = next;
            last_parameter = point.parameter;
            points.push(point);
            step = (step * 1.3).min(budget.step_sizes[2]);
        }
        poll(eval.cancelled)?;
        Ok(AutonomousBranch { points, residual_evaluations: eval.evaluations })
    })();
    eval.problem.set_parameter(last_parameter);
    result
}
