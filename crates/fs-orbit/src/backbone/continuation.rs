//! Conservative families retain both their phase and energy neutral directions.
//! Adding phase and arclength to ALL balance rows produces a rectangular
//! corrector, solved by the original Householder QR owner without regularization.

pub(crate) mod linear;

use crate::{
    AftSum, AftWorkspace, HbAnchor, HbBudget, HbOrbit, HbScales, IslandStateJacobian,
    MAX_AFT_REFINEMENTS, OrbitError, OrbitProblem, Packing, TAU, hb_residual, pack_orbit,
    phase::PhaseCondition, scaling::AdmittedScales, tangent::Workspace,
};
use fs_math::{c64::C64, det};

/// Bounded, fixed-order conservative continuation in the original physical units.
#[derive(Debug, Clone)]
pub struct BackboneContinuationBudget {
    /// Every source state retains DC and this many complex positive harmonics.
    pub harmonics: usize,
    /// Additional accepted points after the two corrected seeds; at most 65,536.
    pub steps: usize,
    /// Corrector iteration limit in 1..=1024.
    pub max_newton: usize,
    /// Attempts per additional point, including step-halving retries, in 1..=64.
    pub max_attempts_per_step: usize,
    /// Initial, minimum and maximum dimensionless Parseval arclength steps.
    pub step_sizes: [f64; 3],
    /// Absolute norm limit on normalized balance, phase and constraint together.
    /// Independent retained-row sampling disagreement has the same limit.
    pub tolerance: f64,
    /// Minimum normalized RMS non-DC waveform norm.
    pub minimum_oscillation: f64,
    /// All real Fourier unknowns plus angular frequency; also hard-capped at 2048.
    pub maximum_unknowns: usize,
    /// Conservative numerical workspace and retained-result payload, excluding
    /// caller-owned models and allocator metadata. Admitted before source work.
    pub maximum_payload_bytes: usize,
    /// Original AFT residual calls, including finite differences and finer grids.
    pub maximum_residual_evaluations: usize,
    /// Maximum estimated 1-norm condition of the column-equilibrated QR R factor.
    /// This finite numerical screen is not a rigorous rank/conditioning bound.
    pub maximum_condition_estimate: f64,
    pub scales: HbScales,
}

/// A physically corrected point with separately measured numerical evidence.
#[derive(Debug, Clone)]
pub struct BackboneContinuationPoint {
    /// Original state order, coefficients and frequency. Native residual fields
    /// contain only unscaled original HB norms, never phase or arclength units.
    pub orbit: HbOrbit,
    /// Peak of the selected coordinate's first harmonic, not a turning point.
    /// Later continuation points do not constrain this value or its direction.
    pub first_harmonic_amplitude: f64,
    /// Zero for seeds; otherwise the accepted normalized arclength step.
    pub step: f64,
    pub scaled_residual: f64,
    pub scaled_sampling_defect: f64,
    pub phase_residual: f64,
    /// Relative amplitude error for seeds, arclength-plane defect thereafter.
    pub constraint_residual: f64,
    pub jacobian_condition_estimate: f64,
}

#[derive(Debug, Clone)]
pub struct BackboneBranch {
    /// Original physical coordinate used to anchor the two seed amplitudes.
    pub coordinate: usize,
    pub points: Vec<BackboneContinuationPoint>,
    pub residual_evaluations: usize,
    /// Complete exact island tangent preparations, separate from residual calls.
    pub tangent_preparations: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub enum BackboneContinuationError {
    Invalid(&'static str),
    Cancelled,
    EvaluationLimit {
        evaluations: usize,
    },
    /// The equilibrated triangular factor exceeded the authored numerical limit.
    IllConditioned {
        estimate: f64,
        maximum: f64,
    },
    Numerical(OrbitError),
}
impl core::fmt::Display for BackboneContinuationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "backbone continuation refused: {self:?}")
    }
}
impl std::error::Error for BackboneContinuationError {}
impl From<OrbitError> for BackboneContinuationError {
    fn from(value: OrbitError) -> Self {
        Self::Numerical(value)
    }
}
type Result<T> = core::result::Result<T, BackboneContinuationError>;
fn invalid(what: &'static str) -> BackboneContinuationError {
    BackboneContinuationError::Invalid(what)
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(BackboneContinuationError::Cancelled)
    } else {
        Ok(())
    }
}
fn length(values: &[f64]) -> f64 {
    crate::scaling::length(values.iter().copied())
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    let mut sum = AftSum::new(0.0);
    for (&a, &b) in a.iter().zip(b) {
        sum.add(a * b);
    }
    sum.total()
}

struct Chart {
    pack: Packing,
    units: AdmittedScales,
    hb: HbBudget,
    coordinate: usize,
}
impl Chart {
    fn new(
        d: usize,
        coordinate: usize,
        seeds: [(f64, &HbOrbit); 2],
        exact: bool,
        b: &BackboneContinuationBudget,
    ) -> Result<Self> {
        let nc = b
            .harmonics
            .checked_mul(2)
            .and_then(|n| n.checked_add(1))
            .and_then(|n| n.checked_mul(d))
            .ok_or_else(|| invalid("harmonic extent overflow"))?;
        let nz = nc
            .checked_add(1)
            .ok_or_else(|| invalid("unknown extent overflow"))?;
        if d == 0
            || coordinate >= d
            || b.harmonics == 0
            || nz > 2048
            || nz > b.maximum_unknowns
            || b.steps > 65_536
            || !(1..=1024).contains(&b.max_newton)
            || !(1..=64).contains(&b.max_attempts_per_step)
            || !positive(b.tolerance)
            || b.tolerance >= 1.0
            || !positive(b.minimum_oscillation)
            || b.step_sizes.iter().any(|&v| !positive(v))
            || b.step_sizes[1] > b.step_sizes[0]
            || b.step_sizes[0] > b.step_sizes[2]
            || !b.maximum_condition_estimate.is_finite()
            || b.maximum_condition_estimate < 1.0
            || b.maximum_residual_evaluations == 0
            || seeds[0].0 == seeds[1].0
        {
            return Err(invalid(
                "finite complete seeds, coordinate and continuation budgets required",
            ));
        }
        for (amplitude, seed) in seeds {
            if !positive(amplitude)
                || !positive(0.5 * amplitude)
                || !positive(seed.omega)
                || !positive(TAU / seed.omega)
                || seed.coeffs.len() != d
                || seed.coeffs.iter().any(|row| {
                    row.len() != b.harmonics + 1
                        || row[0].im != 0.0
                        || row.iter().any(|z| !z.re.is_finite() || !z.im.is_finite())
                })
                || !positive(det::hypot(
                    seed.coeffs[coordinate][1].re,
                    seed.coeffs[coordinate][1].im,
                ))
            {
                return Err(invalid(
                    "seeds require every original coefficient and a selected fundamental",
                ));
            }
        }
        let samples = (4 * b.harmonics + 4).next_power_of_two();
        let finest = (samples as u128 + 1) * (1_u128 << MAX_AFT_REFINEMENTS) - 1;
        let mut fft_samples = samples;
        let mut fft_peak = 0;
        for refinement in 0..=MAX_AFT_REFINEMENTS {
            fft_peak = fft_peak.max(crate::transform::payload(fft_samples, b.harmonics, d)?);
            if refinement < MAX_AFT_REFINEMENTS {
                fft_samples = 2 * fft_samples + 1;
            }
        }
        let size = nz as u128;
        let square = d as u128 * d as u128;
        let orders = 2 * b.harmonics as u128 + 1;
        let tangent = if exact {
            4 * square * orders + 2 * square * (b.harmonics as u128 + 1) + 2 * orders + 2 * square
        } else {
            0
        };
        // Jacobian + QR copy, bounded phase/column/trial work, simultaneous AFT
        // grids/FFT plans, optional analytic tangent, and every retained orbit.
        let bytes = 2 * fft_peak as u128
            + 8 * (3 * (size + 1) * size
                + 80 * size
                + 4 * finest * (d + 1) as u128
                + tangent
                + (b.steps as u128 + 2)
                    * (2 * d as u128 * (b.harmonics + 1) as u128
                        + 3 * d as u128
                        + b.max_newton as u128
                        + 32));
        if bytes > b.maximum_payload_bytes as u128 || bytes > isize::MAX as u128 {
            return Err(invalid(
                "complete backbone continuation payload exceeds its allowance",
            ));
        }
        let hb = HbBudget {
            harmonics: b.harmonics,
            max_newton: b.max_newton,
            tolerance: b.tolerance,
            deflate_equilibrium: false,
        };
        let anchor = HbAnchor::Autonomous {
            omega_guess: seeds[0].1.omega,
        };
        let units =
            AdmittedScales::new(d, anchor, &pack_orbit(seeds[0].1, &anchor), &hb, &b.scales)?;
        Ok(Self {
            pack: Packing {
                d,
                n_harm: b.harmonics,
                omega_slot: Some(nc),
            },
            units,
            hb,
            coordinate,
        })
    }
    fn seed(&self, orbit: &HbOrbit) -> Result<Vec<f64>> {
        let mut z = pack_orbit(
            orbit,
            &HbAnchor::Autonomous {
                omega_guess: orbit.omega,
            },
        );
        for (v, &unit) in z.iter_mut().zip(&self.units.columns) {
            let original = *v;
            *v /= unit;
            if !v.is_finite() || (original != 0.0 && *v == 0.0) {
                return Err(invalid("normalized seed is not representable"));
            }
        }
        Ok(z)
    }
    fn physical(&self, z: &[f64], out: &mut [f64]) -> Result<()> {
        for ((v, &value), &unit) in out.iter_mut().zip(z).zip(&self.units.columns) {
            *v = value * unit;
            if !v.is_finite() || (value != 0.0 && *v == 0.0) {
                return Err(invalid("physical trial coordinates are not representable"));
            }
        }
        if !positive(out[self.pack.coeff_len()]) || !positive(TAU / out[self.pack.coeff_len()]) {
            return Err(invalid("physical trial period must be positive and finite"));
        }
        Ok(())
    }
    fn phase<P: OrbitProblem>(
        &self,
        problem: &P,
        z: &[f64],
        b: &BackboneContinuationBudget,
    ) -> Result<PhaseCondition> {
        if length(&z[self.pack.d..self.pack.coeff_len()]) <= b.minimum_oscillation {
            return Err(OrbitError::TrivialCollapse.into());
        }
        let mut physical = vec![0.0; z.len()];
        self.physical(z, &mut physical)?;
        Ok(PhaseCondition::new(
            problem, &physical, &self.hb, &b.scales,
        )?)
    }
}

// Incomplete source writes and malformed port dimensions must not become stale
// samples or indexing panics in the original AFT owner.
struct Checked<'a, P> {
    source: &'a P,
    dimension: usize,
}
impl<P: OrbitProblem> OrbitProblem for Checked<'_, P> {
    fn dim(&self) -> usize {
        self.dimension
    }
    fn autonomous(&self) -> bool {
        true
    }
    fn island(&self, t: f64, x: &[f64], out: &mut [f64]) {
        out.fill(f64::NAN);
        self.source.island(t, x, out);
    }
    fn port(&self, s: C64) -> Vec<C64> {
        let matrix = self.source.port(s);
        if matrix.len() == self.dimension * self.dimension
            && matrix.iter().all(|z| z.re.is_finite() && z.im.is_finite())
        {
            matrix
        } else {
            vec![C64::new(f64::NAN, f64::NAN); self.dimension * self.dimension]
        }
    }
}

struct Evaluation<'a, P, F> {
    problem: &'a P,
    chart: &'a Chart,
    cancelled: &'a mut F,
    evaluations: usize,
    limit: usize,
    physical: Vec<f64>,
    callback: Option<IslandStateJacobian<P>>,
    tangent: Option<Workspace>,
    tangent_preparations: usize,
}
impl<P: OrbitProblem, F: FnMut() -> bool> Evaluation<'_, P, F> {
    fn eval(&mut self, z: &[f64], aft: &mut AftWorkspace, out: &mut [f64]) -> Result<()> {
        poll(self.cancelled)?;
        if self.evaluations >= self.limit {
            return Err(BackboneContinuationError::EvaluationLimit {
                evaluations: self.evaluations,
            });
        }
        self.chart.physical(z, &mut self.physical)?;
        if self.problem.dim() != self.chart.pack.d || !self.problem.autonomous() {
            return Err(invalid(
                "source must preserve its complete autonomous state chart",
            ));
        }
        self.evaluations += 1;
        hb_residual(
            &Checked {
                source: self.problem,
                dimension: self.chart.pack.d,
            },
            &self.chart.pack,
            &self.physical,
            f64::NAN,
            aft,
            out,
        );
        poll(self.cancelled)?;
        if out.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "original harmonic balance produced nonfinite values",
            ));
        }
        Ok(())
    }
    fn prepare_tangent(&mut self, aft: &AftWorkspace) -> Result<()> {
        if let Some(callback) = self.callback {
            poll(self.cancelled)?;
            if self.tangent.is_none() {
                self.tangent = Some(Workspace::new(&self.chart.pack)?);
            }
            self.tangent_preparations += 1;
            self.tangent.as_mut().expect("admitted tangent").prepare(
                self.problem,
                callback,
                &self.chart.pack,
                self.physical[self.chart.pack.coeff_len()],
                aft,
            )?;
            poll(self.cancelled)?;
        }
        Ok(())
    }
}

enum Constraint<'a> {
    Amplitude(f64),
    Arclength {
        predictor: &'a [f64],
        direction: &'a [f64],
    },
}
impl Constraint<'_> {
    fn value(&self, chart: &Chart, z: &[f64]) -> f64 {
        match self {
            Self::Amplitude(amplitude) => {
                let re = chart.pack.re_index(1, chart.coordinate);
                let im = chart.pack.im_index(1, chart.coordinate);
                det::hypot(
                    z[re] * chart.units.columns[re],
                    z[im] * chart.units.columns[im],
                ) / (0.5 * amplitude)
                    - 1.0
            }
            Self::Arclength {
                predictor,
                direction,
            } => {
                let mut sum = AftSum::new(0.0);
                for ((&v, &reference), &t) in z.iter().zip(*predictor).zip(*direction) {
                    sum.add((v - reference) * t);
                }
                sum.total()
            }
        }
    }
    fn derivative(&self, chart: &Chart, z: &[f64], out: &mut [f64]) -> Result<()> {
        out.fill(0.0);
        match self {
            Self::Amplitude(amplitude) => {
                let re = chart.pack.re_index(1, chart.coordinate);
                let im = chart.pack.im_index(1, chart.coordinate);
                let x = z[re] * chart.units.columns[re];
                let y = z[im] * chart.units.columns[im];
                let radius = det::hypot(x, y);
                if radius == 0.0 {
                    return Err(OrbitError::TrivialCollapse.into());
                }
                out[re] = (x / radius) * (chart.units.columns[re] / (0.5 * amplitude));
                out[im] = (y / radius) * (chart.units.columns[im] / (0.5 * amplitude));
            }
            Self::Arclength { direction, .. } => out.copy_from_slice(direction),
        }
        Ok(())
    }
}

fn augmented<P, F>(
    eval: &Evaluation<'_, P, F>,
    raw: &[f64],
    z: &[f64],
    phase: &PhaseCondition,
    constraint: &Constraint<'_>,
) -> Vec<f64> {
    let mut r: Vec<_> = raw
        .iter()
        .zip(&eval.chart.units.rows)
        .map(|(v, unit)| v / unit)
        .collect();
    r.push(phase.value(&eval.physical));
    r.push(constraint.value(eval.chart, z));
    r
}

fn correct<P: OrbitProblem, F: FnMut() -> bool>(
    eval: &mut Evaluation<'_, P, F>,
    reference: &[f64],
    predictor: &[f64],
    constraint: Constraint<'_>,
    b: &BackboneContinuationBudget,
    step: f64,
) -> Result<(Vec<f64>, BackboneContinuationPoint)> {
    let chart = eval.chart;
    let nc = chart.pack.coeff_len();
    let nz = nc + 1;
    let nr = nc + 2;
    let phase = chart.phase(eval.problem, reference, b)?;
    let mut z = predictor.to_vec();
    let mut aft = AftWorkspace::new(&chart.pack)?;
    let mut raw = vec![0.0; nc];
    let mut plus = vec![0.0; nc];
    let mut minus = vec![0.0; nc];
    let mut jac = vec![0.0; nr * nz];
    let mut trace = Vec::with_capacity(b.max_newton);
    for _ in 0..b.max_newton {
        eval.eval(&z, &mut aft, &mut raw)?;
        let r = augmented(eval, &raw, &z, &phase, &constraint);
        let rn = length(&r);
        trace.push(length(&raw));
        if !rn.is_finite() {
            return Err(invalid("scaled corrector residual overflow"));
        }
        // Rank admission is mandatory even for an already converged seed. An
        // exact residual alone cannot exclude additional neutral directions.
        eval.prepare_tangent(&aft)?;
        constraint.derivative(chart, &z, &mut jac[(nc + 1) * nz..])?;
        for col in 0..nz {
            poll(eval.cancelled)?;
            jac[nc * nz + col] = phase.scaled_derivative(col, chart.units.columns[col]);
            if let Some(workspace) = eval.tangent.as_ref().filter(|_| col < nc) {
                for row in 0..nc {
                    jac[row * nz + col] = (workspace.coefficient_entry(row, col, &chart.pack)
                        / chart.units.rows[row])
                        * chart.units.columns[col];
                }
                continue;
            }
            let old = z[col];
            let h = if col == nc {
                (1e-6 * old.abs().max(1.0)).min(0.25 * old)
            } else {
                1e-6 * old.abs().max(1.0)
            };
            if !positive(h) || old + h == old || old - h == old {
                return Err(invalid("finite-difference column is not representable"));
            }
            z[col] = old + h;
            eval.eval(&z, &mut aft, &mut plus)?;
            z[col] = old - h;
            eval.eval(&z, &mut aft, &mut minus)?;
            z[col] = old;
            for row in 0..nc {
                jac[row * nz + col] = (plus[row] / chart.units.rows[row]
                    - minus[row] / chart.units.rows[row])
                    / (2.0 * h);
            }
        }
        let linear = linear::admit(&mut jac, nr, nz, b.maximum_condition_estimate)?;
        poll(eval.cancelled)?;
        if rn <= b.tolerance {
            let mut finer = aft.finer(&chart.pack)?;
            eval.eval(&z, &mut finer, &mut plus)?;
            let checked = augmented(eval, &plus, &z, &phase, &constraint);
            let defect = chart.units.difference_norm(&plus, &raw);
            let checked_norm = length(&checked);
            let raw_defect = crate::scaling::length(plus.iter().zip(&raw).map(|(a, b)| a - b));
            if checked_norm <= b.tolerance && defect <= b.tolerance {
                if length(&z[chart.pack.d..nc]) <= b.minimum_oscillation {
                    return Err(OrbitError::TrivialCollapse.into());
                }
                let (coeffs, omega) = chart.pack.unpack(&eval.physical);
                let first_harmonic_amplitude = 2.0
                    * det::hypot(
                        coeffs[chart.coordinate][1].re,
                        coeffs[chart.coordinate][1].im,
                    );
                if !first_harmonic_amplitude.is_finite() {
                    return Err(invalid("reported fundamental amplitude overflows"));
                }
                *trace.last_mut().expect("current corrector") = length(&plus);
                return Ok((
                    z,
                    BackboneContinuationPoint {
                        orbit: HbOrbit {
                            omega,
                            coeffs,
                            residual: length(&plus),
                            residual_trace: trace,
                            aft_samples: finer.samples,
                            sampling_defect: raw_defect,
                        },
                        first_harmonic_amplitude,
                        step,
                        scaled_residual: checked_norm,
                        scaled_sampling_defect: defect,
                        phase_residual: checked[nc],
                        constraint_residual: checked[nc + 1],
                        jacobian_condition_estimate: linear.condition,
                    },
                ));
            }
            if finer.refinement >= MAX_AFT_REFINEMENTS {
                return Err(OrbitError::SamplingUnresolved {
                    samples: finer.samples,
                    residual: length(&plus),
                    defect: raw_defect,
                }
                .into());
            }
            aft = finer;
            continue;
        }
        let delta = linear.solve(&r)?;
        let mut accepted = None;
        for backtrack in 0..12_i32 {
            let alpha = 2.0_f64.powi(-backtrack);
            let trial: Vec<_> = z.iter().zip(&delta).map(|(v, dv)| v - alpha * dv).collect();
            if trial[nc] <= 0.0 || trial.iter().any(|v| !v.is_finite()) {
                continue;
            }
            eval.eval(&trial, &mut aft, &mut plus)?;
            if length(&augmented(eval, &plus, &trial, &phase, &constraint)) < rn {
                accepted = Some(trial);
                break;
            }
        }
        match accepted {
            Some(trial) => z = trial,
            None => {
                return Err(OrbitError::NewtonStalled {
                    residual: *trace.last().expect("current corrector"),
                    iterations: trace.len(),
                    trace,
                }
                .into());
            }
        }
    }
    Err(OrbitError::NewtonStalled {
        residual: *trace.last().unwrap_or(&f64::INFINITY),
        iterations: trace.len(),
        trace,
    }
    .into())
}

fn tangent<P: OrbitProblem>(
    problem: &P,
    chart: &Chart,
    old: &[f64],
    now: &[f64],
    b: &BackboneContinuationBudget,
) -> Result<Vec<f64>> {
    let phase = chart.phase(problem, now, b)?;
    let mut direction: Vec<_> = now.iter().zip(old).map(|(a, b)| a - b).collect();
    let unprojected = length(&direction);
    let phase_direction: Vec<_> = chart
        .units
        .columns
        .iter()
        .enumerate()
        .map(|(col, &unit)| phase.scaled_derivative(col, unit))
        .collect();
    let drift = dot(&direction, &phase_direction);
    for (v, p) in direction.iter_mut().zip(&phase_direction) {
        *v -= drift * p;
    }
    let norm = length(&direction);
    if !norm.is_finite() || norm <= 64.0 * f64::EPSILON * unprojected {
        return Err(invalid(
            "seeds or accepted points do not resolve a branch tangent",
        ));
    }
    for value in &mut direction {
        *value /= norm;
    }
    Ok(direction)
}

fn align_seed(chart: &Chart, reference: &[f64], seed: &mut [f64]) {
    let re = chart.pack.re_index(1, chart.coordinate);
    let im = chart.pack.im_index(1, chart.coordinate);
    let shift = det::atan2(reference[im], reference[re]) - det::atan2(seed[im], seed[re]);
    for n in 1..=chart.pack.n_harm {
        let (sine, cosine) = det::sin_cos(n as f64 * shift);
        for k in 0..chart.pack.d {
            let re = chart.pack.re_index(n, k);
            let im = chart.pack.im_index(n, k);
            let (x, y) = (seed[re], seed[im]);
            seed[re] = cosine * x - sine * y;
            seed[im] = sine * x + cosine * y;
        }
    }
}

/// Continue a conservative periodic family with frequency and every original
/// Fourier coefficient free, through turns in frequency or seed amplitude.
///
/// Each seed is (positive first-harmonic peak, complete physical orbit). The
/// chosen original coordinate anchors only these two seeds. Their supplied
/// convergence metadata is ignored; both are corrected with every original
/// balance row, a full-waveform phase condition and the exact amplitude row.
/// All harmonics of the second seed are rotated together before correction.
/// Seed order selects direction. Later points use a normalized pseudo-arclength
/// plane rather than amplitude stepping; no physical parameter is introduced.
///
/// Conservation supplies a dependent balance equation. The rectangular system
/// retains it and uses column-equilibrated Householder QR, with explicit pivot
/// and estimated-conditioning admission even at residual-zero seeds. No row is
/// dropped, damping removed, normal equation formed or Jacobian jitter added.
/// Additional neutral directions, unresolved charts or dissipative sources
/// without a nonzero periodic family may refuse; this does not classify why.
///
/// The original AFT, Parseval physical scaling, phase and exact-island tangent
/// owners are reused. A declared exact tangent cannot fall back to probes.
/// Otherwise coefficient columns use centered probes; frequency always does.
/// Each accepted point passes a separately sampled retained-balance check.
/// Harmonic order stays fixed: this does not certify omitted tails, primitive
/// periods, stability, branch identity, bifurcations or continuous-time error.
///
/// Success contains exactly steps+2 points. Refusal publishes no partial branch.
/// Source callbacks must be deterministic and not hide evolving state. Work and
/// payload are bounded before source calls; cancellation is polled around AFT
/// and exact tangent passes, between columns, around QR, and between steps.
/// Each complete native AFT/tangent/QR call remains noninterruptible.
pub fn continue_backbone_branch<P: OrbitProblem>(
    problem: &P,
    coordinate: usize,
    seeds: [(f64, &HbOrbit); 2],
    budget: &BackboneContinuationBudget,
    mut cancelled: impl FnMut() -> bool,
) -> Result<BackboneBranch> {
    poll(&mut cancelled)?;
    if !problem.autonomous() {
        return Err(invalid("an autonomous problem is required"));
    }
    let callback = problem.island_state_jacobian();
    let chart = Chart::new(problem.dim(), coordinate, seeds, callback.is_some(), budget)?;
    let seed0 = chart.seed(seeds[0].1)?;
    let mut seed1 = chart.seed(seeds[1].1)?;
    let mut eval = Evaluation {
        problem,
        chart: &chart,
        cancelled: &mut cancelled,
        evaluations: 0,
        limit: budget.maximum_residual_evaluations,
        physical: vec![0.0; chart.pack.len()],
        callback,
        tangent: None,
        tangent_preparations: 0,
    };
    let (a, first) = correct(
        &mut eval,
        &seed0,
        &seed0,
        Constraint::Amplitude(seeds[0].0),
        budget,
        0.0,
    )?;
    align_seed(&chart, &a, &mut seed1);
    let (mut z, second) = correct(
        &mut eval,
        &a,
        &seed1,
        Constraint::Amplitude(seeds[1].0),
        budget,
        0.0,
    )?;
    let mut direction = tangent(problem, &chart, &a, &z, budget)?;
    let mut points = Vec::with_capacity(budget.steps + 2);
    points.push(first);
    points.push(second);
    let mut step = budget.step_sizes[0];
    for _ in 0..budget.steps {
        poll(eval.cancelled)?;
        let mut accepted = None;
        for _ in 0..budget.max_attempts_per_step {
            let predictor: Vec<_> = z
                .iter()
                .zip(&direction)
                .map(|(a, t)| a + step * t)
                .collect();
            if predictor[chart.pack.coeff_len()] > 0.0 {
                let constraint = Constraint::Arclength {
                    predictor: &predictor,
                    direction: &direction,
                };
                match correct(&mut eval, &z, &predictor, constraint, budget, step) {
                    Ok(point) => {
                        accepted = Some(point);
                        break;
                    }
                    Err(
                        BackboneContinuationError::Numerical(
                            OrbitError::NewtonStalled { .. }
                            | OrbitError::SingularJacobian
                            | OrbitError::TrivialCollapse,
                        )
                        | BackboneContinuationError::IllConditioned { .. },
                    ) => {}
                    Err(error) => return Err(error),
                }
            }
            if step <= budget.step_sizes[1] {
                break;
            }
            step = (0.5 * step).max(budget.step_sizes[1]);
        }
        let Some((next, point)) = accepted else {
            return Err(OrbitError::ContinuationExhausted {
                steps: points.len(),
            }
            .into());
        };
        let mut next_direction = tangent(problem, &chart, &z, &next, budget)?;
        if dot(&next_direction, &direction) < 0.0 {
            for v in &mut next_direction {
                *v = -*v;
            }
        }
        direction = next_direction;
        z = next;
        points.push(point);
        step = (step * 1.3).min(budget.step_sizes[2]);
    }
    poll(eval.cancelled)?;
    Ok(BackboneBranch {
        coordinate,
        points,
        residual_evaluations: eval.evaluations,
        tangent_preparations: eval.tangent_preparations,
    })
}
