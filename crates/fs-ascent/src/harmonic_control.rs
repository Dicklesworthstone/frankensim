//! Coherent source control through the original BEM energy adjoint and L-BFGS.
//!
//! Only authored actuators vary. Every other source remains at its original
//! full-bank baseline, including uncontrollable disturbance sources.

use crate::{LbfgsError, LbfgsReport, LbfgsState, StopReason, StopRule};
use fs_bem::observation_transfers::{PreparedHarmonicEnergyControl, TransferError};
use fs_exec::CancelGate;
use fs_math::{c64::C64, dd::Dd};

/// Maximum independently controlled complex source coordinates in one study.
pub const MAX_HARMONIC_ACTUATORS: usize = 256;
const MAX_MEMORY: usize = 64;
const MAX_ITERATIONS: usize = 1_000_000;
const MAX_EVALUATIONS: usize = 1_000_000;

/// One actuator in the original physical source order.
///
/// Its complex coefficient is baseline[source] + scale_si*(x_re + i*x_im).
/// The positive scale has the original source's SI unit; different source kinds
/// need distinct explicitly authored scales. It is not a hard amplitude limit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HarmonicActuator {
    pub source: usize,
    pub scale_si: f64,
    /// Nonnegative dimensionless coefficient of 0.5*|(a-reference)/scale|².
    /// This is an authored control-effort penalty, not electrical/mechanical
    /// power, actuator efficiency, saturation, or a material dissipation law.
    pub effort_weight: f64,
    /// Preferred complex source coefficient, in that source's physical SI unit.
    pub reference: C64,
}

/// Whole-study limits and explicit physical-to-numerical normalization.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicControlOptions {
    /// Positive mean-energy scale [J]. The optimizer minimizes J/E_scale plus
    /// the dimensionless authored effort penalty. Physical J is reported too.
    pub energy_scale_j: f64,
    /// Original L-BFGS curvature memory, in 1..=64.
    pub memory: usize,
    /// Cumulative accepted iterations, at most one million, across all resumes.
    pub maximum_iterations: usize,
    /// Cumulative callback attempts, including initialization and refusals,
    /// in 1..=one million. Original StopRule::Budget leaves can tighten this.
    pub maximum_evaluations: usize,
    /// Conservative logical numerical workspace, including the original energy
    /// owner, full source buffers, L-BFGS vectors/curvature and bounded history.
    /// Allocator overhead and the borrowed solved transfer bank are excluded.
    pub maximum_workspace_bytes: usize,
}

/// Original solver attribution with the original physical refusal retained.
pub type HarmonicControlError = LbfgsError<TransferError>;

/// One accepted result, with optimizer stop attribution kept intact.
/// Gradient norms and optimizer objective use the declared scaled coordinates.
#[derive(Clone, Debug)]
pub struct HarmonicControlReport {
    pub optimizer: LbfgsReport,
    /// Original cycle-mean pressure/velocity target mismatch energy [J].
    pub mean_residual_energy_j: f64,
    /// Dimensionless authored source effort, distinct from physical field energy.
    pub effort_penalty: f64,
}

#[derive(Clone, Copy)]
struct AcceptedEnergy {
    physical_j: f64,
    effort: f64,
}

/// A resumable solve tied to one immutable physical operator and source mapping.
///
/// The objective is J(a)/E_scale + sum_i lambda_i*|(a_i-reference_i)/scale_i|²/2.
/// J and its exact complex Cartesian gradient come from the original BEM owner,
/// including actual near-field velocity and every coherent cross term. The
/// selected real/imaginary gradient is scale_i*dJ/da_i/E_scale plus the exact
/// effort derivative; there is no extra conjugation or peak/RMS factor.
///
/// The existing fallible L-BFGS owns every Wolfe probe, curvature pair, budget
/// and stopping reason. Completed accepted iterations can resume bit-for-bit.
/// A refused/cancelled probe retains the accepted coefficients, objective,
/// gradient, history and curvature; its spent attempt remains charged. Neither
/// mutable optimizer state nor replacement physics can be supplied on resume.
///
/// This controls the fixed-frequency, fixed-operator linear model. It adds no
/// geometry/material derivative, hard actuator limit, transient controller,
/// stability/robustness claim, or guarantee that a desired field is realizable.
pub struct HarmonicControlStudy<'a> {
    objective: PreparedHarmonicEnergyControl<'a>,
    baseline: Vec<C64>,
    actuators: Vec<HarmonicActuator>,
    work_coefficients: Vec<C64>,
    accepted_coefficients: Vec<C64>,
    accepted_energy: AcceptedEnergy,
    state: LbfgsState,
    options: HarmonicControlOptions,
    workspace_bytes: usize,
}

impl<'a> HarmonicControlStudy<'a> {
    /// Initialize from real/imaginary decisions interleaved in actuator order.
    /// Zero decisions preserve the complete baseline. At least one distinct
    /// actuator is required; unselected source coefficients never change.
    ///
    /// Shape, finite physical inputs and the full memory bound are admitted
    /// before allocation or the first physical evaluation. The original energy
    /// owner already bounds each forward/adjoint evaluation and receiver storage.
    pub fn new(
        mut objective: PreparedHarmonicEnergyControl<'a>,
        baseline: &[C64],
        actuators: &[HarmonicActuator],
        initial_decisions: &[f64],
        options: HarmonicControlOptions,
        gate: &CancelGate,
    ) -> Result<Self, HarmonicControlError> {
        let mut cancelled = || gate.is_requested();
        poll(&mut cancelled).map_err(LbfgsError::Evaluation)?;
        let workspace_bytes = admit(&objective, baseline, actuators, initial_decisions, options)?;
        let baseline = copy(baseline)?;
        let actuators = copy(actuators)?;
        let mut work_coefficients = copy(&baseline)?;
        let mut accepted_coefficients = copy(&baseline)?;
        let mut accepted_energy = None;
        let state = LbfgsState::try_new(initial_decisions, options.memory, &mut |x| {
            let (value, gradient, energy) = evaluate(
                &mut objective,
                &baseline,
                &actuators,
                x,
                options.energy_scale_j,
                &mut work_coefficients,
                &mut cancelled,
            )?;
            accepted_energy = Some(energy);
            Ok((value, gradient))
        })?;
        accepted_coefficients.copy_from_slice(&work_coefficients);
        Ok(Self {
            objective,
            baseline,
            actuators,
            work_coefficients,
            accepted_coefficients,
            accepted_energy: accepted_energy.expect("successful initialization evaluated physics"),
            state,
            options,
            workspace_bytes,
        })
    }

    /// Immutable original optimizer checkpoint. Cloning it records the complete
    /// accepted numerical state; subsequent runs keep the same bound study.
    pub fn state(&self) -> &LbfgsState {
        &self.state
    }

    /// Accepted physical source coefficients in the ORIGINAL complete order.
    /// Unselected coordinates retain their exact authored baseline bits.
    pub fn coefficients(&self) -> &[C64] {
        &self.accepted_coefficients
    }

    pub fn baseline(&self) -> &[C64] {
        &self.baseline
    }
    pub fn actuators(&self) -> &[HarmonicActuator] {
        &self.actuators
    }
    pub fn objective(&self) -> &PreparedHarmonicEnergyControl<'a> {
        &self.objective
    }
    pub fn options(&self) -> HarmonicControlOptions {
        self.options
    }
    pub fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }
    pub fn mean_residual_energy_j(&self) -> f64 {
        self.accepted_energy.physical_j
    }
    pub fn effort_penalty(&self) -> f64 {
        self.accepted_energy.effort
    }

    /// Request additional accepted iterations under the original stop algebra.
    /// Cumulative study ceilings are never reset on resume. IterationCap, Budget
    /// and Stall remain distinct from the original gradient/objective stops.
    /// Cancellation polls before iterations and inside physical evaluations.
    /// Completed callbacks may be accepted; incomplete trials never replace the
    /// checkpoint. An interrupted search restarts without refunding its probes.
    pub fn run(
        &mut self,
        rule: &StopRule,
        additional_iterations: usize,
        gate: &CancelGate,
    ) -> Result<HarmonicControlReport, HarmonicControlError> {
        self.run_cancellable(rule, additional_iterations, || gate.is_requested())
    }

    /// The same solve with a caller-owned cancellation poll, also propagated
    /// through the original physical forward and Hermitian adjoint kernels.
    pub fn run_cancellable(
        &mut self,
        rule: &StopRule,
        additional_iterations: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HarmonicControlReport, HarmonicControlError> {
        poll(&mut cancelled).map_err(LbfgsError::Evaluation)?;
        let steps = additional_iterations.min(self.options.maximum_iterations - self.state.iters);
        // The zero-step original call performs no callback and preserves its
        // exact stop precedence, including an exhausted evaluation budget.
        let mut report = self.state.try_run(
            &mut |_| -> Result<_, TransferError> { unreachable!("zero iterations") },
            rule,
            0,
            self.options.maximum_evaluations,
        )?;
        if report.reason != StopReason::IterationCap {
            return Ok(self.report(report));
        }
        for _ in 0..steps {
            poll(&mut cancelled).map_err(LbfgsError::Evaluation)?;
            let previous = self.state.iters;
            let mut last_energy = None;
            let objective = &mut self.objective;
            let baseline = &self.baseline;
            let actuators = &self.actuators;
            let work = &mut self.work_coefficients;
            let result = self.state.try_run(
                &mut |x| {
                    let (value, gradient, energy) = evaluate(
                        objective,
                        baseline,
                        actuators,
                        x,
                        self.options.energy_scale_j,
                        work,
                        &mut cancelled,
                    )?;
                    last_energy = Some(energy);
                    Ok((value, gradient))
                },
                rule,
                1,
                self.options.maximum_evaluations,
            );
            if self.state.iters != previous {
                // One original accepted-iteration call ends on its accepted
                // probe. Keep physical diagnostics without an extra evaluation
                // or subtracting a large effort term from a tiny field energy.
                self.accepted_coefficients
                    .copy_from_slice(&self.work_coefficients);
                self.accepted_energy = last_energy.expect("an accepted step evaluated physics");
            }
            report = result?;
            if report.reason != StopReason::IterationCap {
                break;
            }
        }
        Ok(self.report(report))
    }

    fn report(&self, optimizer: LbfgsReport) -> HarmonicControlReport {
        HarmonicControlReport {
            optimizer,
            mean_residual_energy_j: self.accepted_energy.physical_j,
            effort_penalty: self.accepted_energy.effort,
        }
    }
}

fn admit(
    objective: &PreparedHarmonicEnergyControl<'_>,
    baseline: &[C64],
    actuators: &[HarmonicActuator],
    decisions: &[f64],
    options: HarmonicControlOptions,
) -> Result<usize, HarmonicControlError> {
    let k = actuators.len();
    if k == 0
        || k > MAX_HARMONIC_ACTUATORS
        || decisions.len() != 2 * k
        || baseline.len() != objective.sources().len()
        || baseline.iter().any(|&z| !finite(z))
        || decisions.iter().any(|v| !v.is_finite())
        || !options.energy_scale_j.is_finite()
        || options.energy_scale_j <= 0.0
        || !(1..=MAX_MEMORY).contains(&options.memory)
        || options.maximum_iterations > MAX_ITERATIONS
        || !(1..=MAX_EVALUATIONS).contains(&options.maximum_evaluations)
    {
        return Err(invalid("harmonic control shape, scales or solver limits"));
    }
    for (i, actuator) in actuators.iter().enumerate() {
        if actuator.source >= baseline.len()
            || actuators[..i]
                .iter()
                .any(|other| other.source == actuator.source)
            || !actuator.scale_si.is_finite()
            || actuator.scale_si <= 0.0
            || !actuator.effort_weight.is_finite()
            || actuator.effort_weight < 0.0
            || !finite(actuator.reference)
        {
            return Err(invalid(
                "harmonic actuators need distinct source indices and finite physical controls",
            ));
        }
    }
    // Curvature vectors + conservative simultaneous accepted/trial/gradient/
    // direction vectors; factor two bounds Vec's growing objective history.
    let n = 2 * k;
    let scalars = (2 * options.memory + 16)
        .checked_mul(n)
        .and_then(|v| v.checked_add(2 * (options.maximum_iterations + 1)))
        .and_then(|v| v.checked_add(2 * options.memory))
        .ok_or_else(|| invalid("harmonic solver workspace overflow"))?;
    let bytes = scalars
        .checked_mul(core::mem::size_of::<f64>())
        .and_then(|v| v.checked_add(3 * baseline.len() * core::mem::size_of::<C64>()))
        .and_then(|v| v.checked_add(k * core::mem::size_of::<HarmonicActuator>()))
        .and_then(|v| v.checked_add(objective.workspace_bytes()))
        .ok_or_else(|| invalid("harmonic solver workspace overflow"))?;
    if bytes > options.maximum_workspace_bytes {
        return Err(invalid(
            "harmonic study exceeds its complete numerical workspace allowance",
        ));
    }
    Ok(bytes)
}

fn evaluate(
    objective: &mut PreparedHarmonicEnergyControl<'_>,
    baseline: &[C64],
    actuators: &[HarmonicActuator],
    x: &[f64],
    energy_scale: f64,
    coefficients: &mut [C64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(f64, Vec<f64>, AcceptedEnergy), TransferError> {
    poll(cancelled)?;
    coefficients.copy_from_slice(baseline);
    for (i, actuator) in actuators.iter().enumerate() {
        poll(cancelled)?;
        let base = baseline[actuator.source];
        let scale = Dd::from_f64(actuator.scale_si);
        let re = Dd::from_f64(base.re) + scale * Dd::from_f64(x[2 * i]);
        let im = Dd::from_f64(base.im) + scale * Dd::from_f64(x[2 * i + 1]);
        coefficients[actuator.source] = C64::new(scalar(re)?, scalar(im)?);
    }
    let physical = objective.value_gradient(coefficients, &mut *cancelled)?;
    let energy_scale = Dd::from_f64(energy_scale);
    let mut effort = Dd::ZERO;
    let mut gradient = Vec::new();
    gradient
        .try_reserve_exact(x.len())
        .map_err(|_| bad("harmonic decision-gradient allocation refused"))?;
    for actuator in actuators {
        poll(cancelled)?;
        let scale = Dd::from_f64(actuator.scale_si);
        let value = coefficients[actuator.source];
        let source_gradient = physical.source_gradient[actuator.source];
        for (value, reference, derivative) in [
            (value.re, actuator.reference.re, source_gradient.re),
            (value.im, actuator.reference.im, source_gradient.im),
        ] {
            let mut g = (Dd::from_f64(derivative) / energy_scale) * scale;
            if actuator.effort_weight != 0.0 {
                let offset = (Dd::from_f64(value) - Dd::from_f64(reference)) / scale;
                let weight = Dd::from_f64(actuator.effort_weight);
                effort = effort + (Dd::from_f64(0.5) * weight * offset) * offset;
                g = g + weight * offset;
            }
            gradient.push(scalar(g)?);
        }
    }
    let value = scalar(Dd::from_f64(physical.mean_residual_energy_j) / energy_scale + effort)?;
    let effort = scalar(effort)?;
    if value < 0.0 || effort < 0.0 {
        return Err(bad("harmonic objective became negative"));
    }
    poll(cancelled)?;
    Ok((
        value,
        gradient,
        AcceptedEnergy {
            physical_j: physical.mean_residual_energy_j,
            effort,
        },
    ))
}

fn finite(value: C64) -> bool {
    value.re.is_finite() && value.im.is_finite()
}
fn scalar(value: Dd) -> Result<f64, TransferError> {
    if !value.hi.is_finite() || !value.lo.is_finite() || !value.to_f64().is_finite() {
        Err(bad("harmonic source control arithmetic is unrepresentable"))
    } else {
        Ok(value.to_f64())
    }
}
fn copy<T: Clone>(input: &[T]) -> Result<Vec<T>, HarmonicControlError> {
    let mut out = Vec::new();
    out.try_reserve_exact(input.len())
        .map_err(|_| invalid("harmonic study allocation refused"))?;
    out.extend_from_slice(input);
    Ok(out)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), TransferError> {
    if cancelled() {
        Err(TransferError::Cancelled)
    } else {
        Ok(())
    }
}
fn bad(what: &'static str) -> TransferError {
    TransferError::Invalid(what)
}
fn invalid(what: &'static str) -> HarmonicControlError {
    LbfgsError::InvalidInput(what)
}
