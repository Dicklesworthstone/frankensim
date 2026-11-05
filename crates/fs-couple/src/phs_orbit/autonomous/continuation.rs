//! Parameter continuation of self-excited physical-port cycles.
//!
//! The complete immutable pHS supplies its original storage, damping and
//! explicit material memory. Only the authored external feedback parameter
//! varies. fs-orbit owns continuation; the existing pHS owner evaluates every
//! equation and the original physical cycle audit supplies work and heat.

use super::{AutonomousPortFeedback, feedback_problem};
use crate::phs_orbit::{PhsOrbitAudit, PhsOrbitError, PhsOrbitProblem, audit_orbit_with_input};
use crate::phs_orbit::variational_feedback::{FeedbackVariationalProblem, feedback_tangent_problem};
use core::cell::Cell;
use fs_orbit::{ContinuableProblem, HbOrbit, IslandStateJacobian, OrbitProblem};
use fs_orbit::autonomous::{
    AutonomousBranch, AutonomousContinuationBudget, AutonomousContinuationError,
    continue_autonomous_branch,
};
use fs_phs::PortHamiltonian;

/// An instantaneous, read-only physical excitation family with one explicit
/// parameter. State/port units and ordering remain those of the original pHS.
/// Fill every effort in `out`; supplied power is `out dot flows`. Return false
/// outside the source domain. Controller, flow and material memory must be
/// explicit physical states, not mutable history hidden in this callback.
/// This is an externally supplied law, not inferred or validated source physics.
pub trait AutonomousFeedbackFamily {
    fn port_count(&self) -> usize;
    fn input(&self, parameter: f64, state: &[f64], flows: &[f64], out: &mut [f64]) -> bool;
}

/// Complete TOTAL state derivative of one physical feedback-family member.
///
/// At the supplied parameter and original state x, with original conjugate
/// flows y=G^T*grad H(x), fill the row-major m-by-n array
/// d[u(parameter,x,G^T*grad H(x))]/dx. Include explicit state dependence AND
/// the flow/Hessian chain exactly once. This is not a parameter derivative;
/// the continuation owner probes that complete original column separately.
/// Return false outside the source domain. False, unwritten or nonfinite
/// entries refuse; the solve never substitutes finite differences.
pub type PhsFeedbackFamilyStateJacobian<'a> =
    dyn Fn(f64, &[f64], &[f64], &mut [f64]) -> bool + 'a;

/// An immutable view at one physical parameter, also usable by the existing
/// autonomous HB, independent shooting/Floquet and physical-work APIs.
pub struct FixedFeedback<'a> {
    family: &'a dyn AutonomousFeedbackFamily,
    parameter: f64,
}
impl<'a> FixedFeedback<'a> {
    pub fn new(family: &'a dyn AutonomousFeedbackFamily, parameter: f64) -> Result<Self> {
        if !parameter.is_finite() { return Err(invalid("finite feedback parameter required")); }
        Ok(Self { family, parameter })
    }
    pub fn parameter(&self) -> f64 { self.parameter }
}
impl AutonomousPortFeedback for FixedFeedback<'_> {
    fn port_count(&self) -> usize { self.family.port_count() }
    fn input(&self, state: &[f64], flows: &[f64], out: &mut [f64]) -> bool {
        self.family.input(self.parameter, state, flows, out)
    }
}

/// Native continuation controls and independent physical work resolution.
#[derive(Debug, Clone)]
pub struct PhsAutonomousBranchControls {
    /// Payload allowance additionally covers this adapter's retained audits and
    /// source scratch. Physical equation scales have original state/time units.
    pub continuation: AutonomousContinuationBudget,
    /// First physical audit grid; the second uses 2*N+1, up to 65,535 samples.
    pub audit_samples: usize,
    /// Whole-branch physical quadrature sample cap, checked before model calls.
    pub maximum_audit_evaluations: usize,
    /// Each audit's |supply-loss|/mean(|supply|+|loss|), strictly in (0,1).
    pub maximum_relative_power_residual: f64,
    /// Relative change in each mean supply/loss between the two actual grids.
    /// Uses their maximum magnitude, with exact-zero handling and no power floor.
    pub maximum_relative_audit_change: f64,
}

/// Original numerical branch plus physical accounts in the same point order.
#[derive(Debug, Clone)]
pub struct PhsAutonomousBranch {
    pub continuation: AutonomousBranch,
    /// Both original physical audits at every accepted parameter, including seeds.
    pub audits: Vec<[PhsOrbitAudit; 2]>,
    pub audit_evaluations: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub enum PhsAutonomousBranchError {
    Invalid(&'static str),
    Source(PhsOrbitError),
    Continuation(AutonomousContinuationError),
    /// No corrected power or partially admitted branch is substituted.
    PowerUnresolved { point: usize, audits: [PhsOrbitAudit; 2] },
    Cancelled,
}
impl core::fmt::Display for PhsAutonomousBranchError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "physical autonomous branch refused: {self:?}")
    }
}
impl std::error::Error for PhsAutonomousBranchError {}
type Result<T> = core::result::Result<T, PhsAutonomousBranchError>;
fn invalid(what: &'static str) -> PhsAutonomousBranchError {
    PhsAutonomousBranchError::Invalid(what)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() { Err(PhsAutonomousBranchError::Cancelled) } else { Ok(()) }
}

struct VariableFeedback<'a> {
    family: &'a dyn AutonomousFeedbackFamily,
    parameter: &'a Cell<f64>,
}
impl AutonomousPortFeedback for VariableFeedback<'_> {
    fn port_count(&self) -> usize { self.family.port_count() }
    fn input(&self, state: &[f64], flows: &[f64], out: &mut [f64]) -> bool {
        self.family.input(self.parameter.get(), state, flows, out)
    }
}
struct FamilyProblem<'a, P> {
    original: P,
    parameter: &'a Cell<f64>,
}
impl<P: OrbitProblem> OrbitProblem for FamilyProblem<'_, P> {
    fn dim(&self) -> usize { self.original.dim() }
    fn autonomous(&self) -> bool { true }
    fn island(&self, time: f64, state: &[f64], out: &mut [f64]) {
        self.original.island(time, state, out);
    }
    fn island_state_jacobian(&self) -> Option<IslandStateJacobian<Self>> {
        self.original.island_state_jacobian().map(|_| Self::state_jacobian as IslandStateJacobian<Self>)
    }
    fn time_domain_rhs(&self) -> Option<fs_orbit::TimeDomainRhs<Self>> { Some(Self::island) }
}
impl<P: OrbitProblem> FamilyProblem<'_, P> {
    fn state_jacobian(&self, time: f64, state: &[f64], out: &mut [f64]) -> bool {
        self.original.island_state_jacobian()
            .is_some_and(|callback| callback(&self.original, time, state, out))
    }
}
impl<P: OrbitProblem> ContinuableProblem for FamilyProblem<'_, P> {
    fn set_parameter(&mut self, value: f64) { self.parameter.set(value); }
}

fn continue_original<P: OrbitProblem>(
    original: P, parameter: &Cell<f64>, seeds: [(f64, &HbOrbit); 2],
    controls: &AutonomousContinuationBudget, cancelled: &mut impl FnMut() -> bool,
    source_refused: fn(&P) -> bool,
) -> Result<AutonomousBranch> {
    let mut problem = FamilyProblem { original, parameter };
    let result = continue_autonomous_branch(&mut problem, seeds, controls, cancelled);
    if source_refused(&problem.original) {
        return Err(PhsAutonomousBranchError::Source(PhsOrbitError::SourceRefused));
    }
    result.map_err(|e| match e {
        AutonomousContinuationError::Cancelled => PhsAutonomousBranchError::Cancelled,
        other => PhsAutonomousBranchError::Continuation(other),
    })
}

fn relative_change(a: f64, b: f64) -> f64 {
    let scale = a.abs().max(b.abs());
    if scale == 0.0 { 0.0 } else { (a / scale - b / scale).abs() }
}

fn audited(
    system: &PortHamiltonian, family: &dyn AutonomousFeedbackFamily,
    parameter: f64, orbit: &HbOrbit, samples: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<PhsOrbitAudit> {
    poll(cancelled)?;
    let mut interrupted = false;
    let audit = audit_orbit_with_input(system, orbit, samples, |_, state, flows, out| {
        if cancelled() { interrupted = true; return false; }
        family.input(parameter, state, flows, out)
    });
    if interrupted { return Err(PhsAutonomousBranchError::Cancelled); }
    let audit = audit.map_err(PhsAutonomousBranchError::Source)?;
    poll(cancelled)?;
    Ok(audit)
}

/// Continue complete self-excited pHS cycles through folds of one authored
/// feedback parameter. Both seeds and every later point retain ALL original
/// state coordinates, nonlinear storage, resistance and explicit material memory.
/// Their complete Fourier coefficients and independently solved frequency enter
/// the existing autonomous corrector; the pHS equations are not reimplemented.
///
/// Every point is then audited at N and 2*N+1 physical phase samples. Both power
/// balances and the separate changes in mean input power and native loss must
/// pass. No numerical loss is introduced, and externally supplied work is not
/// relabeled as material dissipation. Both system and source family are borrowed
/// immutably; refusal/cancellation publishes no partial physical branch.
///
/// Two finite grids cannot certify arbitrary omitted harmonics or the continuum
/// source law. Fixed-order branch convergence is not a stability, fold-location,
/// fundamental-period or global-hysteresis certificate. Use `FixedFeedback` at
/// an accepted parameter for the existing independent autonomous shooting and
/// Floquet analysis. This path does not vary the system's storage/material cards
/// or support hybrid events or conservative families with extra neutral modes.
/// Cancellation is checked between native residuals and within physical audits;
/// one native AFT residual and dense LU remain the atomic numerical units.
pub fn continue_autonomous_feedback_branch(
    system: &PortHamiltonian,
    family: &dyn AutonomousFeedbackFamily,
    seeds: [(f64, &HbOrbit); 2],
    controls: &PhsAutonomousBranchControls,
    cancelled: impl FnMut() -> bool,
) -> Result<PhsAutonomousBranch> {
    continue_feedback_branch(system, family, None, seeds, controls, cancelled)
}

/// Continue the same complete physical feedback family with exact state
/// tangents, reusing the shared feedback/storage Jacobian and orbit owners.
///
/// The native storage must expose its complete exact quadratic form or
/// continuous Hessian action, including every nonlinear and material-memory
/// coordinate. The authored derivative obeys [`PhsFeedbackFamilyStateJacobian`].
/// The shared owner adds G*du_total/dx to (J-R)*H'', exactly once; this adapter
/// neither reconstructs that tangent nor removes a physical state. All original
/// sources, phase/arclength conditions, independent physical audits and
/// cancellation/refusal semantics remain those of the ordinary branch entry.
/// Frequency and family-parameter columns still probe original residuals.
///
/// Unavailable native storage tangents refuse before source evaluation;
/// declared storage/feedback tangent failures refuse without a probe fallback.
/// Extra storage/feedback derivative buffers count against the same payload
/// allowance. The source and derivative must be deterministic and read-only;
/// the derivative must correspond to the exact family and storage supplied.
pub fn continue_autonomous_feedback_branch_with_tangent(
    system: &PortHamiltonian,
    family: &dyn AutonomousFeedbackFamily,
    derivative: &PhsFeedbackFamilyStateJacobian<'_>,
    seeds: [(f64, &HbOrbit); 2],
    controls: &PhsAutonomousBranchControls,
    cancelled: impl FnMut() -> bool,
) -> Result<PhsAutonomousBranch> {
    continue_feedback_branch(system, family, Some(derivative), seeds, controls, cancelled)
}

fn continue_feedback_branch(
    system: &PortHamiltonian,
    family: &dyn AutonomousFeedbackFamily,
    derivative: Option<&PhsFeedbackFamilyStateJacobian<'_>>,
    seeds: [(f64, &HbOrbit); 2],
    controls: &PhsAutonomousBranchControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<PhsAutonomousBranch> {
    poll(&mut cancelled)?;
    let points = controls.continuation.steps.checked_add(2)
        .ok_or_else(|| invalid("branch extent overflow"))?;
    let samples = controls.audit_samples;
    let finer = samples.checked_mul(2).and_then(|n| n.checked_add(1))
        .ok_or_else(|| invalid("audit extent overflow"))?;
    let evaluations = samples.checked_add(finer).and_then(|n| n.checked_mul(points))
        .ok_or_else(|| invalid("whole-branch audit extent overflow"))?;
    let tolerance = |v: f64| v.is_finite() && v > 0.0 && v < 1.0;
    if system.state_dim() == 0 || family.port_count() != system.port_dim()
        || !(32..=32_767).contains(&samples)
        || evaluations > controls.maximum_audit_evaluations
        || controls.continuation.harmonics > (samples - 1) / 2
        || !tolerance(controls.maximum_relative_power_residual)
        || !tolerance(controls.maximum_relative_audit_change)
    { return Err(invalid("complete physical ports and resolved bounded work quadrature required")); }
    let mut extra = points as u128 * core::mem::size_of::<[PhsOrbitAudit; 2]>() as u128
        + 128 * (system.state_dim() as u128 + system.port_dim() as u128);
    if derivative.is_some() {
        // Bound the native chart before requesting potentially owned storage
        // metadata or rectangular feedback derivative buffers. The orbit owner
        // retains all its finer admission checks and its own tangent payload.
        let n = system.state_dim();
        let unknowns = controls.continuation.harmonics.checked_mul(2)
            .and_then(|v| v.checked_add(1)).and_then(|v| v.checked_mul(n))
            .and_then(|v| v.checked_add(2));
        if unknowns.is_none_or(|v| v > controls.continuation.maximum_unknowns.min(2048)) {
            return Err(invalid("complete analytic branch state exceeds its allowance"));
        }
        let n = n as u128;
        let m = system.port_dim() as u128;
        // Owned exact storage metadata plus its cached vector-field tangent,
        // or the smaller native Hessian-action scratch, and complete m-by-n du.
        extra += 16 * n * n + 8 * m * n;
    }
    if extra > controls.continuation.maximum_payload_bytes as u128 {
        return Err(invalid("physical branch accounts exceed the payload allowance"));
    }
    let mut numerical = controls.continuation.clone();
    numerical.maximum_payload_bytes -= extra as usize;
    let parameter = Cell::new(seeds[0].0);
    let law = VariableFeedback { family, parameter: &parameter };
    let continuation = if let Some(derivative) = derivative {
        let at_parameter = |state: &[f64], flows: &[f64], out: &mut [f64]| {
            derivative(parameter.get(), state, flows, out)
        };
        let original = feedback_tangent_problem(system, &law, &at_parameter, None)
            .map_err(PhsAutonomousBranchError::Source)?;
        if original.island_state_jacobian().is_none() {
            return Err(invalid("complete native storage tangent required for analytic continuation"));
        }
        continue_original(original, &parameter, seeds, &numerical, &mut cancelled,
            FeedbackVariationalProblem::source_refused)?
    } else {
        let original = feedback_problem(system, &law).map_err(PhsAutonomousBranchError::Source)?;
        continue_original(original, &parameter, seeds, &numerical, &mut cancelled,
            PhsOrbitProblem::source_refused)?
    };
    let mut audits = Vec::with_capacity(continuation.points.len());
    for (index, point) in continuation.points.iter().enumerate() {
        let pair = [
            audited(system, family, point.parameter, &point.orbit, samples, &mut cancelled)?,
            audited(system, family, point.parameter, &point.orbit, finer, &mut cancelled)?,
        ];
        let changed = relative_change(pair[0].mean_supplied_power_w, pair[1].mean_supplied_power_w)
            .max(relative_change(pair[0].mean_resistance_loss_w, pair[1].mean_resistance_loss_w));
        if pair.iter().any(|a| a.relative_power_residual > controls.maximum_relative_power_residual)
            || changed > controls.maximum_relative_audit_change
        { return Err(PhsAutonomousBranchError::PowerUnresolved { point: index, audits: pair }); }
        audits.push(pair);
    }
    poll(&mut cancelled)?;
    Ok(PhsAutonomousBranch { continuation, audits, audit_evaluations: evaluations })
}

#[cfg(test)]
mod tests;
