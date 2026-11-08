//! Complete continuous feedback tangents on the shared periodic orbit owners.

use super::{
    AutonomousPortFeedback, InputKind, PeriodicPortInput, PhsOrbitError, PhsOrbitProblem,
    PhsStorageHessianAction, autonomous::feedback_problem, bad, port_flows_into,
    tangent::StorageTangent,
};
use core::cell::RefCell;
use fs_orbit::{
    AdaptiveHbOrbit, HbBudget, HbRefinementBudget, HbScales, IslandStateJacobian, OrbitProblem,
    ScaledHbOrbit, ShootBudget, TimeDomainRhs, TimeDomainStateJacobian, VariationalShootControls,
    VariationalShootError, VariationalShootOrbit, solve_hb_autonomous_adaptive_scaled_seeded,
    solve_hb_autonomous_scaled_seeded, solve_shooting_variational,
};
use fs_phs::PortHamiltonian;

/// Complete TOTAL derivative of an original instantaneous feedback effort law.
///
/// Given original state x and its original conjugate flows y=G^T*grad H(x),
/// write every entry of d[u(x,G^T*grad H(x))]/dx into a row-major m-by-n array,
/// with physical port rows and original state columns. This MUST include both
/// explicit state dependence and the flow/Hessian chain, exactly once. A partial
/// derivative at fixed flow is insufficient. The adapter adds G*du/dx to the
/// original storage tangent; it never adds another flow chain or probes the law.
///
/// False, unwritten or nonfinite entries refuse the entire solve. Callbacks must
/// be deterministic and read-only. Controller, transport and constitutive memory
/// belong in explicit original states; there is no externally prescribed clock.
pub type PhsFeedbackStateJacobian<'a> = dyn Fn(&[f64], &[f64], &mut [f64]) -> bool + 'a;

/// Complete parameter derivative of the SAME original feedback at FIXED state.
///
/// Arguments are original state `x`, conjugate flows `y=G^T grad H`, and their
/// direct parameter derivative `y_p|x=G^T partial_p grad H+G_p^T grad H`.
/// Write every original port's TOTAL `partial_p u(x,y(x,p),p)|x`, including
/// explicit parameter dependence and the flow chain through the supplied third
/// argument exactly once. The state chain through `x_p` belongs to the separate
/// [`PhsFeedbackStateJacobian`] and MUST NOT be included here. External periodic
/// drive derivatives are supplied separately by their authored coefficients.
///
/// The callback must describe the same immutable, deterministic, read-only law
/// and original physical chart. Return false for unsupported source states;
/// incomplete/nonfinite writes refuse. An explicit callback writing zeros can
/// declare a fixed feedback law only when its COMPLETE fixed-state derivative
/// is zero. No missing constitutive/port/feedback chain rule is inferred.
pub type PhsFeedbackParameterAction<'a> =
    dyn Fn(&[f64], &[f64], &[f64], &mut [f64]) -> bool + 'a;

pub(super) struct FeedbackVariationalProblem<'a> {
    pub(super) original: PhsOrbitProblem<'a>,
    pub(super) external: Option<&'a PeriodicPortInput>,
    feedback_jacobian: &'a PhsFeedbackStateJacobian<'a>,
    derivative: RefCell<Vec<f64>>,
}

impl OrbitProblem for FeedbackVariationalProblem<'_> {
    fn dim(&self) -> usize {
        self.original.dim()
    }
    fn autonomous(&self) -> bool {
        self.external.is_none()
    }
    fn island(&self, time: f64, state: &[f64], out: &mut [f64]) {
        self.original
            .island_with_external_input(time, state, out, self.external);
    }
    fn island_state_jacobian(&self) -> Option<IslandStateJacobian<Self>> {
        self.original
            .tangent
            .available()
            .then_some(Self::state_jacobian)
    }
    fn time_domain_rhs(&self) -> Option<TimeDomainRhs<Self>> {
        Some(Self::island)
    }
    fn time_domain_state_jacobian(&self) -> Option<TimeDomainStateJacobian<Self>> {
        self.original
            .tangent
            .available()
            .then_some(Self::state_jacobian)
    }
}

impl FeedbackVariationalProblem<'_> {
    pub(super) fn source_refused(&self) -> bool {
        self.original.source_refused()
    }

    /// The same original source for work-conjugate physical port observation.
    pub(super) fn original_problem(&self) -> &PhsOrbitProblem<'_> {
        &self.original
    }

    fn state_jacobian(&self, time: f64, state: &[f64], out: &mut [f64]) -> bool {
        let n = self.original.system.state_dim();
        let m = self.original.system.port_dim();
        if state.len() != n
            || out.len() != n * n
            || state.iter().any(|v| !v.is_finite())
            || !self.original.state_jacobian(time, state, out)
        {
            return false;
        }
        // Re-evaluate the actual effort at this same state. The derivative
        // callback never depends on a previous RHS stage's scratch contents.
        let mut effort = self.original.effort.borrow_mut();
        effort.fill(f64::NAN);
        self.original.system.storage().gradient(state, &mut effort);
        if effort.iter().any(|v| !v.is_finite()) {
            return false;
        }
        let InputKind::Feedback { flows, .. } = &self.original.input else {
            return false;
        };
        let mut flows = flows.borrow_mut();
        if !port_flows_into(self.original.system, &effort, &mut flows) {
            return false;
        }
        let mut derivative = self.derivative.borrow_mut();
        derivative.fill(f64::NAN);
        if !(self.feedback_jacobian)(state, &flows, &mut derivative)
            || derivative.iter().any(|v| !v.is_finite())
        {
            return false;
        }
        let (_, _, g) = self.original.system.structure();
        for i in 0..n {
            for j in 0..n {
                let at = i * n + j;
                let mut correction = 0.0;
                for p in 0..m {
                    let adjusted = g[i * m + p] * derivative[p * n + j] - correction;
                    let next = out[at] + adjusted;
                    correction = (next - out[at]) - adjusted;
                    out[at] = next;
                }
            }
        }
        out.iter().all(|v| v.is_finite())
    }
}

pub(super) fn feedback_tangent_problem<'a>(
    system: &'a PortHamiltonian,
    feedback: &'a dyn AutonomousPortFeedback,
    feedback_jacobian: &'a PhsFeedbackStateJacobian<'a>,
    hessian: Option<&'a PhsStorageHessianAction<'a>>,
) -> Result<FeedbackVariationalProblem<'a>, PhsOrbitError> {
    let n = system.state_dim();
    let derivative_count = n
        .checked_mul(system.port_dim())
        .ok_or_else(|| bad("complete feedback derivative extent overflows"))?;
    let mut original = feedback_problem(system, feedback)?;
    original.tangent = match hessian {
        Some(hessian) => StorageTangent::hessian(n, hessian),
        None => StorageTangent::from_storage(system),
    };
    Ok(FeedbackVariationalProblem {
        original,
        external: None,
        feedback_jacobian,
        derivative: RefCell::new(vec![0.0; derivative_count]),
    })
}

/// Resolve a self-excited harmonic waveform from a complete physical seed.
///
/// The original pHS flow `(J-R)*grad H + G*u(x,G^T*grad H)` supplies every
/// residual. Its complete analytic derivative `(J-R)*H'' + G*du_total/dx`
/// supplies the shared AFT Newton matrix, once per physical sample rather than
/// by perturbing every Fourier coefficient. Exact quadratic storage metadata or
/// the storage's own continuous tangent supplies H''; a merely quadratic part
/// of nonlinear storage is insufficient. The feedback derivative is TOTAL as
/// specified by [`PhsFeedbackStateJacobian`].
///
/// The packed seed uses [`fs_orbit::HbOrbit::packed_seed`] with an autonomous
/// anchor: all original DC and complex harmonics, then omega in rad/s. No
/// coordinate is removed or permuted. A full-waveform Parseval phase condition
/// uses this seed and the physical state scales, so a static first coordinate
/// or a first coordinate with no fundamental does not suppress the cycle.
/// All state/equation scales, coefficients, material-memory coordinates and
/// source work retain their original physical meaning. Passive R loss remains
/// separate from active supply; no negative material resistance is inserted.
///
/// This is a fixed-order independently sampled balance, not a stability,
/// primitive-period, omitted-tail or excitation-physics certificate. Sources
/// must be deterministic/read-only. No partial orbit escapes a refusal.
///
/// # Errors
/// Invalid original ports, unavailable complete storage tangent, malformed
/// physical seed/scales, source or declared tangent refusal, or native HB
/// numerical failure. Declared tangents never fall back to finite differences;
/// an already converged independently checked seed need not evaluate them.
pub fn solve_autonomous_feedback_harmonic_balance_scaled_seeded(
    system: &PortHamiltonian,
    feedback: &dyn AutonomousPortFeedback,
    feedback_jacobian: &PhsFeedbackStateJacobian<'_>,
    seed: Vec<f64>,
    budget: &HbBudget,
    scales: &HbScales,
) -> Result<ScaledHbOrbit, PhsOrbitError> {
    solve_harmonic(
        system,
        feedback,
        feedback_jacobian,
        seed,
        budget,
        scales,
        None,
    )
}

/// The same full-state scaled feedback solve with an authored complete original
/// storage Hessian action. Its contract is [`PhsStorageHessianAction`]; all
/// original memory cross derivatives and nonlinear stiffness remain present.
/// The TOTAL feedback derivative already includes this same flow/Hessian chain.
pub fn solve_autonomous_feedback_harmonic_balance_analytic_scaled_seeded(
    system: &PortHamiltonian,
    feedback: &dyn AutonomousPortFeedback,
    feedback_jacobian: &PhsFeedbackStateJacobian<'_>,
    seed: Vec<f64>,
    budget: &HbBudget,
    scales: &HbScales,
    hessian: &PhsStorageHessianAction<'_>,
) -> Result<ScaledHbOrbit, PhsOrbitError> {
    solve_harmonic(
        system,
        feedback,
        feedback_jacobian,
        seed,
        budget,
        scales,
        Some(hessian),
    )
}

/// Refine the same complete feedback waveform with independently checked
/// retained and omitted equation harmonics. The original seed's full-state
/// phase reference survives every order extension; no harmonic is discarded
/// merely because the first coordinate is static. Acceptance requires the
/// shared solver's waveform, frequency, sampling and finite-band residual gates.
/// This measured refinement does not bound an infinite tail or certify stability.
pub fn solve_autonomous_feedback_harmonic_balance_adaptive_scaled_seeded(
    system: &PortHamiltonian,
    feedback: &dyn AutonomousPortFeedback,
    feedback_jacobian: &PhsFeedbackStateJacobian<'_>,
    seed: Vec<f64>,
    budget: &HbBudget,
    refinement: &HbRefinementBudget,
    scales: &HbScales,
) -> Result<AdaptiveHbOrbit, PhsOrbitError> {
    solve_harmonic_adaptive(
        system,
        feedback,
        feedback_jacobian,
        seed,
        budget,
        refinement,
        scales,
        None,
    )
}

/// Adaptive full-state feedback HB using an explicitly authored original
/// Hessian. Both exact derivatives and the unchanged original vector field are
/// evaluated at every admitted order; all refusals and physical units retain
/// the fixed-order analytic entry's semantics.
#[allow(clippy::too_many_arguments)]
pub fn solve_autonomous_feedback_harmonic_balance_analytic_adaptive_scaled_seeded(
    system: &PortHamiltonian,
    feedback: &dyn AutonomousPortFeedback,
    feedback_jacobian: &PhsFeedbackStateJacobian<'_>,
    seed: Vec<f64>,
    budget: &HbBudget,
    refinement: &HbRefinementBudget,
    scales: &HbScales,
    hessian: &PhsStorageHessianAction<'_>,
) -> Result<AdaptiveHbOrbit, PhsOrbitError> {
    solve_harmonic_adaptive(
        system,
        feedback,
        feedback_jacobian,
        seed,
        budget,
        refinement,
        scales,
        Some(hessian),
    )
}

fn solve_harmonic(
    system: &PortHamiltonian,
    feedback: &dyn AutonomousPortFeedback,
    feedback_jacobian: &PhsFeedbackStateJacobian<'_>,
    seed: Vec<f64>,
    budget: &HbBudget,
    scales: &HbScales,
    hessian: Option<&PhsStorageHessianAction<'_>>,
) -> Result<ScaledHbOrbit, PhsOrbitError> {
    let problem = feedback_tangent_problem(system, feedback, feedback_jacobian, hessian)?;
    if !problem.original.tangent.available() {
        return Err(bad(
            "complete original storage tangent required for analytic feedback HB",
        ));
    }
    let result = solve_hb_autonomous_scaled_seeded(&problem, seed, budget, scales);
    if problem.original.source_refused() {
        return Err(PhsOrbitError::SourceRefused);
    }
    result.map_err(PhsOrbitError::Orbit)
}

#[allow(clippy::too_many_arguments)]
fn solve_harmonic_adaptive(
    system: &PortHamiltonian,
    feedback: &dyn AutonomousPortFeedback,
    feedback_jacobian: &PhsFeedbackStateJacobian<'_>,
    seed: Vec<f64>,
    budget: &HbBudget,
    refinement: &HbRefinementBudget,
    scales: &HbScales,
    hessian: Option<&PhsStorageHessianAction<'_>>,
) -> Result<AdaptiveHbOrbit, PhsOrbitError> {
    let problem = feedback_tangent_problem(system, feedback, feedback_jacobian, hessian)?;
    if !problem.original.tangent.available() {
        return Err(bad(
            "complete original storage tangent required for analytic feedback HB",
        ));
    }
    let result =
        solve_hb_autonomous_adaptive_scaled_seeded(&problem, seed, budget, refinement, scales);
    if problem.original.source_refused() {
        return Err(PhsOrbitError::SourceRefused);
    }
    result.map_err(PhsOrbitError::Orbit)
}

/// Resolve a self-excited cycle of the COMPLETE original pHS and feedback law.
/// Exact quadratic metadata or [`fs_phs::Storage::continuous_tangent`] supplies
/// the original full Hessian automatically; a quadratic part of nonlinear H is
/// insufficient and returns TangentUnavailable. The explicitly declared feedback
/// Jacobian is required even for a constant law.
///
/// Every RHS remains `(J-R)*grad H + G*u(x,G^T*grad H)`; the only added callback
/// supplies its full continuous Jacobian `(J-R)*H'' + G*du_total/dx`. The shared
/// fs-orbit owner differentiates the SAME RK4 stages, including the unknown-period
/// column. Its full-state phase plane does not require the first coordinate to
/// move. Original material memory and every Floquet multiplier remain present.
/// No finite differences, phase-coordinate permutation or new integration is used.
///
/// Seed/scales are physical state coordinates, and period_guess is in seconds.
/// The result concerns one smooth discrete map: it does not certify temporal
/// convergence, primitive period, global attraction or excitation-law validity.
/// Active feedback work remains physical port supply, separate from original
/// passive R loss; no negative material resistance or heat account is invented.
/// All sources must be deterministic and read-only. A refusal/cancellation
/// publishes no partial cycle and leaves the original sources unchanged.
///
/// Before allocation, the adapter reserves 16*n²+8*n*m+64*n+48*m+256 bytes and
/// 32*n²+(4*maps*steps+1)*(4*n³+32*n²+16*n²*m+32*n*m+32*m) work units from
/// the supplied controls. These include original storage metadata, gradient/
/// flow scratch and full feedback derivatives. Callback internal resources are
/// caller-owned; at most n Hessian actions and one feedback Jacobian are used
/// per native Jacobian evaluation. The remaining caps admit the shared solver.
#[allow(clippy::too_many_arguments)]
pub fn solve_autonomous_feedback_shooting_variational(
    system: &PortHamiltonian,
    feedback: &dyn AutonomousPortFeedback,
    feedback_jacobian: &PhsFeedbackStateJacobian<'_>,
    seed: &[f64],
    period_guess: f64,
    budget: &ShootBudget,
    controls: &VariationalShootControls,
    cancelled: impl FnMut() -> bool,
) -> Result<VariationalShootOrbit, PhsOrbitError> {
    solve(
        system,
        feedback,
        feedback_jacobian,
        seed,
        period_guess,
        budget,
        controls,
        None,
        cancelled,
    )
}

/// The same complete self-excited problem with an authored original nonlinear
/// storage Hessian action. Its contract is [`PhsStorageHessianAction`], including
/// every material-memory coordinate. The feedback Jacobian is still TOTAL: its
/// flow dependence already includes that same storage Hessian. Neither callback
/// may omit entries or silently fall back to numerical differentiation.
#[allow(clippy::too_many_arguments)]
pub fn solve_autonomous_feedback_shooting_variational_analytic(
    system: &PortHamiltonian,
    feedback: &dyn AutonomousPortFeedback,
    feedback_jacobian: &PhsFeedbackStateJacobian<'_>,
    seed: &[f64],
    period_guess: f64,
    budget: &ShootBudget,
    controls: &VariationalShootControls,
    hessian: &PhsStorageHessianAction<'_>,
    cancelled: impl FnMut() -> bool,
) -> Result<VariationalShootOrbit, PhsOrbitError> {
    solve(
        system,
        feedback,
        feedback_jacobian,
        seed,
        period_guess,
        budget,
        controls,
        Some(hessian),
        cancelled,
    )
}

#[allow(clippy::too_many_arguments)]
fn solve(
    system: &PortHamiltonian,
    feedback: &dyn AutonomousPortFeedback,
    feedback_jacobian: &PhsFeedbackStateJacobian<'_>,
    seed: &[f64],
    period_guess: f64,
    budget: &ShootBudget,
    controls: &VariationalShootControls,
    hessian: Option<&PhsStorageHessianAction<'_>>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<VariationalShootOrbit, PhsOrbitError> {
    if cancelled() {
        return Err(PhsOrbitError::Variational(VariationalShootError::Cancelled));
    }
    let (n, m) = (system.state_dim(), system.port_dim());
    if n == 0
        || n > 256
        || n > controls.maximum_dimension
        || feedback.port_count() != m
        || seed.len() != n
        || seed.iter().any(|v| !v.is_finite())
        || controls.state_scales.len() != n
        || controls
            .state_scales
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
        || !period_guess.is_finite()
        || period_guess <= 0.0
        || !controls.period_scale.is_finite()
        || controls.period_scale <= 0.0
        || !budget.tolerance.is_finite()
        || budget.tolerance <= 0.0
        || budget.tolerance >= 1.0
        || !controls.minimum_autonomous_excursion.is_finite()
        || controls.minimum_autonomous_excursion <= budget.tolerance
        || budget.steps_per_period == 0
        || !(1..=4096).contains(&budget.max_newton)
        || controls.maximum_period_maps == 0
    {
        return Err(bad(
            "complete finite original feedback/state data, positive physical scales and bounded autonomous shooting controls required",
        ));
    }
    let (nn, mm) = (n as u128, m as u128);
    let bytes = 16 * nn * nn + 8 * nn * mm + 64 * nn + 48 * mm + 256;
    let work = (controls.maximum_period_maps as u128)
        .checked_mul(budget.steps_per_period as u128)
        .and_then(|v| v.checked_mul(4))
        .and_then(|v| v.checked_add(1))
        .and_then(|v| {
            v.checked_mul(
                4 * nn * nn * nn + 32 * nn * nn + 16 * nn * nn * mm + 32 * nn * mm + 32 * mm,
            )
        })
        .and_then(|v| v.checked_add(32 * nn * nn))
        .ok_or_else(|| bad("autonomous pHS variational work envelope overflows"))?;
    for (quantity, requested, cap) in [
        (
            "autonomous pHS variational payload bytes",
            bytes,
            controls.maximum_payload_bytes as u128,
        ),
        (
            "autonomous pHS variational work",
            work,
            controls.maximum_work as u128,
        ),
    ] {
        if requested > cap {
            return Err(PhsOrbitError::Variational(VariationalShootError::Budget {
                quantity,
                requested,
                cap,
            }));
        }
    }
    let mut admitted = controls.clone();
    admitted.maximum_payload_bytes -= bytes as usize;
    admitted.maximum_work -= work as u64;
    let problem = feedback_tangent_problem(system, feedback, feedback_jacobian, hessian)?;
    let result = solve_shooting_variational(
        &problem,
        seed,
        period_guess,
        budget,
        &admitted,
        &mut cancelled,
    );
    match result {
        Err(VariationalShootError::Cancelled) => {
            Err(PhsOrbitError::Variational(VariationalShootError::Cancelled))
        }
        _ if problem.original.source_refused() => Err(PhsOrbitError::SourceRefused),
        other => other.map_err(PhsOrbitError::Variational),
    }
}
