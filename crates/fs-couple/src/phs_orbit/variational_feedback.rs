//! Exact continuous feedback tangents on the shared autonomous shooting owner.

use super::{
    AutonomousPortFeedback, InputKind, PhsOrbitError, PhsOrbitProblem, PhsStorageHessianAction,
    autonomous::feedback_problem, bad, port_flows_into, tangent::StorageTangent,
};
use core::cell::RefCell;
use fs_orbit::{
    OrbitProblem, ShootBudget, TimeDomainRhs, TimeDomainStateJacobian, VariationalShootControls,
    VariationalShootError, VariationalShootOrbit, solve_shooting_variational,
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

struct FeedbackVariationalProblem<'a> {
    original: PhsOrbitProblem<'a>,
    feedback_jacobian: &'a PhsFeedbackStateJacobian<'a>,
    derivative: RefCell<Vec<f64>>,
}

impl OrbitProblem for FeedbackVariationalProblem<'_> {
    fn dim(&self) -> usize {
        self.original.dim()
    }
    fn autonomous(&self) -> bool {
        true
    }
    fn island(&self, time: f64, state: &[f64], out: &mut [f64]) {
        self.original.island(time, state, out);
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

/// Resolve a self-excited cycle of the COMPLETE original pHS and feedback law.
/// Exact quadratic storage supplies its full Hessian automatically; a quadratic
/// part of nonlinear H is insufficient and returns TangentUnavailable. The
/// explicitly declared feedback Jacobian is required even for a constant law.
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
    let mut original = feedback_problem(system, feedback)?;
    original.tangent = match hessian {
        Some(hessian) => StorageTangent::hessian(n, hessian),
        None => StorageTangent::from_storage(system),
    };
    let problem = FeedbackVariationalProblem {
        original,
        feedback_jacobian,
        derivative: RefCell::new(vec![0.0; n * m]),
    };
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
