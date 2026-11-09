//! Original average-acceleration step with two actual thermal force endpoints.
//!
//! R0 is retained from the accepted source, while R1(q)=K1 q-b_free,1
//! belongs to the immutable candidate source. The original generalized-alpha
//! owner at rho_inf=1 solves M(a0+a1)/2+(R0(q0)+R1(q1))/2=(F0+F1)/2.
//! No stiffness averaging, new time solver, evolving geometric stiffness or
//! numerical-loss-to-heat conversion is introduced.

use fs_math::dd::Dd;
use fs_solver::norm2;
use fs_time::galpha::{ImplicitSolveConfig, ImplicitStepTelemetry, OperatorGeneralizedAlpha,
    SecondOrderOperatorWeights, SecondOrderProblem, SecondOrderState, TimeSolveError};

use super::{Action, CutElasticityRef3, CutSecondOrderProblem3, poll};
use fs_cutfem::elastic3::ElasticityError3;

/// Bounded original midpoint solve and independent endpoint force admission.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CutThermalMidpointSettings3 {
    pub solve: ImplicitSolveConfig,
    /// Total accepted-step and retained step-history cap.
    pub maximum_steps: usize,
    /// Conservative additional scalar storage for cloned state/history,
    /// Newton/FGMRES, native adaptive action scratch and returned fields.
    /// Existing immutable model caches, borrowed inputs, allocator bookkeeping
    /// and constant stack scratch are excluded. A scalar slot is eight bytes.
    pub maximum_workspace_components: usize,
    /// Absolute Euclidean gate on BOTH M a0+R0-F0 and M a1+R1-F1 [N].
    /// The original Newton gate remains separate and cannot replace this check.
    pub maximum_force_residual_n: f64,
}

/// One candidate motion and its independently recomputed endpoint force.
#[derive(Clone, Debug, PartialEq)]
pub struct CutThermalMidpointStep3 {
    pub state: SecondOrderState,
    pub telemetry: ImplicitStepTelemetry,
    /// R1 at the represented accepted endpoint, in original independent native
    /// coordinates. Strong-support rows are zero; these are not support reactions.
    pub native_resistance_n: Vec<f64>,
    pub initial_residual_norm_n: f64,
    pub endpoint_residual_norm_n: f64,
}

fn refused(phase: &'static str, reason: &'static str) -> TimeSolveError {
    TimeSolveError::ModelRefused { phase, reason: reason.into() }
}
fn native_error(phase: &'static str, error: ElasticityError3) -> TimeSolveError {
    match error {
        ElasticityError3::Cancelled => TimeSolveError::Cancelled,
        error => TimeSolveError::ModelRefused { phase, reason: error.to_string() },
    }
}

fn checked_workspace(
    n: usize,
    physical_n: usize,
    native_n: usize,
    state: &SecondOrderState,
    settings: CutThermalMidpointSettings3,
) -> Option<usize> {
    let restart = settings.solve.newton.linear_restart;
    let outer = settings.solve.max_newton_iterations;
    let cycles = settings.solve.newton.max_linear_cycles;
    // Two Arnoldi bases, Hessenberg, Newton/globalization/native residual
    // vectors, bounded inner reports, and retained outer reports.
    let mut slots = restart.checked_mul(2)?.checked_add(64)?.checked_mul(n)?
        .checked_add(restart.checked_add(1)?.checked_mul(restart)?)?
        .checked_add(restart.checked_mul(16)?)?
        .checked_add(cycles.checked_mul(restart.checked_add(1)?)?.checked_mul(4)?)?
        .checked_add(outer.checked_add(1)?.checked_mul(64)?)?
        .checked_add(physical_n.checked_mul(2)?)?
        .checked_add(native_n)?;
    // The original accepted history is borrowed; its candidate clone and the
    // current-step report are charged before any cloning or solver allocation.
    // Three descriptor copies cover the cloned history's possible Vec growth
    // during its final push, including the allocation that it replaces.
    for receipt in &state.history {
        let bytes = core::mem::size_of::<ImplicitStepTelemetry>().checked_mul(3)?
            .checked_add(receipt.newton.history.len()
                .checked_mul(core::mem::size_of::<fs_solver::NewtonIteration>())?)?;
        slots = slots.checked_add(bytes.checked_add(7)? / 8)?;
    }
    Some(slots)
}

/// Private algebraic view of the SAME candidate K1 and original mass.
///
/// Center endpoint loads before adding the resolved displacement increment.
/// Forcing passed to this view is zero; authored F0/F1 are retained in balance.
struct EndpointMidpoint<'a, 'm> {
    current: &'a CutSecondOrderProblem3<'m>,
    base: &'a [f64],
    balance: Vec<f64>,
}
impl SecondOrderProblem for EndpointMidpoint<'_, '_> {
    fn dimension(&self) -> usize { self.current.dimension() }
    fn validate_step(&self, step: f64) -> Result<(), String> {
        self.current.validate_step(step)
    }
    fn validate_state(&self, q: &[f64], v: &[f64]) -> Result<(), String> {
        // Preserve original candidate-source admission at initial state,
        // intermediate stages and endpoint. A mixed old-q/new-T refusal is
        // conservative; this adapter never bypasses a physical source domain.
        self.current.validate_state(q, v)
    }
    fn mass_apply(&self, input: &[f64], output: &mut [f64]) {
        self.current.mass_apply(input, output);
    }
    fn damping_apply(&self, input: &[f64], output: &mut [f64]) {
        self.current.damping_apply(input, output);
    }
    fn internal_force(&self, q: &[f64], output: &mut [f64]) {
        let increment: Vec<_> = q.iter().zip(self.base).map(|(q, base)| q-base).collect();
        self.current.tangent_apply(q, &increment, output);
        for (value, base) in output.iter_mut().zip(&self.balance) {
            *value = (Dd::from_f64(*base) + Dd::from_f64(*value)).to_f64();
        }
    }
    fn tangent_apply(&self, q: &[f64], direction: &[f64], output: &mut [f64]) {
        self.current.tangent_apply(q, direction, output);
    }
    fn resistance_balance(
        &self, q: &[f64], increment: &[f64], fraction: f64,
        _: &[f64], force: &[f64], output: &mut [f64],
    ) {
        self.current.tangent_apply(q, increment, output);
        for ((value, base), force) in output.iter_mut().zip(&self.balance).zip(force) {
            *value = (Dd::from_f64(*base) - Dd::from_f64(*force)
                + Dd::from_f64(fraction) * Dd::from_f64(*value)).to_f64();
        }
    }
    fn preconditioner_apply(
        &self, q: &[f64], weights: SecondOrderOperatorWeights,
        outer: usize, inner: usize, residual: &[f64], output: &mut [f64],
    ) {
        self.current.preconditioner_apply(q, weights, outer, inner, residual, output);
    }
}

impl CutSecondOrderProblem3<'_> {
    /// Advance actual inertia between two source-owned force endpoints.
    ///
    /// This adapter must be constructed with the AFTER eigenstrain/material.
    /// before_native_resistance_n is R0(q0)=K0 q0-b_free,0 on the accepted
    /// original native owner. Supply independent external nodal forces at both
    /// endpoints; do not add either free-strain load to them. All native arrays
    /// use the SAME independent map and fixed homogeneous support assignment.
    /// The caller retains the same reference geometry, physical occupancy and
    /// consistent mass across material refresh. Initial acceleration must come
    /// from the original consistent mass solve, including any authored release.
    ///
    /// The original rho_inf=1 generalized-alpha method gives average-acceleration
    /// Newmark kinematics. Its residual uses the actual endpoint mean resistance,
    /// not new-temperature resistance at midpoint displacement. Both the held
    /// affine free-strain force and changing K are included exactly once.
    /// Frozen reference G is refused: evolving its energy requires a separate
    /// model/account. This first endpoint bridge has zero physical damping.
    ///
    /// Source-domain checks remain in the original integrator at the initial,
    /// intermediate and endpoint coordinates. An old displacement outside the
    /// new source's elastic box can refuse conservatively. No clipping occurs.
    ///
    /// The supplied checkpoint is immutable; all work is performed on a bounded
    /// clone, returned only after true initial/endpoint force and cancellation
    /// gates. The caller owns joint thermal/material publication. Energy/work
    /// must use the REPRESENTED q/v endpoints and the original consistent mass.
    /// A temporal or representability defect is never reversible heat and must
    /// pass the coupled owner's independent energy allowance.
    #[allow(clippy::too_many_arguments)]
    pub fn thermal_midpoint_candidate(
        &self,
        before: &SecondOrderState,
        before_native_resistance_n: &[f64],
        before_native_external_n: &[f64],
        after_native_external_n: &[f64],
        step_s: f64,
        settings: CutThermalMidpointSettings3,
    ) -> Result<CutThermalMidpointStep3, TimeSolveError> {
        poll(self.gate).map_err(|e| native_error("thermal midpoint", e))?;
        if self.eigenstrain.is_none() || self.prestress.is_some() {
            return Err(refused("thermal midpoint", "thermal endpoint stepping requires held eigenstrain and no geometric-force extension"));
        }
        if !settings.maximum_force_residual_n.is_finite()
            || settings.maximum_force_residual_n < 0.0
            || settings.solve.max_newton_iterations == 0
            || before.steps >= settings.maximum_steps
            || before.history.len() >= settings.maximum_steps
        {
            return Err(refused("thermal midpoint controls", "invalid force allowance or exhausted native step/history budget"));
        }
        settings.solve.newton.validate().map_err(TimeSolveError::NewtonSetup)?;
        let next_t = before.t + step_s;
        if !before.t.is_finite() || !step_s.is_finite() || step_s <= 0.0
            || !next_t.is_finite() || next_t <= before.t
        {
            return Err(TimeSolveError::InvalidTimeAdvance {
                t_bits: before.t.to_bits(), h_bits: step_s.to_bits(), next_t_bits: next_t.to_bits(),
            });
        }
        let n = self.dimension();
        for field in [&before.q, &before.v, &before.a] {
            self.admit_vector(field, n).map_err(|e| native_error("thermal midpoint state", e))?;
        }
        for force in [before_native_resistance_n, before_native_external_n, after_native_external_n] {
            self.admit_vector(force, self.mass.n())
                .map_err(|e| native_error("thermal midpoint native force", e))?;
        }
        let physical_n = match self.source {
            CutElasticityRef3::Cartesian(_) => self.mass.n(),
            CutElasticityRef3::Adaptive(operator) => operator.physical_nodes().len()
                .checked_mul(3).ok_or_else(|| refused("thermal midpoint budget", "physical dimension overflow"))?,
        };
        let required = checked_workspace(n, physical_n, self.mass.n(), before, settings)
            .ok_or_else(|| refused("thermal midpoint budget", "workspace count overflow"))?;
        if required > settings.maximum_workspace_components {
            return Err(refused("thermal midpoint budget", "native midpoint workspace allowance exhausted"));
        }
        let old_resistance = self.restrict_force(before_native_resistance_n)
            .map_err(|e| native_error("thermal midpoint old resistance", e))?;
        let old_force = self.restrict_force(before_native_external_n)
            .map_err(|e| native_error("thermal midpoint old force", e))?;
        let next_force = self.restrict_force(after_native_external_n)
            .map_err(|e| native_error("thermal midpoint new force", e))?;
        let mut mass_force = vec![0.0; n];
        self.action(&before.a, &mut mass_force, Action::Mass)
            .map_err(|e| native_error("thermal midpoint initial inertia", e))?;
        let mut residual = vec![0.0; n];
        for i in 0..n {
            if i % 192 == 0 { poll(self.gate).map_err(|e| native_error("thermal midpoint initial force", e))?; }
            residual[i] = (Dd::from_f64(mass_force[i]) + Dd::from_f64(old_resistance[i])
                - Dd::from_f64(old_force[i])).to_f64();
        }
        let initial_residual_norm_n = norm2(&residual);
        if !initial_residual_norm_n.is_finite()
            || initial_residual_norm_n > settings.maximum_force_residual_n
        {
            return Err(refused("thermal midpoint initial force", "initial acceleration does not satisfy the original dynamic force allowance"));
        }
        self.action(&before.q, &mut residual, Action::Resistance)
            .map_err(|e| native_error("thermal midpoint candidate resistance", e))?;
        let mut balance = Vec::with_capacity(n);
        for i in 0..n {
            if i % 192 == 0 { poll(self.gate).map_err(|e| native_error("thermal midpoint centered balance", e))?; }
            balance.push((Dd::from_f64(0.5)
                * ((Dd::from_f64(old_resistance[i]) - Dd::from_f64(old_force[i]))
                    + (Dd::from_f64(residual[i]) - Dd::from_f64(next_force[i])))).to_f64());
        }
        self.admit_vector(&balance, n)
            .map_err(|e| native_error("thermal midpoint centered balance", e))?;
        let problem = EndpointMidpoint { current: self, base: &before.q, balance };
        let method = OperatorGeneralizedAlpha::new(n, step_s, 1.0, settings.solve);
        let mut state = before.clone();
        let zero = vec![0.0; n];
        let telemetry = method.step_controlled(
            &mut state, &problem, &zero, &mut || self.gate.is_requested(),
        )?;
        // Use the represented endpoint, independent of the centered stage
        // residual and its Krylov estimates. No equilibrium re-solve occurs.
        self.action(&state.q, &mut residual, Action::Resistance)
            .map_err(|e| native_error("thermal midpoint endpoint resistance", e))?;
        let native_resistance_n = self.expand_field(&residual)
            .map_err(|e| native_error("thermal midpoint endpoint field", e))?;
        self.action(&state.a, &mut mass_force, Action::Mass)
            .map_err(|e| native_error("thermal midpoint endpoint inertia", e))?;
        for i in 0..n {
            if i % 192 == 0 { poll(self.gate).map_err(|e| native_error("thermal midpoint endpoint force", e))?; }
            residual[i] = (Dd::from_f64(mass_force[i]) + Dd::from_f64(residual[i])
                - Dd::from_f64(next_force[i])).to_f64();
        }
        let endpoint_residual_norm_n = norm2(&residual);
        if !endpoint_residual_norm_n.is_finite()
            || endpoint_residual_norm_n > settings.maximum_force_residual_n
        {
            return Err(refused("thermal midpoint endpoint force", "represented endpoint exceeds the native dynamic force allowance"));
        }
        poll(self.gate).map_err(|e| native_error("thermal midpoint final gate", e))?;
        Ok(CutThermalMidpointStep3 { state, telemetry, native_resistance_n,
            initial_residual_norm_n, endpoint_residual_norm_n })
    }
}
