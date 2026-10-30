//! Exact resisting-flow tangents of the same nonlinear law used in time stepping.
//! The positive resisting flow D enters xdot=(J-R)e-D(x,e)+Gu. Its total
//! state derivative T=D_x+D_e W is retained, without inverting W, dropping
//! cross terms, symmetrizing, or replacing the differential slope by D/e.
use super::*;
use crate::harmonic::C64;

/// Read-only callbacks with the SAME independent state/effort directions as
/// StepWorkspace::step_into_dissipative_controlled. Every output entry must be
/// filled; false, nonfinite or incomplete results refuse the complete image.
/// At a corner the callback must refuse unless a two-sided derivative exists.
/// Hidden evolving constitutive histories must be explicit source coordinates.
#[derive(Clone, Copy)]
pub struct NonlinearDissipation<'a> {
    pub flow: &'a dyn Fn(&[f64], &[f64], &mut [f64]) -> bool,
    pub tangent: &'a dyn Fn(&[f64], &[f64], &[f64], &[f64], &mut [f64]) -> bool,
    /// Exactly n total state-direction actions. No truncation or finite differences.
    /// Callback internals are caller-owned, not a measured FLOP/wall-time budget.
    pub max_tangent_actions: usize,
}

impl EquilibriumLinearization {
    /// Include the original nonlinear resisting flow in BOTH the stationary
    /// balance and its complete tangent. Passive D may have a negative local
    /// slope at a driven bias: incremental instability is retained, not clamped.
    /// Source J/R/G remain constant; D has no independent input dependence.
    pub fn from_dissipative_ode(
        source: &PortHamiltonian, point: EquilibriumPoint<'_>, budget: LinearizationBudget,
        hessian: &dyn Fn(&[f64], &[f64], &mut [f64]) -> bool,
        dissipation: NonlinearDissipation<'_>, mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, LinearizationError> {
        Self::prepare(source.n, source.n, source.m, &source.j, &source.r, &source.g,
            source.storage.as_ref(), point, budget, hessian, Some(dissipation), &mut cancelled)
    }

    /// Same source law with every original algebraic effort/constraint retained.
    /// D and its action cover the FULL state, including algebraic equations.
    /// This is a continuous local model, not a new descriptor time integrator.
    pub fn from_dissipative_descriptor(
        source: &DescriptorPortHamiltonian, point: EquilibriumPoint<'_>, budget: LinearizationBudget,
        hessian: &dyn Fn(&[f64], &[f64], &mut [f64]) -> bool,
        dissipation: NonlinearDissipation<'_>, mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, LinearizationError> {
        Self::prepare(source.n, source.n_diff, source.m, &source.j, &source.r, &source.g,
            source.storage.as_ref(), point, budget, hessian, Some(dissipation), &mut cancelled)
    }

    /// Original D(x0,e0), not its tangent times the biased state. None denotes
    /// the ordinary storage-only path, whose arithmetic is unchanged.
    pub fn equilibrium_nonlinear_resisting_flow(&self) -> Option<&[f64]> {
        self.nonlinear_flow.as_deref()
    }

    /// Actual e0 dot D(x0,e0) [W], separate from signed INCREMENTAL mode/AC work.
    /// Raw roundoff-sized negative values are disclosed, never clipped to heat.
    pub fn equilibrium_nonlinear_dissipated_power_w(&self) -> f64 {
        self.nonlinear_power_w
    }
}

impl NonlinearDissipation<'_> {
    pub(super) fn evaluate(
        &self, state: &[f64], effort: &[f64], cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Vec<f64>, f64), LinearizationError> {
        poll(cancelled)?;
        let mut flow = vec![f64::NAN; state.len()];
        if !(self.flow)(state, effort, &mut flow) || flow.iter().any(|x| !x.is_finite()) {
            return Err(LinearizationError::DissipationRefused { column: None });
        }
        poll(cancelled)?;
        let (mut power, mut correction, mut absolute) = (0.0, 0.0, 0.0);
        for (&e, &d) in effort.iter().zip(&flow) {
            let term = e * d;
            add(&mut power, &mut correction, term);
            absolute += term.abs();
        }
        // Same original-law power admission scale as the prepared time owner.
        // This pointwise check is NOT a global or incremental passivity proof.
        let roundoff = 8.0 * f64::EPSILON * state.len() as f64 * absolute;
        if !power.is_finite() || !absolute.is_finite() || !roundoff.is_finite() {
            return Err(invalid("nonlinear equilibrium work is not representable"));
        }
        if power < -roundoff {
            return Err(LinearizationError::ActiveDissipation { power_w: power, roundoff_bound_w: roundoff });
        }
        Ok((flow, power))
    }

    pub(super) fn total_tangent(
        &self, state: &[f64], effort: &[f64], hessian: &[f64], nd: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<f64>, LinearizationError> {
        let n = state.len();
        let mut result = vec![0.0; n*n];
        let mut dx = vec![0.0; n];
        let mut de = vec![0.0; n];
        let mut column = vec![f64::NAN; n];
        for j in 0..n {
            poll(cancelled)?;
            dx[j] = 1.0;
            for i in 0..n {
                de[i] = if i < nd { hessian[i*n+j] } else if i == j { 1.0 } else { 0.0 };
            }
            column.fill(f64::NAN);
            if !(self.tangent)(state, effort, &dx, &de, &mut column)
                || column.iter().any(|x| !x.is_finite()) {
                return Err(LinearizationError::DissipationRefused { column: Some(j) });
            }
            poll(cancelled)?;
            dx[j] = 0.0;
            for i in 0..n { result[i*n+j] = column[i]; }
        }
        Ok(result)
    }
}

impl HarmonicModel {
    /// Total resisting-flow derivative T=D_x+D_e W in original state order.
    /// Not a resistance matrix: T can be nonsymmetric and need not be PSD.
    pub fn nonlinear_resisting_tangent(&self) -> Option<&[f64]> {
        self.resisting_tangent.as_deref()
    }

    // One physical contraction shared by AC, mode and parameter observations.
    // Callers poll cancellation between rows; no allocation or new factorization.
    pub(in crate::harmonic) fn nonlinear_loss_row(&self, row: usize, state: &[C64]) -> C64 {
        let Some(t) = &self.resisting_tangent else { return C64::ZERO; };
        let mut sum = C64::ZERO;
        let mut correction = C64::ZERO;
        for (j, &x) in state.iter().enumerate() {
            let term = x.scale(t[row*self.n+j]) - correction;
            let next = sum + term;
            correction = (next-sum)-term;
            sum = next;
        }
        sum
    }
}
