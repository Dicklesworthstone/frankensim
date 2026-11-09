//! Cold original material and radiation memory in the same exponential solve.
use super::super::{radiation, relaxation::InitialMemory};
use super::{ExponentialImpactSystem, ImpactError, ImpactSystem, invalid};
use fs_exec::CancelGate;
use fs_phs::{RelaxationBranch, linear_exact::PreparedNonlinearZoh};

fn admit_capacity(inner: &ImpactSystem, additional: usize) -> Result<(), ImpactError> {
    let states = inner
        .x
        .len()
        .checked_add(additional)
        .ok_or_else(|| invalid("complete exponential memory extent overflows"))?;
    if !PreparedNonlinearZoh::fits_reference_capacity(states, inner.modes) {
        return Err(invalid(
            "complete exponential memory exceeds original dense propagation capacity",
        ));
    }
    Ok(())
}

impl ExponentialImpactSystem {
    /// Cold-attach the original Maxwell branches to the complete current image.
    ///
    /// Calls ImpactSystem's existing source admission and fs-phs realization:
    /// H_arm=(sqrt(k)*projection.x-z)^2/2 and R_zz=1/tau. Every cross-mode
    /// projection, external port, Kelvin coordinate and initial history survives.
    /// Equilibrium stiffness must exclude these additional branch increments.
    ///
    /// Preparation propagates that complete quadratic block with the existing
    /// nonlinear contact remainder in one solve. No modal diagonalization,
    /// split relaxation step, pole merging or fitted decay is performed.
    /// The actual combined spectrum and all states must satisfy the original
    /// exponential sampled-band and dense-capacity gates. Nonlinear contact
    /// accuracy still requires refinement. Empty arms preserve the exact image.
    ///
    /// This is a before-first-step operation, not live material replacement.
    /// Nonlinear mapped-strain branches are not inferred from linear projections.
    pub fn with_relaxation_branches(
        self,
        branches: Vec<RelaxationBranch>,
        initial: InitialMemory,
        maximum_branches: usize,
        gate: &CancelGate,
    ) -> Result<Self, ImpactError> {
        if gate.is_requested() {
            return Err(ImpactError::Cancelled);
        }
        if branches.is_empty() {
            return Ok(self);
        }
        admit_capacity(&self.inner, branches.len())?;
        let inner = self
            .inner
            .with_relaxation_branches(branches, initial, maximum_branches)?;
        Self::prepare_inner(inner, gate)
    }

    /// Cold-attach the original complete passive radiation model to the same
    /// exponential equation. The existing signed mode mapping, pole/coupling
    /// admission, initial rest state and work/heat ownership are unchanged.
    ///
    /// Material and acoustic memory may be attached in either order. All
    /// original physical states remain; each load's native positional metadata
    /// still addresses its own memory. A zero-coupling model preserves the
    /// original prepared image, while nonzero loads require complete capacity
    /// and the combined sampled-band gate. This API makes no new claim that an
    /// authored passive pole model represents an acoustic geometry.
    pub fn with_radiation_load(
        mut self,
        model: &radiation::Model,
        source_modes: &[usize],
        maximum_poles: usize,
        gate: &CancelGate,
    ) -> Result<Self, ImpactError> {
        if gate.is_requested() {
            return Err(ImpactError::Cancelled);
        }
        if model
            .poles
            .iter()
            .any(|pole| pole.coupling.iter().any(|v| *v != 0.0))
        {
            let additional = model
                .poles
                .len()
                .checked_mul(2)
                .ok_or_else(|| invalid("complete radiation memory extent overflows"))?;
            admit_capacity(&self.inner, additional)?;
        }
        let before = self.inner.x.len();
        self.inner = self
            .inner
            .with_radiation_load(model, source_modes, maximum_poles)?;
        if gate.is_requested() {
            return Err(ImpactError::Cancelled);
        }
        if self.inner.x.len() == before {
            return Ok(self);
        }
        Self::prepare_inner(self.inner, gate)
    }
}
