//! Chemical endpoint with explicit heat and mechanical boundary work.
//!
//! The original mechanism solves xi_r = dt V_end rate_r(T_end, p_end, n_end).
//! Its material reconstruction uses U_end = U_start + Q - W at V_end, including
//! original species formation energies. The mechanical owner supplies W and
//! must independently match it to its actual force/volume trajectory.
//! Reference: <https://cantera.org/stable/reference/reactors/ideal-gas-reactor.html>.

use super::*;

/// One closed-parcel boundary exchange; positive work leaves the gas.
/// These are proposed extensive data, not an independently solved trajectory.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasReactionBoundary {
    pub final_volume_m3: f64,
    pub heat_j: f64,
    pub work_by_gas_j: f64,
}
impl GasReactionBoundary {
    pub(super) fn validate(self) -> Result<(), GasReactionError> {
        let exchange = self.heat_j - self.work_by_gas_j;
        if !positive(self.final_volume_m3)
            || !self.heat_j.is_finite()
            || !self.work_by_gas_j.is_finite()
            || !exchange.is_finite()
            || (self.heat_j != 0.0 && exchange == -self.work_by_gas_j)
            || (self.work_by_gas_j != 0.0 && exchange == self.heat_j)
        {
            return Err(invalid(
                "positive volume and representable heat/work exchange required",
            ));
        }
        Ok(())
    }
}

/// Original chemical equations and their complete volume/heat/work account.
#[derive(Clone, Debug, PartialEq)]
pub struct BoundaryGasReactionStep {
    pub reaction: GasReactionStep,
    pub initial_volume_m3: f64,
    pub boundary: GasReactionBoundary,
    /// Actual delta U + W - Q [J].
    pub energy_residual_j: f64,
    /// Signed Q/T_end [J/K]. Mechanical work transports no entropy here.
    pub endpoint_heat_entropy_j_k: f64,
    /// Delta S - Q/T_end [J/K], including the temporal discretization's entropy.
    /// This is not isolated chemical production or an external reservoir's entropy.
    pub entropy_residual_j_k: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct BoundaryGasReactionCandidate {
    gas: ClosedGasParcel,
    step: BoundaryGasReactionStep,
}
impl BoundaryGasReactionCandidate {
    pub const fn gas(&self) -> &ClosedGasParcel {
        &self.gas
    }
    pub const fn step(&self) -> &BoundaryGasReactionStep {
        &self.step
    }
    pub fn into_parts(self) -> (ClosedGasParcel, BoundaryGasReactionStep) {
        (self.gas, self.step)
    }
}

impl GasReaction {
    /// Solve the existing implicit mechanism at the actual proposed moving
    /// boundary endpoint. The caller owns the pressure/volume work law and its
    /// reciprocal mechanical debit. No reaction heat is added separately.
    ///
    /// Requires the original coupled implicit-mechanism policy; separate
    /// chemistry subcycling and scalar split policies refuse. This is a
    /// first-order homogeneous closed-gas endpoint, not a flame or shock model.
    /// Final species, source, energy and signed heat-entropy gates all apply.
    pub fn boundary_implicit_candidate(
        &self,
        gas: &ClosedGasParcel,
        duration_s: f64,
        boundary: GasReactionBoundary,
        controls: GasReactionSolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<BoundaryGasReactionCandidate, GasReactionError> {
        self.boundary_candidate(gas, duration_s, boundary, true, controls, cancelled)
    }

    // Only an enclosing native mechanical residual may use this trial path.
    // Its arbitrary pressure/work probes need not satisfy the final entropy
    // inequality. That owner MUST reapply entropy admission after its physical
    // pressure and work residuals close, before publishing any participant.
    pub(crate) fn mechanical_implicit_trial(
        &self,
        gas: &ClosedGasParcel,
        duration_s: f64,
        boundary: GasReactionBoundary,
        controls: GasReactionSolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<BoundaryGasReactionCandidate, GasReactionError> {
        self.boundary_candidate(gas, duration_s, boundary, false, controls, cancelled)
    }

    #[allow(clippy::too_many_arguments)]
    fn boundary_candidate(
        &self,
        gas: &ClosedGasParcel,
        duration_s: f64,
        boundary: GasReactionBoundary,
        admit_entropy: bool,
        controls: GasReactionSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BoundaryGasReactionCandidate, GasReactionError> {
        poll(&mut cancelled)?;
        controls.validate(duration_s)?;
        boundary.validate()?;
        if self.subcycling.is_some() || self.implicit.is_some() {
            return Err(invalid(
                "chemistry and boundary work require one coupled implicit interval",
            ));
        }
        let policy = self.implicit_mechanism.as_ref().ok_or_else(|| {
            invalid("moving-boundary chemistry requires the original implicit mechanism")
        })?;
        let candidate = implicit::coupled::candidate_at_boundary(
            self,
            gas,
            duration_s,
            boundary,
            admit_entropy,
            controls,
            policy,
            &mut cancelled,
        )?;
        let (next, reaction) = candidate.into_parts();
        let energy_residual_j =
            (reaction.stored_energy_change_j + boundary.work_by_gas_j) - boundary.heat_j;
        let endpoint_heat_entropy_j_k = boundary.heat_j / next.state().state().temperature;
        let entropy_residual_j_k = reaction.entropy_change_j_k - endpoint_heat_entropy_j_k;
        poll(&mut cancelled)?;
        Ok(BoundaryGasReactionCandidate {
            gas: next,
            step: BoundaryGasReactionStep {
                reaction,
                initial_volume_m3: gas.volume_m3(),
                boundary,
                energy_residual_j,
                endpoint_heat_entropy_j_k,
                entropy_residual_j_k,
            },
        })
    }
}
