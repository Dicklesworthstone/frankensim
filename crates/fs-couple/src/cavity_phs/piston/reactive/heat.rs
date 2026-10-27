//! Nonadiabatic endpoint chemistry on the original conservative reaction owner.
//!
//! The SAME backward-Euler mechanism equations use the actual endpoint at
//! U_end = U_start + Q, fixed mass and fixed volume. Formation energies remain
//! in U; reaction heat is never added as an independent source. The caller owns
//! the signed heat exchange and its reciprocal thermal participant.
//! Reference: Cantera ideal-gas control-volume energy/species equations,
//! <https://cantera.org/stable/reference/reactors/ideal-gas-reactor.html>.

use super::*;
use fs_material::thermochemical_gas::energy::GasHeatStep;

/// Original chemical equations and the distinct external heat account.
#[derive(Clone, Debug, PartialEq)]
pub struct HeatedGasReactionStep {
    pub reaction: GasReactionStep,
    pub heat: GasHeatStep,
    /// Signed Q/T_end [J/K] at the actual uniform gas endpoint, not the entropy
    /// supplied by an unspecified external reservoir at a different temperature.
    pub endpoint_heat_entropy_j_k: f64,
    /// Delta S - Q/T_end [J/K], independently admitted nonnegative up to the
    /// original tolerance. This includes backward-Euler numerical entropy and
    /// is not an isolated chemical entropy-production measurement.
    pub entropy_residual_j_k: f64,
}

/// Immutable complete gas/chemistry/heat proposal.
#[derive(Clone, Debug, PartialEq)]
pub struct HeatedGasReactionCandidate {
    gas: ClosedGasParcel,
    step: HeatedGasReactionStep,
}
impl HeatedGasReactionCandidate {
    pub const fn gas(&self) -> &ClosedGasParcel {
        &self.gas
    }
    pub const fn step(&self) -> &HeatedGasReactionStep {
        &self.step
    }
    pub fn into_parts(self) -> (ClosedGasParcel, HeatedGasReactionStep) {
        (self.gas, self.step)
    }
}

impl GasReaction {
    /// Solve every declared endpoint rate with a signed heat increment inside
    /// the SAME source-backed mass/composition/energy reconstruction.
    ///
    /// Requires `with_implicit_mechanism`, including for a single channel.
    /// Adiabatic/scalar split selectors and chemistry-only subcycling refuse:
    /// their intermediate states do not represent a coupled thermal trajectory.
    /// The enclosing owner must refine chemistry AND heat on the same interval.
    /// Every rate (including material Arrhenius/collider/Troe response) uses the
    /// actual heated/cooled endpoint T, partial pressures and composition.
    ///
    /// This is first-order homogeneous fixed-volume chemistry, with the original
    /// stoichiometry, source catalog, consumed-fraction and work limits. Heat is
    /// external to this gas account, not additional formation-energy release.
    /// The endpoint entropy residual is checked after subtracting Q/T_end, so
    /// cooling may lower gas entropy without being mislabeled a chemical fault.
    /// No accepted input changes on source, numerical or cancellation refusal.
    pub fn heated_implicit_candidate(
        &self,
        gas: &ClosedGasParcel,
        duration_s: f64,
        heat_j: f64,
        controls: GasReactionSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HeatedGasReactionCandidate, GasReactionError> {
        poll(&mut cancelled)?;
        controls.validate(duration_s)?;
        if !heat_j.is_finite() || self.subcycling.is_some() || self.implicit.is_some() {
            return Err(invalid(
                "finite heat and one shared implicit reaction/thermal interval required",
            ));
        }
        let policy = self.implicit_mechanism.as_ref().ok_or_else(|| {
            invalid("heated chemistry requires the original implicit mechanism policy")
        })?;
        let candidate = implicit::coupled::candidate_with_heat(
            self,
            gas,
            duration_s,
            heat_j,
            controls,
            policy,
            &mut cancelled,
        )?;
        let (gas, reaction) = candidate.into_parts();
        let heat = GasHeatStep {
            supplied_heat_j: heat_j,
            stored_energy_change_j: reaction.stored_energy_change_j,
            residual_j: reaction.stored_energy_change_j - heat_j,
            iterations: reaction
                .implicit_mechanism
                .as_ref()
                .expect("implicit endpoint retains its report")
                .endpoint_energy_iterations,
        };
        let endpoint_heat_entropy_j_k = heat_j / gas.state().state().temperature;
        let entropy_residual_j_k = reaction.entropy_change_j_k - endpoint_heat_entropy_j_k;
        poll(&mut cancelled)?;
        Ok(HeatedGasReactionCandidate {
            gas,
            step: HeatedGasReactionStep {
                reaction,
                heat,
                endpoint_heat_entropy_j_k,
                entropy_residual_j_k,
            },
        })
    }
}
