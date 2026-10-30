//! Chemistry and spatial heat transport at ONE backward-Euler endpoint.
//!
//! Each outer thermal residual solves the original implicit chemical mechanism
//! at its proposed total gas energy. The original gas and every spatial cell
//! then exchange heat at that SAME composition, pressure and temperature. This
//! is nonlinear elimination of the existing chemical equations, not a new
//! reaction/thermal solver or an adiabatic reaction-then-cooling split.

use super::super::transport::GasHeatTransport;
use super::super::{
    Cx, DepositionSimulation, HeatNetworkSample, IsochoricGasParcel, MaterialNetworkSample,
    PreparedHeatNetwork, PreparedMaterialNetwork, SolidEnthalpy, SolidEnthalpyState,
    SolidGasCandidate, SolidGasContact, SolidGasError, SolidGasReport, SolidGasSolve,
};
use crate::cavity_phs::piston::reactive::{
    GasReaction, GasReactionSolve, heat::HeatedGasReactionStep,
};
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;

/// The original accepted thermal candidate and its identical chemical endpoint.
/// There is no second temperature, composition, clock or energy store.
pub struct CoupledReactiveSolidGasCandidate<'m, 't, S = MaterialNetworkSample> {
    coupled: SolidGasCandidate<'m, 't, S>,
    reaction: HeatedGasReactionStep,
}
impl<'m, 't, S> CoupledReactiveSolidGasCandidate<'m, 't, S> {
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        self.coupled.thermal()
    }
    pub const fn gas(&self) -> &IsochoricGasParcel {
        self.coupled.gas()
    }
    pub const fn state(&self) -> &SolidEnthalpyState<'m> {
        self.coupled.state()
    }
    /// The original local cell, patch, gas and combined energy balances.
    pub const fn report(&self) -> &SolidGasReport<S> {
        self.coupled.report()
    }
    /// Actual endpoint channel equations, source evaluations and heat/entropy.
    pub const fn reaction(&self) -> &HeatedGasReactionStep {
        &self.reaction
    }
    pub fn into_parts(self) -> (SolidGasCandidate<'m, 't, S>, HeatedGasReactionStep) {
        (self.coupled, self.reaction)
    }
}

impl<'m> SolidEnthalpy<'m> {
    /// Resolve chemistry, gas heat and local solid enthalpy simultaneously.
    ///
    /// The unknowns remain the original cell energy increments plus total gas
    /// energy increment. At each trial, `heated_implicit_candidate` solves all
    /// declared reaction channels using that final energy; their final T/p/x
    /// determine the original material films and spatial thermal residuals.
    /// Chemical formation energy remains inside gas storage. Only authored
    /// external heat enters the enclosing supply; patch heats cancel once.
    /// Material k(T) and the original solid elasticity/expansion sources are
    /// re-evaluated at every actual proposed cell endpoint.
    ///
    /// Requires the original `with_implicit_mechanism` policy. Chemistry-only
    /// subcycling refuses: refine the complete coupled interval instead. Inner
    /// extent residuals and outer cell/gas energy residuals are independently
    /// admitted. `maximum_trials` bounds all outer probes; each inner solve
    /// retains its original evaluation/caloric caps. This nested solve can be
    /// expensive and provides no uniqueness or unconditional convergence claim.
    ///
    /// The gas remains homogeneous, closed and fixed-volume. Solid geometry and
    /// its reference thermal chart remain fixed; the returned mechanical field
    /// retains the original sequential thermoelastic interpretation. This adds
    /// no flame front, pressure traction, reciprocal thermoelastic heat or moving
    /// interface. Backward Euler remains first order. Any stage/domain/work or
    /// late cancellation refusal publishes no gas, thermal or mechanical state.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_coupled_reacting_isochoric_gas<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        gas: &IsochoricGasParcel,
        reaction: &GasReaction,
        target_s: f64,
        external_gas_heat_j: f64,
        contacts: &[SolidGasContact],
        transport: Option<&PreparedMaterialNetwork<'_, '_>>,
        reaction_controls: GasReactionSolve,
        controls: SolidGasSolve,
        initial_energy_increments_j: &[f64],
        cx: &Cx<'_>,
    ) -> Result<CoupledReactiveSolidGasCandidate<'m, 't>, SolidGasError> {
        self.step_reacting_gas_with_transport(
            thermal,
            gas,
            reaction,
            target_s,
            external_gas_heat_j,
            contacts,
            transport,
            reaction_controls,
            controls,
            initial_energy_increments_j,
            cx,
        )
    }

    /// The same joint implicit chemistry/gas/solid endpoint with the complete
    /// shared heat network. Source contact resistance, integrated material layers
    /// and closed surface radiation respond to the SAME cell temperatures that
    /// exchange energy with the reacting gas. Chemical formation energy stays in
    /// the original gas inventory and is not added again as a heat source.
    ///
    /// All original mechanism, source, numerical and fixed-volume restrictions
    /// remain. Electrical graph bindings refuse; surface radiation assumes a
    /// nonparticipating gas and does not model volumetric radiative chemistry.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_coupled_reacting_isochoric_gas_heat_network<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        gas: &IsochoricGasParcel,
        reaction: &GasReaction,
        target_s: f64,
        external_gas_heat_j: f64,
        contacts: &[SolidGasContact],
        transport: &PreparedHeatNetwork<'_, '_>,
        reaction_controls: GasReactionSolve,
        controls: SolidGasSolve,
        initial_energy_increments_j: &[f64],
        cx: &Cx<'_>,
    ) -> Result<CoupledReactiveSolidGasCandidate<'m, 't, HeatNetworkSample>, SolidGasError> {
        self.step_reacting_gas_with_transport(
            thermal,
            gas,
            reaction,
            target_s,
            external_gas_heat_j,
            contacts,
            transport,
            reaction_controls,
            controls,
            initial_energy_increments_j,
            cx,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_reacting_gas_with_transport<'t, T: GasHeatTransport>(
        &self,
        thermal: &DepositionSimulation<'t>,
        gas: &IsochoricGasParcel,
        reaction: &GasReaction,
        target_s: f64,
        external_gas_heat_j: f64,
        contacts: &[SolidGasContact],
        transport: T,
        reaction_controls: GasReactionSolve,
        controls: SolidGasSolve,
        initial_energy_increments_j: &[f64],
        cx: &Cx<'_>,
    ) -> Result<CoupledReactiveSolidGasCandidate<'m, 't, T::Sample>, SolidGasError> {
        cx.checkpoint().map_err(|_| SolidGasError::Cancelled)?;
        if reaction_controls.gas != controls.gas {
            return Err(SolidGasError::Invalid(
                "coupled chemistry and heat must share one gas energy-inversion policy",
            ));
        }
        let original = ClosedGasParcel::from_isochoric(gas.clone());
        let duration_s = target_s - thermal.time_s();
        let (coupled, reaction) = self.step_gas_with_endpoint(
            thermal,
            gas,
            target_s,
            external_gas_heat_j,
            contacts,
            transport,
            controls,
            initial_energy_increments_j,
            cx,
            |heat_j| {
                let candidate = reaction.heated_implicit_candidate(
                    &original,
                    duration_s,
                    heat_j,
                    reaction_controls,
                    || cx.checkpoint().is_err(),
                )?;
                let (gas, report) = candidate.into_parts();
                Ok((gas.into_isochoric(), report.heat, report))
            },
        )?;
        cx.checkpoint().map_err(|_| SolidGasError::Cancelled)?;
        Ok(CoupledReactiveSolidGasCandidate { coupled, reaction })
    }
}
