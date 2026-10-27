//! Same-endpoint chemistry on the original finite electrical/thermal network.
//!
//! Every outer enthalpy trial resolves the prepared bulk/interface/radiative
//! sources, then solves each reacting chamber's chemistry, pressure work and
//! signed wall heats together. Unlisted chambers use the original midpoint
//! driver. The shared network owns heat deposition and the all-owner balance.

use super::gas_law::{GasSurface, GasSurfaceSet, close_radiation};
use super::*;
use crate::cavity_phs::piston::thermal::ConductanceSample;
use crate::cavity_phs::piston::thermal::boundary::ports::network::{
    Contact, DrivenNetworkCandidate,
    reactive::{self, CoupledBoundaryError, CoupledChamberBoundary, CoupledChamberReaction},
};
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;

/// Original circuit/thermal candidate with each chamber's actual time policy.
pub type CoupledMaterialExchangeCandidate<'m, I = ResolvedMaterialStatePoint> =
    DrivenNetworkCandidate<
        'm,
        ContactThermalSample<MaterialExchangeSample>,
        CoupledChamberBoundary<I>,
    >;

/// Material, reaction, pressure/heat and global-budget failures retain their owner.
pub type CoupledMaterialExchangeError<E = MaterialNetworkError> =
    CircuitEnthalpyError<DrivenNetworkError<MaterialExchangeError, CoupledBoundaryError<E>>>;

impl PreparedMaterialExchangeNetwork<'_, '_> {
    /// Connect same-endpoint chemistry to finite circuit energy, material cells,
    /// thermal/electrical contacts and the closed radiation enclosure.
    ///
    /// Reaction bindings select the coupled policy per original chamber. Their
    /// rates see the actual final composition/T/p after the SAME interval's
    /// material-wall heat and mechanical work. All other chambers retain their
    /// original nonreacting histories. Every numerical trial restarts from the
    /// accepted circuit charge/flux, cell enthalpy and complete gas inventory.
    ///
    /// Source selection, conductivity/resistivity, thermal contacts, emissivity,
    /// surface coverage and receiving-cell heat routes use the existing owners.
    /// A material change therefore reaches chemistry through physical heat and
    /// pressure feedback; no reaction heat is credited separately from gas U.
    #[allow(clippy::too_many_arguments)]
    pub fn coupled_reactive_step_candidate<'m>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        circuit: &MaterialCircuit,
        ports: &[f64],
        chambers: &[Chamber<'_>],
        reactions: &[CoupledChamberReaction<'_>],
        gas_patches: &[MaterialContact],
        controls: NetworkSolve,
        guess_j: &[f64],
    ) -> Result<CoupledMaterialExchangeCandidate<'m>, CoupledMaterialExchangeError> {
        let (surfaces, mapped) = conductive_surfaces(cx, gas_patches, controls)?;
        self.coupled_reactive_with_gas_law_candidate(
            cx,
            thermal,
            target_s,
            circuit,
            ports,
            chambers,
            reactions,
            &surfaces,
            controls,
            guess_j,
            |contact, state, source, gas| {
                self.resolve_conductive_gas_contact(mapped[&contact.id], state, source, gas)
            },
        )
    }

    /// Use an admitted custom gas heat law in the same complete assembly.
    ///
    /// Reacting contacts see the actual endpoint gas; ordinary contacts retain
    /// the native gas-energy midpoint. Both see the same proposed cell material
    /// snapshot. Return total conductance and its original typed source without
    /// advancing/debiting live state. Each evaluated contact heat is deposited
    /// exactly once, including opposing transfers at zero summed gas heat.
    ///
    /// A single result retains all gas endpoints and their original reaction or
    /// midpoint reports. Any source, reaction, pressure, global-work, balance or
    /// cancellation refusal publishes none. Independent walls, fixed circuit
    /// topology/C/L, fixed thermal geometry and homogeneous closed gas remain
    /// the supported scope. Temporal convergence requires whole-interval
    /// refinement; source/event discontinuities require interval boundaries.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn coupled_reactive_with_gas_law_candidate<'m, I, E>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        circuit: &MaterialCircuit,
        ports: &[f64],
        chambers: &[Chamber<'_>],
        reactions: &[CoupledChamberReaction<'_>],
        surfaces: &[GasSurface],
        controls: NetworkSolve,
        guess_j: &[f64],
        resolve: impl Fn(
            Contact,
            &DepositionSimulation<'m>,
            &MaterialExchangeSample,
            &ResolvedThermochemicalGasState,
        ) -> Result<ConductanceSample<I>, E>,
    ) -> Result<CoupledMaterialExchangeCandidate<'m, I>, CoupledMaterialExchangeError<E>> {
        let admitted = GasSurfaceSet::admit(cx, thermal, chambers, surfaces, controls)?;
        let binding = self.binding.electrical;
        let routes = binding
            .heat_destinations(thermal)
            .map_err(|e| coupled_error(e.into()))?;
        let destinations = binding
            .bulk
            .resistors
            .iter()
            .map(|r| r.destination)
            .collect::<Vec<_>>();
        let result = reactive::step_with_contacts_candidate(
            thermal,
            target_s,
            CircuitEnthalpyInput {
                circuit,
                ports,
                destinations: &destinations,
            },
            &routes,
            chambers,
            reactions,
            &admitted.contacts,
            controls,
            guess_j,
            |state| {
                admitted.check_trial(cx, state)?;
                self.resolve(cx, state)
            },
            |contact, state, source, gas| resolve(contact, state, source, gas),
            || cx.checkpoint().is_err(),
        );
        let mut candidate = match result {
            Err(CircuitEnthalpyError::Source(DrivenNetworkError::Material(e))) => {
                return Err(coupled_error(e));
            }
            other => other?,
        };
        admitted.check_samples(
            candidate
                .chambers()
                .iter()
                .flat_map(|chamber| chamber.boundary.source())
                .map(|sample| (sample.port.id, sample.conductance_w_k)),
        )?;
        let radiation_j = candidate
            .material()
            .material
            .radiation
            .as_ref()
            .map_or(0.0, |r| (target_s - thermal.time_s()) * r.net_heat_w());
        close_radiation(cx, &mut candidate.step, radiation_j, controls)?;
        Ok(candidate)
    }
}
