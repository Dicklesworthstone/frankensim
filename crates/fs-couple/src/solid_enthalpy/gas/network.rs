//! Separate finite gas inventories sharing one spatial material wall.
//!
//! All compartments and wall cells use the original joint backward-Euler
//! residual. A wall cell can exchange heat through distinct patches with
//! several gases without advancing its clock or debiting its storage twice.

use super::*;
use crate::cavity_phs::piston::reactive::{
    GasReaction, GasReactionSolve, heat::HeatedGasReactionStep,
};
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;

/// Conservative inter-compartment flow with the same reacting spatial wall.
pub mod exchange;

/// One closed, well-mixed, fixed-volume compartment. Identity labels the
/// physical inventory, not a material name. Different compartments may use
/// different gas catalogs, compositions, volumes and chemical mechanisms.
#[derive(Clone, Copy)]
pub struct SolidGasCompartment<'a> {
    pub id: u64,
    pub gas: &'a IsochoricGasParcel,
    pub external_heat_j: f64,
    /// Empty for an adiabatic reservoir; the complete network must contain at
    /// least one wall patch. Such a reservoir still exchanges species/enthalpy
    /// when the original nozzle-network owner is selected.
    pub contacts: &'a [SolidGasContact],
    /// None freezes composition. Some solves the original implicit mechanism
    /// at this compartment's actual endpoint energy in every thermal trial.
    pub reaction: Option<(&'a GasReaction, GasReactionSolve)>,
}

/// The existing work/energy policy plus an explicit compartment count bound.
#[derive(Clone, Copy, Debug)]
pub struct SolidGasNetworkSolve {
    pub maximum_gases: usize,
    /// maximum_contacts bounds the TOTAL number of patches across all gases.
    pub coupled: SolidGasSolve,
}

/// One accepted inventory, ordered by its stable physical identifier.
pub struct SolidGasEndpoint {
    pub id: u64,
    pub gas: IsochoricGasParcel,
}

/// Separate caloric balance and exact chemical payload for one compartment.
#[derive(Clone, Debug)]
pub struct SolidGasCompartmentReport<P = Option<HeatedGasReactionStep>> {
    pub id: u64,
    pub gas: GasHeatStep,
    pub contacts: Vec<SolidGasContactSample>,
    pub external_heat_j: f64,
    pub solid_to_gas_heat_j: f64,
    pub energy_residual_j: f64,
    pub endpoint: P,
}

/// One wall account, with independent gas balances. Internal patch heat
/// cancels once; formation energy remains inside each original gas inventory.
#[derive(Clone, Debug)]
pub struct SolidGasNetworkReport<S = MaterialNetworkSample, P = Option<HeatedGasReactionStep>> {
    pub thermal: EnergyAudit,
    pub gases: Vec<SolidGasCompartmentReport<P>>,
    pub transport: Option<S>,
    pub closed_radiation_net_j: f64,
    pub external_thermal_heat_j: f64,
    pub combined_energy_residual_j: f64,
    pub maximum_cell_residual_j: f64,
    pub nonlinear: NewtonReport,
    pub trials: usize,
}

/// Detached simultaneous endpoint. No gas, chemistry or wall state is
/// published if any compartment, cell, source, numerical or work gate refuses.
pub struct SolidGasNetworkCandidate<
    'm,
    't,
    S = MaterialNetworkSample,
    P = Option<HeatedGasReactionStep>,
> {
    pub(super) thermal: DepositionSimulation<'t>,
    pub(super) gases: Vec<SolidGasEndpoint>,
    pub(super) state: SolidEnthalpyState<'m>,
    pub(super) report: SolidGasNetworkReport<S, P>,
}
impl<'m, 't, S, P> SolidGasNetworkCandidate<'m, 't, S, P> {
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn gases(&self) -> &[SolidGasEndpoint] {
        &self.gases
    }
    pub fn gas(&self, id: u64) -> Option<&IsochoricGasParcel> {
        self.gases
            .binary_search_by_key(&id, |entry| entry.id)
            .ok()
            .map(|index| &self.gases[index].gas)
    }
    pub const fn state(&self) -> &SolidEnthalpyState<'m> {
        &self.state
    }
    pub const fn report(&self) -> &SolidGasNetworkReport<S, P> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        Vec<SolidGasEndpoint>,
        SolidEnthalpyState<'m>,
        SolidGasNetworkReport<S, P>,
    ) {
        (self.thermal, self.gases, self.state, self.report)
    }
}

/// Canonically ordered inputs to the original gas/solid residual. Source
/// chemistry is provided by its endpoint callback, never by a second solve.
pub(super) struct GasEndpointInput<'a> {
    pub id: u64,
    pub gas: &'a IsochoricGasParcel,
    pub external_heat_j: f64,
    pub contacts: &'a [SolidGasContact],
}

impl<'m> SolidEnthalpy<'m> {
    /// Advance all finite compartments and the shared wall at ONE endpoint.
    ///
    /// Initial guesses are sorted CellId energy increments followed by gas
    /// increments sorted by compartment ID, independent of input enumeration.
    /// Patch IDs must be unique across the whole network; their aggregate area
    /// must fit each cell's exposed area. Geometry must describe disjoint patches.
    /// Every compartment retains its own mass, volume, material catalog and
    /// physical query coordinates. No gas/species exchange is implied.
    ///
    /// Chemical participants require the existing implicit-mechanism policy
    /// and the same caloric inversion controls as the heat solve. Source-dependent
    /// rates and product transport see the actual joint endpoint. The wall's
    /// original structural source field is returned for native deformation.
    /// This remains fixed thermal geometry and sequential thermoelastic loading;
    /// no pressure traction, moving wall, flame, reciprocal thermoelastic heat
    /// or whole-system entropy theorem is introduced.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_isochoric_gases<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        gases: &[SolidGasCompartment<'_>],
        target_s: f64,
        transport: Option<&PreparedMaterialNetwork<'_, '_>>,
        controls: SolidGasNetworkSolve,
        initial_energy_increments_j: &[f64],
        cx: &Cx<'_>,
    ) -> Result<SolidGasNetworkCandidate<'m, 't>, SolidGasError> {
        self.step_gas_network(
            thermal,
            gases,
            target_s,
            transport,
            controls,
            initial_energy_increments_j,
            cx,
        )
    }

    /// The identical joint endpoint with the complete original material heat
    /// network, including thermal contacts/layers and closed surface radiation.
    /// Electrical bindings still refuse because no circuit is advanced here.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_isochoric_gases_heat_network<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        gases: &[SolidGasCompartment<'_>],
        target_s: f64,
        transport: &PreparedHeatNetwork<'_, '_>,
        controls: SolidGasNetworkSolve,
        initial_energy_increments_j: &[f64],
        cx: &Cx<'_>,
    ) -> Result<SolidGasNetworkCandidate<'m, 't, HeatNetworkSample>, SolidGasError> {
        self.step_gas_network(
            thermal,
            gases,
            target_s,
            transport,
            controls,
            initial_energy_increments_j,
            cx,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_gas_network<'t, T: GasHeatTransport>(
        &self,
        thermal: &DepositionSimulation<'t>,
        gases: &[SolidGasCompartment<'_>],
        target_s: f64,
        transport: T,
        controls: SolidGasNetworkSolve,
        initial_energy_increments_j: &[f64],
        cx: &Cx<'_>,
    ) -> Result<SolidGasNetworkCandidate<'m, 't, T::Sample>, SolidGasError> {
        cx.checkpoint().map_err(|_| SolidGasError::Cancelled)?;
        if gases.is_empty() || gases.len() > controls.maximum_gases {
            return Err(SolidGasError::Invalid(
                "bounded nonempty gas network required",
            ));
        }
        let mut ordered: Vec<_> = gases.iter().collect();
        ordered.sort_by_key(|entry| entry.id);
        for entry in &ordered {
            if entry
                .reaction
                .is_some_and(|(_, policy)| policy.gas != controls.coupled.gas)
            {
                return Err(SolidGasError::Invalid(
                    "chemistry and heat must share one gas energy-inversion policy",
                ));
            }
        }
        let inputs: Vec<_> = ordered
            .iter()
            .map(|entry| GasEndpointInput {
                id: entry.id,
                gas: entry.gas,
                external_heat_j: entry.external_heat_j,
                contacts: entry.contacts,
            })
            .collect();
        let dt = target_s - thermal.time_s();
        self.step_gases_with_endpoint(
            thermal,
            &inputs,
            target_s,
            transport,
            controls.coupled,
            initial_energy_increments_j,
            cx,
            |index, heat| {
                let input = ordered[index];
                match input.reaction {
                    Some((reaction, policy)) => {
                        let initial = ClosedGasParcel::from_isochoric(input.gas.clone());
                        let candidate = reaction.heated_implicit_candidate(
                            &initial,
                            dt,
                            heat,
                            policy,
                            || cx.checkpoint().is_err(),
                        )?;
                        let (gas, report) = candidate.into_parts();
                        Ok((gas.into_isochoric(), report.heat, Some(report)))
                    }
                    None => {
                        let mut gas = input.gas.clone();
                        let report =
                            gas.add_heat(heat, controls.coupled.gas, || cx.checkpoint().is_err())?;
                        Ok((gas, report, None))
                    }
                }
            },
        )
    }
}
