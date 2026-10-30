//! Original-snapshot nozzle transport followed by joint chemistry/wall heat.
//!
//! Both stages span ONE interval. Every flow sees the accepted initial gas
//! network; every chemistry/heat trial sees the same fully transferred network.
//! No trial withdraws mass twice, and no transported enthalpy is added as heat.

use super::*;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::network::snapshot::{
    self, SnapshotNozzleError, SnapshotNozzleLink, SnapshotNozzleStep,
};
use fs_material::thermochemical_gas::energy::parcel::transfer::snapshot::{
    GasSnapshotSolve, MAX_SNAPSHOT_PARCELS, MAX_SNAPSHOT_TRANSFERS,
};

/// Existing flow and joint heat admission, plus an enclosing energy allowance.
#[derive(Clone, Copy, Debug)]
pub struct OpenSolidGasSolve {
    pub flow: GasSnapshotSolve,
    pub thermal: SolidGasNetworkSolve,
    pub total_energy_tolerance_j: f64,
}

#[derive(Debug)]
pub enum OpenSolidGasError {
    Invalid(&'static str),
    Flow(SnapshotNozzleError),
    Thermal(SolidGasError),
    Cancelled,
}
impl std::fmt::Display for OpenSolidGasError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "open gas/shared material wall coupling refused: {self:?}"
        )
    }
}
impl std::error::Error for OpenSolidGasError {}
impl From<SnapshotNozzleError> for OpenSolidGasError {
    fn from(error: SnapshotNozzleError) -> Self {
        match error {
            SnapshotNozzleError::Cancelled => Self::Cancelled,
            other => Self::Flow(other),
        }
    }
}
impl From<SolidGasError> for OpenSolidGasError {
    fn from(error: SolidGasError) -> Self {
        match error {
            SolidGasError::Cancelled => Self::Cancelled,
            other => Self::Thermal(other),
        }
    }
}

/// Native stage accounts with a closed boundary around all compartments and
/// the finite wall. The individual compartments are open to each other.
#[derive(Clone, Debug)]
pub struct OpenSolidGasReport<S = MaterialNetworkSample> {
    pub flow: SnapshotNozzleStep,
    /// This stage starts after transfer. Its gas heat contains wall and authored
    /// heat only; incoming enthalpy is already present in the stage's inventory.
    pub thermal: SolidGasNetworkReport<S>,
    /// Fresh sum of actual final-minus-original gas masses [kg].
    pub total_mass_residual_kg: f64,
    /// Enclosing residual from the original datum-cancelled transport and
    /// thermal/chemical energy accounts [J]. Do not reconstruct this by
    /// subtracting rounded large formation-energy or enthalpy totals.
    pub total_energy_residual_j: f64,
}

/// All gas inventories, reaction memory, finite thermal history and local solid
/// sources publish together. No accepted input is changed by either stage.
pub struct OpenSolidGasCandidate<'m, 't, S = MaterialNetworkSample> {
    thermal: DepositionSimulation<'t>,
    gases: Vec<SolidGasEndpoint>,
    state: SolidEnthalpyState<'m>,
    report: OpenSolidGasReport<S>,
}
impl<'m, 't, S> OpenSolidGasCandidate<'m, 't, S> {
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn gases(&self) -> &[SolidGasEndpoint] {
        &self.gases
    }
    pub fn gas(&self, id: u64) -> Option<&IsochoricGasParcel> {
        self.gases
            .binary_search_by_key(&id, |gas| gas.id)
            .ok()
            .map(|i| &self.gases[i].gas)
    }
    pub const fn state(&self) -> &SolidEnthalpyState<'m> {
        &self.state
    }
    pub const fn report(&self) -> &OpenSolidGasReport<S> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        Vec<SolidGasEndpoint>,
        SolidEnthalpyState<'m>,
        OpenSolidGasReport<S>,
    ) {
        (self.thermal, self.gases, self.state, self.report)
    }
}

impl<'m> SolidEnthalpy<'m> {
    /// Transport mass/species/upstream enthalpy through a branching or cyclic
    /// nozzle network, then solve all admitted reactions and spatial wall heat.
    ///
    /// Every original flow is evaluated before ANY inventory changes. All
    /// reconstruction, gross withdrawal, source compatibility and pressure-order
    /// gates remain owned by `snapshot::step_candidate`. The existing shared
    /// wall residual then advances its original thermal history exactly once,
    /// using the actual transferred compositions, capacities and gas transport.
    /// Reactions use their original implicit mechanisms at the joint heat
    /// endpoint; formation energy stays in storage. Guess order remains sorted
    /// CellId followed by sorted gas ID, and guesses are heat-stage increments.
    ///
    /// This is first-order flow -> joint reaction/heat splitting over ONE
    /// duration, not simultaneous implicit flow. Wall/chemical pressure changes
    /// affect rates on the next interval; refine time and split at events.
    /// Compartments are well mixed and fixed-volume, with no external mass
    /// reservoir, jet thrust or spatial flame. Wall geometry and the original
    /// sequential thermoelastic interpretation remain unchanged. All native
    /// stages must succeed before the complete candidate is returned.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_open_gas_network<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        gases: &[SolidGasCompartment<'_>],
        links: &[SnapshotNozzleLink],
        target_s: f64,
        transport: Option<&PreparedMaterialNetwork<'_, '_>>,
        controls: OpenSolidGasSolve,
        initial_energy_increments_j: &[f64],
        cx: &Cx<'_>,
    ) -> Result<OpenSolidGasCandidate<'m, 't>, OpenSolidGasError> {
        self.step_open_gases(
            thermal,
            gases,
            links,
            target_s,
            transport,
            controls,
            initial_energy_increments_j,
            cx,
        )
    }

    /// The same transport/chemistry/wall transaction with material contact
    /// resistance, integrated layers and closed surface radiation through the
    /// original complete heat network. Electrical bindings still require their
    /// circuit owner; surface radiation assumes nonparticipating gas.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_open_gas_heat_network<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        gases: &[SolidGasCompartment<'_>],
        links: &[SnapshotNozzleLink],
        target_s: f64,
        transport: &PreparedHeatNetwork<'_, '_>,
        controls: OpenSolidGasSolve,
        initial_energy_increments_j: &[f64],
        cx: &Cx<'_>,
    ) -> Result<OpenSolidGasCandidate<'m, 't, HeatNetworkSample>, OpenSolidGasError> {
        self.step_open_gases(
            thermal,
            gases,
            links,
            target_s,
            transport,
            controls,
            initial_energy_increments_j,
            cx,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_open_gases<'t, T: GasHeatTransport>(
        &self,
        thermal: &DepositionSimulation<'t>,
        gases: &[SolidGasCompartment<'_>],
        links: &[SnapshotNozzleLink],
        target_s: f64,
        transport: T,
        controls: OpenSolidGasSolve,
        initial_energy_increments_j: &[f64],
        cx: &Cx<'_>,
    ) -> Result<OpenSolidGasCandidate<'m, 't, T::Sample>, OpenSolidGasError> {
        let checkpoint = || cx.checkpoint().map_err(|_| OpenSolidGasError::Cancelled);
        checkpoint()?;
        let dt = target_s - thermal.time_s();
        if !dt.is_finite()
            || dt <= 0.0
            || !target_s.is_finite()
            || gases.is_empty()
            || gases.len() > controls.thermal.maximum_gases
            || gases.len() > MAX_SNAPSHOT_PARCELS
            || links.len() > MAX_SNAPSHOT_TRANSFERS
            || !controls.total_energy_tolerance_j.is_finite()
            || controls.total_energy_tolerance_j < 0.0
            || initial_energy_increments_j.len()
                != thermal.cells().len().saturating_add(gases.len())
        {
            return Err(OpenSolidGasError::Invalid(
                "bounded complete gas/wall graph and one positive interval required",
            ));
        }
        // Let the authoritative timeline reject pending/interior events before
        // spending a flow stage. The heat owner repeats this admission itself.
        fs_thermal::enthalpy::backward_euler::BackwardEulerStep::new(thermal, target_s, || {
            cx.checkpoint().map_err(|_| DepositionError::Cancelled)
        })
        .map_err(SolidGasError::from)?;
        let mut initial = BTreeMap::new();
        let mut endpoints = BTreeMap::new();
        for input in gases {
            checkpoint()?;
            if initial
                .insert(input.id, ClosedGasParcel::from_isochoric(input.gas.clone()))
                .is_some()
            {
                return Err(OpenSolidGasError::Invalid(
                    "unique compartment IDs required",
                ));
            }
            endpoints.insert(input.id, controls.thermal.coupled.gas);
        }
        let flow =
            snapshot::step_candidate(&initial, links, dt, &endpoints, controls.flow, || {
                cx.checkpoint().is_err()
            })?;
        let (transferred, flow) = flow.into_parts();
        let transferred: BTreeMap<_, _> = transferred
            .into_iter()
            .map(|(id, gas)| (id, gas.into_isochoric()))
            .collect();
        let inputs: Vec<_> = gases
            .iter()
            .map(|input| SolidGasCompartment {
                id: input.id,
                gas: &transferred[&input.id],
                external_heat_j: input.external_heat_j,
                contacts: input.contacts,
                reaction: input.reaction,
            })
            .collect();
        let coupled = self.step_gas_network(
            thermal,
            &inputs,
            target_s,
            transport,
            controls.thermal,
            initial_energy_increments_j,
            cx,
        )?;
        let (thermal, gases, state, report) = coupled.into_parts();
        // Both original balances are endpoint-based and cancel conserved
        // species reference terms before rounding. Their intermediate state is
        // exactly the same `transferred` inventory, so those terms telescope.
        let residual = flow.exchange.total_energy_residual_j + report.combined_energy_residual_j;
        let (mut mass, mut correction) = (0.0, 0.0);
        for endpoint in &gases {
            checkpoint()?;
            let initial_temperature = initial[&endpoint.id].state().state().temperature;
            if [
                transferred[&endpoint.id].state().state().temperature,
                endpoint.gas.state().state().temperature,
            ]
            .iter()
            .any(|temperature| {
                (*temperature - initial_temperature).abs()
                    > controls.thermal.coupled.maximum_gas_temperature_increment_k
            }) {
                return Err(OpenSolidGasError::Invalid(
                    "gas transfer or final endpoint exceeds the temperature-increment allowance",
                ));
            }
            let term = endpoint.gas.mass_kg() - initial[&endpoint.id].mass_kg();
            let adjusted = term - correction;
            let next = mass + adjusted;
            correction = (next - mass) - adjusted;
            mass = next;
            if ![term, adjusted, mass, correction]
                .iter()
                .all(|x| x.is_finite())
            {
                return Err(OpenSolidGasError::Invalid(
                    "final gas mass balance is unrepresentable",
                ));
            }
        }
        if !residual.is_finite()
            || residual.abs() > controls.total_energy_tolerance_j
            || mass.abs() > controls.flow.mass_tolerance_kg
        {
            return Err(OpenSolidGasError::Invalid(
                "flow, reacting gases and finite wall fail the enclosing balance",
            ));
        }
        checkpoint()?;
        Ok(OpenSolidGasCandidate {
            thermal,
            gases,
            state,
            report: OpenSolidGasReport {
                flow,
                thermal: report,
                total_mass_residual_kg: mass,
                total_energy_residual_j: residual,
            },
        })
    }
}
