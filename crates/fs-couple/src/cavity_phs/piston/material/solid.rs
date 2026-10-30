//! Source-bound moving chambers on the original shared spatial-wall transaction.
//!
//! Original-snapshot flow and joint chemistry/wall heat precede native midpoint
//! pressure work over ONE interval. The returned moving inventories become the
//! next interval's flow snapshot; no solver or independent gas state is added.

use super::{MaterialGasPiston, PistonError, PistonInput, PistonSolve, PistonStep};
use crate::cavity_phs::piston::reactive::{GasReaction, GasReactionSolve};
use crate::circuit_enthalpy::contacts::heat::{HeatNetworkSample, PreparedHeatNetwork};
use crate::circuit_enthalpy::material::{MaterialNetworkSample, PreparedMaterialNetwork};
use crate::port_adoption::{PortTransfer, StepAccount};
use crate::solid_enthalpy::gas::SolidGasContact;
use crate::solid_enthalpy::gas::network::SolidGasCompartment;
use crate::solid_enthalpy::gas::network::exchange::{
    OpenSolidGasCandidate, OpenSolidGasError, OpenSolidGasReport, OpenSolidGasSolve,
};
use crate::solid_enthalpy::{SolidEnthalpy, SolidEnthalpyState};
use fs_exec::Cx;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::network::snapshot::SnapshotNozzleLink;
use fs_material::thermochemical_gas::energy::parcel::transfer::snapshot::MAX_SNAPSHOT_PARCELS;
use fs_thermal::enthalpy::DepositionSimulation;
use std::collections::BTreeMap;

/// One accepted source-bound chamber and its loads for the next interval.
/// Fixed support represents a stationary vessel using the same native owner.
#[derive(Clone, Copy)]
pub struct SolidPistonCompartment<'a> {
    pub id: u64,
    pub piston: &'a MaterialGasPiston,
    pub external_heat_j: f64,
    pub contacts: &'a [SolidGasContact],
    pub reaction: Option<(&'a GasReaction, GasReactionSolve)>,
    pub exterior_pressure_pa: f64,
    pub external_force_n: f64,
    /// Its caloric controls must equal the shared heat owner's controls.
    /// Duration comes only from the thermal timeline, never from another input.
    pub motion: PistonSolve,
}

#[derive(Debug)]
pub enum MovingSolidGasError {
    Invalid(&'static str),
    Heat(OpenSolidGasError),
    Piston { id: u64, error: PistonError },
    Cancelled,
}
impl std::fmt::Display for MovingSolidGasError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "moving gas/shared spatial wall refused: {self:?}")
    }
}
impl std::error::Error for MovingSolidGasError {}
impl From<OpenSolidGasError> for MovingSolidGasError {
    fn from(error: OpenSolidGasError) -> Self {
        match error {
            OpenSolidGasError::Cancelled => Self::Cancelled,
            other => Self::Heat(other),
        }
    }
}
fn mechanical(id: u64, error: PistonError) -> MovingSolidGasError {
    match error {
        PistonError::Cancelled => MovingSolidGasError::Cancelled,
        error => MovingSolidGasError::Piston { id, error },
    }
}

/// Native pressure-work account plus the complete original-to-final gas change.
#[derive(Clone, Debug)]
pub struct SolidPistonStep {
    pub id: u64,
    pub motion: PistonStep,
    /// Sum of original source-bearing transport, heat/reaction and work storage.
    /// This is diagnostic; use the datum-cancelled total balance for admission.
    pub gas_energy_change_j: f64,
}

#[derive(Clone, Debug)]
pub struct MovingSolidGasReport<S = MaterialNetworkSample> {
    pub heat: OpenSolidGasReport<S>,
    /// Canonical compartment order. The heat stage ends where motion begins.
    pub pistons: Vec<SolidPistonStep>,
    /// Actual final-minus-original masses, after all moving endpoints exist.
    pub total_mass_residual_kg: f64,
    /// Native telescoping balance across the ORIGINAL and FINAL gas, finite
    /// wall and kinetic stores, including exterior pressure and force work.
    /// Transport enthalpy, reaction storage and wall heat are not added again.
    pub total_energy_residual_j: f64,
}

/// All source bodies, accepted gas/velocity/support and port histories are
/// retained. Returning a candidate never changes any accepted input owner.
pub struct MovingSolidGasCandidate<'m, 't, S = MaterialNetworkSample> {
    thermal: DepositionSimulation<'t>,
    pistons: BTreeMap<u64, MaterialGasPiston>,
    state: SolidEnthalpyState<'m>,
    report: MovingSolidGasReport<S>,
}
impl<'m, 't, S> MovingSolidGasCandidate<'m, 't, S> {
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn pistons(&self) -> &BTreeMap<u64, MaterialGasPiston> {
        &self.pistons
    }
    pub fn piston(&self, id: u64) -> Option<&MaterialGasPiston> {
        self.pistons.get(&id)
    }
    pub const fn state(&self) -> &SolidEnthalpyState<'m> {
        &self.state
    }
    pub const fn report(&self) -> &MovingSolidGasReport<S> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        BTreeMap<u64, MaterialGasPiston>,
        SolidEnthalpyState<'m>,
        MovingSolidGasReport<S>,
    ) {
        (self.thermal, self.pistons, self.state, self.report)
    }
}

// Compensated canonical sums combine native accounts, never absolute formation
// energy totals. Each stage begins from the exact preceding physical endpoint.
fn sum(values: impl IntoIterator<Item = f64>) -> Result<f64, MovingSolidGasError> {
    let (mut value, mut correction) = (0.0, 0.0);
    for term in values {
        let adjusted = term - correction;
        let next = value + adjusted;
        correction = (next - value) - adjusted;
        value = next;
        if ![term, adjusted, value, correction]
            .iter()
            .all(|v| v.is_finite())
        {
            return Err(MovingSolidGasError::Invalid(
                "extensive sum is unrepresentable",
            ));
        }
    }
    Ok(value)
}

impl<'m> SolidEnthalpy<'m> {
    /// Advance shared nozzle flow, reactions and finite spatial-wall heat, then
    /// the original material-bound gas/piston momentum and pressure-work owner.
    ///
    /// This is explicit first-order flow -> joint reaction/heat -> motion over
    /// ONE physical duration. Every heat trial uses the same transferred gas;
    /// every mechanical trial uses the accepted heat endpoint. The final volume,
    /// species and pressure drive the NEXT interval's flow and heat. Piston mass
    /// comes from its original source body and never changes during heating.
    ///
    /// Gas-wall contact areas/layers and the structural reference mesh remain
    /// fixed. The translating piston is a separate rigid source body, not a
    /// pressure-loaded/deforming tetrahedron. Returned structural sources retain
    /// sequential thermal-loading semantics; this adds no reciprocal elastic
    /// heat, jet thrust, crank/contact/friction law or spatial combustion model.
    /// Refine time and split at events. All native admissions must pass before
    /// heat history, material pistons and one complete port tick return together.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_material_pistons<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        gases: &[SolidPistonCompartment<'_>],
        links: &[SnapshotNozzleLink],
        target_s: f64,
        transport: Option<&PreparedMaterialNetwork<'_, '_>>,
        controls: OpenSolidGasSolve,
        initial_energy_increments_j: &[f64],
        cx: &Cx<'_>,
    ) -> Result<MovingSolidGasCandidate<'m, 't>, MovingSolidGasError> {
        self.step_material_pistons(thermal, gases, target_s, controls, cx, |inputs| {
            self.step_with_open_gas_network(
                thermal,
                inputs,
                links,
                target_s,
                transport,
                controls,
                initial_energy_increments_j,
                cx,
            )
        })
    }

    /// The same continuing moving assembly with original material contact,
    /// layer and closed surface-radiation heat paths. Electrical bindings still
    /// require their circuit owner; gas is nonparticipating in surface radiation.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_material_pistons_heat_network<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        gases: &[SolidPistonCompartment<'_>],
        links: &[SnapshotNozzleLink],
        target_s: f64,
        transport: &PreparedHeatNetwork<'_, '_>,
        controls: OpenSolidGasSolve,
        initial_energy_increments_j: &[f64],
        cx: &Cx<'_>,
    ) -> Result<MovingSolidGasCandidate<'m, 't, HeatNetworkSample>, MovingSolidGasError> {
        self.step_material_pistons(thermal, gases, target_s, controls, cx, |inputs| {
            self.step_with_open_gas_heat_network(
                thermal,
                inputs,
                links,
                target_s,
                transport,
                controls,
                initial_energy_increments_j,
                cx,
            )
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn step_material_pistons<'t, S>(
        &self,
        thermal: &DepositionSimulation<'t>,
        inputs: &[SolidPistonCompartment<'_>],
        target_s: f64,
        controls: OpenSolidGasSolve,
        cx: &Cx<'_>,
        heat_step: impl FnOnce(
            &[SolidGasCompartment<'_>],
        ) -> Result<OpenSolidGasCandidate<'m, 't, S>, OpenSolidGasError>,
    ) -> Result<MovingSolidGasCandidate<'m, 't, S>, MovingSolidGasError> {
        let checkpoint = || cx.checkpoint().map_err(|_| MovingSolidGasError::Cancelled);
        checkpoint()?;
        let dt = target_s - thermal.time_s();
        if !dt.is_finite()
            || dt <= 0.0
            || !target_s.is_finite()
            || inputs.is_empty()
            || inputs.len() > controls.thermal.maximum_gases
            || inputs.len() > MAX_SNAPSHOT_PARCELS
        {
            return Err(MovingSolidGasError::Invalid(
                "one positive interval and a bounded nonempty chamber graph required",
            ));
        }
        let mut ordered = BTreeMap::new();
        for input in inputs {
            checkpoint()?;
            if input.motion.gas != controls.thermal.coupled.gas
                || !input.exterior_pressure_pa.is_finite()
                || input.exterior_pressure_pa < 0.0
                || !input.external_force_n.is_finite()
                || ordered.insert(input.id, input).is_some()
            {
                return Err(MovingSolidGasError::Invalid(
                    "unique chambers, finite loads and one caloric admission required",
                ));
            }
        }
        let parcels: BTreeMap<_, _> = ordered
            .iter()
            .map(|(&id, input)| (id, input.piston.gas().clone().into_isochoric()))
            .collect();
        let heat_inputs: Vec<_> = ordered
            .iter()
            .map(|(&id, input)| SolidGasCompartment {
                id,
                gas: &parcels[&id],
                external_heat_j: input.external_heat_j,
                contacts: input.contacts,
                reaction: input.reaction,
            })
            .collect();
        let (thermal, heated, state, heat) = heat_step(&heat_inputs)?.into_parts();
        let mut pistons = BTreeMap::new();
        let mut reports = Vec::with_capacity(heated.len());
        let mut mass_changes = Vec::with_capacity(heated.len());
        let mut residuals = vec![heat.total_energy_residual_j];
        for ((endpoint, gas_heat), transfer) in heated
            .into_iter()
            .zip(&heat.thermal.gases)
            .zip(&heat.flow.exchange.parcels)
        {
            checkpoint()?;
            let id = endpoint.id;
            if gas_heat.id != id || transfer.id != id {
                return Err(MovingSolidGasError::Invalid(
                    "native compartment accounts disagree",
                ));
            }
            let input = ordered[&id];
            let old = input.piston;
            let heat_temperature = endpoint.gas.state().state().temperature;
            // Clone the actual accepted mechanical owner; constructing a fresh
            // piston would lose momentum, support and its preceding port clock.
            let mut next = old.clone();
            next.piston.gas = ClosedGasParcel::from_isochoric(endpoint.gas);
            let motion = next
                .step(
                    PistonInput {
                        duration_s: dt,
                        heat_j: 0.0,
                        exterior_pressure_pa: input.exterior_pressure_pa,
                        external_force_n: input.external_force_n,
                    },
                    input.motion,
                    || cx.checkpoint().is_err(),
                )
                .map_err(|e| mechanical(id, e))?;
            if (next.gas().state().state().temperature - old.gas().state().state().temperature)
                .abs()
                > controls.thermal.coupled.maximum_gas_temperature_increment_k
            {
                return Err(MovingSolidGasError::Invalid(
                    "moving gas exceeds temperature-increment allowance",
                ));
            }
            let gas_change = sum([
                transfer.energy_change_j,
                gas_heat.gas.stored_energy_change_j,
                motion.gas.stored_energy_change_j,
            ])?;
            let storage = sum([gas_change, motion.kinetic_energy_change_j])?;
            let dm = next.gas().mass_kg() - old.gas().mass_kg();
            let s0 = old
                .gas()
                .state()
                .specific_entropy_j_kg_k()
                .map_err(|e| mechanical(id, PistonError::Gas(e.into())))?;
            let s1 = next
                .gas()
                .state()
                .specific_entropy_j_kg_k()
                .map_err(|e| mechanical(id, PistonError::Gas(e.into())))?;
            let entropy = sum([next.gas().mass_kg() * (s1 - s0), dm * s0])?;
            // Replace the mechanical substage account with ONE complete interval.
            // The inherited two-port schema cannot certify an open species stream.
            next.piston.ports.accept(StepAccount {
                tick: old.piston.ports.next_tick(), duration_s: dt,
                storage_change_j: Some(storage), loss_j: Some(0.0), included_work_j: 0.0,
                entropy_change_j_k: entropy, tolerance_j: Some(controls.total_energy_tolerance_j),
                unavailable: (transfer.incoming_mass_kg > 0.0 || transfer.outgoing_mass_kg > 0.0)
                    .then_some("gas streams require stream ports; use the enclosing moving spatial-wall balance"),
            }, [PortTransfer::work(motion.external_force_work_j - motion.exterior_pressure_work_j
                    + motion.support_work_j, true),
                PortTransfer::heat(gas_heat.gas.supplied_heat_j, heat_temperature)]);
            mass_changes.push(dm);
            residuals.push(motion.total_energy_residual_j);
            reports.push(SolidPistonStep {
                id,
                motion,
                gas_energy_change_j: gas_change,
            });
            pistons.insert(id, next);
        }
        // Native stage endpoints coincide exactly, so their source-datum-cancelled
        // balances telescope across the whole interval, including actual kinetic
        // change and exterior work. Never subtract rounded formation totals here.
        let mass = sum(mass_changes)?;
        let energy = sum(residuals)?;
        if pistons.len() != ordered.len()
            || heat.thermal.gases.len() != ordered.len()
            || heat.flow.exchange.parcels.len() != ordered.len()
            || mass.abs() > controls.flow.mass_tolerance_kg
            || energy.abs() > controls.total_energy_tolerance_j
        {
            return Err(MovingSolidGasError::Invalid(
                "moving gas, finite wall and kinetic balances do not close",
            ));
        }
        checkpoint()?;
        Ok(MovingSolidGasCandidate {
            thermal,
            pistons,
            state,
            report: MovingSolidGasReport {
                heat,
                pistons: reports,
                total_mass_residual_kg: mass,
                total_energy_residual_j: energy,
            },
        })
    }
}
