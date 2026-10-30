//! Kirchhoff storage and frame materials on the original shared heat residual.
mod contact;

use super::*;
use crate::material_circuit::{
    MaterialCircuit, MaterialCircuitCandidate, MaterialCircuitSolve, MaterialCircuitStep,
};
use fs_material::conductor::ElectricalPropertyBasis;

/// Held electrical source ports and independent native circuit tolerances.
/// Use an empty port slice for a circuit powered only by its finite C/L stores.
#[derive(Clone, Copy)]
pub struct FrameElectricalInput<'p> {
    pub ports: &'p [f64],
    pub solve: MaterialCircuitSolve,
}

pub(super) struct Binding<'n, 'a> {
    contacts: Option<contact::Binding<'n, 'a>>,
    pub(super) circuit: MaterialCircuit,
    destinations: BTreeMap<usize, CellId>,
}

impl<'n, 'a> FrameThermoelasticSystem<'n, 'a> {
    /// Bind a circuit to the SAME material cells and conduction links as the
    /// complete thermoelastic assembly. `source` adds whole-cell resistors to
    /// the original network; its `cells` and `links` must borrow the original
    /// slices, not independent copies or a second set of material laws.
    ///
    /// Every circuit resistor must equal its freshly resolved initial source
    /// and geometry. Conductivity/resistivity bases are explicit per branch.
    /// The one prepared network replaces the thermal-only preparation, so each
    /// residual trial resolves mechanical temperature, R(T), and conduction at
    /// one endpoint. No additional material evaluator or thermal solve is used.
    ///
    /// The original charges and fluxes are retained. This is cold binding only:
    /// the frame sample clock must still be zero, and thermal cells must match
    /// the accepted frame state. Electrical contacts use
    /// [`Self::with_material_contact_circuit`]; material-derived C/L currently
    /// require their reciprocal thermal laws and refuse here. Explicitly
    /// authored fixed C/L retain their native storage accounting.
    /// Fixed geometry means no motional EMF, magnetic force, or changing wire
    /// length is implied. The mechanical/thermal split remains first order.
    #[allow(clippy::too_many_arguments)]
    pub fn with_material_circuit(
        mut self,
        thermal: &DepositionSimulation<'_>,
        source: &'n MaterialNetwork<'a>,
        circuit: MaterialCircuit,
        electrical_bases: &[(usize, ElectricalPropertyBasis)],
        gate: &CancelGate,
    ) -> Result<Self> {
        self.admit_electrical_binding(thermal, source, &circuit, gate)?;
        if !circuit.contacts().is_empty() {
            return Err(invalid("contact-bound frames require their original electrical contact network"));
        }
        let prepared = source.prepare_with_bases(thermal, electrical_bases,
            || gate.is_requested()).map_err(network)?;
        let mut resolved = prepared.resolve(thermal, || gate.is_requested()).map_err(network)?.resistors;
        resolved.sort_unstable_by_key(|resistor| resistor.branch);
        if resolved.as_slice() != circuit.resistors() {
            return Err(invalid("initial circuit must retain every original material resistor and its geometry"));
        }
        let destinations = source.resistors.iter().map(|r|
            (r.destination.branch, r.destination.cell)).collect();
        check(gate)?;
        self.prepared = prepared;
        self.electrical = Some(Binding { circuit, destinations, contacts: None });
        Ok(self)
    }

    /// Accepted native charge, flux, topology, and source material snapshots.
    /// No mutable circuit reference can advance it independently of frame heat.
    pub fn electrical(&self) -> Option<&MaterialCircuit> {
        self.electrical.as_ref().map(|binding| &binding.circuit)
    }

    /// Advance mechanics, circuit storage, and all finite material cells as one
    /// transaction. Every electrical trial starts from the same accepted charge
    /// and flux. Resistance and conduction use endpoint temperature; mechanical
    /// coefficients use the accepted start temperature over this interval.
    /// New source stiffness/expansion/viscosity act in the next interval.
    /// Any electrical, thermal, material, mechanical, or cancellation failure
    /// publishes none of the three histories. Sources are held over the sample.
    #[allow(clippy::too_many_arguments)]
    pub fn step_electrical(
        &mut self,
        thermal: &mut DepositionSimulation<'_>,
        external_modal_forces: &[f64],
        input: FrameElectricalInput<'_>,
        controls: FrameThermalSolve,
        energy_guess_j: &[f64],
        gate: &CancelGate,
    ) -> Result<FrameThermoelasticReport> {
        self.step_inner(thermal, external_modal_forces, Some(input), controls, energy_guess_j, gate)
    }
}

impl Binding<'_, '_> {
    pub(super) fn deposit(
        &self,
        candidate: &MaterialCircuitCandidate,
        dt: f64,
        powers: &mut BTreeMap<CellId, f64>,
        checkpoint: &mut dyn FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<()> {
        for branch in &candidate.step().resistors {
            checkpoint()?;
            let cell = self.destinations.get(&branch.material.branch)
                .ok_or_else(|| invalid("electrical loss has no original frame-network receiver"))?;
            let power = finite(branch.heat_j / dt)?;
            if branch.heat_j > 0.0 && power == 0.0 {
                return Err(invalid("frame electrical heat power underflow"));
            }
            let total = powers.entry(*cell).or_insert(0.0);
            *total = finite(*total + power)?;
        }
        if let Some(binding) = &self.contacts {
            crate::circuit_enthalpy::contacts::deposit(
                &candidate.step().contacts, &binding.routes, dt, powers, checkpoint,
            )?;
        } else if !candidate.step().contacts.is_empty() {
            return Err(invalid("frame contact heat has no original two-sided receivers"));
        }
        Ok(())
    }
}

// Add admitted native driver accounts, not synthetic heat rates. The shared
// caloric owner still independently compares integrated source powers with
// total physical dissipation and closes all supplied work and storage.
pub(super) struct Driver {
    mechanical: MechanicalCandidate,
    pub(super) electrical: Option<MaterialCircuitCandidate>,
}
impl Driver {
    pub(super) fn new(
        mechanical: MechanicalCandidate,
        electrical: Option<MaterialCircuitCandidate>,
    ) -> Result<Self> {
        let driver = Self { mechanical, electrical };
        let a = driver.account();
        finite(a.heat_j)?;
        finite(a.stored_change_j)?;
        finite(a.supplied_j)?;
        Ok(driver)
    }
}
impl ThermalDriver for Driver {
    type Model = Option<MaterialCircuit>;
    type Step = (CoupledModalFrame, Option<MaterialCircuitStep>);

    fn account(&self) -> DriverAccount<'_> {
        let mut account = self.mechanical.account();
        if let Some(candidate) = &self.electrical {
            let electrical = candidate.account();
            account.resistors = electrical.resistors;
            account.contacts = electrical.contacts;
            account.heat_j += electrical.heat_j;
            account.stored_change_j += electrical.stored_change_j;
            account.supplied_j += electrical.supplied_j;
        }
        account
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        let (_, mechanical) = self.mechanical.into_parts();
        match self.electrical {
            Some(candidate) => {
                let (circuit, electrical) = candidate.into_parts();
                (Some(circuit), (mechanical, Some(electrical)))
            }
            None => (None, (mechanical, None)),
        }
    }
}

// One constitutive/source-resolution seam for both bulk-only and contact-bearing
// circuits. Resolve the entire material network once per trial, then use the
// existing native source-continuity and circuit step. No nested thermal solve.
pub(super) fn resolve(
    binding: Option<&Binding<'_, '_>>,
    input: Option<FrameElectricalInput<'_>>,
    prepared: &PreparedMaterialNetwork<'_, '_>,
    trial: &DepositionSimulation<'_>,
    dt: f64,
    gate: &CancelGate,
) -> Result<(
    crate::circuit_enthalpy::CircuitThermalCoefficients<MaterialNetworkSample>,
    Option<MaterialCircuitCandidate>,
)> {
    check(gate)?;
    match (binding, input) {
        (Some(binding), Some(input)) => {
            let (mut bulk, contacts) = match &binding.contacts {
                Some(source) => {
                    let value = source.prepared.resolve(trial, || gate.is_requested())
                        .map_err(contact::source_error)?;
                    (value.bulk, value.contacts)
                }
                None => (prepared.resolve(trial, || gate.is_requested()).map_err(network)?, Vec::new()),
            };
            let candidate = crate::circuit_enthalpy::contacts::advance_circuit(
                &binding.circuit, core::mem::take(&mut bulk.resistors), contacts,
                input.ports, dt, input.solve, || gate.is_requested(),
            )?;
            Ok((bulk, Some(candidate)))
        }
        (None, None) => {
            let bulk = prepared.resolve(trial, || gate.is_requested()).map_err(network)?;
            if !bulk.resistors.is_empty() {
                return Err(invalid("frame electrical sources lack their original coupled owner"));
            }
            Ok((bulk, None))
        }
        _ => Err(invalid("frame electrical sources and held inputs disagree")),
    }
}
