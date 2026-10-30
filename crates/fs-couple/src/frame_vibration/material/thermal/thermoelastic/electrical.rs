//! Kirchhoff storage and frame materials on the original shared heat residual.
mod contact;
mod field;
mod geometry;
pub use field::FrameFieldStorageReport;
pub use geometry::{FrameConductorGeometry, FrameConductorPath};
pub(super) use field::prepare_fields;
use crate::circuit_enthalpy::field_storage::{FieldThermalBinding, PreparedFieldThermal};

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
    fields: Vec<FieldThermalBinding>,
    pub(super) conductor_paths: Vec<geometry::BoundConductorPath>,
    pub(super) conductor_geometry: Vec<FrameConductorGeometry>,
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
    /// [`Self::with_material_contact_circuit`]; material-derived C/L use
    /// [`Self::with_material_field_circuit`] and their reciprocal thermal laws. Explicitly
    /// authored fixed C/L retain their native storage accounting.
    /// Fixed geometry means no motional EMF, magnetic force, or changing wire
    /// length is implied. The mechanical/thermal split remains first order.
    #[allow(clippy::too_many_arguments)]
    pub fn with_material_circuit(
        self,
        thermal: &DepositionSimulation<'_>,
        source: &'n MaterialNetwork<'a>,
        circuit: MaterialCircuit,
        electrical_bases: &[(usize, ElectricalPropertyBasis)],
        gate: &CancelGate,
    ) -> Result<Self> {
        self.bind_material_circuit(thermal, source, circuit, electrical_bases, Vec::new(), gate)
    }

    #[allow(clippy::too_many_arguments)]
    fn bind_material_circuit(
        mut self,
        thermal: &DepositionSimulation<'_>,
        source: &'n MaterialNetwork<'a>,
        circuit: MaterialCircuit,
        electrical_bases: &[(usize, ElectricalPropertyBasis)],
        fields: Vec<FieldThermalBinding>,
        gate: &CancelGate,
    ) -> Result<Self> {
        self.admit_electrical_binding(thermal, source, &circuit, &fields, gate)?;
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
        self.electrical = Some(Binding { circuit, destinations, contacts: None, fields,
            conductor_paths: Vec::new(), conductor_geometry: Vec::new() });
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
        let previous = self.refresh_conductor_geometry(gate)?;
        let result = self.step_inner(
            thermal, external_modal_forces, Some(input), None, controls, energy_guess_j, gate,
        );
        if result.is_err() { self.restore_conductor_geometry(previous); }
        result
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

// Original material resolution and native circuit work, with optional field
// thermodynamics. A separate reversible-source map keeps physical exchange out
// of authored-source entropy; the caller inserts it into the SAME residual.
pub(super) struct ResolvedElectrical {
    pub bulk: crate::circuit_enthalpy::CircuitThermalCoefficients<MaterialNetworkSample>,
    pub circuit: Option<MaterialCircuitCandidate>,
    pub field_storage: Option<FrameFieldStorageReport>,
    pub reversible_sources: BTreeMap<CellId, f64>,
}

#[allow(clippy::too_many_arguments)]
pub(super) fn resolve(
    binding: Option<&Binding<'_, '_>>,
    input: Option<FrameElectricalInput<'_>>,
    prepared: &PreparedMaterialNetwork<'_, '_>,
    fields: Option<&PreparedFieldThermal<'_>>,
    trial: &DepositionSimulation<'_>,
    dt: f64,
    tolerance_j: f64,
    gate: &CancelGate,
) -> Result<ResolvedElectrical> {
    check(gate)?;
    let coefficients = match binding {
        Some(binding) => match &binding.contacts {
            Some(source) => source.prepared.resolve(trial, || gate.is_requested())
                .map_err(contact::source_error)?,
            None => crate::circuit_enthalpy::contacts::ContactThermalCoefficients {
                bulk: prepared.resolve(trial, || gate.is_requested()).map_err(network)?,
                contacts: Vec::new(),
            },
        },
        None => crate::circuit_enthalpy::contacts::ContactThermalCoefficients {
            bulk: prepared.resolve(trial, || gate.is_requested()).map_err(network)?,
            contacts: Vec::new(),
        },
    };
    advance_resolved(binding, input, fields, coefficients, trial, dt, tolerance_j, gate)
}

// Source resolution is separable from native electrical evolution so a complete
// heat network can supply the SAME thermal trial without repeating bulk work.
#[allow(clippy::too_many_arguments)]
pub(super) fn advance_resolved(
    binding: Option<&Binding<'_, '_>>,
    input: Option<FrameElectricalInput<'_>>,
    fields: Option<&PreparedFieldThermal<'_>>,
    coefficients: crate::circuit_enthalpy::contacts::ContactThermalCoefficients<MaterialNetworkSample>,
    trial: &DepositionSimulation<'_>,
    dt: f64,
    tolerance_j: f64,
    gate: &CancelGate,
) -> Result<ResolvedElectrical> {
    check(gate)?;
    let crate::circuit_enthalpy::contacts::ContactThermalCoefficients { mut bulk, contacts } = coefficients;
    if let Some(binding) = binding {
        geometry::apply(&binding.conductor_geometry, &mut bulk.resistors, gate)?;
    }
    match (binding, input) {
        (Some(binding), Some(input)) =>
            field::advance(binding, input, fields, bulk, contacts, trial, dt, tolerance_j, gate),
        (None, None) => {
            if !bulk.resistors.is_empty() || !contacts.is_empty() || fields.is_some() {
                return Err(invalid("frame electrical sources lack their original coupled owner"));
            }
            Ok(ResolvedElectrical { bulk, circuit: None, field_storage: None,
                reversible_sources: BTreeMap::new() })
        }
        _ => Err(invalid("frame electrical sources and held inputs disagree")),
    }
}

// Full heat transport may replace the coefficient evaluator, never the accepted
// electrical specimen or its physical heat partitions. Keep original maps on
// the native owner and compare before any speculative mechanical step.
pub(super) fn admit_heat_sources<S>(
    binding: Option<&Binding<'_, '_>>,
    coefficients: &mut crate::circuit_enthalpy::contacts::ContactThermalCoefficients<S>,
    destinations: &[crate::circuit_enthalpy::ResistorCell],
    routes: &[crate::circuit_enthalpy::contacts::ContactHeatCells],
    thermal: &DepositionSimulation<'_>,
    gate: &CancelGate,
) -> Result<()> {
    check(gate)?;
    match binding {
        Some(binding) => {
            geometry::admit_heat_resistors(
                binding, &mut coefficients.bulk.resistors, gate,
            )?;
            crate::circuit_enthalpy::contacts::admit_sources(
                &binding.circuit, &coefficients.contacts, || gate.is_requested(),
            )?;
            let map: BTreeMap<_, _> = destinations.iter().map(|r| (r.branch, r.cell)).collect();
            if map.len() != destinations.len() || map != binding.destinations {
                return Err(invalid("complete heat network changed original resistor receivers"));
            }
            let routes = crate::circuit_enthalpy::contacts::admit_routes(
                thermal, &binding.circuit, routes, || {
                    if gate.is_requested() { Err(DepositionError::Cancelled) } else { Ok(()) }
                },
            )?;
            let same = match &binding.contacts {
                Some(original) => routes == original.routes,
                None => routes.is_empty(),
            };
            if !same {
                return Err(invalid("complete heat network changed original contact heat partitions"));
            }
        }
        None => {
            if !coefficients.bulk.resistors.is_empty() || !coefficients.contacts.is_empty()
                || !destinations.is_empty() || !routes.is_empty() {
                return Err(invalid("unpowered frame cannot discard heat-network electrical work"));
            }
        }
    }
    check(gate)
}
