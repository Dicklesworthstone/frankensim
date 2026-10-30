//! Original field thermodynamics on the shared mechanical/caloric transaction.
use super::*;
use crate::circuit_enthalpy::{
    CircuitThermalCoefficients,
    contacts::material::ElectricalContactNetwork,
    field_storage::{FieldThermalError, FieldThermalSample},
};
use crate::material_circuit::{MaterialCircuitError, contacts::SeriesElectricalContact};
use std::cell::RefCell;

/// Source-resolved field storage at both actual thermal/electrical endpoints.
/// The native electrical report instead uses isothermal free energy after the
/// coefficient rebind. These physical U/S quantities belong to the material.
#[derive(Clone, Debug, PartialEq)]
pub struct FrameFieldStorageReport {
    pub fields: Vec<FieldThermalSample>,
    pub internal_energy_j: [f64; 2],
    /// Signed reversible transfer into caloric storage; NEVER Joule loss.
    pub reversible_heat_j: f64,
}
impl FrameFieldStorageReport {
    pub(in super::super) fn entropy_changes(&self) -> Result<(f64, f64)> {
        Ok((
            sum(self.fields.iter().map(|f| f.state[1].entropy_j_k - f.state[0].entropy_j_k))?,
            sum(self.fields.iter().map(|f| f.reversible_heat_j / f.temperature_k[1]))?,
        ))
    }
}

// Preserve all source/numerical error variants through the frame's existing
// error vocabulary; a field source refusal is not an invented generic failure.
fn error(e: CircuitEnthalpyError<FieldThermalError>) -> FrameThermalError {
    match e {
        CircuitEnthalpyError::Source(FieldThermalError::Network(e)) => network(e),
        CircuitEnthalpyError::Source(FieldThermalError::Field(e)) =>
            MaterialCircuitError::Storage(e).into(),
        CircuitEnthalpyError::Source(FieldThermalError::Invalid(s)) => invalid(s),
        CircuitEnthalpyError::Invalid(s) => CircuitEnthalpyError::Invalid(s),
        CircuitEnthalpyError::Circuit(e) => CircuitEnthalpyError::Circuit(e),
        CircuitEnthalpyError::Thermal(e) => e.into(),
        CircuitEnthalpyError::Newton(e) => CircuitEnthalpyError::Newton(e),
        CircuitEnthalpyError::NotConverged(e) => CircuitEnthalpyError::NotConverged(e),
        CircuitEnthalpyError::Interrupted(e) => CircuitEnthalpyError::Interrupted(e),
        CircuitEnthalpyError::CircuitBudget => CircuitEnthalpyError::CircuitBudget,
        CircuitEnthalpyError::Cancelled => CircuitEnthalpyError::Cancelled,
    }
}

impl<'n, 'a> FrameThermoelasticSystem<'n, 'a> {
    /// Bind original source-dependent C(T)/L(T) and bulk resistance to the SAME
    /// finite material cells and native thermoelastic assembly. Every storage
    /// branch needs a whole-volume FieldThermalBinding, including constant laws.
    /// The field and frame may occupy the same cell: their capacity contributions
    /// are admitted TOGETHER, with caloric storage counted only once.
    ///
    /// Use step_electrical normally after binding. Field coefficient work and
    /// entropy exchange come from the existing field-storage owner; actual
    /// charge/flux and source laws remain on the original circuit. Changes in
    /// field material data can therefore alter heat, elasticity and later motion.
    /// Geometry is fixed: no motional EMF, Maxwell traction, electrostriction or
    /// field-induced geometry change is inferred. The existing first-order
    /// mechanical/thermal and field/circuit splits still require refinement.
    #[allow(clippy::too_many_arguments)]
    pub fn with_material_field_circuit(
        self,
        thermal: &DepositionSimulation<'_>,
        source: &'n MaterialNetwork<'a>,
        circuit: MaterialCircuit,
        electrical_bases: &[(usize, ElectricalPropertyBasis)],
        fields: Vec<FieldThermalBinding>,
        gate: &CancelGate,
    ) -> Result<Self> {
        check(gate)?;
        if fields.is_empty() {
            return Err(invalid("material field circuit requires complete thermal field bindings"));
        }
        self.bind_material_circuit(thermal, source, circuit, electrical_bases, fields, gate)
    }

    /// Field-storage version of with_material_contact_circuit. The same pinned
    /// interface network resolves both sides at the joint thermal endpoint;
    /// its existing partition deposits each contact loss once. Reversible field
    /// exchange, reversible elastic exchange and Joule/dashpot heat remain
    /// distinct. No additional thermal step or contact law is introduced.
    pub fn with_material_contact_field_circuit(
        self,
        thermal: &DepositionSimulation<'_>,
        source: &'n ElectricalContactNetwork<'a>,
        circuit: MaterialCircuit,
        fields: Vec<FieldThermalBinding>,
        gate: &CancelGate,
    ) -> Result<Self> {
        check(gate)?;
        if fields.is_empty() {
            return Err(invalid("material contact field circuit requires thermal field bindings"));
        }
        self.bind_material_contact_circuit(thermal, source, circuit, fields, gate)
    }
}

// Independent positivity of (caloric + elastic) and each (caloric + field)
 // would spend the same caloric capacity repeatedly. Sum every co-located field
 // contribution first, then intersect that total with the original conservative
 // thermoelastic floor BEFORE any solve. Source-frame layers have disjoint cells.
pub(super) fn admit(
    owner: &FrameThermoelasticSystem<'_, '_>,
    source: &MaterialNetwork<'_>,
    thermal: &DepositionSimulation<'_>,
    circuit: &MaterialCircuit,
    fields: &[FieldThermalBinding],
    gate: &CancelGate,
) -> Result<()> {
    let cancel = RefCell::new(|| gate.is_requested());
    PreparedFieldThermal::new(source, thermal, circuit.into(), fields, &cancel).map_err(error)?;
    let mut mechanical_floors = BTreeMap::new();
    for member in &owner.members {
        for (receiver, &floor) in member.region.receivers.iter().zip(&member.capacity_floor_j_k) {
            check(gate)?;
            mechanical_floors.insert(receiver.binding.cell, floor);
        }
    }
    let mut field_floors = BTreeMap::<CellId, f64>::new();
    for field in fields {
        check(gate)?;
        let field_floor = field.law.minimum_heat_capacity_j_k(field.maximum_abs_coordinate)
            .map_err(|e| FrameThermalError::from(MaterialCircuitError::Storage(e)))?;
        let total = field_floors.entry(field.cell).or_insert(0.0);
        *total = finite(*total + field_floor)?;
    }
    for (cell, field_floor) in field_floors {
        check(gate)?;
        if let Some(&mechanical_floor) = mechanical_floors.get(&cell)
            && finite(mechanical_floor + field_floor)? <= 0.0
        {
            return Err(invalid(
                "combined caloric, elastic and all co-located field heat capacities must stay positive",
            ));
        }
    }
    check(gate)
}

// Prepare once per physical interval, never once per Newton stencil. The
// original helper retains only accepted endpoint coordinates and source laws.
pub(in super::super) fn prepare_fields<'b>(
    binding: Option<&'b Binding<'_, '_>>,
    source: &MaterialNetwork<'_>,
    thermal: &DepositionSimulation<'_>,
    gate: &CancelGate,
) -> Result<Option<PreparedFieldThermal<'b>>> {
    check(gate)?;
    match binding.filter(|b| !b.fields.is_empty()) {
        Some(binding) => {
            let cancel = RefCell::new(|| gate.is_requested());
            PreparedFieldThermal::new(source, thermal, (&binding.circuit).into(),
                &binding.fields, &cancel).map(Some).map_err(error)
        }
        None => Ok(None),
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn advance(
    binding: &Binding<'_, '_>,
    input: FrameElectricalInput<'_>,
    fields: Option<&PreparedFieldThermal<'_>>,
    mut bulk: CircuitThermalCoefficients<MaterialNetworkSample>,
    contacts: Vec<SeriesElectricalContact>,
    trial: &DepositionSimulation<'_>,
    dt: f64,
    tolerance_j: f64,
    gate: &CancelGate,
) -> Result<ResolvedElectrical> {
    check(gate)?;
    if fields.is_some() != !binding.fields.is_empty() {
        return Err(invalid("material field history cannot be silently frozen"));
    }
    let cancel = RefCell::new(|| gate.is_requested());
    let (rebound, samples) = match fields {
        Some(fields) => {
            let (storage, samples) = fields.refresh(trial, &cancel).map_err(error)?;
            let rebound = binding.circuit.rebind_storage_candidate(storage, input.ports,
                || gate.is_requested())?;
            (Some(rebound), samples)
        }
        None => (None, Vec::new()),
    };
    let circuit = rebound.as_ref().map_or(&binding.circuit, |r| r.circuit());
    let candidate = crate::circuit_enthalpy::contacts::advance_circuit(
        circuit, core::mem::take(&mut bulk.resistors), contacts, input.ports, dt,
        input.solve, || gate.is_requested(),
    )?;
    let mut reversible_sources = BTreeMap::new();
    let field_storage = match (fields, rebound) {
        (Some(fields), Some(rebound)) => {
            let source = fields.account(samples, candidate.circuit().into(),
                rebound.parameter_work_j(), dt, bulk.source, &mut reversible_sources,
                tolerance_j, &cancel).map_err(error)?;
            bulk.source = source.network;
            Some(FrameFieldStorageReport { fields: source.fields,
                internal_energy_j: source.field_internal_energy_j,
                reversible_heat_j: source.reversible_heat_j })
        }
        (None, None) => None,
        _ => return Err(invalid("complete field/circuit coefficient candidate is absent")),
    };
    check(gate)?;
    Ok(ResolvedElectrical { bulk, circuit: Some(candidate), field_storage, reversible_sources })
}

#[cfg(test)]
mod tests;
