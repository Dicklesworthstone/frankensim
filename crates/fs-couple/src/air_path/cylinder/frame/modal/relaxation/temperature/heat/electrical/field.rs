//! Existing thermodynamic field storage on the complete relaxation trajectory.
//! The original field owner supplies U=Psi+TS and reversible caloric exchange;
//! Joule and Maxwell loss keep their own independent dissipative accounts.
use super::*;
use crate::circuit_enthalpy::field_storage::{
    FieldThermalError, FieldThermalSource, PreparedFieldThermal,
};
use crate::circuit_enthalpy::contacts::material::ElectricalContactNetwork;
use std::cell::RefCell;

/// Actual material fields at both accepted endpoints; the unused network slot
/// is unit because the original bulk sample remains on RelaxationHeatReport.
pub type RelaxationFieldStorage = FieldThermalSource<()>;

/// Original field entropy and the signed first-order thermal-splitting defect.
/// These are not additional heat or dissipative entropy production.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RelaxationFieldEntropy {
    pub storage_change_j_k: f64,
    pub reversible_transfer_j_k: f64,
    /// Delta S_field + sum(Q_rev/T_endpoint), zero only in the reversible limit.
    pub temporal_defect_j_k: f64,
}

fn error(error: CircuitEnthalpyError<FieldThermalError>) -> RelaxationHeatError {
    match error {
        CircuitEnthalpyError::Source(e) =>
            CircuitEnthalpyError::Source(RelaxationHeatSourceError::Field(e)),
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

impl<'m, 'f, 'n, 'a, F> FrameRelaxationThermalSystem<'m, 'f, 'n, 'a, F> {
    /// Cold-bind source-defined C(T)/L(T) and bulk transport to the SAME finite
    /// material layers and complete Maxwell/fluid trajectory. Every field branch
    /// needs its original whole-volume FieldThermalBinding, even a constant law.
    /// Several independent fields can share a caloric cell: the original field
    /// owner sums their capacity contributions BEFORE admitting that cell once.
    /// Mechanical storage is temperature-independent in this rate-only model.
    ///
    /// All charge/flux, material strains and fluid momenta survive preparation.
    /// Existing step_electrical / step_heat_network advance the resulting owner.
    /// Source C/L are globally affine, geometry/density fixed, and caloric cells
    /// solid under the original field contract. No motional EMF, Maxwell force,
    /// expansion, hysteresis or implicit heat bath is inferred from this binding.
    #[allow(clippy::too_many_arguments)]
    pub fn bind_material_field_circuit(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        source: &'n MaterialNetwork<'a>,
        circuit: MaterialCircuit,
        electrical_bases: &[(usize, ElectricalPropertyBasis)],
        fields: Vec<FieldThermalBinding>,
        gate: &CancelGate,
    ) -> Result<()> {
        if fields.is_empty() { return Err(invalid("material field circuit requires complete thermal field bindings")); }
        self.bind_bulk_circuit(thermal, source, circuit, electrical_bases, fields, gate)
    }

    /// The same field/thermal transaction with original electrical interfaces
    /// and their two-sided Joule partitions. Conservative field exchange is
    /// deposited in its actual field material, never through a contact split.
    pub fn bind_material_contact_field_circuit(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        source: &'n ElectricalContactNetwork<'a>,
        circuit: MaterialCircuit,
        fields: Vec<FieldThermalBinding>,
        gate: &CancelGate,
    ) -> Result<()> {
        if fields.is_empty() { return Err(invalid("contact field circuit requires complete thermal field bindings")); }
        self.bind_contact_circuit(thermal, source, circuit, fields, gate)
    }
}

pub(super) fn admit(source: &MaterialNetwork<'_>, thermal: &DepositionSimulation<'_>,
    circuit: &MaterialCircuit, fields: &[FieldThermalBinding], gate: &CancelGate,
) -> Result<()> {
    let cancel = RefCell::new(|| gate.is_requested());
    PreparedFieldThermal::new(source, thermal, circuit.into(), fields, &cancel).map_err(error)?;
    check(gate)
}

// One accepted-coordinate/source snapshot per physical interval, not per
// Newton stencil. The helper does not borrow mutable electrical histories.
pub(in super::super) fn prepare<'b>(binding: Option<&'b Binding<'_, '_>>,
    source: &MaterialNetwork<'_>, thermal: &DepositionSimulation<'_>, gate: &CancelGate,
) -> Result<Option<PreparedFieldThermal<'b>>> {
    check(gate)?;
    match binding.filter(|binding| !binding.fields.is_empty()) {
        Some(binding) => PreparedFieldThermal::new(source, thermal, (&binding.circuit).into(),
            &binding.fields, &RefCell::new(|| gate.is_requested())).map(Some).map_err(error),
        None => Ok(None),
    }
}

pub(in super::super) struct Trial {
    pub storage: RelaxationFieldStorage,
    // Directly captured before reversible heat is added. A full heat transport
    // source reclassifies its own radiation separately at final acceptance.
    independent_entropy_j_k: f64,
    has_independent_sources: bool,
}

#[allow(clippy::too_many_arguments)]
pub(in super::super) fn advance(
    binding: Option<&Binding<'_, '_>>,
    input: Option<RelaxationElectricalInput<'_>>,
    fields: Option<&PreparedFieldThermal<'_>>,
    resistors: Vec<MaterialResistor>, contacts: Vec<SeriesElectricalContact>,
    trial: &DepositionSimulation<'_>, sources_w: &mut BTreeMap<CellId, f64>,
    dt: f64, tolerance_j: f64, gate: &CancelGate,
) -> Result<(Option<MaterialCircuitCandidate>, Option<Trial>)> {
    check(gate)?;
    if fields.is_none() {
        if binding.is_some_and(|binding| !binding.fields.is_empty()) {
            return Err(invalid("source field history cannot be frozen by omitting its thermal owner"));
        }
        return Ok((super::advance(binding, input, resistors, contacts, dt, gate)?, None));
    }
    let (Some(fields), Some(binding), Some(input)) = (fields, binding, input) else {
        return Err(invalid("complete source field and electrical input required"));
    };
    let cancel = RefCell::new(|| gate.is_requested());
    let (storage, samples) = fields.refresh(trial, &cancel).map_err(error)?;
    let rebound = binding.circuit.rebind_storage_candidate(storage, input.ports,
        || gate.is_requested())?;
    let candidate = crate::circuit_enthalpy::contacts::advance_circuit(
        rebound.circuit(), resistors, contacts, input.ports, dt, input.solve,
        || gate.is_requested(),
    )?;
    let mut independent_entropy_j_k = 0.0;
    let mut has_independent_sources = false;
    for (&id, cell) in trial.cells() {
        check(gate)?;
        let temperature = trial.temperature_k(id)?;
        for power in [cell.boundary.source_w, sources_w.get(&id).copied().unwrap_or(0.0)] {
            if power != 0.0 { has_independent_sources = true; }
            independent_entropy_j_k = finite(independent_entropy_j_k
                + super::super::transport::source_entropy(power, dt, temperature)?)?;
        }
    }
    let storage = fields.account(samples, candidate.circuit().into(),
        rebound.parameter_work_j(), dt, (), sources_w, tolerance_j, &cancel).map_err(error)?;
    check(gate)?;
    Ok((Some(candidate), Some(Trial {storage, independent_entropy_j_k, has_independent_sources})))
}

pub(in super::super) fn finish(trial: Option<Trial>, has_transport: bool,
    radiation_entropy: f64, entropy: &mut FrameHeatEntropy, gate: &CancelGate,
) -> Result<(Option<RelaxationFieldStorage>, Option<RelaxationFieldEntropy>)> {
    let Some(trial) = trial else { return Ok((None, None)); };
    let mut storage_change = 0.0;
    let mut reversible = 0.0;
    for sample in &trial.storage.fields {
        check(gate)?;
        storage_change = finite(storage_change + sample.state[1].entropy_j_k - sample.state[0].entropy_j_k)?;
        reversible = finite(reversible + sample.reversible_heat_j / sample.temperature_k[1])?;
    }
    if !has_transport {
        entropy.unclassified_source_at_cell_j_k = trial.independent_entropy_j_k;
        entropy.has_unclassified_sources = trial.has_independent_sources;
    }
    entropy.accounting_residual_j_k = [entropy.thermal_storage_change_j_k,
        -entropy.reservoir_entropy_input_j_k, -entropy.conduction_generation_j_k,
        -entropy.material_loss_generation_j_k, -entropy.unclassified_source_at_cell_j_k,
        -radiation_entropy, -reversible, -entropy.temporal_defect_j_k,
        -entropy.algebraic_residual_j_k,
    ].into_iter().try_fold(0.0, |sum, term| finite(sum + term))?;
    let field_entropy = RelaxationFieldEntropy { storage_change_j_k: storage_change,
        reversible_transfer_j_k: reversible, temporal_defect_j_k: finite(storage_change + reversible)? };
    check(gate)?;
    Ok((Some(trial.storage), Some(field_entropy)))
}
