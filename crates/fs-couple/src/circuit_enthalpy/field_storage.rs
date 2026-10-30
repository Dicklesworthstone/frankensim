//! Source-bound electric and magnetic storage on the original thermal residual.
//!
//! The same material temperature determines C(T) or L(T), the native circuit
//! effort, and the field's thermodynamic internal energy. At fixed charge q or
//! flux linkage lambda, psi = z^2/(2 K(T)), s = -partial_T psi and U = psi + T s.
//! A frozen-coordinate coefficient update followed by the original midpoint
//! circuit step supplies Q_rev = -delta_psi_parameter - delta(T s) to the
//! material's caloric store. This signed reversible transfer is NEVER Joule loss.
//!
//! This is a first-order endpoint-temperature coupling with fixed geometry,
//! globally affine source coefficients and fixed-density solid caloric charts.
//! It does not model hysteresis, dielectric relaxation, changing field geometry,
//! phase transitions or spatially resolved magnetic-field reassembly.

/// Complete thermal-contact and radiative transport on the same field/circuit step.
pub mod heat;

use super::{
    CircuitEnthalpyError, CircuitEnthalpyInput, CircuitEnthalpyReport, CircuitEnthalpySolve,
    DrivenCoefficients, coupled_driver_candidate, finite,
    material::{MaterialNetwork, MaterialNetworkError, MaterialNetworkSample},
};
use crate::material_circuit::storage::{MaterialStorage, MaterialStorageElement};
use crate::material_circuit::{MaterialCircuit, MaterialCircuitStep};
use fs_material::{
    conductor::ElectricalPropertyBasis,
    field_storage::{
        FieldStorageError,
        thermal::{
            FieldThermalState, ResolvedThermalFieldStorage, ThermalFieldElement,
            ThermalFieldStorage,
        },
    },
    state_point::ResolvedMaterialStatePoint,
};
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use std::{cell::RefCell, collections::BTreeMap};

/// One original field branch and its whole, thermally uniform material volume.
/// Every capacitor and inductor requires a binding, including constant laws.
pub struct FieldThermalBinding {
    pub branch: usize,
    pub cell: CellId,
    /// The original material, context, geometry and continuously admitted T range.
    pub law: ThermalFieldStorage,
    /// Bound on |q| [C] or |lambda| [Wb] for this physical step. Both accepted
    /// endpoints must satisfy it. Together with the entire temperature range,
    /// it admits positive total caloric + field heat capacity before solving.
    pub maximum_abs_coordinate: f64,
}

/// Source and model refusals remain distinct from native solver failures.
#[derive(Debug)]
pub enum FieldThermalError {
    Network(MaterialNetworkError),
    Field(FieldStorageError),
    Invalid(&'static str),
}
impl core::fmt::Display for FieldThermalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material field storage/enthalpy: {self:?}")
    }
}
impl std::error::Error for FieldThermalError {}
type Error = CircuitEnthalpyError<FieldThermalError>;

/// Both endpoints of one original material field, including reversible heat.
#[derive(Clone, Debug, PartialEq)]
pub struct FieldThermalSample {
    pub branch: usize,
    pub cell: CellId,
    pub temperature_k: [f64; 2],
    /// Charge or flux linkage, in the original branch's coordinate convention.
    pub coordinate: [f64; 2],
    pub state: [FieldThermalState; 2],
    /// Actual source-resolved field element at the accepted thermal endpoint.
    pub material: ThermalFieldElement,
    /// Frozen-initial-coordinate free-energy change, NOT external supply [J].
    pub parameter_work_j: f64,
    /// Signed transfer INTO this same material's caloric store [J].
    pub reversible_heat_j: f64,
}

/// Actual endpoint material data and the total physical field-energy account.
#[derive(Clone, Debug, PartialEq)]
pub struct FieldThermalSource<S = MaterialNetworkSample> {
    pub network: S,
    pub fields: Vec<FieldThermalSample>,
    pub field_internal_energy_j: [f64; 2],
    pub reversible_heat_j: f64,
}

/// Independent physical account, using U rather than the isothermal psi alone.
pub struct FieldThermalReport<D = MaterialCircuitStep, S = MaterialNetworkSample> {
    /// Original native circuit and thermal reports. The native circuit's
    /// stored_before_j is psi at the proposed endpoint temperature and the
    /// INITIAL charge/flux. non_electrical_heat_j includes reversible heat.
    pub coupled: CircuitEnthalpyReport<FieldThermalSource<S>, D>,
    /// Authored thermal boundaries/sources, excluding Joule AND reversible heat.
    pub external_heat_j: f64,
    /// delta(U_field) + delta(U_caloric) - electrical work - external heat [J].
    pub total_energy_residual_j: f64,
}

/// One publishable electrical/thermal endpoint; trials mutate no accepted state.
pub struct FieldThermalCandidate<'m, C = MaterialCircuit, D = MaterialCircuitStep, S = MaterialNetworkSample> {
    thermal: DepositionSimulation<'m>,
    circuit: C,
    report: FieldThermalReport<D, S>,
}
impl<'m, C, D, S> FieldThermalCandidate<'m, C, D, S> {
    pub(super) fn new(
        thermal: DepositionSimulation<'m>,
        circuit: C,
        report: FieldThermalReport<D, S>,
    ) -> Self {
        Self {
            thermal,
            circuit,
            report,
        }
    }
    pub const fn thermal(&self) -> &DepositionSimulation<'m> {
        &self.thermal
    }
    pub const fn circuit(&self) -> &C {
        &self.circuit
    }
    pub const fn report(&self) -> &FieldThermalReport<D, S> {
        &self.report
    }
    pub fn into_parts(self) -> (DepositionSimulation<'m>, C, FieldThermalReport<D, S>) {
        (self.thermal, self.circuit, self.report)
    }
}

fn invalid(message: &'static str) -> Error {
    CircuitEnthalpyError::Source(FieldThermalError::Invalid(message))
}
fn field(error: FieldStorageError) -> Error {
    CircuitEnthalpyError::Source(FieldThermalError::Field(error))
}
fn network(error: MaterialNetworkError) -> Error {
    match error {
        MaterialNetworkError::Thermal(DepositionError::Cancelled) => {
            CircuitEnthalpyError::Cancelled
        }
        other => CircuitEnthalpyError::Source(FieldThermalError::Network(other)),
    }
}
fn check(cancel: &RefCell<impl FnMut() -> bool>) -> Result<(), Error> {
    if (*cancel.borrow_mut())() {
        Err(CircuitEnthalpyError::Cancelled)
    } else {
        Ok(())
    }
}
fn close_geometry(a: f64, b: f64) -> bool {
    a.is_finite()
        && b.is_finite()
        && a > 0.0
        && b > 0.0
        && (a - b).abs() <= 64.0 * f64::EPSILON * a.max(b)
}

// A pinned source refresh changes the selection-policy receipt even at the same
// physical state. Compare exact actual sources/coordinates/values, preserving
// both full receipts rather than fabricating policy equality.
fn same_material(a: &ResolvedMaterialStatePoint, b: &ResolvedMaterialStatePoint) -> bool {
    a.card_identity() == b.card_identity()
        && a.material() == b.material()
        && a.axis_quantities() == b.axis_quantities()
        && a.query_point().len() == b.query_point().len()
        && a.query_point()
            .iter()
            .zip(b.query_point())
            .all(|((n, v), (m, w))| n == m && v.to_bits() == w.to_bits())
        && a.properties().len() == b.properties().len()
        && a.properties().iter().zip(b.properties()).all(|(p, q)| {
            p.requirement() == q.requirement()
                && p.source_axes() == q.source_axes()
                && p.answer().receipt.selected == q.answer().receipt.selected
                && p.value_si().to_bits() == q.value_si().to_bits()
        })
}
fn element(value: &ThermalFieldElement) -> MaterialStorageElement {
    match value {
        ThermalFieldElement::Capacitor(c) => MaterialStorageElement::Capacitor(c.clone()),
        ThermalFieldElement::Inductor(l) => MaterialStorageElement::Inductor(l.clone()),
    }
}
fn same_element(a: &MaterialStorageElement, b: &ThermalFieldElement) -> bool {
    match (a, b) {
        (MaterialStorageElement::Capacitor(a), ThermalFieldElement::Capacitor(b)) => {
            a.area_m2().to_bits() == b.area_m2().to_bits()
                && a.gap_m().to_bits() == b.gap_m().to_bits()
                && a.capacitance_f().to_bits() == b.capacitance_f().to_bits()
                && same_material(a.material(), b.material())
        }
        (MaterialStorageElement::Inductor(a), ThermalFieldElement::Inductor(b)) => {
            a.area_m2().to_bits() == b.area_m2().to_bits()
                && a.path_length_m().to_bits() == b.path_length_m().to_bits()
                && a.turns() == b.turns()
                && a.inductance_h().to_bits() == b.inductance_h().to_bits()
                && same_material(a.material(), b.material())
        }
        _ => false,
    }
}
fn volume(value: &ThermalFieldElement) -> f64 {
    match value {
        ThermalFieldElement::Capacitor(c) => c.area_m2() * c.gap_m(),
        ThermalFieldElement::Inductor(l) => l.area_m2() * l.path_length_m(),
    }
}
// Borrow the original native storage; this view owns no electrical history or
// stepping law. Lumped and spatial circuits use the same field thermodynamics.
pub(crate) struct FieldCircuitView<'a> {
    pub graph: &'a CircuitGraph,
    pub storage: &'a [MaterialStorage],
    pub charges_c: &'a [f64],
    pub fluxes_wb: &'a [f64],
}
impl<'a> From<&'a MaterialCircuit> for FieldCircuitView<'a> {
    fn from(circuit: &'a MaterialCircuit) -> Self {
        Self {
            graph: circuit.graph(),
            storage: circuit.storage(),
            charges_c: circuit.charges_c(),
            fluxes_wb: circuit.fluxes_wb(),
        }
    }
}
fn coordinates(circuit: &FieldCircuitView<'_>) -> Result<BTreeMap<usize, f64>, Error> {
    let mut charges = circuit.charges_c.iter();
    let mut fluxes = circuit.fluxes_wb.iter();
    let mut result = BTreeMap::new();
    for (branch, (_, _, kind)) in circuit.graph.branches.iter().enumerate() {
        let value = match kind {
            Branch::Capacitor { .. } => charges.next(),
            Branch::Inductor { .. } => fluxes.next(),
            _ => continue,
        }
        .copied()
        .ok_or_else(|| invalid("native circuit storage coordinate is absent"))?;
        if !value.is_finite() {
            return Err(invalid("nonfinite native field coordinate"));
        }
        result.insert(branch, value);
    }
    if charges.next().is_some() || fluxes.next().is_some() {
        return Err(invalid(
            "native field coordinates do not match their original branches",
        ));
    }
    Ok(result)
}

struct InitialField {
    temperature_k: f64,
    coordinate: f64,
    state: FieldThermalState,
}

fn admit(
    material: &MaterialNetwork<'_>,
    thermal: &DepositionSimulation<'_>,
    circuit: &FieldCircuitView<'_>,
    bindings: &[FieldThermalBinding],
    cancel: &RefCell<impl FnMut() -> bool>,
) -> Result<Vec<InitialField>, Error> {
    let initial = coordinates(circuit)?;
    if bindings.is_empty()
        || bindings.len() != initial.len()
        || bindings.len() != circuit.storage.len()
        || bindings.len() > crate::material_circuit::MAX_MATERIAL_CIRCUIT_BRANCHES
    {
        return Err(invalid(
            "complete bounded thermal binding for all native field storage required",
        ));
    }
    let mut branches = BTreeMap::new();
    let mut cells = BTreeMap::new();
    let mut result = Vec::with_capacity(bindings.len());
    for binding in bindings {
        check(cancel)?;
        let bound = binding.maximum_abs_coordinate;
        if !bound.is_finite()
            || bound < 0.0
            || branches.insert(binding.branch, ()).is_some()
            || cells.insert(binding.cell, ()).is_some()
        {
            return Err(invalid(
                "unique field branches and whole material cells with finite coordinate bounds required",
            ));
        }
        let source = circuit
            .storage
            .iter()
            .find(|s| s.branch == binding.branch)
            .ok_or_else(|| invalid("field branch lacks original source-bound storage"))?;
        let cell = thermal
            .cells()
            .get(&binding.cell)
            .ok_or_else(|| invalid("field thermal cell is inactive"))?;
        let material_cell = material
            .cells
            .iter()
            .find(|s| s.cell == binding.cell)
            .ok_or_else(|| invalid("field thermal cell lacks its original caloric binding"))?;
        let phase = material_cell
            .curve
            .state_at_specific_enthalpy(cell.specific_h_j_kg)
            .map_err(|e| network(MaterialNetworkError::Phase(e)))?;
        let point = material_cell
            .point
            .at(phase)
            .map_err(|e| network(MaterialNetworkError::Material(e)))?;
        let current = binding
            .law
            .at_temperature(phase.temperature_k())
            .map_err(field)?;
        let snapshot = match current.element() {
            ThermalFieldElement::Capacitor(c) => c.material(),
            ThermalFieldElement::Inductor(l) => l.material(),
        };
        if material_cell.point.material.content_hash() != binding.law.card().content_hash()
            || material_cell.point.temperature_axis != binding.law.temperature_axis()
            || point != snapshot.query()
            || !same_element(&source.element, current.element())
            || !close_geometry(volume(current.element()), cell.geometry.volume_m3)
        {
            return Err(invalid(
                "field, caloric and transport owners must use the same card, actual context and whole volume",
            ));
        }
        let range = binding.law.temperature_range_k();
        let knots = material_cell.curve.knots();
        if range[0] < knots[0].temperature_k || range[1] > knots[knots.len() - 1].temperature_k {
            return Err(invalid(
                "the entire field temperature range requires caloric support",
            ));
        }
        let mut minimum_caloric_capacity = f64::INFINITY;
        for knot in knots {
            check(cancel)?;
            if knot.liquid_mass_fraction != 0.0
                || !close_geometry(
                    cell.mass_kg / knot.bulk_density_kg_m3,
                    cell.geometry.volume_m3,
                )
            {
                return Err(invalid(
                    "field thermal geometry requires a fully solid constant-density caloric chart",
                ));
            }
        }
        for pair in knots.windows(2) {
            check(cancel)?;
            let dt = pair[1].temperature_k - pair[0].temperature_k;
            if dt <= 0.0 {
                return Err(invalid(
                    "field caloric chart requires strictly increasing temperature",
                ));
            }
            if pair[1].temperature_k >= range[0] && pair[0].temperature_k <= range[1] {
                let capacity = finite(
                    cell.mass_kg
                        * (pair[1].specific_enthalpy_j_kg - pair[0].specific_enthalpy_j_kg)
                        / dt,
                )?;
                minimum_caloric_capacity = minimum_caloric_capacity.min(capacity);
            }
        }
        let minimum_field_capacity = binding
            .law
            .minimum_heat_capacity_j_k(bound)
            .map_err(field)?;
        if !minimum_caloric_capacity.is_finite()
            || finite(minimum_caloric_capacity + minimum_field_capacity)? <= 0.0
        {
            return Err(invalid(
                "caloric plus field heat capacity must stay positive over the declared state envelope",
            ));
        }
        let coordinate = *initial
            .get(&binding.branch)
            .ok_or_else(|| invalid("thermal field mapping is not an original storage branch"))?;
        if coordinate.abs() > bound {
            return Err(invalid(
                "accepted field coordinate exceeds the thermal capacity envelope",
            ));
        }
        result.push(InitialField {
            temperature_k: phase.temperature_k(),
            coordinate,
            state: current.state(coordinate).map_err(field)?,
        });
    }
    check(cancel)?;
    Ok(result)
}

/// Shared source admission and reversible transfer for original circuit owners.
/// All trials refer to these SAME accepted field coordinates and caloric cells.
pub(crate) struct PreparedFieldThermal<'a> {
    bindings: &'a [FieldThermalBinding],
    initial: Vec<InitialField>,
}
impl<'a> PreparedFieldThermal<'a> {
    pub fn new(
        material: &MaterialNetwork<'_>,
        thermal: &DepositionSimulation<'_>,
        circuit: FieldCircuitView<'_>,
        bindings: &'a [FieldThermalBinding],
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<Self, Error> {
        Ok(Self {
            bindings,
            initial: admit(material, thermal, &circuit, bindings, cancel)?,
        })
    }

    pub fn refresh(
        &self,
        trial: &DepositionSimulation<'_>,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<(Vec<MaterialStorage>, Vec<ResolvedThermalFieldStorage>), Error> {
        let mut sources = Vec::with_capacity(self.bindings.len());
        let mut samples = Vec::with_capacity(self.bindings.len());
        for binding in self.bindings {
            check(cancel)?;
            let sample = binding
                .law
                .at_temperature(trial.temperature_k(binding.cell)?)
                .map_err(field)?;
            sources.push(MaterialStorage {
                branch: binding.branch,
                element: element(sample.element()),
            });
            samples.push(sample);
        }
        Ok((sources, samples))
    }

    #[allow(clippy::too_many_arguments)]
    pub fn account<S>(
        &self,
        samples: Vec<ResolvedThermalFieldStorage>,
        endpoint: FieldCircuitView<'_>,
        native_parameter_work_j: f64,
        dt: f64,
        network: S,
        sources_w: &mut BTreeMap<CellId, f64>,
        tolerance_j: f64,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<FieldThermalSource<S>, Error> {
        if samples.len() != self.bindings.len() || !dt.is_finite() || dt <= 0.0 {
            return Err(invalid(
                "complete field endpoint and positive interval required",
            ));
        }
        let endpoints = coordinates(&endpoint)?;
        let mut fields = Vec::with_capacity(self.bindings.len());
        let mut internal = [0.0; 2];
        let mut parameter_work = 0.0;
        let mut reversible_heat_j = 0.0;
        for ((binding, before), sample) in self.bindings.iter().zip(&self.initial).zip(samples) {
            check(cancel)?;
            let coordinate = *endpoints
                .get(&binding.branch)
                .ok_or_else(|| invalid("native endpoint lost an original field coordinate"))?;
            if coordinate.abs() > binding.maximum_abs_coordinate {
                return Err(invalid(
                    "proposed field coordinate exceeds the thermal capacity envelope",
                ));
            }
            let frozen = sample.state(before.coordinate).map_err(field)?;
            let after = sample.state(coordinate).map_err(field)?;
            let work = finite(frozen.free_energy_j - before.state.free_energy_j)?;
            let delta_ts = finite(
                sample.temperature_k() * after.entropy_j_k
                    - before.temperature_k * before.state.entropy_j_k,
            )?;
            let heat = finite(-work - delta_ts)?;
            let power = finite(heat / dt)?;
            if heat != 0.0 && power == 0.0 {
                return Err(invalid(
                    "field reversible heat disappeared during power conversion",
                ));
            }
            let destination = sources_w.entry(binding.cell).or_insert(0.0);
            *destination = finite(*destination + power)?;
            internal[0] = finite(internal[0] + before.state.internal_energy_j)?;
            internal[1] = finite(internal[1] + after.internal_energy_j)?;
            parameter_work = finite(parameter_work + work)?;
            reversible_heat_j = finite(reversible_heat_j + dt * power)?;
            fields.push(FieldThermalSample {
                branch: binding.branch,
                cell: binding.cell,
                temperature_k: [before.temperature_k, sample.temperature_k()],
                coordinate: [before.coordinate, coordinate],
                state: [before.state, after],
                material: sample.element().clone(),
                parameter_work_j: work,
                reversible_heat_j: heat,
            });
        }
        if finite(parameter_work - native_parameter_work_j)?.abs() > tolerance_j {
            return Err(invalid(
                "original circuit parameter work and sourced thermodynamic field change disagree",
            ));
        }
        Ok(FieldThermalSource {
            network,
            fields,
            field_internal_energy_j: internal,
            reversible_heat_j,
        })
    }
}

pub(super) fn close_report<D, S>(
    coupled: CircuitEnthalpyReport<FieldThermalSource<S>, D>,
    electrical_work_j: f64,
    tolerance_j: f64,
) -> Result<FieldThermalReport<D, S>, Error> {
    close_report_with_internal_heat(coupled, electrical_work_j, 0.0, tolerance_j)
}

// Internal exchanges may have a small nonzero numerical net. Never relabel it
// as supplied heat or redistribute it into the cells to force a zero residual.
fn close_report_with_internal_heat<D, S>(
    coupled: CircuitEnthalpyReport<FieldThermalSource<S>, D>,
    electrical_work_j: f64,
    internal_heat_j: f64,
    tolerance_j: f64,
) -> Result<FieldThermalReport<D, S>, Error> {
    let external_heat_j = finite(
        coupled.non_electrical_heat_j - coupled.source.reversible_heat_j - internal_heat_j,
    )?;
    // Reconstruct physical U at both actual endpoints; the native electrical
    // report instead measures isothermal psi after the coefficient update.
    let total_energy_residual_j = finite(
        coupled.source.field_internal_energy_j[1] - coupled.source.field_internal_energy_j[0]
            + coupled.thermal.stored_energy_change_j
            - electrical_work_j
            - external_heat_j,
    )?;
    if total_energy_residual_j.abs() > tolerance_j {
        return Err(invalid(
            "physical field and caloric internal energy do not close the external supply",
        ));
    }
    Ok(FieldThermalReport {
        coupled,
        external_heat_j,
        total_energy_residual_j,
    })
}

/// Evolve source-dependent C/L, native circuit charge/flux, and caloric enthalpy.
///
/// Geometry/topology stay fixed. The original prepared network owns all
/// resistor and conductive-link source selection; the original thermal
/// backward-Euler/Newton owner solves every field/circuit/material trial from
/// the SAME accepted history. All C/L branches need explicit thermal material
/// cells; no implicit isothermal reservoir pays for changing field entropy.
///
/// Initial and final q/lambda must stay within the caller's capacity envelope.
/// Caloric plus field capacity is admitted throughout its temperature range.
/// Native source-constraint refusals (including an ideal-voltage-clamped
/// capacitor whose q cannot remain fixed during rebinding) propagate.
/// Every external source jump or thermal topology event requires an interval
/// cut as in the original circuit/enthalpy API. This is not an unconditional
/// passivity or temporal-error certificate; refine the physical timestep.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate<'m>(
    material: &MaterialNetwork<'_>,
    thermal: &DepositionSimulation<'m>,
    target_s: f64,
    circuit: &MaterialCircuit,
    ports: &[f64],
    bindings: &[FieldThermalBinding],
    electrical_bases: &[(usize, ElectricalPropertyBasis)],
    controls: CircuitEnthalpySolve,
    initial_energy_increments_j: &[f64],
    cancelled: impl FnMut() -> bool,
) -> Result<FieldThermalCandidate<'m>, Error> {
    let cancel = RefCell::new(cancelled);
    check(&cancel)?;
    super::admit_driver_unknowns(thermal, controls, initial_energy_increments_j)?;
    let prepared = material
        .prepare_with_bases(thermal, electrical_bases, || (*cancel.borrow_mut())())
        .map_err(network)?;
    if !circuit.contacts().is_empty() {
        return Err(invalid("field thermal circuit requires bulk-only resistors"));
    }
    let destinations = prepared.resistor_destinations();
    driven_candidate(
        material, thermal, target_s, circuit, ports, bindings, &destinations, &[],
        controls, initial_energy_increments_j,
        |trial| {
            Ok(super::contacts::ContactThermalCoefficients {
                bulk: prepared.resolve(trial, || (*cancel.borrow_mut())()).map_err(network)?,
                contacts: Vec::new(),
            })
        },
        |_| 0.0,
        || (*cancel.borrow_mut())(),
    )
}

// Share native field admission, conserved-coordinate rebinding and circuit
// advancement across bulk, electrical-contact and complete heat networks.
// Only the source payload/evaluator and declared internal heat net vary.
#[allow(clippy::too_many_arguments)]
pub(crate) fn driven_candidate<'m, S, E: From<FieldThermalError>>(
    material: &MaterialNetwork<'_>,
    thermal: &DepositionSimulation<'m>,
    target_s: f64,
    circuit: &MaterialCircuit,
    ports: &[f64],
    bindings: &[FieldThermalBinding],
    destinations: &[super::ResistorCell],
    routes: &[super::contacts::ContactHeatCells],
    controls: CircuitEnthalpySolve,
    guess_j: &[f64],
    resolve: impl Fn(&DepositionSimulation<'m>) -> Result<
        super::contacts::ContactThermalCoefficients<S>, CircuitEnthalpyError<E>>,
    internal_power_w: impl Fn(&S) -> f64,
    cancelled: impl FnMut() -> bool,
) -> Result<FieldThermalCandidate<'m, MaterialCircuit, MaterialCircuitStep, S>, CircuitEnthalpyError<E>> {
    let cancel = RefCell::new(cancelled);
    check(&cancel).map_err(map_field_error)?;
    super::admit_driver_unknowns(thermal, controls, guess_j)?;
    let mut initial = resolve(thermal)?;
    initial.bulk.resistors.sort_unstable_by_key(|r| r.branch);
    // Source pinning changes policy receipts, not physical material state.
    // Keep the original exact source, complete key, basis and geometry checks.
    if initial.bulk.resistors.len() != circuit.resistors().len()
        || initial.bulk.resistors.iter().zip(circuit.resistors()).any(|(a, b)| {
            let x = &a.conductor;
            let y = &b.conductor;
            a.branch != b.branch || x.source_basis() != y.source_basis()
                || x.source_key() != y.source_key()
                || x.length_m().to_bits() != y.length_m().to_bits()
                || x.area_m2().to_bits() != y.area_m2().to_bits()
                || x.resistance_ohm().to_bits() != y.resistance_ohm().to_bits()
                || !same_material(x.material(), y.material())
        }) {
        return Err(CircuitEnthalpyError::Invalid(
            "field circuit must retain its original bulk sources and geometry",
        ));
    }
    super::contacts::admit_sources(circuit, &initial.contacts, || (*cancel.borrow_mut())())?;
    initial.contacts.sort_unstable_by_key(|contact| contact.id());
    for (actual, expected) in circuit.contacts().iter().zip(&initial.contacts) {
        check(&cancel).map_err(map_field_error)?;
        if actual.id() != expected.id() || actual.branch() != expected.branch()
            || actual.area_m2().to_bits() != expected.area_m2().to_bits()
            || actual.resistance_ohm().to_bits() != expected.resistance_ohm().to_bits()
            || actual.source().resolved().query() != expected.source().resolved().query()
        {
            return Err(CircuitEnthalpyError::Invalid(
                "field circuit and original interface state or geometry disagree",
            ));
        }
    }
    let fields = PreparedFieldThermal::new(material, thermal, circuit.into(), bindings, &cancel)
        .map_err(map_field_error)?;
    let dt = target_s - thermal.time_s();
    let (next, electrical, coupled) = coupled_driver_candidate(
        thermal, target_s, CircuitEnthalpyInput { circuit, ports, destinations },
        routes, controls, guess_j,
        |trial, _remaining| {
            check(&cancel).map_err(map_field_error)?;
            let coefficients = resolve(trial)?;
            let mut bulk = coefficients.bulk;
            let (storage, samples) = fields.refresh(trial, &cancel).map_err(map_field_error)?;
            let rebound = circuit.rebind_storage_candidate(
                storage, ports, || (*cancel.borrow_mut())(),
            )?;
            let candidate = super::contacts::advance_circuit(
                rebound.circuit(), bulk.resistors, coefficients.contacts, ports, dt,
                controls.circuit, || (*cancel.borrow_mut())(),
            )?;
            let source = fields.account(
                samples, candidate.circuit().into(), rebound.parameter_work_j(), dt,
                bulk.source, &mut bulk.sources_w, controls.total_energy_tolerance_j, &cancel,
            ).map_err(map_field_error)?;
            Ok(DrivenCoefficients {
                circuit: candidate, sources_w: bulk.sources_w,
                conductances_w_k: bulk.conductances_w_k, source, circuit_trials: 1,
            })
        },
        || (*cancel.borrow_mut())(),
    )?;
    let supplied = coupled.electrical.native.supplied;
    let internal_heat = finite(dt * internal_power_w(&coupled.source.network))?;
    let report = close_report_with_internal_heat(
        coupled, supplied, internal_heat, controls.total_energy_tolerance_j,
    ).map_err(map_field_error)?;
    check(&cancel).map_err(map_field_error)?;
    Ok(FieldThermalCandidate::new(next, electrical, report))
}

fn map_field_error<E: From<FieldThermalError>>(error: Error) -> CircuitEnthalpyError<E> {
    match error {
        CircuitEnthalpyError::Source(error) => CircuitEnthalpyError::Source(error.into()),
        CircuitEnthalpyError::Invalid(message) => CircuitEnthalpyError::Invalid(message),
        CircuitEnthalpyError::Circuit(error) => CircuitEnthalpyError::Circuit(error),
        CircuitEnthalpyError::Thermal(error) => error.into(),
        CircuitEnthalpyError::Newton(error) => CircuitEnthalpyError::Newton(error),
        CircuitEnthalpyError::NotConverged(report) => CircuitEnthalpyError::NotConverged(report),
        CircuitEnthalpyError::Interrupted(progress) => CircuitEnthalpyError::Interrupted(progress),
        CircuitEnthalpyError::CircuitBudget => CircuitEnthalpyError::CircuitBudget,
        CircuitEnthalpyError::Cancelled => CircuitEnthalpyError::Cancelled,
    }
}
