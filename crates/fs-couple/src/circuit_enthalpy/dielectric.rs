//! Debye polarization on the existing circuit/enthalpy residual.
//!
//! Each native circuit trial starts from the same accepted charge/flux history
//! with source coefficients at the trial thermal endpoint. Native resistor loss
//! heats its explicit receiver once. The dielectric adds the distinct reversible
//! transfer Q_rev=-delta(psi_parameter)-delta(T*S), so the enclosing account
//! closes caloric energy plus field U=psi+T*S against external supply.
//! The original nonlinear profile has temperature-independent potentials and
//! psi=U, S=0, with tau(T) changing relaxation. The explicit reversible-storage
//! entry point additionally admits a source-backed D_inf(E,T)=s(T)D_ref(E)
//! potential and affine polarization increments, including their entropy and
//! fixed-charge field capacity. The linear background is never an energy proxy.
//!
//! This shares the original backward-Euler/Newton and circuit midpoint owners.
//! The coupled method is generally first order. Fixed solid volume, positive
//! combined heat capacity and the admitted affine storage paths are explicit;
//! there is no assumed isothermal reservoir or unaccounted ideal-source impulse.

use std::{cell::RefCell, collections::BTreeMap};

use fs_material::field_storage::{
    FieldStorageError,
    dielectric::reversible::thermal::ThermalReversibleDielectricMaterial,
    dielectric::thermal::{ResolvedThermalDebyeCapacitor, ThermalDebyeCapacitor},
    thermal::FieldThermalState,
};
use fs_phs::circuit::Branch;
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};

use super::{
    CircuitEnthalpyError, CircuitEnthalpyReport, CircuitEnthalpySolve, DrivenCoefficients,
    DriverAccount, ThermalDriver, finite,
    material::{MaterialNetwork, MaterialNetworkError, MaterialNetworkSample},
};
use crate::material_circuit::capacitance::dielectric::MaterialDielectricCapacitor;
use crate::material_circuit::dielectric::{DebyeCircuit, DebyeCircuitCandidate, DebyeCircuitStep};

/// One complete uniform dielectric volume and its original thermal source law.
/// Every polarization process belongs to the same material cell. Charge bounds
/// admit a positive caloric-plus-field heat capacity over the whole T interval.
/// For nonlinear instantaneous storage, this law owns the original pole/rate
/// paths; the explicit reversible-storage entry point separately owns the
/// nonlinear instantaneous potential. The legacy nonlinear path requires zero
/// storage slopes. Neither counts the linear background as a second field store.
pub struct DebyeThermalBinding<'a> {
    pub cell: CellId,
    pub law: &'a ThermalDebyeCapacitor,
    pub maximum_instantaneous_charge_c: f64,
    pub maximum_polarization_charges_c: &'a [f64],
}

#[derive(Debug)]
pub enum DebyeThermalError {
    Material(MaterialNetworkError),
    Field(FieldStorageError),
    Invalid(&'static str),
}
impl core::fmt::Display for DebyeThermalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "dielectric circuit/enthalpy: {self:?}")
    }
}
impl std::error::Error for DebyeThermalError {}
type Error = CircuitEnthalpyError<DebyeThermalError>;

/// Source material and thermodynamic accounts at the accepted endpoints.
pub struct DebyeThermalSource {
    pub network: MaterialNetworkSample,
    pub cell: CellId,
    pub temperature_k: [f64; 2],
    pub dielectric_state: [FieldThermalState; 2],
    /// Signed fixed-charge change of the original Helmholtz potential [J].
    pub parameter_work_j: f64,
    /// Stable change of the field's T*S contribution [J].
    pub thermal_offset_change_j: f64,
    /// Signed internal transfer to the caloric chart, separate from dissipation.
    pub reversible_heat_j: f64,
    /// Stable dielectric internal-energy increment, including charge evolution.
    pub dielectric_internal_energy_change_j: f64,
}

/// The original numerical account plus physical internal-energy closure.
pub struct DebyeThermalReport {
    /// Native isothermal electrical account and original thermal residual.
    /// Its non-electrical source total includes the explicit reversible transfer.
    pub coupled: CircuitEnthalpyReport<DebyeThermalSource, DebyeCircuitStep>,
    /// External heat only: reversible field transfer has been removed [J].
    pub external_heat_j: f64,
    /// All electrical internal storage plus original caloric storage [J].
    pub total_energy_residual_j: f64,
}

/// Both histories publish together; all reported heat is already deposited.
pub struct DebyeThermalCandidate<'m> {
    thermal: DepositionSimulation<'m>,
    circuit: DebyeCircuit,
    report: DebyeThermalReport,
}
impl<'m> DebyeThermalCandidate<'m> {
    pub const fn thermal(&self) -> &DepositionSimulation<'m> {
        &self.thermal
    }
    pub const fn circuit(&self) -> &DebyeCircuit {
        &self.circuit
    }
    pub const fn report(&self) -> &DebyeThermalReport {
        &self.report
    }
    pub fn into_parts(self) -> (DepositionSimulation<'m>, DebyeCircuit, DebyeThermalReport) {
        (self.thermal, self.circuit, self.report)
    }
}

struct Driver(DebyeCircuitCandidate);
impl ThermalDriver for Driver {
    type Model = DebyeCircuit;
    type Step = DebyeCircuitStep;
    fn account(&self) -> DriverAccount<'_> {
        let step = self.0.step();
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: step.native.dissipated,
            numerical_dissipation_j: 0.0,
            stored_change_j: step.native.delta_h,
            supplied_j: step.native.supplied,
        }
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        self.0.into_parts()
    }
}
fn invalid(message: &'static str) -> Error {
    CircuitEnthalpyError::Source(DebyeThermalError::Invalid(message))
}
fn field(error: FieldStorageError) -> Error {
    CircuitEnthalpyError::Source(DebyeThermalError::Field(error))
}
fn material(error: MaterialNetworkError) -> Error {
    match error {
        MaterialNetworkError::Thermal(DepositionError::Cancelled) => {
            CircuitEnthalpyError::Cancelled
        }
        other => CircuitEnthalpyError::Source(DebyeThermalError::Material(other)),
    }
}
fn check(cancel: &RefCell<impl FnMut() -> bool>) -> Result<(), Error> {
    if (*cancel.borrow_mut())() {
        Err(CircuitEnthalpyError::Cancelled)
    } else {
        Ok(())
    }
}
fn close(a: f64, b: f64) -> bool {
    a.is_finite() && b.is_finite() && (a - b).abs() <= 128.0 * f64::EPSILON * a.abs().max(b.abs())
}
fn displacement(charge: f64, area: f64) -> Result<f64, Error> {
    let result = finite(charge / area)?;
    if charge != 0.0 && result == 0.0 {
        return Err(invalid("nonlinear thermal displacement underflows"));
    }
    Ok(result)
}
fn check_charges(
    binding: &DebyeThermalBinding<'_>,
    circuit: &DebyeCircuit,
) -> Result<Vec<f64>, Error> {
    let bound = binding.maximum_instantaneous_charge_c;
    let polarization = circuit.polarization_charges_c();
    if !bound.is_finite()
        || bound < 0.0
        || circuit.instantaneous_charge_c().abs() > bound
        || polarization.len() != binding.maximum_polarization_charges_c.len()
        || polarization
            .iter()
            .zip(binding.maximum_polarization_charges_c)
            .any(|(q, bound)| !bound.is_finite() || *bound < 0.0 || q.abs() > *bound)
    {
        return Err(invalid(
            "dielectric charges must remain within the complete thermal capacity envelope",
        ));
    }
    Ok(polarization)
}

fn instantaneous_at_temperature(
    circuit: &DebyeCircuit,
    axis: &str,
    temperature: f64,
    thermal_law: Option<&ThermalReversibleDielectricMaterial>,
) -> Result<Option<MaterialDielectricCapacitor>, Error> {
    let Some(original) = circuit.reversible_instantaneous() else {
        return Ok(None);
    };
    if let Some(law) = thermal_law {
        let resolved = law.at_temperature(temperature).map_err(field)?;
        return Ok(Some(original.at_resolved_temperature(&resolved)?));
    }
    let point = original
        .material()
        .material()
        .query()
        .with(axis, temperature)
        .map_err(|_| invalid("nonlinear dielectric temperature context is unrepresentable"))?;
    Ok(Some(original.at_temperature_query_point(&point)?))
}

fn physical_state(
    circuit: &DebyeCircuit,
    thermal: &ResolvedThermalDebyeCapacitor,
    polarization: &[f64],
) -> Result<FieldThermalState, Error> {
    if let Some(capacitor) = circuit.reversible_instantaneous()
        && let Some(instantaneous) = capacitor.thermal_material()
    {
        let volume = finite(capacitor.area_m2() * capacitor.gap_m())?;
        let d = displacement(circuit.instantaneous_charge_c().abs(), capacitor.area_m2())?;
        let instant = instantaneous.state(d, volume).map_err(field)?;
        let poles = thermal.state(0.0, polarization).map_err(field)?;
        return Ok(FieldThermalState {
            free_energy_j: finite(instant.free_energy_j + poles.free_energy_j)?,
            entropy_j_k: finite(instant.entropy_j_k + poles.entropy_j_k)?,
            internal_energy_j: finite(instant.internal_energy_j + poles.internal_energy_j)?,
            heat_capacity_j_k: finite(instant.heat_capacity_j_k + poles.heat_capacity_j_k)?,
        });
    }
    if circuit.reversible_instantaneous().is_some() {
        if thermal
            .capacitance_slopes_f_k()
            .iter()
            .any(|&slope| slope != 0.0)
        {
            return Err(invalid(
                "nonlinear dielectric thermal storage requires temperature-independent potentials",
            ));
        }
        let energy = finite(circuit.dielectric_energy_j()?)?;
        Ok(FieldThermalState {
            free_energy_j: energy,
            entropy_j_k: 0.0,
            internal_energy_j: energy,
            heat_capacity_j_k: 0.0,
        })
    } else {
        thermal
            .state(circuit.instantaneous_charge_c(), polarization)
            .map_err(field)
    }
}

fn admit(
    network: &MaterialNetwork<'_>,
    thermal: &DepositionSimulation<'_>,
    circuit: &DebyeCircuit,
    binding: &DebyeThermalBinding<'_>,
    instantaneous_law: Option<&ThermalReversibleDielectricMaterial>,
    cancel: &RefCell<impl FnMut() -> bool>,
) -> Result<ResolvedThermalDebyeCapacitor, Error> {
    check(cancel)?;
    match (circuit.thermal_reversible_storage(), instantaneous_law) {
        (Some(authority), Some(law))
            if authority.polarization() == binding.law && authority.instantaneous() == law => {}
        (None, None) => {}
        _ => {
            return Err(invalid(
                "finite nonlinear thermal feedback requires the circuit's complete retained potential authority",
            ));
        }
    }
    check_charges(binding, circuit)?;
    let cell = thermal
        .cells()
        .get(&binding.cell)
        .ok_or_else(|| invalid("dielectric thermal receiver is inactive"))?;
    let body = network
        .cells
        .iter()
        .find(|body| body.cell == binding.cell)
        .ok_or_else(|| invalid("dielectric receiver lacks its original caloric binding"))?;
    let phase = body
        .curve
        .state_at_specific_enthalpy(cell.specific_h_j_kg)
        .map_err(|e| material(MaterialNetworkError::Phase(e)))?;
    let point = body
        .point
        .at(phase)
        .map_err(|e| material(MaterialNetworkError::Material(e)))?;
    let original = binding.law.at_source(circuit.source()).map_err(field)?;
    let volume = finite(circuit.source().area_m2() * circuit.source().gap_m())?;
    if body.point.material.content_hash() != binding.law.card().content_hash()
        || body.point.temperature_axis != binding.law.temperature_axis()
        || point != circuit.source().material().query()
        || !close(original.temperature_k(), phase.temperature_k())
        || !close(volume, cell.geometry.volume_m3)
    {
        return Err(invalid(
            "dielectric and caloric owners must share the original card, actual context and whole volume",
        ));
    }
    let range = binding.law.temperature_range_k();
    if let Some(instantaneous) = circuit.reversible_instantaneous() {
        if binding.maximum_instantaneous_charge_c > instantaneous.maximum_charge_c()
            || instantaneous_law.is_none()
                && binding
                    .law
                    .capacitance_slopes_f_k()
                    .iter()
                    .any(|&slope| slope != 0.0)
        {
            return Err(invalid(
                "nonlinear dielectric thermal envelope requires supported charge and its admitted storage profile",
            ));
        }
        // The rate/caloric charts alone do not establish validity of the
        // original nonlinear D(E) claim throughout this temperature interval.
        for temperature in range {
            check(cancel)?;
            let refreshed = instantaneous_at_temperature(
                circuit,
                binding.law.temperature_axis(),
                temperature,
                instantaneous_law,
            )?
            .ok_or_else(|| invalid("nonlinear dielectric source disappeared"))?;
            if instantaneous_law.is_none() && refreshed.curve() != instantaneous.curve() {
                return Err(invalid(
                    "nonlinear thermal feedback cannot change the original field potential",
                ));
            }
        }
    }
    let knots = body.curve.knots();
    if range[0] < knots[0].temperature_k || range[1] > knots[knots.len() - 1].temperature_k {
        return Err(invalid(
            "the complete dielectric temperature interval needs caloric support",
        ));
    }
    let mut minimum_capacity = f64::INFINITY;
    for knot in knots {
        check(cancel)?;
        if knot.liquid_mass_fraction != 0.0
            || !close(cell.mass_kg / knot.bulk_density_kg_m3, volume)
        {
            return Err(invalid(
                "dielectric thermal feedback requires fixed fully solid geometry and density",
            ));
        }
    }
    for pair in knots.windows(2) {
        check(cancel)?;
        let dt = pair[1].temperature_k - pair[0].temperature_k;
        if dt <= 0.0 {
            return Err(invalid(
                "dielectric caloric temperature must be strictly increasing",
            ));
        }
        if pair[1].temperature_k >= range[0] && pair[0].temperature_k <= range[1] {
            minimum_capacity = minimum_capacity.min(finite(
                cell.mass_kg * (pair[1].specific_enthalpy_j_kg - pair[0].specific_enthalpy_j_kg)
                    / dt,
            )?);
        }
    }
    let field_capacity = if let Some(law) = instantaneous_law {
        finite(
            law.minimum_heat_capacity_j_k(
                displacement(
                    binding.maximum_instantaneous_charge_c,
                    circuit.source().area_m2(),
                )?,
                volume,
            )
            .map_err(field)?
                + binding
                    .law
                    .minimum_heat_capacity_j_k(0.0, binding.maximum_polarization_charges_c)
                    .map_err(field)?,
        )?
    } else {
        binding
            .law
            .minimum_heat_capacity_j_k(
                binding.maximum_instantaneous_charge_c,
                binding.maximum_polarization_charges_c,
            )
            .map_err(field)?
    };
    if !minimum_capacity.is_finite() || finite(minimum_capacity + field_capacity)? <= 0.0 {
        return Err(invalid(
            "caloric plus all dielectric polarization heat capacities must stay positive",
        ));
    }
    Ok(original)
}

/// Joint dielectric polarization, source R(T), caloric storage and conduction.
///
/// The material network binds every active thermal cell/link and exactly the
/// ORIGINAL non-polarization resistors. The dielectric's internal friction arms
/// belong to its own material cell and must not be supplied as conductor data.
/// Other original C/L coefficients stay fixed; every charge and flux persists.
/// The existing source-selection, fixed geometry and thermal-event rules apply.
///
/// Every trial uses the accepted histories. Parameter work and delta(T*S) are
/// internal reversible transfer; native loss is irreversible heat. The final
/// audit uses physical U plus caloric energy. Ideal clamps that require a charge
/// impulse during the coefficient update refuse. This endpoint coefficient
/// scheme needs timestep refinement and does not certify unconditional entropy
/// production, phase transitions or a changing capacitor geometry.
/// Nonlinear instantaneous D_inf(E) with constant polarization increments uses
/// its actual native energy and zero field entropy. Positive supported tau(T)
/// may vary. This legacy entry point refuses changing nonlinear field or pole
/// storage; use [`step_candidate_with_reversible_storage`] with the circuit's
/// retained thermal potential authority for the admitted temperature-dependent
/// profile.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate<'m>(
    network: &MaterialNetwork<'_>,
    thermal: &DepositionSimulation<'m>,
    target_s: f64,
    circuit: &DebyeCircuit,
    ports: &[f64],
    binding: DebyeThermalBinding<'_>,
    controls: CircuitEnthalpySolve,
    initial_energy_increments_j: &[f64],
    cancelled: impl FnMut() -> bool,
) -> Result<DebyeThermalCandidate<'m>, Error> {
    step_with_storage(
        network,
        thermal,
        target_s,
        circuit,
        ports,
        binding,
        None,
        controls,
        initial_energy_increments_j,
        cancelled,
    )
}

/// Finite heat feedback for an explicitly sourced D_inf(E,T)=s(T)D_ref(E)
/// potential and original affine Debye polarization storage. The circuit must
/// retain these SAME admitted laws. Parameter work and stable delta(T*S) supply
/// reversible heat; physical native polarization/conductor losses stay separate.
/// Uses the existing electrical and caloric residual, with no alternate stepper.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate_with_reversible_storage<'m>(
    network: &MaterialNetwork<'_>,
    thermal: &DepositionSimulation<'m>,
    target_s: f64,
    circuit: &DebyeCircuit,
    ports: &[f64],
    binding: DebyeThermalBinding<'_>,
    instantaneous_law: &ThermalReversibleDielectricMaterial,
    controls: CircuitEnthalpySolve,
    initial_energy_increments_j: &[f64],
    cancelled: impl FnMut() -> bool,
) -> Result<DebyeThermalCandidate<'m>, Error> {
    step_with_storage(
        network,
        thermal,
        target_s,
        circuit,
        ports,
        binding,
        Some(instantaneous_law),
        controls,
        initial_energy_increments_j,
        cancelled,
    )
}

#[allow(clippy::too_many_arguments)]
fn step_with_storage<'m>(
    network: &MaterialNetwork<'_>,
    thermal: &DepositionSimulation<'m>,
    target_s: f64,
    circuit: &DebyeCircuit,
    ports: &[f64],
    binding: DebyeThermalBinding<'_>,
    instantaneous_law: Option<&ThermalReversibleDielectricMaterial>,
    controls: CircuitEnthalpySolve,
    initial_energy_increments_j: &[f64],
    cancelled: impl FnMut() -> bool,
) -> Result<DebyeThermalCandidate<'m>, Error> {
    let cancel = RefCell::new(cancelled);
    check(&cancel)?;
    super::admit_driver_unknowns(thermal, controls, initial_energy_increments_j)?;
    let dt = target_s - thermal.time_s();
    if !dt.is_finite() || dt <= 0.0 {
        return Err(invalid("positive dielectric thermal interval required"));
    }
    let prepared = network
        .prepare(thermal, || (*cancel.borrow_mut())())
        .map_err(material)?;
    let initial_material = prepared
        .resolve(thermal, || (*cancel.borrow_mut())())
        .map_err(material)?;
    let original = admit(
        network,
        thermal,
        circuit,
        &binding,
        instantaneous_law,
        &cancel,
    )?;
    let original_polarization = check_charges(&binding, circuit)?;
    let original_state = physical_state(circuit, &original, &original_polarization)?;
    let nonlinear = circuit.reversible_instantaneous().is_some();
    let original_count = circuit.original_branch_count();
    let resistor_count = circuit.graph().branches[..original_count]
        .iter()
        .filter(|(_, _, b)| matches!(b, Branch::Resistor { .. }))
        .count();
    if initial_material.resistors.len() != resistor_count {
        return Err(invalid(
            "every original resistor needs its material and thermal destination",
        ));
    }
    for resistor in &initial_material.resistors {
        check(&cancel)?;
        match circuit.graph().branches.get(resistor.branch) {
            Some((_, _, Branch::Resistor { ohms }))
                if resistor.branch < original_count
                    && close(*ohms, resistor.conductor.resistance_ohm()) => {}
            _ => {
                return Err(invalid(
                    "accepted resistor and its initial caloric source must agree",
                ));
            }
        }
    }
    let destinations: BTreeMap<_, _> = prepared
        .resistor_destinations()
        .iter()
        .map(|r| (r.branch, r.cell))
        .collect();
    if destinations.len() != resistor_count {
        return Err(invalid(
            "one unique receiver per original resistor required",
        ));
    }
    let (next_thermal, next_circuit, coupled) = super::cell_driver_candidate(
        thermal,
        target_s,
        controls,
        initial_energy_increments_j,
        |trial, _remaining| {
            check(&cancel)?;
            let mut bulk = prepared
                .resolve(trial, || (*cancel.borrow_mut())())
                .map_err(material)?;
            let endpoint = binding
                .law
                .at_temperature(trial.temperature_k(binding.cell)?)
                .map_err(field)?;
            let rebound = if let Some(instantaneous) = instantaneous_at_temperature(
                circuit,
                binding.law.temperature_axis(),
                endpoint.temperature_k(),
                instantaneous_law,
            )? {
                if instantaneous_law.is_some() {
                    circuit.rebind_thermal_reversible_sources_and_resistors_candidate(
                        endpoint.source().clone(),
                        instantaneous,
                        bulk.resistors,
                        ports,
                        || (*cancel.borrow_mut())(),
                    )?
                } else {
                    circuit.rebind_reversible_sources_and_resistors_candidate(
                        endpoint.source().clone(),
                        instantaneous,
                        bulk.resistors,
                        ports,
                        || (*cancel.borrow_mut())(),
                    )?
                }
            } else {
                circuit.rebind_source_and_resistors_candidate(
                    endpoint.source().clone(),
                    bulk.resistors,
                    ports,
                    || (*cancel.borrow_mut())(),
                )?
            };
            let parameter_work_j = rebound.parameter_work_j();
            if nonlinear && instantaneous_law.is_none() && parameter_work_j != 0.0 {
                return Err(invalid(
                    "temperature-independent nonlinear dielectric acquired parameter work",
                ));
            }
            let candidate =
                rebound
                    .circuit()
                    .step_candidate(ports, dt, controls.circuit, || (*cancel.borrow_mut())())?;
            let polarization = check_charges(&binding, candidate.circuit())?;
            let state = physical_state(candidate.circuit(), &endpoint, &polarization)?;
            let thermal_offset_change_j = if instantaneous_law.is_some() {
                let before = circuit.reversible_instantaneous().ok_or_else(|| {
                    invalid("thermal dielectric lost its original instantaneous law")
                })?;
                let after = candidate
                    .circuit()
                    .reversible_instantaneous()
                    .ok_or_else(|| {
                        invalid("thermal dielectric lost its refreshed instantaneous law")
                    })?;
                let initial_field = before
                    .thermal_material()
                    .ok_or_else(|| invalid("original nonlinear field thermal potential absent"))?;
                let final_field = after
                    .thermal_material()
                    .ok_or_else(|| invalid("refreshed nonlinear field thermal potential absent"))?;
                let field_offset = initial_field
                    .thermal_offset_change_j(
                        final_field,
                        displacement(circuit.instantaneous_charge_c(), before.area_m2())?,
                        displacement(
                            candidate.circuit().instantaneous_charge_c(),
                            after.area_m2(),
                        )?,
                        finite(before.area_m2() * before.gap_m())?,
                    )
                    .map_err(field)?;
                let pole_offset = original
                    .thermal_offset_change_j(
                        &endpoint,
                        0.0,
                        &original_polarization,
                        0.0,
                        &polarization,
                    )
                    .map_err(field)?;
                finite(field_offset + pole_offset)?
            } else if nonlinear {
                0.0
            } else {
                original
                    .thermal_offset_change_j(
                        &endpoint,
                        circuit.instantaneous_charge_c(),
                        &original_polarization,
                        candidate.circuit().instantaneous_charge_c(),
                        &polarization,
                    )
                    .map_err(field)?
            };
            let reversible_heat_j = finite(-parameter_work_j - thermal_offset_change_j)?;
            let power = finite(reversible_heat_j / dt)?;
            if reversible_heat_j != 0.0 && power == 0.0 {
                return Err(invalid(
                    "dielectric reversible heat disappeared during power conversion",
                ));
            }
            let supplied = bulk.sources_w.entry(binding.cell).or_insert(0.0);
            *supplied = finite(*supplied + power)?;
            let dielectric_internal_energy_change_j = finite(
                parameter_work_j
                    + candidate.step().dielectric_energy_change_j
                    + thermal_offset_change_j,
            )?;
            let source = DebyeThermalSource {
                network: bulk.source,
                cell: binding.cell,
                temperature_k: [original.temperature_k(), endpoint.temperature_k()],
                dielectric_state: [original_state, state],
                parameter_work_j,
                thermal_offset_change_j,
                reversible_heat_j,
                dielectric_internal_energy_change_j,
            };
            Ok(DrivenCoefficients {
                circuit: Driver(candidate),
                sources_w: bulk.sources_w,
                conductances_w_k: bulk.conductances_w_k,
                source,
                circuit_trials: 1,
            })
        },
        |driver: &Driver, duration, checkpoint| {
            checkpoint()?;
            let step = driver.0.step();
            let mut powers = BTreeMap::new();
            let mut deposit = |cell, heat| -> Result<(), Error> {
                checkpoint()?;
                let power = finite(heat / duration)?;
                if heat != 0.0 && power == 0.0 {
                    return Err(invalid(
                        "native dielectric or resistor heat disappeared during deposition",
                    ));
                }
                let total = powers.entry(cell).or_insert(0.0);
                *total = finite(*total + power)?;
                Ok(())
            };
            deposit(binding.cell, step.polarization_heat_j)?;
            if step.other_resistors.len() != destinations.len() {
                return Err(invalid("native resistor heat lost an original destination"));
            }
            for resistor in &step.other_resistors {
                let cell = *destinations.get(&resistor.branch).ok_or_else(|| {
                    invalid("native resistor heat has no original material receiver")
                })?;
                deposit(cell, resistor.heat_j)?;
            }
            Ok(powers)
        },
        || (*cancel.borrow_mut())(),
    )?;
    let external_heat_j = finite(coupled.non_electrical_heat_j - coupled.source.reversible_heat_j)?;
    let total_energy_residual_j = finite(
        coupled.electrical.native.delta_h
            + coupled.source.parameter_work_j
            + coupled.source.thermal_offset_change_j
            + coupled.thermal.stored_energy_change_j
            - coupled.electrical.native.supplied
            - external_heat_j,
    )?;
    if total_energy_residual_j.abs() > controls.total_energy_tolerance_j {
        return Err(invalid(
            "dielectric internal and caloric energy do not close the external supply",
        ));
    }
    check(&cancel)?;
    Ok(DebyeThermalCandidate {
        thermal: next_thermal,
        circuit: next_circuit,
        report: DebyeThermalReport {
            coupled,
            external_heat_j,
            total_energy_residual_j,
        },
    })
}
