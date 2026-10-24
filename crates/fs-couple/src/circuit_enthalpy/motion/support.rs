//! Thermal expansion and elasticity of the SAME member driving a field actuator.
//!
//! The shared section law owns per-layer Helmholtz energy and entropy. The
//! existing field-motion descriptor supplies actual electrical/gap/momentum
//! evolution, and the existing finite enthalpy residual supplies heat transport.
//! Each thermal trial restarts from the accepted histories. Rebinding the
//! member at fixed initial gap costs parameter energy P; each layer receives
//! Q_rev = -P - delta(T S). Thus caloric energy plus Psi + T S, field energy and
//! kinetic energy close without relabeling reversible exchange as Joule heat.

use super::super::{
    CircuitEnthalpyError, CircuitEnthalpyInput, CircuitEnthalpyReport,
    CircuitEnthalpySolve, DrivenCoefficients, coupled_driver_candidate, finite,
    material::{MaterialNetwork, MaterialNetworkError, MaterialNetworkSample},
};
use crate::material_circuit::motion::{FieldMotionInput, FieldMotionStep, MaterialFieldMotion};
use crate::material_circuit::storage::MaterialStorageElement;
use fs_matdb::EvaluationDecision;
use fs_material::{
    conductor::ElectricalPropertyBasis,
    field_storage::{ELECTRIC_PERMITTIVITY_PROPERTY, MAGNETIC_PERMEABILITY_PROPERTY},
    section_temperature::thermoelastic::{ConcentricThermoelasticLaw, ConcentricThermoelasticState},
    state_point::MaterialStatePointError,
};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use std::{cell::RefCell, collections::BTreeSet};

/// Original constitutive law and exactly one whole thermal cell per layer.
/// Spatial layer order is retained; there is no averaged expansion coefficient.
pub struct SupportThermalBinding<'a> {
    pub law: ConcentricThermoelasticLaw<'a>,
    pub layer_cells: Vec<CellId>,
    /// Absolute geometry agreement allowance [m3], not an inferred fill factor.
    pub volume_tolerance_m3: f64,
}

#[derive(Debug)]
pub enum SupportThermalError {
    Network(MaterialNetworkError),
    Material(MaterialStatePointError),
    Invalid(&'static str),
}
impl core::fmt::Display for SupportThermalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "thermoelastic actuator support: {self:?}")
    }
}
impl std::error::Error for SupportThermalError {}
type Error = CircuitEnthalpyError<SupportThermalError>;

/// Physical endpoints of one layer, using its reference volume A_i L.
#[derive(Clone, Debug, PartialEq)]
pub struct SupportLayerSample {
    pub cell: CellId,
    pub temperature_k: [f64; 2],
    pub free_energy_j: [f64; 2],
    pub entropy_j_k: [f64; 2],
    /// Frozen-initial-gap change in Helmholtz energy [J].
    pub parameter_work_j: f64,
    /// Signed transfer into this layer's caloric store [J].
    pub reversible_heat_j: f64,
}
#[derive(Clone, Debug, PartialEq)]
pub struct SupportThermalSource {
    pub network: MaterialNetworkSample,
    pub layers: Vec<SupportLayerSample>,
    pub support_internal_energy_j: [f64; 2],
    pub reversible_heat_j: f64,
}
pub struct SupportThermalReport {
    /// Original descriptor and finite enthalpy records. The native elastic
    /// starting energy uses the trial endpoint temperatures at the INITIAL gap.
    /// Its non_electrical_heat_j includes the reversible heat supplied above.
    pub coupled: CircuitEnthalpyReport<SupportThermalSource, FieldMotionStep>,
    /// Authored boundaries and sources, excluding internal Joule/reversible heat.
    pub external_heat_j: f64,
    /// Independent delta(U_caloric + U_support + H_field + K) minus supplies.
    pub total_energy_residual_j: f64,
}
pub struct SupportThermalCandidate<'m> {
    thermal: DepositionSimulation<'m>,
    motion: MaterialFieldMotion,
    report: SupportThermalReport,
}
impl<'m> SupportThermalCandidate<'m> {
    pub const fn thermal(&self) -> &DepositionSimulation<'m> { &self.thermal }
    pub const fn motion(&self) -> &MaterialFieldMotion { &self.motion }
    pub const fn report(&self) -> &SupportThermalReport { &self.report }
    /// Publish all histories together. Both internal heat transfers are present.
    pub fn into_parts(self) -> (DepositionSimulation<'m>, MaterialFieldMotion, SupportThermalReport) {
        (self.thermal, self.motion, self.report)
    }
}

fn invalid(message: &'static str) -> Error {
    CircuitEnthalpyError::Source(SupportThermalError::Invalid(message))
}
fn material(error: MaterialStatePointError) -> Error {
    CircuitEnthalpyError::Source(SupportThermalError::Material(error))
}
fn network(error: MaterialNetworkError) -> Error {
    match error {
        MaterialNetworkError::Thermal(DepositionError::Cancelled) => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(SupportThermalError::Network(other)),
    }
}
fn check(cancel: &RefCell<impl FnMut() -> bool>) -> Result<(), Error> {
    if (*cancel.borrow_mut())() { Err(CircuitEnthalpyError::Cancelled) } else { Ok(()) }
}
fn temperatures(thermal: &DepositionSimulation<'_>, cells: &[CellId]) -> Result<Vec<f64>, Error> {
    cells.iter().map(|&id| thermal.temperature_k(id).map_err(Into::into)).collect()
}
fn strain(motion: &MaterialFieldMotion) -> Result<f64, Error> {
    let support = motion.support().ok_or_else(|| invalid("an existing axial member is required"))?;
    Ok(finite((motion.gap_m() - support.reference_gap_m()) / support.length_m())?)
}

fn admit(
    material_network: &MaterialNetwork<'_>, thermal: &DepositionSimulation<'_>,
    motion: &MaterialFieldMotion, binding: &SupportThermalBinding<'_>,
    cancel: &RefCell<impl FnMut() -> bool>,
) -> Result<ConcentricThermoelasticState, Error> {
    check(cancel)?;
    let support = motion.support().ok_or_else(|| invalid("an existing axial member is required"))?;
    if !support.is_thermally_live() || !motion.circuit().contacts().is_empty()
        || binding.layer_cells.len() != support.section().section().layers().len()
        || !binding.volume_tolerance_m3.is_finite() || binding.volume_tolerance_m3 < 0.0
        || binding.layer_cells.iter().copied().collect::<BTreeSet<_>>().len() != binding.layer_cells.len()
    {
        return Err(invalid("live support, unique complete layer cells, bulk-only electrical heat and finite geometry tolerance required"));
    }
    // The moving element already enforces this at construction. The complete
    // U audit also encloses any OTHER source-bound field stores; their omitted
    // thermal entropy cannot be hidden behind an unchanged coefficient snapshot.
    for storage in motion.circuit().storage() {
        check(cancel)?;
        let (state, name) = match &storage.element {
            MaterialStorageElement::Capacitor(c) => (c.material(), ELECTRIC_PERMITTIVITY_PROPERTY),
            MaterialStorageElement::Inductor(l) => (l.material(), MAGNETIC_PERMEABILITY_PROPERTY),
        };
        if state.property(name).is_none_or(|p|
            p.answer().receipt.decision != EvaluationDecision::ConstantWithinValidity)
        {
            return Err(invalid("all enclosed field storage needs temperature-independent coefficients on the support thermal lane"));
        }
    }
    let current = binding.law.at_temperatures_k(&temperatures(thermal, &binding.layer_cells)?,
        strain(motion)?, support.linear_strain_limit()).map_err(material)?;
    if support.section() != &current.section {
        return Err(invalid("accepted support must match its original law and actual receiving-cell temperatures"));
    }
    for (index, (&id, state)) in binding.layer_cells.iter().zip(&current.layers).enumerate() {
        check(cancel)?;
        let source = material_network.cells.iter().find(|s| s.cell == id)
            .ok_or_else(|| invalid("support layer lacks its original material/caloric cell binding"))?;
        let cell = thermal.cells().get(&id).ok_or_else(|| invalid("support layer cell is inactive"))?;
        let layer = &current.section.section().layers()[index];
        let phase = source.curve.state_at_specific_enthalpy(cell.specific_h_j_kg)
            .map_err(|e| network(MaterialNetworkError::Phase(e)))?;
        let query = source.point.at(phase).map_err(|e| network(MaterialNetworkError::Material(e)))?;
        let volume = finite(state.area_m2 * support.length_m())?;
        let volume_error = finite(cell.geometry.volume_m3 - volume)?.abs();
        let mass_error = finite(cell.mass_kg / state.density_kg_m3 - volume)?.abs();
        if source.point.material.content_hash() != layer.material.card_identity()
            || source.point.temperature_axis != "T" || query != layer.material.query()
            || volume_error > binding.volume_tolerance_m3 || mass_error > binding.volume_tolerance_m3
            || volume <= 0.0 || cell.mass_kg <= 0.0
        {
            return Err(invalid("support mechanics and caloric cell must share exact material context and whole reference inventory"));
        }
        let knots = source.curve.knots();
        let mut minimum_capacity = f64::INFINITY;
        for knot in knots {
            check(cancel)?;
            if knot.liquid_mass_fraction != 0.0 || knot.bulk_density_kg_m3 != state.density_kg_m3 {
                return Err(invalid("live support needs a wholly solid constant-reference-density caloric chart"));
            }
        }
        for pair in knots.windows(2) {
            check(cancel)?;
            let dt = finite(pair[1].temperature_k - pair[0].temperature_k)?;
            if dt <= 0.0 {
                return Err(invalid("live support caloric temperature must increase strictly"));
            }
            let capacity = finite(cell.mass_kg *
                (pair[1].specific_enthalpy_j_kg - pair[0].specific_enthalpy_j_kg) / dt)?;
            minimum_capacity = minimum_capacity.min(capacity);
        }
        // Bound -T*Psi_TT over the caloric interval and every admitted elastic
        // strain. The affine formula is used ONLY as a conservative arithmetic
        // bound here; actual queries still require the shared complete source path.
        let bounds = [knots[0].temperature_k, knots[knots.len()-1].temperature_k];
        let modulus = |t| finite(state.young_modulus_derivative_pa_k.mul_add(
            t - state.temperature_k, state.young_modulus_pa));
        let maximum_modulus = modulus(bounds[0])?.max(modulus(bounds[1])?).max(0.0);
        let curvature_bound = finite((1.0 + 2048.0 * f64::EPSILON) * volume
            * (maximum_modulus * state.alpha_per_k.powi(2)
            + 2.0 * (state.young_modulus_derivative_pa_k * state.alpha_per_k).abs()
                * support.linear_strain_limit()))?;
        let capacity = finite(minimum_capacity - bounds[1] * curvature_bound)?;
        if capacity <= 0.0 {
            return Err(invalid("caloric plus thermoelastic capacity must stay positive over the admitted support envelope"));
        }
    }
    Ok(current)
}

impl MaterialNetwork<'_> {
    /// Evolve heating, layer expansion/stiffness, real gap motion and electrical
    /// response on one accepted clock. Each source layer occupies its complete
    /// fixed reference cell; optional conductive links and authored boundaries
    /// retain the original thermal owner. Scalar density/expansion and globally
    /// affine E(T) use the shared section law, never material-name parameters.
    ///
    /// Geometry follows the original linear axial displacement shape. The
    /// member's consistent reference inertia is fixed; it is not another copy
    /// of the guided body's regional mass. No radial work, remeshing, phase
    /// change, interlayer slip, magnetic hysteresis or moving-gap medium heating
    /// is inferred. Electrical paths retain their authored reference geometry;
    /// strain-dependent resistance and distributed Lorentz work are not added.
    /// Field coefficients remain in their original constant lane.
    /// This endpoint-temperature coupling is generally first order in time.
    #[allow(clippy::too_many_arguments)]
    pub fn field_motion_support_candidate<'m>(
        &self, thermal: &DepositionSimulation<'m>, target_s: f64,
        motion: &MaterialFieldMotion, ports: &[f64], external_force_n: f64,
        electrical_bases: &[(usize, ElectricalPropertyBasis)],
        binding: &SupportThermalBinding<'_>, controls: CircuitEnthalpySolve,
        guess_j: &[f64], cancelled: impl FnMut() -> bool,
    ) -> Result<SupportThermalCandidate<'m>, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        super::super::admit_driver_unknowns(thermal, controls, guess_j)?;
        let prepared = self.prepare_with_bases(thermal, electrical_bases,
            || (*cancel.borrow_mut())()).map_err(network)?;
        let initial = admit(self, thermal, motion, binding, &cancel)?;
        let support = motion.support().expect("admitted support");
        let length = support.length_m();
        let initial_strain = strain(motion)?;
        let destinations = self.resistors.iter().map(|r| r.destination).collect::<Vec<_>>();
        let input = CircuitEnthalpyInput { circuit: motion.circuit(), ports, destinations: &destinations };
        let dt = target_s - thermal.time_s();
        let (next_thermal, next_motion, coupled) = coupled_driver_candidate(
            thermal, target_s, input, &[], controls, guess_j,
            |trial, _remaining| {
                check(&cancel)?;
                let mut coefficients = prepared.resolve(trial, || (*cancel.borrow_mut())()).map_err(network)?;
                let temperatures = temperatures(trial, &binding.layer_cells)?;
                let frozen = binding.law.at_temperatures_k(&temperatures, initial_strain,
                    support.linear_strain_limit()).map_err(material)?;
                let native = motion.step_thermal_support_candidate(frozen.section,
                    coefficients.resistors, FieldMotionInput { electrical_ports: ports,
                        external_force_n, duration_s: dt }, controls.circuit,
                    || (*cancel.borrow_mut())())?;
                let after = binding.law.at_temperatures_k(&temperatures, strain(native.model())?,
                    support.linear_strain_limit()).map_err(material)?;
                if native.model().support().expect("retained support").section() != &after.section {
                    return Err(invalid("native mechanical endpoint and thermoelastic source endpoint disagree"));
                }
                let mut layers = Vec::with_capacity(initial.layers.len());
                let mut internal = [0.0; 2];
                let mut reversible_heat = 0.0;
                let mut free = [0.0; 2];
                for (index, &id) in binding.layer_cells.iter().enumerate() {
                    check(&cancel)?;
                    let before = &initial.layers[index];
                    let held = &frozen.layers[index];
                    let end = &after.layers[index];
                    let psi = [finite(length * before.axial_free_energy_per_reference_length_j_m)?,
                        finite(length * end.axial_free_energy_per_reference_length_j_m)?];
                    let entropy = [finite(length * before.entropy_per_reference_length_j_m_k)?,
                        finite(length * end.entropy_per_reference_length_j_m_k)?];
                    let parameter = finite(length * held.axial_free_energy_per_reference_length_j_m - psi[0])?;
                    let ts = [finite(before.temperature_k * entropy[0])?, finite(end.temperature_k * entropy[1])?];
                    let heat = finite(-parameter - (ts[1] - ts[0]))?;
                    let power = finite(heat / dt)?;
                    if heat != 0.0 && power == 0.0 {
                        return Err(invalid("support reversible heat disappeared during power conversion"));
                    }
                    let source = coefficients.sources_w.entry(id).or_insert(0.0);
                    *source = finite(*source + power)?;
                    reversible_heat = finite(reversible_heat + power * dt)?;
                    for i in 0..2 {
                        internal[i] = finite(internal[i] + psi[i] + ts[i])?;
                    }
                    free[0] = finite(free[0] + psi[0] + parameter)?;
                    free[1] = finite(free[1] + psi[1])?;
                    layers.push(SupportLayerSample { cell: id,
                        temperature_k: [before.temperature_k, end.temperature_k],
                        free_energy_j: psi, entropy_j_k: entropy,
                        parameter_work_j: parameter, reversible_heat_j: heat });
                }
                if free.iter().zip(native.step().elastic_energy_j).any(|(&a, b)|
                    (a - b).abs() > controls.total_energy_tolerance_j)
                {
                    return Err(invalid("native support potential and shared per-layer Helmholtz energy disagree"));
                }
                Ok(DrivenCoefficients { circuit: native, sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k, circuit_trials: 1,
                    source: SupportThermalSource { network: coefficients.source, layers,
                        support_internal_energy_j: internal, reversible_heat_j: reversible_heat } })
            }, || (*cancel.borrow_mut())())?;
        let external_heat = finite(coupled.non_electrical_heat_j - coupled.source.reversible_heat_j)?;
        let step = &coupled.electrical;
        let internal = coupled.source.support_internal_energy_j;
        let residual = finite((next_thermal.energy_j() - thermal.energy_j())
            + (internal[1] - internal[0])
            + (step.electrical_energy_j[1] - step.electrical_energy_j[0])
            + (step.kinetic_energy_j[1] - step.kinetic_energy_j[0])
            - step.electrical_supply_j - step.mechanical_supply_j - external_heat)?;
        if residual.abs() > controls.total_energy_tolerance_j {
            return Err(invalid("support internal energy, field, inertia and caloric energy do not close"));
        }
        check(&cancel)?;
        Ok(SupportThermalCandidate { thermal: next_thermal, motion: next_motion,
            report: SupportThermalReport { coupled, external_heat_j: external_heat,
                total_energy_residual_j: residual } })
    }
}
