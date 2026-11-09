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
    CircuitEnthalpySolve, DrivenCoefficients, ThermalDriver, coupled_driver_candidate, finite,
    contacts::{self, ContactHeatCells, ContactThermalCoefficients,
        material::{ElectricalContactNetwork, ElectricalContactNetworkError},
        heat::{HeatAssembly, HeatNetworkError, HeatNetworkSample, PreparedHeatNetwork}},
    material::{MaterialNetwork, MaterialNetworkError, MaterialNetworkSample},
};
use crate::material_circuit::{MaterialCircuitError, MaterialResistor};
use crate::material_circuit::conductivity::{ConductionCoordinate, ConductivityIteration, ConductorResponse};
use crate::material_circuit::conductivity::contact::{ContactCoordinate, ContactResponse};
use crate::material_circuit::contacts::SeriesElectricalContact;
use crate::material_circuit::motion::{FieldMotionCandidate, FieldMotionInput, FieldMotionStep,
    MaterialFieldMotion, conductivity::{FieldConductivityCandidate, FieldConductivityStep}};
use crate::material_circuit::storage::MaterialStorageElement;
use fs_exec::Cx;
use fs_matdb::EvaluationDecision;
use fs_material::{
    conductor::ElectricalPropertyBasis,
    field_storage::{ELECTRIC_PERMITTIVITY_PROPERTY, MAGNETIC_PERMEABILITY_PROPERTY},
    section_temperature::thermoelastic::{ConcentricThermoelasticLaw, ConcentricThermoelasticState},
    section::thermal::ConcentricThermoelasticSection,
    state_point::MaterialStatePointError,
};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
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
pub enum SupportThermalError<E = HeatNetworkError> {
    Network(MaterialNetworkError),
    Contact(ElectricalContactNetworkError),
    Material(MaterialStatePointError),
    /// Original material/transport provider failure, without losing its type.
    Heat(E),
    Invalid(&'static str),
}
impl<E: core::fmt::Debug> core::fmt::Display for SupportThermalError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "thermoelastic actuator support: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for SupportThermalError<E> {}
type Error<E = HeatNetworkError> = CircuitEnthalpyError<SupportThermalError<E>>;

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
pub struct SupportThermalSource<S = MaterialNetworkSample> {
    pub network: S,
    /// Original ordered interface heat partitions, separate from reversible heat.
    pub contact_destinations: Vec<ContactHeatCells>,
    pub layers: Vec<SupportLayerSample>,
    pub support_internal_energy_j: [f64; 2],
    pub reversible_heat_j: f64,
}
pub struct SupportThermalReport<S = MaterialNetworkSample, D = FieldMotionStep> {
    /// Original descriptor and finite enthalpy records. The native elastic
    /// starting energy uses the trial endpoint temperatures at the INITIAL gap.
    /// Its non_electrical_heat_j includes the reversible heat supplied above.
    pub coupled: CircuitEnthalpyReport<SupportThermalSource<S>, D>,
    /// Authored boundaries and sources, excluding internal Joule/reversible heat
    /// and the original closed enclosure's numerical net radiation exchange.
    pub external_heat_j: f64,
    /// Independent delta(U_caloric + U_support + H_field + K) minus supplies.
    pub total_energy_residual_j: f64,
}
pub struct SupportThermalCandidate<'m, S = MaterialNetworkSample, D = FieldMotionStep> {
    thermal: DepositionSimulation<'m>,
    motion: MaterialFieldMotion,
    report: SupportThermalReport<S, D>,
}
impl<'m, S, D> SupportThermalCandidate<'m, S, D> {
    pub const fn thermal(&self) -> &DepositionSimulation<'m> { &self.thermal }
    pub const fn motion(&self) -> &MaterialFieldMotion { &self.motion }
    pub const fn report(&self) -> &SupportThermalReport<S, D> { &self.report }
    /// Publish all histories together. Both internal heat transfers are present.
    pub fn into_parts(self) -> (DepositionSimulation<'m>, MaterialFieldMotion, SupportThermalReport<S, D>) {
        (self.thermal, self.motion, self.report)
    }
}

// Borrow each original driver's physical trajectory without rebuilding its
// candidate or discarding its constitutive samples and native work count.
trait SupportDriver: ThermalDriver<Model = MaterialFieldMotion> {
    fn model(&self) -> &MaterialFieldMotion;
    fn step(&self) -> &Self::Step;
    fn field(step: &Self::Step) -> &FieldMotionStep;
    fn circuit_trials(&self) -> usize;
}
impl SupportDriver for FieldMotionCandidate {
    fn model(&self) -> &MaterialFieldMotion { FieldMotionCandidate::model(self) }
    fn step(&self) -> &Self::Step { FieldMotionCandidate::step(self) }
    fn field(step: &Self::Step) -> &FieldMotionStep { step }
    fn circuit_trials(&self) -> usize { 1 }
}
impl SupportDriver for FieldConductivityCandidate {
    fn model(&self) -> &MaterialFieldMotion { self.motion() }
    fn step(&self) -> &Self::Step { FieldConductivityCandidate::step(self) }
    fn field(step: &Self::Step) -> &FieldMotionStep { &step.field }
    fn circuit_trials(&self) -> usize { self.step().iterations }
}

fn invalid<E>(message: &'static str) -> Error<E> {
    CircuitEnthalpyError::Source(SupportThermalError::Invalid(message))
}
fn material<E>(error: MaterialStatePointError) -> Error<E> {
    CircuitEnthalpyError::Source(SupportThermalError::Material(error))
}
fn network<E>(error: MaterialNetworkError) -> Error<E> {
    match error {
        MaterialNetworkError::Thermal(DepositionError::Cancelled) => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(SupportThermalError::Network(other)),
    }
}
fn contact<E>(error: ElectricalContactNetworkError) -> Error<E> {
    match error {
        ElectricalContactNetworkError::Thermal(DepositionError::Cancelled)
        | ElectricalContactNetworkError::Bulk(MaterialNetworkError::Thermal(DepositionError::Cancelled))
            => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(SupportThermalError::Contact(other)),
    }
}
fn check<E>(cancel: &RefCell<impl FnMut() -> bool>) -> Result<(), Error<E>> {
    if (*cancel.borrow_mut())() { Err(CircuitEnthalpyError::Cancelled) } else { Ok(()) }
}
fn temperatures<E>(thermal: &DepositionSimulation<'_>, cells: &[CellId]) -> Result<Vec<f64>, Error<E>> {
    cells.iter().map(|&id| thermal.temperature_k(id).map_err(Into::into)).collect()
}
fn strain<E>(motion: &MaterialFieldMotion) -> Result<f64, Error<E>> {
    let support = motion.support().ok_or_else(|| invalid("an existing axial member is required"))?;
    Ok(finite((motion.gap_m() - support.reference_gap_m()) / support.length_m())?)
}

fn admit<E>(
    material_network: &MaterialNetwork<'_>, thermal: &DepositionSimulation<'_>,
    motion: &MaterialFieldMotion, binding: &SupportThermalBinding<'_>,
    cancel: &RefCell<impl FnMut() -> bool>,
) -> Result<ConcentricThermoelasticState, Error<E>> {
    check(cancel)?;
    let support = motion.support().ok_or_else(|| invalid("an existing axial member is required"))?;
    if !support.is_thermally_live()
        || binding.layer_cells.len() != support.section().section().layers().len()
        || !binding.volume_tolerance_m3.is_finite() || binding.volume_tolerance_m3 < 0.0
        || binding.layer_cells.iter().copied().collect::<BTreeSet<_>>().len() != binding.layer_cells.len()
    {
        return Err(invalid("live support, unique complete layer cells and finite geometry tolerance required"));
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
        for (index, pair) in knots.windows(2).enumerate() {
            check(cancel)?;
            let dt = finite(pair[1].temperature_k - pair[0].temperature_k)?;
            if dt <= 0.0 {
                return Err(invalid("live support caloric temperature must increase strictly"));
            }
            let capacity = finite(cell.mass_kg *
                source.curve.sensible_heat_capacity_bounds_j_kg_k(index)
                    .map_err(|_| invalid("support caloric heat capacity is unrepresentable"))?[0])?;
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

// The electrical owner may vary E/J and individual interface coordinates. The
// existing support potential and fixed caloric chart do not own those axes.
// Check consumed source claims, including applicability, rather than accepting
// an initially valid zero-field slice of a different thermodynamic model.
fn admit_conductivity_support<E>(
    network: &MaterialNetwork<'_>, binding: &SupportThermalBinding<'_>,
    responses: &[ConductorResponse<'_>], contact_responses: Option<&[ContactResponse<'_>]>,
    routes: &[ContactHeatCells], cancel: &RefCell<impl FnMut() -> bool>,
) -> Result<(), Error<E>> {
    let electrical = |dims| [ConductionCoordinate::ElectricField.dimensions(),
        ConductionCoordinate::CurrentDensity.dimensions(), ContactCoordinate::Voltage.dimensions(),
        ContactCoordinate::Current.dimensions()].contains(&dims);
    let reference = binding.law.reference();
    for (index, layer) in reference.layers().iter().enumerate() {
        check(cancel)?;
        let id = binding.layer_cells.get(index)
            .ok_or_else(|| invalid("support layer lacks its original caloric cell"))?;
        let cell = network.cells.iter().find(|cell| &cell.cell == id)
            .ok_or_else(|| invalid("support layer lacks its original material source"))?;
        let mut properties = reference.elastic_source_properties(index).map_err(material)?;
        properties.push(layer.material.property(
            fs_material::state_point::LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY)
            .ok_or_else(|| invalid("support lost its original expansion source"))?);
        for property in properties {
            let receipt = &property.answer().receipt;
            let claim = cell.point.material.claims().claim(receipt.selected)
                .ok_or_else(|| invalid("original support source claim is absent"))?;
            let field_axis = |axis: &str| [claim.validity.axis_quantities().get(axis),
                receipt.axis_quantities.get(axis), cell.point.context.axis_quantities().get(axis)]
                .into_iter().flatten().any(|quantity| electrical(quantity.dims()));
            if claim.value.axes().any(|(axis, dims)| electrical(dims) || field_axis(axis))
                || claim.validity.bounds().keys().any(|axis| field_axis(axis)) {
                return Err(invalid("support mechanical response and applicability must be independent of electrical coordinates"));
            }
        }
    }
    for cell in network.cells {
        for receipt in cell.curve.source_receipts() {
            check(cancel)?;
            let claim = cell.point.material.claims().claim(receipt.selected)
                .ok_or_else(|| invalid("original caloric source claim is absent"))?;
            let field_axis = |axis: &str| {
                responses.iter().any(|response| response.axes.iter().any(|a| a.name == axis)
                    && network.resistors.iter().any(|resistor|
                        resistor.destination.cell == cell.cell && resistor.destination.branch == response.branch))
                    || contact_responses.into_iter().flatten().any(|response|
                        response.axes.iter().any(|a| a.name == axis)
                        && routes.iter().any(|route| route.id == response.id && route.cells.contains(&cell.cell)))
                    || [claim.validity.axis_quantities().get(axis), receipt.axis_quantities.get(axis),
                        cell.point.context.axis_quantities().get(axis)].into_iter().flatten()
                        .any(|quantity| electrical(quantity.dims()))
            };
            if claim.value.axes().any(|(axis, dims)| electrical(dims) || field_axis(axis))
                || claim.validity.bounds().keys().any(|axis| field_axis(axis)) {
                return Err(invalid("support caloric response and applicability must be independent of electrical coordinates"));
            }
        }
    }
    Ok(())
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
        self.resolved_field_motion_support_candidate(thermal, target_s, motion, ports,
            &[], binding, controls, guess_j, |trial| {
                Ok(ContactThermalCoefficients {
                    bulk: prepared.resolve(trial, || (*cancel.borrow_mut())()).map_err(network)?,
                    contacts: Vec::new(),
                })
            }, |section, seeds, contacts, _| {
                Ok(motion.step_thermal_support_candidate(section, seeds, contacts,
                    FieldMotionInput { electrical_ports: ports, external_force_n,
                        duration_s: target_s - thermal.time_s() }, controls.circuit,
                    || (*cancel.borrow_mut())())?)
            }, |_, _, transport| Ok(transport), |_| 0.0, || (*cancel.borrow_mut())())
    }

    // All public material frontends use one thermoelastic residual and one
    // native field/body trajectory. Only the existing source resolver differs.
    #[allow(clippy::too_many_arguments)]
    fn resolved_field_motion_support_candidate<'m, S, E, D: SupportDriver>(
        &self, thermal: &DepositionSimulation<'m>, target_s: f64,
        motion: &MaterialFieldMotion, ports: &[f64],
        routes: &[ContactHeatCells], binding: &SupportThermalBinding<'_>,
        controls: CircuitEnthalpySolve, guess_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'m>)
            -> Result<ContactThermalCoefficients<S>, Error<E>>,
        advance: impl Fn(ConcentricThermoelasticSection, Vec<MaterialResistor>,
            Vec<SeriesElectricalContact>, usize) -> Result<D, Error<E>>,
        transport: impl Fn(&DepositionSimulation<'m>, &D::Step, StepCoefficients<S>)
            -> Result<StepCoefficients<S>, Error<E>>,
        internal_heat_w: impl Fn(&S) -> f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<SupportThermalCandidate<'m, S, D::Step>, Error<E>> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        let routes = contacts::admit_routes(thermal, motion.circuit(), routes, || {
            if (*cancel.borrow_mut())() { Err(DepositionError::Cancelled) } else { Ok(()) }
        })?.into_values().collect::<Vec<_>>();
        let initial = admit(self, thermal, motion, binding, &cancel)?;
        let support = motion.support().expect("admitted support");
        let length = support.length_m();
        let initial_strain = strain(motion)?;
        let destinations = self.resistors.iter().map(|r| r.destination).collect::<Vec<_>>();
        let input = CircuitEnthalpyInput { circuit: motion.circuit(), ports, destinations: &destinations };
        let dt = target_s - thermal.time_s();
        let (next_thermal, next_motion, coupled) = coupled_driver_candidate(
            thermal, target_s, input, &routes, controls, guess_j,
            |trial, remaining| {
                check(&cancel)?;
                let resolved = resolve(trial)?;
                contacts::admit_sources(motion.circuit(), &resolved.contacts,
                    || (*cancel.borrow_mut())())?;
                let coefficients = resolved.bulk;
                let temperatures = temperatures(trial, &binding.layer_cells)?;
                let frozen = binding.law.at_temperatures_k(&temperatures, initial_strain,
                    support.linear_strain_limit()).map_err(material)?;
                let native = advance(frozen.section, coefficients.resistors,
                    resolved.contacts, remaining)?;
                let mut coefficients = transport(trial, native.step(), StepCoefficients {
                    sources_w: coefficients.sources_w, conductances_w_k: coefficients.conductances_w_k,
                    payload: coefficients.source })?;
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
                if free.iter().zip(D::field(native.step()).elastic_energy_j).any(|(&a, b)|
                    (a - b).abs() > controls.total_energy_tolerance_j)
                {
                    return Err(invalid("native support potential and shared per-layer Helmholtz energy disagree"));
                }
                let circuit_trials = native.circuit_trials();
                Ok(DrivenCoefficients { circuit: native, sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k, circuit_trials,
                    source: SupportThermalSource { network: coefficients.payload,
                        contact_destinations: routes.clone(), layers,
                        support_internal_energy_j: internal, reversible_heat_j: reversible_heat } })
            }, || (*cancel.borrow_mut())())?;
        let internal_heat = finite(dt * internal_heat_w(&coupled.source.network))?;
        let external_heat = finite(coupled.non_electrical_heat_j
            - coupled.source.reversible_heat_j - internal_heat)?;
        let step = D::field(&coupled.electrical);
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

impl ElectricalContactNetwork<'_> {
    /// Advance original electrical contacts, both finite heat receivers, and a
    /// thermally live actuator support in one shared physical transaction.
    /// Pinned interface resistance responds to actual endpoint temperatures;
    /// its declared A/B partition enters the existing enthalpy residual once.
    /// Layer expansion and reversible heat use the same support owner as the
    /// contact-free path. Source refusal, budget failure and cancellation retain
    /// all accepted electrical, mechanical and caloric histories.
    ///
    /// Interface area and partition remain fixed over the interval. This adds
    /// no contact mechanics, opening, arcs, interface storage or strain-dependent
    /// electrical geometry. All original support/field restrictions apply.
    #[allow(clippy::too_many_arguments)]
    pub fn field_motion_support_candidate<'m>(
        &self, thermal: &DepositionSimulation<'m>, target_s: f64,
        motion: &MaterialFieldMotion, ports: &[f64], external_force_n: f64,
        binding: &SupportThermalBinding<'_>, controls: CircuitEnthalpySolve,
        guess_j: &[f64], cancelled: impl FnMut() -> bool,
    ) -> Result<SupportThermalCandidate<'m>, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        super::super::admit_driver_unknowns(thermal, controls, guess_j)?;
        let routes = self.heat_destinations(thermal).map_err(contact)?;
        let prepared = self.prepare(thermal, || (*cancel.borrow_mut())()).map_err(contact)?;
        self.bulk.resolved_field_motion_support_candidate(thermal, target_s, motion, ports,
            &routes, binding, controls, guess_j,
            |trial| prepared.resolve(trial, || (*cancel.borrow_mut())()).map_err(contact),
            |section, seeds, contacts, _| {
                Ok(motion.step_thermal_support_candidate(section, seeds, contacts,
                    FieldMotionInput { electrical_ports: ports, external_force_n,
                        duration_s: target_s - thermal.time_s() }, controls.circuit,
                    || (*cancel.borrow_mut())())?)
            }, |_, _, transport| Ok(transport), |_| 0.0, || (*cancel.borrow_mut())())
    }
}

/// Close source-dependent bulk E/J, individual contact V/I/J and thermally live
/// support motion on the original finite-enthalpy residual. The existing native
/// conductivity owner supplies every material iteration and its actual loss
/// account; the existing support potential supplies reversible heat once.
///
/// `network` retains the original support/caloric inventory and resistor routes.
/// `resolve` must refresh the original pinned bulk/contact bundles at the trial
/// temperatures. `transport` then receives their ACTUAL solved electrical point
/// and refreshes any E/J-dependent thermal coefficients and observations before
/// heat acceptance. It supplies no Joule or reversible support heat and advances
/// no history. `internal_heat_w` is the original closed-enclosure numerical net,
/// or zero when no such internal source is present; it is never balancing heat.
///
/// All native electrical trials, including discarded thermal stencils and final
/// admission, consume the outer circuit budget. Narrowing that budget retains
/// fixed/adaptive source-coordinate relaxation. Applied/evaluated material and
/// contact residuals remain inspectable in the native electrical step report.
/// Source refusal, nonconvergence, exhausted work or cancellation publishes none
/// of the caloric, electrical, gap, momentum or support histories.
///
/// Support and caloric claims must remain independent of electrical coordinates,
/// including their applicability. Original fixed-inventory, stationary-conductor,
/// constant-field-storage and first-order thermal timing restrictions still apply.
#[allow(clippy::too_many_arguments)]
pub fn step_conductivity_candidate<'m, S, E>(
    network: &MaterialNetwork<'_>, thermal: &DepositionSimulation<'m>, target_s: f64,
    motion: &MaterialFieldMotion, ports: &[f64], external_force_n: f64,
    responses: &[ConductorResponse<'_>], contact_responses: Option<&[ContactResponse<'_>]>,
    routes: &[ContactHeatCells], binding: &SupportThermalBinding<'_>,
    controls: CircuitEnthalpySolve, conductivity: impl Into<ConductivityIteration>, guess_j: &[f64],
    resolve: impl Fn(&DepositionSimulation<'m>) -> Result<ContactThermalCoefficients<S>, E>,
    transport: impl Fn(&DepositionSimulation<'m>, &FieldConductivityStep, StepCoefficients<S>)
        -> Result<StepCoefficients<S>, E>,
    internal_heat_w: impl Fn(&S) -> f64, cancelled: impl FnMut() -> bool,
) -> Result<SupportThermalCandidate<'m, S, FieldConductivityStep>, Error<E>> {
    let cancel = RefCell::new(cancelled);
    check(&cancel)?;
    super::super::admit_driver_unknowns(thermal, controls, guess_j)?;
    let conductivity = conductivity.into();
    conductivity.validate()?;
    admit_conductivity_support(network, binding, responses, contact_responses, routes, &cancel)?;
    network.resolved_field_motion_support_candidate(thermal, target_s, motion, ports,
        routes, binding, controls, guess_j,
        |trial| resolve(trial).map_err(|error| CircuitEnthalpyError::Source(SupportThermalError::Heat(error))),
        |section, seeds, contacts, remaining| {
            let bounded = conductivity.bounded(remaining)?;
            motion.conductivity_thermal_support_candidate(section, seeds, contacts,
                responses, contact_responses, FieldMotionInput { electrical_ports: ports,
                    external_force_n, duration_s: target_s - thermal.time_s() },
                controls.circuit, bounded, || (*cancel.borrow_mut())()).map_err(|error| match error {
                    MaterialCircuitError::ConductivityConvergence { .. }
                        if bounded.controls().maximum_iterations < conductivity.controls().maximum_iterations =>
                            CircuitEnthalpyError::CircuitBudget,
                    other => CircuitEnthalpyError::from(other),
                })
        }, |state, electrical, coefficients| transport(state, electrical, coefficients)
            .map_err(|error| CircuitEnthalpyError::Source(SupportThermalError::Heat(error))),
        internal_heat_w, || (*cancel.borrow_mut())())
}

/// Advance a thermally live actuator through its complete original heat assembly.
/// Native thermal contacts, enclosure radiation, source-integrated/directional
/// links and composed coolant/reservoir laws all resolve at the SAME proposed
/// layer temperatures as electrical resistance and support expansion. The
/// original prepared assembly retains every source pin and both electrical
/// heat routes; its complete sample is retained in the accepted support report.
///
/// This uses the same thermoelastic residual, per-layer reversible heat and
/// reciprocal field/body trajectory as the plain material entrypoints. Closed
/// enclosure exchange remains internal, including its numerical net; physical
/// reservoir exchange remains external. All heat enters the existing finite
/// cells once. Source, budget and cancellation failures publish no histories.
///
/// Thermal link and radiating geometry remain authored and fixed over this
/// interval. Support displacement does not infer changing contact area,
/// pressure, view factors or moving coolant volume. All original live-support
/// material, inertia, field-storage and endpoint-timing restrictions still apply.
#[allow(clippy::too_many_arguments)]
pub fn step_heat_assembly_candidate<'m, A: HeatAssembly + ?Sized>(
    assembly: &A, cx: &Cx<'_>, thermal: &DepositionSimulation<'m>, target_s: f64,
    motion: &MaterialFieldMotion, ports: &[f64], external_force_n: f64,
    binding: &SupportThermalBinding<'_>, controls: CircuitEnthalpySolve,
    guess_j: &[f64],
) -> Result<SupportThermalCandidate<'m, A::Sample>, Error<A::Error>> {
    cx.checkpoint().map_err(|_| CircuitEnthalpyError::Cancelled)?;
    super::super::admit_driver_unknowns(thermal, controls, guess_j)?;
    let prepared = assembly.prepared_heat_network();
    let routes = prepared.heat_destinations(thermal).map_err(|error| match error {
        HeatNetworkError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(SupportThermalError::Heat(other.into())),
    })?;
    prepared.field_storage_material_network().resolved_field_motion_support_candidate(
        thermal, target_s, motion, ports, &routes, binding, controls,
        guess_j, |trial| assembly.resolve_heat(cx, trial)
            .map_err(|error| error.map_source(SupportThermalError::Heat)),
        |section, seeds, contacts, _| {
            Ok(motion.step_thermal_support_candidate(section, seeds, contacts,
                FieldMotionInput { electrical_ports: ports, external_force_n,
                    duration_s: target_s - thermal.time_s() }, controls.circuit,
                || cx.checkpoint().is_err())?)
        }, |_, _, transport| Ok(transport),
        |sample| sample.as_ref().radiation_net_power_w(), || cx.checkpoint().is_err())
}

impl PreparedHeatNetwork<'_, '_> {
    /// Continue this source-pinned heat assembly with the original live support.
    /// See [`step_heat_assembly_candidate`] for transport and work ownership.
    #[allow(clippy::too_many_arguments)]
    pub fn field_motion_support_candidate<'m>(
        &self, cx: &Cx<'_>, thermal: &DepositionSimulation<'m>, target_s: f64,
        motion: &MaterialFieldMotion, ports: &[f64], external_force_n: f64,
        binding: &SupportThermalBinding<'_>, controls: CircuitEnthalpySolve,
        guess_j: &[f64],
    ) -> Result<SupportThermalCandidate<'m, HeatNetworkSample>, Error> {
        step_heat_assembly_candidate(self, cx, thermal, target_s, motion, ports,
            external_force_n, binding, controls, guess_j)
    }
}
