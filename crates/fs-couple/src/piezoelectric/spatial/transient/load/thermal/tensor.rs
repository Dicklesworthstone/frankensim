//! Temperature-dependent full-tensor piezoelectric motion with a finite caloric body.
//!
//! One spatial body has one uniform ceramic temperature and one complete caloric
//! inventory. Original per-cell tensors and nonuniform strain/electric fields
//! remain resolved. Each thermal trial reassembles the native operators, solves
//! electrostatics at fixed displacement and free charge, then advances the SAME
//! consistent-mass midpoint system. No mechanical coordinate is reset.
//!
//! With psi the full field Helmholtz energy, the signed reversible transfer is
//! Q_rev = -[psi(T_next,z_old)-psi(T_old,z_old)]-delta(T*S).
//! It enters the original enthalpy residual separately from physical Joule heat.
//! Uniform ceramic temperature makes redistribution of the resolved electric
//! field internal to this one thermal inventory. Independent cell temperatures
//! require a local field-energy transfer model and are not inferred here.
//!
//! Geometry and density stay fixed. Ground electrodes remain at zero; other
//! electrodes retain charge/current control, optionally with physical resistors.
//! Changing ideal voltage constraints, intrinsic memory, thermal expansion,
//! pyroelectric offsets and phase transitions require their own accounted laws.
//! Endpoint-temperature coefficients make this generally first order in time.

use super::{
    LoadThermalError, SpatialPiezoelectricThermalStep, ThermalCeramicCell,
    ThermalElectrodeLoad, bind_loads,
};
use crate::circuit_enthalpy::{
    CircuitEnthalpyError, CircuitEnthalpyReport, CircuitEnthalpySolve, DrivenCoefficients,
    DriverAccount, ThermalDriver, admit_driver_unknowns, cell_driver_candidate,
    material::{MaterialNetworkError, MaterialNetworkSample, PreparedMaterialNetwork},
};
use crate::piezoelectric::spatial::{
    ElectrodeCondition, PiezoelectricElectrode, TetPiezoelectricControls,
    TetPiezoelectricError, TetPiezoelectricModel,
    transient::{FixedBoundary, TetPiezoelectricTransient, TetPiezoelectricTransientState},
};
use crate::piezoelectric::spatial::transient::load::ElectrodeDrive;
use fs_exec::Cx;
use fs_material::{
    field_storage::thermal::FieldThermalState,
    piezoelectric::{
        PiezoelectricError,
        tensor::TensorPiezoelectric,
        thermal::{ResolvedThermalTensorPiezoelectric, ThermalTensorPiezoelectric},
    },
};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation, ThermalLaw};
use std::{collections::BTreeMap, sync::Arc};

type Error = CircuitEnthalpyError<LoadThermalError>;
type Result<T> = core::result::Result<T, Error>;

fn invalid(message: &'static str) -> Error { CircuitEnthalpyError::Invalid(message) }
fn field(error: PiezoelectricError) -> Error {
    spatial(TetPiezoelectricError::Material(error))
}
fn spatial(error: TetPiezoelectricError) -> Error {
    match error {
        TetPiezoelectricError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(LoadThermalError::Piezoelectric(other)),
    }
}
fn network(error: MaterialNetworkError) -> Error {
    match error {
        MaterialNetworkError::Thermal(DepositionError::Cancelled) => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(LoadThermalError::Material(other)),
    }
}
fn check(cx: &Cx<'_>) -> Result<()> {
    cx.checkpoint().map_err(|_| CircuitEnthalpyError::Cancelled)
}
fn finite(value: f64) -> Result<f64> {
    if value.is_finite() { Ok(value) } else { Err(invalid("tensor thermal account is unrepresentable")) }
}
fn norm(values: &[f64]) -> f64 {
    values.iter().fold(0.0_f64, |value, x| value.hypot(*x))
}
fn same_inventory(actual: f64, expected: f64) -> bool {
    actual.is_finite() && actual > 0.0 && expected.is_finite() && expected > 0.0
        && (actual / expected - 1.0).abs() <= 1e-10
}

/// Complete field/kinetic history with immutable geometry and thermal source laws.
/// Accepted material tensors are retained separately from the original model.
/// Native operator views borrow them only during a trial; no self-reference,
/// alternative field solver or hidden permanent-temperature reservoir is used.
#[derive(Clone)]
pub struct ThermalTensorPiezoelectricBody<'model, 'mesh> {
    reference: &'model TetPiezoelectricModel<'mesh>,
    accepted_model: Option<Arc<TetPiezoelectricModel<'mesh>>>,
    boundary: Arc<FixedBoundary>,
    state: TetPiezoelectricTransientState,
    laws: Arc<[ThermalTensorPiezoelectric]>,
    caloric_cell: CellId,
    maximum_displacement_norm_c_m2: f64,
}

/// Actual endpoint field thermodynamics and the original material-network sample.
#[derive(Clone, Debug)]
pub struct TensorThermalSource {
    pub network: MaterialNetworkSample,
    pub caloric_cell: CellId,
    pub temperature_k: [f64; 2],
    /// One pair of integrated local potentials per original tetrahedron.
    pub cells: Vec<[FieldThermalState; 2]>,
    /// Whole ceramic field U=psi+T*S, excluding consistent kinetic energy.
    pub field_internal_energy_j: [f64; 2],
    /// Frozen-displacement/free-charge Helmholtz change, not external work.
    pub parameter_work_j: f64,
    /// Signed transfer into the single finite ceramic caloric cell, not loss.
    pub reversible_heat_j: f64,
}

/// Native midpoint/Newton reports plus the physical U-based combined audit.
pub struct TensorThermalReport {
    /// The common driver's isothermal psi account is retained in this subreport.
    pub coupled: CircuitEnthalpyReport<TensorThermalSource, SpatialPiezoelectricThermalStep>,
    /// Authored boundary/source heat, excluding resistor and reversible transfer.
    pub external_heat_j: f64,
    /// Delta(kinetic+field U+caloric U) minus mechanical/electrical/heat supply.
    pub total_energy_residual_j: f64,
}

/// Transactional field, source-tensor and thermal endpoint.
pub struct TensorThermalCandidate<'thermal, 'model, 'mesh> {
    thermal: DepositionSimulation<'thermal>,
    body: ThermalTensorPiezoelectricBody<'model, 'mesh>,
    report: TensorThermalReport,
}
impl<'thermal, 'model, 'mesh> TensorThermalCandidate<'thermal, 'model, 'mesh> {
    pub const fn thermal(&self) -> &DepositionSimulation<'thermal> { &self.thermal }
    pub const fn body(&self) -> &ThermalTensorPiezoelectricBody<'model, 'mesh> { &self.body }
    pub const fn report(&self) -> &TensorThermalReport { &self.report }
    pub fn into_parts(self) -> (
        DepositionSimulation<'thermal>,
        ThermalTensorPiezoelectricBody<'model, 'mesh>,
        TensorThermalReport,
    ) { (self.thermal, self.body, self.report) }
}

struct TensorDriver<'model, 'mesh> {
    body: ThermalTensorPiezoelectricBody<'model, 'mesh>,
    step: SpatialPiezoelectricThermalStep,
}
impl<'model, 'mesh> ThermalDriver for TensorDriver<'model, 'mesh> {
    type Model = ThermalTensorPiezoelectricBody<'model, 'mesh>;
    type Step = SpatialPiezoelectricThermalStep;
    fn account(&self) -> DriverAccount<'_> {
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: self.step.loads.joule_heat_j,
            numerical_dissipation_j: 0.0,
            stored_change_j: self.step.native.kinetic_energy_change_j
                + self.step.native.stored_energy_change_j,
            supplied_j: self.step.native.mechanical_work_j
                + self.step.loads.electrical_source_work_j,
        }
    }
    fn into_parts(self) -> (Self::Model, Self::Step) { (self.body, self.step) }
}

impl<'model, 'mesh> ThermalTensorPiezoelectricBody<'model, 'mesh> {
    /// Attach one thermal law to each original material slot (uniform or per
    /// tetrahedron), preserving the accepted body. All tensors share the one
    /// finite caloric cell and temperature. The displacement bound is a physical
    /// capacity-admission envelope, enforced at every evaluated field endpoint.
    pub fn new(
        body: &TetPiezoelectricTransient<'model, 'mesh>,
        laws: Vec<ThermalTensorPiezoelectric>,
        caloric_cell: CellId,
        maximum_displacement_norm_c_m2: f64,
        cx: &Cx<'_>,
    ) -> Result<Self> {
        check(cx)?;
        if body.intrinsic.is_some() || laws.len() != body.model.materials.len()
            || laws.is_empty() || !maximum_displacement_norm_c_m2.is_finite()
            || maximum_displacement_norm_c_m2 < 0.0
            || body.state.electrode_conditions.iter().any(|condition|
                matches!(condition, ElectrodeCondition::Voltage(v) if *v != 0.0))
        {
            return Err(invalid(
                "complete tensor thermal laws, fixed ground/charge electrodes and no intrinsic memories required",
            ));
        }
        let mut temperature = None;
        for (law, source) in laws.iter().zip(body.model.materials()) {
            check(cx)?;
            let at = law.at_material(source).map_err(field)?;
            if law.temperature_axis() != laws[0].temperature_axis()
                || temperature.is_some_and(|t: f64| t.to_bits() != at.temperature_k().to_bits())
            {
                return Err(invalid("all ceramic sources must share the same actual temperature and axis"));
            }
            temperature = Some(at.temperature_k());
        }
        if body.state.cells().iter().any(|cell|
            norm(&cell.response.displacement_c_m2) > maximum_displacement_norm_c_m2)
        {
            return Err(invalid("initial electric displacement exceeds the thermal capacity envelope"));
        }
        Ok(Self {
            reference: body.model,
            accepted_model: None,
            boundary: Arc::clone(&body.boundary),
            state: body.state.clone(),
            laws: laws.into(),
            caloric_cell,
            maximum_displacement_norm_c_m2,
        })
    }

    /// Original node/cell order, accepted fields, velocities, charge and clock.
    pub const fn state(&self) -> &TetPiezoelectricTransientState { &self.state }
    /// Actual accepted source tensors, uniform or in original cell order.
    pub fn materials(&self) -> &[TensorPiezoelectric] { self.model().materials() }
    fn model(&self) -> &TetPiezoelectricModel<'mesh> {
        self.accepted_model.as_deref().unwrap_or(self.reference)
    }
    pub const fn caloric_cell(&self) -> CellId { self.caloric_cell }
    pub const fn maximum_displacement_norm_c_m2(&self) -> f64 {
        self.maximum_displacement_norm_c_m2
    }

    fn temperature(
        &self, thermal: &DepositionSimulation<'_>, material: &PreparedMaterialNetwork<'_, '_>,
    ) -> Result<f64> {
        let binding = material.cell_bindings().iter().find(|entry| entry.cell == self.caloric_cell)
            .ok_or_else(|| invalid("ceramic caloric binding is absent"))?;
        let cell = thermal.cells().get(&self.caloric_cell)
            .ok_or_else(|| invalid("ceramic caloric cell is inactive"))?;
        Ok(binding.curve.state_at_specific_enthalpy(cell.specific_h_j_kg)
            .map_err(|error| network(MaterialNetworkError::Phase(error)))?.temperature_k())
    }

    // Existing source/curve identities and whole inventory remain authoritative.
    fn admit(
        &self, thermal: &DepositionSimulation<'_>, material: &PreparedMaterialNetwork<'_, '_>,
        controls: TetPiezoelectricControls, cx: &Cx<'_>,
    ) -> Result<Vec<FieldThermalState>> {
        check(cx)?;
        let binding = material.cell_bindings().iter().find(|entry| entry.cell == self.caloric_cell)
            .ok_or_else(|| invalid("ceramic caloric binding is absent"))?;
        let cell = thermal.cells().get(&self.caloric_cell)
            .ok_or_else(|| invalid("ceramic caloric cell is inactive"))?;
        let curve = binding.curve;
        if cell.material_id != binding.material_id
            || !thermal.uses_material_law(cell.material_id, curve as &dyn ThermalLaw)
        {
            return Err(invalid("ceramic thermal cell must retain its actual borrowed caloric law"));
        }
        let phase = curve.state_at_specific_enthalpy(cell.specific_h_j_kg)
            .map_err(|error| network(MaterialNetworkError::Phase(error)))?;
        let query = binding.point.at(phase)
            .map_err(|error| network(MaterialNetworkError::Material(error)))?;
        let mut range = [0.0_f64, f64::INFINITY];
        let mut samples = Vec::with_capacity(self.laws.len());
        for (law, source) in self.laws.iter().zip(self.materials()) {
            check(cx)?;
            let at = law.at_material(source).map_err(field)?;
            if binding.point.temperature_axis != law.temperature_axis()
                || binding.point.material.content_hash() != law.card().content_hash()
                || curve.material_card_identity() != law.card().content_hash()
                || query != source.material().query()
                || at.temperature_k().to_bits() != phase.temperature_k().to_bits()
            {
                return Err(invalid("ceramic field and caloric owners must share card, complete context and temperature"));
            }
            let density = source.material().property("density")
                .ok_or_else(|| invalid("ceramic density source is missing"))?
                .answer().receipt.selected;
            let receipts = curve.source_receipts();
            if receipts.is_empty() || receipts.len() % 2 != 0
                || receipts.chunks_exact(2).any(|pair|
                    pair[1].property != "density" || pair[1].selected != density)
                || receipts.iter().any(|receipt|
                    receipt.axis_quantities != *source.material().axis_quantities()
                    || !receipt.query_point.iter().filter(|(axis, _)| axis != law.temperature_axis())
                        .eq(source.material().query_point().iter()
                            .filter(|(axis, _)| axis != law.temperature_axis())))
                || curve.knots().iter().any(|knot|
                    knot.liquid_mass_fraction != 0.0
                    || knot.bulk_density_kg_m3.to_bits() != source.density_kg_m3().to_bits())
            {
                return Err(invalid("ceramic caloric chart must preserve original density, context and solid phase"));
            }
            let supported = law.temperature_range_k();
            range[0] = range[0].max(supported[0]);
            range[1] = range[1].min(supported[1]);
            samples.push(at);
        }
        let knots = curve.knots();
        if range[0] >= range[1] || range[0] < knots[0].temperature_k
            || range[1] > knots[knots.len()-1].temperature_k
        {
            return Err(invalid("complete tensor temperature interval requires caloric support"));
        }
        let mut volumes = vec![0.0; self.laws.len()];
        let mut mass = 0.0;
        let mut volume = 0.0;
        for (index, &v) in self.reference.dielectric.volumes_m3().iter().enumerate() {
            check(cx)?;
            let slot = if self.laws.len() == 1 { 0 } else { index };
            volumes[slot] = finite(volumes[slot] + v)?;
            volume = finite(volume + v)?;
            mass = finite(mass + v * self.materials()[slot].density_kg_m3())?;
        }
        if !same_inventory(cell.geometry.volume_m3, volume) || !same_inventory(cell.mass_kg, mass) {
            return Err(invalid("ceramic caloric inventory must equal the whole original spatial volume and mass"));
        }
        let mut floor = 0.0;
        for (law, volume) in self.laws.iter().zip(volumes) {
            floor = finite(floor + law.minimum_heat_capacity_j_k(
                controls.maximum_strain_norm, self.maximum_displacement_norm_c_m2, volume,
            ).map_err(field)?)?;
        }
        let mut caloric = f64::INFINITY;
        for pair in knots.windows(2) {
            check(cx)?;
            let delta = pair[1].temperature_k - pair[0].temperature_k;
            if delta <= 0.0 { return Err(invalid("ceramic caloric chart needs increasing temperature")); }
            if pair[1].temperature_k >= range[0] && pair[0].temperature_k <= range[1] {
                caloric = caloric.min(finite(cell.mass_kg
                    * (pair[1].specific_enthalpy_j_kg-pair[0].specific_enthalpy_j_kg)/delta)?);
            }
        }
        if !caloric.is_finite() || finite(caloric + floor)? <= 0.0 {
            return Err(invalid("whole ceramic caloric plus field capacity must remain positive"));
        }
        self.field_states(&self.state, &samples, controls, cx)
    }

    fn field_states(
        &self, state: &TetPiezoelectricTransientState,
        samples: &[ResolvedThermalTensorPiezoelectric], controls: TetPiezoelectricControls,
        cx: &Cx<'_>,
    ) -> Result<Vec<FieldThermalState>> {
        let mut result = Vec::with_capacity(state.cells().len());
        let mut free = 0.0;
        for (i, cell) in state.cells().iter().enumerate() {
            check(cx)?;
            if norm(&cell.strain_mandel) > controls.maximum_strain_norm
                || norm(&cell.response.displacement_c_m2) > self.maximum_displacement_norm_c_m2
                || norm(&cell.response.electric_field_v_m) > controls.maximum_electric_field_v_m
            {
                return Err(invalid("ceramic field leaves its strain, displacement or electric-field envelope"));
            }
            let sample = &samples[if samples.len() == 1 { 0 } else { i }];
            let value = sample.state(cell.strain_mandel, cell.response.displacement_c_m2,
                self.reference.dielectric.volumes_m3()[i]).map_err(field)?;
            free = finite(free + value.free_energy_j)?;
            result.push(value);
        }
        if (free - state.stored_energy_j()).abs() > controls.energy_tolerance_j {
            return Err(invalid("native tensor storage and sourced thermal Helmholtz energy disagree"));
        }
        Ok(result)
    }

    // Only electrical algebraic coordinates refresh. A zero-velocity initializer
    // prepares the native field view; then the original velocity, consistent
    // kinetic energy, boundary topology and absolute clock are restored exactly.
    fn rebound<'trial>(
        &self, model: &'trial TetPiezoelectricModel<'mesh>,
        controls: TetPiezoelectricControls, cx: &Cx<'_>,
    ) -> Result<TetPiezoelectricTransient<'trial, 'mesh>> {
        let mut supports = Vec::new();
        for dof in 0..3*self.reference.nodes_m.len() {
            check(cx)?;
            if self.boundary.map.index[dof].is_none() {
                supports.push((dof, self.state.displacement_m()[dof/3][dof%3]));
            }
        }
        let electrodes = self.boundary.electrode_nodes.iter()
            .zip(&self.state.electrode_conditions)
            .map(|(nodes, &condition)| PiezoelectricElectrode { nodes, condition })
            .collect::<Vec<_>>();
        let zero = vec![[0.0; 3]; self.reference.nodes_m.len()];
        let mut body = TetPiezoelectricTransient::new(model, &supports, &electrodes,
            self.state.displacement_m(), &zero, &self.state.nodal_free_charge_c, controls, cx,
        ).map_err(spatial)?;
        body.boundary = Arc::clone(&self.boundary);
        body.state.velocity_m_s.clone_from(&self.state.velocity_m_s);
        body.state.kinetic_energy_j = self.state.kinetic_energy_j;
        body.state.time_s = self.state.time_s;
        Ok(body)
    }

    /// Jointly advance source-dependent 3-D fields, motion, shunt heating and
    /// finite caloric/conductive state. Each Newton proposal starts from the
    /// complete accepted history. Supports may move through explicit endpoint
    /// velocities; omitted motion requires stationary supports as in the native
    /// API. All voltage drives must remain the fixed zero-potential reference.
    ///
    /// Tensors are reassembled under the original mesh/assembly budgets.
    /// Reversible field/caloric transfer is signed and kept separate from
    /// source-bound resistor heat. Refusal/cancellation publishes no history.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate<'thermal>(
        &self, thermal: &DepositionSimulation<'thermal>, target_s: f64,
        material: &PreparedMaterialNetwork<'_, '_>, bindings: &[ThermalElectrodeLoad],
        drives: &[ElectrodeDrive], midpoint_nodal_forces_n: &[[f64; 3]],
        endpoint_nodal_free_charge_c: &[f64],
        endpoint_support_velocities_m_s: Option<&[(usize, f64)]>,
        field_controls: TetPiezoelectricControls, thermal_controls: CircuitEnthalpySolve,
        guess_j: &[f64], cx: &Cx<'_>,
    ) -> Result<TensorThermalCandidate<'thermal, 'model, 'mesh>> {
        check(cx)?;
        crate::piezoelectric::spatial::validate_controls(field_controls).map_err(spatial)?;
        admit_driver_unknowns(thermal, thermal_controls, guess_j)?;
        let dt = target_s - thermal.time_s();
        if self.state.time_s() != thermal.time_s() || !dt.is_finite() || dt <= 0.0
            || drives.iter().any(|drive| matches!(drive, ElectrodeDrive::Voltage(v) if *v != 0.0))
        {
            return Err(invalid("shared positive clock and fixed zero-voltage references required"));
        }
        let cancelled = || cx.checkpoint().is_err();
        let original_network = material.resolve(thermal, cancelled).map_err(network)?;
        // Bind accepted resistors too: an endpoint query cannot repair a missing
        // original source. The exact driver connection is checked once here.
        bind_loads(original_network.resistors, bindings)?;
        let before = self.admit(thermal, material, field_controls, cx)?;
        let temperature_before = self.temperature(thermal, material)?;
        let mut routes = material.resistor_destinations().into_iter()
            .map(|route| (route.branch, route)).collect::<BTreeMap<_, _>>();
        if bindings.len() != routes.len() {
            return Err(invalid("every original resistor needs one electrode binding"));
        }
        let destinations = bindings.iter().map(|binding|
            routes.remove(&binding.resistor_branch)
                .ok_or_else(|| invalid("duplicate or foreign thermal electrode binding"))
        ).collect::<Result<Vec<_>>>()?;
        if !routes.is_empty() { return Err(invalid("unconnected material-network resistor")); }
        let ceramic_destinations = (0..self.reference.tetrahedra.len()).map(|piezoelectric_cell|
            ThermalCeramicCell { piezoelectric_cell, thermal_cell: self.caloric_cell }
        ).collect::<Vec<_>>();
        let (thermal_next, body, coupled) = cell_driver_candidate(
            thermal, target_s, thermal_controls, guess_j,
            |trial, _remaining| {
                check(cx)?;
                let mut coefficients = material.resolve(trial, cancelled).map_err(network)?;
                let temperature = self.temperature(trial, material)?;
                let samples = self.laws.iter().map(|law|
                    law.at_temperature(temperature).map_err(field)
                ).collect::<Result<Vec<_>>>()?;
                // Identical-temperature trials reuse the accepted native
                // operators and field history. Changed temperatures own a new
                // complete operator snapshot under the original budgets.
                let unchanged = temperature.to_bits() == temperature_before.to_bits();
                let accepted_model = if unchanged {
                    self.accepted_model.as_ref().map(Arc::clone)
                } else {
                    Some(Arc::new(TetPiezoelectricModel::assemble_owned_materials(
                        self.reference.nodes_m, self.reference.tetrahedra,
                        samples.iter().map(|sample| sample.material().clone()).collect(),
                        self.reference.mesh_frame, self.reference.budget, cx,
                    ).map_err(spatial)?))
                };
                let loads = bind_loads(coefficients.resistors, bindings)?;
                let (mut state, native, loads, midpoint_support_force_n, frozen_states) = {
                    let model = accepted_model.as_deref().unwrap_or(self.reference);
                    let rebound = if unchanged {
                        TetPiezoelectricTransient {
                            model, boundary: Arc::clone(&self.boundary),
                            state: self.state.clone(), intrinsic: None,
                        }
                    } else {
                        self.rebound(model, field_controls, cx)?
                    };
                    let frozen = self.field_states(rebound.state(), &samples, field_controls, cx)?;
                    let (mut candidate, loads) = rebound.advance_loaded_candidate(
                        dt, endpoint_support_velocities_m_s, drives, &loads,
                        midpoint_nodal_forces_n, endpoint_nodal_free_charge_c, field_controls, cx,
                    ).map_err(spatial)?;
                    let reactions = std::mem::take(&mut candidate.midpoint_support_force_n);
                    let (next, native) = candidate.into_parts();
                    (next.state, native, loads, reactions, frozen)
                };
                state.time_s = target_s;
                let after = self.field_states(&state, &samples, field_controls, cx)?;
                let mut parameter_work_j = 0.0;
                let mut delta_ts = 0.0;
                let mut internal = [0.0; 2];
                let mut cells = Vec::with_capacity(before.len());
                for ((old, frozen), new) in before.iter().zip(&frozen_states).zip(&after) {
                    check(cx)?;
                    parameter_work_j = finite(parameter_work_j + frozen.free_energy_j-old.free_energy_j)?;
                    delta_ts = finite(delta_ts + temperature*new.entropy_j_k
                        - temperature_before*old.entropy_j_k)?;
                    internal[0] = finite(internal[0] + old.internal_energy_j)?;
                    internal[1] = finite(internal[1] + new.internal_energy_j)?;
                    cells.push([*old, *new]);
                }
                let reversible_heat_j = finite(-parameter_work_j-delta_ts)?;
                let power = finite(reversible_heat_j/dt)?;
                if reversible_heat_j != 0.0 && power == 0.0 {
                    return Err(invalid("ceramic reversible transfer underflows"));
                }
                let total_source = coefficients.sources_w.entry(self.caloric_cell).or_insert(0.0);
                *total_source = finite(*total_source + power)?;
                let body = Self {
                    reference: self.reference, accepted_model,
                    boundary: Arc::clone(&self.boundary), state, laws: Arc::clone(&self.laws),
                    caloric_cell: self.caloric_cell,
                    maximum_displacement_norm_c_m2: self.maximum_displacement_norm_c_m2,
                };
                Ok(DrivenCoefficients {
                    circuit: TensorDriver { body, step: SpatialPiezoelectricThermalStep {
                        native, loads, midpoint_support_force_n,
                        destinations: destinations.clone(),
                        ceramic_destinations: ceramic_destinations.clone(),
                    }},
                    sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k,
                    source: TensorThermalSource {
                        network: coefficients.source, caloric_cell: self.caloric_cell,
                        temperature_k: [temperature_before, temperature], cells,
                        field_internal_energy_j: internal, parameter_work_j, reversible_heat_j,
                    },
                    circuit_trials: 1,
                })
            },
            |driver, dt, checkpoint| {
                let mut powers = BTreeMap::new();
                for (heat, route) in driver.step.loads.loads.iter().zip(&driver.step.destinations) {
                    checkpoint()?;
                    let power = heat.joule_heat_j/dt;
                    let total = powers.entry(route.cell).or_insert(0.0);
                    *total += power;
                    if !power.is_finite() || power < 0.0 || !total.is_finite()
                        || heat.joule_heat_j != 0.0 && power == 0.0
                    { return Err(DepositionError::Arithmetic.into()); }
                }
                Ok(powers)
            },
            cancelled,
        )?;
        let external_heat_j = finite(coupled.non_electrical_heat_j-coupled.source.reversible_heat_j)?;
        let total_energy_residual_j = finite(
            coupled.source.field_internal_energy_j[1]-coupled.source.field_internal_energy_j[0]
                + coupled.electrical.native.kinetic_energy_change_j
                + coupled.thermal.stored_energy_change_j
                - coupled.electrical.native.mechanical_work_j
                - coupled.electrical.loads.electrical_source_work_j
                - external_heat_j,
        )?;
        if total_energy_residual_j.abs() > thermal_controls.total_energy_tolerance_j {
            return Err(invalid("full tensor field U, kinetic energy, caloric U and external supply disagree"));
        }
        check(cx)?;
        Ok(TensorThermalCandidate {
            thermal: thermal_next, body,
            report: TensorThermalReport { coupled, external_heat_j, total_energy_residual_j },
        })
    }
}
