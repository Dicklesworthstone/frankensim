//! Finite heating, thermal bending and prestress on the original frame timeline.
//!
//! Every trial uses the native implicit caloric residual and native linear-
//! geometry force-based equilibrium or Newmark motion at held nodal loads.
//! The original free-sided laminate potential supplies each layer's endpoint-
//! average stress work and internal-energy change. Their difference enters its original
//! caloric cell once, with either sign. Independently integrated nodal work,
//! cell heat, thermal energy and combined energy must all pass admission.
//!
//! Both thermal endpoint schemes are generally first order in physical
//! time. Refine the thermal cells, member stations and time intervals separately.
//! Lobatto tributary volumes are the explicit thermal spatial approximation;
//! they do not resolve temperature inside a declared uniform layer. Reference
//! geometry, mass and transport paths stay fixed. The explicit inertial route
//! retains the original native mass and undamped average-acceleration Newmark
//! owner. Finite rotations, phase changes, material damping and deformation-
//! dependent transport require their own owners.

use crate::circuit_enthalpy::{
    CircuitEnthalpyError, CircuitEnthalpySolve, DrivenCoefficients, DriverAccount, ThermalDriver,
    cell_driver_candidate,
    material::{
        MaterialNetwork, MaterialNetworkError, MaterialNetworkSample, PreparedMaterialNetwork,
    },
};
use fs_exec::Cx;
use fs_material::laminate::thermal::beam::thermoelastic::{
    LaminateBeamThermoelasticError, LaminateBeamThermoelasticExchange,
};
use fs_material::phase::PhaseStateError;
use fs_math::dd::Dd;
use fs_solid::frame2d::{DynamicHistory, laminate_thermal::dynamics::FrameLaminateNewmarkControls};
use fs_solid::{SolidError, frame2d::Frame2d};
use fs_solver::NewtonReport;
use fs_thermal::enthalpy::{CellId, DepositedCell, DepositionSimulation, EnergyAudit, ThermalLaw};
use std::{collections::BTreeMap, sync::Arc};

/// Declared thermodynamic role of the original caloric chart.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LaminateCaloricReference {
    /// The chart describes this same unloaded solid material path. Its energy
    /// convention is compatible with adding the elastic potential's U=Psi+TS;
    /// maintained-pressure work is already consistently removed or negligible
    /// in the declared small-strain model. A common card alone does not prove
    /// these physical assumptions, and this declaration does not convert h to U.
    StressFreeSolid,
}

/// One original station/layer and its receiving finite caloric inventory.
/// Several points may share a cell only as an explicit uniform-T approximation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct LaminateFrameCell {
    pub element: usize,
    pub station: usize,
    pub layer: usize,
    pub cell: CellId,
}

#[derive(Debug)]
pub enum LaminateFrameSourceError {
    Mechanical(SolidError),
    Potential(LaminateBeamThermoelasticError),
    Caloric(PhaseStateError),
    Network(MaterialNetworkError),
}
impl core::fmt::Display for LaminateFrameSourceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "laminate-frame thermal source: {self:?}")
    }
}
impl std::error::Error for LaminateFrameSourceError {}
pub type LaminateFrameEnthalpyError = CircuitEnthalpyError<LaminateFrameSourceError>;
type Result<T> = core::result::Result<T, LaminateFrameEnthalpyError>;

fn invalid(why: &'static str) -> LaminateFrameEnthalpyError {
    CircuitEnthalpyError::Invalid(why)
}
fn check(cx: &Cx<'_>) -> Result<()> {
    cx.checkpoint().map_err(|_| CircuitEnthalpyError::Cancelled)
}
fn mechanical(error: SolidError) -> LaminateFrameEnthalpyError {
    if matches!(error, SolidError::Cancelled) {
        CircuitEnthalpyError::Cancelled
    } else {
        CircuitEnthalpyError::Source(LaminateFrameSourceError::Mechanical(error))
    }
}
fn network(error: MaterialNetworkError) -> LaminateFrameEnthalpyError {
    if matches!(
        error,
        MaterialNetworkError::Thermal(fs_thermal::enthalpy::DepositionError::Cancelled)
    ) {
        CircuitEnthalpyError::Cancelled
    } else {
        CircuitEnthalpyError::Source(LaminateFrameSourceError::Network(error))
    }
}
fn caloric(error: PhaseStateError) -> LaminateFrameEnthalpyError {
    CircuitEnthalpyError::Source(LaminateFrameSourceError::Caloric(error))
}
fn finite(x: f64) -> Result<f64> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(invalid("unrepresentable laminate-frame thermal account"))
    }
}
fn same_measure(a: f64, b: f64) -> bool {
    a.is_finite()
        && b.is_finite()
        && a > 0.0
        && b > 0.0
        && (a - b).abs() <= 512.0 * f64::EPSILON * a.abs().max(b.abs())
}

#[derive(Clone)]
struct Cell {
    original: DepositedCell,
    source_index: usize,
    minimum_capacity_j_k: Option<f64>,
}
#[derive(Clone)]
struct Station {
    element: usize,
    station: usize,
    weight_m: f64,
    cells: Vec<usize>,
}

/// Native accepted mechanics and immutable original thermal correspondence.
/// Quasistatic construction explicitly equilibrates the frame at its current
/// temperatures and held load. That equilibrium defines the initial state;
/// initialization is not a simulated heating interval. The explicit dynamic
/// constructor instead retains the actual displacement and velocity, including
/// unequilibrated launch states and impulses. These evolution choices remain
/// fixed across accepted steps and retries.
#[derive(Clone)]
pub struct ReciprocalLaminateFrame<'n, 'a> {
    frame: Frame2d,
    source: &'n MaterialNetwork<'a>,
    transport: Arc<PreparedMaterialNetwork<'n, 'a>>,
    cells: Vec<Cell>,
    stations: Vec<Station>,
    held_load: Vec<f64>,
    time_s: f64,
    dynamics: Option<FrameLaminateNewmarkControls>,
}

/// One source-owned exchange per reference length; multiply by `weight_m` once.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LaminateFrameLayerExchange {
    pub location: LaminateFrameCell,
    pub weight_m: f64,
    pub source: LaminateBeamThermoelasticExchange,
}

#[derive(Clone, Debug, PartialEq)]
pub struct LaminateFrameMechanicalExchange {
    pub layers: Vec<LaminateFrameLayerExchange>,
    pub stress_work_j: f64,
    /// Held original nodal loads dotted with the original native displacement
    /// increment (the represented increment for quasistatic equilibrium).
    pub external_work_j: f64,
    /// Layer stress work plus native kinetic change minus nodal-load work.
    pub work_residual_j: f64,
    pub internal_energy_change_j: f64,
    /// Signed exchange into the caloric owner, never numerical dissipation.
    pub reversible_heat_j: f64,
    /// Original native inertial interval, absent for quasistatic equilibrium.
    pub dynamics: Option<LaminateFrameDynamicExchange>,
}

/// Actual original Newmark trajectory and its independently audited work.
#[derive(Clone, Debug, PartialEq)]
pub struct LaminateFrameDynamicExchange {
    /// Times are relative to the beginning of this coupled interval.
    pub history: DynamicHistory,
    /// Stable average(v) dot M delta(v), not absolute kinetic subtraction.
    pub kinetic_energy_change_j: f64,
    pub native_internal_work_j: f64,
    pub native_external_work_j: f64,
    pub native_energy_residual_j: f64,
    /// Zero for this undamped average-acceleration method. A residual is never
    /// relabeled as numerical loss, and numerical loss is never caloric heat.
    pub numerical_dissipation_j: f64,
}

pub struct LaminateFrameEnthalpyReport {
    pub mechanical: LaminateFrameMechanicalExchange,
    pub thermal: EnergyAudit,
    pub source: MaterialNetworkSample,
    pub nonlinear: NewtonReport,
    pub material_trials: usize,
    pub maximum_cell_residual_j: f64,
    /// Authored boundary/source heat excluding reciprocal material exchange.
    pub external_thermal_heat_j: f64,
    pub total_energy_residual_j: f64,
}

/// Both original histories advance together only after independent acceptance.
pub struct LaminateFrameEnthalpyCandidate<'m, 'n, 'a> {
    thermal: DepositionSimulation<'m>,
    frame: ReciprocalLaminateFrame<'n, 'a>,
    report: LaminateFrameEnthalpyReport,
}
impl<'m, 'n, 'a> LaminateFrameEnthalpyCandidate<'m, 'n, 'a> {
    pub fn thermal(&self) -> &DepositionSimulation<'m> {
        &self.thermal
    }
    pub fn frame(&self) -> &ReciprocalLaminateFrame<'n, 'a> {
        &self.frame
    }
    pub fn report(&self) -> &LaminateFrameEnthalpyReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'m>,
        ReciprocalLaminateFrame<'n, 'a>,
        LaminateFrameEnthalpyReport,
    ) {
        (self.thermal, self.frame, self.report)
    }
}

impl<'n, 'a> ReciprocalLaminateFrame<'n, 'a> {
    /// Bind every original frame station/layer exactly once. The network may
    /// include passive thermal cells, but cannot contain unowned resistors.
    /// Mapped cells must equal the complete sum of their original layer volumes
    /// and masses. Source-qualified solid calorics, density claims and fixed
    /// material contexts must agree, including their typed temperature axis.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        frame: Frame2d,
        thermal: &DepositionSimulation<'_>,
        source: &'n MaterialNetwork<'a>,
        mapping: &[LaminateFrameCell],
        held_load: &[f64],
        caloric_reference: LaminateCaloricReference,
        maximum_layer_points: usize,
        cx: &Cx<'_>,
    ) -> Result<Self> {
        Self::bind(
            frame,
            thermal,
            source,
            mapping,
            held_load,
            caloric_reference,
            maximum_layer_points,
            None,
            cx,
        )
    }

    /// Bind actual native motion to finite caloric storage without static
    /// reinitialization or a velocity reset. Every thermal trial starts from
    /// this SAME displacement/velocity and original-temperature force. Native
    /// Newmark computes consistent initial acceleration from the held load;
    /// its endpoint force uses the newly solved layer temperatures. Initial
    /// kinetic energy, original mass and attached source/strength selections
    /// remain physical state. Only the original undamped linear-geometry
    /// laminate Newmark domain is admitted, using its physical work controls.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new_dynamic(
        frame: Frame2d,
        thermal: &DepositionSimulation<'_>,
        source: &'n MaterialNetwork<'a>,
        mapping: &[LaminateFrameCell],
        held_load: &[f64],
        caloric_reference: LaminateCaloricReference,
        maximum_layer_points: usize,
        dynamics: FrameLaminateNewmarkControls,
        cx: &Cx<'_>,
    ) -> Result<Self> {
        Self::bind(
            frame,
            thermal,
            source,
            mapping,
            held_load,
            caloric_reference,
            maximum_layer_points,
            Some(dynamics),
            cx,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn bind(
        mut frame: Frame2d,
        thermal: &DepositionSimulation<'_>,
        source: &'n MaterialNetwork<'a>,
        mapping: &[LaminateFrameCell],
        held_load: &[f64],
        _caloric_reference: LaminateCaloricReference,
        maximum_layer_points: usize,
        dynamics: Option<FrameLaminateNewmarkControls>,
        cx: &Cx<'_>,
    ) -> Result<Self> {
        check(cx)?;
        if maximum_layer_points == 0
            || mapping.is_empty()
            || mapping.len() > maximum_layer_points
            || frame.element_count() == 0
            || !source.resistors.is_empty()
            || held_load.len() != frame.displacements().len()
            || held_load.iter().any(|v| !v.is_finite())
            || frame
                .velocities()
                .iter()
                .chain(frame.accelerations())
                .any(|v| !v.is_finite())
            || dynamics.is_none()
                && frame
                    .velocities()
                    .iter()
                    .chain(frame.accelerations())
                    .any(|&v| v != 0.0)
        {
            return Err(invalid(
                "complete bounded laminate frame, finite held loads/motion and zero quasistatic inertial history required",
            ));
        }
        let transport = Arc::new(
            source
                .prepare(thermal, || cx.checkpoint().is_err())
                .map_err(network)?,
        );
        let mut positions = BTreeMap::new();
        for item in mapping {
            check(cx)?;
            if positions
                .insert((item.element, item.station, item.layer), item.cell)
                .is_some()
            {
                return Err(invalid(
                    "duplicate original laminate station/layer correspondence",
                ));
            }
        }
        let source_indices: BTreeMap<_, _> = source
            .cells
            .iter()
            .enumerate()
            .map(|(i, s)| (s.cell, i))
            .collect();
        let mut cells: Vec<_> = thermal
            .cells()
            .values()
            .map(|cell| Cell {
                original: cell.clone(),
                source_index: source_indices[&cell.id],
                minimum_capacity_j_k: None,
            })
            .collect();
        let indices: BTreeMap<_, _> = cells
            .iter()
            .enumerate()
            .map(|(i, c)| (c.original.id, i))
            .collect();
        let mut volumes = vec![Dd::ZERO; cells.len()];
        let mut masses = vec![Dd::ZERO; cells.len()];
        let mut elastic_floors = vec![Dd::ZERO; cells.len()];
        let mut stations = Vec::new();
        for element in 0..frame.element_count() {
            check(cx)?;
            frame.reference_strain_traces(element).map_err(mechanical)?;
            for station in 0..frame.section_deformations(element).len() {
                let actual = frame
                    .laminate_thermoelastic_response(element, station, cx)
                    .map_err(mechanical)?;
                let section = frame
                    .laminate_thermal_section(element, station)
                    .ok_or_else(|| {
                        invalid("all native members must retain original laminate sources")
                    })?;
                let mut receivers = Vec::with_capacity(actual.layers.len());
                for layer in &actual.layers {
                    check(cx)?;
                    let id = positions
                        .remove(&(element, station, layer.layer))
                        .ok_or_else(|| {
                            invalid("every original laminate station/layer requires one cell")
                        })?;
                    let index = *indices
                        .get(&id)
                        .ok_or_else(|| invalid("laminate caloric cell is inactive"))?;
                    let cell = &cells[index];
                    let material = &source.cells[cell.source_index];
                    let original = &section.source().source().layers()[layer.layer];
                    let input = original.input();
                    let reference = &input.reference;
                    let curve = material.curve;
                    let range = section
                        .thermoelastic_temperature_range_k(layer.layer)
                        .ok_or_else(|| {
                            invalid("original complete thermoelastic interval is absent")
                        })?;
                    let knots = curve.knots();
                    if material.point.material.content_hash() != input.card.content_hash()
                        || material.point.temperature_axis != "T"
                        || material.point.liquid_fraction_axis.is_some()
                        || material
                            .point
                            .context
                            .clone()
                            .with("T", reference.axes()["T"])
                            .ok()
                            .as_ref()
                            != Some(reference)
                        || curve.source_receipts().is_empty()
                        || knots.iter().any(|k| k.liquid_mass_fraction != 0.0)
                        || knots[0].temperature_k > range[0]
                        || knots[knots.len() - 1].temperature_k < range[1]
                        || thermal.temperature_k(id)? != layer.temperature_k
                    {
                        return Err(invalid(
                            "laminate and unloaded caloric source, reference context, interval or temperature disagree",
                        ));
                    }
                    for receipt in curve.source_receipts() {
                        check(cx)?;
                        if receipt.axis_quantities != *reference.axis_quantities()
                            || receipt.query_point.len() != reference.axes().len()
                            || receipt.query_point.iter().any(|(axis, value)| {
                                axis != "T" && reference.axes().get(axis) != Some(value)
                            })
                            || receipt.property == "density"
                                && receipt.selected
                                    != original.reference_density().answer().receipt.selected
                        {
                            return Err(invalid(
                                "original caloric conditions or reference-density claim differ from the laminate",
                            ));
                        }
                    }
                    let volume = finite(
                        actual.reference_length_weight_m
                            * section.prepared().width_m()
                            * input.thickness_m,
                    )?;
                    let mass = finite(volume * original.reference_density_kg_m3())?;
                    if volume <= 0.0 || mass <= 0.0 {
                        return Err(invalid(
                            "positive original laminate volume and mass required",
                        ));
                    }
                    volumes[index] = volumes[index] + Dd::from_f64(volume);
                    masses[index] = masses[index] + Dd::from_f64(mass);
                    elastic_floors[index] = elastic_floors[index]
                        + Dd::from_f64(actual.reference_length_weight_m)
                            * Dd::from_f64(layer.minimum_elastic_heat_capacity_j_m_k);
                    receivers.push(index);
                }
                stations.push(Station {
                    element,
                    station,
                    weight_m: actual.reference_length_weight_m,
                    cells: receivers,
                });
            }
        }
        if !positions.is_empty() {
            return Err(invalid(
                "correspondence names a foreign native station/layer",
            ));
        }
        for (i, cell) in cells.iter_mut().enumerate() {
            check(cx)?;
            if volumes[i] == Dd::ZERO {
                continue;
            }
            if !same_measure(cell.original.geometry.volume_m3, volumes[i].to_f64())
                || !same_measure(cell.original.mass_kg, masses[i].to_f64())
            {
                return Err(invalid(
                    "whole caloric cell must equal its original reference layer inventory",
                ));
            }
            let curve = source.cells[cell.source_index].curve;
            let mut cp_min = f64::INFINITY;
            for segment in 0..curve.knots().len() - 1 {
                check(cx)?;
                cp_min = cp_min.min(
                    curve
                        .sensible_heat_capacity_bounds_j_kg_k(segment)
                        .map_err(caloric)?[0],
                );
            }
            let floor = finite(
                (Dd::from_f64(cell.original.mass_kg) * Dd::from_f64(cp_min) + elastic_floors[i])
                    .to_f64()
                    .next_down(),
            )?;
            if floor <= 0.0 {
                return Err(invalid(
                    "caloric plus elastic capacity must stay positive over the complete source and strain envelope",
                ));
            }
            cell.minimum_capacity_j_k = Some(floor);
        }
        // Dynamic initialization must retain actual motion and prestress.
        if dynamics.is_none() {
            frame
                .static_load_with_cx(held_load, 1, cx)
                .map_err(mechanical)?;
        }
        for element in 0..frame.element_count() {
            frame
                .laminate_thermoelastic_energy(element, cx)
                .map_err(mechanical)?;
        }
        check(cx)?;
        Ok(Self {
            frame,
            source,
            transport,
            cells,
            stations,
            held_load: held_load.to_vec(),
            time_s: thermal.time_s(),
            dynamics,
        })
    }

    pub fn frame(&self) -> &Frame2d {
        &self.frame
    }
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    /// Immutable native evolution policy; None denotes quasistatic equilibrium.
    pub const fn dynamics(&self) -> Option<FrameLaminateNewmarkControls> {
        self.dynamics
    }
    /// Conservative total capacity of each mapped cell at fixed station strains.
    /// Native stable quasistatic mechanical relaxation adds a nonnegative term.
    pub fn minimum_cell_heat_capacities_j_k(&self) -> impl Iterator<Item = (CellId, f64)> + '_ {
        self.cells
            .iter()
            .filter_map(|c| c.minimum_capacity_j_k.map(|floor| (c.original.id, floor)))
    }

    fn check_thermal(
        &self,
        thermal: &DepositionSimulation<'_>,
        accepted: bool,
        cx: &Cx<'_>,
    ) -> Result<()> {
        if thermal.cells().len() != self.cells.len() {
            return Err(invalid("original thermal topology changed"));
        }
        for bound in &self.cells {
            check(cx)?;
            let old = &bound.original;
            let cell = thermal
                .cells()
                .get(&old.id)
                .ok_or_else(|| invalid("original thermal cell disappeared"))?;
            if cell.id != old.id
                || cell.segment_id != old.segment_id
                || cell.body_id != old.body_id
                || cell.material_id != old.material_id
                || cell.geometry != old.geometry
                || cell.mass_kg.to_bits() != old.mass_kg.to_bits()
                || accepted && cell.specific_h_j_kg.to_bits() != old.specific_h_j_kg.to_bits()
                || !thermal.uses_material_law(
                    cell.material_id,
                    self.source.cells[bound.source_index].curve as &dyn ThermalLaw,
                )
            {
                return Err(invalid(
                    "original caloric identity, reference geometry, mass or accepted enthalpy changed",
                ));
            }
        }
        Ok(())
    }

    fn proposal(
        &self,
        thermal: &DepositionSimulation<'_>,
        duration_s: f64,
        tolerance: f64,
        cx: &Cx<'_>,
    ) -> Result<Proposal> {
        self.check_thermal(thermal, false, cx)?;
        let mut field = Vec::with_capacity(self.frame.element_count());
        let mut start = 0;
        while start < self.stations.len() {
            check(cx)?;
            let element = self.stations[start].element;
            let end = start
                + self.stations[start..]
                    .iter()
                    .take_while(|s| s.element == element)
                    .count();
            let temperatures: Vec<Vec<f64>> = self.stations[start..end]
                .iter()
                .map(|station| {
                    station
                        .cells
                        .iter()
                        .map(|&i| {
                            thermal
                                .temperature_k(self.cells[i].original.id)
                                .map_err(Into::into)
                        })
                        .collect()
                })
                .collect::<Result<_>>()?;
            field.push(temperatures);
            start = end;
        }
        let (frame, kinetic, displacement_increment, dynamics) =
            if let Some(controls) = self.dynamics {
                let native = self
                    .frame
                    .laminate_newmark_candidate(&field, &self.held_load, duration_s, controls, cx)
                    .map_err(mechanical)?;
                let kinetic = native.retained_kinetic_energy_change_j();
                let report = LaminateFrameDynamicExchange {
                    history: native.history,
                    kinetic_energy_change_j: native.kinetic_energy_change_j,
                    native_internal_work_j: native.internal_work_j,
                    native_external_work_j: native.external_work_j,
                    native_energy_residual_j: native.energy_balance_residual_j,
                    numerical_dissipation_j: 0.0,
                };
                let increment: Vec<_> = native
                    .displacement_increment
                    .into_iter()
                    .map(Dd::from_f64)
                    .collect();
                (native.frame, kinetic, increment, Some(report))
            } else {
                let mut frame = self
                    .frame
                    .with_laminate_temperature_field(&field, cx)
                    .map_err(mechanical)?;
                frame
                    .static_load_with_cx(&self.held_load, 1, cx)
                    .map_err(mechanical)?;
                let increment = frame
                    .displacements()
                    .iter()
                    .zip(self.frame.displacements())
                    .map(|(&after, &before)| Dd::from_f64(after) - Dd::from_f64(before))
                    .collect();
                (frame, Dd::ZERO, increment, None)
            };
        let mut layers = Vec::new();
        let mut cell_heat = vec![Dd::ZERO; self.cells.len()];
        let (mut work, mut energy, mut heat) = (Dd::ZERO, Dd::ZERO, Dd::ZERO);
        for station in &self.stations {
            check(cx)?;
            let old = self
                .frame
                .laminate_thermal_section(station.element, station.station)
                .and_then(|s| s.thermoelastic())
                .ok_or_else(|| invalid("original laminate potential disappeared"))?;
            let new = frame
                .laminate_thermal_section(station.element, station.station)
                .and_then(|s| s.thermoelastic())
                .ok_or_else(|| invalid("trial laminate potential disappeared"))?;
            let (e0, k0) = self.frame.section_deformations(station.element)[station.station];
            let (e1, k1) = frame.section_deformations(station.element)[station.station];
            let weight = Dd::from_f64(station.weight_m);
            for (layer, &cell) in station.cells.iter().enumerate() {
                check(cx)?;
                let source = old
                    .layer_exchange(new, layer, [e0, k0], [e1, k1])
                    .map_err(|e| {
                        CircuitEnthalpyError::Source(LaminateFrameSourceError::Potential(e))
                    })?;
                work = work + weight * Dd::from_f64(source.mechanical_work_j_m);
                energy = energy + weight * Dd::from_f64(source.internal_energy_change_j_m);
                let q = weight * Dd::from_f64(source.reversible_heat_j_m);
                heat = heat + q;
                cell_heat[cell] = cell_heat[cell] + q;
                layers.push(LaminateFrameLayerExchange {
                    location: LaminateFrameCell {
                        element: station.element,
                        station: station.station,
                        layer,
                        cell: self.cells[cell].original.id,
                    },
                    weight_m: station.weight_m,
                    source,
                });
            }
        }
        let mut external = Dd::ZERO;
        for (&p, &increment) in self.held_load.iter().zip(&displacement_increment) {
            check(cx)?;
            external = external + Dd::from_f64(p) * increment;
        }
        if let Some(native) = &dynamics {
            for difference in [
                work - Dd::from_f64(native.native_internal_work_j),
                external - Dd::from_f64(native.native_external_work_j),
            ] {
                if finite(difference.to_f64())?.abs() > tolerance {
                    return Err(invalid(
                        "layer stress or nodal-load work differs from the original Newmark account",
                    ));
                }
            }
        }
        let residual = finite((work + kinetic - external).to_f64())?;
        if residual.abs() > tolerance {
            return Err(invalid(
                "original layer stress work and kinetic change do not match held native nodal-load work",
            ));
        }
        let account = LaminateFrameMechanicalExchange {
            layers,
            stress_work_j: finite(work.to_f64())?,
            external_work_j: finite(external.to_f64())?,
            work_residual_j: residual,
            internal_energy_change_j: finite(energy.to_f64())?,
            reversible_heat_j: finite(heat.to_f64())?,
            dynamics,
        };
        let cell_heat = self
            .cells
            .iter()
            .zip(cell_heat)
            .filter(|(c, _)| c.minimum_capacity_j_k.is_some())
            .map(|(c, q)| Ok((c.original.id, finite(q.to_f64())?)))
            .collect::<Result<_>>()?;
        Ok(Proposal {
            frame,
            account,
            cell_heat,
            stored_change_j: finite((energy + kinetic).to_f64())?,
            retained: energy + kinetic - external,
        })
    }

    /// Advance finite heat and compatible shape, prestress or native motion on
    /// one physical clock, using the immutable construction-time evolution.
    /// The existing implicit energy-coordinate Newton solve owns all stencils,
    /// work limits and acceptance. `maximum_circuit_trials` bounds mechanical/
    /// material proposals here; there is no electrical state. Original native
    /// transport/event restrictions apply. Each trial restarts both accepted
    /// states, including source-pinned strength, and failed work publishes none.
    pub fn step_candidate<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        controls: CircuitEnthalpySolve,
        guess_j: &[f64],
        cx: &Cx<'_>,
    ) -> Result<LaminateFrameEnthalpyCandidate<'m, 'n, 'a>> {
        check(cx)?;
        if thermal.time_s() != self.time_s {
            return Err(invalid("native mechanical and caloric clocks disagree"));
        }
        if self.dynamics.is_some_and(|d| {
            !d.maximum_energy_defect_j.is_finite()
                || d.maximum_energy_defect_j < 0.0
                || d.maximum_energy_defect_j > controls.total_energy_tolerance_j
        }) {
            return Err(invalid(
                "native mechanical energy allowance must fit the combined thermal allowance",
            ));
        }
        self.check_thermal(thermal, true, cx)?;
        let (thermal, frame, report) = cell_driver_candidate(
            thermal,
            target_s,
            controls,
            guess_j,
            |state, _remaining| {
                let transport = self
                    .transport
                    .resolve(state, || cx.checkpoint().is_err())
                    .map_err(network)?;
                let circuit = self.proposal(
                    state,
                    target_s - self.time_s,
                    controls.total_energy_tolerance_j,
                    cx,
                )?;
                Ok(DrivenCoefficients {
                    circuit,
                    sources_w: transport.sources_w,
                    conductances_w_k: transport.conductances_w_k,
                    source: transport.source,
                    circuit_trials: 1,
                })
            },
            |proposal: &Proposal, dt, poll| {
                let mut powers = BTreeMap::new();
                for (&id, &q) in &proposal.cell_heat {
                    poll()?;
                    let power = finite(q / dt)?;
                    if q != 0.0 && power == 0.0 {
                        return Err(invalid(
                            "reversible material heat disappeared during power conversion",
                        ));
                    }
                    powers.insert(id, power);
                }
                Ok(powers)
            },
            || cx.checkpoint().is_err(),
        )?;
        let mut next = Self {
            frame,
            source: self.source,
            transport: Arc::clone(&self.transport),
            cells: self.cells.clone(),
            stations: self.stations.clone(),
            held_load: self.held_load.clone(),
            time_s: thermal.time_s(),
            dynamics: self.dynamics,
        };
        for cell in &mut next.cells {
            cell.original = thermal.cells()[&cell.original.id].clone();
        }
        check(cx)?;
        Ok(LaminateFrameEnthalpyCandidate {
            thermal,
            frame: next,
            report: LaminateFrameEnthalpyReport {
                mechanical: report.electrical,
                thermal: report.thermal,
                source: report.source,
                nonlinear: report.nonlinear,
                material_trials: report.circuit_trials,
                maximum_cell_residual_j: report.maximum_cell_residual_j,
                external_thermal_heat_j: report.non_electrical_heat_j,
                total_energy_residual_j: report.total_energy_residual_j,
            },
        })
    }
}

struct Proposal {
    frame: Frame2d,
    account: LaminateFrameMechanicalExchange,
    cell_heat: BTreeMap<CellId, f64>,
    retained: Dd,
    stored_change_j: f64,
}
impl ThermalDriver for Proposal {
    type Model = Frame2d;
    type Step = LaminateFrameMechanicalExchange;
    fn account(&self) -> DriverAccount<'_> {
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: self.account.reversible_heat_j,
            numerical_dissipation_j: 0.0,
            stored_change_j: self.stored_change_j,
            supplied_j: self.account.external_work_j,
        }
    }
    fn retained_energy_balance_j(&self) -> Dd {
        self.retained
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        (self.frame, self.account)
    }
}
