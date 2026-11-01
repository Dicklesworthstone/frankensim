//! Pinned contact sources on the actual two native half-cell surfaces.

use std::sync::Arc;

use fs_exec::CancelGate;
use fs_matdb::{InterfaceSystemCard, QueryPoint};
use fs_material::state_point::{MaterialPropertySelection, ResolvedInterfaceStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, is_fixed_external_coordinate,
    resolve_interface_state_point};
use fs_material::thermal::interface::{AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as RD,
    AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as RP};
use fs_qty::{QuantitySpec, Temperature, semantic::{QuantityKind, SemanticType, ValueForm}};

use super::{ChtError, EnergySolution, FaceContactResistance, SolidMaterial, ThermalSetup,
    Voxel, VoxelDomain};
use super::super::conductance::{ConductivityCurve, interface_response};

/// One immutable interface source assigned to an ordered Cartesian cell face.
/// The caller assigns the card's A/B surface roles to `cells` in that order.
/// Native material indices do not infer a surface finish or source identity.
#[derive(Debug, Clone, PartialEq)]
pub struct SourceFaceContact {
    cells: [usize; 2],
    card: Arc<InterfaceSystemCard>,
    axes: [String; 2],
    reference: ResolvedInterfaceStatePoint,
}

impl SourceFaceContact {
    /// Select one positive static resistance, retaining its original context.
    /// Only the two distinct, explicit absolute-temperature axes evolve.
    /// Pressure, medium, finish and history remain at their supplied values;
    /// no phase, mechanical pressure or contact-area evolution is inferred.
    /// Geometry, occupancy and duplicate owners are checked by the native solve.
    pub fn new(
        cells: [usize; 2], card: Arc<InterfaceSystemCard>, point: &QueryPoint,
        temperature_axes: [&str; 2], selection: MaterialPropertySelection,
    ) -> Result<Self, ChtError> {
        if cells[0] == cells[1] || temperature_axes[0] == temperature_axes[1]
            || temperature_axes.iter().any(|axis| axis.trim().is_empty())
        { return Err(invalid("source contact needs distinct cells and temperature axes")); }
        let absolute = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature, ValueForm::Static));
        for axis in temperature_axes {
            if !point.axes().get(axis).is_some_and(|t| t.is_finite() && *t > 0.0)
                || !point.axis_quantities().get(axis).is_some_and(|q|
                    *q == absolute || *q == QuantitySpec::dimensional(Temperature::DIMS))
            { return Err(invalid("source contact requires typed positive absolute temperatures")); }
        }
        let requirement = ScalarPropertyRequirement::try_new(RP, RD, ScalarAdmissibility::StrictlyPositive)
            .map_err(|e| invalid(&format!("contact property requirement refused: {e}")))?;
        let reference = resolve_interface_state_point(&card, point, &[requirement], selection)
            .map_err(|e| invalid(&format!("contact source refused: {e}")))?;
        let property = reference.property(RP).ok_or_else(|| invalid("contact source omitted its resistance"))?;
        for (axis, dims) in property.source_axes() {
            let thermal = temperature_axes.contains(&axis.as_str()) && *dims == Temperature::DIMS;
            let quantity = point.axis_quantities().get(axis).copied();
            if !thermal && !(quantity.is_some_and(|q| q.dims() == *dims)
                && is_fixed_external_coordinate(axis, quantity, point.axes().get(axis).copied()))
            { return Err(invalid("contact response needs owned thermal or typed fixed external coordinates")); }
        }
        if point.axes().iter().any(|(axis, value)| {
            (axis == "omega" || point.axis_quantities().get(axis).and_then(|q| q.semantic_type())
                .is_some_and(|s| matches!(s.kind(), QuantityKind::Frequency(_)))) && *value != 0.0
        }) || property.angular_frequency_band_rad_s()
            .map_err(|e| invalid(&format!("contact static applicability refused: {e}")))?
            .is_some_and(|(low, _)| low > 0.0)
        { return Err(invalid("native thermal contact requires static source applicability")); }
        Ok(Self { cells, card, axes: temperature_axes.map(str::to_owned), reference })
    }

    /// Ordered cell assignment, matching the original card's A/B roles.
    #[must_use]
    pub const fn cells(&self) -> [usize; 2] { self.cells }
    /// Original source selection, quantities and full reference conditions.
    #[must_use]
    pub const fn reference(&self) -> &ResolvedInterfaceStatePoint { &self.reference }

    fn at(&self, surfaces: [f64; 2]) -> Result<ResolvedInterfaceStatePoint, ChtError> {
        let mut point = self.reference.query();
        for (axis, temperature) in self.axes.iter().zip(surfaces) {
            if !temperature.is_finite() || temperature <= 0.0 {
                return Err(invalid("contact surface temperature must be finite and positive"));
            }
            point = point.with(axis, temperature)
                .map_err(|e| invalid(&format!("contact surface coordinate refused: {e}")))?;
        }
        self.reference.at_query_point(&self.card, &point)
            .map_err(|e| invalid(&format!("contact source refused at its physical surfaces: {e}")))
    }
}

/// Numerical controls for the shared nonlinear source-contact endpoint.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SourceContactControls {
    /// Maximum complete fixed-resistance thermal solves per endpoint.
    pub max_iterations: usize,
    /// Symmetric relative difference between applied and freshly queried R.
    pub resistance_relative_tolerance: f64,
    /// Absolute contact-jump residual allowance, K, before roundoff allowance.
    pub temperature_tolerance_k: f64,
    /// Convex resistance update fraction in (0,1].
    pub relaxation: f64,
}
impl Default for SourceContactControls {
    fn default() -> Self { Self { max_iterations: 64, resistance_relative_tolerance: 1e-10,
        temperature_tolerance_k: 1e-9, relaxation: 1.0 } }
}

/// Physical source query at the returned cell field and its integrated surfaces.
/// Applied/source differences are retained and bounded by the requested controls.
#[derive(Debug, Clone, PartialEq)]
pub struct SourceContactState {
    /// Original A/B cell order; it is never replaced by a sorted source order.
    pub cells: [usize; 2],
    /// Solved half-cell surface temperatures in A/B order, K.
    pub surface_temperature_k: [f64; 2],
    /// Recomputed common heat transfer from first to second cell, W.
    pub heat_rate_first_to_second_w: f64,
    /// Fresh source resistance at the reported surface temperatures, m² K/W.
    pub resistance_m2_k_w: f64,
    /// Resistance used by the returned fixed-resistance thermal solve, m² K/W.
    pub applied_resistance_m2_k_w: f64,
    /// Symmetric relative applied/source resistance discrepancy.
    pub resistance_relative_residual: f64,
    /// Signed (T_surface_A - T_surface_B) - R_source*q, K.
    pub contact_temperature_residual_k: f64,
    /// Exact original source, selected claim and physical surface query receipt.
    pub source: ResolvedInterfaceStatePoint,
}

fn resistance(source: &ResolvedInterfaceStatePoint) -> Result<f64, ChtError> {
    source.property(RP).map(|p| p.value_si())
        .filter(|r| r.is_finite() && *r > 0.0)
        .ok_or_else(|| invalid("source contact resistance must remain finite and positive"))
}

/// One setup clone and one existing thermal solve per coupling pass. Every
/// caller's physical old state stays captured by `fixed_solve`; constitutive
/// iterations never become time steps. A refusal publishes no partial result.
pub(super) fn solve(
    domain: &VoxelDomain, solids: &[SolidMaterial], setup: &ThermalSetup,
    native_tables: bool, gate: &CancelGate,
    mut fixed_solve: impl FnMut(&ThermalSetup) -> Result<EnergySolution, ChtError>,
) -> Result<EnergySolution, ChtError> {
    super::super::poll(gate)?;
    let controls = setup.source_contact_controls;
    if controls.max_iterations == 0
        || !controls.resistance_relative_tolerance.is_finite() || controls.resistance_relative_tolerance <= 0.0
        || !controls.temperature_tolerance_k.is_finite() || controls.temperature_tolerance_k <= 0.0
        || !controls.relaxation.is_finite() || controls.relaxation <= 0.0 || controls.relaxation > 1.0
        || setup.source_face_contacts.len() > domain.cell_count().saturating_mul(3)
    { return Err(invalid("source-contact controls or face count are inadmissible")); }
    let mut fixed = setup.clone();
    fixed.source_face_contacts.clear();
    let start = fixed.face_contacts.len();
    for (index, law) in setup.source_face_contacts.iter().enumerate() {
        if index.is_multiple_of(256) { super::super::poll(gate)?; }
        fixed.face_contacts.push(FaceContactResistance {
            cells: law.cells, resistance_m2_k_w: resistance(&law.reference)?,
        });
    }
    let mut maximum_residual = f64::INFINITY;
    for iteration in 1..=controls.max_iterations {
        super::super::poll(gate)?;
        // Native face admission owns adjacency, duplicate contacts, compact
        // exclusion and zero through-flow, including all newly appended faces.
        let mut next = fixed_solve(&fixed)?;
        let mut states = Vec::with_capacity(setup.source_face_contacts.len());
        maximum_residual = 0.0;
        for (index, law) in setup.source_face_contacts.iter().enumerate() {
            if index.is_multiple_of(256) { super::super::poll(gate)?; }
            let [a, b] = law.cells;
            let axis = (0..6).find(|&face| domain.neighbor(a, face) == Some(b))
                .ok_or_else(|| invalid("source-contact face lost its admitted adjacency"))? / 2;
            let curve = |cell: usize| {
                if native_tables && let Voxel::Solid(m) = domain.voxel_at(cell) {
                    let material = &solids[usize::from(m)];
                    if !material.conductivity_table.is_empty() {
                        return ConductivityCurve::tabulated(material.conductivity_w_m_k,
                            material.axis_conductivity()[axis], &material.conductivity_table);
                    }
                }
                ConductivityCurve::constant(next.conductivity[cell][axis])
            };
            let applied = fixed.face_contacts[start + index].resistance_m2_k_w;
            let response = interface_response(curve(a), next.temperature[a], 0.5 * domain.widths(a)[axis],
                curve(b), next.temperature[b], 0.5 * domain.widths(b)[axis], applied)
                .ok_or_else(|| invalid("source-contact half-cell response is unrepresentable"))?;
            let surfaces = response.surface_temperature_k;
            let source = law.at(surfaces)?;
            let actual = resistance(&source)?;
            let flux = response.conductance_w_m2_k * (next.temperature[a] - next.temperature[b]);
            let heat = domain.face_area(a, axis) * flux;
            let residual = (surfaces[0] - surfaces[1]) - actual * flux;
            let relative = (actual - applied).abs() / actual.max(applied);
            let roundoff = 64.0 * f64::EPSILON * surfaces[0].abs()
                + 64.0 * f64::EPSILON * surfaces[1].abs()
                + 64.0 * f64::EPSILON * (actual * flux).abs();
            if [flux, heat, residual, relative, roundoff].iter().any(|v| !v.is_finite()) {
                return Err(invalid("source-contact heat or consistency residual is unrepresentable"));
            }
            maximum_residual = maximum_residual.max(relative / controls.resistance_relative_tolerance)
                .max((residual.abs() - roundoff).max(0.0) / controls.temperature_tolerance_k);
            states.push(SourceContactState { cells: law.cells, surface_temperature_k: surfaces,
                heat_rate_first_to_second_w: heat, resistance_m2_k_w: actual,
                applied_resistance_m2_k_w: applied, resistance_relative_residual: relative,
                contact_temperature_residual_k: residual, source });
        }
        if maximum_residual <= 1.0 {
            super::super::poll(gate)?;
            next.source_contact_states = states;
            next.source_contact_iterations = iteration;
            return Ok(next);
        }
        for (index, (applied, state)) in fixed.face_contacts[start..].iter_mut().zip(states).enumerate() {
            if index.is_multiple_of(256) { super::super::poll(gate)?; }
            let r = (1.0 - controls.relaxation) * applied.resistance_m2_k_w
                + controls.relaxation * state.resistance_m2_k_w;
            if !r.is_finite() || r <= 0.0 {
                return Err(invalid("relaxed contact resistance is unrepresentable"));
            }
            applied.resistance_m2_k_w = r;
        }
    }
    Err(ChtError::SolverNotConverged { system: "source thermal contacts",
        iterations: controls.max_iterations, relative_residual: maximum_residual, tolerance: 1.0 })
}

fn invalid(reason: &str) -> ChtError {
    ChtError::InvalidInput { field: "thermal.source_face_contacts", reason: reason.into() }
}
