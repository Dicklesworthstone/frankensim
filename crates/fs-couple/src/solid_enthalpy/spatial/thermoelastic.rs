//! Reciprocal heat on the original nodal enthalpy and tetrahedral equilibrium.
//!
//! The existing P1 heat owner retains every caloric chart, reference nodal mass,
//! conductivity, boundary and contact. Each native tetrahedron uses the existing
//! enthalpy-centroid material quadrature: T = curve^{-1}(sum(h_a)/4). Its original
//! source potential supplies endpoint stress work W and internal-energy change
//! dU. The signed exchange W-dU is assembled as one quarter at each original
//! vertex, in joules, and enters the original heat residual exactly once.
//!
//! This is a quasistatic, generally first-order endpoint coupling on stationary
//! reference geometry. Spatial and temporal refinement remain necessary. The
//! centroid material temperature is not a claim of resolved within-element
//! temperature or exact nodal entropy production. The physical material heat
//! is separate from every numerical residual. Native pointwise radiation stays
//! with its original boundary owner. The dynamic child retains the original
//! inertial trajectory; phase changes and moving thermal geometry remain outside
//! this small-strain, stationary-reference model.

/// Reciprocal nodal heat and the original full tetrahedral inertia.
pub mod dynamic;

use super::*;
use crate::solid_enthalpy::thermoelastic::{MAX_THERMOELASTIC_ITERATIONS, ThermoelasticSolve};
use fs_conduction::transient::enthalpy::nodal_heat::{NodalHeatEndpoint, NodalHeatStep};
use fs_conduction::interface::temperature::{
    TemperatureContacts, enthalpy::ContactEnthalpyNewtonConfig,
};
use fs_material::elastic::thermal::{
    AffineThermoelastic, ThermoelasticError, ThermoelasticState, ThermoelasticSymmetry,
    exchange::ThermoelasticExchange,
};
use fs_math::dd::Dd;
use std::sync::Arc;

/// Explicit thermodynamic role; a common source card alone cannot prove it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SpatialCaloricReference {
    /// The original specific-enthalpy chart is the same unloaded solid path.
    /// Its reference/work convention is compatible with adding Psi+T*S;
    /// maintained-pressure work is consistently removed or negligible in this
    /// small-strain model. This declaration does not convert arbitrary h to U.
    StressFreeSolid,
}

#[derive(Debug)]
pub enum SpatialThermoelasticError {
    Spatial(SpatialSolidEnthalpyError),
    Potential(ThermoelasticError),
    Invalid(&'static str),
    NotConverged {
        iterations: usize,
        heat_defect_j: f64,
        maximum_nodal_residual_j: f64,
        combined_residual_j: f64,
    },
    Cancelled,
}
impl core::fmt::Display for SpatialThermoelasticError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "reciprocal nodal thermoelasticity: {self:?}")
    }
}
impl std::error::Error for SpatialThermoelasticError {}
impl From<SpatialSolidEnthalpyError> for SpatialThermoelasticError {
    fn from(error: SpatialSolidEnthalpyError) -> Self {
        if matches!(error, SpatialSolidEnthalpyError::Cancelled) {
            Self::Cancelled
        } else {
            Self::Spatial(error)
        }
    }
}
impl From<EnthalpyError> for SpatialThermoelasticError {
    fn from(error: EnthalpyError) -> Self {
        SpatialSolidEnthalpyError::from(error).into()
    }
}
impl From<ThermoelasticError> for SpatialThermoelasticError {
    fn from(error: ThermoelasticError) -> Self { Self::Potential(error) }
}
type Result<T> = core::result::Result<T, SpatialThermoelasticError>;
fn check(cx: &Cx<'_>) -> Result<()> {
    cx.checkpoint().map_err(|_| SpatialThermoelasticError::Cancelled)
}
fn invalid(message: &'static str) -> SpatialThermoelasticError {
    SpatialThermoelasticError::Invalid(message)
}
fn finite(value: Dd) -> Result<f64> {
    let rounded = value.to_f64();
    if value.hi.is_finite() && value.lo.is_finite() && rounded.is_finite() {
        Ok(rounded)
    } else {
        Err(invalid("unrepresentable reciprocal nodal account"))
    }
}
fn dd(value: f64) -> Dd { Dd::from_f64(value) }
// Compare the retained physical account, not only its rounded report value.
fn within(value: Dd, tolerance: f64) -> Result<bool> {
    finite(value)?;
    Ok(!dd(tolerance).lt(value.abs()))
}
fn admit_controls(config: ThermoelasticSolve) -> Result<()> {
    if config.maximum_iterations == 0
        || config.maximum_iterations > MAX_THERMOELASTIC_ITERATIONS
        || !config.heat_relaxation.is_finite()
        || config.heat_relaxation <= 0.0 || config.heat_relaxation > 1.0
        || [config.heat_tolerance_j, config.energy_tolerance_j, config.stress_tolerance_pa]
            .iter().any(|v| !v.is_finite() || *v <= 0.0)
    {
        return Err(invalid("bounded iterations, positive tolerances and relaxation in (0,1] required"));
    }
    Ok(())
}
// One dispatch for both native mechanical trajectories. Live contact sets
// already include their constant interfaces, so a second interface set cannot
// silently occupy the same faces or apply the transfer twice.
#[derive(Clone, Copy)]
enum HeatBoundary<'b, 'mesh, 'card> {
    Fixed(Option<&'b ThermalInterfaces>, EnthalpyStepConfig),
    Temperature(&'b TemperatureContacts<'mesh, 'card>, ContactEnthalpyNewtonConfig),
}
impl HeatBoundary<'_, '_, '_> {
    fn step_config(self) -> EnthalpyStepConfig {
        match self {
            Self::Fixed(_, config) => config,
            Self::Temperature(_, config) => config.step,
        }
    }

    fn admits_contact_residual(self, residual_norm_j: f64) -> bool {
        match self {
            Self::Fixed(..) => true,
            Self::Temperature(_, config) => residual_norm_j <= config.residual_tolerance_j,
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn advance(
        self, source: &SpatialSolidEnthalpy<'_, '_, '_>, boundary: &ThermalBoundary,
        volumetric_source: &ScalarField, old_h: &[f64], dt_s: f64, heat_j: &[Dd],
        radiation: EnthalpyRadiation<'_>, cx: &Cx<'_>,
    ) -> Result<NodalHeatStep> {
        let problem = source.thermal_problem(boundary, volumetric_source)?;
        Ok(match self {
            Self::Fixed(interfaces, config) => source.thermal.advance_with_nodal_heat_and_radiation(
                cx, problem, interfaces, old_h, dt_s, config, heat_j, radiation,
            )?,
            Self::Temperature(contacts, config) => source.thermal.advance_with_nodal_heat_and_temperature_contacts(
                cx, problem, contacts, old_h, dt_s, config, heat_j, radiation,
            )?,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn evaluate(
        self, source: &SpatialSolidEnthalpy<'_, '_, '_>, boundary: &ThermalBoundary,
        volumetric_source: &ScalarField, old_h: &[f64], candidate_h: &[f64],
        dt_s: f64, heat_j: &[Dd], radiation: EnthalpyRadiation<'_>, cx: &Cx<'_>,
    ) -> Result<NodalHeatEndpoint> {
        let problem = source.thermal_problem(boundary, volumetric_source)?;
        Ok(match self {
            Self::Fixed(interfaces, config) => source.thermal.evaluate_endpoint_with_nodal_heat_and_radiation(
                cx, problem, interfaces, old_h, candidate_h, dt_s, config, heat_j, radiation,
            )?,
            Self::Temperature(contacts, config) => source.thermal.evaluate_endpoint_with_nodal_heat_and_temperature_contacts(
                cx, problem, contacts, old_h, candidate_h, dt_s, config, heat_j, radiation,
            )?,
        })
    }
}

fn stiffness(model: ThermalStiffnessModel) -> ThermoelasticSymmetry {
    match model {
        ThermalStiffnessModel::Isotropic => ThermoelasticSymmetry::Isotropic,
        ThermalStiffnessModel::IsotropicBulkShear => ThermoelasticSymmetry::IsotropicBulkShear,
        ThermalStiffnessModel::IsotropicYoungShear => ThermoelasticSymmetry::IsotropicYoungShear,
        ThermalStiffnessModel::OrientedOrthotropic { principal_to_mesh, .. } =>
            ThermoelasticSymmetry::OrientedOrthotropic { principal_to_target: principal_to_mesh },
    }
}
fn expansion(model: ThermalExpansionModel) -> ThermoelasticSymmetry {
    match model {
        ThermalExpansionModel::Isotropic => ThermoelasticSymmetry::Isotropic,
        ThermalExpansionModel::OrientedOrthotropic { principal_to_mesh, .. } =>
            ThermoelasticSymmetry::OrientedOrthotropic { principal_to_target: principal_to_mesh },
    }
}

/// Accepted original nodal heat, compatible equilibrium and immutable source laws.
/// Initialization explicitly establishes the supplied equilibrium; that operation
/// is not a simulated heat interval. Subsequent rejected candidates change none
/// of these histories.
pub struct ReciprocalSpatialSolid<'s, 'a, 'm, 'c> {
    source: &'s SpatialSolidEnthalpy<'a, 'm, 'c>,
    state: SpatialSolidState<'s>,
    equilibrium: SolidEquilibrium,
    laws: Arc<[AffineThermoelastic]>,
    specific_h: Vec<f64>,
    time_s: f64,
    fixed_dofs: Vec<usize>,
    nodal_forces_n: Vec<(usize, f64)>,
    solve: SolidEquilibriumConfig,
    minimum_nodal_capacity_j_k: Arc<[f64]>,
}

/// Every element exchange is in original tetrahedron order and original joules.
#[derive(Debug)]
pub struct SpatialThermoelasticReport {
    pub iterations: usize,
    pub elements: Vec<ThermoelasticExchange>,
    pub nodal_reversible_heat_j: Vec<Dd>,
    pub stress_work_j: f64,
    pub external_work_j: f64,
    pub mechanical_residual_j: f64,
    pub internal_energy_change_j: f64,
    pub reversible_heat_j: f64,
    /// Includes authored volume/boundary heat and prescribed-T reactions.
    pub external_thermal_heat_j: f64,
    pub heat_defect_j: f64,
    pub total_energy_residual_j: f64,
}

/// The native heat endpoint is independently re-evaluated with the final actual
/// material heat. Its Newton report belongs to the final held-source native solve.
pub struct SpatialThermoelasticCandidate<'s, 'a, 'm, 'c> {
    pub solid: ReciprocalSpatialSolid<'s, 'a, 'm, 'c>,
    pub thermal: NodalHeatStep,
    pub report: SpatialThermoelasticReport,
}

impl<'s, 'a: 's, 'm: 'a, 'c: 'a> ReciprocalSpatialSolid<'s, 'a, 'm, 'c> {
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        source: &'s SpatialSolidEnthalpy<'a, 'm, 'c>,
        specific_h: &[f64],
        time_s: f64,
        loads: &TetLoadCase3<'_, '_>,
        fixed_dofs: &[usize],
        solve: SolidEquilibriumConfig,
        config: ThermoelasticSolve,
        _caloric_reference: SpatialCaloricReference,
        cx: &Cx<'_>,
    ) -> Result<Self> {
        check(cx)?;
        admit_controls(config)?;
        let state = source.resolve_at(specific_h, time_s, cx)?;
        let mut laws = Vec::with_capacity(source.sources.len());
        let mut minimum_cp = Vec::with_capacity(source.sources.len());
        for bound in &source.sources {
            check(cx)?;
            let input = &bound.source;
            let phase = input.curve
                .state_at_specific_enthalpy(input.reference_specific_enthalpy_j_kg)
                .map_err(SpatialSolidEnthalpyError::from)?;
            let (elastic, alpha, _) = bound.resolve(
                phase, source.mesh_frame, source.maximum_absolute_mandel_strain,
            )?;
            let law = AffineThermoelastic::new(
                input.point.material, &elastic, alpha.reference(),
                stiffness(input.stiffness_model), expansion(input.expansion_model),
                input.temperature_interval_k, source.maximum_absolute_mandel_strain,
            )?;
            let range = input.temperature_interval_k;
            let knots = input.curve.knots();
            if range[0] < knots[0].temperature_k || range[1] > knots[knots.len()-1].temperature_k {
                return Err(invalid("original caloric chart must cover the complete source interval"));
            }
            let mut cp = f64::INFINITY;
            for (index, pair) in knots.windows(2).enumerate() {
                check(cx)?;
                if pair[1].temperature_k == pair[0].temperature_k {
                    if pair[0].temperature_k >= range[0] && pair[0].temperature_k <= range[1] {
                        return Err(invalid("reciprocal solid interval intersects a latent plateau"));
                    }
                    continue;
                }
                if pair[0].temperature_k >= range[1] || pair[1].temperature_k <= range[0] { continue; }
                if pair.iter().any(|knot| knot.liquid_mass_fraction != 0.0) {
                    return Err(invalid("the complete reciprocal source interval must remain fully solid"));
                }
                cp = cp.min(input.curve.sensible_heat_capacity_bounds_j_kg_k(index)
                    .map_err(SpatialSolidEnthalpyError::from)?[0]);
            }
            if !cp.is_finite() || cp <= 0.0 {
                return Err(invalid("positive original caloric capacity required over the source interval"));
            }
            minimum_cp.push(cp);
            laws.push(law);
        }
        // In h coordinates one centroid element contributes C_el/(16*Cp)
        // times the four-node all-ones matrix. This bound makes its negative
        // part diagonally dominated by the ORIGINAL mass. Native heterogeneous
        // preparation forbids different caloric charts at a shared vertex.
        // Positive elastic capacity and positive surface storage need no credit.
        let mut floors = Vec::with_capacity(specific_h.len());
        for (&mass, &id) in source.thermal.reference_nodal_masses_kg().iter()
            .zip(source.thermal.vertex_material_ids())
        {
            check(cx)?;
            floors.push(dd(mass) * dd(minimum_cp[id]));
        }
        for (element, tet) in source.tetrahedra.iter().enumerate() {
            check(cx)?;
            let sample = &state.samples[element];
            let lower = laws[sample.material_index]
                .minimum_heat_capacity_j_k(sample.reference_volume_m3)?.min(0.0);
            for &node in tet { floors[node] = floors[node] + dd(0.25) * dd(lower); }
        }
        let mut minimum_nodal_capacity_j_k = Vec::with_capacity(floors.len());
        for floor in floors {
            check(cx)?;
            let floor = finite(floor)?.next_down();
            if floor <= 0.0 {
                return Err(invalid("caloric plus centroid elastic capacity must remain positive at every original node"));
            }
            minimum_nodal_capacity_j_k.push(floor);
        }
        let equilibrium = state.solve(loads, fixed_dofs, solve, cx)?;
        for (element, sample) in state.samples.iter().enumerate() {
            check(cx)?;
            let law = laws[sample.material_index].at_materials(
                &sample.elastic, sample.expansion.current(),
            )?;
            let native = &equilibrium.solution.stress.elements()[element];
            let potential = law.state(*native.total_strain_mandel(), sample.reference_volume_m3)?;
            admit_native(&potential, native, source.maximum_absolute_mandel_strain, config)?;
        }
        check(cx)?;
        Ok(Self {
            source, state, equilibrium, laws: laws.into(), specific_h: specific_h.to_vec(),
            time_s, fixed_dofs: fixed_dofs.to_vec(),
            nodal_forces_n: loads.nodal_forces_n().to_vec(), solve,
            minimum_nodal_capacity_j_k: minimum_nodal_capacity_j_k.into(),
        })
    }

    pub fn state(&self) -> &SpatialSolidState<'s> { &self.state }
    pub fn equilibrium(&self) -> &SolidEquilibrium { &self.equilibrium }
    pub fn specific_enthalpy_j_kg(&self) -> &[f64] { &self.specific_h }
    pub const fn time_s(&self) -> f64 { self.time_s }
    /// Conservative fixed-strain capacity admission; not the equilibrium
    /// trajectory's effective heat capacity or an inferred dissipation.
    pub fn minimum_nodal_heat_capacities_j_k(&self) -> &[f64] {
        &self.minimum_nodal_capacity_j_k
    }

    /// Advance the SAME native nodal heat inventory and loaded equilibrium.
    /// Endpoint forces define trapezoidal mechanical work; supports and source
    /// identities remain fixed. Each iteration restarts the accepted h field.
    /// Final acceptance checks every native free heat row, native global heat,
    /// actual nodal/element mechanical work, and combined storage independently.
    #[allow(clippy::too_many_arguments)]
    pub fn advance(
        &self,
        dt_s: f64,
        boundary: &ThermalBoundary,
        interfaces: Option<&ThermalInterfaces>,
        volumetric_source: &ScalarField,
        loads: &TetLoadCase3<'_, '_>,
        thermal_config: EnthalpyStepConfig,
        config: ThermoelasticSolve,
        cx: &Cx<'_>,
    ) -> Result<SpatialThermoelasticCandidate<'s, 'a, 'm, 'c>> {
        self.advance_with_radiation(
            dt_s, boundary, interfaces, volumetric_source, loads,
            thermal_config, config, EnthalpyRadiation::None, cx,
        )
    }

    /// The same joint candidate with an original pointwise radiative boundary.
    /// Emission/reflection and its signed external heat remain with native
    /// conduction; they are never added as a second Robin or nodal source.
    /// Unsupported mean-patch approximations refuse in the native heat owner.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_radiation(
        &self,
        dt_s: f64,
        boundary: &ThermalBoundary,
        interfaces: Option<&ThermalInterfaces>,
        volumetric_source: &ScalarField,
        loads: &TetLoadCase3<'_, '_>,
        thermal_config: EnthalpyStepConfig,
        config: ThermoelasticSolve,
        radiation: EnthalpyRadiation<'_>,
        cx: &Cx<'_>,
    ) -> Result<SpatialThermoelasticCandidate<'s, 'a, 'm, 'c>> {
        self.advance_heat_boundary(
            dt_s, boundary, HeatBoundary::Fixed(interfaces, thermal_config),
            volumetric_source, loads, config, radiation, cx,
        )
    }

    /// Evolve live source-backed contact and this same compatible deformation.
    /// Both independent thermal traces are queried at each native Newton trial;
    /// the actual final contact sources remain on the returned thermal endpoint.
    ///
    /// The contact owner includes its fixed interfaces and retains original
    /// source claims, coordinates, geometry, pressure, gap, finish and medium.
    /// Only its explicitly declared thermal/material coordinates evolve. Fixed
    /// mechanical supports do not introduce an unstated pressure or gap law.
    /// Contact transfer is internal heat and is never added as external input.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_temperature_contacts(
        &self,
        dt_s: f64,
        boundary: &ThermalBoundary,
        contacts: &TemperatureContacts<'_, '_>,
        volumetric_source: &ScalarField,
        loads: &TetLoadCase3<'_, '_>,
        contact_config: ContactEnthalpyNewtonConfig,
        config: ThermoelasticSolve,
        radiation: EnthalpyRadiation<'_>,
        cx: &Cx<'_>,
    ) -> Result<SpatialThermoelasticCandidate<'s, 'a, 'm, 'c>> {
        self.advance_heat_boundary(
            dt_s, boundary, HeatBoundary::Temperature(contacts, contact_config),
            volumetric_source, loads, config, radiation, cx,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_heat_boundary(
        &self,
        dt_s: f64,
        boundary: &ThermalBoundary,
        heat_boundary: HeatBoundary<'_, '_, '_>,
        volumetric_source: &ScalarField,
        loads: &TetLoadCase3<'_, '_>,
        config: ThermoelasticSolve,
        radiation: EnthalpyRadiation<'_>,
        cx: &Cx<'_>,
    ) -> Result<SpatialThermoelasticCandidate<'s, 'a, 'm, 'c>> {
        check(cx)?;
        let thermal_config = heat_boundary.step_config();
        admit_controls(config)?;
        let end_s = self.time_s + dt_s;
        if !dt_s.is_finite() || dt_s <= 0.0 || !end_s.is_finite() || end_s <= self.time_s
            || !thermal_config.energy_tolerance_j.is_finite()
            || thermal_config.energy_tolerance_j < 0.0
            || thermal_config.energy_tolerance_j > config.energy_tolerance_j
        {
            return Err(invalid("representable positive physical interval and compatible native heat tolerance required"));
        }
        let source = self.source;
        let mut applied = vec![Dd::ZERO; self.specific_h.len()];
        let mut defects = [f64::INFINITY; 3];
        for iteration in 1..=config.maximum_iterations {
            check(cx)?;
            let mut thermal = heat_boundary.advance(
                source, boundary, volumetric_source, &self.specific_h, dt_s, &applied, radiation, cx,
            )?;
            let state = source.resolve_at(&thermal.specific_enthalpy_j_kg, end_s, cx)?;
            let equilibrium = state.solve(loads, &self.fixed_dofs, self.solve, cx)?;
            let ElementExchange { elements, nodal_heat, work, internal, heat } = element_exchange(
                source, &self.laws, &self.state, &self.equilibrium.solution.stress,
                &state, &equilibrium.solution.stress, config, cx,
            )?;
            let before = self.equilibrium.solution.displacement.displacement_m();
            let after = equilibrium.solution.displacement.displacement_m();
            let mut external = Dd::ZERO;
            for &(dof, force) in self.nodal_forces_n.iter().chain(loads.nodal_forces_n()) {
                check(cx)?;
                external = external + dd(0.5) * dd(force)
                    * (dd(after[dof / 3][dof % 3]) - dd(before[dof / 3][dof % 3]));
            }
            let mechanical_residual = work - external;
            let mechanical_residual_j = finite(mechanical_residual)?;
            if !within(mechanical_residual, config.energy_tolerance_j)? {
                return Err(invalid("native nodal work and original element stress work disagree"));
            }
            let endpoint = heat_boundary.evaluate(
                source, boundary, volumetric_source, &self.specific_h,
                &thermal.specific_enthalpy_j_kg, dt_s, &nodal_heat, radiation, cx,
            )?;
            let external_heat = endpoint.energy.volumetric_source_j
                + endpoint.energy.dirichlet_in_j - endpoint.energy.neumann_out_j
                - endpoint.energy.robin_out_j - endpoint.energy.radiation_out_j;
            let total = endpoint.energy.stored_energy_change_j + internal - external - external_heat;
            let mut heat_defect = Dd::ZERO;
            for (&actual, &used) in nodal_heat.iter().zip(&applied) {
                heat_defect = heat_defect + (actual - used).abs();
            }
            let heat_defect_j = finite(heat_defect)?;
            let total_energy_residual_j = finite(total)?;
            defects = [heat_defect_j, endpoint.maximum_free_residual_j, total_energy_residual_j];
            if within(heat_defect, config.heat_tolerance_j)?
                && heat_boundary.admits_contact_residual(endpoint.residual_norm_j)
                && endpoint.maximum_free_residual_j <= thermal_config.energy_tolerance_j
                && within(endpoint.energy.residual_j, thermal_config.energy_tolerance_j)?
                && within(total, config.energy_tolerance_j)?
            {
                thermal.endpoint = endpoint;
                let report = SpatialThermoelasticReport {
                    iterations: iteration, elements, nodal_reversible_heat_j: nodal_heat,
                    stress_work_j: finite(work)?, external_work_j: finite(external)?,
                    mechanical_residual_j, internal_energy_change_j: finite(internal)?,
                    reversible_heat_j: finite(heat)?, external_thermal_heat_j: finite(external_heat)?,
                    heat_defect_j, total_energy_residual_j,
                };
                let next = Self {
                    source, state, equilibrium, laws: Arc::clone(&self.laws),
                    specific_h: thermal.specific_enthalpy_j_kg.clone(), time_s: end_s,
                    fixed_dofs: self.fixed_dofs.clone(), nodal_forces_n: loads.nodal_forces_n().to_vec(),
                    solve: self.solve, minimum_nodal_capacity_j_k: Arc::clone(&self.minimum_nodal_capacity_j_k),
                };
                check(cx)?;
                return Ok(SpatialThermoelasticCandidate { solid: next, thermal, report });
            }
            for (used, actual) in applied.iter_mut().zip(nodal_heat) {
                *used = *used + dd(config.heat_relaxation) * (actual - *used);
                finite(*used)?;
            }
        }
        Err(SpatialThermoelasticError::NotConverged {
            iterations: config.maximum_iterations,
            heat_defect_j: defects[0], maximum_nodal_residual_j: defects[1],
            combined_residual_j: defects[2],
        })
    }
}

struct ElementExchange {
    elements: Vec<ThermoelasticExchange>,
    nodal_heat: Vec<Dd>,
    work: Dd,
    internal: Dd,
    heat: Dd,
}

// Both equilibrium and inertia use the same original quadrature, exact total
// strains, source/native point checks and retained finite material increments.
#[allow(clippy::too_many_arguments)]
fn element_exchange(
    source: &SpatialSolidEnthalpy<'_, '_, '_>,
    laws: &[AffineThermoelastic],
    before_state: &SpatialSolidState<'_>,
    before_stress: &fs_solid::linear3::TetThermalStressReport,
    after_state: &SpatialSolidState<'_>,
    after_stress: &fs_solid::linear3::TetThermalStressReport,
    config: ThermoelasticSolve,
    cx: &Cx<'_>,
) -> Result<ElementExchange> {
    let mut nodal_heat = vec![Dd::ZERO; source.thermal.mesh().vertex_count()];
    let mut elements = Vec::with_capacity(source.tetrahedra.len());
    let (mut work, mut internal, mut heat) = (Dd::ZERO, Dd::ZERO, Dd::ZERO);
    for (element, tet) in source.tetrahedra.iter().enumerate() {
        check(cx)?;
        let before = &before_stress.elements()[element];
        let after = &after_stress.elements()[element];
        let sample = &after_state.samples[element];
        let exchange = laws[sample.material_index].exchange(
            before_state.samples[element].phase.temperature_k(),
            *before.total_strain_mandel(),
            sample.phase.temperature_k(),
            *after.total_strain_mandel(),
            sample.reference_volume_m3,
        )?;
        admit_native(&exchange.before, before, source.maximum_absolute_mandel_strain, config)?;
        admit_native(&exchange.after, after, source.maximum_absolute_mandel_strain, config)?;
        work = work + exchange.retained_mechanical_work_j();
        internal = internal + exchange.retained_internal_energy_change_j();
        heat = heat + exchange.retained_reversible_heat_j();
        let quarter = dd(0.25) * exchange.retained_reversible_heat_j();
        for &node in tet { nodal_heat[node] = nodal_heat[node] + quarter; }
        elements.push(exchange);
    }
    Ok(ElementExchange { elements, nodal_heat, work, internal, heat })
}

fn admit_native(
    potential: &ThermoelasticState,
    native: &fs_solid::linear3::TetElementEigenstrainStress,
    maximum_strain: f64,
    config: ThermoelasticSolve,
) -> Result<()> {
    let strain_tolerance = 4096.0 * f64::EPSILON * maximum_strain;
    if potential.elastic_strain_mandel.iter().zip(native.elastic_strain_mandel())
        .any(|(a, b)| (a - b).abs() > strain_tolerance)
        || potential.stress_mandel_pa.iter().zip(native.stress_mandel_pa())
            .any(|(a, b)| (a - b).abs() > config.stress_tolerance_pa)
        || !within(dd(potential.free_energy_j) - dd(native.elastic_energy_j()), config.energy_tolerance_j)?
    {
        return Err(invalid("actual native strain, stress or energy differs from the original source potential"));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn physical_account_gate_retains_the_low_component_at_the_tolerance() {
        assert!(!within(Dd::from_pair(1.0, 1e-30), 1.0).unwrap());
        assert!(within(Dd::from_pair(1.0, -1e-30), 1.0).unwrap());
        assert!(!within(Dd::from_pair(-1.0, -1e-30), 1.0).unwrap());
        assert!(within(Dd::from_pair(-1.0, 1e-30), 1.0).unwrap());
    }
}
