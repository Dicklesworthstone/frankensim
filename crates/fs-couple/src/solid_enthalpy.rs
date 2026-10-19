//! Local material enthalpy drives the existing three-dimensional solid owner.
//!
//! One active thermal control volume is explicitly assigned to each reference
//! P1 tetrahedron. Its original caloric law supplies temperature and phase;
//! source-pinned elasticity and integrated expansion supply the native material
//! and eigenstrain fields. The reference mesh and delivered mass stay fixed.
//!
//! This is a sequential, quasistatic thermal-loading approximation. The thermal
//! audit covers caloric storage and heat; it does not include elastic energy or
//! infer reciprocal thermoelastic heat, plastic flow, phase-support laws, or
//! deformation-dependent thermal conductances. Those require their own closures.
//! Partial melt and excessive displacement/strain refuse the whole candidate.
//! The sibling [`plastic`] path couples actual loaded J2 tetrahedra and spatial
//! physical-heat feedback under its separately stated energetic assumptions.

use std::collections::{BTreeMap, BTreeSet};

/// Loaded J2 tetrahedra exchange physical plastic heat with spatial enthalpy.
pub mod plastic;

use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_matdb::QueryPoint;
use fs_material::phase::{
    EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState, PhaseStateError,
};
use fs_material::state_point::{
    IntegratedIsotropicThermalExpansion, IsotropicElasticStatePoint,
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, MaterialPropertySelection,
    MaterialStatePointError, integrate_isotropic_thermal_expansion,
    resolve_isotropic_elastic_state_point,
};
use fs_material::thermal::source::{ThermalMaterialError, ThermalMaterialPoint};
use fs_solid::linear3::{
    TetAssemblyBudget, TetDeformedMesh, TetElasticError, TetElasticMaterial,
    TetLinearElasticProblem, TetMaterialField, TetStaticSolveConfig,
    TetThermalStrainField, TetThermalStrainState,
};
use fs_solid::loads3::{LoadError3, TetLoadCase3, TetLoadGeometry3, TetLoadSolution3};
use fs_thermal::enthalpy::{
    CellId, DepositedCell, DepositionError, DepositionObserver, DepositionSimulation,
    EnergyAudit, Event, StepConfig, ThermalLaw,
};

#[derive(Debug)]
pub enum SolidEnthalpyError {
    Invalid(&'static str),
    Thermal(DepositionError),
    Phase(PhaseStateError),
    Material(MaterialStatePointError),
    Source(ThermalMaterialError),
    Solid(TetElasticError),
    Loads(LoadError3),
    Cancelled,
}
impl std::fmt::Display for SolidEnthalpyError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "material/enthalpy solid coupling refused: {self:?}")
    }
}
impl std::error::Error for SolidEnthalpyError {}
impl From<DepositionError> for SolidEnthalpyError {
    fn from(e: DepositionError) -> Self { Self::Thermal(e) }
}
impl From<PhaseStateError> for SolidEnthalpyError {
    fn from(e: PhaseStateError) -> Self { Self::Phase(e) }
}
impl From<MaterialStatePointError> for SolidEnthalpyError {
    fn from(e: MaterialStatePointError) -> Self { Self::Material(e) }
}
impl From<ThermalMaterialError> for SolidEnthalpyError {
    fn from(e: ThermalMaterialError) -> Self { Self::Source(e) }
}
impl From<TetElasticError> for SolidEnthalpyError {
    fn from(e: TetElasticError) -> Self { Self::Solid(e) }
}
impl From<LoadError3> for SolidEnthalpyError {
    fn from(e: LoadError3) -> Self { Self::Loads(e) }
}
fn poll(cx: &Cx<'_>) -> Result<(), SolidEnthalpyError> {
    cx.checkpoint().map_err(|_| SolidEnthalpyError::Cancelled)
}
fn same_measure(a: f64, b: f64) -> bool {
    a.is_finite() && b.is_finite() && a > 0.0 && b > 0.0
        && (a-b).abs() <= 512.0*f64::EPSILON*a.abs().max(b.abs())
}

/// Ordered like tetrahedron connectivity; cell IDs are never sorted into mesh
/// order. Different elements can use different cards and caloric charts.
#[derive(Clone, Debug)]
pub struct TetCellMaterial<'a> {
    pub cell: CellId,
    pub curve: &'a EquilibriumEnthalpyPhaseCurve,
    pub point: ThermalMaterialPoint<'a>,
    /// Requires density, Young modulus and Poisson ratio at one common point.
    pub mechanical_selection: MaterialPropertySelection,
    pub expansion_selection: MaterialPropertySelection,
}

struct BoundCell<'a> {
    source: TetCellMaterial<'a>,
    reference_cell: DepositedCell,
    reference_point: QueryPoint,
    reference_elastic: IsotropicElasticStatePoint,
    reference_volume_m3: f64,
}

/// A fixed, explicitly mapped reference specimen. The existing load-geometry
/// owner admits the mesh and supplies volumes; no second geometry kernel exists.
pub struct SolidEnthalpy<'a> {
    nodes_m: &'a [[f64; 3]],
    tetrahedra: &'a [[usize; 4]],
    geometry: TetLoadGeometry3<'a>,
    cells: Vec<BoundCell<'a>>,
    budget: TetAssemblyBudget,
    maximum_absolute_mandel_strain: f64,
}

impl<'a> SolidEnthalpy<'a> {
    /// Pin every material claim at the accepted initial thermal state, which is
    /// the stress-free reference. Each cell volume must match its assigned tet;
    /// delivered mass must match reference material density times that volume
    /// within roundoff. Extra thermal cells may be reservoirs or other bodies.
    ///
    /// The literal T axis is required by the existing integrated-expansion law.
    /// External pressure/composition/processing coordinates remain unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        thermal: &DepositionSimulation<'_>,
        nodes_m: &'a [[f64; 3]],
        tetrahedra: &'a [[usize; 4]],
        sources: Vec<TetCellMaterial<'a>>,
        budget: TetAssemblyBudget,
        maximum_absolute_mandel_strain: f64,
        cx: &Cx<'_>,
    ) -> Result<Self, SolidEnthalpyError> {
        poll(cx)?;
        if sources.len() != tetrahedra.len()
            || !maximum_absolute_mandel_strain.is_finite()
            || maximum_absolute_mandel_strain <= 0.0
            || maximum_absolute_mandel_strain >= 1.0 {
            return Err(SolidEnthalpyError::Invalid(
                "one material cell per tet and a small-strain bound in (0,1) required"));
        }
        let entries = tetrahedra.len().checked_mul(4)
            .ok_or(SolidEnthalpyError::Invalid("face count overflow"))?;
        let geometry = TetLoadGeometry3::new(nodes_m, tetrahedra, budget, entries, cx)?;
        let mut used = BTreeSet::new();
        let mut cells = Vec::with_capacity(sources.len());
        for (element, mut source) in sources.into_iter().enumerate() {
            poll(cx)?;
            if !used.insert(source.cell) || source.point.temperature_axis != "T" {
                return Err(SolidEnthalpyError::Invalid(
                    "unique active cell bindings and the canonical T coordinate required"));
            }
            let cell = thermal.cells().get(&source.cell)
                .ok_or(SolidEnthalpyError::Invalid("structural cell is inactive"))?;
            let phase = phase_at(thermal, cell, &source)?;
            let reference_point = source.point.at(phase)?;
            let elastic = resolve_isotropic_elastic_state_point(
                source.point.material, &reference_point, source.mechanical_selection.clone())?;
            for property in elastic.resolved().properties() {
                source.point.admit_static(property)?;
            }
            let expansion = integrate_isotropic_thermal_expansion(
                source.point.material, &reference_point, &reference_point,
                source.expansion_selection.clone())?;
            for property in expansion.current().resolved().properties() {
                source.point.admit_static(property)?;
            }
            TetThermalStrainState::try_from_isotropic_expansion(
                &elastic, &expansion, maximum_absolute_mandel_strain)?;
            let volume = geometry.element_volumes_m3()[element];
            if !same_measure(cell.geometry.volume_m3, volume)
                || !same_measure(cell.mass_kg, elastic.density_kg_m3()*volume) {
                return Err(SolidEnthalpyError::Invalid(
                    "thermal inventory must match the tetrahedron reference volume and source density"));
            }
            source.mechanical_selection = MaterialPropertySelection::PinnedByProperty(
                elastic.resolved().properties().iter().map(|p| (
                    p.requirement().name().to_owned(), p.answer().receipt.selected,
                )).collect());
            source.expansion_selection = MaterialPropertySelection::PinnedByProperty(vec![(
                LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY.to_owned(), expansion.selected_claim(),
            )]);
            cells.push(BoundCell { source, reference_cell: cell.clone(), reference_point,
                reference_elastic: elastic, reference_volume_m3: volume });
        }
        poll(cx)?;
        Ok(Self { nodes_m, tetrahedra, geometry, cells, budget, maximum_absolute_mandel_strain })
    }

    /// Use the retained geometry for native pressure/traction/body load assembly.
    pub const fn geometry(&self) -> &TetLoadGeometry3<'a> { &self.geometry }

    /// Refresh every local source from the exact current thermal inventory.
    /// Expansion always integrates from the original reference, never from the
    /// previous observation. Source support cannot silently switch to another
    /// same-name claim as temperature evolves.
    pub fn resolve(
        &self, thermal: &DepositionSimulation<'_>, cx: &Cx<'_>,
    ) -> Result<SolidEnthalpyState<'a>, SolidEnthalpyError> {
        poll(cx)?;
        let mut materials = Vec::with_capacity(self.cells.len());
        let mut strains = Vec::with_capacity(self.cells.len());
        let mut samples = Vec::with_capacity(self.cells.len());
        for bound in &self.cells {
            poll(cx)?;
            let source = &bound.source;
            let cell = thermal.cells().get(&source.cell)
                .ok_or(SolidEnthalpyError::Invalid("structural cell is inactive"))?;
            let reference = &bound.reference_cell;
            if cell.id != reference.id || cell.material_id != reference.material_id
                || cell.segment_id != reference.segment_id || cell.body_id != reference.body_id
                || cell.geometry != reference.geometry
                || cell.mass_kg.to_bits() != reference.mass_kg.to_bits() {
                return Err(SolidEnthalpyError::Invalid(
                    "structural material assignment, geometry and reference mass must remain unchanged"));
            }
            let phase = phase_at(thermal, cell, source)?;
            let point = source.point.at(phase)?;
            let elastic = resolve_isotropic_elastic_state_point(
                source.point.material, &point, source.mechanical_selection.clone())?;
            for property in elastic.resolved().properties() { source.point.admit_static(property)?; }
            let expansion = integrate_isotropic_thermal_expansion(
                source.point.material, &bound.reference_point, &point,
                source.expansion_selection.clone())?;
            for property in expansion.current().resolved().properties() {
                source.point.admit_static(property)?;
            }
            let free = TetThermalStrainState::try_from_isotropic_expansion(
                &elastic, &expansion, self.maximum_absolute_mandel_strain)?;

            // Inertia and weight use delivered reference mass. Current rho(T)
            // remains in the source receipt; it must not create/remove matter on
            // the unchanged reference mesh. Bind this lowering to BOTH states.
            let reference_density = reference.mass_kg/bound.reference_volume_m3;
            let mut identity = DomainHasher::new("org.frankensim.fs-couple.solid-enthalpy-material.v1");
            identity.update(elastic.resolved().identity().as_bytes());
            identity.update(bound.reference_elastic.resolved().identity().as_bytes());
            identity.update(phase.identity().as_bytes());
            identity.update(&reference.id.0.to_le_bytes());
            identity.update(&reference.mass_kg.to_bits().to_le_bytes());
            identity.update(&bound.reference_volume_m3.to_bits().to_le_bytes());
            let identity = *identity.finalize().as_bytes();
            materials.push(TetElasticMaterial::try_new(reference_density,
                elastic.young_modulus_pa(), elastic.poisson_ratio(), identity)?);
            strains.push(TetThermalStrainState::try_new(
                free.temperature_k(), free.reference_temperature_k(), free.free_strain_mandel(),
                identity, free.thermal_law_identity(), self.maximum_absolute_mandel_strain)?);
            samples.push(TetMaterialSample { cell: source.cell, phase, elastic, expansion,
                reference_mass_kg: reference.mass_kg, reference_volume_m3: bound.reference_volume_m3 });
        }
        poll(cx)?;
        Ok(SolidEnthalpyState { nodes_m: self.nodes_m, tetrahedra: self.tetrahedra,
            budget: self.budget, time_s: thermal.time_s(), materials, strains, samples,
            maximum_absolute_mandel_strain: self.maximum_absolute_mandel_strain })
    }

    /// Advance the existing heat timeline and recover compatible structural
    /// equilibrium as one immutable candidate. Any thermal, material, solve,
    /// small-strain, geometric or cancellation refusal publishes neither state.
    /// The existing substep/event observer checks source and solid-phase support
    /// throughout the actual discrete thermal path, including heat/cool cycles.
    ///
    /// The thermal audit is the original CALORIC balance. The returned native
    /// stress report records elastic energy separately, without claiming that
    /// this sequential approximation supplies reversible thermoelastic feedback.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_to<'t>(
        &self, thermal: &DepositionSimulation<'t>, target_s: f64, thermal_config: StepConfig,
        loads: &TetLoadCase3<'_, '_>, fixed_dofs: &[usize], solve: SolidEquilibriumConfig,
        cx: &Cx<'_>,
    ) -> Result<SolidEnthalpyCandidate<'a, 't>, SolidEnthalpyError> {
        poll(cx)?;
        // Validate the accepted starting inventory before any proposed aging.
        self.resolve(thermal, cx)?;
        let observer = SolidSourceObserver { source: self, cx };
        let (thermal, _, thermal_audit) = thermal.advance_observed(
            &observer, target_s, thermal_config, ||
                cx.checkpoint().map_err(|_| DepositionError::Cancelled))?;
        let state = self.resolve(&thermal, cx)?;
        let equilibrium = state.solve(loads, fixed_dofs, solve, cx)?;
        poll(cx)?;
        Ok(SolidEnthalpyCandidate { thermal, thermal_audit, state, equilibrium })
    }
}

// The thermal owner already visits every accepted substep and individual event.
// Reuse that path to prevent unobserved melting/source-support exits followed by
// re-solidification/re-entry. No second heat integrator or mechanical solver runs.
#[derive(Clone, Copy)]
struct SolidSourceObserver<'s, 'm, 'cx> {
    source: &'s SolidEnthalpy<'m>,
    cx: &'s Cx<'cx>,
}
impl DepositionObserver for SolidSourceObserver<'_, '_, '_> {
    type Error = SolidEnthalpyError;
    fn event(&mut self, _: &Event, thermal: &DepositionSimulation<'_>)
        -> Result<(), Self::Error> {
        self.source.resolve(thermal, self.cx).map(|_| ())
    }
    fn step(&mut self, _: f64, _: f64, _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        checkpoint()?;
        self.source.resolve(thermal, self.cx).map(|_| ())
    }
}

fn phase_at(
    thermal: &DepositionSimulation<'_>, cell: &DepositedCell, source: &TetCellMaterial<'_>,
) -> Result<EquilibriumPhaseState, SolidEnthalpyError> {
    if !thermal.uses_material_law(cell.material_id, source.curve as &dyn ThermalLaw) {
        return Err(SolidEnthalpyError::Invalid("cell must retain its actual borrowed caloric law"));
    }
    let phase = source.curve.state_at_specific_enthalpy(cell.specific_h_j_kg)?;
    if !phase.is_fully_solid() {
        return Err(SolidEnthalpyError::Invalid(
            "solid-liquid support/flow closure required before partially molten material can carry elastic stress"));
    }
    // The exact law instance is checked above; its temperature is not a second
    // independently supplied mechanical temperature.
    Ok(phase)
}

/// Original constitutive receipts plus the independent reference inventory.
#[derive(Clone, Debug)]
pub struct TetMaterialSample {
    pub cell: CellId,
    pub phase: EquilibriumPhaseState,
    pub elastic: IsotropicElasticStatePoint,
    pub expansion: IntegratedIsotropicThermalExpansion,
    pub reference_mass_kg: f64,
    pub reference_volume_m3: f64,
}

/// A complete local state ready for the native load, support-motion, stress,
/// prestress/vibration and transient-solid APIs. Those owners remain unchanged.
#[derive(Clone, Debug)]
pub struct SolidEnthalpyState<'a> {
    nodes_m: &'a [[f64; 3]],
    tetrahedra: &'a [[usize; 4]],
    budget: TetAssemblyBudget,
    time_s: f64,
    materials: Vec<TetElasticMaterial>,
    strains: Vec<TetThermalStrainState>,
    samples: Vec<TetMaterialSample>,
    maximum_absolute_mandel_strain: f64,
}
impl SolidEnthalpyState<'_> {
    pub const fn time_s(&self) -> f64 { self.time_s }
    pub fn samples(&self) -> &[TetMaterialSample] { &self.samples }
    pub fn materials(&self) -> &[TetElasticMaterial] { &self.materials }
    pub fn thermal_strains(&self) -> TetThermalStrainField<'_> {
        TetThermalStrainField::PerElement(&self.strains)
    }
    /// Construct the original solver input with the SAME refreshed material
    /// fields. Reference density also reaches acceleration loads and inertia.
    pub fn problem<'s>(&'s self, fixed_dofs: &'s [usize]) -> TetLinearElasticProblem<'s> {
        TetLinearElasticProblem { nodes_m: self.nodes_m, tetrahedra: self.tetrahedra,
            materials: TetMaterialField::PerElement(&self.materials),
            fixed_dofs, budget: self.budget }
    }
    pub fn solve(
        &self, loads: &TetLoadCase3<'_, '_>, fixed_dofs: &[usize],
        config: SolidEquilibriumConfig, cx: &Cx<'_>,
    ) -> Result<SolidEquilibrium, SolidEnthalpyError> {
        poll(cx)?;
        let problem = self.problem(fixed_dofs);
        let solution = loads.solve(&problem, self.thermal_strains(), config.static_solve,
            config.mesh_frame, self.maximum_absolute_mandel_strain, cx)?;
        let deformed_mesh = problem.update_geometry_from_displacement(
            &solution.displacement, config.maximum_nodal_displacement_m, cx)?;
        poll(cx)?;
        Ok(SolidEquilibrium { solution, deformed_mesh })
    }
}

#[derive(Clone, Copy, Debug)]
pub struct SolidEquilibriumConfig {
    pub static_solve: TetStaticSolveConfig,
    pub mesh_frame: ContentHash,
    pub maximum_nodal_displacement_m: f64,
}
#[derive(Debug)]
pub struct SolidEquilibrium {
    pub solution: TetLoadSolution3,
    pub deformed_mesh: TetDeformedMesh,
}
pub struct SolidEnthalpyCandidate<'m, 't> {
    pub thermal: DepositionSimulation<'t>,
    pub thermal_audit: EnergyAudit,
    pub state: SolidEnthalpyState<'m>,
    pub equilibrium: SolidEquilibrium,
}

