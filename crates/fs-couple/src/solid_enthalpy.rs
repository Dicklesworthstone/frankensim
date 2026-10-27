//! Local material enthalpy drives the existing three-dimensional solid owner.
//!
//! One active thermal control volume is explicitly assigned to each reference
//! P1 tetrahedron. Its original caloric law supplies temperature and phase;
//! source-pinned elasticity and integrated expansion supply the native material
//! and eigenstrain fields. The reference mesh and delivered mass stay fixed.
//!
//! This is a sequential thermal-loading approximation. The default solve is
//! quasistatic; resolved states also drive the existing inertial Newmark owner
//! through `advance_transient_candidate` with refreshed stiffness and fixed
//! reference mass. Expansion and stiffness work are explicit. The thermal
//! audit covers caloric storage and heat; it does not include elastic energy or
//! infer reciprocal thermoelastic heat, plastic flow, phase-support laws, or
//! deformation-dependent thermal conductances. Those require their own closures.
//! Partial melt and excessive displacement/strain refuse the whole candidate.
//! The [`thermoelastic`] path closes reciprocal elastic/finite-heat feedback.
//! The sibling [`plastic`] path couples actual loaded J2 tetrahedra and spatial
//! physical-heat feedback under its separately stated energetic assumptions.

use std::collections::{BTreeMap, BTreeSet};

/// Loaded J2 tetrahedra exchange physical plastic heat with spatial enthalpy.
pub mod plastic;
/// Reciprocal source thermodynamics on the native elastic equilibrium/heat history.
pub mod thermoelastic;
/// One finite source gas exchanges heat with distinct spatial solid cells.
pub mod gas;
/// Accepted local heat histories drive the existing pressure-valve reduction.
pub mod valve;

use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_matdb::QueryPoint;
use fs_material::phase::{
    EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState, PhaseStateError,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, MaterialStatePointError,
    ResolvedMaterialStatePoint,
};
use fs_material::thermal::source::{ThermalMaterialError, ThermalMaterialPoint};
use fs_solid::linear3::{
    TetAssemblyBudget, TetDeformedMesh, TetElasticError, TetElasticMaterial,
    TetLinearElasticProblem, TetMaterialField, TetStaticSolveConfig,
    TetThermalStrainField, TetThermalStrainState,
};
pub use fs_solid::linear3::source::{ThermalExpansionModel, ThermalExpansionSource, ThermalStiffnessModel};
use fs_solid::linear3::source::SolidMaterialError;
use fs_solid::loads3::{LoadError3, TetLoadCase3, TetLoadGeometry3, TetLoadSolution3};
use fs_solid::loads3::transient::{TetTransientRun3, TransientError3};
use fs_solid::loads3::transient::events::{MaterialRebind3, MaterialRebindConfig3};
use fs_solid::loads3::transient::thermal::ThermalTransientStep3;
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
    Transient(TransientError3),
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
impl From<SolidMaterialError> for SolidEnthalpyError {
    fn from(error: SolidMaterialError) -> Self {
        match error {
            SolidMaterialError::Material(error) => Self::Material(error),
            SolidMaterialError::Solid(error) => Self::Solid(error),
        }
    }
}
impl From<LoadError3> for SolidEnthalpyError {
    fn from(e: LoadError3) -> Self { Self::Loads(e) }
}
impl From<TransientError3> for SolidEnthalpyError {
    fn from(e: TransientError3) -> Self {
        if e == TransientError3::Cancelled { Self::Cancelled } else { Self::Transient(e) }
    }
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
    /// Reference density and the chosen model's complete elastic schema.
    /// Current stiffness refreshes the pinned elastic claims without re-querying density.
    pub mechanical_selection: MaterialPropertySelection,
    /// Elastic symmetry and its fixed principal-to-mesh orientation.
    pub stiffness_model: ThermalStiffnessModel,
    pub expansion_selection: MaterialPropertySelection,
    /// Expansion symmetry and fixed material axes, independent of stiffness.
    pub expansion_model: ThermalExpansionModel,
}

struct BoundCell<'a> {
    source: TetCellMaterial<'a>,
    reference_cell: DepositedCell,
    reference_point: QueryPoint,
    reference_elastic: ResolvedMaterialStatePoint,
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
    mesh_frame: ContentHash,
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
    /// The nonzero mesh frame binds every material orientation and must also be
    /// used by subsequent load/stress/vibration consumers.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        thermal: &DepositionSimulation<'_>,
        nodes_m: &'a [[f64; 3]],
        tetrahedra: &'a [[usize; 4]],
        sources: Vec<TetCellMaterial<'a>>,
        budget: TetAssemblyBudget,
        maximum_absolute_mandel_strain: f64,
        mesh_frame: ContentHash,
        cx: &Cx<'_>,
    ) -> Result<Self, SolidEnthalpyError> {
        poll(cx)?;
        if mesh_frame == ContentHash([0; 32]) {
            return Err(SolidEnthalpyError::Invalid("nonzero reference mesh frame required"));
        }
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
            let elastic = source.stiffness_model.resolve_reference(
                source.point.material, &reference_point, source.mechanical_selection.clone())?;
            for property in elastic.properties() {
                source.point.admit_static(property)?;
            }
            let density = elastic.property(DENSITY_PROPERTY)
                .ok_or(SolidEnthalpyError::Invalid("reference density required"))?.value_si();
            let material = source.stiffness_model.material(&elastic, density,
                elastic.identity(), maximum_absolute_mandel_strain)?;
            let expansion = source.expansion_model.resolve(
                source.point.material, &reference_point, &reference_point,
                source.expansion_selection.clone(), &material, mesh_frame,
                maximum_absolute_mandel_strain)?;
            for property in expansion.current().properties() {
                source.point.admit_static(property)?;
            }
            let volume = geometry.element_volumes_m3()[element];
            if !same_measure(cell.geometry.volume_m3, volume)
                || !same_measure(cell.mass_kg, density*volume) {
                return Err(SolidEnthalpyError::Invalid(
                    "thermal inventory must match the tetrahedron reference volume and source density"));
            }
            source.mechanical_selection = MaterialPropertySelection::PinnedByProperty(
                elastic.properties().iter().filter(|p| p.requirement().name() != DENSITY_PROPERTY).map(|p| (
                    p.requirement().name().to_owned(), p.answer().receipt.selected,
                )).collect());
            source.expansion_selection = MaterialPropertySelection::PinnedByProperty(
                expansion.reference().properties().iter().map(|p| (
                    p.requirement().name().to_owned(), p.answer().receipt.selected,
                )).collect());
            cells.push(BoundCell { source, reference_cell: cell.clone(), reference_point,
                reference_elastic: elastic, reference_volume_m3: volume });
        }
        poll(cx)?;
        Ok(Self { nodes_m, tetrahedra, geometry, cells, budget, mesh_frame, maximum_absolute_mandel_strain })
    }

    /// Use the retained geometry for native pressure/traction/body load assembly.
    pub const fn geometry(&self) -> &TetLoadGeometry3<'a> { &self.geometry }
    /// Coordinate frame in which the reference nodes and all admitted tensors live.
    pub const fn mesh_frame(&self) -> ContentHash { self.mesh_frame }

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
            let elastic = source.stiffness_model.resolve(
                source.point.material, &point, source.mechanical_selection.clone())?;
            for property in elastic.properties() { source.point.admit_static(property)?; }

            // Inertia and weight use the exact delivered reference inventory.
            // No current density is required and no heating step creates matter.
            let reference_density = reference.mass_kg/bound.reference_volume_m3;
            let mut identity = DomainHasher::new("org.frankensim.fs-couple.solid-enthalpy-material.v2");
            identity.update(elastic.identity().as_bytes());
            identity.update(bound.reference_elastic.identity().as_bytes());
            identity.update(phase.identity().as_bytes());
            identity.update(&reference.id.0.to_le_bytes());
            identity.update(&reference.mass_kg.to_bits().to_le_bytes());
            identity.update(&bound.reference_volume_m3.to_bits().to_le_bytes());
            identity.update(self.mesh_frame.as_bytes());
            let material = source.stiffness_model.material(&elastic,
                reference_density, identity.finalize(), self.maximum_absolute_mandel_strain)?;
            let expansion = source.expansion_model.resolve(
                source.point.material, &bound.reference_point, &point,
                source.expansion_selection.clone(), &material, self.mesh_frame,
                self.maximum_absolute_mandel_strain)?;
            for property in expansion.current().properties() { source.point.admit_static(property)?; }
            strains.push(*expansion.thermal_strain());
            materials.push(material);
            samples.push(TetMaterialSample { cell: source.cell, phase, elastic, expansion,
                reference_elastic: bound.reference_elastic.clone(),
                stiffness_model: source.stiffness_model, expansion_model: source.expansion_model,
                reference_mass_kg: reference.mass_kg, reference_volume_m3: bound.reference_volume_m3 });
        }
        poll(cx)?;
        Ok(SolidEnthalpyState { nodes_m: self.nodes_m, tetrahedra: self.tetrahedra,
            budget: self.budget, mesh_frame: self.mesh_frame, time_s: thermal.time_s(), materials, strains, samples,
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
    /// Current elastic constants; density belongs to the original inventory.
    pub elastic: ResolvedMaterialStatePoint,
    /// Original density and elastic claims at the stress-free reference state.
    pub reference_elastic: ResolvedMaterialStatePoint,
    pub stiffness_model: ThermalStiffnessModel,
    pub expansion: ThermalExpansionSource,
    pub expansion_model: ThermalExpansionModel,
    pub reference_mass_kg: f64,
    pub reference_volume_m3: f64,
}

/// One explicitly replaced cell; paired values are ordered [before, after].
/// The two caloric energies retain their OWN chart references. Their difference
/// is not declared heat or work: replacement may exchange matter, and distinct
/// materials need not share an enthalpy zero or chemical-energy convention.
#[derive(Clone, Debug, PartialEq)]
pub struct SolidMaterialRebindCell {
    pub cell: CellId,
    pub reference_mass_kg: [f64; 2],
    /// Complete local h/T/phase/card/chart bindings, without a second heat state.
    pub phase: [EquilibriumPhaseState; 2],
    /// The extensive m*h inventories on their respective admitted charts [J].
    pub caloric_chart_energy_j: [f64; 2],
}

/// Native mechanical transfer plus the actual exchanged material inventories.
/// The external replacement owner must account for matter and its chosen
/// caloric/chemical reference; this report does not silently debit either heat
/// graph, infer a phase transition, or treat mechanical replacement work as heat.
#[derive(Clone, Debug, PartialEq)]
pub struct SolidMaterialRebind {
    pub mechanical: MaterialRebind3,
    /// Original tetrahedron order, with unchanged persistent cell IDs.
    pub cells: Vec<SolidMaterialRebindCell>,
}

/// A complete local state ready for the native load, support-motion, stress,
/// prestress/vibration and transient-solid APIs. Those owners remain unchanged.
#[derive(Clone, Debug)]
pub struct SolidEnthalpyState<'a> {
    nodes_m: &'a [[f64; 3]],
    tetrahedra: &'a [[usize; 4]],
    budget: TetAssemblyBudget,
    mesh_frame: ContentHash,
    time_s: f64,
    materials: Vec<TetElasticMaterial>,
    strains: Vec<TetThermalStrainState>,
    samples: Vec<TetMaterialSample>,
    maximum_absolute_mandel_strain: f64,
}
impl SolidEnthalpyState<'_> {
    pub const fn time_s(&self) -> f64 { self.time_s }
    pub const fn mesh_frame(&self) -> ContentHash { self.mesh_frame }
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

    fn matches_problem(&self, problem: &TetLinearElasticProblem<'_>) -> bool {
        self.nodes_m == problem.nodes_m && self.tetrahedra == problem.tetrahedra
            && match problem.materials {
                TetMaterialField::Uniform(material) => self.materials.iter().all(|m| m == material),
                TetMaterialField::PerElement(materials) => materials == self.materials,
            }
    }

    fn matches_transient(&self, run: &TetTransientRun3<'_, '_, '_, '_>) -> bool {
        let thermal_matches = match run.thermal_strain_field() {
            Some(TetThermalStrainField::Uniform(state)) => self.strains.iter().all(|s| s == state),
            Some(TetThermalStrainField::PerElement(states)) => states == self.strains,
            None => false,
        };
        run.state().time_s.to_bits() == self.time_s.to_bits()
            && self.matches_problem(run.material_problem()) && thermal_matches
    }

    /// Replace a running specimen with this source-resolved material state at
    /// the SAME physical time. Construct both states through their actual
    /// `SolidEnthalpy`/`DepositionSimulation` owners; supply this state's
    /// `problem()` and an explicitly assembled replacement load, including any
    /// changed weight. No density, elastic constant or thermal strain is copied
    /// into a separate coefficient API.
    ///
    /// Displacement, time and step count survive. The native rebind solves the
    /// replacement velocity from retained reduced canonical momentum and then
    /// equilibrates acceleration with the new mass/stiffness/thermal/load data.
    /// Reference mesh, zero supports, frame and cell correspondence stay fixed.
    /// The new binding supplies its own stress-free reference and caloric chart.
    ///
    /// This is an explicit external replacement, not the ordinary heating path
    /// or a melting/mass-flow closure. Mechanical work is bounded by `controls`;
    /// chart-specific mass/enthalpy inventories are reported separately because
    /// unrelated caloric zeros cannot be subtracted as physical energy. Accept
    /// the new thermal owner, material binding, resolved state and returned run
    /// together. Every input remains unchanged on any refusal.
    #[allow(clippy::too_many_arguments)]
    pub fn rebind_transient_candidate<'g, 'm, 'p, 's>(
        &self, previous: &Self, run: &TetTransientRun3<'g, 'm, '_, '_>,
        problem: &'p TetLinearElasticProblem<'s>, loads: &TetLoadCase3<'g, 'm>,
        controls: MaterialRebindConfig3, cx: &Cx<'_>,
    ) -> Result<(TetTransientRun3<'g, 'm, 'p, 's>, SolidMaterialRebind), SolidEnthalpyError> {
        poll(cx)?;
        if !previous.matches_transient(run) || !self.matches_problem(problem)
            || self.time_s.to_bits() != previous.time_s.to_bits()
            || self.mesh_frame != previous.mesh_frame
            || self.nodes_m != previous.nodes_m || self.tetrahedra != previous.tetrahedra
            || self.samples.len() != previous.samples.len()
            || !controls.maximum_strain.is_finite()
            || controls.maximum_strain > self.maximum_absolute_mandel_strain
            || controls.maximum_strain > previous.maximum_absolute_mandel_strain {
            return Err(SolidEnthalpyError::Invalid(
                "material replacement must bind the accepted motion and the exact new sources at one clock, mesh, frame and admitted strain limit"));
        }
        let mut cells = Vec::with_capacity(self.samples.len());
        for (old, next) in previous.samples.iter().zip(&self.samples) {
            poll(cx)?;
            if old.cell != next.cell
                || old.reference_volume_m3.to_bits() != next.reference_volume_m3.to_bits() {
                return Err(SolidEnthalpyError::Invalid(
                    "material replacement must retain each reference tetrahedron and cell correspondence"));
            }
            let energies = [old.reference_mass_kg*old.phase.specific_enthalpy_j_kg(),
                next.reference_mass_kg*next.phase.specific_enthalpy_j_kg()];
            if !energies.iter().all(|e| e.is_finite()) {
                return Err(SolidEnthalpyError::Invalid(
                    "replacement caloric inventories must be representable on both original charts"));
            }
            cells.push(SolidMaterialRebindCell { cell: old.cell,
                reference_mass_kg: [old.reference_mass_kg,next.reference_mass_kg],
                phase: [old.phase,next.phase], caloric_chart_energy_j: energies });
        }
        let (candidate, mechanical) = run.rebind_material_candidate(problem,
            Some(self.thermal_strains()), loads, self.mesh_frame, controls, cx)?;
        poll(cx)?;
        Ok((candidate, SolidMaterialRebind { mechanical, cells }))
    }

    /// Use actual local enthalpy temperatures to drive the SAME inertial state.
    /// `previous` must exactly match the accepted run's material/thermal fields
    /// and clock. `problem` must come from this endpoint's `problem()`; it is
    /// borrowed by the returned run so refreshed source identities remain alive.
    /// Keep both resolved states and their problems until the candidate is accepted.
    ///
    /// The native owner retains reference mass, mesh and zero supports while
    /// rebuilding changed stiffness. Motion/history continues, and expansion
    /// and stiffness work are audited separately, even in fully clamped elements.
    /// Changing stiffness admits mass-proportional damping only. This is
    /// prescribed one-way thermal loading: signed material work is NOT subtracted
    /// from caloric storage. Reciprocal thermoelastic heat needs its own closure.
    /// This endpoint adapter does not inspect intervening thermal events; admit
    /// that path, resolve heating events and verify temporal refinement explicitly.
    pub fn advance_transient_candidate<'g, 'm, 'p, 's>(
        &self, previous: &Self, run: &TetTransientRun3<'g, 'm, '_, '_>,
        problem: &'p TetLinearElasticProblem<'s>, loads: &TetLoadCase3<'_, '_>, cx: &Cx<'_>,
    ) -> Result<(TetTransientRun3<'g, 'm, 'p, 's>, ThermalTransientStep3), SolidEnthalpyError> {
        poll(cx)?;
        if !previous.matches_transient(run)
            || !self.matches_problem(problem) || self.mesh_frame != previous.mesh_frame
            || self.maximum_absolute_mandel_strain.to_bits()
                != previous.maximum_absolute_mandel_strain.to_bits()
            || self.samples.len() != previous.samples.len() {
            return Err(SolidEnthalpyError::Invalid("transient endpoints must retain their exact material, thermal, frame and clock bindings"));
        }
        for (old, next) in previous.samples.iter().zip(&self.samples) {
            poll(cx)?;
            if old.cell != next.cell || old.reference_elastic != next.reference_elastic
                || old.expansion.reference() != next.expansion.reference()
                || old.stiffness_model != next.stiffness_model || old.expansion_model != next.expansion_model
                || old.reference_mass_kg.to_bits() != next.reference_mass_kg.to_bits()
                || old.reference_volume_m3.to_bits() != next.reference_volume_m3.to_bits() {
                return Err(SolidEnthalpyError::Invalid("transient heating must preserve the original specimen, sources and reference inventory"));
            }
        }
        let h_s = self.time_s - previous.time_s;
        if !h_s.is_finite() || h_s <= 0.0
            || (run.state().time_s + h_s).to_bits() != self.time_s.to_bits() {
            return Err(SolidEnthalpyError::Invalid("thermal and mechanical endpoints must share one distinct representable clock"));
        }
        Ok(run.advance_thermal_candidate(problem, self.thermal_strains(), loads, h_s,
            self.mesh_frame, self.maximum_absolute_mandel_strain, cx)?)
    }

    pub fn solve(
        &self, loads: &TetLoadCase3<'_, '_>, fixed_dofs: &[usize],
        config: SolidEquilibriumConfig, cx: &Cx<'_>,
    ) -> Result<SolidEquilibrium, SolidEnthalpyError> {
        poll(cx)?;
        if config.mesh_frame != self.mesh_frame {
            return Err(SolidEnthalpyError::Invalid(
                "equilibrium frame must match the reference mesh and admitted material tensors"));
        }
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
