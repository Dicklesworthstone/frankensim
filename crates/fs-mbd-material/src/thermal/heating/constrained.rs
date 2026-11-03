//! One source-bound heat timeline and compatible constrained FEM configuration.
//!
//! The original thermal owner evolves enthalpy; the original solid owner solves
//! quasistatic displacement/stress at every accepted thermal substep. Geometry
//! and reference mass use the same tetrahedral inventory, not independent free
//! expansion. This is one-way thermal-to-mechanical coupling: elastic energy
//! is reported separately, never mistaken for a total thermomechanical balance.

use std::collections::BTreeMap;
use fs_exec::Cx;
use fs_material::state_point::MaterialPropertySelection;
use fs_thermal::enthalpy::{CellId, DepositionBudget, DepositionError, DepositionObserver,
    DepositionSimulation, EnergyAudit, Event, StepConfig, TimedEvent};
use fs_thermal::observers::{CoupledObservers, ObserverCouplingError};
use crate::deformation::TetrahedralMass;
use crate::deformation::thermoelastic::{ThermalBody, ThermalBodyError, ThermalCell, ThermalSolveConfig, ThermalStiffnessModel, ThermalExpansionModel};
use super::{HeatingError, ThermalRegion, initialize};

/// One original thermal/source description and its current stiffness selection.
#[derive(Clone, Debug)]
pub struct ConstrainedRegion<'a> {
    /// Reuse the existing storage, original card, caloric chart, alpha selection,
    /// physical cell IDs, boundary and reference heat geometry. The descriptor's
    /// `fixed_point_m` is NOT a FEM constraint: it applies only to the separate
    /// free-expansion API. Here supports determine compatible displacement.
    pub storage: ThermalRegion<'a>,
    /// Current elastic-constant claim selection; density stays reference-bound.
    pub stiffness_selection: MaterialPropertySelection,
    /// Same material symmetry and axes used in every accepted thermal substep.
    pub stiffness_model: ThermalStiffnessModel,
    /// Same thermal symmetry and fixed principal axes throughout this history.
    pub expansion_model: ThermalExpansionModel,
}

/// Original failing owner; no partly updated thermal or geometric state escapes.
#[derive(Clone, Debug, PartialEq)]
pub enum ConstrainedHeatError {
    /// Source/storage admission, enthalpy-domain or thermal step failure.
    Heating(HeatingError),
    /// Constitutive, equilibrium, strain, geometry or conserved-inertia failure.
    Shape(ThermalBodyError),
}
impl From<HeatingError> for ConstrainedHeatError {
    fn from(e: HeatingError) -> Self { Self::Heating(e) }
}
impl From<DepositionError> for ConstrainedHeatError {
    fn from(e: DepositionError) -> Self { Self::Heating(e.into()) }
}
impl From<ThermalBodyError> for ConstrainedHeatError {
    fn from(e: ThermalBodyError) -> Self { Self::Shape(e) }
}
impl core::fmt::Display for ConstrainedHeatError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self { Self::Heating(e) => write!(f, "constrained heating: {e}"),
            Self::Shape(e) => write!(f, "constrained heating shape: {e}") }
    }
}
impl std::error::Error for ConstrainedHeatError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self { Self::Heating(e) => Some(e), Self::Shape(e) => Some(e) }
    }
}

/// Heat balance and actual number of constrained configurations in this advance.
#[derive(Clone, Debug, PartialEq)]
pub struct ConstrainedHeatReport {
    /// Existing thermal audit; generated heat is already included exactly once.
    /// This does not include the separately reported elastic potential energy.
    pub thermal: EnergyAudit,
    /// Actual successful shape solves, one per thermal substep, not output row.
    pub shape_solves: usize,
}

/// A restartable complete heat/shape state over one immutable reference inventory.
/// Mutable trajectories and source pointers cannot be replaced independently.
#[derive(Clone)]
pub struct ConstrainedHeatedBody<'a> {
    mesh: TetrahedralMass<'a>,
    regions: Vec<ConstrainedRegion<'a>>,
    fixed_dofs: Vec<usize>,
    nodal_forces_n: Vec<(usize, f64)>,
    prescribed_dofs_m: Vec<(usize, f64)>,
    shape_limits: ThermalSolveConfig,
    thermal: DepositionSimulation<'a>,
    configuration: ThermalBody<'a>,
}
impl<'a> ConstrainedHeatedBody<'a> {
    /// Bind the SAME source cards, mass inventory and caloric charts to heat and FEM.
    /// Uses the original single-body/mechanism storage initializer and its source,
    /// state-context, cell-volume, mass and complete-assignment checks. Initialize
    /// the original schedule at zero, then solve the initial physical shape once.
    /// All subsequent cell births refuse. Supports remain prescribed zero DOFs.
    ///
    /// No independent free-expansion geometry is built. Temperature support for
    /// Cp alone is insufficient: current elasticity and the whole alpha path must also
    /// support every physical endpoint. Thermal contact geometry stays reference
    /// geometry; mechanical deformation does not silently change heat transfer.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        mesh: TetrahedralMass<'a>, regions: Vec<ConstrainedRegion<'a>>,
        fixed_dofs: Vec<usize>, shape_limits: ThermalSolveConfig,
        events: Vec<TimedEvent>, budget: DepositionBudget, config: StepConfig, cx: &Cx<'_>,
    ) -> Result<Self, ConstrainedHeatError> {
        Self::new_with_nodal_forces(mesh,regions,fixed_dofs,Vec::new(),shape_limits,
            events,budget,config,cx)
    }

    /// Keep the same explicit mechanical dead loads in every constrained thermal
    /// substep. This composes heating with gravity, applied force or previously
    /// integrated pressure/traction through the original solid solve.
    ///
    /// Forces are sorted, unique full xyz DOFs in newtons, fixed in the reference
    /// frame. They remain immutable with this history and its source inventory.
    /// The exact accepted force vector also drives each state's vibration pencil.
    /// Empty forces preserve `new` behavior. Dynamic/follower loads and their
    /// mechanical work are not part of the existing heat-only energy audit.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_nodal_forces(
        mesh: TetrahedralMass<'a>, regions: Vec<ConstrainedRegion<'a>>,
        fixed_dofs: Vec<usize>, nodal_forces_n: Vec<(usize, f64)>,
        shape_limits: ThermalSolveConfig, events: Vec<TimedEvent>,
        budget: DepositionBudget, config: StepConfig, cx: &Cx<'_>,
    ) -> Result<Self, ConstrainedHeatError> {
        Self::new_with_prescribed_dofs(mesh, regions, fixed_dofs, nodal_forces_n,
            Vec::new(), shape_limits, events, budget, config, cx)
    }

    /// Heat a specimen held at explicit, fixed support offsets from its
    /// reference geometry. The same native Dirichlet lift supplies the initial
    /// prestress and every later equilibrium as temperature/material laws evolve.
    /// Offsets remain absolute: they are neither reset nor accumulated by steps,
    /// events, restart or cooling. All increments of the vibration pencil use
    /// the same fixed support components about the accepted equilibrium.
    ///
    /// Entries are sorted unique full xyz DOFs already present in `fixed_dofs`,
    /// in metres, bounded by `shape_limits.maximum_nodal_displacement_m`.
    /// Heating and generated-source histories publish only after shape accepts.
    /// Supports do not move during this timeline: no actuator work is inserted
    /// into the heat audit, and reciprocal thermoelastic heat is not inferred.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_prescribed_dofs(
        mesh: TetrahedralMass<'a>, mut regions: Vec<ConstrainedRegion<'a>>,
        fixed_dofs: Vec<usize>, nodal_forces_n: Vec<(usize, f64)>,
        prescribed_dofs_m: Vec<(usize, f64)>, shape_limits: ThermalSolveConfig,
        events: Vec<TimedEvent>, budget: DepositionBudget, config: StepConfig, cx: &Cx<'_>,
    ) -> Result<Self, ConstrainedHeatError> {
        poll(cx)?;
        if regions.len() != mesh.reference().regions().len()
            || fixed_dofs.len() > mesh.nodes_m().len().saturating_mul(3)
            || nodal_forces_n.len() > mesh.nodes_m().len().saturating_mul(3)
            || prescribed_dofs_m.len() > fixed_dofs.len() {
            return Err(HeatingError::Binding("complete thermal regions and bounded supports required").into());
        }
        regions.sort_by_key(|r| r.storage.key);
        let storage: Vec<_> = regions.iter().map(|r| r.storage.clone()).collect();
        let thermal = initialize(&[(mesh.reference(), &storage)], events, budget, config,
            &mut || poll(cx))?;
        let configuration = configure(&mesh, &regions, &fixed_dofs, &nodal_forces_n,
            &prescribed_dofs_m, shape_limits, &thermal, cx)?;
        poll(cx)?;
        Ok(Self { mesh, regions, fixed_dofs, nodal_forces_n, prescribed_dofs_m,
            shape_limits, thermal, configuration })
    }

    /// Original shared-node material mass map, never replaced while heating.
    pub const fn mesh(&self) -> &TetrahedralMass<'a> { &self.mesh }
    /// Actual thermal time, enthalpies, links and original source laws, read-only.
    pub const fn thermal(&self) -> &DepositionSimulation<'a> { &self.thermal }
    /// Latest compatible displacement, reactions, stress, geometry and full inertia.
    pub const fn configuration(&self) -> &ThermalBody<'a> { &self.configuration }
    /// Immutable original source/storage descriptions in canonical region order.
    pub fn regions(&self) -> &[ConstrainedRegion<'a>] { &self.regions }
    /// Immutable reference-frame mechanical dead loads shared by all endpoints.
    pub fn nodal_forces_n(&self) -> &[(usize, f64)] { &self.nodal_forces_n }
    /// Immutable absolute support offsets, shared by heat, shape and vibration.
    pub fn prescribed_dofs_m(&self) -> &[(usize, f64)] { &self.prescribed_dofs_m }

    /// Advance actual heat transport with a constrained configuration at EACH
    /// accepted substep, before endpoint events. A bad intermediate state cannot
    /// hide behind a valid final temperature. `max_shape_solves` tightens the
    /// thermal substep budget before work; it does not renew at schedule events.
    /// Zero admits a same-time query/events, but no positive-duration heat step.
    pub fn advance_to(
        &self, target_s: f64, config: StepConfig, max_shape_solves: usize, cx: &Cx<'_>,
    ) -> Result<(Self, ConstrainedHeatReport), ConstrainedHeatError> {
        let observer = ShapeObserver { body: self, cx, configuration: self.configuration.clone(), calls: 0 };
        let (thermal, observer, audit) = self.thermal.advance_observed_limited(
            &observer, target_s, config, max_shape_solves, || poll(cx))?;
        let mut next = self.clone();
        next.thermal = thermal;
        next.configuration = observer.configuration;
        poll(cx)?;
        Ok((next, ConstrainedHeatReport { thermal: audit, shape_solves: observer.calls }))
    }

    /// Independent generated heat and conductance owners on the SAME thermal
    /// substeps and original events. Their accepted histories publish only with
    /// successful constrained mechanics and the original final heat audit.
    /// Observer Clone must be rollback-safe. No source is re-added to boundaries.
    pub fn advance_with<O: DepositionObserver>(
        &self, observer: &O, target_s: f64, config: StepConfig,
        max_shape_solves: usize, cx: &Cx<'_>,
    ) -> Result<(Self, O, ConstrainedHeatReport), ObserverCouplingError<O::Error, ConstrainedHeatError>> {
        let shape = ShapeObserver { body: self, cx, configuration: self.configuration.clone(), calls: 0 };
        let pair = CoupledObservers::new(observer.clone(), shape);
        let (thermal, pair, audit) = self.thermal.advance_observed_limited(
            &pair, target_s, config, max_shape_solves, || poll(cx))?;
        let (observer, shape) = pair.into_parts();
        let mut next = self.clone();
        next.thermal = thermal;
        next.configuration = shape.configuration;
        poll(cx)?;
        Ok((next, observer, ConstrainedHeatReport { thermal: audit, shape_solves: shape.calls }))
    }
}

#[derive(Clone)]
struct ShapeObserver<'b, 'a, 'cx> {
    body: &'b ConstrainedHeatedBody<'a>,
    cx: &'b Cx<'cx>,
    configuration: ThermalBody<'a>,
    calls: usize,
}
impl DepositionObserver for ShapeObserver<'_, '_, '_> {
    type Error = ConstrainedHeatError;
    fn event(&mut self, event: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        if matches!(event, Event::Birth(_)) {
            return Err(HeatingError::Binding("fixed material mesh cannot accept new mass").into());
        }
        Ok(())
    }
    fn step(&mut self, _: f64, _: f64, _: &BTreeMap<CellId, f64>, thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        checkpoint()?;
        let configuration = configure(&self.body.mesh, &self.body.regions, &self.body.fixed_dofs,
            &self.body.nodal_forces_n, &self.body.prescribed_dofs_m,
            self.body.shape_limits, thermal, self.cx)?;
        checkpoint()?;
        self.calls = self.calls.checked_add(1).ok_or(DepositionError::Budget("shape solves"))?;
        self.configuration = configuration;
        Ok(())
    }
}

#[allow(clippy::too_many_arguments)]
fn configure<'a>(
    mesh: &TetrahedralMass<'a>, regions: &[ConstrainedRegion<'_>], fixed_dofs: &[usize],
    nodal_forces_n: &[(usize, f64)], prescribed_dofs_m: &[(usize, f64)], limits: ThermalSolveConfig,
    thermal: &DepositionSimulation<'_>, cx: &Cx<'_>,
) -> Result<ThermalBody<'a>, ConstrainedHeatError> {
    poll(cx)?;
    let mut cells = Vec::with_capacity(regions.len());
    for region in regions {
        poll(cx)?;
        let source = &region.storage;
        let cell = thermal.cells().get(&source.cell_id)
            .ok_or(HeatingError::Binding("mapped constrained thermal cell is absent"))?;
        let state = source.caloric.curve().state_at_specific_enthalpy(cell.specific_h_j_kg)
            .map_err(HeatingError::from)?;
        if !state.is_fully_solid() {
            return Err(HeatingError::Binding("non-solid heat state cannot enter small-strain solid mechanics").into());
        }
        cells.push(ThermalCell { key: source.key, card: source.card,
            temperature_k: state.temperature_k(), stiffness_selection: region.stiffness_selection.clone(),
            stiffness_model: region.stiffness_model,
            expansion_model: region.expansion_model,
            expansion_selection: source.expansion_selection.clone() });
    }
    Ok(mesh.solve_thermal_with_prescribed_dofs(&cells, fixed_dofs, nodal_forces_n,
        prescribed_dofs_m, limits, cx)?)
}

fn poll(cx: &Cx<'_>) -> Result<(), DepositionError> {
    cx.checkpoint().map_err(|_| DepositionError::Cancelled)
}
