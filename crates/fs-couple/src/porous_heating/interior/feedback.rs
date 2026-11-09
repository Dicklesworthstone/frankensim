//! Accepted-temperature material feedback in the original full pressure space.
//!
//! The same original isobaric gas/liquid and molecular laws used by the porous
//! stack owner refresh each tetrahedron at its native caloric temperature.
//! Changed coefficients participate in one complete FE pressure/wall solve.
//! Unchanged temperatures, including latent plateaus, reuse the immutable field.
//!
//! This is quasisteady cycle-mean heating: acoustic settling must be fast
//! compared with thermal evolution. Accepted-start watts drive each native
//! thermal step; an endpoint solve admits the next state without retroactively
//! changing its heat. Acoustic startup, parametric stored-energy work and moving
//! geometry are not inferred. Refine the thermal timestep for this feedback;
//! the native conduction CFL alone does not bound its temporal error.
//!
//! Pressure/composition are prescribed background conditions of the ORIGINAL
//! fluid source, not a finite sealed gas inventory. Thermal mass is unchanged,
//! never replaced by rho(T)*volume. The native caloric law must justify the
//! thermally lumped specimen storage and its shared fluid/skeleton temperature.
//! Mean mass transport, pressure evolution and separate constituent temperatures
//! require their own owners. Mechanical coefficients and Robin data are held
//! fixed explicitly, although wall motion and absorption respond to the field.

use super::heat::{InteriorHeatAccount, InteriorHeatState};
use super::{
    InteriorAcousticHeating, InteriorCellHeating, InteriorHeatCell, InteriorHeatingError,
    InteriorHeatingStep, InteriorWallHeatCell, InteriorWallHeating, MAX_INTERIOR_HEATING_CELLS,
};
use crate::cavity_modal::interior::heterogeneous::porous::AcousticVolumeRegionState;
use crate::porous_heating::PorousHeatingError;
use crate::porous_heating::feedback::{IsobaricPorousGas, PorousTemperatureFluid, fluid};
use fs_feec::acoustics::helmholtz::HelmholtzFe;
use fs_feec::acoustics::helmholtz::heterogeneous::time::mechanical::MechanicalBoundary;
use fs_feec::acoustics::helmholtz::heterogeneous::{
    BoundaryAdmittance, EquivalentFluidControls, EquivalentFluidPressureField, HeterogeneousError,
    PreparedEquivalentFluid, VelocityTrace,
};
use fs_material::MaterialError;
use fs_material::porous::PorousFluidLaw;
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionObserver, DepositionSimulation, Event, StepConfig,
};
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

/// One complete original tetrahedron and its thermally lumped specimen law.
#[derive(Clone, Copy)]
pub struct TemperatureInteriorCell<'t, 'm> {
    pub destination: InteriorHeatCell<'t>,
    /// Original fixed porous skeleton. None selects the source's free-volume
    /// response, retaining molecular reservoirs when that source supplies them.
    pub porous: Option<&'m PorousFluidLaw>,
}

/// Fixed frequency, boundaries and coherent physical velocity drive.
///
/// All original pressure and mechanical coordinates remain in the solve. These
/// mechanical coefficients are explicitly held; temperature-dependent material
/// wall laws must be refreshed by their own source owner, outside this closure.
#[derive(Clone, Debug)]
pub struct InteriorTemperatureAcoustics {
    pub frequency_hz: f64,
    pub inward_velocity: Vec<VelocityTrace>,
    pub admittances: Vec<BoundaryAdmittance>,
    pub mechanical_walls: Vec<MechanicalBoundary>,
    pub pressure_release_faces: Vec<[u32; 3]>,
    pub controls: EquivalentFluidControls,
}

/// Original source refusal with its spatial identity, or native solve/heat error.
#[derive(Debug)]
pub enum InteriorFeedbackError {
    Heating(InteriorHeatingError),
    Fluid {
        element: usize,
        source: PorousHeatingError,
    },
    Material {
        element: usize,
        source: MaterialError,
    },
}

impl core::fmt::Display for InteriorFeedbackError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Heating(error) => error.fmt(f),
            Self::Fluid { element, source } => {
                write!(f, "interior element {element} fluid: {source}")
            }
            Self::Material { element, source } => {
                write!(f, "interior element {element} material: {source}")
            }
        }
    }
}
impl std::error::Error for InteriorFeedbackError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Heating(error) => Some(error),
            Self::Fluid { source, .. } => Some(source),
            Self::Material { source, .. } => Some(source),
        }
    }
}
impl From<InteriorHeatingError> for InteriorFeedbackError {
    fn from(value: InteriorHeatingError) -> Self {
        Self::Heating(value)
    }
}
impl From<DepositionError> for InteriorFeedbackError {
    fn from(value: DepositionError) -> Self {
        Self::Heating(value.into())
    }
}
impl From<HeterogeneousError> for InteriorFeedbackError {
    fn from(value: HeterogeneousError) -> Self {
        Self::Heating(value.into())
    }
}

/// Work performed during one advance, excluding initial preparation.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct InteriorFeedbackWork {
    pub fluid_evaluations: usize,
    pub porous_evaluations: usize,
    /// Complete global FE preparations/solves, including changed endpoints.
    pub field_solves: usize,
}

fn increment(value: &mut usize) -> Result<(), InteriorFeedbackError> {
    *value = value.checked_add(1).ok_or(DepositionError::Budget(
        "interior acoustic feedback evaluations",
    ))?;
    Ok(())
}

struct InteriorTemperatureSnapshot<'f, 'c, 't> {
    prepared: PreparedEquivalentFluid<'f, 'c>,
    field: EquivalentFluidPressureField,
    heat: Arc<InteriorHeatState<'t>>,
}

/// Source-retaining temperature feedback with immutable native field snapshots.
///
/// The common original fluid source can be the existing gas or evidence-bearing
/// liquid or molecular adapter, or another implementation declaring the same
/// volume-response contract. Each cell retains its complete evaluated fluid,
/// reservoir and porous state. Original molecular laws reevaluate reservoirs;
/// RelaxationState snapshots are never silently temperature-rescaled. A custom
/// source without a declared volume response refuses instead of losing its
/// dispersive constitutive behavior through an equilibrium-only fallback.
///
/// Snapshot clones share the immutable factorization and physical heat state.
/// During a refresh the accepted factor remains live alongside the new native
/// preparation. Its retained storage is additional to the native per-preparation
/// work/memory allowance; native thermal substep limits bound the solve count.
#[derive(Clone)]
pub struct TemperatureInteriorHeating<
    'f,
    'c,
    't,
    'm,
    F: PorousTemperatureFluid = IsobaricPorousGas<'m>,
> {
    geometry: &'f HelmholtzFe<'c>,
    cells: Arc<[TemperatureInteriorCell<'t, 'm>]>,
    wall_destinations: Arc<[InteriorWallHeatCell<'t>]>,
    fluid: F,
    problem: Arc<InteriorTemperatureAcoustics>,
    maximum_cells: usize,
    states: Vec<F::State>,
    responses: Vec<AcousticVolumeRegionState>,
    snapshot: Arc<InteriorTemperatureSnapshot<'f, 'c, 't>>,
}

/// One jointly accepted native thermal candidate and refreshed acoustic owner.
pub struct InteriorTemperatureStep<
    'f,
    'c,
    't,
    'm,
    'h,
    F: PorousTemperatureFluid = IsobaricPorousGas<'m>,
> {
    pub heating: TemperatureInteriorHeating<'f, 'c, 't, 'm, F>,
    /// Original thermal audit and independent four-mechanism heat receipts.
    pub step: InteriorHeatingStep<'h>,
    pub work: InteriorFeedbackWork,
}

fn acoustic<T>(
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    solve: impl FnOnce(&mut dyn FnMut() -> bool) -> Result<T, InteriorFeedbackError>,
) -> Result<T, InteriorFeedbackError> {
    let mut refusal = None;
    let result = solve(&mut || match checkpoint() {
        Ok(()) => false,
        Err(error) => {
            refusal = Some(error);
            true
        }
    });
    match refusal {
        Some(error) => Err(error.into()),
        None => result,
    }
}

fn material<F: PorousTemperatureFluid>(
    fluid: &F,
    element: usize,
    state: &F::State,
    law: Option<&PorousFluidLaw>,
    omega: f64,
) -> Result<AcousticVolumeRegionState, InteriorFeedbackError> {
    fluid
        .volume_response(state, law, omega)
        .map_err(|source| InteriorFeedbackError::Material { element, source })
}

#[allow(clippy::too_many_arguments)]
fn prepare<'f, 'c, 't>(
    geometry: &'f HelmholtzFe<'c>,
    problem: &InteriorTemperatureAcoustics,
    cells: &[TemperatureInteriorCell<'t, '_>],
    wall_destinations: &[InteriorWallHeatCell<'t>],
    maximum_cells: usize,
    thermal: &DepositionSimulation<'_>,
    response_at: &mut impl FnMut(usize) -> Result<AcousticVolumeRegionState, InteriorFeedbackError>,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<InteriorTemperatureSnapshot<'f, 'c, 't>, InteriorFeedbackError> {
    let prepared = acoustic(checkpoint, |cancelled| {
        PreparedEquivalentFluid::new_with_coupled_mechanical_walls(
            geometry,
            problem.frequency_hz,
            &problem.admittances,
            &problem.mechanical_walls,
            &problem.pressure_release_faces,
            problem.controls,
            &mut |element, _| Ok::<_, InteriorFeedbackError>(response_at(element)?.coefficients()),
            &mut || cancelled(),
        )
    })?;
    let field = acoustic(checkpoint, |cancelled| {
        prepared
            .solve(&problem.inward_velocity, &mut || cancelled())
            .map_err(Into::into)
    })?;
    let destinations: Vec<_> = cells.iter().map(|cell| cell.destination).collect();
    let heat = {
        let frozen = InteriorAcousticHeating::try_new_with_mechanical_walls(
            &prepared,
            &field.pressure_pa,
            &field.walls,
            thermal,
            &destinations,
            wall_destinations,
            maximum_cells,
            &mut *checkpoint,
        )?;
        Arc::clone(&frozen.heat)
    };
    checkpoint()?;
    Ok(InteriorTemperatureSnapshot {
        prepared,
        field,
        heat,
    })
}

impl<'f, 'c, 't, 'm, F: PorousTemperatureFluid> TemperatureInteriorHeating<'f, 'c, 't, 'm, F> {
    /// Bind every original element and mechanical wall to its existing receiver.
    ///
    /// Spatial correspondence and the thermally lumped specimen approximation
    /// have the same explicit requirements as InteriorAcousticHeating. Each
    /// element is homogeneous at its original native temperature; no averaging
    /// of distinct cell temperatures or rescaling of porosity is introduced.
    /// Native FE work/storage admission precedes the first acoustic material
    /// evaluation; native caloric temperatures are read before that preparation.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        geometry: &'f HelmholtzFe<'c>,
        thermal: &DepositionSimulation<'_>,
        cells: &[TemperatureInteriorCell<'t, 'm>],
        wall_destinations: &[InteriorWallHeatCell<'t>],
        fluid: F,
        problem: InteriorTemperatureAcoustics,
        maximum_cells: usize,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, InteriorFeedbackError> {
        checkpoint()?;
        let count = geometry.tetrahedra().len();
        let total = count.checked_add(problem.mechanical_walls.len());
        if maximum_cells == 0
            || maximum_cells > MAX_INTERIOR_HEATING_CELLS
            || total.is_none_or(|total| total > maximum_cells)
        {
            return Err(DepositionError::Budget("interior acoustic feedback cells").into());
        }
        if count == 0
            || cells.len() != count
            || wall_destinations.len() != problem.mechanical_walls.len()
        {
            return Err(InteriorHeatingError::Invalid(
                "interior feedback needs every original element and wall receiver",
            )
            .into());
        }
        let mut ids = BTreeSet::new();
        for id in cells
            .iter()
            .map(|cell| cell.destination.cell)
            .chain(wall_destinations.iter().map(|cell| cell.cell))
        {
            checkpoint()?;
            if !ids.insert(id) {
                return Err(InteriorHeatingError::Invalid(
                    "interior feedback destinations must be distinct",
                )
                .into());
            }
        }
        let mut temperatures = Vec::with_capacity(count);
        for cell in cells {
            checkpoint()?;
            temperatures.push(thermal.temperature_k(cell.destination.cell)?);
        }
        let omega = core::f64::consts::TAU * problem.frequency_hz;
        let mut states = vec![None; count];
        let mut responses = vec![None; count];
        let snapshot = prepare(
            geometry,
            &problem,
            cells,
            wall_destinations,
            maximum_cells,
            thermal,
            &mut |element| {
                if let Some(response) = responses[element] {
                    return Ok(response);
                }
                let state = fluid::evaluate(&fluid, temperatures[element])
                    .map_err(|source| InteriorFeedbackError::Fluid { element, source })?;
                let response =
                    material::<F>(&fluid, element, &state, cells[element].porous, omega)?;
                states[element] = Some(state);
                responses[element] = Some(response);
                Ok(response)
            },
            &mut checkpoint,
        )?;
        let states =
            states
                .into_iter()
                .collect::<Option<Vec<_>>>()
                .ok_or(InteriorHeatingError::Invalid(
                    "original interior element was not evaluated",
                ))?;
        let responses = responses.into_iter().collect::<Option<Vec<_>>>().ok_or(
            InteriorHeatingError::Invalid("original interior material was not evaluated"),
        )?;
        checkpoint()?;
        Ok(Self {
            geometry,
            cells: cells.to_vec().into(),
            wall_destinations: wall_destinations.to_vec().into(),
            fluid,
            problem: Arc::new(problem),
            maximum_cells,
            states,
            responses,
            snapshot: Arc::new(snapshot),
        })
    }

    pub fn source(&self) -> &F {
        &self.fluid
    }
    pub fn fluid_states(&self) -> &[F::State] {
        &self.states
    }
    pub fn material_responses(&self) -> &[AcousticVolumeRegionState] {
        &self.responses
    }
    pub fn cells(&self) -> &[TemperatureInteriorCell<'t, 'm>] {
        &self.cells
    }
    pub fn wall_destinations(&self) -> &[InteriorWallHeatCell<'t>] {
        &self.wall_destinations
    }
    pub fn problem(&self) -> &InteriorTemperatureAcoustics {
        &self.problem
    }
    pub fn prepared(&self) -> &PreparedEquivalentFluid<'f, 'c> {
        &self.snapshot.prepared
    }
    pub fn field(&self) -> &EquivalentFluidPressureField {
        &self.snapshot.field
    }
    pub fn volume_heating(&self) -> impl ExactSizeIterator<Item = &InteriorCellHeating> {
        self.snapshot.heat.cells()
    }
    pub fn wall_heating(&self) -> impl ExactSizeIterator<Item = &InteriorWallHeating> {
        self.snapshot.heat.walls()
    }

    fn refresh(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        work: &mut InteriorFeedbackWork,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), InteriorFeedbackError> {
        self.snapshot.heat.validate(thermal, checkpoint)?;
        let omega = core::f64::consts::TAU * self.problem.frequency_hz;
        let mut changed = Vec::new();
        for (element, cell) in self.cells.iter().enumerate() {
            checkpoint()?;
            let temperature = thermal.temperature_k(cell.destination.cell)?;
            if temperature.to_bits() == self.fluid.temperature_k(&self.states[element]).to_bits() {
                continue;
            }
            let state = fluid::evaluate(&self.fluid, temperature)
                .map_err(|source| InteriorFeedbackError::Fluid { element, source })?;
            increment(&mut work.fluid_evaluations)?;
            let response = material::<F>(&self.fluid, element, &state, cell.porous, omega)?;
            if cell.porous.is_some() {
                increment(&mut work.porous_evaluations)?;
            }
            changed.push((element, state, response));
        }
        if changed.is_empty() {
            return Ok(());
        }
        let mut responses = self.responses.clone();
        for &(element, _, response) in &changed {
            responses[element] = response;
        }
        let snapshot = prepare(
            self.geometry,
            &self.problem,
            &self.cells,
            &self.wall_destinations,
            self.maximum_cells,
            thermal,
            &mut |element| Ok(responses[element]),
            checkpoint,
        )?;
        increment(&mut work.field_solves)?;
        checkpoint()?;
        // Commit the coherent material/field cache only after every original
        // source, solve, physical heat form and native binding admits it.
        for (element, state, _) in changed {
            self.states[element] = state;
        }
        self.responses = responses;
        self.snapshot = Arc::new(snapshot);
        Ok(())
    }

    /// Advance material feedback and the native thermal timeline atomically.
    /// The returned owner is already evaluated at the final native temperatures
    /// and can continue without an unnecessary initial solve. Refusal publishes
    /// neither its field cache nor its thermal cells or accumulated heat.
    pub fn advance_to<'h>(
        &self,
        initial: &DepositionSimulation<'h>,
        target_s: f64,
        config: StepConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<InteriorTemperatureStep<'f, 'c, 't, 'm, 'h, F>, InteriorFeedbackError> {
        checkpoint()?;
        let mut observer = FeedbackObserver {
            source: self.clone(),
            account: InteriorHeatAccount::new(&self.snapshot.heat),
            work: InteriorFeedbackWork::default(),
        };
        observer
            .source
            .refresh(initial, &mut observer.work, &mut checkpoint)?;
        let (thermal, mut observer, heat_report) =
            initial.advance_observed(&observer, target_s, config, &mut checkpoint)?;
        observer
            .source
            .refresh(&thermal, &mut observer.work, &mut checkpoint)?;
        let step = observer.account.finish(
            &observer.source.snapshot.heat,
            thermal,
            heat_report,
            &mut checkpoint,
        )?;
        checkpoint()?;
        Ok(InteriorTemperatureStep {
            heating: observer.source,
            step,
            work: observer.work,
        })
    }
}

#[derive(Clone)]
struct FeedbackObserver<'f, 'c, 't, 'm, F: PorousTemperatureFluid> {
    source: TemperatureInteriorHeating<'f, 'c, 't, 'm, F>,
    account: InteriorHeatAccount,
    work: InteriorFeedbackWork,
}

impl<F: PorousTemperatureFluid> DepositionObserver for FeedbackObserver<'_, '_, '_, '_, F> {
    type Error = InteriorFeedbackError;
    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.source.refresh(thermal, &mut self.work, checkpoint)?;
        Ok(self.source.snapshot.heat.sources())
    }
    fn event(&mut self, _: &Event, thermal: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        self.source
            .snapshot
            .heat
            .validate(thermal, &mut || Ok(()))
            .map_err(Into::into)
    }
    fn step(
        &mut self,
        _: f64,
        dt: f64,
        _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.source.snapshot.heat.validate(thermal, checkpoint)?;
        self.account
            .record(&self.source.snapshot.heat, dt, checkpoint)?;
        self.source.refresh(thermal, &mut self.work, checkpoint)
    }
}
