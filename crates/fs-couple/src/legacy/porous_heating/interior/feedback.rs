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
//! require their own owners. Original material walls can refresh their stiffness,
//! viscosity and complete relaxation spectrum at their own native temperatures.
//! Authored mechanical coefficients and Robin data otherwise remain held.
//! Original quadratic boundary motion and applied mechanical traction remain
//! coherent fixed drives as the fluid and wall materials refresh. Their signed
//! source work is observed separately; only native positive physical losses
//! feed the existing thermal receivers.
//! Distributed volume injection retains its original spatial profile and is
//! projected again on each changed-temperature owner, without normalization.

use super::heat::{InteriorHeatAccount, InteriorHeatState};
use super::{
    InteriorAcousticHeating, InteriorCellHeating, InteriorHeatCell, InteriorHeatingError,
    InteriorHeatingStep, InteriorWallHeatCell, InteriorWallHeating, MAX_INTERIOR_HEATING_CELLS,
};
use crate::bernoulli_aperture::wall::material::LocalMaterialWall;
use crate::cavity_modal::interior::heterogeneous::mechanical::MaterialMechanicalHarmonicError;
use crate::cavity_modal::interior::heterogeneous::porous::AcousticVolumeRegionState;
use crate::porous_heating::PorousHeatingError;
use crate::porous_heating::feedback::{IsobaricPorousGas, PorousTemperatureFluid, fluid};
use fs_feec::acoustics::helmholtz::HelmholtzFe;
use fs_feec::acoustics::helmholtz::heterogeneous::time::mechanical::{
    MechanicalBoundary, traction::WallTractionProfile,
};
use fs_feec::acoustics::helmholtz::heterogeneous::{
    BoundaryAdmittance, EquivalentFluidControls, EquivalentFluidPressureField, HeterogeneousError,
    MAX_BOUNDARY_SOURCES, MAX_VOLUME_SOURCES, PreparedEquivalentFluid, VelocityTrace,
};
use fs_math::c64::C64;
use fs_material::MaterialError;
use fs_material::porous::PorousFluidLaw;
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionObserver, DepositionSimulation, Event, StepConfig,
};
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

mod wall;
pub use wall::TemperatureInteriorWall;

/// One complete original tetrahedron and its thermally lumped specimen law.
#[derive(Clone, Copy)]
pub struct TemperatureInteriorCell<'t, 'm> {
    pub destination: InteriorHeatCell<'t>,
    /// Original fixed porous skeleton. None selects the source's free-volume
    /// response, retaining molecular reservoirs when that source supplies them.
    pub porous: Option<&'m PorousFluidLaw>,
}

/// One prescribed quadratic inward-normal boundary shape.
///
/// Each original triangle carries dimensionless values in native order:
/// corner0, corner1, corner2, edge01, edge12, edge20. A complex peak velocity
/// [m/s] on the exp(+i omega t) clock multiplies the complete signed projection.
/// No area normalization or replacement by three corner values is performed.
#[derive(Clone, Debug)]
pub struct InteriorTemperatureBoundaryDrive {
    pub traces: Vec<([u32; 3], [f64; 6])>,
    /// Original boundary-source quadrature, from pressure order + 2 through 16.
    pub quadrature_order: usize,
    pub amplitude_m_s: C64,
}

/// Original immutable spatial volume-injection shape [m^-3].
///
/// Arguments are the original tetrahedron index and physical position [m].
/// Values must remain the same for the same arguments throughout the owner's
/// lifetime: changing captured state is not a cache-invalidation mechanism.
/// Signed and zero-integral profiles are valid. Callback failures propagate
/// through the original acoustic error without publishing a partial snapshot.
pub type InteriorTemperatureVolumeProfile =
    dyn Fn(usize, [f64; 3]) -> Result<f64, HeterogeneousError> + Send + Sync;

/// Distributed continuity injection on the original full pressure geometry.
///
/// The peak exp(+i omega t) flow [m^3/s] multiplies the spatial shape [m^-3].
/// Neither its integral nor its pressure-basis projection is normalized. No
/// effective density, porosity, mass/enthalpy injection or heat-release model
/// is inferred. Clone shares the same immutable profile; Debug never calls it.
#[derive(Clone)]
pub struct InteriorTemperatureVolumeDrive {
    pub density_per_m3: Arc<InteriorTemperatureVolumeProfile>,
    /// Independent native source quadrature, from pressure order + 2 to 16.
    pub quadrature_order: usize,
    pub amplitude_m3_s: C64,
}

impl core::fmt::Debug for InteriorTemperatureVolumeDrive {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("InteriorTemperatureVolumeDrive")
            .field("density_per_m3", &"<retained spatial profile>")
            .field("quadrature_order", &self.quadrature_order)
            .field("amplitude_m3_s", &self.amplitude_m3_s)
            .finish()
    }
}

/// Applied mechanical pressure on the original wall's full trace space.
///
/// The native profile uses original global hierarchical coordinates. A complex
/// peak pressure [Pa] acts outward from the fluid and drives wall momentum;
/// wall motion, its fluid reaction and every memory coordinate remain solved.
#[derive(Clone, Debug)]
pub struct InteriorTemperatureTractionDrive {
    pub profile: WallTractionProfile,
    pub amplitude_pa: C64,
}

/// Fixed frequency, boundary geometry and coherent physical source drives.
///
/// All original pressure and mechanical coordinates remain in the solve.
/// `try_new` holds the authored mechanical coefficients. The material-wall
/// constructor supplies and refreshes them through their original source laws;
/// it requires an empty input `mechanical_walls` instead of replacing a law.
#[derive(Clone, Debug)]
pub struct InteriorTemperatureAcoustics {
    pub frequency_hz: f64,
    pub inward_velocity: Vec<VelocityTrace>,
    /// Original quadratic boundary motion, simultaneous with the linear traces.
    pub boundary_drives: Vec<InteriorTemperatureBoundaryDrive>,
    /// Original distributed injection, coherent with every boundary/wall drive.
    pub volume_drives: Vec<InteriorTemperatureVolumeDrive>,
    /// Independent applied wall pressures, simultaneous with all fluid drives.
    pub traction_drives: Vec<InteriorTemperatureTractionDrive>,
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
    Wall(MaterialMechanicalHarmonicError),
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
            Self::Wall(error) => error.fmt(f),
        }
    }
}
impl std::error::Error for InteriorFeedbackError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Heating(error) => Some(error),
            Self::Fluid { source, .. } => Some(source),
            Self::Material { source, .. } => Some(source),
            Self::Wall(error) => Some(error),
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
impl From<MaterialMechanicalHarmonicError> for InteriorFeedbackError {
    fn from(value: MaterialMechanicalHarmonicError) -> Self {
        Self::Wall(value)
    }
}

/// Work performed during one advance, excluding initial preparation.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct InteriorFeedbackWork {
    pub fluid_evaluations: usize,
    pub porous_evaluations: usize,
    /// Changed wall-temperature evaluations of the original material bundles.
    pub wall_evaluations: usize,
    /// Complete global FE preparations/solves, including changed endpoints.
    pub field_solves: usize,
}

fn increment(value: &mut usize) -> Result<(), InteriorFeedbackError> {
    *value = value.checked_add(1).ok_or(DepositionError::Budget(
        "interior acoustic feedback evaluations",
    ))?;
    Ok(())
}

fn admit_count(elements: usize, walls: usize, maximum: usize) -> Result<(), InteriorFeedbackError> {
    if maximum == 0
        || maximum > MAX_INTERIOR_HEATING_CELLS
        || elements
            .checked_add(walls)
            .is_none_or(|total| total > maximum)
    {
        return Err(DepositionError::Budget("interior acoustic feedback cells").into());
    }
    Ok(())
}

fn admit_drives(
    problem: &InteriorTemperatureAcoustics,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<(), InteriorFeedbackError> {
    for (requested, maximum, what) in [
        (
            problem.boundary_drives.len(),
            MAX_BOUNDARY_SOURCES,
            "interior feedback quadratic boundary sources",
        ),
        (
            problem.volume_drives.len(),
            MAX_VOLUME_SOURCES,
            "interior feedback distributed volume sources",
        ),
        (
            problem.traction_drives.len(),
            MAX_VOLUME_SOURCES,
            "interior feedback mechanical traction sources",
        ),
    ] {
        checkpoint()?;
        if requested > maximum {
            return Err(HeterogeneousError::Limit {
                what,
                requested,
                maximum,
            }
            .into());
        }
    }
    for amplitude in problem
        .boundary_drives
        .iter()
        .map(|drive| drive.amplitude_m_s)
        .chain(problem.volume_drives.iter().map(|drive| drive.amplitude_m3_s))
        .chain(problem.traction_drives.iter().map(|drive| drive.amplitude_pa))
    {
        checkpoint()?;
        if !amplitude.re.is_finite() || !amplitude.im.is_finite() {
            return Err(InteriorHeatingError::Invalid(
                "interior feedback needs finite physical source amplitudes",
            )
            .into());
        }
    }
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
    wall_bindings: Arc<[TemperatureInteriorWall<'t, 'm>]>,
    wall_states: Vec<wall::State>,
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
        // Each refreshed factor owns its sources. Rebinding keeps the original
        // geometry, signed shapes and physical amplitudes; an old material/wall
        // preparation is never accepted by the new solve. Unchanged temperature
        // snapshots bypass this entire preparation and reuse their Arc.
        let mut boundary = Vec::with_capacity(problem.boundary_drives.len());
        for drive in &problem.boundary_drives {
            boundary.push(prepared.quadratic_velocity_source_with_quadrature(
                &drive.traces,
                drive.quadrature_order,
                &mut || cancelled(),
            )?);
        }
        let mut volume = Vec::with_capacity(problem.volume_drives.len());
        for drive in &problem.volume_drives {
            volume.push(prepared.volume_velocity_source(
                drive.quadrature_order,
                &mut |element, position| (drive.density_per_m3)(element, position),
                &mut || cancelled(),
            )?);
        }
        let mut traction = Vec::with_capacity(problem.traction_drives.len());
        for drive in &problem.traction_drives {
            traction.push(prepared.prepare_wall_traction(&drive.profile, &mut || cancelled())?);
        }
        let boundary_values = boundary
            .iter()
            .zip(&problem.boundary_drives)
            .map(|(source, drive)| (source, drive.amplitude_m_s))
            .collect::<Vec<_>>();
        let traction_values = traction
            .iter()
            .zip(&problem.traction_drives)
            .map(|(source, drive)| (source, drive.amplitude_pa))
            .collect::<Vec<_>>();
        let volume_values = volume
            .iter()
            .zip(&problem.volume_drives)
            .map(|(source, drive)| (source, drive.amplitude_m3_s))
            .collect::<Vec<_>>();
        prepared
            .solve_with_sources(
                &problem.inward_velocity,
                &boundary_values,
                &volume_values,
                &traction_values,
                &mut || cancelled(),
            )
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
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, InteriorFeedbackError> {
        Self::try_new_inner(
            geometry,
            thermal,
            cells,
            wall_destinations,
            &[],
            fluid,
            problem,
            maximum_cells,
            checkpoint,
        )
    }

    /// Bind each wall to its original temperature-dependent material source.
    ///
    /// `problem.mechanical_walls` must be empty: every supplied material binding
    /// contributes its complete original mechanical boundary. Cold preparation
    /// evaluates the wall at its receiver's native temperature, then performs
    /// one complete fluid/wall solve. Later wall-only changes also trigger that
    /// solve, without reevaluating unchanged fluid sources.
    ///
    /// The same fixed-mass, geometry, caloric and accepted-start heat bindings
    /// apply. Constitutive refresh is quasisteady; finite acoustic settling,
    /// transient Maxwell-history transfer and parametric stored-energy work
    /// require a transient owner and are not represented by harmonic snapshots.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new_with_material_walls(
        geometry: &'f HelmholtzFe<'c>,
        thermal: &DepositionSimulation<'_>,
        cells: &[TemperatureInteriorCell<'t, 'm>],
        wall_bindings: &[TemperatureInteriorWall<'t, 'm>],
        fluid: F,
        problem: InteriorTemperatureAcoustics,
        maximum_cells: usize,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, InteriorFeedbackError> {
        checkpoint()?;
        admit_count(
            geometry.tetrahedra().len(),
            wall_bindings.len(),
            maximum_cells,
        )?;
        if !problem.mechanical_walls.is_empty() {
            return Err(InteriorHeatingError::Invalid(
                "material wall feedback cannot replace authored mechanical boundaries",
            )
            .into());
        }
        let destinations: Vec<_> = wall_bindings.iter().map(|wall| wall.destination).collect();
        Self::try_new_inner(
            geometry,
            thermal,
            cells,
            &destinations,
            wall_bindings,
            fluid,
            problem,
            maximum_cells,
            checkpoint,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn try_new_inner(
        geometry: &'f HelmholtzFe<'c>,
        thermal: &DepositionSimulation<'_>,
        cells: &[TemperatureInteriorCell<'t, 'm>],
        wall_destinations: &[InteriorWallHeatCell<'t>],
        wall_bindings: &[TemperatureInteriorWall<'t, 'm>],
        fluid: F,
        mut problem: InteriorTemperatureAcoustics,
        maximum_cells: usize,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, InteriorFeedbackError> {
        checkpoint()?;
        admit_drives(&problem, &mut checkpoint)?;
        let count = geometry.tetrahedra().len();
        let wall_count = if wall_bindings.is_empty() {
            problem.mechanical_walls.len()
        } else {
            wall_bindings.len()
        };
        admit_count(count, wall_count, maximum_cells)?;
        if count == 0 || cells.len() != count || wall_destinations.len() != wall_count {
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
        let mut wall_states = Vec::with_capacity(wall_bindings.len());
        for (index, binding) in wall_bindings.iter().enumerate() {
            checkpoint()?;
            let temperature = thermal.temperature_k(binding.destination.cell)?;
            let (state, boundary) = acoustic(&mut checkpoint, |cancelled| {
                wall::evaluate(
                    binding,
                    index,
                    None,
                    temperature,
                    geometry,
                    problem.frequency_hz,
                    &mut || cancelled(),
                )
                .map_err(Into::into)
            })?;
            wall_states.push(state);
            problem.mechanical_walls.push(boundary);
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
            wall_bindings: wall_bindings.to_vec().into(),
            wall_states,
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
    /// Original wall cards, reference material states and physical bindings.
    /// Empty for a problem with explicitly held mechanical coefficients.
    pub fn wall_bindings(&self) -> &[TemperatureInteriorWall<'t, 'm>] {
        &self.wall_bindings
    }
    /// Complete wall material snapshots at the accepted native temperatures.
    pub fn wall_materials(&self) -> impl ExactSizeIterator<Item = &LocalMaterialWall> {
        self.wall_states.iter().map(|state| state.material.as_ref())
    }
    /// Original drive/controls and the current accepted mechanical coefficients.
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
        let mut changed_walls = Vec::new();
        for (index, binding) in self.wall_bindings.iter().enumerate() {
            checkpoint()?;
            let temperature = thermal.temperature_k(binding.destination.cell)?;
            if temperature.to_bits() == self.wall_states[index].temperature_k.to_bits() {
                continue;
            }
            let (state, boundary) = acoustic(checkpoint, |cancelled| {
                wall::evaluate(
                    binding,
                    index,
                    Some(&self.wall_states[index]),
                    temperature,
                    self.geometry,
                    self.problem.frequency_hz,
                    &mut || cancelled(),
                )
                .map_err(Into::into)
            })?;
            increment(&mut work.wall_evaluations)?;
            changed_walls.push((index, state, boundary));
        }
        if changed.is_empty() && changed_walls.is_empty() {
            return Ok(());
        }
        let mut responses = self.responses.clone();
        for &(element, _, response) in &changed {
            responses[element] = response;
        }
        let mut next_wall_states = None;
        let problem = if changed_walls.is_empty() {
            Arc::clone(&self.problem)
        } else {
            let mut problem = (*self.problem).clone();
            let mut states = self.wall_states.clone();
            for (index, state, boundary) in changed_walls {
                states[index] = state;
                problem.mechanical_walls[index] = boundary;
            }
            next_wall_states = Some(states);
            Arc::new(problem)
        };
        let snapshot = prepare(
            self.geometry,
            &problem,
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
        if let Some(states) = next_wall_states {
            self.wall_states = states;
        }
        self.problem = problem;
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
