//! Original species diffusion, Fourier heat and complete volumetric stress.
//!
//! The shared three-dimensional grid owns the conservative face divergence
//! and boundary quadrature. The original gas source owns Maxwell--Stefan,
//! partial enthalpies, Fourier heat and Newtonian stress/work. Prepared field,
//! gradient and face buffers compose these owners without a second transport
//! law or an additional viscous-heating source.

use super::reservoir::Reservoirs;
use super::{
    Boundaries3d, BoundaryCondition, BoundaryImpulse, BoundarySide, CartesianGrid,
    MAX_MIXTURE_WORKSPACE_BYTES, MixtureEuler3d, Primitive, Totals, WIDTH_LIMIT, allocated,
    boundary_view, cartesian, cell, checked_totals, copy_values, stepping, sum_rows, tangents,
    update,
};
use crate::compressible::AxisBoundaries;
use crate::compressible::mixture::molecular::interval::{self, Action, Stage, parabolic};
use crate::compressible::stress;
pub use fs_gas::multicomponent::dissipation::{
    DissipationControl, DissipationError, SpatialTransport, StabilizedDissipationControl,
    TransportLaw, WallThermalCondition,
};
use fs_gas::multicomponent::{MixtureGas, SpatialPrimitive};
use fs_qty::Time;

mod adaptive;
mod coupled;
mod driving;
mod inertia;
mod reservoir;
mod stabilized;
pub(super) mod timed;
pub use adaptive::{AdaptiveDissipationControl, AdaptiveMolecularReport};
pub use coupled::{
    MolecularEuler3d, MolecularStepControl, MolecularStepReport3d, StabilizedMolecularStepControl,
    StabilizedMolecularStepReport3d,
};
pub use driving::WallAcceleration3d;
pub(super) use driving::WallDrive3d;
use inertia::WallInertia;
pub(super) use inertia::WallWork;
pub use inertia::{
    InertialWall3d, WallInertiaControl3d, WallInertiaExchange3d, WallInertiaReport3d,
};
pub use stabilized::StabilizedMolecularReport;

#[cfg(feature = "experimental-sgs")]
pub mod sgs;
#[cfg(feature = "experimental-sgs")]
use sgs::SubgridTransport;
#[cfg(feature = "experimental-sgs")]
mod reactive_lagrangian;

/// Molecular behavior of one physical reflecting wall, in three dimensions.
/// The original public path retains the shared planar/spatial wall condition.
pub type MolecularWall3d = crate::compressible::molecular::MolecularWall<3>;

/// Explicit molecular closure on source-owned Euler pressure outlets.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum PressureOutletPolicy {
    /// Pressure alone does not specify molecular exchange; refuse the outlet.
    #[default]
    Refuse,
    /// Set every species, viscous momentum and diffusive energy flux to zero
    /// on the outlet face. Interior molecular transport remains active.
    ///
    /// This assumes zero normal species diffusion, viscous traction and heat
    /// transfer at the selected truncation boundary. Advective species,
    /// pressure/momentum and enthalpy still use the original Euler outlet.
    /// No bath is inferred, no incoming composition is supplied, and no general
    /// nonreflecting or backflow-stable boundary claim follows.
    ZeroDiffusiveFlux,
}

/// Complete molecular interval; physical time remains with its composing owner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolecularReport {
    /// Complete represented physical duration; the Euler clock is unchanged.
    pub elapsed: Time,
    /// Accepted SSP RK2 substeps or opt-in RKL2 macrosteps.
    pub accepted_steps: usize,
    /// Rejected internal trials before successful substeps.
    pub rejected_steps: usize,
    /// Largest source/stencil parabolic estimate over accepted stages.
    pub maximum_parabolic_number: f64,
    /// Largest original Maxwell--Stefan row residual [1/m].
    pub maximum_gradient_residual_per_m: f64,
    /// Largest absolute net diffusive species flux [kg/(m² s)].
    pub maximum_mass_frame_residual_kg_m2_s: f64,
    /// Complete initial extensive state.
    pub before: Totals,
    /// Complete accepted extensive state.
    pub after: Totals,
    /// Inward exchange integrated over all accepted stages and substeps.
    pub boundary_change: Totals,
    /// `after-before-boundary_change` in original component order.
    pub conservation_residual: Totals,
    /// Original reference transform of the conservative energy residual [J].
    pub absolute_energy_residual_j: f64,
    /// Accepted finite tangential carrier storage and reciprocal wall exchange.
    pub wall_inertia: Option<WallInertiaReport3d>,
}

struct Workspace {
    current: Vec<f64>,
    stage: Vec<f64>,
    candidate: Vec<f64>,
    primitive: Vec<Primitive>,
    /// Row-major derivative axis, then velocity component: du_i/dx_j.
    gradient: Vec<[[f64; 3]; 3]>,
    viscosity: Vec<stress::Viscosity>,
    flux: [Vec<f64>; 3],
    boundary_step: Vec<f64>,
    boundary_trial: Vec<f64>,
    boundary_correction: Vec<f64>,
}

/// Prepared molecular transport on the original evolving-species 3-D field.
///
/// One stored face flux serves both neighbors, including exact periodic copies.
/// Reflecting faces default to stationary adiabatic slip walls and may be
/// configured as impermeable no-slip or finite Navier-slip walls with prescribed
/// tangential velocity and fixed-temperature or adiabatic conductive transfer.
/// Navier walls explicitly partition interfacial friction heat between the
/// fluid and prescribed wall. A transmissive end needs
/// the same complete, source-admitted reservoir used by the Euler face. Its
/// species diffusion, heat and stress use the physical half-cell distance and
/// the bath's enthalpy/velocity for boundary energy transfer. Uniform and
/// face-resolved reservoirs remain held throughout the requested interval.
/// Unbound open ends refuse. Pressure-only outlets also refuse unless the
/// caller explicitly selects [`PressureOutletPolicy::ZeroDiffusiveFlux`]. All
/// successful field stages reuse flat buffers; no independent reservoir
/// inventory is introduced.
/// Compatible cell-weighted molecular stress applies on periodic,
/// zero-diffusive-outlet and reflecting domains. No-slip, finite Navier and
/// natural free-slip walls share one strain and physical friction law; the
/// same boundary traction supplies actual wall work and Navier heat retains
/// its authored fluid/wall partition. Source Fourier heat remains separate.
/// Baths and authored equilibrium SGS wall laws retain the original full
/// discretization.
/// Its sixteen-byte/cell cache is included in the combined workspace budget.
/// Interior SGS stays a separate face model without a new energy guarantee.
pub struct MolecularTransport3d<L: TransportLaw> {
    source: SpatialTransport<L>,
    grid: CartesianGrid,
    boundaries: Boundaries3d,
    pressure_outlets: PressureOutletPolicy,
    walls: [MolecularWall3d; 6],
    inertia: Option<WallInertia>,
    bound_flow_workspace_bytes: usize,
    maximum_workspace_bytes: usize,
    boundary_impulses: Vec<f64>,
    last_step: Option<MolecularReport>,
    work: Workspace,
    stabilized: Option<stabilized::Workspace>,
    adaptive: Option<adaptive::Workspace>,
    #[cfg(feature = "experimental-sgs")]
    sgs: Option<SubgridTransport>,
}
impl<L: TransportLaw> MolecularTransport3d<L> {
    /// Combined numerical bytes for the supplied flow and prepared operator.
    /// Includes source/MS work, all stages/faces/gradients/boundary buffers and
    /// input state excess capacity. Immutable source/coefficient data, bounded
    /// stack scratch and arbitrary law storage/evaluation allocations are
    /// outside this numerical-buffer budget.
    pub fn required_workspace_bytes(
        flow: &MixtureEuler3d,
        source: &SpatialTransport<L>,
    ) -> Result<usize, DissipationError> {
        let checked = || -> Option<usize> {
            let width = flow.gas.species_count().checked_add(4)?;
            let faces = flow
                .grid
                .face_len(0)
                .checked_add(flow.grid.face_len(1))?
                .checked_add(flow.grid.face_len(2))?;
            let scalars = flow
                .grid
                .len()
                .checked_mul(3)?
                .checked_add(faces)?
                .checked_add(flow.grid.boundary_len().checked_mul(4)?)?
                .checked_mul(width)?;
            scalars
                .checked_mul(size_of::<f64>())?
                .checked_add(flow.grid.len().checked_mul(
                    size_of::<Primitive>()
                        + size_of::<[[f64; 3]; 3]>()
                        + size_of::<stress::Viscosity>(),
                )?)?
                .checked_add(source.workspace_bytes())?
                .checked_add(flow.workspace_bytes())
        };
        checked().ok_or(DissipationError::Invalid(
            "volumetric molecular workspace size overflow",
        ))
    }

    /// Bind identical original source/domain, physical grid and boundaries.
    /// The limit covers this operator together with the supplied flow. Pressure
    /// outlets retain the strict [`PressureOutletPolicy::Refuse`] default.
    pub fn new(
        flow: &MixtureEuler3d,
        source: SpatialTransport<L>,
        maximum_workspace_bytes: usize,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DissipationError> {
        Self::new_with_pressure_outlet_policy(
            flow,
            source,
            PressureOutletPolicy::Refuse,
            maximum_workspace_bytes,
            cancelled,
        )
    }

    /// Prepare the original molecular operator with an explicit outlet closure.
    ///
    /// The policy applies only to `BoundaryCondition::PressureOutlet` faces.
    /// Transmissive faces still require complete source-owned bath traces;
    /// periodic and reflecting faces retain their existing laws. The selected
    /// closure follows every SSP/RKL stage, retry and coupled composition.
    /// It adds no field storage and does not change Euler outlet admission,
    /// acoustic bounds, backflow refusal or conservative advective transfer.
    pub fn new_with_pressure_outlet_policy(
        flow: &MixtureEuler3d,
        source: SpatialTransport<L>,
        pressure_outlets: PressureOutletPolicy,
        maximum_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DissipationError> {
        poll(&mut cancelled)?;
        if source.gas() != flow.gas() {
            return Err(DissipationError::Invalid(
                "volumetric molecular source or caloric path differs",
            ));
        }
        reservoir::admit(flow.boundaries, &flow.reservoirs, pressure_outlets)?;
        admit_workspace(
            Self::required_workspace_bytes(flow, &source)?,
            maximum_workspace_bytes,
        )?;
        let grid = flow.grid;
        let width = flow.gas.species_count() + 4;
        let count = grid.len() * width;
        let boundary = grid.boundary_len() * width;
        let result = Self {
            source,
            grid,
            boundaries: flow.boundaries,
            pressure_outlets,
            walls: [MolecularWall3d::SlipAdiabatic; 6],
            inertia: None,
            bound_flow_workspace_bytes: flow.workspace_bytes(),
            maximum_workspace_bytes,
            boundary_impulses: allocated(boundary, 0.0)?,
            last_step: None,
            stabilized: None,
            adaptive: None,
            #[cfg(feature = "experimental-sgs")]
            sgs: None,
            work: Workspace {
                current: allocated(count, 0.0)?,
                stage: allocated(count, 0.0)?,
                candidate: allocated(count, 0.0)?,
                primitive: allocated(grid.len(), Primitive::default())?,
                gradient: allocated(grid.len(), [[0.0; 3]; 3])?,
                viscosity: allocated(grid.len(), stress::Viscosity::default())?,
                flux: [
                    allocated(grid.face_len(0) * width, 0.0)?,
                    allocated(grid.face_len(1) * width, 0.0)?,
                    allocated(grid.face_len(2) * width, 0.0)?,
                ],
                boundary_step: allocated(boundary, 0.0)?,
                boundary_trial: allocated(boundary, 0.0)?,
                boundary_correction: allocated(boundary, 0.0)?,
            },
        };
        admit_workspace(
            result.workspace_bytes() + flow.workspace_bytes(),
            maximum_workspace_bytes,
        )?;
        poll(&mut cancelled)?;
        Ok(result)
    }

    /// Declared molecular closure on this owner's Euler pressure outlets.
    pub const fn pressure_outlet_policy(&self) -> PressureOutletPolicy {
        self.pressure_outlets
    }

    /// Whether the complete domain admits the compatible molecular stress
    /// identity. Prescribed wall/bath work and any modeled SGS contribution
    /// retain their original qualifications.
    pub fn uses_compatible_molecular_stress(&self) -> bool {
        compatible_stress_domain(
            self.boundaries,
            &self.walls,
            #[cfg(feature = "experimental-sgs")]
            self.sgs.as_ref(),
        )
    }
    /// Configure one physical reflecting side. Periodic axes have no wall.
    /// Normal wall motion is refused until conservative moving geometry exists.
    pub fn with_wall(
        mut self,
        side: BoundarySide,
        wall: MolecularWall3d,
    ) -> Result<Self, DissipationError> {
        if self
            .inertia
            .is_some_and(|model| model.masses[side_index(side)] != 0.0)
        {
            return Err(DissipationError::Invalid(
                "finite-wall model and accepted velocity cannot be reconfigured",
            ));
        }
        let Some(condition) = physical_condition(self.boundaries, side) else {
            return Err(DissipationError::Invalid(
                "periodic molecular side has no physical wall",
            ));
        };
        if condition != BoundaryCondition::Reflecting {
            return Err(DissipationError::Invalid(
                "molecular wall model requires a reflecting Euler side",
            ));
        }
        wall.validate(
            side.axis(),
            self.gas().domain().minimum_temperature,
            self.gas().domain().maximum_temperature,
        )
        .map_err(DissipationError::Invalid)?;
        #[cfg(feature = "experimental-sgs")]
        if let Some(sgs) = &self.sgs {
            sgs.admit_wall(side, wall)?;
        }
        self.walls[side_index(side)] = wall;
        Ok(self)
    }

    /// Configured molecular model for one side, including periodic placeholders.
    pub fn wall(&self, side: BoundarySide) -> MolecularWall3d {
        self.walls[side_index(side)]
    }

    /// Exact original caloric source, species and admitted temperature interval.
    pub fn gas(&self) -> &MixtureGas {
        self.source.gas()
    }
    /// Original transport law and its source/material identities.
    pub fn law(&self) -> &L {
        self.source.law()
    }
    /// Owned numerical capacities plus the bounded finite-wall solver reserve.
    /// External flow, allocator metadata and arbitrary provider data are excluded.
    pub fn workspace_bytes(&self) -> usize {
        let bytes = self.source.workspace_bytes()
            + (self.boundary_impulses.capacity()
                + self.work.current.capacity()
                + self.work.stage.capacity()
                + self.work.candidate.capacity()
                + self.work.flux.iter().map(Vec::capacity).sum::<usize>()
                + self.work.boundary_step.capacity()
                + self.work.boundary_trial.capacity()
                + self.work.boundary_correction.capacity())
                * size_of::<f64>()
            + self.work.primitive.capacity() * size_of::<Primitive>()
            + self.work.gradient.capacity() * size_of::<[[f64; 3]; 3]>()
            + self.work.viscosity.capacity() * size_of::<stress::Viscosity>()
            + self
                .inertia
                .map_or(0, |_| WallInertiaControl3d::solver_workspace_bytes())
            + self
                .stabilized
                .as_ref()
                .map_or(0, stabilized::Workspace::workspace_bytes)
            + self
                .adaptive
                .as_ref()
                .map_or(0, adaptive::Workspace::workspace_bytes);
        #[cfg(feature = "experimental-sgs")]
        let bytes = bytes
            + self
                .sgs
                .as_ref()
                .map_or(0, SubgridTransport::workspace_bytes);
        bytes
    }
    /// Last completely published requested-duration interval.
    pub const fn last_step(&self) -> Option<&MolecularReport> {
        self.last_step.as_ref()
    }
    /// Outward molecular exchange on each face of the last complete interval.
    pub fn boundary_impulses(
        &self,
        side: BoundarySide,
    ) -> impl ExactSizeIterator<Item = BoundaryImpulse<'_>> {
        boundary_view(
            &self.boundary_impulses,
            &self.grid,
            self.gas().species_count(),
            side,
        )
    }

    /// Advance the complete requested interval atomically. Euler time and its
    /// advection report/loads continue to describe their original interval.
    pub fn advance(
        &mut self,
        flow: &mut MixtureEuler3d,
        duration: Time,
        controls: DissipationControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularReport, DissipationError> {
        let report = self.advance_unpublished(
            flow,
            duration,
            controls,
            self.maximum_workspace_bytes,
            &mut cancelled,
        )?;
        self.publish(report, None);
        Ok(report)
    }

    /// Advance with prescribed tangential wall acceleration and native wall work.
    ///
    /// Every SSP stage uses the wall velocity at its physical stage time.
    /// The accepted endpoint wall velocities publish with the complete field
    /// and molecular loads. Refusal preserves all of them. The Euler clock,
    /// geometry and wall thermal conditions retain their original semantics.
    pub fn advance_with_wall_acceleration(
        &mut self,
        flow: &mut MixtureEuler3d,
        duration: Time,
        controls: DissipationControl,
        accelerations: &[WallAcceleration3d],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularReport, DissipationError> {
        poll(&mut cancelled)?;
        let drive = self.prepare_wall_drive(accelerations, &mut cancelled)?;
        let endpoint = drive.as_ref().map(|d| d.at(duration.value())).transpose()?;
        let report = self.advance_driven_unpublished(
            flow,
            duration,
            controls,
            self.maximum_workspace_bytes,
            drive.as_ref().map(|d| (d, 0.0)),
            None,
            &WallWork::default(),
            &mut cancelled,
        )?;
        self.publish(report, endpoint);
        Ok(report)
    }

    pub(super) fn admit_time_resolved_profiles(
        &self,
        flow: &MixtureEuler3d,
        clock: timed::Clock,
        duration: f64,
    ) -> Result<(), DissipationError> {
        timed::admit(flow, clock, duration, self.inertia.is_some())
    }

    pub(super) fn prepare_wall_drive(
        &self,
        accelerations: &[WallAcceleration3d],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Option<WallDrive3d>, DissipationError> {
        if let Some(model) = self.inertia {
            for acceleration in accelerations {
                if model.masses[side_index(acceleration.side)] != 0.0
                    && acceleration.acceleration.iter().any(|value| *value != 0.0)
                {
                    return Err(DissipationError::Invalid(
                        "a finite wall cannot also have prescribed acceleration",
                    ));
                }
            }
        }
        WallDrive3d::new(self.walls, self.boundaries, accelerations, cancelled)
    }

    pub(super) fn trial_impulses(&self) -> &[f64] {
        &self.work.boundary_trial
    }
    pub(super) fn publish(
        &mut self,
        report: MolecularReport,
        endpoint_walls: Option<[MolecularWall3d; 6]>,
    ) {
        if let Some(walls) = endpoint_walls {
            self.walls = walls;
        }
        if let Some(inertia) = &report.wall_inertia {
            inertia.apply(&mut self.walls);
        }
        std::mem::swap(&mut self.boundary_impulses, &mut self.work.boundary_trial);
        self.last_step = Some(report);
    }
    pub(super) fn advance_unpublished(
        &mut self,
        flow: &mut MixtureEuler3d,
        duration: Time,
        controls: DissipationControl,
        maximum_workspace_bytes: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MolecularReport, DissipationError> {
        self.advance_driven_unpublished(
            flow,
            duration,
            controls,
            maximum_workspace_bytes,
            None,
            None,
            &WallWork::default(),
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn advance_driven_unpublished(
        &mut self,
        flow: &mut MixtureEuler3d,
        duration: Time,
        controls: DissipationControl,
        maximum_workspace_bytes: usize,
        drive: Option<(&WallDrive3d, f64)>,
        previous: Option<&WallInertiaReport3d>,
        wall_work: &WallWork,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MolecularReport, DissipationError> {
        self.advance_driven_with_profiles_unpublished(
            flow,
            duration,
            controls,
            maximum_workspace_bytes,
            drive,
            previous,
            wall_work,
            None,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn advance_driven_with_profiles_unpublished(
        &mut self,
        flow: &mut MixtureEuler3d,
        duration: Time,
        controls: DissipationControl,
        maximum_workspace_bytes: usize,
        drive: Option<(&WallDrive3d, f64)>,
        previous: Option<&WallInertiaReport3d>,
        wall_work: &WallWork,
        profiles: Option<(&mut super::reservoir_stages::Context<'_>, timed::Clock)>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MolecularReport, DissipationError> {
        poll(cancelled)?;
        controls.validate()?;
        if let Some((_, clock)) = profiles.as_ref() {
            timed::admit(flow, *clock, duration.value(), self.inertia.is_some())?;
        }
        let timed = profiles.is_some();
        let before = self.prepare_interval(flow, duration, maximum_workspace_bytes, cancelled)?;
        let mut reservoirs = match profiles {
            Some((context, clock)) => timed::Sampling::Timed {
                reservoirs: &mut flow.reservoirs,
                context,
                clock,
                duration: duration.value(),
            },
            None => timed::Sampling::Held(&flow.reservoirs),
        };
        let mut current_walls = self.walls;
        if let Some(previous) = previous {
            previous.apply(&mut current_walls);
        }
        let initial_walls = current_walls;
        let mut mechanical_work = [0.0; 6];
        let mut work_correction = [0.0; 6];
        let interval = interval::advance(
            duration.value(),
            controls,
            cancelled,
            |action, elapsed, cancelled| {
                let result = (|| -> Result<Stage, DissipationError> {
                    let (stage, endpoint, exchange) = match action {
                        Action::First => {
                            let walls = self.stage_walls(drive, elapsed, Some(&current_walls))?;
                            wall_work.field(self.inertia)?;
                            let trace =
                                reservoirs.at(self.source.gas(), &self.grid, elapsed, cancelled)?;
                            return spatial(
                                &mut self.source,
                                &self.grid,
                                self.boundaries,
                                self.pressure_outlets,
                                trace,
                                &walls,
                                &self.work.current,
                                &mut self.work.primitive,
                                &mut self.work.gradient,
                                &mut self.work.viscosity,
                                &mut self.work.flux,
                                #[cfg(feature = "experimental-sgs")]
                                self.sgs.as_mut(),
                                controls,
                                cancelled,
                            );
                        }
                        Action::Trial { dt, first, cached } => {
                            wall_work.trial(self.inertia)?;
                            if self.inertia.is_some() {
                                self.solve_wall_trial(
                                    reservoirs.active(),
                                    dt,
                                    controls,
                                    drive,
                                    elapsed,
                                    &current_walls,
                                    wall_work,
                                    cancelled,
                                )?
                            } else {
                                let stage = self.ssp_trial_sampled(
                                    &mut reservoirs,
                                    dt,
                                    controls,
                                    drive,
                                    elapsed,
                                    None,
                                    cached.then_some(first),
                                    wall_work,
                                    cancelled,
                                )?;
                                (stage, current_walls, [0.0; 6])
                            }
                        }
                    };
                    for (index, ((total, correction), change)) in self
                        .work
                        .boundary_trial
                        .iter_mut()
                        .zip(&mut self.work.boundary_correction)
                        .zip(&self.work.boundary_step)
                        .enumerate()
                    {
                        if index % 4096 == 0 {
                            poll(cancelled)?;
                        }
                        let adjusted = change - *correction;
                        let next = *total + adjusted;
                        *correction = (next - *total) - adjusted;
                        *total = next;
                    }
                    for index in 0..6 {
                        let adjusted = exchange[index] - work_correction[index];
                        let next = mechanical_work[index] + adjusted;
                        work_correction[index] = (next - mechanical_work[index]) - adjusted;
                        mechanical_work[index] = next;
                    }
                    current_walls = endpoint;
                    std::mem::swap(&mut self.work.current, &mut self.work.candidate);
                    Ok(stage)
                })();
                result.map_err(|error| reservoirs.failure(error))
            },
        )
        .map_err(|failure| failure.error)?;
        let inertia =
            self.inertia_report(&initial_walls, &current_walls, mechanical_work, wall_work)?;
        let report =
            self.complete_interval(flow, duration, before, interval, inertia, cancelled)?;
        if timed && duration.value() > 0.0 {
            flow.reservoirs
                .profiles
                .as_mut()
                .expect("prepared stage profiles")
                .commit();
        }
        Ok(report)
    }

    /// Original shared admission and private state/load initialization.
    fn prepare_interval(
        &mut self,
        flow: &MixtureEuler3d,
        duration: Time,
        maximum_workspace_bytes: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Totals, DissipationError> {
        if self.gas() != flow.gas() || self.grid != flow.grid || self.boundaries != flow.boundaries
        {
            return Err(DissipationError::Invalid(
                "volumetric molecular field ownership differs",
            ));
        }
        reservoir::admit(flow.boundaries, &flow.reservoirs, self.pressure_outlets)?;
        let bytes = self
            .workspace_bytes()
            .checked_add(flow.workspace_bytes())
            .ok_or(DissipationError::Invalid(
                "volumetric molecular workspace size overflow",
            ))?;
        admit_workspace(bytes, maximum_workspace_bytes)?;
        if !duration.value().is_finite() || duration.value() < 0.0 {
            return Err(DissipationError::Invalid(
                "nonnegative finite molecular duration required",
            ));
        }
        let species = self.gas().species_count();
        let before = checked_totals(&flow.states, &self.grid, species, cancelled)?;
        copy_values(&flow.states, &mut self.work.current, cancelled)?;
        for values in [
            &mut self.work.boundary_trial,
            &mut self.work.boundary_correction,
        ] {
            for tile in values.chunks_mut(4096) {
                poll(cancelled)?;
                tile.fill(0.0);
            }
        }
        Ok(before)
    }

    /// One original physical ledger and atomic field commit for both integrators.
    fn complete_interval(
        &mut self,
        flow: &mut MixtureEuler3d,
        duration: Time,
        before: Totals,
        interval: interval::Interval,
        wall_inertia: Option<WallInertiaReport3d>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MolecularReport, DissipationError> {
        let species = self.gas().species_count();
        let width = species + 4;
        let after = checked_totals(&self.work.current, &self.grid, species, cancelled)?;
        let outward = sum_rows(&self.work.boundary_trial, width, 1.0, cancelled)?;
        let boundary_change = Totals::from_row(&outward.map(|value| -value), species);
        let residual = std::array::from_fn::<_, WIDTH_LIMIT, _>(|d| {
            if d < width {
                after.component(d) - before.component(d) - boundary_change.component(d)
            } else {
                0.0
            }
        });
        if !residual.iter().all(|value| value.is_finite()) {
            return Err(DissipationError::Invalid(
                "volumetric molecular balance overflow",
            ));
        }
        let conservation_residual = Totals::from_row(&residual[..width], species);
        let report = MolecularReport {
            elapsed: duration,
            accepted_steps: interval.accepted_steps,
            rejected_steps: interval.rejected_steps,
            maximum_parabolic_number: interval.maximum_parabolic_number,
            maximum_gradient_residual_per_m: interval.maximum_gradient_residual_per_m,
            maximum_mass_frame_residual_kg_m2_s: interval.maximum_mass_frame_residual_kg_m2_s,
            before,
            after,
            boundary_change,
            conservation_residual,
            absolute_energy_residual_j: conservation_residual.absolute_energy(self.gas())?,
            wall_inertia,
        };
        // No copied prefix can escape a late cancellation. Excess input
        // capacity follows the actual allocation across this O(1) exchange.
        let bytes = flow.workspace_bytes - flow.states.capacity() * size_of::<f64>()
            + self.work.current.capacity() * size_of::<f64>();
        poll(cancelled)?;
        std::mem::swap(&mut flow.states, &mut self.work.current);
        flow.workspace_bytes = bytes;
        Ok(report)
    }
    /// Sample prescribed histories, then override only finite wall coordinates.
    fn stage_walls(
        &self,
        drive: Option<(&WallDrive3d, f64)>,
        elapsed: f64,
        finite: Option<&[MolecularWall3d; 6]>,
    ) -> Result<[MolecularWall3d; 6], DissipationError> {
        let mut walls = drive
            .map(|(d, offset)| d.at(offset + elapsed))
            .transpose()?
            .unwrap_or(self.walls);
        if let (Some(model), Some(finite)) = (self.inertia, finite) {
            for index in 0..6 {
                if model.masses[index] != 0.0 {
                    walls[index] = finite[index];
                }
            }
        }
        Ok(walls)
    }

    /// The single original SSP RK2 field trial. Probes only overwrite prepared
    /// stage/candidate buffers; the caller owns accepted-state/load publication.
    #[allow(clippy::too_many_arguments)]
    fn ssp_trial(
        &mut self,
        reservoirs: &Reservoirs,
        dt: f64,
        controls: DissipationControl,
        drive: Option<(&WallDrive3d, f64)>,
        elapsed: f64,
        midpoint: Option<&[MolecularWall3d; 6]>,
        cached: Option<Stage>,
        wall_work: &WallWork,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Stage, DissipationError> {
        self.ssp_trial_sampled(
            &mut timed::Sampling::Held(reservoirs),
            dt,
            controls,
            drive,
            elapsed,
            midpoint,
            cached,
            wall_work,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn ssp_trial_sampled(
        &mut self,
        reservoirs: &mut timed::Sampling<'_, '_>,
        dt: f64,
        controls: DissipationControl,
        drive: Option<(&WallDrive3d, f64)>,
        elapsed: f64,
        midpoint: Option<&[MolecularWall3d; 6]>,
        cached: Option<Stage>,
        wall_work: &WallWork,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Stage, DissipationError> {
        let width = self.gas().species_count() + 4;
        let walls = self.stage_walls(drive, elapsed, midpoint)?;
        let first = if let Some(first) = cached {
            first
        } else {
            wall_work.field(self.inertia)?;
            let trace = reservoirs.at(self.source.gas(), &self.grid, elapsed, cancelled)?;
            spatial(
                &mut self.source,
                &self.grid,
                self.boundaries,
                self.pressure_outlets,
                trace,
                &walls,
                &self.work.current,
                &mut self.work.primitive,
                &mut self.work.gradient,
                &mut self.work.viscosity,
                &mut self.work.flux,
                #[cfg(feature = "experimental-sgs")]
                self.sgs.as_mut(),
                controls,
                cancelled,
            )?
        };
        parabolic(first.rate * dt, controls.parabolic_cfl)?;
        update(
            self.source.gas(),
            &self.grid,
            &self.work.current,
            &self.work.current,
            &self.work.flux,
            &mut self.work.stage,
            dt,
            false,
            None,
            cancelled,
        )?;
        // Every probe starts a new SSP quadrature, never accumulating a
        // perturbed candidate into the enclosing physical interval.
        capture(
            &self.grid,
            width,
            &self.work.flux,
            0.5 * dt,
            false,
            &mut self.work.boundary_step,
            cancelled,
        )?;
        let walls = self.stage_walls(drive, elapsed + dt, midpoint)?;
        wall_work.field(self.inertia)?;
        let trace = reservoirs.at_stage(self.source.gas(), &self.grid, elapsed, dt, cancelled)?;
        let second = spatial(
            &mut self.source,
            &self.grid,
            self.boundaries,
            self.pressure_outlets,
            trace,
            &walls,
            &self.work.stage,
            &mut self.work.primitive,
            &mut self.work.gradient,
            &mut self.work.viscosity,
            &mut self.work.flux,
            #[cfg(feature = "experimental-sgs")]
            self.sgs.as_mut(),
            controls,
            cancelled,
        )?;
        parabolic(second.rate * dt, controls.parabolic_cfl)?;
        update(
            self.source.gas(),
            &self.grid,
            &self.work.current,
            &self.work.stage,
            &self.work.flux,
            &mut self.work.candidate,
            dt,
            true,
            None,
            cancelled,
        )?;
        capture(
            &self.grid,
            width,
            &self.work.flux,
            0.5 * dt,
            true,
            &mut self.work.boundary_step,
            cancelled,
        )?;
        Ok(first.merge(second))
    }
}

#[allow(clippy::too_many_arguments)]
fn spatial<'a, L: TransportLaw>(
    source: &mut SpatialTransport<L>,
    grid: &CartesianGrid,
    boundaries: Boundaries3d,
    pressure_outlets: PressureOutletPolicy,
    reservoirs: impl Into<reservoir::View<'a>>,
    walls: &[MolecularWall3d; 6],
    states: &[f64],
    primitive: &mut [Primitive],
    gradient: &mut [[[f64; 3]; 3]],
    viscosity: &mut [stress::Viscosity],
    fluxes: &mut [Vec<f64>; 3],
    #[cfg(feature = "experimental-sgs")] mut sgs: Option<&mut SubgridTransport>,
    controls: DissipationControl,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Stage, DissipationError> {
    let reservoirs = reservoirs.into();
    let species = source.gas().species_count();
    let width = species + 4;
    let compatible = compatible_stress_domain(
        boundaries,
        walls,
        #[cfg(feature = "experimental-sgs")]
        sgs.as_deref(),
    );
    let mut minimum_density = f64::INFINITY;
    let mut maximum_cell_stress_scale: f64 = 0.0;
    for (row, value) in states.chunks_exact(width).zip(primitive.iter_mut()) {
        *value = cell(row, species).primitive(source.gas(), &mut *cancelled)?;
        minimum_density = minimum_density.min(value.density);
    }
    if compatible {
        for (index, value) in viscosity.iter_mut().enumerate() {
            poll(cancelled)?;
            let p = primitive[index];
            let coefficients = source.coefficients_from_primitive(
                SpatialPrimitive {
                    partial_densities: &states[index * width..index * width + species],
                    normal_velocity: p.velocity[0],
                    transverse_velocity: [p.velocity[1], p.velocity[2]],
                    pressure: p.pressure,
                },
                &mut *cancelled,
            )?;
            *value = stress::Viscosity {
                shear: coefficients.dynamic_viscosity.value(),
                bulk: coefficients.bulk_viscosity.value(),
            };
            if !value.valid() || !value.scale::<3>().is_finite() {
                return Err(DissipationError::Invalid(
                    "nonfinite or negative cell viscosity",
                ));
            }
            maximum_cell_stress_scale = maximum_cell_stress_scale.max(value.scale::<3>());
        }
    }
    // A single nine-derivative cache supplies both tangential directions of
    // every face and the SGS strain. Compatible walls solve the coupled slip
    // traces locally; bath and equilibrium-wall domains keep their native ghosts.
    for (index, row) in gradient.iter_mut().enumerate() {
        poll(cancelled)?;
        if compatible {
            let value = stress::velocity_gradient(
                grid.spacing,
                primitive[index].velocity,
                |axis, side| {
                    let (neighbor, reflected) =
                        cartesian::neighbour(grid, boundaries, index, axis, side == 1);
                    if reflected {
                        let (velocity, slip, _) = walls[2 * axis + side].stress_boundary();
                        stress::GradientSample::Wall { velocity, slip }
                    } else {
                        stress::GradientSample::Cell(primitive[neighbor].velocity)
                    }
                },
            )
            .map_err(DissipationError::Invalid)?;
            for (axis, derivative) in row.iter_mut().enumerate() {
                for (component, entry) in derivative.iter_mut().enumerate() {
                    *entry = value[component][axis];
                }
            }
            continue;
        }
        for (axis, derivative) in row.iter_mut().enumerate() {
            let low = reservoir::velocity_neighbour(
                grid, boundaries, reservoirs, primitive, index, axis, false, walls,
            )?;
            let high = reservoir::velocity_neighbour(
                grid, boundaries, reservoirs, primitive, index, axis, true, walls,
            )?;
            for (component, value) in derivative.iter_mut().enumerate() {
                *value = (high[component] - low[component]) / (2.0 * grid.spacing[axis]);
            }
        }
    }
    #[cfg(feature = "experimental-sgs")]
    if let Some(sgs) = sgs.as_mut() {
        sgs.prepare(source.gas(), states, primitive, gradient, walls, cancelled)?;
    }
    let normal = |index: usize, axis: usize| {
        let [a, b] = tangents(axis);
        SpatialPrimitive {
            partial_densities: &states[index * width..index * width + species],
            normal_velocity: primitive[index].velocity[axis],
            transverse_velocity: [primitive[index].velocity[a], primitive[index].velocity[b]],
            pressure: primitive[index].pressure,
        }
    };
    let mut stats = Stage::default();
    let mut mechanical_faces = false;
    cartesian::faces::<DissipationError, _>(grid, boundaries, cancelled, |face, cancelled| {
        let axis = face.axis;
        let [a, b] = tangents(axis);
        let flux = &mut fluxes[axis];
        if let cartesian::FaceKind::PeriodicCopy { first } = face.kind {
            flux.copy_within(first * width..(first + 1) * width, face.index * width);
            return Ok(());
        }
        mechanical_faces |= matches!(
            face.kind,
            cartesian::FaceKind::Interior { .. }
                | cartesian::FaceKind::Boundary {
                    condition: BoundaryCondition::Reflecting,
                    ..
                }
        );
        let output = &mut flux[face.index * width..(face.index + 1) * width];
        let receipt = match face.kind {
            cartesian::FaceKind::Interior { left, right } => {
                let tangential = [a, b].map(|derivative| {
                    [axis, a, b].map(|component| {
                        0.5 * gradient[left][derivative][component]
                            + 0.5 * gradient[right][derivative][component]
                    })
                });
                let mut receipt = source.face_from_primitives(
                    normal(left, axis),
                    normal(right, axis),
                    grid.spacing[axis],
                    tangential,
                    controls.diffusion,
                    output,
                    &mut *cancelled,
                )?;
                let face_viscosity = stress::Viscosity {
                    shear: receipt.dynamic_viscosity_pa_s,
                    bulk: receipt.bulk_viscosity_pa_s,
                };
                if compatible
                    && (viscosity[left] != face_viscosity || viscosity[right] != face_viscosity)
                {
                    let traction = stress::traction(
                        axis,
                        grid.spacing[axis],
                        [primitive[left].velocity, primitive[right].velocity],
                        [viscosity[left], viscosity[right]],
                        |side, component, derivative| {
                            gradient[[left, right][side]][derivative][component]
                        },
                    );
                    let velocity = [axis, a, b].map(|component| {
                        0.5 * primitive[left].velocity[component]
                            + 0.5 * primitive[right].velocity[component]
                    });
                    output[species] = -traction[axis];
                    output[species + 1] = -traction[a];
                    output[species + 2] = -traction[b];
                    // Source enthalpy, Fourier and Dufour terms are retained
                    // BEFORE stress work; subtracting a huge old face work
                    // would irreversibly erase a finite scalar-energy flux.
                    output[species + 3] = receipt.nonmechanical_energy_flux_w_m2
                        - velocity[0] * traction[axis]
                        - velocity[1] * traction[a]
                        - velocity[2] * traction[b];
                    if !output.iter().all(|v| v.is_finite()) {
                        return Err(DissipationError::Invalid(
                            "compatible molecular face overflow",
                        ));
                    }
                }
                if compatible {
                    // Flux arithmetic may bypass replacement for matching
                    // samples; its mechanical admission still uses the actual
                    // cell coefficients and global minimum stage density.
                    receipt.diffusivity_m2_s = receipt.nonmechanical_diffusivity_m2_s.max(
                        viscosity[left]
                            .scale::<3>()
                            .max(viscosity[right].scale::<3>())
                            / minimum_density,
                    );
                    if !receipt.diffusivity_m2_s.is_finite() {
                        return Err(DissipationError::Invalid(
                            "compatible molecular face scale overflow",
                        ));
                    }
                }
                #[cfg(feature = "experimental-sgs")]
                let receipt = if let Some(sgs) = sgs.as_ref() {
                    sgs.add_face(
                        source.gas(),
                        states,
                        primitive,
                        gradient,
                        left,
                        right,
                        axis,
                        grid.spacing[axis],
                        receipt,
                        controls.diffusion,
                        output,
                        cancelled,
                    )?
                } else {
                    receipt
                };
                receipt
            }
            cartesian::FaceKind::Boundary {
                cell,
                sign,
                condition: BoundaryCondition::Reflecting,
            } => {
                let side = side_for(axis, sign > 0.0);
                let wall = walls[side_index(side)];
                let mut receipt = match wall {
                    MolecularWall3d::SlipAdiabatic => source.wall_from_primitive(
                        normal(cell, axis),
                        -sign * 2.0 * primitive[cell].velocity[axis] / grid.spacing[axis],
                        [gradient[cell][a][a], gradient[cell][b][b]],
                        output,
                        &mut *cancelled,
                    )?,
                    MolecularWall3d::NoSlip { velocity, thermal }
                    | MolecularWall3d::NavierSlip {
                        velocity,
                        thermal,
                        slip_length_m: 0.0,
                        ..
                    } => source.no_slip_wall_from_primitive(
                        normal(cell, axis),
                        0.5 * grid.spacing[axis],
                        sign,
                        [velocity[a], velocity[b]],
                        thermal,
                        -sign * 2.0 * primitive[cell].velocity[axis] / grid.spacing[axis],
                        // One constant prescribed wall velocity has zero
                        // tangential divergence at the physical face.
                        [0.0; 2],
                        output,
                        &mut *cancelled,
                    )?,
                    MolecularWall3d::NavierSlip { thermal, .. } => {
                        let distance = 0.5 * grid.spacing[axis];
                        let p = primitive[cell];
                        let ghost = wall
                            .velocity_ghost(axis, distance, p.velocity)
                            .ok_or(DissipationError::Invalid("missing Navier wall fluid trace"))?;
                        let trace = [a, b].map(|t| 0.5 * p.velocity[t] + 0.5 * ghost[t]);
                        let fraction = wall.tangential_trace_fraction(distance);
                        // The original source admits impermeable species,
                        // normal stress and Fourier heat at the actual local
                        // composition and cell/wall temperature. Evaluate it
                        // once, using the resolved fluid trace and its retained
                        // tangential derivatives, not the solid's velocity.
                        let receipt = source.no_slip_wall_from_primitive(
                            normal(cell, axis),
                            distance,
                            sign,
                            trace,
                            thermal,
                            -sign * 2.0 * p.velocity[axis] / grid.spacing[axis],
                            [
                                fraction * gradient[cell][a][a],
                                fraction * gradient[cell][b][b],
                            ],
                            output,
                            &mut *cancelled,
                        )?;
                        if !compatible {
                            // Baths and authored equilibrium wall laws keep
                            // their original complete finite-slip mechanics.
                            let (shear, mechanical_energy) = wall.shear_and_energy(
                                axis,
                                sign,
                                p.velocity,
                                p.temperature,
                                distance,
                                receipt.dynamic_viscosity_pa_s,
                                0.0,
                            );
                            output[species + 1] = -shear[a];
                            output[species + 2] = -shear[b];
                            output[species + 3] =
                                mechanical_energy + receipt.reduced_heat_flux_w_m2;
                            if !output.iter().all(|value| value.is_finite()) {
                                return Err(DissipationError::Invalid(
                                    "Navier molecular wall flux overflow",
                                ));
                            }
                        }
                        receipt
                    }
                };
                if compatible {
                    let (velocity, slip, fluid_heat_fraction) = wall.stress_boundary();
                    let mechanics = stress::wall_stress(
                        axis,
                        sign,
                        0.5 * grid.spacing[axis],
                        primitive[cell].velocity,
                        velocity,
                        slip,
                        viscosity[cell],
                        |component, derivative| gradient[cell][derivative][component],
                    )
                    .map_err(DissipationError::Invalid)?;
                    let traction = mechanics.traction;
                    output[species] = -traction[axis];
                    output[species + 1] = -traction[a];
                    output[species + 2] = -traction[b];
                    output[species + 3] = receipt.nonmechanical_energy_flux_w_m2
                        - velocity[0] * traction[0]
                        - velocity[1] * traction[1]
                        - velocity[2] * traction[2]
                        + sign * (1.0 - fluid_heat_fraction) * mechanics.slip_dissipation;
                    // Wall-temperature sampling still owns Fourier heat.
                    // The compatible mechanical operator uses cell values,
                    // not the unused midpoint wall viscosity's rate bound.
                    receipt.diffusivity_m2_s = receipt
                        .nonmechanical_diffusivity_m2_s
                        .max(viscosity[cell].scale::<3>() / minimum_density);
                    if !output.iter().all(|value| value.is_finite())
                        || !receipt.diffusivity_m2_s.is_finite()
                    {
                        return Err(DissipationError::Invalid(
                            "compatible molecular wall flux overflow",
                        ));
                    }
                }
                #[cfg(feature = "experimental-sgs")]
                let receipt = if let Some(sgs) = sgs.as_ref() {
                    sgs.add_wall(
                        source.gas(),
                        &states[cell * width..cell * width + species],
                        side,
                        wall,
                        primitive[cell],
                        0.5 * grid.spacing[axis],
                        receipt,
                        output,
                        cancelled,
                    )?
                } else {
                    receipt
                };
                receipt
            }
            cartesian::FaceKind::Boundary {
                cell,
                sign,
                condition: BoundaryCondition::Transmissive,
            } => {
                let side = side_for(axis, sign > 0.0);
                let face = reservoir::face_index(grid, cell, axis);
                let exterior = reservoirs
                    .normal_face(axis, sign, face, species)
                    .ok_or(DissipationError::Invalid("unbound molecular reservoir"))?;
                let interior = normal(cell, axis);
                let (left, right) = if sign < 0.0 {
                    (exterior, interior)
                } else {
                    (interior, exterior)
                };
                let tangential = reservoir::tangential_gradient(
                    grid,
                    boundaries,
                    reservoirs,
                    side,
                    cell,
                    &gradient[cell],
                )?;
                // The original source uses bath enthalpies and velocity for
                // the physical boundary energy effort. Modeled SGS flux is
                // not extrapolated across this externally maintained bath.
                source.reservoir_face_from_primitives(
                    left,
                    right,
                    0.5 * grid.spacing[axis],
                    tangential,
                    sign < 0.0,
                    controls.diffusion,
                    output,
                    &mut *cancelled,
                )?
            }
            cartesian::FaceKind::Boundary {
                condition: BoundaryCondition::PressureOutlet { .. },
                ..
            } => {
                if pressure_outlets != PressureOutletPolicy::ZeroDiffusiveFlux {
                    return Err(DissipationError::Invalid(
                        "pressure alone does not define a species/thermal molecular reservoir",
                    ));
                }
                // This is the declared physical diffusive boundary flux,
                // including stress work and species enthalpy. Copying the
                // interior velocity supplies the existing normal gradient
                // extrapolation; no bath or additional wall source is created.
                // Euler separately owns all advective outlet transfer.
                output.fill(0.0);
                return Ok(());
            }
            cartesian::FaceKind::PeriodicCopy { .. } => unreachable!(),
        };
        // Restore the field's global momentum order from the normal/tangent
        // frame without changing species or referenced total-energy fluxes.
        let momentum = [output[species], output[species + 1], output[species + 2]];
        output[species + axis] = momentum[0];
        output[species + a] = momentum[1];
        output[species + b] = momentum[2];
        stats.add(receipt);
        Ok(())
    })?;
    if compatible && mechanical_faces {
        // Also covers singleton/all-wall grids without an interior face. A
        // domain made solely of true zero-flux faces retains its zero rate.
        stats.rate = stats.rate.max(maximum_cell_stress_scale / minimum_density);
    }
    // Covers all three axes, full cross-stress stencils and doubled wall
    // gradients. Both SSP states also undergo original source admission.
    // This engineering interval bound is not a positivity/entropy theorem.
    stats.rate *= 4.0 * grid.spacing.iter().map(|h| h.recip().powi(2)).sum::<f64>();
    if !stats.rate.is_finite() || stats.rate < 0.0 {
        return Err(DissipationError::Invalid(
            "volumetric molecular parabolic scale overflow",
        ));
    }
    Ok(stats)
}

fn side_for(axis: usize, upper: bool) -> BoundarySide {
    match (axis, upper) {
        (0, false) => BoundarySide::XMin,
        (0, true) => BoundarySide::XMax,
        (1, false) => BoundarySide::YMin,
        (1, true) => BoundarySide::YMax,
        (2, false) => BoundarySide::ZMin,
        (2, true) => BoundarySide::ZMax,
        _ => unreachable!("Cartesian axis"),
    }
}

fn compatible_stress_domain(
    boundaries: Boundaries3d,
    walls: &[MolecularWall3d; 6],
    #[cfg(feature = "experimental-sgs")] subgrid: Option<&SubgridTransport>,
) -> bool {
    #[cfg(feature = "experimental-sgs")]
    if subgrid.is_some_and(SubgridTransport::has_equilibrium_wall) {
        return false;
    }
    walls.iter().all(|wall| wall.supports_compatible_stress())
        && (0..3).all(|axis| match boundaries.axis(axis) {
            AxisBoundaries::Periodic => true,
            AxisBoundaries::Ends { lower, upper } => [lower, upper].into_iter().all(|condition| {
                matches!(
                    condition,
                    BoundaryCondition::Reflecting | BoundaryCondition::PressureOutlet { .. }
                )
            }),
        })
}
fn side_index(side: BoundarySide) -> usize {
    match side {
        BoundarySide::XMin => 0,
        BoundarySide::XMax => 1,
        BoundarySide::YMin => 2,
        BoundarySide::YMax => 3,
        BoundarySide::ZMin => 4,
        BoundarySide::ZMax => 5,
    }
}
fn physical_condition(boundaries: Boundaries3d, side: BoundarySide) -> Option<BoundaryCondition> {
    let AxisBoundaries::Ends { lower, upper } = boundaries.axis(side.axis()) else {
        return None;
    };
    Some(if side.sign() < 0.0 { lower } else { upper })
}
#[allow(clippy::too_many_arguments)]
fn wall_velocity_sample(
    primitive: &[Primitive],
    neighbour: usize,
    cell: usize,
    axis: usize,
    upper: bool,
    reflected: bool,
    walls: &[MolecularWall3d; 6],
    wall_distance: f64,
) -> Result<[f64; 3], DissipationError> {
    if !reflected {
        return Ok(primitive[neighbour].velocity);
    }
    let side = side_for(axis, upper);
    let cell_velocity = primitive[cell].velocity;
    let velocity = walls[side_index(side)]
        .velocity_ghost(axis, wall_distance, cell_velocity)
        .unwrap_or_else(|| {
            std::array::from_fn(|component| {
                cell_velocity[component] * if component == axis { -1.0 } else { 1.0 }
            })
        });
    if !velocity.iter().all(|value| value.is_finite()) {
        return Err(DissipationError::Invalid(
            "molecular wall ghost velocity overflow",
        ));
    }
    Ok(velocity)
}

#[allow(clippy::too_many_arguments)]
fn capture(
    grid: &CartesianGrid,
    width: usize,
    fluxes: &[Vec<f64>; 3],
    duration: f64,
    add: bool,
    output: &mut [f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), DissipationError> {
    cartesian::capture::<DissipationError>(
        grid,
        width,
        duration,
        add,
        |axis, face, d| fluxes[axis][face * width + d],
        |face, d, value, add| {
            let out = &mut output[face * width + d];
            *out = if add { *out + value } else { value };
        },
        cancelled,
    )
}
fn admit_workspace(bytes: usize, maximum: usize) -> Result<(), DissipationError> {
    if bytes > maximum || bytes > MAX_MIXTURE_WORKSPACE_BYTES {
        Err(DissipationError::Invalid(
            "volumetric molecular total workspace budget",
        ))
    } else {
        Ok(())
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), DissipationError> {
    stepping::poll::<DissipationError>(cancelled)
}
