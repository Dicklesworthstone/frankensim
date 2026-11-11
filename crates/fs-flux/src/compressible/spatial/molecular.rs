//! Source-owned Fourier conduction and full three-dimensional Newtonian stress.
//!
//! The original transport law supplies actual temperature, pressure and frozen
//! composition coefficients. Compatible molecular stress samples cell values
//! on periodic, zero-flux and reflecting domains, with the full compatible
//! no-slip/Navier/free-slip wall trace. Fourier heat keeps its native face
//! sampling. Stress work and Fourier heat share one conservative total-energy
//! flux. Bulk viscous loss remains internal energy; interfacial slip friction
//! follows its authored fluid/wall heat partition, without an additional volume
//! heating source. These operators do not
//! resolve a physical shock thickness without an adequate grid and source data.

pub use super::mixture::molecular::MolecularWall3d;
use super::*;
use crate::compressible::molecular::parabolic;
pub use crate::compressible::molecular::{
    MolecularControl, MolecularError, MolecularStabilization, MolecularStepControl,
};
use fs_gas::multicomponent::dissipation::TransportLaw;
pub use fs_gas::multicomponent::dissipation::WallThermalCondition;
use fs_qty::{Pressure, Temperature};

mod coupled;
pub use coupled::{
    InertialWallLoad, InertialWallsControl, MolecularInertialWallsError,
    MolecularInertialWallsStep, Piston, PistonError,
};
pub use coupled::{
    MolecularEuler3d, MolecularStepReport3d, MovingMolecularStepReport3d,
    StabilizedMolecularStepReport3d,
};
pub(in crate::compressible::spatial) use coupled::{admit_planes, bounded_duration, displacement};
mod stabilized;
pub use stabilized::StabilizedMolecularReport3d;
mod reservoir;
mod timed;
mod vorticity;
pub use vorticity::{StressVorticityBudget3d, StressVorticityWorkspace3d};
mod stress_observation;
pub use stress_observation::StressField3d;

/// Explicit local or periodic dynamic subgrid transport for low-Mach studies.
#[cfg(feature = "experimental-sgs")]
pub mod sgs;

/// One complete molecular interval, without advancing the enclosing flow clock.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolecularReport3d {
    /// Complete represented duration [s].
    pub elapsed: Time,
    /// Accepted SSP RK2 substeps or opt-in RKL2 macrosteps.
    pub accepted_steps: usize,
    /// Rejected internal attempts.
    pub rejected_steps: usize,
    /// Largest accepted-stage parabolic estimate.
    pub maximum_parabolic_number: f64,
    /// Original extensive state.
    pub before: Totals,
    /// Accepted extensive state.
    pub after: Totals,
    /// Complete inward molecular boundary transfer.
    pub boundary_change: Totals,
    /// `after - before - boundary_change`.
    pub conservation_residual: Totals,
}

struct MolecularWorkspace {
    current: Vec<State>,
    stage: Vec<State>,
    candidate: Vec<State>,
    primitive: Vec<Primitive>,
    // gradient[velocity component][derivative axis], prepared once per cell.
    gradient: Vec<[[f64; 3]; 3]>,
    viscosity: Vec<crate::compressible::stress::Viscosity>,
    flux: [Vec<SpatialFaceFlux>; 3],
    binary: Vec<f64>,
    boundary_step: Vec<BoundaryImpulse>,
    boundary_trial: Vec<BoundaryImpulse>,
}

/// Prepared three-dimensional molecular operator bound to the exact gas source.
///
/// Periodic faces share one flux. Reflecting faces are stationary adiabatic
/// slip walls by default. Explicit no-slip or finite Navier-slip walls prescribe
/// tangential velocity and adiabatic or fixed-temperature conductive transfer.
/// Navier slip has an explicit interfacial-friction heat partition; its energy
/// flux includes wall mechanical work and friction heat sent to the wall.
/// Fourier flux uses the physical half-cell distance. SGS wall flux vanishes
/// by default; an explicit smooth equilibrium
/// stress model can supply adiabatic wall traction on a coarse SGS mesh, or
/// explicitly include integrated heat transfer at an isothermal wall.
/// Other walls require a resolved molecular wall layer. Transmissive ends
/// have zero molecular flux by default. [`Self::with_reservoir_exchange`]
/// enables heat and stress exchange through complete source-admitted uniform
/// or face-resolved reservoirs. All workspaces are preallocated.
/// Compatible cell-weighted molecular stress retains sixteen bytes of
/// viscosity scratch per cell, including prescribed wall motion, physical
/// Navier-slip dissipation, and the corresponding natural free-slip limit.
/// Enabled molecular bath exchange or an authored equilibrium SGS wall keeps
/// the original full mechanical discretization.
/// SGS retains its authored face closure and has no new dissipation guarantee.
pub struct MolecularTransport3d<L> {
    gas: FrozenGas,
    grid: CartesianGrid,
    boundaries: Boundaries3d,
    walls: [MolecularWall3d; 6],
    reservoir_exchange: bool,
    law: L,
    boundary_impulses: Vec<BoundaryImpulse>,
    work: MolecularWorkspace,
    stabilized: Option<stabilized::Workspace>,
    #[cfg(feature = "experimental-sgs")]
    sgs: Option<sgs::SubgridTransport>,
}

impl<L: TransportLaw> MolecularTransport3d<L> {
    /// Bind complete original source models, fixed composition and geometry.
    pub fn new(
        flow: &Euler3d,
        law: L,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MolecularError> {
        poll(&mut cancelled)?;
        if law.source().models() != flow.gas.source().components() {
            return Err(MolecularError::Invalid(
                "transport and Euler caloric sources differ",
            ));
        }
        let n = flow.states.len();
        let species = law.source().models().len();
        let binary = species.checked_mul(species).ok_or(MolecularError::Invalid(
            "molecular species workspace overflow",
        ))?;
        let faces = flow.grid.boundary_len();
        let result = Self {
            gas: flow.gas.clone(),
            grid: flow.grid,
            boundaries: flow.boundaries,
            walls: [MolecularWall3d::SlipAdiabatic; 6],
            reservoir_exchange: false,
            law,
            stabilized: None,
            #[cfg(feature = "experimental-sgs")]
            sgs: None,
            boundary_impulses: allocated(faces, Totals::default())?,
            work: MolecularWorkspace {
                current: allocated(n, State::default())?,
                stage: allocated(n, State::default())?,
                candidate: allocated(n, State::default())?,
                primitive: allocated(n, Primitive::default())?,
                gradient: allocated(n, [[0.0; 3]; 3])?,
                viscosity: allocated(n, crate::compressible::stress::Viscosity::default())?,
                flux: [
                    allocated(flow.grid.face_len(0), SpatialFaceFlux::default())?,
                    allocated(flow.grid.face_len(1), SpatialFaceFlux::default())?,
                    allocated(flow.grid.face_len(2), SpatialFaceFlux::default())?,
                ],
                binary: allocated(binary, 0.0)?,
                boundary_step: allocated(faces, Totals::default())?,
                boundary_trial: allocated(faces, Totals::default())?,
            },
        };
        poll(&mut cancelled)?;
        Ok(result)
    }

    /// Original immutable material/collision source and its reusable workspace.
    pub const fn law(&self) -> &L {
        &self.law
    }

    /// Owned numerical capacities: the state, primitive, gradient, viscosity,
    /// face-flux, binary-diffusivity and boundary-impulse workspaces, the RKL2
    /// buffers once prepared and, with `experimental-sgs`, the attached SGS
    /// closure. The external flow, allocator metadata and the transport law's own
    /// data are excluded. Attaching walls or a wall model reuses these capacities
    /// and leaves the value unchanged.
    pub fn workspace_bytes(&self) -> usize {
        let work = &self.work;
        let bytes = (work.current.capacity() + work.stage.capacity() + work.candidate.capacity())
            * size_of::<State>()
            + work.primitive.capacity() * size_of::<Primitive>()
            + work.gradient.capacity() * size_of::<[[f64; 3]; 3]>()
            + work.viscosity.capacity() * size_of::<crate::compressible::stress::Viscosity>()
            + work.flux.iter().map(Vec::capacity).sum::<usize>() * size_of::<SpatialFaceFlux>()
            + work.binary.capacity() * size_of::<f64>()
            + (self.boundary_impulses.capacity()
                + work.boundary_step.capacity()
                + work.boundary_trial.capacity())
                * size_of::<BoundaryImpulse>()
            + self.stabilized_workspace_bytes();
        #[cfg(feature = "experimental-sgs")]
        let bytes = bytes
            + self
                .sgs
                .as_ref()
                .map_or(0, sgs::SubgridTransport::workspace_bytes);
        bytes
    }

    /// Whether the complete domain admits the compatible molecular stress
    /// identity. This describes the molecular contribution; modeled SGS stress
    /// retains its own physical and numerical qualifications.
    pub fn uses_compatible_molecular_stress(&self) -> bool {
        #[cfg(feature = "experimental-sgs")]
        if self
            .sgs
            .as_ref()
            .is_some_and(sgs::SubgridTransport::has_equilibrium_wall)
        {
            return false;
        }
        !self.reservoir_exchange
            && self
                .walls
                .iter()
                .all(|wall| wall.supports_compatible_stress())
    }

    /// Enable molecular heat and full Newtonian stress at complete reservoirs.
    ///
    /// Each stage borrows the original Euler owner's held, source-admitted
    /// trace at the physical face. Normal gradients use half a cell; mechanical
    /// work uses the prescribed bath velocity. The original transport law is
    /// evaluated at the cell/bath midpoint temperature and pressure. Modeled
    /// bath flux requires a separate local SGS reservoir selection. No additional
    /// numerical storage is introduced. Unbound transmissive sides and
    /// pressure-only outlets retain zero molecular flux.
    /// This does not sample a time-dependent inlet during molecular stages.
    pub fn with_reservoir_exchange(mut self) -> Self {
        self.reservoir_exchange = true;
        self
    }

    /// Whether complete reservoirs supply the optional molecular bath law.
    pub const fn reservoir_exchange_enabled(&self) -> bool {
        self.reservoir_exchange
    }

    /// Configure an impermeable physical wall. Normal wall motion requires
    /// conservative moving geometry and is not represented by this operator.
    pub fn with_wall(
        mut self,
        side: BoundarySide,
        wall: MolecularWall3d,
    ) -> Result<Self, MolecularError> {
        let AxisBoundaries::Ends { lower, upper } = self.boundaries.axis(side.axis()) else {
            return Err(MolecularError::Invalid(
                "periodic side has no physical wall",
            ));
        };
        let condition = if side.sign() < 0.0 { lower } else { upper };
        if condition != BoundaryCondition::Reflecting {
            return Err(MolecularError::Invalid(
                "wall requires a reflecting Euler side",
            ));
        }
        #[cfg(feature = "experimental-sgs")]
        if !matches!(wall, MolecularWall3d::NoSlip { .. })
            && self
                .sgs
                .as_ref()
                .is_some_and(sgs::SubgridTransport::requires_no_slip_wall)
        {
            return Err(MolecularError::Invalid(
                "channel dynamic transport requires no-slip walls",
            ));
        }
        #[cfg(feature = "experimental-sgs")]
        if let Some(sgs) = &self.sgs {
            sgs.admit_wall(side, wall)?;
        }
        let domain = self.gas.domain();
        wall.validate(
            side.axis(),
            domain.minimum_temperature,
            domain.maximum_temperature,
        )
        .map_err(MolecularError::Invalid)?;
        self.walls[2 * side.axis() + usize::from(side.sign() > 0.0)] = wall;
        Ok(self)
    }

    /// Authored wall condition, including the slip default on unused sides.
    pub fn wall(&self, side: BoundarySide) -> MolecularWall3d {
        self.walls[2 * side.axis() + usize::from(side.sign() > 0.0)]
    }

    /// Outward molecular transfer during the last complete accepted interval.
    pub fn boundary_impulses(&self, side: BoundarySide) -> &[BoundaryImpulse] {
        &self.boundary_impulses[self.grid.boundary_range(side)]
    }

    /// Advance the whole molecular duration atomically, leaving time unchanged.
    ///
    /// The composing owner owns physical time. Existing Euler step metadata
    /// continues to describe its original advective interval; the returned
    /// molecular report and this operator's loads describe this update.
    pub fn advance(
        &mut self,
        flow: &mut Euler3d,
        duration: Time,
        controls: MolecularControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularReport3d, MolecularError> {
        self.advance_impl(
            flow,
            duration,
            controls,
            timed::Clock::interval(0.0, duration.value()),
            None,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_impl(
        &mut self,
        flow: &mut Euler3d,
        duration: Time,
        controls: MolecularControl,
        clock: timed::Clock,
        profiles: Option<&mut timed::TimedSource<'_>>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularReport3d, MolecularError> {
        self.advance_at_geometry(
            flow, duration, controls, clock, profiles, None, true, cancelled,
        )
        .map(|(report, _)| report)
    }

    // The coupled moving owner passes only its private trial field. Geometry
    // belongs to that particular half step, while this operator remains bound
    // to the accepted box until the complete D/2--ALE--D/2 transaction succeeds.
    #[allow(clippy::too_many_arguments)]
    fn advance_at_geometry(
        &mut self,
        flow: &mut Euler3d,
        duration: Time,
        controls: MolecularControl,
        clock: timed::Clock,
        mut profiles: Option<&mut timed::TimedSource<'_>>,
        motion: Option<MeshMotion>,
        publish: bool,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(MolecularReport3d, WallTransfer), MolecularError> {
        poll(&mut cancelled)?;
        controls.validate()?;
        if self.gas != flow.gas
            || self.boundaries != flow.boundaries
            || self.grid.shape != flow.grid.shape
            || (motion.is_none() && self.grid != flow.grid)
        {
            return Err(MolecularError::Invalid(
                "molecular operator and flow differ",
            ));
        }
        if !positive(duration.value()) {
            return Err(MolecularError::Invalid("molecular duration"));
        }
        if profiles.is_some() {
            clock.admit(duration.value())?;
        }
        let grid = flow.grid;
        let before = checked_totals(&flow.states, grid.cell_volume(), &mut cancelled)?;
        self.work.current.copy_from_slice(&flow.states);
        self.work.boundary_trial.fill(Totals::default());
        let mut wall_power = WallTransfer::default();
        let mut wall_transfer = WallTransfer::default();
        let mut elapsed = 0.0;
        let mut accepted = 0;
        let mut rejected = 0;
        let mut maximum_parabolic: f64 = 0.0;
        while elapsed < duration.value() {
            poll(&mut cancelled)?;
            if accepted + rejected >= controls.maximum_steps {
                return Err(MolecularError::WorkLimit {
                    elapsed_s: elapsed,
                    accepted,
                    rejected,
                });
            }
            let reservoirs = timed::stage(
                &mut flow.reservoirs,
                profiles.as_deref_mut(),
                &self.gas,
                &grid,
                clock.at(elapsed, duration.value()),
                self.reservoir_exchange,
                &mut cancelled,
            )?;
            let first_rate = spatial_with_motion(
                &self.gas,
                &grid,
                self.boundaries,
                &self.walls,
                motion,
                motion.map(|_| &mut wall_power),
                reservoirs,
                &self.work.current,
                &mut self.law,
                &mut self.work.primitive,
                &mut self.work.gradient,
                &mut self.work.viscosity,
                &mut self.work.flux,
                &mut self.work.binary,
                #[cfg(feature = "experimental-sgs")]
                self.sgs.as_mut(),
                &mut cancelled,
            )?;
            let remaining = duration.value() - elapsed;
            let mut dt = remaining.min(controls.maximum_step_s);
            if first_rate > 0.0 {
                dt = dt.min(controls.parabolic_cfl / first_rate);
            }
            for halving in 0..=controls.maximum_halvings {
                poll(&mut cancelled)?;
                if accepted + rejected >= controls.maximum_steps {
                    return Err(MolecularError::WorkLimit {
                        elapsed_s: elapsed,
                        accepted,
                        rejected,
                    });
                }
                if !positive(dt)
                    || elapsed + dt <= elapsed
                    || (dt < controls.minimum_step_s && dt != remaining)
                {
                    return Err(MolecularError::Invalid("unrepresentable molecular substep"));
                }
                let mut wall_step = WallTransfer::default();
                let trial = (|| {
                    // First-trial faces are already prepared for the CFL. Only
                    // a retry needs to reconstruct them after stage two reused
                    // the buffers; this avoids one whole source/face pass.
                    let r1 = if halving == 0 {
                        first_rate
                    } else {
                        let reservoirs = timed::stage(
                            &mut flow.reservoirs,
                            profiles.as_deref_mut(),
                            &self.gas,
                            &grid,
                            clock.at(elapsed, duration.value()),
                            self.reservoir_exchange,
                            &mut cancelled,
                        )?;
                        spatial_with_motion(
                            &self.gas,
                            &grid,
                            self.boundaries,
                            &self.walls,
                            motion,
                            motion.map(|_| &mut wall_power),
                            reservoirs,
                            &self.work.current,
                            &mut self.law,
                            &mut self.work.primitive,
                            &mut self.work.gradient,
                            &mut self.work.viscosity,
                            &mut self.work.flux,
                            &mut self.work.binary,
                            #[cfg(feature = "experimental-sgs")]
                            self.sgs.as_mut(),
                            &mut cancelled,
                        )?
                    };
                    parabolic(r1 * dt, controls.parabolic_cfl)?;
                    wall_step.add_scaled(wall_power, 0.5 * dt)?;
                    evolve::update(
                        &self.gas,
                        &grid,
                        &self.work.current,
                        &self.work.current,
                        &self.work.flux,
                        &mut self.work.stage,
                        dt,
                        false,
                        None,
                        None, // Gravity belongs to the native advective interval, not diffusion.
                        &mut cancelled,
                    )?;
                    evolve::capture_boundaries(
                        &grid,
                        &self.work.flux,
                        0.5 * dt,
                        &mut self.work.boundary_step,
                        false,
                        &mut cancelled,
                    )?;
                    let endpoint = if dt == remaining {
                        duration.value()
                    } else {
                        elapsed + dt
                    };
                    let reservoirs = timed::stage(
                        &mut flow.reservoirs,
                        profiles.as_deref_mut(),
                        &self.gas,
                        &grid,
                        clock.at(endpoint, duration.value()),
                        self.reservoir_exchange,
                        &mut cancelled,
                    )?;
                    let r2 = spatial_with_motion(
                        &self.gas,
                        &grid,
                        self.boundaries,
                        &self.walls,
                        motion,
                        motion.map(|_| &mut wall_power),
                        reservoirs,
                        &self.work.stage,
                        &mut self.law,
                        &mut self.work.primitive,
                        &mut self.work.gradient,
                        &mut self.work.viscosity,
                        &mut self.work.flux,
                        &mut self.work.binary,
                        #[cfg(feature = "experimental-sgs")]
                        self.sgs.as_mut(),
                        &mut cancelled,
                    )?;
                    parabolic(r2 * dt, controls.parabolic_cfl)?;
                    wall_step.add_scaled(wall_power, 0.5 * dt)?;
                    evolve::update(
                        &self.gas,
                        &grid,
                        &self.work.current,
                        &self.work.stage,
                        &self.work.flux,
                        &mut self.work.candidate,
                        dt,
                        true,
                        None,
                        None,
                        &mut cancelled,
                    )?;
                    evolve::capture_boundaries(
                        &grid,
                        &self.work.flux,
                        0.5 * dt,
                        &mut self.work.boundary_step,
                        true,
                        &mut cancelled,
                    )?;
                    Ok::<f64, MolecularError>(r1.max(r2) * dt)
                })();
                match trial {
                    Ok(number) => {
                        wall_transfer.add_scaled(wall_step, 1.0)?;
                        maximum_parabolic = maximum_parabolic.max(number);
                        for (i, (total, step)) in self
                            .work
                            .boundary_trial
                            .iter_mut()
                            .zip(&self.work.boundary_step)
                            .enumerate()
                        {
                            if i % 64 == 0 {
                                poll(&mut cancelled)?;
                            }
                            *total = add_totals(*total, *step);
                        }
                        std::mem::swap(&mut self.work.current, &mut self.work.candidate);
                        elapsed = if dt == remaining {
                            duration.value()
                        } else {
                            elapsed + dt
                        };
                        accepted += 1;
                        break;
                    }
                    Err(MolecularError::Cancelled) => return Err(MolecularError::Cancelled),
                    Err(error) if profiles.as_deref().is_some_and(timed::TimedSource::failed) => {
                        return Err(error);
                    }
                    Err(error) => {
                        rejected += 1;
                        if halving == controls.maximum_halvings
                            || dt * 0.5 < controls.minimum_step_s
                        {
                            return Err(MolecularError::StepFailure {
                                elapsed_s: elapsed,
                                step_s: dt,
                                last: Box::new(error),
                            });
                        }
                        dt *= 0.5;
                    }
                }
            }
        }
        let after = checked_totals(&self.work.current, grid.cell_volume(), &mut cancelled)?;
        let outward = sum_arrays(
            self.work.boundary_trial.iter().map(|q| q.array()),
            &mut cancelled,
        )?;
        let boundary_change = Totals::from_array(outward.map(|v| -v));
        let conservation_residual = balance(before, after, boundary_change)?;
        poll(&mut cancelled)?;
        flow.states.copy_from_slice(&self.work.current);
        if publish {
            self.boundary_impulses
                .copy_from_slice(&self.work.boundary_trial);
        }
        if profiles.is_some() {
            flow.reservoirs.profiles.publish();
        }
        Ok((
            MolecularReport3d {
                elapsed: duration,
                accepted_steps: accepted,
                rejected_steps: rejected,
                maximum_parabolic_number: maximum_parabolic,
                before,
                after,
                boundary_change,
                conservation_residual,
            },
            wall_transfer,
        ))
    }
}

fn add_totals(a: Totals, b: Totals) -> Totals {
    let a = a.array();
    let b = b.array();
    Totals::from_array(std::array::from_fn(|k| a[k] + b[k]))
}

fn balance(before: Totals, after: Totals, boundary: Totals) -> Result<Totals, MolecularError> {
    let a = after.array();
    let b = before.array();
    let c = boundary.array();
    let result = Totals::from_array(std::array::from_fn(|k| a[k] - b[k] - c[k]));
    if !a
        .into_iter()
        .chain(b)
        .chain(c)
        .chain(result.array())
        .all(f64::is_finite)
    {
        return Err(MolecularError::Invalid(
            "unrepresentable molecular extensive balance",
        ));
    }
    Ok(result)
}

/// Separate inward mechanical and thermal transfers, in boundary-side order.
/// Used as power during a face pass and as energy after SSP quadrature.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
struct WallTransfer {
    work: [f64; 6],
    heat: [f64; 6],
}

impl WallTransfer {
    fn add_scaled(&mut self, rhs: Self, scale: f64) -> Result<(), MolecularError> {
        for side in 0..6 {
            self.work[side] += scale * rhs.work[side];
            self.heat[side] += scale * rhs.heat[side];
        }
        if !self.work.into_iter().chain(self.heat).all(f64::is_finite) {
            return Err(MolecularError::Invalid("nonfinite moving wall transfer"));
        }
        Ok(())
    }
}

#[allow(clippy::too_many_arguments)]
fn spatial<L: TransportLaw>(
    gas: &FrozenGas,
    grid: &CartesianGrid,
    boundaries: Boundaries3d,
    walls: &[MolecularWall3d; 6],
    reservoirs: Option<super::reservoir::View<'_>>,
    states: &[State],
    law: &mut L,
    primitive: &mut [Primitive],
    gradient: &mut [[[f64; 3]; 3]],
    viscosity: &mut [crate::compressible::stress::Viscosity],
    fluxes: &mut [Vec<SpatialFaceFlux>; 3],
    binary: &mut [f64],
    #[cfg(feature = "experimental-sgs")] subgrid: Option<&mut sgs::SubgridTransport>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, MolecularError> {
    spatial_with_motion(
        gas,
        grid,
        boundaries,
        walls,
        None,
        None,
        reservoirs,
        states,
        law,
        primitive,
        gradient,
        viscosity,
        fluxes,
        binary,
        #[cfg(feature = "experimental-sgs")]
        subgrid,
        cancelled,
    )
}

#[allow(clippy::too_many_arguments)]
fn spatial_with_motion<L: TransportLaw>(
    gas: &FrozenGas,
    grid: &CartesianGrid,
    boundaries: Boundaries3d,
    walls: &[MolecularWall3d; 6],
    motion: Option<MeshMotion>,
    mut wall_power: Option<&mut WallTransfer>,
    reservoirs: Option<super::reservoir::View<'_>>,
    states: &[State],
    law: &mut L,
    primitive: &mut [Primitive],
    gradient: &mut [[[f64; 3]; 3]],
    viscosity: &mut [crate::compressible::stress::Viscosity],
    fluxes: &mut [Vec<SpatialFaceFlux>; 3],
    binary: &mut [f64],
    #[cfg(feature = "experimental-sgs")] mut subgrid: Option<&mut sgs::SubgridTransport>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, MolecularError> {
    if let Some(power) = wall_power.as_deref_mut() {
        *power = WallTransfer::default();
    }
    if reservoirs.is_some_and(|view| view.has_stagnation()) {
        return Err(MolecularError::Invalid(
            "compact stagnation reservoirs do not specify a molecular face bath",
        ));
    }
    let compatible =
        reservoirs.is_none() && walls.iter().all(|wall| wall.supports_compatible_stress());
    #[cfg(feature = "experimental-sgs")]
    let compatible = compatible
        && subgrid
            .as_ref()
            .is_none_or(|sgs| !sgs.has_equilibrium_wall());
    if motion.is_some() && !compatible {
        return Err(MolecularError::Invalid(
            "moving molecular walls require the compatible molecular stress",
        ));
    }
    let mut minimum_density = f64::INFINITY;
    let mut minimum_cv = f64::INFINITY;
    for (q, p) in states.iter().zip(primitive.iter_mut()) {
        *p = q.primitive(gas, &mut *cancelled)?;
        let props = gas
            .source()
            .evaluate_heat_capacities(Temperature::new(p.temperature))
            .map_err(Error::from)?;
        let volumetric_cv = p.density * props.mass_specific_cv().value();
        if !positive(volumetric_cv) {
            return Err(MolecularError::Invalid(
                "molecular volumetric heat capacity",
            ));
        }
        minimum_density = minimum_density.min(p.density);
        minimum_cv = minimum_cv.min(volumetric_cv);
    }
    if compatible {
        for (p, value) in primitive.iter().zip(viscosity.iter_mut()) {
            poll(cancelled)?;
            let coefficients = law.evaluate(
                Temperature::new(p.temperature),
                Pressure::new(p.pressure),
                gas.source().mole_composition().fractions(),
                binary,
                cancelled,
            )?;
            *value = crate::compressible::stress::Viscosity {
                shear: coefficients.dynamic_viscosity.value(),
                bulk: coefficients.bulk_viscosity.value(),
            };
            if !value.valid() || !value.scale::<3>().is_finite() {
                return Err(MolecularError::Invalid(
                    "nonfinite or negative cell viscosity",
                ));
            }
        }
    }
    // Cell central gradients are reused by the six adjacent face traces. This
    // retains the same centered stencil without redoing every neighbor lookup.
    for (cell, derivatives) in gradient.iter_mut().enumerate() {
        if cell % 64 == 0 {
            poll(cancelled)?;
        }
        let c = grid.coordinates(cell);
        if compatible {
            *derivatives = crate::compressible::stress::velocity_gradient(
                grid.spacing,
                primitive[cell].velocity,
                |axis, side| {
                    let outside = if side == 0 {
                        c[axis] == 0
                    } else {
                        c[axis] + 1 == grid.shape[axis]
                    };
                    if outside
                        && let AxisBoundaries::Ends { lower, upper } = boundaries.axis(axis)
                        && [lower, upper][side] == BoundaryCondition::Reflecting
                    {
                        let (mut velocity, slip, _) = walls[2 * axis + side].stress_boundary();
                        if let Some(motion) = motion {
                            velocity[axis] = if side == 0 {
                                motion.lower_velocity[axis]
                            } else {
                                motion.upper_velocity[axis]
                            };
                        }
                        return crate::compressible::stress::GradientSample::Wall {
                            velocity,
                            slip,
                        };
                    }
                    let sample = evolve::neighbour(primitive, grid, boundaries, c, axis, side == 1);
                    crate::compressible::stress::GradientSample::Cell([
                        sample[1], sample[2], sample[3],
                    ])
                },
            )
            .map_err(MolecularError::Invalid)?;
            continue;
        }
        for axis in 0..3 {
            let mut lower = evolve::neighbour(primitive, grid, boundaries, c, axis, false);
            let mut upper = evolve::neighbour(primitive, grid, boundaries, c, axis, true);
            // Mirror about the actual fluid trace: the prescribed velocity
            // for no slip, the resolved Navier trace for finite slip.
            for (outside, sample, wall, sign) in [
                (c[axis] == 0, &mut lower, walls[2 * axis], -1.0),
                (
                    c[axis] + 1 == grid.shape[axis],
                    &mut upper,
                    walls[2 * axis + 1],
                    1.0,
                ),
            ] {
                if outside
                    && let Some(ghost) = wall.velocity_ghost(
                        axis,
                        0.5 * grid.spacing[axis],
                        primitive[cell].velocity,
                    )
                {
                    sample[1..4].copy_from_slice(&ghost);
                } else if outside
                    && let Some(exterior) =
                        reservoirs.and_then(|r| r.get_at(grid, axis, sign, cell))
                {
                    // Mirror about the physical bath trace, half a cell away.
                    let ghost: [f64; 3] = std::array::from_fn(|d| {
                        2.0 * exterior.primitive().velocity[d] - primitive[cell].velocity[d]
                    });
                    if !ghost.iter().all(|v| v.is_finite()) {
                        return Err(MolecularError::Invalid(
                            "molecular reservoir ghost velocity overflow",
                        ));
                    }
                    sample[1..4].copy_from_slice(&ghost);
                }
            }
            for component in 0..3 {
                derivatives[component][axis] =
                    0.5 * (upper[component + 1] - lower[component + 1]) / grid.spacing[axis];
            }
        }
    }
    #[cfg(feature = "experimental-sgs")]
    if let Some(sgs) = subgrid.as_mut() {
        sgs.prepare(primitive, walls, cancelled)?;
    }
    let mut maximum_diffusivity: f64 = 0.0;
    for axis in 0..3 {
        let count = grid.shape[axis];
        let pair = boundaries.axis(axis);
        let tangential = tangents(axis);
        let flux = &mut fluxes[axis];
        let mut shape = grid.shape;
        shape[axis] += 1;
        for index in 0..flux.len() {
            poll(cancelled)?;
            let mut c = [
                index % shape[0],
                (index / shape[0]) % shape[1],
                index / (shape[0] * shape[1]),
            ];
            let face = c[axis];
            if face == count && pair == AxisBoundaries::Periodic {
                c[axis] = 0;
                flux[index] = flux[grid.cell_face(axis, c)];
                continue;
            }
            let exterior = face == 0 || face == count;
            let wall = if let AxisBoundaries::Ends { lower, upper } = pair {
                exterior && (if face == 0 { lower } else { upper }) == BoundaryCondition::Reflecting
            } else {
                false
            };
            let sign = if face == 0 { -1.0 } else { 1.0 };
            let mut adjacent = c;
            adjacent[axis] = if face == 0 { 0 } else { count - 1 };
            let adjacent = grid.index(adjacent);
            let bath = if exterior && pair != AxisBoundaries::Periodic && !wall {
                reservoirs.and_then(|r| r.get_at(grid, axis, sign, adjacent))
            } else {
                None
            };
            if exterior && pair != AxisBoundaries::Periodic && !wall && bath.is_none() {
                flux[index] = SpatialFaceFlux::default();
                continue;
            }
            let mut left = c;
            let mut right = c;
            left[axis] = if face == 0 {
                if wall || bath.is_some() { 0 } else { count - 1 }
            } else {
                face - 1
            };
            right[axis] = if face == count { count - 1 } else { face };
            let left = grid.index(left);
            let right = grid.index(right);
            let pl = if face == 0 {
                bath.map_or(primitive[left], |r| r.primitive())
            } else {
                primitive[left]
            };
            let pr = if face == count {
                bath.map_or(primitive[right], |r| r.primitive())
            } else {
                primitive[right]
            };
            let distance = if bath.is_some() {
                0.5 * grid.spacing[axis]
            } else {
                grid.spacing[axis]
            };
            let wall_model = wall.then(|| walls[2 * axis + usize::from(face == count)]);
            let temperature = match wall_model {
                Some(model) => model.face_temperature(pl.temperature),
                _ => 0.5 * pl.temperature + 0.5 * pr.temperature,
            };
            let pressure = 0.5 * pl.pressure + 0.5 * pr.pressure;
            let coefficients = law.evaluate(
                Temperature::new(temperature),
                Pressure::new(pressure),
                gas.source().mole_composition().fractions(),
                binary,
                cancelled,
            )?;
            let mu = coefficients.dynamic_viscosity.value();
            let bulk = coefficients.bulk_viscosity.value();
            let molecular_viscosity = crate::compressible::stress::Viscosity { shear: mu, bulk };
            let replace_molecular = compatible
                && !wall
                && (viscosity[left] != molecular_viscosity
                    || viscosity[right] != molecular_viscosity);
            let replace_wall = compatible && wall;
            let conductivity = coefficients.thermal_conductivity.value();
            if ![mu, bulk, conductivity]
                .into_iter()
                .all(|v| v.is_finite() && v >= 0.0)
            {
                return Err(MolecularError::Invalid(
                    "nonnegative finite molecular coefficients required",
                ));
            }
            let bath_gradient = if bath.is_some() {
                Some(reservoir::tangential_gradient(
                    grid,
                    boundaries,
                    reservoirs.ok_or(MolecularError::Invalid("missing molecular reservoir"))?,
                    axis,
                    sign,
                    adjacent,
                    &gradient[adjacent],
                )?)
            } else {
                None
            };
            #[cfg(feature = "experimental-sgs")]
            let (mu, conductivity, modeled_mu) =
                if let Some(sgs) = subgrid.as_ref().filter(|_| !wall) {
                    if let Some(derivatives) = bath_gradient {
                        sgs.reservoir_face_coefficients(
                            gas,
                            axis,
                            index,
                            pl,
                            pr,
                            distance,
                            derivatives,
                            temperature,
                            pressure,
                            mu,
                            conductivity,
                        )?
                    } else {
                        sgs.face_coefficients(
                            gas,
                            left,
                            right,
                            axis,
                            index,
                            grid.spacing[axis],
                            primitive,
                            gradient,
                            temperature,
                            pressure,
                            mu,
                            conductivity,
                        )?
                    }
                } else {
                    (mu, conductivity, 0.0)
                };
            #[cfg(not(feature = "experimental-sgs"))]
            let modeled_mu = 0.0;
            let longitudinal = (4.0 / 3.0) * mu + bulk;
            let cross = bulk - (2.0 / 3.0) * mu;
            // The wall law replaces only tangential traction. Bulk and normal
            // stresses retain the source molecular viscosity; the nonlinear
            // wall-traction Jacobian supplies its own diffusion bound.
            #[cfg(feature = "experimental-sgs")]
            let wall_viscosity = match (subgrid.as_ref(), wall_model) {
                (Some(sgs), Some(model)) => sgs.wall_viscosity(
                    2 * axis + usize::from(face == count),
                    model,
                    pl,
                    0.5 * grid.spacing[axis],
                    mu,
                )?,
                _ => sgs::wall::WallViscosity {
                    secant: mu,
                    tangent: mu,
                },
            };
            #[cfg(feature = "experimental-sgs")]
            let (wall_mu, shear_diffusion) = (wall_viscosity.secant, wall_viscosity.tangent);
            #[cfg(feature = "experimental-sgs")]
            let conductivity = match (subgrid.as_ref(), wall_model) {
                (Some(sgs), Some(model)) => sgs.wall_conductivity(
                    gas,
                    2 * axis + usize::from(face == count),
                    model,
                    pl,
                    0.5 * grid.spacing[axis],
                    mu,
                    conductivity,
                    temperature,
                    pressure,
                    cancelled,
                )?,
                _ => conductivity,
            };
            #[cfg(not(feature = "experimental-sgs"))]
            let (wall_mu, shear_diffusion) = (mu, mu);
            // Cover both tangential cross-derivative stencils, variable face
            // coefficients and doubled half-cell wall or reservoir gradients.
            let stress_scale = if replace_molecular || replace_wall {
                let modeled = crate::compressible::stress::Viscosity {
                    shear: modeled_mu,
                    bulk: 0.0,
                };
                viscosity[left]
                    .scale::<3>()
                    .max(viscosity[right].scale::<3>())
                    + modeled.scale::<3>()
            } else {
                longitudinal + 2.0 * cross.abs() + 4.0 * shear_diffusion
            };
            let diffusivity = (stress_scale / minimum_density).max(conductivity / minimum_cv);
            if !diffusivity.is_finite() {
                return Err(MolecularError::Invalid(
                    "molecular diffusion scale overflow",
                ));
            }
            maximum_diffusivity = maximum_diffusivity.max(diffusivity);
            let normal_derivative = if wall {
                let sign = if face == 0 { 1.0 } else { -1.0 };
                let wall_speed = motion.map_or(0.0, |m| m.velocity(grid, axis, face));
                if wall_speed == 0.0 {
                    sign * 2.0 * pl.velocity[axis] / grid.spacing[axis]
                } else {
                    sign * 2.0 * (pl.velocity[axis] - wall_speed) / grid.spacing[axis]
                }
            } else {
                (pr.velocity[axis] - pl.velocity[axis]) / distance
            };
            let wall_trace_fraction = wall_model.map_or(1.0, |model| {
                model.tangential_trace_fraction(0.5 * grid.spacing[axis])
            });
            let transverse_divergence: f64 = if wall_trace_fraction == 0.0 {
                // No-slip wall velocity is spatially constant. Finite-slip
                // fluid traces retain the fraction b/(h+b) of the adjacent
                // tangential gradient, including its normal-stress effect.
                0.0
            } else if let Some(derivatives) = bath_gradient {
                derivatives[0][1] + derivatives[1][2]
            } else {
                tangential
                    .into_iter()
                    .map(|t| 0.5 * gradient[left][t][t] + 0.5 * gradient[right][t][t])
                    .sum::<f64>()
                    * wall_trace_fraction
            };
            let mut normal_stress =
                longitudinal * normal_derivative + cross * transverse_divergence;
            let (mut shear_stress, mut energy) = if let Some(model) = wall_model {
                if compatible {
                    let (mut wall_velocity, slip, heat_fraction) = model.stress_boundary();
                    if let Some(motion) = motion {
                        wall_velocity[axis] = motion.velocity(grid, axis, face);
                    }
                    let stress = crate::compressible::stress::wall_stress(
                        axis,
                        sign,
                        0.5 * grid.spacing[axis],
                        pl.velocity,
                        wall_velocity,
                        slip,
                        viscosity[left],
                        |component, derivative| gradient[left][component][derivative],
                    )
                    .map_err(MolecularError::Invalid)?;
                    let traction = stress.traction;
                    normal_stress = traction[axis];
                    let shear = tangential.map(|t| traction[t]);
                    let heat = model.conductive_energy_flux(
                        sign,
                        pl.temperature,
                        0.5 * grid.spacing[axis],
                        conductivity,
                    );
                    let tangential_work = -shear[0] * wall_velocity[tangential[0]]
                        - shear[1] * wall_velocity[tangential[1]];
                    let mechanical_flux = if wall_velocity[axis] == 0.0 {
                        tangential_work
                    } else {
                        -normal_stress * wall_velocity[axis] + tangential_work
                    };
                    let thermal_flux =
                        sign * (1.0 - heat_fraction) * stress.slip_dissipation + heat;
                    if let Some(power) = wall_power.as_deref_mut() {
                        let side = 2 * axis + usize::from(face == count);
                        let inward_area = -sign * grid.face_area(BoundarySide::ALL[side]);
                        power.work[side] += inward_area * mechanical_flux;
                        power.heat[side] += inward_area * thermal_flux;
                    }
                    // Keep the stationary expression's evaluation order.
                    let mut energy = -shear[0] * wall_velocity[tangential[0]]
                        - shear[1] * wall_velocity[tangential[1]]
                        + sign * (1.0 - heat_fraction) * stress.slip_dissipation
                        + heat;
                    if wall_velocity[axis] != 0.0 {
                        energy -= normal_stress * wall_velocity[axis];
                    }
                    (shear, energy)
                } else {
                    let (shear, energy) = model.shear_and_energy(
                        axis,
                        sign,
                        pl.velocity,
                        pl.temperature,
                        0.5 * grid.spacing[axis],
                        wall_mu,
                        conductivity,
                    );
                    ([shear[tangential[0]], shear[tangential[1]]], energy)
                }
            } else {
                let shear: [f64; 2] = std::array::from_fn(|i| {
                    let t = tangential[i];
                    if let Some(derivatives) = bath_gradient {
                        mu * ((pr.velocity[t] - pl.velocity[t]) / distance + derivatives[i][0])
                    } else {
                        mu * ((pr.velocity[t] - pl.velocity[t]) / distance
                            + 0.5 * gradient[left][axis][t]
                            + 0.5 * gradient[right][axis][t])
                    }
                });
                let energy = if let Some(bath) = bath {
                    let velocity = bath.primitive().velocity;
                    -normal_stress * velocity[axis]
                        - shear[0] * velocity[tangential[0]]
                        - shear[1] * velocity[tangential[1]]
                        - conductivity * (pr.temperature - pl.temperature) / distance
                } else {
                    -normal_stress * (0.5 * pl.velocity[axis] + 0.5 * pr.velocity[axis])
                        - shear[0]
                            * (0.5 * pl.velocity[tangential[0]] + 0.5 * pr.velocity[tangential[0]])
                        - shear[1]
                            * (0.5 * pl.velocity[tangential[1]] + 0.5 * pr.velocity[tangential[1]])
                        - conductivity * (pr.temperature - pl.temperature) / grid.spacing[axis]
                };
                (shear, energy)
            };
            if replace_molecular {
                // Form molecular traction directly, retaining modeled viscosity
                // separately. Subtracting a huge old midpoint-law stress would
                // lose finite compatible stress, SGS and conductive heat.
                let traction = crate::compressible::stress::traction(
                    axis,
                    distance,
                    [pl.velocity, pr.velocity],
                    [viscosity[left], viscosity[right]],
                    |side, component, derivative| {
                        gradient[[left, right][side]][component][derivative]
                    },
                );
                normal_stress = traction[axis]
                    + modeled_mu
                        * ((4.0 / 3.0) * normal_derivative - (2.0 / 3.0) * transverse_divergence);
                for tangent in 0..2 {
                    let t = tangential[tangent];
                    shear_stress[tangent] = traction[t]
                        + modeled_mu
                            * ((pr.velocity[t] - pl.velocity[t]) / distance
                                + 0.5 * gradient[left][axis][t]
                                + 0.5 * gradient[right][axis][t]);
                }
                energy = -normal_stress * (0.5 * pl.velocity[axis] + 0.5 * pr.velocity[axis])
                    - shear_stress[0]
                        * (0.5 * pl.velocity[tangential[0]] + 0.5 * pr.velocity[tangential[0]])
                    - shear_stress[1]
                        * (0.5 * pl.velocity[tangential[1]] + 0.5 * pr.velocity[tangential[1]])
                    - conductivity * (pr.temperature - pl.temperature) / grid.spacing[axis];
            }
            if ![normal_stress, shear_stress[0], shear_stress[1], energy]
                .into_iter()
                .all(f64::is_finite)
            {
                return Err(MolecularError::Invalid("nonfinite molecular face flux"));
            }
            flux[index] = SpatialFaceFlux {
                mass: 0.0,
                normal_momentum: -normal_stress,
                transverse_momentum: shear_stress.map(|v| -v),
                energy,
                signal_speed: 0.0,
                used_hll_fallback: false,
            };
        }
    }
    // A conservative step estimate, not a fully discrete entropy theorem.
    // Original source admission is checked again after every conservative stage.
    let rate = if maximum_diffusivity == 0.0 {
        0.0
    } else {
        4.0 * maximum_diffusivity
            * grid
                .spacing
                .into_iter()
                .map(|h| h.recip().powi(2))
                .sum::<f64>()
    };
    if !rate.is_finite() {
        return Err(MolecularError::Invalid("parabolic scale overflow"));
    }
    Ok(rate)
}

#[cfg(test)]
#[path = "molecular_tests.rs"]
mod tests;
