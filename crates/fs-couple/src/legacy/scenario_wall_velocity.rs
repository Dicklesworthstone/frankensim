//! Analytic scenario-frame drives for the native fixed-topology moving walls.
//!
//! Each native wall link is an Eulerian location. Its velocity is evaluated at
//! that fixed WORLD point, rather than at a material point advected from the
//! initial frame. One analytic frame sample serves every link at a given time.
//! The original LBM owner admits tangency, discrete impermeability and Mach
//! limits, applies its moving bounce-back law, and publishes the complete step.

use fs_exec::{CancelGate, Cx, KernelRunner};
use fs_ga::Vec3;
use fs_lbm::lattice3::moving_walls::MovingWallMomentum3;
use fs_lbm::lattice3::wall_loads::WallMomentumProbe3;
use fs_lbm::lattice3::{Grid3, GridError3, Lattice3};
use fs_lbm::material::{MaterialLbmError, MaterialScaling};
use fs_material::fluid_transport::FluidTransport;
use fs_scenario::frame::{FrameKinematics, FrameKinematicsError, PreparedFrameKinematics};
use fs_scenario::signal::DEFAULT_MAXIMUM_CYCLES_PER_STEP;
use std::ops::ControlFlow;

/// Refused physical wall-velocity sampling; no fluid step has begun.
#[derive(Debug)]
pub enum FrameWallVelocityError {
    /// Invalid fixed grid origin, buffer size, or position arithmetic.
    Invalid(&'static str),
    /// Cancellation during sampling.
    Cancelled,
    /// The original frame owner refused its history or kinematics.
    Frame(FrameKinematicsError),
    /// SI-to-lattice conversion failed.
    Scaling(MaterialLbmError),
    /// A complete wall velocity exceeds the supplied scaling's speed bound.
    SpeedBound {
        /// Index in the original wall probe order.
        link: usize,
        /// Euclidean speed in lattice units.
        actual: f64,
        /// Admitted scenario speed in lattice units.
        maximum: f64,
    },
    /// The selected frame source cannot be resolved by the fixed lattice step.
    /// Rebuild the material scaling and grid with a smaller physical timestep.
    FrequencyResolution {
        /// Caller-selected physical sample time [s].
        time_s: f64,
        /// Frozen physical duration of one native lattice step [s].
        timestep_s: f64,
    },
}

impl core::fmt::Display for FrameWallVelocityError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        if let Self::FrequencyResolution { time_s, timestep_s } = self {
            return write!(
                f,
                "wall frame source at {time_s} s is unresolved by the fixed lattice timestep {timestep_s} s; rebuild material scaling and the grid with a smaller timestep"
            );
        }
        write!(f, "scenario wall velocity refused: {self:?}")
    }
}
impl std::error::Error for FrameWallVelocityError {}

/// Refusal before the original fluid owner publishes its complete step.
#[derive(Debug)]
pub enum FrameWallStepError<E> {
    /// Sampling failed before entering the native operator.
    Velocity(FrameWallVelocityError),
    /// Native topology, tangency, Mach, collision, receipt or cancellation refusal.
    Native(GridError3<E>),
}

impl<E: core::fmt::Debug> core::fmt::Display for FrameWallStepError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        if let Self::Velocity(error) = self {
            return write!(f, "scenario-driven wall step refused: {error}");
        }
        write!(f, "scenario-driven wall step refused: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for FrameWallStepError<E> {}

/// Borrowed physical drive, prepared frame path and frozen material scaling.
///
/// The grid's Cartesian axes are WORLD axes. `grid_origin_world_m` locates its
/// zero cell-coordinate corner; the scaling must be the one used for the grid.
/// No histories, wall topology or material laws are cloned or re-evaluated.
/// The caller provides one scratch velocity per selected native wall link.
pub struct FrameWallVelocity3<'a, 'frame, 'material, F: FluidTransport + ?Sized> {
    frame: &'a PreparedFrameKinematics<'frame>,
    scaling: &'a MaterialScaling<'material, F>,
    grid_origin_world_m: [f64; 3],
}

impl<'a, 'frame, 'material, F: FluidTransport + ?Sized>
    FrameWallVelocity3<'a, 'frame, 'material, F>
{
    /// Bind one analytic drive to an explicitly located and scaled grid.
    pub fn try_new(
        frame: &'a PreparedFrameKinematics<'frame>,
        scaling: &'a MaterialScaling<'material, F>,
        grid_origin_world_m: [f64; 3],
    ) -> Result<Self, FrameWallVelocityError> {
        if grid_origin_world_m.iter().any(|v| !v.is_finite()) {
            return Err(FrameWallVelocityError::Invalid("nonfinite grid origin"));
        }
        Ok(Self {
            frame,
            scaling,
            grid_origin_world_m,
        })
    }

    /// Fill caller-owned scratch in the native probe order, in lattice units.
    ///
    /// Frame sampling occurs once; link work is allocation-free and polls every
    /// 64 links. An error may leave scratch partly filled, but never changes
    /// the grid, its populations, its targets or its step count. Scratch is not
    /// an accepted control. The native step remains responsible for geometric
    /// compatibility, local mass conservation and its own collision limits.
    pub fn sample_velocities<L: Lattice3<Q>, const Q: usize>(
        &self,
        time_s: f64,
        probe: &WallMomentumProbe3<L, Q>,
        scratch: &mut [[f64; 3]],
        checkpoint: &mut impl FnMut() -> ControlFlow<()>,
    ) -> Result<FrameKinematics, FrameWallVelocityError> {
        poll(checkpoint)?;
        if scratch.len() != probe.measured_links() {
            return Err(FrameWallVelocityError::Invalid(
                "wall velocity scratch length",
            ));
        }
        let frame = self
            .frame
            .sample(time_s, checkpoint)
            .map_err(FrameWallVelocityError::Frame)?;
        let dx = self.scaling.cell_size_m();
        let maximum = self.scaling.flow().u_lattice;
        for (index, (link, velocity)) in probe.wall_links().zip(scratch.iter_mut()).enumerate() {
            if index % 64 == 0 {
                poll(checkpoint)?;
            }
            let point: [f64; 3] =
                core::array::from_fn(|a| link.position[a].mul_add(dx, self.grid_origin_world_m[a]));
            if point.iter().any(|v| !v.is_finite()) {
                return Err(FrameWallVelocityError::Invalid(
                    "wall point overflow in SI coordinates",
                ));
            }
            let si = frame
                .velocity_at_world_point(Vec3::new(point[0], point[1], point[2]))
                .map_err(FrameWallVelocityError::Frame)?;
            let mut lattice = [0.0; 3];
            for (out, value) in lattice.iter_mut().zip([si.x, si.y, si.z]) {
                *out = self
                    .scaling
                    .velocity_lattice(value)
                    .map_err(FrameWallVelocityError::Scaling)?;
            }
            let actual = lattice[0].hypot(lattice[1]).hypot(lattice[2]);
            if actual > maximum {
                return Err(FrameWallVelocityError::SpeedBound {
                    link: index,
                    actual,
                    maximum,
                });
            }
            *velocity = lattice;
        }
        poll(checkpoint)?;
        Ok(frame)
    }

    fn admit_step_frequency(
        &self,
        time_s: f64,
        checkpoint: &mut impl FnMut() -> ControlFlow<()>,
    ) -> Result<(), FrameWallVelocityError> {
        poll(checkpoint)?;
        let timestep_s = self.scaling.timestep_s();
        if !self
            .frame
            .admits_sampled_step(
                time_s,
                timestep_s,
                DEFAULT_MAXIMUM_CYCLES_PER_STEP,
                checkpoint,
            )
            .map_err(FrameWallVelocityError::Frame)?
        {
            return Err(FrameWallVelocityError::FrequencyResolution { time_s, timestep_s });
        }
        Ok(())
    }

    /// Sample the held drive and perform one original serial lattice step.
    ///
    /// Populations and step count publish only after the native moving-wall
    /// receipt is accepted. The returned receipt retains its original lattice
    /// units, including wall impulse and mechanical work. This prescribes
    /// surface motion on fixed topology; normal/swept geometry still refuses.
    /// Choose the sample time explicitly (for example the interval midpoint).
    /// Unresolved native chirps in the prepared frame path refuse before scratch
    /// sampling. Use a smaller material-scaling timestep to resolve the source;
    /// the fluid timestep and transport scaling remain frozen during the step.
    pub fn step_serial<L: Lattice3<Q>, const Q: usize>(
        &self,
        grid: &mut Grid3<L, Q>,
        cx: &Cx<'_>,
        probe: &WallMomentumProbe3<L, Q>,
        scratch: &mut [[f64; 3]],
        time_s: f64,
        collisions: Option<&[L::Collision]>,
    ) -> Result<MovingWallMomentum3, FrameWallStepError<L::Error>> {
        let mut checkpoint = || {
            if cx.checkpoint().is_err() {
                ControlFlow::Break(())
            } else {
                ControlFlow::Continue(())
            }
        };
        self.admit_step_frequency(time_s, &mut checkpoint)
            .map_err(FrameWallStepError::Velocity)?;
        self.sample_velocities(time_s, probe, scratch, &mut checkpoint)
            .map_err(FrameWallStepError::Velocity)?;
        grid.step_serial_with_moving_wall_momentum(cx, probe, scratch, collisions)
            .map_err(FrameWallStepError::Native)
    }

    /// Executor counterpart using the same sampling and native publication law.
    pub fn step_pooled<L: Lattice3<Q>, const Q: usize, P: KernelRunner>(
        &self,
        grid: &mut Grid3<L, Q>,
        pool: &P,
        gate: &CancelGate,
        probe: &WallMomentumProbe3<L, Q>,
        scratch: &mut [[f64; 3]],
        time_s: f64,
        collisions: Option<&[L::Collision]>,
    ) -> Result<MovingWallMomentum3, FrameWallStepError<L::Error>> {
        let mut checkpoint = || {
            if gate.is_requested() {
                ControlFlow::Break(())
            } else {
                ControlFlow::Continue(())
            }
        };
        self.admit_step_frequency(time_s, &mut checkpoint)
            .map_err(FrameWallStepError::Velocity)?;
        self.sample_velocities(time_s, probe, scratch, &mut checkpoint)
            .map_err(FrameWallStepError::Velocity)?;
        grid.step_pooled_with_moving_wall_momentum(pool, gate, probe, scratch, collisions)
            .map_err(FrameWallStepError::Native)
    }
}

fn poll(checkpoint: &mut impl FnMut() -> ControlFlow<()>) -> Result<(), FrameWallVelocityError> {
    if checkpoint().is_break() {
        Err(FrameWallVelocityError::Cancelled)
    } else {
        Ok(())
    }
}
