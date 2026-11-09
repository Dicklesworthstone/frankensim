//! Prescribed acceleration on the original native collision/transport owner.

use super::{ForcedGrid3, ForcedLattice3, admit_forces};
use crate::lattice3::moving_walls::{MovingWallMomentum3, WallMotion3};
use crate::lattice3::wall_loads::WallMomentumProbe3;
use crate::lattice3::{EquilibriumLattice3, Grid3, GridError3};
use fs_exec::{CancelGate, Cx, KernelRunner};

pub(super) enum EndpointForcing<'a> {
    Density(&'a [[f64; 3]]),
    Acceleration(&'a mut [[f64; 3]]),
}

impl<'a> EndpointForcing<'a> {
    pub(super) fn values(&self) -> &[[f64; 3]] {
        match self {
            Self::Density(values) => values,
            Self::Acceleration(values) => values,
        }
    }

    pub(super) fn split(self) -> (Option<&'a [[f64; 3]]>, Option<&'a mut [[f64; 3]]>) {
        match self {
            Self::Density(values) => (Some(values), None),
            Self::Acceleration(values) => (None, Some(values)),
        }
    }
}

impl<L: EquilibriumLattice3<Q>, const Q: usize> ForcedGrid3<L, Q> {
    /// Initialize a spatial acceleration field without an artificial half kick.
    /// Values are lattice acceleration (`a_SI * dt² / dx`) in x-fastest order;
    /// solid entries must be zero. The initial force is the actual local `rho*a`.
    /// The original native source transform preserves velocity and stress.
    pub fn try_new_with_body_acceleration_preserving_velocity(
        grid: Grid3<L, Q>,
        cx: &Cx<'_>,
        mut acceleration: Vec<[f64; 3]>,
    ) -> Result<Self, GridError3<L::Error>> {
        let mut poll = || cx.checkpoint().map_err(|_| GridError3::Cancelled);
        admit_forces(&grid, &acceleration, None, &mut poll)?;
        for (cell, value) in acceleration.iter_mut().enumerate() {
            if cell % 64 == 0 {
                poll()?;
            }
            if !grid.solid[cell] {
                let rho: f64 = grid.f[cell].iter().sum();
                *value = value.map(|a| rho * a);
            }
        }
        Self::try_new_preserving_velocity(grid, cx, acceleration)
    }
}

impl<L: ForcedLattice3<Q>, const Q: usize> ForcedGrid3<L, Q> {
    /// Advance with prescribed spatial acceleration at the next lattice time.
    /// Collision retains `F_n`; the completed native stream determines
    /// `F_(n+1)(x) = rho_(n+1)(x) * a_(n+1)(x)`. Physical momentum therefore
    /// receives the trapezoidal body impulse plus native boundary transport.
    /// Supply each endpoint explicitly: ordinary `step_serial` holds force
    /// density, not acceleration. No same-time rebind or extra kick occurs.
    pub fn advance_serial_with_body_acceleration(
        &mut self,
        cx: &Cx<'_>,
        collisions: &[L::Collision],
        mut acceleration_at_end: Vec<[f64; 3]>,
    ) -> Result<(), GridError3<L::Error>> {
        super::serial_with_wall_motion(
            &mut self.grid,
            &self.forces,
            EndpointForcing::Acceleration(&mut acceleration_at_end),
            collisions,
            cx,
            None,
            |_| Ok(()),
        )?;
        self.forces = acceleration_at_end;
        Ok(())
    }

    /// Pooled counterpart; the same transported density and canonical endpoint
    /// admission produce identical populations and source history for any pool.
    pub fn advance_pooled_with_body_acceleration<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        collisions: &[L::Collision],
        mut acceleration_at_end: Vec<[f64; 3]>,
    ) -> Result<(), GridError3<L::Error>> {
        super::pooled_with_wall_motion(
            &mut self.grid,
            &self.forces,
            EndpointForcing::Acceleration(&mut acceleration_at_end),
            collisions,
            pool,
            gate,
            None,
            |_| Ok(()),
        )?;
        self.forces = acceleration_at_end;
        Ok(())
    }

    /// Couple the endpoint acceleration to prescribed native wall velocities.
    /// The receipt measures the actual corrected stream; wall impulse and work
    /// retain their native meaning. Source, populations and clock publish only
    /// after endpoint Mach, wall receipt and final cancellation admission.
    pub fn advance_serial_with_body_acceleration_and_moving_wall_momentum(
        &mut self,
        cx: &Cx<'_>,
        probe: &WallMomentumProbe3<L, Q>,
        velocities: &[[f64; 3]],
        collisions: &[L::Collision],
        mut acceleration_at_end: Vec<[f64; 3]>,
    ) -> Result<MovingWallMomentum3, GridError3<L::Error>> {
        let poll = || cx.checkpoint().map_err(|_| GridError3::Cancelled);
        let motion = WallMotion3 { probe, velocities };
        let moving = motion.admit(&self.grid, poll)?;
        let receipt = super::serial_with_wall_motion(
            &mut self.grid,
            &self.forces,
            EndpointForcing::Acceleration(&mut acceleration_at_end),
            collisions,
            cx,
            moving.then_some(&motion),
            |grid| motion.measure(grid, poll),
        )?;
        self.forces = acceleration_at_end;
        Ok(receipt)
    }

    /// Executor counterpart using the same two native passes and wall receipt.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_pooled_with_body_acceleration_and_moving_wall_momentum<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        probe: &WallMomentumProbe3<L, Q>,
        velocities: &[[f64; 3]],
        collisions: &[L::Collision],
        mut acceleration_at_end: Vec<[f64; 3]>,
    ) -> Result<MovingWallMomentum3, GridError3<L::Error>> {
        let poll = || {
            if gate.is_requested() {
                Err(GridError3::Cancelled)
            } else {
                Ok(())
            }
        };
        let motion = WallMotion3 { probe, velocities };
        let moving = motion.admit(&self.grid, poll)?;
        let receipt = super::pooled_with_wall_motion(
            &mut self.grid,
            &self.forces,
            EndpointForcing::Acceleration(&mut acceleration_at_end),
            collisions,
            pool,
            gate,
            moving.then_some(&motion),
            |grid| motion.measure(grid, poll),
        )?;
        self.forces = acceleration_at_end;
        Ok(receipt)
    }
}
