//! Prescribed walls on the original spatial/time-dependent force transaction.
//!
//! Incoming collision uses F_n and endpoint admission uses F_(n+1). The wall
//! receipt reads the actual corrected stream before populations, source history
//! or clock publish. Consequently delta physical momentum is the trapezoidal
//! body impulse plus the receipt's fluid_population_impulse. Wall-relative mass
//! terms retain their native meaning; wall_work is work ON the prescribed wall.
//! Fixed topology, tangency and local impermeability use the native wall owner.

use super::{ForcedGrid3, ForcedLattice3};
use crate::lattice3::moving_walls::{MovingWallMomentum3, WallMotion3};
use crate::lattice3::wall_loads::{StationaryWallLink3, WallMomentumProbe3};
use crate::lattice3::{Grid3, GridError3};
use fs_exec::{CancelGate, Cx, KernelRunner};

impl<L: ForcedLattice3<Q>, const Q: usize> ForcedGrid3<L, Q> {
    /// Compile native wall links without exposing populations without their
    /// accepted force frame. Geometry and selection match Grid3 exactly.
    pub fn wall_momentum_probe(
        &self,
        cx: &Cx<'_>,
        origin: [f64; 3],
        select: impl FnMut(StationaryWallLink3) -> bool,
    ) -> Result<WallMomentumProbe3<L, Q>, GridError3<L::Error>> {
        self.grid.wall_momentum_probe(cx, origin, select)
    }

    /// Hold current body forces while applying prescribed wall velocities and
    /// the original per-cell molecular or prepared SGS collision field.
    pub fn step_serial_with_moving_wall_momentum(
        &mut self,
        cx: &Cx<'_>,
        probe: &WallMomentumProbe3<L, Q>,
        velocities: &[[f64; 3]],
        collisions: &[L::Collision],
    ) -> Result<MovingWallMomentum3, GridError3<L::Error>> {
        serial(
            &mut self.grid, &self.forces, &self.forces, collisions, cx, probe, velocities,
        )
    }

    /// Advance the body-force endpoint and wall-driven fluid together.
    /// No accepted state changes if source, wall, stream or receipt admission
    /// fails. Wall velocities are held over this one lattice step.
    pub fn advance_serial_with_moving_wall_momentum(
        &mut self,
        cx: &Cx<'_>,
        probe: &WallMomentumProbe3<L, Q>,
        velocities: &[[f64; 3]],
        collisions: &[L::Collision],
        forces_at_end: Vec<[f64; 3]>,
    ) -> Result<MovingWallMomentum3, GridError3<L::Error>> {
        let receipt = serial(
            &mut self.grid, &self.forces, &forces_at_end, collisions, cx, probe, velocities,
        )?;
        self.forces = forces_at_end;
        Ok(receipt)
    }

    /// Pooled native moving-wall step with the accepted spatial force held.
    pub fn step_pooled_with_moving_wall_momentum<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        probe: &WallMomentumProbe3<L, Q>,
        velocities: &[[f64; 3]],
        collisions: &[L::Collision],
    ) -> Result<MovingWallMomentum3, GridError3<L::Error>> {
        pooled(
            &mut self.grid, &self.forces, &self.forces, collisions, pool, gate, probe, velocities,
        )
    }

    /// Pooled endpoint advance; source history publishes only after the same
    /// candidate stream and canonical-order wall receipt have both succeeded.
    pub fn advance_pooled_with_moving_wall_momentum<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        probe: &WallMomentumProbe3<L, Q>,
        velocities: &[[f64; 3]],
        collisions: &[L::Collision],
        forces_at_end: Vec<[f64; 3]>,
    ) -> Result<MovingWallMomentum3, GridError3<L::Error>> {
        let receipt = pooled(
            &mut self.grid, &self.forces, &forces_at_end, collisions, pool, gate, probe, velocities,
        )?;
        self.forces = forces_at_end;
        Ok(receipt)
    }
}

fn serial<L: ForcedLattice3<Q>, const Q: usize>(
    grid: &mut Grid3<L, Q>,
    start: &[[f64; 3]],
    end: &[[f64; 3]],
    collisions: &[L::Collision],
    cx: &Cx<'_>,
    probe: &WallMomentumProbe3<L, Q>,
    velocities: &[[f64; 3]],
) -> Result<MovingWallMomentum3, GridError3<L::Error>> {
    let poll = || cx.checkpoint().map_err(|_| GridError3::Cancelled);
    let motion = WallMotion3 { probe, velocities };
    let moving = motion.admit(grid, poll)?;
    super::serial_with_wall_motion(
        grid, start, end, collisions, cx, moving.then_some(&motion),
        |candidate| motion.measure(candidate, poll),
    )
}

#[allow(clippy::too_many_arguments)]
fn pooled<L: ForcedLattice3<Q>, const Q: usize, P: KernelRunner>(
    grid: &mut Grid3<L, Q>,
    start: &[[f64; 3]],
    end: &[[f64; 3]],
    collisions: &[L::Collision],
    pool: &P,
    gate: &CancelGate,
    probe: &WallMomentumProbe3<L, Q>,
    velocities: &[[f64; 3]],
) -> Result<MovingWallMomentum3, GridError3<L::Error>> {
    let poll = || {
        if gate.is_requested() {
            Err(GridError3::Cancelled)
        } else {
            Ok(())
        }
    };
    let motion = WallMotion3 { probe, velocities };
    let moving = motion.admit(grid, poll)?;
    super::pooled_with_wall_motion(
        grid, start, end, collisions, pool, gate, moving.then_some(&motion),
        |candidate| motion.measure(candidate, poll),
    )
}
