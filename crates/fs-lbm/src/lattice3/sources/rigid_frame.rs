//! Full prescribed rigid-frame evolution on the original native source owner.

use super::acceleration::EndpointForcing;
use super::{ForcedGrid3, ForcedLattice3, admit_forces};
use crate::d3q19::{Coriolis3, RigidFrame3};
use crate::lattice3::moving_walls::{MovingWallMomentum3, WallMotion3};
use crate::lattice3::wall_loads::WallMomentumProbe3;
use crate::lattice3::{EquilibriumLattice3, Grid3, GridError3, coordinates, macro_cell};
use fs_exec::{CancelGate, Cx, KernelRunner};

/// One admitted skew solve for local Coriolis or full rigid-frame endpoints.
#[derive(Clone, Copy)]
pub(in crate::lattice3) struct FrameEndpoint3 {
    frame: RigidFrame3,
    denominator: f64,
    spatial: bool,
}

impl FrameEndpoint3 {
    pub(super) fn coriolis(rotation: Coriolis3) -> Result<Option<Self>, &'static str> {
        let frame = rotation.frame();
        let denominator = frame.midpoint_denominator([0.0; 3])?;
        Ok((rotation.angular_velocity != [0.0; 3]).then_some(Self {
            frame,
            denominator,
            spatial: false,
        }))
    }

    pub(super) fn rigid(
        frame: RigidFrame3,
        periodic: [bool; 3],
    ) -> Result<Option<Self>, &'static str> {
        let denominator = frame.midpoint_denominator([0.0; 3])?;
        for (axis, periodic) in periodic.into_iter().enumerate() {
            if periodic {
                frame.admit_periodic_axis(axis)?;
            }
        }
        Ok((frame.angular_velocity != [0.0; 3]
            || frame.angular_acceleration != [0.0; 3]
            || frame.origin_acceleration != [0.0; 3])
            .then_some(Self {
                frame,
                denominator,
                spatial: true,
            }))
    }

    fn base(self, cell: usize, dims: [usize; 3], acceleration: [f64; 3]) -> [f64; 3] {
        if self.spatial {
            let position = coordinates(cell, dims).map(|coordinate| coordinate as f64 + 0.5);
            self.frame.position_acceleration(position, acceleration)
        } else {
            acceleration
        }
    }

    pub(super) fn acceleration(self, velocity: [f64; 3], base: [f64; 3]) -> [f64; 3] {
        let omega = self.frame.angular_velocity;
        core::array::from_fn(|a| {
            let b = (a + 1) % 3;
            let c = (a + 2) % 3;
            base[a] - 2.0 * omega[b].mul_add(velocity[c], -omega[c] * velocity[b])
        })
    }

    pub(in crate::lattice3) fn force(
        self,
        rho: f64,
        momentum: [f64; 3],
        acceleration: [f64; 3],
        cell: usize,
        dims: [usize; 3],
    ) -> Option<[f64; 3]> {
        let base = self.base(cell, dims, acceleration);
        // rho*u = m_raw + rho*(base - 2*Omega cross u)/2, hence
        // (I + [Omega]x)u = m_raw/rho + base/2. The existing analytic
        // skew inverse closes this in the actual streamed density and velocity.
        let effective = self.frame.source_acceleration::<true>(
            [0.0; 3],
            momentum.map(|value| value / rho),
            base,
            self.denominator,
        );
        let force = effective.map(|value| rho * value);
        let velocity = core::array::from_fn(|a| (momentum[a] + 0.5 * force[a]) / rho);
        let actual = self.acceleration(velocity, base);
        let omega = self.frame.angular_velocity;
        for a in 0..3 {
            let b = (a + 1) % 3;
            let c = (a + 2) % 3;
            let scale = force[a].abs()
                + rho
                    * (base[a].abs()
                        + 2.0 * ((omega[b] * velocity[c]).abs() + (omega[c] * velocity[b]).abs()));
            if !force[a].is_finite()
                || !actual[a].is_finite()
                || !scale.is_finite()
                || (force[a] - rho * actual[a]).abs() > 128.0 * f64::EPSILON * scale
            {
                return None;
            }
        }
        Some(force)
    }
}

impl<L: EquilibriumLattice3<Q>, const Q: usize> ForcedGrid3<L, Q> {
    /// Attach the full instantaneous rigid-frame force without a startup kick.
    ///
    /// The input zero-force grid already contains velocities RELATIVE to the
    /// declared frame. Physical body acceleration, including gravity, is also
    /// resolved in those frame axes and supplied in x-fastest cell order.
    /// Solid entries must be zero. All quantities use lattice units, dx=dt=1.
    /// At cell centers r = [x+0.5,y+0.5,z+0.5] - frame.origin, this attaches
    /// F = rho*(g - a_origin - alpha cross r - Omega cross (Omega cross r)
    ///          - 2*Omega cross u).
    ///
    /// The original source transform preserves physical velocity and corrected
    /// stress. It neither converts an inertial initial condition nor advances
    /// time. Periodic axes obey the same admission as RigidFrame3::prepare_step.
    pub fn try_new_with_rigid_frame_preserving_velocity(
        grid: Grid3<L, Q>,
        cx: &Cx<'_>,
        frame: RigidFrame3,
        mut acceleration: Vec<[f64; 3]>,
    ) -> Result<Self, GridError3<L::Error>> {
        cx.checkpoint().map_err(|_| GridError3::Cancelled)?;
        let Some(law) = FrameEndpoint3::rigid(frame, grid.periodic).map_err(GridError3::Input)?
        else {
            return Self::try_new_with_body_acceleration_preserving_velocity(grid, cx, acceleration);
        };
        let mut poll = || cx.checkpoint().map_err(|_| GridError3::Cancelled);
        admit_forces(&grid, &acceleration, None, &mut poll)?;
        for (cell, value) in acceleration.iter_mut().enumerate() {
            if cell % 64 == 0 {
                poll()?;
            }
            if !grid.solid[cell] {
                let (rho, velocity) =
                    macro_cell::<L, Q>(&grid.f[cell], [0.0; 3], cell, grid.budget.maximum_mach)?;
                let base = law.base(cell, grid.dims, *value);
                *value = law.acceleration(velocity, base).map(|a| rho * a);
            }
        }
        Self::try_new_preserving_velocity(grid, cx, acceleration)
    }
}

impl<L: ForcedLattice3<Q>, const Q: usize> ForcedGrid3<L, Q> {
    /// Advance the complete prescribed rigid-frame law in the native stream.
    ///
    /// Collision retains the accepted F_n. After transport and wall reflection,
    /// each original tile evaluates the next position-dependent terms at its
    /// cell centers and closes Coriolis force with the actual rho_next and
    /// half-force velocity. The resulting F_next replaces the supplied vector
    /// in place. There is no extra population buffer, source sweep or iteration.
    ///
    /// Supply the frame motion and physical body acceleration at the NEXT
    /// endpoint. The caller owns their time history and consistent frame-axis
    /// components, including alpha=dOmega/dt. No inertial coordinate remap is
    /// performed. Grid and wall velocities are relative to this moving frame.
    /// Stationary native walls therefore describe a vessel at rest in it.
    ///
    /// A periodic direction requires both Omega and alpha parallel to that
    /// axis; two periodic directions require zero rotation and angular
    /// acceleration. Translation is compatible with every periodic direction.
    /// Zero frame motion follows the original body-acceleration path exactly.
    ///
    /// Physical momentum receives the endpoint trapezoidal body impulse plus
    /// native transport. Coriolis endpoint power is zero; the other frame terms
    /// exchange work. This is not an inertial/full-flow energy balance or a
    /// finite-inertia vessel model. Resolve the prescribed motion and retain
    /// the original low-Mach/collision regime. Failure publishes no population,
    /// force history or clock. Ordinary step_serial holds accepted force density.
    pub fn advance_serial_with_rigid_frame(
        &mut self,
        cx: &Cx<'_>,
        collisions: &[L::Collision],
        frame_at_end: RigidFrame3,
        mut acceleration_at_end: Vec<[f64; 3]>,
    ) -> Result<(), GridError3<L::Error>> {
        cx.checkpoint().map_err(|_| GridError3::Cancelled)?;
        let endpoint = EndpointForcing::with_rigid_frame(
            &mut acceleration_at_end, frame_at_end, self.grid.periodic,
        )
        .map_err(GridError3::Input)?;
        super::serial_with_wall_motion(
            &mut self.grid, &self.forces, endpoint, collisions, cx, None, |_| Ok(()),
        )?;
        self.forces = acceleration_at_end;
        Ok(())
    }

    /// The same full frame closure in the original executor stream tiles.
    pub fn advance_pooled_with_rigid_frame<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        collisions: &[L::Collision],
        frame_at_end: RigidFrame3,
        mut acceleration_at_end: Vec<[f64; 3]>,
    ) -> Result<(), GridError3<L::Error>> {
        if gate.is_requested() {
            return Err(GridError3::Cancelled);
        }
        let endpoint = EndpointForcing::with_rigid_frame(
            &mut acceleration_at_end, frame_at_end, self.grid.periodic,
        )
        .map_err(GridError3::Input)?;
        super::pooled_with_wall_motion(
            &mut self.grid, &self.forces, endpoint, collisions, pool, gate, None, |_| Ok(()),
        )?;
        self.forces = acceleration_at_end;
        Ok(())
    }

    /// Apply full frame acceleration and relative moving-wall reflection in one
    /// transaction. The receipt is the original population impulse, torque and
    /// relative wall work; no fictitious frame force is counted as wall recoil.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_serial_with_rigid_frame_and_moving_wall_momentum(
        &mut self,
        cx: &Cx<'_>,
        probe: &WallMomentumProbe3<L, Q>,
        velocities: &[[f64; 3]],
        collisions: &[L::Collision],
        frame_at_end: RigidFrame3,
        mut acceleration_at_end: Vec<[f64; 3]>,
    ) -> Result<MovingWallMomentum3, GridError3<L::Error>> {
        let poll = || cx.checkpoint().map_err(|_| GridError3::Cancelled);
        let motion = WallMotion3 { probe, velocities };
        let moving = motion.admit(&self.grid, poll)?;
        let endpoint = EndpointForcing::with_rigid_frame(
            &mut acceleration_at_end, frame_at_end, self.grid.periodic,
        )
        .map_err(GridError3::Input)?;
        let receipt = super::serial_with_wall_motion(
            &mut self.grid,
            &self.forces,
            endpoint,
            collisions,
            cx,
            moving.then_some(&motion),
            |grid| motion.measure(grid, poll),
        )?;
        self.forces = acceleration_at_end;
        Ok(receipt)
    }

    /// Pooled full-frame counterpart with identical source and wall ownership.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_pooled_with_rigid_frame_and_moving_wall_momentum<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        probe: &WallMomentumProbe3<L, Q>,
        velocities: &[[f64; 3]],
        collisions: &[L::Collision],
        frame_at_end: RigidFrame3,
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
        let endpoint = EndpointForcing::with_rigid_frame(
            &mut acceleration_at_end, frame_at_end, self.grid.periodic,
        )
        .map_err(GridError3::Input)?;
        let receipt = super::pooled_with_wall_motion(
            &mut self.grid,
            &self.forces,
            endpoint,
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

#[cfg(test)]
mod tests;
