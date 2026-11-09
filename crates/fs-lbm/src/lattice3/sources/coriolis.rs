//! Velocity-dependent rotating-flow sources on the original native stream.

use super::acceleration::EndpointForcing;
use super::{ForcedGrid3, ForcedLattice3, admit_forces};
use crate::d3q19::Coriolis3;
use crate::lattice3::moving_walls::{MovingWallMomentum3, WallMotion3};
use crate::lattice3::wall_loads::WallMomentumProbe3;
use crate::lattice3::{EquilibriumLattice3, Grid3, GridError3, macro_cell};
use fs_exec::{CancelGate, Cx, KernelRunner};

/// Coefficients are admitted once; every native tile uses the same skew solve.
#[derive(Clone, Copy)]
pub(in crate::lattice3) struct CoriolisEndpoint3 {
    rotation: Coriolis3,
    denominator: f64,
}

impl CoriolisEndpoint3 {
    pub(super) fn new(rotation: Coriolis3) -> Result<Option<Self>, &'static str> {
        let denominator = rotation.frame().midpoint_denominator([0.0; 3])?;
        Ok((rotation.angular_velocity != [0.0; 3]).then_some(Self {
            rotation,
            denominator,
        }))
    }

    fn acceleration(self, velocity: [f64; 3], base: [f64; 3]) -> [f64; 3] {
        let omega = self.rotation.angular_velocity;
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
    ) -> Option<[f64; 3]> {
        // rho*u = m_raw + rho*(a - 2*Omega cross u)/2. Thus
        // (I + [Omega]x)u = m_raw/rho + a/2: exactly the existing analytic
        // midpoint skew inverse, here applied to the streamed raw velocity.
        let effective = self.rotation.frame().source_acceleration::<true>(
            [0.0; 3],
            momentum.map(|value| value / rho),
            acceleration,
            self.denominator,
        );
        let force = effective.map(|value| rho * value);
        let velocity = core::array::from_fn(|a| (momentum[a] + 0.5 * force[a]) / rho);
        let actual = self.acceleration(velocity, acceleration);
        let omega = self.rotation.angular_velocity;
        for a in 0..3 {
            let b = (a + 1) % 3;
            let c = (a + 2) % 3;
            let scale = force[a].abs()
                + rho
                    * (acceleration[a].abs()
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
    /// Attach `F = rho*(a - 2*Omega cross u)` while preserving the input
    /// physical velocity and native force-corrected stress. The input is an
    /// unstepped zero-force grid. `acceleration` is the non-Coriolis lattice
    /// acceleration in x-fastest order; solid entries must be zero.
    /// This changes only the population representation, without a half kick.
    pub fn try_new_with_coriolis_preserving_velocity(
        grid: Grid3<L, Q>,
        cx: &Cx<'_>,
        rotation: Coriolis3,
        mut acceleration: Vec<[f64; 3]>,
    ) -> Result<Self, GridError3<L::Error>> {
        cx.checkpoint().map_err(|_| GridError3::Cancelled)?;
        let Some(law) = CoriolisEndpoint3::new(rotation).map_err(GridError3::Input)? else {
            return Self::try_new_with_body_acceleration_preserving_velocity(
                grid,
                cx,
                acceleration,
            );
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
                *value = law.acceleration(velocity, *value).map(|a| rho * a);
            }
        }
        Self::try_new_preserving_velocity(grid, cx, acceleration)
    }
}

impl<L: ForcedLattice3<Q>, const Q: usize> ForcedGrid3<L, Q> {
    /// Evolve local rotating flow in one original collision/stream transaction.
    ///
    /// Collision uses the accepted `F_n`. After native transport, an analytic
    /// skew solve closes `F_next = rho_next*(a_next - 2*Omega cross u_next)` and
    /// the physical half-force velocity together. No stale velocity, reference
    /// density, explicit Coriolis kick or zero-time force reframe is used.
    /// `acceleration_at_end` becomes the next accepted force storage in place.
    /// Serial and pooled entries retain the original two native passes.
    ///
    /// For fixed Omega, isolated constant-density flow has the Cayley rotation
    /// angle `-2*atan(|Omega|)` per tick and preserves speed to roundoff. The
    /// Coriolis endpoint power is zero; this is not a transported fluid energy
    /// balance. Resolve the inertial period to control phase error. A changing
    /// Omega alone is not a complete accelerating frame: Euler, centrifugal
    /// and translational terms must be represented in the prescribed `a`.
    /// Velocities use the declared grid axes, without an inertial/frame remap.
    ///
    /// Initialize the accepted force with the same physical law, then supply
    /// each next endpoint explicitly. Ordinary `step_serial` holds force density.
    /// Zero spin retains the original body-acceleration arithmetic. Source
    /// representation, collision/Mach or cancellation refusal changes no
    /// accepted populations, force history or step count.
    pub fn advance_serial_with_coriolis(
        &mut self,
        cx: &Cx<'_>,
        collisions: &[L::Collision],
        rotation: Coriolis3,
        mut acceleration_at_end: Vec<[f64; 3]>,
    ) -> Result<(), GridError3<L::Error>> {
        cx.checkpoint().map_err(|_| GridError3::Cancelled)?;
        let endpoint = EndpointForcing::with_coriolis(&mut acceleration_at_end, rotation)
            .map_err(GridError3::Input)?;
        super::serial_with_wall_motion(
            &mut self.grid,
            &self.forces,
            endpoint,
            collisions,
            cx,
            None,
            |_| Ok(()),
        )?;
        self.forces = acceleration_at_end;
        Ok(())
    }

    /// Same implicit physical endpoint on the original executor stream tiles.
    pub fn advance_pooled_with_coriolis<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        collisions: &[L::Collision],
        rotation: Coriolis3,
        mut acceleration_at_end: Vec<[f64; 3]>,
    ) -> Result<(), GridError3<L::Error>> {
        if gate.is_requested() {
            return Err(GridError3::Cancelled);
        }
        let endpoint = EndpointForcing::with_coriolis(&mut acceleration_at_end, rotation)
            .map_err(GridError3::Input)?;
        super::pooled_with_wall_motion(
            &mut self.grid,
            &self.forces,
            endpoint,
            collisions,
            pool,
            gate,
            None,
            |_| Ok(()),
        )?;
        self.forces = acceleration_at_end;
        Ok(())
    }

    /// Close the velocity-dependent force AFTER the actual moving reflection.
    /// Native wall impulse, angular impulse and work remain distinct from the
    /// rotating-flow source. All source, wall and final Mach checks precede
    /// publication. Wall velocities are relative to the same grid axes.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_serial_with_coriolis_and_moving_wall_momentum(
        &mut self,
        cx: &Cx<'_>,
        probe: &WallMomentumProbe3<L, Q>,
        velocities: &[[f64; 3]],
        collisions: &[L::Collision],
        rotation: Coriolis3,
        mut acceleration_at_end: Vec<[f64; 3]>,
    ) -> Result<MovingWallMomentum3, GridError3<L::Error>> {
        let poll = || cx.checkpoint().map_err(|_| GridError3::Cancelled);
        let motion = WallMotion3 { probe, velocities };
        let moving = motion.admit(&self.grid, poll)?;
        let endpoint = EndpointForcing::with_coriolis(&mut acceleration_at_end, rotation)
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

    /// Pooled rotating-flow/moving-wall counterpart with identical cell order.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_pooled_with_coriolis_and_moving_wall_momentum<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        probe: &WallMomentumProbe3<L, Q>,
        velocities: &[[f64; 3]],
        collisions: &[L::Collision],
        rotation: Coriolis3,
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
        let endpoint = EndpointForcing::with_coriolis(&mut acceleration_at_end, rotation)
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
