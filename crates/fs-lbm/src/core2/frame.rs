//! Rigid-frame accelerations on the original D2Q9 force/collision path.
//!
//! In frame coordinates, `a = g - a_origin - alpha cross r
//! - omega cross (omega cross r) - 2 omega cross u`. All quantities use
//! lattice units, `dx = dt = 1`; positive spin is counterclockwise.
//! The velocity-dependent Coriolis term uses a two-by-two midpoint solve.
//! Its isolated source step changes direction without creating kinetic energy.
//! This does not make the full spatial/free-surface discretization energy exact.

use super::{Cell, Grid};
use crate::Q;

/// Prescribed rigid-frame motion, expressed in the moving frame.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct RigidFrame2 {
    /// Fixed rotation origin in lattice coordinates (cell centres are integers).
    pub origin: [f64; 2],
    /// Counterclockwise angular velocity, radians per lattice time.
    pub angular_velocity: f64,
    /// Counterclockwise angular acceleration, radians per lattice time squared.
    pub angular_acceleration: f64,
    /// Inertial acceleration of the frame origin, resolved in frame coordinates.
    pub origin_acceleration: [f64; 2],
}

/// Refusal before changing populations or forces.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FrameForcingError2 {
    /// Motion or prescribed acceleration is not finite/representable.
    InvalidMotion,
    /// Public grid storage does not cover its nonempty dimensions.
    InvalidGrid,
    /// A live cell has invalid density, velocity, force or transformed populations.
    InvalidCell {
        /// Row-major cell index.
        cell: usize,
    },
}

impl core::fmt::Display for FrameForcingError2 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidMotion => write!(f, "invalid D2Q9 rigid-frame motion"),
            Self::InvalidGrid => write!(f, "invalid D2Q9 rigid-frame grid storage"),
            Self::InvalidCell { cell } => write!(f, "invalid D2Q9 rigid-frame cell {cell}"),
        }
    }
}
impl std::error::Error for FrameForcingError2 {}

impl RigidFrame2 {
    /// Replace the owned force field for the next native lattice step.
    ///
    /// `acceleration` is the physical body acceleration (e.g. inertial gravity)
    /// already rotated into this frame. Motion and acceleration are held over
    /// one unit lattice step; a schedule should supply their midpoint values.
    /// Pressure, viscosity, wall reactions and free-surface traction remain in
    /// the original flow operator and must not be included here a second time.
    ///
    /// This owner replaces both `grid.g` and `grid.fext`. Each live cell gets
    /// the total force density, including translation, centrifugal, Euler and
    /// midpoint Coriolis terms. Updating Guo's half-source representation keeps
    /// the current physical velocity and stress unchanged: preparation itself
    /// advances no time and applies no impulse. Call the native step afterwards.
    /// Refresh after every step; topology conversion retains the original
    /// free-surface owner's force-aware velocity initialization.
    ///
    /// Two linear passes provide atomic admission without a grid-sized scratch
    /// allocation. The original low-Mach/resolution requirements still apply.
    /// This is a prescribed-frame model, not finite-inertia vessel dynamics.
    ///
    /// # Errors
    /// Invalid motion, grid storage or any nonrepresentable live-cell state.
    pub fn prepare_step(
        self,
        grid: &mut Grid,
        acceleration: [f64; 2],
    ) -> Result<(), FrameForcingError2> {
        if !self
            .origin
            .into_iter()
            .chain(self.origin_acceleration)
            .chain([self.angular_velocity, self.angular_acceleration])
            .chain(acceleration)
            .all(f64::is_finite)
            || !(1.0 + self.angular_velocity * self.angular_velocity).is_finite()
        {
            return Err(FrameForcingError2::InvalidMotion);
        }
        let count = grid
            .nx
            .checked_mul(grid.ny)
            .filter(|&count| count > 0)
            .ok_or(FrameForcingError2::InvalidGrid)?;
        if grid.f.len() != count
            || grid.flags.len() != count
            || grid.fext.len() != count
            || grid.tau.len() != count
            || !grid.g.into_iter().all(f64::is_finite)
        {
            return Err(FrameForcingError2::InvalidGrid);
        }
        for cell in 0..count {
            if matches!(grid.flags[cell], Cell::Fluid | Cell::Interface) {
                self.cell_force(grid, cell, acceleration)?;
            }
        }
        for cell in 0..count {
            if matches!(grid.flags[cell], Cell::Fluid | Cell::Interface) {
                let (force, populations) = self.cell_force(grid, cell, acceleration)?;
                grid.f[cell] = populations;
                grid.fext[cell] = force;
            } else {
                grid.fext[cell] = [0.0; 2];
            }
        }
        grid.g = [0.0; 2];
        Ok(())
    }

    fn cell_force(
        self,
        grid: &Grid,
        cell: usize,
        acceleration: [f64; 2],
    ) -> Result<([f64; 2], [f64; Q]), FrameForcingError2> {
        let invalid = FrameForcingError2::InvalidCell { cell };
        let rho: f64 = grid.f[cell].iter().sum();
        if !rho.is_finite()
            || rho <= 0.0
            || !grid.f[cell]
                .into_iter()
                .chain(grid.fext[cell])
                .all(f64::is_finite)
        {
            return Err(invalid);
        }
        let velocity = grid.moments(cell).u;
        let previous = core::array::from_fn(|d| grid.g[d].mul_add(rho, grid.fext[cell][d]));
        let position = [(cell % grid.nx) as f64, (cell / grid.nx) as f64];
        let r = [position[0] - self.origin[0], position[1] - self.origin[1]];
        let omega = self.angular_velocity;
        let alpha = self.angular_acceleration;
        let base = [
            acceleration[0] - self.origin_acceleration[0] + omega * (omega * r[0]) + alpha * r[1],
            acceleration[1] - self.origin_acceleration[1] + omega * (omega * r[1]) - alpha * r[0],
        ];
        // (I + omega J) u_mid = u_start + a_non_coriolis / 2.
        // No iteration or velocity clipping; the skew source does no midpoint work.
        let b = [velocity[0] + 0.5 * base[0], velocity[1] + 0.5 * base[1]];
        let denominator = 1.0 + omega * omega;
        let midpoint = [
            (b[0] + omega * b[1]) / denominator,
            (b[1] - omega * b[0]) / denominator,
        ];
        let force = [
            rho * (base[0] + 2.0 * omega * midpoint[1]),
            rho * (base[1] - 2.0 * omega * midpoint[0]),
        ];
        if !velocity
            .into_iter()
            .chain(previous)
            .chain(force)
            .all(f64::is_finite)
        {
            return Err(invalid);
        }
        let mut populations = grid.f[cell];
        if previous != force {
            crate::thermal::forcing::reframe2(&mut populations, velocity, previous, force);
        }
        if !populations.into_iter().all(f64::is_finite) {
            return Err(invalid);
        }
        Ok((force, populations))
    }
}
