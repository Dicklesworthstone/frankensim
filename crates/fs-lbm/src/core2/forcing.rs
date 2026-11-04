//! Held-force changes and time-resolved D2Q9 momentum sources.

use super::{Cell, CollisionModel2, Grid, raw_population_moments};
use crate::Q;
use fs_exec::CancelGate;

const FORCE_STEP_CHUNK: usize = 256;

/// A prescribed force-history step was refused without changing the grid.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ForceStepError2 {
    /// Grid dimensions, field dimensions or collision parameters are invalid.
    Invalid(&'static str),
    /// An active cell, imposed force or resulting population is unrepresentable.
    Unphysical {
        /// Row-major cell index.
        cell: usize,
    },
    /// Cancellation was observed before publishing the next fluid state.
    Cancelled,
}

impl core::fmt::Display for ForceStepError2 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "D2Q9 force-history step: {self:?}")
    }
}
impl std::error::Error for ForceStepError2 {}

/// A force-frame proposal was refused before any accepted grid field changed.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ForceRebindError2 {
    /// Invalid grid/field dimensions or a nonfinite acceleration.
    Invalid(&'static str),
    /// A nonfinite force or an unrepresentable active-cell state.
    Unphysical {
        /// Row-major cell index.
        cell: usize,
    },
}

impl core::fmt::Display for ForceRebindError2 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "D2Q9 force rebinding: {self:?}")
    }
}
impl std::error::Error for ForceRebindError2 {}

impl Grid {
    /// Advance one lattice time step with prescribed endpoint acceleration and
    /// external force density, using the native collision and stationary pull.
    ///
    /// The accepted force is `F_n = rho_n*g_n + fext_n`. Collision uses this
    /// force and the resulting populations are interpreted at
    /// `F_(n+1) = rho_(n+1)*gravity + external`. Therefore physical momentum
    /// receives the trapezoidal source impulse `(F_n + F_(n+1))/2`, in addition
    /// to the native momentum flux and stationary-wall reaction. Density in
    /// the endpoint acceleration is the streamed density, not the old density.
    ///
    /// This is the smooth, explicitly prescribed source convention. Switching
    /// a piecewise-held load without advancing time instead uses
    /// [`Self::rebind_forcing`]. Velocity-dependent drag and rotating-frame
    /// loads require their own implicit source solve. Like [`Self::step`], this
    /// method keeps geometry fixed and performs no free-surface mass exchange.
    ///
    /// `post` and `next` are reusable collision and streaming buffers. Refusal
    /// or cancellation may modify them, but publishes neither populations nor
    /// force fields. Cancellation is polled at most every 256 cells of the
    /// numerical traversals. With unchanged forces, accepted populations are
    /// bit-identical to [`Self::step`] for every native collision family.
    pub fn step_with_endpoint_forcing(
        &mut self,
        gravity: [f64; 2],
        external: &[[f64; 2]],
        post: &mut Vec<[f64; Q]>,
        next: &mut Vec<[f64; Q]>,
        gate: &CancelGate,
    ) -> Result<(), ForceStepError2> {
        let checkpoint = || {
            if gate.is_requested() {
                Err(ForceStepError2::Cancelled)
            } else {
                Ok(())
            }
        };
        checkpoint()?;
        let count = self
            .nx
            .checked_mul(self.ny)
            .filter(|&n| n > 0)
            .ok_or(ForceStepError2::Invalid("invalid grid dimensions"))?;
        if self.flags.len() != count
            || self.f.len() != count
            || self.tau.len() != count
            || self.fext.len() != count
            || external.len() != count
        {
            return Err(ForceStepError2::Invalid(
                "populations, flags, relaxation times and force fields must cover the grid",
            ));
        }
        if !gravity.iter().chain(&self.g).all(|g| g.is_finite()) {
            return Err(ForceStepError2::Invalid(
                "old and endpoint acceleration must be finite",
            ));
        }
        if let CollisionModel2::Trt { magic } = self.collision {
            if !magic.is_finite() || magic <= 0.0 {
                return Err(ForceStepError2::Invalid(
                    "TRT magic parameter must be finite and positive",
                ));
            }
        }
        for (cell, next_external) in external.iter().enumerate() {
            if cell % FORCE_STEP_CHUNK == 0 {
                checkpoint()?;
            }
            if !self.fext[cell]
                .iter()
                .chain(next_external)
                .all(|f| f.is_finite())
            {
                return Err(ForceStepError2::Unphysical { cell });
            }
            if !matches!(self.flags[cell], Cell::Fluid | Cell::Interface) {
                continue;
            }
            if !self.tau[cell].is_finite() || self.tau[cell] <= 0.5 {
                return Err(ForceStepError2::Invalid(
                    "active relaxation times must be finite and greater than 0.5",
                ));
            }
            validate_forced_cell(&self.f[cell], self.g, self.fext[cell], cell)?;
        }

        post.clear();
        post.resize(count, [0.0; Q]);
        for (chunk, cells) in post.chunks_mut(FORCE_STEP_CHUNK).enumerate() {
            checkpoint()?;
            self.collide_cells_into(chunk * FORCE_STEP_CHUNK, cells);
        }
        next.clear();
        next.resize(count, [0.0; Q]);
        for (cell, populations) in next.iter_mut().enumerate() {
            if cell % FORCE_STEP_CHUNK == 0 {
                checkpoint()?;
            }
            if !matches!(self.flags[cell], Cell::Fluid | Cell::Interface) {
                *populations = self.f[cell];
                continue;
            }
            let (x, y) = (cell % self.nx, cell / self.nx);
            for (q, population) in populations.iter_mut().enumerate() {
                *population = self.stationary_pull_population(post, cell, q, self.source(x, y, q));
            }
            validate_forced_cell(populations, gravity, external[cell], cell)?;
        }
        checkpoint()?;
        core::mem::swap(&mut self.f, next);
        self.fext.copy_from_slice(external);
        self.g = gravity;
        Ok(())
    }

    /// Replace a held lattice acceleration, preserving physical fluid state.
    ///
    /// At each fluid/interface cell the force is `rho*g + fext`. Replacing its
    /// Guo half-source at the accepted physical velocity preserves density,
    /// momentum and the force-corrected second moment to population roundoff.
    /// This advances no time and supplies no impulse or work; the next native
    /// collision applies the new force. The external force field is retained.
    ///
    /// `scratch` reuses the caller's collision buffer. A refusal may change
    /// scratch, but leaves all accepted grid fields unchanged. Gas and wall
    /// populations are unchanged; identical acceleration preserves all bits.
    pub fn rebind_gravity(
        &mut self,
        gravity: [f64; 2],
        scratch: &mut Vec<[f64; Q]>,
    ) -> Result<(), ForceRebindError2> {
        self.rebind_forcing_inner(gravity, None, scratch)
    }

    /// Replace held acceleration and per-cell external force density together.
    ///
    /// This has the same physical-state-preserving convention as
    /// [`Grid::rebind_gravity`]. The supplied field includes every cell,
    /// including gas and walls that may become wet after a topology change.
    /// Those inactive populations remain unchanged. Use this operation for
    /// piecewise-held loads, switch-on/off or rotation of gravity. Direct force
    /// assignment instead retains the old stored populations and changes their
    /// physical half-source interpretation; it is a different time convention.
    ///
    /// This does not initialize or equilibrate flow. Existing nonequilibrium
    /// stress, cell flags, relaxation times and boundary conditions survive.
    pub fn rebind_forcing(
        &mut self,
        gravity: [f64; 2],
        external: &[[f64; 2]],
        scratch: &mut Vec<[f64; Q]>,
    ) -> Result<(), ForceRebindError2> {
        self.rebind_forcing_inner(gravity, Some(external), scratch)
    }

    fn rebind_forcing_inner(
        &mut self,
        gravity: [f64; 2],
        external: Option<&[[f64; 2]]>,
        scratch: &mut Vec<[f64; Q]>,
    ) -> Result<(), ForceRebindError2> {
        let count = self
            .nx
            .checked_mul(self.ny)
            .filter(|&n| n > 0)
            .ok_or(ForceRebindError2::Invalid("invalid grid dimensions"))?;
        if self.flags.len() != count
            || self.f.len() != count
            || self.fext.len() != count
            || external.is_some_and(|field| field.len() != count)
        {
            return Err(ForceRebindError2::Invalid(
                "populations, flags and force fields must cover the grid",
            ));
        }
        if !gravity.iter().chain(&self.g).all(|g| g.is_finite()) {
            return Err(ForceRebindError2::Invalid(
                "old and new acceleration must be finite",
            ));
        }
        let next_external = external.unwrap_or(&self.fext);
        for (cell, (old, next)) in self.fext.iter().zip(next_external).enumerate() {
            if !old.iter().chain(next).all(|f| f.is_finite()) {
                return Err(ForceRebindError2::Unphysical { cell });
            }
        }
        if gravity == self.g && next_external == self.fext {
            return Ok(());
        }

        scratch.clear();
        scratch.extend_from_slice(&self.f);
        for (cell, populations) in scratch.iter_mut().enumerate() {
            if !matches!(self.flags[cell], Cell::Fluid | Cell::Interface) {
                continue;
            }
            let (rho, momentum) = raw_population_moments(populations);
            if !rho.is_finite()
                || rho <= 0.0
                || !populations.iter().chain(&momentum).all(|v| v.is_finite())
            {
                return Err(ForceRebindError2::Unphysical { cell });
            }
            let previous = core::array::from_fn(|d| self.g[d].mul_add(rho, self.fext[cell][d]));
            let next = core::array::from_fn(|d| gravity[d].mul_add(rho, next_external[cell][d]));
            if !previous.iter().chain(&next).all(|f| f.is_finite()) {
                return Err(ForceRebindError2::Unphysical { cell });
            }
            if previous == next {
                continue;
            }
            let velocity = core::array::from_fn(|d| (momentum[d] + 0.5 * previous[d]) / rho);
            if !velocity.iter().all(|u| u.is_finite()) {
                return Err(ForceRebindError2::Unphysical { cell });
            }
            crate::thermal::forcing::reframe2(populations, velocity, previous, next);
            let (rho_after, momentum_after) = raw_population_moments(populations);
            if !rho_after.is_finite()
                || rho_after <= 0.0
                || !populations.iter().all(|f| f.is_finite())
                || (0..2).any(|d| !((momentum_after[d] + 0.5 * next[d]) / rho_after).is_finite())
            {
                return Err(ForceRebindError2::Unphysical { cell });
            }
        }

        core::mem::swap(&mut self.f, scratch);
        if let Some(field) = external {
            self.fext.copy_from_slice(field);
        }
        self.g = gravity;
        Ok(())
    }
}

fn validate_forced_cell(
    populations: &[f64; Q],
    gravity: [f64; 2],
    external: [f64; 2],
    cell: usize,
) -> Result<(), ForceStepError2> {
    let (rho, momentum) = raw_population_moments(populations);
    if !rho.is_finite() || rho <= 0.0 || !populations.iter().chain(&momentum).all(|f| f.is_finite())
    {
        return Err(ForceStepError2::Unphysical { cell });
    }
    for d in 0..2 {
        let force = gravity[d].mul_add(rho, external[d]);
        if !force.is_finite() || !((momentum[d] + 0.5 * force) / rho).is_finite() {
            return Err(ForceStepError2::Unphysical { cell });
        }
    }
    Ok(())
}
