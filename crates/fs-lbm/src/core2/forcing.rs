//! Changes of held D2Q9 forces without an artificial half-source impulse.

use super::{Cell, Grid, raw_population_moments};
use crate::Q;

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
