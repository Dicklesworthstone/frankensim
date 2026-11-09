//! Bounded characteristic tracing for the existing signed SGS histories.
//!
//! The opt-in scheme solves `dX/ds = -u(X)` in the CURRENT frozen physical
//! velocity field by explicit midpoint RK2. For each authored subinterval h:
//! `X_half = X - h u(X)/2`, `X_next = X - h u(X_half)`.
//! The original positive trilinear stencil transports accepted history only
//! once, at the final departure. Substeps refine trajectories, not the history
//! remap, and therefore do not add repeated interpolation diffusion.
//!
//! Characteristic integration and history interpolation are separate errors;
//! see Diamantakis, ECMWF Seminar 2013, §§2.1--2.2:
//! <https://www.ecmwf.int/sites/default/files/elibrary/2014/9054-semi-lagrangian-technique-atmospheric-modelling-current-status-and-future-challenges.pdf>.
//! Unlike the temporally extrapolated winds discussed there, this owner has
//! only the current field. RK2 gives second-order trajectory convergence for
//! a smooth FROZEN field; it does not upgrade unsteady velocity sampling,
//! frozen-source relaxation, trilinear history accuracy, or the whole solver
//! to second order. No large-Courant or DNS qualification is supplied.

use super::*;

/// One replayable policy for every history in a trajectory. The original
/// Euler policy preserves its arithmetic and cancellation behavior exactly.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum LagrangianAdvection {
    #[default]
    Euler,
    /// Positive authored substep count, at most 255. Each substep uses two
    /// velocity evaluations, with the initial cell value reused. Every stage
    /// and the accumulated per-axis travel must meet the existing departure
    /// bound; substepping never bypasses that admission or a physical wall.
    Midpoint { substeps: u8 },
}

impl LagrangianAdvection {
    pub fn validate(self) -> Result<(), DynamicSgsError> {
        if self == (Self::Midpoint { substeps: 0 }) {
            return Err(DynamicSgsError::Invalid(
                "midpoint characteristic tracing needs at least one substep",
            ));
        }
        Ok(())
    }
}

impl LagrangianDynamicSgs {
    /// Select characteristic tracing before the first accepted sample.
    /// This fixed-size policy adds no numerical field storage. It also governs
    /// the advected contractions used by an adaptive momentum-memory timescale.
    pub fn with_advection(mut self, policy: LagrangianAdvection) -> Result<Self, DynamicSgsError> {
        policy.validate()?;
        if self.accepted_time.is_some() {
            return Err(DynamicSgsError::Invalid(
                "select Lagrangian advection before the first accepted sample",
            ));
        }
        self.advection = policy;
        Ok(self)
    }

    pub const fn advection(&self) -> LagrangianAdvection {
        self.advection
    }
}

impl HistoryUpdate {
    pub(super) fn with_advection(
        mut self,
        policy: LagrangianAdvection,
        walls: [[f64; 3]; 2],
    ) -> Self {
        self.advection = policy;
        self.wall_velocity = walls;
        self
    }

    pub(super) fn pair_with_field(
        &self,
        cell: usize,
        velocity: &[[f64; 3]],
        history: &[[f64; 2]],
        contraction: &[[f64; 2]],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[f64; 2], DynamicSgsError> {
        if self.advection == LagrangianAdvection::Euler {
            return self.pair(cell, velocity[cell], history, contraction);
        }
        self.pair_using(cell, contraction, || {
            self.advected_with_field(cell, velocity, history, cancelled)
        })
    }

    pub(super) fn advected_with_field(
        &self,
        cell: usize,
        velocity: &[[f64; 3]],
        history: &[[f64; 2]],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[f64; 2], DynamicSgsError> {
        if self.advection == LagrangianAdvection::Euler || self.elapsed == 0.0 {
            return self.advected(cell, velocity[cell], history);
        }
        let departure = self.midpoint_departure(cell, velocity, cancelled)?;
        sample_history(self.geometry, self.shape, cell, departure, history)
    }

    fn midpoint_departure(
        &self,
        cell: usize,
        velocity: &[[f64; 3]],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[f64; 3], DynamicSgsError> {
        self.advection.validate()?;
        let LagrangianAdvection::Midpoint { substeps } = self.advection else {
            return Err(DynamicSgsError::Invalid("midpoint policy required"));
        };
        let interval = self.elapsed / f64::from(substeps);
        if self.elapsed > 0.0 && interval == 0.0 {
            return Err(DynamicSgsError::Invalid(
                "unrepresentable midpoint characteristic subinterval",
            ));
        }
        let time_per_cell = self.spacing.map(|h| interval / h);
        let mut departure = [0.0; 3];
        let mut travel = [0.0; 3];
        for step in 0..substeps {
            // The authored work ceiling cannot hide cancellation inside a
            // long trajectory: poll at most eight substeps apart per cell.
            if step % 8 == 0 && cancelled() {
                return Err(DynamicSgsError::Cancelled);
            }
            let start = if step == 0 {
                velocity[cell]
            } else {
                self.velocity_at(cell, departure, velocity)?
            };
            let half = core::array::from_fn::<_, 3, _>(|a| 0.5 * (start[a] * time_per_cell[a]));
            let midpoint = core::array::from_fn(|a| departure[a] + half[a]);
            self.admit_path(
                cell,
                midpoint,
                core::array::from_fn(|a| travel[a] + half[a].abs()),
            )?;
            let middle = self.velocity_at(cell, midpoint, velocity)?;
            let increment = core::array::from_fn::<_, 3, _>(|a| middle[a] * time_per_cell[a]);
            for a in 0..3 {
                departure[a] += increment[a];
                travel[a] += increment[a].abs();
            }
            self.admit_path(cell, departure, travel)?;
        }
        Ok(departure)
    }

    fn admit_path(
        &self,
        cell: usize,
        departure: [f64; 3],
        travel: [f64; 3],
    ) -> Result<(), DynamicSgsError> {
        if (0..3).any(|a| {
            !departure[a].is_finite()
                || !travel[a].is_finite()
                || departure[a].abs() > self.maximum_departure_cells
                || travel[a] > self.maximum_departure_cells
        }) {
            return Err(DynamicSgsError::Cell {
                cell,
                reason: "midpoint characteristic travel exceeds declared cell limit",
            });
        }
        if let FilterGeometry::Channel(normal) = self.geometry {
            let position = coordinates(self.shape, cell)[normal] as f64 - departure[normal];
            if position < -0.5 || position > self.shape[normal] as f64 - 0.5 {
                return Err(DynamicSgsError::Cell {
                    cell,
                    reason: "midpoint characteristic crosses impermeable channel wall",
                });
            }
        }
        Ok(())
    }

    fn velocity_at(
        &self,
        cell: usize,
        departure: [f64; 3],
        velocity: &[[f64; 3]],
    ) -> Result<[f64; 3], DynamicSgsError> {
        let mut stencil = [(0usize, 0.0); 8];
        let mut next = 0;
        each_history_sample(self.geometry, self.shape, cell, departure, |j, weight| {
            stencil[next] = (j, weight);
            next += 1;
        })?;
        let mut value = [0.0; 3];
        for a in 0..3 {
            // Component normalization avoids overflowing differences between
            // finite velocities. Anchoring preserves a uniform field exactly.
            let scale = stencil.iter().fold(velocity[cell][a].abs(), |s, &(j, _)| {
                s.max(velocity[j][a].abs())
            });
            if scale == 0.0 {
                continue;
            }
            let reference = velocity[cell][a] / scale;
            let mut normalized = reference;
            for &(j, weight) in &stencil {
                normalized += weight * (velocity[j][a] / scale - reference);
            }
            value[a] = normalized * scale;
        }
        if let FilterGeometry::Channel(normal) = self.geometry {
            let position = coordinates(self.shape, cell)[normal] as f64 - departure[normal];
            let last = (self.shape[normal] - 1) as f64;
            let side_and_weight = if position < 0.0 {
                Some((0, -2.0 * position))
            } else if position > last {
                Some((1, 2.0 * (position - last)))
            } else {
                None
            };
            if let Some((side, weight)) = side_and_weight {
                // History retains its original nearest-center extension;
                // velocity instead resolves the actual physical half-cell
                // wall, including its authored tangential translation.
                for a in 0..3 {
                    let wall = self.wall_velocity[side][a];
                    if wall != value[a] {
                        value[a] = (1.0 - weight) * value[a] + weight * wall;
                    }
                }
            }
        }
        if value.iter().any(|v| !v.is_finite()) {
            return Err(DynamicSgsError::Cell {
                cell,
                reason: "unrepresentable characteristic velocity interpolation",
            });
        }
        Ok(value)
    }
}

#[cfg(test)]
#[path = "advection/tests.rs"]
mod tests;
