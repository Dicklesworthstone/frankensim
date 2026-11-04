//! Shared bounded endpoint search for two- and three-coordinate native flow coupling.
//!
//! Physical owners supply dimensional residual allowances and actual represented
//! rates. This layer measures their coupled response, counts no hidden flow
//! solves, and never substitutes a requested rate for native stored momentum.

use fs_la::{LuWorkspace, LuWorkspaceError};
use fs_math::dd::Dd;

#[derive(Clone, Copy)]
pub(super) struct Observation<const N: usize> {
    pub velocity: [f64; N],
    pub residual: [f64; N],
    pub limits: [f64; N],
}
pub(super) enum SearchError<const N: usize> {
    Invalid(&'static str),
    Stalled { residual: [f64; N] },
    LinearSolve(LuWorkspaceError),
}
#[derive(Clone, Copy)]
enum Phase {
    Base,
    Probe(usize),
    Descent,
}

pub(super) struct EndpointSearch<const N: usize> {
    pub requested: [f64; N],
    base: Observation<N>,
    probes: [Observation<N>; N],
    phase: Phase,
    direction: [f64; N],
    fraction: f64,
    bounds: [[f64; 2]; N],
    steps: [f64; N],
    free_axes: [bool; N],
    // The only owners have two or three coordinates; keep their small
    // response matrix on the stack, including the existing two-axis lane.
    matrix: [f64; 9],
}
impl<const N: usize> EndpointSearch<N> {
    pub fn new(velocity: [f64; N], bounds: [[f64; 2]; N], steps: [f64; N]) -> Self {
        Self::with_free_axes(velocity, bounds, steps, [true; N])
    }

    /// Locked coordinates are already stationary in the physical owner. Keep
    /// them exact and omit their fluid probes; the small fixed-size linear
    /// solve uses identity rows/columns for those inactive coordinates.
    pub fn with_free_axes(
        velocity: [f64; N],
        bounds: [[f64; 2]; N],
        steps: [f64; N],
        free_axes: [bool; N],
    ) -> Self {
        assert!(
            N == 2 || N == 3,
            "native endpoint search has two or three coordinates"
        );
        let empty = Observation {
            velocity,
            residual: [0.0; N],
            limits: [1.0; N],
        };
        Self {
            requested: core::array::from_fn(|i| {
                if free_axes[i] {
                    velocity[i].clamp(bounds[i][0], bounds[i][1])
                } else {
                    0.0
                }
            }),
            base: empty,
            probes: [empty; N],
            phase: Phase::Base,
            direction: [0.0; N],
            fraction: 1.0,
            bounds,
            steps,
            free_axes,
            matrix: [0.0; 9],
        }
    }
    fn stalled(&self) -> SearchError<N> {
        SearchError::Stalled {
            residual: self.base.residual,
        }
    }
    fn next_free_axis(&self, first: usize) -> Option<usize> {
        (first..N).find(|&i| self.free_axes[i])
    }
    fn first_probe(&mut self) -> Result<(), SearchError<N>> {
        let axis = self.next_free_axis(0).ok_or_else(|| self.stalled())?;
        self.probe(axis)
    }
    fn probe(&mut self, axis: usize) -> Result<(), SearchError<N>> {
        self.requested = self.base.velocity;
        let [lo, hi] = self.bounds[axis];
        let v = self.requested[axis];
        let (room, sign) = if hi - v >= v - lo {
            (hi - v, 1.0)
        } else {
            (v - lo, -1.0)
        };
        self.requested[axis] = (v + sign * self.steps[axis].min(room)).clamp(lo, hi);
        if self.requested[axis] == v {
            return Err(self.stalled());
        }
        self.phase = Phase::Probe(axis);
        Ok(())
    }
    fn line_point(&mut self) -> Result<(), SearchError<N>> {
        self.requested = core::array::from_fn(|i| {
            if self.free_axes[i] {
                let [lo, hi] = self.bounds[i];
                self.direction[i]
                    .mul_add(self.fraction, self.base.velocity[i])
                    .clamp(lo, hi)
            } else {
                0.0
            }
        });
        if !self.fraction.is_finite()
            || self.fraction <= 0.0
            || self.requested == self.base.velocity
        {
            return Err(self.stalled());
        }
        Ok(())
    }
    pub fn observe(
        &mut self,
        observation: Observation<N>,
        schur: &mut LuWorkspace,
    ) -> Result<(), SearchError<N>> {
        match self.phase {
            Phase::Base => {
                self.base = observation;
                self.first_probe()?;
            }
            Phase::Probe(axis) => {
                self.probes[axis] = observation;
                if let Some(next) = self.next_free_axis(axis + 1) {
                    return self.probe(next);
                }
                // Normalize using actual native probe changes, including their
                // cross-coordinate rounding in the reconstructed direction.
                let h: [f64; N] =
                    core::array::from_fn(|j| self.probes[j].velocity[j] - self.base.velocity[j]);
                if (0..N).any(|j| self.free_axes[j] && (!h[j].is_finite() || h[j] == 0.0)) {
                    return Err(self.stalled());
                }
                for row in 0..N {
                    for col in 0..N {
                        self.matrix[N * row + col] = if self.free_axes[row] && self.free_axes[col] {
                            (self.probes[col].residual[row] - self.base.residual[row])
                                / h[col]
                                / self.base.limits[row]
                        } else if row == col {
                            1.0
                        } else {
                            0.0
                        };
                    }
                }
                let rhs: [f64; N] = core::array::from_fn(|i| {
                    if self.free_axes[i] {
                        -self.base.residual[i] / self.base.limits[i]
                    } else {
                        0.0
                    }
                });
                if !self.matrix.iter().copied().chain(rhs).all(f64::is_finite) {
                    return Err(SearchError::Invalid("unrepresentable coupled response"));
                }
                let mut correction = [0.0; N];
                schur
                    .solve_into(&self.matrix[..N * N], &rhs, &mut correction)
                    .map_err(SearchError::LinearSolve)?;
                self.direction = core::array::from_fn(|i| {
                    if self.free_axes[i] {
                        sum(core::array::from_fn::<_, N, _>(|j| {
                            if self.free_axes[j] {
                                ((self.probes[j].velocity[i] - self.base.velocity[i]) / h[j])
                                    * correction[j]
                            } else {
                                0.0
                            }
                        }))
                    } else {
                        0.0
                    }
                });
                if !self.direction.into_iter().all(f64::is_finite) {
                    return Err(SearchError::Invalid(
                        "unrepresentable coupled velocity correction",
                    ));
                }
                self.fraction = 1.0;
                for i in 0..N {
                    let d = self.direction[i];
                    let [lo, hi] = self.bounds[i];
                    let room = if d > 0.0 {
                        hi - self.base.velocity[i]
                    } else {
                        lo - self.base.velocity[i]
                    };
                    if d != 0.0 {
                        self.fraction = self.fraction.min(room / d);
                    }
                }
                self.phase = Phase::Descent;
                self.line_point()?;
            }
            Phase::Descent => {
                // Keep base physical scaling fixed throughout backtracking.
                let merit = |r: [f64; N]| {
                    maximum(core::array::from_fn::<_, N, _>(|i| {
                        r[i] / self.base.limits[i]
                    }))
                };
                let before = merit(self.base.residual);
                let after = merit(observation.residual);
                if before.is_finite() && after <= (1.0 - 1e-4 * self.fraction) * before {
                    self.base = observation;
                    self.first_probe()?;
                } else {
                    self.fraction *= 0.5;
                    self.line_point()?;
                }
            }
        }
        Ok(())
    }
}

pub(super) fn sum<const N: usize>(values: [f64; N]) -> f64 {
    values
        .into_iter()
        .fold(Dd::ZERO, |s, x| s + Dd::from_f64(x))
        .to_f64()
}
pub(super) fn maximum<const N: usize>(values: [f64; N]) -> f64 {
    values.into_iter().map(f64::abs).fold(0.0, f64::max)
}
