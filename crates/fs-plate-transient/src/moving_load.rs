//! Prescribed finite-footprint loading through the original plate dynamics.
//!
//! A fixed quadrature footprint translates along a piecewise-linear path. The
//! three total physical loads at each time knot are the design parameters;
//! geometry, knot times and quadrature weights are held fixed by the adjoint.
//! Location, load transfer and observation share PlatePatchStencil. Power uses
//! material velocity, not the derivative of a probe moving over a displaced
//! surface. This is a prescribed applied load, not a contact law.

use crate::{PlateLoad, PlateModel};
use fs_plate::{
    PlateMesh,
    loading::{
        PlateAppliedLoad, PlateLoadError, PlatePatchBudget, PlatePatchPoint, PlatePatchStencil,
        PlatePointSample,
    },
};

/// One physical-time knot of a translating patch and its total applied load.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlateLoadKeyframe {
    /// Physical time in seconds, nonnegative and strictly increasing.
    pub time_s: f64,
    /// Translation of the footprint in the undeformed xy plane, in metres.
    pub center_m: [f64; 2],
    /// Total transverse force and Cartesian bending couples [Fz, Mx, My],
    /// in N, N m, N m, distributed in proportion to quadrature area.
    pub load: [f64; 3],
}

/// Explicit bounds for admission and every subsequent location query.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MovingPlateLoadBudget {
    /// Maximum number of time/path/load knots inspected at admission.
    pub max_keyframes: usize,
    /// Footprint size and per-point spatial-location bounds.
    pub patch: PlatePatchBudget,
}

/// A moving load refuses before publishing a partial force or derivative.
#[derive(Debug, Clone, PartialEq)]
pub enum MovingPlateLoadError {
    /// Caller-supplied storage/work bounds cannot cover the supplied data.
    Budget(&'static str),
    /// Invalid path, time or load data.
    InvalidInput(&'static str),
    /// Evaluation would extrapolate the explicitly supplied time horizon.
    OutsideTime {
        /// Requested time in seconds.
        time_s: f64,
        /// Admitted first and last knot times in seconds.
        horizon_s: [f64; 2],
    },
    /// The original spatial transfer owner's refusal.
    Transfer(PlateLoadError),
}

impl core::fmt::Display for MovingPlateLoadError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-PLATE-MOVING-LOAD: {self:?}")
    }
}
impl std::error::Error for MovingPlateLoadError {}
impl From<PlateLoadError> for MovingPlateLoadError {
    fn from(value: PlateLoadError) -> Self {
        Self::Transfer(value)
    }
}

/// Immutable moving-patch load, usable directly by the original PlateDynamics.
///
/// Footprint point entries are offsets from each keyframe's center_m; their
/// area_m2 entries are physical quadrature measures. Forces/couples interpolate
/// linearly between knots, independently of position. Evaluation is stateless:
/// rejected integrator trials, replay and out-of-order queries cannot advance a
/// load cursor or debit work twice. A changed mesh/support map requires rebinding.
pub struct MovingPlateLoad<'a> {
    mesh: &'a PlateMesh,
    model: &'a PlateModel,
    footprint: &'a [PlatePatchPoint],
    keyframes: &'a [PlateLoadKeyframe],
    budget: MovingPlateLoadBudget,
    parameters: usize,
}

impl<'a> MovingPlateLoad<'a> {
    /// Admit finite data, bounded dimensions and all knot locations.
    ///
    /// Paths through nonconvex meshes may leave the domain between valid knots;
    /// every evaluation locates the actual translated footprint and refuses if
    /// any quadrature point is outside. No clipping, extrapolation or silent
    /// reduction of the footprint is performed. The caller owns quadrature
    /// accuracy and path resolution; this does not certify the swept area.
    ///
    /// # Errors
    /// Invalid knots or footprint, exhausted bounds, or spatial-location refusal.
    pub fn new(
        mesh: &'a PlateMesh,
        model: &'a PlateModel,
        footprint: &'a [PlatePatchPoint],
        keyframes: &'a [PlateLoadKeyframe],
        budget: MovingPlateLoadBudget,
    ) -> Result<Self, MovingPlateLoadError> {
        if keyframes.len() > budget.max_keyframes || footprint.len() > budget.patch.max_points {
            return Err(MovingPlateLoadError::Budget(
                "keyframes or footprint points",
            ));
        }
        let parameters = keyframes
            .len()
            .checked_mul(3)
            .ok_or(MovingPlateLoadError::Budget(
                "load parameter count overflow",
            ))?;
        if keyframes.len() < 2 || footprint.is_empty() {
            return Err(MovingPlateLoadError::InvalidInput(
                "at least two time knots and a nonempty footprint required",
            ));
        }
        for (i, knot) in keyframes.iter().enumerate() {
            if !knot.time_s.is_finite()
                || knot.time_s < 0.0
                || knot
                    .center_m
                    .iter()
                    .chain(&knot.load)
                    .any(|v| !v.is_finite())
                || (i > 0 && knot.time_s <= keyframes[i - 1].time_s)
            {
                return Err(MovingPlateLoadError::InvalidInput(
                    "finite positions/loads and nonnegative strictly increasing times required",
                ));
            }
        }
        let load = Self {
            mesh,
            model,
            footprint,
            keyframes,
            budget,
            parameters,
        };
        for knot in keyframes {
            load.locate(knot.center_m)?;
        }
        Ok(load)
    }

    /// Closed physical-time horizon. Evaluations outside it refuse.
    #[must_use]
    pub fn horizon_s(&self) -> [f64; 2] {
        [
            self.keyframes[0].time_s,
            self.keyframes[self.keyframes.len() - 1].time_s,
        ]
    }

    /// Resolve the same finite footprint and physical load used by the solver.
    ///
    /// # Errors
    /// Refuses times outside the horizon or a footprint outside the spatial domain.
    pub fn stencil(
        &self,
        time_s: f64,
    ) -> Result<(PlatePatchStencil, [f64; 3]), MovingPlateLoadError> {
        let (left, right, fraction) = self.interval(time_s)?;
        let a = &self.keyframes[left];
        let b = &self.keyframes[right];
        let center = std::array::from_fn(|i| interpolate(a.center_m[i], b.center_m[i], fraction));
        let load = std::array::from_fn(|i| interpolate(a.load[i], b.load[i], fraction));
        Ok((self.locate(center)?, load))
    }

    /// Overwrite the reduced load, returning global resultants and the portions
    /// applied on fixed DOFs. Those portions are not complete support reactions.
    /// The output remains unchanged on any refusal.
    ///
    /// # Errors
    /// Invalid time, location, dimensions or unrepresentable load accumulation.
    pub fn applied_load(
        &self,
        time_s: f64,
        output: &mut [f64],
    ) -> Result<PlateAppliedLoad, MovingPlateLoadError> {
        self.dimension(output, self.model.free)?;
        let (stencil, load) = self.stencil(time_s)?;
        let mut candidate = vec![0.0; self.model.free];
        let applied = stencil.add_load(load, &mut candidate)?;
        output.copy_from_slice(&candidate);
        Ok(applied)
    }

    /// Area-averaged displacement/slopes, material velocity and instantaneous
    /// supplied power in watts. Observing does not accumulate or debit energy.
    ///
    /// # Errors
    /// Invalid time, location, state dimensions or nonfinite observations.
    pub fn sample(
        &self,
        time_s: f64,
        displacement: &[f64],
        velocity: &[f64],
    ) -> Result<(PlatePointSample, f64), MovingPlateLoadError> {
        let (stencil, load) = self.stencil(time_s)?;
        Ok((
            stencil.sample(displacement, velocity)?,
            stencil.power(load, velocity)?,
        ))
    }

    fn locate(&self, center: [f64; 2]) -> Result<PlatePatchStencil, MovingPlateLoadError> {
        let points: Vec<_> = self
            .footprint
            .iter()
            .map(|p| PlatePatchPoint {
                point: [center[0] + p.point[0], center[1] + p.point[1]],
                area_m2: p.area_m2,
            })
            .collect();
        Ok(PlatePatchStencil::locate(
            self.mesh,
            self.model,
            &points,
            self.budget.patch,
        )?)
    }

    fn interval(&self, time_s: f64) -> Result<(usize, usize, f64), MovingPlateLoadError> {
        let [first, last] = self.horizon_s();
        if !time_s.is_finite() || time_s < first || time_s > last {
            return Err(MovingPlateLoadError::OutsideTime {
                time_s,
                horizon_s: [first, last],
            });
        }
        let right = self.keyframes.partition_point(|knot| knot.time_s < time_s);
        if self.keyframes[right].time_s == time_s {
            return Ok((right, right, 0.0));
        }
        let left = right - 1;
        let fraction = (time_s - self.keyframes[left].time_s)
            / (self.keyframes[right].time_s - self.keyframes[left].time_s);
        Ok((left, right, fraction))
    }

    fn dimension(&self, values: &[f64], expected: usize) -> Result<(), MovingPlateLoadError> {
        if values.len() != expected {
            return Err(PlateLoadError::Dimension {
                expected,
                actual: values.len(),
            }
            .into());
        }
        Ok(())
    }

    fn pullback(
        &self,
        time_s: f64,
        seed: &[f64],
        output: &mut [f64],
    ) -> Result<(), MovingPlateLoadError> {
        self.dimension(output, self.parameters)?;
        let (left, right, fraction) = self.interval(time_s)?;
        let (stencil, _) = self.stencil(time_s)?;
        let physical = stencil.load_vjp(seed)?;
        let mut candidate = vec![0.0; self.parameters];
        for i in 0..3 {
            candidate[3 * left + i] = (1.0 - fraction) * physical[i];
            candidate[3 * right + i] += fraction * physical[i];
        }
        output.copy_from_slice(&candidate);
        Ok(())
    }
}

impl PlateLoad for MovingPlateLoad<'_> {
    fn parameter_count(&self) -> usize {
        self.parameters
    }

    fn forcing(&self, time: f64, output: &mut [f64]) -> Result<(), String> {
        self.applied_load(time, output)
            .map(|_| ())
            .map_err(|error| error.to_string())
    }

    fn forcing_vjp(&self, time: f64, seed: &[f64], output: &mut [f64]) -> Result<(), String> {
        self.pullback(time, seed, output)
            .map_err(|error| error.to_string())
    }
}

// Avoid subtracting opposite-sign extreme endpoints, and preserve exact knots.
fn interpolate(a: f64, b: f64, fraction: f64) -> f64 {
    if fraction == 0.0 {
        a
    } else if fraction == 1.0 {
        b
    } else if a.is_sign_negative() == b.is_sign_negative() {
        fraction.mul_add(b - a, a)
    } else {
        fraction.mul_add(b, (1.0 - fraction) * a)
    }
}
