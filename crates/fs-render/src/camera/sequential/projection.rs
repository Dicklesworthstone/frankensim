//! Local inverse of one native spectral lens/pupil realization.
//!
//! The original Snell/coating/ghost traversal defines the image. A bounded
//! damped solve changes only the physical sensor coordinate. Its finite
//! differences are solver directions, not certified optical derivatives.
//! Every accepted projection is checked against an actual native outgoing ray.

use super::{CameraError, LensSample, PhysicalCamera};
use crate::charts::Ray;
use fs_exec::Cx;
use fs_geom::{Point3, Vec3};
use fs_math::{dd::Dd, det};

/// Explicit local branch seed, physical residual tolerance and native work cap.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SequentialProjectionOptions {
    /// Initial sensor tangent coordinates, in the original camera convention.
    pub initial_tangent_xy: [f64; 2],
    /// Maximum transverse miss distance at the target's axial plane [m].
    pub position_tolerance_m: f64,
    /// Relative sensor-coordinate difference used for local solve directions.
    pub difference_step: f64,
    /// Maximum Euclidean change of sensor tangent coordinates per iteration.
    pub maximum_tangent_step: f64,
    /// Counts all original ray traces, including clipped probes and backtracks.
    pub maximum_trace_evaluations: u32,
}

impl Default for SequentialProjectionOptions {
    fn default() -> Self {
        Self {
            initial_tangent_xy: [0.0; 2],
            position_tolerance_m: 1e-8,
            difference_step: 1e-5,
            maximum_tangent_step: 0.1,
            maximum_trace_evaluations: 128,
        }
    }
}

/// A converged conditional lens projection; it does not establish visibility.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SequentialProjection {
    tangent_xy: [f64; 2],
    ndc_xy: [f64; 2],
    depth_m: f64,
    residual_m: f64,
    trace_evaluations: u32,
}

impl SequentialProjection {
    pub const fn tangent_xy(self) -> [f64; 2] {
        self.tangent_xy
    }
    /// +x right and +y up. Off-screen values remain valid geometric results.
    pub const fn ndc_xy(self) -> [f64; 2] {
        self.ndc_xy
    }
    pub const fn depth_m(self) -> f64 {
        self.depth_m
    }
    pub const fn residual_m(self) -> f64 {
        self.residual_m
    }
    pub const fn trace_evaluations(self) -> u32 {
        self.trace_evaluations
    }
}

/// Local inverse failure, kept distinct from a projected or invisible point.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SequentialProjectionError {
    InvalidOptions,
    /// A finite scene point must lie beyond the complete front lens cap.
    TargetOutsideModel,
    /// The caller's selected branch has no admitted ray at its initial seed.
    SeedClipped,
    /// Neither side of a bounded local stencil supplies an admitted ray.
    DerivativeUnavailable,
    SingularJacobian,
    NoDescent,
    NumericalRange,
    WorkLimit {
        trace_evaluations: u32,
    },
    Camera(CameraError),
}

impl From<CameraError> for SequentialProjectionError {
    fn from(error: CameraError) -> Self {
        Self::Camera(error)
    }
}

impl core::fmt::Display for SequentialProjectionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "sequential projection refused: {self:?}")
    }
}

impl core::error::Error for SequentialProjectionError {}

type Result<T> = core::result::Result<T, SequentialProjectionError>;

impl PhysicalCamera {
    /// Project a scene point through one original wavelength and pupil sample.
    ///
    /// The same lens sample retains the native direct/ghost path selection and
    /// conditional pupil coordinate at every trial. Source bands, coatings,
    /// absorption, stops and TIR use the original steady ray generator; this
    /// query never requests group timing. The local seed selects a root basin.
    /// A solution is accepted only when an actual outgoing ray approaches the
    /// target within the requested physical tolerance in its axial plane.
    ///
    /// This is a local conditional image, not a unique projection of a blurred
    /// point, a search for all ghost images, or a visibility test. A clipped
    /// seed, unresolved inverse, exhausted work or cancellation returns an
    /// error, never a fabricated image coordinate.
    #[allow(clippy::too_many_arguments)]
    pub fn project_sequential_point(
        &self,
        cx: &Cx<'_>,
        world_point: Point3,
        aspect_ratio: f64,
        wavelength_nm: f64,
        lens_sample: LensSample,
        options: SequentialProjectionOptions,
    ) -> Result<SequentialProjection> {
        cx.checkpoint().map_err(|_| CameraError::Cancelled)?;
        if !options.initial_tangent_xy.into_iter().all(f64::is_finite)
            || ![
                options.position_tolerance_m,
                options.difference_step,
                options.maximum_tangent_step,
            ]
            .into_iter()
            .all(|v| v.is_finite() && v > 0.0)
            || options.maximum_trace_evaluations == 0
            || !aspect_ratio.is_finite()
            || aspect_ratio <= 0.0
        {
            return Err(SequentialProjectionError::InvalidOptions);
        }
        if ![world_point.x, world_point.y, world_point.z]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(CameraError::NonFiniteInput.into());
        }
        let lens = self
            .sequential_lens()
            .ok_or(CameraError::InvalidSequentialLens)?;
        let depth_m = dot(difference(world_point, self.eye()), self.forward()).to_f64();
        if !depth_m.is_finite()
            || depth_m <= lens.surfaces().last().expect("admitted lens").bounds()[1]
        {
            return Err(SequentialProjectionError::TargetOutsideModel);
        }
        let mut trace = Trace {
            camera: self,
            cx,
            target: world_point,
            wavelength_nm,
            lens_sample,
            count: 0,
            maximum: options.maximum_trace_evaluations,
        };
        let mut tangent = options.initial_tangent_xy;
        let mut residual = trace
            .evaluate(tangent)?
            .ok_or(SequentialProjectionError::SeedClipped)?;
        loop {
            let norm = residual_norm(residual)?;
            if norm <= options.position_tolerance_m {
                let half = self.raster_vertical_half_tan();
                let ndc_xy = [tangent[0] / aspect_ratio / half, tangent[1] / half];
                if !ndc_xy.into_iter().all(f64::is_finite) {
                    return Err(SequentialProjectionError::NumericalRange);
                }
                return Ok(SequentialProjection {
                    tangent_xy: tangent,
                    ndc_xy,
                    depth_m,
                    residual_m: norm,
                    trace_evaluations: trace.count,
                });
            }
            let x = trace.derivative(tangent, residual, 0, options.difference_step)?;
            let y = trace.derivative(tangent, residual, 1, options.difference_step)?;
            let mut step = solve_step([x, y], residual)?;
            let length = residual_norm(step)?;
            if length > options.maximum_tangent_step {
                let scale = options.maximum_tangent_step / length;
                step = step.map(|value| value * scale);
            }
            let mut accepted = None;
            let mut scale = 1.0;
            for _ in 0..12 {
                let candidate = [tangent[0] + scale * step[0], tangent[1] + scale * step[1]];
                if candidate == tangent {
                    break;
                }
                if candidate.into_iter().all(f64::is_finite) {
                    if let Some(value) = trace.evaluate(candidate)? {
                        if residual_norm(value)? < norm {
                            accepted = Some((candidate, value));
                            break;
                        }
                    }
                }
                scale *= 0.5;
            }
            let Some((next, value)) = accepted else {
                return Err(SequentialProjectionError::NoDescent);
            };
            tangent = next;
            residual = value;
        }
    }
}

struct Trace<'a, 'cx> {
    camera: &'a PhysicalCamera,
    cx: &'a Cx<'cx>,
    target: Point3,
    wavelength_nm: f64,
    lens_sample: LensSample,
    count: u32,
    maximum: u32,
}

impl Trace<'_, '_> {
    fn evaluate(&mut self, tangent: [f64; 2]) -> Result<Option<[f64; 2]>> {
        if self.count == self.maximum {
            return Err(SequentialProjectionError::WorkLimit {
                trace_evaluations: self.count,
            });
        }
        self.count += 1;
        let [sample] = self
            .camera
            .generate_sequential_spectral_rays_from_tangent_offsets(
                self.cx,
                tangent[0],
                tangent[1],
                self.lens_sample,
                &[self.wavelength_nm],
            )?;
        sample.map(|sample| self.residual(sample.ray)).transpose()
    }

    fn residual(&self, ray: Ray) -> Result<[f64; 2]> {
        let denominator = ray.dir.dot(self.camera.forward());
        if !denominator.is_finite() || denominator <= 0.0 {
            return Err(SequentialProjectionError::NumericalRange);
        }
        let delta = difference(ray.origin, self.target);
        let distance = -dot(delta, self.camera.forward()) / Dd::from_f64(denominator);
        if !distance.hi.is_finite() || !distance.lo.is_finite() || distance.lt(Dd::ZERO) {
            return Err(SequentialProjectionError::TargetOutsideModel);
        }
        let residual = [self.camera.right(), self.camera.up()]
            .map(|axis| (dot(delta, axis) + distance * Dd::from_f64(ray.dir.dot(axis))).to_f64());
        if !residual.into_iter().all(f64::is_finite) {
            return Err(SequentialProjectionError::NumericalRange);
        }
        Ok(residual)
    }

    fn derivative(
        &mut self,
        tangent: [f64; 2],
        residual: [f64; 2],
        axis: usize,
        relative_step: f64,
    ) -> Result<[f64; 2]> {
        let mut h = relative_step * tangent[axis].abs().max(1.0);
        for _ in 0..4 {
            let mut plus = tangent;
            let mut minus = tangent;
            plus[axis] += h;
            minus[axis] -= h;
            if !plus[axis].is_finite()
                || !minus[axis].is_finite()
                || plus[axis] == tangent[axis]
                || minus[axis] == tangent[axis]
            {
                return Err(SequentialProjectionError::NumericalRange);
            }
            let positive = self.evaluate(plus)?;
            let negative = self.evaluate(minus)?;
            let (a, b, distance) = match (positive, negative) {
                (Some(a), Some(b)) => (a, b, plus[axis] - minus[axis]),
                (Some(a), None) => (a, residual, plus[axis] - tangent[axis]),
                (None, Some(b)) => (residual, b, tangent[axis] - minus[axis]),
                (None, None) => {
                    h *= 0.5;
                    continue;
                }
            };
            let result = [
                ((Dd::from_f64(a[0]) - Dd::from_f64(b[0])) / Dd::from_f64(distance)).to_f64(),
                ((Dd::from_f64(a[1]) - Dd::from_f64(b[1])) / Dd::from_f64(distance)).to_f64(),
            ];
            if !result.into_iter().all(f64::is_finite) {
                return Err(SequentialProjectionError::NumericalRange);
            }
            return Ok(result);
        }
        Err(SequentialProjectionError::DerivativeUnavailable)
    }
}

// Columns of the native local Jacobian. Scaling precedes the DD determinant
// so a weak but resolved determinant is not confused with arithmetic overflow.
fn solve_step(columns: [[f64; 2]; 2], residual: [f64; 2]) -> Result<[f64; 2]> {
    let scale = columns
        .into_iter()
        .flatten()
        .map(f64::abs)
        .fold(0.0, f64::max);
    if !scale.is_finite() || scale == 0.0 {
        return Err(SequentialProjectionError::SingularJacobian);
    }
    let [[a, c], [b, d]] =
        columns.map(|column| column.map(|v| Dd::from_f64(v) / Dd::from_f64(scale)));
    let [r, s] = residual.map(|v| -Dd::from_f64(v) / Dd::from_f64(scale));
    let determinant = a * d - b * c;
    if determinant == Dd::ZERO {
        return Err(SequentialProjectionError::SingularJacobian);
    }
    let result = [
        ((d * r - b * s) / determinant).to_f64(),
        ((a * s - c * r) / determinant).to_f64(),
    ];
    if !result.into_iter().all(f64::is_finite) {
        return Err(SequentialProjectionError::NumericalRange);
    }
    Ok(result)
}

fn difference(a: Point3, b: Point3) -> [Dd; 3] {
    [
        Dd::from_f64(a.x) - Dd::from_f64(b.x),
        Dd::from_f64(a.y) - Dd::from_f64(b.y),
        Dd::from_f64(a.z) - Dd::from_f64(b.z),
    ]
}

fn dot(v: [Dd; 3], axis: Vec3) -> Dd {
    v[0] * Dd::from_f64(axis.x) + v[1] * Dd::from_f64(axis.y) + v[2] * Dd::from_f64(axis.z)
}

fn residual_norm(value: [f64; 2]) -> Result<f64> {
    let norm = det::hypot(value[0], value[1]);
    if norm.is_finite() {
        Ok(norm)
    } else {
        Err(SequentialProjectionError::NumericalRange)
    }
}
