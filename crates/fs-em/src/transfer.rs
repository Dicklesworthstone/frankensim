//! Explicit geometry and frame correspondence into existing field consumers.
//!
//! This is point sampling of the computed P1 field, not a conservative mesh
//! projection. In particular it does not certify a different voxel scheme's
//! discrete divergence. Each sample names its source triangle, so a material
//! interface is never crossed or averaged silently. Output vectors can be
//! supplied directly to the existing moving-conductor or magnetothermal LBM
//! owner in the SAME global frame as its velocity and cell-centre positions.

use crate::planar::{MagneticError, PlanarField, gradients, poll};

/// Right-handed orthonormal placement of a planar problem in a 3-D SI frame.
/// Global position = origin + x*e_x + y*e_y + z*(e_x cross e_y).
/// The physical extrusion is 0 <= z <= the problem's planar thickness.
#[derive(Clone, Copy, Debug)]
pub struct PlanarFrame {
    origin_m: [f64; 3],
    axes: [[f64; 3]; 3],
}
impl PlanarFrame {
    /// Admit finite unit perpendicular x/y directions (roundoff tolerance
    /// 1e-12). Their cross product fixes the positive A_z/current direction.
    /// Reflections/scales are not silently absorbed into magnetic components.
    pub fn new(
        origin_m: [f64; 3],
        x_direction: [f64; 3],
        y_direction: [f64; 3],
    ) -> Result<Self, MagneticError> {
        if !origin_m
            .iter()
            .chain(&x_direction)
            .chain(&y_direction)
            .all(|v| v.is_finite())
            || (dot(x_direction, x_direction) - 1.0).abs() > 1e-12
            || (dot(y_direction, y_direction) - 1.0).abs() > 1e-12
            || dot(x_direction, y_direction).abs() > 1e-12
        {
            return Err(MagneticError::Invalid(
                "field frame must be finite and orthonormal",
            ));
        }
        let [a, b, c] = x_direction;
        let [d, e, f] = y_direction;
        let normal = [b * f - c * e, c * d - a * f, a * e - b * d];
        Ok(Self {
            origin_m,
            axes: [x_direction, y_direction, normal],
        })
    }

    pub(crate) fn vector_to_global(&self, local: [f64; 3]) -> [f64; 3] {
        [0, 1, 2].map(|a| {
            local[0].mul_add(
                self.axes[0][a],
                local[1].mul_add(self.axes[1][a], local[2] * self.axes[2][a]),
            )
        })
    }

    pub(crate) fn point_to_global(&self, local: [f64; 3]) -> [f64; 3] {
        let vector = self.vector_to_global(local);
        [0, 1, 2].map(|a| self.origin_m[a] + vector[a])
    }
}

/// Caller-owned correspondence between a consumer point and source triangle.
#[derive(Clone, Copy, Debug)]
pub struct CellSample {
    /// Original mesh triangle index, including an explicit side at interfaces.
    pub triangle: usize,
    /// Actual consumer position in metres in the frame's global coordinates.
    pub position_m: [f64; 3],
}

impl PlanarField<'_> {
    /// Check every geometry correspondence, then express B [T] in the consumer
    /// frame, in sample order. The position tolerance is the caller's admitted
    /// absolute coordinate tolerance; unsupported/outside points refuse.
    /// No source lookup, extrapolation, nearest-cell guess or field averaging
    /// occurs. Cancellation/refusal publishes no partially transferred vector.
    pub fn sample_cells(
        &self,
        frame: &PlanarFrame,
        samples: &[CellSample],
        position_tolerance_m: f64,
        max_samples: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<[f64; 3]>, MagneticError> {
        poll(cancelled)?;
        if samples.len() > max_samples
            || !position_tolerance_m.is_finite()
            || position_tolerance_m < 0.0
        {
            return Err(MagneticError::Invalid(
                "field sample budget or position tolerance is invalid",
            ));
        }
        let mesh = self.problem().mesh();
        let mut result = Vec::with_capacity(samples.len());
        for sample in samples {
            poll(cancelled)?;
            let face = mesh
                .faces()
                .get(sample.triangle)
                .ok_or(MagneticError::Invalid(
                    "sample triangle is outside the source mesh",
                ))?;
            if !sample.position_m.iter().all(|v| v.is_finite()) {
                return Err(MagneticError::Invalid(
                    "field sample position must be finite",
                ));
            }
            let delta = [0, 1, 2].map(|a| sample.position_m[a] - frame.origin_m[a]);
            let local = frame.axes.map(|axis| dot(delta, axis));
            if !local.iter().all(|v| v.is_finite())
                || local[2] < -position_tolerance_m
                || local[2] > self.thickness_m() + position_tolerance_m
            {
                return Err(MagneticError::Invalid(
                    "sample is outside the planar extrusion",
                ));
            }
            let points = face.map(|v| mesh.vertices()[v as usize]);
            let shape = gradients(points)?;
            let offset = [local[0] - points[0][0], local[1] - points[0][1]];
            for (i, gradient) in shape.iter().enumerate() {
                let barycentric = gradient[0].mul_add(offset[0], gradient[1] * offset[1])
                    + if i == 0 { 1.0 } else { 0.0 };
                let margin = position_tolerance_m * (gradient[0].abs() + gradient[1].abs());
                if !barycentric.is_finite()
                    || !margin.is_finite()
                    || barycentric < -margin
                    || barycentric > 1.0 + margin
                {
                    return Err(MagneticError::Invalid(
                        "consumer position is outside its named triangle",
                    ));
                }
            }
            let b = self.magnetic_flux_density_t()[sample.triangle];
            let global = [0, 1, 2].map(|a| b[0].mul_add(frame.axes[0][a], b[1] * frame.axes[1][a]));
            if !global.iter().all(|v| v.is_finite()) {
                return Err(MagneticError::Invalid(
                    "transferred magnetic field is not representable",
                ));
            }
            result.push(global);
        }
        poll(cancelled)?;
        Ok(result)
    }
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0].mul_add(b[0], a[1].mul_add(b[1], a[2] * b[2]))
}
