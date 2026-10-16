//! The voxel chart: an fs-geom `Chart` over an occupancy field.
//! Inside/outside comes from occupancy and distance magnitude comes from
//! exact center-to-center Euclidean DTs (both polarities). The resulting
//! representative is piecewise constant over voxel cells: it is not a
//! continuous signed-distance field. Its error model is therefore an
//! estimate-grade full-voxel-diagonal resolution band, never an enclosure.

use crate::VoxelError;
use crate::dt::{CheckedBox, DistanceField, active_bounds, checked_dense_box, euclidean_dt};
use crate::field::OccupancyField;
use fs_evidence::NumericalCertificate;
use fs_exec::Cx;
use fs_geom::{Aabb, Chart, ChartSample, Differentiability, Point3};

/// An occupancy-backed chart with precomputed distance transforms.
pub struct OccupancyChart {
    field: OccupancyField,
    /// Distance to the solid (positive side).
    to_solid: Option<DistanceField>,
    /// Distance to the void (negative side, from the complement within
    /// a one-voxel-padded bounding box).
    to_void: Option<DistanceField>,
    /// Conservative geometric scale for the two center substitutions: the
    /// query point to its voxel center and the target cube to its center.
    /// This is reported only as an estimate because the center metric is not
    /// itself a continuous realization of the union-of-cubes boundary.
    resolution_band: f64,
}

impl OccupancyChart {
    /// Build from a field and precompute both distance transforms.
    ///
    /// `max_voxels` bounds each dense transform and also bounds the
    /// padded complement scan performed during construction.
    ///
    /// # Errors
    /// Returns a structured coordinate, volume, or budget error before
    /// the inadmissible scan or allocation begins. Empty fields are
    /// refused because they cannot produce a finite chart sample.
    pub fn try_new(field: OccupancyField, max_voxels: usize) -> Result<Self, VoxelError> {
        if field.active() == 0 {
            return Err(VoxelError::EmptyOccupancy {
                operation: "occupancy chart construction",
            });
        }
        let (min, max) = active_bounds(&field).ok_or(VoxelError::EmptyOccupancy {
            operation: "occupancy chart construction",
        })?;
        let complement_box =
            checked_dense_box(min, max, 1, max_voxels, "occupancy complement halo")?;
        let to_solid = euclidean_dt(&field, max_voxels)?;
        let to_void = complement_dt(&field, max_voxels, complement_box)?;
        let resolution_band = fs_math::det::sqrt(3.0) * field.voxel_size();
        Ok(OccupancyChart {
            field,
            to_solid,
            to_void,
            resolution_band,
        })
    }

    /// The underlying field.
    #[must_use]
    pub fn field(&self) -> &OccupancyField {
        &self.field
    }
}

/// DT of the field's complement over a one-voxel-padded active box
/// (distance to the nearest EMPTY voxel; used inside the solid).
fn complement_dt(
    field: &OccupancyField,
    max_voxels: usize,
    padded: CheckedBox,
) -> Result<Option<DistanceField>, VoxelError> {
    let mut complement = OccupancyField::new(field.voxel_size(), field.origin())?;
    for x in padded.min[0]..=padded.max[0] {
        for y in padded.min[1]..=padded.max[1] {
            for z in padded.min[2]..=padded.max[2] {
                if !field.is_solid([x, y, z]) {
                    complement.set([x, y, z]);
                }
            }
        }
    }
    euclidean_dt(&complement, max_voxels)
}

impl Chart for OccupancyChart {
    fn eval(&self, x: Point3, _cx: &Cx<'_>) -> ChartSample {
        let Ok(coord) = self.field.voxel_of([x.x, x.y, x.z]) else {
            return ChartSample {
                signed_distance: f64::NAN,
                gradient: None,
                lipschitz: None,
                error: NumericalCertificate::no_claim(),
            };
        };
        let inside = self.field.is_solid(coord);
        let magnitude = if inside {
            self.to_void
                .as_ref()
                .and_then(|dt| dt.distance(coord))
                .unwrap_or(self.field.voxel_size())
        } else {
            // Outside the DT's bounding box, fall back to an exact scan
            // over active-voxel centers (same center-to-center metric).
            self.to_solid
                .as_ref()
                .and_then(|dt| dt.distance(coord))
                .unwrap_or_else(|| {
                    let q = self.field.center(coord);
                    let mut best = f64::INFINITY;
                    for (c, _) in self.field.grid.iter_active() {
                        let p = self.field.center(c);
                        let d = fs_math::det::sqrt(
                            (p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2),
                        );
                        best = best.min(d);
                    }
                    best
                })
        };
        let signed = if inside { -magnitude } else { magnitude };
        let error = if signed.is_finite() && self.resolution_band.is_finite() {
            let lo = signed - self.resolution_band;
            let hi = signed + self.resolution_band;
            if lo.is_finite() && hi.is_finite() {
                NumericalCertificate::estimate(lo, hi)
            } else {
                NumericalCertificate::no_claim()
            }
        } else {
            NumericalCertificate::no_claim()
        };
        ChartSample {
            signed_distance: signed,
            gradient: None,
            // Center selection changes discontinuously at voxel faces, so a
            // finite local Lipschitz claim would be false there.
            lipschitz: None,
            error,
        }
    }

    fn support(&self) -> Aabb {
        let mut min = [f64::INFINITY; 3];
        let mut max = [f64::NEG_INFINITY; 3];
        for (c, _) in self.field.grid.iter_active() {
            let center = self.field.center(c);
            let half_voxel = 0.5 * self.field.voxel_size();
            for k in 0..3 {
                min[k] = min[k].min(center[k] - half_voxel);
                max[k] = max[k].max(center[k] + half_voxel);
            }
        }
        if min[0] > max[0] {
            return Aabb::new(Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 0.0));
        }
        Aabb::new(
            Point3::new(min[0], min[1], min[2]),
            Point3::new(max[0], max[1], max[2]),
        )
    }

    fn name(&self) -> &'static str {
        "voxel-occupancy"
    }

    fn differentiability(&self) -> Differentiability {
        Differentiability::Unknown
    }
}
