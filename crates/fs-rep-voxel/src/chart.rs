//! The voxel chart: an fs-geom `Chart` over an occupancy field.
//! Inside/outside from occupancy; distance magnitude from the exact
//! Euclidean DT (both polarities); the declared error model is HONEST —
//! an Enclosure of half a voxel diagonal (resolution error), never
//! "exact".

use crate::dt::{DistanceField, euclidean_dt};
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
    half_diagonal: f64,
}

impl OccupancyChart {
    /// Build from a field (precomputes both distance transforms).
    #[must_use]
    pub fn new(field: OccupancyField) -> Self {
        let to_solid = euclidean_dt(&field);
        let to_void = complement_dt(&field);
        let half_diagonal = 0.5 * fs_math::det::sqrt(3.0) * field.voxel_size;
        OccupancyChart {
            field,
            to_solid,
            to_void,
            half_diagonal,
        }
    }

    /// The underlying field.
    #[must_use]
    pub fn field(&self) -> &OccupancyField {
        &self.field
    }
}

/// DT of the field's complement over a one-voxel-padded active box
/// (distance to the nearest EMPTY voxel; used inside the solid).
fn complement_dt(field: &OccupancyField) -> Option<DistanceField> {
    let mut min = [i32::MAX; 3];
    let mut max = [i32::MIN; 3];
    for (c, _) in field.grid.iter_active() {
        for k in 0..3 {
            min[k] = min[k].min(c[k]);
            max[k] = max[k].max(c[k]);
        }
    }
    if min[0] == i32::MAX {
        return None;
    }
    let mut complement = OccupancyField::new(field.voxel_size, field.origin).ok()?;
    for x in (min[0] - 1)..=(max[0] + 1) {
        for y in (min[1] - 1)..=(max[1] + 1) {
            for z in (min[2] - 1)..=(max[2] + 1) {
                if !field.is_solid([x, y, z]) {
                    complement.set([x, y, z]);
                }
            }
        }
    }
    euclidean_dt(&complement)
}

impl Chart for OccupancyChart {
    fn eval(&self, x: Point3, _cx: &Cx<'_>) -> ChartSample {
        let coord = self.field.voxel_of([x.x, x.y, x.z]);
        let inside = self.field.is_solid(coord);
        let magnitude = if inside {
            self.to_void
                .as_ref()
                .and_then(|dt| dt.distance(coord))
                .unwrap_or(self.field.voxel_size)
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
        ChartSample {
            signed_distance: signed,
            gradient: None,
            lipschitz: Some(1.0),
            error: NumericalCertificate::enclosure(
                signed - self.half_diagonal,
                signed + self.half_diagonal,
            ),
        }
    }

    fn support(&self) -> Aabb {
        let mut min = [f64::INFINITY; 3];
        let mut max = [f64::NEG_INFINITY; 3];
        for (c, _) in self.field.grid.iter_active() {
            let center = self.field.center(c);
            for k in 0..3 {
                min[k] = min[k].min(center[k] - self.field.voxel_size);
                max[k] = max[k].max(center[k] + self.field.voxel_size);
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
        Differentiability::C0
    }
}
