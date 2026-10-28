//! Moving geometry on the original partially saturated D2Q9 operator.
//!
//! The disk half-drifts, exchanges a coupled impulse at that configuration,
//! then half-drifts with the accepted velocity. Each free drift is parallel
//! to its own linear momentum, so it cannot create orbital angular momentum.
//! Cell fractions are analytic circle/square intersection areas; populations
//! remain on the original grid and are never invented, removed or remapped.

use super::{Cell, Grid, PartialSaturationBody2, PartialSaturationBodyExchange2};
use crate::Q;
use fs_math::det;

/// One moving disk and reusable geometry/collision scratch, in lattice units.
///
/// `body.mass` and `body.inertia` describe the **additional solid degree of
/// freedom** in the partial-saturation model. Fluid populations, including
/// those inside the disk, retain their mass and momentum. These parameters
/// are not silently interpreted as physical solid density times disk volume;
/// a displaced/interior-fluid inertia closure is not supplied here.
///
/// One disk, fixed radius and inertia, force-free BGK fluid, and lattice
/// `dx = dt = 1` are admitted. Its swept bounding boxes must contain only
/// fluid cells and retain a one-cell domain halo. Collision/contact with
/// permanent walls, periodic seam crossing, multiple disks and free surfaces
/// therefore refuse rather than inventing a remapping or contact law.
/// Outer fluid boundaries retain the grid's original streaming behavior.
#[derive(Debug, Clone)]
pub struct MovingDisk2 {
    body: PartialSaturationBody2,
    radius: f64,
    angle: f64,
    dims: [usize; 2],
    fractions: Vec<f64>,
    raster_bounds: Option<Bounds>,
    post: Vec<[f64; Q]>,
}

/// The native exchange and its actual application point for one lattice tick.
#[derive(Debug, Clone, Copy, PartialEq)]
#[must_use]
pub struct MovingDiskStep2 {
    /// Original finite-inertia collision receipt, with its couple about
    /// `collision_center`. Its work uses the solved midpoint velocities.
    pub exchange: PartialSaturationBodyExchange2,
    /// Body centre after the first half-drift, where the impulse is applied.
    pub collision_center: [f64; 2],
    /// Applied external angular impulse about the fixed lattice origin:
    /// `collision_center cross external_impulse + external_couple`.
    pub external_angular_impulse_about_origin: f64,
}

#[derive(Debug, Clone, Copy)]
struct Bounds {
    lower: [usize; 2],
    upper: [usize; 2],
}

impl MovingDisk2 {
    /// Bind a finite disk to a grid's dimensions, without changing the fluid.
    ///
    /// The angle records the disk's material orientation. It does not alter
    /// circular geometry; spin still enters every solid point velocity and
    /// the original reciprocal hydrodynamic torque.
    ///
    /// # Panics
    /// Invalid finite state/inertia/radius, malformed grid dimensions, or an
    /// initial swept box that touches nonfluid cells or the domain halo.
    #[must_use]
    pub fn new(grid: &Grid, radius: f64, body: PartialSaturationBody2, angle: f64) -> Self {
        assert!(
            grid.nx >= 3 && grid.ny >= 3,
            "a disk needs an interior fluid domain"
        );
        let count = grid
            .nx
            .checked_mul(grid.ny)
            .expect("disk grid size overflowed");
        assert_eq!(
            grid.flags.len(),
            count,
            "disk grid flags must cover its dimensions"
        );
        assert_eq!(
            grid.f.len(),
            count,
            "disk populations must cover its dimensions"
        );
        admit_radius(radius);
        assert!(
            body.mass.is_finite()
                && body.mass > 0.0
                && body.inertia.is_finite()
                && body.inertia > 0.0
                && body
                    .center
                    .into_iter()
                    .chain(body.velocity)
                    .chain([body.angular_velocity, angle])
                    .all(f64::is_finite),
            "disk state and positive inertia must be finite"
        );
        let disk = Self {
            body,
            radius,
            angle,
            dims: [grid.nx, grid.ny],
            fractions: vec![0.0; count],
            raster_bounds: None,
            post: Vec::new(),
        };
        disk.admit_drift(grid, body.center, body.center);
        disk
    }

    /// Accepted centre, velocities and added-solid inertia.
    #[must_use]
    pub fn body(&self) -> PartialSaturationBody2 {
        self.body
    }

    /// Accepted material orientation, counterclockwise radians.
    #[must_use]
    pub fn angle(&self) -> f64 {
        self.angle
    }

    /// Fixed disk radius, in lattice lengths.
    #[must_use]
    pub fn radius(&self) -> f64 {
        self.radius
    }

    /// Analytic area fraction of a unit cell at the accepted body position.
    /// This observation neither advances nor changes cached geometry.
    ///
    /// # Panics
    /// An index outside the bound grid dimensions.
    #[must_use]
    pub fn solid_fraction(&self, x: usize, y: usize) -> f64 {
        assert!(
            x < self.dims[0] && y < self.dims[1],
            "disk cell index out of range"
        );
        cell_fraction(
            self.radius,
            [
                x as f64 - self.body.center[0],
                y as f64 - self.body.center[1],
            ],
        )
    }

    /// Advance both body pose/velocity and actual fluid populations one tick.
    ///
    /// External linear impulse acts at the collision centre; the angular
    /// impulse is an additional couple. The original finite-inertia owner
    /// performs one 3-by-3 response solve and one collision evaluation. This
    /// adapter rasterizes only the disk's bounding box and clears only the
    /// previous box. It does not allocate another full fluid candidate.
    ///
    /// Body and fluid publish only after the final drift is admitted. Scratch
    /// may change on refusal; subsequent calls rebuild it from accepted state.
    /// Native stationary boundaries may exchange further momentum and periodic
    /// fluid seams transport angular momentum; the returned exchange covers
    /// the immersed collision, not those distinct outer-boundary fluxes.
    /// No order claim for the complete moving-boundary fluid discretization
    /// or calibrated sedimentation/drag prediction is made by this adapter.
    ///
    /// # Panics
    /// The original partial-saturation admission failures, changed dimensions,
    /// invalid impulses, or either half-drift leaving the admitted fluid region.
    pub fn step(
        &mut self,
        grid: &mut Grid,
        external_impulse: [f64; 2],
        external_couple: f64,
    ) -> MovingDiskStep2 {
        assert_eq!(
            [grid.nx, grid.ny],
            self.dims,
            "disk grid dimensions changed"
        );
        assert_eq!(
            grid.flags.len(),
            self.fractions.len(),
            "disk grid flags changed size"
        );
        assert!(
            external_impulse
                .into_iter()
                .chain([external_couple])
                .all(f64::is_finite),
            "disk external impulses must be finite"
        );
        let mut candidate = self.body;
        let collision_center =
            core::array::from_fn(|axis| self.body.center[axis] + 0.5 * self.body.velocity[axis]);
        self.admit_drift(grid, self.body.center, collision_center);
        self.rasterize(collision_center);
        candidate.center = collision_center;
        let exchange = grid.collide_into_with_partial_saturation_body(
            &mut self.post,
            &self.fractions,
            &mut candidate,
            external_impulse,
            external_couple,
        );
        let endpoint =
            core::array::from_fn(|axis| collision_center[axis] + 0.5 * candidate.velocity[axis]);
        let angle =
            self.angle + 0.5 * self.body.angular_velocity + 0.5 * candidate.angular_velocity;
        let angular_impulse = collision_center[0].mul_add(
            external_impulse[1],
            -collision_center[1] * external_impulse[0],
        ) + external_couple;
        assert!(
            angle.is_finite() && angular_impulse.is_finite(),
            "disk pose or angular impulse overflowed"
        );
        self.admit_drift(grid, collision_center, endpoint);
        // The original collision validated all dimensions and populations.
        // Streaming is an admitted copy/bounce pass with no later refusal.
        grid.stream_from(&self.post);
        candidate.center = endpoint;
        self.body = candidate;
        self.angle = angle;
        MovingDiskStep2 {
            exchange,
            collision_center,
            external_angular_impulse_about_origin: angular_impulse,
        }
    }

    fn admit_drift(&self, grid: &Grid, from: [f64; 2], to: [f64; 2]) {
        for axis in 0..2 {
            assert!(
                from[axis].is_finite()
                    && to[axis].is_finite()
                    && from[axis].min(to[axis]) - self.radius >= 0.5
                    && from[axis].max(to[axis]) + self.radius <= self.dims[axis] as f64 - 1.5,
                "moving disk left the interior lattice domain"
            );
        }
        let bounds = Bounds {
            lower: core::array::from_fn(|axis| {
                (from[axis].min(to[axis]) - self.radius + 0.5).floor() as usize
            }),
            upper: core::array::from_fn(|axis| {
                (from[axis].max(to[axis]) + self.radius + 0.5).ceil() as usize
            }),
        };
        for y in bounds.lower[1]..bounds.upper[1] {
            for x in bounds.lower[0]..bounds.upper[0] {
                assert_eq!(
                    grid.flags[grid.idx(x, y)],
                    Cell::Fluid,
                    "moving disk swept box requires fluid cells; solid contact is not modeled"
                );
            }
        }
    }

    fn rasterize(&mut self, center: [f64; 2]) {
        if let Some(old) = self.raster_bounds {
            for y in old.lower[1]..old.upper[1] {
                let row = y * self.dims[0];
                self.fractions[row + old.lower[0]..row + old.upper[0]].fill(0.0);
            }
        }
        let bounds = Bounds {
            lower: core::array::from_fn(|axis| (center[axis] - self.radius + 0.5).floor() as usize),
            upper: core::array::from_fn(|axis| (center[axis] + self.radius + 0.5).ceil() as usize),
        };
        // Record the scratch support before filling so a numerical refusal
        // still lets the next call erase the entire attempted raster.
        self.raster_bounds = Some(bounds);
        for y in bounds.lower[1]..bounds.upper[1] {
            for x in bounds.lower[0]..bounds.upper[0] {
                self.fractions[y * self.dims[0] + x] =
                    cell_fraction(self.radius, [x as f64 - center[0], y as f64 - center[1]]);
            }
        }
    }
}

/// Area of a disk inside the unit square centred at `cell_offset` from its
/// centre. Interior/exterior cells return exactly one/zero; boundary cells
/// integrate clipped triangle and circular-sector areas analytically.
///
/// # Panics
/// Nonfinite offsets or a radius whose positive square is not representable.
#[must_use]
pub fn disk_cell_fraction(radius: f64, cell_offset: [f64; 2]) -> f64 {
    admit_radius(radius);
    assert!(
        cell_offset.into_iter().all(f64::is_finite),
        "disk cell offset must be finite"
    );
    cell_fraction(radius, cell_offset)
}

fn admit_radius(radius: f64) {
    assert!(
        radius.is_finite()
            && radius > 0.0
            && (radius * radius).is_finite()
            && radius * radius > 0.0,
        "disk radius and its square must be finite and positive"
    );
}

fn cell_fraction(radius: f64, [x, y]: [f64; 2]) -> f64 {
    let r2 = radius * radius;
    let nearest = [(x.abs() - 0.5).max(0.0), (y.abs() - 0.5).max(0.0)];
    if nearest[0].mul_add(nearest[0], nearest[1] * nearest[1]) >= r2 {
        return 0.0;
    }
    let farthest = [x.abs() + 0.5, y.abs() + 0.5];
    if farthest[0].mul_add(farthest[0], farthest[1] * farthest[1]) <= r2 {
        return 1.0;
    }
    let corners = [
        [x - 0.5, y - 0.5],
        [x + 0.5, y - 0.5],
        [x + 0.5, y + 0.5],
        [x - 0.5, y + 0.5],
    ];
    let mut area = 0.0;
    for edge in 0..4 {
        let a = corners[edge];
        let b = corners[(edge + 1) % 4];
        let axis = usize::from(a[0] == b[0]);
        let fixed = a[1 - axis];
        let delta = b[axis] - a[axis];
        let mut split = [0.0, 1.0, 0.0, 0.0];
        let mut count = 2;
        if fixed.abs() < radius {
            let crossing = det::sqrt((radius - fixed.abs()) * (radius + fixed.abs()));
            for coordinate in [-crossing, crossing] {
                let t = (coordinate - a[axis]) / delta;
                if t > 0.0 && t < 1.0 {
                    split[count] = t;
                    count += 1;
                }
            }
        }
        split[..count].sort_by(f64::total_cmp);
        let d = [b[0] - a[0], b[1] - a[1]];
        let edge_cross = a[0] * d[1] - a[1] * d[0];
        for segment in split[..count].windows(2) {
            let u = [a[0] + segment[0] * d[0], a[1] + segment[0] * d[1]];
            let v = [a[0] + segment[1] * d[0], a[1] + segment[1] * d[1]];
            let middle = [0.5 * (u[0] + v[0]), 0.5 * (u[1] + v[1])];
            let cross = edge_cross * (segment[1] - segment[0]);
            // A tangent edge has its midpoint on the circle but otherwise
            // lies outside it; only a strictly interior segment is triangular.
            area += if middle[0].mul_add(middle[0], middle[1] * middle[1]) < r2 {
                0.5 * cross
            } else {
                0.5 * r2 * det::atan2(cross, u[0].mul_add(v[0], u[1] * v[1]))
            };
        }
    }
    let tolerance = 128.0 * f64::EPSILON * r2.max(1.0);
    assert!(
        area.is_finite() && area >= -tolerance && area <= 1.0 + tolerance,
        "analytic disk/cell area is numerically unresolved"
    );
    area.clamp(0.0, 1.0)
}
