//! SI solid geometry sampled into the original stationary NSF cell mask.
//!
//! The existing F-rep owner evaluates boxes and hard half-space intersections.
//! Coordinates start at the low domain wall; cell centers are `(i+1/2)*dx`,
//! as in the original reservoir observer. Input order applies additions and
//! subtractions. Boundary centers belong to the declared shape. The result is
//! the existing halfway bounce-back staircase, not a cut-cell/off-lattice wall
//! or a certified geometry approximation. Keep these SI shapes fixed while
//! refining dx; a shape selecting no cell is refused instead of disappearing.

use fs_exec::CancelGate;
use fs_geom::{Aabb, Point3, Vec3};
use fs_rep_frep::{BoolOp, BoolStyle, FrepBuilder};

const MAX_SHAPES: usize = 256;
const MAX_CELL_EVALUATIONS: usize = 16_777_216;

#[derive(Clone, Copy, Debug)]
pub(super) struct Primitive {
    pub add: bool,
    pub bounds: Aabb,
    /// Unit-normal distance from the origin. F-rep normalizes the supplied
    /// direction; its offset remains a physical distance in metres.
    pub clip: Option<(Vec3, f64)>,
}

pub(super) fn admit_work(
    dims: [usize; 3],
    dx_m: f64,
    count: usize,
    maximum_cell_evaluations: usize,
) -> Result<usize, &'static str> {
    let cells = dims.iter().try_fold(1usize, |n, &d| n.checked_mul(d));
    if dims.contains(&0)
        || !dx_m.is_finite()
        || dx_m <= 0.0
        || dims.iter().any(|&n| !(n as f64 * dx_m).is_finite())
        || !(1..=MAX_SHAPES).contains(&count)
        || maximum_cell_evaluations == 0
        || maximum_cell_evaluations > MAX_CELL_EVALUATIONS
        || cells
            .and_then(|n| n.checked_mul(count))
            .is_none_or(|n| n > maximum_cell_evaluations)
    {
        return Err(
            "SI geometry requires finite dimensions, 1..256 shapes and a complete bounded cell-evaluation allowance",
        );
    }
    Ok(cells.expect("admitted positive bounded work"))
}

pub(super) fn rasterize(
    dims: [usize; 3],
    dx_m: f64,
    primitives: &[Primitive],
    maximum_cell_evaluations: usize,
    gate: &CancelGate,
) -> Result<Vec<bool>, &'static str> {
    let cells = admit_work(dims, dx_m, primitives.len(), maximum_cell_evaluations)?;
    if gate.is_requested() {
        return Err("SI geometry construction cancelled");
    }
    if !primitives[0].add {
        return Err("SI solid geometry must start with an addition");
    }
    let extent = dims.map(|n| n as f64 * dx_m);
    let mut fields = Vec::with_capacity(primitives.len());
    // Admit every source shape before allocating or mutating the cell mask.
    for primitive in primitives {
        let Aabb { min, max } = primitive.bounds;
        let lo = [min.x, min.y, min.z];
        let hi = [max.x, max.y, max.z];
        if (0..3).any(|a| {
            !lo[a].is_finite()
                || !hi[a].is_finite()
                || lo[a] < 0.0
                || hi[a] > extent[a]
                || lo[a] >= hi[a]
        }) {
            return Err("SI solid bounds must have positive finite spans inside the physical grid");
        }
        let mut builder = FrepBuilder::new();
        let center = Point3::new(
            f64::midpoint(min.x, max.x),
            f64::midpoint(min.y, max.y),
            f64::midpoint(min.z, max.z),
        );
        let half = Vec3::new(
            0.5 * (max.x - min.x),
            0.5 * (max.y - min.y),
            0.5 * (max.z - min.z),
        );
        let mut root = builder
            .box_prim(center, half)
            .map_err(|_| "original F-rep refused the SI solid box")?;
        if let Some((normal, offset)) = primitive.clip {
            let plane = builder
                .half_space(normal, offset)
                .map_err(|_| "original F-rep refused the SI wedge normal or offset")?;
            root = builder
                .boolean(BoolOp::Intersect, BoolStyle::Hard, root, plane)
                .map_err(|_| "original F-rep refused the SI wedge intersection")?;
        }
        fields.push(
            builder
                .finish(root)
                .map_err(|_| "original F-rep refused the SI solid geometry")?,
        );
    }
    let mut solid = vec![false; cells];
    for (primitive, field) in primitives.iter().zip(fields) {
        let mut selected = 0usize;
        for (i, value) in solid.iter_mut().enumerate() {
            if i % 256 == 0 && gate.is_requested() {
                return Err("SI geometry construction cancelled");
            }
            let x = [i % dims[0], i / dims[0] % dims[1], i / (dims[0] * dims[1])];
            let p = x.map(|n| (n as f64 + 0.5) * dx_m);
            let distance = field.value(Point3::new(p[0], p[1], p[2]));
            if !distance.is_finite() {
                return Err("SI solid field evaluation must remain finite");
            }
            if distance <= 0.0 {
                selected += 1;
                *value = primitive.add;
            }
        }
        if selected == 0 {
            return Err(
                "SI solid shape selects no cell centers; refine the physical grid or correct its bounds",
            );
        }
    }
    if solid.iter().all(|&value| value) {
        return Err("SI solid geometry must leave actual fluid cells");
    }
    Ok(solid)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn shape(add: bool, min: [f64; 3], max: [f64; 3]) -> Primitive {
        Primitive {
            add,
            bounds: Aabb {
                min: Point3::new(min[0], min[1], min[2]),
                max: Point3::new(max[0], max[1], max[2]),
            },
            clip: None,
        }
    }

    #[test]
    fn g0_original_frep_solids_and_subtracted_slit_use_physical_cell_centers() {
        let solids = [
            shape(true, [1.0, 0.0, 0.0], [3.0, 4.0, 2.0]),
            shape(false, [1.0, 1.0, 0.0], [3.0, 2.0, 2.0]),
        ];
        let mask = rasterize([4, 4, 2], 1.0, &solids, 64, &CancelGate::new_clock_free()).unwrap();
        for z in 0..2 {
            for y in 0..4 {
                for x in 0..4 {
                    assert_eq!(mask[(z * 4 + y) * 4 + x], (1..3).contains(&x) && y != 1);
                }
            }
        }
    }

    #[test]
    fn g3_same_si_labium_converges_under_grid_refinement_without_changing_shape() {
        let mut wedge = shape(true, [1.0, 0.0, 0.0], [3.0, 2.0, 1.0]);
        // x + 0.37*y <= 2.7; analytic prism volume is 2.66 m^3.
        wedge.clip = Some((
            Vec3::new(1.0, 0.37, 0.0),
            2.7 / (1.0_f64 + 0.37 * 0.37).sqrt(),
        ));
        let volume = |dims, dx| {
            let mask = rasterize(dims, dx, &[wedge], 4096, &CancelGate::new_clock_free()).unwrap();
            mask.iter().filter(|&&value| value).count() as f64 * dx * dx * dx
        };
        let coarse = volume([8, 4, 2], 0.5);
        let fine = volume([16, 8, 4], 0.25);
        assert!((coarse - 2.5).abs() < 1e-12);
        assert!((fine - 2.6875).abs() < 1e-12);
        assert!((fine - 2.66).abs() < (coarse - 2.66).abs());
    }

    #[test]
    fn g0_g4_si_geometry_refuses_missing_resolution_budget_domain_and_cancellation() {
        let gate = CancelGate::new_clock_free();
        let solid = shape(true, [1.0, 0.0, 0.0], [3.0, 4.0, 2.0]);
        assert!(rasterize([4, 4, 2], 1.0, &[solid], 31, &gate).is_err());
        assert!(admit_work([usize::MAX, 2, 1], 1.0, 1, usize::MAX).is_err());
        assert!(
            rasterize(
                [4, 4, 2],
                1.0,
                &[shape(true, [0.0; 3], [5.0, 4.0, 2.0])],
                32,
                &gate
            )
            .is_err()
        );
        assert!(
            rasterize(
                [4, 4, 2],
                1.0,
                &[shape(true, [1.1, 0.0, 0.0], [1.2, 4.0, 2.0])],
                32,
                &gate
            )
            .is_err()
        );
        let mut malformed = solid;
        malformed.clip = Some((Vec3::new(0.0, 0.0, 0.0), 0.0));
        assert!(rasterize([4, 4, 2], 1.0, &[malformed], 32, &gate).is_err());
        gate.request();
        assert!(rasterize([4, 4, 2], 1.0, &[solid], 32, &gate).is_err());
    }
}
