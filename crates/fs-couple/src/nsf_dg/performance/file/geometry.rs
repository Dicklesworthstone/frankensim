//! SI solid geometry sampled into the original stationary NSF cell mask.
//!
//! The existing F-rep owner evaluates boxes, half-space intersections, spheres
//! and capped cylinders. Cylinders retain their authored physical endpoints;
//! the original rigid transform maps a capped +z primitive into that frame.
//! Coordinates start at the low domain wall; cell centers are `(i+1/2)*dx`,
//! as in the original reservoir observer. Input order applies additions and
//! subtractions. Boundary centers belong to the declared shape. The result is
//! the existing halfway bounce-back staircase, not a cut-cell/off-lattice wall
//! or a certified geometry approximation. Keep these SI shapes fixed while
//! refining dx; a shape selecting no cell is refused instead of disappearing.

use fs_exec::CancelGate;
use fs_geom::{Aabb, Point3, Vec3};
use fs_math::det;
use fs_rep_frep::{BoolOp, BoolStyle, Frep, FrepBuilder};

const MAX_SHAPES: usize = 256;
const MAX_CELL_EVALUATIONS: usize = 16_777_216;

/// Sample one physical face rectangle on the two increasing Cartesian axes
/// other than `axis`, exactly as the original characteristic-patch owner.
/// Bounds are half-open in metres from the grid's low wall. Face tiles are
/// selected by their centers, so each represented edge is within dx/2 of its
/// authored position. This is a staircase aperture, not subcell quadrature.
pub(super) fn face_rectangle(
    dims: [usize; 3],
    axis: usize,
    dx_m: f64,
    bounds: [f64; 4],
    gate: &CancelGate,
) -> Result<([usize; 2], [usize; 2]), &'static str> {
    admit_work(dims, dx_m, 1, MAX_CELL_EVALUATIONS)?;
    if axis >= 3 {
        return Err("physical aperture requires a Cartesian normal axis");
    }
    let mut start = [0; 2];
    let mut end = [0; 2];
    let mut tangent = 0;
    for (a, count) in dims.into_iter().enumerate() {
        if a == axis {
            continue;
        }
        let low = bounds[2 * tangent];
        let high = bounds[2 * tangent + 1];
        if !low.is_finite()
            || !high.is_finite()
            || low < 0.0
            || high > count as f64 * dx_m
            || low >= high
        {
            return Err(
                "physical aperture bounds must have positive finite spans inside the grid face",
            );
        }
        let mut first = None;
        for i in 0..count {
            if i % 256 == 0 && gate.is_requested() {
                return Err("physical aperture construction cancelled");
            }
            let center = (i as f64 + 0.5) * dx_m;
            if center >= low && center < high {
                first.get_or_insert(i);
                end[tangent] = i + 1;
            }
        }
        start[tangent] = first.ok_or(
            "physical aperture selects no face tiles; refine the grid or correct its bounds",
        )?;
        tangent += 1;
    }
    Ok((start, end))
}

#[derive(Clone, Copy, Debug)]
pub(super) struct Primitive {
    pub add: bool,
    pub shape: Shape,
}

#[derive(Clone, Copy, Debug)]
pub(super) enum Shape {
    Box {
        bounds: Aabb,
        /// Unit-normal distance from the origin. F-rep normalizes the supplied
        /// direction; its offset remains a physical distance in metres.
        clip: Option<(Vec3, f64)>,
    },
    Cylinder {
        start: Point3,
        end: Point3,
        radius_m: f64,
    },
    Sphere {
        center: Point3,
        radius_m: f64,
    },
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

fn field(shape: Shape, extent: [f64; 3]) -> Result<Frep, &'static str> {
    let mut builder = FrepBuilder::new();
    let (bounds, root) = match shape {
        Shape::Box { bounds, clip } => {
            let Aabb { min, max } = bounds;
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
            if let Some((normal, offset)) = clip {
                let plane = builder
                    .half_space(normal, offset)
                    .map_err(|_| "original F-rep refused the SI wedge normal or offset")?;
                root = builder
                    .boolean(BoolOp::Intersect, BoolStyle::Hard, root, plane)
                    .map_err(|_| "original F-rep refused the SI wedge intersection")?;
            }
            (bounds, root)
        }
        Shape::Sphere { center, radius_m } => {
            let root = builder
                .sphere(center, radius_m)
                .map_err(|_| "original F-rep refused the SI sphere center or radius")?;
            (
                Aabb {
                    min: Point3::new(
                        center.x - radius_m,
                        center.y - radius_m,
                        center.z - radius_m,
                    ),
                    max: Point3::new(
                        center.x + radius_m,
                        center.y + radius_m,
                        center.z + radius_m,
                    ),
                },
                root,
            )
        }
        Shape::Cylinder {
            start,
            end,
            radius_m,
        } => {
            let a = [start.x, start.y, start.z];
            let b = [end.x, end.y, end.z];
            let delta = [end.x - start.x, end.y - start.y, end.z - start.z];
            let length = det::hypot(det::hypot(delta[0], delta[1]), delta[2]);
            if a.iter().chain(&b).any(|v| !v.is_finite())
                || !length.is_finite()
                || length <= 0.0
                || !radius_m.is_finite()
                || radius_m <= 0.0
            {
                return Err(
                    "SI cylinder requires distinct finite endpoints and a positive finite radius",
                );
            }
            let direction = delta.map(|v| v / length);
            // Projection of the circular cross-section onto each Cartesian
            // axis. Endpoint +/- radius would unnecessarily reject an axial
            // end cap exactly on a physical domain boundary.
            let radial: [f64; 3] = core::array::from_fn(|i| {
                radius_m * det::hypot(direction[(i + 1) % 3], direction[(i + 2) % 3])
            });
            let low: [f64; 3] = core::array::from_fn(|i| a[i].min(b[i]) - radial[i]);
            let high: [f64; 3] = core::array::from_fn(|i| a[i].max(b[i]) + radial[i]);
            let origin = Point3::new(0.0, 0.0, 0.0);
            let cylinder = builder
                .cylinder(origin, radius_m)
                .map_err(|_| "original F-rep refused the SI cylinder")?;
            let caps = builder
                .box_prim(origin, Vec3::new(radius_m, radius_m, 0.5 * length))
                .map_err(|_| "original F-rep refused the SI cylinder caps")?;
            let mut root = builder
                .boolean(BoolOp::Intersect, BoolStyle::Hard, cylinder, caps)
                .map_err(|_| "original F-rep refused the capped SI cylinder")?;
            let transverse = det::hypot(direction[0], direction[1]);
            // A finite cylinder is unchanged by reversing its axis, so either
            // axial orientation already has the required support and caps.
            if transverse > 0.0 {
                root = builder
                    .rotate(
                        root,
                        Vec3::new(-direction[1] / transverse, direction[0] / transverse, 0.0),
                        det::atan2(transverse, direction[2]),
                    )
                    .map_err(|_| "original F-rep refused the SI cylinder orientation")?;
            }
            root = builder
                .translate(
                    root,
                    Vec3::new(
                        f64::midpoint(start.x, end.x),
                        f64::midpoint(start.y, end.y),
                        f64::midpoint(start.z, end.z),
                    ),
                )
                .map_err(|_| "original F-rep refused the SI cylinder placement")?;
            (
                Aabb {
                    min: Point3::new(low[0], low[1], low[2]),
                    max: Point3::new(high[0], high[1], high[2]),
                },
                root,
            )
        }
    };
    let Aabb { min, max } = bounds;
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
    builder
        .finish(root)
        .map_err(|_| "original F-rep refused the SI solid geometry")
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
        fields.push(field(primitive.shape, extent)?);
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
            shape: Shape::Box {
                bounds: Aabb {
                    min: Point3::new(min[0], min[1], min[2]),
                    max: Point3::new(max[0], max[1], max[2]),
                },
                clip: None,
            },
        }
    }

    #[test]
    fn g3_physical_face_rectangles_keep_area_axes_and_half_open_edges_under_refinement() {
        let gate = CancelGate::new_clock_free();
        for axis in 0..3 {
            let bounds = [0.375, 0.875, 0.125, 0.625];
            let (a, b) = face_rectangle([8; 3], axis, 0.25, bounds, &gate).unwrap();
            assert_eq!((a, b), ([1, 0], [3, 2]));
            let (c, d) = face_rectangle([16; 3], axis, 0.125, bounds, &gate).unwrap();
            assert_eq!((c, d), ([3, 1], [7, 5]));
            let area = |start: [usize; 2], end: [usize; 2], dx: f64| {
                (end[0] - start[0]) as f64 * (end[1] - start[1]) as f64 * dx * dx
            };
            assert_eq!(area(a, b, 0.25), 0.25);
            assert_eq!(area(c, d, 0.125), 0.25);
            for (start, end, dx) in [(a, b, 0.25), (c, d, 0.125)] {
                for j in 0..2 {
                    assert!((start[j] as f64 * dx - bounds[2 * j]).abs() <= 0.5 * dx);
                    assert!((end[j] as f64 * dx - bounds[2 * j + 1]).abs() <= 0.5 * dx);
                }
            }
        }
        // Axis 1 uses x,z in increasing Cartesian order, not cyclic z,x.
        assert_eq!(
            face_rectangle([12, 8, 6], 1, 0.25, [2.0, 2.5, 1.0, 1.25], &gate).unwrap(),
            ([8, 4], [10, 5])
        );
    }

    #[test]
    fn g0_g4_physical_apertures_refuse_lost_resolution_domain_and_cancellation() {
        let gate = CancelGate::new_clock_free();
        for bounds in [
            [0.2, 0.3, 0.0, 1.0],
            [1.0, 1.0, 0.0, 1.0],
            [-0.1, 1.0, 0.0, 1.0],
            [0.0, 5.0, 0.0, 1.0],
            [f64::NAN, 1.0, 0.0, 1.0],
        ] {
            assert!(face_rectangle([4; 3], 0, 1.0, bounds, &gate).is_err());
        }
        assert!(face_rectangle([4; 3], 3, 1.0, [0.0, 1.0, 0.0, 1.0], &gate).is_err());
        gate.request();
        assert!(face_rectangle([4; 3], 0, 1.0, [0.0, 1.0, 0.0, 1.0], &gate).is_err());
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
    fn g0_capped_curved_solids_match_independent_physical_membership() {
        let start = Point3::new(0.83, 1.07, 0.71);
        let end = Point3::new(3.11, 2.89, 2.93);
        let delta = [end.x - start.x, end.y - start.y, end.z - start.z];
        let length_squared: f64 = delta.iter().map(|v| v * v).sum();
        let gate = CancelGate::new_clock_free();
        // Endpoint reversal, arbitrary orientation and finite caps must all
        // describe the same physical cylinder, not an infinite tube or capsule.
        for (a, b) in [(start, end), (end, start)] {
            let solid = Primitive {
                add: true,
                shape: Shape::Cylinder {
                    start: a,
                    end: b,
                    radius_m: 0.41,
                },
            };
            let mask = rasterize([16; 3], 0.25, &[solid], 4096, &gate).unwrap();
            for (i, &actual) in mask.iter().enumerate() {
                let cell = [i % 16, i / 16 % 16, i / 256];
                let p = cell.map(|n| (n as f64 + 0.5) * 0.25);
                let v = [p[0] - start.x, p[1] - start.y, p[2] - start.z];
                let dot: f64 = (0..3).map(|j| v[j] * delta[j]).sum();
                let radial_squared =
                    v.iter().map(|v| v * v).sum::<f64>() - dot * dot / length_squared;
                let expected =
                    (0.0..=length_squared).contains(&dot) && radial_squared <= 0.41_f64.powi(2);
                assert_eq!(actual, expected, "cell {cell:?}");
            }
        }

        let outer = Shape::Cylinder {
            start: Point3::new(0.75, 2.0, 2.0),
            end: Point3::new(3.25, 2.0, 2.0),
            radius_m: 0.8,
        };
        let inner = Shape::Cylinder {
            start: Point3::new(0.75, 2.0, 2.0),
            end: Point3::new(3.25, 2.0, 2.0),
            radius_m: 0.4,
        };
        let sphere = Shape::Sphere {
            center: Point3::new(2.0, 2.0, 2.0),
            radius_m: 0.55,
        };
        let primitives = [
            Primitive {
                add: true,
                shape: outer,
            },
            Primitive {
                add: false,
                shape: inner,
            },
            Primitive {
                add: true,
                shape: sphere,
            },
        ];
        let mask = rasterize([16; 3], 0.25, &primitives, 3 * 4096, &gate).unwrap();
        for (i, &actual) in mask.iter().enumerate() {
            let cell = [i % 16, i / 16 % 16, i / 256];
            let p = cell.map(|n| (n as f64 + 0.5) * 0.25);
            let radial = (p[1] - 2.0).hypot(p[2] - 2.0);
            let tube = (0.75..=3.25).contains(&p[0]) && radial <= 0.8 && radial > 0.4;
            let cap = radial.hypot(p[0] - 2.0) <= 0.55;
            assert_eq!(actual, tube || cap, "ordered curved CSG at {cell:?}");
        }
    }

    #[test]
    fn g3_curved_physical_volumes_refine_without_changing_authored_geometry() {
        for (shape, exact) in [
            (
                Shape::Sphere {
                    center: Point3::new(2.0, 2.0, 2.0),
                    radius_m: 0.91,
                },
                4.0 * core::f64::consts::PI * 0.91_f64.powi(3) / 3.0,
            ),
            (
                Shape::Cylinder {
                    start: Point3::new(0.75, 2.0, 2.0),
                    end: Point3::new(3.25, 2.0, 2.0),
                    radius_m: 0.73,
                },
                core::f64::consts::PI * 0.73_f64.powi(2) * 2.5,
            ),
        ] {
            let mut errors = Vec::new();
            for n in [16, 32, 64] {
                let dx = 4.0 / n as f64;
                let mask = rasterize(
                    [n; 3],
                    dx,
                    &[Primitive { add: true, shape }],
                    n * n * n,
                    &CancelGate::new_clock_free(),
                )
                .unwrap();
                let volume = mask.iter().filter(|&&v| v).count() as f64 * dx.powi(3);
                errors.push((volume - exact).abs());
            }
            assert!(
                errors.windows(2).all(|pair| pair[1] < pair[0]),
                "{shape:?}: {errors:?}"
            );
            assert!(errors[2] / exact < 0.01, "{shape:?}: {errors:?}");
        }
    }

    #[test]
    fn g0_g4_curved_solids_refuse_invalid_domain_radius_resolution_and_cancellation() {
        let gate = CancelGate::new_clock_free();
        let sphere = |radius_m| Shape::Sphere {
            center: Point3::new(2.0, 2.0, 2.0),
            radius_m,
        };
        let cylinder = |start, end, radius_m| Shape::Cylinder {
            start,
            end,
            radius_m,
        };
        for shape in [
            sphere(0.0),
            sphere(-1.0),
            sphere(f64::NAN),
            sphere(3.0),
            sphere(0.01),
            cylinder(Point3::new(2.0, 2.0, 1.0), Point3::new(2.0, 2.0, 1.0), 0.5),
            cylinder(Point3::new(2.0, 2.0, 1.0), Point3::new(2.0, 2.0, 3.0), -0.5),
            cylinder(
                Point3::new(2.0, 2.0, 1.0),
                Point3::new(2.0, 2.0, f64::INFINITY),
                0.5,
            ),
            cylinder(Point3::new(0.1, 0.1, 1.0), Point3::new(0.2, 0.2, 3.0), 0.5),
        ] {
            assert!(
                rasterize(
                    [16; 3],
                    0.25,
                    &[Primitive { add: true, shape }],
                    4096,
                    &gate
                )
                .is_err(),
                "{shape:?}"
            );
        }
        // End caps may meet a domain boundary; a spurious +/- radius bound
        // along the cylinder axis must not crop or reject them.
        for (start, end) in [
            (Point3::new(2.0, 2.0, 0.0), Point3::new(2.0, 2.0, 4.0)),
            (Point3::new(2.0, 2.0, 4.0), Point3::new(2.0, 2.0, 0.0)),
        ] {
            assert!(
                rasterize(
                    [16; 3],
                    0.25,
                    &[Primitive {
                        add: true,
                        shape: cylinder(start, end, 0.5)
                    }],
                    4096,
                    &gate
                )
                .is_ok()
            );
        }
        gate.request();
        assert!(
            rasterize(
                [16; 3],
                0.25,
                &[Primitive {
                    add: true,
                    shape: sphere(0.9)
                }],
                4096,
                &gate
            )
            .is_err()
        );
    }

    #[test]
    fn g3_same_si_labium_converges_under_grid_refinement_without_changing_shape() {
        let mut wedge = shape(true, [1.0, 0.0, 0.0], [3.0, 2.0, 1.0]);
        // x + 0.37*y <= 2.7; analytic prism volume is 2.66 m^3.
        let Shape::Box { clip, .. } = &mut wedge.shape else {
            unreachable!()
        };
        *clip = Some((
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
        let Shape::Box { clip, .. } = &mut malformed.shape else {
            unreachable!()
        };
        *clip = Some((Vec3::new(0.0, 0.0, 0.0), 0.0));
        assert!(rasterize([4, 4, 2], 1.0, &[malformed], 32, &gate).is_err());
        gate.request();
        assert!(rasterize([4, 4, 2], 1.0, &[solid], 32, &gate).is_err());
    }
}
