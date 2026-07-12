//! `frankensim-gmik` vector CutFEM elasticity acceptance battery.
//!
//! - cte-001 G0 affine patch law on exactly represented closed polygons;
//! - cte-002 G1 curved-domain MMS convergence;
//! - cte-003 cut-fraction-independent conditioning under ghost penalty.
//!
//! The feature-required cte-004 VJP gate lives in `elasticity_adjoint.rs`.

use fs_cutfem::{
    Circle, CutElasticity, CutElasticityOperator, CutFemError, CutSdf, HalfPlane, Quadtree,
    condition_estimate,
};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use std::f64::consts::PI;
use std::fmt::Write as _;

fn verdict(name: &str, pass: bool, details: &str) {
    println!(
        "{{\"test\":\"{name}\",\"verdict\":\"{}\",{details}}}",
        if pass { "pass" } else { "fail" }
    );
    assert!(pass, "{name} failed: {details}");
}

fn json_number(value: f64) -> String {
    if value.is_finite() {
        format!("{value:.6e}")
    } else if value.is_nan() {
        "\"nan\"".to_string()
    } else if value.is_sign_positive() {
        "\"+inf\"".to_string()
    } else {
        "\"-inf\"".to_string()
    }
}

fn material() -> IsotropicElastic {
    // The largest MMS strain component is below 2*pi. Keep every acceptance
    // fixture inside the material model card rather than testing a numerically
    // correct solve in a scientifically out-of-domain regime.
    IsotropicElastic::new(1.0, 0.3, 10.0).expect("fixture material")
}

fn problem<'a>(
    grid: &'a Quadtree,
    sdf: &'a dyn CutSdf,
    material: &'a IsotropicElastic,
    ghost_gamma: f64,
) -> CutElasticity<'a> {
    CutElasticity {
        grid,
        sdf,
        material,
        nitsche_beta: 20.0,
        ghost_gamma,
        quad_depth: 3,
        clamp: None,
        boundary_traction: None,
        traction_free_interface: false,
        solver_tol: 1e-13,
        solver_max_iters: 60_000,
    }
}

/// Convex piecewise-linear level set. Vertices are placed on background-grid
/// nodes while its edges cross cells at non-axis-aligned, unequal fractions.
/// Consequently every interface segment and its normal are represented exactly
/// by the reused cut quadrature; there is no curved/corner geometry error to
/// contaminate the algebraic patch law.
#[derive(Debug, Clone)]
struct ConvexPolygon {
    /// Outward normal and offset for `normal dot x - offset <= 0`.
    /// `cut_cell_rules` normalizes the gradient at interface points.
    planes: Vec<([f64; 2], f64)>,
}

impl ConvexPolygon {
    fn from_ccw(vertices: &[[f64; 2]]) -> Self {
        assert!(vertices.len() >= 3, "a polygon needs three vertices");
        let signed_double_area: f64 = vertices
            .iter()
            .zip(vertices.iter().cycle().skip(1))
            .take(vertices.len())
            .map(|(a, b)| a[0] * b[1] - a[1] * b[0])
            .sum();
        assert!(
            signed_double_area > 0.0,
            "vertices must be counter-clockwise"
        );
        let planes = vertices
            .iter()
            .zip(vertices.iter().cycle().skip(1))
            .take(vertices.len())
            .map(|(a, b)| {
                let edge = [b[0] - a[0], b[1] - a[1]];
                assert!(edge[0] != 0.0 || edge[1] != 0.0, "duplicate polygon vertex");
                // Keep dyadic coefficients dyadic. Normalization here would
                // introduce irrational roundoff and move a nominal grid-node
                // vertex off the exact zero set; the shared chord rule already
                // normalizes gradients before using them as physical normals.
                let normal = [edge[1], -edge[0]];
                let offset = normal[0] * a[0] + normal[1] * a[1];
                (normal, offset)
            })
            .collect();
        ConvexPolygon { planes }
    }

    fn plane_values(&self, point: [f64; 2]) -> impl Iterator<Item = (f64, [f64; 2])> + '_ {
        self.planes.iter().map(move |&(normal, offset)| {
            (normal[0] * point[0] + normal[1] * point[1] - offset, normal)
        })
    }
}

impl CutSdf for ConvexPolygon {
    fn value(&self, point: [f64; 2]) -> f64 {
        self.plane_values(point)
            .map(|(value, _)| value)
            .fold(f64::NEG_INFINITY, f64::max)
    }

    fn gradient(&self, point: [f64; 2]) -> [f64; 2] {
        self.plane_values(point)
            .max_by(|a, b| a.0.total_cmp(&b.0))
            .expect("polygon planes")
            .1
    }

    fn enclose(&self, lo: [f64; 2], hi: [f64; 2]) -> Interval {
        let x = Interval::new(lo[0], hi[0]);
        let y = Interval::new(lo[1], hi[1]);
        let planes: Vec<Interval> = self
            .planes
            .iter()
            .map(|&(normal, offset)| {
                Interval::point(normal[0]) * x + Interval::point(normal[1]) * y
                    - Interval::point(offset)
            })
            .collect();
        // max(I_i) is enclosed by [max lo(I_i), max hi(I_i)]. The
        // plane intervals above are already outward-rounded by fs-ivl.
        let lower = planes
            .iter()
            .map(|interval| interval.lo())
            .fold(f64::NEG_INFINITY, f64::max);
        let upper = planes
            .iter()
            .map(|interval| interval.hi())
            .fold(f64::NEG_INFINITY, f64::max);
        Interval::new(lower, upper)
    }
}

type Displacement = fn(f64, f64) -> [f64; 2];
type DisplacementGradient = fn(f64, f64) -> [[f64; 2]; 2];

fn tx(_: f64, _: f64) -> [f64; 2] {
    [0.01, 0.0]
}
fn ty(_: f64, _: f64) -> [f64; 2] {
    [0.0, 0.01]
}
fn ux(x: f64, _: f64) -> [f64; 2] {
    [0.01 * x, 0.0]
}
fn uy(_: f64, y: f64) -> [f64; 2] {
    [0.01 * y, 0.0]
}
fn vx(x: f64, _: f64) -> [f64; 2] {
    [0.0, 0.01 * x]
}
fn vy(_: f64, y: f64) -> [f64; 2] {
    [0.0, 0.01 * y]
}
fn grad_zero(_: f64, _: f64) -> [[f64; 2]; 2] {
    [[0.0; 2]; 2]
}
fn grad_ux(_: f64, _: f64) -> [[f64; 2]; 2] {
    [[0.01, 0.0], [0.0, 0.0]]
}
fn grad_uy(_: f64, _: f64) -> [[f64; 2]; 2] {
    [[0.0, 0.01], [0.0, 0.0]]
}
fn grad_vx(_: f64, _: f64) -> [[f64; 2]; 2] {
    [[0.0, 0.0], [0.01, 0.0]]
}
fn grad_vy(_: f64, _: f64) -> [[f64; 2]; 2] {
    [[0.0, 0.0], [0.0, 0.01]]
}

#[test]
fn cte_001_constant_strain_patch_on_arbitrary_cuts() {
    let grid = Quadtree::uniform(4);
    let mat = material();
    let polygons = [
        ConvexPolygon::from_ccw(&[
            [3.0 / 16.0, 3.0 / 16.0],
            [13.0 / 16.0, 4.0 / 16.0],
            [11.0 / 16.0, 13.0 / 16.0],
        ]),
        ConvexPolygon::from_ccw(&[
            [2.0 / 16.0, 5.0 / 16.0],
            [11.0 / 16.0, 2.0 / 16.0],
            [14.0 / 16.0, 11.0 / 16.0],
            [6.0 / 16.0, 14.0 / 16.0],
        ]),
        ConvexPolygon::from_ccw(&[
            [4.0 / 16.0, 2.0 / 16.0],
            [14.0 / 16.0, 6.0 / 16.0],
            [10.0 / 16.0, 14.0 / 16.0],
            [2.0 / 16.0, 11.0 / 16.0],
        ]),
    ];
    let cases: [(&str, Displacement, DisplacementGradient); 6] = [
        ("tx", tx, grad_zero),
        ("ty", ty, grad_zero),
        ("ux", ux, grad_ux),
        ("uy", uy, grad_uy),
        ("vx", vx, grad_vx),
        ("vy", vy, grad_vy),
    ];
    let mut rows = String::new();
    let mut worst_l2 = 0.0f64;
    let mut worst_h1 = 0.0f64;
    let mut worst_affine_residual = 0.0f64;
    let mut worst_solve_residual = 0.0f64;
    let mut all_finite = true;
    for (polygon_index, sdf) in polygons.iter().enumerate() {
        for &(case, exact, gradient) in &cases {
            let cut = problem(&grid, sdf, &mat, 0.5);
            let operator = cut
                .assemble(&|_, _| [0.0, 0.0], &exact)
                .expect("patch operator");
            let mut exact_coefficients = vec![0.0; operator.dof_count()];
            for (&node, &id) in operator.node_ids() {
                let point = grid.node_pos(node);
                let value = exact(point[0], point[1]);
                exact_coefficients[2 * id] = value[0];
                exact_coefficients[2 * id + 1] = value[1];
            }
            let applied = operator.apply_vec(&exact_coefficients);
            let scale = applied
                .iter()
                .chain(operator.rhs())
                .map(|value| value.abs())
                .fold(f64::MIN_POSITIVE, f64::max);
            let affine_residual = applied
                .iter()
                .zip(operator.rhs())
                .map(|(lhs, rhs)| (lhs - rhs).abs())
                .fold(0.0f64, f64::max)
                / scale;
            let solution = cut.solve(&|_, _| [0.0, 0.0], &exact).expect("patch solve");
            let (l2, h1) = cut.l2_h1_error(&solution, &exact, &gradient);
            all_finite &= applied
                .iter()
                .chain(operator.rhs())
                .chain(exact_coefficients.iter())
                .all(|value| value.is_finite())
                && scale.is_finite()
                && affine_residual.is_finite()
                && solution.rel_residual.is_finite()
                && l2.is_finite()
                && h1.is_finite();
            worst_l2 = worst_l2.max(l2);
            worst_h1 = worst_h1.max(h1);
            worst_affine_residual = worst_affine_residual.max(affine_residual);
            worst_solve_residual = worst_solve_residual.max(solution.rel_residual);
            let _ = write!(
                rows,
                "{{\"polygon\":{polygon_index},\"case\":\"{case}\",\"affine_residual\":{affine_residual:.3e},\"solve_residual\":{:.3e},\"l2\":{l2:.3e},\"h1\":{h1:.3e},\"iters\":{}}},",
                solution.rel_residual, solution.iters
            );
        }
    }
    // The six fields span vector-valued affine displacements. The polygon
    // interface is represented exactly. Gate the discrete affine law itself
    // at roundoff separately from the condition-amplified CG forward error.
    let pass = all_finite
        && worst_affine_residual < 2e-12
        && worst_solve_residual <= 1.1e-13
        && worst_l2 < 2e-8
        && worst_h1 < 2e-7;
    verdict(
        "cte-001",
        pass,
        &format!(
            "\"detail\":\"roundoff affine residual plus absolute solver forward-error gates across arbitrary cut fractions\",\
             \"rows\":[{}],\"worst_affine_residual\":{worst_affine_residual:.3e},\
             \"worst_solve_residual\":{worst_solve_residual:.3e},\
             \"worst_l2\":{worst_l2:.3e},\"worst_h1\":{worst_h1:.3e}",
            rows.trim_end_matches(',')
        ),
    );
}

fn mms_u(x: f64, y: f64) -> [f64; 2] {
    [
        (PI * x).sin() * (PI * y).sin(),
        (2.0 * PI * x).sin() * (PI * y).sin(),
    ]
}

fn mms_grad(x: f64, y: f64) -> [[f64; 2]; 2] {
    [
        [
            PI * (PI * x).cos() * (PI * y).sin(),
            PI * (PI * x).sin() * (PI * y).cos(),
        ],
        [
            2.0 * PI * (2.0 * PI * x).cos() * (PI * y).sin(),
            PI * (2.0 * PI * x).sin() * (PI * y).cos(),
        ],
    ]
}

fn mms_f(x: f64, y: f64) -> [f64; 2] {
    let mat = material();
    let (lambda, mu) = mat.lame();
    let u = mms_u(x, y);
    let pi2 = PI * PI;
    [
        (lambda + 3.0 * mu) * pi2 * u[0]
            - 2.0 * (lambda + mu) * pi2 * (2.0 * PI * x).cos() * (PI * y).cos(),
        (lambda + 6.0 * mu) * pi2 * u[1] - (lambda + mu) * pi2 * (PI * x).cos() * (PI * y).cos(),
    ]
}

#[test]
fn cte_002_curved_sdf_mms_order_gate() {
    let sdf = Circle {
        center: [0.5, 0.5],
        radius: 0.3,
    };
    let mat = material();
    let mut errors = Vec::new();
    let mut rows = String::new();
    for level in [4u32, 5, 6] {
        let grid = Quadtree::uniform(level);
        let cut = problem(&grid, &sdf, &mat, 0.5);
        let solution = cut.solve(&mms_f, &mms_u).expect("MMS solve");
        let error = cut.l2_h1_error(&solution, &mms_u, &mms_grad);
        let _ = write!(
            rows,
            "{{\"level\":{level},\"l2\":{:.3e},\"h1\":{:.3e},\"iters\":{}}},",
            error.0, error.1, solution.iters
        );
        errors.push(error);
    }
    let l2_orders = [
        (errors[0].0 / errors[1].0).log2(),
        (errors[1].0 / errors[2].0).log2(),
    ];
    let h1_orders = [
        (errors[0].1 / errors[1].1).log2(),
        (errors[1].1 / errors[2].1).log2(),
    ];
    let pass = errors
        .iter()
        .all(|(l2, h1)| l2.is_finite() && h1.is_finite())
        && l2_orders
            .iter()
            .all(|order| order.is_finite() && (1.8..=2.2).contains(order))
        && h1_orders
            .iter()
            .all(|order| order.is_finite() && (0.8..=1.2).contains(order))
        && errors
            .windows(2)
            .all(|pair| pair[1].0 < pair[0].0 && pair[1].1 < pair[0].1)
        && errors[2].0 < 3e-3;
    verdict(
        "cte-002",
        pass,
        &format!(
            "\"detail\":\"G1 vector MMS on curved SDF\",\"rows\":[{}],\
             \"l2_orders\":[{:.2},{:.2}],\"h1_orders\":[{:.2},{:.2}]",
            rows.trim_end_matches(','),
            l2_orders[0],
            l2_orders[1],
            h1_orders[0],
            h1_orders[1]
        ),
    );
}

fn elasticity_condition(epsilon: f64, ghost_gamma: f64) -> f64 {
    // Full-spectrum cyclic Jacobi is cubic in the vector-system dimension.
    // Level 3 retains the same five-decade cut-fraction sweep while keeping
    // this acceptance test bounded in debug CI.
    let grid = Quadtree::uniform(3);
    let h = 1.0 / 8.0;
    let sdf = HalfPlane {
        normal: [1.0, 0.0],
        offset: 0.5 + epsilon * h,
    };
    let mat = material();
    let clamp_all = |_: f64, _: f64| true;
    let cut = CutElasticity {
        clamp: Some(&clamp_all),
        ..problem(&grid, &sdf, &mat, ghost_gamma)
    };
    let operator = cut
        .assemble(&|_, _| [0.0, 0.0], &|_, _| [0.0, 0.0])
        .expect("conditioning operator");
    assert_bit_symmetric(&operator);
    condition_estimate(operator.matrix()).cond
}

fn assert_bit_symmetric(operator: &CutElasticityOperator) {
    let matrix = operator.matrix();
    for row in 0..matrix.nrows() {
        let (columns, values) = matrix.row(row);
        for (&column, &value) in columns.iter().zip(values) {
            assert!(
                value.is_finite(),
                "elasticity matrix contains a non-finite entry at ({row}, {column})"
            );
            assert_eq!(
                value.to_bits(),
                matrix.get(column, row).to_bits(),
                "elasticity matrix lost exact symmetry at ({row}, {column})"
            );
        }
    }
}

#[test]
fn cte_003_ghost_penalty_bounds_degenerate_cut_conditioning() {
    let epsilons = [0.5, 1e-2, 1e-4, 1e-6, 1e-8];
    let mut with_ghost = Vec::new();
    let mut without_ghost = Vec::new();
    let mut rows = String::new();
    for &epsilon in &epsilons {
        with_ghost.push(elasticity_condition(epsilon, 0.5));
        without_ghost.push(elasticity_condition(epsilon, 0.0));
    }
    for (index, epsilon) in epsilons.iter().enumerate() {
        let ghost = json_number(with_ghost[index]);
        let bare = json_number(without_ghost[index]);
        let _ = write!(
            rows,
            "{{\"epsilon\":{epsilon:.1e},\"ghost\":{ghost},\"bare\":{bare}}},"
        );
    }
    let ghost_min = with_ghost.iter().copied().fold(f64::INFINITY, f64::min);
    let ghost_max = with_ghost.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let ghost_ratio = ghost_max / ghost_min;
    let bare_blowup = without_ghost[4] / without_ghost[0];
    let degenerate_improvement = without_ghost[4] / with_ghost[4];
    let pass = with_ghost.iter().all(|value| value.is_finite())
        && without_ghost.iter().all(|value| value.is_finite())
        && ghost_ratio.is_finite()
        && bare_blowup.is_finite()
        && degenerate_improvement.is_finite()
        && ghost_ratio < 100.0
        && bare_blowup > 100.0
        && degenerate_improvement > 100.0;
    verdict(
        "cte-003",
        pass,
        &format!(
            "\"detail\":\"cut-independent Nitsche with ghost-controlled slivers\",\
             \"rows\":[{}],\"ghost_max_over_min\":{},\"bare_blowup\":{},\
             \"degenerate_bare_over_ghost\":{}",
            rows.trim_end_matches(','),
            json_number(ghost_ratio),
            json_number(bare_blowup),
            json_number(degenerate_improvement)
        ),
    );
}

#[test]
fn elasticity_inputs_fail_closed() {
    let grid = Quadtree::uniform(2);
    let sdf = Circle {
        center: [0.5, 0.5],
        radius: 0.3,
    };
    let mat = material();
    let invalid = CutElasticity {
        ghost_gamma: f64::NAN,
        ..problem(&grid, &sdf, &mat, 0.5)
    };
    let error = invalid
        .assemble(&|_, _| [0.0, 0.0], &|_, _| [0.0, 0.0])
        .expect_err("NaN stabilization must refuse");
    assert!(matches!(error, CutFemError::InvalidElasticityInput { .. }));

    let malformed_material = IsotropicElastic {
        youngs: 1.0,
        poisson: 0.75,
        strain_limit: 10.0,
    };
    let malformed = CutElasticity {
        material: &malformed_material,
        ..problem(&grid, &sdf, &mat, 0.5)
    };
    assert!(matches!(
        malformed.assemble(&|_, _| [0.0, 0.0], &|_, _| [0.0, 0.0]),
        Err(CutFemError::InvalidElasticityInput { .. })
    ));

    // Lock both sides of the advertised compressible-regime boundary.  For
    // plane strain, nu = 1/3 gives (lambda + 2*mu)/mu = 4 exactly.
    let boundary_poisson = 1.0f64 / 3.0;
    let boundary_material =
        IsotropicElastic::new(1.0, boundary_poisson, 10.0).expect("boundary material");
    let boundary_problem = CutElasticity {
        material: &boundary_material,
        ..problem(&grid, &sdf, &mat, 0.5)
    };
    boundary_problem
        .assemble(&|_, _| [0.0, 0.0], &|_, _| [0.0, 0.0])
        .expect("stiffness-ratio limit itself remains admissible");

    let above_boundary_poisson = f64::from_bits(boundary_poisson.to_bits() + 1);
    let above_boundary_material =
        IsotropicElastic::new(1.0, above_boundary_poisson, 10.0).expect("adjacent material");
    let above_boundary_problem = CutElasticity {
        material: &above_boundary_material,
        ..problem(&grid, &sdf, &mat, 0.5)
    };
    assert!(matches!(
        above_boundary_problem.assemble(&|_, _| [0.0, 0.0], &|_, _| [0.0, 0.0]),
        Err(CutFemError::InvalidElasticityInput { .. })
    ));

    let near_incompressible =
        IsotropicElastic::new(1.0, 0.49, 10.0).expect("base material accepts nu < 0.5");
    let unsupported_regime = CutElasticity {
        material: &near_incompressible,
        ..problem(&grid, &sdf, &mat, 0.5)
    };
    let error = unsupported_regime
        .assemble(&|_, _| [0.0, 0.0], &|_, _| [0.0, 0.0])
        .expect_err("mu-scaled v1 stabilization must refuse near incompressibility");
    assert!(
        matches!(&error, CutFemError::InvalidElasticityInput { what } if what.contains("near-incompressible")),
        "unexpected refusal: {error}"
    );

    // Assembly-only users and VJP consumers do not pay for irrelevant solver
    // validation; a natural interface likewise does not inspect Nitsche data.
    let assembly_only = CutElasticity {
        traction_free_interface: true,
        nitsche_beta: f64::NAN,
        solver_tol: f64::NAN,
        ..problem(&grid, &sdf, &mat, 0.5)
    };
    assembly_only
        .assemble(&|_, _| [0.0, 0.0], &|_, _| [f64::NAN; 2])
        .expect("unused solve and Dirichlet settings must not poison assembly");
    assert!(matches!(
        assembly_only.solve(&|_, _| [0.0, 0.0], &|_, _| [0.0, 0.0]),
        Err(CutFemError::InvalidElasticityInput { .. })
    ));

    // The current boundary-load contract is fail-closed until certified 1-D
    // clipping exists; never quantize a partially cut edge by sample masking.
    let half_plane = HalfPlane {
        normal: [1.0, 0.0],
        offset: 0.51,
    };
    let traction = |_: f64, _: f64| [1.0, 0.0];
    let clamp_all = |_: f64, _: f64| true;
    let cut_loaded_edge = CutElasticity {
        boundary_traction: Some(&traction),
        traction_free_interface: true,
        clamp: Some(&clamp_all),
        ..problem(&grid, &half_plane, &mat, 0.5)
    };
    assert!(matches!(
        cut_loaded_edge.assemble(&|_, _| [0.0, 0.0], &|_, _| [0.0, 0.0]),
        Err(CutFemError::InvalidElasticityInput { .. })
    ));
}
