//! `frankensim-gmik` vector CutFEM elasticity acceptance battery.
//!
//! - cte-000/000b clamped and unclamped uniform-operator whole-evidence goldens;
//! - cte-001 G0 affine patch law on exactly represented closed polygons;
//! - cte-002 G1 curved-domain MMS convergence;
//! - cte-003 cut-fraction-independent conditioning under ghost penalty.
//! - cte-005 G0/G3 graded componentwise constraint reduction and replay.
//!
//! The feature-required cte-004 VJP gate lives in `elasticity_adjoint.rs`.

use fs_cutfem::{
    BoundaryTraction, Circle, CutElasticity, CutElasticityOperator, CutElasticitySolution,
    CutFemError, CutSdf, DesignBoxEdge, EdgeBand, HalfPlane, NodeKey, Quadtree, condition_estimate,
};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use std::collections::{BTreeMap, BTreeSet};
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

fn hash_bytes(hash: &mut u64, bytes: &[u8]) {
    for byte in bytes {
        *hash ^= u64::from(*byte);
        *hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
    }
}

fn hash_u64(hash: &mut u64, value: u64) {
    hash_bytes(hash, &value.to_le_bytes());
}

fn hash_usize(hash: &mut u64, value: usize) {
    hash_u64(hash, u64::try_from(value).expect("fixture index fits u64"));
}

fn hash_cell(hash: &mut u64, cell: (u32, u32, u32)) {
    hash_u64(hash, u64::from(cell.0));
    hash_u64(hash, u64::from(cell.1));
    hash_u64(hash, u64::from(cell.2));
}

#[allow(clippy::too_many_lines)]
fn uniform_operator_golden_hash(operator: &CutElasticityOperator) -> u64 {
    let mut hash = 0xcbf2_9ce4_8422_2325;
    hash_bytes(&mut hash, b"fs-cutfem/uniform-elasticity-golden/v1\0");
    let matrix = operator.matrix();
    hash_bytes(&mut hash, b"M");
    hash_usize(&mut hash, matrix.nrows());
    hash_usize(&mut hash, matrix.ncols());
    hash_usize(&mut hash, matrix.nnz());
    for row in 0..matrix.nrows() {
        let (columns, values) = matrix.row(row);
        hash_usize(&mut hash, row);
        hash_usize(&mut hash, columns.len());
        for (&column, &value) in columns.iter().zip(values) {
            hash_usize(&mut hash, column);
            hash_u64(&mut hash, value.to_bits());
        }
    }
    hash_bytes(&mut hash, b"R");
    hash_usize(&mut hash, operator.rhs().len());
    for &value in operator.rhs() {
        hash_u64(&mut hash, value.to_bits());
    }
    hash_bytes(&mut hash, b"N");
    hash_usize(&mut hash, operator.node_ids().len());
    for (&(i, j), &id) in operator.node_ids() {
        hash_u64(&mut hash, u64::from(i));
        hash_u64(&mut hash, u64::from(j));
        hash_usize(&mut hash, id);
    }
    hash_bytes(&mut hash, b"C");
    hash_usize(&mut hash, operator.clamped_dofs().len());
    for &clamped in operator.clamped_dofs() {
        hash_bytes(&mut hash, &[u8::from(clamped)]);
    }
    hash_bytes(&mut hash, b"A");
    hash_usize(&mut hash, operator.active_cells().len());
    for &cell in operator.active_cells() {
        hash_cell(&mut hash, cell);
    }
    hash_bytes(&mut hash, b"Q");
    hash_usize(&mut hash, operator.cut_rules().len());
    for (&cell, rule) in operator.cut_rules() {
        hash_cell(&mut hash, cell);
        hash_usize(&mut hash, rule.bulk.len());
        for &(point, weight) in &rule.bulk {
            hash_u64(&mut hash, point[0].to_bits());
            hash_u64(&mut hash, point[1].to_bits());
            hash_u64(&mut hash, weight.to_bits());
        }
        hash_usize(&mut hash, rule.iface.len());
        for &(point, weight, normal) in &rule.iface {
            hash_u64(&mut hash, point[0].to_bits());
            hash_u64(&mut hash, point[1].to_bits());
            hash_u64(&mut hash, weight.to_bits());
            hash_u64(&mut hash, normal[0].to_bits());
            hash_u64(&mut hash, normal[1].to_bits());
        }
    }
    hash_bytes(&mut hash, b"G");
    hash_usize(&mut hash, operator.ghost_faces().len());
    for &(left, right) in operator.ghost_faces() {
        hash_cell(&mut hash, left);
        hash_cell(&mut hash, right);
    }
    hash_bytes(&mut hash, b"D");
    hash_usize(&mut hash, operator.dropped_cut_cells());
    hash
}

#[test]
fn cte_000_uniform_operator_bits_are_frozen() {
    let grid = Quadtree::uniform(1);
    let sdf = HalfPlane {
        normal: [1.0, 0.0],
        offset: 0.75,
    };
    let mat = material();
    let clamp_box = |_: f64, _: f64| true;
    let cut = CutElasticity {
        quad_depth: 0,
        clamp: Some(&clamp_box),
        ..problem(&grid, &sdf, &mat, 0.5)
    };
    let operator = cut
        .assemble(&|_, _| [1.0, -2.0], &|x, y| [0.125 * x, -0.25 * y])
        .expect("uniform golden operator");
    assert_eq!(operator.dof_count(), 18);
    assert_eq!(operator.matrix().nnz(), 20);
    assert_eq!(
        operator
            .clamped_dofs()
            .iter()
            .filter(|&&value| value)
            .count(),
        16
    );
    assert_eq!(operator.active_cells().len(), 4);
    assert_eq!(operator.cut_rules().len(), 2);
    assert_eq!(operator.ghost_faces().len(), 3);
    assert_eq!(operator.dropped_cut_cells(), 0);
    assert_eq!(
        operator.node_ids(),
        &BTreeMap::from([
            ((0, 0), 0),
            ((0, 1), 3),
            ((0, 2), 5),
            ((1, 0), 1),
            ((1, 1), 2),
            ((1, 2), 4),
            ((2, 0), 6),
            ((2, 1), 7),
            ((2, 2), 8),
        ])
    );
    let actual = uniform_operator_golden_hash(&operator);
    assert_eq!(
        actual, 0xeaff_c0cc_edce_3c66,
        "uniform operator golden changed; reviewed replacement = {actual:#018x}"
    );
}

#[test]
fn cte_000b_unclamped_uniform_operator_bits_are_frozen() {
    let grid = Quadtree::uniform(1);
    let sdf = HalfPlane {
        normal: [1.0, 0.0],
        offset: 0.75,
    };
    let mat = material();
    let cut = CutElasticity {
        quad_depth: 0,
        ..problem(&grid, &sdf, &mat, 0.5)
    };
    let operator = cut
        .assemble(&|_, _| [1.0, -2.0], &|x, y| [0.125 * x, -0.25 * y])
        .expect("unclamped uniform golden operator");
    assert_eq!(operator.dof_count(), 18);
    assert_eq!(operator.matrix().nnz(), 240);
    assert!(operator.clamped_dofs().iter().all(|value| !value));
    let row_lengths: Vec<usize> = (0..operator.matrix().nrows())
        .map(|row| operator.matrix().row(row).0.len())
        .collect();
    assert_eq!(
        row_lengths,
        [
            10, 10, 14, 14, 18, 18, 15, 15, 14, 14, 10, 10, 12, 12, 15, 15, 12, 12
        ]
    );
    let actual = uniform_operator_golden_hash(&operator);
    assert_eq!(
        actual, 0x3ec9_48f9_76c7_36b8,
        "unclamped uniform operator golden changed; reviewed replacement = {actual:#018x}"
    );
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

#[allow(clippy::too_many_lines)] // one focused conformance narrative over the support surface
fn assert_support_api(
    grid: &Quadtree,
    operator: &CutElasticityOperator,
    solution: &CutElasticitySolution,
    exact: Displacement,
    exact_gradient: DisplacementGradient,
) {
    assert_eq!(solution.dof_count(), operator.dof_count());
    assert_eq!(solution.active_cells(), operator.active_cells());

    // The operator was assembled independently from the assembly internal to
    // `solve`; deterministic `b^T x` must therefore agree bit-for-bit.
    let independently_evaluated = operator
        .algebraic_compliance(solution.coefficients())
        .expect("independent algebraic compliance");
    assert!(solution.compliance().is_finite());
    assert_eq!(
        solution.compliance().to_bits(),
        independently_evaluated.to_bits(),
        "solution retains the exact deterministic assembled-load dot product"
    );

    assert!(!solution.cut_rules().is_empty(), "fixture has cut cells");
    assert_eq!(
        solution.cut_rules().len(),
        operator.cut_rules().len(),
        "solve retains every assembled cut rule"
    );
    for ((operator_cell, operator_rule), (solution_cell, solution_rule)) in
        operator.cut_rules().iter().zip(solution.cut_rules())
    {
        assert_eq!(operator_cell, solution_cell);
        assert_eq!(operator_rule.bulk.len(), solution_rule.bulk.len());
        assert_eq!(operator_rule.iface.len(), solution_rule.iface.len());
        for ((operator_point, operator_weight), (solution_point, solution_weight)) in
            operator_rule.bulk.iter().zip(&solution_rule.bulk)
        {
            assert!(
                operator_point
                    .iter()
                    .zip(solution_point)
                    .all(|(a, b)| a.to_bits() == b.to_bits())
            );
            assert_eq!(operator_weight.to_bits(), solution_weight.to_bits());
        }
        for (
            (operator_point, operator_weight, operator_normal),
            (solution_point, solution_weight, solution_normal),
        ) in operator_rule.iface.iter().zip(&solution_rule.iface)
        {
            assert!(
                operator_point
                    .iter()
                    .zip(solution_point)
                    .chain(operator_normal.iter().zip(solution_normal))
                    .all(|(a, b)| a.to_bits() == b.to_bits())
            );
            assert_eq!(operator_weight.to_bits(), solution_weight.to_bits());
        }
    }

    assert!(
        !solution.ghost_faces().is_empty(),
        "ghost stabilization is active"
    );
    assert_eq!(solution.ghost_faces(), operator.ghost_faces());
    assert!(
        solution
            .ghost_faces()
            .iter()
            .all(|(left, right)| left < right),
        "every face pair is canonical"
    );
    assert!(
        solution
            .ghost_faces()
            .windows(2)
            .all(|pair| pair[0] < pair[1]),
        "canonical faces are strictly ordered and unique"
    );

    let assert_affine_value = |cell, point| {
        let (value, gradient) = solution
            .value_gradient(grid, cell, point)
            .expect("active-cell field evaluation");
        let expected_value = exact(point[0], point[1]);
        let expected_gradient = exact_gradient(point[0], point[1]);
        for component in 0..2 {
            assert!(
                (value[component] - expected_value[component]).abs() < 5e-7,
                "affine value mismatch on {cell:?} at {point:?}"
            );
            for axis in 0..2 {
                assert!(
                    (gradient[component][axis] - expected_gradient[component][axis]).abs() < 5e-6,
                    "affine gradient mismatch on {cell:?} at {point:?}"
                );
            }
        }
    };

    let (&cut_cell, cut_rule) = solution.cut_rules().iter().next().expect("cut cell");
    let cut_point = cut_rule.bulk.first().expect("kept cut-cell bulk point").0;
    assert_affine_value(cut_cell, cut_point);
    let inside_cell = solution
        .active_cells()
        .iter()
        .copied()
        .find(|cell| !solution.cut_rules().contains_key(cell))
        .expect("certified-interior cell");
    let (inside_lo, inside_hi) = grid.rect(inside_cell);
    assert_affine_value(
        inside_cell,
        [
            f64::midpoint(inside_lo[0], inside_hi[0]),
            f64::midpoint(inside_lo[1], inside_hi[1]),
        ],
    );

    let inactive = grid
        .leaves()
        .find(|cell| solution.active_cells().binary_search(cell).is_err())
        .expect("fixture has an inactive cell");
    let (inactive_lo, inactive_hi) = grid.rect(inactive);
    assert!(matches!(
        solution.value_gradient(
            grid,
            inactive,
            [
                f64::midpoint(inactive_lo[0], inactive_hi[0]),
                f64::midpoint(inactive_lo[1], inactive_hi[1]),
            ],
        ),
        Err(CutFemError::InvalidElasticityInput { .. })
    ));
    assert!(matches!(
        solution.value_gradient(grid, cut_cell, [f64::NAN, cut_point[1]]),
        Err(CutFemError::InvalidElasticityInput { .. })
    ));
}

fn cte_polygons() -> [ConvexPolygon; 3] {
    [
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
    ]
}

fn graded_affine(x: f64, y: f64) -> [f64; 2] {
    [0.01 + 0.02 * x - 0.01 * y, -0.02 + 0.015 * x + 0.025 * y]
}

fn graded_affine_gradient(_: f64, _: f64) -> [[f64; 2]; 2] {
    [[0.02, -0.01], [0.015, 0.025]]
}

fn active_nodes(grid: &Quadtree, operator: &CutElasticityOperator) -> BTreeSet<NodeKey> {
    operator
        .active_cells()
        .iter()
        .flat_map(|&cell| grid.corner_nodes(cell))
        .collect()
}

fn expand_expected_node(
    node: NodeKey,
    constraints: &BTreeMap<NodeKey, [(NodeKey, f64); 2]>,
    node_ids: &BTreeMap<NodeKey, usize>,
    memo: &mut BTreeMap<NodeKey, Vec<(usize, f64)>>,
    stack: &mut BTreeSet<NodeKey>,
) {
    if memo.contains_key(&node) {
        return;
    }
    assert!(stack.insert(node), "fixture constraint cycle at {node:?}");
    let expansion = if let Some(terms) = constraints.get(&node) {
        let mut composed = BTreeMap::<usize, f64>::new();
        for &(child, weight) in terms {
            expand_expected_node(child, constraints, node_ids, memo, stack);
            for &(id, child_weight) in &memo[&child] {
                *composed.entry(id).or_insert(0.0) += weight * child_weight;
            }
        }
        composed.into_iter().collect()
    } else {
        vec![(node_ids[&node], 1.0)]
    };
    stack.remove(&node);
    memo.insert(node, expansion);
}

fn independent_full_domain_load(
    grid: &Quadtree,
    operator: &CutElasticityOperator,
    body: [f64; 2],
    traction: Option<[f64; 2]>,
) -> Vec<f64> {
    let nodes = active_nodes(grid, operator);
    let active: BTreeSet<_> = operator.active_cells().iter().copied().collect();
    let constraints: BTreeMap<_, _> = grid
        .hanging_constraints(&active, &nodes)
        .into_iter()
        .collect();
    let mut expansions = BTreeMap::new();
    for &node in &nodes {
        expand_expected_node(
            node,
            &constraints,
            operator.node_ids(),
            &mut expansions,
            &mut BTreeSet::new(),
        );
    }
    let mut expected = vec![0.0; operator.dof_count()];
    let mut scatter = |node: NodeKey, component: usize, load: f64| {
        for &(id, weight) in &expansions[&node] {
            let dof = 2 * id + component;
            if !operator.clamped_dofs()[dof] {
                expected[dof] += weight * load;
            }
        }
    };
    for &cell in operator.active_cells() {
        let h = grid.cell_h(cell);
        let corners = grid.corner_nodes(cell);
        for node in corners {
            for (component, value) in body.iter().enumerate() {
                scatter(node, component, 0.25 * h * h * value);
            }
        }
        if let Some(traction) = traction {
            let (level, i, j) = cell;
            let nmax = 1u32 << level;
            for (on_boundary, endpoints) in [
                (j == 0, [corners[0], corners[1]]),
                (i + 1 == nmax, [corners[1], corners[2]]),
                (j + 1 == nmax, [corners[2], corners[3]]),
                (i == 0, [corners[3], corners[0]]),
            ] {
                if on_boundary {
                    for node in endpoints {
                        for (component, value) in traction.iter().enumerate() {
                            scatter(node, component, 0.5 * h * value);
                        }
                    }
                }
            }
        }
    }
    expected
}

fn assert_matrix_bits_eq(left: &CutElasticityOperator, right: &CutElasticityOperator) {
    assert_eq!(left.matrix().nrows(), right.matrix().nrows());
    assert_eq!(left.matrix().ncols(), right.matrix().ncols());
    for row in 0..left.matrix().nrows() {
        let (left_columns, left_values) = left.matrix().row(row);
        let (right_columns, right_values) = right.matrix().row(row);
        assert_eq!(left_columns, right_columns);
        assert!(
            left_values
                .iter()
                .zip(right_values)
                .all(|(a, b)| a.to_bits() == b.to_bits()),
            "matrix row {row} moved at the bit level"
        );
    }
}

fn assert_operator_evidence_bits_eq(left: &CutElasticityOperator, right: &CutElasticityOperator) {
    assert_matrix_bits_eq(left, right);
    assert!(
        left.rhs()
            .iter()
            .zip(right.rhs())
            .all(|(a, b)| a.to_bits() == b.to_bits())
    );
    assert_eq!(left.node_ids(), right.node_ids());
    assert_eq!(left.clamped_dofs(), right.clamped_dofs());
    assert_eq!(left.active_cells(), right.active_cells());
    assert_eq!(left.ghost_faces(), right.ghost_faces());
    assert_eq!(left.dropped_cut_cells(), right.dropped_cut_cells());
    assert_eq!(left.cut_rules().len(), right.cut_rules().len());
    for ((left_cell, left_rule), (right_cell, right_rule)) in
        left.cut_rules().iter().zip(right.cut_rules())
    {
        assert_eq!(left_cell, right_cell);
        assert_eq!(left_rule.bulk.len(), right_rule.bulk.len());
        assert_eq!(left_rule.iface.len(), right_rule.iface.len());
        assert!(left_rule.bulk.iter().zip(&right_rule.bulk).all(
            |((left_point, left_weight), (right_point, right_weight))| {
                left_point
                    .iter()
                    .zip(right_point)
                    .all(|(a, b)| a.to_bits() == b.to_bits())
                    && left_weight.to_bits() == right_weight.to_bits()
            }
        ));
        assert!(left_rule.iface.iter().zip(&right_rule.iface).all(
            |(
                (left_point, left_weight, left_normal),
                (right_point, right_weight, right_normal),
            )| {
                left_point
                    .iter()
                    .zip(right_point)
                    .chain(left_normal.iter().zip(right_normal))
                    .all(|(a, b)| a.to_bits() == b.to_bits())
                    && left_weight.to_bits() == right_weight.to_bits()
            }
        ));
    }
}

fn assert_public_expansion_basis(
    operator: &CutElasticityOperator,
    constraints: &[(NodeKey, [(NodeKey, f64); 2])],
) {
    for dof in 0..operator.dof_count() {
        let mut basis = vec![0.0; operator.dof_count()];
        basis[dof] = 1.0;
        let nodal = operator.nodal_values(&basis);
        for (&node, &id) in operator.node_ids() {
            for component in 0..2 {
                let expected: f64 = if dof == 2 * id + component { 1.0 } else { 0.0 };
                assert_eq!(
                    nodal[&node][component].to_bits(),
                    expected.to_bits(),
                    "terminal basis reconstruction moved at {node:?}, component {component}"
                );
            }
        }
        for &(midpoint, endpoints) in constraints {
            for component in 0..2 {
                let expected = endpoints[0].1 * nodal[&endpoints[0].0][component]
                    + endpoints[1].1 * nodal[&endpoints[1].0][component];
                assert_eq!(
                    nodal[&midpoint][component].to_bits(),
                    expected.to_bits(),
                    "terminal basis violates midpoint constraint at {midpoint:?}, component {component}"
                );
            }
        }
    }
}

fn test_q1_gradients(lo: [f64; 2], hi: [f64; 2], point: [f64; 2]) -> [[f64; 2]; 4] {
    let hx = hi[0] - lo[0];
    let hy = hi[1] - lo[1];
    let xi = (point[0] - lo[0]) / hx;
    let eta = (point[1] - lo[1]) / hy;
    [
        [-(1.0 - eta) / hx, -(1.0 - xi) / hy],
        [(1.0 - eta) / hx, -xi / hy],
        [eta / hx, xi / hy],
        [-eta / hx, (1.0 - xi) / hy],
    ]
}

fn test_elastic_traction(
    gradient: [[f64; 2]; 2],
    lambda: f64,
    mu: f64,
    normal: [f64; 2],
) -> [f64; 2] {
    let trace = gradient[0][0] + gradient[1][1];
    let sigma_xx = lambda * trace + 2.0 * mu * gradient[0][0];
    let sigma_yy = lambda * trace + 2.0 * mu * gradient[1][1];
    let sigma_xy = mu * (gradient[0][1] + gradient[1][0]);
    [
        sigma_xx * normal[0] + sigma_xy * normal[1],
        sigma_xy * normal[0] + sigma_yy * normal[1],
    ]
}

fn test_vector_field(
    grid: &Quadtree,
    cell: (u32, u32, u32),
    point: [f64; 2],
    nodal: &BTreeMap<NodeKey, [f64; 2]>,
) -> ([f64; 2], [[f64; 2]; 2]) {
    let (lo, hi) = grid.rect(cell);
    let hx = hi[0] - lo[0];
    let hy = hi[1] - lo[1];
    let xi = (point[0] - lo[0]) / hx;
    let eta = (point[1] - lo[1]) / hy;
    let shape = [
        (1.0 - xi) * (1.0 - eta),
        xi * (1.0 - eta),
        xi * eta,
        (1.0 - xi) * eta,
    ];
    let gradients = test_q1_gradients(lo, hi, point);
    let values = grid.corner_nodes(cell).map(|node| nodal[&node]);
    let mut value = [0.0; 2];
    let mut gradient = [[0.0; 2]; 2];
    for corner in 0..4 {
        for component in 0..2 {
            value[component] += shape[corner] * values[corner][component];
            gradient[component][0] += gradients[corner][0] * values[corner][component];
            gradient[component][1] += gradients[corner][1] * values[corner][component];
        }
    }
    (value, gradient)
}

#[allow(clippy::too_many_arguments)]
fn independent_nitsche_forms(
    grid: &Quadtree,
    operator: &CutElasticityOperator,
    test_nodal: &BTreeMap<NodeKey, [f64; 2]>,
    trial_nodal: &BTreeMap<NodeKey, [f64; 2]>,
    data: &dyn Fn(f64, f64) -> [f64; 2],
    nitsche_beta: f64,
    lambda: f64,
    mu: f64,
) -> (f64, f64) {
    let mut bilinear = 0.0;
    let mut linear = 0.0;
    for (&cell, rules) in operator.cut_rules() {
        let penalty = nitsche_beta * mu / grid.cell_h(cell);
        for &(point, weight, normal) in &rules.iface {
            let (test_value, test_gradient) = test_vector_field(grid, cell, point, test_nodal);
            let (trial_value, trial_gradient) = test_vector_field(grid, cell, point, trial_nodal);
            let test_traction = test_elastic_traction(test_gradient, lambda, mu, normal);
            let trial_traction = test_elastic_traction(trial_gradient, lambda, mu, normal);
            let boundary_data = data(point[0], point[1]);
            bilinear += weight
                * (penalty * (test_value[0] * trial_value[0] + test_value[1] * trial_value[1])
                    - test_traction[0] * trial_value[0]
                    - test_traction[1] * trial_value[1]
                    - test_value[0] * trial_traction[0]
                    - test_value[1] * trial_traction[1]);
            linear += weight
                * (penalty * (test_value[0] * boundary_data[0] + test_value[1] * boundary_data[1])
                    - test_traction[0] * boundary_data[0]
                    - test_traction[1] * boundary_data[1]);
        }
    }
    (bilinear, linear)
}

fn independent_ghost_energy(
    grid: &Quadtree,
    operator: &CutElasticityOperator,
    coefficients: &[f64],
    gamma: f64,
    mu: f64,
) -> f64 {
    let nodal = operator.nodal_values(coefficients);
    let mut energy = 0.0;
    for &(cell_a, cell_b) in operator.ghost_faces() {
        let (lo_a, hi_a) = grid.rect(cell_a);
        let (lo_b, hi_b) = grid.rect(cell_b);
        let h = grid.cell_h(cell_a);
        let axis = usize::from(cell_a.1 == cell_b.1);
        let (t0, t1) = if axis == 0 {
            (lo_a[1], hi_a[1])
        } else {
            (lo_a[0], hi_a[0])
        };
        let normal = if axis == 0 { [1.0, 0.0] } else { [0.0, 1.0] };
        let face_coordinate = if axis == 0 { hi_a[0] } else { hi_a[1] };
        let corners_a = grid.corner_nodes(cell_a);
        let corners_b = grid.corner_nodes(cell_b);
        let gauss = 0.5 / 3.0f64.sqrt();
        let weight = 0.5 * (t1 - t0);
        for t in [0.5 - gauss, 0.5 + gauss] {
            let varying = t0 + t * (t1 - t0);
            let point = if axis == 0 {
                [face_coordinate, varying]
            } else {
                [varying, face_coordinate]
            };
            let gradients_a = test_q1_gradients(lo_a, hi_a, point);
            let gradients_b = test_q1_gradients(lo_b, hi_b, point);
            for component in 0..2 {
                let mut jump = 0.0;
                for corner in 0..4 {
                    let derivative_a =
                        gradients_a[corner][0] * normal[0] + gradients_a[corner][1] * normal[1];
                    let derivative_b =
                        gradients_b[corner][0] * normal[0] + gradients_b[corner][1] * normal[1];
                    jump += derivative_a * nodal[&corners_a[corner]][component]
                        - derivative_b * nodal[&corners_b[corner]][component];
                }
                energy += gamma * mu * h * weight * jump * jump;
            }
        }
    }
    energy
}

#[test]
#[allow(clippy::too_many_lines)]
fn cte_005_graded_componentwise_patch_reconstructs_midpoints_and_replays() {
    let sdf = cte_polygons().into_iter().next().expect("polygon fixture");
    let mut grid = Quadtree::with_room(3, 5);
    grid.refine_toward_interface(&sdf, 5);
    let mat = material();
    let cut = problem(&grid, &sdf, &mat, 0.5);
    let operator = cut
        .assemble(&|_, _| [0.0, 0.0], &graded_affine)
        .expect("graded affine operator");
    assert_bit_symmetric(&operator);

    let nodes = active_nodes(&grid, &operator);
    let active: BTreeSet<_> = operator.active_cells().iter().copied().collect();
    let constraints = grid.hanging_constraints(&active, &nodes);
    assert!(
        !constraints.is_empty(),
        "fixture must contain hanging nodes"
    );
    assert!(
        operator.node_ids().len() < nodes.len(),
        "hanging nodes must not own algebraic terminal blocks"
    );
    assert_public_expansion_basis(&operator, &constraints);

    let ghost_free = CutElasticity {
        ghost_gamma: 0.0,
        ..cut
    }
    .assemble(&|_, _| [0.0, 0.0], &graded_affine)
    .expect("matching ghost-free graded operator");
    assert_eq!(operator.node_ids(), ghost_free.node_ids());
    let probe: Vec<f64> = (0..operator.dof_count())
        .map(|dof| {
            let residue = (dof * 37 + 11) % 101;
            (f64::from(u32::try_from(residue).expect("small residue")) - 50.0) * 0.001
        })
        .collect();
    let ghost_on_apply = operator.apply_vec(&probe);
    let ghost_off_apply = ghost_free.apply_vec(&probe);
    let actual_ghost_energy: f64 = probe
        .iter()
        .zip(ghost_on_apply.iter().zip(ghost_off_apply))
        .map(|(coefficient, (on, off))| coefficient * (on - off))
        .sum();
    let (_, mu) = mat.lame();
    let expected_ghost_energy = independent_ghost_energy(&grid, &operator, &probe, 0.5, mu);
    assert!(expected_ghost_energy > 0.0 && expected_ghost_energy.is_finite());
    assert!(
        (actual_ghost_energy - expected_ghost_energy).abs()
            <= 2e-11 * expected_ghost_energy.abs().max(1.0),
        "equal-level ghost energy mismatch: actual={actual_ghost_energy:e}, expected={expected_ghost_energy:e}"
    );

    let mut exact_coefficients = vec![0.0; operator.dof_count()];
    for (&node, &id) in operator.node_ids() {
        let point = grid.node_pos(node);
        let value = graded_affine(point[0], point[1]);
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
    assert!(
        affine_residual < 3e-12,
        "graded affine residual {affine_residual:e}"
    );

    let nodal = operator.nodal_values(&exact_coefficients);
    for (&node, value) in &nodal {
        let point = grid.node_pos(node);
        let expected = graded_affine(point[0], point[1]);
        for component in 0..2 {
            assert!(
                (value[component] - expected[component]).abs() < 4e-17,
                "affine reconstruction mismatch at {node:?}, component {component}"
            );
        }
    }
    for (midpoint, endpoints) in constraints {
        for component in 0..2 {
            let expected = endpoints[0].1 * nodal[&endpoints[0].0][component]
                + endpoints[1].1 * nodal[&endpoints[1].0][component];
            assert!(
                (nodal[&midpoint][component] - expected).abs() < 4e-17,
                "midpoint reconstruction mismatch at {midpoint:?}, component {component}"
            );
        }
    }

    let replay = cut
        .assemble(&|_, _| [0.0, 0.0], &graded_affine)
        .expect("graded replay operator");
    assert_operator_evidence_bits_eq(&operator, &replay);
    let replay_nodal = replay.nodal_values(&exact_coefficients);
    assert!(
        nodal
            .iter()
            .zip(&replay_nodal)
            .all(|((left_node, left), (right_node, right))| {
                left_node == right_node
                    && left
                        .iter()
                        .zip(right)
                        .all(|(a, b)| a.to_bits() == b.to_bits())
            })
    );

    let solution = cut
        .solve(&|_, _| [0.0, 0.0], &graded_affine)
        .expect("graded affine solve");
    let (l2, h1) = cut.l2_h1_error(&solution, &graded_affine, &graded_affine_gradient);
    assert!(solution.rel_residual <= 1.1e-13);
    assert!(l2 < 3e-8, "graded affine L2 error {l2:e}");
    assert!(h1 < 3e-7, "graded affine H1 error {h1:e}");
}

#[test]
#[allow(clippy::too_many_lines)]
fn graded_nitsche_reduction_matches_reconstructed_field_energy() {
    // Split only the lower-right base cell. The vertical interface x=0.6
    // cuts two fine cells whose left corners include the midpoint hanging on
    // the coarse lower-left cell. Ghost stabilization is disabled because its
    // mixed-level cut-band redesign is deliberately outside v6dn.2.1.
    let mut grid = Quadtree::with_room(1, 2);
    grid.split((1, 1, 0));
    let sdf = HalfPlane {
        normal: [1.0, 0.0],
        offset: 0.6,
    };
    let mat = material();
    let data = |x: f64, y: f64| [0.2 + 0.1 * x, -0.3 + 0.05 * y];
    let zero = |_: f64, _: f64| [0.0, 0.0];
    let nitsche_problem = problem(&grid, &sdf, &mat, 0.0);
    let nitsche = nitsche_problem
        .assemble(&zero, &data)
        .expect("mixed-level Nitsche operator");
    let natural = CutElasticity {
        traction_free_interface: true,
        ..nitsche_problem
    }
    .assemble(&zero, &data)
    .expect("matching natural-interface operator");
    assert_eq!(nitsche.node_ids(), natural.node_ids());

    let nodes = active_nodes(&grid, &nitsche);
    let active: BTreeSet<_> = nitsche.active_cells().iter().copied().collect();
    let constraints = grid.hanging_constraints(&active, &nodes);
    let hanging: BTreeSet<_> = constraints.iter().map(|(node, _)| *node).collect();
    let constrained_cut_cells: BTreeSet<_> = nitsche
        .cut_rules()
        .keys()
        .copied()
        .filter(|&cell| {
            grid.corner_nodes(cell)
                .into_iter()
                .any(|node| hanging.contains(&node))
        })
        .collect();
    assert!(
        !constrained_cut_cells.is_empty(),
        "fixture must put a non-identity hanging transform on Nitsche support"
    );

    let constraint_map: BTreeMap<_, _> = constraints.iter().copied().collect();
    let mut expansions = BTreeMap::new();
    for &node in &nodes {
        expand_expected_node(
            node,
            &constraint_map,
            nitsche.node_ids(),
            &mut expansions,
            &mut BTreeSet::new(),
        );
    }
    assert!(
        constrained_cut_cells.iter().any(|&cell| {
            grid.corner_nodes(cell)
                .into_iter()
                .any(|node| expansions[&node].len() > 1)
        }),
        "constrained Nitsche support must carry a genuinely non-identity expansion"
    );

    let inverse_ids: BTreeMap<_, _> = nitsche
        .node_ids()
        .iter()
        .map(|(&node, &id)| (id, node))
        .collect();
    let mut transformed_basis: Vec<BTreeMap<NodeKey, [f64; 2]>> =
        Vec::with_capacity(nitsche.dof_count());
    let mut naive_basis: Vec<BTreeMap<NodeKey, [f64; 2]>> = Vec::with_capacity(nitsche.dof_count());
    for dof in 0..nitsche.dof_count() {
        let terminal_id = dof / 2;
        let component = dof % 2;
        let transformed = expansions
            .iter()
            .map(|(&node, expansion)| {
                let mut value = [0.0; 2];
                value[component] = expansion
                    .iter()
                    .find_map(|&(id, weight)| (id == terminal_id).then_some(weight))
                    .unwrap_or(0.0);
                (node, value)
            })
            .collect();
        let mut naive: BTreeMap<_, _> = nodes.iter().map(|&node| (node, [0.0; 2])).collect();
        naive
            .get_mut(&inverse_ids[&terminal_id])
            .expect("terminal node belongs to the active nodal set")[component] = 1.0;
        transformed_basis.push(transformed);
        naive_basis.push(naive);
    }

    assert_bit_symmetric(&nitsche);
    assert_bit_symmetric(&natural);
    let (lambda, mu) = mat.lame();
    let mut max_naive_matrix_error = 0.0f64;
    let mut max_naive_rhs_error = 0.0f64;
    for (column, (trial_transformed, trial_naive)) in
        transformed_basis.iter().zip(&naive_basis).enumerate()
    {
        let mut coefficient = vec![0.0; nitsche.dof_count()];
        coefficient[column] = 1.0;
        let with_nitsche = nitsche.apply_vec(&coefficient);
        let without_nitsche = natural.apply_vec(&coefficient);
        for (row, (test_transformed, test_naive)) in
            transformed_basis.iter().zip(&naive_basis).enumerate()
        {
            let (expected, _) = independent_nitsche_forms(
                &grid,
                &nitsche,
                test_transformed,
                trial_transformed,
                &data,
                nitsche_problem.nitsche_beta,
                lambda,
                mu,
            );
            let (naive, _) = independent_nitsche_forms(
                &grid,
                &nitsche,
                test_naive,
                trial_naive,
                &data,
                nitsche_problem.nitsche_beta,
                lambda,
                mu,
            );
            let actual = with_nitsche[row] - without_nitsche[row];
            let scale = actual.abs().max(expected.abs()).max(1.0);
            assert!(
                (actual - expected).abs() <= 2e-12 * scale,
                "reduced Nitsche matrix mismatch at ({row}, {column}): actual={actual:e}, expected={expected:e}"
            );
            max_naive_matrix_error = max_naive_matrix_error.max((expected - naive).abs());
        }
    }
    for (row, (test_transformed, test_naive)) in
        transformed_basis.iter().zip(&naive_basis).enumerate()
    {
        let (_, expected) = independent_nitsche_forms(
            &grid,
            &nitsche,
            test_transformed,
            test_transformed,
            &data,
            nitsche_problem.nitsche_beta,
            lambda,
            mu,
        );
        let (_, naive) = independent_nitsche_forms(
            &grid,
            &nitsche,
            test_naive,
            test_naive,
            &data,
            nitsche_problem.nitsche_beta,
            lambda,
            mu,
        );
        let actual = nitsche.rhs()[row] - natural.rhs()[row];
        let scale = actual.abs().max(expected.abs()).max(1.0);
        assert!(
            (actual - expected).abs() <= 2e-12 * scale,
            "reduced Nitsche RHS mismatch at {row}: actual={actual:e}, expected={expected:e}"
        );
        max_naive_rhs_error = max_naive_rhs_error.max((expected - naive).abs());
    }
    assert!(
        max_naive_matrix_error > 1e-8,
        "fixture did not distinguish T^T K T from a naive terminal submatrix"
    );
    assert!(
        max_naive_rhs_error > 1e-8,
        "fixture did not distinguish T^T f from a naive terminal subvector"
    );
}

#[test]
fn graded_body_reduces_and_outer_traction_conserves_terminal_loads() {
    let mut grid = Quadtree::with_room(1, 3);
    grid.refine_where(3, &|lo, _| lo[0] < 0.5 && lo[1] < 0.5);
    let sdf = HalfPlane {
        normal: [1.0, 0.0],
        offset: 2.0,
    };
    let mat = material();
    let clamp_left = |x: f64, _: f64| x == 0.0;
    let body_value = [1.0, -2.0];
    let traction_value = [0.25, 0.75];
    let traction = |_: f64, _: f64| traction_value;
    let body = |_: f64, _: f64| body_value;
    let zero = |_: f64, _: f64| [0.0, 0.0];
    let base = CutElasticity {
        ghost_gamma: 0.0,
        clamp: Some(&clamp_left),
        boundary_traction: Some(&traction),
        traction_free_interface: true,
        ..problem(&grid, &sdf, &mat, 0.0)
    };
    let combined = base.assemble(&body, &zero).expect("combined graded load");
    let body_only = CutElasticity {
        boundary_traction: None,
        ..base
    }
    .assemble(&body, &zero)
    .expect("graded body load");
    let traction_only = base
        .assemble(&zero, &zero)
        .expect("graded outer traction load");
    assert_matrix_bits_eq(&combined, &body_only);
    assert_matrix_bits_eq(&combined, &traction_only);
    let nodes = active_nodes(&grid, &combined);
    assert!(combined.node_ids().len() < nodes.len());
    let extent = grid.node_extent();
    for &cell in traction_only.active_cells() {
        for node in grid.corner_nodes(cell) {
            if node.0 == 0 || node.0 == extent || node.1 == 0 || node.1 == extent {
                assert!(
                    traction_only.node_ids().contains_key(&node),
                    "outer-box traction node {node:?} must be an unconstrained terminal"
                );
            }
        }
    }
    assert!(body_only.rhs().iter().any(|value| *value != 0.0));
    assert!(traction_only.rhs().iter().any(|value| *value != 0.0));
    let expected_body = independent_full_domain_load(&grid, &body_only, body_value, None);
    let expected_traction =
        independent_full_domain_load(&grid, &traction_only, [0.0; 2], Some(traction_value));
    for (actual, expected) in body_only.rhs().iter().zip(expected_body) {
        assert!((actual - expected).abs() <= 32.0 * f64::EPSILON * expected.abs().max(1.0));
    }
    for (actual, expected) in traction_only.rhs().iter().zip(expected_traction) {
        assert!((actual - expected).abs() <= 32.0 * f64::EPSILON * expected.abs().max(1.0));
    }
    for ((combined, body), traction) in combined
        .rhs()
        .iter()
        .zip(body_only.rhs())
        .zip(traction_only.rhs())
    {
        let scale = combined.abs().max(body.abs()).max(traction.abs()).max(1.0);
        assert!((combined - body - traction).abs() <= 32.0 * f64::EPSILON * scale);
    }
}

#[test]
fn typed_boundary_edge_bands_are_checked() {
    for (start, end) in [
        (f64::NAN, 0.5),
        (0.5, f64::INFINITY),
        (-f64::EPSILON, 0.5),
        (0.5, 1.0 + f64::EPSILON),
        (0.75, 0.25),
    ] {
        assert!(matches!(
            EdgeBand::new(DesignBoxEdge::Right, start, end),
            Err(CutFemError::InvalidElasticityInput { .. })
        ));
    }

    let point = EdgeBand::new(DesignBoxEdge::Top, 0.5, 0.5)
        .expect("measure-zero closed support remains representable");
    assert_eq!(point.edge(), DesignBoxEdge::Top);
    assert_eq!(point.start().to_bits(), 0.5f64.to_bits());
    assert_eq!(point.end().to_bits(), 0.5f64.to_bits());
}

#[test]
fn typed_right_edge_band_certifies_only_its_exact_support() {
    let grid = Quadtree::uniform(2);
    let sdf = HalfPlane {
        normal: [0.0, 1.0],
        offset: 0.5,
    };
    let mat = material();
    let cut = CutElasticity {
        ghost_gamma: 0.0,
        traction_free_interface: true,
        ..problem(&grid, &sdf, &mat, 0.0)
    };
    let traction = |_: f64, _: f64| [0.0, 1.0];
    let zero = |_: f64, _: f64| [0.0, 0.0];
    let disjoint_band =
        EdgeBand::new(DesignBoxEdge::Right, 0.0, 0.25).expect("valid lower-right support");
    let operator = cut
        .assemble_with_boundary_traction(
            &zero,
            &zero,
            BoundaryTraction::EdgeBand {
                support: disjoint_band,
                value: &traction,
            },
        )
        .expect("support disjoint from the right-edge SDF cut must assemble");
    let total_x: f64 = operator.rhs().iter().step_by(2).sum();
    let total_y: f64 = operator.rhs().iter().skip(1).step_by(2).sum();
    assert_eq!(total_x.to_bits(), 0.0f64.to_bits());
    assert!((total_y - (disjoint_band.end() - disjoint_band.start())).abs() <= f64::EPSILON);

    let cut_band =
        EdgeBand::new(DesignBoxEdge::Right, 0.25, 0.75).expect("valid straddling support");
    assert!(matches!(
        cut.assemble_with_boundary_traction(
            &zero,
            &zero,
            BoundaryTraction::EdgeBand {
                support: cut_band,
                value: &traction,
            },
        ),
        Err(CutFemError::InvalidElasticityInput { .. })
    ));

    let endpoint_contact = EdgeBand::new(DesignBoxEdge::Right, 0.5, 0.75)
        .expect("valid support touching the SDF crossing");
    assert!(matches!(
        cut.assemble_with_boundary_traction(
            &zero,
            &zero,
            BoundaryTraction::EdgeBand {
                support: endpoint_contact,
                value: &traction,
            },
        ),
        Err(CutFemError::InvalidElasticityInput { .. })
    ));

    let zero_traction = |_: f64, _: f64| [0.0, 0.0];
    assert!(matches!(
        cut.assemble_with_boundary_traction(
            &zero,
            &zero,
            BoundaryTraction::Uncertified(&zero_traction),
        ),
        Err(CutFemError::InvalidElasticityInput { .. })
    ));
}

#[test]
fn aligned_typed_edge_band_matches_legacy_callback_bits() {
    let grid = Quadtree::uniform(2);
    let sdf = HalfPlane {
        normal: [1.0, 0.0],
        offset: 2.0,
    };
    let mat = material();
    let value = [0.25, -0.75];
    let legacy_traction = |x: f64, y: f64| {
        if x == 1.0 && (0.25..=0.75).contains(&y) {
            value
        } else {
            [0.0, 0.0]
        }
    };
    let supported_traction = |_: f64, _: f64| value;
    let zero = |_: f64, _: f64| [0.0, 0.0];
    let legacy = CutElasticity {
        ghost_gamma: 0.0,
        boundary_traction: Some(&legacy_traction),
        traction_free_interface: true,
        ..problem(&grid, &sdf, &mat, 0.0)
    }
    .assemble(&zero, &zero)
    .expect("legacy aligned callback");
    let typed_problem = CutElasticity {
        ghost_gamma: 0.0,
        traction_free_interface: true,
        ..problem(&grid, &sdf, &mat, 0.0)
    };
    let typed = typed_problem
        .assemble_with_boundary_traction(
            &zero,
            &zero,
            BoundaryTraction::EdgeBand {
                support: EdgeBand::new(DesignBoxEdge::Right, 0.25, 0.75).expect("aligned support"),
                value: &supported_traction,
            },
        )
        .expect("typed aligned callback");
    assert_operator_evidence_bits_eq(&legacy, &typed);

    let occupied = CutElasticity {
        boundary_traction: Some(&legacy_traction),
        ..typed_problem
    };
    assert!(matches!(
        occupied.assemble_with_boundary_traction(
            &zero,
            &zero,
            BoundaryTraction::EdgeBand {
                support: EdgeBand::new(DesignBoxEdge::Right, 0.25, 0.75).expect("aligned support"),
                value: &supported_traction,
            },
        ),
        Err(CutFemError::InvalidElasticityInput { .. })
    ));
    assert!(matches!(
        occupied.solve_with_boundary_traction(
            &zero,
            &zero,
            BoundaryTraction::EdgeBand {
                support: EdgeBand::new(DesignBoxEdge::Right, 0.25, 0.75).expect("aligned support"),
                value: &supported_traction,
            },
        ),
        Err(CutFemError::InvalidElasticityInput { .. })
    ));
}

#[test]
fn cte_001_constant_strain_patch_on_arbitrary_cuts() {
    let grid = Quadtree::uniform(4);
    let mat = material();
    let polygons = cte_polygons();
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
            if polygon_index == 0 && case == "ux" {
                assert_support_api(&grid, &operator, &solution, exact, gradient);
            }
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
    // Level 4 is still visibly pre-asymptotic for this curved cut (its first
    // L2 slope is about 2.33). Gate the theoretical orders on a fixed
    // three-level log-log fit over the asymptotic ladder. Unfitted curved
    // boundaries move relative to the background lattice at each level, so
    // individual adjacent slopes oscillate slightly even when the ladder's
    // convergence trend is stable; retain those adjacent values in the
    // evidence row and require strict error monotonicity.
    for level in [5u32, 6, 7] {
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
    let adjacent_l2_orders = [
        (errors[0].0 / errors[1].0).log2(),
        (errors[1].0 / errors[2].0).log2(),
    ];
    let adjacent_h1_orders = [
        (errors[0].1 / errors[1].1).log2(),
        (errors[1].1 / errors[2].1).log2(),
    ];
    // With equally spaced refinement levels, the least-squares log-log slope
    // is the endpoint slope divided by the two level intervals.
    let l2_fit_order = 0.5 * (errors[0].0 / errors[2].0).log2();
    let h1_fit_order = 0.5 * (errors[0].1 / errors[2].1).log2();
    let pass = errors
        .iter()
        .all(|(l2, h1)| *l2 > 0.0 && l2.is_finite() && *h1 > 0.0 && h1.is_finite())
        && adjacent_l2_orders
            .iter()
            .chain(&adjacent_h1_orders)
            .all(|order| order.is_finite() && *order > 0.0)
        && (1.8..=2.2).contains(&l2_fit_order)
        && (0.8..=1.2).contains(&h1_fit_order)
        && l2_fit_order.is_finite()
        && h1_fit_order.is_finite()
        && errors
            .windows(2)
            .all(|pair| pair[1].0 < pair[0].0 && pair[1].1 < pair[0].1)
        && errors[2].0 < 3e-3;
    verdict(
        "cte-002",
        pass,
        &format!(
            "\"detail\":\"G1 vector MMS on curved SDF\",\"rows\":[{}],\
             \"adjacent_l2_orders\":[{:.3},{:.3}],\
             \"adjacent_h1_orders\":[{:.3},{:.3}],\
             \"l2_fit_order\":{l2_fit_order:.4},\"h1_fit_order\":{h1_fit_order:.4}",
            rows.trim_end_matches(','),
            adjacent_l2_orders[0],
            adjacent_l2_orders[1],
            adjacent_h1_orders[0],
            adjacent_h1_orders[1]
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
        // Conservative full-element trace constant for this vector Q1,
        // stiffness-ratio<=4 material family. It stays fixed across the
        // entire cut-fraction sweep; only the ghost term controls slivers.
        nitsche_beta: 100.0,
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
        // An indefinite bare matrix is reported as +infinity and is valid
        // evidence of lost coercivity. NaN, negative infinity, and malformed
        // finite values remain fail-closed.
        && without_ghost
            .iter()
            .all(|value| !value.is_nan() && *value > 0.0)
        && ghost_ratio.is_finite()
        && !bare_blowup.is_nan()
        && !degenerate_improvement.is_nan()
        && ghost_ratio < 100.0
        && bare_blowup > 100.0
        && degenerate_improvement > 100.0;
    verdict(
        "cte-003",
        pass,
        &format!(
            "\"detail\":\"cut-independent Nitsche with ghost-controlled slivers\",\
             \"nitsche_beta\":100.0,\
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
