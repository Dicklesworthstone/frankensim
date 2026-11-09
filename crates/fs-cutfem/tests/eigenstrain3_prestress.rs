//! Signed CutFEM initial stress, native constraints and geometric Nitsche flux.
use std::cell::Cell;
use std::ops::ControlFlow;

use fs_blake3::ContentHash;
use fs_cutfem::elastic3::adaptive::AdaptiveElasticity3;
use fs_cutfem::elastic3::{
    CutElasticity3, ElasticMaterialField3, ElasticityError3, ElasticityFormulation3,
    ElasticityOptions3, PreparedPrestress3,
};
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_material::elastic::acoustic::ElasticAcousticState;
use fs_material::tensor::{from_mandel, to_mandel};
use fs_solver::op::LinearOp;

const CUT: f64 = 0.73;
const YOUNG: f64 = 30.0;
const SCALE: f64 = 0.6;

struct Slab(usize);
impl CutSdf3 for Slab {
    fn value(&self, p: [f64; 3]) -> f64 {
        p[self.0] - CUT
    }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        Interval::new(lo[self.0], hi[self.0]) - Interval::new(CUT, CUT)
    }
    fn derivative_enclose(&self, _: [f64; 3], _: [f64; 3], axis: HeightAxis) -> Interval {
        let index = match axis {
            HeightAxis::X => 0,
            HeightAxis::Y => 1,
            HeightAxis::Z => 2,
        };
        let value = if self.0 == index { 1.0 } else { 0.0 };
        Interval::new(value, value)
    }
}
fn poll() -> ControlFlow<()> {
    ControlFlow::Continue(())
}
fn domain() -> HexCell {
    HexCell::try_new([0.0; 3], [1.0; 3]).unwrap()
}
fn control<T>(f: impl FnOnce(&mut QuadratureControl3<'_>) -> T) -> T {
    let mut checkpoint = |_| poll();
    let mut q = QuadratureControl3::new(
        QuadratureOptions3 {
            depth: 0,
            ..Default::default()
        },
        &mut checkpoint,
    )
    .unwrap();
    f(&mut q)
}
fn options(formulation: ElasticityFormulation3) -> ElasticityOptions3 {
    ElasticityOptions3 {
        formulation,
        ..Default::default()
    }
}
fn cartesian(
    formulation: ElasticityFormulation3,
    axis: usize,
    counts: [usize; 3],
    embedded: bool,
    poisson: f64,
    limit: f64,
) -> CutElasticity3 {
    control(|q| {
        let law = IsotropicElastic::new(YOUNG, poisson, limit).unwrap();
        let clamp = |p| p == [0.0; 3];
        let mut op = if embedded {
            CutElasticity3::build_with_embedded_dirichlet(
                domain(),
                counts,
                &Slab(axis),
                &law,
                &clamp,
                &|_, _| true,
                options(formulation),
                Default::default(),
                Default::default(),
                q,
            )
            .unwrap()
        } else {
            CutElasticity3::build(
                domain(),
                counts,
                &Slab(axis),
                &law,
                &clamp,
                options(formulation),
                q,
            )
            .unwrap()
        };
        op.set_scales(&vec![SCALE; op.cells()]).unwrap();
        op
    })
}
fn adaptive() -> AdaptiveElasticity3 {
    let tree = Octree3::uniform(1, 3, 256).unwrap();
    let tree = tree
        .refined(&[*tree.leaves().iter().next().unwrap()], poll)
        .unwrap();
    control(|q| {
        AdaptiveElasticity3::build(
            domain(),
            &tree,
            &Slab(2),
            &IsotropicElastic::new(YOUNG, 0.0, 0.1).unwrap(),
            &|p| p == [0.0; 3],
            options(ElasticityFormulation3::Displacement),
            q,
        )
        .unwrap()
    })
}
fn close(a: f64, b: f64, tol: f64) {
    assert!(
        a.is_finite() && b.is_finite() && (a - b).abs() <= tol * a.abs().max(b.abs()).max(1.0),
        "{a:e} != {b:e}"
    );
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn sigma_grad(s: [f64; 6], g: [f64; 3]) -> [f64; 3] {
    [
        s[0] * g[0] + s[3] * g[1] + s[5] * g[2],
        s[3] * g[0] + s[1] * g[1] + s[4] * g[2],
        s[5] * g[0] + s[4] * g[1] + s[2] * g[2],
    ]
}
fn action(g: &PreparedPrestress3, x: &[f64]) -> Vec<f64> {
    let mut out = vec![0.0; g.n()];
    g.apply_into(x, &mut out, poll).unwrap();
    out
}
fn prepare(op: &CutElasticity3, free: [f64; 6], u: &[f64]) -> PreparedPrestress3 {
    let field = op
        .prepare_eigenstrain(&|_, _| Ok(free), 100_000, poll)
        .unwrap();
    op.prepare_eigenstrain_prestress(&field, u, 0.1, 100_000, 100_000, poll)
        .unwrap()
}

#[test]
fn compatible_expansion_vanishes_while_compression_and_tension_keep_their_signs() {
    let expansion = 0.002;
    for formulation in [
        ElasticityFormulation3::Displacement,
        ElasticityFormulation3::BBar,
    ] {
        let op = cartesian(formulation, 2, [1; 3], false, 0.0, 0.1);
        let u: Vec<_> = op
            .nodes()
            .iter()
            .flat_map(|p| [expansion * p[0], 0.0, 0.0])
            .collect();
        let mode: Vec<_> = op.nodes().iter().flat_map(|p| [0.0, p[0], 0.0]).collect();
        let zero = prepare(&op, [expansion, 0.0, 0.0, 0.0, 0.0, 0.0], &u);
        close(zero.bilinear_work(&mode, &mode, poll).unwrap(), 0.0, 1e-12);
        for entry in zero.diagonal(poll).unwrap() {
            close(entry, 0.0, 1e-12);
        }
        let state = vec![0.0; op.n()];
        let compression = prepare(&op, [expansion, 0.0, 0.0, 0.0, 0.0, 0.0], &state);
        let tension = prepare(&op, [-expansion, 0.0, 0.0, 0.0, 0.0, 0.0], &state);
        close(
            compression.bilinear_work(&mode, &mode, poll).unwrap(),
            -SCALE * YOUNG * expansion * CUT,
            1e-11,
        );
        close(
            tension.bilinear_work(&mode, &mode, poll).unwrap(),
            SCALE * YOUNG * expansion * CUT,
            1e-11,
        );
        assert!(
            compression
                .diagonal(poll)
                .unwrap()
                .iter()
                .any(|v| *v < -1e-5)
        );
        assert!(tension.diagonal(poll).unwrap().iter().any(|v| *v > 1e-5));
        assert_eq!(compression.fixed(), op.fixed());
        let output = action(&compression, &mode);
        for (node, fixed) in op.fixed().iter().enumerate() {
            if *fixed {
                assert_eq!(&output[3 * node..3 * node + 3], &[0.0; 3]);
            }
        }
    }
}

#[test]
fn original_full_tensor_stress_matches_affine_geometric_contraction_and_symmetric_actions() {
    let diagonal = [1200.0, 900.0, 1500.0, 800.0, 650.0, 700.0];
    let mut tensor = [[0.0; 6]; 6];
    for i in 0..6 {
        tensor[i][i] = diagonal[i];
    }
    tensor[0][3] = 80.0;
    tensor[3][0] = 80.0;
    tensor[1][2] = 150.0;
    tensor[2][1] = 150.0;
    let identity = ContentHash([93; 32]);
    let material = ElasticAcousticState::try_new_mandel(tensor, 2700.0, identity, 0.1).unwrap();
    let mut op = control(|q| {
        CutElasticity3::build_with_material_field(
            domain(),
            [1; 3],
            &Slab(2),
            ElasticMaterialField3::Uniform(&material),
            &|p| p == [0.0; 3],
            options(ElasticityFormulation3::Displacement),
            q,
        )
        .unwrap()
    });
    op.set_scales(&[SCALE]).unwrap();
    let free = [0.001, -0.0004, 0.0007, 0.0003, -0.0002, 0.0005];
    let g = prepare(&op, free, &vec![0.0; op.n()]);
    let mandel = to_mandel(&free);
    let sigma = from_mandel(&std::array::from_fn(|i| {
        -SCALE * (0..6).map(|j| tensor[i][j] * mandel[j]).sum::<f64>()
    }));
    let left_gradient = [[1.0, 2.0, -0.5], [-0.4, 0.7, 0.2], [0.1, -0.2, 0.3]];
    let right_gradient = [[-0.5, 0.3, 1.0], [0.8, -0.2, 0.4], [0.6, 0.1, -0.7]];
    let affine = |gradient: [[f64; 3]; 3]| -> Vec<f64> {
        op.nodes()
            .iter()
            .flat_map(|p| gradient.map(|row| dot(&row, p)))
            .collect()
    };
    let left = affine(left_gradient);
    let right = affine(right_gradient);
    let expected = CUT
        * (0..3)
            .map(|c| dot(&left_gradient[c], &sigma_grad(sigma, right_gradient[c])))
            .sum::<f64>();
    close(
        g.bilinear_work(&left, &right, poll).unwrap(),
        expected,
        1e-10,
    );
    close(
        g.bilinear_work(&right, &left, poll).unwrap(),
        expected,
        1e-10,
    );
    close(dot(&left, &action(&g, &right)), expected, 1e-10);
    let diagonal = g.diagonal(poll).unwrap();
    for i in 0..g.n() {
        let mut basis = vec![0.0; g.n()];
        basis[i] = 1.0;
        close(action(&g, &basis)[i], diagonal[i], 1e-12);
    }
    assert_eq!(
        op.cell_material_state(0).unwrap().material_state_identity(),
        identity
    );
    assert_eq!(g.scalar_pairs().len(), g.scalar_entries());
}

#[test]
fn spatial_stress_is_correlated_with_actual_gradients_including_bbar_recovery() {
    let expansion = 0.002;
    for formulation in [
        ElasticityFormulation3::Displacement,
        ElasticityFormulation3::BBar,
    ] {
        let poisson = if formulation == ElasticityFormulation3::BBar {
            0.49
        } else {
            0.0
        };
        let op = cartesian(formulation, 2, [1; 3], false, poisson, 0.1);
        let free = op
            .prepare_eigenstrain(
                &|_, p| Ok([expansion * p[2], 0.0, 0.0, 0.0, 0.0, 0.0]),
                100_000,
                poll,
            )
            .unwrap();
        let u = vec![0.0; op.n()];
        let g = op
            .prepare_eigenstrain_prestress(&free, &u, 0.1, 100_000, 100_000, poll)
            .unwrap();
        let direction: Vec<_> = op
            .nodes()
            .iter()
            .flat_map(|p| [0.0, p[0] * p[2], 0.0])
            .collect();
        let response = op
            .eigenstrain_response(&free, &u, 0.1, 100_000, poll)
            .unwrap();
        let mut expected = 0.0;
        let mut mean = [0.0; 6];
        let volume: f64 = response.iter().map(|p| p.weight).sum();
        for point in &response {
            let gradient = [point.position[2], 0.0, point.position[0]];
            expected += point.weight * dot(&gradient, &sigma_grad(point.stress, gradient));
            for i in 0..6 {
                mean[i] += point.weight / volume * point.stress[i];
            }
        }
        let averaged: f64 = response
            .iter()
            .map(|point| {
                let gradient = [point.position[2], 0.0, point.position[0]];
                point.weight * dot(&gradient, &sigma_grad(mean, gradient))
            })
            .sum();
        close(
            g.bilinear_work(&direction, &direction, poll).unwrap(),
            expected,
            1e-10,
        );
        assert!(
            (expected - averaged).abs() > 1e-7,
            "element-mean stress must not replace original stress/gradient correlation"
        );
        if formulation == ElasticityFormulation3::Displacement {
            close(
                expected,
                -SCALE * YOUNG * expansion * CUT.powi(4) / 4.0,
                1e-11,
            );
        }
    }
}

// Independent physical-field integration for sigma_xx on a rectangular slab.
// Only d/dx enters; x is integrated analytically and y/z use exact Gauss rules.
fn physical_axial_work(
    nodes: &[[f64; 3]],
    cells: &[[usize; 8]],
    field: &[f64],
    scales: &[f64],
    sigma: f64,
) -> f64 {
    let gauss = [
        (-0.7745966692414834, 5.0 / 9.0),
        (0.0, 8.0 / 9.0),
        (0.7745966692414834, 5.0 / 9.0),
    ];
    let mut work = 0.0;
    for (cell, scale) in cells.iter().zip(scales) {
        let lo: [f64; 3] = std::array::from_fn(|i| {
            cell.iter()
                .map(|&node| nodes[node][i])
                .fold(f64::INFINITY, f64::min)
        });
        let hi: [f64; 3] = std::array::from_fn(|i| {
            cell.iter()
                .map(|&node| nodes[node][i])
                .fold(f64::NEG_INFINITY, f64::max)
        });
        let top = hi[2].min(CUT);
        if top <= lo[2] {
            continue;
        }
        for (eta, wy) in gauss {
            let y = 0.5 * (lo[1] + hi[1]) + 0.5 * (hi[1] - lo[1]) * eta;
            for (zeta, wz) in gauss {
                let z = 0.5 * (lo[2] + top) + 0.5 * (top - lo[2]) * zeta;
                let mut dx = [0.0; 3];
                for &node in cell {
                    let p = nodes[node];
                    let nx = if p[0] == lo[0] { -1.0 } else { 1.0 } / (hi[0] - lo[0]);
                    let ny = if p[1] == lo[1] { hi[1] - y } else { y - lo[1] } / (hi[1] - lo[1]);
                    let nz = if p[2] == lo[2] { hi[2] - z } else { z - lo[2] } / (hi[2] - lo[2]);
                    for c in 0..3 {
                        dx[c] += nx * ny * nz * field[3 * node + c];
                    }
                }
                work += scale
                    * sigma
                    * dot(&dx, &dx)
                    * (hi[0] - lo[0])
                    * (hi[1] - lo[1])
                    * (top - lo[2])
                    * 0.25
                    * wy
                    * wz;
            }
        }
    }
    work
}

#[test]
fn adaptive_cache_equals_original_physical_reconstruction_and_cartesian_affine_work() {
    let expansion = 0.001;
    let mut op = adaptive();
    assert!(op.physical_nodes().len() > op.nodes().len());
    let scales: Vec<_> = (0..op.cells()).map(|i| 0.4 + 0.01 * i as f64).collect();
    op.set_scales(&scales).unwrap();
    let free = op
        .prepare_eigenstrain(
            &|_, _| Ok([expansion, 0.0, 0.0, 0.0, 0.0, 0.0]),
            100_000,
            poll,
        )
        .unwrap();
    let u = vec![0.0; op.n()];
    let g = op
        .prepare_eigenstrain_prestress(&free, &u, 0.1, 100_000, 100_000, poll)
        .unwrap();
    let direction: Vec<_> = (0..op.n()).map(|i| 0.1 * ((i % 7) as f64 - 3.0)).collect();
    let physical = op.physical_displacements(&direction).unwrap();
    let expected = physical_axial_work(
        op.physical_nodes(),
        &op.cell_nodes(),
        &physical,
        &scales,
        -YOUNG * expansion,
    );
    close(
        g.bilinear_work(&direction, &direction, poll).unwrap(),
        expected,
        1e-10,
    );
    close(dot(&direction, &action(&g, &direction)), expected, 1e-10);
    op.set_scales(&vec![SCALE; op.cells()]).unwrap();
    assert!(op.validate_prestress(&g).is_err());
    let g = op
        .prepare_eigenstrain_prestress(&free, &u, 0.1, 100_000, 100_000, poll)
        .unwrap();
    let direction: Vec<_> = op
        .nodes()
        .iter()
        .flat_map(|p| [0.0, p[0] * p[2], 0.0])
        .collect();
    let adaptive_work = g.bilinear_work(&direction, &direction, poll).unwrap();
    let cart = cartesian(
        ElasticityFormulation3::Displacement,
        2,
        [4; 3],
        false,
        0.0,
        0.1,
    );
    let gc = prepare(
        &cart,
        [expansion, 0.0, 0.0, 0.0, 0.0, 0.0],
        &vec![0.0; cart.n()],
    );
    let dc: Vec<_> = cart
        .nodes()
        .iter()
        .flat_map(|p| [0.0, p[0] * p[2], 0.0])
        .collect();
    close(
        adaptive_work,
        gc.bilinear_work(&dc, &dc, poll).unwrap(),
        1e-10,
    );
    close(
        adaptive_work,
        -SCALE * YOUNG * expansion * CUT.powi(3) / 3.0,
        1e-10,
    );
    assert!(cart.validate_prestress(&g).is_err());
}

#[test]
fn nitsche_geometric_flux_contains_both_symmetric_boundary_terms() {
    let expansion = 0.001;
    for formulation in [
        ElasticityFormulation3::Displacement,
        ElasticityFormulation3::BBar,
    ] {
        let op = cartesian(formulation, 0, [1; 3], true, 0.0, 0.1);
        let queries = Cell::new(0);
        let free = op
            .prepare_eigenstrain(
                &|_, _| {
                    queries.set(queries.get() + 1);
                    Ok([expansion, 0.0, 0.0, 0.0, 0.0, 0.0])
                },
                100_000,
                poll,
            )
            .unwrap();
        let samples = queries.get();
        let g = op
            .prepare_eigenstrain_prestress(&free, &vec![0.0; op.n()], 0.1, 100_000, 100_000, poll)
            .unwrap();
        assert_eq!(queries.get(), samples);
        let left: Vec<_> = op.nodes().iter().flat_map(|p| [0.0, p[0], 0.0]).collect();
        let right: Vec<_> = op
            .nodes()
            .iter()
            .flat_map(|p| [0.0, p[0] * p[1], 0.0])
            .collect();
        let sigma = -SCALE * YOUNG * expansion;
        // Volume gives sigma*CUT; -2 integral v sigma grad(v).n gives -2*sigma*CUT.
        close(
            g.bilinear_work(&left, &left, poll).unwrap(),
            -sigma * CUT,
            1e-10,
        );
        close(
            g.bilinear_work(&left, &right, poll).unwrap(),
            -0.5 * sigma * CUT,
            1e-10,
        );
        close(
            g.bilinear_work(&right, &left, poll).unwrap(),
            -0.5 * sigma * CUT,
            1e-10,
        );
        assert!(
            op.prepare_eigenstrain_prestress(
                &free,
                &vec![0.0; op.n()],
                0.1,
                free.bulk_points(),
                100_000,
                poll,
            )
            .is_err(),
            "surface recovery belongs to the point allowance"
        );
    }
}

#[test]
fn selected_surface_strain_is_admitted_even_when_bulk_strain_is_zero() {
    let op = cartesian(
        ElasticityFormulation3::Displacement,
        0,
        [1; 3],
        true,
        0.0,
        0.001,
    );
    let free = op
        .prepare_eigenstrain(
            &|_, p| Ok([if p[0] > 0.7 { 0.01 } else { 0.0 }, 0.0, 0.0, 0.0, 0.0, 0.0]),
            100_000,
            poll,
        )
        .unwrap();
    let u = vec![0.0; op.n()];
    assert!(
        op.eigenstrain_response(&free, &u, 0.1, 100_000, poll)
            .is_ok()
    );
    assert!(
        op.prepare_eigenstrain_prestress(&free, &u, 0.1, 100_000, 100_000, poll)
            .is_err()
    );
}

#[test]
fn frozen_owner_scales_budgets_and_cancellation_are_enforced() {
    let mut op = cartesian(
        ElasticityFormulation3::Displacement,
        2,
        [1; 3],
        false,
        0.0,
        0.1,
    );
    let free = op
        .prepare_eigenstrain(&|_, _| Ok([0.001, 0.0, 0.0, 0.0, 0.0, 0.0]), 100_000, poll)
        .unwrap();
    let u = vec![0.0; op.n()];
    let g = op
        .prepare_eigenstrain_prestress(&free, &u, 0.1, 100_000, 100_000, poll)
        .unwrap();
    op.validate_prestress(&g).unwrap();
    op.set_scales(&[f64::from_bits(SCALE.to_bits() + 1)])
        .unwrap();
    assert!(
        op.validate_prestress(&g).is_err(),
        "even one scale bit changes the frozen tangent"
    );
    op.set_scales(&[SCALE]).unwrap();
    op.validate_prestress(&g).unwrap();
    let other = cartesian(
        ElasticityFormulation3::Displacement,
        2,
        [1; 3],
        false,
        0.0,
        0.1,
    );
    assert!(other.validate_prestress(&g).is_err());
    assert!(
        other
            .prepare_eigenstrain_prestress(&free, &u, 0.1, 100_000, 100_000, poll)
            .is_err()
    );
    assert!(
        op.prepare_eigenstrain_prestress(&free, &u, 0.1, 0, 100_000, poll)
            .is_err()
    );
    assert!(
        op.prepare_eigenstrain_prestress(&free, &u, 0.1, 100_000, 0, poll)
            .is_err()
    );
    assert!(
        op.prepare_eigenstrain_prestress(&free, &u, 1.0, 100_000, 100_000, poll)
            .is_err()
    );
    assert!(matches!(
        op.prepare_eigenstrain_prestress(&free, &u, 0.1, 100_000, 100_000, || ControlFlow::Break(
            ()
        ),),
        Err(ElasticityError3::Cancelled)
    ));
    let mut output = vec![5.0; g.n()];
    let mut invalid = u.clone();
    invalid[0] = f64::NAN;
    assert!(g.apply_into(&invalid, &mut output, poll).is_err());
    assert!(output.iter().all(|v| *v == 5.0));
    let mut count = 0;
    g.apply_into(&u, &mut output, || {
        count += 1;
        poll()
    })
    .unwrap();
    let reference: Vec<_> = g.scalar_pairs().collect();
    let mut current = 0;
    assert!(matches!(
        g.apply_into(&u, &mut output, || {
            current += 1;
            if current == count {
                ControlFlow::Break(())
            } else {
                poll()
            }
        }),
        Err(ElasticityError3::Cancelled)
    ));
    assert_eq!(g.scalar_pairs().collect::<Vec<_>>(), reference);
    assert!(action(&g, &u).iter().all(|v| *v == 0.0));
}
