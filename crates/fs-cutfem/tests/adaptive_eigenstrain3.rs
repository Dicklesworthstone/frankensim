//! G1/G3: thermal mechanics and work transposes survive hanging-node elimination.

use std::ops::ControlFlow;

use fs_cutfem::elastic3::adaptive::AdaptiveElasticity3;
use fs_cutfem::elastic3::{ElasticityError3, ElasticityFormulation3, ElasticityOptions3};
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_solver::op::LinearOp;

struct Halfspace(f64);
impl CutSdf3 for Halfspace {
    fn value(&self, p: [f64; 3]) -> f64 {
        p[2] - self.0
    }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        Interval::new(lo[2] - self.0, hi[2] - self.0)
    }
    fn derivative_enclose(&self, _: [f64; 3], _: [f64; 3], a: HeightAxis) -> Interval {
        let d = if a == HeightAxis::Z { 1.0 } else { 0.0 };
        Interval::new(d, d)
    }
}
fn poll() -> ControlFlow<()> {
    ControlFlow::Continue(())
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn close(a: f64, b: f64, relative: f64, absolute: f64) {
    assert!(
        (a - b).abs() <= absolute + relative * a.abs().max(b.abs()),
        "{a} != {b}"
    );
}
fn build(formulation: ElasticityFormulation3) -> AdaptiveElasticity3 {
    let tree = Octree3::uniform(1, 4, 1000).unwrap();
    let mark = *tree
        .leaves()
        .iter()
        .find(|cell| cell.index() == [0, 0, 1])
        .unwrap();
    let tree = tree.refined(&[mark], poll).unwrap();
    let mut callback = |_| poll();
    let mut control =
        QuadratureControl3::new(QuadratureOptions3::default(), &mut callback).unwrap();
    AdaptiveElasticity3::build(
        HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(),
        &tree,
        &Halfspace(0.73),
        &IsotropicElastic::new(3.0, 0.0, 0.1).unwrap(),
        &|p| p[0] == 0.0,
        ElasticityOptions3 {
            formulation,
            ..Default::default()
        },
        &mut control,
    )
    .unwrap()
}

#[test]
fn free_expansion_on_mixed_octree_is_stress_free_with_original_material_scales() {
    for formulation in [
        ElasticityFormulation3::Displacement,
        ElasticityFormulation3::BBar,
    ] {
        let mut op = build(formulation);
        assert!(op.physical_nodes().len() > op.nodes().len());
        let base: Vec<_> = (0..op.cells()).map(|i| 0.7 + 0.02 * i as f64).collect();
        op.set_scales(&base).unwrap();
        let alpha = 0.002;
        let prepared = op
            .prepare_eigenstrain(&|_, _| Ok([alpha, 0.0, 0.0, 0.0, 0.0, 0.0]), 100_000, poll)
            .unwrap();
        let rhs = op.eigenstrain_load(&prepared, poll).unwrap();
        let solution = op
            .solve_controlled(&rhs, 1e-11, 30_000, 32, |_| poll())
            .unwrap();
        let full = op.physical_displacements(solution.coefficients()).unwrap();
        for (position, displacement) in op.physical_nodes().iter().zip(full.chunks_exact(3)) {
            close(displacement[0], alpha * position[0], 2e-7, 2e-11);
            assert!(displacement[1].abs() < 2e-11 && displacement[2].abs() < 2e-11);
        }
        for point in op
            .eigenstrain_response(&prepared, solution.coefficients(), 0.05, 100_000, poll)
            .unwrap()
        {
            assert!(point.stress.iter().all(|s| s.abs() < 2e-10));
            assert!(point.energy_density < 1e-19);
        }
        let dual: Vec<_> = (0..op.n()).map(|i| 0.1 * ((i % 9) as f64 - 4.0)).collect();
        let pullback = op
            .eigenstrain_load_scale_pullback(&prepared, &dual, poll)
            .unwrap();
        close(dot(&rhs, &dual), dot(&base, &pullback), 2e-13, 1e-14);
    }
}

#[test]
fn adaptive_thermal_load_and_stress_derivatives_use_the_same_master_transpose() {
    let mut op = build(ElasticityFormulation3::BBar);
    let prepared = op
        .prepare_eigenstrain(
            &|cell, p| {
                Ok([
                    0.001 * (1.0 + p[0]),
                    -0.0002 * p[2],
                    0.0001 * cell as f64,
                    0.0003 * p[1],
                    -0.0001 * p[0],
                    0.0002 * p[2],
                ])
            },
            100_000,
            poll,
        )
        .unwrap();
    let base: Vec<_> = (0..op.cells()).map(|i| 0.7 + 0.02 * i as f64).collect();
    op.set_scales(&base).unwrap();
    let u: Vec<_> = op
        .nodes()
        .iter()
        .flat_map(|&[x, y, z]| [0.003 * x + 0.001 * x * y, 0.002 * x * z, -0.001 * x * y])
        .collect();
    let points = op
        .eigenstrain_response(&prepared, &u, 0.05, 100_000, poll)
        .unwrap();
    let derivatives: Vec<_> = points
        .iter()
        .map(|p| [1.0 + p.position[0], -0.2, 0.4, 0.3, -0.7, 0.9].map(|v| p.weight * v))
        .collect();
    let pb = op
        .eigenstrain_stress_pullback(&prepared, &u, 0.05, &derivatives, poll)
        .unwrap();
    let observed = |op: &AdaptiveElasticity3, u: &[f64]| -> f64 {
        op.eigenstrain_response(&prepared, u, 0.05, points.len(), poll)
            .unwrap()
            .iter()
            .zip(&derivatives)
            .map(|(point, d)| dot(&point.stress, d))
            .sum()
    };
    let direction: Vec<_> = (0..op.n())
        .map(|i| 0.0001 * ((i % 5) as f64 - 2.0))
        .collect();
    let h = 1e-4;
    let plus: Vec<_> = u.iter().zip(&direction).map(|(u, d)| u + h * d).collect();
    let minus: Vec<_> = u.iter().zip(&direction).map(|(u, d)| u - h * d).collect();
    close(
        (observed(&op, &plus) - observed(&op, &minus)) / (2.0 * h),
        dot(&pb.displacement, &direction),
        2e-8,
        1e-12,
    );
    let ds: Vec<_> = (0..op.cells()).map(|i| 0.1 * (1 + i % 3) as f64).collect();
    let load_pb = op
        .eigenstrain_load_scale_pullback(&prepared, &direction, poll)
        .unwrap();
    op.set_scales(
        &base
            .iter()
            .zip(&ds)
            .map(|(s, d)| s + h * d)
            .collect::<Vec<_>>(),
    )
    .unwrap();
    let jp = observed(&op, &u);
    let lp = dot(&op.eigenstrain_load(&prepared, poll).unwrap(), &direction);
    op.set_scales(
        &base
            .iter()
            .zip(&ds)
            .map(|(s, d)| s - h * d)
            .collect::<Vec<_>>(),
    )
    .unwrap();
    let jm = observed(&op, &u);
    let lm = dot(&op.eigenstrain_load(&prepared, poll).unwrap(), &direction);
    close((jp - jm) / (2.0 * h), dot(&pb.scales, &ds), 2e-8, 1e-12);
    close((lp - lm) / (2.0 * h), dot(&load_pb, &ds), 2e-8, 1e-14);
    op.set_scales(&base).unwrap();
    let other = build(ElasticityFormulation3::BBar);
    assert!(other.validate_eigenstrain(&prepared).is_err());
    assert!(
        op.eigenstrain_load_scale_pullback(&prepared, &[0.0], poll)
            .is_err()
    );
    assert!(matches!(
        op.eigenstrain_load(&prepared, || ControlFlow::Break(())),
        Err(ElasticityError3::Cancelled)
    ));
    assert!(matches!(
        op.eigenstrain_response(&prepared, &u, 0.05, points.len(), || ControlFlow::Break(())),
        Err(ElasticityError3::Cancelled)
    ));
    assert_eq!(
        points,
        op.eigenstrain_response(&prepared, &u, 0.05, points.len(), poll)
            .unwrap()
    );
    assert_eq!(op.scales(), base);
}
