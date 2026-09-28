//! G0 analytical stress/VJP laws and G4 bounded observation refusals.
use fs_cutfem::elastic3::adaptive::AdaptiveElasticity3;
use fs_cutfem::elastic3::stress::BulkStressPoint3;
use fs_cutfem::elastic3::{CutElasticity3, ElasticityError3, ElasticityOptions3};
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_solver::op::LinearOp;
use std::ops::ControlFlow;

struct Slab;
impl CutSdf3 for Slab {
    fn value(&self, p: [f64; 3]) -> f64 {
        p[2] - 0.73
    }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        Interval::new(lo[2], hi[2]) - Interval::new(0.73, 0.73)
    }
    fn derivative_enclose(&self, _: [f64; 3], _: [f64; 3], a: HeightAxis) -> Interval {
        let d = if a == HeightAxis::Z { 1.0 } else { 0.0 };
        Interval::new(d, d)
    }
}
fn poll() -> ControlFlow<()> {
    ControlFlow::Continue(())
}
fn material() -> IsotropicElastic {
    IsotropicElastic::new(10.0, 0.25, 1.0).unwrap()
}
fn domain() -> HexCell {
    HexCell::try_new([0.0; 3], [1.0; 3]).unwrap()
}
fn cart() -> CutElasticity3 {
    let mut callback = |_| poll();
    let mut control = QuadratureControl3::new(
        QuadratureOptions3 {
            depth: 0,
            ..Default::default()
        },
        &mut callback,
    )
    .unwrap();
    // One origin clamp admits the observation operator and preserves affine
    // fields through that point; this fixture makes no unique-solve claim.
    CutElasticity3::build(
        domain(),
        [2; 3],
        &Slab,
        &material(),
        &|p| p == [0.0; 3],
        ElasticityOptions3::default(),
        &mut control,
    )
    .unwrap()
}
fn adaptive() -> AdaptiveElasticity3 {
    let tree = Octree3::uniform(1, 3, 128).unwrap();
    let mark = *tree
        .leaves()
        .iter()
        .find(|c| c.index() == [1, 0, 0])
        .unwrap();
    let tree = tree.refined(&[mark], poll).unwrap();
    let mut callback = |_| poll();
    let mut control = QuadratureControl3::new(
        QuadratureOptions3 {
            depth: 0,
            ..Default::default()
        },
        &mut callback,
    )
    .unwrap();
    AdaptiveElasticity3::build(
        domain(),
        &tree,
        &Slab,
        &material(),
        &|p| p == [0.0; 3],
        ElasticityOptions3::default(),
        &mut control,
    )
    .unwrap()
}
fn affine(nodes: &[[f64; 3]]) -> Vec<f64> {
    nodes
        .iter()
        .flat_map(|&[x, y, z]| {
            [
                0.2 * x + 0.3 * y - 0.1 * z,
                -0.4 * x + 0.1 * y + 0.2 * z,
                0.6 * x - 0.2 * y + 0.4 * z,
            ]
        })
        .collect()
}
fn cotangents(points: &[BulkStressPoint3]) -> Vec<[f64; 6]> {
    points
        .iter()
        .map(|p| {
            let [x, y, z] = p.position;
            [
                1.0 + x,
                -0.2 + y,
                0.4 - z,
                0.3 + x * y,
                -0.7 + y * z,
                0.9 + x * z,
            ]
            .map(|v| p.weight * v)
        })
        .collect()
}
fn observed(points: &[BulkStressPoint3], d: &[[f64; 6]], reference: bool) -> f64 {
    points
        .iter()
        .zip(d)
        .map(|(p, d)| {
            let s = if reference {
                &p.reference_stress
            } else {
                &p.stress
            };
            s.iter().zip(d).map(|(a, b)| a * b).sum::<f64>()
        })
        .sum()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * a.abs().max(b.abs()).max(1.0),
        "{a} != {b}"
    );
}

#[test]
fn affine_cut_stress_uses_original_tensor_weights_and_current_density() {
    let mut op = cart();
    let scales: Vec<_> = (0..op.cells()).map(|i| 0.5 + 0.03 * i as f64).collect();
    op.set_scales(&scales).unwrap();
    let u = affine(op.nodes());
    let points = op.bulk_stress(&u, 10_000, poll).unwrap();
    let expected = [4.4, 3.6, 6.0, -0.4, 0.0, 2.0];
    let mut volume = vec![0.0; op.cells()];
    for p in &points {
        assert!(p.weight > 0.0 && p.position[2] < 0.73);
        volume[p.cell] += p.weight;
        for (i, value) in expected.iter().enumerate() {
            close(p.reference_stress[i], *value, 2e-14);
            close(p.stress[i], scales[p.cell] * value, 2e-14);
        }
    }
    assert_eq!(volume, op.volumes());
    close(volume.iter().sum(), 0.73, 1e-9);
    assert_eq!(points, op.bulk_stress(&u, points.len(), poll).unwrap());
    assert_eq!(op.scales(), scales);
}

#[test]
fn hydrostatic_and_rigid_rotation_have_zero_deviatoric_stress() {
    let op = cart();
    let hydro: Vec<_> = op.nodes().iter().flat_map(|p| p.map(|v| 0.2 * v)).collect();
    let rotation: Vec<_> = op
        .nodes()
        .iter()
        .flat_map(|&[x, y, _]| [-y, x, 0.0])
        .collect();
    for p in op.bulk_stress(&hydro, 10_000, poll).unwrap() {
        for i in 0..3 {
            close(p.stress[i], 4.0, 2e-14);
        }
        for i in 3..6 {
            close(p.stress[i], 0.0, 2e-14);
        }
        let s = p.stress;
        let vm2 = 0.5 * ((s[0] - s[1]).powi(2) + (s[1] - s[2]).powi(2) + (s[2] - s[0]).powi(2))
            + 3.0 * (s[3].powi(2) + s[4].powi(2) + s[5].powi(2));
        assert!(vm2 < 1e-26);
    }
    for p in op.bulk_stress(&rotation, 10_000, poll).unwrap() {
        for s in p.stress {
            close(s, 0.0, 2e-14);
        }
    }
}

#[test]
fn physical_and_reference_pullbacks_match_state_and_scale_differences() {
    let mut op = cart();
    let scales: Vec<_> = (0..op.cells()).map(|i| 0.5 + 0.02 * i as f64).collect();
    op.set_scales(&scales).unwrap();
    let u: Vec<_> = (0..op.n()).map(|i| (i % 11) as f64 / 17.0).collect();
    let v: Vec<_> = (0..op.n()).map(|i| (i % 7) as f64 / 13.0 - 0.2).collect();
    let points = op.bulk_stress(&u, 10_000, poll).unwrap();
    let d = cotangents(&points);
    let pb = op.bulk_stress_pullback(&u, &d, poll).unwrap();
    let reference = op.reference_bulk_stress_pullback(&d, poll).unwrap();
    close(
        reference.iter().zip(&u).map(|(a, b)| a * b).sum(),
        observed(&points, &d, true),
        2e-13,
    );
    let h = 1e-5;
    let plus: Vec<_> = u.iter().zip(&v).map(|(u, v)| u + h * v).collect();
    let minus: Vec<_> = u.iter().zip(&v).map(|(u, v)| u - h * v).collect();
    let fd = (observed(&op.bulk_stress(&plus, 10_000, poll).unwrap(), &d, false)
        - observed(&op.bulk_stress(&minus, 10_000, poll).unwrap(), &d, false))
        / (2.0 * h);
    close(
        fd,
        pb.displacement.iter().zip(&v).map(|(a, b)| a * b).sum(),
        2e-9,
    );
    let ds: Vec<_> = (0..op.cells())
        .map(|i| 0.1 + (i % 3) as f64 / 10.0)
        .collect();
    let plus: Vec<_> = scales.iter().zip(&ds).map(|(s, d)| s + h * d).collect();
    op.set_scales(&plus).unwrap();
    let jp = observed(&op.bulk_stress(&u, 10_000, poll).unwrap(), &d, false);
    let minus: Vec<_> = scales.iter().zip(&ds).map(|(s, d)| s - h * d).collect();
    op.set_scales(&minus).unwrap();
    let jm = observed(&op.bulk_stress(&u, 10_000, poll).unwrap(), &d, false);
    close(
        (jp - jm) / (2.0 * h),
        pb.scales.iter().zip(&ds).map(|(a, b)| a * b).sum(),
        2e-9,
    );
    op.set_scales(&scales).unwrap();
    let energy = op.scale_bilinear_forms(&v, &u, poll).unwrap();
    let mut ku = vec![0.0; op.n()];
    op.apply(&u, &mut ku);
    let actual: f64 = v
        .iter()
        .zip(&ku)
        .enumerate()
        .filter(|(i, _)| !op.fixed()[i / 3])
        .map(|(_, (a, b))| a * b)
        .sum();
    close(
        energy.iter().zip(&scales).map(|(a, b)| a * b).sum(),
        actual,
        2e-13,
    );
    for (i, &fixed) in op.fixed().iter().enumerate() {
        if fixed {
            assert_eq!(&pb.displacement[3 * i..3 * i + 3], &[0.0; 3]);
        }
    }
}

#[test]
fn adaptive_stress_preserves_affine_fields_and_exact_hanging_transpose() {
    let mut op = adaptive();
    assert!(op.physical_nodes().len() > op.nodes().len());
    let scales: Vec<_> = (0..op.cells()).map(|i| 0.4 + 0.01 * i as f64).collect();
    op.set_scales(&scales).unwrap();
    let affine = affine(op.nodes());
    for p in op.bulk_stress(&affine, 10_000, poll).unwrap() {
        for (a, b) in p
            .reference_stress
            .iter()
            .zip([4.4, 3.6, 6.0, -0.4, 0.0, 2.0])
        {
            close(*a, b, 4e-14);
        }
    }
    let u: Vec<_> = (0..op.n()).map(|i| (i % 13) as f64 / 19.0).collect();
    let v: Vec<_> = (0..op.n()).map(|i| (i % 5) as f64 / 7.0 - 0.2).collect();
    let points = op.bulk_stress(&u, 10_000, poll).unwrap();
    let d = cotangents(&points);
    let pb = op.bulk_stress_pullback(&u, &d, poll).unwrap();
    let reference = op.reference_bulk_stress_pullback(&d, poll).unwrap();
    close(
        reference.iter().zip(&u).map(|(a, b)| a * b).sum(),
        observed(&points, &d, true),
        4e-13,
    );
    let h = 1e-5;
    let plus: Vec<_> = u.iter().zip(&v).map(|(u, v)| u + h * v).collect();
    let minus: Vec<_> = u.iter().zip(&v).map(|(u, v)| u - h * v).collect();
    let fd = (observed(&op.bulk_stress(&plus, 10_000, poll).unwrap(), &d, false)
        - observed(&op.bulk_stress(&minus, 10_000, poll).unwrap(), &d, false))
        / (2.0 * h);
    close(
        fd,
        pb.displacement.iter().zip(&v).map(|(a, b)| a * b).sum(),
        2e-9,
    );
    let energy = op.scale_bilinear_forms(&v, &u, poll).unwrap();
    let mut ku = vec![0.0; op.n()];
    op.apply(&u, &mut ku);
    let actual: f64 = v
        .iter()
        .zip(&ku)
        .enumerate()
        .filter(|(i, _)| !op.fixed()[i / 3])
        .map(|(_, (a, b))| a * b)
        .sum();
    close(
        energy.iter().zip(&scales).map(|(a, b)| a * b).sum(),
        actual,
        5e-13,
    );
    assert_eq!(points, op.bulk_stress(&u, 10_000, poll).unwrap());
}

#[test]
fn bounded_stress_observation_refuses_invalid_and_interrupted_work_atomically() {
    let op = adaptive();
    let u = affine(op.nodes());
    let points = op.bulk_stress(&u, 10_000, poll).unwrap();
    let d = cotangents(&points);
    let before = op.scales().to_vec();
    assert!(matches!(
        op.bulk_stress(&u, points.len() - 1, poll),
        Err(ElasticityError3::Invalid(
            "bulk stress point allowance exhausted"
        ))
    ));
    assert!(op.bulk_stress(&u[..u.len() - 1], 10_000, poll).is_err());
    let mut bad = u.clone();
    bad[0] = f64::NAN;
    assert!(op.bulk_stress(&bad, 10_000, poll).is_err());
    assert!(
        op.bulk_stress_pullback(&u, &d[..d.len() - 1], poll)
            .is_err()
    );
    let mut bad = d.clone();
    bad[0][0] = f64::NAN;
    assert!(op.reference_bulk_stress_pullback(&bad, poll).is_err());
    let mut count = 0;
    assert!(matches!(
        op.bulk_stress(&u, 10_000, || {
            count += 1;
            if count > 100 {
                ControlFlow::Break(())
            } else {
                poll()
            }
        }),
        Err(ElasticityError3::Cancelled)
    ));
    assert!(matches!(
        op.bulk_stress_pullback(&u, &d, || ControlFlow::Break(())),
        Err(ElasticityError3::Cancelled)
    ));
    assert_eq!(op.scales(), before);
    assert_eq!(points, op.bulk_stress(&u, 10_000, poll).unwrap());
}
