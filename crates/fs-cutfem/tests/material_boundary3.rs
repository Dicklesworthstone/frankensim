//! G0/G1/G3: source-tensor stresses, transposes, embedded motion and reactions.
use std::ops::ControlFlow;

use fs_blake3::ContentHash;
use fs_cutfem::elastic3::stress::BulkStressPoint3;
use fs_cutfem::elastic3::{CutElasticity3, ElasticMaterialField3, ElasticityError3};
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_ivl::Interval;
use fs_material::elastic::acoustic::ElasticAcousticState;
use fs_solver::LinearOp;

struct Slab;
impl CutSdf3 for Slab {
    fn value(&self, p: [f64; 3]) -> f64 {
        (p[0] - 0.17) * (p[0] - 0.83)
    }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        let x = Interval::new(lo[0], hi[0]);
        (x - Interval::new(0.17, 0.17)) * (x - Interval::new(0.83, 0.83))
    }
    fn derivative_enclose(&self, lo: [f64; 3], hi: [f64; 3], axis: HeightAxis) -> Interval {
        if axis == HeightAxis::X {
            Interval::new(2.0, 2.0) * Interval::new(lo[0], hi[0]) - Interval::new(1.0, 1.0)
        } else {
            Interval::new(0.0, 0.0)
        }
    }
}
fn poll() -> ControlFlow<()> {
    ControlFlow::Continue(())
}
fn domain() -> HexCell {
    HexCell::try_new([0.0; 3], [1.0; 3]).unwrap()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * a.abs().max(b.abs()).max(1.0),
        "{a} != {b}"
    );
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn state(c: [[f64; 6]; 6], ceiling: f64) -> ElasticAcousticState {
    ElasticAcousticState::try_new_mandel(c, 3.0, ContentHash([71; 32]), ceiling).unwrap()
}
fn observed(points: &[BulkStressPoint3], d: &[[f64; 6]], reference: bool) -> f64 {
    points
        .iter()
        .zip(d)
        .map(|(p, d)| {
            dot(
                if reference {
                    &p.reference_stress
                } else {
                    &p.stress
                },
                d,
            )
        })
        .sum()
}

#[test]
fn coupled_tensor_stress_and_vjps_share_the_source_law_and_strain_domain() {
    // Strict diagonal dominance gives an SPD law with normal/shear coupling.
    let material = state(
        [
            [11.0, 2.0, 1.0, 0.4, 0.2, -0.3],
            [2.0, 9.0, 1.5, -0.2, 0.3, 0.1],
            [1.0, 1.5, 8.0, 0.1, -0.3, 0.2],
            [0.4, -0.2, 0.1, 5.0, 0.25, -0.15],
            [0.2, 0.3, -0.3, 0.25, 4.0, 0.35],
            [-0.3, 0.1, 0.2, -0.15, 0.35, 6.0],
        ],
        1.0,
    );
    let mut callback = |_| poll();
    let mut control =
        QuadratureControl3::new(QuadratureOptions3::default(), &mut callback).unwrap();
    // One origin clamp preserves affine fields; this test makes no solve claim.
    let mut op = CutElasticity3::build_with_material_field_and_surface(
        domain(),
        [2; 3],
        &Slab,
        ElasticMaterialField3::Uniform(&material),
        &|p| p == [0.0; 3],
        Default::default(),
        Default::default(),
        &mut control,
    )
    .unwrap();
    let pressure = op.pressure_load(&|_| 0.3, poll).unwrap().rhs;
    let scales: Vec<_> = (0..op.cells()).map(|i| 0.5 + 0.02 * i as f64).collect();
    op.set_scales(&scales).unwrap();
    assert_eq!(pressure, op.pressure_load(&|_| 0.3, poll).unwrap().rhs);
    let u: Vec<_> = op
        .nodes()
        .iter()
        .flat_map(|&[x, y, z]| {
            [
                0.2 * x + 0.3 * y - 0.1 * z,
                -0.4 * x + 0.1 * y + 0.2 * z,
                0.6 * x - 0.2 * y + 0.4 * z,
            ]
        })
        .collect();
    let q = core::f64::consts::FRAC_1_SQRT_2;
    let expected = [
        2.8 - 0.19 * q,
        1.9 + 0.07 * q,
        3.55 + 0.09 * q,
        0.1 * q - 0.2875,
        -0.05 * q + 0.075,
        0.03 * q + 1.5075,
    ];
    let points = op.bulk_stress(&u, 100_000, poll).unwrap();
    for p in &points {
        for k in 0..6 {
            close(p.reference_stress[k], expected[k], 2e-13);
            close(p.stress[k], scales[p.cell] * expected[k], 2e-13);
        }
    }
    // The source domain applies after homogeneous coefficient masking.
    let mut masked = u.clone();
    for (node, fixed) in op.fixed().iter().enumerate() {
        if *fixed {
            masked[3 * node..3 * node + 3].fill(1e100);
        }
    }
    assert_eq!(points, op.bulk_stress(&masked, points.len(), poll).unwrap());
    let d: Vec<_> = points
        .iter()
        .map(|p| [1.0 + p.position[0], -0.2, 0.4, 0.3, -0.7, 0.9].map(|v| p.weight * v))
        .collect();
    let pb = op.bulk_stress_pullback(&u, &d, poll).unwrap();
    let reference = op.reference_bulk_stress_pullback(&d, poll).unwrap();
    close(dot(&reference, &u), observed(&points, &d, true), 2e-12);
    let direction: Vec<_> = (0..op.n()).map(|i| 0.01 * ((i % 7) as f64 - 3.0)).collect();
    let h = 1e-5;
    let plus: Vec<_> = u.iter().zip(&direction).map(|(u, d)| u + h * d).collect();
    let minus: Vec<_> = u.iter().zip(&direction).map(|(u, d)| u - h * d).collect();
    let fd = (observed(
        &op.bulk_stress(&plus, points.len(), poll).unwrap(),
        &d,
        false,
    ) - observed(
        &op.bulk_stress(&minus, points.len(), poll).unwrap(),
        &d,
        false,
    )) / (2.0 * h);
    close(fd, dot(&pb.displacement, &direction), 2e-9);
    let ds: Vec<_> = (0..op.cells()).map(|i| 0.02 * (1 + i % 3) as f64).collect();
    op.set_scales(
        &scales
            .iter()
            .zip(&ds)
            .map(|(s, d)| s + h * d)
            .collect::<Vec<_>>(),
    )
    .unwrap();
    let jp = observed(&op.bulk_stress(&u, points.len(), poll).unwrap(), &d, false);
    op.set_scales(
        &scales
            .iter()
            .zip(&ds)
            .map(|(s, d)| s - h * d)
            .collect::<Vec<_>>(),
    )
    .unwrap();
    let jm = observed(&op.bulk_stress(&u, points.len(), poll).unwrap(), &d, false);
    close((jp - jm) / (2.0 * h), dot(&pb.scales, &ds), 2e-9);
    op.set_scales(&scales).unwrap();
    let outside: Vec<_> = u.iter().map(|v| 20.0 * v).collect();
    assert!(matches!(
        op.bulk_stress(&outside, points.len(), poll),
        Err(ElasticityError3::Invalid(_))
    ));
    assert!(matches!(
        op.bulk_stress_pullback(&outside, &d, poll),
        Err(ElasticityError3::Invalid(_))
    ));
    assert!(matches!(
        op.bulk_stress(&u, points.len(), || ControlFlow::Break(())),
        Err(ElasticityError3::Cancelled)
    ));
    assert_eq!(op.scales(), scales);
    assert_eq!(points, op.bulk_stress(&u, points.len(), poll).unwrap());
}

#[test]
fn anisotropic_embedded_motion_reaction_and_lifting_have_matching_derivatives() {
    let diagonal = [10.0, 7.0, 6.0, 4.0, 3.0, 8.0];
    let material = state(
        std::array::from_fn(|i| std::array::from_fn(|j| if i == j { diagonal[i] } else { 0.0 })),
        0.2,
    );
    let mut callback = |_| poll();
    let mut control =
        QuadratureControl3::new(QuadratureOptions3::default(), &mut callback).unwrap();
    let mut op = CutElasticity3::build_with_material_field_and_embedded_dirichlet(
        domain(),
        [2; 3],
        &Slab,
        ElasticMaterialField3::Uniform(&material),
        &|_| false,
        &|_, _| true,
        Default::default(),
        Default::default(),
        Default::default(),
        &mut control,
    )
    .unwrap();
    let scales = vec![0.6; op.cells()];
    op.set_scales(&scales).unwrap();
    let motion = |p: [f64; 3], _: [f64; 3]| [0.02 * (p[0] - 0.17), 0.0, 0.0];
    let right = |_: [f64; 3], n: [f64; 3]| {
        if n[0] > 0.0 {
            [1.0, 0.0, 0.0]
        } else {
            [0.0; 3]
        }
    };
    let load = op.prescribed_displacement_load(&motion, poll).unwrap();
    let solved = op
        .solve_controlled(&load, 1e-11, 20_000, 32, |_| poll())
        .unwrap();
    let u = solved.coefficients();
    for (p, value) in op.nodes().iter().zip(u.chunks_exact(3)) {
        close(value[0], 0.02 * (p[0] - 0.17), 2e-8);
        close(value[1], 0.0, 2e-8);
        close(value[2], 0.0, 2e-8);
    }
    let reaction = op
        .embedded_reaction(u, Some(&motion), &right, poll)
        .unwrap();
    close(reaction.value, 0.6 * 10.0 * 0.02, 2e-8);
    for p in op.bulk_stress(u, 100_000, poll).unwrap() {
        close(p.stress[0], 0.12, 2e-8);
        for s in &p.stress[1..] {
            close(*s, 0.0, 2e-8);
        }
    }
    let direction: Vec<_> = (0..op.n()).map(|i| 0.01 * ((i % 7) as f64 - 3.0)).collect();
    let h = 1e-5;
    let plus: Vec<_> = u.iter().zip(&direction).map(|(u, d)| u + h * d).collect();
    let minus: Vec<_> = u.iter().zip(&direction).map(|(u, d)| u - h * d).collect();
    let fd = (op
        .embedded_reaction(&plus, Some(&motion), &right, poll)
        .unwrap()
        .value
        - op.embedded_reaction(&minus, Some(&motion), &right, poll)
            .unwrap()
            .value)
        / (2.0 * h);
    close(fd, dot(&reaction.displacement_gradient, &direction), 2e-9);
    let lifting = op
        .prescribed_displacement_scale_work(&motion, &direction, poll)
        .unwrap();
    let ds: Vec<_> = (0..op.cells()).map(|i| 0.02 * (1 + i % 3) as f64).collect();
    op.set_scales(
        &scales
            .iter()
            .zip(&ds)
            .map(|(s, d)| s + h * d)
            .collect::<Vec<_>>(),
    )
    .unwrap();
    let rp = op
        .embedded_reaction(u, Some(&motion), &right, poll)
        .unwrap()
        .value;
    let bp = op.prescribed_displacement_load(&motion, poll).unwrap();
    op.set_scales(
        &scales
            .iter()
            .zip(&ds)
            .map(|(s, d)| s - h * d)
            .collect::<Vec<_>>(),
    )
    .unwrap();
    let rm = op
        .embedded_reaction(u, Some(&motion), &right, poll)
        .unwrap()
        .value;
    let bm = op.prescribed_displacement_load(&motion, poll).unwrap();
    close(
        (rp - rm) / (2.0 * h),
        dot(&reaction.scale_gradient, &ds),
        2e-9,
    );
    close(
        (dot(&bp, &direction) - dot(&bm, &direction)) / (2.0 * h),
        dot(&lifting, &ds),
        2e-9,
    );
    op.set_scales(&scales).unwrap();
    assert_eq!(
        load,
        op.prescribed_displacement_load(&motion, poll).unwrap()
    );
}
