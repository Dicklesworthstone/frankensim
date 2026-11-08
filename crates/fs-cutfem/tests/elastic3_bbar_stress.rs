//! B-bar stress, its exact pullbacks, and embedded-support consistency.
use std::ops::ControlFlow;

use fs_blake3::ContentHash;
use fs_cutfem::elastic3::adaptive::AdaptiveElasticity3;
use fs_cutfem::elastic3::stress::BulkStressPoint3;
use fs_cutfem::elastic3::{
    CutElasticity3, ElasticMaterialField3, ElasticityError3, ElasticityFormulation3,
    ElasticityOptions3,
};
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_material::elastic::acoustic::ElasticAcousticState;
use fs_solver::op::LinearOp;

struct HalfSpace;
impl CutSdf3 for HalfSpace {
    fn value(&self, p: [f64; 3]) -> f64 {
        p[2] - 0.73
    }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        Interval::new(lo[2], hi[2]) - Interval::new(0.73, 0.73)
    }
    fn derivative_enclose(&self, _: [f64; 3], _: [f64; 3], axis: HeightAxis) -> Interval {
        let d = if axis == HeightAxis::Z { 1.0 } else { 0.0 };
        Interval::new(d, d)
    }
}

struct Strip;
impl CutSdf3 for Strip {
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
fn options() -> ElasticityOptions3 {
    ElasticityOptions3 {
        formulation: ElasticityFormulation3::BBar,
        ghost_gamma: 0.0,
        ..Default::default()
    }
}
fn source(limit: f64) -> ElasticAcousticState {
    ElasticAcousticState::from_bulk_and_shear(1000.0, 100_000.0, 3.0, ContentHash([73; 32]), limit)
        .unwrap()
}
fn cartesian(count: usize, limit: f64, formulation: ElasticityFormulation3) -> CutElasticity3 {
    let mut checkpoint = |_| poll();
    let mut control = QuadratureControl3::new(
        QuadratureOptions3 {
            depth: 0,
            ..Default::default()
        },
        &mut checkpoint,
    )
    .unwrap();
    CutElasticity3::build_with_material_field(
        domain(),
        [count; 3],
        &HalfSpace,
        ElasticMaterialField3::Uniform(&source(limit)),
        &|p| p == [0.0; 3],
        ElasticityOptions3 {
            formulation,
            ..options()
        },
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
    let mut checkpoint = |_| poll();
    let mut control = QuadratureControl3::new(
        QuadratureOptions3 {
            depth: 0,
            ..Default::default()
        },
        &mut checkpoint,
    )
    .unwrap();
    AdaptiveElasticity3::build_with_material_field(
        domain(),
        &tree,
        &HalfSpace,
        ElasticMaterialField3::Uniform(&source(0.1)),
        &|p| p == [0.0; 3],
        options(),
        &mut control,
    )
    .unwrap()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(a.is_finite() && b.is_finite());
    assert!(
        (a - b).abs() <= tolerance * a.abs().max(b.abs()).max(1.0),
        "{a:e} != {b:e}"
    );
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn observed(points: &[BulkStressPoint3], dual: &[[f64; 6]], reference: bool) -> f64 {
    points
        .iter()
        .zip(dual)
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
fn bilinear_displacement_retains_local_deviator_and_cell_average_pressure() {
    let op = cartesian(1, 0.1, ElasticityFormulation3::BBar);
    let amplitude = 0.02;
    let u: Vec<_> = op
        .nodes()
        .iter()
        .flat_map(|&[x, y, _]| [amplitude * x * y, 0.0, 0.0])
        .collect();
    let points = op.bulk_stress(&u, 1000, poll).unwrap();
    let (bulk, shear) = (100_000.0, 3.0);
    let mut work = 0.0;
    for p in &points {
        let [x, y, _] = p.position;
        let pressure = bulk * amplitude * 0.5;
        let expected = [
            pressure + 4.0 * shear * amplitude * y / 3.0,
            pressure - 2.0 * shear * amplitude * y / 3.0,
            pressure - 2.0 * shear * amplitude * y / 3.0,
            shear * amplitude * x,
            0.0,
            0.0,
        ];
        for k in 0..6 {
            close(p.stress[k], expected[k], 2e-12);
        }
        close(p.stress[..3].iter().sum::<f64>() / 3.0, pressure, 2e-12);
        // Independent virtual work with the ACTUAL displacement gradient.
        // Its cell-average divergence is the same as the projected divergence.
        work += p.weight * amplitude * (y * p.stress[0] + x * p.stress[3]);
    }
    let mut ku = vec![0.0; op.n()];
    op.apply(&u, &mut ku);
    close(dot(&u, &ku), work, 2e-12);

    let observed = op.bulk_strain_stress(&u, points.len(), poll).unwrap();
    let mut storage = 0.0;
    for (point, stress) in observed.iter().zip(&points) {
        assert_eq!(&point.response, stress);
        let [x, y, _] = stress.position;
        let strain = [amplitude * y, 0.0, 0.0, 0.5 * amplitude * x, 0.0, 0.0];
        let correction = amplitude * (0.5 - y) / 3.0;
        for k in 0..6 {
            close(point.strain[k], strain[k], 2e-12);
            close(
                point.constitutive_strain[k],
                strain[k] + if k < 3 { correction } else { 0.0 },
                2e-12,
            );
        }
        storage += stress.weight * point.energy_density;
    }
    close(2.0 * storage, work, 2e-12);
    let beyond_domain: Vec<_> = u.iter().map(|value| 100.0 * value).collect();
    assert!(op.bulk_strain_stress(&beyond_domain, points.len(), poll).is_err());

    // Translations are admissible on this observation operator because it has
    // no strong zero clamp that would change the translated nodal field.
    let mut checkpoint = |_| poll();
    let mut control = QuadratureControl3::new(
        QuadratureOptions3 {
            depth: 0,
            ..Default::default()
        },
        &mut checkpoint,
    )
    .unwrap();
    let movable = CutElasticity3::build_with_material_field_and_embedded_dirichlet(
        domain(),
        [1; 3],
        &HalfSpace,
        ElasticMaterialField3::Uniform(&source(0.1)),
        &|_| false,
        &|_, _| true,
        options(),
        Default::default(),
        Default::default(),
        &mut control,
    )
    .unwrap();
    let dyadic: Vec<_> = movable
        .nodes()
        .iter()
        .flat_map(|&[x, y, _]| [x * y / 32.0, 0.0, 0.0])
        .collect();
    let shift = [2.0_f64.powi(40), -2.0_f64.powi(39), 2.0_f64.powi(38)];
    let translated: Vec<_> = dyadic
        .iter()
        .enumerate()
        .map(|(i, u)| u + shift[i % 3])
        .collect();
    // Every local nodal difference is exactly representable before and after
    // this shift, so the stress and both strain-domain queries must agree bitwise.
    assert_eq!(
        movable.bulk_stress(&dyadic, 1000, poll).unwrap(),
        movable.bulk_stress(&translated, 1000, poll).unwrap(),
    );
    movable.admit_material_domain(&translated, poll).unwrap();
    let rigid: Vec<_> = (0..movable.n())
        .map(|i| [2.0_f64.powi(80), -2.0_f64.powi(79), 2.0_f64.powi(78)][i % 3])
        .collect();
    movable.admit_material_domain(&rigid, poll).unwrap();
    for point in movable.bulk_stress(&rigid, 1000, poll).unwrap() {
        assert_eq!(point.stress, [0.0; 6]);
    }
}

#[test]
fn all_strain_and_stress_observations_remove_common_translation() {
    for formulation in [
        ElasticityFormulation3::Displacement,
        ElasticityFormulation3::BBar,
    ] {
        let mut checkpoint = |_| poll();
        let mut control = QuadratureControl3::new(
            QuadratureOptions3 {
                depth: 0,
                ..Default::default()
            },
            &mut checkpoint,
        )
        .unwrap();
        let state = ElasticAcousticState::from_bulk_and_shear(
            1000.0, 12.0, 3.0, ContentHash([74; 32]), 0.1,
        )
        .unwrap();
        let op = CutElasticity3::build_with_material_field_and_embedded_dirichlet(
            domain(),
            [1; 3],
            &HalfSpace,
            ElasticMaterialField3::Uniform(&state),
            &|_| false,
            &|_, _| true,
            ElasticityOptions3 { formulation, ..options() },
            Default::default(),
            Default::default(),
            &mut control,
        )
        .unwrap();
        let u: Vec<_> = op
            .nodes()
            .iter()
            .flat_map(|&[x, y, _]| [x * y / 32.0, 0.0, 0.0])
            .collect();
        let shift = [2.0_f64.powi(40), -2.0_f64.powi(39), 2.0_f64.powi(38)];
        let shifted: Vec<_> = u.iter().enumerate().map(|(i, u)| u + shift[i % 3]).collect();
        // Dyadic nodal differences survive the shift exactly. Removing it
        // before differentiation therefore preserves every response bit.
        assert_eq!(
            op.bulk_strain_stress(&u, 1000, poll).unwrap(),
            op.bulk_strain_stress(&shifted, 1000, poll).unwrap()
        );
        assert_eq!(
            op.bulk_stress(&u, 1000, poll).unwrap(),
            op.bulk_stress(&shifted, 1000, poll).unwrap()
        );
        op.admit_material_domain(&shifted, poll).unwrap();
        let rigid: Vec<_> = (0..op.n())
            .map(|i| [2.0_f64.powi(80), -2.0_f64.powi(79), 2.0_f64.powi(78)][i % 3])
            .collect();
        op.admit_material_domain(&rigid, poll).unwrap();
        for point in op.bulk_strain_stress(&rigid, 1000, poll).unwrap() {
            assert_eq!(point.strain, [0.0; 6]);
            assert_eq!(point.constitutive_strain, [0.0; 6]);
            assert_eq!(point.response.stress, [0.0; 6]);
            assert_eq!(point.energy_density, 0.0);
        }
    }
}

// Cartesian and adaptive APIs must share exactly the same constitutive map;
// the latter additionally exercises hanging-node expansion and its transpose.
macro_rules! check_pullbacks {
    ($operator:expr) => {{
        let mut op = $operator;
        let scales: Vec<_> = (0..op.cells()).map(|i| 0.4 + 0.01 * i as f64).collect();
        op.set_scales(&scales).unwrap();
        let u: Vec<_> = (0..op.n())
            .map(|i| 0.0003 * ((i % 11) as f64 - 5.0))
            .collect();
        let v: Vec<_> = (0..op.n())
            .map(|i| 0.0002 * ((i % 7) as f64 - 3.0))
            .collect();
        let points = op.bulk_stress(&u, 10_000, poll).unwrap();
        let dual: Vec<_> = points
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
            .collect();
        let pb = op.bulk_stress_pullback(&u, &dual, poll).unwrap();
        let reference = op.reference_bulk_stress_pullback(&dual, poll).unwrap();
        close(dot(&reference, &u), observed(&points, &dual, true), 2e-12);
        let h = 1e-4;
        let plus: Vec<_> = u.iter().zip(&v).map(|(u, v)| u + h * v).collect();
        let minus: Vec<_> = u.iter().zip(&v).map(|(u, v)| u - h * v).collect();
        let fd = (observed(&op.bulk_stress(&plus, 10_000, poll).unwrap(), &dual, false)
            - observed(&op.bulk_stress(&minus, 10_000, poll).unwrap(), &dual, false))
            / (2.0 * h);
        close(fd, dot(&pb.displacement, &v), 3e-9);
        let ds: Vec<_> = (0..op.cells())
            .map(|i| 0.1 + 0.03 * (i % 3) as f64)
            .collect();
        let plus: Vec<_> = scales.iter().zip(&ds).map(|(s, d)| s + h * d).collect();
        let minus: Vec<_> = scales.iter().zip(&ds).map(|(s, d)| s - h * d).collect();
        op.set_scales(&plus).unwrap();
        let jp = observed(&op.bulk_stress(&u, 10_000, poll).unwrap(), &dual, false);
        op.set_scales(&minus).unwrap();
        let jm = observed(&op.bulk_stress(&u, 10_000, poll).unwrap(), &dual, false);
        close((jp - jm) / (2.0 * h), dot(&pb.scales, &ds), 3e-9);
        op.set_scales(&scales).unwrap();
        let family = op
            .affine_bulk_stress(&u, &[&v], 2 * points.len(), poll)
            .unwrap();
        assert_eq!(family.baseline, points);
        let direction = op.bulk_stress(&v, 10_000, poll).unwrap();
        for (basis, p) in family.basis_reference_stresses[0].iter().zip(direction) {
            assert_eq!(basis, &p.reference_stress);
        }
        let work = op.scale_bilinear_forms(&v, &u, poll).unwrap();
        let mut ku = vec![0.0; op.n()];
        op.apply(&u, &mut ku);
        for (node, fixed) in op.fixed().iter().enumerate() {
            if *fixed {
                ku[3 * node..3 * node + 3].fill(0.0);
            }
        }
        close(dot(&work, &scales), dot(&v, &ku), 2e-12);
    }};
}

#[test]
fn cartesian_and_hanging_stress_pullbacks_match_directional_work() {
    check_pullbacks!(cartesian(2, 0.1, ElasticityFormulation3::BBar));
    let op = adaptive();
    assert!(op.physical_nodes().len() > op.nodes().len());
    check_pullbacks!(op);
}

#[test]
fn source_admission_checks_actual_and_projected_strain_for_entire_affine_family() {
    let op = cartesian(1, 0.1, ElasticityFormulation3::BBar);
    let standard = cartesian(1, 0.1, ElasticityFormulation3::Displacement);
    // All raw components stay below 0.1 at retained points, but the projection
    // adds over 0.01 to exx on the left side of the cell.
    let projected_exit: Vec<_> = op
        .nodes()
        .iter()
        .flat_map(|&[x, y, _]| [0.09 * x, 0.15 * (x - 0.5) * y, 0.0])
        .collect();
    standard.bulk_stress(&projected_exit, 1000, poll).unwrap();
    assert!(matches!(
        op.bulk_stress(&projected_exit, 1000, poll),
        Err(ElasticityError3::Invalid(
            "projected bulk strain exceeds retained material domain"
        ))
    ));
    // Conversely this field's projected strain fits, while actual exx exceeds
    // the calibrated ceiling near y=1. B-bar must not hide that excursion.
    let raw_exit: Vec<_> = op
        .nodes()
        .iter()
        .flat_map(|&[x, y, _]| [0.13 * x * y, 0.0, 0.0])
        .collect();
    assert!(matches!(
        op.bulk_stress(&raw_exit, 1000, poll),
        Err(ElasticityError3::Invalid(
            "bulk strain exceeds retained material domain"
        ))
    ));
    let zero = vec![0.0; op.n()];
    for exit in [&projected_exit, &raw_exit] {
        assert!(op.admit_material_domain(exit, poll).is_err());
        assert!(op.admit_load_space_material_domain(&[exit], poll).is_err());
        assert!(op.affine_bulk_stress(&zero, &[exit], 2000, poll).is_err());
        let half: Vec<_> = exit.iter().map(|x| 0.5 * x).collect();
        op.admit_affine_material_domain(&half, &[&half], poll)
            .unwrap_err();
        op.admit_material_domain(&half, poll).unwrap();
        op.affine_bulk_stress(&zero, &[&half], 2000, poll).unwrap();
    }
}

macro_rules! check_embedded_extension {
    ($operator:expr, $youngs:expr, $poisson:expr) => {{
        let op = $operator;
        let (youngs, poisson) = ($youngs, $poisson);
        let strain = 0.01;
        let motion = |p: [f64; 3], _: [f64; 3]| {
            [
                0.003 + strain * (p[0] - 0.17),
                -poisson * strain * p[1],
                -poisson * strain * p[2],
            ]
        };
        let traction = |_: [f64; 3], n: [f64; 3]| {
            if n[0] > 0.0 {
                [youngs * strain, 0.0, 0.0]
            } else {
                [0.0; 3]
            }
        };
        let f = op.surface_load(&traction, poll).unwrap().rhs;
        let mut rhs = op.prescribed_displacement_load(&motion, poll).unwrap();
        for (b, f) in rhs.iter_mut().zip(f) {
            *b += f;
        }
        let exact: Vec<_> = op
            .nodes()
            .iter()
            .flat_map(|&p| motion(p, [0.0; 3]))
            .collect();
        let mut applied = vec![0.0; op.n()];
        op.apply(&exact, &mut applied);
        for (a, b) in applied.iter().zip(&rhs) {
            close(*a, *b, 1e-9);
        }
        let solved = op
            .solve_controlled(&rhs, 1e-12, 20_000, 32, |_| poll())
            .unwrap();
        for (actual, expected) in solved.coefficients().iter().zip(&exact) {
            close(*actual, *expected, 2e-7);
        }
        for point in op.bulk_stress(&exact, 100_000, poll).unwrap() {
            close(point.stress[0], youngs * strain, 1e-9);
            for stress in &point.stress[1..] {
                close(*stress, 0.0, 1e-9);
            }
        }
        let reaction = op
            .embedded_reaction(&exact, Some(&motion), &|_, _| [1.0, 0.0, 0.0], poll)
            .unwrap();
        close(reaction.value, -youngs * strain, 1e-9);
    }};
}

#[test]
fn nearly_incompressible_embedded_extension_preserves_affine_motion_and_traction() {
    let (youngs, poisson) = (2.0, 0.4999);
    let material = IsotropicElastic::new(youngs, poisson, 0.1).unwrap();
    let mut checkpoint = |_| poll();
    let mut control = QuadratureControl3::new(
        QuadratureOptions3 {
            depth: 1,
            ..Default::default()
        },
        &mut checkpoint,
    )
    .unwrap();
    let cartesian = CutElasticity3::build_with_embedded_dirichlet(
        domain(),
        [2; 3],
        &Strip,
        &material,
        &|_| false,
        &|_, n| n[0] < 0.0,
        options(),
        Default::default(),
        Default::default(),
        &mut control,
    )
    .unwrap();
    // A weighted reaction distinguishes projected traces from ordinary traces;
    // an affine patch alone would not, because both have the same affine stress.
    let amplitude = 0.02;
    let motion = |p: [f64; 3], _: [f64; 3]| [amplitude * p[0] * p[1], 0.0, 0.0];
    let u: Vec<_> = cartesian
        .nodes()
        .iter()
        .flat_map(|&p| motion(p, [0.0; 3]))
        .collect();
    let reaction = cartesian
        .embedded_reaction(&u, Some(&motion), &|p, _| [p[1], 0.0, 0.0], poll)
        .unwrap();
    let shear = youngs / (2.0 * (1.0 + poisson));
    let bulk = youngs / (3.0 * (1.0 - 2.0 * poisson));
    // Integral y*mean_cell(y) dy = 5/16 on the two equal y intervals.
    close(
        reaction.value,
        -amplitude * (4.0 * shear / 9.0 + 5.0 * bulk / 16.0),
        2e-10,
    );
    check_embedded_extension!(cartesian, youngs, poisson);
    let tree = Octree3::uniform(1, 3, 128).unwrap();
    let tree = tree
        .refined(&[*tree.leaves().iter().next().unwrap()], poll)
        .unwrap();
    let adaptive = AdaptiveElasticity3::build_with_embedded_dirichlet(
        domain(),
        &tree,
        &Strip,
        &material,
        &|_| false,
        &|_, n| n[0] < 0.0,
        options(),
        Default::default(),
        Default::default(),
        &mut control,
    )
    .unwrap();
    assert!(adaptive.physical_nodes().len() > adaptive.nodes().len());
    check_embedded_extension!(adaptive, youngs, poisson);
}

#[test]
fn complete_q1_material_domain_detects_corner_strain_missed_by_bulk_samples() {
    // u_x = a*x*y is represented exactly by Q1. epsilon_xx = a*y reaches
    // 0.101 at y=1, although all interior Gauss points stay below 0.1. The
    // B-bar projection is also inside the ceiling: actual strain must still gate.
    for formulation in [
        ElasticityFormulation3::Displacement,
        ElasticityFormulation3::BBar,
    ] {
        let op = cartesian(1, 0.1, formulation);
        let field = |a: f64| {
            op.nodes()
                .iter()
                .flat_map(|p| [a * p[0] * p[1], 0.0, 0.0])
                .collect::<Vec<_>>()
        };
        let excessive = field(0.101);
        op.bulk_stress(&excessive, 100_000, poll).unwrap();
        assert!(op.admit_material_domain(&excessive, poll).is_err());
        op.admit_material_domain(&field(0.099), poll).unwrap();

        // Each listed endpoint is admissible; the combined coefficient unit
        // ball reaches 0.04 + sqrt(2)*0.045 at the corner and must refuse.
        let offset = field(0.04);
        let direction = field(0.045);
        op.admit_affine_material_domain(&offset, &[&direction], poll)
            .unwrap();
        op.affine_bulk_stress(&offset, &[&direction, &direction], 100_000, poll)
            .unwrap();
        assert!(
            op.admit_affine_material_domain(&offset, &[&direction, &direction], poll)
                .is_err()
        );
        let direction = field(0.075);
        op.admit_load_space_material_domain(&[&direction], poll)
            .unwrap();
        assert!(
            op.admit_load_space_material_domain(&[&direction, &direction], poll)
                .is_err()
        );
    }
}

#[test]
fn adaptive_q1_material_domain_bounds_the_hanging_field_between_bulk_points() {
    let op = adaptive();
    assert!(op.physical_nodes().len() > op.nodes().len());
    let field = |a: f64| {
        op.nodes()
            .iter()
            .flat_map(|p| [a * p[0] * p[1], 0.0, 0.0])
            .collect::<Vec<_>>()
    };
    let excessive = field(0.101);
    op.bulk_stress(&excessive, 100_000, poll).unwrap();
    assert!(op.admit_material_domain(&excessive, poll).is_err());
    op.admit_material_domain(&field(0.099), poll).unwrap();
    let offset = field(0.04);
    let direction = field(0.045);
    op.affine_bulk_stress(&offset, &[&direction, &direction], 300_000, poll)
        .unwrap();
    assert!(
        op.admit_affine_material_domain(&offset, &[&direction, &direction], poll)
            .is_err()
    );
    op.admit_affine_material_domain(&offset, &[&direction], poll)
        .unwrap();
    let direction = field(0.075);
    assert!(
        op.admit_load_space_material_domain(&[&direction, &direction], poll)
            .is_err()
    );
    op.admit_load_space_material_domain(&[&direction], poll)
        .unwrap();
}
