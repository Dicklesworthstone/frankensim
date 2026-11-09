//! G0/G3: thermal design gradients retain load dependence and full tensor stress.

use std::ops::ControlFlow;

use fs_blake3::ContentHash;
use fs_cutfem::elastic3::eigenstrain::EigenstrainResponse3;
use fs_cutfem::elastic3::{
    CutElasticity3, ElasticMaterialField3, ElasticityError3, ElasticityFormulation3,
    ElasticityOptions3,
};
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_material::elastic::acoustic::ElasticAcousticState;
use fs_solver::op::LinearOp;

struct Slab;
impl CutSdf3 for Slab {
    fn value(&self, p: [f64; 3]) -> f64 {
        (p[0] - 0.17) * (p[0] - 0.83)
    }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        let x = Interval::new(lo[0], hi[0]);
        (x - Interval::new(0.17, 0.17)) * (x - Interval::new(0.83, 0.83))
    }
    fn derivative_enclose(&self, lo: [f64; 3], hi: [f64; 3], a: HeightAxis) -> Interval {
        if a == HeightAxis::X {
            Interval::new(2.0 * lo[0] - 1.0, 2.0 * hi[0] - 1.0)
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
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn close(a: f64, b: f64, relative: f64, absolute: f64) {
    assert!(
        (a - b).abs() <= absolute + relative * a.abs().max(b.abs()),
        "{a} != {b}"
    );
}
fn thermal(_: usize, p: [f64; 3]) -> Result<[f64; 6], ElasticityError3> {
    Ok([
        0.001 * (1.0 + p[0]),
        -0.0002 * p[2],
        0.0003 * p[1],
        0.0001 * p[1],
        0.0002 * p[0],
        -0.0001 * p[2],
    ])
}

#[test]
fn complete_thermal_compliance_gradient_requires_both_load_work_terms() {
    for formulation in [
        ElasticityFormulation3::Displacement,
        ElasticityFormulation3::BBar,
    ] {
        let mut callback = |_| poll();
        let mut control =
            QuadratureControl3::new(QuadratureOptions3::default(), &mut callback).unwrap();
        let poisson = if formulation == ElasticityFormulation3::BBar {
            0.49
        } else {
            0.2
        };
        let mut op = CutElasticity3::build_with_embedded_dirichlet(
            domain(),
            [2; 3],
            &Slab,
            &IsotropicElastic::new(20.0, poisson, 0.1).unwrap(),
            &|_| false,
            &|_, normal| normal[0] < 0.0,
            ElasticityOptions3 {
                formulation,
                ..Default::default()
            },
            Default::default(),
            Default::default(),
            &mut control,
        )
        .unwrap();
        let prepared = op.prepare_eigenstrain(&thermal, 100_000, poll).unwrap();
        let base: Vec<_> = (0..op.cells()).map(|i| 0.7 + 0.02 * i as f64).collect();
        op.set_scales(&base).unwrap();
        let external = op.body_load(&|_| [0.0, 0.0001, -0.0002], poll).unwrap();
        let solve = |op: &CutElasticity3| {
            let free = op.eigenstrain_load(&prepared, poll).unwrap();
            let rhs: Vec<_> = external.iter().zip(free).map(|(f, b)| f + b).collect();
            op.solve_controlled(&rhs, 1e-11, 30_000, 32, |_| poll())
                .unwrap()
        };
        let solution = solve(&op);
        op.eigenstrain_response(&prepared, solution.coefficients(), 0.05, 100_000, poll)
            .unwrap();
        let load = op
            .eigenstrain_load_scale_pullback(&prepared, solution.coefficients(), poll)
            .unwrap();
        let stiffness = op.scale_quadratic_forms(solution.coefficients()).unwrap();
        let gradient: Vec<_> = load
            .iter()
            .zip(&stiffness)
            .map(|(b, k)| 2.0 * b - k)
            .collect();
        let direction: Vec<_> = (0..op.cells())
            .map(|i| 0.3 + 0.1 * (i % 3) as f64)
            .collect();
        let h = 1e-4;
        op.set_scales(
            &base
                .iter()
                .zip(&direction)
                .map(|(s, d)| s + h * d)
                .collect::<Vec<_>>(),
        )
        .unwrap();
        let plus = solve(&op).compliance();
        op.set_scales(
            &base
                .iter()
                .zip(&direction)
                .map(|(s, d)| s - h * d)
                .collect::<Vec<_>>(),
        )
        .unwrap();
        let minus = solve(&op).compliance();
        let finite_difference = (plus - minus) / (2.0 * h);
        close(finite_difference, dot(&gradient, &direction), 5e-5, 1e-11);
        assert!(
            (finite_difference + dot(&stiffness, &direction)).abs()
                > 0.01 * finite_difference.abs().max(1e-12)
        );
        op.set_scales(&base).unwrap();
        assert_eq!(solution.coefficients(), solve(&op).coefficients());
    }
}

fn tensor_operator() -> CutElasticity3 {
    let mut c = [[0.0; 6]; 6];
    for (i, row) in c.iter_mut().enumerate() {
        row[i] = 12.0 + i as f64;
    }
    for (i, j, value) in [(0, 3, 1.2), (1, 4, -0.7), (2, 5, 0.8), (3, 5, 0.4)] {
        c[i][j] = value;
        c[j][i] = value;
    }
    // The total free expansion is much larger than the admissible elastic strain.
    let material =
        ElasticAcousticState::try_new_mandel(c, 2.0, ContentHash([83; 32]), 0.002).unwrap();
    let mut callback = |_| poll();
    let mut control =
        QuadratureControl3::new(QuadratureOptions3::default(), &mut callback).unwrap();
    CutElasticity3::build_with_material_field(
        domain(),
        [2; 3],
        &Slab,
        ElasticMaterialField3::Uniform(&material),
        &|p| p == [0.0; 3],
        Default::default(),
        &mut control,
    )
    .unwrap()
}

fn observed(points: &[EigenstrainResponse3], derivatives: &[[f64; 6]]) -> f64 {
    points
        .iter()
        .zip(derivatives)
        .map(|(p, d)| dot(&p.stress, d))
        .sum()
}

#[test]
fn full_tensor_thermal_stress_transposes_admit_elastic_strain_and_match_finite_differences() {
    let mut op = tensor_operator();
    let prepared = op
        .prepare_eigenstrain(&|_, _| Ok([0.01, 0.01, 0.01, 0.0, 0.0, 0.0]), 100_000, poll)
        .unwrap();
    let base: Vec<_> = (0..op.cells()).map(|i| 0.7 + 0.02 * i as f64).collect();
    op.set_scales(&base).unwrap();
    let u: Vec<_> = op
        .nodes()
        .iter()
        .flat_map(|&[x, y, z]| {
            [
                0.0104 * x + 0.0002 * y,
                0.0097 * y + 0.0001 * z,
                0.0102 * z - 0.0001 * x,
            ]
        })
        .collect();
    let points = op
        .eigenstrain_response(&prepared, &u, 0.05, 100_000, poll)
        .unwrap();
    assert!(op.bulk_stress(&u, 100_000, poll).is_err());
    let derivatives: Vec<_> = points
        .iter()
        .map(|p| [1.0 + p.position[0], -0.2, 0.4, 0.3, -0.7, 0.9].map(|v| p.weight * v))
        .collect();
    let pullback = op
        .eigenstrain_stress_pullback(&prepared, &u, 0.05, &derivatives, poll)
        .unwrap();
    let direction: Vec<_> = (0..op.n())
        .map(|i| 0.0001 * ((i % 7) as f64 - 3.0))
        .collect();
    let h = 1e-4;
    let plus: Vec<_> = u.iter().zip(&direction).map(|(u, d)| u + h * d).collect();
    let minus: Vec<_> = u.iter().zip(&direction).map(|(u, d)| u - h * d).collect();
    let jp = observed(
        &op.eigenstrain_response(&prepared, &plus, 0.05, points.len(), poll)
            .unwrap(),
        &derivatives,
    );
    let jm = observed(
        &op.eigenstrain_response(&prepared, &minus, 0.05, points.len(), poll)
            .unwrap(),
        &derivatives,
    );
    close(
        (jp - jm) / (2.0 * h),
        dot(&pullback.displacement, &direction),
        2e-8,
        1e-12,
    );
    let ds: Vec<_> = (0..op.cells()).map(|i| 0.1 * (1 + i % 3) as f64).collect();
    op.set_scales(
        &base
            .iter()
            .zip(&ds)
            .map(|(s, d)| s + h * d)
            .collect::<Vec<_>>(),
    )
    .unwrap();
    let jp = observed(
        &op.eigenstrain_response(&prepared, &u, 0.05, points.len(), poll)
            .unwrap(),
        &derivatives,
    );
    op.set_scales(
        &base
            .iter()
            .zip(&ds)
            .map(|(s, d)| s - h * d)
            .collect::<Vec<_>>(),
    )
    .unwrap();
    let jm = observed(
        &op.eigenstrain_response(&prepared, &u, 0.05, points.len(), poll)
            .unwrap(),
        &derivatives,
    );
    close(
        (jp - jm) / (2.0 * h),
        dot(&pullback.scales, &ds),
        2e-8,
        1e-12,
    );
    op.set_scales(&base).unwrap();

    let dual: Vec<_> = (0..op.n())
        .map(|i| 1000.0 * ((i % 7) as f64 - 3.0))
        .collect();
    let load_pullback = op
        .eigenstrain_load_scale_pullback(&prepared, &dual, poll)
        .unwrap();
    close(
        dot(&dual, &op.eigenstrain_load(&prepared, poll).unwrap()),
        dot(&base, &load_pullback),
        2e-13,
        1e-10,
    );
    let other = tensor_operator();
    assert!(
        other
            .eigenstrain_load_scale_pullback(&prepared, &dual, poll)
            .is_err()
    );
    assert!(
        op.eigenstrain_stress_pullback(
            &prepared,
            &u,
            0.05,
            &derivatives[..derivatives.len() - 1],
            poll
        )
        .is_err()
    );
    assert!(matches!(
        op.eigenstrain_load_scale_pullback(&prepared, &dual, || ControlFlow::Break(())),
        Err(ElasticityError3::Cancelled)
    ));
    assert!(matches!(
        op.eigenstrain_stress_pullback(
            &prepared,
            &u,
            0.05,
            &derivatives,
            || ControlFlow::Break(())
        ),
        Err(ElasticityError3::Cancelled)
    ));
    assert_eq!(op.scales(), base);
    assert_eq!(
        points,
        op.eigenstrain_response(&prepared, &u, 0.05, points.len(), poll)
            .unwrap()
    );
}
