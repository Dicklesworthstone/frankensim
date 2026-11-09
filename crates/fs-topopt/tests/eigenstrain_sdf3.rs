//! Fixed thermal fields through real CutFEM equilibrium and density derivatives.
use std::ops::ControlFlow;

use fs_blake3::ContentHash;
use fs_cutfem::elastic3::adaptive::AdaptiveElasticity3;
use fs_cutfem::elastic3::adaptive::enrichment::precondition::AdaptiveSolveSpace3;
use fs_cutfem::elastic3::eigenstrain::PreparedEigenstrain3;
use fs_cutfem::elastic3::{CutElasticity3, ElasticMaterialField3};
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_material::elastic::acoustic::ElasticAcousticState;
use fs_solver::LinearOp;
use fs_topopt::pipeline::LoadCase;
use fs_topopt::sdf3::eigenstrain::{CutEigenstrainCase3, CutEigenstrainPolicy3};
use fs_topopt::sdf3::{CutDensityStudy3, PhysicalRegion3, Sdf3Elasticity};
use fs_topopt::{EvaluationStop, SimpParams, SolveBudget, SolveControl, SolveProgress};

struct Slab;
impl CutSdf3 for Slab {
    fn value(&self, p: [f64; 3]) -> f64 {
        p[2] - 0.73
    }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        Interval::new(lo[2], hi[2]) - Interval::new(0.73, 0.73)
    }
    fn derivative_enclose(&self, _: [f64; 3], _: [f64; 3], axis: HeightAxis) -> Interval {
        let value = if axis == HeightAxis::Z { 1.0 } else { 0.0 };
        Interval::new(value, value)
    }
}
fn proceed() -> ControlFlow<()> {
    ControlFlow::Continue(())
}
fn policy() -> CutEigenstrainPolicy3 {
    CutEigenstrainPolicy3 {
        maximum_kinematic_strain: 0.1,
        max_points_per_case: 100_000,
    }
}
fn params() -> SimpParams {
    SimpParams {
        penal: 3.0,
        beta: 2.0,
        eta: 0.5,
        e_min: 0.02,
    }
}
fn material(limit: f64) -> ElasticAcousticState {
    ElasticAcousticState::from_isotropic(
        &IsotropicElastic::new(20.0, 0.25, limit).unwrap(),
        3.0,
        ContentHash([4; 32]),
    )
    .unwrap()
}
fn raw(counts: [usize; 3], limit: f64, embedded: bool) -> CutElasticity3 {
    let mut callback = |_| proceed();
    let mut control = QuadratureControl3::new(
        QuadratureOptions3 {
            depth: 1,
            ..Default::default()
        },
        &mut callback,
    )
    .unwrap();
    let domain = HexCell::try_new([0.0; 3], [1.0; 3]).unwrap();
    let material = material(limit);
    if embedded {
        CutElasticity3::build_with_material_field_and_embedded_dirichlet(
            domain,
            counts,
            &Slab,
            ElasticMaterialField3::Uniform(&material),
            &|p| p[0] == 0.0,
            &|_, _| true,
            Default::default(),
            Default::default(),
            Default::default(),
            &mut control,
        )
        .unwrap()
    } else {
        CutElasticity3::build_with_material_field(
            domain,
            counts,
            &Slab,
            ElasticMaterialField3::Uniform(&material),
            &|p| p[0] == 0.0,
            Default::default(),
            &mut control,
        )
        .unwrap()
    }
}
fn adaptive() -> AdaptiveElasticity3 {
    let tree = Octree3::uniform(1, 4, 1000).unwrap();
    let mark = *tree
        .leaves()
        .iter()
        .find(|cell| cell.index() == [0, 0, 1])
        .unwrap();
    let tree = tree.refined(&[mark], proceed).unwrap();
    let mut callback = |_| proceed();
    let mut control = QuadratureControl3::new(
        QuadratureOptions3 {
            depth: 1,
            ..Default::default()
        },
        &mut callback,
    )
    .unwrap();
    AdaptiveElasticity3::build_with_material_field_and_embedded_dirichlet(
        HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(),
        &tree,
        &Slab,
        ElasticMaterialField3::Uniform(&material(0.1)),
        &|p| p[0] == 0.0,
        &|_, _| true,
        Default::default(),
        Default::default(),
        Default::default(),
        &mut control,
    )
    .unwrap()
}
fn free(_: usize, p: [f64; 3]) -> Result<[f64; 6], fs_cutfem::elastic3::ElasticityError3> {
    Ok([
        0.002 * (1.0 + p[0]),
        -0.0005 * p[2],
        0.0003,
        0.0002 * p[1],
        0.0,
        0.0,
    ])
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * a.abs().max(b.abs()).max(1e-12),
        "{a:e} != {b:e}"
    );
}

fn check_chain<O: Sdf3Elasticity>(
    operator: O,
    first: PreparedEigenstrain3,
    second: PreparedEigenstrain3,
) {
    let body = operator
        .prepare_material_body_load(&|p| [0.0001 * p[1], -0.0002, 0.0], 100_000, proceed)
        .unwrap();
    let force = vec![0.0; operator.n()];
    let cases = [
        CutEigenstrainCase3 {
            force: &force,
            weight: 0.3,
            eigenstrain: &first,
        },
        CutEigenstrainCase3 {
            force: &force,
            weight: 0.7,
            eigenstrain: &second,
        },
    ];
    let mut study = CutDensityStudy3::new(operator, 0.12, params())
        .with_material_body_load(body)
        .unwrap();
    let mut regions = vec![PhysicalRegion3::Design; study.cells()];
    regions[0] = PhysicalRegion3::Solid;
    study = study.with_physical_regions(regions).unwrap();
    let rho: Vec<_> = (0..study.cells())
        .map(|i| 0.45 + 0.01 * (i % 9) as f64)
        .collect();
    let mut callback = |_| proceed();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let base = study
        .evaluate_eigenstrain(&rho, &cases, policy(), &mut control)
        .unwrap();
    assert!(base.objective.compliance > 0.0);
    close(
        base.objective.compliance,
        0.3 * base.objective.case_compliances[0] + 0.7 * base.objective.case_compliances[1],
        1e-14,
    );
    assert_ne!(
        base.objective.displacements[0],
        base.objective.displacements[1]
    );
    for direction in [
        vec![1.0; rho.len()],
        (0..rho.len()).map(|i| (i % 3) as f64 - 1.0).collect(),
    ] {
        let h = 1e-4;
        let plus: Vec<_> = rho.iter().zip(&direction).map(|(r, d)| r + h * d).collect();
        let minus: Vec<_> = rho.iter().zip(&direction).map(|(r, d)| r - h * d).collect();
        let a = study
            .evaluate_eigenstrain(&plus, &cases, policy(), &mut control)
            .unwrap();
        let b = study
            .evaluate_eigenstrain(&minus, &cases, policy(), &mut control)
            .unwrap();
        close(
            (a.objective.compliance - b.objective.compliance) / (2.0 * h),
            dot(&base.objective.gradient, &direction),
            4e-4,
        );
        close(
            (a.volume_fraction - b.volume_fraction) / (2.0 * h),
            dot(&base.volume_gradient, &direction),
            5e-5,
        );
    }
}

#[test]
fn independent_thermal_cases_and_reference_density_load_have_the_complete_raw_gradient() {
    let op = raw([2, 1, 2], 0.1, true);
    let a = op.prepare_eigenstrain(&free, 100_000, proceed).unwrap();
    let b = op
        .prepare_eigenstrain(
            &|cell, p| free(cell, p).map(|e| e.map(|v| -0.7 * v)),
            100_000,
            proceed,
        )
        .unwrap();
    check_chain(op, a, b);
}

#[test]
fn adaptive_and_preconditioned_thermal_cases_use_the_same_master_space_chain() {
    let op = adaptive();
    let a = op.prepare_eigenstrain(&free, 100_000, proceed).unwrap();
    let b = op
        .prepare_eigenstrain(
            &|cell, p| free(cell, p).map(|e| e.map(|v| -0.7 * v)),
            100_000,
            proceed,
        )
        .unwrap();
    check_chain(op, a, b);
    let op = adaptive();
    let a = op.prepare_eigenstrain(&free, 100_000, proceed).unwrap();
    let b = op
        .prepare_eigenstrain(
            &|cell, p| free(cell, p).map(|e| e.map(|v| -0.7 * v)),
            100_000,
            proceed,
        )
        .unwrap();
    check_chain(AdaptiveSolveSpace3::jacobi(op, 100_000_000), a, b);
}

#[test]
fn free_expansion_has_positive_work_gradient_and_retains_its_tight_elastic_domain() {
    let op = raw([1; 3], 1e-6, false);
    let free = [0.003, 0.0, 0.0, 0.0, 0.0, 0.0];
    let field = op
        .prepare_eigenstrain(&|_, _| Ok(free), 100_000, proceed)
        .unwrap();
    let force = vec![0.0; op.n()];
    let cases = [CutEigenstrainCase3 {
        force: &force,
        weight: 1.0,
        eigenstrain: &field,
    }];
    let mut study = CutDensityStudy3::new(op, 0.0, params());
    let mut callback = |_| proceed();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let base = study
        .evaluate_eigenstrain(&[0.6], &cases, policy(), &mut control)
        .unwrap();
    assert!(
        base.objective.gradient[0] > 0.0,
        "the thermal load derivative must reverse the fixed-load-only sign"
    );
    for (node, u) in study
        .operator()
        .nodes()
        .iter()
        .zip(base.objective.displacements[0].chunks_exact(3))
    {
        assert!((u[0] - free[0] * node[0]).abs() < 1e-11);
        assert!(u[1].abs() < 1e-11 && u[2].abs() < 1e-11);
    }
    let h = 1e-4;
    let a = study
        .evaluate_eigenstrain(&[0.6 + h], &cases, policy(), &mut control)
        .unwrap();
    let b = study
        .evaluate_eigenstrain(&[0.6 - h], &cases, policy(), &mut control)
        .unwrap();
    close(
        (a.objective.compliance - b.objective.compliance) / (2.0 * h),
        base.objective.gradient[0],
        1e-6,
    );
    let accepted = study.operator().scales().to_vec();
    let tight = CutEigenstrainPolicy3 {
        maximum_kinematic_strain: 0.001,
        ..policy()
    };
    assert!(
        study
            .evaluate_eigenstrain(&[0.5], &cases, tight, &mut control)
            .is_err()
    );
    assert_eq!(study.operator().scales(), accepted);
}

#[test]
fn zero_free_strain_preserves_the_original_cold_functional_and_gradient() {
    let op = raw([2, 1, 2], 0.1, false);
    let field = op
        .prepare_eigenstrain(&|_, _| Ok([0.0; 6]), 100_000, proceed)
        .unwrap();
    let force = op.body_load(&|_| [0.0, -0.0001, 0.0], proceed).unwrap();
    let cases = [CutEigenstrainCase3 {
        force: &force,
        weight: 1.0,
        eigenstrain: &field,
    }];
    let loads = [LoadCase {
        force: &force,
        weight: 1.0,
    }];
    let mut study = CutDensityStudy3::new(op, 0.12, params());
    let rho = vec![0.55; study.cells()];
    let mut callback = |_| proceed();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let cold = study.evaluate(&rho, &loads, &mut control).unwrap();
    let thermal = study
        .evaluate_eigenstrain(&rho, &cases, policy(), &mut control)
        .unwrap();
    assert_eq!(cold.objective.compliance, thermal.objective.compliance);
    assert_eq!(
        cold.objective.displacements,
        thermal.objective.displacements
    );
    assert_eq!(cold.objective.gradient, thermal.objective.gradient);
}

#[test]
fn foreign_preparation_point_allowance_and_cancelled_trials_preserve_the_operator() {
    let op = raw([2, 1, 2], 0.1, true);
    let field = op.prepare_eigenstrain(&free, 100_000, proceed).unwrap();
    let foreign = raw([2, 1, 2], 0.1, true)
        .prepare_eigenstrain(&free, 100_000, proceed)
        .unwrap();
    let force = vec![0.0; op.n()];
    let cases = [CutEigenstrainCase3 {
        force: &force,
        weight: 1.0,
        eigenstrain: &field,
    }];
    let mut study = CutDensityStudy3::new(op, 0.1, params());
    let rho = vec![0.6; study.cells()];
    let mut callback = |_| proceed();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    study
        .evaluate_eigenstrain(&rho, &cases, policy(), &mut control)
        .unwrap();
    let accepted = study.operator().scales().to_vec();
    let wrong = [CutEigenstrainCase3 {
        eigenstrain: &foreign,
        ..cases[0]
    }];
    assert!(
        study
            .evaluate_eigenstrain(&rho, &wrong, policy(), &mut control)
            .is_err()
    );
    let small = CutEigenstrainPolicy3 {
        max_points_per_case: 1,
        ..policy()
    };
    assert!(
        study
            .evaluate_eigenstrain(&rho, &cases, small, &mut control)
            .is_err()
    );
    assert_eq!(study.operator().scales(), accepted);
    for stage in [
        "sdf3-eigenstrain-load",
        "sdf3-eigenstrain-domain",
        "sdf3-eigenstrain-pullback",
        "sdf3-evaluation-publish",
    ] {
        let mut callback = |progress: SolveProgress| {
            if progress.stage == stage {
                ControlFlow::Break(())
            } else {
                proceed()
            }
        };
        let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
        assert!(matches!(
            study.evaluate_eigenstrain(&vec![0.4; study.cells()], &cases, policy(), &mut control),
            Err(EvaluationStop::Cancelled)
        ));
        assert_eq!(study.operator().scales(), accepted);
    }
}
