//! G0/G3/G4: material-aware 3-D modal fields and sensitivities.
//! Small complete search spaces make the numerical comparisons independent of
//! the rate of convergence of block inverse iteration.
use std::ops::ControlFlow;

use fs_blake3::ContentHash;
use fs_cutfem::elastic3::adaptive::AdaptiveElasticity3;
use fs_cutfem::elastic3::{CutElasticity3, ElasticMaterialField3, MaterialMass3};
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_material::elastic::acoustic::ElasticAcousticState;
use fs_solver::LinearOp;
use fs_topopt::modal::{MatrixFreeEigenError, MatrixFreeEigenOptions};
use fs_topopt::sdf3::modal::{
    Sdf3EigenfrequencyEvaluation3, Sdf3EigenfrequencyOptions3, controlled_sdf3_eigenpairs,
};
use fs_topopt::sdf3::{CutDensityStudy3, PhysicalRegion3, Sdf3Elasticity};
use fs_topopt::{EvaluationStop, SimpParams, SolveBudget, SolveControl, SolveProgress};

struct Top;
impl CutSdf3 for Top {
    fn value(&self, p: [f64; 3]) -> f64 {
        p[2] - 0.53
    }

    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        Interval::new(lo[2], hi[2]) - Interval::new(0.53, 0.53)
    }

    fn derivative_enclose(&self, _: [f64; 3], _: [f64; 3], axis: HeightAxis) -> Interval {
        let d = if axis == HeightAxis::Z { 1.0 } else { 0.0 };
        Interval::new(d, d)
    }
}

fn poll() -> ControlFlow<()> {
    ControlFlow::Continue(())
}

fn domain() -> HexCell {
    HexCell::try_new([0.0; 3], [1.6, 0.9, 0.7]).unwrap()
}

fn params() -> SimpParams {
    SimpParams {
        e_min: 0.02,
        penal: 3.0,
        beta: 2.0,
        eta: 0.5,
    }
}

fn material(density: f64, stiffness: f64) -> ElasticAcousticState {
    let diagonal = [20.0, 17.0, 29.0, 8.0, 13.0, 11.0];
    let tensor = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            if i == j {
                stiffness * diagonal[i]
            } else {
                0.0
            }
        })
    });
    // A normalized eigenvector is an arbitrary algebraic amplitude, not a
    // physical displacement experiment. This tiny source strain ceiling must
    // therefore not reject an otherwise admissible vibration calculation.
    ElasticAcousticState::try_new_mandel(
        tensor,
        density,
        ContentHash([(density as u8) * 3 + stiffness as u8; 32]),
        1e-12,
    )
    .unwrap()
}

fn cartesian(densities: [f64; 2], stiffness: f64) -> CutElasticity3 {
    let mut callback = |_| poll();
    let mut quadrature = QuadratureControl3::new(
        QuadratureOptions3 {
            depth: 1,
            ..Default::default()
        },
        &mut callback,
    )
    .unwrap();
    let field = |cell: HexCell| {
        Ok(material(
            densities[usize::from(cell.lo()[0] >= 0.8)],
            stiffness,
        ))
    };
    CutElasticity3::build_with_material_field(
        domain(),
        [2, 2, 1],
        &Top,
        ElasticMaterialField3::ByCell(&field),
        &|p| p[0] == 0.0,
        Default::default(),
        &mut quadrature,
    )
    .unwrap()
}

fn options() -> Sdf3EigenfrequencyOptions3 {
    Sdf3EigenfrequencyOptions3 {
        solver: MatrixFreeEigenOptions {
            count: 2,
            oversampling: 60,
            max_iterations: 8,
            relative_tolerance: 1e-9,
            linear_tolerance: 1e-12,
            seed: 17,
        },
        aggregation_beta: 30.0,
    }
}

fn evaluate<O: Sdf3Elasticity>(
    study: &mut CutDensityStudy3<O>,
    mass: &MaterialMass3,
    rho: &[f64],
) -> Sdf3EigenfrequencyEvaluation3 {
    let mut callback = |_| poll();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    study
        .evaluate_eigenfrequency(mass, rho, options(), &mut control)
        .unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    let scale = actual.abs().max(expected.abs()).max(1e-10);
    assert!(
        actual.is_finite()
            && expected.is_finite()
            && (actual - expected).abs() <= tolerance * scale,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}",
    );
}

fn dot(left: &[f64], right: &[f64]) -> f64 {
    left.iter().zip(right).map(|(a, b)| a * b).sum()
}

fn norm(values: &[f64]) -> f64 {
    values.iter().fold(0.0_f64, |n, v| n.hypot(*v))
}

#[test]
fn source_density_and_tensor_scale_the_same_native_vibration_pencil() {
    let mut light = CutDensityStudy3::new(cartesian([2.0; 2], 1.0), 0.12, params());
    let mut heavy = CutDensityStudy3::new(cartesian([4.0; 2], 1.0), 0.12, params());
    let mut stiff = CutDensityStudy3::new(cartesian([2.0; 2], 2.0), 0.12, params());
    let light_mass = light
        .operator()
        .prepare_material_mass(100_000, poll)
        .unwrap();
    let heavy_mass = heavy
        .operator()
        .prepare_material_mass(100_000, poll)
        .unwrap();
    let stiff_mass = stiff
        .operator()
        .prepare_material_mass(100_000, poll)
        .unwrap();
    let rho = [0.47, 0.56, 0.63, 0.71];
    let first = evaluate(&mut light, &light_mass, &rho);
    let second = evaluate(&mut heavy, &heavy_mass, &rho);
    let third = evaluate(&mut stiff, &stiff_mass, &rho);
    for i in 0..options().solver.count {
        close(
            second.objective.modal.values[i],
            0.5 * first.objective.modal.values[i],
            2e-8,
        );
        close(
            third.objective.modal.values[i],
            2.0 * first.objective.modal.values[i],
            2e-8,
        );
    }
    close(second.physical_mass_kg, 2.0 * first.physical_mass_kg, 1e-12);
    close(third.physical_mass_kg, first.physical_mass_kg, 1e-12);
    close(second.volume_fraction, first.volume_fraction, 1e-12);
    for (heavy, light) in second
        .physical_mass_gradient
        .iter()
        .zip(&first.physical_mass_gradient)
    {
        close(*heavy, 2.0 * light, 1e-11);
    }
    assert!(
        light
            .operator()
            .admit_material_domain(&first.objective.modal.modes[0], poll)
            .is_err(),
        "the same arbitrary normalized amplitude exceeds the physical source domain",
    );
    assert!(first.next_eigenvalue.unwrap() > first.objective.modal.values[1]);
    assert!(
        first.objective.modal.relative_residuals.iter().all(|r| *r <= 1e-9)
    );
}

#[test]
fn heterogeneous_frequency_derivatives_include_mass_filter_and_protected_regions() {
    // Exercise both mass-interpolation branches, away from their nondifferentiable
    // switch. Protect one cell only in the regular-density experiment so the
    // low-density case still probes the full rho^6 inertia derivative.
    for (rho, protected) in [
        ([0.43, 0.58, 0.66, 0.51], true),
        ([0.07, 0.09, 0.11, 0.10], false),
    ] {
        let operator = cartesian([2.0, 5.0], 1.0);
        let mass = operator.prepare_material_mass(100_000, poll).unwrap();
        let mut study = CutDensityStudy3::new(operator, 0.12, params());
        if protected {
            let mut regions = vec![PhysicalRegion3::Design; study.cells()];
            regions[0] = PhysicalRegion3::Solid;
            study = study.with_physical_regions(regions).unwrap();
        }
        let original = study.operator().scales().to_vec();
        let result = evaluate(&mut study, &mass, &rho);
        if protected {
            assert_eq!(result.projected_rho[0], 1.0);
        } else {
            assert!(result.projected_rho.iter().all(|r| *r > 0.0 && *r < 0.1));
        }
        for i in 0..rho.len() {
            let h = 1e-5;
            let mut plus = rho;
            let mut minus = rho;
            plus[i] += h;
            minus[i] -= h;
            let plus = evaluate(&mut study, &mass, &plus);
            let minus = evaluate(&mut study, &mass, &minus);
            close(
                (plus.objective.aggregate - minus.objective.aggregate) / (2.0 * h),
                result.objective.gradient[i],
                8e-4,
            );
            close(
                (plus.physical_mass_kg - minus.physical_mass_kg) / (2.0 * h),
                result.physical_mass_gradient[i],
                2e-6,
            );
            close(
                (plus.volume_fraction - minus.volume_fraction) / (2.0 * h),
                result.volume_gradient[i],
                2e-6,
            );
        }
        assert_eq!(study.operator().scales(), original);
    }
}

#[test]
fn adaptive_modes_use_the_reduced_hanging_space_and_consistent_mass_metric() {
    let tree = Octree3::uniform(1, 3, 128).unwrap();
    let mark = *tree
        .leaves()
        .iter()
        .find(|leaf| leaf.index() == [1, 0, 0])
        .unwrap();
    let tree = tree.refined(&[mark], poll).unwrap();
    let material = material(3.0, 1.0);
    let mut callback = |_| poll();
    let mut quadrature = QuadratureControl3::new(
        QuadratureOptions3 { depth: 1, ..Default::default() },
        &mut callback,
    )
    .unwrap();
    // Coordinate-plane clamps have exactly clamped parent traces, including
    // hanging nodes. They keep this complete spectral fixture small without
    // clamping an interior hanging node against an unconstrained parent.
    let operator = AdaptiveElasticity3::build_with_material_field(
        domain(),
        &tree,
        &Top,
        ElasticMaterialField3::Uniform(&material),
        &|p| p[0] == 0.0 || p[1] == 0.0 || p[2] == 0.0,
        Default::default(),
        &mut quadrature,
    )
    .unwrap();
    assert!(operator.physical_nodes().len() > operator.nodes().len());
    assert!(3 * operator.fixed().iter().filter(|fixed| !**fixed).count() <= 64);
    let mass = operator.prepare_material_mass(100_000, poll).unwrap();
    let occupancy = vec![0.8; operator.cells()];
    let mut callback = |_| poll();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let report = controlled_sdf3_eigenpairs(
        &operator, &mass, &occupancy, options().solver, &mut control,
    )
    .unwrap();
    for (i, mode) in report.modes.iter().enumerate() {
        assert_eq!(mode.len(), operator.n());
        for (node, &fixed) in operator.fixed().iter().enumerate() {
            if fixed {
                assert_eq!(&mode[3 * node..3 * node + 3], &[0.0; 3]);
            }
        }
        let physical = operator.physical_displacements(mode).unwrap();
        assert!(physical.len() > mode.len());
        let inertia = mass.apply(&occupancy, mode, poll).unwrap();
        close(dot(mode, &inertia), 1.0, 2e-9);
        for earlier in &report.modes[..i] {
            assert!(dot(earlier, &inertia).abs() < 2e-8);
        }
        let mut stiffness = vec![0.0; operator.n()];
        operator.apply(mode, &mut stiffness);
        let lambda = report.values[i];
        let residual: Vec<_> = stiffness
            .iter()
            .zip(&inertia)
            .map(|(k, m)| k - lambda * m)
            .collect();
        assert!(
            norm(&residual) / (norm(&stiffness) + lambda.abs() * norm(&inertia)) <= 2e-9
        );
    }
}

#[test]
fn absent_or_foreign_inertia_and_cancelled_modal_trials_publish_no_result() {
    let mut callback = |_| poll();
    let mut quadrature =
        QuadratureControl3::new(QuadratureOptions3::default(), &mut callback).unwrap();
    let legacy = CutElasticity3::build(
        domain(),
        [2, 2, 1],
        &Top,
        &IsotropicElastic::new(20.0, 0.2, 1.0).unwrap(),
        &|p| p[0] == 0.0,
        Default::default(),
        &mut quadrature,
    )
    .unwrap();
    assert!(legacy.prepare_material_mass(100_000, poll).is_err());
    let mut operator = cartesian([2.0, 5.0], 1.0);
    let mass = operator.prepare_material_mass(100_000, poll).unwrap();
    let incoming = vec![0.37; operator.cells()];
    operator.set_scales(&incoming).unwrap();
    let mut study = CutDensityStudy3::new(operator, 0.12, params());
    let rho = [0.43, 0.58, 0.66, 0.51];
    let mut stop = |progress: SolveProgress| {
        if progress.stage == "sdf3-eigen-setup" { ControlFlow::Break(()) } else { poll() }
    };
    let mut control = SolveControl::new(SolveBudget::default(), &mut stop);
    assert!(matches!(
        study.evaluate_eigenfrequency(&mass, &rho, options(), &mut control),
        Err(MatrixFreeEigenError::Operator(EvaluationStop::Cancelled)),
    ));
    assert!(control.work().linear_solves > 0);
    assert_eq!(study.operator().scales(), incoming);
    let accepted = evaluate(&mut study, &mass, &rho);
    assert_eq!(study.operator().scales(), incoming);
    assert!(accepted.objective.aggregate.is_finite());

    let mut other = CutDensityStudy3::new(cartesian([2.0, 5.0], 1.0), 0.12, params());
    let other_scales = other.operator().scales().to_vec();
    let mut callback = |_| poll();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    assert!(other.evaluate_eigenfrequency(&mass, &rho, options(), &mut control).is_err());
    assert_eq!(other.operator().scales(), other_scales);
}

#[test]
fn both_design_goals_use_real_modes_and_recheck_the_original_physical_limit() {
    use fs_ascent::projected_al::ProjectedAlOptions;
    use fs_topopt::sdf3::modal::design::{
        Sdf3EigenfrequencyDesign3, Sdf3EigenfrequencyDesignOptions3, Sdf3EigenfrequencyGoal3,
    };
    use fs_topopt::sdf3::response::limits::design::MinimumMaterialObjective3;

    for minimum_mass in [false, true] {
        let operator = cartesian([2.0, 5.0], 1.0);
        let mass = operator.prepare_material_mass(100_000, poll).unwrap();
        let mut study = CutDensityStudy3::new(operator, 0.12, params());
        let incoming = study.operator().scales().to_vec();
        let rho = [0.64, 0.71, 0.68, 0.74];
        let baseline = evaluate(&mut study, &mass, &rho);
        assert!(baseline.objective.aggregate > 0.0);
        let reference = 0.45 * baseline.objective.aggregate;
        let goal = if minimum_mass {
            Sdf3EigenfrequencyGoal3::MinimumMaterial {
                objective: MinimumMaterialObjective3::MassKg {
                    scale_kg: baseline.physical_mass_kg,
                },
                volume_scale: 1.0,
            }
        } else {
            Sdf3EigenfrequencyGoal3::MaximumFrequency { volume_limit: 0.97 }
        };
        let settings = Sdf3EigenfrequencyDesignOptions3 {
            frequency: options(),
            eigenvalue_reference: reference,
            density_floor: 0.3,
            optimizer: ProjectedAlOptions {
                max_spectral_step: 0.01,
                ..Default::default()
            },
            ..Default::default()
        };
        let mut callback = |_| poll();
        let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
        let best = {
            let mut design = Sdf3EigenfrequencyDesign3::new(
                &mut study, &mass, &rho, goal, settings, &mut control,
            )
            .unwrap();
            assert!(design.gradient_check().passed());
            design.run(3).unwrap();
            assert!(design.history().len() > 1);
            design.recheck_best().unwrap().unwrap()
        };
        if minimum_mass {
            assert!(
                best.evaluation.physical_mass_kg
                    < baseline.physical_mass_kg * (1.0 - 1e-7)
            );
            assert!(best.evaluation.objective.aggregate >= reference);
            assert!(best.evaluation.objective.modal.values[0] >= reference);
        } else {
            assert!(
                best.evaluation.objective.aggregate
                    > baseline.objective.aggregate * (1.0 + 1e-7)
            );
            assert!(best.evaluation.volume_fraction <= 0.97);
        }
        assert_eq!(study.operator().scales(), incoming);
    }
}

#[path = "material_modal_sdf3/adaptive.rs"]
mod adaptive;
