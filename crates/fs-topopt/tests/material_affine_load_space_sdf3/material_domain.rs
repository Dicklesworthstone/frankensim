//! G1/G3/G4: material-domain admission is independent of requested objectives.
use super::*;
use fs_topopt::sdf3::Sdf3Elasticity;

fn material() -> ElasticAcousticState {
    let tensor = std::array::from_fn(|i| std::array::from_fn(|j| if i == j { 10.0 } else { 0.0 }));
    ElasticAcousticState::try_new_mandel(tensor, 3.0, ContentHash([93; 32]), 0.1).unwrap()
}

fn adaptive(tree: &Octree3) -> AdaptiveSolveSpace3 {
    let mut callback = |_| poll();
    let mut control = QuadratureControl3::new(
        QuadratureOptions3 {
            depth: 1,
            ..Default::default()
        },
        &mut callback,
    )
    .unwrap();
    let op = AdaptiveElasticity3::build_with_material_field(
        HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(),
        tree,
        &Top,
        ElasticMaterialField3::Uniform(&material()),
        &|p| p[0] == 0.0,
        Default::default(),
        &mut control,
    )
    .unwrap();
    AdaptiveSolveSpace3::jacobi(op, 100_000_000)
}

fn check_objective_independent_admission<O: Sdf3Elasticity>(mut op: O, nodes: &[[f64; 3]]) {
    // Manufacture equilibrium loads through the actual native operator. The
    // two coordinate endpoints each have 0.09 axial strain, but the combined
    // unit-ball direction reaches 0.06 + sqrt(2)*0.03 > the source ceiling 0.1.
    let offset: Vec<_> = nodes.iter().flat_map(|p| [0.06 * p[0], 0.0, 0.0]).collect();
    let direction: Vec<_> = nodes.iter().flat_map(|p| [0.03 * p[0], 0.0, 0.0]).collect();
    let endpoint: Vec<_> = offset.iter().zip(&direction).map(|(a, b)| a + b).collect();
    op.bulk_stress(&endpoint, 100_000, poll).unwrap();
    let outside: Vec<_> = offset
        .iter()
        .zip(&direction)
        .map(|(a, b)| a + std::f64::consts::SQRT_2 * b)
        .collect();
    assert!(op.bulk_stress(&outside, 100_000, poll).is_err());
    let mut dead = vec![0.0; op.n()];
    let mut load = vec![0.0; op.n()];
    op.apply(&offset, &mut dead);
    op.apply(&direction, &mut load);
    let scales: Vec<_> = (0..op.cells()).map(|i| 0.4 + 0.01 * i as f64).collect();
    op.set_scales(&scales).unwrap();
    let mut study = CutDensityStudy3::new(op, 0.0, SimpParams::default());
    let rho = vec![1.0; study.cells()];
    let single = [&load[..]];
    let pair = [&load[..], &load[..]];
    let single_policy = AffineLoadSpace3::new(1, 1.0).unwrap();
    let pair_policy = AffineLoadSpace3::new(2, 1.0).unwrap();
    let mut callback = |_| poll();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let admitted = study
        .evaluate_robust_affine_load_space(
            &rho,
            &dead,
            &single,
            RobustOptions3::default(),
            &single_policy,
            &mut control,
        )
        .unwrap();
    close_vector(&admitted.offset_displacement, &offset);
    close_vector(&admitted.basis_displacements[0], &direction);
    assert_eq!(study.operator().scales(), scales);

    let limits = [ResponseLimit3 {
        response: BoundedResponse3::Displacement(&load),
        sense: ResponseLimitSense3::AtMost,
        bound: 100.0,
        scale: 100.0,
    }];
    let responses = AffineResponseLimits3 {
        limits: &limits,
        options: AffineResponseOptions3::default(),
    };
    for with_responses in [false, true] {
        let mut callback = |_| poll();
        let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
        let result = if with_responses {
            study.evaluate_robust_affine_response_limits(
                &rho,
                &dead,
                &pair,
                RobustOptions3::default(),
                &pair_policy,
                responses,
                &mut control,
            )
        } else {
            study.evaluate_robust_affine_load_space(
                &rho,
                &dead,
                &pair,
                RobustOptions3::default(),
                &pair_policy,
                &mut control,
            )
        };
        assert_eq!(
            result.unwrap_err(),
            EvaluationStop::Breakdown {
                stage: "sdf3-affine-material-domain"
            }
        );
        // One filter, one volume pullback, and three primals. No response
        // adjoints or objective-gradient pullbacks run after material refusal.
        assert_eq!(control.work().linear_solves, 5);
        assert_eq!(study.operator().scales(), scales);
    }

    let mut stop = |progress: SolveProgress| {
        if progress.stage == "sdf3-affine-material-domain" {
            ControlFlow::Break(())
        } else {
            poll()
        }
    };
    let mut control = SolveControl::new(SolveBudget::default(), &mut stop);
    assert_eq!(
        study
            .evaluate_robust_affine_load_space(
                &rho,
                &dead,
                &single,
                RobustOptions3::default(),
                &single_policy,
                &mut control,
            )
            .unwrap_err(),
        EvaluationStop::Cancelled
    );
    assert_eq!(study.operator().scales(), scales);
    let mut callback = |_| poll();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let replay = study
        .evaluate_robust_affine_load_space(
            &rho,
            &dead,
            &single,
            RobustOptions3::default(),
            &single_policy,
            &mut control,
        )
        .unwrap();
    close_vector(&replay.worst_displacement, &admitted.worst_displacement);
    assert_eq!(study.operator().scales(), scales);
}

#[test]
fn affine_compliance_and_response_studies_admit_the_complete_material_strain_ball() {
    let mut callback = |_| poll();
    let mut control =
        QuadratureControl3::new(QuadratureOptions3::default(), &mut callback).unwrap();
    let uniform = CutElasticity3::build_with_material_field(
        HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(),
        [2; 3],
        &Top,
        ElasticMaterialField3::Uniform(&material()),
        &|p| p[0] == 0.0,
        Default::default(),
        &mut control,
    )
    .unwrap();
    let nodes = uniform.nodes().to_vec();
    check_objective_independent_admission(uniform, &nodes);

    let tree = Octree3::uniform(1, 4, 4096).unwrap();
    let mark = *tree.leaves().iter().next().unwrap();
    let graded = tree.refined(&[mark], poll).unwrap();
    let adaptive = adaptive(&graded);
    let nodes = adaptive.elasticity().nodes().to_vec();
    check_objective_independent_admission(adaptive, &nodes);
}

#[test]
fn affine_refinement_refuses_material_escape_on_the_enriched_grid() {
    let coarse = Octree3::uniform(0, 4, 4096).unwrap();
    let fine = Octree3::uniform(1, 4, 4096).unwrap();
    let mut op = adaptive(&coarse);
    let scales = vec![0.37; op.cells()];
    op.set_scales(&scales).unwrap();
    let mut study = CutDensityStudy3::new(op, 0.0, SimpParams::default());
    let rho = vec![1.0; study.cells()];
    // Orthogonal to the coarse Q1 x modes. The fine grid resolves the same
    // pure physical load law; its displacement exceeds the retained ceiling.
    let polynomial = |p: [f64; 3]| [100.0 * (6.0 * p[0] * p[0] - 6.0 * p[0] + 1.0), 0.0, 0.0];
    let loads = [ReferenceLoad3::body(&polynomial)];
    let coarse_force = study
        .operator()
        .elasticity()
        .reference_load(loads[0], poll)
        .unwrap();
    assert!(coarse_force.iter().all(|f| f.abs() < 1e-10));
    let zero = vec![0.0; coarse_force.len()];
    let policy = AffineLoadSpace3::new(1, 1.0).unwrap();
    let mut callback = |_| poll();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let baseline = study
        .evaluate_robust_affine_load_space(
            &rho,
            &zero,
            &[&coarse_force],
            RobustOptions3::default(),
            &policy,
            &mut control,
        )
        .unwrap();
    assert!(baseline.worst_displacement.iter().all(|u| u.abs() < 1e-8));
    assert_eq!(study.operator().scales(), scales);
    let error = study
        .estimate_affine_load_space_enrichment(
            adaptive(&fine),
            &rho,
            ReferenceLoad3::default(),
            &loads,
            RobustOptions3::default(),
            &policy,
            AffineLoadSpaceRefinementOptions3::default(),
            &mut control,
        )
        .unwrap_err();
    assert_eq!(
        error,
        GoalRefinementError3::Evaluation(EvaluationStop::Breakdown {
            stage: "sdf3-affine-material-domain"
        })
    );
    assert_eq!(study.operator().scales(), scales);
}

#[test]
fn zero_centered_uncertainty_checks_combinations_beyond_valid_independent_cases() {
    use fs_topopt::pipeline::LoadCase;
    use fs_topopt::sdf3::robust::envelope::load_space::LoadSpace3;

    let tree = Octree3::uniform(1, 4, 4096).unwrap();
    let mut op = adaptive(&tree);
    let field: Vec<_> = op
        .elasticity()
        .nodes()
        .iter()
        .flat_map(|p| [0.08 * p[0], 0.0, 0.0])
        .collect();
    op.bulk_stress(&field, 100_000, poll).unwrap();
    let mut force = vec![0.0; op.n()];
    op.apply(&field, &mut force);
    let scales = vec![0.43; op.cells()];
    op.set_scales(&scales).unwrap();
    let mut study = CutDensityStudy3::new(op, 0.0, SimpParams::default());
    let rho = vec![1.0; study.cells()];
    let mut callback = |_| poll();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    // Independent experiments do not combine their loads. Both fields remain
    // valid, even though the larger continuous uncertainty set below does not.
    let independent = study
        .evaluate_robust(
            &rho,
            &[
                LoadCase {
                    force: &force,
                    weight: 0.5,
                },
                LoadCase {
                    force: &force,
                    weight: 0.5,
                },
            ],
            RobustOptions3::default(),
            &mut control,
        )
        .unwrap();
    for u in &independent.eroded.objective.displacements {
        close_vector(u, &field);
    }
    assert_eq!(study.operator().scales(), scales);
    let single = study
        .evaluate_robust_load_space(
            &rho,
            &[&force],
            RobustOptions3::default(),
            &LoadSpace3::new(1, 1.0, 0.02).unwrap(),
            &mut control,
        )
        .unwrap();
    close_vector(&single.worst_displacement, &field);
    assert_eq!(study.operator().scales(), scales);
    assert_eq!(
        study
            .evaluate_robust_load_space(
                &rho,
                &[&force, &force],
                RobustOptions3::default(),
                &LoadSpace3::new(2, 1.0, 0.02).unwrap(),
                &mut control,
            )
            .unwrap_err(),
        EvaluationStop::Breakdown {
            stage: "sdf3-affine-material-domain"
        }
    );
    assert_eq!(study.operator().scales(), scales);
}

#[test]
fn zero_centered_refinement_checks_the_combined_fine_material_domain() {
    use fs_topopt::pipeline::LoadCase;
    use fs_topopt::sdf3::robust::envelope::load_space::LoadSpace3;
    use fs_topopt::sdf3::robust::envelope::load_space::design::limits::LoadSpaceCap3;
    use fs_topopt::sdf3_goal::{GoalPreconditioner3, GoalRefinementOptions3};

    let coarse = Octree3::uniform(0, 4, 4096).unwrap();
    let fine = Octree3::uniform(1, 4, 4096).unwrap();
    let polynomial = |p: [f64; 3]| [6.0 * p[0] * p[0] - 6.0 * p[0] + 1.0, 0.0, 0.0];
    let mut fine_study = CutDensityStudy3::new(adaptive(&fine), 0.0, SimpParams::default());
    let fine_rho = vec![1.0; fine_study.cells()];
    let small_force = fine_study
        .operator()
        .elasticity()
        .reference_load(ReferenceLoad3::body(&polynomial), poll)
        .unwrap();
    let mut callback = |_| poll();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let small = fine_study
        .evaluate(
            &fine_rho,
            &[LoadCase {
                force: &small_force,
                weight: 1.0,
            }],
            &mut control,
        )
        .unwrap();
    // For this authored tensor C_Mandel=10 I, convert actual native physical
    // stresses back to Mandel strain. Calibrate inside the original domain;
    // no source ceiling is enlarged or inferred from an elastic modulus.
    let strain = fine_study
        .operator()
        .bulk_stress(&small.objective.displacements[0], 100_000, poll)
        .unwrap()
        .iter()
        .flat_map(|p| {
            p.reference_stress.iter().enumerate().map(|(i, stress)| {
                stress.abs() / 10.0 * if i < 3 { 1.0 } else { std::f64::consts::SQRT_2 }
            })
        })
        .fold(0.0_f64, f64::max);
    assert!(strain > 0.0 && strain < 0.1);
    let amplitude = 0.08 / strain;
    let load = |p: [f64; 3]| polynomial(p).map(|v| amplitude * v);
    let law = ReferenceLoad3::body(&load);
    let fine_force = fine_study
        .operator()
        .elasticity()
        .reference_load(law, poll)
        .unwrap();
    fine_study
        .evaluate_robust_load_space(
            &fine_rho,
            &[&fine_force],
            RobustOptions3::default(),
            &LoadSpace3::new(1, 1.0, 0.02).unwrap(),
            &mut control,
        )
        .unwrap();

    let mut op = adaptive(&coarse);
    let scales = vec![0.37; op.cells()];
    op.set_scales(&scales).unwrap();
    let mut study = CutDensityStudy3::new(op, 0.0, SimpParams::default());
    let rho = vec![1.0; study.cells()];
    let coarse_force = study
        .operator()
        .elasticity()
        .reference_load(law, poll)
        .unwrap();
    assert!(coarse_force.iter().all(|f| f.abs() < 1e-10));
    let params = study.params();
    let error = study
        .estimate_load_space_cap_enrichment(
            adaptive(&fine).into_elasticity(),
            &[law, law],
            &rho,
            params,
            RobustOptions3::default(),
            &LoadSpaceCap3::new(2, 1e6, 0.02).unwrap(),
            GoalRefinementOptions3 {
                preconditioner: GoalPreconditioner3::Jacobi {
                    max_contributions: 100_000_000,
                },
                ..Default::default()
            },
            &mut control,
        )
        .unwrap_err();
    assert_eq!(
        error,
        GoalRefinementError3::Evaluation(EvaluationStop::Breakdown {
            stage: "sdf3-affine-material-domain"
        })
    );
    let restored = study.params();
    assert_eq!(
        [restored.e_min, restored.penal, restored.beta, restored.eta],
        [params.e_min, params.penal, params.beta, params.eta]
    );
    assert_eq!(study.operator().scales(), scales);
}
