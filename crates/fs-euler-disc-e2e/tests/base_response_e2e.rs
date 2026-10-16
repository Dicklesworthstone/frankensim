use fs_euler_disc_e2e::base_response::{
    ReducedBasePort, ReducedBasePortIdentity, ReducedBaseStepInput,
};
use fs_euler_disc_e2e::{
    BaseGeometryScope, BaseResponseError, BaseResponseInput, ContactLoadScope, LevelSupportInput,
    MAX_BASE_RESPONSE_STEPS, MovingContactLoad, refine_reduced_base_response,
    run_reduced_base_response,
};
use fs_solid::{
    AssemblyBudget, DampingModel, OperatorDiagnostics, ShellIdentity, ShellMaterial, ShellNode,
    ShellPlate, ShellSupport,
};

fn input(damping: DampingModel) -> BaseResponseInput {
    let support = ShellSupport {
        node_indices: [0, 1, 2],
        normal: [0.0, 0.0, 1.0],
    };
    BaseResponseInput {
        plate: ShellPlate {
            nodes: vec![
                ShellNode {
                    position_m: [-0.10, -0.08, 0.0],
                },
                ShellNode {
                    position_m: [0.10, -0.08, 0.0],
                },
                ShellNode {
                    position_m: [0.0, 0.12, 0.0],
                },
                ShellNode {
                    position_m: [0.0, 0.0, 0.0],
                },
            ],
            triangles: vec![[0, 1, 3], [1, 2, 3], [2, 0, 3]],
            thickness_m: 0.004,
            material: ShellMaterial {
                youngs_modulus_pa: 70.0e9,
                poisson_ratio: 0.33,
                density_kg_m3: 2_700.0,
            },
            identity: ShellIdentity {
                model_id: "e2e/base-response-v1".into(),
                source_id: "synthetic/flat-tripod-plate".into(),
                state_id: "initial".into(),
            },
            support: Some(support),
            damping,
            budget: AssemblyBudget::default(),
        },
        level_support: LevelSupportInput {
            support,
            level_normal: [0.0, 0.0, 1.0],
            maximum_tilt_rad: 1.0e-6,
        },
        geometry_scope: BaseGeometryScope::FlatSinglePatch,
        contact_scope: ContactLoadScope::NodalNormalLoad,
        load: MovingContactLoad {
            start_node: 3,
            end_node: 1,
            normal_force_n: 1.0,
        },
        initial_modal_displacement_m: 0.0,
        initial_modal_velocity_m_per_s: 0.0,
        timestep_s: 1.0e-6,
        steps: 200,
    }
}

fn independently_reconstructed_scaled_residual(
    request: &BaseResponseInput,
    static_shape: &[f64],
) -> f64 {
    let assembly = request.plate.assemble().expect("fixture plate assembles");
    let mut translation_scale_m = 0.0_f64;
    for (index, node) in request.plate.nodes.iter().enumerate() {
        for other in request.plate.nodes.iter().skip(index + 1) {
            let delta = [
                node.position_m[0] - other.position_m[0],
                node.position_m[1] - other.position_m[1],
                node.position_m[2] - other.position_m[2],
            ];
            translation_scale_m = translation_scale_m
                .max((delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2]).sqrt());
        }
    }
    let mut full_unit_load = vec![0.0; request.plate.nodes.len() * 6];
    for (node, weight) in [(request.load.start_node, 0.5), (request.load.end_node, 0.5)] {
        for component in 0..3 {
            full_unit_load[node * 6 + component] -=
                weight * request.level_support.level_normal[component];
        }
    }
    let reduced_shape: Vec<_> = assembly
        .free_dofs
        .iter()
        .map(|&dof| static_shape[dof])
        .collect();
    let reduced_load: Vec<_> = assembly
        .free_dofs
        .iter()
        .map(|&dof| full_unit_load[dof])
        .collect();
    let scales: Vec<_> = assembly
        .free_dofs
        .iter()
        .map(|dof| {
            if *dof % 6 < 3 {
                translation_scale_m
            } else {
                1.0
            }
        })
        .collect();
    let residual_squared: f64 = (0..assembly.stiffness.dimension())
        .map(|row| {
            let residual = assembly.stiffness.values()
                [row * assembly.stiffness.dimension()..(row + 1) * assembly.stiffness.dimension()]
                .iter()
                .zip(&reduced_shape)
                .map(|(coefficient, value)| coefficient * value)
                .sum::<f64>()
                - reduced_load[row];
            let scaled = scales[row] * residual;
            scaled * scaled
        })
        .sum();
    let load_squared: f64 = reduced_load
        .iter()
        .zip(&scales)
        .map(|(force, scale)| (force * scale) * (force * scale))
        .sum();
    residual_squared.sqrt() / load_squared.sqrt()
}

fn moving_input(damping: DampingModel) -> BaseResponseInput {
    let support = ShellSupport {
        node_indices: [0, 1, 2],
        normal: [0.0, 0.0, 1.0],
    };
    BaseResponseInput {
        plate: ShellPlate {
            nodes: vec![
                ShellNode {
                    position_m: [-0.12, -0.10, 0.0],
                },
                ShellNode {
                    position_m: [0.12, -0.10, 0.0],
                },
                ShellNode {
                    position_m: [0.0, 0.14, 0.0],
                },
                ShellNode {
                    position_m: [-0.025, -0.015, 0.0],
                },
                ShellNode {
                    position_m: [0.035, 0.025, 0.0],
                },
            ],
            triangles: vec![[0, 1, 3], [1, 4, 3], [1, 2, 4], [2, 0, 3], [2, 3, 4]],
            thickness_m: 0.004,
            material: ShellMaterial {
                youngs_modulus_pa: 70.0e9,
                poisson_ratio: 0.33,
                density_kg_m3: 2_700.0,
            },
            identity: ShellIdentity {
                model_id: "e2e/moving-base-response-v1".into(),
                source_id: "synthetic/flat-tripod-two-free-nodes".into(),
                state_id: "initial".into(),
            },
            support: Some(support),
            damping,
            budget: AssemblyBudget::default(),
        },
        level_support: LevelSupportInput {
            support,
            level_normal: [0.0, 0.0, 1.0],
            maximum_tilt_rad: 1.0e-6,
        },
        geometry_scope: BaseGeometryScope::FlatSinglePatch,
        contact_scope: ContactLoadScope::NodalNormalLoad,
        load: MovingContactLoad {
            start_node: 3,
            end_node: 4,
            normal_force_n: 1.0,
        },
        initial_modal_displacement_m: 0.0,
        initial_modal_velocity_m_per_s: 0.0,
        timestep_s: 1.0e-6,
        steps: 200,
    }
}

#[test]
fn e2e_reduced_flexible_base_retains_modal_energy_damping_work_support_reaction_and_conditioning() {
    let run = run_reduced_base_response(&input(DampingModel::Rayleigh {
        mass_proportional_per_s: 0.2,
        stiffness_proportional_s: 1.0e-6,
    }))
    .expect("flat production plate reduces");
    assert_eq!(run.samples.len(), 201);
    assert!(run.diagnostics.modal_mass_kg > 0.0);
    assert!(run.diagnostics.modal_stiffness_n_per_m > 0.0);
    assert!(run.diagnostics.modal_damping_n_s_per_m > 0.0);
    assert!(
        run.diagnostics.reduced_solve_scaled_residual
            <= run.diagnostics.reduced_solve_scaled_residual_limit
    );
    assert_eq!(run.diagnostics.level_tilt_rad, 0.0);
    let OperatorDiagnostics::Computed {
        raw_mass_min_eigenvalue,
        raw_stiffness_eigenvalue_spread,
        symmetry_residual,
        ..
    } = &run.diagnostics.operator
    else {
        panic!("fixture remains below the production conditioning budget");
    };
    assert!(*raw_mass_min_eigenvalue > 0.0);
    assert!(raw_stiffness_eigenvalue_spread.is_finite());
    assert!(*symmetry_residual <= f64::EPSILON);
    let terminal = run
        .final_sample()
        .expect("integrator retains initial sample");
    assert!(terminal.elastic_energy_j > 0.0);
    assert!(terminal.damping_work_j > 0.0);
    assert!(terminal.support_reaction_norm_n > 0.0);
    assert!(run.energy_closure_residual_j.abs() < 2.0e-6);
    assert!(run.normalized_energy_closure_residual < 1.0e-4);
    let reconstructed = independently_reconstructed_scaled_residual(
        &input(DampingModel::Rayleigh {
            mass_proportional_per_s: 0.2,
            stiffness_proportional_s: 1.0e-6,
        }),
        &run.diagnostics.supported_static_shape,
    );
    assert!(
        (reconstructed - run.diagnostics.reduced_solve_scaled_residual).abs() < 1.0e-12,
        "independent residual {reconstructed:e}"
    );
}

#[test]
fn e2e_reduced_flexible_base_moving_load_and_dimensioned_scaling_are_consistent() {
    let unit = moving_input(DampingModel::None);
    let run = run_reduced_base_response(&unit).expect("two-free-node moving response");
    assert_ne!(
        run.samples.first().expect("initial sample").modal_force_n,
        run.samples.last().expect("terminal sample").modal_force_n,
        "distinct free nodes exercise a moving projected load"
    );
    assert!(run.diagnostics.modal_shape_translation_scale_m > 0.0);
    assert!(run.energy_closure_residual_j.abs() < 2.0e-6);
    assert!(run.normalized_energy_closure_residual < 1.0e-4);

    let mut doubled = unit.clone();
    doubled.load.normal_force_n = 2.0;
    let doubled_run = run_reduced_base_response(&doubled).expect("scaled moving response");
    assert_eq!(
        run.diagnostics.modal_mass_kg,
        doubled_run.diagnostics.modal_mass_kg
    );
    assert_eq!(
        run.diagnostics.modal_stiffness_n_per_m,
        doubled_run.diagnostics.modal_stiffness_n_per_m
    );
    assert_eq!(
        run.diagnostics.modal_damping_n_s_per_m,
        doubled_run.diagnostics.modal_damping_n_s_per_m
    );
    assert!(
        (doubled_run
            .final_sample()
            .expect("terminal")
            .modal_displacement_m
            - 2.0 * run.final_sample().expect("terminal").modal_displacement_m)
            .abs()
            < 1.0e-12
    );
    assert!(
        (doubled_run
            .final_sample()
            .expect("terminal")
            .elastic_energy_j
            - 4.0 * run.final_sample().expect("terminal").elastic_energy_j)
            .abs()
            < 1.0e-16
    );
}

#[test]
fn e2e_reduced_flexible_base_zero_damping_preserves_free_modal_energy() {
    let mut request = input(DampingModel::None);
    request.load.normal_force_n = 0.0;
    request.initial_modal_velocity_m_per_s = 0.01;
    let run = run_reduced_base_response(&request).expect("undamped free response");
    assert_eq!(
        run.final_sample()
            .expect("integrator retains initial sample")
            .damping_work_j,
        0.0
    );
    assert!(run.energy_closure_residual_j.abs() < 1.0e-8);
    assert!(run.normalized_energy_closure_residual < 1.0e-5);
}

#[test]
fn e2e_reduced_flexible_base_timestep_refinement_is_retained() {
    let mut request = input(DampingModel::Rayleigh {
        mass_proportional_per_s: 0.1,
        stiffness_proportional_s: 5.0e-7,
    });
    request.steps = 100;
    let evidence = refine_reduced_base_response(&request).expect("refinement pair");
    assert!(evidence.terminal_displacement_difference_m < 1.0e-6);
    assert!(evidence.terminal_elastic_energy_difference_j < 1.0e-6);
    assert!(evidence.coarse.energy_closure_residual_j.abs() < 2.0e-6);
    assert!(evidence.fine.energy_closure_residual_j.abs() < 2.0e-6);
    assert!(evidence.coarse.normalized_energy_closure_residual < 1.0e-4);
    assert!(evidence.medium.normalized_energy_closure_residual < 1.0e-4);
    assert!(evidence.fine.normalized_energy_closure_residual < 1.0e-4);
    assert!(evidence.terminal_normalized_energy_difference < 1.0e-4);
    assert!(evidence.displacement_refinement_improved);
    assert!(evidence.energy_refinement_improved);
}

#[test]
fn e2e_public_base_response_api_refuses_unresolved_damping_and_long_synchronous_runs() {
    let unstable = input(DampingModel::Rayleigh {
        mass_proportional_per_s: 1.0e6,
        stiffness_proportional_s: 0.0,
    });
    assert!(matches!(
        run_reduced_base_response(&unstable),
        Err(BaseResponseError::TimestepOutsideResolution { .. })
    ));

    let mut too_long = input(DampingModel::None);
    too_long.steps = MAX_BASE_RESPONSE_STEPS + 1;
    assert!(matches!(
        run_reduced_base_response(&too_long),
        Err(BaseResponseError::StepBudgetExceeded)
    ));
}

#[test]
fn e2e_reduced_flexible_base_length_rescaling_preserves_scaled_solve_admission() {
    let original_request = input(DampingModel::None);
    let original = run_reduced_base_response(&original_request).expect("reference scale admits");
    let mut scaled_request = original_request.clone();
    let factor = 10.0;
    for node in &mut scaled_request.plate.nodes {
        for coordinate in &mut node.position_m {
            *coordinate *= factor;
        }
    }
    scaled_request.plate.thickness_m *= factor;
    scaled_request.timestep_s *= 0.1;
    let scaled = run_reduced_base_response(&scaled_request).expect("rescaled plate admits");
    for (run, request) in [(&original, &original_request), (&scaled, &scaled_request)] {
        assert!(run.diagnostics.reduced_solve_scaled_residual.is_finite());
        assert!(
            run.diagnostics.reduced_solve_scaled_residual
                <= run.diagnostics.reduced_solve_scaled_residual_limit
        );
        let reconstructed = independently_reconstructed_scaled_residual(
            request,
            &run.diagnostics.supported_static_shape,
        );
        assert!((reconstructed - run.diagnostics.reduced_solve_scaled_residual).abs() < 1.0e-12);
    }
}

#[test]
fn e2e_reduced_flexible_base_refuses_level_scope_and_resolved_contact_overreach() {
    let mut curved = input(DampingModel::None);
    curved.geometry_scope = BaseGeometryScope::CurvedShell;
    assert!(matches!(
        run_reduced_base_response(&curved),
        Err(BaseResponseError::UnsupportedScope {
            scope: "curved shell"
        })
    ));
    let mut multipatch = input(DampingModel::None);
    multipatch.geometry_scope = BaseGeometryScope::MultiPatch;
    assert!(matches!(
        run_reduced_base_response(&multipatch),
        Err(BaseResponseError::UnsupportedScope {
            scope: "multi-patch shell"
        })
    ));
    let mut resolved = input(DampingModel::None);
    resolved.contact_scope = ContactLoadScope::ResolvedFinitePatch;
    assert!(matches!(
        run_reduced_base_response(&resolved),
        Err(BaseResponseError::UnsupportedScope {
            scope: "resolved finite-patch contact"
        })
    ));
    let mut tilted = input(DampingModel::None);
    tilted.level_support.level_normal = [0.01, 0.0, (1.0 - 0.01_f64.powi(2)).sqrt()];
    assert!(matches!(
        run_reduced_base_response(&tilted),
        Err(BaseResponseError::SupportMismatch | BaseResponseError::InvalidInput { .. })
    ));
}

#[test]
fn e2e_reduced_base_port_replays_the_existing_implicit_midpoint_trajectory() {
    let request = moving_input(DampingModel::Rayleigh {
        mass_proportional_per_s: 0.2,
        stiffness_proportional_s: 1.0e-6,
    });
    let reference = run_reduced_base_response(&request).expect("reference trajectory");
    let port = ReducedBasePort::build(
        ReducedBasePortIdentity {
            model_id: "e2e/reduced-base-port-v1".into(),
            configuration_id: "moving-rayleigh-fixture-v1".into(),
        },
        request.clone(),
        u64::from(request.steps),
    )
    .expect("same flat nodal model prepares once");
    assert_eq!(
        port.diagnostics().modal_mass_kg,
        reference.diagnostics.modal_mass_kg,
        "the port reuses the prepared one-mode reduction"
    );

    let mut checkpoint = port.initial_checkpoint();
    let mut preceding_lineage_root = checkpoint.accepted_step_lineage_root();
    for step_index in 1..=request.steps {
        let proposal = port
            .propose(
                &checkpoint,
                &ReducedBaseStepInput {
                    // The version is part of the idempotency tuple, so a
                    // stable descriptive label does not require retaining an
                    // ever-growing set inside an audio-rate checkpoint.
                    step_id: "mechanics-substep".into(),
                    expected_version: checkpoint.accepted_version(),
                    duration_s: request.timestep_s,
                    compressive_normal_force_on_base_n: request.load.normal_force_n,
                    load_progress_start: f64::from(step_index - 1) / f64::from(request.steps),
                    load_progress_end: f64::from(step_index) / f64::from(request.steps),
                },
            )
            .expect("prepare interval");
        assert_eq!(proposal.receipt().parent_version, u64::from(step_index - 1));
        assert_eq!(proposal.receipt().next_version, u64::from(step_index));
        assert!(proposal.receipt().end_support_reaction_norm_n.is_finite());
        assert!(proposal.receipt().energy_closure_residual_j.abs() < 1.0e-14);
        if step_index == 1 {
            assert_eq!(
                port.refuse(&checkpoint, &proposal)
                    .expect("refuse leaves state intact"),
                checkpoint,
                "a refused mechanics candidate cannot advance base state"
            );
        }
        checkpoint = port.accept(&checkpoint, proposal).expect("commit interval");
        assert_ne!(
            checkpoint.accepted_step_lineage_root(),
            preceding_lineage_root
        );
        preceding_lineage_root = checkpoint.accepted_step_lineage_root();
    }
    let terminal = reference.final_sample().expect("terminal reference sample");
    assert_eq!(checkpoint.accepted_version(), u64::from(request.steps));
    for (label, actual, expected) in [
        ("elapsed time", checkpoint.elapsed_time_s(), terminal.time_s),
        (
            "modal displacement",
            checkpoint.modal_displacement_m(),
            terminal.modal_displacement_m,
        ),
        (
            "modal velocity",
            checkpoint.modal_velocity_m_per_s(),
            terminal.modal_velocity_m_per_s,
        ),
        (
            "damping work",
            checkpoint.cumulative_damping_work_j(),
            terminal.damping_work_j,
        ),
        (
            "external work",
            checkpoint.cumulative_external_work_j(),
            terminal.external_work_j,
        ),
    ] {
        assert!(
            (actual - expected).abs() <= 2.0e-15,
            "{label}: port={actual:e}, legacy={expected:e}"
        );
    }
}

#[test]
fn e2e_reduced_base_port_refuses_stale_replay_capacity_and_higher_fidelity_claims() {
    let request = input(DampingModel::None);
    let identity = ReducedBasePortIdentity {
        model_id: "e2e/reduced-base-port-v1".into(),
        configuration_id: "refusal-fixture-v1".into(),
    };
    let port = ReducedBasePort::build(identity, request.clone(), 1).expect("one interval budget");
    let checkpoint = port.initial_checkpoint();
    let accepted = port
        .propose(
            &checkpoint,
            &ReducedBaseStepInput {
                step_id: "accepted-step".into(),
                expected_version: 0,
                duration_s: 0.5 * request.timestep_s,
                compressive_normal_force_on_base_n: 1.0,
                load_progress_start: 0.0,
                load_progress_end: 0.01,
            },
        )
        .expect("first proposal");
    assert_eq!(accepted.receipt().compressive_normal_force_on_base_n, 1.0);
    assert_eq!(
        accepted.receipt().normal_reaction_on_disc_world_n,
        [0.0, 0.0, 1.0]
    );
    assert_eq!(accepted.receipt().timestep_s, 0.5 * request.timestep_s);
    let checkpoint = port.accept(&checkpoint, accepted).expect("first accept");
    assert!(matches!(
        port.propose(
            &checkpoint,
            &ReducedBaseStepInput {
                step_id: "stale-step".into(),
                expected_version: 0,
                duration_s: request.timestep_s,
                compressive_normal_force_on_base_n: 1.0,
                load_progress_start: 0.01,
                load_progress_end: 0.02,
            },
        ),
        Err(BaseResponseError::PortVersionMismatch { .. })
    ));
    assert!(matches!(
        port.propose(
            &checkpoint,
            &ReducedBaseStepInput {
                step_id: "over-budget-step".into(),
                expected_version: 1,
                duration_s: request.timestep_s,
                compressive_normal_force_on_base_n: 1.0,
                load_progress_start: 0.01,
                load_progress_end: 0.02,
            },
        ),
        Err(BaseResponseError::PortStepBudgetExceeded)
    ));
    let mut resolved = request;
    resolved.contact_scope = ContactLoadScope::ResolvedFinitePatch;
    assert!(matches!(
        ReducedBasePort::build(
            ReducedBasePortIdentity {
                model_id: "e2e/reduced-base-port-v1".into(),
                configuration_id: "resolved-finite-patch-v1".into(),
            },
            resolved,
            1,
        ),
        Err(BaseResponseError::UnsupportedScope {
            scope: "resolved finite-patch contact"
        })
    ));
    let mut as_built = input(DampingModel::None);
    as_built.geometry_scope = BaseGeometryScope::AsBuiltSurface;
    assert!(matches!(
        ReducedBasePort::build(
            ReducedBasePortIdentity {
                model_id: "e2e/reduced-base-port-v1".into(),
                configuration_id: "as-built-base-v1".into(),
            },
            as_built,
            1,
        ),
        Err(BaseResponseError::UnsupportedScope {
            scope: "as-built base surface"
        })
    ));
}
