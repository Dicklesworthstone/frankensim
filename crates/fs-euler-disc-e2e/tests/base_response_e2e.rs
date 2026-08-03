#[path = "../src/base_response.rs"]
mod base_response;

use base_response::{
    BaseGeometryScope, BaseResponseError, BaseResponseInput, ContactLoadScope, LevelSupportInput,
    MovingContactLoad, refine_reduced_base_response, run_reduced_base_response,
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
    assert_eq!(run.diagnostics.level_tilt_rad, 0.0);
    let OperatorDiagnostics::Computed {
        mass_min_eigenvalue,
        stiffness_condition_number,
        symmetry_residual,
        ..
    } = &run.diagnostics.operator
    else {
        panic!("fixture remains below the production conditioning budget");
    };
    assert!(*mass_min_eigenvalue > 0.0);
    assert!(stiffness_condition_number.is_finite());
    assert!(*symmetry_residual <= f64::EPSILON);
    let terminal = run
        .final_sample()
        .expect("integrator retains initial sample");
    assert!(terminal.elastic_energy_j > 0.0);
    assert!(terminal.damping_work_j > 0.0);
    assert!(terminal.support_reaction_norm_n > 0.0);
    assert!(run.energy_closure_residual_j.abs() < 2.0e-6);
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
}

#[test]
fn e2e_reduced_flexible_base_timestep_refinement_is_retained() {
    let evidence = refine_reduced_base_response(&input(DampingModel::Rayleigh {
        mass_proportional_per_s: 0.1,
        stiffness_proportional_s: 5.0e-7,
    }))
    .expect("refinement pair");
    assert!(evidence.terminal_displacement_difference_m < 1.0e-6);
    assert!(evidence.terminal_elastic_energy_difference_j < 1.0e-6);
    assert!(evidence.coarse.energy_closure_residual_j.abs() < 2.0e-6);
    assert!(evidence.fine.energy_closure_residual_j.abs() < 2.0e-6);
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
