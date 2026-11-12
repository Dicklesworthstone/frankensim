//! G1 physical fixtures through the original damaged-body midpoint owner.
use super::{
    AuditControls, NaturalMapScales, SustainedContact, SustainedControls, step_candidate_controlled,
};
use crate::cohesive_contact::impact::TraceTerm;
use crate::cohesive_contact::trajectory::BodyLoad;
use fs_contact::sustained::solve::SolveControls;
use fs_material::cohesive::BilinearCohesive;
use fs_solid::cohesive::dynamics::{CohesiveMotion, DynamicsConfig};
use fs_solid::cohesive::{Budget, CohesivePatch, CohesiveSystem, Config, JumpTerm};
use fs_solver::{LinearOp, NewtonKrylovConfig};

const START_S: f64 = 3.0;
const STEP_S: f64 = 0.125;
const HELD_OPENING_M: f64 = 0.002;
const HELD_FORCE_N: f64 = 10.0 * HELD_OPENING_M + 0.01 * (0.004 - HELD_OPENING_M) / 0.003;

/// Original four-DOF trajectory fixture: free translations and a held,
/// already-damaged cohesive opening. Contact never replaces the fourth DOF.
struct HeldElastic;
impl LinearOp for HeldElastic {
    fn n(&self) -> usize {
        4
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        y.fill(0.0);
        y[3] = 10.0 * x[3];
    }
    // Diagonal, so its transpose is itself.
    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
        self.apply(x, y);
    }
}
struct Mass(f64);
impl LinearOp for Mass {
    fn n(&self) -> usize {
        4
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        for i in 0..4 {
            y[i] = self.0 * x[i];
        }
    }
    // Diagonal, so its transpose is itself.
    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
        self.apply(x, y);
    }
}
fn newton() -> NewtonKrylovConfig {
    NewtonKrylovConfig {
        absolute_tolerance: 1e-11,
        relative_tolerance: 1e-12,
        linear_restart: 12,
        forcing_maximum: 0.05,
        ..NewtonKrylovConfig::default()
    }
}
fn dynamics() -> DynamicsConfig {
    DynamicsConfig {
        solve: Config {
            newton: newton(),
            max_newton_iterations: 60,
            max_operator_evaluations: 4000,
            force_tolerance_n: 1e-9,
            local_energy_tolerance_j_m2: 1e-10,
            total_cohesive_energy_tolerance_j: 1e-10,
            maximum_displacement_m: 2.0,
        },
        difference_step_m: 1.0 / 4096.0,
        maximum_step_s: STEP_S,
        maximum_velocity_m_s: 20.0,
        maximum_steps: 1000,
        absolute_energy_tolerance_j: 1e-10,
        relative_energy_tolerance: 1e-8,
    }
}
fn controls() -> SustainedControls {
    SustainedControls {
        solve: SolveControls {
            newton: newton(),
            difference_step_n: 0.125,
            force_residual_tolerance_n: 1e-8,
            maximum_response_evaluations: 1000,
            ..SolveControls::default()
        },
        audit: AuditControls {
            force_tolerance_n: 1e-8,
            gap_tolerance_m: 1e-10,
            slip_tolerance_m: 1e-10,
            velocity_tolerance_m_s: 1e-8,
            energy_tolerance_j: 1e-8,
        },
        maximum_bodies: 4,
        maximum_trace_terms: 64,
        maximum_body_trials: 4000,
        maximum_residual_evaluations: 1_000_000,
        maximum_native_newton_iterations: 100_000,
        maximum_contact_force_n: 100.0,
        maximum_workspace_entries: 1_000_000,
    }
}
fn scale() -> NaturalMapScales {
    NaturalMapScales {
        normal_n_per_m: 256.0,
        tangent_n_per_m: 256.0,
    }
}
fn system<'a>(elastic: &'a HeldElastic, law: &'a BilinearCohesive) -> CohesiveSystem<'a> {
    CohesiveSystem::new(
        elastic,
        vec![
            CohesivePatch::new(
                41,
                0.01,
                law,
                vec![JumpTerm {
                    dof: 3,
                    local: [1.0, 0.0, 0.0],
                }],
            )
            .unwrap(),
        ],
        Budget {
            max_dofs: 4,
            max_patches: 1,
            max_jump_terms: 1,
            max_krylov_basis_entries: 1000,
        },
    )
    .unwrap()
}
fn motion<'m, 'a>(
    system: &'m CohesiveSystem<'a>,
    mass: &'m Mass,
    velocity: [f64; 3],
) -> CohesiveMotion<'m, 'a> {
    let state = system
        .bonded_state()
        .solve_controlled(
            &[0.0, 0.0, 0.0, HELD_FORCE_N],
            &[[0.0; 3]],
            dynamics().solve,
            &mut || false,
        )
        .unwrap()
        .state;
    assert!(state.points()[0].damage() > 0.0 && state.points()[0].damage() < 1.0);
    close(state.displacements_m()[3], HELD_OPENING_M, 1e-11);
    state
        .into_dynamics(
            mass,
            vec![velocity[0], velocity[1], velocity[2], 0.0],
            START_S,
            dynamics(),
        )
        .unwrap()
}
fn wall_trace() -> [[TraceTerm; 1]; 3] {
    core::array::from_fn(|dof| {
        [TraceTerm {
            body: 0,
            dof,
            coefficient: 1.0,
        }]
    })
}
fn contact<'t>(trace: [&'t [TraceTerm]; 3], prescribed: [f64; 3]) -> SustainedContact<'t> {
    SustainedContact {
        id: 9,
        gap_m: 0.0,
        activation_tolerance_m: 1e-8,
        friction_coefficient: 0.25,
        prescribed_velocity_m_s: prescribed,
        trace,
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}
fn unchanged(actual: &CohesiveMotion<'_, '_>, original: &CohesiveMotion<'_, '_>) {
    assert_eq!(actual.velocity_m_s(), original.velocity_m_s());
    assert_eq!(
        actual.state().displacements_m(),
        original.state().displacements_m()
    );
    assert_eq!(actual.state().points(), original.state().points());
    assert_eq!(actual.time_s().to_bits(), original.time_s().to_bits());
    assert_eq!(actual.steps(), original.steps());
}
fn retained_history(actual: &CohesiveMotion<'_, '_>, original: &CohesiveMotion<'_, '_>) {
    assert!(core::ptr::eq(
        actual.state().system(),
        original.state().system()
    ));
    assert!(core::ptr::eq(
        actual.mass_operator(),
        original.mass_operator()
    ));
    close(
        actual.state().displacements_m()[3],
        original.state().displacements_m()[3],
        1e-10,
    );
    close(
        actual.state().points()[0].damage(),
        original.state().points()[0].damage(),
        1e-8,
    );
    assert!(actual.state().points()[0].damage() >= original.state().points()[0].damage());
}

#[test]
fn already_damaged_bodies_remain_supported_with_static_friction_across_steps() {
    let elastic = HeldElastic;
    let masses = [Mass(1.0), Mass(2.0)];
    let law = BilinearCohesive::new(1000.0, 1.0, 0.002, 1000.0, 0.02).unwrap();
    let system = system(&elastic, &law);
    let left = motion(&system, &masses[0], [0.0; 3]);
    let right = motion(&system, &masses[1], [0.0; 3]);
    let original = [left.clone(), right.clone()];
    let rows: [[TraceTerm; 2]; 3] = core::array::from_fn(|dof| {
        [
            TraceTerm {
                body: 0,
                dof,
                coefficient: -1.0,
            },
            TraceTerm {
                body: 1,
                dof,
                coefficient: 1.0,
            },
        ]
    });
    let contacts = [contact([&rows[0], &rows[1], &rows[2]], [0.0; 3])];
    let forces = [
        [20.0, 3.0, 0.0, HELD_FORCE_N],
        [-20.0, -3.0, 0.0, HELD_FORCE_N],
    ];
    let offsets = [[0.0; 3]];
    let loads: [BodyLoad<'_>; 2] = core::array::from_fn(|body| BodyLoad {
        mean_force_n: &forces[body],
        endpoint_prescribed_jumps_m: &offsets,
    });
    let first = step_candidate_controlled(
        &[&left, &right],
        &contacts,
        &loads,
        &[scale()],
        &[[20.0, 3.0, 0.0]],
        STEP_S,
        controls(),
        &mut || false,
    )
    .unwrap();
    let second = step_candidate_controlled(
        &[&first.bodies[0].motion, &first.bodies[1].motion],
        &contacts,
        &loads,
        &[scale()],
        &first.report.contact.forces_n,
        STEP_S,
        controls(),
        &mut || false,
    )
    .unwrap();
    for (step, elapsed) in [(&first, STEP_S), (&second, 2.0 * STEP_S)] {
        assert_eq!(step.report.contact_ids, vec![9]);
        assert_eq!(step.report.contact.forces_n, vec![[20.0, 3.0, 0.0]]);
        assert_eq!(step.report.contact.paths_m, vec![[0.0; 3]]);
        close(step.report.audit.friction_loss_j, 0.0, 1e-12);
        close(step.report.audit.normal_constraint_work_j, 0.0, 1e-12);
        close(step.report.audit.energy_residual_j, 0.0, 1e-12);
        for body in 0..2 {
            let actual = &step.bodies[body].motion;
            retained_history(actual, &original[body]);
            assert_eq!(actual.velocity_m_s(), original[body].velocity_m_s());
            assert_eq!(
                actual.state().displacements_m(),
                original[body].state().displacements_m()
            );
            assert_eq!(actual.state().points(), original[body].state().points());
            assert_eq!(actual.time_s().to_bits(), (START_S + elapsed).to_bits());
            assert_eq!(actual.steps(), (elapsed / STEP_S) as usize);
        }
    }
    unchanged(&left, &original[0]);
    unchanged(&right, &original[1]);
}

#[test]
fn cold_contact_solve_decelerates_actual_native_sliding_and_accounts_for_work() {
    let elastic = HeldElastic;
    let mass = Mass(2.0);
    let law = BilinearCohesive::new(1000.0, 1.0, 0.002, 1000.0, 0.02).unwrap();
    let system = system(&elastic, &law);
    let before = motion(&system, &mass, [0.0, 3.0, 0.0]);
    let original = before.clone();
    let rows = wall_trace();
    let force = [-20.0, 0.0, 0.0, HELD_FORCE_N];
    let offsets = [[0.0; 3]];
    let candidate = step_candidate_controlled(
        &[&before],
        &[contact([&rows[0], &rows[1], &rows[2]], [0.0; 3])],
        &[BodyLoad {
            mean_force_n: &force,
            endpoint_prescribed_jumps_m: &offsets,
        }],
        &[scale()],
        &[[0.0; 3]],
        STEP_S,
        controls(),
        &mut || false,
    )
    .unwrap();
    let after = &candidate.bodies[0].motion;
    let report = &candidate.report;
    assert!(report.contact.nonlinear.iterations > 0);
    for (actual, expected) in report.contact.forces_n[0].iter().zip([20.0, -5.0, 0.0]) {
        close(*actual, expected, 1e-8);
    }
    assert!(report.contact.paths_m[0][0] >= 0.0);
    assert!(after.state().displacements_m()[0] >= 0.0);
    close(after.velocity_m_s()[0], 0.0, 1e-8);
    // m=2, h=1/8, N=20, mu=1/4: dv=-5/16, s=91/256.
    close(after.velocity_m_s()[1], 43.0 / 16.0, 1e-8);
    close(after.state().displacements_m()[1], 91.0 / 256.0, 1e-9);
    close(report.audit.friction_loss_j, 455.0 / 256.0, 1e-8);
    close(report.body_work[0].kinetic_change_j, -455.0 / 256.0, 1e-8);
    close(report.audit.contact_body_work_j, -455.0 / 256.0, 1e-8);
    close(report.body_work[0].authored_work_j, 0.0, 1e-8);
    close(report.audit.surface_work_j, 0.0, 1e-12);
    close(report.audit.energy_residual_j, 0.0, 1e-8);
    close(report.audit.contact_work_residual_j, 0.0, 1e-8);
    assert_eq!(after.time_s().to_bits(), (START_S + STEP_S).to_bits());
    assert_eq!(after.steps(), before.steps() + 1);
    retained_history(after, &before);
    unchanged(&before, &original);
}

#[test]
fn moving_surface_supplies_native_kinetic_energy_and_relative_friction_loss() {
    let elastic = HeldElastic;
    let mass = Mass(2.0);
    let law = BilinearCohesive::new(1000.0, 1.0, 0.002, 1000.0, 0.02).unwrap();
    let system = system(&elastic, &law);
    let before = motion(&system, &mass, [0.0; 3]);
    let original = before.clone();
    let rows = wall_trace();
    let force = [-20.0, 0.0, 0.0, HELD_FORCE_N];
    let offsets = [[0.0; 3]];
    let candidate = step_candidate_controlled(
        &[&before],
        &[contact([&rows[0], &rows[1], &rows[2]], [0.0, -2.0, 0.0])],
        &[BodyLoad {
            mean_force_n: &force,
            endpoint_prescribed_jumps_m: &offsets,
        }],
        &[scale()],
        &[[20.0, 5.0, 0.0]],
        STEP_S,
        controls(),
        &mut || false,
    )
    .unwrap();
    let after = &candidate.bodies[0].motion;
    let report = &candidate.report;
    // Wall speed 2 m/s: body moves 5/256 m while the wall moves 1/4 m.
    close(after.velocity_m_s()[1], 5.0 / 16.0, 1e-8);
    close(after.state().displacements_m()[1], 5.0 / 256.0, 1e-9);
    close(report.contact.paths_m[0][1], -59.0 / 256.0, 1e-9);
    close(report.audit.surface_work_j, 5.0 / 4.0, 1e-8);
    close(report.body_work[0].kinetic_change_j, 25.0 / 256.0, 1e-8);
    close(report.audit.contact_body_work_j, 25.0 / 256.0, 1e-8);
    close(report.audit.friction_loss_j, 295.0 / 256.0, 1e-8);
    close(report.audit.normal_constraint_work_j, 0.0, 1e-10);
    close(report.audit.energy_residual_j, 0.0, 1e-8);
    retained_history(after, &before);
    unchanged(&before, &original);
}
