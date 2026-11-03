//! Original cohesive bodies across contact modes, with independent mechanics.
use super::*;
use crate::cohesive_contact::impact::TraceTerm;
use crate::cohesive_contact::impact::frictional::{
    FrictionalImpactControls, FrictionalResponseControls, NormalVelocityLaw,
};
use crate::cohesive_contact::trajectory::arrival::ArrivalControls;
use crate::cohesive_contact::trajectory::sustained::{StoppingControls, SustainedControls};
use fs_contact::frictional_multicontact::FrictionalContactConfig;
use fs_contact::sustained::solve::SolveControls;
use fs_material::cohesive::BilinearCohesive;
use fs_solid::cohesive::dynamics::DynamicsConfig;
use fs_solid::cohesive::dynamics::impulse::ImpulseControls;
use fs_solid::cohesive::{Budget, CohesivePatch, CohesiveSystem, Config, JumpTerm};
use fs_solver::{CheckedCgConfig, LinearOp, NewtonKrylovConfig};
use fs_sparse::precond::IdentityPrecond;

const START_S: f64 = 3.0;
const HELD_OPENING_M: f64 = 0.002;
const HELD_FORCE_N: f64 = 10.0 * HELD_OPENING_M + 0.01 * (0.004 - HELD_OPENING_M) / 0.003;

// Three free translations and the original independent damaged cohesive DOF.
struct HeldElastic;
impl LinearOp for HeldElastic {
    fn n(&self) -> usize {
        4
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        y.fill(0.0);
        y[3] = 10.0 * x[3];
    }
}
struct Mass;
impl LinearOp for Mass {
    fn n(&self) -> usize {
        4
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        for i in 0..4 {
            y[i] = 2.0 * x[i];
        }
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
fn dynamics(maximum_step_s: f64) -> DynamicsConfig {
    DynamicsConfig {
        solve: Config {
            newton: newton(),
            max_newton_iterations: 60,
            max_operator_evaluations: 4000,
            force_tolerance_n: 1e-9,
            local_energy_tolerance_j_m2: 1e-10,
            total_cohesive_energy_tolerance_j: 1e-10,
            maximum_displacement_m: 100.0,
        },
        difference_step_m: 1.0 / 4096.0,
        maximum_step_s,
        maximum_velocity_m_s: 100.0,
        maximum_steps: 1000,
        absolute_energy_tolerance_j: 1e-10,
        relative_energy_tolerance: 1e-8,
    }
}
fn controls(maximum_step_s: f64) -> EvolutionControls {
    let cg = CheckedCgConfig {
        tolerance: 1e-13,
        max_iterations: 64,
        max_corrections: 2,
    };
    EvolutionControls {
        trajectory: TrajectoryControls {
            impact: FrictionalImpactControls {
                response: FrictionalResponseControls {
                    mass_solve: cg,
                    contact_solve: FrictionalContactConfig::default(),
                    maximum_bodies: 4,
                    maximum_trace_terms: 64,
                    maximum_workspace_entries: 1_000_000,
                    velocity_tolerance_m_s: 1e-9,
                    energy_tolerance_j: 1e-8,
                },
                impulse: ImpulseControls {
                    solve: cg,
                    max_workspace_entries: 8192,
                    absolute_momentum_tolerance_ns: 1e-12,
                    relative_momentum_tolerance: 1e-12,
                    absolute_work_tolerance_j: 1e-12,
                    relative_work_tolerance: 1e-12,
                },
            },
            maximum_step_s,
            minimum_step_s: 1e-6,
            maximum_accepted_steps: 64,
            maximum_body_trials: 20_000,
            maximum_refinements_per_step: 64,
            maximum_impacts: 4,
            maximum_residual_evaluations: 10_000_000,
            maximum_newton_iterations: 1_000_000,
            maximum_workspace_entries: 20_000_000,
            closing_velocity_tolerance_m_s: 1e-8,
            energy_tolerance_j: 1e-7,
        },
        arrival: ArrivalControls {
            minimum_step_s: 1e-6,
            time_tolerance_s: 1e-10,
            gap_tolerance_m: 1e-10,
            maximum_trials: 64,
        },
        sustained: SegmentControls {
            step: SustainedControls {
                solve: SolveControls {
                    // The response is itself an original nonlinear solve.
                    // Resolve its attainable outer force accuracy while staying
                    // ten times tighter than the independent physical gate.
                    newton: NewtonKrylovConfig {
                        absolute_tolerance: 1e-9,
                        ..newton()
                    },
                    difference_step_n: 0.125,
                    force_residual_tolerance_n: 1e-8,
                    maximum_response_evaluations: 1000,
                    ..SolveControls::default()
                },
                audit: super::super::sustained::AuditControls {
                    force_tolerance_n: 1e-8,
                    gap_tolerance_m: 1e-10,
                    slip_tolerance_m: 1e-10,
                    velocity_tolerance_m_s: 1e-8,
                    energy_tolerance_j: 1e-8,
                },
                maximum_bodies: 4,
                maximum_trace_terms: 64,
                maximum_body_trials: 4000,
                maximum_residual_evaluations: 4_000_000,
                maximum_native_newton_iterations: 400_000,
                maximum_contact_force_n: 100.0,
                maximum_workspace_entries: 1_000_000,
            },
            stopping: StoppingControls {
                minimum_step_s: 1e-6,
                time_tolerance_s: 1e-10,
                velocity_tolerance_m_s: 1e-8,
                maximum_trials: 64,
            },
            maximum_contact_solves: 128,
            maximum_response_evaluations: 4096,
        },
        maximum_contact_solves: 1024,
        maximum_response_evaluations: 16384,
        maximum_stopping_events: 16,
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
    maximum_step_s: f64,
) -> CohesiveMotion<'m, 'a> {
    let state = system
        .bonded_state()
        .solve_controlled(
            &[0.0, 0.0, 0.0, HELD_FORCE_N],
            &[[0.0; 3]],
            dynamics(maximum_step_s).solve,
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
            dynamics(maximum_step_s),
        )
        .unwrap()
}
fn evolve<'m, 'a>(
    motion: &CohesiveMotion<'m, 'a>,
    force: [f64; 4],
    initial_force: [f64; 3],
    duration_s: f64,
    maximum_step_s: f64,
) -> EvolutionCandidate<'m, 'a> {
    let rows: [[TraceTerm; 1]; 3] = core::array::from_fn(|dof| {
        [TraceTerm {
            body: 0,
            dof,
            coefficient: 1.0,
        }]
    });
    let contact = FrictionalContact {
        id: 9,
        gap_m: 0.0,
        activation_tolerance_m: 1e-8,
        friction_coefficient: 0.25,
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.0 },
        prescribed_velocity_m_s: [0.0; 3],
        trace: [&rows[0], &rows[1], &rows[2]],
    };
    let offsets = [[0.0; 3]];
    let preconditioner = IdentityPrecond;
    advance_candidate_controlled(
        &[motion],
        &[contact],
        &[BodyLoad {
            mean_force_n: &force,
            endpoint_prescribed_jumps_m: &offsets,
        }],
        &[&preconditioner],
        &[NaturalMapScales {
            normal_n_per_m: 256.0,
            tangent_n_per_m: 256.0,
        }],
        &[initial_force],
        duration_s,
        controls(maximum_step_s),
        &mut || false,
    )
    .unwrap()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}
fn unchanged(actual: &CohesiveMotion<'_, '_>, before: &CohesiveMotion<'_, '_>) {
    assert_eq!(
        actual.state().displacements_m(),
        before.state().displacements_m()
    );
    assert_eq!(actual.state().points(), before.state().points());
    assert_eq!(actual.velocity_m_s(), before.velocity_m_s());
    assert_eq!(actual.time_s().to_bits(), before.time_s().to_bits());
    assert_eq!(actual.steps(), before.steps());
}
fn retained(after: &CohesiveMotion<'_, '_>, before: &CohesiveMotion<'_, '_>) {
    assert!(core::ptr::eq(
        after.state().system(),
        before.state().system()
    ));
    assert!(core::ptr::eq(after.mass_operator(), before.mass_operator()));
    close(
        after.state().displacements_m()[3],
        before.state().displacements_m()[3],
        1e-10,
    );
    close(
        after.state().points()[0].damage(),
        before.state().points()[0].damage(),
        1e-8,
    );
    assert!(after.state().points()[0].damage() >= before.state().points()[0].damage());
}

#[test]
fn native_sliding_stops_once_then_remains_supported_for_the_whole_window() {
    let elastic = HeldElastic;
    let mass = Mass;
    let law = BilinearCohesive::new(1000.0, 1.0, 0.002, 1000.0, 0.02).unwrap();
    let system = system(&elastic, &law);
    let before = motion(&system, &mass, [0.0, 2.5, 0.0], 1.5);
    let original = before.clone();
    let candidate = evolve(
        &before,
        [-20.0, 0.0, 0.0, HELD_FORCE_N],
        [20.0, -5.0, 0.0],
        2.0,
        1.5,
    );
    let after = &candidate.motions[0];
    let report = &candidate.report;
    close(after.velocity_m_s()[1], 0.0, 1e-8);
    close(after.state().displacements_m()[1], 1.25, 1e-9);
    assert_eq!(after.time_s().to_bits(), (START_S + 2.0).to_bits());
    assert_eq!(report.stopping_events.len(), 1);
    close(report.stopping_events[0].time_s, START_S + 1.0, 1e-9);
    assert_eq!(report.stopping_events[0].contact_ids, vec![9]);
    assert_eq!(report.sustained_steps, 2);
    assert_eq!(report.trajectory.accepted_steps, 2);
    assert!(report.trajectory.events.is_empty());
    close(report.contact_friction_loss_j, 6.25, 1e-8);
    close(report.trajectory.kinetic_change_j, -6.25, 1e-8);
    close(report.trajectory.external_work_j, 0.0, 1e-8);
    close(report.trajectory.internal_friction_loss_j, 0.0, 1e-12);
    close(report.trajectory.energy_residual_j, 0.0, 1e-8);
    close(report.warm_start_forces_n[0][0], 20.0, 1e-8);
    close(report.warm_start_forces_n[0][1], 0.0, 1e-8);
    retained(after, &before);
    unchanged(&before, &original);
}

#[test]
fn authored_force_drives_breakaway_in_the_opposite_direction_after_the_stop() {
    let elastic = HeldElastic;
    let mass = Mass;
    let law = BilinearCohesive::new(1000.0, 1.0, 0.002, 1000.0, 0.02).unwrap();
    let system = system(&elastic, &law);
    let before = motion(&system, &mass, [0.0, 3.0, 0.0], 0.75);
    let original = before.clone();
    let candidate = evolve(
        &before,
        [-20.0, -7.0, 0.0, HELD_FORCE_N],
        [20.0, -5.0, 0.0],
        1.0,
        0.75,
    );
    let after = &candidate.motions[0];
    let report = &candidate.report;
    // First half: a=-6, s=.75, D=3.75. Second half: a=-1,
    // s=-.125, D=.625. Force changes from -5 to +5 N after actual rest.
    close(after.velocity_m_s()[1], -0.5, 1e-8);
    close(after.state().displacements_m()[1], 0.625, 1e-9);
    assert_eq!(report.stopping_events.len(), 1);
    close(report.stopping_events[0].time_s, START_S + 0.5, 1e-9);
    close(report.warm_start_forces_n[0][1], 5.0, 1e-8);
    close(report.contact_friction_loss_j, 4.375, 1e-8);
    close(report.trajectory.external_work_j, -4.375, 1e-8);
    close(report.trajectory.kinetic_change_j, -8.75, 1e-8);
    close(report.trajectory.energy_residual_j, 0.0, 1e-8);
    assert_eq!(after.time_s().to_bits(), (START_S + 1.0).to_bits());
    retained(after, &before);
    unchanged(&before, &original);
}

#[test]
fn actual_zero_gap_plastic_impact_hands_off_to_loaded_rest_without_chatter() {
    let elastic = HeldElastic;
    let mass = Mass;
    let law = BilinearCohesive::new(1000.0, 1.0, 0.002, 1000.0, 0.02).unwrap();
    let system = system(&elastic, &law);
    let before = motion(&system, &mass, [-1.0, 0.0, 0.0], 0.125);
    let original = before.clone();
    let candidate = evolve(
        &before,
        [-20.0, 0.0, 0.0, HELD_FORCE_N],
        [0.0; 3],
        0.25,
        0.125,
    );
    let after = &candidate.motions[0];
    let report = &candidate.report;
    assert_eq!(report.trajectory.events.len(), 1);
    assert_eq!(
        report.trajectory.events[0].time_s.to_bits(),
        START_S.to_bits()
    );
    assert!(report.stopping_events.is_empty());
    assert_eq!(report.sustained_steps, 2);
    assert_eq!(after.time_s().to_bits(), (START_S + 0.25).to_bits());
    for &velocity in &after.velocity_m_s()[..3] {
        close(velocity, 0.0, 1e-8);
    }
    close(after.state().displacements_m()[0], 0.0, 1e-10);
    close(report.warm_start_forces_n[0][0], 20.0, 1e-8);
    close(report.trajectory.impact_loss_j, 1.0, 1e-8);
    close(report.contact_friction_loss_j, 0.0, 1e-10);
    close(report.trajectory.kinetic_change_j, -1.0, 1e-8);
    close(report.trajectory.energy_residual_j, 0.0, 1e-8);
    retained(after, &before);
    unchanged(&before, &original);
}

#[test]
fn unloaded_contact_releases_into_original_free_flight() {
    let elastic = HeldElastic;
    let mass = Mass;
    let law = BilinearCohesive::new(1000.0, 1.0, 0.002, 1000.0, 0.02).unwrap();
    let system = system(&elastic, &law);
    let before = motion(&system, &mass, [0.0; 3], 0.125);
    let original = before.clone();
    let candidate = evolve(
        &before,
        [20.0, 0.0, 0.0, HELD_FORCE_N],
        [0.0; 3],
        0.5,
        0.125,
    );
    let after = &candidate.motions[0];
    let report = &candidate.report;
    assert_eq!(report.sustained_steps, 1);
    assert_eq!(report.trajectory.accepted_steps, 4);
    assert!(report.trajectory.events.is_empty());
    assert!(report.stopping_events.is_empty());
    close(after.velocity_m_s()[0], 5.0, 1e-8);
    close(after.state().displacements_m()[0], 1.25, 1e-9);
    close(report.trajectory.final_gaps_m[0], 1.25, 1e-9);
    close(report.trajectory.external_work_j, 25.0, 1e-8);
    close(report.trajectory.kinetic_change_j, 25.0, 1e-8);
    close(report.contact_friction_loss_j, 0.0, 1e-10);
    close(report.trajectory.energy_residual_j, 0.0, 1e-8);
    retained(after, &before);
    unchanged(&before, &original);
}

#[path = "tests/arrival.rs"]
mod arrival;

#[path = "tests/transitions.rs"]
mod transitions;
