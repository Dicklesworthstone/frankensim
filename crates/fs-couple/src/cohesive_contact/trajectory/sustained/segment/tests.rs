//! G1 stopping and continuation through the original damaged cohesive bodies.
use crate::cohesive_contact::impact::TraceTerm;
use crate::cohesive_contact::trajectory::BodyLoad;
use crate::cohesive_contact::trajectory::sustained::{
    AuditControls, NaturalMapScales, StoppingControls, SustainedContact, SustainedControls,
};
use crate::cohesive_contact::trajectory::sustained::{
    SegmentControls, segment_candidate_controlled,
};
use fs_contact::sustained::solve::SolveControls;
use fs_material::cohesive::BilinearCohesive;
use fs_solid::cohesive::dynamics::{CohesiveMotion, DynamicsConfig};
use fs_solid::cohesive::{Budget, CohesivePatch, CohesiveSystem, Config, JumpTerm};
use fs_solver::{LinearOp, NewtonKrylovConfig};

const START: f64 = 3.0;
const REQUEST: f64 = 3.0 / 32.0;
const STOP: f64 = 1.0 / 16.0;
const HELD_OPENING: f64 = 0.002;
const HELD_FORCE: f64 = 10.0 * HELD_OPENING + 0.01 * (0.004 - HELD_OPENING) / 0.003;

// Three translations and the original independent damaged cohesive opening.
struct Elastic;
impl LinearOp for Elastic {
    fn n(&self) -> usize {
        4
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        y.fill(0.0);
        y[3] = 10.0 * x[3];
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
        maximum_step_s: 0.125,
        maximum_velocity_m_s: 20.0,
        maximum_steps: 1000,
        absolute_energy_tolerance_j: 1e-10,
        relative_energy_tolerance: 1e-8,
    }
}
fn controls() -> SegmentControls {
    SegmentControls {
        step: SustainedControls {
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
        },
        stopping: StoppingControls {
            minimum_step_s: 1e-4,
            time_tolerance_s: 1e-9,
            velocity_tolerance_m_s: 1e-8,
            maximum_trials: 64,
        },
        maximum_contact_solves: 128,
        maximum_response_evaluations: 5000,
    }
}
fn scale() -> NaturalMapScales {
    NaturalMapScales {
        normal_n_per_m: 1024.0,
        tangent_n_per_m: 1024.0,
    }
}
fn system<'a>(elastic: &'a Elastic, law: &'a BilinearCohesive) -> CohesiveSystem<'a> {
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
    tangent: f64,
) -> CohesiveMotion<'m, 'a> {
    let state = system
        .bonded_state()
        .solve_controlled(
            &[0.0, 0.0, 0.0, HELD_FORCE],
            &[[0.0; 3]],
            dynamics().solve,
            &mut || false,
        )
        .unwrap()
        .state;
    assert!(state.points()[0].damage() > 0.0 && state.points()[0].damage() < 1.0);
    close(state.displacements_m()[3], HELD_OPENING);
    state
        .into_dynamics(mass, vec![0.0, tangent, 0.0, 0.0], START, dynamics())
        .unwrap()
}
fn rows() -> [[TraceTerm; 1]; 3] {
    core::array::from_fn(|dof| {
        [TraceTerm {
            body: 0,
            dof,
            coefficient: 1.0,
        }]
    })
}
fn contact<'a>(trace: [&'a [TraceTerm]; 3], surface_speed: f64) -> SustainedContact<'a> {
    SustainedContact {
        id: 9,
        gap_m: 0.0,
        activation_tolerance_m: 1e-8,
        friction_coefficient: 0.25,
        prescribed_velocity_m_s: [0.0, -surface_speed, 0.0],
        trace,
    }
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-8,
        "{actual:e} != {expected:e}"
    );
}
fn retained(actual: &CohesiveMotion<'_, '_>, before: &CohesiveMotion<'_, '_>) {
    assert!(core::ptr::eq(
        actual.state().system(),
        before.state().system()
    ));
    assert!(core::ptr::eq(
        actual.mass_operator(),
        before.mass_operator()
    ));
    assert_eq!(actual.state().points(), before.state().points());
    close(
        actual.state().displacements_m()[3],
        before.state().displacements_m()[3],
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

#[test]
fn actual_sliding_stops_without_projection_then_remains_supported() {
    let elastic = Elastic;
    let mass = Mass(2.0);
    let law = BilinearCohesive::new(1000.0, 1.0, 0.002, 1000.0, 0.02).unwrap();
    let system = system(&elastic, &law);
    let before = motion(&system, &mass, 5.0 / 32.0);
    let original = before.clone();
    let rows = rows();
    let contacts = [contact([&rows[0], &rows[1], &rows[2]], 0.0)];
    let force = [-20.0, 0.0, 0.0, HELD_FORCE];
    let offsets = [[0.0; 3]];
    let loads = [BodyLoad {
        mean_force_n: &force,
        endpoint_prescribed_jumps_m: &offsets,
    }];
    let first = segment_candidate_controlled(
        &[&before],
        &contacts,
        &loads,
        &[scale()],
        &[[20.0, -5.0, 0.0]],
        REQUEST,
        controls(),
        &mut || false,
    )
    .unwrap();
    close(first.report.accepted_duration_s, STOP);
    assert_eq!(first.report.stopped_contact_ids, vec![9]);
    let event = first.report.stopping.as_ref().unwrap();
    close(event.elapsed_s, STOP);
    close(first.step.bodies[0].motion.velocity_m_s()[1], 0.0);
    close(
        first.step.bodies[0].motion.state().displacements_m()[1],
        5.0 / 1024.0,
    );
    close(first.step.report.audit.friction_loss_j, 25.0 / 1024.0);
    close(first.step.report.audit.energy_residual_j, 0.0);
    retained(&first.step.bodies[0].motion, &before);
    let second = segment_candidate_controlled(
        &[&first.step.bodies[0].motion],
        &contacts,
        &loads,
        &[scale()],
        &first.step.report.contact.forces_n,
        REQUEST - first.report.accepted_duration_s,
        controls(),
        &mut || false,
    )
    .unwrap();
    assert!(second.report.stopping.is_none());
    close(second.step.bodies[0].motion.time_s(), START + REQUEST);
    close(second.step.bodies[0].motion.velocity_m_s()[1], 0.0);
    close(
        second.step.bodies[0].motion.state().displacements_m()[1],
        5.0 / 1024.0,
    );
    close(second.step.report.contact.forces_n[0][1], 0.0);
    close(second.step.report.audit.friction_loss_j, 0.0);
    retained(&second.step.bodies[0].motion, &before);
    let exact = segment_candidate_controlled(
        &[&before],
        &contacts,
        &loads,
        &[scale()],
        &[[20.0, -5.0, 0.0]],
        STOP,
        controls(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(exact.report.stopped_contact_ids, vec![9]);
    assert_eq!(exact.report.contact_solves, 1);
    let endpoint = exact.report.stopping.unwrap();
    assert_eq!(endpoint.bracket_s, [STOP, STOP]);
    assert_eq!(endpoint.response_evaluations, 0);
    unchanged(&before, &original);
}

#[test]
fn a_moving_surface_stops_relative_slip_at_nonzero_world_speed() {
    let elastic = Elastic;
    let mass = Mass(2.0);
    let law = BilinearCohesive::new(1000.0, 1.0, 0.002, 1000.0, 0.02).unwrap();
    let system = system(&elastic, &law);
    let before = motion(&system, &mass, 27.0 / 32.0);
    let original = before.clone();
    let rows = rows();
    let contacts = [contact([&rows[0], &rows[1], &rows[2]], 1.0)];
    let force = [-20.0, 0.0, 0.0, HELD_FORCE];
    let offsets = [[0.0; 3]];
    let loads = [BodyLoad {
        mean_force_n: &force,
        endpoint_prescribed_jumps_m: &offsets,
    }];
    let first = segment_candidate_controlled(
        &[&before],
        &contacts,
        &loads,
        &[scale()],
        &[[20.0, 5.0, 0.0]],
        REQUEST,
        controls(),
        &mut || false,
    )
    .unwrap();
    close(first.report.accepted_duration_s, STOP);
    assert_eq!(first.report.stopped_contact_ids, vec![9]);
    close(first.step.bodies[0].motion.velocity_m_s()[1], 1.0);
    close(
        first.step.bodies[0].motion.state().displacements_m()[1],
        59.0 / 1024.0,
    );
    close(first.step.report.audit.surface_work_j, 320.0 / 1024.0);
    close(
        first.step.report.body_work[0].kinetic_change_j,
        295.0 / 1024.0,
    );
    close(first.step.report.audit.friction_loss_j, 25.0 / 1024.0);
    let second = segment_candidate_controlled(
        &[&first.step.bodies[0].motion],
        &contacts,
        &loads,
        &[scale()],
        &first.step.report.contact.forces_n,
        REQUEST - first.report.accepted_duration_s,
        controls(),
        &mut || false,
    )
    .unwrap();
    close(second.step.bodies[0].motion.velocity_m_s()[1], 1.0);
    close(second.step.bodies[0].motion.time_s(), START + REQUEST);
    close(second.step.report.audit.friction_loss_j, 0.0);
    retained(&second.step.bodies[0].motion, &before);
    unchanged(&before, &original);
}

#[test]
fn an_opposing_load_can_reverse_motion_only_after_the_resolved_stop() {
    let elastic = Elastic;
    let mass = Mass(2.0);
    let law = BilinearCohesive::new(1000.0, 1.0, 0.002, 1000.0, 0.02).unwrap();
    let system = system(&elastic, &law);
    let before = motion(&system, &mass, 15.0 / 32.0);
    let original = before.clone();
    let rows = rows();
    let contacts = [contact([&rows[0], &rows[1], &rows[2]], 0.0)];
    let force = [-20.0, -10.0, 0.0, HELD_FORCE];
    let offsets = [[0.0; 3]];
    let loads = [BodyLoad {
        mean_force_n: &force,
        endpoint_prescribed_jumps_m: &offsets,
    }];
    let first = segment_candidate_controlled(
        &[&before],
        &contacts,
        &loads,
        &[scale()],
        &[[20.0, -5.0, 0.0]],
        REQUEST,
        controls(),
        &mut || false,
    )
    .unwrap();
    close(first.report.accepted_duration_s, STOP);
    assert_eq!(first.report.stopped_contact_ids, vec![9]);
    close(first.step.bodies[0].motion.velocity_m_s()[1], 0.0);
    close(first.step.report.audit.friction_loss_j, 75.0 / 1024.0);
    close(
        first.step.report.body_work[0].authored_work_j,
        -150.0 / 1024.0,
    );
    let second = segment_candidate_controlled(
        &[&first.step.bodies[0].motion],
        &contacts,
        &loads,
        &[scale()],
        &first.step.report.contact.forces_n,
        REQUEST - first.report.accepted_duration_s,
        controls(),
        &mut || false,
    )
    .unwrap();
    close(second.step.bodies[0].motion.velocity_m_s()[1], -5.0 / 64.0);
    close(second.step.report.contact.forces_n[0][1], 5.0);
    close(second.step.report.audit.friction_loss_j, 25.0 / 4096.0);
    close(
        second.step.report.body_work[0].authored_work_j,
        50.0 / 4096.0,
    );
    close(second.step.report.audit.energy_residual_j, 0.0);
    retained(&second.step.bodies[0].motion, &before);
    unchanged(&before, &original);
}

#[test]
fn two_original_bodies_stop_relative_slip_and_preserve_total_momentum() {
    let elastic = Elastic;
    let masses = [Mass(1.0), Mass(2.0)];
    let law = BilinearCohesive::new(1000.0, 1.0, 0.002, 1000.0, 0.02).unwrap();
    let system = system(&elastic, &law);
    let left = motion(&system, &masses[0], 21.0 / 16.0);
    let right = motion(&system, &masses[1], 27.0 / 32.0);
    let originals = [left.clone(), right.clone()];
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
    let contacts = [contact([&rows[0], &rows[1], &rows[2]], 0.0)];
    let forces = [[20.0, 0.0, 0.0, HELD_FORCE], [-20.0, 0.0, 0.0, HELD_FORCE]];
    let offsets = [[0.0; 3]];
    let loads: [BodyLoad<'_>; 2] = core::array::from_fn(|i| BodyLoad {
        mean_force_n: &forces[i],
        endpoint_prescribed_jumps_m: &offsets,
    });
    let step = segment_candidate_controlled(
        &[&left, &right],
        &contacts,
        &loads,
        &[scale()],
        &[[20.0, 5.0, 0.0]],
        REQUEST,
        controls(),
        &mut || false,
    )
    .unwrap();
    close(step.report.accepted_duration_s, STOP);
    assert_eq!(step.report.stopped_contact_ids, vec![9]);
    close(step.step.bodies[0].motion.velocity_m_s()[1], 1.0);
    close(step.step.bodies[1].motion.velocity_m_s()[1], 1.0);
    close(
        step.step.bodies[0].motion.velocity_m_s()[1]
            + 2.0 * step.step.bodies[1].motion.velocity_m_s()[1],
        3.0,
    );
    close(step.step.report.audit.friction_loss_j, 75.0 / 1024.0);
    close(step.step.report.audit.surface_work_j, 0.0);
    close(step.step.report.audit.energy_residual_j, 0.0);
    for i in 0..2 {
        retained(&step.step.bodies[i].motion, &originals[i]);
    }
    unchanged(&left, &originals[0]);
    unchanged(&right, &originals[1]);
}
