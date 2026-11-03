//! Real cohesive midpoint histories, timed Coulomb impacts, and atomic trajectory refusal.
use fs_contact::frictional_multicontact::FrictionalContactConfig;
use fs_couple::cohesive_contact::impact::TraceTerm;
use fs_couple::cohesive_contact::impact::frictional::{
    FrictionalContact, FrictionalImpactControls, FrictionalResponseControls, NormalVelocityLaw,
};
use fs_couple::cohesive_contact::trajectory::{
    BodyLoad, TrajectoryControls, TrajectoryError, advance_candidate_controlled,
};
use fs_material::cohesive::BilinearCohesive;
use fs_solid::cohesive::dynamics::impulse::ImpulseControls;
use fs_solid::cohesive::dynamics::{CohesiveMotion, DynamicsConfig};
use fs_solid::cohesive::{Budget, CohesivePatch, CohesiveSystem, Config, JumpTerm};
use fs_solver::{CheckedCgConfig, LinearOp, NewtonKrylovConfig};
use fs_sparse::precond::{IdentityPrecond, Precond};

const INITIAL_TIME_S: f64 = 3.0;
const DURATION_S: f64 = 0.01;
const HELD_DISPLACEMENT_M: f64 = 0.002;
const HELD_FORCE_N: f64 = 10.0 * HELD_DISPLACEMENT_M + 0.01 * (0.004 - HELD_DISPLACEMENT_M) / 0.003;

/// Three free translations and an independent, held cohesive opening.
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
fn dynamics() -> DynamicsConfig {
    DynamicsConfig {
        solve: Config {
            newton: NewtonKrylovConfig {
                absolute_tolerance: 1e-11,
                relative_tolerance: 1e-12,
                linear_restart: 12,
                forcing_maximum: 0.05,
                ..NewtonKrylovConfig::default()
            },
            max_newton_iterations: 60,
            max_operator_evaluations: 4000,
            force_tolerance_n: 1e-9,
            local_energy_tolerance_j_m2: 1e-10,
            total_cohesive_energy_tolerance_j: 1e-10,
            maximum_displacement_m: 0.04,
        },
        difference_step_m: 1e-8,
        maximum_step_s: DURATION_S,
        maximum_velocity_m_s: 4.0,
        maximum_steps: 1000,
        absolute_energy_tolerance_j: 1e-10,
        relative_energy_tolerance: 1e-8,
    }
}
fn controls() -> TrajectoryControls {
    let solve = CheckedCgConfig {
        tolerance: 1e-13,
        max_iterations: 32,
        max_corrections: 2,
    };
    TrajectoryControls {
        impact: FrictionalImpactControls {
            response: FrictionalResponseControls {
                mass_solve: solve,
                contact_solve: FrictionalContactConfig::default(),
                maximum_bodies: 4,
                maximum_trace_terms: 64,
                maximum_workspace_entries: 1_000_000,
                velocity_tolerance_m_s: 1e-9,
                energy_tolerance_j: 1e-9,
            },
            impulse: ImpulseControls {
                solve,
                max_workspace_entries: 8192,
                absolute_momentum_tolerance_ns: 1e-12,
                relative_momentum_tolerance: 1e-12,
                absolute_work_tolerance_j: 1e-12,
                relative_work_tolerance: 1e-12,
            },
        },
        maximum_step_s: DURATION_S,
        // Keep native displacement-difference probes inside their velocity bound.
        minimum_step_s: 1e-6,
        maximum_accepted_steps: 64,
        maximum_body_trials: 1000,
        maximum_refinements_per_step: 64,
        maximum_impacts: 4,
        maximum_residual_evaluations: 1_000_000,
        maximum_newton_iterations: 100_000,
        maximum_workspace_entries: 20_000_000,
        closing_velocity_tolerance_m_s: 1e-8,
        energy_tolerance_j: 1e-8,
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
    close(state.displacements_m()[3], HELD_DISPLACEMENT_M, 1e-11);
    state
        .into_dynamics(
            mass,
            vec![velocity[0], velocity[1], velocity[2], 0.0],
            INITIAL_TIME_S,
            dynamics(),
        )
        .unwrap()
}
fn pair_trace() -> [[TraceTerm; 2]; 3] {
    core::array::from_fn(|dof| {
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
    })
}
fn contact(rows: &[[TraceTerm; 2]; 3], gap_m: f64, mu: f64) -> FrictionalContact<'_> {
    FrictionalContact {
        id: 9,
        gap_m,
        activation_tolerance_m: 1e-8,
        friction_coefficient: mu,
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.5 },
        prescribed_velocity_m_s: [0.0; 3],
        trace: [&rows[0], &rows[1], &rows[2]],
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}
fn unchanged(actual: &CohesiveMotion<'_, '_>, expected: &CohesiveMotion<'_, '_>) {
    assert_eq!(actual.velocity_m_s(), expected.velocity_m_s());
    assert_eq!(
        actual.state().displacements_m(),
        expected.state().displacements_m()
    );
    assert_eq!(actual.state().points(), expected.state().points());
    assert_eq!(actual.time_s().to_bits(), expected.time_s().to_bits());
    assert_eq!(actual.steps(), expected.steps());
}

#[test]
fn damaged_bodies_fly_collide_and_continue_at_their_physical_clock() {
    let elastic = HeldElastic;
    let masses = [Mass(1.0), Mass(2.0)];
    let law = BilinearCohesive::new(1000.0, 1.0, 0.002, 1000.0, 0.02).unwrap();
    let system = system(&elastic, &law);
    let left = motion(&system, &masses[0], [1.0, 0.5, 0.0]);
    let right = motion(&system, &masses[1], [-0.5, -0.25, 0.0]);
    let original = [left.clone(), right.clone()];
    let rows = pair_trace();
    let force = [0.0, 0.0, 0.0, HELD_FORCE_N];
    let offsets = [[0.0; 3]];
    let loads: [BodyLoad<'_>; 2] = core::array::from_fn(|_| BodyLoad {
        mean_force_n: &force,
        endpoint_prescribed_jumps_m: &offsets,
    });
    let preconditioners: [&dyn Precond; 2] = [&IdentityPrecond, &IdentityPrecond];
    for (mu, tangential_impulse, loss) in [(0.0, 0.0, 0.5625), (0.25, 0.375, 0.73828125)] {
        let candidate = advance_candidate_controlled(
            &[&left, &right],
            &[contact(&rows, 0.006, mu)],
            &loads,
            &preconditioners,
            DURATION_S,
            controls(),
            &mut || false,
        )
        .unwrap();
        let report = &candidate.report;
        assert_eq!(report.events.len(), 1);
        assert!(report.accepted_steps >= 2);
        assert!(
            report.rejected_trials >= 1,
            "the first complete trial crosses contact"
        );
        close(report.elapsed_s, DURATION_S, 1e-14);
        close(report.events[0].time_s, INITIAL_TIME_S + 0.004, 1e-8);
        close(report.events[0].contact.impulses[0][0], 1.5, 1e-8);
        close(
            report.events[0].contact.impulses[0][1],
            tangential_impulse,
            1e-8,
        );
        close(report.events[0].contact.impulses[0][2], 0.0, 1e-10);
        let expected = [
            [-0.5, 0.5 - tangential_impulse, 0.0],
            [0.25, -0.25 + 0.5 * tangential_impulse, 0.0],
        ];
        for body in 0..2 {
            let actual = &candidate.motions[body];
            for axis in 0..3 {
                close(actual.velocity_m_s()[axis], expected[body][axis], 1e-8);
                let before = original[body].velocity_m_s()[axis];
                close(
                    actual.state().displacements_m()[axis],
                    before * 0.004 + expected[body][axis] * 0.006,
                    5e-8,
                );
            }
            assert_eq!(
                actual.time_s().to_bits(),
                (INITIAL_TIME_S + DURATION_S).to_bits()
            );
            assert_eq!(actual.steps(), report.accepted_steps);
            assert_eq!(actual.state().points(), original[body].state().points());
            assert_eq!(
                actual.state().displacements_m()[3],
                original[body].state().displacements_m()[3]
            );
            assert!(core::ptr::eq(actual.state().system(), &system));
        }
        for axis in 0..3 {
            close(
                candidate.motions[0].velocity_m_s()[axis]
                    + 2.0 * candidate.motions[1].velocity_m_s()[axis],
                0.0,
                1e-8,
            );
        }
        close(report.initial_gaps_m[0], 0.006, 1e-14);
        close(report.final_gaps_m[0], 0.0045, 5e-8);
        assert!(report.minimum_gap_m >= 0.0);
        close(report.impact_loss_j, loss, 1e-8);
        close(report.kinetic_change_j, -loss, 1e-8);
        for value in [
            report.external_work_j,
            report.prescribed_jump_work_j,
            report.contact_prescribed_work_j,
            report.bulk_stored_change_j,
            report.cohesive_stored_change_j,
            report.fracture_loss_j,
            report.viscous_loss_j,
            report.internal_friction_loss_j,
            report.energy_residual_j,
            report.absolute_smooth_energy_defect_j,
            report.absolute_impact_energy_defect_j,
        ] {
            close(value, 0.0, 1e-8);
        }
    }
    unchanged(&left, &original[0]);
    unchanged(&right, &original[1]);
}

#[test]
fn loaded_rest_refuses_a_finite_interval_without_inventing_sustained_contact_force() {
    let elastic = HeldElastic;
    let masses = [Mass(1.0), Mass(2.0)];
    let law = BilinearCohesive::new(1000.0, 1.0, 0.002, 1000.0, 0.02).unwrap();
    let system = system(&elastic, &law);
    let left = motion(&system, &masses[0], [0.0; 3]);
    let right = motion(&system, &masses[1], [0.0; 3]);
    let original = [left.clone(), right.clone()];
    let rows = pair_trace();
    let forces = [
        [1.0, 0.0, 0.0, HELD_FORCE_N],
        [-1.0, 0.0, 0.0, HELD_FORCE_N],
    ];
    let offsets = [[0.0; 3]];
    let loads: [BodyLoad<'_>; 2] = core::array::from_fn(|body| BodyLoad {
        mean_force_n: &forces[body],
        endpoint_prescribed_jumps_m: &offsets,
    });
    let preconditioners: [&dyn Precond; 2] = [&IdentityPrecond, &IdentityPrecond];
    assert!(matches!(
        advance_candidate_controlled(
            &[&left, &right],
            &[contact(&rows, 0.0, 0.25)],
            &loads,
            &preconditioners,
            DURATION_S,
            controls(),
            &mut || false,
        ),
        Err(TrajectoryError::NoProgress { .. })
    ));
    unchanged(&left, &original[0]);
    unchanged(&right, &original[1]);
}

#[test]
fn body_work_budget_and_late_cancellation_publish_no_partial_timeline() {
    let elastic = HeldElastic;
    let masses = [Mass(1.0), Mass(2.0)];
    let law = BilinearCohesive::new(1000.0, 1.0, 0.002, 1000.0, 0.02).unwrap();
    let system = system(&elastic, &law);
    let left = motion(&system, &masses[0], [1.0, 0.5, 0.0]);
    let right = motion(&system, &masses[1], [-0.5, -0.25, 0.0]);
    let original = [left.clone(), right.clone()];
    let rows = pair_trace();
    let contacts = [contact(&rows, 0.006, 0.25)];
    let force = [0.0, 0.0, 0.0, HELD_FORCE_N];
    let offsets = [[0.0; 3]];
    let loads: [BodyLoad<'_>; 2] = core::array::from_fn(|_| BodyLoad {
        mean_force_n: &force,
        endpoint_prescribed_jumps_m: &offsets,
    });
    let preconditioners: [&dyn Precond; 2] = [&IdentityPrecond, &IdentityPrecond];
    let mut complete_polls = 0;
    let expected = advance_candidate_controlled(
        &[&left, &right],
        &contacts,
        &loads,
        &preconditioners,
        DURATION_S,
        controls(),
        &mut || {
            complete_polls += 1;
            false
        },
    )
    .unwrap();
    let mut limited = controls();
    limited.maximum_body_trials = 1;
    assert!(matches!(
        advance_candidate_controlled(
            &[&left, &right],
            &contacts,
            &loads,
            &preconditioners,
            DURATION_S,
            limited,
            &mut || false,
        ),
        Err(TrajectoryError::Budget(_))
    ));
    let mut polls = 0;
    assert!(matches!(
        advance_candidate_controlled(
            &[&left, &right],
            &contacts,
            &loads,
            &preconditioners,
            DURATION_S,
            controls(),
            &mut || {
                polls += 1;
                polls == complete_polls
            },
        ),
        Err(TrajectoryError::Cancelled)
    ));
    let retry = advance_candidate_controlled(
        &[&left, &right],
        &contacts,
        &loads,
        &preconditioners,
        DURATION_S,
        controls(),
        &mut || false,
    )
    .unwrap();
    for body in 0..2 {
        unchanged(&retry.motions[body], &expected.motions[body]);
    }
    assert_eq!(retry.report.final_gaps_m, expected.report.final_gaps_m);
    assert_eq!(retry.report.impact_loss_j, expected.report.impact_loss_j);
    unchanged(&left, &original[0]);
    unchanged(&right, &original[1]);
}
