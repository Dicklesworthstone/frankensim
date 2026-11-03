//! G0: the native contact law updates the original damaged motion owners.

use fs_contact::multicontact::MultiImpactConfig;
use fs_couple::cohesive_contact::impact::{
    ImpactControls, ImpactError, NormalContact, ResponseControls, TraceTerm,
    normal_impacts_candidate_controlled,
};
use fs_material::cohesive::BilinearCohesive;
use fs_solid::cohesive::dynamics::impulse::{ImpulseControls, ImpulseError};
use fs_solid::cohesive::dynamics::{CohesiveMotion, DynamicsConfig};
use fs_solid::cohesive::{Budget, CohesivePatch, CohesiveSystem, Config, JumpTerm};
use fs_solver::{CheckedCgConfig, LinearOp, NewtonKrylovConfig};
use fs_sparse::precond::{IdentityPrecond, Precond};

struct Scalar(f64);
impl LinearOp for Scalar {
    fn n(&self) -> usize {
        1
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        y[0] = self.0 * x[0];
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
            maximum_displacement_m: 0.02,
        },
        difference_step_m: 1e-8,
        maximum_step_s: 0.02,
        maximum_velocity_m_s: 2.0,
        maximum_steps: 1000,
        absolute_energy_tolerance_j: 1e-10,
        relative_energy_tolerance: 1e-8,
    }
}

fn controls() -> ImpactControls {
    let solve = CheckedCgConfig {
        tolerance: 1e-13,
        max_iterations: 32,
        max_corrections: 2,
    };
    ImpactControls {
        response: ResponseControls {
            mass_solve: solve,
            contact_solve: MultiImpactConfig::default(),
            maximum_bodies: 4,
            maximum_trace_terms: 24,
            maximum_workspace_entries: 100_000,
            velocity_tolerance_m_s: 1e-9,
            energy_tolerance_j: 1e-9,
        },
        impulse: ImpulseControls {
            solve,
            max_workspace_entries: 2048,
            absolute_momentum_tolerance_ns: 1e-12,
            relative_momentum_tolerance: 1e-12,
            absolute_work_tolerance_j: 1e-12,
            relative_work_tolerance: 1e-12,
        },
    }
}

fn system<'a>(elastic: &'a Scalar, law: &'a BilinearCohesive) -> CohesiveSystem<'a> {
    CohesiveSystem::new(
        elastic,
        vec![
            CohesivePatch::new(
                41,
                0.01,
                law,
                vec![JumpTerm {
                    dof: 0,
                    local: [1.0, 0.0, 0.0],
                }],
            )
            .unwrap(),
        ],
        Budget {
            max_dofs: 1,
            max_patches: 1,
            max_jump_terms: 1,
            max_krylov_basis_entries: 1000,
        },
    )
    .unwrap()
}

fn damaged<'m, 'a>(
    system: &'m CohesiveSystem<'a>,
    mass: &'m Scalar,
    velocity: f64,
    config: DynamicsConfig,
) -> CohesiveMotion<'m, 'a> {
    let u = 0.002;
    let force = 10.0 * u + 0.01 * (0.004 - u) / 0.003;
    let state = system
        .bonded_state()
        .solve_controlled(&[force], &[[0.0; 3]], config.solve, &mut || false)
        .unwrap()
        .state;
    assert!(state.points()[0].damage() > 0.0 && state.points()[0].damage() < 1.0);
    state
        .into_dynamics(mass, vec![velocity], 3.0, config)
        .unwrap()
}

fn contact(trace: &[TraceTerm], restitution: f64) -> NormalContact<'_> {
    NormalContact {
        id: 9,
        gap_m: 0.0,
        activation_tolerance_m: 1e-9,
        restitution,
        prescribed_normal_velocity_m_s: 0.0,
        trace,
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn unequal_damaged_bodies_restitute_with_original_masses_and_continue_natively() {
    let elastic = Scalar(10.0);
    let masses = [Scalar(1.0), Scalar(2.0)];
    let law = BilinearCohesive::new(1000.0, 1.0, 0.002, 1000.0, 0.02).unwrap();
    let system = system(&elastic, &law);
    let left = damaged(&system, &masses[0], 1.0, dynamics());
    let right = damaged(&system, &masses[1], -0.5, dynamics());
    let old_point = left.state().points()[0];
    let old_u = left.state().displacements_m()[0];
    let trace = [
        TraceTerm {
            body: 0,
            dof: 0,
            coefficient: -1.0,
        },
        TraceTerm {
            body: 1,
            dof: 0,
            coefficient: 1.0,
        },
    ];
    let preconditioners: [&dyn Precond; 2] = [&IdentityPrecond, &IdentityPrecond];
    for restitution in [0.0, 0.5, 1.0] {
        let event = normal_impacts_candidate_controlled(
            &[&left, &right],
            &[contact(&trace, restitution)],
            &preconditioners,
            controls(),
            &mut || false,
        )
        .unwrap();
        // Independent two-mass oracle: W=1/1+1/2, approach speed=1.5.
        let impulse = 1.0 + restitution;
        close(
            event.report.response.contact.normal_impulses[0],
            impulse,
            1e-10,
        );
        close(event.motions[0].velocity_m_s()[0], 1.0 - impulse, 1e-10);
        close(
            event.motions[1].velocity_m_s()[0],
            -0.5 + impulse / 2.0,
            1e-10,
        );
        close(
            event.motions[0].velocity_m_s()[0] + 2.0 * event.motions[1].velocity_m_s()[0],
            0.0,
            1e-10,
        );
        close(
            event.report.response.contact.dissipated_energy,
            0.75 * (1.0 - restitution * restitution),
            1e-10,
        );
        close(
            event.report.response.kinetic_change_j,
            -0.75 * (1.0 - restitution * restitution),
            1e-10,
        );
        for motion in &event.motions {
            assert_eq!(motion.state().points(), &[old_point]);
            assert_eq!(motion.state().displacements_m(), &[old_u]);
            assert_eq!(motion.time_s(), 3.0);
            assert_eq!(motion.steps(), 0);
            assert!(core::ptr::eq(motion.state().system(), &system));
        }
        if restitution == 0.0 {
            let next = event.motions[0]
                .advance_controlled(&[0.0], &[[0.0; 3]], 1e-4, &mut || false)
                .unwrap();
            assert_eq!(next.motion.steps(), 1);
            assert_eq!(next.motion.state().points()[0].damage(), old_point.damage());
            assert_eq!(next.cohesive_energy_j.dissipated, 0.0);
            close(next.energy_residual_j, 0.0, 1e-10);
        }
    }
    assert_eq!(left.velocity_m_s(), &[1.0]);
    assert_eq!(right.velocity_m_s(), &[-0.5]);
    assert_eq!(left.state().points(), &[old_point]);
}

#[test]
fn second_body_refusal_and_late_cancellation_leave_the_entire_event_unpublished() {
    let elastic = Scalar(10.0);
    let masses = [Scalar(1.0), Scalar(2.0)];
    let law = BilinearCohesive::new(1000.0, 1.0, 0.002, 1000.0, 0.02).unwrap();
    let system = system(&elastic, &law);
    let left = damaged(&system, &masses[0], 1.0, dynamics());
    let mut limited = dynamics();
    limited.maximum_velocity_m_s = 0.5;
    let right = damaged(&system, &masses[1], 0.0, limited);
    let trace = [
        TraceTerm {
            body: 0,
            dof: 0,
            coefficient: -1.0,
        },
        TraceTerm {
            body: 1,
            dof: 0,
            coefficient: 1.0,
        },
    ];
    let contacts = [contact(&trace, 1.0)];
    let preconditioners: [&dyn Precond; 2] = [&IdentityPrecond, &IdentityPrecond];
    assert!(matches!(
        normal_impacts_candidate_controlled(
            &[&left, &right],
            &contacts,
            &preconditioners,
            controls(),
            &mut || false,
        ),
        Err(ImpactError::Body {
            body: 1,
            error: ImpulseError::Velocity { .. }
        })
    ));
    assert_eq!(left.velocity_m_s(), &[1.0]);
    assert_eq!(right.velocity_m_s(), &[0.0]);
    assert_eq!(right.time_s(), 3.0);
    assert_eq!(right.steps(), 0);

    let admitted_right = damaged(&system, &masses[1], 0.0, dynamics());
    let mut complete_polls = 0;
    let expected = normal_impacts_candidate_controlled(
        &[&left, &admitted_right],
        &contacts,
        &preconditioners,
        controls(),
        &mut || {
            complete_polls += 1;
            false
        },
    )
    .unwrap();
    let mut polls = 0;
    assert!(matches!(
        normal_impacts_candidate_controlled(
            &[&left, &admitted_right],
            &contacts,
            &preconditioners,
            controls(),
            &mut || {
                polls += 1;
                polls == complete_polls
            },
        ),
        Err(ImpactError::Cancelled)
    ));
    let retry = normal_impacts_candidate_controlled(
        &[&left, &admitted_right],
        &contacts,
        &preconditioners,
        controls(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(
        retry.motions[0].velocity_m_s(),
        expected.motions[0].velocity_m_s()
    );
    assert_eq!(
        retry.motions[1].velocity_m_s(),
        expected.motions[1].velocity_m_s()
    );
    assert_eq!(left.velocity_m_s(), &[1.0]);
    assert_eq!(admitted_right.velocity_m_s(), &[0.0]);
    assert_eq!(left.state().points(), admitted_right.state().points());
}
