//! Original damaged material states survive coupled normal/tangential impacts.

use super::*;
use fs_contact::frictional_multicontact::FrictionalContactConfig;
use fs_couple::cohesive_contact::impact::frictional::{
    FrictionalContact, FrictionalImpactControls, FrictionalResponseControls, NormalVelocityLaw,
    impacts_candidate_controlled,
};

struct Isotropic3(f64);
impl LinearOp for Isotropic3 {
    fn n(&self) -> usize {
        3
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        for i in 0..3 {
            y[i] = self.0 * x[i];
        }
    }
}

fn friction_controls() -> FrictionalImpactControls {
    let original = controls();
    FrictionalImpactControls {
        impulse: original.impulse,
        response: FrictionalResponseControls {
            mass_solve: original.response.mass_solve,
            contact_solve: FrictionalContactConfig::default(),
            maximum_bodies: 4,
            maximum_trace_terms: 24,
            maximum_workspace_entries: 1_000_000,
            velocity_tolerance_m_s: 1e-8,
            energy_tolerance_j: 1e-8,
        },
    }
}

fn system3<'a>(elastic: &'a Isotropic3, law: &'a BilinearCohesive) -> CohesiveSystem<'a> {
    let terms = (0..3)
        .map(|dof| {
            let mut local = [0.0; 3];
            local[dof] = 1.0;
            JumpTerm { dof, local }
        })
        .collect();
    CohesiveSystem::new(
        elastic,
        vec![CohesivePatch::new(41, 0.01, law, terms).unwrap()],
        Budget {
            max_dofs: 3,
            max_patches: 1,
            max_jump_terms: 3,
            max_krylov_basis_entries: 1000,
        },
    )
    .unwrap()
}

fn damaged3<'m, 'a>(
    system: &'m CohesiveSystem<'a>,
    mass: &'m Isotropic3,
    velocity: [f64; 3],
    config: DynamicsConfig,
) -> CohesiveMotion<'m, 'a> {
    let u = 0.002;
    let force = 10.0 * u + 0.01 * (0.004 - u) / 0.003;
    let state = system
        .bonded_state()
        .solve_controlled(&[force, 0.0, 0.0], &[[0.0; 3]], config.solve, &mut || false)
        .unwrap()
        .state;
    assert!(state.points()[0].damage() > 0.0 && state.points()[0].damage() < 1.0);
    state
        .into_dynamics(mass, velocity.to_vec(), 3.0, config)
        .unwrap()
}

fn paired_trace() -> [[TraceTerm; 2]; 3] {
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
fn contact3(
    trace: &[[TraceTerm; 2]; 3],
    friction_coefficient: f64,
    coefficient: f64,
) -> FrictionalContact<'_> {
    FrictionalContact {
        id: 71,
        gap_m: 0.0,
        activation_tolerance_m: 1e-9,
        friction_coefficient,
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient },
        prescribed_velocity_m_s: [0.0; 3],
        trace: [&trace[0], &trace[1], &trace[2]],
    }
}

#[test]
fn sliding_and_sticking_impacts_retain_actual_damaged_owners_and_momentum() {
    let elastic = Isotropic3(10.0);
    let masses = [Isotropic3(1.0), Isotropic3(2.0)];
    let law = BilinearCohesive::new(1000.0, 1.0, 0.002, 1000.0, 0.02).unwrap();
    let system = system3(&elastic, &law);
    let left = damaged3(&system, &masses[0], [1.0, 0.6, 0.0], dynamics());
    let right = damaged3(&system, &masses[1], [-0.5, -0.3, 0.0], dynamics());
    let old_point = left.state().points()[0];
    let old_u = left.state().displacements_m().to_vec();
    let trace = paired_trace();
    let preconditioners: [&dyn Precond; 2] = [&IdentityPrecond, &IdentityPrecond];
    for (mu, tangent_impulse, expected_loss) in [(0.25, 0.25, 0.928125), (0.8, 0.6, 1.02)] {
        let event = impacts_candidate_controlled(
            &[&left, &right],
            &[contact3(&trace, mu, 0.0)],
            &preconditioners,
            friction_controls(),
            &mut || false,
        )
        .unwrap();
        let impulse = event.report.response.contact.impulses[0];
        close(impulse[0], 1.0, 1e-8);
        close(impulse[1], tangent_impulse, 1e-8);
        close(impulse[2], 0.0, 1e-8);
        let a = event.motions[0].velocity_m_s();
        let b = event.motions[1].velocity_m_s();
        close(a[0], 0.0, 1e-8);
        close(b[0], 0.0, 1e-8);
        close(a[1], 0.6 - tangent_impulse, 1e-8);
        close(b[1], -0.3 + 0.5 * tangent_impulse, 1e-8);
        for d in 0..3 {
            close(a[d] + 2.0 * b[d], 0.0, 1e-8);
        }
        close(
            event.report.response.contact.dissipated_energy,
            expected_loss,
            1e-8,
        );
        close(event.report.response.kinetic_change_j, -expected_loss, 1e-8);
        close(event.report.response.dissipation_residual_j, 0.0, 1e-8);
        for motion in &event.motions {
            assert_eq!(motion.state().points(), &[old_point]);
            assert_eq!(motion.state().displacements_m(), old_u);
            assert_eq!(motion.time_s(), 3.0);
            assert_eq!(motion.steps(), 0);
            assert!(core::ptr::eq(motion.state().system(), &system));
        }
        if mu == 0.8 {
            // The next original material step unloads the normal opening;
            // the impact has neither erased damage nor consumed physical time.
            let next = event.motions[0]
                .advance_controlled(&[0.0; 3], &[[0.0; 3]], 1e-4, &mut || false)
                .unwrap();
            assert_eq!(next.motion.steps(), 1);
            assert_eq!(next.motion.state().points()[0].damage(), old_point.damage());
            assert_eq!(next.cohesive_energy_j.dissipated, 0.0);
            close(next.energy_residual_j, 0.0, 1e-9);
        }
    }
    assert_eq!(left.velocity_m_s(), &[1.0, 0.6, 0.0]);
    assert_eq!(right.velocity_m_s(), &[-0.5, -0.3, 0.0]);
    assert_eq!(left.state().points(), &[old_point]);
}

#[test]
fn frictional_event_retains_all_inputs_when_a_later_body_or_publication_refuses() {
    let elastic = Isotropic3(10.0);
    let masses = [Isotropic3(1.0), Isotropic3(2.0)];
    let law = BilinearCohesive::new(1000.0, 1.0, 0.002, 1000.0, 0.02).unwrap();
    let system = system3(&elastic, &law);
    let left = damaged3(&system, &masses[0], [1.0, 1.0, 0.0], dynamics());
    let mut limited = dynamics();
    limited.maximum_velocity_m_s = 0.5;
    let right = damaged3(&system, &masses[1], [0.0; 3], limited);
    let trace = paired_trace();
    let contacts = [contact3(&trace, 0.8, 1.0)];
    let preconditioners: [&dyn Precond; 2] = [&IdentityPrecond, &IdentityPrecond];
    assert!(matches!(
        impacts_candidate_controlled(
            &[&left, &right],
            &contacts,
            &preconditioners,
            friction_controls(),
            &mut || false
        ),
        Err(ImpactError::Body {
            body: 1,
            error: ImpulseError::Velocity { .. }
        })
    ));
    assert_eq!(left.velocity_m_s(), &[1.0, 1.0, 0.0]);
    assert_eq!(right.velocity_m_s(), &[0.0; 3]);
    assert_eq!(right.steps(), 0);

    let allowed = damaged3(&system, &masses[1], [0.0; 3], dynamics());
    let mut total_polls = 0;
    let expected = impacts_candidate_controlled(
        &[&left, &allowed],
        &contacts,
        &preconditioners,
        friction_controls(),
        &mut || {
            total_polls += 1;
            false
        },
    )
    .unwrap();
    let mut polls = 0;
    assert!(matches!(
        impacts_candidate_controlled(
            &[&left, &allowed],
            &contacts,
            &preconditioners,
            friction_controls(),
            &mut || {
                polls += 1;
                polls == total_polls
            }
        ),
        Err(ImpactError::Cancelled)
    ));
    let retry = impacts_candidate_controlled(
        &[&left, &allowed],
        &contacts,
        &preconditioners,
        friction_controls(),
        &mut || false,
    )
    .unwrap();
    for i in 0..2 {
        assert_eq!(
            retry.motions[i].velocity_m_s(),
            expected.motions[i].velocity_m_s()
        );
        assert_eq!(
            retry.motions[i].state().points(),
            expected.motions[i].state().points()
        );
    }
    assert_eq!(left.velocity_m_s(), &[1.0, 1.0, 0.0]);
    assert_eq!(allowed.velocity_m_s(), &[0.0; 3]);
    assert_eq!(allowed.time_s(), 3.0);
}
