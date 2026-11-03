//! G0: instantaneous momentum preserves the original dissipative owners.
use super::*;
use fs_solid::cohesive::dynamics::impulse::ImpulseControls;
use fs_solver::CheckedCgConfig;
use fs_sparse::precond::IdentityPrecond;

fn impulse_controls() -> ImpulseControls {
    ImpulseControls {
        solve: CheckedCgConfig {
            tolerance: 1e-12,
            max_iterations: 16,
            max_corrections: 1,
        },
        max_workspace_entries: 1024,
        absolute_momentum_tolerance_ns: 1e-13,
        relative_momentum_tolerance: 1e-10,
        absolute_work_tolerance_j: 1e-13,
        relative_work_tolerance: 1e-10,
    }
}

#[test]
fn impacted_damaged_motion_retains_sourced_viscosity_and_regularized_friction() {
    let (_, source) = sourced([0.6; 2], [0.3; 2], true);
    let elastic = Scalar(10.0);
    let mass = Scalar(0.1);
    let system = system(&elastic, &source, true);
    let motion = compressed(&system, 0.002)
        .into_dynamics(&mass, vec![0.02], 2.0, mechanical())
        .unwrap()
        .with_damaged_area_friction_controlled(&[Some(1e-4)], &mut || false)
        .unwrap()
        .advance_controlled(&[0.0], &[[-COMPRESSION, 0.0, 0.0]], 0.001, &mut || false)
        .unwrap()
        .motion;
    let old_point = motion.state().points()[0];
    let old_friction = motion.damaged_area_friction().unwrap().to_vec();
    let old_u = motion.state().displacements_m()[0];
    let old_v = motion.velocity_m_s()[0];
    let event = motion
        .impulse_candidate_controlled(&[0.001], &IdentityPrecond, impulse_controls(), &mut || {
            false
        })
        .unwrap();
    let kicked_v = old_v + 0.001 / mass.0;
    close(event.motion.velocity_m_s()[0], kicked_v, 1e-14);
    assert_eq!(event.motion.state().points(), &[old_point]);
    assert_eq!(event.motion.state().displacements_m(), &[old_u]);
    assert_eq!(event.motion.damaged_area_friction().unwrap(), old_friction);
    assert_eq!(event.motion.time_s(), motion.time_s());
    assert_eq!(event.motion.steps(), motion.steps());
    assert!(core::ptr::eq(event.motion.state().system(), &system));
    assert!(core::ptr::eq(
        system.patches()[0].source().unwrap(),
        &source
    ));

    // Positive-duration continuation, not the event, dissipates the retained
    // surviving-bond dashpot and the original damaged-area friction work.
    let dt = 0.001;
    let resistance = AREA * 0.3 * old_point.damage() * KC * COMPRESSION;
    let stiffness = elastic.0 + AREA * 1000.0 * (1.0 - old_point.damage());
    let viscosity = AREA * 4.0 * (1.0 - old_point.damage());
    let du = (2.0 * mass.0 * kicked_v / dt - stiffness * old_u - resistance)
        / (2.0 * mass.0 / (dt * dt) + 0.5 * stiffness + viscosity / dt);
    assert!(du / dt > 1e-4);
    let next = event
        .motion
        .advance_controlled(&[0.0], &[[-COMPRESSION, 0.0, 0.0]], dt, &mut || false)
        .unwrap();
    close(next.motion.state().displacements_m()[0], old_u + du, 1e-12);
    close(
        next.motion.velocity_m_s()[0],
        2.0 * du / dt - kicked_v,
        1e-10,
    );
    close(next.friction_dissipated_j, resistance * du, 1e-13);
    close(next.viscous_dissipated_j, viscosity * du * du / dt, 1e-13);
    assert!(next.friction_dissipated_j > 0.0);
    assert!(next.viscous_dissipated_j > 0.0);
    assert_eq!(next.cohesive_energy_j.dissipated, 0.0);
    assert_eq!(next.motion.state().points()[0].damage(), old_point.damage());
    assert_eq!(
        next.motion.state().points()[0].dissipated_j_m2(),
        old_point.dissipated_j_m2()
    );
    close(
        next.kinetic_before_j,
        event.report.kinetic_energy_after_j,
        1e-15,
    );
    close(next.energy_residual_j, 0.0, 1e-12);
    assert_eq!(next.motion.damaged_area_friction().unwrap(), old_friction);
    assert_eq!(next.motion.steps(), motion.steps() + 1);
    assert_eq!(motion.state().points(), &[old_point]);
    assert_eq!(motion.velocity_m_s(), &[old_v]);
}

#[test]
fn an_impulse_stopping_a_sliding_bond_allows_original_static_capacity_to_recapture() {
    let (_, source) = sourced([0.6; 2], [0.3; 2], true);
    let elastic = Scalar(10.0);
    let mass = Scalar(0.1);
    let system = system(&elastic, &source, true);
    let motion = compressed(&system, 0.002)
        .into_dynamics(&mass, vec![0.02], 2.0, mechanical())
        .unwrap()
        .with_damaged_area_coulomb_controlled(
            &[Some(CoulombActivation::Sliding)],
            controls(),
            &mut || false,
        )
        .unwrap()
        .advance_controlled(&[0.0], &[[-COMPRESSION, 0.0, 0.0]], 0.001, &mut || false)
        .unwrap()
        .motion;
    let old_point = motion.state().points()[0];
    let old_u = motion.state().displacements_m()[0];
    let old_v = motion.velocity_m_s()[0];
    let old_contact = motion.damaged_area_coulomb_states().unwrap().to_vec();
    assert_eq!(old_contact[0].occupancy, CoulombRegime::Sliding);
    assert!(old_contact[0].last_traction_pa[0] > 0.0);
    let event = motion
        .impulse_candidate_controlled(
            &[-mass.0 * old_v],
            &IdentityPrecond,
            impulse_controls(),
            &mut || false,
        )
        .unwrap();
    close(event.motion.velocity_m_s()[0], 0.0, 1e-14);
    close(
        event.report.applied_impulse_work_j,
        -0.5 * mass.0 * old_v * old_v,
        1e-15,
    );
    assert_eq!(event.motion.state().points(), &[old_point]);
    assert_eq!(event.motion.state().displacements_m(), &[old_u]);
    assert_eq!(
        event.motion.damaged_area_coulomb_states().unwrap(),
        old_contact
    );
    assert_eq!(event.motion.time_s(), motion.time_s());
    assert_eq!(event.motion.steps(), motion.steps());
    assert!(core::ptr::eq(event.motion.state().system(), &system));

    // A finite-duration owner re-evaluates capture using the actual post-event
    // velocity. Reusing the old Sliding occupancy would select only mu_k and
    // make this same stopped contact accelerate under a supportable load.
    let kinetic_capacity = AREA * 0.3 * old_point.damage() * KC * COMPRESSION;
    let static_capacity = AREA * 0.6 * old_point.damage() * KC * COMPRESSION;
    let reaction = 0.5 * (kinetic_capacity + static_capacity);
    let stiffness = elastic.0 + AREA * 1000.0 * (1.0 - old_point.damage());
    let force = stiffness * old_u + reaction;
    let held = event
        .motion
        .advance_controlled(&[force], &[[-COMPRESSION, 0.0, 0.0]], 0.001, &mut || false)
        .unwrap();
    close(held.motion.state().displacements_m()[0], old_u, 1e-12);
    close(held.motion.velocity_m_s()[0], 0.0, 1e-10);
    let contact = held.motion.damaged_area_coulomb_states().unwrap()[0];
    assert_eq!(contact.occupancy, CoulombRegime::Sticking);
    close(AREA * contact.last_traction_pa[0], reaction, 1e-10);
    assert_eq!(held.friction_dissipated_j, 0.0);
    close(held.viscous_dissipated_j, 0.0, 1e-15);
    assert_eq!(held.cohesive_energy_j.dissipated, 0.0);
    assert_eq!(held.motion.state().points()[0].damage(), old_point.damage());
    assert_eq!(
        held.motion.state().points()[0].dissipated_j_m2(),
        old_point.dissipated_j_m2()
    );
    close(held.energy_residual_j, 0.0, 1e-12);
    assert_eq!(held.motion.steps(), motion.steps() + 1);
    assert_eq!(motion.damaged_area_coulomb_states().unwrap(), old_contact);
    assert_eq!(motion.velocity_m_s(), &[old_v]);
}
