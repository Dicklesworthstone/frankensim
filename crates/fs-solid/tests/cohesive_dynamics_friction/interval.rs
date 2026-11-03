//! Complete native intervals across unaligned stops; no velocity projection.
use super::*;
use fs_solid::cohesive::dynamics::interval::IntervalConfig;

fn interval_controls() -> IntervalConfig {
    IntervalConfig {
        minimum_step_s: 1e-12,
        reversal_time_tolerance_s: 1e-11,
        maximum_trials: 160,
        maximum_substeps: 30,
        maximum_total_residual_evaluations: 200_000,
        maximum_total_newton_iterations: 6000,
    }
}

#[test]
fn interval_crosses_an_unaligned_stop_then_holds_and_breaks_away() {
    let (_, source) = sourced([0.6; 2], [0.3; 2], false);
    let elastic = Scalar(0.0);
    let mass = Scalar(0.1);
    let system = system(&elastic, &source, true);
    let v0 = 0.00037;
    let motion = compressed(&system, 0.005)
        .into_dynamics(&mass, vec![v0], 1.0, mechanical())
        .unwrap()
        .with_damaged_area_coulomb_controlled(
            &[Some(CoulombActivation::Sliding)],
            controls(),
            &mut || false,
        )
        .unwrap();
    let offsets = [[-COMPRESSION, 0.0, 0.0]];
    let stop_time = mass.0 * v0 / 0.006;
    let stop_distance = mass.0 * v0 * v0 / (2.0 * 0.006);
    assert!(stop_time < 0.02 && stop_time != 0.01);
    let stopped = motion
        .advance_interval_controlled(&[0.0], &offsets, 0.02, interval_controls(), &mut || false)
        .unwrap();
    assert_eq!(stopped.motion.time_s(), 1.02);
    assert!(stopped.accepted_substeps >= 2);
    assert!(stopped.attempted_substeps > stopped.accepted_substeps);
    assert!(stopped.motion.velocity_m_s()[0].abs() <= controls().velocity_tolerance_m_s);
    close(
        stopped.motion.state().displacements_m()[0],
        stop_distance,
        1e-12,
    );
    close(stopped.friction_dissipated_j, 0.5 * mass.0 * v0 * v0, 1e-15);
    close(
        stopped.friction_dissipated_j,
        stopped.patch_dissipation[0].friction_j,
        1e-20,
    );
    assert_eq!(stopped.patch_dissipation[0].fracture_j, 0.0);
    assert_eq!(stopped.patch_dissipation[0].viscous_j, 0.0);
    close(stopped.energy_residual_j, 0.0, 1e-13);
    assert_eq!(
        stopped.motion.damaged_area_coulomb_states().unwrap()[0].occupancy,
        CoulombRegime::Sticking
    );
    assert_eq!(stopped.motion.state().points()[0].damage(), 1.0);
    let held = stopped
        .motion
        .advance_interval_controlled(&[0.009], &offsets, 0.012, interval_controls(), &mut || {
            false
        })
        .unwrap();
    close(
        held.motion.state().displacements_m()[0],
        stop_distance,
        1e-12,
    );
    assert_eq!(held.friction_dissipated_j, 0.0);
    let incoming = held.motion.velocity_m_s()[0];
    let sliding = held
        .motion
        .advance_interval_controlled(&[0.02], &offsets, 0.01, interval_controls(), &mut || false)
        .unwrap();
    let du = incoming * 0.01 + (0.02 - 0.006) * 0.01_f64.powi(2) / (2.0 * mass.0);
    close(
        sliding.motion.state().displacements_m()[0] - held.motion.state().displacements_m()[0],
        du,
        1e-12,
    );
    close(sliding.friction_dissipated_j, 0.006 * du, 1e-13);
    assert_eq!(
        sliding.motion.damaged_area_coulomb_states().unwrap()[0].occupancy,
        CoulombRegime::Sliding
    );
    assert_eq!(motion.velocity_m_s(), &[v0]);
    assert_eq!(motion.time_s(), 1.0);
    assert_eq!(motion.steps(), 0);
}

#[test]
fn interval_sums_original_fracture_viscosity_and_friction_on_prescribed_paths() {
    let (_, source) = sourced([0.6; 2], [0.3; 2], true);
    let elastic = Scalar(10.0);
    let mass = Scalar(0.1);
    let system = system(&elastic, &source, false);
    let start = [-COMPRESSION, 0.0005, 0.0];
    let end = [-COMPRESSION, 0.003, 0.0];
    let state = system
        .bonded_state()
        .solve_over_time_controlled(&[0.0], &[start], 1.0, mechanical().solve, &mut || false)
        .unwrap()
        .state;
    let motion = state
        .into_dynamics(&mass, vec![0.0], 0.0, mechanical())
        .unwrap()
        .with_damaged_area_coulomb_controlled(
            &[Some(CoulombActivation::AtRest)],
            controls(),
            &mut || false,
        )
        .unwrap();
    let integrated = motion
        .advance_interval_controlled(&[0.0], &[end], 0.04, interval_controls(), &mut || false)
        .unwrap();
    assert_eq!(integrated.accepted_substeps, 2);
    let mid = [-COMPRESSION, (start[1] + end[1]) * 0.5, 0.0];
    let a = motion
        .advance_controlled(&[0.0], &[mid], 0.02, &mut || false)
        .unwrap();
    let b = a
        .motion
        .advance_controlled(&[0.0], &[end], 0.02, &mut || false)
        .unwrap();
    assert_eq!(
        integrated.motion.state().points(),
        b.motion.state().points()
    );
    for (got, expected) in [
        (
            integrated.friction_dissipated_j,
            a.friction_dissipated_j + b.friction_dissipated_j,
        ),
        (
            integrated.viscous_dissipated_j,
            a.viscous_dissipated_j + b.viscous_dissipated_j,
        ),
        (
            integrated.cohesive_energy_j.dissipated,
            a.cohesive_energy_j.dissipated + b.cohesive_energy_j.dissipated,
        ),
        (
            integrated.prescribed_jump_work_j,
            a.prescribed_jump_work_j + b.prescribed_jump_work_j,
        ),
    ] {
        close(got, expected, 1e-15);
        assert!(got > 0.0);
    }
    close(
        integrated.patch_dissipation[0].fracture_j,
        integrated.cohesive_energy_j.dissipated,
        1e-15,
    );
    close(
        integrated.patch_dissipation[0].viscous_j,
        integrated.viscous_dissipated_j,
        1e-15,
    );
    close(
        integrated.patch_dissipation[0].friction_j,
        integrated.friction_dissipated_j,
        1e-15,
    );
    close(integrated.energy_residual_j, 0.0, 1e-12);
}

#[test]
fn interval_refusal_or_late_cancel_preserves_all_original_motion_and_history() {
    let (_, source) = sourced([0.6; 2], [0.3; 2], false);
    let elastic = Scalar(0.0);
    let mass = Scalar(0.1);
    let system = system(&elastic, &source, true);
    let motion = compressed(&system, 0.005)
        .into_dynamics(&mass, vec![0.00037], 0.0, mechanical())
        .unwrap()
        .with_damaged_area_coulomb_controlled(
            &[Some(CoulombActivation::Sliding)],
            controls(),
            &mut || false,
        )
        .unwrap();
    let offsets = [[-COMPRESSION, 0.0, 0.0]];
    let point = motion.state().points()[0];
    let state = motion.damaged_area_coulomb_states().unwrap()[0];
    for limited in [
        IntervalConfig {
            maximum_trials: 1,
            ..interval_controls()
        },
        IntervalConfig {
            maximum_substeps: 1,
            ..interval_controls()
        },
        IntervalConfig {
            maximum_total_residual_evaluations: 1,
            ..interval_controls()
        },
        IntervalConfig {
            maximum_total_newton_iterations: 1,
            ..interval_controls()
        },
    ] {
        assert!(
            motion
                .advance_interval_controlled(&[0.0], &offsets, 0.02, limited, &mut || false)
                .is_err()
        );
    }
    let mut calls = 0;
    let complete = motion
        .advance_interval_controlled(&[0.0], &offsets, 0.02, interval_controls(), &mut || {
            calls += 1;
            false
        })
        .unwrap();
    assert!(
        complete.charged_residual_evaluations
            <= interval_controls().maximum_total_residual_evaluations
    );
    assert!(
        complete.charged_newton_iterations <= interval_controls().maximum_total_newton_iterations
    );
    let mut count = 0;
    assert!(matches!(
        motion.advance_interval_controlled(
            &[0.0],
            &offsets,
            0.02,
            interval_controls(),
            &mut || {
                count += 1;
                count == calls
            },
        ),
        Err(DynamicsError::Model(Error::Cancelled))
    ));
    assert_eq!(motion.state().points()[0], point);
    assert_eq!(motion.damaged_area_coulomb_states().unwrap()[0], state);
    assert_eq!(motion.velocity_m_s(), &[0.00037]);
    assert_eq!(motion.state().displacements_m(), &[0.0]);
    assert_eq!(motion.time_s(), 0.0);
    assert_eq!(motion.steps(), 0);
}
