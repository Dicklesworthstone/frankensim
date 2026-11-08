//! Physical momentum jumps, unchanged material history and original source work.
use fs_phs::RelaxationBranch;
use fs_phs::linear_exact::{LinearSplitBudget, PreparedExactZoh};
use fs_phs::mechanical::{LinearMechanicalSystem, MechanicalImpulseControls};

fn controls() -> MechanicalImpulseControls {
    MechanicalImpulseControls {
        absolute_momentum_tolerance: 1e-14,
        relative_momentum_tolerance: 1e-13,
        work_tolerance: 1e-12,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}
fn bits(state: &[f64]) -> Vec<u64> {
    state.iter().map(|v| v.to_bits()).collect()
}
fn coupled() -> LinearMechanicalSystem {
    LinearMechanicalSystem::new(
        &[2.0, 0.4, 0.4, 1.5],
        &[0.8, 1.55, -1.05, 0.6],
        &[7.0, -1.0, -1.0, 4.0],
        &[1.0, 0.2, -0.3, 0.8],
        2,
        2,
        vec![
            RelaxationBranch {
                projection: vec![1.0, -0.6],
                stiffness: 3.0,
                relaxation_time_s: 0.08,
            },
            RelaxationBranch {
                projection: vec![0.3, 1.0],
                stiffness: 2.0,
                relaxation_time_s: 0.6,
            },
        ],
    )
    .unwrap()
}

#[test]
fn consistent_mass_impulse_and_braking_retain_gyro_and_all_maxwell_history() {
    let model = coupled();
    let initial = model
        .state(&[0.2, -0.15], &[0.1, 0.25], &[0.05, -0.1])
        .unwrap();
    let snapshot = bits(&initial);
    let impulse = [0.3, -0.7];
    let jump = model.apply_impulse(&initial, &impulse, controls()).unwrap();
    let delta = [
        (1.5 * 0.3 + 0.4 * 0.7) / 2.84,
        (-0.4 * 0.3 - 2.0 * 0.7) / 2.84,
    ];
    let mut velocity = [0.0; 2];
    model
        .velocities_into(&jump.state_after, &mut velocity)
        .unwrap();
    let mut expected_work = 0.0;
    for i in 0..2 {
        close(jump.velocity_change[i], delta[i], 3e-16);
        close(velocity[i], [0.1, 0.25][i] + delta[i], 3e-16);
        close(jump.momentum_change[i], impulse[i], 2e-16);
        close(jump.momentum_residual[i], 0.0, 2e-16);
        expected_work += impulse[i] * ([0.1, 0.25][i] + 0.5 * delta[i]);
    }
    close(jump.supplied_work, expected_work, 3e-16);
    close(jump.kinetic_energy_change, expected_work, 3e-16);
    close(jump.work_residual, 0.0, 3e-16);
    assert_eq!(jump.dissipated, 0.0);
    assert_eq!(bits(&jump.state_after[..2]), snapshot[..2]);
    assert_eq!(bits(&jump.state_after[4..]), snapshot[4..]);
    assert_eq!(
        model.branch_stresses(&initial).unwrap(),
        model.branch_stresses(&jump.state_after).unwrap()
    );
    assert_eq!(bits(&initial), snapshot);

    let stopped = model
        .apply_impulse(&initial, &[-initial[2], -initial[3]], controls())
        .unwrap();
    assert_eq!(model.velocities(&stopped.state_after).unwrap(), [0.0, 0.0]);
    let kinetic = 0.5 * (2.0 * 0.1_f64.powi(2) + 0.8 * 0.1 * 0.25 + 1.5 * 0.25_f64.powi(2));
    close(stopped.supplied_work, -kinetic, 1e-16);
    close(stopped.kinetic_energy_change, -kinetic, 1e-16);
    assert_eq!(stopped.dissipated, 0.0);
    assert_eq!(bits(&stopped.state_after[4..]), snapshot[4..]);
}

#[test]
fn port_impulses_are_simultaneous_with_signed_source_work_and_cancellation() {
    let model =
        LinearMechanicalSystem::new(&[2.0], &[0.0], &[0.0], &[1.0, 1.0], 1, 2, vec![]).unwrap();
    let initial = model.state(&[0.2], &[0.3], &[]).unwrap();
    let jump = model
        .apply_port_impulse(&initial, &[0.5, -0.2], controls())
        .unwrap();
    close(model.velocities(&jump.state_after).unwrap()[0], 0.45, 2e-16);
    close(jump.port_work[0], 0.5 * 0.375, 2e-16);
    close(jump.port_work[1], -0.2 * 0.375, 2e-16);
    close(jump.supplied_work, 0.1125, 2e-16);
    close(jump.kinetic_energy_change, 0.1125, 2e-16);
    let cancel = model
        .apply_port_impulse(&initial, &[2.0, -2.0], controls())
        .unwrap();
    assert_eq!(bits(&cancel.state_after), bits(&initial));
    close(cancel.port_work[0], 0.6, 2e-16);
    close(cancel.port_work[1], -0.6, 2e-16);
    assert_eq!(cancel.supplied_work, 0.0);
    assert_eq!(cancel.kinetic_energy_change, 0.0);
    assert_eq!(cancel.dissipated, 0.0);

    let large = 2.0_f64.powi(54);
    let source =
        LinearMechanicalSystem::new(&[1.0], &[0.0], &[0.0], &[large, 1.0, -large], 1, 3, vec![])
            .unwrap();
    let jump = source
        .apply_port_impulse(&[0.0, 0.0], &[1.0; 3], controls())
        .unwrap();
    assert_eq!(jump.state_after, [0.0, 1.0]);
    assert_eq!(jump.port_work, [0.5 * large, 0.5, -0.5 * large]);
    assert_eq!(jump.supplied_work, 0.5);
    assert_eq!(jump.kinetic_energy_change, 0.5);
    assert_eq!(jump.work_residual, 0.0);
}

#[test]
fn weak_work_survives_opposing_momenta_and_large_fixed_elastic_memory_energy() {
    let identity = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    let model =
        LinearMechanicalSystem::new(&identity, &[0.0; 9], &[0.0; 9], &[], 3, 0, vec![]).unwrap();
    let large = 2.0_f64.powi(52);
    let initial = model.state(&[0.0; 3], &[large, 0.5, -large], &[]).unwrap();
    let jump = model
        .apply_impulse(&initial, &[1.0; 3], controls())
        .unwrap();
    assert_eq!(jump.momentum_change, [1.0; 3]);
    assert_eq!(jump.supplied_work, 2.0);
    assert_eq!(jump.kinetic_energy_change, 2.0);
    assert_eq!(jump.work_residual, 0.0);

    let loaded = LinearMechanicalSystem::new(
        &[2.0],
        &[3.0],
        &[8.0],
        &[],
        1,
        0,
        vec![RelaxationBranch {
            projection: vec![1.0],
            stiffness: 3.0,
            relaxation_time_s: 0.4,
        }],
    )
    .unwrap();
    let initial = loaded.state(&[1e12], &[0.2], &[-1e12]).unwrap();
    let jump = loaded.apply_impulse(&initial, &[0.6], controls()).unwrap();
    // Subtraction of total endpoint energies loses the entire kinetic event.
    assert_eq!(
        loaded.system().hamiltonian(&jump.state_after),
        loaded.system().hamiltonian(&initial)
    );
    close(jump.supplied_work, 0.21, 2e-16);
    close(jump.kinetic_energy_change, 0.21, 2e-16);
    assert_eq!(bits(&jump.state_after[2..]), bits(&initial[2..]));
    assert_eq!(jump.dissipated, 0.0);
}

#[test]
fn impulse_continues_existing_prepared_oscillator_without_advancing_time() {
    let model = LinearMechanicalSystem::new(&[2.0], &[0.0], &[8.0], &[], 1, 0, vec![]).unwrap();
    let duration = 0.7;
    let flow = PreparedExactZoh::new(
        model.system(),
        duration,
        LinearSplitBudget {
            max_modes: 4,
            max_omega_dt: 100.0,
        },
    )
    .unwrap();
    let initial = model.state(&[0.4], &[-0.3], &[]).unwrap();
    let jump = model.apply_impulse(&initial, &[0.6], controls()).unwrap();
    assert_eq!(jump.state_after[0], 0.4);
    assert_eq!(jump.state_after[1], 0.0);
    let next = flow.step(&jump.state_after, &[]).unwrap();
    close(next.x[0], 0.4 * (2.0 * duration).cos(), 2e-13);
    close(
        model.velocities(&next.x).unwrap()[0],
        -0.8 * (2.0 * duration).sin(),
        2e-13,
    );
    assert_eq!(jump.dissipated, 0.0);
    assert_eq!(next.dissipated, 0.0);
}

#[test]
fn physical_initialization_keeps_cancelling_inertia_and_weak_negative_storage_refuses() {
    let mass = [4.0, 1.0, 1.0, 1.0, 2.0, 0.0, 1.0, 0.0, 2.0];
    let model =
        LinearMechanicalSystem::new(&mass, &[0.0; 9], &[0.0; 9], &[], 3, 0, vec![]).unwrap();
    let large = 2.0_f64.powi(54);
    let state = model
        .state(&[0.0; 3], &[large / 4.0, 1.0, -large], &[])
        .unwrap();
    assert_eq!(state[3], 1.0);
    let identity = [1.0, 0.0, 0.0, 1.0];
    let weak_negative = [1.0, 0.0, 0.0, -1e-300];
    assert!(
        LinearMechanicalSystem::new(&identity, &[0.0; 4], &weak_negative, &[], 2, 0, vec![])
            .is_err()
    );
    assert!(
        LinearMechanicalSystem::new(&identity, &weak_negative, &[0.0; 4], &[], 2, 0, vec![])
            .is_err()
    );
}

#[test]
fn zero_event_preserves_bits_and_refusals_leave_original_state_unchanged() {
    let model = LinearMechanicalSystem::new(&[1.0], &[0.0], &[0.0], &[1.0], 1, 1, vec![]).unwrap();
    let signed = [-0.0, -0.0];
    assert_eq!(
        bits(
            &model
                .apply_impulse(&signed, &[0.0], controls())
                .unwrap()
                .state_after
        ),
        bits(&signed)
    );
    let initial = [0.2, 2.0_f64.powi(54)];
    let snapshot = bits(&initial);
    for impulse in [&[][..], &[1.0], &[f64::NAN], &[f64::MAX]] {
        assert!(model.apply_impulse(&initial, impulse, controls()).is_err());
        assert_eq!(bits(&initial), snapshot);
    }
    assert!(
        model
            .apply_impulse(&initial[..1], &[0.0], controls())
            .is_err()
    );
    assert!(model.apply_port_impulse(&initial, &[], controls()).is_err());
    assert!(
        model
            .apply_port_impulse(&initial, &[f64::INFINITY], controls())
            .is_err()
    );
    for invalid in [
        MechanicalImpulseControls {
            work_tolerance: -1.0,
            ..controls()
        },
        MechanicalImpulseControls {
            relative_momentum_tolerance: 1.0,
            ..controls()
        },
        MechanicalImpulseControls {
            absolute_momentum_tolerance: f64::NAN,
            ..controls()
        },
    ] {
        assert!(model.apply_impulse(&initial, &[0.0], invalid).is_err());
    }
    assert!(model.velocities_into(&initial, &mut []).is_err());
    assert!(model.velocities_into(&[0.0, f64::NAN], &mut [0.0]).is_err());
    assert_eq!(bits(&initial), snapshot);
}
