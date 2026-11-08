//! G1/G3 canonical parameter changes against original physical dynamics.
use fs_phs::RelaxationBranch;
use fs_phs::linear_exact::{LinearSplitBudget, PreparedExactZoh};
use fs_phs::mechanical::{
    LinearMechanicalSystem, MechanicalCoordinateMap, MechanicalTransferControls,
    MechanicalTransferError,
};

fn controls() -> MechanicalTransferControls {
    MechanicalTransferControls {
        maximum_coordinates: 8,
        maximum_payload_bytes: 1 << 20,
        maximum_condition_number: 1e8,
        maximum_relative_mapping_residual: 2e-13,
        maximum_relative_strain_residual: 2e-13,
        maximum_absolute_parameter_work_j: 10.0,
        energy_tolerance_j: 2e-12,
    }
}
fn budget() -> LinearSplitBudget {
    LinearSplitBudget {
        max_modes: 32,
        max_omega_dt: 100.0,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}; tolerance={tolerance:e}"
    );
}
fn arm(projection: &[f64], stiffness: f64, time: f64) -> RelaxationBranch {
    RelaxationBranch {
        projection: projection.to_vec(),
        stiffness,
        relaxation_time_s: time,
    }
}
fn scalar(mass: f64, stiffness: f64, branch: f64, time: f64) -> LinearMechanicalSystem {
    LinearMechanicalSystem::new(
        &[mass],
        &[0.3],
        &[stiffness],
        &[1.0],
        1,
        1,
        vec![arm(&[1.0], branch, time)],
    )
    .unwrap()
}

#[test]
fn changed_inertia_stiffness_and_relaxation_preserve_canonical_motion_and_viscous_strain() {
    let previous = scalar(2.0, 3.0, 4.0, 0.5);
    let target = scalar(8.0, 7.0, 16.0, 0.1);
    let state = previous.state(&[0.3], &[0.4], &[0.1]).unwrap();
    let plan = previous
        .prepare_state_transfer(
            &target,
            MechanicalCoordinateMap::Unchanged,
            &[0],
            controls(),
            || false,
        )
        .unwrap();
    assert!(core::ptr::eq(plan.previous(), &previous));
    assert!(core::ptr::eq(plan.target(), &target));
    let transferred = plan.transfer(&state, || false).unwrap();
    assert_eq!(transferred.state_after[..2], state[..2]);
    close(transferred.state_after[2], 0.4, 1e-16);
    close(
        target.velocities(&transferred.state_after).unwrap()[0],
        0.1,
        1e-16,
    );
    close(
        target.branch_stresses(&transferred.state_after).unwrap()[0],
        3.2,
        1e-14,
    );
    let report = transferred.report;
    close(report.kinetic_parameter_work_j, -0.12, 1e-15);
    close(report.elastic_parameter_work_j, 0.18, 1e-15);
    close(report.memory_parameter_work_j, 0.24, 1e-15);
    close(report.parameter_work_j, 0.30, 1e-15);
    close(report.storage_change_j, 0.30, 1e-15);
    close(report.energy_residual_j, 0.0, 1e-15);
    close(report.previous_stored_energy_j, 0.375, 1e-15);
    close(report.target_stored_energy_j, 0.675, 1e-15);

    let returned = target
        .prepare_state_transfer(
            &previous,
            MechanicalCoordinateMap::Unchanged,
            &[0],
            controls(),
            || false,
        )
        .unwrap()
        .transfer(&transferred.state_after, || false)
        .unwrap();
    assert_eq!(returned.state_after, state);
    close(returned.report.parameter_work_j, -0.30, 1e-15);

    // Independent physical ODE [q,v,viscous strain,heat,force work]. All new
    // coefficients act immediately; neither velocity nor the old z is copied.
    let derivative = |x: [f64; 5]| {
        let stress = 16.0 * (x[0] - x[2]);
        [
            x[1],
            (1.2 - 0.3 * x[1] - 7.0 * x[0] - stress) / 8.0,
            (x[0] - x[2]) / 0.1,
            0.3 * x[1] * x[1] + stress * stress / (16.0 * 0.1),
            1.2 * x[1],
        ]
    };
    let mut reference = [0.3, 0.1, 0.1, 0.0, 0.0];
    let dt = 0.2 / 4096.0;
    for _ in 0..4096 {
        let a = derivative(reference);
        let b = derivative(core::array::from_fn(|i| reference[i] + 0.5 * dt * a[i]));
        let c = derivative(core::array::from_fn(|i| reference[i] + 0.5 * dt * b[i]));
        let d = derivative(core::array::from_fn(|i| reference[i] + dt * c[i]));
        for i in 0..5 {
            reference[i] += dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    let flow = PreparedExactZoh::new(target.system(), 0.2, budget()).unwrap();
    let evolved = flow.step(&transferred.state_after, &[1.2]).unwrap();
    close(evolved.x[0], reference[0], 2e-13);
    close(
        target.velocities(&evolved.x).unwrap()[0],
        reference[1],
        2e-13,
    );
    close(evolved.x[2] / 4.0, reference[2], 2e-13);
    close(evolved.dissipated, reference[3], 2e-12);
    close(evolved.supplied, reference[4], 2e-13);
}

#[test]
fn degenerate_rotation_and_branch_permutation_commute_with_full_damped_evolution() {
    // The relaxed spectrum is repeated. Equal-time arms still carry separate
    // histories, including the completely dark zero-projection branch.
    let previous = LinearMechanicalSystem::new(
        &[2.0, 0.0, 0.0, 2.0],
        &[0.4, 0.3, -0.1, 0.2],
        &[9.0, 0.0, 0.0, 9.0],
        &[1.0, -0.2],
        2,
        1,
        vec![
            arm(&[1.0, 0.0], 4.0, 0.2),
            arm(&[0.0, 0.0], 9.0, 0.2),
            arm(&[0.0, 1.0], 1.0, 0.7),
        ],
    )
    .unwrap();
    // A*C*A^T retains the original symmetric loss and skew gyroscopic part.
    let a = [0.6, -0.8, 0.8, 0.6];
    let target = LinearMechanicalSystem::new(
        &[2.0, 0.0, 0.0, 2.0],
        &[0.176, 0.268, -0.132, 0.424],
        &[9.0, 0.0, 0.0, 9.0],
        &[0.76, 0.68],
        2,
        1,
        vec![
            arm(&[-0.8, 0.6], 1.0, 0.7),
            arm(&[0.6, 0.8], 4.0, 0.2),
            arm(&[0.0, 0.0], 9.0, 0.2),
        ],
    )
    .unwrap();
    let plan = previous
        .prepare_state_transfer(
            &target,
            MechanicalCoordinateMap::TargetFromPrevious(&a),
            &[2, 0, 1],
            controls(),
            || false,
        )
        .unwrap();
    let state = previous
        .state(&[0.3, -0.2], &[0.1, 0.4], &[0.12, 0.07, -0.1])
        .unwrap();
    let transferred = plan.transfer(&state, || false).unwrap();
    close(transferred.state_after[0], 0.34, 2e-16);
    close(transferred.state_after[1], 0.12, 2e-16);
    close(transferred.state_after[2], -0.52, 2e-16);
    close(transferred.state_after[3], 0.64, 2e-16);
    assert_eq!(transferred.state_after[4..], [state[6], state[4], state[5]]);
    close(transferred.report.parameter_work_j, 0.0, 2e-15);
    close(transferred.report.energy_residual_j, 0.0, 2e-15);

    let before_flow = PreparedExactZoh::new(previous.system(), 0.13, budget()).unwrap();
    let after_flow = PreparedExactZoh::new(target.system(), 0.13, budget()).unwrap();
    let before_step = before_flow.step(&state, &[0.7]).unwrap();
    let after_step = after_flow.step(&transferred.state_after, &[0.7]).unwrap();
    let mapped_step = plan.transfer(&before_step.x, || false).unwrap();
    for (got, want) in after_step.x.iter().zip(&mapped_step.state_after) {
        close(*got, *want, 2e-13);
    }
    close(after_step.supplied, before_step.supplied, 2e-13);
    close(after_step.dissipated, before_step.dissipated, 2e-13);
    close(after_step.x[6], state[5] * (-0.13_f64 / 0.2).exp(), 2e-13);
}

#[test]
fn nonorthogonal_coordinates_preserve_virtual_work_full_inertia_and_history() {
    let previous = LinearMechanicalSystem::new(
        &[1.0, 0.0, 0.0, 1.0],
        &[0.0; 4],
        &[4.0, 0.0, 0.0, 9.0],
        &[1.0, -0.25],
        2,
        1,
        vec![arm(&[1.0, 0.5], 4.0, 0.3)],
    )
    .unwrap();
    // A=[2,1/2;0,1/4], B=A^-T=[1/2,0;-1,4]. Exact binary coefficients.
    let target = LinearMechanicalSystem::new(
        &[0.25, -0.5, -0.5, 17.0],
        &[0.0; 4],
        &[1.0, -2.0, -2.0, 148.0],
        &[0.5, -2.0],
        2,
        1,
        vec![arm(&[0.5, 1.0], 4.0, 0.3)],
    )
    .unwrap();
    let plan = previous
        .prepare_state_transfer(
            &target,
            MechanicalCoordinateMap::TargetFromPrevious(&[2.0, 0.5, 0.0, 0.25]),
            &[0],
            controls(),
            || false,
        )
        .unwrap();
    let state = previous.state(&[0.1, -0.2], &[0.3, 0.4], &[-0.05]).unwrap();
    let result = plan.transfer(&state, || false).unwrap();
    for (got, want) in result.state_after.iter().zip([0.1, -0.05, 0.15, 1.3, -0.1]) {
        close(*got, want, 2e-15);
    }
    let velocity = target.velocities(&result.state_after).unwrap();
    close(velocity[0], 0.8, 2e-15);
    close(velocity[1], 0.1, 2e-15);
    close(
        result.state_after[0] * result.state_after[2]
            + result.state_after[1] * result.state_after[3],
        state[0] * state[2] + state[1] * state[3],
        2e-16,
    );
    close(result.report.parameter_work_j, 0.0, 2e-15);
    close(result.report.storage_change_j, 0.0, 2e-15);
    let mut ill_conditioned = controls();
    ill_conditioned.maximum_condition_number = 1.0;
    assert!(
        previous
            .prepare_state_transfer(
                &target,
                MechanicalCoordinateMap::TargetFromPrevious(&[2.0, 0.5, 0.0, 0.25]),
                &[0],
                ill_conditioned,
                || false,
            )
            .is_err()
    );
}

#[test]
fn strong_coordinate_scaling_is_preserved_without_normalization_or_rank_loss() {
    let a = 2.0_f64.powi(200);
    let inverse = 2.0_f64.powi(-200);
    let previous = LinearMechanicalSystem::new(
        &[1.0],
        &[0.0],
        &[4.0],
        &[1.0],
        1,
        1,
        vec![arm(&[1.0], 4.0, 0.3)],
    )
    .unwrap();
    let target = LinearMechanicalSystem::new(
        &[inverse * inverse],
        &[0.0],
        &[4.0 * inverse * inverse],
        &[inverse],
        1,
        1,
        vec![arm(&[inverse], 4.0, 0.3)],
    )
    .unwrap();
    let plan = previous
        .prepare_state_transfer(
            &target,
            MechanicalCoordinateMap::TargetFromPrevious(&[a]),
            &[0],
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(plan.condition_number(), 1.0);
    let state = previous.state(&[0.25], &[0.5], &[0.125]).unwrap();
    let result = plan.transfer(&state, || false).unwrap();
    assert_eq!(result.state_after, [0.25 * a, 0.5 * inverse, 0.25]);
    assert_eq!(result.report.parameter_work_j, 0.0);
    assert_eq!(result.report.storage_change_j, 0.0);
    assert_eq!(result.report.energy_residual_j, 0.0);
}

#[test]
fn weak_parameter_work_survives_large_unchanged_storage() {
    let previous = LinearMechanicalSystem::new(
        &[1.0, 0.0, 0.0, 1.0],
        &[0.0; 4],
        &[1e100, 0.0, 0.0, 2.0],
        &[],
        2,
        0,
        vec![],
    )
    .unwrap();
    let target = LinearMechanicalSystem::new(
        &[1.0, 0.0, 0.0, 1.0],
        &[0.0; 4],
        &[1e100, 0.0, 0.0, 4.0],
        &[],
        2,
        0,
        vec![],
    )
    .unwrap();
    let state = previous.state(&[1.0, 1.0], &[0.0, 0.0], &[]).unwrap();
    let result = previous
        .prepare_state_transfer(
            &target,
            MechanicalCoordinateMap::Unchanged,
            &[],
            controls(),
            || false,
        )
        .unwrap()
        .transfer(&state, || false)
        .unwrap();
    assert_eq!(
        result.report.previous_stored_energy_j,
        result.report.target_stored_energy_j
    );
    assert_eq!(result.report.elastic_parameter_work_j, 1.0);
    assert_eq!(result.report.parameter_work_j, 1.0);
    assert_eq!(result.report.storage_change_j, 1.0);
    assert_eq!(result.report.energy_residual_j, 0.0);
}

#[test]
fn cancellation_and_refused_laws_preserve_output_and_workspace_can_resume() {
    let previous = scalar(2.0, 3.0, 4.0, 0.5);
    let target = scalar(8.0, 7.0, 16.0, 0.1);
    let state = previous.state(&[0.3], &[0.4], &[0.1]).unwrap();
    let before = state.clone();
    let plan = previous
        .prepare_state_transfer(
            &target,
            MechanicalCoordinateMap::Unchanged,
            &[0],
            controls(),
            || false,
        )
        .unwrap();
    let mut scratch = plan.workspace();
    let mut output = vec![11.0; 3];
    let mut calls = 0;
    scratch
        .transfer_into(&state, &mut output, || {
            calls += 1;
            false
        })
        .unwrap();
    let accepted = output.clone();
    output.fill(11.0);
    let mut count = 0;
    assert_eq!(
        scratch
            .transfer_into(&state, &mut output, || {
                count += 1;
                count == calls
            })
            .unwrap_err(),
        MechanicalTransferError::Cancelled
    );
    assert_eq!(output, [11.0; 3]);
    scratch
        .transfer_into(&state, &mut output, || false)
        .unwrap();
    assert_eq!(output, accepted);
    assert_eq!(state, before);
    assert!(matches!(
        previous.prepare_state_transfer(
            &target,
            MechanicalCoordinateMap::Unchanged,
            &[0],
            controls(),
            || true
        ),
        Err(MechanicalTransferError::Cancelled)
    ));

    let mut low_work = controls();
    low_work.maximum_absolute_parameter_work_j = 0.29;
    let bounded = previous
        .prepare_state_transfer(
            &target,
            MechanicalCoordinateMap::Unchanged,
            &[0],
            low_work,
            || false,
        )
        .unwrap();
    output.fill(11.0);
    assert!(
        bounded
            .workspace()
            .transfer_into(&state, &mut output, || false)
            .is_err()
    );
    assert_eq!(output, [11.0; 3]);
    let mut low_memory = controls();
    low_memory.maximum_payload_bytes = 0;
    assert!(
        previous
            .prepare_state_transfer(
                &target,
                MechanicalCoordinateMap::Unchanged,
                &[0],
                low_memory,
                || false
            )
            .is_err()
    );
    assert!(
        previous
            .prepare_state_transfer(
                &target,
                MechanicalCoordinateMap::TargetFromPrevious(&[0.0]),
                &[0],
                controls(),
                || false
            )
            .is_err()
    );
    assert!(
        previous
            .prepare_state_transfer(
                &target,
                MechanicalCoordinateMap::Unchanged,
                &[1],
                controls(),
                || false
            )
            .is_err()
    );
    let changed_strain = LinearMechanicalSystem::new(
        &[2.0],
        &[0.3],
        &[3.0],
        &[1.0],
        1,
        1,
        vec![arm(&[2.0], 4.0, 0.5)],
    )
    .unwrap();
    assert!(
        previous
            .prepare_state_transfer(
                &changed_strain,
                MechanicalCoordinateMap::Unchanged,
                &[0],
                controls(),
                || false
            )
            .is_err()
    );
}

#[test]
fn pure_relaxation_change_retains_signed_zero_and_independent_dark_history() {
    let body = |time| {
        LinearMechanicalSystem::new(
            &[1.0],
            &[0.0],
            &[0.0],
            &[],
            1,
            0,
            vec![arm(&[0.0], 4.0, time), arm(&[0.0], 9.0, time)],
        )
        .unwrap()
    };
    let previous = body(0.2);
    let target = body(0.03);
    let state = [-0.0, -0.0, 0.14, -0.18];
    let plan = previous
        .prepare_state_transfer(
            &target,
            MechanicalCoordinateMap::Unchanged,
            &[0, 1],
            controls(),
            || false,
        )
        .unwrap();
    let result = plan.transfer(&state, || false).unwrap();
    assert!(
        result
            .state_after
            .iter()
            .zip(state)
            .all(|(a, b)| a.to_bits() == b.to_bits())
    );
    assert_eq!(result.report.parameter_work_j, 0.0);
    let evolved = PreparedExactZoh::new(target.system(), 0.09, budget())
        .unwrap()
        .step(&result.state_after, &[])
        .unwrap();
    close(evolved.x[2], state[2] * (-3.0_f64).exp(), 2e-14);
    close(evolved.x[3], state[3] * (-3.0_f64).exp(), 2e-14);
    assert!(
        previous
            .prepare_state_transfer(
                &target,
                MechanicalCoordinateMap::Unchanged,
                &[0, 0],
                controls(),
                || false
            )
            .is_err()
    );
}
