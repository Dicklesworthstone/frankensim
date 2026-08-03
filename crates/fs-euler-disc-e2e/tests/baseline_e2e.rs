use std::process::Command;

use fs_euler_disc_e2e::{
    BaselineDynamicsClass, BaselineRefusalReason, BaselineRunOutput, BaselineTerminal,
    SquatDiscInput, run_ideal_conservative_baseline,
};

#[test]
fn ideal_conservative_production_slice_emits_deterministic_trajectory_and_ledger() {
    let input = SquatDiscInput::nominal();
    let first = run_ideal_conservative_baseline(input.clone());
    let second = run_ideal_conservative_baseline(input);
    assert_eq!(first, second);

    let BaselineRunOutput::Completed(trajectory) = first else {
        panic!("nominal production input must be admitted");
    };
    assert_eq!(trajectory.samples.len(), 101);
    assert_eq!(
        trajectory.terminal,
        BaselineTerminal::TimeHorizonReached {
            completed_steps: 100
        }
    );
    let terminal = trajectory.samples.last().expect("terminal sample");
    assert!(terminal.support.plane_gap_m.abs() < 1.0e-12);
    assert!(terminal.support.no_slip_residual_mps < 1.0e-12);
    assert!(
        terminal.support.required_tangential_force_n
            <= terminal.support.available_static_friction_n
    );
    assert!(terminal.energy.residual_from_initial_j.abs() < 1.0e-9);
    assert!(
        trajectory
            .equilibrium
            .precession_balance_residual_s_inv2
            .abs()
            < 1.0e-7
    );
    assert!(
        trajectory.equilibrium.small_angle_energy_residual_j.abs()
            < 0.03 * trajectory.equilibrium.small_angle_energy_oracle_j
    );
}

#[test]
fn perturbing_precession_refuses_the_dynamic_steady_rolling_oracle() {
    let mut input = SquatDiscInput::nominal();
    input.precession_rad_s *= 1.01;
    assert_eq!(
        run_ideal_conservative_baseline(input),
        BaselineRunOutput::Refused(fs_euler_disc_e2e::BaselineRefusal {
            reason: BaselineRefusalReason::DynamicOracleViolation,
        })
    );
}

#[test]
fn non_oracle_rates_are_explicitly_labeled_as_prescribed_kinematics() {
    let mut input = SquatDiscInput::nominal();
    input.dynamics_class = BaselineDynamicsClass::PrescribedKinematicPath;
    input.precession_rad_s *= 1.01;
    input.reset_supported_initial_state();
    let BaselineRunOutput::Completed(trajectory) = run_ideal_conservative_baseline(input) else {
        panic!("prescribed path should remain a kinematic, not dynamic, output");
    };
    assert_eq!(
        trajectory.model_id,
        "prescribed_no_slip_rolling_disc_kinematics"
    );
}

#[test]
fn executable_runs_the_production_baseline_and_reports_structured_output() {
    let binary = env!("CARGO_BIN_EXE_ideal_conservative_baseline");
    let output = Command::new(binary)
        .output()
        .expect("execute Euler baseline binary");
    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).expect("UTF-8 structured output");
    assert!(stdout.contains("\"model\":\"ideal_conservative_no_slip_rolling_disc\""));
    assert!(stdout.contains("\"disposition\":\"completed\""));
    assert!(stdout.contains("\"samples\":101"));
}
