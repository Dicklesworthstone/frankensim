//! Independent winding/Cayley, recoil, added-inertia and source-energy oracles.
use super::*;
use fs_flyer::partitioned::articulated::actuation::{
    JointVoltageDrive, VoltageDrivenState, step_with_voltage_drives,
};
use fs_mbd::rotordynamics::transient::drive::voltage::{RotorVoltageControl, RotorWindingMotor};
use fs_mbd::rotordynamics::transient::drive::{RotorCurrentDrive, RotorDrivePowerControl};

fn drive(joint: usize, voltage: f64, inductance: f64, resistance: f64) -> JointVoltageDrive {
    JointVoltageDrive {
        motor: RotorWindingMotor {
            drive: RotorCurrentDrive {
                shaft_coordinate: joint,
                torque_constant_nm_per_a: 0.5,
                winding_resistance_ohm: resistance,
                maximum_current_a: 100.0,
                maximum_terminal_voltage_v: 20.0,
                minimum_speed_rad_s: -100.0,
                maximum_speed_rad_s: 100.0,
                motoring_efficiency: 1.0,
                regenerative_efficiency: 1.0,
            },
            inductance_h: inductance,
        },
        command_voltage_v: voltage,
        control: RotorVoltageControl {
            source: RotorDrivePowerControl {
                maximum_source_input_w: 10_000.0,
                maximum_source_return_w: 10_000.0,
                absolute_power_tolerance_w: 1e-10,
                absolute_voltage_tolerance_v: 1e-10,
            },
            absolute_current_tolerance_a: 1e-12,
        },
    }
}

fn single_rotor() -> ArticulatedModel {
    ArticulatedModel::new(vec![
        Link::new(
            "carrier",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(4.0, [0.0; 3], [3.0, 3.0, 4.0]),
        ),
        Link::new(
            "hinge",
            Some(0),
            Se3::identity(),
            JointModel::revolute(vector([0.0, 0.0, 1.0]), None).unwrap(),
            inertia(2.0, [0.0; 3], [0.8, 0.8, 1.0]),
        ),
    ])
    .unwrap()
}

fn single_state(current: f64, rate: f64) -> VoltageDrivenState {
    let mut mechanical = rotor_state();
    mechanical.motion.joint_position = vec![0.2];
    mechanical.motion.joint_velocity = vec![rate];
    mechanical.motion.base.twist_body.angular.z = -rate / 5.0;
    VoltageDrivenState {
        mechanical,
        currents_a: vec![current],
    }
}

fn one_config() -> ArticulatedStepConfig<'static> {
    ArticulatedStepConfig {
        joint_velocity_tolerances: &[1e-12],
        ..config()
    }
}

#[test]
fn winding_cayley_retains_airframe_recoil_and_stiff_electrical_storage() {
    let model = single_rotor();
    // Ordinary powered winding, very stiff short circuit, and lossless LC-like
    // mechanical/magnetic exchange. No native mass query supplies this oracle.
    for (l, r, voltage, dt) in [
        (0.35, 0.4, 0.8, 0.04),
        (1e-6, 20.0, 0.0, 0.02),
        (0.15, 0.0, 0.0, 0.1),
    ] {
        let mut state = single_state(0.3, 0.6);
        let motor = drive(0, voltage, l, r);
        let k = motor.motor.drive.torque_constant_nm_per_a;
        let reduced_inertia = 4.0 / 5.0;
        let h = dt / 2.0;
        for tick in 0..12 {
            let w0 = state.mechanical.motion.joint_velocity[0];
            let i0 = state.currents_a[0];
            let im =
                (l * i0 + h * voltage - h * k * w0) / (l + h * r + h * h * k * k / reduced_inertia);
            let w1 = w0 + dt * k * im / reduced_inertia;
            let i1 = 2.0 * im - i0;
            let out = step_with_voltage_drives(
                &model,
                &[ZERO_IMPEDANCE],
                &state,
                &[],
                &[motor],
                tick as f64 * dt,
                dt,
                one_config(),
                |_, _, map| loads(map, &[], 1.0),
                || false,
            )
            .unwrap();
            close(out.mechanical.state.motion.joint_velocity[0], w1, 3e-12);
            close(
                out.mechanical.state.motion.base.twist_body.angular.z,
                -w1 / 5.0,
                3e-12,
            );
            close(
                out.mechanical.state.motion.joint_position[0],
                state.mechanical.motion.joint_position[0] + h * (w0 + w1),
                3e-12,
            );
            let receipt = &out.drives[0];
            close(receipt.current_after_a, i1, 3e-12);
            close(receipt.midpoint_current_a, im, 3e-12);
            close(
                receipt.joint_work_j,
                0.5 * reduced_inertia * (w1 * w1 - w0 * w0),
                3e-12,
            );
            close(
                receipt.magnetic_energy_change_j,
                0.5 * l * (i1 * i1 - i0 * i0),
                3e-12,
            );
            close(
                receipt.terminal_work_j,
                receipt.joint_work_j + receipt.magnetic_energy_change_j + receipt.winding_heat_j,
                3e-12,
            );
            state = out.into_state();
        }
    }
}

#[test]
fn two_windings_share_native_fluid_inertia_and_carrier_reaction() {
    let model = rotor_tree(None);
    let rho = 1.2;
    let dimensions = [(1.3_f64, 1.8_f64), (0.9_f64, 1.2_f64)];
    let geometry: Vec<_> = dimensions
        .iter()
        .enumerate()
        .map(|(i, (c, s))| {
            (
                i + 1,
                plate(*c, *s, [0.0; 3], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]),
            )
        })
        .collect();
    let strips: Vec<_> = geometry
        .iter()
        .map(|(i, p)| LinkedStrip::new(*i, p.clone()).unwrap())
        .collect();
    let j: [f64; 2] = core::array::from_fn(|i| {
        [1.0, 0.5][i]
            + rho * std::f64::consts::PI * dimensions[i].0.powi(4) * dimensions[i].1 / 128.0
    });
    let total = 4.0 + j[0] + j[1];
    let a = j[0] - j[0] * j[0] / total;
    let b = -j[0] * j[1] / total;
    let c = j[1] - j[1] * j[1] / total;
    let determinant = a * c - b * b;
    let inverse = [
        [c / determinant, -b / determinant],
        [-b / determinant, a / determinant],
    ];
    let motors = [drive(0, 1.3, 0.2, 0.3), drive(1, -0.6, 0.1, 0.2)];
    let mut state = VoltageDrivenState {
        mechanical: rotor_state(),
        currents_a: vec![0.0; 2],
    };
    let dt = 0.06;
    let h = dt / 2.0;
    for tick in 0..8 {
        let before = complete_invariants(&model, &state.mechanical, &geometry, rho);
        // Eliminate the two mechanical rates in the independent simultaneous
        // circuit equations, retaining the off-diagonal electrical feedback.
        let aa = 0.2 + h * 0.3 + h * h * 0.25 * inverse[0][0];
        let bb = h * h * 0.25 * inverse[0][1];
        let cc = 0.1 + h * 0.2 + h * h * 0.25 * inverse[1][1];
        let rhs = [
            0.2 * state.currents_a[0] + h * 1.3
                - h * 0.5 * state.mechanical.motion.joint_velocity[0],
            0.1 * state.currents_a[1]
                - h * 0.6
                - h * 0.5 * state.mechanical.motion.joint_velocity[1],
        ];
        let det = aa * cc - bb * bb;
        let im = [
            (cc * rhs[0] - bb * rhs[1]) / det,
            (aa * rhs[1] - bb * rhs[0]) / det,
        ];
        let want: [f64; 2] = core::array::from_fn(|i| {
            state.mechanical.motion.joint_velocity[i]
                + dt * 0.5 * (inverse[i][0] * im[0] + inverse[i][1] * im[1])
        });
        let out = step_with_voltage_drives(
            &model,
            &[ZERO_IMPEDANCE; 2],
            &state,
            &[],
            &motors,
            tick as f64 * dt,
            dt,
            config(),
            |_, _, map| loads(map, &strips, rho),
            || false,
        )
        .unwrap();
        for i in 0..2 {
            close(
                out.mechanical.state.motion.joint_velocity[i],
                want[i],
                4e-12,
            );
            close(out.drives[i].midpoint_current_a, im[i], 4e-12);
        }
        close(
            out.mechanical.state.motion.base.twist_body.angular.z,
            -(j[0] * want[0] + j[1] * want[1]) / total,
            4e-12,
        );
        let after = complete_invariants(&model, &out.mechanical.state, &geometry, rho);
        close(
            after.0 - before.0,
            out.drives.iter().map(|r| r.joint_work_j).sum(),
            4e-12,
        );
        close((after.2 - before.2).norm(), 0.0, 4e-12);
        state = out.into_state();
    }
}

#[test]
fn actual_link_load_and_voltage_saturation_enter_the_same_accepted_stage() {
    let model = single_rotor();
    let state = single_state(0.2, 0.6);
    let mut motor = drive(0, 10.0, 0.3, 0.4);
    motor.motor.drive.maximum_terminal_voltage_v = 2.0;
    let (dt, h, k, root_load, hinge_load) = (0.08, 0.04, 0.5, -0.2, -0.45);
    let im = (0.3 * 0.2 + h * 2.0 - h * k * (0.6 + h * (hinge_load - root_load / 4.0)))
        / (0.3 + h * 0.4 + h * h * k * k / 0.8);
    let want_base = -0.6 / 5.0 + dt * (root_load - k * im) / 4.0;
    let want_rate = 0.6 + dt * (k * im / 0.8 + hinge_load - root_load / 4.0);
    let calls = Cell::new(0);
    let out = step_with_voltage_drives(
        &model,
        &[ZERO_IMPEDANCE],
        &state,
        &[],
        &[motor],
        0.0,
        dt,
        one_config(),
        |_, _, map| {
            calls.set(calls.get() + 1);
            let mut result = loads(map, &[], 1.0)?;
            result.world_wrenches[0].torque.z = root_load;
            result.world_wrenches[1].torque.z = hinge_load;
            Ok(result)
        },
        || false,
    )
    .unwrap();
    assert_eq!(calls.get(), out.mechanical.iterations);
    close(
        out.mechanical.state.motion.base.twist_body.angular.z,
        want_base,
        3e-12,
    );
    close(
        out.mechanical.state.motion.joint_velocity[0],
        want_rate,
        3e-12,
    );
    let receipt = &out.drives[0];
    assert_eq!(receipt.applied_voltage_v, 2.0);
    close(receipt.midpoint_current_a, im, 3e-12);
    let base_mid = 0.5 * (-0.6 / 5.0 + want_base);
    let rotor_mid = base_mid + 0.5 * (0.6 + want_rate);
    let before = complete_invariants(&model, &state.mechanical, &[], 1.0).0;
    let after = complete_invariants(&model, &out.mechanical.state, &[], 1.0).0;
    close(
        after - before,
        receipt.joint_work_j + dt * (root_load * base_mid + hinge_load * rotor_mid),
        3e-12,
    );
}

#[test]
fn reconstructed_current_reversal_keeps_both_source_transfers_and_bridge_heat() {
    let model = single_rotor();
    let state = single_state(-1.0, 0.3);
    let mut motor = drive(0, 4.0, 0.2, 0.1);
    motor.motor.drive.motoring_efficiency = 0.8;
    motor.motor.drive.regenerative_efficiency = 0.6;
    let dt = 0.2;
    let out = step_with_voltage_drives(
        &model,
        &[ZERO_IMPEDANCE],
        &state,
        &[],
        &[motor],
        0.0,
        dt,
        one_config(),
        |_, _, map| loads(map, &[], 1.0),
        || false,
    )
    .unwrap();
    let receipt = &out.drives[0];
    let p0 = 4.0 * state.currents_a[0];
    let p1 = 4.0 * receipt.current_after_a;
    assert!(p0 < 0.0 && p1 > 0.0);
    // Exact signed triangle areas for the declared affine current segment.
    let positive = 0.5 * p1 * p1 / (p1 - p0);
    let negative = 0.5 * p0 * p0 / (p1 - p0);
    close(receipt.electrical_draw_j, dt * positive / 0.8, 3e-12);
    close(receipt.electrical_return_j, dt * negative * 0.6, 3e-12);
    close(
        receipt.bridge_heat_j,
        dt * (positive * 0.25 + negative * 0.4),
        3e-12,
    );
    assert!(receipt.bridge_heat_j > 0.0);
    close(receipt.electrical_balance_residual_j, 0.0, 3e-12);
    close(receipt.winding_residual_v, 0.0, 3e-12);
}

#[test]
fn infeasible_current_source_and_authored_impedance_refuse_without_partial_state() {
    let model = single_rotor();
    let state = single_state(0.0, 0.0);
    let original = state.clone();
    let motor = drive(0, 2.0, 0.2, 0.1);
    for source_limit in [false, true] {
        let mut bounded = motor;
        if source_limit {
            bounded.control.source.maximum_source_input_w = 0.01;
        } else {
            bounded.motor.drive.maximum_current_a = 0.01;
        }
        let err = step_with_voltage_drives(
            &model,
            &[ZERO_IMPEDANCE],
            &state,
            &[],
            &[bounded],
            0.0,
            0.1,
            one_config(),
            |_, _, map| loads(map, &[], 1.0),
            || false,
        )
        .unwrap_err();
        assert_eq!(
            err.code,
            if source_limit {
                "articulated-drive-power-limit"
            } else {
                "articulated-drive-current-limit"
            }
        );
        assert_eq!(state, original);
    }
    let invalid = [JointImpedance {
        damping: -0.001,
        ..ZERO_IMPEDANCE
    }];
    let calls = Cell::new(0);
    let err = step_with_voltage_drives(
        &model,
        &invalid,
        &state,
        &[],
        &[motor],
        0.0,
        0.1,
        one_config(),
        |_, _, map| {
            calls.set(calls.get() + 1);
            loads(map, &[], 1.0)
        },
        || false,
    )
    .unwrap_err();
    assert_eq!(err.code, "articulated-drive-invalid");
    assert_eq!(calls.get(), 0);
    let cancel = Cell::new(false);
    assert!(
        step_with_voltage_drives(
            &model,
            &[ZERO_IMPEDANCE],
            &state,
            &[],
            &[motor],
            0.0,
            0.1,
            one_config(),
            |_, _, map| {
                cancel.set(true);
                loads(map, &[], 1.0)
            },
            || cancel.get()
        )
        .is_err()
    );
    assert_eq!(state, original);
}
