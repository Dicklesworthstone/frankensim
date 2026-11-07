//! Independent original-energy and reaction-torque checks for nonlinear motors.
//! These certify the native discrete quadrature separately from the analytic
//! continuous RL transient; a balanced discrete step is not an exact trajectory.

use fs_phs::circuit::inductance::InductorFluxCurve;
use fs_phs::circuit::motor::nonlinear::{
    NonlinearDcMotor, NonlinearDcMotorParameters, NonlinearMotorControls,
};
use fs_phs::{
    InterconnectionBudget, PortHamiltonian, PreparedStepError, QuadraticStorage, SolveBudget,
    StepWorkspace, interconnect_mapped,
};

fn curve() -> InductorFluxCurve {
    InductorFluxCurve::try_new(vec![
        (-4.0, -16.0),
        (-1.0, -1.0),
        (0.0, 0.0),
        (1.0, 1.0),
        (4.0, 16.0),
    ])
    .unwrap()
}
fn parameters() -> NonlinearDcMotorParameters {
    NonlinearDcMotorParameters {
        winding_resistance_ohm: 0.4,
        rotor_inertia_kg_m2: 2.0,
        stator_inertia_kg_m2: 3.0,
        coupling_nm_per_a: 1.25,
        bearing_drag_nm_s_per_rad: 0.1,
    }
}
// Independent integral of the fixture's affine current pieces, with H(0)=0.
fn magnetic_energy(flux: f64) -> f64 {
    let f = flux.abs();
    if f <= 1.0 {
        0.5 * f * f
    } else {
        2.5 * f * f - 4.0 * f + 2.0
    }
}
fn endpoint_current(flux: f64) -> f64 {
    if flux.abs() <= 1.0 {
        flux
    } else {
        flux.signum() * (5.0 * flux.abs() - 4.0)
    }
}
fn close(got: f64, want: f64, tolerance: f64) {
    assert!(
        got.is_finite() && want.is_finite() && (got - want).abs() <= tolerance,
        "got={got:e}, want={want:e}, tolerance={tolerance:e}"
    );
}

#[test]
fn crossing_saturation_retains_reciprocal_work_both_shafts_and_positive_heat() {
    let p = parameters();
    let h = 0.3;
    for sign in [-1.0, 1.0] {
        for (start, end) in [(0.5, 1.5), (1.5, 0.5)] {
            let (flux0, flux1) = (sign * start, sign * end);
            let (wr, ws) = (sign * 1.2, sign * -0.4);
            let (tr, ts) = (sign * 0.2, sign * -0.1);
            let dh = magnetic_energy(flux1) - magnetic_energy(flux0);
            let current = dh / (flux1 - flux0);
            // Algebraic elimination of the two independent midpoint momentum
            // equations supplies an exact endpoint fixture, not a second solver.
            let inverse_inertia = 1.0 / p.rotor_inertia_kg_m2 + 1.0 / p.stator_inertia_kg_m2;
            let drag = 0.5 * h * p.bearing_drag_nm_s_per_rad * inverse_inertia;
            let relative1 = ((wr - ws) * (1.0 - drag)
                + h * (p.coupling_nm_per_a * current * inverse_inertia
                    + tr / p.rotor_inertia_kg_m2
                    - ts / p.stator_inertia_kg_m2))
                / (1.0 + drag);
            let relative_mean = 0.5 * (wr - ws + relative1);
            let internal =
                p.coupling_nm_per_a * current - p.bearing_drag_nm_s_per_rad * relative_mean;
            let speeds1 = [
                wr + h * (internal + tr) / p.rotor_inertia_kg_m2,
                ws + h * (-internal + ts) / p.stator_inertia_kg_m2,
            ];
            let means = [0.5 * (wr + speeds1[0]), 0.5 * (ws + speeds1[1])];
            let voltage = (flux1 - flux0) / h
                + p.winding_resistance_ohm * current
                + p.coupling_nm_per_a * relative_mean;
            let motor = NonlinearDcMotor::new(p, curve()).unwrap();
            let state = motor.state(endpoint_current(flux0), [wr, ws]).unwrap();
            let mut prepared = motor
                .prepare_controlled(h, NonlinearMotorControls::default(), || false)
                .unwrap();
            let (mut next, mut output) = ([0.0; 3], [0.0; 3]);
            let record = prepared
                .step_controlled(&state, &[voltage, tr, ts], &mut next, &mut output, || false)
                .unwrap();
            close(next[0], flux1, 2e-11);
            close(next[1], p.rotor_inertia_kg_m2 * speeds1[0], 2e-11);
            close(next[2], p.stator_inertia_kg_m2 * speeds1[1], 2e-11);
            close(output[0], current, 2e-11);
            assert!(
                (output[0] - 0.5 * (endpoint_current(flux0) + endpoint_current(flux1))).abs() > 0.4
            );
            close(record.magnetic_energy_change_j, dh, 3e-11);
            close(
                record.discrete_back_emf_v,
                p.coupling_nm_per_a * relative_mean,
                3e-11,
            );
            close(
                record.copper_heat_j,
                h * p.winding_resistance_ohm * current * current,
                3e-11,
            );
            close(
                record.bearing_heat_j,
                h * p.bearing_drag_nm_s_per_rad * relative_mean.powi(2),
                3e-11,
            );
            assert!(record.copper_heat_j > 0.0 && record.bearing_heat_j > 0.0);
            close(record.electrical_work_j, h * voltage * current, 3e-11);
            assert_eq!(record.electrical_work_j < 0.0, end < start);
            for i in 0..2 {
                let inertia = [p.rotor_inertia_kg_m2, p.stator_inertia_kg_m2][i];
                let speed0 = [wr, ws][i];
                let reaction = if i == 0 { 1.0 } else { -1.0 };
                let external = h * [tr, ts][i] * means[i];
                let electromagnetic = reaction * h * p.coupling_nm_per_a * current * means[i];
                let bearing =
                    -reaction * h * p.bearing_drag_nm_s_per_rad * relative_mean * means[i];
                let kinetic = 0.5 * inertia * (speeds1[i] - speed0) * (speeds1[i] + speed0);
                close(output[i + 1], means[i], 3e-11);
                close(record.shaft_angle_increment_rad[i], h * means[i], 3e-11);
                close(record.shaft_work_j[i], external, 3e-11);
                close(record.electromagnetic_work_j[i], electromagnetic, 3e-11);
                close(record.bearing_work_j[i], bearing, 3e-11);
                close(record.kinetic_energy_change_j[i], kinetic, 3e-11);
                close(kinetic, external + electromagnetic + bearing, 3e-11);
            }
            close(
                next[1] + next[2] - state[1] - state[2],
                h * (tr + ts),
                3e-11,
            );
            close(record.angular_momentum_residual_nms, 0.0, 3e-11);
            close(record.native.balance_residual(), 0.0, 3e-11);
        }
    }
}

fn driven_rl(steps: usize) -> f64 {
    let p = NonlinearDcMotorParameters {
        winding_resistance_ohm: 1.0,
        coupling_nm_per_a: 0.0,
        bearing_drag_nm_s_per_rad: 0.0,
        ..parameters()
    };
    let motor = NonlinearDcMotor::new(p, curve()).unwrap();
    // Huge unchanged kinetic energy may neither hide the small magnetic work
    // nor set the physical tolerance for its independent winding boundary.
    let mut state = motor.state(0.0, [1e7, -2e7]).unwrap();
    let h = 1.0 / steps as f64;
    let mut prepared = motor
        .prepare_controlled(h, NonlinearMotorControls::default(), || false)
        .unwrap();
    let (mut source, mut heat) = (0.0, 0.0);
    for _ in 0..steps {
        let (mut next, mut output) = ([0.0; 3], [0.0; 3]);
        let record = prepared
            .step_controlled(&state, &[3.0, 0.0, 0.0], &mut next, &mut output, || false)
            .unwrap();
        close(next[0] - state[0], h * (3.0 - output[0]), 2e-11);
        close(
            record.magnetic_energy_change_j,
            magnetic_energy(next[0]) - magnetic_energy(state[0]),
            2e-11,
        );
        close(record.electrical_work_j, h * 3.0 * output[0], 2e-11);
        close(record.copper_heat_j, h * output[0] * output[0], 2e-11);
        assert_eq!(&next[1..], &state[1..]);
        assert_eq!(record.kinetic_energy_change_j, [0.0, 0.0]);
        assert_eq!(record.bearing_heat_j, 0.0);
        source += record.electrical_work_j;
        heat += record.copper_heat_j;
        state = next;
    }
    close(magnetic_energy(state[0]) + heat, source, 2e-9);
    state[0]
}

#[test]
fn continuous_rl_saturation_refines_while_discrete_work_already_closes() {
    // λ'=3-I(λ): first 3(1-exp(-t)), then 1.4-.4exp(-5(t-ln(1.5))).
    let crossing = 1.5_f64.ln();
    let exact = 1.4 - 0.4 * (-5.0 * (1.0 - crossing)).exp();
    let errors = [40, 80, 160].map(|steps| (driven_rl(steps) - exact).abs());
    assert!(errors[2] < 1e-5, "{errors:?}");
    for pair in errors.windows(2) {
        assert!((3.7..4.3).contains(&(pair[0] / pair[1])), "{errors:?}");
    }
}

#[test]
fn isolated_nonlinear_motor_exchanges_energy_with_finite_stator_reaction() {
    let p = NonlinearDcMotorParameters {
        winding_resistance_ohm: 0.0,
        rotor_inertia_kg_m2: 1.0,
        stator_inertia_kg_m2: 3.0,
        coupling_nm_per_a: 2.0,
        bearing_drag_nm_s_per_rad: 0.0,
    };
    let motor = NonlinearDcMotor::new(p, curve()).unwrap();
    let mut state = motor.state(endpoint_current(1.4), [0.0, 0.0]).unwrap();
    let initial_energy = magnetic_energy(state[0]);
    let mut prepared = motor
        .prepare_controlled(0.01, NonlinearMotorControls::default(), || false)
        .unwrap();
    let mut negative_flux = false;
    let mut moving_stator = false;
    for _ in 0..400 {
        let (mut next, mut output) = ([0.0; 3], [0.0; 3]);
        let record = prepared
            .step_controlled(&state, &[0.0; 3], &mut next, &mut output, || false)
            .unwrap();
        close(next[1] + next[2], 0.0, 2e-10);
        let energy = magnetic_energy(next[0])
            + 0.5 * next[1] * next[1] / p.rotor_inertia_kg_m2
            + 0.5 * next[2] * next[2] / p.stator_inertia_kg_m2;
        close(energy, initial_energy, 2e-9);
        close(
            record.magnetic_energy_change_j + record.kinetic_energy_change_j.iter().sum::<f64>(),
            0.0,
            2e-10,
        );
        assert_eq!(record.electrical_work_j, 0.0);
        assert_eq!(record.copper_heat_j, 0.0);
        assert_eq!(record.bearing_heat_j, 0.0);
        negative_flux |= next[0] < -1.0;
        moving_stator |= output[2].abs() > 0.1;
        state = next;
    }
    assert!(negative_flux && moving_stator);
}

#[test]
fn finite_capacitor_bus_retains_nonlinear_magnetic_work_after_native_composition() {
    let capacitance = 2.0;
    let bus = PortHamiltonian::new(
        1,
        1,
        vec![0.0],
        vec![0.0],
        vec![1.0],
        Box::new(QuadraticStorage::new(vec![1.0 / capacitance], 1).unwrap()),
    )
    .unwrap();
    let p = NonlinearDcMotorParameters {
        winding_resistance_ohm: 0.1,
        coupling_nm_per_a: 0.5,
        bearing_drag_nm_s_per_rad: 0.02,
        ..parameters()
    };
    let motor = NonlinearDcMotor::new(p, curve()).unwrap();
    let system = interconnect_mapped(
        bus,
        motor.into_system(),
        &[0],
        &[0],
        &[1.0],
        InterconnectionBudget {
            max_states: 4,
            max_total_ports: 4,
            max_dense_bytes: 1 << 20,
            max_coupling_multiply_terms: 1_000_000,
        },
    )
    .unwrap();
    assert_eq!(system.state_dim(), 4);
    assert_eq!(system.port_dim(), 2);
    let mut workspace = StepWorkspace::new(&system).unwrap();
    workspace.set_solve_budget(SolveBudget::default()).unwrap();
    let h = 0.02;
    let mut state = [4.0, 0.0, 0.0, 0.0];
    let mut heat = 0.0;
    let mut saturation = false;
    let mut regeneration = false;
    for _ in 0..400 {
        let (mut next, mut outputs) = ([0.0; 4], [0.0; 2]);
        let record = workspace
            .step_into_controlled(
                &system,
                &state,
                &[0.0; 2],
                h,
                &mut next,
                &mut outputs,
                || false,
            )
            .unwrap();
        let magnetic = magnetic_energy(next[1]) - magnetic_energy(state[1]);
        let current = if next[1] == state[1] {
            endpoint_current(state[1])
        } else {
            magnetic / (next[1] - state[1])
        };
        let voltage = 0.5 * (state[0] + next[0]) / capacitance;
        let relative = outputs[0] - outputs[1];
        let bus_change = 0.5 * (next[0] - state[0]) * (next[0] + state[0]) / capacitance;
        let copper = h * p.winding_resistance_ohm * current * current;
        let bearing = h * p.bearing_drag_nm_s_per_rad * relative * relative;
        close(next[0] - state[0], -h * current, 2e-9);
        close(
            next[1] - state[1],
            h * (voltage - p.winding_resistance_ohm * current - p.coupling_nm_per_a * relative),
            2e-9,
        );
        close(bus_change, -h * voltage * current, 2e-9);
        close(record.dissipated, copper + bearing, 2e-9);
        close(next[2] + next[3], 0.0, 2e-9);
        assert_eq!(record.supplied, 0.0);
        heat += copper + bearing;
        let energy = 0.5 * next[0] * next[0] / capacitance
            + magnetic_energy(next[1])
            + 0.5 * next[2] * next[2] / p.rotor_inertia_kg_m2
            + 0.5 * next[3] * next[3] / p.stator_inertia_kg_m2;
        close(energy + heat, 4.0, 2e-8);
        saturation |= next[1].abs() > 1.0;
        regeneration |= bus_change > 1e-5;
        state = next;
    }
    assert!(saturation && regeneration);
}

#[test]
fn invalid_source_budget_and_late_cancellation_publish_neither_output() {
    let p = parameters();
    let motor = || NonlinearDcMotor::new(p, curve()).unwrap();
    assert!(
        NonlinearDcMotor::new(
            NonlinearDcMotorParameters {
                stator_inertia_kg_m2: 0.0,
                ..p
            },
            curve()
        )
        .is_err()
    );
    assert!(motor().state(17.0, [0.0; 2]).is_err());
    assert!(
        motor()
            .prepare_controlled(0.0, NonlinearMotorControls::default(), || false)
            .is_err()
    );
    assert!(
        motor()
            .prepare_controlled(
                0.1,
                NonlinearMotorControls {
                    absolute_energy_tolerance_j: f64::NAN,
                    ..NonlinearMotorControls::default()
                },
                || false
            )
            .is_err()
    );
    assert!(
        motor()
            .prepare_controlled(0.1, NonlinearMotorControls::default(), || true)
            .is_err()
    );
    let state = motor().state(0.5, [0.2, -0.1]).unwrap();
    let input = [2.0, 0.1, -0.05];
    let mut limited = motor()
        .prepare_controlled(
            0.1,
            NonlinearMotorControls {
                solve: SolveBudget {
                    max_iterations: 0,
                    ..SolveBudget::default()
                },
                ..NonlinearMotorControls::default()
            },
            || false,
        )
        .unwrap();
    let sentinel = ([7.0, 8.0, 9.0], [10.0, 11.0, 12.0]);
    let (mut next, mut output) = sentinel;
    assert!(
        limited
            .step_controlled(&state, &input, &mut next, &mut output, || false)
            .is_err()
    );
    assert_eq!((next, output), sentinel);
    let mut prepared = motor()
        .prepare_controlled(0.1, NonlinearMotorControls::default(), || false)
        .unwrap();
    // Deliberately leave the bounded constitutive source during a native probe.
    assert!(
        prepared
            .step_controlled(&state, &[1e6, 0.0, 0.0], &mut next, &mut output, || false)
            .is_err()
    );
    assert_eq!((next, output), sentinel);
    let mut polls = 0;
    let expected_record = prepared
        .step_controlled(&state, &input, &mut next, &mut output, || {
            polls += 1;
            false
        })
        .unwrap();
    let expected = (next, output);
    assert!(polls > 1);
    (next, output) = sentinel;
    let mut count = 0;
    let cancelled = prepared.step_controlled(&state, &input, &mut next, &mut output, || {
        count += 1;
        count == polls
    });
    assert!(matches!(cancelled, Err(PreparedStepError::Cancelled)));
    assert_eq!((next, output), sentinel);
    let retry = prepared
        .step_controlled(&state, &input, &mut next, &mut output, || false)
        .unwrap();
    assert_eq!((next, output), expected);
    assert_eq!(retry.native, expected_record.native);
    assert_eq!(
        retry.magnetic_energy_change_j,
        expected_record.magnetic_energy_change_j
    );
    assert_eq!(retry.copper_heat_j, expected_record.copper_heat_j);
    assert_eq!(retry.bearing_heat_j, expected_record.bearing_heat_j);
}
