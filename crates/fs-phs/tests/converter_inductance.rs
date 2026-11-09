//! Nonlinear converter current, finite-bus energy and physical conduction heat.
//! Discrete work identities and continuous-time refinement are checked
//! separately; neither fixture claims PWM ripple or semiconductor switching loss.

use fs_phs::circuit::converter::inductance::{
    NonlinearAveragedBuckBoost, NonlinearConverterControls,
};
use fs_phs::circuit::converter::{AveragedBuckBoost, BridgeModulation, BuckBoostParameters};
use fs_phs::circuit::inductance::InductorFluxCurve;
use fs_phs::linear_exact::LinearSplitBudget;
use fs_phs::{PhsError, PreparedStepError, SolveBudget, StepWorkspace};

#[path = "converter_inductance/power_load.rs"]
mod power_load;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        actual.is_finite() && expected.is_finite() && (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

fn parameters() -> BuckBoostParameters {
    BuckBoostParameters {
        inductance_h: 1.0,
        bus_capacitance_f: [2.0, 3.0],
        winding_resistance_ohm: 0.03,
        input_on_resistance_ohm: [0.02, 0.04],
        output_on_resistance_ohm: [0.015, 0.025],
    }
}

fn modulation(a: f64, b: f64) -> BridgeModulation {
    BridgeModulation {
        input_high_fraction: a,
        output_high_fraction: b,
    }
}

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

fn nonlinear(p: BuckBoostParameters, m: BridgeModulation) -> NonlinearAveragedBuckBoost {
    AveragedBuckBoost::new(p, m)
        .unwrap()
        .with_inductor_curve(curve())
        .unwrap()
}

// Independent antiderivative of I(lambda)=lambda inside |lambda|<=1,
// and sign(lambda)*(5|lambda|-4) outside the saturation knee.
fn magnetic_energy(flux: f64) -> f64 {
    let x = flux.abs();
    if x <= 1.0 {
        0.5 * x * x
    } else {
        2.5 * x * x - 4.0 * x + 2.0
    }
}

fn current(flux: f64) -> f64 {
    if flux.abs() <= 1.0 {
        flux
    } else {
        flux.signum() * (5.0 * flux.abs() - 4.0)
    }
}

fn chord_current(before: f64, after: f64) -> f64 {
    if before == after {
        current(before)
    } else {
        (magnetic_energy(after) - magnetic_energy(before)) / (after - before)
    }
}

fn resistances(p: BuckBoostParameters, m: BridgeModulation) -> [f64; 5] {
    let (a, b) = (m.input_high_fraction, m.output_high_fraction);
    [
        p.winding_resistance_ohm,
        a * p.input_on_resistance_ohm[0],
        (1.0 - a) * p.input_on_resistance_ohm[1],
        b * p.output_on_resistance_ohm[0],
        (1.0 - b) * p.output_on_resistance_ohm[1],
    ]
}

fn energy(x: [f64; 3], p: BuckBoostParameters) -> f64 {
    magnetic_energy(x[0])
        + 0.5 * x[1] * x[1] / p.bus_capacitance_f[0]
        + 0.5 * x[2] * x[2] / p.bus_capacitance_f[1]
}

fn exact_budget() -> LinearSplitBudget {
    LinearSplitBudget {
        max_modes: 3,
        max_omega_dt: 0.04,
    }
}

#[test]
fn saturation_and_regeneration_retain_reciprocal_bus_work_and_all_five_losses() {
    let (p, m, h) = (parameters(), modulation(0.6, 0.8), 0.2);
    let resistance = resistances(p, m);
    let total_r: f64 = resistance.iter().sum();
    for sign in [-1.0, 1.0] {
        let (flux0, flux1) = (sign * 0.5, sign * 1.5);
        let i_bar = chord_current(flux0, flux1);
        let input = [sign * 0.2, sign * -0.1];
        let dq = [
            h * (input[0] - m.input_high_fraction * i_bar),
            h * (input[1] + m.output_high_fraction * i_bar),
        ];
        let input_v0 = 12.0;
        let input_v_bar = input_v0 + dq[0] / (2.0 * p.bus_capacitance_f[0]);
        // Manufacture an endpoint using the independent capacitor balances
        // and original winding voltage law, with both bus voltages positive.
        let output_v_bar =
            (m.input_high_fraction * input_v_bar - (flux1 - flux0) / h - total_r * i_bar)
                / m.output_high_fraction;
        let voltage0 = [
            input_v0,
            output_v_bar - dq[1] / (2.0 * p.bus_capacitance_f[1]),
        ];
        assert!(voltage0.iter().all(|v| *v > 0.0));
        let model = nonlinear(p, m);
        let initial = model.state(current(flux0), voltage0).unwrap();
        let mut prepared = model
            .prepare_controlled(h, NonlinearConverterControls::default(), || false)
            .unwrap();
        let (mut next, mut voltage) = ([0.0; 3], [0.0; 2]);
        let record = prepared
            .step_controlled(&initial, &input, &mut next, &mut voltage, || false)
            .unwrap();
        close(next[0], flux1, 2e-10);
        close(record.discrete_current_a, i_bar, 2e-10);
        close(
            record.magnetic_energy_change_j,
            magnetic_energy(flux1) - magnetic_energy(flux0),
            2e-10,
        );
        // Crossing the knee makes this 1.5 A chord current distinct from the
        // 2 A arithmetic mean of endpoint currents, and from a frozen L=1 H.
        assert!((i_bar - 0.5 * (current(flux0) + current(flux1))).abs() > 0.4);
        for j in 0..2 {
            close(next[j + 1] - initial[j + 1], dq[j], 2e-10);
            let expected_voltage = [input_v_bar, output_v_bar][j];
            close(voltage[j], expected_voltage, 2e-10);
            close(
                record.bus_energy_change_j[j],
                dq[j] * expected_voltage,
                2e-10,
            );
            close(
                record.terminal_work_j[j],
                h * input[j] * expected_voltage,
                2e-10,
            );
        }
        for (heat, r) in record.conduction.joules.iter().zip(resistance) {
            assert!(*heat > 0.0);
            close(*heat, h * r * i_bar * i_bar, 2e-10);
        }
        close(
            record.native.delta_h,
            energy(next, p) - energy(initial, p),
            3e-10,
        );
        close(record.native.dissipated, h * total_r * i_bar * i_bar, 2e-10);
        close(record.native.balance_residual(), 0.0, 3e-10);
        assert_eq!(record.terminal_work_j[0] > 0.0, sign > 0.0);
        assert_eq!(record.terminal_work_j[1] > 0.0, sign < 0.0);
        if sign < 0.0 {
            assert!(next[1] > initial[1] && next[2] < initial[2]);
        } else {
            let linear = AveragedBuckBoost::new(p, m).unwrap();
            let mut linear = linear
                .prepare_controlled(h, exact_budget(), || false)
                .unwrap();
            let (mut linear_next, mut linear_voltage) = ([0.0; 3], [0.0; 2]);
            linear
                .step_controlled(
                    &initial,
                    &input,
                    &mut linear_next,
                    &mut linear_voltage,
                    || false,
                )
                .unwrap();
            assert!(current(next[0]) > linear_next[0] / p.inductance_h + 1.5);
        }
    }
}

#[test]
fn closed_finite_buses_exchange_saturated_magnetic_energy_without_hidden_supply() {
    let (p, m, h) = (parameters(), modulation(0.7, 0.6), 0.01);
    let model = nonlinear(p, m);
    let mut state = model.state(0.0, [6.0, 0.0]).unwrap();
    let initial_energy = energy(state, p);
    let conserved_charge = m.output_high_fraction * state[1] + m.input_high_fraction * state[2];
    // Rounding the first accepted input charge near 12 C contributes about
    // 5e-15 J: the zero-absolute default correctly refuses a finer requested
    // per-bus resolution. It must not publish either candidate output.
    let mut strict = nonlinear(p, m)
        .prepare_controlled(h, NonlinearConverterControls::default(), || false)
        .unwrap();
    let sentinel = ([7.0, 8.0, 9.0], [10.0, 11.0]);
    let (mut next, mut voltage) = sentinel;
    assert!(matches!(
        strict.step_controlled(&state, &[0.0; 2], &mut next, &mut voltage, || false),
        Err(PreparedStepError::Solver(PhsError::Dimension { what }))
            if what.contains("accounts do not close")
    ));
    assert_eq!((next, voltage), sentinel);
    // Declare two initial charge quanta at the input voltage, rather than
    // scaling an allowance by unchanged stored energy or relaxing the gate.
    let energy_quantum =
        2.0 * (state[1].next_up() - state[1]) * (state[1] / p.bus_capacitance_f[0]);
    let mut prepared = model
        .prepare_controlled(
            h,
            NonlinearConverterControls {
                absolute_energy_tolerance_j: energy_quantum,
                ..NonlinearConverterControls::default()
            },
            || false,
        )
        .unwrap();
    let resistance = resistances(p, m);
    let mut heat = [0.0; 5];
    let mut crossed_knee = false;
    for iteration in 0..160 {
        let (mut next, mut voltage) = ([0.0; 3], [0.0; 2]);
        let record = prepared
            .step_controlled(&state, &[0.0; 2], &mut next, &mut voltage, || false)
            .unwrap_or_else(|error| {
                let mut probe = StepWorkspace::new(prepared.system()).unwrap();
                let raw = probe.step_into_controlled(
                    prepared.system(),
                    &state,
                    &[0.0; 2],
                    h,
                    &mut next,
                    &mut voltage,
                    || false,
                );
                panic!(
                    "iteration={iteration}, initial={state:?}, candidate={next:?}, \
                     voltages={voltage:?}, raw={raw:?}, refusal={error:?}"
                );
            });
        let i_bar = chord_current(state[0], next[0]);
        close(
            next[1] - state[1],
            -h * m.input_high_fraction * i_bar,
            3e-10,
        );
        close(
            next[2] - state[2],
            h * m.output_high_fraction * i_bar,
            3e-10,
        );
        for j in 0..5 {
            let expected_heat = h * resistance[j] * i_bar * i_bar;
            close(record.conduction.joules[j], expected_heat, 3e-10);
            heat[j] += record.conduction.joules[j];
        }
        assert_eq!(record.terminal_work_j, [0.0; 2]);
        assert_eq!(record.native.supplied, 0.0);
        close(
            energy(next, p) + heat.iter().sum::<f64>(),
            initial_energy,
            2e-8,
        );
        close(
            m.output_high_fraction * next[1] + m.input_high_fraction * next[2],
            conserved_charge,
            2e-9,
        );
        crossed_knee |= next[0] > 1.0;
        state = next;
    }
    assert!(crossed_knee);
    assert!(heat.iter().all(|value| *value > 0.0));
}

#[test]
fn linear_flux_curve_refines_to_the_existing_exact_converter_trajectory() {
    let p = BuckBoostParameters {
        inductance_h: 0.02,
        bus_capacitance_f: [0.2, 0.3],
        ..parameters()
    };
    let (m, duration, input) = (modulation(0.63, 0.48), 0.05, [0.3, -0.4]);
    let linear = AveragedBuckBoost::new(p, m).unwrap();
    let initial = linear.state(1.7, [24.0, 8.0]).unwrap();
    let mut exact = linear
        .prepare_controlled(duration, exact_budget(), || false)
        .unwrap();
    let (mut expected, mut expected_voltage) = ([0.0; 3], [0.0; 2]);
    exact
        .step_controlled(
            &initial,
            &input,
            &mut expected,
            &mut expected_voltage,
            || false,
        )
        .unwrap();
    let errors = [8, 16, 32].map(|steps| {
        let source =
            InductorFluxCurve::try_new(vec![(-2.0, -100.0), (0.0, 0.0), (2.0, 100.0)]).unwrap();
        let model = AveragedBuckBoost::new(p, m)
            .unwrap()
            .with_inductor_curve(source)
            .unwrap();
        let mut prepared = model
            .prepare_controlled(
                duration / steps as f64,
                NonlinearConverterControls::default(),
                || false,
            )
            .unwrap();
        let mut state = initial;
        for _ in 0..steps {
            let (mut next, mut voltage) = ([0.0; 3], [0.0; 2]);
            prepared
                .step_controlled(&state, &input, &mut next, &mut voltage, || false)
                .unwrap();
            state = next;
        }
        let scales = [
            p.inductance_h,
            p.bus_capacitance_f[0],
            p.bus_capacitance_f[1],
        ];
        (0..3)
            .map(|j| ((state[j] - expected[j]) / scales[j]).abs())
            .fold(0.0_f64, f64::max)
    });
    assert!(errors[2] < 0.005, "errors={errors:?}");
    for pair in errors.windows(2) {
        assert!(
            (3.8..4.2).contains(&(pair[0] / pair[1])),
            "errors={errors:?}"
        );
    }
}

#[test]
fn into_system_retains_saturation_when_physical_bus_loads_are_composed() {
    let (p, m, h) = (parameters(), modulation(0.7, 0.6), 0.05);
    let model = nonlinear(p, m);
    let initial = model.state(2.0, [6.0, 3.0]).unwrap();
    let load = [0.1, 0.2];
    let system = model
        .into_system()
        .with_port_conductance(&[0, 1], &[load[0], 0.0, 0.0, load[1]])
        .unwrap();
    assert!(system.storage().exact_quadratic_matrix().is_none());
    close(system.effort(&initial)[0], 2.0, 1e-14);
    close(system.hamiltonian(&initial), energy(initial, p), 1e-13);
    let mut workspace = StepWorkspace::new(&system).unwrap();
    let (mut next, mut voltage) = ([0.0; 3], [0.0; 2]);
    let record = workspace
        .step_into_controlled(
            &system,
            &initial,
            &[0.0; 2],
            h,
            &mut next,
            &mut voltage,
            || false,
        )
        .unwrap();
    let i_bar = chord_current(initial[0], next[0]);
    let total_r: f64 = resistances(p, m).iter().sum();
    close(
        next[0] - initial[0],
        h * (m.input_high_fraction * voltage[0]
            - m.output_high_fraction * voltage[1]
            - total_r * i_bar),
        3e-10,
    );
    close(
        next[1] - initial[1],
        h * (-m.input_high_fraction * i_bar - load[0] * voltage[0]),
        3e-10,
    );
    close(
        next[2] - initial[2],
        h * (m.output_high_fraction * i_bar - load[1] * voltage[1]),
        3e-10,
    );
    let heat = h
        * (total_r * i_bar * i_bar
            + load[0] * voltage[0] * voltage[0]
            + load[1] * voltage[1] * voltage[1]);
    close(record.dissipated, heat, 3e-10);
    close(energy(next, p) - energy(initial, p) + heat, 0.0, 3e-10);
    assert_eq!(record.supplied, 0.0);
    assert!(current(next[0]) > next[0] / p.inductance_h + 0.5);
}

#[test]
fn source_budget_and_late_cancellation_refusals_publish_neither_output() {
    let (p, m, h) = (parameters(), modulation(0.6, 0.8), 0.1);
    let model = || nonlinear(p, m);
    assert!(model().state(17.0, [2.0, 1.0]).is_err());
    let initial = model().state(0.5, [2.0, 1.0]).unwrap();
    let input = [0.2, -0.1];
    let sentinel = ([7.0, 8.0, 9.0], [10.0, 11.0]);
    let (mut next, mut voltage) = sentinel;
    let mut limited = model()
        .prepare_controlled(
            h,
            NonlinearConverterControls {
                solve: SolveBudget {
                    max_iterations: 0,
                    ..SolveBudget::default()
                },
                ..NonlinearConverterControls::default()
            },
            || false,
        )
        .unwrap();
    assert!(
        limited
            .step_controlled(&initial, &input, &mut next, &mut voltage, || false)
            .is_err()
    );
    assert_eq!((next, voltage), sentinel);
    let mut prepared = model()
        .prepare_controlled(h, NonlinearConverterControls::default(), || false)
        .unwrap();
    let unsupported = model().state(15.5, [1000.0, 0.0]).unwrap();
    let held = unsupported;
    assert!(
        prepared
            .step_controlled(&unsupported, &[0.0; 2], &mut next, &mut voltage, || false)
            .is_err()
    );
    assert_eq!((next, voltage), sentinel);
    assert_eq!(unsupported, held);
    let mut polls = 0;
    let expected_record = prepared
        .step_controlled(&initial, &input, &mut next, &mut voltage, || {
            polls += 1;
            false
        })
        .unwrap();
    let expected = (next, voltage);
    assert!(polls > 1);
    (next, voltage) = sentinel;
    let mut count = 0;
    let cancelled = prepared.step_controlled(&initial, &input, &mut next, &mut voltage, || {
        count += 1;
        count == polls
    });
    assert!(matches!(cancelled, Err(PreparedStepError::Cancelled)));
    assert_eq!((next, voltage), sentinel);
    let retry = prepared
        .step_controlled(&initial, &input, &mut next, &mut voltage, || false)
        .unwrap();
    assert_eq!((next, voltage), expected);
    assert_eq!(retry.native, expected_record.native);
    assert_eq!(
        retry.magnetic_energy_change_j,
        expected_record.magnetic_energy_change_j
    );
    assert_eq!(
        retry.bus_energy_change_j,
        expected_record.bus_energy_change_j
    );
    assert_eq!(retry.conduction, expected_record.conduction);
}
