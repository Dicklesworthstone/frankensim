//! Finite chemical energy through a real DC/DC stage and physical demand.
use fs_phs::circuit::battery::{BatteryOcvCurve, BatteryParameters, BatterySource};
use fs_phs::circuit::converter::battery_fed::BatteryFedBuckBoost;
use fs_phs::circuit::converter::power_load::{
    PowerLoadControls, PowerLoadMode, RegulatedPowerLoad,
};
use fs_phs::circuit::converter::{AveragedBuckBoost, BridgeModulation, BuckBoostParameters};
use fs_phs::{PreparedStepError, SolveBudget};

const LEAD_R: f64 = 0.04;
fn battery_parameters() -> BatteryParameters {
    BatteryParameters {
        series_resistance_ohm: 0.08,
        polarization_resistance_ohm: 0.1,
        polarization_capacitance_f: 20.0,
        bus_capacitance_f: 0.3,
    }
}
fn converter_parameters() -> BuckBoostParameters {
    BuckBoostParameters {
        inductance_h: 0.03,
        bus_capacitance_f: [0.2, 0.25],
        winding_resistance_ohm: 0.03,
        input_on_resistance_ohm: [0.01, 0.02],
        output_on_resistance_ohm: [0.01, 0.02],
    }
}
fn model(capacity: f64, curve: &[(f64, f64)]) -> BatteryFedBuckBoost {
    let source = BatterySource::new(
        BatteryOcvCurve::new(capacity, curve).unwrap(),
        battery_parameters(),
    )
    .unwrap();
    AveragedBuckBoost::new(
        converter_parameters(),
        BridgeModulation {
            input_high_fraction: 0.5,
            output_high_fraction: 1.0,
        },
    )
    .unwrap()
    .with_battery_source(source, LEAD_R)
    .unwrap()
}
fn controls() -> PowerLoadControls {
    PowerLoadControls {
        solve: SolveBudget {
            relative: 1e-13,
            absolute: 1e-13,
            max_iterations: 40,
        },
        backtracking_halvings: 8,
        absolute_energy_tolerance_j: 2e-10,
        relative_energy_tolerance: 2e-11,
    }
}
fn load(power: f64, limit: f64, foldback: f64) -> RegulatedPowerLoad {
    RegulatedPowerLoad::try_new(power, limit, foldback, 100.0).unwrap()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn plateau_battery_depletes_at_the_analytic_converter_operating_point() {
    let model = model(100.0, &[(0.0, 12.0), (1.0, 12.0)]);
    assert_eq!(model.system().state_dim(), 6);
    let bp = battery_parameters();
    let cp = converter_parameters();
    let converter_current = 2.0;
    let cell_current = 0.5 * converter_current;
    let battery_v =
        12.0 - (bp.series_resistance_ohm + bp.polarization_resistance_ohm) * cell_current;
    let input_v = battery_v - LEAD_R * cell_current;
    let converter_r = cp.winding_resistance_ohm
        + 0.5 * cp.input_on_resistance_ohm[0]
        + 0.5 * cp.input_on_resistance_ohm[1]
        + cp.output_on_resistance_ohm[0];
    let output_v = 0.5 * input_v - converter_r * converter_current;
    let demand_w = output_v * converter_current;
    let initial = model
        .state(
            0.6,
            bp.polarization_resistance_ohm * cell_current,
            battery_v,
            converter_current,
            [input_v, output_v],
        )
        .unwrap();
    let initial_energy = model.system().hamiltonian(&initial);
    let dt = 0.01;
    let mut prepared = model.prepare_controlled(dt, controls(), || false).unwrap();
    let expected_heat_w = [
        cp.winding_resistance_ohm * 4.0,
        0.5 * cp.input_on_resistance_ohm[0] * 4.0,
        0.5 * cp.input_on_resistance_ohm[1] * 4.0,
        cp.output_on_resistance_ohm[0] * 4.0,
        0.0,
        bp.series_resistance_ohm,
        bp.polarization_resistance_ohm,
        LEAD_R,
    ];
    close(
        demand_w + expected_heat_w.iter().sum::<f64>(),
        12.0 * cell_current,
        2e-14,
    );
    let mut state = initial;
    let mut exported = 0.0;
    let mut heat = 0.0;
    for _ in 0..100 {
        let mut next = [0.0; 6];
        let step = prepared
            .step_controlled(
                &state,
                0.0,
                Some(load(demand_w, 3.0, 3.0)),
                &mut next,
                || false,
            )
            .unwrap();
        assert_eq!(step.load.mode, PowerLoadMode::ConstantPower);
        assert_eq!(step.external_work_j, 0.0);
        close(next[0] - state[0], -dt * cell_current, 2e-12);
        for i in 1..6 {
            close(next[i], initial[i], 2e-11);
        }
        close(
            step.battery.chemical_energy_change_j,
            -12.0 * cell_current * dt,
            2e-11,
        );
        close(step.load_work_j, demand_w * dt, 2e-12);
        for (actual, expected) in step.heat_joules().into_iter().zip(expected_heat_w) {
            close(actual, dt * expected, 2e-12);
        }
        close(step.native.supplied, -step.load_work_j, 2e-12);
        close(step.balance_residual_j(), 0.0, 2e-11);
        exported += step.load_work_j;
        heat += step.heat_joules().iter().sum::<f64>();
        state = next;
    }
    close(state[0], initial[0] - 1.0, 2e-11);
    close(
        prepared.system().hamiltonian(&state) - initial_energy + exported + heat,
        0.0,
        2e-9,
    );
}

#[test]
fn changing_power_current_limit_and_foldback_share_the_same_finite_energy_states() {
    let model = model(100.0, &[(0.0, 10.0), (0.5, 12.0), (1.0, 13.0)]);
    let initial = model.state(0.75, 0.0, 12.5, 0.0, [12.5, 6.25]).unwrap();
    let curve = model.battery_curve().clone();
    let initial_energy = model.system().hamiltonian(&initial);
    let bp = battery_parameters();
    let cp = converter_parameters();
    let dt = 0.002;
    let mut prepared = model.prepare_controlled(dt, controls(), || false).unwrap();
    let mut state = initial;
    let mut exported = 0.0;
    let mut total_heat = 0.0;
    for (demand, mode) in [
        (load(2.0, 2.0, 3.0), PowerLoadMode::ConstantPower),
        (load(20.0, 1.0, 3.0), PowerLoadMode::CurrentLimited),
        (load(8.0, 2.0, 8.0), PowerLoadMode::VoltageFoldback),
    ] {
        for _ in 0..200 {
            let mut next = [0.0; 6];
            let step = prepared
                .step_controlled(&state, 0.0, Some(demand), &mut next, || false)
                .unwrap();
            assert_eq!(step.load.mode, mode);
            let mid = core::array::from_fn::<_, 6, _>(|i| f64::midpoint(state[i], next[i]));
            let vb = mid[2] / bp.bus_capacitance_f;
            let vi = mid[4] / cp.bus_capacitance_f[0];
            let vo = mid[5] / cp.bus_capacitance_f[1];
            let lead_current = (vb - vi) / LEAD_R;
            let winding_current = mid[3] / cp.inductance_h;
            let chemical = curve.energy_change_j(state[0], next[0]).unwrap();
            let stored = chemical
                + (next[1] - state[1]) * mid[1] / bp.polarization_capacitance_f
                + (next[2] - state[2]) * vb
                + (next[3] - state[3]) * winding_current
                + (next[4] - state[4]) * vi
                + (next[5] - state[5]) * vo;
            let delivered = dt * vo * step.load.current_a;
            let heat = step.heat_joules().iter().sum::<f64>();
            close(stored + heat + delivered, 0.0, 3e-10);
            close(step.native.dissipated, heat, 3e-11);
            close(step.load_work_j, delivered, 3e-12);
            close(
                step.lead_heat_j,
                dt * LEAD_R * lead_current * lead_current,
                3e-12,
            );
            close(
                -step.battery.terminal_work_j - step.converter_terminal_work_j[0],
                step.lead_heat_j,
                3e-12,
            );
            close(
                (next[0] - state[0]) + (next[2] - state[2]),
                -dt * lead_current,
                3e-11,
            );
            assert!(step.heat_joules().iter().all(|value| *value >= 0.0));
            if mode != PowerLoadMode::ConstantPower {
                assert!(step.curtailed_load_energy_j > 0.0);
            }
            exported += delivered;
            total_heat += heat;
            state = next;
        }
    }
    assert!(state[0] < initial[0] - 0.1);
    assert!(exported > 3.0 && total_heat > 0.0);
    close(
        prepared.system().hamiltonian(&state) - initial_energy + exported + total_heat,
        0.0,
        3e-9,
    );
}

#[test]
fn output_generation_returns_energy_through_the_converter_into_chemical_storage() {
    let model = model(100.0, &[(0.0, 10.0), (0.5, 12.0), (1.0, 13.0)]);
    let initial = model.state(0.4, 0.0, 11.6, 0.0, [11.6, 5.8]).unwrap();
    let initial_energy = model.system().hamiltonian(&initial);
    let mut prepared = model
        .prepare_controlled(0.005, controls(), || false)
        .unwrap();
    let mut state = initial;
    let mut supplied = 0.0;
    let mut heat = 0.0;
    let mut charging_steps = 0;
    for _ in 0..1000 {
        let mut next = [0.0; 6];
        let step = prepared
            .step_controlled(&state, 1.0, None, &mut next, || false)
            .unwrap();
        assert_eq!(step.load_work_j, 0.0);
        close(step.balance_residual_j(), 0.0, 3e-10);
        if step.battery.discrete_cell_current_a < 0.0 && step.lead_current_a < 0.0 {
            assert!(step.battery.terminal_work_j > 0.0);
            assert!(step.converter_terminal_work_j[0] < 0.0);
            charging_steps += 1;
        }
        supplied += step.external_work_j;
        heat += step.heat_joules().iter().sum::<f64>();
        state = next;
    }
    assert!(charging_steps > 800);
    assert!(state[0] > initial[0] + 1.0);
    assert!(state[3] < 0.0);
    assert!(supplied > 20.0 && heat > 0.0);
    close(
        prepared.system().hamiltonian(&state) - initial_energy + heat,
        supplied,
        3e-9,
    );
}

#[test]
fn unsupported_demand_exhaustion_and_cancellation_leave_all_six_states_unpublished() {
    let assembly = model(100.0, &[(0.0, 12.0), (1.0, 12.0)]);
    let initial = assembly.state(0.5, 0.0, 12.0, 0.0, [12.0, 6.0]).unwrap();
    let mut prepared = assembly
        .prepare_controlled(0.01, controls(), || false)
        .unwrap();
    let mut next = [91.0; 6];
    let unsupported = RegulatedPowerLoad::try_new(1.0, 1.0, 0.5, 1.0).unwrap();
    assert!(
        prepared
            .step_controlled(&initial, 0.0, Some(unsupported), &mut next, || false)
            .is_err()
    );
    assert_eq!(next, [91.0; 6]);
    assert!(matches!(
        prepared.step_controlled(&initial, 0.0, Some(load(2.0, 2.0, 3.0)), &mut next, || true),
        Err(PreparedStepError::Cancelled)
    ));
    assert_eq!(next, [91.0; 6]);
    let first = prepared
        .step_controlled(&initial, 0.0, Some(load(2.0, 2.0, 3.0)), &mut next, || {
            false
        })
        .unwrap();
    let accepted = next;
    let retry = prepared
        .step_controlled(&initial, 0.0, Some(load(2.0, 2.0, 3.0)), &mut next, || {
            false
        })
        .unwrap();
    assert_eq!(next, accepted);
    assert_eq!(first.native, retry.native);
    assert_eq!(first.heat_joules(), retry.heat_joules());

    let small = model(0.02, &[(0.0, 12.0), (1.0, 12.0)]);
    let initial = small.state(0.5, 0.1, 11.82, 2.0, [11.78, 5.78]).unwrap();
    let mut prepared = small
        .prepare_controlled(0.02, controls(), || false)
        .unwrap();
    next = [91.0; 6];
    assert!(
        prepared
            .step_controlled(
                &initial,
                0.0,
                Some(load(11.56, 3.0, 3.0)),
                &mut next,
                || false
            )
            .is_err()
    );
    assert_eq!(next, [91.0; 6]);
}
