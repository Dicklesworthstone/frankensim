//! Regulated demand and nonlinear magnetic storage share one physical step.

use super::{
    NonlinearConverterControls, chord_current, close, current, energy, magnetic_energy, modulation,
    nonlinear, parameters, resistances,
};
use fs_phs::circuit::converter::power_load::{PowerLoadMode, RegulatedPowerLoad};
use fs_phs::{PreparedStepError, SolveBudget};

fn load(power: f64, maximum_current: f64, foldback: f64) -> RegulatedPowerLoad {
    RegulatedPowerLoad::try_new(power, maximum_current, foldback, 100.0).unwrap()
}

fn controls() -> NonlinearConverterControls {
    NonlinearConverterControls {
        solve: SolveBudget {
            relative: 1e-14,
            ..SolveBudget::default()
        },
        backtracking_halvings: 8,
        ..NonlinearConverterControls::default()
    }
}

#[test]
fn regulated_load_export_is_separate_from_saturated_storage_and_five_heat_owners() {
    let (p, m, h) = (parameters(), modulation(0.6, 0.8), 0.2);
    let (flux0, flux1) = (0.5, 1.5);
    let i_bar = chord_current(flux0, flux1);
    let resistance = resistances(p, m);
    let total_r: f64 = resistance.iter().sum();
    let dq = [-0.1, 0.15];
    let input_v0 = 12.0;
    let vin = input_v0 + dq[0] / (2.0 * p.bus_capacitance_f[0]);
    let vout = (m.input_high_fraction * vin - (flux1 - flux0) / h - total_r * i_bar)
        / m.output_high_fraction;
    let initial_voltage = [input_v0, vout - dq[1] / (2.0 * p.bus_capacitance_f[1])];
    let requested_power = [2.0, 3.0];
    let loads = [
        Some(load(requested_power[0], 10.0, 1.0)),
        Some(load(requested_power[1], 10.0, 1.0)),
    ];
    // Prescribe source currents using the wanted capacitor increments and the
    // load evaluated at the independently manufactured discrete voltage.
    let input = [
        dq[0] / h + m.input_high_fraction * i_bar + requested_power[0] / vin,
        dq[1] / h - m.output_high_fraction * i_bar + requested_power[1] / vout,
    ];
    let owner = nonlinear(p, m);
    let initial = owner.state(current(flux0), initial_voltage).unwrap();
    let mut prepared = owner.prepare_controlled(h, controls(), || false).unwrap();
    let (mut next, mut voltage) = ([0.0; 3], [0.0; 2]);
    let record = prepared
        .step_power_loaded_controlled(&initial, &input, &loads, &mut next, &mut voltage, || false)
        .unwrap();
    close(next[0], flux1, 3e-10);
    close(record.converter.discrete_current_a, 1.5, 3e-10);
    close(record.converter.magnetic_energy_change_j, 1.5, 3e-10);
    assert!(current(next[0]) > next[0] / p.inductance_h + 1.5);
    for j in 0..2 {
        let expected_voltage = [vin, vout][j];
        close(next[j + 1] - initial[j + 1], dq[j], 3e-10);
        close(voltage[j], expected_voltage, 3e-10);
        assert_eq!(record.loads[j].mode, PowerLoadMode::ConstantPower);
        assert!(record.loads[j].incremental_conductance_s < 0.0);
        close(
            record.loads[j].current_a,
            requested_power[j] / expected_voltage,
            3e-10,
        );
        close(
            record.source_work_j[j],
            h * input[j] * expected_voltage,
            3e-10,
        );
        close(record.load_work_j[j], h * requested_power[j], 3e-10);
        close(
            record.converter.bus_energy_change_j[j],
            dq[j] * expected_voltage,
            3e-10,
        );
        close(
            record.converter.terminal_work_j[j],
            record.source_work_j[j] - record.load_work_j[j],
            3e-10,
        );
        close(record.curtailed_load_energy_j[j], 0.0, 3e-10);
    }
    // A previous-voltage command is detectably different from the solved load.
    assert!((record.loads[1].current_a - requested_power[1] / initial_voltage[1]).abs() > 1e-3);
    for (heat, r) in record.converter.conduction.joules.iter().zip(resistance) {
        assert!(*heat > 0.0);
        close(*heat, h * r * i_bar * i_bar, 3e-10);
    }
    let heat = h * total_r * i_bar * i_bar;
    let exported: f64 = record.load_work_j.iter().sum();
    let supplied: f64 = record.source_work_j.iter().sum();
    close(record.converter.native.dissipated, heat, 3e-10);
    close(record.converter.native.supplied, supplied - exported, 3e-10);
    close(
        energy(next, p) - energy(initial, p) + heat + exported,
        supplied,
        5e-10,
    );
    assert!(exported > 10.0 * heat);
}

#[test]
fn bipolar_foldback_and_current_limits_discharge_real_buses_beside_saturated_winding_loss() {
    let (p, m) = (parameters(), modulation(0.0, 0.0));
    let resistance = resistances(p, m);
    let total_r: f64 = resistance.iter().sum();
    for sign in [-1.0, 1.0] {
        for limited in [false, true] {
            let (h, voltage0, demand, mode) = if limited {
                (
                    0.1,
                    10.0 * sign,
                    load(100.0, 2.0, 1.0),
                    PowerLoadMode::CurrentLimited,
                )
            } else {
                (
                    0.2,
                    sign,
                    load(12.0, 100.0, 6.0),
                    PowerLoadMode::VoltageFoldback,
                )
            };
            let owner = nonlinear(p, m);
            let initial = owner.state(sign * 2.0, [0.0, voltage0]).unwrap();
            let mut prepared = owner.prepare_controlled(h, controls(), || false).unwrap();
            let (mut next, mut voltage) = ([0.0; 3], [0.0; 2]);
            let record = prepared
                .step_power_loaded_controlled(
                    &initial,
                    &[0.0; 2],
                    &[None, Some(demand)],
                    &mut next,
                    &mut voltage,
                    || false,
                )
                .unwrap();
            let capacitance = p.bus_capacitance_f[1];
            let expected_voltage = if limited {
                voltage0 - h * 2.0 * sign / capacitance
            } else {
                let conductance = 12.0 / 36.0;
                voltage0 * (capacitance - 0.5 * h * conductance)
                    / (capacitance + 0.5 * h * conductance)
            };
            // The winding remains beyond its knee: I=5lambda-4sign. Its
            // independent affine midpoint RL solution is simultaneous with
            // the independent capacitive demand response.
            let expected_flux = (initial[0] * (1.0 - 2.5 * h * total_r) + 4.0 * sign * h * total_r)
                / (1.0 + 2.5 * h * total_r);
            close(next[0], expected_flux, 3e-10);
            assert!(next[0].abs() > 1.0);
            close(next[2] / capacitance, expected_voltage, 3e-10);
            close(voltage[1], 0.5 * (voltage0 + expected_voltage), 3e-10);
            assert_eq!(record.loads[1].mode, mode);
            assert_eq!(record.loads[1].current_a.is_sign_positive(), sign > 0.0);
            let load_current = if limited {
                2.0 * sign
            } else {
                voltage[1] / 3.0
            };
            close(record.loads[1].current_a, load_current, 3e-10);
            let exported =
                0.5 * capacitance * (voltage0 - expected_voltage) * (voltage0 + expected_voltage);
            let winding_heat = magnetic_energy(initial[0]) - magnetic_energy(next[0]);
            let i_bar = chord_current(initial[0], next[0]);
            close(record.load_work_j[1], exported, 3e-10);
            close(exported, h * voltage[1] * load_current, 3e-10);
            assert!(exported > 0.0 && winding_heat > 0.0);
            close(record.converter.native.dissipated, winding_heat, 3e-10);
            for (heat, r) in record.converter.conduction.joules.iter().zip(resistance) {
                close(*heat, h * r * i_bar * i_bar, 3e-10);
            }
            close(
                record.curtailed_load_energy_j[1],
                h * demand.requested_power_w() - exported,
                3e-10,
            );
            close(
                record.converter.native.delta_h,
                -winding_heat - exported,
                3e-10,
            );
            close(record.converter.native.supplied, -exported, 3e-10);
            assert_eq!(record.source_work_j, [0.0; 2]);
        }
    }
}

#[test]
fn absent_and_zero_demand_preserve_the_ordinary_nonlinear_converter_step() {
    let (p, m, h) = (parameters(), modulation(0.7, 0.6), 0.05);
    let owner = nonlinear(p, m);
    let initial = owner.state(2.0, [6.0, 3.0]).unwrap();
    let input = [0.2, -0.1];
    let mut ordinary = owner.prepare_controlled(h, controls(), || false).unwrap();
    let (mut expected, mut expected_voltage) = ([0.0; 3], [0.0; 2]);
    let expected_record = ordinary
        .step_controlled(
            &initial,
            &input,
            &mut expected,
            &mut expected_voltage,
            || false,
        )
        .unwrap();
    for loads in [[None, None], [Some(load(0.0, 1.0, 1.0)); 2]] {
        let mut prepared = nonlinear(p, m)
            .prepare_controlled(h, controls(), || false)
            .unwrap();
        let (mut next, mut voltage) = ([0.0; 3], [0.0; 2]);
        let record = prepared
            .step_power_loaded_controlled(&initial, &input, &loads, &mut next, &mut voltage, || {
                false
            })
            .unwrap();
        for (actual, expected) in next.iter().zip(expected) {
            close(*actual, expected, 3e-10);
        }
        for j in 0..2 {
            close(voltage[j], expected_voltage[j], 3e-10);
            assert_eq!(record.loads[j].mode, PowerLoadMode::Disabled);
            close(
                record.source_work_j[j],
                expected_record.terminal_work_j[j],
                3e-10,
            );
            close(
                record.converter.bus_energy_change_j[j],
                expected_record.bus_energy_change_j[j],
                3e-10,
            );
        }
        assert_eq!(record.load_work_j, [0.0; 2]);
        assert_eq!(record.curtailed_load_energy_j, [0.0; 2]);
        close(
            record.converter.magnetic_energy_change_j,
            expected_record.magnetic_energy_change_j,
            3e-10,
        );
        close(
            record.converter.native.dissipated,
            expected_record.native.dissipated,
            3e-10,
        );
        for (actual, expected) in record
            .converter
            .conduction
            .joules
            .iter()
            .zip(expected_record.conduction.joules)
        {
            close(*actual, expected, 3e-10);
        }
    }
}

#[test]
fn final_load_domain_source_budget_and_late_cancellation_refuse_transactionally() {
    let p = parameters();
    let (m, h) = (modulation(0.6, 0.8), 0.1);
    let owner = nonlinear(p, m);
    let initial = owner.state(0.5, [6.0, 3.0]).unwrap();
    let input = [0.2, 0.1];
    let loads = [None, Some(load(2.0, 10.0, 1.0))];
    let mut prepared = owner.prepare_controlled(h, controls(), || false).unwrap();
    let sentinel = ([7.0, 8.0, 9.0], [10.0, 11.0]);
    let (mut next, mut voltage) = sentinel;
    let outside_source = [4.01, initial[1], initial[2]];
    assert!(
        prepared
            .step_power_loaded_controlled(
                &outside_source,
                &input,
                &loads,
                &mut next,
                &mut voltage,
                || false,
            )
            .is_err()
    );
    assert_eq!((next, voltage), sentinel);

    let mut limited_controls = controls();
    limited_controls.solve.max_iterations = 0;
    let mut limited = nonlinear(p, m)
        .prepare_controlled(h, limited_controls, || false)
        .unwrap();
    assert!(
        limited
            .step_power_loaded_controlled(&initial, &input, &loads, &mut next, &mut voltage, || {
                false
            },)
            .is_err()
    );
    assert_eq!((next, voltage), sentinel);

    // A disabled authored load still owns its voltage support. The initial
    // 3.9 V and discrete 4 V are supported; the requested final 4.1 V is not.
    let isolated = nonlinear(p, modulation(0.0, 0.0));
    let charged = isolated.state(2.0, [0.0, 3.9]).unwrap();
    let mut isolated = isolated
        .prepare_controlled(h, controls(), || false)
        .unwrap();
    let bounded_zero = [
        None,
        Some(RegulatedPowerLoad::try_new(0.0, 1.0, 1.0, 4.0).unwrap()),
    ];
    assert!(
        isolated
            .step_power_loaded_controlled(
                &charged,
                &[0.0, 6.0],
                &bounded_zero,
                &mut next,
                &mut voltage,
                || false,
            )
            .is_err()
    );
    assert_eq!((next, voltage), sentinel);

    let mut polls = 0;
    let expected_record = prepared
        .step_power_loaded_controlled(&initial, &input, &loads, &mut next, &mut voltage, || {
            polls += 1;
            false
        })
        .unwrap();
    let expected = (next, voltage);
    assert!(polls > 1);
    (next, voltage) = sentinel;
    let mut count = 0;
    let cancelled = prepared.step_power_loaded_controlled(
        &initial,
        &input,
        &loads,
        &mut next,
        &mut voltage,
        || {
            count += 1;
            count == polls
        },
    );
    assert!(matches!(cancelled, Err(PreparedStepError::Cancelled)));
    assert_eq!((next, voltage), sentinel);
    let retry = prepared
        .step_power_loaded_controlled(&initial, &input, &loads, &mut next, &mut voltage, || false)
        .unwrap();
    assert_eq!((next, voltage), expected);
    assert_eq!(retry.converter.native, expected_record.converter.native);
    assert_eq!(retry.load_work_j, expected_record.load_work_j);
    assert_eq!(retry.source_work_j, expected_record.source_work_j);
    assert_eq!(
        retry.converter.conduction,
        expected_record.converter.conduction
    );
}
