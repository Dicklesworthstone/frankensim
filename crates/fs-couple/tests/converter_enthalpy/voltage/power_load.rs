//! Independent source/load calorimetry and material feedback at both buses.
use super::*;
use fs_couple::circuit_enthalpy::converter::voltage::power_load::*;
use fs_phs::SolveBudget;
use fs_phs::circuit::converter::power_load::{
    PowerLoadControls, PowerLoadMode, RegulatedPowerLoad,
};

fn power_controls() -> PowerLoadControls {
    PowerLoadControls {
        solve: SolveBudget {
            relative: 1e-13,
            absolute: 1e-12,
            max_iterations: 40,
        },
        backtracking_halvings: 8,
        absolute_energy_tolerance_j: 1e-10,
        relative_energy_tolerance: 1e-12,
    }
}
fn all_r() -> [f64; 7] {
    [R0[0], R0[1], R0[2], R0[3], R0[4], 1.0, 1.0]
}
fn fixed_paths() -> [ConverterThermalChannel<'static>; 7] {
    core::array::from_fn(|i| ConverterThermalChannel {
        cell: ALL_IDS[i],
        resistance: ConverterResistance::Fixed,
    })
}
fn load(power: f64) -> RegulatedPowerLoad {
    RegulatedPowerLoad::try_new(power, 10.0, 1.0, 50.0).unwrap()
}

#[test]
fn two_loads_share_one_finite_body_and_keep_each_source_heat_separate() {
    let cp = [1.0; 7];
    let (cards, curves) = sources(all_r(), cp, 0.0);
    let context = point(300.0);
    let cells = bindings(&cards, &curves, &context);
    let channels = fixed_paths();
    let network = VoltageDrivenConverterThermalNetwork {
        cells: &cells,
        links: &[],
        channels: &channels,
    };
    let thermal = body(&curves, cp);
    let modulation = modulation(0.0, 0.0);
    let native = AveragedBuckBoost::new(parameters(), modulation).unwrap();
    let initial = native.state(0.0, [8.0, 9.0]).unwrap();
    let converter = network
        .bind(&thermal, &native, [1.0; 2], initial, || false)
        .unwrap();
    let demands = [Some(load(16.0)), Some(load(27.0))];
    let dt = 0.02;
    let run = |receivers| {
        network
            .step_power_loaded_candidate(
                &thermal,
                dt,
                &converter,
                modulation,
                &[10.0, 12.0],
                &demands,
                receivers,
                power_controls(),
                controls(),
                &[0.0; 7],
                || false,
            )
            .unwrap()
    };
    let exported = run([None; 2]);
    let partial = run([Some(ALL_IDS[0]), None]);
    let heated = run([Some(ALL_IDS[0]); 2]);
    let report = heated.report();
    // Two uncoupled steady buses: (10-8)/1=2A and (12-9)/1=3A.
    // Loads use16+27W; source resistors heat4+9W; ideal supplies deliver20+36W.
    for port in 0..2 {
        near(report.electrical.source_current_a[port], [2.0, 3.0][port]);
        near(
            report.electrical.source_work_j[port],
            [20.0, 36.0][port] * dt,
        );
        near(
            report.electrical.terminal_resistor_heat_j[port],
            [4.0, 9.0][port] * dt,
        );
        near(report.electrical.load_work_j[port], [16.0, 27.0][port] * dt);
        let route = report.source.load_heat[port].unwrap();
        assert_eq!(route.cell, ALL_IDS[0]);
        near(route.energy_j, [16.0, 27.0][port] * dt);
    }
    near(report.deposited_joule_j, 56.0 * dt);
    near(exported.report().deposited_joule_j, 13.0 * dt);
    near(partial.report().deposited_joule_j, 29.0 * dt);
    near(report.electrical.converter.native.delta_h, 0.0);
    near(report.total_energy_residual_j, 0.0);
    near(report.non_electrical_heat_j, 0.0);
    assert_eq!(heated.converter().state(), exported.converter().state());
    assert_eq!(exported.report().source.load_heat, [None; 2]);
    for (i, expected_w) in [43.0, 0.0, 0.0, 0.0, 0.0, 4.0, 9.0].into_iter().enumerate() {
        near(
            heated.thermal().temperature_k(ALL_IDS[i]).unwrap() - 300.0,
            expected_w * dt,
        );
    }
    near(heated.thermal().energy_j() - thermal.energy_j(), 56.0 * dt);

    let (thermal, converter, _) = heated.into_parts();
    let stopped = network
        .step_power_loaded_candidate(
            &thermal,
            2.0 * dt,
            &converter,
            modulation,
            &[10.0, 12.0],
            &[None; 2],
            [Some(ALL_IDS[0]); 2],
            power_controls(),
            controls(),
            &[0.0; 7],
            || false,
        )
        .unwrap();
    assert_eq!(stopped.report().source.load_heat[0].unwrap().energy_j, 0.0);
    assert_eq!(stopped.report().source.load_heat[1].unwrap().energy_j, 0.0);
    assert_eq!(
        stopped.thermal().temperature_k(ALL_IDS[0]).unwrap(),
        thermal.temperature_k(ALL_IDS[0]).unwrap()
    );
}

#[test]
fn routed_load_heat_changes_its_source_resistance_and_original_bus_trajectory() {
    let cp = [0.01; 7];
    let alpha = 0.0039;
    let (cards, curves) = sources(all_r(), cp, alpha);
    let context = point(300.0);
    let cells = bindings(&cards, &curves, &context);
    let key = key();
    let channels = paths(ALL_IDS, &key);
    let network = VoltageDrivenConverterThermalNetwork {
        cells: &cells,
        links: &[],
        channels: &channels,
    };
    let thermal = body(&curves, cp);
    let modulation = modulation(0.0, 0.0);
    let native = AveragedBuckBoost::new(parameters(), modulation).unwrap();
    let initial = native.state(0.0, [0.0, 9.0]).unwrap();
    let converter = network
        .bind(&thermal, &native, [1.0; 2], initial, || false)
        .unwrap();
    let demands = [None, Some(load(27.0))];
    let dt = 0.02;
    let run = |receiver| {
        network
            .step_power_loaded_candidate(
                &thermal,
                dt,
                &converter,
                modulation,
                &[0.0, 12.0],
                &demands,
                [None, receiver],
                power_controls(),
                controls(),
                &[0.0; 7],
                || false,
            )
            .unwrap()
    };
    let exported = run(None);
    let heated = run(Some(ALL_IDS[6]));
    let theta = heated.thermal().temperature_k(ALL_IDS[6]).unwrap() - 300.0;
    let r = 1.0 + alpha * theta;
    let v = f64::midpoint(initial[2], heated.converter().state()[2])
        / parameters().bus_capacitance_f[1];
    let current = (12.0 - v) / r;
    near(heated.converter().terminal_resistance_ohm()[1], r);
    near(heated.report().electrical.source_current_a[1], current);
    near(cp[6] * theta, dt * (current * current * r + 27.0));
    near(heated.report().electrical.load_work_j[1], dt * 27.0);
    near(heated.report().total_energy_residual_j, 0.0);
    assert!(theta > 40.0);
    assert!(heated.converter().state()[2] < exported.converter().state()[2] - 1e-4);
}

#[test]
fn regeneration_keeps_negative_source_work_and_positive_material_heat() {
    let cp = [1.0; 7];
    let (cards, curves) = sources(all_r(), cp, 0.0);
    let context = point(300.0);
    let cells = bindings(&cards, &curves, &context);
    let channels = fixed_paths();
    let network = VoltageDrivenConverterThermalNetwork {
        cells: &cells,
        links: &[],
        channels: &channels,
    };
    let thermal = body(&curves, cp);
    let modulation = modulation(0.0, 0.0);
    let native = AveragedBuckBoost::new(parameters(), modulation).unwrap();
    let before = [6.0, 4.0];
    let source = [2.0, 0.0];
    let initial = native.state(0.0, before).unwrap();
    let converter = network
        .bind(&thermal, &native, [1.0; 2], initial, || false)
        .unwrap();
    let dt = 0.02;
    let result = network
        .step_power_loaded_candidate(
            &thermal,
            dt,
            &converter,
            modulation,
            &source,
            &[None; 2],
            [None; 2],
            power_controls(),
            controls(),
            &[0.0; 7],
            || false,
        )
        .unwrap();
    let mut work = 0.0;
    for port in 0..2 {
        let c = parameters().bus_capacitance_f[port];
        // Independent midpoint RC equation, with physical series resistance1.
        let after =
            (before[port] + dt / c * (source[port] - 0.5 * before[port])) / (1.0 + 0.5 * dt / c);
        let current = source[port] - 0.5 * (before[port] + after);
        let heat = dt * current * current;
        let supplied = dt * source[port] * current;
        near(result.converter().state()[port + 1], c * after);
        near(
            result.report().electrical.terminal_resistor_heat_j[port],
            heat,
        );
        near(
            result.thermal().temperature_k(ALL_IDS[5 + port]).unwrap() - 300.0,
            heat,
        );
        near(result.report().electrical.source_work_j[port], supplied);
        work += supplied;
    }
    assert!(result.report().electrical.source_work_j[0] < 0.0);
    assert!(result.report().deposited_joule_j > 0.0);
    near(
        result.report().electrical.converter.native.delta_h + result.thermal().energy_j()
            - thermal.energy_j(),
        work,
    );
    near(result.report().total_energy_residual_j, 0.0);
}

#[test]
fn foldback_and_current_limiting_heat_actual_delivery_only() {
    let cp = [1.0; 7];
    let (cards, curves) = sources(all_r(), cp, 0.0);
    let context = point(300.0);
    let cells = bindings(&cards, &curves, &context);
    let channels = fixed_paths();
    let network = VoltageDrivenConverterThermalNetwork {
        cells: &cells,
        links: &[],
        channels: &channels,
    };
    let thermal = body(&curves, cp);
    let modulation = modulation(0.0, 0.0);
    let native = AveragedBuckBoost::new(parameters(), modulation).unwrap();
    let initial = native.state(0.0, [0.0, 4.0]).unwrap();
    let converter = network
        .bind(&thermal, &native, [1.0; 2], initial, || false)
        .unwrap();
    let dt = 0.02;
    for (demand, mode) in [
        (
            RegulatedPowerLoad::try_new(4.0, 10.0, 10.0, 50.0).unwrap(),
            PowerLoadMode::VoltageFoldback,
        ),
        (
            RegulatedPowerLoad::try_new(40.0, 0.1, 1.0, 50.0).unwrap(),
            PowerLoadMode::CurrentLimited,
        ),
    ] {
        let result = network
            .step_power_loaded_candidate(
                &thermal,
                dt,
                &converter,
                modulation,
                &[0.0, 6.0],
                &[None, Some(demand)],
                [None, Some(ALL_IDS[0])],
                power_controls(),
                controls(),
                &[0.0; 7],
                || false,
            )
            .unwrap();
        let v = f64::midpoint(initial[2], result.converter().state()[2])
            / parameters().bus_capacitance_f[1];
        let actual = match mode {
            PowerLoadMode::VoltageFoldback => 4.0 * v * v / 100.0,
            _ => 0.1 * v,
        };
        assert_eq!(result.report().electrical.loads[1].mode, mode);
        near(
            result.report().source.load_heat[1].unwrap().energy_j,
            dt * actual,
        );
        near(
            result.thermal().temperature_k(ALL_IDS[0]).unwrap() - 300.0,
            dt * actual,
        );
        near(
            result.report().electrical.curtailed_load_energy_j[1],
            dt * (demand.requested_power_w() - actual),
        );
        assert!(actual < 0.5 * demand.requested_power_w());
        near(result.report().total_energy_residual_j, 0.0);
    }
}

#[test]
fn unsupported_receiver_calorimetry_and_late_cancel_preserve_both_histories() {
    let cp = [1.0; 7];
    let (cards, curves) = sources(all_r(), cp, 0.0);
    let context = point(300.0);
    let cells = bindings(&cards, &curves, &context);
    let channels = fixed_paths();
    let network = VoltageDrivenConverterThermalNetwork {
        cells: &cells,
        links: &[],
        channels: &channels,
    };
    let thermal = body(&curves, cp);
    let native = AveragedBuckBoost::new(parameters(), modulation(0.0, 0.0)).unwrap();
    let initial = native.state(0.0, [8.0, 9.0]).unwrap();
    let converter = network
        .bind(&thermal, &native, [1.0; 2], initial, || false)
        .unwrap();
    let energy = thermal.energy_j();
    let run =
        |loads: &[Option<RegulatedPowerLoad>; 2], receivers, cancel: &mut dyn FnMut() -> bool| {
            network.step_power_loaded_candidate(
                &thermal,
                0.02,
                &converter,
                modulation(0.0, 0.0),
                &[10.0, 12.0],
                loads,
                receivers,
                power_controls(),
                controls(),
                &[0.0; 7],
                cancel,
            )
        };
    assert!(run(&[None; 2], [None, Some(CellId(9999))], &mut || false).is_err());
    let demands = [Some(load(16.0)), Some(load(27.0))];
    let calls = Cell::new(0usize);
    let expected = run(&demands, [Some(ALL_IDS[0]); 2], &mut || {
        calls.set(calls.get() + 1);
        false
    })
    .unwrap();
    let stop = calls.get() / 2;
    assert!(stop > 10);
    calls.set(0);
    assert!(
        run(&demands, [Some(ALL_IDS[0]); 2], &mut || {
            calls.set(calls.get() + 1);
            calls.get() >= stop
        })
        .is_err()
    );
    assert_eq!(thermal.energy_j(), energy);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(converter.state(), &initial);
    let retry = run(&demands, [Some(ALL_IDS[0]); 2], &mut || false).unwrap();
    assert_eq!(retry.converter().state(), expected.converter().state());
    assert_eq!(retry.thermal().energy_j(), expected.thermal().energy_j());

    let mut tiny_cp = [1.0; 7];
    tiny_cp[0] = 1e-5;
    let (cards, curves) = sources(all_r(), tiny_cp, 0.0);
    let cells = bindings(&cards, &curves, &context);
    let network = VoltageDrivenConverterThermalNetwork {
        cells: &cells,
        links: &[],
        channels: &channels,
    };
    let thermal = body(&curves, tiny_cp);
    let converter = network
        .bind(&thermal, &native, [1.0; 2], initial, || false)
        .unwrap();
    let energy = thermal.energy_j();
    // The steady shared load supplies0.86J, beyond this cell's original chart.
    assert!(
        network
            .step_power_loaded_candidate(
                &thermal,
                0.02,
                &converter,
                modulation(0.0, 0.0),
                &[10.0, 12.0],
                &demands,
                [Some(ALL_IDS[0]); 2],
                power_controls(),
                controls(),
                &[0.0; 7],
                || false,
            )
            .is_err()
    );
    assert_eq!(thermal.energy_j(), energy);
    assert_eq!(converter.state(), &initial);
}
