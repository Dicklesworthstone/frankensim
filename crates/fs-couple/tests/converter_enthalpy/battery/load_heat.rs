//! G1/G3: actual regulated-load work heats its finite receiver exactly once.
use super::*;

// Retain all original physical resistances, but share the lead/polarization
// receiver. The last finite material cell then represents only the load.
fn isolated_load_channels() -> [ConverterThermalChannel<'static>; 8] {
    core::array::from_fn(|i| ConverterThermalChannel {
        cell: ALL_IDS[if i == 7 { 6 } else { i }],
        resistance: ConverterResistance::Fixed,
    })
}

#[test]
fn finite_load_calorimetry_closes_chemical_depletion_without_double_counting() {
    let fixture = Fixture::new([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0], 0.0);
    let context = point(300.0);
    let cells = fixture.cells(&context);
    let channels = isolated_load_channels();
    let network = BatteryFedConverterThermalNetwork {
        cells: &cells,
        links: &[],
        channels: &channels,
    };
    let body = fixture.body();
    let (native, state, load, heat_w) = steady();
    let converter = network.bind(&body, &native, state, || false).unwrap();
    let dt = 0.02;
    let run = |receiver| {
        network
            .step_with_load_heat_candidate(
                &body,
                dt,
                &converter,
                native.modulation(),
                0.0,
                Some(load),
                receiver,
                electrical_controls(),
                controls(),
                &[0.0; 8],
                || false,
            )
            .unwrap()
    };
    let exported = run(None);
    let heated = run(Some(ALL_IDS[7]));
    let report = heated.report();
    // The independently derived steady load is 8.031 W and total resistive
    // loss is 0.369 W. With flat 12 V OCV, chemical depletion is 8.4 W.
    near(load.requested_power_w(), 8.031);
    near(report.electrical.battery.chemical_energy_change_j, -8.4 * dt);
    near(report.electrical.load_work_j, 8.031 * dt);
    near(report.deposited_joule_j, 8.4 * dt);
    near(report.non_electrical_heat_j, 0.0);
    near(report.total_energy_residual_j, 0.0);
    assert_eq!(exported.converter().state(), heated.converter().state());
    assert_eq!(exported.report().source.load_heat, None);
    assert_eq!(
        report.source.load_heat,
        Some(RegulatedLoadHeat {
            cell: ALL_IDS[7],
            energy_j: report.electrical.load_work_j,
        })
    );
    for i in 0..8 {
        let power = match i {
            6 => heat_w[6] + heat_w[7],
            7 => 8.031,
            _ => heat_w[i],
        };
        let measured = fixture.cp[i]
            * (heated.thermal().temperature_k(ALL_IDS[i]).unwrap() - 300.0);
        assert!((measured - power * dt).abs() < 2e-9, "cell {i}");
        near(report.electrical.heat_joules()[i], heat_w[i] * dt);
    }
    near(exported.thermal().temperature_k(ALL_IDS[7]).unwrap(), 300.0);
    near(
        heated.thermal().energy_j() - exported.thermal().energy_j(),
        report.electrical.load_work_j,
    );
    near(
        report.electrical.native.delta_h + heated.thermal().energy_j() - body.energy_j(),
        0.0,
    );

    // Publishing both histories and turning the command off must not deposit
    // the previous interval's load work again.
    let (body, converter, _) = heated.into_parts();
    let disabled = network
        .step_with_load_heat_candidate(
            &body,
            2.0 * dt,
            &converter,
            native.modulation(),
            0.0,
            None,
            Some(ALL_IDS[7]),
            electrical_controls(),
            controls(),
            &[0.0; 8],
            || false,
        )
        .unwrap();
    assert_eq!(disabled.report().source.load_heat.unwrap().energy_j, 0.0);
    assert_eq!(
        disabled.thermal().temperature_k(ALL_IDS[7]).unwrap(),
        body.temperature_k(ALL_IDS[7]).unwrap()
    );
}

#[test]
fn foldback_and_current_limit_heat_only_actual_delivered_power() {
    let fixture = Fixture::new([1.0; 8], 0.0);
    let context = point(300.0);
    let cells = fixture.cells(&context);
    let channels = isolated_load_channels();
    let network = BatteryFedConverterThermalNetwork {
        cells: &cells,
        links: &[],
        channels: &channels,
    };
    let body = fixture.body();
    let (native, before, _, _) = steady();
    let converter = network.bind(&body, &native, before, || false).unwrap();
    let dt = 0.02;
    let exported_current_a = -0.2;
    for (load, mode) in [
        (
            RegulatedPowerLoad::try_new(4.0, 10.0, 20.0, 50.0).unwrap(),
            PowerLoadMode::VoltageFoldback,
        ),
        (
            RegulatedPowerLoad::try_new(40.0, 0.1, 1.0, 50.0).unwrap(),
            PowerLoadMode::CurrentLimited,
        ),
    ] {
        let candidate = network
            .step_with_load_heat_candidate(
                &body,
                dt,
                &converter,
                native.modulation(),
                exported_current_a,
                Some(load),
                Some(ALL_IDS[7]),
                electrical_controls(),
                controls(),
                &[0.0; 8],
                || false,
            )
            .unwrap();
        let report = candidate.report();
        let voltage = f64::midpoint(before[5], candidate.converter().state()[5])
            / parameters().bus_capacitance_f[1];
        let actual_power = match mode {
            PowerLoadMode::VoltageFoldback => 4.0 * voltage * voltage / 400.0,
            _ => 0.1 * voltage,
        };
        assert_eq!(report.electrical.load.mode, mode);
        near(report.source.load_heat.unwrap().energy_j, dt * actual_power);
        near(
            candidate.thermal().temperature_k(ALL_IDS[7]).unwrap() - 300.0,
            dt * actual_power,
        );
        near(
            report.electrical.curtailed_load_energy_j,
            dt * (load.requested_power_w() - actual_power),
        );
        // An independent output-current load stays outside this receiver.
        // Only the explicitly routed regulated demand becomes internal heat.
        near(report.electrical.external_work_j, dt * voltage * exported_current_a);
        near(
            report.electrical.native.delta_h + candidate.thermal().energy_j() - body.energy_j(),
            dt * voltage * exported_current_a,
        );
        assert!(actual_power < 0.5 * load.requested_power_w());
        near(report.total_energy_residual_j, 0.0);
    }
}

#[test]
fn shared_load_heat_feeds_back_into_original_material_resistance_and_current() {
    let fixture = Fixture::new([0.01; 8], 0.0039);
    let context = point(300.0);
    let cells = fixture.cells(&context);
    let resistance = key();
    let channels = paths(&resistance);
    let network = BatteryFedConverterThermalNetwork {
        cells: &cells,
        links: &[],
        channels: &channels,
    };
    let body = fixture.body();
    let (native, before, load, _) = steady();
    let converter = network.bind(&body, &native, before, || false).unwrap();
    let run = |receiver| {
        network
            .step_with_load_heat_candidate(
                &body,
                0.02,
                &converter,
                native.modulation(),
                0.0,
                Some(load),
                receiver,
                electrical_controls(),
                controls(),
                &[0.0; 8],
                || false,
            )
            .unwrap()
    };
    let exported = run(None);
    let heated = run(Some(ALL_IDS[0]));
    let temperature = heated.thermal().temperature_k(ALL_IDS[0]).unwrap();
    let resistance = R0[0] * (1.0 + 0.0039 * (temperature - 300.0));
    let current = f64::midpoint(before[3], heated.converter().state()[3])
        / parameters().inductance_h;
    near(heated.converter().channels()[0].resistance_ohm(), resistance);
    near(
        0.01 * (temperature - 300.0),
        0.02 * resistance * current * current + heated.report().electrical.load_work_j,
    );
    assert!(
        temperature > exported.thermal().temperature_k(ALL_IDS[0]).unwrap() + 10.0
    );
    assert!((heated.converter().state()[3] - exported.converter().state()[3]).abs() > 1e-6);
    near(heated.report().total_energy_residual_j, 0.0);
}

#[test]
fn invalid_receivers_and_late_cancellation_preserve_both_histories() {
    let fixture = Fixture::new([1.0; 8], 0.0);
    let context = point(300.0);
    let cells = fixture.cells(&context);
    let channels = isolated_load_channels();
    let network = BatteryFedConverterThermalNetwork {
        cells: &cells,
        links: &[],
        channels: &channels,
    };
    let body = fixture.body();
    let energy = body.energy_j();
    let (native, before, load, _) = steady();
    let converter = network.bind(&body, &native, before, || false).unwrap();
    let run = |load, cell, cancel: &mut dyn FnMut() -> bool| {
        network.step_with_load_heat_candidate(
            &body,
            0.02,
            &converter,
            native.modulation(),
            0.0,
            load,
            Some(cell),
            electrical_controls(),
            controls(),
            &[0.0; 8],
            cancel,
        )
    };
    // Even a dormant command cannot authorize a nonexistent thermal receiver.
    assert!(run(None, CellId(9999), &mut || false).is_err());
    let calls = Cell::new(0usize);
    let expected = run(Some(load), ALL_IDS[7], &mut || {
        calls.set(calls.get() + 1);
        false
    })
    .unwrap();
    let stop = calls.get() / 2;
    assert!(stop > 10);
    calls.set(0);
    assert!(run(Some(load), ALL_IDS[7], &mut || {
        calls.set(calls.get() + 1);
        calls.get() >= stop
    })
    .is_err());
    assert_eq!(body.energy_j(), energy);
    assert_eq!(body.time_s(), 0.0);
    assert_eq!(converter.state(), &before);
    assert_eq!(converter.time_s(), 0.0);
    let retried = run(Some(load), ALL_IDS[7], &mut || false).unwrap();
    assert_eq!(retried.converter().state(), expected.converter().state());
    assert_eq!(retried.thermal().energy_j(), expected.thermal().energy_j());
}

#[test]
fn dissipative_load_cannot_heat_beyond_its_material_caloric_support() {
    let mut capacities = [1.0; 8];
    capacities[7] = 1e-5;
    let fixture = Fixture::new(capacities, 0.0);
    let context = point(300.0);
    let cells = fixture.cells(&context);
    let channels = isolated_load_channels();
    let network = BatteryFedConverterThermalNetwork {
        cells: &cells,
        links: &[],
        channels: &channels,
    };
    let body = fixture.body();
    let energy = body.energy_j();
    let (native, before, load, _) = steady();
    let converter = network.bind(&body, &native, before, || false).unwrap();
    // The fixed 0.16062 J load needs 16062 K of temperature rise, outside the
    // original [200,2000] K chart. No extrapolated endpoint may be published.
    assert!(
        network
            .step_with_load_heat_candidate(
                &body,
                0.02,
                &converter,
                native.modulation(),
                0.0,
                Some(load),
                Some(ALL_IDS[7]),
                electrical_controls(),
                controls(),
                &[0.0; 8],
                || false,
            )
            .is_err()
    );
    assert_eq!(body.energy_j(), energy);
    assert_eq!(body.time_s(), 0.0);
    assert_eq!(converter.state(), &before);
}
