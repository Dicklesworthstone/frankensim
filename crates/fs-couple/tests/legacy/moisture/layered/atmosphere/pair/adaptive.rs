//! Four-store continuous-time oracle for a heterogeneous finite barrier.
use super::*;
use fs_couple::moisture::adaptive::MoistureTimeControls;

fn time_controls(maximum: f64, tolerance: f64) -> MoistureTimeControls {
    MoistureTimeControls {
        initial_step_s: maximum,
        minimum_step_s: 1e-7,
        maximum_step_s: maximum,
        maximum_evaluations: 5000,
        moisture_absolute_tolerance: tolerance,
        activity_absolute_tolerance: tolerance,
        reservoir_water_absolute_tolerance_kg: tolerance * 1e-3,
        relative_tolerance: 0.0,
        water_balance_tolerance_kg: 1e-10,
    }
}

/// Independent series for exp(t A) applied to [RH_left,a_left,a_right,RH_right].
/// A is the continuous four-inventory generator, not a native implicit update.
/// Here ||t A||∞ <= 8, so 64 terms leave far below double-precision tail error.
fn continuous(
    initial: [f64; 4],
    gas_capacity: [f64; 2],
    boundary: [f64; 2],
    join: f64,
    duration: f64,
) -> [f64; 4] {
    let rate = |state: [f64; 4]| {
        let left = boundary[0] * (state[0] - state[1]);
        let middle = join * (state[1] - state[2]);
        let right = boundary[1] * (state[2] - state[3]);
        [
            -left / gas_capacity[0],
            (left - middle) / CAPACITY[0],
            (middle - right) / CAPACITY[1],
            right / gas_capacity[1],
        ]
    };
    let mut result = initial;
    let mut term = initial;
    for order in 1..=64 {
        term = rate(term).map(|value| value * duration / order as f64);
        for (value, increment) in result.iter_mut().zip(term) {
            *value += increment;
        }
    }
    result
}

#[test]
fn g1_paired_layered_refinement_and_continuation_match_four_native_inventories() {
    with_cx(|cx| {
        let ports = [
            Some(MoistureAtmospherePort::First),
            Some(MoistureAtmospherePort::Second),
        ];
        let gas = [air::atmosphere(0.85, 0.3), air::atmosphere(0.1, 0.7)];
        let saved = gas.clone();
        let capacity = gas
            .each_ref()
            .map(|g| g.vapor_mass_kg().unwrap() / g.relative_humidity().unwrap());
        let initial = [0.85, 0.3, 0.65, 0.1];
        for finite_transfer in [false, true] {
            let mut old = stack([initial[1], initial[2]], cx);
            let (boundary, join) = if finite_transfer {
                old = old
                    .with_transfer(
                        [
                            MoistureActivityTransfer::ConductanceKgM2S(40.0),
                            MoistureActivityTransfer::ConductanceKgM2S(60.0),
                        ],
                        vec![MoistureActivityTransfer::ConductanceKgM2S(30.0)],
                        cx,
                    )
                    .unwrap();
                (
                    [
                        1.0 / (1.0 / 0.016 + 1.0 / 0.040),
                        1.0 / (1.0 / 0.048 + 1.0 / 0.060),
                    ],
                    1.0 / (1.0 / 0.012 + 1.0 / 0.030),
                )
            } else {
                ([0.016, 0.048], 0.012)
            };
            let exact = continuous(initial, capacity, boundary, join, 0.4);
            let run = |time| {
                gas[0]
                    .advance_paired_layered_adaptive_to(
                        &gas[1],
                        &old,
                        0.4,
                        &ports,
                        config(),
                        nonlinear(),
                        controls(),
                        time,
                        cx,
                    )
                    .unwrap()
            };
            let (loose_gas, loose, loose_report) = run(time_controls(0.4, 2e-3));
            let (fine_gas, fine, report) = run(time_controls(0.4, 2e-5));
            let error = |air: &[fs_couple::moisture::atmosphere::FiniteMoistureAtmosphere; 2],
                         field: &MoistureLayeredSlab| {
                [
                    air[0].relative_humidity().unwrap(),
                    field.activity()[0],
                    field.activity()[1],
                    air[1].relative_humidity().unwrap(),
                ]
                .into_iter()
                .zip(exact)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0_f64, f64::max)
            };
            assert!(error(&fine_gas, &fine) < 0.3 * error(&loose_gas, &loose));
            let summary = &report.transport.summary;
            assert!(summary.evaluations > loose_report.transport.summary.evaluations);
            assert!(summary.rejected_intervals > 0);
            assert!(summary.maximum_accepted_error_ratio <= 1.0);
            assert_eq!(summary.interval_s, [0.0, 0.4]);
            assert_eq!(summary.activity, None);
            assert_eq!(fine.time_s(), 0.4);
            assert_eq!(fine.definition_identity(), old.definition_identity());
            assert_eq!(report.transport.layers.len(), 2);
            assert_eq!(
                report.transport.layers[0].boundary_water_kg[1],
                -report.transport.layers[1].boundary_water_kg[0]
            );
            for i in 0..2 {
                let layer = &report.transport.layers[i];
                near(
                    layer.water_change_kg,
                    layer.boundary_water_kg.iter().sum(),
                    1e-10,
                );
                near(
                    layer.water_change_kg,
                    fine.layers()[i].slab().water_mass_kg()
                        - old.layers()[i].slab().water_mass_kg(),
                    1e-10,
                );
                assert!(layer.maximum_accepted_error_ratio <= 1.0);
                let chamber = &report.atmospheres[i];
                near(
                    chamber.water_to_solids_kg,
                    summary.boundary_water_kg[i],
                    0.0,
                );
                near(chamber.water_balance_residual_kg, 0.0, 1e-10);
                near(chamber.total_mass_residual_kg, 0.0, 1e-10);
                near(
                    air::carrier_mass(&fine_gas[i]),
                    air::carrier_mass(&gas[i]),
                    1e-12,
                );
                assert_eq!(
                    summary.last_substep.transport.sorption_sources[i],
                    old.layers()[i].sorption().identity()
                );
            }
            near(summary.water_balance_residual_kg, 0.0, 1e-10);
            near(summary.total_mass_residual_kg.unwrap(), 0.0, 1e-10);
            near(
                fine_gas[0].vapor_mass_kg().unwrap()
                    + fine_gas[1].vapor_mass_kg().unwrap()
                    + fine.water_mass_kg(),
                gas[0].vapor_mass_kg().unwrap()
                    + gas[1].vapor_mass_kg().unwrap()
                    + old.water_mass_kg(),
                1e-10,
            );

            let (continued_gas, continued, next_report) = fine_gas[0]
                .advance_paired_layered_adaptive_to(
                    &fine_gas[1],
                    &fine,
                    0.8,
                    &ports,
                    config(),
                    nonlinear(),
                    controls(),
                    time_controls(0.4, 2e-5),
                    cx,
                )
                .unwrap();
            let exact = continuous(initial, capacity, boundary, join, 0.8);
            for (actual, expected) in [
                continued_gas[0].relative_humidity().unwrap(),
                continued.activity()[0],
                continued.activity()[1],
                continued_gas[1].relative_humidity().unwrap(),
            ]
            .into_iter()
            .zip(exact)
            {
                near(actual, expected, 2e-3);
            }
            assert_eq!(next_report.transport.summary.interval_s, [0.4, 0.8]);
            assert_eq!(continued.time_s(), 0.8);
            assert_eq!(old.time_s(), 0.0);
            near(old.activity()[0], 0.3, 1e-14);
            near(old.activity()[1], 0.65, 1e-14);
        }
        // Actual authored sealing overrides an otherwise assigned second port.
        // The complete-window debit must follow the native effective mapping.
        let sealed = stack([0.3, 0.65], cx)
            .with_transfer(
                [
                    MoistureActivityTransfer::IdealContact,
                    MoistureActivityTransfer::Sealed,
                ],
                vec![MoistureActivityTransfer::IdealContact],
                cx,
            )
            .unwrap();
        let (air, _, report) = gas[0]
            .advance_paired_layered_adaptive_to(
                &gas[1],
                &sealed,
                0.4,
                &ports,
                config(),
                nonlinear(),
                controls(),
                time_controls(0.4, 2e-4),
                cx,
            )
            .unwrap();
        assert_eq!(air[1], gas[1]);
        assert_eq!(report.atmospheres[1].water_to_solids_kg, 0.0);
        assert_eq!(
            report.transport.summary.last_substep.face_atmospheres,
            [ports[0], None]
        );
        assert_eq!(gas, saved);
    });
}
