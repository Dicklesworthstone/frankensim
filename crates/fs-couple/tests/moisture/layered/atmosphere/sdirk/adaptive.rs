//! Original duration control over coupled native finite-air stages.
use super::*;
use fs_couple::moisture::adaptive::{
    LayeredMoistureAdvanceReport, MoistureAdvanceError, MoistureTimeControls,
};
use fs_couple::moisture::atmosphere::AtmosphereMoistureStep;
use fs_couple::moisture::layered::MoistureLayeredSdirk2Report;

type Report = LayeredMoistureAdvanceReport<AtmosphereMoistureStep<MoistureLayeredSdirk2Report>>;

fn time(maximum: f64, tolerance: f64) -> MoistureTimeControls {
    MoistureTimeControls {
        initial_step_s: maximum,
        minimum_step_s: 1e-7,
        maximum_step_s: maximum,
        maximum_evaluations: 5000,
        moisture_absolute_tolerance: tolerance,
        activity_absolute_tolerance: tolerance,
        reservoir_water_absolute_tolerance_kg: 1.0,
        relative_tolerance: 0.0,
        water_balance_tolerance_kg: 1e-10,
    }
}

fn account(report: &Report, capacity: [f64; 2], initial: [f64; 2], next: &[f64], duration: f64) {
    assert_eq!(report.summary.interval_s, [0.0, duration]);
    assert_eq!(report.layers.len(), 2);
    assert_eq!(
        report.layers[0].boundary_water_kg[1],
        -report.layers[1].boundary_water_kg[0]
    );
    assert_eq!(
        report.summary.accepted_substeps,
        2 * report.summary.accepted_intervals
    );
    assert!(report.summary.maximum_accepted_error_ratio <= 1.0);
    assert_eq!(report.summary.evaluations % 2, 0);
    for i in 0..2 {
        let change = capacity[i] * (next[i] - initial[i]);
        near(report.layers[i].water_change_kg, change, 1e-10);
        near(
            report.layers[i].boundary_water_kg.iter().sum(),
            change,
            1e-10,
        );
        assert_eq!(
            report.summary.boundary_water_kg[i],
            report.layers[i].boundary_water_kg[i]
        );
    }
    near(report.summary.water_balance_residual_kg, 0.0, 1e-10);
    near(report.summary.transport_balance_residual_kg, 0.0, 1e-10);
    near(
        report.summary.maximum_carrier_mass_residual_kg.unwrap(),
        0.0,
        1e-10,
    );
    near(report.summary.total_mass_residual_kg.unwrap(), 0.0, 1e-10);
    near(
        report.summary.water_change_kg,
        report.summary.boundary_water_kg.iter().sum(),
        1e-10,
    );
}

#[test]
fn finite_air_duration_accuracy_refines_toward_the_continuous_closed_system() {
    with_cx(|cx| {
        for symmetry in [
            None,
            Some(RadialMoistureSymmetry::Sphere),
            Some(RadialMoistureSymmetry::Cylinder { length_m: 0.1 }),
        ] {
            let (initial, capacity, links, volume, duration) = if let Some(symmetry) = symmetry {
                let p = curved(symmetry);
                (
                    [0.1, 0.2, 0.8],
                    p.capacity,
                    vec![(0, 1, p.join), (1, 2, p.outer)],
                    0.5,
                    8.0,
                )
            } else {
                ([0.2, 0.1, 0.8], CAPACITY, SLAB_LINKS.to_vec(), 1.0, 1.0)
            };
            let gas = air::atmosphere(initial[2], volume);
            let saved = gas.clone();
            let full_capacity = [
                capacity[0],
                capacity[1],
                gas.vapor_mass_kg().unwrap() / initial[2],
            ];
            let rate = rate_matrix(full_capacity, &links);
            let exact = exact_three(full_capacity, &rate, initial, duration);
            let mut errors = Vec::new();
            let mut evaluations = Vec::new();
            for tolerance in [1e-3, 1e-6] {
                let (next_gas, activity, final_water, old_wet, next_wet, report) =
                    if let Some(symmetry) = symmetry {
                        let old = coated(symmetry, [initial[0], initial[1]], cx);
                        let (next_gas, next, report) = gas
                            .advance_layered_radial_sdirk2_adaptive_to(
                                &old,
                                duration,
                                [false, true],
                                config(),
                                nonlinear(),
                                air::controls(),
                                time(duration, tolerance),
                                cx,
                            )
                            .unwrap();
                        assert_eq!(old.time_s(), 0.0);
                        assert_eq!(next.time_s(), duration);
                        (
                            next_gas,
                            next.activity().to_vec(),
                            next.water_mass_kg(),
                            old.wet_mass_kg(),
                            next.wet_mass_kg(),
                            report,
                        )
                    } else {
                        let old = stack([initial[0], initial[1]], cx);
                        let (next_gas, next, report) = gas
                            .advance_layered_sdirk2_adaptive_to(
                                &old,
                                duration,
                                [true; 2],
                                config(),
                                nonlinear(),
                                air::controls(),
                                time(duration, tolerance),
                                cx,
                            )
                            .unwrap();
                        assert_eq!(old.time_s(), 0.0);
                        assert_eq!(next.time_s(), duration);
                        (
                            next_gas,
                            next.activity().to_vec(),
                            next.water_mass_kg(),
                            old.wet_mass_kg(),
                            next.wet_mass_kg(),
                            report,
                        )
                    };
                account(
                    &report,
                    capacity,
                    [initial[0], initial[1]],
                    &activity,
                    duration,
                );
                if symmetry.is_some() {
                    assert_eq!(report.summary.boundary_water_kg[0], 0.0);
                }
                let original_water = capacity[0] * initial[0] + capacity[1] * initial[1];
                near(
                    next_gas.vapor_mass_kg().unwrap() + final_water,
                    gas.vapor_mass_kg().unwrap() + original_water,
                    1e-10,
                );
                near(
                    next_gas.gas().mass_kg() + next_wet,
                    gas.gas().mass_kg() + old_wet,
                    1e-10,
                );
                near(air::carrier_mass(&next_gas), air::carrier_mass(&gas), 1e-10);
                assert_eq!(
                    report.summary.reservoir_water_kg.unwrap(),
                    [
                        gas.vapor_mass_kg().unwrap(),
                        next_gas.vapor_mass_kg().unwrap()
                    ]
                );
                let state = [
                    activity[0],
                    activity[1],
                    next_gas.relative_humidity().unwrap(),
                ];
                errors.push(
                    (0..3)
                        .map(|i| (state[i] - exact[i]).abs())
                        .fold(0.0_f64, f64::max),
                );
                evaluations.push(report.summary.evaluations);
            }
            assert!(errors[1] < 0.1 * errors[0], "{symmetry:?}: {errors:?}");
            assert!(
                evaluations[1] > evaluations[0],
                "{symmetry:?}: {evaluations:?}"
            );
            assert_eq!(gas, saved);
        }
    });
}

#[test]
fn loose_finite_air_accuracy_matches_native_halves_without_debiting_stage_diagnostics() {
    with_cx(|cx| {
        let initial = [0.2, 0.1];
        for symmetry in [
            None,
            Some(RadialMoistureSymmetry::Sphere),
            Some(RadialMoistureSymmetry::Cylinder { length_m: 0.1 }),
        ] {
            let duration = if symmetry.is_some() { 2.0 } else { 0.25 };
            let gas = air::atmosphere(0.8, 1.0);
            let (actual_gas, actual, report, expected_gas, expected, receipts, capacity) =
                if let Some(symmetry) = symmetry {
                    let old = coated(symmetry, initial, cx);
                    let (half_gas, half, first) = gas
                        .advance_layered_radial_sdirk2(
                            &old,
                            0.5 * duration,
                            [false, true],
                            config(),
                            nonlinear(),
                            air::controls(),
                            cx,
                        )
                        .unwrap()
                        .into_parts();
                    let (expected_gas, expected, second) = half_gas
                        .advance_layered_radial_sdirk2(
                            &half,
                            0.5 * duration,
                            [false, true],
                            config(),
                            nonlinear(),
                            air::controls(),
                            cx,
                        )
                        .unwrap()
                        .into_parts();
                    let (actual_gas, actual, report) = gas
                        .advance_layered_radial_sdirk2_adaptive_to(
                            &old,
                            duration,
                            [false, true],
                            config(),
                            nonlinear(),
                            air::controls(),
                            time(duration, 1.0),
                            cx,
                        )
                        .unwrap();
                    (
                        actual_gas,
                        actual.activity().to_vec(),
                        report,
                        expected_gas,
                        expected.activity().to_vec(),
                        [first, second],
                        curved(symmetry).capacity,
                    )
                } else {
                    let old = stack(initial, cx);
                    let (half_gas, half, first) = gas
                        .advance_layered_sdirk2(
                            &old,
                            0.5 * duration,
                            [true; 2],
                            config(),
                            nonlinear(),
                            air::controls(),
                            cx,
                        )
                        .unwrap()
                        .into_parts();
                    let (expected_gas, expected, second) = half_gas
                        .advance_layered_sdirk2(
                            &half,
                            0.5 * duration,
                            [true; 2],
                            config(),
                            nonlinear(),
                            air::controls(),
                            cx,
                        )
                        .unwrap()
                        .into_parts();
                    let (actual_gas, actual, report) = gas
                        .advance_layered_sdirk2_adaptive_to(
                            &old,
                            duration,
                            [true; 2],
                            config(),
                            nonlinear(),
                            air::controls(),
                            time(duration, 1.0),
                            cx,
                        )
                        .unwrap();
                    (
                        actual_gas,
                        actual.activity().to_vec(),
                        report,
                        expected_gas,
                        expected.activity().to_vec(),
                        [first, second],
                        CAPACITY,
                    )
                };
            assert_eq!(actual, expected);
            assert_eq!(actual_gas, expected_gas);
            assert_eq!(report.summary.evaluations, 6);
            assert_eq!(report.summary.accepted_intervals, 1);
            for face in 0..2 {
                let exterior: f64 = receipts
                    .iter()
                    .map(|r| r.transport.endpoint.boundary_water_kg[face])
                    .sum();
                near(report.summary.boundary_water_kg[face], exterior, 2e-14);
            }
            let interface: f64 = receipts
                .iter()
                .map(|r| r.transport.endpoint.interfaces[0].water_left_to_right_kg)
                .sum();
            near(report.layers[0].boundary_water_kg[1], -interface, 2e-14);
            near(report.layers[1].boundary_water_kg[0], interface, 2e-14);
            account(&report, capacity, initial, &actual, duration);
        }
    });
}

#[test]
fn finite_air_sdirk_work_limit_discards_an_accepted_private_prefix() {
    with_cx(|cx| {
        let gas = air::atmosphere(0.8, 1.0);
        let saved = gas.clone();
        let initial = [0.2, 0.1];
        let step = 1.0 / 64.0;
        let mut limited = time(step, 1.0);
        limited.maximum_evaluations = 6;
        let slab = stack(initial, cx);
        let slab_activity = slab.activity().to_vec();
        assert!(matches!(
            gas.advance_layered_sdirk2_adaptive_to(
                &slab, 2.0 * step, [true; 2], config(), nonlinear(),
                air::controls(), limited, cx,
            ),
            Err(MoistureAdvanceError::WorkLimit {
                time_s, accepted_intervals: 1, evaluations: 6,
            }) if time_s == step
        ));
        assert_eq!(slab.time_s(), 0.0);
        assert_eq!(slab.activity(), slab_activity);
        for symmetry in [
            RadialMoistureSymmetry::Sphere,
            RadialMoistureSymmetry::Cylinder { length_m: 0.1 },
        ] {
            let body = coated(symmetry, initial, cx);
            let body_activity = body.activity().to_vec();
            assert!(matches!(
                gas.advance_layered_radial_sdirk2_adaptive_to(
                    &body, 2.0 * step, [false, true], config(), nonlinear(),
                    air::controls(), limited, cx,
                ),
                Err(MoistureAdvanceError::WorkLimit {
                    time_s, accepted_intervals: 1, evaluations: 6,
                }) if time_s == step
            ));
            assert_eq!(body.time_s(), 0.0);
            assert_eq!(body.activity(), body_activity);
        }
        assert_eq!(gas, saved);
    });
}
