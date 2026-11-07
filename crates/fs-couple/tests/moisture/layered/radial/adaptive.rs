//! Physical unequal-shell error scales, exact exchange and whole-duration ledgers.
use super::*;
use crate::atmosphere as air;
use fs_couple::moisture::adaptive::{
    LayeredMoistureAdvanceReport, MoistureAdvanceError, MoistureTimeControls,
};

const RHO: [f64; 2] = [300.0, 900.0];
const SLOPE: [f64; 2] = [0.2, 0.4];
const D: [f64; 2] = [0.001, 0.0005];
const INITIAL: [f64; 2] = [0.8, 0.1];
const CONTACT: f64 = 0.03;

fn coated(symmetry: RadialMoistureSymmetry, cx: &Cx<'_>) -> MoistureLayeredRadial {
    MoistureLayeredRadial::new(
        (0..2)
            .map(|i| {
                radial_layer(
                    linear_card(RHO[i], SLOPE[i], D[i]),
                    symmetry,
                    if i == 0 {
                        vec![0.0, 0.03]
                    } else {
                        vec![0.03, 0.05]
                    },
                    vec![INITIAL[i]],
                    cx,
                )
            })
            .collect(),
        2,
        cx,
    )
    .unwrap()
    .with_transfer(
        [MoistureActivityTransfer::IdealContact; 2],
        vec![MoistureActivityTransfer::ConductanceKgM2S(CONTACT)],
        cx,
    )
    .unwrap()
}

fn time(step: f64, tolerance: f64) -> MoistureTimeControls {
    MoistureTimeControls {
        initial_step_s: step,
        minimum_step_s: 1e-7,
        maximum_step_s: step,
        maximum_evaluations: 5000,
        moisture_absolute_tolerance: tolerance,
        activity_absolute_tolerance: tolerance,
        reservoir_water_absolute_tolerance_kg: 1.0,
        relative_tolerance: 0.0,
        water_balance_tolerance_kg: 1e-10,
    }
}

fn parameters(symmetry: RadialMoistureSymmetry) -> ([f64; 2], f64) {
    let capacity = [
        RHO[0] * SLOPE[0] * measures(symmetry, 0.0, 0.03).0,
        RHO[1] * SLOPE[1] * measures(symmetry, 0.03, 0.05).0,
    ];
    let resistance = measures(symmetry, 0.015, 0.03).1 / (RHO[0] * D[0] * SLOPE[0])
        + 1.0 / (area(symmetry, 0.03) * CONTACT)
        + measures(symmetry, 0.03, 0.04).1 / (RHO[1] * D[1] * SLOPE[1]);
    (
        capacity,
        (1.0 / capacity[0] + 1.0 / capacity[1]) / resistance,
    )
}

fn mode(capacity: [f64; 2], decay: f64) -> [f64; 2] {
    let mean = (capacity[0] * INITIAL[0] + capacity[1] * INITIAL[1]) / (capacity[0] + capacity[1]);
    INITIAL.map(|a| mean + (a - mean) * decay)
}

fn decay(z: f64, second_order: bool) -> f64 {
    if second_order {
        let gamma = 1.0 - core::f64::consts::FRAC_1_SQRT_2;
        (1.0 + (2.0 * gamma - 1.0) * z) / (1.0 + gamma * z).powi(2)
    } else {
        1.0 / (1.0 + z)
    }
}

fn audit<R>(
    old: &MoistureLayeredRadial,
    next: &MoistureLayeredRadial,
    report: &LayeredMoistureAdvanceReport<R>,
) {
    assert_eq!(report.layers.len(), 2);
    assert_eq!(report.layers[0].id, 0);
    assert_eq!(report.layers[1].id, 1);
    assert_eq!(
        report.layers[0].boundary_water_kg[1],
        -report.layers[1].boundary_water_kg[0]
    );
    assert_eq!(
        report.summary.accepted_substeps,
        2 * report.summary.accepted_intervals
    );
    assert!(report.summary.maximum_accepted_error_ratio <= 1.0);
    for (i, layer) in report.layers.iter().enumerate() {
        let field = next.layers()[i].radial();
        let original = old.layers()[i].radial();
        let change: f64 = original
            .moisture()
            .iter()
            .zip(field.moisture())
            .zip(original.dry_masses_kg())
            .map(|((&a, &b), &dry)| dry * (b - a))
            .sum();
        near(layer.water_change_kg, change, 1e-11);
        near(change, layer.boundary_water_kg.iter().sum(), 1e-11);
        assert!(layer.maximum_accepted_error_ratio <= 1.0);
    }
    near(report.summary.transport_balance_residual_kg, 0.0, 1e-10);
}

#[test]
fn curved_layered_accuracy_refines_both_time_methods_against_continuous_exchange() {
    with_cx(|cx| {
        for symmetry in [
            RadialMoistureSymmetry::Sphere,
            RadialMoistureSymmetry::Cylinder { length_m: 0.1 },
        ] {
            let original = coated(symmetry, cx);
            let (capacity, rate) = parameters(symmetry);
            let duration = 0.8 / rate;
            let exact = mode(capacity, (-rate * duration).exp());
            for second_order in [false, true] {
                let mut errors = Vec::new();
                let mut work = Vec::new();
                for tolerance in [1e-3, 1e-5] {
                    let (next, evaluations) = if second_order {
                        let (next, report) = original
                            .advance_sdirk2_adaptive_to(
                                duration,
                                &sealed(),
                                config(),
                                nonlinear(),
                                time(duration, tolerance),
                                cx,
                            )
                            .unwrap();
                        audit(&original, &next, &report);
                        assert_eq!(report.summary.boundary_water_kg, [0.0; 2]);
                        (next, report.summary.evaluations)
                    } else {
                        let (next, report) = original
                            .advance_adaptive_to(
                                duration,
                                &sealed(),
                                config(),
                                nonlinear(),
                                time(duration, tolerance),
                                cx,
                            )
                            .unwrap();
                        audit(&original, &next, &report);
                        assert_eq!(report.summary.boundary_water_kg, [0.0; 2]);
                        (next, report.summary.evaluations)
                    };
                    assert_eq!(next.time_s(), duration);
                    near(next.water_mass_kg(), original.water_mass_kg(), 1e-10);
                    errors.push(
                        next.activity()
                            .iter()
                            .zip(exact)
                            .map(|(&a, b)| (a - b).abs())
                            .fold(0.0_f64, f64::max),
                    );
                    work.push(evaluations);
                }
                assert!(
                    errors[1] < 0.3 * errors[0],
                    "{symmetry:?} second_order={second_order}: {errors:?}"
                );
                assert!(work[1] > work[0]);
            }
            assert_eq!(original.time_s(), 0.0);
            for (&a, b) in original.activity().iter().zip(INITIAL) {
                near(a, b, 1e-15);
            }
        }
    });
}

#[test]
fn layered_radial_water_error_cannot_cancel_across_materials() {
    with_cx(|cx| {
        let symmetry = RadialMoistureSymmetry::Sphere;
        let original = coated(symmetry, cx);
        let (capacity, rate) = parameters(symmetry);
        assert!(capacity[1] > 20.0 * capacity[0]);
        let duration = 0.4 / rate;
        for second_order in [false, true] {
            let coarse = mode(capacity, decay(rate * duration, second_order));
            let fine = mode(capacity, decay(0.5 * rate * duration, second_order).powi(2));
            let discrepancy = (0..2)
                .map(|i| capacity[i] * (coarse[i] - fine[i]).abs())
                .sum::<f64>();
            let fine_factor = if second_order { 1.0 / 3.0 } else { 1.0 };
            let mut controls = time(duration, 1.0);
            controls.reservoir_water_absolute_tolerance_kg = 2.0 * fine_factor * discrepancy;
            if second_order {
                let (next, report) = original
                    .advance_sdirk2_adaptive_to(
                        duration,
                        &sealed(),
                        config(),
                        nonlinear(),
                        controls,
                        cx,
                    )
                    .unwrap();
                audit(&original, &next, &report);
                near(report.summary.maximum_accepted_error_ratio, 0.5, 3e-8);
                assert_eq!(report.summary.accepted_intervals, 1);
                assert_eq!(report.summary.evaluations, 6);
                for (&a, b) in next.activity().iter().zip(fine) {
                    near(a, b, 4e-12);
                }
            } else {
                let (next, report) = original
                    .advance_adaptive_to(duration, &sealed(), config(), nonlinear(), controls, cx)
                    .unwrap();
                audit(&original, &next, &report);
                near(report.summary.maximum_accepted_error_ratio, 0.5, 3e-8);
                assert_eq!(report.summary.accepted_intervals, 1);
                assert_eq!(report.summary.evaluations, 3);
                for (&a, b) in next.activity().iter().zip(fine) {
                    near(a, b, 4e-12);
                }
            }
        }
    });
}

#[test]
fn layered_radial_range_retry_and_accepted_prefix_work_refusal_are_atomic() {
    with_cx(|cx| {
        let original = coated(RadialMoistureSymmetry::Sphere, cx);
        let (_, rate) = parameters(RadialMoistureSymmetry::Sphere);
        let duration = 10.0 / rate;
        assert!(matches!(
            original.advance_sdirk2(duration, &sealed(), config(), nonlinear(), cx),
            Err(MoistureError::ActivityStageRange { .. })
        ));
        let (next, report) = original
            .advance_sdirk2_adaptive_to(
                duration,
                &sealed(),
                config(),
                nonlinear(),
                time(duration, 1e-4),
                cx,
            )
            .unwrap();
        audit(&original, &next, &report);
        assert!(report.summary.rejected_intervals > 0);
        assert_eq!(next.time_s(), duration);
        near(next.water_mass_kg(), original.water_mass_kg(), 1e-10);

        let step = 1.0 / 64.0;
        for second_order in [false, true] {
            let budget = if second_order { 6 } else { 3 };
            let mut controls = time(step, 1.0);
            controls.maximum_evaluations = budget;
            if second_order {
                assert!(matches!(original.advance_sdirk2_adaptive_to(
                    2.0*step,&sealed(),config(),nonlinear(),controls,cx,
                ),Err(MoistureAdvanceError::WorkLimit{time_s,accepted_intervals:1,evaluations})
                    if time_s==step && evaluations==budget));
            } else {
                assert!(matches!(original.advance_adaptive_to(
                    2.0*step,&sealed(),config(),nonlinear(),controls,cx,
                ),Err(MoistureAdvanceError::WorkLimit{time_s,accepted_intervals:1,evaluations})
                    if time_s==step && evaluations==budget));
            }
            assert_eq!(original.time_s(), 0.0);
            for (&a, b) in original.activity().iter().zip(INITIAL) {
                near(a, b, 1e-15);
            }
        }
    });
}

#[test]
fn adaptive_coated_body_and_finite_air_close_every_layer_and_gas_inventory() {
    with_cx(|cx| {
        for symmetry in [
            RadialMoistureSymmetry::Sphere,
            RadialMoistureSymmetry::Cylinder { length_m: 0.1 },
        ] {
            let original = coated(symmetry, cx);
            let atmosphere = air::atmosphere(0.8, 0.5);
            let total_water = atmosphere.vapor_mass_kg().unwrap() + original.water_mass_kg();
            let (next_air, next, report) = atmosphere
                .advance_layered_radial_adaptive_to(
                    &original,
                    2.0,
                    [false, true],
                    config(),
                    nonlinear(),
                    air::controls(),
                    time(0.5, 1e-4),
                    cx,
                )
                .unwrap();
            audit(&original, &next, &report);
            assert_eq!(next.time_s(), 2.0);
            assert!(report.summary.accepted_intervals >= 4);
            assert!(report.summary.water_change_kg > 0.0);
            assert!(
                next_air.relative_humidity().unwrap() < atmosphere.relative_humidity().unwrap()
            );
            assert_eq!(
                report.summary.reservoir_water_kg.unwrap(),
                [
                    atmosphere.vapor_mass_kg().unwrap(),
                    next_air.vapor_mass_kg().unwrap()
                ]
            );
            near(
                next_air.vapor_mass_kg().unwrap() + next.water_mass_kg(),
                total_water,
                1e-10,
            );
            near(
                air::carrier_mass(&next_air),
                air::carrier_mass(&atmosphere),
                1e-12,
            );
            near(report.summary.total_mass_residual_kg.unwrap(), 0.0, 1e-10);
            near(
                report.summary.maximum_carrier_mass_residual_kg.unwrap(),
                0.0,
                1e-10,
            );
            near(
                report.summary.water_change_kg,
                report.summary.boundary_water_kg[1],
                1e-10,
            );
            assert_eq!(report.summary.boundary_water_kg[0], 0.0);
            assert_eq!(original.time_s(), 0.0);
        }
    });
}
