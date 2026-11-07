//! Independent continuous partition exchange under the native order-two driver.
use super::*;

#[test]
fn layered_sdirk_accuracy_refines_the_physical_exchange_and_each_material_ledger() {
    with_cx(|cx| {
        let initial = [0.9, 0.1];
        let contact = 0.008;
        let original = stack(initial, [0.2, 0.4], cx)
            .with_transfer(
                [MoistureActivityTransfer::IdealContact; 2],
                vec![MoistureActivityTransfer::ConductanceKgM2S(contact / 0.001)],
                cx,
            )
            .unwrap();
        let conductance = 1.0 / (1.0 / HALF[0] + 1.0 / contact + 1.0 / HALF[1]);
        let mean =
            (CAPACITY[0] * initial[0] + CAPACITY[1] * initial[1]) / (CAPACITY[0] + CAPACITY[1]);
        let exact = initial.map(|a| {
            mean + (a - mean) * (-conductance * (1.0 / CAPACITY[0] + 1.0 / CAPACITY[1])).exp()
        });
        let mut errors = Vec::new();
        let mut work = Vec::new();
        for tolerance in [2e-3, 2e-5] {
            let (next, report) = original
                .advance_sdirk2_adaptive_to(
                    1.0,
                    &sealed(),
                    config(),
                    nonlinear(),
                    controls(1.0, tolerance),
                    cx,
                )
                .unwrap();
            errors.push(
                next.activity()
                    .iter()
                    .zip(exact)
                    .map(|(&a, b)| (a - b).abs())
                    .fold(0.0_f64, f64::max),
            );
            work.push(report.summary.evaluations);
            assert_eq!(next.time_s(), 1.0);
            assert_eq!(report.summary.interval_s, [0.0, 1.0]);
            assert_eq!(report.summary.boundary_water_kg, [0.0; 2]);
            assert_eq!(
                report.summary.accepted_substeps,
                2 * report.summary.accepted_intervals
            );
            assert!(report.summary.maximum_accepted_error_ratio <= 1.0);
            assert_eq!(report.summary.evaluations % 2, 0);
            assert_eq!(report.layers.len(), 2);
            assert_eq!(
                report.layers[0].boundary_water_kg[1],
                -report.layers[1].boundary_water_kg[0]
            );
            for layer in &report.layers {
                near(
                    layer.water_change_kg,
                    layer.boundary_water_kg.iter().sum(),
                    1e-11,
                );
                assert!(layer.maximum_accepted_error_ratio <= 1.0);
            }
            near(next.water_mass_kg(), original.water_mass_kg(), 1e-11);
            near(report.summary.transport_balance_residual_kg, 0.0, 1e-11);
        }
        assert!(errors[1] < 0.3 * errors[0], "{errors:?}");
        assert!(work[1] > work[0]);
        assert_eq!(original.time_s(), 0.0);
        for (a, b) in original.activity().iter().zip(initial) {
            near(*a, b, 1e-15);
        }
    });
}
