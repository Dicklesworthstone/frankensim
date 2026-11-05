//! Independent continuous-time throughflow, native spatial parity and refusal.
use super::*;
use fs_couple::moisture::adaptive::{MoistureAdvanceError, MoistureTimeControls};

fn time_controls(maximum: f64) -> MoistureTimeControls {
    MoistureTimeControls {
        initial_step_s: maximum,
        minimum_step_s: 1e-6,
        maximum_step_s: maximum,
        maximum_evaluations: 5000,
        // A uniform specimen and constant total vapor must not decide accuracy
        // for the independent chamber-to-chamber exchange in these fixtures.
        moisture_absolute_tolerance: 1.0,
        activity_absolute_tolerance: 1.0,
        reservoir_water_absolute_tolerance_kg: 1.0,
        relative_tolerance: 0.0,
        water_balance_tolerance_kg: 1e-10,
    }
}

#[test]
fn g1_two_chamber_accuracy_resolves_throughflow_hidden_by_constant_total_storage() {
    with_cx(|cx| {
        let gas = [atmosphere(0.8, 1.0), atmosphere(0.2, 1.0)];
        let old = slab(card([0.01; 2], 0.3), 1.0, vec![0.125], cx);
        let faces = [
            port(MoistureAtmospherePort::First, 0.25),
            port(MoistureAtmospherePort::Second, 0.25),
        ];
        let target = 2.0;
        let capacity = gas[0].vapor_mass_kg().unwrap() / 0.8;
        // Symmetric chambers around M/s=0.5 keep M and total vapor constant.
        // C*dr/dt = -m_dry*(2D/L²)*s*(r-0.5), independent of backward Euler.
        let rate = old.dry_mass_kg() * 2.0 * 0.01 * 0.25 / capacity;
        let exact = 0.5 + 0.3 * (-rate * target).exp();
        let run = |time| {
            gas[0]
                .advance_paired_adaptive_to(
                    &gas[1],
                    &old,
                    target,
                    &faces,
                    config(),
                    pair_controls(),
                    time,
                    cx,
                )
                .unwrap()
        };
        let (loose_air, _, loose) = run(time_controls(target));
        let loose_error = (loose_air[0].relative_humidity().unwrap() - exact).abs();
        assert_eq!(loose.transport.evaluations, 3);
        for water_only in [false, true] {
            let mut controls = time_controls(target);
            if water_only {
                controls.reservoir_water_absolute_tolerance_kg = capacity * 1e-5;
            } else {
                controls.activity_absolute_tolerance = 1e-5;
            }
            let (air, field, report) = run(controls);
            let summary = &report.transport;
            let error = (air[0].relative_humidity().unwrap() - exact).abs();
            assert!(error < 0.2 * loose_error);
            assert!(summary.evaluations > loose.transport.evaluations);
            assert!(summary.rejected_intervals > 0);
            assert!(summary.maximum_accepted_error_ratio <= 1.0);
            assert_eq!(summary.accepted_substeps, 2 * summary.accepted_intervals);
            assert_eq!(summary.interval_s, [0.0, target]);
            assert_eq!(field.time_s(), target);
            assert_eq!(summary.activity, None);
            near(field.moisture()[0], 0.125, 2e-11);
            near(
                summary.water_change_kg,
                0.0,
                controls.water_balance_tolerance_kg,
            );
            let [before, after] = summary.reservoir_water_kg.unwrap();
            near(after, before, controls.water_balance_tolerance_kg);
            assert!(report.atmospheres[0].water_to_solids_kg > 1e-4);
            assert!(report.atmospheres[1].water_to_solids_kg < -1e-4);
            near(
                report.atmospheres[0].water_to_solids_kg,
                -report.atmospheres[1].water_to_solids_kg,
                controls.water_balance_tolerance_kg,
            );
            for i in 0..2 {
                let account = &report.atmospheres[i];
                near(
                    account.water_to_solids_kg,
                    summary.boundary_water_kg[i],
                    0.0,
                );
                near(
                    account.vapor_mass_kg[1] - account.vapor_mass_kg[0]
                        + account.water_to_solids_kg,
                    0.0,
                    controls.water_balance_tolerance_kg,
                );
                near(
                    account.water_balance_residual_kg,
                    0.0,
                    controls.water_balance_tolerance_kg,
                );
                near(
                    account.total_mass_residual_kg,
                    0.0,
                    controls.water_balance_tolerance_kg,
                );
                near(carrier_mass(&air[i]), carrier_mass(&gas[i]), 2e-12);
                assert_eq!(air[i].source_catalog(), gas[i].source_catalog());
            }
            let (repeat_air, repeat, repeated) = run(controls);
            assert_eq!(repeat_air, air);
            assert_eq!(repeat.moisture(), field.moisture());
            assert_eq!(repeated.transport.evaluations, summary.evaluations);
            assert_eq!(
                repeated.transport.boundary_water_kg,
                summary.boundary_water_kg
            );

            // A genuine 3-D field with the same axial cells/area and sealed
            // transverse faces must resolve the same two finite inventories.
            if !water_only {
                let spatial = MoistureCartesian::new(
                    card([0.01; 2], 0.3),
                    point(0.1),
                    ["moisture_diffusivity"; 3],
                    MaterialPropertySelection::SingleClaimOnly,
                    MoistureCartesianGeometry {
                        lengths_m: [1.0, 0.02, 0.05],
                        cells: [1; 3],
                        frame: ContentHash([4; 32]),
                        maximum_moisture: 0.3,
                        maximum_cells: 1,
                    },
                    vec![0.125],
                    cx,
                )
                .unwrap();
                let spatial_faces = [faces[0].clone(), faces[1].clone(), None, None, None, None];
                let (spatial_air, spatial_next, spatial_report) = gas[0]
                    .advance_paired_cartesian_adaptive_to(
                        &gas[1],
                        &spatial,
                        target,
                        &spatial_faces,
                        config(),
                        nonlinear(),
                        pair_controls(),
                        controls,
                        cx,
                    )
                    .unwrap();
                near(spatial_next.moisture()[0], field.moisture()[0], 2e-10);
                for i in 0..2 {
                    near(
                        spatial_air[i].relative_humidity().unwrap(),
                        air[i].relative_humidity().unwrap(),
                        2e-8,
                    );
                    near(
                        spatial_report.atmospheres[i].water_to_solids_kg,
                        report.atmospheres[i].water_to_solids_kg,
                        2e-10,
                    );
                }
                assert_eq!(&spatial_report.transport.boundary_water_kg[2..], &[0.0; 4]);
            }
        }
        assert_eq!(old.moisture(), &[0.125]);
        assert_eq!(old.time_s(), 0.0);
    });
}

#[test]
fn g4_paired_duration_keeps_sealed_gas_and_discards_late_work_or_cancelled_proposals() {
    with_cx(|cx| {
        let gas = [atmosphere(0.8, 1.0), atmosphere(0.2, 0.7)];
        let saved = gas.clone();
        let old = slab(card([0.01; 2], 0.3), 1.0, vec![0.1; 2], cx);
        let faces = [port(MoistureAtmospherePort::First, 0.25), None];
        let mut controls = time_controls(1.0);
        let (air, field, report) = gas[0]
            .advance_paired_adaptive_to(
                &gas[1],
                &old,
                2.0,
                &faces,
                config(),
                pair_controls(),
                controls,
                cx,
            )
            .unwrap();
        assert_eq!(air[1], gas[1]);
        assert_eq!(report.atmospheres[1].water_to_solids_kg, 0.0);
        assert_eq!(report.transport.boundary_water_kg[1], 0.0);
        assert_eq!(field.time_s(), 2.0);
        controls.maximum_evaluations = 3;
        assert!(matches!(
            gas[0].advance_paired_adaptive_to(
                &gas[1],
                &old,
                2.0,
                &faces,
                config(),
                pair_controls(),
                controls,
                cx,
            ),
            Err(MoistureAdvanceError::WorkLimit {
                time_s,
                accepted_intervals: 1,
                evaluations: 3,
            }) if time_s == 1.0
        ));
        assert_eq!(gas, saved);
        assert_eq!(old.moisture(), &[0.1; 2]);
        assert_eq!(old.time_s(), 0.0);
    });
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 3,
                kernel_id: 17,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let gas = [atmosphere(0.8, 1.0), atmosphere(0.2, 0.7)];
        let saved = gas.clone();
        let field = slab(card([0.01; 2], 0.3), 1.0, vec![0.1], &cx);
        gate.request();
        assert!(matches!(
            gas[0].advance_paired_adaptive_to(
                &gas[1],
                &field,
                2.0,
                &[port(MoistureAtmospherePort::First, 0.25), None],
                config(),
                pair_controls(),
                time_controls(1.0),
                &cx,
            ),
            Err(MoistureAdvanceError::Cancelled)
        ));
        assert_eq!(gas, saved);
        assert_eq!(field.moisture(), &[0.1]);
        assert_eq!(field.time_s(), 0.0);
    });
}
