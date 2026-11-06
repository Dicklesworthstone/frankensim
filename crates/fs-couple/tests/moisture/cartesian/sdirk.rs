//! G1 discrete modes and analytic nonlinear rates; G4 physical-stage refusal.
use super::*;
use fs_couple::moisture::adaptive::{MoistureAdvanceError, MoistureTimeControls};

fn decay(z: f64) -> f64 {
    let gamma = 1.0 - core::f64::consts::FRAC_1_SQRT_2;
    (1.0 - (1.0 - 2.0 * gamma) * z) / (1.0 + gamma * z).powi(2)
}

#[test]
fn g1_cartesian_sdirk_has_second_order_anisotropic_mode_decay() {
    with_cx(|cx| {
        let geometry = grid_geometry([1.0, 0.8, 0.6], [4, 3, 2]);
        let d = [0.02, 0.03, 0.01];
        let [nx, ny, nz] = geometry.cells;
        let wave = |i: usize| {
            [i % nx, (i / nx) % ny, i / (nx * ny)]
                .into_iter()
                .zip(geometry.cells)
                .map(|(j, n)| (core::f64::consts::PI * (j as f64 + 0.5) / n as f64).cos())
                .product::<f64>()
        };
        let initial = MoistureCartesian::new(
            directional_card(d),
            point(0.1),
            ["diffusion_x", "diffusion_y", "diffusion_z"],
            MaterialPropertySelection::SingleClaimOnly,
            geometry,
            (0..nx * ny * nz).map(|i| 0.15 + 0.03 * wave(i)).collect(),
            cx,
        ).unwrap();
        let lambda = (0..3).map(|k| {
            let dx = geometry.lengths_m[k] / geometry.cells[k] as f64;
            4.0 * d[k] / (dx * dx)
                * (core::f64::consts::PI / (2.0 * geometry.cells[k] as f64)).sin().powi(2)
        }).sum::<f64>();
        let duration = 2.0;
        let errors = [8, 16, 32].map(|steps| {
            let mut state = initial.clone();
            for _ in 0..steps {
                let dt = duration / steps as f64;
                let (next, report) = state.advance_sdirk2(
                    dt, &sealed(), config(), nonlinear(), cx,
                ).unwrap();
                assert_eq!(report.endpoint.boundary_water_kg, [0.0; 6]);
                assert!(report.endpoint.nonlinear_relative_residual < 1e-12);
                near(next.water_mass_kg(), initial.water_mass_kg(), 2e-11);
                state = next;
            }
            assert_eq!(state.time_s(), duration);
            state.moisture().iter().enumerate().map(|(i, &m)| {
                (m - (0.15 + 0.03 * (-lambda * duration).exp() * wave(i))).abs()
            }).fold(0.0, f64::max)
        });
        assert!(errors[0] / errors[1] > 3.8, "{errors:?}");
        assert!(errors[1] / errors[2] > 3.8, "{errors:?}");
        assert!(errors[1] / errors[2] < 4.5, "{errors:?}");
        assert_eq!(initial.time_s(), 0.0);
    });
}

#[test]
fn g1_cartesian_sdirk_accounts_for_six_distinct_surface_fluxes() {
    with_cx(|cx| {
        let geometry = grid_geometry([1.0, 2.0, 0.5], [1; 3]);
        let d = [0.01, 0.02, 0.04];
        let initial = MoistureCartesian::new(
            directional_card(d),
            point(0.1),
            ["diffusion_x", "diffusion_y", "diffusion_z"],
            MaterialPropertySelection::SingleClaimOnly,
            geometry,
            vec![0.1],
            cx,
        ).unwrap();
        let reservoir = [0.08, 0.21, 0.12, 0.17, 0.09, 0.24];
        let rates: [f64; 6] = std::array::from_fn(|side| {
            2.0 * d[side / 2] / geometry.lengths_m[side / 2].powi(2)
        });
        let lambda: f64 = rates.iter().sum();
        let equilibrium: f64 =
            rates.iter().zip(reservoir).map(|(a, m)| a * m).sum::<f64>() / lambda;
        let gamma = 1.0 - core::f64::consts::FRAC_1_SQRT_2;
        let dt = 2.0;
        let first = equilibrium + (0.1 - equilibrium) / (1.0 + gamma * dt * lambda);
        let expected = equilibrium + (0.1 - equilibrium) * decay(dt * lambda);
        let boundaries = reservoir.map(MoistureFace::PrescribedDryBasis);
        let (next, report) = initial.advance_sdirk2(
            dt, &boundaries, config(), nonlinear(), cx,
        ).unwrap();
        near(next.moisture()[0], expected, 3e-13);
        assert_eq!(report.endpoint.interval_s, [0.0, dt]);
        assert_eq!(report.first_stage.interval_s, [0.0, gamma * dt]);
        for side in 0..6 {
            let water = initial.dry_mass_kg() * dt * rates[side]
                * ((1.0 - gamma) * (reservoir[side] - first)
                    + gamma * (reservoir[side] - expected));
            near(report.endpoint.boundary_water_kg[side], water, 2e-10);
            assert_eq!(report.endpoint.surface_moisture[side], vec![(0, reservoir[side])]);
        }
        near(report.endpoint.water_balance_residual_kg, 0.0, 1e-11);
        // The optional method does not replace ordinary backward Euler.
        let ordinary = next.advance(0.5, &boundaries, config(), nonlinear(), cx).unwrap().0;
        near(ordinary.moisture()[0],
            equilibrium + (expected - equilibrium) / (1.0 + 0.5 * lambda), 3e-13);
        let replay = initial.advance_sdirk2(
            dt, &boundaries, config(), nonlinear(), cx,
        ).unwrap();
        assert_eq!(replay.0.moisture(), next.moisture());
    });
}

#[test]
fn g1_cartesian_sdirk_requeries_nonlinear_diffusion_at_both_stages() {
    with_cx(|cx| {
        let initial = grid(card([0.01, 0.04], 0.3),
            grid_geometry([1.0; 3], [1; 3]), vec![0.1], cx);
        let gamma = 1.0 - core::f64::consts::FRAC_1_SQRT_2;
        let dt = 2.0;
        // Six independent half-cell faces: F(M)=12*(.01+.1*M)*(.2-M).
        let root = |base: f64, h: f64| {
            let a = 1.2 * h;
            let b = 1.0 - 0.12 * h;
            let c = base + 0.024 * h;
            2.0 * c / (b + (b * b + 4.0 * a * c).sqrt())
        };
        let first = root(0.1, gamma * dt);
        let rate = 12.0 * (0.01 + 0.1 * first) * (0.2 - first);
        let expected = root(0.1 + (1.0 - gamma) * dt * rate, gamma * dt);
        let boundaries = std::array::from_fn(|_| MoistureFace::PrescribedDryBasis(0.2));
        let (next, report) = initial.advance_sdirk2(
            dt, &boundaries, config(), nonlinear(), cx,
        ).unwrap();
        near(next.moisture()[0], expected, 2e-12);
        near(report.first_stage.transport_end[0].property("moisture_diffusivity")
            .unwrap().value_si(), 0.01 + 0.1 * first, 2e-12);
        near(report.endpoint.transport_end[0].property("moisture_diffusivity")
            .unwrap().value_si(), 0.01 + 0.1 * expected, 2e-12);
        near(report.endpoint.water_balance_residual_kg, 0.0, 1e-11);
        assert_eq!(initial.moisture(), &[0.1]);
    });
}

#[test]
fn g1_cartesian_sdirk_kirchhoff_wetting_matches_analytic_dry_front() {
    with_cx(|cx| {
        let initial = MoistureCartesian::new_kirchhoff(
            card([0.0, 0.03], 0.3),
            point(0.1),
            ["moisture_diffusivity"; 3],
            MaterialPropertySelection::SingleClaimOnly,
            grid_geometry([1.0; 3], [1; 3]),
            vec![0.0],
            cx,
        ).unwrap();
        let boundaries = std::array::from_fn(|_| MoistureFace::PrescribedDryBasis(0.2));
        // D(M)=.1*M, six Kirchhoff faces: M'=.6*(.2^2-M^2).
        // The exact one-cell ODE solution is .2*tanh(.12*t).
        let errors = [4, 8, 16].map(|steps| {
            let mut state = initial.clone();
            for _ in 0..steps {
                let (next, report) = state.advance_sdirk2(
                    3.0 / steps as f64, &boundaries, config(), nonlinear(), cx,
                ).unwrap();
                assert!(next.moisture()[0] > state.moisture()[0]);
                assert!(next.moisture()[0] < 0.2);
                assert!(report.endpoint.kirchhoff_sources.is_some());
                near(report.endpoint.water_balance_residual_kg, 0.0, 1e-11);
                state = next;
            }
            (state.moisture()[0] - 0.2 * (0.36_f64).tanh()).abs()
        });
        assert!(errors[0] / errors[1] > 3.5, "{errors:?}");
        assert!(errors[1] / errors[2] > 3.5, "{errors:?}");
        assert_eq!(initial.moisture(), &[0.0]);
    });
}

#[test]
fn g4_cartesian_sdirk_dormant_reservoirs_cannot_authorize_overshoot() {
    with_cx(|cx| {
        let initial = MoistureCartesian::new(
            directional_card([0.01, 0.0, 0.0]),
            point(0.1),
            ["diffusion_x", "diffusion_y", "diffusion_z"],
            MaterialPropertySelection::SingleClaimOnly,
            grid_geometry([1.0; 3], [1; 3]),
            vec![0.1],
            cx,
        ).unwrap();
        let mut boundaries = sealed();
        boundaries[0] = MoistureFace::PrescribedDryBasis(0.2);
        boundaries[1] = MoistureFace::SourceSorption {
            law: super::super::surface::law([0.0; 2]),
            relative_humidity: 1.0,
        };
        boundaries[2] = MoistureFace::PrescribedDryBasis(0.3);
        // z=3 from the sole active face produces an inadmissible SDIRK
        // overshoot above .2; neither dormant .25 nor .3 may admit it.
        let error = initial.advance_sdirk2(
            150.0, &boundaries, config(), nonlinear(), cx,
        ).unwrap_err();
        assert!(matches!(error,
            MoistureError::StageRange { bounds, .. } if bounds == [0.1, 0.2]
        ), "{error:?}");
        assert_eq!(initial.moisture(), &[0.1]);
        assert_eq!(initial.time_s(), 0.0);
        let (next, _) = initial.advance_sdirk2(
            25.0, &boundaries, config(), nonlinear(), cx,
        ).unwrap();
        assert!(next.moisture()[0] > 0.1 && next.moisture()[0] < 0.2);
    });
}

#[test]
fn g0_cartesian_sdirk_preserves_stiff_equilibrium_and_cancellation() {
    with_cx(|cx| {
        let initial = grid(card([0.01; 2], 0.3),
            grid_geometry([1.0; 3], [2; 3]), vec![0.15; 8], cx);
        let (next, report) = initial.advance_sdirk2(
            1e8, &sealed(), config(), nonlinear(), cx,
        ).unwrap();
        assert_eq!(next.moisture(), initial.moisture());
        assert_eq!(report.endpoint.boundary_water_kg, [0.0; 6]);
        assert_eq!(report.endpoint.linear_iterations, 0);
        let gate = CancelGate::new();
        gate.request();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            let cancelled = Cx::new(&gate, arena, StreamKey {
                seed: 3, kernel_id: 17, tile: 0, iteration: 0,
            }, Budget::INFINITE, ExecMode::Deterministic);
            assert!(matches!(initial.advance_sdirk2(
                1.0, &sealed(), config(), nonlinear(), &cancelled,
            ), Err(MoistureError::Cancelled)));
        });
        assert_eq!(initial.time_s(), 0.0);
    });
}

fn time_controls(dt: f64, absolute: f64) -> MoistureTimeControls {
    MoistureTimeControls {
        initial_step_s: dt,
        minimum_step_s: dt / 4096.0,
        maximum_step_s: dt,
        maximum_evaluations: 100_000,
        moisture_absolute_tolerance: absolute,
        activity_absolute_tolerance: 1e-8,
        reservoir_water_absolute_tolerance_kg: 1e-9,
        relative_tolerance: 0.0,
        water_balance_tolerance_kg: 1e-9,
    }
}

#[test]
fn g1_cartesian_sdirk_adaptive_uses_cubic_error_and_weighted_water() {
    with_cx(|cx| {
        let initial = grid(card([0.01; 2], 0.3),
            grid_geometry([1.0; 3], [1; 3]), vec![0.1], cx);
        let boundaries = std::array::from_fn(|_| MoistureFace::PrescribedDryBasis(0.2));
        // Six half-cell faces give lambda=12D/L^2=.12.
        let coarse = 0.2 - 0.1 * decay(0.12);
        let fine = 0.2 - 0.1 * decay(0.06).powi(2);
        let mut controls = time_controls(1.0, (coarse - fine).abs() / (3.0 * 0.2));
        controls.maximum_step_s = 2.0;
        let (next, report) = initial.advance_sdirk2_adaptive_to(
            1.0, &boundaries, config(), nonlinear(), controls, cx,
        ).unwrap();
        near(next.moisture()[0], fine, 3e-13);
        near(report.maximum_accepted_error_ratio, 0.2, 2e-8);
        assert_eq!(report.comparisons, 1);
        assert_eq!(report.accepted_intervals, 1);
        assert_eq!(report.accepted_substeps, 2);
        assert_eq!(report.evaluations, 6);
        assert_eq!(report.rejected_intervals, 0);
        assert_eq!(report.suggested_step_s, 1.0);
        assert_eq!(report.last_substep.endpoint.interval_s, [0.5, 1.0]);
        let water = initial.dry_mass_kg() * (fine - 0.1);
        for transfer in report.boundary_water_kg {
            near(transfer, water / 6.0, 2e-10);
        }
        near(report.water_change_kg, water, 2e-10);
        near(report.water_balance_residual_kg, 0.0, 1e-9);
        controls.moisture_absolute_tolerance *= 2.0;
        let (_, growing) = initial.advance_sdirk2_adaptive_to(
            1.0, &boundaries, config(), nonlinear(), controls, cx,
        ).unwrap();
        assert_eq!(growing.suggested_step_s, 2.0);
        assert_eq!(initial.time_s(), 0.0);
    });
}

#[test]
fn g4_cartesian_sdirk_adaptive_recovers_stiff_drying_without_clipping() {
    with_cx(|cx| {
        let initial = grid(card([0.01; 2], 0.3),
            grid_geometry([1.0; 3], [1; 3]), vec![0.2], cx);
        let boundaries = std::array::from_fn(|_| MoistureFace::PrescribedDryBasis(0.0));
        let mut controls = time_controls(50.0, 1e-5);
        let (next, report) = initial.advance_sdirk2_adaptive_to(
            50.0, &boundaries, config(), nonlinear(), controls, cx,
        ).unwrap();
        near(next.moisture()[0], 0.2 * (-6.0_f64).exp(), 1e-4);
        assert!(next.moisture()[0] > 0.0);
        assert_eq!(next.time_s(), 50.0);
        assert!(report.rejected_intervals > 0);
        assert!(report.evaluations < 6 * report.comparisons);
        assert_eq!(report.evaluations % 2, 0);
        assert_eq!(report.accepted_substeps, 2 * report.accepted_intervals);
        near(report.boundary_water_kg.iter().sum(),
            initial.dry_mass_kg() * (next.moisture()[0] - 0.2), 1e-9);
        controls.minimum_step_s = 50.0;
        assert!(matches!(initial.advance_sdirk2_adaptive_to(
            50.0, &boundaries, config(), nonlinear(), controls, cx,
        ), Err(MoistureAdvanceError::Step {
            evaluations: 2, error: MoistureError::StageRange { .. }, ..
        })));
        assert_eq!(initial.moisture(), &[0.2]);
        assert_eq!(initial.time_s(), 0.0);
    });
}

#[test]
fn g4_cartesian_sdirk_adaptive_charges_complete_stages_and_discards_prefixes() {
    with_cx(|cx| {
        let initial = grid(card([0.01; 2], 0.3),
            grid_geometry([1.0; 3], [2; 3]), vec![0.1; 8], cx);
        let mut controls = time_controls(1.0, 1e-4);
        controls.maximum_evaluations = 5;
        assert!(matches!(initial.advance_sdirk2_adaptive_to(
            1.0, &sealed(), config(), nonlinear(), controls, cx,
        ), Err(MoistureAdvanceError::WorkLimit {
            evaluations: 4, accepted_intervals: 0, ..
        })));
        controls.maximum_evaluations = 11;
        assert!(matches!(initial.advance_sdirk2_adaptive_to(
            2.0, &sealed(), config(), nonlinear(), controls, cx,
        ), Err(MoistureAdvanceError::WorkLimit {
            evaluations: 10, accepted_intervals: 1, ..
        })));
        controls.maximum_evaluations = 12;
        let (next, report) = initial.advance_sdirk2_adaptive_to(
            2.0, &sealed(), config(), nonlinear(), controls, cx,
        ).unwrap();
        assert_eq!(report.evaluations, 12);
        assert_eq!(report.accepted_intervals, 2);
        assert_eq!(report.accepted_substeps, 4);
        assert_eq!(next.moisture(), initial.moisture());
        assert_eq!(next.time_s(), 2.0);
        assert_eq!(initial.time_s(), 0.0);
    });
}
