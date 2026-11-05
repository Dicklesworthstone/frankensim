//! Independent scalar/eigenmode references for conservative second-order time.
use super::*;

fn nonlinear() -> MoistureNonlinearConfig {
    MoistureNonlinearConfig {
        relative_tolerance: 1e-12,
        maximum_iterations: 128,
        relaxation: 1.0,
    }
}

#[test]
fn g1_sdirk2_neumann_mode_has_second_order_decay_and_closed_water() {
    with_cx(|cx| {
        let n = 16;
        let diffusivity = 0.02;
        let total = 8.0;
        let wave = |i: usize| (core::f64::consts::PI * (i as f64 + 0.5) / n as f64).cos();
        let original = slab(
            card([diffusivity; 2], 0.3),
            1.0,
            (0..n).map(|i| 0.15 + 0.05 * wave(i)).collect(),
            cx,
        );
        // Exact eigenvalue of the independent centered, sealed FV Laplacian.
        let lambda = 4.0
            * diffusivity
            * (n * n) as f64
            * (core::f64::consts::PI / (2.0 * n as f64)).sin().powi(2);
        let exact_amplitude = 0.05 * (-lambda * total).exp();
        let errors = [4, 8, 16].map(|steps| {
            let mut state = original.clone();
            for _ in 0..steps {
                let (next, report) = state
                    .advance_sdirk2(
                        total / steps as f64,
                        &[MoistureFace::Sealed, MoistureFace::Sealed],
                        config(),
                        nonlinear(),
                        cx,
                    )
                    .unwrap();
                assert_eq!(report.endpoint.boundary_water_kg, [0.0; 2]);
                near(next.water_mass_kg(), original.water_mass_kg(), 2e-13);
                assert!(report.endpoint.nonlinear.unwrap().relative_residual < 1e-12);
                state = next;
            }
            assert_eq!(state.time_s(), total);
            state
                .moisture()
                .iter()
                .enumerate()
                .map(|(i, &m)| (m - (0.15 + exact_amplitude * wave(i))).abs())
                .fold(0.0, f64::max)
        });
        assert!(errors[0] / errors[1] > 3.8, "{errors:?}");
        assert!(errors[1] / errors[2] > 3.8, "{errors:?}");
        assert!(errors[1] / errors[2] < 4.5, "{errors:?}");
        assert_eq!(original.time_s(), 0.0);
    });
}

#[test]
fn g1_sdirk2_nonlinear_stages_use_original_rates_and_weighted_surface_water() {
    with_cx(|cx| {
        let initial = slab(card([0.01, 0.04], 0.3), 1.0, vec![0.1], cx);
        let dt = 2.0;
        let gamma = 1.0 - core::f64::consts::FRAC_1_SQRT_2;
        // Independent positive root of x-base=h*4*(.01+.1*x)*(.2-x).
        let root = |base: f64, h: f64| {
            let a = 0.4 * h;
            let b = 1.0 - 0.04 * h;
            let c = base + 0.008 * h;
            2.0 * c / (b + (b * b + 4.0 * a * c).sqrt())
        };
        let first = root(0.1, gamma * dt);
        let rate = |m: f64| 4.0 * (0.01 + 0.1 * m) * (0.2 - m);
        let expected = root(0.1 + (1.0 - gamma) * dt * rate(first), gamma * dt);
        let (next, report) = initial
            .advance_sdirk2(dt, &faces(0.2), config(), nonlinear(), cx)
            .unwrap();
        near(next.moisture()[0], expected, 2e-12);
        near(report.first_stage.surface_moisture[0].unwrap(), 0.2, 1e-15);
        let first_water = initial.dry_mass_kg() * gamma * dt * rate(first);
        near(
            report.first_stage.boundary_water_kg.iter().sum(),
            first_water,
            2e-12,
        );
        let total_water =
            initial.dry_mass_kg() * dt * ((1.0 - gamma) * rate(first) + gamma * rate(expected));
        for actual in report.endpoint.boundary_water_kg {
            near(actual, 0.5 * total_water, 2e-12);
        }
        near(report.endpoint.water_balance_residual_kg, 0.0, 2e-12);
        assert_eq!(report.endpoint.interval_s, [0.0, dt]);
        assert_eq!(report.first_stage.interval_s, [0.0, gamma * dt]);
        // Calling this explicit method must not silently change ordinary BE.
        let ordinary = next.advance(0.1, &faces(0.2), config(), cx).unwrap();
        assert!(ordinary.1.nonlinear.is_none());
    });
}

#[test]
fn g1_sdirk2_kirchhoff_dry_front_converges_to_independent_nonlinear_ode() {
    with_cx(|cx| {
        let initial = MoistureSlab::new_kirchhoff(
            card([0.0, 0.03], 0.3),
            point(0.1),
            "moisture_diffusivity",
            MaterialPropertySelection::SingleClaimOnly,
            geometry(1.0),
            vec![0.0],
            nonlinear(),
            cx,
        )
        .unwrap();
        // D(M)=.1*M, so two exact half-cell potentials give
        // dM/dt=.2*(.2²-M²), with M(0)=0 and M=.2*tanh(.04*t).
        let exact = 0.2 * (0.4_f64).tanh();
        let errors = [4, 8, 16].map(|steps| {
            let mut state = initial.clone();
            for _ in 0..steps {
                let (next, report) = state
                    .advance_sdirk2(10.0 / steps as f64, &faces(0.2), config(), nonlinear(), cx)
                    .unwrap();
                assert!(next.moisture()[0] > state.moisture()[0]);
                assert!(next.moisture()[0] < 0.2);
                assert!(report.endpoint.kirchhoff_source.is_some());
                near(report.endpoint.water_balance_residual_kg, 0.0, 2e-12);
                state = next;
            }
            (state.moisture()[0] - exact).abs()
        });
        assert!(errors[0] / errors[1] > 3.5, "{errors:?}");
        assert!(errors[1] / errors[2] > 3.5, "{errors:?}");
    });
}

#[test]
fn g4_sdirk2_refuses_nonphysical_stiff_drying_without_clipping_or_mutation() {
    with_cx(|cx| {
        let initial = slab(card([0.01; 2], 0.3), 1.0, vec![0.2], cx);
        let error = initial
            .advance_sdirk2(100.0, &faces(0.0), config(), nonlinear(), cx)
            .unwrap_err();
        assert!(
            matches!(error, MoistureError::StageRange { .. }),
            "{error:?}"
        );
        assert_eq!(initial.moisture(), &[0.2]);
        assert_eq!(initial.time_s(), 0.0);
        let (next, report) = initial
            .advance_sdirk2(10.0, &faces(0.0), config(), nonlinear(), cx)
            .unwrap();
        assert!(next.moisture()[0] > 0.0 && next.moisture()[0] < 0.2);
        let gamma = 1.0 - core::f64::consts::FRAC_1_SQRT_2;
        let z = 0.4; // dt*4D/L² for two equal prescribed faces.
        let decay = (1.0 - (1.0 - 2.0 * gamma) * z) / (1.0 + gamma * z).powi(2);
        near(next.moisture()[0], 0.2 * decay, 2e-13);
        assert!(report.endpoint.boundary_water_kg.iter().all(|&m| m < 0.0));
        let replay = initial
            .advance_sdirk2(10.0, &faces(0.0), config(), nonlinear(), cx)
            .unwrap();
        assert_eq!(replay.0.moisture(), next.moisture());
    });
}

#[test]
fn g4_zero_transfer_reservoir_cannot_authorize_sdirk_overshoot() {
    with_cx(|cx| {
        let initial = slab(card([0.01; 2], 0.3), 1.0, vec![0.1], cx);
        let boundaries = [
            MoistureFace::PrescribedDryBasis(0.2),
            MoistureFace::SourceSorption {
                law: super::surface::law([0.0; 2]),
                relative_humidity: 1.0, // M_eq=.25, but this face is sealed.
            },
        ];
        // One active face: lambda=2D/L², so z=3. The unconstrained
        // second-order endpoint exceeds .2 while remaining below .25.
        let error = initial
            .advance_sdirk2(150.0, &boundaries, config(), nonlinear(), cx)
            .unwrap_err();
        assert!(
            matches!(error,
                MoistureError::StageRange { bounds, .. } if bounds == [0.1,0.2]
            ),
            "{error:?}"
        );
        let (next, report) = initial
            .advance_sdirk2(25.0, &boundaries, config(), nonlinear(), cx)
            .unwrap();
        assert!(next.moisture()[0] > 0.1 && next.moisture()[0] < 0.2);
        assert_eq!(report.endpoint.boundary_water_kg[1], 0.0);
        assert_eq!(initial.moisture(), &[0.1]);
        assert_eq!(initial.time_s(), 0.0);
    });
}

#[test]
fn g4_sdirk2_cancelled_or_invalid_controls_publish_no_state() {
    with_cx(|cx| {
        let initial = slab(card([0.01; 2], 0.3), 1.0, vec![0.1; 8], cx);
        let mut invalid = nonlinear();
        invalid.maximum_iterations = 0;
        assert!(
            initial
                .advance_sdirk2(1.0, &faces(0.2), config(), invalid, cx)
                .is_err()
        );
        let gate = CancelGate::new();
        gate.request();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            let cancelled = Cx::new(
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
            assert!(matches!(
                initial.advance_sdirk2(1.0, &faces(0.2), config(), nonlinear(), &cancelled,),
                Err(MoistureError::Cancelled)
            ));
        });
        assert_eq!(initial.moisture(), &[0.1; 8]);
        assert_eq!(initial.time_s(), 0.0);
    });
}
