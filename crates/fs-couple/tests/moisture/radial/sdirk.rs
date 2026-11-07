//! Independent radial semigroups and nonlinear ODEs for second-order time.

use super::*;

#[test]
fn radial_sdirk_unequal_shell_mode_has_second_order_decay_and_closed_water() {
    with_cx(|cx| {
        let radii: [f64; 3] = [0.02, 0.035, 0.1];
        let centers = [0.5 * (radii[0] + radii[1]), 0.5 * (radii[1] + radii[2])];
        let diffusion = 0.001;
        let duration = 0.25;
        let original = [0.22, 0.08];
        for symmetry in [
            RadialMoistureSymmetry::Cylinder { length_m: 0.3 },
            RadialMoistureSymmetry::Sphere,
        ] {
            // Independent physical volumes/resistance of two unequal shells.
            let (dimension, factor, resistance) = match symmetry {
                RadialMoistureSymmetry::Cylinder { length_m } => {
                    let factor = 2.0 * core::f64::consts::PI * length_m;
                    (2, factor, (centers[1] / centers[0]).ln() / factor)
                }
                RadialMoistureSymmetry::Sphere => {
                    let factor = 4.0 * core::f64::consts::PI;
                    (3, factor, (1.0 / centers[0] - 1.0 / centers[1]) / factor)
                }
            };
            let mass: [f64; 2] = core::array::from_fn(|i| {
                500.0 * factor * (radii[i + 1].powi(dimension) - radii[i].powi(dimension))
                    / f64::from(dimension)
            });
            let total = mass[0] + mass[1];
            let average = (mass[0] * original[0] + mass[1] * original[1]) / total;
            let rate = 500.0 * diffusion / resistance * (1.0 / mass[0] + 1.0 / mass[1]);
            let difference = (original[0] - original[1]) * (-rate * duration).exp();
            let exact = [
                average + mass[1] / total * difference,
                average - mass[0] / total * difference,
            ];
            let initial = specimen(
                symmetry,
                radii.to_vec(),
                original.to_vec(),
                [diffusion; 2],
                cx,
            );
            let errors = [4, 8, 16].map(|steps| {
                let mut state = initial.clone();
                for _ in 0..steps {
                    let (next, report) = state
                        .advance_sdirk2(
                            duration / f64::from(steps),
                            &[MoistureFace::Sealed, MoistureFace::Sealed],
                            config(),
                            nonlinear(),
                            cx,
                        )
                        .unwrap();
                    assert_eq!(report.endpoint.boundary_water_kg, [0.0; 2]);
                    assert_eq!(report.first_stage.boundary_water_kg, [0.0; 2]);
                    assert!(report.endpoint.nonlinear.unwrap().relative_residual < 1e-12);
                    near(next.water_mass_kg(), initial.water_mass_kg(), 2e-12);
                    assert!(next.moisture()[0] < state.moisture()[0]);
                    assert!(next.moisture()[1] > state.moisture()[1]);
                    state = next;
                }
                assert_eq!(state.time_s(), duration);
                state
                    .moisture()
                    .iter()
                    .zip(exact)
                    .map(|(&m, expected)| (m - expected).abs())
                    .fold(0.0, f64::max)
            });
            assert!(errors[0] / errors[1] > 3.8, "{symmetry:?}: {errors:?}");
            assert!(errors[1] / errors[2] > 3.8, "{symmetry:?}: {errors:?}");
            assert!(errors[1] / errors[2] < 4.5, "{symmetry:?}: {errors:?}");
            assert!(errors[2] < 9e-6, "{symmetry:?}: {errors:?}");
            assert_eq!(initial.moisture(), original);
            assert_eq!(initial.time_s(), 0.0);
        }
    });
}

#[test]
fn radial_sdirk_nonlinear_stages_use_actual_face_rates_and_weighted_water() {
    with_cx(|cx| {
        let radius: f64 = 0.1;
        let equilibrium: f64 = 0.2;
        let dt = 0.04;
        let gamma = 1.0 - core::f64::consts::FRAC_1_SQRT_2;
        let initial = specimen(
            RadialMoistureSymmetry::Sphere,
            vec![0.0, radius],
            vec![0.0],
            [0.0, 0.03],
            cx,
        );
        // D(M)=.1*M and center R/2 give M'=b*(Meq^2-M^2), b=.15/R^2.
        let b = 0.15 / radius.powi(2);
        let rate = |m: f64| b * (equilibrium.powi(2) - m * m);
        let root = |base: f64, h: f64| {
            let a = b * h;
            let c = base + a * equilibrium.powi(2);
            2.0 * c / (1.0 + (1.0 + 4.0 * a * c).sqrt())
        };
        let first = root(0.0, gamma * dt);
        let expected = root((1.0 - gamma) * dt * rate(first), gamma * dt);
        let boundaries = [
            MoistureFace::Sealed,
            MoistureFace::PrescribedDryBasis(equilibrium),
        ];
        let (next, report) = initial
            .advance_sdirk2(dt, &boundaries, config(), nonlinear(), cx)
            .unwrap();
        near(next.moisture()[0], expected, 3e-13);
        near(
            report.first_stage.stored_water_kg[1] / initial.dry_mass_kg(),
            first,
            3e-13,
        );
        near(
            report.first_stage.boundary_water_kg[1],
            initial.dry_mass_kg() * gamma * dt * rate(first),
            2e-12,
        );
        let total =
            initial.dry_mass_kg() * dt * ((1.0 - gamma) * rate(first) + gamma * rate(expected));
        assert_eq!(report.endpoint.boundary_water_kg[0], 0.0);
        near(report.endpoint.boundary_water_kg[1], total, 2e-12);
        near(report.endpoint.water_balance_residual_kg, 0.0, 2e-12);
        assert_eq!(report.endpoint.interval_s, [0.0, dt]);
        assert_eq!(report.first_stage.interval_s, [0.0, gamma * dt]);
        assert!(report.endpoint.kirchhoff_source.is_some());
        assert_eq!(initial.moisture(), &[0.0]);
        assert_eq!(initial.time_s(), 0.0);
        let replay = initial
            .advance_sdirk2(dt, &boundaries, config(), nonlinear(), cx)
            .unwrap();
        assert_eq!(next.moisture(), replay.0.moisture());
        assert_eq!(
            report.endpoint.boundary_water_kg,
            replay.1.endpoint.boundary_water_kg
        );
    });
}

#[test]
fn radial_sdirk_dry_front_converges_to_nonlinear_continuous_time_solution() {
    with_cx(|cx| {
        let initial = specimen(
            RadialMoistureSymmetry::Sphere,
            vec![0.0, 0.1],
            vec![0.0],
            [0.0, 0.03],
            cx,
        );
        let duration = 0.25;
        // Exact ODE from the original one-shell nonlinear radial resistance.
        let exact = 0.2 * (3.0_f64 * duration).tanh();
        let boundaries = [MoistureFace::Sealed, MoistureFace::PrescribedDryBasis(0.2)];
        let errors = [4, 8, 16].map(|steps| {
            let mut state = initial.clone();
            for _ in 0..steps {
                let (next, report) = state
                    .advance_sdirk2(
                        duration / f64::from(steps),
                        &boundaries,
                        config(),
                        nonlinear(),
                        cx,
                    )
                    .unwrap();
                assert!(next.moisture()[0] > state.moisture()[0] && next.moisture()[0] < 0.2);
                near(
                    report.endpoint.boundary_water_kg.iter().sum(),
                    next.water_mass_kg() - state.water_mass_kg(),
                    2e-12,
                );
                state = next;
            }
            (state.moisture()[0] - exact).abs()
        });
        assert!(errors[0] / errors[1] > 3.8, "{errors:?}");
        assert!(errors[1] / errors[2] > 3.8, "{errors:?}");
        assert!(errors[1] / errors[2] < 4.5, "{errors:?}");
        assert!(errors[2] < 2e-6, "{errors:?}");
    });
}

#[test]
fn radial_sdirk_finite_surface_exchange_uses_curved_area_and_final_surface_content() {
    with_cx(|cx| {
        let radius: f64 = 0.2;
        let h = 0.004;
        let diffusion = 0.01;
        let dt = 2.0;
        let equilibrium = 0.2;
        let initial = specimen(
            RadialMoistureSymmetry::Sphere,
            vec![0.0, radius],
            vec![0.1],
            [diffusion; 2],
            cx,
        );
        let boundaries = [
            MoistureFace::Sealed,
            MoistureFace::SourceSorption {
                law: surface::law([h; 2]),
                relative_humidity: 0.8,
            },
        ];
        let speed = diffusion / radius;
        let effective = 1.0 / (1.0 / h + 1.0 / speed);
        let z = dt * 3.0 * effective / radius;
        let gamma = 1.0 - core::f64::consts::FRAC_1_SQRT_2;
        let decay = (1.0 - (1.0 - 2.0 * gamma) * z) / (1.0 + gamma * z).powi(2);
        let expected = equilibrium + (0.1 - equilibrium) * decay;
        let (next, report) = initial
            .advance_sdirk2(dt, &boundaries, config(), nonlinear(), cx)
            .unwrap();
        near(next.moisture()[0], expected, 2e-13);
        near(
            report.endpoint.surface_moisture[1].unwrap(),
            (speed * expected + h * equilibrium) / (speed + h),
            2e-13,
        );
        near(
            report.endpoint.boundary_water_kg[1],
            initial.dry_mass_kg() * (expected - 0.1),
            2e-12,
        );
        assert!(report.endpoint.surface_sorption[1].is_some());
        assert!(report.endpoint.surface_transfer[1].is_some());
    });
}

#[test]
fn radial_sdirk_refuses_stiff_nonphysical_stages_without_clipping_or_method_fallback() {
    with_cx(|cx| {
        let initial = specimen(
            RadialMoistureSymmetry::Sphere,
            vec![0.0, 0.1],
            vec![0.2],
            [0.01; 2],
            cx,
        );
        let boundaries = [MoistureFace::Sealed, MoistureFace::PrescribedDryBasis(0.0)];
        let error = initial
            .advance_sdirk2(1.0, &boundaries, config(), nonlinear(), cx)
            .unwrap_err();
        assert!(
            matches!(error, MoistureError::StageRange { .. }),
            "{error:?}"
        );
        assert_eq!(initial.moisture(), &[0.2]);
        assert_eq!(initial.time_s(), 0.0);
        // Backward Euler has a different positive endpoint at this same dt.
        // An SDIRK failure must never publish that first-order substitute.
        let ordinary = initial
            .advance(1.0, &boundaries, config(), nonlinear(), cx)
            .unwrap();
        near(ordinary.0.moisture()[0], 0.05, 2e-13);
        let (next, _) = initial
            .advance_sdirk2(0.1, &boundaries, config(), nonlinear(), cx)
            .unwrap();
        assert!(next.moisture()[0] > 0.0 && next.moisture()[0] < 0.2);
    });
}

#[test]
fn radial_sdirk_dormant_reservoir_cannot_authorize_overshoot() {
    with_cx(|cx| {
        let initial = specimen(
            RadialMoistureSymmetry::Sphere,
            vec![0.05, 0.1],
            vec![0.1],
            [0.01; 2],
            cx,
        );
        let boundaries = [
            MoistureFace::SourceSorption {
                law: surface::law([0.0; 2]),
                relative_humidity: 1.0,
            },
            MoistureFace::PrescribedDryBasis(0.2),
        ];
        // Only the outer surface is active. z=3 overshoots .2 while staying
        // below the dormant inner reservoir's .25 equilibrium.
        let volume = 4.0 * core::f64::consts::PI / 3.0 * (0.1_f64.powi(3) - 0.05_f64.powi(3));
        let resistance = (1.0 / 0.075 - 1.0 / 0.1) / (4.0 * core::f64::consts::PI);
        let dt = 3.0 * resistance * volume / 0.01;
        let error = initial
            .advance_sdirk2(dt, &boundaries, config(), nonlinear(), cx)
            .unwrap_err();
        assert!(
            matches!(error, MoistureError::StageRange { bounds, .. } if bounds == [0.1, 0.2]),
            "{error:?}"
        );
        assert_eq!(initial.moisture(), &[0.1]);
        assert_eq!(initial.time_s(), 0.0);
    });
}
