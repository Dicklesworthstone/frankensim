use super::*;
use crate::bulk_viscous_heating::finite_window::*;

fn fine() -> BulkViscousWindowLimits {
    BulkViscousWindowLimits {
        maximum_panel_exponent: 0.2,
        maximum_diagonal_relative_error: 1e-11,
        ..BulkViscousWindowLimits::default()
    }
}
fn source<'a>(
    waves: &'a [BulkTractionAbsorptionResponse],
    initial: &DepositionSimulation<'_>,
    ends: &[f64],
    limits: BulkViscousWindowLimits,
) -> BulkViscousWindowHeating<'a> {
    BulkViscousWindowHeating::try_new(
        waves,
        "original/absolute-seconds",
        initial,
        &routes(ends),
        limits,
        1.0,
        || Ok(()),
    )
    .unwrap()
}

// Independent tensor Simpson integral of the real stress gradient. This does
// not use the new gradient field, the spatial Gauss rule, harmonic kernels,
// or the prepared quadratic spectrum. All tones interfere before squaring.
fn real_field_heat(
    waves: &[BulkTractionAbsorptionResponse],
    from: f64,
    to: f64,
    start: f64,
    duration: f64,
) -> f64 {
    let (nx, nt) = (192, 1024);
    let eta = 4.0 * 0.003 / 3.0 + 0.002;
    let dx = (to - from) / nx as f64;
    let dt = duration / nt as f64;
    let mut total = 0.0;
    for x in 0..=nx {
        let position = from + x as f64 * dx;
        let gradients: Vec<_> = waves
            .iter()
            .map(|wave| {
                let field = wave.at_distance(position).unwrap();
                (field.pressure_pa - field.normal_traction_pa) / C64::from_re(eta)
            })
            .collect();
        let mut time = 0.0;
        for t in 0..=nt {
            let seconds = start + t as f64 * dt;
            let strain: f64 = gradients
                .iter()
                .zip(waves)
                .map(|(g, wave)| {
                    let phase = wave.angular_frequency_rad_s() * seconds;
                    g.re * phase.cos() + g.im * phase.sin()
                })
                .sum();
            let weight = if t == 0 || t == nt {
                1.0
            } else if t % 2 == 0 {
                2.0
            } else {
                4.0
            };
            time += weight * strain * strain;
        }
        let weight = if x == 0 || x == nx {
            1.0
        } else if x % 2 == 0 {
            2.0
        } else {
            4.0
        };
        total += weight * time * dt / 3.0;
    }
    total * eta * area() * dx / 3.0
}

#[test]
fn finite_single_carrier_retains_the_actual_sum_frequency_heat() {
    let law = Caloric::new(8.0, 0);
    let waves = [response(0.0, 0.01, 1000.0, C64::new(2e-5, -1e-5), true)];
    let initial = initial(events(&[0.01], 0.0), &law);
    let source = source(&waves, &initial, &[0.01], fine());
    let cell = &source.cells()[0];
    let cycle = waves[0]
        .absorption_between(0.0, 0.01)
        .unwrap()
        .viscous_heat_w;
    near(cell.spectrum().dc_dd().to_f64(), cycle, 1e-12);
    assert_eq!(cell.spectrum().tones().len(), 1);
    let (start, duration) = (0.00031, 0.00017);
    let actual = cell.mean_heat_w(start, duration).unwrap() * duration;
    let independent = real_field_heat(&waves, 0.0, 0.01, start, duration);
    near(actual, independent, 2e-10);
    assert!((actual / (cycle * duration) - 1.0).abs() > 0.2);
    let field = waves[0].at_distance(0.004).unwrap();
    let eta = 4.0 * 0.003 / 3.0 + 0.002;
    let gradient = (field.pressure_pa - field.normal_traction_pa) / C64::from_re(eta);
    assert!((gradient - field.velocity_gradient_per_s).abs() < 2e-13 * gradient.abs());
}

#[test]
fn close_tone_beats_match_real_fields_and_are_not_independent_mean_powers() {
    let law = Caloric::new(8.0, 0);
    let waves = [
        response(0.0, 0.01, 1000.0, C64::from_re(2e-5), true),
        response(
            0.0,
            0.01,
            1030.0,
            C64::from_re(-2e-5 * 1000.0 / 1030.0),
            true,
        ),
    ];
    let initial = initial(events(&[0.01], 0.0), &law);
    let both = source(&waves, &initial, &[0.01], fine());
    let a = source(&waves[..1], &initial, &[0.01], fine());
    let b = source(&waves[1..], &initial, &[0.01], fine());
    for (start, maximum_ratio, minimum_ratio) in [(0.0, 0.02, 0.0), (0.1, 2.1, 1.8)] {
        let duration = 0.002;
        let coherent = both.cells()[0].mean_heat_w(start, duration).unwrap();
        let independent = a.cells()[0].mean_heat_w(start, duration).unwrap()
            + b.cells()[0].mean_heat_w(start, duration).unwrap();
        assert!(
            coherent < maximum_ratio * independent && coherent > minimum_ratio * independent,
            "coherent/independent = {}",
            coherent / independent
        );
        near(
            coherent * duration,
            real_field_heat(&waves, 0.0, 0.01, start, duration),
            3e-9,
        );
    }
    assert!(
        both.cells()[0]
            .spectrum()
            .tones()
            .iter()
            .any(|tone| tone.angular_frequency_dd().to_f64() == 30.0)
    );
}

#[test]
fn positive_spatial_quadrature_has_a_live_three_rung_convergence_axis() {
    let law = Caloric::new(8.0, 0);
    let waves = [response(0.0, 0.04, 1000.0, C64::new(2e-5, -1e-5), false)];
    let initial = initial(events(&[0.04], 0.0), &law);
    let exact = waves[0]
        .absorption_between(0.0, 0.04)
        .unwrap()
        .viscous_heat_w;
    let mut errors = Vec::new();
    let mut points = Vec::new();
    for exponent in [4.0, 2.0, 1.0] {
        let limits = BulkViscousWindowLimits {
            maximum_panel_exponent: exponent,
            maximum_diagonal_relative_error: 0.5,
            ..BulkViscousWindowLimits::default()
        };
        let source = source(&waves, &initial, &[0.04], limits);
        let cell = &source.cells()[0];
        let error = (cell.spectrum().dc_dd().to_f64() / exact - 1.0).abs();
        errors.push(error);
        points.push(cell.spatial_points);
        near(cell.maximum_diagonal_relative_error, error, 1e-6);
    }
    assert!(points[0] < points[1] && points[1] < points[2], "{points:?}");
    assert!(
        errors[0] > 1e-4 && errors[1] < errors[0] / 50.0 && errors[2] < errors[1] / 100.0,
        "positive Gauss refinement errors {errors:?}"
    );
    let strict = BulkViscousWindowLimits {
        maximum_panel_exponent: 4.0,
        maximum_diagonal_relative_error: 1e-10,
        ..BulkViscousWindowLimits::default()
    };
    assert!(matches!(
        BulkViscousWindowHeating::try_new(
            &waves,
            "clock",
            &initial,
            &routes(&[0.04]),
            strict,
            1.0,
            || Ok(())
        ),
        Err(BulkViscousWindowError::SpatialResolution { .. })
    ));
}

#[test]
fn equal_frequency_cancellation_is_exact_and_small_retained_fields_survive() {
    let law = Caloric::new(8.0, 0);
    let drive = C64::new(2e-5, -1e-5);
    let large = drive.scale(2.0_f64.powi(40));
    let waves = [
        response(0.0, 0.04, 1000.0, large, false),
        response(0.0, 0.04, 1000.0, -large, false),
    ];
    let initial = initial(events(&[0.04], 0.0), &law);
    let dark = source(&waves, &initial, &[0.04], fine());
    assert_eq!(dark.cells()[0].spectrum().dc_dd().to_f64(), 0.0);
    assert!(dark.cells()[0].spectrum().tones().is_empty());
    assert_eq!(dark.cells()[0].mean_heat_w(1e300, 0.01).unwrap(), 0.0);
    let step = dark
        .advance_to(&initial, 0.02, config(0.003), || Ok(()))
        .unwrap();
    assert_eq!(step.viscous_heat_j, 0.0);
    assert_eq!(step.thermal.energy_j(), 0.0);
    let small = response(0.0, 0.04, 1000.0, drive.scale(1e-9), false);
    let residual = [waves[0], small, waves[1]];
    let actual = source(&residual, &initial, &[0.04], fine());
    let only_small = [small];
    let expected = source(&only_small, &initial, &[0.04], fine());
    near(
        actual.cells()[0].mean_heat_w(0.003, 0.011).unwrap(),
        expected.cells()[0].mean_heat_w(0.003, 0.011).unwrap(),
        1e-10,
    );
}

#[test]
fn native_interval_energy_is_additive_and_uses_one_absolute_phase_clock() {
    let law = Caloric::new(8.0, 0);
    let frequencies = [300.0, 517.0];
    let drives = [C64::new(2e-5, -1e-5), C64::new(-1e-5, 1.7e-5)];
    let waves =
        std::array::from_fn::<_, 2, _>(|i| response(0.0, 0.04, frequencies[i], drives[i], false));
    let ends = [0.013, 0.04];
    let initial = initial(events(&ends, 0.0), &law);
    let source = source(&waves, &initial, &ends, fine());
    assert_eq!(source.source_clock(), "original/absolute-seconds");
    assert!(core::ptr::eq(source.responses().as_ptr(), waves.as_ptr()));
    let whole = source
        .advance_to(&initial, 0.02, config(0.02), || Ok(()))
        .unwrap();
    let first = source
        .advance_to(&initial, 0.013, config(0.003), || Ok(()))
        .unwrap();
    let second = source
        .advance_to(&first.thermal, 0.02, config(0.002), || Ok(()))
        .unwrap();
    near(
        first.viscous_heat_j + second.viscous_heat_j,
        whole.viscous_heat_j,
        2e-14,
    );
    near(second.thermal.energy_j(), whole.thermal.energy_j(), 2e-14);
    near(
        whole.heat_report.external_heat_j,
        whole.viscous_heat_j,
        2e-15,
    );
    near(
        whole.heat_report.stored_energy_change_j,
        whole.viscous_heat_j,
        2e-15,
    );
    near(
        whole.viscous_heat_j,
        real_field_heat(&waves, 0.0, 0.04, 0.0, 0.02),
        3e-9,
    );
    for cell in whole.thermal.cells().values() {
        assert_eq!(cell.boundary.source_w, 0.0);
    }
    let shift = 0.017;
    let shifted = std::array::from_fn::<_, 2, _>(|i| {
        let phase = frequencies[i] * shift;
        response(
            0.0,
            0.04,
            frequencies[i],
            drives[i] * C64::new(phase.cos(), -phase.sin()),
            false,
        )
    });
    let shifted_source = self::source(&shifted, &initial, &ends, fine());
    for (a, b) in source.cells().iter().zip(shifted_source.cells()) {
        near(
            a.mean_heat_w(shift, 0.009).unwrap(),
            b.mean_heat_w(0.0, 0.009).unwrap(),
            3e-14,
        );
    }
}

#[test]
fn native_events_and_conduction_receive_the_actual_window_heat_once() {
    let law = Caloric::new(8.0, 0);
    let waves = [
        response(0.0, 0.04, 1000.0, C64::from_re(2e-5), false),
        response(0.0, 0.04, 1030.0, C64::from_re(-1.5e-5), false),
    ];
    let ends = [0.02, 0.04];
    let mut schedule = events(&ends, 0.01);
    schedule.push(TimedEvent {
        time_s: 0.013,
        event: Event::Boundary(id(0), boundary(2e-5)),
    });
    let initial = initial(schedule, &law);
    let source = source(&waves, &initial, &ends, fine());
    let step = source
        .advance_to(&initial, 0.04, config(0.03), || Ok(()))
        .unwrap();
    let mut exact_interval = 0.0;
    for (cell, (_, actual_j)) in source.cells().iter().zip(&step.cell_energy_j) {
        let expected = cell.mean_heat_w(0.0, 0.04).unwrap() * 0.04;
        near(*actual_j, expected, 4e-14);
        exact_interval += expected;
    }
    let authored = 2e-5 * (0.04 - 0.013);
    near(step.viscous_heat_j, exact_interval, 3e-14);
    near(
        step.heat_report.external_heat_j,
        exact_interval + authored,
        3e-14,
    );
    near(step.thermal.energy_j(), exact_interval + authored, 3e-14);
    assert_eq!(step.heat_report.events, 1);
    assert_eq!(step.thermal.links()[&1].conductance_w_k, 0.01);
    assert_eq!(step.thermal.cells()[&id(0)].boundary.source_w, 2e-5);
}

#[test]
fn foreign_physics_resolution_budgets_and_late_refusals_are_transactional() {
    let law = Caloric::new(8.0, 0);
    let base = response(0.0, 0.04, 1000.0, C64::from_re(2e-5), false);
    let initial = initial(events(&[0.04], 0.0), &law);
    for other in [
        response(0.08, 0.04, 1030.0, C64::from_re(2e-5), false),
        response(0.0, 0.03, 1030.0, C64::from_re(2e-5), false),
        response(0.0, 0.04, 1000.0, C64::from_re(2e-5), true),
    ] {
        assert!(
            BulkViscousWindowHeating::try_new(
                &[base, other],
                "clock",
                &initial,
                &routes(&[0.04]),
                fine(),
                1.0,
                || Ok(())
            )
            .is_err()
        );
    }
    let waves = [base];
    let exhausted = BulkViscousWindowLimits {
        maximum_spatial_points: 1,
        ..fine()
    };
    assert!(
        BulkViscousWindowHeating::try_new(
            &waves,
            "clock",
            &initial,
            &routes(&[0.04]),
            exhausted,
            1.0,
            || Ok(())
        )
        .is_err()
    );
    let source = source(&waves, &initial, &[0.04], fine());
    let mut calls = 0;
    let refusal = source.advance_to(&initial, 0.02, config(0.002), || {
        calls += 1;
        if calls == 25 {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    });
    assert!(matches!(
        refusal,
        Err(BulkViscousWindowError::Heating(
            BulkViscousHeatingError::Thermal(DepositionError::Cancelled)
        ))
    ));
    let guarded = BulkViscousWindowHeating::try_new(
        &waves,
        "clock",
        &initial,
        &routes(&[0.04]),
        fine(),
        1e-8,
        || Ok(()),
    )
    .unwrap();
    assert!(matches!(
        guarded.advance_to(&initial, 0.02, config(0.002), || Ok(())),
        Err(BulkViscousWindowError::Heating(
            BulkViscousHeatingError::TemperatureExcursion { .. }
        ))
    ));
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.energy_j(), 0.0);
}
