//! Independent time/refinement and exact schedule tests on the native film.
use super::*;
use crate::film::evolution::adaptive::*;
use std::convert::Infallible;

fn time_controls(transport: f64, energy: f64) -> ThinFilmAdaptiveControls {
    ThinFilmAdaptiveControls {
        initial_dt_s: 1.0,
        minimum_dt_s: 1e-7,
        maximum_dt_s: 1.0,
        maximum_trials: 16_384,
        maximum_accepted_steps: 32_768,
        height_absolute_tolerance_m: 1e-7,
        transport_absolute_tolerance_m2: transport,
        energy_absolute_tolerance_j_per_m: energy,
        relative_tolerance: 0.0,
    }
}
fn open_drive(acceleration: f64, traction: f64, temperature: Option<f64>) -> ThinFilmDrive {
    let end = ThinFilmOpenEnd {
        pressure_pa: 0.3,
        surface_shear_pa: traction,
    };
    ThinFilmDrive {
        acceleration_m_s2: acceleration,
        surface_shear_pa: vec![traction; 3],
        open_ends: Some(ThinFilmOpenEnds {
            left: end,
            right: end,
        }),
        face_temperatures_k: temperature.map(|t| vec![t; 5]),
    }
}
fn open_newtonian() -> ThinFilmLine {
    ThinFilmLine::newtonian(
        material(0.04),
        geometry(ThinFilmBoundary::Open, 1.5),
        vec![0.1; 4],
        0.0,
    )
    .unwrap()
}
fn open_carreau() -> ThinFilmLine {
    ThinFilmLine::carreau(
        source(0.4),
        &[300.0; 5],
        2.0,
        0.04,
        20.0,
        geometry(ThinFilmBoundary::Open, 1.5),
        vec![0.1; 4],
        0.0,
    )
    .unwrap()
}

#[test]
fn g1_uniform_height_still_refines_actual_throughflow_and_counterflow_heat() {
    // F=0.4+0.6t; direct time integrals, with no reference to reports.
    let integral_force = 0.4 + 0.6 / 2.0;
    let integral_force_squared = 0.4_f64.powi(2) + 0.4 * 0.6 + 0.6_f64.powi(2) / 3.0;
    let (h, mu, length) = (0.1_f64, 0.5, 2.0);
    for traction_ratio in [0.0, -2.0 / 3.0, -0.4] {
        // The middle case has exactly zero net Q but nonzero counterflow loss.
        // The last case has opposing body/surface work; compare both channels.
        let volume = h.powi(3) / mu * (1.0 / 3.0 + traction_ratio / 2.0) * integral_force;
        let heat = length * h.powi(3) / mu
            * (1.0 / 3.0 + traction_ratio + traction_ratio * traction_ratio)
            * integral_force_squared;
        let mut errors = Vec::new();
        for allowance in [2e-7, 2e-9] {
            let mut line = open_newtonian();
            let report = line
                .advance_adaptive(
                    1.0,
                    &[],
                    limits(),
                    controls(),
                    time_controls(allowance, allowance),
                    |query| {
                        let acceleration = 0.2 + 0.3 * query.time_s;
                        Ok::<_, Infallible>(open_drive(
                            acceleration,
                            traction_ratio * 2.0 * acceleration * h,
                            None,
                        ))
                    },
                    || false,
                )
                .unwrap();
            assert!(report.accuracy_rejections > 0);
            assert!(report.accepted_comparisons > 1);
            assert!(report.maximum_accepted_normalized_error <= 1.0);
            assert_eq!(report.accepted_steps, 2 * report.accepted_comparisons);
            for &height in line.heights_m() {
                near(height, h, 3e-12);
            }
            for face in &report.face_integrals {
                // All physical profiles are uniform, including both half-cells.
                near(
                    face.transported_volume_m2,
                    report.face_integrals[0].transported_volume_m2,
                    2e-13,
                );
            }
            if traction_ratio == -2.0 / 3.0 {
                near(report.face_integrals[0].transported_volume_m2, 0.0, 2e-13);
                assert!(report.viscous_heat_j_per_m > 0.0);
            }
            if traction_ratio == -0.4 {
                assert!(report.body_work_j_per_m > 0.0);
                assert!(report.surface_work_j_per_m < 0.0);
            }
            near(report.volume_change_m2, 0.0, 1e-11);
            near(report.volume_balance_residual_m2, 0.0, 1e-11);
            near(report.energy_residual_j_per_m, 0.0, 1e-11);
            errors.push((
                (report.face_integrals[0].transported_volume_m2 - volume).abs(),
                (report.viscous_heat_j_per_m - heat).abs(),
            ));
        }
        if traction_ratio != -2.0 / 3.0 {
            assert!(errors[1].0 < 0.25 * errors[0].0, "{errors:?}");
        }
        assert!(errors[1].1 < 0.25 * errors[0].1, "{errors:?}");
    }
}

#[test]
fn g1_adaptive_height_decay_converges_to_the_independent_discrete_fourier_reference() {
    let (n, h, amplitude, dx) = (12, 0.1, 1e-6, 0.5);
    let shape: Vec<_> = (0..n)
        .map(|i| (std::f64::consts::TAU * i as f64 / n as f64).cos())
        .collect();
    let initial: Vec<_> = shape.iter().map(|&x| h + amplitude * x).collect();
    let k2 = 4.0 * (std::f64::consts::PI / n as f64).sin().powi(2) / (dx * dx);
    let rate = h * h * h / (3.0 * 0.5) * (3.0 * k2 + 0.04 * k2 * k2);
    let target = 0.5 / rate;
    let expected = amplitude * (-0.5_f64).exp();
    let mut errors = Vec::new();
    for tolerance in [2e-8, 2e-10] {
        let mut line = ThinFilmLine::newtonian(
            material(0.04),
            geometry(ThinFilmBoundary::Periodic, 1.5),
            initial.clone(),
            0.0,
        )
        .unwrap();
        let report = line
            .advance_adaptive(
                target,
                &[],
                limits(),
                controls(),
                ThinFilmAdaptiveControls {
                    initial_dt_s: target,
                    minimum_dt_s: target / 1e6,
                    maximum_dt_s: target,
                    height_absolute_tolerance_m: tolerance,
                    // This case independently isolates the evolving-height estimate.
                    transport_absolute_tolerance_m2: 1.0,
                    energy_absolute_tolerance_j_per_m: 1.0,
                    ..time_controls(1.0, 1.0)
                },
                |_| {
                    Ok::<_, Infallible>(ThinFilmDrive {
                        acceleration_m_s2: 0.0,
                        surface_shear_pa: vec![0.0; n],
                        open_ends: None,
                        face_temperatures_k: None,
                    })
                },
                || false,
            )
            .unwrap();
        let modal = 2.0 / n as f64
            * line
                .heights_m()
                .iter()
                .zip(&shape)
                .map(|(&v, &s)| (v - h) * s)
                .sum::<f64>();
        errors.push((modal - expected).abs());
        assert!(report.stored_energy_change_j_per_m < 0.0);
        assert!(report.viscous_heat_j_per_m > 0.0);
        assert!(report.backward_euler_energy_defect_j_per_m > 0.0);
        near(report.volume_balance_residual_m2, 0.0, 1e-11);
        near(report.energy_residual_j_per_m, 0.0, 1e-11);
    }
    assert!(errors[1] < 0.3 * errors[0], "{errors:?}");
    assert!(errors[1] < 0.005 * amplitude, "{errors:?}");
}

#[test]
fn g3_decimal_events_keep_exact_left_limits_and_original_carreau_sources_and_accounts() {
    let events = [0.037, 0.074, 0.111];
    let temperatures = [300.0, 320.0, 340.0];
    let acceleration = [0.2, -0.15, 0.3];
    let mut line = open_carreau();
    let mut queries = Vec::new();
    let report = line
        .advance_adaptive(
            0.111,
            &events,
            limits(),
            controls(),
            ThinFilmAdaptiveControls {
                initial_dt_s: 0.037,
                maximum_dt_s: 0.037,
                ..time_controls(1e-10, 1e-10)
            },
            |query| {
                queries.push(query);
                assert!(query.segment < 3);
                let lower = if query.segment == 0 {
                    0.0
                } else {
                    events[query.segment - 1]
                };
                assert!(query.time_s > lower && query.time_s <= events[query.segment]);
                Ok::<_, Infallible>(open_drive(
                    acceleration[query.segment],
                    0.0,
                    Some(temperatures[query.segment]),
                ))
            },
            || false,
        )
        .unwrap();
    assert_eq!(line.accepted_time_s(), 0.111);
    assert_eq!(report.events_reached, 3);
    assert_eq!(report.accepted_comparisons, 3);
    assert_eq!(report.accuracy_rejections, 0);
    for (segment, &event) in events.iter().enumerate() {
        assert!(queries.contains(&ThinFilmDriveQuery {
            time_s: event,
            segment
        }));
    }
    for i in 0..5 {
        assert_eq!(line.face_temperature_k(i), Some(340.0));
    }
    let mut volume = 0.0;
    let mut heat = 0.0;
    let mut start = 0.0;
    for segment in 0..3 {
        let response = original_profile(
            source(0.4),
            temperatures[segment],
            0.1,
            2.0 * acceleration[segment],
            0.0,
        );
        let dt = events[segment] - start;
        volume += dt * response[0];
        heat += dt * 2.0 * response[2];
        start = events[segment];
    }
    near(
        report.face_integrals[0].transported_volume_m2,
        volume,
        2e-11,
    );
    near(report.viscous_heat_j_per_m, heat, 2e-11);
    near(report.energy_residual_j_per_m, 0.0, 1e-11);
    // Comparison trajectories must not appear in accepted heat/volume.
    assert!(report.endpoint_evaluations > report.accepted_steps);
}

#[test]
fn g4_late_source_work_and_cancellation_refusals_preserve_all_heights_temperatures_and_time() {
    let target = 0.111;
    let events = [0.037, 0.074, 0.111];
    let adaptive = ThinFilmAdaptiveControls {
        initial_dt_s: 0.037,
        maximum_dt_s: 0.037,
        ..time_controls(1e-8, 1e-8)
    };
    let drive = |query: ThinFilmDriveQuery| {
        Ok::<_, &'static str>(open_drive(0.2, 0.0, Some(300.0 + 20.0 * query.time_s)))
    };
    let mut reference = open_carreau();
    let mut polls = 0;
    let expected = reference
        .advance_adaptive(
            target,
            &events,
            limits(),
            controls(),
            adaptive,
            drive,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let assert_opening = |line: &ThinFilmLine| {
        assert_eq!(line.heights_m(), &[0.1; 4]);
        assert_eq!(line.accepted_time_s(), 0.0);
        for i in 0..5 {
            assert_eq!(line.face_temperature_k(i), Some(300.0));
        }
    };
    let mut victim = open_carreau();
    assert!(
        victim
            .set_face_temperatures(&[300.0, 300.0, 300.0, 300.0, 361.0], || false)
            .is_err()
    );
    assert_opening(&victim);
    assert!(matches!(
        victim.advance_adaptive(
            target,
            &events,
            limits(),
            controls(),
            adaptive,
            |query| if query.segment == 2 {
                Err("original temperature source refused")
            } else {
                drive(query)
            },
            || false,
        ),
        Err(ThinFilmAdaptiveError::Source(
            "original temperature source refused"
        ))
    ));
    assert_opening(&victim);
    assert!(matches!(
        victim.advance_adaptive(
            target,
            &events,
            limits(),
            controls(),
            adaptive,
            |query| {
                let mut prescribed = drive(query)?;
                if query.segment == 2 {
                    prescribed.face_temperatures_k.as_mut().unwrap()[4] = 361.0;
                }
                Ok::<_, &'static str>(prescribed)
            },
            || false,
        ),
        Err(ThinFilmAdaptiveError::Film(ThinFilmError::Material(_)))
    ));
    assert_opening(&victim);
    assert!(matches!(
        victim.advance_adaptive(
            target,
            &events,
            limits(),
            ThinFilmControls {
                maximum_endpoint_evaluations: expected.endpoint_evaluations - 1,
                ..controls()
            },
            adaptive,
            drive,
            || false,
        ),
        Err(ThinFilmAdaptiveError::Film(ThinFilmError::EndpointBudget))
    ));
    assert_opening(&victim);
    let mut calls = 0;
    assert!(matches!(
        victim.advance_adaptive(
            target,
            &events,
            limits(),
            controls(),
            adaptive,
            drive,
            || {
                calls += 1;
                calls == polls
            },
        ),
        Err(ThinFilmAdaptiveError::Film(ThinFilmError::Cancelled))
    ));
    assert_eq!(calls, polls);
    assert_opening(&victim);
    let actual = victim
        .advance_adaptive(
            target,
            &events,
            limits(),
            controls(),
            adaptive,
            drive,
            || false,
        )
        .unwrap();
    assert_eq!(victim.heights_m(), reference.heights_m());
    assert_eq!(actual.face_integrals, expected.face_integrals);
    assert_eq!(actual.comparisons, expected.comparisons);
    assert_eq!(
        actual.viscous_heat_j_per_m.to_bits(),
        expected.viscous_heat_j_per_m.to_bits()
    );
    for i in 0..5 {
        assert_eq!(
            victim.face_temperature_k(i),
            reference.face_temperature_k(i)
        );
    }

    let mut unresolved = open_newtonian();
    let strict = ThinFilmAdaptiveControls {
        minimum_dt_s: 1.0,
        initial_dt_s: 1.0,
        maximum_dt_s: 1.0,
        ..time_controls(1e-12, 1e-12)
    };
    assert!(matches!(
        unresolved.advance_adaptive(
            1.0,
            &[],
            limits(),
            controls(),
            strict,
            |query| Ok::<_, Infallible>(open_drive(0.2 + 0.3 * query.time_s, 0.0, None)),
            || false,
        ),
        Err(ThinFilmAdaptiveError::TimeAccuracy { .. })
    ));
    assert_eq!(unresolved.heights_m(), &[0.1; 4]);
    assert_eq!(unresolved.accepted_time_s(), 0.0);
}
