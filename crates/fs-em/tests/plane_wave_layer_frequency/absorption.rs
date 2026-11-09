//! Independent carrier, material-heat and Airy oracles for absorption tangents.
use super::*;

fn near_real(actual: f64, expected: f64, tolerance: f64) {
    near(C64::from_re(actual), C64::from_re(expected), tolerance);
}

#[test]
fn coherent_film_heat_slopes_follow_native_fields_and_actual_incident_power() {
    let ambient = material(2.0 * E, M);
    let exit_material = debye(1.4 * E);
    let exit = DielectricResponse::Debye(&exit_material);
    let first = debye(3.0 * E);
    let second = debye(7.0 * E);
    let films = [
        DielectricLayer {
            response: DielectricResponse::Debye(&first),
            thickness_m: 0.14,
        },
        DielectricLayer {
            response: DielectricResponse::Debye(&second),
            thickness_m: 0.0,
        },
        DielectricLayer {
            response: DielectricResponse::Debye(&second),
            thickness_m: 0.09,
        },
    ];
    let omega = 0.7 * W;
    let cosine = 0.73;
    let step = 1e-4 * omega;
    let normal_admittance = det::sqrt(2.0 * E / M);
    let stack = PlanarDielectricLayers::try_new(&ambient, exit, omega, &films).unwrap();
    for polarization in [Polarization::Te, Polarization::Tm] {
        for constraint in [
            Constraint::FixedAngle,
            Constraint::FixedTransverseWavenumber,
        ] {
            let jet = stack
                .frequency_absorption(cosine, polarization, constraint)
                .unwrap();
            assert_eq!(jet.layer_count(), films.len());
            // Every frequency uses the ORIGINAL source law, native interior
            // field owner and actual incidence angle. No jet enters this oracle.
            // Each pair is [absorbed fraction, physical heat for unit peak E].
            let samples: [[[f64; 2]; 3]; 5] = core::array::from_fn(|sample_index| {
                let frequency = omega + (sample_index as f64 - 2.0) * step;
                let angle = match constraint {
                    Constraint::FixedAngle => cosine,
                    Constraint::FixedTransverseWavenumber => {
                        (1.0 - (1.0 - cosine * cosine) * (omega / frequency).powi(2)).sqrt()
                    }
                };
                let incident_power = 0.5
                    * normal_admittance
                    * match polarization {
                        Polarization::Te => angle,
                        Polarization::Tm => 1.0 / angle,
                    };
                let native =
                    PlanarDielectricLayers::try_new(&ambient, exit, frequency, &films).unwrap();
                core::array::from_fn(|layer_index| {
                    let heat = native
                        .prepare_layer_field(layer_index, angle, polarization)
                        .unwrap()
                        .absorbed_power_w_m2();
                    [heat / incident_power, heat]
                })
            });
            let incident_power = 0.5
                * normal_admittance
                * match polarization {
                    Polarization::Te => cosine,
                    Polarization::Tm => 1.0 / cosine,
                };
            let mut fraction_sum = 0.0;
            let mut slope_sum = 0.0;
            let mut heat_sum = 0.0;
            let mut heat_slope_sum = 0.0;
            for layer_index in 0..films.len() {
                let fraction = jet.layer_absorption(layer_index).unwrap();
                let heat = jet.layer_power_density(layer_index, C64::ONE).unwrap();
                let slope: [f64; 2] = core::array::from_fn(|component| {
                    (samples[0][layer_index][component] - 8.0 * samples[1][layer_index][component]
                        + 8.0 * samples[3][layer_index][component]
                        - samples[4][layer_index][component])
                        * omega
                        / (12.0 * step)
                });
                near_real(fraction.fraction(), samples[2][layer_index][0], 3e-12);
                near_real(fraction.scaled_derivative(), slope[0], 3e-8);
                near_real(
                    heat.power_w_m2(),
                    samples[2][layer_index][1],
                    3e-12 * incident_power,
                );
                near_real(
                    heat.scaled_derivative_w_m2(),
                    slope[1],
                    3e-8 * incident_power,
                );
                near_real(
                    heat.derivative_w_s_m2() * omega,
                    heat.scaled_derivative_w_m2(),
                    1e-15 * incident_power,
                );
                if layer_index == 1 {
                    assert_eq!(fraction.fraction(), 0.0);
                    assert_eq!(fraction.scaled_derivative(), 0.0);
                    assert_eq!(heat.power_w_m2(), 0.0);
                    assert_eq!(heat.scaled_derivative_w_m2(), 0.0);
                } else {
                    assert!(fraction.fraction() > 0.0);
                    assert!(heat.power_w_m2() > 0.0);
                }
                fraction_sum += fraction.fraction();
                slope_sum += fraction.scaled_derivative();
                heat_sum += heat.power_w_m2();
                heat_slope_sum += heat.scaled_derivative_w_m2();
            }
            let total = jet.total_absorption();
            let heat = jet.total_power_density(C64::ONE).unwrap();
            near_real(total.fraction(), fraction_sum, 3e-13);
            near_real(total.scaled_derivative(), slope_sum, 3e-13);
            near_real(heat.power_w_m2(), heat_sum, 3e-13 * incident_power);
            near_real(
                heat.scaled_derivative_w_m2(),
                heat_slope_sum,
                3e-13 * incident_power,
            );
            assert_eq!(
                total.fraction(),
                jet.response().carrier().absorbed_power_fraction
            );
            let native_flux_slopes = power_slope(
                &ambient,
                exit,
                &films,
                omega,
                cosine,
                polarization,
                constraint,
            );
            near_real(
                native_flux_slopes[0] + native_flux_slopes[1] + total.scaled_derivative(),
                0.0,
                6e-8,
            );
            assert!(jet.power_balance_relative_residual() < 4096.0 * f64::EPSILON);

            let incident_log_slope = match constraint {
                Constraint::FixedAngle => 0.0,
                Constraint::FixedTransverseWavenumber => {
                    let magnitude = (1.0 - cosine * cosine) / (cosine * cosine);
                    if polarization == Polarization::Te {
                        magnitude
                    } else {
                        -magnitude
                    }
                }
            };
            near_real(
                heat.scaled_derivative_w_m2(),
                incident_power
                    * (total.scaled_derivative() + incident_log_slope * total.fraction()),
                3e-13 * incident_power,
            );
            if constraint == Constraint::FixedTransverseWavenumber {
                assert!(
                    (heat.scaled_derivative_w_m2() - incident_power * total.scaled_derivative())
                        .abs()
                        > 1e-5 * incident_power
                );
            }
            let scaled = jet.total_power_density(C64::new(2.0, -3.0)).unwrap();
            near_real(
                scaled.power_w_m2(),
                13.0 * heat.power_w_m2(),
                3e-13 * incident_power,
            );
            near_real(
                scaled.scaled_derivative_w_m2(),
                13.0 * heat.scaled_derivative_w_m2(),
                3e-13 * incident_power,
            );
            let zero = jet.total_power_density(C64::ZERO).unwrap();
            assert_eq!(zero.power_w_m2(), 0.0);
            assert_eq!(zero.scaled_derivative_w_m2(), 0.0);
            assert!(jet.total_power_density(C64::new(f64::NAN, 0.0)).is_err());
            assert!(jet.layer_absorption(3).is_none());
            assert!(jet.layer_power_density(3, C64::ONE).is_err());
        }
    }
}

#[test]
fn positive_sub_epsilon_heat_critical_zero_and_opaque_limits_remain_distinct() {
    let ambient = material(E, M);
    let dielectric = debye(2.0 * E);
    let films = [
        DielectricLayer {
            response: DielectricResponse::Debye(&dielectric),
            thickness_m: 1e-20,
        },
        DielectricLayer {
            response: DielectricResponse::Debye(&dielectric),
            thickness_m: 0.0,
        },
    ];
    let stack = PlanarDielectricLayers::try_new(
        &ambient,
        DielectricResponse::Lossless(&ambient),
        W,
        &films,
    )
    .unwrap();
    let jet = stack
        .frequency_absorption(1.0, Polarization::Te, Constraint::FixedAngle)
        .unwrap();
    // At omega*tau=1, -Im(epsilon)=2E and D[-Im(epsilon)]=0.
    // The complete thin-film DA therefore equals A to leading order.
    let expected = W * (2.0 * E) * 1e-20 / det::sqrt(E / M);
    let absorbed = jet.layer_absorption(0).unwrap();
    assert!(absorbed.fraction() > 0.0 && absorbed.fraction() < 1e-16);
    near_real(absorbed.fraction(), expected, 5e-13 * expected);
    near_real(absorbed.scaled_derivative(), expected, 5e-13 * expected);
    let heat = jet.layer_power_density(0, C64::ONE).unwrap();
    let expected_heat = 0.5 * W * (2.0 * E) * 1e-20;
    near_real(heat.power_w_m2(), expected_heat, 5e-13 * expected_heat);
    near_real(
        heat.scaled_derivative_w_m2(),
        expected_heat,
        5e-13 * expected_heat,
    );
    near_real(
        heat.power_w_m2(),
        stack
            .prepare_layer_field(0, 1.0, Polarization::Te)
            .unwrap()
            .absorbed_power_w_m2(),
        5e-13 * expected_heat,
    );
    assert_eq!(jet.layer_absorption(1).unwrap().fraction(), 0.0);
    assert_eq!(jet.layer_absorption(1).unwrap().scaled_derivative(), 0.0);

    let critical_ambient = material(4.0 * E, M);
    let critical = material(3.0 * E, M);
    let critical_films = [
        DielectricLayer {
            response: DielectricResponse::Lossless(&critical),
            thickness_m: 0.2,
        },
        DielectricLayer {
            response: DielectricResponse::Debye(&dielectric),
            thickness_m: 0.0,
        },
    ];
    let critical_stack = PlanarDielectricLayers::try_new(
        &critical_ambient,
        DielectricResponse::Lossless(&critical_ambient),
        W,
        &critical_films,
    )
    .unwrap();
    for polarization in [Polarization::Te, Polarization::Tm] {
        for constraint in [
            Constraint::FixedAngle,
            Constraint::FixedTransverseWavenumber,
        ] {
            let lossless = critical_stack
                .frequency_absorption(0.5, polarization, constraint)
                .unwrap();
            for layer_index in 0..2 {
                let fraction = lossless.layer_absorption(layer_index).unwrap();
                assert_eq!(fraction.fraction(), 0.0);
                assert_eq!(fraction.scaled_derivative(), 0.0);
            }
            assert_eq!(lossless.total_absorption().fraction(), 0.0);
            assert_eq!(lossless.total_absorption().scaled_derivative(), 0.0);
        }
    }
    let opaque_films = [DielectricLayer {
        response: DielectricResponse::Debye(&dielectric),
        thickness_m: 1e6,
    }];
    let opaque = PlanarDielectricLayers::try_new(
        &ambient,
        DielectricResponse::Lossless(&ambient),
        W,
        &opaque_films,
    )
    .unwrap();
    let native = opaque.response(0.8, Polarization::Tm).unwrap();
    assert!(native.absorbed_power_fraction > 0.0);
    assert_eq!(native.transmitted_tangential_e, C64::ZERO);
    assert!(
        opaque
            .frequency_absorption(0.8, Polarization::Tm, Constraint::FixedAngle)
            .is_err()
    );
}

// Independent symmetric single-film Airy law. Subtracting R and T is safe for
// this order-one TEST oracle; it is deliberately not the product loss owner.
fn airy_absorption(omega: f64, thickness: f64) -> f64 {
    let epsilon = C64::from_re(64.0 * E) + C64::from_re(4.0 * E) / C64::new(1.0, omega / W);
    let admittance_ratio = epsilon.scale(1.0 / E).sqrt();
    let phase = epsilon.scale(M).sqrt().scale(omega * thickness);
    let (sine, cosine) = phase.re.sin_cos();
    let propagator = C64::new(cosine, -sine).scale(phase.im.exp());
    let front = (C64::ONE - admittance_ratio) / (C64::ONE + admittance_ratio);
    let denominator = C64::ONE - front * front * propagator * propagator;
    let reflected = front * (C64::ONE - propagator * propagator) / denominator;
    let transmitted = (C64::ONE - front * front) * propagator / denominator;
    1.0 - reflected.norm_sq() - transmitted.norm_sq()
}

#[test]
fn dispersive_absorption_peak_and_half_height_match_independent_airy_law() {
    let ambient = material(E, M);
    let dielectric = debye(64.0 * E);
    let thickness = core::f64::consts::PI / (W * det::sqrt(66.0 * E * M));
    let films = [DielectricLayer {
        response: DielectricResponse::Debye(&dielectric),
        thickness_m: thickness,
    }];
    let stack = PlanarDielectricLayers::try_new(
        &ambient,
        DielectricResponse::Lossless(&ambient),
        W,
        &films,
    )
    .unwrap();
    let interval = (0.6 * W, 1.4 * W);
    let (mut left, mut right) = interval;
    for _ in 0..90 {
        let first = left + (right - left) / 3.0;
        let second = right - (right - left) / 3.0;
        if airy_absorption(first, thickness) < airy_absorption(second, thickness) {
            left = first;
        } else {
            right = second;
        }
    }
    let expected_peak = (left + right) * 0.5;
    let expected_power = airy_absorption(expected_peak, thickness);
    let half_height = expected_power * 0.5;
    let crossing = |mut a: f64, mut b: f64, rising: bool| {
        for _ in 0..80 {
            let middle = (a + b) * 0.5;
            if (airy_absorption(middle, thickness) < half_height) == rising {
                a = middle;
            } else {
                b = middle;
            }
        }
        (a + b) * 0.5
    };
    let lower = crossing(interval.0, expected_peak, true);
    let upper = crossing(expected_peak, interval.1, false);
    let incident_power = 0.5 * det::sqrt(E / M);
    let mut reference = None;
    for (observation, scale) in [
        (Observation::AbsorbedPowerFraction, 1.0),
        (
            Observation::LayerAbsorbedPowerFraction { layer_index: 0 },
            1.0,
        ),
        (Observation::AbsorbedPowerDensity, incident_power),
        (
            Observation::LayerAbsorbedPowerDensity { layer_index: 0 },
            incident_power,
        ),
    ] {
        let actual = stack
            .driven_resonance(
                1.0,
                Polarization::Te,
                Constraint::FixedAngle,
                observation,
                interval,
                controls(true),
                || false,
            )
            .unwrap();
        // A value-only maximum loses first-order resolution at a smooth peak.
        // This oracle tolerance is independent of the much tighter solver bracket.
        near_real(actual.peak.angular_frequency_rad_s, expected_peak, 20.0);
        assert!(actual.peak.bracket_rad_s.1 - actual.peak.bracket_rad_s.0 <= 0.1);
        near_real(
            actual.observation.norm_sq(),
            expected_power * scale,
            3e-12 * scale,
        );
        assert_eq!(actual.observation.im, 0.0);
        let band = actual.half_power.unwrap();
        near_real(band.lower.angular_frequency_rad_s, lower, 1.0);
        near_real(band.upper.angular_frequency_rad_s, upper, 1.0);
        near_real(band.bandwidth_rad_s, upper - lower, 1.5);
        near_real(
            band.response_quality_factor,
            expected_peak / (upper - lower),
            2e-7,
        );
        assert!(actual.evaluations <= controls(true).maximum_evaluations);
        if let Some((previous_peak, previous_width)) = reference {
            near_real(actual.peak.angular_frequency_rad_s, previous_peak, 0.3);
            near_real(band.bandwidth_rad_s, previous_width, 0.5);
        } else {
            reference = Some((actual.peak.angular_frequency_rad_s, band.bandwidth_rad_s));
        }
    }
    assert!(matches!(
        stack.driven_resonance(
            1.0,
            Polarization::Te,
            Constraint::FixedAngle,
            Observation::AbsorbedPowerFraction,
            interval,
            controls(false),
            || true
        ),
        Err(ResonanceSearchError::Cancelled)
    ));
    for observation in [
        Observation::LayerAbsorbedPowerFraction { layer_index: 1 },
        Observation::LayerAbsorbedPowerDensity { layer_index: 1 },
    ] {
        assert!(matches!(
            stack.driven_resonance(
                1.0,
                Polarization::Te,
                Constraint::FixedAngle,
                observation,
                interval,
                controls(false),
                || false
            ),
            Err(ResonanceSearchError::Source(_))
        ));
    }
}
