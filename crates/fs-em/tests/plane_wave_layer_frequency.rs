#![cfg(feature = "plane-wave")]
//! Source-bound film tangents and driven amplitude resonances.
//! Numerical high-contrast fixtures are not material qualification.
use fs_em::plane_wave::{
    DielectricResponse, PlanarDielectricInterface, Polarization,
    dispersion::IncidenceFrequencyConstraint as Constraint,
    layers::{DielectricLayer, LayeredResonanceObservation as Observation, PlanarDielectricLayers},
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    field_storage::{
        ELECTRIC_PERMITTIVITY_DIMS as EPS, ELECTRIC_PERMITTIVITY_PROPERTY as EPS_NAME,
        MAGNETIC_PERMEABILITY_DIMS as MU, MAGNETIC_PERMEABILITY_PROPERTY as MU_NAME,
        dielectric::{DebyePoleProperties, ResolvedDebyeCapacitor},
    },
    state_point::{
        MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
        ScalarPropertyRequirement, resolve_material_state_point,
    },
};
use fs_math::{c64::C64, dd::Dd, det};
use fs_phs::harmonic::resonance::{ResonanceControls, ResonanceSearchError};
use fs_qty::{Dims, Time};

#[path = "plane_wave_layer_frequency/absorption.rs"]
mod absorption;
#[path = "plane_wave_layer_frequency/antiresonance.rs"]
mod antiresonance;

const E: f64 = 1e-11;
const M: f64 = 1e-6;
const W: f64 = 1e9;

fn material(epsilon: f64, mu: f64) -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    let properties: [(&str, Dims, f64); 4] = [
        (EPS_NAME, EPS, epsilon),
        (MU_NAME, MU, mu),
        ("delta", EPS, 4.0 * E),
        ("tau", Time::DIMS, 1.0 / W),
    ];
    for (name, dims, value) in properties {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained(),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic coherent film frequency fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
                .unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic dielectric".into(),
            phase: "solid".into(),
            process: "coherent tangent fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    resolve_material_state_point(
        &card,
        &QueryPoint::new(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn debye(epsilon: f64) -> ResolvedDebyeCapacitor {
    ResolvedDebyeCapacitor::try_from_material(
        material(epsilon, M),
        EPS_NAME,
        &[DebyePoleProperties {
            permittivity_increment: "delta".into(),
            relaxation_time: "tau".into(),
        }],
        1.0,
        1.0,
    )
    .unwrap()
}
fn near(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:?} != {expected:?}, tolerance {tolerance}"
    );
}
fn controls(bandwidth: bool) -> ResonanceControls {
    ResonanceControls {
        absolute_frequency_tolerance_rad_s: 0.1,
        relative_frequency_tolerance: 0.0,
        maximum_evaluations: 160,
        measure_half_power_bandwidth: bandwidth,
    }
}

#[test]
fn homogeneous_transit_and_split_film_keep_both_incidence_constraints() {
    let ambient = material(E, M);
    let response = DielectricResponse::Lossless(&ambient);
    let layer = [DielectricLayer {
        response,
        thickness_m: 0.13,
    }];
    let split = [
        DielectricLayer {
            response,
            thickness_m: 0.04,
        },
        DielectricLayer {
            response,
            thickness_m: 0.0,
        },
        DielectricLayer {
            response,
            thickness_m: 0.09,
        },
    ];
    let stack = PlanarDielectricLayers::try_new(&ambient, response, W, &layer).unwrap();
    let split = PlanarDielectricLayers::try_new(&ambient, response, W, &split).unwrap();
    let cosine = 0.8;
    for polarization in [Polarization::Te, Polarization::Tm] {
        for constraint in [
            Constraint::FixedAngle,
            Constraint::FixedTransverseWavenumber,
        ] {
            let jet = stack
                .frequency_dispersion(cosine, polarization, constraint)
                .unwrap();
            let other = split
                .frequency_dispersion(cosine, polarization, constraint)
                .unwrap();
            let phase = jet.carrier().incident_normal_wavenumber_rad_m * 0.13;
            let slope = match constraint {
                Constraint::FixedAngle => phase,
                Constraint::FixedTransverseWavenumber => phase / (cosine * cosine),
            };
            let amplitude = jet.transmitted_tangential_e().unwrap();
            near(
                amplitude.scaled_derivative(),
                amplitude.amplitude() * C64::new(0.0, -slope),
                2e-13,
            );
            near(
                jet.reflected_tangential_e().scaled_derivative(),
                C64::ZERO,
                2e-13,
            );
            near(
                other
                    .transmitted_tangential_e()
                    .unwrap()
                    .scaled_derivative(),
                amplitude.scaled_derivative(),
                2e-13,
            );
            assert!((amplitude.group_delay_s().unwrap() - slope / W).abs() < 2e-22);
        }
    }
}

// Independent five-point derivative of native carriers. Production queries
// never perturb frequency and never use this oracle.
fn native_derivatives(
    ambient: &ResolvedMaterialStatePoint,
    exit: DielectricResponse<'_>,
    films: &[DielectricLayer<'_>],
    omega: f64,
    cosine: f64,
    polarization: Polarization,
    constraint: Constraint,
) -> [C64; 3] {
    let step = 1e-4 * omega;
    let samples: [[C64; 3]; 5] = core::array::from_fn(|index| {
        let frequency = omega + (index as f64 - 2.0) * step;
        let angle = match constraint {
            Constraint::FixedAngle => cosine,
            Constraint::FixedTransverseWavenumber => {
                let sine_squared = (1.0 - cosine * cosine) * (omega / frequency).powi(2);
                (1.0 - sine_squared).sqrt()
            }
        };
        let value = PlanarDielectricLayers::try_new(ambient, exit, frequency, films)
            .unwrap()
            .response(angle, polarization)
            .unwrap();
        [
            value.reflected_tangential_e,
            value.transmitted_tangential_e,
            value.transmitted_conjugate_h,
        ]
    });
    core::array::from_fn(|component| {
        (samples[0][component] - samples[1][component].scale(8.0)
            + samples[3][component].scale(8.0)
            - samples[4][component])
            .scale(omega / (12.0 * step))
    })
}

#[test]
fn material_dispersion_loss_and_magnetic_contrast_match_native_complete_stack() {
    let ambient = material(2.0 * E, M);
    let exit = material(1.4 * E, M);
    let dielectric = debye(3.0 * E);
    let second = material(7.0 * E, 1.3 * M);
    let films = [
        DielectricLayer {
            response: DielectricResponse::Debye(&dielectric),
            thickness_m: 0.14,
        },
        DielectricLayer {
            response: DielectricResponse::Lossless(&second),
            thickness_m: 0.09,
        },
    ];
    let exit = DielectricResponse::Lossless(&exit);
    let omega = 0.7 * W;
    let stack = PlanarDielectricLayers::try_new(&ambient, exit, omega, &films).unwrap();
    for polarization in [Polarization::Te, Polarization::Tm] {
        for constraint in [
            Constraint::FixedAngle,
            Constraint::FixedTransverseWavenumber,
        ] {
            let jet = stack
                .frequency_dispersion(0.73, polarization, constraint)
                .unwrap();
            assert!(jet.carrier().absorbed_power_fraction > 0.0);
            let expected = native_derivatives(
                &ambient,
                exit,
                &films,
                omega,
                0.73,
                polarization,
                constraint,
            );
            let actual = [
                jet.reflected_tangential_e(),
                jet.transmitted_tangential_e().unwrap(),
                jet.transmitted_conjugate_h().unwrap(),
            ];
            for i in 0..3 {
                near(
                    actual[i].scaled_derivative(),
                    expected[i],
                    3e-8 * expected[i].abs().max(1.0),
                );
            }
            assert!(jet.relative_equation_residual() <= 512.0 * f64::EPSILON);
        }
    }
}

#[test]
fn critical_finite_film_retains_nonzero_transverse_frequency_tangent() {
    let ambient = material(4.0 * E, M);
    let critical = material(3.0 * E, M);
    let exit = DielectricResponse::Lossless(&ambient);
    let films = [DielectricLayer {
        response: DielectricResponse::Lossless(&critical),
        thickness_m: 0.2,
    }];
    let stack = PlanarDielectricLayers::try_new(&ambient, exit, W, &films).unwrap();
    for polarization in [Polarization::Te, Polarization::Tm] {
        let jet = stack
            .frequency_dispersion(0.5, polarization, Constraint::FixedAngle)
            .unwrap();
        let depth = jet.carrier().incident_normal_wavenumber_rad_m / 0.5 * 0.2;
        let coupling = depth
            * if polarization == Polarization::Te {
                0.5
            } else {
                0.375
            };
        let denominator = C64::new(2.0, coupling);
        let derivative = (C64::new(0.0, -2.0 * coupling) / denominator) / denominator;
        near(
            jet.transmitted_tangential_e().unwrap().scaled_derivative(),
            derivative,
            2e-13,
        );
        let constrained = stack
            .frequency_dispersion(0.5, polarization, Constraint::FixedTransverseWavenumber)
            .unwrap();
        let expected = native_derivatives(
            &ambient,
            exit,
            &films,
            W,
            0.5,
            polarization,
            Constraint::FixedTransverseWavenumber,
        );
        near(
            constrained
                .transmitted_tangential_e()
                .unwrap()
                .scaled_derivative(),
            expected[1],
            2e-8,
        );
        assert!(
            (constrained
                .transmitted_tangential_e()
                .unwrap()
                .scaled_derivative()
                - derivative)
                .abs()
                > 0.01
        );
    }
}

#[test]
fn high_finesse_etalon_uses_original_peak_and_amplitude_half_height() {
    let ambient = material(E, M);
    let index: f64 = 1024.0;
    let film = material(E * index * index, M);
    let exit = DielectricResponse::Lossless(&ambient);
    let wave_number = W * det::sqrt(E) * det::sqrt(M);
    let thickness = core::f64::consts::PI / (index * wave_number);
    let films = [DielectricLayer {
        response: DielectricResponse::Lossless(&film),
        thickness_m: thickness,
    }];
    let stack = PlanarDielectricLayers::try_new(&ambient, exit, W, &films).unwrap();
    let result = stack
        .driven_resonance(
            1.0,
            Polarization::Te,
            Constraint::FixedAngle,
            Observation::TransmittedElectricAmplitude,
            (0.99 * W, 1.01 * W),
            controls(true),
            || false,
        )
        .unwrap();
    assert!((result.peak.angular_frequency_rad_s - W).abs() < 0.2);
    assert!((result.observation.abs() - 1.0).abs() < 1e-12);
    let airy = 0.5 * (index - 1.0 / index);
    let expected_bandwidth = 2.0 * W / core::f64::consts::PI * (1.0 / airy).asin();
    let measured = result.half_power.unwrap();
    assert!((measured.bandwidth_rad_s - expected_bandwidth).abs() < 0.3);
    assert!((measured.response_quality_factor - W / expected_bandwidth).abs() < 3e-4);
    assert!(result.evaluations <= 160);
    assert!(matches!(
        stack.driven_resonance(
            1.0,
            Polarization::Te,
            Constraint::FixedAngle,
            Observation::TransmittedElectricAmplitude,
            (0.99 * W, 1.01 * W),
            controls(false),
            || true,
        ),
        Err(ResonanceSearchError::Cancelled)
    ));
}

#[test]
fn extreme_finesse_tangent_preserves_native_double_double_phase() {
    let ambient = material(E, M);
    let index: f64 = 1099511627776.0; // 2^40: exact index contrast, no fit.
    let film = material(E * index * index, M);
    let wave_number = W * det::sqrt(E) * det::sqrt(M);
    let thickness = core::f64::consts::PI / (index * wave_number);
    let films = [DielectricLayer {
        response: DielectricResponse::Lossless(&film),
        thickness_m: thickness,
    }];
    let stack = PlanarDielectricLayers::try_new(
        &ambient,
        DielectricResponse::Lossless(&ambient),
        W,
        &films,
    )
    .unwrap();
    let result = stack
        .frequency_dispersion(1.0, Polarization::Te, Constraint::FixedAngle)
        .unwrap();
    let phase = Dd::from_f64(index) * (Dd::from_f64(wave_number) * Dd::from_f64(thickness));
    let (sine, cosine) = det::sin_cos_dd(phase);
    let (sine, cosine) = (sine.to_f64(), cosine.to_f64());
    let airy = 0.5 * (index + 1.0 / index);
    let denominator = C64::new(cosine, airy * sine);
    let numerator_slope = C64::new(-sine, airy * cosine).scale(phase.to_f64());
    let expected = -(numerator_slope / denominator) / denominator;
    let jet = result.transmitted_tangential_e().unwrap();
    near(jet.amplitude(), C64::ONE / denominator, 3e-12);
    near(jet.scaled_derivative(), expected, 3e-12 * expected.abs());
}

#[test]
fn opaque_native_reflection_keeps_source_slope_without_invented_transmission_delay() {
    let ambient = material(E, M);
    let dielectric = debye(2.0 * E);
    let response = DielectricResponse::Debye(&dielectric);
    let films = [DielectricLayer {
        response,
        thickness_m: 1e6,
    }];
    let stack = PlanarDielectricLayers::try_new(
        &ambient,
        DielectricResponse::Lossless(&ambient),
        W,
        &films,
    )
    .unwrap();
    let jet = stack
        .frequency_dispersion(0.8, Polarization::Tm, Constraint::FixedAngle)
        .unwrap();
    let front = PlanarDielectricInterface::try_new(&ambient, response)
        .unwrap()
        .frequency_dispersion(W, 0.8, Polarization::Tm, Constraint::FixedAngle)
        .unwrap();
    near(
        jet.reflected_tangential_e().scaled_derivative(),
        front.reflected_tangential_e().scaled_derivatives()[0],
        1e-15,
    );
    assert_eq!(jet.carrier().transmitted_tangential_e, C64::ZERO);
    assert!(jet.transmitted_tangential_e().is_none());
    assert!(jet.transmitted_conjugate_h().is_none());
    assert!(matches!(
        stack.driven_resonance(
            0.8,
            Polarization::Tm,
            Constraint::FixedAngle,
            Observation::TransmittedElectricAmplitude,
            (0.8 * W, 1.2 * W),
            controls(false),
            || false,
        ),
        Err(ResonanceSearchError::Source(_))
    ));
}

fn power_slope(
    ambient: &ResolvedMaterialStatePoint,
    exit: DielectricResponse<'_>,
    films: &[DielectricLayer<'_>],
    omega: f64,
    cosine: f64,
    polarization: Polarization,
    constraint: Constraint,
) -> [f64; 2] {
    let step = 1e-4 * omega;
    let samples: [[f64; 2]; 4] = core::array::from_fn(|i| {
        let frequency = omega + [-2.0, -1.0, 1.0, 2.0][i] * step;
        let angle = match constraint {
            Constraint::FixedAngle => cosine,
            Constraint::FixedTransverseWavenumber => {
                (1.0 - (1.0 - cosine * cosine) * (omega / frequency).powi(2)).sqrt()
            }
        };
        let sample = PlanarDielectricLayers::try_new(ambient, exit, frequency, films)
            .unwrap()
            .response(angle, polarization)
            .unwrap();
        [
            sample.reflected_power_fraction,
            sample.transmitted_power_fraction,
        ]
    });
    core::array::from_fn(|i| {
        (samples[0][i] - 8.0 * samples[1][i] + 8.0 * samples[2][i] - samples[3][i]) * omega
            / (12.0 * step)
    })
}

#[test]
fn original_interface_and_lossy_stack_flux_slopes_include_both_port_admittances() {
    let ambient = material(2.0 * E, M);
    let exit_material = debye(1.4 * E);
    let exit = DielectricResponse::Debye(&exit_material);
    let film = debye(3.0 * E);
    let second = material(7.0 * E, 1.3 * M);
    let films = [
        DielectricLayer {
            response: DielectricResponse::Debye(&film),
            thickness_m: 0.14,
        },
        DielectricLayer {
            response: DielectricResponse::Lossless(&second),
            thickness_m: 0.09,
        },
    ];
    let omega = 0.7 * W;
    let stack = PlanarDielectricLayers::try_new(&ambient, exit, omega, &films).unwrap();
    let interface = PlanarDielectricInterface::try_new(&ambient, exit).unwrap();
    let mut amplitude_and_flux_differ = false;
    for polarization in [Polarization::Te, Polarization::Tm] {
        for constraint in [
            Constraint::FixedAngle,
            Constraint::FixedTransverseWavenumber,
        ] {
            let boundary = interface
                .frequency_dispersion(omega, 0.73, polarization, constraint)
                .unwrap();
            let r = boundary.reflected_power_fraction().unwrap();
            let t = boundary.transmitted_power_fraction().unwrap();
            assert!((r.fraction() + t.fraction() - 1.0).abs() < 2e-13);
            assert!((r.scaled_derivative() + t.scaled_derivative()).abs() < 2e-13);
            let expected = power_slope(&ambient, exit, &[], omega, 0.73, polarization, constraint);
            assert!((t.scaled_derivative() - expected[1]).abs() < 3e-8);
            let field = stack
                .frequency_dispersion(0.73, polarization, constraint)
                .unwrap();
            let r = field.reflected_power_fraction().unwrap();
            let t = field.transmitted_power_fraction().unwrap().unwrap();
            let second_order = stack
                .second_order_frequency_dispersion(0.73, polarization, constraint)
                .unwrap();
            let second_r = second_order.reflected_power_fraction().unwrap();
            let second_t = second_order.transmitted_power_fraction().unwrap().unwrap();
            assert!((second_r.scaled_derivative() - r.scaled_derivative()).abs() < 2e-12);
            assert!((second_t.scaled_derivative() - t.scaled_derivative()).abs() < 2e-12);
            assert_eq!(r.fraction(), field.carrier().reflected_power_fraction);
            assert_eq!(t.fraction(), field.carrier().transmitted_power_fraction);
            let expected = power_slope(
                &ambient,
                exit,
                &films,
                omega,
                0.73,
                polarization,
                constraint,
            );
            assert!((r.scaled_derivative() - expected[0]).abs() < 3e-8);
            assert!((t.scaled_derivative() - expected[1]).abs() < 3e-8);
            let electric = field.transmitted_tangential_e().unwrap();
            let amplitude_slope =
                2.0 * (electric.amplitude().conj() * electric.scaled_derivative()).re;
            amplitude_and_flux_differ |= (amplitude_slope - t.scaled_derivative()).abs() > 1e-4;
        }
    }
    assert!(amplitude_and_flux_differ);
}

// Independent lossless single-film TE Airy law at fixed transverse wave number.
// The incident index is1, film8, exit2, carrier incidence cosine0.6, phase pi.
fn asymmetric_airy(ratio: f64) -> (f64, f64) {
    let u2 = 0.64 / (ratio * ratio);
    let cosine = (1.0 - u2).sqrt();
    let q = (64.0 - u2).sqrt();
    let qo = (4.0 - u2).sqrt();
    let g = qo / cosine;
    let dg = g * (u2 / (qo * qo) - u2 / (cosine * cosine));
    let a = q / cosine;
    let da = a * (u2 / (q * q) - u2 / (cosine * cosine));
    let b = qo / q;
    let db = b * (u2 / (qo * qo) - u2 / (q * q));
    let phase = core::f64::consts::PI * ratio * q / 63.36_f64.sqrt();
    let dphase = phase * (1.0 + u2 / (q * q));
    let (sn, cs) = phase.sin_cos();
    let left = 1.0 + g;
    let right = a + b;
    let den = cs * cs * left * left + sn * sn * right * right;
    let dden = 2.0 * cs * cs * left * dg
        + 2.0 * sn * sn * right * (da + db)
        + 2.0 * sn * cs * dphase * (right * right - left * left);
    let power = 4.0 * g / den;
    (power, power * (dg / g - dden / den))
}

#[test]
fn transmitted_power_resonance_and_half_height_follow_actual_asymmetric_flux() {
    let ambient = material(E, M);
    let film = material(64.0 * E, M);
    let exit = material(4.0 * E, M);
    let wave_number = W * det::sqrt(E) * det::sqrt(M);
    let thickness = core::f64::consts::PI / (63.36_f64.sqrt() * wave_number);
    let films = [DielectricLayer {
        response: DielectricResponse::Lossless(&film),
        thickness_m: thickness,
    }];
    let stack =
        PlanarDielectricLayers::try_new(&ambient, DielectricResponse::Lossless(&exit), W, &films)
            .unwrap();
    let interval = (0.84 * W, 1.2 * W);
    let actual = stack
        .driven_resonance(
            0.6,
            Polarization::Te,
            Constraint::FixedTransverseWavenumber,
            Observation::TransmittedPowerFraction,
            interval,
            controls(true),
            || false,
        )
        .unwrap();
    let electric = stack
        .driven_resonance(
            0.6,
            Polarization::Te,
            Constraint::FixedTransverseWavenumber,
            Observation::TransmittedElectricAmplitude,
            interval,
            controls(false),
            || false,
        )
        .unwrap();
    let (mut lower, mut upper) = (0.84, 1.2);
    for _ in 0..64 {
        let middle = (lower + upper) * 0.5;
        if asymmetric_airy(middle).1 > 0.0 {
            lower = middle;
        } else {
            upper = middle;
        }
    }
    let peak = (lower + upper) * 0.5;
    let target = asymmetric_airy(peak).0 * 0.5;
    let crossing = |mut a: f64, mut b: f64, rising: bool| {
        for _ in 0..64 {
            let middle = (a + b) * 0.5;
            if (asymmetric_airy(middle).0 < target) == rising {
                a = middle;
            } else {
                b = middle;
            }
        }
        (a + b) * 0.5
    };
    let width = (crossing(peak, 1.2, false) - crossing(0.84, peak, true)) * W;
    assert!((actual.peak.angular_frequency_rad_s - peak * W).abs() < 0.3);
    assert!((actual.observation.norm_sq() - 2.0 * target).abs() < 2e-12);
    assert_eq!(actual.observation.im, 0.0);
    assert!((actual.half_power.unwrap().bandwidth_rad_s - width).abs() < 1.0);
    assert!(
        (actual.peak.angular_frequency_rad_s - electric.peak.angular_frequency_rad_s).abs() > 1e6
    );
    assert!(actual.evaluations <= 160);
    assert!(matches!(
        stack.driven_resonance(
            0.6,
            Polarization::Te,
            Constraint::FixedTransverseWavenumber,
            Observation::TransmittedPowerFraction,
            interval,
            controls(false),
            || true
        ),
        Err(ResonanceSearchError::Cancelled)
    ));
}
