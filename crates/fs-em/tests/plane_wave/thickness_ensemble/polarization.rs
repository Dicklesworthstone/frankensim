//! G1/G3 joint-polarization ensemble checks against rowwise Jones/Airy optics.
use super::*;

fn pure_coherency(a: [C64; 2]) -> [[C64; 2]; 2] {
    [
        [C64::from_re(a[0].norm_sq()), a[0] * a[1].conj()],
        [a[1] * a[0].conj(), C64::from_re(a[1].norm_sq())],
    ]
}
fn matrix_near(a: [[C64; 2]; 2], b: [[C64; 2]; 2], tolerance: f64) {
    for i in 0..2 {
        for j in 0..2 {
            complex_near(a[i][j], b[i][j], tolerance);
        }
    }
}
fn trace(a: [[C64; 2]; 2]) -> f64 {
    a[0][0].re + a[1][1].re
}

#[test]
fn paired_channels_retain_native_response_and_single_row_has_no_covariance() {
    let outside = material(1.0, 1.0);
    let film = material(9.0, 1.0);
    let stack = PlanarDielectricLayers::try_new(
        &outside,
        DielectricResponse::Lossless(&outside),
        1.0,
        &[DielectricLayer {
            response: DielectricResponse::Lossless(&film),
            thickness_m: 0.7,
        }],
    )
    .unwrap();
    let geometry = [0.7];
    let rows = [LayerThicknessRealization {
        thicknesses_m: &geometry,
        probability: 1.0,
    }];
    let prepared = stack
        .prepare_thickness_ensemble(&rows, budget(1, 1), || false)
        .unwrap();
    let paired = prepared.evaluate_polarized(0.6, 2, || false).unwrap();
    for (index, polarization) in [Polarization::Te, Polarization::Tm].into_iter().enumerate() {
        assert_eq!(
            paired.channels()[index],
            prepared.evaluate(0.6, polarization, || false).unwrap()
        );
    }
    assert_eq!(paired.fluctuating_response_moment(), [[C64::ZERO; 2]; 2]);
    matrix_near(
        paired.reflected_response_moment(),
        paired.coherent_response_moment(),
        2e-16,
    );
    let r = paired.channels().map(|x| x.mean_reflected_tangential_e);
    for input in [
        [C64::ONE, C64::ZERO],
        [C64::ZERO, C64::ONE],
        [C64::from_re(0.6), C64::new(0.0, 0.8)],
    ] {
        let incident = pure_coherency(input);
        let result = paired.reflect_coherency(incident).unwrap();
        matrix_near(
            result.reflected_coherency,
            pure_coherency([r[0] * input[0], r[1] * input[1]]),
            3e-16,
        );
        near(
            trace(result.reflected_coherency) + result.transmitted_power + result.absorbed_power,
            trace(incident),
            2e-14,
        );
    }
}

#[test]
fn same_geometry_preserves_complex_cross_correlation_and_real_depolarization() {
    let outside = material(1.0, 1.0);
    let state = source(
        &[
            (EPS_NAME, EPS, 2.0),
            (MU_NAME, MU, 1.0),
            ("delta", EPS, 2.0),
            ("tau", Time::DIMS, 1.0),
        ],
        ValidityDomain::unconstrained(),
        QueryPoint::new(),
    );
    let dielectric = ResolvedDebyeCapacitor::try_from_material(
        state,
        EPS_NAME,
        &[DebyePoleProperties {
            permittivity_increment: "delta".into(),
            relaxation_time: "tau".into(),
        }],
        1.0,
        1.0,
    )
    .unwrap();
    let stack = PlanarDielectricLayers::try_new(
        &outside,
        DielectricResponse::Lossless(&outside),
        1.0,
        &[DielectricLayer {
            response: DielectricResponse::Debye(&dielectric),
            thickness_m: 0.7,
        }],
    )
    .unwrap();
    let geometry = [[0.0], [0.2], [0.7], [2.0]];
    let probabilities = [0.1, 0.2, 0.3, 0.4];
    let rows: Vec<_> = geometry
        .iter()
        .zip(probabilities)
        .map(|(d, probability)| LayerThicknessRealization {
            thicknesses_m: d,
            probability,
        })
        .collect();
    let prepared = stack
        .prepare_thickness_ensemble(&rows, budget(4, 1), || false)
        .unwrap();
    let result = prepared.evaluate_polarized(0.6, 8, || false).unwrap();
    let epsilon = C64::new(3.0, -1.0);
    let amplitudes: Vec<_> = geometry
        .iter()
        .map(|d| {
            [
                airy(epsilon, d[0], 0.6, Polarization::Te).0,
                airy(epsilon, d[0], 0.6, Polarization::Tm).0,
            ]
        })
        .collect();
    let mut moment = [[C64::ZERO; 2]; 2];
    for (r, probability) in amplitudes.iter().zip(probabilities) {
        let dyad = pure_coherency(*r);
        for i in 0..2 {
            for j in 0..2 {
                moment[i][j] = moment[i][j] + dyad[i][j].scale(probability);
            }
        }
    }
    matrix_near(result.reflected_response_moment(), moment, 3e-14);
    assert!(
        moment[0][1].im.abs() > 1e-4,
        "relative TE/TM phase must be retained"
    );
    assert!((moment[0][1] - result.coherent_response_moment()[0][1]).abs() > 1e-4);
    let input = pure_coherency([C64::from_re(0.6), C64::new(0.0, 0.8)]);
    let out = result.reflect_coherency(input).unwrap();
    let mut rowwise = [[C64::ZERO; 2]; 2];
    for (r, probability) in amplitudes.iter().zip(probabilities) {
        let dyad = pure_coherency([r[0].scale(0.6), r[1] * C64::new(0.0, 0.8)]);
        for i in 0..2 {
            for j in 0..2 {
                rowwise[i][j] = rowwise[i][j] + dyad[i][j].scale(probability);
            }
        }
    }
    matrix_near(out.reflected_coherency, rowwise, 3e-14);
    let j = out.reflected_coherency;
    assert!(
        j[0][0].re * j[1][1].re - j[0][1].norm_sq() > 1e-7,
        "the ensemble depolarizes even though each row is deterministic Jones optics"
    );
    for i in 0..2 {
        for j in 0..2 {
            complex_near(
                out.coherent_reflected_coherency[i][j] + out.fluctuating_reflected_coherency[i][j],
                out.reflected_coherency[i][j],
                3e-16,
            );
        }
    }
    assert!(out.absorbed_power > 0.2);
    near(
        trace(out.reflected_coherency) + out.transmitted_power + out.absorbed_power,
        1.0,
        3e-14,
    );
    let stokes = [
        trace(input),
        input[0][0].re - input[1][1].re,
        2.0 * input[0][1].re,
        2.0 * input[0][1].im,
    ];
    let y = result
        .reflected_mueller()
        .map(|row| row.iter().zip(stokes).map(|(a, b)| a * b).sum::<f64>());
    near(y[0], trace(j), 3e-16);
    near(y[1], j[0][0].re - j[1][1].re, 3e-16);
    near(y[2], 2.0 * j[0][1].re, 3e-16);
    near(y[3], 2.0 * j[0][1].im, 3e-16);
}

#[test]
fn normal_incidence_fluctuation_does_not_invent_depolarization_or_heat() {
    let outside = material(1.0, 1.0);
    let film = material(9.0, 1.0);
    let stack = PlanarDielectricLayers::try_new(
        &outside,
        DielectricResponse::Lossless(&outside),
        1.0,
        &[DielectricLayer {
            response: DielectricResponse::Lossless(&film),
            thickness_m: 0.0,
        }],
    )
    .unwrap();
    let geometry: Vec<_> = (0..64)
        .map(|i| [std::f64::consts::PI * (i as f64 + 0.5) / 192.0])
        .collect();
    let rows: Vec<_> = geometry
        .iter()
        .map(|d| LayerThicknessRealization {
            thicknesses_m: d,
            probability: 1.0 / 64.0,
        })
        .collect();
    let prepared = stack
        .prepare_thickness_ensemble(&rows, budget(64, 1), || false)
        .unwrap();
    let response = prepared.evaluate_polarized(1.0, 128, || false).unwrap();
    let input = pure_coherency([C64::from_re(0.6), C64::new(0.0, 0.8)]);
    let output = response.reflect_coherency(input).unwrap();
    matrix_near(
        output.reflected_coherency,
        input.map(|row| row.map(|z| z.scale(0.4))),
        2e-14,
    );
    matrix_near(
        output.coherent_reflected_coherency,
        input.map(|row| row.map(|z| z.scale(0.25))),
        2e-14,
    );
    matrix_near(
        output.fluctuating_reflected_coherency,
        input.map(|row| row.map(|z| z.scale(0.15))),
        2e-14,
    );
    assert_eq!(output.absorbed_power, 0.0);
    near(output.transmitted_power, 0.6, 2e-14);
}

#[test]
fn weak_channels_and_tiny_geometry_covariance_remain_representable() {
    let outside = material(1.0, 1.0);
    let glass = material(2.25, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&outside, DielectricResponse::Lossless(&glass), 1.0, &[])
            .unwrap();
    let rows = [LayerThicknessRealization {
        thicknesses_m: &[],
        probability: 1.0,
    }];
    let response = stack
        .prepare_thickness_ensemble(&rows, budget(1, 0), || false)
        .unwrap()
        .evaluate_polarized(1.0 / 3.25_f64.sqrt() + 1e-12, 2, || false)
        .unwrap();
    let p = response.channels()[1].mean_reflected_power_fraction;
    assert!(p > 0.0 && p < 1e-20);
    let pure_p = [[C64::ZERO, C64::ZERO], [C64::ZERO, C64::ONE]];
    assert_eq!(
        response
            .reflect_coherency(pure_p)
            .unwrap()
            .reflected_coherency[1][1]
            .re,
        p
    );
    let layer = [DielectricLayer {
        response: DielectricResponse::Lossless(&glass),
        thickness_m: 0.3,
    }];
    let stack = PlanarDielectricLayers::try_new(
        &outside,
        DielectricResponse::Lossless(&outside),
        1.0,
        &layer,
    )
    .unwrap();
    let geometry = [[0.3], [0.3 + 1e-10]];
    let rows: Vec<_> = geometry
        .iter()
        .map(|d| LayerThicknessRealization {
            thicknesses_m: d,
            probability: 0.5,
        })
        .collect();
    let response = stack
        .prepare_thickness_ensemble(&rows, budget(2, 1), || false)
        .unwrap()
        .evaluate_polarized(0.7, 4, || false)
        .unwrap();
    let covariance = response.fluctuating_response_moment();
    for index in 0..2 {
        assert!(covariance[index][index].re > 0.0 && covariance[index][index].re < 1e-18);
    }
    assert!(covariance[0][1].abs() > 0.0);
}

#[test]
fn joint_budget_cancellation_replay_and_polarization_admission() {
    let outside = material(1.0, 1.0);
    let glass = material(2.25, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&outside, DielectricResponse::Lossless(&glass), 1.0, &[])
            .unwrap();
    let rows = [LayerThicknessRealization {
        thicknesses_m: &[],
        probability: 1.0,
    }];
    let prepared = stack
        .prepare_thickness_ensemble(&rows, budget(1, 0), || false)
        .unwrap();
    let mut polls = 0;
    assert!(
        prepared
            .evaluate_polarized(0.7, 1, || {
                polls += 1;
                false
            })
            .is_err()
    );
    assert_eq!(polls, 1, "budget refuses before solving either channel");
    let expected = prepared.evaluate_polarized(0.7, 2, || false).unwrap();
    for stop in 1..=4 {
        let mut polls = 0;
        assert!(
            prepared
                .evaluate_polarized(0.7, 2, || {
                    polls += 1;
                    polls == stop
                })
                .is_err()
        );
        assert_eq!(
            prepared.evaluate_polarized(0.7, 2, || false).unwrap(),
            expected
        );
    }
    let zero = expected.reflect_coherency([[C64::ZERO; 2]; 2]).unwrap();
    assert_eq!(zero.reflected_coherency, [[C64::ZERO; 2]; 2]);
    assert_eq!(zero.transmitted_power, 0.0);
    assert_eq!(zero.absorbed_power, 0.0);
    for invalid in [
        [[C64::from_re(-1.0), C64::ZERO], [C64::ZERO, C64::ONE]],
        [[C64::ONE, C64::from_re(2.0)], [C64::from_re(2.0), C64::ONE]],
        [
            [C64::ONE, C64::new(0.0, 1.0)],
            [C64::new(0.0, 1.0), C64::ONE],
        ],
        [[C64::new(1.0, 0.1), C64::ZERO], [C64::ZERO, C64::ONE]],
        [[C64::from_re(f64::NAN), C64::ZERO], [C64::ZERO, C64::ONE]],
    ] {
        assert!(expected.reflect_coherency(invalid).is_err());
    }
    // The positive covariance bound must not overflow/underflow.
    let extreme = [
        [C64::from_re(1e300), C64::ONE],
        [C64::ONE, C64::from_re(1e-300)],
    ];
    let out = expected.reflect_coherency(extreme).unwrap();
    assert!(out.reflected_coherency[0][0].re.is_finite());
    assert!(out.reflected_coherency[1][1].re > 0.0);
}
