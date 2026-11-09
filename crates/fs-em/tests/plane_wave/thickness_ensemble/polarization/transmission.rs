//! G1/G3 common-plane Jones transport and positive downstream heat.
use super::*;

fn debye(base: f64, delta: f64) -> ResolvedDebyeCapacitor {
    let state = source(
        &[
            (EPS_NAME, EPS, base),
            (MU_NAME, MU, 1.0),
            ("delta", EPS, delta),
            ("tau", Time::DIMS, 1.0),
        ],
        ValidityDomain::unconstrained(),
        QueryPoint::new(),
    );
    ResolvedDebyeCapacitor::try_from_material(
        state,
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

fn propagation(q: C64, distance: f64) -> C64 {
    C64::new((q.re * distance).cos(), -(q.re * distance).sin()).scale((q.im * distance).exp())
}

// Independent one-film Airy fields, with the actual unequal exit admittance.
fn transmitted_airy(
    film: C64,
    exit: C64,
    thickness: f64,
    z: f64,
    cosine: f64,
    polarization: Polarization,
) -> (C64, f64) {
    let tangent_squared = (1.0 - cosine) * (1.0 + cosine);
    let qf = (film - C64::from_re(tangent_squared)).sqrt();
    let qe = (exit - C64::from_re(tangent_squared)).sqrt();
    let (yi, yf, ye) = match polarization {
        Polarization::Te => (cosine, qf, qe),
        Polarization::Tm => (1.0 / cosine, film / qf, exit / qe),
    };
    let y = C64::from_re(yi);
    let r01 = (y - yf) / (y + yf);
    let r12 = (yf - ye) / (yf + ye);
    let p = propagation(qf, thickness);
    let t = (C64::ONE + r01) * (C64::ONE + r12) * p / (C64::ONE + r01 * r12 * p * p);
    let power = t.norm_sq() * ye.re / yi;
    let loss = power * -(2.0 * qe.im * (z - thickness)).exp_m1();
    (
        t.scale((ye.re / yi).sqrt()) * propagation(qe, z - thickness),
        loss,
    )
}

#[test]
fn common_plane_removes_false_fluctuations_from_homogeneous_exit_motion() {
    let outside = material(1.0, 1.0);
    let stack = PlanarDielectricLayers::try_new(
        &outside,
        DielectricResponse::Lossless(&outside),
        1.0,
        &[DielectricLayer {
            response: DielectricResponse::Lossless(&outside),
            thickness_m: 0.0,
        }],
    )
    .unwrap();
    let geometry = [[0.0], [0.2], [0.8]];
    let probabilities = [0.2, 0.3, 0.5];
    let rows: Vec<_> = geometry
        .iter()
        .zip(probabilities)
        .map(|(d, probability)| LayerThicknessRealization {
            thicknesses_m: d,
            probability,
        })
        .collect();
    let ensemble = stack
        .prepare_thickness_ensemble(&rows, budget(3, 1), || false)
        .unwrap();
    let out = ensemble
        .evaluate_polarized_at_plane(0.6, 1.3, 6, || false)
        .unwrap();
    assert_eq!(
        out.reflected,
        ensemble.evaluate_polarized(0.6, 6, || false).unwrap()
    );
    assert_eq!(out.transmitted.reference_plane_m(), 1.3);
    let phase = propagation(C64::from_re(0.6), 1.3);
    for t in out.transmitted.mean_transmitted_power_wave() {
        complex_near(t, phase, 8e-16);
    }
    matrix_near(
        out.transmitted.transmitted_response_moment(),
        [[C64::ONE; 2]; 2],
        1e-15,
    );
    assert!(
        out.transmitted
            .fluctuating_response_moment()
            .iter()
            .flatten()
            .all(|z| z.abs() < 1e-28)
    );
    assert_eq!(
        out.transmitted.exit_medium_absorbed_power_fractions(),
        [0.0; 2]
    );
    let naive = geometry
        .iter()
        .zip(probabilities)
        .fold(C64::ZERO, |sum, (d, w)| {
            sum + propagation(C64::from_re(0.6), d[0]).scale(w)
        });
    assert!(
        naive.norm_sq() < 0.97,
        "exit-face field averaging invents a fluctuation"
    );
    // The common plane may be much farther away than the spread in exit
    // positions. Preserve that spread before reducing the carrier phase.
    let far = ensemble
        .evaluate_polarized_at_plane(0.6, 1e16, 6, || false)
        .unwrap()
        .transmitted;
    let phase =
        fs_math::det::sin_cos_dd(fs_math::dd::Dd::from_f64(0.6) * fs_math::dd::Dd::from_f64(1e16));
    for t in far.mean_transmitted_power_wave() {
        complex_near(t, C64::new(phase.1.to_f64(), -phase.0.to_f64()), 1e-15);
    }
    assert!(
        far.fluctuating_response_moment()
            .iter()
            .flatten()
            .all(|z| z.abs() < 1e-28)
    );
}

#[test]
fn lossy_exit_and_film_match_rowwise_airy_and_keep_two_heat_regions_distinct() {
    let outside = material(1.0, 1.0);
    let film = debye(2.0, 2.0);
    let exit = debye(1.5, 1.0);
    let stack = PlanarDielectricLayers::try_new(
        &outside,
        DielectricResponse::Debye(&exit),
        1.0,
        &[DielectricLayer {
            response: DielectricResponse::Debye(&film),
            thickness_m: 0.0,
        }],
    )
    .unwrap();
    let geometry = [[0.0], [0.3], [0.9]];
    let probabilities = [0.2, 0.3, 0.5];
    let rows: Vec<_> = geometry
        .iter()
        .zip(probabilities)
        .map(|(d, probability)| LayerThicknessRealization {
            thicknesses_m: d,
            probability,
        })
        .collect();
    let ensemble = stack
        .prepare_thickness_ensemble(&rows, budget(3, 1), || false)
        .unwrap();
    let response = ensemble
        .evaluate_polarized_at_plane(0.6, 1.4, 6, || false)
        .unwrap();
    let samples: Vec<_> = geometry
        .iter()
        .map(|d| {
            [Polarization::Te, Polarization::Tm].map(|pol| {
                transmitted_airy(
                    C64::new(3.0, -1.0),
                    C64::new(2.0, -0.5),
                    d[0],
                    1.4,
                    0.6,
                    pol,
                )
            })
        })
        .collect();
    let mut mean = [C64::ZERO; 2];
    let mut moment = [[C64::ZERO; 2]; 2];
    let mut loss = [0.0; 2];
    for (row, w) in samples.iter().zip(probabilities) {
        let dyad = pure_coherency(row.map(|x| x.0));
        for i in 0..2 {
            mean[i] = mean[i] + row[i].0.scale(w);
            loss[i] += w * row[i].1;
            for j in 0..2 {
                moment[i][j] = moment[i][j] + dyad[i][j].scale(w);
            }
        }
    }
    let transmitted = response.transmitted;
    for i in 0..2 {
        complex_near(transmitted.mean_transmitted_power_wave()[i], mean[i], 3e-14);
        near(
            transmitted.exit_medium_absorbed_power_fractions()[i],
            loss[i],
            3e-14,
        );
    }
    matrix_near(transmitted.transmitted_response_moment(), moment, 3e-14);
    matrix_near(
        transmitted.coherent_response_moment(),
        pure_coherency(mean),
        3e-14,
    );
    let input = pure_coherency([C64::from_re(0.6), C64::new(0.0, 0.8)]);
    let r = response.reflected.reflect_coherency(input).unwrap();
    let t = transmitted.transmit_coherency(input).unwrap();
    assert!(r.absorbed_power > 0.1);
    assert!(t.exit_medium_absorbed_power > 0.1);
    near(
        trace(r.reflected_coherency)
            + trace(t.transmitted_coherency)
            + r.absorbed_power
            + t.exit_medium_absorbed_power,
        1.0,
        3e-14,
    );
    let j = t.transmitted_coherency;
    assert!(j[0][0].re * j[1][1].re - j[0][1].norm_sq() > 1e-7);
    for i in 0..2 {
        for k in 0..2 {
            complex_near(
                t.coherent_transmitted_coherency[i][k] + t.fluctuating_transmitted_coherency[i][k],
                j[i][k],
                3e-15,
            );
        }
    }
    let input_stokes = [
        trace(input),
        input[0][0].re - input[1][1].re,
        2.0 * input[0][1].re,
        2.0 * input[0][1].im,
    ];
    let out_stokes = transmitted.transmitted_mueller().map(|row| {
        row.into_iter()
            .zip(input_stokes)
            .map(|(a, b)| a * b)
            .sum::<f64>()
    });
    for (a, b) in out_stokes.into_iter().zip([
        trace(j),
        j[0][0].re - j[1][1].re,
        2.0 * j[0][1].re,
        2.0 * j[0][1].im,
    ]) {
        near(a, b, 2e-15);
    }
    let later = ensemble
        .evaluate_polarized_at_plane(0.6, 2.1, 6, || false)
        .unwrap()
        .transmitted;
    let q = (C64::new(2.0, -0.5) - C64::from_re(0.64)).sqrt();
    let shift = propagation(q, 0.7);
    for i in 0..2 {
        complex_near(
            later.mean_transmitted_power_wave()[i],
            transmitted.mean_transmitted_power_wave()[i] * shift,
            3e-14,
        );
        near(
            later.exit_medium_absorbed_power_fractions()[i],
            loss[i] + moment[i][i].re * -(2.0 * q.im * 0.7).exp_m1(),
            3e-14,
        );
    }
    matrix_near(
        later.transmitted_response_moment(),
        moment.map(|r| r.map(|z| z.scale(shift.norm_sq()))),
        3e-14,
    );
    let opaque = ensemble
        .evaluate_polarized_at_plane(0.6, 1e8, 6, || false)
        .unwrap();
    assert_eq!(
        opaque.transmitted.transmitted_response_moment(),
        [[C64::ZERO; 2]; 2]
    );
    for i in 0..2 {
        near(
            opaque.transmitted.exit_medium_absorbed_power_fractions()[i],
            opaque.reflected.channels()[i].mean_transmitted_power_fraction,
            2e-15,
        );
    }
}

#[test]
fn single_row_flux_normalization_and_tir_closed_channel() {
    let outside = material(4.0, 1.0);
    let glass = material(1.0, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&outside, DielectricResponse::Lossless(&glass), 1.0, &[])
            .unwrap();
    let rows = [LayerThicknessRealization {
        thicknesses_m: &[],
        probability: 1.0,
    }];
    let ensemble = stack
        .prepare_thickness_ensemble(&rows, budget(1, 0), || false)
        .unwrap();
    for cosine in [1.0, 0.95, 0.5] {
        let result = ensemble
            .evaluate_polarized_at_plane(cosine, 0.2, 2, || false)
            .unwrap();
        assert_eq!(
            result.transmitted.fluctuating_response_moment(),
            [[C64::ZERO; 2]; 2]
        );
        for (i, pol) in [Polarization::Te, Polarization::Tm].into_iter().enumerate() {
            let native = stack.response(cosine, pol).unwrap();
            near(
                result.transmitted.transmitted_response_moment()[i][i].re,
                native.transmitted_power_fraction,
                1e-15,
            );
            if cosine == 0.5 {
                assert_eq!(native.transmitted_power_fraction, 0.0);
            } else {
                let phase = native
                    .transmitted_tangential_e
                    .scale(1.0 / native.transmitted_tangential_e.abs());
                complex_near(
                    result.transmitted.mean_transmitted_power_wave()[i],
                    phase.scale(native.transmitted_power_fraction.sqrt())
                        * propagation(native.exit_normal_wavenumber_rad_m, 0.2),
                    2e-15,
                );
            }
        }
    }
}

#[test]
fn tiny_geometry_covariance_and_fixed_partition_keep_native_precision() {
    let outside = material(1.0, 1.0);
    let film = material(2.25, 1.0);
    let stack = PlanarDielectricLayers::try_new(
        &outside,
        DielectricResponse::Lossless(&outside),
        1.0,
        &[DielectricLayer {
            response: DielectricResponse::Lossless(&film),
            thickness_m: 0.3,
        }],
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
    let result = stack
        .prepare_thickness_ensemble(&rows, budget(2, 1), || false)
        .unwrap()
        .evaluate_polarized_at_plane(0.7, 1.0, 4, || false)
        .unwrap()
        .transmitted;
    let covariance = result.fluctuating_response_moment();
    for i in 0..2 {
        assert!(covariance[i][i].re > 0.0 && covariance[i][i].re < 1e-18);
    }
    assert!(covariance[0][1].abs() > 0.0);
    let same_total = [[0.2, 0.6], [0.5, 0.3]];
    let films = [DielectricLayer {
        response: DielectricResponse::Lossless(&film),
        thickness_m: 0.4,
    }; 2];
    let stack = PlanarDielectricLayers::try_new(
        &outside,
        DielectricResponse::Lossless(&outside),
        1.0,
        &films,
    )
    .unwrap();
    let rows: Vec<_> = same_total
        .iter()
        .map(|d| LayerThicknessRealization {
            thicknesses_m: d,
            probability: 0.5,
        })
        .collect();
    let result = stack
        .prepare_thickness_ensemble(&rows, budget(2, 2), || false)
        .unwrap()
        .evaluate_polarized_at_plane(0.7, 1.0, 8, || false)
        .unwrap()
        .transmitted;
    assert!(
        result
            .fluctuating_response_moment()
            .iter()
            .flatten()
            .all(|z| z.abs() < 1e-28)
    );
}

#[test]
fn common_plane_admission_budget_cancellation_and_dark_rows_are_transactional() {
    let outside = material(1.0, 1.0);
    let stack = PlanarDielectricLayers::try_new(
        &outside,
        DielectricResponse::Lossless(&outside),
        1.0,
        &[DielectricLayer {
            response: DielectricResponse::Lossless(&outside),
            thickness_m: 0.0,
        }],
    )
    .unwrap();
    let geometry = [[0.2], [1.0]];
    let rows = [
        LayerThicknessRealization {
            thicknesses_m: &geometry[0],
            probability: 1.0,
        },
        LayerThicknessRealization {
            thicknesses_m: &geometry[1],
            probability: 0.0,
        },
    ];
    let ensemble = stack
        .prepare_thickness_ensemble(&rows, budget(2, 1), || false)
        .unwrap();
    for plane in [f64::NAN, f64::INFINITY, -1.0, 0.8] {
        assert!(
            ensemble
                .evaluate_polarized_at_plane(0.7, plane, 4, || false)
                .is_err()
        );
    }
    let mut polls = 0;
    assert!(
        ensemble
            .evaluate_polarized_at_plane(0.7, 1.0, 3, || {
                polls += 1;
                false
            })
            .is_err()
    );
    assert_eq!(polls, 1);
    let expected = ensemble
        .evaluate_polarized_at_plane(0.7, 1.0, 4, || false)
        .unwrap();
    for stop in 1..=7 {
        let mut polls = 0;
        assert!(
            ensemble
                .evaluate_polarized_at_plane(0.7, 1.0, 4, || {
                    polls += 1;
                    polls == stop
                })
                .is_err()
        );
        assert_eq!(
            ensemble
                .evaluate_polarized_at_plane(0.7, 1.0, 4, || false)
                .unwrap(),
            expected
        );
    }
    let zero = expected
        .transmitted
        .transmit_coherency([[C64::ZERO; 2]; 2])
        .unwrap();
    assert_eq!(zero.transmitted_coherency, [[C64::ZERO; 2]; 2]);
    assert_eq!(zero.exit_medium_absorbed_power, 0.0);
    assert!(
        expected
            .transmitted
            .transmit_coherency([[C64::ONE, C64::from_re(2.0)], [C64::from_re(2.0), C64::ONE]])
            .is_err()
    );
}
