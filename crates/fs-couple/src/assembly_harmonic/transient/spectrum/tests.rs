use super::super::joint::PreparedJointFlow;
use super::super::tests::{budget, maxwell, source, system};
use super::*;
use fs_vfit::{mimo::MimoModel, model::PoleTerm};

fn controls() -> AssemblyResonanceControls {
    AssemblyResonanceControls {
        numerical: DenseSpectrumOptions {
            maximum_dimension: 16,
            maximum_dense_bytes: 4 << 20,
            maximum_work: 10_000_000_000,
            qr_max_sweeps_per_window: 160,
            svd_max_sweeps: 100,
            svd_tolerance: 2e-14,
            maximum_relative_residual: 1e-11,
            minimum_eigenbasis_reciprocal_condition: 1e-9,
        },
        reference_angular_frequency_rad_s: 3.,
        maximum_payload_bytes: 1 << 20,
        maximum_back_transform_terms: 1_000_000,
        maximum_original_relative_residual: 1e-10,
    }
}
fn near(a: C64, b: C64) {
    assert!((a - b).abs() < 3e-11 * (1. + b.abs()), "{a:?} != {b:?}");
}
fn expected(values: &[AssemblyResonancePole], want: &[C64]) {
    assert_eq!(values.len(), want.len());
    for value in want {
        assert!(
            values
                .iter()
                .any(|p| (p.laplace_pole_per_s - *value).abs() < 3e-11 * (1. + value.abs())),
            "missing {value:?}"
        );
    }
}

#[test]
fn g1_loaded_maxwell_cubic_frequency_decay_and_fourier_sign_survive_units_and_clock() {
    let sys = system(&[6.], &[1.], &[0.], 1);
    let b = 35_f64.sqrt() / 2.;
    for dt in [0.02, 0.8] {
        let load = maxwell(dt);
        let joint = PreparedJointFlow::new(
            &sys,
            &[6.],
            &[1.],
            &[0.],
            &[1.],
            &load,
            budget(),
            &mut || false,
        )
        .unwrap();
        let scales = [6_f64.sqrt(), 1., joint.extra_scales[0]];
        let (poles, _) = analyze_generator(
            &joint.generator_per_s,
            1,
            &scales,
            [0.01, 100_000.],
            controls(),
            &mut || false,
        )
        .unwrap();
        expected(
            &poles,
            &[C64::from_re(-2.), C64::new(-0.5, -b), C64::new(-0.5, b)],
        );
        for pole in &poles {
            let s = pole.laplace_pole_per_s;
            near(pole.angular_frequency_rad_s, C64::new(-s.im, s.re));
            assert!(pole.angular_frequency_rad_s.im < 0.);
            assert!(
                pole.kinematic_relative_residual < 1e-11
                    && pole.momentum_relative_residual < 1e-11
                    && pole.load_relative_residual < 1e-11
            );
            if s.im.abs() > 1. {
                assert!((pole.frequency_hz() - b / core::f64::consts::TAU).abs() < 1e-11);
                assert!((pole.decay_rate_per_s() - 0.5).abs() < 1e-11);
                assert!((pole.quality_factor().unwrap() - b).abs() < 1e-10);
            }
        }
        // Canonical unit change q_new=1e-3*q, p_new=1e3*p plus a load-state
        // chart rescaling. This is a similarity, not altered mass or damping.
        let units = [1e-3, 1e3, 1e-2];
        let mut changed = joint.generator_per_s.clone();
        for i in 0..3 {
            for j in 0..3 {
                changed[i * 3 + j] *= units[i] / units[j];
            }
        }
        let rescaled = [
            scales[0] / units[0],
            scales[1] / units[1],
            scales[2] / units[2],
        ];
        let (changed_poles, _) = analyze_generator(
            &changed,
            1,
            &rescaled,
            [0.01, 100_000.],
            controls(),
            &mut || false,
        )
        .unwrap();
        expected(
            &changed_poles,
            &poles
                .iter()
                .map(|p| p.laplace_pole_per_s)
                .collect::<Vec<_>>(),
        );
    }
}

#[test]
fn g3_high_q_memory_loading_resolves_the_physical_decay_without_a_loss_floor() {
    // The complete characteristic polynomial is exactly
    // (s+a)*(s^2+2*sigma*s+omega_n^2); the small damping comes from memory,
    // with no structural viscous term or direct feedthrough resistance.
    let (a, sigma, omega_n_squared) = (2., 1e-6, 9.);
    let alpha = a + 2. * sigma;
    let stiffness = a * omega_n_squared / alpha;
    let residue = omega_n_squared + 2. * a * sigma - stiffness;
    let model = MimoModel {
        ports: 1,
        poles: vec![PoleTerm::Real {
            pole: -alpha,
            residue: 0.,
        }],
        residues: vec![vec![C64::from_re(residue)]],
        d: vec![0.],
        e: vec![0.],
        reciprocal: true,
    };
    let load = source(&model, 0.1);
    let sys = system(&[stiffness], &[1.], &[0.], 1);
    let joint = PreparedJointFlow::new(
        &sys,
        &[stiffness],
        &[1.],
        &[0.],
        &[1.],
        &load,
        budget(),
        &mut || false,
    )
    .unwrap();
    let scales = [stiffness.sqrt(), 1., joint.extra_scales[0]];
    let (poles, _) = analyze_generator(
        &joint.generator_per_s,
        1,
        &scales,
        [0.01, 100_000.],
        controls(),
        &mut || false,
    )
    .unwrap();
    let oscillation = (omega_n_squared - sigma * sigma).sqrt();
    expected(
        &poles,
        &[
            C64::from_re(-a),
            C64::new(-sigma, -oscillation),
            C64::new(-sigma, oscillation),
        ],
    );
    for pole in poles
        .iter()
        .filter(|pole| pole.laplace_pole_per_s.im.abs() > 1.)
    {
        assert!((pole.decay_rate_per_s() / sigma - 1.).abs() < 2e-8);
        assert!((pole.quality_factor().unwrap() / (oscillation / (2. * sigma)) - 1.).abs() < 2e-8);
        assert!(pole.momentum_relative_residual < 1e-11 && pole.load_relative_residual < 1e-11);
    }
}

#[test]
fn g3_full_mass_glue_and_nonreciprocal_mimo_poles_reproduce_the_exact_transient() {
    let m = [2., 0.2, 0.2, 1.];
    let k = [8., -2., -2., 5.];
    let c = [0.26, -0.06, -0.06, 0.2];
    let determinant = m[0] * m[3] - m[1] * m[2];
    let inverse = [
        m[3] / determinant,
        -m[1] / determinant,
        -m[2] / determinant,
        m[0] / determinant,
    ];
    let b = [1., 0.3, -0.2, 0.7];
    let model = MimoModel {
        ports: 2,
        poles: vec![PoleTerm::Real {
            pole: -3.,
            residue: 0.,
        }],
        residues: vec![[5., 1., 1., 2.].map(C64::from_re).to_vec()],
        d: vec![0.3, 0.08, -0.02, 0.2],
        e: vec![0.; 4],
        reciprocal: false,
    };
    let dt = 0.35;
    let load = source(&model, dt);
    let history = load.pole_states().to_vec();
    let sys = system(&k, &inverse, &c, 2);
    let joint =
        PreparedJointFlow::new(&sys, &k, &inverse, &c, &b, &load, budget(), &mut || false).unwrap();
    let mut scales = vec![
        k[0].sqrt(),
        k[3].sqrt(),
        inverse[0].sqrt(),
        inverse[3].sqrt(),
    ];
    scales.extend_from_slice(&joint.extra_scales);
    let (poles, _) = analyze_generator(
        &joint.generator_per_s,
        2,
        &scales,
        [0.01, 100_000.],
        controls(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(poles.len(), 6);
    for pole in poles {
        let real: Vec<_> = pole.state.iter().map(|v| v.re).collect();
        let imag: Vec<_> = pole.state.iter().map(|v| v.im).collect();
        let re = joint.flow.evaluate(&real, &[0.; 2]).unwrap();
        let im = joint.flow.evaluate(&imag, &[0.; 2]).unwrap();
        let (sin, cos) = (pole.laplace_pole_per_s.im * dt).sin_cos();
        let factor = C64::new(cos, sin).scale((pole.laplace_pole_per_s.re * dt).exp());
        for i in 0..6 {
            near(C64::new(re.state[i], im.state[i]), factor * pole.state[i]);
        }
    }
    assert_eq!(load.pole_states(), history);
    assert_eq!(load.source_model(), &model);
}

#[test]
fn g0_dark_memory_and_out_of_band_poles_remain_complete_and_refusals_are_atomic() {
    let load = source(
        &MimoModel {
            ports: 1,
            poles: vec![
                PoleTerm::Real {
                    pole: -3.,
                    residue: 0.,
                },
                PoleTerm::Real {
                    pole: -7.,
                    residue: 0.,
                },
            ],
            residues: vec![vec![C64::from_re(5.)], vec![C64::ZERO]],
            d: vec![0.],
            e: vec![0.],
            reciprocal: true,
        },
        0.1,
    );
    let sys = system(&[6.], &[1.], &[0.], 1);
    let joint = PreparedJointFlow::new(
        &sys,
        &[6.],
        &[1.],
        &[0.],
        &[1.],
        &load,
        budget(),
        &mut || false,
    )
    .unwrap();
    let scales = [
        6_f64.sqrt(),
        1.,
        joint.extra_scales[0],
        joint.extra_scales[1],
    ];
    let (poles, _) = analyze_generator(
        &joint.generator_per_s,
        1,
        &scales,
        [10., 100.],
        controls(),
        &mut || false,
    )
    .unwrap();
    let b = 35_f64.sqrt() / 2.;
    expected(
        &poles,
        &[
            C64::from_re(-7.),
            C64::from_re(-2.),
            C64::new(-0.5, -b),
            C64::new(-0.5, b),
        ],
    );
    assert!(
        poles
            .iter()
            .all(|p| !p.oscillation_frequency_in_source_band)
    );
    let dark = poles
        .iter()
        .find(|p| (p.laplace_pole_per_s - C64::from_re(-7.)).abs() < 1e-10)
        .unwrap();
    assert!(dark.state[..3].iter().all(|v| v.abs() < 1e-10));
    assert!(dark.state[3].abs() > 0.99);
    for cap in [
        AssemblyResonanceControls {
            maximum_payload_bytes: 1,
            ..controls()
        },
        AssemblyResonanceControls {
            maximum_back_transform_terms: 0,
            ..controls()
        },
        AssemblyResonanceControls {
            reference_angular_frequency_rad_s: 0.,
            ..controls()
        },
    ] {
        assert!(
            analyze_generator(
                &joint.generator_per_s,
                1,
                &scales,
                [1., 100.],
                cap,
                &mut || false
            )
            .is_err()
        );
    }
    let mut calls = 0;
    assert!(matches!(
        analyze_generator(
            &joint.generator_per_s,
            1,
            &scales,
            [1., 100.],
            controls(),
            &mut || {
                calls += 1;
                calls >= 20
            }
        ),
        Err(AssemblyResonanceError::Cancelled)
    ));
    assert!(load.pole_states().iter().all(|z| *z == C64::ZERO));
}
