use super::*;

fn sweep(stack: &PlanarDielectricLayers<'_>, p: f64, polarization: Polarization) {
    let ki = stack
        .response(1.0, polarization)
        .unwrap()
        .incident_normal_wavenumber_rad_m;
    let kx = p * ki;
    let jet = stack.transverse_dispersion(kx, polarization).unwrap();
    assert!(std::ptr::eq(jet.source(), stack));
    assert_eq!(jet.transverse_wavenumber_rad_m(), kx);
    assert_eq!(
        jet.carrier(),
        stack
            .response(jet.incidence_cosine(), polarization)
            .unwrap()
    );
    let step = 1e-4 * ki;
    let samples: [[C64; 3]; 5] = core::array::from_fn(|index| {
        let transverse = (kx + (index as f64 - 2.0) * step) / ki;
        let cosine = det::sqrt((1.0 - transverse.abs()) * (1.0 + transverse.abs()));
        let value = stack.response(cosine, polarization).unwrap();
        [
            value.reflected_tangential_e,
            value.transmitted_tangential_e,
            value.transmitted_conjugate_h,
        ]
    });
    for (index, amplitude) in [
        jet.reflected_tangential_e(),
        jet.transmitted_tangential_e().unwrap(),
        jet.transmitted_conjugate_h().unwrap(),
    ]
    .into_iter()
    .enumerate()
    {
        let [a, b, c, d, e] = samples.map(|value| value[index]);
        let first = (a - b.scale(8.0) + d.scale(8.0) - e).scale(1.0 / (12.0 * step));
        let second = (-a + b.scale(16.0) - c.scale(30.0) + d.scale(16.0) - e)
            .scale(1.0 / (12.0 * step * step));
        complex_near(amplitude.derivative_m(), first, 4e-8 / ki);
        complex_near(amplitude.second_derivative_m2(), second, 4e-6 / ki / ki);
    }
}

#[test]
fn g1_transverse_tir_recovers_signed_artmann_shift_and_curvature() {
    let incident = material(1.0, 1.0);
    let substrate = material(0.25, 1.0);
    for ki in [0.01, 1.0, 1e9] {
        let stack = PlanarDielectricLayers::try_new(
            &incident,
            DielectricResponse::Lossless(&substrate),
            ki,
            &[],
        )
        .unwrap();
        for p in [0.55_f64, 0.75, 0.95] {
            let c = (1.0 - p * p).sqrt();
            let gamma = (p * p - 0.25).sqrt();
            let first = 2.0 * p / (c * gamma);
            let second = 2.0 / (c * gamma) + 2.0 * p * p / (c.powi(3) * gamma)
                - 2.0 * p * p / (c * gamma.powi(3));
            for polarization in [Polarization::Te, Polarization::Tm] {
                let plus = stack.transverse_dispersion(ki * p, polarization).unwrap();
                let minus = stack.transverse_dispersion(-ki * p, polarization).unwrap();
                let a = plus.reflected_tangential_e();
                let b = minus.reflected_tangential_e();
                assert_eq!(a.amplitude(), b.amplitude());
                complex_near(a.derivative_m(), -b.derivative_m(), 4e-14 / ki);
                complex_near(
                    a.second_derivative_m2(),
                    b.second_derivative_m2(),
                    4e-13 / ki / ki,
                );
                let (slope, curvature) = if polarization == Polarization::Te {
                    (first, second)
                } else {
                    let denominator = 1.25 * p * p - 0.25;
                    let factor = 0.25 / denominator;
                    (
                        first * factor,
                        second * factor - first * 0.625 * p / denominator.powi(2),
                    )
                };
                let phase = a.logarithmic_dispersion().unwrap();
                near(phase.lateral_shift_m, slope / ki, 8e-13 * slope.abs() / ki);
                near(
                    phase.phase_curvature_m2,
                    curvature / ki / ki,
                    2e-12 * curvature.abs().max(1.0) / ki / ki,
                );
                near(phase.log_amplitude_derivative_m, 0.0, 8e-13 / ki);
                near(
                    phase.log_amplitude_second_derivative_m2,
                    0.0,
                    8e-12 / ki / ki,
                );
            }
        }
    }
}

#[test]
fn g1_transverse_lossy_coatings_and_internal_critical_film_follow_native_sweeps() {
    let incident = material(E, M);
    let film = material(4.0 * E, 1.2 * M);
    let relaxing = crate::finite_layers::debye(2.0 * E, 5.0 * E, 0.01, 0.001);
    let layers = [
        DielectricLayer {
            response: DielectricResponse::Debye(&relaxing),
            thickness_m: 0.17,
        },
        DielectricLayer {
            response: DielectricResponse::Lossless(&film),
            thickness_m: 0.22,
        },
        DielectricLayer {
            response: DielectricResponse::Debye(&relaxing),
            thickness_m: 0.09,
        },
    ];
    let stack = PlanarDielectricLayers::try_new(
        &incident,
        DielectricResponse::Debye(&relaxing),
        0.7e9,
        &layers,
    )
    .unwrap();
    for p in [-0.8, -0.3, 0.0, 0.6] {
        for polarization in [Polarization::Te, Polarization::Tm] {
            sweep(&stack, p, polarization);
        }
    }

    let incident = material(1.0, 1.0);
    // The admitted c=.8 gives exactly critical n²=(1-c)(1+c), including
    // the native material-product residual. The p sweep crosses both branches.
    let critical = material(1.0 - 0.8, 1.0 + 0.8);
    let substrate = material(2.25, 1.0);
    let stack = PlanarDielectricLayers::try_new(
        &incident,
        DielectricResponse::Lossless(&substrate),
        1.0,
        &[DielectricLayer {
            response: DielectricResponse::Lossless(&critical),
            thickness_m: 2.3,
        }],
    )
    .unwrap();
    for polarization in [Polarization::Te, Polarization::Tm] {
        sweep(&stack, 0.6, polarization);
    }
}

#[test]
fn g1_transverse_homogeneous_film_keeps_geometrical_exit_displacement() {
    let incident = material(1.0, 1.0);
    let depth = 3.7;
    let ki = 2.3;
    let layers = [DielectricLayer {
        response: DielectricResponse::Lossless(&incident),
        thickness_m: depth,
    }];
    let stack = PlanarDielectricLayers::try_new(
        &incident,
        DielectricResponse::Lossless(&incident),
        ki,
        &layers,
    )
    .unwrap();
    for p in [-0.7_f64, 0.0, 0.7] {
        let c = (1.0 - p * p).sqrt();
        for polarization in [Polarization::Te, Polarization::Tm] {
            let jet = stack.transverse_dispersion(ki * p, polarization).unwrap();
            let t = jet.transmitted_tangential_e().unwrap();
            complex_near(t.amplitude(), phase(-ki * c * depth), 8e-15);
            let phase = t.logarithmic_dispersion().unwrap();
            near(phase.lateral_shift_m, depth * p / c, 8e-14);
            near(phase.phase_curvature_m2, depth / (ki * c.powi(3)), 3e-13);
            near(phase.log_amplitude_derivative_m, 0.0, 2e-14);
            for derivative in jet.reflected_tangential_e().scaled_derivatives() {
                complex_near(derivative, C64::ZERO, 2e-13);
            }
        }
    }
}

#[test]
fn g0_transverse_zero_critical_exit_and_opaque_limits_remain_explicit() {
    let incident = material(1.0, 1.0);
    let matched = PlanarDielectricLayers::try_new(
        &incident,
        DielectricResponse::Lossless(&incident),
        1.0,
        &[],
    )
    .unwrap();
    let zero = matched
        .transverse_dispersion(0.6, Polarization::Te)
        .unwrap()
        .reflected_tangential_e();
    assert_eq!(zero.amplitude(), C64::ZERO);
    assert!(zero.logarithmic_dispersion().is_none());
    for k in [-1.0, 1.0, 1.1, f64::NAN, f64::INFINITY] {
        assert!(matched.transverse_dispersion(k, Polarization::Te).is_err());
    }
    // Native binary64 Brewster pair at c=.8: q=.5=epsilon_ratio*c.
    let substrate = material(0.625, 0.976);
    let brewster = PlanarDielectricLayers::try_new(
        &incident,
        DielectricResponse::Lossless(&substrate),
        1.0,
        &[],
    )
    .unwrap();
    let zero = brewster
        .transverse_dispersion(0.6, Polarization::Tm)
        .unwrap()
        .reflected_tangential_e();
    assert_eq!(zero.amplitude(), C64::ZERO);
    assert!(zero.scaled_derivatives()[0].abs() > 0.1);
    assert!(zero.logarithmic_dispersion().is_none());

    let critical = material(1.0 - 0.8, 1.0 + 0.8);
    let exit = PlanarDielectricLayers::try_new(
        &incident,
        DielectricResponse::Lossless(&critical),
        1.0,
        &[],
    )
    .unwrap();
    assert_eq!(
        exit.response(0.8, Polarization::Te)
            .unwrap()
            .exit_normal_wavenumber_rad_m,
        C64::ZERO
    );
    assert!(exit.transverse_dispersion(0.6, Polarization::Te).is_err());

    let incident = material(E, M);
    let lossy = crate::finite_layers::debye(2.0 * E, 5.0 * E, 0.01, 0.001);
    let opaque = PlanarDielectricLayers::try_new(
        &incident,
        DielectricResponse::Lossless(&incident),
        1e9,
        &[DielectricLayer {
            response: DielectricResponse::Debye(&lossy),
            thickness_m: 1e5,
        }],
    )
    .unwrap();
    let front =
        PlanarDielectricLayers::try_new(&incident, DielectricResponse::Debye(&lossy), 1e9, &[])
            .unwrap();
    for polarization in [Polarization::Te, Polarization::Tm] {
        let a = opaque.transverse_dispersion(1.0, polarization).unwrap();
        let b = front.transverse_dispersion(1.0, polarization).unwrap();
        assert!(a.transmitted_tangential_e().is_none());
        assert!(a.transmitted_conjugate_h().is_none());
        for (a, b) in a
            .reflected_tangential_e()
            .scaled_derivatives()
            .into_iter()
            .zip(b.reflected_tangential_e().scaled_derivatives())
        {
            complex_near(a, b, 2e-13);
        }
    }
}

#[test]
fn g1_existing_finite_angular_spectrum_beam_converges_to_the_local_tir_shift() {
    let incident = material(1.0, 1.0);
    let substrate = material(0.25, 1.0);
    let stack = PlanarDielectricLayers::try_new(
        &incident,
        DielectricResponse::Lossless(&substrate),
        1.0,
        &[],
    )
    .unwrap();
    let carrier = 0.75;
    let predicted = stack
        .transverse_dispersion(carrier, Polarization::Te)
        .unwrap()
        .reflected_tangential_e()
        .logarithmic_dispersion()
        .unwrap()
        .lateral_shift_m;
    let n = 1024;
    let dk = 0.001;
    let spacing = core::f64::consts::TAU / (n as f64 * dk);
    let grid = AngularSpectrumGrid {
        width: n,
        height: 1,
        spacing_m: [spacing, 1.0],
    };
    let mut previous_error = f64::INFINITY;
    for sigma in [0.01, 0.005] {
        let mut modes = vec![[C64::ZERO; 2]; n];
        for bin in -100_i32..=100 {
            let offset = bin as f64 * dk;
            let sign = if bin % 2 == 0 { 1.0 } else { -1.0 };
            let value = sign * det::exp(-offset * offset / (4.0 * sigma * sigma));
            modes[bin.rem_euclid(n as i32) as usize][1] = C64::from_re(value);
        }
        let beam = PlanarAngularSpectrum::prepare_modes_with_carrier(
            &stack,
            grid,
            [carrier, 0.0],
            &modes,
            limits(),
            || false,
        )
        .unwrap();
        let centroid = |side| {
            let field = beam.field(side, 0.0, || false).unwrap();
            let mut sum = 0.0;
            let mut moment = 0.0;
            for (index, value) in field.electric_v_per_m().iter().enumerate() {
                // Scalar TE intensity centroid on the entrance plane, not a
                // paraxial ray position or a beam-axis-normal displacement.
                let intensity = value[1].norm_sq();
                sum += intensity;
                moment += index as f64 * spacing * intensity;
            }
            moment / sum
        };
        let displacement = centroid(Side::Reflected) - centroid(Side::Incident);
        let error = (displacement - predicted).abs();
        assert!(displacement > 0.0);
        assert!(error < 0.35 * previous_error);
        assert!(error < 0.004 * predicted);
        previous_error = error;
    }
}
