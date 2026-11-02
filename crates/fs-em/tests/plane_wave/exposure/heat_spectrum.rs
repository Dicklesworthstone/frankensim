//! G1/G3/G4 native film heat spectra, source identity and finite-time loads.
use super::*;
use fs_em::plane_wave::angular_spectrum::{
    AngularSpectrumExposureHeatChannel as Channel, AngularSpectrumExposureHeatSpectrumLimits,
    AngularSpectrumFactorStatistics as Statistics, AngularSpectrumIlluminationLimits,
    PreparedAngularSpectrumIllumination,
};
use fs_fft::harmonic_window::HarmonicWindowProfile;

fn spectrum_limits() -> AngularSpectrumExposureHeatSpectrumLimits {
    AngularSpectrumExposureHeatSpectrumLimits {
        maximum_products: 4096,
        maximum_generated_terms: 8192,
        maximum_workspace_bytes: 32 << 20,
    }
}
fn illuminate<'b, 's, 'm>(
    left: &'b PlanarAngularSpectrum<'s, 'm>,
    right: [C64; 2],
) -> PreparedAngularSpectrumIllumination<'b, 's, 'm> {
    left.prepare_illumination(
        &[right; 4],
        AngularSpectrumIlluminationLimits {
            maximum_layer_mode_pairs: 4096,
            maximum_workspace_bytes: 8 << 20,
            maximum_fft_work: 1 << 20,
        },
        || false,
    )
    .unwrap()
}
fn instant_heat(
    fields: &[AngularSpectrumLayerField],
    omega: &[f64],
    time: f64,
    coefficient: C64,
) -> [f64; 3] {
    let mut e = [0.0; 3];
    let mut pd = [0.0; 3];
    for (field, &omega) in fields.iter().zip(omega) {
        let (s, c) = (omega * time).sin_cos();
        let real = |a: C64| a.re * c - a.im * s;
        for axis in 0..3 {
            let a = field.field().electric_v_per_m()[0][axis] * coefficient;
            let p = a.scale(0.7) / C64::new(1.0, omega * 0.4);
            e[axis] += real(a);
            pd[axis] += real(C64::new(-p.im, p.re).scale(omega));
        }
    }
    let ohmic = 0.3 * e.iter().map(|x| x * x).sum::<f64>();
    let polarization = (0.4 / 0.7) * pd.iter().map(|x| x * x).sum::<f64>();
    [ohmic, polarization, ohmic + polarization]
}

#[test]
fn g1_original_internal_heat_spectra_match_physical_fields_and_shutter_observations() {
    let medium = material(1.0, 1.0);
    let (dielectric, conductor) = heat_material();
    let response = DielectricResponse::DebyeWithOhmic {
        dielectric: &dielectric,
        conductor: &conductor,
    };
    let omega = [2.0, 2.3];
    let stacks = omega.map(|omega| {
        PlanarDielectricLayers::try_new(
            &medium,
            DielectricResponse::Lossless(&medium),
            omega,
            &[DielectricLayer {
                response,
                thickness_m: 0.31,
            }],
        )
        .unwrap()
    });
    let lines = [
        constant(&stacks[0], [C64::new(0.8, 0.2), C64::ZERO]),
        constant(&stacks[1], [C64::new(-0.3, 0.6), C64::ZERO]),
    ];
    let exposure = AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap();
    let layer = exposure.prepare_layer(0, heat_limits(), || false).unwrap();
    let spectra = [Channel::Ohmic, Channel::Polarization, Channel::Total].map(|channel| {
        layer
            .heat_spectrum(0.155, 0, channel, spectrum_limits(), || false)
            .unwrap()
    });
    let pole = layer
        .heat_spectrum(
            0.155,
            0,
            Channel::DebyePole { index: 0 },
            spectrum_limits(),
            || false,
        )
        .unwrap();
    let native = lines
        .iter()
        .map(|line| {
            line.prepare_layer(
                0,
                AngularSpectrumLayerLimits {
                    maximum_layer_mode_pairs: 4096,
                    maximum_workspace_bytes: 16 << 20,
                    maximum_fft_work: 1 << 20,
                },
                || false,
            )
            .unwrap()
            .field(0.155, || false)
            .unwrap()
        })
        .collect::<Vec<_>>();
    for time in [-0.3, 0.0, 0.17, 0.49] {
        let expected = instant_heat(&native, &omega, time, C64::ONE);
        for (spectrum, expected) in spectra.iter().zip(expected) {
            near(spectrum.evaluate(time).unwrap(), expected, 3e-12);
        }
        near(pole.evaluate(time).unwrap(), expected[1], 3e-12);
    }
    assert!(core::ptr::eq(spectra[2].layer(), &layer));
    assert_eq!(spectra[2].grid(), grid());
    assert_eq!(spectra[2].depth_m(), 0.155);
    assert_eq!(spectra[2].point_index(), 0);
    assert_eq!(spectra[2].channel(), Channel::Total);
    for profile in [
        HarmonicWindowProfile::Uniform,
        HarmonicWindowProfile::Triangular,
    ] {
        let window = HarmonicWindow::with_profile(0.17, 0.12, profile).unwrap();
        let original = layer.field(0.155, window, || false).unwrap();
        let q = original.dissipation()[0];
        for (spectrum, expected) in spectra.iter().zip([
            q.ohmic_heat_w_m3,
            q.polarization_heat_w_m3,
            q.total_heat_w_m3,
        ]) {
            near(spectrum.observe_window(window).unwrap(), expected, 3e-12);
        }
        // This fixture returns reversible stored energy while still heating.
        assert!(q.total_heat_w_m3 > 0.0 && q.storage_rate_w_m3 < -1.0);
        assert!(
            (spectra[2].observe_window(window).unwrap()
                - (q.total_heat_w_m3 + q.storage_rate_w_m3))
                .abs()
                > 1.0
        );
    }
}

#[test]
fn g3_joint_face_mixed_heat_spectra_match_the_exact_eight_state_ensemble() {
    let medium = material(1.0, 1.0);
    let (dielectric, conductor) = heat_material();
    let response = DielectricResponse::DebyeWithOhmic {
        dielectric: &dielectric,
        conductor: &conductor,
    };
    let omega = [2.0, 2.3];
    let stacks = omega.map(|omega| {
        PlanarDielectricLayers::try_new(
            &medium,
            DielectricResponse::Lossless(&medium),
            omega,
            &[DielectricLayer {
                response,
                thickness_m: 0.31,
            }],
        )
        .unwrap()
    });
    let left = [
        constant(&stacks[0], [C64::new(0.8, 0.2), C64::new(-0.1, 0.3)]),
        constant(&stacks[1], [C64::new(-0.3, 0.6), C64::new(0.4, 0.1)]),
    ];
    let lines = [
        illuminate(&left[0], [C64::new(-0.2, 0.4), C64::new(0.3, -0.2)]),
        illuminate(&left[1], [C64::new(0.7, -0.1), C64::new(-0.1, 0.2)]),
    ];
    let mean = AngularSpectrumExposure::try_new_illumination(&left[0], &lines, limits(), || false)
        .unwrap();
    let factors = [
        AngularSpectrumExposure::try_new_illumination(&left[0], &lines, limits(), || false)
            .unwrap(),
        AngularSpectrumExposure::try_new_illumination(&left[0], &lines, limits(), || false)
            .unwrap(),
    ];
    let statistics = [Statistics::Real, Statistics::ProperComplex];
    let ensemble = AngularSpectrumExposureEnsemble::try_new_with_statistics(
        &left[0],
        Some(&mean),
        &factors,
        &statistics,
        ensemble_limits(),
        || false,
    )
    .unwrap();
    let layer = ensemble.prepare_layer(0, heat_limits(), || false).unwrap();
    assert_eq!(layer.preparation_layer_mode_pairs(), 0);
    let spectra = [Channel::Ohmic, Channel::Polarization, Channel::Total].map(|channel| {
        layer
            .heat_spectrum(0.155, 0, channel, spectrum_limits(), || false)
            .unwrap()
    });
    let native = lines
        .iter()
        .map(|line| line.layer_field(0, 0.155, || false).unwrap())
        .collect::<Vec<_>>();
    for time in [-0.17, 0.0, 0.23, 0.49] {
        let mut expected = [0.0; 3];
        for sign in [-1.0, 1.0] {
            for phase in [C64::ONE, -C64::ONE, C64::new(0.0, 1.0), C64::new(0.0, -1.0)] {
                let coefficient = C64::from_re(1.0 + sign) + phase;
                for (sum, value) in
                    expected
                        .iter_mut()
                        .zip(instant_heat(&native, &omega, time, coefficient))
                {
                    *sum += value / 8.0;
                }
            }
        }
        for (spectrum, expected) in spectra.iter().zip(expected) {
            near(spectrum.evaluate(time).unwrap(), expected, 5e-12);
        }
    }
    for profile in [
        HarmonicWindowProfile::Uniform,
        HarmonicWindowProfile::Triangular,
    ] {
        let window = HarmonicWindow::with_profile(0.17, 0.12, profile).unwrap();
        let original = layer.field(0.155, window, || false).unwrap();
        near(
            spectra[2].observe_window(window).unwrap(),
            original.dissipation()[0].total_heat_w_m3,
            5e-12,
        );
    }
}

#[test]
fn g4_coherent_residual_heat_keeps_original_owner_after_resource_and_cancellation_refusals() {
    let medium = material(1.0, 1.0);
    let (dielectric, conductor) = heat_material();
    let stack = PlanarDielectricLayers::try_new(
        &medium,
        DielectricResponse::Lossless(&medium),
        2.0,
        &[DielectricLayer {
            response: DielectricResponse::DebyeWithOhmic {
                dielectric: &dielectric,
                conductor: &conductor,
            },
            thickness_m: 0.31,
        }],
    )
    .unwrap();
    let lines = [1e16, 1.0, -1e16].map(|a| constant(&stack, [C64::from_re(a), C64::ZERO]));
    let bank = AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap();
    let layer = bank.prepare_layer(0, heat_limits(), || false).unwrap();
    let spectrum = layer
        .heat_spectrum(0.155, 0, Channel::Total, spectrum_limits(), || false)
        .unwrap();
    let window = HarmonicWindow::new(0.17, 0.12).unwrap();
    let before = layer.field(0.155, window, || false).unwrap();
    near(
        spectrum.observe_window(window).unwrap(),
        before.dissipation()[0].total_heat_w_m3,
        3e-12,
    );
    let exact = AngularSpectrumExposureHeatSpectrumLimits {
        maximum_products: 6,
        maximum_generated_terms: 12,
        maximum_workspace_bytes: spectrum.preparation_workspace_bytes(),
    };
    let exact_spectrum = layer
        .heat_spectrum(0.155, 0, Channel::Total, exact, || false)
        .unwrap();
    assert_eq!(
        spectrum.evaluate(0.17).unwrap(),
        exact_spectrum.evaluate(0.17).unwrap()
    );
    for limits in [
        AngularSpectrumExposureHeatSpectrumLimits {
            maximum_products: 5,
            ..exact
        },
        AngularSpectrumExposureHeatSpectrumLimits {
            maximum_generated_terms: 11,
            ..exact
        },
        AngularSpectrumExposureHeatSpectrumLimits {
            maximum_workspace_bytes: exact.maximum_workspace_bytes - 1,
            ..exact
        },
    ] {
        assert!(
            layer
                .heat_spectrum(0.155, 0, Channel::Total, limits, || false)
                .is_err()
        );
    }
    for (depth, point, channel) in [
        (-0.1, 0, Channel::Total),
        (0.32, 0, Channel::Total),
        (0.155, 4, Channel::Total),
        (0.155, 0, Channel::DebyePole { index: 1 }),
    ] {
        assert!(
            layer
                .heat_spectrum(depth, point, channel, spectrum_limits(), || false)
                .is_err()
        );
    }
    let mut polls = 0;
    assert!(matches!(
        layer.heat_spectrum(0.155, 0, Channel::Total, spectrum_limits(), || {
            polls += 1;
            polls == 7
        }),
        Err(AngularSpectrumError::Cancelled)
    ));
    assert!(matches!(
        spectrum.evaluate_cancellable(0.17, || true),
        Err(AngularSpectrumError::Cancelled)
    ));
    assert!(matches!(
        spectrum.observe_window_cancellable(window, || true),
        Err(AngularSpectrumError::Cancelled)
    ));
    assert!(spectrum.evaluate(f64::NAN).is_err());
    assert_eq!(
        before.dissipation(),
        layer.field(0.155, window, || false).unwrap().dissipation()
    );
}
