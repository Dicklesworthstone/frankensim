//! G1 spectral Fourier references and G3 passive narrowband interference.

use super::*;
use fs_render::camera::{
    CoherentDiffractionInterferometer, CoherentJonesArm, CoherentJonesPupilLayer,
    InterferometerTemporalCoherence,
};

fn receiver(powers: [f64; 2], phase: f64) -> CoherentDiffractionInterferometer {
    CoherentDiffractionInterferometer::try_new([[0.0; 2]; 2], 0.03, powers, 0.0, phase).unwrap()
}

#[test]
fn g1_line_profiles_match_independent_spectral_fourier_integrals() {
    let camera = focused();
    let waves = [wave([0.0; 2], C64::ONE)];
    let query = camera
        .prepare_coherent_diffraction(
            LAMBDA_NM,
            CoherentPupilIllumination {
                mean: &waves,
                factors: &[],
            },
            LIMITS,
            || false,
        )
        .unwrap();
    let tau = core::f64::consts::TAU;
    for q in [0.0, 0.2, 0.7, 2.0] {
        // Independent normalized POWER spectra in x=(nu-nu0)/FWHM. Integrate
        // their Fourier transforms rather than reusing the production envelope.
        let gaussian = integrate_reference(-6.0, 6.0, |x| {
            2.0 * (core::f64::consts::LN_2 / core::f64::consts::PI).sqrt()
                * (-4.0 * core::f64::consts::LN_2 * x * x).exp()
                * (tau * x * q).cos()
        });
        // Lorentzian truncation leaves an absolute tail below 1/(pi*100).
        let lorentzian = integrate_reference(-100.0, 100.0, |x| {
            (tau * x * q).cos() / (tau * (x * x + 0.25))
        });
        for (coherence, expected, tolerance) in [
            (
                InterferometerTemporalCoherence::try_gaussian(1e6, q / 1e6).unwrap(),
                gaussian,
                2e-12,
            ),
            (
                InterferometerTemporalCoherence::try_lorentzian(1e6, q / 1e6).unwrap(),
                lorentzian,
                0.0032,
            ),
        ] {
            let measured = query
                .observe_interferometer_with_linewidth(
                    receiver([0.5; 2], 0.0),
                    coherence,
                    2,
                    || false,
                )
                .unwrap();
            assert!((measured.fringe_visibility - expected).abs() < tolerance);
            let total = measured.arm_power_density_w.iter().sum::<f64>();
            close(
                measured.output_power_density_w.iter().sum(),
                total,
                3e-13,
            );
            close(
                measured.mutual_coherence_density_w.re / total,
                0.5 * measured.fringe_visibility,
                3e-13,
            );
        }
    }
}

#[test]
fn g1_linewidth_keeps_phase_energy_unequal_arms_and_weak_dark_power() {
    let camera = focused();
    let waves = [wave([0.0; 2], C64::ONE)];
    let query = camera
        .prepare_coherent_diffraction(
            LAMBDA_NM,
            CoherentPupilIllumination {
                mean: &waves,
                factors: &[],
            },
            LIMITS,
            || false,
        )
        .unwrap();
    let peak = query.observe([0.0; 2], || false).unwrap().power_density_w;
    let contrast = 2.0 * (0.25_f64 * 0.75).sqrt();
    for delay in [-0.4e-6_f64, 0.4e-6] {
        let retained = (-core::f64::consts::PI * 1e6 * delay.abs()).exp();
        let coherence = InterferometerTemporalCoherence::try_lorentzian(1e6, delay).unwrap();
        for phase in [0.0_f64, 0.7, core::f64::consts::PI] {
            let measured = query
                .observe_interferometer_with_linewidth(
                    receiver([0.25, 0.75], phase),
                    coherence,
                    2,
                    || false,
                )
                .unwrap();
            close(measured.fringe_visibility, contrast * retained, 3e-13);
            close(
                measured.output_power_density_w[0] / peak,
                0.5 + 0.5 * contrast * retained * phase.cos(),
                3e-13,
            );
            close(
                measured.output_power_density_w[1] / peak,
                0.5 - 0.5 * contrast * retained * phase.cos(),
                3e-13,
            );
            close(
                measured.output_power_density_w.iter().sum(),
                peak,
                3e-13,
            );
        }
    }
    for coherence in [
        InterferometerTemporalCoherence::try_lorentzian(0.0, 1e300).unwrap(),
        InterferometerTemporalCoherence::try_gaussian(1e6, 0.0).unwrap(),
    ] {
        assert_eq!(
            query
                .observe_interferometer_with_linewidth(
                    receiver([0.5; 2], 0.2),
                    coherence,
                    2,
                    || false,
                )
                .unwrap(),
            query
                .observe_interferometer(receiver([0.5; 2], 0.2), 2, || false)
                .unwrap(),
        );
    }

    // The coherence rounds to one in binary64, but the dark-port loss remains
    // positive and agrees with the first nonzero analytic linewidth term.
    let x = core::f64::consts::PI * 1e-20;
    let weak = query
        .observe_interferometer_with_linewidth(
            receiver([0.5; 2], 0.0),
            InterferometerTemporalCoherence::try_lorentzian(1e-8, 1e-12).unwrap(),
            2,
            || false,
        )
        .unwrap();
    assert!(weak.output_power_density_w[1] > 0.0);
    close(weak.output_power_density_w[1] / (0.5 * peak * x), 1.0, 3e-13);

    let incoherent = query
        .observe_interferometer_with_linewidth(
            receiver([0.5; 2], 0.0),
            InterferometerTemporalCoherence::try_gaussian(1e6, 1e308).unwrap(),
            2,
            || false,
        )
        .unwrap();
    assert_eq!(incoherent.fringe_visibility, 0.0);
    for port in incoherent.output_power_density_w {
        close(port, 0.5 * peak, 3e-13);
    }
}

#[test]
fn g1_vacuum_path_mismatch_retains_carrier_fringe_and_loses_temporal_visibility() {
    let camera = focused();
    let waves = [wave([0.0; 2], C64::ONE)];
    let query = camera
        .prepare_coherent_diffraction(
            512.0,
            CoherentPupilIllumination {
                mean: &waves,
                factors: &[],
            },
            LIMITS,
            || false,
        )
        .unwrap();
    // An integer metre is exactly 1,953,125 vacuum wavelengths at 512 nm.
    // The carrier is constructive while finite linewidth washes out over L/c.
    for path in [0.0, 1.0, 100.0, 1000.0] {
        let detector =
            CoherentDiffractionInterferometer::try_new([[0.0; 2]; 2], 0.03, [0.5; 2], path, 0.0)
                .unwrap();
        let delay = path / 299_792_458.0;
        let expected = (-core::f64::consts::PI * 1e6 * delay).exp();
        let measured = query
            .observe_interferometer_with_linewidth(
                detector,
                InterferometerTemporalCoherence::try_lorentzian(1e6, delay).unwrap(),
                2,
                || false,
            )
            .unwrap();
        close(measured.fringe_visibility, expected, 4e-13);
        let total = measured.arm_power_density_w.iter().sum::<f64>();
        close(
            measured.output_power_density_w[0] / total,
            0.5 * (1.0 + expected),
            3e-13,
        );
        assert_eq!(
            query
                .observe_interferometer(detector, 2, || false)
                .unwrap()
                .output_power_density_w[1],
            0.0,
        );
    }
}

#[test]
fn g3_linewidth_preserves_independent_spatial_factors_and_polarized_coherencies() {
    let camera = focused();
    let x = [jones_wave([0.0; 2], [C64::ONE, C64::ZERO])];
    let y = [jones_wave([0.0; 2], [C64::ZERO, C64::ONE])];
    let together = [jones_wave([0.0; 2], [C64::ONE; 2])];
    let factors = [
        JonesPupilCoherenceFactor { waves: &x },
        JonesPupilCoherenceFactor { waves: &y },
    ];
    let arms = [[1.0, 0.0], [0.0, 1.0]].map(|powers| {
        CoherentJonesArm::try_new(&[
            CoherentJonesPupilLayer::try_new(powers, [0.0; 2], 0.0, [0.0; 2]).unwrap(),
        ])
        .unwrap()
    });
    let coherence = InterferometerTemporalCoherence::try_gaussian(1e6, 0.4e-6).unwrap();
    let retained =
        (-(core::f64::consts::PI * 0.4).powi(2) / (4.0 * core::f64::consts::LN_2)).exp();
    let phase = 0.3_f64;
    for illumination in [
        CoherentJonesPupilIllumination {
            mean: &together,
            factors: &[],
        },
        CoherentJonesPupilIllumination {
            mean: &[],
            factors: &factors,
        },
    ] {
        let query = camera
            .prepare_polarized_coherent_diffraction(LAMBDA_NM, illumination, LIMITS, || false)
            .unwrap();
        let measured = query
            .observe_interferometer_with_linewidth(
                receiver([0.5; 2], phase),
                arms,
                coherence,
                2 * query.wave_count(),
                || false,
            )
            .unwrap();
        let total = measured.arm_power_density_w.iter().sum::<f64>();
        assert_eq!(measured.fringe_visibility, 0.0);
        for (port, output) in measured.output_coherency_density_w.into_iter().enumerate() {
            let stokes = output.stokes();
            close(stokes.i, 0.5 * total, 3e-13);
            assert!(stokes.q.abs() <= 3e-13 * total);
            if illumination.mean.is_empty() {
                assert_eq!(stokes.u, 0.0);
                assert_eq!(stokes.v, 0.0);
            } else {
                let sign = if port == 0 { 1.0 } else { -1.0 };
                close(stokes.u, sign * 0.5 * total * retained * phase.cos(), 3e-13);
                close(stokes.v, -sign * 0.5 * total * retained * phase.sin(), 3e-13);
            }
        }
    }
}

#[test]
fn g1_finite_pixels_resolve_localized_temporal_coherence_before_integration() {
    let camera =
        base(Aperture::try_regular_polygon(RADIUS, 4, core::f64::consts::FRAC_PI_4).unwrap())
            .with_in_focus_diffraction(FocusedDiffraction::try_new(0.0).unwrap())
            .unwrap();
    let waves = [wave([0.0; 2], C64::ONE)];
    let query = camera
        .prepare_coherent_diffraction(
            LAMBDA_NM,
            CoherentPupilIllumination {
                mean: &waves,
                factors: &[],
            },
            LIMITS,
            || false,
        )
        .unwrap();
    let coherence = InterferometerTemporalCoherence::try_gaussian(1e6, -3.0 * 0.37 / 1e6)
        .unwrap()
        .try_with_relay_delay_gradient([3.0 / (1e6 * SCALE), 0.0])
        .unwrap();
    let footprint = pixel([-3.0, -1.0], [3.0, 1.0]);
    let profile = |q: f64| sinc_reference(q * core::f64::consts::FRAC_1_SQRT_2).powi(2);
    let y = integrate_reference(-1.0, 1.0, profile) / (2.0 * core::f64::consts::PI.powi(2));
    let power = y * integrate_reference(-3.0, 3.0, profile);
    let correlated = y
        * integrate_reference(-3.0, 3.0, |q| {
            profile(q)
                * (-(core::f64::consts::PI * 3.0 * (q - 0.37)).powi(2)
                    / (4.0 * core::f64::consts::LN_2))
                    .exp()
        });
    let bare = query
        .integrate_interferometer_pixel(receiver([0.5; 2], 0.0), footprint, PIXEL_LIMITS, || {
            false
        })
        .unwrap();
    for phase_span in [1.0, 0.5] {
        let measured = query
            .integrate_interferometer_pixel_with_linewidth(
                receiver([0.5; 2], 0.0),
                coherence,
                footprint,
                CoherentDiffractionPixelLimits {
                    maximum_phase_span_rad: phase_span,
                    ..PIXEL_LIMITS
                },
                || false,
            )
            .unwrap();
        assert!(measured.subdivisions[0] > bare.subdivisions[0]);
        close(measured.output_power_w[0], 0.5 * (power + correlated), 3e-10);
        close(measured.output_power_w[1], 0.5 * (power - correlated), 3e-10);
        close(measured.fringe_visibility, correlated / power, 3e-10);
        close(measured.mutual_coherence_w.re, 0.5 * correlated, 3e-10);
        assert_eq!(measured.pupil_evaluations, 2 * measured.quadrature_samples);
    }
}

#[test]
fn g0_g4_linewidth_admission_budget_and_cancellation_leave_original_query_reusable() {
    for invalid in [-1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            InterferometerTemporalCoherence::try_lorentzian(invalid, 0.0),
            Err(CoherentDiffractionError::InvalidTemporalCoherence),
        );
        assert_eq!(
            InterferometerTemporalCoherence::try_gaussian(invalid, 0.0),
            Err(CoherentDiffractionError::InvalidTemporalCoherence),
        );
    }
    assert_eq!(
        InterferometerTemporalCoherence::try_gaussian(1e6, f64::INFINITY),
        Err(CoherentDiffractionError::InvalidTemporalCoherence),
    );
    let coherence = InterferometerTemporalCoherence::try_lorentzian(1e6, 1e-6).unwrap();
    assert_eq!(
        coherence.try_with_relay_delay_gradient([f64::NAN, 0.0]),
        Err(CoherentDiffractionError::InvalidTemporalCoherence),
    );
    let camera = focused();
    let waves = [wave([0.0; 2], C64::ONE)];
    let query = camera
        .prepare_coherent_diffraction(
            LAMBDA_NM,
            CoherentPupilIllumination {
                mean: &waves,
                factors: &[],
            },
            LIMITS,
            || false,
        )
        .unwrap();
    let detector = receiver([0.5; 2], 0.0);
    let expected = query.observe_interferometer(detector, 2, || false).unwrap();
    assert_eq!(
        query.observe_interferometer_with_linewidth(detector, coherence, 2, || true),
        Err(CoherentDiffractionError::Cancelled),
    );
    assert!(matches!(
        query.observe_interferometer_with_linewidth(detector, coherence, 1, || false),
        Err(CoherentDiffractionError::Limit { .. }),
    ));
    let broad = InterferometerTemporalCoherence::try_gaussian(
        0.02 * 299_792_458.0 / LAMBDA,
        0.0,
    )
    .unwrap();
    assert_eq!(
        query.observe_interferometer_with_linewidth(detector, broad, 2, || false),
        Err(CoherentDiffractionError::InvalidTemporalCoherence),
    );
    assert!(matches!(
        query.integrate_interferometer_pixel_with_linewidth(
            detector,
            coherence,
            pixel([-1.0; 2], [1.0; 2]),
            CoherentDiffractionPixelLimits {
                maximum_samples: 1,
                ..PIXEL_LIMITS
            },
            || false,
        ),
        Err(CoherentDiffractionError::Limit { .. }),
    ));
    assert_eq!(
        query.observe_interferometer(detector, 2, || false).unwrap(),
        expected,
    );
}
