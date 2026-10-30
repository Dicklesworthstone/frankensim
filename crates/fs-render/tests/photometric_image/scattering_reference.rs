//! Independent airlight integrals and the limits of a refinement discrepancy.
//! These are analytic transport checks, not CIE 171 or fixture measurements.

use super::*;

const FIXTURES: [([f64; 3], f64); 3] = [
    ([0.75, 0.0, 1.25], 120.0),
    ([-2.0, 1.5, 5.0], 75.0),
    ([0.0, -1.3, -2.0], 180.0),
];

fn fixtures(scale: f64) -> Vec<PhotometricSource> {
    FIXTURES
        .iter()
        .enumerate()
        .map(|(index, &(position, intensity))| {
            // A gray medium transports equal candela identically regardless
            // of spectrum. Include narrow red/blue bins as well as white.
            let spectrum = if index == 0 {
                SourceSpectrum::equal_energy()
            } else {
                let mut bins = [0.0; SPECTRAL_BINS];
                bins[if index == 1 { 14 } else { 60 }] = 1.0;
                SourceSpectrum::try_from_bins(bins).unwrap()
            };
            source(
                position.map(|coordinate| coordinate * scale),
                intensity * scale * scale,
                spectrum,
            )
        })
        .collect()
}

// Camera: (0,0,s), 0 <= s <= D. Source: radial offset b, axial position c.
// L = I*sigma_s * integral p((c-s)/r) exp[-sigma_t*(s+r)] / r² ds.
// Substitute theta=atan((s-c)/b): ds/r²=dtheta/b, r=b/cos(theta).
// This independently evaluates that smooth angle-space integral with Simpson
// weights, rather than the renderer's distance-space midpoint/Beer panels.
// No production phase, transmittance, geometry or spectral function supplies
// the expected value. The source's declared intensity is already in candela.
fn airlight_reference(phase: PhaseFunction, absorption: f64, scattering: f64) -> f64 {
    let panels = 4096;
    FIXTURES
        .iter()
        .map(|&(position, intensity)| {
            let b = position[0].hypot(position[1]);
            let c = position[2];
            let low = (-c / b).atan();
            let high = ((4.0 - c) / b).atan();
            let step = (high - low) / f64::from(panels);
            let value = |theta: f64| {
                let cosine = -theta.sin();
                let phase = match phase {
                    PhaseFunction::Isotropic => 1.0 / (4.0 * PI),
                    PhaseFunction::HenyeyGreenstein { g } => {
                        (1.0 - g * g) / (4.0 * PI * (1.0 + g * g - 2.0 * g * cosine).powf(1.5))
                    }
                    PhaseFunction::Rayleigh => 3.0 * (1.0 + cosine * cosine) / (16.0 * PI),
                };
                let s = c + b * theta.tan();
                let r = b / theta.cos();
                phase * (-(absorption + scattering) * (s + r)).exp()
            };
            let mut sum = value(low) + value(high);
            for panel in 1..panels {
                sum +=
                    if panel % 2 == 0 { 2.0 } else { 4.0 } * value(low + f64::from(panel) * step);
            }
            intensity * scattering * sum * step / (3.0 * b)
        })
        .sum()
}

#[test]
fn g1_noncollinear_multi_fixture_airlight_matches_independent_reference() {
    with_cx(|cx, _| {
        let sources = fixtures(1.0);
        for phase in [
            PhaseFunction::Isotropic,
            PhaseFunction::HenyeyGreenstein { g: 0.7 },
            PhaseFunction::HenyeyGreenstein { g: -0.45 },
            PhaseFunction::Rayleigh,
        ] {
            let medium = HomogeneousMedium::try_new(0.04, 0.08, phase).unwrap();
            let scene =
                PhotometricScene::try_new(&[], &sources, medium, TransportOptions::default(), cx)
                    .unwrap();
            let result = scene
                .luminance_along_ray(
                    cx,
                    &forward_ray(),
                    4.0,
                    ScatterSettings {
                        initial_panels: 32,
                        max_panels: 16_384,
                        absolute_tolerance_cd_m2: 0.0,
                        relative_tolerance: 1.0e-8,
                    },
                )
                .unwrap();
            let expected = airlight_reference(phase, 0.04, 0.08);
            close(result.luminance_cd_m2, expected, 5.0e-8);
            assert_eq!(result.surface_cd_m2, 0.0);
            assert!(result.single_scatter.tolerance_met);
            assert!(result.total_error_bound_cd_m2.is_none());
            assert!(
                result
                    .single_scatter
                    .omitted_multiple_scattering_bound_cd_m2
                    .is_none()
            );
            if phase == PhaseFunction::Isotropic {
                // Removing attenuation gives a closed-form upper envelope.
                // Since s+r grows along this ray, its endpoint values bound
                // extinction. These are checks on this single-scatter model,
                // not a bound on the omitted higher scattering orders.
                let (mut lower, mut upper) = (0.0, 0.0);
                for &(position, intensity) in &FIXTURES {
                    let b = position[0].hypot(position[1]);
                    let c = position[2];
                    let unattenuated = intensity * 0.08 / (4.0 * PI * b)
                        * (((4.0 - c) / b).atan() - (-c / b).atan());
                    lower += unattenuated * (-0.12 * (4.0 + b.hypot(4.0 - c))).exp();
                    upper += unattenuated * (-0.12 * b.hypot(c)).exp();
                }
                assert!(result.luminance_cd_m2 >= lower);
                assert!(result.luminance_cd_m2 <= upper);
            }
        }
    });
}

#[test]
fn g3_airlight_preserves_si_scaling_of_lengths_coefficients_and_intensity() {
    with_cx(|cx, _| {
        let mut reference = None;
        for scale in [1.0, 1.0 / 1024.0, 1024.0] {
            let sources = fixtures(scale);
            let medium = HomogeneousMedium::try_new(
                0.04 / scale,
                0.08 / scale,
                PhaseFunction::HenyeyGreenstein { g: 0.7 },
            )
            .unwrap();
            let scene = PhotometricScene::try_new(
                &[],
                &sources,
                medium,
                TransportOptions {
                    ray_epsilon_m: 1.0e-6 * scale,
                    max_distance_m: 32.0 * scale,
                    ..TransportOptions::default()
                },
                cx,
            )
            .unwrap();
            let result = scene
                .luminance_along_ray(
                    cx,
                    &forward_ray(),
                    4.0 * scale,
                    ScatterSettings {
                        initial_panels: 256,
                        max_panels: 256,
                        absolute_tolerance_cd_m2: 0.0,
                        relative_tolerance: 0.0,
                    },
                )
                .unwrap();
            // I grows as length², sigma shrinks as 1/length. Thus irradiance,
            // optical depth and the complete first-order luminance stay fixed.
            close(
                result.luminance_cd_m2,
                *reference.get_or_insert(result.luminance_cd_m2),
                5.0e-13,
            );
            assert!(!result.single_scatter.tolerance_met);
            assert!(
                result
                    .single_scatter
                    .estimated_quadrature_error_cd_m2
                    .is_none()
            );
            assert!(result.total_error_bound_cd_m2.is_none());
        }
    });
}

#[test]
fn g1_zero_refinement_discrepancy_does_not_bound_an_unsampled_beam() {
    with_cx(|cx, _| {
        let half_angle = 0.01_f64;
        let intensity = 100.0;
        // The fixture at (-1,0,1/16) aims along +X. Its narrow cone crosses
        // the camera ray near z=1/16, between every 2- and 4-panel midpoint.
        let q = core::f64::consts::FRAC_1_SQRT_2;
        let sources = [PhotometricSource::try_new(
            RigidTransform::try_new([0.0, q, 0.0, q], [-1.0, 0.0, 1.0 / 16.0]).unwrap(),
            intensity * 2.0 * PI * (1.0 - half_angle.cos()),
            1.0,
            AngularDistribution::uniform_cone(half_angle).unwrap(),
            SourceSpectrum::equal_energy(),
        )
        .unwrap()];
        let medium = HomogeneousMedium::try_new(0.05, 0.1, PhaseFunction::Isotropic).unwrap();
        let scene =
            PhotometricScene::try_new(&[], &sources, medium, TransportOptions::default(), cx)
                .unwrap();
        let settings = ScatterSettings {
            initial_panels: 2,
            max_panels: 4096,
            absolute_tolerance_cd_m2: 0.0,
            relative_tolerance: 0.0,
        };
        let aliased = scene
            .luminance_along_ray(cx, &forward_ray(), 1.0, settings)
            .unwrap();
        assert_eq!(aliased.luminance_cd_m2, 0.0);
        assert_eq!(aliased.single_scatter.panels, 4);
        assert_eq!(
            aliased.single_scatter.estimated_quadrature_error_cd_m2,
            Some(0.0)
        );
        assert!(aliased.single_scatter.tolerance_met);
        assert!(aliased.total_error_bound_cd_m2.is_none());
        assert!(
            aliased
                .single_scatter
                .omitted_multiple_scattering_bound_cd_m2
                .is_none()
        );

        // The actual cone intersection is nonempty. Integrating 1/r² over
        // it gives 2*half_angle (b=1 m); using the longest s+r gives a strictly
        // positive lower envelope, independently of quadrature midpoints.
        let longest_path = 1.0 / 16.0 + half_angle.tan() + 1.0 / half_angle.cos();
        let lower = intensity * 0.1 * half_angle / (2.0 * PI) * (-0.15 * longest_path).exp();
        assert!(lower > 0.01);
        let resolved = scene
            .luminance_along_ray(
                cx,
                &forward_ray(),
                1.0,
                ScatterSettings {
                    initial_panels: 4096,
                    ..settings
                },
            )
            .unwrap();
        assert!(resolved.luminance_cd_m2 > 0.95 * lower);
        assert!(resolved.total_error_bound_cd_m2.is_none());
    });
}
