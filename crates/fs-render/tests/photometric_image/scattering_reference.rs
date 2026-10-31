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
fn airlight_reference(phase: &PhaseFunction, absorption: f64, scattering: f64) -> f64 {
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
                    PhaseFunction::Tabulated(_)
                    | PhaseFunction::RayleighDepolarized { .. }
                    | PhaseFunction::SpectralTable(_) => {
                        unreachable!("analytic reference fixtures only")
                    }
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
            let medium = HomogeneousMedium::try_new(0.04, 0.08, phase.clone()).unwrap();
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
            let expected = airlight_reference(&phase, 0.04, 0.08);
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

const CROSSING_M: f64 = 1.0 / 16.0;
const PEAK_CANDELA: f64 = 100.0;

fn narrow_settings() -> ScatterSettings {
    ScatterSettings {
        initial_panels: 2,
        max_panels: 4096,
        absolute_tolerance_cd_m2: 0.0,
        relative_tolerance: 1.0e-8,
    }
}

// Source distance to the camera line is exactly 1 m and its axial coordinate
// is c. With t=atan(s-c), ds/r²=dt and s+r=c+tan(t)+sec(t). The supplied
// profile and integration intervals below come from the declared table or
// cone, not from production source evaluation or the partition helper.
fn narrow_reference(intervals: &[(f64, f64)], profile: impl Fn(f64) -> f64) -> f64 {
    const PANELS: u32 = 1024;
    let value = |t: f64| profile(t) * (-0.15 * (CROSSING_M + t.tan() + 1.0 / t.cos())).exp();
    let integral: f64 = intervals
        .iter()
        .map(|&(low, high)| {
            let step = (high - low) / f64::from(PANELS);
            let mut sum = value(low) + value(high);
            for i in 1..PANELS {
                sum += if i % 2 == 0 { 2.0 } else { 4.0 } * value(low + f64::from(i) * step);
            }
            sum * step / 3.0
        })
        .sum();
    PEAK_CANDELA * 0.1 * integral / (4.0 * PI)
}

fn narrow_source_frame() -> RigidTransform {
    let q = core::f64::consts::FRAC_1_SQRT_2;
    RigidTransform::try_new([0.0, q, 0.0, q], [-1.0, 0.0, CROSSING_M]).unwrap()
}

fn check_narrow_result(result: &LuminanceEstimate, expected: f64) {
    // A fixed oracle tolerance, not a confidence interval or a numerical bound.
    close(result.luminance_cd_m2, expected, 2.0e-7);
    assert_eq!(result.surface_cd_m2, 0.0);
    // A cap can prevent individual intervals from meeting their requested
    // discrepancy even when this independent whole-integral check passes.
    assert!(result.single_scatter.panels <= narrow_settings().max_panels);
    assert!(result.single_scatter.evaluations < 4 * u64::from(narrow_settings().max_panels));
    assert!(result.total_error_bound_cd_m2.is_none());
    assert!(
        result
            .single_scatter
            .omitted_multiple_scattering_bound_cd_m2
            .is_none()
    );
}

#[test]
fn g1_scene_partition_resolves_narrow_cone_in_rigid_frames() {
    with_cx(|cx, _| {
        let half_angle = 0.01_f64;
        let expected = narrow_reference(&[(-half_angle, half_angle)], |_| 1.0);
        // All illuminated points lie between the original 2/4-panel nodes.
        // The second case rotates the complete scene +90 degrees about Y,
        // then translates it. Its scalar reference is unchanged, with no
        // production transform or source-boundary lookup used in the oracle.
        for (frame, ray) in [
            (narrow_source_frame(), forward_ray()),
            (
                RigidTransform::try_new([0.0, 1.0, 0.0, 0.0], [2.0 + CROSSING_M, -3.0, 5.0])
                    .unwrap(),
                Ray {
                    origin: Point3::new(2.0, -3.0, 4.0),
                    dir: Vec3::new(1.0, 0.0, 0.0),
                },
            ),
        ] {
            let sources = [PhotometricSource::try_new(
                frame,
                PEAK_CANDELA * 2.0 * PI * (1.0 - half_angle.cos()),
                1.0,
                AngularDistribution::uniform_cone(half_angle).unwrap(),
                SourceSpectrum::equal_energy(),
            )
            .unwrap()];
            let scene = PhotometricScene::try_new(
                &[],
                &sources,
                HomogeneousMedium::try_new(0.05, 0.1, PhaseFunction::Isotropic).unwrap(),
                TransportOptions::default(),
                cx,
            )
            .unwrap();
            let result = scene
                .luminance_along_ray(cx, &ray, 1.0, narrow_settings())
                .unwrap();
            check_narrow_result(&result, expected);
            assert!(result.single_scatter.tolerance_met);
            // Extinction over the longest path supplies an independent
            // positive lower envelope for this single-scattering model.
            let longest = CROSSING_M + half_angle.tan() + 1.0 / half_angle.cos();
            let lower = PEAK_CANDELA * 0.1 * half_angle / (2.0 * PI) * (-0.15 * longest).exp();
            assert!(result.luminance_cd_m2 >= lower);
        }
    });
}

#[test]
fn g1_scene_partition_resolves_narrow_tabulated_polar_peak() {
    with_cx(|cx, _| {
        let (low, peak, high) = (0.01_f64, 0.012_f64, 0.014_f64);
        let (mu_low, mu_peak, mu_high) = (low.cos(), peak.cos(), high.cos());
        // An off-axis annular peak: zero on axis and outside [low,high].
        // Integrating its two linear-in-cos triangles gives pi*(mu_low-mu_high)
        // steradians, independently of the location of the middle knot.
        let distribution = AngularDistribution::tabulated(
            vec![0.0, low, peak, high, PI],
            vec![0.0, 2.0 * PI],
            vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
        )
        .unwrap();
        let sources = [PhotometricSource::try_new(
            narrow_source_frame(),
            PEAK_CANDELA * PI * (mu_low - mu_high),
            1.0,
            distribution,
            SourceSpectrum::equal_energy(),
        )
        .unwrap()];
        let scene = PhotometricScene::try_new(
            &[],
            &sources,
            HomogeneousMedium::try_new(0.05, 0.1, PhaseFunction::Isotropic).unwrap(),
            TransportOptions::default(),
            cx,
        )
        .unwrap();
        let expected = narrow_reference(
            &[(-high, -peak), (-peak, -low), (low, peak), (peak, high)],
            |t| {
                if t.abs() <= peak {
                    (mu_low - t.cos()) / (mu_low - mu_peak)
                } else {
                    (t.cos() - mu_high) / (mu_peak - mu_high)
                }
            },
        );
        let result = scene
            .luminance_along_ray(cx, &forward_ray(), 1.0, narrow_settings())
            .unwrap();
        check_narrow_result(&result, expected);
    });
}

#[test]
fn g1_scene_partition_resolves_tabulated_azimuth_seam() {
    with_cx(|cx, _| {
        let half_width = 0.01_f64;
        // The fixture's optical axis is world +Y. On the entire camera ray,
        // local theta=pi/2 and phi=-atan(s-c), traversing the 0/2pi seam.
        // The table is (1-|cos(theta)|) times a triangular azimuth profile;
        // its exact represented solid-angle integral is half_width sr.
        let q = core::f64::consts::FRAC_1_SQRT_2;
        let distribution = AngularDistribution::tabulated(
            vec![0.0, PI / 2.0, PI],
            vec![0.0, half_width, 2.0 * PI - half_width, 2.0 * PI],
            vec![0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
        )
        .unwrap();
        let sources = [PhotometricSource::try_new(
            RigidTransform::try_new([-q, 0.0, 0.0, q], [-1.0, 0.0, CROSSING_M]).unwrap(),
            PEAK_CANDELA * half_width,
            1.0,
            distribution,
            SourceSpectrum::equal_energy(),
        )
        .unwrap()];
        let scene = PhotometricScene::try_new(
            &[],
            &sources,
            HomogeneousMedium::try_new(0.05, 0.1, PhaseFunction::Isotropic).unwrap(),
            TransportOptions::default(),
            cx,
        )
        .unwrap();
        let expected = narrow_reference(&[(-half_width, 0.0), (0.0, half_width)], |t| {
            1.0 - t.abs() / half_width
        });
        let result = scene
            .luminance_along_ray(cx, &forward_ray(), 1.0, narrow_settings())
            .unwrap();
        check_narrow_result(&result, expected);
    });
}

#[test]
fn g4_scattering_partition_budget_is_global_and_repeated_knots_are_shared() {
    with_cx(|cx, _| {
        let half_angle = 0.01_f64;
        let emitter = PhotometricSource::try_new(
            narrow_source_frame(),
            PEAK_CANDELA * 2.0 * PI * (1.0 - half_angle.cos()),
            1.0,
            AngularDistribution::uniform_cone(half_angle).unwrap(),
            SourceSpectrum::equal_energy(),
        )
        .unwrap();
        let sources = vec![emitter; 4];
        let scene = PhotometricScene::try_new(
            &[],
            &sources,
            HomogeneousMedium::try_new(0.05, 0.1, PhaseFunction::Isotropic).unwrap(),
            TransportOptions::default(),
            cx,
        )
        .unwrap();
        let settings = ScatterSettings {
            max_panels: 2,
            ..narrow_settings()
        };
        assert!(matches!(
            scene.luminance_along_ray(cx, &forward_ray(), 1.0, settings),
            Err(PhotometryError::ScatteringPartitionBudgetExceeded { max_panels: 2 })
        ));
        // Four coincident fixtures share the same two boundary locations and
        // therefore only three mandatory intervals, not four copies of them.
        let result = scene
            .luminance_along_ray(
                cx,
                &forward_ray(),
                1.0,
                ScatterSettings {
                    initial_panels: 1,
                    max_panels: 3,
                    ..settings
                },
            )
            .unwrap();
        assert!(result.single_scatter.panels <= 3);
        let expected = 4.0 * narrow_reference(&[(-half_angle, half_angle)], |_| 1.0);
        close(result.luminance_cd_m2, expected, 1.0e-3);
        assert!(
            result
                .single_scatter
                .estimated_quadrature_error_cd_m2
                .is_none()
        );
        assert!(!result.single_scatter.tolerance_met);
        assert!(result.total_error_bound_cd_m2.is_none());
    });
}

#[test]
fn g1_raw_callback_discrepancy_does_not_bound_an_unsampled_beam() {
    with_cx(|cx, _| {
        let half_angle = 0.01_f64;
        let medium = HomogeneousMedium::try_new(0.05, 0.1, PhaseFunction::Isotropic).unwrap();
        // The generic callback API has no declared source profile or support
        // boundaries. Its refinement discrepancy retains its heuristic meaning.
        let aliased = fs_render::photometry::medium::single_scattered_luminance(
            &medium,
            &forward_ray(),
            1.0,
            narrow_settings(),
            cx,
            |point, _| {
                let axial = point.z - CROSSING_M;
                let radius_squared = 1.0 + axial * axial;
                Ok::<_, ()>(if axial.abs() <= half_angle.tan() {
                    PEAK_CANDELA * (-0.15 * radius_squared.sqrt()).exp()
                        / (4.0 * PI * radius_squared)
                } else {
                    0.0
                })
            },
        )
        .unwrap();
        assert_eq!(aliased.luminance_cd_m2, 0.0);
        assert_eq!(aliased.panels, 4);
        assert_eq!(aliased.estimated_quadrature_error_cd_m2, Some(0.0));
        assert!(aliased.tolerance_met);
        assert!(aliased.omitted_multiple_scattering_bound_cd_m2.is_none());
        assert!(narrow_reference(&[(-half_angle, half_angle)], |_| 1.0) > 0.01);
    });
}
