use super::*;
use crate::charts::TriMesh;
use crate::dielectric::{
    BeerLambertAbsorption, CauchyIor, DielectricGlass, DielectricSurface, GlassProvenance,
};
use crate::gradient_index::{
    GradientIndexError, GradientIndexField, GradientIndexOptions, GradientIndexSample,
    event::GradientEventOptions,
};
use crate::instances::{GeometryInstance, InstanceScene, RigidTransform, SharedGeometry};
use crate::photometry::source::{
    AngularDistribution, BIN_WIDTH_NM, LUMINOUS_EFFICACY_LM_W, SourceSpectrum, photopic_weight_nm,
};
use crate::polarization::dielectric::scene::curved::SceneMedium;
use crate::polarization::dielectric::scene::{Interaction, SurfaceBinding};
use crate::polarization::photometry::SpectralOptics;
use fs_blake3::hash_domain;
use fs_exec::{CancelGate, ExecMode, StreamKey};
use fs_geom::Point3;

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 93,
                kernel_id: 8,
                tile: 0,
                iteration: 0,
            },
            asupersync::types::Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate);
    });
}

fn source(distribution: AngularDistribution) -> PhotometricSource {
    PhotometricSource::try_new(
        RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [0.0, 0.0, -1.0]).unwrap(),
        1200.0,
        0.75,
        distribution,
        SourceSpectrum::equal_energy(),
    )
    .unwrap()
}

fn glass(index: f64, absorption: f64) -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(index).unwrap(),
        BeerLambertAbsorption::try_constant(absorption).unwrap(),
        GlassProvenance::Custom,
    )
}

fn square(id: u64, z: f64, half_width: f64, positive: bool) -> GeometryInstance {
    let triangles = if positive {
        vec![[0, 1, 2], [0, 2, 3]]
    } else {
        vec![[0, 2, 1], [0, 3, 2]]
    };
    GeometryInstance::try_new(
        id,
        hash_domain("finite-detector-test", &id.to_le_bytes()),
        SharedGeometry::mesh(TriMesh::new(
            vec![
                [-half_width, -half_width, z],
                [half_width, -half_width, z],
                [half_width, half_width, z],
                [-half_width, half_width, z],
            ],
            triangles,
        )),
        RigidTransform::identity(),
    )
    .unwrap()
}

fn detector() -> SurfaceBinding {
    SurfaceBinding {
        object_id: 9,
        interaction: Interaction::Detector {
            reference_axis: Vec3::new(0.0, 1.0, 0.0),
        },
    }
}

fn request(samples: u32) -> FiniteDetectorRequest {
    FiniteDetectorRequest {
        detector_object_id: 9,
        source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
        initial_medium: 0,
        samples,
        seed: 0x7318_2a93_1c74_b6e1,
        trace: CurvedTraceSettings {
            max_surface_events: 16,
            max_distance_m: 12.0,
            max_total_steps: 16_384,
            max_query_work: 1_000_000,
            event: GradientEventOptions {
                integration: GradientIndexOptions {
                    position_tolerance_m: 1e-10,
                    direction_tolerance: 1e-11,
                    max_step_m: 0.5,
                    ..GradientIndexOptions::default()
                },
                surface_tolerance_m: 1e-7,
                ..GradientEventOptions::default()
            },
        },
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

#[test]
fn g1_square_solid_angle_counts_misses_and_retains_spectral_covariance() {
    with_cx(|cx, _| {
        let source = source(AngularDistribution::isotropic());
        let geometry = InstanceScene::try_new(vec![square(9, 0.0, 1.0, true)]).unwrap();
        let air = glass(1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&air)];
        let scene = CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![detector()], 100_000)
            .unwrap();
        let options = request(128);
        let result = source.measure_detector(cx, &scene, options).unwrap();
        assert_eq!(result.sample_count(), options.samples);
        assert_eq!(
            result.launched_paths(),
            u64::from(options.samples) * SPECTRAL_BINS as u64
        );
        // Nondispersive vacuum and a common source direction give the same hit
        // event at every wavelength. The detector subtends Ω=2π/3, or 1/6 sphere.
        assert_eq!(result.detected_paths() % SPECTRAL_BINS as u64, 0);
        let hits = result.detected_paths() / SPECTRAL_BINS as u64;
        assert!(hits > 0 && hits < u64::from(options.samples));
        assert_eq!(u64::from(result.received_samples()), hits);
        let fraction = hits as f64 / f64::from(options.samples);
        let bernoulli_sem = (fraction * (1.0 - fraction) / f64::from(options.samples - 1)).sqrt();
        let solid_angle = 4.0 * (1.0 / 3.0_f64.sqrt()).atan();
        let exact_fraction = solid_angle / (4.0 * core::f64::consts::PI);
        close(fraction, exact_fraction, 5.0 * bernoulli_sem);

        let mut emitted_w = 0.0;
        let mut emitted_lm = 0.0;
        for bin in 0..SPECTRAL_BINS {
            let power = source.flux_lm() * source.spectrum().bins_w_per_lumen_nm()[bin];
            close(
                result.spectral_power_w_nm()[bin].i,
                power * fraction,
                2e-13 * power,
            );
            close(
                result.spectral_standard_error_w_nm()[bin][0],
                power * bernoulli_sem,
                2e-13 * power,
            );
            emitted_w += power * BIN_WIDTH_NM;
            emitted_lm += power * LUMINOUS_EFFICACY_LM_W * photopic_weight_nm(bin);
        }
        close(
            result.radiant_stokes_w().i,
            emitted_w * fraction,
            2e-12 * emitted_w,
        );
        close(
            result.photopic_stokes_lm().i,
            emitted_lm * fraction,
            2e-12 * emitted_lm,
        );
        // Summing variances as if the 80 bins were independent would fail by a
        // large factor. Integrate each whole sampled spectrum BEFORE its variance.
        close(
            result.radiant_standard_error_w()[0],
            emitted_w * bernoulli_sem,
            2e-12 * emitted_w,
        );
        close(
            result.photopic_standard_error_lm()[0],
            emitted_lm * bernoulli_sem,
            2e-12 * emitted_lm,
        );
        // This detector's actual square area is 4 m²; averaging power over it
        // must scale uncertainty and illuminance units along with the mean.
        let average = result.area_average(4.0).unwrap();
        close(
            average.radiant_w_m2.i,
            emitted_w * fraction / 4.0,
            2e-12 * emitted_w,
        );
        close(
            average.photopic_lux.i,
            emitted_lm * fraction / 4.0,
            2e-12 * emitted_lm,
        );
        close(
            average.photopic_standard_error_lux[0],
            emitted_lm * bernoulli_sem / 4.0,
            2e-12 * emitted_lm,
        );
        assert!(result.area_average(0.0).is_err());
        assert!(result.max_contact_resolution_m() > 0.0);
    });
}

#[test]
fn g1_full_cone_capture_preserves_calibration_and_filtered_polarized_power_once() {
    with_cx(|cx, _| {
        let original = source(AngularDistribution::uniform_cone(0.2).unwrap());
        let transmission = std::array::from_fn(|bin| if bin % 2 == 0 { 0.25 } else { 0.75 });
        let polarized = original
            .clone()
            .with_unpolarized_optics()
            .through(&SpectralOptics::neutral(transmission).unwrap())
            .unwrap()
            .through(&SpectralOptics::linear_polarizer(0.0).unwrap())
            .unwrap();
        let geometry = InstanceScene::try_new(vec![square(9, 0.0, 2.0, true)]).unwrap();
        let air = glass(1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&air)];
        let scene = CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![detector()], 100_000)
            .unwrap();
        let result = polarized.measure_detector(cx, &scene, request(4)).unwrap();
        assert_eq!(result.detected_paths(), 4 * SPECTRAL_BINS as u64);
        for bin in 0..SPECTRAL_BINS {
            let expected = 0.5
                * transmission[bin]
                * original.flux_lm()
                * original.spectrum().bins_w_per_lumen_nm()[bin];
            let state = result.spectral_power_w_nm()[bin];
            close(state.i, expected, 2e-12 * expected);
            close(state.q, expected, 2e-12 * expected);
            close(state.u, 0.0, 2e-12 * expected);
            close(state.v, 0.0, 2e-12 * expected);
        }
        close(
            result.photopic_stokes_lm().i,
            polarized.intensity_source().flux_lm(),
            2e-9,
        );
        assert!(result.photopic_standard_error_lm()[0] < 1e-9);
    });
}

#[test]
fn g1_absorbing_refractive_slab_matches_incoherent_fresnel_series() {
    with_cx(|cx, _| {
        let source = source(AngularDistribution::uniform_cone(0.001).unwrap());
        let geometry = InstanceScene::try_new(vec![
            square(1, 0.0, 2.0, false),
            square(2, 1.0, 2.0, true),
            square(9, 2.0, 2.0, true),
        ])
        .unwrap();
        let air = glass(1.0, 0.0);
        let material = glass(1.5, 0.4);
        let media = [
            SceneMedium::Homogeneous(&air),
            SceneMedium::Homogeneous(&material),
        ];
        let bindings = vec![
            SurfaceBinding {
                object_id: 1,
                interaction: Interaction::Dielectric {
                    outside: 0,
                    inside: 1,
                    surface: DielectricSurface::SMOOTH,
                },
            },
            SurfaceBinding {
                object_id: 2,
                interaction: Interaction::Dielectric {
                    outside: 0,
                    inside: 1,
                    surface: DielectricSurface::SMOOTH,
                },
            },
            detector(),
        ];
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings, 100_000).unwrap();
        let result = source.measure_detector(cx, &scene, request(128)).unwrap();
        // At normal incidence R=((n-1)/(n+1))². Sum every even internal
        // reflection before transmission, with Beer loss on the actual length.
        // The 1 mrad cone changes this limit by less than 2e-6 relative here.
        let reflectance = ((1.5_f64 - 1.0) / (1.5 + 1.0)).powi(2);
        let beer = (-0.4_f64).exp();
        let expected_fraction =
            (1.0 - reflectance).powi(2) * beer / (1.0 - reflectance.powi(2) * beer.powi(2));
        close(
            result.photopic_stokes_lm().i,
            source.flux_lm() * expected_fraction,
            5.0 * result.photopic_standard_error_lm()[0] + 2e-6 * source.flux_lm(),
        );
        assert!(result.detected_paths() < result.launched_paths());
        assert!(result.photopic_stokes_lm().i < source.flux_lm());
    });
}

#[test]
fn g1_finite_detector_at_true_grin_focus_retains_finite_total_power() {
    // n(r)=n0 sech(g*r), with g=1/m. Translational symmetry conserves
    // p_z=n0*cos(theta); in each meridional plane the exact ray satisfies
    // sinh(g*r)=tan(theta)*sin(g*(z-z_source)). Thus all forward rays from
    // the axis meet again at z-z_source=pi/g, even away from the paraxial
    // limit. Point irradiance is singular; finite captured flux is not.
    struct FocusingField;
    impl GradientIndexField for FocusingField {
        fn sample(
            &self,
            position: Point3,
            _: f64,
        ) -> Result<GradientIndexSample, GradientIndexError> {
            let radius = position.x.hypot(position.y);
            let index = 1.5 / radius.cosh();
            let radial_derivative = if radius == 0.0 {
                0.0
            } else {
                -index * radius.tanh() / radius
            };
            GradientIndexSample::try_new(
                index,
                Vec3::new(
                    radial_derivative * position.x,
                    radial_derivative * position.y,
                    0.0,
                ),
                0.0,
                Some(index),
            )
        }
    }
    with_cx(|cx, _| {
        let source = source(AngularDistribution::uniform_cone(0.2).unwrap());
        let geometry =
            InstanceScene::try_new(vec![square(9, core::f64::consts::PI - 1.0, 0.01, true)])
                .unwrap();
        let field = FocusingField;
        // This field intentionally supplies no Hessian or bundle derivatives.
        let media = [SceneMedium::GradientIndex {
            field: &field,
            identity: hash_domain("finite-detector-focus", b"n=1.5sech(r),g=1/m"),
        }];
        let scene = CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![detector()], 100_000)
            .unwrap();
        let mut options = request(4);
        options.trace.event.integration.max_step_m = 0.05;
        let result = source.measure_detector(cx, &scene, options).unwrap();
        assert_eq!(result.detected_paths(), 4 * SPECTRAL_BINS as u64);
        close(result.photopic_stokes_lm().i, source.flux_lm(), 2e-9);
        assert!(result.radiant_stokes_w().i.is_finite());
        assert!(result.photopic_standard_error_lm()[0] < 1e-8);
        // The result remains a finite aperture average at the geometric focus.
        let average = result.area_average(4.0 * 0.01_f64.powi(2)).unwrap();
        close(average.photopic_lux.i, source.flux_lm() / 0.0004, 1e-5);
    });
}

#[test]
fn g4_distance_truncation_never_becomes_a_zero_power_measurement() {
    with_cx(|cx, _| {
        let source = source(AngularDistribution::uniform_cone(0.1).unwrap());
        let geometry = InstanceScene::try_new(vec![square(9, 0.0, 2.0, true)]).unwrap();
        let air = glass(1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&air)];
        let scene = CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![detector()], 100_000)
            .unwrap();
        let mut options = request(2);
        options.trace.max_distance_m = 0.1;
        assert!(matches!(
            source.measure_detector(cx, &scene, options),
            Err(FiniteDetectorError::IncompletePath {
                sample: 0,
                bin: 0,
                termination: PathTermination::DistanceLimit
            })
        ));
        options.samples = 0;
        assert!(matches!(
            source.measure_detector(cx, &scene, options),
            Err(FiniteDetectorError::InvalidRequest)
        ));
    });
}
