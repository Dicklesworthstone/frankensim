use super::*;
use crate::charts::{Ray, TriMesh};
use crate::dielectric::{
    BeerLambertAbsorption, CauchyIor, DielectricGlass, DielectricIor, DielectricIorSample,
    DielectricSurface, GlassProvenance, TabulatedIor,
};
use crate::gradient_index::{GradientIndexOptions, event::GradientEventOptions};
use crate::instances::{GeometryInstance, InstanceScene, RigidTransform, SharedGeometry};
use crate::photometry::source::bundle::PlanarIrradianceReceiver;
use crate::photometry::source::{
    AngularDistribution, BIN_WIDTH_NM, LUMINOUS_EFFICACY_LM_W, SourceSpectrum, photopic_weight_nm,
};
use crate::polarization::dielectric::scene::curved::bundle::CurvedBundleSettings;
use crate::polarization::dielectric::scene::curved::bundle::connection::{
    BundleBranchPolicy, BundleConnectionOptions,
    images::{BundleImageOptions, BundleImageSeed, BundleImageSeedOutcome},
};
use crate::polarization::dielectric::scene::curved::{CurvedTraceSettings, SceneMedium};
use crate::polarization::dielectric::scene::{Interaction, SurfaceBinding};
use crate::polarization::photometry::SpectralOptics;
use fs_blake3::hash_domain;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3, fixtures::SphereChart};

const TARGET: Point3 = Point3::new(0.03, 0.0, 2.5);

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 101,
                kernel_id: 9,
                tile: 0,
                iteration: 0,
            },
            asupersync::types::Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate);
    });
}

fn source() -> PhotometricSource {
    PhotometricSource::try_new(
        RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [0.0, 0.0, -3.0]).unwrap(),
        1200.0,
        0.75,
        AngularDistribution::cosine_power(2.0, core::f64::consts::FRAC_PI_2).unwrap(),
        SourceSpectrum::equal_energy(),
    )
    .unwrap()
}

fn glass(ior: impl Into<DielectricIor>, absorption: f64) -> DielectricGlass {
    DielectricGlass::new(
        ior,
        BeerLambertAbsorption::try_constant(absorption).unwrap(),
        GlassProvenance::Custom,
    )
}

fn geometry() -> InstanceScene {
    let instance = |id: u64, shape| {
        GeometryInstance::try_new(
            id,
            hash_domain("multiple-image-photometry-test", &id.to_le_bytes()),
            shape,
            RigidTransform::identity(),
        )
        .unwrap()
    };
    InstanceScene::try_new(vec![
        instance(
            1,
            SharedGeometry::chart(SphereChart {
                center: Point3::new(0.0, 0.0, 0.0),
                radius: 1.0,
            }),
        ),
        instance(
            9,
            SharedGeometry::mesh(TriMesh::new(
                vec![[-0.6, -0.6, 2.5], [0.6, -0.6, 2.5], [0.0, 0.6, 2.5]],
                vec![[0, 1, 2]],
            )),
        ),
    ])
    .unwrap()
}

fn bindings() -> Vec<SurfaceBinding> {
    vec![
        SurfaceBinding {
            object_id: 1,
            interaction: Interaction::Dielectric {
                outside: 1,
                inside: 0,
                surface: DielectricSurface::SMOOTH,
            },
        },
        SurfaceBinding {
            object_id: 9,
            interaction: Interaction::Detector {
                reference_axis: Vec3::new(0.0, 1.0, 0.0),
            },
        },
    ]
}

fn request() -> ImageIrradianceRequest {
    ImageIrradianceRequest {
        detector_object_id: 9,
        receiver: PlanarIrradianceReceiver::try_new(TARGET, Vec3::new(0.0, 0.0, -1.0), 2e-8)
            .unwrap(),
        source_transverse_axis: Vec3::new(1.0, 0.0, 0.0),
        initial_medium: 1,
        discovery: BundleImageOptions::default(),
        max_continuation_step_rad: 0.02,
        connection: BundleConnectionOptions {
            scene: CurvedBundleSettings {
                trace: CurvedTraceSettings {
                    max_surface_events: 8,
                    max_distance_m: 8.0,
                    max_total_steps: 16_384,
                    max_query_work: 1_000_000,
                    event: GradientEventOptions {
                        integration: GradientIndexOptions {
                            position_tolerance_m: 1e-12,
                            direction_tolerance: 1e-13,
                            optical_path_tolerance_m: 1e-12,
                            optical_depth_tolerance: 1e-13,
                            relative_tolerance: 1e-12,
                            ..GradientIndexOptions::default()
                        },
                        surface_tolerance_m: 1e-9,
                        ..GradientEventOptions::default()
                    },
                },
                ..CurvedBundleSettings::default()
            },
            ..BundleConnectionOptions::default()
        },
    }
}

fn direction(theta: f64) -> Vec3 {
    Vec3::new(theta.sin(), 0.0, theta.cos())
}

fn seeds() -> [BundleImageSeed<'static>; 3] {
    [-0.20, 0.06, 0.14].map(|theta| BundleImageSeed {
        initial_direction: direction(theta),
        policy: BundleBranchPolicy::TransmitOrReflectTir,
    })
}

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-15),
        "{actual:.16e} != {expected:.16e}"
    );
}

fn unit(v: Vec3) -> Vec3 {
    v.scale(1.0 / v.dot(v).sqrt())
}
fn add(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(a.x + b.x, a.y + b.y, a.z + b.z)
}
fn cross(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    )
}

fn snell(direction: Vec3, outward: Vec3, ratio: f64) -> Vec3 {
    let normal = if direction.dot(outward) < 0.0 {
        outward
    } else {
        outward.scale(-1.0)
    };
    let ci = -direction.dot(normal);
    let ct = (1.0 - ratio * ratio * (1.0 - ci * ci)).sqrt();
    unit(add(direction.scale(ratio), normal.scale(ratio * ci - ct)))
}

// Independent quadratic sphere intersections, Snell directions and Fresnel
// channel powers. Production tracing and bundle derivatives are not used.
fn reference(ray: Ray, index: f64) -> (Point3, f64, [f64; 2]) {
    let origin = Vec3::new(ray.origin.x, ray.origin.y, ray.origin.z);
    let b = origin.dot(ray.dir);
    let entry = ray.at(-b - (b * b - origin.dot(origin) + 1.0).sqrt());
    let normal = unit(Vec3::new(entry.x, entry.y, entry.z));
    let inside = snell(ray.dir, normal, 1.0 / index);
    let chord = -2.0 * inside.dot(Vec3::new(entry.x, entry.y, entry.z));
    let exit = entry.offset(inside.scale(chord));
    let outgoing = snell(inside, unit(Vec3::new(exit.x, exit.y, exit.z)), index);
    let ci = -ray.dir.dot(normal);
    let ct = -inside.dot(normal);
    let ts = 1.0 - ((ci - index * ct) / (ci + index * ct)).powi(2);
    let tp = 1.0 - ((index * ci - ct) / (index * ci + ct)).powi(2);
    (
        exit.offset(outgoing.scale((TARGET.z - exit.z) / outgoing.z)),
        chord,
        [tp * tp, ts * ts],
    )
}

fn reference_ray(index: f64, image: usize) -> Ray {
    let (mut lo, mut hi) = [(-0.22, -0.16), (0.04, 0.10), (0.11, 0.16)][image];
    let origin = Point3::new(0.0, 0.0, -3.0);
    let residual = |theta| {
        reference(
            Ray {
                origin,
                dir: direction(theta),
            },
            index,
        )
        .0
        .x - TARGET.x
    };
    assert!(residual(lo) * residual(hi) < 0.0);
    for _ in 0..64 {
        let mid = 0.5 * (lo + hi);
        if residual(lo) * residual(mid) <= 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    Ray {
        origin,
        dir: direction(0.5 * (lo + hi)),
    }
}

fn reference_area(ray: Ray, index: f64) -> f64 {
    let first = unit(add(Vec3::new(1.0, 0.0, 0.0), ray.dir.scale(-ray.dir.x)));
    let axes = [first, cross(ray.dir, first)];
    let h = 1e-6;
    let columns = axes.map(|axis| {
        let endpoint = |angle| {
            reference(
                Ray {
                    dir: unit(add(ray.dir, axis.scale(angle))),
                    ..ray
                },
                index,
            )
            .0
        };
        endpoint(h).delta_from(endpoint(-h)).scale(0.5 / h)
    });
    cross(columns[0], columns[1]).z.abs()
}

#[test]
fn g1_three_same_history_images_match_independent_spectral_snell_fresnel_and_area() {
    with_cx(|cx, _| {
        let source = source();
        let geometry = geometry();
        let air = glass(CauchyIor::try_constant(1.0).unwrap(), 0.0);
        let lens = glass(CauchyIor::try_new(1.498, 0.0005, 0.0).unwrap(), 0.15);
        let media = [
            SceneMedium::Homogeneous(&lens),
            SceneMedium::Homogeneous(&air),
        ];
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let result = source
            .discover_irradiance_images(cx, &scene, request(), &seeds())
            .unwrap();
        assert_eq!(result.images().len(), 3);
        let signature = result.images()[0].branch_signature();
        let mut expected_readings = [[0.0; 2]; 3];
        for image in result.images() {
            assert_eq!(image.branch_signature(), signature);
            assert!(
                (image.launch_directions()[0].x - image.launch_directions()[79].x).abs() > 1e-3
            );
        }
        for bin in 0..SPECTRAL_BINS {
            let wavelength = wavelength_nm(bin);
            let index = 1.498 + 0.0005 / (wavelength * 1e-3).powi(2);
            let mut total = 0.0;
            for (image, expected_reading) in expected_readings.iter_mut().enumerate() {
                let ray = reference_ray(index, image);
                let (endpoint, chord, channels) = reference(ray, index);
                close(endpoint.x, TARGET.x, 1e-10);
                let connected = &result.images()[image];
                let native_direction = connected.launch_directions()[bin];
                assert!(add(native_direction, ray.dir.scale(-1.0)).norm() < 2e-6);
                let area = reference_area(ray, index);
                assert!(area > 0.03, "the three images must be regular");
                let expected = source
                    .spectral_intensity_w_sr_nm(ray.dir, wavelength)
                    .unwrap()
                    * 0.5
                    * (channels[0] + channels[1])
                    * (-0.15 * chord).exp()
                    / area;
                close(connected.spectral_irradiance_w_m2_nm()[bin], expected, 2e-4);
                expected_reading[0] += expected * BIN_WIDTH_NM;
                expected_reading[1] += expected * photopic_weight_nm(bin) * LUMINOUS_EFFICACY_LM_W;
                total += expected;
            }
            close(result.spectral_irradiance_w_m2_nm()[bin], total, 2e-4);
        }
        for (image, expected) in result.images().iter().zip(expected_readings) {
            close(image.reading().radiant_irradiance_w_m2, expected[0], 2e-4);
            close(image.reading().illuminance_lux, expected[1], 2e-4);
        }
        close(
            result.reading().radiant_irradiance_w_m2,
            expected_readings.iter().map(|value| value[0]).sum(),
            2e-4,
        );
        close(
            result.reading().illuminance_lux,
            expected_readings.iter().map(|value| value[1]).sum(),
            2e-4,
        );
        assert!(result.trace_evaluations() > (3 * SPECTRAL_BINS) as u64);

        // Repeated seeds and a backward ray must not inject additional light.
        let family = seeds();
        let repeated = [
            family[0],
            family[1],
            family[2],
            family[1],
            BundleImageSeed {
                initial_direction: Vec3::new(0.0, 0.0, -1.0),
                policy: BundleBranchPolicy::TransmitOrReflectTir,
            },
        ];
        let replay = source
            .discover_irradiance_images(cx, &scene, request(), &repeated)
            .unwrap();
        assert_eq!(replay.images().len(), 3);
        assert_eq!(
            replay.spectral_irradiance_w_m2_nm(),
            result.spectral_irradiance_w_m2_nm()
        );
        assert!(matches!(
            replay.discovery().seed_results()[3].outcome(),
            BundleImageSeedOutcome::Duplicate { image: 1 }
        ));
        assert!(matches!(
            replay.discovery().seed_results()[4].outcome(),
            BundleImageSeedOutcome::Miss(_)
        ));
        assert!(replay.trace_evaluations() > result.trace_evaluations());
    });
}

#[test]
fn g1_polarized_multiple_images_preserve_each_channel_and_upstream_filter_once() {
    with_cx(|cx, _| {
        let source = source();
        let geometry = geometry();
        let air = glass(CauchyIor::try_constant(1.0).unwrap(), 0.0);
        let lens = glass(CauchyIor::try_constant(1.5).unwrap(), 0.15);
        let media = [
            SceneMedium::Homogeneous(&lens),
            SceneMedium::Homogeneous(&air),
        ];
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let parallel = core::array::from_fn(|bin| 0.2 + 0.6 * bin as f64 / 79.0);
        let perpendicular = [0.1; SPECTRAL_BINS];
        let neutral = core::array::from_fn(|bin| 0.6 - 0.2 * bin as f64 / 79.0);
        let filtered = source
            .clone()
            .with_unpolarized_optics()
            .through(&SpectralOptics::linear_diattenuator(0.0, parallel, perpendicular).unwrap())
            .unwrap()
            .through(&SpectralOptics::neutral(neutral).unwrap())
            .unwrap();
        let result = filtered
            .discover_irradiance_images(cx, &scene, request(), &seeds())
            .unwrap();
        assert_eq!(result.images().len(), 3);
        let references = [0, 1, 2].map(|image| {
            let ray = reference_ray(1.5, image);
            let (_, chord, channels) = reference(ray, 1.5);
            (
                ray.dir,
                channels,
                (-0.15 * chord).exp() / reference_area(ray, 1.5),
            )
        });
        let mut expected_lux = 0.0;
        for bin in 0..SPECTRAL_BINS {
            let mut total = 0.0;
            for (image, (direction, channels, attenuation)) in references.iter().enumerate() {
                // Projected world-x is p polarization for every xz-plane ray;
                // all source filtering and both native interfaces apply once.
                let expected = source
                    .spectral_intensity_w_sr_nm(*direction, wavelength_nm(bin))
                    .unwrap()
                    * 0.5
                    * neutral[bin]
                    * (parallel[bin] * channels[0] + perpendicular[bin] * channels[1])
                    * attenuation;
                close(
                    result.images()[image].spectral_irradiance_w_m2_nm()[bin],
                    expected,
                    2e-4,
                );
                total += expected;
            }
            close(result.spectral_irradiance_w_m2_nm()[bin], total, 2e-4);
            expected_lux += total * photopic_weight_nm(bin) * LUMINOUS_EFFICACY_LM_W;
        }
        close(result.reading().illuminance_lux, expected_lux, 2e-4);
        assert!(filtered.intensity_source().flux_lm() < 0.3 * source.flux_lm());
    });
}

#[test]
fn g4_spectral_fold_disappearance_refuses_partial_or_duplicate_image_power() {
    with_cx(|cx, _| {
        let source = source();
        let geometry = geometry();
        let air = glass(CauchyIor::try_constant(1.0).unwrap(), 0.0);
        // This declared increasing phase law crosses the ball lens's fold.
        // Three roots exist in the first bin; the two positive-angle images
        // coalesce and disappear while the negative-angle image survives.
        let law = TabulatedIor::try_new(core::array::from_fn(|k| {
            DielectricIorSample::try_new(380.0 + 50.0 * k as f64, 1.5 + 0.00125 * k as f64).unwrap()
        }))
        .unwrap();
        let lens = glass(law, 0.15);
        let media = [
            SceneMedium::Homogeneous(&lens),
            SceneMedium::Homogeneous(&air),
        ];
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let error = source
            .discover_irradiance_images(cx, &scene, request(), &seeds())
            .unwrap_err();
        let bin = match error {
            ImageIrradianceError::Connection { bin, .. }
            | ImageIrradianceError::ImageCollision { bin, .. }
            | ImageIrradianceError::AmbiguousImages { bin, .. }
            | ImageIrradianceError::SpectralParityChanged { bin, .. }
            | ImageIrradianceError::ContinuationJump { bin, .. } => bin,
            other => panic!("unexpected spectral fold refusal: {other:?}"),
        };
        assert!(
            bin > 0 && bin < SPECTRAL_BINS,
            "first-bin discovery must succeed"
        );
    });
}
