use super::*;
use crate::charts::TriMesh;
use crate::dielectric::{
    BeerLambertAbsorption, CauchyIor, DielectricGlass, DielectricSurface, GlassProvenance,
};
use crate::gradient_index::{GradientIndexOptions, event::GradientEventOptions};
use crate::instances::{GeometryInstance, InstanceScene, RigidTransform, SharedGeometry};
use crate::photometry::source::{
    AngularDistribution, BIN_WIDTH_NM, LUMINOUS_EFFICACY_LM_W, SourceSpectrum, photopic_weight_nm,
};
use crate::polarization::dielectric::scene::curved::bundle::CurvedBundleSettings;
use crate::polarization::dielectric::scene::curved::{CurvedTraceSettings, SceneMedium};
use crate::polarization::dielectric::scene::{Interaction, SurfaceBinding};
use crate::polarization::photometry::SpectralOptics;
use fs_blake3::hash_domain;
use fs_exec::{CancelGate, ExecMode, StreamKey};
use fs_geom::{Point3, fixtures::SphereChart};

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 81,
                kernel_id: 7,
                tile: 0,
                iteration: 0,
            },
            asupersync::types::Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx);
    });
}

fn fixture(position: [f64; 3]) -> PhotometricSource {
    PhotometricSource::try_new(
        RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], position).unwrap(),
        1200.0,
        0.75,
        AngularDistribution::cosine_power(2.0, core::f64::consts::FRAC_PI_2).unwrap(),
        SourceSpectrum::equal_energy(),
    )
    .unwrap()
}
fn glass(ior: CauchyIor, absorption: f64) -> DielectricGlass {
    DielectricGlass::new(
        ior,
        BeerLambertAbsorption::try_constant(absorption).unwrap(),
        GlassProvenance::Custom,
    )
}
fn sphere() -> GeometryInstance {
    GeometryInstance::try_new(
        1,
        hash_domain("connected-photometry-test", b"unit-ball"),
        SharedGeometry::chart(SphereChart {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        }),
        RigidTransform::identity(),
    )
    .unwrap()
}
fn plane(id: u64, point: Point3, radius: f64) -> GeometryInstance {
    GeometryInstance::try_new(
        id,
        hash_domain("connected-photometry-test", &id.to_le_bytes()),
        SharedGeometry::mesh(TriMesh::new(
            vec![
                [point.x - radius, point.y - radius, point.z],
                [point.x + radius, point.y - radius, point.z],
                [point.x, point.y + radius, point.z],
            ],
            vec![[0, 1, 2]],
        )),
        RigidTransform::identity(),
    )
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
fn request(point: Point3, initial_direction: Vec3) -> ConnectedIrradianceRequest {
    ConnectedIrradianceRequest {
        detector_object_id: 9,
        receiver: PlanarIrradianceReceiver::try_new(point, Vec3::new(0.0, 0.0, -1.0), 2e-7)
            .unwrap(),
        initial_direction,
        source_transverse_axis: Vec3::new(1.0, 0.0, 0.0),
        initial_medium: 1,
        connection: BundleConnectionOptions {
            scene: CurvedBundleSettings {
                trace: CurvedTraceSettings {
                    max_surface_events: 8,
                    max_distance_m: 8.0,
                    max_total_steps: 16_384,
                    max_query_work: 1_000_000,
                    event: GradientEventOptions {
                        integration: GradientIndexOptions {
                            position_tolerance_m: 1e-11,
                            direction_tolerance: 1e-12,
                            optical_path_tolerance_m: 1e-11,
                            optical_depth_tolerance: 1e-12,
                            relative_tolerance: 1e-11,
                            ..GradientIndexOptions::default()
                        },
                        surface_tolerance_m: 1e-8,
                        ..GradientEventOptions::default()
                    },
                },
                ..CurvedBundleSettings::default()
            },
            ..BundleConnectionOptions::default()
        },
    }
}
fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-10),
        "{actual:.16e} != {expected:.16e}"
    );
}

#[test]
fn g1_connected_eighty_bin_ball_lens_matches_abcd_fresnel_beer_and_cie() {
    with_cx(|cx| {
        let source = fixture([0.0, 0.0, -3.0]);
        let receiver = Point3::new(0.0, 0.0, 2.5);
        let geometry = InstanceScene::try_new(vec![sphere(), plane(9, receiver, 0.2)]).unwrap();
        let air = glass(CauchyIor::try_constant(1.0).unwrap(), 0.0);
        let lens = glass(CauchyIor::try_new(1.45, 0.012, 0.0).unwrap(), 0.15);
        let media = [
            SceneMedium::Homogeneous(&lens),
            SceneMedium::Homogeneous(&air),
        ];
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let controls = request(receiver, Vec3::new(0.0, 0.0, 1.0));
        let result = source.connect_irradiance(cx, &scene, controls).unwrap();
        let mut expected_w = 0.0;
        let mut expected_lux = 0.0;
        for bin in 0..SPECTRAL_BINS {
            let wavelength = wavelength_nm(bin);
            let index = 1.45 + 0.012 / (wavelength * 1e-3).powi(2);
            let footprint = 15.0 / index - 9.5;
            let transmission = 4.0 * index / (1.0 + index).powi(2);
            let expected = source
                .spectral_intensity_w_sr_nm(controls.initial_direction, wavelength)
                .unwrap()
                * transmission.powi(2)
                * (-0.3_f64).exp()
                / footprint.powi(2);
            close(result.spectral_irradiance_w_m2_nm()[bin], expected, 5e-5);
            assert_eq!(result.launch_directions()[bin], controls.initial_direction);
            expected_w += expected * BIN_WIDTH_NM;
            expected_lux += expected * photopic_weight_nm(bin) * LUMINOUS_EFFICACY_LM_W;
        }
        close(result.reading().radiant_irradiance_w_m2, expected_w, 5e-5);
        close(result.reading().illuminance_lux, expected_lux, 5e-5);
        assert!(result.max_endpoint_error_m() <= controls.receiver.position_tolerance_m());
        assert!(result.max_contact_resolution_m() > 0.0);
        assert_eq!(result.trace_evaluations(), SPECTRAL_BINS as u64);
        assert!(result.field_evaluations() >= 5 * SPECTRAL_BINS as u64);
        assert!(result.query_work() > result.field_evaluations());
    });
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

// Quadratic intersections and vector Snell, independent of the native event
// walk, bundle derivatives, connection solver, and production Fresnel routine.
fn ball_reference(ray: Ray, index: f64) -> (Point3, f64, f64, f64) {
    let origin = Vec3::new(ray.origin.x, ray.origin.y, ray.origin.z);
    let b = origin.dot(ray.dir);
    let before = -b - (b * b - origin.dot(origin) + 1.0).sqrt();
    let entry = ray.at(before);
    let normal = unit(Vec3::new(entry.x, entry.y, entry.z));
    let inside = snell(ray.dir, normal, 1.0 / index);
    let chord = -2.0 * inside.dot(Vec3::new(entry.x, entry.y, entry.z));
    let exit = entry.offset(inside.scale(chord));
    let outgoing = snell(inside, unit(Vec3::new(exit.x, exit.y, exit.z)), index);
    let after = (2.5 - exit.z) / outgoing.z;
    let ci = -ray.dir.dot(normal);
    let ct = -inside.dot(normal);
    let ts = 1.0 - ((ci - index * ct) / (ci + index * ct)).powi(2);
    let tp = 1.0 - ((index * ci - ct) / (index * ci + ct)).powi(2);
    (exit.offset(outgoing.scale(after)), chord, ts, tp)
}

fn reference_area(ray: Ray, index: f64) -> f64 {
    let first = unit(add(Vec3::new(1.0, 0.0, 0.0), ray.dir.scale(-ray.dir.x)));
    let axes = [first, cross(ray.dir, first)];
    let h = 1e-6;
    let columns = axes.map(|axis| {
        let endpoint = |angle| {
            ball_reference(
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
fn g1_connected_off_axis_spectrum_solves_each_wavelength_against_independent_snell() {
    with_cx(|cx| {
        let source = fixture([0.03, 0.0, -3.0]);
        let reference_ray = Ray {
            origin: source.position_m(),
            dir: unit(Vec3::new(0.01, 0.0, 1.0)),
        };
        let target = ball_reference(reference_ray, 1.45 + 0.012 / (0.55_f64).powi(2)).0;
        let geometry = InstanceScene::try_new(vec![sphere(), plane(9, target, 0.2)]).unwrap();
        let air = glass(CauchyIor::try_constant(1.0).unwrap(), 0.0);
        let lens = glass(CauchyIor::try_new(1.45, 0.012, 0.0).unwrap(), 0.15);
        let media = [
            SceneMedium::Homogeneous(&lens),
            SceneMedium::Homogeneous(&air),
        ];
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let result = source
            .connect_irradiance(cx, &scene, request(target, reference_ray.dir))
            .unwrap();
        assert!(result.trace_evaluations() > SPECTRAL_BINS as u64);
        assert!(
            (result.launch_directions()[0].x - result.launch_directions()[SPECTRAL_BINS - 1].x)
                .abs()
                > 1e-4
        );
        for bin in 0..SPECTRAL_BINS {
            let wavelength = wavelength_nm(bin);
            let index = 1.45 + 0.012 / (wavelength * 1e-3).powi(2);
            let ray = Ray {
                dir: result.launch_directions()[bin],
                ..reference_ray
            };
            let (endpoint, inside, ts, tp) = ball_reference(ray, index);
            let miss = endpoint.delta_from(target);
            assert!(
                miss.dot(miss).sqrt() < 1e-6,
                "independent spectral miss at bin {bin}: {miss:?}"
            );
            let area = reference_area(ray, index);
            let expected = source
                .spectral_intensity_w_sr_nm(ray.dir, wavelength)
                .unwrap()
                * 0.5
                * (ts * ts + tp * tp)
                * (-0.15 * inside).exp()
                / area;
            close(result.spectral_irradiance_w_m2_nm()[bin], expected, 2e-4);
        }
    });
}

#[test]
fn g1_connected_polarized_fixture_preserves_upstream_losses_and_native_sp_transfer() {
    with_cx(|cx| {
        let source = fixture([0.25, 0.0, -3.0]);
        let ray = Ray {
            origin: source.position_m(),
            dir: unit(Vec3::new(0.05, 0.0, 1.0)),
        };
        let (target, inside, ts, tp) = ball_reference(ray, 1.5);
        let area = reference_area(ray, 1.5);
        let geometry = InstanceScene::try_new(vec![sphere(), plane(9, target, 0.2)]).unwrap();
        let air = glass(CauchyIor::try_constant(1.0).unwrap(), 0.0);
        let lens = glass(CauchyIor::try_constant(1.5).unwrap(), 0.15);
        let media = [
            SceneMedium::Homogeneous(&lens),
            SceneMedium::Homogeneous(&air),
        ];
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let parallel =
            core::array::from_fn(|bin| 0.2 + 0.6 * bin as f64 / (SPECTRAL_BINS - 1) as f64);
        let perpendicular = [0.1; SPECTRAL_BINS];
        let neutral =
            core::array::from_fn(|bin| 0.6 - 0.2 * bin as f64 / (SPECTRAL_BINS - 1) as f64);
        // The projected world-x axis is p polarization for these xz-plane
        // rays. Surviving source power remains relative to the original lamp.
        let filtered = source
            .clone()
            .with_unpolarized_optics()
            .through(&SpectralOptics::linear_diattenuator(0.0, parallel, perpendicular).unwrap())
            .unwrap()
            .through(&SpectralOptics::neutral(neutral).unwrap())
            .unwrap();
        let result = filtered
            .connect_irradiance(cx, &scene, request(target, ray.dir))
            .unwrap();
        let mut expected_lux = 0.0;
        for bin in 0..SPECTRAL_BINS {
            let power =
                0.5 * neutral[bin] * (parallel[bin] * tp * tp + perpendicular[bin] * ts * ts);
            let expected = source
                .spectral_intensity_w_sr_nm(ray.dir, wavelength_nm(bin))
                .unwrap()
                * power
                * (-0.15 * inside).exp()
                / area;
            close(result.spectral_irradiance_w_m2_nm()[bin], expected, 1e-4);
            close(
                filtered.relative_stokes()[bin].i,
                0.5 * neutral[bin] * (parallel[bin] + perpendicular[bin]),
                1e-13,
            );
            expected_lux += expected * photopic_weight_nm(bin) * LUMINOUS_EFFICACY_LM_W;
        }
        close(result.reading().illuminance_lux, expected_lux, 1e-4);
        assert!(filtered.intensity_source().flux_lm() < 0.3 * source.flux_lm());
        assert!(
            (tp * tp - ts * ts).abs() > 0.01,
            "the lens must distinguish source polarization"
        );
    });
}

#[test]
fn g4_failed_later_spectral_bin_returns_no_partial_connected_reading() {
    use crate::gradient_index::bundle::{IndexHessian, IndexJet, SecondOrderIndexField};
    use crate::gradient_index::{GradientIndexError, GradientIndexField, GradientIndexSample};
    use std::sync::atomic::{AtomicU32, Ordering};
    struct LimitedSpectrum {
        successful: AtomicU32,
        failed: AtomicU32,
    }
    impl GradientIndexField for LimitedSpectrum {
        fn sample(
            &self,
            _: Point3,
            wavelength: f64,
        ) -> Result<GradientIndexSample, GradientIndexError> {
            if wavelength >= wavelength_nm(17) {
                self.failed.fetch_add(1, Ordering::Relaxed);
                return Err(GradientIndexError::OutsideFieldDomain);
            }
            self.successful.fetch_add(1, Ordering::Relaxed);
            GradientIndexSample::try_new(1.0, Vec3::new(0.0, 0.0, 0.0), 0.0, Some(1.0))
        }
        fn is_spatially_uniform(&self) -> bool {
            true
        }
    }
    impl SecondOrderIndexField for LimitedSpectrum {
        fn sample_second_order(
            &self,
            p: Point3,
            wavelength: f64,
        ) -> Result<IndexJet, GradientIndexError> {
            Ok(IndexJet::new(
                self.sample(p, wavelength)?,
                IndexHessian::ZERO,
            ))
        }
    }
    with_cx(|cx| {
        let source = fixture([0.0; 3]);
        let target = Point3::new(0.0, 0.0, 2.5);
        let geometry = InstanceScene::try_new(vec![plane(9, target, 1.0)]).unwrap();
        let field = LimitedSpectrum {
            successful: AtomicU32::new(0),
            failed: AtomicU32::new(0),
        };
        let media = [SceneMedium::Differentiable {
            field: &field,
            identity: hash_domain("connected-photometry-test", b"limited-spectrum"),
        }];
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![bindings()[1]], 100_000)
                .unwrap();
        let mut controls = request(target, Vec3::new(0.0, 0.0, 1.0));
        controls.initial_medium = 0;
        assert!(matches!(
            source.connect_irradiance(cx, &scene, controls),
            Err(ConnectedIrradianceError::Connection {
                bin: 17,
                error: BundleConnectionError::Scene(_)
            })
        ));
        assert!(
            field.successful.load(Ordering::Relaxed) >= 34,
            "earlier bins really completed their field evaluations"
        );
        assert_eq!(field.failed.load(Ordering::Relaxed), 1);
        assert_eq!(source.flux_lm(), 900.0);
    });
}
