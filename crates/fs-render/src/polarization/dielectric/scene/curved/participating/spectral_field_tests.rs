use super::*;
use crate::charts::{Ray, TriMesh};
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricGlass, GlassProvenance};
use crate::gradient_index::event::GradientEventOptions;
use crate::gradient_index::{
    GradientIndexError, GradientIndexField, GradientIndexOptions, GradientIndexSample,
};
use crate::instances::{GeometryInstance, InstanceScene, RigidTransform, SharedGeometry};
use crate::polarization::Stokes;
use crate::polarization::dielectric::scene::curved::camera::participating::CurvedRayleighCameraScene;
use crate::polarization::dielectric::scene::curved::camera::{
    CameraRadianceError, CurvedCameraScene, EmissionSides, RadianceSpectrum, SurfaceEmitter,
};
use crate::polarization::dielectric::scene::{Interaction, PathTermination, SurfaceBinding};
use crate::polarization::dielectric::transport::PolarizationFrame;
use fs_blake3::hash_domain;
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_geom::Vec3;
use std::sync::atomic::{AtomicUsize, Ordering};

const FIELD_ID: ContentHash = ContentHash([143; 32]);

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 31,
                kernel_id: 7,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate);
    });
}

fn geometry(z: f64) -> InstanceScene {
    InstanceScene::try_new(vec![
        GeometryInstance::try_new(
            9,
            hash_domain(
                "org.frankensim.test.curved-spectral-molecules",
                &z.to_le_bytes(),
            ),
            SharedGeometry::mesh(TriMesh::new(
                vec![[-10.0, -10.0, z], [10.0, -10.0, z], [0.0, 10.0, z]],
                vec![[0, 1, 2]],
            )),
            RigidTransform::identity(),
        )
        .unwrap(),
    ])
    .unwrap()
}

fn bindings() -> Vec<SurfaceBinding> {
    vec![SurfaceBinding {
        object_id: 9,
        interaction: Interaction::Opaque,
    }]
}

fn emitters() -> [SurfaceEmitter; 1] {
    [SurfaceEmitter::new(
        9,
        RadianceSpectrum::try_constant(4.0).unwrap(),
        EmissionSides::Both,
    )]
}

fn ray() -> Ray {
    Ray {
        origin: Point3::new(0.0, 0.0, 0.0),
        dir: Vec3::new(0.0, 0.0, 1.0),
    }
}

fn launch(scene: &CurvedPolarizedScene<'_>, wavelength: f64) -> PolarizedBeam {
    scene
        .launch_beam(
            ray().origin,
            PolarizationFrame::try_new(ray().dir, Vec3::new(0.0, 1.0, 0.0)).unwrap(),
            wavelength,
            0,
            Stokes::unpolarized(1.0),
        )
        .unwrap()
}

fn settings() -> CurvedRayleighTraceSettings {
    CurvedRayleighTraceSettings {
        surfaces: CurvedTraceSettings {
            max_surface_events: 8,
            max_distance_m: 4.0,
            max_total_steps: 16_384,
            max_query_work: 1_000_000,
            event: GradientEventOptions {
                integration: GradientIndexOptions {
                    position_tolerance_m: 1e-11,
                    direction_tolerance: 1e-12,
                    ..GradientIndexOptions::default()
                },
                surface_tolerance_m: 1e-7,
                ..GradientEventOptions::default()
            },
        },
        max_volume_events: 8,
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

struct Index {
    gradient: f64,
    absorption: f64,
}
impl GradientIndexField for Index {
    fn sample(&self, position: Point3, _: f64) -> Result<GradientIndexSample, GradientIndexError> {
        let index = 1.2 + self.gradient * position.x;
        GradientIndexSample::try_new(
            index,
            Vec3::new(self.gradient, 0.0, 0.0),
            self.absorption,
            Some(index),
        )
    }
}

#[derive(Clone, Copy)]
struct LocalLaw {
    identity: ContentHash,
    majorant: f64,
    local: f64,
    delta: f64,
}
impl CurvedRayleighField for LocalLaw {
    fn source_field_identity(&self) -> ContentHash {
        self.identity
    }
    fn scattering_majorant_per_m(&self, _: f64) -> Result<f64, VolumeTransportError> {
        Ok(self.majorant)
    }
    fn sample_scattering(
        &self,
        _: Point3,
        _: f64,
    ) -> Result<CurvedRayleighSample, VolumeTransportError> {
        Ok(CurvedRayleighSample {
            scattering_per_m: self.local,
            depolarization_ratio: self.delta,
        })
    }
}

#[test]
fn g0_spectral_scattering_requires_the_original_refractive_snapshot() {
    with_cx(|cx, _| {
        let geometry = geometry(1.0);
        let field = Index {
            gradient: 0.0,
            absorption: 0.0,
        };
        let media = [SceneMedium::GradientIndex {
            field: &field,
            identity: FIELD_ID,
        }];
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let mut law = LocalLaw {
            identity: FIELD_ID,
            majorant: 1.0,
            local: 1.0,
            delta: 0.03,
        };
        assert!(CurvedRayleighScene::try_new_with_fields(cx, &scene, &[]).is_err());
        let volumes = [Some(CurvedRayleighBinding::Field(&law))];
        assert!(CurvedRayleighScene::try_new_with_fields(cx, &scene, &volumes).is_ok());
        for identity in [ContentHash([0; 32]), ContentHash([144; 32])] {
            law.identity = identity;
            let volumes = [Some(CurvedRayleighBinding::Field(&law))];
            assert!(CurvedRayleighScene::try_new_with_fields(cx, &scene, &volumes).is_err());
        }
        law.identity = FIELD_ID;
        let material = DielectricGlass::new(
            CauchyIor::try_constant(1.2).unwrap(),
            BeerLambertAbsorption::try_constant(0.0).unwrap(),
            GlassProvenance::Custom,
        );
        let homogeneous = [SceneMedium::Homogeneous(&material)];
        let scene = CurvedPolarizedScene::try_new(cx, &geometry, &homogeneous, bindings(), 100_000)
            .unwrap();
        let volumes = [Some(CurvedRayleighBinding::Field(&law))];
        assert!(matches!(
            CurvedRayleighScene::try_new_with_fields(cx, &scene, &volumes),
            Err(CurvedSceneError::Scene(SceneError::InvalidBindings))
        ));
    });
}

#[test]
fn g3_mixed_binding_preserves_original_density_transport_bit_for_bit() {
    with_cx(|cx, _| {
        let geometry = geometry(5.0);
        let field = Index {
            gradient: 0.4,
            absorption: 0.3,
        };
        let media = [SceneMedium::GradientIndex {
            field: &field,
            identity: FIELD_ID,
        }];
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let spectrum = MassScatteringSpectrum::rayleigh(550.0, 1.0).unwrap();
        let phase = PhaseFunction::rayleigh_depolarized(0.03).unwrap();
        let law = CurvedRayleighMedium::try_new(&spectrum, 2.0, &phase).unwrap();
        let original = [Some(law)];
        let mixed = [Some(CurvedRayleighBinding::Density(law))];
        let before = CurvedRayleighScene::try_new(cx, &scene, &original).unwrap();
        let after = CurvedRayleighScene::try_new_with_fields(cx, &scene, &mixed).unwrap();
        let mut limited = settings();
        limited.max_volume_events = 1;
        for wavelength in [400.0, 550.0, 700.0] {
            let beam = launch(&scene, wavelength);
            let old = before.trace(cx, beam, 0, limited, 31, 42).unwrap();
            let new = after.trace(cx, beam, 0, limited, 31, 42).unwrap();
            assert_eq!(old, new);
        }
    });
}

#[test]
fn g3_zero_field_majorant_preserves_camera_bits_and_consumes_no_flight() {
    with_cx(|cx, _| {
        let geometry = geometry(1.0);
        let field = Index {
            gradient: 0.2,
            absorption: 0.3,
        };
        let media = [SceneMedium::GradientIndex {
            field: &field,
            identity: FIELD_ID,
        }];
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let sources = emitters();
        let camera = CurvedCameraScene::try_new(cx, &scene, &sources, None).unwrap();
        let law = LocalLaw {
            identity: FIELD_ID,
            majorant: 0.0,
            local: 0.0,
            delta: 0.0,
        };
        let volumes = [Some(CurvedRayleighBinding::Field(&law))];
        let transport = CurvedRayleighScene::try_new_with_fields(cx, &scene, &volumes).unwrap();
        let participating = CurvedRayleighCameraScene::try_new(cx, &camera, &transport).unwrap();
        let before = camera
            .sample_ray_with_samples(cx, ray(), 550.0, 0, settings().surfaces, |_| {
                panic!("opaque source draws no surface branch")
            })
            .unwrap();
        let after = participating
            .sample_ray_with_samples(
                cx,
                ray(),
                550.0,
                0,
                settings(),
                |_| panic!("opaque source draws no surface branch"),
                |_| panic!("zero majorant draws no flight"),
            )
            .unwrap();
        assert_eq!(*after.radiance_sample(), before);
        assert_eq!(after.free_flight_draws(), 0);
        assert_eq!(after.volume_events(), 0);
    });
}

struct DispersiveLaw;
impl CurvedRayleighField for DispersiveLaw {
    fn source_field_identity(&self) -> ContentHash {
        FIELD_ID
    }
    fn scattering_majorant_per_m(&self, wavelength: f64) -> Result<f64, VolumeTransportError> {
        let ratio = 550.0 / wavelength;
        Ok(2.0 * ratio.powi(4) * (1.0 + 0.1 * ratio * ratio))
    }
    fn sample_scattering(
        &self,
        position: Point3,
        wavelength: f64,
    ) -> Result<CurvedRayleighSample, VolumeTransportError> {
        Ok(CurvedRayleighSample {
            scattering_per_m: self.scattering_majorant_per_m(wavelength)?,
            depolarization_ratio: 0.03 + 0.02 * position.z + 0.01 * (550.0 / wavelength).powi(2),
        })
    }
}

#[test]
fn g1_spectral_rate_and_local_phase_are_sampled_on_the_analytic_curved_arc() {
    with_cx(|cx, _| {
        let geometry = geometry(5.0);
        let field = Index {
            gradient: 0.4,
            absorption: 0.3,
        };
        let media = [SceneMedium::GradientIndex {
            field: &field,
            identity: FIELD_ID,
        }];
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let volumes = [Some(CurvedRayleighBinding::Field(&DispersiveLaw))];
        let transport = CurvedRayleighScene::try_new_with_fields(cx, &scene, &volumes).unwrap();
        let mut limited = settings();
        limited.max_volume_events = 1;
        for wavelength in [400.0_f64, 550.0, 700.0] {
            let ratio = 550.0 / wavelength;
            let rate = 2.0 * ratio.powi(4) * (1.0 + 0.1 * ratio * ratio);
            let arc = 1.4 / rate;
            let n = (1.2_f64.powi(2) + (0.4 * arc).powi(2)).sqrt();
            let x = (n - 1.2) / 0.4;
            let z = 3.0 * (arc / 3.0).asinh();
            let delta = 0.03 + 0.02 * z + 0.01 * ratio * ratio;
            let isotropic = 3.0 * delta / (2.0 + delta);
            let right_angle = isotropic + 0.5 * (1.0 - isotropic);
            let output = transport
                .trace_with_tracking_samples(
                    cx,
                    launch(&scene, wavelength),
                    0,
                    limited,
                    |_| panic!("collision precedes the surface"),
                    |_| [1.0 - (-1.4_f64).exp(), right_angle, 0.0, 0.9],
                )
                .unwrap();
            let path = output.path();
            assert_eq!(path.termination(), PathTermination::VolumeLimit);
            assert_eq!(path.medium_snapshot_identity(), Some(FIELD_ID));
            close(path.traveled_m(), arc, 2e-10);
            close(path.beam().ray().origin.x, x, 2e-9);
            close(path.beam().ray().origin.z, z, 2e-9);
            close(path.optical_depth(), 0.3 * arc, 2e-10);
            let power = path.beam().power();
            close(power.i, (-0.3 * arc).exp(), 3e-9);
            close(power.q / power.i, (1.0 - delta) / (1.0 + delta), 3e-9);
            assert_eq!(output.volume_events(), 1);
            assert_eq!(output.free_flight_draws(), 1);
        }
    });
}

struct PositionPhase;
impl CurvedRayleighField for PositionPhase {
    fn source_field_identity(&self) -> ContentHash {
        FIELD_ID
    }
    fn scattering_majorant_per_m(&self, _: f64) -> Result<f64, VolumeTransportError> {
        Ok(1.0)
    }
    fn sample_scattering(
        &self,
        position: Point3,
        _: f64,
    ) -> Result<CurvedRayleighSample, VolumeTransportError> {
        Ok(CurvedRayleighSample {
            scattering_per_m: 1.0,
            depolarization_ratio: 0.02 + 0.12 * position.x,
        })
    }
}

#[test]
fn g1_changing_local_molecular_phase_reaches_camera_radiance_with_full_stokes() {
    with_cx(|cx, _| {
        let geometry = geometry(1.0);
        let field = Index {
            gradient: 0.0,
            absorption: 0.2,
        };
        let media = [SceneMedium::GradientIndex {
            field: &field,
            identity: FIELD_ID,
        }];
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let sources = emitters();
        let camera = CurvedCameraScene::try_new(cx, &scene, &sources, None).unwrap();
        let volumes = [Some(CurvedRayleighBinding::Field(&PositionPhase))];
        let transport = CurvedRayleighScene::try_new_with_fields(cx, &scene, &volumes).unwrap();
        let participating = CurvedRayleighCameraScene::try_new(cx, &camera, &transport).unwrap();
        let right_angle = |delta| {
            let isotropic = 3.0 * delta / (2.0 + delta);
            isotropic + 0.5 * (1.0 - isotropic) + f64::EPSILON
        };
        let output = participating
            .sample_ray_with_tracking_samples(
                cx,
                ray(),
                550.0,
                0,
                settings(),
                |_| panic!("no dielectric branch"),
                |flight| match flight {
                    0 => [1.0 - (-0.25_f64).exp(), right_angle(0.02), 0.0, 0.9],
                    1 => [1.0 - (-0.5_f64).exp(), right_angle(0.08), 0.5, 0.9],
                    2 => [1.0 - 1e-12, 0.5, 0.0, 0.9],
                    _ => panic!("finite source precedes another collision"),
                },
            )
            .unwrap();
        // Two common-plane 90-degree collisions with DIFFERENT local phases:
        // scalar response is 1 + rho_1 rho_2, not the product of their M00s.
        let rho1 = (1.0 - 0.02) / (1.0 + 0.02);
        let rho2 = (1.0 - 0.08) / (1.0 + 0.08);
        close(
            output.radiance_sample().radiance_w_m2_sr_nm(),
            4.0 * (1.0 + rho1 * rho2) * (-0.3_f64).exp(),
            4e-6,
        );
        close(output.radiance_sample().path().traveled_m(), 1.5, 2e-7);
        assert_eq!(output.volume_events(), 2);
        assert_eq!(output.free_flight_draws(), 3);
    });
}

struct CountedLaw {
    law: LocalLaw,
    majorant_queries: AtomicUsize,
    local_queries: AtomicUsize,
}
impl CurvedRayleighField for CountedLaw {
    fn source_field_identity(&self) -> ContentHash {
        self.law.identity
    }
    fn scattering_majorant_per_m(&self, wavelength: f64) -> Result<f64, VolumeTransportError> {
        self.majorant_queries.fetch_add(1, Ordering::Relaxed);
        self.law.scattering_majorant_per_m(wavelength)
    }
    fn sample_scattering(
        &self,
        position: Point3,
        wavelength: f64,
    ) -> Result<CurvedRayleighSample, VolumeTransportError> {
        self.local_queries.fetch_add(1, Ordering::Relaxed);
        self.law.sample_scattering(position, wavelength)
    }
}

#[test]
fn g1_null_candidates_reuse_the_sampled_majorant_without_extra_camera_weights() {
    with_cx(|cx, _| {
        let geometry = geometry(1.0);
        let field = Index {
            gradient: 0.0,
            absorption: 0.2,
        };
        let media = [SceneMedium::GradientIndex {
            field: &field,
            identity: FIELD_ID,
        }];
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let sources = emitters();
        let camera = CurvedCameraScene::try_new(cx, &scene, &sources, None).unwrap();
        let law = CountedLaw {
            law: LocalLaw {
                identity: FIELD_ID,
                majorant: 2.0,
                local: 0.5,
                delta: 0.03,
            },
            majorant_queries: AtomicUsize::new(0),
            local_queries: AtomicUsize::new(0),
        };
        let volumes = [Some(CurvedRayleighBinding::Field(&law))];
        let transport = CurvedRayleighScene::try_new_with_fields(cx, &scene, &volumes).unwrap();
        let participating = CurvedRayleighCameraScene::try_new(cx, &camera, &transport).unwrap();
        let output = participating
            .sample_ray_with_tracking_samples(
                cx,
                ray(),
                550.0,
                0,
                settings(),
                |_| panic!("no dielectric branch"),
                |flight| match flight {
                    0 | 1 => [1.0 - (-0.4_f64).exp(), 0.5, 0.0, 0.9],
                    2 => [1.0 - 1e-12, 0.5, 0.0, 0.9],
                    _ => panic!("finite source precedes another candidate"),
                },
            )
            .unwrap();
        close(
            output.radiance_sample().radiance_w_m2_sr_nm(),
            4.0 * (-0.2_f64).exp(),
            2e-7,
        );
        assert_eq!(output.volume_events(), 0);
        assert_eq!(output.free_flight_draws(), 3);
        assert_eq!(law.majorant_queries.load(Ordering::Relaxed), 1);
        assert_eq!(law.local_queries.load(Ordering::Relaxed), 2);
    });
}

#[test]
fn g0_invalid_spectral_rates_phases_and_missing_independent_draws_refuse() {
    with_cx(|cx, _| {
        let geometry = geometry(1.0);
        let field = Index {
            gradient: 0.0,
            absorption: 0.0,
        };
        let media = [SceneMedium::GradientIndex {
            field: &field,
            identity: FIELD_ID,
        }];
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let beam = launch(&scene, 550.0);
        for (majorant, local, delta, error) in [
            (-1.0, 0.0, 0.03, VolumeTransportError::ScatteringSpectrum),
            (
                f64::NAN,
                0.0,
                0.03,
                VolumeTransportError::ScatteringSpectrum,
            ),
            (
                f64::INFINITY,
                0.0,
                0.03,
                VolumeTransportError::ScatteringSpectrum,
            ),
            (1.0, -0.1, 0.03, VolumeTransportError::ScatteringSpectrum),
            (1.0, 1.01, 0.03, VolumeTransportError::ScatteringSpectrum),
            (
                1.0,
                f64::NAN,
                0.03,
                VolumeTransportError::ScatteringSpectrum,
            ),
            (1.0, 0.0, f64::NAN, VolumeTransportError::Phase),
            (1.0, 1.0, 6.0 / 7.0, VolumeTransportError::Phase),
        ] {
            let law = LocalLaw {
                identity: FIELD_ID,
                majorant,
                local,
                delta,
            };
            let volumes = [Some(CurvedRayleighBinding::Field(&law))];
            let transport = CurvedRayleighScene::try_new_with_fields(cx, &scene, &volumes).unwrap();
            assert_eq!(
                transport
                    .trace_with_tracking_samples(
                        cx,
                        beam,
                        0,
                        settings(),
                        |_| [0.5; 3],
                        |_| [0.1, 0.5, 0.0, 0.0],
                    )
                    .unwrap_err(),
                volume_error(error)
            );
        }
        let law = LocalLaw {
            identity: FIELD_ID,
            majorant: 1.0,
            local: 1.0,
            delta: 0.03,
        };
        let volumes = [Some(CurvedRayleighBinding::Field(&law))];
        let transport = CurvedRayleighScene::try_new_with_fields(cx, &scene, &volumes).unwrap();
        assert_eq!(
            transport
                .trace_with_samples(
                    cx,
                    beam,
                    0,
                    settings(),
                    |_| panic!("refusal precedes geometry"),
                    |_| panic!("three variates cannot sample a local field")
                )
                .unwrap_err(),
            volume_error(VolumeTransportError::Sample)
        );
        assert_eq!(
            transport
                .trace_with_tracking_samples(
                    cx,
                    beam,
                    0,
                    settings(),
                    |_| [0.5; 3],
                    |_| [0.1, 0.5, 0.0, f64::NAN]
                )
                .unwrap_err(),
            volume_error(VolumeTransportError::Sample)
        );
    });
}

#[test]
fn g0_native_volume_steps_distance_and_cancellation_remain_camera_refusals() {
    with_cx(|cx, gate| {
        let geometry = geometry(1.0);
        let field = Index {
            gradient: 0.0,
            absorption: 0.0,
        };
        let media = [SceneMedium::GradientIndex {
            field: &field,
            identity: FIELD_ID,
        }];
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let sources = emitters();
        let camera = CurvedCameraScene::try_new(cx, &scene, &sources, None).unwrap();
        let law = LocalLaw {
            identity: FIELD_ID,
            majorant: 1.0,
            local: 1.0,
            delta: 0.03,
        };
        let volumes = [Some(CurvedRayleighBinding::Field(&law))];
        let transport = CurvedRayleighScene::try_new_with_fields(cx, &scene, &volumes).unwrap();
        let participating = CurvedRayleighCameraScene::try_new(cx, &camera, &transport).unwrap();
        let mut limited = settings();
        limited.max_volume_events = 1;
        assert!(matches!(
            participating.sample_ray_with_tracking_samples(
                cx,
                ray(),
                550.0,
                0,
                limited,
                |_| [0.5; 3],
                |_| [0.1, 0.5, 0.0, 0.0],
            ),
            Err(CameraRadianceError::Unresolved(
                PathTermination::VolumeLimit
            ))
        ));
        limited = settings();
        limited.surfaces.max_total_steps = 1;
        assert!(matches!(
            participating.sample_ray_with_tracking_samples(
                cx,
                ray(),
                550.0,
                0,
                limited,
                |_| [0.5; 3],
                |_| [0.99, 0.5, 0.0, 0.0],
            ),
            Err(CameraRadianceError::Scene(CurvedSceneError::Propagation(
                GradientIndexError::StepBudgetExceeded
            )))
        ));
        let vacuum = LocalLaw { local: 0.0, ..law };
        let volumes = [Some(CurvedRayleighBinding::Field(&vacuum))];
        let vacuum_transport =
            CurvedRayleighScene::try_new_with_fields(cx, &scene, &volumes).unwrap();
        let vacuum_camera =
            CurvedRayleighCameraScene::try_new(cx, &camera, &vacuum_transport).unwrap();
        let away = Ray {
            dir: Vec3::new(0.0, 0.0, -1.0),
            ..ray()
        };
        assert!(matches!(
            vacuum_camera.sample_ray_with_tracking_samples(
                cx,
                away,
                550.0,
                0,
                settings(),
                |_| [0.5; 3],
                |_| [0.9, 0.5, 0.0, 0.0],
            ),
            Err(CameraRadianceError::Unresolved(
                PathTermination::DistanceLimit
            ))
        ));
        gate.request();
        assert!(matches!(
            participating.sample_ray_with_tracking_samples(
                cx,
                ray(),
                550.0,
                0,
                settings(),
                |_| panic!("cancelled surface sample"),
                |_| panic!("cancelled flight sample"),
            ),
            Err(CameraRadianceError::Scene(CurvedSceneError::Scene(
                SceneError::Cancelled
            )))
        ));
    });
}
