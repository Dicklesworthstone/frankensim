//! Analytic native diffuse/refraction fixtures; no measured-material claim.

use super::{CurvedPolarizedScene, CurvedSceneError, CurvedTraceSettings, SceneMedium};
use crate::charts::{Ray, TriMesh};
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricGlass, GlassProvenance};
use crate::gradient_index::event::GradientEventOptions;
use crate::gradient_index::{
    GradientIndexError, GradientIndexField, GradientIndexOptions, GradientIndexSample,
};
use crate::instances::{GeometryInstance, InstanceScene, RigidTransform, SharedGeometry};
use crate::polarization::Stokes;
use crate::polarization::dielectric::scene::curved::bundle::{
    CurvedBundleSceneError, CurvedBundleSettings,
};
use crate::polarization::dielectric::scene::curved::camera::{
    CameraRadianceError, CurvedCameraScene, EmissionSides, RadianceSpectrum, SurfaceEmitter,
};
use crate::polarization::dielectric::scene::{
    Interaction, PathTermination, PolarizedScene, SceneError, SurfaceBinding, TraceSettings,
};
use crate::polarization::dielectric::transport::{BeamMedium, PolarizationFrame, PolarizedBeam};
use crate::spectral::SpectralReflectance;
use fs_blake3::{ContentHash, hash_domain};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 47,
                kernel_id: 11,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate);
    });
}

fn glass(index: f64, absorption: f64) -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(index).unwrap(),
        BeerLambertAbsorption::try_constant(absorption).unwrap(),
        GlassProvenance::Custom,
    )
}

fn plane(id: u64, z: f64, center_x: f64, radius: f64) -> GeometryInstance {
    GeometryInstance::try_new(
        id,
        hash_domain("org.frankensim.test.curved-diffuse", &id.to_le_bytes()),
        SharedGeometry::mesh(TriMesh::new(
            vec![
                [center_x - radius, -radius, z],
                [center_x + radius, -radius, z],
                [center_x, radius, z],
            ],
            vec![[0, 1, 2]],
        )),
        RigidTransform::identity(),
    )
    .unwrap()
}

fn geometry() -> InstanceScene {
    InstanceScene::try_new(vec![plane(1, 1.0, 0.0, 10.0), plane(9, -1.0, 0.0, 10.0)]).unwrap()
}

fn bindings(receiver: Interaction) -> Vec<SurfaceBinding> {
    vec![
        SurfaceBinding {
            object_id: 1,
            interaction: Interaction::DepolarizingDiffuse {
                medium: 0,
                reflectance: 0,
            },
        },
        SurfaceBinding {
            object_id: 9,
            interaction: receiver,
        },
    ]
}

fn ray() -> Ray {
    Ray {
        origin: Point3::new(0.0, 0.0, 0.5),
        dir: Vec3::new(0.0, 0.0, 1.0),
    }
}

fn settings() -> CurvedTraceSettings {
    CurvedTraceSettings {
        max_surface_events: 8,
        max_distance_m: 4.0,
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
    }
}

fn source() -> [SurfaceEmitter; 1] {
    [SurfaceEmitter::new(
        9,
        RadianceSpectrum::try_constant(8.0).unwrap(),
        EmissionSides::Both,
    )]
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}, tolerance {tolerance:.3e}"
    );
}

#[test]
fn g1_gray_white_and_exact_black_have_the_declared_reflected_source_radiance() {
    with_cx(|cx, _| {
        let geometry = geometry();
        let material = glass(1.3, 0.2);
        let media = [SceneMedium::Homogeneous(&material)];
        let sources = source();
        for rho in [0.0, 0.37, 1.0] {
            let reflectances = [SpectralReflectance::constant(rho).unwrap()];
            let scene = CurvedPolarizedScene::try_new_with_materials(
                cx,
                &geometry,
                &media,
                &[],
                &[],
                &reflectances,
                bindings(Interaction::Opaque),
                100_000,
            )
            .unwrap();
            let camera = CurvedCameraScene::try_new(cx, &scene, &sources, None).unwrap();
            // u=0 samples the reflected geometric normal: 0.5 m to the
            // diffuse plane and 2 m to the emitting plane behind the camera.
            let sample = camera
                .sample_ray_with_samples(cx, ray(), 550.0, 0, settings(), |_| {
                    assert!(rho > 0.0, "exact black must not sample a continuation");
                    [0.0, 0.25, 0.5]
                })
                .unwrap();
            close(
                sample.radiance_w_m2_sr_nm(),
                8.0 * rho * (-0.5_f64).exp(),
                2e-6,
            );
            assert_eq!(sample.path().termination(), PathTermination::Opaque);
            assert_eq!(sample.path().medium_index(), 0);
            assert_eq!(
                sample.path().last_contact().unwrap().object_id,
                if rho == 0.0 { 1 } else { 9 }
            );
            if rho == 0.0 {
                assert_eq!(sample.radiance_w_m2_sr_nm(), 0.0);
                assert_eq!(sample.path().surface_events(), 1);
            } else {
                close(sample.path().traveled_m(), 2.5, 2e-7);
                close(sample.path().optical_depth(), 0.5, 5e-8);
                close(sample.path().optical_path_length_m(), 3.25, 3e-7);
                close(
                    sample.path().group_optical_path_length_m().unwrap(),
                    3.25,
                    3e-7,
                );
                assert_eq!(sample.path().surface_events(), 2);
            }
        }
    });
}

#[test]
fn g1_spectral_diffuse_depolarizes_forward_power_in_both_native_scene_owners() {
    with_cx(|cx, _| {
        let geometry = geometry();
        let materials = [glass(1.3, 0.2)];
        let media = [SceneMedium::Homogeneous(&materials[0])];
        let reflectances = [
            SpectralReflectance::table(&[(380.0, 0.2), (550.0, 0.6), (780.0, 0.9)]).unwrap(),
        ];
        let receiver = Interaction::Detector {
            reference_axis: Vec3::new(1.0, 0.0, 0.0),
        };
        let curved = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &[],
            &[],
            &reflectances,
            bindings(receiver),
            100_000,
        )
        .unwrap();
        let straight = PolarizedScene::try_new_with_materials(
            &geometry,
            &materials,
            &[],
            &[],
            &reflectances,
            bindings(receiver),
        )
        .unwrap();
        let frame = PolarizationFrame::try_new(ray().dir, Vec3::new(1.0, 0.0, 0.0)).unwrap();
        for wavelength in [400.0, 550.0, 700.0] {
            let input = PolarizedBeam::try_new(
                ray().origin,
                frame,
                wavelength,
                materials[0],
                Stokes::try_new(3.0, 1.0, -1.0, 1.0).unwrap(),
            )
            .unwrap();
            let curved = curved
                .trace_with_samples(cx, input, 0, settings(), |_| [0.0, 0.25, 0.5])
                .unwrap();
            let straight = straight
                .trace_with_samples(
                    cx,
                    input,
                    0,
                    TraceSettings {
                        max_surface_events: 8,
                        max_distance_m: 4.0,
                        surface_epsilon_m: 1e-9,
                        query_offset_m: 1e-7,
                    },
                    |_| [0.0, 0.25, 0.5],
                )
                .unwrap();
            let expected = 3.0 * reflectances[0].at(wavelength).unwrap() * (-0.5_f64).exp();
            for (termination, traveled, beam) in [
                (curved.termination(), curved.traveled_m(), curved.beam()),
                (
                    straight.termination(),
                    straight.traveled_m(),
                    straight.beam(),
                ),
            ] {
                assert_eq!(termination, PathTermination::Detected);
                close(traveled, 2.5, 2e-7);
                close(beam.power().i, expected, 5e-7);
                assert_eq!([beam.power().q, beam.power().u, beam.power().v], [0.0; 3]);
                assert_eq!(beam.medium(), input.medium());
                close(beam.ray().dir.z, -1.0, 1e-14);
            }
        }
    });
}

struct DispersiveAffineIndex;

impl GradientIndexField for DispersiveAffineIndex {
    fn sample(
        &self,
        p: Point3,
        wavelength: f64,
    ) -> Result<GradientIndexSample, GradientIndexError> {
        let dispersion = 0.035 * (550.0 / wavelength).powi(2);
        let phase = 1.165 + 0.2 * p.x + dispersion;
        // n_g = n - lambda dn/dlambda for this declared Cauchy-like field.
        GradientIndexSample::try_new(
            phase,
            Vec3::new(0.2, 0.0, 0.0),
            0.2,
            Some(phase + 2.0 * dispersion),
        )
    }
}

#[test]
fn g1_diffuse_reflection_keeps_curved_arc_absorption_phase_ratio_and_group_delay() {
    with_cx(|cx, _| {
        let n0 = 1.2_f64;
        let gradient = 0.2_f64;
        let h1 = 1.0_f64;
        let h2 = 1.25_f64;
        // Each normal launch through n(x)=n_start+g*x has
        // n(z)=n_start*cosh(g*z/n_start), s=n_start/g*sinh(g*z/n_start).
        let n1 = n0 * (gradient * h1 / n0).cosh();
        let n2 = n1 * (gradient * h2 / n1).cosh();
        let x1 = (n1 - n0) / gradient;
        let x2 = (n2 - n0) / gradient;
        let length1 = n0 / gradient * (gradient * h1 / n0).sinh();
        let length2 = n1 / gradient * (gradient * h2 / n1).sinh();
        let length = length1 + length2;
        let optical_segment = |n: f64, h: f64| {
            n * h / 2.0 + n * n / (4.0 * gradient) * (2.0 * gradient * h / n).sinh()
        };
        let optical = optical_segment(n0, h1) + optical_segment(n1, h2);
        let geometry =
            InstanceScene::try_new(vec![plane(1, 1.25, x1, 0.02), plane(9, 0.0, x2, 0.02)])
                .unwrap();
        let media = [SceneMedium::GradientIndex {
            field: &DispersiveAffineIndex,
            identity: ContentHash([83; 32]),
        }];
        let rho = 0.4;
        let reflectances = [SpectralReflectance::constant(rho).unwrap()];
        let scene = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &[],
            &[],
            &reflectances,
            bindings(Interaction::Opaque),
            100_000,
        )
        .unwrap();
        let sources = source();
        let camera = CurvedCameraScene::try_new(cx, &scene, &sources, None).unwrap();
        let launch = Ray {
            origin: Point3::new(0.0, 0.0, 0.25),
            ..ray()
        };
        let sample = camera
            .sample_ray_with_samples(cx, launch, 550.0, 0, settings(), |_| [0.0, 0.25, 0.5])
            .unwrap();
        let path = sample.path();
        let expected_power = rho * (-0.2 * length).exp();
        close(path.beam().power().i, expected_power, 3e-7);
        close(
            sample.radiance_w_m2_sr_nm(),
            8.0 * expected_power * (n0 / n2).powi(2),
            3e-6,
        );
        close(path.traveled_m(), length, 5e-7);
        close(path.optical_depth(), 0.2 * length, 1e-7);
        close(path.optical_path_length_m(), optical, 1e-6);
        close(
            path.group_optical_path_length_m().unwrap(),
            optical + 0.07 * length,
            1e-6,
        );
        close(path.beam().ray().origin.x, x2, 5e-7);
        assert_eq!(path.medium_snapshot_identity(), Some(ContentHash([83; 32])));
        let BeamMedium::GradientIndex(local) = path.beam().medium() else {
            panic!("lost local index field")
        };
        close(local.sample().phase_index(), n2, 1e-7);
        close(local.sample().group_index().unwrap(), n2 + 0.07, 1e-7);
        assert_eq!(path.surface_events(), 2);
        assert_eq!(path.last_contact().unwrap().object_id, 9);
        assert!(
            x1 > 0.02 && x2 - x1 > 0.02,
            "straight segments miss both finite targets"
        );
    });
}

#[test]
fn g0_diffuse_continuations_preserve_cumulative_limits_and_cancellation() {
    with_cx(|cx, gate| {
        let geometry = geometry();
        let material = glass(1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&material)];
        let reflectances = [SpectralReflectance::constant(0.8).unwrap()];
        let scene = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &[],
            &[],
            &reflectances,
            bindings(Interaction::Opaque),
            100_000,
        )
        .unwrap();
        let sources = source();
        let camera = CurvedCameraScene::try_new(cx, &scene, &sources, None).unwrap();
        let trace = |limits| {
            camera.sample_ray_with_samples(cx, ray(), 550.0, 0, limits, |_| [0.0, 0.25, 0.5])
        };
        let mut limits = settings();
        limits.max_surface_events = 1;
        assert!(matches!(
            trace(limits),
            Err(CameraRadianceError::Unresolved(
                PathTermination::SurfaceLimit
            ))
        ));
        limits = settings();
        limits.max_distance_m = 1.0;
        assert!(matches!(
            trace(limits),
            Err(CameraRadianceError::Unresolved(
                PathTermination::DistanceLimit
            ))
        ));
        limits = settings();
        limits.max_total_steps = 1;
        assert!(matches!(
            trace(limits),
            Err(CameraRadianceError::Scene(CurvedSceneError::Propagation(
                GradientIndexError::StepBudgetExceeded
            )))
        ));
        limits = settings();
        limits.max_query_work = 1;
        assert!(
            trace(limits).is_err(),
            "diffuse traversal bypassed native geometry work budget"
        );
        gate.request();
        assert!(matches!(
            trace(settings()),
            Err(CameraRadianceError::Scene(CurvedSceneError::Scene(
                SceneError::Cancelled
            )))
        ));
    });
}

#[test]
fn g0_diffuse_surface_cannot_replace_the_arriving_medium() {
    with_cx(|cx, _| {
        let geometry = geometry();
        let first = glass(1.0, 0.0);
        let other = glass(1.4, 0.0);
        let media = [
            SceneMedium::Homogeneous(&first),
            SceneMedium::Homogeneous(&other),
        ];
        let reflectances = [SpectralReflectance::constant(0.5).unwrap()];
        for (medium, reflectance) in [(0, 1), (2, 0)] {
            let mut surfaces = bindings(Interaction::Opaque);
            surfaces[0].interaction = Interaction::DepolarizingDiffuse {
                medium,
                reflectance,
            };
            assert!(matches!(
                CurvedPolarizedScene::try_new_with_materials(
                    cx,
                    &geometry,
                    &media,
                    &[],
                    &[],
                    &reflectances,
                    surfaces.clone(),
                    100_000,
                ),
                Err(CurvedSceneError::Scene(SceneError::InvalidBindings))
            ));
            assert!(matches!(
                PolarizedScene::try_new_with_materials(
                    &geometry,
                    &[first, other],
                    &[],
                    &[],
                    &reflectances,
                    surfaces,
                ),
                Err(SceneError::InvalidBindings)
            ));
        }
        let mut surfaces = bindings(Interaction::Opaque);
        surfaces[0].interaction = Interaction::DepolarizingDiffuse {
            medium: 1,
            reflectance: 0,
        };
        let scene = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &[],
            &[],
            &reflectances,
            surfaces,
            100_000,
        )
        .unwrap();
        let camera = CurvedCameraScene::try_new(cx, &scene, &[], None).unwrap();
        assert!(matches!(
            camera.sample_ray_with_samples(cx, ray(), 550.0, 0, settings(), |_| panic!(
                "medium mismatch reached diffuse sampling"
            )),
            Err(CameraRadianceError::Scene(CurvedSceneError::Scene(
                SceneError::MediumMismatch {
                    object_id: 1,
                    expected: 1,
                    actual: 0,
                }
            )))
        ));
    });
}

#[test]
fn g0_a_diffuse_direction_has_no_invented_smooth_bundle_derivative() {
    with_cx(|cx, _| {
        let geometry = geometry();
        let material = glass(1.2, 0.0);
        let media = [SceneMedium::Homogeneous(&material)];
        let reflectances = [SpectralReflectance::constant(0.5).unwrap()];
        let scene = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &[],
            &[],
            &reflectances,
            bindings(Interaction::Opaque),
            100_000,
        )
        .unwrap();
        let bundle = scene
            .launch_bundle(
                ray(),
                Vec3::new(1.0, 0.0, 0.0),
                550.0,
                Stokes::unpolarized(1.0),
                0,
            )
            .unwrap();
        assert!(matches!(
            scene.trace_bundle_transmitted(
                cx,
                bundle,
                0,
                CurvedBundleSettings {
                    trace: settings(),
                    ..CurvedBundleSettings::default()
                }
            ),
            Err(CurvedBundleSceneError::UnsupportedInteraction { object_id: 1 })
        ));
    });
}
