use super::*;
use crate::charts::{Ray, TriMesh};
use crate::dielectric::{
    BeerLambertAbsorption, CauchyIor, DielectricGlass, DielectricSurface, GlassProvenance,
};
use crate::gradient_index::{GradientIndexOptions, event::GradientEventOptions};
use crate::lighting::EnvironmentMap;
use crate::photometry::scene::{
    NativeGeometryBinding, NativeOpticalBinding, NativeOpticalMedium, NativeOpticalRole,
};
use crate::polarization::dielectric::scene::curved::CurvedTraceSettings;
use crate::polarization::dielectric::scene::curved::camera::{
    RadianceEndpoint, area_lights::CurvedAreaEmitter, raster::CurvedCameraRaster,
};
use crate::spectral::{LiftedSpectrum, SpectralReflectance, cie_x, cie_y, cie_z, y_integral};
use crate::tracer::{Camera, Material, Primitive, RectLight, Scene, Shape};
use fs_blake3::hash_domain;
use fs_exec::{CancelGate, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3, fixtures::SphereChart};

const SPECTRUM: LiftedSpectrum = LiftedSpectrum {
    c: [3.0, -2.0, -0.25],
};

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
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

fn panel(z: f64, rho: f64, strength: f64) -> Primitive {
    Primitive {
        shape: Shape::Mesh(TriMesh::new(
            vec![
                [-4.0, -4.0, z],
                [4.0, -4.0, z],
                [4.0, 4.0, z],
                [-4.0, 4.0, z],
            ],
            vec![[0, 1, 2], [0, 2, 3]],
        )),
        material: Material::SpectralDiffuse {
            reflectance: SpectralReflectance::constant(rho).unwrap(),
        },
        emission: Some((SPECTRUM, strength)),
    }
}

fn scene(primitives: Vec<Primitive>) -> Scene {
    Scene {
        primitives,
        lights: Vec::new(),
        environment: None,
        camera: Camera {
            eye: Point3::new(0.0, 0.0, -3.0),
            forward: Vec3::new(0.0, 0.0, 1.0),
            up: Vec3::new(0.0, 1.0, 0.0),
            half_tan: 0.5,
        },
    }
}

fn binding(id: u64) -> NativeOpticalBinding {
    NativeOpticalBinding {
        geometry: NativeGeometryBinding::Unplaced {
            object_id: id,
            geometry_identity: hash_domain("native-camera-test", &id.to_le_bytes()),
        },
        role: NativeOpticalRole::Material {
            outside: NativeOpticalMedium::Ambient,
        },
    }
}

fn admit(cx: &Cx<'_>, scene: Scene, ambient: DielectricGlass) -> NativeOpticalScene {
    let bindings: Vec<_> = (0..scene.primitives.len())
        .map(|i| binding(i as u64 + 1))
        .collect();
    NativeOpticalScene::try_from_scene(cx, scene, ambient, &bindings, bindings.len()).unwrap()
}

fn settings() -> CurvedTraceSettings {
    CurvedTraceSettings {
        max_surface_events: 12,
        max_distance_m: 12.0,
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

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-10),
        "{actual:.16e} != {expected:.16e}"
    );
}

#[test]
fn g1_original_source_through_glass_retains_beer_and_immersed_camera_index() {
    with_cx(|cx, _| {
        let mut emitter = panel(2.0, 0.0, 2.0);
        // A single original triangle keeps the axial arrival inside one face.
        emitter.shape = Shape::Mesh(TriMesh::new(
            vec![[-4.0, -4.0, 2.0], [4.0, -4.0, 2.0], [0.0, 4.0, 2.0]],
            vec![[0, 1, 2]],
        ));
        let owner = admit(
            cx,
            scene(vec![
                Primitive {
                    shape: Shape::Chart(Box::new(SphereChart {
                        center: Point3::new(0.0, 0.0, 0.0),
                        radius: 1.0,
                    })),
                    material: Material::Dielectric {
                        glass: glass(1.5, 0.2),
                        surface: DielectricSurface::SMOOTH,
                    },
                    emission: None,
                },
                emitter,
            ]),
            glass(1.0, 0.1),
        );
        let calibration = NativeRadianceCalibration::try_new(3.0, 0.0).unwrap();
        let ambient = owner.medium_index(NativeOpticalMedium::Ambient).unwrap();
        let interior = owner
            .medium_index(NativeOpticalMedium::InteriorOf(0))
            .unwrap();
        owner
            .with_curved_camera(cx, 100_000, calibration, |native, camera| {
                assert_eq!(camera.eye(), owner.render_scene().camera.eye);
                assert_eq!(camera.forward(), owner.render_scene().camera.forward);
                assert_eq!(camera.up(), owner.render_scene().camera.up);
                assert_eq!(
                    camera.projection().vertical_half_tan().to_bits(),
                    0.5_f64.to_bits()
                );
                assert!(camera.aperture().is_pinhole());
                for wavelength in [401.25, 550.0, 699.125] {
                    for (origin, medium, depth, index, events) in [
                        (camera.eye(), ambient, 0.7_f64, 1.0_f64, 3),
                        (Point3::new(0.0, 0.0, 0.0), interior, 0.3, 1.5, 2),
                    ] {
                        let mut draws = 0;
                        let sample = native
                            .sample_ray_with_samples(
                                cx,
                                Ray {
                                    origin,
                                    dir: camera.forward(),
                                },
                                wavelength,
                                medium,
                                settings(),
                                |_| {
                                    draws += 1;
                                    [0.0, 0.0, 0.99]
                                },
                            )
                            .unwrap();
                        // Conditional transmission divides its Fresnel proposal;
                        // Beer follows actual lengths and reduced radiance supplies
                        // the launch/source index ratio once, including immersion.
                        let expected =
                            SPECTRUM.eval(wavelength) * 6.0 * (-depth).exp() * index * index;
                        close(sample.radiance_w_m2_sr_nm(), expected, 2e-7);
                        close(sample.path().optical_depth(), depth, 2e-7);
                        assert_eq!(sample.launch_phase_index(), index);
                        assert_eq!(sample.path().surface_events(), events);
                        assert_eq!(draws, events - 1);
                        assert_eq!(
                            sample.endpoint(),
                            RadianceEndpoint::Surface { object_id: 2 }
                        );
                    }
                }
                assert!(
                    native
                        .sample_ray_with_samples(
                            cx,
                            Ray {
                                origin: Point3::new(0.0, 0.0, 0.0),
                                dir: camera.forward()
                            },
                            550.0,
                            ambient,
                            settings(),
                            |_| [0.0, 0.0, 0.99],
                        )
                        .is_err(),
                    "camera containment must not silently replace the declared medium"
                );
            })
            .unwrap();
    });
}

#[test]
fn g1_original_emission_precedes_diffuse_scattering_and_environment_raster_is_replayable() {
    with_cx(|cx, _| {
        for rho in [0.0, 0.4] {
            let mut original = scene(vec![panel(0.0, rho, 4.0)]);
            original.environment = Some(
                EnvironmentMap::try_from_linear_srgb(1, 1, vec![[0.2, 0.7, 1.3]], 0.0).unwrap(),
            );
            let owner = admit(cx, original, glass(1.0, 0.0));
            let calibration = NativeRadianceCalibration::try_new(2.0, 3.0).unwrap();
            assert_eq!(calibration.surface_w_m2_sr_nm_per_unit(), 2.0);
            assert_eq!(calibration.environment_w_m2_sr_nm_per_unit(), 3.0);
            let background = owner
                .render_scene()
                .environment
                .as_ref()
                .unwrap()
                .evaluate(Vec3::new(0.0, 0.0, -1.0))
                .unwrap()
                .emission;
            let expected = |wavelength| {
                SPECTRUM.eval(wavelength) * 8.0
                    + rho * background.0.eval(wavelength) * background.1 * 3.0
            };
            owner
                .with_curved_camera(cx, 100_000, calibration, |native, camera| {
                    // Both sides emit the same original radiance; reflectance
                    // multiplies only the continuation, never the local source.
                    for side in [-1.0, 1.0] {
                        let mut draws = 0;
                        let sample = native
                            .sample_ray_with_samples(
                                cx,
                                Ray {
                                    origin: Point3::new(0.25, -0.25, side * 3.0),
                                    dir: Vec3::new(0.0, 0.0, -side),
                                },
                                611.75,
                                0,
                                settings(),
                                |_| {
                                    draws += 1;
                                    [0.0; 3]
                                },
                            )
                            .unwrap();
                        close(sample.radiance_w_m2_sr_nm(), expected(611.75), 2e-13);
                        assert_eq!(draws, usize::from(rho > 0.0));
                        assert_eq!(
                            sample.endpoint(),
                            if rho > 0.0 {
                                RadianceEndpoint::Environment
                            } else {
                                RadianceEndpoint::Surface { object_id: 1 }
                            }
                        );
                    }
                    let raster = CurvedCameraRaster::try_new(camera, 2, 2).unwrap();
                    for pixel in 0..4 {
                        let sample = raster
                            .sample_pixel(cx, native, 0, settings(), 93, pixel, 7)
                            .unwrap();
                        assert_eq!(
                            sample,
                            raster
                                .sample_pixel(cx, native, 0, settings(), 93, pixel, 7)
                                .unwrap()
                        );
                        let mut xyz = [0.0; 3];
                        for (&wavelength, spectral) in sample
                            .wavelengths_nm()
                            .iter()
                            .zip(sample.spectral_samples())
                        {
                            let radiance = expected(wavelength);
                            close(spectral.radiance_w_m2_sr_nm(), radiance, 2e-13);
                            for (sum, observer) in xyz.iter_mut().zip([
                                cie_x(wavelength),
                                cie_y(wavelength),
                                cie_z(wavelength),
                            ]) {
                                *sum += radiance * (100.0 / y_integral()) * observer;
                            }
                        }
                        for (actual, expected) in sample.xyz().into_iter().zip(xyz) {
                            close(actual, expected, 2e-13);
                        }
                    }
                })
                .unwrap();
        }
    });
}

fn rectangle() -> RectLight {
    RectLight {
        corner: Point3::new(-4.0, -4.0, 0.0),
        edge_u: Vec3::new(8.0, 0.0, 0.0),
        edge_v: Vec3::new(0.0, 8.0, 0.0),
        prim: 0,
        emission: (SPECTRUM, 2.0),
    }
}

#[test]
fn g1_rectangle_proposal_reuses_original_camera_area_source_without_double_emission() {
    with_cx(|cx, _| {
        let mut without_proposal = None;
        for metadata in [false, true] {
            let mut original = scene(vec![panel(0.0, 0.0, 2.0)]);
            if metadata {
                original.lights.push(rectangle());
            }
            let owner = admit(cx, original, glass(1.3, 0.2));
            let sample = owner
                .with_curved_camera(
                    cx,
                    100_000,
                    NativeRadianceCalibration::try_new(3.0, 0.0).unwrap(),
                    |native, _| {
                        let area = CurvedAreaEmitter::from_camera(cx, native, 1, 2).unwrap();
                        assert!(core::ptr::eq(
                            area.source_instance(),
                            &owner.geometry().instances()[0]
                        ));
                        assert_eq!(area.total_area_m2(), 64.0);
                        let source = area.source_emitter().unwrap();
                        assert_eq!(source.object_id(), 1);
                        assert_eq!(source.sides(), EmissionSides::Both);
                        assert_eq!(source.spectrum().lifted_source(), Some((SPECTRUM, 6.0)));
                        native
                            .sample_ray_with_samples(
                                cx,
                                Ray {
                                    origin: Point3::new(0.25, -0.25, -2.0),
                                    dir: Vec3::new(0.0, 0.0, 1.0),
                                },
                                444.125,
                                0,
                                settings(),
                                |_| panic!("black source does not scatter"),
                            )
                            .unwrap()
                    },
                )
                .unwrap();
            close(
                sample.radiance_w_m2_sr_nm(),
                SPECTRUM.eval(444.125) * 6.0 * (-0.4_f64).exp(),
                2e-7,
            );
            if let Some(original) = without_proposal {
                assert_eq!(sample, original);
            } else {
                without_proposal = Some(sample);
            }
        }
    });
}

#[test]
fn g0_original_camera_source_ownership_calibration_and_budgets_refuse_before_observation() {
    for bad in [-1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            NativeRadianceCalibration::try_new(bad, 1.0),
            Err(NativeCameraError::InvalidCalibration)
        );
        assert_eq!(
            NativeRadianceCalibration::try_new(1.0, bad),
            Err(NativeCameraError::InvalidCalibration)
        );
    }
    with_cx(|cx, gate| {
        let calibration = NativeRadianceCalibration::try_new(2.0, 3.0).unwrap();
        for bad in [-1.0, f64::NAN, f64::INFINITY, f64::MAX] {
            let owner = admit(cx, scene(vec![panel(0.0, 0.0, bad)]), glass(1.0, 0.0));
            assert_eq!(
                owner.with_curved_camera(cx, 100_000, calibration, |_, _| ()),
                Err(NativeCameraError::InvalidEmission { primitive: 0 })
            );
        }
        let owner = admit(
            cx,
            scene(vec![panel(0.0, 0.0, f64::from_bits(1))]),
            glass(1.0, 0.0),
        );
        assert_eq!(
            owner.with_curved_camera(
                cx,
                100_000,
                NativeRadianceCalibration::try_new(0.25, 1.0).unwrap(),
                |_, _| ()
            ),
            Err(NativeCameraError::InvalidEmission { primitive: 0 })
        );
        let mut invalid = scene(vec![panel(0.0, 0.0, 2.0)]);
        invalid.primitives[0].emission.as_mut().unwrap().0.c[1] = f64::NAN;
        let owner = admit(cx, invalid, glass(1.0, 0.0));
        assert_eq!(
            owner.with_curved_camera(cx, 100_000, calibration, |_, _| ()),
            Err(NativeCameraError::InvalidEmission { primitive: 0 })
        );
        for mismatch in 0..3 {
            let mut original = scene(vec![panel(0.0, 0.0, 2.0)]);
            let mut light = rectangle();
            match mismatch {
                0 => light.prim = 9,
                1 => light.emission.1 = 3.0,
                _ => light.corner.x += 1.0,
            }
            let primitive = light.prim;
            original.lights.push(light);
            let owner = admit(cx, original, glass(1.0, 0.0));
            assert_eq!(
                owner.with_curved_camera(cx, 100_000, calibration, |_, _| ()),
                Err(NativeCameraError::LightPrimitiveMismatch {
                    light: 0,
                    primitive
                })
            );
        }
        let mut original = scene(vec![panel(0.0, 0.0, 2.0)]);
        original.lights = vec![rectangle(), rectangle()];
        let owner = admit(cx, original, glass(1.0, 0.0));
        assert_eq!(
            owner.with_curved_camera(cx, 100_000, calibration, |_, _| ()),
            Err(NativeCameraError::LightBudgetExceeded {
                lights: 2,
                objects: 1
            })
        );
        let mut original = scene(vec![panel(0.0, 0.0, 2.0)]);
        original.camera.up = original.camera.forward;
        let owner = admit(cx, original, glass(1.0, 0.0));
        assert_eq!(
            owner.with_curved_camera(cx, 100_000, calibration, |_, _| ()),
            Err(NativeCameraError::Camera(
                CameraError::InvalidOrientationBasis
            ))
        );
        let owner = NativeOpticalScene::try_from_scene(
            cx,
            scene(vec![panel(0.0, 0.0, 2.0)]),
            glass(1.0, 0.0),
            &[NativeOpticalBinding {
                role: NativeOpticalRole::Detector {
                    reference_axis: Vec3::new(1.0, 0.0, 0.0),
                },
                ..binding(1)
            }],
            1,
        )
        .unwrap();
        assert_eq!(
            owner.with_curved_camera(cx, 100_000, calibration, |_, _| ()),
            Err(NativeCameraError::Radiance(
                CameraRadianceError::UnsupportedSurface { object_id: 1 }
            ))
        );
        let owner = admit(cx, scene(vec![panel(0.0, 0.0, 0.0)]), glass(1.0, 0.0));
        assert!(matches!(
            owner.with_curved_camera(cx, 0, calibration, |_, _| panic!("no partial view")),
            Err(NativeCameraError::OpticalScene(
                NativeOpticalSceneError::Curved(_)
            ))
        ));
        gate.request();
        assert_eq!(
            owner.with_curved_camera(cx, 100_000, calibration, |_, _| ()),
            Err(NativeCameraError::OpticalScene(
                NativeOpticalSceneError::Cancelled
            ))
        );
    });
}
