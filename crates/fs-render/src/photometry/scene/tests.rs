use super::*;
use crate::charts::TriMesh;
use crate::conductor::{CoherentDielectricCoating, ConductorSurface};
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricSurface, GlassProvenance};
use crate::gradient_index::{GradientIndexOptions, event::GradientEventOptions};
use crate::photometry::source::bundle::PlanarIrradianceReceiver;
use crate::photometry::source::bundle::connection::ConnectedIrradianceRequest;
use crate::photometry::source::{
    AngularDistribution, BIN_WIDTH_NM, LUMINOUS_EFFICACY_LM_W, PhotometricSource, SPECTRAL_BINS,
    SourceSpectrum, photopic_weight_nm, wavelength_nm,
};
use crate::polarization::Stokes;
use crate::polarization::dielectric::scene::curved::CurvedTraceSettings;
use crate::polarization::dielectric::scene::curved::bundle::CurvedBundleSettings;
use crate::polarization::dielectric::scene::curved::bundle::connection::BundleConnectionOptions;
use crate::polarization::dielectric::scene::{PathTermination, SceneError};
use crate::polarization::dielectric::transport::PolarizationFrame;
use crate::spectral::LiftedSpectrum;
use crate::tracer::{Camera, RectLight};
use fs_blake3::hash_domain;
use fs_exec::{CancelGate, ExecMode, StreamKey};
use fs_geom::{Point3, fixtures::SphereChart};

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

fn diffuse(rho: f64) -> Material {
    Material::SpectralDiffuse {
        reflectance: SpectralReflectance::constant(rho).unwrap(),
    }
}

fn sphere(radius: f64, material: Material) -> Primitive {
    Primitive {
        shape: Shape::Chart(Box::new(SphereChart {
            center: Point3::new(0.0, 0.0, 0.0),
            radius,
        })),
        material,
        emission: None,
    }
}

fn plane(z: f64, material: Material) -> Primitive {
    Primitive {
        shape: Shape::Mesh(TriMesh::new(
            vec![[-4.0, -4.0, z], [4.0, -4.0, z], [0.0, 4.0, z]],
            vec![[0, 1, 2]],
        )),
        material,
        emission: None,
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

fn raw(id: u64, role: NativeOpticalRole) -> NativeOpticalBinding {
    NativeOpticalBinding {
        geometry: NativeGeometryBinding::Unplaced {
            object_id: id,
            geometry_identity: hash_domain("native-optical-scene-test", &id.to_le_bytes()),
        },
        role,
    }
}

fn material(outside: NativeOpticalMedium) -> NativeOpticalRole {
    NativeOpticalRole::Material { outside }
}
fn detector() -> NativeOpticalRole {
    NativeOpticalRole::Detector {
        reference_axis: Vec3::new(1.0, 0.0, 0.0),
    }
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

fn frame() -> PolarizationFrame {
    PolarizationFrame::try_new(Vec3::new(0.0, 0.0, 1.0), Vec3::new(1.0, 0.0, 0.0)).unwrap()
}

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-10),
        "{actual:.16e} != {expected:.16e}"
    );
}

#[test]
fn g1_original_chart_ball_lens_connects_calibrated_spectrum_without_geometry_substitution() {
    with_cx(|cx, _| {
        let lens = DielectricGlass::new(
            CauchyIor::try_new(1.45, 0.012, 0.0).unwrap(),
            BeerLambertAbsorption::try_constant(0.15).unwrap(),
            GlassProvenance::Custom,
        );
        let owner = NativeOpticalScene::try_from_scene(
            cx,
            scene(vec![
                sphere(
                    1.0,
                    Material::Dielectric {
                        glass: lens,
                        surface: DielectricSurface::SMOOTH,
                    },
                ),
                plane(2.5, diffuse(0.3)),
            ]),
            glass(1.0, 0.0),
            &[
                raw(1, material(NativeOpticalMedium::Ambient)),
                raw(9, detector()),
            ],
            2,
        )
        .unwrap();
        let source = PhotometricSource::try_new(
            RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [0.0, 0.0, -3.0]).unwrap(),
            1200.0,
            0.75,
            AngularDistribution::cosine_power(2.0, core::f64::consts::FRAC_PI_2).unwrap(),
            SourceSpectrum::equal_energy(),
        )
        .unwrap();
        let request = ConnectedIrradianceRequest {
            detector_object_id: owner.surface_binding(1).unwrap().object_id,
            receiver: PlanarIrradianceReceiver::try_new(
                Point3::new(0.0, 0.0, 2.5),
                Vec3::new(0.0, 0.0, -1.0),
                2e-7,
            )
            .unwrap(),
            initial_direction: Vec3::new(0.0, 0.0, 1.0),
            source_transverse_axis: Vec3::new(1.0, 0.0, 0.0),
            initial_medium: owner.medium_index(NativeOpticalMedium::Ambient).unwrap(),
            connection: BundleConnectionOptions {
                scene: CurvedBundleSettings {
                    trace: settings(),
                    ..CurvedBundleSettings::default()
                },
                ..BundleConnectionOptions::default()
            },
        };
        let result = owner
            .with_curved_scene(cx, 100_000, |native| {
                source.connect_irradiance(cx, native, request)
            })
            .unwrap()
            .unwrap();
        // Independent on-axis ray-transfer derivative for a radius-one ball:
        // source is two metres before the first vertex, detector 1.5 after exit.
        let mut radiant = 0.0;
        let mut illuminance = 0.0;
        for bin in 0..SPECTRAL_BINS {
            let wavelength = wavelength_nm(bin);
            let index = 1.45 + 0.012 / (wavelength * 1e-3).powi(2);
            let footprint = 15.0 / index - 9.5;
            let transmission = 4.0 * index / (1.0 + index).powi(2);
            let expected = source
                .spectral_intensity_w_sr_nm(request.initial_direction, wavelength)
                .unwrap()
                * transmission.powi(2)
                * (-0.3_f64).exp()
                / footprint.powi(2);
            close(result.spectral_irradiance_w_m2_nm()[bin], expected, 5e-5);
            radiant += expected * BIN_WIDTH_NM;
            illuminance += expected * photopic_weight_nm(bin) * LUMINOUS_EFFICACY_LM_W;
        }
        close(result.reading().radiant_irradiance_w_m2, radiant, 5e-5);
        close(result.reading().illuminance_lux, illuminance, 5e-5);
        assert!(result.query_work() > 0);
        assert_eq!(result.trace_evaluations(), SPECTRAL_BINS as u64);
        assert!(matches!(
            owner.geometry().instances()[0].geometry(),
            SharedGeometry::Chart(_)
        ));
    });
}

fn nested(cx: &Cx<'_>, inner_outside: NativeOpticalMedium) -> NativeOpticalScene {
    NativeOpticalScene::try_from_scene(
        cx,
        scene(vec![
            sphere(
                2.0,
                Material::Dielectric {
                    glass: glass(1.3, 0.2),
                    surface: DielectricSurface::SMOOTH,
                },
            ),
            sphere(
                1.0,
                Material::Dielectric {
                    glass: glass(1.3, 0.7),
                    surface: DielectricSurface::SMOOTH,
                },
            ),
            plane(3.0, diffuse(0.4)),
        ]),
        glass(1.3, 0.05),
        &[
            raw(3, material(NativeOpticalMedium::Ambient)),
            raw(2, material(inner_outside)),
            raw(1, detector()),
        ],
        3,
    )
    .unwrap()
}

#[test]
fn g1_embedded_source_and_equal_index_nested_media_keep_distinct_absorption_and_topology() {
    with_cx(|cx, _| {
        let owner = nested(cx, NativeOpticalMedium::InteriorOf(0));
        let inner = owner
            .medium_index(NativeOpticalMedium::InteriorOf(1))
            .unwrap();
        let outer = owner
            .medium_index(NativeOpticalMedium::InteriorOf(0))
            .unwrap();
        assert_ne!(inner, outer);
        assert_ne!(
            outer,
            owner.medium_index(NativeOpticalMedium::Ambient).unwrap()
        );
        let result = owner
            .with_curved_scene(cx, 100_000, |native| {
                let beam = native
                    .launch_beam(
                        Point3::new(0.0, 0.0, 0.0),
                        frame(),
                        550.0,
                        inner,
                        Stokes::try_new(1.0, 0.2, 0.1, 0.3).unwrap(),
                    )
                    .unwrap();
                native.trace_with_samples(cx, beam, inner, settings(), |_| [0.4, 0.6, 0.9])
            })
            .unwrap()
            .unwrap();
        assert_eq!(result.termination(), PathTermination::Detected);
        assert_eq!(result.medium_index(), 0);
        assert_eq!(result.last_contact().unwrap().object_id, 1);
        close(result.traveled_m(), 3.0, 2e-7);
        close(result.optical_depth(), 0.7 + 0.2 + 0.05, 2e-7);
        close(result.beam().power().i, (-0.95_f64).exp(), 2e-7);
        close(result.beam().power().v, 0.3 * (-0.95_f64).exp(), 2e-7);

        // Acyclic declarations alone cannot prove geometry containment. An
        // incorrect exterior must still fail at the actual outer boundary,
        // even though all three numerical phase indices are exactly equal.
        let wrong = nested(cx, NativeOpticalMedium::Ambient);
        let result = wrong
            .with_curved_scene(cx, 100_000, |native| {
                let beam = native
                    .launch_beam(
                        Point3::new(0.0, 0.0, 0.0),
                        frame(),
                        550.0,
                        inner,
                        Stokes::try_new(1.0, 0.0, 0.0, 0.0).unwrap(),
                    )
                    .unwrap();
                native.trace_with_samples(cx, beam, inner, settings(), |_| [0.4, 0.6, 0.9])
            })
            .unwrap();
        assert!(
            matches!(result, Err(CurvedSceneError::Scene(SceneError::MediumMismatch { object_id: 3, expected, actual: 0 })) if expected == outer)
        );
    });
}

#[test]
fn g1_original_black_occluder_blocks_curved_transport_without_drawing_a_scatter() {
    with_cx(|cx, _| {
        let owner = NativeOpticalScene::try_from_scene(
            cx,
            scene(vec![plane(1.0, diffuse(0.0)), plane(2.0, diffuse(0.5))]),
            glass(1.0, 0.0),
            &[
                raw(7, material(NativeOpticalMedium::Ambient)),
                raw(9, detector()),
            ],
            2,
        )
        .unwrap();
        let result = owner
            .with_curved_scene(cx, 100_000, |native| {
                let beam = native
                    .launch_beam(
                        Point3::new(0.0, 0.0, 0.0),
                        frame(),
                        550.0,
                        0,
                        Stokes::try_new(1.0, 0.0, 0.0, 0.0).unwrap(),
                    )
                    .unwrap();
                native.trace_with_samples(cx, beam, 0, settings(), |_| {
                    panic!("exact black consumes no random draw")
                })
            })
            .unwrap()
            .unwrap();
        assert_eq!(result.termination(), PathTermination::Opaque);
        assert_eq!(result.last_contact().unwrap().object_id, 7);
        close(result.traveled_m(), 1.0, 2e-7);
    });
}

#[test]
fn g1_original_instance_placement_geometry_optics_and_renderer_emission_are_retained() {
    with_cx(|cx, _| {
        let geometry = SharedGeometry::chart(SphereChart {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 0.3,
        });
        let pose = RigidTransform::try_new([0.0, 0.0, 0.6, 0.8], [2.0, 1.0, 0.0]).unwrap();
        let instance = GeometryInstance::try_new(
            91,
            hash_domain("native-optical-scene-test", b"placed"),
            geometry.clone(),
            pose,
        )
        .unwrap();
        let identity = instance.frame_identity();
        let bare = ConductorOptics::representative_tungsten();
        let coating = CoherentDielectricCoating::try_new(
            CauchyIor::try_constant(1.45).unwrap(),
            125e-9,
            bare.source(),
        )
        .unwrap();
        let optics = bare.with_coating(coating);
        let original_material = Material::Conductor {
            optics,
            surface: ConductorSurface::try_rough(0.25).unwrap(),
        };
        let emission = (
            LiftedSpectrum {
                c: [0.1, -0.2, 0.5],
            },
            13.0,
        );
        let light = RectLight {
            corner: Point3::new(-1.0, -1.0, 3.0),
            edge_u: Vec3::new(2.0, 0.0, 0.0),
            edge_v: Vec3::new(0.0, 2.0, 0.0),
            prim: 1,
            emission,
        };
        let mut original = scene(vec![
            Primitive {
                shape: Shape::Instance(instance),
                material: original_material,
                emission: Some(emission),
            },
            plane(3.0, diffuse(0.5)),
        ]);
        original.lights.push(light);
        let owner = NativeOpticalScene::try_from_scene(
            cx,
            original,
            glass(1.0, 0.0),
            &[
                NativeOpticalBinding {
                    geometry: NativeGeometryBinding::OriginalInstance,
                    role: material(NativeOpticalMedium::Ambient),
                },
                raw(5, detector()),
            ],
            2,
        )
        .unwrap();
        let retained = owner
            .geometry()
            .instances()
            .iter()
            .find(|shape| shape.object_id() == 91)
            .unwrap();
        assert!(retained.geometry().ptr_eq(&geometry));
        assert_eq!(retained.frame_identity(), identity);
        assert_eq!(retained.transform(), pose);
        assert_eq!(owner.conductors[0], optics);
        assert!(
            matches!(owner.surface_binding(0).unwrap().interaction, Interaction::RoughConductor { optics: 0, medium: 0, surface } if surface == ConductorSurface::try_rough(0.25).unwrap())
        );
        for (index, primitive) in owner.render_scene().primitives.iter().enumerate() {
            let Shape::Instance(render_instance) = &primitive.shape else {
                panic!("static shape should be shared")
            };
            let optical = owner
                .geometry()
                .instances()
                .iter()
                .find(|shape| shape.object_id() == owner.surface_binding(index).unwrap().object_id)
                .unwrap();
            assert!(render_instance.geometry().ptr_eq(optical.geometry()));
            assert_eq!(render_instance.frame_identity(), optical.frame_identity());
        }
        assert_eq!(
            owner.render_scene().primitives[0].material,
            original_material
        );
        assert_eq!(owner.render_scene().primitives[0].emission, Some(emission));
        assert_eq!(owner.render_scene().lights, [light]);
        assert_eq!(owner.render_scene().camera.half_tan, 0.5);
        let restored = owner.into_scene();
        assert_eq!(restored.primitives[0].emission, Some(emission));
        assert_eq!(restored.lights, [light]);
    });
}

#[test]
fn g0_invalid_topology_materials_and_budgets_refuse_before_observations() {
    with_cx(|cx, gate| {
        let air = glass(1.0, 0.0);
        let dielectric = Material::Dielectric {
            glass: glass(1.5, 0.0),
            surface: DielectricSurface::SMOOTH,
        };
        for topology in [
            [
                NativeOpticalMedium::InteriorOf(0),
                NativeOpticalMedium::Ambient,
            ],
            [
                NativeOpticalMedium::InteriorOf(1),
                NativeOpticalMedium::InteriorOf(0),
            ],
        ] {
            assert!(matches!(
                NativeOpticalScene::try_from_scene(
                    cx,
                    scene(vec![sphere(2.0, dielectric), sphere(1.0, dielectric)]),
                    air,
                    &[raw(1, material(topology[0])), raw(2, material(topology[1]))],
                    2,
                ),
                Err(NativeOpticalSceneError::CyclicMedia { .. })
            ));
        }
        assert!(matches!(
            NativeOpticalScene::try_from_scene(
                cx,
                scene(vec![sphere(1.0, dielectric), plane(2.0, diffuse(0.2))]),
                air,
                &[
                    raw(1, material(NativeOpticalMedium::InteriorOf(1))),
                    raw(2, detector())
                ],
                2,
            ),
            Err(NativeOpticalSceneError::InvalidMediumReference {
                medium: NativeOpticalMedium::InteriorOf(1)
            })
        ));
        assert!(matches!(
            NativeOpticalScene::try_from_scene(
                cx,
                scene(vec![plane(
                    1.0,
                    Material::Ggx {
                        reflectance: LiftedSpectrum { c: [1.0, 0.0, 0.0] },
                        alpha: 0.2,
                    }
                )]),
                air,
                &[raw(1, material(NativeOpticalMedium::Ambient))],
                1,
            ),
            Err(NativeOpticalSceneError::UnsupportedMaterial { primitive: 0, .. })
        ));
        assert!(matches!(
            NativeOpticalScene::try_from_scene(
                cx,
                scene(vec![plane(1.0, diffuse(0.3)), plane(2.0, diffuse(0.3))]),
                air,
                &[
                    raw(1, material(NativeOpticalMedium::Ambient)),
                    raw(2, detector())
                ],
                1,
            ),
            Err(NativeOpticalSceneError::ObjectBudgetExceeded {
                objects: 2,
                max_objects: 1
            })
        ));
        let owner = NativeOpticalScene::try_from_scene(
            cx,
            scene(vec![plane(
                1.0,
                Material::Lambertian {
                    reflectance: LiftedSpectrum { c: [0.0, 0.0, 0.3] },
                },
            )]),
            air,
            &[raw(1, material(NativeOpticalMedium::Ambient))],
            1,
        )
        .unwrap();
        assert!(owner.reflectances.is_empty());
        assert_eq!(
            owner.surface_binding(0).unwrap().interaction,
            Interaction::DepolarizingLiftedDiffuse {
                medium: 0,
                reflectance: LiftedSpectrum { c: [0.0, 0.0, 0.3] }
            }
        );
        assert!(matches!(
            owner.with_curved_scene(cx, 0, |_| panic!(
                "build budget must refuse before observation"
            )),
            Err(NativeOpticalSceneError::Curved(_))
        ));
        gate.request();
        assert_eq!(
            owner.with_curved_scene(cx, 100_000, |_| ()),
            Err(NativeOpticalSceneError::Cancelled)
        );
    });
}

#[test]
fn g1_original_nonconstant_lambertian_is_shared_without_a_spectral_table_fit() {
    with_cx(|cx, _| {
        let reflectance = LiftedSpectrum {
            c: [2.3, -1.7, 0.2],
        };
        let original_material = Material::Lambertian { reflectance };
        let owner = NativeOpticalScene::try_from_scene(
            cx,
            scene(vec![plane(1.0, original_material)]),
            glass(1.0, 0.0),
            &[raw(1, material(NativeOpticalMedium::Ambient))],
            1,
        )
        .unwrap();
        assert!(owner.reflectances.is_empty());
        assert_eq!(
            owner.render_scene().primitives[0].material,
            original_material
        );
        assert_eq!(
            owner.surface_binding(0).unwrap().interaction,
            Interaction::DepolarizingLiftedDiffuse {
                medium: 0,
                reflectance
            }
        );
        owner
            .with_curved_scene(cx, 100_000, |native| {
                for wavelength in [387.5, 552.5, 772.5] {
                    let beam = native
                        .launch_beam(
                            Point3::new(0.0, 0.0, 0.0),
                            frame(),
                            wavelength,
                            0,
                            Stokes::try_new(1.0, 0.2, -0.1, 0.3).unwrap(),
                        )
                        .unwrap();
                    let mut limits = settings();
                    limits.max_surface_events = 1;
                    let path = native
                        .trace_with_samples(cx, beam, 0, limits, |_| [0.0; 3])
                        .unwrap();
                    assert_eq!(path.termination(), PathTermination::SurfaceLimit);
                    close(path.beam().power().i, reflectance.eval(wavelength), 2e-12);
                    assert_eq!(path.beam().power().q, 0.0);
                    assert_eq!(path.beam().power().u, 0.0);
                    assert_eq!(path.beam().power().v, 0.0);
                }
            })
            .unwrap();
        assert!(matches!(
            NativeOpticalScene::try_from_scene(
                cx,
                scene(vec![plane(
                    1.0,
                    Material::Lambertian {
                        reflectance: LiftedSpectrum {
                            c: [f64::NAN, 0.0, 0.0]
                        }
                    }
                )]),
                glass(1.0, 0.0),
                &[raw(1, material(NativeOpticalMedium::Ambient))],
                1,
            ),
            Err(NativeOpticalSceneError::InvalidMaterial { primitive: 0 })
        ));
    });
}
