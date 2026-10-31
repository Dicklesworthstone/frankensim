use super::*;
use crate::charts::TriMesh;
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricSurface, GlassProvenance};
use crate::gradient_index::{GradientIndexOptions, GradientIndexSample};
use crate::instances::{GeometryInstance, RigidTransform, SharedGeometry};
use asupersync::types::Budget;
use fs_blake3::hash_domain;
use fs_exec::{CancelGate, ExecMode, StreamKey};
use fs_geom::Vec3;

mod curved_charts;

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new();
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

fn glass(n: f64, absorption: f64) -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(n).unwrap(),
        BeerLambertAbsorption::try_constant(absorption).unwrap(),
        GlassProvenance::Custom,
    )
}

fn plane(id: u64, z: f64, center_x: f64, radius: f64) -> GeometryInstance {
    GeometryInstance::try_new(
        id,
        hash_domain("curved-scene-test", &id.to_le_bytes()),
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
                ..GradientIndexOptions::default()
            },
            surface_tolerance_m: 1e-7,
            ..GradientEventOptions::default()
        },
    }
}

fn bindings(detector: u64) -> Vec<SurfaceBinding> {
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
            object_id: detector,
            interaction: Interaction::Detector {
                reference_axis: Vec3::new(0.0, 1.0, 0.0),
            },
        },
    ]
}

fn frame() -> PolarizationFrame {
    PolarizationFrame::try_new(Vec3::new(0.0, 0.0, 1.0), Vec3::new(0.0, 1.0, 0.0)).unwrap()
}

struct SlabField;
impl GradientIndexField for SlabField {
    fn sample(&self, p: Point3, _: f64) -> Result<GradientIndexSample, GradientIndexError> {
        GradientIndexSample::try_new(
            1.2 + 0.2 * p.x,
            Vec3::new(0.2, 0.0, 0.0),
            0.3,
            Some(1.2 + 0.2 * p.x),
        )
    }
}

#[test]
fn g1_curved_slab_refracts_to_finite_native_detector_and_closes_optical_lengths() {
    with_cx(|cx, _| {
        // Independent solution for n(x)=n0+g*x with launch along +z:
        // z=(n0/g) asinh(g*s/n0), x=(sqrt(n0²+g²s²)-n0)/g.
        let (n0, g, d, after) = (1.2_f64, 0.2_f64, 1.0_f64, 1.0_f64);
        let s = n0 / g * (g * d / n0).sinh();
        let n_exit = n0 * (g * d / n0).cosh();
        let x_exit = (n_exit - n0) / g;
        let air_x = g * s;
        let air_z = (1.0 - air_x * air_x).sqrt();
        let receiver_x = x_exit + after * air_x / air_z;
        let geometry = InstanceScene::try_new(vec![
            plane(1, d, 0.0, 10.0),
            plane(9, d + after, receiver_x, 0.04),
        ])
        .unwrap();
        let air = glass(1.0, 0.0);
        let field = SlabField;
        let media = [
            SceneMedium::GradientIndex {
                field: &field,
                identity: ContentHash([7; 32]),
            },
            SceneMedium::Homogeneous(&air),
        ];
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(9), 100_000).unwrap();
        let input = scene
            .launch_beam(
                Point3::new(0.0, 0.0, 0.0),
                frame(),
                550.0,
                0,
                Stokes::unpolarized(1.0),
            )
            .unwrap();
        let mut dielectric_events = 0;
        let result = scene
            .trace_with_samples(cx, input, 0, settings(), |_| {
                dielectric_events += 1;
                [0.3, 0.7, 0.999_999]
            })
            .unwrap_or_else(|error| {
                panic!("{error:?} after {dielectric_events} dielectric events")
            });
        assert_eq!(result.termination(), PathTermination::Detected);
        assert_eq!(result.last_contact().unwrap().object_id, 9);
        assert_eq!(result.surface_events(), 2);
        assert_eq!(result.medium_index(), 1);
        assert!((result.beam().ray().origin.x - receiver_x).abs() < 2e-6);
        assert!((result.beam().ray().origin.z - (d + after)).abs() < 2e-6);
        assert!((result.beam().ray().dir.x - air_x).abs() < 2e-6);
        assert!((result.traveled_m() - (s + after / air_z)).abs() < 3e-6);
        let opl = 0.5 * (s * n_exit + n0 * d) + after / air_z;
        assert!((result.optical_path_length_m() - opl).abs() < 4e-6);
        assert_eq!(
            result.group_optical_path_length_m(),
            Some(result.optical_path_length_m())
        );
        assert!((result.optical_depth() - 0.3 * s).abs() < 1e-6);
        assert!((result.beam().power().i - (-0.3 * s).exp()).abs() < 2e-6);
        assert!(result.max_contact_resolution_m() > 0.0);
        // The finite receiver is actually missed by the unrefracted ray.
        assert!(
            geometry.instances()[1]
                .intersect(
                    cx,
                    &crate::charts::Ray {
                        origin: Point3::new(0.0, 0.0, 1.1),
                        dir: Vec3::new(0.0, 0.0, 1.0)
                    },
                    2.0,
                    1e-9
                )
                .unwrap()
                .is_none()
        );
    });
}

#[test]
fn g1_reflection_departs_into_same_field_and_counts_every_metre_of_absorption() {
    with_cx(|cx, _| {
        let geometry =
            InstanceScene::try_new(vec![plane(1, 1.0, 0.0, 10.0), plane(9, -1.0, 0.0, 10.0)])
                .unwrap();
        let medium = glass(1.5, 0.4);
        let uniform = HomogeneousIndexField::new(&medium);
        let air = glass(1.0, 0.0);
        let media = [
            SceneMedium::GradientIndex {
                field: &uniform,
                identity: ContentHash([7; 32]),
            },
            SceneMedium::Homogeneous(&air),
        ];
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(9), 100_000).unwrap();
        let input = scene
            .launch_beam(
                Point3::new(0.0, 0.0, 0.0),
                frame(),
                550.0,
                0,
                Stokes::unpolarized(1.0),
            )
            .unwrap();
        let result = scene
            .trace_with_samples(cx, input, 0, settings(), |_| [0.3, 0.7, 0.0])
            .unwrap();
        assert_eq!(result.termination(), PathTermination::Detected);
        assert_eq!(
            result.medium_snapshot_identity(),
            Some(ContentHash([7; 32]))
        );
        assert_eq!(result.medium_index(), 0);
        assert!((result.traveled_m() - 3.0).abs() < 2e-6);
        assert!((result.optical_depth() - 0.4 * result.traveled_m()).abs() < 1e-12);
        assert!((result.beam().power().i - (-result.optical_depth()).exp()).abs() < 1e-12);
        assert!((result.optical_path_length_m() - 1.5 * result.traveled_m()).abs() < 1e-12);
    });
}

#[test]
fn g1_curved_metal_mirror_uses_local_index_phase_and_entire_absorbing_path() {
    use crate::conductor::{ConductorDataStatus, ConductorIorSample, ConductorSource};
    with_cx(|cx, _| {
        let (n0, g, d) = (1.2_f64, 0.2_f64, 1.0_f64);
        // Unfolding reflection about z=d gives the affine-index solution at 2d.
        let length = n0 / g * (2.0 * g * d / n0).sinh();
        let receiver_x = n0 / g * ((2.0 * g * d / n0).cosh() - 1.0);
        let geometry = InstanceScene::try_new(vec![
            plane(1, d, 0.0, 10.0),
            plane(9, 0.0, receiver_x, 0.04),
        ])
        .unwrap();
        let field = SlabField;
        let identity = ContentHash([8; 32]);
        let media = [SceneMedium::GradientIndex {
            field: &field,
            identity,
        }];
        let mut bindings = bindings(9);
        bindings[0].interaction = Interaction::SmoothConductor {
            optics: 0,
            medium: 0,
        };
        let source = ConductorSource::try_new(
            hash_domain("curved-mirror-test", b"synthetic-n0.2-k3"),
            ConductorDataStatus::Representative,
        )
        .unwrap();
        let metals = [ConductorOptics::try_new(
            core::array::from_fn(|i| {
                ConductorIorSample::try_new(380.0 + 50.0 * i as f64, 0.2, 3.0).unwrap()
            }),
            source,
        )
        .unwrap()];
        assert!(matches!(
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings.clone(), 100_000),
            Err(CurvedSceneError::Scene(SceneError::InvalidBindings))
        ));
        let scene = CurvedPolarizedScene::try_new_with_conductors(
            cx, &geometry, &media, &metals, bindings, 100_000,
        )
        .unwrap();
        let input = scene
            .launch_beam(
                Point3::new(0.0, 0.0, 0.0),
                frame(),
                550.0,
                0,
                Stokes::try_new(1.0, 0.2, 0.6, 0.3).unwrap(),
            )
            .unwrap();
        let result = scene
            .trace_with_samples(cx, input, 0, settings(), |_| {
                panic!("smooth conductor reflection needs no random draw")
            })
            .unwrap();
        assert_eq!(result.termination(), PathTermination::Detected);
        assert_eq!(result.last_contact().unwrap().object_id, 9);
        assert_eq!(result.surface_events(), 2);
        assert_eq!(result.medium_index(), 0);
        assert_eq!(result.medium_snapshot_identity(), Some(identity));
        assert!((result.beam().ray().origin.x - receiver_x).abs() < 2e-6);
        assert!(result.beam().ray().origin.z.abs() < 2e-6);
        assert!((result.beam().ray().dir.x - (2.0 * g * d / n0).tanh()).abs() < 2e-6);
        assert!((result.beam().ray().dir.z + 1.0 / (2.0 * g * d / n0).cosh()).abs() < 2e-6);
        assert!((result.traveled_m() - length).abs() < 2e-6);
        assert!((result.optical_depth() - 0.3 * length).abs() < 1e-6);
        let BeamMedium::GradientIndex(local) = result.beam().medium() else {
            panic!("mirror lost its gradient-index medium")
        };
        let final_index = n0 * (2.0 * g * d / n0).cosh();
        assert!((local.sample().phase_index() - final_index).abs() < 1e-6);
        assert!(
            (result.optical_path_length_m() - 0.5 * (length * final_index + 2.0 * n0 * d)).abs()
                < 3e-6
        );
        // Independent complex-field Fresnel evaluation for N=0.2+3i,
        // n_m=n0*cosh(g*d/n0), c_m=sech(g*d/n0):
        // q=sqrt(N²-n_m²(1-c_m²)), r_s=(n_m*c_m-q)/(n_m*c_m+q),
        // r_p=(N²*c_m-n_m*q)/(N²*c_m+n_m*q). Apply these amplitudes
        // to input (1,.2,.6,.3), then exp(-0.3*length) on every component.
        // These numbers use the actual mirror index 1.2167052826537559;
        // substituting vacuum or the launch index 1.2 fails the comparison.
        let expected = [
            0.494_865_469_422_422_5,
            0.099_675_321_446_841_12,
            -0.300_063_902_369_324_7,
            -0.141_765_970_292_339_26,
        ];
        let power = result.beam().power();
        for (actual, expected) in [power.i, power.q, power.u, power.v]
            .into_iter()
            .zip(expected)
        {
            assert!((actual - expected).abs() < 2e-6, "{actual} != {expected}");
        }
    });
}

#[test]
fn g0_curved_path_limits_and_stale_medium_refuse_without_darkening_survivors() {
    with_cx(|cx, _| {
        let geometry =
            InstanceScene::try_new(vec![plane(1, 1.0, 0.0, 10.0), plane(9, 2.0, 0.0, 10.0)])
                .unwrap();
        let field = SlabField;
        let air = glass(1.0, 0.0);
        let media = [
            SceneMedium::GradientIndex {
                field: &field,
                identity: ContentHash([7; 32]),
            },
            SceneMedium::Homogeneous(&air),
        ];
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(9), 100_000).unwrap();
        let input = scene
            .launch_beam(
                Point3::new(0.0, 0.0, 0.0),
                frame(),
                550.0,
                0,
                Stokes::unpolarized(1.0),
            )
            .unwrap();
        let mut limits = settings();
        limits.max_distance_m = 0.5;
        let result = scene.trace(cx, input, 0, limits, 1, 1).unwrap();
        assert_eq!(result.termination(), PathTermination::DistanceLimit);
        assert!((result.traveled_m() - 0.5).abs() < 1e-14);
        assert!(result.beam().power().i > 0.0);
        assert_eq!(result.surface_events(), 0);
        assert!(matches!(
            scene.trace(cx, input, 1, limits, 1, 1),
            Err(CurvedSceneError::Scene(SceneError::InvalidInitialMedium))
        ));
        limits = settings();
        limits.max_total_steps = 1;
        assert!(scene.trace(cx, input, 0, limits, 1, 1).is_err());
    });
}

#[test]
fn g4_cancellation_after_interface_discards_the_complete_curved_path() {
    with_cx(|cx, gate| {
        let geometry =
            InstanceScene::try_new(vec![plane(1, 1.0, 0.0, 10.0), plane(9, 2.0, 0.0, 10.0)])
                .unwrap();
        let medium = glass(1.2, 0.0);
        let air = glass(1.0, 0.0);
        let media = [
            SceneMedium::Homogeneous(&medium),
            SceneMedium::Homogeneous(&air),
        ];
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(9), 100_000).unwrap();
        let input = scene
            .launch_beam(
                Point3::new(0.0, 0.0, 0.0),
                frame(),
                550.0,
                0,
                Stokes::unpolarized(1.0),
            )
            .unwrap();
        let result = scene.trace_with_samples(cx, input, 0, settings(), |_| {
            gate.request();
            [0.3, 0.7, 0.999_999]
        });
        assert!(matches!(
            result,
            Err(CurvedSceneError::Scene(SceneError::Cancelled))
        ));
        assert_eq!(input.power(), Stokes::unpolarized(1.0));
    });
}

#[cfg(feature = "photometry")]
#[test]
fn g1_curved_native_detector_feeds_original_spectral_meter() {
    use crate::photometry::source::{SPECTRAL_BINS, wavelength_nm};
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![
            plane(1, 1.0, 0.0, 10.0),
            plane(5, -1.0, 0.0, 10.0),
            plane(9, 2.0, 0.0, 10.0),
        ])
        .unwrap();
        let field = SlabField;
        let air = glass(1.0, 0.0);
        let media = [
            SceneMedium::GradientIndex {
                field: &field,
                identity: ContentHash([7; 32]),
            },
            SceneMedium::Homogeneous(&air),
        ];
        let mut surfaces = bindings(9);
        surfaces.push(SurfaceBinding {
            object_id: 5,
            interaction: Interaction::Opaque,
        });
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, surfaces, 100_000).unwrap();
        let mut meter = scene.power_meter(9).unwrap();
        let s = 1.2 / 0.2 * (0.2_f64 / 1.2).sinh();
        let n_exit = 1.2 * (0.2_f64 / 1.2).cosh();
        let ci = 1.2 / n_exit;
        let ct = (1.0 - (0.2 * s).powi(2)).sqrt();
        let rs = ((n_exit * ci - ct) / (n_exit * ci + ct)).powi(2);
        let rp = ((ci - n_exit * ct) / (ci + n_exit * ct)).powi(2);
        let reflection = 0.5 * (rs + rp);
        let transmission = 1.0 - reflection;
        for bin in 0..SPECTRAL_BINS {
            // Two explicit branch strata. Their weights include the branch's
            // probability mass and the factor two canceled by the meter's
            // all-launch mean; reflected launches are recorded too.
            for (event, weight) in [(0.0, 2.0 * reflection), (0.999_999, 2.0 * transmission)] {
                let input = scene
                    .launch_beam(
                        Point3::new(0.0, 0.0, 0.0),
                        frame(),
                        wavelength_nm(bin),
                        0,
                        Stokes::unpolarized(weight),
                    )
                    .unwrap();
                let path = scene
                    .trace_with_samples(cx, input, 0, settings(), |_| [0.3, 0.7, event])
                    .unwrap();
                assert_eq!(
                    path.termination(),
                    if event == 0.0 {
                        PathTermination::Opaque
                    } else {
                        PathTermination::Detected
                    }
                );
                meter.record_curved(bin, &path).unwrap();
            }
            assert!(
                (meter.spectral_mean_w_nm(bin).unwrap().i - transmission * (-0.3 * s).exp()).abs()
                    < 2e-6
            );
        }
        let reading = meter.reading().unwrap().require_complete().unwrap();
        let expected_i = 400.0 * transmission * (-0.3 * s).exp();
        assert!((reading.radiant_stokes_w().i - expected_i).abs() < 0.001);
        let expected_q = expected_i * (rp - rs) / (2.0 - rs - rp);
        assert!((reading.radiant_stokes_w().q - expected_q).abs() < 0.001);
        assert!(reading.max_contact_resolution_m() > 0.0);
        assert_eq!(
            reading
                .linear_analyzer(0.0)
                .unwrap()
                .max_contact_resolution_m,
            reading.max_contact_resolution_m()
        );
        assert_eq!(reading.detected_paths(), SPECTRAL_BINS as u64);
        assert_eq!(reading.launched_paths(), 2 * SPECTRAL_BINS as u64);
    });
}
