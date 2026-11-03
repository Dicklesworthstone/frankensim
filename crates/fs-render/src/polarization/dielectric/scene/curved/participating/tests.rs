use super::super::{CurvedPolarizedScene, CurvedSceneError, CurvedTraceSettings, SceneMedium};
use super::{
    CurvedDensityField, CurvedRayleighBinding, CurvedRayleighEvent, CurvedRayleighField,
    CurvedRayleighMedium, CurvedRayleighSample, CurvedRayleighScene, CurvedRayleighTraceSettings,
    CurvedTraceEvent,
};
use crate::charts::TriMesh;
use crate::dielectric::{
    BeerLambertAbsorption, CauchyIor, DielectricGlass, DielectricSurface, GlassProvenance,
};
use crate::gradient_index::event::GradientEventOptions;
use crate::gradient_index::{
    GradientIndexError, GradientIndexField, GradientIndexOptions, GradientIndexSample,
};
use crate::instances::{GeometryInstance, InstanceScene, RigidTransform, SharedGeometry};
use crate::polarization::Stokes;
use crate::polarization::dielectric::scene::{
    Interaction, PathTermination, SceneError, SurfaceBinding,
};
use crate::polarization::dielectric::transport::{BeamMedium, PolarizationFrame};
use crate::volumes::transport::scattering::{MassScatteringSpectrum, PhaseFunction};
use fs_blake3::{ContentHash, hash_domain};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(
            &Cx::new(
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
            ),
            &gate,
        );
    });
}

fn glass(index: f64, absorption: f64) -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(index).unwrap(),
        BeerLambertAbsorption::try_constant(absorption).unwrap(),
        GlassProvenance::Custom,
    )
}

fn plane(id: u64, axis: usize, position: f64) -> GeometryInstance {
    let vertices = if axis == 0 {
        vec![
            [position, -10.0, -10.0],
            [position, 10.0, -10.0],
            [position, 0.0, 10.0],
        ]
    } else {
        vec![
            [-10.0, -10.0, position],
            [10.0, -10.0, position],
            [0.0, 10.0, position],
        ]
    };
    GeometryInstance::try_new(
        id,
        hash_domain("org.frankensim.test.curved-scattering", &id.to_le_bytes()),
        SharedGeometry::mesh(TriMesh::new(vertices, vec![[0, 1, 2]])),
        RigidTransform::identity(),
    )
    .unwrap()
}

fn detector(object_id: u64) -> SurfaceBinding {
    SurfaceBinding {
        object_id,
        interaction: Interaction::Detector {
            reference_axis: Vec3::new(0.0, 1.0, 0.0),
        },
    }
}

fn frame() -> PolarizationFrame {
    PolarizationFrame::try_new(Vec3::new(0.0, 0.0, 1.0), Vec3::new(0.0, 1.0, 0.0)).unwrap()
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
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.3e}"
    );
}

struct LinearIndex;
impl GradientIndexField for LinearIndex {
    fn sample(&self, position: Point3, _: f64) -> Result<GradientIndexSample, GradientIndexError> {
        let index = 1.2 + 0.4 * position.x;
        GradientIndexSample::try_new(index, Vec3::new(0.4, 0.0, 0.0), 0.3, Some(index))
    }
}

#[test]
fn g1_spectral_free_flights_follow_analytic_curved_arc_and_physical_stokes_power() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(9, 2, 5.0)]).unwrap();
        let field = LinearIndex;
        let media = [SceneMedium::GradientIndex {
            field: &field,
            identity: ContentHash([71; 32]),
        }];
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![detector(9)], 100_000)
                .unwrap();
        let spectrum = MassScatteringSpectrum::rayleigh(550.0, 1.0).unwrap();
        let delta = 0.03;
        let phase = PhaseFunction::rayleigh_depolarized(delta).unwrap();
        let volumes = [Some(
            CurvedRayleighMedium::try_new(&spectrum, 2.0, &phase).unwrap(),
        )];
        let scattering = CurvedRayleighScene::try_new(cx, &scene, &volumes).unwrap();
        let mut limit = settings();
        limit.max_volume_events = 1;
        for wavelength in [400.0_f64, 550.0, 700.0] {
            // Independent exponential CDF, exact lambda^-4 rate, and analytic
            // ray for n(x)=n0+g*x. Collision distance is ARC length, not z or chord.
            let sigma = 2.0 * (550.0 / wavelength).powi(4);
            let arc = 1.4 / sigma;
            let (n0, g) = (1.2_f64, 0.4_f64);
            let local_index = (n0 * n0 + g * g * arc * arc).sqrt();
            let x = (local_index - n0) / g;
            let z = n0 / g * (g * arc / n0).asinh();
            let opl = 0.5 * (arc * local_index + n0 * n0 / g * (g * arc / n0).asinh());
            let isotropic_fraction = 3.0 * delta / (2.0 + delta);
            let right_angle_sample = isotropic_fraction + 0.5 * (1.0 - isotropic_fraction);
            let input = scene
                .launch_beam(
                    Point3::new(0.0, 0.0, 0.0),
                    frame(),
                    wavelength,
                    0,
                    Stokes::unpolarized(1.0),
                )
                .unwrap();
            let output = scattering
                .trace_with_samples(
                    cx,
                    input,
                    0,
                    limit,
                    |_| panic!("collision occurs before a surface"),
                    |flight| {
                        assert_eq!(flight, 0);
                        [1.0 - (-1.4_f64).exp(), right_angle_sample, 0.0]
                    },
                )
                .unwrap();
            let path = output.path();
            assert_eq!(path.termination(), PathTermination::VolumeLimit);
            assert_eq!(output.volume_events(), 1);
            assert_eq!(output.free_flight_draws(), 1);
            assert_eq!(path.surface_events(), 0);
            assert!(path.last_contact().is_none());
            close(path.traveled_m(), arc, 2e-10);
            close(path.beam().ray().origin.x, x, 2e-9);
            close(path.beam().ray().origin.y, 0.0, 2e-10);
            close(path.beam().ray().origin.z, z, 2e-9);
            close(path.optical_path_length_m(), opl, 3e-9);
            close(path.group_optical_path_length_m().unwrap(), opl, 3e-9);
            close(path.optical_depth(), 0.3 * arc, 2e-10);
            // Mueller/PDF has I=1, Q=(1-d)/(1+d) at 90 degrees for unpolarized
            // input. Only absorption attenuates: no n^2 or extra exp(-sigma_s*s).
            let power = (-0.3 * arc).exp();
            close(path.beam().power().i, power, 3e-9);
            close(
                path.beam().power().q,
                power * (1.0 - delta) / (1.0 + delta),
                3e-9,
            );
            close(path.beam().power().u, 0.0, 2e-10);
            close(path.beam().power().v, 0.0, 2e-10);
            close(path.beam().ray().dir.x, n0 / local_index, 3e-9);
            close(path.beam().ray().dir.y, 0.0, 2e-10);
            close(path.beam().ray().dir.z, -g * arc / local_index, 3e-9);
            let BeamMedium::GradientIndex(local) = path.beam().medium() else {
                panic!("GRIN identity lost")
            };
            assert_eq!(local.position(), path.beam().ray().origin);
            close(local.sample().phase_index(), local_index, 3e-9);
            assert_eq!(path.medium_snapshot_identity(), Some(ContentHash([71; 32])));
        }
    });
}

#[test]
fn g1_scattering_can_return_from_a_straight_miss_and_retains_dipole_extinction() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(9, 0, 1.0)]).unwrap();
        let material = glass(1.4, 0.2);
        let media = [SceneMedium::Homogeneous(&material)];
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![detector(9)], 100_000)
                .unwrap();
        let spectrum = MassScatteringSpectrum::constant(0.7).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let volumes = [Some(
            CurvedRayleighMedium::try_new(&spectrum, 1.0, &phase).unwrap(),
        )];
        let scattering = CurvedRayleighScene::try_new(cx, &scene, &volumes).unwrap();
        let input = scene
            .launch_beam(
                Point3::new(0.0, 0.0, 0.0),
                frame(),
                550.0,
                0,
                Stokes::unpolarized(1.0),
            )
            .unwrap();
        assert_eq!(
            scene
                .trace(cx, input, 0, settings().surfaces, 31, 42)
                .unwrap()
                .termination(),
            PathTermination::Escaped
        );
        let output = scattering
            .trace_with_samples(
                cx,
                input,
                0,
                settings(),
                |_| panic!("detector consumes no boundary sample"),
                |flight| match flight {
                    0 => [1.0 - (-0.35_f64).exp(), 0.5, 0.0],
                    1 => [1.0 - 1e-12, 0.5, 0.0],
                    _ => panic!("receiver must precede another flight"),
                },
            )
            .unwrap();
        let path = output.path();
        assert_eq!(path.termination(), PathTermination::Detected);
        assert_eq!(path.last_contact().unwrap().object_id, 9);
        assert_eq!(path.surface_events(), 1);
        assert_eq!(output.volume_events(), 1);
        assert_eq!(output.free_flight_draws(), 2);
        close(path.traveled_m(), 1.5, 2e-7);
        close(path.beam().ray().origin.x, 1.0, 2e-7);
        close(path.beam().ray().origin.z, 0.5, 2e-9);
        close(path.optical_path_length_m(), 1.4 * 1.5, 3e-7);
        close(path.beam().power().i, (-0.3_f64).exp(), 3e-7);
        close(path.beam().power().q, (-0.3_f64).exp(), 3e-7);
        close(path.beam().power().u, 0.0, 1e-10);
        close(path.beam().power().v, 0.0, 1e-10);
        // The y-directed dipole after the first event cannot emit along y.
        let dark = scattering
            .trace_with_samples(
                cx,
                input,
                0,
                settings(),
                |_| panic!("no boundary before dipole extinction"),
                |flight| match flight {
                    0 => [1.0 - (-0.35_f64).exp(), 0.5, 0.0],
                    1 => [1.0 - (-0.175_f64).exp(), 0.5, 0.25],
                    _ => panic!("extinguished dipole cannot launch another flight"),
                },
            )
            .unwrap();
        assert_eq!(dark.path().termination(), PathTermination::Extinguished);
        assert_eq!(dark.volume_events(), 2);
        assert_eq!(dark.path().surface_events(), 0);
        close(dark.path().traveled_m(), 0.75, 2e-9);
        close(dark.path().beam().power().i, 0.0, 1e-12);
    });
}

#[test]
fn g1_scattering_inside_glass_returns_through_native_refraction_with_circular_handedness() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(1, 2, 0.5), plane(9, 2, -0.5)]).unwrap();
        let air = glass(1.0, 0.0);
        let lens = glass(1.5, 0.2);
        let media = [
            SceneMedium::Homogeneous(&air),
            SceneMedium::Homogeneous(&lens),
        ];
        let scene = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![
                SurfaceBinding {
                    object_id: 1,
                    interaction: Interaction::Dielectric {
                        inside: 0,
                        outside: 1,
                        surface: DielectricSurface::SMOOTH,
                    },
                },
                detector(9),
            ],
            100_000,
        )
        .unwrap();
        let spectrum = MassScatteringSpectrum::constant(1.0).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let volumes = [
            None,
            Some(CurvedRayleighMedium::try_new(&spectrum, 1.0, &phase).unwrap()),
        ];
        let scattering = CurvedRayleighScene::try_new(cx, &scene, &volumes).unwrap();
        let input = scene
            .launch_beam(
                Point3::new(0.0, 0.0, 0.0),
                frame(),
                550.0,
                0,
                Stokes::try_new(1.0, 0.0, 0.0, 1.0).unwrap(),
            )
            .unwrap();
        let mut surface_draws = 0;
        let output = scattering
            .trace_with_samples(
                cx,
                input,
                0,
                settings(),
                |event| {
                    assert_eq!(event, surface_draws);
                    surface_draws += 1;
                    [0.3, 0.7, 0.999999]
                },
                |flight| match flight {
                    0 => [1.0 - (-0.4_f64).exp(), 0.0, 0.0],
                    1 => [1.0 - 1e-12, 0.5, 0.0],
                    _ => panic!("clear air must consume no flight samples"),
                },
            )
            .unwrap();
        let path = output.path();
        assert_eq!(surface_draws, 2);
        assert_eq!(path.surface_events(), 3);
        assert_eq!(output.volume_events(), 1);
        assert_eq!(output.free_flight_draws(), 2);
        assert_eq!(path.termination(), PathTermination::Detected);
        assert_eq!(path.medium_index(), 0);
        assert_eq!(path.last_contact().unwrap().object_id, 9);
        close(path.beam().ray().origin.z, -0.5, 3e-7);
        close(path.traveled_m(), 2.3, 5e-7);
        close(path.optical_path_length_m(), 2.7, 6e-7);
        close(path.optical_depth(), 0.16, 3e-7);
        // At normal incidence conditional Fresnel/PDF transmission has weight
        // one. Glass absorption applies twice; backscatter reverses V once.
        let expected = (-0.16_f64).exp();
        close(path.beam().power().i, expected, 3e-7);
        close(path.beam().power().q, 0.0, 1e-10);
        close(path.beam().power().u, 0.0, 1e-10);
        close(path.beam().power().v, -expected, 3e-7);
        assert_eq!(
            scattering.trace(cx, input, 0, settings(), 31, 42),
            scattering.trace(cx, input, 0, settings(), 31, 42)
        );
        // A collision inside the admitted numerical departure collar cannot
        // safely choose which side of the dielectric supplied the scatterer.
        assert!(matches!(
            scattering.trace_with_samples(
                cx,
                input,
                0,
                settings(),
                |_| [0.3, 0.7, 0.999999],
                |_| [1e-12, 0.5, 0.0]
            ),
            Err(CurvedSceneError::Propagation(
                GradientIndexError::BoundaryResolutionLimit
            ))
        ));
    });
}

#[test]
fn g3_zero_scattering_is_bitwise_native_curved_transport_without_new_draws() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(1, 2, 1.0), plane(9, 2, 2.0)]).unwrap();
        let field = LinearIndex;
        let air = glass(1.0, 0.0);
        let media = [
            SceneMedium::GradientIndex {
                field: &field,
                identity: ContentHash([72; 32]),
            },
            SceneMedium::Homogeneous(&air),
        ];
        let scene = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![
                SurfaceBinding {
                    object_id: 1,
                    interaction: Interaction::Dielectric {
                        inside: 0,
                        outside: 1,
                        surface: DielectricSurface::SMOOTH,
                    },
                },
                detector(9),
            ],
            100_000,
        )
        .unwrap();
        let zero = MassScatteringSpectrum::constant(0.0).unwrap();
        let positive = MassScatteringSpectrum::constant(1.0).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let input = scene
            .launch_beam(
                Point3::new(0.0, 0.0, 0.0),
                frame(),
                550.0,
                0,
                Stokes::try_new(1.0, 0.2, -0.3, 0.4).unwrap(),
            )
            .unwrap();
        for volume in [
            None,
            Some(CurvedRayleighMedium::try_new(&zero, 1.0, &phase).unwrap()),
            Some(CurvedRayleighMedium::try_new(&positive, 0.0, &phase).unwrap()),
        ] {
            let volumes = [volume, None];
            let scattering = CurvedRayleighScene::try_new(cx, &scene, &volumes).unwrap();
            for seed in [1, 31, 97] {
                let output = scattering
                    .trace(cx, input, 0, settings(), seed, 42)
                    .unwrap();
                assert_eq!(
                    *output.path(),
                    scene
                        .trace(cx, input, 0, settings().surfaces, seed, 42)
                        .unwrap()
                );
                assert_eq!(output.volume_events(), 0);
                assert_eq!(output.free_flight_draws(), 0);
            }
            let output = scattering
                .trace_with_samples(
                    cx,
                    input,
                    0,
                    settings(),
                    |_| [0.3, 0.7, 0.999999],
                    |_| panic!("zero rate consumes no sample"),
                )
                .unwrap();
            let native = scene
                .trace_with_samples(cx, input, 0, settings().surfaces, |_| [0.3, 0.7, 0.999999])
                .unwrap();
            assert_eq!(*output.path(), native);
            assert_eq!(output.path().termination(), PathTermination::Detected);
        }
    });
}

#[test]
fn g0_invalid_limits_samples_and_cancellation_never_become_dark_measurements() {
    with_cx(|cx, gate| {
        let geometry = InstanceScene::try_new(vec![plane(9, 2, 2.0)]).unwrap();
        let material = glass(1.33, 0.1);
        let media = [SceneMedium::Homogeneous(&material)];
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![detector(9)], 100_000)
                .unwrap();
        let spectrum = MassScatteringSpectrum::constant(1.0).unwrap();
        let phase = PhaseFunction::Rayleigh;
        assert!(CurvedRayleighMedium::try_new(&spectrum, f64::NAN, &phase).is_err());
        assert!(CurvedRayleighMedium::try_new(&spectrum, -1.0, &phase).is_err());
        assert!(CurvedRayleighMedium::try_new(&spectrum, 1.0, &PhaseFunction::Isotropic).is_err());
        let volumes = [Some(
            CurvedRayleighMedium::try_new(&spectrum, 1.0, &phase).unwrap(),
        )];
        assert!(CurvedRayleighScene::try_new(cx, &scene, &[]).is_err());
        let scattering = CurvedRayleighScene::try_new(cx, &scene, &volumes).unwrap();
        let input = scene
            .launch_beam(
                Point3::new(0.0, 0.0, 0.0),
                frame(),
                550.0,
                0,
                Stokes::unpolarized(1.0),
            )
            .unwrap();
        let mut invalid = settings();
        invalid.max_volume_events = 0;
        assert!(matches!(
            scattering.trace_with_samples(
                cx,
                input,
                0,
                invalid,
                |_| panic!("invalid limits consume no draw"),
                |_| panic!("invalid limits consume no draw")
            ),
            Err(CurvedSceneError::Scene(SceneError::InvalidSettings))
        ));
        for value in [f64::NAN, -0.1, 1.0] {
            assert!(
                scattering
                    .trace_with_samples(
                        cx,
                        input,
                        0,
                        settings(),
                        |_| [0.5; 3],
                        |_| [value, 0.5, 0.0]
                    )
                    .is_err()
            );
        }
        // Exactly competing detector and collision events cannot be assigned
        // an arbitrary order inside the native contact-resolution allowance.
        assert!(matches!(
            scattering.trace_with_samples(
                cx,
                input,
                0,
                settings(),
                |_| [0.5; 3],
                |_| [1.0 - (-2.0_f64).exp(), 0.5, 0.0]
            ),
            Err(CurvedSceneError::Propagation(
                GradientIndexError::BoundaryResolutionLimit
            ))
        ));
        gate.request();
        assert!(matches!(
            scattering.trace_with_samples(
                cx,
                input,
                0,
                settings(),
                |_| panic!("cancelled path consumes no draw"),
                |_| panic!("cancelled path consumes no draw")
            ),
            Err(CurvedSceneError::Scene(SceneError::Cancelled))
        ));
    });
}


struct LinearDensity;
impl CurvedDensityField for LinearDensity {
    fn density_kg_m3(
        &self,
        position: Point3,
    ) -> Result<f64, crate::volumes::transport::VolumeTransportError> {
        if !position.x.is_finite() || !position.y.is_finite() || !position.z.is_finite() {
            return Err(crate::volumes::transport::VolumeTransportError::Density { index: 0 });
        }
        Ok((1.0 + position.z).clamp(0.0, 2.0))
    }
}

struct ViolatingDensity;
impl CurvedDensityField for ViolatingDensity {
    fn density_kg_m3(
        &self,
        _: Point3,
    ) -> Result<f64, crate::volumes::transport::VolumeTransportError> {
        Ok(3.0)
    }
}

#[test]
fn g1_delta_tracking_resolves_heterogeneous_density_on_the_actual_path() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(9, 2, 5.0)]).unwrap();
        let material = glass(1.4, 0.0);
        let media = [SceneMedium::Homogeneous(&material)];
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![detector(9)], 100_000)
                .unwrap();
        let spectrum = MassScatteringSpectrum::constant(1.0).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let density = LinearDensity;
        let volumes = [Some(
            CurvedRayleighMedium::try_new_heterogeneous(
                &spectrum,
                &density,
                2.0,
                ContentHash([91; 32]),
                &phase,
            )
            .unwrap(),
        )];
        let scattering = CurvedRayleighScene::try_new(cx, &scene, &volumes).unwrap();
        let input = scene
            .launch_beam(
                Point3::new(0.0, 0.0, 0.0),
                frame(),
                550.0,
                0,
                Stokes::unpolarized(1.0),
            )
            .unwrap();
        let mut limit = settings();
        limit.max_volume_events = 1;
        let output = scattering
            .trace_with_tracking_samples(
                cx,
                input,
                0,
                limit,
                |_| panic!("candidate scattering precedes the detector"),
                |candidate| match candidate {
                    // majorant=2/m: first candidate at z=0.2, rho=1.2 => p=0.6.
                    0 => [1.0 - (-0.4_f64).exp(), 0.5, 0.0, 0.9],
                    // Null collision leaves state unchanged. Next 0.3 m reaches
                    // z=0.5, rho=1.5 => p=0.75 and is accepted.
                    1 => [1.0 - (-0.6_f64).exp(), 0.5, 0.0, 0.5],
                    _ => panic!("accepted event reaches the volume limit"),
                },
            )
            .unwrap();
        assert_eq!(output.path().termination(), PathTermination::VolumeLimit);
        assert_eq!(output.volume_events(), 1);
        assert_eq!(output.free_flight_draws(), 2);
        close(output.path().traveled_m(), 0.5, 2e-12);
        close(output.path().beam().ray().origin.z, 0.5, 2e-12);
        close(output.path().optical_path_length_m(), 0.7, 2e-12);
        close(output.path().beam().power().i, 1.0, 2e-12);
        assert_eq!(
            volumes[0].unwrap().density_identity(),
            Some(ContentHash([91; 32]))
        );
        assert!(volumes[0].unwrap().is_heterogeneous());

        // The legacy three-variate custom sampler cannot silently correlate
        // acceptance and phase sampling in a heterogeneous field.
        assert!(
            scattering
                .trace_with_samples(cx, input, 0, limit, |_| [0.5; 3], |_| [0.5; 3])
                .is_err()
        );
    });
}

#[test]
fn g0_delta_tracking_refuses_density_majorant_violations() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(9, 2, 5.0)]).unwrap();
        let material = glass(1.33, 0.0);
        let media = [SceneMedium::Homogeneous(&material)];
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![detector(9)], 100_000)
                .unwrap();
        let spectrum = MassScatteringSpectrum::constant(1.0).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let density = ViolatingDensity;
        assert!(
            CurvedRayleighMedium::try_new_heterogeneous(
                &spectrum,
                &density,
                0.0,
                ContentHash([92; 32]),
                &phase,
            )
            .is_err()
        );
        assert!(
            CurvedRayleighMedium::try_new_heterogeneous(
                &spectrum,
                &density,
                2.0,
                ContentHash([0; 32]),
                &phase,
            )
            .is_err()
        );
        let volumes = [Some(
            CurvedRayleighMedium::try_new_heterogeneous(
                &spectrum,
                &density,
                2.0,
                ContentHash([92; 32]),
                &phase,
            )
            .unwrap(),
        )];
        let scattering = CurvedRayleighScene::try_new(cx, &scene, &volumes).unwrap();
        let input = scene
            .launch_beam(
                Point3::new(0.0, 0.0, 0.0),
                frame(),
                550.0,
                0,
                Stokes::unpolarized(1.0),
            )
            .unwrap();
        assert!(
            scattering
                .trace(cx, input, 0, settings(), 19, 7)
                .is_err()
        );
    });
}

#[test]
fn g3_native_event_observer_retains_density_and_field_null_tracking() {
    struct Medium;
    impl GradientIndexField for Medium {
        fn sample(&self, _: Point3, _: f64) -> Result<GradientIndexSample, GradientIndexError> {
            GradientIndexSample::try_new(1.4, Vec3::new(0.0, 0.0, 0.0), 0.2, Some(1.4))
        }
    }
    impl CurvedDensityField for Medium {
        fn density_kg_m3(
            &self,
            _: Point3,
        ) -> Result<f64, crate::volumes::transport::VolumeTransportError> {
            Ok(0.25)
        }
    }
    impl CurvedRayleighField for Medium {
        fn source_field_identity(&self) -> ContentHash {
            ContentHash([183; 32])
        }
        fn scattering_majorant_per_m(
            &self,
            _: f64,
        ) -> Result<f64, crate::volumes::transport::VolumeTransportError> {
            Ok(1.0)
        }
        fn sample_scattering(
            &self,
            _: Point3,
            _: f64,
        ) -> Result<CurvedRayleighSample, crate::volumes::transport::VolumeTransportError> {
            Ok(CurvedRayleighSample {
                scattering_per_m: 0.25,
                depolarization_ratio: 0.0,
            })
        }
    }
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(9, 0, 1.0)]).unwrap();
        let field = Medium;
        let media = [SceneMedium::GradientIndex {
            field: &field,
            identity: ContentHash([183; 32]),
        }];
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![detector(9)], 100_000)
                .unwrap();
        let spectrum = MassScatteringSpectrum::constant(1.0).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let density = CurvedRayleighMedium::try_new_heterogeneous(
            &spectrum,
            &field,
            1.0,
            ContentHash([184; 32]),
            &phase,
        )
        .unwrap();
        let input = scene
            .launch_beam(
                Point3::new(0.0, 0.0, 0.0),
                frame(),
                550.0,
                0,
                Stokes::unpolarized(1.0),
            )
            .unwrap();
        let draws = |flight| match flight {
            0 => [1.0 - (-0.2_f64).exp(), 0.5, 0.0, 0.8], // null at z=0.2
            1 => [1.0 - (-0.3_f64).exp(), 0.5, 0.0, 0.1], // real at z=0.5
            2 => [1.0 - 1e-12, 0.5, 0.0, 0.0],
            _ => panic!("detector interrupts the third flight"),
        };
        for law in [
            CurvedRayleighBinding::Density(density),
            CurvedRayleighBinding::Field(&field),
        ] {
            let volumes = [Some(law)];
            let scattering = CurvedRayleighScene::try_new_with_fields(cx, &scene, &volumes)
                .unwrap();
            let expected = scattering
                .trace_with_tracking_samples(
                    cx,
                    input,
                    0,
                    settings(),
                    |_| panic!("detector consumes no surface variate"),
                    draws,
                )
                .unwrap();
            let mut events = Vec::new();
            let mut arriving = None;
            let actual = scattering
                .trace_with_tracking_samples_events_observed(
                    cx,
                    input,
                    0,
                    settings(),
                    |_| panic!("observer must not resample a surface"),
                    draws,
                    &[],
                    |event| {
                        match event {
                            CurvedRayleighEvent::VolumeArrival { path, phase } => {
                                events.push(0);
                                close(path.beam().ray().origin.z, 0.5, 2e-9);
                                close(path.beam().ray().origin.x, 0.0, 2e-9);
                                close(
                                    phase.eval_cos(0.0).unwrap(),
                                    3.0 / (16.0 * core::f64::consts::PI),
                                    1e-14,
                                );
                                arriving = Some(*path.beam());
                            }
                            CurvedRayleighEvent::VolumeDeparture {
                                path,
                                incoming_beam,
                                phase,
                                pdf,
                            } => {
                                events.push(1);
                                assert_eq!(Some(incoming_beam), arriving);
                                let native = phase.sample_polarized(incoming_beam, [0.5, 0.0])
                                    .unwrap();
                                assert_eq!(*path.beam(), native.beam);
                                assert_eq!(pdf, native.pdf);
                                close(pdf, 3.0 / (16.0 * core::f64::consts::PI), 1e-14);
                                close(path.beam().ray().dir.x, 1.0, 1e-14);
                                close(path.beam().power().i, (-0.1_f64).exp(), 2e-9);
                                close(path.beam().power().q, (-0.1_f64).exp(), 2e-9);
                            }
                            CurvedRayleighEvent::Surface(CurvedTraceEvent::Arrival(path)) => {
                                events.push(2);
                                assert_eq!(path.last_contact().unwrap().object_id, 9);
                            }
                            CurvedRayleighEvent::Surface(CurvedTraceEvent::Departure { .. }) => {
                                panic!("detector has no continuation")
                            }
                        }
                        Ok(())
                    },
                )
                .unwrap();
            assert_eq!(events, [0, 1, 2]);
            assert_eq!(actual, expected);
            assert_eq!(actual.volume_events(), 1);
            assert_eq!(actual.free_flight_draws(), 3);
            close(actual.path().traveled_m(), 1.5, 3e-7);
            close(actual.path().beam().power().i, (-0.3_f64).exp(), 3e-7);
        }
    });
}

#[test]
fn g3_native_event_observer_retains_refractive_surface_proposals() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(1, 2, 0.5), plane(9, 2, 1.0)]).unwrap();
        let first = glass(1.0, 0.0);
        let second = glass(1.5, 0.2);
        let media = [
            SceneMedium::Homogeneous(&first),
            SceneMedium::Homogeneous(&second),
        ];
        let scene = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![
                SurfaceBinding {
                    object_id: 1,
                    interaction: Interaction::Dielectric {
                        inside: 0,
                        outside: 1,
                        surface: DielectricSurface::SMOOTH,
                    },
                },
                detector(9),
            ],
            100_000,
        )
        .unwrap();
        let volumes = [None, None];
        let scattering = CurvedRayleighScene::try_new(cx, &scene, &volumes).unwrap();
        let input = scene
            .launch_beam(
                Point3::new(0.0, 0.0, 0.0),
                frame(),
                550.0,
                0,
                Stokes::unpolarized(1.0),
            )
            .unwrap();
        let expected = scattering
            .trace_with_tracking_samples(
                cx,
                input,
                0,
                settings(),
                |_| [0.3, 0.7, 0.999999],
                |_| panic!("clear media consume no flight draw"),
            )
            .unwrap();
        let mut events = Vec::new();
        let mut incoming = None;
        let actual = scattering
            .trace_with_tracking_samples_events_observed(
                cx,
                input,
                0,
                settings(),
                |event| {
                    assert_eq!(event, 0);
                    [0.3, 0.7, 0.999999]
                },
                |_| panic!("event observation must not introduce flight draws"),
                &[],
                |event| {
                    match event {
                        CurvedRayleighEvent::Surface(CurvedTraceEvent::Arrival(path)) => {
                            events.push(path.last_contact().unwrap().object_id);
                            if path.surface_events() == 1 {
                                incoming = Some(*path.beam());
                                assert_eq!(path.medium_index(), 0);
                            }
                        }
                        CurvedRayleighEvent::Surface(CurvedTraceEvent::Departure {
                            path,
                            incoming_beam,
                            event,
                            pdf,
                            delta,
                        }) => {
                            events.push(2);
                            assert_eq!(Some(incoming_beam), incoming);
                            assert_eq!(path.medium_index(), 1);
                            assert_eq!(event, crate::dielectric::DielectricEvent::Transmission);
                            assert!(delta);
                            close(pdf, 0.96, 1e-14);
                        }
                        _ => panic!("clear media produce no molecular event"),
                    }
                    Ok(())
                },
            )
            .unwrap();
        assert_eq!(events, [1, 2, 9]);
        assert_eq!(actual, expected);
        assert_eq!(actual.volume_events(), 0);
        assert_eq!(actual.free_flight_draws(), 0);
    });
}
