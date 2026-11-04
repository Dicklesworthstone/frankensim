//! Multi-emitter fixtures through the existing SI direct and spectral paths.
//! Analytic controls concern the represented point-source model. Large-rig
//! cases exercise integration, not measured throughput or physical error bounds.

#![cfg(feature = "photometry")]

use core::f64::consts::{FRAC_1_SQRT_2, PI};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use fs_render::camera::{Aperture, CameraProjection, PhysicalCamera};
use fs_render::charts::{Ray, TriMesh};
use fs_render::instances::{GeometryInstance, RigidTransform, SharedGeometry};
use fs_render::photometry::area::{
    EmissionSidedness, GoboTransmission, PhotometricAreaEmitter, PhotometricEmitters,
    ThinLensProjector,
};
use fs_render::photometry::fixture::{
    EmitterId, FixtureAssemblyLimits, FixtureEmitter, FixtureEmitterIdentity, FixtureError,
    FixtureGeometry, FixtureId, FixtureOpticalState, FixtureSurfaceEmitter, PhotometricFixture,
    PhotometricRig,
};
use fs_render::photometry::image::LuminanceJob;
use fs_render::photometry::medium::{HomogeneousMedium, PhaseFunction, ScatterSettings};
use fs_render::photometry::path::{
    IlluminanceSampling, PathIlluminanceJob, PathLuminanceEstimate, PathLuminanceJob, PathMedium,
};
use fs_render::photometry::source::{
    AngularDistribution, PhotometricSource, SPECTRAL_BINS, SourceSpectrum, wavelength_nm,
};
use fs_render::photometry::{
    IlluminanceJob, IlluminanceProbe, PhotometricScene, PhotometryError, TransportOptions,
};
use fs_render::spectral::LiftedSpectrum;
use fs_render::tracer::{
    Camera, DirectStrategy, Material, Primitive, Sampler, Scene, Settings, Shape,
};
use fs_render::volumes::VolumeGrid;
use fs_render::volumes::transport::scattering::{MassScatteringSpectrum, ScatteringGrid};
use fs_render::volumes::transport::{AbsorbingGrid, AbsorbingGridBudget, MassAbsorptionSpectrum};

#[path = "photometric_fixtures/lens_controls.rs"]
mod lens_controls;
#[path = "photometric_fixtures/optical_components.rs"]
mod optical_components;

fn with_cx<R>(run: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 173,
                kernel_id: 0x5048_4649,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        run(&cx, &gate)
    })
}

fn placement(position: [f64; 3]) -> RigidTransform {
    RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], position).unwrap()
}

fn isotropic(position: [f64; 3], candela: f64, dimmer: f64) -> PhotometricSource {
    PhotometricSource::try_new(
        placement(position),
        4.0 * PI * candela,
        dimmer,
        AngularDistribution::isotropic(),
        SourceSpectrum::equal_energy(),
    )
    .unwrap()
}

fn fixture(
    id: u64,
    transform: RigidTransform,
    sources: Vec<PhotometricSource>,
) -> PhotometricFixture {
    PhotometricFixture {
        id: FixtureId(id),
        transform,
        optics: FixtureOpticalState::default(),
        emitters: sources
            .into_iter()
            .enumerate()
            .map(|(i, source)| FixtureEmitter {
                id: EmitterId(i as u64 + 1),
                source,
            })
            .collect(),
        geometry: Vec::new(),
        surfaces: Vec::new(),
    }
}

fn assemble(fixtures: &[PhotometricFixture], cx: &Cx<'_>) -> PhotometricRig {
    PhotometricRig::try_new(fixtures, FixtureAssemblyLimits::default(), cx).unwrap()
}

fn probe(position: [f64; 3], normal: [f64; 3]) -> IlluminanceProbe {
    IlluminanceProbe::try_new(
        Point3::new(position[0], position[1], position[2]),
        Vec3::new(normal[0], normal[1], normal[2]),
    )
    .unwrap()
}

fn close(actual: f64, expected: f64) {
    assert!(actual.is_finite());
    assert!(
        (actual - expected).abs() <= 2.0e-12 * expected.abs().max(1.0),
        "observed={actual:.17e}, expected={expected:.17e}"
    );
}

fn gray() -> Material {
    Material::Lambertian {
        reflectance: LiftedSpectrum { c: [0.0; 3] },
    }
}

fn square(half_width: f64) -> TriMesh {
    TriMesh::new(
        vec![
            [-half_width, -half_width, 0.0],
            [half_width, -half_width, 0.0],
            [half_width, half_width, 0.0],
            [-half_width, half_width, 0.0],
        ],
        vec![[0, 1, 2], [0, 2, 3]],
    )
}

fn native_scene(primitives: Vec<Primitive>) -> Scene {
    Scene {
        primitives,
        lights: Vec::new(),
        environment: None,
        camera: Camera {
            eye: Point3::new(0.0, 0.0, 2.0),
            forward: Vec3::new(0.0, 0.0, -1.0),
            up: Vec3::new(0.0, 1.0, 0.0),
            half_tan: 0.0,
        },
    }
}

fn settings(spp: u32) -> Settings {
    Settings {
        width: 1,
        height: 1,
        spp,
        max_depth: 3,
        sampler: Sampler::Iid,
        strategy: DirectStrategy::Mis,
        seed: 173,
    }
}

fn receiver_scene() -> Scene {
    native_scene(vec![Primitive {
        shape: Shape::Mesh(square(10.0)),
        material: gray(),
        emission: None,
    }])
}

fn direct<'a>(
    rig: &'a PhotometricRig,
    primitives: &'a [Primitive],
    cx: &Cx<'_>,
) -> PhotometricScene<'a> {
    PhotometricScene::try_from_rig(
        primitives,
        rig,
        HomogeneousMedium::vacuum(),
        TransportOptions::default(),
        cx,
    )
    .unwrap()
}

fn paths<'scene, 'field>(
    rig: &'scene PhotometricRig,
    scene: &'scene Scene,
    medium: PathMedium<'scene, 'field>,
    spp: u32,
    cx: &Cx<'_>,
) -> PathLuminanceJob<'scene, 'field> {
    PathLuminanceJob::new_with_rig(
        scene,
        rig,
        &[],
        medium,
        TransportOptions::default(),
        settings(spp),
        cx,
    )
    .unwrap()
}

fn path_bits(value: PathLuminanceEstimate) -> ([u64; 4], u32, Option<u64>) {
    (
        [
            value.xyz[0].to_bits(),
            value.xyz[1].to_bits(),
            value.xyz[2].to_bits(),
            value.luminance_cd_m2.to_bits(),
        ],
        value.samples,
        value.estimated_standard_error_cd_m2.map(f64::to_bits),
    )
}

fn finite_fixture(id: u64, transform: RigidTransform) -> PhotometricFixture {
    let mut item = fixture(id, transform, Vec::new());
    item.optics = FixtureOpticalState::try_new(0.5, [0.4; SPECTRAL_BINS]).unwrap();
    let mesh = SharedGeometry::mesh(square(0.25));
    let mut red = [0.0; SPECTRAL_BINS];
    red[46] = 1.0;
    for (index, (lumens, dimmer, spectrum)) in [
        (100.0 * PI, 0.5, SourceSpectrum::equal_energy()),
        (40.0 * PI, 0.25, SourceSpectrum::try_from_bins(red).unwrap()),
    ]
    .into_iter()
    .enumerate()
    {
        item.geometry.push(FixtureGeometry {
            instance: GeometryInstance::try_new(
                id * 10 + index as u64,
                fs_blake3::hash_domain(
                    "org.frankensim.test.fixture-aperture",
                    b"square-quarter-metre",
                ),
                mesh.clone(),
                placement([index as f64 * 1.5 - 0.75, 0.0, 0.0]),
            )
            .unwrap(),
            material: gray(),
        });
        item.surfaces.push(FixtureSurfaceEmitter {
            id: EmitterId(index as u64 + 11),
            source: PhotometricAreaEmitter::try_new(
                index,
                lumens,
                dimmer,
                spectrum,
                EmissionSidedness::Front,
            )
            .unwrap(),
        });
    }
    item
}

fn world_mesh(instance: &GeometryInstance) -> TriMesh {
    let SharedGeometry::Mesh(mesh) = instance.geometry() else {
        unreachable!()
    };
    let vertices = mesh
        .vertices
        .iter()
        .map(|&[x, y, z]| {
            let point = instance.transform().transform_point(Point3::new(x, y, z));
            [point.x, point.y, point.z]
        })
        .collect();
    TriMesh::new(vertices, mesh.triangles.clone())
}

#[test]
fn transformed_emitters_retain_offsets_axes_calibration_and_distinct_identities() {
    with_cx(|cx, _| {
        let turn =
            RigidTransform::try_new([0.0, FRAC_1_SQRT_2, 0.0, FRAC_1_SQRT_2], [10.0, 2.0, 3.0])
                .unwrap();
        let solid_angle = 2.0 * PI * (1.0 - 0.25_f64.cos());
        let forward = PhotometricSource::try_new(
            placement([0.0, 0.0, 1.0]),
            120.0 * solid_angle,
            0.5,
            AngularDistribution::uniform_cone(0.25).unwrap(),
            SourceSpectrum::equal_energy(),
        )
        .unwrap();
        let backward = PhotometricSource::try_new(
            RigidTransform::try_new([1.0, 0.0, 0.0, 0.0], [0.0, 1.0, -1.0]).unwrap(),
            80.0 * solid_angle,
            0.75,
            AngularDistribution::uniform_cone(0.25).unwrap(),
            SourceSpectrum::equal_energy(),
        )
        .unwrap();
        let mut authored = fixture(41, turn, vec![forward, backward]);
        authored.optics = FixtureOpticalState::try_new(0.5, [0.5; SPECTRAL_BINS]).unwrap();
        let rig = assemble(&[authored.clone()], cx);
        let scene = direct(&rig, &[], cx);
        for (index, (position, direction)) in [
            ([15.0, 2.0, 3.0], [-1.0, 0.0, 0.0]),
            ([5.0, 3.0, 3.0], [1.0, 0.0, 0.0]),
        ]
        .into_iter()
        .enumerate()
        {
            let result = scene
                .illuminance_breakdown(cx, &probe(position, direction))
                .unwrap();
            close(result.total.lux, 15.0 / 16.0);
            close(result.emitters[index].estimate.lux, 15.0 / 16.0);
            assert_eq!(result.emitters[1 - index].estimate.lux, 0.0);
            assert_eq!(
                result.emitters[index].identity,
                Some(FixtureEmitterIdentity {
                    fixture: FixtureId(41),
                    emitter: EmitterId(index as u64 + 1)
                })
            );
            assert_eq!(result.fixtures[0].fixture, FixtureId(41));
            assert_eq!(
                result.total.lux.to_bits(),
                scene
                    .illuminance(cx, &probe(position, direction))
                    .unwrap()
                    .lux
                    .to_bits()
            );
        }
        close(rig.fixtures()[0].flux_lm, 30.0 * solid_angle);
        close(authored.emitters[0].source.flux_lm(), 60.0 * solid_angle);
        assert_eq!(scene.emitter_identities(), Some(rig.emitter_identities()));
        assert_eq!(rig.fixtures()[0].sources, 0..2);
    });
}

#[test]
fn independent_spectra_and_controls_survive_shared_passive_optics() {
    with_cx(|cx, _| {
        let bins = [14, 46];
        let sources = bins.map(|bin| {
            let mut power = [0.0; SPECTRAL_BINS];
            power[bin] = 1.0;
            PhotometricSource::try_new(
                placement([0.0; 3]),
                120.0,
                0.25,
                AngularDistribution::isotropic(),
                SourceSpectrum::try_from_bins(power).unwrap(),
            )
            .unwrap()
        });
        let mut transmission = [0.0; SPECTRAL_BINS];
        transmission[bins[0]] = 0.8;
        transmission[bins[1]] = 0.2;
        let mut authored = fixture(55, placement([0.0, 0.0, 3.0]), sources.to_vec());
        authored.optics = FixtureOpticalState::try_new(0.5, transmission).unwrap();
        let lit = assemble(&[authored.clone()], cx);
        for (index, expected_flux) in [12.0, 3.0].into_iter().enumerate() {
            close(lit.sources()[index].flux_lm(), expected_flux);
            close(
                lit.sources()[index]
                    .spectral_intensity_w_sr_nm(
                        Vec3::new(0.0, 0.0, -1.0),
                        wavelength_nm(bins[index]),
                    )
                    .unwrap(),
                sources[index]
                    .spectral_intensity_w_sr_nm(
                        Vec3::new(0.0, 0.0, -1.0),
                        wavelength_nm(bins[index]),
                    )
                    .unwrap()
                    * 0.5
                    * transmission[bins[index]],
            );
            assert_eq!(
                lit.sources()[index]
                    .spectral_intensity_w_sr_nm(
                        Vec3::new(0.0, 0.0, -1.0),
                        wavelength_nm(bins[1 - index])
                    )
                    .unwrap(),
                0.0
            );
        }
        authored.emitters[0].source = PhotometricSource::try_new(
            placement([0.0; 3]),
            120.0,
            0.0,
            AngularDistribution::isotropic(),
            sources[0].spectrum().clone(),
        )
        .unwrap();
        let dimmed = assemble(&[authored], cx);
        assert_eq!(dimmed.emitter_identities(), lit.emitter_identities());
        assert_eq!(dimmed.sources()[0].flux_lm(), 0.0);
        assert_eq!(
            dimmed.sources()[1].flux_lm().to_bits(),
            lit.sources()[1].flux_lm().to_bits()
        );
        let scene = direct(&dimmed, &[], cx);
        let result = scene
            .illuminance_breakdown(cx, &probe([0.0; 3], [0.0, 0.0, 1.0]))
            .unwrap();
        close(result.total.lux, 3.0 / (4.0 * PI * 9.0));
        assert_eq!(result.emitters.len(), 2);
        assert_eq!(result.emitters[0].estimate.lux, 0.0);
        assert_eq!(
            result.total.numerical.kind,
            fs_evidence::NumericalKind::NoClaim
        );
        assert!(result.total.total_error_bound_lux.is_none());
        let native = receiver_scene();
        let mut lit_image = paths(&lit, &native, PathMedium::Vacuum, 128, cx);
        let mut dimmed_image = paths(&dimmed, &native, PathMedium::Vacuum, 128, cx);
        lit_image.advance(cx, usize::MAX).unwrap();
        dimmed_image.advance(cx, usize::MAX).unwrap();
        let lit_xyz = lit_image.pixel_estimate(0).unwrap().xyz;
        let dimmed_xyz = dimmed_image.pixel_estimate(0).unwrap().xyz;
        assert!(lit_xyz[1] > dimmed_xyz[1] && dimmed_xyz[1] > 0.0);
        assert!(lit_xyz[2] / lit_xyz[1] > dimmed_xyz[2] / dimmed_xyz[1]);
    });
}

#[test]
fn transformed_body_geometry_casts_native_shadows_without_duplicating_emitters() {
    with_cx(|cx, _| {
        let turn =
            RigidTransform::try_new([0.0, FRAC_1_SQRT_2, 0.0, FRAC_1_SQRT_2], [3.0, -4.0, 2.0])
                .unwrap();
        let mut authored = fixture(
            9,
            turn,
            vec![
                isotropic([0.0, 0.0, 1.0], 100.0, 1.0),
                isotropic([2.0, 0.0, 1.0], 100.0, 1.0),
            ],
        );
        let body = SharedGeometry::mesh(square(0.5));
        authored.geometry.push(FixtureGeometry {
            instance: GeometryInstance::try_new(
                90,
                fs_blake3::hash_domain("org.frankensim.test.fixture-body", b"square"),
                body.clone(),
                RigidTransform::identity(),
            )
            .unwrap(),
            material: gray(),
        });
        let rig = assemble(&[authored], cx);
        assert!(body.ptr_eq(rig.geometry()[0].instance.geometry()));
        let mut native = native_scene(Vec::new());
        assert_eq!(
            rig.append_geometry_to(&mut native.primitives, cx).unwrap(),
            0..1
        );
        assert_eq!(native.primitives.len(), 1);
        assert!(matches!(
            rig.append_geometry_to(&mut native.primitives, cx),
            Err(FixtureError::DuplicateGeometryId(90))
        ));
        assert_eq!(native.primitives.len(), 1);
        let point = turn.transform_point(Point3::new(0.0, 0.0, -1.0));
        let receiver =
            IlluminanceProbe::try_new(point, turn.transform_vector(Vec3::new(0.0, 0.0, 1.0)))
                .unwrap();
        let scene = direct(&rig, &native.primitives, cx);
        let result = scene.illuminance_breakdown(cx, &receiver).unwrap();
        assert_eq!(result.emitters[0].estimate.lux, 0.0);
        close(
            result.emitters[1].estimate.lux,
            100.0 / (8.0 * 2.0_f64.sqrt()),
        );
        let receivers = [receiver];
        let sampling = IlluminanceSampling {
            samples_per_probe: 8,
            max_depth: 2,
            sampler: Sampler::Iid,
            seed: 173,
        };
        // A calibrated finite emitter may use the same native body primitive.
        // It remains an independent binding beside the grouped point sources.
        let surfaces = [PhotometricAreaEmitter::try_new(
            0,
            100.0,
            1.0,
            SourceSpectrum::equal_energy(),
            EmissionSidedness::Both,
        )
        .unwrap()];
        let mut paths = PathIlluminanceJob::new_with_rig(
            &native,
            &rig,
            &surfaces,
            PathMedium::Vacuum,
            &receivers,
            TransportOptions::default(),
            sampling,
            cx,
        )
        .unwrap();
        let mut flat = PathIlluminanceJob::new_with_emitters(
            &native,
            PhotometricEmitters {
                points: rig.sources(),
                surfaces: &surfaces,
            },
            PathMedium::Vacuum,
            &receivers,
            TransportOptions::default(),
            sampling,
            cx,
        )
        .unwrap();
        paths.advance(cx, usize::MAX).unwrap();
        flat.advance(cx, usize::MAX).unwrap();
        assert_eq!(paths.emitter_identities(), Some(rig.emitter_identities()));
        close(
            paths.probe_estimate(0).unwrap().direct_lux,
            result.total.lux,
        );
        assert!(paths.probe_estimate(0).unwrap().continuous_lux > 0.0);
        assert_eq!(paths.probe_estimate(0), flat.probe_estimate(0));
    });
}

#[test]
fn physical_camera_and_spectral_paths_preserve_lux_and_luminance_units() {
    with_cx(|cx, _| {
        let mut authored = fixture(
            7,
            placement([0.0, 0.0, 3.0]),
            vec![
                isotropic([0.0; 3], 90.0, 0.5),
                isotropic([0.0, 0.0, 1.0], 160.0, 0.25),
            ],
        );
        authored.optics = FixtureOpticalState::try_new(0.8, [0.5; SPECTRAL_BINS]).unwrap();
        let rig = assemble(&[authored], cx);
        let native = receiver_scene();
        let scene = direct(&rig, &native.primitives, cx);
        close(
            scene
                .illuminance(cx, &probe([0.0; 3], [0.0, 0.0, 1.0]))
                .unwrap()
                .lux,
            3.0,
        );
        let camera = PhysicalCamera::try_look_at(
            Point3::new(0.0, 0.0, 2.0),
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            CameraProjection::try_vertical_fov(PI / 6.0).unwrap(),
            2.0,
            Aperture::try_circular(0.0).unwrap(),
        )
        .unwrap();
        let mut image =
            LuminanceJob::new(&scene, &camera, 1, 1, 8.0, ScatterSettings::default()).unwrap();
        image.advance(cx, 1).unwrap();
        close(image.results()[0].surface_cd_m2, 1.5 / PI);
        assert_eq!(image.results()[0].single_scatter.luminance_cd_m2, 0.0);
        assert!(image.results()[0].total_error_bound_cd_m2.is_none());
        let mut assembled = PathLuminanceJob::new_with_rig_camera(
            &native,
            &camera,
            &rig,
            &[],
            PathMedium::Vacuum,
            TransportOptions::default(),
            settings(32),
            cx,
        )
        .unwrap();
        let explicit = [
            isotropic([0.0, 0.0, 3.0], 18.0, 1.0),
            isotropic([0.0, 0.0, 4.0], 16.0, 1.0),
        ];
        let mut flat = PathLuminanceJob::new_with_camera(
            &native,
            &camera,
            PhotometricEmitters::points_only(&explicit),
            PathMedium::Vacuum,
            TransportOptions::default(),
            settings(32),
            cx,
        )
        .unwrap();
        assembled.advance(cx, usize::MAX).unwrap();
        flat.advance(cx, usize::MAX).unwrap();
        let observed = assembled.pixel_estimate(0).unwrap();
        assert_eq!(
            assembled.emitter_identities(),
            Some(rig.emitter_identities())
        );
        assert!(flat.emitter_identities().is_none());
        assert_eq!(image.emitter_identities(), Some(rig.emitter_identities()));
        let reference = flat.pixel_estimate(0).unwrap();
        for (actual, expected) in observed.xyz.into_iter().zip(reference.xyz) {
            close(actual, expected);
        }
        assert!(observed.luminance_cd_m2 > 0.0);
        assert_eq!(observed.numerical.kind, fs_evidence::NumericalKind::NoClaim);
        assert!(observed.total_error_bound_cd_m2.is_none());
    });
}

#[test]
fn assembled_rig_reuses_single_scattering_and_native_spectral_volume_transport() {
    with_cx(|cx, _| {
        let rig = assemble(
            &[fixture(
                3,
                placement([0.0, 0.0, 2.0]),
                vec![
                    isotropic([-0.3, 0.0, 0.0], 50.0, 1.0),
                    isotropic([0.3, 0.0, 0.0], 75.0, 0.5),
                ],
            )],
            cx,
        );
        let homogeneous = HomogeneousMedium::try_new(0.2, 0.3, PhaseFunction::Isotropic).unwrap();
        let grouped = PhotometricScene::try_from_rig(
            &[],
            &rig,
            homogeneous.clone(),
            TransportOptions::default(),
            cx,
        )
        .unwrap();
        let flat = PhotometricScene::try_new(
            &[],
            rig.sources(),
            homogeneous,
            TransportOptions::default(),
            cx,
        )
        .unwrap();
        let ray = Ray {
            origin: Point3::new(0.0, 0.0, 0.0),
            dir: Vec3::new(0.0, 0.0, 1.0),
        };
        let light = grouped
            .luminance_along_ray(cx, &ray, 1.0, ScatterSettings::default())
            .unwrap();
        let reference = flat
            .luminance_along_ray(cx, &ray, 1.0, ScatterSettings::default())
            .unwrap();
        assert!(light.single_scatter.luminance_cd_m2 > 0.0);
        assert_eq!(
            light.luminance_cd_m2.to_bits(),
            reference.luminance_cd_m2.to_bits()
        );
        assert!(light.total_error_bound_cd_m2.is_none());
        let density = [1.0];
        let grid = ScatteringGrid::try_new(
            cx,
            AbsorbingGrid::try_new(
                cx,
                VolumeGrid::new([1; 3], &density, [-1.0, -1.0, 0.0], [2.0, 2.0, 1.0]),
                MassAbsorptionSpectrum::constant(0.2).unwrap(),
                None,
                AbsorbingGridBudget {
                    max_grid_cells: 16,
                    max_steps_per_segment: 64,
                },
            )
            .unwrap(),
            MassScatteringSpectrum::constant(0.3).unwrap(),
            PhaseFunction::Isotropic,
        )
        .unwrap();
        let mut native = native_scene(Vec::new());
        native.camera.eye = ray.origin;
        native.camera.forward = ray.dir;
        let mut paths = paths(&rig, &native, PathMedium::Scattering(&grid), 64, cx);
        let mut flat_paths = PathLuminanceJob::new(
            &native,
            rig.sources(),
            PathMedium::Scattering(&grid),
            TransportOptions::default(),
            settings(64),
            cx,
        )
        .unwrap();
        paths.advance(cx, usize::MAX).unwrap();
        flat_paths.advance(cx, usize::MAX).unwrap();
        assert!(paths.pixel_estimate(0).unwrap().luminance_cd_m2 > 0.0);
        assert_eq!(
            path_bits(paths.pixel_estimate(0).unwrap()),
            path_bits(flat_paths.pixel_estimate(0).unwrap())
        );
    });
}

#[test]
fn cancellation_retains_completed_receivers_and_paths_and_resumes_identically() {
    with_cx(|cx, gate| {
        let authored = [fixture(
            81,
            placement([0.0, 0.0, 3.0]),
            vec![
                isotropic([-0.5, 0.0, 0.0], 50.0, 1.0),
                isotropic([0.5, 0.0, 0.0], 30.0, 0.4),
            ],
        )];
        let rig = assemble(&authored, cx);
        let native = receiver_scene();
        let scene = direct(&rig, &native.primitives, cx);
        let receivers = [
            probe([0.0; 3], [0.0, 0.0, 1.0]),
            probe([1.0, 0.0, 0.0], [0.0, 0.0, 1.0]),
        ];
        let mut full_direct = IlluminanceJob::new(&scene, &receivers).unwrap();
        full_direct.advance(cx, usize::MAX).unwrap();
        let mut split_direct = IlluminanceJob::new(&scene, &receivers).unwrap();
        split_direct.advance(cx, 1).unwrap();
        let mut full = paths(&rig, &native, PathMedium::Vacuum, 16, cx);
        full.advance(cx, usize::MAX).unwrap();
        let mut split = paths(&rig, &native, PathMedium::Vacuum, 16, cx);
        split.advance(cx, 3).unwrap();
        let prefix = path_bits(split.pixel_estimate(0).unwrap());
        gate.request();
        assert!(matches!(
            PhotometricRig::try_new(&authored, FixtureAssemblyLimits::default(), cx),
            Err(FixtureError::Cancelled)
        ));
        assert!(matches!(
            split_direct.advance(cx, 1),
            Err(PhotometryError::Cancelled)
        ));
        assert!(matches!(
            split.advance(cx, 8),
            Err(PhotometryError::Cancelled)
        ));
        assert_eq!(split_direct.results(), &full_direct.results()[..1]);
        assert_eq!(path_bits(split.pixel_estimate(0).unwrap()), prefix);
        assert_eq!(split.completed_path_samples(), 3);
        assert_eq!(split.emitter_identities(), Some(rig.emitter_identities()));
        with_cx(|resumed, _| {
            split_direct.advance(resumed, usize::MAX).unwrap();
            while !split.is_complete() {
                split.advance(resumed, 5).unwrap();
            }
        });
        assert_eq!(split_direct.results(), full_direct.results());
        assert_eq!(
            path_bits(split.pixel_estimate(0).unwrap()),
            path_bits(full.pixel_estimate(0).unwrap())
        );
        assert_eq!(
            split_direct.emitter_identities(),
            Some(rig.emitter_identities())
        );
    });
}

#[test]
fn three_hundred_fixtures_with_twenty_four_emitters_preserve_flux_and_direct_lux() {
    with_cx(|cx, _| {
        let mut expected_lux = 0.0;
        let mut expected_lumens = 0.0;
        let fixtures: Vec<_> = (0..300)
            .map(|fixture_index| {
                let x = f64::from(fixture_index % 20) * 0.5 - 4.75;
                let y = f64::from(fixture_index / 20) * 0.5 - 3.5;
                let sources = (0..24)
                    .map(|emitter_index| {
                        let offset = f64::from(emitter_index) * 0.03 - 0.345;
                        let dimmer = [0.0, 0.25, 0.5, 1.0][emitter_index as usize % 4];
                        let intensity = (1.0 + f64::from(emitter_index % 3)) * dimmer * 0.4;
                        let r2 = (x + offset).powi(2) + y * y + 36.0;
                        expected_lux += intensity * 6.0 / (r2 * r2.sqrt());
                        expected_lumens += 4.0 * PI * intensity;
                        isotropic(
                            [offset, 0.0, 0.0],
                            1.0 + f64::from(emitter_index % 3),
                            dimmer,
                        )
                    })
                    .collect();
                let mut item =
                    fixture(fixture_index as u64 + 1000, placement([x, y, 6.0]), sources);
                item.optics = FixtureOpticalState::try_new(0.8, [0.5; SPECTRAL_BINS]).unwrap();
                item
            })
            .collect();
        let rig = assemble(&fixtures, cx);
        let replay = assemble(&fixtures, cx);
        assert_eq!(rig.sources().len(), 7200);
        assert_eq!(rig.fixtures().len(), 300);
        assert_eq!(rig.fixtures()[299].sources, 7176..7200);
        assert_eq!(rig.emitter_identities(), replay.emitter_identities());
        assert_eq!(rig.fixtures(), replay.fixtures());
        close(
            rig.fixtures().iter().map(|item| item.flux_lm).sum(),
            expected_lumens,
        );
        let direct = direct(&rig, &[], cx);
        let receiver = probe([0.0; 3], [0.0, 0.0, 1.0]);
        let result = direct.illuminance_breakdown(cx, &receiver).unwrap();
        close(result.total.lux, expected_lux);
        assert_eq!(result.emitters.len(), 7200);
        assert_eq!(result.fixtures.len(), 300);
        let native = native_scene(Vec::new());
        let receivers = [receiver];
        let mut paths = PathIlluminanceJob::new_with_rig(
            &native,
            &rig,
            &[],
            PathMedium::Vacuum,
            &receivers,
            TransportOptions::default(),
            IlluminanceSampling {
                samples_per_probe: 1,
                max_depth: 1,
                sampler: Sampler::Iid,
                seed: 173,
            },
            cx,
        )
        .unwrap();
        paths.advance(cx, 1).unwrap();
        let estimate = paths.probe_estimate(0).unwrap();
        assert_eq!(estimate.direct_lux.to_bits(), result.total.lux.to_bits());
        assert_eq!(estimate.continuous_lux, 0.0);
        assert!(estimate.total_error_bound_lux.is_none());
    });
}

#[test]
fn finite_fixture_controls_and_frames_reach_visible_camera_and_probe_transport() {
    with_cx(|cx, _| {
        let turn =
            RigidTransform::try_new([0.0, FRAC_1_SQRT_2, 0.0, FRAC_1_SQRT_2], [4.0, -3.0, 2.0])
                .unwrap();
        let authored = finite_fixture(101, turn);
        let rig = assemble(&[authored.clone()], cx);
        assert!(
            rig.sources().is_empty(),
            "finite emitters must not also add point-source power"
        );
        assert_eq!(rig.fixtures()[0].surfaces, 0..2);
        close(rig.fixtures()[0].flux_lm, 12.0 * PI);
        close(rig.surface_sources()[0].flux_lm(), 10.0 * PI);
        close(rig.surface_sources()[1].flux_lm(), 2.0 * PI);
        close(authored.surfaces[0].source.flux_lm(), 50.0 * PI);

        // Existing geometry precedes the rig: its local slots 0/1 must bind to
        // native primitive indices 1/2 for both source hits and light sampling.
        let mut native = receiver_scene();
        assert_eq!(
            rig.append_geometry_to(&mut native.primitives, cx).unwrap(),
            1..3
        );
        let mut explicit = receiver_scene();
        let mut flat_sources = Vec::new();
        for (index, body) in authored.geometry.iter().enumerate() {
            let instance = GeometryInstance::try_new(
                body.instance.object_id(),
                body.instance.geometry_identity(),
                body.instance.geometry().clone(),
                turn.compose(body.instance.transform()).unwrap(),
            )
            .unwrap();
            let primitive_index = explicit.primitives.len();
            explicit.primitives.push(Primitive {
                shape: Shape::Mesh(world_mesh(&instance)),
                material: body.material,
                emission: None,
            });
            let source = &authored.surfaces[index].source;
            flat_sources.push(
                PhotometricAreaEmitter::try_new(
                    primitive_index,
                    source.clear_flux_lm(),
                    0.5,
                    source.spectrum().clone(),
                    source.sidedness(),
                )
                .unwrap()
                .filtered([0.4; SPECTRAL_BINS])
                .unwrap(),
            );
            assert_eq!(rig.surface_sources()[index].primitive_index(), index);
            assert_eq!(
                rig.surface_emitter_identities()[index],
                FixtureEmitterIdentity {
                    fixture: FixtureId(101),
                    emitter: EmitterId(index as u64 + 11),
                }
            );
        }

        let policy = Settings {
            max_depth: 1,
            sampler: Sampler::OwenSobol,
            ..settings(8192)
        };
        let mut colors = Vec::new();
        for body in &authored.geometry {
            let frame = turn.compose(body.instance.transform()).unwrap();
            let centre = frame.transform_point(Point3::new(0.0, 0.0, 0.0));
            let camera = PhysicalCamera::try_look_at(
                centre.offset(frame.transform_vector(Vec3::new(0.0, 0.0, 2.0))),
                centre,
                frame.transform_vector(Vec3::new(0.0, 1.0, 0.0)),
                CameraProjection::try_vertical_fov(0.04).unwrap(),
                2.0,
                Aperture::try_circular(0.02).unwrap(),
            )
            .unwrap();
            let mut grouped = PathLuminanceJob::new_with_rig_camera(
                &native,
                &camera,
                &rig,
                &[],
                PathMedium::Vacuum,
                TransportOptions::default(),
                policy,
                cx,
            )
            .unwrap();
            let mut flat = PathLuminanceJob::new_with_camera(
                &explicit,
                &camera,
                PhotometricEmitters {
                    points: &[],
                    surfaces: &flat_sources,
                },
                PathMedium::Vacuum,
                TransportOptions::default(),
                policy,
                cx,
            )
            .unwrap();
            grouped.advance(cx, usize::MAX).unwrap();
            flat.advance(cx, usize::MAX).unwrap();
            let value = grouped.pixel_estimate(0).unwrap();
            for (actual, expected) in value
                .xyz
                .into_iter()
                .zip(flat.pixel_estimate(0).unwrap().xyz)
            {
                close(actual, expected);
            }
            assert_eq!(
                grouped.surface_emitter_identities(),
                Some(rig.surface_emitter_identities())
            );
            assert!(value.total_error_bound_cd_m2.is_none());
            colors.push(value.xyz);
        }
        // Native finite-aperture luminance uses actual area: 10*pi lm over
        // 0.25 m² emits 40 cd/m². This is a bounded sampling tolerance.
        assert!((colors[0][1] - 40.0).abs() < 0.002 * 40.0);
        assert!(colors[1][1] > 0.0);
        assert!(colors[1][0] / colors[1][1] > colors[0][0] / colors[0][1]);

        let probes = [IlluminanceProbe::try_new(
            turn.transform_point(Point3::new(0.0, 0.0, 2.0)),
            turn.transform_vector(Vec3::new(0.0, 0.0, -1.0)),
        )
        .unwrap()];
        let sampling = IlluminanceSampling {
            samples_per_probe: 512,
            max_depth: 1,
            sampler: Sampler::OwenSobol,
            seed: 173,
        };
        let mut grouped = PathIlluminanceJob::new_with_rig(
            &native,
            &rig,
            &[],
            PathMedium::Vacuum,
            &probes,
            TransportOptions::default(),
            sampling,
            cx,
        )
        .unwrap();
        let mut flat = PathIlluminanceJob::new_with_emitters(
            &explicit,
            PhotometricEmitters {
                points: &[],
                surfaces: &flat_sources,
            },
            PathMedium::Vacuum,
            &probes,
            TransportOptions::default(),
            sampling,
            cx,
        )
        .unwrap();
        grouped.advance(cx, usize::MAX).unwrap();
        flat.advance(cx, usize::MAX).unwrap();
        let result = grouped.probe_estimate(0).unwrap();
        assert_eq!(result.direct_lux, 0.0);
        assert!(result.continuous_lux > 0.0);
        close(
            result.illuminance_lux,
            flat.probe_estimate(0).unwrap().illuminance_lux,
        );
        assert_eq!(
            grouped.surface_emitter_identities(),
            Some(rig.surface_emitter_identities())
        );
    });
}

#[test]
fn fixture_projector_preserves_local_optical_frame_clear_flux_and_mask_loss() {
    with_cx(|cx, _| {
        let fixture_frame =
            RigidTransform::try_new([0.0, FRAC_1_SQRT_2, 0.0, FRAC_1_SQRT_2], [5.0, -2.0, 4.0])
                .unwrap();
        let body_frame =
            RigidTransform::try_new([0.1_f64.sin(), 0.0, 0.0, 0.1_f64.cos()], [0.4, 0.0, 0.2])
                .unwrap();
        let optical_frame = RigidTransform::try_new(
            [0.0, 0.0, 0.2_f64.sin(), 0.2_f64.cos()],
            [0.02, -0.01, 0.04],
        )
        .unwrap();
        let mut pupil = square(0.05);
        let vertices = pupil
            .vertices
            .iter()
            .map(|&[x, y, z]| {
                let point = optical_frame.transform_point(Point3::new(x, y, z));
                [point.x, point.y, point.z]
            })
            .collect();
        // Construct a fresh native BVH after transforming the aperture vertices.
        pupil = TriMesh::new(vertices, pupil.triangles.clone());
        let projector = ThinLensProjector::try_new(
            optical_frame,
            0.1,
            0.1 * 3.0 / 3.1,
            GoboTransmission::try_new(0.04, 0.02, 2, 1, vec![1.0, 0.25]).unwrap(),
        )
        .unwrap();
        let mut authored = fixture(202, fixture_frame, Vec::new());
        authored.optics = FixtureOpticalState::try_new(0.5, [0.4; SPECTRAL_BINS]).unwrap();
        authored.geometry.push(FixtureGeometry {
            instance: GeometryInstance::try_new(
                2021,
                fs_blake3::hash_domain("org.frankensim.test.fixture-projector", b"offset-pupil"),
                SharedGeometry::mesh(pupil),
                body_frame,
            )
            .unwrap(),
            material: gray(),
        });
        authored.surfaces.push(FixtureSurfaceEmitter {
            id: EmitterId(21),
            source: PhotometricAreaEmitter::try_new(
                0,
                240.0,
                0.5,
                SourceSpectrum::equal_energy(),
                EmissionSidedness::Front,
            )
            .unwrap()
            .with_projector(projector)
            .unwrap(),
        });
        let rig = assemble(&[authored], cx);
        let source = &rig.surface_sources()[0];
        close(source.clear_flux_lm(), 24.0);
        close(source.flux_lm(), 15.0);
        close(rig.fixtures()[0].flux_lm, 15.0);
        assert_eq!(source.projector().unwrap().frame(), optical_frame);
        let mut scene = receiver_scene();
        assert_eq!(
            rig.append_geometry_to(&mut scene.primitives, cx).unwrap(),
            1..2
        );
        let world_optical = fixture_frame
            .compose(body_frame)
            .unwrap()
            .compose(optical_frame)
            .unwrap();
        let policy = Settings {
            max_depth: 1,
            sampler: Sampler::OwenSobol,
            ..settings(8192)
        };
        for (slope, transmission) in [(0.06_f64, 1.0), (-0.06, 0.25)] {
            let camera = PhysicalCamera::try_look_at(
                world_optical.transform_point(Point3::new(3.0 * slope, 0.0, 3.0)),
                world_optical.transform_point(Point3::new(0.0, 0.0, 0.0)),
                world_optical.transform_vector(Vec3::new(0.0, 1.0, 0.0)),
                CameraProjection::try_vertical_fov(1.0e-6).unwrap(),
                3.0,
                Aperture::try_circular(0.0).unwrap(),
            )
            .unwrap();
            let mut job = PathLuminanceJob::new_with_rig_camera(
                &scene,
                &camera,
                &rig,
                &[],
                PathMedium::Vacuum,
                TransportOptions::default(),
                policy,
                cx,
            )
            .unwrap();
            job.advance(cx, usize::MAX).unwrap();
            // Independent native-pupil law: clear_flux/A * s²/G * M/cos⁴.
            let expected =
                24.0 / 0.01 * 0.01 / (0.04 * 0.02) * transmission * (1.0 + slope * slope).powi(2);
            let value = job.pixel_estimate(0).unwrap();
            assert!((value.luminance_cd_m2 - expected).abs() < 0.002 * expected);
            assert!(value.total_error_bound_cd_m2.is_none());
        }
    });
}

#[test]
fn finite_fixture_mount_refusals_reordering_and_cancel_resume_preserve_results() {
    with_cx(|cx, gate| {
        let rig = assemble(&[finite_fixture(303, placement([0.0, 0.0, 1.0]))], cx);
        let mut native = receiver_scene();
        assert!(matches!(
            PathLuminanceJob::new_with_rig(
                &native,
                &rig,
                &[],
                PathMedium::Vacuum,
                TransportOptions::default(),
                settings(16),
                cx,
            ),
            Err(PhotometryError::Fixture(
                FixtureError::SurfaceNotMounted { .. }
            ))
        ));
        rig.append_geometry_to(&mut native.primitives, cx).unwrap();
        assert!(matches!(
            PhotometricScene::try_from_rig(
                &native.primitives,
                &rig,
                HomogeneousMedium::vacuum(),
                TransportOptions::default(),
                cx,
            ),
            Err(PhotometryError::UnsupportedEmitters { .. })
        ));

        let original = rig.geometry()[0].instance.clone();
        native.primitives[1].shape = Shape::Instance(
            GeometryInstance::try_new(
                original.object_id(),
                original.geometry_identity(),
                original.geometry().clone(),
                placement([0.0, 0.0, 0.1])
                    .compose(original.transform())
                    .unwrap(),
            )
            .unwrap(),
        );
        assert!(matches!(
            PathLuminanceJob::new_with_rig(
                &native,
                &rig,
                &[],
                PathMedium::Vacuum,
                TransportOptions::default(),
                settings(16),
                cx,
            ),
            Err(PhotometryError::Fixture(
                FixtureError::SurfaceGeometryMismatch { .. }
            ))
        ));
        native.primitives[1].shape = Shape::Instance(original.clone());
        native.primitives.push(Primitive {
            shape: Shape::Instance(original),
            material: gray(),
            emission: None,
        });
        assert!(matches!(
            PathLuminanceJob::new_with_rig(
                &native,
                &rig,
                &[],
                PathMedium::Vacuum,
                TransportOptions::default(),
                settings(16),
                cx,
            ),
            Err(PhotometryError::Fixture(FixtureError::DuplicateGeometryId(
                _
            )))
        ));
        native.primitives.pop().unwrap();

        let duplicate = [PhotometricAreaEmitter::try_new(
            1,
            1.0,
            1.0,
            SourceSpectrum::equal_energy(),
            EmissionSidedness::Front,
        )
        .unwrap()];
        assert!(matches!(
            PathLuminanceJob::new_with_rig(
                &native,
                &rig,
                &duplicate,
                PathMedium::Vacuum,
                TransportOptions::default(),
                settings(16),
                cx,
            ),
            Err(PhotometryError::InvalidOptions)
        ));
        native.primitives.swap(1, 2);
        native.camera.eye = Point3::new(-0.75, 0.0, 3.0);
        let mut full = paths(&rig, &native, PathMedium::Vacuum, 16, cx);
        full.advance(cx, usize::MAX).unwrap();
        assert!(full.pixel_estimate(0).unwrap().luminance_cd_m2 > 0.0);
        let mut split = paths(&rig, &native, PathMedium::Vacuum, 16, cx);
        split.advance(cx, 3).unwrap();
        let accepted = path_bits(split.pixel_estimate(0).unwrap());
        gate.request();
        assert!(matches!(
            split.advance(cx, 4),
            Err(PhotometryError::Cancelled)
        ));
        assert_eq!(path_bits(split.pixel_estimate(0).unwrap()), accepted);
        assert_eq!(
            split.surface_emitter_identities(),
            Some(rig.surface_emitter_identities())
        );
        with_cx(|resumed, _| {
            while !split.is_complete() {
                split.advance(resumed, 5).unwrap();
            }
        });
        assert_eq!(
            path_bits(split.pixel_estimate(0).unwrap()),
            path_bits(full.pixel_estimate(0).unwrap())
        );
    });
}
