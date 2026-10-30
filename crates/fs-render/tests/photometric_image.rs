//! G1/G3/G4/G5 checks for physical surface luminance and single-scattered haze.
//! Analytic references use metres, candela, lux and cd/m² without exposure or
//! a display transform. Refinement discrepancies are not physical error bounds.

#![cfg(feature = "photometry")]

use core::f64::consts::PI;

use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use fs_render::camera::{Aperture, CameraProjection, PhysicalCamera};
use fs_render::charts::{Ray, TriMesh};
use fs_render::instances::RigidTransform;
use fs_render::photometry::image::{LuminanceEstimate, LuminanceJob};
use fs_render::photometry::medium::{HomogeneousMedium, PhaseFunction, ScatterSettings};
use fs_render::photometry::source::{
    AngularDistribution, PhotometricSource, SPECTRAL_BINS, SourceSpectrum, wavelength_nm,
};
use fs_render::photometry::{PhotometricScene, PhotometryError, TransportOptions};
use fs_render::spectral::LiftedSpectrum;
use fs_render::tracer::{Material, Primitive, Shape};

fn with_cx<R>(operation: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 7,
                kernel_id: 0x5048_494d,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        operation(&cx, &gate)
    })
}

fn source(position: [f64; 3], candela: f64, spectrum: SourceSpectrum) -> PhotometricSource {
    PhotometricSource::try_new(
        RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], position).unwrap(),
        4.0 * PI * candela,
        1.0,
        AngularDistribution::isotropic(),
        spectrum,
    )
    .unwrap()
}

fn floor(z: f64, reflectance: LiftedSpectrum, facing_up: bool) -> Primitive {
    Primitive {
        shape: Shape::Mesh(TriMesh::new(
            vec![
                [-10.0, -10.0, z],
                [10.0, -10.0, z],
                [10.0, 10.0, z],
                [-10.0, 10.0, z],
            ],
            if facing_up {
                vec![[0, 1, 2], [0, 2, 3]]
            } else {
                vec![[0, 2, 1], [0, 3, 2]]
            },
        )),
        material: Material::Lambertian { reflectance },
        emission: None,
    }
}

fn camera(height: f64) -> PhysicalCamera {
    PhysicalCamera::try_look_at(
        Point3::new(0.0, 0.0, height),
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        CameraProjection::try_vertical_fov(PI / 3.0).unwrap(),
        height,
        Aperture::try_circular(0.0).unwrap(),
    )
    .unwrap()
}

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(actual.is_finite());
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1.0e-12),
        "observed={actual:.17e}, expected={expected:.17e}"
    );
}

fn forward_ray() -> Ray {
    Ray {
        origin: Point3::new(0.0, 0.0, 0.0),
        dir: Vec3::new(0.0, 0.0, 1.0),
    }
}

fn result_bits(result: &LuminanceEstimate) -> [u64; 3] {
    [
        result.luminance_cd_m2.to_bits(),
        result.surface_cd_m2.to_bits(),
        result.single_scatter.luminance_cd_m2.to_bits(),
    ]
}

#[test]
fn g1_center_pixel_reports_lambertian_rho_e_over_pi_in_cd_per_square_metre() {
    with_cx(|cx, _| {
        let surfaces = [floor(0.0, LiftedSpectrum { c: [0.0; 3] }, true)];
        let sources = [source(
            [0.0, 0.0, 3.0],
            90.0,
            SourceSpectrum::equal_energy(),
        )];
        let scene = PhotometricScene::try_new(
            &surfaces,
            &sources,
            HomogeneousMedium::vacuum(),
            TransportOptions::default(),
            cx,
        )
        .unwrap();
        let camera = camera(2.0);
        let mut image =
            LuminanceJob::new(&scene, &camera, 1, 1, 8.0, ScatterSettings::default()).unwrap();
        assert_eq!(image.advance(cx, 1).unwrap(), 1);
        assert!(image.is_complete());
        let pixel = &image.results()[0];
        // 90 cd at 3 m gives 10 lx; a rho=0.5 Lambertian returns 5/pi cd/m².
        close(pixel.luminance_cd_m2, 5.0 / PI, 2.0e-12);
        close(pixel.surface_cd_m2, 5.0 / PI, 2.0e-12);
        assert_eq!(pixel.single_scatter.luminance_cd_m2, 0.0);
        assert!(pixel.total_error_bound_cd_m2.is_none());
    });
}

#[test]
fn g1_surface_reflection_uses_emitter_spectrum_and_existing_material_spectrum() {
    with_cx(|cx, _| {
        let surfaces = [floor(
            0.0,
            LiftedSpectrum {
                c: [0.0, 2.0, -0.5],
            },
            true,
        )];
        let ray = Ray {
            origin: Point3::new(0.0, 0.0, 2.0),
            dir: Vec3::new(0.0, 0.0, -1.0),
        };
        let mut values = Vec::new();
        for bin in [12, 60] {
            let mut spectrum = [0.0; SPECTRAL_BINS];
            spectrum[bin] = 1.0;
            let sources = [source(
                [0.0, 0.0, 3.0],
                90.0,
                SourceSpectrum::try_from_bins(spectrum).unwrap(),
            )];
            let scene = PhotometricScene::try_new(
                &surfaces,
                &sources,
                HomogeneousMedium::vacuum(),
                TransportOptions::default(),
                cx,
            )
            .unwrap();
            let observed = scene
                .luminance_along_ray(cx, &ray, 8.0, ScatterSettings::default())
                .unwrap();
            let normalized_lambda = (wavelength_nm(bin) - 380.0) / 400.0;
            let polynomial = 2.0 * normalized_lambda - 0.5;
            let reflectance = 0.5 + 0.5 * polynomial / (1.0 + polynomial * polynomial).sqrt();
            close(observed.surface_cd_m2, 10.0 * reflectance / PI, 2.0e-12);
            values.push(observed.surface_cd_m2);
        }
        assert!(
            values[1] > 1.8 * values[0],
            "equal lux must retain spectral reflectance contrast"
        );
    });
}

#[test]
fn g1_absorption_attenuates_the_light_and_camera_legs_once_each() {
    with_cx(|cx, _| {
        let surfaces = [floor(0.0, LiftedSpectrum { c: [0.0; 3] }, true)];
        let sources = [source(
            [0.0, 0.0, 3.0],
            90.0,
            SourceSpectrum::equal_energy(),
        )];
        let medium = HomogeneousMedium::try_new(0.2, 0.0, PhaseFunction::Isotropic).unwrap();
        let scene =
            PhotometricScene::try_new(&surfaces, &sources, medium, TransportOptions::default(), cx)
                .unwrap();
        for height in [2.0, 4.0] {
            let ray = Ray {
                origin: Point3::new(0.0, 0.0, height),
                dir: Vec3::new(0.0, 0.0, -1.0),
            };
            let result = scene
                .luminance_along_ray(cx, &ray, 10.0, ScatterSettings::default())
                .unwrap();
            let expected = (5.0 / PI) * (-0.2_f64 * (3.0 + height)).exp();
            close(result.luminance_cd_m2, expected, 2.0e-12);
            assert_eq!(result.single_scatter.luminance_cd_m2, 0.0);
        }
    });
}

#[test]
fn g3_opaque_mesh_shadows_the_fixture_to_haze_connection() {
    with_cx(|cx, _| {
        // The x=1 wall hides a source at x=2 from the entire x=0 camera ray,
        // but does not itself intersect that ray.
        let blockers = [Primitive {
            shape: Shape::Mesh(TriMesh::new(
                vec![
                    [1.0, -10.0, -10.0],
                    [1.0, 10.0, -10.0],
                    [1.0, 10.0, 10.0],
                    [1.0, -10.0, 10.0],
                ],
                vec![[0, 1, 2], [0, 2, 3]],
            )),
            material: Material::Lambertian {
                reflectance: LiftedSpectrum { c: [0.0; 3] },
            },
            emission: None,
        }];
        let sources = [source(
            [2.0, 0.0, 2.0],
            100.0,
            SourceSpectrum::equal_energy(),
        )];
        let medium = HomogeneousMedium::try_new(0.05, 0.1, PhaseFunction::Isotropic).unwrap();
        let clear =
            PhotometricScene::try_new(&[], &sources, medium, TransportOptions::default(), cx)
                .unwrap();
        let shadow =
            PhotometricScene::try_new(&blockers, &sources, medium, TransportOptions::default(), cx)
                .unwrap();
        let visible = clear
            .luminance_along_ray(cx, &forward_ray(), 4.0, ScatterSettings::default())
            .unwrap();
        let hidden = shadow
            .luminance_along_ray(cx, &forward_ray(), 4.0, ScatterSettings::default())
            .unwrap();
        assert!(visible.single_scatter.luminance_cd_m2 > 0.1);
        assert_eq!(hidden.luminance_cd_m2, 0.0);
        assert_eq!(hidden.single_scatter.luminance_cd_m2, 0.0);
        assert!(
            hidden
                .single_scatter
                .omitted_multiple_scattering_bound_cd_m2
                .is_none()
        );
        assert!(hidden.total_error_bound_cd_m2.is_none());
    });
}

#[test]
fn g1_henyey_greenstein_uses_propagation_directions_for_forward_scattering() {
    with_cx(|cx, _| {
        let medium =
            HomogeneousMedium::try_new(0.05, 0.1, PhaseFunction::HenyeyGreenstein { g: 0.5 })
                .unwrap();
        let settings = ScatterSettings {
            initial_panels: 1,
            max_panels: 1,
            absolute_tolerance_cd_m2: 0.0,
            relative_tolerance: 0.0,
        };
        let mut observed = Vec::new();
        // The one quadrature point is z=2; both fixtures are exactly 4 m away.
        // From z=6 the incoming and outgoing propagation vectors are both -Z.
        for z in [6.0, -2.0] {
            let sources = [source([0.0, 0.0, z], 100.0, SourceSpectrum::equal_energy())];
            let scene =
                PhotometricScene::try_new(&[], &sources, medium, TransportOptions::default(), cx)
                    .unwrap();
            let result = scene
                .luminance_along_ray(cx, &forward_ray(), 4.0, settings)
                .unwrap();
            assert_eq!(result.single_scatter.panels, 1);
            assert!(!result.single_scatter.tolerance_met);
            assert!(
                result
                    .single_scatter
                    .estimated_quadrature_error_cd_m2
                    .is_none()
            );
            observed.push(result.single_scatter.luminance_cd_m2);
        }
        // HG p(+1)/p(-1) = ((1+g)/(1-g))^3 = 27 for g=0.5.
        close(observed[0] / observed[1], 27.0, 3.0e-13);
        let forward_phase = 1.5 / (4.0 * PI * 0.5 * 0.5);
        let expected = (100.0 / 16.0)
            * (-0.15_f64 * 4.0).exp()
            * forward_phase
            * (0.1 / 0.15)
            * -(-0.15_f64 * 4.0).exp_m1();
        close(observed[0], expected, 2.0e-12);
    });
}

#[test]
fn g3_first_opaque_surface_truncates_haze_at_the_real_hit_distance() {
    with_cx(|cx, _| {
        let surfaces = [floor(2.0, LiftedSpectrum { c: [0.0; 3] }, false)];
        let sources = [source(
            [2.0, 0.0, 1.0],
            100.0,
            SourceSpectrum::equal_energy(),
        )];
        let medium = HomogeneousMedium::try_new(0.02, 0.08, PhaseFunction::Isotropic).unwrap();
        let clear =
            PhotometricScene::try_new(&[], &sources, medium, TransportOptions::default(), cx)
                .unwrap();
        let wall =
            PhotometricScene::try_new(&surfaces, &sources, medium, TransportOptions::default(), cx)
                .unwrap();
        let settings = ScatterSettings::default();
        let reference = clear
            .luminance_along_ray(cx, &forward_ray(), 2.0, settings)
            .unwrap();
        let observed = wall
            .luminance_along_ray(cx, &forward_ray(), 4.0, settings)
            .unwrap();
        let untruncated = clear
            .luminance_along_ray(cx, &forward_ray(), 4.0, settings)
            .unwrap();
        assert_eq!(
            observed.single_scatter.luminance_cd_m2.to_bits(),
            reference.single_scatter.luminance_cd_m2.to_bits()
        );
        assert!(
            untruncated.single_scatter.luminance_cd_m2 > observed.single_scatter.luminance_cd_m2
        );
        assert!(observed.surface_cd_m2 > 0.0);
        close(
            observed.luminance_cd_m2,
            observed.surface_cd_m2 + observed.single_scatter.luminance_cd_m2,
            2.0e-14,
        );
    });
}

#[test]
fn g3_haze_refuses_unsampled_point_source_singularities_on_the_camera_segment() {
    with_cx(|cx, _| {
        let settings = ScatterSettings {
            initial_panels: 1,
            max_panels: 2,
            absolute_tolerance_cd_m2: 1.0e9,
            relative_tolerance: 1.0,
        };
        let medium = HomogeneousMedium::try_new(0.0, 0.1, PhaseFunction::Isotropic).unwrap();
        // None of the midpoint panels samples z=0.3. Even a permissive
        // quadrature tolerance must not turn its divergent integral finite.
        for x in [0.0, 0.5e-6] {
            let sources = [source([x, 0.0, 0.3], 100.0, SourceSpectrum::equal_energy())];
            let scene =
                PhotometricScene::try_new(&[], &sources, medium, TransportOptions::default(), cx)
                    .unwrap();
            assert!(matches!(
                scene.luminance_along_ray(cx, &forward_ray(), 1.0, settings),
                Err(PhotometryError::SourceTooClose { source: 0 })
            ));
        }
        // No scattered-light singularity exists in vacuum, and an opaque
        // surface terminates the visible segment before a hidden source.
        let sources = [source(
            [0.0, 0.0, 0.3],
            100.0,
            SourceSpectrum::equal_energy(),
        )];
        let vacuum = PhotometricScene::try_new(
            &[],
            &sources,
            HomogeneousMedium::vacuum(),
            TransportOptions::default(),
            cx,
        )
        .unwrap();
        assert_eq!(
            vacuum
                .luminance_along_ray(cx, &forward_ray(), 1.0, settings)
                .unwrap()
                .luminance_cd_m2,
            0.0
        );
        let surfaces = [floor(0.1, LiftedSpectrum { c: [0.0; 3] }, false)];
        let shadow =
            PhotometricScene::try_new(&surfaces, &sources, medium, TransportOptions::default(), cx)
                .unwrap();
        assert_eq!(
            shadow
                .luminance_along_ray(cx, &forward_ray(), 1.0, settings)
                .unwrap()
                .luminance_cd_m2,
            0.0
        );
    });
}

#[test]
fn g4_g5_progressive_luminance_cancel_and_resume_preserve_pixel_bits() {
    with_cx(|cx, gate| {
        let surfaces = [floor(0.0, LiftedSpectrum { c: [0.0; 3] }, true)];
        let sources = [source(
            [0.0, 0.0, 3.0],
            90.0,
            SourceSpectrum::equal_energy(),
        )];
        let scene = PhotometricScene::try_new(
            &surfaces,
            &sources,
            HomogeneousMedium::vacuum(),
            TransportOptions::default(),
            cx,
        )
        .unwrap();
        let camera = camera(2.0);
        let settings = ScatterSettings::default();
        let mut full = LuminanceJob::new(&scene, &camera, 3, 2, 10.0, settings).unwrap();
        assert_eq!(full.advance(cx, 6).unwrap(), 6);
        let reference: Vec<_> = full.results().iter().map(result_bits).collect();
        let mut job = LuminanceJob::new(&scene, &camera, 3, 2, 10.0, settings).unwrap();
        assert_eq!(job.advance(cx, 2).unwrap(), 2);
        let prefix: Vec<_> = job.results().iter().map(result_bits).collect();
        gate.request();
        assert!(matches!(
            job.advance(cx, 3),
            Err(PhotometryError::Cancelled)
        ));
        assert_eq!(
            job.results().iter().map(result_bits).collect::<Vec<_>>(),
            prefix
        );
        with_cx(|resumed, _| {
            while !job.is_complete() {
                job.advance(resumed, 1).unwrap();
            }
        });
        assert_eq!(
            job.results().iter().map(result_bits).collect::<Vec<_>>(),
            reference
        );
    });
}
