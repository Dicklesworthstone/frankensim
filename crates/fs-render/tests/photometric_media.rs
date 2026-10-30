//! Calibrated measurements through the original physical density/spectrum grid.
#![cfg(feature = "photometry")]

use core::f64::consts::PI;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use fs_render::camera::{Aperture, CameraProjection, PhysicalCamera};
use fs_render::charts::{Ray, TriMesh};
use fs_render::instances::RigidTransform;
use fs_render::photometry::image::LuminanceJob;
use fs_render::photometry::medium::{HomogeneousMedium, PhaseFunction, ScatterSettings};
use fs_render::photometry::source::{
    AngularDistribution, PhotometricSource, SPECTRAL_BINS, SourceSpectrum, photopic_weight_nm,
    wavelength_nm,
};
use fs_render::photometry::transport::PhotometricMedium;
use fs_render::photometry::{
    IlluminanceJob, IlluminanceProbe, PhotometricScene, PhotometryError, TransportOptions,
};
use fs_render::spectral::LiftedSpectrum;
use fs_render::tracer::{Material, Primitive, Shape};
use fs_render::volumes::VolumeGrid;
use fs_render::volumes::transport::scattering::{
    MassScatteringSpectrum, PhaseFunction as PhysicalPhase, ScatteringGrid,
};
use fs_render::volumes::transport::{
    AbsorbingGrid, AbsorbingGridBudget, MassAbsorptionSpectrum, VolumeTransportError,
};

fn with_cx<R>(operation: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 3,
                kernel_id: 0x5048_4d45,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        operation(&cx, &gate)
    })
}

fn fixture(position: [f64; 3], spectrum: SourceSpectrum) -> PhotometricSource {
    PhotometricSource::try_new(
        RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], position).unwrap(),
        400.0 * PI,
        1.0,
        AngularDistribution::isotropic(),
        spectrum,
    )
    .unwrap()
}

fn receiver(z: f64) -> IlluminanceProbe {
    IlluminanceProbe::try_new(Point3::new(0.0, 0.0, z), Vec3::new(0.0, 0.0, 1.0)).unwrap()
}

fn floor(reflectance: LiftedSpectrum) -> Primitive {
    Primitive {
        shape: Shape::Mesh(TriMesh::new(
            vec![
                [-10.0, -10.0, 0.0],
                [10.0, -10.0, 0.0],
                [10.0, 10.0, 0.0],
                [-10.0, 10.0, 0.0],
            ],
            vec![[0, 1, 2], [0, 2, 3]],
        )),
        material: Material::Lambertian { reflectance },
        emission: None,
    }
}

fn camera_ray() -> Ray {
    Ray {
        origin: Point3::new(0.0, 0.0, 3.0),
        dir: Vec3::new(0.0, 0.0, -1.0),
    }
}

fn budget(steps: u32) -> AbsorbingGridBudget {
    AbsorbingGridBudget {
        max_grid_cells: 64,
        max_steps_per_segment: steps,
    }
}

fn close(observed: f64, expected: f64) {
    assert!(observed.is_finite());
    assert!(
        (observed - expected).abs() <= 3.0e-12 * expected.abs().max(1.0),
        "observed={observed:.17e} expected={expected:.17e}"
    );
}

#[test]
fn g1_two_density_layers_match_column_mass_extinction_in_lux() {
    with_cx(|cx, _| {
        let density = [1.0, 3.0];
        let grid = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1, 1, 2], &density, [-1.0, -1.0, 1.0], [2.0, 2.0, 0.5]),
            MassAbsorptionSpectrum::constant(0.4).unwrap(),
            None,
            budget(16),
        )
        .unwrap();
        let sources = [fixture([0.0, 0.0, 3.0], SourceSpectrum::equal_energy())];
        let scene = PhotometricScene::try_with_medium(
            &[],
            &sources,
            PhotometricMedium::Absorbing(&grid),
            TransportOptions::default(),
            cx,
        )
        .unwrap();
        // Two half-metre cells: column mass = 1*0.5 + 3*0.5 = 2 kg/m².
        close(
            scene.illuminance(cx, &receiver(0.0)).unwrap().lux,
            (100.0 / 9.0) * (-0.8_f64).exp(),
        );
        assert!(
            scene
                .illuminance(cx, &receiver(0.0))
                .unwrap()
                .total_error_bound_lux
                .is_none()
        );
    });
}

#[test]
fn g1_oblique_connection_integrates_only_the_bounded_grid_chord() {
    with_cx(|cx, _| {
        let density = [1.0];
        let grid = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1; 3], &density, [1.0, -1.0, 1.0], [1.0, 2.0, 1.0]),
            MassAbsorptionSpectrum::constant(0.5).unwrap(),
            None,
            budget(16),
        )
        .unwrap();
        let sources = [fixture([3.0, 0.0, 3.0], SourceSpectrum::equal_energy())];
        let scene = PhotometricScene::try_with_medium(
            &[],
            &sources,
            PhotometricMedium::Absorbing(&grid),
            TransportOptions::default(),
            cx,
        )
        .unwrap();
        // The diagonal segment is 3*sqrt(2) long; only sqrt(2) lies in the box.
        close(
            scene.illuminance(cx, &receiver(0.0)).unwrap().lux,
            (100.0 / 18.0) / 2.0_f64.sqrt() * (-0.5 * 2.0_f64.sqrt()).exp(),
        );
        let offset =
            IlluminanceProbe::try_new(Point3::new(3.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0))
                .unwrap();
        close(scene.illuminance(cx, &offset).unwrap().lux, 100.0 / 9.0);
        // Rotate fixture and bounded box together by 90 degrees around Z.
        // This permutation is exactly representable by the native grid axes.
        let rotated_grid = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1; 3], &density, [-1.0, 1.0, 1.0], [2.0, 1.0, 1.0]),
            MassAbsorptionSpectrum::constant(0.5).unwrap(),
            None,
            budget(16),
        )
        .unwrap();
        let rotated_sources = [fixture([0.0, 3.0, 3.0], SourceSpectrum::equal_energy())];
        let rotated = PhotometricScene::try_with_medium(
            &[],
            &rotated_sources,
            PhotometricMedium::Absorbing(&rotated_grid),
            TransportOptions::default(),
            cx,
        )
        .unwrap();
        assert_eq!(
            scene.illuminance(cx, &receiver(0.0)).unwrap().lux.to_bits(),
            rotated
                .illuminance(cx, &receiver(0.0))
                .unwrap()
                .lux
                .to_bits()
        );
    });
}

#[test]
fn g1_spectral_grid_filters_equal_lux_sources_before_photopic_integration() {
    with_cx(|cx, _| {
        let density = [2.0];
        let grid = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1; 3], &density, [-1.0, -1.0, 1.0], [2.0, 2.0, 1.0]),
            MassAbsorptionSpectrum::table(&[(380.0, 0.1), (780.0, 1.1)]).unwrap(),
            None,
            budget(16),
        )
        .unwrap();
        let mut measured = Vec::new();
        for bin in [12, 60] {
            let mut bins = [0.0; SPECTRAL_BINS];
            bins[bin] = 1.0;
            let sources = [fixture(
                [0.0, 0.0, 3.0],
                SourceSpectrum::try_from_bins(bins).unwrap(),
            )];
            let scene = PhotometricScene::try_with_medium(
                &[],
                &sources,
                PhotometricMedium::Absorbing(&grid),
                TransportOptions::default(),
                cx,
            )
            .unwrap();
            let kappa = 0.1 + (wavelength_nm(bin) - 380.0) / 400.0;
            let observed = scene.illuminance(cx, &receiver(0.0)).unwrap().lux;
            close(observed, (100.0 / 9.0) * (-2.0 * kappa).exp());
            measured.push(observed);
        }
        assert!(measured[0] > 3.0 * measured[1]);
    });
}

#[test]
fn g1_surface_response_retains_spectral_attenuation_on_both_legs() {
    with_cx(|cx, _| {
        let density = [2.0];
        let grid = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1; 3], &density, [-1.0, -1.0, 1.0], [2.0, 2.0, 1.0]),
            MassAbsorptionSpectrum::table(&[(380.0, 0.1), (780.0, 1.1)]).unwrap(),
            None,
            budget(16),
        )
        .unwrap();
        let spectrum = SourceSpectrum::equal_energy();
        let sources = [fixture([0.0, 0.0, 3.0], spectrum.clone())];
        let surfaces = [floor(LiftedSpectrum {
            c: [0.0, 2.0, -0.5],
        })];
        let scene = PhotometricScene::try_with_medium(
            &surfaces,
            &sources,
            PhotometricMedium::Absorbing(&grid),
            TransportOptions::default(),
            cx,
        )
        .unwrap();
        let response: f64 = (0..SPECTRAL_BINS)
            .map(|bin| {
                let x = (wavelength_nm(bin) - 380.0) / 400.0;
                let polynomial = 2.0 * x - 0.5;
                let rho = 0.5 + 0.5 * polynomial / (1.0 + polynomial * polynomial).sqrt();
                spectrum.bins_w_per_lumen_nm()[bin]
                    * photopic_weight_nm(bin)
                    * rho
                    * (-4.0 * (0.1 + x)).exp()
            })
            .sum();
        let expected = (100.0 / 9.0) * 683.0 * response / PI;
        let component = scene
            .direct_surface_luminance_along_ray(cx, &camera_ray(), 5.0)
            .unwrap();
        let image = scene
            .luminance_along_ray(cx, &camera_ray(), 5.0, ScatterSettings::default())
            .unwrap();
        close(component.surface_cd_m2, expected);
        assert_eq!(
            component.surface_cd_m2.to_bits(),
            image.luminance_cd_m2.to_bits()
        );
        assert_eq!(image.single_scatter.luminance_cd_m2, 0.0);
        assert!(image.total_error_bound_cd_m2.is_none());
    });
}

#[test]
fn g5_vacuum_and_gray_grids_preserve_homogeneous_fast_path_bits() {
    with_cx(|cx, _| {
        let density = [1.0];
        let sources = [fixture([0.0, 0.0, 3.0], SourceSpectrum::equal_energy())];
        let surfaces = [floor(LiftedSpectrum { c: [0.0; 3] })];
        for absorption in [0.0, 0.2] {
            let grid = AbsorbingGrid::try_new(
                cx,
                VolumeGrid::new([1; 3], &density, [-1.0, -1.0, 0.0], [2.0, 2.0, 3.0]),
                MassAbsorptionSpectrum::constant(absorption).unwrap(),
                None,
                budget(16),
            )
            .unwrap();
            let physical = PhotometricScene::try_with_medium(
                &surfaces,
                &sources,
                PhotometricMedium::Absorbing(&grid),
                TransportOptions::default(),
                cx,
            )
            .unwrap();
            let homogeneous = PhotometricScene::try_new(
                &surfaces,
                &sources,
                HomogeneousMedium::try_new(absorption, 0.0, PhaseFunction::Isotropic).unwrap(),
                TransportOptions::default(),
                cx,
            )
            .unwrap();
            assert_eq!(
                physical
                    .illuminance(cx, &receiver(0.0))
                    .unwrap()
                    .lux
                    .to_bits(),
                homogeneous
                    .illuminance(cx, &receiver(0.0))
                    .unwrap()
                    .lux
                    .to_bits()
            );
            let settings = ScatterSettings::default();
            assert_eq!(
                physical
                    .luminance_along_ray(cx, &camera_ray(), 5.0, settings)
                    .unwrap()
                    .luminance_cd_m2
                    .to_bits(),
                homogeneous
                    .luminance_along_ray(cx, &camera_ray(), 5.0, settings)
                    .unwrap()
                    .luminance_cd_m2
                    .to_bits()
            );
        }
    });
}

#[test]
fn g1_scattering_grid_attenuates_direct_components_and_refuses_missing_haze() {
    with_cx(|cx, _| {
        let density = [2.0];
        let absorption = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1; 3], &density, [-1.0, -1.0, 1.0], [2.0, 2.0, 1.0]),
            MassAbsorptionSpectrum::constant(0.2).unwrap(),
            None,
            budget(16),
        )
        .unwrap();
        let grid = ScatteringGrid::try_new(
            cx,
            absorption,
            MassScatteringSpectrum::constant(0.3).unwrap(),
            PhysicalPhase::henyey_greenstein(0.4).unwrap(),
        )
        .unwrap();
        let sources = [fixture([0.0, 0.0, 3.0], SourceSpectrum::equal_energy())];
        let surfaces = [floor(LiftedSpectrum { c: [0.0; 3] })];
        let scene = PhotometricScene::try_with_medium(
            &surfaces,
            &sources,
            PhotometricMedium::Scattering(&grid),
            TransportOptions::default(),
            cx,
        )
        .unwrap();
        close(
            scene.illuminance(cx, &receiver(0.0)).unwrap().lux,
            (100.0 / 9.0) * (-1.0_f64).exp(),
        );
        close(
            scene
                .direct_surface_luminance_along_ray(cx, &camera_ray(), 5.0)
                .unwrap()
                .surface_cd_m2,
            (100.0 / 9.0) * 0.5 / PI * (-2.0_f64).exp(),
        );
        assert!(matches!(
            scene.luminance_along_ray(cx, &camera_ray(), 5.0, ScatterSettings::default()),
            Err(PhotometryError::UnsupportedMedium { .. })
        ));
        let camera = PhysicalCamera::try_look_at(
            camera_ray().origin,
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            CameraProjection::try_vertical_fov(PI / 3.0).unwrap(),
            3.0,
            Aperture::try_circular(0.0).unwrap(),
        )
        .unwrap();
        assert!(matches!(
            LuminanceJob::new(&scene, &camera, 1, 1, 5.0, ScatterSettings::default()),
            Err(PhotometryError::UnsupportedMedium { .. })
        ));
    });
}

#[test]
fn g3_thermal_grid_emission_is_not_silently_dropped_from_images() {
    with_cx(|cx, _| {
        let density = [1.0];
        let temperature = [1800.0];
        let grid = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1; 3], &density, [-1.0, -1.0, 1.0], [2.0, 2.0, 1.0]),
            MassAbsorptionSpectrum::constant(0.2).unwrap(),
            Some(&temperature),
            budget(16),
        )
        .unwrap();
        let sources = [fixture([0.0, 0.0, 3.0], SourceSpectrum::equal_energy())];
        let surfaces = [floor(LiftedSpectrum { c: [0.0; 3] })];
        let scene = PhotometricScene::try_with_medium(
            &surfaces,
            &sources,
            PhotometricMedium::Absorbing(&grid),
            TransportOptions::default(),
            cx,
        )
        .unwrap();
        close(
            scene.illuminance(cx, &receiver(0.0)).unwrap().lux,
            (100.0 / 9.0) * (-0.2_f64).exp(),
        );
        close(
            scene
                .direct_surface_luminance_along_ray(cx, &camera_ray(), 5.0)
                .unwrap()
                .surface_cd_m2,
            (100.0 / 9.0) * 0.5 / PI * (-0.4_f64).exp(),
        );
        assert!(matches!(
            scene.luminance_along_ray(cx, &camera_ray(), 5.0, ScatterSettings::default()),
            Err(PhotometryError::UnsupportedMedium { .. })
        ));
    });
}

#[test]
fn g4_grid_budget_refusal_retains_the_completed_probe_prefix() {
    with_cx(|cx, _| {
        let density = [1.0, 2.0];
        let grid = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1, 1, 2], &density, [-1.0, -1.0, 0.0], [2.0, 2.0, 1.0]),
            MassAbsorptionSpectrum::constant(0.2).unwrap(),
            None,
            budget(1),
        )
        .unwrap();
        let sources = [fixture([0.0, 0.0, 3.0], SourceSpectrum::equal_energy())];
        let scene = PhotometricScene::try_with_medium(
            &[],
            &sources,
            PhotometricMedium::Absorbing(&grid),
            TransportOptions::default(),
            cx,
        )
        .unwrap();
        let probes = [receiver(2.5), receiver(0.0)];
        let mut job = IlluminanceJob::new(&scene, &probes).unwrap();
        assert!(matches!(
            job.advance(cx, 2),
            Err(PhotometryError::PhysicalMedium(
                VolumeTransportError::Budget
            ))
        ));
        assert_eq!(job.results().len(), 1);
        close(job.results()[0].lux, 400.0);
    });
}

#[test]
fn g4_g5_physical_grid_jobs_cancel_and_resume_without_changing_bits() {
    with_cx(|cx, gate| {
        let density = [1.0, 2.0];
        let grid = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1, 1, 2], &density, [-1.0, -1.0, 0.0], [2.0, 2.0, 1.0]),
            MassAbsorptionSpectrum::table(&[(380.0, 0.1), (780.0, 0.7)]).unwrap(),
            None,
            budget(16),
        )
        .unwrap();
        let sources = [fixture([0.0, 0.0, 3.0], SourceSpectrum::equal_energy())];
        let scene = PhotometricScene::try_with_medium(
            &[],
            &sources,
            PhotometricMedium::Absorbing(&grid),
            TransportOptions::default(),
            cx,
        )
        .unwrap();
        let probes = [receiver(0.0), receiver(0.5), receiver(1.5)];
        let reference: Vec<_> = probes
            .iter()
            .map(|probe| scene.illuminance(cx, probe).unwrap().lux.to_bits())
            .collect();
        let mut job = IlluminanceJob::new(&scene, &probes).unwrap();
        assert_eq!(job.advance(cx, 1).unwrap(), 1);
        gate.request();
        assert!(matches!(
            job.advance(cx, 2),
            Err(PhotometryError::Cancelled)
        ));
        assert_eq!(job.results().len(), 1);
        with_cx(|fresh, _| {
            while !job.is_complete() {
                job.advance(fresh, 1).unwrap();
            }
        });
        assert_eq!(
            job.results()
                .iter()
                .map(|value| value.lux.to_bits())
                .collect::<Vec<_>>(),
            reference
        );
    });
}
