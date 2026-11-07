//! Aspheric geometry through the original spectral, Stokes and group-time path.
#![cfg(feature = "tracer")]

use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use fs_render::camera::{
    Aperture, CameraError, CameraProjection, LensSample, LensSurface, PhysicalCamera,
    SequentialLens,
};
use fs_render::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricGlass, GlassProvenance};

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 53,
                kernel_id: 16,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx);
    });
}

fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}

fn camera(aspheric: bool, absorption: f64) -> PhysicalCamera {
    let glass = DielectricGlass::new(
        CauchyIor::try_new(1.46, 0.004, 0.0).unwrap(),
        BeerLambertAbsorption::try_constant(absorption).unwrap(),
        GlassProvenance::Custom,
    );
    let surfaces = if aspheric {
        [
            LensSurface::try_new_aspheric(
                0.08,
                20.0,
                -1.0,
                [-1e4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                0.012,
                Some(glass),
            )
            .unwrap(),
            LensSurface::try_new_aspheric(
                0.09,
                -20.0,
                -0.4,
                [2e3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                0.012,
                None,
            )
            .unwrap(),
        ]
    } else {
        [
            LensSurface::try_new(0.08, 20.0, 0.012, Some(glass)).unwrap(),
            LensSurface::try_new(0.09, -20.0, 0.012, None).unwrap(),
        ]
    };
    PhysicalCamera::try_look_at(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.0, 0.0, 1.0),
        Vec3::new(0.0, 1.0, 0.0),
        CameraProjection::try_focal_sensor(0.05, 0.024).unwrap(),
        1.0,
        Aperture::try_circular(0.001).unwrap(),
    )
    .unwrap()
    .with_sequential_lens(SequentialLens::try_new(550.0, &surfaces).unwrap())
    .unwrap()
}

#[test]
fn zero_aspheric_departure_preserves_original_surfaces_and_spectral_beams() {
    let original = camera(false, 0.2);
    let lens = original.sequential_lens().unwrap();
    let surfaces: Vec<_> = lens
        .surfaces()
        .iter()
        .map(|s| {
            let zero = LensSurface::try_new_aspheric(
                s.vertex_z_m(),
                s.curvature_per_m(),
                0.0,
                [0.0; 7],
                s.clear_radius_m(),
                s.medium_after().copied(),
            )
            .unwrap();
            assert_eq!(*s, zero);
            assert_eq!(zero.conic_constant(), 0.0);
            assert_eq!(zero.even_aspheric_coefficients(), &[0.0; 7]);
            zero
        })
        .collect();
    let same = original
        .clone()
        .with_sequential_lens(SequentialLens::try_new(550.0, &surfaces).unwrap())
        .unwrap();
    with_cx(|cx| {
        for sample in [LensSample::CENTER, LensSample::try_new(0.75, 0.5).unwrap()] {
            let a = original
                .generate_sequential_spectral_timed_beams_from_tangent_offsets(
                    cx,
                    0.025,
                    -0.013,
                    sample,
                    &[450.0, 650.0],
                )
                .unwrap();
            let b = same
                .generate_sequential_spectral_timed_beams_from_tangent_offsets(
                    cx,
                    0.025,
                    -0.013,
                    sample,
                    &[450.0, 650.0],
                )
                .unwrap();
            assert_eq!(a, b);
        }
    });
}

#[test]
fn aspheric_spectral_beams_keep_actual_bulk_paths_and_opt_in_group_transit() {
    let clear = camera(true, 0.0);
    let absorbed = camera(true, 2.0);
    let sphere = camera(false, 0.0);
    assert_eq!(
        clear.sequential_lens().unwrap().reference_focus().unwrap(),
        sphere.sequential_lens().unwrap().reference_focus().unwrap()
    );
    with_cx(|cx| {
        let sample = LensSample::try_new(0.75, 0.5).unwrap();
        let wavelengths = [450.0, 650.0];
        let steady = clear
            .generate_sequential_spectral_beams_from_tangent_offsets(
                cx,
                0.025,
                -0.013,
                sample,
                &wavelengths,
            )
            .unwrap();
        let timed = clear
            .generate_sequential_spectral_timed_beams_from_tangent_offsets(
                cx,
                0.025,
                -0.013,
                sample,
                &wavelengths,
            )
            .unwrap();
        let attenuated = absorbed
            .generate_sequential_spectral_timed_beams_from_tangent_offsets(
                cx,
                0.025,
                -0.013,
                sample,
                &wavelengths,
            )
            .unwrap();
        let spherical = sphere
            .generate_sequential_spectral_beams_from_tangent_offsets(
                cx,
                0.025,
                -0.013,
                sample,
                &wavelengths,
            )
            .unwrap();
        let x = 0.006_f64;
        let rear = Point3::new(-x, 0.0, 0.08 + 10.0 * x * x - 1e4 * x.powi(4));
        let sensor = Point3::new(0.025 * 0.08, 0.013 * 0.08, 0.0);
        let sensor_length = rear.delta_from(sensor).norm();
        for (j, wavelength) in wavelengths.into_iter().enumerate() {
            let (a, b, dark) = (
                steady[j].unwrap(),
                timed[j].unwrap(),
                attenuated[j].unwrap(),
            );
            assert_eq!(a.ray(), b.ray());
            assert_eq!(a.frame(), b.frame());
            assert_eq!(a.power(), b.power());
            assert_eq!(a.geometry_weight(), b.geometry_weight());
            assert_eq!(a.group_delay_s(), None);
            assert_eq!(b.ray(), dark.ray());
            assert_eq!(b.group_delay_s(), dark.group_delay_s());
            let glass_length = b.ray().origin.delta_from(rear).norm();
            let group_index = 1.46 + 3.0 * 0.004 / (wavelength * 1e-3).powi(2);
            near(
                b.group_delay_s().unwrap().to_f64(),
                (sensor_length + group_index * glass_length) / 299_792_458.0,
                4e-24,
            );
            let attenuation = (-2.0 * glass_length).exp();
            let components = |s: fs_render::polarization::Stokes| [s.i, s.q, s.u, s.v];
            for (lossy, lossless) in components(dark.power())
                .into_iter()
                .zip(components(b.power()))
            {
                near(lossy, attenuation * lossless, 3e-14);
            }
            assert_ne!(b.ray().dir, spherical[j].unwrap().ray().dir);
        }
        assert_ne!(timed[0].unwrap().ray().dir, timed[1].unwrap().ray().dir);
    });
}

#[test]
fn complete_aspheric_bounds_refuse_overlapping_lens_surfaces() {
    let rear = LensSurface::try_new_aspheric(
        0.08,
        0.0,
        0.0,
        [1e5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        0.012,
        None,
    )
    .unwrap();
    let front = LensSurface::try_new(0.081, 0.0, 0.012, None).unwrap();
    assert_eq!(
        SequentialLens::try_new(550.0, &[rear, front]),
        Err(CameraError::InvalidSequentialLens)
    );
    assert!(
        LensSurface::try_new_aspheric(
            1e-4,
            0.0,
            0.0,
            [-1e5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            0.012,
            None
        )
        .is_err()
    );
}
