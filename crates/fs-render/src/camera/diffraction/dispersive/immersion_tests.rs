use super::*;
use crate::camera::{Aperture, CameraProjection, DefocusedDiffraction, LayeredDiffraction};
use crate::dielectric::{DielectricIorSample, SellmeierIor, TabulatedIor};
use fs_exec::{CancelGate, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}

fn constant(index: f64) -> DielectricIor {
    SellmeierIor::try_new_in_band([index * index - 1.0, 0.0, 0.0], [0.0; 3], (8000.0, 14000.0))
        .unwrap()
        .into()
}

fn exterior(scale: f64, low: f64, high: f64) -> DielectricIor {
    let gap = (high - low) / 8.0;
    TabulatedIor::try_new_in_band(
        core::array::from_fn(|i| {
            let wavelength = low + gap * i as f64;
            DielectricIorSample::try_new_in_band(
                wavelength,
                scale * (1.2 + (11000.0 - wavelength) * 1e-5),
            )
            .unwrap()
        }),
        gap,
    )
    .unwrap()
    .into()
}

fn lens(scale: f64) -> DispersiveThinLens {
    DispersiveThinLens::try_new_in_medium(
        11000.0,
        &[ThinLensElement::try_new(constant(2.0 * scale), 20.0, -20.0).unwrap()],
        exterior(scale, 8000.0, 14000.0),
    )
    .unwrap()
}

fn camera(lens: &DispersiveThinLens, aperture: Aperture) -> PhysicalCamera {
    PhysicalCamera::try_look_at(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.0, 0.0, 2.0),
        Vec3::new(0.0, 1.0, 0.0),
        CameraProjection::try_focal_sensor(lens.reference_focal_length_m(), 0.002).unwrap(),
        2.0,
        aperture,
    )
    .unwrap()
}

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 21,
                kernel_id: 3,
                tile: 0,
                iteration: 0,
            },
            asupersync::types::Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}

#[test]
fn g1_immersed_lensmaker_power_agrees_with_two_independent_snell_surfaces() {
    for (glass, ambient, front, back) in [
        (1.5, 4.0 / 3.0, 20.0, -20.0),
        (1.5, 2.0, -20.0, 20.0),
        (4.0, 1.5, 6.0, -6.0),
    ] {
        let lens = DispersiveThinLens::try_new_in_medium(
            11000.0,
            &[ThinLensElement::try_new(constant(glass), front, back).unwrap()],
            constant(ambient),
        )
        .unwrap();
        // Two exact angular Snell transitions at a common small ray height;
        // their thin paraxial limit is independent of the implementation.
        let height = 1e-8;
        let first_normal = -(height * front).asin();
        let inside = first_normal + ((ambient / glass) * (-first_normal).sin()).asin();
        let second_normal = -(height * back).asin();
        let outside = second_normal + ((glass / ambient) * (inside - second_normal).sin()).asin();
        close(lens.reference_power_per_m(), -outside.tan() / height, 1e-10);
        close(
            lens.sensor_distance_m(2.0).unwrap(),
            1.0 / (((glass / ambient) - 1.0) * (front - back) - 0.5),
            2e-15,
        );
    }
    let matched = ThinLensElement::try_new(constant(1.5), 20.0, -20.0).unwrap();
    assert!(DispersiveThinLens::try_new_in_medium(11000.0, &[matched], constant(1.5)).is_err());
}

#[test]
fn g1_dispersive_exterior_alone_changes_each_native_pupil_ray() {
    with_cx(|cx| {
        let lens = lens(1.0);
        for aperture in [
            Aperture::try_circular(0.001).unwrap(),
            Aperture::try_annular(0.0015, 0.0002).unwrap(),
            Aperture::try_regular_polygon(0.001, 5, 0.31).unwrap(),
        ] {
            let plain = camera(&lens, aperture);
            let camera = plain
                .clone()
                .with_geometric_dispersive_lens(lens.clone())
                .unwrap();
            let wavelengths = [8000.0, 11000.0, 14000.0];
            let sample = LensSample::try_new(0.83, 0.27).unwrap();
            let rays = camera
                .generate_geometric_spectral_rays_from_tangent_offsets(
                    cx,
                    0.013,
                    -0.007,
                    sample,
                    &wavelengths,
                )
                .unwrap();
            let reference = plain
                .generate_ray_from_tangent_offsets(cx, 0.013, -0.007, sample)
                .unwrap();
            assert_eq!(rays[1], reference);
            assert_ne!(rays[0].dir, rays[2].dir);
            for (ray, wavelength) in rays.iter().zip(wavelengths) {
                let ambient = 1.2 + (11000.0 - wavelength) * 1e-5;
                let power = (2.0 / ambient - 1.0) * 40.0;
                let inverse_sensor = (2.0 / 1.2 - 1.0) * 40.0 - 0.5;
                let pupil = ray.origin.delta_from(camera.eye());
                let axial = ray.dir.dot(camera.forward());
                close(
                    ray.dir.dot(camera.right()) / axial,
                    0.013 + pupil.dot(camera.right()) * (inverse_sensor - power),
                    2e-15,
                );
                close(
                    ray.dir.dot(camera.up()) / axial,
                    -0.007 + pupil.dot(camera.up()) * (inverse_sensor - power),
                    2e-15,
                );
            }
        }
    });
}

#[test]
fn g3_common_index_scaling_preserves_relative_power_and_geometric_rays() {
    with_cx(|cx| {
        let original = lens(1.0);
        let scaled = lens(2.0);
        let aperture = Aperture::try_circular(0.001).unwrap();
        let a = camera(&original, aperture.clone())
            .with_geometric_dispersive_lens(original.clone())
            .unwrap();
        let b = camera(&scaled, aperture)
            .with_geometric_dispersive_lens(scaled.clone())
            .unwrap();
        let wavelengths = [8000.0, 9375.0, 11000.0, 14000.0];
        for wavelength in wavelengths {
            assert_eq!(
                original.power_per_m(wavelength).unwrap(),
                scaled.power_per_m(wavelength).unwrap()
            );
            assert_eq!(
                original.power_difference_per_m(wavelength).unwrap(),
                scaled.power_difference_per_m(wavelength).unwrap()
            );
        }
        let sample = LensSample::try_new(0.72, 0.31).unwrap();
        assert_eq!(
            a.generate_geometric_spectral_rays_from_tangent_offsets(
                cx,
                0.01,
                -0.02,
                sample,
                &wavelengths
            )
            .unwrap(),
            b.generate_geometric_spectral_rays_from_tangent_offsets(
                cx,
                0.01,
                -0.02,
                sample,
                &wavelengths
            )
            .unwrap(),
        );
    });
}

#[test]
fn g1_weak_relative_chromatic_power_survives_equal_rounded_index_ratios() {
    let delta = f64::EPSILON;
    let table = |initial: f64| {
        TabulatedIor::try_new_in_band(
            core::array::from_fn(|i| {
                DielectricIorSample::try_new_in_band(
                    8000.0 + 750.0 * i as f64,
                    if i == 0 { initial } else { initial + delta },
                )
                .unwrap()
            }),
            750.0,
        )
        .unwrap()
    };
    let lens = DispersiveThinLens::try_new_in_medium(
        8000.0,
        &[ThinLensElement::try_new(table(1.5).into(), 20.0, -20.0).unwrap()],
        table(1.25).into(),
    )
    .unwrap();
    assert_eq!(1.5 / 1.25 - (1.5 + delta) / (1.25 + delta), 0.0);
    let expected = 40.0 * delta * 0.25 / (1.25 * (1.25 + delta));
    close(
        lens.power_difference_per_m(14000.0).unwrap(),
        expected,
        2e-30,
    );
    assert!(lens.power_difference_per_m(14000.0).unwrap() > 0.0);
}

#[test]
fn g0_immersed_source_support_scene_match_and_wave_policy_are_explicit() {
    let law = exterior(1.0, 9000.0, 13000.0);
    let lens = DispersiveThinLens::try_new_in_medium(
        11000.0,
        &[ThinLensElement::try_new(constant(2.0), 20.0, -20.0).unwrap()],
        law,
    )
    .unwrap();
    assert_eq!(lens.wavelength_bounds_nm(), (9000.0, 13000.0));
    assert_eq!(lens.exterior_ior(), Some(&law));
    for outside in [8999.0, 13001.0, 0.0, f64::NAN] {
        assert!(lens.power_per_m(outside).is_err());
        assert!(lens.power_difference_per_m(outside).is_err());
    }
    let actual = law.eval(11000.0).unwrap();
    lens.admit_exterior_phase_index(11000.0, actual).unwrap();
    assert!(
        lens.admit_exterior_phase_index(11000.0, actual.next_up())
            .is_err()
    );
    assert!(lens.admit_exterior_phase_index(11000.0, f64::NAN).is_err());
    let plain = camera(&lens, Aperture::try_circular(0.001).unwrap());
    let geometric = plain
        .clone()
        .with_geometric_dispersive_lens(lens.clone())
        .unwrap();
    let defocus = DefocusedDiffraction::try_new(Some(2.0), 0.01).unwrap();
    let layers = LayeredDiffraction::try_new(&[2.0], 0.01).unwrap();
    assert!(
        plain
            .clone()
            .with_defocused_diffraction(defocus)
            .unwrap()
            .with_dispersive_lens(lens.clone())
            .is_err()
    );
    assert!(
        plain
            .with_layered_diffraction(layers)
            .unwrap()
            .with_dispersive_lens(lens)
            .is_err()
    );
    assert!(
        geometric
            .clone()
            .with_defocused_diffraction(defocus)
            .is_err()
    );
    assert!(geometric.with_layered_diffraction(layers).is_err());
}

#[test]
fn g1_geometric_opening_uses_medium_wavelength_without_changing_vacuum_source_coordinates() {
    with_cx(|cx| {
        let immersed = DispersiveThinLens::try_new_in_medium(
            11000.0,
            &[ThinLensElement::try_new(constant(5.0), 20.0, -20.0).unwrap()],
            constant(4.0),
        )
        .unwrap();
        let vacuum = DispersiveThinLens::try_new(
            11000.0,
            &[ThinLensElement::try_new(constant(5.0), 1.25, -1.25).unwrap()],
        )
        .unwrap();
        assert_eq!(
            immersed.reference_focal_length_m(),
            vacuum.reference_focal_length_m()
        );
        let aperture = Aperture::try_circular(0.0003).unwrap();
        let a = camera(&immersed, aperture.clone())
            .with_geometric_dispersive_lens(immersed)
            .unwrap();
        let b = camera(&vacuum, aperture)
            .with_geometric_dispersive_lens(vacuum)
            .unwrap();
        let wavelengths = [8000.0, 11000.0, 14000.0];
        assert!(
            a.generate_geometric_spectral_rays_from_tangent_offsets(
                cx,
                0.0,
                0.0,
                LensSample::CENTER,
                &wavelengths
            )
            .is_ok()
        );
        assert!(
            b.generate_geometric_spectral_rays_from_tangent_offsets(
                cx,
                0.0,
                0.0,
                LensSample::CENTER,
                &wavelengths
            )
            .is_err()
        );
    });
}
