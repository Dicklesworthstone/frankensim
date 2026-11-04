//! Source-bounded phase power must change the actual geometric camera rays.

use super::*;
use fs_render::camera::{DispersiveThinLens, ThinLensElement};
use fs_render::dielectric::{DielectricIorSample, TabulatedIor};

fn infrared_lens() -> DispersiveThinLens {
    let law = TabulatedIor::try_new_in_band(
        core::array::from_fn(|i| {
            let wavelength = 8000.0 + 750.0 * i as f64;
            DielectricIorSample::try_new_in_band(wavelength, 1.5 + (11000.0 - wavelength) * 2.0e-5)
                .unwrap()
        }),
        750.0,
    )
    .unwrap();
    DispersiveThinLens::try_new(
        11000.0,
        &[ThinLensElement::try_new(law.into(), 20.0, -20.0).unwrap()],
    )
    .unwrap()
}

fn infrared_camera(aperture: Aperture) -> PhysicalCamera {
    PhysicalCamera::try_look_at(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.0, 0.0, -1.0),
        Vec3::new(0.0, 1.0, 0.0),
        CameraProjection::try_focal_sensor(0.05, 0.002).unwrap(),
        2.0,
        aperture,
    )
    .unwrap()
}

#[test]
fn g1_infrared_phase_lensmaker_law_moves_rays_through_every_original_pupil() {
    with_cx(|_, cx| {
        let lens = infrared_lens();
        assert_eq!(lens.wavelength_bounds_nm(), (8000.0, 14000.0));
        for pupil in [
            Aperture::try_circular(0.002).unwrap(),
            Aperture::try_annular(0.002, 0.0005).unwrap(),
            Aperture::try_regular_polygon(0.002, 5, 0.31).unwrap(),
        ] {
            let original = infrared_camera(pupil);
            let camera = original
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
            let reference = original
                .generate_ray_from_tangent_offsets(cx, 0.013, -0.007, sample)
                .unwrap();
            assert_eq!(rays[1], reference);
            assert_ne!(rays[0].dir, rays[1].dir);
            assert_ne!(rays[1].dir, rays[2].dir);
            for (ray, wavelength) in rays.iter().zip(wavelengths) {
                assert_eq!(ray.origin, reference.origin);
                let [px, py] = lens_coordinates(&camera, ray.origin);
                // Independent lensmaker/conjugate calculation in SI metres.
                let phase_index = 1.5 + (11000.0 - wavelength) * 2.0e-5;
                let power = 40.0 * (phase_index - 1.0);
                let inverse_sensor = 20.0 - 0.5;
                let axial = ray.dir.dot(camera.forward());
                assert_close(
                    ray.dir.dot(camera.right()) / axial,
                    0.013 + px * (inverse_sensor - power),
                    2e-15,
                    "horizontal phase ray",
                );
                assert_close(
                    ray.dir.dot(camera.up()) / axial,
                    -0.007 + py * (inverse_sensor - power),
                    2e-15,
                    "vertical phase ray",
                );
                let intercept = ray.origin.offset(ray.dir.scale(2.0 / axial));
                assert_close(
                    intercept.delta_from(camera.eye()).dot(camera.right()),
                    2.0 * 0.013 + 2.0 * px * (20.0 - power),
                    5e-15,
                    "chromatic focal blur",
                );
            }
        }
    });
}

#[test]
fn g0_geometric_dispersion_refuses_source_exits_shared_rays_wave_mixing_and_cancellation() {
    with_cx(|gate, cx| {
        let original = infrared_camera(Aperture::try_circular(0.001).unwrap());
        let camera = original
            .clone()
            .with_geometric_dispersive_lens(infrared_lens())
            .unwrap();
        for wavelength in [7999.0, 14001.0, 0.0, f64::NAN] {
            assert!(
                camera
                    .generate_geometric_spectral_rays_from_tangent_offsets(
                        cx,
                        0.0,
                        0.0,
                        LensSample::CENTER,
                        &[wavelength],
                    )
                    .is_err()
            );
        }
        assert_eq!(
            camera
                .generate_ray_from_tangent_offsets(cx, 0.0, 0.0, LensSample::CENTER,)
                .unwrap_err(),
            CameraError::DispersionRequiresSpectralRays
        );
        assert_eq!(
            camera
                .generate_spectral_ray_from_tangent_offsets(
                    cx,
                    0.0,
                    0.0,
                    LensSample::CENTER,
                    LensSample::CENTER,
                    &[9000.0, 13000.0],
                )
                .unwrap_err(),
            CameraError::DispersionRequiresSpectralRays
        );
        assert!(
            camera
                .generate_geometric_spectral_rays_from_tangent_offsets(
                    cx,
                    0.1,
                    0.1,
                    LensSample::CENTER,
                    &[11000.0],
                )
                .is_err()
        );
        let wave = original
            .with_defocused_diffraction(DefocusedDiffraction::try_new(Some(2.0), 0.0).unwrap())
            .unwrap();
        assert!(
            wave.clone()
                .with_geometric_dispersive_lens(infrared_lens())
                .is_err()
        );
        assert_eq!(
            wave.generate_geometric_spectral_rays_from_tangent_offsets(
                cx,
                0.0,
                0.0,
                LensSample::CENTER,
                &[11000.0],
            )
            .unwrap_err(),
            CameraError::DiffractionRequiresSpectralPath
        );
        // The existing wave owner still accepts the same source-bound lens.
        assert!(
            wave.with_dispersive_lens(infrared_lens())
                .unwrap()
                .defocused_diffraction_density([0.0; 2], 11000.0)
                .is_ok()
        );
        let undersized = infrared_camera(Aperture::try_circular(1e-5).unwrap())
            .with_geometric_dispersive_lens(infrared_lens())
            .unwrap();
        assert!(
            undersized
                .generate_geometric_spectral_rays_from_tangent_offsets(
                    cx,
                    0.0,
                    0.0,
                    LensSample::CENTER,
                    &[11000.0],
                )
                .is_err()
        );
        gate.request();
        assert_eq!(
            camera
                .generate_geometric_spectral_rays_from_tangent_offsets(
                    cx,
                    0.0,
                    0.0,
                    LensSample::CENTER,
                    &[11000.0],
                )
                .unwrap_err(),
            CameraError::Cancelled
        );
    });
}
