use super::*;
use crate::polarization::dielectric::scene::curved::camera::SpatialVolumeEmitter;

fn temperature(p: Point3) -> Result<f64, crate::gradient_index::GradientIndexError> {
    let radiance = 1.0 + 2.0 * p.z;
    let h = 6.626_070_15e-34;
    let c = 299_792_458.0_f64;
    let k = 1.380_649e-23;
    let wavelength = 550e-9_f64;
    Ok((h * c / (wavelength * k))
        / (2.0 * h * c * c / wavelength.powi(5) * 1e-9 / radiance).ln_1p())
}
fn integral(length: f64, initial: f64, slope: f64) -> f64 {
    let survival = (-0.2 * length).exp();
    let absorbed = -(-0.2 * length).exp_m1();
    initial * absorbed + slope * (absorbed / 0.2 - length * survival)
}

#[test]
fn g1_nonisothermal_emission_keeps_polarized_collision_weights_and_native_survival() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(2, 1.0)]).unwrap();
        let material = glass(1.4, 0.2);
        let media = [SceneMedium::Homogeneous(&material)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, binding(), 100_000).unwrap();
        let source =
            SpatialVolumeEmitter::try_new(0, &temperature, ContentHash([110; 32]), 1e-12).unwrap();
        let camera = CurvedCameraScene::try_new(cx, &native, &[], None)
            .unwrap()
            .with_spatial_volume_emitters(cx, &[source])
            .unwrap();
        let spectrum = MassScatteringSpectrum::constant(1.0).unwrap();
        for delta in [0.0, 0.03] {
            let phase = PhaseFunction::rayleigh_depolarized(delta).unwrap();
            let laws = [Some(
                CurvedRayleighMedium::try_new(&spectrum, 1.0, &phase).unwrap(),
            )];
            let scattering = CurvedRayleighScene::try_new(cx, &native, &laws).unwrap();
            let molecular = CurvedRayleighCameraScene::try_new(cx, &camera, &scattering).unwrap();
            let isotropic = 3.0 * delta / (2.0 + delta);
            let right_angle = isotropic + 0.5 * (1.0 - isotropic) + f64::EPSILON;
            let flights = |ordinal| match ordinal {
                0 => [1.0 - (-0.25_f64).exp(), right_angle, 0.0],
                1 => [1.0 - (-0.5_f64).exp(), right_angle, 0.5],
                2 => [1.0 - 1e-12, 0.5, 0.0],
                _ => panic!("opaque boundary precedes next collision"),
            };
            let sample = molecular
                .sample_ray_with_samples(
                    cx,
                    ray(),
                    550.0,
                    0,
                    settings(),
                    |_| panic!("opaque surface uses no sample"),
                    flights,
                )
                .unwrap();
            let polarization = (1.0 - delta) / (1.0 + delta);
            // Along +z,.25m: B=1+2s; then +x,.5m: B=1.5;
            // finally +z,.75m: B=1.5+2s, with the twice-polarized prefix.
            let expected = integral(0.25, 1.0, 2.0)
                + (-0.05_f64).exp() * integral(0.5, 1.5, 0.0)
                + (-0.15_f64).exp()
                    * (1.0 + polarization * polarization)
                    * integral(0.75, 1.5, 2.0);
            close(
                sample.radiance_sample().radiance_w_m2_sr_nm() / 1.4_f64.powi(2),
                expected,
                1e-6,
            );
            assert_eq!(sample.volume_events(), 2);
            assert_eq!(sample.free_flight_draws(), 3);
            let tracked = molecular
                .sample_ray_with_tracking_samples(
                    cx,
                    ray(),
                    550.0,
                    0,
                    settings(),
                    |_| panic!("opaque surface uses no sample"),
                    |ordinal| {
                        let draw = flights(ordinal);
                        [draw[0], draw[1], draw[2], 0.5]
                    },
                )
                .unwrap();
            assert_eq!(sample, tracked);
        }
    });
}

#[test]
fn g1_null_candidates_and_zero_scattering_preserve_varying_source_integrals() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(2, 1.0)]).unwrap();
        let material = glass(1.2, 0.2);
        let media = [SceneMedium::Homogeneous(&material)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, binding(), 100_000).unwrap();
        let source =
            SpatialVolumeEmitter::try_new(0, &temperature, ContentHash([111; 32]), 1e-12).unwrap();
        let camera = CurvedCameraScene::try_new(cx, &native, &[], None)
            .unwrap()
            .with_spatial_volume_emitters(cx, &[source])
            .unwrap();
        let spectrum = MassScatteringSpectrum::constant(1.0).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let laws = [Some(
            CurvedRayleighMedium::try_new_heterogeneous(
                &spectrum,
                &VacuumDensity,
                1.0,
                ContentHash([79; 32]),
                &phase,
            )
            .unwrap(),
        )];
        let scattering = CurvedRayleighScene::try_new(cx, &native, &laws).unwrap();
        let molecular = CurvedRayleighCameraScene::try_new(cx, &camera, &scattering).unwrap();
        let sample = molecular
            .sample_ray_with_tracking_samples(
                cx,
                ray(),
                550.0,
                0,
                settings(),
                |_| panic!("no surface sample"),
                |ordinal| match ordinal {
                    0 => [1.0 - (-0.25_f64).exp(), 0.5, 0.0, 0.5],
                    1 => [0.99, 0.5, 0.0, 0.5],
                    _ => panic!("opaque boundary precedes the next candidate"),
                },
            )
            .unwrap();
        close(
            sample.radiance_sample().radiance_w_m2_sr_nm() / 1.2_f64.powi(2),
            integral(1.0, 1.0, 2.0),
            1e-6,
        );
        assert_eq!(sample.volume_events(), 0);
        assert_eq!(sample.free_flight_draws(), 2);
        let clear_laws = [Some(
            CurvedRayleighMedium::try_new(&spectrum, 0.0, &phase).unwrap(),
        )];
        let clear = CurvedRayleighScene::try_new(cx, &native, &clear_laws).unwrap();
        let molecular = CurvedRayleighCameraScene::try_new(cx, &camera, &clear).unwrap();
        let ordinary = camera
            .sample_ray(cx, ray(), 550.0, 0, settings().surfaces, 41, 3)
            .unwrap();
        let actual = molecular
            .sample_ray(cx, ray(), 550.0, 0, settings(), 41, 3)
            .unwrap();
        assert_eq!(*actual.radiance_sample(), ordinary);
        assert_eq!(actual.volume_events(), 0);
        assert_eq!(actual.free_flight_draws(), 0);
    });
}
