use super::*;
use crate::polarization::dielectric::scene::curved::camera::IsothermalVolumeEmitter;
use crate::thermal_surface::planck_radiance_w_m2_sr_nm;

#[test]
fn g1_lte_emission_uses_each_pre_collision_mueller_prefix_without_extra_scattering_extinction() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(2, 1.0)]).unwrap();
        let material = glass(1.4, 0.2);
        let media = [SceneMedium::Homogeneous(&material)];
        let native = CurvedPolarizedScene::try_new(cx, &geometry, &media, binding(), 100_000).unwrap();
        let emitter = IsothermalVolumeEmitter::try_new(0, 2400.0, ContentHash([93; 32])).unwrap();
        let camera = CurvedCameraScene::try_new(cx, &native, &[], None).unwrap()
            .with_isothermal_volume_emitters(cx, &[emitter]).unwrap();
        let spectrum = MassScatteringSpectrum::constant(1.0).unwrap();
        for delta in [0.0, 0.03] {
            let phase = PhaseFunction::rayleigh_depolarized(delta).unwrap();
            let laws = [Some(CurvedRayleighMedium::try_new(&spectrum, 1.0, &phase).unwrap())];
            let scattering = CurvedRayleighScene::try_new(cx, &native, &laws).unwrap();
            let molecular = CurvedRayleighCameraScene::try_new(cx, &camera, &scattering).unwrap();
            let isotropic = 3.0 * delta / (2.0 + delta);
            let right_angle = isotropic + 0.5 * (1.0 - isotropic) + f64::EPSILON;
            let flights = |ordinal| match ordinal {
                0 => [1.0 - (-0.25_f64).exp(), right_angle, 0.0],
                1 => [1.0 - (-0.5_f64).exp(), right_angle, 0.5],
                2 => [1.0 - 1e-12, 0.5, 0.0],
                _ => panic!("opaque boundary precedes the next collision"),
            };
            let sample = molecular.sample_ray_with_samples(cx, ray(), 550.0, 0, settings(),
                |_| panic!("no surface sample"), flights).unwrap();
            // First two legs are .25m and .5m with unit phase weight. The
            // second polarized 90-degree event multiplies subsequent power
            // by 1+p², so the last .75m has a different thermal prefix.
            let polarization = (1.0 - delta) / (1.0 + delta);
            let transmission = (-0.15_f64).exp();
            let expected = (1.0 - transmission)
                * (1.0 + (1.0 + polarization * polarization) * transmission);
            let blackbody = planck_radiance_w_m2_sr_nm(2400.0, 550.0).unwrap();
            close(sample.radiance_sample().radiance_w_m2_sr_nm() / (1.4_f64.powi(2) * blackbody),
                expected, 3e-7);
            assert_eq!(sample.volume_events(), 2);
            assert_eq!(sample.free_flight_draws(), 3);
            let tracked = molecular.sample_ray_with_tracking_samples(cx, ray(), 550.0, 0, settings(),
                |_| panic!("no surface sample"), |ordinal| {
                    let draw = flights(ordinal);
                    [draw[0], draw[1], draw[2], 0.5]
                }).unwrap();
            assert_eq!(sample, tracked);
        }
    });
}

#[test]
fn g1_null_candidates_preserve_the_full_lte_absorption_integral() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(2, 1.0)]).unwrap();
        let material = glass(1.2, 0.2);
        let media = [SceneMedium::Homogeneous(&material)];
        let native = CurvedPolarizedScene::try_new(cx, &geometry, &media, binding(), 100_000).unwrap();
        let emitter = IsothermalVolumeEmitter::try_new(0, 2400.0, ContentHash([94; 32])).unwrap();
        let camera = CurvedCameraScene::try_new(cx, &native, &[], None).unwrap()
            .with_isothermal_volume_emitters(cx, &[emitter]).unwrap();
        let spectrum = MassScatteringSpectrum::constant(1.0).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let laws = [Some(CurvedRayleighMedium::try_new_heterogeneous(&spectrum,
            &VacuumDensity, 1.0, ContentHash([79; 32]), &phase).unwrap())];
        let scattering = CurvedRayleighScene::try_new(cx, &native, &laws).unwrap();
        let molecular = CurvedRayleighCameraScene::try_new(cx, &camera, &scattering).unwrap();
        let sample = molecular.sample_ray_with_tracking_samples(cx, ray(), 550.0, 0, settings(),
            |_| panic!("no surface sample"), |ordinal| match ordinal {
                0 => [1.0 - (-0.25_f64).exp(), 0.5, 0.0, 0.5],
                1 => [0.99, 0.5, 0.0, 0.5],
                _ => panic!("opaque boundary precedes the next candidate"),
            }).unwrap();
        let blackbody = planck_radiance_w_m2_sr_nm(2400.0, 550.0).unwrap();
        close(sample.radiance_sample().radiance_w_m2_sr_nm() / (1.2_f64.powi(2) * blackbody),
            -(-0.2_f64).exp_m1(), 2e-7);
        assert_eq!(sample.volume_events(), 0);
        assert_eq!(sample.free_flight_draws(), 2);
    });
}
