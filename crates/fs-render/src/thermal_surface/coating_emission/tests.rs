use super::*;
use crate::conductor::{
    CoatingIorSample, ConductorDataStatus, ConductorIorSample, ConductorSource,
};
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricGlass, GlassProvenance};
use core::f64::consts::PI;

mod infrared;

fn source(tag: u8) -> ConductorSource {
    ConductorSource::try_new(ContentHash([tag; 32]), ConductorDataStatus::Representative).unwrap()
}

fn film(eta: f64, k: f64, thickness_m: f64, tag: u8) -> CoherentAbsorbingCoating {
    CoherentAbsorbingCoating::try_new(
        core::array::from_fn(|i| {
            CoatingIorSample::try_new(380.0 + 50.0 * i as f64, eta, k).unwrap()
        }),
        thickness_m,
        source(tag),
    )
    .unwrap()
}

fn substrate(eta: f64, k: f64) -> ConductorOptics {
    ConductorOptics::try_new(
        core::array::from_fn(|i| {
            ConductorIorSample::try_new(380.0 + 50.0 * i as f64, eta, k).unwrap()
        }),
        source(4),
    )
    .unwrap()
}

fn vacuum() -> DielectricIor {
    CauchyIor::try_constant(1.0).unwrap().into()
}

fn close(actual: f64, expected: f64) {
    let scale = expected.abs().max(f64::MIN_POSITIVE);
    assert!(
        (actual - expected).abs() <= 3.0e-11 * scale,
        "actual={actual:e}, expected={expected:e}"
    );
}

fn planck_si(temperature_k: f64, wavelength_nm: f64) -> f64 {
    // Independent direct SI expression in a representable temperature range.
    // Production uses the native logarithmic Planck owner instead.
    let wavelength_m = wavelength_nm * 1.0e-9;
    let h = 6.626_070_15e-34;
    let c = 299_792_458.0;
    let k = 1.380_649e-23;
    2.0 * h * c * c / (wavelength_m.powi(5) * (h * c / (wavelength_m * k * temperature_k)).exp_m1())
        * 1.0e-9
}

#[test]
fn g1_depth_resolved_slab_emits_original_local_temperatures() {
    let glass = DielectricGlass::new(
        CauchyIor::try_constant(1.5).unwrap(),
        BeerLambertAbsorption::try_constant(40.0).unwrap(),
        GlassProvenance::Custom,
    );
    let depths = [0.005, 0.02];
    let temperatures = [1500.0, 2500.0];
    let emitter = CoatingThermalEmission::try_depth_resolved_incoherent_slab(
        glass,
        0.02,
        &depths,
        &temperatures,
        vacuum(),
        ContentHash([93; 32]),
    )
    .unwrap();
    let sample = emitter.evaluate(550.0, 1.0).unwrap();
    assert_eq!(emitter.film_count(), 2);
    assert_eq!(emitter.film_temperature_k(1), Some(2500.0));
    assert_eq!(emitter.film_temperature_k(2), None);
    // Independent normal-incidence Fresnel and forward/backward Beer series.
    let reflectance = (0.5_f64 / 2.5).powi(2);
    let entering =
        (1.0 - reflectance) / (1.0 - reflectance.powi(2) * (-2.0_f64 * 40.0 * 0.02).exp());
    let mut expected = 0.0;
    let mut start = 0.0;
    for (i, &end) in depths.iter().enumerate() {
        let absorbed = entering
            * -(-40.0 * (end - start)).exp_m1()
            * ((-40.0 * start).exp() + reflectance * (-40.0 * (0.04 - end)).exp());
        let radiance = absorbed * planck_si(temperatures[i], 550.0);
        for mode in sample.film_radiance_w_m2_sr_nm()[i] {
            close(mode, radiance / 2.0);
        }
        expected += radiance;
        start = end;
    }
    close(sample.stokes().i, expected);
    let average = CoatingThermalEmission::try_incoherent_slab(
        glass,
        0.02,
        2250.0,
        vacuum(),
        ContentHash([93; 32]),
    )
    .unwrap();
    assert!(
        (sample.stokes().i / average.evaluate(550.0, 1.0).unwrap().stokes().i - 1.0).abs() > 0.1
    );
    // Reversing the actual hot/cold profile changes directional emission,
    // while the homogeneous optical R/T/A remains unchanged.
    let reverse = CoatingThermalEmission::try_depth_resolved_incoherent_slab(
        glass,
        0.02,
        &[0.015, 0.02],
        &[2500.0, 1500.0],
        vacuum(),
        ContentHash([93; 32]),
    )
    .unwrap()
    .evaluate(550.0, 1.0)
    .unwrap();
    assert!(reverse.stokes().i > sample.stokes().i);
    assert_eq!(
        reverse.reciprocal_absorbed_fraction(),
        sample.reciprocal_absorbed_fraction()
    );
}

#[test]
fn g3_depth_partition_retains_isothermal_emission_and_temperature_only_cache_identity() {
    let glass = DielectricGlass::new(
        CauchyIor::try_constant(1.5).unwrap(),
        BeerLambertAbsorption::try_constant(40.0).unwrap(),
        GlassProvenance::Custom,
    );
    let other: DielectricIor = CauchyIor::try_constant(1.33).unwrap().into();
    let prepare = |depths: &[f64], temperatures: &[f64]| {
        CoatingThermalEmission::try_depth_resolved_incoherent_slab(
            glass,
            0.02,
            depths,
            temperatures,
            other,
            ContentHash([94; 32]),
        )
        .unwrap()
    };
    let depths = [0.001, 0.008, 0.02];
    let partitioned = prepare(&depths, &[1500.0; 3]);
    let uniform = CoatingThermalEmission::try_incoherent_slab(
        glass,
        0.02,
        1500.0,
        other,
        ContentHash([94; 32]),
    )
    .unwrap();
    for cosine in [1.0, 0.5, 0.01] {
        let result = partitioned.evaluate(550.0, cosine).unwrap();
        let reference = uniform.evaluate(550.0, cosine).unwrap();
        for pol in 0..2 {
            close(
                result.polarized_radiance_w_m2_sr_nm()[pol],
                reference.polarized_radiance_w_m2_sr_nm()[pol],
            );
        }
        let work = partitioned
            .reciprocal_layer_absorption(550.0, cosine)
            .unwrap();
        for pol in 0..2 {
            close(
                work.films[..3].iter().map(|value| value[pol]).sum(),
                work.total[pol],
            );
        }
    }
    let warmer = prepare(&depths, &[1400.0, 1600.0, 1700.0]);
    assert_eq!(warmer.optical_identity(), partitioned.optical_identity());
    assert_ne!(warmer.identity(), partitioned.identity());
    assert_ne!(
        prepare(&[0.002, 0.008, 0.02], &[1500.0; 3]).optical_identity(),
        partitioned.optical_identity()
    );
}

#[test]
fn g0_depth_emission_refuses_incomplete_or_unsupported_profiles() {
    let glass = DielectricGlass::new(
        CauchyIor::try_constant(1.5).unwrap(),
        BeerLambertAbsorption::try_constant(0.0).unwrap(),
        GlassProvenance::Custom,
    );
    let prepare = |depths: &[f64], temperatures: &[f64]| {
        CoatingThermalEmission::try_depth_resolved_incoherent_slab(
            glass,
            0.02,
            depths,
            temperatures,
            vacuum(),
            ContentHash([95; 32]),
        )
    };
    for (depths, temperatures) in [
        (vec![], vec![]),
        (vec![0.01, 0.02], vec![1500.0]),
        (vec![0.0, 0.02], vec![1500.0; 2]),
        (vec![0.01], vec![1500.0]),
        (vec![0.02, 0.01], vec![1500.0; 2]),
        (vec![f64::NAN], vec![1500.0]),
        (vec![0.02], vec![f64::NAN]),
        (vec![0.02], vec![0.0]),
        (
            vec![0.02; MAX_COHERENT_STACK_LAYERS + 1],
            vec![1500.0; MAX_COHERENT_STACK_LAYERS + 1],
        ),
    ] {
        assert!(prepare(&depths, &temperatures).is_err());
    }
    let transparent = prepare(&[0.01, 0.02], &[1500.0, 2500.0]).unwrap();
    assert_eq!(transparent.evaluate(550.0, 0.7).unwrap().stokes().i, 0.0);
    assert!(transparent.evaluate(300.0, 0.7).is_err());
}

#[test]
fn incoherent_slab_retains_one_physical_receiver_and_directional_kirchhoff_emission() {
    let glass = DielectricGlass::new(
        CauchyIor::try_constant(1.5).unwrap(),
        BeerLambertAbsorption::try_constant(40.0).unwrap(),
        GlassProvenance::Custom,
    );
    let temperature = 1500.0;
    let emitter = CoatingThermalEmission::try_incoherent_slab(
        glass,
        0.02,
        temperature,
        CauchyIor::try_constant(1.33).unwrap().into(),
        ContentHash([91; 32]),
    )
    .unwrap();
    assert_eq!(emitter.film_count(), 1);
    assert_eq!(emitter.film_temperature_k(0), Some(temperature));
    assert_eq!(emitter.film_temperature_k(1), None);
    assert_eq!(emitter.substrate_temperature_k(), None);
    assert!(!emitter.is_vacuum_transmitting());
    for cosine in [1.0, 0.6, 0.1] {
        let sample = emitter.evaluate(550.0, cosine).unwrap();
        let optical = glass
            .incoherent_slab_response_between(550.0, 0.02, cosine, 1.0, 1.33)
            .unwrap();
        let powers = [optical.perpendicular, optical.parallel];
        assert_eq!(
            sample.reciprocal_absorbed_fraction(),
            powers.map(|p| p.absorptance)
        );
        assert_eq!(
            sample.reciprocal_transmitted_fraction(),
            powers.map(|p| p.transmittance)
        );
        for pol in 0..2 {
            close(
                sample.polarized_radiance_w_m2_sr_nm()[pol],
                0.5 * powers[pol].absorptance * planck_si(temperature, 550.0),
            );
        }
        assert_eq!(
            sample.film_radiance_w_m2_sr_nm(),
            &[sample.polarized_radiance_w_m2_sr_nm()]
        );
        for opposite in [false, true] {
            assert!(
                emitter
                    .reciprocal_vacuum_port_power(550.0, cosine, opposite)
                    .unwrap()
                    .is_none()
            );
        }
    }
    let clone = emitter.clone();
    assert_eq!(clone.identity(), emitter.identity());
    assert_eq!(
        clone.evaluate(550.0, 0.6).unwrap(),
        emitter.evaluate(550.0, 0.6).unwrap()
    );
    let warmer = CoatingThermalEmission::try_incoherent_slab(
        glass,
        0.02,
        1700.0,
        CauchyIor::try_constant(1.33).unwrap().into(),
        ContentHash([91; 32]),
    )
    .unwrap();
    assert_eq!(warmer.optical_identity(), emitter.optical_identity());
    assert_ne!(warmer.identity(), emitter.identity());
    assert!(
        warmer.radiance_w_m2_sr_nm(550.0, 0.6).unwrap()
            > emitter.radiance_w_m2_sr_nm(550.0, 0.6).unwrap()
    );
}

#[test]
fn slab_emission_keeps_transparent_thin_and_black_limits_and_two_vacuum_ports() {
    use crate::thermal_surface::{
        CoatingExchangeQuadrature, CoatingTransmittingPairBudget, CoatingTransmittingPairResponse,
    };
    for alpha in [0.0, 1.0e-20, 1.0e5] {
        let glass = DielectricGlass::new(
            CauchyIor::try_constant(1.0).unwrap(),
            BeerLambertAbsorption::try_constant(alpha).unwrap(),
            GlassProvenance::Custom,
        );
        let emitter = CoatingThermalEmission::try_incoherent_slab(
            glass,
            0.02,
            1500.0,
            vacuum(),
            ContentHash([92; 32]),
        )
        .unwrap();
        let sample = emitter.evaluate(550.0, 0.6).unwrap();
        let absorption = -(-alpha * 0.02 / 0.6_f64).exp_m1();
        close(sample.stokes().i, absorption * planck_si(1500.0, 550.0));
        assert_eq!(sample.film_radiance_w_m2_sr_nm().len(), 1);
        if alpha == 0.0 {
            assert_eq!(sample.stokes().i, 0.0);
        } else {
            assert!(sample.stokes().i > 0.0);
        }
        let front = emitter
            .reciprocal_vacuum_port_power(550.0, 0.6, false)
            .unwrap()
            .unwrap();
        let rear = emitter
            .reciprocal_vacuum_port_power(550.0, 0.6, true)
            .unwrap()
            .unwrap();
        assert_eq!(front.absorbed, rear.absorbed);
        assert_eq!(front.transmitted, rear.transmitted);
        let response = CoatingTransmittingPairResponse::from_emitters(
            &emitter,
            &emitter,
            0.1,
            CoatingExchangeQuadrature {
                wavelength_bins: 4,
                cosine_bins: 4,
            },
            CoatingTransmittingPairBudget {
                max_channel_spectral_nodes: 24,
                max_optical_evaluations: 64,
            },
            || Ok(()),
        )
        .unwrap();
        let equilibrium = response
            .evaluate(1500.0, 1500.0, Some(1500.0), Some(1500.0), || Ok(()))
            .unwrap();
        assert_eq!(equilibrium.body_outward_power_w(), [0.0; 2]);
        assert_eq!(equilibrium.bath_outward_power_w(), [0.0; 2]);
        if alpha == 0.0 {
            let dark = response
                .evaluate(1200.0, 1800.0, Some(900.0), Some(1900.0), || Ok(()))
                .unwrap();
            assert_eq!(dark.body_outward_power_w(), [0.0; 2]);
        }
    }
}

#[test]
fn unequal_exterior_media_use_their_own_reciprocal_absorption_and_radiance_measure() {
    let glass = DielectricGlass::new(
        CauchyIor::try_constant(1.5).unwrap(),
        BeerLambertAbsorption::try_constant(40.0).unwrap(),
        GlassProvenance::Custom,
    );
    let temperature = 1500.0;
    let front_cosine = 0.6_f64;
    let rear_cosine = (1.0 - (1.0 - front_cosine.powi(2)) / 1.33_f64.powi(2)).sqrt();
    let front = glass
        .incoherent_slab_thermal_response_between(550.0, 0.02, front_cosine, 1.0, 1.33, temperature)
        .unwrap();
    let rear = glass
        .incoherent_slab_thermal_response_between(550.0, 0.02, rear_cosine, 1.33, 1.0, temperature)
        .unwrap();
    for (response, index) in [(front, 1.0_f64), (rear, 1.33_f64)] {
        for (pol, power) in [response.optical.perpendicular, response.optical.parallel]
            .into_iter()
            .enumerate()
        {
            close(
                response.polarized_radiance_w_m2_sr_nm[pol],
                0.5 * index.powi(2) * power.absorptance * planck_si(temperature, 550.0),
            );
            // In equal-temperature surroundings, each mode remains its local
            // blackbody radiance after reflection, reciprocal transmission and emission.
            close(
                index.powi(2)
                    * planck_si(temperature, 550.0)
                    * (power.reflectance + power.transmittance + power.absorptance),
                index.powi(2) * planck_si(temperature, 550.0),
            );
        }
    }
    close(
        front.optical.perpendicular.transmittance,
        rear.optical.perpendicular.transmittance,
    );
    assert_ne!(
        front.polarized_radiance_w_m2_sr_nm,
        rear.polarized_radiance_w_m2_sr_nm
    );
}

#[test]
fn nonisothermal_matched_layers_follow_independent_beer_lambert_and_planck_sources() {
    // Same complex index in both finite films and the half space: only the
    // entrance reflects. Each finite interval has a closed-form positive loss.
    // Different temperatures must weight those depths separately.
    let eta = 1.5;
    let k = 0.3;
    let thickness = [80.0e-9, 120.0e-9];
    let temperatures = [900.0, 1600.0];
    let substrate_temperature = 2300.0;
    let emitter = {
        let films = [film(eta, k, thickness[0], 1), film(eta, k, thickness[1], 2)];
        CoatingThermalEmission::try_opaque(
            &[
                CoatingLayer::Absorbing(&films[0]),
                CoatingLayer::Absorbing(&films[1]),
            ],
            &temperatures,
            substrate(eta, k),
            substrate_temperature,
            ContentHash([9; 32]),
        )
        .unwrap()
    }; // Both original borrowed tables are gone; the prepared emitter owns them.
    assert_eq!(emitter.film_count(), 2);
    assert_eq!(emitter.film_temperature_k(1), Some(1600.0));
    assert_eq!(emitter.substrate_temperature_k(), Some(2300.0));
    for wavelength in [380.0, 550.0, 780.0] {
        let entrance = 4.0 * eta / ((1.0 + eta) * (1.0 + eta) + k * k);
        let alpha = 4.0 * PI * k / (wavelength * 1.0e-9);
        let optical_depth = thickness.map(|d| alpha * d);
        let absorption = [
            entrance * -(-optical_depth[0]).exp_m1(),
            entrance * (-optical_depth[0]).exp() * -(-optical_depth[1]).exp_m1(),
        ];
        let substrate_absorption = entrance * (-(optical_depth[0] + optical_depth[1])).exp();
        let sample = emitter.evaluate(wavelength, 1.0).unwrap();
        let mut expected_total = 0.0;
        for index in 0..2 {
            let expected = absorption[index] * planck_si(temperatures[index], wavelength);
            expected_total += expected;
            close(sample.film_radiance_w_m2_sr_nm()[index][0], 0.5 * expected);
            close(sample.film_radiance_w_m2_sr_nm()[index][1], 0.5 * expected);
        }
        let substrate_expected =
            substrate_absorption * planck_si(substrate_temperature, wavelength);
        close(
            sample.substrate_radiance_w_m2_sr_nm().unwrap()[0],
            0.5 * substrate_expected,
        );
        expected_total += substrate_expected;
        close(sample.stokes().i, expected_total);
        assert_eq!(sample.stokes().q, 0.0);
        assert_eq!(sample.stokes().u, 0.0);
        assert_eq!(sample.stokes().v, 0.0);
        close(sample.reciprocal_absorbed_fraction()[0], entrance);
        assert_eq!(sample.reciprocal_transmitted_fraction(), [0.0; 2]);
        let average_temperature = (temperatures[0] + temperatures[1] + substrate_temperature) / 3.0;
        let averaged_source = entrance * planck_si(average_temperature, wavelength);
        assert!((sample.stokes().i / averaged_source - 1.0).abs() > 0.1);
    }
}

#[test]
fn isothermal_stacks_recover_polarized_kirchhoff_with_original_dispersive_optics() {
    let absorbing = [film(2.1, 0.27, 93.0e-9, 1), film(1.3, 0.08, 170.0e-9, 2)];
    let transparent = CoherentDielectricCoating::try_new(
        CauchyIor::try_new(1.45, 0.006, 0.0).unwrap(),
        110.0e-9,
        source(3),
    )
    .unwrap();
    let layers = [
        CoatingLayer::Absorbing(&absorbing[0]),
        CoatingLayer::Dielectric(&transparent),
        CoatingLayer::Absorbing(&absorbing[1]),
    ];
    let metal = substrate(0.8, 3.0);
    let opposite: DielectricIor = CauchyIor::try_new(1.5, 0.004, 0.0).unwrap().into();
    let temperature = 2000.0;
    let opaque = CoatingThermalEmission::try_opaque(
        &layers,
        &[temperature; 3],
        metal,
        temperature,
        ContentHash([9; 32]),
    )
    .unwrap();
    let mut nonzero_polarization = false;
    for wavelength in [380.0, 532.0, 780.0] {
        let blackbody = planck_si(temperature, wavelength);
        for cosine in [1.0, 0.66, 0.1] {
            let original = CoherentCoatingStack::try_new(&layers)
                .unwrap()
                .layer_power_balance(&metal, wavelength, 1.0, cosine)
                .unwrap();
            let sample = opaque.evaluate(wavelength, cosine).unwrap();
            let [a_s, a_p] = original.total.absorbed;
            close(sample.stokes().i, 0.5 * (a_s + a_p) * blackbody);
            if cosine == 1.0 {
                assert_eq!(sample.stokes().q, 0.0);
            } else {
                close(sample.stokes().q, 0.5 * (a_s - a_p) * blackbody);
                nonzero_polarization |= sample.stokes().q != 0.0;
            }
            assert!(sample.stokes().degree() <= 1.0);
            assert_eq!(sample.film_radiance_w_m2_sr_nm()[1], [0.0; 2]);
            for reverse in [false, true] {
                // For reverse observation this describes the declared geometry
                // with vacuum at its authored last-film side.
                let emitter = CoatingThermalEmission::try_transmitting(
                    &layers,
                    &[temperature; 3],
                    opposite,
                    reverse,
                    ContentHash([9; 32]),
                )
                .unwrap();
                let work = DielectricCoatingStack::try_mixed(&layers)
                    .unwrap()
                    .layer_power_balance(
                        wavelength,
                        reverse,
                        cosine,
                        1.0,
                        opposite.eval(wavelength).unwrap(),
                    )
                    .unwrap();
                let sample = emitter.evaluate(wavelength, cosine).unwrap();
                close(
                    sample.stokes().i,
                    0.5 * (work.absorbed[0] + work.absorbed[1]) * blackbody,
                );
                assert_eq!(sample.substrate_radiance_w_m2_sr_nm(), None);
                for pol in 0..2 {
                    close(
                        sample.reciprocal_absorbed_fraction()[pol],
                        work.absorbed[pol],
                    );
                    close(
                        sample.reciprocal_reflected_fraction()[pol]
                            + sample.reciprocal_transmitted_fraction()[pol]
                            + sample.reciprocal_absorbed_fraction()[pol],
                        1.0,
                    );
                    for (index, fraction) in work.layer_absorbed().iter().enumerate() {
                        close(
                            sample.film_radiance_w_m2_sr_nm()[index][pol],
                            0.5 * fraction[pol] * blackbody,
                        );
                    }
                }
            }
        }
    }
    assert!(nonzero_polarization);
}

#[test]
fn reverse_view_preserves_each_physical_temperature_and_weak_loss_emits_without_flux_subtraction() {
    let films = [film(2.2, 0.15, 100.0e-9, 1), film(1.4, 0.03, 270.0e-9, 2)];
    let layers = [
        CoatingLayer::Absorbing(&films[0]),
        CoatingLayer::Absorbing(&films[1]),
    ];
    let temperatures = [1000.0, 2100.0];
    let mut output = [0.0; 2];
    for (side, reverse) in [false, true].into_iter().enumerate() {
        let emitter = CoatingThermalEmission::try_transmitting(
            &layers,
            &temperatures,
            vacuum(),
            reverse,
            ContentHash([9; 32]),
        )
        .unwrap();
        let sample = emitter.evaluate(550.0, 0.63).unwrap();
        let work = DielectricCoatingStack::try_mixed(&layers)
            .unwrap()
            .layer_power_balance(550.0, reverse, 0.63, 1.0, 1.0)
            .unwrap();
        for (index, fraction) in work.layer_absorbed().iter().enumerate() {
            for pol in 0..2 {
                close(
                    sample.film_radiance_w_m2_sr_nm()[index][pol],
                    0.5 * fraction[pol] * planck_si(temperatures[index], 550.0),
                );
            }
        }
        output[side] = sample.stokes().i;
    }
    assert!((output[0] / output[1] - 1.0).abs() > 0.01);

    let k = 1.0e-20;
    let thickness = 100.0e-9;
    let weak = film(1.0, k, thickness, 7);
    let emitter = CoatingThermalEmission::try_transmitting(
        &[CoatingLayer::Absorbing(&weak)],
        &[1500.0],
        vacuum(),
        false,
        ContentHash([9; 32]),
    )
    .unwrap();
    let sample = emitter.evaluate(550.0, 1.0).unwrap();
    let expected = 4.0 * PI * k * thickness / 550.0e-9 * planck_si(1500.0, 550.0);
    assert!(sample.stokes().i > 0.0);
    close(sample.stokes().i, expected);
    assert_eq!(
        sample.reciprocal_reflected_fraction()[0] + sample.reciprocal_transmitted_fraction()[0],
        1.0
    );
    assert!(sample.reciprocal_absorbed_fraction()[0] > 0.0);
}

#[test]
fn dark_stacks_and_admission_keep_visible_vacuum_model_boundaries() {
    let dispersive = CoherentDielectricCoating::try_new(
        CauchyIor::try_new(1.5, 0.005, 0.0).unwrap(),
        100.0e-9,
        source(1),
    )
    .unwrap();
    let zero_depth = film(2.0, 0.7, 0.0, 2);
    let zero_loss = film(1.8, 0.0, 80.0e-9, 3);
    let layers = [
        CoatingLayer::Dielectric(&dispersive),
        CoatingLayer::Absorbing(&zero_depth),
        CoatingLayer::Absorbing(&zero_loss),
    ];
    let source_id = ContentHash([9; 32]);
    let emitter = CoatingThermalEmission::try_transmitting(
        &layers,
        &[1100.0, 2200.0, 3300.0],
        vacuum(),
        false,
        source_id,
    )
    .unwrap();
    for wavelength in [380.0, 600.0, 780.0] {
        let sample = emitter.evaluate(wavelength, 0.4).unwrap();
        assert_eq!(sample.stokes(), Stokes::unpolarized(0.0));
        assert_eq!(sample.film_radiance_w_m2_sr_nm(), &[[0.0; 2]; 3]);
        assert_eq!(sample.reciprocal_absorbed_fraction(), [0.0; 2]);
        assert!(sample.reciprocal_reflected_fraction()[0] > 0.0);
        assert!(sample.reciprocal_transmitted_fraction()[0] > 0.0);
    }
    let empty = CoatingThermalEmission::try_transmitting(
        &[],
        &[],
        CauchyIor::try_constant(1.5).unwrap().into(),
        false,
        source_id,
    )
    .unwrap();
    assert_eq!(
        empty.evaluate(550.0, 1.0).unwrap().stokes(),
        Stokes::unpolarized(0.0)
    );
    assert_eq!(empty.film_count(), 0);
    assert_eq!(empty.substrate_temperature_k(), None);
    for wavelength in [379.0, 781.0, 0.0, f64::NAN] {
        assert_eq!(
            empty.evaluate(wavelength, 1.0),
            Err(CoatingEmissionError::Thermal(
                ThermalSurfaceError::Wavelength
            ))
        );
    }
    for cosine in [0.0, -0.5, 1.1, f64::NAN] {
        assert_eq!(
            empty.evaluate(550.0, cosine),
            Err(CoatingEmissionError::EmissionDirection)
        );
    }
    for temperature in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(
            CoatingThermalEmission::try_transmitting(
                &layers[..1],
                &[temperature],
                vacuum(),
                false,
                source_id,
            )
            .is_err()
        );
    }
    assert!(matches!(
        CoatingThermalEmission::try_transmitting(&layers, &[1000.0], vacuum(), false, source_id,),
        Err(CoatingEmissionError::LayerTemperatures)
    ));
    assert!(
        CoatingThermalEmission::try_transmitting(&[], &[], vacuum(), false, ContentHash([0; 32]),)
            .is_err()
    );
    assert!(
        CoatingThermalEmission::try_transmitting(
            &[layers[0]; MAX_COHERENT_STACK_LAYERS + 1],
            &[1000.0; MAX_COHERENT_STACK_LAYERS + 1],
            vacuum(),
            false,
            source_id,
        )
        .is_err()
    );
    assert_ne!(
        empty.identity(),
        CoatingThermalEmission::try_transmitting(&[], &[], vacuum(), false, source_id,)
            .unwrap()
            .identity()
    );
    assert_ne!(
        emitter.identity(),
        CoatingThermalEmission::try_transmitting(
            &layers,
            &[1100.0, 2200.0, 3300.0],
            vacuum(),
            true,
            source_id,
        )
        .unwrap()
        .identity()
    );
}
