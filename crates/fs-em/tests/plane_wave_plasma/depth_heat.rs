//! Conservative cell heat from the original admitted Maxwell/carrier fields.
use super::*;
use fs_em::plane_wave::layers::{
    IonizedGasLayerExcitation, IonizedGasLayerHeating, PreparedIonizedGasLayerIllumination,
};
use fs_material::conductor::ResolvedConductorState;

fn frequency(gas: &IonizedGasConductivity) -> f64 {
    (gas.carriers()
        .iter()
        .map(|c| {
            let q = c.charge_number() as f64 * fixture::E;
            c.number_density_per_m3() * q * q / c.mass_kg()
        })
        .sum::<f64>()
        / EPSILON)
        .sqrt()
}
fn phase(k: C64, depth: f64) -> C64 {
    C64::new((k.re * depth).cos(), -(k.re * depth).sin()).scale((k.im * depth).exp())
}
fn values(h: IonizedGasLayerHeating) -> [f64; 4] {
    [
        h.total_w_m2,
        h.electron_w_m2,
        h.heavy_w_m2,
        h.exchange_correction_to_electrons_w_m2,
    ]
}

#[test]
fn coherent_cells_integrate_original_vector_loss_and_retain_sub_roundoff_heat() {
    let gas = gas();
    let transport = gas.prepare_harmonic(|| false).unwrap();
    let background = background(&gas);
    let omega = 2.0 * frequency(&gas);
    let law = IonizedGasWaveResponse::try_new(&transport, &background, omega, || false).unwrap();
    let epsilon = law.effective_permittivity_f_m();
    let ki = omega * (EPSILON * MU).sqrt();
    let d = 0.7 / ki;
    let c: f64 = 0.8;
    let kp = ki * (1.0 - c * c).sqrt();
    let kz = (epsilon.scale(omega * omega * MU) - C64::from_re(kp * kp)).sqrt();
    let input = [C64::new(0.2, 0.05), C64::new(-0.07, 0.12)];
    let films = [DielectricLayer {
        response: DielectricResponse::IonizedGas(law),
        thickness_m: d,
    }];
    for pol in [Polarization::Te, Polarization::Tm] {
        let wave = PreparedIonizedGasLayerIllumination::try_new(
            DielectricResponse::Lossless(&background),
            DielectricResponse::Lossless(&background),
            omega,
            &films,
            IonizedGasLayerExcitation {
                incidence_cosine: c,
                polarization: pol,
                incoming_tangential_e_v_m: input,
            },
            || false,
        )
        .unwrap();
        let yi = (EPSILON / MU).sqrt() * if pol == Polarization::Te { c } else { 1.0 / c };
        let y = if pol == Polarization::Te {
            kz.scale(1.0 / (omega * MU))
        } else {
            epsilon.scale(omega) / kz
        };
        let r = (C64::from_re(yi) - y) / (C64::from_re(yi) + y);
        let t = C64::from_re(2.0 * yi) / (C64::from_re(yi) + y);
        let p = phase(kz, d);
        let denominator = C64::ONE - r * r * p * p;
        let f = t * (input[0] - r * p * input[1]) / denominator;
        let b = t * (input[1] - r * p * input[0]) / denominator;
        let normal = if pol == Polarization::Tm {
            (C64::from_re(kp) / kz).norm_sq()
        } else {
            0.0
        };
        let alpha = -kz.im;
        let mut sum = 0.0;
        for ends in [0.0, 0.07 * d, 0.41 * d, 0.8 * d, d].windows(2) {
            let (start, end) = (ends[0], ends[1]);
            let width = end - start;
            // Analytic two-wave volume integral, including the TM normal E
            // and the original complex left/right source cross term.
            let diagonal = -(-2.0 * alpha * width).exp_m1() / (2.0 * alpha)
                * (f.norm_sq() * (-2.0 * alpha * start).exp()
                    + b.norm_sq() * (-2.0 * alpha * (d - end)).exp());
            let cross = 2.0 * (-alpha * d).exp() * (kz.re * width).sin() / kz.re
                * (f * b.conj() * phase(C64::from_re(-kz.re), d - start - end)).re;
            let expected = 0.5
                * law.conductivity_s_m().re
                * ((1.0 + normal) * diagonal + (1.0 - normal) * cross);
            let cell = wave
                .layer_absorbed_power_in_depth_interval_w_m2(0, start, end)
                .unwrap();
            assert!(cell > 0.0);
            near(cell, expected, 4e-10);
            sum += cell;
        }
        near(sum, wave.layer_absorbed_power_w_m2(0).unwrap(), 4e-11);
        assert_eq!(
            wave.layer_absorbed_power_in_depth_interval_w_m2(0, 0.0, d)
                .unwrap(),
            wave.layer_absorbed_power_w_m2(0).unwrap()
        );
        let tiny = 1e-20 * d;
        let heat = wave
            .layer_absorbed_power_in_depth_interval_w_m2(0, 0.0, tiny)
            .unwrap();
        let density = wave
            .field(0, 0.0)
            .unwrap()
            .fields()
            .statistics
            .absorbed_power_w_m3;
        assert!(heat > 0.0 && heat < f64::EPSILON * wave.boundary().power.incoming_w_m2[0]);
        near(heat, density * tiny, 3e-11);
    }
}

#[test]
fn passive_exit_cells_preserve_original_coulomb_heat_partitions() {
    let original = gas();
    let gas = original
        .clone()
        .with_coulomb_drag(fixture::coulomb_model(original.state()), || false)
        .unwrap();
    let transport = gas.prepare_harmonic(|| false).unwrap();
    let background = background(&gas);
    let omega = 2.0 * frequency(&gas);
    let law = IonizedGasWaveResponse::try_new(&transport, &background, omega, || false).unwrap();
    let d = 0.7 / (omega * (EPSILON * MU).sqrt());
    let films = [DielectricLayer {
        response: DielectricResponse::IonizedGas(law),
        thickness_m: d,
    }];
    for (ratio, loss) in [(0.75, 0.0), (1.4, 0.3)] {
        let state = gas.state();
        let material = background_with_properties(
            state.heavy_temperature().value(),
            state.electron_temperature().value(),
            state.pressure().value(),
            true,
            ratio * EPSILON,
            (loss > 0.0).then_some(loss * omega * EPSILON),
        );
        let conductor = (loss > 0.0)
            .then(|| ResolvedConductorState::try_from_material(material.clone()).unwrap());
        let exit = conductor.as_ref().map_or(
            DielectricResponse::Lossless(&material),
            DielectricResponse::Ohmic,
        );
        for pol in [Polarization::Te, Polarization::Tm] {
            let wave = PreparedIonizedGasLayerIllumination::try_new_one_sided(
                DielectricResponse::Lossless(&background),
                exit,
                omega,
                &films,
                IonizedGasLayerExcitation {
                    incidence_cosine: 0.5,
                    polarization: pol,
                    incoming_tangential_e_v_m: [C64::new(0.01, 0.003), C64::ZERO],
                },
                || false,
            )
            .unwrap();
            let local = wave.field(0, 0.3 * d).unwrap();
            let response = transport
                .response(omega, local.fields().electric_v_m, || false)
                .unwrap();
            let original = &response.mean_heating_w_m3;
            let fractions = [
                original.total_power_density_w_m3(),
                original.electron_power_density_w_m3(),
                original.heavy_power_density_w_m3(),
                original.thermal_exchange_correction_to_electrons_w_m3(),
            ];
            let mut summed = [0.0; 4];
            for ends in [0.0, 0.2 * d, 0.6 * d, d].windows(2) {
                let total = wave
                    .layer_absorbed_power_in_depth_interval_w_m2(0, ends[0], ends[1])
                    .unwrap();
                let heat = wave
                    .layer_carrier_heating_in_depth_interval(0, ends[0], ends[1])
                    .unwrap()
                    .unwrap();
                assert_eq!(heat.total_w_m2, total);
                assert!(heat.exchange_correction_to_electrons_w_m2 > 0.0);
                for (i, value) in values(heat).into_iter().enumerate() {
                    near(value, total * fractions[i] / fractions[0], 4e-11);
                    summed[i] += value;
                }
            }
            let full = wave.layer_carrier_heating(0).unwrap();
            for (sum, whole) in summed.into_iter().zip(values(full)) {
                near(sum, whole, 5e-11);
            }
            assert_eq!(
                wave.layer_carrier_heating_in_depth_interval(0, 0.0, d)
                    .unwrap(),
                Some(full)
            );
        }
    }
}

#[test]
fn opaque_cells_retain_skin_heat_and_validate_zero_volume_and_lossless_ranges() {
    let gas = gas();
    let transport = gas.prepare_harmonic(|| false).unwrap();
    let background = background(&gas);
    let omega = 0.6 * frequency(&gas);
    let law = IonizedGasWaveResponse::try_new(&transport, &background, omega, || false).unwrap();
    let kz = law
        .effective_permittivity_f_m()
        .scale(omega * omega * MU)
        .sqrt();
    let attenuation = -2.0 * kz.im;
    let d = 4000.0 / attenuation;
    let films = [
        DielectricLayer {
            response: DielectricResponse::IonizedGas(law),
            thickness_m: d,
        },
        DielectricLayer {
            response: DielectricResponse::IonizedGas(law),
            thickness_m: 0.0,
        },
        DielectricLayer {
            response: DielectricResponse::Lossless(&background),
            thickness_m: d,
        },
    ];
    let wave = PreparedIonizedGasLayerIllumination::try_new_one_sided(
        DielectricResponse::Lossless(&background),
        DielectricResponse::Lossless(&background),
        omega,
        &films,
        IonizedGasLayerExcitation {
            incidence_cosine: 1.0,
            polarization: Polarization::Te,
            incoming_tangential_e_v_m: [C64::from_re(0.1), C64::ZERO],
        },
        || false,
    )
    .unwrap();
    let density = wave
        .field(0, 0.0)
        .unwrap()
        .fields()
        .statistics
        .absorbed_power_w_m3;
    let mut sum = 0.0;
    for ends in [
        0.0,
        0.2 / attenuation,
        1.0 / attenuation,
        5.0 / attenuation,
        50.0 / attenuation,
        d,
    ]
    .windows(2)
    {
        let expected = density / attenuation
            * (-attenuation * ends[0]).exp()
            * -(-attenuation * (ends[1] - ends[0])).exp_m1();
        let cell = wave
            .layer_absorbed_power_in_depth_interval_w_m2(0, ends[0], ends[1])
            .unwrap();
        near(cell, expected, 4e-11);
        sum += cell;
    }
    near(sum, wave.layer_absorbed_power_w_m2(0).unwrap(), 4e-11);
    assert_eq!(
        wave.field(0, 0.5 * d).unwrap().fields().electric_v_m,
        [C64::ZERO; 3]
    );
    assert_eq!(
        wave.layer_carrier_heating_in_depth_interval(1, 0.0, 0.0)
            .unwrap(),
        Some(IonizedGasLayerHeating::default())
    );
    assert_eq!(
        wave.layer_carrier_heating_in_depth_interval(2, 0.0, d)
            .unwrap(),
        None
    );
    assert_eq!(
        wave.layer_absorbed_power_in_depth_interval_w_m2(2, 0.0, d)
            .unwrap(),
        0.0
    );
    assert_eq!(
        wave.layer_carrier_heating_in_depth_interval(0, 0.3 * d, 0.3 * d)
            .unwrap(),
        Some(IonizedGasLayerHeating::default())
    );
    for index in 0..3 {
        for (start, end) in [
            (-1.0, 0.0),
            (0.0, d * 1.01),
            (1.0, 0.0),
            (f64::NAN, 0.0),
            (0.0, f64::INFINITY),
        ] {
            assert!(
                wave.layer_absorbed_power_in_depth_interval_w_m2(index, start, end)
                    .is_err()
            );
            assert!(
                wave.layer_carrier_heating_in_depth_interval(index, start, end)
                    .is_err()
            );
        }
    }
    assert!(
        wave.layer_absorbed_power_in_depth_interval_w_m2(3, 0.0, 0.0)
            .is_err()
    );
    assert!(
        wave.layer_carrier_heating_in_depth_interval(3, 0.0, 0.0)
            .is_err()
    );
}
