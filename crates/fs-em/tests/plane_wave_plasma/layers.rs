//! Original carrier physics in finite coherent refracting films.
use super::*;
use fs_em::plane_wave::layers::{IonizedGasLayerExcitation, PreparedIonizedGasLayerIllumination};

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

fn excitation(c: f64, p: Polarization, input: [C64; 2]) -> IonizedGasLayerExcitation {
    IonizedGasLayerExcitation {
        incidence_cosine: c,
        polarization: p,
        incoming_tangential_e_v_m: input,
    }
}

fn phase(k: C64, d: f64) -> C64 {
    C64::new((k.re * d).cos(), -(k.re * d).sin()).scale((k.im * d).exp())
}

fn prepare<'a>(
    background: &'a ResolvedMaterialStatePoint,
    omega: f64,
    films: &[DielectricLayer<'a>],
    drive: IonizedGasLayerExcitation,
) -> Result<PreparedIonizedGasLayerIllumination<'a>, ReflectionError> {
    PreparedIonizedGasLayerIllumination::try_new(
        DielectricResponse::Lossless(background),
        DielectricResponse::Lossless(background),
        omega,
        films,
        drive,
        || false,
    )
}

#[test]
fn plasma_films_match_independent_bidirectional_fabry_perot_and_volume_heat() {
    let gas = gas();
    let transport = gas.prepare_harmonic(|| false).unwrap();
    let background = background(&gas);
    let input = [C64::new(0.2, 0.05), C64::new(-0.07, 0.12)];
    for ratio in [0.6, 1.0, 2.0] {
        let omega = ratio * frequency(&gas);
        let law =
            IonizedGasWaveResponse::try_new(&transport, &background, omega, || false).unwrap();
        let epsilon = law.effective_permittivity_f_m();
        let ki = omega * (EPSILON * MU).sqrt();
        let d = 0.7 / ki;
        for polarization in [Polarization::Te, Polarization::Tm] {
            let c: f64 = 0.8;
            let kp = ki * (1.0 - c * c).sqrt();
            let kz = (epsilon.scale(omega * omega * MU) - C64::from_re(kp * kp)).sqrt();
            let y = match polarization {
                Polarization::Te => kz.scale(1.0 / (omega * MU)),
                Polarization::Tm => epsilon.scale(omega) / kz,
            };
            let yi = (EPSILON / MU).sqrt()
                * match polarization {
                    Polarization::Te => c,
                    Polarization::Tm => 1.0 / c,
                };
            let r = (C64::from_re(yi) - y) / (C64::from_re(yi) + y);
            let t = C64::from_re(2.0 * yi) / (C64::from_re(yi) + y);
            let p = phase(kz, d);
            let denom = C64::ONE - r * r * p * p;
            let f = t * (input[0] - r * p * input[1]) / denom;
            let b = t * (input[1] - r * p * input[0]) / denom;
            let reflected = r * (C64::ONE - p * p) / denom;
            let transmitted = (C64::ONE - r * r) * p / denom;
            let films = [DielectricLayer {
                response: DielectricResponse::IonizedGas(law),
                thickness_m: d,
            }];
            let wave = PreparedIonizedGasLayerIllumination::try_new(
                DielectricResponse::Lossless(&background),
                DielectricResponse::Lossless(&background),
                omega,
                &films,
                excitation(c, polarization, input),
                || false,
            )
            .unwrap();
            let boundary = wave.boundary();
            complex_near(
                boundary.outgoing_tangential_e_v_m[0],
                reflected * input[0] + transmitted * input[1],
                2e-11,
            );
            complex_near(
                boundary.outgoing_tangential_e_v_m[1],
                transmitted * input[0] + reflected * input[1],
                2e-11,
            );
            near(
                boundary.power.incoming_w_m2.iter().sum::<f64>(),
                boundary.power.outgoing_w_m2.iter().sum::<f64>() + boundary.power.absorbed_w_m2,
                2e-11,
            );
            for fraction in [0.0, 0.17, 0.5, 0.83, 1.0] {
                let depth = d * fraction;
                let forward = f * phase(kz, depth);
                let backward = b * phase(kz, d - depth);
                let et = forward + backward;
                let ht = y * (forward - backward);
                let (e, h) = match polarization {
                    Polarization::Te => (
                        [C64::ZERO, et, C64::ZERO],
                        [-ht, C64::ZERO, et.scale(kp / (omega * MU))],
                    ),
                    Polarization::Tm => (
                        [et, C64::ZERO, -ht.scale(kp / omega) / epsilon],
                        [C64::ZERO, ht, C64::ZERO],
                    ),
                };
                let point = wave.field(0, depth).unwrap();
                for axis in 0..3 {
                    complex_near(point.fields().electric_v_m[axis], e[axis], 3e-11);
                    complex_near(point.fields().magnetic_a_m[axis], h[axis], 3e-11);
                }
                let original = transport
                    .response(omega, point.fields().electric_v_m, || false)
                    .unwrap();
                for axis in 0..3 {
                    complex_near(
                        point.current_density_a_m2[axis],
                        original.current_density_a_m2[axis],
                        3e-13,
                    );
                }
                for carrier in 0..gas.carriers().len() {
                    let actual = point.carrier_velocity_m_s(carrier).unwrap();
                    for axis in 0..3 {
                        complex_near(
                            actual[axis],
                            original.carrier_velocity_m_s[carrier][axis],
                            4e-12,
                        );
                    }
                }
                near(
                    point.mean_carrier_heating_w_m3().unwrap(),
                    original.mean_heating_w_m3.total_power_density_w_m3(),
                    4e-12,
                );
                near(
                    point.fields().statistics.absorbed_power_w_m3,
                    point.mean_carrier_heating_w_m3().unwrap(),
                    4e-12,
                );
                near(
                    point.mean_drift_kinetic_energy_j_m3().unwrap(),
                    original.mean_drift_kinetic_energy_j_m3,
                    4e-12,
                );
                near(
                    point.reactive_carrier_power_var_m3().unwrap(),
                    original.reactive_power_var_m3,
                    4e-12,
                );
                let storage = 0.25
                    * (EPSILON * e.iter().map(|z| z.norm_sq()).sum::<f64>()
                        + MU * h.iter().map(|z| z.norm_sq()).sum::<f64>());
                near(
                    point.mean_background_field_energy_j_m3().unwrap(),
                    storage,
                    5e-11,
                );
                let norm = e.iter().map(|z| z.norm_sq()).sum::<f64>().sqrt();
                assert!(norm <= wave.layer_electric_envelope_v_m(0).unwrap() * (1.0 + 3e-11));
            }
            // Closed-form integral of the independent two-wave E field.
            // Opposing waves share constant cross-term attenuation exp(-alpha*d).
            let alpha = -kz.im;
            let diagonal = -(-2.0 * alpha * d).exp_m1() / (2.0 * alpha);
            let cross = 2.0 * (-alpha * d).exp() * (kz.re * d).sin() / kz.re * (f * b.conj()).re;
            let normal = if polarization == Polarization::Tm {
                (C64::from_re(kp) / kz).norm_sq()
            } else {
                0.0
            };
            let integral =
                (1.0 + normal) * (f.norm_sq() + b.norm_sq()) * diagonal + (1.0 - normal) * cross;
            let expected_heat = 0.5 * law.conductivity_s_m().re * integral;
            near(
                wave.layer_absorbed_power_w_m2(0).unwrap(),
                expected_heat,
                4e-10,
            );
            near(
                wave.layer_carrier_heating(0).unwrap().total_w_m2,
                expected_heat,
                4e-10,
            );
        }
    }
}

#[test]
fn multilayer_coulomb_fields_retain_original_heat_partition_and_zero_thickness_identity() {
    let original = gas();
    let gas = original
        .clone()
        .with_coulomb_drag(fixture::coulomb_model(original.state()), || false)
        .unwrap();
    let transport = gas.prepare_harmonic(|| false).unwrap();
    let background = background(&gas);
    let omega = 2.0 * frequency(&gas);
    let ki = omega * (EPSILON * MU).sqrt();
    let law = IonizedGasWaveResponse::try_new(&transport, &background, omega, || false).unwrap();
    let films = [
        DielectricLayer {
            response: DielectricResponse::IonizedGas(law),
            thickness_m: 0.0,
        },
        DielectricLayer {
            response: DielectricResponse::IonizedGas(law),
            thickness_m: 0.8 / ki,
        },
        DielectricLayer {
            response: DielectricResponse::Lossless(&background),
            thickness_m: 0.3 / ki,
        },
        DielectricLayer {
            response: DielectricResponse::IonizedGas(law),
            thickness_m: 0.4 / ki,
        },
    ];
    let drive = excitation(
        0.7,
        Polarization::Tm,
        [C64::new(0.2, 0.03), C64::new(0.04, -0.08)],
    );
    let wave = prepare(&background, omega, &films, drive).unwrap();
    let without_identity = prepare(&background, omega, &films[1..], drive).unwrap();
    assert_eq!(wave.layer_count(), 4);
    assert_eq!(wave.omega_rad_s(), omega);
    for port in 0..2 {
        complex_near(
            wave.boundary().outgoing_tangential_e_v_m[port],
            without_identity.boundary().outgoing_tangential_e_v_m[port],
            4e-12,
        );
    }
    assert_eq!(wave.layer_absorbed_power_w_m2(0), Some(0.0));
    assert_eq!(wave.layer_carrier_heating(0).unwrap().total_w_m2, 0.0);
    assert!(wave.field(0, 0.0).is_err());
    assert!(wave.layer_electric_envelope_v_m(0).is_none());
    assert!(wave.layer_carrier_heating(2).is_none());
    assert_eq!(
        wave.field(2, 0.1 / ki).unwrap().current_density_a_m2,
        [C64::ZERO; 3]
    );
    for index in [1, 3] {
        let point = wave.field(index, 0.4 * films[index].thickness_m).unwrap();
        let original = transport
            .response(omega, point.fields().electric_v_m, || false)
            .unwrap();
        near(
            point.mean_electron_heating_w_m3().unwrap(),
            original.mean_heating_w_m3.electron_power_density_w_m3(),
            3e-11,
        );
        near(
            point.mean_heavy_heating_w_m3().unwrap(),
            original.mean_heating_w_m3.heavy_power_density_w_m3(),
            3e-11,
        );
        near(
            point.mean_exchange_correction_to_electrons_w_m3().unwrap(),
            original
                .mean_heating_w_m3
                .thermal_exchange_correction_to_electrons_w_m3(),
            3e-11,
        );
        let heat = wave.layer_carrier_heating(index).unwrap();
        near(heat.electron_w_m2 + heat.heavy_w_m2, heat.total_w_m2, 3e-12);
        near(
            heat.electron_w_m2 / heat.total_w_m2,
            point.mean_electron_heating_w_m3().unwrap()
                / point.mean_carrier_heating_w_m3().unwrap(),
            3e-11,
        );
        assert!(heat.exchange_correction_to_electrons_w_m2 > 0.0);
    }
    let sum = (0..4)
        .map(|i| wave.layer_absorbed_power_w_m2(i).unwrap())
        .sum::<f64>();
    near(sum, wave.boundary().power.absorbed_w_m2, 3e-12);
    assert!(wave.field(4, 0.0).is_err());
    assert!(wave.layer_source(4).is_none());
}

#[test]
fn simultaneous_fields_are_admitted_throughout_the_original_gas_and_cancel_atomically() {
    let gas = gas();
    let transport = gas.prepare_harmonic(|| false).unwrap();
    let background = background(&gas);
    let omega = 2.0 * frequency(&gas);
    let ki = omega * (EPSILON * MU).sqrt();
    let law = IonizedGasWaveResponse::try_new(&transport, &background, omega, || false).unwrap();
    let unit = transport
        .response(omega, [C64::ONE, C64::ZERO, C64::ZERO], || false)
        .unwrap();
    assert!(gas.coulomb_pairs().is_empty());
    let max_field = gas
        .carriers()
        .iter()
        .enumerate()
        .map(|(i, c)| {
            (1.380649e-23 * c.temperature_k() / c.mass_kg()).sqrt()
                * gas.limits().maximum_drift_thermal_ratio
                / unit.carrier_velocity_m_s[i][0].abs()
        })
        .fold(f64::INFINITY, f64::min);
    let amplitude = 0.55 * max_field;
    assert!(
        transport
            .response(
                omega,
                [C64::from_re(amplitude), C64::ZERO, C64::ZERO],
                || false
            )
            .is_ok()
    );
    let films = [DielectricLayer {
        response: DielectricResponse::IonizedGas(law),
        thickness_m: 1e-6 / ki,
    }];
    let prepare = |input, frequency, cancel: &mut dyn FnMut() -> bool| {
        PreparedIonizedGasLayerIllumination::try_new(
            DielectricResponse::Lossless(&background),
            DielectricResponse::Lossless(&background),
            frequency,
            &films,
            excitation(1.0, Polarization::Te, input),
            cancel,
        )
    };
    assert!(prepare([C64::from_re(amplitude); 2], omega, &mut || false).is_err());
    assert!(prepare([C64::from_re(amplitude), C64::ZERO], omega, &mut || false).is_ok());
    let silence = prepare([C64::ZERO; 2], omega, &mut || false).unwrap();
    assert_eq!(silence.boundary().power.absorbed_w_m2, 0.0);
    let point = silence.field(0, 0.5 * films[0].thickness_m).unwrap();
    assert_eq!(point.fields().electric_v_m, [C64::ZERO; 3]);
    assert_eq!(point.mean_carrier_heating_w_m3(), Some(0.0));
    assert_eq!(point.carrier_velocity_m_s(0), Some([C64::ZERO; 3]));
    assert!(point.carrier_velocity_m_s(gas.carriers().len()).is_none());
    for depth in [-1.0, f64::NAN, f64::INFINITY, 2.0 * films[0].thickness_m] {
        assert!(silence.field(0, depth).is_err());
    }
    assert!(prepare([C64::ONE, C64::ZERO], omega * 1.01, &mut || false).is_err());
    assert!(prepare([C64::new(f64::NAN, 0.0), C64::ZERO], omega, &mut || false).is_err());
    let mut polls = 0;
    assert!(matches!(
        prepare([C64::ONE, C64::ZERO], omega, &mut || {
            polls += 1;
            polls >= 4
        }),
        Err(ReflectionError::IonizedGas(IonizedGasError::Cancelled))
    ));
    // Generic normalized films remain unable to apply arbitrary unadmitted fields.
    assert!(
        PlanarDielectricLayers::try_new(
            &background,
            DielectricResponse::Lossless(&background),
            omega,
            &films
        )
        .is_err()
    );
    assert!(
        PreparedIonizedGasLayerIllumination::try_new(
            DielectricResponse::Lossless(&background),
            DielectricResponse::IonizedGas(law),
            omega,
            &films,
            excitation(1.0, Polarization::Te, [C64::ONE, C64::ZERO]),
            || false
        )
        .is_err()
    );
}
