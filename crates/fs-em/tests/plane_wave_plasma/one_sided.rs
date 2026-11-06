//! Original finite gas fields with passive absorbing, evanescent and critical exits.
use super::*;
use fs_em::plane_wave::layers::{IonizedGasLayerExcitation, PreparedIonizedGasLayerIllumination};
use fs_material::conductor::ResolvedConductorState;

fn frequency(gas: &IonizedGasConductivity) -> f64 {
    (gas.carriers()
        .iter()
        .map(|c| {
            let charge = c.charge_number() as f64 * fixture::E;
            c.number_density_per_m3() * charge * charge / c.mass_kg()
        })
        .sum::<f64>()
        / EPSILON)
        .sqrt()
}

fn drive(c: f64, p: Polarization, amplitude: C64) -> IonizedGasLayerExcitation {
    IonizedGasLayerExcitation {
        incidence_cosine: c,
        polarization: p,
        incoming_tangential_e_v_m: [amplitude, C64::ZERO],
    }
}

fn medium(
    gas: &IonizedGasConductivity,
    ratio: f64,
    sigma: Option<f64>,
) -> ResolvedMaterialStatePoint {
    let s = gas.state();
    background_with_properties(
        s.heavy_temperature().value(),
        s.electron_temperature().value(),
        s.pressure().value(),
        true,
        ratio * EPSILON,
        sigma,
    )
}

fn phase(k: C64, d: f64) -> C64 {
    C64::new((k.re * d).cos(), -(k.re * d).sin()).scale((k.im * d).exp())
}

#[test]
fn one_sided_gas_films_retain_original_two_port_fields_currents_and_heat() {
    let gas = gas();
    let transport = gas.prepare_harmonic(|| false).unwrap();
    let background = background(&gas);
    for ratio in [0.6, 2.0] {
        let omega = ratio * frequency(&gas);
        let ki = omega * (EPSILON * MU).sqrt();
        let law =
            IonizedGasWaveResponse::try_new(&transport, &background, omega, || false).unwrap();
        let films = [
            DielectricLayer {
                response: DielectricResponse::IonizedGas(law),
                thickness_m: 0.0,
            },
            DielectricLayer {
                response: DielectricResponse::IonizedGas(law),
                thickness_m: 0.7 / ki,
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
        for pol in [Polarization::Te, Polarization::Tm] {
            let excitation = drive(0.8, pol, C64::new(0.2, -0.07));
            let one = PreparedIonizedGasLayerIllumination::try_new_one_sided(
                DielectricResponse::Lossless(&background),
                DielectricResponse::Lossless(&background),
                omega,
                &films,
                excitation,
                || false,
            )
            .unwrap();
            let two = PreparedIonizedGasLayerIllumination::try_new(
                DielectricResponse::Lossless(&background),
                DielectricResponse::Lossless(&background),
                omega,
                &films,
                excitation,
                || false,
            )
            .unwrap();
            for face in 0..2 {
                complex_near(
                    one.boundary().outgoing_tangential_e_v_m[face],
                    two.boundary().outgoing_tangential_e_v_m[face],
                    3e-11,
                );
                complex_near(
                    one.boundary().boundary_conjugate_h_a_m[face],
                    two.boundary().boundary_conjugate_h_a_m[face],
                    3e-11,
                );
                near(
                    one.boundary().power.outgoing_w_m2[face],
                    two.boundary().power.outgoing_w_m2[face],
                    3e-11,
                );
            }
            for (i, film) in films.iter().enumerate() {
                near(
                    one.layer_absorbed_power_w_m2(i).unwrap(),
                    two.layer_absorbed_power_w_m2(i).unwrap(),
                    4e-11,
                );
                if film.thickness_m == 0.0 {
                    assert!(one.field(i, 0.0).is_err());
                    assert_eq!(one.layer_carrier_heating(i).unwrap().total_w_m2, 0.0);
                    continue;
                }
                if let Some(heat) = one.layer_carrier_heating(i) {
                    let old = two.layer_carrier_heating(i).unwrap();
                    near(heat.electron_w_m2, old.electron_w_m2, 4e-11);
                    near(heat.heavy_w_m2, old.heavy_w_m2, 4e-11);
                }
                for fraction in [0.0, 0.37, 1.0] {
                    let a = one.field(i, fraction * film.thickness_m).unwrap();
                    let b = two.field(i, fraction * film.thickness_m).unwrap();
                    for axis in 0..3 {
                        complex_near(
                            a.fields().electric_v_m[axis],
                            b.fields().electric_v_m[axis],
                            4e-11,
                        );
                        complex_near(
                            a.fields().magnetic_a_m[axis],
                            b.fields().magnetic_a_m[axis],
                            4e-11,
                        );
                        complex_near(
                            a.current_density_a_m2[axis],
                            b.current_density_a_m2[axis],
                            4e-11,
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn passive_substrates_match_dimensional_maxwell_fields_and_absorb_original_exit_power() {
    let gas = gas();
    let transport = gas.prepare_harmonic(|| false).unwrap();
    let background = background(&gas);
    let omega = 2.0 * frequency(&gas);
    let law = IonizedGasWaveResponse::try_new(&transport, &background, omega, || false).unwrap();
    let epsilon = law.effective_permittivity_f_m();
    let ki = omega * (EPSILON * MU).sqrt();
    let c: f64 = 0.5;
    let kp = ki * (1.0 - c * c).sqrt();
    let kz = (epsilon.scale(omega * omega * MU) - C64::from_re(kp * kp)).sqrt();
    let thickness = 0.7 / ki;
    let films = [DielectricLayer {
        response: DielectricResponse::IonizedGas(law),
        thickness_m: thickness,
    }];
    let amplitude = C64::new(0.2, 0.07);
    for (ratio, loss) in [(1.4, 0.3), (0.5, 0.0), (0.75, 0.0)] {
        let sigma = loss * omega * EPSILON;
        let material = medium(&gas, ratio, (sigma > 0.0).then_some(sigma));
        let conductor = (sigma > 0.0)
            .then(|| ResolvedConductorState::try_from_material(material.clone()).unwrap());
        let exit = conductor.as_ref().map_or(
            DielectricResponse::Lossless(&material),
            DielectricResponse::Ohmic,
        );
        let exit_epsilon = C64::new(ratio * EPSILON, -sigma / omega);
        // The authored ratio 3/4 and cosine 1/2 give kz=0 exactly. Otherwise
        // use independent dimensional dispersion and the outgoing passive root.
        let mut ke = if ratio == 0.75 {
            C64::ZERO
        } else {
            (exit_epsilon.scale(omega * omega * MU) - C64::from_re(kp * kp)).sqrt()
        };
        if ke.im > 0.0 {
            ke = -ke;
        }
        for pol in [Polarization::Te, Polarization::Tm] {
            let wave = PreparedIonizedGasLayerIllumination::try_new_one_sided(
                DielectricResponse::Lossless(&background),
                exit,
                omega,
                &films,
                drive(c, pol, amplitude),
                || false,
            )
            .unwrap();
            let yi = (EPSILON / MU).sqrt() * if pol == Polarization::Te { c } else { 1.0 / c };
            let y = if pol == Polarization::Te {
                kz.scale(1.0 / (omega * MU))
            } else {
                epsilon.scale(omega) / kz
            };
            // Match physical E/H directly. Choosing an H-normalized TM exit
            // remains finite at criticality, where tangential E alone is zero.
            let [ee, hh] = if pol == Polarization::Te {
                [C64::ONE, ke.scale(1.0 / (omega * MU))]
            } else {
                [ke / exit_epsilon.scale(omega), C64::ONE]
            };
            let r = (C64::from_re(yi) - y) / (C64::from_re(yi) + y);
            let tail = if ke == C64::ZERO {
                if pol == Polarization::Te {
                    C64::ONE
                } else {
                    -C64::ONE
                }
            } else {
                (y * ee - hh) / (y * ee + hh)
            };
            let p = phase(kz, thickness);
            let f = C64::from_re(2.0 * yi) / (C64::from_re(yi) + y) * amplitude
                / (C64::ONE + r * tail * p * p);
            let b = tail * f * p;
            let u = y.scale(2.0) * f * p / (y * ee + hh);
            let et = ee * u;
            let ht = hh * u;
            complex_near(
                wave.boundary().outgoing_tangential_e_v_m[0],
                f + b * p - amplitude,
                5e-11,
            );
            complex_near(wave.boundary().boundary_tangential_e_v_m[1], et, 5e-11);
            complex_near(wave.boundary().boundary_conjugate_h_a_m[1], ht, 5e-11);
            near(
                wave.boundary().power.incoming_w_m2[0],
                wave.boundary().power.outgoing_w_m2.iter().sum::<f64>()
                    + wave.boundary().power.absorbed_w_m2,
                5e-11,
            );
            for fraction in [0.0, 0.37, 1.0] {
                let depth = fraction * thickness;
                let field = wave.field(0, depth).unwrap();
                let ft = f * phase(kz, depth);
                let bt = b * phase(kz, thickness - depth);
                let (electric, magnetic) = (ft + bt, y * (ft - bt));
                let (e, h) = if pol == Polarization::Te {
                    (
                        [C64::ZERO, electric, C64::ZERO],
                        [-magnetic, C64::ZERO, electric.scale(kp / (omega * MU))],
                    )
                } else {
                    (
                        [electric, C64::ZERO, -magnetic.scale(kp / omega) / epsilon],
                        [C64::ZERO, magnetic, C64::ZERO],
                    )
                };
                for axis in 0..3 {
                    complex_near(field.fields().electric_v_m[axis], e[axis], 6e-11);
                    complex_near(field.fields().magnetic_a_m[axis], h[axis], 6e-11);
                }
            }
            let (exit_e, exit_h) = if pol == Polarization::Te {
                (
                    [C64::ZERO, et, C64::ZERO],
                    [-ht, C64::ZERO, et.scale(kp / (omega * MU))],
                )
            } else {
                (
                    [et, C64::ZERO, -ht.scale(kp / omega) / exit_epsilon],
                    [C64::ZERO, ht, C64::ZERO],
                )
            };
            for depth in [0.0, 0.35 / ki, 1e-14 / ki] {
                let point = wave.exit_field(depth).unwrap();
                for axis in 0..3 {
                    complex_near(
                        point.electric_v_m[axis],
                        exit_e[axis] * phase(ke, depth),
                        6e-11,
                    );
                    complex_near(
                        point.magnetic_a_m[axis],
                        exit_h[axis] * phase(ke, depth),
                        6e-11,
                    );
                }
                let heat = wave.absorbed_exit_power_to_depth_w_m2(depth).unwrap();
                if sigma > 0.0 {
                    let density = 0.5 * sigma * exit_e.iter().map(|e| e.norm_sq()).sum::<f64>();
                    let expected = density * (2.0 * ke.im * depth).exp_m1() / (2.0 * ke.im);
                    near(heat, expected, 6e-11);
                    near(
                        point.statistics.absorbed_power_w_m3,
                        density * (2.0 * ke.im * depth).exp(),
                        6e-11,
                    );
                } else {
                    assert_eq!(heat, 0.0);
                }
            }
            if sigma == 0.0 {
                assert_eq!(wave.boundary().power.outgoing_w_m2[1], 0.0);
            }
            if ratio == 0.75 && pol == Polarization::Tm {
                let critical = wave.exit_field(1.0 / ki).unwrap();
                assert_eq!(critical.electric_v_m[0], C64::ZERO);
                assert!(critical.electric_v_m[2].abs() > 0.0);
                assert!(critical.magnetic_a_m[1].abs() > 0.0);
            }
            if sigma > 0.0 {
                let depth = 1000.0 / -ke.im;
                assert_eq!(wave.exit_field(depth).unwrap().electric_v_m, [C64::ZERO; 3]);
                near(
                    wave.absorbed_exit_power_to_depth_w_m2(depth).unwrap(),
                    wave.boundary().power.outgoing_w_m2[1],
                    6e-11,
                );
            }
        }
    }
}

#[test]
fn one_sided_whole_film_drift_admission_and_passive_boundary_domain_remain_atomic() {
    let gas = gas();
    let transport = gas.prepare_harmonic(|| false).unwrap();
    let background = background(&gas);
    let critical = medium(&gas, 0.75, None);
    let omega = 2.0 * frequency(&gas);
    let ki = omega * (EPSILON * MU).sqrt();
    let law = IonizedGasWaveResponse::try_new(&transport, &background, omega, || false).unwrap();
    let films = [DielectricLayer {
        response: DielectricResponse::IonizedGas(law),
        thickness_m: 1e-6 / ki,
    }];
    let prepare = |excitation, cancel: &mut dyn FnMut() -> bool| {
        PreparedIonizedGasLayerIllumination::try_new_one_sided(
            DielectricResponse::Lossless(&background),
            DielectricResponse::Lossless(&critical),
            omega,
            &films,
            excitation,
            cancel,
        )
    };
    let unit = transport
        .response(omega, [C64::ONE, C64::ZERO, C64::ZERO], || false)
        .unwrap();
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
    let amplitude = 0.6 * max_field;
    assert!(
        transport
            .response(
                omega,
                [C64::from_re(amplitude), C64::ZERO, C64::ZERO],
                || false
            )
            .is_ok()
    );
    // A thin film against the TE critical exit has nearly twice its incoming
    // field. The original gas is admitted against that standing wave envelope.
    assert!(
        prepare(
            drive(0.5, Polarization::Te, C64::from_re(amplitude)),
            &mut || false
        )
        .is_err()
    );
    assert!(
        prepare(
            drive(0.5, Polarization::Te, C64::from_re(0.2 * max_field)),
            &mut || false
        )
        .is_ok()
    );
    let zero = prepare(drive(0.5, Polarization::Te, C64::ZERO), &mut || false).unwrap();
    assert_eq!(
        zero.field(0, 0.0).unwrap().mean_carrier_heating_w_m3(),
        Some(0.0)
    );
    assert_eq!(zero.exit_field(0.0).unwrap().electric_v_m, [C64::ZERO; 3]);
    for depth in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(zero.exit_field(depth).is_err());
        assert!(zero.absorbed_exit_power_to_depth_w_m2(depth).is_err());
    }
    let mut invalid = drive(0.5, Polarization::Te, C64::ONE);
    invalid.incoming_tangential_e_v_m[1] = C64::ONE;
    assert!(prepare(invalid, &mut || false).is_err());
    assert!(matches!(
        prepare(drive(0.5, Polarization::Te, C64::ONE), &mut || true),
        Err(ReflectionError::IonizedGas(IonizedGasError::Cancelled))
    ));
    assert!(
        PlanarDielectricLayers::try_new(
            &background,
            DielectricResponse::Lossless(&critical),
            omega,
            &films
        )
        .is_err()
    );
    // The existing coherent two-port owner still refuses a nonpropagating port.
    assert!(
        PreparedIonizedGasLayerIllumination::try_new(
            DielectricResponse::Lossless(&background),
            DielectricResponse::Lossless(&critical),
            omega,
            &films,
            drive(0.5, Polarization::Te, C64::ONE),
            || false
        )
        .is_err()
    );
}
