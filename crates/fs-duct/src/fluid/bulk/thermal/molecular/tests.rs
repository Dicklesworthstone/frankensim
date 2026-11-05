use super::*;
use crate::DuctError;
use crate::compose::{ComposeOptions, redheffer_star};
use ThermalEndCondition::{Adiabatic, Isothermal};
use fs_material::gas::{ConductivityModel, GasSpec, GasState};
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw, RelaxingMedium,
    SpeciesRelaxation, VibrationalMode,
};
use fs_math::{c64::C64, det};
use fs_thermochem::SpeciesId;

const RADIUS: f64 = 0.02;
const LIMIT: f64 = 1e12;

mod regional;

// Authored transport and reservoir data isolate the coupled equations; these
// are not measured air properties. Strong transport makes BOTH spatial modes
// resolvable in a short independent integration without a growing long-span
// transfer oracle. Gas construction, calorics and admission are all native.
fn medium(count: usize, time_scale: f64) -> RelaxationState {
    let gas = GasState::try_new(
        &GasSpec {
            sutherland_beta: 1.458e-1,
            conductivity: ConductivityModel::Eucken,
            ..GasSpec::dry_air_ussa1976()
        },
        293.15,
        101325.0,
    )
    .unwrap();
    state_with_gas(gas, count, time_scale, "authored")
}

fn state_with_gas(gas: GasState, count: usize, time_scale: f64, name: &str) -> RelaxationState {
    let id = SpeciesId::new(name).unwrap();
    let processes = (0..count)
        .map(|i| RelaxationProcess {
            reservoir: InternalReservoir::Vibrational(vec![VibrationalMode {
                theta_k: 800.0 + 700.0 * i as f64,
                degeneracy: 1,
            }]),
            rate: RelaxationRateLaw::PartnerIndependent {
                p_tau_ref_pa_s: gas.pressure * [0.0007, 0.003][i] * time_scale,
                t_ref_k: gas.temperature,
                exponent: 0.0,
                t_min_k: 200.0,
                t_max_k: 500.0,
            },
        })
        .collect();
    RelaxationAbsorption::new(vec![SpeciesRelaxation {
        species: id.clone(),
        processes,
    }])
    .unwrap()
    .evaluate(&RelaxingMedium::new(gas, vec![(id, 1.0)]).unwrap())
    .unwrap()
}

fn prepare(
    state: &RelaxationState,
    length: f64,
    omega: f64,
    ends: [ThermalEndCondition; 2],
) -> BulkThermalCylinder {
    BulkThermalCylinder::from_relaxation(state, RADIUS, length, omega, ends, LIMIT).unwrap()
}

fn near(a: C64, b: C64, scale: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * scale,
        "{a:?} != {b:?}, scale {scale}"
    );
}

#[test]
fn g3_empty_reservoirs_are_the_original_complete_nsf_cylinder() {
    let state = medium(0, 1.0);
    let material = BulkThermoviscousEquilibrium::from_material(state.gas(), 0.0).unwrap();
    for omega in [0.1, 100.0, 10000.0] {
        for ends in [
            [Isothermal, Isothermal],
            [Isothermal, Adiabatic],
            [Adiabatic, Adiabatic],
        ] {
            let original =
                BulkThermalCylinder::new(&material, RADIUS, 0.04, omega, ends, LIMIT).unwrap();
            let molecular = prepare(&state, 0.04, omega, ends);
            assert_eq!(molecular.scattering(), original.scattering());
            assert_eq!(molecular.references(), original.references());
            assert_eq!(
                molecular.traction_scattering(),
                original.traction_scattering()
            );
            for load in [
                None,
                Some(C64::new(1.2, 0.3).scale(original.references()[0])),
            ] {
                let a = molecular.response(load, ComposeOptions::default()).unwrap();
                let b = original.response(load, ComposeOptions::default()).unwrap();
                assert_eq!(
                    a.mechanical.inlet_traction_impedance_pa_s_m3,
                    b.mechanical.inlet_traction_impedance_pa_s_m3
                );
                assert_eq!(
                    a.mechanical.inlet_pressure_over_flow_pa_s_m3,
                    b.mechanical.inlet_pressure_over_flow_pa_s_m3
                );
                assert_eq!(
                    a.mechanical.outlet_traction_over_inlet_flow_pa_s_m3,
                    b.mechanical.outlet_traction_over_inlet_flow_pa_s_m3
                );
                assert_eq!(
                    a.mechanical.outlet_pressure_over_inlet_flow_pa_s_m3,
                    b.mechanical.outlet_pressure_over_inlet_flow_pa_s_m3
                );
                assert_eq!(
                    a.mechanical.outlet_over_inlet_flow,
                    b.mechanical.outlet_over_inlet_flow
                );
                assert_eq!(
                    a.mechanical.condition_estimate,
                    b.mechanical.condition_estimate
                );
                assert_eq!(
                    a.inlet_temperature_over_inlet_flow_k_s_m3,
                    b.inlet_temperature_over_inlet_flow_k_s_m3
                );
                assert_eq!(
                    a.outlet_temperature_over_inlet_flow_k_s_m3,
                    b.outlet_temperature_over_inlet_flow_k_s_m3
                );
                assert_eq!(
                    a.inlet_heat_power_over_inlet_flow_pa,
                    b.inlet_heat_power_over_inlet_flow_pa
                );
                assert_eq!(
                    a.outlet_heat_power_over_inlet_flow_pa,
                    b.outlet_heat_power_over_inlet_flow_pa
                );
            }
        }
    }
}

#[test]
fn g1_acoustic_mode_matches_the_original_coupled_nsf_dispersion() {
    let state = medium(2, 1.0);
    let gas = state.gas();
    let area = core::f64::consts::PI * RADIUS * RADIUS;
    let length = 0.04;
    for omega in [100.0, 1000.0, 10000.0] {
        let cylinder = prepare(&state, length, omega, [Isothermal, Adiabatic]);
        // This is the material's coupled NSF branch. complex_wavenumber adds
        // weak Stokes--Kirchhoff attenuation and is NOT this exact oracle.
        let k = state.causal_wavenumber(omega).unwrap();
        let traction = C64::from_re(omega * gas.density / area) / k;
        let pressure = traction + C64::new(0.0, 4.0 * gas.dynamic_viscosity / (3.0 * area)) * k;
        let theta = pressure.scale(1.0 / gas.pressure) - k.scale(1.0 / (omega * area));
        let heat = C64::new(0.0, -gas.thermal_conductivity * area * gas.temperature) * k * theta;
        let inlet = [traction, theta, C64::ONE, heat];
        let (sin, cos) = det::sin_cos(k.re * length);
        let propagation = C64::new(cos, sin).scale(det::exp(-k.im * length));
        let outlet = inlet.map(|v| v * propagation);
        let roots = cylinder.references().map(det::sqrt);
        let mut incoming = [C64::ZERO; 4];
        let mut outgoing = incoming;
        for i in 0..2 {
            incoming[i] =
                (inlet[i].scale(1.0 / roots[i]) + inlet[i + 2].scale(roots[i])).scale(0.5);
            outgoing[i] =
                (inlet[i].scale(1.0 / roots[i]) - inlet[i + 2].scale(roots[i])).scale(0.5);
            incoming[i + 2] =
                (outlet[i].scale(1.0 / roots[i]) - outlet[i + 2].scale(roots[i])).scale(0.5);
            outgoing[i + 2] =
                (outlet[i].scale(1.0 / roots[i]) + outlet[i + 2].scale(roots[i])).scale(0.5);
        }
        let scale = incoming
            .iter()
            .chain(&outgoing)
            .map(|v| v.abs())
            .fold(0.0, f64::max);
        let scattering = cylinder.scattering();
        for i in 0..4 {
            let actual = (0..4).fold(C64::ZERO, |sum, j| {
                sum + scattering[4 * i + j] * incoming[j]
            });
            near(actual, outgoing[i], scale, 2e-11);
        }
    }
}

// Primitive dimensional mass, momentum, Fourier and caloric equations in
// [traction, T'/T0, U, Q]. Reservoir temperatures are solved independently as
// theta_j=theta/(1-i*omega*tau_j); their positive heat is evaluated separately.
fn rhs(state: &RelaxationState, y: [C64; 4], area: f64, omega: f64) -> ([C64; 4], [f64; 3]) {
    let gas = state.gas();
    let eta = 4.0 * gas.dynamic_viscosity / 3.0;
    let pressure =
        (y[0] + C64::new(0.0, -omega * eta) * y[1]) / C64::new(1.0, -omega * eta / gas.pressure);
    let ux = C64::new(0.0, omega) * (pressure.scale(1.0 / gas.pressure) - y[1]);
    let theta_x = y[3].scale(-1.0 / (gas.thermal_conductivity * area * gas.temperature));
    let mut caloric_temperature = y[1].scale(state.frozen_cv_j_kg_k() + gas.specific_gas_constant);
    let mut molecular = 0.0;
    for reservoir in state.processes() {
        let internal = y[1] / C64::new(1.0, -omega * reservoir.relaxation_time_s);
        caloric_temperature = caloric_temperature + internal.scale(reservoir.heat_capacity_j_kg_k);
        molecular += reservoir.heat_capacity_j_kg_k / reservoir.relaxation_time_s
            * (y[1] - internal).norm_sq();
    }
    (
        [
            C64::new(0.0, omega * gas.density / area) * y[2],
            theta_x,
            ux.scale(area),
            C64::new(0.0, omega * area)
                * (caloric_temperature.scale(gas.density * gas.temperature) - pressure),
        ],
        [
            0.5 * area * eta * ux.norm_sq(),
            0.5 * area * gas.thermal_conductivity * gas.temperature * theta_x.norm_sq(),
            0.5 * area * gas.density * gas.temperature * molecular,
        ],
    )
}

#[test]
fn g1_complete_end_fields_close_independent_nsf_and_each_heat_source() {
    let state = medium(2, 1.0);
    let gas = state.gas();
    let area = core::f64::consts::PI * RADIUS * RADIUS;
    let length = 0.04;
    for omega in [200.0, 2000.0, 10000.0] {
        for ends in [
            [Isothermal, Isothermal],
            [Isothermal, Adiabatic],
            [Adiabatic, Adiabatic],
        ] {
            let cylinder = prepare(&state, length, omega, ends);
            let result = cylinder
                .response(
                    Some(C64::new(1.1, 0.25).scale(cylinder.references()[0])),
                    ComposeOptions::default(),
                )
                .unwrap();
            let m = result.mechanical;
            let inlet = [
                m.inlet_traction_impedance_pa_s_m3,
                result
                    .inlet_temperature_over_inlet_flow_k_s_m3
                    .scale(1.0 / gas.temperature),
                C64::ONE,
                result.inlet_heat_power_over_inlet_flow_pa,
            ];
            let outlet = [
                m.outlet_traction_over_inlet_flow_pa_s_m3,
                result
                    .outlet_temperature_over_inlet_flow_k_s_m3
                    .scale(1.0 / gas.temperature),
                m.outlet_over_inlet_flow,
                result.outlet_heat_power_over_inlet_flow_pa,
            ];
            for (i, y) in [inlet, outlet].iter().enumerate() {
                match ends[i] {
                    Isothermal => assert_eq!(y[1], C64::ZERO),
                    Adiabatic => assert_eq!(y[3], C64::ZERO),
                }
                let eta = 4.0 * gas.dynamic_viscosity / 3.0;
                let pressure = (y[0] + C64::new(0.0, -omega * eta) * y[1])
                    / C64::new(1.0, -omega * eta / gas.pressure);
                near(
                    [
                        m.inlet_pressure_over_flow_pa_s_m3,
                        m.outlet_pressure_over_inlet_flow_pa_s_m3,
                    ][i],
                    pressure,
                    pressure.abs(),
                    3e-12,
                );
            }
            let steps = 2048;
            let dx = length / steps as f64;
            let mut y = inlet;
            let mut heat = [0.0; 3];
            for step in 0..=steps {
                let (k1, source) = rhs(&state, y, area, omega);
                let weight = if step == 0 || step == steps {
                    1.0
                } else if step % 2 == 0 {
                    2.0
                } else {
                    4.0
                };
                for i in 0..3 {
                    heat[i] += weight * source[i];
                }
                if step == steps {
                    break;
                }
                let advance = |k: [C64; 4], h: f64| core::array::from_fn(|i| y[i] + k[i].scale(h));
                let k2 = rhs(&state, advance(k1, dx * 0.5), area, omega).0;
                let k3 = rhs(&state, advance(k2, dx * 0.5), area, omega).0;
                let k4 = rhs(&state, advance(k3, dx), area, omega).0;
                y = core::array::from_fn(|i| {
                    y[i] + (k1[i] + k2[i].scale(2.0) + k3[i].scale(2.0) + k4[i]).scale(dx / 6.0)
                });
            }
            for i in 0..4 {
                let scale = inlet[i].abs().max(outlet[i].abs()).max(1.0);
                near(y[i], outlet[i], scale, 2e-9);
            }
            heat.iter_mut().for_each(|v| *v *= dx / 3.0);
            assert!(
                heat.iter().all(|v| *v > 0.0),
                "inactive physical source: {heat:?}"
            );
            let flux = |z: [C64; 4]| 0.5 * (z[0] * z[2].conj() + z[1] * z[3].conj()).re;
            let loss = flux(inlet) - flux(outlet);
            let total: f64 = heat.iter().sum();
            assert!(
                (loss / total - 1.0).abs() < 2e-9,
                "work {loss}, sources {heat:?}"
            );
            assert!(
                heat[2] > 1e-4 * total,
                "reservoirs must affect physical work"
            );
        }
    }
}

#[test]
fn g1_closed_adiabatic_compression_retains_equilibrium_and_frozen_limits() {
    let omega = 1000.0;
    let length = 1e-6;
    let volume = core::f64::consts::PI * RADIUS * RADIUS * length;
    for time_scale in [1e-8, 1.0, 1e8] {
        let state = medium(2, time_scale);
        let result = prepare(&state, length, omega, [Adiabatic, Adiabatic])
            .response(None, ComposeOptions::default())
            .unwrap();
        let modulus =
            result.mechanical.inlet_pressure_over_flow_pa_s_m3 * C64::new(0.0, -omega * volume);
        let expected = state.caloric_response(omega).unwrap().bulk_modulus_pa;
        near(modulus, expected, expected.abs(), 2e-7);
        assert!(modulus.im <= 0.0);
        if time_scale != 1.0 {
            let gas = state.gas();
            let speed = if time_scale < 1.0 {
                gas.sound_speed
            } else {
                state.frozen_sound_speed_m_s()
            };
            near(
                modulus,
                C64::from_re(gas.density * speed * speed),
                expected.abs(),
                2e-7,
            );
        }
    }
}

#[test]
fn g3_molecular_full_scattering_retains_opaque_composition_and_passivity() {
    let state = medium(2, 1.0);
    for length in [0.04, 10000.0] {
        let full = prepare(&state, length, 1000.0, [Isothermal, Adiabatic]);
        let half = prepare(&state, length * 0.5, 1000.0, [Isothermal, Adiabatic]);
        let (split, _) =
            redheffer_star(&half.scattering(), &half.scattering(), 2, LIMIT, || false).unwrap();
        let s = full.scattering();
        for i in 0..16 {
            near(s[i], split[i], 1.0, 3e-12);
        }
        for i in 0..4 {
            assert!((0..4).map(|j| s[4 * j + i].norm_sq()).sum::<f64>() <= 1.0 + 3e-11);
            for j in 0..4 {
                near(s[4 * i + j], s[4 * j + i], 1.0, 3e-12);
            }
        }
        assert!(
            full.response(
                Some(C64::from_re(full.references()[0])),
                ComposeOptions::default()
            )
            .unwrap()
            .mechanical
            .inlet_traction_impedance_pa_s_m3
            .re > 0.0
        );
    }
}

#[test]
fn g4_molecular_preparation_refuses_invalid_or_cancelled_work_atomically() {
    let state = medium(2, 1.0);
    let before = state.clone();
    let ends = [Isothermal, Adiabatic];
    for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(
            BulkThermalCylinder::from_relaxation(&state, invalid, 0.04, 1000.0, ends, LIMIT)
                .is_err()
        );
        assert!(
            BulkThermalCylinder::from_relaxation(&state, RADIUS, invalid, 1000.0, ends, LIMIT)
                .is_err()
        );
        assert!(
            BulkThermalCylinder::from_relaxation(&state, RADIUS, 0.04, invalid, ends, LIMIT)
                .is_err()
        );
    }
    assert!(BulkThermalCylinder::from_relaxation(&state, RADIUS, 0.04, 1000.0, ends, 1.0).is_err());
    let mut polls = 0;
    BulkThermalCylinder::from_relaxation_controlled(
        &state,
        RADIUS,
        10000.0,
        1000.0,
        ends,
        LIMIT,
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    assert!(polls > 5);
    for stop in [1, 3, polls - 1, polls] {
        let mut count = 0;
        let result = BulkThermalCylinder::from_relaxation_controlled(
            &state,
            RADIUS,
            10000.0,
            1000.0,
            ends,
            LIMIT,
            || {
                count += 1;
                count == stop
            },
        );
        assert!(matches!(
            result,
            Err(BulkThermalError::Duct(DuctError::Cancelled { .. }))
        ));
    }
    assert_eq!(state, before);
}
