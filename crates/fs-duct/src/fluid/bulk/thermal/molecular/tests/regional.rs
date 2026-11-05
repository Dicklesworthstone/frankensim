use super::*;
use crate::fluid::bulk::thermal::{BulkThermalRegion, BulkThermalResponse};

// A second independently authored ideal gas changes density, sound speed,
// shear/conductive transport and reservoir spectrum, at the SAME equilibrium
// pressure and temperature. No fitted attenuation or average gas enters.
fn other_medium(count: usize, temperature_k: f64) -> RelaxationState {
    let gas = GasState::try_new(
        &GasSpec {
            molar_mass: 0.044,
            gamma: 1.3,
            sutherland_beta: 0.28,
            sutherland_s: 180.0,
            conductivity: ConductivityModel::Eucken,
        },
        temperature_k,
        101325.0,
    )
    .unwrap();
    state_with_gas(gas, count, 3.5, "second-authored")
}

fn region(state: &RelaxationState, length_m: f64) -> BulkThermalRelaxationRegion<'_> {
    BulkThermalRelaxationRegion {
        material: state,
        length_m,
    }
}

fn chain(
    regions: &[BulkThermalRelaxationRegion<'_>],
    omega: f64,
    ends: [ThermalEndCondition; 2],
) -> BulkThermalCylinder {
    BulkThermalCylinder::from_relaxation_regions(
        regions,
        RADIUS,
        omega,
        ends,
        ComposeOptions::default(),
    )
    .unwrap()
}

fn fields(response: BulkThermalResponse, temperature: f64) -> [[C64; 4]; 2] {
    let m = response.mechanical;
    [
        [
            m.inlet_traction_impedance_pa_s_m3,
            response
                .inlet_temperature_over_inlet_flow_k_s_m3
                .scale(1.0 / temperature),
            C64::ONE,
            response.inlet_heat_power_over_inlet_flow_pa,
        ],
        [
            m.outlet_traction_over_inlet_flow_pa_s_m3,
            response
                .outlet_temperature_over_inlet_flow_k_s_m3
                .scale(1.0 / temperature),
            m.outlet_over_inlet_flow,
            response.outlet_heat_power_over_inlet_flow_pa,
        ],
    ]
}

fn pressure(state: &RelaxationState, y: [C64; 4], omega: f64) -> C64 {
    let gas = state.gas();
    let eta = 4.0 * gas.dynamic_viscosity / 3.0;
    (y[0] + C64::new(0.0, -omega * eta) * y[1]) / C64::new(1.0, -omega * eta / gas.pressure)
}

#[test]
fn g3_regional_molecular_single_and_empty_sources_retain_original_arithmetic() {
    let options = ComposeOptions::default();
    for count in [0, 2] {
        let a = medium(count, 1.0);
        let b = other_medium(count, a.gas().temperature);
        let regions = [region(&a, 0.025), region(&b, 0.03)];
        for omega in [0.1, 2000.0, 10000.0] {
            for ends in [[Isothermal, Adiabatic], [Adiabatic, Adiabatic]] {
                let one = chain(&regions[..1], omega, ends);
                let native = prepare(&a, regions[0].length_m, omega, ends);
                assert_eq!(one.scattering(), native.scattering());
                assert_eq!(one.references(), native.references());
                assert_eq!(one.traction_scattering(), native.traction_scattering());
                if count != 0 {
                    continue;
                }
                let classical: Vec<_> = regions
                    .iter()
                    .map(|r| BulkThermalRegion {
                        material: BulkThermoviscousEquilibrium::from_material(
                            r.material.gas(),
                            0.0,
                        )
                        .unwrap(),
                        length_m: r.length_m,
                    })
                    .collect();
                let original =
                    BulkThermalCylinder::from_regions(&classical, RADIUS, omega, ends, options)
                        .unwrap();
                let molecular = chain(&regions, omega, ends);
                assert_eq!(molecular.scattering(), original.scattering());
                assert_eq!(molecular.references(), original.references());
                assert_eq!(
                    molecular.traction_scattering(),
                    original.traction_scattering()
                );
                let load = Some(C64::new(1.2, 0.3).scale(original.references()[0]));
                let x = molecular.response(load, options).unwrap();
                let y = original.response(load, options).unwrap();
                assert_eq!(
                    fields(x, a.gas().temperature),
                    fields(y, a.gas().temperature)
                );
                assert_eq!(
                    x.mechanical.inlet_pressure_over_flow_pa_s_m3,
                    y.mechanical.inlet_pressure_over_flow_pa_s_m3
                );
                assert_eq!(
                    x.mechanical.outlet_pressure_over_inlet_flow_pa_s_m3,
                    y.mechanical.outlet_pressure_over_inlet_flow_pa_s_m3
                );
                assert_eq!(
                    x.mechanical.condition_estimate,
                    y.mechanical.condition_estimate
                );
            }
        }
    }
}

#[test]
fn g1_molecular_contacts_preserve_full_nsf_fields_and_separate_source_heat() {
    let a = medium(2, 1.0);
    let b = other_medium(1, a.gas().temperature);
    let area = core::f64::consts::PI * RADIUS * RADIUS;
    for regions in [
        [region(&a, 0.025), region(&b, 0.03)],
        [region(&b, 0.03), region(&a, 0.025)],
    ] {
        for omega in [200.0, 2000.0, 10000.0] {
            for ends in [[Isothermal, Adiabatic], [Adiabatic, Adiabatic]] {
                let cylinder = chain(&regions, omega, ends);
                let result = cylinder
                    .response(
                        Some(C64::new(1.2, 0.3).scale(cylinder.references()[0])),
                        ComposeOptions::default(),
                    )
                    .unwrap();
                let [inlet, outlet] = fields(result, a.gas().temperature);
                for (state, y, actual) in [
                    (
                        regions[0].material,
                        inlet,
                        result.mechanical.inlet_pressure_over_flow_pa_s_m3,
                    ),
                    (
                        regions[1].material,
                        outlet,
                        result.mechanical.outlet_pressure_over_inlet_flow_pa_s_m3,
                    ),
                ] {
                    let expected = pressure(state, y, omega);
                    near(actual, expected, expected.abs(), 3e-12);
                }
                let mut y = inlet;
                let mut scale = inlet.map(|z| z.abs());
                let mut total_heat = 0.0;
                for (index, r) in regions.iter().enumerate() {
                    // Integrate the independent DIMENSIONAL conservation
                    // equations and each original reservoir's own slip heat.
                    // The interface carries this same y into the second gas;
                    // no local thermal termination or material averaging.
                    let steps = 2048;
                    let dx = r.length_m / steps as f64;
                    let mut heat = [0.0; 3];
                    for step in 0..=steps {
                        let (k1, source) = rhs(r.material, y, area, omega);
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
                        for i in 0..4 {
                            scale[i] = scale[i].max(y[i].abs());
                        }
                        if step == steps {
                            break;
                        }
                        let advance =
                            |k: [C64; 4], h: f64| core::array::from_fn(|i| y[i] + k[i].scale(h));
                        let k2 = rhs(r.material, advance(k1, dx * 0.5), area, omega).0;
                        let k3 = rhs(r.material, advance(k2, dx * 0.5), area, omega).0;
                        let k4 = rhs(r.material, advance(k3, dx), area, omega).0;
                        y = core::array::from_fn(|i| {
                            y[i] + (k1[i] + k2[i].scale(2.0) + k3[i].scale(2.0) + k4[i])
                                .scale(dx / 6.0)
                        });
                    }
                    heat.iter_mut().for_each(|v| *v *= dx / 3.0);
                    assert!(
                        heat.iter().all(|v| *v > 0.0),
                        "inactive region {index}: {heat:?}"
                    );
                    assert!(heat[2] > 1e-5 * heat.iter().sum::<f64>());
                    total_heat += heat.iter().sum::<f64>();
                    if index == 0 {
                        let left = pressure(regions[0].material, y, omega);
                        let right = pressure(regions[1].material, y, omega);
                        assert!((left - right).abs() > 1e-7 * left.abs());
                        assert!(y[1].abs() > 0.0 && y[3].abs() > 0.0);
                    }
                }
                for i in 0..4 {
                    near(y[i], outlet[i], scale[i].max(1.0), 3e-9);
                }
                let flux = |z: [C64; 4]| 0.5 * (z[0] * z[2].conj() + z[1] * z[3].conj()).re;
                let lost_work = flux(inlet) - flux(outlet);
                assert!(
                    (lost_work / total_heat - 1.0).abs() < 3e-9,
                    "boundary work {lost_work}, regional source heat {total_heat}"
                );
            }
        }
    }
}

#[test]
fn g3_molecular_regions_keep_reciprocal_opaque_composition_and_split_identity() {
    let a = medium(2, 1.0);
    let b = other_medium(1, a.gas().temperature);
    let omega = 2000.0;
    let ends = [Adiabatic, Adiabatic];
    for length in [0.04, 10000.0] {
        let original = prepare(&a, length, omega, ends);
        let split = chain(&[region(&a, length / 4.0); 4], omega, ends);
        for (x, y) in original.scattering().iter().zip(split.scattering()) {
            near(*x, y, 1.0, 3e-11);
        }
        let half_regions = [region(&a, length * 0.25), region(&b, length * 0.25)];
        let half = chain(&half_regions, omega, [Isothermal, Adiabatic]);
        let full = chain(
            &[
                half_regions[0],
                half_regions[1],
                half_regions[0],
                half_regions[1],
            ],
            omega,
            ends,
        );
        let (joined, _) =
            redheffer_star(&half.scattering(), &half.scattering(), 2, LIMIT, || false).unwrap();
        let s = full.scattering();
        for i in 0..4 {
            for j in 0..4 {
                near(s[4 * i + j], joined[4 * i + j], 1.0, 3e-11);
                near(s[4 * i + j], s[4 * j + i], 1.0, 3e-11);
            }
        }
        for drive in 0..16 {
            let input: [C64; 4] = core::array::from_fn(|i| {
                C64::new(
                    (i + 1) as f64,
                    if drive & (1 << i) == 0 { 0.7 } else { -0.7 },
                )
            });
            let incoming: f64 = input.iter().map(|z| z.norm_sq()).sum();
            let outgoing: f64 = (0..4)
                .map(|i| {
                    (0..4)
                        .fold(C64::ZERO, |v, j| v + s[4 * i + j] * input[j])
                        .norm_sq()
                })
                .sum();
            assert!(outgoing <= incoming * (1.0 + 5e-11));
        }
    }
}

#[test]
fn g4_regional_molecular_admission_and_late_cancellation_preserve_sources() {
    let a = medium(2, 1.0);
    let b = other_medium(1, a.gas().temperature);
    let before = [a.clone(), b.clone()];
    let regions = [region(&a, 0.025), region(&b, 0.03)];
    let ends = [Isothermal, Adiabatic];
    let options = ComposeOptions {
        maximum_blocks: 2,
        ..ComposeOptions::default()
    };
    let build = |r: &[BulkThermalRelaxationRegion<'_>], o| {
        BulkThermalCylinder::from_relaxation_regions(r, RADIUS, 2000.0, ends, o)
    };
    assert!(build(&regions, options).is_ok());
    assert!(build(&[], options).is_err());
    assert!(
        build(
            &regions,
            ComposeOptions {
                maximum_blocks: 1,
                ..options
            }
        )
        .is_err()
    );
    assert!(
        build(
            &regions,
            ComposeOptions {
                maximum_condition: 1.0,
                ..options
            }
        )
        .is_err()
    );
    for length in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(build(&[regions[0], region(&b, length)], options).is_err());
    }
    let warm = other_medium(1, 305.0);
    assert!(build(&[regions[0], region(&warm, 0.03)], options).is_err());
    let mut polls = 0;
    BulkThermalCylinder::from_relaxation_regions_controlled(
        &regions,
        RADIUS,
        2000.0,
        ends,
        options,
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    assert!(polls > 5);
    for stop in [1, polls / 2, polls - 1, polls] {
        let mut count = 0;
        let result = BulkThermalCylinder::from_relaxation_regions_controlled(
            &regions,
            RADIUS,
            2000.0,
            ends,
            options,
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
    assert_eq!([a, b], before);
}
