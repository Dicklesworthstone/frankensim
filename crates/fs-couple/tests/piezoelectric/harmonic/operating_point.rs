//! G1/G3/G4: DC thermal bias and AC increments retain distinct physical roles.
//! Uniform axial fields give exact capacitor, strain-envelope and acoustic
//! power oracles without inferring a measured ceramic or a continuum accuracy.
#[path = "operating_point/relaxation.rs"]
mod relaxation;

use super::*;
use fs_couple::piezoelectric::spatial::harmonic::acoustic::{
    HarmonicAcousticPort, PlaneWavePortControls,
};
use fs_couple::piezoelectric::spatial::transient::load::reactive::ReactiveElectrodeLoad;
use fs_material::fluid_acoustics::FluidAcousticState;
use std::collections::BTreeSet;

const COUPLING: f64 = 12.0;
const ALPHA: f64 = 2e-6;
const PYRO: f64 = -3e-5;
const REFERENCE: f64 = 300.0;
const HOT: f64 = 310.0;
const AREA: f64 = LENGTH * LENGTH;
const VOLUME: f64 = AREA * LENGTH;

fn thermal_model<'a>(cube: &'a Cube, temperature: f64, cx: &Cx<'_>) -> Model<'a> {
    let (_, law) = crate::thermal_bias::source();
    Model::assemble_thermal_materials(
        &cube.nodes,
        &cube.tets,
        vec![law.at_temperature(temperature).unwrap()],
        FRAME,
        budget(),
        cx,
    )
    .unwrap()
}

fn fixed(cube: &Cube, strain: f64) -> Vec<(usize, f64)> {
    cube.nodes
        .iter()
        .enumerate()
        .flat_map(|(node, position)| {
            let u = [0.0, 0.0, strain * position[2]];
            (0..3).map(move |axis| (3 * node + axis, u[axis]))
        })
        .collect()
}

fn oscillating(cube: &Cube, strain: C64) -> Vec<(usize, C64)> {
    cube.nodes
        .iter()
        .enumerate()
        .flat_map(|(node, position)| {
            let u = [ZERO, ZERO, scale(strain, position[2])];
            (0..3).map(move |axis| (3 * node + axis, u[axis]))
        })
        .collect()
}

fn top_faces(cube: &Cube) -> Vec<[usize; 3]> {
    let mut faces = BTreeSet::new();
    for tet in &cube.tets {
        for omit in 0..4 {
            let face = (0..4)
                .filter(|&i| i != omit)
                .map(|i| tet[i])
                .collect::<Vec<_>>();
            if face.iter().all(|node| cube.top.contains(node)) {
                let mut face = [face[0], face[1], face[2]];
                face.sort_unstable();
                faces.insert(face);
            }
        }
    }
    faces.into_iter().collect()
}

#[test]
fn g1_g3_hot_open_and_short_operating_points_share_the_unbiased_ac_tangent() {
    with_cx(|cx, _| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let zero = vec![[0.0; 3]; n];
        let charge = vec![0.0; n];
        let ac_zero = vec![[ZERO; 3]; n];
        let ac_charge = vec![ZERO; n];
        let dc_supports = fixed(&cube, 0.0);
        let ac_supports = oscillating(&cube, ZERO);
        let omega = 2e6;
        let capacitance = EPSILON[2] * LENGTH;
        let voltage = C64::new(3.0, -4.0);
        let current = mul(C64::new(0.0, -omega * capacitance), voltage);
        let ac_energy = 0.25 * capacitance * voltage.norm_sq();
        for temperature in [REFERENCE, HOT] {
            let model = thermal_model(&cube, temperature, cx);
            let tangents = [model.materials()[0].clone()];
            let tangent =
                Model::assemble(&cube.nodes, &cube.tets, &tangents, FRAME, budget(), cx).unwrap();
            let expansion = ALPHA * (temperature - REFERENCE);
            let polarization = PYRO * (temperature - REFERENCE);
            for open in [false, true] {
                let dc = if open {
                    Condition::Charge(0.0)
                } else {
                    Condition::Voltage(0.0)
                };
                let ac = if open {
                    HCondition::Current(current)
                } else {
                    HCondition::Voltage(voltage)
                };
                let dc_field = if open {
                    -polarization / EPSILON[2]
                } else {
                    0.0
                };
                let dc_energy = 0.5
                    * VOLUME
                    * (STIFFNESS[2] * expansion * expansion + EPSILON[2] * dc_field * dc_field);
                let point = model
                    .harmonic_operating_point(
                        &dc_supports,
                        &electrodes(&cube, 0.0, dc),
                        &zero,
                        &charge,
                        controls(),
                        cx,
                    )
                    .unwrap();
                let response = point
                    .solve_harmonic(
                        omega,
                        &[],
                        &ac_supports,
                        &helectrodes(&cube, ZERO, ac),
                        &[],
                        &ac_zero,
                        &ac_charge,
                        hc(),
                        cx,
                    )
                    .unwrap();
                let reference = tangent
                    .solve_harmonic(
                        omega,
                        &ac_supports,
                        &helectrodes(&cube, ZERO, ac),
                        &[],
                        &ac_zero,
                        &ac_charge,
                        hc(),
                        cx,
                    )
                    .unwrap();
                near(point.state().stored_energy_j, dc_energy, 2e-17, 3e-8);
                near(
                    point.state().electrodes[1].voltage_v,
                    -LENGTH * dc_field,
                    2e-8,
                    3e-8,
                );
                let harmonic = &response.harmonic;
                cnear(harmonic.electrodes[1].voltage_v, voltage, 2e-8, 3e-8);
                cnear(
                    harmonic.electrodes[1].charge_c,
                    scale(voltage, capacitance),
                    2e-18,
                    3e-8,
                );
                cnear(
                    harmonic.electrodes[1].source_current_a,
                    current,
                    2e-12,
                    3e-8,
                );
                for (cell, unbiased) in harmonic.cells.iter().zip(&reference.cells) {
                    cnear(cell.strain_mandel[2], ZERO, 1e-14, 0.0);
                    cnear(
                        cell.electric_field_v_m[2],
                        scale(voltage, -1.0 / LENGTH),
                        2e-5,
                        3e-8,
                    );
                    cnear(
                        cell.displacement_c_m2[2],
                        scale(voltage, -EPSILON[2] / LENGTH),
                        2e-12,
                        3e-8,
                    );
                    cnear(
                        cell.stress_mandel_pa[2],
                        scale(voltage, COUPLING / LENGTH),
                        2e-3,
                        3e-8,
                    );
                    cnear(
                        cell.stress_mandel_pa[2],
                        unbiased.stress_mandel_pa[2],
                        2e-3,
                        3e-8,
                    );
                    cnear(
                        cell.displacement_c_m2[2],
                        unbiased.displacement_c_m2[2],
                        2e-12,
                        3e-8,
                    );
                }
                near(harmonic.average_stored_energy_j, ac_energy, 2e-18, 3e-8);
                near(
                    harmonic.average_stored_energy_j,
                    reference.average_stored_energy_j,
                    2e-18,
                    3e-8,
                );
                near(
                    response.average_total_stored_energy_j,
                    dc_energy + ac_energy,
                    2e-17,
                    3e-8,
                );
                near(response.maximum_total_strain_norm, 0.0, 1e-14, 0.0);
                near(
                    response.maximum_elastic_strain_norm,
                    expansion.abs(),
                    1e-12,
                    3e-8,
                );
                near(
                    response.maximum_electric_field_v_m,
                    dc_field.abs() + voltage.norm_sq().sqrt() / LENGTH,
                    2e-4,
                    3e-8,
                );
                near(harmonic.average_electrical_power_w, 0.0, 1e-12, 0.0);
                near(harmonic.power_balance_residual_w, 0.0, 1e-12, 0.0);
                assert_eq!(harmonic.average_material_heat_w, 0.0);
                assert_eq!(response.dc_load_heat_w, 0.0);
            }
        }
    });
}

#[test]
fn g1_g4_full_cycle_limits_include_dc_total_elastic_and_electric_bias_between_quadratures() {
    with_cx(|cx, gate| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let model = thermal_model(&cube, HOT, cx);
        let zero = vec![[0.0; 3]; n];
        let charge = vec![0.0; n];
        let ac_zero = vec![[ZERO; 3]; n];
        let ac_charge = vec![ZERO; n];
        let strain = C64::new(3e-4, 3e-4);
        let amplitude = strain.norm_sq().sqrt();
        let expansion = ALPHA * (HOT - REFERENCE);
        let ac_supports = oscillating(&cube, strain);
        let ac_electrodes = helectrodes(&cube, ZERO, HCondition::Voltage(ZERO));
        let mut generous = hc();
        generous.field.maximum_strain_norm = 2e-3;
        for dc_strain in [5.9e-4_f64, -5.7e-4] {
            // All four quadratures pass 1e-3. The true extrema at pi/4 and
            // 5*pi/4 violate TOTAL strain in one case and ELASTIC strain in
            // the other; checking only the AC phasor cannot catch either.
            assert!(dc_strain.abs() + strain.re.abs() < 1e-3);
            assert!((dc_strain - expansion).abs() + strain.re.abs() < 1e-3);
            let point = model
                .harmonic_operating_point(
                    &fixed(&cube, dc_strain),
                    &electrodes(&cube, 0.0, Condition::Voltage(0.0)),
                    &zero,
                    &charge,
                    controls(),
                    cx,
                )
                .unwrap();
            let accepted = point
                .solve_harmonic(
                    1e5,
                    &[],
                    &ac_supports,
                    &ac_electrodes,
                    &[],
                    &ac_zero,
                    &ac_charge,
                    generous,
                    cx,
                )
                .unwrap();
            near(
                accepted.maximum_total_strain_norm,
                dc_strain.abs() + amplitude,
                2e-12,
                3e-8,
            );
            near(
                accepted.maximum_elastic_strain_norm,
                (dc_strain - expansion).abs() + amplitude,
                2e-12,
                3e-8,
            );
            assert!(
                point
                    .solve_harmonic(
                        1e5,
                        &[],
                        &ac_supports,
                        &ac_electrodes,
                        &[],
                        &ac_zero,
                        &ac_charge,
                        hc(),
                        cx
                    )
                    .is_err()
            );
            assert!(
                point
                    .solve_harmonic(
                        1e5,
                        &[],
                        &ac_supports[1..],
                        &ac_electrodes,
                        &[],
                        &ac_zero,
                        &ac_charge,
                        generous,
                        cx
                    )
                    .is_err()
            );
            near(
                point.state().cells[0].strain_mandel[2],
                dc_strain,
                1e-12,
                3e-8,
            );
        }

        let dc_field = 1.18e6;
        let ac_field = C64::new(6e5, 6e5);
        let ac_supports = oscillating(&cube, ZERO);
        let ac_electrodes = helectrodes(&cube, ZERO, HCondition::Voltage(scale(ac_field, -LENGTH)));
        let point = model
            .harmonic_operating_point(
                &fixed(&cube, 0.0),
                &electrodes(&cube, 0.0, Condition::Voltage(-LENGTH * dc_field)),
                &zero,
                &charge,
                controls(),
                cx,
            )
            .unwrap();
        generous.field.maximum_electric_field_v_m = 3e6;
        let accepted = point
            .solve_harmonic(
                1e5,
                &[],
                &ac_supports,
                &ac_electrodes,
                &[],
                &ac_zero,
                &ac_charge,
                generous,
                cx,
            )
            .unwrap();
        assert!(dc_field + ac_field.re.abs() < hc().field.maximum_electric_field_v_m);
        near(
            accepted.maximum_electric_field_v_m,
            dc_field + ac_field.norm_sq().sqrt(),
            2e-3,
            3e-8,
        );
        assert!(
            point
                .solve_harmonic(
                    1e5,
                    &[],
                    &ac_supports,
                    &ac_electrodes,
                    &[],
                    &ac_zero,
                    &ac_charge,
                    hc(),
                    cx
                )
                .is_err()
        );
        let before = point.state().stored_energy_j;
        gate.request();
        assert!(matches!(
            point.solve_harmonic(
                1e5,
                &[],
                &ac_supports,
                &ac_electrodes,
                &[],
                &ac_zero,
                &ac_charge,
                generous,
                cx
            ),
            Err(TetPiezoelectricError::Cancelled)
        ));
        assert_eq!(point.state().stored_energy_j, before);
    });
}

#[test]
fn g1_acoustic_support_power_contains_only_motion_about_the_hot_dc_state() {
    with_cx(|cx, _| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let model = thermal_model(&cube, HOT, cx);
        let zero = vec![[0.0; 3]; n];
        let charge = vec![0.0; n];
        let point = model
            .harmonic_operating_point(
                &fixed(&cube, 0.0),
                &electrodes(&cube, 0.0, Condition::Charge(0.0)),
                &zero,
                &charge,
                controls(),
                cx,
            )
            .unwrap();
        let trace = top_faces(&cube);
        let port = model
            .plane_wave_acoustic_port(
                &trace,
                AREA,
                PlaneWavePortControls {
                    maximum_surface_dofs: 64,
                    coplanarity_tolerance: 1e-10,
                    maximum_wavenumber_rad_m: core::f64::consts::PI / LENGTH,
                    maximum_relative_amplitude: 1e-2,
                },
                cx,
            )
            .unwrap();
        let (density, sound_speed, omega) = (1000.0, 1480.0, 1e5);
        let medium = FluidAcousticState::try_new(density, sound_speed).unwrap();
        let strain = C64::new(1e-5, -2e-5);
        let velocity = mul(C64::new(0.0, -omega * LENGTH), strain);
        let flow = scale(velocity, AREA);
        let pressure = scale(velocity, density * sound_speed);
        let power = 0.5 * density * sound_speed * AREA * velocity.norm_sq();
        let energy =
            0.25 * VOLUME * (STIFFNESS[2] + COUPLING * COUPLING / EPSILON[2]) * strain.norm_sq();
        let response = point
            .solve_harmonic(
                omega,
                &[HarmonicAcousticPort {
                    port: &port,
                    medium: &medium,
                    incident_pressure_pa: ZERO,
                }],
                &oscillating(&cube, strain),
                &helectrodes(&cube, ZERO, HCondition::Current(ZERO)),
                &[],
                &vec![[ZERO; 3]; n],
                &vec![ZERO; n],
                hc(),
                cx,
            )
            .unwrap();
        let harmonic = &response.harmonic;
        cnear(
            harmonic.electrodes[1].voltage_v,
            scale(strain, COUPLING * LENGTH / EPSILON[2]),
            2e-7,
            3e-8,
        );
        cnear(harmonic.electrodes[1].charge_c, ZERO, 2e-18, 0.0);
        let acoustic = &harmonic.acoustic_ports[0];
        cnear(acoustic.volume_velocity_m3_s, flow, 1e-18, 3e-8);
        cnear(acoustic.outgoing_pressure_pa, pressure, 1e-6, 3e-8);
        near(acoustic.average_outgoing_power_w, power, 2e-13, 3e-8);
        near(harmonic.average_acoustic_power_w, power, 2e-13, 3e-8);
        near(harmonic.average_mechanical_power_w, power, 2e-13, 3e-8);
        near(harmonic.average_electrical_power_w, 0.0, 2e-13, 0.0);
        near(harmonic.average_stored_energy_j, energy, 2e-17, 3e-8);
        near(
            harmonic.average_kinetic_energy_j,
            7500.0 * VOLUME * LENGTH * LENGTH * omega * omega * strain.norm_sq() / 12.0,
            2e-18,
            3e-8,
        );
        near(
            response.average_total_stored_energy_j,
            point.state().stored_energy_j + energy,
            2e-17,
            3e-8,
        );
        near(harmonic.power_balance_residual_w, 0.0, 2e-13, 0.0);
        assert_eq!(acoustic.average_incident_power_w, 0.0);
        assert_eq!(harmonic.average_material_heat_w, 0.0);
        assert_eq!(response.average_total_load_heat_w, 0.0);
    });
}

#[test]
fn g1_dc_shunt_supplies_storage_and_heat_are_separate_from_ac_power() {
    use crate::transient_load::reactive::{capacitor, inductor};

    with_cx(|cx, _| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let model = thermal_model(&cube, HOT, cx);
        let (gauge, dc_voltage, resistance, winding_resistance, inductance) =
            (1000.0, 4.0, 1e5, 2e5, 2e-3);
        let capacitance = 2.0 * EPSILON[2] * LENGTH;
        let point = model
            .harmonic_operating_point(
                &fixed(&cube, 0.0),
                &electrodes(&cube, gauge, Condition::Voltage(gauge + dc_voltage)),
                &vec![[0.0; 3]; n],
                &vec![0.0; n],
                controls(),
                cx,
            )
            .unwrap();
        let loads = [
            HLoad::Resistive(crate::transient_load::resistor(resistance, [1, 0])),
            HLoad::Reactive(ReactiveElectrodeLoad::SeriesRl {
                electrodes: [1, 0],
                inductor: inductor(inductance),
                conductor: Some(
                    crate::transient_load::resistor(winding_resistance, [1, 0]).conductor,
                ),
            }),
            HLoad::Reactive(ReactiveElectrodeLoad::Capacitor {
                electrodes: [1, 0],
                capacitor: capacitor(capacitance),
            }),
        ];
        let omega = 2e6;
        let voltage = C64::new(3.0, -4.0);
        let response = point
            .solve_harmonic(
                omega,
                &[],
                &oscillating(&cube, ZERO),
                &helectrodes(&cube, ZERO, HCondition::Voltage(voltage)),
                &loads,
                &vec![[ZERO; 3]; n],
                &vec![ZERO; n],
                hc(),
                cx,
            )
            .unwrap();
        let winding_current = dc_voltage / winding_resistance;
        let total_current = dc_voltage / resistance + winding_current;
        let dc_heat = dc_voltage * total_current;
        let capacitor_energy = 0.5 * capacitance * dc_voltage * dc_voltage;
        let inductor_energy = 0.5 * inductance * winding_current * winding_current;
        assert_eq!(response.dc_loads.len(), 3);
        for branch in &response.dc_loads {
            near(branch.voltage_v, dc_voltage, 1e-12, 0.0);
        }
        near(
            response.dc_loads[0].current_a,
            dc_voltage / resistance,
            1e-14,
            3e-8,
        );
        near(
            response.dc_loads[0].joule_power_w,
            dc_voltage * dc_voltage / resistance,
            1e-14,
            3e-8,
        );
        near(response.dc_loads[1].current_a, winding_current, 1e-14, 3e-8);
        near(
            response.dc_loads[1].stored_energy_j,
            inductor_energy,
            1e-20,
            3e-8,
        );
        near(
            response.dc_loads[2].charge_c,
            capacitance * dc_voltage,
            1e-20,
            3e-8,
        );
        near(
            response.dc_loads[2].stored_energy_j,
            capacitor_energy,
            1e-20,
            3e-8,
        );
        assert_eq!(response.dc_loads[2].current_a, 0.0);
        assert_eq!(response.dc_loads[2].joule_power_w, 0.0);
        near(
            response.dc_electrode_source_current_a[0],
            -total_current,
            1e-14,
            3e-8,
        );
        near(
            response.dc_electrode_source_current_a[1],
            total_current,
            1e-14,
            3e-8,
        );
        near(response.dc_load_heat_w, dc_heat, 1e-14, 3e-8);
        near(response.dc_electrical_supply_w, dc_heat, 1e-14, 3e-8);
        near(
            response.dc_external_stored_energy_j,
            capacitor_energy + inductor_energy,
            1e-20,
            3e-8,
        );
        let ac_winding = div(voltage, C64::new(winding_resistance, -omega * inductance));
        let ac_heat =
            0.5 * (voltage.norm_sq() / resistance + winding_resistance * ac_winding.norm_sq());
        near(response.harmonic.average_load_heat_w, ac_heat, 1e-14, 3e-8);
        near(
            response.harmonic.average_electrical_power_w,
            ac_heat,
            1e-14,
            3e-8,
        );
        near(
            response.average_total_load_heat_w,
            dc_heat + ac_heat,
            1e-14,
            3e-8,
        );
        near(
            response.average_total_external_stored_energy_j,
            capacitor_energy
                + inductor_energy
                + 0.25 * (capacitance * voltage.norm_sq() + inductance * ac_winding.norm_sq()),
            1e-20,
            3e-8,
        );
        let lossless = [HLoad::Reactive(ReactiveElectrodeLoad::SeriesRl {
            electrodes: [1, 0],
            inductor: inductor(inductance),
            conductor: None,
        })];
        // A nonzero DC terminal voltage cannot coexist with a steady ideal
        // inductor. Its AC admittance alone must not admit that impossible bias.
        assert!(
            point
                .solve_harmonic_with_inductor_bias(
                    omega,
                    &[],
                    &[(0, 0.02)],
                    &oscillating(&cube, ZERO),
                    &helectrodes(&cube, ZERO, HCondition::Voltage(voltage)),
                    &lossless,
                    &vec![[ZERO; 3]; n],
                    &vec![ZERO; n],
                    hc(),
                    cx
                )
                .is_err()
        );
        // At zero DC voltage, a lossless winding can instead retain an
        // explicitly supplied persistent current. Zero voltage alone must not
        // silently choose zero magnetic storage for this independent state.
        let point = model
            .harmonic_operating_point(
                &fixed(&cube, 0.0),
                &electrodes(&cube, gauge, Condition::Voltage(gauge)),
                &vec![[0.0; 3]; n],
                &vec![0.0; n],
                controls(),
                cx,
            )
            .unwrap();
        assert!(
            point
                .solve_harmonic(
                    omega,
                    &[],
                    &oscillating(&cube, ZERO),
                    &helectrodes(&cube, ZERO, HCondition::Voltage(voltage)),
                    &lossless,
                    &vec![[ZERO; 3]; n],
                    &vec![ZERO; n],
                    hc(),
                    cx
                )
                .is_err()
        );
        let persistent = 0.02;
        let response = point
            .solve_harmonic_with_inductor_bias(
                omega,
                &[],
                &[(0, persistent)],
                &oscillating(&cube, ZERO),
                &helectrodes(&cube, ZERO, HCondition::Voltage(voltage)),
                &lossless,
                &vec![[ZERO; 3]; n],
                &vec![ZERO; n],
                hc(),
                cx,
            )
            .unwrap();
        let magnetic = 0.5 * inductance * persistent * persistent;
        let ac_current = div(voltage, C64::new(0.0, -omega * inductance));
        near(response.dc_loads[0].current_a, persistent, 1e-14, 3e-8);
        near(response.dc_loads[0].stored_energy_j, magnetic, 1e-20, 3e-8);
        near(
            response.dc_electrode_source_current_a[0],
            -persistent,
            1e-14,
            3e-8,
        );
        near(
            response.dc_electrode_source_current_a[1],
            persistent,
            1e-14,
            3e-8,
        );
        near(
            response.average_total_external_stored_energy_j,
            magnetic + 0.25 * inductance * ac_current.norm_sq(),
            1e-20,
            3e-8,
        );
        assert_eq!(response.dc_electrical_supply_w, 0.0);
        assert_eq!(response.average_total_load_heat_w, 0.0);
    });
}
