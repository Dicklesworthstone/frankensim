//! G1/G3/G4: reciprocal waveguide loading, sensing and mechanical work.
//! Synthetic source values specify implementation references, not measured devices.
use super::*;
use fs_couple::piezoelectric::spatial::harmonic::acoustic::{
    HarmonicAcousticPort, PlaneWavePortControls,
};
use fs_material::fluid_acoustics::FluidAcousticState;
use std::collections::BTreeSet;

fn port_controls() -> PlaneWavePortControls {
    PlaneWavePortControls {
        maximum_surface_dofs: 64,
        coplanarity_tolerance: 1e-10,
        // First nonuniform mode of the authored square rigid guide.
        maximum_wavenumber_rad_m: core::f64::consts::PI / LENGTH,
        maximum_relative_amplitude: 1e-2,
    }
}
fn face_trace(cube: &Cube, nodes: &[usize]) -> Vec<[usize; 3]> {
    let mut faces = BTreeSet::new();
    for tet in &cube.tets {
        for omit in 0..4 {
            let face: Vec<_> = (0..4).filter(|&i| i != omit).map(|i| tet[i]).collect();
            if face.iter().all(|i| nodes.contains(i)) {
                let mut f = [face[0], face[1], face[2]];
                f.sort_unstable();
                faces.insert(f);
            }
        }
    }
    faces.into_iter().collect()
}
// Integrate the selected top vertex's P1 shape directly from each triangle.
fn top_weight(cube: &Cube, faces: &[[usize; 3]], node: usize) -> f64 {
    faces
        .iter()
        .filter(|f| f.contains(&node))
        .map(|f| {
            let [a, b, c] = f.map(|i| cube.nodes[i]);
            ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])).abs() / 6.0
        })
        .sum()
}
fn diagonal(matrix: &fs_sparse::Csr, dof: usize) -> f64 {
    let (columns, values) = matrix.row(dof);
    values[columns.iter().position(|&i| i == dof).unwrap()]
}
fn electrode_coupling(model: &Model<'_>, cube: &Cube, dof: usize) -> f64 {
    let (columns, values) = model.coupling().row(dof);
    columns
        .iter()
        .zip(values)
        .filter(|(i, _)| cube.top.contains(i))
        .map(|(_, v)| v)
        .sum()
}

#[test]
fn g1_voltage_actuator_loads_the_actual_fluid_and_supplies_outgoing_wave_power() {
    with_cx(|cx, _| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let materials = [material([0.0, 0.0, 12.0, 0.0, 0.0, 0.0])];
        let model =
            Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let faces = face_trace(&cube, &cube.top);
        let area = LENGTH * LENGTH;
        let port = model
            .plane_wave_acoustic_port(&faces, area, port_controls(), cx)
            .unwrap();
        let node = *cube.top.last().unwrap();
        let dof = 3 * node + 2;
        let b = top_weight(&cube, &faces, node);
        near(port.moving_area_m2(), area, 1e-20, 1e-14);
        near(
            port.weights_m2().iter().find(|(i, _)| *i == dof).unwrap().1,
            b,
            1e-20,
            1e-14,
        );
        let supports: Vec<_> = (0..3 * n)
            .filter(|&i| i != dof)
            .map(|i| (i, ZERO))
            .collect();
        let k = diagonal(&model.solid().stiffness, dof);
        let m = diagonal(&model.solid().mass, dof);
        let g = electrode_coupling(&model, &cube, dof);
        let voltage = C64::new(-2.0, 1.0);
        for medium in [
            FluidAcousticState::try_new(1.2, 343.0).unwrap(),
            FluidAcousticState::try_new(1000.0, 1480.0).unwrap(),
        ] {
            let omega = 0.5 * medium.sound_speed_m_s() * port_controls().maximum_wavenumber_rad_m;
            let z = medium.characteristic_impedance_pa_s_m() / area;
            let dynamic = C64::new(k - omega * omega * m, -omega * z * b * b);
            let motion = div(scale(voltage, -g), dynamic);
            let flow = scale(mul(C64::new(0.0, -omega), motion), b);
            let pressure = scale(flow, z);
            let expected_power = 0.5 * z * flow.norm_sq();
            let response = model
                .solve_harmonic_with_acoustic_ports(
                    omega,
                    &[HarmonicAcousticPort {
                        port: &port,
                        medium: &medium,
                        incident_pressure_pa: ZERO,
                    }],
                    None,
                    &supports,
                    &helectrodes(&cube, ZERO, HCondition::Voltage(voltage)),
                    &[],
                    &vec![[ZERO; 3]; n],
                    &vec![ZERO; n],
                    hc(),
                    cx,
                )
                .unwrap();
            cnear(response.displacement_m[node][2], motion, 1e-18, 5e-8);
            let acoustic = &response.acoustic_ports[0];
            assert_eq!(acoustic.medium, medium);
            cnear(acoustic.volume_velocity_m3_s, flow, 1e-20, 5e-8);
            cnear(acoustic.pressure_pa, pressure, 1e-10, 5e-8);
            cnear(acoustic.outgoing_pressure_pa, pressure, 1e-10, 5e-8);
            assert_eq!(acoustic.average_incident_power_w, 0.0);
            near(
                acoustic.average_outgoing_power_w,
                expected_power,
                1e-20,
                2e-7,
            );
            near(
                response.average_acoustic_power_w,
                expected_power,
                1e-20,
                2e-7,
            );
            near(
                response.average_electrical_power_w,
                expected_power,
                2e-18,
                2e-7,
            );
            near(response.average_mechanical_power_w, 0.0, 1e-20, 0.0);
            near(response.average_load_heat_w, 0.0, 1e-20, 0.0);
            near(response.average_material_heat_w, 0.0, 1e-20, 0.0);
            assert!(expected_power > 0.0);
        }
    });
}

#[test]
fn g1_incident_sound_drives_sensor_voltage_and_actual_shunt_heat_with_backreaction() {
    with_cx(|cx, _| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let materials = [material([0.0, 0.0, 12.0, 0.0, 0.0, 0.0])];
        let model =
            Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let faces = face_trace(&cube, &cube.top);
        let area = LENGTH * LENGTH;
        let port = model
            .plane_wave_acoustic_port(&faces, area, port_controls(), cx)
            .unwrap();
        let medium = FluidAcousticState::try_new(1000.0, 1480.0).unwrap();
        let node = *cube.top.last().unwrap();
        let dof = 3 * node + 2;
        let b = top_weight(&cube, &faces, node);
        let fixed: Vec<_> = (0..3 * n)
            .filter(|&i| i != dof)
            .map(|i| (i, ZERO))
            .collect();
        let k = diagonal(&model.solid().stiffness, dof);
        let m = diagonal(&model.solid().mass, dof);
        let g = electrode_coupling(&model, &cube, dof);
        let omega = 0.4 * medium.sound_speed_m_s() * port_controls().maximum_wavenumber_rad_m;
        let capacitance = EPSILON[2] * LENGTH;
        let resistance = 1.0 / (omega * capacitance);
        let z = medium.characteristic_impedance_pa_s_m() / area;
        let incident = C64::new(3000.0, -1000.0);
        let voltage_per_motion = div(
            C64::new(0.0, -omega * g),
            C64::new(1.0 / resistance, -omega * capacitance),
        );
        let dynamic = add(
            C64::new(k - omega * omega * m, -omega * z * b * b),
            scale(voltage_per_motion, g),
        );
        let motion = div(scale(incident, -2.0 * b), dynamic);
        let voltage = mul(voltage_per_motion, motion);
        let flow = scale(mul(C64::new(0.0, -omega), motion), b);
        let outgoing = add(incident, scale(flow, z));
        let heat = 0.5 * voltage.norm_sq() / resistance;
        let loads = [HLoad::Resistive(super::super::transient_load::resistor(
            resistance,
            [1, 0],
        ))];
        let response = model
            .solve_harmonic_with_acoustic_ports(
                omega,
                &[HarmonicAcousticPort {
                    port: &port,
                    medium: &medium,
                    incident_pressure_pa: incident,
                }],
                None,
                &fixed,
                &helectrodes(&cube, ZERO, HCondition::Current(ZERO)),
                &loads,
                &vec![[ZERO; 3]; n],
                &vec![ZERO; n],
                hc(),
                cx,
            )
            .unwrap();
        cnear(response.displacement_m[node][2], motion, 1e-18, 5e-8);
        cnear(response.electrodes[1].voltage_v, voltage, 1e-10, 5e-8);
        cnear(
            response.acoustic_ports[0].outgoing_pressure_pa,
            outgoing,
            1e-7,
            5e-8,
        );
        near(response.average_load_heat_w, heat, 1e-20, 2e-7);
        near(response.average_acoustic_power_w, -heat, 2e-17, 2e-7);
        near(response.average_electrical_power_w, 0.0, 1e-18, 0.0);
        near(response.average_mechanical_power_w, 0.0, 1e-18, 0.0);
        assert!(
            response.acoustic_ports[0].average_outgoing_power_w
                < response.acoustic_ports[0].average_incident_power_w
        );
        assert!(heat > 0.0);
    });
}

#[test]
fn g3_rotated_pressure_trace_keeps_flow_reaction_work_and_rigid_reflection() {
    with_cx(|cx, _| {
        let original = cube(1);
        let faces = face_trace(&original, &original.top);
        let medium = FluidAcousticState::try_new(1000.0, 1480.0).unwrap();
        let omega = 1e5;
        let displacement = C64::new(1e-10, 3e-11);
        let area = LENGTH * LENGTH;
        let expected_flow = scale(mul(C64::new(0.0, -omega), displacement), area);
        let z = medium.characteristic_impedance_pa_s_m() / area;
        let expected_power = 0.5 * z * expected_flow.norm_sq();
        for angle in [0.0, 0.61] {
            let q = rotation([1.0 / 3.0_f64.sqrt(); 3], angle);
            let frame = ContentHash([101; 32]);
            let nodes: Vec<_> = original.nodes.iter().map(|p| apply(&q, p)).collect();
            let materials = [material([0.0, 0.0, 12.0, 0.0, 0.0, 0.0])
                .rotated(frame, q)
                .unwrap()];
            let model =
                Model::assemble(&nodes, &original.tets, &materials, frame, budget(), cx).unwrap();
            let port = model
                .plane_wave_acoustic_port(&faces, area, port_controls(), cx)
                .unwrap();
            let mut reordered = faces.clone();
            reordered.reverse();
            for f in &mut reordered {
                f.swap(0, 2);
            }
            let reordered_port = model
                .plane_wave_acoustic_port(&reordered, area, port_controls(), cx)
                .unwrap();
            assert_eq!(port.weights_m2(), reordered_port.weights_m2());
            for moving in [false, true] {
                let supports: Vec<_> = (0..3 * nodes.len())
                    .map(|dof| {
                        (
                            dof,
                            if moving {
                                scale(displacement, q[dof % 3][2])
                            } else {
                                ZERO
                            },
                        )
                    })
                    .collect();
                let incident = if moving { ZERO } else { C64::new(200.0, -50.0) };
                let response = model
                    .solve_harmonic_with_acoustic_ports(
                        omega,
                        &[HarmonicAcousticPort {
                            port: &port,
                            medium: &medium,
                            incident_pressure_pa: incident,
                        }],
                        None,
                        &supports,
                        &helectrodes(&original, ZERO, HCondition::Current(ZERO)),
                        &[],
                        &vec![[ZERO; 3]; nodes.len()],
                        &vec![ZERO; nodes.len()],
                        hc(),
                        cx,
                    )
                    .unwrap();
                let a = &response.acoustic_ports[0];
                if moving {
                    cnear(a.volume_velocity_m3_s, expected_flow, 1e-22, 5e-12);
                    near(
                        response.average_acoustic_power_w,
                        expected_power,
                        1e-20,
                        5e-11,
                    );
                    near(
                        response.average_mechanical_power_w,
                        expected_power,
                        1e-18,
                        5e-8,
                    );
                } else {
                    cnear(a.volume_velocity_m3_s, ZERO, 1e-30, 0.0);
                    cnear(a.outgoing_pressure_pa, incident, 1e-12, 1e-14);
                    cnear(a.pressure_pa, scale(incident, 2.0), 1e-12, 1e-14);
                    near(
                        a.average_outgoing_power_w,
                        a.average_incident_power_w,
                        1e-25,
                        1e-14,
                    );
                    near(response.average_acoustic_power_w, 0.0, 1e-25, 0.0);
                }
            }
        }
    });
}

#[test]
fn g0_g4_acoustic_ports_refuse_foreign_duplicate_or_out_of_domain_loading() {
    with_cx(|cx, gate| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let materials = [material([0.0, 0.0, 12.0, 0.0, 0.0, 0.0])];
        let model =
            Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let other =
            Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let faces = face_trace(&cube, &cube.top);
        let area = LENGTH * LENGTH;
        let pc = port_controls();
        let port = model
            .plane_wave_acoustic_port(&faces, area, pc, cx)
            .unwrap();
        let medium = FluidAcousticState::try_new(1000.0, 1480.0).unwrap();
        let input = HarmonicAcousticPort {
            port: &port,
            medium: &medium,
            incident_pressure_pa: ZERO,
        };
        let fixed: Vec<_> = (0..3 * n).map(|i| (i, ZERO)).collect();
        assert!(
            model
                .plane_wave_acoustic_port(&[faces[0], faces[0]], area, pc, cx)
                .is_err()
        );
        assert!(
            model
                .plane_wave_acoustic_port(&faces, area / 2.0, pc, cx)
                .is_err()
        );
        assert!(
            model
                .plane_wave_acoustic_port(
                    &faces,
                    area,
                    PlaneWavePortControls {
                        maximum_surface_dofs: 1,
                        ..pc
                    },
                    cx
                )
                .is_err()
        );
        let mut opposed = faces.clone();
        opposed.extend(face_trace(&cube, &cube.bottom));
        assert!(
            model
                .plane_wave_acoustic_port(&opposed, 2.0 * area, pc, cx)
                .is_err()
        );
        let solve = |owner: &Model<'_>, ports: &[HarmonicAcousticPort<'_>], omega| {
            owner.solve_harmonic_with_acoustic_ports(
                omega,
                ports,
                None,
                &fixed,
                &helectrodes(&cube, ZERO, HCondition::Current(ZERO)),
                &[],
                &vec![[ZERO; 3]; n],
                &vec![ZERO; n],
                hc(),
                cx,
            )
        };
        assert!(solve(&other, &[input], 1e5).is_err());
        assert!(solve(&model, &[input, input], 1e5).is_err());
        assert!(
            solve(
                &model,
                &[input],
                2.0 * pc.maximum_wavenumber_rad_m * medium.sound_speed_m_s()
            )
            .is_err()
        );
        // A finite but excessive incoming pressure cannot bypass linear acoustics.
        assert!(
            solve(
                &model,
                &[HarmonicAcousticPort {
                    incident_pressure_pa: C64::new(medium.adiabatic_bulk_modulus_pa(), 0.0),
                    ..input
                }],
                1e5
            )
            .is_err()
        );
        solve(&model, &[input], 1e5).unwrap();
        gate.request();
        assert!(matches!(
            solve(&model, &[input], 1e5),
            Err(TetPiezoelectricError::Cancelled)
        ));
    });
}
