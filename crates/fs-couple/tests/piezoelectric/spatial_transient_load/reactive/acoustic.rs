//! G1/G3/G4: incident sound, native motion and retained R/L/C storage.
use super::*;
use fs_couple::piezoelectric::spatial::transient::acoustic::{
    PlaneWavePortControls, TransientAcousticPort,
};
use fs_material::fluid_acoustics::FluidAcousticState;
use std::collections::BTreeSet;

fn faces(cube: &Cube) -> Vec<[usize; 3]> {
    let mut faces = BTreeSet::new();
    for tet in &cube.tets {
        for omit in 0..4 {
            let nodes: Vec<_> = (0..4).filter(|&i| i != omit).map(|i| tet[i]).collect();
            if nodes.iter().all(|i| cube.top.contains(i)) {
                let mut face = [nodes[0], nodes[1], nodes[2]];
                face.sort_unstable();
                faces.insert(face);
            }
        }
    }
    faces.into_iter().collect()
}

fn port_controls() -> PlaneWavePortControls {
    PlaneWavePortControls {
        maximum_surface_dofs: 64,
        coplanarity_tolerance: 1e-10,
        maximum_wavenumber_rad_m: core::f64::consts::PI / LENGTH,
        maximum_relative_amplitude: 1e-2,
    }
}

fn ac() -> Controls {
    Controls { energy_tolerance_j: 2e-18, ..controls() }
}

#[test]
fn g1_incident_sound_drives_reactive_storage_and_returns_as_radiation_and_heat() {
    with_cx(|cx, _| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let materials = [material([0.0, 0.0, 12.0, 0.0, 0.0, 0.0])];
        let model = Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let trace = faces(&cube);
        let area = LENGTH * LENGTH;
        let port = model.plane_wave_acoustic_port(&trace, area, port_controls(), cx).unwrap();
        // Authored fluid keeps the retained mechanical pole inside the guide band.
        let medium = FluidAcousticState::try_new(1000.0, 14800.0).unwrap();
        let z = medium.characteristic_impedance_pa_s_m() / area;
        let node = *cube.top.last().unwrap();
        let dof = 3 * node + 2;
        let fixed: Vec<_> = (0..3 * n).filter(|&i| i != dof).map(|i| (i, 0.0)).collect();
        let zero = vec![[0.0; 3]; n];
        let charge = vec![0.0; n];
        let body = Body::new(&model, &fixed, &electrodes(&cube, 0.0, Condition::Charge(0.0)),
            &zero, &zero, &charge, ac(), cx).unwrap();
        // Independent P1 triangle integration, not the prepared acoustic weights.
        let b: f64 = trace.iter().filter(|f| f.contains(&node)).map(|f| {
            let [a, b, c] = f.map(|i| cube.nodes[i]);
            ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])).abs() / 6.0
        }).sum();
        let diagonal = |matrix: &fs_sparse::Csr| {
            let (cols, values) = matrix.row(dof);
            values[cols.iter().position(|&i| i == dof).unwrap()]
        };
        let (m, k) = (diagonal(&model.solid().mass), diagonal(&model.solid().stiffness));
        let (cols, values) = model.coupling().row(dof);
        let g: f64 = cols.iter().zip(values).filter(|(i, _)| cube.top.contains(i))
            .map(|(_, v)| v).sum();
        let ceramic_c = EPSILON[2] * LENGTH;
        let cap = capacitor(ceramic_c);
        let coil = inductor(1e-5);
        let series = resistor(100.0, [1, 0]);
        let parallel = resistor(1e4, [1, 0]);
        let (c, l, r, conductance) = (ceramic_c + cap.capacitance_f(), coil.inductance_h(),
            series.conductor.resistance_ohm(), 1.0 / parallel.conductor.resistance_ohm());
        let mut coupled = Shunted::new(body, vec![
            Reactive::SeriesRl { electrodes: [1, 0], inductor: coil, conductor: Some(series.conductor) },
            Reactive::Capacitor { electrodes: [1, 0], capacitor: cap },
        ], vec![State::FluxLinkage(0.0), State::Charge(0.0)], vec![parallel], ac(), cx).unwrap();
        let (dt, h) = (2e-8, 1e-8);
        let (mut x, mut speed, mut voltage, mut flux) = (0.0, 0.0, 0.0, 0.0);
        let mut all_heat = 0.0;
        let mut all_outward = 0.0;
        for incident in [1e5, -5e4, 0.0, 0.0] {
            // Eliminate the one mechanical mode, electrode voltage and inductor
            // current analytically, retaining acoustic reaction and BOTH loads.
            let electric = c / h + conductance + h / (l + h * r);
            let rhs_electric = c * voltage / h - flux / (l + h * r);
            let midpoint_speed = (m * speed / h - k * x - 2.0 * b * incident
                - g * rhs_electric / electric)
                / (m / h + h * k + z * b * b + g * g / electric);
            let midpoint_voltage = (rhs_electric + g * midpoint_speed) / electric;
            let current = (flux + h * midpoint_voltage) / (l + h * r);
            let (xn, un, vn, ln) = (x + dt * midpoint_speed, 2.0 * midpoint_speed - speed,
                2.0 * midpoint_voltage - voltage, 2.0 * l * current - flux);
            let flow = b * midpoint_speed;
            let pressure = 2.0 * incident + z * flow;
            let heat = dt * (r * current * current + conductance * midpoint_voltage.powi(2));
            let outward = dt * pressure * flow;
            let before_time = coupled.body().state().time_s();
            let before_states = coupled.states().to_vec();
            let input = TransientAcousticPort { port: &port, medium: &medium,
                midpoint_incident_pressure_pa: incident,
                maximum_angular_frequency_rad_s: 0.5 * medium.sound_speed_m_s()
                    * port_controls().maximum_wavenumber_rad_m };
            let candidate = coupled.step_with_acoustic_ports_candidate(dt, &[input], None,
                &[Drive::Voltage(0.0), Drive::Current(0.0)], &zero, &charge, ac(), cx).unwrap();
            let state = candidate.model().body().state();
            let report = candidate.step();
            let sound = &report.acoustic_ports[0];
            near(state.displacement_m()[node][2], xn, 1e-18, 3e-8);
            near(state.velocity_m_s()[node][2], un, 1e-10, 3e-8);
            near(state.electrodes()[1].voltage_v, vn, 1e-10, 3e-8);
            near(state.electrodes()[1].charge_c, ceramic_c * vn - g * xn, 2e-18, 3e-8);
            near(coordinate(candidate.model().states()[0]), ln, 1e-20, 3e-8);
            near(coordinate(candidate.model().states()[1]), (c - ceramic_c) * vn, 1e-20, 3e-8);
            near(report.reactive[0].midpoint_current_a, current, 1e-12, 3e-8);
            near(report.reactive_energy_j[1], 0.5 * ln * ln / l + 0.5 * (c - ceramic_c) * vn * vn,
                2e-20, 3e-8);
            near(report.joule_heat_j, heat, 2e-20, 3e-8);
            near(sound.midpoint_volume_velocity_m3_s, flow, 1e-18, 3e-8);
            near(sound.midpoint_pressure_pa, pressure, 1e-7, 3e-8);
            near(sound.net_outward_work_j, outward, 2e-20, 3e-8);
            near(report.body.acoustic_outward_work_j, outward, 2e-20, 3e-8);
            near(report.energy_residual_j, 0.0, 2e-18, 0.0);
            near(report.electrical_source_work_j, 0.0, 1e-22, 0.0);
            assert_eq!(report.body.mechanical_work_j, 0.0);
            assert_eq!(report.body.material_heat_j, 0.0);
            assert!(report.joule_heat_j > 0.0);
            if incident == 0.0 { assert!(sound.outgoing_energy_j > 0.0); }
            all_heat += heat;
            all_outward += outward;
            near(state.total_energy_j() + report.reactive_energy_j[1] + all_heat + all_outward,
                0.0, 2e-18, 0.0);
            assert_eq!(coupled.body().state().time_s(), before_time);
            assert_eq!(coupled.states(), before_states);
            coupled = candidate.into_parts().0;
            (x, speed, voltage, flux) = (xn, un, vn, ln);
        }
        assert!(all_heat > 0.0 && all_outward < 0.0);
    });
}

#[test]
fn g3_g4_acoustic_shunt_preserves_empty_port_behavior_and_complete_moving_support_work() {
    with_cx(|cx, gate| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let materials = [material([0.0, 0.0, 12.0, 0.0, 0.0, 0.0])];
        let model = Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let area = LENGTH * LENGTH;
        let port = model.plane_wave_acoustic_port(&faces(&cube), area, port_controls(), cx).unwrap();
        let medium = FluidAcousticState::try_new(1000.0, 1480.0).unwrap();
        let fixed: Vec<_> = (0..3 * n).map(|i| (i, 0.0)).collect();
        let supports: Vec<_> = (0..3 * n).map(|i| (i, if i % 3 == 2 { 0.01 } else { 0.0 })).collect();
        let zero = vec![[0.0; 3]; n];
        let charge = vec![0.0; n];
        let voltage = 3.0;
        let c = EPSILON[2] * LENGTH;
        let body = Body::new(&model, &fixed, &electrodes(&cube, 0.0, Condition::Charge(c * voltage)),
            &zero, &zero, &charge, ac(), cx).unwrap();
        let coupled = Shunted::new(body, vec![
            Reactive::SeriesRl { electrodes: [1, 0], inductor: inductor(2e-3),
                conductor: Some(resistor(100.0, [1, 0]).conductor) },
            Reactive::Capacitor { electrodes: [1, 0], capacitor: capacitor(c) },
        ], vec![State::FluxLinkage(0.0), State::Charge(c * voltage)], vec![], ac(), cx).unwrap();
        let drives = [Drive::Voltage(0.0), Drive::Current(0.0)];
        let dt = 2e-7;
        let ordinary = coupled.step_candidate(dt, &drives, &zero, &charge, ac(), cx).unwrap();
        let empty = coupled.step_with_acoustic_ports_candidate(dt, &[], None,
            &drives, &zero, &charge, ac(), cx).unwrap();
        assert_eq!(ordinary.model().states(), empty.model().states());
        assert_eq!(ordinary.model().body().state().displacement_m(), empty.model().body().state().displacement_m());
        assert_eq!(ordinary.model().body().state().velocity_m_s(), empty.model().body().state().velocity_m_s());
        assert_eq!(ordinary.step().energy_residual_j, empty.step().energy_residual_j);
        assert!(empty.step().acoustic_ports.is_empty());
        let input = TransientAcousticPort { port: &port, medium: &medium,
            midpoint_incident_pressure_pa: 0.0,
            maximum_angular_frequency_rad_s: 0.5 * medium.sound_speed_m_s()
                * port_controls().maximum_wavenumber_rad_m };
        let baseline = coupled.step_with_support_motion_candidate(dt, &supports,
            &drives, &zero, &charge, ac(), cx).unwrap();
        let acoustic = coupled.step_with_acoustic_ports_candidate(dt, &[input], Some(&supports),
            &drives, &zero, &charge, ac(), cx).unwrap();
        let flow = area * 0.005;
        let pressure = medium.characteristic_impedance_pa_s_m() * 0.005;
        let emitted = dt * pressure * flow;
        let kinetic = 0.5 * 7500.0 * LENGTH.powi(3) * 0.01_f64.powi(2);
        near(acoustic.step().body.support_work_j, kinetic + emitted, 2e-18, 3e-9);
        near(acoustic.step().body.acoustic_outward_work_j, emitted, 2e-18, 3e-9);
        near(acoustic.step().body.support_work_j - baseline.step().body.support_work_j, emitted, 2e-18, 3e-9);
        near(acoustic.step().acoustic_ports[0].outgoing_energy_j, emitted, 2e-18, 3e-9);
        assert_eq!(acoustic.model().states(), baseline.model().states());
        near(acoustic.step().energy_residual_j, 0.0, 2e-18, 0.0);
        assert!(acoustic.step().joule_heat_j > 0.0);
        assert!(coupled.step_with_acoustic_ports_candidate(dt, &[input, input], Some(&supports),
            &drives, &zero, &charge, ac(), cx).is_err());
        gate.request();
        assert!(matches!(coupled.step_with_acoustic_ports_candidate(dt, &[input], Some(&supports),
            &drives, &zero, &charge, ac(), cx), Err(TetPiezoelectricError::Cancelled)));
        assert_eq!(coupled.body().state().time_s(), 0.0);
        assert_eq!(coupled.states(), &[State::FluxLinkage(0.0), State::Charge(c * voltage)]);
    });
}
