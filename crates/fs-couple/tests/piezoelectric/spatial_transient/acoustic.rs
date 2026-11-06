//! G1/G3/G4 time-domain wave/solid/charge work with authored material coefficients.
use super::*;
use fs_couple::piezoelectric::spatial::transient::acoustic::{
    PlaneWaveAcousticPort, PlaneWavePortControls, TransientAcousticPort,
};
use fs_couple::piezoelectric::spatial::transient::load::ElectrodeDrive as Drive;
use fs_material::fluid_acoustics::FluidAcousticState;
use std::collections::BTreeSet;

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
fn faces(cube: &Cube) -> Vec<[usize; 3]> {
    let mut faces = BTreeSet::new();
    for tet in &cube.tets {
        for omit in 0..4 {
            let face: Vec<_> = (0..4).filter(|&i| i != omit).map(|i| tet[i]).collect();
            if face.iter().all(|i| cube.top.contains(i)) {
                let mut face = [face[0], face[1], face[2]];
                face.sort_unstable();
                faces.insert(face);
            }
        }
    }
    faces.into_iter().collect()
}
fn wave<'a>(
    port: &'a PlaneWaveAcousticPort<'a>,
    medium: &'a FluidAcousticState,
    pressure: f64,
) -> TransientAcousticPort<'a> {
    TransientAcousticPort {
        port,
        medium,
        midpoint_incident_pressure_pa: pressure,
        maximum_angular_frequency_rad_s:
            0.5 * medium.sound_speed_m_s() * port_controls().maximum_wavenumber_rad_m,
    }
}
fn diagonal(matrix: &fs_sparse::Csr, dof: usize) -> f64 {
    let (cols, values) = matrix.row(dof);
    values[cols.iter().position(|&i| i == dof).unwrap()]
}

#[test]
fn g1_incident_sound_couples_motion_shunt_heat_and_intrinsic_memory_without_lost_wave_work() {
    for memory in [false, true] {
        with_cx(|cx, _| {
            let cube = cube(1);
            let n = cube.nodes.len();
            let (base, loss) = crate::harmonic::relaxation::fixture(
                Some((0.5, 0.7e-6)), Some((0.8 * EPSILON[2], 1e-6)), FRAME,
            );
            let materials = [base];
            let model = Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx)
                .unwrap();
            let trace = faces(&cube);
            let port = model.plane_wave_acoustic_port(
                &trace, LENGTH * LENGTH, port_controls(), cx,
            ).unwrap();
            // Authored fast-fluid coefficients keep the retained mechanical
            // pole below the declared square-guide transverse cutoff.
            let medium = FluidAcousticState::try_new(1000.0, 14800.0).unwrap();
            let node = *cube.top.last().unwrap();
            let dof = 3 * node + 2;
            let fixed: Vec<_> = (0..3 * n).filter(|&i| i != dof).map(|i| (i, 0.0)).collect();
            let (zero, charge) = zero_source(n);
            let mut body = Body::new(
                &model, &fixed, &electrodes(&cube, 0.0, Condition::Charge(0.0)),
                &zero, &zero, &charge, ac(), cx,
            ).unwrap();
            if memory {
                body = body.with_relaxation(vec![loss], 256, ac(), cx).unwrap();
            }
            // Independent integral of the selected vertex's P1 shape over
            // horizontal triangles; no native acoustic weight enters the oracle.
            let b: f64 = trace.iter().filter(|f| f.contains(&node)).map(|f| {
                let [a, b, c] = f.map(|i| cube.nodes[i]);
                ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])).abs() / 6.0
            }).sum();
            let k = diagonal(&model.solid().stiffness, dof);
            let m = diagonal(&model.solid().mass, dof);
            let (cols, values) = model.coupling().row(dof);
            let g: f64 = cols.iter().zip(values).filter(|(i, _)| cube.top.contains(i))
                .map(|(_, v)| v).sum();
            let c = EPSILON[2] * LENGTH;
            let cp = if memory { 0.8 * c } else { 0.0 };
            let dk = if memory { 0.5 * k } else { 0.0 };
            let (tau_m, tau_p, dt, resistance) = (0.7e-6, 1e-6, 2e-7, 1e4);
            let h = 0.5 * dt;
            let z = medium.characteristic_impedance_pa_s_m() / (LENGTH * LENGTH);
            let loads = [crate::transient_load::resistor(resistance, [1, 0])];
            let (mut x, mut v, mut voltage, mut polarization, mut maxwell_force) =
                (0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64);
            for incident in [1e5, -5e4, 0.0, 0.0] {
                // Eliminate one original mechanical coordinate, electrode,
                // Maxwell force and Debye charge; retain actual guide damping.
                let electric = c / h + 1.0 / resistance + cp / (tau_p + h);
                let electric_rhs = c * voltage / h + polarization / (tau_p + h);
                let vm = (m * v / h - k * x - 2.0 * b * incident
                    - tau_m * maxwell_force / (tau_m + h) - g * electric_rhs / electric)
                    / (m / h + k * h + z * b * b + h * dk * tau_m / (tau_m + h)
                        + g * g / electric);
                let mid_voltage = (electric_rhs + g * vm) / electric;
                let pm = (tau_p * polarization + h * cp * mid_voltage) / (tau_p + h);
                let sm = (tau_m * maxwell_force + h * dk * tau_m * vm) / (tau_m + h);
                let (xn, vn) = (x + dt * vm, 2.0 * vm - v);
                let (en, pn, sn) = (2.0 * mid_voltage - voltage, 2.0 * pm - polarization,
                    2.0 * sm - maxwell_force);
                let ceramic_heat = if memory {
                    dt * sm * sm / (dk * tau_m)
                        + tau_p * (pn - polarization).powi(2) / (cp * dt)
                } else { 0.0 };
                let memory_energy = if memory {
                    0.5 * sn * sn / dk + 0.5 * pn * pn / cp
                } else { 0.0 };
                let resistor_heat = dt * mid_voltage * mid_voltage / resistance;
                let flow = b * vm;
                let outgoing = incident + z * flow;
                let pressure = outgoing + incident;
                let outward = dt * pressure * flow;
                let before = body.state().total_energy_j();
                let before_time = body.state().time_s();
                let (candidate, loaded) = body.step_with_acoustic_loads_candidate(
                    dt, &[wave(&port, &medium, incident)], None,
                    &[Drive::Voltage(0.0), Drive::Current(0.0)], &loads,
                    &zero, &charge, ac(), cx,
                ).unwrap();
                let sound = &candidate.acoustic_ports()[0];
                let step = candidate.step();
                near(candidate.body().state().displacement_m()[node][2], xn, 1e-17, 3e-8);
                near(candidate.body().state().velocity_m_s()[node][2], vn, 1e-10, 3e-8);
                near(candidate.body().state().electrodes()[1].voltage_v, en, 1e-9, 3e-8);
                near(candidate.body().state().electrodes()[1].charge_c, c * en + pn - g * xn, 2e-18, 3e-8);
                near(candidate.body().state().total_energy_j(),
                    0.5 * m * vn * vn + 0.5 * k * xn * xn + 0.5 * c * en * en + memory_energy,
                    2e-20, 3e-8);
                near(loaded.joule_heat_j, resistor_heat, 2e-20, 3e-8);
                near(step.material_heat_j, ceramic_heat, 2e-20, 3e-8);
                near(sound.midpoint_volume_velocity_m3_s, flow, 1e-18, 3e-8);
                near(sound.midpoint_pressure_pa, pressure, 1e-7, 3e-8);
                near(sound.midpoint_outgoing_pressure_pa, outgoing, 1e-7, 3e-8);
                near(sound.incident_energy_j, dt * incident * incident / z, 2e-20, 3e-8);
                near(sound.outgoing_energy_j, dt * outgoing * outgoing / z, 2e-20, 3e-8);
                near(sound.net_outward_work_j, outward, 2e-20, 3e-8);
                near(step.acoustic_outward_work_j, outward, 2e-20, 3e-8);
                near(sound.wave_energy_residual_j, 0.0, 2e-18, 0.0);
                near(step.energy_residual_j, 0.0, 2e-18, 0.0);
                near(loaded.energy_residual_j, 0.0, 2e-18, 0.0);
                assert_eq!(step.mechanical_work_j, 0.0);
                assert_eq!(loaded.electrical_source_work_j, 0.0);
                if incident == 0.0 {
                    assert!(sound.outgoing_energy_j > 0.0);
                    assert!(candidate.body().state().total_energy_j() < before);
                }
                if !memory { assert_eq!(step.material_heat_j, 0.0); }
                assert_eq!(body.state().time_s(), before_time);
                assert_eq!(body.state().total_energy_j(), before);
                body = candidate.into_parts().0;
                (x, v, voltage, polarization, maxwell_force) = (xn, vn, en, pn, sn);
            }
        });
    }
}

#[test]
fn g3_rotated_moving_interface_pays_exact_inertial_and_outgoing_wave_work() {
    with_cx(|cx, _| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let trace = faces(&cube);
        let area = LENGTH * LENGTH;
        let mass = 7500.0 * LENGTH.powi(3);
        let medium = FluidAcousticState::try_new(1000.0, 1480.0).unwrap();
        let z = medium.characteristic_impedance_pa_s_m() / area;
        let dt = 2e-7;
        for angle in [0.0, 0.61] {
            let q = rotation([1.0 / 3.0_f64.sqrt(); 3], angle);
            let frame = ContentHash([109; 32]);
            let nodes: Vec<_> = cube.nodes.iter().map(|p| apply(&q, p)).collect();
            let materials = [material([0.0, 0.0, 12.0, 0.0, 0.0, 0.0]).rotated(frame, q).unwrap()];
            let model = Model::assemble(&nodes, &cube.tets, &materials, frame, budget(), cx).unwrap();
            let port = model.plane_wave_acoustic_port(&trace, area, port_controls(), cx).unwrap();
            let fixed: Vec<_> = (0..3 * n).map(|i| (i, 0.0)).collect();
            let (zero, charge) = zero_source(n);
            let conditions = [Condition::Voltage(0.0), Condition::Charge(0.0)];
            let mut body = Body::new(
                &model, &fixed, &electrodes(&cube, 0.0, Condition::Charge(0.0)),
                &zero, &zero, &charge, ac(), cx,
            ).unwrap();
            let mut old_velocity = 0.0;
            for velocity in [0.01, 0.02, 0.0] {
                let supports: Vec<_> = (0..3 * n)
                    .map(|dof| (dof, velocity * q[dof % 3][2])).collect();
                let vm = 0.5 * (old_velocity + velocity);
                let flow = area * vm;
                let pressure = z * flow;
                let emitted = dt * z * flow * flow;
                let kinetic_change = 0.5 * mass * (velocity * velocity - old_velocity * old_velocity);
                let candidate = body.step_with_acoustic_support_motion_candidate(
                    dt, &[wave(&port, &medium, 0.0)], &supports, &conditions,
                    &zero, &charge, ac(), cx,
                ).unwrap();
                let reaction: f64 = candidate.midpoint_support_force_n().iter()
                    .map(|f| (0..3).map(|axis| f[axis] * q[axis][2]).sum::<f64>()).sum();
                near(reaction, mass * (velocity - old_velocity) / dt + area * pressure, 1e-9, 3e-9);
                near(candidate.step().support_work_j, kinetic_change + emitted, 2e-20, 3e-9);
                near(candidate.step().acoustic_outward_work_j, emitted, 2e-20, 3e-9);
                near(candidate.acoustic_ports()[0].outgoing_energy_j, emitted, 2e-20, 3e-9);
                near(candidate.body().state().kinetic_energy_j(), 0.5 * mass * velocity * velocity, 2e-20, 3e-9);
                assert_eq!(candidate.step().material_heat_j, 0.0);
                assert_eq!(candidate.step().applied_force_work_j, 0.0);
                body = candidate.into_parts().0;
                old_velocity = velocity;
            }
            // A stopped interface reflects pressure with unit amplitude and
            // double total pressure, but supplies no mechanical work or heat.
            let reflected = body.step_with_acoustic_ports_candidate(
                dt, &[wave(&port, &medium, 1250.0)], &conditions,
                &zero, &charge, ac(), cx,
            ).unwrap();
            let sound = &reflected.acoustic_ports()[0];
            near(sound.midpoint_pressure_pa, 2500.0, 1e-10, 0.0);
            near(sound.midpoint_outgoing_pressure_pa, 1250.0, 1e-10, 0.0);
            assert_eq!(sound.incident_energy_j, sound.outgoing_energy_j);
            assert_eq!(sound.net_outward_work_j, 0.0);
            assert_eq!(reflected.step().mechanical_work_j, 0.0);
        }
    });
}

#[test]
fn g4_acoustic_source_domains_and_cancellation_leave_the_accepted_body_untouched() {
    with_cx(|cx, gate| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let materials = [material([0.0, 0.0, 12.0, 0.0, 0.0, 0.0])];
        let model = Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let other = Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let trace = faces(&cube);
        let area = LENGTH * LENGTH;
        let port = model.plane_wave_acoustic_port(&trace, area, port_controls(), cx).unwrap();
        let foreign = other.plane_wave_acoustic_port(&trace, area, port_controls(), cx).unwrap();
        let medium = FluidAcousticState::try_new(1000.0, 1480.0).unwrap();
        let fixed: Vec<_> = (0..3 * n).map(|i| (i, 0.0)).collect();
        let (zero, charge) = zero_source(n);
        let body = Body::new(
            &model, &fixed, &electrodes(&cube, 0.0, Condition::Charge(0.0)),
            &zero, &zero, &charge, ac(), cx,
        ).unwrap();
        let run = |ports: &[TransientAcousticPort<'_>]| body.step_with_acoustic_ports_candidate(
            2e-7, ports, &[Condition::Voltage(0.0), Condition::Charge(0.0)],
            &zero, &charge, ac(), cx,
        );
        let input = wave(&port, &medium, 1000.0);
        assert!(run(&[wave(&foreign, &medium, 1000.0)]).is_err());
        assert!(run(&[input, input]).is_err());
        assert!(run(&[TransientAcousticPort {
            maximum_angular_frequency_rad_s: 4.0 * input.maximum_angular_frequency_rad_s,
            ..input
        }]).is_err());
        assert!(run(&[TransientAcousticPort {
            midpoint_incident_pressure_pa: medium.adiabatic_bulk_modulus_pa(),
            ..input
        }]).is_err());
        run(&[input]).unwrap();
        assert_eq!(body.state().time_s(), 0.0);
        assert_eq!(body.state().total_energy_j(), 0.0);
        gate.request();
        assert!(matches!(run(&[input]), Err(TetPiezoelectricError::Cancelled)));
        assert_eq!(body.state().time_s(), 0.0);
        assert_eq!(body.state().displacement_m(), zero.as_slice());
    });
}
