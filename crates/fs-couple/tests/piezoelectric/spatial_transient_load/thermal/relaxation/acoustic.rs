//! G1/G3: actual sound work, intrinsic memory and temperature-dependent load heat.
use super::*;
use fs_couple::piezoelectric::spatial::transient::acoustic::{
    PlaneWaveAcousticPort, PlaneWavePortControls, TransientAcousticPort,
};
use fs_material::fluid_acoustics::FluidAcousticState;
use std::collections::BTreeSet;

fn faces(cube: &Cube) -> Vec<[usize; 3]> {
    let mut result = BTreeSet::new();
    for tet in &cube.tets {
        for omit in 0..4 {
            let face: Vec<_> = (0..4).filter(|&i| i != omit).map(|i| tet[i]).collect();
            if face.iter().all(|i| cube.top.contains(i)) {
                let mut face = [face[0], face[1], face[2]];
                face.sort_unstable();
                result.insert(face);
            }
        }
    }
    result.into_iter().collect()
}
fn port_controls() -> PlaneWavePortControls {
    PlaneWavePortControls {
        maximum_surface_dofs: 64,
        coplanarity_tolerance: 1e-10,
        maximum_wavenumber_rad_m: core::f64::consts::PI / LENGTH,
        maximum_relative_amplitude: 1e-2,
    }
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
fn g1_sound_drives_motion_memory_and_distinct_finite_heat_cells_with_resistor_feedback() {
    with_cx(|cx, _| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let (card, curve, base, memory) =
            ceramic_with_mechanical_memory(false, Some((0.5, 0.7e-6)));
        let materials = [base];
        let model = Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let trace = faces(&cube);
        let port = model.plane_wave_acoustic_port(&trace, AREA, port_controls(), cx).unwrap();
        let medium = FluidAcousticState::try_new(1000.0, 14800.0).unwrap();
        let z = medium.characteristic_impedance_pa_s_m() / AREA;
        let node = *cube.top.last().unwrap();
        let dof = 3 * node + 2;
        let fixed: Vec<_> = (0..3 * n).filter(|&i| i != dof).map(|i| (i, 0.0)).collect();
        let (zero, charge) = (vec![[0.0; 3]; n], vec![0.0; n]);
        let mut body = Body::new(
            &model, &fixed, &electrodes(&cube, 0.0, Condition::Charge(0.0)),
            &zero, &zero, &charge, field_controls(), cx,
        ).unwrap().with_relaxation(vec![memory], MEMORY_BUDGET, field_controls(), cx).unwrap();
        let alpha = 0.02;
        let (load_card, load_curve) = thermal_material(2.0, 4.0, alpha, 400.0);
        let cards = [card, load_card];
        let curves = [curve, load_curve];
        let context = point(300.0);
        let cells = caloric_cells(&cards, &curves, &context);
        let resistors = [CellResistor {
            destination: ResistorCell { branch: LOAD.resistor_branch, cell: CELLS[1] },
            selection: MaterialPropertySelection::SingleClaimOnly,
            length_m: LENGTH,
            area_m2: AREA,
        }];
        let network = MaterialNetwork { cells: &cells, resistors: &resistors, links: &[] };
        let mut thermal = isolated_cells(&curves);
        let prepared = network.prepare(&thermal, || false).unwrap();
        let routes = ceramic_routes(cube.tets.len());
        let m = diagonal(&model.solid().mass, dof);
        let k = diagonal(&model.solid().stiffness, dof);
        let (cols, values) = model.coupling().row(dof);
        let g: f64 = cols.iter().zip(values).filter(|(i, _)| cube.top.contains(i))
            .map(|(_, v)| v).sum();
        // Independent triangular P1 shape integral, not the native port weights.
        let b: f64 = trace.iter().filter(|face| face.contains(&node)).map(|face| {
            let [a, b, c] = face.map(|i| cube.nodes[i]);
            ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])).abs() / 6.0
        }).sum();
        let (cp, dk, tau_m, dt) = (POLE_INCREMENT_F_M * LENGTH, 0.5 * k, 0.7e-6, 2e-7);
        let h = 0.5 * dt;
        let (mut x, mut v, mut voltage, mut polarization, mut maxwell) = (0.0, 0.0, 0.0, 0.0, 0.0);
        let (mut ceramic_heat, mut load_heat, mut sound_work) = (0.0, 0.0, 0.0);
        let initial_energy = thermal.energy_j() + body.state().total_energy_j();
        for incident in [1e5, -5e4, 0.0, 0.0] {
            let old_temperature = 300.0 + load_heat / (2.0 * VOLUME);
            // Independently solve one mechanical/electrical coordinate and
            // two memory ODEs, then the scalar finite-resistor heat equation.
            let midpoint = |temperature: f64| {
                let resistance = RESISTANCE * (1.0 + alpha * (temperature - 300.0));
                let electric = CAPACITANCE / h + 1.0 / resistance + cp / (POLE_TIME_S + h);
                let electric_rhs = CAPACITANCE * voltage / h + polarization / (POLE_TIME_S + h);
                let vm = (m * v / h - k * x - 2.0 * b * incident
                    - tau_m * maxwell / (tau_m + h) - g * electric_rhs / electric)
                    / (m / h + k * h + z * b * b + h * dk * tau_m / (tau_m + h)
                        + g * g / electric);
                let em = (electric_rhs + g * vm) / electric;
                let pm = (POLE_TIME_S * polarization + h * cp * em) / (POLE_TIME_S + h);
                let sm = (tau_m * maxwell + h * dk * tau_m * vm) / (tau_m + h);
                [vm, em, pm, sm, dt * em * em / resistance]
            };
            let (mut low, mut high) = (old_temperature, 400.0);
            for _ in 0..64 {
                let t = 0.5 * (low + high);
                if t - old_temperature > midpoint(t)[4] / (2.0 * VOLUME) { high = t; }
                else { low = t; }
            }
            let [vm, em, pm, sm, joule] = midpoint(0.5 * (low + high));
            let (xn, vn, en, pn, sn) =
                (x + dt * vm, 2.0 * vm - v, 2.0 * em - voltage,
                 2.0 * pm - polarization, 2.0 * sm - maxwell);
            let heat = dt * sm * sm / (dk * tau_m)
                + POLE_TIME_S * (pn - polarization).powi(2) / (cp * dt);
            let outward = dt * (2.0 * incident + z * b * vm) * b * vm;
            let before_time = body.state().time_s();
            let before_heat = thermal.energy_j();
            let candidate = body.step_with_acoustic_thermal_loads_candidate(
                &thermal, thermal.time_s() + dt, &[wave(&port, &medium, incident)], None,
                &prepared, &[LOAD], &routes, "T", &[Drive::Voltage(0.0), Drive::Current(0.0)],
                &zero, &charge, field_controls(), thermal_controls(), &[0.0; 2], cx,
            ).unwrap();
            let report = candidate.report();
            near(candidate.body().state().displacement_m()[node][2], xn, 1e-17, 3e-8);
            near(candidate.body().state().velocity_m_s()[node][2], vn, 1e-10, 3e-8);
            near(candidate.body().state().electrodes()[1].voltage_v, en, 1e-8, 3e-8);
            near(candidate.body().state().electrodes()[1].charge_c,
                CAPACITANCE * en + pn - g * xn, 2e-18, 3e-8);
            near(report.electrical.native.material_heat_j, heat, 2e-18, 3e-8);
            near(report.electrical.loads.joule_heat_j, joule, 2e-18, 3e-8);
            near(report.electrical.native.acoustic_outward_work_j, outward, 2e-18, 3e-8);
            near(report.electrical.acoustic_ports[0].net_outward_work_j, outward, 2e-18, 3e-8);
            near(report.deposited_joule_j, heat + joule, 2e-18, 3e-8);
            near(report.non_electrical_heat_j, 0.0, 2e-16, 0.0);
            near(report.total_energy_residual_j, 0.0, 2e-16, 0.0);
            assert_eq!(report.electrical.native.mechanical_work_j, 0.0);
            assert_eq!(report.electrical.loads.electrical_source_work_j, 0.0);
            ceramic_heat += heat;
            load_heat += joule;
            sound_work += outward;
            near(candidate.thermal().temperature_k(CELLS[0]).unwrap(),
                300.0 + ceramic_heat / (CERAMIC_DENSITY_KG_M3 * VOLUME * CERAMIC_CP_J_KG_K),
                2e-8, 0.0);
            near(candidate.thermal().temperature_k(CELLS[1]).unwrap(),
                300.0 + load_heat / (2.0 * VOLUME), 2e-8, 0.0);
            near(candidate.thermal().energy_j() + candidate.body().state().total_energy_j()
                + sound_work, initial_energy, 4e-16, 0.0);
            assert_eq!(body.state().time_s(), before_time);
            assert_eq!(thermal.energy_j(), before_heat);
            (thermal, body, _) = candidate.into_parts();
            (x, v, voltage, polarization, maxwell) = (xn, vn, en, pn, sn);
        }
        assert!(ceramic_heat > 0.0 && load_heat > 0.0);
    });
}

#[test]
fn g3_acoustic_support_work_is_not_deposited_as_heat_and_refusals_leave_both_histories() {
    with_cx(|cx, gate| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let (card, curve, base, memory) = ceramic(false);
        let materials = [base];
        let model = Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let trace = faces(&cube);
        let port = model.plane_wave_acoustic_port(&trace, AREA, port_controls(), cx).unwrap();
        let medium = FluidAcousticState::try_new(1000.0, 1480.0).unwrap();
        let fixed: Vec<_> = (0..3 * n).map(|i| (i, 0.0)).collect();
        let (zero, charge) = (vec![[0.0; 3]; n], vec![0.0; n]);
        let body = Body::new(
            &model, &fixed, &electrodes(&cube, 0.0, Condition::Charge(0.0)),
            &zero, &zero, &charge, field_controls(), cx,
        ).unwrap().with_relaxation(vec![memory], MEMORY_BUDGET, field_controls(), cx).unwrap();
        let (other_card, other_curve) = thermal_material(2.0, 4.0, 0.0, 400.0);
        let cards = [card, other_card];
        let curves = [curve, other_curve];
        let context = point(300.0);
        let cells = caloric_cells(&cards, &curves, &context);
        let network = MaterialNetwork { cells: &cells, resistors: &[], links: &[] };
        let thermal = isolated_cells(&curves);
        let prepared = network.prepare(&thermal, || false).unwrap();
        let routes = ceramic_routes(cube.tets.len());
        let dt = 2e-7;
        let speed = 0.02;
        let support: Vec<_> = (0..3 * n).map(|i| (i, if i % 3 == 2 { speed } else { 0.0 })).collect();
        let drives = [Drive::Voltage(0.0), Drive::Current(0.0)];
        let solve = |ports: &[TransientAcousticPort<'_>]| {
            body.step_with_acoustic_thermal_loads_candidate(
                &thermal, dt, ports, Some(&support), &prepared, &[], &routes, "T", &drives,
                &zero, &charge, field_controls(), thermal_controls(), &[0.0; 2], cx,
            )
        };
        let input = wave(&port, &medium, 0.0);
        let candidate = solve(&[input]).unwrap();
        let outward = dt * medium.characteristic_impedance_pa_s_m() * AREA * (speed * 0.5).powi(2);
        let kinetic = 0.5 * CERAMIC_DENSITY_KG_M3 * VOLUME * speed * speed;
        near(candidate.report().electrical.native.support_work_j, kinetic + outward, 2e-18, 3e-8);
        near(candidate.report().electrical.native.acoustic_outward_work_j, outward, 2e-18, 3e-8);
        near(candidate.report().deposited_joule_j, 0.0, 2e-18, 0.0);
        near(candidate.thermal().energy_j(), thermal.energy_j(), 2e-16, 0.0);
        let empty = solve(&[]).unwrap();
        let original = body.step_with_material_thermal_support_motion_candidate(
            &thermal, dt, &support, &prepared, &[], &routes, "T", &drives,
            &zero, &charge, field_controls(), thermal_controls(), &[0.0; 2], cx,
        ).unwrap();
        assert!(empty.report().electrical.acoustic_ports.is_empty());
        assert_eq!(empty.body().state().displacement_m(), original.body().state().displacement_m());
        assert_eq!(empty.thermal().energy_j(), original.thermal().energy_j());
        assert!(matches!(solve(&[input, input]), Err(CircuitEnthalpyError::Source(
            LoadThermalError::Piezoelectric(TetPiezoelectricError::Invalid(_))))));
        gate.request();
        assert!(matches!(solve(&[input]), Err(CircuitEnthalpyError::Cancelled)));
        assert_eq!(body.state().time_s(), 0.0);
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(body.state().material_memory_energy_j(), 0.0);
    });
}
