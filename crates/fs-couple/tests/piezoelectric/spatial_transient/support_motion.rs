//! G1/G3/G4 imposed-support motion on the original coupled field trajectory.
use super::*;
use fs_couple::piezoelectric::spatial::transient::load::ElectrodeDrive;
use fs_couple::piezoelectric::spatial::transient::load::reactive::{
    ReactiveElectrodeLoad, ReactiveElectrodeState, ReactivePiezoelectricTransient,
};

#[test]
fn g1_rigid_support_acceleration_retains_consistent_inertia_and_reaction_work() {
    with_cx(|cx, _| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let materials = [material([0.0, 0.0, 12.0, 0.0, 0.0, 0.0])];
        let model =
            Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let (zero, charge) = zero_source(n);
        let fixed = (0..3 * n).map(|i| (i, 0.0)).collect::<Vec<_>>();
        let body = Body::new(
            &model,
            &fixed,
            &electrodes(&cube, 0.0, Condition::Voltage(0.0)),
            &zero,
            &zero,
            &charge,
            controls(),
            cx,
        )
        .unwrap();
        let dt = 1e-7_f64;
        let acceleration = [2e4, -1e4, 5e3];
        let supports = (0..3 * n)
            .rev()
            .map(|dof| (dof, dt * acceleration[dof % 3]))
            .collect::<Vec<_>>();
        let mut force = zero.clone();
        force[0] = [0.007, -0.002, 0.003];
        let candidate = body
            .step_with_support_motion_candidate(
                dt,
                &supports,
                &[Condition::Voltage(0.0); 2],
                &force,
                &charge,
                controls(),
                cx,
            )
            .unwrap();
        let state = candidate.body().state();
        let cell_mass = 7500.0 * LENGTH.powi(3) / 6.0;
        for node in 0..n {
            // Independent integral N_i dV for a uniform rigid acceleration.
            let nodal_mass = 0.25
                * cell_mass
                * cube.tets.iter().filter(|tet| tet.contains(&node)).count() as f64;
            for axis in 0..3 {
                near(
                    state.displacement_m()[node][axis],
                    0.5 * dt * dt * acceleration[axis],
                    1e-22,
                    1e-12,
                );
                assert_eq!(state.velocity_m_s()[node][axis], dt * acceleration[axis]);
                near(
                    candidate.midpoint_support_force_n()[node][axis],
                    nodal_mass * acceleration[axis] - force[node][axis],
                    1e-12,
                    1e-10,
                );
            }
        }
        let kinetic = 0.5
            * 7500.0
            * LENGTH.powi(3)
            * acceleration.iter().map(|a| (dt * a).powi(2)).sum::<f64>();
        near(state.kinetic_energy_j(), kinetic, 1e-25, 1e-12);
        assert!(state.stored_energy_j() < 1e-25);
        assert!(candidate.step().applied_force_work_j > 0.0);
        near(candidate.step().mechanical_work_j, kinetic, 1e-25, 1e-12);
        near(
            candidate.step().support_work_j,
            kinetic - candidate.step().applied_force_work_j,
            1e-25,
            1e-12,
        );
        assert!(candidate.step().energy_residual_j.abs() < 1e-23);
        let moving = candidate.into_parts().0;
        assert!(
            moving
                .step_candidate(
                    dt,
                    &[Condition::Voltage(0.0); 2],
                    &zero,
                    &charge,
                    controls(),
                    cx,
                )
                .is_err()
        );
        let stopped = moving
            .step_with_support_motion_candidate(
                dt,
                &fixed,
                &[Condition::Voltage(0.0); 2],
                &zero,
                &charge,
                controls(),
                cx,
            )
            .unwrap();
        assert_eq!(stopped.body().state().velocity_m_s(), &zero);
        near(stopped.step().support_work_j, -kinetic, 1e-24, 1e-11);
        assert_eq!(body.state().displacement_m(), &zero);
        assert_eq!(body.state().time_s(), 0.0);
    });
}

fn entry(matrix: &fs_sparse::Csr, row: usize, col: usize) -> f64 {
    let (columns, values) = matrix.row(row);
    columns
        .iter()
        .position(|&i| i == col)
        .map_or(0.0, |i| values[i])
}

#[test]
fn g1_base_excited_sensor_matches_condensed_scalar_inertia_and_open_charge() {
    with_cx(|cx, _| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let materials = [material([0.0, 0.0, 12.0, 0.0, 0.0, 0.0])];
        let model =
            Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let (zero, charge) = zero_source(n);
        let free = 3 * *cube.top.last().unwrap() + 2;
        let driven = 3 * cube.bottom[0] + 2;
        let fixed = (0..3 * n)
            .filter(|&dof| dof != free)
            .map(|dof| (dof, 0.0))
            .collect::<Vec<_>>();
        let coupling = |dof| {
            cube.top
                .iter()
                .map(|&node| entry(model.coupling(), dof, node))
                .sum::<f64>()
        };
        let (g_free, g_driven) = (coupling(free), coupling(driven));
        let capacitance = EPSILON[2] * LENGTH;
        let mass = &model.solid().mass;
        let stiffness = &model.solid().stiffness;
        assert!(entry(mass, free, driven) > 0.0);
        assert_ne!(g_free, 0.0);
        assert_ne!(g_driven, 0.0);
        let body = Body::new(
            &model,
            &fixed,
            &electrodes(&cube, 0.0, Condition::Charge(0.0)),
            &zero,
            &zero,
            &charge,
            controls(),
            cx,
        )
        .unwrap();
        let dt = 1e-7_f64;
        let velocity = 0.03;
        let supports = fixed
            .iter()
            .map(|&(dof, _)| (dof, if dof == driven { velocity } else { 0.0 }))
            .collect::<Vec<_>>();
        // Independent charge condensation and scalar mechanical elimination.
        // The mixed consistent-mass term is essential to base excitation.
        let factor = 4.0 / (dt * dt);
        let half_base = 0.25 * dt * velocity;
        let half_free = -half_base
            * (entry(stiffness, free, driven)
                + g_free * g_driven / capacitance
                + factor * entry(mass, free, driven))
            / (entry(stiffness, free, free)
                + g_free * g_free / capacitance
                + factor * entry(mass, free, free));
        let midpoint_voltage = (g_free * half_free + g_driven * half_base) / capacitance;
        let reaction = (entry(stiffness, driven, free)
            + factor * entry(mass, driven, free))
            * half_free
            + (entry(stiffness, driven, driven) + factor * entry(mass, driven, driven))
                * half_base
            + g_driven * midpoint_voltage;
        let candidate = body
            .step_with_support_motion_candidate(
                dt,
                &supports,
                &[Condition::Voltage(0.0), Condition::Charge(0.0)],
                &zero,
                &charge,
                controls(),
                cx,
            )
            .unwrap();
        let state = candidate.body().state();
        near(
            state.displacement_m()[free / 3][2],
            2.0 * half_free,
            1e-20,
            1e-9,
        );
        near(
            state.electrodes()[1].voltage_v,
            2.0 * midpoint_voltage,
            1e-12,
            1e-9,
        );
        near(state.electrodes()[1].charge_c, 0.0, 1e-20, 0.0);
        near(
            candidate.midpoint_support_force_n()[driven / 3][2],
            reaction,
            1e-10,
            1e-9,
        );
        assert_eq!(candidate.midpoint_support_force_n()[free / 3][2], 0.0);
        near(
            candidate.step().support_work_j,
            reaction * 2.0 * half_base,
            1e-21,
            1e-9,
        );
        near(
            state.total_energy_j(),
            candidate.step().support_work_j,
            1e-21,
            1e-9,
        );
        assert!(state.electrodes()[1].voltage_v.abs() > 1e-4);
    });
}

#[test]
fn g1_imposed_strain_supplies_harvesting_heat_and_external_capacitor_storage() {
    with_cx(|cx, _| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let materials = [material([0.0, 0.0, 12.0, 0.0, 0.0, 0.0])];
        let model =
            Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let (zero, charge) = zero_source(n);
        let fixed = (0..3 * n).map(|i| (i, 0.0)).collect::<Vec<_>>();
        let body = Body::new(
            &model,
            &fixed,
            &electrodes(&cube, 0.0, Condition::Charge(0.0)),
            &zero,
            &zero,
            &charge,
            controls(),
            cx,
        )
        .unwrap();
        let dt = 1e-7_f64;
        let rate = 20.0;
        let supports = (0..3 * n)
            .map(|dof| {
                (
                    dof,
                    if dof % 3 == 2 {
                        rate * cube.nodes[dof / 3][2]
                    } else {
                        0.0
                    },
                )
            })
            .collect::<Vec<_>>();
        let resistance = 1e4;
        let load = super::super::transient_load::resistor(resistance, [1, 0]);
        let capacitance = EPSILON[2] * LENGTH;
        let midpoint_strain = 0.25 * dt * rate;
        let piezo_charge = 12.0 * LENGTH.powi(2) * midpoint_strain;
        let midpoint_voltage = piezo_charge / (capacitance + 0.5 * dt / resistance);
        let expected_heat = dt * midpoint_voltage.powi(2) / resistance;
        let drives = [ElectrodeDrive::Voltage(0.0), ElectrodeDrive::Current(0.0)];
        let (candidate, report) = body
            .step_with_loaded_support_motion_candidate(
                dt,
                &supports,
                &drives,
                std::slice::from_ref(&load),
                &zero,
                &charge,
                controls(),
                cx,
            )
            .unwrap();
        near(
            candidate.body().state().electrodes()[1].voltage_v,
            2.0 * midpoint_voltage,
            1e-12,
            1e-10,
        );
        near(
            report.loads[0].midpoint_current_a,
            midpoint_voltage / resistance,
            1e-17,
            1e-10,
        );
        near(report.joule_heat_j, expected_heat, 1e-22, 1e-10);
        assert!(report.joule_heat_j > 0.0);
        near(report.electrical_source_work_j, 0.0, 1e-22, 0.0);
        assert_eq!(candidate.step().applied_force_work_j, 0.0);
        near(
            candidate.step().support_work_j,
            candidate.body().state().total_energy_j() + expected_heat,
            1e-21,
            1e-10,
        );
        near(report.energy_residual_j, 0.0, 1e-21, 0.0);

        // The same imposed strain now drives a real source-resolved capacitor
        // in parallel with the conductor. Its storage must not become heat.
        let external = super::super::transient_load::reactive::capacitor(2.0 * capacitance);
        let coupled = ReactivePiezoelectricTransient::new(
            body.clone(),
            vec![ReactiveElectrodeLoad::Capacitor {
                electrodes: [1, 0],
                capacitor: external.clone(),
            }],
            vec![ReactiveElectrodeState::Charge(0.0)],
            vec![load],
            controls(),
            cx,
        )
        .unwrap();
        let candidate = coupled
            .step_with_support_motion_candidate(
                dt,
                &supports,
                &drives,
                &zero,
                &charge,
                controls(),
                cx,
            )
            .unwrap();
        let midpoint_voltage =
            piezo_charge / (capacitance + external.capacitance_f() + 0.5 * dt / resistance);
        let heat = dt * midpoint_voltage.powi(2) / resistance;
        let storage = 0.5 * external.capacitance_f() * (2.0 * midpoint_voltage).powi(2);
        let report = candidate.step();
        near(
            candidate.model().body().state().electrodes()[1].voltage_v,
            2.0 * midpoint_voltage,
            1e-12,
            1e-10,
        );
        near(report.joule_heat_j, heat, 1e-22, 1e-10);
        near(report.reactive_energy_j[1], storage, 1e-22, 1e-10);
        near(
            report.body.support_work_j,
            candidate.model().body().state().total_energy_j() + storage + heat,
            1e-21,
            1e-10,
        );
        near(report.electrical_source_work_j, 0.0, 1e-22, 0.0);
        near(report.energy_residual_j, 0.0, 1e-21, 0.0);
        assert_eq!(report.midpoint_support_force_n.len(), n);
        assert_eq!(body.state().time_s(), 0.0);
    });
}

#[test]
fn g3_small_kinetic_exchange_survives_rigid_motion_through_native_and_loaded_owners() {
    with_cx(|cx, _| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let materials = [material([0.0, 0.0, 12.0, 0.0, 0.0, 0.0])];
        let model = Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let (zero, charge) = zero_source(n);
        let fixed = (0..3 * n).map(|dof| (dof, 0.0)).collect::<Vec<_>>();
        let original = Body::new(&model, &fixed,
            &electrodes(&cube, 0.0, Condition::Charge(0.0)),
            &zero, &zero, &charge, controls(), cx).unwrap();
        let dt = 2.0_f64.powi(-20);
        let rate = 0.25;
        let expected_kinetic = 7500.0 * LENGTH.powi(5) * rate * rate / 6.0;
        let tight = Controls { energy_tolerance_j: 1e-18, ..controls() };
        let conditions = [Condition::Voltage(0.0), Condition::Charge(0.0)];
        let drives = [ElectrodeDrive::Voltage(0.0), ElectrodeDrive::Current(0.0)];
        let resistor = super::super::transient_load::resistor(1e4, [1, 0]);
        let capacitor = super::super::transient_load::reactive::capacitor(2.0 * EPSILON[2] * LENGTH);
        let mut reference = None;
        for boost in [0.0, 2.0_f64.powi(16)] {
            // A dyadic startup supplies the complete rigid kinetic inventory.
            // The subsequent transverse strain has its own tiny inertia/work;
            // the common translation must neither supply heat nor enlarge a gate.
            let startup = (0..3 * n)
                .map(|dof| (dof, if dof % 3 == 0 { boost } else { 0.0 }))
                .collect::<Vec<_>>();
            let body = original.step_with_support_motion_candidate(
                dt, &startup, &conditions, &zero, &charge, controls(), cx,
            ).unwrap().into_parts().0;
            let supports = (0..3 * n).map(|dof| (dof, match dof % 3 {
                0 => boost,
                2 => rate * cube.nodes[dof / 3][2],
                _ => 0.0,
            })).collect::<Vec<_>>();
            let native = body.step_with_support_motion_candidate(
                dt, &supports, &conditions, &zero, &charge, tight, cx,
            ).unwrap();
            let (loaded, shunt) = body.step_with_loaded_support_motion_candidate(
                dt, &supports, &drives, std::slice::from_ref(&resistor),
                &zero, &charge, tight, cx,
            ).unwrap();
            let reactive = ReactivePiezoelectricTransient::new(
                body.clone(),
                vec![ReactiveElectrodeLoad::Capacitor {
                    electrodes: [1, 0], capacitor: capacitor.clone(),
                }],
                vec![ReactiveElectrodeState::Charge(0.0)],
                vec![resistor.clone()], tight, cx,
            ).unwrap().step_with_support_motion_candidate(
                dt, &supports, &drives, &zero, &charge, tight, cx,
            ).unwrap();
            for step in [native.step(), loaded.step(), &reactive.step().body] {
                near(step.kinetic_energy_change_j, expected_kinetic, 1e-24, 1e-10);
                assert!(step.energy_residual_j.abs() < tight.energy_tolerance_j);
                assert_eq!(step.material_heat_j, 0.0);
            }
            assert!(shunt.joule_heat_j > 0.0);
            assert!(reactive.step().reactive_energy_j[1] > 0.0);
            assert!(shunt.energy_residual_j.abs() < tight.energy_tolerance_j);
            assert!(reactive.step().energy_residual_j.abs() < tight.energy_tolerance_j);
            let measured = [native.step().support_work_j, loaded.step().support_work_j,
                shunt.joule_heat_j, loaded.body().state().electrodes()[1].voltage_v,
                reactive.step().body.support_work_j, reactive.step().joule_heat_j,
                reactive.step().reactive[0].stored_energy_change_j];
            if let Some(reference) = reference {
                for (actual, expected) in measured.into_iter().zip(reference) {
                    near(actual, expected, 1e-24, 1e-9);
                }
                assert!(body.state().kinetic_energy_j() > 1e17 * expected_kinetic);
            } else {
                reference = Some(measured);
            }
            assert_eq!(body.state().time_s(), dt);
            assert!(loaded.body().state().velocity_m_s().iter().all(|v| v[0] == boost));
        }
    });
}

#[test]
fn g4_support_topology_and_cancelled_motion_refuse_without_mutation() {
    with_cx(|cx, gate| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let materials = [material([0.0, 0.0, 12.0, 0.0, 0.0, 0.0])];
        let model =
            Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let (zero, charge) = zero_source(n);
        let body = Body::new(
            &model,
            &cube.supports,
            &electrodes(&cube, 0.0, Condition::Charge(0.0)),
            &zero,
            &zero,
            &charge,
            controls(),
            cx,
        )
        .unwrap();
        let attempt = |supports: &[(usize, f64)]| {
            body.step_with_support_motion_candidate(
                1e-7,
                supports,
                &[Condition::Voltage(0.0), Condition::Charge(0.0)],
                &zero,
                &charge,
                controls(),
                cx,
            )
        };
        assert!(attempt(&cube.supports[..cube.supports.len() - 1]).is_err());
        let mut invalid = cube.supports.clone();
        invalid[1] = invalid[0];
        assert!(attempt(&invalid).is_err());
        invalid[1] = (3 * *cube.top.last().unwrap() + 2, 0.0);
        assert!(attempt(&invalid).is_err());
        invalid[1] = (3 * n, 0.0);
        assert!(attempt(&invalid).is_err());
        invalid[1] = (cube.supports[1].0, f64::NAN);
        assert!(attempt(&invalid).is_err());
        gate.request();
        assert!(matches!(
            attempt(&cube.supports),
            Err(TetPiezoelectricError::Cancelled)
        ));
        assert_eq!(body.state().displacement_m(), &zero);
        assert_eq!(body.state().velocity_m_s(), &zero);
        assert_eq!(body.state().time_s(), 0.0);
    });
}
