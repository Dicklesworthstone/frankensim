//! Original dense-gas, finite thermal, electrical and acoustic storage together.

use super::*;
use crate::cavity_phs::piston::double_acting::real_gas::{real_controls, real_pair};
use crate::cavity_phs::piston::material::tests::body;
use crate::dg_port::sum_terminal_work;
use fs_ports::PortedModule;

const DENSE_DT: f64 = 2e-6;
// Small additional acoustic faces keep the fixture's 1 mm ducts in their
// perturbation regime even for the light specimen's faster second interval.
const DENSE_AREAS: [f64; 2] = [5e-5, -3e-5];

fn dense_drive() -> ElectromagneticNetworkDrive {
    let mut drive = super::drive();
    drive.gas.input.duration_s = DENSE_DT;
    drive.gas.piston = real_controls();
    for face in &mut drive.gas.exchange.faces {
        face.heat_bracket_j = [-50.0 * DENSE_DT, 50.0 * DENSE_DT];
    }
    drive.gas.exchange.maximum_joint_trials = 128;
    drive.initial_midpoint_velocity_m_s = 0.4;
    drive
}

#[test]
fn dense_network_and_native_acoustic_replay_close_every_independent_store() {
    let (a, ca) = fixture::card(0.2, 2.0, 0.002);
    let (b, cb) = fixture::card(0.4, 8.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let context = fixture::point(400.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = fixture::cells(&cards, &curves, &context, &key);
    let resistors = fixture::resistor();
    let binding = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &[],
    };
    let thermal = fixture::initial(&curves);
    let circuit = fixture::circuit(&binding, &thermal);
    let piston = real_pair();
    let loads = vec![
        acoustic_fixture::load(0.15, DENSE_DT),
        acoustic_fixture::load(0.85, DENSE_DT),
    ];
    let drive = dense_drive();
    let candidate = piston
        .electroacoustic_network_candidate(
            &faces(&loads, DENSE_AREAS),
            &binding,
            &thermal,
            DENSE_DT,
            &circuit,
            &[0.0],
            drive,
            fixture::controls(),
            acoustic(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    let source = &candidate.native().network().report().source;
    let step = &source.boundary.step().piston;
    let displacement = step.displacement_m;
    let velocity = displacement / DENSE_DT;
    acoustic_fixture::near(
        source.motion.back_emf_v,
        drive.transducer.factor_n_a * velocity,
        1e-7,
    );
    acoustic_fixture::near(
        source.motion.force_n,
        drive.transducer.factor_n_a * source.motion.current_a,
        1e-12,
    );
    let observed_force: f64 = candidate
        .motion()
        .faces
        .iter()
        .map(|face| -face.signed_area_m2 * face.mean_pressure_pa)
        .sum();
    acoustic_fixture::near(
        piston.primary().mass_kg() * (step.velocity_m_s[1] - step.velocity_m_s[0]),
        DENSE_DT
            * (piston.primary().area_m2() * step.midpoint[0].pressure_pa()
                - piston.opposing_area_m2() * step.midpoint[1].pressure_pa()
                + drive.gas.input.external_force_n
                + source.motion.force_n
                + observed_force),
        real_controls().piston.momentum_tolerance_n_s,
    );
    for i in 0..2 {
        let mut replay = loads[i].clone();
        let output = replay
            .preview_step(DENSE_AREAS[i] * velocity)
            .unwrap()
            .accept();
        assert_eq!(candidate.motion().faces[i].duct, output);
        acoustic_fixture::assert_load_same(&candidate.loads()[i], &replay);
        assert_eq!(candidate.loads()[i].state().step_count(), 1);
        assert_eq!(
            candidate.loads()[i].duct().gas_field(),
            loads[i].duct().gas_field()
        );
    }
    let unloaded = piston
        .electromagnetic_network_candidate(
            &binding,
            &thermal,
            DENSE_DT,
            &circuit,
            &[0.0],
            drive,
            fixture::controls(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    assert!(
        candidate.native().piston().primary().velocity_m_s()
            < unloaded.piston().primary().velocity_m_s()
    );
    assert!(source.motion.current_a > unloaded.network().report().source.motion.current_a);
    assert!(
        candidate.native().network().report().electrical.heat_j
            > unloaded.network().report().electrical.heat_j
    );
    let (native, next_loads, motion, total) = candidate.into_parts();
    let (next_thermal, next_circuit, next_piston, report) = native.into_parts();
    // Reconstruct all finite changes independently. Original PR increments
    // avoid subtracting absolute formation energies; neither contact heat nor
    // resistor heat appears as an external supply in this enclosing balance.
    let (q0, q1) = (circuit.charges_c()[0], next_circuit.charges_c()[0]);
    let (v0, v1) = (
        piston.primary().velocity_m_s(),
        next_piston.primary().velocity_m_s(),
    );
    let mut energy = (q1 - q0) * (q1 + q0) / (2.0 * 0.01)
        + 0.5 * piston.primary().mass_kg() * (v1 - v0) * (v1 + v0);
    for (old, new) in [
        (piston.primary().gas(), next_piston.primary().gas()),
        (piston.opposing(), next_piston.opposing()),
    ] {
        assert_eq!(new.mass_kg(), old.mass_kg());
        assert_eq!(new.model(), old.model());
        assert_eq!(new.point().axes().get("specimen-lot"), Some(&42.0));
        assert_eq!(new.density_kg_m3(), new.mass_kg() / new.volume_m3());
        energy += old.mass_kg()
            * old
                .model()
                .specific_internal_energy_change_between_densities(
                    old.density_kg_m3(),
                    old.state().temperature_k(),
                    new.density_kg_m3(),
                    new.state().temperature_k(),
                )
                .unwrap();
    }
    for (id, old) in thermal.cells() {
        energy += old.mass_kg * (next_thermal.cells()[id].specific_h_j_kg - old.specific_h_j_kg);
    }
    for (old, new) in loads.iter().zip(&next_loads) {
        energy += new.stored_energy_j() - old.stored_energy_j() + new.cumulative_dissipated_j()
            - old.cumulative_dissipated_j()
            + sum_terminal_work(new.last_terminal_frames());
    }
    let supplied = drive.gas.input.heat_j.iter().sum::<f64>()
        + drive.gas.input.external_force_n * displacement;
    acoustic_fixture::near(
        energy,
        supplied,
        fixture::controls().total_energy_tolerance_j,
    );
    acoustic_fixture::near(total, 0.0, acoustic().energy_tolerance_j);
    acoustic_fixture::near(
        motion.transfer_residual_j,
        0.0,
        acoustic().energy_tolerance_j,
    );
    acoustic_fixture::near(
        report.motion.transfer_residual_j,
        0.0,
        drive.motion.energy_tolerance_j,
    );
    acoustic_fixture::near(
        report.network.deposited_joule_j,
        report.network.electrical.heat_j,
        1e-12,
    );
    acoustic_fixture::near(report.motion.external_electrical_work_j, 0.0, 1e-15);
    acoustic_fixture::near(report.coupled.external_network_heat_j, 0.0, 1e-15);
    next_piston.audit_last_window().unwrap();
}

#[test]
fn source_bound_dense_actuator_keeps_all_clocks_and_final_cancellation_is_atomic() {
    let (a, ca) = fixture::card(0.2, 2.0, 0.002);
    let (b, cb) = fixture::card(0.4, 8.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let context = fixture::point(400.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = fixture::cells(&cards, &curves, &context, &key);
    let resistors = fixture::resistor();
    let binding = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &[],
    };
    let mut thermal = fixture::initial(&curves);
    let mut circuit = fixture::circuit(&binding, &thermal);
    let gas = real_pair();
    let specimen = body(1000.0);
    let mut piston = MaterialDoubleActingGasPiston::new(
        gas.primary().gas().clone(),
        gas.opposing().clone(),
        specimen.clone(),
        [0.01, 0.008],
        0.4,
    )
    .unwrap();
    assert_ne!(piston.piston().primary().mass_kg(), gas.primary().mass_kg());
    let mut loads = vec![
        acoustic_fixture::load(0.15, DENSE_DT),
        acoustic_fixture::load(0.85, DENSE_DT),
    ];
    for index in 0..2u64 {
        // These first two decimal intervals retain exactly the propagator's dt.
        let target = thermal.time_s() + DENSE_DT;
        assert_eq!(target - thermal.time_s(), DENSE_DT);
        let old_piston = piston.clone();
        let old_circuit = circuit.clone();
        let old_cells = thermal.cells().clone();
        let old_loads = loads.clone();
        let mut polls = 0;
        let candidate = piston
            .electroacoustic_network_candidate(
                &faces(&loads, DENSE_AREAS),
                &binding,
                &thermal,
                target,
                &circuit,
                &[0.0],
                dense_drive(),
                fixture::controls(),
                acoustic(),
                &[0.0; 2],
                || {
                    polls += 1;
                    false
                },
            )
            .unwrap();
        assert_eq!(candidate.native().piston().body(), &specimen);
        assert_eq!(
            candidate.native().piston().piston().primary().mass_kg(),
            specimen.inertia().mass()
        );
        assert_eq!(
            candidate.native().piston().piston().ports.next_tick(),
            index + 2
        );
        if index == 0 {
            let native = piston
                .piston()
                .electroacoustic_network_candidate(
                    &faces(&loads, DENSE_AREAS),
                    &binding,
                    &thermal,
                    target,
                    &circuit,
                    &[0.0],
                    dense_drive(),
                    fixture::controls(),
                    acoustic(),
                    &[0.0; 2],
                    || false,
                )
                .unwrap();
            assert_eq!(
                candidate.native().piston().piston(),
                native.native().piston()
            );
            assert_eq!(
                candidate.native().circuit(),
                native.native().network().circuit()
            );
            assert_eq!(
                candidate.native().thermal().cells(),
                native.native().network().thermal().cells()
            );
            assert_eq!(candidate.motion(), native.motion());
            for (actual, expected) in candidate.loads().iter().zip(native.loads()) {
                acoustic_fixture::assert_load_same(actual, expected);
            }
        } else {
            let mut count = 0;
            assert!(matches!(
                piston.electroacoustic_network_candidate(
                    &faces(&loads, DENSE_AREAS),
                    &binding,
                    &thermal,
                    target,
                    &circuit,
                    &[0.0],
                    dense_drive(),
                    fixture::controls(),
                    acoustic(),
                    &[0.0; 2],
                    || {
                        count += 1;
                        count == polls
                    },
                ),
                Err(AcousticMotionError::Cancelled)
            ));
            assert_eq!(count, polls);
        }
        assert_eq!(piston, old_piston);
        assert_eq!(circuit, old_circuit);
        assert_eq!(thermal.cells(), &old_cells);
        for (actual, expected) in loads.iter().zip(&old_loads) {
            acoustic_fixture::assert_load_same(actual, expected);
        }
        let displacement = candidate.native().report().boundary.piston.displacement_m;
        for i in 0..2 {
            let mut replay = loads[i].clone();
            replay
                .preview_step(DENSE_AREAS[i] * (displacement / DENSE_DT))
                .unwrap()
                .accept();
            acoustic_fixture::assert_load_same(&candidate.loads()[i], &replay);
            assert_eq!(candidate.loads()[i].state().step_count(), index + 1);
        }
        let (native, next_loads, _, _) = candidate.into_parts();
        let (next_thermal, next_circuit, next_piston, _) = native.into_parts();
        thermal = next_thermal;
        circuit = next_circuit;
        piston = next_piston;
        loads = next_loads;
        piston.piston().audit_last_window().unwrap();
    }
    assert_eq!(piston.body(), &specimen);
    assert_eq!(piston.piston().ports.next_tick(), 3);
    assert_eq!(thermal.time_s(), 2.0 * DENSE_DT);
    assert_eq!(piston.gases()[0].model(), gas.primary().gas().model());
    assert_eq!(piston.gases()[1].model(), gas.opposing().model());
}
