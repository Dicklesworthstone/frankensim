//! Dense independently composed chambers retain their native electrical/thermal owners.

use super::*;
use crate::cavity_phs::piston::double_acting::real_gas::{
    real_controls as mechanical_controls, real_pair,
};
use crate::cavity_phs::piston::material::tests::body;
use fs_material::thermochemical_gas::nonideal::energy::parcel::IsochoricRealGasParcel;
use fs_ports::PortedModule;

fn dense_drive(duration_s: f64) -> ElectromagneticNetworkDrive {
    let mut drive = super::drive(0.01);
    drive.gas.input.duration_s = duration_s;
    drive.gas.piston = mechanical_controls();
    // The finite contacts exchange about 1e-4 J at 20 us. Bracket that
    // physical scale so neither nested scalar search wastes its joint budget.
    let heat_bound = 50.0 * duration_s;
    for face in &mut drive.gas.exchange.faces {
        face.heat_bracket_j = [-heat_bound, heat_bound];
    }
    drive.gas.exchange.maximum_joint_trials = 128;
    drive.initial_midpoint_velocity_m_s = 0.4;
    drive
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.16e} != {expected:.16e} (allowance {tolerance:.3e})"
    );
}

fn gas_energy(old: &IsochoricRealGasParcel, new: &IsochoricRealGasParcel) -> f64 {
    old.mass_kg()
        * old
            .model()
            .specific_internal_energy_change_between_densities(
                old.density_kg_m3(),
                old.state().temperature_k(),
                new.density_kg_m3(),
                new.state().temperature_k(),
            )
            .unwrap()
}

fn whole_energy_change(
    before: &DoubleActingGasPiston<IsochoricRealGasParcel>,
    after: &DoubleActingGasPiston<IsochoricRealGasParcel>,
    thermal_before: &DepositionSimulation<'_>,
    thermal_after: &DepositionSimulation<'_>,
    circuit_before: &MaterialCircuit,
    circuit_after: &MaterialCircuit,
) -> f64 {
    let q0 = circuit_before.charges_c()[0];
    let q1 = circuit_after.charges_c()[0];
    let v0 = before.primary().velocity_m_s();
    let v1 = after.primary().velocity_m_s();
    // Independent original-law increments, not the adapter's total residual.
    // Joule loss and both contact heats are internal to these finite stores.
    (q1 - q0) * (q1 + q0) / (2.0 * 0.01)
        + thermal_before
            .cells()
            .iter()
            .map(|(id, old)| {
                old.mass_kg * (thermal_after.cells()[id].specific_h_j_kg - old.specific_h_j_kg)
            })
            .sum::<f64>()
        + gas_energy(before.primary().gas(), after.primary().gas())
        + gas_energy(before.opposing(), after.opposing())
        + 0.5 * before.primary().mass_kg() * (v1 - v0) * (v1 + v0)
}

#[test]
fn dense_motion_closes_native_back_emf_joule_heat_and_every_finite_store() {
    let (a, ca) = card(0.2, 2.0, 0.002);
    let (b, cb) = card(0.4, 8.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let context = point(400.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(&cards, &curves, &context, &key);
    let resistors = resistor();
    let binding = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &[],
    };
    let thermal = initial(&curves);
    let circuit = circuit(&binding, &thermal);
    let piston = real_pair();
    let drive = dense_drive(20e-6);
    let candidate = piston
        .electromagnetic_network_candidate(
            &binding,
            &thermal,
            drive.gas.input.duration_s,
            &circuit,
            &[0.0],
            drive,
            controls(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    let (next_thermal, next_circuit, next_piston, report) = candidate.into_parts();
    close(
        whole_energy_change(
            &piston,
            &next_piston,
            &thermal,
            &next_thermal,
            &circuit,
            &next_circuit,
        ),
        0.0,
        controls().total_energy_tolerance_j,
    );
    let motion = &report.motion;
    let paired = &report.boundary.piston;
    close(
        motion.back_emf_v,
        drive.transducer.factor_n_a * f64::midpoint(paired.velocity_m_s[0], paired.velocity_m_s[1]),
        2.0 * drive.motion.velocity_tolerance_m_s,
    );
    close(
        motion.force_n,
        drive.transducer.factor_n_a * motion.current_a,
        1e-12,
    );
    close(
        motion.mechanical_work_j,
        motion.force_n * paired.displacement_m,
        1e-12,
    );
    close(
        motion.electrical_work_j + motion.mechanical_work_j,
        0.0,
        drive.motion.energy_tolerance_j,
    );
    close(motion.external_electrical_work_j, 0.0, 1e-15);
    let resistor = &report.network.electrical.resistors[0];
    let endpoint_r = 2.0 * (1.0 + 0.002 * (next_thermal.temperature_k(CellId(0)).unwrap() - 400.0));
    close(
        resistor.material.conductor.resistance_ohm(),
        endpoint_r,
        1e-12,
    );
    close(resistor.current_a, motion.current_a, 1e-10);
    close(
        resistor.heat_j,
        drive.gas.input.duration_s * endpoint_r * motion.current_a.powi(2),
        1e-12,
    );
    close(report.network.deposited_joule_j, resistor.heat_j, 1e-12);
    assert!(resistor.heat_j > 0.0 && motion.mechanical_work_j > 0.0);
    assert_ne!(resistor.material.conductor.resistance_ohm(), 2.0);
    let original = [piston.primary().gas(), piston.opposing()];
    let next = [next_piston.primary().gas(), next_piston.opposing()];
    for i in 0..2 {
        let midpoint = &paired.midpoint[i];
        let transport = midpoint.transport_state().unwrap();
        let geometry = drive.geometry[i];
        let conductance = geometry.area_m2
            / (geometry.gas_distance_m / transport.thermal_conductivity_w_m_k()
                + geometry.cell_distance_m / [0.2, 0.4][i]);
        let temperature = next_thermal.temperature_k(CellId(i as u64)).unwrap();
        close(report.boundary.conductance_w_k[i], conductance, 1e-14);
        close(
            report.boundary.exchanged_heat_j[i],
            drive.gas.input.duration_s * conductance * (temperature - midpoint.temperature_k()),
            drive.gas.exchange.faces[i].source_tolerance_j,
        );
        assert!(report.boundary.exchanged_heat_j[i] > 0.0);
        assert_eq!(transport.density_kg_m3(), midpoint.density_kg_m3());
        assert_ne!(transport.density_kg_m3(), midpoint.source().state().density);
        assert_eq!(next[i].model(), original[i].model());
        assert_eq!(next[i].mass_kg(), original[i].mass_kg());
        assert_eq!(next[i].point().axes().get("specimen-lot"), Some(&42.0));
        let sample = report
            .network
            .source
            .cells
            .iter()
            .find(|source| source.cell == CellId(i as u64))
            .unwrap();
        assert_eq!(sample.conductivity.as_ref(), Some(&report.interfaces[i]));
        let deposited = if i == 0 { resistor.heat_j } else { 0.0 };
        let old_cell = &thermal.cells()[&CellId(i as u64)];
        let new_cell = &next_thermal.cells()[&CellId(i as u64)];
        close(
            old_cell.mass_kg * (new_cell.specific_h_j_kg - old_cell.specific_h_j_kg),
            deposited - report.boundary.exchanged_heat_j[i],
            controls().network.cell_energy_tolerance_j,
        );
    }
    close(
        piston.primary().mass_kg() * (paired.velocity_m_s[1] - paired.velocity_m_s[0]),
        drive.gas.input.duration_s
            * (piston.primary().area_m2() * paired.midpoint[0].pressure_pa()
                - piston.opposing_area_m2() * paired.midpoint[1].pressure_pa()
                + motion.force_n),
        mechanical_controls().piston.momentum_tolerance_n_s,
    );
    next_piston.audit_last_window().unwrap();
}

#[test]
fn sourced_solid_density_changes_current_and_both_real_gas_histories() {
    let (a, ca) = card(0.2, 2.0, 0.002);
    let (b, cb) = card(0.4, 8.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let context = point(400.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(&cards, &curves, &context, &key);
    let resistors = resistor();
    let binding = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &[],
    };
    let initial_thermal = initial(&curves);
    let initial_circuit = circuit(&binding, &initial_thermal);
    let gas = real_pair();
    let run = |density| {
        let specimen = body(density);
        let mut piston = MaterialDoubleActingGasPiston::new(
            gas.primary().gas().clone(),
            gas.opposing().clone(),
            specimen.clone(),
            [gas.primary().area_m2(), gas.opposing_area_m2()],
            gas.primary().velocity_m_s(),
        )
        .unwrap();
        let original = piston.clone();
        let mut thermal = initial_thermal.clone();
        let mut circuit = initial_circuit.clone();
        let mut last_electrical = [0.0; 2];
        for index in 0..2 {
            let target = thermal.time_s() + 2e-6;
            let drive = dense_drive(target - thermal.time_s());
            let candidate = piston
                .electromagnetic_network_candidate(
                    &binding,
                    &thermal,
                    target,
                    &circuit,
                    &[0.0],
                    drive,
                    controls(),
                    &[0.0; 2],
                    || false,
                )
                .unwrap();
            if index == 0 && density == 1000.0 {
                let native = piston
                    .piston()
                    .electromagnetic_network_candidate(
                        &binding,
                        &thermal,
                        target,
                        &circuit,
                        &[0.0],
                        drive,
                        controls(),
                        &[0.0; 2],
                        || false,
                    )
                    .unwrap();
                assert_eq!(candidate.piston().piston(), native.piston());
                assert_eq!(candidate.circuit(), native.network().circuit());
                assert_eq!(
                    candidate.report().boundary,
                    *native.network().report().source.boundary.step()
                );
                assert_eq!(
                    candidate.piston().piston().last_window_transfers(),
                    native.piston().last_window_transfers()
                );
            }
            let (next_thermal, next_circuit, next_piston, report) = candidate.into_parts();
            thermal = next_thermal;
            circuit = next_circuit;
            piston = next_piston;
            last_electrical = [report.motion.back_emf_v, report.motion.current_a];
            assert_eq!(piston.body(), &specimen);
            assert_eq!(
                piston.piston().primary().mass_kg(),
                specimen.inertia().mass()
            );
            for (new, old) in piston.gases().into_iter().zip(original.gases()) {
                assert_eq!(new.mass_kg(), old.mass_kg());
                assert_eq!(new.model(), old.model());
                assert_eq!(new.density_kg_m3(), new.mass_kg() / new.volume_m3());
            }
            close(
                whole_energy_change(
                    original.piston(),
                    piston.piston(),
                    &initial_thermal,
                    &thermal,
                    &initial_circuit,
                    &circuit,
                ),
                0.0,
                (index + 1) as f64 * controls().total_energy_tolerance_j,
            );
            piston.piston().audit_last_window().unwrap();
        }
        (piston, last_electrical)
    };
    let (light, light_electrical) = run(1000.0);
    let (dense, dense_electrical) = run(8000.0);
    assert!(light.piston().primary().velocity_m_s() > dense.piston().primary().velocity_m_s());
    assert!(light_electrical[0] > dense_electrical[0]);
    assert!(light_electrical[1] < dense_electrical[1]);
    assert!(light.gases()[0].density_kg_m3() < dense.gases()[0].density_kg_m3());
    assert!(light.gases()[1].density_kg_m3() > dense.gases()[1].density_kg_m3());
}

#[test]
fn missing_contact_source_and_final_cancel_preserve_all_accepted_owners() {
    let (a, ca) = card(0.2, 2.0, 0.002);
    let (b, cb) = card(0.4, 8.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let context = point(400.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(&cards, &curves, &context, &key);
    let resistors = resistor();
    let binding = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &[],
    };
    let thermal = initial(&curves);
    let circuit = circuit(&binding, &thermal);
    let gas = real_pair();
    let piston = MaterialDoubleActingGasPiston::new(
        gas.primary().gas().clone(),
        gas.opposing().clone(),
        body(8000.0),
        [0.01, 0.008],
        0.4,
    )
    .unwrap();
    let (thermal, circuit, piston, _) = piston
        .electromagnetic_network_candidate(
            &binding,
            &thermal,
            2e-6,
            &circuit,
            &[0.0],
            dense_drive(2e-6),
            controls(),
            &[0.0; 2],
            || false,
        )
        .unwrap()
        .into_parts();
    let saved_piston = piston.clone();
    let saved_circuit = circuit.clone();
    let saved_cells = thermal.cells().clone();
    let target = thermal.time_s() + 2e-6;
    let drive = dense_drive(target - thermal.time_s());
    let mut missing = super::cells(&cards, &curves, &context, &key);
    missing[1].conductivity = None;
    let missing_binding = MaterialNetwork {
        cells: &missing,
        resistors: &resistors,
        links: &[],
    };
    assert!(matches!(
        piston.electromagnetic_network_candidate(
            &missing_binding,
            &thermal,
            target,
            &circuit,
            &[0.0],
            drive,
            controls(),
            &[0.0; 2],
            || false,
        ),
        Err(CircuitEnthalpyError::Source(DrivenSourceError::Motion(
            ElectricalMotionError::Mechanical(ThermalPistonError::Interface(_))
        )))
    ));
    let mut polls = 0;
    piston
        .electromagnetic_network_candidate(
            &binding,
            &thermal,
            target,
            &circuit,
            &[0.0],
            drive,
            controls(),
            &[0.0; 2],
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut count = 0;
    assert!(matches!(
        piston.electromagnetic_network_candidate(
            &binding,
            &thermal,
            target,
            &circuit,
            &[0.0],
            drive,
            controls(),
            &[0.0; 2],
            || {
                count += 1;
                count == polls
            },
        ),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    assert_eq!(count, polls);
    assert_eq!(piston, saved_piston);
    assert_eq!(circuit, saved_circuit);
    assert_eq!(thermal.cells(), &saved_cells);
    assert_eq!(thermal.time_s(), 2e-6);
    assert_eq!(
        piston.piston().last_window_transfers(),
        saved_piston.piston().last_window_transfers()
    );
    piston.piston().audit_last_window().unwrap();
}
