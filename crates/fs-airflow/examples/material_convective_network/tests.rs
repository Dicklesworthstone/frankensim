use super::*;
use fs_airflow::graph::thermal::medium::convection::evaluate_forced_from_material;
use fs_convection::ThermalDirection;
use fs_couple::cavity_phs::piston::thermal::ConductanceSample;
use fs_couple::cavity_phs::piston::thermal::boundary::ports::network as original;
use fs_couple::circuit_enthalpy::{CircuitEnthalpyError, CircuitEnthalpyInput};

#[test]
fn shared_adapter_matches_original_coupled_solver_and_preserves_each_source_point() {
    fixture(|binding, thermal, circuit| {
        let gas = gases(0.75); let drives = chambers(&gas, DT); let patches = contacts(1.0);
        let a = network::step_candidate(binding, &thermal, DT, &circuit, &[], &BASES,
            &drives, &patches, controls(), &[0.0; 2], || false).unwrap();
        let prepared = binding.prepare_with_bases(&thermal, &BASES, || false).unwrap();
        let destinations = binding.resistors.iter().map(|r| r.destination).collect::<Vec<_>>();
        let ports = patches.iter().map(|p| p.contact).collect::<Vec<_>>();
        let b = original::step_candidate(&thermal, DT,
            CircuitEnthalpyInput { circuit: &circuit, ports: &[], destinations: &destinations },
            &drives, &ports, controls(), &[0.0; 2], |t| prepared.resolve(t, || false),
            |contact, _, source, gas| {
                let patch = patches.iter().find(|p| p.contact == contact).unwrap();
                let cell = source.cells.iter().find(|s| s.cell == contact.cell).unwrap();
                let mut settings = patch.surface.settings;
                settings.direction = if cell.phase.temperature_k() > gas.state().temperature {
                    ThermalDirection::HeatingFluid
                } else { ThermalDirection::CoolingFluid };
                let resolved = evaluate_forced_from_material(patch.surface.correlation, gas,
                    patch.surface.speed_m_s, patch.surface.characteristic_length, settings)?;
                Ok::<_, fs_airflow::graph::thermal::medium::convection::MaterialConvectionError>(
                    ConductanceSample { conductance_w_k: patch.surface.area_m2*resolved.coefficient().value.value(), source: () })
            }, || false).unwrap();
        assert_eq!(a.thermal().cells(), b.thermal().cells());
        assert_eq!(a.circuit(), b.circuit());
        assert_eq!(a.material(), b.material());
        for (one, two) in a.chambers().iter().zip(b.chambers()) {
            assert_eq!(one.id, two.id);
            assert_eq!(one.boundary.piston(), two.boundary.piston());
            assert_eq!(one.boundary.step(), two.boundary.step());
            for surface in one.boundary.source() {
                let (contact, sample) = &surface.source;
                assert_eq!(sample.gas, one.boundary.step().piston.midpoint);
                assert_eq!(sample.body.temperature_k(), a.thermal().temperature_k(contact.cell).unwrap());
                assert_eq!(sample.gas.source_packs().len(), 2);
                assert_eq!(surface.conductance_w_k, sample.area_m2*sample.evaluation.as_ref().unwrap().coefficient.value.value());
            }
        }
        let mut reverse = patches; reverse.reverse();
        let reordered = network::step_candidate(binding, &thermal, DT, &circuit, &[], &BASES,
            &drives, &reverse, controls(), &[0.0; 2], || false).unwrap();
        assert_eq!(a.thermal().cells(), reordered.thermal().cells());
        assert_eq!(a.circuit(), reordered.circuit());
    });
}

#[test]
fn gas_material_data_reaches_other_cells_charge_motion_and_sound_over_time() {
    let a = run(0.75, 1.0, 24, false).unwrap();
    let b = run(0.25, 1.0, 24, false).unwrap();
    assert!((a.cell_k-b.cell_k).abs() > 1e-6);
    assert!((a.charge_c-b.charge_c).abs() > 1e-10);
    assert!((a.velocity_m_s-b.velocity_m_s).abs() > 1e-6);
    assert!((a.sound_m_s-b.sound_m_s).abs() > 100.0);
    assert!(a.energy_defect_j.abs() < 2e-5 && b.energy_defect_j.abs() < 2e-5);
}

#[test]
fn shared_cell_receives_both_signed_gas_transfers_and_closes_total_energy() {
    fixture(|binding, thermal, circuit| {
        let mut gas = gases(0.75);
        gas.insert(20, GasPiston::new(ClosedGasParcel::new(gas_source::model(0.5),
            gas_source::point(250.0, 100000.0), 0.001).unwrap(), 1.0, 0.01, 0.0).unwrap());
        let initial = storage(&thermal, &circuit, &gas);
        let drives = chambers(&gas, DT);
        let patches = [contacts(1.0)[0], contacts(1.0)[2]];
        let next = network::step_candidate(binding, &thermal, DT, &circuit, &[], &BASES,
            &drives, &patches, controls(), &[0.0; 2], || false).unwrap();
        let hot = &next.chambers()[0].boundary;
        let cold = &next.chambers()[1].boundary;
        assert!(hot.source()[0].heat_j < 0.0 && cold.source()[0].heat_j > 0.0);
        assert_eq!(hot.source()[0].source.0.cell, cold.source()[0].source.0.cell);
        let final_gas = next.chambers().iter().map(|c| (c.id, c.boundary.piston().clone())).collect();
        let work: f64 = next.chambers().iter().map(|c| c.boundary.step().piston.exterior_pressure_work_j).sum();
        assert!((storage(next.thermal(), next.circuit(), &final_gas)+work-initial).abs() < 4e-6);
        assert!(next.step().total_energy_residual_j.abs() <= controls().total_energy_tolerance_j);
    });
}

#[test]
fn aggregate_surface_overlap_unknown_cells_and_incomplete_bindings_refuse() {
    fixture(|binding, thermal, circuit| {
        let gas = gases(0.75); let drives = chambers(&gas, DT);
        let patches = contacts(1.0);
        let attempt = |patches: &[_]| network::step_candidate(binding, &thermal, DT, &circuit,
            &[], &BASES, &drives, patches, controls(), &[0.0; 2], || false);
        // Existing contact covers 1 m2 of each 6 m2 cell: two 3 m2 patches do not fit.
        let mut overlap = [patches[0], patches[2]];
        for patch in &mut overlap { patch.surface.area_m2 = 3.0; }
        assert!(attempt(&overlap).is_err());
        assert!(attempt(&[patches[0], patches[0]]).is_err());
        let mut unknown = patches[0]; unknown.contact.cell = CellId(999);
        assert!(attempt(&[unknown]).is_err());
        let mut regime = patches[0]; regime.surface.speed_m_s = 1e20;
        assert!(attempt(&[regime]).is_err());
        assert!(network::step_candidate(binding, &thermal, DT, &circuit, &[], &[],
            &drives, &patches, controls(), &[0.0; 2], || false).is_err());
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(circuit.charges_c(), &[10.0]);
    });
}

#[test]
fn failed_late_trial_and_cancellation_publish_no_subset_and_closed_faces_are_exact() {
    fixture(|binding, thermal, circuit| {
        let gas = gases(0.75); let saved = gas.clone(); let drives = chambers(&gas, DT);
        let cells = thermal.cells().clone(); let electrical = circuit.clone(); let patches = contacts(1.0);
        let mut limited = controls(); limited.maximum_boundary_solves = 1;
        assert!(network::step_candidate(binding, &thermal, DT, &circuit, &[], &BASES,
            &drives, &patches, limited, &[0.0; 2], || false).is_err());
        let mut polls = 0;
        network::step_candidate(binding, &thermal, DT, &circuit, &[], &BASES,
            &drives, &patches, controls(), &[0.0; 2], || { polls += 1; false }).unwrap();
        let mut at = 0;
        assert!(matches!(network::step_candidate(binding, &thermal, DT, &circuit, &[], &BASES,
            &drives, &patches, controls(), &[0.0; 2], || { at += 1; at == polls }),
            Err(CircuitEnthalpyError::Cancelled)));
        assert_eq!(gas, saved); assert_eq!(thermal.cells(), &cells); assert_eq!(circuit, electrical);
        let mut closed = patches;
        for patch in &mut closed { patch.surface.area_m2 = 0.0; patch.surface.speed_m_s = 0.0; }
        let a = network::step_candidate(binding, &thermal, DT, &circuit, &[], &BASES,
            &drives, &closed, controls(), &[0.0; 2], || false).unwrap();
        let b = network::step_candidate(binding, &thermal, DT, &circuit, &[], &BASES,
            &drives, &[], controls(), &[0.0; 2], || false).unwrap();
        assert_eq!(a.thermal().cells(), b.thermal().cells()); assert_eq!(a.circuit(), b.circuit());
        for (one, two) in a.chambers().iter().zip(b.chambers()) {
            assert_eq!(one.boundary.piston(), two.boundary.piston());
            assert!(one.boundary.source().iter().all(|s| s.source.1.evaluation.is_none()));
        }
    });
}

#[test]
fn authored_natural_convection_drives_the_same_electrical_and_gas_trajectory() {
    let forced = run_with_law(0.75, heat_law("1").unwrap(), 12, false).unwrap();
    let first = run_with_law(0.75, heat_law("natural").unwrap(), 12, false).unwrap();
    let second = run_with_law(0.25, heat_law("natural").unwrap(), 12, false).unwrap();
    let low_gravity = run_with_law(0.75, heat_law("natural:3.71").unwrap(), 12, false).unwrap();
    for outcome in [&forced, &first, &second, &low_gravity] {
        assert!(outcome.energy_defect_j.abs() < 2e-5);
        assert!(outcome.cell_k > 300.0);
    }
    assert_ne!(first.cell_k, forced.cell_k);
    assert_ne!(first.cell_k, second.cell_k);
    assert_ne!(first.charge_c, second.charge_c);
    assert_ne!(first.velocity_m_s, second.velocity_m_s);
    assert_ne!(first.sound_m_s, second.sound_m_s);
    assert!(first.cell_k > low_gravity.cell_k);
    for invalid in ["natural:0", "natural:-1", "natural:NaN", "0", "-1"] {
        assert!(heat_law(invalid).is_err());
    }
}
