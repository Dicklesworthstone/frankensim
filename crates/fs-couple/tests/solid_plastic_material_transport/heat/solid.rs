//! Original 3D solve, source binding and closed accounting through one entry point.
#[path = "solid/dynamics.rs"]
mod dynamics;

use super::*;
use fs_couple::circuit_enthalpy::contacts::heat::solid::SolidHeatError;
use fs_couple::circuit_enthalpy::material::MaterialCell;
use fs_material::thermal::source::ThermalMaterialPoint;

#[test]
fn direct_heat_entry_preserves_native_states_and_closes_the_original_energy_tolerance() {
    with_heat(0.7, 0.1, |cx, _, _, thermal, solid, prepared| {
        let record = Recorder::default();
        let (loaded, audit) = solid.advance_heat_network_candidate(thermal, prepared, &record,
            1.0, &supports(SHEAR), &applied(), IncrementSettings::default(), thermal_config(),
            control(), cx).unwrap();
        let observer = prepared.observe(cx, record.clone(), control().energy_tolerance_j).unwrap();
        let manual = solid.advance(thermal, &observer, 1.0, &supports(SHEAR), &applied(),
            IncrementSettings::default(), thermal_config(), control(), cx).unwrap();
        assert_eq!(loaded.thermal().cells(), manual.thermal().cells());
        assert_eq!(loaded.solid().solution(), manual.solid().solution());
        assert_eq!(loaded.observer().second, manual.observer().second);
        assert_eq!(loaded.observer().first.sample().unwrap().radiation,
            manual.observer().first.sample().unwrap().radiation);
        assert_eq!(audit.external_heat_j,
            loaded.report().external_heat_j - audit.closed_radiation_net_j);
        assert_eq!(audit.combined_residual_j,
            loaded.report().combined_residual_j + audit.closed_radiation_net_j);
        assert!(audit.combined_residual_j.abs() <= control().energy_tolerance_j);
        assert!(audit.thermal_residual_j.abs() <= control().energy_tolerance_j);
        close(audit.external_heat_j, 0.0, 1e-10);
        let history = loaded.solid().solution().element_states.clone();
        // Keep only the caller's original history; the direct entry point adds
        // the one material heat owner, not another copy of the accepted pair.
        let (unloaded, second) = loaded.solid().advance_heat_network_candidate(
            loaded.thermal(), prepared, &loaded.observer().second, 2.0,
            &supports(SHEAR - 0.0005), &[], IncrementSettings::default(), thermal_config(),
            control(), cx).unwrap();
        assert_eq!(unloaded.solid().solution().element_states, history);
        assert!(second.combined_residual_j.abs() <= control().energy_tolerance_j);
        close(unloaded.observer().second.durations.iter().sum(), 2.0, 1e-12);
        assert_eq!(solid.time_s(), 0.0); assert_eq!(thermal.time_s(), 0.0);
        assert!(record.durations.is_empty());
    });
}

#[test]
fn same_card_and_caloric_law_do_not_excuse_different_mechanical_source_conditions() {
    with_heat(0.7, 0.1, |cx, gate, f, thermal, solid, original| {
        let keys = vec![Some((f.conductivity.clone(), MaterialPropertySelection::SingleClaimOnly)); 2];
        let sources = solid.thermal_material_cells(&keys).unwrap();
        // Both pressures are inside the conductivity claim. The independently
        // prepared heat network is valid, but describes different body conditions.
        let different = f.point.clone().with_quantity("p", pressure_quantity(), PRESSURE_PA + 0.5).unwrap();
        let cells = sources.iter().map(|s| MaterialCell {
            cell: s.cell, material_id: s.material_id, curve: s.curve,
            point: ThermalMaterialPoint { context: &different, ..s.point },
            conductivity: s.conductivity.clone(),
        }).collect::<Vec<_>>();
        let links = links();
        let bulk = MaterialNetwork { cells: &cells, resistors: &[], links: &links };
        let electrical = ElectricalContactNetwork { bulk: &bulk, contacts: &[], electrical_bases: &[] };
        let network = MaterialHeatNetwork { electrical: &electrical, thermal_contacts: &[], radiation: None };
        let prepared = network.prepare(cx, thermal).unwrap();
        assert!(prepared.resolve(cx, thermal).is_ok());
        let saved = thermal.cells().clone(); let history = solid.solution().clone();
        assert!(matches!(solid.advance_heat_network_candidate(thermal, &prepared, &(),
            1.0, &supports(SHEAR), &applied(), IncrementSettings::default(), thermal_config(),
            control(), cx), Err(SolidHeatError::Invalid("heat and solid material conditions disagree"))));
        let moving = fs_couple::solid_enthalpy::plastic::dynamics::ThermoplasticDynamics::try_new(
            solid.clone(), vec![[0.0; 3]; NODES.len()], cx).unwrap();
        let velocity = moving.velocity_mps().to_vec();
        assert!(matches!(moving.advance_heat_network_candidate(thermal, &prepared, &(),
            1.0, &supports(SHEAR), &applied(), IncrementSettings::default(), thermal_config(),
            control(), cx), Err(SolidHeatError::Invalid("heat and solid material conditions disagree"))));
        assert_eq!(moving.velocity_mps(), velocity.as_slice());
        assert_eq!(thermal.cells(), &saved); assert_eq!(solid.solution(), &history);
        gate.request();
        assert!(matches!(solid.advance_heat_network_candidate(thermal, original, &(),
            1.0, &supports(SHEAR), &applied(), IncrementSettings::default(), thermal_config(),
            control(), cx), Err(SolidHeatError::Cancelled)));
        assert!(matches!(moving.advance_heat_network_candidate(thermal, original, &(),
            1.0, &supports(SHEAR), &applied(), IncrementSettings::default(), thermal_config(),
            control(), cx), Err(SolidHeatError::Cancelled)));
        assert_eq!(moving.velocity_mps(), velocity.as_slice());
        assert_eq!(thermal.cells(), &saved); assert_eq!(solid.solution(), &history);
    });
}
