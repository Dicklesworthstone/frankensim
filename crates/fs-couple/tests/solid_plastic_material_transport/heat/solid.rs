//! Original 3D solve, source binding and closed accounting through one entry point.
#[path = "solid/dynamics.rs"]
mod dynamics;

use super::*;
use fs_couple::circuit_enthalpy::contacts::heat::{
    HeatAssembly, HeatAssemblyAdvanceError, HeatNetworkSample, solid::SolidHeatError,
};
use fs_couple::circuit_enthalpy::material::MaterialCell;
use fs_material::thermal::source::ThermalMaterialPoint;

#[test]
fn direct_heat_entry_preserves_native_states_and_closes_the_original_energy_tolerance() {
    with_heat(0.7, 0.1, |cx, _, _, thermal, solid, prepared| {
        let record = Recorder::default();
        let (loaded, audit) = solid.advance_heat_network_candidate(thermal, prepared, &record,
            1.0, &supports(SHEAR), &applied(), IncrementSettings::default(), thermal_config(),
            control(), cx).unwrap();
        // Exercise the public unsized provider path with the exact same native owner.
        let assembly: &dyn HeatAssembly<Sample = HeatNetworkSample, Error = HeatNetworkError> = prepared;
        let (assembled, assembly_audit) = solid.advance_heat_assembly_candidate(
            thermal, assembly, &record, 1.0, &supports(SHEAR), &applied(),
            IncrementSettings::default(), thermal_config(), control(), cx).unwrap();
        assert_eq!(assembled.thermal().cells(), loaded.thermal().cells());
        assert_eq!(assembled.solid().solution(), loaded.solid().solution());
        assert_eq!(assembled.observer().second, loaded.observer().second);
        assert_eq!(assembled.observer().first.sample().unwrap().radiation,
            loaded.observer().first.sample().unwrap().radiation);
        assert_eq!(assembly_audit, audit);
        assert!(assembled.report().physical_heat_j > 0.0);
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
fn complete_assembly_preserves_actual_inertial_plastic_state_and_native_heat_account() {
    with_heat(0.7, 0.1, |cx, _, _, thermal, solid, prepared| {
        let mut velocity = vec![[0.0; 3]; NODES.len()];
        velocity[3][0] = 1.0;
        let moving = fs_couple::solid_enthalpy::plastic::dynamics::ThermoplasticDynamics::try_new(
            solid.clone(), velocity.clone(), cx,
        ).unwrap();
        let record = Recorder::default();
        let (native, native_audit) = moving.advance_heat_network_candidate(
            thermal, prepared, &record, 1.0, &supports(SHEAR), &applied(),
            IncrementSettings::default(), thermal_config(), control(), cx,
        ).unwrap();
        let assembly: &dyn HeatAssembly<Sample = HeatNetworkSample, Error = HeatNetworkError> = prepared;
        let (complete, audit) = moving.advance_heat_assembly_candidate(
            thermal, assembly, &record, 1.0, &supports(SHEAR), &applied(),
            IncrementSettings::default(), thermal_config(), control(), cx,
        ).unwrap();
        assert_eq!(complete.thermal().cells(), native.thermal().cells());
        assert_eq!(complete.dynamics().solid().solution(), native.dynamics().solid().solution());
        assert_eq!(complete.dynamics().velocity_mps(), native.dynamics().velocity_mps());
        assert_eq!(complete.dynamics().time_s(), native.dynamics().time_s());
        assert_eq!(complete.observer().second, native.observer().second);
        assert_eq!(complete.observer().first.sample().unwrap().radiation,
            native.observer().first.sample().unwrap().radiation);
        assert_eq!(audit, native_audit);
        let actual = complete.report();
        let expected = native.report();
        assert_eq!(actual.inertia.kinetic_before_j, expected.inertia.kinetic_before_j);
        assert_eq!(actual.inertia.kinetic_after_j, expected.inertia.kinetic_after_j);
        assert_eq!(actual.inertia.numerical_loss_j, expected.inertia.numerical_loss_j);
        assert_eq!(actual.material.physical_heat_j, expected.material.physical_heat_j);
        assert!(actual.inertia.kinetic_before_j > 0.0);
        assert!(actual.inertia.kinetic_after_j > 0.0);
        assert!(actual.inertia.numerical_loss_j > 0.0);
        assert!(actual.material.physical_heat_j > 0.0);
        assert!(complete.dynamics().solid().solution().element_states[0].alpha > 0.0);
        close(
            actual.material.thermal_audit.stored_energy_change_j,
            actual.material.physical_heat_j + actual.material.reversible_heat_j + audit.external_heat_j,
            control().energy_tolerance_j,
        );
        assert_eq!(moving.velocity_mps(), velocity.as_slice());
        assert_eq!(moving.time_s(), 0.0);
        assert_eq!(thermal.time_s(), 0.0);
        assert!(record.durations.is_empty());
    });
}

#[test]
fn cyclic_history_reaches_complete_material_heat_networks_and_unsized_assemblies() {
    use fs_material::plastic::chaboche::source::BackstressProperties;
    use fs_material::plastic::voce::source::{
        VOCE_SATURATION_INCREMENT_PROPERTY as Q, VOCE_SATURATION_RATE_PROPERTY as B,
    };
    let schema = || [
        BackstressProperties::new("fast_backstress_C", "fast_backstress_gamma"),
        BackstressProperties::new("slow_backstress_C", "slow_backstress_gamma"),
    ];
    let f = fixture_with_properties(0.7, &[
        scalar(Q, Pressure::DIMS, 1e6),
        scalar(B, Dims([0; 6]), 200.0),
        scalar("fast_backstress_C", Pressure::DIMS, 2e8),
        scalar("fast_backstress_gamma", Dims([0; 6]), 80.0),
        scalar("slow_backstress_C", Pressure::DIMS, 5e7),
        scalar("slow_backstress_gamma", Dims([0; 6]), 10.0),
    ]);
    with_heat_fixture(f, 0.1, |cx, _, _, thermal, original, network| {
        let prepared = network.prepare(cx, thermal).unwrap();
        let cyclic = original.clone()
            .with_voce_hardening(&vec![Some(MaterialPropertySelection::SingleClaimOnly); 2], cx)
            .unwrap()
            .with_chaboche_branches(&[0, 1].map(|_| (
                MaterialPropertySelection::SingleClaimOnly, schema(),
            )), cx).unwrap();
        let assembly: &dyn HeatAssembly<Sample = HeatNetworkSample, Error = HeatNetworkError> = &prepared;
        let observer = prepared.observe(cx, (), control().energy_tolerance_j).unwrap();
        let native = cyclic.advance(thermal, &observer, 1.0, &supports(SHEAR), &applied(),
            IncrementSettings::default(), thermal_config(), control(), cx).unwrap();
        let (direct, direct_audit) = cyclic.advance_heat_network_candidate(
            thermal, &prepared, &(), 1.0, &supports(SHEAR), &applied(),
            IncrementSettings::default(), thermal_config(), control(), cx).unwrap();
        let (complete, audit) = cyclic.advance_heat_assembly_candidate(
            thermal, assembly, &(), 1.0, &supports(SHEAR), &applied(),
            IncrementSettings::default(), thermal_config(), control(), cx).unwrap();
        assert_eq!(native.thermal().cells(), direct.thermal().cells());
        assert_eq!(native.solid().solution(), direct.solid().solution());
        assert_eq!(complete.thermal().cells(), direct.thermal().cells());
        assert_eq!(complete.solid().solution(), direct.solid().solution());
        assert_eq!(audit, direct_audit);
        let loaded = &complete.solid().solution().element_states[0];
        assert!(loaded.plastic.alpha > 0.0);
        assert!(loaded.backstress_pa.iter().all(|branch| branch.iter().any(|x| *x != 0.0)));
        assert!(complete.report().physical_heat_j > 0.0);
        assert!(audit.combined_residual_j.abs() <= control().energy_tolerance_j);
        assert_eq!(cyclic.time_s(), 0.0);

        let mut velocity = vec![[0.0; 3]; NODES.len()];
        velocity[3][0] = 1.0;
        let moving = fs_couple::solid_enthalpy::plastic::dynamics::ThermoplasticDynamics::try_new(
            cyclic, velocity, cx).unwrap();
        let (direct, direct_audit) = moving.advance_heat_network_candidate(
            thermal, &prepared, &(), 1.0, &supports(SHEAR), &applied(),
            IncrementSettings::default(), thermal_config(), control(), cx).unwrap();
        let (complete, audit) = moving.advance_heat_assembly_candidate(
            thermal, assembly, &(), 1.0, &supports(SHEAR), &applied(),
            IncrementSettings::default(), thermal_config(), control(), cx).unwrap();
        assert_eq!(complete.thermal().cells(), direct.thermal().cells());
        assert_eq!(complete.dynamics().solid().solution(), direct.dynamics().solid().solution());
        assert_eq!(complete.dynamics().velocity_mps(), direct.dynamics().velocity_mps());
        assert_eq!(complete.report().inertia, direct.report().inertia);
        assert_eq!(audit, direct_audit);
        assert!(complete.report().inertia.kinetic_before_j > 0.0);
        assert!(complete.report().material.physical_heat_j > 0.0);
        let history = &complete.dynamics().solid().solution().element_states[0];
        assert!(history.plastic.alpha > 0.0);
        assert!(history.backstress_pa.iter().all(|branch| branch.iter().any(|x| *x != 0.0)));
        close(complete.report().material.thermal_audit.stored_energy_change_j,
            complete.report().material.physical_heat_j + complete.report().material.reversible_heat_j
                + audit.external_heat_j, control().energy_tolerance_j);
        assert!(audit.combined_residual_j.abs() <= control().energy_tolerance_j);
        assert_eq!(moving.time_s(), 0.0);
        assert_eq!(thermal.time_s(), 0.0);
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
        assert!(matches!(solid.advance_heat_assembly_candidate(thermal, &prepared, &(),
            1.0, &supports(SHEAR), &applied(), IncrementSettings::default(), thermal_config(),
            control(), cx), Err(HeatAssemblyAdvanceError::Binding(SolidHeatError::Invalid(
                "heat and solid material conditions disagree")))));
        let moving = fs_couple::solid_enthalpy::plastic::dynamics::ThermoplasticDynamics::try_new(
            solid.clone(), vec![[0.0; 3]; NODES.len()], cx).unwrap();
        let velocity = moving.velocity_mps().to_vec();
        assert!(matches!(moving.advance_heat_network_candidate(thermal, &prepared, &(),
            1.0, &supports(SHEAR), &applied(), IncrementSettings::default(), thermal_config(),
            control(), cx), Err(SolidHeatError::Invalid("heat and solid material conditions disagree"))));
        assert!(matches!(moving.advance_heat_assembly_candidate(thermal, &prepared, &(),
            1.0, &supports(SHEAR), &applied(), IncrementSettings::default(), thermal_config(),
            control(), cx), Err(HeatAssemblyAdvanceError::Binding(SolidHeatError::Invalid(
                "heat and solid material conditions disagree")))));
        assert_eq!(moving.velocity_mps(), velocity.as_slice());
        assert_eq!(thermal.cells(), &saved); assert_eq!(solid.solution(), &history);
        gate.request();
        assert!(matches!(solid.advance_heat_network_candidate(thermal, original, &(),
            1.0, &supports(SHEAR), &applied(), IncrementSettings::default(), thermal_config(),
            control(), cx), Err(SolidHeatError::Cancelled)));
        assert!(matches!(moving.advance_heat_network_candidate(thermal, original, &(),
            1.0, &supports(SHEAR), &applied(), IncrementSettings::default(), thermal_config(),
            control(), cx), Err(SolidHeatError::Cancelled)));
        assert!(matches!(solid.advance_heat_assembly_candidate(thermal, original, &(),
            1.0, &supports(SHEAR), &applied(), IncrementSettings::default(), thermal_config(),
            control(), cx), Err(HeatAssemblyAdvanceError::Binding(SolidHeatError::Cancelled))));
        assert!(matches!(moving.advance_heat_assembly_candidate(thermal, original, &(),
            1.0, &supports(SHEAR), &applied(), IncrementSettings::default(), thermal_config(),
            control(), cx), Err(HeatAssemblyAdvanceError::Binding(SolidHeatError::Cancelled))));
        assert_eq!(moving.velocity_mps(), velocity.as_slice());
        assert_eq!(thermal.cells(), &saved); assert_eq!(solid.solution(), &history);
    });
}
