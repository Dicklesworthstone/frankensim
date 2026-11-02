//! The complete source/caloric/thermocapillary owner selects native conservation.

use super::*;
use fs_lbm::d3q19::freesurface3::ConversionRedistribution3;

#[test]
fn g1_conservative_thermocapillary_trajectory_transports_original_heat_and_source_state() {
    for slope in [1e-6, -1e-6] {
        with_graph(slope, |original, gate| {
            let mut thermal = original
                .clone()
                .with_conversion_redistribution(ConversionRedistribution3::ConservativeLocal);
            let populations = original.flow().lattice().grid().state_bits();
            let mass = original.flow().ledger_mass_kg().unwrap();
            let enthalpy = original.total_enthalpy_j().unwrap();
            assert_eq!(thermal.flow().lattice().grid().state_bits(), populations);
            assert_eq!(thermal.total_enthalpy_j().unwrap(), enthalpy);
            assert_eq!(thermal.physical_time_s(), original.physical_time_s());
            assert_eq!(
                thermal.height_graph_thermocapillary_surface(),
                Some(graph())
            );
            for step in 1..=3 {
                let report = thermal.step(gate).unwrap();
                near(report.mass_after_kg, mass, 2e-12);
                near(report.enthalpy_after_j, enthalpy, 2e-11);
                assert!(report.mass_routes > 0);
                assert_eq!(report.external_heat_j, 0.0);
                assert_eq!(report.viscous_heat_j, 0.0);
                assert_eq!(report.body_heat_j, 0.0);
                assert_eq!(thermal.accepted_steps(), step);
                assert_eq!(thermal.physical_time_s(), step as f64 * DT);
                let conversion = thermal
                    .flow()
                    .lattice()
                    .conversion_momentum_report()
                    .unwrap();
                assert!(conversion.kinetic_energy_change.is_finite());
                assert!(
                    conversion
                        .absolute_momentum_defect
                        .iter()
                        .all(|p| *p < 1e-10)
                );
                assert_eq!(
                    thermal
                        .material_at([3, 3, 2])
                        .unwrap()
                        .material()
                        .card_identity(),
                    original
                        .material_at([3, 3, 2])
                        .unwrap()
                        .material()
                        .card_identity(),
                );
            }
            assert_ne!(thermal.flow().lattice().grid().state_bits(), populations);
        });
    }
}

#[test]
fn g4_conservative_thermal_refusal_retains_the_accepted_conversion_receipt_and_clock() {
    with_graph(1e-6, |original, _| {
        let mut thermal = original
            .clone()
            .with_conversion_redistribution(ConversionRedistribution3::ConservativeLocal);
        let gate = CancelGate::new_clock_free();
        thermal.step(&gate).unwrap();
        let populations = thermal.flow().lattice().grid().state_bits();
        let mass = thermal.flow().ledger_mass_kg().unwrap().to_bits();
        let enthalpy = thermal.total_enthalpy_j().unwrap().to_bits();
        let time = thermal.physical_time_s();
        let steps = thermal.accepted_steps();
        let receipt = thermal
            .flow()
            .lattice()
            .conversion_momentum_report()
            .copied();
        assert!(
            thermal
                .step_with_heat(
                    &[FreeSurfaceCellHeat {
                        coordinates: [3, 3, 3],
                        energy_j: 1000.0,
                    }],
                    &gate
                )
                .is_err()
        );
        assert_eq!(thermal.flow().lattice().grid().state_bits(), populations);
        assert_eq!(thermal.flow().ledger_mass_kg().unwrap().to_bits(), mass);
        assert_eq!(thermal.total_enthalpy_j().unwrap().to_bits(), enthalpy);
        assert_eq!(
            thermal.flow().lattice().conversion_momentum_report(),
            receipt.as_ref()
        );
        gate.request();
        assert!(thermal.step(&gate).is_err());
        assert_eq!(thermal.flow().lattice().grid().state_bits(), populations);
        assert_eq!(thermal.total_enthalpy_j().unwrap().to_bits(), enthalpy);
        assert_eq!(
            thermal.flow().lattice().conversion_momentum_report(),
            receipt.as_ref()
        );
        assert_eq!(thermal.physical_time_s(), time);
        assert_eq!(thermal.accepted_steps(), steps);
    });
}
