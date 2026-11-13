//! G1/G3 moving phase sources enter the actual finite gas/wall thermal owner.
//! G4 late refusal preserves all accepted storage and momentum inventories.
use super::*;
use fs_couple::liquid_gas::evaporation::momentum::{
    MovingGasParcel,
    conduction::{
        advance_liquid_surfaces_with_moving_gas_and_conduction,
        advance_liquid_surfaces_with_moving_gas_and_radiative_conduction,
    },
    multiple::advance_liquid_surfaces_with_moving_gas,
};
use fs_lbm::material::free_surface::enthalpy::conduction::FreeSurfaceRadiationFace;

#[test]
fn opposing_moving_phase_streams_retain_momentum_through_finite_wall_gas_and_radiation() {
    let cold = [3, 2, 5];
    let warm = [3, 4, 2];
    with_wall(|body, contact| {
        with_liquid_interfaces(
            0.0,
            &[(cold, 0.4, 295.0), (warm, 0.4, 305.0)],
            |liquid, gate| {
                let source = source(-13.4e6);
                let gas = MovingGasParcel::new(gas(&source, 300.0, 3500.0), [0.025, -0.015, 0.005])
                    .unwrap();
                let original_gas = gas.clone();
                let original_bits = liquid.flow().lattice().grid().state_bits();
                let binding = binding(&source, 0.0);
                let surfaces = [cold, warm].map(|liquid| LiquidGasEvaporationSurface {
                    liquid,
                    gas: film(liquid).gas,
                    ..surface(&Phase, &binding)
                });
                let bodies = [body];
                let initial_wall = bodies[0].enthalpy_j();
                let phase = advance_liquid_surfaces_with_moving_gas(
                    liquid,
                    &gas,
                    &source,
                    &surfaces,
                    joint_controls(),
                    None,
                    gate,
                )
                .unwrap();
                let next = advance_liquid_surfaces_with_moving_gas_and_radiative_conduction(
                    liquid,
                    &gas,
                    &source,
                    &surfaces,
                    &bodies,
                    &[FreeSurfaceConductionFace {
                        cells: [cold, [2, 2, 5]],
                        area_m2: 1e-6,
                        distances_m: [0.0005; 2],
                    }],
                    &[FreeSurfaceCellHeat {
                        coordinates: [1, 3, 3],
                        energy_j: 2e-6,
                    }],
                    &[contact],
                    &[],
                    &[film(cold), film(warm)],
                    &[FreeSurfaceRadiationFace {
                        liquid: warm,
                        gas: film(warm).gas,
                        area_m2: 1e-6,
                        emissivity: 0.6.into(),
                        surroundings_k: 450.0,
                    }],
                    joint_controls(),
                    thermal_controls(),
                    controls().phase.gas.gas,
                    None,
                    true,
                    gate,
                )
                .unwrap();
                let report = next.report();
                let native = &report.native;
                let motion = &report.momentum;
                let thermal = &native.thermal.conduction;
                assert!(motion.condensed_mass_kg > 0.0);
                assert!(motion.evaporated_mass_kg > 0.0);
                assert!(motion.liquid_mixing_heat_j > 0.0);
                assert!(motion.gas_mixing_heat_j > 0.0);
                assert_eq!(next.gas().momentum_kg_m_s(), phase.gas().momentum_kg_m_s());
                assert_eq!(
                    next.gas().parcel().mass_kg().to_bits(),
                    phase.gas().parcel().mass_kg().to_bits()
                );
                assert_eq!(
                    next.gas().kinetic_energy_j().to_bits(),
                    phase.gas().kinetic_energy_j().to_bits()
                );
                assert_eq!(
                    native.thermal_gas_temperature_k[0],
                    phase.gas().parcel().state().state().temperature
                );
                assert_eq!(
                    native.thermal_gas_temperature_k[1],
                    next.gas().parcel().state().state().temperature
                );
                assert_eq!(native.phase.flow, thermal.flow);
                assert_eq!(next.liquid().accepted_steps(), 1);
                assert_eq!(next.liquid().physical_time_s(), DT);
                assert!(thermal.thermal.substeps >= 4);
                assert!(thermal.radiation.as_ref().unwrap().heat_to_liquid_j > 0.0);
                assert_ne!(next.bodies()[0].enthalpy_j(), initial_wall);
                for d in 0..3 {
                    near(
                        native.phase.mass_exchange.momentum_change_kg_m_s[d]
                            + next.gas().momentum_kg_m_s()[d]
                            - gas.momentum_kg_m_s()[d],
                        0.0,
                        controls().native.momentum_tolerance_kg_m_s,
                    );
                }
                let dm = next.liquid().flow().ledger_mass_kg().unwrap()
                    - liquid.flow().ledger_mass_kg().unwrap();
                let dh =
                    next.liquid().total_enthalpy_j().unwrap() - liquid.total_enthalpy_j().unwrap();
                let du = next.gas().parcel().mass_kg()
                    * next.gas().parcel().state().specific_internal_energy_j_kg()
                    - gas.parcel().mass_kg() * gas.parcel().state().specific_internal_energy_j_kg();
                let dw = next.bodies()[0].enthalpy_j() - initial_wall;
                let dk = next.gas().kinetic_energy_j() - gas.kinetic_energy_j();
                near(
                    dm + next.gas().parcel().mass_kg() - gas.parcel().mass_kg(),
                    0.0,
                    controls().mass_tolerance_kg,
                );
                near(dk, motion.gas_kinetic_energy_change_j, 1e-17);
                near(
                    dh + binding.liquid_enthalpy_offset_j_kg * dm
                        + du
                        + dw
                        + dk
                        + motion.liquid_kinetic_energy_change_j,
                    2e-6 + thermal.flow.viscous_heat_j
                        + thermal.radiation.as_ref().unwrap().heat_to_liquid_j
                        + thermal.body_radiation.as_ref().unwrap().heat_to_solids_j,
                    controls().energy_tolerance_j,
                );
                assert!(native.absolute_energy_defect_j <= joint_controls().energy_tolerance_j);
                assert_eq!(gas, original_gas);
                assert_eq!(liquid.flow().lattice().grid().state_bits(), original_bits);
                assert_eq!(bodies[0].enthalpy_j(), initial_wall);
            },
        );
    });
}

#[test]
fn moving_gas_pressure_uses_post_thermal_eos_without_a_second_source_impulse() {
    with_liquid_pressure(0.0, 0.4, Some(100000.0), |liquid, gate| {
        let source = source(-13.4e6);
        let initial = gas(&source, 320.0, 6000.0);
        let gas = MovingGasParcel::new(
            ClosedGasParcel::new(initial.model().clone(), initial.point().clone(), 3.686e-7)
                .unwrap(),
            [0.025, -0.015, 0.005],
        )
        .unwrap();
        let binding = binding(&source, 0.0);
        let surfaces = [surface(&Phase, &binding)];
        let phase = advance_liquid_surfaces_with_moving_gas(
            liquid,
            &gas,
            &source,
            &surfaces,
            joint_controls(),
            None,
            gate,
        )
        .unwrap();
        let advance = |pressure| {
            advance_liquid_surfaces_with_moving_gas_and_conduction(
                liquid,
                &gas,
                &source,
                &surfaces,
                &[],
                &[],
                &[],
                &[],
                &[],
                &[film(INTERFACE)],
                joint_controls(),
                thermal_controls(),
                controls().phase.gas.gas,
                pressure,
                false,
                gate,
            )
            .unwrap()
        };
        let passive = advance(None);
        let loaded = advance(Some(0.05));
        let report = &loaded.report().native;
        let pressure = report
            .thermal
            .conduction
            .finite_gas
            .as_ref()
            .unwrap()
            .pressure
            .unwrap();
        assert_eq!(
            pressure.initial_pressure_pa,
            phase.gas().parcel().state().state().pressure
        );
        assert_eq!(
            pressure.applied_pressure_pa,
            loaded.gas().parcel().state().state().pressure
        );
        assert_eq!(report.phase.pressure, Some(pressure));
        assert_eq!(loaded.gas(), passive.gas());
        assert_eq!(
            loaded.gas().momentum_kg_m_s(),
            phase.gas().momentum_kg_m_s()
        );
        assert_eq!(
            loaded.gas().kinetic_energy_j(),
            phase.gas().kinetic_energy_j()
        );
        assert_ne!(
            loaded.liquid().flow().lattice().grid().state_bits(),
            passive.liquid().flow().lattice().grid().state_bits()
        );
        assert_eq!(loaded.liquid().accepted_steps(), 1);
        assert_eq!(loaded.liquid().physical_time_s(), DT);
    });
}

#[test]
fn late_thermal_refusal_discards_mixed_momentum_heat_and_finite_walls_together() {
    with_wall(|body, contact| {
        with_liquid(0.0, 0.4, |liquid, gate| {
            let source = source(-13.4e6);
            let gas =
                MovingGasParcel::new(gas(&source, 320.0, 6000.0), [0.025, -0.015, 0.005]).unwrap();
            let original_gas = gas.clone();
            let original_bits = liquid.flow().lattice().grid().state_bits();
            let original_h = liquid.total_enthalpy_j().unwrap();
            let binding = binding(&source, 0.0);
            let bodies = [body];
            let initial_wall = bodies[0].enthalpy_j();
            let mut exhausted = thermal_controls();
            exhausted.maximum_substeps = 1;
            let advance = |thermal| {
                advance_liquid_surfaces_with_moving_gas_and_conduction(
                    liquid,
                    &gas,
                    &source,
                    &[surface(&Phase, &binding)],
                    &bodies,
                    &[],
                    &[],
                    &[contact],
                    &[],
                    &[film(INTERFACE)],
                    joint_controls(),
                    thermal,
                    controls().phase.gas.gas,
                    None,
                    false,
                    gate,
                )
            };
            assert!(matches!(
                advance(exhausted),
                Err(EvaporatingLiquidGasError::Liquid(_))
            ));
            assert_eq!(gas, original_gas);
            assert_eq!(liquid.flow().lattice().grid().state_bits(), original_bits);
            assert_eq!(liquid.total_enthalpy_j().unwrap(), original_h);
            assert_eq!(liquid.accepted_steps(), 0);
            assert_eq!(bodies[0].enthalpy_j(), initial_wall);
            let next = advance(thermal_controls()).unwrap();
            assert_eq!(next.liquid().accepted_steps(), 1);
            assert_ne!(next.gas().momentum_kg_m_s(), gas.momentum_kg_m_s());
            assert_ne!(next.bodies()[0].enthalpy_j(), initial_wall);
        });
    });
}
