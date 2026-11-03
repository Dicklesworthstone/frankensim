//! Original cell-temperature radiation joins phase, finite gas and native flow.
use super::*;
use fs_couple::liquid_gas::evaporation::conduction::advance_liquid_surfaces_with_gas_and_radiative_conduction;
use fs_lbm::material::free_surface::enthalpy::conduction::{
    FreeSurfaceEmissivity, FreeSurfaceRadiationFace,
};

fn radiating_face(emissivity: FreeSurfaceEmissivity) -> FreeSurfaceRadiationFace {
    FreeSurfaceRadiationFace {
        liquid: INTERFACE,
        gas: [4, 3, 3],
        area_m2: 1e-6,
        emissivity,
        surroundings_k: 450.0,
    }
}

#[test]
fn g1_liquid_radiation_joins_phase_gas_walls_and_actual_endpoint_pressure_once() {
    for maximum_density_excursion in [None, Some(0.3)] {
        with_wall(|body, contact| {
            with_liquid_pressure(
                0.0,
                0.4,
                maximum_density_excursion.map(|_| 100000.0),
                |liquid, gate| {
                    let source = source(-13.4e6);
                    let gas = gas(&source, 300.0, 1000.0);
                    let original_gas = gas.clone();
                    let original_bits = liquid.flow().lattice().grid().state_bits();
                    let binding = binding(&source, 0.0);
                    let surfaces = [surface(&Phase, &binding)];
                    let bodies = [body];
                    let initial_wall_h = bodies[0].enthalpy_j();
                    let radiation = [radiating_face(0.6.into())];
                    let heat = [FreeSurfaceCellHeat {
                        coordinates: [1, 3, 3],
                        energy_j: 2e-6,
                    }];
                    let next = advance_liquid_surfaces_with_gas_and_radiative_conduction(
                        liquid,
                        &gas,
                        &source,
                        &surfaces,
                        &bodies,
                        &[],
                        &heat,
                        &[contact],
                        &[],
                        &[film(INTERFACE)],
                        &radiation,
                        joint_controls(),
                        thermal_controls(),
                        controls().phase.gas.gas,
                        maximum_density_excursion,
                        false,
                        gate,
                    )
                    .unwrap();
                    let report = next.report();
                    let thermal = &report.thermal.conduction;
                    let radiative = thermal.radiation.as_ref().unwrap();
                    let wall_radiation = thermal.body_radiation.as_ref().unwrap();
                    assert!(report.phase.phase.transferred_mass_kg > 0.0);
                    assert!(radiative.heat_to_liquid_j > 0.0);
                    assert!(wall_radiation.heat_to_solids_j > 0.0);
                    assert_ne!(thermal.environment.as_ref().unwrap().heat_to_liquid_j, 0.0);
                    assert_eq!(radiative.last_applied[0].face, radiation[0]);
                    let sample = &radiative.last_applied[0];
                    let expected_power = 0.6
                        * 5.670374419e-8
                        * 1e-6
                        * (450.0_f64.powi(4) - sample.liquid_temperature_k.powi(4));
                    near(sample.heat_into_liquid_w, expected_power, 1e-16);
                    inventory(
                        liquid,
                        &gas,
                        &binding,
                        &bodies,
                        &next,
                        heat[0].energy_j
                            + radiative.heat_to_liquid_j
                            + wall_radiation.heat_to_solids_j,
                    );
                    assert_eq!(next.liquid().accepted_steps(), 1);
                    assert_eq!(next.liquid().physical_time_s(), DT);
                    let finite = thermal.finite_gas.as_ref().unwrap();
                    match maximum_density_excursion {
                        Some(maximum) => {
                            let pressure = finite.pressure.unwrap();
                            assert_eq!(
                                pressure.applied_pressure_pa,
                                next.gas().state().state().pressure
                            );
                            assert_eq!(pressure.maximum_density_excursion, maximum);
                            assert_eq!(report.phase.pressure, Some(pressure));
                        }
                        None => assert!(finite.pressure.is_none()),
                    }
                    assert_eq!(gas, original_gas);
                    assert_eq!(liquid.flow().lattice().grid().state_bits(), original_bits);
                    assert_eq!(liquid.accepted_steps(), 0);
                    assert_eq!(bodies[0].enthalpy_j(), initial_wall_h);
                },
            );
        });
    }
}

#[test]
fn g1_radiative_supply_can_heat_liquid_without_an_invented_finite_gas_debit() {
    with_liquid(0.0, 0.4, |liquid, gate| {
        let source = source(-13.4e6);
        let closed = gas(&source, 300.0, 1000.0);
        let mut gas = closed.clone().into_isochoric();
        let mut next = liquid.clone();
        let before_h = next.total_enthalpy_j().unwrap();
        let radiation = [radiating_face(0.6.into())];
        let report = next
            .step_with_radiative_finite_gas_environment(
                &[],
                &[],
                thermal_controls(),
                &mut [],
                &[],
                &[],
                &mut gas,
                &[],
                &radiation,
                controls().phase.gas.gas,
                None,
                false,
                gate,
            )
            .unwrap();
        let thermal = &report.conduction;
        assert_eq!(thermal.thermal.substeps, 4);
        // Independently integrate the existing explicit fourth-power law.
        let capacity_j_k = 0.4 * 1000.0 * 0.001_f64.powi(3) * CP;
        let mut temperature = 300.0_f64;
        let mut supplied = 0.0;
        for _ in 0..4 {
            let heat = (DT / 4.0)
                * 0.6
                * 5.670374419e-8
                * 1e-6
                * (450.0_f64.powi(4) - temperature.powi(4));
            supplied += heat;
            temperature += heat / capacity_j_k;
        }
        let radiant = thermal.radiation.as_ref().unwrap();
        near(radiant.heat_to_liquid_j, supplied, 1e-14);
        near(
            radiant.thermal_endpoint[0].liquid_temperature_k,
            temperature,
            1e-10,
        );
        near(next.total_enthalpy_j().unwrap() - before_h, supplied, 1e-11);
        assert_eq!(ClosedGasParcel::from_isochoric(gas), closed);
        assert_eq!(
            thermal
                .finite_gas
                .as_ref()
                .unwrap()
                .internal_energy_change_j,
            0.0
        );
        assert!(thermal.environment.is_none());
        assert_eq!(next.accepted_steps(), 1);
        assert_eq!(liquid.accepted_steps(), 0);
    });
}

#[test]
fn g3_zero_radiation_preserves_the_old_route_and_missing_original_emissivity_rolls_back() {
    with_wall(|body, contact| {
        with_liquid(0.0, 0.4, |liquid, gate| {
            let source = source(-13.4e6);
            let gas = gas(&source, 300.0, 1000.0);
            let original_gas = gas.clone();
            let binding = binding(&source, 0.0);
            let surfaces = [surface(&Phase, &binding)];
            let bodies = [body];
            let initial_wall_h = bodies[0].enthalpy_j();
            let initial_h = liquid.total_enthalpy_j().unwrap();
            let bits = liquid.flow().lattice().grid().state_bits();
            let old = advance_liquid_surfaces_with_gas_and_conduction(
                liquid,
                &gas,
                &source,
                &surfaces,
                &bodies,
                &[],
                &[],
                &[contact],
                &[],
                &[film(INTERFACE)],
                joint_controls(),
                thermal_controls(),
                controls().phase.gas.gas,
                None,
                false,
                gate,
            )
            .unwrap();
            let run = |emissivity| {
                advance_liquid_surfaces_with_gas_and_radiative_conduction(
                    liquid,
                    &gas,
                    &source,
                    &surfaces,
                    &bodies,
                    &[],
                    &[],
                    &[contact],
                    &[],
                    &[film(INTERFACE)],
                    &[radiating_face(emissivity)],
                    joint_controls(),
                    thermal_controls(),
                    controls().phase.gas.gas,
                    None,
                    false,
                    gate,
                )
            };
            let reflecting = run(0.0.into()).unwrap();
            assert_eq!(reflecting.gas(), old.gas());
            assert_eq!(
                reflecting.liquid().flow().lattice().grid().state_bits(),
                old.liquid().flow().lattice().grid().state_bits()
            );
            assert_eq!(
                reflecting.liquid().total_enthalpy_j().unwrap(),
                old.liquid().total_enthalpy_j().unwrap()
            );
            assert_eq!(
                reflecting.bodies()[0].enthalpy_j(),
                old.bodies()[0].enthalpy_j()
            );
            assert_eq!(
                reflecting
                    .report()
                    .thermal
                    .conduction
                    .radiation
                    .as_ref()
                    .unwrap()
                    .heat_to_liquid_j,
                0.0
            );
            // The original bundle supplies no emissivity. The new connection
            // must retain that late native source refusal after phase transfer.
            assert!(matches!(
                run(FreeSurfaceEmissivity::Material),
                Err(EvaporatingLiquidGasError::Liquid(_))
            ));
            assert_eq!(gas, original_gas);
            assert_eq!(liquid.flow().lattice().grid().state_bits(), bits);
            assert_eq!(liquid.total_enthalpy_j().unwrap(), initial_h);
            assert_eq!(liquid.accepted_steps(), 0);
            assert_eq!(bodies[0].enthalpy_j(), initial_wall_h);
            let retry = run(0.0.into()).unwrap();
            assert_eq!(retry.gas(), reflecting.gas());
            assert_eq!(
                retry.liquid().flow().lattice().grid().state_bits(),
                reflecting.liquid().flow().lattice().grid().state_bits()
            );
        });
    });
}
