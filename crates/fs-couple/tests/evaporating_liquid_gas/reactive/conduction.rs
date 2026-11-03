//! Chemistry, phase transfer, finite walls and native thermal/flow publish together.
use super::*;
use crate::conduction::{film, joint_controls, thermal_controls, with_wall};
use fs_couple::liquid_gas::evaporation::conduction::advance_liquid_surfaces_with_gas_and_radiative_conduction;
use fs_lbm::material::free_surface::enthalpy::{
    FreeSurfaceCellHeat,
    conduction::{FreeSurfaceEmissivity, FreeSurfaceRadiationFace},
};

#[test]
fn g1_reactive_phase_and_radiative_conduction_match_original_owners_on_one_clock() {
    for maximum_density_excursion in [None, Some(0.3)] {
        with_liquid_pressure(
            0.0,
            0.4,
            maximum_density_excursion.map(|_| 100000.0),
            |liquid, gate| {
                with_wall(|body, contact| {
                    let law = reaction(500.0, 75.0);
                    let source = law.catalog();
                    let gas = gas(source, 300.0, 1000.0);
                    let original_gas = gas.clone();
                    let original_bits = liquid.flow().lattice().grid().state_bits();
                    let original_h = liquid.total_enthalpy_j().unwrap();
                    let binding = binding(source, 0.0);
                    let surfaces = [surface(&Phase, &binding)];
                    let exposed = [film(INTERFACE)];
                    let radiation = [FreeSurfaceRadiationFace {
                        liquid: INTERFACE,
                        gas: exposed[0].gas,
                        area_m2: 1e-6,
                        emissivity: FreeSurfaceEmissivity::Prescribed(0.6),
                        surroundings_k: 450.0,
                    }];
                    let bodies = [body];
                    let initial_wall_h = bodies[0].enthalpy_j();
                    let heat = [FreeSurfaceCellHeat {
                        coordinates: INTERFACE,
                        energy_j: 1e-6,
                    }];
                    let chemical = law
                        .parcel_candidate(&gas, DT, reactive_controls().reaction, || false)
                        .unwrap();
                    let direct = advance_liquid_surfaces_with_gas_and_radiative_conduction(
                        liquid,
                        chemical.gas(),
                        source,
                        &surfaces,
                        &bodies,
                        &[],
                        &heat,
                        &[contact],
                        &[],
                        &exposed,
                        &radiation,
                        joint_controls(),
                        thermal_controls(),
                        reactive_controls().reaction.gas,
                        maximum_density_excursion,
                        false,
                        gate,
                    )
                    .unwrap();
                    let next = law
                        .radiative_conducting_evaporating_liquid_candidate(
                            liquid,
                            &gas,
                            &surfaces,
                            &bodies,
                            &[],
                            &heat,
                            &[contact],
                            &[],
                            &exposed,
                            &radiation,
                            joint_controls(),
                            thermal_controls(),
                            maximum_density_excursion,
                            false,
                            reactive_controls(),
                            gate,
                        )
                        .unwrap();
                    let report = next.report();
                    assert!(report.reaction.extent_mol > 0.0);
                    assert_eq!(&report.reaction, chemical.step());
                    assert_eq!(next.gas(), direct.gas());
                    assert_eq!(
                        next.liquid().flow().lattice().grid().state_bits(),
                        direct.liquid().flow().lattice().grid().state_bits(),
                    );
                    assert_eq!(
                        next.liquid().total_enthalpy_j().unwrap(),
                        direct.liquid().total_enthalpy_j().unwrap()
                    );
                    assert_eq!(
                        next.bodies()[0].enthalpy_j(),
                        direct.bodies()[0].enthalpy_j()
                    );
                    assert_ne!(next.bodies()[0].enthalpy_j(), initial_wall_h);
                    assert!(core::ptr::eq(next.catalog(), source));
                    assert_eq!(report.duration_s, DT);
                    assert_eq!(next.liquid().accepted_steps(), 1);
                    assert_eq!(next.liquid().physical_time_s(), DT);
                    assert_eq!(next.gas().volume_m3(), gas.volume_m3());

                    let thermal = &report.phase.thermal.conduction;
                    let finite_gas = thermal.finite_gas.as_ref().unwrap();
                    assert_ne!(thermal.environment.as_ref().unwrap().heat_to_liquid_j, 0.0);
                    assert_ne!(thermal.bodies.as_ref().unwrap().heat_from_gas_j, 0.0);
                    let radiant_supply = thermal.body_radiation.as_ref().unwrap().heat_to_solids_j;
                    assert!(radiant_supply > 0.0);
                    let liquid_radiation = thermal.radiation.as_ref().unwrap().heat_to_liquid_j;
                    assert!(liquid_radiation > 0.0);
                    assert_eq!(
                        report.phase.thermal_gas_temperature_k[1],
                        next.gas().state().state().temperature,
                    );
                    if let Some(maximum) = maximum_density_excursion {
                        let pressure = finite_gas.pressure.unwrap();
                        assert_eq!(
                            pressure.applied_pressure_pa,
                            next.gas().state().state().pressure
                        );
                        assert_eq!(pressure.maximum_density_excursion, maximum);
                        assert_ne!(pressure.applied_pressure_pa, gas.state().state().pressure);
                    } else {
                        assert!(finite_gas.pressure.is_none());
                    }

                    // Independently close the whole original inventory: gas/liquid
                    // and wall/liquid heat are internal. Only the authored cell
                    // source and both prescribed radiative boundaries enter this case.
                    let dm = next.liquid().flow().ledger_mass_kg().unwrap()
                        - liquid.flow().ledger_mass_kg().unwrap();
                    let dh = next.liquid().total_enthalpy_j().unwrap() - original_h;
                    let du = next.gas().mass_kg()
                        * next.gas().state().specific_internal_energy_j_kg()
                        - gas.mass_kg() * gas.state().specific_internal_energy_j_kg();
                    let dw = next.bodies()[0].enthalpy_j() - initial_wall_h;
                    near(
                        dm + next.gas().mass_kg() - gas.mass_kg(),
                        0.0,
                        reactive_controls().mass_tolerance_kg,
                    );
                    near(
                        dh + binding.liquid_enthalpy_offset_j_kg * dm + du + dw,
                        heat[0].energy_j + radiant_supply + liquid_radiation,
                        reactive_controls().energy_tolerance_j,
                    );
                    let volatile = report
                        .reaction
                        .species
                        .iter()
                        .find(|row| row.species == water())
                        .unwrap();
                    near(
                        dm + next.gas().species_mass_kg(&water()).unwrap()
                            - gas.species_mass_kg(&water()).unwrap()
                            - M * volatile.stoichiometric_change_mol,
                        0.0,
                        reactive_controls().mass_tolerance_kg,
                    );
                    assert!(
                        report.absolute_mass_defect_kg <= reactive_controls().mass_tolerance_kg
                    );
                    assert!(
                        report.absolute_energy_defect_j <= reactive_controls().energy_tolerance_j
                    );

                    assert_eq!(gas, original_gas);
                    assert_eq!(liquid.flow().lattice().grid().state_bits(), original_bits);
                    assert_eq!(liquid.total_enthalpy_j().unwrap(), original_h);
                    assert_eq!(liquid.accepted_steps(), 0);
                    assert_eq!(bodies[0].enthalpy_j(), initial_wall_h);
                    let (_, _, catalog, accepted_bodies, accepted_report) = next.into_parts();
                    assert!(core::ptr::eq(catalog, source));
                    assert_eq!(accepted_bodies.len(), 1);
                    assert_eq!(accepted_report.duration_s, DT);
                })
            },
        );
    }
}

#[test]
fn g3_late_thermal_refusal_discards_chemistry_phase_and_all_finite_owners() {
    with_liquid(0.0, 0.4, |liquid, gate| {
        with_wall(|body, contact| {
            let law = reaction(500.0, 75.0);
            let source = law.catalog();
            let gas = gas(source, 300.0, 1000.0);
            let initial_gas = gas.clone();
            let initial_bits = liquid.flow().lattice().grid().state_bits();
            let initial_h = liquid.total_enthalpy_j().unwrap();
            let binding = binding(source, 0.0);
            let surfaces = [surface(&Phase, &binding)];
            let bodies = [body];
            let initial_wall_h = bodies[0].enthalpy_j();
            let exposed = [film(INTERFACE)];
            let run = |heat: &[FreeSurfaceCellHeat]| {
                law.conducting_evaporating_liquid_candidate(
                    liquid,
                    &gas,
                    &surfaces,
                    &bodies,
                    &[],
                    heat,
                    &[contact],
                    &[],
                    &exposed,
                    joint_controls(),
                    thermal_controls(),
                    None,
                    false,
                    reactive_controls(),
                    gate,
                )
            };
            let accepted = run(&[]).unwrap();
            assert!(accepted.report().reaction.extent_mol > 0.0);
            assert!(
                accepted
                    .report()
                    .phase
                    .phase
                    .surfaces
                    .iter()
                    .any(|row| row.phase.transferred_mass_kg != 0.0)
            );
            // Rejection belongs to the full native thermal owner, after chemistry
            // and phase have both created private candidates.
            let invalid = [FreeSurfaceCellHeat {
                coordinates: INTERFACE,
                energy_j: f64::NAN,
            }];
            assert!(matches!(
                run(&invalid),
                Err(ReactiveLiquidPhaseError::Phase(_))
            ));
            assert_eq!(gas, initial_gas);
            assert_eq!(liquid.flow().lattice().grid().state_bits(), initial_bits);
            assert_eq!(liquid.total_enthalpy_j().unwrap(), initial_h);
            assert_eq!(liquid.accepted_steps(), 0);
            assert_eq!(bodies[0].enthalpy_j(), initial_wall_h);
            let repeated = run(&[]).unwrap();
            assert_eq!(accepted.gas(), repeated.gas());
            assert_eq!(
                accepted.bodies()[0].enthalpy_j(),
                repeated.bodies()[0].enthalpy_j()
            );
            assert_eq!(
                accepted.liquid().flow().lattice().grid().state_bits(),
                repeated.liquid().flow().lattice().grid().state_bits()
            );
        })
    });
}
