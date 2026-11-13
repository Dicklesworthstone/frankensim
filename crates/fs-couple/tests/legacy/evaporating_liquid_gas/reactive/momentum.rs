//! G1/G3 original chemistry and moving phase/thermal owners share one clock.
//! G4 temperature, late thermal and cancellation refusals publish no owner.
use super::*;
use crate::conduction::{film, joint_controls, thermal_controls, with_wall};
use fs_couple::liquid_gas::evaporation::momentum::{
    MovingGasParcel, advance_liquid_with_moving_gas,
    conduction::advance_liquid_surfaces_with_moving_gas_and_radiative_conduction,
    multiple::{
        advance_liquid_surfaces_with_moving_gas, advance_liquid_surfaces_with_moving_gas_and_heat,
    },
};
use fs_lbm::material::free_surface::enthalpy::{
    FreeSurfaceCellHeat, conduction::FreeSurfaceRadiationFace,
};

// Powers of two make the public constructor's mass/velocity reconstruction
// exact, so direct original-owner chains can compare retained momentum bits.
const DRIFT: [f64; 3] = [0.03125, -0.015625, 0.0078125];

#[test]
fn chemistry_preserves_drift_then_native_condensation_uses_that_actual_donor() {
    with_liquid(0.0, 0.4, |liquid, gate| {
        let law = reaction(2000.0, 75.0);
        let source = law.catalog();
        let gas = MovingGasParcel::new(gas(source, 300.0, 1000.0), DRIFT).unwrap();
        let original_gas = gas.clone();
        let bits = liquid.flow().lattice().grid().state_bits();
        let binding = binding(source, 0.0);
        let patch = surface(&Phase, &binding);
        let chemical = law
            .parcel_candidate(gas.parcel(), DT, reactive_controls().reaction, || false)
            .unwrap();
        assert_eq!(
            chemical.gas().mass_kg().to_bits(),
            gas.parcel().mass_kg().to_bits()
        );
        let reacted = MovingGasParcel::new(chemical.gas().clone(), gas.velocity_m_s()).unwrap();
        assert_eq!(reacted.momentum_kg_m_s(), gas.momentum_kg_m_s());
        assert_eq!(reacted.kinetic_energy_j(), gas.kinetic_energy_j());
        let direct = advance_liquid_with_moving_gas(
            liquid,
            &reacted,
            source,
            patch,
            controls(),
            None,
            None,
            gate,
        )
        .unwrap();
        let next = law
            .moving_gas_evaporating_liquid_candidate(
                liquid,
                &gas,
                patch,
                controls(),
                None,
                None,
                reactive_controls(),
                gate,
            )
            .unwrap();
        let report = next.report();
        let phase = &report.phase.native;
        let motion = &report.phase.momentum;
        assert_eq!(&report.reaction, chemical.step());
        assert!(report.reaction.extent_mol > 0.0);
        assert!(phase.phase.transferred_mass_kg < 0.0);
        assert_eq!(motion.gas_momentum_before_kg_m_s, gas.momentum_kg_m_s());
        assert_eq!(motion.donor_velocity_m_s, DRIFT);
        assert_eq!(next.gas(), direct.gas());
        assert_eq!(
            next.liquid().flow().lattice().grid().state_bits(),
            direct.liquid().flow().lattice().grid().state_bits()
        );
        assert_eq!(
            next.liquid().total_enthalpy_j().unwrap(),
            direct.liquid().total_enthalpy_j().unwrap()
        );
        let dm = next.liquid().flow().ledger_mass_kg().unwrap()
            - liquid.flow().ledger_mass_kg().unwrap();
        let dh = next.liquid().total_enthalpy_j().unwrap() - liquid.total_enthalpy_j().unwrap();
        let du = next.gas().parcel().mass_kg()
            * next.gas().parcel().state().specific_internal_energy_j_kg()
            - gas.parcel().mass_kg() * gas.parcel().state().specific_internal_energy_j_kg();
        let dk = next.gas().kinetic_energy_j() - gas.kinetic_energy_j();
        near(
            dh + binding.liquid_enthalpy_offset_j_kg * dm
                + du
                + dk
                + motion.liquid_kinetic_energy_change_j,
            phase.flow.external_heat_j + phase.flow.viscous_heat_j + phase.flow.body_heat_j,
            reactive_controls().energy_tolerance_j,
        );
        near(
            dm + next.gas().parcel().mass_kg() - gas.parcel().mass_kg(),
            0.0,
            reactive_controls().mass_tolerance_kg,
        );
        let volatile = report
            .reaction
            .species
            .iter()
            .find(|row| row.species == water())
            .unwrap();
        near(
            dm + next.gas().parcel().species_mass_kg(&water()).unwrap()
                - gas.parcel().species_mass_kg(&water()).unwrap()
                - M * volatile.stoichiometric_change_mol,
            0.0,
            reactive_controls().mass_tolerance_kg,
        );
        for d in 0..3 {
            near(
                phase.mass_exchange.momentum_change_kg_m_s[d] + next.gas().momentum_kg_m_s()[d]
                    - gas.momentum_kg_m_s()[d],
                0.0,
                controls().native.momentum_tolerance_kg_m_s,
            );
        }
        assert!(report.absolute_energy_defect_j <= reactive_controls().energy_tolerance_j);
        assert_eq!(report.duration_s, DT);
        assert_eq!(next.liquid().accepted_steps(), 1);
        assert_eq!(next.liquid().physical_time_s(), DT);
        assert_eq!(gas, original_gas);
        assert_eq!(liquid.flow().lattice().grid().state_bits(), bits);
    });
}

#[test]
fn reacting_distributed_phase_and_sensible_heat_keep_one_shared_moving_inventory() {
    let coordinates = [[3, 2, 5], [3, 4, 2]];
    with_liquid_interfaces(
        0.0,
        &[(coordinates[0], 0.4, 295.0), (coordinates[1], 0.4, 305.0)],
        |liquid, gate| {
            let law = reaction(500.0, 75.0);
            let source = law.catalog();
            let gas = MovingGasParcel::new(gas(source, 300.0, 2500.0), DRIFT).unwrap();
            let binding = binding(source, 0.0);
            let surfaces = coordinates.map(|liquid| LiquidGasEvaporationSurface {
                liquid,
                gas: film(liquid).gas,
                ..surface(&Phase, &binding)
            });
            let chemical = law
                .parcel_candidate(gas.parcel(), DT, reactive_controls().reaction, || false)
                .unwrap();
            let reacted = MovingGasParcel::new(chemical.gas().clone(), gas.velocity_m_s()).unwrap();
            let direct = advance_liquid_surfaces_with_moving_gas(
                liquid,
                &reacted,
                source,
                &surfaces,
                joint_controls(),
                None,
                gate,
            )
            .unwrap();
            let next = law
                .moving_gas_evaporating_liquid_surfaces_candidate(
                    liquid,
                    &gas,
                    &surfaces,
                    joint_controls(),
                    None,
                    reactive_controls(),
                    gate,
                )
                .unwrap();
            assert_eq!(next.gas(), direct.gas());
            assert_eq!(
                next.liquid().flow().lattice().grid().state_bits(),
                direct.liquid().flow().lattice().grid().state_bits()
            );
            assert!(next.report().phase.momentum.condensed_mass_kg > 0.0);
            assert!(next.report().phase.momentum.evaporated_mass_kg > 0.0);
            let heated = surfaces.map(|surface| LiquidGasHeatSurface {
                surface,
                heat_transfer_coefficient_w_m2_k: 1000.0,
                heat_controls: heat_controls(),
            });
            let heat = GasLiquidChamberHeatSolve {
                maximum_iterations: 128,
            };
            let direct = advance_liquid_surfaces_with_moving_gas_and_heat(
                liquid,
                &reacted,
                source,
                &heated,
                joint_controls(),
                heat,
                None,
                gate,
            )
            .unwrap();
            let next = law
                .moving_gas_evaporating_liquid_surfaces_heat_candidate(
                    liquid,
                    &gas,
                    &heated,
                    joint_controls(),
                    heat,
                    None,
                    reactive_controls(),
                    gate,
                )
                .unwrap();
            let permuted = law
                .moving_gas_evaporating_liquid_surfaces_heat_candidate(
                    liquid,
                    &gas,
                    &[heated[1], heated[0]],
                    joint_controls(),
                    heat,
                    None,
                    reactive_controls(),
                    gate,
                )
                .unwrap();
            assert_eq!(next.gas(), direct.gas());
            assert_eq!(next.gas(), permuted.gas());
            assert_eq!(
                next.liquid().flow().lattice().grid().state_bits(),
                direct.liquid().flow().lattice().grid().state_bits()
            );
            assert_eq!(
                next.liquid().flow().lattice().grid().state_bits(),
                permuted.liquid().flow().lattice().grid().state_bits()
            );
            assert_eq!(
                next.report().phase.momentum.gas_momentum_before_kg_m_s,
                gas.momentum_kg_m_s()
            );
            assert_eq!(next.liquid().accepted_steps(), 1);
            assert_eq!(next.liquid().physical_time_s(), DT);
            assert!(next.report().absolute_mass_defect_kg <= reactive_controls().mass_tolerance_kg);
            assert!(
                next.report().absolute_energy_defect_j <= reactive_controls().energy_tolerance_j
            );
        },
    );
}

#[test]
fn reacting_moving_gas_radiative_walls_and_pressure_match_original_owners_on_one_clock() {
    for maximum_density_excursion in [None, Some(0.3)] {
        with_liquid_pressure(
            0.0,
            0.4,
            maximum_density_excursion.map(|_| 100000.0),
            |liquid, gate| {
                with_wall(|body, contact| {
                    let law = reaction(500.0, 75.0);
                    let source = law.catalog();
                    let gas = MovingGasParcel::new(gas(source, 300.0, 1000.0), DRIFT).unwrap();
                    let original_gas = gas.clone();
                    let bits = liquid.flow().lattice().grid().state_bits();
                    let binding = binding(source, 0.0);
                    let surfaces = [surface(&Phase, &binding)];
                    let bodies = [body];
                    let initial_wall = bodies[0].enthalpy_j();
                    let radiation = [FreeSurfaceRadiationFace {
                        liquid: INTERFACE,
                        gas: film(INTERFACE).gas,
                        area_m2: 1e-6,
                        emissivity: 0.6.into(),
                        surroundings_k: 450.0,
                    }];
                    let heat = [FreeSurfaceCellHeat {
                        coordinates: INTERFACE,
                        energy_j: 1e-6,
                    }];
                    let chemical = law
                        .parcel_candidate(gas.parcel(), DT, reactive_controls().reaction, || false)
                        .unwrap();
                    let reacted =
                        MovingGasParcel::new(chemical.gas().clone(), gas.velocity_m_s()).unwrap();
                    let direct = advance_liquid_surfaces_with_moving_gas_and_radiative_conduction(
                        liquid,
                        &reacted,
                        source,
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
                        reactive_controls().reaction.gas,
                        maximum_density_excursion,
                        false,
                        gate,
                    )
                    .unwrap();
                    let next = law
                        .moving_gas_radiative_conducting_evaporating_liquid_candidate(
                            liquid,
                            &gas,
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
                            maximum_density_excursion,
                            false,
                            reactive_controls(),
                            gate,
                        )
                        .unwrap();
                    let report = next.report();
                    let motion = &report.phase.momentum;
                    let thermal = &report.phase.native.thermal.conduction;
                    assert_eq!(&report.reaction, chemical.step());
                    assert_eq!(motion.gas_momentum_before_kg_m_s, gas.momentum_kg_m_s());
                    assert_eq!(next.gas(), direct.gas());
                    assert_eq!(
                        next.liquid().flow().lattice().grid().state_bits(),
                        direct.liquid().flow().lattice().grid().state_bits()
                    );
                    assert_eq!(
                        next.liquid().total_enthalpy_j().unwrap(),
                        direct.liquid().total_enthalpy_j().unwrap()
                    );
                    assert_eq!(
                        next.bodies()[0].enthalpy_j(),
                        direct.bodies()[0].enthalpy_j()
                    );
                    assert_ne!(next.bodies()[0].enthalpy_j(), initial_wall);
                    assert_eq!(
                        next.gas().momentum_kg_m_s(),
                        motion.gas_momentum_after_kg_m_s
                    );
                    assert_eq!(report.duration_s, DT);
                    assert_eq!(next.liquid().accepted_steps(), 1);
                    assert_eq!(next.liquid().physical_time_s(), DT);
                    let finite = thermal.finite_gas.as_ref().unwrap();
                    match maximum_density_excursion {
                        Some(maximum) => {
                            let pressure = finite.pressure.unwrap();
                            assert_eq!(
                                pressure.applied_pressure_pa,
                                next.gas().parcel().state().state().pressure
                            );
                            assert_eq!(pressure.maximum_density_excursion, maximum);
                            assert_eq!(report.phase.native.phase.pressure, Some(pressure));
                        }
                        None => assert!(finite.pressure.is_none()),
                    }
                    let dm = next.liquid().flow().ledger_mass_kg().unwrap()
                        - liquid.flow().ledger_mass_kg().unwrap();
                    let dh = next.liquid().total_enthalpy_j().unwrap()
                        - liquid.total_enthalpy_j().unwrap();
                    let du = next.gas().parcel().mass_kg()
                        * next.gas().parcel().state().specific_internal_energy_j_kg()
                        - gas.parcel().mass_kg()
                            * gas.parcel().state().specific_internal_energy_j_kg();
                    let dw = next.bodies()[0].enthalpy_j() - initial_wall;
                    near(
                        dh + binding.liquid_enthalpy_offset_j_kg * dm
                            + du
                            + dw
                            + next.gas().kinetic_energy_j()
                            - gas.kinetic_energy_j()
                            + motion.liquid_kinetic_energy_change_j,
                        heat[0].energy_j
                            + thermal.radiation.as_ref().unwrap().heat_to_liquid_j
                            + thermal.body_radiation.as_ref().unwrap().heat_to_solids_j,
                        reactive_controls().energy_tolerance_j,
                    );
                    assert!(
                        report.absolute_energy_defect_j <= reactive_controls().energy_tolerance_j
                    );
                    assert_eq!(gas, original_gas);
                    assert_eq!(liquid.flow().lattice().grid().state_bits(), bits);
                    assert_eq!(bodies[0].enthalpy_j(), initial_wall);
                });
            },
        );
    }
}

#[test]
fn reaction_temperature_late_thermal_and_cancellation_refusals_preserve_all_moving_owners() {
    with_wall(|body, contact| {
        with_liquid(0.0, 0.4, |liquid, gate| {
            let law = reaction(500.0, 75.0);
            let gas = MovingGasParcel::new(gas(law.catalog(), 300.0, 1000.0), DRIFT).unwrap();
            let original_gas = gas.clone();
            let bits = liquid.flow().lattice().grid().state_bits();
            let original_h = liquid.total_enthalpy_j().unwrap();
            let binding = binding(law.catalog(), 0.0);
            let bodies = [body];
            let initial_wall = bodies[0].enthalpy_j();
            let mut narrow = reactive_controls();
            narrow.maximum_gas_temperature_increment_k = 1e-12;
            let mut exhausted = thermal_controls();
            exhausted.maximum_substeps = 1;
            for (reaction_controls, thermal) in [
                (narrow, thermal_controls()),
                (reactive_controls(), exhausted),
            ] {
                let result = law.moving_gas_conducting_evaporating_liquid_candidate(
                    liquid,
                    &gas,
                    &[surface(&Phase, &binding)],
                    &bodies,
                    &[],
                    &[],
                    &[contact],
                    &[],
                    &[film(INTERFACE)],
                    joint_controls(),
                    thermal,
                    None,
                    false,
                    reaction_controls,
                    gate,
                );
                assert!(result.is_err());
                assert_eq!(gas, original_gas);
                assert_eq!(liquid.flow().lattice().grid().state_bits(), bits);
                assert_eq!(liquid.total_enthalpy_j().unwrap(), original_h);
                assert_eq!(liquid.accepted_steps(), 0);
                assert_eq!(bodies[0].enthalpy_j(), initial_wall);
            }
            // Original phase-source cancellation is observed after chemistry,
            // without publishing its changed composition or the preserved P.
            struct Cancel<'a>(&'a CancelGate);
            impl EvaporationLaw for Cancel<'_> {
                fn temperature_range_k(&self) -> [f64; 2] {
                    [280.0, 340.0]
                }
                fn properties(
                    &self,
                    temperature: f64,
                ) -> Result<fs_thermal::evaporation::EvaporationProperties, EvaporationError>
                {
                    self.0.request();
                    Phase.properties(temperature)
                }
            }
            let result = law.moving_gas_evaporating_liquid_candidate(
                liquid,
                &gas,
                surface(&Cancel(gate), &binding),
                controls(),
                None,
                None,
                reactive_controls(),
                gate,
            );
            assert!(matches!(result, Err(ReactiveLiquidPhaseError::Cancelled)));
            assert!(gate.is_requested());
            assert_eq!(gas, original_gas);
            assert_eq!(liquid.flow().lattice().grid().state_bits(), bits);
            assert_eq!(liquid.total_enthalpy_j().unwrap(), original_h);
            assert_eq!(liquid.accepted_steps(), 0);
            assert_eq!(bodies[0].enthalpy_j(), initial_wall);
        });
    });
}
