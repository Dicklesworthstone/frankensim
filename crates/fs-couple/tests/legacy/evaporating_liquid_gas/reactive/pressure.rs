//! Chemistry, phase transfer and native pressure share one accepted endpoint.
use super::*;
use fs_couple::liquid_gas::evaporation::multiple::{
    advance_liquid_surfaces_with_pressurized_gas,
    advance_liquid_surfaces_with_pressurized_gas_and_heat,
};
use fs_couple::liquid_gas::evaporation::{
    advance_liquid_with_pressurized_gas, advance_liquid_with_pressurized_gas_and_heat,
};

const PRESSURE_REFERENCE: f64 = 100000.0;
const DENSITY_ALLOWANCE: f64 = 0.3;

#[test]
fn g1_reaction_phase_and_heat_load_flow_with_the_actual_final_gas_pressure() {
    for sensible_heat in [false, true] {
        with_liquid_pressure(0.0, 0.4, Some(PRESSURE_REFERENCE), |liquid, gate| {
            let law = reaction(500.0, 75.0);
            let source = law.catalog();
            let gas = gas(source, if sensible_heat { 300.25 } else { 300.0 }, 1000.0);
            let original_gas = gas.clone();
            let bits = liquid.flow().lattice().grid().state_bits();
            let binding = binding(source, 0.0);
            let patch = surface(&Phase, &binding);
            let chemical = law
                .parcel_candidate(&gas, DT, reactive_controls().reaction, || false)
                .unwrap();
            let (direct, next, unloaded) = if sensible_heat {
                (
                    advance_liquid_with_pressurized_gas_and_heat(
                        liquid,
                        chemical.gas(),
                        source,
                        patch,
                        500.0,
                        DENSITY_ALLOWANCE,
                        controls(),
                        heat_controls(),
                        gate,
                    )
                    .unwrap(),
                    law.pressurized_evaporating_liquid_heat_candidate(
                        liquid,
                        &gas,
                        patch,
                        500.0,
                        DENSITY_ALLOWANCE,
                        controls(),
                        heat_controls(),
                        reactive_controls(),
                        gate,
                    )
                    .unwrap(),
                    law.evaporating_liquid_heat_candidate(
                        liquid,
                        &gas,
                        patch,
                        500.0,
                        controls(),
                        heat_controls(),
                        reactive_controls(),
                        gate,
                    )
                    .unwrap(),
                )
            } else {
                (
                    advance_liquid_with_pressurized_gas(
                        liquid,
                        chemical.gas(),
                        source,
                        patch,
                        DENSITY_ALLOWANCE,
                        controls(),
                        gate,
                    )
                    .unwrap(),
                    law.pressurized_evaporating_liquid_candidate(
                        liquid,
                        &gas,
                        patch,
                        DENSITY_ALLOWANCE,
                        controls(),
                        reactive_controls(),
                        gate,
                    )
                    .unwrap(),
                    law.evaporating_liquid_candidate(
                        liquid,
                        &gas,
                        patch,
                        controls(),
                        reactive_controls(),
                        gate,
                    )
                    .unwrap(),
                )
            };
            let pressure = next.report().phase.pressure.unwrap();
            assert_eq!(
                pressure.initial_pressure_pa,
                chemical.gas().state().state().pressure
            );
            assert_eq!(
                pressure.applied_pressure_pa,
                next.gas().state().state().pressure
            );
            assert_ne!(pressure.applied_pressure_pa, pressure.initial_pressure_pa);
            assert_ne!(pressure.applied_pressure_pa, gas.state().state().pressure);
            assert_eq!(pressure.volume_m3, gas.volume_m3());
            assert_eq!(pressure.maximum_density_excursion, DENSITY_ALLOWANCE);
            assert_eq!(&next.report().reaction, chemical.step());
            assert_eq!(next.gas(), direct.gas());
            assert_eq!(
                next.gas(),
                unloaded.gas(),
                "traction adds no extra gas heat or step"
            );
            assert_eq!(
                next.liquid().flow().lattice().grid().state_bits(),
                direct.liquid().flow().lattice().grid().state_bits()
            );
            assert_ne!(
                next.liquid().flow().lattice().grid().state_bits(),
                unloaded.liquid().flow().lattice().grid().state_bits()
            );
            assert_eq!(
                next.liquid().total_enthalpy_j().unwrap(),
                direct.liquid().total_enthalpy_j().unwrap()
            );
            assert!(core::ptr::eq(next.catalog(), source));
            assert_eq!(next.liquid().physical_time_s(), DT);
            let flow = next.report().phase.flow;
            complete_inventory(
                liquid,
                &gas,
                &binding,
                &next,
                flow.external_heat_j + flow.viscous_heat_j + flow.body_heat_j,
            );
            assert_eq!(gas, original_gas);
            assert_eq!(liquid.flow().lattice().grid().state_bits(), bits);
            assert_eq!(liquid.accepted_steps(), 0);
        });
    }
}

#[test]
fn g3_reactive_joint_phase_pressure_uses_one_gas_and_canonical_surface_order() {
    let coordinates = [INTERFACE, [3, 4, 3]];
    with_liquid_interfaces_pressure(
        0.0,
        &[(coordinates[0], 0.4, 298.0), (coordinates[1], 0.4, 302.0)],
        Some(PRESSURE_REFERENCE),
        |liquid, gate| {
            let law = reaction(500.0, 75.0);
            let source = law.catalog();
            let gas = gas(source, 300.0, 2500.0);
            let binding = binding(source, 0.0);
            let surfaces = coordinates.map(|liquid| LiquidGasEvaporationSurface {
                liquid,
                gas: [4, liquid[1], liquid[2]],
                mass_transfer_coefficient_m_s: 0.002,
                ..surface(&Phase, &binding)
            });
            let single = controls();
            let joint = DistributedLiquidGasSolve {
                chamber: GasLiquidChamberSolve {
                    surface: single.phase,
                    maximum_surfaces: 2,
                },
                native: single.native,
                mass_tolerance_kg: single.mass_tolerance_kg,
                energy_tolerance_j: single.energy_tolerance_j,
            };
            let chemical = law
                .parcel_candidate(&gas, DT, reactive_controls().reaction, || false)
                .unwrap();
            let direct = advance_liquid_surfaces_with_pressurized_gas(
                liquid,
                chemical.gas(),
                source,
                &surfaces,
                DENSITY_ALLOWANCE,
                joint,
                gate,
            )
            .unwrap();
            let next = law
                .pressurized_evaporating_liquid_surfaces_candidate(
                    liquid,
                    &gas,
                    &surfaces,
                    DENSITY_ALLOWANCE,
                    joint,
                    reactive_controls(),
                    gate,
                )
                .unwrap();
            assert_eq!(next.gas(), direct.gas());
            assert_eq!(
                next.liquid().flow().lattice().grid().state_bits(),
                direct.liquid().flow().lattice().grid().state_bits()
            );
            assert_eq!(
                next.report().phase.pressure.unwrap().applied_pressure_pa,
                next.gas().state().state().pressure
            );
            assert!(next.report().phase.surfaces[0].phase.transferred_mass_kg < 0.0);
            assert!(next.report().phase.surfaces[1].phase.transferred_mass_kg > 0.0);
            let flow = next.report().phase.flow;
            complete_inventory(
                liquid,
                &gas,
                &binding,
                &next,
                flow.external_heat_j + flow.viscous_heat_j + flow.body_heat_j,
            );

            let heated = surfaces.map(|surface| LiquidGasHeatSurface {
                surface,
                heat_transfer_coefficient_w_m2_k: 50.0,
                heat_controls: heat_controls(),
            });
            let heat = GasLiquidChamberHeatSolve {
                maximum_iterations: 128,
            };
            let direct = advance_liquid_surfaces_with_pressurized_gas_and_heat(
                liquid,
                chemical.gas(),
                source,
                &heated,
                DENSITY_ALLOWANCE,
                joint,
                heat,
                gate,
            )
            .unwrap();
            let next = law
                .pressurized_evaporating_liquid_surfaces_heat_candidate(
                    liquid,
                    &gas,
                    &heated,
                    DENSITY_ALLOWANCE,
                    joint,
                    heat,
                    reactive_controls(),
                    gate,
                )
                .unwrap();
            let permuted = law
                .pressurized_evaporating_liquid_surfaces_heat_candidate(
                    liquid,
                    &gas,
                    &[heated[1], heated[0]],
                    DENSITY_ALLOWANCE,
                    joint,
                    heat,
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
            let pressure = next.report().phase.pressure.unwrap();
            assert_eq!(
                pressure.initial_pressure_pa,
                chemical.gas().state().state().pressure
            );
            assert_eq!(
                pressure.applied_pressure_pa,
                next.gas().state().state().pressure
            );
            assert_eq!(
                pressure.applied_pressure_pa,
                permuted
                    .report()
                    .phase
                    .pressure
                    .unwrap()
                    .applied_pressure_pa
            );
            assert_eq!(next.liquid().physical_time_s(), DT);
            let flow = next.report().phase.flow;
            complete_inventory(
                liquid,
                &gas,
                &binding,
                &next,
                flow.external_heat_j + flow.viscous_heat_j + flow.body_heat_j,
            );
        },
    );
}

#[test]
fn g4_reactive_phase_pressure_refusals_preserve_inputs_and_exact_retry() {
    with_liquid_pressure(0.0, 0.4, Some(PRESSURE_REFERENCE), |liquid, gate| {
        let law = reaction(500.0, 75.0);
        let source = law.catalog();
        let gas = gas(source, 300.0, 1000.0);
        let original_gas = gas.clone();
        let bits = liquid.flow().lattice().grid().state_bits();
        let original_h = liquid.total_enthalpy_j().unwrap();
        let original_mass = liquid.flow().ledger_mass_kg().unwrap();
        let binding = binding(source, 0.0);
        let patch = surface(&Phase, &binding);
        let accepted = law
            .pressurized_evaporating_liquid_candidate(
                liquid,
                &gas,
                patch,
                DENSITY_ALLOWANCE,
                controls(),
                reactive_controls(),
                gate,
            )
            .unwrap();
        let refused = law.pressurized_evaporating_liquid_candidate(
            liquid,
            &gas,
            patch,
            1e-12,
            controls(),
            reactive_controls(),
            gate,
        );
        assert!(matches!(refused, Err(ReactiveLiquidPhaseError::Phase(_))));
        assert_eq!(liquid.flow().lattice().grid().state_bits(), bits);
        assert_eq!(liquid.total_enthalpy_j().unwrap(), original_h);
        assert_eq!(liquid.flow().ledger_mass_kg().unwrap(), original_mass);
        assert_eq!(liquid.accepted_steps(), 0);
        assert_eq!(liquid.physical_time_s(), 0.0);
        assert_eq!(gas, original_gas);
        let retry = law
            .pressurized_evaporating_liquid_candidate(
                liquid,
                &gas,
                patch,
                DENSITY_ALLOWANCE,
                controls(),
                reactive_controls(),
                gate,
            )
            .unwrap();
        assert_eq!(retry.gas(), accepted.gas());
        assert_eq!(
            retry.liquid().flow().lattice().grid().state_bits(),
            accepted.liquid().flow().lattice().grid().state_bits()
        );
        assert_eq!(
            retry.liquid().total_enthalpy_j().unwrap(),
            accepted.liquid().total_enthalpy_j().unwrap()
        );
    });
    with_liquid_pressure(0.0, 0.4, None, |liquid, gate| {
        let law = reaction(500.0, 75.0);
        let gas = gas(law.catalog(), 300.0, 1000.0);
        let binding = binding(law.catalog(), 0.0);
        let refused = law.pressurized_evaporating_liquid_candidate(
            liquid,
            &gas,
            surface(&Phase, &binding),
            DENSITY_ALLOWANCE,
            controls(),
            reactive_controls(),
            gate,
        );
        assert!(matches!(refused, Err(ReactiveLiquidPhaseError::Phase(_))));
        assert_eq!(liquid.accepted_steps(), 0);
    });
}
