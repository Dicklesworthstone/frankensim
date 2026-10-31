//! Joint sensible heat reaches each native h before phase and one flow step.

use fs_couple::liquid_gas::evaporation::multiple::{
    LiquidGasHeatSurface, advance_liquid_surfaces_with_gas_and_heat,
};
use fs_material::evaporation::gas::chamber::heat::GasLiquidChamberHeatSolve;

use super::*;

fn heated_patch<'a>(
    coordinates: [usize; 3],
    coefficient: f64,
    binding: &'a GasLiquidBinding,
) -> LiquidGasHeatSurface<'a> {
    LiquidGasHeatSurface {
        surface: patch(coordinates, &Phase, binding),
        heat_transfer_coefficient_w_m2_k: coefficient,
        heat_controls: heat_controls(),
    }
}

fn joint_heat_controls() -> GasLiquidChamberHeatSolve {
    GasLiquidChamberHeatSolve {
        maximum_iterations: 128,
    }
}

#[test]
fn g1_joint_native_reciprocal_heat_precedes_phase_with_one_common_gas_and_flow_step() {
    with_liquid_interfaces(0.0, &INITIAL, |liquid, gate| {
        let source = source(-13.4e6);
        let gas = gas(&source, 300.0, 3500.0);
        let binding = binding(&source, 0.0);
        let surfaces = [
            heated_patch(COLD, 10000.0, &binding),
            heated_patch(WARM, 11000.0, &binding),
        ];
        // Strong heat moves each film by about 0.03 K and produces a resolved
        // gas temperature change. Heated-mass carry exceeds the accepted defect.
        // Omitting that carry from the native ledger must therefore refuse.
        let mut controls = joint_controls();
        controls.chamber.surface.gas.energy_tolerance_j = 1e-13;
        controls.energy_tolerance_j = 1e-12;
        let next = advance_liquid_surfaces_with_gas_and_heat(
            liquid,
            &gas,
            &source,
            &surfaces,
            controls,
            joint_heat_controls(),
            gate,
        )
        .unwrap();
        let report = next.report();
        let heat = report.heat.as_ref().unwrap();
        let cold = report.surfaces[0].heat.as_ref().unwrap();
        let warm = report.surfaces[1].heat.as_ref().unwrap();
        assert!(cold.heat_into_liquid_j > 0.0);
        assert!(warm.heat_into_liquid_j < 0.0);
        assert!(cold.temperature_k[1] > cold.temperature_k[0]);
        assert!(warm.temperature_k[1] < warm.temperature_k[0]);
        assert!(report.surfaces[0].phase.transferred_mass_kg < 0.0);
        assert!(report.surfaces[1].phase.transferred_mass_kg > 0.0);
        for row in &report.surfaces {
            let heat = row.heat.as_ref().unwrap();
            assert_eq!(row.phase.id, heat.id);
            assert_eq!(
                row.phase.surface_temperature_k[0],
                liquid
                    .liquid_temperature_at_specific_enthalpy_k(
                        heat.phase_initial_liquid.specific_enthalpy_j_kg,
                    )
                    .unwrap(),
            );
            assert_eq!(row.native.mass_change_kg, -row.phase.transferred_mass_kg);
            let heated_delta_h = heat.phase_initial_liquid.specific_enthalpy_j_kg
                - liquid.specific_enthalpy_at(row.coordinates).unwrap();
            assert!(
                row.phase.transferred_mass_kg.abs() * heated_delta_h.abs()
                    > 5.0 * controls.energy_tolerance_j
            );
        }
        assert_eq!(report.phase.gas_temperature_k[0], heat.gas_temperature_k[1]);
        assert_eq!(
            report.exchange.gas_temperature_k[0],
            gas.state().state().temperature
        );
        assert_eq!(report.mass_exchange.physical_time_s, 0.0);
        assert_eq!(next.liquid().accepted_steps(), 1);
        assert_eq!(next.liquid().physical_time_s(), DT);
        near(
            cold.liquid_enthalpy_change_j
                + warm.liquid_enthalpy_change_j
                + heat.gas_internal_energy_change_j,
            heat.energy_residual_j,
            controls.chamber.surface.gas.energy_tolerance_j,
        );
        complete_inventory(liquid, &gas, &binding, &next);

        let reversed = advance_liquid_surfaces_with_gas_and_heat(
            liquid,
            &gas,
            &source,
            &[surfaces[1], surfaces[0]],
            controls,
            joint_heat_controls(),
            gate,
        )
        .unwrap();
        assert_eq!(next.gas(), reversed.gas());
        assert_eq!(report.exchange, reversed.report().exchange);
        assert_eq!(report.mass_exchange, reversed.report().mass_exchange);
        assert_eq!(
            next.liquid().flow().lattice().grid().state_bits(),
            reversed.liquid().flow().lattice().grid().state_bits(),
        );
        assert_eq!(
            next.liquid().total_enthalpy_j().unwrap(),
            reversed.liquid().total_enthalpy_j().unwrap(),
        );
        assert_eq!(liquid.physical_time_s(), 0.0);
    });
}

#[test]
fn g3_joint_native_zero_sensible_heat_preserves_phase_only_endpoints() {
    with_liquid_interfaces(0.0, &INITIAL, |liquid, gate| {
        let source = source(-13.4e6);
        let gas = gas(&source, 300.0, 3500.0);
        let binding = binding(&source, 0.0);
        let surfaces = [
            heated_patch(COLD, 0.0, &binding),
            heated_patch(WARM, 0.0, &binding),
        ];
        let legacy = advance_liquid_surfaces_with_gas(
            liquid,
            &gas,
            &source,
            &[surfaces[0].surface, surfaces[1].surface],
            joint_controls(),
            gate,
        )
        .unwrap();
        let heated = advance_liquid_surfaces_with_gas_and_heat(
            liquid,
            &gas,
            &source,
            &surfaces,
            joint_controls(),
            joint_heat_controls(),
            gate,
        )
        .unwrap();
        assert!(legacy.report().heat.is_none());
        assert_eq!(
            heated.report().heat.as_ref().unwrap().heat_into_liquid_j,
            0.0
        );
        assert_eq!(legacy.gas(), heated.gas());
        assert_eq!(legacy.report().phase, heated.report().phase);
        assert_eq!(legacy.report().exchange, heated.report().exchange);
        assert_eq!(legacy.report().mass_exchange, heated.report().mass_exchange);
        assert_eq!(legacy.report().flow, heated.report().flow);
        assert_eq!(
            legacy.liquid().flow().lattice().grid().state_bits(),
            heated.liquid().flow().lattice().grid().state_bits(),
        );
        assert_eq!(
            legacy.liquid().total_enthalpy_j().unwrap(),
            heated.liquid().total_enthalpy_j().unwrap(),
        );
        for row in &heated.report().surfaces {
            let heat = row.heat.as_ref().unwrap();
            assert_eq!(heat.heat_into_liquid_j, 0.0);
            assert_eq!(
                heat.phase_initial_liquid.specific_enthalpy_j_kg,
                liquid.specific_enthalpy_at(row.coordinates).unwrap(),
            );
        }
        complete_inventory(liquid, &gas, &binding, &heated);
    });
}

#[test]
fn g1_joint_native_heat_only_transports_opposing_sensible_energy_without_mass() {
    with_liquid_interfaces(0.0, &INITIAL, |liquid, gate| {
        let source = source(-13.4e6);
        let gas = gas(&source, 300.0, 3500.0);
        let original_gas = gas.clone();
        let binding = binding(&source, 0.0);
        let mut surfaces = [
            heated_patch(COLD, 50.0, &binding),
            heated_patch(WARM, 150.0, &binding),
        ];
        for surface in &mut surfaces {
            surface.surface.mass_transfer_coefficient_m_s = 0.0;
        }
        let next = advance_liquid_surfaces_with_gas_and_heat(
            liquid,
            &gas,
            &source,
            &surfaces,
            joint_controls(),
            joint_heat_controls(),
            gate,
        )
        .unwrap();
        let report = next.report();
        assert_eq!(report.exchange.transferred_mass_kg, 0.0);
        assert_eq!(report.phase.latent_energy_j, 0.0);
        assert_eq!(report.mass_exchange.mass_change_kg, 0.0);
        assert!(report.surfaces[0].native.heat_j > 0.0);
        assert!(report.surfaces[1].native.heat_j < 0.0);
        assert!(next.gas().state().state().temperature > gas.state().state().temperature);
        for row in &report.surfaces {
            assert_eq!(row.native.mass_change_kg, 0.0);
            assert_eq!(row.phase.transferred_mass_kg, 0.0);
            assert_eq!(row.phase.liquid_enthalpy_change_j, 0.0);
        }
        assert_eq!(next.liquid().accepted_steps(), 1);
        assert_eq!(next.liquid().physical_time_s(), DT);
        complete_inventory(liquid, &gas, &binding, &next);
        assert_eq!(gas, original_gas);
    });
}
