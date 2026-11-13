//! Distributed finite-gas drift uses the actual original native phase owners.

use super::*;
use fs_couple::liquid_gas::evaporation::momentum::{
    MovingGasParcel,
    multiple::{
        advance_liquid_surfaces_with_moving_gas, advance_liquid_surfaces_with_moving_gas_and_heat,
    },
};
use fs_couple::liquid_gas::evaporation::multiple::{
    LiquidGasHeatSurface, advance_liquid_surfaces_with_gas_and_heat,
};
use fs_material::evaporation::gas::GasEvaporationError;
use fs_material::evaporation::gas::chamber::heat::GasLiquidChamberHeatSolve;

#[test]
fn opposing_native_phase_streams_close_finite_gas_momentum_and_total_energy() {
    with_liquid_interfaces(0.0, &INITIAL, |liquid, gate| {
        let source = source(-13.4e6);
        let gas =
            MovingGasParcel::new(gas(&source, 300.0, 3500.0), [0.025, -0.015, 0.005]).unwrap();
        let original_gas = gas.clone();
        let original_bits = liquid.flow().lattice().grid().state_bits();
        let binding = binding(&source, 0.0);
        let faces = [patch(COLD, &Phase, &binding), patch(WARM, &Phase, &binding)];
        let next = advance_liquid_surfaces_with_moving_gas(
            liquid,
            &gas,
            &source,
            &faces,
            joint_controls(),
            None,
            gate,
        )
        .unwrap();
        let report = next.report();
        let motion = &report.momentum;
        assert!(motion.evaporated_mass_kg > 0.0);
        assert!(motion.condensed_mass_kg > 0.0);
        assert!(motion.gas_mixing_heat_j > 0.0);
        assert!(motion.liquid_mixing_heat_j > 0.0);
        assert_eq!(report.native.surfaces.len(), 2);
        let mut expected = gas.momentum_kg_m_s();
        for row in &report.native.surfaces {
            let mass = row.phase.transferred_mass_kg;
            let [x, y, z] = row.coordinates;
            let donor = if mass > 0.0 {
                liquid
                    .flow()
                    .sample_at_si(x, y, z)
                    .unwrap()
                    .unwrap()
                    .velocity_m_s
            } else {
                gas.velocity_m_s()
            };
            for d in 0..3 {
                expected[d] += mass * donor[d];
            }
            assert_eq!(row.native.mass_change_kg, -mass);
            if mass > 0.0 {
                assert_eq!(row.native.mixing_heat_j, 0.0);
            } else {
                assert!(row.native.mixing_heat_j > 0.0);
            }
        }
        for d in 0..3 {
            near(next.gas().momentum_kg_m_s()[d], expected[d], 1e-22);
            near(
                report.native.mass_exchange.momentum_change_kg_m_s[d]
                    + next.gas().momentum_kg_m_s()[d]
                    - gas.momentum_kg_m_s()[d],
                0.0,
                joint_controls().native.momentum_tolerance_kg_m_s,
            );
        }
        near(
            report
                .native
                .surfaces
                .iter()
                .map(|s| s.native.mixing_heat_j)
                .sum(),
            motion.liquid_mixing_heat_j,
            1e-26,
        );
        let liquid_dm = next.liquid().flow().ledger_mass_kg().unwrap()
            - liquid.flow().ledger_mass_kg().unwrap();
        let dh = next.liquid().total_enthalpy_j().unwrap() - liquid.total_enthalpy_j().unwrap();
        let before = gas.parcel();
        let after = next.gas().parcel();
        let du = after.mass_kg() * after.state().specific_internal_energy_j_kg()
            - before.mass_kg() * before.state().specific_internal_energy_j_kg();
        let dk = next.gas().kinetic_energy_j() - gas.kinetic_energy_j();
        near(dk, motion.gas_kinetic_energy_change_j, 1e-18);
        near(
            liquid_dm + after.mass_kg() - before.mass_kg(),
            0.0,
            joint_controls().mass_tolerance_kg,
        );
        // The source kinetic account precedes native flow. Its subsequent
        // mechanical truncation is never manufactured as thermal energy.
        let flow = &report.native.flow;
        near(
            dh + binding.liquid_enthalpy_offset_j_kg * liquid_dm
                + du
                + motion.liquid_kinetic_energy_change_j
                + dk,
            flow.external_heat_j + flow.viscous_heat_j + flow.body_heat_j,
            joint_controls().energy_tolerance_j,
        );
        assert!(report.native.absolute_energy_defect_j <= joint_controls().energy_tolerance_j);
        let reversed = advance_liquid_surfaces_with_moving_gas(
            liquid,
            &gas,
            &source,
            &[faces[1], faces[0]],
            joint_controls(),
            None,
            gate,
        )
        .unwrap();
        assert_eq!(next.gas(), reversed.gas());
        assert_eq!(
            report.native.mass_exchange,
            reversed.report().native.mass_exchange
        );
        assert_eq!(
            next.liquid().flow().lattice().grid().state_bits(),
            reversed.liquid().flow().lattice().grid().state_bits()
        );
        assert_eq!(next.liquid().accepted_steps(), 1);
        assert_eq!(next.liquid().physical_time_s(), DT);
        assert_eq!(gas, original_gas);
        assert_eq!(liquid.flow().lattice().grid().state_bits(), original_bits);
    });
}

#[test]
fn thermal_only_interfaces_preserve_gas_drift_and_advance_reciprocal_heat_once() {
    with_liquid_interfaces(0.0, &INITIAL, |liquid, gate| {
        let source = source(-13.4e6);
        let gas = MovingGasParcel::new(gas(&source, 320.0, 6000.0), [0.025, 0.01, 0.0]).unwrap();
        let binding = binding(&source, 0.0);
        let faces = [COLD, WARM].map(|coordinates| {
            let mut surface = patch(coordinates, &Phase, &binding);
            surface.mass_transfer_coefficient_m_s = 0.0;
            LiquidGasHeatSurface {
                surface,
                heat_transfer_coefficient_w_m2_k: 50.0,
                heat_controls: heat_controls(),
            }
        });
        let next = advance_liquid_surfaces_with_moving_gas_and_heat(
            liquid,
            &gas,
            &source,
            &faces,
            joint_controls(),
            GasLiquidChamberHeatSolve {
                maximum_iterations: 128,
            },
            None,
            gate,
        )
        .unwrap();
        let report = next.report();
        assert_eq!(report.native.mass_exchange.samples.len(), 2);
        assert!(report.native.heat.as_ref().unwrap().heat_into_liquid_j > 0.0);
        assert_eq!(report.momentum.evaporated_mass_kg, 0.0);
        assert_eq!(report.momentum.condensed_mass_kg, 0.0);
        assert_eq!(report.momentum.gas_mixing_heat_j, 0.0);
        assert_eq!(report.momentum.liquid_mixing_heat_j, 0.0);
        assert_eq!(next.gas().momentum_kg_m_s(), gas.momentum_kg_m_s());
        assert!(
            next.gas().parcel().state().state().temperature
                < gas.parcel().state().state().temperature
        );
        assert_eq!(next.liquid().accepted_steps(), 1);
        assert_eq!(next.liquid().physical_time_s(), DT);
        assert!(report.native.absolute_energy_defect_j <= joint_controls().energy_tolerance_j);
        assert_eq!(liquid.physical_time_s(), 0.0);
    });
}

#[test]
fn heat_plus_receiver_mixing_cannot_reset_the_original_temperature_increment_budget() {
    let initial = [(COLD, 0.4, 305.0), (WARM, 0.4, 305.0)];
    with_liquid_interfaces(0.0, &initial, |liquid, gate| {
        let source = source(-13.4e6);
        let gas = MovingGasParcel::new(gas(&source, 300.0, 0.0), [300.0, -100.0, 50.0]).unwrap();
        let original = gas.clone();
        let bits = liquid.flow().lattice().grid().state_bits();
        let binding = binding(&source, 0.0);
        let faces = [COLD, WARM].map(|coordinates| LiquidGasHeatSurface {
            surface: patch(coordinates, &Phase, &binding),
            heat_transfer_coefficient_w_m2_k: 50.0,
            heat_controls: heat_controls(),
        });
        let heat = GasLiquidChamberHeatSolve {
            maximum_iterations: 128,
        };
        let mut controls = joint_controls();
        controls.chamber.surface.maximum_temperature_increment_k = 10.0;
        let baseline = advance_liquid_surfaces_with_moving_gas_and_heat(
            liquid, &gas, &source, &faces, controls, heat, None, gate,
        )
        .unwrap();
        let native = &baseline.report().native;
        assert!(native.heat.as_ref().unwrap().gas_temperature_k[1] > 300.0);
        assert!(baseline.report().momentum.gas_mixing_heat_j > 0.0);
        let final_gas = baseline.gas().parcel().state().state().temperature;
        let original_increment = final_gas - gas.parcel().state().state().temperature;
        let phase_increment = (native.exchange.gas_temperature_k[1] - 300.0).abs();
        let reset_increment = (final_gas - native.phase.gas_temperature_k[0]).abs();
        let liquid_increment = native
            .surfaces
            .iter()
            .map(|s| (s.phase.surface_temperature_k[1] - 305.0).abs())
            .fold(0.0_f64, f64::max);
        let before_mixing_limit = phase_increment.max(reset_increment).max(liquid_increment);
        assert!(original_increment > before_mixing_limit + 1e-6);
        controls.chamber.surface.maximum_temperature_increment_k =
            0.5 * (before_mixing_limit + original_increment);
        // The original heat+phase transaction fits this cap. A check measured
        // from its post-heat state would also pass; only the full accumulated
        // heat+phase+receiver-mixing increment must refuse.
        advance_liquid_surfaces_with_gas_and_heat(
            liquid,
            gas.parcel(),
            &source,
            &faces,
            controls,
            heat,
            gate,
        )
        .unwrap();
        let error = advance_liquid_surfaces_with_moving_gas_and_heat(
            liquid, &gas, &source, &faces, controls, heat, None, gate,
        )
        .err()
        .expect("combined original-to-final temperature increment must refuse");
        assert!(matches!(
            error,
            EvaporatingLiquidGasError::Phase(GasEvaporationError::StepSizeLimit {
                reason: "phase plus receiver mixing exceeds the original temperature increment bound",
                ..
            })
        ));
        assert_eq!(gas, original);
        assert_eq!(liquid.flow().lattice().grid().state_bits(), bits);
        assert_eq!(liquid.physical_time_s(), 0.0);
    });
}
