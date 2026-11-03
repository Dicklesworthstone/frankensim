//! Phase and sensible heat feed the original native gas-pressure reconstruction.

use super::*;
use fs_couple::liquid_gas::evaporation::{
    advance_liquid_with_pressurized_gas, advance_liquid_with_pressurized_gas_and_heat,
    multiple::{
        DistributedLiquidGasSolve, LiquidGasHeatSurface, advance_liquid_surfaces_with_gas,
        advance_liquid_surfaces_with_gas_and_heat, advance_liquid_surfaces_with_pressurized_gas,
        advance_liquid_surfaces_with_pressurized_gas_and_heat,
    },
};
use fs_lbm::material::free_surface::enthalpy::FreeSurfaceEnthalpyReport;
use fs_material::evaporation::gas::chamber::{
    GasLiquidChamberSolve, heat::GasLiquidChamberHeatSolve,
};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;

const REFERENCE: f64 = 100000.0;
const MAXIMUM: f64 = 0.05;

fn atmosphere(source: &ThermochemicalGasModel, t: f64, vapor: f64, volume: f64) -> ClosedGasParcel {
    let initial = gas(source, t, vapor);
    ClosedGasParcel::new(initial.model().clone(), initial.point().clone(), volume).unwrap()
}

fn pressure_from_inventory(gas: &ClosedGasParcel) -> f64 {
    let amount = gas.species_mass_kg(&carrier()).unwrap() / 0.028
        + gas.species_mass_kg(&water()).unwrap() / M;
    amount * UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * gas.state().state().temperature / gas.volume_m3()
}

fn normal_velocity(liquid: &MaterialFreeSurfaceEnthalpy3<'_, '_>) -> f64 {
    liquid
        .flow()
        .sample_at_si(3, 3, 3)
        .unwrap()
        .unwrap()
        .velocity_m_s[0]
}

fn inventory(
    liquid: &MaterialFreeSurfaceEnthalpy3<'_, '_>,
    gas: &ClosedGasParcel,
    binding: &GasLiquidBinding,
    next_liquid: &MaterialFreeSurfaceEnthalpy3<'_, '_>,
    next_gas: &ClosedGasParcel,
    flow: FreeSurfaceEnthalpyReport,
) {
    let dm = next_liquid.flow().ledger_mass_kg().unwrap() - liquid.flow().ledger_mass_kg().unwrap();
    let dh = next_liquid.total_enthalpy_j().unwrap() - liquid.total_enthalpy_j().unwrap();
    let du = next_gas.mass_kg() * next_gas.state().specific_internal_energy_j_kg()
        - gas.mass_kg() * gas.state().specific_internal_energy_j_kg();
    near(
        dm + next_gas.mass_kg() - gas.mass_kg(),
        0.0,
        controls().mass_tolerance_kg,
    );
    near(
        dh + binding.liquid_enthalpy_offset_j_kg * dm + du,
        flow.external_heat_j + flow.viscous_heat_j + flow.body_heat_j,
        controls().energy_tolerance_j,
    );
    assert_eq!(
        flow.external_heat_j, 0.0,
        "gas pressure is not a second caloric source"
    );
    assert_eq!(next_gas.volume_m3(), gas.volume_m3());
    assert_eq!(next_gas.point().axes().get("specimen-lot"), Some(&27.0));
}

#[test]
fn g1_evaporation_condensation_and_sensible_heat_load_the_actual_final_gas_pressure() {
    for (temperature, vapor) in [(300.0, 0.0), (320.0, 6000.0)] {
        for with_heat in [false, true] {
            with_liquid_pressure(0.0, 0.4, Some(REFERENCE), |liquid, gate| {
                let source = source(-13.4e6);
                // Initially 512 cells minus 143 full and one 0.4-filled liquid
                // cell, each 1 mm cubed. This physical volume is then held fixed.
                let original = atmosphere(&source, temperature, vapor, 3.686e-7);
                let binding = binding(&source, 0.0);
                let face = surface(&Phase, &binding);
                let next = if with_heat {
                    advance_liquid_with_pressurized_gas_and_heat(
                        liquid,
                        &original,
                        &source,
                        face,
                        100.0,
                        MAXIMUM,
                        controls(),
                        heat_controls(),
                        gate,
                    )
                } else {
                    advance_liquid_with_pressurized_gas(
                        liquid,
                        &original,
                        &source,
                        face,
                        MAXIMUM,
                        controls(),
                        gate,
                    )
                }
                .unwrap();
                let unloaded = if with_heat {
                    advance_liquid_with_gas_and_heat(
                        liquid,
                        &original,
                        &source,
                        face,
                        100.0,
                        controls(),
                        heat_controls(),
                        gate,
                    )
                } else {
                    advance_liquid_with_gas(liquid, &original, &source, face, controls(), gate)
                }
                .unwrap();
                let report = next.report();
                let pressure = report.pressure.unwrap();
                assert_eq!(next.gas(), unloaded.gas());
                assert_eq!(report.phase, unloaded.report().phase);
                assert_eq!(report.exchange, unloaded.report().exchange);
                assert_eq!(report.mass_exchange, unloaded.report().mass_exchange);
                assert!(unloaded.report().pressure.is_none());
                assert_eq!(
                    pressure.initial_pressure_pa,
                    original.state().state().pressure
                );
                assert_eq!(
                    pressure.applied_pressure_pa,
                    next.gas().state().state().pressure
                );
                assert_eq!(pressure.volume_m3, original.volume_m3());
                assert_eq!(pressure.maximum_density_excursion, MAXIMUM);
                assert!(
                    report.phase.transferred_mass_kg * (pressure.applied_pressure_pa - REFERENCE)
                        > 0.0
                );
                assert_eq!(report.heat.is_some(), with_heat);

                // Independent ideal-gas inventory and native flat-front moment
                // oracle. The five missing links have total weight 1/6; their
                // normal population/momentum increment gives -d/(1+d), with
                // d = delta_p / (rho * (dx/dt)^2). Here rho=1000 and dx/dt=1.
                // The fixture's tangential linear shear cancels in that sum.
                let expected_pressure = pressure_from_inventory(next.gas());
                near(pressure.applied_pressure_pa, expected_pressure, 2e-10);
                let d = (expected_pressure - REFERENCE) / 1000.0;
                near(normal_velocity(next.liquid()), -d / (1.0 + d), 2e-11);
                near(normal_velocity(unloaded.liquid()), 0.0, 2e-13);
                assert_ne!(
                    next.liquid().flow().lattice().grid().state_bits(),
                    unloaded.liquid().flow().lattice().grid().state_bits()
                );
                assert_eq!(next.liquid().accepted_steps(), 1);
                assert_eq!(next.liquid().physical_time_s(), DT);
                assert_eq!(report.mass_exchange.physical_time_s, 0.0);
                inventory(
                    liquid,
                    &original,
                    &binding,
                    next.liquid(),
                    next.gas(),
                    report.flow,
                );
                assert_eq!(liquid.accepted_steps(), 0);
            });
        }
    }
}

#[test]
fn g3_opposing_surfaces_and_optional_heat_share_one_post_phase_pressure_in_any_order() {
    let cells = [[3, 2, 5], [3, 4, 2]];
    for with_heat in [false, true] {
        with_liquid_interfaces_pressure(
            0.0,
            &[(cells[0], 0.4, 295.0), (cells[1], 0.4, 305.0)],
            Some(REFERENCE),
            |liquid, gate| {
                let source = source(-13.4e6);
                let original = atmosphere(&source, 300.0, 3500.0, 3.692e-7);
                let binding = binding(&source, 0.0);
                let patches = cells.map(|cell| LiquidGasEvaporationSurface {
                    liquid: cell,
                    gas: [4, cell[1], cell[2]],
                    ..surface(&Phase, &binding)
                });
                let heated = core::array::from_fn::<_, 2, _>(|i| LiquidGasHeatSurface {
                    surface: patches[i],
                    heat_transfer_coefficient_w_m2_k: [50.0, 150.0][i],
                    heat_controls: heat_controls(),
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
                let thermal = GasLiquidChamberHeatSolve {
                    maximum_iterations: 128,
                };
                let next = if with_heat {
                    advance_liquid_surfaces_with_pressurized_gas_and_heat(
                        liquid, &original, &source, &heated, MAXIMUM, joint, thermal, gate,
                    )
                } else {
                    advance_liquid_surfaces_with_pressurized_gas(
                        liquid, &original, &source, &patches, MAXIMUM, joint, gate,
                    )
                }
                .unwrap();
                let unloaded = if with_heat {
                    advance_liquid_surfaces_with_gas_and_heat(
                        liquid, &original, &source, &heated, joint, thermal, gate,
                    )
                } else {
                    advance_liquid_surfaces_with_gas(
                        liquid, &original, &source, &patches, joint, gate,
                    )
                }
                .unwrap();
                let reversed = if with_heat {
                    advance_liquid_surfaces_with_pressurized_gas_and_heat(
                        liquid,
                        &original,
                        &source,
                        &[heated[1], heated[0]],
                        MAXIMUM,
                        joint,
                        thermal,
                        gate,
                    )
                } else {
                    advance_liquid_surfaces_with_pressurized_gas(
                        liquid,
                        &original,
                        &source,
                        &[patches[1], patches[0]],
                        MAXIMUM,
                        joint,
                        gate,
                    )
                }
                .unwrap();
                let report = next.report();
                assert!(report.surfaces[0].phase.transferred_mass_kg < 0.0);
                assert!(report.surfaces[1].phase.transferred_mass_kg > 0.0);
                assert_eq!(next.gas(), unloaded.gas());
                assert_eq!(report.exchange, unloaded.report().exchange);
                assert_eq!(report.mass_exchange, unloaded.report().mass_exchange);
                assert_eq!(next.gas(), reversed.gas());
                assert_eq!(report.pressure, reversed.report().pressure);
                assert_eq!(report.flow, reversed.report().flow);
                assert_eq!(
                    next.liquid().flow().lattice().grid().state_bits(),
                    reversed.liquid().flow().lattice().grid().state_bits()
                );
                assert_eq!(next.liquid().accepted_steps(), 1);
                assert_eq!(next.liquid().physical_time_s(), DT);
                let pressure = report.pressure.unwrap();
                assert_eq!(
                    pressure.applied_pressure_pa,
                    next.gas().state().state().pressure
                );
                let d = (pressure_from_inventory(next.gas()) - REFERENCE) / 1000.0;
                near(normal_velocity(next.liquid()), -d / (1.0 + d), 2e-11);
                assert!(normal_velocity(next.liquid()).abs() > 1e-8);
                inventory(
                    liquid,
                    &original,
                    &binding,
                    next.liquid(),
                    next.gas(),
                    report.flow,
                );
                assert_eq!(report.heat.is_some(), with_heat);
            },
        );
    }
}

#[test]
fn g4_pressure_reference_budget_and_late_native_refusal_preserve_phase_inventory() {
    with_liquid_pressure(0.0, 0.4, Some(REFERENCE), |liquid, gate| {
        let source = source(-13.4e6);
        let original = atmosphere(&source, 300.0, 0.0, 3.686e-7);
        let preserved = original.clone();
        let bits = liquid.flow().lattice().grid().state_bits();
        let enthalpy = liquid.total_enthalpy_j().unwrap();
        let binding = binding(&source, 0.0);
        let face = surface(&Phase, &binding);
        for maximum in [0.0, 1.0, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                advance_liquid_with_pressurized_gas(
                    liquid,
                    &original,
                    &source,
                    face,
                    maximum,
                    controls(),
                    gate
                ),
                Err(EvaporatingLiquidGasError::Invalid(_))
            ));
        }
        let expected = advance_liquid_with_pressurized_gas(
            liquid,
            &original,
            &source,
            face,
            MAXIMUM,
            controls(),
            gate,
        )
        .unwrap();
        // This positive bound passes admission, then rejects the actual new
        // EOS traction AFTER the original gas phase and native mass/h updates.
        assert!(matches!(
            advance_liquid_with_pressurized_gas(
                liquid,
                &original,
                &source,
                face,
                1e-12,
                controls(),
                gate
            ),
            Err(EvaporatingLiquidGasError::Liquid(_))
        ));
        assert_eq!(original, preserved);
        assert_eq!(liquid.flow().lattice().grid().state_bits(), bits);
        assert_eq!(liquid.total_enthalpy_j().unwrap(), enthalpy);
        assert_eq!(liquid.physical_time_s(), 0.0);
        let retry = advance_liquid_with_pressurized_gas(
            liquid,
            &original,
            &source,
            face,
            MAXIMUM,
            controls(),
            gate,
        )
        .unwrap();
        assert_eq!(retry.gas(), expected.gas());
        assert_eq!(
            retry.liquid().flow().lattice().grid().state_bits(),
            expected.liquid().flow().lattice().grid().state_bits()
        );
        gate.request();
        assert!(matches!(
            advance_liquid_with_pressurized_gas(
                liquid,
                &original,
                &source,
                face,
                MAXIMUM,
                controls(),
                gate
            ),
            Err(EvaporatingLiquidGasError::Cancelled)
        ));
        assert_eq!(liquid.flow().lattice().grid().state_bits(), bits);
        assert_eq!(original, preserved);
    });
    with_liquid(0.0, 0.4, |liquid, gate| {
        let source = source(-13.4e6);
        let original = atmosphere(&source, 300.0, 0.0, 3.686e-7);
        let binding = binding(&source, 0.0);
        assert!(matches!(
            advance_liquid_with_pressurized_gas(
                liquid,
                &original,
                &source,
                surface(&Phase, &binding),
                MAXIMUM,
                controls(),
                gate
            ),
            Err(EvaporatingLiquidGasError::Invalid(_))
        ));
        assert_eq!(liquid.accepted_steps(), 0);
    });
}
