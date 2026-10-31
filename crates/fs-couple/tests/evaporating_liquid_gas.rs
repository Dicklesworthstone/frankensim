//! Finite original gas and phase sources drive real native mass, h and flow.

use fs_couple::liquid_gas::evaporation::{
    EvaporatingLiquidGasCandidate, EvaporatingLiquidGasError, EvaporatingLiquidGasSolve,
    LiquidGasEvaporationSurface, advance_liquid_with_gas, advance_liquid_with_gas_and_heat,
};
use fs_exec::CancelGate;
use fs_lbm::material::free_surface::enthalpy::{
    FreeSurfaceMassExchangeControls, MaterialFreeSurfaceEnthalpy3,
};
use fs_material::evaporation::gas::{
    GasEvaporationSolve, GasFilmHeatSolve, GasFilmHeatStep, GasLiquidBinding,
};
use fs_material::thermochemical_gas::{
    ThermochemicalGasModel,
    energy::{
        GasEnergySolve, parcel::ClosedGasParcel,
        parcel::transfer::species::boundary::GasSpeciesBoundarySolve,
    },
};
use fs_thermal::evaporation::{EvaporationError, EvaporationLaw, GasFilm, surface_flux};

#[path = "evaporating_liquid_gas/fixture.rs"]
mod fixture;
use fixture::*;

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, tolerance {tolerance:e}"
    );
}

fn heat_controls() -> GasFilmHeatSolve {
    GasFilmHeatSolve {
        pure_liquid_enthalpy_range_j_kg: [CP * (280.0 - 300.0), CP * (340.0 - 300.0)],
        temperature_tolerance_k: 1e-8,
        maximum_iterations: 128,
    }
}

fn check_complete_inventory(
    liquid: &MaterialFreeSurfaceEnthalpy3<'_, '_>,
    gas: &ClosedGasParcel,
    binding: &GasLiquidBinding,
    next: &EvaporatingLiquidGasCandidate<'_, '_, '_>,
) {
    let liquid_mass_change =
        next.liquid().flow().ledger_mass_kg().unwrap() - liquid.flow().ledger_mass_kg().unwrap();
    let liquid_h_change =
        next.liquid().total_enthalpy_j().unwrap() - liquid.total_enthalpy_j().unwrap();
    let gas_u_change = next.gas().mass_kg() * next.gas().state().specific_internal_energy_j_kg()
        - gas.mass_kg() * gas.state().specific_internal_energy_j_kg();
    let flow = next.report().flow;
    near(
        liquid_mass_change + next.gas().mass_kg() - gas.mass_kg(),
        0.0,
        controls().mass_tolerance_kg,
    );
    near(
        liquid_h_change + binding.liquid_enthalpy_offset_j_kg * liquid_mass_change + gas_u_change,
        flow.external_heat_j + flow.viscous_heat_j + flow.body_heat_j,
        controls().energy_tolerance_j,
    );
    near(
        liquid_h_change,
        next.report().native_liquid_enthalpy_change_j,
        1e-12,
    );
}

#[test]
fn g1_condensation_credits_actual_vof_and_upstream_vapor_heat_on_one_clock() {
    with_liquid(0.0, 0.4, |liquid, gate| {
        let source = source(-13.4e6);
        let gas = gas(&source, 320.0, 6000.0);
        let binding = binding(&source, 0.0);
        let face = surface(&Phase, &binding);
        let original_mass = liquid.flow().ledger_mass_kg().unwrap();
        let original_gas = gas.clone();
        let original_bits = liquid.flow().lattice().grid().state_bits();
        let flux = surface_flux(
            &Phase,
            300.0,
            GasFilm {
                total_pressure_pa: 100000.0,
                vapor_partial_pressure_pa: 6000.0,
                temperature_k: 320.0,
                mass_transfer_coefficient_m_s: face.mass_transfer_coefficient_m_s,
            },
        )
        .unwrap();
        let next = advance_liquid_with_gas(liquid, &gas, &source, face, controls(), gate).unwrap();
        let report = next.report();
        let dm = report.phase.transferred_mass_kg;
        assert!(dm < 0.0);
        near(
            report.phase.requested_transferred_mass_kg,
            flux.mass_flux_kg_m2_s * face.area_m2 * DT,
            1e-26,
        );
        assert_eq!(report.mass_exchange.mass_change_kg, -dm);
        assert_eq!(report.mass_exchange.physical_time_s, 0.0);
        assert_eq!(next.liquid().accepted_steps(), 1);
        assert_eq!(next.liquid().physical_time_s(), DT);
        assert!(report.phase.surface_temperature_k[1] > 300.0);
        assert!(report.phase.gas_temperature_k[1] < 320.0);
        // Incoming vapor retains its original upstream sensible enthalpy.
        let expected_heat = -dm * (LATENT + 3.5 * R * 20.0);
        near(report.mass_exchange.heat_j, expected_heat, 1e-15);
        assert!(report.mass_exchange.heat_j > -report.phase.latent_energy_j);
        near(
            next.gas().species_mass_kg(&carrier()).unwrap(),
            gas.species_mass_kg(&carrier()).unwrap(),
            1e-22,
        );
        near(
            next.liquid().flow().ledger_mass_kg().unwrap() + next.gas().mass_kg(),
            original_mass + gas.mass_kg(),
            1e-15,
        );
        near(report.mass_residual_kg, 0.0, controls().mass_tolerance_kg);
        near(report.energy_residual_j, 0.0, controls().energy_tolerance_j);
        // Independent WHOLE native inventory after advection, with the same
        // constant gas/liquid datum. The interface alone no longer holds all
        // of the deposited heat once native transport has advanced.
        let native_delta_h =
            next.liquid().total_enthalpy_j().unwrap() - liquid.total_enthalpy_j().unwrap();
        let native_delta_mass = next.liquid().flow().ledger_mass_kg().unwrap() - original_mass;
        let gas_delta_u = next.gas().mass_kg() * next.gas().state().specific_internal_energy_j_kg()
            - gas.mass_kg() * gas.state().specific_internal_energy_j_kg();
        let supplied =
            report.flow.external_heat_j + report.flow.viscous_heat_j + report.flow.body_heat_j;
        near(
            native_delta_h + binding.liquid_enthalpy_offset_j_kg * native_delta_mass + gas_delta_u,
            supplied,
            controls().energy_tolerance_j,
        );
        near(
            native_delta_h,
            report.native_liquid_enthalpy_change_j,
            1e-12,
        );
        assert_eq!(next.gas().point().axes().get("specimen-lot"), Some(&27.0));
        assert_eq!(next.gas().volume_m3(), gas.volume_m3());
        assert_ne!(next.gas().state(), gas.state());
        assert_ne!(
            next.liquid().flow().lattice().grid().state_bits(),
            original_bits
        );
        assert_eq!(gas, original_gas);
        assert_eq!(liquid.flow().lattice().grid().state_bits(), original_bits);
        assert_eq!(liquid.physical_time_s(), 0.0);
    });
}

#[test]
fn g1_dry_gas_evaporates_native_liquid_and_continuation_uses_finite_backpressure() {
    with_liquid(0.0, 0.4, |liquid, gate| {
        let source = source(-13.4e6);
        let gas = gas(&source, 300.0, 0.0);
        let binding = binding(&source, 0.0);
        let face = surface(&Phase, &binding);
        let first = advance_liquid_with_gas(liquid, &gas, &source, face, controls(), gate).unwrap();
        assert!(first.report().phase.transferred_mass_kg > 0.0);
        assert!(first.report().mass_exchange.heat_j < 0.0);
        assert!(first.report().phase.surface_temperature_k[1] < 300.0);
        assert!(first.gas().species_partial_pressure_pa(&water()).unwrap() > 0.0);
        let second = advance_liquid_with_gas(
            first.liquid(),
            first.gas(),
            first.catalog(),
            face,
            controls(),
            gate,
        )
        .unwrap();
        assert_eq!(second.liquid().accepted_steps(), 2);
        assert_eq!(second.liquid().physical_time_s(), 2.0 * DT);
        assert!(
            second.gas().species_partial_pressure_pa(&water()).unwrap()
                > first.gas().species_partial_pressure_pa(&water()).unwrap()
        );
        // Compare the SAME second-step liquid against a reset dry atmosphere.
        let reset =
            advance_liquid_with_gas(first.liquid(), &gas, &source, face, controls(), gate).unwrap();
        assert!(
            second.report().phase.transferred_mass_kg < reset.report().phase.transferred_mass_kg
        );
        near(
            second.report().energy_residual_j,
            0.0,
            controls().energy_tolerance_j,
        );
    });
}

#[test]
fn g3_common_gas_formation_shift_preserves_native_phase_temperature_and_transfer() {
    with_liquid(0.0, 0.4, |liquid, gate| {
        let run = |reference| {
            let source = source(reference);
            let gas = gas(&source, 320.0, 6000.0);
            let binding = binding(&source, 0.0);
            let next = advance_liquid_with_gas(
                liquid,
                &gas,
                &source,
                surface(&Phase, &binding),
                controls(),
                gate,
            )
            .unwrap();
            (
                next.report().phase.transferred_mass_kg,
                next.report().phase.surface_temperature_k[1],
                next.gas().state().state().temperature,
                next.report().energy_residual_j,
            )
        };
        let a = run(-13.4e6);
        let b = run(-13.4e6 + 1e9);
        assert_eq!(a.0, b.0);
        near(a.1, b.1, 1e-10);
        near(a.2, b.2, 1e-8);
        near(a.3, b.3, 1e-10);
    });
}

#[test]
fn g4_geometry_late_source_and_cancellation_preserve_both_original_owners() {
    with_liquid(0.0, 0.4, |liquid, gate| {
        let source = source(-13.4e6);
        let gas = gas(&source, 320.0, 6000.0);
        let binding = binding(&source, 0.0);
        let face = surface(&Phase, &binding);
        let original_bits = liquid.flow().lattice().grid().state_bits();
        let original_gas = gas.clone();
        let mut wall = face;
        wall.gas = [2, 3, 3];
        assert!(advance_liquid_with_gas(liquid, &gas, &source, wall, controls(), gate).is_err());
        let mut excessive_area = face;
        excessive_area.area_m2 = 2e-6;
        assert!(
            advance_liquid_with_gas(liquid, &gas, &source, excessive_area, controls(), gate)
                .is_err()
        );
        let mut wrong_datum = binding.clone();
        wrong_datum.liquid_enthalpy_offset_j_kg += 1e5;
        assert!(
            advance_liquid_with_gas(
                liquid,
                &gas,
                &source,
                surface(&Phase, &wrong_datum),
                controls(),
                gate
            )
            .is_err()
        );
        struct EndpointRefusal;
        impl EvaporationLaw for EndpointRefusal {
            fn temperature_range_k(&self) -> [f64; 2] {
                Phase.temperature_range_k()
            }
            fn properties(
                &self,
                t: f64,
            ) -> Result<fs_thermal::evaporation::EvaporationProperties, EvaporationError>
            {
                if t > 300.0 {
                    Err(EvaporationError::Invalid("endpoint source unavailable"))
                } else {
                    Phase.properties(t)
                }
            }
        }
        // Initial flux and native mass selection are supported, but the
        // condensation-heated endpoint leaves the declared source support.
        assert!(
            advance_liquid_with_gas(
                liquid,
                &gas,
                &source,
                surface(&EndpointRefusal, &binding),
                controls(),
                gate
            )
            .is_err()
        );
        gate.request();
        assert!(matches!(
            advance_liquid_with_gas(liquid, &gas, &source, face, controls(), gate),
            Err(EvaporatingLiquidGasError::Cancelled)
        ));
        assert_eq!(gas, original_gas);
        assert_eq!(liquid.flow().lattice().grid().state_bits(), original_bits);
        assert_eq!(liquid.physical_time_s(), 0.0);
    });
}

#[test]
fn g4_condensation_cannot_overfill_a_native_interface_or_publish_only_its_gas() {
    with_liquid(0.0, 1.0 - 1e-8, |liquid, gate| {
        let source = source(-13.4e6);
        let gas = gas(&source, 320.0, 6000.0);
        let binding = binding(&source, 0.0);
        let before = liquid.flow().ledger_mass_kg().unwrap();
        let original_gas = gas.clone();
        assert!(matches!(
            advance_liquid_with_gas(
                liquid,
                &gas,
                &source,
                surface(&Phase, &binding),
                controls(),
                gate
            ),
            Err(EvaporatingLiquidGasError::Liquid(_))
        ));
        assert_eq!(liquid.flow().ledger_mass_kg().unwrap(), before);
        assert_eq!(gas, original_gas);
    });
}

#[test]
fn g1_reciprocal_film_heat_cools_or_warms_before_phase_and_one_native_flow_step() {
    for (temperature, vapor_pressure) in [(290.0, 1000.0), (320.0, 6000.0)] {
        with_liquid(0.0, 0.4, |liquid, gate| {
            let source = source(-13.4e6);
            let gas = gas(&source, temperature, vapor_pressure);
            let binding = binding(&source, 0.0);
            let face = surface(&Phase, &binding);
            let legacy =
                advance_liquid_with_gas(liquid, &gas, &source, face, controls(), gate).unwrap();
            let next = advance_liquid_with_gas_and_heat(
                liquid,
                &gas,
                &source,
                face,
                50.0,
                controls(),
                heat_controls(),
                gate,
            )
            .unwrap();
            let report = next.report();
            let heat = report.heat.unwrap();
            assert!(heat.heat_into_liquid_j * (temperature - 300.0) > 0.0);
            assert!(report.phase.transferred_mass_kg != 0.0);
            assert_ne!(
                report.phase.requested_transferred_mass_kg,
                legacy.report().phase.requested_transferred_mass_kg
            );
            assert_eq!(report.exchange.surface_temperature_k[0], 300.0);
            assert_eq!(report.exchange.gas_temperature_k[0], temperature);
            assert_ne!(report.phase.surface_temperature_k[0], 300.0);
            // The reciprocal equation refers to the postheat/prephase states,
            // not the states after phase exchange or native advection.
            let conductance_time = 50.0 * face.area_m2 * DT;
            near(
                heat.heat_into_liquid_j,
                conductance_time
                    * (report.phase.gas_temperature_k[0] - report.phase.surface_temperature_k[0]),
                conductance_time * heat_controls().temperature_tolerance_k + 1e-18,
            );
            near(
                heat.liquid_enthalpy_change_j + heat.gas_internal_energy_change_j,
                0.0,
                controls().phase.gas.energy_tolerance_j,
            );
            near(
                report.exchange.gas_internal_energy_change_j,
                report.phase.gas_internal_energy_change_j + heat.gas_internal_energy_change_j,
                1e-15,
            );
            near(
                report.exchange.liquid_enthalpy_change_j,
                report.phase.liquid_enthalpy_change_j + heat.liquid_enthalpy_change_j,
                1e-15,
            );
            assert_eq!(
                report.mass_exchange.mass_change_kg,
                -report.exchange.transferred_mass_kg
            );
            assert_eq!(report.phase.duration_s, DT);
            assert_eq!(report.exchange.duration_s, DT);
            assert_eq!(next.liquid().accepted_steps(), 1);
            assert_eq!(next.liquid().physical_time_s(), DT);
            near(report.energy_residual_j, 0.0, controls().energy_tolerance_j);
            check_complete_inventory(liquid, &gas, &binding, &next);
        });
    }
}

#[test]
fn g3_zero_sensible_coefficient_retains_the_exact_legacy_native_and_gas_endpoints() {
    with_liquid(0.0, 0.4, |liquid, gate| {
        let source = source(-13.4e6);
        let gas = gas(&source, 320.0, 6000.0);
        let binding = binding(&source, 0.0);
        let face = surface(&Phase, &binding);
        let legacy =
            advance_liquid_with_gas(liquid, &gas, &source, face, controls(), gate).unwrap();
        let coupled = advance_liquid_with_gas_and_heat(
            liquid,
            &gas,
            &source,
            face,
            0.0,
            controls(),
            heat_controls(),
            gate,
        )
        .unwrap();
        assert!(legacy.report().heat.is_none());
        assert_eq!(coupled.report().heat, Some(GasFilmHeatStep::default()));
        assert_eq!(legacy.report().phase, legacy.report().exchange);
        assert_eq!(coupled.report().phase, legacy.report().phase);
        assert_eq!(coupled.report().exchange, legacy.report().exchange);
        assert_eq!(
            coupled.report().mass_exchange,
            legacy.report().mass_exchange
        );
        assert_eq!(coupled.gas(), legacy.gas());
        assert_eq!(
            coupled.liquid().flow().lattice().grid().state_bits(),
            legacy.liquid().flow().lattice().grid().state_bits()
        );
        assert_eq!(
            coupled.liquid().total_enthalpy_j().unwrap(),
            legacy.liquid().total_enthalpy_j().unwrap()
        );
        assert_eq!(coupled.liquid().physical_time_s(), DT);
    });
}

#[test]
fn g1_zero_mass_transfer_keeps_reciprocal_heat_and_matches_finite_capacity_reference() {
    with_liquid(0.0, 0.4, |liquid, gate| {
        let source = source(-13.4e6);
        let gas = gas(&source, 320.0, 6000.0);
        let binding = binding(&source, 0.0);
        let mut face = surface(&Phase, &binding);
        face.mass_transfer_coefficient_m_s = 0.0;
        let next = advance_liquid_with_gas_and_heat(
            liquid,
            &gas,
            &source,
            face,
            50.0,
            controls(),
            heat_controls(),
            gate,
        )
        .unwrap();
        let report = next.report();
        let heat = report.heat.unwrap();
        // Both authored gases have Cp/R=3.5 and the native chart has constant
        // cp. This independent two-capacity backward-Euler result is exact for
        // this fixture, without replacing the production caloric owners.
        let liquid_capacity = liquid
            .project_interface_mass_exchange(INTERFACE, 0.0)
            .unwrap()
            .initial_mass_kg
            * CP;
        let gas_capacity = gas.mass_kg() * 2.5 * gas.state().state().specific_gas_constant;
        let conductance_time = 50.0 * face.area_m2 * DT;
        let expected_heat = conductance_time * 20.0
            / (1.0 + conductance_time * (1.0 / liquid_capacity + 1.0 / gas_capacity));
        near(heat.heat_into_liquid_j, expected_heat, 2e-14);
        assert_eq!(report.phase.transferred_mass_kg, 0.0);
        assert_eq!(report.mass_exchange.mass_change_kg, 0.0);
        assert_eq!(report.phase.latent_energy_j, 0.0);
        assert_eq!(next.gas().mass_kg(), gas.mass_kg());
        near(report.mass_exchange.heat_j, heat.heat_into_liquid_j, 1e-15);
        near(
            report.exchange.gas_internal_energy_change_j,
            -expected_heat,
            2e-14,
        );
        assert_eq!(next.liquid().accepted_steps(), 1);
        assert_eq!(next.liquid().physical_time_s(), DT);
        check_complete_inventory(liquid, &gas, &binding, &next);
    });
}
