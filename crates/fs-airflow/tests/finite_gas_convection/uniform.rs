//! G1/G3/G4: the original finite-gas/enthalpy residual consumes live transport.
use super::{near, source, surface};
use fs_airflow::graph::thermal::medium::convection::{
    evaluate_forced_from_material, finite_gas::ForcedGasSurface,
};
use fs_conduction::lumped::{LumpedThermalEnvironment, LumpedThermalTransport};
use fs_couple::cavity_phs::piston::thermal::{HeatExchangeSolve, ThermalPistonError};
use fs_couple::finite_gas_thermal::{
    ConvectionCoefficient, ConvectionError, FiniteGasThermalEnvironment,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode};
use fs_material::fluid_transport::FluidTransport;
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, UniformEnthalpyStepInput};
use fs_material::thermochemical_gas::energy::{GasEnergySolve, IsochoricGasParcel};

fn scope(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 0,
                kernel_id: 0x464743,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate);
    });
}

fn gas(fraction: f64, pressure: f64) -> IsochoricGasParcel {
    IsochoricGasParcel::new(
        source::model_transport(fraction, true),
        source::point(400.0, pressure),
        0.003,
    )
    .unwrap()
}
fn body(curve: &EquilibriumEnthalpyPhaseCurve, enthalpy: f64) -> UniformEnthalpyStepInput<'_> {
    UniformEnthalpyStepInput {
        curve,
        initial: curve.state_at_specific_enthalpy(enthalpy).unwrap(),
        mass_kg: 0.001,
        volume_m3: 1e-6,
        surface_area_m2: 0.02,
        duration_s: 0.2,
        internal_heat_j: 0.0,
    }
}
fn boundary() -> LumpedThermalEnvironment {
    LumpedThermalEnvironment {
        temperature_k: 12345.0,
        radiation_temperature_k: 300.0,
        // Both are replaced by the actual gas and material-derived film.
        convection_w_per_m2_k: f64::NAN,
        transport: LumpedThermalTransport::try_declared(100.0, 0.0).unwrap(),
        enthalpy_tolerance_j_kg: 1e-9,
        maximum_thermal_residual_j: 1e-7,
    }
}
fn environment<'a>(
    gas: &'a IsochoricGasParcel,
    body: &'a LumpedThermalEnvironment,
) -> FiniteGasThermalEnvironment<'a> {
    FiniteGasThermalEnvironment {
        gas,
        body,
        gas_heat_j: 0.0,
        gas_solve: GasEnergySolve {
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 900.0,
            energy_tolerance_j_kg: 1e-7,
            maximum_iterations: 128,
        },
        exchange: HeatExchangeSolve {
            heat_bracket_j: [-25.0, 25.0],
            source_tolerance_j: 1e-9,
            energy_tolerance_j: 1e-7,
            maximum_iterations: 128,
        },
        maximum_gas_temperature_increment_k: 100.0,
    }
}

#[test]
fn live_gas_transport_changes_the_actual_heat_solution_and_closes_storage() {
    scope(|cx, _| {
        let chart = source::body_curve();
        let thermal = body(&chart, 260000.0); // 500 K
        let gas = gas(0.75, 100000.0);
        let original = gas.clone();
        let boundary = boundary();
        let environment = environment(&gas, &boundary);
        let next = surface()
            .isochoric_body_candidate(cx, &environment, thermal)
            .unwrap();
        let report = &next.report;
        let sample = report.convection.source.as_ref().unwrap();
        let end_gas = report.next_gas.state();
        assert_eq!(sample.transport, end_gas.transport_state().unwrap());
        let independent = evaluate_forced_from_material(
            surface().correlation,
            end_gas,
            surface().speed_m_s,
            surface().characteristic_length,
            surface().settings,
        )
        .unwrap();
        assert_eq!(sample.coefficient.value, independent.coefficient().value);
        let q = report.solid_to_gas_heat_j;
        assert!(q > 0.0);
        near(
            q,
            thermal.duration_s
                * thermal.surface_area_m2
                * sample.coefficient.value.value()
                * (next.state.temperature_k() - end_gas.state().temperature),
            1e-9,
        );
        let du = gas.mass_kg()
            * (end_gas.specific_internal_energy_j_kg()
                - gas.state().specific_internal_energy_j_kg());
        let dh = thermal.mass_kg
            * (next.state.specific_enthalpy_j_kg() - thermal.initial.specific_enthalpy_j_kg());
        near(du, q, 1e-7);
        near(dh, -q, 1e-7);
        near(du + dh, 0.0, 1e-7);
        assert_eq!(gas, original);
        assert_eq!(end_gas.source_packs(), gas.state().source_packs());

        let initial = evaluate_forced_from_material(
            surface().correlation,
            gas.state(),
            surface().speed_m_s,
            surface().characteristic_length,
            surface().settings,
        )
        .unwrap();
        let frozen = LumpedThermalEnvironment {
            convection_w_per_m2_k: initial.coefficient().value.value(),
            ..boundary.clone()
        };
        let frozen_result = FiniteGasThermalEnvironment {
            body: &frozen,
            ..environment
        }
        .advance(cx, thermal)
        .unwrap();
        assert!((q - frozen_result.report.solid_to_gas_heat_j).abs() > 1e-5);
    });
}

#[test]
fn data_only_pressure_and_species_changes_reach_finite_body_temperature() {
    scope(|cx, _| {
        let chart = source::body_curve();
        let thermal = body(&chart, 260000.0);
        let boundary = boundary();
        let run = |fraction, pressure| {
            let gas = gas(fraction, pressure);
            surface()
                .isochoric_body_candidate(cx, &environment(&gas, &boundary), thermal)
                .unwrap()
        };
        let a = run(0.75, 100000.0);
        let b = run(0.75, 200000.0);
        let c = run(0.25, 100000.0);
        assert!(b.report.solid_to_gas_heat_j > a.report.solid_to_gas_heat_j);
        assert!(b.state.temperature_k() < a.state.temperature_k());
        assert!((c.state.temperature_k() - a.state.temperature_k()).abs() > 0.01);
        for next in [a, b, c] {
            assert!(next.report.combined_energy_residual_j.abs() < 1e-7);
            assert_eq!(
                next.report.convection.source.unwrap().transport,
                next.report.next_gas.state().transport_state().unwrap()
            );
        }
    });
}

#[test]
fn heat_direction_selects_the_actual_correlation_branch_inside_the_residual() {
    use fs_convection::{CorrelationId, ThermalDirection};
    scope(|cx, _| {
        let chart = source::body_curve();
        let gas = gas(0.999, 100000.0);
        let boundary = boundary();
        let environment = environment(&gas, &boundary);
        for (enthalpy, heating) in [(110000.0, true), (90000.0, false)] {
            let mut patch = surface();
            patch.correlation = CorrelationId::DittusBoelter;
            patch.speed_m_s = 10.0;
            patch.settings.length_over_hydraulic_diameter = Some(100.0);
            // Deliberately stale direction: the converged transfer owns it.
            patch.settings.direction = if heating {
                ThermalDirection::CoolingFluid
            } else {
                ThermalDirection::HeatingFluid
            };
            let next = patch
                .isochoric_body_candidate(cx, &environment, body(&chart, enthalpy))
                .unwrap();
            assert_eq!(next.report.solid_to_gas_heat_j > 0.0, heating);
            let mut actual = patch.settings;
            actual.direction = if heating {
                ThermalDirection::HeatingFluid
            } else {
                ThermalDirection::CoolingFluid
            };
            let independent = evaluate_forced_from_material(
                patch.correlation,
                next.report.next_gas.state(),
                patch.speed_m_s,
                patch.characteristic_length,
                actual,
            )
            .unwrap();
            assert_eq!(
                next.report.convection.coefficient_w_per_m2_k,
                independent.coefficient().value.value()
            );
        }
    });
}

#[test]
fn zero_heat_seed_does_not_let_a_stale_direction_falsely_fail_biot() {
    use fs_conduction::lumped::LUMPED_BIOT_CEILING;
    use fs_convection::{CorrelationId, ThermalDirection};
    scope(|cx, _| {
        let chart = source::body_curve();
        let gas = gas(0.999, 100000.0);
        let thermal = UniformEnthalpyStepInput {
            duration_s: 1e-4,
            ..body(&chart, 110000.0) // 410 K body heats the 400 K gas.
        };
        let mut patch = surface();
        patch.correlation = CorrelationId::DittusBoelter;
        patch.speed_m_s = 10.0;
        patch.settings.length_over_hydraulic_diameter = Some(100.0);
        patch.settings.direction = ThermalDirection::CoolingFluid;
        let coefficient = |direction| {
            let mut settings = patch.settings;
            settings.direction = direction;
            evaluate_forced_from_material(
                patch.correlation,
                gas.state(),
                patch.speed_m_s,
                patch.characteristic_length,
                settings,
            )
            .unwrap()
            .coefficient()
            .value
            .value()
        };
        let heating = coefficient(ThermalDirection::HeatingFluid);
        let cooling = coefficient(ThermalDirection::CoolingFluid);
        // This gas has Pr < 1. Put the two physical branch coefficients on
        // opposite sides of the original owner's Biot ceiling.
        assert!(cooling > heating);
        let length = thermal.volume_m3 / thermal.surface_area_m2;
        let conductivity = 0.5 * (heating + cooling) * length / LUMPED_BIOT_CEILING;
        assert!(heating * length / conductivity < LUMPED_BIOT_CEILING);
        assert!(cooling * length / conductivity > LUMPED_BIOT_CEILING);
        let boundary = LumpedThermalEnvironment {
            transport: LumpedThermalTransport::try_declared(conductivity, 0.0).unwrap(),
            ..boundary()
        };
        let stale = LumpedThermalEnvironment {
            temperature_k: gas.state().state().temperature,
            convection_w_per_m2_k: cooling,
            ..boundary.clone()
        };
        assert!(stale.advance(cx, thermal).is_err());
        let mut environment = environment(&gas, &boundary);
        environment.exchange.heat_bracket_j = [-0.01, 0.01];
        let next = patch
            .isochoric_body_candidate(cx, &environment, thermal)
            .unwrap();
        assert!(next.report.solid_to_gas_heat_j > 0.0);
        assert!(next.report.body.maximum_biot() < LUMPED_BIOT_CEILING);
        let mut actual = patch.settings;
        actual.direction = ThermalDirection::HeatingFluid;
        let independent = evaluate_forced_from_material(
            patch.correlation,
            next.report.next_gas.state(),
            patch.speed_m_s,
            patch.characteristic_length,
            actual,
        )
        .unwrap();
        assert_eq!(
            next.report.convection.coefficient_w_per_m2_k,
            independent.coefficient().value.value()
        );
        patch.settings.direction = ThermalDirection::HeatingFluid;
        let correctly_authored = patch
            .isochoric_body_candidate(cx, &environment, thermal)
            .unwrap();
        assert_eq!(next.state, correctly_authored.state);
        assert_eq!(
            next.report.solid_to_gas_heat_j,
            correctly_authored.report.solid_to_gas_heat_j
        );
    });
}

#[test]
fn gas_cooling_debits_latent_heat_and_insulation_does_not_query_a_correlation() {
    scope(|cx, _| {
        let chart = source::body_curve();
        let thermal = body(&chart, 180000.0); // 450 K, half melted
        let gas = gas(0.75, 100000.0);
        let boundary = boundary();
        let environment = environment(&gas, &boundary);
        let next = surface()
            .isochoric_body_candidate(cx, &environment, thermal)
            .unwrap();
        near(next.state.temperature_k(), 450.0, 1e-12);
        assert!(next.state.liquid_mass_fraction() < 0.5);
        near(
            thermal.mass_kg * 60000.0 * (0.5 - next.state.liquid_mass_fraction()),
            next.report.solid_to_gas_heat_j,
            1e-7,
        );
        let closed = ForcedGasSurface {
            area_m2: 0.0,
            speed_m_s: 0.0,
            ..surface()
        };
        let isolated = closed
            .isochoric_body_candidate(cx, &environment, thermal)
            .unwrap();
        near(
            isolated.state.specific_enthalpy_j_kg(),
            thermal.initial.specific_enthalpy_j_kg(),
            boundary.enthalpy_tolerance_j_kg,
        );
        let zero_boundary = LumpedThermalEnvironment {
            convection_w_per_m2_k: 0.0,
            ..boundary.clone()
        };
        let explicit_zero = FiniteGasThermalEnvironment {
            body: &zero_boundary,
            ..environment.clone()
        }
        .advance(cx, thermal)
        .unwrap();
        assert_eq!(isolated.state, explicit_zero.state);
        assert_eq!(isolated.report.body, explicit_zero.report.body);
        assert_eq!(isolated.report.next_gas, gas);
        assert_eq!(isolated.report.solid_to_gas_heat_j, 0.0);
        assert!(isolated.report.convection.source.is_none());
    });
}

#[test]
fn material_domain_biot_and_late_callback_cancellation_preserve_both_inputs() {
    scope(|cx, gate| {
        let chart = source::body_curve();
        let thermal = body(&chart, 260000.0);
        let gas = gas(0.75, 100000.0);
        let saved = gas.clone();
        let boundary = boundary();
        let environment = environment(&gas, &boundary);
        for surface in [
            ForcedGasSurface {
                area_m2: 0.01,
                ..surface()
            },
            ForcedGasSurface {
                speed_m_s: 1e20,
                ..surface()
            },
        ] {
            assert!(
                surface
                    .isochoric_body_candidate(cx, &environment, thermal)
                    .is_err()
            );
        }
        let insulating_body = LumpedThermalEnvironment {
            transport: LumpedThermalTransport::try_declared(0.001, 0.0).unwrap(),
            ..boundary.clone()
        };
        assert!(
            surface()
                .isochoric_body_candidate(
                    cx,
                    &FiniteGasThermalEnvironment {
                        body: &insulating_body,
                        ..environment.clone()
                    },
                    thermal
                )
                .is_err()
        );
        // Cancellation occurs after the source has actually seen a trial gas,
        // before the original body solver can publish any result.
        let result = environment.advance_with_convection_feedback(cx, thermal, |_, _| {
            gate.request();
            Ok::<_, core::convert::Infallible>(ConvectionCoefficient {
                coefficient_w_per_m2_k: 10.0,
                source: (),
            })
        });
        assert!(matches!(
            result,
            Err(ConvectionError::Thermal(ThermalPistonError::Cancelled))
        ));
        assert_eq!(gas, saved);
        assert_eq!(
            thermal.initial,
            chart.state_at_specific_enthalpy(260000.0).unwrap()
        );
    });
}
