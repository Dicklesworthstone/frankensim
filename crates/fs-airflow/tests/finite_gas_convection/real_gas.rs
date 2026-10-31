//! Dense homogeneous gas through the original reciprocal thermal/convection path.
//! Manufactured source packs and invariants; no empirical validation is implied.

mod natural;

use super::{body, near, scope, source, surface, uniform_boundary};
use fs_airflow::graph::thermal::medium::convection::{
    evaluate_forced_from_material, finite_gas::ForcedGasSurface,
};
use fs_conduction::lumped::LumpedThermalEnvironment;
use fs_couple::cavity_phs::piston::thermal::{HeatExchangeSolve, ThermalPistonError};
use fs_couple::cavity_phs::piston::{GasPiston, PistonInput, PistonSolve};
use fs_couple::finite_gas_thermal::{
    ConvectionCoefficient, ConvectionError, FiniteGasThermalEnvironment, WallConvectionSolve,
};
use fs_material::fluid_acoustics::FluidAcoustics;
use fs_material::fluid_transport::FluidTransport;
use fs_material::phase::UniformEnthalpyStepInput;
use fs_material::thermochemical_gas::energy::{GasEnergySolve, IsochoricGasParcel};
use fs_material::thermochemical_gas::nonideal::{
    PengRobinsonSpeciesBinding, ThermochemicalPengRobinsonModel,
    energy::parcel::IsochoricRealGasParcel,
};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::{SpeciesId, eos::peng_robinson::PengRobinsonPhaseV1};

fn gas(temperature_k: f64, energy_reference_k: f64) -> IsochoricRealGasParcel {
    let model = ThermochemicalPengRobinsonModel::new(
        source::model_transport_with_reference(0.75, true, energy_reference_k),
        vec![
            PengRobinsonSpeciesBinding {
                species: SpeciesId::new("A").unwrap(),
                critical_temperature: Temperature::new(190.564),
                critical_pressure: Pressure::new(4.5992e6),
                acentric_factor: 0.01142,
            },
            PengRobinsonSpeciesBinding {
                species: SpeciesId::new("B").unwrap(),
                critical_temperature: Temperature::new(150.687),
                critical_pressure: Pressure::new(4.863e6),
                acentric_factor: -0.00219,
            },
        ],
        vec![0.0, 0.03, 0.03, 0.0],
    )
    .unwrap();
    IsochoricRealGasParcel::new(
        model,
        source::point(temperature_k, 5.0e6),
        0.0001,
        PengRobinsonPhaseV1::Vapor,
    )
    .unwrap()
}

fn environment<'a>(
    gas: &'a IsochoricRealGasParcel,
    body: &'a LumpedThermalEnvironment,
) -> FiniteGasThermalEnvironment<'a, IsochoricRealGasParcel> {
    FiniteGasThermalEnvironment {
        gas,
        body,
        gas_heat_j: 0.0,
        gas_solve: GasEnergySolve {
            minimum_temperature_k: 280.0,
            maximum_temperature_k: 450.0,
            energy_tolerance_j_kg: 1e-7,
            maximum_iterations: 128,
        },
        exchange: HeatExchangeSolve {
            heat_bracket_j: [-60.0, 100.0],
            source_tolerance_j: 1e-9,
            energy_tolerance_j: 1e-7,
            maximum_iterations: 128,
        },
        maximum_gas_temperature_increment_k: 100.0,
    }
}

#[test]
fn dense_gas_heat_exchange_uses_real_density_calorics_and_endpoint_transport() {
    scope(|cx, _| {
        let curve = source::body_curve();
        let boundary = uniform_boundary();
        let surface = ForcedGasSurface {
            speed_m_s: 0.5,
            ..surface()
        };
        for (temperature_k, body_enthalpy) in [(300.0, 260000.0), (400.0, 50000.0)] {
            let old = gas(temperature_k, 0.0);
            let saved = old.clone();
            let environment = environment(&old, &boundary);
            let thermal = UniformEnthalpyStepInput {
                duration_s: 0.2,
                ..body(&curve, body_enthalpy)
            };
            let next = surface
                .isochoric_body_candidate(cx, &environment, thermal)
                .unwrap();
            let report = &next.report;
            let endpoint = report.next_gas.state();
            let sample = report.convection.source.as_ref().unwrap();
            let q = report.solid_to_gas_heat_j;
            assert_eq!(q > 0.0, thermal.initial.temperature_k() > temperature_k);
            assert_eq!(sample.transport, endpoint.transport_state().unwrap());
            assert_eq!(sample.transport.density_kg_m3(), old.density_kg_m3());
            assert_eq!(
                sample.transport.specific_heat_j_kg_k(),
                endpoint.specific_heat_cp_j_kg_k()
            );
            assert_ne!(
                sample.transport.density_kg_m3(),
                endpoint.source().state().density
            );
            let independently_resolved = evaluate_forced_from_material(
                surface.correlation,
                endpoint,
                surface.speed_m_s,
                surface.characteristic_length,
                surface.settings,
            )
            .unwrap();
            assert_eq!(
                sample.coefficient.value,
                independently_resolved.coefficient().value
            );
            near(
                q,
                thermal.duration_s
                    * thermal.surface_area_m2
                    * sample.coefficient.value.value()
                    * (next.state.temperature_k() - endpoint.temperature_k()),
                1e-9,
            );
            let gas_storage = old.mass_kg()
                * old
                    .model()
                    .specific_internal_energy_change_at_density(
                        old.density_kg_m3(),
                        temperature_k,
                        endpoint.temperature_k(),
                    )
                    .unwrap();
            let body_storage = thermal.mass_kg
                * (next.state.specific_enthalpy_j_kg() - thermal.initial.specific_enthalpy_j_kg());
            near(gas_storage, q, 1e-7);
            near(body_storage, -q, 1e-7);
            near(gas_storage + body_storage, 0.0, 1e-7);
            near(
                report.entropy.storage_change_j_k[1],
                old.mass_kg()
                    * old
                        .model()
                        .specific_entropy_change_between_densities(
                            old.density_kg_m3(),
                            temperature_k,
                            old.density_kg_m3(),
                            endpoint.temperature_k(),
                        )
                        .unwrap(),
                1e-12,
            );
            assert!(report.entropy.convection_production_j_k > 0.0);
            assert!(report.entropy.balance_residual_j_k.abs() < 1e-12);
            assert_eq!(report.next_gas.mass_kg(), old.mass_kg());
            assert_eq!(report.next_gas.volume_m3(), old.volume_m3());
            assert_eq!(report.next_gas.density_kg_m3(), old.density_kg_m3());
            assert!(
                (endpoint.pressure_pa() - 5e6 * endpoint.temperature_k() / temperature_k).abs()
                    > 100.0
            );
            assert_eq!(
                endpoint.acoustic_state().unwrap().density_kg_m3(),
                old.density_kg_m3()
            );
            assert_eq!(old, saved);

            // Identical source chemistry, volume and PT with the ideal EOS is
            // a different physical inventory and therefore a different heat solution.
            let ideal = IsochoricGasParcel::new(
                old.model().source().clone(),
                source::point(temperature_k, 5e6),
                old.volume_m3(),
            )
            .unwrap();
            let ideal_environment = FiniteGasThermalEnvironment {
                gas: &ideal,
                body: &boundary,
                gas_heat_j: 0.0,
                gas_solve: environment.gas_solve,
                exchange: environment.exchange,
                maximum_gas_temperature_increment_k: 100.0,
            };
            let ideal_next = surface
                .isochoric_body_candidate(cx, &ideal_environment, thermal)
                .unwrap();
            assert!((ideal_next.report.solid_to_gas_heat_j - q).abs() > 0.01);
        }
    });
}

#[test]
fn reciprocal_storage_does_not_lose_heat_beside_large_formation_references() {
    scope(|cx, _| {
        let curve = source::body_curve();
        let thermal = UniformEnthalpyStepInput {
            duration_s: 0.02,
            ..body(&curve, 260000.0)
        };
        let boundary = uniform_boundary();
        let baseline = gas(300.0, 0.0);
        let shifted = gas(300.0, 1e24);
        let surface = ForcedGasSurface {
            speed_m_s: 0.5,
            ..surface()
        };
        let a = surface
            .isochoric_body_candidate(cx, &environment(&baseline, &boundary), thermal)
            .unwrap();
        let b = surface
            .isochoric_body_candidate(cx, &environment(&shifted, &boundary), thermal)
            .unwrap();
        assert!(a.report.solid_to_gas_heat_j > 1.0);
        // Absolute subtraction really has erased the increment in this fixture.
        assert_eq!(
            b.report.next_gas.state().specific_internal_energy_j_kg()
                - shifted.state().specific_internal_energy_j_kg(),
            0.0
        );
        assert_eq!(a.state, b.state);
        assert_eq!(
            a.report.next_gas.state().temperature_k(),
            b.report.next_gas.state().temperature_k()
        );
        assert_eq!(a.report.solid_to_gas_heat_j, b.report.solid_to_gas_heat_j);
        assert_eq!(a.report.gas, b.report.gas);
        assert_eq!(a.report.entropy, b.report.entropy);
        assert!(b.report.combined_energy_residual_j.abs() < 1e-7);

        // The backward-compatible ideal specialization uses the same storage
        // audit and must retain the same small transfer for shifted references.
        let run_ideal = |source: &IsochoricRealGasParcel| {
            let gas = IsochoricGasParcel::new(
                source.model().source().clone(),
                source::point(300.0, 5e6),
                source.volume_m3(),
            )
            .unwrap();
            let controls = environment(source, &boundary);
            let environment = FiniteGasThermalEnvironment {
                gas: &gas,
                body: &boundary,
                gas_heat_j: 0.0,
                gas_solve: controls.gas_solve,
                exchange: controls.exchange,
                maximum_gas_temperature_increment_k: 100.0,
            };
            surface
                .isochoric_body_candidate(cx, &environment, thermal)
                .unwrap()
        };
        let ideal_a = run_ideal(&baseline);
        let ideal_b = run_ideal(&shifted);
        assert_eq!(ideal_a.state, ideal_b.state);
        assert_eq!(ideal_a.report.gas, ideal_b.report.gas);
        assert_eq!(ideal_a.report.entropy, ideal_b.report.entropy);
        assert!(ideal_b.report.combined_energy_residual_j.abs() < 1e-7);
    });
}

#[test]
fn real_gas_uses_the_existing_wall_feedback_and_entropy_account() {
    scope(|cx, _| {
        let curve = source::body_curve();
        let thermal = UniformEnthalpyStepInput {
            duration_s: 0.2,
            ..body(&curve, 260000.0)
        };
        let old = gas(300.0, 0.0);
        let boundary = uniform_boundary();
        let environment = environment(&old, &boundary);
        let next = environment
            .advance_with_wall_convection(
                cx,
                thermal,
                WallConvectionSolve {
                    temperature_increment_bracket_k: [-100.0, 20.0],
                    temperature_tolerance_k: 1e-7,
                    maximum_iterations: 128,
                },
                |gas, wall_k| {
                    // Manufactured wall-dependent law to exercise the same inner and
                    // outer residuals; its density is the actual nonideal endpoint.
                    Ok::<_, core::convert::Infallible>(ConvectionCoefficient {
                        coefficient_w_per_m2_k: 0.5
                            * gas.density_kg_m3()
                            * (wall_k / gas.temperature_k()).sqrt(),
                        source: (gas.density_kg_m3(), wall_k),
                    })
                },
            )
            .unwrap();
        let sample = &next.report.convection.source;
        assert_eq!(sample.source.0, old.density_kg_m3());
        assert_eq!(sample.endpoint_source.1, next.state.temperature_k());
        assert!(sample.temperature_residual_k.abs() <= 1e-7);
        assert!(sample.boundary_residual_j.abs() <= environment.exchange.source_tolerance_j);
        assert!(next.report.combined_energy_residual_j.abs() < 1e-7);
        assert!(next.report.entropy.balance_residual_j_k.abs() < 1e-12);
    });
}

#[test]
fn real_gas_source_refusal_and_callback_cancellation_preserve_both_inputs() {
    scope(|cx, gate| {
        let curve = source::body_curve();
        let thermal = UniformEnthalpyStepInput {
            duration_s: 0.2,
            ..body(&curve, 260000.0)
        };
        let old = gas(300.0, 0.0);
        let saved = old.clone();
        let boundary = uniform_boundary();
        let environment = environment(&old, &boundary);
        let unsupported = FiniteGasThermalEnvironment {
            gas_solve: GasEnergySolve {
                maximum_temperature_k: 1100.0,
                ..environment.gas_solve
            },
            ..environment.clone()
        };
        assert!(
            surface()
                .isochoric_body_candidate(cx, &unsupported, thermal)
                .is_err()
        );
        let result = environment.advance_with_convection_feedback(cx, thermal, |state, _| {
            assert_eq!(state.density_kg_m3(), old.density_kg_m3());
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
        assert_eq!(old, saved);
        assert_eq!(
            thermal.initial,
            curve.state_at_specific_enthalpy(260000.0).unwrap()
        );
    });
}

fn moving_controls() -> PistonSolve {
    PistonSolve {
        volume_bracket_m3: [0.000095, 0.000105],
        gas: GasEnergySolve {
            minimum_temperature_k: 280.0,
            maximum_temperature_k: 450.0,
            ..super::controls().gas
        },
        ..super::controls()
    }
}

fn moving_input() -> PistonInput {
    // The real 0.1-litre inventory and 2 kg wall need this shorter acoustic
    // clock; no resolution gate is relaxed to admit the dense-gas example.
    PistonInput {
        duration_s: 1e-5,
        heat_j: 0.0,
        exterior_pressure_pa: 4.5e6,
        external_force_n: 0.0,
    }
}

#[test]
fn moving_real_gas_convection_retains_actual_midpoint_and_reciprocal_energy() {
    let curve = source::body_curve();
    let surface = ForcedGasSurface {
        speed_m_s: 0.5,
        ..surface()
    };
    for (temperature_k, enthalpy) in [(300.0, 260000.0), (400.0, 50000.0)] {
        let old = GasPiston::new(gas(temperature_k, 0.0), 2.0, 0.01, 0.15).unwrap();
        let saved = old.clone();
        let thermal = UniformEnthalpyStepInput {
            duration_s: moving_input().duration_s,
            ..body(&curve, enthalpy)
        };
        let candidate = surface
            .finite_body_candidate(
                &old,
                moving_input(),
                thermal,
                moving_controls(),
                super::exchange(),
                || false,
            )
            .unwrap();
        let next = candidate.piston();
        let report = candidate.step();
        let sample = &candidate.thermal().report;
        let midpoint = &report.piston.midpoint;
        let evaluation = sample.evaluation.as_ref().unwrap();
        assert_eq!(sample.gas, *midpoint);
        assert_eq!(sample.body, report.thermal_midpoint);
        assert_eq!(evaluation.transport, midpoint.transport_state().unwrap());
        assert_eq!(
            evaluation.transport.density_kg_m3(),
            midpoint.density_kg_m3()
        );
        assert_ne!(
            evaluation.transport.density_kg_m3(),
            midpoint.source().state().density
        );
        assert_eq!(
            evaluation.transport.specific_heat_j_kg_k(),
            midpoint.specific_heat_cp_j_kg_k()
        );
        assert_eq!(
            sample.gas.source().source_packs(),
            old.gas().model().source().source_packs()
        );
        assert_eq!(sample.gas.source().source_packs().len(), 2);
        let checked = evaluate_forced_from_material(
            surface.correlation,
            midpoint,
            surface.speed_m_s,
            surface.characteristic_length,
            surface.settings,
        )
        .unwrap();
        assert_eq!(
            evaluation.coefficient.value.value(),
            checked.coefficient().value.value()
        );
        assert_eq!(
            report.conductance_w_k,
            surface.area_m2 * checked.coefficient().value.value()
        );
        let heat = report.exchanged_heat_j;
        assert_eq!(heat > 0.0, thermal.initial.temperature_k() > temperature_k);
        near(
            heat,
            moving_input().duration_s
                * report.conductance_w_k
                * (sample.body.temperature_k() - sample.gas.temperature_k()),
            super::exchange().source_tolerance_j,
        );
        assert_eq!(candidate.thermal().external_heat_j, -heat);

        let gas_energy = old.gas().mass_kg()
            * (next.gas().state().specific_internal_energy_j_kg()
                - old.gas().state().specific_internal_energy_j_kg());
        let body_energy = thermal.mass_kg
            * (candidate.thermal().state.specific_enthalpy_j_kg()
                - thermal.initial.specific_enthalpy_j_kg());
        let kinetic = 0.5
            * old.mass_kg()
            * (next.velocity_m_s() - old.velocity_m_s())
            * (next.velocity_m_s() + old.velocity_m_s());
        let dv = next.gas().volume_m3() - old.gas().volume_m3();
        near(
            gas_energy + body_energy + kinetic + moving_input().exterior_pressure_pa * dv,
            0.0,
            super::exchange().energy_tolerance_j,
        );
        near(
            body_energy + heat,
            0.0,
            super::exchange().energy_tolerance_j,
        );
        assert!(report.total_energy_residual_j.abs() <= super::exchange().energy_tolerance_j);
        assert!(next.gas().volume_m3() > old.gas().volume_m3());
        assert_eq!(next.gas().mass_kg(), old.gas().mass_kg());
        assert_eq!(
            next.gas().density_kg_m3(),
            next.gas().mass_kg() / next.gas().volume_m3()
        );
        assert_eq!(old, saved);
        assert_eq!(
            thermal.initial,
            curve.state_at_specific_enthalpy(enthalpy).unwrap()
        );
    }
}

#[test]
fn moving_real_gas_correlation_refusal_and_late_cancellation_publish_neither_state() {
    let curve = source::body_curve();
    let old = GasPiston::new(gas(300.0, 0.0), 2.0, 0.01, 0.15).unwrap();
    let saved = old.clone();
    let thermal = UniformEnthalpyStepInput {
        duration_s: moving_input().duration_s,
        ..body(&curve, 260000.0)
    };
    let initial_body = thermal.initial;
    let surface = ForcedGasSurface {
        speed_m_s: 0.5,
        ..surface()
    };
    let unsupported = ForcedGasSurface {
        speed_m_s: 1e20,
        ..surface
    };
    assert!(matches!(
        unsupported.finite_body_candidate(
            &old,
            moving_input(),
            thermal,
            moving_controls(),
            super::exchange(),
            || false
        ),
        Err(ThermalPistonError::Interface(_))
    ));
    let mut total = 0;
    surface
        .finite_body_candidate(
            &old,
            moving_input(),
            thermal,
            moving_controls(),
            super::exchange(),
            || {
                total += 1;
                false
            },
        )
        .unwrap();
    let mut polls = 0;
    assert!(matches!(
        surface.finite_body_candidate(
            &old,
            moving_input(),
            thermal,
            moving_controls(),
            super::exchange(),
            || {
                polls += 1;
                polls == total
            }
        ),
        Err(ThermalPistonError::Cancelled)
    ));
    assert_eq!(old, saved);
    assert_eq!(thermal.initial, initial_body);
}
