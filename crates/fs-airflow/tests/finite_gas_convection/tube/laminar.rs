//! G1/G3: whole-tube laminar mean on the original finite gas/body heat owners.
use super::*;
use fs_couple::finite_gas_thermal::WallConvectionSolve;
use fs_material::fluid_transport::FluidTransportState;

fn laminar_surface(old: &GasPiston) -> ForcedTubeGasSurface<'_> {
    let base = tube(old);
    let heated_length = Length::new(0.2);
    ForcedTubeGasSurface {
        area_m2: core::f64::consts::PI * base.diameter.value() * heated_length.value(),
        heated_length,
        speed_m_s: base.speed_m_s / 20.0, // Re=1000 in the original bulk state.
        ..base
    }
}

// Independent source formula, not the production correlation or tube helper.
fn reference(
    transport: FluidTransportState,
    speed: f64,
    diameter: f64,
    length: f64,
    ratio: f64,
) -> f64 {
    let re = transport.density_kg_m3() * speed * diameter / transport.dynamic_viscosity_pa_s();
    let nu = 1.86 * (re * transport.prandtl() * diameter / length).cbrt() * ratio.powf(0.14);
    nu * transport.thermal_conductivity_w_m_k() / diameter
}

#[test]
fn g1_g0_laminar_mean_keeps_original_wall_viscosity_whole_length_and_entry_admission() {
    let old = gas(400.0);
    let saved = old.clone();
    let bulk = old.gas().state();
    let curve = source::body_curve();
    let bank = source::model_transport(0.25, true);
    let surface = ForcedTubeGasSurface {
        model: &bank,
        ..laminar_surface(&old)
    };
    let transport = bulk.transport_state().unwrap();
    for enthalpy in [80000.0, 100000.0, 120000.0] {
        let wall = body(&curve, enthalpy).initial;
        assert!(
            surface.resolve(bulk, wall).is_err(),
            "selection must be explicit"
        );
        let sample = surface.laminar().resolve(bulk, wall).unwrap();
        let evaluated = sample.source.interface.evaluation.as_ref().unwrap();
        let at_wall = sample.source.surface_gas.as_ref().unwrap();
        let ratio = 400.0 / wall.temperature_k();
        let expected = reference(transport, surface.speed_m_s, 0.01, 0.2, ratio);
        assert_eq!(
            evaluated.nusselt.card().id,
            CorrelationId::SiederTateLaminar
        );
        assert_eq!(evaluated.transport, transport);
        assert_eq!(sample.source.interface.gas, *bulk);
        assert_eq!(sample.source.interface.body, wall);
        assert_eq!(at_wall.state().temperature, wall.temperature_k());
        assert_eq!(at_wall.state().pressure, bulk.state().pressure);
        assert_eq!(at_wall.source_packs(), bulk.source_packs());
        assert_eq!(
            at_wall.thermodynamics().receipt().mole_fraction_bits(),
            bulk.thermodynamics().receipt().mole_fraction_bits()
        );
        near(evaluated.nusselt.groups()["Re"], 1000.0, 1e-11);
        near(evaluated.nusselt.groups()["L_over_Dh"], 20.0, 1e-13);
        near(
            evaluated.nusselt.groups()["mu_bulk_over_surface"],
            ratio,
            1e-13,
        );
        near(
            evaluated.coefficient.value.value(),
            expected,
            2e-12 * expected,
        );
        near(
            sample.conductance_w_k,
            surface.area_m2 * expected,
            2e-12 * expected,
        );
        assert!(
            evaluated.section.is_none()
                && evaluated.tube_section.is_none()
                && evaluated.rectangular_section.is_none(),
            "a whole mean is not a local profile"
        );
    }
    let wall = body(&curve, 100000.0).initial;
    assert!(
        ForcedTubeGasSurface {
            area_m2: surface.area_m2.next_up(),
            ..surface
        }
        .laminar()
        .resolve(bulk, wall)
        .is_ok(),
        "geometry roundoff is admitted"
    );
    for invalid in [
        ForcedTubeGasSurface {
            area_m2: surface.area_m2 * 0.5,
            ..surface
        },
        ForcedTubeGasSurface {
            speed_m_s: 0.0,
            ..surface
        },
        ForcedTubeGasSurface {
            speed_m_s: surface.speed_m_s * 3.0,
            ..surface
        },
        ForcedTubeGasSurface {
            heated_length: Length::new(2.0),
            area_m2: core::f64::consts::PI * 0.01 * 2.0,
            ..surface
        },
    ] {
        assert!(invalid.laminar().resolve(bulk, wall).is_err());
    }
    let wrong_bank = source::model_transport(1.0, false);
    assert!(
        ForcedTubeGasSurface {
            model: &wrong_bank,
            ..surface
        }
        .laminar()
        .resolve(bulk, wall)
        .is_err()
    );
    let insulated = ForcedTubeGasSurface {
        area_m2: 0.0,
        speed_m_s: 0.0,
        ..surface
    }
    .laminar()
    .resolve(bulk, wall)
    .unwrap();
    assert_eq!(insulated.conductance_w_k, 0.0);
    assert!(insulated.source.surface_gas.is_none());
    assert!(insulated.source.interface.evaluation.is_none());
    assert_eq!(old, saved);
}

#[test]
fn g3_laminar_mean_closes_original_piston_body_and_latent_energy() {
    let curve = source::body_curve();
    for (temperature, enthalpy) in [(400.0, 80000.0), (400.0, 120000.0), (470.0, 180000.0)] {
        let old = gas(temperature);
        let saved = old.clone();
        let thermal = body(&curve, enthalpy);
        let surface = laminar_surface(&old);
        let next = surface
            .laminar()
            .finite_body_candidate(&old, input(), thermal, controls(), exchange(), || false)
            .unwrap();
        let sample = &next.thermal().report;
        let report = next.step();
        assert_eq!(sample.interface.gas, report.piston.midpoint);
        assert_eq!(sample.interface.body, report.thermal_midpoint);
        assert_eq!(
            sample
                .interface
                .evaluation
                .as_ref()
                .unwrap()
                .nusselt
                .card()
                .id,
            CorrelationId::SiederTateLaminar
        );
        near(
            sample.surface_gas.as_ref().unwrap().state().temperature,
            report.thermal_midpoint.temperature_k(),
            1e-12,
        );
        near(
            report.exchanged_heat_j,
            DT * report.conductance_w_k
                * (sample.interface.body.temperature_k()
                    - sample.interface.gas.state().temperature),
            exchange().source_tolerance_j,
        );
        assert_eq!(
            report.exchanged_heat_j > 0.0,
            thermal.initial.temperature_k() > temperature
        );
        let gas_heat = old.gas().mass_kg()
            * (next.piston().gas().state().specific_internal_energy_j_kg()
                - old.gas().state().specific_internal_energy_j_kg());
        let body_heat = thermal.mass_kg
            * (next.thermal().state.specific_enthalpy_j_kg()
                - thermal.initial.specific_enthalpy_j_kg());
        let kinetic = 0.5 * old.mass_kg() * next.piston().velocity_m_s().powi(2);
        let outside = input().exterior_pressure_pa
            * (next.piston().gas().volume_m3() - old.gas().volume_m3());
        near(gas_heat + body_heat + kinetic + outside, 0.0, 1e-7);
        near(body_heat + report.exchanged_heat_j, 0.0, 1e-7);
        if enthalpy == 180000.0 {
            assert_eq!(next.thermal().state.temperature_k(), 450.0);
            assert!(next.thermal().state.liquid_mass_fraction() > 0.5);
        }
        let insulated = ForcedTubeGasSurface {
            area_m2: 0.0,
            speed_m_s: 0.0,
            ..surface
        }
        .laminar()
        .finite_body_candidate(&old, input(), thermal, controls(), exchange(), || false)
        .unwrap();
        assert_eq!(insulated.step().exchanged_heat_j, 0.0);
        assert_eq!(insulated.thermal().state, thermal.initial);
        assert!(matches!(
            surface.laminar().finite_body_candidate(
                &old,
                input(),
                thermal,
                controls(),
                exchange(),
                || true,
            ),
            Err(ThermalPistonError::Cancelled)
        ));
        assert_eq!(old, saved);
    }
}

#[test]
fn g3_laminar_isochoric_wall_feedback_closes_both_stores_and_preserves_bulk_source() {
    scope(|cx, _| {
        let old = IsochoricGasParcel::new(
            source::model_transport(1.0, true),
            source::point(400.0, 100000.0),
            0.003,
        )
        .unwrap();
        let saved = old.clone();
        let curve = source::body_curve();
        let boundary = uniform_boundary();
        let environment = uniform_environment(&old, &boundary);
        let transport = old.state().transport_state().unwrap();
        let surface = ForcedTubeGasSurface {
            model: old.model(),
            area_m2: core::f64::consts::PI * 0.01 * 0.2,
            diameter: Length::new(0.01),
            heated_length: Length::new(0.2),
            speed_m_s: 1000.0 * transport.dynamic_viscosity_pa_s()
                / (transport.density_kg_m3() * 0.01),
        };
        let solve = WallConvectionSolve {
            temperature_increment_bracket_k: [-15.0, 15.0],
            temperature_tolerance_k: 1e-8,
            maximum_iterations: 128,
        };
        for enthalpy in [80000.0, 120000.0] {
            let thermal = UniformEnthalpyStepInput {
                duration_s: 0.2,
                surface_area_m2: surface.area_m2,
                ..body(&curve, enthalpy)
            };
            let next = surface
                .laminar()
                .isochoric_body_candidate(cx, &environment, thermal, solve)
                .unwrap();
            let report = &next.report;
            let wall = &report.convection.source;
            let sample = &wall.source;
            let bulk = report.next_gas.state();
            let at_wall = sample.surface_gas.as_ref().unwrap();
            let evaluated = sample.interface.evaluation.as_ref().unwrap();
            let transport = bulk.transport_state().unwrap();
            assert_eq!(&sample.interface.gas, bulk);
            assert_eq!(evaluated.transport, transport);
            assert_eq!(
                evaluated.nusselt.card().id,
                CorrelationId::SiederTateLaminar
            );
            assert_eq!(at_wall.state().temperature, wall.temperature_k);
            assert_eq!(at_wall.state().pressure, bulk.state().pressure);
            assert_eq!(at_wall.source_packs(), bulk.source_packs());
            let ratio = transport.dynamic_viscosity_pa_s()
                / at_wall.transport_state().unwrap().dynamic_viscosity_pa_s();
            let expected = reference(transport, surface.speed_m_s, 0.01, 0.2, ratio);
            near(
                report.convection.coefficient_w_per_m2_k,
                expected,
                2e-12 * expected,
            );
            near(
                wall.endpoint_source
                    .surface_gas
                    .as_ref()
                    .unwrap()
                    .state()
                    .temperature,
                next.state.temperature_k(),
                1e-12,
            );
            near(
                report.solid_to_gas_heat_j,
                thermal.duration_s
                    * thermal.surface_area_m2
                    * expected
                    * (next.state.temperature_k() - bulk.state().temperature),
                environment.exchange.source_tolerance_j,
            );
            let gas_heat = old.mass_kg()
                * (bulk.specific_internal_energy_j_kg()
                    - old.state().specific_internal_energy_j_kg());
            let body_heat = thermal.mass_kg
                * (next.state.specific_enthalpy_j_kg() - thermal.initial.specific_enthalpy_j_kg());
            near(
                gas_heat,
                report.solid_to_gas_heat_j + report.external_gas_heat_j,
                1e-7,
            );
            near(
                body_heat,
                thermal.internal_heat_j + report.radiation_heat_j - report.solid_to_gas_heat_j,
                1e-7,
            );
            near(report.combined_energy_residual_j, 0.0, 1e-7);
            assert!(report.entropy.convection_production_j_k >= 0.0);
            near(report.entropy.balance_residual_j_k, 0.0, 1e-12);
            assert!(wall.temperature_residual_k.abs() <= solve.temperature_tolerance_k);
            assert!(wall.boundary_residual_j.abs() <= environment.exchange.source_tolerance_j);
            assert_eq!(report.solid_to_gas_heat_j > 0.0, enthalpy > 100000.0);
        }
        assert_eq!(old, saved);
    });
}
