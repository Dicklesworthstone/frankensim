//! Physical turbulent entrance history on the original finite-gas heat owners.
use super::*;
use fs_airflow::graph::thermal::medium::convection::{
    evaluate_turbulent_tube_section_from_material,
    finite_gas::film::{FilmGasHeatLaw, FilmGasSurface, ForcedTurbulentTubeSectionGasSurface},
};
use fs_couple::finite_gas_thermal::WallConvectionSolve;
use fs_material::fluid_transport::FluidTransport;

fn gas(temperature_k: f64) -> GasPiston {
    GasPiston::new(
        ClosedGasParcel::new(source::model_transport(1.0, true),
            source::point(temperature_k, 100000.0), 0.001).unwrap(),
        1.0, 0.01, 0.0,
    ).unwrap()
}

fn section(fluid: &impl FluidTransport, start: f64, end: f64, area_m2: f64)
    -> ForcedTurbulentTubeSectionGasSurface
{
    let transport = fluid.transport_state().unwrap();
    ForcedTurbulentTubeSectionGasSurface {
        area_m2, diameter: Length::new(0.01), start: Length::new(start), end: Length::new(end),
        speed_m_s: 20000.0 * transport.dynamic_viscosity_pa_s()
            / (transport.density_kg_m3() * 0.01),
    }
}

#[test]
fn split_tube_keeps_bulk_sources_and_recovers_the_original_entrance_integral() {
    let actual = source::model_transport(1.0, true);
    let bulk = actual.evaluate(&source::point(400.0, 200000.0)).unwrap();
    // Extra species and unrelated initial fractions cannot replace the actual gas.
    let bank = source::model_transport(0.25, true);
    let curve = source::body_curve();
    let wall = body(&curve, 120000.0).initial;
    let resolve = |patch| FilmGasSurface {
        model: &bank, law: FilmGasHeatLaw::BulkForcedTubeSection(patch),
    }.resolve(&bulk, wall).unwrap();
    let whole_patch = section(&bulk, 0.0, 1.0, 0.02);
    let whole = resolve(whole_patch);
    let first = resolve(section(&bulk, 0.0, 0.5, 0.01));
    let second = resolve(section(&bulk, 0.5, 1.0, 0.01));
    let mean = FilmGasSurface {
        model: &bank,
        law: FilmGasHeatLaw::BulkForced(ForcedGasSurface {
            area_m2: whole_patch.area_m2, characteristic_length: whole_patch.diameter,
            speed_m_s: whole_patch.speed_m_s, correlation: CorrelationId::GnielinskiDeveloping,
            settings: ForcedConvectionSettings {
                length_over_hydraulic_diameter: Some(100.0), ..Default::default()
            },
        }),
    }.resolve(&bulk, wall).unwrap();
    assert_eq!(whole.conductance_w_k, mean.conductance_w_k);
    near((first.conductance_w_k + second.conductance_w_k) / whole.conductance_w_k, 1.0, 5e-12);
    assert!(2.0 * first.conductance_w_k > whole.conductance_w_k);

    let transport = bulk.transport_state().unwrap();
    let re = 20000.0_f64;
    let pr = transport.prandtl();
    let f = (1.8 * re.log10() - 1.5).powi(-2);
    let developed = (f / 8.0) * (re - 1000.0) * pr
        / (1.0 + 12.7 * (f / 8.0).sqrt() * (pr.powf(2.0 / 3.0) - 1.0));
    // Independently integrate Nu_mean(x) = Nu_developed [1 + (x/D)^(-2/3)].
    let expected_second = developed
        * (1.0 + (100.0_f64.cbrt() - 50.0_f64.cbrt()) / 50.0);
    near(second.source.interface.evaluation.as_ref().unwrap().effective_nusselt(), expected_second, 2e-10);
    for result in [&whole, &first, &second] {
        let sample = &result.source;
        let property = sample.property_gas.as_ref().unwrap();
        let evaluation = sample.interface.evaluation.as_ref().unwrap();
        let tube = evaluation.tube_section.as_ref().unwrap();
        assert!(evaluation.section.is_none());
        assert_eq!(sample.interface.gas, bulk);
        assert_eq!(property.state().temperature, bulk.state().temperature);
        assert_eq!(property.state().pressure, bulk.state().pressure);
        assert_eq!(property.source_packs(), bulk.source_packs());
        assert_eq!(property.thermodynamics().receipt().mole_fraction_bits(),
            bulk.thermodynamics().receipt().mole_fraction_bits());
        assert_eq!(evaluation.transport, transport);
        assert_eq!(evaluation.nusselt.card().id, CorrelationId::GnielinskiDeveloping);
        near(evaluation.nusselt.groups()["Re"], re, 1e-10);
        assert_eq!(evaluation.nusselt.groups()["L_over_Dh"], tube.end.value() / tube.diameter.value());
        near(evaluation.coefficient.value.value(), evaluation.effective_nusselt()
            * transport.thermal_conductivity_w_m_k() / tube.diameter.value(), 1e-11);
        assert_eq!(tube.upstream.is_some(), tube.start.value() > 0.0);
    }
    // K=1 retains bulk transport in heating, cooling and exact equilibrium.
    for enthalpy in [80000.0, 100000.0, 120000.0] {
        let result = FilmGasSurface { model: &bank,
            law: FilmGasHeatLaw::BulkForcedTubeSection(whole_patch) }
            .resolve(&bulk, body(&curve, enthalpy).initial).unwrap();
        assert_eq!(result.conductance_w_k, whole.conductance_w_k);
        assert_eq!(result.source.property_gas.unwrap().state().temperature, 400.0);
    }
}

#[test]
fn tube_sections_use_native_body_gas_and_piston_energy_in_both_directions() {
    let curve = source::body_curve();
    for (temperature, enthalpy) in [(400.0, 80000.0), (400.0, 120000.0), (470.0, 180000.0)] {
        let old = gas(temperature);
        let saved = old.clone();
        let thermal = body(&curve, enthalpy);
        let patch = section(old.gas().state(), 0.2, 1.0, thermal.surface_area_m2);
        let next = FilmGasSurface { model: old.gas().model(),
            law: FilmGasHeatLaw::BulkForcedTubeSection(patch) }
            .finite_body_candidate(&old, input(), thermal, controls(), exchange(), || false).unwrap();
        let sample = &next.thermal().report;
        let report = next.step();
        let evaluation = sample.interface.evaluation.as_ref().unwrap();
        let tube = evaluation.tube_section.as_ref().unwrap();
        assert_eq!((tube.start, tube.end, tube.diameter), (patch.start, patch.end, patch.diameter));
        assert_eq!(sample.interface.gas, report.piston.midpoint);
        assert_eq!(sample.interface.body, report.thermal_midpoint);
        assert_eq!(sample.property_gas.as_ref().unwrap().state().temperature,
            report.piston.midpoint.state().temperature);
        near(report.exchanged_heat_j, DT * report.conductance_w_k
            * (report.thermal_midpoint.temperature_k() - report.piston.midpoint.state().temperature),
            exchange().source_tolerance_j);
        assert_eq!(report.exchanged_heat_j > 0.0, thermal.initial.temperature_k() > temperature);
        let gas_energy = old.gas().mass_kg()
            * (next.piston().gas().state().specific_internal_energy_j_kg()
                - old.gas().state().specific_internal_energy_j_kg());
        let body_energy = thermal.mass_kg * (next.thermal().state.specific_enthalpy_j_kg()
            - thermal.initial.specific_enthalpy_j_kg());
        let kinetic = 0.5 * old.mass_kg() * next.piston().velocity_m_s().powi(2);
        let outside = input().exterior_pressure_pa * (next.piston().gas().volume_m3() - old.gas().volume_m3());
        near(gas_energy + body_energy + kinetic + outside, 0.0, 1e-7);
        near(body_energy + report.exchanged_heat_j, 0.0, 1e-7);
        if enthalpy == 180000.0 {
            assert_eq!(next.thermal().state.temperature_k(), 450.0);
            assert!(next.thermal().state.liquid_mass_fraction() > 0.5);
        }
        assert_eq!(old, saved);
    }
}

#[test]
fn fixed_volume_tube_sections_refresh_the_actual_endpoint_and_balance_both_stores() {
    scope(|cx, _| {
        let old = IsochoricGasParcel::new(source::model_transport(1.0, true),
            source::point(400.0, 100000.0), 0.003).unwrap();
        let saved = old.clone();
        let curve = source::body_curve();
        let thermal = UniformEnthalpyStepInput { duration_s: 0.01, internal_heat_j: 0.1,
            ..body(&curve, 120000.0) };
        let boundary = uniform_boundary();
        let mut environment = uniform_environment(&old, &boundary);
        environment.gas_heat_j = 0.2;
        let patch = section(old.state(), 0.2, 1.0, thermal.surface_area_m2);
        let next = FilmGasSurface { model: old.model(),
            law: FilmGasHeatLaw::BulkForcedTubeSection(patch) }
            .isochoric_body_candidate(cx, &environment, thermal, WallConvectionSolve {
                temperature_increment_bracket_k: [-15.0, 15.0],
                temperature_tolerance_k: 1e-8, maximum_iterations: 128,
            }).unwrap();
        let report = &next.report;
        let endpoint = report.next_gas.state();
        let sample = &report.convection.source.source;
        assert_eq!(&sample.interface.gas, endpoint);
        assert_eq!(sample.property_gas.as_ref().unwrap().state().temperature, endpoint.state().temperature);
        let evaluated = sample.interface.evaluation.as_ref().unwrap();
        let fresh = evaluate_turbulent_tube_section_from_material(endpoint, patch.speed_m_s,
            patch.diameter, patch.start, patch.end).unwrap();
        assert_eq!(evaluated.coefficient.value.value(), fresh.coefficient().value.value());
        assert!(evaluated.tube_section.as_ref().unwrap().upstream.is_some());
        near(report.solid_to_gas_heat_j, thermal.duration_s * thermal.surface_area_m2
            * report.convection.coefficient_w_per_m2_k
            * (next.state.temperature_k() - endpoint.state().temperature), environment.exchange.source_tolerance_j);
        let gas_energy = old.mass_kg() * (endpoint.specific_internal_energy_j_kg()
            - old.state().specific_internal_energy_j_kg());
        let body_energy = thermal.mass_kg * (next.state.specific_enthalpy_j_kg()
            - thermal.initial.specific_enthalpy_j_kg());
        near(gas_energy, report.solid_to_gas_heat_j + environment.gas_heat_j, 1e-7);
        near(body_energy, thermal.internal_heat_j - report.solid_to_gas_heat_j, 1e-7);
        near(gas_energy + body_energy, environment.gas_heat_j + thermal.internal_heat_j, 1e-7);
        assert_eq!(old, saved);
    });
}

#[test]
fn tube_section_source_and_geometry_refusals_preserve_the_original_state() {
    let old = gas(400.0);
    let saved = old.clone();
    let curve = source::body_curve();
    let thermal = body(&curve, 100000.0);
    let valid = section(old.gas().state(), 0.2, 1.0, thermal.surface_area_m2);
    // The two source cards agree numerically at 400 K but have different claims.
    let wrong = source::model_transport(1.0, false);
    assert!(FilmGasSurface { model: &wrong, law: FilmGasHeatLaw::BulkForcedTubeSection(valid) }
        .resolve(old.gas().state(), thermal.initial).is_err());
    for patch in [
        ForcedTurbulentTubeSectionGasSurface { diameter: Length::new(0.0), ..valid },
        ForcedTurbulentTubeSectionGasSurface { start: Length::new(-0.1), ..valid },
        ForcedTurbulentTubeSectionGasSurface { start: valid.end, ..valid },
        ForcedTurbulentTubeSectionGasSurface { start: Length::new(0.05), ..valid },
        ForcedTurbulentTubeSectionGasSurface { end: Length::new(f64::INFINITY), ..valid },
        ForcedTurbulentTubeSectionGasSurface { area_m2: 0.1, ..valid },
        ForcedTurbulentTubeSectionGasSurface { speed_m_s: 0.0, ..valid },
        ForcedTurbulentTubeSectionGasSurface { speed_m_s: valid.speed_m_s / 10.0, ..valid },
    ] {
        assert!(FilmGasSurface { model: old.gas().model(),
            law: FilmGasHeatLaw::BulkForcedTubeSection(patch) }
            .finite_body_candidate(&old, input(), thermal, controls(), exchange(), || false).is_err());
        assert_eq!(old, saved);
    }
    let rectangular_mean = ForcedGasSurface {
        area_m2: valid.area_m2, characteristic_length: valid.diameter, speed_m_s: valid.speed_m_s,
        correlation: CorrelationId::GnielinskiDeveloping,
        settings: ForcedConvectionSettings {
            length_over_hydraulic_diameter: Some(100.0), aspect_ratio: Some(0.5),
            ..Default::default()
        },
    };
    assert!(FilmGasSurface { model: old.gas().model(),
        law: FilmGasHeatLaw::BulkForced(rectangular_mean) }
        .resolve(old.gas().state(), thermal.initial).is_err());
    let insulated = FilmGasSurface { model: &wrong,
        law: FilmGasHeatLaw::BulkForcedTubeSection(ForcedTurbulentTubeSectionGasSurface {
            area_m2: 0.0, speed_m_s: 0.0, ..valid
        }) }
        .finite_body_candidate(&old, input(), thermal, controls(), exchange(), || false).unwrap();
    let closed = old.thermal_candidate(input(), thermal, controls(), exchange(),
        |_, _| Ok::<_, std::convert::Infallible>(ConductanceSample { conductance_w_k: 0.0, source: () }),
        || false).unwrap();
    assert_eq!(insulated.piston(), closed.piston());
    assert_eq!(insulated.thermal().state, closed.thermal().state);
    assert_eq!(insulated.step(), closed.step());
    assert!(insulated.thermal().report.property_gas.is_none());
    assert_eq!(old, saved);
}
