//! Physical laminar entrances on the original finite-gas/body energy owners.
use super::*;
use fs_airflow::graph::thermal::medium::convection::{
    evaluate_laminar_tube_section_from_material, evaluate_rectangular_duct_section_from_material,
    finite_gas::film::{
        FilmGasHeatLaw, FilmGasSurface, ForcedRectangularDuctSectionGasSurface,
        ForcedTubeSectionGasSurface,
    },
};
use fs_couple::finite_gas_thermal::WallConvectionSolve;
use fs_material::fluid_transport::FluidTransport;

fn tube(
    fluid: &impl FluidTransport,
    start: f64,
    end: f64,
    area_m2: f64,
) -> ForcedTubeSectionGasSurface {
    let transport = fluid.transport_state().unwrap();
    ForcedTubeSectionGasSurface {
        area_m2,
        diameter: Length::new(0.01),
        start: Length::new(start),
        end: Length::new(end),
        speed_m_s: 1000.0 * transport.kinematic_viscosity_m2_s() / 0.01,
    }
}

fn rectangular(
    fluid: &impl FluidTransport,
    start: f64,
    end: f64,
) -> ForcedRectangularDuctSectionGasSurface {
    let transport = fluid.transport_state().unwrap();
    ForcedRectangularDuctSectionGasSurface {
        area_m2: 4.5 * 0.01 * (end - start),
        hydraulic_diameter: Length::new(0.01),
        aspect_ratio: 0.5,
        start: Length::new(start),
        end: Length::new(end),
        speed_m_s: 100.0 * transport.kinematic_viscosity_m2_s() / 0.01,
    }
}

fn pr072_pressure() -> f64 {
    fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * 400.0
}

#[test]
fn laminar_circular_patches_preserve_bulk_sources_and_the_original_entrance_integral() {
    let actual = source::model_transport(1.0, true);
    let bulk = actual.evaluate(&source::point(400.0, 100000.0)).unwrap();
    let bank = source::model_transport(0.25, true);
    let curve = source::body_curve();
    let wall = body(&curve, 120000.0).initial;
    let resolve = |patch| {
        FilmGasSurface {
            model: &bank,
            law: FilmGasHeatLaw::BulkForcedLaminarTubeSection(patch),
        }
        .resolve(&bulk, wall)
        .unwrap()
    };
    let patch = tube(&bulk, 0.0, 1.0, 0.02);
    let whole = resolve(patch);
    let first = resolve(tube(&bulk, 0.0, 0.2, 0.004));
    let second = resolve(tube(&bulk, 0.2, 1.0, 0.016));
    let original = FilmGasSurface {
        model: &bank,
        law: FilmGasHeatLaw::BulkForced(ForcedGasSurface {
            area_m2: patch.area_m2,
            characteristic_length: patch.diameter,
            speed_m_s: patch.speed_m_s,
            correlation: CorrelationId::CircularDuctHausen,
            settings: ForcedConvectionSettings {
                length_over_hydraulic_diameter: Some(100.0),
                ..Default::default()
            },
        }),
    }
    .resolve(&bulk, wall)
    .unwrap();
    assert_eq!(whole.conductance_w_k, original.conductance_w_k);
    near(
        (first.conductance_w_k + second.conductance_w_k) / whole.conductance_w_k,
        1.0,
        3e-12,
    );
    assert!(first.conductance_w_k / 0.004 > second.conductance_w_k / 0.016);
    let reset = resolve(tube(&bulk, 0.0, 0.8, 0.016));
    assert!(first.conductance_w_k + reset.conductance_w_k > 1.05 * whole.conductance_w_k);
    let pr = bulk.transport_state().unwrap().prandtl();
    let mean = |length: f64| {
        let gz = 1000.0 * pr / (length / 0.01);
        3.66 + 0.0668 * gz / (1.0 + 0.04 * gz.powf(2.0 / 3.0))
    };
    near(
        second
            .source
            .interface
            .evaluation
            .as_ref()
            .unwrap()
            .effective_nusselt(),
        (mean(1.0) - 0.2 * mean(0.2)) / 0.8,
        3e-12,
    );
    for result in [&whole, &first, &second] {
        let sample = &result.source;
        let property = sample.property_gas.as_ref().unwrap();
        let evaluation = sample.interface.evaluation.as_ref().unwrap();
        let section = evaluation.tube_section.as_ref().unwrap();
        assert!(evaluation.section.is_none() && evaluation.rectangular_section.is_none());
        assert_eq!(sample.interface.gas, bulk);
        assert_eq!(property.state().temperature, bulk.state().temperature);
        assert_eq!(property.state().pressure, bulk.state().pressure);
        assert_eq!(property.source_packs(), bulk.source_packs());
        assert_eq!(
            property.thermodynamics().receipt().mole_fraction_bits(),
            bulk.thermodynamics().receipt().mole_fraction_bits()
        );
        assert_eq!(evaluation.transport, bulk.transport_state().unwrap());
        assert_eq!(
            evaluation.nusselt.card().id,
            CorrelationId::CircularDuctHausen
        );
        assert_eq!(section.upstream.is_some(), section.start.value() > 0.0);
        assert_eq!(
            evaluation.coefficient.model.discrepancy_rel.is_infinite(),
            section.start.value() > 0.0
        );
        near(
            evaluation.coefficient.qoi,
            evaluation.effective_nusselt() * evaluation.transport.thermal_conductivity_w_m_k()
                / section.diameter.value(),
            1e-12,
        );
    }
}

#[test]
fn rectangular_patches_preserve_four_wall_heat_across_the_table_and_declared_bridge() {
    let model = source::model_pr072();
    let bulk = model
        .evaluate(&source::point(400.0, pr072_pressure()))
        .unwrap();
    assert_eq!(bulk.transport_state().unwrap().prandtl(), 0.72);
    let curve = source::body_curve();
    let wall = body(&curve, 120000.0).initial;
    let resolve = |patch| {
        FilmGasSurface {
            model: &model,
            law: FilmGasHeatLaw::BulkForcedRectangularDuctSection(patch),
        }
        .resolve(&bulk, wall)
        .unwrap()
    };
    let whole = resolve(rectangular(&bulk, 0.0, 0.12));
    let first = resolve(rectangular(&bulk, 0.0, 0.04));
    let second = resolve(rectangular(&bulk, 0.04, 0.12));
    near(
        (first.conductance_w_k + second.conductance_w_k) / whole.conductance_w_k,
        1.0,
        3e-12,
    );
    // Independent table values: Nu(10)=4.20, Nu(20)=4.79; CWT limit at a=0.5.
    let developed = 3.388736875;
    let first_mean = 4.20 + 0.8 * (4.79 - 4.20); // Gz=18.
    let whole_mean = developed + 0.6 * (4.20 - developed); // Gz=6 bridge.
    let expected_section = (0.12 * whole_mean - 0.04 * first_mean) / 0.08;
    let evaluation = second.source.interface.evaluation.as_ref().unwrap();
    near(evaluation.effective_nusselt(), expected_section, 3e-12);
    let section = evaluation.rectangular_section.as_ref().unwrap();
    assert_eq!(section.aspect_ratio, 0.5);
    assert_eq!((section.start.value(), section.end.value()), (0.04, 0.12));
    assert_eq!(section.hydraulic_diameter.value(), 0.01);
    assert_eq!(
        evaluation.nusselt.card().id,
        CorrelationId::RectangularDuctLaminarCwtDevelopingPr072
    );
    assert_eq!(
        section.upstream.as_ref().unwrap().source_region(),
        fs_convection::EvaluationSourceRegion::TableSliceInterpolation
    );
    assert_eq!(
        evaluation.nusselt.source_region(),
        fs_convection::EvaluationSourceRegion::DeclaredEngineeringBridge
    );
    assert!(evaluation.tube_section.is_none() && evaluation.section.is_none());
    assert!(evaluation.coefficient.model.discrepancy_rel.is_infinite());
    assert_eq!(
        second
            .source
            .property_gas
            .as_ref()
            .unwrap()
            .state()
            .temperature,
        400.0
    );
    assert_eq!(
        second.source.property_gas.as_ref().unwrap().source_packs(),
        bulk.source_packs()
    );
}

#[test]
fn laminar_tube_heat_uses_native_body_gas_and_piston_energy_in_both_directions() {
    let curve = source::body_curve();
    for (temperature, enthalpy) in [(400.0, 80000.0), (400.0, 120000.0), (470.0, 180000.0)] {
        let old = GasPiston::new(
            ClosedGasParcel::new(
                source::model_transport(1.0, true),
                source::point(temperature, 100000.0),
                0.001,
            )
            .unwrap(),
            1.0,
            0.01,
            0.0,
        )
        .unwrap();
        let saved = old.clone();
        let thermal = body(&curve, enthalpy);
        let patch = tube(old.gas().state(), 0.2, 1.0, thermal.surface_area_m2);
        let next = FilmGasSurface {
            model: old.gas().model(),
            law: FilmGasHeatLaw::BulkForcedLaminarTubeSection(patch),
        }
        .finite_body_candidate(&old, input(), thermal, controls(), exchange(), || false)
        .unwrap();
        let report = next.step();
        let sample = &next.thermal().report;
        assert_eq!(sample.interface.gas, report.piston.midpoint);
        assert_eq!(sample.interface.body, report.thermal_midpoint);
        assert_eq!(
            sample.property_gas.as_ref().unwrap().state().temperature,
            report.piston.midpoint.state().temperature
        );
        near(
            report.exchanged_heat_j,
            DT * report.conductance_w_k
                * (report.thermal_midpoint.temperature_k()
                    - report.piston.midpoint.state().temperature),
            exchange().source_tolerance_j,
        );
        assert_eq!(
            report.exchanged_heat_j > 0.0,
            thermal.initial.temperature_k() > temperature
        );
        let gas_energy = old.gas().mass_kg()
            * (next.piston().gas().state().specific_internal_energy_j_kg()
                - old.gas().state().specific_internal_energy_j_kg());
        let body_energy = thermal.mass_kg
            * (next.thermal().state.specific_enthalpy_j_kg()
                - thermal.initial.specific_enthalpy_j_kg());
        let kinetic = 0.5 * old.mass_kg() * next.piston().velocity_m_s().powi(2);
        let outside = input().exterior_pressure_pa
            * (next.piston().gas().volume_m3() - old.gas().volume_m3());
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
fn laminar_sections_refresh_fixed_volume_endpoints_and_balance_both_stores() {
    scope(|cx, _| {
        for is_rectangular in [false, true] {
            let model = if is_rectangular {
                source::model_pr072()
            } else {
                source::model_transport(1.0, true)
            };
            let pressure = if is_rectangular {
                pr072_pressure()
            } else {
                100000.0
            };
            let old =
                IsochoricGasParcel::new(model, source::point(400.0, pressure), 0.003).unwrap();
            let saved = old.clone();
            let law = if is_rectangular {
                FilmGasHeatLaw::BulkForcedRectangularDuctSection(rectangular(
                    old.state(),
                    0.04,
                    0.12,
                ))
            } else {
                FilmGasHeatLaw::BulkForcedLaminarTubeSection(tube(old.state(), 0.2, 1.0, 0.02))
            };
            let area = match law {
                FilmGasHeatLaw::BulkForcedRectangularDuctSection(surface) => surface.area_m2,
                FilmGasHeatLaw::BulkForcedLaminarTubeSection(surface) => surface.area_m2,
                _ => unreachable!(),
            };
            let curve = source::body_curve();
            let thermal = UniformEnthalpyStepInput {
                surface_area_m2: area,
                duration_s: 0.02,
                internal_heat_j: 0.01,
                ..body(&curve, 120000.0)
            };
            let boundary = uniform_boundary();
            let mut environment = uniform_environment(&old, &boundary);
            environment.gas_heat_j = 0.02;
            environment.exchange.heat_bracket_j = [-0.5, 0.5];
            let next = FilmGasSurface {
                model: old.model(),
                law,
            }
            .isochoric_body_candidate(
                cx,
                &environment,
                thermal,
                WallConvectionSolve {
                    temperature_increment_bracket_k: [-15.0, 15.0],
                    temperature_tolerance_k: 1e-8,
                    maximum_iterations: 128,
                },
            )
            .unwrap();
            let report = &next.report;
            let endpoint = report.next_gas.state();
            let sample = &report.convection.source.source;
            assert_eq!(&sample.interface.gas, endpoint);
            assert_eq!(
                sample.property_gas.as_ref().unwrap().state().temperature,
                endpoint.state().temperature
            );
            assert_ne!(
                endpoint.state().temperature,
                old.state().state().temperature
            );
            let fresh = match law {
                FilmGasHeatLaw::BulkForcedRectangularDuctSection(surface) => {
                    evaluate_rectangular_duct_section_from_material(
                        endpoint,
                        surface.speed_m_s,
                        surface.hydraulic_diameter,
                        surface.aspect_ratio,
                        surface.start,
                        surface.end,
                    )
                    .unwrap()
                    .coefficient()
                    .qoi
                }
                FilmGasHeatLaw::BulkForcedLaminarTubeSection(surface) => {
                    evaluate_laminar_tube_section_from_material(
                        endpoint,
                        surface.speed_m_s,
                        surface.diameter,
                        surface.start,
                        surface.end,
                    )
                    .unwrap()
                    .coefficient()
                    .qoi
                }
                _ => unreachable!(),
            };
            assert_eq!(
                sample
                    .interface
                    .evaluation
                    .as_ref()
                    .unwrap()
                    .coefficient
                    .qoi,
                fresh
            );
            near(
                report.solid_to_gas_heat_j,
                thermal.duration_s
                    * area
                    * fresh
                    * (next.state.temperature_k() - endpoint.state().temperature),
                environment.exchange.source_tolerance_j,
            );
            let gas_energy = old.mass_kg()
                * (endpoint.specific_internal_energy_j_kg()
                    - old.state().specific_internal_energy_j_kg());
            let body_energy = thermal.mass_kg
                * (next.state.specific_enthalpy_j_kg() - thermal.initial.specific_enthalpy_j_kg());
            near(
                gas_energy,
                report.solid_to_gas_heat_j + environment.gas_heat_j,
                1e-7,
            );
            near(
                body_energy,
                thermal.internal_heat_j - report.solid_to_gas_heat_j,
                1e-7,
            );
            near(
                gas_energy + body_energy,
                environment.gas_heat_j + thermal.internal_heat_j,
                1e-7,
            );
            assert_eq!(old, saved);
        }
    });
}

#[test]
fn laminar_sections_refuse_unsupported_sources_geometry_and_regimes_without_mutation() {
    let curve = source::body_curve();
    let wall = body(&curve, 120000.0).initial;
    let model = source::model_transport(1.0, true);
    let bulk = model.evaluate(&source::point(400.0, 100000.0)).unwrap();
    let saved = bulk.clone();
    let valid = tube(&bulk, 0.2, 1.0, 0.02);
    let wrong_bank = source::model_transport(1.0, false);
    assert!(
        FilmGasSurface {
            model: &wrong_bank,
            law: FilmGasHeatLaw::BulkForcedLaminarTubeSection(valid)
        }
        .resolve(&bulk, wall)
        .is_err()
    );
    for patch in [
        ForcedTubeSectionGasSurface {
            speed_m_s: 0.0,
            ..valid
        },
        ForcedTubeSectionGasSurface {
            speed_m_s: 10.0 * valid.speed_m_s,
            ..valid
        },
        ForcedTubeSectionGasSurface {
            area_m2: 0.1,
            ..valid
        },
        ForcedTubeSectionGasSurface {
            start: Length::new(1e-12),
            ..valid
        },
        ForcedTubeSectionGasSurface {
            start: valid.end,
            ..valid
        },
    ] {
        assert!(
            FilmGasSurface {
                model: &model,
                law: FilmGasHeatLaw::BulkForcedLaminarTubeSection(patch)
            }
            .resolve(&bulk, wall)
            .is_err()
        );
    }
    let isolated = FilmGasSurface {
        model: &wrong_bank,
        law: FilmGasHeatLaw::BulkForcedLaminarTubeSection(ForcedTubeSectionGasSurface {
            area_m2: 0.0,
            speed_m_s: 0.0,
            ..valid
        }),
    }
    .resolve(&bulk, wall)
    .unwrap();
    assert_eq!(isolated.conductance_w_k, 0.0);
    assert!(
        isolated.source.property_gas.is_none() && isolated.source.interface.evaluation.is_none()
    );
    assert_eq!(bulk, saved);

    let rectangular_model = source::model_pr072();
    let rectangular_bulk = rectangular_model
        .evaluate(&source::point(400.0, pr072_pressure()))
        .unwrap();
    let valid = rectangular(&rectangular_bulk, 0.04, 0.12);
    for patch in [
        ForcedRectangularDuctSectionGasSurface {
            area_m2: 0.5 * valid.area_m2,
            ..valid
        },
        ForcedRectangularDuctSectionGasSurface {
            area_m2: 1.5 * valid.area_m2,
            ..valid
        },
        ForcedRectangularDuctSectionGasSurface {
            hydraulic_diameter: Length::new(0.0),
            ..valid
        },
        ForcedRectangularDuctSectionGasSurface {
            speed_m_s: 0.0,
            ..valid
        },
        rectangular(&rectangular_bulk, 1e-5, 0.12), // Unsupported positive entrance Gz.
    ] {
        assert!(
            FilmGasSurface {
                model: &rectangular_model,
                law: FilmGasHeatLaw::BulkForcedRectangularDuctSection(patch)
            }
            .resolve(&rectangular_bulk, wall)
            .is_err()
        );
    }
    let different_aspect = ForcedRectangularDuctSectionGasSurface {
        aspect_ratio: 0.4,
        area_m2: 0.01 * 1.4 * 1.4 / 0.4 * 0.08,
        ..valid
    };
    assert!(
        FilmGasSurface {
            model: &rectangular_model,
            law: FilmGasHeatLaw::BulkForcedRectangularDuctSection(different_aspect)
        }
        .resolve(&rectangular_bulk, wall)
        .is_err()
    );
    assert!(
        FilmGasSurface {
            model: &model,
            law: FilmGasHeatLaw::BulkForcedRectangularDuctSection(rectangular(&bulk, 0.04, 0.12))
        }
        .resolve(&bulk, wall)
        .is_err()
    ); // Actual Pr is not the table's Pr=0.72.
    let insulated = FilmGasSurface {
        model: &wrong_bank,
        law: FilmGasHeatLaw::BulkForcedRectangularDuctSection(
            ForcedRectangularDuctSectionGasSurface {
                area_m2: 0.0,
                speed_m_s: 0.0,
                ..valid
            },
        ),
    }
    .resolve(&rectangular_bulk, wall)
    .unwrap();
    assert_eq!(insulated.conductance_w_k, 0.0);
    assert!(insulated.source.property_gas.is_none());
}
