//! G1/G3/G4: adaptive physical time with actual source wall and conducting surface.
use super::*;
use fs_couple::finite_gas_thermal::window::{FiniteGasWindowError, FiniteGasWindowSolve};

fn time_controls(duration: f64) -> FiniteGasWindowSolve {
    FiniteGasWindowSolve {
        maximum_step_s: duration,
        minimum_step_s: 1e-6,
        maximum_attempts: 2048,
        energy_absolute_tolerance_j: 1e-3,
        temperature_absolute_tolerance_k: 1e-3,
        phase_fraction_absolute_tolerance: 1e-4,
        relative_tolerance: 0.0,
        energy_balance_tolerance_j: 1e-5,
    }
}

fn wall_controls() -> WallConvectionSolve {
    WallConvectionSolve {
        temperature_increment_bracket_k: [-5.0, 5.0],
        temperature_tolerance_k: 1e-9,
        maximum_iterations: 128,
    }
}

#[test]
fn adaptive_source_film_refreshes_wall_mass_flux_and_spends_each_prescribed_joule_once() {
    scope(|cx, gate| {
        let card = solid(true, 1.0);
        let curve = caloric(&card);
        let old = IsochoricGasParcel::new(
            source::model_transport(1.0, true),
            source::point(350.0, 100000.0),
            0.003,
        )
        .unwrap();
        let saved = old.clone();
        let boundary = uniform_boundary();
        let mut environment = uniform_environment(&old, &boundary);
        environment.gas_heat_j = 0.3;
        let body = UniformEnthalpyStepInput {
            duration_s: 0.2,
            internal_heat_j: 1.0,
            ..thermal(&curve, 420.0)
        };
        let surface = plate(old.model())
            .with_prescribed_mass_flow(RATE, FLOW_AREA)
            .unwrap();
        let controls = time_controls(body.duration_s);
        let next = environment
            .advance_adaptive_with_wall_convection(
                cx,
                body,
                controls,
                wall_controls(),
                |gas, temperature| {
                    surface.resolve_isochoric_coefficient(
                        gas,
                        body.surface_area_m2,
                        Temperature::new(temperature),
                    )
                },
            )
            .unwrap();
        let report = &next.report;
        assert!(report.rejected_comparisons > 0);
        assert!(report.samples.len() > 2);
        assert!(report.maximum_accepted_normalized_error <= 1.0);
        assert!(report.attempted_intervals <= controls.maximum_attempts);
        let mut previous = 0.0;
        let mut property_temperatures = Vec::new();
        for sample in &report.samples {
            let dt = sample.time_s - previous;
            assert!(dt > 0.0);
            let feedback = &sample.convection.source;
            let film = &feedback.source;
            let bulk = &film.interface.gas;
            let property = film.property_gas.as_ref().unwrap();
            let flow = film.prescribed_mass_flow.unwrap();
            assert_eq!(bulk.state().temperature, sample.gas_temperature_k);
            assert_eq!(film.interface.body.value(), feedback.temperature_k);
            assert_eq!(
                feedback.endpoint_source.interface.body.value(),
                sample.body_state.temperature_k()
            );
            near(
                property.state().temperature,
                0.5 * (feedback.temperature_k + sample.gas_temperature_k),
                1e-12,
            );
            assert_eq!(property.state().pressure, bulk.state().pressure);
            assert_eq!(property.source_packs(), bulk.source_packs());
            assert_eq!(flow.bulk_density_kg_m3, bulk.state().density);
            near(flow.bulk_density_kg_m3, old.state().state().density, 2e-12);
            assert!((property.state().density / bulk.state().density - 1.0).abs() > 0.05);
            near(
                flow.speed_m_s * flow.bulk_density_kg_m3 * FLOW_AREA,
                RATE,
                2e-15,
            );
            let h = plate_reference(property, bulk.state().density);
            near(sample.convection.coefficient_w_per_m2_k, h, 2e-12 * h);
            near(
                sample.solid_to_gas_heat_j,
                dt * AREA * h * (sample.body_state.temperature_k() - sample.gas_temperature_k),
                environment.exchange.source_tolerance_j,
            );
            assert!(
                feedback.temperature_residual_k.abs() <= wall_controls().temperature_tolerance_k
            );
            assert!(feedback.boundary_residual_j.abs() <= environment.exchange.source_tolerance_j);
            property_temperatures.push(property.state().temperature);
            previous = sample.time_s;
        }
        assert_eq!(previous, body.duration_s);
        assert!((property_temperatures.last().unwrap() - property_temperatures[0]).abs() > 0.1);
        assert!((next.state.temperature_k() - body.initial.temperature_k()).abs() > 0.05);
        near(report.internal_body_heat_j, body.internal_heat_j, 1e-14);
        near(report.external_gas_heat_j, environment.gas_heat_j, 1e-14);
        near(report.body_heat_apportionment_residual_j, 0.0, 1e-14);
        near(report.gas_heat_apportionment_residual_j, 0.0, 1e-14);
        let dh = body.mass_kg
            * (next.state.specific_enthalpy_j_kg() - body.initial.specific_enthalpy_j_kg());
        let du = old.mass_kg()
            * (report.next_gas.state().specific_internal_energy_j_kg()
                - old.state().specific_internal_energy_j_kg());
        near(
            dh + du,
            body.internal_heat_j + environment.gas_heat_j,
            controls.energy_balance_tolerance_j,
        );
        near(
            dh,
            body.internal_heat_j - report.solid_to_gas_heat_j,
            controls.energy_balance_tolerance_j,
        );
        assert!(report.sum_absolute_energy_defect_j <= controls.energy_balance_tolerance_j);
        assert_eq!(report.radiation_heat_j, 0.0);
        assert_eq!(report.next_gas.mass_kg(), old.mass_kg());
        near(report.entropy.balance_residual_j_k, 0.0, 1e-10);

        let exhausted = environment.advance_adaptive_with_wall_convection(
            cx,
            body,
            FiniteGasWindowSolve {
                maximum_attempts: 1,
                ..controls
            },
            wall_controls(),
            |gas, temperature| {
                surface.resolve_isochoric_coefficient(
                    gas,
                    body.surface_area_m2,
                    Temperature::new(temperature),
                )
            },
        );
        assert!(matches!(
            exhausted,
            Err(FiniteGasWindowError::WorkLimit { .. })
        ));
        let cancelled = environment.advance_adaptive_with_wall_convection(
            cx,
            body,
            controls,
            wall_controls(),
            |gas, temperature| {
                let value = surface.resolve_isochoric_coefficient(
                    gas,
                    body.surface_area_m2,
                    Temperature::new(temperature),
                );
                gate.request();
                value
            },
        );
        assert!(matches!(cancelled, Err(FiniteGasWindowError::Cancelled)));
        assert_eq!(old, saved);
        assert_eq!(body.initial.temperature_k(), 420.0);
    });
}

#[test]
fn adaptive_conducting_surface_reuses_native_intervals_and_refuses_actual_source_exit() {
    scope(|cx, _| {
        let card = solid(true, 1.0);
        let curve = caloric(&card);
        let wall = layer(&card, &curve, [300.0, 440.0]);
        let geometry = ReservoirLayerGeometry::planar(AREA, DEPTH).unwrap();
        let old = IsochoricGasParcel::new(
            source::model_transport(1.0, true),
            source::point(350.0, 100000.0),
            0.003,
        )
        .unwrap();
        let saved = old.clone();
        let boundary = uniform_boundary();
        let mut environment = uniform_environment(&old, &boundary);
        environment.gas_heat_j = 0.2;
        let body = UniformEnthalpyStepInput {
            duration_s: 0.2,
            internal_heat_j: 0.6,
            ..thermal(&curve, 420.0)
        };
        let owner = plate(old.model())
            .with_prescribed_mass_flow(RATE, FLOW_AREA)
            .unwrap()
            .with_conducting_layer(&wall, geometry, config())
            .unwrap();
        let controls = time_controls(body.duration_s);
        let next = owner
            .adaptive_isochoric_body_candidate(cx, &environment, body, controls, wall_controls())
            .unwrap();
        let report = &next.report;
        assert!(report.rejected_comparisons > 0);
        assert!(report.maximum_accepted_normalized_error <= 1.0);
        let mut previous = 0.0;
        for sample in &report.samples {
            let dt = sample.time_s - previous;
            assert!(dt > 0.0);
            let feedback = &sample.convection.source;
            let endpoint = &feedback.endpoint_source;
            let response = endpoint.layer.as_ref().unwrap();
            let film = &response.film;
            let property = film.property_gas.as_ref().unwrap();
            let bulk = &film.interface.gas;
            assert_eq!(
                endpoint.body_temperature_k,
                sample.body_state.temperature_k()
            );
            assert_eq!(bulk.state().temperature, sample.gas_temperature_k);
            assert_eq!(film.interface.body.value(), response.surface_temperature_k);
            assert!(endpoint.body_temperature_k > response.surface_temperature_k + 0.1);
            near(
                property.state().temperature,
                0.5 * (response.surface_temperature_k + sample.gas_temperature_k),
                1e-12,
            );
            assert_eq!(property.source_packs(), bulk.source_packs());
            assert!(Arc::ptr_eq(
                &endpoint.conductivity_source,
                wall.source_line()
            ));
            let h = plate_reference(property, bulk.state().density);
            let film_power = AREA * h * (response.surface_temperature_k - sample.gas_temperature_k);
            let solid_power = (AREA / DEPTH)
                * 0.001
                * (endpoint.body_temperature_k - response.surface_temperature_k)
                * (endpoint.body_temperature_k + response.surface_temperature_k);
            near(response.power_to_bath_w, film_power, 2e-10);
            near(response.solid_power_to_bath_w, solid_power, 2e-9);
            near(
                sample.solid_to_gas_heat_j,
                dt * film_power,
                2.0 * environment.exchange.source_tolerance_j,
            );
            for source in [&feedback.source, endpoint] {
                assert!(
                    (dt * source.layer.as_ref().unwrap().heat_residual_w).abs()
                        <= environment.exchange.source_tolerance_j
                );
            }
            previous = sample.time_s;
        }
        assert_eq!(previous, body.duration_s);
        near(report.internal_body_heat_j, body.internal_heat_j, 1e-14);
        near(report.external_gas_heat_j, environment.gas_heat_j, 1e-14);
        near(
            report.body_enthalpy_change_j + report.gas_internal_energy_change_j,
            body.internal_heat_j + environment.gas_heat_j,
            controls.energy_balance_tolerance_j,
        );
        assert!(report.sum_absolute_energy_defect_j <= controls.energy_balance_tolerance_j);
        assert_eq!(report.next_gas.mass_kg(), old.mass_kg());

        // A loose one-comparison window must return the exact original two
        // half transactions, including half-window heater/gas source energies.
        let loose = owner
            .adaptive_isochoric_body_candidate(
                cx,
                &environment,
                body,
                FiniteGasWindowSolve {
                    energy_absolute_tolerance_j: 100.0,
                    temperature_absolute_tolerance_k: 100.0,
                    phase_fraction_absolute_tolerance: 1.0,
                    ..controls
                },
                wall_controls(),
            )
            .unwrap();
        assert_eq!(loose.report.attempted_intervals, 3);
        assert_eq!(loose.report.samples.len(), 2);
        let half = UniformEnthalpyStepInput {
            duration_s: 0.5 * body.duration_s,
            internal_heat_j: 0.5 * body.internal_heat_j,
            ..body
        };
        let half_environment = FiniteGasThermalEnvironment {
            gas_heat_j: 0.5 * environment.gas_heat_j,
            ..environment.clone()
        };
        let first = owner
            .isochoric_body_candidate(cx, &half_environment, half, wall_controls())
            .unwrap();
        let second_environment = FiniteGasThermalEnvironment {
            gas: &first.report.next_gas,
            ..half_environment
        };
        let second = owner
            .isochoric_body_candidate(
                cx,
                &second_environment,
                UniformEnthalpyStepInput {
                    initial: first.state,
                    ..half
                },
                wall_controls(),
            )
            .unwrap();
        assert_eq!(loose.state, second.state);
        assert_eq!(loose.report.next_gas, second.report.next_gas);
        near(
            loose.report.solid_to_gas_heat_j,
            first.report.solid_to_gas_heat_j + second.report.solid_to_gas_heat_j,
            1e-12,
        );

        // Insulation skips unused gas-property roots, but the coarse proposal's
        // actual final native body still leaves the pinned conductivity source.
        let limited = layer(&card, &curve, [300.0, 400.0]);
        let insulated = FilmGasSurface {
            model: old.model(),
            law: FilmGasHeatLaw::Forced(ForcedGasSurface {
                area_m2: 0.0,
                speed_m_s: 0.0,
                ..surface()
            }),
        }
        .with_conducting_layer(&limited, geometry, config())
        .unwrap();
        let warming = UniformEnthalpyStepInput {
            duration_s: 0.1,
            internal_heat_j: 4.0,
            ..thermal(&curve, 399.9)
        };
        assert!(matches!(
            insulated.adaptive_isochoric_body_candidate(
                cx,
                &uniform_environment(&old, &boundary),
                warming,
                time_controls(warming.duration_s),
                wall_controls(),
            ),
            Err(FiniteGasWindowError::Step(ConvectionError::Convection(
                WallConvectionError::Source(ConjugateGasError::Layer(_)),
            )))
        ));
        let cloned = curve.clone();
        assert!(matches!(
            owner.adaptive_isochoric_body_candidate(
                cx,
                &environment,
                UniformEnthalpyStepInput {
                    curve: &cloned,
                    ..body
                },
                controls,
                wall_controls(),
            ),
            Err(FiniteGasWindowError::Step(ConvectionError::Convection(
                WallConvectionError::Source(ConjugateGasError::Invalid(_)),
            )))
        ));
        assert_eq!(old, saved);
        near(warming.initial.temperature_k(), 399.9, 1e-12);
    });
}
