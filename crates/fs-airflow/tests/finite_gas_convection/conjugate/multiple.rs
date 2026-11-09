//! G1/G3/G4: distinct original conducting bodies sharing one finite gas window.
use super::*;
use fs_airflow::graph::thermal::medium::convection::finite_gas::conjugate::multiple::{
    ConjugateGasBody, adaptive_isochoric_bodies_candidate,
};
use fs_couple::finite_gas_thermal::{
    ConvectionCoefficient,
    bodies::{
        FiniteGasBodiesEnvironment, FiniteGasBodiesError, FiniteGasBody,
        window::{FiniteGasBodiesWindowError, FiniteGasBodiesWindowSolve},
    },
    window::FiniteGasWindowSolve,
};
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;

fn window_controls(duration: f64) -> FiniteGasBodiesWindowSolve {
    FiniteGasBodiesWindowSolve {
        time: FiniteGasWindowSolve {
            maximum_step_s: duration,
            minimum_step_s: 1e-6,
            maximum_attempts: 2048,
            energy_absolute_tolerance_j: 1e-3,
            temperature_absolute_tolerance_k: 1e-3,
            phase_fraction_absolute_tolerance: 1e-4,
            relative_tolerance: 0.0,
            energy_balance_tolerance_j: 1e-5,
        },
        maximum_retained_body_steps: 8192,
    }
}

fn wall_controls() -> WallConvectionSolve {
    WallConvectionSolve {
        temperature_increment_bracket_k: [-5.0, 5.0],
        temperature_tolerance_k: 1e-9,
        maximum_iterations: 128,
    }
}

fn shared_environment<'a>(
    gas: &'a IsochoricGasParcel,
    boundary: &'a LumpedThermalEnvironment,
) -> FiniteGasBodiesEnvironment<'a> {
    let single = uniform_environment(gas, boundary);
    FiniteGasBodiesEnvironment {
        gas,
        gas_heat_j: single.gas_heat_j,
        gas_solve: single.gas_solve,
        exchange: single.exchange,
        maximum_gas_temperature_increment_k: single.maximum_gas_temperature_increment_k,
        maximum_bodies: 4,
    }
}

#[test]
fn distinct_conducting_sources_share_each_actual_gas_endpoint_and_one_energy_account() {
    scope(|cx, _| {
        let hot_card = solid(true, 1.0);
        let cold_card = solid(true, 2.0);
        let hot_curve = caloric(&hot_card);
        let cold_curve = caloric(&cold_card);
        let hot_layer = layer(&hot_card, &hot_curve, [280.0, 440.0]);
        let cold_layer = layer(&cold_card, &cold_curve, [280.0, 440.0]);
        let geometry = ReservoirLayerGeometry::planar(AREA, DEPTH).unwrap();
        let old = IsochoricGasParcel::new(
            source::model_transport(1.0, true),
            source::point(350.0, 100000.0),
            0.003,
        )
        .unwrap();
        let saved = old.clone();
        let boundary = uniform_boundary();
        let mut environment = shared_environment(&old, &boundary);
        environment.gas_heat_j = 0.25;
        let hot = UniformEnthalpyStepInput {
            duration_s: 0.2,
            internal_heat_j: 0.6,
            ..thermal(&hot_curve, 420.0)
        };
        let cold = UniformEnthalpyStepInput {
            duration_s: 0.2,
            internal_heat_j: -0.2,
            ..thermal(&cold_curve, 315.0)
        };
        let film = plate(old.model())
            .with_prescribed_mass_flow(RATE, FLOW_AREA)
            .unwrap();
        let hot_surface = film
            .with_conducting_layer(&hot_layer, geometry, config())
            .unwrap();
        let cold_surface = film
            .with_conducting_layer(&cold_layer, geometry, config())
            .unwrap();
        let bindings = [
            ConjugateGasBody {
                body: FiniteGasBody {
                    id: 9,
                    environment: &boundary,
                    input: hot,
                },
                surface: &hot_surface,
            },
            ConjugateGasBody {
                body: FiniteGasBody {
                    id: 2,
                    environment: &boundary,
                    input: cold,
                },
                surface: &cold_surface,
            },
        ];
        let controls = window_controls(hot.duration_s);
        let next = adaptive_isochoric_bodies_candidate(
            cx,
            &environment,
            &bindings,
            controls,
            wall_controls(),
        )
        .unwrap();
        let report = next.report();
        assert!(report.rejected_comparisons > 0);
        assert!(report.maximum_accepted_normalized_error <= 1.0);
        assert_eq!(
            next.bodies()
                .iter()
                .map(|b| b.report.id)
                .collect::<Vec<_>>(),
            [2, 9]
        );
        let mut previous = 0.0;
        let mut total_heat = [0.0; 2];
        for sample in &report.samples {
            let dt = sample.time_s - previous;
            assert!(dt > 0.0);
            assert_eq!(sample.bodies.len(), 2);
            let mut layer_work = 0.0;
            for (slot, body) in sample.bodies.iter().enumerate() {
                let feedback = &body.report.convection.source;
                let endpoint = &feedback.endpoint_source;
                let response = endpoint.layer.as_ref().unwrap();
                let original = if body.report.id == 9 {
                    &hot_layer
                } else {
                    &cold_layer
                };
                let k_scale = if body.report.id == 9 { 1.0 } else { 2.0 };
                assert!(Arc::ptr_eq(
                    &endpoint.conductivity_source,
                    original.source_line()
                ));
                assert_eq!(endpoint.body_temperature_k, body.state.temperature_k());
                assert_eq!(endpoint.bulk_temperature_k, sample.gas_temperature_k);
                let film = &response.film;
                let gas = &film.interface.gas;
                let property = film.property_gas.as_ref().unwrap();
                assert_eq!(gas.state().temperature, sample.gas_temperature_k);
                assert_eq!(property.source_packs(), gas.source_packs());
                near(
                    property.state().temperature,
                    0.5 * (response.surface_temperature_k + sample.gas_temperature_k),
                    1e-12,
                );
                let h = plate_reference(property, gas.state().density);
                let film_power =
                    AREA * h * (response.surface_temperature_k - sample.gas_temperature_k);
                // Independent integral of the manufactured original k(T)=scale*0.002*T.
                let solid_power = (AREA / DEPTH)
                    * k_scale
                    * 0.001
                    * (endpoint.body_temperature_k - response.surface_temperature_k)
                    * (endpoint.body_temperature_k + response.surface_temperature_k);
                near(response.power_to_bath_w, film_power, 2e-10);
                near(response.solid_power_to_bath_w, solid_power, 2e-9);
                near(
                    body.report.solid_to_gas_heat_j,
                    dt * film_power,
                    2.0 * environment.exchange.source_tolerance_j,
                );
                assert!(
                    feedback.boundary_residual_j.abs()
                        <= environment.exchange.source_tolerance_j / 2.0
                );
                layer_work +=
                    dt * (response.heat_residual_w.abs() + endpoint.representation_error_w.abs());
                total_heat[slot] += body.report.solid_to_gas_heat_j;
            }
            assert!(layer_work <= environment.exchange.source_tolerance_j);
            near(
                sample.shared.heat_from_bodies_j,
                sample
                    .bodies
                    .iter()
                    .map(|b| b.report.solid_to_gas_heat_j)
                    .sum(),
                1e-12,
            );
            previous = sample.time_s;
        }
        assert_eq!(previous, hot.duration_s);
        assert!(total_heat[0] < 0.0 && total_heat[1] > 0.0);
        let cold_storage = cold.mass_kg
            * 1000.0
            * (next.bodies()[0].state.temperature_k() - cold.initial.temperature_k());
        let hot_storage = hot.mass_kg
            * 1000.0
            * (next.bodies()[1].state.temperature_k() - hot.initial.temperature_k());
        let gas_storage = old.mass_kg()
            * (next.gas().state().specific_internal_energy_j_kg()
                - old.state().specific_internal_energy_j_kg());
        near(
            hot_storage + cold_storage + gas_storage,
            hot.internal_heat_j + cold.internal_heat_j + environment.gas_heat_j,
            controls.time.energy_balance_tolerance_j,
        );
        for (index, body) in next.bodies().iter().enumerate() {
            near(body.report.solid_to_gas_heat_j, total_heat[index], 1e-12);
            near(body.report.heat_apportionment_residual_j, 0.0, 1e-14);
        }
        near(report.external_gas_heat_j, environment.gas_heat_j, 1e-14);
        near(report.entropy.balance_residual_j_k, 0.0, 1e-10);
        assert!(report.sum_absolute_energy_defect_j <= controls.time.energy_balance_tolerance_j);
        assert_eq!(next.gas().mass_kg(), old.mass_kg());
        let reversed = adaptive_isochoric_bodies_candidate(
            cx,
            &environment,
            &[bindings[1], bindings[0]],
            controls,
            wall_controls(),
        )
        .unwrap();
        assert_eq!(next.gas(), reversed.gas());
        assert_eq!(next.bodies(), reversed.bodies());
        assert_eq!(old, saved);
    });
}

#[test]
fn shared_conjugate_window_agrees_with_two_original_simultaneous_half_transactions() {
    scope(|cx, _| {
        let card = solid(true, 1.0);
        let curve = caloric(&card);
        let material = layer(&card, &curve, [280.0, 440.0]);
        let old = IsochoricGasParcel::new(
            source::model_transport(1.0, true),
            source::point(350.0, 100000.0),
            0.003,
        )
        .unwrap();
        let boundary = uniform_boundary();
        let mut environment = shared_environment(&old, &boundary);
        environment.gas_heat_j = 0.2;
        let surface = plate(old.model())
            .with_prescribed_mass_flow(RATE, FLOW_AREA)
            .unwrap()
            .with_conducting_layer(
                &material,
                ReservoirLayerGeometry::planar(AREA, DEPTH).unwrap(),
                config(),
            )
            .unwrap();
        let bodies = [
            FiniteGasBody {
                id: 7,
                environment: &boundary,
                input: UniformEnthalpyStepInput {
                    duration_s: 0.1,
                    internal_heat_j: 0.3,
                    ..thermal(&curve, 420.0)
                },
            },
            FiniteGasBody {
                id: 3,
                environment: &boundary,
                input: UniformEnthalpyStepInput {
                    duration_s: 0.1,
                    internal_heat_j: -0.1,
                    ..thermal(&curve, 315.0)
                },
            },
        ];
        let bindings = bodies.map(|body| ConjugateGasBody {
            body,
            surface: &surface,
        });
        let controls = window_controls(0.1);
        let loose = FiniteGasBodiesWindowSolve {
            time: FiniteGasWindowSolve {
                energy_absolute_tolerance_j: 100.0,
                temperature_absolute_tolerance_k: 100.0,
                phase_fraction_absolute_tolerance: 1.0,
                ..controls.time
            },
            ..controls
        };
        let next = adaptive_isochoric_bodies_candidate(
            cx,
            &environment,
            &bindings,
            loose,
            wall_controls(),
        )
        .unwrap();
        assert_eq!(next.report().attempted_intervals, 3);
        assert_eq!(next.report().samples.len(), 2);
        let mut half = bodies.map(|body| FiniteGasBody {
            input: UniformEnthalpyStepInput {
                duration_s: 0.05,
                internal_heat_j: 0.5 * body.input.internal_heat_j,
                ..body.input
            },
            ..body
        });
        let half_environment = FiniteGasBodiesEnvironment {
            gas_heat_j: 0.5 * environment.gas_heat_j,
            ..environment.clone()
        };
        // Use the original shared native owner and original public material
        // resolver, without the new collection adapter or adaptive controller.
        let resolve = |body: &FiniteGasBody<'_>, gas: &ResolvedThermochemicalGasState, t: f64| {
            let state = body
                .input
                .curve
                .state_at_specific_enthalpy(1000.0 * (t - 200.0))
                .unwrap();
            let sample = surface.resolve(gas, state, || Ok(())).unwrap();
            Ok::<_, core::convert::Infallible>(ConvectionCoefficient {
                coefficient_w_per_m2_k: sample.conductance_w_k / body.input.surface_area_m2,
                source: (),
            })
        };
        let first = half_environment
            .advance_with_wall_convection(cx, &half, wall_controls(), resolve)
            .unwrap();
        for body in &mut half {
            body.input.initial = first
                .bodies()
                .iter()
                .find(|candidate| candidate.report.id == body.id)
                .unwrap()
                .state;
        }
        let second_environment = FiniteGasBodiesEnvironment {
            gas: first.gas(),
            ..half_environment
        };
        let second = second_environment
            .advance_with_wall_convection(cx, &half, wall_controls(), resolve)
            .unwrap();
        near(
            next.gas().state().state().temperature,
            second.gas().state().state().temperature,
            1e-10,
        );
        for (actual, expected) in next.bodies().iter().zip(second.bodies()) {
            assert_eq!(actual.report.id, expected.report.id);
            near(
                actual.state.temperature_k(),
                expected.state.temperature_k(),
                1e-10,
            );
            near(
                actual.state.specific_enthalpy_j_kg(),
                expected.state.specific_enthalpy_j_kg(),
                1e-7,
            );
        }
        near(
            next.report().heat_from_bodies_j,
            first.report().heat_from_bodies_j + second.report().heat_from_bodies_j,
            2e-9,
        );
    });
}

#[test]
fn every_shared_body_keeps_original_source_identity_and_atomic_refusal() {
    scope(|cx, gate| {
        let card = solid(true, 1.0);
        let curve = caloric(&card);
        let cloned = curve.clone();
        let limited = layer(&card, &curve, [300.0, 400.0]);
        let old = IsochoricGasParcel::new(
            source::model_transport(1.0, true),
            source::point(350.0, 100000.0),
            0.003,
        )
        .unwrap();
        let saved = old.clone();
        let boundary = uniform_boundary();
        let environment = shared_environment(&old, &boundary);
        let insulated = FilmGasSurface {
            model: old.model(),
            law: FilmGasHeatLaw::Forced(ForcedGasSurface {
                area_m2: 0.0,
                speed_m_s: 0.0,
                ..surface()
            }),
        }
        .with_conducting_layer(
            &limited,
            ReservoirLayerGeometry::planar(AREA, DEPTH).unwrap(),
            config(),
        )
        .unwrap();
        let safe = FiniteGasBody {
            id: 1,
            environment: &boundary,
            input: UniformEnthalpyStepInput {
                duration_s: 0.1,
                ..thermal(&curve, 350.0)
            },
        };
        let warming = FiniteGasBody {
            id: 8,
            environment: &boundary,
            input: UniformEnthalpyStepInput {
                duration_s: 0.1,
                internal_heat_j: 4.0,
                ..thermal(&curve, 399.9)
            },
        };
        let bindings = [
            ConjugateGasBody {
                body: safe,
                surface: &insulated,
            },
            ConjugateGasBody {
                body: warming,
                surface: &insulated,
            },
        ];
        let controls = window_controls(0.1);
        assert!(matches!(
            adaptive_isochoric_bodies_candidate(
                cx,
                &environment,
                &bindings,
                controls,
                wall_controls()
            ),
            Err(FiniteGasBodiesWindowError::Step(
                FiniteGasBodiesError::Body {
                    id: 8,
                    error: ConvectionError::Convection(WallConvectionError::Source(
                        ConjugateGasError::Layer(_)
                    ))
                }
            ))
        ));
        let mut wrong_identity = bindings;
        wrong_identity[1].body.input = UniformEnthalpyStepInput {
            curve: &cloned,
            internal_heat_j: 0.0,
            ..warming.input
        };
        assert!(matches!(
            adaptive_isochoric_bodies_candidate(
                cx,
                &environment,
                &wrong_identity,
                controls,
                wall_controls()
            ),
            Err(FiniteGasBodiesWindowError::Step(
                FiniteGasBodiesError::Body {
                    id: 8,
                    error: ConvectionError::Convection(WallConvectionError::Source(
                        ConjugateGasError::Invalid(_)
                    ))
                }
            ))
        ));
        assert!(matches!(
            adaptive_isochoric_bodies_candidate(
                cx,
                &environment,
                &[bindings[0], bindings[0]],
                controls,
                wall_controls()
            ),
            Err(FiniteGasBodiesWindowError::Invalid(_))
        ));
        let mut corrected = bindings;
        corrected[1].body.input.internal_heat_j = 0.0;
        let after_refusal = adaptive_isochoric_bodies_candidate(
            cx,
            &environment,
            &corrected,
            controls,
            wall_controls(),
        )
        .unwrap();
        assert_eq!(after_refusal.gas(), &old);
        assert_eq!(after_refusal.bodies()[0].state, safe.input.initial);
        assert_eq!(after_refusal.bodies()[1].state, warming.input.initial);
        gate.request();
        assert!(matches!(
            adaptive_isochoric_bodies_candidate(
                cx,
                &environment,
                &corrected,
                controls,
                wall_controls()
            ),
            Err(FiniteGasBodiesWindowError::Cancelled)
        ));
        assert_eq!(old, saved);
        near(warming.input.initial.temperature_k(), 399.9, 1e-12);
    });
}
