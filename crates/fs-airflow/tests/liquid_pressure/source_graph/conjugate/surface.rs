//! Independent film formulas against the original finite graph/material path.
use super::*;
use fs_airflow::graph::thermal::medium::convection::NaturalConvectionGeometry;
use fs_airflow::graph::thermal::transport::transient::material::convection::PlateFlowDirection;

fn source_binding<'a, 'm>(
    layer: &'a ResolvedConductionLayer<'m>,
    kind: usize,
) -> GraphConjugateLinkBinding<'a, 'm> {
    let mut value = binding(layer);
    value.convection.property_temperature = ConvectionPropertyTemperature::FilmMean;
    if kind == 1 {
        value.convection.property_temperature = ConvectionPropertyTemperature::Reservoir;
        value.convection.law = GraphConvectionLaw::ForcedSphere { branch: 0 };
    } else if kind == 2 {
        value.convection.characteristic_length = Length::new(0.05);
        value.convection.law = GraphConvectionLaw::ForcedPlateSection {
            branch: 0,
            correlation: CorrelationId::FlatPlateLaminarAverage,
            start: Length::new(0.02),
            direction: PlateFlowDirection::AlongBranch,
        };
    }
    value
}

// Manufactured original laws: k_liquid(T)=0.005T-1; mu(T)=0.005-1e-5T.
// These equations do not call the production correlations or surface solver.
fn film_conductance(kind: usize, surface: f64) -> f64 {
    let property = if kind == 1 {
        300.0
    } else {
        0.5 * (surface + 300.0)
    };
    let k = 0.005 * property - 1.0;
    let mu = 0.005 - 1e-5 * property;
    let pr = mu * 2000.0 / k;
    let h = match kind {
        0 => 3.66 * k / 0.03,
        1 => {
            let re: f64 = 1000.0 * 0.01 * 0.03 / mu;
            let surface_mu = 0.005 - 1e-5 * surface;
            (2.0 + (0.4 * re.sqrt() + 0.06 * re.powf(2.0 / 3.0))
                * pr.powf(0.4)
                * (mu / surface_mu).powf(0.25))
                * k
                / 0.03
        }
        _ => {
            0.664 * k * (1000.0 * 0.01 / mu).sqrt() * pr.powf(1.0 / 3.0)
                / (0.02_f64.sqrt() + 0.05_f64.sqrt())
        }
    };
    AREA * h
}

#[test]
fn film_duct_sphere_and_plate_use_the_solved_surface_and_original_sources() {
    context(|cx| {
        let card = fluid_card();
        let (law, caloric) = source_pair(&card);
        let solid = solid(true);
        let wall = wall_caloric(&solid);
        let point = QueryPoint::new().with("T", 300.0).unwrap();
        let layer = layer(&solid, &wall, &point, [200.0, 500.0]);
        for kind in 0..3 {
            for tw in [400.0_f64, 280.0, 300.0] {
                // Whitaker's original liquid viscosity ratio excludes cooling.
                if kind == 1 && tw < 300.0 {
                    continue;
                }
                let initial = initial(&caloric, &wall, tw);
                let graph = graph();
                let inventory = bind_inventory(cx, &graph, &initial);
                let laws = laws(2300.0);
                let owner = LiquidGraphSources::try_new(
                    cx,
                    &law,
                    &caloric,
                    &initial,
                    SourceGraphBinding {
                        graph: &graph,
                        inventory: &inventory,
                        laws: &laws,
                        acceleration_m_s2: [0.0; 3],
                    },
                    caps(),
                )
                .unwrap();
                let flows = [1e-5, 0.0, 1e-5];
                let source = owner.resolve(cx, &initial, &flows).unwrap();
                let report = native_report(
                    cx,
                    &source.graph,
                    &inventory,
                    flows,
                    [1.0, 0.0, 1.0],
                    flows,
                    0.1,
                );
                let prepared = PreparedGraphConjugateLinks::new_with_surface_feedback(
                    &initial,
                    vec![source_binding(&layer, kind)],
                    1,
                    || Ok(()),
                )
                .unwrap();
                let heat = owner
                    .evaluate_conjugate_links(
                        cx,
                        &initial,
                        &source.sources,
                        &report,
                        &prepared,
                        config(),
                    )
                    .unwrap();
                let sample = &heat.samples()[0];
                let ts = sample.layer.surface_temperature_k;
                let film = &sample.layer.film;
                let mut lo = tw.min(300.0);
                let mut hi = tw.max(300.0);
                for _ in 0..80 {
                    let t = 0.5 * (lo + hi);
                    let solid_q = AREA / DISTANCE * 0.001 * (tw - t) * (tw + t);
                    if solid_q - film_conductance(kind, t) * (t - 300.0) > 0.0 {
                        lo = t;
                    } else {
                        hi = t;
                    }
                }
                near(ts, 0.5 * (lo + hi));
                near(film.wall_temperature_k, ts);
                near(sample.wall_temperature_k, tw);
                near(film.conductance_w_k, film_conductance(kind, ts));
                assert!((sample.layer.heat_residual_w).abs() < 2e-10);
                assert!(
                    (sample.applied_heat_into_fluid_w
                        - AREA / DISTANCE * 0.001 * (tw - ts) * (tw + ts))
                        .abs()
                        < 2e-10
                );
                assert_eq!(
                    film.material.material().card_identity(),
                    card.content_hash()
                );
                assert!(Arc::ptr_eq(&sample.source_line, layer.source_line()));
                if kind == 1 {
                    near(film.property_temperature_k, 300.0);
                    let surface = film.surface.as_ref().unwrap();
                    near(surface.thermal.temperature_k(), ts);
                    assert_eq!(
                        surface.material.material().card_identity(),
                        card.content_hash()
                    );
                } else {
                    near(film.property_temperature_k, 0.5 * (ts + 300.0));
                }
                if tw != 300.0 {
                    assert_ne!(ts, tw);
                }
                let mut coefficients = BTreeMap::new();
                heat.apply_to(&mut coefficients, || Ok(())).unwrap();
                let before = coefficients.clone();
                assert!(heat.apply_to(&mut coefficients, || Ok(())).is_err());
                assert_eq!(coefficients, before);
                let mut cancelled = BTreeMap::new();
                assert!(
                    heat.apply_to(&mut cancelled, || Err(DepositionError::Cancelled))
                        .is_err()
                );
                assert!(cancelled.is_empty());
                let trial = BackwardEulerStep::new(&initial, 0.1, || Ok(()))
                    .unwrap()
                    .evaluate(
                        &[0.0; 7],
                        |_| {
                            Ok::<_, DepositionError>(StepCoefficients {
                                conductances_w_k: coefficients,
                                sources_w: BTreeMap::new(),
                                payload: (),
                            })
                        },
                        || Ok(()),
                    )
                    .unwrap();
                near(
                    trial.residual_j()[2],
                    -0.1 * sample.applied_heat_into_fluid_w,
                );
                near(
                    trial.residual_j()[6],
                    0.1 * sample.applied_heat_into_fluid_w,
                );
                near(trial.residual_j().iter().sum(), 0.0);
                if kind == 0 && tw == 400.0 {
                    let source = owner.resolve(cx, &initial, &[1e-5, 1e-5, -1e-5]).unwrap();
                    let flow = source
                        .graph
                        .solve(
                            &[
                                FixedPressure {
                                    node: 0,
                                    pressure: Pressure::new(1.0),
                                },
                                FixedPressure {
                                    node: 1,
                                    pressure: Pressure::new(0.0),
                                },
                            ],
                            GraphSolveConfig {
                                max_sweeps: 256,
                                max_node_iterations: 96,
                                absolute_flow_tolerance: VolumetricFlowRate::new(1e-14),
                                relative_flow_tolerance: 1e-12,
                            },
                            cx,
                        )
                        .unwrap();
                    let steady = owner
                        .evaluate_quasi_steady_conjugate_links(
                            cx,
                            &initial,
                            &source.sources,
                            &flow,
                            0.1,
                            &prepared,
                            config(),
                        )
                        .unwrap();
                    near(
                        steady.samples()[0].applied_heat_into_fluid_w,
                        sample.applied_heat_into_fluid_w,
                    );
                }
            }
        }
    });
}

#[test]
fn source_limits_and_unclosed_surface_refuse_without_mutating_native_history() {
    context(|cx| {
        for source_limit in [false, true] {
            let card = if source_limit {
                card_with_transport_ceiling(true, 1.0, 330.0)
            } else {
                fluid_card()
            };
            let (law, caloric) = source_pair(&card);
            let solid = solid(true);
            let wall = wall_caloric(&solid);
            let point = QueryPoint::new().with("T", 300.0).unwrap();
            let layer = layer(&solid, &wall, &point, [200.0, 500.0]);
            let initial = initial(&caloric, &wall, 400.0);
            let graph = graph();
            let inventory = bind_inventory(cx, &graph, &initial);
            let laws = laws(2300.0);
            let owner = LiquidGraphSources::try_new(
                cx,
                &law,
                &caloric,
                &initial,
                SourceGraphBinding {
                    graph: &graph,
                    inventory: &inventory,
                    laws: &laws,
                    acceleration_m_s2: [0.0; 3],
                },
                caps(),
            )
            .unwrap();
            let flows = [1e-5, 0.0, 1e-5];
            let source = owner.resolve(cx, &initial, &flows).unwrap();
            let report = native_report(
                cx,
                &source.graph,
                &inventory,
                flows,
                [1.0, 0.0, 1.0],
                flows,
                0.1,
            );
            let prepared = PreparedGraphConjugateLinks::new_with_surface_feedback(
                &initial,
                vec![source_binding(&layer, 0)],
                1,
                || Ok(()),
            )
            .unwrap();
            let controls = LayerConductionConfig {
                max_iterations: 1,
                absolute_heat_tolerance_w: 0.0,
                relative_heat_tolerance: 0.0,
            };
            let failure = owner
                .evaluate_conjugate_links(
                    cx,
                    &initial,
                    &source.sources,
                    &report,
                    &prepared,
                    controls,
                )
                .unwrap_err();
            if source_limit {
                assert!(matches!(failure, GraphConjugateError::Convection(_)));
            } else {
                assert!(matches!(
                    failure,
                    GraphConjugateError::SurfaceBalance { link: 9, .. }
                ));
                // The old bulk-only callback entrypoint cannot bypass the
                // explicit surface feedback path by reusing this preparation.
                let called = std::cell::Cell::new(false);
                assert!(
                    prepared
                        .resolve(
                            &initial,
                            config(),
                            |candidate, links| {
                                called.set(true);
                                owner.evaluate_links(
                                    cx,
                                    candidate,
                                    &source.sources,
                                    &report,
                                    links,
                                    1,
                                )
                            },
                            || Ok(())
                        )
                        .is_err()
                );
                assert!(!called.get());
                let narrow = super::layer(&solid, &wall, &point, [350.0, 450.0]);
                let restricted = PreparedGraphConjugateLinks::new_with_surface_feedback(
                    &initial,
                    vec![source_binding(&narrow, 0)],
                    1,
                    || Ok(()),
                )
                .unwrap();
                assert!(matches!(
                    owner.evaluate_conjugate_links(
                        cx,
                        &initial,
                        &source.sources,
                        &report,
                        &restricted,
                        config()
                    ),
                    Err(GraphConjugateError::Thermal(
                        DepositionError::MaterialDomain
                    ))
                ));
            }
            near(initial.temperature_k(CellId(6)).unwrap(), 400.0);
            assert_eq!(initial.time_s(), 0.0);
            assert_eq!(initial.links()[&9].conductance_w_k, 0.0);
        }
    });
}

#[test]
fn natural_and_aiding_surface_feedback_retain_buoyancy_sign_and_native_flow_envelope() {
    context(|cx| {
        let card = buoyant_card();
        let (law, caloric) = source_pair(&card);
        let law = law
            .with_density_expansion(Temperature::new(250.0), Temperature::new(450.0))
            .unwrap();
        let solid = solid(true);
        let wall = wall_caloric(&solid);
        let point = QueryPoint::new().with("T", 300.0).unwrap();
        let layer = layer(&solid, &wall, &point, [200.0, 500.0]);
        for natural in [true, false] {
            for tw in [310.0_f64, 290.0, 300.0] {
                let initial = initial(&caloric, &wall, tw);
                let graph = graph();
                let inventory = bind_inventory(cx, &graph, &initial);
                let laws = laws(2300.0);
                let owner = LiquidGraphSources::try_new(
                    cx,
                    &law,
                    &caloric,
                    &initial,
                    SourceGraphBinding {
                        graph: &graph,
                        inventory: &inventory,
                        laws: &laws,
                        acceleration_m_s2: [0.0, 0.0, -9.81],
                    },
                    caps(),
                )
                .unwrap();
                let mut b = binding(&layer);
                b.convection.property_temperature = ConvectionPropertyTemperature::FilmMean;
                b.convection.characteristic_length = Length::new(0.025);
                b.convection.law = if natural {
                    GraphConvectionLaw::Natural {
                        geometry: NaturalConvectionGeometry::VerticalPlate,
                    }
                } else {
                    GraphConvectionLaw::AidingVerticalPlate {
                        branch: 0,
                        plate_up: [0.0, 0.0, 1.0],
                        maximum_relative_density_departure: 0.1,
                    }
                };
                let prepared = PreparedGraphConjugateLinks::new_with_surface_feedback(
                    &initial,
                    vec![b],
                    1,
                    || Ok(()),
                )
                .unwrap();
                let q = if natural || tw == 300.0 {
                    0.0
                } else {
                    (tw - 300.0).signum() * 1e-5
                };
                let flows = [q, 0.0, q];
                let source = owner.resolve(cx, &initial, &flows).unwrap();
                let report = native_report(
                    cx,
                    &source.graph,
                    &inventory,
                    flows,
                    [1e5 * q, 0.0, 1e5 * q],
                    flows,
                    0.1,
                );
                let heat = owner
                    .evaluate_conjugate_links(
                        cx,
                        &initial,
                        &source.sources,
                        &report,
                        &prepared,
                        config(),
                    )
                    .unwrap();
                let sample = &heat.samples()[0];
                let ts = sample.layer.surface_temperature_k;
                let film = &sample.layer.film;
                let tf = 0.5 * (ts + 300.0);
                let rho = 1000.0 - 0.5 * (tf - 300.0);
                near(film.property_temperature_k, tf);
                near(
                    film.thermal.volumetric_expansion_per_k().unwrap(),
                    0.5 / rho,
                );
                if tw == 300.0 {
                    assert_eq!(sample.applied_heat_into_fluid_w, 0.0);
                    assert_eq!(sample.layer.conductance_w_k, 0.0);
                    assert!(film.evaluation.is_none());
                } else {
                    let evaluation = film.evaluation.as_ref().unwrap();
                    let ra = (9.81 * 0.5 / rho * (ts - 300.0) * 0.025_f64.powi(3)
                        / ((0.002 / rho) * (0.5 / (rho * 2000.0))))
                        .abs();
                    near(evaluation.nusselt.groups()["Ra"], ra);
                    if !natural {
                        near(
                            evaluation.nusselt.groups()["Re"],
                            rho * q.abs() / PIPE_AREA * 0.025 / 0.002,
                        );
                    }
                    assert!(sample.applied_heat_into_fluid_w * (tw - 300.0) > 0.0);
                    assert!(sample.layer.heat_residual_w.abs() < 2e-10);
                    // The source is evaluated at the actual surface, not at
                    // the unresisted stored wall temperature.
                    assert!((ts - 300.0).abs() < (tw - 300.0).abs());
                }
                // A nonzero natural envelope or opposing mixed stream cannot
                // be legitimized by the conducting wall's smaller heat drop.
                if tw != 300.0 {
                    let bad_q = if natural { 1e-5 } else { -q };
                    let bad_flows = [bad_q, 0.0, bad_q];
                    let bad = owner.resolve(cx, &initial, &bad_flows).unwrap();
                    let bad_report = native_report(
                        cx,
                        &bad.graph,
                        &inventory,
                        bad_flows,
                        [1e5 * bad_q, 0.0, 1e5 * bad_q],
                        bad_flows,
                        0.1,
                    );
                    assert!(
                        owner
                            .evaluate_conjugate_links(
                                cx,
                                &initial,
                                &bad.sources,
                                &bad_report,
                                &prepared,
                                config()
                            )
                            .is_err()
                    );
                }
                assert_eq!(initial.time_s(), 0.0);
                assert_eq!(initial.links()[&9].conductance_w_k, 0.0);
            }
        }
    });
}
