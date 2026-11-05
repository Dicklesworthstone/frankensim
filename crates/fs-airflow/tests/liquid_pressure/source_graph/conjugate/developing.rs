//! Original duct-section evidence, material surface balance and reciprocal heat.
use super::*;
use fs_airflow::graph::thermal::transport::transient::material::convection::developing::{
    DevelopingDuct, DevelopingWall, GraphDevelopingConvectionError, GraphDevelopingSection,
    conjugate::{GraphDevelopingConjugateLinkBinding, PreparedGraphDevelopingConjugateLinks},
};

fn dimensions(rectangular: bool) -> (DevelopingDuct, f64, f64) {
    if rectangular {
        let width = (2.0 * PIPE_AREA).sqrt();
        let height = width / 2.0;
        let perimeter = 2.0 * (width + height);
        (
            DevelopingDuct::Rectangular {
                width: Length::new(width),
                height: Length::new(height),
            },
            4.0 * PIPE_AREA / perimeter,
            perimeter * 0.5,
        )
    } else {
        let diameter = (4.0 * PIPE_AREA / core::f64::consts::PI).sqrt();
        (
            DevelopingDuct::Circular {
                diameter: Length::new(diameter),
            },
            diameter,
            core::f64::consts::PI * diameter * 0.5,
        )
    }
}

fn fixture<'a>(
    fluid: &'a ResolvedHeatCapacityCurve,
    wall: &'a EquilibriumEnthalpyPhaseCurve,
    temperature: f64,
    area: f64,
) -> DepositionSimulation<'a> {
    let base = initial(fluid, wall, temperature);
    let mut events: Vec<_> = base
        .cells()
        .values()
        .cloned()
        .map(|mut cell| {
            if cell.id == CellId(6) {
                cell.mass_kg = 1000.0 * area * 2.0 * DISTANCE;
                cell.geometry.volume_m3 = area * 2.0 * DISTANCE;
                cell.geometry.surface_area_m2 = 2.0 * area + 0.1;
            }
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(cell),
            }
        })
        .collect();
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            area_m2: area,
            ..base.links()[&9].clone()
        }),
    });
    DepositionSimulation::new(
        events,
        &[
            (MaterialId(1), fluid.curve() as &dyn ThermalLaw),
            (MaterialId(2), wall as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 16,
            max_links: 16,
            max_events: 32,
            max_substeps: 16,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-8,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

fn binding<'a, 'm>(
    layer: &'a ResolvedConductionLayer<'m>,
    rectangular: bool,
    property_temperature: ConvectionPropertyTemperature,
) -> GraphDevelopingConjugateLinkBinding<'a, 'm> {
    let (duct, _, area) = dimensions(rectangular);
    GraphDevelopingConjugateLinkBinding {
        section: GraphDevelopingSection {
            fluid_cell: CellId(2),
            wall: DevelopingWall::Link(9),
            duct,
            forward_start: Length::new(0.0),
            reverse_start: Length::new(0.5),
            property_temperature,
        },
        wall_layer: layer,
        geometry: ReservoirLayerGeometry::planar(area, DISTANCE).unwrap(),
    }
}

// Independent cumulative Hausen equation. No production section or root
// implementation participates in the expected heat or surface temperature.
fn circular_h(rate: f64, conductivity: f64) -> f64 {
    let (_, diameter, _) = dimensions(false);
    let start = if rate > 0.0 { 0.0 } else { 0.5 };
    let mean = |length: f64| {
        let gz = 1000.0 * rate.abs() / PIPE_AREA * 2000.0 * diameter * diameter
            / (conductivity * length);
        3.66 + 0.0668 * gz / (1.0 + 0.04 * gz.powf(2.0 / 3.0))
    };
    ((start + 0.5) * mean(start + 0.5)
        - if start == 0.0 {
            0.0
        } else {
            start * mean(start)
        })
        / 0.5
        * conductivity
        / diameter
}

fn affine_surface(cell: f64, bulk: f64, h: f64) -> f64 {
    let a = 0.001 / DISTANCE; // original wall k(T)=0.002*T
    let rhs = a * cell * cell + h * bulk;
    2.0 * rhs / (h + (h * h + 4.0 * a * rhs).sqrt())
}

#[test]
fn g1_developing_wall_integral_matches_analytic_surface_and_reciprocal_native_heat() {
    context(|cx| {
        let card = card(false, 1.0);
        let (law, caloric) = source_pair(&card);
        for variable in [false, true] {
            let solid = solid(variable);
            let wall = wall_caloric(&solid);
            let point = QueryPoint::new().with("T", 300.0).unwrap();
            let layer = layer(&solid, &wall, &point, [200.0, 500.0]);
            let (_, _, area) = dimensions(false);
            for tw in [400.0, 280.0, 300.0] {
                let initial = fixture(&caloric, &wall, tw, area);
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
                let prepared = PreparedGraphDevelopingConjugateLinks::new(
                    &initial,
                    vec![binding(
                        &layer,
                        false,
                        ConvectionPropertyTemperature::Reservoir,
                    )],
                    1,
                    || Ok(()),
                )
                .unwrap();
                for sign in [-1.0, 1.0] {
                    let q = sign * 1e-5;
                    let flows = [q, 0.0, q];
                    let source = owner.resolve(cx, &initial, &flows).unwrap();
                    let report = native_report(
                        cx,
                        &source.graph,
                        &inventory,
                        flows,
                        [sign, 0.0, sign],
                        flows,
                        0.1,
                    );
                    let heat = owner
                        .evaluate_developing_conjugate_links(
                            cx,
                            &initial,
                            &source.sources,
                            &report,
                            &prepared,
                            config(),
                        )
                        .unwrap();
                    let row = &heat.samples()[0];
                    let h = circular_h(q, 0.5);
                    let ts = if tw == 300.0 {
                        300.0
                    } else if variable {
                        affine_surface(tw, 300.0, h)
                    } else {
                        (0.6 / DISTANCE * tw + h * 300.0) / (0.6 / DISTANCE + h)
                    };
                    let power = area * h * (ts - 300.0);
                    near(row.layer.surface_temperature_k, ts);
                    near(row.layer.film.wall_temperature_k, ts);
                    near(row.wall_temperature_k, tw);
                    near(row.layer.film.property_temperature_k, 300.0);
                    near(row.layer.film.conductance_w_k, area * h);
                    near(row.applied_heat_into_fluid_w, power);
                    assert_eq!(row.layer.film.section.upstream().is_some(), sign < 0.0);
                    assert!(row.layer.heat_residual_w.abs() <= 2e-10);
                    assert!(Arc::ptr_eq(&row.source_line, layer.source_line()));
                    assert_eq!(
                        row.layer.film.material.material().card_identity(),
                        card.content_hash()
                    );
                    let mut coefficients = BTreeMap::new();
                    heat.apply_to(&mut coefficients, || Ok(())).unwrap();
                    let saved = coefficients.clone();
                    assert!(heat.apply_to(&mut coefficients, || Ok(())).is_err());
                    assert_eq!(coefficients, saved);
                    let mut empty = BTreeMap::new();
                    assert!(
                        heat.apply_to(&mut empty, || Err(DepositionError::Cancelled))
                            .is_err()
                    );
                    assert!(empty.is_empty());
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
                    near(trial.residual_j()[2], -0.1 * power);
                    near(trial.residual_j()[6], 0.1 * power);
                    near(trial.residual_j().iter().sum(), 0.0);
                    if variable && tw == 400.0 && sign > 0.0 {
                        let source = owner.resolve(cx, &initial, &[q, q, -q]).unwrap();
                        let flow = source
                            .graph
                            .solve(
                                &[
                                    FixedPressure {
                                        node: 0,
                                        pressure: Pressure::new(sign),
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
                            .evaluate_quasi_steady_developing_conjugate_links(
                                cx,
                                &initial,
                                &source.sources,
                                &flow,
                                0.1,
                                &prepared,
                                config(),
                            )
                            .unwrap();
                        near(steady.samples()[0].applied_heat_into_fluid_w, power);
                    }
                }
            }
        }
    });
}

#[test]
fn g1_g3_actual_surface_source_retains_rectangular_table_and_circular_film_feedback() {
    context(|cx| {
        for rectangular in [false, true] {
            let card = if rectangular {
                card(false, 0.09)
            } else {
                fluid_card()
            };
            let (law, caloric) = source_pair(&card);
            let solid = solid(true);
            let wall = wall_caloric(&solid);
            let point = QueryPoint::new().with("T", 300.0).unwrap();
            let layer = layer(&solid, &wall, &point, [200.0, 500.0]);
            let (_, diameter, area) = dimensions(rectangular);
            let q = if rectangular {
                40.0 * 0.5 * PIPE_AREA / (1000.0 * 2000.0 * diameter * diameter)
            } else {
                1e-5
            };
            for sign in [-1.0, 1.0] {
                let initial = fixture(&caloric, &wall, 400.0, area);
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
                let prepared = PreparedGraphDevelopingConjugateLinks::new(
                    &initial,
                    vec![binding(
                        &layer,
                        rectangular,
                        ConvectionPropertyTemperature::FilmMean,
                    )],
                    1,
                    || Ok(()),
                )
                .unwrap();
                let flows = [sign * q, 0.0, sign * q];
                let source = owner.resolve(cx, &initial, &flows).unwrap();
                let pressure = A * if rectangular { 0.00018 } else { 0.002 } * sign * q;
                let report = native_report(
                    cx,
                    &source.graph,
                    &inventory,
                    flows,
                    [pressure, 0.0, pressure],
                    flows,
                    0.1,
                );
                let heat = owner
                    .evaluate_developing_conjugate_links(
                        cx,
                        &initial,
                        &source.sources,
                        &report,
                        &prepared,
                        config(),
                    )
                    .unwrap();
                let row = &heat.samples()[0];
                let ts = row.layer.surface_temperature_k;
                let film = &row.layer.film;
                near(film.property_temperature_k, 0.5 * (ts + 300.0));
                near(film.thermal.temperature_k(), film.property_temperature_k);
                if rectangular {
                    // Original tabulated cumulative means: Gz=80 ->6.88,
                    // Gz=40 ->5.61; downstream section Nu=2*5.61-6.88.
                    let nu = if sign > 0.0 { 6.88 } else { 4.34 };
                    near(film.section.nusselt(), nu);
                    near(ts, affine_surface(400.0, 300.0, nu * 0.5 / diameter));
                    near(film.section.downstream().groups()["Pr"], 0.72);
                } else {
                    let mut lo = 300.0;
                    let mut hi = 400.0;
                    for _ in 0..80 {
                        let t = 0.5 * (lo + hi);
                        let k = 0.005 * (0.5 * (t + 300.0)) - 1.0;
                        let residual = 0.001 / DISTANCE * (400.0 - t) * (400.0 + t)
                            - circular_h(sign * q, k) * (t - 300.0);
                        if residual > 0.0 {
                            lo = t;
                        } else {
                            hi = t;
                        }
                    }
                    near(ts, 0.5 * (lo + hi));
                    near(
                        film.conductance_w_k,
                        area * circular_h(sign * q, 0.005 * film.property_temperature_k - 1.0),
                    );
                }
                near(
                    row.applied_heat_into_fluid_w,
                    area / DISTANCE * 0.001 * (400.0 - ts) * (400.0 + ts),
                );
                assert!(row.layer.heat_residual_w.abs() <= 2e-10);
                assert_eq!(
                    film.material.material().card_identity(),
                    card.content_hash()
                );
            }
        }
    });
}

#[test]
fn g0_g4_developing_surface_budget_geometry_and_crossing_fail_without_advancing_heat() {
    context(|cx| {
        let card = fluid_card();
        let (law, caloric) = source_pair(&card);
        let solid = solid(true);
        let wall = wall_caloric(&solid);
        let point = QueryPoint::new().with("T", 300.0).unwrap();
        let layer = layer(&solid, &wall, &point, [200.0, 500.0]);
        let initial = fixture(&caloric, &wall, 400.0, dimensions(false).2);
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
        let make = || binding(&layer, false, ConvectionPropertyTemperature::FilmMean);
        let prepared =
            PreparedGraphDevelopingConjugateLinks::new(&initial, vec![make()], 1, || Ok(()))
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
        assert!(matches!(
            owner.evaluate_developing_conjugate_links(
                cx,
                &initial,
                &source.sources,
                &report,
                &prepared,
                LayerConductionConfig {
                    max_iterations: 1,
                    absolute_heat_tolerance_w: 0.0,
                    relative_heat_tolerance: 0.0
                }
            ),
            Err(GraphConjugateError::SurfaceBalance { link: 9, .. })
        ));
        let mut wrong = make();
        wrong.section.duct = DevelopingDuct::Circular {
            diameter: Length::new(1.01 * dimensions(false).1),
        };
        let wrong = PreparedGraphDevelopingConjugateLinks::new(&initial, vec![wrong], 1, || Ok(()))
            .unwrap();
        assert!(matches!(
            owner.evaluate_developing_conjugate_links(
                cx,
                &initial,
                &source.sources,
                &report,
                &wrong,
                config()
            ),
            Err(GraphConjugateError::Developing(_))
        ));
        let mut reservoir = make();
        reservoir.section.wall = DevelopingWall::Reservoir;
        assert!(matches!(
            PreparedGraphDevelopingConjugateLinks::new(&initial, vec![reservoir], 1, || Ok(())),
            Err(GraphConjugateError::Invalid(_))
        ));
        assert!(
            PreparedGraphDevelopingConjugateLinks::new(
                &initial,
                vec![make(), make()],
                2,
                || Ok(())
            )
            .is_err()
        );
        assert!(matches!(
            PreparedGraphDevelopingConjugateLinks::new(&initial, vec![make()], 1, || Err(
                DepositionError::Cancelled
            )),
            Err(GraphConjugateError::Thermal(DepositionError::Cancelled))
        ));
        let average = [5e-6, 0.0, 5e-6];
        let source = owner.resolve(cx, &initial, &average).unwrap();
        let crossing = native_report(
            cx,
            &source.graph,
            &inventory,
            [-1e-5, 0.0, -1e-5],
            [300.5, 0.0, 300.5],
            average,
            0.1,
        );
        assert!(matches!(
            owner.evaluate_developing_conjugate_links(
                cx,
                &initial,
                &source.sources,
                &crossing,
                &prepared,
                config()
            ),
            Err(GraphConjugateError::Developing(
                GraphDevelopingConvectionError::Convection(_)
            ))
        ));
        assert_eq!(initial.time_s(), 0.0);
        assert_eq!(initial.links()[&9].conductance_w_k, 0.0);
        near(initial.temperature_k(CellId(6)).unwrap(), 400.0);
    });
}
