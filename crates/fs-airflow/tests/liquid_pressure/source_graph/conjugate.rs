//! Original material layer, original graph convection, and one native heat owner.
use super::*;
use fs_airflow::graph::thermal::transport::transient::{
    inertial::{
        InertialGraphHeat, InertialThermalControls, InertialThermalProblem, ResolvedInertialGraph,
        advance_inertial_graph_distributed,
    },
    material::convection::{
        GraphConvectionError,
        conjugate::{GraphConjugateError, GraphConjugateLinkBinding, PreparedGraphConjugateLinks},
    },
};
use fs_matdb::ClaimSelection;
use fs_material::{
    phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve},
    state_point::THERMAL_CONDUCTIVITY_DIMS,
    thermal::{
        conduction::{
            LayerConductionConfig, ReservoirLayerGeometry, source::ResolvedConductionLayer,
        },
        source::ThermalMaterialPoint,
    },
};
use fs_solver::NewtonKrylovConfig;
use std::sync::Arc;

const AREA: f64 = 0.002;
const DISTANCE: f64 = 0.05;
#[path = "conjugate/surface.rs"]
mod surface;
#[path = "conjugate/developing.rs"]
mod developing;
#[path = "conjugate/native_trace.rs"]
mod native_trace;
fn config() -> LayerConductionConfig {
    LayerConductionConfig {
        max_iterations: 96,
        absolute_heat_tolerance_w: 1e-10,
        relative_heat_tolerance: 1e-12,
    }
}
fn solid(variable: bool) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: [200.0, 350.0, 500.0]
                    .map(|t| (t, if variable { 0.002 * t } else { 0.6 }))
                    .to_vec(),
                dims: THERMAL_CONDUCTIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained().with("T", 200.0, 500.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured conjugate solid".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic-solid".into(),
            phase: "solid".into(),
            process: "numerical test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn wall_caloric(card: &MaterialCard) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_single_phase(
        card.content_hash(),
        SolidLiquidPhase::Solid,
        [200.0, 500.0]
            .map(|t| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 100.0 * (t - 200.0),
                temperature_k: t,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 1000.0,
            })
            .to_vec(),
    )
    .unwrap()
}
fn layer<'a>(
    card: &'a MaterialCard,
    caloric: &'a EquilibriumEnthalpyPhaseCurve,
    context: &'a QueryPoint,
    span: [f64; 2],
) -> ResolvedConductionLayer<'a> {
    ResolvedConductionLayer::new(
        caloric,
        ThermalMaterialPoint {
            material: card,
            context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        &PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        SolidLiquidPhase::Solid,
        span,
        || Ok(()),
    )
    .unwrap()
}
fn fluid_card() -> MaterialCard {
    let base = card(true, 1.0);
    let mut claims = ClaimSet::new();
    for (_, claim) in base.claims().claims_ordered() {
        let mut claim = claim.clone();
        if claim.key.name() == "thermal_conductivity" {
            claim.value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(250.0, 0.25), (450.0, 1.25)],
                dims: THERMAL_CONDUCTIVITY_DIMS,
            };
            claim.interpolation = InterpolationPolicy::LinearInside;
        }
        claims.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap()
}
fn initial<'a>(
    fluid: &'a ResolvedHeatCapacityCurve,
    wall: &'a EquilibriumEnthalpyPhaseCurve,
    wall_temperature: f64,
) -> DepositionSimulation<'a> {
    let base = thermal(fluid, [300.0; 6], None, false);
    let mut events: Vec<_> = base
        .cells()
        .values()
        .cloned()
        .map(|cell| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell),
        })
        .collect();
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Birth(DepositedCell {
            id: CellId(6),
            segment_id: 6,
            body_id: 1,
            material_id: MaterialId(2),
            mass_kg: 1000.0 * AREA * 2.0 * DISTANCE,
            specific_h_j_kg: 100.0 * (wall_temperature - 200.0),
            geometry: CellGeometry {
                start_m: [0.0; 3],
                end_m: [2.0 * DISTANCE, 0.0, 0.0],
                volume_m3: AREA * 2.0 * DISTANCE,
                surface_area_m2: 0.1,
            },
            boundary: Boundary {
                ambient_k: 300.0,
                htc_w_m2_k: 0.0,
                reservoir: None,
                source_w: 0.0,
            },
        }),
    });
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 9,
            a: CellId(2),
            b: CellId(6),
            area_m2: AREA,
            conductance_w_k: 0.0,
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
fn binding<'a, 'm>(layer: &'a ResolvedConductionLayer<'m>) -> GraphConjugateLinkBinding<'a, 'm> {
    GraphConjugateLinkBinding {
        convection: GraphConvectiveLink {
            link: 9,
            fluid_cell: CellId(2),
            characteristic_length: Length::new(0.03),
            property_temperature: ConvectionPropertyTemperature::Reservoir,
            law: GraphConvectionLaw::Forced {
                branch: 0,
                correlation: CorrelationId::CircularDuctLaminarCwt,
                settings: developed_settings(),
            },
        },
        wall_layer: layer,
        geometry: ReservoirLayerGeometry::planar(AREA, DISTANCE).unwrap(),
    }
}

#[test]
fn actual_surface_closes_constant_and_variable_solid_flux_in_both_directions_and_equilibrium() {
    context(|cx| {
        let fluid_card = fluid_card();
        let (law, caloric) = source_pair(&fluid_card);
        for variable in [false, true] {
            let solid = solid(variable);
            let wall = wall_caloric(&solid);
            let point = QueryPoint::new().with("T", 300.0).unwrap();
            let layer = layer(&solid, &wall, &point, [200.0, 500.0]);
            for tw in [400.0, 280.0, 300.0] {
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
                let prepared = PreparedGraphConjugateLinks::new(
                    &initial,
                    vec![binding(&layer)],
                    16,
                    || Ok(()),
                )
                .unwrap();
                let heat = prepared
                    .resolve(
                        &initial,
                        config(),
                        |candidate, links| {
                            owner.evaluate_links(cx, candidate, &source.sources, &report, links, 16)
                        },
                        || Ok(()),
                    )
                    .unwrap();
                let sample = &heat.samples()[0];
                let film = &heat.fluid_samples()[0];
                let ts = sample.layer.material_face_temperature_k;
                let g_film = AREA * 3.66 * 0.5 / 0.03;
                near(film.property_temperature_k, 300.0);
                near(film.conductance_w_k, g_film);
                assert!((tw.min(300.0)..=tw.max(300.0)).contains(&ts));
                assert!(Arc::ptr_eq(&sample.source_line, layer.source_line()));
                assert_eq!(sample.source_line.knots.len(), 3);
                let solid_power = if variable {
                    AREA / DISTANCE * 0.001 * (tw - ts) * (tw + ts)
                } else {
                    AREA / DISTANCE * 0.6 * (tw - ts)
                };
                assert!((sample.applied_heat_into_fluid_w - solid_power).abs() < 2e-10);
                assert!((sample.applied_heat_into_fluid_w - g_film * (ts - 300.0)).abs() < 2e-10);
                assert!(sample.layer.heat_residual_w.abs() < 2e-10);
                if !variable {
                    let g = 1.0 / (1.0 / g_film + DISTANCE / (AREA * 0.6));
                    near(sample.conductance_w_k, g);
                    near(sample.applied_heat_into_fluid_w, g * (tw - 300.0));
                }
                if tw == 300.0 {
                    assert_eq!(sample.applied_heat_into_fluid_w, 0.0);
                    near(ts, 300.0);
                } else {
                    assert!(sample.applied_heat_into_fluid_w * (tw - 300.0) > 0.0);
                }
                let mut coefficients = BTreeMap::new();
                heat.apply_to(&mut coefficients, || Ok(())).unwrap();
                let before = coefficients.clone();
                assert!(heat.apply_to(&mut coefficients, || Ok(())).is_err());
                assert_eq!(coefficients, before);
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
            }
        }
    });
}

#[test]
fn foreign_wall_area_film_stale_callback_source_domain_and_cancellation_refuse_atomically() {
    context(|cx| {
        let fluid_card = fluid_card();
        let (law, caloric) = source_pair(&fluid_card);
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
        let mut wrong_area = binding(&layer);
        wrong_area.geometry = ReservoirLayerGeometry::planar(2.0 * AREA, DISTANCE).unwrap();
        assert!(
            PreparedGraphConjugateLinks::new(&initial, vec![wrong_area], 16, || Ok(())).is_err()
        );
        for film in [false, true] {
            let mut unsupported = binding(&layer);
            if film {
                unsupported.convection.property_temperature =
                    ConvectionPropertyTemperature::FilmMean;
            } else {
                unsupported.convection.law = GraphConvectionLaw::ForcedSphere { branch: 0 };
            }
            assert!(
                PreparedGraphConjugateLinks::new(&initial, vec![unsupported], 16, || Ok(()))
                    .is_err()
            );
        }
        let foreign_wall = wall_caloric(&solid);
        let foreign_layer = self::layer(&solid, &foreign_wall, &point, [200.0, 500.0]);
        assert!(
            PreparedGraphConjugateLinks::new(
                &initial,
                vec![binding(&foreign_layer)],
                16,
                || Ok(())
            )
            .is_err()
        );
        let prepared =
            PreparedGraphConjugateLinks::new(&initial, vec![binding(&layer)], 16, || Ok(()))
                .unwrap();
        let other = self::initial(&caloric, &wall, 390.0);
        assert!(
            prepared
                .resolve(
                    &other,
                    config(),
                    |_, links| owner.evaluate_links(
                        cx,
                        &initial,
                        &source.sources,
                        &report,
                        links,
                        16
                    ),
                    || Ok(())
                )
                .is_err()
        );
        assert!(
            prepared
                .resolve(
                    &initial,
                    config(),
                    |candidate, links| {
                        let mut changed = links.to_vec();
                        changed[0].property_temperature = ConvectionPropertyTemperature::FilmMean;
                        owner.evaluate_links(cx, candidate, &source.sources, &report, &changed, 16)
                    },
                    || Ok(())
                )
                .is_err()
        );
        let narrow = self::layer(&solid, &wall, &point, [350.0, 450.0]);
        let constrained =
            PreparedGraphConjugateLinks::new(&initial, vec![binding(&narrow)], 16, || Ok(()))
                .unwrap();
        assert!(
            constrained
                .resolve(
                    &initial,
                    config(),
                    |candidate, links| owner.evaluate_links(
                        cx,
                        candidate,
                        &source.sources,
                        &report,
                        links,
                        16
                    ),
                    || Ok(())
                )
                .is_err()
        );
        let called = std::cell::Cell::new(false);
        assert!(
            prepared
                .resolve(
                    &initial,
                    config(),
                    |candidate, links| {
                        called.set(true);
                        owner.evaluate_links(cx, candidate, &source.sources, &report, links, 16)
                    },
                    || Err(DepositionError::Cancelled)
                )
                .is_err()
        );
        assert!(!called.get());
        assert_eq!(initial.time_s(), 0.0);
        near(initial.temperature_k(CellId(6)).unwrap(), 400.0);
        assert_eq!(initial.links()[&9].conductance_w_k, 0.0);
    });
}

#[test]
fn native_coupled_momentum_advection_and_surface_film_use_one_endpoint_heat() {
    context(|cx| {
        let fluid_card = fluid_card();
        let (law, caloric) = source_pair(&fluid_card);
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
        let mut surface_binding = binding(&layer);
        surface_binding.convection.property_temperature = ConvectionPropertyTemperature::FilmMean;
        let prepared = PreparedGraphConjugateLinks::new_with_surface_feedback(
            &initial,
            vec![surface_binding],
            16,
            || Ok(()),
        )
        .unwrap();
        let initial_flows = [1e-5, 0.0, 1e-5];
        let first = owner.resolve(cx, &initial, &initial_flows).unwrap();
        let held_flow = native_report(
            cx,
            &first.graph,
            &inventory,
            initial_flows,
            [1.0, 0.0, 1.0],
            initial_flows,
            1.0,
        );
        let held = owner
            .evaluate_conjugate_links(
                cx,
                &initial,
                &first.sources,
                &held_flow,
                &prepared,
                config(),
            )
            .unwrap();
        let branches: Vec<_> = inventory
            .inertial_branches()
            .iter()
            .map(|b| b.unwrap())
            .collect();
        let boundaries = vec![FixedPressure {
            node: 0,
            pressure: Pressure::new(0.0),
        }];
        let momentum = InertialGraph::try_new(
            cx,
            &first.graph,
            Density::new(1000.0),
            &branches,
            &boundaries,
            &initial_flows,
            0.0,
            32,
        )
        .unwrap();
        let problem = InertialThermalProblem {
            momentum: &momentum,
            node_cells: inventory.node_cells(),
            maximum_connections: 32,
        };
        let controls = InertialThermalControls {
            solver: NewtonKrylovConfig {
                absolute_tolerance: 1e-9,
                relative_tolerance: 0.0,
                ..NewtonKrylovConfig::default()
            },
            maximum_cells: 16,
            maximum_newton_iterations: 32,
            maximum_endpoint_evaluations: 2048,
            energy_scale_j: 1.0,
            flow_scale_m3_s: 1e-5,
            finite_difference_step: 1e-4,
            cell_tolerance_j: 1e-6,
            total_thermal_tolerance_j: 1e-6,
            impulse_tolerance_pa_s: 1e-6,
            mechanical_tolerance_j: 1e-8,
            average_flow_tolerance_m3_s: 1e-12,
            momentum: InertialGraphControls {
                maximum_drag_step: 2.0,
                maximum_neglected_power_w: 100.0,
                maximum_neglected_energy_j: 100.0,
            },
        };
        let (next_momentum, next, report) = advance_inertial_graph_distributed(
            cx,
            &problem,
            &inventory,
            &initial,
            1.0,
            controls,
            &[0.0; 7],
            &initial_flows,
            |cx, candidate, flows| {
                let resolved = owner.resolve(cx, candidate, flows).map_err(|e| {
                    GraphConjugateError::Convection(GraphConvectionError::Source(e))
                })?;
                Ok::<_, GraphConjugateError>(ResolvedInertialGraph {
                    graph: resolved.graph,
                    pressure_boundaries: boundaries.clone(),
                    driving_pressure_pa: vec![1.0, 0.0, 1.0],
                    payload: resolved.sources,
                })
            },
            |cx, candidate, sources, hydraulics| {
                let heat = owner.evaluate_conjugate_links(
                    cx,
                    candidate,
                    sources,
                    hydraulics,
                    &prepared,
                    config(),
                )?;
                let mut conductances = BTreeMap::new();
                heat.apply_to(&mut conductances, || {
                    cx.checkpoint().map_err(|_| DepositionError::Cancelled)
                })?;
                Ok::<_, GraphConjugateError>(InertialGraphHeat {
                    mass_boundaries: vec![],
                    heat: StepCoefficients {
                        sources_w: BTreeMap::new(),
                        conductances_w_k: conductances,
                        payload: heat,
                    }
                    .into(),
                })
            },
        )
        .unwrap();
        near(next.time_s(), 1.0);
        near(next_momentum.accepted_time_s(), 1.0);
        let sample = &report.heat.samples()[0];
        let film = &sample.layer.film;
        let tf = next.temperature_k(CellId(2)).unwrap();
        let tw = next.temperature_k(CellId(6)).unwrap();
        assert!(tf > 300.0 && tw < 400.0);
        assert!(sample.layer.surface_temperature_k > tf && sample.layer.surface_temperature_k < tw);
        near(sample.fluid_temperature_k, tf);
        near(sample.wall_temperature_k, tw);
        near(
            film.property_temperature_k,
            0.5 * (sample.layer.surface_temperature_k + tf),
        );
        assert!(
            (sample.layer.conductance_w_k - held.samples()[0].layer.conductance_w_k).abs() > 1e-9
        );
        assert!((film.conductance_w_k - held.samples()[0].layer.film.conductance_w_k).abs() > 1e-9);
        let wall_change = next.cells()[&CellId(6)].mass_kg
            * (next.cells()[&CellId(6)].specific_h_j_kg
                - initial.cells()[&CellId(6)].specific_h_j_kg);
        assert!((wall_change + sample.applied_heat_into_fluid_w).abs() < 1e-6);
        assert!(report.maximum_cell_residual_j <= controls.cell_tolerance_j);
        assert!(report.hydraulics.maximum_impulse_residual_pa_s <= controls.impulse_tolerance_pa_s);
        assert!(report.thermal.residual_j.abs() <= controls.total_thermal_tolerance_j);
        assert!(report.thermal.stored_energy_change_j.abs() < 1e-6);
        near(report.thermal.external_heat_j, 0.0);
        assert!(Arc::ptr_eq(&sample.source_line, layer.source_line()));
        assert_eq!(initial.time_s(), 0.0);
        assert_eq!(initial.links()[&9].conductance_w_k, 0.0);
    });
}
