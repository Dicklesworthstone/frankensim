//! Axial Fourier transport and advection share the original finite-cell energy.
use super::*;
use fs_airflow::graph::thermal::transport::transient::material::conduction::{
    GraphConductionError, GraphConductiveLink,
};

fn conductivity_card(scale: f64) -> MaterialCard {
    let base = card(false, 1.0);
    let mut claims = ClaimSet::new();
    for (_, claim) in base.claims().claims_ordered() {
        let mut claim = claim.clone();
        if claim.key.name() == "thermal_conductivity" {
            claim.value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(250.0, 0.25 * scale), (450.0, 1.25 * scale)],
                dims: Dims([1, 1, -3, -1, 0, 0]),
            };
            claim.interpolation = InterpolationPolicy::LinearInside;
        }
        claims.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap()
}

fn link(id: u64, a: u64, b: u64, area: f64) -> ThermalLink {
    ThermalLink {
        id,
        a: CellId(a),
        b: CellId(b),
        area_m2: area,
        conductance_w_k: 0.0,
    }
}

fn linked<'a>(
    caloric: &'a ResolvedHeatCapacityCurve,
    temperatures: [f64; 6],
    stepped: bool,
    links: &[ThermalLink],
) -> DepositionSimulation<'a> {
    let initial = thermal(caloric, temperatures, None, false);
    let mut schedule = Vec::new();
    for mut cell in initial.cells().values().cloned() {
        if stepped && cell.id == CellId(2) {
            cell.geometry.end_m[2] = 0.25;
            cell.geometry.volume_m3 = PIPE_AREA * 0.25;
            cell.mass_kg = 1000.0 * cell.geometry.volume_m3;
        } else if stepped && cell.id == CellId(3) {
            cell.geometry.start_m[2] = 0.25;
            cell.geometry.volume_m3 = 2.0 * PIPE_AREA * 0.75;
            cell.mass_kg = 1000.0 * cell.geometry.volume_m3;
        }
        schedule.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell),
        });
    }
    schedule.extend(links.iter().cloned().map(|link| TimedEvent {
        time_s: 0.0,
        event: Event::Contact(link),
    }));
    DepositionSimulation::new(
        schedule,
        &[(MaterialId(1), caloric.curve() as &dyn ThermalLaw)],
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

#[test]
fn original_local_k_and_actual_half_areas_control_stagnant_and_reversed_pipe_heat() {
    context(|cx| {
        for scale in [1.0, 2.0] {
            let card = conductivity_card(scale);
            let (law, caloric) = source_pair(&card);
            for reversed in [false, true] {
                let [a, b] = if reversed { [3, 2] } else { [2, 3] };
                let initial = linked(
                    &caloric,
                    [300.0, 300.0, 300.0, 400.0, 300.0, 300.0],
                    true,
                    &[link(20, a, b, PIPE_AREA)],
                );
                let graph = graph();
                let inventory = bind_inventory(cx, &graph, &initial);
                let mut laws = laws(1e6);
                laws[0][0].resistance = MaterialFlowResistance::try_new(
                    A * 0.25,
                    0.0,
                    Area::new(PIPE_AREA),
                    Length::new(0.03),
                    (0.0, 1e6),
                )
                .unwrap();
                laws[0][1].resistance = MaterialFlowResistance::try_new(
                    A * 0.75,
                    0.0,
                    Area::new(2.0 * PIPE_AREA),
                    Length::new(0.03),
                    (0.0, 1e6),
                )
                .unwrap();
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
                for flow in [0.0, -1e-6, 1e-6] {
                    let resolved = owner.resolve(cx, &initial, &[flow, 0.0, flow]).unwrap();
                    let heat = owner
                        .evaluate_conduction_links(
                            cx,
                            &initial,
                            &resolved.sources,
                            &[GraphConductiveLink {
                                link: 20,
                                branch: 0,
                            }],
                            16,
                        )
                        .unwrap();
                    let sample = &heat.samples()[0];
                    let expected =
                        scale / (0.125 / (0.5 * PIPE_AREA) + 0.375 / (1.0 * 2.0 * PIPE_AREA));
                    near(sample.conductance_w_k, expected);
                    near(
                        sample.heat_a_to_b_w,
                        if reversed {
                            100.0 * expected
                        } else {
                            -100.0 * expected
                        },
                    );
                    assert_eq!(sample.cells, [CellId(a), CellId(b)]);
                    let index = if reversed { 1 } else { 0 };
                    near(
                        sample.transport[index].thermal_conductivity_w_m_k(),
                        0.5 * scale,
                    );
                    near(
                        sample.transport[1 - index].thermal_conductivity_w_m_k(),
                        scale,
                    );
                    near(sample.half_resistances_k_w[index], 250.0 / scale);
                    near(sample.half_resistances_k_w[1 - index], 187.5 / scale);
                    assert!(std::ptr::eq(
                        &heat.source_sections()[0][0],
                        &resolved.sources.sections()[0][0]
                    ));
                    assert_eq!(
                        heat.source_sections()[0][1]
                            .momentum
                            .medium
                            .material()
                            .card_identity(),
                        law.reference_state().material().card_identity()
                    );
                }
            }
        }
    });
}

#[test]
fn refined_original_cells_converge_to_the_variable_conductivity_fourier_flux() {
    context(|cx| {
        let card = conductivity_card(1.0);
        let (law, caloric) = source_pair(&card);
        let graph = graph();
        let mut errors = Vec::new();
        for n in [8, 16, 32] {
            let template = thermal(&caloric, [300.0; 6], None, false);
            let dx = 1.0 / n as f64;
            let mut cells: Vec<_> = [CellId(0), CellId(1)]
                .map(|id| template.cells()[&id].clone())
                .into_iter()
                .collect();
            let path: Vec<_> = (0..n).map(|i| CellId(2 + i as u64)).collect();
            for (i, &id) in path.iter().enumerate() {
                let mut cell = template.cells()[&CellId(2)].clone();
                cell.id = id;
                cell.segment_id = id.0;
                cell.geometry.start_m[2] = i as f64 * dx;
                cell.geometry.end_m[2] = (i + 1) as f64 * dx;
                cell.geometry.volume_m3 = PIPE_AREA * dx;
                cell.mass_kg = 1000.0 * cell.geometry.volume_m3;
                // k(T)=0.005*T-1; integrating k*dT/dx=-75 gives this
                // independent steady profile between 400 K and 300 K.
                let x = (i as f64 + 0.5) * dx;
                let t = ((1.0 - 0.75 * x).sqrt() + 1.0) / 0.005;
                cell.specific_h_j_kg = 2000.0 * (t - 250.0);
                cells.push(cell);
            }
            for (old, id) in [
                (CellId(4), CellId(n as u64 + 2)),
                (CellId(5), CellId(n as u64 + 3)),
            ] {
                let mut cell = template.cells()[&old].clone();
                cell.id = id;
                cell.segment_id = id.0;
                cells.push(cell);
            }
            let mut schedule: Vec<_> = cells
                .into_iter()
                .map(|cell| TimedEvent {
                    time_s: 0.0,
                    event: Event::Birth(cell),
                })
                .collect();
            let bindings: Vec<_> = (0..n - 1)
                .map(|i| {
                    let id = 100 + i as u64;
                    schedule.push(TimedEvent {
                        time_s: 0.0,
                        event: Event::Contact(link(id, path[i].0, path[i + 1].0, PIPE_AREA)),
                    });
                    GraphConductiveLink {
                        link: id,
                        branch: 0,
                    }
                })
                .collect();
            let initial = DepositionSimulation::new(
                schedule,
                &[(MaterialId(1), caloric.curve() as &dyn ThermalLaw)],
                DepositionBudget {
                    max_cells: 64,
                    max_links: 64,
                    max_events: 128,
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
            .0;
            let inventory = DistributedGraphInventory::new(
                cx,
                &graph,
                Density::new(1000.0),
                &[CellId(0), CellId(1)],
                &[path, vec![CellId(n as u64 + 2)], vec![CellId(n as u64 + 3)]],
                &initial,
                64,
            )
            .unwrap();
            let mut laws = laws(1e6);
            let mut section = laws[0][0].clone();
            section.resistance = MaterialFlowResistance::try_new(
                A * dx,
                0.0,
                Area::new(PIPE_AREA),
                Length::new(0.03),
                (0.0, 1e6),
            )
            .unwrap();
            laws[0] = vec![section; n];
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
                SourceGraphLimits {
                    maximum_cells: 64,
                    ..caps()
                },
            )
            .unwrap();
            let resolved = owner.resolve(cx, &initial, &[0.0; 3]).unwrap();
            let conduction = owner
                .evaluate_conduction_links(cx, &initial, &resolved.sources, &bindings, 64)
                .unwrap();
            let expected = 75.0 * PIPE_AREA;
            let error = conduction
                .samples()
                .iter()
                .map(|sample| (sample.heat_a_to_b_w - expected).abs() / expected)
                .fold(0.0_f64, f64::max);
            errors.push(error);
        }
        assert!(errors[1] < errors[0] / 2.8);
        assert!(errors[2] < errors[1] / 2.8);
        assert!(errors[2] < 5e-4);
    });
}

#[test]
fn stagnant_original_cells_equilibrate_through_one_native_backward_euler_link() {
    context(|cx| {
        let card = card(false, 1.0);
        let (law, caloric) = source_pair(&card);
        let initial = linked(
            &caloric,
            [300.0, 300.0, 400.0, 300.0, 300.0, 300.0],
            false,
            &[link(20, 2, 3, PIPE_AREA)],
        );
        let graph = graph();
        let inventory = bind_inventory(cx, &graph, &initial);
        let laws = laws(1e6);
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
        let dt = 100_000.0;
        let conductance = 0.5 * PIPE_AREA / 0.5;
        let capacity = 0.5 * 2000.0;
        // Independent two-body backward-Euler balance, C*dT = +/-dt*G*(Ta-Tb).
        let exchanged = dt * conductance * 100.0 / (1.0 + 2.0 * dt * conductance / capacity);
        let step = BackwardEulerStep::new(&initial, dt, || Ok(())).unwrap();
        let trial = step
            .evaluate(
                &[0.0, 0.0, -exchanged, exchanged, 0.0, 0.0],
                |candidate| {
                    let resolved = owner.resolve(cx, candidate, &[0.0; 3]).unwrap();
                    let conduction = owner
                        .evaluate_conduction_links(
                            cx,
                            candidate,
                            &resolved.sources,
                            &[GraphConductiveLink {
                                link: 20,
                                branch: 0,
                            }],
                            16,
                        )
                        .unwrap();
                    let mut map = BTreeMap::new();
                    conduction.apply_to(cx, &mut map).unwrap();
                    Ok::<_, DepositionError>(StepCoefficients {
                        sources_w: BTreeMap::new(),
                        conductances_w_k: map,
                        payload: (resolved.sources, conduction),
                    })
                },
                || Ok(()),
            )
            .unwrap();
        let (next, (sources, conduction), audit) = trial.accept(1e-7, 1e-7, || Ok(())).unwrap();
        near(
            next.temperature_k(CellId(2)).unwrap(),
            400.0 - exchanged / capacity,
        );
        near(
            next.temperature_k(CellId(3)).unwrap(),
            300.0 + exchanged / capacity,
        );
        near(dt * conduction.samples()[0].heat_a_to_b_w, exchanged);
        assert!(std::ptr::eq(
            &sources.sections()[0][0],
            &conduction.source_sections()[0][0]
        ));
        assert!(audit.residual_j.abs() < 1e-7);
        near(audit.external_heat_j, 0.0);
        near(initial.time_s(), 0.0);
        near(initial.temperature_k(CellId(2)).unwrap(), 400.0);
        let old = owner.resolve(cx, &initial, &[0.0; 3]).unwrap();
        assert!(matches!(
            owner.evaluate_conduction_links(
                cx,
                &next,
                &old.sources,
                &[GraphConductiveLink {
                    link: 20,
                    branch: 0
                }],
                16
            ),
            Err(GraphConductionError::Source(SourceGraphError::Invalid(_)))
        ));
        let mut occupied = BTreeMap::from([(20, 0.123), (99, 0.456)]);
        let before = occupied.clone();
        assert!(conduction.apply_to(cx, &mut occupied).is_err());
        assert_eq!(occupied, before);
    });
}

#[test]
fn native_directed_advection_and_reciprocal_axial_conduction_add_once_through_reversal() {
    context(|cx| {
        let card = card(false, 1.0);
        let (law, caloric) = source_pair(&card);
        let initial = linked(
            &caloric,
            [300.0, 300.0, 400.0, 300.0, 300.0, 300.0],
            false,
            &[link(20, 2, 3, PIPE_AREA)],
        );
        let graph = graph();
        let inventory = bind_inventory(cx, &graph, &initial);
        let laws = laws(1e6);
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
        let dt = 0.1;
        for q in [-1e-9_f64, 0.0, 1e-9] {
            let step = BackwardEulerStep::new(&initial, dt, || Ok(())).unwrap();
            let trial = step
                .evaluate_with_transport(
                    &[0.0; 6],
                    |candidate| {
                        let resolved = owner.resolve(cx, candidate, &[q, 0.0, q]).unwrap();
                        let transport = inventory
                            .bind(cx, &resolved.graph, &[q, 0.0, q], candidate, vec![], 32)
                            .unwrap();
                        let mut coefficients = transport.endpoint_transport(cx, candidate).unwrap();
                        let conduction = owner
                            .evaluate_conduction_links(
                                cx,
                                candidate,
                                &resolved.sources,
                                &[GraphConductiveLink {
                                    link: 20,
                                    branch: 0,
                                }],
                                16,
                            )
                            .unwrap();
                        conduction
                            .apply_to(cx, &mut coefficients.heat.bulk.conductances_w_k)
                            .unwrap();
                        Ok::<_, DepositionError>(coefficients)
                    },
                    || Ok(()),
                )
                .unwrap();
            let advective = 1000.0 * q.abs() * 2000.0 * 100.0;
            let conductive = 0.1;
            near(trial.residual_j()[2], dt * (advective + conductive));
            near(
                trial.residual_j()[3],
                -dt * (conductive + if q > 0.0 { advective } else { 0.0 }),
            );
            near(
                trial.residual_j()[0],
                if q < 0.0 { -dt * advective } else { 0.0 },
            );
            assert!(trial.residual_j().iter().sum::<f64>().abs() < 1e-12);
            near(initial.time_s(), 0.0);
        }
    });
}

#[test]
fn partial_apertures_junctions_and_duplicate_whole_faces_cannot_borrow_axial_geometry() {
    context(|cx| {
        let card = card(false, 1.0);
        let (law, caloric) = source_pair(&card);
        let graph = graph();
        let laws = laws(1e6);
        for contacts in [
            vec![link(20, 2, 3, 0.5 * PIPE_AREA)],
            vec![link(20, 2, 3, 2.0 * PIPE_AREA)],
            vec![link(20, 0, 2, PIPE_AREA)],
            vec![link(20, 2, 4, PIPE_AREA)],
            vec![link(20, 2, 3, PIPE_AREA), link(21, 3, 2, PIPE_AREA)],
        ] {
            let initial = linked(&caloric, [300.0; 6], false, &contacts);
            let inventory = bind_inventory(cx, &graph, &initial);
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
            let resolved = owner.resolve(cx, &initial, &[0.0; 3]).unwrap();
            assert!(matches!(
                owner.evaluate_conduction_links(
                    cx,
                    &initial,
                    &resolved.sources,
                    &[GraphConductiveLink {
                        link: 20,
                        branch: 0
                    }],
                    16
                ),
                Err(GraphConductionError::Link { .. })
            ));
            assert_eq!(initial.time_s(), 0.0);
            assert!(
                initial
                    .links()
                    .values()
                    .all(|link| link.conductance_w_k == 0.0)
            );
        }
    });
}
