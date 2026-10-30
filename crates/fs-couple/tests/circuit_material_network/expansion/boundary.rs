//! Actual source-expanded areas feed the native boundary heat/entropy owner.
use super::*;
use fs_material::thermal::transport::expansion::surface::PreparedExpansionSurface;

fn timeline(
    caloric: &ResolvedHeatCapacityCurve,
    boundary: Boundary,
    mut later: Vec<TimedEvent>,
) -> DepositionSimulation<'_> {
    let mut events = vec![TimedEvent {
        time_s: 0.0,
        event: Event::Birth(DepositedCell {
            id: CELL,
            segment_id: 0,
            body_id: 0,
            material_id: MaterialId(0),
            mass_kg: 1.0,
            specific_h_j_kg: 10.0 * CP,
            geometry: CellGeometry {
                start_m: [0.0; 3],
                end_m: [1.0, 0.0, 0.0],
                volume_m3: 1.0,
                surface_area_m2: 6.0,
            },
            boundary,
        }),
    }];
    events.append(&mut later);
    DepositionSimulation::new(
        events,
        &[(MaterialId(0), caloric.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 4,
            max_substeps: 32,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-7,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

fn faces(
    path: &PreparedFreeExpansion,
    branch: usize,
    exposed: [f64; 3],
    reservoir: Option<[f64; 3]>,
) -> ExpandingBoundary {
    ExpandingBoundary {
        branch,
        exposed: PreparedExpansionSurface::try_new(path.clone(), exposed).unwrap(),
        reservoir: reservoir
            .map(|areas| PreparedExpansionSurface::try_new(path.clone(), areas).unwrap()),
    }
}
fn area(reference: [f64; 3], alpha: [f64; 3], t: f64) -> f64 {
    let s = alpha.map(|a| 1.0 + a * (t - 300.0));
    (0..3)
        .map(|i| reference[i] * s[(i + 1) % 3] * s[(i + 2) % 3])
        .sum()
}
fn driven(sources: Vec<fs_couple::material_circuit::MaterialResistor>) -> MaterialCircuit {
    MaterialCircuit::new(
        CircuitGraph {
            node_count: 2,
            branches: vec![
                (1, 0, Branch::Resistor { ohms: f64::NAN }),
                (1, 0, Branch::VoltageSource { port: 0 }),
            ],
            transformers: vec![],
        },
        sources,
        vec![],
        vec![],
        &[4.0],
    )
    .unwrap()
}
fn environment() -> Boundary {
    Boundary {
        ambient_k: 295.0,
        htc_w_m2_k: 0.2,
        reservoir: Some(Reservoir {
            area_m2: 3.0,
            conductance_w_k: 0.4,
            temperature_k: 320.0,
        }),
        source_w: 0.0,
    }
}

#[test]
fn g1_source_area_and_current_share_the_endpoint_heat_and_reservoir_entropy() {
    let alpha = [0.006, 0.003, 0.009];
    let (card, caloric, keys) = material_expansion(ElectricalPropertyBasis::Conductivity, alpha);
    let context = temperature(300.0);
    let cells = [cell(&card, &caloric, &context)];
    let resistors = [resistor(0)];
    let network = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &[],
    };
    let mut heat = timeline(&caloric, environment(), vec![]);
    let reference = heat.cells()[&CELL].clone();
    let source = path(&card, &keys, 380.0);
    let expanded = network
        .prepare_with_bases(&heat, &[(0, ElectricalPropertyBasis::Conductivity)], || {
            false
        })
        .unwrap()
        .with_free_expansion_and_boundaries(
            &heat,
            &[ConductorFreeExpansion {
                branch: 0,
                source: source.clone(),
                longitudinal_axis: 0,
            }],
            &[faces(&source, 0, [2.0, 1.0, 0.0], Some([0.0, 1.0, 2.0]))],
            || false,
        )
        .unwrap();
    assert!(
        expanded.resolve(&heat, || false).is_err(),
        "surface heat cannot be silently omitted"
    );
    let mut circuit = driven(
        expanded
            .resolve_with_boundary_conductances(&heat, || false)
            .unwrap()
            .bulk
            .resistors,
    );
    for step in 1..=3 {
        let old_t = heat.temperature_k(CELL).unwrap();
        let dt = 0.2;
        let ga = |t| 0.2 * area([2.0, 1.0, 0.0], alpha, t);
        let gr = |t| 0.4 * area([0.0, 1.0, 2.0], alpha, t) / 3.0;
        let expected = root(old_t, 379.0, |t| {
            CP * (t - old_t)
                - dt * (16.0 / resistance(alpha, 0, t) + ga(t) * (295.0 - t) + gr(t) * (320.0 - t))
        });
        let (candidate, entropy) = expanded
            .step_candidate_with_entropy(
                &heat,
                f64::from(step) * dt,
                &circuit,
                &[4.0],
                config(),
                &[0.0],
                || false,
            )
            .unwrap();
        let t = candidate.thermal().temperature_k(CELL).unwrap();
        assert!((t - expected).abs() < 1e-7);
        let report = candidate.report();
        let boundary = &report.source.boundaries[0];
        near(
            boundary.exposed.total_area_m2(),
            area([2.0, 1.0, 0.0], alpha, t),
        );
        near(
            boundary.reservoir.unwrap().total_area_m2(),
            area([0.0, 1.0, 2.0], alpha, t),
        );
        near(boundary.ambient_conductance_w_k, ga(t));
        near(boundary.reservoir_conductance_w_k, gr(t));
        near(boundary.ambient_heat_into_cell_w, ga(t) * (295.0 - t));
        near(boundary.reservoir_heat_into_cell_w, gr(t) * (320.0 - t));
        let qa = dt * ga(t) * (295.0 - t);
        let qr = dt * gr(t) * (320.0 - t);
        near(report.non_electrical_heat_j, qa + qr);
        near(
            report.thermal.stored_energy_change_j,
            report.deposited_joule_j + qa + qr,
        );
        near(entropy.reservoir_entropy_input_j_k, qa / 295.0 + qr / 320.0);
        near(
            entropy.conduction_generation_j_k,
            qa * (1.0 / t - 1.0 / 295.0) + qr * (1.0 / t - 1.0 / 320.0),
        );
        near(entropy.thermal_storage_change_j_k, CP * (t / old_t).ln());
        near(
            entropy.electrical_loss_generation_j_k,
            report.deposited_joule_j / t,
        );
        near(entropy.unclassified_source_at_cell_j_k, 0.0);
        assert!(!entropy.has_unclassified_sources);
        near(entropy.accounting_residual_j_k, 0.0);
        assert_eq!(
            candidate.thermal().cells()[&CELL].geometry,
            reference.geometry
        );
        assert_eq!(
            candidate.thermal().cells()[&CELL].boundary,
            reference.boundary
        );
        assert_eq!(
            candidate.thermal().cells()[&CELL].mass_kg,
            reference.mass_kg
        );
        let (next_heat, next_circuit, _) = candidate.into_parts();
        heat = next_heat;
        circuit = next_circuit;
    }
}

#[test]
fn g3_zero_strain_matches_original_boundary_and_native_events_change_actual_heat() {
    let (card, caloric, keys) = material_expansion(ElectricalPropertyBasis::Resistivity, [0.0; 3]);
    let context = temperature(300.0);
    let cells = [cell(&card, &caloric, &context)];
    let resistors = [resistor(0)];
    let network = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &[],
    };
    let changed = Boundary {
        ambient_k: 330.0,
        htc_w_m2_k: 0.6,
        reservoir: Some(Reservoir {
            temperature_k: 290.0,
            conductance_w_k: 0.7,
            ..environment().reservoir.unwrap()
        }),
        source_w: 0.0,
    };
    let heat = timeline(
        &caloric,
        environment(),
        vec![TimedEvent {
            time_s: 0.2,
            event: Event::Boundary(CELL, changed),
        }],
    );
    let source = path(&card, &keys, 380.0);
    let expanded = network
        .prepare(&heat, || false)
        .unwrap()
        .with_free_expansion_and_boundaries(
            &heat,
            &[ConductorFreeExpansion {
                branch: 0,
                source: source.clone(),
                longitudinal_axis: 0,
            }],
            &[faces(&source, 0, [2.0, 1.0, 0.0], Some([0.0, 1.0, 2.0]))],
            || false,
        )
        .unwrap();
    let circuit = driven(
        network
            .prepare(&heat, || false)
            .unwrap()
            .resolve(&heat, || false)
            .unwrap()
            .resistors,
    );
    let fixed = network
        .step_candidate(&heat, 0.2, &circuit, &[4.0], config(), &[0.0], || false)
        .unwrap();
    let current = expanded
        .step_candidate(&heat, 0.2, &circuit, &[4.0], config(), &[0.0], || false)
        .unwrap();
    assert_eq!(current.thermal().cells(), fixed.thermal().cells());
    assert_eq!(current.circuit(), fixed.circuit());
    assert_eq!(current.report().thermal, fixed.report().thermal);
    assert_eq!(
        current.thermal().cells()[&CELL].boundary,
        environment(),
        "endpoint event remains pending"
    );
    let (left, circuit, _) = current.into_parts();
    let (updated, _) = left
        .advance_to(
            0.2,
            StepConfig {
                max_dt_s: 0.1,
                energy_tolerance_j: 1e-7,
            },
            || Ok(()),
        )
        .unwrap();
    assert_eq!(updated.cells()[&CELL].boundary, changed);
    let next = expanded
        .step_candidate(&updated, 0.3, &circuit, &[4.0], config(), &[0.0], || false)
        .unwrap();
    let sample = &next.report().source.boundaries[0];
    near(
        sample.ambient_heat_into_cell_w,
        1.8 * (330.0 - sample.temperature_k),
    );
    near(
        sample.reservoir_heat_into_cell_w,
        0.7 * (290.0 - sample.temperature_k),
    );
    near(
        next.report().non_electrical_heat_j,
        0.1 * (sample.ambient_heat_into_cell_w + sample.reservoir_heat_into_cell_w),
    );
}

#[test]
fn g3_surface_direction_data_changes_cooling_current_and_magnetic_motion() {
    let run = |normal| {
        let (card, caloric, keys) =
            material_expansion(ElectricalPropertyBasis::Resistivity, [0.04, 0.0, 0.0]);
        let context = temperature(300.0);
        let cells = [cell(&card, &caloric, &context)];
        let resistors = [resistor(1)];
        let network = MaterialNetwork {
            cells: &cells,
            resistors: &resistors,
            links: &[],
        };
        let mut heat = timeline(
            &caloric,
            Boundary {
                htc_w_m2_k: 1.0,
                reservoir: Some(Reservoir {
                    area_m2: 2.0,
                    conductance_w_k: 0.0,
                    temperature_k: 300.0,
                }),
                ..boundary()
            },
            vec![],
        );
        let source = path(&card, &keys, 380.0);
        // The same cube's complete reference surface is [2,2,2]. Only which
        // face pair is insulated changes; both partitions remain disjoint.
        let covered = if normal == 0 {
            [0.0, 2.0, 0.0]
        } else {
            [2.0, 0.0, 0.0]
        };
        let areas = std::array::from_fn(|i| 2.0 - covered[i]);
        let expanded = network
            .prepare(&heat, || false)
            .unwrap()
            .with_free_expansion_and_boundaries(
                &heat,
                &[ConductorFreeExpansion {
                    branch: 1,
                    source: source.clone(),
                    longitudinal_axis: 0,
                }],
                &[faces(&source, 1, areas, Some(covered))],
                || false,
            )
            .unwrap();
        let magnetic = scalar_card(
            MAGNETIC_PERMEABILITY_PROPERTY,
            MAGNETIC_PERMEABILITY_DIMS,
            1.0,
        );
        let coil = MaterialStorage {
            branch: 0,
            element: MaterialStorageElement::Inductor(
                ResolvedInductor::resolve(
                    &magnetic,
                    &QueryPoint::new(),
                    MaterialPropertySelection::SingleClaimOnly,
                    1.0,
                    1.0,
                    1,
                )
                .unwrap(),
            ),
        };
        let circuit = MaterialCircuit::new_with_storage(
            CircuitGraph {
                node_count: 2,
                branches: vec![
                    (1, 0, Branch::Inductor { henries: f64::NAN }),
                    (1, 0, Branch::Resistor { ohms: f64::NAN }),
                ],
                transformers: vec![],
            },
            expanded
                .resolve_with_boundary_conductances(&heat, || false)
                .unwrap()
                .bulk
                .resistors,
            vec![coil],
            vec![4.0],
            vec![],
            &[],
        )
        .unwrap();
        let mut moving = MaterialFieldMotion::new(circuit, 0, body(), 0.0, [0.5, 1.5]).unwrap();
        let stored = |motion: &MaterialFieldMotion| {
            motion
                .inductor()
                .unwrap()
                .stored_energy_j(motion.circuit().fluxes_wb()[0])
                .unwrap()
                + 0.5 * motion.body().inertia().mass() * motion.velocity_m_s().powi(2)
        };
        for step in 1..=6 {
            let energy = heat.energy_j() + stored(&moving);
            let next = expanded
                .field_motion_candidate(
                    &heat,
                    f64::from(step) * 0.1,
                    &moving,
                    &[],
                    0.0,
                    config(),
                    &[0.0],
                    || false,
                )
                .unwrap();
            near(
                next.thermal().energy_j() + stored(next.motion()) - energy,
                next.report().non_electrical_heat_j,
            );
            near(
                next.report().non_electrical_heat_j,
                0.1 * next.report().source.boundaries[0].ambient_heat_into_cell_w,
            );
            let (t, m, _) = next.into_parts();
            heat = t;
            moving = m;
        }
        (
            heat.temperature_k(CELL).unwrap(),
            moving.circuit().fluxes_wb()[0],
            moving.gap_m(),
        )
    };
    let unchanged_area = run(0);
    let growing_area = run(2);
    assert!(
        growing_area.0 < unchanged_area.0
            && growing_area.1 > unchanged_area.1
            && growing_area.2 < unchanged_area.2,
        "only reference face direction changed: greater area cools the resistor and strengthens actual field motion"
    );
}

#[test]
fn g4_surface_source_coverage_changes_and_late_cancellation_preserve_histories() {
    let (card, caloric, keys) = material_expansion(ElectricalPropertyBasis::Resistivity, [0.01; 3]);
    let context = temperature(300.0);
    let cells = [cell(&card, &caloric, &context)];
    let resistors = [resistor(0)];
    let network = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &[],
    };
    let changed = Boundary {
        reservoir: Some(Reservoir {
            area_m2: 2.0,
            ..environment().reservoir.unwrap()
        }),
        ..environment()
    };
    let heat = timeline(
        &caloric,
        environment(),
        vec![TimedEvent {
            time_s: 0.2,
            event: Event::Boundary(CELL, changed),
        }],
    );
    let source = path(&card, &keys, 380.0);
    let expansion = [ConductorFreeExpansion {
        branch: 0,
        source: source.clone(),
        longitudinal_axis: 0,
    }];
    let source_faces = faces(&source, 0, [2.0, 1.0, 0.0], Some([0.0, 1.0, 2.0]));
    assert!(
        network
            .prepare(&heat, || false)
            .unwrap()
            .with_free_expansion_and_boundaries(
                &heat,
                &expansion,
                &[faces(&source, 0, [3.0, 1.0, 0.0], Some([0.0, 1.0, 2.0]))],
                || false
            )
            .is_err()
    );
    let (foreign, _, foreign_keys) =
        material_expansion(ElectricalPropertyBasis::Resistivity, [0.011; 3]);
    let foreign_path = path(&foreign, &foreign_keys, 380.0);
    assert!(
        network
            .prepare(&heat, || false)
            .unwrap()
            .with_free_expansion_and_boundaries(
                &heat,
                &expansion,
                &[faces(
                    &foreign_path,
                    0,
                    [2.0, 1.0, 0.0],
                    Some([0.0, 1.0, 2.0])
                )],
                || false
            )
            .is_err()
    );
    let expanded = network
        .prepare(&heat, || false)
        .unwrap()
        .with_free_expansion_and_boundaries(&heat, &expansion, &[source_faces], || false)
        .unwrap();
    let circuit = driven(
        expanded
            .resolve_with_boundary_conductances(&heat, || false)
            .unwrap()
            .bulk
            .resistors,
    );
    let saved = heat.cells().clone();
    let saved_circuit = circuit.clone();
    let mut polls = 0;
    let expected = expanded
        .step_candidate_with_entropy(&heat, 0.2, &circuit, &[4.0], config(), &[0.0], || {
            polls += 1;
            false
        })
        .unwrap()
        .0;
    let mut count = 0;
    assert!(matches!(
        expanded.step_candidate_with_entropy(
            &heat,
            0.2,
            &circuit,
            &[4.0],
            config(),
            &[0.0],
            || {
                count += 1;
                count == polls
            }
        ),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    assert_eq!(heat.cells(), &saved);
    assert_eq!(circuit, saved_circuit);
    let replay = expanded
        .step_candidate(&heat, 0.2, &circuit, &[4.0], config(), &[0.0], || false)
        .unwrap();
    assert_eq!(replay.thermal().cells(), expected.thermal().cells());
    assert_eq!(replay.circuit(), expected.circuit());
    let (left, circuit, _) = expected.into_parts();
    let (updated, _) = left
        .advance_to(
            0.2,
            StepConfig {
                max_dt_s: 0.1,
                energy_tolerance_j: 1e-7,
            },
            || Ok(()),
        )
        .unwrap();
    let saved = updated.cells().clone();
    let saved_circuit = circuit.clone();
    assert!(
        expanded
            .step_candidate(&updated, 0.3, &circuit, &[4.0], config(), &[0.0], || false)
            .is_err()
    );
    assert_eq!(updated.cells(), &saved);
    assert_eq!(circuit, saved_circuit);
}
