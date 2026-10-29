//! The existing integrated material flux must reach the actual circuit history.
use super::*;
use fs_material::phase::SolidLiquidPhase;
use fs_material::thermal::conduction::{
    LayerConductionConfig,
    network::{LayerLinkBinding, PreparedLayerConduction},
    source::ResolvedConductionLayer,
};

fn layer_controls() -> LayerConductionConfig {
    LayerConductionConfig {
        max_iterations: 80,
        absolute_heat_tolerance_w: 1e-10,
        relative_heat_tolerance: 1e-12,
    }
}

fn with_layers<R>(
    scale: f64,
    run: impl FnOnce(
        &Cx<'_>,
        &CancelGate,
        &MaterialHeatNetwork<'_>,
        &PreparedLayerConduction<'_, '_>,
        DepositionSimulation<'_>,
        MaterialCircuit,
    ) -> R,
) -> R {
    let cards = [source(0.2, 1.0), source(0.6, 2.0)].map(|original| {
        let mut claims = ClaimSet::new();
        for (_, claim) in original.claims().claims_ordered() {
            let mut claim = claim.clone();
            if claim.key.name() == "thermal_conductivity" {
                let PropertyValue::Curve { knots, .. } = &mut claim.value else {
                    unreachable!()
                };
                for (_, value) in knots {
                    *value *= scale;
                }
            }
            claims.insert_claim(claim).unwrap();
        }
        MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap()
    });
    let curves = [curve(&cards[0], 100.0), curve(&cards[1], 100.0)];
    let context = context();
    let ids = [CellId(30), CellId(10)];
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = [0, 1].map(|i| MaterialCell {
        cell: ids[i],
        material_id: MaterialId(i as u64),
        curve: &curves[i],
        point: binding(&cards[i], &context),
        conductivity: Some((&key, MaterialPropertySelection::SingleClaimOnly)),
    });
    let resistors = [0, 1].map(|i| CellResistor {
        destination: ResistorCell {
            branch: i + 1,
            cell: ids[i],
        },
        selection: MaterialPropertySelection::SingleClaimOnly,
        length_m: 1.0,
        area_m2: 1.0,
    });
    let links = [ConductionLink {
        link: 7,
        distance_a_m: 0.5,
        distance_b_m: 0.5,
    }];
    let bulk = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &links,
    };
    let bases = [
        (1, ElectricalPropertyBasis::Resistivity),
        (2, ElectricalPropertyBasis::Resistivity),
    ];
    let electrical = ElectricalContactNetwork {
        bulk: &bulk,
        contacts: &[],
        electrical_bases: &bases,
    };
    let assembly = MaterialHeatNetwork {
        electrical: &electrical,
        thermal_contacts: &[],
        radiation: None,
    };
    let mut events = [0, 1]
        .map(|i| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: ids[i],
                segment_id: i as u64,
                body_id: i as u64,
                material_id: MaterialId(i as u64),
                mass_kg: 1.0,
                specific_h_j_kg: if i == 0 { -80.0 } else { -20.0 },
                geometry: CellGeometry {
                    start_m: [i as f64, 0.0, 0.0],
                    end_m: [i as f64 + 1.0, 0.0, 0.0],
                    volume_m3: 1.0,
                    surface_area_m2: 2.5,
                },
                boundary: Boundary {
                    ambient_k: 350.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        })
        .to_vec();
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 7,
            a: ids[0],
            b: ids[1],
            area_m2: 0.1,
            conductance_w_k: 999.0,
        }),
    });
    let thermal = DepositionSimulation::new(
        events,
        &[
            (MaterialId(0), &curves[0] as &dyn ThermalLaw),
            (MaterialId(1), &curves[1] as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.001,
            energy_tolerance_j: 1e-7,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    let resolved = electrical.resolve(&thermal, || false).unwrap();
    let circuit = MaterialCircuit::new(
        CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, Branch::Capacitor { farads: 0.01 }),
                (1, 2, Branch::Resistor { ohms: 0.0 }),
                (2, 0, Branch::Resistor { ohms: 0.0 }),
            ],
            transformers: vec![],
        },
        resolved.bulk.resistors,
        vec![],
        vec![0.5],
        &[],
    )
    .unwrap();
    let layers = [0, 1].map(|i| {
        ResolvedConductionLayer::new(
            &curves[i],
            binding(&cards[i], &context),
            &key,
            ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
            SolidLiquidPhase::Solid,
            [300.0, 400.0],
            || Ok(()),
        )
        .unwrap()
    });
    let integrated = PreparedLayerConduction::new(
        &thermal,
        vec![LayerLinkBinding {
            link: 7,
            layers: [&layers[0], &layers[1]],
            distances_m: [0.5; 2],
            contact_area_resistance_m2_k_w: 0.0,
        }],
        1,
        || Ok(()),
    )
    .unwrap();
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 73,
                kernel_id: 52,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        run(&cx, &gate, &assembly, &integrated, thermal, circuit)
    })
}

#[test]
fn integrated_source_curve_changes_native_heat_and_retains_the_joint_circuit_balance() {
    with_layers(1.0, |cx, _, assembly, layers, mut thermal, mut circuit| {
        let prepared = assembly
            .prepare(cx, &thermal)
            .unwrap()
            .with_integrated_conduction(cx, &thermal, layers, layer_controls())
            .unwrap();
        let initial = prepared.resolve(cx, &thermal).unwrap();
        // For identical linear k(T) layers, integrating the actual temperature
        // path gives the endpoint arithmetic mean, independently of the split.
        near(initial.bulk.conductances_w_k[&7], 0.1 * (1.1 + 1.4) / 2.0);
        assert!((initial.bulk.conductances_w_k[&7] - 0.1 / (0.5 / 1.1 + 0.5 / 1.4)).abs() > 1e-3);
        let repeated = prepared.resolve(cx, &thermal).unwrap();
        assert_eq!(
            initial.bulk.conductances_w_k,
            repeated.bulk.conductances_w_k
        );
        let total = thermal.energy_j() + 12.5;
        for step in 1..=20 {
            let candidate = prepared
                .step_candidate(
                    cx,
                    &thermal,
                    step as f64 * 0.0005,
                    &circuit,
                    &[],
                    controls(),
                    &[0.0; 2],
                )
                .unwrap();
            let source = &candidate.report().source.material;
            let path = &source.integrated_conduction[0];
            let ta = candidate.thermal().temperature_k(CellId(30)).unwrap();
            let tb = candidate.thermal().temperature_k(CellId(10)).unwrap();
            let g = 0.1 * (1.0 + ((ta + tb) * 0.5 - 300.0) / 200.0);
            near(path.response.conductance_w_k, g);
            near(path.applied_power_a_to_b_w, g * (ta - tb));
            near(source.material.links[0].conductance_w_k, g);
            near(
                candidate.thermal().energy_j() + candidate.report().electrical.stored_after_j,
                total,
            );
            let (next_thermal, next_circuit, _) = candidate.into_parts();
            thermal = next_thermal;
            circuit = next_circuit;
        }
    });
}

#[test]
fn data_only_layer_conductivity_change_reaches_temperature_and_electrical_history() {
    let run = |scale| {
        with_layers(
            scale,
            |cx, _, assembly, layers, mut thermal, mut circuit| {
                let prepared = assembly
                    .prepare(cx, &thermal)
                    .unwrap()
                    .with_integrated_conduction(cx, &thermal, layers, layer_controls())
                    .unwrap();
                for step in 1..=20 {
                    let candidate = prepared
                        .step_candidate(
                            cx,
                            &thermal,
                            step as f64 * 0.0005,
                            &circuit,
                            &[],
                            controls(),
                            &[0.0; 2],
                        )
                        .unwrap();
                    let (t, c, _) = candidate.into_parts();
                    thermal = t;
                    circuit = c;
                }
                [
                    thermal.temperature_k(CellId(30)).unwrap(),
                    thermal.temperature_k(CellId(10)).unwrap(),
                    circuit.charges_c()[0],
                ]
            },
        )
    };
    let original = run(1.0);
    let changed = run(4.0);
    assert!(changed[0] > original[0] && changed[1] < original[1]);
    assert!((changed[2] - original[2]).abs() > 1e-8);
}

#[test]
fn layer_geometry_context_contact_ownership_and_cancellation_refuse() {
    with_layers(1.0, |cx, gate, assembly, layers, thermal, circuit| {
        let original = layers.bindings().next().unwrap();
        let mismatched = PreparedLayerConduction::new(
            &thermal,
            vec![LayerLinkBinding {
                link: 7,
                layers: original.layers,
                distances_m: [0.25, 0.75],
                contact_area_resistance_m2_k_w: 0.0,
            }],
            1,
            || Ok(()),
        )
        .unwrap();
        assert!(
            assembly
                .prepare(cx, &thermal)
                .unwrap()
                .with_integrated_conduction(cx, &thermal, &mismatched, layer_controls())
                .is_err()
        );

        let cell = &assembly.electrical.bulk.cells[0];
        let context = cell
            .point
            .context
            .clone()
            .with("pressure", 105000.0)
            .unwrap();
        let other = ResolvedConductionLayer::new(
            cell.curve,
            fs_material::thermal::source::ThermalMaterialPoint {
                context: &context,
                ..cell.point
            },
            cell.conductivity.as_ref().unwrap().0,
            ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
            SolidLiquidPhase::Solid,
            [300.0, 400.0],
            || Ok(()),
        )
        .unwrap();
        let other_context = PreparedLayerConduction::new(
            &thermal,
            vec![LayerLinkBinding {
                link: 7,
                layers: [&other, original.layers[1]],
                distances_m: [0.5; 2],
                contact_area_resistance_m2_k_w: 0.0,
            }],
            1,
            || Ok(()),
        )
        .unwrap();
        assert!(
            assembly
                .prepare(cx, &thermal)
                .unwrap()
                .with_integrated_conduction(cx, &thermal, &other_context, layer_controls())
                .is_err()
        );

        let cells = assembly.electrical.bulk.cells;
        let cards = [
            cells[0].point.material.clone(),
            cells[1].point.material.clone(),
        ];
        let card = interface(&cards, 0.2);
        let contact_context = QueryPoint::new()
            .with("aT", 320.0)
            .unwrap()
            .with("bT", 380.0)
            .unwrap()
            .with("aM", 0.0)
            .unwrap()
            .with("bM", 0.0)
            .unwrap()
            .with("pressure", 100000.0)
            .unwrap();
        let contacts = [MaterialContactLink {
            link: 7,
            point: ThermalInterfacePoint {
                interface: &card,
                materials: [cells[0].point.material, cells[1].point.material],
                context: &contact_context,
                temperature_axes: ["aT", "bT"],
                liquid_fraction_axes: [Some("aM"), Some("bM")],
                mean_temperature_axis: None,
            },
            sides: [0, 1].map(|i| ContactCell {
                cell: cells[i].cell,
                material_id: cells[i].material_id,
                curve: cells[i].curve,
            }),
            selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        }];
        let duplicate = MaterialHeatNetwork {
            electrical: assembly.electrical,
            thermal_contacts: &contacts,
            radiation: None,
        };
        assert!(
            duplicate
                .prepare(cx, &thermal)
                .unwrap()
                .with_integrated_conduction(cx, &thermal, layers, layer_controls())
                .is_err()
        );

        let prepared = assembly
            .prepare(cx, &thermal)
            .unwrap()
            .with_integrated_conduction(cx, &thermal, layers, layer_controls())
            .unwrap();
        assert!(
            prepared
                .extend_fixed_circuit(cx, assembly, &thermal, &[])
                .is_err()
        );
        let accepted = thermal.cells().clone();
        gate.request();
        assert!(matches!(
            prepared.step_candidate(cx, &thermal, 0.0005, &circuit, &[], controls(), &[0.0; 2]),
            Err(circuit_enthalpy::CircuitEnthalpyError::Cancelled)
        ));
        assert_eq!(thermal.cells(), &accepted);
        assert_eq!(circuit.charges_c(), &[0.5]);
    });
}
