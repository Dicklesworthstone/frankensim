use super::super::affine::{admit_boundaries, node_mass_streams};
use super::*;

fn bind_affine(
    cx: &Cx<'_>,
    graph: &LossGraph,
    old: &DepositionSimulation<'_>,
    initial: &[f64],
    mean: &[f64],
    final_flow: &[f64],
    boundaries: Vec<MassBoundary>,
    maximum: usize,
) -> Result<GraphEnthalpyTransport, TransientTransportError> {
    GraphEnthalpyTransport::from_affine_volume_flows(
        cx,
        graph,
        initial,
        mean,
        final_flow,
        Density::new(1.0),
        &[CellId(0), CellId(1)],
        old,
        boundaries,
        maximum,
    )
}

#[test]
fn g1_zero_mean_reversal_transports_both_original_donors_and_conserves_heat() {
    with_cx(|cx| {
        let graph =
            LossGraph::new(2, vec![edge("out", 0, 1, 1.0), edge("back", 1, 0, 1.0)]).unwrap();
        for datum in [0.0, 1_099_511_627_776.0] {
            let law = Law {
                datum,
                plateau: false,
            };
            let old = thermal(&law, [1.0; 2], [0.0, 100.0], false);
            let transport = bind_affine(
                cx,
                &graph,
                &old,
                &[1.0; 2],
                &[0.0; 2],
                &[-1.0; 2],
                vec![],
                4,
            )
            .unwrap();
            assert_eq!(transport.advection.flows().len(), 4);
            // Each branch sweeps 1/4 kg in each direction over one second.
            // Two branches give the independent balances 1.5 h0-.5 h1=0,
            // 1.5 h1-.5 h0=100, even though both signed means are zero.
            let (next, audit) = transport
                .advance_advection_implicit(cx, &old, 1.0, 1e-8, 1e-8)
                .unwrap();
            near(next.cells()[&CellId(0)].specific_h_j_kg - datum, 25.0, 1e-8);
            near(next.cells()[&CellId(1)].specific_h_j_kg - datum, 75.0, 1e-8);
            near(audit.stored_energy_change_j, 0.0, 1e-8);
            assert_eq!(audit.external_heat_j, 0.0);
            assert_eq!(old.time_s(), 0.0);
        }
    });
}

fn reversing_streams(datum: f64) -> Vec<MassBoundary> {
    vec![
        MassBoundary::Inlet {
            id: 10,
            cell: CellId(0),
            kg_s: 0.5,
            material_id: MaterialId(1),
            specific_h_j_kg: datum + 100.0,
        },
        MassBoundary::Outlet {
            id: 11,
            cell: CellId(0),
            kg_s: 0.5,
        },
        MassBoundary::Inlet {
            id: 12,
            cell: CellId(1),
            kg_s: 0.5,
            material_id: MaterialId(1),
            specific_h_j_kg: datum,
        },
        MassBoundary::Outlet {
            id: 13,
            cell: CellId(1),
            kg_s: 0.5,
        },
    ]
}

#[test]
fn g1_open_reversal_retains_distinct_inlet_enthalpies_and_gross_supply() {
    with_cx(|cx| {
        let law = Law {
            datum: 0.0,
            plateau: false,
        };
        let old = thermal(&law, [1.0; 2], [0.0; 2], false);
        let graph = LossGraph::new(2, vec![edge("open", 0, 1, 1.0)]).unwrap();
        let expected =
            node_mass_streams(cx, &graph, Density::new(1.0), &[2.0], &[0.0], &[-2.0]).unwrap();
        for node in expected {
            assert_eq!(node.inlet_kg_s, 0.5);
            assert_eq!(node.outlet_kg_s, 0.5);
        }
        let transport = bind_affine(
            cx,
            &graph,
            &old,
            &[2.0],
            &[0.0],
            &[-2.0],
            reversing_streams(0.0),
            6,
        )
        .unwrap();
        let (next, audit) = transport
            .advance_advection_implicit(cx, &old, 1.0, 1e-9, 1e-9)
            .unwrap();
        // Actual endpoint donor balances: 2 h0-.5 h1=50; 2 h1-.5 h0=0.
        near(next.cells()[&CellId(0)].specific_h_j_kg, 80.0 / 3.0, 1e-9);
        near(next.cells()[&CellId(1)].specific_h_j_kg, 20.0 / 3.0, 1e-9);
        near(audit.external_heat_j, 100.0 / 3.0, 1e-9);
        near(audit.stored_energy_change_j, audit.external_heat_j, 1e-9);
        // Native explicit stepping sees the GROSS outgoing mass, not zero
        // signed mean, and retains both independently authored supply accounts.
        let (_, _, advected) = transport
            .advance_to(
                cx,
                &old,
                1.0,
                StepConfig {
                    max_dt_s: 1.0,
                    energy_tolerance_j: 1e-8,
                },
            )
            .unwrap();
        near(advected.inlet_mass_kg, 1.0, 1e-12);
        near(advected.outlet_mass_kg, 1.0, 1e-12);
        // The native gross transfer count includes all internal, inlet and
        // outlet connections: the same mass can cross more than one of them.
        near(advected.transported_mass_kg, 3.0, 1e-12);
        assert!(advected.maximum_courant > 0.0);
    });
}

#[test]
fn g3_no_reversal_parity_and_through_junction_has_no_fictitious_reservoir() {
    with_cx(|cx| {
        let law = Law {
            datum: 0.0,
            plateau: false,
        };
        let old = thermal(&law, [1.0; 2], [10.0, 40.0], false);
        let graph =
            LossGraph::new(2, vec![edge("out", 0, 1, 1.0), edge("back", 1, 0, 1.0)]).unwrap();
        let original = GraphEnthalpyTransport::from_volume_flows(
            cx,
            &graph,
            &[1.0; 2],
            Density::new(1.0),
            &[CellId(0), CellId(1)],
            &old,
            vec![],
            2,
        )
        .unwrap();
        let affine =
            bind_affine(cx, &graph, &old, &[0.0; 2], &[1.0; 2], &[2.0; 2], vec![], 2).unwrap();
        assert_eq!(original.advection.flows(), affine.advection.flows());
        let (expected, _) = original
            .advance_advection_implicit(cx, &old, 0.7, 1e-9, 1e-9)
            .unwrap();
        let (actual, _) = affine
            .advance_advection_implicit(cx, &old, 0.7, 1e-9, 1e-9)
            .unwrap();
        for (&id, cell) in expected.cells() {
            assert_eq!(
                cell.specific_h_j_kg.to_bits(),
                actual.cells()[&id].specific_h_j_kg.to_bits()
            );
        }
        let series =
            LossGraph::new(3, vec![edge("in", 0, 1, 1.0), edge("out", 1, 2, 1.0)]).unwrap();
        let expected = node_mass_streams(
            cx,
            &series,
            Density::new(1.0),
            &[2.0; 2],
            &[0.0; 2],
            &[-2.0; 2],
        )
        .unwrap();
        assert_eq!(expected[1].inlet_kg_s, 0.0);
        assert_eq!(expected[1].outlet_kg_s, 0.0);
        assert_eq!(expected[0].inlet_kg_s, 0.5);
        assert_eq!(expected[2].inlet_kg_s, 0.5);
    });
}

#[test]
fn g0_missing_gross_supply_budget_collision_and_donor_domain_refuse_without_mutation() {
    with_cx(|cx| {
        let law = Law {
            datum: 0.0,
            plateau: false,
        };
        let old = thermal(&law, [1.0; 2], [0.0; 2], false);
        let graph = LossGraph::new(2, vec![edge("open", 0, 1, 1.0)]).unwrap();
        assert!(bind_affine(cx, &graph, &old, &[2.0], &[0.0], &[-2.0], vec![], 6).is_err());
        assert!(
            bind_affine(
                cx,
                &graph,
                &old,
                &[2.0],
                &[0.0],
                &[-2.0],
                reversing_streams(0.0),
                5
            )
            .is_err()
        );
        let mut collision = reversing_streams(0.0);
        if let MassBoundary::Inlet { id, .. } = &mut collision[0] {
            *id = 1;
        }
        assert!(bind_affine(cx, &graph, &old, &[2.0], &[0.0], &[-2.0], collision, 6).is_err());
        let mut bad_donor = reversing_streams(0.0);
        if let MassBoundary::Inlet {
            specific_h_j_kg, ..
        } = &mut bad_donor[2]
        {
            *specific_h_j_kg = 2000.0;
        }
        assert!(bind_affine(cx, &graph, &old, &[2.0], &[0.0], &[-2.0], bad_donor, 6).is_err());
        assert_eq!(old.time_s(), 0.0);
        assert_eq!(old.cells()[&CellId(0)].specific_h_j_kg, 0.0);
        assert!(
            bind_affine(
                cx,
                &graph,
                &old,
                &[2.0],
                &[0.0],
                &[-2.0],
                reversing_streams(0.0),
                6
            )
            .is_ok()
        );
    });
}

#[test]
fn g0_large_through_flow_cannot_hide_small_external_reversal_or_invent_a_supply() {
    with_cx(|cx| {
        let graph = LossGraph::new(
            4,
            vec![
                edge("large-in", 0, 1, 1.0),
                edge("large-out", 1, 2, 1.0),
                edge("large-return", 2, 0, 1.0),
                edge("small-reversal", 1, 3, 1.0),
            ],
        )
        .unwrap();
        let initial = [1e16, 1e16, 1e16, 1.0];
        let mean = [1e16, 1e16, 1e16, 0.0];
        let final_flow = [1e16, 1e16, 1e16, -1.0];
        let nodes = [CellId(0), CellId(1), CellId(2), CellId(3)];
        let expected =
            node_mass_streams(cx, &graph, Density::new(1.0), &initial, &mean, &final_flow).unwrap();
        for node in [1, 3] {
            assert_eq!(expected[node].inlet_kg_s, 0.25);
            assert_eq!(expected[node].outlet_kg_s, 0.25);
        }
        for node in [0, 2] {
            assert_eq!(expected[node].inlet_kg_s, 0.0);
            assert_eq!(expected[node].outlet_kg_s, 0.0);
        }
        let mut correct = Vec::new();
        for node in [1, 3] {
            correct.push(MassBoundary::Inlet {
                id: 20 + 2 * node,
                cell: CellId(node),
                kg_s: 0.25,
                material_id: MaterialId(1),
                specific_h_j_kg: 10.0 * node as f64,
            });
            correct.push(MassBoundary::Outlet {
                id: 21 + 2 * node,
                cell: CellId(node),
                kg_s: 0.25,
            });
        }
        assert!(admit_boundaries(cx, &nodes, &expected, &correct).is_ok());
        // Both missing streams at node 1 leave NET continuity unchanged, but
        // lose their original entering heat. The large circulation is irrelevant.
        assert!(admit_boundaries(cx, &nodes, &expected, &correct[2..]).is_err());
        let mut wrong = correct.clone();
        if let MassBoundary::Inlet { kg_s, .. } = &mut wrong[0] {
            *kg_s = 0.2;
        }
        assert!(admit_boundaries(cx, &nodes, &expected, &wrong).is_err());
        let mut invented = correct;
        invented.push(MassBoundary::Inlet {
            id: 40,
            cell: CellId(0),
            kg_s: 1e-30,
            material_id: MaterialId(1),
            specific_h_j_kg: 100.0,
        });
        invented.push(MassBoundary::Outlet {
            id: 41,
            cell: CellId(0),
            kg_s: 1e-30,
        });
        assert!(admit_boundaries(cx, &nodes, &expected, &invented).is_err());
    });
}
