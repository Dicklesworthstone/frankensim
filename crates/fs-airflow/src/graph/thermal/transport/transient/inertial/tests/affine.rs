use super::*;

#[test]
fn g1_original_inertial_reversal_advects_finite_heat_at_zero_signed_mean() {
    with_cx(|cx| {
        let graph = network(2, &[(0, 1), (1, 0)], &[1.0; 2]);
        let momentum = owner(cx, &graph, &[1.0; 2], &[1.0; 2]);
        let law = Caloric;
        let old = inventory(&law, &[300.0, 400.0], &[0.0; 2]);
        let nodes = [CellId(0), CellId(1)];
        let problem = InertialThermalProblem {
            momentum: &momentum,
            node_cells: &nodes,
            maximum_connections: 4,
        };
        let (next_momentum, next, report) = advance_inertial_graph_affine(
            cx,
            &problem,
            None,
            &old,
            1.0,
            controls(),
            &[25_000.0, -25_000.0],
            &[0.0; 2],
            |_, _, _| Ok::<_, &'static str>(response(&graph, &[-2.0; 2])),
            |_, _, _, _| Ok(heat(())),
        )
        .unwrap();
        // Independent momentum: q1-q0 + (q1+q0)/2 = -2 gives q1=-1.
        // Each of two branches sweeps 1/4 kg each way, giving the two-cell
        // endpoint thermal solution (325 K, 375 K), not the old (300,400).
        for flow in &report.hydraulics.average_volume_flows_m3_s {
            near(*flow, 0.0, 2e-12);
        }
        for flow in next_momentum.volume_flows_m3_s() {
            near(flow, -1.0, 2e-12);
        }
        near(next.temperature_k(CellId(0)).unwrap(), 325.0, 1e-8);
        near(next.temperature_k(CellId(1)).unwrap(), 375.0, 1e-8);
        assert_eq!(report.advection_rule, InertialAdvectionRule::AffineSwept);
        near(report.thermal.stored_energy_change_j, 0.0, 1e-6);
        assert_eq!(report.thermal.external_heat_j, 0.0);
        assert_eq!(momentum.accepted_time_s(), 0.0);
        assert_eq!(old.time_s(), 0.0);
    });
}

#[test]
fn g4_original_source_refusal_during_affine_coupling_preserves_both_owners() {
    with_cx(|cx| {
        let graph = network(2, &[(0, 1), (1, 0)], &[1.0; 2]);
        let momentum = owner(cx, &graph, &[1.0; 2], &[1.0; 2]);
        let law = Caloric;
        let old = inventory(&law, &[300.0, 400.0], &[0.0; 2]);
        let nodes = [CellId(0), CellId(1)];
        let problem = InertialThermalProblem {
            momentum: &momentum,
            node_cells: &nodes,
            maximum_connections: 4,
        };
        let calls = Cell::new(0usize);
        let failed = advance_inertial_graph_affine(
            cx,
            &problem,
            None,
            &old,
            1.0,
            controls(),
            &[25_000.0, -25_000.0],
            &[0.0; 2],
            |_, _, _| Ok::<_, &'static str>(response(&graph, &[-2.0; 2])),
            |_, _, _, _| {
                calls.set(calls.get() + 1);
                if calls.get() >= 2 {
                    Err("original film source declined")
                } else {
                    Ok(heat(()))
                }
            },
        );
        assert!(matches!(
            failed,
            Err(InertialThermalError::Exchange(
                "original film source declined"
            ))
        ));
        assert_eq!(old.time_s(), 0.0);
        assert_eq!(momentum.accepted_time_s(), 0.0);
        assert_eq!(momentum.volume_flows_m3_s(), vec![1.0; 2]);
        near(old.temperature_k(CellId(1)).unwrap(), 400.0, 0.0);
        assert!(
            advance_inertial_graph_affine(
                cx,
                &problem,
                None,
                &old,
                1.0,
                controls(),
                &[25_000.0, -25_000.0],
                &[0.0; 2],
                |_, _, _| Ok::<_, &'static str>(response(&graph, &[-2.0; 2])),
                |_, _, _, _| Ok(heat(()))
            )
            .is_ok()
        );
    });
}
