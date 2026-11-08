//! Original-law pump work in sealed and reservoir-connected radial networks.
use super::*;

fn pumped_loop() -> TransientCarreauGraph {
    let graph = CarreauPipeGraph::new(
        3,
        vec![
            branch("pump", 0, 1, 0.2, 0.01, true),
            branch("transfer", 1, 2, 0.5, 0.01, true),
            branch("return", 2, 0, 0.3, 0.01, true),
        ],
    )
    .unwrap();
    TransientCarreauGraph::new(graph, 1000.0, &[16; 3], vec![storage(2e-9), None, None]).unwrap()
}

#[test]
fn closed_pumped_loop_matches_original_full_profile_and_signed_work() {
    let model = pumped_loop();
    let mut state =
        with_cx(|cx| model.rest_state(&[Some(Pressure::new(52.0)), None, None], cx)).unwrap();
    let mut workspace = model.workspace().unwrap();
    let pipe = RadialCarreauPipe::new(source(1.0, 0.01, true), 1000.0, 16).unwrap();
    let mut reference = pipe.rest_state().unwrap();
    let mut reference_work = pipe.workspace().unwrap();
    let dt = 0.003;
    let mut recovered_work = false;
    for total_head in [1000.0, 1000.0, -100.0, -1000.0, 0.0] {
        let expected = pipe
            .step(
                &mut reference,
                &mut reference_work,
                total_head,
                dt,
                controls().radial,
                || false,
            )
            .unwrap();
        let report = with_cx(|cx| {
            model.step_with_heads(
                &mut state,
                &mut workspace,
                &[],
                &[],
                &[
                    Pressure::new(total_head),
                    Pressure::new(0.0),
                    Pressure::new(0.0),
                ],
                dt,
                controls(),
                cx,
            )
        })
        .unwrap();
        balance(&report);
        assert_eq!(report.reservoir_work_j, 0.0);
        assert_eq!(report.injection_work_j, 0.0);
        assert!(report.reservoir_outflows.is_empty());
        near(
            state.storage_pressures()[0].unwrap().value(),
            52.0,
            1e-9,
            1e-8,
        );
        near(report.storage_energy_change_j, 0.0, 0.0, 2e-15);
        near(
            report.branch_head_work_j,
            expected.pressure_work_j,
            2e-8,
            2e-15,
        );
        near(report.viscous_heat_j, expected.viscous_heat_j, 2e-8, 2e-15);
        near(
            report.kinetic_energy_change_j,
            expected.kinetic_energy_change_j,
            2e-8,
            2e-15,
        );
        for (edge, length) in [0.2, 0.5, 0.3].into_iter().enumerate() {
            for (&actual, &expected) in state.branches()[edge]
                .velocities_m_s()
                .iter()
                .zip(reference.velocities_m_s())
            {
                near(actual, expected, 2e-8, 1e-11);
            }
            near(
                report.branches[edge].viscous_heat_j,
                length * expected.viscous_heat_j,
                2e-8,
                2e-15,
            );
            near(
                report.branches[edge].midpoint_volume_flow_m3_s,
                expected.midpoint_volume_flow_m3_s,
                2e-8,
                2e-15,
            );
        }
        near(
            report.midpoint_pressures[1].value(),
            52.0 + 0.8 * total_head,
            2e-8,
            1e-7,
        );
        near(
            report.midpoint_pressures[2].value(),
            52.0 + 0.3 * total_head,
            2e-8,
            1e-7,
        );
        recovered_work |= report.branch_head_work_j < 0.0;
    }
    assert!(
        recovered_work,
        "retained fluid inertia can return work to an opposing source"
    );
}

#[test]
fn separate_reservoir_and_head_work_match_the_original_total_pressure_trajectory() {
    let graph = CarreauPipeGraph::new(2, vec![branch("pipe", 0, 1, 0.5, 0.01, true)]).unwrap();
    let model = TransientCarreauGraph::new(graph, 1000.0, &[16], vec![None; 2]).unwrap();
    let mut state = with_cx(|cx| model.rest_state(&[None; 2], cx)).unwrap();
    let mut ordinary = state.clone();
    let mut workspace = model.workspace().unwrap();
    let mut ordinary_workspace = model.workspace().unwrap();
    for (left, right, head) in [
        (200.0, 50.0, 100.0),
        (200.0, 50.0, -180.0),
        (-50.0, 100.0, 150.0),
        (0.0, 0.0, 0.0),
    ] {
        let report = with_cx(|cx| {
            model.step_with_heads(
                &mut state,
                &mut workspace,
                &[boundary(0, left), boundary(1, right)],
                &[],
                &[Pressure::new(head)],
                0.004,
                controls(),
                cx,
            )
        })
        .unwrap();
        let reference = with_cx(|cx| {
            model.step(
                &mut ordinary,
                &mut ordinary_workspace,
                &[boundary(0, left - right + head), boundary(1, 0.0)],
                &[],
                0.004,
                controls(),
                cx,
            )
        })
        .unwrap();
        balance(&report);
        let volume = 0.004 * report.branches[0].midpoint_volume_flow_m3_s;
        near(
            report.reservoir_work_j,
            (left - right) * volume,
            1e-12,
            1e-16,
        );
        near(report.branch_head_work_j, head * volume, 1e-12, 1e-16);
        near(
            report.reservoir_work_j + report.branch_head_work_j,
            reference.reservoir_work_j,
            2e-9,
            2e-15,
        );
        near(report.viscous_heat_j, reference.viscous_heat_j, 2e-9, 2e-15);
        assert_eq!(snapshot(&state), snapshot(&ordinary));
    }
}

#[test]
fn invalid_or_unsupported_head_preserves_full_accepted_history_and_retry() {
    let model = pumped_loop();
    let mut state =
        with_cx(|cx| model.rest_state(&[Some(Pressure::new(52.0)), None, None], cx)).unwrap();
    let mut workspace = model.workspace().unwrap();
    let before = snapshot(&state);
    for heads in [
        vec![Pressure::new(0.0); 2],
        vec![Pressure::new(f64::NAN); 3],
        vec![Pressure::new(f64::MAX); 3],
    ] {
        assert!(
            with_cx(|cx| model.step_with_heads(
                &mut state,
                &mut workspace,
                &[],
                &[],
                &heads,
                0.003,
                controls(),
                cx,
            ))
            .is_err()
        );
        assert_eq!(snapshot(&state), before);
    }
    let mut clean = state.clone();
    let mut clean_workspace = model.workspace().unwrap();
    let heads = [Pressure::new(500.0), Pressure::new(0.0), Pressure::new(0.0)];
    with_cx(|cx| {
        model.step_with_heads(
            &mut state,
            &mut workspace,
            &[],
            &[],
            &heads,
            0.003,
            controls(),
            cx,
        )
    })
    .unwrap();
    with_cx(|cx| {
        model.step_with_heads(
            &mut clean,
            &mut clean_workspace,
            &[],
            &[],
            &heads,
            0.003,
            controls(),
            cx,
        )
    })
    .unwrap();
    assert_eq!(snapshot(&state), snapshot(&clean));
}
