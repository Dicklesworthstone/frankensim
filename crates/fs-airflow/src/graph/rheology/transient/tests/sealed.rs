//! Closed hydraulic circuits are anchored by their actual finite storage.
use super::*;

// Independent simultaneous Newtonian momentum/volume system. There is no
// prescribed-pressure row: two finite stores enclose an ideal middle junction.
fn reference_step(
    old_u: &[Vec<f64>],
    old_p: [f64; 2],
    injection: [f64; 2],
    dt: f64,
) -> (Vec<Vec<f64>>, [f64; 3]) {
    let cells = old_u[0].len();
    let extent = 2 * cells + 3;
    let mut matrix = vec![vec![0.0; extent]; extent];
    let mut rhs = vec![0.0; extent];
    let radius = 0.01;
    let dr = radius / cells as f64;
    for (edge, length) in [0.6, 0.4].into_iter().enumerate() {
        let offset = edge * cells;
        for cell in 0..cells {
            let area = core::f64::consts::PI * dr * dr * (2 * cell + 1) as f64;
            let mass = 950.0 * area * length;
            matrix[offset + cell][offset + cell] = mass;
            rhs[offset + cell] = mass * old_u[edge][cell];
            for (node, sign) in [(edge, 1.0), (edge + 1, -1.0)] {
                matrix[offset + cell][2 * cells + node] -= dt * area * sign;
                matrix[2 * cells + node][offset + cell] += 0.5 * sign * area;
                rhs[2 * cells + node] -= 0.5 * sign * area * old_u[edge][cell];
            }
        }
        for cell in 0..cells {
            let distance = if cell + 1 == cells { 0.5 * dr } else { dr };
            let stiffness =
                0.2 * 2.0 * core::f64::consts::PI * (cell + 1) as f64 * dr * length / distance;
            matrix[offset + cell][offset + cell] += 0.5 * dt * stiffness;
            rhs[offset + cell] -= 0.5 * dt * stiffness * old_u[edge][cell];
            if cell + 1 < cells {
                matrix[offset + cell][offset + cell + 1] -= 0.5 * dt * stiffness;
                matrix[offset + cell + 1][offset + cell] -= 0.5 * dt * stiffness;
                matrix[offset + cell + 1][offset + cell + 1] += 0.5 * dt * stiffness;
                rhs[offset + cell] += 0.5 * dt * stiffness * old_u[edge][cell + 1];
                rhs[offset + cell + 1] +=
                    0.5 * dt * stiffness * (old_u[edge][cell] - old_u[edge][cell + 1]);
            }
        }
    }
    for (index, (node, compliance)) in [(0, 2e-9), (2, 3e-9)].into_iter().enumerate() {
        matrix[2 * cells + node][2 * cells + node] = 2.0 * compliance / dt;
        rhs[2 * cells + node] += injection[index] + 2.0 * compliance / dt * old_p[index];
    }
    let answer = dense(matrix, rhs);
    (
        vec![answer[..cells].to_vec(), answer[cells..2 * cells].to_vec()],
        [
            answer[2 * cells],
            answer[2 * cells + 1],
            answer[2 * cells + 2],
        ],
    )
}

#[test]
fn closed_stores_exchange_pressure_and_full_radial_history_without_a_reservoir() {
    let graph = CarreauPipeGraph::new(
        3,
        vec![
            branch("first", 0, 1, 0.6, 0.01, false),
            branch("second", 1, 2, 0.4, 0.01, false),
        ],
    )
    .unwrap();
    let model = TransientCarreauGraph::new(
        graph,
        950.0,
        &[6, 6],
        vec![storage(2e-9), None, storage(3e-9)],
    )
    .unwrap();
    let mut state = with_cx(|cx| {
        model.rest_state(
            &[Some(Pressure::new(200.0)), None, Some(Pressure::new(-20.0))],
            cx,
        )
    })
    .unwrap();
    let mut workspace = model.workspace().unwrap();
    let mut old_u = vec![vec![0.0; 6]; 2];
    let mut old_p = [200.0, -20.0];
    for (dt, inputs) in [
        (0.003, [0.0, 0.0]),
        (0.007, [0.0, 0.0]),
        (0.004, [1e-7, -2e-8]),
        (0.002, [-3e-8, 8e-8]),
        (0.005, [0.0, 0.0]),
    ] {
        let (expected_u, mid) = reference_step(&old_u, old_p, inputs, dt);
        let report = with_cx(|cx| {
            model.step(
                &mut state,
                &mut workspace,
                &[],
                &[injection(0, inputs[0]), injection(2, inputs[1])],
                dt,
                controls(),
                cx,
            )
        })
        .unwrap();
        balance(&report);
        assert!(report.reservoir_outflows.is_empty());
        assert_eq!(report.reservoir_work_j, 0.0);
        for node in 0..3 {
            near(
                report.midpoint_pressures[node].value(),
                mid[node],
                2e-9,
                1e-8,
            );
        }
        let mut stored_volume = 0.0;
        let mut stored_energy = 0.0;
        for (index, (node, compliance)) in [(0, 2e-9), (2, 3e-9)].into_iter().enumerate() {
            let next = 2.0 * mid[node] - old_p[index];
            near(
                state.storage_pressures()[node].unwrap().value(),
                next,
                2e-9,
                1e-8,
            );
            stored_volume += compliance * (next - old_p[index]);
            stored_energy += 0.5 * compliance * (next + old_p[index]) * (next - old_p[index]);
            old_p[index] = next;
        }
        near(stored_volume, dt * (inputs[0] + inputs[1]), 1e-8, 2e-17);
        near(
            report.storage_volume_changes_m3.iter().sum(),
            stored_volume,
            1e-8,
            2e-17,
        );
        near(report.storage_energy_change_j, stored_energy, 2e-9, 1e-14);
        near(
            report.injection_work_j,
            dt * (mid[0] * inputs[0] + mid[2] * inputs[1]),
            2e-9,
            1e-15,
        );
        for (actual, expected) in state.branches().iter().zip(&expected_u) {
            for (&actual, &expected) in actual.velocities_m_s().iter().zip(expected) {
                near(actual, expected, 2e-9, 1e-12);
            }
        }
        assert!(report.max_final_rigid_node_imbalance.value() < 1e-13);
        old_u = expected_u;
    }
    assert!(old_p[0] < 200.0 && old_p[1] > -20.0);
}

#[test]
fn disconnected_closed_and_reservoir_components_keep_independent_physical_storage() {
    let original = CarreauPipeGraph::new(2, vec![branch("closed", 0, 1, 0.5, 0.01, true)]).unwrap();
    let standalone =
        TransientCarreauGraph::new(original, 1000.0, &[12], vec![storage(2e-9), storage(3e-9)])
            .unwrap();
    let graph = CarreauPipeGraph::new(
        4,
        vec![
            branch("closed", 0, 1, 0.5, 0.01, true),
            branch("open", 2, 3, 0.3, 0.01, true),
        ],
    )
    .unwrap();
    let model = TransientCarreauGraph::new(
        graph,
        1000.0,
        &[12, 12],
        vec![storage(2e-9), storage(3e-9), None, storage(4e-9)],
    )
    .unwrap();
    let pressures = [Some(Pressure::new(300.0)), Some(Pressure::new(-50.0))];
    let mut isolated = with_cx(|cx| standalone.rest_state(&pressures, cx)).unwrap();
    let mut state = with_cx(|cx| {
        model.rest_state(
            &[
                pressures[0],
                pressures[1],
                None,
                Some(Pressure::new(1000.0)),
            ],
            cx,
        )
    })
    .unwrap();
    let mut a = standalone.workspace().unwrap();
    let mut b = model.workspace().unwrap();
    for pressure in [800.0, 1600.0, 500.0] {
        let first =
            with_cx(|cx| standalone.step(&mut isolated, &mut a, &[], &[], 0.003, controls(), cx))
                .unwrap();
        let combined = with_cx(|cx| {
            model.step(
                &mut state,
                &mut b,
                &[boundary(2, pressure)],
                &[],
                0.003,
                controls(),
                cx,
            )
        })
        .unwrap();
        balance(&combined);
        for node in 0..2 {
            near(
                state.storage_pressures()[node].unwrap().value(),
                isolated.storage_pressures()[node].unwrap().value(),
                1e-8,
                1e-8,
            );
        }
        for (&actual, &expected) in state.branches()[0]
            .velocities_m_s()
            .iter()
            .zip(isolated.branches()[0].velocities_m_s())
        {
            near(actual, expected, 1e-8, 1e-11);
        }
        near(
            combined.branches[0].viscous_heat_j,
            first.viscous_heat_j,
            1e-8,
            1e-15,
        );
        near(
            combined.storage_volume_changes_m3[0] + combined.storage_volume_changes_m3[1],
            0.0,
            0.0,
            3e-17,
        );
        near(
            combined.reservoir_outflows[0].value() * 0.003,
            combined.storage_volume_changes_m3[3],
            1e-8,
            3e-17,
        );
    }
}

#[test]
fn storage_anchors_do_not_pin_pressure_or_admit_floating_rigid_components() {
    let graph = CarreauPipeGraph::new(
        4,
        vec![
            branch("stored", 0, 1, 0.5, 0.01, false),
            branch("floating", 2, 3, 0.5, 0.01, false),
        ],
    )
    .unwrap();
    let model = TransientCarreauGraph::new(
        graph,
        1000.0,
        &[8, 8],
        vec![storage(2e-9), storage(3e-9), None, None],
    )
    .unwrap();
    let mut state = with_cx(|cx| {
        model.rest_state(
            &[
                Some(Pressure::new(20.0)),
                Some(Pressure::new(20.0)),
                None,
                None,
            ],
            cx,
        )
    })
    .unwrap();
    let mut workspace = model.workspace().unwrap();
    let before = snapshot(&state);
    assert!(
        with_cx(|cx| model.step(&mut state, &mut workspace, &[], &[], 0.003, controls(), cx,))
            .is_err()
    );
    assert_eq!(snapshot(&state), before);
    assert!(
        with_cx(|cx| model.step(
            &mut state,
            &mut workspace,
            &[boundary(0, 20.0), boundary(2, 0.0)],
            &[],
            0.003,
            controls(),
            cx,
        ))
        .is_err(),
        "a physical store cannot simultaneously be an ideal pressure source"
    );
    assert_eq!(snapshot(&state), before);
    assert!(
        with_cx(|cx| model.graph().solve(
            &[boundary(2, 0.0)],
            super::super::super::CarreauGraphControls::default(),
            cx,
        ))
        .is_err(),
        "steady flow still needs actual pressure data in every component"
    );
    let report = with_cx(|cx| {
        model.step(
            &mut state,
            &mut workspace,
            &[boundary(2, 0.0)],
            &[],
            0.003,
            controls(),
            cx,
        )
    })
    .unwrap();
    balance(&report);
    for pressure in state.storage_pressures().iter().flatten() {
        assert_eq!(pressure.value(), 20.0);
    }
    assert_eq!(report.viscous_heat_j, 0.0);
}
