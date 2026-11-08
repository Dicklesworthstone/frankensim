//! Synthetic G1/G3 coupled-flow fixtures, not calibrated material data.

mod sealed;

use super::super::CarreauPipeBranch;
use super::*;
use crate::graph::tests::with_cx;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_material::fluid_transport::rheology::pipe::LaminarCarreauPipe;
use fs_material::fluid_transport::rheology::{
    CarreauYasuda, ThermalCarreauYasuda, ViscosityShiftMode,
};
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};

fn near(a: f64, b: f64, relative: f64, absolute: f64) {
    assert!(a.is_finite() && b.is_finite());
    assert!(
        (a - b).abs() <= absolute + relative * a.abs().max(b.abs()),
        "{a:.16e} != {b:.16e}"
    );
}

fn source(length: f64, radius: f64, nonlinear: bool) -> LaminarCarreauPipe {
    let law = ThermalCarreauYasuda::new(
        CarreauYasuda::new(
            0.2,
            if nonlinear { 0.01 } else { 0.2 },
            0.8,
            0.5,
            1.2,
            [0.0, 1e7],
        )
        .unwrap(),
        TemperatureShift::new(
            TemperatureShiftModel::Arrhenius {
                activation_energy_j_mol: 9000.0,
            },
            300.0,
            [280.0, 360.0],
        )
        .unwrap(),
        ViscosityShiftMode::Vertical,
    );
    LaminarCarreauPipe::new(law, 300.0, radius, length, 1e5).unwrap()
}

fn branch(
    name: &str,
    from: usize,
    to: usize,
    length: f64,
    radius: f64,
    nonlinear: bool,
) -> CarreauPipeBranch {
    CarreauPipeBranch {
        name: name.into(),
        from,
        to,
        pipe: source(length, radius, nonlinear),
    }
}

fn boundary(node: usize, pressure: f64) -> FixedPressure {
    FixedPressure {
        node,
        pressure: Pressure::new(pressure),
    }
}

fn injection(node: usize, flow: f64) -> FlowInjection {
    FlowInjection {
        node,
        flow: VolumetricFlowRate::new(flow),
    }
}

fn storage(compliance: f64) -> Option<LinearHydraulicStorage> {
    Some(LinearHydraulicStorage {
        compliance_m3_pa: compliance,
        minimum_pressure: Pressure::new(-1e5),
        maximum_pressure: Pressure::new(1e5),
    })
}

fn controls() -> TransientCarreauControls {
    TransientCarreauControls {
        absolute_flow_tolerance: VolumetricFlowRate::new(2e-15),
        relative_flow_tolerance: 1e-11,
        radial: RadialCarreauControls {
            absolute_velocity_tolerance_m_s: 1e-14,
            relative_velocity_tolerance: 1e-12,
            ..RadialCarreauControls::default()
        },
        ..TransientCarreauControls::default()
    }
}

fn snapshot(state: &TransientCarreauState) -> (u64, Vec<Option<u64>>, Vec<Vec<u64>>) {
    (
        state.accepted_time_s.to_bits(),
        state
            .storage_pressures
            .iter()
            .map(|p| p.map(|p| p.value().to_bits()))
            .collect(),
        state
            .branches
            .iter()
            .map(|b| b.velocities_m_s().iter().map(|x| x.to_bits()).collect())
            .collect(),
    )
}

fn balance(step: &TransientCarreauStep) {
    assert!(step.viscous_heat_j >= 0.0);
    assert!(step.accounting_residual_j.abs() <= step.accounting_tolerance_j);
    let heat: f64 = step.branches.iter().map(|b| b.viscous_heat_j).sum();
    let kinetic: f64 = step
        .branches
        .iter()
        .map(|b| b.kinetic_energy_change_j)
        .sum();
    near(step.viscous_heat_j, heat, 2e-15, 1e-22);
    near(step.kinetic_energy_change_j, kinetic, 2e-15, 1e-22);
    near(
        step.physical_work_residual_j,
        step.radial_work_residual_j + step.storage_midpoint_work_residual_j
            - step.node_residual_work_j,
        1e-6,
        step.accounting_tolerance_j + 1e-22,
    );
    assert!(
        step.branches
            .iter()
            .all(|b| b.midpoint_pressure_conductance_m3_s_pa > 0.0)
    );
}

#[test]
fn g1_nonlinear_series_retains_the_single_original_pipe_trajectory() {
    let graph = CarreauPipeGraph::new(
        3,
        vec![
            branch("short", 0, 1, 0.2, 0.01, true),
            branch("long", 1, 2, 0.8, 0.01, true),
        ],
    )
    .unwrap();
    let model = TransientCarreauGraph::new(graph, 1000.0, &[24, 24], vec![None; 3]).unwrap();
    let mut state = with_cx(|cx| model.rest_state(&[None; 3], cx)).unwrap();
    let mut work = model.workspace().unwrap();
    let whole = RadialCarreauPipe::new(source(1.0, 0.01, true), 1000.0, 24).unwrap();
    let mut original = whole.rest_state().unwrap();
    let mut original_work = whole.workspace().unwrap();
    for pressure in [4000.0, 4000.0, 1200.0, -1000.0, -2000.0, 0.0] {
        let report = with_cx(|cx| {
            model.step(
                &mut state,
                &mut work,
                &[boundary(0, pressure), boundary(2, 0.0)],
                &[],
                0.003,
                controls(),
                cx,
            )
        })
        .unwrap();
        let reference = whole
            .step(
                &mut original,
                &mut original_work,
                pressure,
                0.003,
                controls().radial,
                || false,
            )
            .unwrap();
        balance(&report);
        near(
            report.midpoint_pressures[1].value(),
            0.8 * pressure,
            5e-8,
            1e-8,
        );
        for actual in state.branches() {
            for (&a, &b) in actual
                .velocities_m_s()
                .iter()
                .zip(original.velocities_m_s())
            {
                near(a, b, 2e-8, 1e-11);
            }
        }
        near(report.viscous_heat_j, reference.viscous_heat_j, 2e-8, 1e-17);
        near(
            report.kinetic_energy_change_j,
            reference.kinetic_energy_change_j,
            2e-8,
            1e-17,
        );
        assert!(report.max_final_rigid_node_imbalance.value() < 3e-15);
    }
}

// Independent test-only direct solve of the simultaneous Newtonian annular
// velocity/storage equations. Production uses neither this dense matrix nor
// its Gaussian elimination; it eliminates the original nonlinear radial laws.
fn dense(mut matrix: Vec<Vec<f64>>, mut rhs: Vec<f64>) -> Vec<f64> {
    let n = rhs.len();
    for k in 0..n {
        let pivot = (k..n)
            .max_by(|&a, &b| matrix[a][k].abs().total_cmp(&matrix[b][k].abs()))
            .unwrap();
        matrix.swap(k, pivot);
        rhs.swap(k, pivot);
        assert!(matrix[k][k].abs() > 1e-25);
        for i in k + 1..n {
            let factor = matrix[i][k] / matrix[k][k];
            for j in k..n {
                matrix[i][j] -= factor * matrix[k][j];
            }
            rhs[i] -= factor * rhs[k];
        }
    }
    let mut answer = vec![0.0; n];
    for i in (0..n).rev() {
        let known: f64 = (i + 1..n).map(|j| matrix[i][j] * answer[j]).sum();
        answer[i] = (rhs[i] - known) / matrix[i][i];
    }
    answer
}

fn dense_storage_step(
    old_u: &[Vec<f64>],
    old_p: [f64; 2],
    reservoir: f64,
    dt: f64,
) -> (Vec<Vec<f64>>, [f64; 2]) {
    let cells = old_u[0].len();
    let unknowns = 2 * cells + 2;
    let mut a = vec![vec![0.0; unknowns]; unknowns];
    let mut b = vec![0.0; unknowns];
    for (edge, (length, radius)) in [(0.6, 0.01), (0.4, 0.013)].into_iter().enumerate() {
        let dr = radius / cells as f64;
        let offset = edge * cells;
        for i in 0..cells {
            let area = core::f64::consts::PI * dr * dr * (2 * i + 1) as f64;
            let mass = 950.0 * area * length;
            a[offset + i][offset + i] += mass;
            b[offset + i] += mass * old_u[edge][i];
            for (node, sign) in [(edge, 1.0), (edge + 1, -1.0)] {
                if node == 0 {
                    b[offset + i] += dt * area * sign * reservoir;
                } else {
                    a[offset + i][2 * cells + node - 1] -= dt * area * sign;
                    a[2 * cells + node - 1][offset + i] += 0.5 * sign * area;
                    b[2 * cells + node - 1] -= 0.5 * sign * area * old_u[edge][i];
                }
            }
        }
        for i in 0..cells {
            let distance = if i + 1 == cells { 0.5 * dr } else { dr };
            let stiffness =
                0.2 * 2.0 * core::f64::consts::PI * (i + 1) as f64 * dr * length / distance;
            a[offset + i][offset + i] += 0.5 * dt * stiffness;
            b[offset + i] -= 0.5 * dt * stiffness * old_u[edge][i];
            if i + 1 < cells {
                a[offset + i][offset + i + 1] -= 0.5 * dt * stiffness;
                a[offset + i + 1][offset + i] -= 0.5 * dt * stiffness;
                a[offset + i + 1][offset + i + 1] += 0.5 * dt * stiffness;
                b[offset + i] += 0.5 * dt * stiffness * old_u[edge][i + 1];
                b[offset + i + 1] += 0.5 * dt * stiffness * (old_u[edge][i] - old_u[edge][i + 1]);
            }
        }
    }
    for (i, c) in [2e-9, 3e-9].into_iter().enumerate() {
        a[2 * cells + i][2 * cells + i] = 2.0 * c / dt;
        b[2 * cells + i] += 2.0 * c / dt * old_p[i];
    }
    let solved = dense(a, b);
    (
        vec![solved[..cells].to_vec(), solved[cells..2 * cells].to_vec()],
        [solved[2 * cells], solved[2 * cells + 1]],
    )
}

#[test]
fn g1_two_storage_nodes_match_independent_coupled_radial_matrix() {
    let graph = CarreauPipeGraph::new(
        3,
        vec![
            branch("feed", 0, 1, 0.6, 0.01, false),
            branch("transfer", 1, 2, 0.4, 0.013, false),
        ],
    )
    .unwrap();
    let model = TransientCarreauGraph::new(
        graph,
        950.0,
        &[6, 6],
        vec![None, storage(2e-9), storage(3e-9)],
    )
    .unwrap();
    let mut state = with_cx(|cx| {
        model.rest_state(
            &[None, Some(Pressure::new(12.0)), Some(Pressure::new(-7.0))],
            cx,
        )
    })
    .unwrap();
    let mut work = model.workspace().unwrap();
    let mut old_u = vec![vec![0.0; 6]; 2];
    let mut old_p = [12.0, -7.0];
    for pressure in [200.0, 50.0, -100.0, 0.0, 150.0] {
        let (expected_u, mid_p) = dense_storage_step(&old_u, old_p, pressure, 0.003);
        let report = with_cx(|cx| {
            model.step(
                &mut state,
                &mut work,
                &[boundary(0, pressure)],
                &[],
                0.003,
                controls(),
                cx,
            )
        })
        .unwrap();
        balance(&report);
        let mut expected_storage_change = 0.0;
        for i in 0..2 {
            near(
                report.midpoint_pressures[i + 1].value(),
                mid_p[i],
                2e-10,
                1e-9,
            );
            let next = 2.0 * mid_p[i] - old_p[i];
            near(
                state.storage_pressures()[i + 1].unwrap().value(),
                next,
                2e-10,
                1e-9,
            );
            expected_storage_change += 0.5 * [2e-9, 3e-9][i] * (next * next - old_p[i] * old_p[i]);
            old_p[i] = next;
            for (&actual, &expected) in state.branches()[i]
                .velocities_m_s()
                .iter()
                .zip(&expected_u[i])
            {
                near(actual, expected, 2e-10, 1e-13);
            }
        }
        near(
            report.storage_energy_change_j,
            expected_storage_change,
            2e-10,
            1e-15,
        );
        let storage_volume: f64 = report.storage_volume_changes_m3.iter().sum();
        near(
            report.reservoir_outflows[0].value() * 0.003,
            storage_volume,
            1e-9,
            2e-17,
        );
        old_u = expected_u;
    }
}

fn loop_model() -> TransientCarreauGraph {
    let graph = CarreauPipeGraph::new(
        4,
        vec![
            branch("upper-feed", 0, 1, 0.3, 0.01, true),
            branch("upper-return", 1, 3, 0.5, 0.012, true),
            branch("lower-feed", 0, 2, 0.7, 0.013, true),
            branch("lower-return", 2, 3, 0.4, 0.01, true),
            branch("cross", 1, 2, 0.25, 0.007, true),
        ],
    )
    .unwrap();
    TransientCarreauGraph::new(graph, 1000.0, &[16; 5], vec![None; 4]).unwrap()
}

#[test]
fn g3_nonlinear_loop_conserves_endpoints_and_reverses_with_its_radial_history() {
    let model = loop_model();
    let mut positive = with_cx(|cx| model.rest_state(&[None; 4], cx)).unwrap();
    let mut negative = positive.clone();
    let mut a = model.workspace().unwrap();
    let mut b = model.workspace().unwrap();
    for pressure in [5000.0, 5000.0, 1000.0, -2000.0, 0.0] {
        let forward = with_cx(|cx| {
            model.step(
                &mut positive,
                &mut a,
                &[boundary(0, pressure), boundary(3, 0.0)],
                &[],
                0.004,
                controls(),
                cx,
            )
        })
        .unwrap();
        let reverse = with_cx(|cx| {
            model.step(
                &mut negative,
                &mut b,
                &[boundary(0, -pressure), boundary(3, 0.0)],
                &[],
                0.004,
                controls(),
                cx,
            )
        })
        .unwrap();
        balance(&forward);
        balance(&reverse);
        near(forward.viscous_heat_j, reverse.viscous_heat_j, 2e-10, 1e-18);
        near(
            forward.reservoir_work_j,
            reverse.reservoir_work_j,
            2e-10,
            1e-18,
        );
        for (u, v) in positive.branches().iter().zip(negative.branches()) {
            for (&u, &v) in u.velocities_m_s().iter().zip(v.velocities_m_s()) {
                near(u, -v, 2e-9, 1e-12);
            }
        }
        for node in [1, 2] {
            assert!(forward.node_outflows[node].value().abs() < 2e-14);
        }
        assert!(forward.max_final_rigid_node_imbalance.value() < 2e-14);
    }
}

#[test]
fn g3_storage_and_fixed_node_injections_preserve_separate_work_and_volume() {
    let graph = CarreauPipeGraph::new(2, vec![branch("line", 0, 1, 0.5, 0.01, true)]).unwrap();
    let model =
        TransientCarreauGraph::new(graph, 1000.0, &[20], vec![None, storage(1e-8)]).unwrap();
    let mut a = with_cx(|cx| model.rest_state(&[None, Some(Pressure::new(0.0))], cx)).unwrap();
    let mut b = a.clone();
    let mut wa = model.workspace().unwrap();
    let mut wb = model.workspace().unwrap();
    for _ in 0..5 {
        let ordinary = with_cx(|cx| {
            model.step(
                &mut a,
                &mut wa,
                &[boundary(0, 50.0)],
                &[injection(1, 1e-6)],
                0.002,
                controls(),
                cx,
            )
        })
        .unwrap();
        let ports = with_cx(|cx| {
            model.step(
                &mut b,
                &mut wb,
                &[boundary(0, 50.0)],
                &[injection(1, 1e-6), injection(0, 4e-6)],
                0.002,
                controls(),
                cx,
            )
        })
        .unwrap();
        assert_eq!(snapshot(&a), snapshot(&b));
        assert_eq!(ordinary.branches, ports.branches);
        balance(&ordinary);
        balance(&ports);
        near(
            ordinary.storage_volume_changes_m3[1],
            0.002 * (ordinary.reservoir_outflows[0].value() + 1e-6),
            1e-9,
            1e-17,
        );
        near(
            ports.reservoir_work_j - ordinary.reservoir_work_j,
            -50.0 * 4e-6 * 0.002,
            1e-11,
            1e-18,
        );
        near(
            ports.injection_work_j - ordinary.injection_work_j,
            50.0 * 4e-6 * 0.002,
            1e-11,
            1e-18,
        );
    }
}

fn profile_flow(pipe: &RadialCarreauPipe, u: &[f64]) -> f64 {
    let area = core::f64::consts::PI * pipe.spacing_m() * pipe.spacing_m();
    u.iter()
        .enumerate()
        .map(|(i, &u)| area * (2 * i + 1) as f64 * u)
        .sum()
}

#[test]
fn g3_held_rigid_node_injection_requires_and_preserves_consistent_history() {
    let graph =
        CarreauPipeGraph::new(2, vec![branch("controlled", 0, 1, 0.5, 0.01, true)]).unwrap();
    let model = TransientCarreauGraph::new(graph, 1000.0, &[20], vec![None; 2]).unwrap();
    let u: Vec<_> = (0..20)
        .map(|i| 0.03 * (1.0 - ((i as f64 + 0.5) / 20.0).powi(2)))
        .collect();
    let flow = profile_flow(&model.pipes()[0], &u);
    let mut state =
        with_cx(|cx| model.state_from_velocities(vec![u], &[None; 2], 0.2, cx)).unwrap();
    let mut work = model.workspace().unwrap();
    for _ in 0..4 {
        let report = with_cx(|cx| {
            model.step(
                &mut state,
                &mut work,
                &[boundary(1, 0.0)],
                &[injection(0, flow)],
                0.001,
                controls(),
                cx,
            )
        })
        .unwrap();
        near(report.branches[0].final_volume_flow_m3_s, flow, 1e-9, 2e-15);
        balance(&report);
    }
    let before = snapshot(&state);
    assert!(matches!(
        with_cx(|cx| model.step(
            &mut state,
            &mut work,
            &[boundary(1, 0.0)],
            &[injection(0, flow * 1.1)],
            0.001,
            controls(),
            cx
        )),
        Err(TransientCarreauError::InconsistentInitialFlow { node: 0, .. })
    ));
    assert_eq!(snapshot(&state), before);
}

#[test]
fn g3_equal_mean_flow_histories_keep_distinct_network_relaxation() {
    let graph = CarreauPipeGraph::new(
        3,
        vec![
            branch("left", 0, 1, 0.3, 0.01, true),
            branch("right", 1, 2, 0.7, 0.01, true),
        ],
    )
    .unwrap();
    let model = TransientCarreauGraph::new(graph, 1000.0, &[16, 16], vec![None; 3]).unwrap();
    let u = vec![0.05; 16];
    let mut v = u.clone();
    v[0] += 0.04;
    v[1] -= 0.04 / 3.0;
    near(
        profile_flow(&model.pipes()[0], &u),
        profile_flow(&model.pipes()[0], &v),
        1e-15,
        1e-20,
    );
    let mut a =
        with_cx(|cx| model.state_from_velocities(vec![u.clone(), u], &[None; 3], 0.0, cx)).unwrap();
    let mut b =
        with_cx(|cx| model.state_from_velocities(vec![v.clone(), v], &[None; 3], 0.0, cx)).unwrap();
    let mut wa = model.workspace().unwrap();
    let mut wb = model.workspace().unwrap();
    let ra = with_cx(|cx| {
        model.step(
            &mut a,
            &mut wa,
            &[boundary(0, 0.0), boundary(2, 0.0)],
            &[],
            0.1,
            controls(),
            cx,
        )
    })
    .unwrap();
    let rb = with_cx(|cx| {
        model.step(
            &mut b,
            &mut wb,
            &[boundary(0, 0.0), boundary(2, 0.0)],
            &[],
            0.1,
            controls(),
            cx,
        )
    })
    .unwrap();
    balance(&ra);
    balance(&rb);
    assert!(
        (ra.branches[0].final_volume_flow_m3_s - rb.branches[0].final_volume_flow_m3_s).abs()
            > 1e-12
    );
    assert_eq!(ra.reservoir_work_j, 0.0);
    assert_eq!(rb.reservoir_work_j, 0.0);
    assert!(ra.kinetic_energy_change_j < 0.0 && rb.kinetic_energy_change_j < 0.0);
}

#[test]
fn g3_nonlinear_pressure_and_storage_schur_tangent_matches_directional_difference() {
    let graph = CarreauPipeGraph::new(
        3,
        vec![
            branch("first", 0, 1, 0.3, 0.01, true),
            branch("second", 1, 2, 0.7, 0.012, true),
        ],
    )
    .unwrap();
    let model = TransientCarreauGraph::new(
        graph,
        1000.0,
        &[20, 20],
        vec![None, storage(2e-9), storage(3e-9)],
    )
    .unwrap();
    let state = with_cx(|cx| {
        model.rest_state(
            &[None, Some(Pressure::new(100.0)), Some(Pressure::new(150.0))],
            cx,
        )
    })
    .unwrap();
    let mut work = model.workspace().unwrap();
    with_cx(|cx| {
        let layout = model
            .graph
            .boundary_layout(&[boundary(0, 900.0)], &[], cx)
            .unwrap();
        let p = [900.0, 300.0, 100.0];
        let direction = [0.0, 1.2, -0.7];
        let dt = 0.01;
        let mut count = 0;
        model
            .evaluate(
                &state,
                &p,
                &layout,
                dt,
                controls(),
                &mut work.radial,
                &mut work.outputs,
                &mut work.current,
                &mut count,
                cx,
            )
            .unwrap();
        let slopes = work.current.conductance.clone();
        let eps = 1e-2;
        let plus: Vec<_> = p.iter().zip(direction).map(|(p, d)| p + eps * d).collect();
        model
            .evaluate(
                &state,
                &plus,
                &layout,
                dt,
                controls(),
                &mut work.radial,
                &mut work.outputs,
                &mut work.trial,
                &mut count,
                cx,
            )
            .unwrap();
        let positive = work.trial.residuals.clone();
        let minus: Vec<_> = p.iter().zip(direction).map(|(p, d)| p - eps * d).collect();
        model
            .evaluate(
                &state,
                &minus,
                &layout,
                dt,
                controls(),
                &mut work.radial,
                &mut work.outputs,
                &mut work.trial,
                &mut count,
                cx,
            )
            .unwrap();
        for node in [1, 2] {
            let mut expected =
                2.0 * model.storage[node].unwrap().compliance_m3_pa / dt * direction[node];
            for &edge in &model.graph.adjacency[node] {
                let other = model.graph.other(edge, node);
                expected += slopes[edge] * (direction[node] - direction[other]);
            }
            near(
                (positive[node] - work.trial.residuals[node]) / (2.0 * eps),
                expected,
                2e-7,
                1e-14,
            );
        }
    });
}

#[test]
fn g4_partial_branch_budget_and_storage_refusal_leave_all_accepted_state_unchanged() {
    let graph = CarreauPipeGraph::new(
        3,
        vec![
            branch("first", 0, 1, 0.3, 0.01, true),
            branch("second", 1, 2, 0.7, 0.01, true),
        ],
    )
    .unwrap();
    let model = TransientCarreauGraph::new(graph, 1000.0, &[16, 16], vec![None; 3]).unwrap();
    let mut state = with_cx(|cx| model.rest_state(&[None; 3], cx)).unwrap();
    let mut clean = state.clone();
    let mut work = model.workspace().unwrap();
    let mut clean_work = model.workspace().unwrap();
    let before = snapshot(&state);
    let limited = TransientCarreauControls {
        maximum_branch_solves: 1,
        ..controls()
    };
    assert!(matches!(
        with_cx(|cx| model.step(
            &mut state,
            &mut work,
            &[boundary(0, 1000.0), boundary(2, 0.0)],
            &[],
            0.003,
            limited,
            cx
        )),
        Err(TransientCarreauError::BranchSolveBudget)
    ));
    assert_eq!(snapshot(&state), before);
    for bad_dt in [0.0, -0.1, f64::NAN] {
        assert!(
            with_cx(|cx| model.step(
                &mut state,
                &mut work,
                &[boundary(0, 1000.0), boundary(2, 0.0)],
                &[],
                bad_dt,
                controls(),
                cx
            ))
            .is_err()
        );
        assert_eq!(snapshot(&state), before);
    }
    let retried = with_cx(|cx| {
        model.step(
            &mut state,
            &mut work,
            &[boundary(0, 1000.0), boundary(2, 0.0)],
            &[],
            0.003,
            controls(),
            cx,
        )
    })
    .unwrap();
    let reference = with_cx(|cx| {
        model.step(
            &mut clean,
            &mut clean_work,
            &[boundary(0, 1000.0), boundary(2, 0.0)],
            &[],
            0.003,
            controls(),
            cx,
        )
    })
    .unwrap();
    assert_eq!(snapshot(&state), snapshot(&clean));
    assert_eq!(retried.branches, reference.branches);
    let mut constrained = model.graph().branches().to_vec();
    constrained[1].pipe =
        LaminarCarreauPipe::new(*constrained[1].pipe.law(), 300.0, 0.01, 0.7, 1e-12).unwrap();
    let constrained = TransientCarreauGraph::new(
        CarreauPipeGraph::new(3, constrained).unwrap(),
        1000.0,
        &[16, 16],
        vec![None; 3],
    )
    .unwrap();
    let mut constrained_state = with_cx(|cx| constrained.rest_state(&[None; 3], cx)).unwrap();
    let mut constrained_work = constrained.workspace().unwrap();
    let constrained_before = snapshot(&constrained_state);
    assert!(matches!(
        with_cx(|cx| constrained.step(
            &mut constrained_state,
            &mut constrained_work,
            &[boundary(0, 1000.0), boundary(2, 0.0)],
            &[],
            0.003,
            controls(),
            cx
        )),
        Err(TransientCarreauError::Radial { branch: 1, .. })
    ));
    assert_eq!(snapshot(&constrained_state), constrained_before);
    let graph =
        CarreauPipeGraph::new(2, vec![branch("bounded-storage", 0, 1, 0.5, 0.01, true)]).unwrap();
    let model = TransientCarreauGraph::new(
        graph,
        1000.0,
        &[16],
        vec![
            None,
            Some(LinearHydraulicStorage {
                compliance_m3_pa: 1e-10,
                minimum_pressure: Pressure::new(-1.0),
                maximum_pressure: Pressure::new(1.0),
            }),
        ],
    )
    .unwrap();
    let mut state = with_cx(|cx| model.rest_state(&[None, Some(Pressure::new(0.0))], cx)).unwrap();
    let mut work = model.workspace().unwrap();
    let before = snapshot(&state);
    assert!(
        with_cx(|cx| model.step(
            &mut state,
            &mut work,
            &[boundary(0, 1000.0)],
            &[injection(1, 1e-4)],
            0.01,
            controls(),
            cx
        ))
        .is_err()
    );
    assert_eq!(snapshot(&state), before);
    assert!(
        with_cx(|cx| model.step(
            &mut state,
            &mut work,
            &[boundary(0, 0.0), boundary(1, 0.0)],
            &[],
            0.01,
            controls(),
            cx
        ))
        .is_err()
    );
    assert_eq!(snapshot(&state), before);
}

#[test]
fn g4_context_cancellation_and_foreign_histories_publish_nothing() {
    let model = loop_model();
    let other = loop_model();
    let mut state = with_cx(|cx| model.rest_state(&[None; 4], cx)).unwrap();
    let mut work = model.workspace().unwrap();
    let before = snapshot(&state);
    let gate = CancelGate::new();
    gate.request();
    let refused = ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 17,
                kernel_id: 711,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        model.step(
            &mut state,
            &mut work,
            &[boundary(0, 100.0), boundary(3, 0.0)],
            &[],
            0.01,
            controls(),
            &cx,
        )
    });
    assert!(matches!(
        refused,
        Err(TransientCarreauError::Graph(CarreauGraphError::Graph(
            GraphError::Interrupted
        )))
    ));
    assert_eq!(snapshot(&state), before);
    assert!(matches!(
        with_cx(|cx| other.step(
            &mut state,
            &mut work,
            &[boundary(0, 100.0), boundary(3, 0.0)],
            &[],
            0.01,
            controls(),
            cx
        )),
        Err(TransientCarreauError::StateMismatch)
    ));
    assert_eq!(snapshot(&state), before);
    assert!(
        with_cx(|cx| model.step(&mut state, &mut work, &[], &[], 0.01, controls(), cx)).is_err()
    );
    assert_eq!(snapshot(&state), before);
    let report = with_cx(|cx| {
        model.step(
            &mut state,
            &mut work,
            &[boundary(0, 100.0), boundary(3, 0.0)],
            &[],
            0.01,
            controls(),
            cx,
        )
    })
    .unwrap();
    balance(&report);
}
