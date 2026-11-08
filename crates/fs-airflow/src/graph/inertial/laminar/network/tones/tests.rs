use super::*;
use crate::graph::inertial::laminar::LaminarPipeSpec;
use crate::graph::tests::with_cx;
use fs_material::fluid_transport::FluidTransportState;
use fs_phs::ZkDiffusiveBasis;
use fs_qty::{Length, VolumetricFlowRate};

fn budget() -> DenseExponentialBudget {
    DenseExponentialBudget {
        max_dimension: 256,
        max_dense_bytes: 128 * 1024 * 1024,
        max_multiply_terms: 10_000_000_000,
    }
}

fn pipe(cx: &Cx<'_>, viscosity: f64) -> LaminarPipe<FluidTransportState> {
    LaminarPipe::try_new(
        cx,
        FluidTransportState::try_new(1000.0, 4200.0, 0.6, viscosity).unwrap(),
        LaminarPipeSpec {
            length: Length::new(0.4),
            radius: Length::new(0.0005),
            maximum_reynolds: 1000.0,
            band_rad_s: [0.01, 10.0],
            maximum_relative_wall_error: 1e-4,
        },
        &ZkDiffusiveBasis::new(4).unwrap(),
    )
    .unwrap()
}

fn fork(cx: &Cx<'_>) -> LaminarNetwork<FluidTransportState> {
    LaminarNetwork::try_new(
        cx,
        LaminarNetworkTopology {
            node_count: 4,
            branch_nodes: vec![[0, 1], [1, 2], [1, 3]],
            pressure_nodes: vec![0, 2, 3],
        },
        vec![pipe(cx, 0.001), pipe(cx, 0.002), pipe(cx, 0.003)],
    )
    .unwrap()
}

fn tone(c: f64, s: f64) -> LaminarPressureTone {
    LaminarPressureTone {
        cosine: Pressure::new(c),
        sine: Pressure::new(s),
    }
}

fn source() -> [LaminarPressureTone; 12] {
    [
        tone(20.0, 0.0),
        tone(10.0, 2.0),
        tone(-5.0, 0.0),
        tone(-3.0, 7.0),
        tone(2.0, 0.0),
        tone(5.0, -4.0),
        tone(8.0, 0.0),
        tone(1.0, 0.0),
        tone(-12.0, 0.0),
        tone(0.0, -6.0),
        tone(4.0, 0.0),
        tone(-2.0, 3.0),
    ]
}

fn near(a: f64, b: f64, scale: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * scale,
        "{a:e} != {b:e}, scale={scale:e}"
    );
}

#[test]
fn single_branch_recovers_original_continuous_pipe_and_profile() {
    with_cx(|cx| {
        let network = LaminarNetwork::try_new(
            cx,
            LaminarNetworkTopology {
                node_count: 2,
                branch_nodes: vec![[0, 1]],
                pressure_nodes: vec![0, 1],
            },
            vec![pipe(cx, 0.001)],
        )
        .unwrap();
        let single = &network.branches()[0];
        let mut state = network.rest_state(0.0).unwrap();
        let mut original = single
            .steady_state(VolumetricFlowRate::new(0.0), 0.0)
            .unwrap();
        let mut pressures = [
            tone(20.0, 0.0),
            tone(10.0, 2.0),
            tone(-5.0, 0.0),
            tone(-3.0, 7.0),
            tone(8.0, 0.0),
            tone(1.0, 0.0),
        ];
        let mut difference = [tone(33.0, 0.0), tone(14.0, -5.0)];
        let mut prepared = network
            .prepare_tones(cx, &[0.0, 5.0], 0.02, budget())
            .unwrap();
        let reference = single
            .prepare_tones(cx, &[0.0, 5.0], 0.02, budget())
            .unwrap()
            .advance(cx, &mut original, &mut difference)
            .unwrap();
        let result = prepared.advance(cx, &mut state, &mut pressures).unwrap();
        let branch = prepared.branch_steps()[0];
        near(
            branch.final_volume_flow_m3_s,
            reference.final_volume_flow_m3_s,
            1e-9,
            2e-10,
        );
        near(
            branch.displaced_volume_m3,
            reference.displaced_volume_m3,
            1e-11,
            2e-10,
        );
        near(
            branch.viscous_heat_j,
            reference.native.dissipated,
            reference.native.dissipated,
            2e-9,
        );
        near(
            result.native.supplied,
            reference.native.supplied,
            reference.native.supplied.abs(),
            2e-9,
        );
        let extracted = network.branch_state(&state, 0).unwrap();
        for (&a, &b) in extracted.coordinates().iter().zip(original.coordinates()) {
            near(a, b, original.kinetic_energy_j().sqrt(), 2e-10);
        }
        for k in 0..2 {
            near(
                pressures[k].cosine.value() - pressures[2 + k].cosine.value()
                    + pressures[4 + k].cosine.value(),
                difference[k].cosine.value(),
                40.0,
                1e-14,
            );
            near(
                pressures[k].sine.value() - pressures[2 + k].sine.value()
                    + pressures[4 + k].sine.value(),
                difference[k].sine.value(),
                40.0,
                1e-14,
            );
        }
    });
}

// Independent small-step RK4 of the original complete generator, simultaneous
// physical port powers, positive branch source squares and displaced volumes.
// It does not use the prepared exponential, power cache or endpoint averages.
fn reference(
    network: &LaminarNetwork<FluidTransportState>,
    initial: &[f64],
    pressures: &[LaminarPressureTone],
    duration: f64,
) -> Vec<f64> {
    let (n, m, e) = (initial.len(), network.port_count(), network.branches.len());
    let (_, resistance, g) = network.system.structure();
    let rhs = |time: f64, values: &[f64]| {
        let (sn, cs) = (5.0 * time).sin_cos();
        let mut input: Vec<_> = (0..m)
            .map(|p| {
                pressures[2 * p].cosine.value()
                    + pressures[2 * p + 1].cosine.value() * cs
                    + pressures[2 * p + 1].sine.value() * sn
            })
            .collect();
        let absolute = input.clone();
        for p in 0..network.topology.pressure_nodes.len() {
            input[p] -= absolute[network.basis.pressure_reference[p]];
        }
        let x = &values[..n];
        let mut result = vec![0.0; n + 2 * e + m];
        for i in 0..n {
            result[i] =
                -dot(&resistance[i * n..(i + 1) * n], x) + dot(&g[i * m..(i + 1) * m], &input);
        }
        for edge in 0..e {
            result[n + edge] = network.source_rows[edge]
                .chunks_exact(n)
                .map(|row| dot(row, x).powi(2))
                .sum();
            result[n + e + m + edge] = dot(&network.flow_rows[edge], x);
        }
        for port in 0..m {
            let output: f64 = (0..n).map(|i| g[i * m + port] * x[i]).sum();
            result[n + e + port] = input[port] * output;
        }
        result
    };
    let steps = 4096;
    let h = duration / steps as f64;
    let mut values = vec![0.0; n + 2 * e + m];
    values[..n].copy_from_slice(initial);
    let shifted = |x: &[f64], k: &[f64], dt: f64| -> Vec<f64> {
        x.iter().zip(k).map(|(x, k)| x + dt * k).collect()
    };
    for step in 0..steps {
        let t = step as f64 * h;
        let a = rhs(t, &values);
        let b = rhs(t + h / 2.0, &shifted(&values, &a, h / 2.0));
        let c = rhs(t + h / 2.0, &shifted(&values, &b, h / 2.0));
        let d = rhs(t + h, &shifted(&values, &c, h));
        for i in 0..values.len() {
            values[i] += h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    values
}

#[test]
fn fork_matches_independent_trajectory_each_heat_and_each_source_work() {
    with_cx(|cx| {
        let network = fork(cx);
        let mut state = network.rest_state(0.0).unwrap();
        let mut pressures = source();
        let expected = reference(&network, state.coordinates(), &pressures, 0.02);
        let mut prepared = network
            .prepare_tones(cx, &[0.0, 5.0], 0.02, budget())
            .unwrap();
        let result = prepared.advance(cx, &mut state, &mut pressures).unwrap();
        let n = state.coordinates().len();
        for (&a, &b) in state.coordinates().iter().zip(&expected[..n]) {
            near(a, b, state.kinetic_energy_j().sqrt(), 2e-8);
        }
        for (edge, branch) in prepared.branch_steps().iter().enumerate() {
            near(
                branch.viscous_heat_j,
                expected[n + edge],
                result.native.dissipated,
                2e-8,
            );
            near(
                branch.displaced_volume_m3,
                expected[n + 3 + 6 + edge],
                1e-11,
                2e-8,
            );
            assert!(branch.viscous_heat_j >= 0.0);
        }
        for (port, &work) in prepared.port_work_j().iter().enumerate() {
            near(
                work,
                expected[n + 3 + port],
                result.native.supplied.abs(),
                2e-8,
            );
        }
        near(result.maximum_internal_imbalance_m3_s, 0.0, 1e-9, 1e-12);
        near(
            result.maximum_component_boundary_imbalance_m3_s,
            0.0,
            1e-9,
            1e-12,
        );
        near(
            result.native.balance_residual(),
            0.0,
            result.native.supplied.abs(),
            2e-10,
        );
        near(
            result.reservoir_work_j + result.branch_head_work_j,
            result.native.supplied,
            result.native.supplied.abs(),
            2e-10,
        );
    });
}

#[test]
fn subdivision_preserves_coupled_history_phase_volume_and_all_accounts() {
    with_cx(|cx| {
        let network = fork(cx);
        let mut whole_state = network.rest_state(2.0).unwrap();
        let mut divided_state = whole_state.clone();
        let mut whole_source = source();
        let mut divided_source = whole_source;
        let mut whole = network
            .prepare_tones(cx, &[0.0, 5.0], 0.02, budget())
            .unwrap();
        let mut divided = network
            .prepare_tones(cx, &[0.0, 5.0], 0.005, budget())
            .unwrap();
        let report = whole
            .advance(cx, &mut whole_state, &mut whole_source)
            .unwrap();
        let (mut heat, mut volume, mut work) = ([0.0; 3], [0.0; 3], [0.0; 6]);
        for _ in 0..4 {
            divided
                .advance(cx, &mut divided_state, &mut divided_source)
                .unwrap();
            for (edge, branch) in divided.branch_steps().iter().enumerate() {
                heat[edge] += branch.viscous_heat_j;
                volume[edge] += branch.displaced_volume_m3;
            }
            for (sum, value) in work.iter_mut().zip(divided.port_work_j()) {
                *sum += value;
            }
        }
        for (&a, &b) in whole_state
            .coordinates()
            .iter()
            .zip(divided_state.coordinates())
        {
            near(a, b, whole_state.kinetic_energy_j().sqrt(), 2e-9);
        }
        for edge in 0..3 {
            near(
                heat[edge],
                whole.branch_steps()[edge].viscous_heat_j,
                report.native.dissipated,
                2e-8,
            );
            near(
                volume[edge],
                whole.branch_steps()[edge].displaced_volume_m3,
                1e-11,
                2e-9,
            );
        }
        for (a, b) in work.iter().zip(whole.port_work_j()) {
            near(*a, *b, report.native.supplied.abs(), 2e-8);
        }
        for (a, b) in whole_source.iter().zip(&divided_source) {
            near(a.cosine.value(), b.cosine.value(), 30.0, 1e-14);
            near(a.sine.value(), b.sine.value(), 30.0, 1e-14);
        }
    });
}

#[test]
fn dc_recovers_held_network_and_component_gauges_do_not_drive_flow() {
    with_cx(|cx| {
        let network = fork(cx);
        let mut held = network.rest_state(0.0).unwrap();
        let mut continuous = held.clone();
        let pressures = [Pressure::new(20.0), Pressure::new(-5.0), Pressure::new(2.0)];
        let heads = [Pressure::new(8.0), Pressure::new(-12.0), Pressure::new(4.0)];
        let reference = network
            .prepare(cx, 0.02)
            .unwrap()
            .advance_with_heads(cx, &mut held, &pressures, &heads)
            .unwrap();
        let mut tones: Vec<_> = pressures
            .iter()
            .chain(&heads)
            .map(|p| tone(p.value(), 0.0))
            .collect();
        let mut prepared = network.prepare_tones(cx, &[0.0], 0.02, budget()).unwrap();
        let result = prepared.advance(cx, &mut continuous, &mut tones).unwrap();
        for (&a, &b) in held.coordinates().iter().zip(continuous.coordinates()) {
            near(a, b, held.kinetic_energy_j().sqrt(), 2e-9);
        }
        near(
            result.native.dissipated,
            reference.native.dissipated,
            reference.native.dissipated,
            2e-8,
        );
        near(
            result.reservoir_work_j,
            reference.reservoir_work_j,
            reference.native.supplied.abs(),
            2e-8,
        );
        near(
            result.branch_head_work_j,
            reference.branch_head_work_j,
            reference.native.supplied.abs(),
            2e-8,
        );

        let disconnected = LaminarNetwork::try_new(
            cx,
            LaminarNetworkTopology {
                node_count: 4,
                branch_nodes: vec![[0, 1], [2, 3]],
                pressure_nodes: vec![0, 1, 2, 3],
            },
            vec![pipe(cx, 0.001), pipe(cx, 0.002)],
        )
        .unwrap();
        let mut rest = disconnected.rest_state(0.0).unwrap();
        let mut gauge = [
            tone(1e8, 0.0),
            tone(2e7, -3e7),
            tone(1e8, 0.0),
            tone(2e7, -3e7),
            tone(-2e8, 0.0),
            tone(-4e7, 5e7),
            tone(-2e8, 0.0),
            tone(-4e7, 5e7),
            tone(0.0, 0.0),
            tone(0.0, 0.0),
            tone(0.0, 0.0),
            tone(0.0, 0.0),
        ];
        let mut prepared = disconnected
            .prepare_tones(cx, &[0.0, 5.0], 0.02, budget())
            .unwrap();
        let report = prepared.advance(cx, &mut rest, &mut gauge).unwrap();
        assert!(rest.coordinates().iter().all(|v| *v == 0.0));
        assert_eq!(report.native.dissipated, 0.0);
        assert!(prepared.port_work_j().iter().all(|v| *v == 0.0));
    });
}

#[test]
fn hidden_pressure_peaks_and_bad_history_leave_all_accepted_outputs_untouched() {
    with_cx(|cx| {
        let network = fork(cx);
        let mut state = network.rest_state(0.0).unwrap();
        let mut pressures = source();
        let mut prepared = network
            .prepare_tones(cx, &[0.0, 5.0], 0.02, budget())
            .unwrap();
        prepared.advance(cx, &mut state, &mut pressures).unwrap();
        let accepted = state.clone();
        let work = prepared.port_work_j().to_vec();
        let flow = prepared.boundary_volume_flows_m3_s().to_vec();
        let heat: Vec<_> = prepared
            .branch_steps()
            .iter()
            .map(|b| b.viscous_heat_j)
            .collect();
        for invalid in [tone(0.0, 1e12), tone(f64::NAN, 0.0)] {
            let mut bad = pressures;
            bad[7] = invalid;
            let before: Vec<_> = bad
                .iter()
                .map(|p| (p.cosine.value().to_bits(), p.sine.value().to_bits()))
                .collect();
            assert!(prepared.advance(cx, &mut state, &mut bad).is_err());
            assert_eq!(state.coordinates(), accepted.coordinates());
            assert_eq!(state.accepted_time_s(), accepted.accepted_time_s());
            assert_eq!(
                bad.iter()
                    .map(|p| (p.cosine.value().to_bits(), p.sine.value().to_bits()))
                    .collect::<Vec<_>>(),
                before
            );
            assert_eq!(prepared.port_work_j(), work);
            assert_eq!(prepared.boundary_volume_flows_m3_s(), flow);
            assert_eq!(
                prepared
                    .branch_steps()
                    .iter()
                    .map(|b| b.viscous_heat_j)
                    .collect::<Vec<_>>(),
                heat
            );
        }
        let foreign = fork(cx);
        let mut foreign_state = foreign.rest_state(0.0).unwrap();
        assert!(
            prepared
                .advance(cx, &mut foreign_state, &mut pressures)
                .is_err()
        );
        let mut bad_dc = pressures;
        bad_dc[0] = tone(20.0, 1.0);
        assert!(prepared.advance(cx, &mut state, &mut bad_dc).is_err());
        assert!(
            prepared
                .advance(cx, &mut state, &mut pressures[..11])
                .is_err()
        );
        prepared.advance(cx, &mut state, &mut pressures).unwrap();
        assert!(state.accepted_time_s() > accepted.accepted_time_s());
        assert!(network.prepare_tones(cx, &[10.1], 0.02, budget()).is_err());
        assert!(
            network
                .prepare_tones(
                    cx,
                    &[0.0, 5.0],
                    0.02,
                    DenseExponentialBudget {
                        max_dense_bytes: 1,
                        ..budget()
                    }
                )
                .is_err()
        );
    });
}
