use super::*;
use crate::graph::inertial::laminar::LaminarPipeSpec;
use crate::graph::tests::with_cx;
use fs_material::fluid_transport::FluidTransportState;
use fs_phs::ZkDiffusiveBasis;
use fs_qty::Length;

fn pipe(cx: &Cx<'_>, length: f64, viscosity: f64) -> LaminarPipe<FluidTransportState> {
    pipe_with_reynolds(cx, length, viscosity, 1000.0)
}

fn pipe_with_reynolds(
    cx: &Cx<'_>,
    length: f64,
    viscosity: f64,
    maximum_reynolds: f64,
) -> LaminarPipe<FluidTransportState> {
    LaminarPipe::try_new(
        cx,
        FluidTransportState::try_new(1000.0, 4200.0, 0.6, viscosity).unwrap(),
        LaminarPipeSpec {
            length: Length::new(length),
            radius: Length::new(0.001),
            maximum_reynolds,
            band_rad_s: [0.01, 900.0],
            maximum_relative_wall_error: 1.0e-4,
        },
        &ZkDiffusiveBasis::new(16).unwrap(),
    )
    .unwrap()
}

fn topology(nodes: usize, edges: &[[usize; 2]], ports: &[usize]) -> LaminarNetworkTopology {
    LaminarNetworkTopology {
        node_count: nodes,
        branch_nodes: edges.to_vec(),
        pressure_nodes: ports.to_vec(),
    }
}

fn budget() -> HarmonicBudget {
    HarmonicBudget {
        max_dimension: 128,
        max_dense_bytes: 64 << 20,
        max_dense_work: 200_000_000,
        relative_residual_tolerance: 1.0e-10,
    }
}

fn near(actual: C64, expected: C64) {
    let scale = actual.abs() + expected.abs();
    assert!((actual - expected).abs() <= 3.0e-9 * scale.max(1.0e-300));
}

#[test]
fn single_and_series_networks_retain_original_complex_impedance() {
    with_cx(|cx| {
        let omega = 17.0;
        let first = pipe(cx, 0.15, 0.001);
        let z0 = first.impedance_at(omega).unwrap().conj();
        let single = LaminarNetwork::try_new(
            cx,
            topology(2, &[[0, 1]], &[0, 1]),
            vec![first],
        )
        .unwrap();
        let pressure = C64::new(80.0, -25.0);
        let response = single
            .harmonic_response(cx, omega, &[pressure, C64::ZERO], budget())
            .unwrap();
        near(response.output[2], pressure / z0);
        near(response.output[0], response.output[2]);
        near(response.output[1], response.output[2].scale(-1.0));

        let a = pipe(cx, 0.15, 0.001);
        let b = pipe(cx, 0.25, 0.003);
        let impedance =
            a.impedance_at(omega).unwrap().conj() + b.impedance_at(omega).unwrap().conj();
        let series = LaminarNetwork::try_new(
            cx,
            topology(3, &[[0, 1], [1, 2]], &[0, 2]),
            vec![a, b],
        )
        .unwrap();
        let response = series
            .harmonic_response(cx, omega, &[pressure, C64::ZERO], budget())
            .unwrap();
        near(response.output[2], pressure / impedance);
        near(response.output[3], response.output[2]);
        assert!(
            response.power_residual_w.abs()
                < 2.0e-9 * response.supplied_power_w.abs().max(1.0e-300)
        );
    });
}

#[test]
fn explicit_head_drives_a_closed_loop_without_reservoir_flow() {
    with_cx(|cx| {
        let omega = 9.0;
        let a = pipe(cx, 0.10, 0.001);
        let b = pipe(cx, 0.20, 0.002);
        let c = pipe(cx, 0.12, 0.0015);
        let impedance = a.impedance_at(omega).unwrap().conj()
            + b.impedance_at(omega).unwrap().conj()
            + c.impedance_at(omega).unwrap().conj();
        let network = LaminarNetwork::try_new(
            cx,
            topology(3, &[[0, 1], [1, 2], [2, 0]], &[0]),
            vec![a, b, c],
        )
        .unwrap();
        let head = C64::new(120.0, 30.0);
        let response = network
            .harmonic_response_with_heads(
                cx,
                omega,
                &[C64::new(1.0e9, -2.0e8)],
                &[head, C64::ZERO, C64::ZERO],
                budget(),
            )
            .unwrap();
        assert!(response.output[0].abs() < 1.0e-18);
        for flow in &response.output[1..] {
            near(*flow, head / impedance);
        }
        assert!(response.supplied_power_w > 0.0);
    });
}

#[test]
fn batch_reuses_one_factor_and_matches_individual_checked_responses() {
    with_cx(|cx| {
        let network = LaminarNetwork::try_new(
            cx,
            topology(3, &[[0, 1], [1, 2]], &[0, 2]),
            vec![pipe(cx, 0.15, 0.001), pipe(cx, 0.25, 0.002)],
        )
        .unwrap();
        let cases = vec![
            vec![C64::new(20.0, -4.0), C64::ZERO],
            vec![
                C64::new(1.0e8 - 8.0, -3.0e7 + 13.0),
                C64::new(1.0e8, -3.0e7),
            ],
        ];
        let batch = network
            .harmonic_responses_batch(cx, 23.0, &cases, budget())
            .unwrap();
        assert_eq!(batch.len(), cases.len());
        for (actual, input) in batch.iter().zip(&cases) {
            let individual = network
                .harmonic_response(cx, 23.0, input, budget())
                .unwrap();
            assert_eq!(actual, &individual);
        }
    });
}

#[test]
fn harmonic_peak_reynolds_is_admitted_and_batch_refusal_is_atomic() {
    with_cx(|cx| {
        let network = LaminarNetwork::try_new(
            cx,
            topology(2, &[[0, 1]], &[0, 1]),
            vec![pipe_with_reynolds(cx, 0.2, 0.001, 10.0)],
        )
        .unwrap();
        let gentle = vec![C64::new(0.01, 0.0), C64::ZERO];
        assert!(
            network
                .harmonic_response(cx, 10.0, &gentle, budget())
                .is_ok()
        );
        let excessive = vec![C64::new(1000.0, 0.0), C64::ZERO];
        assert!(
            network
                .harmonic_response(cx, 10.0, &excessive, budget())
                .is_err()
        );
        let cases = vec![gentle, excessive];
        assert!(
            network
                .harmonic_responses_batch(cx, 10.0, &cases, budget())
                .is_err()
        );
    });
}

#[test]
fn malformed_ports_batches_and_out_of_band_frequency_refuse() {
    with_cx(|cx| {
        let network = LaminarNetwork::try_new(
            cx,
            topology(2, &[[0, 1]], &[0, 1]),
            vec![pipe(cx, 0.2, 0.001)],
        )
        .unwrap();
        assert!(
            network
                .harmonic_response(cx, 10.0, &[C64::ZERO], budget())
                .is_err()
        );
        assert!(
            network
                .harmonic_response_with_heads(
                    cx,
                    10.0,
                    &[C64::ZERO; 2],
                    &[C64::ZERO; 2],
                    budget(),
                )
                .is_err()
        );
        assert!(
            network
                .harmonic_response(cx, 1000.0, &[C64::ONE, C64::ZERO], budget())
                .is_err()
        );
        assert!(
            network
                .harmonic_responses_batch(cx, 10.0, &[], budget())
                .is_err()
        );
        let pressure_cases = vec![vec![C64::ZERO; 2], vec![C64::ZERO; 2]];
        let head_cases = vec![vec![C64::ZERO; 1]];
        assert!(
            network
                .harmonic_responses_with_heads_batch(
                    cx,
                    10.0,
                    &pressure_cases,
                    &head_cases,
                    budget(),
                )
                .is_err()
        );
    });
}
