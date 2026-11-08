use super::*;
use crate::graph::inertial::laminar::LaminarPipeSpec;
use crate::graph::tests::with_cx;
use fs_material::fluid_transport::FluidTransportState;
use fs_math::det;
use fs_phs::ZkDiffusiveBasis;
use fs_phs::linear_exact::DenseExponentialBudget;
use fs_qty::Length;

fn pipe(cx: &Cx<'_>, viscosity: f64, maximum_reynolds: f64) -> LaminarPipe<FluidTransportState> {
    LaminarPipe::try_new(
        cx,
        FluidTransportState::try_new(1000.0, 4200.0, 0.6, viscosity).unwrap(),
        LaminarPipeSpec {
            length: Length::new(0.4),
            radius: Length::new(0.0005),
            maximum_reynolds,
            band_rad_s: [0.01, 10.0],
            maximum_relative_wall_error: 1.0e-4,
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
        vec![
            pipe(cx, 0.001, 1000.0),
            pipe(cx, 0.002, 1000.0),
            pipe(cx, 0.003, 1000.0),
        ],
    )
    .unwrap()
}

fn budget() -> HarmonicBudget {
    HarmonicBudget {
        max_dimension: 128,
        max_dense_bytes: 64 << 20,
        max_dense_work: 200_000_000,
        relative_residual_tolerance: 1.0e-10,
    }
}

fn tone(cosine: f64, sine: f64) -> LaminarPressureTone {
    LaminarPressureTone {
        cosine: Pressure::new(cosine),
        sine: Pressure::new(sine),
    }
}

#[test]
fn multifrequency_settled_state_advances_without_a_free_transient() {
    with_cx(|cx| {
        let network = fork(cx);
        let frequencies = [0.0, 3.0, det::sqrt(17.0)];
        let mut sources = vec![
            tone(20.0, 0.0),
            tone(10.0, 2.0),
            tone(-4.0, 5.0),
            tone(-5.0, 0.0),
            tone(-3.0, 7.0),
            tone(2.0, -1.0),
            tone(2.0, 0.0),
            tone(5.0, -4.0),
            tone(6.0, 3.0),
            tone(8.0, 0.0),
            tone(1.0, 0.0),
            tone(2.0, -5.0),
            tone(-12.0, 0.0),
            tone(0.0, -6.0),
            tone(-3.0, 1.0),
            tone(4.0, 0.0),
            tone(-2.0, 3.0),
            tone(1.0, 4.0),
        ];
        let mut state = network
            .settled_tone_state(cx, &frequencies, &sources, 12.25, budget())
            .unwrap();
        assert!(
            state.coordinates()[network.basis.columns..]
                .iter()
                .any(|v| *v != 0.0)
        );
        let mut prepared = network
            .prepare_tones(
                cx,
                &frequencies,
                0.08,
                DenseExponentialBudget {
                    max_dimension: 256,
                    max_dense_bytes: 128 << 20,
                    max_multiply_terms: 10_000_000_000,
                },
            )
            .unwrap();
        for _ in 0..2 {
            prepared.advance(cx, &mut state, &mut sources).unwrap();
            let expected = network
                .settled_tone_state(
                    cx,
                    &frequencies,
                    &sources,
                    state.accepted_time_s(),
                    budget(),
                )
                .unwrap();
            let difference: Vec<_> = state
                .coordinates()
                .iter()
                .zip(expected.coordinates())
                .map(|(actual, target)| actual - target)
                .collect();
            let scale = norm(state.coordinates()).max(norm(expected.coordinates()));
            assert!(
                norm(&difference) <= 2.0e-8 * scale.max(1.0e-300),
                "settled phase mismatch: {} / {scale}",
                norm(&difference),
            );
        }
    });
}

#[test]
fn dc_matches_cached_steady_state_and_invalid_authored_sources_refuse() {
    with_cx(|cx| {
        let network = fork(cx);
        let sources = [
            tone(20.0, 0.0),
            tone(-5.0, 0.0),
            tone(2.0, 0.0),
            tone(8.0, 0.0),
            tone(-12.0, 0.0),
            tone(4.0, 0.0),
        ];
        let expected = network
            .steady_state_with_heads(
                &sources[..3].iter().map(|s| s.cosine).collect::<Vec<_>>(),
                &sources[3..].iter().map(|s| s.cosine).collect::<Vec<_>>(),
                9.0,
            )
            .unwrap();
        let state = network
            .settled_tone_state(cx, &[0.0], &sources, 9.0, budget())
            .unwrap();
        assert_eq!(state.coordinates(), expected.coordinates());
        assert_eq!(state.accepted_time_s(), 9.0);
        assert!(network.branch_state(&state, 2).is_ok());

        let mut invalid_dc = sources;
        invalid_dc[4].sine = Pressure::new(1.0);
        assert!(
            network
                .settled_tone_state(cx, &[0.0], &invalid_dc, 9.0, budget())
                .is_err()
        );
        assert!(
            network
                .settled_tone_state(cx, &[0.0], &sources[..5], 9.0, budget())
                .is_err()
        );
        assert!(
            network
                .settled_tone_state(cx, &[0.0], &sources, -1.0, budget())
                .is_err()
        );
        for frequency in [f64::NAN, -1.0, 11.0] {
            assert!(
                network
                    .settled_tone_state(cx, &[frequency], &sources, 0.0, budget())
                    .is_err()
            );
        }
    });
}

#[test]
fn combined_dc_and_tone_peak_envelope_enforces_original_reynolds_ceiling() {
    with_cx(|cx| {
        let original = pipe(cx, 0.001, 10.0);
        // Re = 2*rho*|Q|/(pi*radius*mu). Choose two quadrature flows:
        // each has 0.4 of the original ceiling, with only DC flowing at tau=0.
        let limit = 10.0 * core::f64::consts::PI * 0.0005 * 0.001 / (2.0 * 1000.0);
        let frequencies = [0.0, 3.0, det::sqrt(17.0)];
        let dc = 0.3 * limit * original.steady_resistance_pa_s_m3();
        let positive: Vec<_> = frequencies[1..]
            .iter()
            .map(|&omega| original.impedance_at(omega).unwrap().conj() * C64::new(0.0, 0.4 * limit))
            .collect();
        let network = LaminarNetwork::try_new(
            cx,
            LaminarNetworkTopology {
                node_count: 2,
                branch_nodes: vec![[0, 1]],
                pressure_nodes: vec![0, 1],
            },
            vec![original],
        )
        .unwrap();
        let mut sources = vec![tone(0.0, 0.0); 9];
        sources[0] = tone(dc, 0.0);
        for (index, &phasor) in positive.iter().enumerate() {
            sources[index + 1] = tone(phasor.re, phasor.im);
            assert!(
                network
                    .harmonic_response(cx, frequencies[index + 1], &[phasor, C64::ZERO], budget())
                    .is_ok()
            );
        }
        assert!(matches!(
            network.settled_tone_state(cx, &frequencies, &sources, 0.0, budget()),
            Err(HarmonicError::Invalid {
                what: "settled simultaneous branch peaks exceed original laminar admission",
            })
        ));
        // Removing only DC admits the same two oscillatory responses.
        sources[0] = tone(0.0, 0.0);
        assert!(
            network
                .settled_tone_state(cx, &frequencies, &sources, 0.0, budget())
                .is_ok()
        );
    });
}

#[test]
fn independent_component_pressure_gauges_preserve_settled_heads_and_memory() {
    with_cx(|cx| {
        let network = LaminarNetwork::try_new(
            cx,
            LaminarNetworkTopology {
                node_count: 4,
                branch_nodes: vec![[0, 1], [2, 3]],
                pressure_nodes: vec![0, 1, 2, 3],
            },
            vec![pipe(cx, 0.001, 1000.0), pipe(cx, 0.002, 1000.0)],
        )
        .unwrap();
        let frequencies = [0.0, 3.0, det::sqrt(17.0)];
        let mut sources = vec![tone(0.0, 0.0); 18];
        sources[12..].copy_from_slice(&[
            tone(20.0, 0.0),
            tone(5.0, -2.0),
            tone(-3.0, 4.0),
            tone(-10.0, 0.0),
            tone(1.0, 6.0),
            tone(2.0, -1.0),
        ]);
        let reference = network
            .settled_tone_state(cx, &frequencies, &sources, 4.5, budget())
            .unwrap();
        for (port, gauge) in [1.0e100, 1.0e100, -2.0e100, -2.0e100]
            .into_iter()
            .enumerate()
        {
            sources[3 * port] = tone(gauge, 0.0);
            sources[3 * port + 1] = tone(-3.0 * gauge, 4.0 * gauge);
            sources[3 * port + 2] = tone(2.0 * gauge, -5.0 * gauge);
        }
        let shifted = network
            .settled_tone_state(cx, &frequencies, &sources, 4.5, budget())
            .unwrap();
        assert_eq!(shifted.coordinates(), reference.coordinates());
    });
}
