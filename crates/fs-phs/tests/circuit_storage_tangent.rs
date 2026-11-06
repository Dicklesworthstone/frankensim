//! Exact continuous storage tangents connect biased nonlinear circuits to
//! existing stationary and harmonic analysis, without finite-difference laws.
use fs_math::c64::C64;
use fs_phs::circuit::{
    Branch, CircuitDae, CircuitGraph, MutualInductance,
    capacitance::{CapacitanceBlock, CapacitorChargeCurve, CapacitorCurve},
    inductance::{
        InductorCurve, InductorFluxCurve, assemble_circuit_with_inductor_curves,
        assemble_circuit_with_storage_curves,
    },
};
use fs_phs::harmonic::{
    HarmonicBudget,
    equilibrium::{
        EquilibriumLinearization, EquilibriumPoint, LinearizationBudget, LinearizationError,
    },
};

fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}
fn curve(scale: f64) -> InductorFluxCurve {
    InductorFluxCurve::try_new(
        [
            (-2.0, -6.0),
            (-1.0, -1.0),
            (0.0, 0.0),
            (1.0, 1.0),
            (2.0, 6.0),
        ]
        .map(|(flux, current)| (scale * flux, current))
        .to_vec(),
    )
    .unwrap()
}
fn rl(scale: f64) -> CircuitDae {
    assemble_circuit_with_inductor_curves(
        &CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, Branch::VoltageSource { port: 0 }),
                (1, 2, Branch::Resistor { ohms: 1.0 }),
                (2, 0, Branch::Inductor { henries: scale }),
            ],
            transformers: vec![],
        },
        &[InductorCurve {
            branch: 2,
            curve: curve(scale),
        }],
        &[],
    )
    .unwrap()
}
fn budget() -> LinearizationBudget {
    LinearizationBudget {
        harmonic: HarmonicBudget {
            max_dimension: 32,
            max_dense_bytes: 4_000_000,
            max_dense_work: 10_000_000,
            relative_residual_tolerance: 1e-11,
        },
        max_hessian_actions: 32,
    }
}

#[test]
fn biased_saturating_rl_ac_uses_differential_inductance_at_both_polarities() {
    for scale in [1.0, 1e-6] {
        let source = rl(scale);
        let storage = source.system.storage();
        let tangent = storage
            .continuous_tangent()
            .expect("original continuous law");
        for bias in [-1.5_f64, -0.5, 0.0, 0.5, 1.5] {
            let mut state = vec![0.0; source.system.state_dim()];
            state[source.flux_index[0]] = scale * bias;
            let input = [curve(scale).current_a(scale * bias).unwrap()];
            let seed = state.clone();
            state = source.consistent_initial_state(&seed, &input).unwrap();
            // Native initialization returns a discrete-gradient endpoint.
            // Continuous equilibrium takes the physical algebraic efforts,
            // which are the initializer's seed/endpoint midpoint.
            for i in source.system.differential_dim()..state.len() {
                state[i] = f64::midpoint(seed[i], state[i]);
            }
            let held = state.clone();
            let scales = vec![1.0; state.len()];
            let image = EquilibriumLinearization::from_descriptor(
                &source.system,
                EquilibriumPoint {
                    state: &state,
                    input: &input,
                    balance_scales: &scales,
                    relative_tolerance: 1e-10,
                },
                budget(),
                &|x, d, out| tangent.apply(x, d, out),
                || false,
            )
            .unwrap();
            let l = scale / if bias.abs() > 1.0 { 5.0 } else { 1.0 };
            let index = source.flux_index[0];
            near(
                image.tangent_storage_matrix()[index * state.len() + index],
                1.0 / l,
                2e-9,
            );
            for normalized_omega in [0.1, 1.0, 3.5] {
                let omega = normalized_omega / scale;
                let drive = C64::new(0.7, -0.2);
                let response = image
                    .model()
                    .factor(omega, || false)
                    .unwrap()
                    .response(&[drive], || false)
                    .unwrap();
                // exp(-i omega t): delta V=(R-i omega L_diff)*delta I.
                let current = drive / C64::new(1.0, -omega * l);
                assert!((response.effort[index] - current).abs() < 3e-11);
                assert!((response.state[index] - current.scale(l)).abs() < 3e-11 * scale);
                assert!((response.output[0] - current).abs() < 3e-11);
                let heat = 0.5 * current.abs().powi(2);
                near(response.dissipated_power_w, heat, 3e-11);
                near(response.supplied_power_w, heat, 3e-11);
                near(
                    response.mean_stored_energy_j,
                    0.25 * l * current.abs().powi(2),
                    3e-11 * scale,
                );
            }
            assert_eq!(state, held);
        }
    }
}

#[test]
fn continuous_actions_keep_nonlinear_l_c_mutual_blocks_and_algebraic_zero_storage() {
    let graph = CircuitGraph {
        node_count: 7,
        branches: vec![
            (1, 0, Branch::Inductor { henries: 1.0 }),
            (1, 0, Branch::Capacitor { farads: 1.0 }),
            (2, 0, Branch::Inductor { henries: 2.0 }),
            (3, 0, Branch::Inductor { henries: 3.0 }),
            (4, 0, Branch::Capacitor { farads: 1.0 }),
            (5, 0, Branch::Capacitor { farads: 1.0 }),
            (6, 0, Branch::Capacitor { farads: 1.0 }),
        ],
        transformers: vec![],
    };
    let source = assemble_circuit_with_storage_curves(
        &graph,
        &[CapacitanceBlock::try_new(vec![4, 5], vec![2.0, 1.0, 1.0, 3.0]).unwrap()],
        &[CapacitorCurve {
            branch: 1,
            curve: CapacitorChargeCurve::try_new(curve(1.0).knots().to_vec()).unwrap(),
        }],
        &[InductorCurve {
            branch: 0,
            curve: curve(1.0),
        }],
        &[MutualInductance {
            branches: (2, 3),
            henries: 1.0,
        }],
    )
    .unwrap();
    let n = source.system.state_dim();
    let mut state = vec![0.0; n];
    state[source.flux_index[0]] = 1.5;
    state[source.charge_index[0]] = -0.25;
    state[source.charge_index[3]] = 1e8;
    let mut direction = vec![0.0; n];
    for (index, value) in [
        (source.flux_index[0], 0.3),
        (source.flux_index[1], 0.2),
        (source.flux_index[2], -0.4),
        (source.charge_index[0], -0.7),
        (source.charge_index[1], 0.4),
        (source.charge_index[2], 0.1),
    ] {
        direction[index] = value;
    }
    let mut expected = vec![0.0; n];
    for (index, value) in [
        (source.flux_index[0], 1.5),
        (source.flux_index[1], 0.2),
        (source.flux_index[2], -0.2),
        (source.charge_index[0], -0.7),
        (source.charge_index[1], 0.22),
        (source.charge_index[2], -0.04),
    ] {
        expected[index] = value;
    }
    let storage = source.system.storage();
    let shared = storage.shared_clone().unwrap();
    for owner in [storage, shared.as_ref()] {
        let tangent = owner.continuous_tangent().unwrap();
        let mut actual = vec![f64::NAN; n];
        assert!(tangent.apply(&state, &direction, &mut actual));
        for (&a, &b) in actual.iter().zip(&expected) {
            near(a, b, 3e-15);
        }
        // Independent differentiation of the original gradient, away from kinks.
        let epsilon = 1e-5;
        let plus: Vec<_> = state
            .iter()
            .zip(&direction)
            .map(|(x, d)| x + epsilon * d)
            .collect();
        let minus: Vec<_> = state
            .iter()
            .zip(&direction)
            .map(|(x, d)| x - epsilon * d)
            .collect();
        let (mut ep, mut em) = (vec![0.0; n], vec![0.0; n]);
        owner.gradient(&plus, &mut ep);
        owner.gradient(&minus, &mut em);
        for ((&p, &m), &a) in ep.iter().zip(&em).zip(&actual) {
            near((p - m) / (2.0 * epsilon), a, 5e-10);
        }
        assert!(
            actual[source.system.differential_dim()..]
                .iter()
                .all(|&x| x == 0.0)
        );
    }
}

#[test]
fn source_edges_and_genuine_kinks_refuse_local_ac_while_original_transients_remain_valid() {
    let source = rl(1.0);
    let tangent = source.system.storage().continuous_tangent().unwrap();
    let n = source.system.state_dim();
    let mut state = vec![0.0; n];
    let direction = vec![0.0; n];
    let mut out = vec![0.0; n];
    for bias in [-2.0, -1.0, 1.0, 2.0, 2.01, f64::NAN] {
        state[source.flux_index[0]] = bias;
        // Refuse even if this requested direction happens to vanish: the
        // complete continuous derivative of the source does not exist here.
        assert!(!tangent.apply(&state, &direction, &mut out));
    }
    state[source.flux_index[0]] = 0.0;
    assert!(tangent.apply(&state, &direction, &mut out));
    assert!(!tangent.apply(&state[..n - 1], &direction, &mut out));
    assert!(!tangent.apply(&state, &direction[..n - 1], &mut out));
    assert!(!tangent.apply(&state, &direction, &mut out[..n - 1]));
    let mut nonfinite = direction.clone();
    nonfinite[0] = f64::INFINITY;
    assert!(!tangent.apply(&state, &nonfinite, &mut out));

    // Time stepping retains its valid piecewise integral and one-sided Newton
    // representative at the same kink; this feature does not alter that owner.
    state[source.flux_index[0]] = 1.0;
    let seed = state.clone();
    let mut state = source.consistent_initial_state(&seed, &[1.0]).unwrap();
    for i in source.system.differential_dim()..state.len() {
        state[i] = f64::midpoint(seed[i], state[i]);
    }
    assert!(source.step_audited(&state, &[1.2], 0.01).is_ok());
    let scales = vec![1.0; n];
    assert!(matches!(
        EquilibriumLinearization::from_descriptor(
            &source.system,
            EquilibriumPoint {
                state: &state,
                input: &[1.0],
                balance_scales: &scales,
                relative_tolerance: 1e-10
            },
            budget(),
            &|x, d, out| tangent.apply(x, d, out),
            || false,
        ),
        Err(LinearizationError::HessianRefused { .. })
    ));
}

#[test]
fn collinear_source_knots_admit_a_genuine_two_sided_continuous_derivative() {
    let source = assemble_circuit_with_inductor_curves(
        &CircuitGraph {
            node_count: 2,
            branches: vec![(1, 0, Branch::Inductor { henries: 2.0 })],
            transformers: vec![],
        },
        &[InductorCurve {
            branch: 0,
            curve: InductorFluxCurve::try_new(vec![
                (-4.0, -2.0),
                (-2.0, -1.0),
                (0.0, 0.0),
                (2.0, 1.0),
                (4.0, 2.0),
            ])
            .unwrap(),
        }],
        &[],
    )
    .unwrap();
    let n = source.system.state_dim();
    let tangent = source.system.storage().continuous_tangent().unwrap();
    let mut state = vec![0.0; n];
    let mut direction = vec![0.0; n];
    direction[source.flux_index[0]] = 3.0;
    let mut out = vec![0.0; n];
    for bias in [-2.0, 0.0, 2.0] {
        state[source.flux_index[0]] = bias;
        assert!(tangent.apply(&state, &direction, &mut out));
        assert_eq!(out[source.flux_index[0]], 1.5);
    }
}
