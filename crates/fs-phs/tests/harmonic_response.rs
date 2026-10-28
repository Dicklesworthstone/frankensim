//! G1 analytic resonance/constraint oracles and G3 scaling/coupling checks.
use fs_math::c64::C64;
use fs_phs::circuit::{
    Branch, CircuitGraph, MutualInductance, TransformerLink, assemble_circuit,
    assemble_circuit_with_mutual_inductance,
};
use fs_phs::harmonic::{HarmonicBudget, HarmonicError, HarmonicModel};
use fs_phs::{PortHamiltonian, QuadraticStorage, RelaxationBranch, ScalarLaw, SeparableStorage};

fn budget() -> HarmonicBudget {
    HarmonicBudget {
        max_dimension: 128,
        max_dense_bytes: 2_000_000,
        max_dense_work: 4_000_000,
        relative_residual_tolerance: 1e-12,
    }
}
fn near(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-200),
        "actual {actual:?}, expected {expected:?}"
    );
}

#[test]
fn damped_resonance_has_analytic_amplitude_phase_and_independent_power() {
    let (mass, stiffness, damping) = (2.5_f64, 1600.0_f64, 0.6_f64);
    let system = fs_phs::mass_spring_damper(mass, stiffness, damping).unwrap();
    let model = HarmonicModel::from_ode(&system, budget(), || false).unwrap();
    let omega0 = (stiffness / mass).sqrt();
    let force = C64::new(3.0, -0.7);
    for ratio in [0.25, 0.999, 1.0, 1.001, 2.0] {
        let omega = ratio * omega0;
        let factor = model.factor(omega, || false).unwrap();
        let result = factor.response(&[force], || false).unwrap();
        let displacement = force / C64::new(stiffness - mass * omega * omega, -damping * omega);
        let velocity = C64::new(0.0, -omega) * displacement;
        near(result.state[0], displacement, 2e-12);
        near(result.state[1], velocity.scale(mass), 2e-12);
        near(result.output[0], velocity, 2e-12);
        let heat = 0.5 * damping * velocity.norm_sq();
        let storage = 0.25 * (stiffness * displacement.norm_sq() + mass * velocity.norm_sq());
        assert!((result.dissipated_power_w / heat - 1.0).abs() < 5e-12);
        assert!((result.mean_stored_energy_j / storage - 1.0).abs() < 5e-12);
        assert!(result.power_residual_w.abs() < 2e-12 * heat);
        assert!(result.relative_residual < 1e-13);
        // The LU is reused for a different complex excitation, preserving phase.
        let multiplier = C64::new(-0.4, 0.9);
        let second = factor.response(&[force * multiplier], || false).unwrap();
        near(second.state[0], displacement * multiplier, 2e-12);
    }
}

#[test]
fn hereditary_internal_states_produce_the_original_frequency_dependent_loss() {
    let (mass, equilibrium, branch, tau) = (0.4_f64, 110.0, 300.0, 0.025);
    let system = fs_phs::mass_spring_damper(mass, equilibrium, 0.0)
        .unwrap()
        .with_relaxation_branches(vec![RelaxationBranch {
            projection: vec![1.0, 0.0],
            stiffness: branch,
            relaxation_time_s: tau,
        }])
        .unwrap();
    let model = HarmonicModel::from_ode(&system, budget(), || false).unwrap();
    for omega in [0.1, 10.0, 30.0, 100.0, 1000.0] {
        let result = model
            .factor(omega, || false)
            .unwrap()
            .response(&[C64::ONE], || false)
            .unwrap();
        let st = C64::new(0.0, -omega * tau);
        let dynamic = C64::from_re(equilibrium - mass * omega * omega)
            + st / (C64::ONE + st).scale(1.0 / branch);
        let expected = C64::ONE / dynamic;
        near(result.state[0], expected, 4e-12);
        near(
            result.state[2],
            expected.scale(branch.sqrt()) / (C64::ONE + st),
            4e-12,
        );
        let loss =
            0.5 * branch * tau * omega * omega / (1.0 + (omega * tau).powi(2)) * expected.norm_sq();
        assert!((result.dissipated_power_w / loss - 1.0).abs() < 5e-11);
        assert!(result.power_residual_w.abs() < 5e-11 * loss);
    }
}

fn series_rlc(resistance: f64, inductance: f64, capacitance: f64) -> CircuitGraph {
    CircuitGraph {
        node_count: 4,
        branches: vec![
            (1, 0, Branch::VoltageSource { port: 0 }),
            (
                1,
                2,
                Branch::Inductor {
                    henries: inductance,
                },
            ),
            (
                2,
                3,
                Branch::Capacitor {
                    farads: capacitance,
                },
            ),
            (3, 0, Branch::Resistor { ohms: resistance }),
        ],
        transformers: vec![],
    }
}

#[test]
fn rlc_resonance_retains_kirchhoff_constraints_and_source_current() {
    let (r, l, c) = (0.2_f64, 0.1_f64, 2.2e-7_f64);
    let circuit = assemble_circuit(&series_rlc(r, l, c)).unwrap();
    let model = HarmonicModel::from_descriptor(&circuit.system, budget(), || false).unwrap();
    let voltage = C64::new(10.0, -3.0);
    for ratio in [0.1, 0.999, 1.0, 1.001, 3.0] {
        let omega = ratio / (l * c).sqrt();
        let result = model
            .factor(omega, || false)
            .unwrap()
            .response(&[voltage], || false)
            .unwrap();
        let current = voltage / C64::new(r, 1.0 / (omega * c) - omega * l);
        near(result.output[0], current, 2e-9);
        near(
            result.state[circuit.flux_index[0]].scale(1.0 / l),
            current,
            2e-9,
        );
        near(
            result.state[circuit.node_potential_index[0]],
            voltage,
            2e-12,
        );
        near(
            result.state[circuit.node_potential_index[2]],
            current.scale(r),
            2e-9,
        );
        let capacitor_voltage = current / C64::new(0.0, -omega * c);
        near(
            result.state[circuit.charge_index[0]].scale(1.0 / c),
            capacitor_voltage,
            2e-9,
        );
        let heat = 0.5 * r * current.norm_sq();
        assert!((result.dissipated_power_w / heat - 1.0).abs() < 4e-9);
        assert!(result.power_residual_w.abs() < 3e-11 * heat);
    }
}

#[test]
fn mutual_inductance_splits_resonance_and_preserves_winding_polarity() {
    let (r, l, c, mutual) = (2.0_f64, 0.1_f64, 1e-6_f64, 0.06_f64);
    let mut graph = series_rlc(r, l, c);
    graph.node_count = 7;
    graph.branches.extend([
        (4, 0, Branch::VoltageSource { port: 1 }),
        (4, 5, Branch::Inductor { henries: l }),
        (5, 6, Branch::Capacitor { farads: c }),
        (6, 0, Branch::Resistor { ohms: r }),
    ]);
    for winding_sign in [-1.0, 1.0] {
        let m = mutual * winding_sign;
        let circuit = assemble_circuit_with_mutual_inductance(
            &graph,
            &[MutualInductance {
                branches: (1, 5),
                henries: m,
            }],
        )
        .unwrap();
        let model = HarmonicModel::from_descriptor(&circuit.system, budget(), || false).unwrap();
        for symmetry in [-1.0, 1.0] {
            let omega = 1.0 / ((l + symmetry * m) * c).sqrt();
            let voltage = C64::new(1.0, 0.3);
            let result = model
                .factor(omega, || false)
                .unwrap()
                .response(&[voltage, voltage.scale(symmetry)], || false)
                .unwrap();
            // Each collective resonance has exactly resistive source impedance.
            near(result.output[0], voltage.scale(1.0 / r), 2e-12);
            near(result.output[1], voltage.scale(symmetry / r), 2e-12);
            let heat = voltage.norm_sq() / r;
            assert!((result.dissipated_power_w / heat - 1.0).abs() < 3e-12);
            assert!(result.power_residual_w.abs() < heat * 2e-12);
        }
    }
}

#[test]
fn algebraic_transformer_is_not_mistaken_for_zero_storage_zero_output() {
    let (turns, resistance) = (3.0, 12.0);
    let graph = CircuitGraph {
        node_count: 3,
        branches: vec![
            (1, 0, Branch::VoltageSource { port: 0 }),
            (2, 0, Branch::Resistor { ohms: resistance }),
        ],
        transformers: vec![TransformerLink {
            primary: (1, 0),
            secondary: (2, 0),
            ratio: turns,
        }],
    };
    let circuit = assemble_circuit(&graph).unwrap();
    assert_eq!(circuit.system.differential_dim(), 0);
    let model = HarmonicModel::from_descriptor(&circuit.system, budget(), || false).unwrap();
    let voltage = C64::new(12.0, -4.0);
    let result = model
        .factor(600.0, || false)
        .unwrap()
        .response(&[voltage], || false)
        .unwrap();
    near(
        result.output[0],
        voltage.scale(1.0 / (turns * turns * resistance)),
        1e-14,
    );
    near(
        result.state[circuit.node_potential_index[1]],
        voltage.scale(1.0 / turns),
        1e-14,
    );
    assert_eq!(result.mean_stored_energy_j, 0.0);
    let power = 0.5 * voltage.norm_sq() / (turns * turns * resistance);
    assert!((result.supplied_power_w / power - 1.0).abs() < 1e-14);
    assert!(result.power_residual_w.abs() < power * 1e-14);
}

#[test]
fn mixed_coordinate_scales_do_not_change_the_resonance_response() {
    for scale in [1e-120, 1.0, 1e120] {
        let system = fs_phs::mass_spring_damper(scale, 9.0 * scale, 0.02 * scale).unwrap();
        let model = HarmonicModel::from_ode(&system, budget(), || false).unwrap();
        let result = model
            .factor(3.0, || false)
            .unwrap()
            .response(&[C64::from_re(0.06 * scale)], || false)
            .unwrap();
        near(result.state[0], C64::new(0.0, 1.0), 1e-12);
        near(result.output[0], C64::from_re(3.0), 1e-12);
        assert!((result.mean_stored_energy_j / (4.5 * scale) - 1.0).abs() < 1e-12);
        assert!(result.power_residual_w.abs() < 1e-12 * scale);
    }
}

#[test]
fn true_lossless_pole_refuses_and_cancellation_does_not_poison_reusable_factor() {
    let lossless = fs_phs::mass_spring_damper(1.0, 4.0, 0.0).unwrap();
    let model = HarmonicModel::from_ode(&lossless, budget(), || false).unwrap();
    assert!(matches!(
        model.factor(2.0, || false),
        Err(HarmonicError::Singular { .. })
    ));
    assert!(matches!(
        model.factor(1.0, || true),
        Err(HarmonicError::Cancelled)
    ));
    let factor = model.factor(1.0, || false).unwrap();
    let before = factor.response(&[C64::ONE], || false).unwrap();
    let mut calls = 0;
    assert!(matches!(
        factor.response(&[C64::ONE], || {
            calls += 1;
            calls == 5
        }),
        Err(HarmonicError::Cancelled)
    ));
    assert_eq!(before, factor.response(&[C64::ONE], || false).unwrap());
    assert!(matches!(
        HarmonicModel::from_ode(&lossless, budget(), || true),
        Err(HarmonicError::Cancelled)
    ));
}

#[test]
fn nonlinear_storage_malformed_inputs_and_insufficient_budgets_refuse() {
    let nonlinear = PortHamiltonian::new(
        1,
        1,
        vec![0.0],
        vec![1.0],
        vec![1.0],
        Box::new(SeparableStorage::new(vec![ScalarLaw {
            h: |x| x.powi(4) * 0.25,
            dh: |x| x.powi(3),
        }])),
    )
    .unwrap();
    assert!(HarmonicModel::from_ode(&nonlinear, budget(), || false).is_err());
    let system = fs_phs::mass_spring_damper(1.0, 4.0, 0.2).unwrap();
    for bad in [
        HarmonicBudget {
            max_dimension: 1,
            ..budget()
        },
        HarmonicBudget {
            max_dense_bytes: 1,
            ..budget()
        },
        HarmonicBudget {
            max_dense_work: 1,
            ..budget()
        },
        HarmonicBudget {
            relative_residual_tolerance: f64::NAN,
            ..budget()
        },
    ] {
        assert!(HarmonicModel::from_ode(&system, bad, || false).is_err());
    }
    let model = HarmonicModel::from_ode(&system, budget(), || false).unwrap();
    for omega in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(model.factor(omega, || false).is_err());
    }
    let factor = model.factor(2.0, || false).unwrap();
    assert!(factor.response(&[], || false).is_err());
    assert!(
        factor
            .response(&[C64::new(1.0, f64::NAN)], || false)
            .is_err()
    );
    // Raw structural misuse must produce a typed refusal, never index a matrix.
    let malformed = PortHamiltonian::from_raw_parts(
        2,
        1,
        vec![],
        vec![],
        vec![],
        Box::new(QuadraticStorage::new(vec![1.0], 1).unwrap()),
    );
    assert!(HarmonicModel::from_ode(&malformed, budget(), || false).is_err());
}
