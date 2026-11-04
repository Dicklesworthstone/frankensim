//! G1 continuous-response tangents and G3 original-coordinate/constraint checks.
use fs_math::c64::C64;
use fs_phs::circuit::{Branch, CircuitGraph, TransformerLink, assemble_circuit};
use fs_phs::harmonic::{HarmonicBudget, HarmonicError, HarmonicModel};
use fs_phs::{PortHamiltonian, QuadraticStorage, RelaxationBranch};

fn budget() -> HarmonicBudget {
    HarmonicBudget {
        max_dimension: 64,
        max_dense_bytes: 2_000_000,
        max_dense_work: 4_000_000,
        relative_residual_tolerance: 1e-12,
    }
}

fn near(actual: C64, expected: C64, tolerance: f64) {
    let scale = if expected == C64::ZERO {
        1.0
    } else {
        expected.abs()
    };
    assert!(
        (actual - expected).abs() <= tolerance * scale,
        "actual {actual:?}, expected {expected:?}"
    );
}

fn near_real(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "actual {actual:e}, expected {expected:e}"
    );
}

#[test]
fn oscillator_derivatives_resolve_high_q_amplitude_phase_and_physical_energy() {
    let (mass, stiffness) = (2.0, 8.0);
    for damping in [0.6, 1e-7] {
        let system = fs_phs::mass_spring_damper(mass, stiffness, damping).unwrap();
        let model = HarmonicModel::from_ode(&system, budget(), || false).unwrap();
        for omega in [1.0, 1.999, 2.0, 2.001, 4.0] {
            let factor = model.factor(omega, || false).unwrap();
            let result = factor
                .response_with_frequency_derivative(&[C64::ONE], None, || false)
                .unwrap();
            assert_eq!(
                result.response,
                factor.response(&[C64::ONE], || false).unwrap()
            );
            let dynamic = C64::new(stiffness - mass * omega * omega, -damping * omega);
            let dynamic_derivative = C64::new(-2.0 * mass * omega, -damping);
            let q = C64::ONE / dynamic;
            let dq = -dynamic_derivative / (dynamic * dynamic);
            let velocity = C64::new(0.0, -omega) * q;
            let dv = C64::new(0.0, -1.0) * (q + dq.scale(omega));
            near(result.state_derivative[0], dq, 1e-8);
            near(result.state_derivative[1], dv.scale(mass), 1e-8);
            near(result.output_derivative[0], dv, 1e-8);
            near(result.effort_derivative[0], dq.scale(stiffness), 1e-8);
            let dh = 0.5 * stiffness * (q.conj() * dq).re + 0.5 * mass * (velocity.conj() * dv).re;
            let dloss = damping * (velocity.conj() * dv).re;
            near_real(result.mean_stored_energy_derivative, dh, 2e-8);
            near_real(result.dissipated_power_derivative, dloss, 2e-8);
            near_real(result.supplied_power_derivative, 0.5 * dv.re, 2e-8);
            assert!(result.relative_derivative_residual < 1e-12);
            if omega == 2.0 {
                let slope = result.output_logarithmic_derivative(0).unwrap().unwrap();
                near(slope, C64::new(0.0, 2.0 * mass / damping), 1e-8);
            }
        }
    }
}

fn series_rlc(r: f64, l: f64, c: f64) -> CircuitGraph {
    CircuitGraph {
        node_count: 4,
        branches: vec![
            (1, 0, Branch::VoltageSource { port: 0 }),
            (1, 2, Branch::Inductor { henries: l }),
            (2, 3, Branch::Capacitor { farads: c }),
            (3, 0, Branch::Resistor { ohms: r }),
        ],
        transformers: vec![],
    }
}

#[test]
fn rlc_tangent_keeps_kirchhoff_reactions_and_frequency_dependent_source_work() {
    let (r, l, c) = (0.02, 0.5, 0.125);
    let circuit = assemble_circuit(&series_rlc(r, l, c)).unwrap();
    let model = HarmonicModel::from_descriptor(&circuit.system, budget(), || false).unwrap();
    let (voltage, du) = (C64::new(3.0, -0.7), C64::new(0.2, 0.5));
    for omega in [1.0, 3.99, 4.0, 4.01, 7.0] {
        let result = model
            .factor(omega, || false)
            .unwrap()
            .response_with_frequency_derivative(&[voltage], Some(&[du]), || false)
            .unwrap();
        let impedance = C64::new(r, 1.0 / (omega * c) - omega * l);
        let dz = C64::new(0.0, -1.0 / (omega * omega * c) - l);
        let current = voltage / impedance;
        let di = (du - dz * current) / impedance;
        let charge = C64::new(0.0, 1.0 / omega) * current;
        let dq = C64::new(0.0, 1.0 / omega) * (di - current.scale(1.0 / omega));
        near(result.output_derivative[0], di, 2e-10);
        near(
            result.state_derivative[circuit.flux_index[0]],
            di.scale(l),
            2e-10,
        );
        near(result.state_derivative[circuit.charge_index[0]], dq, 2e-10);
        near(
            result.state_derivative[circuit.node_potential_index[0]],
            du,
            2e-10,
        );
        near(
            result.state_derivative[circuit.node_potential_index[2]],
            di.scale(r),
            2e-10,
        );
        let dh = 0.5 * l * (current.conj() * di).re + 0.5 / c * (charge.conj() * dq).re;
        let dp = 0.5 * (du.conj() * current + voltage.conj() * di).re;
        near_real(result.mean_stored_energy_derivative, dh, 2e-9);
        near_real(result.supplied_power_derivative, dp, 2e-9);
        near_real(
            result.dissipated_power_derivative,
            r * (current.conj() * di).re,
            2e-9,
        );
        assert!(result.power_residual_derivative.abs() < 2e-9 * dp.abs().max(1.0));
    }
}

#[test]
fn maxwell_memory_derivative_retains_dispersion_internal_motion_and_heat() {
    let (mass, equilibrium, branch, tau) = (0.4, 110.0, 300.0_f64, 0.025);
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
            .response_with_frequency_derivative(&[C64::ONE], None, || false)
            .unwrap();
        let st = C64::new(0.0, -omega * tau);
        let dst = C64::new(0.0, -tau);
        let memory = C64::ONE + st;
        let dynamic =
            C64::from_re(equilibrium - mass * omega * omega) + (st / memory).scale(branch);
        let ddynamic = C64::from_re(-2.0 * mass * omega) + (dst / (memory * memory)).scale(branch);
        let q = C64::ONE / dynamic;
        let dq = -ddynamic / (dynamic * dynamic);
        let dz = (dq / memory - q * dst / (memory * memory)).scale(branch.sqrt());
        near(result.state_derivative[0], dq, 2e-11);
        near(result.state_derivative[2], dz, 2e-11);
        let denominator = 1.0 + (omega * tau).powi(2);
        let dloss = branch * tau * omega / (denominator * denominator) * q.norm_sq()
            + branch * tau * omega * omega / denominator * (q.conj() * dq).re;
        near_real(result.dissipated_power_derivative, dloss, 2e-12);
        assert!(result.power_residual_derivative.abs() < 2e-12);
    }
}

#[test]
fn purely_algebraic_circuit_has_no_spurious_frequency_motion() {
    let (turns, resistance) = (3.0, 12.0);
    let circuit = assemble_circuit(&CircuitGraph {
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
    })
    .unwrap();
    assert_eq!(circuit.system.differential_dim(), 0);
    let model = HarmonicModel::from_descriptor(&circuit.system, budget(), || false).unwrap();
    let factor = model.factor(600.0, || false).unwrap();
    let (u, du) = (C64::new(12.0, -4.0), C64::new(-0.3, 0.9));
    let fixed = factor
        .response_with_frequency_derivative(&[u], None, || false)
        .unwrap();
    assert!(fixed.state_derivative.iter().all(|&x| x == C64::ZERO));
    assert_eq!(
        fixed.output_logarithmic_derivative(0).unwrap(),
        Some(C64::ZERO)
    );
    let driven = factor
        .response_with_frequency_derivative(&[u], Some(&[du]), || false)
        .unwrap();
    near(
        driven.output_derivative[0],
        du.scale(1.0 / (turns * turns * resistance)),
        1e-14,
    );
    near(
        driven.state_derivative[circuit.node_potential_index[1]],
        du.scale(1.0 / turns),
        1e-14,
    );
    near_real(
        driven.supplied_power_derivative,
        (u.conj() * du).re / (turns * turns * resistance),
        1e-14,
    );
    assert_eq!(driven.mean_stored_energy_derivative, 0.0);
}

#[test]
fn exact_antiresonance_retains_finite_tangent_without_a_phase_floor() {
    let mut q = vec![0.0; 16];
    for (i, value) in [2.0, 1.0, 6.0, 1.0].into_iter().enumerate() {
        q[5 * i] = value;
    }
    let mut j = vec![0.0; 16];
    j[1] = 1.0;
    j[4] = -1.0;
    j[11] = 1.0;
    j[14] = -1.0;
    let system = PortHamiltonian::new(
        4,
        1,
        j,
        vec![0.0; 16],
        vec![0.0, 1.0, 0.0, 1.0],
        Box::new(QuadraticStorage::new(q, 4).unwrap()),
    )
    .unwrap();
    let model = HarmonicModel::from_ode(&system, budget(), || false).unwrap();
    let factor = model.factor(2.0, || false).unwrap();
    let result = factor
        .response_with_frequency_derivative(&[C64::ONE], None, || false)
        .unwrap();
    near(result.response.output[0], C64::ZERO, 2e-15);
    near(result.output_derivative[0], C64::new(0.0, -4.0), 2e-14);
    if result.response.output[0] == C64::ZERO {
        assert_eq!(result.output_logarithmic_derivative(0).unwrap(), None);
    }
    // Exactly zero excitation, including a nonzero source tangent, has no
    // output phase. The tangent itself remains a defined physical response.
    let zero = factor
        .response_with_frequency_derivative(&[C64::ZERO], Some(&[C64::ONE]), || false)
        .unwrap();
    assert_eq!(zero.output_logarithmic_derivative(0).unwrap(), None);
    assert!(zero.output_logarithmic_derivative(1).is_err());
}

#[test]
fn dimensional_scaling_preserves_the_same_resonant_slope() {
    for scale in [1e-120, 1.0, 1e120] {
        let system = fs_phs::mass_spring_damper(scale, 9.0 * scale, 0.02 * scale).unwrap();
        let model = HarmonicModel::from_ode(&system, budget(), || false).unwrap();
        let result = model
            .factor(3.0, || false)
            .unwrap()
            .response_with_frequency_derivative(&[C64::from_re(0.06 * scale)], None, || false)
            .unwrap();
        near(result.output_derivative[0], C64::new(0.0, 300.0), 2e-11);
        near(
            result.output_logarithmic_derivative(0).unwrap().unwrap(),
            C64::new(0.0, 100.0),
            2e-11,
        );
    }
}

#[test]
fn derivative_admission_and_cancellation_leave_the_frequency_factor_reusable() {
    let system = fs_phs::mass_spring_damper(1.0, 4.0, 0.2).unwrap();
    let model = HarmonicModel::from_ode(&system, budget(), || false).unwrap();
    let factor = model.factor(2.0, || false).unwrap();
    let before = factor.response(&[C64::ONE], || false).unwrap();
    let mut calls = 0;
    let expected = factor
        .response_with_frequency_derivative(&[C64::ONE], None, || {
            calls += 1;
            false
        })
        .unwrap();
    for stop in [1, calls / 2, calls - 1] {
        let mut count = 0;
        assert!(matches!(
            factor.response_with_frequency_derivative(&[C64::ONE], None, || {
                count += 1;
                count == stop
            }),
            Err(HarmonicError::Cancelled)
        ));
    }
    for invalid in [vec![], vec![C64::from_re(f64::NAN)]] {
        assert!(
            factor
                .response_with_frequency_derivative(&[C64::ONE], Some(&invalid), || false)
                .is_err()
        );
    }
    assert_eq!(factor.response(&[C64::ONE], || false).unwrap(), before);
    assert_eq!(
        factor
            .response_with_frequency_derivative(&[C64::ONE], None, || false)
            .unwrap(),
        expected
    );
    for tight in [
        HarmonicBudget {
            max_dense_work: 100,
            ..budget()
        },
        HarmonicBudget {
            max_dense_bytes: 1000,
            ..budget()
        },
    ] {
        let limited = HarmonicModel::from_ode(&system, tight, || false).unwrap();
        let factor = limited.factor(2.0, || false).unwrap();
        factor.response(&[C64::ONE], || false).unwrap();
        assert!(matches!(
            factor.response_with_frequency_derivative(&[C64::ONE], None, || false),
            Err(HarmonicError::Invalid { .. })
        ));
    }
}
