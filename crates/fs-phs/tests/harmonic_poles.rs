//! Physical poles of the original continuous pHS and descriptor models.
use fs_la::nonlinear_eig::NonlinearEigError;
use fs_math::c64::C64;
use fs_phs::RelaxationBranch;
use fs_phs::circuit::{Branch, CircuitGraph, assemble_circuit};
use fs_phs::harmonic::{HarmonicBudget, HarmonicModel, HarmonicPoleError, PoleOptions};
use std::cell::Cell;

fn model_budget() -> HarmonicBudget {
    HarmonicBudget {
        max_dimension: 32,
        max_dense_bytes: 512_000,
        max_dense_work: 1_000_000,
        relative_residual_tolerance: 1e-12,
    }
}

fn pole_options(spectral_scale: f64) -> PoleOptions {
    PoleOptions {
        spectral_scale,
        step_tolerance: 1e-13,
        residual_tolerance: 1e-13,
        max_iterations: 40,
        max_dimension: 32,
        max_dense_bytes: 512_000,
        max_work: 1_000_000,
    }
}

#[test]
fn damped_mechanical_pole_retains_original_coordinates_output_and_decay() {
    let (mass, stiffness, damping) = (2.5_f64, 1600.0_f64, 0.6_f64);
    let system = fs_phs::mass_spring_damper(mass, stiffness, damping).unwrap();
    let model = HarmonicModel::from_ode(&system, model_budget(), || false).unwrap();
    let gamma = damping / (2.0 * mass);
    let expected = C64::new(-gamma, -(stiffness / mass - gamma * gamma).sqrt());
    let pole = model
        .resonance_pole(
            expected + C64::new(0.03, 0.5),
            &[C64::ONE, C64::new(-0.25, -60.0)],
            &pole_options((stiffness / mass).sqrt()),
            || false,
        )
        .unwrap();

    assert!(
        (pole.laplace_pole_per_s - expected).abs() < 1e-10,
        "{pole:?}"
    );
    assert!((pole.state[1] - pole.laplace_pole_per_s.scale(mass) * pole.state[0]).abs() < 1e-11);
    assert!((pole.port_output[0] - pole.state[1].scale(1.0 / mass)).abs() < 1e-12);
    assert_eq!(
        pole.angular_frequency_rad_s,
        C64::new(-pole.laplace_pole_per_s.im, pole.laplace_pole_per_s.re)
    );
    assert!(pole.stored_energy_j > 0.0 && pole.dissipated_power_w > 0.0);
    assert!(
        pole.decay_power_residual_w.abs() < 1e-11 * pole.dissipated_power_w,
        "{pole:?}"
    );
    assert!(pole.relative_residual < 1e-12);
}

#[test]
fn singular_descriptor_rlc_pole_retains_kirchhoff_constraints_and_decay() {
    let (resistance, inductance, capacitance) = (0.2_f64, 0.1_f64, 2.2e-7_f64);
    let circuit = assemble_circuit(&CircuitGraph {
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
    })
    .unwrap();
    let model = HarmonicModel::from_descriptor(&circuit.system, model_budget(), || false).unwrap();
    let omega0 = 1.0 / (inductance * capacitance).sqrt();
    let seed = model
        .factor(omega0, || false)
        .unwrap()
        .response(&[C64::ONE], || false)
        .unwrap()
        .state;
    let gamma = resistance / (2.0 * inductance);
    let expected = C64::new(-gamma, -(omega0 * omega0 - gamma * gamma).sqrt());
    let pole = model
        .resonance_pole(
            expected + C64::new(0.1, 2.0),
            &seed,
            &pole_options(omega0),
            || false,
        )
        .unwrap();

    assert!(
        (pole.laplace_pole_per_s - expected).abs() < 1e-7,
        "{pole:?}"
    );
    // A free mode keeps the ideal voltage source's homogeneous constraint.
    assert!(pole.state[circuit.node_potential_index[0]].abs() < 1e-12);
    let current = pole.state[circuit.flux_index[0]].scale(1.0 / inductance);
    assert!(
        (pole.state[circuit.node_potential_index[2]] - current.scale(resistance)).abs() < 1e-10
    );
    assert!(pole.relative_residual < 1e-12);
    assert!(pole.decay_power_residual_w.abs() < 1e-9 * pole.dissipated_power_w);
}

#[test]
fn hereditary_pole_retains_frequency_dependent_storage_loss_and_memory() {
    let (mass, equilibrium, branch, tau) = (0.4_f64, 110.0_f64, 300.0_f64, 0.025_f64);
    let system = fs_phs::mass_spring_damper(mass, equilibrium, 0.0)
        .unwrap()
        .with_relaxation_branches(vec![RelaxationBranch {
            projection: vec![1.0, 0.0],
            stiffness: branch,
            relaxation_time_s: tau,
        }])
        .unwrap();
    let model = HarmonicModel::from_ode(&system, model_budget(), || false).unwrap();
    let seed = model
        .factor(22.0, || false)
        .unwrap()
        .response(&[C64::ONE], || false)
        .unwrap()
        .state;
    let pole = model
        .resonance_pole(C64::new(-11.0, -22.0), &seed, &pole_options(32.0), || false)
        .unwrap();

    assert_eq!(pole.state.len(), 3);
    let s = pole.laplace_pole_per_s;
    // Eliminate the Maxwell coordinate analytically, independently of the
    // prepared descriptor generator: M*s^2 + K0 + Kb*s*tau/(1+s*tau)=0.
    let dispersion = (s * s).scale(mass)
        + C64::from_re(equilibrium)
        + s.scale(tau * branch) / (C64::ONE + s.scale(tau));
    assert!(
        dispersion.abs() < 1e-9,
        "dispersion={dispersion:?}, pole={pole:?}"
    );
    assert!(
        (pole.state[2] - pole.state[0].scale(branch.sqrt()) / (C64::ONE + s.scale(tau))).abs()
            < 1e-11
    );
    let exact_loss = (pole.state[2] - pole.state[0].scale(branch.sqrt())).norm_sq() / tau;
    assert!((pole.dissipated_power_w / exact_loss - 1.0).abs() < 1e-11);
    assert!(pole.decay_power_residual_w.abs() < 1e-10 * exact_loss);
    assert!(s.re < 0.0 && s.im < 0.0);
}

#[test]
fn pole_cancellation_and_budget_refusal_preserve_replay_and_reusable_ac_factor() {
    let system = fs_phs::mass_spring_damper(2.5, 1600.0, 0.6).unwrap();
    let model = HarmonicModel::from_ode(&system, model_budget(), || false).unwrap();
    let omega = (1600.0_f64 / 2.5).sqrt();
    let factor = model.factor(omega, || false).unwrap();
    let before = factor.response(&[C64::new(2.0, -0.3)], || false).unwrap();
    let start = C64::new(-0.09, -omega + 0.5);
    let seed = [C64::ONE, C64::new(-0.25, -60.0)];
    let options = pole_options(omega);

    let polls = Cell::new(0_usize);
    let reference = model
        .resonance_pole(start, &seed, &options, || {
            polls.set(polls.get() + 1);
            false
        })
        .unwrap();
    let cut = polls.get() / 2;
    assert!(cut > 1, "exercise cancellation during refinement");
    let interrupted = Cell::new(0_usize);
    let cancelled = model.resonance_pole(start, &seed, &options, || {
        interrupted.set(interrupted.get() + 1);
        interrupted.get() >= cut
    });
    assert!(matches!(
        cancelled,
        Err(HarmonicPoleError::Solve(NonlinearEigError::Cancelled))
    ));

    let mut insufficient = options;
    insufficient.max_dense_bytes = 0;
    assert!(matches!(
        model.resonance_pole(start, &seed, &insufficient, || false),
        Err(HarmonicPoleError::Solve(NonlinearEigError::BadInput { .. }))
    ));

    let replay = model
        .resonance_pole(start, &seed, &options, || false)
        .unwrap();
    assert_eq!(reference, replay);
    let after = factor.response(&[C64::new(2.0, -0.3)], || false).unwrap();
    assert_eq!(before.state, after.state);
    assert_eq!(before.output, after.output);
    assert_eq!(
        before.dissipated_power_w.to_bits(),
        after.dissipated_power_w.to_bits()
    );
}
