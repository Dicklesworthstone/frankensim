//! G1 physical spectra, G3 hidden-mode completeness evidence, G4 cancellation.
use fs_math::c64::C64;
use fs_phs::circuit::{Branch, CircuitGraph, assemble_circuit};
use fs_phs::harmonic::{BeynOptions, HarmonicBudget, HarmonicModel, PoleOptions, ResonanceRegion};
use fs_phs::{PortHamiltonian, QuadraticStorage, RelaxationBranch};
use std::cell::Cell;

fn budget() -> HarmonicBudget {
    HarmonicBudget {
        max_dimension: 64,
        max_dense_bytes: 2_000_000,
        max_dense_work: 10_000_000,
        relative_residual_tolerance: 1e-12,
    }
}
fn options(n: usize, scale: f64) -> BeynOptions {
    BeynOptions {
        nodes: 128,
        probe_columns: n,
        rank_relative_tolerance: 1e-10,
        count_tolerance: 1e-6,
        quadrature_tolerance: 1e-6,
        contour_min_reciprocal_condition: 1e-14,
        linear_residual_tolerance: 1e-12,
        svd_tolerance: 1e-13,
        svd_max_sweeps: 100,
        eigen_max_sweeps: 120,
        max_dimension: 64,
        max_dense_bytes: 2_000_000,
        max_work: 2_000_000_000,
        refinement: PoleOptions {
            spectral_scale: scale,
            step_tolerance: 1e-13,
            residual_tolerance: 1e-13,
            max_iterations: 32,
            max_dimension: 64,
            max_dense_bytes: 2_000_000,
            max_work: 100_000_000,
        },
    }
}

#[test]
fn complete_mechanical_spectrum_retains_an_unobservable_resonance() {
    // Two physical oscillators; the only applied/observed port acts on the
    // first. Discovery must retain the second oscillator without a response seed.
    let j = vec![
        0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0,
    ];
    let mut r = vec![0.0; 16];
    r[10] = 0.2;
    r[15] = 0.4;
    let mut q = vec![0.0; 16];
    q[0] = 100.0;
    q[5] = 225.0;
    q[10] = 1.0;
    q[15] = 1.0;
    let system = PortHamiltonian::new(
        4,
        1,
        j,
        r,
        vec![0.0, 0.0, 1.0, 0.0],
        Box::new(QuadraticStorage::new(q, 4).unwrap()),
    )
    .unwrap();
    let model = HarmonicModel::from_ode(&system, budget(), || false).unwrap();
    let region = ResonanceRegion {
        center: C64::new(-0.15, -12.5),
        real_radius: 3.0,
        imaginary_radius: 5.0,
    };
    let result = model
        .resonance_spectrum(region, &options(4, 15.0), || false)
        .unwrap();
    assert_eq!(
        (
            result.argument_count,
            result.resolved_rank,
            result.poles.len()
        ),
        (2, 2, 2)
    );
    for (gamma, stiffness, dark) in [(0.1_f64, 100.0_f64, false), (0.2, 225.0, true)] {
        let target = C64::new(-gamma, -(stiffness - gamma * gamma).sqrt());
        let pole = result
            .poles
            .iter()
            .find(|p| (p.laplace_pole_per_s - target).abs() < 1e-10)
            .unwrap();
        assert_eq!(pole.state.len(), 4);
        assert!(pole.stored_energy_j > 0.0 && pole.dissipated_power_w > 0.0);
        assert!(pole.decay_power_residual_w.abs() < 1e-10 * pole.dissipated_power_w);
        if dark {
            assert!(pole.port_output[0].abs() < 1e-13);
        } else {
            assert!(pole.port_output[0].abs() > 0.1);
        }
    }
}

#[test]
fn complete_hereditary_spectrum_contains_relaxation_and_both_oscillatory_modes() {
    let (mass, k0, kb, tau) = (0.4_f64, 110.0_f64, 300.0_f64, 0.025_f64);
    let system = fs_phs::mass_spring_damper(mass, k0, 0.0)
        .unwrap()
        .with_relaxation_branches(vec![RelaxationBranch {
            projection: vec![1.0, 0.0],
            stiffness: kb,
            relaxation_time_s: tau,
        }])
        .unwrap();
    let model = HarmonicModel::from_ode(&system, budget(), || false).unwrap();
    let region = ResonanceRegion {
        center: C64::from_re(-15.0),
        real_radius: 30.0,
        imaginary_radius: 40.0,
    };
    let result = model
        .resonance_spectrum(region, &options(3, 40.0), || false)
        .unwrap();
    assert_eq!(result.argument_count, 3);
    assert_eq!(
        result
            .poles
            .iter()
            .filter(|p| p.laplace_pole_per_s.im.abs() < 1e-9)
            .count(),
        1
    );
    for pole in &result.poles {
        let s = pole.laplace_pole_per_s;
        assert!(s.re < 0.0);
        let original =
            (s * s).scale(mass) + C64::from_re(k0) + s.scale(tau * kb) / (C64::ONE + s.scale(tau));
        assert!(original.abs() < 1e-9, "{pole:?}");
        assert!(
            (pole.state[2] - pole.state[0].scale(kb.sqrt()) / (C64::ONE + s.scale(tau))).abs()
                < 1e-11
        );
        assert!(pole.decay_power_residual_w.abs() < 1e-10 * pole.dissipated_power_w);
    }
}

#[test]
fn descriptor_spectrum_counts_only_finite_rlc_modes_and_retains_kirchhoff_rows() {
    let (r, l, c) = (0.2_f64, 0.1_f64, 2.2e-7_f64);
    let circuit = assemble_circuit(&CircuitGraph {
        node_count: 4,
        branches: vec![
            (1, 0, Branch::VoltageSource { port: 0 }),
            (1, 2, Branch::Inductor { henries: l }),
            (2, 3, Branch::Capacitor { farads: c }),
            (3, 0, Branch::Resistor { ohms: r }),
        ],
        transformers: vec![],
    })
    .unwrap();
    let n = circuit.system.state_dim();
    let model = HarmonicModel::from_descriptor(&circuit.system, budget(), || false).unwrap();
    let omega0 = 1.0 / (l * c).sqrt();
    let region = ResonanceRegion {
        center: C64::ZERO,
        real_radius: omega0 * 1.5,
        imaginary_radius: omega0 * 1.5,
    };
    let result = model
        .resonance_spectrum(region, &options(n, omega0), || false)
        .unwrap();
    assert_eq!((result.argument_count, result.poles.len()), (2, 2));
    for pole in &result.poles {
        let s = pole.laplace_pole_per_s;
        assert!((s * s + C64::from_re(1.0 / (l * c)) + s.scale(r / l)).abs() < 1e-6);
        assert_eq!(pole.state.len(), n);
        assert!(pole.state[circuit.node_potential_index[0]].abs() < 1e-12);
        let current = pole.state[circuit.flux_index[0]].scale(1.0 / l);
        assert!((pole.state[circuit.node_potential_index[2]] - current.scale(r)).abs() < 1e-10);
        assert!(pole.relative_residual < 1e-12);
    }
}

#[test]
fn cancelled_spectrum_preserves_reusable_ac_factor_and_replays() {
    let system = fs_phs::mass_spring_damper(2.5, 1600.0, 0.6).unwrap();
    let model = HarmonicModel::from_ode(&system, budget(), || false).unwrap();
    let factor = model.factor(25.0, || false).unwrap();
    let before = factor.response(&[C64::ONE], || false).unwrap();
    let region = ResonanceRegion {
        center: C64::ZERO,
        real_radius: 40.0,
        imaginary_radius: 40.0,
    };
    let polls = Cell::new(0);
    let interrupted = model.resonance_spectrum(region, &options(2, 40.0), || {
        polls.set(polls.get() + 1);
        polls.get() > 100
    });
    assert!(interrupted.is_err());
    let first = model
        .resonance_spectrum(region, &options(2, 40.0), || false)
        .unwrap();
    let second = model
        .resonance_spectrum(region, &options(2, 40.0), || false)
        .unwrap();
    assert_eq!(first.argument_count, 2);
    assert_eq!(first.poles, second.poles);
    let after = factor.response(&[C64::ONE], || false).unwrap();
    assert_eq!(before.state, after.state);
    assert_eq!(before.output, after.output);
    assert_eq!(
        before.dissipated_power_w.to_bits(),
        after.dissipated_power_w.to_bits()
    );
}
