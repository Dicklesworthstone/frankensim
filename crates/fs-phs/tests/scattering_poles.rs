//! Independent original-equation and energy oracles for outgoing natural modes.
use fs_la::nonlinear_eig::NonlinearEigError;
use fs_math::c64::C64;
use fs_phs::harmonic::scattering::{ScatteringBudget, ScatteringModel, ScatteringPole};
use fs_phs::harmonic::{
    BeynOptions, HarmonicBudget, HarmonicModel, HarmonicPoleError, PoleOptions, ResonanceRegion,
};
use fs_phs::{PortHamiltonian, QuadraticStorage, RelaxationBranch};

fn budget() -> HarmonicBudget {
    HarmonicBudget {
        max_dimension: 64,
        max_dense_bytes: 4_000_000,
        max_dense_work: 8_000_000,
        relative_residual_tolerance: 2e-11,
    }
}
fn wave_budget() -> ScatteringBudget {
    ScatteringBudget {
        maximum_ports: 8,
        maximum_payload_bytes: 1_000_000,
        maximum_work: 2_000_000,
    }
}
fn prepare(system: &PortHamiltonian, ports: &[(usize, f64)]) -> ScatteringModel {
    ScatteringModel::new(
        HarmonicModel::from_ode(system, budget(), || false).unwrap(),
        ports,
        wave_budget(),
        || false,
    )
    .unwrap()
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
            max_iterations: 40,
            max_dimension: 64,
            max_dense_bytes: 2_000_000,
            max_work: 100_000_000,
        },
    }
}
fn close(a: C64, b: C64, floor: f64) {
    assert!(
        (a - b).abs() < 3e-10 * a.abs().max(b.abs()).max(floor),
        "{a:?} != {b:?}"
    );
}
fn balance(p: &ScatteringPole) {
    let scale = p.internal_dissipated_power_w.abs() + p.outgoing_power_w + 1e-24;
    assert!(p.relative_residual < 2e-11, "{p:?}");
    assert!(p.decay_power_residual_w.abs() < 3e-10 * scale, "{p:?}");
    assert!(p.physical_power_residual_w.abs() < 3e-10 * scale, "{p:?}");
    for w in &p.waves {
        assert!(w.power_residual_w.abs() < 3e-10 * scale, "{p:?}");
        close(
            w.physical_input,
            w.outgoing.scale(-1.0 / w.reference),
            1e-15,
        );
    }
}

#[test]
fn matched_boundary_causes_decay_without_manufacturing_material_heat() {
    let (mass, stiffness, reference) = (2.5_f64, 1600.0_f64, 0.5_f64);
    for damping in [0.0, 0.6] {
        let system = fs_phs::mass_spring_damper(mass, stiffness, damping).unwrap();
        let model = prepare(&system, &[(0, reference)]);
        let gamma = (damping + 1.0 / reference) / (2.0 * mass);
        let target = C64::new(-gamma, -(stiffness / mass - gamma * gamma).sqrt());
        let pole = model
            .resonance_pole(
                target + C64::new(0.05, 0.3),
                &[C64::ONE, target.scale(mass)],
                &options(2, 30.0).refinement,
                || false,
            )
            .unwrap();
        close(pole.laplace_pole_per_s, target, 1e-12);
        let v = pole.state[1].scale(1.0 / mass);
        close(v, target * pole.state[0], 1e-12);
        let heat = damping * v.norm_sq();
        let outgoing = v.norm_sq() / reference;
        assert!((pole.internal_dissipated_power_w - heat).abs() < 1e-12);
        assert!((pole.outgoing_power_w - outgoing).abs() < 1e-12);
        if damping == 0.0 {
            assert_eq!(pole.internal_dissipated_power_w, 0.0);
            assert!(pole.outgoing_power_w > 0.0);
        }
        assert!(
            (pole.stored_energy_j
                - 0.5 * (stiffness * pole.state[0].norm_sq() + mass * v.norm_sq()))
            .abs()
                < 1e-12
        );
        balance(&pole);
    }
}

#[test]
fn hereditary_spectrum_retains_relaxation_and_distinct_heat_and_outgoing_channels() {
    let (mass, k0, kb, tau, reference) = (0.4_f64, 110.0_f64, 300.0_f64, 0.025_f64, 0.2_f64);
    let system = fs_phs::mass_spring_damper(mass, k0, 0.0)
        .unwrap()
        .with_relaxation_branches(vec![RelaxationBranch {
            projection: vec![1.0, 0.0],
            stiffness: kb,
            relaxation_time_s: tau,
        }])
        .unwrap();
    let model = prepare(&system, &[(0, reference)]);
    let region = ResonanceRegion {
        center: C64::from_re(-25.0),
        real_radius: 40.0,
        imaginary_radius: 60.0,
    };
    let mut previous: Vec<C64> = Vec::new();
    for nodes in [64, 128] {
        let mut opts = options(3, 60.0);
        opts.nodes = nodes;
        let spectrum = model.resonance_spectrum(region, &opts, || false).unwrap();
        assert_eq!(
            (
                spectrum.argument_count,
                spectrum.resolved_rank,
                spectrum.poles.len()
            ),
            (3, 3, 3)
        );
        assert_eq!(
            spectrum
                .poles
                .iter()
                .filter(|p| p.laplace_pole_per_s.im.abs() < 1e-9)
                .count(),
            1
        );
        for p in &spectrum.poles {
            let s = p.laplace_pole_per_s;
            let polynomial = (s * s).scale(mass)
                + s.scale(1.0 / reference)
                + C64::from_re(k0)
                + s.scale(tau * kb) / (C64::ONE + s.scale(tau));
            assert!(polynomial.abs() < 1e-8, "{p:?}");
            close(
                p.state[2],
                p.state[0].scale(kb.sqrt()) / (C64::ONE + s.scale(tau)),
                1e-12,
            );
            let memory_effort = p.state[2] - p.state[0].scale(kb.sqrt());
            assert!((p.internal_dissipated_power_w - memory_effort.norm_sq() / tau).abs() < 1e-10);
            assert!((p.outgoing_power_w - (s * p.state[0]).norm_sq() / reference).abs() < 1e-10);
            assert_eq!(p.state.len(), 3);
            balance(p);
            if !previous.is_empty() {
                assert!(previous.iter().any(|&old| (old - s).abs() < 1e-9));
            }
        }
        previous = spectrum
            .poles
            .iter()
            .map(|p| p.laplace_pole_per_s)
            .collect();
    }
}

#[test]
fn mutual_wave_feedback_keeps_dark_modes_and_homogeneous_ordinary_ports() {
    // The same physical wave port acts on both momenta. Its off-diagonal
    // feedback damps the symmetric mode while the antisymmetric mode is dark.
    let j = vec![
        0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0,
    ];
    let mut q = vec![0.0; 16];
    q[0] = 100.0;
    q[5] = 100.0;
    q[10] = 1.0;
    q[15] = 1.0;
    let mut r = vec![0.0; 16];
    r[10] = 0.2;
    r[15] = 0.2;
    let system = PortHamiltonian::new(
        4,
        2,
        j,
        r,
        vec![0.0, 0.0, 0.0, 0.0, 1.0, 0.7, 1.0, -0.2],
        Box::new(QuadraticStorage::new(q, 4).unwrap()),
    )
    .unwrap();
    let model = prepare(&system, &[(0, 1.0)]);
    let spectrum = model
        .resonance_spectrum(
            ResonanceRegion {
                center: C64::new(-0.6, -10.0),
                real_radius: 3.0,
                imaginary_radius: 3.0,
            },
            &options(4, 10.0),
            || false,
        )
        .unwrap();
    assert_eq!((spectrum.argument_count, spectrum.poles.len()), (2, 2));
    for (gamma, dark) in [(0.1_f64, true), (1.1_f64, false)] {
        let expected = C64::new(-gamma, -(100.0 - gamma * gamma).sqrt());
        let p = spectrum
            .poles
            .iter()
            .find(|p| (p.laplace_pole_per_s - expected).abs() < 1e-9)
            .unwrap();
        assert_eq!(p.physical_inputs[1], C64::ZERO);
        assert_eq!(p.state.len(), 4);
        if dark {
            assert!(p.port_output[0].abs() < 1e-12);
            assert!(p.port_output[1].abs() > 0.1);
            assert!(p.outgoing_power_w < 1e-24);
        } else {
            assert!(p.outgoing_power_w > 0.1);
        }
        balance(p);
    }
}

#[test]
fn descriptor_rlc_outgoing_modes_preserve_actual_boundary_and_kirchhoff_reactions() {
    use fs_phs::circuit::{Branch, CircuitGraph, assemble_circuit};
    let (r, l, c, z) = (0.2_f64, 0.1_f64, 2.2e-7_f64, 2.0_f64);
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
    let model = ScatteringModel::new(
        HarmonicModel::from_descriptor(&circuit.system, budget(), || false).unwrap(),
        &[(0, z)],
        wave_budget(),
        || false,
    )
    .unwrap();
    let omega0 = 1.0 / (l * c).sqrt();
    let seed = model
        .factor(omega0, || false)
        .unwrap()
        .response(&[C64::ONE], || false)
        .unwrap()
        .physical
        .state;
    let gamma = (r + 1.0 / z) / (2.0 * l);
    let expected = C64::new(-gamma, -(omega0 * omega0 - gamma * gamma).sqrt());
    let p = model
        .resonance_pole(
            expected + C64::new(0.1, 2.0),
            &seed,
            &options(model.dimension(), omega0).refinement,
            || false,
        )
        .unwrap();
    close(p.laplace_pole_per_s, expected, 1e-12);
    let current = p.state[circuit.flux_index[0]].scale(1.0 / l);
    close(p.port_output[0], current, 1e-12);
    close(
        p.state[circuit.node_potential_index[0]],
        current.scale(-1.0 / z),
        1e-12,
    );
    close(
        p.state[circuit.node_potential_index[2]],
        current.scale(r),
        1e-12,
    );
    assert!((p.internal_dissipated_power_w - r * current.norm_sq()).abs() < 1e-11);
    assert!((p.outgoing_power_w - current.norm_sq() / z).abs() < 1e-11);
    assert_eq!(p.state.len(), circuit.system.state_dim());
    balance(&p);
}

#[test]
fn spectral_refusal_and_cancellation_do_not_change_reusable_ac_factor() {
    let system = fs_phs::mass_spring_damper(1.0, 4.0, 0.2).unwrap();
    let model = prepare(&system, &[(0, 2.0)]);
    let factor = model.factor(2.0, || false).unwrap();
    let before = factor.response(&[C64::ONE], || false).unwrap();
    let seed = &before.physical.state;
    let opts = options(2, 2.0);
    for limited in [
        PoleOptions {
            max_dense_bytes: 0,
            ..opts.refinement
        },
        PoleOptions {
            max_work: 0,
            ..opts.refinement
        },
    ] {
        assert!(
            model
                .resonance_pole(C64::new(-0.3, -2.0), seed, &limited, || false)
                .is_err()
        );
    }
    let mut polls = 0;
    assert!(matches!(
        model.resonance_pole(C64::new(-0.3, -2.0), seed, &opts.refinement, || {
            polls += 1;
            polls == 8
        }),
        Err(HarmonicPoleError::Solve(NonlinearEigError::Cancelled))
    ));
    let region = ResonanceRegion {
        center: C64::ZERO,
        real_radius: 4.0,
        imaginary_radius: 4.0,
    };
    assert!(model.resonance_spectrum(region, &opts, || true).is_err());
    for limited in [
        BeynOptions {
            max_dense_bytes: 0,
            ..opts
        },
        BeynOptions {
            max_work: 0,
            ..opts
        },
    ] {
        assert!(
            model
                .resonance_spectrum(region, &limited, || false)
                .is_err()
        );
    }
    let spectrum = model.resonance_spectrum(region, &opts, || false).unwrap();
    assert_eq!(spectrum.poles.len(), 2);
    assert_eq!(factor.response(&[C64::ONE], || false).unwrap(), before);
}
