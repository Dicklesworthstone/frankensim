//! G1 independent physical oracles for original-state wave feedback.
use fs_math::c64::C64;
use fs_phs::harmonic::scattering::{ScatteringBudget, ScatteringModel};
use fs_phs::harmonic::{HarmonicBudget, HarmonicError, HarmonicModel};
use fs_phs::{PortHamiltonian, QuadraticStorage, RelaxationBranch};

fn budget() -> HarmonicBudget {
    HarmonicBudget {
        max_dimension: 128,
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
fn close(a: C64, b: C64, scale: f64) {
    assert!(
        (a - b).abs() < 3e-10 * scale.max(a.abs()).max(b.abs()),
        "{a:?} != {b:?}"
    );
}
fn power(r: &fs_phs::harmonic::scattering::ScatteringResponse) {
    let scale = r
        .waves
        .iter()
        .map(|p| p.incident_power_w + p.outgoing_power_w)
        .sum::<f64>()
        + r.physical.dissipated_power_w.abs()
        + 1e-20;
    assert!(r.physical.power_residual_w.abs() < 2e-10 * scale, "{r:?}");
    assert!(r.wave_power_residual_w.abs() < 2e-10 * scale, "{r:?}");
    for p in &r.waves {
        assert!(p.power_residual_w.abs() < 2e-10 * scale);
    }
}

#[test]
fn bare_resonance_is_closed_before_factorization_and_reference_is_not_heat() {
    let incident = C64::new(0.7, -0.3);
    for damping in [0.0, 0.4] {
        let system = fs_phs::mass_spring_damper(1.0, 4.0, damping).unwrap();
        if damping == 0.0 {
            let bare = HarmonicModel::from_ode(&system, budget(), || false).unwrap();
            assert!(matches!(
                bare.factor(2.0, || false),
                Err(HarmonicError::Singular { .. })
            ));
        }
        let model = prepare(&system, &[(0, 0.25)]);
        for omega in [0.1, 1.999, 2.0, 2.001, 10.0] {
            let factor = model.factor(omega, || false).unwrap();
            let r = factor.response(&[incident], || false).unwrap();
            let displacement =
                incident.scale(8.0) / C64::new(4.0 - omega * omega, -omega * (damping + 4.0));
            let velocity = C64::new(0.0, -omega) * displacement;
            close(r.physical.state[0], displacement, 1e-12);
            close(r.physical.output[0], velocity, 1e-12);
            close(r.waves[0].outgoing, velocity - incident, 1e-12);
            close(
                r.physical_inputs[0],
                (incident.scale(2.0) - velocity).scale(4.0),
                1e-12,
            );
            let heat = 0.5 * damping * velocity.norm_sq();
            assert!((r.physical.dissipated_power_w - heat).abs() < 1e-12);
            if damping == 0.0 {
                assert_eq!(r.physical.dissipated_power_w, 0.0);
            }
            power(&r);
            let phase = C64::new(-0.2, 0.6);
            let second = factor.response(&[incident * phase], || false).unwrap();
            close(second.physical.state[0], displacement * phase, 1e-12);
        }
    }
}

#[test]
fn full_mutual_feedback_and_ordinary_forcing_match_an_independent_two_body_solve() {
    let n = 4;
    let mut j = vec![0.0; n * n];
    j[2] = 1.0;
    j[7] = 1.0;
    j[8] = -1.0;
    j[13] = -1.0;
    let mut q = vec![0.0; n * n];
    q[0] = 9.0;
    q[5] = 25.0;
    q[10] = 0.5;
    q[15] = 1.0 / 3.0;
    let mut r = vec![0.0; n * n];
    r[10] = 0.5;
    r[11] = 0.1;
    r[14] = 0.1;
    r[15] = 0.8;
    let g = vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.2, 0.7, 0.4, 1.0, -0.2];
    let system = PortHamiltonian::new(
        n,
        3,
        j,
        r,
        g,
        Box::new(QuadraticStorage::new(q, n).unwrap()),
    )
    .unwrap();
    let model = prepare(&system, &[(1, 0.35), (0, 0.2)]);
    assert_eq!(model.port_count(), 3);
    assert_eq!(model.dimension(), 4);
    let omega = 2.7;
    let factor = model.factor(omega, || false).unwrap();
    let input = [C64::new(0.3, 0.1), C64::new(-0.2, 0.4), C64::new(0.7, -0.1)];
    let response = factor.response(&input, || false).unwrap();
    let c00 = 0.5 + 1.0 / 0.2 + 0.04 / 0.35;
    let c01 = 0.1 + 0.4 / 0.2 + 0.2 / 0.35;
    let c11 = 0.8 + 0.16 / 0.2 + 1.0 / 0.35;
    let a = C64::new(9.0 - 2.0 * omega * omega, -omega * c00);
    let b = C64::new(0.0, -omega * c01);
    let d = C64::new(25.0 - 3.0 * omega * omega, -omega * c11);
    let f0 = input[0].scale(2.0 / 0.2) + input[1].scale(0.4 / 0.35) + input[2].scale(0.7);
    let f1 = input[0].scale(0.8 / 0.2) + input[1].scale(2.0 / 0.35) - input[2].scale(0.2);
    let det = a * d - b * b;
    close(response.physical.state[0], (d * f0 - b * f1) / det, 1e-12);
    close(response.physical.state[1], (a * f1 - b * f0) / det, 1e-12);
    assert_eq!(response.physical_inputs[2], input[2]);
    power(&response);
    let first = factor
        .response(&[C64::ONE, C64::ZERO, C64::ZERO], || false)
        .unwrap();
    let second = factor
        .response(&[C64::ZERO, C64::ONE, C64::ZERO], || false)
        .unwrap();
    // Power-normalized reciprocity with unequal port references and reordered slots.
    close(
        first.waves[0].outgoing.scale((0.2_f64 / 0.35).sqrt()),
        second.waves[1].outgoing.scale((0.35_f64 / 0.2).sqrt()),
        1e-12,
    );
}

#[test]
fn hereditary_storage_and_heat_survive_wave_closure_without_frequency_fitting() {
    let (mass, eq, branch, tau, z) = (0.4_f64, 110.0_f64, 300.0_f64, 0.025_f64, 0.2_f64);
    let system = fs_phs::mass_spring_damper(mass, eq, 0.0)
        .unwrap()
        .with_relaxation_branches(vec![RelaxationBranch {
            projection: vec![1.0, 0.0],
            stiffness: branch,
            relaxation_time_s: tau,
        }])
        .unwrap();
    let model = prepare(&system, &[(0, z)]);
    assert_eq!(model.dimension(), 3);
    for omega in [0.1, 10.0, 30.0, 100.0, 1000.0] {
        let result = model
            .factor(omega, || false)
            .unwrap()
            .response(&[C64::ONE], || false)
            .unwrap();
        let st = C64::new(0.0, -omega * tau);
        let dynamic =
            C64::new(eq - mass * omega * omega, -omega / z) + st.scale(branch) / (C64::ONE + st);
        let expected = C64::from_re(2.0 / z) / dynamic;
        close(result.physical.state[0], expected, 1e-20);
        close(
            result.physical.state[2],
            expected.scale(branch.sqrt()) / (C64::ONE + st),
            1e-20,
        );
        let heat =
            0.5 * branch * tau * omega * omega / (1.0 + (omega * tau).powi(2)) * expected.norm_sq();
        assert!((result.physical.dissipated_power_w / heat - 1.0).abs() < 1e-8);
        power(&result);
    }
}

#[test]
fn descriptor_rlc_keeps_original_constraints_and_physical_resistor_loss() {
    use fs_phs::circuit::{Branch, CircuitGraph, assemble_circuit};
    let (r, l, c, z) = (0.2_f64, 0.1_f64, 2.2e-7_f64, 2.0_f64);
    let graph = CircuitGraph {
        node_count: 4,
        branches: vec![
            (1, 0, Branch::VoltageSource { port: 0 }),
            (1, 2, Branch::Inductor { henries: l }),
            (2, 3, Branch::Capacitor { farads: c }),
            (3, 0, Branch::Resistor { ohms: r }),
        ],
        transformers: vec![],
    };
    let circuit = assemble_circuit(&graph).unwrap();
    let original = HarmonicModel::from_descriptor(&circuit.system, budget(), || false).unwrap();
    let model = ScatteringModel::new(original, &[(0, z)], wave_budget(), || false).unwrap();
    let wave = C64::new(0.8, -0.1);
    for ratio in [0.3, 1.0, 2.0] {
        let omega = ratio / (l * c).sqrt();
        let result = model
            .factor(omega, || false)
            .unwrap()
            .response(&[wave], || false)
            .unwrap();
        let current = wave.scale(2.0 / z) / C64::new(r + 1.0 / z, 1.0 / (omega * c) - omega * l);
        close(result.physical.output[0], current, 1e-12);
        close(
            result.physical.state[circuit.flux_index[0]].scale(1.0 / l),
            current,
            1e-12,
        );
        assert!((result.physical.dissipated_power_w - 0.5 * r * current.norm_sq()).abs() < 1e-10);
        power(&result);
    }
}

#[test]
fn admission_cancellation_and_reusable_factor_remain_explicit() {
    let system = fs_phs::mass_spring_damper(1.0, 4.0, 0.0).unwrap();
    let original = || HarmonicModel::from_ode(&system, budget(), || false).unwrap();
    for ports in [
        vec![],
        vec![(0, 0.0)],
        vec![(0, -1.0)],
        vec![(0, f64::NAN)],
        vec![(1, 1.0)],
        vec![(0, 1.0), (0, 2.0)],
    ] {
        assert!(ScatteringModel::new(original(), &ports, wave_budget(), || false).is_err());
    }
    for limited in [
        ScatteringBudget {
            maximum_ports: 0,
            ..wave_budget()
        },
        ScatteringBudget {
            maximum_payload_bytes: 0,
            ..wave_budget()
        },
        ScatteringBudget {
            maximum_work: 0,
            ..wave_budget()
        },
    ] {
        assert!(ScatteringModel::new(original(), &[(0, 1.0)], limited, || false).is_err());
    }
    assert!(matches!(
        ScatteringModel::new(original(), &[(0, 1.0)], wave_budget(), || true),
        Err(HarmonicError::Cancelled)
    ));
    let model = prepare(&system, &[(0, 1.0)]);
    assert!(matches!(
        model.factor(2.0, || true),
        Err(HarmonicError::Cancelled)
    ));
    let factor = model.factor(2.0, || false).unwrap();
    let first = factor.response(&[C64::ONE], || false).unwrap();
    assert!(factor.response(&[], || false).is_err());
    assert!(matches!(
        factor.response(&[C64::ONE], || true),
        Err(HarmonicError::Cancelled)
    ));
    assert_eq!(factor.response(&[C64::ONE], || false).unwrap(), first);
}
