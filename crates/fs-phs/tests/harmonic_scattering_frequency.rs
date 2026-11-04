//! G1 loaded physical tangents and G3 power-conjugate/source invariance.
use fs_math::c64::C64;
use fs_phs::circuit::{Branch, CircuitGraph, assemble_circuit};
use fs_phs::harmonic::scattering::{
    ScatteringBudget, ScatteringFrequencyResponse, ScatteringModel,
};
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
fn waves() -> ScatteringBudget {
    ScatteringBudget {
        maximum_ports: 8,
        maximum_payload_bytes: 1_000_000,
        maximum_work: 2_000_000,
    }
}
fn model(source: &PortHamiltonian, ports: &[(usize, f64)]) -> ScatteringModel {
    ScatteringModel::new(
        HarmonicModel::from_ode(source, budget(), || false).unwrap(),
        ports,
        waves(),
        || false,
    )
    .unwrap()
}
fn close(actual: C64, expected: C64) {
    assert!(
        (actual - expected).abs() <= 5e-10 * expected.abs().max(1e-12),
        "{actual:?} != {expected:?}"
    );
}
fn scalar(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 5e-10 * expected.abs().max(1.0),
        "{actual:e} != {expected:e}"
    );
}
fn balances(value: &ScatteringFrequencyResponse) {
    // A fixed incident wave reflected by a lossless system has exactly zero
    // power derivative while its phase derivative is finite. Admit floating
    // contraction error against the uncancelled amplitude products, not the
    // already cancelled real power slopes.
    let scale = value
        .response
        .waves
        .iter()
        .zip(&value.wave_derivatives)
        .map(|(w, d)| {
            (w.incident.abs() * d.incident.abs() + w.outgoing.abs() * d.outgoing.abs())
                / w.reference
        })
        .sum::<f64>()
        + (0..value.output_derivative.len())
            .map(|p| {
                0.5 * (value.physical_input_derivative[p].abs()
                    * value.response.physical.output[p].abs()
                    + value.response.physical_inputs[p].abs() * value.output_derivative[p].abs())
            })
            .sum::<f64>()
        + value.dissipated_power_derivative.abs()
        + 1e-30;
    assert!(value.power_residual_derivative.abs() <= 1e-9 * scale);
    assert!(value.wave_power_residual_derivative.abs() <= 1e-9 * scale);
    for wave in &value.wave_derivatives {
        assert!(wave.power_residual_w.abs() <= 1e-9 * scale);
    }
    assert!(value.relative_derivative_residual <= 2e-11);
}

#[test]
fn loaded_bare_pole_retains_reaction_derivatives_and_never_reports_boundary_heat() {
    let reference = 0.25;
    let a = C64::new(0.7, -0.3);
    for damping in [0.0, 0.4, 4.0] {
        let source = fs_phs::mass_spring_damper(1.0, 4.0, damping).unwrap();
        let loaded = model(&source, &[(0, reference)]);
        for omega in [0.1, 1.999, 2.0, 2.001, 10.0] {
            for da in [C64::ZERO, C64::new(0.2, 0.5)] {
                let factor = loaded.factor(omega, || false).unwrap();
                let value = factor
                    .response_with_frequency_derivative(&[a], Some(&[da]), || false)
                    .unwrap();
                assert_eq!(value.response, factor.response(&[a], || false).unwrap());
                let dynamic = C64::new(4.0 - omega * omega, -omega * (damping + 4.0));
                let ddynamic = C64::new(-2.0 * omega, -(damping + 4.0));
                let q = a.scale(8.0) / dynamic;
                let dq = (da.scale(8.0) - ddynamic * q) / dynamic;
                let velocity = C64::new(0.0, -omega) * q;
                let dv = C64::new(0.0, -1.0) * (q + dq.scale(omega));
                let u = (a.scale(2.0) - velocity).scale(4.0);
                let du = (da.scale(2.0) - dv).scale(4.0);
                close(value.state_derivative[0], dq);
                close(value.output_derivative[0], dv);
                close(value.physical_input_derivative[0], du);
                close(value.wave_derivatives[0].outgoing, dv - da);
                scalar(
                    value.dissipated_power_derivative,
                    damping * (velocity.conj() * dv).re,
                );
                scalar(
                    value.supplied_power_derivative,
                    0.5 * (du.conj() * velocity + u.conj() * dv).re,
                );
                scalar(
                    value.wave_derivatives[0].incident_power_w,
                    (a.conj() * da).re / reference,
                );
                if damping == 0.0 {
                    assert_eq!(value.dissipated_power_derivative, 0.0);
                }
                balances(&value);
            }
        }
    }
}

#[test]
fn mutual_ports_and_ordinary_drive_retain_the_full_two_body_derivative() {
    let mut j = vec![0.0; 16];
    j[2] = 1.0;
    j[7] = 1.0;
    j[8] = -1.0;
    j[13] = -1.0;
    let mut q = vec![0.0; 16];
    q[0] = 9.0;
    q[5] = 25.0;
    q[10] = 0.5;
    q[15] = 1.0 / 3.0;
    let mut r = vec![0.0; 16];
    r[10] = 0.5;
    r[11] = 0.1;
    r[14] = 0.1;
    r[15] = 0.8;
    let source = PortHamiltonian::new(
        4,
        3,
        j,
        r,
        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.2, 0.7, 0.4, 1.0, -0.2],
        Box::new(QuadraticStorage::new(q, 4).unwrap()),
    )
    .unwrap();
    let loaded = model(&source, &[(1, 0.35), (0, 0.2)]);
    let input = [C64::new(0.3, 0.1), C64::new(-0.2, 0.4), C64::new(0.7, -0.1)];
    let derivative = [C64::new(-0.1, 0.4), C64::new(0.5, -0.2), C64::new(0.3, 0.1)];
    let w = 2.7;
    let value = loaded
        .factor(w, || false)
        .unwrap()
        .response_with_frequency_derivative(&input, Some(&derivative), || false)
        .unwrap();
    let c00 = 0.5 + 1.0 / 0.2 + 0.04 / 0.35;
    let c01 = 0.1 + 0.4 / 0.2 + 0.2 / 0.35;
    let c11 = 0.8 + 0.16 / 0.2 + 1.0 / 0.35;
    let a = C64::new(9.0 - 2.0 * w * w, -w * c00);
    let b = C64::new(0.0, -w * c01);
    let d = C64::new(25.0 - 3.0 * w * w, -w * c11);
    let da = C64::new(-4.0 * w, -c00);
    let db = C64::new(0.0, -c01);
    let dd = C64::new(-6.0 * w, -c11);
    let forces = |u: &[C64]| {
        [
            u[0].scale(10.0) + u[1].scale(0.4 / 0.35) + u[2].scale(0.7),
            u[0].scale(4.0) + u[1].scale(2.0 / 0.35) - u[2].scale(0.2),
        ]
    };
    let f = forces(&input);
    let df = forces(&derivative);
    let det = a * d - b * b;
    let q0 = (d * f[0] - b * f[1]) / det;
    let q1 = (a * f[1] - b * f[0]) / det;
    let df0 = df[0] - da * q0 - db * q1;
    let df1 = df[1] - db * q0 - dd * q1;
    let dq0 = (d * df0 - b * df1) / det;
    let dq1 = (a * df1 - b * df0) / det;
    close(value.state_derivative[0], dq0);
    close(value.state_derivative[1], dq1);
    let v0 = C64::new(0.0, -w) * q0;
    let v1 = C64::new(0.0, -w) * q1;
    let dv0 = C64::new(0.0, -1.0) * (q0 + dq0.scale(w));
    let dv1 = C64::new(0.0, -1.0) * (q1 + dq1.scale(w));
    let dy = [
        dv0 + dv1.scale(0.4),
        dv0.scale(0.2) + dv1,
        dv0.scale(0.7) - dv1.scale(0.2),
    ];
    for p in 0..3 {
        close(value.output_derivative[p], dy[p]);
    }
    close(
        value.physical_input_derivative[0],
        (derivative[0].scale(2.0) - dy[0]).scale(5.0),
    );
    close(
        value.physical_input_derivative[1],
        (derivative[1].scale(2.0) - dy[1]).scale(1.0 / 0.35),
    );
    assert_eq!(value.physical_input_derivative[2], derivative[2]);
    assert_eq!(
        value
            .wave_derivatives
            .iter()
            .map(|x| x.original_port)
            .collect::<Vec<_>>(),
        vec![1, 0]
    );
    scalar(
        value.dissipated_power_derivative,
        (v0.conj() * (dv0.scale(0.5) + dv1.scale(0.1))
            + v1.conj() * (dv0.scale(0.1) + dv1.scale(0.8)))
        .re,
    );
    balances(&value);
}

#[test]
fn descriptor_reaction_derivative_uses_actual_terminal_voltage_at_series_resonance() {
    let (r, l, c, z) = (0.2, 0.5, 0.125, 2.0);
    let circuit = assemble_circuit(&CircuitGraph {
        node_count: 4,
        transformers: vec![],
        branches: vec![
            (1, 0, Branch::VoltageSource { port: 0 }),
            (1, 2, Branch::Inductor { henries: l }),
            (2, 3, Branch::Capacitor { farads: c }),
            (3, 0, Branch::Resistor { ohms: r }),
        ],
    })
    .unwrap();
    let loaded = ScatteringModel::new(
        HarmonicModel::from_descriptor(&circuit.system, budget(), || false).unwrap(),
        &[(0, z)],
        waves(),
        || false,
    )
    .unwrap();
    let (a, da) = (C64::new(0.8, -0.1), C64::new(-0.2, 0.3));
    for w in [1.2, 4.0, 8.0] {
        let value = loaded
            .factor(w, || false)
            .unwrap()
            .response_with_frequency_derivative(&[a], Some(&[da]), || false)
            .unwrap();
        let dynamic = C64::new(r + 1.0 / z, 1.0 / (w * c) - w * l);
        let dz = C64::new(0.0, -1.0 / (w * w * c) - l);
        let current = a.scale(2.0 / z) / dynamic;
        let di = (da.scale(2.0 / z) - dz * current) / dynamic;
        let du = (da.scale(2.0) - di).scale(1.0 / z);
        close(value.output_derivative[0], di);
        close(value.physical_input_derivative[0], du);
        close(value.state_derivative[circuit.node_potential_index[0]], du);
        close(value.state_derivative[circuit.flux_index[0]], di.scale(l));
        close(
            value.state_derivative[circuit.charge_index[0]],
            C64::new(0.0, 1.0 / w) * (di - current.scale(1.0 / w)),
        );
        scalar(
            value.dissipated_power_derivative,
            r * (current.conj() * di).re,
        );
        balances(&value);
    }
}

#[test]
fn loaded_maxwell_tangent_preserves_material_dispersion_history_and_its_own_heat() {
    let (mass, eq, branch, tau, z) = (0.4, 110.0, 300.0_f64, 0.025, 0.2);
    let source = fs_phs::mass_spring_damper(mass, eq, 0.0)
        .unwrap()
        .with_relaxation_branches(vec![RelaxationBranch {
            projection: vec![1.0, 0.0],
            stiffness: branch,
            relaxation_time_s: tau,
        }])
        .unwrap();
    let loaded = model(&source, &[(0, z)]);
    for w in [0.1, 10.0, 30.0, 100.0, 1000.0] {
        let value = loaded
            .factor(w, || false)
            .unwrap()
            .response_with_frequency_derivative(&[C64::ONE], None, || false)
            .unwrap();
        let st = C64::new(0.0, -w * tau);
        let dst = C64::new(0.0, -tau);
        let memory = C64::ONE + st;
        let dynamic = C64::new(eq - mass * w * w, -w / z) + (st / memory).scale(branch);
        let ddynamic =
            C64::new(-2.0 * mass * w, -1.0 / z) + (dst / (memory * memory)).scale(branch);
        let q = C64::from_re(2.0 / z) / dynamic;
        let dq = -ddynamic * q / dynamic;
        close(value.state_derivative[0], dq);
        close(
            value.state_derivative[2],
            (dq / memory - q * dst / (memory * memory)).scale(branch.sqrt()),
        );
        let d = 1.0 + (w * tau).powi(2);
        let heat_derivative = branch * tau * w / (d * d) * q.norm_sq()
            + branch * tau * w * w / d * (q.conj() * dq).re;
        scalar(value.dissipated_power_derivative, heat_derivative);
        balances(&value);
    }
}

#[test]
fn power_conjugate_port_units_preserve_physical_loaded_frequency_tangents() {
    let evaluate = |scale: f64| {
        let source = PortHamiltonian::new(
            2,
            1,
            vec![0.0, 1.0, -1.0, 0.0],
            vec![0.0, 0.0, 0.0, 0.2],
            vec![0.0, 1.0 / scale],
            Box::new(QuadraticStorage::new(vec![8.0, 0.0, 0.0, 0.5], 2).unwrap()),
        )
        .unwrap();
        let loaded = model(&source, &[(0, 0.5 / (scale * scale))]);
        loaded
            .factor(1.7, || false)
            .unwrap()
            .response_with_frequency_derivative(
                &[C64::new(0.7, -0.2).scale(1.0 / scale)],
                Some(&[C64::new(0.1, 0.3).scale(1.0 / scale)]),
                || false,
            )
            .unwrap()
    };
    let reference = evaluate(1.0);
    for scale in [1e-40, 1e40] {
        let value = evaluate(scale);
        for i in 0..2 {
            close(value.state_derivative[i], reference.state_derivative[i]);
        }
        close(
            value.output_derivative[0].scale(scale),
            reference.output_derivative[0],
        );
        close(
            value.physical_input_derivative[0].scale(1.0 / scale),
            reference.physical_input_derivative[0],
        );
        close(
            value.wave_derivatives[0].outgoing.scale(scale),
            reference.wave_derivatives[0].outgoing,
        );
        scalar(
            value.dissipated_power_derivative,
            reference.dissipated_power_derivative,
        );
        scalar(
            value.supplied_power_derivative,
            reference.supplied_power_derivative,
        );
        balances(&value);
    }
}

#[test]
fn loaded_tangent_caps_cancellation_and_refusal_preserve_the_factor() {
    let source = fs_phs::mass_spring_damper(1.0, 4.0, 0.2).unwrap();
    let loaded = model(&source, &[(0, 0.25)]);
    let factor = loaded.factor(2.0, || false).unwrap();
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
    for bad in [
        vec![],
        vec![C64::from_re(f64::NAN)],
        vec![C64::from_re(f64::MAX)],
    ] {
        assert!(
            factor
                .response_with_frequency_derivative(&[C64::ONE], Some(&bad), || false)
                .is_err()
        );
    }
    assert_eq!(before, factor.response(&[C64::ONE], || false).unwrap());
    assert_eq!(
        expected,
        factor
            .response_with_frequency_derivative(&[C64::ONE], None, || false)
            .unwrap()
    );
    for tight in [
        ScatteringBudget {
            maximum_payload_bytes: 336,
            ..waves()
        },
        ScatteringBudget {
            maximum_work: 51,
            ..waves()
        },
    ] {
        let limited = ScatteringModel::new(
            HarmonicModel::from_ode(&source, budget(), || false).unwrap(),
            &[(0, 0.25)],
            tight,
            || false,
        )
        .unwrap();
        let factor = limited.factor(2.0, || false).unwrap();
        factor.response(&[C64::ONE], || false).unwrap();
        assert!(matches!(
            factor.response_with_frequency_derivative(&[C64::ONE], None, || false),
            Err(HarmonicError::Invalid { .. })
        ));
    }
}
