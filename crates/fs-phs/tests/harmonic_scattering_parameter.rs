//! G1 physical loaded tuning and G3 conjugate-unit/reference invariance.
use fs_math::c64::C64;
use fs_phs::circuit::{Branch, CircuitGraph, assemble_circuit};
use fs_phs::harmonic::scattering::{
    ScatteringBudget, ScatteringModel, ScatteringParameterResponse, ScatteringResponse,
};
use fs_phs::harmonic::{
    HarmonicBudget, HarmonicError, HarmonicModel, HarmonicParameterDirection as Direction,
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
fn loaded(source: &PortHamiltonian, ports: &[(usize, f64)]) -> ScatteringModel {
    ScatteringModel::new(
        HarmonicModel::from_ode(source, budget(), || false).unwrap(),
        ports,
        wave_budget(),
        || false,
    )
    .unwrap()
}
fn oscillator(mass: f64, stiffness: f64, damping: f64, coupling: f64) -> PortHamiltonian {
    PortHamiltonian::new(
        2,
        1,
        vec![0., 1., -1., 0.],
        vec![0., 0., 0., damping],
        vec![0., coupling],
        Box::new(QuadraticStorage::new(vec![stiffness, 0., 0., 1. / mass], 2).unwrap()),
    )
    .unwrap()
}
fn near(actual: C64, expected: C64, tol: f64) {
    assert!(
        (actual - expected).abs() <= tol * expected.abs().max(1e-25),
        "{actual:?} != {expected:?}"
    );
}
fn real(actual: f64, expected: f64, tol: f64) {
    assert!(
        (actual - expected).abs() <= tol * expected.abs().max(1e-12),
        "{actual:e} != {expected:e}"
    );
}
fn balance(value: &ScatteringParameterResponse) {
    let scale = value
        .response
        .waves
        .iter()
        .zip(&value.wave_derivatives)
        .map(|(w, d)| {
            (w.incident.abs() * d.incident.abs() + w.outgoing.abs() * d.outgoing.abs())
                / w.reference
                + (w.incident_power_w + w.outgoing_power_w)
                    * (d.reference_derivative / w.reference).abs()
        })
        .sum::<f64>()
        + value
            .response
            .physical
            .output
            .iter()
            .enumerate()
            .map(|(p, &y)| {
                0.5 * (y.abs() * value.physical_input_derivative[p].abs()
                    + value.response.physical_inputs[p].abs() * value.output_derivative[p].abs())
            })
            .sum::<f64>()
        + value.dissipated_power_derivative.abs()
        + 1e-25;
    assert!(
        value.power_residual_derivative.abs() <= 2e-9 * scale,
        "physical power {} / {scale}",
        value.power_residual_derivative
    );
    assert!(value.wave_power_residual_derivative.abs() <= 2e-9 * scale);
    for wave in &value.wave_derivatives {
        assert!(wave.power_residual_w.abs() <= 2e-9 * scale);
    }
    assert!(value.relative_derivative_residual <= 2e-11);
}

#[test]
fn loaded_physical_mass_loss_coupling_reference_and_drive_match_exact_mechanics() {
    let (mass, stiffness, coupling, reference) = (2.0, 8.0, 0.7, 0.35);
    let a = C64::new(0.4, -0.2);
    let da = C64::new(0.03, 0.02);
    for damping in [0.0, 0.2, 1e-5] {
        let system = loaded(
            &oscillator(mass, stiffness, damping, coupling),
            &[(0, reference)],
        );
        for omega in [0.5, 2.0, 3.5] {
            let s = C64::new(0.0, -omega);
            let dynamic = C64::new(
                stiffness - mass * omega * omega,
                -omega * (damping + coupling * coupling / reference),
            );
            let q = a.scale(2.0 * coupling / reference) / dynamic;
            let v = s * q;
            for (dm, dk, dc, dg, dz) in [
                (1., 0., 0., 0., 0.),
                (0., 1., 0., 0., 0.),
                (0., 0., 1., 0., 0.),
                (0., 0., 0., 1., 0.),
                (0., 0., 0., 0., -1.),
                (0.2, -0.3, 0.1, -0.4, 0.2),
            ] {
                let dq_matrix = [dk, 0., 0., -dm / (mass * mass)];
                let dr = [0., 0., 0., dc];
                let d_port = [0., dg];
                let factor = system.factor(omega, || false).unwrap();
                let result = factor
                    .response_with_parameter_derivative(
                        &[a],
                        Direction {
                            storage: Some(&dq_matrix),
                            dissipation: Some(&dr),
                            ports: Some(&d_port),
                            ..Default::default()
                        },
                        Some(&[da]),
                        Some(&[dz]),
                        || false,
                    )
                    .unwrap();
                assert_eq!(result.response, factor.response(&[a], || false).unwrap());
                let d_dynamic = C64::new(
                    dk - dm * omega * omega,
                    -omega
                        * (dc + 2. * coupling * dg / reference
                            - coupling * coupling * dz / (reference * reference)),
                );
                let d_force = da.scale(2. * coupling / reference)
                    + a.scale(2. * dg / reference - 2. * coupling * dz / (reference * reference));
                let dq = (d_force - d_dynamic * q) / dynamic;
                let dv = s * dq;
                let y = v.scale(coupling);
                let dy = dv.scale(coupling) + v.scale(dg);
                let u = (a.scale(2.) - y).scale(1. / reference);
                let du = (da.scale(2.) - dy - u.scale(dz)).scale(1. / reference);
                near(result.state_derivative[0], dq, 2e-10);
                near(
                    result.state_derivative[1],
                    dv.scale(mass) + v.scale(dm),
                    2e-10,
                );
                near(
                    result.effort_derivative[0],
                    dq.scale(stiffness) + q.scale(dk),
                    2e-10,
                );
                near(result.effort_derivative[1], dv, 2e-10);
                near(result.output_derivative[0], dy, 2e-10);
                // At a lossless bare pole the actual force and its derivative
                // can vanish by exact boundary cancellation. Compare the two
                // independent formulas against their uncancelled wave scale.
                let force_roundoff =
                    128. * f64::EPSILON * (2. * da.abs() + dy.abs() + u.abs() * dz.abs())
                        / reference;
                assert!(
                    (result.physical_input_derivative[0] - du).abs()
                        <= 2e-10 * du.abs() + force_roundoff
                );
                near(result.wave_derivatives[0].outgoing, dy - da, 2e-10);
                real(
                    result.mean_stored_energy_derivative,
                    0.25 * (dk * q.norm_sq() + dm * v.norm_sq())
                        + 0.5 * (stiffness * (q.conj() * dq).re + mass * (v.conj() * dv).re),
                    2e-9,
                );
                real(
                    result.dissipated_power_derivative,
                    0.5 * dc * v.norm_sq() + damping * (v.conj() * dv).re,
                    2e-9,
                );
                real(
                    result.wave_derivatives[0].incident_power_w,
                    (a.conj() * da).re / reference
                        - 0.5 * a.norm_sq() * dz / (reference * reference),
                    2e-10,
                );
                balance(&result);
            }
        }
    }
}

fn compare_fd(
    value: &ScatteringParameterResponse,
    lo: &ScatteringResponse,
    hi: &ScatteringResponse,
    step: f64,
) {
    // Perturbation oracles amplify primal rounding by 1/step, including when
    // a physical power derivative cancels exactly. Scale that uncertainty in
    // the original quantity's units, independently of the computed tangent.
    let scalar_fd = |d: f64, a: f64, b: f64| {
        let expected = (b - a) * 0.5 / step;
        let rounding = 32. * f64::EPSILON * (a.abs() + b.abs()) / step;
        assert!(
            (d - expected).abs() <= 3e-6 * d.abs().max(expected.abs()) + rounding,
            "{d:e} != {expected:e}; rounding {rounding:e}"
        );
    };
    for (actual, a, b) in [
        (
            &value.state_derivative,
            &lo.physical.state,
            &hi.physical.state,
        ),
        (
            &value.effort_derivative,
            &lo.physical.effort,
            &hi.physical.effort,
        ),
        (
            &value.output_derivative,
            &lo.physical.output,
            &hi.physical.output,
        ),
        (
            &value.physical_input_derivative,
            &lo.physical_inputs,
            &hi.physical_inputs,
        ),
    ] {
        for ((&d, &a), &b) in actual.iter().zip(a).zip(b) {
            near(d, (b - a).scale(0.5 / step), 3e-6);
        }
    }
    for (d, a, b) in [
        (
            value.mean_stored_energy_derivative,
            lo.physical.mean_stored_energy_j,
            hi.physical.mean_stored_energy_j,
        ),
        (
            value.dissipated_power_derivative,
            lo.physical.dissipated_power_w,
            hi.physical.dissipated_power_w,
        ),
        (
            value.supplied_power_derivative,
            lo.physical.supplied_power_w,
            hi.physical.supplied_power_w,
        ),
    ] {
        scalar_fd(d, a, b);
    }
    for ((d, a), b) in value.wave_derivatives.iter().zip(&lo.waves).zip(&hi.waves) {
        near(
            d.outgoing,
            (b.outgoing - a.outgoing).scale(0.5 / step),
            3e-6,
        );
        scalar_fd(d.incident_power_w, a.incident_power_w, b.incident_power_w);
        scalar_fd(d.outgoing_power_w, a.outgoing_power_w, b.outgoing_power_w);
        scalar_fd(d.physical_power_w, a.physical_power_w, b.physical_power_w);
    }
    balance(value);
}

#[test]
fn mutual_materials_interconnection_and_mixed_source_slots_keep_the_entire_chain_rule() {
    let j = [
        0., 0.2, 1., 0., -0.2, 0., 0., 1., -1., 0., 0., 0.3, 0., -1., -0.3, 0.,
    ];
    let r = [
        0.1, 0., 0., 0., 0., 0.2, 0., 0., 0., 0., 0.5, 0.1, 0., 0., 0.1, 0.8,
    ];
    let q = [
        9.,
        0.4,
        0.,
        0.,
        0.4,
        25.,
        0.,
        0.,
        0.,
        0.,
        0.5,
        0.02,
        0.,
        0.,
        0.02,
        1. / 3.,
    ];
    let g = [0.1, 0., 0., 0., 0., 0.1, 1., 0.2, 0.7, 0.4, 1., -0.2];
    let dj = [
        0., 0.1, 0.2, 0., -0.1, 0., 0., -0.15, -0.2, 0., 0., 0.04, 0., 0.15, -0.04, 0.,
    ];
    let dr = [
        0.02, 0., 0., 0., 0., -0.03, 0., 0., 0., 0., 0.2, -0.04, 0., 0., -0.04, -0.1,
    ];
    let dq = [
        0.3, -0.1, 0., 0., -0.1, -0.4, 0., 0., 0., 0., 0.02, 0.01, 0., 0., 0.01, -0.03,
    ];
    let dg = [0., 0.1, 0., 0., 0., -0.1, 0.3, -0.2, 0.1, -0.1, 0.4, 0.2];
    let input = [C64::new(0.3, 0.1), C64::new(-0.2, 0.4), C64::new(0.7, -0.1)];
    let di = [
        C64::new(-0.1, 0.2),
        C64::new(0.03, -0.04),
        C64::new(0.2, 0.3),
    ];
    let ports = [(1, 0.35), (0, 0.2)];
    let dz = [0.1, -0.04];
    let at_parameter = |theta: f64| {
        let shift = |a: &[f64], d: &[f64]| a.iter().zip(d).map(|(a, d)| a + theta * d).collect();
        let source = PortHamiltonian::new(
            4,
            3,
            shift(&j, &dj),
            shift(&r, &dr),
            shift(&g, &dg),
            Box::new(QuadraticStorage::new(shift(&q, &dq), 4).unwrap()),
        )
        .unwrap();
        let selected: Vec<_> = ports
            .iter()
            .zip(&dz)
            .map(|(&(p, z), &dz)| (p, z + theta * dz))
            .collect();
        loaded(&source, &selected)
    };
    let point = at_parameter(0.);
    let result = point
        .factor(2.7, || false)
        .unwrap()
        .response_with_parameter_derivative(
            &input,
            Direction {
                interconnection: Some(&dj),
                dissipation: Some(&dr),
                storage: Some(&dq),
                ports: Some(&dg),
                resisting_flow: None,
            },
            Some(&di),
            Some(&dz),
            || false,
        )
        .unwrap();
    let h = 1e-5;
    let sample = |theta: f64| {
        at_parameter(theta)
            .factor(2.7, || false)
            .unwrap()
            .response(
                &input
                    .iter()
                    .zip(&di)
                    .map(|(&a, &d)| a + d.scale(theta))
                    .collect::<Vec<_>>(),
                || false,
            )
            .unwrap()
    };
    compare_fd(&result, &sample(-h), &sample(h), h);
    assert_eq!(result.physical_input_derivative[2], di[2]);
    assert_eq!(result.wave_derivatives[0].original_port, 1);
    assert_eq!(result.wave_derivatives[0].reference_derivative, dz[0]);
}

#[test]
fn loaded_rlc_parameters_retain_descriptor_reactions_and_real_resistor_heat() {
    let (resistance, inductance, capacitance, reference) = (0.02, 0.5, 0.125, 0.4);
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
    let n = circuit.system.state_dim();
    let flux = circuit.flux_index[0];
    let charge = circuit.charge_index[0];
    let resistor_node = circuit.node_potential_index[2];
    let loaded = ScatteringModel::new(
        HarmonicModel::from_descriptor(&circuit.system, budget(), || false).unwrap(),
        &[(0, reference)],
        wave_budget(),
        || false,
    )
    .unwrap();
    let a = C64::new(0.3, -0.7);
    let da = C64::new(0.02, 0.03);
    for omega in [1., 4., 7.] {
        for (dl, dc, dr, dz) in [
            (1., 0., 0., 0.),
            (0., 1., 0., 0.),
            (0., 0., 1., 0.),
            (0., 0., 0., 1.),
            (0.1, 0.2, -0.03, 0.04),
        ] {
            let mut dq = vec![0.; n * n];
            dq[flux * n + flux] = -dl / (inductance * inductance);
            dq[charge * n + charge] = -dc / (capacitance * capacitance);
            let mut d_r = vec![0.; n * n];
            d_r[resistor_node * n + resistor_node] = -dr / (resistance * resistance);
            let result = loaded
                .factor(omega, || false)
                .unwrap()
                .response_with_parameter_derivative(
                    &[a],
                    Direction {
                        storage: Some(&dq),
                        dissipation: Some(&d_r),
                        ..Default::default()
                    },
                    Some(&[da]),
                    Some(&[dz]),
                    || false,
                )
                .unwrap();
            let z = C64::new(resistance, 1. / (omega * capacitance) - omega * inductance);
            let dz_circuit = C64::new(dr, -omega * dl - dc / (omega * capacitance * capacitance));
            let denominator = C64::ONE + z.scale(reference);
            let current = a.scale(2.) / denominator;
            let di = (da.scale(2.) - current * (z.scale(dz) + dz_circuit.scale(reference)))
                / denominator;
            let voltage = (a.scale(2.) - current).scale(1. / reference);
            let dv = (da.scale(2.) - di - voltage.scale(dz)).scale(1. / reference);
            near(result.output_derivative[0], di, 4e-10);
            near(result.physical_input_derivative[0], dv, 4e-10);
            near(
                result.state_derivative[flux],
                di.scale(inductance) + current.scale(dl),
                4e-10,
            );
            near(
                result.state_derivative[charge],
                C64::new(0., 1. / omega) * di,
                4e-10,
            );
            near(
                result.state_derivative[circuit.node_potential_index[0]],
                dv,
                4e-10,
            );
            near(
                result.state_derivative[resistor_node],
                di.scale(resistance) + current.scale(dr),
                4e-10,
            );
            real(
                result.dissipated_power_derivative,
                0.5 * dr * current.norm_sq() + resistance * (current.conj() * di).re,
                2e-9,
            );
            balance(&result);
        }
    }
}

#[test]
fn loaded_maxwell_strength_and_time_keep_internal_history_storage_and_heat() {
    let (mass, k0, damping, branch, tau, reference) = (0.4, 110., 0.2, 300.0_f64, 0.025, 0.3);
    let source = fs_phs::mass_spring_damper(mass, k0, damping)
        .unwrap()
        .with_relaxation_branches(vec![RelaxationBranch {
            projection: vec![1., 0.],
            stiffness: branch,
            relaxation_time_s: tau,
        }])
        .unwrap();
    let loaded = loaded(&source, &[(0, reference)]);
    let a = C64::new(0.7, -0.3);
    for omega in [1., 30., 100.] {
        let s = C64::new(0., -omega);
        let memory = C64::ONE + s.scale(tau);
        let gain = s.scale(tau) / memory;
        let dynamic = C64::from_re(k0 - mass * omega * omega)
            + s.scale(damping + 1. / reference)
            + gain.scale(branch);
        let q = a.scale(2. / reference) / dynamic;
        let v = s * q;
        let history = q.scale(branch.sqrt()) / memory;
        let mismatch = q.scale(branch.sqrt()) - history;
        for (dk, dtau) in [(1., 0.), (0., 1.)] {
            let cross = -0.5 * dk / branch.sqrt();
            let dq_matrix = [dk, 0., cross, 0., 0., 0., cross, 0., 0.];
            let dr = [0., 0., 0., 0., 0., 0., 0., 0., -dtau / (tau * tau)];
            let result = loaded
                .factor(omega, || false)
                .unwrap()
                .response_with_parameter_derivative(
                    &[a],
                    Direction {
                        storage: Some(&dq_matrix),
                        dissipation: Some(&dr),
                        ..Default::default()
                    },
                    None,
                    None,
                    || false,
                )
                .unwrap();
            let d_dynamic = gain.scale(dk) + s.scale(branch * dtau) / (memory * memory);
            let dq = -d_dynamic * q / dynamic;
            let dv = s * dq;
            let dz = (dq.scale(branch.sqrt()) + q.scale(0.5 * dk / branch.sqrt())) / memory
                - q.scale(branch.sqrt() * dtau) * s / (memory * memory);
            let dmismatch = dq.scale(branch.sqrt()) + q.scale(0.5 * dk / branch.sqrt()) - dz;
            near(result.state_derivative[0], dq, 4e-10);
            near(result.state_derivative[2], dz, 4e-10);
            near(result.effort_derivative[2], -dmismatch, 4e-10);
            real(
                result.mean_stored_energy_derivative,
                0.5 * (k0 * (q.conj() * dq).re
                    + mass * (v.conj() * dv).re
                    + (mismatch.conj() * dmismatch).re),
                2e-9,
            );
            real(
                result.dissipated_power_derivative,
                damping * (v.conj() * dv).re + (mismatch.conj() * dmismatch).re / tau
                    - 0.5 * dtau * mismatch.norm_sq() / (tau * tau),
                2e-9,
            );
            balance(&result);
        }
    }
}

#[test]
fn conjugate_port_unit_change_preserves_all_physical_states_and_watts() {
    for scale in [1e-60, 1., 1e60] {
        let reference = 0.3 * scale * scale;
        let source = oscillator(2., 8., 0.2, scale);
        let loaded = loaded(&source, &[(0, reference)]);
        let a = C64::new(0.4 * scale, -0.2 * scale);
        let dg = [0., scale];
        let result = loaded
            .factor(2., || false)
            .unwrap()
            .response_with_parameter_derivative(
                &[a],
                Direction {
                    ports: Some(&dg),
                    ..Default::default()
                },
                Some(&[a]),
                Some(&[2. * reference]),
                || false,
            )
            .unwrap();
        for (d, x) in result
            .state_derivative
            .iter()
            .zip(&result.response.physical.state)
        {
            assert!(d.abs() <= 3e-12 * x.abs());
        }
        for (d, e) in result
            .effort_derivative
            .iter()
            .zip(&result.response.physical.effort)
        {
            assert!(d.abs() <= 3e-12 * e.abs());
        }
        near(
            result.output_derivative[0],
            result.response.physical.output[0],
            3e-12,
        );
        near(
            result.physical_input_derivative[0],
            -result.response.physical_inputs[0],
            3e-12,
        );
        near(
            result.wave_derivatives[0].outgoing,
            result.response.waves[0].outgoing,
            3e-12,
        );
        let norm =
            result.response.waves[0].incident_power_w + result.response.waves[0].outgoing_power_w;
        for power in [
            result.dissipated_power_derivative,
            result.supplied_power_derivative,
            result.wave_derivatives[0].incident_power_w,
            result.wave_derivatives[0].outgoing_power_w,
        ] {
            assert!(power.abs() <= 4e-12 * norm);
        }
        balance(&result);
    }
}

#[test]
fn zero_drive_source_slopes_and_refusals_preserve_the_reusable_closed_factor() {
    let source = oscillator(1., 4., 0.2, 1.);
    let loaded = loaded(&source, &[(0, 0.3)]);
    let factor = loaded.factor(2., || false).unwrap();
    let drive = C64::new(0.3, -0.2);
    let zero = factor
        .response_with_parameter_derivative(
            &[C64::ZERO],
            Direction::default(),
            Some(&[drive]),
            Some(&[0.1]),
            || false,
        )
        .unwrap();
    let driven = factor.response(&[drive], || false).unwrap();
    assert_eq!(zero.state_derivative, driven.physical.state);
    assert_eq!(zero.output_derivative, driven.physical.output);
    assert_eq!(zero.physical_input_derivative, driven.physical_inputs);
    assert_eq!(zero.mean_stored_energy_derivative, 0.);
    assert_eq!(zero.dissipated_power_derivative, 0.);
    assert_eq!(zero.wave_derivatives[0].incident_power_w, 0.);
    for references in [&[][..], &[f64::NAN][..], &[1., 2.][..]] {
        assert!(
            factor
                .response_with_parameter_derivative(
                    &[drive],
                    Direction::default(),
                    None,
                    Some(references),
                    || false
                )
                .is_err()
        );
    }
    assert!(
        factor
            .response_with_parameter_derivative(
                &[drive],
                Direction {
                    ports: Some(&[1.]),
                    ..Default::default()
                },
                None,
                None,
                || false
            )
            .is_err()
    );
    assert!(matches!(
        factor.response_with_parameter_derivative(
            &[drive],
            Direction::default(),
            None,
            None,
            || true
        ),
        Err(HarmonicError::Cancelled)
    ));
    let mut polls = 0;
    assert!(matches!(
        factor.response_with_parameter_derivative(
            &[drive],
            Direction::default(),
            None,
            None,
            || {
                polls += 1;
                polls == 15
            }
        ),
        Err(HarmonicError::Cancelled)
    ));
    assert_eq!(factor.response(&[drive], || false).unwrap(), driven);
    let limited = ScatteringModel::new(
        HarmonicModel::from_ode(&source, budget(), || false).unwrap(),
        &[(0, 0.3)],
        ScatteringBudget {
            maximum_ports: 1,
            maximum_payload_bytes: 400,
            maximum_work: 100,
        },
        || false,
    )
    .unwrap();
    assert!(
        limited
            .factor(2., || false)
            .unwrap()
            .response_with_parameter_derivative(&[drive], Direction::default(), None, None, || {
                false
            })
            .is_err()
    );
}
