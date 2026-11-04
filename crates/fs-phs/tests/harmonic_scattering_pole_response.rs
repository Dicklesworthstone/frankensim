//! Original loaded-pole participation and physical parameter directions.
use fs_la::pole_response::PoleResponseOptions;
use fs_math::c64::C64;
use fs_phs::circuit::{Branch, CircuitGraph, assemble_circuit};
use fs_phs::harmonic::equilibrium::{
    EquilibriumLinearization, EquilibriumPoint, LinearizationBudget, NonlinearDissipation,
};
use fs_phs::harmonic::response::{HarmonicPoleResponseError, HarmonicPoleResponseOptions};
use fs_phs::harmonic::scattering::{
    ScatteringBudget, ScatteringModel, ScatteringPole, ScatteringPoleDerivative,
    ScatteringPoleParameterDirection as Direction,
};
use fs_phs::harmonic::{HarmonicBudget, HarmonicModel, HarmonicParameterDirection, PoleOptions};
use fs_phs::{PortHamiltonian, QuadraticStorage, RelaxationBranch};
use std::cell::Cell;

fn budget() -> HarmonicBudget {
    HarmonicBudget {
        max_dimension: 32,
        max_dense_bytes: 2_000_000,
        max_dense_work: 10_000_000,
        relative_residual_tolerance: 1e-10,
    }
}
fn options(scale: f64) -> HarmonicPoleResponseOptions {
    HarmonicPoleResponseOptions {
        numerical: PoleResponseOptions {
            spectral_scale: scale,
            pole_residual_tolerance: 1e-12,
            pole_step_tolerance: 1e-12,
            linear_residual_tolerance: 1e-11,
            max_dimension: 32,
            max_vectors: 8,
            max_dense_bytes: 2_000_000,
            max_work: 10_000_000,
        },
        max_physical_bytes: 2_000_000,
        max_physical_work: 10_000_000,
    }
}
fn close(a: C64, b: C64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * (1.0 + b.abs()),
        "{a:?} != {b:?}"
    );
}
fn scalar(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * (1.0 + b.abs()),
        "{a:e} != {b:e}"
    );
}
fn loaded(model: HarmonicModel, ports: &[(usize, f64)]) -> ScatteringModel {
    ScatteringModel::new(
        model,
        ports,
        ScatteringBudget {
            maximum_ports: 8,
            maximum_payload_bytes: 1_000_000,
            maximum_work: 2_000_000,
        },
        || false,
    )
    .unwrap()
}
fn ode(source: &PortHamiltonian, ports: &[(usize, f64)]) -> ScatteringModel {
    loaded(
        HarmonicModel::from_ode(source, budget(), || false).unwrap(),
        ports,
    )
}
fn resolve(model: &ScatteringModel, s: C64, seed: &[C64], scale: f64) -> ScatteringPole {
    model
        .resonance_pole(
            s,
            seed,
            &PoleOptions {
                spectral_scale: scale,
                step_tolerance: 1e-13,
                residual_tolerance: 1e-13,
                max_iterations: 40,
                max_dimension: 32,
                max_dense_bytes: 2_000_000,
                max_work: 10_000_000,
            },
            || false,
        )
        .unwrap()
}
fn mechanical(
    m: f64,
    k: f64,
    c: f64,
    g: [f64; 3],
    z: [f64; 2],
) -> (ScatteringModel, ScatteringPole) {
    let source = PortHamiltonian::new(
        2,
        3,
        vec![0., 1., -1., 0.],
        vec![0., 0., 0., c],
        vec![0., 0., 0., g[0], g[1], g[2]],
        Box::new(QuadraticStorage::new(vec![k, 0., 0., 1. / m], 2).unwrap()),
    )
    .unwrap();
    // Two selected boundaries in reversed source order; the middle port is ordinary.
    let model = ode(&source, &[(2, z[1]), (0, z[0])]);
    let gamma = c + g[0] * g[0] / z[0] + g[2] * g[2] / z[1];
    let decay = gamma / (2. * m);
    let s = C64::new(-decay, -(k / m - decay * decay).sqrt());
    let pole = resolve(&model, s, &[C64::ONE, s.scale(m)], (k / m).sqrt());
    (model, pole)
}
fn balance(d: &ScatteringPoleDerivative) {
    let scale = 1. + d.internal_dissipated_power_w.abs() + d.outgoing_power_w.abs();
    assert!(d.relative_residual <= 1e-10);
    assert!(d.linear_residual <= 1e-11);
    assert!(d.decay_power_residual_w.abs() <= 2e-9 * scale, "{d:?}");
    assert!(d.physical_power_residual_w.abs() <= 2e-9 * scale, "{d:?}");
    for wave in &d.waves {
        assert!(wave.power_residual_w.abs() <= 2e-9 * scale);
    }
}

#[test]
fn mixed_incident_and_ordinary_forcing_has_exact_physical_laurent_residues() {
    let (m, k, c, g, z) = (2.5, 1600., 0.6, [0.7, 0.2, -0.4], [0.5, 0.8]);
    let (model, pole) = mechanical(m, k, c, g, z);
    let prepared = model
        .prepare_pole_response(&pole, options(26.), || false)
        .unwrap();
    let s = prepared.laplace_pole_per_s();
    let gamma = c + g[0] * g[0] / z[0] + g[2] * g[2] / z[1];
    let input = [C64::new(0.8, -0.2), C64::new(-0.3, 0.4), C64::new(0.2, 0.1)];
    for scale in [1e-50, 1., 1e50] {
        let source = input.map(|v| v.scale(scale));
        let force = source[0].scale(2. * g[0] / z[0])
            + source[1].scale(g[1])
            + source[2].scale(2. * g[2] / z[1]);
        let q = force / (s.scale(2. * m) + C64::from_re(gamma));
        let result = prepared
            .port_residues(&[&source], || false)
            .unwrap()
            .remove(0);
        close(
            result.physical.state[0].scale(1. / scale),
            q.scale(1. / scale),
            2e-11,
        );
        close(
            result.physical.state[1].scale(1. / scale),
            (q * s).scale(m / scale),
            2e-11,
        );
        for p in 0..3 {
            close(
                result.physical.output[p].scale(1. / scale),
                (q * s).scale(g[p] / scale),
                2e-11,
            );
        }
        close(
            result.physical_inputs[0],
            result.physical.output[0].scale(-1. / z[0]),
            2e-11,
        );
        close(
            result.physical_inputs[2],
            result.physical.output[2].scale(-1. / z[1]),
            2e-11,
        );
        assert_eq!(result.physical_inputs[1], C64::ZERO);
        assert_eq!(
            result.outgoing,
            vec![result.physical.output[2], result.physical.output[0]]
        );
        for (value, mode) in result.physical.state.iter().zip(prepared.right_state()) {
            close(
                value.scale(1. / scale),
                (*mode * result.physical.participation).scale(1. / scale),
                2e-11,
            );
        }
        assert!(result.relative_residual <= 1e-10);
    }
    // Analytic regular incident feedthrough must disappear in (s-s0)*response.
    let residue = prepared
        .port_residues(&[&input], || false)
        .unwrap()
        .remove(0);
    let force =
        input[0].scale(2. * g[0] / z[0]) + input[1].scale(g[1]) + input[2].scale(2. * g[2] / z[1]);
    let mut errors = Vec::new();
    for h in [0.2, 0.1, 0.05] {
        let delta = C64::new(h, 0.4 * h);
        let at = s + delta;
        let q = force / ((at * at).scale(m) + at.scale(gamma) + C64::from_re(k));
        let actual_input = (input[0].scale(2.) - (at * q).scale(g[0])).scale(1. / z[0]);
        let outgoing = (at * q).scale(g[2]) - input[2];
        errors.push(
            (delta * actual_input - residue.physical_inputs[0]).abs()
                + (delta * outgoing - residue.outgoing[0]).abs(),
        );
    }
    assert!(errors[0] / errors[1] > 1.9 && errors[1] / errors[2] > 1.9);
}

#[test]
fn original_material_port_and_reference_directions_tune_poles_without_boundary_heat() {
    let (m, k, g, z) = (2.0, 200.0, [0.8, 0.3, 0.6], [0.4, 1.2]);
    for c in [0., 4.] {
        let (model, pole) = mechanical(m, k, c, g, z);
        let prepared = model
            .prepare_pole_response(&pole, options(10.), || false)
            .unwrap();
        assert_eq!(prepared.normalization_pivot(), 1);
        let s = prepared.laplace_pole_per_s();
        let q = (s.scale(m)).recip();
        let gamma = c + g[0] * g[0] / z[0] + g[2] * g[2] / z[1];
        let mode = prepared.mode_observation();
        scalar(
            mode.stored_energy_j,
            0.5 * (k * q.norm_sq() + 1. / m),
            1e-11,
        );
        scalar(mode.internal_dissipated_power_w, c / (m * m), 1e-11);
        scalar(mode.outgoing_power_w, (gamma - c) / (m * m), 1e-11);
        assert!(mode.decay_power_residual_w.abs() < 1e-10);
        for (dm, dk, dc, dg, dz) in [
            (1., 0., 0., [0., 0., 0.], [0., 0.]),
            (0., 1., 0., [0., 0., 0.], [0., 0.]),
            (0., 0., 1., [0., 0., 0.], [0., 0.]),
            (0., 0., 0., [1., 0., 0.], [0., 0.]),
            (0., 0., 0., [0., 1., 0.], [0., 0.]),
            (0., 0., 0., [0., 0., 0.], [1., 0.]),
            (0.2, -0.3, 0.1, [-0.4, 0.2, 0.3], [0.2, -0.1]),
        ] {
            let dq_matrix = [dk, 0., 0., -dm / (m * m)];
            let dr = [0., 0., 0., dc];
            let dp = [0., 0., 0., dg[0], dg[1], dg[2]];
            let references = [dz[1], dz[0]];
            let derivative = prepared
                .parameter_derivatives(
                    &[Direction {
                        physical: HarmonicParameterDirection {
                            storage: Some(&dq_matrix),
                            dissipation: Some(&dr),
                            ports: Some(&dp),
                            ..Default::default()
                        },
                        reference_derivative: Some(&references),
                    }],
                    || false,
                )
                .unwrap()
                .remove(0);
            let dgamma = dc + 2. * g[0] * dg[0] / z[0] - g[0] * g[0] * dz[0] / (z[0] * z[0])
                + 2. * g[2] * dg[2] / z[1]
                - g[2] * g[2] * dz[1] / (z[1] * z[1]);
            let ds = -((s * s).scale(dm) + s.scale(dgamma) + C64::from_re(dk))
                / (s.scale(2. * m) + C64::from_re(gamma));
            let dq = -q * (C64::from_re(dm / m) + ds / s);
            close(derivative.laplace_pole_per_s, ds, 1e-10);
            close(
                derivative.angular_frequency_rad_s,
                C64::new(-ds.im, ds.re),
                1e-10,
            );
            close(derivative.state[0], dq, 1e-10);
            assert_eq!(derivative.state[1], C64::ZERO);
            close(derivative.effort[0], dq.scale(k) + q.scale(dk), 1e-10);
            close(derivative.effort[1], C64::from_re(-dm / (m * m)), 1e-10);
            for p in 0..3 {
                close(
                    derivative.output[p],
                    C64::from_re(dg[p] / m - g[p] * dm / (m * m)),
                    1e-10,
                );
            }
            scalar(
                derivative.stored_energy_j,
                k * (dq.conj() * q).re + 0.5 * dk * q.norm_sq() - 0.5 * dm / (m * m),
                1e-10,
            );
            scalar(
                derivative.internal_dissipated_power_w,
                dc / (m * m) - 2. * c * dm / (m * m * m),
                1e-10,
            );
            scalar(
                derivative.outgoing_power_w,
                (dgamma - dc) / (m * m) - 2. * (gamma - c) * dm / (m * m * m),
                1e-10,
            );
            assert_eq!(derivative.physical_inputs[1], C64::ZERO);
            for (at, p) in [2, 0].into_iter().enumerate() {
                let reference = if p == 0 { z[0] } else { z[1] };
                close(
                    derivative.physical_inputs[p],
                    -derivative.output[p].scale(1. / reference)
                        - mode.physical_inputs[p].scale(references[at] / reference),
                    1e-10,
                );
                assert_eq!(derivative.waves[at].original_port, p);
            }
            if c == 0. && dc == 0. {
                assert_eq!(derivative.internal_dissipated_power_w, 0.);
            }
            balance(&derivative);
        }
    }
}

#[test]
fn descriptor_rlc_retains_voltage_reactions_and_original_l_c_r_reference_chain_rule() {
    let (r, l, c, z) = (0.2_f64, 0.5_f64, 0.02_f64, 2.0_f64);
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
    let model = loaded(
        HarmonicModel::from_descriptor(&circuit.system, budget(), || false).unwrap(),
        &[(0, z)],
    );
    let gamma = r + 1. / z;
    let s = C64::new(
        -gamma / (2. * l),
        -(1. / (l * c) - (gamma / (2. * l)).powi(2)).sqrt(),
    );
    let seed = model
        .factor(10., || false)
        .unwrap()
        .response(&[C64::ONE], || false)
        .unwrap()
        .physical
        .state;
    let pole = resolve(&model, s, &seed, 10.);
    let prepared = model
        .prepare_pole_response(&pole, options(10.), || false)
        .unwrap();
    let a = C64::new(0.7, -0.2);
    let residue = prepared.port_residues(&[&[a]], || false).unwrap().remove(0);
    let current = a.scale(2. / z) * s / (s.scale(2. * l) + C64::from_re(gamma));
    close(residue.physical.output[0], current, 1e-10);
    close(
        residue.physical.state[circuit.node_potential_index[0]],
        residue.physical_inputs[0],
        1e-10,
    );
    let (dl, dc, dr, dz) = (0.2, 0.003, 0.1, 0.3);
    let n = seed.len();
    let mut dq = vec![0.; n * n];
    dq[circuit.flux_index[0] * n + circuit.flux_index[0]] = -dl / (l * l);
    dq[circuit.charge_index[0] * n + circuit.charge_index[0]] = -dc / (c * c);
    let loss: Vec<_> = circuit
        .system
        .structure()
        .1
        .iter()
        .map(|v| -v * dr / r)
        .collect();
    let result = prepared
        .parameter_derivatives(
            &[Direction {
                physical: HarmonicParameterDirection {
                    storage: Some(&dq),
                    dissipation: Some(&loss),
                    ..Default::default()
                },
                reference_derivative: Some(&[dz]),
            }],
            || false,
        )
        .unwrap()
        .remove(0);
    let expected = -((s * s).scale(dl) + s.scale(dr - dz / (z * z)) - C64::from_re(dc / (c * c)))
        / (s.scale(2. * l) + C64::from_re(gamma));
    close(result.laplace_pole_per_s, expected, 2e-9);
    close(
        result.state[circuit.node_potential_index[0]],
        result.physical_inputs[0],
        1e-9,
    );
    balance(&result);
    // Algebraic reactions cannot acquire storage through a direction.
    let algebraic = circuit.node_potential_index[0];
    dq[algebraic * n + algebraic] = 1.;
    assert!(
        prepared
            .parameter_derivatives(
                &[Direction {
                    physical: HarmonicParameterDirection {
                        storage: Some(&dq),
                        ..Default::default()
                    },
                    ..Default::default()
                }],
                || false
            )
            .is_err()
    );
}

#[test]
fn maxwell_oscillatory_and_relaxation_poles_keep_complete_causal_residues_and_tangents() {
    let (kb, tau, z) = (50_000.0_f64, 1. / 300., 0.5);
    let source = fs_phs::mass_spring_damper(1., 60_000., 0.)
        .unwrap()
        .with_relaxation_branches(vec![RelaxationBranch {
            projection: vec![1., 0.],
            stiffness: kb,
            relaxation_time_s: tau,
        }])
        .unwrap();
    let model = ode(&source, &[(0, z)]);
    let dq = [
        1.,
        0.,
        -0.5 / kb.sqrt(),
        0.,
        0.,
        0.,
        -0.5 / kb.sqrt(),
        0.,
        0.,
    ];
    let mut dr = [0.; 9];
    dr[8] = -1. / (tau * tau);
    for seed in [C64::new(-50., -87_500.0_f64.sqrt()), C64::from_re(-200.)] {
        let pole = resolve(
            &model,
            seed,
            &[
                C64::ONE,
                seed,
                C64::from_re(kb.sqrt()) / (C64::ONE + seed.scale(tau)),
            ],
            400.,
        );
        let prepared = model
            .prepare_pole_response(&pole, options(400.), || false)
            .unwrap();
        let s = prepared.laplace_pole_per_s();
        let arm = C64::ONE + s.scale(tau);
        let slope = s.scale(2.) + C64::from_re(1. / z) + C64::from_re(kb * tau) / (arm * arm);
        let residue = prepared
            .port_residues(&[&[C64::ONE]], || false)
            .unwrap()
            .remove(0);
        let q = C64::from_re(2. / z) / slope;
        close(residue.physical.state[0], q, 1e-10);
        close(residue.physical.state[2], q.scale(kb.sqrt()) / arm, 1e-10);
        close(residue.physical.output[0], s * q, 1e-10);
        for (direction, expected) in [
            (
                Direction {
                    physical: HarmonicParameterDirection {
                        storage: Some(&dq),
                        ..Default::default()
                    },
                    ..Default::default()
                },
                -s.scale(tau) / arm / slope,
            ),
            (
                Direction {
                    physical: HarmonicParameterDirection {
                        dissipation: Some(&dr),
                        ..Default::default()
                    },
                    ..Default::default()
                },
                -s.scale(kb) / (arm * arm) / slope,
            ),
            (
                Direction {
                    reference_derivative: Some(&[1.]),
                    ..Default::default()
                },
                s.scale(1. / (z * z)) / slope,
            ),
        ] {
            let d = prepared
                .parameter_derivatives(&[direction], || false)
                .unwrap()
                .remove(0);
            close(d.laplace_pole_per_s, expected, 2e-8);
            assert_eq!(d.state.len(), 3);
            balance(&d);
        }
    }
}

#[test]
fn dark_internal_mode_has_zero_residue_but_a_physical_visibility_tangent() {
    let source = PortHamiltonian::new(
        4,
        1,
        vec![
            0., 0., 1., 0., 0., 0., 0., 1., -1., 0., 0., 0., 0., -1., 0., 0.,
        ],
        vec![
            0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.2, 0., 0., 0., 0., 0.4,
        ],
        vec![0., 0., 1., 0.],
        Box::new(
            QuadraticStorage::new(
                vec![
                    100., 0., 0., 0., 0., 225., 0., 0., 0., 0., 1., 0., 0., 0., 0., 1.,
                ],
                4,
            )
            .unwrap(),
        ),
    )
    .unwrap();
    let model = ode(&source, &[(0, 1.)]);
    let s = C64::new(-0.2, -(225.0_f64 - 0.04).sqrt());
    let pole = resolve(&model, s, &[C64::ZERO, C64::ONE, C64::ZERO, s], 15.);
    let prepared = model
        .prepare_pole_response(&pole, options(15.), || false)
        .unwrap();
    assert_eq!(prepared.mode_output(), &[C64::ZERO]);
    let r = prepared
        .port_residues(&[&[C64::ONE]], || false)
        .unwrap()
        .remove(0);
    assert_eq!(r.physical.participation, C64::ZERO);
    assert!(r.physical.state.iter().all(|v| *v == C64::ZERO));
    let d = prepared
        .parameter_derivatives(
            &[Direction {
                physical: HarmonicParameterDirection {
                    ports: Some(&[0., 0., 0., 1.]),
                    ..Default::default()
                },
                ..Default::default()
            }],
            || false,
        )
        .unwrap()
        .remove(0);
    close(d.laplace_pole_per_s, C64::ZERO, 1e-12);
    let bright = s * s + s.scale(1.2) + C64::from_re(100.);
    close(d.output[0], C64::ONE - s / bright, 1e-11);
    assert!(d.output[0].abs() > 0.5);
    balance(&d);
}

#[test]
fn biased_nonlinear_drag_keeps_total_dq_dt_and_incremental_material_work() {
    let (mass, k, linear, beta, p0, z) = (2.0_f64, 200.0_f64, 0.1_f64, 0.7_f64, 0.4_f64, 1.5_f64);
    let source = PortHamiltonian::new(
        2,
        2,
        vec![0., 1., -1., 0.],
        vec![0., 0., 0., linear],
        vec![1., 0., 0., 1.],
        Box::new(QuadraticStorage::new(vec![k, 0., 0., 1. / mass], 2).unwrap()),
    )
    .unwrap();
    let v0 = p0 / mass;
    let hessian = |_: &[f64], v: &[f64], out: &mut [f64]| {
        out[0] = k * v[0];
        out[1] = v[1] / mass;
        true
    };
    let flow = |_: &[f64], e: &[f64], out: &mut [f64]| {
        out[0] = 0.;
        out[1] = beta * e[1].powi(3);
        true
    };
    let tangent = |_: &[f64], e: &[f64], _: &[f64], de: &[f64], out: &mut [f64]| {
        out[0] = 0.;
        out[1] = 3. * beta * e[1] * e[1] * de[1];
        true
    };
    let local = EquilibriumLinearization::from_dissipative_ode(
        &source,
        EquilibriumPoint {
            state: &[0., p0],
            input: &[-v0, linear * v0 + beta * v0.powi(3)],
            balance_scales: &[1., 1.],
            relative_tolerance: 1e-12,
        },
        LinearizationBudget {
            harmonic: budget(),
            max_hessian_actions: 2,
        },
        &hessian,
        NonlinearDissipation {
            flow: &flow,
            tangent: &tangent,
            max_tangent_actions: 2,
        },
        || false,
    )
    .unwrap();
    let model = loaded(local.into_model(), &[(1, z)]);
    let material = linear + 3. * beta * v0 * v0;
    let gamma = material + 1. / z;
    let s = C64::new(
        -gamma / (2. * mass),
        -(k / mass - (gamma / (2. * mass)).powi(2)).sqrt(),
    );
    let pole = resolve(&model, s, &[C64::ONE, s.scale(mass)], 10.);
    let prepared = model
        .prepare_pole_response(&pole, options(10.), || false)
        .unwrap();
    for (dm, dbeta) in [(1., 0.), (0., 1.)] {
        let dq = [0., 0., 0., -dm / (mass * mass)];
        let dt = [
            0.,
            0.,
            0.,
            3. * v0 * v0 * dbeta / mass - 9. * beta * p0 * p0 * dm / mass.powi(4),
        ];
        let dmaterial = 3. * v0 * v0 * dbeta - 6. * beta * p0 * p0 * dm / mass.powi(3);
        let d = prepared
            .parameter_derivatives(
                &[Direction {
                    physical: HarmonicParameterDirection {
                        storage: Some(&dq),
                        resisting_flow: Some(&dt),
                        ..Default::default()
                    },
                    ..Default::default()
                }],
                || false,
            )
            .unwrap()
            .remove(0);
        close(
            d.laplace_pole_per_s,
            -((s * s).scale(dm) + s.scale(dmaterial)) / (s.scale(2. * mass) + C64::from_re(gamma)),
            1e-10,
        );
        scalar(
            d.internal_dissipated_power_w,
            dmaterial / (mass * mass) - 2. * material * dm / mass.powi(3),
            1e-10,
        );
        balance(&d);
    }
}

#[test]
fn stale_poles_batch_caps_and_cancelled_actions_leave_the_loaded_factor_reusable() {
    let (model, pole) = mechanical(2., 200., 0.2, [0.8, 0.3, 0.6], [0.4, 1.2]);
    let mut opt = options(10.);
    opt.numerical.max_vectors = 2;
    let mut stale = pole.clone();
    stale.laplace_pole_per_s = stale.laplace_pole_per_s + C64::new(0.1, 0.2);
    assert!(model.prepare_pole_response(&stale, opt, || false).is_err());
    let mut cached = pole.clone();
    cached.effort.fill(C64::new(f64::NAN, 0.));
    cached.stored_energy_j = f64::NAN;
    let prepared = model.prepare_pole_response(&cached, opt, || false).unwrap();
    assert!(prepared.mode_observation().stored_energy_j.is_finite());
    for (bytes, work) in [(1, 10_000_000), (2_000_000, 1)] {
        let limited = HarmonicPoleResponseOptions {
            max_physical_bytes: bytes,
            max_physical_work: work,
            ..opt
        };
        assert!(
            model
                .prepare_pole_response(&pole, limited, || false)
                .is_err()
        );
    }
    let input = [C64::ONE, C64::ZERO, C64::ZERO];
    assert!(
        prepared
            .port_residues(&[&input, &input, &input], || false)
            .is_err()
    );
    assert!(
        prepared
            .parameter_derivatives(&[Direction::default(); 3], || false)
            .is_err()
    );
    for references in [&[][..], &[f64::NAN, 0.][..]] {
        assert!(
            prepared
                .parameter_derivatives(
                    &[Direction {
                        reference_derivative: Some(references),
                        ..Default::default()
                    }],
                    || false
                )
                .is_err()
        );
    }
    let direction = Direction {
        reference_derivative: Some(&[0.3, -0.2]),
        ..Default::default()
    };
    let count = Cell::new(0);
    let before = prepared
        .parameter_derivatives(&[direction, direction], || {
            count.set(count.get() + 1);
            false
        })
        .unwrap();
    let total = count.get();
    for cut in [1, total / 2, total] {
        let count = Cell::new(0);
        assert!(
            prepared
                .parameter_derivatives(&[direction, direction], || {
                    count.set(count.get() + 1);
                    count.get() >= cut
                })
                .is_err()
        );
    }
    let after = prepared
        .parameter_derivatives(&[direction], || false)
        .unwrap();
    assert_eq!(before[0].state, after[0].state);
    assert_eq!(before[0].waves, after[0].waves);
    assert!(matches!(
        prepared.port_residues(&[&input], || true),
        Err(HarmonicPoleResponseError::Cancelled)
    ));
    assert!(prepared.port_residues(&[&input], || false).is_ok());
}
