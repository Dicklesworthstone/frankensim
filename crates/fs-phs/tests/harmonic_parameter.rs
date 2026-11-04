//! Physical parameter tangents of complete forced responses, not pole proxies.
use fs_math::c64::C64;
use fs_phs::circuit::{Branch, CircuitGraph, assemble_circuit};
use fs_phs::harmonic::equilibrium::{
    EquilibriumLinearization, EquilibriumPoint, LinearizationBudget, NonlinearDissipation,
};
use fs_phs::harmonic::{
    HarmonicBudget, HarmonicError, HarmonicModel, HarmonicParameterDirection as Direction,
};
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
fn mass_stiffness_and_loss_tangents_match_the_driven_high_q_oscillator() {
    let (mass, stiffness) = (2.0, 8.0);
    let (force, du) = (C64::new(0.7, -0.3), C64::new(0.03, -0.02));
    for damping in [0.6, 1e-5] {
        let system = fs_phs::mass_spring_damper(mass, stiffness, damping).unwrap();
        let model = HarmonicModel::from_ode(&system, budget(), || false).unwrap();
        for omega in [1.0, 2.0, 3.0] {
            let factor = model.factor(omega, || false).unwrap();
            let dynamic = C64::new(stiffness - mass * omega * omega, -damping * omega);
            let q = force / dynamic;
            let s = C64::new(0.0, -omega);
            let v = s * q;
            for (dm, dk, dc) in [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, -1.0)] {
                // Negative loss is a signed LOCAL direction, not a replacement R.
                let d_q = [dk, 0.0, 0.0, -dm / (mass * mass)];
                let d_r = [0.0, 0.0, 0.0, dc];
                let result = factor
                    .response_with_parameter_derivative(
                        &[force],
                        Direction {
                            storage: Some(&d_q),
                            dissipation: Some(&d_r),
                            ..Default::default()
                        },
                        Some(&[du]),
                        || false,
                    )
                    .unwrap();
                assert_eq!(
                    result.response,
                    factor.response(&[force], || false).unwrap()
                );
                let d_dynamic = C64::new(dk - dm * omega * omega, -dc * omega);
                let dq = (du - d_dynamic * q) / dynamic;
                let dv = s * dq;
                near(result.state_derivative[0], dq, 3e-9);
                near(
                    result.state_derivative[1],
                    dv.scale(mass) + v.scale(dm),
                    3e-9,
                );
                near(
                    result.effort_derivative[0],
                    dq.scale(stiffness) + q.scale(dk),
                    3e-9,
                );
                near(result.effort_derivative[1], dv, 3e-9);
                near(result.output_derivative[0], dv, 3e-9);
                let dh = 0.25 * (dk * q.norm_sq() + dm * v.norm_sq())
                    + 0.5 * (stiffness * (q.conj() * dq).re + mass * (v.conj() * dv).re);
                let heat = 0.5 * dc * v.norm_sq() + damping * (v.conj() * dv).re;
                let work = 0.5 * (du.conj() * v + force.conj() * dv).re;
                near_real(result.mean_stored_energy_derivative, dh, 2e-8);
                near_real(result.dissipated_power_derivative, heat, 2e-8);
                near_real(result.supplied_power_derivative, work, 2e-8);
                near(
                    result.output_logarithmic_derivative(0).unwrap().unwrap(),
                    dv / v,
                    3e-9,
                );
                assert!(result.relative_derivative_residual < 1e-12);
            }
        }
    }
}

#[test]
fn rlc_material_tangents_retain_the_full_kirchhoff_and_source_chain_rule() {
    let (resistance, inductance, capacitance) = (0.02, 0.5, 0.125);
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
    let model = HarmonicModel::from_descriptor(&circuit.system, budget(), || false).unwrap();
    let n = circuit.system.state_dim();
    let flux = circuit.flux_index[0];
    let charge = circuit.charge_index[0];
    let resistor_node = circuit.node_potential_index[2];
    let (voltage, du) = (C64::new(0.3, -0.7), C64::new(0.02, 0.03));
    for omega in [1.0, 4.0, 7.0] {
        for (dl, dc, dr) in [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)] {
            let mut d_q = vec![0.0; n * n];
            d_q[flux * n + flux] = -dl / (inductance * inductance);
            d_q[charge * n + charge] = -dc / (capacitance * capacitance);
            let mut d_r = vec![0.0; n * n];
            d_r[resistor_node * n + resistor_node] = -dr / (resistance * resistance);
            let result = model
                .factor(omega, || false)
                .unwrap()
                .response_with_parameter_derivative(
                    &[voltage],
                    Direction {
                        storage: Some(&d_q),
                        dissipation: Some(&d_r),
                        ..Default::default()
                    },
                    Some(&[du]),
                    || false,
                )
                .unwrap();
            let impedance = C64::new(resistance, 1.0 / (omega * capacitance) - omega * inductance);
            let d_impedance = C64::new(dr, -omega * dl - dc / (omega * capacitance * capacitance));
            let current = voltage / impedance;
            let di = (du - d_impedance * current) / impedance;
            let q = C64::new(0.0, 1.0 / omega) * current;
            let dq = C64::new(0.0, 1.0 / omega) * di;
            near(result.output_derivative[0], di, 3e-10);
            near(
                result.state_derivative[flux],
                di.scale(inductance) + current.scale(dl),
                3e-10,
            );
            near(result.state_derivative[charge], dq, 3e-10);
            near(result.effort_derivative[flux], di, 3e-10);
            near(
                result.state_derivative[circuit.node_potential_index[0]],
                du,
                3e-10,
            );
            near(
                result.state_derivative[resistor_node],
                di.scale(resistance) + current.scale(dr),
                3e-10,
            );
            near(
                result.state_derivative[circuit.node_potential_index[1]],
                di.scale(resistance) + current.scale(dr) + dq.scale(1.0 / capacitance)
                    - q.scale(dc / (capacitance * capacitance)),
                3e-10,
            );
            let dh = 0.25
                * (dl * current.norm_sq() - dc / (capacitance * capacitance) * q.norm_sq())
                + 0.5 * (inductance * (current.conj() * di).re + (q.conj() * dq).re / capacitance);
            near_real(result.mean_stored_energy_derivative, dh, 3e-9);
            near_real(
                result.dissipated_power_derivative,
                0.5 * dr * current.norm_sq() + resistance * (current.conj() * di).re,
                3e-9,
            );
            near_real(
                result.supplied_power_derivative,
                0.5 * (du.conj() * current + voltage.conj() * di).re,
                3e-9,
            );
            assert!(result.relative_derivative_residual < 1e-12);
        }
    }
    let mut bad = vec![0.0; n * n];
    bad[resistor_node * n + resistor_node] = 1.0;
    assert!(
        model
            .factor(4.0, || false)
            .unwrap()
            .response_with_parameter_derivative(
                &[voltage],
                Direction {
                    storage: Some(&bad),
                    ..Default::default()
                },
                None,
                || false
            )
            .is_err()
    );
}

#[test]
fn maxwell_strength_and_time_tangents_keep_real_internal_storage_and_heat() {
    let (mass, equilibrium, damping, branch, tau) = (0.4, 110.0, 0.2, 300.0_f64, 0.025);
    let source = fs_phs::mass_spring_damper(mass, equilibrium, damping)
        .unwrap()
        .with_relaxation_branches(vec![RelaxationBranch {
            projection: vec![1.0, 0.0],
            stiffness: branch,
            relaxation_time_s: tau,
        }])
        .unwrap();
    let model = HarmonicModel::from_ode(&source, budget(), || false).unwrap();
    let force = C64::new(0.7, -0.3);
    for omega in [0.1, 10.0, 30.0, 100.0, 1000.0] {
        let s = C64::new(0.0, -omega);
        let memory = C64::ONE + s.scale(tau);
        let gain = s.scale(tau) / memory;
        let dynamic = C64::from_re(equilibrium - mass * omega * omega)
            + s.scale(damping)
            + gain.scale(branch);
        let q = force / dynamic;
        let v = s * q;
        let z = q.scale(branch.sqrt()) / memory;
        let mismatch = q.scale(branch.sqrt()) - z;
        for (dk, dtau) in [(1.0, 0.0), (0.0, 1.0)] {
            let cross = -0.5 * dk / branch.sqrt();
            let d_q = [dk, 0.0, cross, 0.0, 0.0, 0.0, cross, 0.0, 0.0];
            let mut d_r = [0.0; 9];
            d_r[8] = -dtau / (tau * tau);
            let result = model
                .factor(omega, || false)
                .unwrap()
                .response_with_parameter_derivative(
                    &[force],
                    Direction {
                        storage: Some(&d_q),
                        dissipation: Some(&d_r),
                        ..Default::default()
                    },
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
            near(result.state_derivative[0], dq, 3e-10);
            near(result.state_derivative[2], dz, 3e-10);
            near(result.effort_derivative[2], -dmismatch, 3e-10);
            near(
                result.effort_derivative[0],
                dq.scale(equilibrium)
                    + dmismatch.scale(branch.sqrt())
                    + mismatch.scale(0.5 * dk / branch.sqrt()),
                3e-10,
            );
            let dh = 0.5
                * (equilibrium * (q.conj() * dq).re
                    + mass * (v.conj() * dv).re
                    + (mismatch.conj() * dmismatch).re);
            let heat = damping * (v.conj() * dv).re + (mismatch.conj() * dmismatch).re / tau
                - 0.5 * dtau * mismatch.norm_sq() / (tau * tau);
            near_real(result.mean_stored_energy_derivative, dh, 2e-11);
            near_real(result.dissipated_power_derivative, heat, 2e-11);
            near_real(
                result.supplied_power_derivative,
                0.5 * (force.conj() * dv).re,
                2e-11,
            );
            assert!(result.power_residual_derivative.abs() < 2e-11 * heat.abs().max(1.0));
        }
    }
}

fn coupled(j: Vec<f64>, r: Vec<f64>, q: Vec<f64>, g: Vec<f64>) -> HarmonicModel {
    let source = PortHamiltonian::new(
        4,
        2,
        j,
        r,
        g,
        Box::new(QuadraticStorage::new(q, 4).unwrap()),
    )
    .unwrap();
    HarmonicModel::from_ode(&source, budget(), || false).unwrap()
}

#[test]
fn full_mutual_storage_loss_interconnection_and_ports_match_independent_reassembly() {
    let j = vec![
        0., 0., 1., 0., 0., 0., 0., 1., -1., 0., 0., 0.2, 0., -1., -0.2, 0.,
    ];
    let mut r = vec![0.; 16];
    r[10] = 0.4;
    r[15] = 0.6;
    r[11] = 0.05;
    r[14] = 0.05;
    let mut q = vec![0.; 16];
    q[0] = 100.;
    q[5] = 225.;
    q[10] = 1.;
    q[15] = 1.;
    q[1] = -15.;
    q[4] = -15.;
    let g = vec![0.1, 0., 0., -0.2, 1., 0.3, 0.4, 1.];
    let mut dj = vec![0.; 16];
    dj[11] = 0.3;
    dj[14] = -0.3;
    let mut dr = vec![0.; 16];
    dr[11] = 0.1;
    dr[14] = 0.1;
    let mut dq = vec![0.; 16];
    dq[1] = 7.;
    dq[4] = 7.;
    dq[11] = -0.04;
    dq[14] = -0.04;
    let dg = vec![0.2, -0.1, 0., 0.1, 0.3, 0.1, -0.1, 0.2];
    let direction = Direction {
        interconnection: Some(&dj),
        dissipation: Some(&dr),
        storage: Some(&dq),
        ports: Some(&dg),
        ..Default::default()
    };
    let input = [C64::new(0.7, 0.2), C64::new(-0.3, 0.1)];
    let du = [C64::new(0.03, -0.01), C64::new(-0.02, 0.04)];
    let model = coupled(j.clone(), r.clone(), q.clone(), g.clone());
    for omega in [9.0, 10.4, 15.0] {
        let result = model
            .factor(omega, || false)
            .unwrap()
            .response_with_parameter_derivative(&input, direction, Some(&du), || false)
            .unwrap();
        let step = 1e-5;
        let perturb = |a: &[f64], d: &[f64], t: f64| -> Vec<_> {
            a.iter().zip(d).map(|(a, d)| a + t * d).collect()
        };
        let evaluate = |t| {
            let changed = coupled(
                perturb(&j, &dj, t),
                perturb(&r, &dr, t),
                perturb(&q, &dq, t),
                perturb(&g, &dg, t),
            );
            let forcing: Vec<_> = input.iter().zip(du).map(|(&u, d)| u + d.scale(t)).collect();
            changed
                .factor(omega, || false)
                .unwrap()
                .response(&forcing, || false)
                .unwrap()
        };
        let lo = evaluate(-step);
        let hi = evaluate(step);
        for (actual, a, b) in [
            (&result.state_derivative, &lo.state, &hi.state),
            (&result.effort_derivative, &lo.effort, &hi.effort),
            (&result.output_derivative, &lo.output, &hi.output),
        ] {
            for ((&value, &lo), &hi) in actual.iter().zip(a).zip(b) {
                near(value, (hi - lo).scale(0.5 / step), 2e-7);
            }
        }
        for (actual, a, b) in [
            (
                result.mean_stored_energy_derivative,
                lo.mean_stored_energy_j,
                hi.mean_stored_energy_j,
            ),
            (
                result.supplied_power_derivative,
                lo.supplied_power_w,
                hi.supplied_power_w,
            ),
            (
                result.dissipated_power_derivative,
                lo.dissipated_power_w,
                hi.dissipated_power_w,
            ),
        ] {
            near_real(actual, (b - a) * 0.5 / step, 2e-8);
        }
        for ((&actual, &a), &b) in result
            .port_power_derivative
            .iter()
            .zip(&lo.port_power_w)
            .zip(&hi.port_power_w)
        {
            near_real(actual, (b - a) * 0.5 / step, 2e-8);
        }
        assert!(result.power_residual_derivative.abs() < 1e-11);
    }
}

#[test]
fn changing_a_power_conjugate_port_retains_both_source_and_observation_derivatives() {
    for scale in [1e-120, 1.0, 1e120] {
        let source = PortHamiltonian::new(
            2,
            1,
            vec![0., 1., -1., 0.],
            vec![0., 0., 0., 0.2],
            vec![0., scale],
            Box::new(QuadraticStorage::new(vec![9., 0., 0., 1.], 2).unwrap()),
        )
        .unwrap();
        let model = HarmonicModel::from_ode(&source, budget(), || false).unwrap();
        let factor = model.factor(3.0, || false).unwrap();
        let input = C64::new(0.3 / scale, -0.2 / scale);
        let dg = [0., 0.5 * scale];
        let direction = Direction {
            ports: Some(&dg),
            ..Default::default()
        };
        let changed = factor
            .response_with_parameter_derivative(&[input], direction, None, || false)
            .unwrap();
        for (&derivative, &state) in changed.state_derivative.iter().zip(&changed.response.state) {
            near(derivative, state.scale(0.5), 2e-12);
        }
        // One factor of G drives, the other observes: dy/dtheta=y here.
        near(
            changed.output_derivative[0],
            changed.response.output[0],
            2e-12,
        );
        near_real(
            changed.mean_stored_energy_derivative,
            changed.response.mean_stored_energy_j,
            2e-12,
        );
        near_real(
            changed.dissipated_power_derivative,
            changed.response.dissipated_power_w,
            2e-12,
        );
        let unchanged = factor
            .response_with_parameter_derivative(
                &[input],
                direction,
                Some(&[input.scale(-0.5)]),
                || false,
            )
            .unwrap();
        assert!(unchanged.state_derivative.iter().all(|&v| v == C64::ZERO));
        near(
            unchanged.output_derivative[0],
            unchanged.response.output[0].scale(0.5),
            2e-12,
        );
        assert_eq!(unchanged.mean_stored_energy_derivative, 0.0);
        assert_eq!(unchanged.dissipated_power_derivative, 0.0);
        near_real(unchanged.supplied_power_derivative, 0.0, 2e-14);
    }
}

#[test]
fn nonlinear_resisting_flow_uses_the_complete_mass_chain_rule_and_incremental_work() {
    let (mass, stiffness, viscous, linear, cubic, bias_p) = (2.0, 4.0, 0.1, 0.2, 0.6, 0.6);
    let source = PortHamiltonian::new(
        2,
        2,
        vec![0., 1., -1., 0.],
        vec![0., 0., 0., viscous],
        vec![1., 0., 0., 1.],
        Box::new(QuadraticStorage::new(vec![stiffness, 0., 0., 1. / mass], 2).unwrap()),
    )
    .unwrap();
    let bias = [0.2, bias_p];
    let v0: f64 = bias_p / mass;
    let dc = [
        -v0,
        stiffness * bias[0] + (viscous + linear) * v0 + cubic * v0.powi(3),
    ];
    let hessian = |_: &[f64], d: &[f64], o: &mut [f64]| {
        o[0] = stiffness * d[0];
        o[1] = d[1] / mass;
        true
    };
    let flow = |_: &[f64], e: &[f64], o: &mut [f64]| {
        o[0] = 0.;
        o[1] = linear * e[1] + cubic * e[1].powi(3);
        true
    };
    let tangent = |_: &[f64], e: &[f64], _: &[f64], de: &[f64], o: &mut [f64]| {
        o[0] = 0.;
        o[1] = (linear + 3. * cubic * e[1] * e[1]) * de[1];
        true
    };
    let local = EquilibriumLinearization::from_dissipative_ode(
        &source,
        EquilibriumPoint {
            state: &bias,
            input: &dc,
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
    let input = [C64::new(0.3, -0.2), C64::new(0.7, 0.1)];
    let d_q = [0., 0., 0., -1. / (mass * mass)];
    let t = (linear + 3. * cubic * v0 * v0) / mass;
    let dt = -linear / (mass * mass) - 9. * cubic * bias_p * bias_p / mass.powi(4);
    let d_t = [0., 0., 0., dt];
    for omega in [0.5, 1.4, 3.0] {
        let result = local
            .model()
            .factor(omega, || false)
            .unwrap()
            .response_with_parameter_derivative(
                &input,
                Direction {
                    storage: Some(&d_q),
                    resisting_flow: Some(&d_t),
                    ..Default::default()
                },
                None,
                || false,
            )
            .unwrap();
        let s = C64::new(0., -omega);
        let d = s + C64::from_re(viscous / mass + t);
        let determinant = s * d + C64::from_re(stiffness / mass);
        let q = (d * input[0] + input[1].scale(1. / mass)) / determinant;
        let p = (s * input[1] - input[0].scale(stiffness)) / determinant;
        let rhs0 = p.scale(-1. / (mass * mass));
        let rhs1 = p.scale(viscous / (mass * mass) - dt);
        let dq = (d * rhs0 + rhs1.scale(1. / mass)) / determinant;
        let dp = (s * rhs1 - rhs0.scale(stiffness)) / determinant;
        let e = p.scale(1. / mass);
        let de = dp.scale(1. / mass) - p.scale(1. / (mass * mass));
        near(result.state_derivative[0], dq, 2e-11);
        near(result.state_derivative[1], dp, 2e-11);
        near(result.effort_derivative[1], de, 2e-11);
        let loss = viscous * (e.conj() * de).re
            + 0.5 * (de.conj() * p.scale(t) + e.conj() * (p.scale(dt) + dp.scale(t))).re;
        near_real(result.dissipated_power_derivative, loss, 2e-11);
        near_real(
            result.mean_stored_energy_derivative,
            0.5 * stiffness * (q.conj() * dq).re + 0.5 / mass * (p.conj() * dp).re
                - 0.25 / (mass * mass) * p.norm_sq(),
            2e-11,
        );
        assert!(result.power_residual_derivative.abs() < 2e-11);
    }
    near_real(
        local.equilibrium_nonlinear_dissipated_power_w(),
        linear * v0 * v0 + cubic * v0.powi(4),
        1e-14,
    );
}

#[test]
fn parameter_tangents_preserve_exact_zero_observations_and_nonzero_source_slopes() {
    let source = PortHamiltonian::new(
        4,
        1,
        vec![
            0., 1., 0., 0., -1., 0., 0., 0., 0., 0., 0., 1., 0., 0., -1., 0.,
        ],
        vec![0.; 16],
        vec![0., 1., 0., 1.],
        Box::new(
            QuadraticStorage::new(
                vec![
                    2., 0., 0., 0., 0., 1., 0., 0., 0., 0., 6., 0., 0., 0., 0., 1.,
                ],
                4,
            )
            .unwrap(),
        ),
    )
    .unwrap();
    let model = HarmonicModel::from_ode(&source, budget(), || false).unwrap();
    let factor = model.factor(2.0, || false).unwrap();
    let mut d_q = [0.; 16];
    d_q[0] = 1.;
    let direction = Direction {
        storage: Some(&d_q),
        ..Default::default()
    };
    let result = factor
        .response_with_parameter_derivative(&[C64::ONE], direction, None, || false)
        .unwrap();
    near(result.response.output[0], C64::ZERO, 2e-15);
    near(result.output_derivative[0], C64::new(0., 0.5), 2e-14);
    if result.response.output[0] == C64::ZERO {
        assert_eq!(result.output_logarithmic_derivative(0).unwrap(), None);
    }
    let zero = factor
        .response_with_parameter_derivative(&[C64::ZERO], direction, Some(&[C64::ONE]), || false)
        .unwrap();
    assert_eq!(
        zero.state_derivative,
        factor.response(&[C64::ONE], || false).unwrap().state
    );
    assert_eq!(zero.mean_stored_energy_derivative, 0.0);
    assert_eq!(zero.supplied_power_derivative, 0.0);
    assert_eq!(zero.output_logarithmic_derivative(0).unwrap(), None);
}

#[test]
fn parameter_admission_and_cancellation_preserve_the_reusable_factor() {
    let source = fs_phs::mass_spring_damper(1., 4., 0.2).unwrap();
    let model = HarmonicModel::from_ode(&source, budget(), || false).unwrap();
    let factor = model.factor(2., || false).unwrap();
    let d_q = [1., 0., 0., 0.];
    let direction = Direction {
        storage: Some(&d_q),
        ..Default::default()
    };
    let before = factor.response(&[C64::ONE], || false).unwrap();
    let mut polls = 0;
    let expected = factor
        .response_with_parameter_derivative(&[C64::ONE], direction, None, || {
            polls += 1;
            false
        })
        .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut count = 0;
        assert!(matches!(
            factor.response_with_parameter_derivative(&[C64::ONE], direction, None, || {
                count += 1;
                count >= stop
            }),
            Err(HarmonicError::Cancelled)
        ));
    }
    for bad in [
        Direction {
            interconnection: Some(&[1., 0., 0., 0.]),
            ..Default::default()
        },
        Direction {
            storage: Some(&[0., 1., 0., 0.]),
            ..Default::default()
        },
        Direction {
            dissipation: Some(&[0.]),
            ..Default::default()
        },
        Direction {
            ports: Some(&[0., f64::NAN]),
            ..Default::default()
        },
        Direction {
            resisting_flow: Some(&[0., 0., 0., f64::INFINITY]),
            ..Default::default()
        },
    ] {
        assert!(
            factor
                .response_with_parameter_derivative(&[C64::ONE], bad, None, || false)
                .is_err()
        );
    }
    for bad in [vec![], vec![C64::new(0., f64::INFINITY)]] {
        assert!(
            factor
                .response_with_parameter_derivative(&[C64::ONE], direction, Some(&bad), || false)
                .is_err()
        );
    }
    assert!(expected.output_logarithmic_derivative(1).is_err());
    assert_eq!(before, factor.response(&[C64::ONE], || false).unwrap());
    assert_eq!(
        expected,
        factor
            .response_with_parameter_derivative(&[C64::ONE], direction, None, || false)
            .unwrap()
    );
    for tight in [
        HarmonicBudget {
            max_dense_work: 100,
            ..budget()
        },
        HarmonicBudget {
            max_dense_bytes: 1100,
            ..budget()
        },
    ] {
        let limited = HarmonicModel::from_ode(&source, tight, || false).unwrap();
        let factor = limited.factor(2., || false).unwrap();
        factor.response(&[C64::ONE], || false).unwrap();
        assert!(matches!(
            factor.response_with_parameter_derivative(&[C64::ONE], direction, None, || false),
            Err(HarmonicError::Invalid { .. })
        ));
    }
}
