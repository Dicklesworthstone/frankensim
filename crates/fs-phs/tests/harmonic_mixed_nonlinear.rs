//! Independent original-coordinate jets of a nonlinear biased resistance.
//! This is a local incremental response, not the biased source's total heat.

use fs_math::c64::C64;
use fs_phs::harmonic::equilibrium::{
    EquilibriumLinearization, EquilibriumPoint, LinearizationBudget, NonlinearDissipation,
};
use fs_phs::harmonic::{HarmonicBudget, HarmonicParameterDirection};
use fs_phs::{PortHamiltonian, QuadraticStorage};

// Value, omega, omega-omega, mass, omega-mass. These are scalar closed-form
// numerator/denominator derivatives, independent of native state solves.
type Jet = [C64; 5];

fn quotient(n: Jet, d: Jet) -> Jet {
    let value = n[0] / d[0];
    let omega = (n[1] - value * d[1]) / d[0];
    let second = (n[2] - omega * d[1].scale(2.0) - value * d[2]) / d[0];
    let parameter = (n[3] - value * d[3]) / d[0];
    let mixed = (n[4] - omega * d[3] - parameter * d[1] - value * d[4]) / d[0];
    [value, omega, second, parameter, mixed]
}

fn norm_square(z: Jet) -> [f64; 5] {
    [
        z[0].norm_sq(),
        2.0 * (z[0].conj() * z[1]).re,
        2.0 * ((z[0].conj() * z[2]).re + z[1].norm_sq()),
        2.0 * (z[0].conj() * z[3]).re,
        2.0 * ((z[3].conj() * z[1]).re + (z[0].conj() * z[4]).re),
    ]
}

fn near(actual: C64, expected: C64, field: &str, omega: f64, order: usize) {
    assert!(
        (actual - expected).abs() < 4e-11 * expected.abs().max(1.0),
        "{field}, omega={omega}, order={order}: {actual:?} != {expected:?}"
    );
}

fn real_near(actual: f64, expected: f64, field: &str, omega: f64, order: usize) {
    near(
        C64::from_re(actual),
        C64::from_re(expected),
        field,
        omega,
        order,
    );
}

#[test]
fn total_nonlinear_resisting_flow_mass_jet_matches_original_rational_fields_and_signed_work() {
    let mass: f64 = 2.0;
    let stiffness: f64 = 4.0;
    let linear: f64 = 0.2;
    let budget = HarmonicBudget {
        max_dimension: 8,
        max_dense_bytes: 100_000,
        max_dense_work: 100_000,
        relative_residual_tolerance: 1e-12,
    };
    let input = [C64::new(0.3, -0.2), C64::new(0.7, 0.1)];
    // The second local law still dissipates at its positive-velocity bias,
    // while its differential resistance is negative. Original viscous loss
    // keeps the complete incremental oscillator dissipative. Do not clip T
    // or mislabel incremental work as total physical biased-source heat.
    for (cubic, viscous, bias_p) in [(0.6_f64, 0.1_f64, 0.6_f64), (-0.6, 0.6, 1.0)] {
        let source = PortHamiltonian::new(
            2,
            2,
            vec![0.0, 1.0, -1.0, 0.0],
            vec![0.0, 0.0, 0.0, viscous],
            vec![1.0, 0.0, 0.0, 1.0],
            Box::new(QuadraticStorage::new(vec![stiffness, 0.0, 0.0, 1.0 / mass], 2).unwrap()),
        )
        .unwrap();
        let bias = [0.2, bias_p];
        let v0 = bias_p / mass;
        let dc = [
            -v0,
            stiffness * bias[0] + (viscous + linear) * v0 + cubic * v0.powi(3),
        ];
        let hessian = |_: &[f64], direction: &[f64], out: &mut [f64]| {
            out[0] = stiffness * direction[0];
            out[1] = direction[1] / mass;
            true
        };
        let flow = |_: &[f64], effort: &[f64], out: &mut [f64]| {
            out[0] = 0.0;
            out[1] = linear * effort[1] + cubic * effort[1].powi(3);
            true
        };
        let tangent = |_: &[f64], effort: &[f64], _: &[f64], de: &[f64], out: &mut [f64]| {
            out[0] = 0.0;
            out[1] = (linear + 3.0 * cubic * effort[1].powi(2)) * de[1];
            true
        };
        let local = EquilibriumLinearization::from_dissipative_ode(
            &source,
            EquilibriumPoint {
                state: &bias,
                input: &dc,
                balance_scales: &[1.0, 1.0],
                relative_tolerance: 1e-12,
            },
            LinearizationBudget {
                harmonic: budget,
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

        // Mass varies at fixed original bias momentum. T=D_e*Q therefore
        // contains both the inverse-mass map and the changing bias velocity.
        // The authored total dT includes all of that chain rule.
        let t = (linear + 3.0 * cubic * v0 * v0) / mass;
        let dt = -linear / mass.powi(2) - 9.0 * cubic * bias_p * bias_p / mass.powi(4);
        let dq = [0.0, 0.0, 0.0, -1.0 / mass.powi(2)];
        let d_t = [0.0, 0.0, 0.0, dt];
        let alpha = viscous / mass + t;
        let alpha_mass = -viscous / mass.powi(2) + dt;
        let loss_coefficient = viscous + linear + 3.0 * cubic * v0 * v0;
        let loss_coefficient_mass = -6.0 * cubic * bias_p * bias_p / mass.powi(3);
        for omega in [0.2, 1.0, 2.0_f64.sqrt(), 2.8] {
            let result = local
                .model()
                .factor(omega, || false)
                .unwrap()
                .response_with_frequency_parameter_derivative(
                    &input,
                    HarmonicParameterDirection {
                        storage: Some(&dq),
                        resisting_flow: Some(&d_t),
                        ..Default::default()
                    },
                    None,
                    None,
                    None,
                    None,
                    || false,
                )
                .unwrap();
            // Original equations: s*q=p/m+u0 and
            // s*p=-k*q-(c/m+T)*p+u1. No native matrix/factor/jet oracle.
            let s = C64::new(0.0, -omega);
            let denominator = [
                s * (s + C64::from_re(alpha)) + C64::from_re(stiffness / mass),
                C64::new(-2.0 * omega, -alpha),
                C64::from_re(-2.0),
                s.scale(alpha_mass) - C64::from_re(stiffness / mass.powi(2)),
                C64::new(0.0, -alpha_mass),
            ];
            let q = quotient(
                [
                    (s + C64::from_re(alpha)) * input[0] + input[1].scale(1.0 / mass),
                    C64::new(0.0, -1.0) * input[0],
                    C64::ZERO,
                    input[0].scale(alpha_mass) - input[1].scale(1.0 / mass.powi(2)),
                    C64::ZERO,
                ],
                denominator,
            );
            let p = quotient(
                [
                    s * input[1] - input[0].scale(stiffness),
                    C64::new(0.0, -1.0) * input[1],
                    C64::ZERO,
                    C64::ZERO,
                    C64::ZERO,
                ],
                denominator,
            );
            let velocity = [
                p[0].scale(1.0 / mass),
                p[1].scale(1.0 / mass),
                p[2].scale(1.0 / mass),
                p[3].scale(1.0 / mass) - p[0].scale(1.0 / mass.powi(2)),
                p[4].scale(1.0 / mass) - p[1].scale(1.0 / mass.powi(2)),
            ];
            let second = &result.frequency_response;
            let first = &second.frequency_response;
            let base = &first.response;
            let parameter = &result.parameter_response;
            let states = [
                &base.state,
                &first.state_derivative,
                &second.state_second_derivative,
                &parameter.state_derivative,
                &result.state_mixed_derivative,
            ];
            let efforts = [
                &base.effort,
                &first.effort_derivative,
                &second.effort_second_derivative,
                &parameter.effort_derivative,
                &result.effort_mixed_derivative,
            ];
            let outputs = [
                &base.output,
                &first.output_derivative,
                &second.output_second_derivative,
                &parameter.output_derivative,
                &result.output_mixed_derivative,
            ];
            let energy = [
                base.mean_stored_energy_j,
                first.mean_stored_energy_derivative,
                second.mean_stored_energy_second_derivative,
                parameter.mean_stored_energy_derivative,
                result.mean_stored_energy_mixed_derivative,
            ];
            let supply = [
                base.supplied_power_w,
                first.supplied_power_derivative,
                second.supplied_power_second_derivative,
                parameter.supplied_power_derivative,
                result.supplied_power_mixed_derivative,
            ];
            let loss = [
                base.dissipated_power_w,
                first.dissipated_power_derivative,
                second.dissipated_power_second_derivative,
                parameter.dissipated_power_derivative,
                result.dissipated_power_mixed_derivative,
            ];
            let q_square = norm_square(q);
            let p_square = norm_square(p);
            let v_square = norm_square(velocity);
            for order in 0..5 {
                near(states[order][0], q[order], "q", omega, order);
                near(states[order][1], p[order], "p", omega, order);
                let expected_effort = [q[order].scale(stiffness), velocity[order]];
                for port in 0..2 {
                    near(
                        efforts[order][port],
                        expected_effort[port],
                        "effort",
                        omega,
                        order,
                    );
                    near(
                        outputs[order][port],
                        expected_effort[port],
                        "output",
                        omega,
                        order,
                    );
                    let port_power = [
                        base.port_power_w[port],
                        first.port_power_derivative[port],
                        second.port_power_second_derivative[port],
                        parameter.port_power_derivative[port],
                        result.port_power_mixed_derivative[port],
                    ];
                    real_near(
                        port_power[order],
                        0.5 * (input[port].conj() * expected_effort[port]).re,
                        "original port work",
                        omega,
                        order,
                    );
                }
                let mut expected_energy =
                    0.25 * (stiffness * q_square[order] + p_square[order] / mass);
                let mut expected_loss = 0.5 * loss_coefficient * v_square[order];
                if order >= 3 {
                    let lower = if order == 3 { 0 } else { 1 };
                    expected_energy -= 0.25 * p_square[lower] / mass.powi(2);
                    expected_loss += 0.5 * loss_coefficient_mass * v_square[lower];
                }
                let expected_supply = 0.5
                    * (input[0].conj() * expected_effort[0] + input[1].conj() * expected_effort[1])
                        .re;
                real_near(
                    energy[order],
                    expected_energy,
                    "original storage",
                    omega,
                    order,
                );
                real_near(
                    loss[order],
                    expected_loss,
                    "independent R/T work",
                    omega,
                    order,
                );
                real_near(supply[order], expected_supply, "source work", omega, order);
                real_near(
                    supply[order],
                    loss[order],
                    "incremental balance",
                    omega,
                    order,
                );
            }
            if cubic < 0.0 {
                assert!(t < 0.0);
                assert!(base.dissipated_power_w < 0.5 * viscous * v_square[0]);
            }
            assert!(result.relative_mixed_derivative_residual < 1e-12);
            assert!(second.relative_second_derivative_residual < 1e-12);
        }
        real_near(
            local.equilibrium_nonlinear_dissipated_power_w(),
            linear * v0 * v0 + cubic * v0.powi(4),
            "separate physical bias heat",
            0.0,
            0,
        );
        assert!(local.equilibrium_nonlinear_dissipated_power_w() > 0.0);
    }
}
