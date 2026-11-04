//! Loaded derivatives of actual biased nonlinear sources, in fixed coordinates.
use fs_math::c64::C64;
use fs_phs::harmonic::equilibrium::{
    EquilibriumLinearization, EquilibriumPoint, LinearizationBudget, NonlinearDissipation,
};
use fs_phs::harmonic::scattering::{ScatteringBudget, ScatteringModel};
use fs_phs::harmonic::{HarmonicBudget, HarmonicParameterDirection};
use fs_phs::{PortHamiltonian, Storage};

const STIFFNESS: f64 = 4.0;
const VISCOUS: f64 = 0.1;
const LINEAR: f64 = 0.2;
const POSITION_DRAG: f64 = 0.7;
const BIAS: [f64; 2] = [0.2, 0.6];
const REFERENCE: f64 = 0.4;

#[derive(Clone, Copy)]
struct Coefficients {
    mass: f64,
    hardening: f64,
    cubic_drag: f64,
}

struct Spring {
    mass: f64,
    hardening: f64,
}
impl Storage for Spring {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        0.5 * STIFFNESS * x[0].powi(2)
            + 0.25 * self.hardening * x[0].powi(4)
            + 0.5 * x[1].powi(2) / self.mass
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        out[0] = STIFFNESS * x[0] + self.hardening * x[0].powi(3);
        out[1] = x[1] / self.mass;
    }
}

fn budget() -> LinearizationBudget {
    LinearizationBudget {
        harmonic: HarmonicBudget {
            max_dimension: 8,
            max_dense_bytes: 2_000_000,
            max_dense_work: 4_000_000,
            relative_residual_tolerance: 1e-11,
        },
        max_hessian_actions: 2,
    }
}
fn waves() -> ScatteringBudget {
    ScatteringBudget {
        maximum_ports: 2,
        maximum_payload_bytes: 1_000_000,
        maximum_work: 2_000_000,
    }
}
fn close(a: C64, b: C64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * a.abs().max(b.abs()).max(1e-12),
        "{a:?} != {b:?}"
    );
}
fn scalar(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * a.abs().max(b.abs()).max(1e-16),
        "{a:e} != {b:e}"
    );
}

fn local(c: Coefficients) -> EquilibriumLinearization {
    let source = PortHamiltonian::new(
        2,
        2,
        vec![0.0, 1.0, -1.0, 0.0],
        vec![0.0, 0.0, 0.0, VISCOUS],
        vec![1.0, 0.0, 0.0, 1.0],
        Box::new(Spring {
            mass: c.mass,
            hardening: c.hardening,
        }),
    )
    .unwrap();
    let v0 = BIAS[1] / c.mass;
    let factor = 1.0 + POSITION_DRAG * BIAS[0].powi(2);
    let dc = [
        -v0,
        STIFFNESS * BIAS[0]
            + c.hardening * BIAS[0].powi(3)
            + VISCOUS * v0
            + factor * (LINEAR * v0 + c.cubic_drag * v0.powi(3)),
    ];
    let hessian = |x: &[f64], dx: &[f64], out: &mut [f64]| {
        out[0] = (STIFFNESS + 3.0 * c.hardening * x[0].powi(2)) * dx[0];
        out[1] = dx[1] / c.mass;
        true
    };
    // Positive original heat for every state, with a nonsymmetric total tangent.
    let flow = |x: &[f64], e: &[f64], out: &mut [f64]| {
        out[0] = 0.0;
        out[1] =
            (1.0 + POSITION_DRAG * x[0].powi(2)) * (LINEAR * e[1] + c.cubic_drag * e[1].powi(3));
        true
    };
    let tangent = |x: &[f64], e: &[f64], dx: &[f64], de: &[f64], out: &mut [f64]| {
        out[0] = 0.0;
        out[1] = 2.0 * POSITION_DRAG * x[0] * (LINEAR * e[1] + c.cubic_drag * e[1].powi(3)) * dx[0]
            + (1.0 + POSITION_DRAG * x[0].powi(2))
                * (LINEAR + 3.0 * c.cubic_drag * e[1].powi(2))
                * de[1];
        true
    };
    EquilibriumLinearization::from_dissipative_ode(
        &source,
        EquilibriumPoint {
            state: &BIAS,
            input: &dc,
            balance_scales: &[1.0; 2],
            relative_tolerance: 1e-12,
        },
        budget(),
        &hessian,
        NonlinearDissipation {
            flow: &flow,
            tangent: &tangent,
            max_tangent_actions: 2,
        },
        || false,
    )
    .unwrap()
}

fn loaded(c: Coefficients) -> ScatteringModel {
    ScatteringModel::new(local(c).into_model(), &[(1, REFERENCE)], waves(), || false).unwrap()
}

#[test]
fn loaded_nonlinear_mass_and_constitutive_tangents_match_original_source_reassembly() {
    let c = Coefficients {
        mass: 2.0,
        hardening: 3.0,
        cubic_drag: 0.6,
    };
    let input = [C64::new(0.003, -0.002), C64::new(0.008, 0.004)];
    let local = local(c);
    let bias_heat = local.equilibrium_nonlinear_dissipated_power_w();
    let v0 = BIAS[1] / c.mass;
    let multiplier = 1.0 + POSITION_DRAG * BIAS[0].powi(2);
    scalar(
        bias_heat,
        multiplier * (LINEAR * v0.powi(2) + c.cubic_drag * v0.powi(4)),
        1e-13,
    );
    // Consuming the admitted model preserves its exact local physics, without
    // rebuilding from the nonlinear source or losing the nonsymmetric T.
    let bare = local
        .model()
        .factor(0.7, || false)
        .unwrap()
        .response(&input, || false)
        .unwrap();
    let original_t = local
        .model()
        .nonlinear_resisting_tangent()
        .unwrap()
        .to_vec();
    let owned = local.into_model();
    assert_eq!(owned.nonlinear_resisting_tangent().unwrap(), original_t);
    assert_eq!(
        owned
            .factor(0.7, || false)
            .unwrap()
            .response(&input, || false)
            .unwrap(),
        bare
    );
    let model = ScatteringModel::new(owned, &[(1, REFERENCE)], waves(), || false).unwrap();
    let (tq, tp) = (original_t[2], original_t[3]);
    assert!(tq > 0.0 && original_t[1] == 0.0);
    for (dm, dh, db) in [(1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (0.4, 0.3, -0.2)] {
        let dq = [3.0 * BIAS[0].powi(2) * dh, 0.0, 0.0, -dm / c.mass.powi(2)];
        // TOTAL dT includes the change of biased effort v0=p0/m and of Q_pp.
        let dtq = -2.0 * POSITION_DRAG * BIAS[0] * (LINEAR + 3.0 * c.cubic_drag * v0.powi(2)) * v0
            / c.mass
            * dm
            + 2.0 * POSITION_DRAG * BIAS[0] * v0.powi(3) * db;
        let dtp = multiplier
            * (-LINEAR / c.mass.powi(2) - 9.0 * c.cubic_drag * BIAS[1].powi(2) / c.mass.powi(4))
            * dm
            + multiplier * 3.0 * v0.powi(2) / c.mass * db;
        let dt = [0.0, 0.0, dtq, dtp];
        for omega in [0.5, 1.48, 3.0] {
            let factor = model.factor(omega, || false).unwrap();
            let result = factor
                .response_with_parameter_derivative(
                    &input,
                    HarmonicParameterDirection {
                        storage: Some(&dq),
                        resisting_flow: Some(&dt),
                        ..Default::default()
                    },
                    None,
                    None,
                    || false,
                )
                .unwrap();
            assert_eq!(result.response, factor.response(&input, || false).unwrap());
            let step = 1e-5;
            let evaluate = |amount: f64| {
                loaded(Coefficients {
                    mass: c.mass + amount * dm,
                    hardening: c.hardening + amount * dh,
                    cubic_drag: c.cubic_drag + amount * db,
                })
                .factor(omega, || false)
                .unwrap()
                .response(&input, || false)
                .unwrap()
            };
            let (lo, hi) = (evaluate(-step), evaluate(step));
            for (actual, low, high) in [
                (
                    &result.state_derivative,
                    &lo.physical.state,
                    &hi.physical.state,
                ),
                (
                    &result.effort_derivative,
                    &lo.physical.effort,
                    &hi.physical.effort,
                ),
                (
                    &result.output_derivative,
                    &lo.physical.output,
                    &hi.physical.output,
                ),
                (
                    &result.physical_input_derivative,
                    &lo.physical_inputs,
                    &hi.physical_inputs,
                ),
            ] {
                for ((&value, &a), &b) in actual.iter().zip(low).zip(high) {
                    close(value, (b - a).scale(0.5 / step), 3e-7);
                }
            }
            for (value, a, b) in [
                (
                    result.mean_stored_energy_derivative,
                    lo.physical.mean_stored_energy_j,
                    hi.physical.mean_stored_energy_j,
                ),
                (
                    result.dissipated_power_derivative,
                    lo.physical.dissipated_power_w,
                    hi.physical.dissipated_power_w,
                ),
                (
                    result.supplied_power_derivative,
                    lo.physical.supplied_power_w,
                    hi.physical.supplied_power_w,
                ),
                (
                    result.wave_derivatives[0].outgoing_power_w,
                    lo.waves[0].outgoing_power_w,
                    hi.waves[0].outgoing_power_w,
                ),
            ] {
                scalar(value, (b - a) * 0.5 / step, 3e-7);
            }
            let response = &result.response.physical;
            let (q, p, v) = (response.state[0], response.state[1], response.effort[1]);
            let (q_prime, p_prime, v_prime) = (
                result.state_derivative[0],
                result.state_derivative[1],
                result.effort_derivative[1],
            );
            let resisting = q.scale(tq) + p.scale(tp);
            let d_resisting = q.scale(dtq) + q_prime.scale(tq) + p.scale(dtp) + p_prime.scale(tp);
            let independent_loss = VISCOUS * (v.conj() * v_prime).re
                + 0.5 * (v_prime.conj() * resisting + v.conj() * d_resisting).re;
            scalar(result.dissipated_power_derivative, independent_loss, 2e-11);
            let wave = &result.response.waves[0];
            let independent_supply = -(wave.outgoing.conj() * result.output_derivative[1]).re
                / REFERENCE
                + 0.5 * (input[0].conj() * result.output_derivative[0]).re;
            scalar(result.supplied_power_derivative, independent_supply, 2e-11);
            let scale = independent_loss.abs() + independent_supply.abs() + 1e-20;
            assert!(result.power_residual_derivative.abs() <= 1e-9 * scale);
            assert!(result.wave_power_residual_derivative.abs() <= 1e-9 * scale);
            assert!(result.relative_derivative_residual <= 1e-11);
            assert!(response.dissipated_power_w.abs() < 0.01 * bias_heat);
        }
    }
}

fn saturating(gain: f64) -> EquilibriumLinearization {
    let (mass, viscous, velocity) = (2.0, 0.01, 3.0_f64);
    let source = PortHamiltonian::new(
        2,
        2,
        vec![0.0, 1.0, -1.0, 0.0],
        vec![0.0, 0.0, 0.0, viscous],
        vec![1.0, 0.0, 0.0, 1.0],
        Box::new(Spring {
            mass,
            hardening: 0.0,
        }),
    )
    .unwrap();
    let state = [0.2, mass * velocity];
    let input = [
        -velocity,
        STIFFNESS * state[0] + viscous * velocity + gain * velocity / (1.0 + velocity.powi(2)),
    ];
    let hessian = |_: &[f64], dx: &[f64], out: &mut [f64]| {
        out[0] = STIFFNESS * dx[0];
        out[1] = dx[1] / mass;
        true
    };
    let flow = |_: &[f64], e: &[f64], out: &mut [f64]| {
        out[0] = 0.0;
        out[1] = gain * e[1] / (1.0 + e[1].powi(2));
        true
    };
    let tangent = |_: &[f64], e: &[f64], _: &[f64], de: &[f64], out: &mut [f64]| {
        out[0] = 0.0;
        out[1] = gain * (1.0 - e[1].powi(2)) / (1.0 + e[1].powi(2)).powi(2) * de[1];
        true
    };
    EquilibriumLinearization::from_dissipative_ode(
        &source,
        EquilibriumPoint {
            state: &state,
            input: &input,
            balance_scales: &[4.0; 2],
            relative_tolerance: 1e-12,
        },
        budget(),
        &hessian,
        NonlinearDissipation {
            flow: &flow,
            tangent: &tangent,
            max_tangent_actions: 2,
        },
        || false,
    )
    .unwrap()
}

#[test]
fn passive_biased_law_retains_negative_incremental_work_and_its_loaded_derivative() {
    let local = saturating(1.0);
    scalar(local.equilibrium_nonlinear_dissipated_power_w(), 0.9, 1e-13);
    let reference = 0.5;
    let model =
        ScatteringModel::new(local.into_model(), &[(1, reference)], waves(), || false).unwrap();
    let incident = C64::new(0.01, 0.004);
    let input = [C64::ZERO, incident];
    let omega = 2.0_f64.sqrt();
    let dt = [0.0, 0.0, 0.0, -0.04];
    let result = model
        .factor(omega, || false)
        .unwrap()
        .response_with_parameter_derivative(
            &input,
            HarmonicParameterDirection {
                resisting_flow: Some(&dt),
                ..Default::default()
            },
            None,
            None,
            || false,
        )
        .unwrap();
    // At mechanical resonance, c_incremental=0.01-0.08. The positive wave
    // termination stabilizes the response but must not be added to material heat.
    let damping = -0.07;
    let denominator = 1.0 + reference * damping;
    let velocity = incident.scale(2.0 / denominator);
    let derivative = incident.scale(2.0 * reference * 0.08 / denominator.powi(2));
    close(result.response.physical.output[1], velocity, 1e-11);
    close(result.output_derivative[1], derivative, 1e-11);
    scalar(
        result.response.physical.dissipated_power_w,
        0.5 * damping * velocity.norm_sq(),
        1e-11,
    );
    let loss_derivative = -0.04 * velocity.norm_sq() + damping * (velocity.conj() * derivative).re;
    scalar(result.dissipated_power_derivative, loss_derivative, 1e-11);
    assert!(result.response.physical.dissipated_power_w < 0.0);
    assert!(result.dissipated_power_derivative < 0.0);
    let wave = &result.response.waves[0];
    assert!(wave.outgoing_power_w > wave.incident_power_w);
    scalar(
        result.supplied_power_derivative,
        -(wave.outgoing.conj() * derivative).re / reference,
        1e-11,
    );
    let step = 1e-5;
    let evaluate = |gain| {
        ScatteringModel::new(
            saturating(gain).into_model(),
            &[(1, reference)],
            waves(),
            || false,
        )
        .unwrap()
        .factor(omega, || false)
        .unwrap()
        .response(&input, || false)
        .unwrap()
    };
    let (lo, hi) = (evaluate(1.0 - step), evaluate(1.0 + step));
    close(
        result.output_derivative[1],
        (hi.physical.output[1] - lo.physical.output[1]).scale(0.5 / step),
        3e-7,
    );
    scalar(
        result.dissipated_power_derivative,
        (hi.physical.dissipated_power_w - lo.physical.dissipated_power_w) * 0.5 / step,
        3e-7,
    );
    assert!(result.relative_derivative_residual <= 1e-11);
    assert!(result.wave_power_residual_derivative.abs() <= 1e-10 * loss_derivative.abs());
}
