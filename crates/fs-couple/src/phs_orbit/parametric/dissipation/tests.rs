//! G0/G1: original nonlinear losses, pumped memory and complete response cycles.

use super::*;
use fs_orbit::sensitivity::{ForcedSensitivityControls, ForcedSensitivityError};
use fs_phs::harmonic::equilibrium::LinearizationError;

const BETA: f64 = 0.35;
const SIGMA: f64 = 0.2;
const ETA: f64 = 0.1;

fn flow(x: &[f64], effort: &[f64], out: &mut [f64]) -> bool {
    let w = effort[1] + ETA * effort[2];
    let d = BETA * (1.0 + SIGMA * x[0] * x[0]) * w.powi(3);
    out.copy_from_slice(&[0.0, d, ETA * d]);
    true
}
fn tangent(x: &[f64], e: &[f64], dx: &[f64], de: &[f64], out: &mut [f64]) -> bool {
    let w = e[1] + ETA * e[2];
    let d = BETA
        * (2.0 * SIGMA * x[0] * dx[0] * w.powi(3)
            + 3.0 * (1.0 + SIGMA * x[0] * x[0]) * w * w * (de[1] + ETA * de[2]));
    out.copy_from_slice(&[0.0, d, ETA * d]);
    true
}
fn law() -> NonlinearDissipation<'static> {
    NonlinearDissipation {
        flow: &flow,
        tangent: &tangent,
        max_tangent_actions: 3,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        a.is_finite() && (a - b).abs() <= tolerance * (1.0 + b.abs()),
        "{a:e} != {b:e}"
    );
}
fn parameters() -> ForcedSensitivityControls {
    ForcedSensitivityControls {
        samples_per_period: 128,
        maximum_unknowns: 256,
        maximum_payload_bytes: 16 << 20,
        maximum_work: 500_000_000,
        maximum_base_residual: 1e-9,
        maximum_tangent_sampling_defect: 1e-9,
        maximum_condition_estimate: 1e10,
        maximum_relative_residual: 1e-11,
        maximum_derivative_residual: 1e-9,
    }
}

#[test]
fn pumped_effort_and_complete_loss_chain_preserve_original_memory_work() {
    let omega = 0.8;
    let system = memory_system();
    let input = drive(omega, 0.2, -0.2);
    let mut matrices = vec![vec![C64::ZERO; 9]; 3];
    matrices[0][0] = C64::from_re(0.07);
    matrices[1][2] = complex(0.06, 0.3);
    matrices[1][6] = matrices[1][2];
    matrices[2][0] = complex(0.125, 0.4);
    let pump = PeriodicQuadraticModulation::new(omega, 3, matrices).unwrap();
    let problem = ParametricHamiltonian::new(&system, &input, &pump)
        .unwrap()
        .with_dissipation(law())
        .unwrap();
    let third = problem.over_periods(3).unwrap();
    let sixth = third.over_periods(2).unwrap();
    assert_eq!(sixth.period_multiplier(), 6);
    for (time, x) in [(0.73, [0.31, -0.27, 0.08]), (1.31, [-0.24, 0.43, 0.51])] {
        let cross = 0.12 * (omega * time + 0.3).cos();
        let stiffness = 0.07 + 0.25 * (2.0 * omega * time + 0.4).cos();
        let e = [
            (2.6 + stiffness) * x[0] + 0.6 * x[0].powi(3) + (-0.8 + cross) * x[2],
            x[1] / 1.2,
            0.8 * x[2] + (-0.8 + cross) * x[0],
        ];
        let w = e[1] + ETA * e[2];
        let d = BETA * (1.0 + SIGMA * x[0] * x[0]) * w.powi(3);
        let expected = [
            e[1],
            -e[0] - 0.32 * e[1] - d + 0.2 * (omega * time - 0.2).cos(),
            -0.7 * e[2] - ETA * d,
        ];
        let mut actual = [0.0; 3];
        problem.vector_field_into(time, &x, &mut actual).unwrap();
        for j in 0..3 {
            close(actual[j], expected[j], 2e-14);
        }
        assert!((actual[2] + 0.7 * e[2]).abs() > 1e-7);
        let h = [
            2.6 + stiffness + 1.8 * x[0] * x[0],
            0.0,
            -0.8 + cross,
            0.0,
            1.0 / 1.2,
            0.0,
            -0.8 + cross,
            0.0,
            0.8,
        ];
        let mut jacobian = [0.0; 9];
        problem
            .state_jacobian_into(time, &x, &mut jacobian)
            .unwrap();
        for j in 0..3 {
            let dd = BETA
                * (if j == 0 {
                    2.0 * SIGMA * x[0] * w.powi(3)
                } else {
                    0.0
                } + 3.0 * (1.0 + SIGMA * x[0] * x[0]) * w * w * (h[3 + j] + ETA * h[6 + j]));
            close(jacobian[j], h[3 + j], 3e-14);
            close(jacobian[3 + j], -h[j] - 0.32 * h[3 + j] - dd, 3e-14);
            close(jacobian[6 + j], -0.7 * h[6 + j] - ETA * dd, 3e-14);
        }
        let power = problem.power(time, &x).unwrap();
        close(
            power.nonlinear_power,
            BETA * (1.0 + SIGMA * x[0] * x[0]) * w.powi(4),
            2e-14,
        );
        close(
            power.energy_rate,
            power.pump_power + power.port_power - power.resistance_power - power.nonlinear_power,
            3e-14,
        );
        let mut other = [0.0; 3];
        let mut other_jacobian = [0.0; 9];
        sixth.vector_field_into(time, &x, &mut other).unwrap();
        sixth
            .state_jacobian_into(time, &x, &mut other_jacobian)
            .unwrap();
        assert_eq!(actual, other);
        assert_eq!(jacobian, other_jacobian);
        let other_power = sixth.power(time, &x).unwrap();
        assert_eq!(power.nonlinear_power, other_power.nonlinear_power);
        assert_eq!(power.pump_power, other_power.pump_power);
        assert_eq!(power.energy_rate, other_power.energy_rate);
    }
}

const OMEGA: f64 = 2.1;
const MEMORY: f64 = 0.08;
const LOSS: f64 = 0.06;
const RELAXATION: f64 = 0.4;
const BIAS: f64 = 0.03;
const PUMP: f64 = 0.6;
const ORDER: usize = 12;
struct PeriodStorage;
impl Storage for PeriodStorage {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        0.5 * (x[0] * x[0] + x[1] * x[1] + MEMORY * (x[2] - x[0]).powi(2)) + 0.25 * x[0].powi(4)
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        out[0] = x[0] + x[0].powi(3) + MEMORY * (x[0] - x[2]);
        out[1] = x[1];
        out[2] = MEMORY * (x[2] - x[0]);
    }
    fn continuous_tangent(&self) -> Option<&dyn ContinuousStorageTangent> {
        Some(self)
    }
}
impl ContinuousStorageTangent for PeriodStorage {
    fn apply(&self, x: &[f64], dx: &[f64], out: &mut [f64]) -> bool {
        out[0] = (1.0 + 3.0 * x[0] * x[0] + MEMORY) * dx[0] - MEMORY * dx[2];
        out[1] = dx[1];
        out[2] = MEMORY * (dx[2] - dx[0]);
        true
    }
}
fn period_components() -> (
    PortHamiltonian,
    PeriodicPortInput,
    PeriodicQuadraticModulation,
) {
    let system = PortHamiltonian::new(
        3,
        1,
        vec![0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        vec![0.0, 0.0, 0.0, 0.0, LOSS, 0.0, 0.0, 0.0, RELAXATION / MEMORY],
        vec![0.0, 1.0, 0.0],
        Box::new(PeriodStorage),
    )
    .unwrap();
    let input = PeriodicPortInput::new(OMEGA, vec![vec![C64::from_re(BIAS)]]).unwrap();
    let mut q = vec![vec![C64::ZERO; 9]; 2];
    q[1][0] = complex(0.5 * PUMP, 0.2);
    (
        system,
        input,
        PeriodicQuadraticModulation::new(OMEGA, 3, q).unwrap(),
    )
}
fn original_rhs(time: f64, x: &[f64; 3]) -> [f64; 3] {
    let w = x[1] + ETA * MEMORY * (x[2] - x[0]);
    let d = BETA * (1.0 + SIGMA * x[0] * x[0]) * w.powi(3);
    [
        x[1],
        -(1.0 + PUMP * (OMEGA * time + 0.2).cos()) * x[0]
            - x[0].powi(3)
            - MEMORY * (x[0] - x[2])
            - LOSS * x[1]
            - d
            + BIAS,
        RELAXATION * (x[0] - x[2]) - ETA * d,
    ]
}
fn original_cycle(steps: usize) -> (HbOrbit, [f64; 5]) {
    let response = OMEGA / 2.0;
    let dt = TAU / (response * steps as f64);
    let mut state = [0.2, 0.1, 0.0];
    for sample in 0..100 * steps {
        state = rk4(&original_rhs, sample as f64 * dt, state, dt);
    }
    let mut orbit = HbOrbit {
        omega: response,
        coeffs: vec![vec![C64::ZERO; ORDER + 1]; 3],
        residual_trace: vec![],
        residual: f64::NAN,
        aft_samples: 0,
        sampling_defect: f64::NAN,
    };
    let mut means = [0.0; 5];
    for sample in 0..steps {
        let time = (100 * steps + sample) as f64 * dt;
        let theta = TAU * sample as f64 / steps as f64;
        for (row, x) in orbit.coeffs.iter_mut().zip(state) {
            for (n, z) in row.iter_mut().enumerate() {
                *z = *z + complex(x / steps as f64, -(n as f64) * theta);
            }
        }
        let w = state[1] + ETA * MEMORY * (state[2] - state[0]);
        let nonlinear = BETA * (1.0 + SIGMA * state[0] * state[0]);
        let values = [
            -0.5 * PUMP * OMEGA * (OMEGA * time + 0.2).sin() * state[0] * state[0],
            BIAS * state[1],
            LOSS * state[1] * state[1] + RELAXATION * MEMORY * (state[2] - state[0]).powi(2),
            nonlinear * w.powi(4),
            -LOSS - RELAXATION - 3.0 * nonlinear * w * w * (1.0 + ETA * ETA * MEMORY),
        ];
        for (mean, value) in means.iter_mut().zip(values) {
            *mean += value / steps as f64;
        }
        state = rk4(&original_rhs, time, state, dt);
    }
    (orbit, means)
}

#[test]
fn nonlinear_loss_saturates_pumped_half_frequency_motion_with_original_heat_and_floquet_tangent() {
    let (system, input, pump) = period_components();
    let source = ParametricHamiltonian::new(&system, &input, &pump)
        .unwrap()
        .with_dissipation(law())
        .unwrap();
    let view = source.over_periods(2).unwrap();
    let omega = view.angular_frequency();
    let (coarse, _) = original_cycle(256);
    let mut seed = coarse.clone();
    for row in &mut seed.coeffs {
        row.truncate(5);
    }
    let orbit = view
        .solve_harmonic_balance_adaptive_scaled_seeded(
            seed.packed_seed(&HbAnchor::Forced { omega }),
            &HbBudget {
                harmonics: 4,
                max_newton: 28,
                tolerance: 2e-12,
                deflate_equilibrium: false,
            },
            &HbRefinementBudget {
                max_harmonics: 32,
                residual_tolerances: vec![2e-8; 3],
                state_tolerances: vec![2e-8; 3],
                frequency_tolerance: 1e-10,
            },
            &scales(3, omega),
        )
        .unwrap();
    assert!(orbit.refinement_trace.len() >= 2);
    let orbit = orbit.orbit;
    view.require_subharmonic_separation(&orbit, &[1e-6; 3])
        .unwrap();
    assert!(orbit.coeffs[0][1].abs() > 0.15);
    assert!(orbit.coeffs[0][3].abs() > 1e-3 && orbit.coeffs[0][5].abs() > 1e-5);
    assert!(orbit.coeffs[2][1].abs() > 0.03);
    let (middle, _) = original_cycle(512);
    let (fine, expected_means) = original_cycle(1024);
    let errors = [&coarse, &middle, &fine]
        .map(|reference| coefficient_error(&orbit.coeffs, &reference.coeffs));
    assert!(
        errors[1] < 0.15 * errors[0] && errors[2] < 0.2 * errors[1] && errors[2] < 3e-9,
        "{errors:?}"
    );
    let mut means = [0.0; 4];
    for sample in 0..1024 {
        let theta = TAU * sample as f64 / 1024.0;
        let x = core::array::from_fn::<_, 3, _>(|j| orbit.sample(j, theta));
        let power = view.power(theta / omega, &x).unwrap();
        close(
            power.energy_rate,
            power.pump_power + power.port_power - power.resistance_power - power.nonlinear_power,
            2e-13,
        );
        for (mean, value) in means.iter_mut().zip([
            power.pump_power,
            power.port_power,
            power.resistance_power,
            power.nonlinear_power,
        ]) {
            *mean += value / 1024.0;
        }
    }
    for j in 0..4 {
        close(means[j], expected_means[j], 2e-9);
    }
    assert!(means[0] > 0.01 && means[2] > 0.005 && means[3] > 0.005);
    close(means[0] + means[1], means[2] + means[3], 2e-10);
    let shot = view
        .solve_shooting_variational(
            &(0..3).map(|j| orbit.sample(j, 0.0)).collect::<Vec<_>>(),
            &ShootBudget {
                steps_per_period: 1024,
                max_newton: 8,
                tolerance: 1e-10,
            },
            &VariationalShootControls {
                state_scales: vec![1.0; 3],
                period_scale: TAU / omega,
                minimum_autonomous_excursion: 1e-5,
                maximum_dimension: 3,
                maximum_payload_bytes: 2_000_000,
                maximum_work: 2_000_000_000,
                maximum_period_maps: 32,
            },
            || false,
        )
        .unwrap();
    for j in 0..3 {
        close(shot.orbit.x0[j], orbit.sample(j, 0.0), 3e-8);
    }
    let a = &shot.monodromy;
    let determinant = a[0] * (a[4] * a[8] - a[5] * a[7]) - a[1] * (a[3] * a[8] - a[5] * a[6])
        + a[2] * (a[3] * a[7] - a[4] * a[6]);
    close(
        determinant,
        (expected_means[4] * shot.orbit.period).exp(),
        5e-9,
    );
    // The original dissipative tangent now reaches the shared sensitivity
    // owner on the complete two-source-period branch.
    let mut parameter_controls = parameters();
    // The native dealiased sensitivity grid needs at least 4*(h+1) samples;
    // the actual h32 branch therefore requires the next power of two.
    parameter_controls.samples_per_period = 256;
    let _sensitivity = view
        .prepare_sensitivity(&orbit, parameter_controls, &scales(3, omega), || false)
        .unwrap();
    assert!(matches!(
        view.prepare_sensitivity(&orbit, parameters(), &scales(3, omega), || true),
        Err(ParametricError::Sensitivity(
            ForcedSensitivityError::Cancelled
        ))
    ));
}

#[test]
fn nonlinear_loss_tangent_reaches_signed_pump_sidebands() {
    let (omega, carrier, resistance, linear_loss) = (0.7, 0.23, 0.3, 0.4);
    let make = |r| {
        PortHamiltonian::new(
            1,
            1,
            vec![0.0],
            vec![r],
            vec![1.0],
            Box::new(fs_phs::QuadraticStorage::new(vec![1.0], 1).unwrap()),
        )
        .unwrap()
    };
    let system = make(resistance);
    let equivalent = make(resistance + linear_loss);
    let input = PeriodicPortInput::new(omega, vec![vec![C64::ZERO]]).unwrap();
    let pump =
        PeriodicQuadraticModulation::new(omega, 1, vec![vec![C64::ZERO], vec![complex(0.04, 0.3)]])
            .unwrap();
    let d = |x: &[f64], e: &[f64], out: &mut [f64]| {
        out[0] = (linear_loss + 0.2 * x[0] * x[0]) * e[0] + 0.6 * e[0].powi(3);
        true
    };
    let dd = |x: &[f64], e: &[f64], dx: &[f64], de: &[f64], out: &mut [f64]| {
        out[0] = 0.4 * x[0] * dx[0] * e[0]
            + (linear_loss + 0.2 * x[0] * x[0] + 1.8 * e[0] * e[0]) * de[0];
        true
    };
    let source = ParametricHamiltonian::new(&system, &input, &pump)
        .unwrap()
        .with_dissipation(NonlinearDissipation {
            flow: &d,
            tangent: &dd,
            max_tangent_actions: 1,
        })
        .unwrap();
    let reference = ParametricHamiltonian::new(&equivalent, &input, &pump).unwrap();
    let order = 6;
    let orbit = hb(&source, 1, order);
    let mut forcing = vec![C64::ZERO; 2 * order + 1];
    forcing[order] = C64::new(0.13, -0.07);
    let scales = SidebandScales {
        states: vec![1.0],
        equations: vec![1.0],
    };
    let response = source
        .solve_port_sideband_response(
            &orbit,
            carrier,
            &forcing,
            sideband_controls(order),
            &scales,
            || false,
        )
        .unwrap();
    let expected = reference
        .solve_port_sideband_response(
            &orbit,
            carrier,
            &forcing,
            sideband_controls(order),
            &scales,
            || false,
        )
        .unwrap();
    for (a, b) in response.coefficients.iter().zip(expected.coefficients) {
        assert!((*a - b).abs() < 3e-13);
    }
    assert!(response.coefficients[order - 1].abs() > 1e-3);
    assert!(response.coefficients[order + 1].abs() > 1e-3);
}

#[test]
fn incomplete_or_active_loss_refuses_and_cancelled_sources_remain_reusable() {
    let (system, input, pump) = period_components();
    assert!(matches!(
        ParametricHamiltonian::new(&system, &input, &pump)
            .unwrap()
            .with_dissipation(NonlinearDissipation {
                max_tangent_actions: 2,
                ..law()
            }),
        Err(ParametricError::Invalid(_))
    ));
    let calls = core::cell::Cell::new(0);
    let refused = core::cell::Cell::new(true);
    let guarded = |x: &[f64], e: &[f64], out: &mut [f64]| {
        calls.set(calls.get() + 1);
        if refused.get() {
            out[0] = 0.0;
            true
        } else {
            flow(x, e, out)
        }
    };
    let source = ParametricHamiltonian::new(&system, &input, &pump)
        .unwrap()
        .with_dissipation(NonlinearDissipation {
            flow: &guarded,
            tangent: &tangent,
            max_tangent_actions: 3,
        })
        .unwrap();
    let budget = HbBudget {
        harmonics: 4,
        max_newton: 24,
        tolerance: 2e-12,
        deflate_equilibrium: false,
    };
    assert!(matches!(
        source.solve_harmonic_balance_scaled_seeded(vec![0.0; 27], &budget, &scales(3, OMEGA)),
        Err(ParametricError::Dissipation(
            LinearizationError::DissipationRefused { column: None }
        ))
    ));
    refused.set(false);
    assert!(
        source
            .solve_harmonic_balance_scaled_seeded(vec![0.0; 27], &budget, &scales(3, OMEGA))
            .is_ok()
    );
    let before = calls.get();
    let orbit = hb(&source, 3, 4);
    let baseline_calls = calls.get();
    assert!(baseline_calls > before);
    assert!(matches!(
        source.solve_port_sideband_response(
            &orbit,
            0.2,
            &[C64::ZERO; 9],
            sideband_controls(4),
            &SidebandScales {
                states: vec![1.0; 3],
                equations: vec![1.0; 3]
            },
            || true
        ),
        Err(ParametricError::Sideband(SidebandError::Cancelled))
    ));
    assert_eq!(calls.get(), baseline_calls);
    let mut field = [0.0; 3];
    source
        .vector_field_into(0.7, &[0.3, 0.4, 0.1], &mut field)
        .unwrap();
    let active = |_: &[f64], e: &[f64], out: &mut [f64]| {
        for (v, e) in out.iter_mut().zip(e) {
            *v = -*e;
        }
        true
    };
    let bad = ParametricHamiltonian::new(&system, &input, &pump)
        .unwrap()
        .with_dissipation(NonlinearDissipation {
            flow: &active,
            tangent: &tangent,
            max_tangent_actions: 3,
        })
        .unwrap();
    assert!(matches!(
        bad.power(0.3, &[0.4, 0.2, 0.1]),
        Err(ParametricError::Dissipation(
            LinearizationError::ActiveDissipation { .. }
        ))
    ));
}
