//! Real finite-amplitude period doubling, unchanged source laws and full work.

use super::*;
use crate::phs_orbit::PhsForcedParameterDirection;
use fs_orbit::sensitivity::{
    ForcedMeanObjectiveControls, ForcedSensitivityControls, ForcedSensitivityDirection,
    ForcedSensitivityError,
};

const ORDER: usize = 12;
const LOSS: f64 = 0.08;
const MEMORY: f64 = 0.08;
const RELAXATION: f64 = 0.4;
const BIAS: f64 = 0.03;

struct PeriodMemory;
impl Storage for PeriodMemory {
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
impl ContinuousStorageTangent for PeriodMemory {
    fn apply(&self, x: &[f64], v: &[f64], out: &mut [f64]) -> bool {
        out[0] = (1.0 + 3.0 * x[0] * x[0] + MEMORY) * v[0] - MEMORY * v[2];
        out[1] = v[1];
        out[2] = MEMORY * (v[2] - v[0]);
        true
    }
}
fn components(
    omega: f64,
) -> (
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
        Box::new(PeriodMemory),
    )
    .unwrap();
    let input = PeriodicPortInput::new(omega, vec![vec![C64::from_re(BIAS)]]).unwrap();
    let mut q = vec![vec![C64::ZERO; 9]; 2];
    q[1][0] = complex(0.225, 0.2);
    (
        system,
        input,
        PeriodicQuadraticModulation::new(omega, 3, q).unwrap(),
    )
}
fn original_rhs(time: f64, x: &[f64; 3], omega: f64) -> [f64; 3] {
    [
        x[1],
        -(1.0 + 0.45 * (omega * time + 0.2).cos()) * x[0]
            - x[0].powi(3)
            - MEMORY * (x[0] - x[2])
            - LOSS * x[1]
            + BIAS,
        RELAXATION * (x[0] - x[2]),
    ]
}
fn empty_orbit(n: usize, order: usize, omega: f64) -> HbOrbit {
    HbOrbit {
        omega,
        coeffs: vec![vec![C64::ZERO; order + 1]; n],
        residual_trace: vec![],
        residual: f64::NAN,
        aft_samples: 0,
        sampling_defect: f64::NAN,
    }
}
// Explicit test-only physical ODE; a nonzero initial state selects one of the
// saturated branches. No production field, HB, pump sampler or tangent is used.
fn timed_cycle(omega: f64, steps: usize) -> (HbOrbit, [f64; 3]) {
    let response = omega / 2.0;
    let dt = TAU / (response * steps as f64);
    let rhs = |time: f64, x: &[f64; 3]| original_rhs(time, x, omega);
    let mut state = [0.2, 0.1, 0.0];
    for sample in 0..100 * steps {
        state = rk4(&rhs, sample as f64 * dt, state, dt);
    }
    let mut orbit = empty_orbit(3, ORDER, response);
    let mut powers = [0.0; 3];
    for sample in 0..steps {
        let theta = TAU * sample as f64 / steps as f64;
        let time = (100 * steps + sample) as f64 * dt;
        for (row, x) in orbit.coeffs.iter_mut().zip(state) {
            for (n, z) in row.iter_mut().enumerate() {
                *z = *z + complex(x / steps as f64, -(n as f64) * theta);
            }
        }
        powers[0] +=
            -0.5 * omega * 0.45 * (omega * time + 0.2).sin() * state[0] * state[0] / steps as f64;
        powers[1] += BIAS * state[1] / steps as f64;
        powers[2] += (LOSS * state[1] * state[1]
            + RELAXATION * MEMORY * (state[2] - state[0]).powi(2))
            / steps as f64;
        state = rk4(&rhs, time, state, dt);
    }
    (orbit, powers)
}
fn budget(order: usize) -> HbBudget {
    HbBudget {
        harmonics: order,
        max_newton: 24,
        tolerance: 2e-12,
        deflate_equilibrium: false,
    }
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
fn objectives() -> ForcedMeanObjectiveControls {
    ForcedMeanObjectiveControls {
        objective_scale: 1.0,
        maximum_value_sampling_defect: 1e-10,
        maximum_gradient_sampling_defect: 1e-9,
    }
}
fn empty_direction<'a>() -> ParametricParameterDirection<'a> {
    ParametricParameterDirection {
        physical: PhsForcedParameterDirection {
            structure: Default::default(),
            storage: None,
            input: None,
        },
        modulation: None,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        a.is_finite() && (a - b).abs() <= tolerance * (1.0 + b.abs()),
        "{a:.17e} != {b:.17e}"
    );
}

#[test]
fn saturated_half_frequency_response_keeps_memory_and_full_cycle_pump_work() {
    let omega = 2.1;
    let (system, input, pump) = components(omega);
    let source = ParametricHamiltonian::new(&system, &input, &pump).unwrap();
    let view = source.over_periods(2).unwrap();
    assert_eq!(view.angular_frequency(), omega / 2.0);
    assert_eq!(view.required_harmonics(), 2);
    let (seed, _) = timed_cycle(omega, 256);
    let anchor = HbAnchor::Forced {
        omega: view.angular_frequency(),
    };
    let orbit = view
        .solve_harmonic_balance_scaled_seeded(
            seed.packed_seed(&anchor),
            &budget(ORDER),
            &scales(3, view.angular_frequency()),
        )
        .unwrap()
        .orbit;
    view.require_subharmonic_separation(&orbit, &[1e-6; 3])
        .unwrap();
    assert!(orbit.coeffs[0][1].abs() > 0.2);
    assert!(orbit.coeffs[2][1].abs() > 0.05);
    // DC bias removes inversion symmetry; the full original tangent need not
    // repeat on the pump period even though the modulation itself does.
    assert!(orbit.coeffs[0][0].re.abs() > 1e-3);
    let (fine, means) = timed_cycle(omega, 512);
    let coarse_error = coefficient_error(&orbit.coeffs, &seed.coeffs);
    let fine_error = coefficient_error(&orbit.coeffs, &fine.coeffs);
    assert!(
        fine_error < 2e-8 && fine_error < 0.25 * coarse_error,
        "coarse={coarse_error:e}, fine={fine_error:e}"
    );
    let mut prepared = view
        .prepare_sensitivity(
            &orbit,
            parameters(),
            &scales(3, view.angular_frequency()),
            || false,
        )
        .unwrap();
    for (index, objective) in [
        ParametricPowerObjective::MeanPumpPower,
        ParametricPowerObjective::MeanPortPower,
        ParametricPowerObjective::MeanResistancePower,
        ParametricPowerObjective::PumpEnergyPerCycle,
        ParametricPowerObjective::PortEnergyPerCycle,
        ParametricPowerObjective::ResistanceEnergyPerCycle,
    ]
    .into_iter()
    .enumerate()
    {
        let value = prepared
            .power_objective(objective, objectives(), || false)
            .unwrap()
            .value;
        close(
            value,
            means[index % 3]
                * if index >= 3 {
                    TAU / view.angular_frequency()
                } else {
                    1.0
                },
            3e-8,
        );
    }
    assert!(means[0] > 0.01 && means[2] > 0.01);
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
                period_scale: TAU / view.angular_frequency(),
                minimum_autonomous_excursion: 1e-5,
                maximum_dimension: 3,
                maximum_payload_bytes: 2_000_000,
                maximum_work: 2_000_000_000,
                maximum_period_maps: 32,
            },
            || false,
        )
        .unwrap();
    close(shot.orbit.period, TAU / view.angular_frequency(), 1e-15);
    for j in 0..3 {
        close(shot.orbit.x0[j], orbit.sample(j, 0.0), 3e-8);
    }
    let a = &shot.monodromy;
    let determinant = a[0] * (a[4] * a[8] - a[5] * a[7]) - a[1] * (a[3] * a[8] - a[5] * a[6])
        + a[2] * (a[3] * a[7] - a[4] * a[6]);
    close(
        determinant,
        (-(LOSS + RELAXATION) * shot.orbit.period).exp(),
        3e-9,
    );
    assert!(
        shot.orbit
            .multipliers
            .iter()
            .all(|&(re, im)| re.hypot(im) < 0.99)
    );
    // The original one-pump-period branch also fits the long chart. Its
    // genuine original-source solve must NOT count as a period-doubled result.
    let shorter = hb(&source, 3, 5);
    assert!(
        source
            .require_subharmonic_separation(&shorter, &[1e-6; 3])
            .is_err()
    );
    let mut repeated = empty_orbit(3, ORDER, view.angular_frequency());
    for (row, original) in repeated.coeffs.iter_mut().zip(&shorter.coeffs) {
        for (n, &z) in original.iter().enumerate() {
            row[2 * n] = z;
        }
    }
    assert!(
        view.require_subharmonic_separation(&repeated, &[1e-6; 3])
            .is_err()
    );
}

#[test]
fn period_views_preserve_exact_physical_laws_and_authored_tangent() {
    let omega = 0.8;
    let system = memory_system();
    let input = drive(omega, 0.2, -0.2);
    let pump = stiffness(omega, 3, 0.25, 0.4);
    let calls = Cell::new(0);
    let hessian = |x: &[f64], v: &[f64], out: &mut [f64]| {
        calls.set(calls.get() + 1);
        MemoryStorage.apply(x, v, out)
    };
    let source = ParametricHamiltonian::new_analytic(&system, &input, &pump, &hessian).unwrap();
    let triple = source.over_periods(3).unwrap();
    let nested = triple.over_periods(2).unwrap();
    assert_eq!(nested.period_multiplier(), 6);
    assert_eq!(nested.angular_frequency(), omega / 6.0);
    assert_eq!(nested.required_harmonics(), 12);
    assert_eq!(calls.get(), 0);
    let bits = |p: ParametricPower| {
        [
            p.energy,
            p.pump_power,
            p.port_power,
            p.resistance_power,
            p.energy_rate,
        ]
        .map(f64::to_bits)
    };
    let state = [0.31, -0.27, 0.08];
    for time in [0.0, 0.43, 13.17, 1e6] {
        let mut original = [0.0; 3];
        source
            .vector_field_into(time, &state, &mut original)
            .unwrap();
        let original_power = bits(source.power(time, &state).unwrap());
        let mut jacobian = [0.0; 9];
        source
            .state_jacobian_into(time, &state, &mut jacobian)
            .unwrap();
        for view in [&triple, &nested] {
            let mut field = [0.0; 3];
            view.vector_field_into(time, &state, &mut field).unwrap();
            assert_eq!(original.map(f64::to_bits), field.map(f64::to_bits));
            assert_eq!(original_power, bits(view.power(time, &state).unwrap()));
            let before = calls.get();
            let mut actual = [0.0; 9];
            view.state_jacobian_into(time, &state, &mut actual).unwrap();
            assert_eq!(calls.get(), before + 3);
            assert_eq!(jacobian.map(f64::to_bits), actual.map(f64::to_bits));
        }
    }
    assert!(source.over_periods(0).is_err());
    assert!(source.over_periods(257).is_err());
    assert!(source.over_periods(129).is_err()); // actual pump harmonic 2 maps above 256
    assert!(triple.over_periods(100).is_err());
    // Full-state separation counts explicit memory, ignores huge DC, and
    // cannot promote a tiny odd tail or exact shorter synchronized repeat.
    let mut candidate = empty_orbit(3, 12, triple.angular_frequency());
    candidate.coeffs[0][0] = C64::from_re(1e6);
    candidate.coeffs[0][3] = C64::from_re(0.4);
    candidate.coeffs[2][1] = C64::new(0.03, 0.04);
    triple
        .require_subharmonic_separation(&candidate, &[1e-8; 3])
        .unwrap();
    candidate.coeffs[2][1] = C64::from_re(1e-12);
    assert!(
        triple
            .require_subharmonic_separation(&candidate, &[1e-8; 3])
            .is_err()
    );
}

fn mechanical_source(
    omega: f64,
) -> (
    PortHamiltonian,
    PeriodicPortInput,
    PeriodicQuadraticModulation,
) {
    (
        mass_spring_damper(1.3, 2.7, 0.35).unwrap(),
        drive(omega, 0.23, -0.24),
        stiffness(omega, 2, 0.55, 0.41),
    )
}
fn mechanical_mean_pump(omega: f64, orbit: &HbOrbit) -> f64 {
    (0..1024)
        .map(|sample| {
            let theta = TAU * sample as f64 / 1024.0;
            let time = theta / orbit.omega;
            -0.55 * omega * (2.0 * omega * time + 0.41).sin() * orbit.sample(0, theta).powi(2)
                / 1024.0
        })
        .sum()
}

#[test]
fn response_clock_gradients_and_direction_bandwidth_use_the_long_chart() {
    let omega = 0.9;
    let (system, input, pump) = mechanical_source(omega);
    let source = ParametricHamiltonian::new(&system, &input, &pump).unwrap();
    let view = source.over_periods(2).unwrap();
    let order = 18;
    let orbit = view
        .solve_harmonic_balance_scaled_seeded(
            vec![0.0; 2 * (2 * order + 1)],
            &budget(order),
            &scales(2, view.angular_frequency()),
        )
        .unwrap()
        .orbit;
    let mut prepared = view
        .prepare_sensitivity(
            &orbit,
            parameters(),
            &scales(2, view.angular_frequency()),
            || false,
        )
        .unwrap();
    let chart = ForcedSensitivityDirection {
        parameter: 0.0,
        parameter_scale: 0.4,
        angular_frequency_derivative: 0.1,
    };
    let response = prepared
        .solve_parameter(chart, empty_direction(), || false)
        .unwrap();
    let mut packed = vec![0.0; 2 * (2 * order + 1)];
    for state in 0..2 {
        packed[state] = response.coeffs[state][0].re;
    }
    for n in 1..=order {
        for state in 0..2 {
            let re = 2 + (n - 1) * 4 + state;
            packed[re] = response.coeffs[state][n].re;
            packed[re + 2] = response.coeffs[state][n].im;
        }
    }
    let epsilon = 2e-5;
    let neighbor = |change: f64| {
        // The parameter derivative denotes the RESPONSE clock: physical
        // input and pump clocks change by twice that amount on this view.
        let w = omega + 2.0 * chart.angular_frequency_derivative * change;
        let (system, input, pump) = mechanical_source(w);
        let source = ParametricHamiltonian::new(&system, &input, &pump).unwrap();
        let view = source.over_periods(2).unwrap();
        let anchor = HbAnchor::Forced {
            omega: view.angular_frequency(),
        };
        let solved = view
            .solve_harmonic_balance_scaled_seeded(
                orbit.packed_seed(&anchor),
                &budget(order),
                &scales(2, view.angular_frequency()),
            )
            .unwrap()
            .orbit;
        (w, mechanical_mean_pump(w, &solved))
    };
    let (lo, hi) = (neighbor(-epsilon), neighbor(epsilon));
    for cycle in [false, true] {
        let objective = if cycle {
            ParametricPowerObjective::PumpEnergyPerCycle
        } else {
            ParametricPowerObjective::MeanPumpPower
        };
        let observation = prepared
            .power_objective(objective, objectives(), || false)
            .unwrap();
        let chain: f64 = observation
            .gradient
            .iter()
            .zip(&packed)
            .map(|(g, x)| g * x)
            .sum();
        let explicit = if cycle {
            0.0
        } else {
            chart.angular_frequency_derivative / view.angular_frequency() * observation.value
        };
        let mut adjoint = prepared
            .prepare_power_adjoint(objective, objectives(), || false)
            .unwrap();
        let result = adjoint
            .gradient_parameter(chart, empty_direction(), || false)
            .unwrap();
        close(result.derivative - chain, explicit, 3e-12);
        let expected = if cycle {
            2.0 * TAU * (hi.1 / hi.0 - lo.1 / lo.0) / (2.0 * epsilon)
        } else {
            (hi.1 - lo.1) / (2.0 * epsilon)
        };
        close(result.derivative, expected, 2e-7);
    }
    let mut high = vec![vec![C64::ZERO; 4]; 11];
    high[10][0] = C64::from_re(1e-12);
    let high = PeriodicQuadraticModulation::new(omega, 2, high).unwrap();
    // Physical harmonic 10 is below retained order 18 numerically, but its
    // RESPONSE-chart index is 20 and must refuse, even for this tiny amplitude.
    assert!(matches!(
        prepared.solve_parameter(
            chart,
            ParametricParameterDirection {
                modulation: Some(&high),
                ..empty_direction()
            },
            || false
        ),
        Err(ParametricError::Invalid(_))
    ));
    let wrong =
        PeriodicQuadraticModulation::new(view.angular_frequency(), 2, vec![vec![C64::ZERO; 4]])
            .unwrap();
    assert!(
        prepared
            .solve_parameter(
                chart,
                ParametricParameterDirection {
                    modulation: Some(&wrong),
                    ..empty_direction()
                },
                || false
            )
            .is_err()
    );
    assert!(matches!(
        prepared.solve_parameter(chart, empty_direction(), || true),
        Err(ParametricError::Sensitivity(
            ForcedSensitivityError::Cancelled
        ))
    ));
    let retry = prepared
        .solve_parameter(chart, empty_direction(), || false)
        .unwrap();
    assert_eq!(response.coeffs, retry.coeffs);
}
