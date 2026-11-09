use super::*;
use fs_phs::{ContinuousStorageTangent, Storage, mass_spring_damper};

#[path = "sensitivity/tests.rs"]
mod sensitivity;

#[path = "periods/tests.rs"]
mod periods;

#[path = "dissipation/tests.rs"]
mod dissipation;

const TAU: f64 = core::f64::consts::TAU;

fn complex(amplitude: f64, phase: f64) -> C64 {
    C64::new(amplitude * phase.cos(), amplitude * phase.sin())
}
fn drive(omega: f64, peak: f64, phase: f64) -> PeriodicPortInput {
    PeriodicPortInput::new(omega, vec![vec![C64::ZERO, complex(0.5 * peak, phase)]]).unwrap()
}
fn stiffness(omega: f64, n: usize, peak: f64, phase: f64) -> PeriodicQuadraticModulation {
    let mut matrices = vec![vec![C64::ZERO; n * n]; 3];
    matrices[2][0] = complex(0.5 * peak, phase);
    PeriodicQuadraticModulation::new(omega, n, matrices).unwrap()
}
fn scales(n: usize, omega: f64) -> HbScales {
    HbScales {
        states: vec![1.0; n],
        equations: vec![1.0; n],
        angular_frequency: omega,
    }
}
fn hb(problem: &ParametricHamiltonian<'_>, n: usize, order: usize) -> HbOrbit {
    problem
        .solve_harmonic_balance_scaled_seeded(
            vec![0.0; n * (2 * order + 1)],
            &HbBudget {
                harmonics: order,
                max_newton: 18,
                tolerance: 2e-12,
                deflate_equilibrium: false,
            },
            &scales(n, problem.modulation.omega),
        )
        .unwrap()
        .orbit
}

// Independent test-only time-domain quadrature of explicitly written physical
// equations. It does not call the adapter, HB, shooting or original field owner.
fn rk4<const D: usize>(
    rhs: &impl Fn(f64, &[f64; D]) -> [f64; D],
    time: f64,
    x: [f64; D],
    dt: f64,
) -> [f64; D] {
    let k1 = rhs(time, &x);
    let k2 = rhs(
        time + 0.5 * dt,
        &core::array::from_fn(|j| x[j] + 0.5 * dt * k1[j]),
    );
    let k3 = rhs(
        time + 0.5 * dt,
        &core::array::from_fn(|j| x[j] + 0.5 * dt * k2[j]),
    );
    let k4 = rhs(time + dt, &core::array::from_fn(|j| x[j] + dt * k3[j]));
    core::array::from_fn(|j| x[j] + dt * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]) / 6.0)
}
fn time_coefficients<const D: usize>(
    omega: f64,
    steps: usize,
    order: usize,
    rhs: impl Fn(f64, &[f64; D]) -> [f64; D],
) -> Vec<Vec<C64>> {
    let dt = TAU / (omega * steps as f64);
    let mut x = [0.0; D];
    for i in 0..120 * steps {
        x = rk4(&rhs, i as f64 * dt, x, dt);
    }
    let mut result = vec![vec![C64::ZERO; order + 1]; D];
    for i in 0..steps {
        let time = (120 * steps + i) as f64 * dt;
        for n in 0..=order {
            let phase = n as f64 * omega * time;
            for j in 0..D {
                result[j][n] = result[j][n] + complex(x[j] / steps as f64, -phase);
            }
        }
        x = rk4(&rhs, time, x, dt);
    }
    result
}
fn coefficient_error(a: &[Vec<C64>], b: &[Vec<C64>]) -> f64 {
    a.iter()
        .zip(b)
        .flat_map(|(a, b)| a.iter().zip(b))
        .map(|(a, b)| (*a - *b).abs())
        .fold(0.0, f64::max)
}

#[test]
fn mathieu_harmonics_agree_with_independent_refined_time_domain_motion() {
    let (m, k, c, omega, delta, phase, force, drive_phase) =
        (1.3, 2.7, 0.35, 0.9, 0.55, 0.41, 0.23, -0.24);
    let system = mass_spring_damper(m, k, c).unwrap();
    let input = drive(omega, force, drive_phase);
    let pump = stiffness(omega, 2, delta, phase);
    let problem = ParametricHamiltonian::new(&system, &input, &pump).unwrap();
    let orbit = hb(&problem, 2, 9);
    let rhs = |time: f64, x: &[f64; 2]| {
        [
            x[1] / m,
            -(k + delta * (2.0 * omega * time + phase).cos()) * x[0] - c * x[1] / m
                + force * (omega * time + drive_phase).cos(),
        ]
    };
    let coarse = time_coefficients(omega, 512, 9, rhs);
    let fine = time_coefficients(omega, 1024, 9, rhs);
    let ec = coefficient_error(&orbit.coeffs, &coarse);
    let ef = coefficient_error(&orbit.coeffs, &fine);
    assert!(ef < 3e-9 && ef < 0.15 * ec, "coarse={ec:e}, fine={ef:e}");
    assert!(orbit.coeffs[0][3].abs() > 1e-4);
    assert!(orbit.coeffs[0][5].abs() > 1e-7);
    assert!(orbit.coeffs[0][2].abs() < 1e-12);
}

#[test]
fn explicit_pump_work_closes_the_physical_energy_balance() {
    let (m, k, c, omega, delta, phase, force, drive_phase) =
        (1.3, 2.7, 0.35, 0.9, 0.55, 0.41, 0.23, -0.24);
    let system = mass_spring_damper(m, k, c).unwrap();
    let input = drive(omega, force, drive_phase);
    let pump = stiffness(omega, 2, delta, phase);
    let problem = ParametricHamiltonian::new(&system, &input, &pump).unwrap();
    let rhs = |time: f64, x: &[f64; 3]| {
        let velocity = x[1] / m;
        let pump_power = -delta * omega * (2.0 * omega * time + phase).sin() * x[0] * x[0];
        let supplied = force * (omega * time + drive_phase).cos() * velocity;
        [
            velocity,
            -(k + delta * (2.0 * omega * time + phase).cos()) * x[0] - c * velocity
                + force * (omega * time + drive_phase).cos(),
            pump_power + supplied - c * velocity * velocity,
        ]
    };
    let energy = |time: f64, x: &[f64; 3]| {
        0.5 * ((k + delta * (2.0 * omega * time + phase).cos()) * x[0] * x[0] + x[1] * x[1] / m)
    };
    let mut errors = Vec::new();
    for steps in [256, 512] {
        let dt = TAU / (omega * steps as f64);
        let mut state = [0.4, -0.15, 0.0];
        let initial = energy(0.0, &state);
        for j in 0..steps {
            let time = j as f64 * dt;
            let power = problem.power(time, &state[..2]).unwrap();
            assert!((power.energy - energy(time, &state)).abs() < 2e-14);
            let exact_pump =
                -delta * omega * (2.0 * omega * time + phase).sin() * state[0] * state[0];
            assert!((power.pump_power - exact_pump).abs() < 2e-14);
            assert!(
                (power.energy_rate - power.pump_power - power.port_power + power.resistance_power)
                    .abs()
                    < 3e-14
            );
            state = rk4(&rhs, time, state, dt);
        }
        errors.push((energy(TAU / omega, &state) - initial - state[2]).abs());
    }
    assert!(
        errors[1] < 0.1 * errors[0] && errors[1] < 2e-9,
        "{errors:?}"
    );
}

// A complete nonlinear original storage with an explicit relaxation coordinate.
// The memory z is neither eliminated nor folded into a fitted damping constant.
struct MemoryStorage;
impl Storage for MemoryStorage {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        0.5 * (1.8 * x[0] * x[0] + x[1] * x[1] / 1.2 + 0.8 * (x[2] - x[0]).powi(2))
            + 0.6 * x[0].powi(4) / 4.0
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        out[0] = 1.8 * x[0] + 0.6 * x[0].powi(3) + 0.8 * (x[0] - x[2]);
        out[1] = x[1] / 1.2;
        out[2] = 0.8 * (x[2] - x[0]);
    }
    fn continuous_tangent(&self) -> Option<&dyn ContinuousStorageTangent> {
        Some(self)
    }
}
impl ContinuousStorageTangent for MemoryStorage {
    fn apply(&self, x: &[f64], v: &[f64], out: &mut [f64]) -> bool {
        if x.len() != 3 || v.len() != 3 || out.len() != 3 {
            return false;
        }
        out[0] = (2.6 + 1.8 * x[0] * x[0]) * v[0] - 0.8 * v[2];
        out[1] = v[1] / 1.2;
        out[2] = 0.8 * (v[2] - v[0]);
        true
    }
}
fn memory_system() -> PortHamiltonian {
    PortHamiltonian::new(
        3,
        1,
        vec![0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        vec![0.0, 0.0, 0.0, 0.0, 0.32, 0.0, 0.0, 0.0, 0.7],
        vec![0.0, 1.0, 0.0],
        Box::new(MemoryStorage),
    )
    .unwrap()
}

#[test]
fn nonlinear_storage_memory_and_mixed_energy_tangents_remain_complete() {
    let omega = 0.8;
    let system = memory_system();
    let input = drive(omega, 0.2, -0.2);
    let mut matrices = vec![vec![C64::ZERO; 9]; 3];
    matrices[0][0] = C64::new(0.07, 0.0);
    matrices[1][2] = complex(0.06, 0.3);
    matrices[1][6] = matrices[1][2];
    matrices[2][0] = complex(0.125, 0.4);
    let pump = PeriodicQuadraticModulation::new(omega, 3, matrices).unwrap();
    let problem = ParametricHamiltonian::new(&system, &input, &pump).unwrap();
    let (time, state) = (0.73, [0.31, -0.27, 0.08]);
    let cross = 0.12 * (omega * time + 0.3).cos();
    let cross_rate = -0.12 * omega * (omega * time + 0.3).sin();
    let stiffness = 0.07 + 0.25 * (2.0 * omega * time + 0.4).cos();
    let stiffness_rate = -0.5 * omega * (2.0 * omega * time + 0.4).sin();
    let e0 = (1.8 + stiffness) * state[0]
        + 0.6 * state[0].powi(3)
        + 0.8 * (state[0] - state[2])
        + cross * state[2];
    let e1 = state[1] / 1.2;
    let e2 = 0.8 * (state[2] - state[0]) + cross * state[0];
    let expected = [
        e1,
        -e0 - 0.32 * e1 + 0.2 * (omega * time - 0.2).cos(),
        -0.7 * e2,
    ];
    let mut actual = [0.0; 3];
    problem
        .vector_field_into(time, &state, &mut actual)
        .unwrap();
    for j in 0..3 {
        assert!((actual[j] - expected[j]).abs() < 2e-14);
    }
    let mut jacobian = [0.0; 9];
    problem
        .state_jacobian_into(time, &state, &mut jacobian)
        .unwrap();
    for column in 0..3 {
        let mut minus = state;
        let mut plus = state;
        minus[column] -= 1e-6;
        plus[column] += 1e-6;
        let mut fm = [0.0; 3];
        let mut fp = [0.0; 3];
        problem.vector_field_into(time, &minus, &mut fm).unwrap();
        problem.vector_field_into(time, &plus, &mut fp).unwrap();
        for row in 0..3 {
            assert!((jacobian[row * 3 + column] - (fp[row] - fm[row]) / 2e-6).abs() < 2e-9);
        }
    }
    let p = problem.power(time, &state).unwrap();
    assert!(
        (p.pump_power
            - 0.5 * stiffness_rate * state[0] * state[0]
            - cross_rate * state[0] * state[2])
            .abs()
            < 2e-14
    );
    assert!((p.energy_rate - p.pump_power - p.port_power + p.resistance_power).abs() < 2e-14);
}

#[test]
fn adaptive_modulated_nonlinear_response_preserves_relaxation_harmonics() {
    let omega = 0.8;
    let system = memory_system();
    let input = drive(omega, 0.2, -0.2);
    let pump = stiffness(omega, 3, 0.25, 0.4);
    let problem = ParametricHamiltonian::new(&system, &input, &pump).unwrap();
    let result = problem
        .solve_harmonic_balance_adaptive_scaled_seeded(
            vec![0.0; 3 * 7],
            &HbBudget {
                harmonics: 3,
                max_newton: 20,
                tolerance: 2e-12,
                deflate_equilibrium: false,
            },
            &HbRefinementBudget {
                max_harmonics: 15,
                residual_tolerances: vec![2e-8; 3],
                state_tolerances: vec![2e-8; 3],
                frequency_tolerance: 1e-10,
            },
            &scales(3, omega),
        )
        .unwrap();
    assert!(result.refinement_trace.len() >= 2);
    let order = result.orbit.coeffs[0].len() - 1;
    let direct = time_coefficients(omega, 1024, order, |time, x: &[f64; 3]| {
        let q = x[0];
        let p = x[1];
        let z = x[2];
        [
            p / 1.2,
            -(1.8 + 0.25 * (2.0 * omega * time + 0.4).cos()) * q
                - 0.6 * q * q * q
                - 0.8 * (q - z)
                - 0.32 * p / 1.2
                + 0.2 * (omega * time - 0.2).cos(),
            -0.56 * (z - q),
        ]
    });
    assert!(coefficient_error(&result.orbit.coeffs, &direct) < 3e-9);
    assert!(result.orbit.coeffs.len() == 3 && result.orbit.coeffs[0][3].abs() > 1e-5);
    for n in [1, 3, 5] {
        let expected =
            C64::new(0.56, 0.0) * result.orbit.coeffs[0][n] / C64::new(0.56, n as f64 * omega);
        assert!((result.orbit.coeffs[2][n] - expected).abs() < 2e-11);
    }
}

#[test]
fn same_parametric_tangent_drives_variational_shooting_and_liouville_identity() {
    let (m, k, c, omega) = (1.3, 2.7, 0.35, 0.9);
    let system = mass_spring_damper(m, k, c).unwrap();
    let input = drive(omega, 0.23, -0.24);
    let pump = stiffness(omega, 2, 0.55, 0.41);
    let problem = ParametricHamiltonian::new(&system, &input, &pump).unwrap();
    let orbit = hb(&problem, 2, 9);
    let shoot = problem
        .solve_shooting_variational(
            &[0.0, 0.0],
            &ShootBudget {
                steps_per_period: 1024,
                max_newton: 8,
                tolerance: 1e-10,
            },
            &VariationalShootControls {
                state_scales: vec![1.0; 2],
                period_scale: TAU / omega,
                minimum_autonomous_excursion: 1e-5,
                maximum_dimension: 2,
                maximum_payload_bytes: 2_000_000,
                maximum_work: 2_000_000_000,
                maximum_period_maps: 32,
            },
            || false,
        )
        .unwrap();
    for j in 0..2 {
        assert!((shoot.orbit.x0[j] - orbit.sample(j, 0.0)).abs() < 2e-8);
    }
    let a = &shoot.monodromy;
    let determinant = a[0] * a[3] - a[1] * a[2];
    // Abel/Liouville: trace of the complete continuous field Jacobian is -c/m.
    assert!((determinant - (-c * TAU / (m * omega)).exp()).abs() < 3e-10);
    assert!(shoot.jacobian_evaluations > 0);
}

fn sideband_controls(order: usize) -> SidebandControls {
    SidebandControls {
        harmonics: order,
        samples_per_period: 128,
        maximum_unknowns: 128,
        maximum_payload_bytes: 8_000_000,
        maximum_work: 2_000_000_000,
        maximum_base_residual: 1e-10,
        maximum_tangent_sampling_defect: 1e-10,
        maximum_relative_residual: 1e-10,
        maximum_omitted_residual: 1e-9,
    }
}
#[test]
fn pump_converts_physical_port_probes_into_both_signed_sidebands() {
    let (m, k, c, omega, delta, phase, carrier) = (1.3, 4.0, 0.3, 0.9, 0.001, 0.41, 0.25);
    let system = mass_spring_damper(m, k, c).unwrap();
    let input = drive(omega, 0.0, 0.0);
    let pump = stiffness(omega, 2, delta, phase);
    let problem = ParametricHamiltonian::new(&system, &input, &pump).unwrap();
    let order = 6;
    let orbit = hb(&problem, 2, order);
    let mut forcing = vec![C64::ZERO; 2 * order + 1];
    forcing[order] = C64::new(0.13, -0.07);
    let controls = sideband_controls(order);
    let response = problem
        .solve_port_sideband_response(
            &orbit,
            carrier,
            &forcing,
            controls,
            &SidebandScales {
                states: vec![1.0; 2],
                equations: vec![1.0; 2],
            },
            || false,
        )
        .unwrap();
    let dynamic = |frequency: f64| C64::new(k - m * frequency * frequency, c * frequency);
    let unpumped = forcing[order] / dynamic(carrier);
    for sign in [-1, 1] {
        let band = (order as isize + 2 * sign) as usize;
        let frequency = carrier + 2.0 * sign as f64 * omega;
        let expected = -complex(0.5 * delta, sign as f64 * phase) * unpumped / dynamic(frequency);
        let q = response.coefficients[2 * band];
        let p = response.coefficients[2 * band + 1];
        assert!(q.abs() > 1e-7);
        assert!((q - expected).abs() < 2e-6 * expected.abs());
        assert!((p - C64::new(0.0, m * frequency) * q).abs() < 2e-13);
    }
    assert!(response.relative_residual < 1e-10 && response.omitted_residual < 1e-9);
    let mut tiny = controls;
    tiny.maximum_payload_bytes = 1;
    assert!(matches!(
        problem.solve_port_sideband_response(
            &orbit,
            carrier,
            &forcing,
            tiny,
            &SidebandScales {
                states: vec![1.0; 2],
                equations: vec![1.0; 2]
            },
            || false
        ),
        Err(ParametricError::Sideband(SidebandError::Budget { .. }))
    ));
}

#[test]
fn physical_fourier_conventions_bandwidth_and_declared_tangent_refusals_are_explicit() {
    let system = mass_spring_damper(1.0, 2.0, 0.2).unwrap();
    let input = drive(1.0, 0.2, 0.0);
    let mut matrix = vec![C64::ZERO; 4];
    matrix[1] = C64::new(0.1, 0.2);
    matrix[2] = matrix[1];
    assert!(
        PeriodicQuadraticModulation::new(1.0, 2, vec![vec![C64::ZERO; 4], matrix.clone()]).is_ok()
    );
    matrix[2] = matrix[1].conj();
    assert!(matches!(
        PeriodicQuadraticModulation::new(1.0, 2, vec![vec![C64::ZERO; 4], matrix]),
        Err(ParametricError::Invalid(_))
    ));
    assert!(PeriodicQuadraticModulation::new(1.0, 1, vec![vec![C64::new(0.0, 1e-200)]]).is_err());
    let mut coefficients = vec![vec![C64::ZERO; 4]; 9];
    coefficients[8][0] = C64::new(1e-200, 0.0);
    let wide = PeriodicQuadraticModulation::new(1.0, 2, coefficients).unwrap();
    let problem = ParametricHamiltonian::new(&system, &input, &wide).unwrap();
    assert!(matches!(
        problem.solve_harmonic_balance_scaled_seeded(
            vec![0.0; 30],
            &HbBudget {
                harmonics: 7,
                ..HbBudget::default()
            },
            &scales(2, 1.0)
        ),
        Err(ParametricError::Invalid(_))
    ));
    let wrong_clock = stiffness(2.0, 2, 0.1, 0.0);
    assert!(ParametricHamiltonian::new(&system, &input, &wrong_clock).is_err());
    let zero = stiffness(1.0, 2, 0.0, 0.0);
    let calls = Cell::new(0);
    let refusing = |_: &[f64], _: &[f64], _: &mut [f64]| {
        calls.set(calls.get() + 1);
        false
    };
    let problem = ParametricHamiltonian::new_analytic(&system, &input, &zero, &refusing).unwrap();
    assert!(matches!(
        problem.solve_harmonic_balance_scaled_seeded(
            vec![0.0; 14],
            &HbBudget {
                harmonics: 3,
                ..HbBudget::default()
            },
            &scales(2, 1.0)
        ),
        Err(ParametricError::Orbit(OrbitError::TangentRefused { .. }))
    ));
    assert!(calls.get() > 0);
}
