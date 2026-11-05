use super::*;
use crate::phs_orbit::{
    audit_forced_orbit, solve_forced_harmonic_balance_scaled_seeded,
    solve_forced_shooting_variational,
};
use core::cell::Cell;
use fs_math::{c64::C64, det};
use fs_nlmodal::{SosModalStorage, StressChannel};
use fs_orbit::{OrbitError, OrbitProblem, ShootBudget, VariationalShootControls};
use fs_phs::{QuadraticStorage, Storage};
use std::rc::Rc;

fn scales(omega: f64) -> HbScales {
    HbScales {
        states: vec![1.0, 2.0],
        equations: vec![2.0, 4.0],
        angular_frequency: omega,
    }
}

struct CountedQuadratic {
    storage: QuadraticStorage,
    gradient_calls: Rc<Cell<usize>>,
    declare: bool,
}
impl Storage for CountedQuadratic {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        self.storage.hamiltonian(x)
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        self.gradient_calls.set(self.gradient_calls.get() + 1);
        self.storage.gradient(x, out);
    }
    fn exact_quadratic_matrix(&self) -> Option<&[f64]> {
        self.declare.then(|| self.storage.matrix())
    }
}

#[test]
fn declared_original_quadratic_tangent_avoids_fourier_column_source_sweeps() {
    let make = |declare| {
        let count = Rc::new(Cell::new(0));
        let storage = CountedQuadratic {
            storage: QuadraticStorage::new(vec![4.0, 0.0, 0.0, 0.5], 2).unwrap(),
            gradient_calls: count.clone(),
            declare,
        };
        let system = PortHamiltonian::new(
            2,
            1,
            vec![0.0, 1.0, -1.0, 0.0],
            vec![0.0, 0.0, 0.0, 0.3],
            vec![0.0, 1.0],
            Box::new(storage),
        )
        .unwrap();
        count.set(0);
        (system, count)
    };
    let (exact, exact_calls) = make(true);
    let (difference, difference_calls) = make(false);
    let omega = 1.3;
    let input = PeriodicPortInput::new(
        omega,
        vec![vec![
            C64::from_re(0.2),
            C64::new(0.3, -0.2),
            C64::ZERO,
            C64::new(-0.07, 0.03),
        ]],
    )
    .unwrap();
    let budget = HbBudget {
        harmonics: 7,
        tolerance: 1e-11,
        ..HbBudget::default()
    };
    let units = scales(omega);
    let a =
        solve_forced_harmonic_balance_scaled_seeded(&exact, &input, vec![0.0; 30], &budget, &units)
            .unwrap();
    let b = solve_forced_harmonic_balance_scaled_seeded(
        &difference,
        &input,
        vec![0.0; 30],
        &budget,
        &units,
    )
    .unwrap();
    for (n, force) in input.coefficients[0].iter().enumerate() {
        let w = n as f64 * omega;
        let q = *force / C64::new(4.0 - 2.0 * w * w, 0.3 * w);
        let p = C64::new(0.0, 2.0 * w) * q;
        for (state, expected) in [(0, q), (1, p)] {
            assert!((a.orbit.coeffs[state][n] - expected).abs() < 1e-10);
            assert!((b.orbit.coeffs[state][n] - expected).abs() < 1e-10);
        }
    }
    assert!(
        exact_calls.get() * 5 < difference_calls.get(),
        "original gradient evaluations: analytic {}, finite-difference {}",
        exact_calls.get(),
        difference_calls.get()
    );
}

struct CubicStorage {
    stiffness: f64,
    cubic: f64,
    mass: f64,
}
impl Storage for CubicStorage {
    fn quadratic_part(&self) -> Option<fs_phs::linear_exact::QuadraticPart> {
        Some(fs_phs::linear_exact::QuadraticPart {
            matrix: vec![self.stiffness, 0.0, 0.0, 1.0 / self.mass],
            exact: self.cubic == 0.0,
        })
    }
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        0.5 * self.stiffness * x[0] * x[0]
            + 0.25 * self.cubic * x[0].powi(4)
            + 0.5 * x[1] * x[1] / self.mass
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        out[0] = self.stiffness * x[0] + self.cubic * x[0].powi(3);
        out[1] = x[1] / self.mass;
    }
}

#[test]
fn original_nonlinear_hessian_recovers_manufactured_multiharmonic_force_and_work() {
    let (mass, stiffness, cubic, damping, omega, amplitude) = (2.0, 4.0, 10.0, 0.3, 1.3, 0.4);
    let system = PortHamiltonian::new(
        2,
        1,
        vec![0.0, 1.0, -1.0, 0.0],
        vec![0.0, 0.0, 0.0, damping],
        vec![0.0, 1.0],
        Box::new(CubicStorage {
            stiffness,
            cubic,
            mass,
        }),
    )
    .unwrap();
    let (sin_phase, cos_phase) = det::sin_cos(0.37);
    let q1 = C64::new(0.5 * amplitude * cos_phase, 0.5 * amplitude * sin_phase);
    let p1 = C64::new(0.0, mass * omega) * q1;
    // Independently substitute q(t)=A cos(omega*t+phase) into
    // m q'' + d q' + k q + alpha q^3 = F(t).
    let f1 = q1
        * C64::new(
            stiffness - mass * omega * omega + 0.75 * cubic * amplitude * amplitude,
            damping * omega,
        );
    let f3 = (q1 * q1 * q1).scale(cubic);
    let input = PeriodicPortInput::new(omega, vec![vec![C64::ZERO, f1, C64::ZERO, f3]]).unwrap();
    assert!(
        forced_problem(&system, &input)
            .unwrap()
            .island_state_jacobian()
            .is_none(),
        "the quadratic part cannot stand in for the nonlinear Hessian"
    );
    let calls = Cell::new(0);
    let hessian = |x: &[f64], v: &[f64], out: &mut [f64]| {
        calls.set(calls.get() + 1);
        out[0] = (stiffness + 3.0 * cubic * x[0] * x[0]) * v[0];
        out[1] = v[1] / mass;
        true
    };
    let budget = HbBudget {
        harmonics: 3,
        tolerance: 1e-11,
        ..HbBudget::default()
    };
    let units = scales(omega);
    let mut seed = vec![0.0; 14];
    seed[2] = 0.8 * q1.re;
    seed[3] = 0.8 * p1.re;
    seed[4] = 0.8 * q1.im;
    seed[5] = 0.8 * p1.im;
    let fixed = solve_forced_harmonic_balance_analytic_scaled_seeded(
        &system,
        &input,
        seed.clone(),
        &budget,
        &units,
        &hessian,
    )
    .unwrap();
    let adaptive = solve_forced_harmonic_balance_analytic_adaptive_scaled_seeded(
        &system,
        &input,
        seed,
        &budget,
        &HbRefinementBudget {
            max_harmonics: 7,
            residual_tolerances: vec![1e-8; 2],
            state_tolerances: vec![1e-8; 2],
            frequency_tolerance: 1e-10,
        },
        &units,
        &hessian,
    )
    .unwrap();
    assert!(calls.get() > 0);
    assert_eq!(
        adaptive
            .refinement_trace
            .iter()
            .map(|x| x.harmonics)
            .collect::<Vec<_>>(),
        [3, 7]
    );
    let exact_loss = 0.5 * damping * omega * omega * amplitude * amplitude;
    for orbit in [&fixed.orbit, &adaptive.orbit] {
        assert!((orbit.coeffs[0][1] - q1).abs() < 1e-9);
        assert!((orbit.coeffs[1][1] - p1).abs() < 1e-9);
        for state in &orbit.coeffs {
            for (n, coefficient) in state.iter().enumerate() {
                if n != 1 {
                    assert!(coefficient.abs() < 1e-9);
                }
            }
        }
        let audit = audit_forced_orbit(&system, &input, orbit, 256).unwrap();
        assert!((audit.mean_resistance_loss_w / exact_loss - 1.0).abs() < 1e-8);
        assert!(audit.relative_power_residual < 1e-8);
    }
}

#[test]
fn assembled_modal_storage_supplies_hb_and_variational_tangents_without_a_second_law() {
    let (stiffness, cubic, damping, omega, amplitude): (f64, f64, f64, f64, f64) =
        (4.0, 10.0, 0.3, 1.3, 0.4);
    let storage = SosModalStorage {
        omegas: vec![stiffness.sqrt()],
        channels: vec![StressChannel {
            coefficient: cubic,
            coupling: vec![1.0],
        }],
    };
    let system =
        fs_nlmodal::assemble(storage, &[damping / (2.0 * stiffness.sqrt())], &[1.0]).unwrap();
    let (sin_phase, cos_phase) = det::sin_cos(0.37);
    let q1 = C64::new(0.5 * amplitude * cos_phase, 0.5 * amplitude * sin_phase);
    let p1 = C64::new(0.0, omega) * q1;
    // Substituting q=A*cos(omega*t+phase) into the ORIGINAL unit-mass law
    // requires both fundamental and third-harmonic physical forcing.
    let f1 = q1
        * C64::new(
            stiffness - omega * omega + 0.75 * cubic * amplitude * amplitude,
            damping * omega,
        );
    let f3 = (q1 * q1 * q1).scale(cubic);
    let input = PeriodicPortInput::new(omega, vec![vec![C64::ZERO, f1, C64::ZERO, f3]]).unwrap();
    let problem = forced_problem(&system, &input).unwrap();
    let action = problem
        .island_state_jacobian()
        .expect("boxed native modal tangent");
    let mut tangent = [f64::NAN; 4];
    assert!(action(&problem, 0.0, &[0.3, -0.2], &mut tangent));
    for (actual, expected) in tangent.iter().zip([0.0, 1.0, -6.7, -damping]) {
        assert!((actual - expected).abs() < 1e-14);
    }
    let budget = HbBudget {
        harmonics: 3,
        tolerance: 1e-11,
        ..HbBudget::default()
    };
    let mut seed = vec![0.0; 14];
    seed[2] = 0.8 * q1.re;
    seed[3] = 0.8 * p1.re;
    seed[4] = 0.8 * q1.im;
    seed[5] = 0.8 * p1.im;
    let response =
        solve_forced_harmonic_balance_scaled_seeded(&system, &input, seed, &budget, &scales(omega))
            .unwrap();
    for (state, expected) in [(0, q1), (1, p1)] {
        assert!((response.orbit.coeffs[state][1] - expected).abs() < 1e-9);
        for (n, coefficient) in response.orbit.coeffs[state].iter().enumerate() {
            if n != 1 {
                assert!(coefficient.abs() < 1e-9);
            }
        }
    }
    let audit = audit_forced_orbit(&system, &input, &response.orbit, 256).unwrap();
    let loss = 0.5 * damping * omega * omega * amplitude * amplitude;
    assert!((audit.mean_resistance_loss_w / loss - 1.0).abs() < 1e-8);
    assert!(audit.relative_power_residual < 1e-8);

    // The ordinary variational owner discovers the SAME boxed action and
    // propagates all state directions at its actual RK4 stages.
    let period = 2.0 * core::f64::consts::PI / omega;
    let expected = [2.0 * q1.re, 2.0 * p1.re];
    let shot = solve_forced_shooting_variational(
        &system,
        &input,
        &expected,
        &ShootBudget {
            steps_per_period: 512,
            max_newton: 12,
            tolerance: 1e-11,
        },
        &VariationalShootControls {
            state_scales: vec![1.0, 2.0],
            period_scale: period,
            minimum_autonomous_excursion: 0.01,
            maximum_dimension: 2,
            maximum_payload_bytes: 100_000,
            maximum_work: 100_000_000,
            maximum_period_maps: 32,
        },
        || false,
    )
    .unwrap();
    for (actual, expected) in shot.orbit.x0.iter().zip(expected) {
        assert!((actual - expected).abs() < 2e-7);
    }
    assert!(shot.jacobian_evaluations > 0);
    let m = &shot.monodromy;
    let determinant = m[0] * m[3] - m[1] * m[2];
    // The continuous vector-field divergence is exactly -damping even
    // though the quartic stiffness varies throughout the original orbit.
    assert!((determinant / (-damping * period).exp() - 1.0).abs() < 3e-7);
}

struct RefusingNativeStorage {
    incomplete: bool,
}
impl Storage for RefusingNativeStorage {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        0.5 * (x[0] * x[0] + x[1] * x[1])
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        out.copy_from_slice(x);
    }
    fn continuous_tangent(&self) -> Option<&dyn ContinuousStorageTangent> {
        Some(self)
    }
}
impl ContinuousStorageTangent for RefusingNativeStorage {
    fn apply(&self, _: &[f64], _: &[f64], out: &mut [f64]) -> bool {
        if self.incomplete {
            out[0] = 1.0;
        }
        self.incomplete
    }
}

#[test]
fn native_storage_tangent_refusal_cannot_fall_back_to_numerical_columns() {
    let input = PeriodicPortInput::new(1.3, vec![vec![C64::ZERO, C64::from_re(0.5)]]).unwrap();
    let budget = HbBudget {
        harmonics: 1,
        ..HbBudget::default()
    };
    for incomplete in [false, true] {
        let system = PortHamiltonian::new(
            2,
            1,
            vec![0.0, 1.0, -1.0, 0.0],
            vec![0.0, 0.0, 0.0, 0.3],
            vec![0.0, 1.0],
            Box::new(RefusingNativeStorage { incomplete }),
        )
        .unwrap();
        assert!(matches!(
            solve_forced_harmonic_balance_scaled_seeded(
                &system,
                &input,
                vec![0.0; 6],
                &budget,
                &scales(1.3),
            ),
            Err(PhsOrbitError::Orbit(OrbitError::TangentRefused { .. }))
        ));
    }
}

#[test]
fn original_relaxation_coordinates_supply_the_complete_analytic_tangent() {
    let (mass, stiffness, damping, arm_stiffness, tau) = (2.0, 4.0, 0.3, 3.0, 0.4);
    let system = fs_phs::mass_spring_damper(mass, stiffness, damping)
        .unwrap()
        .with_relaxation_branches(vec![fs_phs::RelaxationBranch {
            projection: vec![1.0, 0.0],
            stiffness: arm_stiffness,
            relaxation_time_s: tau,
        }])
        .unwrap();
    let input = PeriodicPortInput::new(1.3, vec![vec![C64::ZERO, C64::from_re(0.5)]]).unwrap();
    let problem = forced_problem(&system, &input).unwrap();
    let callback = problem
        .island_state_jacobian()
        .expect("complete original quadratic memory law");
    let mut actual = vec![f64::NAN; 9];
    assert!(callback(&problem, 0.0, &[0.2, -0.1, 0.05], &mut actual));
    // The original memory variable stores H_arm=(sqrt(k_arm)*q-z)^2/2.
    // q'=p/m, p'=-(k+k_arm)q-d*p/m+sqrt(k_arm)z,
    // z'=(sqrt(k_arm)q-z)/tau; all three coordinates stay in the tangent.
    let coupling = arm_stiffness.sqrt();
    let expected = [
        0.0,
        1.0 / mass,
        0.0,
        -(stiffness + arm_stiffness),
        -damping / mass,
        coupling,
        coupling / tau,
        0.0,
        -1.0 / tau,
    ];
    for (actual, expected) in actual.iter().zip(expected) {
        assert!((actual - expected).abs() < 1e-13);
    }
}

#[test]
fn declared_storage_hessian_refusal_and_incomplete_output_cannot_fall_back() {
    let system = fs_phs::mass_spring_damper(2.0, 4.0, 0.3).unwrap();
    let input = PeriodicPortInput::new(1.3, vec![vec![C64::ZERO, C64::from_re(0.5)]]).unwrap();
    let budget = HbBudget {
        harmonics: 1,
        ..HbBudget::default()
    };
    let units = scales(1.3);
    for mode in 0..3 {
        let action = |_: &[f64], _: &[f64], out: &mut [f64]| {
            if mode == 0 {
                return false;
            }
            if mode == 1 {
                out.fill(f64::NAN);
            } else {
                out[0] = 1.0;
            }
            true
        };
        assert!(matches!(
            solve_forced_harmonic_balance_analytic_scaled_seeded(
                &system,
                &input,
                vec![0.0; 6],
                &budget,
                &units,
                &action,
            ),
            Err(PhsOrbitError::Orbit(OrbitError::TangentRefused { .. }))
        ));
    }
}
