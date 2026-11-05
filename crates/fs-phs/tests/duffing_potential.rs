//! Physical spring work must survive unrelated stored energy and composition.
use fs_phs::{
    PortHamiltonian, QuadraticStorage, StepWorkspace, discrete_gradient, duffing_oscillator,
    interconnect, step,
};

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}, tolerance {tolerance:e}");
}

#[test]
fn spring_force_and_work_survive_large_kinetic_energy() {
    let system = duffing_oscillator(1.0, 2.0, 3.0, 0.0).unwrap();
    let (a, b) = ([1.0, 1e9], [1.25, 1e9]);
    // The 5e17 J kinetic background swallows this potential increment if
    // total endpoint energies are subtracted. These dyadic references follow
    // directly from V(q) = q^2 + 3q^4/4 at q=1 and q=5/4.
    assert_eq!(system.hamiltonian(&a), system.hamiltonian(&b));
    let effort = discrete_gradient(system.storage(), &a, &b);
    assert_eq!(effort, [6.57421875, 1e9]);
    assert_eq!(effort[0] * (b[0] - a[0]), 1.6435546875);
    assert_eq!(discrete_gradient(system.storage(), &b, &a), effort);

    let mid = [1.125, 1e9];
    let remainder = system.storage().midpoint_energy_remainder(&a, &b, &mid).unwrap();
    let change = system.effort(&mid)[0] * (b[0] - a[0]) + remainder;
    assert_eq!(change, 1.6435546875);
    assert_eq!(system.storage().midpoint_energy_remainder(&b, &a, &mid), Some(-remainder));
}

#[test]
fn separate_spring_and_momentum_work_and_coincident_efforts_are_preserved() {
    for k3 in [0.0, 3.0] {
        let system = duffing_oscillator(2.0, 2.0, k3, 0.0).unwrap();
        for (a, b) in [
            ([0.5, -1.25], [1.75, 0.75]),
            ([-1.5, 0.25], [1.5, -0.75]),
            ([1.0, 0.5], [1.0, -1.0]),
        ] {
            let effort = discrete_gradient(system.storage(), &a, &b);
            let potential = |q: f64| q * q + 0.25 * k3 * q.powi(4);
            near(effort[0] * (b[0] - a[0]), potential(b[0]) - potential(a[0]), 2e-14);
            near(effort[1] * (b[1] - a[1]), 0.25 * (b[1] * b[1] - a[1] * a[1]), 2e-14);
            assert_eq!(discrete_gradient(system.storage(), &b, &a), effort);
            assert_eq!(discrete_gradient(system.storage(), &a, &a), system.effort(&a));
        }
    }
}

#[test]
fn energy_remainder_retains_the_rounded_midpoint_offset() {
    let system = duffing_oscillator(4.0, 2.0, 3.0, 0.0).unwrap();
    let epsilon = f64::EPSILON;
    let (a, b) = ([1.0, 0.0], [1.0 + epsilon, 0.0]);
    let mid = [f64::midpoint(a[0], b[0]), 0.0];
    assert_eq!(mid, a);
    // Binomial expansion of V(1+eps)-V(1)-V'(1)*eps.
    let expected = 5.5 * epsilon.powi(2) + 3.0 * epsilon.powi(3) + 0.75 * epsilon.powi(4);
    near(system.storage().midpoint_energy_remainder(&a, &b, &mid).unwrap(), expected, 1e-45);

    let p = 9_007_199_254_740_992.0;
    assert_eq!(f64::midpoint(p, p + 2.0), p);
    // Kinetic Taylor remainder is (delta p)^2 / (2m), even though
    // the represented midpoint has rounded all the way to the left state.
    assert_eq!(system.storage().midpoint_energy_remainder(
        &[0.0, p], &[0.0, p + 2.0], &[0.0, p]), Some(0.5));
}

#[test]
fn analytic_tangents_differentiate_the_original_force_and_the_owned_secant() {
    let system = duffing_oscillator(0.7, 2.0, 3.0, 0.0).unwrap();
    let tangent = system.storage().continuous_tangent().unwrap();
    let rule = system.storage().discrete_gradient_rule().unwrap();
    let a = [-0.7, 0.3];
    for b in [[1.2, -0.4], a, [0.7, -0.3]] {
        for direction in [[0.2, -0.7], [-0.9, 0.4]] {
            let epsilon = 1e-6;
            let plus = [b[0] + epsilon * direction[0], b[1] + epsilon * direction[1]];
            let minus = [b[0] - epsilon * direction[0], b[1] - epsilon * direction[1]];
            let mut continuous = [f64::NAN; 2];
            let mut endpoint = [f64::NAN; 2];
            assert!(tangent.apply(&b, &direction, &mut continuous));
            assert!(rule.end_jacobian_action(&a, &b, &direction, &mut endpoint));
            let (force_plus, force_minus) = (system.effort(&plus), system.effort(&minus));
            let (secant_plus, secant_minus) = (
                discrete_gradient(system.storage(), &a, &plus),
                discrete_gradient(system.storage(), &a, &minus),
            );
            for i in 0..2 {
                near(continuous[i], (force_plus[i] - force_minus[i]) / (2.0 * epsilon), 3e-8);
                near(endpoint[i], (secant_plus[i] - secant_minus[i]) / (2.0 * epsilon), 3e-8);
                if b == a { near(endpoint[i], 0.5 * continuous[i], 2e-15); }
            }
        }
    }
    let mut out = [0.0; 2];
    assert!(!tangent.apply(&a[..1], &[1.0, 0.0], &mut out));
    assert!(!tangent.apply(&a, &[f64::NAN, 0.0], &mut out));
    assert!(!rule.end_jacobian_action(&a, &a[..1], &[1.0, 0.0], &mut out));
    assert!(!rule.end_jacobian_action(&a, &a, &[1.0, f64::INFINITY], &mut out));
}

fn potential(q: f64) -> f64 { q * q + 0.75 * q.powi(4) }

// Independent scalar solve after eliminating momentum. Evaluate work from
// the endpoint polynomial itself; do not call storage, its secant, or Newton.
fn scalar_reference(a: [f64; 2], force: f64, dt: f64) -> [f64; 2] {
    let secant = |q: f64| {
        if q == a[0] { 2.0 * q + 3.0 * q.powi(3) }
        else { (potential(q) - potential(a[0])) / (q - a[0]) }
    };
    let residual = |q: f64| {
        (1.4 + 0.15 * dt) * (q - a[0]) + dt * dt * secant(q)
            - 2.0 * dt * a[1] - dt * dt * force
    };
    let (mut lo, mut hi) = (-2.0, 2.0);
    assert!(residual(lo) < 0.0 && residual(hi) > 0.0);
    for _ in 0..80 {
        let q = f64::midpoint(lo, hi);
        if residual(q) < 0.0 { lo = q; } else { hi = q; }
    }
    let q = f64::midpoint(lo, hi);
    [q, a[1] + dt * (force - 0.15 * (q - a[0]) / dt - secant(q))]
}

#[test]
fn native_and_prepared_steps_match_independent_motion_spring_work_and_heat() {
    let system = duffing_oscillator(0.7, 2.0, 3.0, 0.15).unwrap();
    let mut prepared = StepWorkspace::new(&system).unwrap();
    prepared.set_iteration_limit(8).unwrap();
    let (force, dt) = (0.3, 0.03);
    let mut a = [0.8, -0.25];
    for _ in 0..20 {
        let expected = scalar_reference(a, force, dt);
        let ordinary = step(&system, &a, &[force], dt).unwrap();
        let (mut next, mut output) = ([0.0; 2], [0.0]);
        // A declared secant must use its own endpoint derivative, even when
        // an explicitly analytic caller supplies no Hamiltonian Hessian.
        let record = prepared.step_into_analytic(&system, &a, &[force], dt,
            &mut next, &mut output, &|_, _, _| panic!("owned secant has its own tangent"))
            .unwrap();
        assert!(record.newton_iters <= 8);
        for i in 0..2 {
            near(ordinary.x[i], expected[i], 2e-12);
            near(next[i], expected[i], 2e-12);
        }
        let velocity = (a[1] + next[1]) / 1.4;
        let energy_change = potential(next[0]) - potential(a[0])
            + (next[1] * next[1] - a[1] * a[1]) / 1.4;
        near(output[0], velocity, 2e-13);
        near(record.delta_h, energy_change, 2e-13);
        near(record.dissipated, dt * 0.15 * velocity * velocity, 2e-13);
        near(record.supplied, dt * force * velocity, 2e-13);
        near(record.balance_residual(), 0.0, 2e-13);
        a = next;
    }
}

#[test]
fn composition_keeps_small_spring_energy_changes_beside_frozen_storage() {
    let spring = duffing_oscillator(0.7, 2.0, 3.0, 0.15).unwrap();
    let frozen = PortHamiltonian::new(1, 0, vec![0.0], vec![0.0], vec![],
        Box::new(QuadraticStorage::new(vec![1.0], 1).unwrap())).unwrap();
    let system = interconnect(spring, frozen, &[]).unwrap();
    let a = [0.8, -0.25, 1e9];
    let expected = scalar_reference([a[0], a[1]], 0.3, 0.03);
    let ordinary = step(&system, &a, &[0.3], 0.03).unwrap();
    let mut prepared = StepWorkspace::new(&system).unwrap();
    prepared.set_iteration_limit(8).unwrap();
    let (mut next, mut output) = ([0.0; 3], [0.0]);
    let record = prepared.step_into(&system, &a, &[0.3], 0.03, &mut next, &mut output).unwrap();
    assert_eq!(next[2], a[2]);
    assert_eq!(system.hamiltonian(&a), system.hamiltonian(&next));
    assert!(record.delta_h.abs() > 1e-3);
    for i in 0..2 {
        near(next[i], expected[i], 2e-12);
        near(ordinary.x[i], expected[i], 2e-12);
    }
    let change = potential(next[0]) - potential(a[0])
        + (next[1] * next[1] - a[1] * a[1]) / 1.4;
    near(record.delta_h, change, 2e-13);
    near(ordinary.delta_h, change, 2e-13);
    near(record.balance_residual(), 0.0, 2e-13);
}

#[test]
fn nonfinite_coefficients_refuse_and_linear_limits_remain_admitted() {
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        for index in 0..4 {
            let mut coefficients = [1.0, 2.0, 3.0, 0.1];
            coefficients[index] = bad;
            assert!(duffing_oscillator(coefficients[0], coefficients[1], coefficients[2],
                coefficients[3]).is_err());
        }
    }
    assert!(duffing_oscillator(1.0, 0.0, 0.0, 0.0).is_ok());
    assert!(duffing_oscillator(1.0, 0.0, 3.0, 0.0).is_ok());
}
