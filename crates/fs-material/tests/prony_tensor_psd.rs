//! Physical Mandel work for correlated passive tensors across component scales.
use fs_material::visco::tensor::{Mandel6, PronyBranch, TensorProny, TensorPronyState};

fn identity(scale: f64) -> Mandel6 {
    core::array::from_fn(|i| core::array::from_fn(|j| if i == j { scale } else { 0.0 }))
}

fn correlated(extra_weak: bool) -> (Mandel6, f64) {
    let a = 2.0_f64.powi(-80);
    let mut c = [[0.0; 6]; 6];
    c[0][0] = 1.0;
    c[0][5] = a;
    c[5][0] = a;
    c[5][5] = if extra_weak { 2.0 * a * a } else { a * a };
    (c, a)
}

fn strain(a: f64, amplitude: f64) -> [f64; 6] {
    [-a * amplitude, 0.0, 0.0, 0.0, 0.0, amplitude]
}

fn dot(a: &[f64; 6], b: &[f64; 6]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn stress(c: &Mandel6, x: &[f64; 6]) -> [f64; 6] {
    core::array::from_fn(|i| dot(&c[i], x))
}

fn quadratic(c: &Mandel6, x: &[f64; 6]) -> f64 {
    dot(x, &stress(c, x))
}

fn close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 3e-11 * scale.max(1e-300),
        "actual={actual:e}, expected={expected:e}, scale={scale:e}"
    );
}

// Independent quadrature of q(t)^T C q(t)/tau in the ORIGINAL coordinates.
// No admitted factor or Prony heat routine participates in this oracle.
fn physical_heat(c: &Mandel6, q0: &[f64; 6], delta: &[f64; 6], dt: f64, tau: f64) -> f64 {
    let intervals = 512;
    let mut integral = 0.0;
    for k in 0..=intervals {
        let time = dt * f64::from(k) / f64::from(intervals);
        let decay = (-time / tau).exp();
        let ramp = -(-time / tau).exp_m1() / (dt / tau);
        let q = core::array::from_fn(|i| decay * q0[i] + ramp * delta[i]);
        let weight = if k == 0 || k == intervals {
            1.0
        } else if k % 2 == 0 {
            2.0
        } else {
            4.0
        };
        integral += weight * quadratic(c, &q) / tau;
    }
    integral * (dt / (3.0 * f64::from(intervals)))
}

fn physical_step(law: &TensorProny, state: &mut TensorPronyState, next: [f64; 6], dt: f64) {
    let before = state.clone();
    let delta = core::array::from_fn(|i| next[i] - before.strain[i]);
    let midpoint = core::array::from_fn(|i| before.strain[i] + 0.5 * delta[i]);
    let mut expected_stress = stress(law.c_eq(), &next);
    let mut expected_mean = stress(law.c_eq(), &midpoint);
    let mut expected_heat = 0.0;
    let mut expected_stored = 0.5 * quadratic(law.c_eq(), &next);
    let mut expected_q = Vec::new();
    if let Some(kv) = law.kv_pa_s() {
        let rate = delta.map(|d| d / dt);
        let viscous = stress(kv, &rate);
        for i in 0..6 {
            expected_stress[i] += viscous[i];
            expected_mean[i] += viscous[i];
        }
        expected_heat += quadratic(kv, &delta) / dt;
    }
    for (branch, q0) in law.branches().iter().zip(&before.q) {
        let r = dt / branch.tau_s;
        let decay = (-r).exp();
        let b = -(-r).exp_m1() / r;
        let c = (1.0 - b) / r;
        let q = core::array::from_fn(|i| decay * q0[i] + b * delta[i]);
        let mean_q = core::array::from_fn(|i| b * q0[i] + c * delta[i]);
        let endpoint = stress(&branch.c_pa, &q);
        let mean = stress(&branch.c_pa, &mean_q);
        for i in 0..6 {
            expected_stress[i] += endpoint[i];
            expected_mean[i] += mean[i];
        }
        expected_heat += physical_heat(&branch.c_pa, q0, &delta, dt, branch.tau_s);
        expected_stored += 0.5 * quadratic(&branch.c_pa, &q);
        expected_q.push(q);
    }
    let step = law.step_with_mean(state, next, dt).unwrap();
    for i in 0..6 {
        close(step.stress[i], expected_stress[i], expected_stress[i].abs());
        close(
            step.mean_stress[i],
            expected_mean[i],
            expected_mean[i].abs(),
        );
    }
    for (actual, expected) in state.q.iter().zip(expected_q) {
        for i in 0..6 {
            close(
                actual[i],
                expected[i],
                expected[i].abs() + before.q[0][i].abs(),
            );
        }
    }
    close(state.stored, expected_stored, expected_stored.abs());
    close(step.step_dissipated, expected_heat, expected_heat.abs());
    let expected_work = dot(&expected_mean, &delta);
    close(step.step_work, expected_work, expected_work.abs());
    close(
        state.stored - before.stored + step.step_dissipated,
        step.step_work,
        state.stored + before.stored + step.step_dissipated + step.step_work.abs(),
    );
}

#[test]
fn correlated_weak_prony_motion_retains_original_strain_energy_heat_and_work() {
    let (c, a) = correlated(true);
    let amplitude = 2.0_f64.powi(-10);
    // C = [1,a]^T[1,a] + diag(0,a^2). This strain annihilates the
    // strong channel, leaving the original weak loss a^2*amplitude^2.
    assert_eq!(
        quadratic(&c, &strain(a, amplitude)),
        a * a * amplitude * amplitude
    );
    for cached in [false, true] {
        let branch = if cached {
            PronyBranch::new(c, 0.7).unwrap()
        } else {
            PronyBranch::unvalidated(c, 0.7)
        };
        // The exact zero arm stays in source order and keeps its strain history.
        let law = TensorProny::new(
            identity(a * a),
            None,
            vec![branch, PronyBranch::unvalidated([[0.0; 6]; 6], 0.11)],
            (0.0, 1e6),
        )
        .unwrap();
        assert_eq!(law.branches()[0].c_pa, c);
        assert_eq!(law.branches()[0].tau_s, 0.7);
        assert_eq!(law.branches()[1].tau_s, 0.11);
        let mut state = law.state();
        for (gain, dt) in [(1.0, 0.3), (0.4, 0.04), (-0.7, 0.8), (0.0, 0.12)] {
            physical_step(&law, &mut state, strain(a, amplitude * gain), dt);
        }
        assert!(state.dissipated > 0.0);
    }
}

#[test]
fn rank_one_prony_and_kelvin_tensors_preserve_the_exact_physical_null() {
    let (c, a) = correlated(false);
    let amplitude = 2.0_f64.powi(-10);
    let law = TensorProny::new(
        identity(a * a),
        Some(c),
        vec![PronyBranch::new(c, 0.7).unwrap()],
        (0.0, 1e6),
    )
    .unwrap();
    let mut state = law.state();
    for (gain, dt) in [(1.0, 0.3), (-0.5, 0.2), (0.0, 0.4)] {
        physical_step(&law, &mut state, strain(a, amplitude * gain), dt);
        assert_eq!(
            state.dissipated, 0.0,
            "source null direction must not acquire heat"
        );
    }
    assert_eq!(state.stored, 0.0);
}

#[test]
fn kelvin_heat_uses_the_same_correlated_mandel_tensor_on_both_attachment_paths() {
    let (c, a) = correlated(true);
    for attached in [false, true] {
        let law = if attached {
            TensorProny::new(identity(a * a), None, vec![], (0.0, 1e6))
                .unwrap()
                .with_kelvin_viscosity_mandel(c)
                .unwrap()
        } else {
            TensorProny::new(identity(a * a), Some(c), vec![], (0.0, 1e6)).unwrap()
        };
        assert_eq!(law.kv_pa_s(), Some(&c));
        let mut state = law.state();
        physical_step(&law, &mut state, strain(a, 2.0_f64.powi(-10)), 0.3);
        assert!(state.dissipated > 0.0);
    }
}

#[test]
fn representable_large_and_small_passive_tensors_reach_actual_ramp_physics() {
    for exponent in [600, -600] {
        let scale = 2.0_f64.powi(exponent);
        let amplitude = if exponent > 0 {
            2.0_f64.powi(-300)
        } else {
            2.0_f64.powi(-10)
        };
        let mut c = [[0.0; 6]; 6];
        c[3][3] = scale;
        let law = TensorProny::new(
            identity(1.0),
            Some(c),
            vec![PronyBranch::unvalidated(c, 0.7)],
            (0.0, 1e6),
        )
        .unwrap();
        let mut state = law.state();
        let mut next = [0.0; 6];
        next[3] = amplitude;
        physical_step(&law, &mut state, next, 0.3);
        physical_step(&law, &mut state, [0.0; 6], 0.2);
        // At zero total strain, this is the independently computed branch
        // energy alone, including the uniformly tiny admitted source tensor.
        assert!(state.stored > 0.0 && state.dissipated > 0.0);
        assert_eq!(law.branches()[0].c_pa, c);
    }
}

#[test]
fn active_cross_terms_cannot_hide_behind_a_strong_unrelated_scale() {
    let (mut active, a) = correlated(false);
    active[0][5] = 2.0 * a;
    active[5][0] = 2.0 * a;
    assert!(quadratic(&active, &strain(a, 1.0)) < 0.0);
    assert!(PronyBranch::new(active, 0.7).is_err());
    assert!(
        TensorProny::new(
            identity(1.0),
            None,
            vec![PronyBranch::unvalidated(active, 0.7)],
            (0.0, 1e6)
        )
        .is_err()
    );
    assert!(TensorProny::new(identity(1.0), Some(active), vec![], (0.0, 1e6)).is_err());
}
