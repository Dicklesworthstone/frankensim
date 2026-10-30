//! G0/G1/G3: full nonlinear storage with exact linear and material-memory flow.
//! Heat uses the original resistance on the method's discrete effort; only the
//! quadratic limit is an exact continuous physical trajectory.
use fs_phs::linear_exact::{
    LinearSplitBudget, PreparedExactZoh, PreparedNonlinearZoh, QuadraticPart,
};
use fs_phs::{PortHamiltonian, QuadraticStorage, RelaxationBranch, Storage, duffing_oscillator};

fn budget(n: usize) -> LinearSplitBudget {
    LinearSplitBudget {
        max_modes: n,
        max_omega_dt: 0.04,
    }
}

struct QuarticStorage {
    base: QuadraticStorage,
    cubic: [f64; 2],
    datum: f64,
}
impl Storage for QuarticStorage {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        self.datum
            + self.base.hamiltonian(x)
            + 0.25 * (self.cubic[0] * x[0].powi(4) + self.cubic[1] * x[2].powi(4))
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        self.base.gradient(x, out);
        out[0] += self.cubic[0] * x[0].powi(3);
        out[2] += self.cubic[1] * x[2].powi(3);
    }
    fn quadratic_part(&self) -> Option<QuadraticPart> {
        let mut part = self.base.quadratic_part()?;
        part.exact = self.cubic == [0.0, 0.0];
        Some(part)
    }
    fn midpoint_energy_remainder(&self, a: &[f64], b: &[f64], mid: &[f64]) -> Option<f64> {
        let mut remainder = self.base.midpoint_energy_remainder(a, b, mid)?;
        for (index, cubic) in [(0, self.cubic[0]), (2, self.cubic[1])] {
            let d = b[index] - a[index];
            let m = mid[index];
            let offset = 0.5 * ((a[index] - m) + (b[index] - m));
            remainder += cubic
                * d
                * (offset * (3.0 * m * m + 3.0 * m * offset + offset * offset)
                    + 0.25 * (m + offset) * d * d);
        }
        Some(remainder)
    }
}

fn body(cubic: [f64; 2], datum: f64) -> PortHamiltonian {
    let mut j = vec![0.0; 16];
    j[1] = 1.0;
    j[4] = -1.0;
    j[11] = 1.0;
    j[14] = -1.0;
    let mut r = vec![0.0; 16];
    // Positive definite but not diagonally dominant: all cross-mode damping
    // terms must survive the exponential and independently integrated heat.
    r[5] = 0.4;
    r[7] = 0.5;
    r[13] = 0.5;
    r[15] = 0.8;
    let base = QuadraticStorage::new(
        vec![
            9.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 25.0, 0.0, 0.0, 0.0, 0.0, 1.0,
        ],
        4,
    )
    .unwrap();
    PortHamiltonian::new(
        4,
        2,
        j,
        r,
        vec![0.0, 0.0, 1.0, 0.3, 0.0, 0.0, -0.4, 0.9],
        Box::new(QuarticStorage { base, cubic, datum }),
    )
    .unwrap()
    .with_relaxation_branches(vec![RelaxationBranch {
        projection: vec![1.0, 0.0, -0.5, 0.0],
        stiffness: 4.0,
        relaxation_time_s: 1e-4,
    }])
    .unwrap()
}

fn rk4<const N: usize>(
    mut x: [f64; N],
    duration: f64,
    steps: usize,
    rhs: impl Fn(&[f64; N]) -> [f64; N],
) -> [f64; N] {
    let h = duration / steps as f64;
    for _ in 0..steps {
        let a = rhs(&x);
        let b = rhs(&core::array::from_fn(|i| x[i] + 0.5 * h * a[i]));
        let c = rhs(&core::array::from_fn(|i| x[i] + 0.5 * h * b[i]));
        let d = rhs(&core::array::from_fn(|i| x[i] + h * c[i]));
        for i in 0..N {
            x[i] += h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    x
}

fn close(a: f64, b: f64, scale: f64) {
    assert!(
        (a - b).abs() <= 3e-10 * scale.max(f64::MIN_POSITIVE),
        "{a:.17e} != {b:.17e}, scale={scale:.17e}"
    );
}

#[test]
fn exact_quadratic_limit_preserves_full_damping_and_stiff_warm_memory() {
    let original = body([0.0, 0.0], 0.0);
    let expected_structure = {
        let (j, r, g) = original.structure();
        (j.to_vec(), r.to_vec(), g.to_vec())
    };
    let dt = 0.01; // One hundred original Maxwell relaxation times.
    let exact = PreparedExactZoh::new_owned(original, dt, budget(5)).unwrap();
    let mut nonlinear =
        PreparedNonlinearZoh::new_owned(body([0.0, 0.0], 0.0), dt, budget(5)).unwrap();
    assert_eq!(
        nonlinear.system().structure(),
        (
            expected_structure.0.as_slice(),
            expected_structure.1.as_slice(),
            expected_structure.2.as_slice()
        )
    );
    let initial = [0.3, -0.15, -0.2, 0.1, 0.65];
    let input = [0.23, -0.17];
    let reference = exact.step(&initial, &input).unwrap();
    let actual = nonlinear.step(&initial, &input).unwrap();
    for (a, b) in actual.x.iter().zip(&reference.x) {
        close(*a, *b, 1.0);
    }
    for (a, b) in actual.y.iter().zip(&reference.y) {
        close(*a, *b, 1.0);
    }
    close(actual.delta_h, reference.delta_h, 1.0);
    close(actual.dissipated, reference.dissipated, 1.0);
    close(actual.supplied, reference.supplied, 1.0);
    assert!(actual.dissipated > 0.01);
}

#[test]
fn coupled_quartic_modes_and_stiff_memory_refine_to_independent_continuous_physics() {
    let initial = [0.3, -0.15, -0.2, 0.1, 0.65];
    let input = [0.23, -0.17];
    let duration = 0.08;
    // This oracle evaluates the original continuous polynomial and Maxwell
    // equations directly, including the full physical damping tensor.
    let reference = rk4(
        [
            initial[0], initial[1], initial[2], initial[3], initial[4], 0.0, 0.0,
        ],
        duration,
        65_536,
        |x| {
            let arm = 2.0 * x[0] - x[2] - x[4];
            let heat = 0.4 * x[1] * x[1] + x[1] * x[3] + 0.8 * x[3] * x[3] + arm * arm / 1e-4;
            let work = input[0] * (x[1] - 0.4 * x[3]) + input[1] * (0.3 * x[1] + 0.9 * x[3]);
            [
                x[1],
                -9.0 * x[0] - 7.0 * x[0].powi(3) - 2.0 * arm - 0.4 * x[1] - 0.5 * x[3]
                    + input[0]
                    + 0.3 * input[1],
                x[3],
                -25.0 * x[2] - 3.0 * x[2].powi(3) + arm - 0.5 * x[1] - 0.8 * x[3] - 0.4 * input[0]
                    + 0.9 * input[1],
                arm / 1e-4,
                heat,
                work,
            ]
        },
    );
    let mut errors = Vec::new();
    for steps in [8, 16, 32] {
        let dt = duration / f64::from(steps);
        let mut runtime =
            PreparedNonlinearZoh::new_owned(body([7.0, 3.0], 0.0), dt, budget(5)).unwrap();
        let before = runtime.system().hamiltonian(&initial);
        let mut x = initial;
        let mut next = [0.0; 5];
        let mut output = [0.0; 2];
        let (mut heat, mut work) = (0.0, 0.0);
        for _ in 0..steps {
            let step = runtime
                .step_into(&x, &input, &mut next, &mut output)
                .unwrap();
            assert!(step.dissipated >= 0.0);
            close(
                step.supplied,
                dt * (input[0] * output[0] + input[1] * output[1]),
                1.0,
            );
            heat += step.dissipated;
            work += step.supplied;
            x = next;
        }
        close(
            runtime.system().hamiltonian(&x) - before + heat,
            work,
            before,
        );
        let weights = [9.0, 1.0, 25.0, 1.0, 1.0];
        errors.push(
            (0..5)
                .map(|i| weights[i] * (x[i] - reference[i]).powi(2))
                .sum::<f64>()
                .sqrt(),
        );
        if steps == 32 {
            assert!((heat - reference[5]).abs() < 3e-6 * reference[5].abs());
            assert!((work - reference[6]).abs() < 3e-6 * reference[6].abs());
        }
    }
    assert!(errors[2] < 2e-5, "continuous full-state errors {errors:?}");
    for pair in errors.windows(2) {
        assert!(
            pair[0] > 2.8 * pair[1],
            "nonlinear refinement errors {errors:?}"
        );
    }
}

#[test]
fn original_resistance_heat_and_port_work_match_independent_discrete_effort_trajectory() {
    let (mass, stiffness, cubic, damping, dt, input) = (1.0, 4.0, 80.0, 0.2, 0.1, 0.15);
    let initial = [0.3, 0.2];
    let mut runtime = PreparedNonlinearZoh::new_owned(
        duffing_oscillator(mass, stiffness, cubic, damping).unwrap(),
        dt,
        budget(2),
    )
    .unwrap();
    let step = runtime.step(&initial, &[input]).unwrap();
    let dq = step.x[0] - initial[0];
    let dp = step.x[1] - initial[1];
    let mid = 0.5 * (step.x[0] + initial[0]);
    // Independent Gonzalez gradient of U(q)=cubic*q^4/4. The correction uses
    // the complete original state increment, including momentum.
    let correction = 0.25 * cubic * mid * dq.powi(3) / (dq * dq + dp * dp);
    let held = [cubic * mid.powi(3) + correction * dq, correction * dp];
    let reference = rk4([initial[0], initial[1], 0.0, 0.0, 0.0], dt, 4096, |x| {
        let e = [stiffness * x[0] + held[0], x[1] / mass + held[1]];
        [
            e[1],
            -e[0] - damping * e[1] + input,
            damping * e[1] * e[1],
            input * e[1],
            e[1],
        ]
    });
    close(step.x[0], reference[0], 1.0);
    close(step.x[1], reference[1], 1.0);
    close(step.dissipated, reference[2], 1.0);
    close(step.supplied, reference[3], 1.0);
    close(step.y[0], reference[4] / dt, 1.0);
    assert!(step.dissipated > dt * damping * step.y[0] * step.y[0]);
    close(step.delta_h + step.dissipated, step.supplied, 1.0);
}

#[test]
fn passive_port_feedback_closes_the_same_nonlinear_endpoint_and_energy_window() {
    let make = || duffing_oscillator(1.0, 4.0, 80.0, 0.2).unwrap();
    let mut runtime = PreparedNonlinearZoh::new_owned(make(), 0.025, budget(2)).unwrap();
    let mut x = [0.3, 0.2];
    let before = runtime.system().hamiltonian(&x);
    let mut input = [0.0];
    let mut next = [0.0; 2];
    let mut output = [0.0; 1];
    let (mut heat, mut boundary_loss) = (0.0, 0.0);
    let mut calls = 0;
    for _ in 0..32 {
        let step = runtime
            .step_with_port_into(&x, &mut input, 0, &mut next, &mut output, |free, beta| {
                calls += 1;
                Ok(-free / (2.0 + beta))
            })
            .unwrap();
        close(output[0], -2.0 * input[0], 1.0);
        close(step.supplied, input[0] * output[0] * 0.025, 1.0);
        heat += step.dissipated;
        boundary_loss += 2.0 * input[0] * input[0] * 0.025;
        x = next;
    }
    assert!(
        calls > 32,
        "nonlinear trials must update feedback, not reuse a lagged pressure"
    );
    assert!(heat > 0.0 && boundary_loss > 0.0);
    close(
        runtime.system().hamiltonian(&x) - before + heat + boundary_loss,
        0.0,
        before,
    );
}

#[test]
fn arbitrary_energy_datum_does_not_change_motion_heat_or_work() {
    let mut original =
        PreparedNonlinearZoh::new_owned(body([7.0, 3.0], 0.0), 0.003, budget(5)).unwrap();
    let mut shifted =
        PreparedNonlinearZoh::new_owned(body([7.0, 3.0], 1e24), 0.003, budget(5)).unwrap();
    let initial = [0.3, -0.15, -0.2, 0.1, 0.65];
    let input = [0.23, -0.17];
    let a = original.step(&initial, &input).unwrap();
    let b = shifted.step(&initial, &input).unwrap();
    assert_eq!(a.x, b.x);
    assert_eq!(a.y, b.y);
    assert_eq!(a.delta_h, b.delta_h);
    assert_eq!(a.dissipated, b.dissipated);
    assert_eq!(a.supplied, b.supplied);
}

#[test]
fn unresolved_or_cancelled_candidates_publish_no_state_output_or_feedback() {
    let mut runtime = PreparedNonlinearZoh::new_owned(
        duffing_oscillator(1.0, 4.0, 80.0, 0.2).unwrap(),
        0.1,
        budget(2),
    )
    .unwrap();
    let initial = [0.3, 0.2];
    let mut next = [17.0, 19.0];
    let mut output = [23.0];
    runtime.set_max_iterations(1).unwrap();
    assert!(
        runtime
            .step_into(&initial, &[0.15], &mut next, &mut output)
            .is_err()
    );
    assert_eq!(next, [17.0, 19.0]);
    assert_eq!(output, [23.0]);
    runtime.set_max_iterations(50).unwrap();
    assert!(
        runtime
            .step_into_controlled(&initial, &[0.15], &mut next, &mut output, || true)
            .is_err()
    );
    assert_eq!(next, [17.0, 19.0]);
    assert_eq!(output, [23.0]);
    let mut input = [29.0];
    assert!(
        runtime
            .step_with_port_into(&initial, &mut input, 0, &mut next, &mut output, |_, _| Err(
                fs_phs::PhsError::LinearSplitBudget {
                    what: "deliberately refused physical load"
                }
            ))
            .is_err()
    );
    assert_eq!(input, [29.0]);
    assert_eq!(next, [17.0, 19.0]);
    assert_eq!(output, [23.0]);
    let accepted = runtime
        .step_into(&initial, &[0.15], &mut next, &mut output)
        .unwrap();
    assert!(next != [17.0, 19.0] && output != [23.0]);
    assert!(accepted.dissipated > 0.0);
}
