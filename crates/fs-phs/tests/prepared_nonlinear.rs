//! G0/G1/G3: full nonlinear storage with exact linear and material-memory flow.
//! Heat uses the original resistance on the method's discrete effort; only the
//! quadratic limit is an exact continuous physical trajectory.
use fs_phs::linear_exact::{
    DenseExponentialBudget, LinearSplitBudget, PreparedExactTones, PreparedExactZoh,
    PreparedNonlinearTones, PreparedNonlinearZoh, QuadraticPart, split_quadratic,
};
use fs_phs::{
    PortHamiltonian, QuadraticStorage, RelaxationBranch, ResistanceWitness, Storage,
    duffing_oscillator,
};

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


#[test]
fn complete_reference_boundary_retains_every_mode_port_and_zero_energy_coordinate() {
    // The actual propagation matrix reaches 256 coordinates with one output.
    // This avoids preparing 51 separate output integrals just to test the same
    // matrix ceiling; the coupled KC regression covers its many-port layout.
    let modes = 63;
    let n = 2 * modes + 1;
    let dt = 0.01;
    let input = 0.2;
    let mut j = vec![0.0; n * n];
    let r = vec![0.0; n * n];
    let mut g = vec![0.0; n];
    let mut q = vec![0.0; n * n];
    let mut initial = vec![0.0; n];
    for mode in 0..modes {
        let position = 2 * mode;
        let momentum = position + 1;
        let omega = 3.0 + 0.1 * mode as f64;
        j[position * n + momentum] = 1.0;
        j[momentum * n + position] = -1.0;
        q[position * n + position] = omega * omega;
        q[momentum * n + momentum] = 1.0;
        g[momentum] = 0.1 + 0.01 * mode as f64;
        initial[position] = if mode % 2 == 0 { 0.01 } else { -0.01 };
        initial[momentum] = 0.003 * ((mode % 3) as f64 - 1.0);
    }
    // This original zero-Q coordinate still has a nonzero held drive. It
    // carries no stored energy, but its state must not disappear in preparation.
    initial[n - 1] = 0.19;
    g[n - 1] = 0.2;
    let system = PortHamiltonian::new(
        n,
        1,
        j.clone(),
        r.clone(),
        g.clone(),
        Box::new(QuadraticStorage::new(q, n).unwrap()),
    )
    .unwrap();
    assert!(PreparedNonlinearZoh::fits_reference_capacity(n, 1));
    // The ordinary split continues to reserve its legacy doubled ZOH block.
    assert!(split_quadratic(&system, budget(n)).is_err());
    let before = system.hamiltonian(&initial);
    let mut runtime = PreparedNonlinearZoh::new_owned(system, dt, budget(n)).unwrap();
    assert_eq!(runtime.system().state_dim(), n);
    assert_eq!(runtime.system().port_dim(), 1);
    assert_eq!(
        runtime.system().structure(),
        (j.as_slice(), r.as_slice(), g.as_slice())
    );
    let step = runtime.step(&initial, &[input]).unwrap();
    let mut reference_work = 0.0;
    for mode in 0..modes {
        let position = 2 * mode;
        let momentum = position + 1;
        let omega = 3.0 + 0.1 * mode as f64;
        let (sin, cos) = (omega * dt).sin_cos();
        let equilibrium = g[momentum] * input / (omega * omega);
        let relative = initial[position] - equilibrium;
        let expected_q = equilibrium + relative * cos + initial[momentum] * sin / omega;
        let expected_p = initial[momentum] * cos - omega * relative * sin;
        close(step.x[position], expected_q, 0.1);
        close(step.x[momentum], expected_p, 0.1);
        reference_work += input * g[momentum] * (expected_q - initial[position]);
    }
    close(step.x[n - 1], initial[n - 1] + dt * g[n - 1] * input, 1.0);
    close(step.y[0], reference_work / (input * dt), 1.0);
    close(step.supplied, reference_work, before);
    close(step.delta_h, runtime.system().hamiltonian(&step.x) - before, before);
    close(step.delta_h, reference_work, before);
    assert_eq!(step.dissipated, 0.0);
    assert!((step.x[n - 3] - initial[n - 3]).abs() > 1e-6);
}

#[test]
fn reference_capacity_is_overflow_safe_and_refuses_before_large_allocation() {
    // Complete bare KC layouts use two states and one force port per mode.
    for modes in [25, 26, 43, 51] {
        assert!(PreparedNonlinearZoh::fits_reference_capacity(2 * modes, modes));
    }
    assert!(!PreparedNonlinearZoh::fits_reference_capacity(104, 52));
    assert!(PreparedNonlinearZoh::fits_reference_capacity(127, 1));
    for (n, m) in [(0, 0), (127, 2), (128, 0), (usize::MAX, 0), (1, usize::MAX)] {
        assert!(!PreparedNonlinearZoh::fits_reference_capacity(n, m));
        // Deliberately omit all large matrices: capacity must refuse before
        // shape access, storage probing or any dimension-dependent allocation.
        let system = PortHamiltonian::from_raw_parts(
            n,
            m,
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Box::new(QuadraticStorage::new(vec![1.0], 1).unwrap()),
        );
        assert!(PreparedNonlinearZoh::new_owned(system, 0.01, budget(usize::MAX)).is_err());
    }
    let make = || duffing_oscillator(1.0, 4.0, 3.0, 0.2).unwrap();
    assert!(PreparedNonlinearZoh::fits_reference_capacity(2, 1));
    assert!(PreparedNonlinearZoh::new_owned(make(), 0.01, budget(1)).is_err());
    assert!(PreparedNonlinearZoh::new_owned(make(), 0.0, budget(2)).is_err());
    assert!(PreparedNonlinearZoh::new_owned(make(), f64::NAN, budget(2)).is_err());
    assert!(
        PreparedNonlinearZoh::new_owned(
            make(),
            0.01,
            LinearSplitBudget { max_modes: 2, max_omega_dt: f64::NAN },
        )
        .is_err()
    );
}

// Preparation must use the explicit decomposition, never a sampled guess for H.
struct PreparationOnlyStorage(Option<Vec<f64>>);
impl Storage for PreparationOnlyStorage {
    fn hamiltonian(&self, _x: &[f64]) -> f64 {
        panic!("preparation must not probe the storage value");
    }
    fn gradient(&self, _x: &[f64], _out: &mut [f64]) {
        panic!("preparation must not probe the storage gradient");
    }
    fn quadratic_part(&self) -> Option<QuadraticPart> {
        self.0.as_ref().map(|matrix| QuadraticPart { matrix: matrix.clone(), exact: true })
    }
}

#[test]
fn prepared_capacity_path_keeps_raw_structure_and_storage_admission() {
    let j = vec![0.0, 1.0, -1.0, 0.0];
    let r = vec![0.0; 4];
    let g = vec![0.0, 1.0];
    let q = Some(vec![4.0, 0.0, 0.0, 1.0]);
    let invalid = [
        ("J shape", vec![0.0; 3], r.clone(), g.clone(), q.clone()),
        ("J finite", vec![0.0, f64::NAN, -1.0, 0.0], r.clone(), g.clone(), q.clone()),
        ("J skew", vec![1.0, 1.0, -1.0, 0.0], r.clone(), g.clone(), q.clone()),
        ("R symmetry", j.clone(), vec![0.0, 1.0, 0.0, 0.0], g.clone(), q.clone()),
        ("R PSD", j.clone(), vec![-1.0, 0.0, 0.0, 0.0], g.clone(), q.clone()),
        ("G shape", j.clone(), r.clone(), vec![0.0], q.clone()),
        ("G finite", j.clone(), r.clone(), vec![0.0, f64::NAN], q.clone()),
        ("Q absent", j.clone(), r.clone(), g.clone(), None),
        ("Q shape", j.clone(), r.clone(), g.clone(), Some(vec![1.0; 3])),
        ("Q finite", j.clone(), r.clone(), g.clone(), Some(vec![f64::NAN, 0.0, 0.0, 1.0])),
        ("Q symmetry", j.clone(), r.clone(), g.clone(), Some(vec![4.0, 0.1, 0.0, 1.0])),
        ("Q PSD", j, r, g, Some(vec![-1.0, 0.0, 0.0, 1.0])),
    ];
    for (label, j, r, g, q) in invalid {
        let system = PortHamiltonian::from_raw_parts(
            2, 1, j, r, g, Box::new(PreparationOnlyStorage(q)),
        );
        assert!(
            PreparedNonlinearZoh::new_owned(system, 0.01, budget(2)).is_err(),
            "{label} must retain its original admission",
        );
    }
}

fn tone_budget() -> DenseExponentialBudget {
    DenseExponentialBudget {
        max_dimension: 256, max_dense_bytes: usize::MAX, max_multiply_terms: u64::MAX,
    }
}

struct PureQuartic;
impl Storage for PureQuartic {
    fn hamiltonian(&self, x: &[f64]) -> f64 { 0.25 * x[0].powi(4) }
    fn gradient(&self, x: &[f64], out: &mut [f64]) { out[0] = x[0].powi(3); }
    fn quadratic_part(&self) -> Option<QuadraticPart> {
        Some(QuadraticPart { matrix: vec![0.0], exact: false })
    }
    fn midpoint_energy_remainder(&self, a: &[f64], b: &[f64], mid: &[f64]) -> Option<f64> {
        let secant = 0.25 * (a[0] + b[0]) * (a[0] * a[0] + b[0] * b[0]);
        Some((b[0] - a[0]) * (secant - mid[0].powi(3)))
    }
}

#[test]
fn nonlinear_tones_integrate_original_quartic_work_and_heat_analytically() {
    let dt = 0.08;
    let resistance = 0.4;
    let original = PortHamiltonian::new(
        1, 1, vec![0.0], vec![resistance], vec![1.0], Box::new(PureQuartic),
    ).unwrap();
    let frequencies = [0.0, 3.0, 3.0, 3.0, 7.1];
    let mut owner = PreparedNonlinearTones::new_owned(
        original, &frequencies, dt, budget(1), tone_budget(),
    ).unwrap();
    let large = 2.0_f64.powi(54);
    let cosine = [0.2, large, 0.7, -large, -0.3];
    let sine = [9.0, 0.1, 0.2, 0.3, 0.4]; // The DC sine is identically zero.
    let initial = [0.7];
    let step = owner.step(&initial, &cosine, &sine).unwrap();
    let end = step.x[0];
    let f = 0.25 * (initial[0] + end) * (initial[0].powi(2) + end * end);
    let integral = 0.2 * dt
        + 0.7 * (3.0 * dt).sin() / 3.0 + 0.6 * (1.0 - (3.0 * dt).cos()) / 3.0
        - 0.3 * (7.1 * dt).sin() / 7.1 + 0.4 * (1.0 - (7.1 * dt).cos()) / 7.1;
    // With Q=0 the accepted discrete effort is constant, even though the
    // original storage is quartic and the applied source varies continuously.
    close(end, initial[0] - dt * resistance * f + integral, 1.0);
    close(step.y[0], f, 1.0);
    close(step.dissipated, dt * resistance * f * f, 1.0);
    close(step.supplied, f * integral, 1.0);
    close(step.delta_h, 0.25 * (end.powi(4) - initial[0].powi(4)), 1.0);
    assert!(step.dissipated > 0.0 && step.newton_iters > 0);
    assert_eq!(owner.angular_frequencies_rad_s(), &frequencies);
    assert_eq!((owner.state_dim(), owner.port_dim(), owner.tone_count()), (1, 1, 5));
    let (c, s) = owner.advance_force_phase(4, cosine[4], sine[4]).unwrap();
    close(c, -0.3 * (7.1 * dt).cos() + 0.4 * (7.1 * dt).sin(), 1.0);
    close(s, 0.4 * (7.1 * dt).cos() + 0.3 * (7.1 * dt).sin(), 1.0);
}

#[test]
fn nonlinear_tone_dc_and_quadratic_limits_reuse_the_original_physics() {
    let make = || duffing_oscillator(1.0, 4.0, 7.0, 0.2).unwrap();
    let initial = [0.31, -0.11];
    let mut held = PreparedNonlinearZoh::new_owned(make(), 0.02, budget(2)).unwrap();
    let mut dc = PreparedNonlinearTones::new_owned(
        make(), &[0.0], 0.02, budget(2), tone_budget(),
    ).unwrap();
    let reference = held.step(&initial, &[0.4]).unwrap();
    let actual = dc.step(&initial, &[0.4], &[100.0]).unwrap();
    for (a, b) in actual.x.iter().zip(&reference.x) { close(*a, *b, 1.0); }
    close(actual.y[0], reference.y[0], 1.0);
    close(actual.supplied, reference.supplied, 1.0);
    close(actual.dissipated, reference.dissipated, 1.0);

    let make = || duffing_oscillator(1.0, 4.0, 0.0, 0.2).unwrap();
    let frequencies = [0.0, 3.0, 3.00001, 7.0];
    let mut exact = PreparedExactTones::new_owned(
        make(), &frequencies, 0.1, budget(2), tone_budget(),
    ).unwrap();
    let mut nonlinear = PreparedNonlinearTones::new_owned(
        make(), &frequencies, 0.1, budget(2), tone_budget(),
    ).unwrap();
    let (cosine, sine) = ([0.2, 0.8, -0.4, 0.3], [0.0, -0.1, 0.2, 0.4]);
    let (mut x, mut y) = ([0.0; 2], [0.0]);
    let reference = exact.step_into(&initial, &cosine, &sine, &mut x, &mut y).unwrap();
    let actual = nonlinear.step(&initial, &cosine, &sine).unwrap();
    for (a, b) in actual.x.iter().zip(&x) { close(*a, *b, 1.0); }
    close(actual.y[0], y[0], 1.0);
    close(actual.supplied, reference.supplied, 1.0);
    close(actual.dissipated, reference.dissipated, 1.0);
    close(actual.delta_h, reference.delta_h, 1.0);
    assert_eq!(actual.newton_iters, 0);
}

#[test]
fn continuous_tones_refine_original_nonlinear_motion_work_and_heat() {
    let initial = [0.7, -0.2];
    let duration = 0.4;
    // Independent continuous Duffing equations, with original heat and work
    // integrated alongside the state. Time is the fifth reference coordinate.
    let reference = rk4([initial[0], initial[1], 0.0, 0.0, 0.0], duration, 8192, |x| {
        let t = x[4];
        let u = 0.8 * (2.3 * t).cos() + 0.3 * (2.3 * t).sin()
            - 0.25 * (4.1 * t).cos() + 0.2 * (4.1 * t).sin();
        [x[1], -4.0 * x[0] - 7.0 * x[0].powi(3) - 0.2 * x[1] + u,
            0.2 * x[1] * x[1], u * x[1], 1.0]
    });
    let mut errors = Vec::new();
    for steps in [8, 16, 32] {
        let mut owner = PreparedNonlinearTones::new_owned(
            duffing_oscillator(1.0, 4.0, 7.0, 0.2).unwrap(), &[2.3, 4.1],
            duration / f64::from(steps), budget(2), tone_budget(),
        ).unwrap();
        let (mut cosine, mut sine) = ([0.8, -0.25], [0.3, 0.2]);
        let (mut x, mut next, mut y) = (initial, [0.0; 2], [0.0]);
        let (mut heat, mut work) = (0.0, 0.0);
        for _ in 0..steps {
            let record = owner.step_into(&x, &cosine, &sine, &mut next, &mut y).unwrap();
            heat += record.dissipated;
            work += record.supplied;
            x = next;
            for tone in 0..2 {
                (cosine[tone], sine[tone]) = owner.advance_force_phase(
                    tone, cosine[tone], sine[tone],
                ).unwrap();
            }
        }
        close(owner.system().hamiltonian(&x) - owner.system().hamiltonian(&initial) + heat,
            work, 1.0);
        errors.push([x[0] - reference[0], x[1] - reference[1],
            heat - reference[2], work - reference[3]].map(f64::abs));
    }
    for component in 0..4 {
        assert!(errors[1][component] < 0.35 * errors[0][component], "{errors:?}");
        assert!(errors[2][component] < 0.35 * errors[1][component], "{errors:?}");
    }
}

#[test]
fn nonlinear_tone_capacity_covers_work_and_refusals_are_transactional() {
    assert!(PreparedNonlinearTones::fits_reference_capacity(126, 1, 1));
    assert!(PreparedNonlinearTones::fits_reference_capacity(1, 85, 1));
    for (n, m, k) in [(127, 1, 1), (1, 86, 1), (0, 1, 1), (1, 0, 1),
        (1, 1, 0), (usize::MAX, 1, 1), (1, usize::MAX, 1), (1, 1, usize::MAX)]
    {
        assert!(!PreparedNonlinearTones::fits_reference_capacity(n, m, k));
    }
    let make = || duffing_oscillator(1.0, 4.0, 7.0, 0.2).unwrap();
    for denied in [
        DenseExponentialBudget { max_dimension: 6, ..tone_budget() },
        DenseExponentialBudget { max_dense_bytes: 0, ..tone_budget() },
        DenseExponentialBudget { max_multiply_terms: 0, ..tone_budget() },
    ] {
        assert!(PreparedNonlinearTones::new_owned(make(), &[3.0], 0.1, budget(2), denied).is_err());
    }
    // Repeated authored terms occupy one rotating pair, not 80 pairs.
    assert!(PreparedNonlinearTones::new_owned(
        make(), &[3.0; 80], 0.1, budget(2),
        DenseExponentialBudget { max_dimension: 7, ..tone_budget() },
    ).is_ok());
    let mut owner = PreparedNonlinearTones::new_owned(
        make(), &[3.0], 0.1, budget(2), tone_budget(),
    ).unwrap();
    let (mut next, mut output) = ([41.0; 2], [42.0]);
    assert!(owner.step_into(&[0.8, 0.3], &[f64::NAN], &[0.2], &mut next, &mut output).is_err());
    let mut polls = 0;
    assert!(matches!(owner.step_into_controlled(
        &[0.8, 0.3], &[0.7], &[0.2], &mut next, &mut output,
        || { polls += 1; polls >= 6 },
    ), Err(fs_phs::PreparedStepError::Cancelled)));
    owner.set_max_iterations(0).unwrap();
    assert!(owner.step_into(&[0.8, 0.3], &[0.7], &[0.2], &mut next, &mut output).is_err());
    assert_eq!(next, [41.0; 2]);
    assert_eq!(output, [42.0]);
    owner.set_max_iterations(50).unwrap();
    owner.step_into(&[0.8, 0.3], &[0.7], &[0.2], &mut next, &mut output).unwrap();
    assert!(next.iter().all(|v| v.is_finite()) && output[0].is_finite());
}

fn affine_force_integral(omega: f64, c0: f64, s0: f64, c1: f64, s1: f64, h: f64) -> f64 {
    if omega == 0.0 { return 0.5 * h * (c0 + c1); }
    let (sin, cos) = (omega * h).sin_cos();
    c0 * sin / omega + s0 * (1.0 - cos) / omega
        + (c1 - c0) * (sin / omega + (cos - 1.0) / (h * omega * omega))
        + (s1 - s0) * (-cos / omega + sin / (h * omega * omega))
}

#[test]
fn affine_nonlinear_tones_integrate_original_quartic_source_and_channel_heat() {
    let (dt, resistance) = (0.08, 0.4);
    let original = PortHamiltonian::new(
        1, 1, vec![0.0], vec![resistance], vec![1.0], Box::new(PureQuartic),
    ).unwrap();
    // The last observation overlaps the first two and exceeds total R. Its
    // account must stay independent; these are not assumed to partition R.
    let losses = [0.1_f64, 0.3, 0.8];
    let channels: Vec<_> = losses.iter().map(|&r| {
        ResistanceWitness::try_from_source_rows(1, vec![r], vec![r.sqrt()]).unwrap()
    }).collect();
    let frequencies = [0.0, 3.0, 3.0, 3.0, 7.1];
    let mut owner = PreparedNonlinearTones::new_owned_affine_with_resistance_channels(
        original, &frequencies, dt, budget(1), tone_budget(), &channels,
    ).unwrap();
    let large = 2.0_f64.powi(54);
    let (c0, s0) = ([0.2, large, 0.7, -large, -0.3], [9.0, 0.1, 0.2, 0.3, 0.4]);
    let (c1, s1) = ([-0.1, large, 0.9, -large, 0.2], [-14.0, 0.2, -0.1, 0.1, -0.25]);
    let initial = [0.7];
    let (mut next, mut mean, mut heat) = ([0.0], [0.0], [0.0; 3]);
    let record = owner.step_affine_observed_into(
        &initial, &c0, &s0, &c1, &s1, &mut next, &mut mean, &mut heat,
    ).unwrap();
    let f = 0.25 * (initial[0] + next[0]) * (initial[0] * initial[0] + next[0] * next[0]);
    let integral = affine_force_integral(0.0, 0.2, 0.0, -0.1, 0.0, dt)
        + affine_force_integral(3.0, 0.7, 0.6, 0.9, 0.2, dt)
        + affine_force_integral(7.1, -0.3, 0.4, 0.2, -0.25, dt);
    // Q=0 leaves a constant accepted discrete effort. Integrals of the
    // prescribed affine sinusoids give an independent complete energy oracle.
    close(next[0], initial[0] - dt * resistance * f + integral, 1.0);
    close(mean[0], f, 1.0);
    close(record.supplied, f * integral, 1.0);
    close(record.dissipated, resistance * dt * f * f, 1.0);
    close(record.delta_h, 0.25 * (next[0].powi(4) - initial[0].powi(4)), 1.0);
    for (&actual, &r) in heat.iter().zip(&losses) { close(actual, r * dt * f * f, 1.0); }
    assert!(owner.is_affine() && record.newton_iters > 0);
    assert_eq!(owner.resistance_channel_count(), 3);
    assert!(heat[2] > record.dissipated);
}

#[test]
fn nonlinear_affine_and_channel_limits_preserve_complete_quadratic_memory() {
    let n = 5;
    let mut row = vec![0.0; n];
    row[1] = 0.3;
    row[3] = -0.2;
    let mut matrix = vec![0.0; n * n];
    for i in 0..n { for j in 0..n { matrix[i * n + j] = row[i] * row[j]; } }
    let channels = [ResistanceWitness::try_from_source_rows(n, matrix, row).unwrap()];
    let frequencies = [0.0, 2.7, 3.01];
    let mut exact = PreparedExactTones::new_owned_affine_with_resistance_channels(
        body([0.0; 2], 0.0), &frequencies, 0.003, budget(n), tone_budget(), &channels,
    ).unwrap();
    let mut nonlinear = PreparedNonlinearTones::new_owned_affine_with_resistance_channels(
        body([0.0; 2], 0.0), &frequencies, 0.003, budget(n), tone_budget(), &channels,
    ).unwrap();
    let initial = [0.07, -0.03, -0.04, 0.02, 0.025];
    let (c0, s0) = ([0.1, 0.3, -0.2, 0.2, -0.1, 0.4], [0.0, 0.2, 0.3, 0.0, 0.1, -0.2]);
    let (c1, s1) = ([0.2, -0.4, 0.1, -0.3, 0.2, 0.1], [0.0, -0.1, 0.1, 0.0, 0.3, 0.2]);
    let (mut x, mut y, mut heat) = ([0.0; 5], [0.0; 2], [0.0]);
    let expected = exact.step_affine_observed_into(
        &initial, &c0, &s0, &c1, &s1, &mut x, &mut y, &mut heat,
    ).unwrap();
    let (mut next, mut mean, mut loss) = ([0.0; 5], [0.0; 2], [0.0]);
    let actual = nonlinear.step_affine_observed_into(
        &initial, &c0, &s0, &c1, &s1, &mut next, &mut mean, &mut loss,
    ).unwrap();
    for (a, b) in next.iter().chain(&mean).chain(&loss).zip(x.iter().chain(&y).chain(&heat)) {
        close(*a, *b, 1.0);
    }
    close(actual.supplied, expected.supplied, 1.0);
    close(actual.dissipated, expected.dissipated, 1.0);
    close(actual.delta_h, expected.delta_h, 1.0);
    assert_eq!(actual.newton_iters, 0);

    let make = || duffing_oscillator(1.0, 4.0, 7.0, 0.2).unwrap();
    let channels = [ResistanceWitness::try_from_source_rows(
        2, vec![0.0, 0.0, 0.0, 0.09], vec![0.0, 0.3],
    ).unwrap()];
    let mut fixed = PreparedNonlinearTones::new_owned_with_resistance_channels(
        make(), &[2.7], 0.03, budget(2), tone_budget(), &channels,
    ).unwrap();
    let mut affine = PreparedNonlinearTones::new_owned_affine_with_resistance_channels(
        make(), &[2.7], 0.03, budget(2), tone_budget(), &channels,
    ).unwrap();
    let (mut a, mut ya, mut ha) = ([0.0; 2], [0.0], [0.0]);
    let (mut b, mut yb, mut hb) = ([0.0; 2], [0.0], [0.0]);
    let ra = fixed.step_observed_into(&[0.2, -0.1], &[0.6], &[-0.2], &mut a, &mut ya, &mut ha).unwrap();
    let rb = affine.step_affine_observed_into(
        &[0.2, -0.1], &[0.6], &[-0.2], &[0.6], &[-0.2], &mut b, &mut yb, &mut hb,
    ).unwrap();
    for (x, y) in a.iter().chain(&ya).chain(&ha).zip(b.iter().chain(&yb).chain(&hb)) { close(*x, *y, 1.0); }
    close(ra.supplied, rb.supplied, 1.0);
    close(ra.dissipated, rb.dissipated, 1.0);
    close(ha[0], 0.45 * ra.dissipated, 1.0);
}

#[test]
fn nonlinear_affine_sources_refine_state_work_and_channels_with_cached_phase() {
    let (initial, duration) = ([0.7, -0.2], 0.4);
    let reference = rk4([initial[0], initial[1], 0.0, 0.0, 0.0], duration, 8192, |x| {
        let t = x[4];
        let a = t / duration;
        let u = (0.6 + 0.9 * a) * (2.3 * t).cos() + (-0.2 + 0.3 * a) * (2.3 * t).sin()
            + (0.1 - 0.4 * a) * (4.1 * t).cos() + (0.3 - 0.2 * a) * (4.1 * t).sin();
        [x[1], -4.0 * x[0] - 7.0 * x[0].powi(3) - 0.2 * x[1] + u,
            0.2 * x[1] * x[1], u * x[1], 1.0]
    });
    let channel = ResistanceWitness::try_from_source_rows(
        2, vec![0.0, 0.0, 0.0, 0.09], vec![0.0, 0.3],
    ).unwrap();
    let mut errors = Vec::new();
    for steps in [8, 16, 32] {
        let mut owner = PreparedNonlinearTones::new_owned_affine_with_resistance_channels(
            duffing_oscillator(1.0, 4.0, 7.0, 0.2).unwrap(), &[2.3, 4.1],
            duration / f64::from(steps), budget(2), tone_budget(), core::slice::from_ref(&channel),
        ).unwrap();
        let (mut cosine, mut sine) = ([0.6, 0.1], [-0.2, 0.3]);
        let (mut dc, mut ds) = ([0.9, -0.4], [0.3, -0.2]);
        let (mut x, mut next, mut mean, mut loss) = (initial, [0.0; 2], [0.0], [0.0]);
        let (mut heat, mut work, mut observed) = (0.0, 0.0, 0.0);
        for _ in 0..steps {
            let c1 = core::array::from_fn::<_, 2, _>(|i| cosine[i] + dc[i] / f64::from(steps));
            let s1 = core::array::from_fn::<_, 2, _>(|i| sine[i] + ds[i] / f64::from(steps));
            let record = owner.step_affine_observed_into(
                &x, &cosine, &sine, &c1, &s1, &mut next, &mut mean, &mut loss,
            ).unwrap();
            x = next;
            heat += record.dissipated;
            work += record.supplied;
            observed += loss[0];
            for tone in 0..2 {
                (cosine[tone], sine[tone]) = owner.advance_force_phase(tone, c1[tone], s1[tone]).unwrap();
                (dc[tone], ds[tone]) = owner.advance_force_phase(tone, dc[tone], ds[tone]).unwrap();
            }
        }
        close(owner.system().hamiltonian(&x) - owner.system().hamiltonian(&initial) + heat, work, 1.0);
        close(observed, 0.45 * heat, 1.0);
        errors.push([x[0] - reference[0], x[1] - reference[1],
            heat - reference[2], work - reference[3], observed - 0.45 * reference[2]].map(f64::abs));
    }
    for component in 0..5 {
        assert!(errors[1][component] < 0.35 * errors[0][component], "{errors:?}");
        assert!(errors[2][component] < 0.35 * errors[1][component], "{errors:?}");
    }
}

#[test]
fn nonlinear_affine_channel_admission_and_publication_are_transactional() {
    assert!(PreparedNonlinearTones::fits_affine_reference_capacity(125, 1, 1));
    assert!(PreparedNonlinearTones::fits_affine_reference_capacity(1, 51, 1));
    for (n, m, k) in [(126, 1, 1), (1, 52, 1), (0, 1, 1), (1, 0, 1), (1, 1, 0),
        (usize::MAX, 1, 1), (1, usize::MAX, 1), (1, 1, usize::MAX)]
    {
        assert!(!PreparedNonlinearTones::fits_affine_reference_capacity(n, m, k));
    }
    let make = || duffing_oscillator(1.0, 4.0, 7.0, 0.2).unwrap();
    let channel = ResistanceWitness::try_from_source_rows(
        2, vec![0.0, 0.0, 0.0, 0.09], vec![0.0, 0.3],
    ).unwrap();
    // The old fixed-source dimension and a budget admitting no channel cannot
    // stand in for the actual four-quadrature/channel payload and cold work.
    for denied in [
        DenseExponentialBudget { max_dimension: 8, ..tone_budget() },
        DenseExponentialBudget { max_dense_bytes: 60_000, ..tone_budget() },
        DenseExponentialBudget { max_multiply_terms: 2_300_000, ..tone_budget() },
    ] {
        assert!(PreparedNonlinearTones::new_owned_affine_with_resistance_channels(
            make(), &[2.7], 0.03, budget(2), denied, core::slice::from_ref(&channel),
        ).is_err());
    }
    let wrong = ResistanceWitness::try_from_source_rows(1, vec![1.0], vec![1.0]).unwrap();
    assert!(PreparedNonlinearTones::new_owned_with_resistance_channels(
        make(), &[2.7], 0.03, budget(2), tone_budget(), &[wrong],
    ).is_err());
    assert!(PreparedNonlinearTones::new_owned_affine(
        make(), &[2.7], f64::from_bits(1), budget(2), tone_budget(),
    ).is_err());
    let mut owner = PreparedNonlinearTones::new_owned_affine_with_resistance_channels(
        make(), &[2.7], 0.03, budget(2), tone_budget(), &[channel],
    ).unwrap();
    let (mut next, mut mean, mut heat) = ([41.0; 2], [42.0], [43.0]);
    assert!(owner.step_into(&[0.7, -0.2], &[0.4], &[0.1], &mut next, &mut mean).is_err());
    assert!(owner.step_affine_observed_into(
        &[0.7, -0.2], &[0.4], &[0.1], &[f64::NAN], &[0.2], &mut next, &mut mean, &mut heat,
    ).is_err());
    assert!(owner.step_affine_observed_into(
        &[0.7, -0.2], &[0.4], &[0.1], &[0.7], &[0.2], &mut next, &mut mean, &mut [],
    ).is_err());
    let mut polls = 0;
    assert!(matches!(owner.step_affine_observed_into_controlled(
        &[0.7, -0.2], &[0.4], &[0.1], &[0.7], &[0.2], &mut next, &mut mean, &mut heat,
        || { polls += 1; polls >= 6 },
    ), Err(fs_phs::PreparedStepError::Cancelled)));
    owner.set_max_iterations(0).unwrap();
    assert!(owner.step_affine_observed_into(
        &[0.7, -0.2], &[0.4], &[0.1], &[0.7], &[0.2], &mut next, &mut mean, &mut heat,
    ).is_err());
    assert_eq!((next, mean, heat), ([41.0; 2], [42.0], [43.0]));
    owner.set_max_iterations(50).unwrap();
    owner.step_affine_observed_into(
        &[0.7, -0.2], &[0.4], &[0.1], &[0.7], &[0.2], &mut next, &mut mean, &mut heat,
    ).unwrap();
    assert!(heat[0] > 0.0);
}
