//! G1/G3: positive, independently integrated physical heat at very small motion.
use fs_phs::linear_exact::{
    LinearSplitBudget, PreparedExactZoh, integrated_resistance_loss, split_quadratic, step_render,
};
use fs_phs::{PortHamiltonian, QuadraticStorage, RelaxationBranch, mass_spring_damper};

fn budget(n: usize) -> LinearSplitBudget {
    LinearSplitBudget {
        max_modes: n,
        max_omega_dt: 0.04,
    }
}

fn relative(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        actual.is_finite() && actual >= 0.0,
        "physical heat must be finite and nonnegative: {actual:e}"
    );
    assert!(
        (actual / expected - 1.0).abs() <= tolerance,
        "physical heat {actual:.17e} vs {expected:.17e}"
    );
}

fn rk4<const N: usize>(
    mut x: [f64; N],
    duration: f64,
    steps: usize,
    rhs: impl Fn(&[f64; N]) -> [f64; N],
) -> [f64; N] {
    let dt = duration / steps as f64;
    for _ in 0..steps {
        let a = rhs(&x);
        let b = rhs(&core::array::from_fn(|i| x[i] + 0.5 * dt * a[i]));
        let c = rhs(&core::array::from_fn(|i| x[i] + 0.5 * dt * b[i]));
        let d = rhs(&core::array::from_fn(|i| x[i] + dt * c[i]));
        for i in 0..N {
            x[i] += dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    x
}

#[test]
fn held_equilibrium_is_silent_and_tiny_vibration_heat_remains_resolved() {
    for damping in [1e-200, 1e-100, 1e-8, 0.02, 2.0, 20.0] {
        let sys = mass_spring_damper(1.0, 1.0, damping).unwrap();
        let split = split_quadratic(&sys, budget(2)).unwrap();
        let prepared = PreparedExactZoh::new(&sys, 1.0, budget(2)).unwrap();
        assert_eq!(prepared.step(&[1.0, 0.0], &[1.0]).unwrap().dissipated, 0.0);
        assert_eq!(
            step_render(&sys, &split, &[1.0, 0.0], &[1.0], 1.0)
                .unwrap()
                .dissipated,
            0.0
        );
        // Integrate the independently written physical perturbation equation.
        // Its true q perturbation is the represented subtraction, not 1e-10.
        let q = 1.0 + 1e-10;
        let velocity = 0.0;
        let reference = rk4([q - 1.0, velocity, 0.0], 1.0, 4096, |x| {
            [x[1], -x[0] - damping * x[1], damping * x[1] * x[1]]
        });
        let direct = step_render(&sys, &split, &[q, velocity], &[1.0], 1.0).unwrap();
        let cached = prepared.step(&[q, velocity], &[1.0]).unwrap();
        relative(direct.dissipated, reference[2], 3e-11);
        assert_eq!(direct.dissipated.to_bits(), cached.dissipated.to_bits());
    }
}

#[test]
fn high_quality_factor_and_stiff_relaxation_keep_independent_positive_heat() {
    // One full lossless-period limit: integral sin(t)^2 dt = pi. The energy
    // subtraction is zero at this tiny damping, but the physical heat is not.
    for damping in [1e-24, 1e-100, 1e-200] {
        let sys = mass_spring_damper(1.0, 1.0, damping).unwrap();
        let split = split_quadratic(&sys, budget(2)).unwrap();
        let step = step_render(&sys, &split, &[0.1, 0.0], &[0.0], core::f64::consts::TAU).unwrap();
        relative(
            step.dissipated,
            damping * 0.01 * core::f64::consts::PI,
            2e-12,
        );
    }
    for rate in [1e-12, 0.25, 1.0, 100.0, 1e6] {
        let sys = PortHamiltonian::new(
            1,
            1,
            vec![0.0],
            vec![rate],
            vec![0.0],
            Box::new(QuadraticStorage::new(vec![1.0], 1).unwrap()),
        )
        .unwrap();
        let split = split_quadratic(&sys, budget(1)).unwrap();
        let prepared = PreparedExactZoh::new(&sys, 1.0, budget(1)).unwrap();
        let expected = -0.5 * (-2.0 * rate).exp_m1();
        relative(
            prepared.step(&[1.0], &[0.0]).unwrap().dissipated,
            expected,
            8e-12,
        );
        relative(
            step_render(&sys, &split, &[1.0], &[0.0], 1.0)
                .unwrap()
                .dissipated,
            expected,
            8e-12,
        );
        relative(
            integrated_resistance_loss(&sys, &split, &[1.0], &[0.0], 1.0, &[0.25 * rate]).unwrap(),
            0.25 * expected,
            8e-12,
        );
    }
}

fn relative_damper() -> PortHamiltonian {
    let mut j = vec![0.0; 16];
    j[1] = 1.0;
    j[4] = -1.0;
    j[11] = 1.0;
    j[14] = -1.0;
    let mut r = vec![0.0; 16];
    r[5] = 0.3;
    r[7] = -0.3;
    r[13] = -0.3;
    r[15] = 0.3;
    let q = vec![
        4.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    ];
    PortHamiltonian::new(
        4,
        1,
        j,
        r,
        vec![0.0, 1.0, 0.0, 1.0],
        Box::new(QuadraticStorage::new(q, 4).unwrap()),
    )
    .unwrap()
}

#[test]
fn a_damping_nullspace_and_its_small_relative_motion_are_not_replaced_by_heat_noise() {
    let sys = relative_damper();
    let prepared = PreparedExactZoh::new(&sys, 0.3, budget(4)).unwrap();
    let initial = [0.4, -0.3, 0.4, -0.3];
    let common = prepared.step(&initial, &[0.7]).unwrap();
    assert!(common.dissipated >= 0.0 && common.dissipated < 1e-28);
    let mut moving = initial;
    moving[3] += 1e-10;
    let reference = rk4([0.0, moving[1] - moving[3], 0.0], 0.3, 2048, |x| {
        [x[1], -4.0 * x[0] - 0.6 * x[1], 0.3 * x[1] * x[1]]
    });
    relative(
        prepared.step(&moving, &[0.7]).unwrap().dissipated,
        reference[2],
        2e-6,
    );
    let lossless = mass_spring_damper(1.0, 4.0, 0.0).unwrap();
    assert_eq!(
        PreparedExactZoh::new(&lossless, 0.3, budget(2))
            .unwrap()
            .step(&[0.4, -0.3], &[0.7])
            .unwrap()
            .dissipated,
        0.0
    );
}

fn material_body() -> PortHamiltonian {
    let mut j = vec![0.0; 16];
    j[1] = 1.0;
    j[4] = -1.0;
    j[11] = 1.0;
    j[14] = -1.0;
    let mut r = vec![0.0; 16];
    r[5] = 0.4;
    r[7] = -0.1;
    r[13] = -0.1;
    r[15] = 0.2;
    let q = vec![
        4.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 9.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    ];
    PortHamiltonian::new(
        4,
        1,
        j,
        r,
        vec![0.0, 1.0, 0.0, -0.3],
        Box::new(QuadraticStorage::new(q, 4).unwrap()),
    )
    .unwrap()
    .with_relaxation_branches(vec![RelaxationBranch {
        projection: vec![1.0, 0.0, -0.5, 0.0],
        stiffness: 4.0,
        relaxation_time_s: 0.15,
    }])
    .unwrap()
}

#[test]
fn complete_material_heat_matches_physical_rhs_and_clock_subdivision_with_atomic_outputs() {
    let sys = material_body();
    let initial = [0.02, -0.03, 0.01, 0.04, 0.003];
    let force = 0.7;
    let duration = 0.3;
    let reference = rk4(
        [
            initial[0], initial[1], initial[2], initial[3], initial[4], 0.0, 0.0,
        ],
        duration,
        4096,
        |x| {
            let arm = 2.0 * x[0] - x[2] - x[4];
            [
                x[1],
                -4.0 * x[0] - 0.4 * x[1] + 0.1 * x[3] - 2.0 * arm + force,
                x[3],
                -9.0 * x[2] + 0.1 * x[1] - 0.2 * x[3] + arm - 0.3 * force,
                arm / 0.15,
                0.4 * x[1] * x[1] - 0.2 * x[1] * x[3] + 0.2 * x[3] * x[3] + arm * arm / 0.15,
                force * (x[1] - 0.3 * x[3]),
            ]
        },
    );
    let mut full = PreparedExactZoh::new(&sys, duration, budget(5)).unwrap();
    let record = full.step(&initial, &[force]).unwrap();
    for (actual, expected) in record.x.iter().zip(reference) {
        assert!((actual - expected).abs() < 2e-12, "{actual} vs {expected}");
    }
    relative(record.dissipated, reference[5], 3e-11);
    assert!((record.supplied - reference[6]).abs() < 2e-12);
    for count in [2, 4] {
        let sub = PreparedExactZoh::new(&sys, duration / count as f64, budget(5)).unwrap();
        let mut x = initial.to_vec();
        let mut heat = 0.0;
        for _ in 0..count {
            let r = sub.step(&x, &[force]).unwrap();
            x = r.x;
            heat += r.dissipated;
        }
        relative(heat, record.dissipated, 3e-12);
        for (a, b) in x.iter().zip(&record.x) {
            assert!((a - b).abs() < 2e-12);
        }
    }
    let mut next = vec![123.0; 5];
    let mut output = vec![456.0];
    assert!(
        full.step_into(&initial, &[f64::NAN], &mut next, &mut output)
            .is_err()
    );
    assert_eq!(next, vec![123.0; 5]);
    assert_eq!(output, vec![456.0]);
    let accepted = full
        .step_into(&initial, &[force], &mut next, &mut output)
        .unwrap();
    assert_eq!(next, record.x);
    assert_eq!(accepted.dissipated, record.dissipated);
}
