use super::*;
use crate::hb_residual;

fn coefficient(values: &[f64], pack: &Packing, state: usize, order: isize) -> C64 {
    let n = order.unsigned_abs();
    if n > pack.n_harm {
        return C64::ZERO;
    }
    let z = C64::new(
        values[pack.re_index(n, state)],
        if n == 0 {
            0.0
        } else {
            values[pack.im_index(n, state)]
        },
    );
    if order < 0 { z.conj() } else { z }
}

fn action<P: OrbitProblem>(problem: &P, pack: &Packing, u: &[f64], direction: &[f64]) -> Vec<f64> {
    let mut aft = AftWorkspace::new(pack).unwrap();
    let mut residual = vec![0.0; pack.coeff_len()];
    hb_residual(problem, pack, u, 1.0, &mut aft, &mut residual);
    let mut workspace = Workspace::new(pack).unwrap();
    workspace
        .prepare(
            problem,
            problem.island_state_jacobian().unwrap(),
            pack,
            1.0,
            &aft,
        )
        .unwrap();
    (0..pack.coeff_len())
        .map(|row| {
            direction
                .iter()
                .enumerate()
                .map(|(col, v)| workspace.coefficient_entry(row, col, pack) * v)
                .sum()
        })
        .collect()
}

#[test]
fn cubic_aft_derivative_matches_independent_full_convolution() {
    struct Cubic;
    impl OrbitProblem for Cubic {
        fn dim(&self) -> usize {
            2
        }
        fn island(&self, t: f64, x: &[f64], out: &mut [f64]) {
            out[0] = det::cos(t) - x[0];
            out[1] = x[0].powi(3) - x[1];
        }
        fn port(&self, s: C64) -> Vec<C64> {
            vec![
                s + C64::from_re(0.2),
                C64::from_re(-0.3),
                C64::from_re(0.4),
                s.scale(2.0) + C64::from_re(0.1),
            ]
        }
        fn autonomous(&self) -> bool {
            false
        }
        fn island_state_jacobian(&self) -> Option<IslandStateJacobian<Self>> {
            Some(|_, _, x, out| {
                out.copy_from_slice(&[-1.0, 0.0, 3.0 * x[0] * x[0], -1.0]);
                true
            })
        }
    }
    let pack = Packing {
        d: 2,
        n_harm: 3,
        omega_slot: None,
    };
    let u: Vec<_> = (0..pack.coeff_len())
        .map(|i| 0.06 * (i as f64 - 3.5))
        .collect();
    let direction: Vec<_> = (0..pack.coeff_len())
        .map(|i| 0.03 * (4.0 - i as f64))
        .collect();
    let measured = action(&Cubic, &pack, &u, &direction);
    for n in 0..=3 {
        let dx = coefficient(&direction, &pack, 0, n);
        let dy = coefficient(&direction, &pack, 1, n);
        let mut cubic = C64::ZERO;
        for a in -3..=3 {
            for b in -3..=3 {
                cubic = cubic
                    + coefficient(&direction, &pack, 0, a)
                        * coefficient(&u, &pack, 0, b)
                        * coefficient(&u, &pack, 0, n - a - b);
            }
        }
        let expected = [
            C64::new(1.2, n as f64) * dx - dy.scale(0.3),
            dx.scale(0.4) + C64::new(1.1, 2.0 * n as f64) * dy - cubic.scale(3.0),
        ];
        for (state, expected) in expected.iter().enumerate() {
            let actual = coefficient(&measured, &pack, state, n);
            assert!(
                (actual - *expected).abs() < 2.0e-13,
                "n={n}, state={state}: {actual:?} != {expected:?}"
            );
        }
    }
}

#[test]
fn explicit_periodic_jacobian_retains_cross_harmonic_phase_and_dc() {
    struct Periodic;
    impl OrbitProblem for Periodic {
        fn dim(&self) -> usize {
            1
        }
        fn island(&self, t: f64, x: &[f64], out: &mut [f64]) {
            out[0] = (-1.0 + 0.4 * det::cos(2.0 * t) + 0.3 * det::sin(t)) * x[0];
        }
        fn autonomous(&self) -> bool {
            false
        }
        fn island_state_jacobian(&self) -> Option<IslandStateJacobian<Self>> {
            Some(|_, t, _, out| {
                out[0] = -1.0 + 0.4 * det::cos(2.0 * t) + 0.3 * det::sin(t);
                true
            })
        }
    }
    let pack = Packing {
        d: 1,
        n_harm: 3,
        omega_slot: None,
    };
    let direction = vec![0.3, 0.2, -0.4, -0.1, 0.15, 0.07, 0.05];
    let measured = action(&Periodic, &pack, &[0.0; 7], &direction);
    // Independently authored Fourier coefficients of the time-local tangent.
    let a = |order: isize| match order {
        0 => C64::from_re(-1.0),
        1 => C64::new(0.0, -0.15),
        -1 => C64::new(0.0, 0.15),
        2 | -2 => C64::from_re(0.2),
        _ => C64::ZERO,
    };
    for n in 0..=3 {
        let mut expected = C64::new(0.0, n as f64) * coefficient(&direction, &pack, 0, n);
        for k in -3..=3 {
            expected = expected - a(n - k) * coefficient(&direction, &pack, 0, k);
        }
        let actual = coefficient(&measured, &pack, 0, n);
        assert!(
            (actual - expected).abs() < 2.0e-14,
            "n={n}: {actual:?} != {expected:?}"
        );
    }
}

#[test]
fn constant_jacobian_has_exactly_zero_dc_to_ac_coupling() {
    struct Constant;
    impl OrbitProblem for Constant {
        fn dim(&self) -> usize { 1 }
        fn island(&self, _: f64, x: &[f64], out: &mut [f64]) { out[0] = 1.0e12 * x[0]; }
        fn autonomous(&self) -> bool { false }
        fn island_state_jacobian(&self) -> Option<IslandStateJacobian<Self>> {
            Some(|_, _, _, out| { out[0] = 1.0e12; true })
        }
    }
    let pack = Packing { d: 1, n_harm: 5, omega_slot: None };
    for samples in [32, 65] {
        let mut aft = AftWorkspace::with_samples(&pack, samples, 0).unwrap();
        let mut residual = vec![0.0; pack.coeff_len()];
        hb_residual(&Constant, &pack, &vec![0.0; pack.coeff_len()], 1.0, &mut aft, &mut residual);
        let mut workspace = Workspace::new(&pack).unwrap();
        workspace.prepare(&Constant, Constant.island_state_jacobian().unwrap(), &pack, 1.0, &aft).unwrap();
        assert_eq!(workspace.coefficient_entry(0, 0, &pack), -1.0e12);
        for row in 1..pack.coeff_len() {
            assert_eq!(workspace.coefficient_entry(row, 0, &pack), 0.0);
        }
    }
}

#[test]
fn centered_nonlinear_tangent_matches_sample_products_and_finite_differences() {
    struct Coupled;
    impl OrbitProblem for Coupled {
        fn dim(&self) -> usize { 2 }
        fn island(&self, t: f64, x: &[f64], out: &mut [f64]) {
            out[0] = 1234.0 + (1.0 + 0.4 * det::cos(t)) * x[0]
                + 0.3 * x[0] * x[1] + 0.1 * x[1].powi(3);
            out[1] = -567.0 + (2.0 + 0.2 * det::sin(t)) * x[1]
                - 0.7 * x[0] * x[0] + 0.25 * x[0] * x[1];
        }
        fn port(&self, s: C64) -> Vec<C64> {
            vec![C64::from_re(1.0) + s, C64::from_re(0.2), C64::from_re(-0.4), C64::from_re(2.0) + s.scale(0.5)]
        }
        fn autonomous(&self) -> bool { false }
        fn island_state_jacobian(&self) -> Option<IslandStateJacobian<Self>> {
            Some(|_, t, x, out| {
                out.copy_from_slice(&[
                    1.0 + 0.4 * det::cos(t) + 0.3 * x[1],
                    0.3 * x[0] + 0.3 * x[1] * x[1],
                    -1.4 * x[0] + 0.25 * x[1],
                    2.0 + 0.2 * det::sin(t) + 0.25 * x[0],
                ]);
                true
            })
        }
    }
    let pack = Packing { d: 2, n_harm: 3, omega_slot: None };
    let u: Vec<_> = (0..pack.coeff_len()).map(|i| 0.025 * (i as f64 - 4.5)).collect();
    let mut aft = AftWorkspace::with_samples(&pack, 33, 0).unwrap();
    let mut residual = vec![0.0; pack.coeff_len()];
    hb_residual(&Coupled, &pack, &u, 1.0, &mut aft, &mut residual);
    let callback = Coupled.island_state_jacobian().unwrap();
    let mut workspace = Workspace::new(&pack).unwrap();
    workspace.prepare(&Coupled, callback, &pack, 1.0, &aft).unwrap();
    let jacobians: Vec<_> = (0..aft.samples).map(|j| {
        let mut jacobian = [0.0; 4];
        assert!(callback(&Coupled, TAU * j as f64 / aft.samples as f64, &aft.xs[2*j..2*j+2], &mut jacobian));
        jacobian
    }).collect();
    let basis = |j: usize, n: usize, imaginary: bool| {
        if n == 0 { 1.0 }
        else if imaginary { -2.0 * aft.sine[j * pack.n_harm + n - 1] }
        else { 2.0 * aft.cosine[j * pack.n_harm + n - 1] }
    };
    for col in 0..pack.coeff_len() {
        let (k, input, imaginary_col) = coordinate(col, pack.d);
        for row in 0..pack.coeff_len() {
            let (n, output, imaginary_row) = coordinate(row, pack.d);
            let entry = output * pack.d + input;
            let reference = jacobians[0][entry] * basis(0, k, imaginary_col);
            let mut sampled = C64::from_re(if n == 0 { reference } else { 0.0 });
            for (j, jacobian) in jacobians.iter().enumerate().skip(1) {
                let weight = if n == 0 { C64::from_re(1.0) }
                    else { C64::new(aft.cosine[j * pack.n_harm + n - 1], -aft.sine[j * pack.n_harm + n - 1]) };
                sampled = sampled + weight.scale((jacobian[entry] * basis(j, k, imaginary_col) - reference) / aft.samples as f64);
            }
            let mut linear = if n == k { Coupled.port(C64::new(0.0, n as f64))[entry] } else { C64::ZERO };
            if imaginary_col { linear = C64::new(-linear.im, linear.re); }
            let expected = linear - sampled;
            let expected = if imaginary_row { expected.im } else { expected.re };
            assert!((workspace.coefficient_entry(row, col, &pack) - expected).abs() < 2.0e-12,
                "direct sampled derivative at ({row}, {col})");
        }
    }
    let h = 1.0e-4;
    for col in 0..pack.coeff_len() {
        let mut plus = u.clone();
        let mut minus = u.clone();
        plus[col] += h;
        minus[col] -= h;
        let mut rp = vec![0.0; pack.coeff_len()];
        let mut rm = rp.clone();
        hb_residual(&Coupled, &pack, &plus, 1.0, &mut aft, &mut rp);
        hb_residual(&Coupled, &pack, &minus, 1.0, &mut aft, &mut rm);
        for row in 0..pack.coeff_len() {
            let difference = (rp[row] - rm[row]) / (2.0 * h);
            assert!((workspace.coefficient_entry(row, col, &pack) - difference).abs() < 3.0e-8,
                "finite-difference derivative at ({row}, {col})");
        }
    }
}
