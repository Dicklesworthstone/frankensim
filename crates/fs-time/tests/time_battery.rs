//! fs-time battery (tfz.12): symplectic energy boundedness vs RK4's
//! secular drift, the variational-integrator equivalence, Lie-group
//! norm preservation + gyroscope physics, generalized-α spectral
//! behavior against Chung–Hulbert theory, IMEX stiff stability +
//! second-order convergence, exponential-Euler exactness on linears,
//! embedded-pair adaptivity with bitwise resumability, the discrete
//! adjoint gradcheck, and the cross-ISA golden hash.
//!
//! Test-side ORACLES may use std/platform libm (disjoint-path rule);
//! everything feeding the golden hash flows through solver code only.

use fs_time::{
    AdaptiveState, ExpEuler, GeneralizedAlpha, Imex2, PiController, So3, Vec3, galpha_step,
    imex2_step, rigid_body_step, rk45_adaptive, so3_body_exp_step, verlet_adjoint, verlet_step,
};

fn log(case: &str, verdict: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-time\",\"case\":\"{case}\",\"verdict\":\"{verdict}\",\"detail\":\"{detail}\"}}"
    );
}

// ---------------------------------------------------------------- symplectic

/// Harmonic-oscillator energy for unit mass/stiffness.
fn ho_energy(q: &[f64], p: &[f64]) -> f64 {
    f64::midpoint(q[0] * q[0], p[0] * p[0])
}

/// Classic RK4 on the first-order form of q̈ = −q (the drift comparator).
fn rk4_ho_step(q: &mut f64, p: &mut f64, h: f64) {
    let f = |q: f64, p: f64| (p, -q);
    let (k1q, k1p) = f(*q, *p);
    let (k2q, k2p) = f(*q + 0.5 * h * k1q, *p + 0.5 * h * k1p);
    let (k3q, k3p) = f(*q + 0.5 * h * k2q, *p + 0.5 * h * k2p);
    let (k4q, k4p) = f(*q + h * k3q, *p + h * k3p);
    *q += h / 6.0 * (k1q + 2.0 * k2q + 2.0 * k3q + k4q);
    *p += h / 6.0 * (k1p + 2.0 * k2p + 2.0 * k3p + k4p);
}

#[test]
fn verlet_energy_bounded_1e6_steps_vs_rk4_drift() {
    let h = 0.1;
    let steps = 1_000_000usize;
    let force = |q: &[f64], out: &mut [f64]| out[0] = -q[0];
    let (mut q, mut p) = (vec![1.0f64], vec![0.0f64]);
    let mut scratch = vec![0.0f64; 1];
    let e0 = ho_energy(&q, &p);
    let (mut max_dev_first, mut max_dev_second) = (0.0f64, 0.0f64);
    for k in 0..steps {
        verlet_step(&mut q, &mut p, h, &force, &mut scratch);
        let dev = (ho_energy(&q, &p) - e0).abs();
        if k < steps / 2 {
            max_dev_first = max_dev_first.max(dev);
        } else {
            max_dev_second = max_dev_second.max(dev);
        }
    }
    // Symplectic ⇒ BOUNDED oscillating energy error: the second half of
    // the run must not exceed the first half beyond roundoff slack.
    assert!(
        max_dev_second <= 1.05 * max_dev_first + 1e-12,
        "secular drift in Verlet: first-half max {max_dev_first:.3e}, second {max_dev_second:.3e}"
    );
    assert!(
        max_dev_first < 2e-3,
        "Verlet energy deviation too large: {max_dev_first:.3e}"
    );
    // RK4 comparator: same h, same span — visible SECULAR energy decay.
    let (mut qr, mut pr) = (1.0f64, 0.0f64);
    for _ in 0..steps {
        rk4_ho_step(&mut qr, &mut pr, h);
    }
    let rk4_drift = (f64::midpoint(qr * qr, pr * pr) - e0).abs();
    assert!(
        rk4_drift > 5.0 * max_dev_second,
        "RK4 should drift visibly: {rk4_drift:.3e} vs Verlet bound {max_dev_second:.3e}"
    );
    log(
        "verlet-vs-rk4",
        "pass",
        &format!("verlet_bound={max_dev_second:.3e} rk4_drift={rk4_drift:.3e}"),
    );
}

#[test]
fn verlet_energy_bounded_kepler_eccentric() {
    // e = 0.6 Kepler orbit, ~16 revolutions: energy bounded, no drift.
    let h = 1e-3;
    let steps = 100_000usize;
    let force = |q: &[f64], out: &mut [f64]| {
        let r2 = q[0] * q[0] + q[1] * q[1];
        let r3 = r2 * r2.sqrt();
        out[0] = -q[0] / r3;
        out[1] = -q[1] / r3;
    };
    let energy = |q: &[f64], p: &[f64]| {
        f64::midpoint(p[0] * p[0], p[1] * p[1]) - 1.0 / (q[0] * q[0] + q[1] * q[1]).sqrt()
    };
    let (mut q, mut p) = (vec![0.4f64, 0.0], vec![0.0f64, 2.0]);
    let mut scratch = vec![0.0f64; 2];
    let e0 = energy(&q, &p);
    let mut max_dev = 0.0f64;
    for _ in 0..steps {
        verlet_step(&mut q, &mut p, h, &force, &mut scratch);
        max_dev = max_dev.max((energy(&q, &p) - e0).abs());
    }
    assert!(max_dev < 5e-4, "Kepler energy deviation {max_dev:.3e}");
    log(
        "verlet-kepler",
        "pass",
        &format!("e=0.6 max_energy_dev={max_dev:.3e}"),
    );
}

#[test]
fn verlet_is_the_variational_integrator() {
    // Marsden–West: extremizing Σ h·[½‖(q_{k+1}−q_k)/h‖² − ½(V_k+V_{k+1})]
    // gives the position two-step q_{k+1} = 2q_k − q_{k−1} + h²F(q_k).
    // The kick–drift–kick positions must satisfy it (same method, two
    // derivations) to accumulated-roundoff accuracy on a NONLINEAR V.
    let h = 0.05;
    let force = |q: &[f64], out: &mut [f64]| out[0] = -q[0].sin(); // pendulum
    let (mut q, mut p) = (vec![1.2f64], vec![0.3f64]);
    let mut scratch = vec![0.0f64; 1];
    let mut positions = vec![q[0]];
    for _ in 0..200 {
        verlet_step(&mut q, &mut p, h, &force, &mut scratch);
        positions.push(q[0]);
    }
    let mut worst = 0.0f64;
    for k in 1..200 {
        let predicted = 2.0 * positions[k] - positions[k - 1] - h * h * positions[k].sin();
        worst = worst.max((positions[k + 1] - predicted).abs());
    }
    assert!(
        worst < 1e-12,
        "discrete Euler–Lagrange residual {worst:.3e}"
    );
    log(
        "verlet-variational",
        "pass",
        &format!("max_dEL_residual={worst:.3e}"),
    );
}

// G0/G1/G3/G4: physical momentum, coupled mass, exact discrete pullbacks and
// refusal without partially advancing a caller's state.
mod momentum_verlet {
    use fs_solver::op::LinearOp;
    use fs_time::{
        MomentumVerletError, verlet_adjoint, verlet_momentum_adjoint, verlet_momentum_step,
        verlet_step,
    };
    use std::cell::Cell;

    struct Diagonal<const N: usize>([f64; N]);
    impl<const N: usize> LinearOp for Diagonal<N> {
        fn n(&self) -> usize {
            N
        }
        fn apply(&self, x: &[f64], out: &mut [f64]) {
            for i in 0..N {
                out[i] = self.0[i] * x[i];
            }
        }
    }
    struct Matrix([[f64; 2]; 2]);
    impl LinearOp for Matrix {
        fn n(&self) -> usize {
            2
        }
        fn apply(&self, x: &[f64], out: &mut [f64]) {
            for i in 0..2 {
                out[i] = self.0[i][0] * x[0] + self.0[i][1] * x[1];
            }
        }
        fn apply_transpose(&self, x: &[f64], out: &mut [f64]) {
            for i in 0..2 {
                out[i] = self.0[0][i] * x[0] + self.0[1][i] * x[1];
            }
        }
    }

    #[test]
    fn unequal_mass_internal_forces_preserve_physical_momenta_and_reverse() {
        let inverse_mass = Diagonal([0.5, 0.5, 1.0 / 3.0, 1.0 / 3.0]);
        let force = |q: &[f64], out: &mut [f64]| {
            for axis in 0..2 {
                out[axis] = q[axis + 2] - q[axis];
                out[axis + 2] = -out[axis];
            }
        };
        let initial = ([-0.6, 0.2, 0.4, -0.3], [0.8, -0.5, -0.2, 1.1]);
        let (mut q, mut p) = initial;
        let mut scratch = [0.0; 12];
        let angular = |q: &[f64], p: &[f64]| q[0] * p[1] - q[1] * p[0] + q[2] * p[3] - q[3] * p[2];
        let energy = |q: &[f64], p: &[f64]| {
            p.iter()
                .zip(inverse_mass.0)
                .map(|(p, b)| 0.5 * p * p * b)
                .sum::<f64>()
                + 0.5 * ((q[2] - q[0]).powi(2) + (q[3] - q[1]).powi(2))
        };
        let l0 = angular(&q, &p);
        let e0 = energy(&q, &p);
        for _ in 0..4000 {
            verlet_momentum_step(
                &mut q,
                &mut p,
                0.01,
                &inverse_mass,
                &force,
                &mut scratch,
                &mut || false,
            )
            .unwrap();
            assert!((p[0] + p[2] - 0.6).abs() < 2e-12);
            assert!((p[1] + p[3] - 0.6).abs() < 2e-12);
            assert!((angular(&q, &p) - l0).abs() < 2e-10);
            assert!((energy(&q, &p) - e0).abs() < 5e-5);
        }
        for axis in 0..2 {
            let com_change = (2.0 * (q[axis] - initial.0[axis])
                + 3.0 * (q[axis + 2] - initial.0[axis + 2]))
                / 5.0;
            assert!((com_change - 40.0 * 0.6 / 5.0).abs() < 2e-10);
        }
        for _ in 0..4000 {
            verlet_momentum_step(
                &mut q,
                &mut p,
                -0.01,
                &inverse_mass,
                &force,
                &mut scratch,
                &mut || false,
            )
            .unwrap();
        }
        for i in 0..4 {
            assert!((q[i] - initial.0[i]).abs() < 1e-9);
            assert!((p[i] - initial.1[i]).abs() < 1e-9);
        }
    }

    #[test]
    fn coupled_consistent_mass_matches_independent_oscillator_and_second_order() {
        // M = [[2,1],[1,2]], K = M: both coordinate modes have frequency one.
        let inverse_mass = Matrix([[2.0 / 3.0, -1.0 / 3.0], [-1.0 / 3.0, 2.0 / 3.0]]);
        let force = |q: &[f64], out: &mut [f64]| {
            out[0] = -2.0 * q[0] - q[1];
            out[1] = -q[0] - 2.0 * q[1];
        };
        let (s, c) = 1.0f64.sin_cos();
        let exact_q = [0.7 * c + 0.2 * s, -0.4 * c + 0.3 * s];
        let exact_v = [-0.7 * s + 0.2 * c, 0.4 * s + 0.3 * c];
        let exact_p = [2.0 * exact_v[0] + exact_v[1], exact_v[0] + 2.0 * exact_v[1]];
        let mut errors = Vec::new();
        for steps in [25, 50, 100] {
            let (mut q, mut p) = ([0.7, -0.4], [0.7, 0.8]);
            let mut scratch = [0.0; 6];
            for _ in 0..steps {
                verlet_momentum_step(
                    &mut q,
                    &mut p,
                    1.0 / f64::from(steps),
                    &inverse_mass,
                    &force,
                    &mut scratch,
                    &mut || false,
                )
                .unwrap();
            }
            errors.push(
                (0..2)
                    .map(|i| (q[i] - exact_q[i]).powi(2) + (p[i] - exact_p[i]).powi(2))
                    .sum::<f64>()
                    .sqrt(),
            );
        }
        for pair in errors.windows(2) {
            assert!((3.9..4.1).contains(&(pair[0] / pair[1])), "{errors:?}");
        }
        assert!(errors[2] < 2e-5, "{errors:?}");
    }

    #[test]
    fn identity_mass_retains_original_step_and_adjoint_bits() {
        let identity = Diagonal([1.0, 1.0]);
        let force = |q: &[f64], out: &mut [f64]| {
            out[0] = -q[0].sin();
            out[1] = -2.0 * q[1];
        };
        let vjp = |q: &[f64], seed: &[f64], out: &mut [f64]| {
            out[0] = -q[0].cos() * seed[0];
            out[1] = -2.0 * seed[1];
        };
        let initial = ([0.8, -0.4], [-0.2, 0.3]);
        let (mut q, mut p) = initial;
        let (mut old_q, mut old_p) = initial;
        for _ in 0..31 {
            verlet_step(&mut old_q, &mut old_p, 0.03, &force, &mut [0.0; 2]);
            verlet_momentum_step(
                &mut q,
                &mut p,
                0.03,
                &identity,
                &force,
                &mut [0.0; 6],
                &mut || false,
            )
            .unwrap();
            assert_eq!(q.map(f64::to_bits), old_q.map(f64::to_bits));
            assert_eq!(p.map(f64::to_bits), old_p.map(f64::to_bits));
        }
        let seed = (&[0.7, -0.2][..], &[-0.5, 0.9][..]);
        let old = verlet_adjoint(&initial.0, &initial.1, 0.03, 31, &force, &vjp, seed);
        let new = verlet_momentum_adjoint(
            &initial.0,
            &initial.1,
            0.03,
            31,
            &identity,
            &force,
            &vjp,
            seed,
            &mut || false,
        )
        .unwrap();
        assert_eq!(
            old.0.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
            new.0.iter().map(|x| x.to_bits()).collect::<Vec<_>>()
        );
        assert_eq!(
            old.1.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
            new.1.iter().map(|x| x.to_bits()).collect::<Vec<_>>()
        );
    }

    #[test]
    fn adjoint_uses_both_true_transposes_and_matches_finite_difference() {
        // Deliberately nonsymmetric actions test the numerical map's transpose;
        // this fixture makes no physical positive-mass or symplectic claim.
        let inverse_mass = Matrix([[0.7, 0.2], [-0.1, 0.5]]);
        let force = |q: &[f64], out: &mut [f64]| {
            out[0] = -q[0] - 0.2 * q[0].powi(3) + 0.3 * q[1];
            out[1] = -0.4 * q[0] - 1.6 * q[1];
        };
        let vjp = |q: &[f64], seed: &[f64], out: &mut [f64]| {
            out[0] = (-1.0 - 0.6 * q[0].powi(2)) * seed[0] - 0.4 * seed[1];
            out[1] = 0.3 * seed[0] - 1.6 * seed[1];
        };
        let initial = ([0.8, -0.4], [-0.2, 0.3]);
        let seed = ([0.7, -0.2], [-0.5, 0.9]);
        let cost = |mut q: [f64; 2], mut p: [f64; 2]| {
            for _ in 0..17 {
                verlet_momentum_step(
                    &mut q,
                    &mut p,
                    0.03,
                    &inverse_mass,
                    &force,
                    &mut [0.0; 6],
                    &mut || false,
                )
                .unwrap();
            }
            (0..2)
                .map(|i| seed.0[i] * q[i] + seed.1[i] * p[i])
                .sum::<f64>()
        };
        let actual = verlet_momentum_adjoint(
            &initial.0,
            &initial.1,
            0.03,
            17,
            &inverse_mass,
            &force,
            &vjp,
            (&seed.0, &seed.1),
            &mut || false,
        )
        .unwrap();
        let epsilon = 1e-5;
        for i in 0..2 {
            let (mut high, mut low) = (initial.0, initial.0);
            high[i] += epsilon;
            low[i] -= epsilon;
            let expected = (cost(high, initial.1) - cost(low, initial.1)) / (2.0 * epsilon);
            assert!((actual.0[i] - expected).abs() < 5e-9);
            let (mut high, mut low) = (initial.1, initial.1);
            high[i] += epsilon;
            low[i] -= epsilon;
            let expected = (cost(initial.0, high) - cost(initial.0, low)) / (2.0 * epsilon);
            assert!((actual.1[i] - expected).abs() < 5e-9);
        }
        let empty = verlet_momentum_adjoint(
            &initial.0,
            &initial.1,
            0.03,
            0,
            &inverse_mass,
            &force,
            &vjp,
            (&seed.0, &seed.1),
            &mut || false,
        )
        .unwrap();
        assert_eq!(empty, (seed.0.to_vec(), seed.1.to_vec()));
    }

    #[test]
    fn invalid_actions_and_late_cancellation_never_publish_partial_state() {
        let identity = Diagonal([1.0, 1.0]);
        let initial = ([1.0, -0.5], [0.2, 0.3]);
        let (mut q, mut p) = initial;
        let force = |q: &[f64], out: &mut [f64]| {
            out[0] = -q[0];
            out[1] = -q[1];
        };
        for (h, scratch_length) in [(0.0, 6), (0.1, 5)] {
            assert!(matches!(
                verlet_momentum_step(
                    &mut q,
                    &mut p,
                    h,
                    &identity,
                    &force,
                    &mut vec![0.0; scratch_length],
                    &mut || false
                ),
                Err(MomentumVerletError::InvalidInput(_))
            ));
            assert_eq!((q, p), initial);
        }
        assert!(matches!(
            verlet_momentum_step(
                &mut q,
                &mut p,
                0.1,
                &Diagonal([1.0]),
                &force,
                &mut [0.0; 6],
                &mut || false
            ),
            Err(MomentumVerletError::InvalidInput(_))
        ));
        let incomplete_force = |_: &[f64], out: &mut [f64]| {
            out[0] = 0.0;
        };
        assert!(matches!(
            verlet_momentum_step(
                &mut q,
                &mut p,
                0.1,
                &identity,
                &incomplete_force,
                &mut [0.0; 6],
                &mut || false
            ),
            Err(MomentumVerletError::NonFinite(_))
        ));
        let calls = Cell::new(0usize);
        let counted = |q: &[f64], out: &mut [f64]| {
            calls.set(calls.get() + 1);
            force(q, out);
        };
        assert_eq!(
            verlet_momentum_step(
                &mut q,
                &mut p,
                0.1,
                &identity,
                &counted,
                &mut [0.0; 6],
                &mut || calls.get() == 2
            ),
            Err(MomentumVerletError::Cancelled)
        );
        assert_eq!((q, p), initial);
        calls.set(0);
        let vjp = |_: &[f64], seed: &[f64], out: &mut [f64]| {
            out[0] = -seed[0];
            out[1] = -seed[1];
        };
        // A pullback must admit the complete primal step, including the final
        // force kick, even though its reverse arithmetic only needs q_next.
        let invalid_endpoint = |position: &[f64], out: &mut [f64]| {
            out.fill(if position[0] == initial.0[0] {
                0.0
            } else {
                f64::NAN
            });
        };
        assert_eq!(
            verlet_momentum_adjoint(
                &q,
                &p,
                0.1,
                1,
                &identity,
                &invalid_endpoint,
                &vjp,
                (&q, &p),
                &mut || false
            ),
            Err(MomentumVerletError::NonFinite("endpoint force"))
        );
        assert_eq!(
            verlet_momentum_adjoint(
                &q,
                &p,
                0.1,
                1 << 20,
                &identity,
                &counted,
                &vjp,
                (&q, &p),
                &mut || calls.get() >= 3
            ),
            Err(MomentumVerletError::Cancelled)
        );
        assert_eq!(
            calls.get(),
            3,
            "cancellation must stop checkpoint replay immediately"
        );
    }

    mod parameters {
        use super::*;
        use fs_ad::{Dual64, Real};
        use fs_time::{
            MomentumVerletGradient, MomentumVerletParameterVjp, verlet_momentum_parameter_adjoint,
        };

        struct SpringDesign([f64; 5]);

        impl LinearOp for SpringDesign {
            fn n(&self) -> usize {
                2
            }
            fn apply(&self, p: &[f64], out: &mut [f64]) {
                let [m1, m2, coupling, _, _] = self.0;
                let determinant = m1 * m2 - coupling * coupling;
                out[0] = (m2 * p[0] - coupling * p[1]) / determinant;
                out[1] = (m1 * p[1] - coupling * p[0]) / determinant;
            }
        }

        impl SpringDesign {
            fn force(&self, q: &[f64], out: &mut [f64]) {
                let extension = q[0] - q[1];
                out[0] = -self.0[3] * extension - self.0[4] * extension.powi(3);
                out[1] = -out[0];
            }
            fn state_vjp(&self, q: &[f64], seed: &[f64], out: &mut [f64]) {
                let tangent = self.0[3] + 3.0 * self.0[4] * (q[0] - q[1]).powi(2);
                out[0] = tangent * (seed[1] - seed[0]);
                out[1] = -out[0];
            }
            fn trajectory(&self) -> ([f64; 2], [f64; 2]) {
                let (mut q, mut p) = ([0.7, -0.25], [0.2, -0.35]);
                for _ in 0..33 {
                    verlet_momentum_step(
                        &mut q,
                        &mut p,
                        0.025,
                        self,
                        &|q, f| self.force(q, f),
                        &mut [0.0; 6],
                        &mut || false,
                    )
                    .unwrap();
                }
                (q, p)
            }
            fn cost(&self) -> f64 {
                let (q, p) = self.trajectory();
                0.5 * ((q[0] - 0.1).powi(2) + (q[1] + 0.2).powi(2))
                    + 0.05 * (p[0] * p[0] + p[1] * p[1])
            }
            fn gradient(&self) -> MomentumVerletGradient {
                let (q, p) = self.trajectory();
                verlet_momentum_parameter_adjoint(
                    &[0.7, -0.25],
                    &[0.2, -0.35],
                    0.025,
                    33,
                    self,
                    &|q, f| self.force(q, f),
                    &|q, seed, out| self.state_vjp(q, seed, out),
                    self,
                    (&[q[0] - 0.1, q[1] + 0.2], &[0.1 * p[0], 0.1 * p[1]]),
                    &mut || false,
                )
                .unwrap()
            }
        }

        impl MomentumVerletParameterVjp for SpringDesign {
            fn parameter_count(&self) -> usize {
                5
            }
            fn force_vjp(
                &self,
                q: &[f64],
                seed: &[f64],
                out: &mut [f64],
            ) -> Result<(), MomentumVerletError> {
                out.fill(0.0);
                let extension = q[0] - q[1];
                out[3] = extension * (seed[1] - seed[0]);
                out[4] = extension.powi(3) * (seed[1] - seed[0]);
                Ok(())
            }
            fn mass_vjp(
                &self,
                v: &[f64],
                seed: &[f64],
                out: &mut [f64],
            ) -> Result<(), MomentumVerletError> {
                out.copy_from_slice(&[
                    seed[0] * v[0],
                    seed[1] * v[1],
                    seed[0] * v[1] + seed[1] * v[0],
                    0.0,
                    0.0,
                ]);
                Ok(())
            }
        }

        #[test]
        fn coupled_mass_and_nonlinear_force_gradients_match_independent_duals_and_differences() {
            let model = SpringDesign([2.0, 3.0, 0.4, 4.0, 0.7]);
            let actual = model.gradient();
            // Independently differentiate the scalar two-body equations forward,
            // including the explicit 2x2 inverse and all nine initial/design inputs.
            type D = Dual64<9>;
            let mut q = [D::variable(0.7, 0), D::variable(-0.25, 1)];
            let mut p = [D::variable(0.2, 2), D::variable(-0.35, 3)];
            let [m1, m2, c, k, a] = std::array::from_fn(|i| D::variable(model.0[i], i + 4));
            let determinant = m1 * m2 - c * c;
            let half = D::constant(0.0125);
            let h = D::constant(0.025);
            for _ in 0..33 {
                let x = q[0] - q[1];
                let f = -(k * x + a * x * x * x);
                p[0] = half.mul_add(f, p[0]);
                p[1] = half.mul_add(-f, p[1]);
                let v = [
                    (m2 * p[0] - c * p[1]) / determinant,
                    (m1 * p[1] - c * p[0]) / determinant,
                ];
                q[0] = h.mul_add(v[0], q[0]);
                q[1] = h.mul_add(v[1], q[1]);
                let x = q[0] - q[1];
                let f = -(k * x + a * x * x * x);
                p[0] = half.mul_add(f, p[0]);
                p[1] = half.mul_add(-f, p[1]);
            }
            let e0 = q[0] - D::constant(0.1);
            let e1 = q[1] + D::constant(0.2);
            let objective = D::constant(0.5) * (e0 * e0 + e1 * e1)
                + D::constant(0.05) * (p[0] * p[0] + p[1] * p[1]);
            assert!((objective.re - model.cost()).abs() < 2e-14);
            let all = actual
                .initial_position
                .iter()
                .chain(&actual.initial_momentum)
                .chain(&actual.parameters);
            for (index, (&reverse, forward)) in all.zip(objective.eps).enumerate() {
                assert!(
                    (reverse - forward).abs() < 3e-12,
                    "component {index}: reverse={reverse:e}, dual={forward:e}"
                );
            }
            for i in 0..5 {
                let epsilon = 1e-5;
                let (mut plus, mut minus) = (model.0, model.0);
                plus[i] += epsilon;
                minus[i] -= epsilon;
                let finite_difference =
                    (SpringDesign(plus).cost() - SpringDesign(minus).cost()) / (2.0 * epsilon);
                assert!((actual.parameters[i] - finite_difference).abs() < 2e-9);
            }
        }

        #[test]
        fn mass_and_force_design_updates_reduce_the_actual_trajectory_objective() {
            let mut model = SpringDesign([2.0, 3.0, 0.4, 4.0, 0.7]);
            let initial = model.cost();
            for _ in 0..6 {
                let before = model.cost();
                let gradient = model.gradient().parameters;
                let norm2 = gradient.iter().map(|g| g * g).sum::<f64>();
                let mut accepted = None;
                for backtrack in 0..8 {
                    let rate = 0.5 * 0.5_f64.powi(backtrack);
                    let next = std::array::from_fn(|i| model.0[i] - rate * gradient[i]);
                    assert!(next[0] > 0.0 && next[0] * next[1] > next[2] * next[2]);
                    assert!(next[3] > 0.0 && next[4] > 0.0);
                    let trial = SpringDesign(next);
                    if trial.cost() < before - 1e-4 * rate * norm2 {
                        accepted = Some(trial);
                        break;
                    }
                }
                model = accepted.expect("physical gradient must produce descent");
            }
            println!(
                "momentum design: initial={initial:e}, final={:e}",
                model.cost()
            );
            assert!(model.cost() < 0.95 * initial);
        }

        struct Ballistic {
            mass: f64,
            force: f64,
        }
        impl LinearOp for Ballistic {
            fn n(&self) -> usize {
                1
            }
            fn apply(&self, p: &[f64], out: &mut [f64]) {
                out[0] = p[0] / self.mass;
            }
        }
        impl MomentumVerletParameterVjp for Ballistic {
            fn parameter_count(&self) -> usize {
                2
            }
            fn force_vjp(
                &self,
                _: &[f64],
                seed: &[f64],
                out: &mut [f64],
            ) -> Result<(), MomentumVerletError> {
                out.copy_from_slice(&[0.0, seed[0]]);
                Ok(())
            }
            fn mass_vjp(
                &self,
                v: &[f64],
                seed: &[f64],
                out: &mut [f64],
            ) -> Result<(), MomentumVerletError> {
                out.copy_from_slice(&[v[0] * seed[0], 0.0]);
                Ok(())
            }
        }

        #[test]
        fn ballistic_mass_and_thrust_gradients_use_fixed_canonical_momentum_in_both_time_directions()
         {
            let model = Ballistic {
                mass: 2.5,
                force: 0.7,
            };
            let initial_momentum = 1.3;
            for h in [-0.03, 0.03] {
                for steps in [0, 1, 17] {
                    let result = verlet_momentum_parameter_adjoint(
                        &[0.2],
                        &[initial_momentum],
                        h,
                        steps,
                        &model,
                        &|_, out| out[0] = model.force,
                        &|_, _, out| out.fill(0.0),
                        &model,
                        (&[0.8], &[-0.4]),
                        &mut || false,
                    )
                    .unwrap();
                    let time = h * steps as f64;
                    let dm = -0.8 * (time * initial_momentum + 0.5 * time * time * model.force)
                        / model.mass.powi(2);
                    let df = 0.8 * 0.5 * time * time / model.mass - 0.4 * time;
                    assert!((result.parameters[0] - dm).abs() < 2e-14);
                    assert!((result.parameters[1] - df).abs() < 2e-14);
                    assert_eq!(result.initial_position, [0.8]);
                    assert!(
                        (result.initial_momentum[0] - (-0.4 + 0.8 * time / model.mass)).abs()
                            < 2e-14
                    );
                    // Fixed-velocity initialization needs this additional chain.
                    let fixed_velocity_dm =
                        dm + result.initial_momentum[0] * initial_momentum / model.mass;
                    let expected = -0.8 * 0.5 * time * time * model.force / model.mass.powi(2)
                        - 0.4 * initial_momentum / model.mass;
                    assert!((fixed_velocity_dm - expected).abs() < 2e-14);
                }
            }
        }

        #[test]
        fn parameter_pullback_retains_initial_state_bits_and_adds_no_inverse_solve() {
            struct Counted<'a> {
                model: &'a SpringDesign,
                applications: Cell<usize>,
            }
            impl LinearOp for Counted<'_> {
                fn n(&self) -> usize {
                    2
                }
                fn apply(&self, x: &[f64], out: &mut [f64]) {
                    self.applications.set(self.applications.get() + 1);
                    self.model.apply(x, out);
                }
            }
            let model = SpringDesign([2.0, 3.0, 0.4, 4.0, 0.7]);
            let mass = Counted {
                model: &model,
                applications: Cell::new(0),
            };
            let force = |q: &[f64], f: &mut [f64]| model.force(q, f);
            let state_vjp =
                |q: &[f64], seed: &[f64], out: &mut [f64]| model.state_vjp(q, seed, out);
            let seed = (&[0.5, -0.3][..], &[-0.2, 0.8][..]);
            let original = verlet_momentum_adjoint(
                &[0.7, -0.25],
                &[0.2, -0.35],
                0.025,
                33,
                &mass,
                &force,
                &state_vjp,
                seed,
                &mut || false,
            )
            .unwrap();
            let original_count = mass.applications.replace(0);
            let expanded = verlet_momentum_parameter_adjoint(
                &[0.7, -0.25],
                &[0.2, -0.35],
                0.025,
                33,
                &mass,
                &force,
                &state_vjp,
                &model,
                seed,
                &mut || false,
            )
            .unwrap();
            assert_eq!(mass.applications.get(), original_count);
            for (a, b) in original.0.iter().chain(&original.1).zip(
                expanded
                    .initial_position
                    .iter()
                    .chain(&expanded.initial_momentum),
            ) {
                assert_eq!(a.to_bits(), b.to_bits());
            }
        }

        #[test]
        fn parameter_derivative_errors_incomplete_writes_and_cancellation_return_no_gradient() {
            struct Faulty<'a> {
                model: &'a SpringDesign,
                mode: usize,
                calls: Cell<usize>,
            }
            impl MomentumVerletParameterVjp for Faulty<'_> {
                fn parameter_count(&self) -> usize {
                    match self.mode {
                        4 => usize::MAX,
                        5 => 0,
                        _ => 5,
                    }
                }
                fn force_vjp(
                    &self,
                    q: &[f64],
                    seed: &[f64],
                    out: &mut [f64],
                ) -> Result<(), MomentumVerletError> {
                    self.calls.set(self.calls.get() + 1);
                    if self.mode == 0 {
                        return Err(MomentumVerletError::InvalidInput("derivative domain"));
                    }
                    self.model.force_vjp(q, seed, out)?;
                    if self.mode == 2 {
                        out[0] = f64::INFINITY;
                    }
                    Ok(())
                }
                fn mass_vjp(
                    &self,
                    v: &[f64],
                    seed: &[f64],
                    out: &mut [f64],
                ) -> Result<(), MomentumVerletError> {
                    self.calls.set(self.calls.get() + 1);
                    if self.mode != 1 {
                        self.model.mass_vjp(v, seed, out)?;
                    }
                    Ok(())
                }
            }
            let model = SpringDesign([2.0, 3.0, 0.4, 4.0, 0.7]);
            for mode in 0..6 {
                let parameters = Faulty {
                    model: &model,
                    mode,
                    calls: Cell::new(0),
                };
                let result = verlet_momentum_parameter_adjoint(
                    &[0.7, -0.25],
                    &[0.2, -0.35],
                    0.025,
                    1,
                    &model,
                    &|q, f| model.force(q, f),
                    &|q, seed, out| model.state_vjp(q, seed, out),
                    &parameters,
                    (&[0.5, -0.3], &[-0.2, 0.8]),
                    &mut || mode == 3 && parameters.calls.get() > 0,
                );
                match mode {
                    0 => assert_eq!(
                        result,
                        Err(MomentumVerletError::InvalidInput("derivative domain"))
                    ),
                    1 => assert_eq!(
                        result,
                        Err(MomentumVerletError::NonFinite("mass parameter derivative"))
                    ),
                    2 => assert_eq!(
                        result,
                        Err(MomentumVerletError::NonFinite(
                            "endpoint force parameter derivative"
                        ))
                    ),
                    3 => assert_eq!(result, Err(MomentumVerletError::Cancelled)),
                    4 => assert!(matches!(result, Err(MomentumVerletError::InvalidInput(_)))),
                    5 => assert!(result.unwrap().parameters.is_empty()),
                    _ => unreachable!(),
                }
                if mode >= 4 {
                    assert_eq!(parameters.calls.get(), 0);
                }
            }
        }
    }
}

// ----------------------------------------------------------------- Lie/SO(3)

#[test]
fn so3_group_invariant_preserved_1e5_steps() {
    // Time-varying ω from a deterministic recurrence; exp-map updates
    // must keep ‖q‖ = 1 to a roundoff random walk (~√N·ε ≈ 4e-14).
    let mut rotation = So3::identity();
    let mut angular = Vec3::new(0.3, -0.2, 0.9);
    for k in 0..100_000usize {
        rotation = so3_body_exp_step(rotation, angular, 0.01).expect("finite SO(3) step");
        let s = if k % 2 == 0 { 1.0 } else { -1.0 };
        angular = Vec3::new(angular.y, angular.z, angular.x + s * 1e-4);
    }
    let q = rotation.as_quat();
    let norm = (q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z).sqrt();
    // Measured ≈ 1e-12 on M4 Pro: the per-step exp-map is exact to a
    // few ulps but the det::sin/cos roundoff walk is mildly biased, so
    // drift runs ~30× above the ideal √N·ε random walk. Still 10⁻¹²
    // after 10⁵ steps with NO renormalization — the point of the test.
    assert!(
        (norm - 1.0).abs() < 5e-12,
        "quaternion norm drift: {:.3e}",
        (norm - 1.0).abs()
    );
    log(
        "quat-norm",
        "pass",
        &format!("drift={:.3e}", (norm - 1.0).abs()),
    );
}

#[test]
fn rigid_body_gyroscope_physics() {
    // Symmetric top I = (2, 2, 1), ω = (0.1, 0, 5): Euler's equations
    // give ω₃ EXACTLY constant and (ω₁, ω₂) precessing in the body
    // frame at Ω = (I₃ − I₁)/I₁·ω₃ = −2.5 rad/s. Second-order method ⇒
    // both hold to O(h²) over T = 10.
    let inertia = Vec3::new(2.0, 2.0, 1.0);
    let h = 1e-3;
    let steps = 10_000usize;
    let mut rotation = So3::identity();
    let mut angular = Vec3::new(0.1, 0.0, 5.0);
    let l0 = rotation
        .rotate(Vec3::new(
            inertia.x * angular.x,
            inertia.y * angular.y,
            inertia.z * angular.z,
        ))
        .expect("finite spatial momentum");
    let e0 = 0.5
        * (inertia.x * angular.x * angular.x
            + inertia.y * angular.y * angular.y
            + inertia.z * angular.z * angular.z);
    for _ in 0..steps {
        let (next_rotation, next_angular) =
            rigid_body_step(rotation, angular, inertia, h).expect("finite rigid-body step");
        rotation = next_rotation;
        angular = next_angular;
    }
    // (a) ω₃ constant.
    assert!(
        (angular.z - 5.0).abs() < 1e-9,
        "omega3 drift {:.3e}",
        (angular.z - 5.0).abs()
    );
    // (b) body-frame precession phase after T = 10: Ω·T = −25 rad
    // (phase(t) = −2.5t since d/dt(ω₁ + iω₂) = −2.5i(ω₁ + iω₂)).
    let two_pi = 2.0 * std::f64::consts::PI;
    let expected = (-25.0f64).rem_euclid(two_pi);
    let measured = angular.y.atan2(angular.x).rem_euclid(two_pi);
    let mut diff = (measured - expected).abs();
    if diff > std::f64::consts::PI {
        diff = two_pi - diff;
    }
    assert!(diff < 1e-3, "precession phase error {diff:.3e}");
    // (c) energy and SPATIAL angular momentum conserved to O(h²).
    let e1 = 0.5
        * (inertia.x * angular.x * angular.x
            + inertia.y * angular.y * angular.y
            + inertia.z * angular.z * angular.z);
    assert!(
        (e1 - e0).abs() / e0 < 1e-4,
        "energy drift {:.3e}",
        (e1 - e0).abs() / e0
    );
    let l1 = rotation
        .rotate(Vec3::new(
            inertia.x * angular.x,
            inertia.y * angular.y,
            inertia.z * angular.z,
        ))
        .expect("finite spatial momentum");
    let delta = l1 - l0;
    let ldev = delta.x.abs().max(delta.y.abs()).max(delta.z.abs());
    assert!(ldev < 1e-3, "spatial L drift {ldev:.3e}");
    log(
        "gyroscope",
        "pass",
        &format!(
            "phase_err={diff:.3e} L_dev={ldev:.3e} E_rel={:.3e}",
            (e1 - e0).abs() / e0
        ),
    );
}

// ------------------------------------------------------------- generalized-α

#[test]
fn galpha_high_frequency_dissipation_matches_rho_inf() {
    // 1-DOF, ωh = 10³ (deep in the high-frequency limit): the numerical
    // amplification's asymptotic per-step contraction must approach ρ∞.
    // 200 steps so the (defective-pair) k·ρᵏ transient factor 1 + 1/k
    // decays under the tolerance.
    for &(rho, annihilates) in &[(0.0f64, true), (0.5, false), (0.9, false)] {
        let ga = GeneralizedAlpha::new(&[1.0], &[0.0], &[1.0e6], 1, 1.0, rho);
        let (mut q, mut v, mut a) = (vec![1.0f64], vec![0.0f64], vec![-1.0e6f64]);
        let f = vec![0.0f64];
        let mut prev_norm = 0.0f64;
        let mut ratio = 0.0f64;
        for k in 0..200 {
            galpha_step(&ga, &mut q, &mut v, &mut a, &f);
            let norm = (q[0] * q[0] + (v[0] / 1.0e3) * (v[0] / 1.0e3)).sqrt();
            if k == 199 && prev_norm > 0.0 {
                ratio = norm / prev_norm;
            }
            prev_norm = norm;
        }
        if annihilates {
            // Asymptotic annihilation: state collapses towards zero.
            assert!(
                prev_norm < 1e-30,
                "rho=0 should annihilate: {prev_norm:.3e}"
            );
        } else {
            assert!(
                (ratio - rho).abs() < 0.02,
                "spectral radius {ratio:.4} vs rho_inf {rho} at omega*h=1e3"
            );
        }
        log(
            "galpha-spectral",
            "pass",
            &format!("rho_inf={rho} measured={ratio:.4}"),
        );
    }
}

#[test]
fn galpha_no_dissipation_at_rho_one_and_order_two() {
    // ρ∞ = 1 (trapezoidal limit): free vibration amplitude preserved.
    let ga = GeneralizedAlpha::new(&[1.0], &[0.0], &[1.0], 1, 0.1, 1.0);
    let (mut q, mut v, mut a) = (vec![1.0f64], vec![0.0f64], vec![-1.0f64]);
    let f = vec![0.0f64];
    for _ in 0..10_000 {
        galpha_step(&ga, &mut q, &mut v, &mut a, &f);
    }
    let energy = f64::midpoint(q[0] * q[0], v[0] * v[0]);
    assert!(
        (energy - 0.5).abs() < 1e-3,
        "rho=1 energy drift {:.3e}",
        (energy - 0.5).abs()
    );
    // O(h²): error over one period, halving h → error / 4 (±30%).
    // METRIC MATTERS: max over (q, v) errors — at t = 2π, cos′ = 0, so
    // a q-only error measures the phase error quadratically and fakes
    // order ≈ 4 (diagnosed by the galpha_probe sweep, kept as its own
    // regression).
    let err_at = |h: f64| -> f64 {
        let steps = (2.0 * std::f64::consts::PI / h).round() as usize;
        let hh = 2.0 * std::f64::consts::PI / steps as f64;
        let ga = GeneralizedAlpha::new(&[1.0], &[0.0], &[1.0], 1, hh, 0.8);
        let (mut q, mut v, mut a) = (vec![1.0f64], vec![0.0f64], vec![-1.0f64]);
        let f = vec![0.0f64];
        for _ in 0..steps {
            galpha_step(&ga, &mut q, &mut v, &mut a, &f);
        }
        (q[0] - 1.0).abs().max(v[0].abs())
    };
    let (e1, e2) = (err_at(0.02), err_at(0.01));
    let order = (e1 / e2).log2();
    assert!(
        (order - 2.0).abs() < 0.4,
        "generalized-alpha order {order:.2} (errors {e1:.3e}/{e2:.3e})"
    );
    log("galpha-order", "pass", &format!("order={order:.2}"));
}

// --------------------------------------------------------------------- stiff

#[test]
fn imex_stiff_stability() {
    // hλ = −100 on the implicit part: explicit methods explode, the
    // ARS(2,2,2) step contracts monotonically. Mild monotone nonlinearity.
    let im = Imex2::new(&[-1.0e4], 1, 0.01);
    let nonlin = |u: &[f64], out: &mut [f64]| out[0] = -u[0] * u[0] * u[0];
    let mut u = vec![1.0f64];
    let mut prev = u[0].abs();
    for _ in 0..200 {
        imex2_step(&im, &mut u, &nonlin);
        assert!(
            u[0].abs() <= prev + 1e-15,
            "IMEX not contracting: {} -> {}",
            prev,
            u[0]
        );
        prev = u[0].abs();
    }
    assert!(u[0].abs() < 1e-8, "stiff mode not damped: {:.3e}", u[0]);
    log("imex-stability", "pass", &format!("u_final={:.3e}", u[0]));
}

#[test]
fn imex_second_order_on_logistic() {
    // u′ = −u + u², u₀ = ½: exact u(t) = u₀ / (u₀ + (1−u₀)eᵗ).
    // L = −1 implicit, N(u) = u² explicit — order 2 requires the ARS
    // (δ, 1−δ) weights (trapezoidal ½,½ degrades to order 1 HERE).
    let exact = |t: f64| 0.5 / (0.5 + 0.5 * t.exp());
    let err_at = |h: f64| -> f64 {
        let steps = (1.0 / h).round() as usize;
        let im = Imex2::new(&[-1.0], 1, h);
        let nonlin = |u: &[f64], out: &mut [f64]| out[0] = u[0] * u[0];
        let mut u = vec![0.5f64];
        for _ in 0..steps {
            imex2_step(&im, &mut u, &nonlin);
        }
        (u[0] - exact(1.0)).abs()
    };
    let (e1, e2, e3) = (err_at(0.1), err_at(0.05), err_at(0.025));
    let (o1, o2) = ((e1 / e2).log2(), (e2 / e3).log2());
    assert!(
        (o1 - 2.0).abs() < 0.35 && (o2 - 2.0).abs() < 0.35,
        "IMEX order ratios {o1:.2}, {o2:.2} (errors {e1:.3e}/{e2:.3e}/{e3:.3e})"
    );
    log("imex-order", "pass", &format!("orders={o1:.2},{o2:.2}"));
}

#[test]
fn exp_euler_exact_on_linear() {
    // A = QᵀDQ (rotation θ = 0.3, D = diag(−2, −3)); with N ≡ 0 the
    // exponential integrator must reproduce u(t) = Qᵀe^{tD}Qu₀ to
    // roundoff over 100 steps — the oracle path (std exp + hand
    // rotation) is DISJOINT from the solver path (Jacobi eigenbasis).
    let (c, s) = (0.3f64.cos(), 0.3f64.sin());
    // A = Qᵀ D Q with Q = [[c, s], [−s, c]].
    let d = [-2.0f64, -3.0];
    let a = [
        c * d[0] * c + s * d[1] * s,
        c * d[0] * s - s * d[1] * c,
        s * d[0] * c - c * d[1] * s,
        s * d[0] * s + c * d[1] * c,
    ];
    let h = 0.05;
    let ee = ExpEuler::new(&a, 2, h);
    let zero = |_: &[f64], out: &mut [f64]| out.fill(0.0);
    let mut u = vec![0.7f64, -0.4];
    for _ in 0..100 {
        ee.step(&mut u, &zero);
    }
    let t = h * 100.0;
    let u0 = [0.7f64, -0.4];
    let w = [c * u0[0] + s * u0[1], -s * u0[0] + c * u0[1]];
    let w = [w[0] * (d[0] * t).exp(), w[1] * (d[1] * t).exp()];
    let expect = [c * w[0] - s * w[1], s * w[0] + c * w[1]];
    let err = (u[0] - expect[0]).abs().max((u[1] - expect[1]).abs());
    assert!(err < 1e-13, "exp-Euler not exact on linear: {err:.3e}");
    log("expeuler-exact", "pass", &format!("err={err:.3e}"));
}

#[test]
fn exp_euler_first_order_with_nonlinearity() {
    // ETD1 is order 1 in N; verify the convergence rate on the logistic
    // problem (same exact solution as the IMEX test).
    let exact = |t: f64| 0.5 / (0.5 + 0.5 * t.exp());
    let err_at = |h: f64| -> f64 {
        let steps = (1.0 / h).round() as usize;
        let ee = ExpEuler::new(&[-1.0], 1, h);
        let nonlin = |u: &[f64], out: &mut [f64]| out[0] = u[0] * u[0];
        let mut u = vec![0.5f64];
        for _ in 0..steps {
            ee.step(&mut u, &nonlin);
        }
        (u[0] - exact(1.0)).abs()
    };
    let (e1, e2) = (err_at(0.01), err_at(0.005));
    let order = (e1 / e2).log2();
    assert!((order - 1.0).abs() < 0.25, "ETD1 order {order:.2}");
    log("expeuler-order", "pass", &format!("order={order:.2}"));
}

// ------------------------------------------------------------------ adaptive

fn ho_rhs(_t: f64, u: &[f64], out: &mut [f64]) {
    out[0] = u[1];
    out[1] = -u[0];
}

#[test]
fn rk45_accuracy_tracks_tolerance() {
    let pi = PiController::default();
    let mut st = AdaptiveState::new(0.0, &[1.0, 0.0], 0.1);
    rk45_adaptive(&mut st, &ho_rhs, 10.0, 1e-8, 1e-10, &pi, 100_000);
    assert!((st.t - 10.0).abs() < 1e-12, "did not reach t_end: {}", st.t);
    let err = (st.u[0] - 10.0f64.cos())
        .abs()
        .max((st.u[1] + 10.0f64.sin()).abs());
    assert!(err < 1e-6, "RK45 error {err:.3e} at rtol=1e-8");
    assert!(
        st.accepted > 20 && st.accepted < 2_000,
        "step count off: {}",
        st.accepted
    );
    log(
        "rk45-accuracy",
        "pass",
        &format!(
            "err={err:.3e} accepted={} rejected={}",
            st.accepted, st.rejected
        ),
    );
}

#[test]
fn rk45_rejection_recovers_from_huge_h0() {
    let pi = PiController::default();
    let mut st = AdaptiveState::new(0.0, &[1.0, 0.0], 50.0); // absurd h₀
    rk45_adaptive(&mut st, &ho_rhs, 10.0, 1e-8, 1e-10, &pi, 100_000);
    assert!(
        st.rejected >= 1,
        "expected at least one rejection from h0=50"
    );
    let err = (st.u[0] - 10.0f64.cos())
        .abs()
        .max((st.u[1] + 10.0f64.sin()).abs());
    assert!(err < 1e-6, "post-rejection accuracy {err:.3e}");
    log(
        "rk45-reject",
        "pass",
        &format!("rejected={} err={err:.3e}", st.rejected),
    );
}

#[test]
fn rk45_resumable_split_run_bitwise() {
    // The P7 contract: state is COMPLETE. Interrupting after any number
    // of attempts (max_steps), checkpointing by clone, and resuming
    // must be bitwise-identical to the straight run — controller memory
    // (err_prev) and counters included.
    let pi = PiController::default();
    let mut straight = AdaptiveState::new(0.0, &[1.0, 0.0], 0.1);
    rk45_adaptive(&mut straight, &ho_rhs, 10.0, 1e-9, 1e-12, &pi, 100_000);
    for &cut in &[1usize, 7, 23, 61] {
        let mut first = AdaptiveState::new(0.0, &[1.0, 0.0], 0.1);
        rk45_adaptive(&mut first, &ho_rhs, 10.0, 1e-9, 1e-12, &pi, cut);
        let mut resumed = first.clone(); // checkpoint = clone
        rk45_adaptive(&mut resumed, &ho_rhs, 10.0, 1e-9, 1e-12, &pi, 100_000);
        assert_eq!(
            resumed.t.to_bits(),
            straight.t.to_bits(),
            "t differs at cut {cut}"
        );
        assert_eq!(
            resumed.h.to_bits(),
            straight.h.to_bits(),
            "h differs at cut {cut}"
        );
        assert_eq!(
            resumed.err_prev.to_bits(),
            straight.err_prev.to_bits(),
            "controller memory differs at cut {cut}"
        );
        for i in 0..2 {
            assert_eq!(
                resumed.u[i].to_bits(),
                straight.u[i].to_bits(),
                "u[{i}] at cut {cut}"
            );
        }
        assert_eq!(
            (resumed.accepted, resumed.rejected),
            (straight.accepted, straight.rejected),
            "counters at cut {cut}"
        );
    }
    log(
        "rk45-resume",
        "pass",
        "4 split points bitwise == straight run",
    );
}

// ------------------------------------------------------------------- adjoint

#[test]
fn verlet_adjoint_gradcheck_vs_central_fd() {
    // F(q) = −q − 0.3‖q‖²q (symmetric Jacobian, as the adjoint's
    // Jᵀv = Jv shortcut requires); J = ½‖q_N‖² + ½‖p_N‖².
    let n = 3usize;
    let (h, steps) = (0.01f64, 64usize);
    let force = |q: &[f64], out: &mut [f64]| {
        let r2: f64 = q.iter().map(|x| x * x).sum();
        for i in 0..q.len() {
            out[i] = -q[i] - 0.3 * r2 * q[i];
        }
    };
    let force_jvp = |q: &[f64], v: &[f64], out: &mut [f64]| {
        let r2: f64 = q.iter().map(|x| x * x).sum();
        let qv: f64 = q.iter().zip(v).map(|(a, b)| a * b).sum();
        for i in 0..q.len() {
            out[i] = -v[i] - 0.3 * (r2 * v[i] + 2.0 * qv * q[i]);
        }
    };
    let q0 = vec![0.4f64, -0.7, 0.2];
    let p0 = vec![0.1f64, 0.3, -0.5];
    let cost = |q0: &[f64], p0: &[f64]| -> f64 {
        let (mut q, mut p) = (q0.to_vec(), p0.to_vec());
        let mut scratch = vec![0.0f64; n];
        for _ in 0..steps {
            verlet_step(&mut q, &mut p, h, &force, &mut scratch);
        }
        f64::midpoint(
            q.iter().map(|x| x * x).sum::<f64>(),
            p.iter().map(|x| x * x).sum::<f64>(),
        )
    };
    // Terminal cotangent: (q_N, p_N).
    let (mut qn, mut pn) = (q0.clone(), p0.clone());
    let mut scratch = vec![0.0f64; n];
    for _ in 0..steps {
        verlet_step(&mut qn, &mut pn, h, &force, &mut scratch);
    }
    let (bar_q0, bar_p0) = verlet_adjoint(&q0, &p0, h, steps, &force, &force_jvp, (&qn, &pn));
    let eps = 1e-6;
    let mut worst = 0.0f64;
    for i in 0..n {
        let mut qp = q0.clone();
        qp[i] += eps;
        let mut qm = q0.clone();
        qm[i] -= eps;
        let fd = (cost(&qp, &p0) - cost(&qm, &p0)) / (2.0 * eps);
        worst = worst.max((fd - bar_q0[i]).abs() / fd.abs().max(1.0));
        let mut pp = p0.clone();
        pp[i] += eps;
        let mut pm = p0.clone();
        pm[i] -= eps;
        let fd = (cost(&q0, &pp) - cost(&q0, &pm)) / (2.0 * eps);
        worst = worst.max((fd - bar_p0[i]).abs() / fd.abs().max(1.0));
    }
    assert!(worst < 1e-7, "adjoint gradcheck worst rel err {worst:.3e}");
    log("verlet-adjoint", "pass", &format!("worst_rel={worst:.3e}"));
}

// --------------------------------------------------------------- golden hash

// Canonical fs-ga So3 migration: the rigid-body lane now canonicalizes the
// quaternion double cover at every accepted group operation. The other five
// lanes and the physical SO(3) trajectory are unchanged; the representative
// bits, and therefore this whole-suite digest, intentionally changed.
const GOLDEN_HASH: u64 = 0x41c7_740c_1a69_00dd;

#[test]
fn time_golden_hash() {
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |v: f64| {
        for byte in v.to_bits().to_le_bytes() {
            acc ^= u64::from(byte);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    // Verlet: 10k pendulum steps.
    let force = |q: &[f64], out: &mut [f64]| out[0] = -fs_math::det::sin(q[0]);
    let (mut q, mut p) = (vec![1.2f64], vec![0.3f64]);
    let mut scratch = vec![0.0f64; 1];
    for _ in 0..10_000 {
        verlet_step(&mut q, &mut p, 0.01, &force, &mut scratch);
    }
    feed(q[0]);
    feed(p[0]);
    // Lie: asymmetric rigid body, 10k steps.
    let mut rotation = So3::identity();
    let mut angular = Vec3::new(0.3, 1.5, -0.4);
    for _ in 0..10_000 {
        let (next_rotation, next_angular) =
            rigid_body_step(rotation, angular, Vec3::new(1.0, 2.0, 3.0), 1e-3)
                .expect("finite golden rigid-body step");
        rotation = next_rotation;
        angular = next_angular;
    }
    let quat = rotation.as_quat();
    for v in [quat.w, quat.x, quat.y, quat.z] {
        feed(v);
    }
    for v in [angular.x, angular.y, angular.z] {
        feed(v);
    }
    // Generalized-α: damped 2-DOF, 1k steps.
    let m = [2.0f64, 0.0, 0.0, 1.0];
    let c = [0.1f64, 0.0, 0.0, 0.1];
    let k = [4.0f64, -1.0, -1.0, 3.0];
    let ga = GeneralizedAlpha::new(&m, &c, &k, 2, 0.05, 0.8);
    let (mut gq, mut gv, mut ga_acc) = (vec![1.0f64, -0.5], vec![0.0f64; 2], vec![0.0f64; 2]);
    // Consistent initial acceleration a = M⁻¹(f − Cv − Kq), f = 0.
    ga_acc[0] = -(k[0] * gq[0] + k[1] * gq[1]) / m[0];
    ga_acc[1] = -(k[2] * gq[0] + k[3] * gq[1]) / m[3];
    let f = vec![0.0f64; 2];
    for _ in 0..1_000 {
        galpha_step(&ga, &mut gq, &mut gv, &mut ga_acc, &f);
    }
    for v in gq.iter().chain(gv.iter()) {
        feed(*v);
    }
    // IMEX: 2-D stiff/nonstiff mix, 500 steps.
    let im = Imex2::new(&[-50.0, 1.0, 0.0, -0.5], 2, 0.02);
    let nl = |u: &[f64], out: &mut [f64]| {
        out[0] = -u[0] * u[1];
        out[1] = u[0] * u[0];
    };
    let mut iu = vec![1.0f64, 0.2];
    for _ in 0..500 {
        imex2_step(&im, &mut iu, &nl);
    }
    feed(iu[0]);
    feed(iu[1]);
    // Exponential Euler: symmetric A, 500 steps.
    let ee = ExpEuler::new(&[-2.0, 0.7, 0.7, -3.0], 2, 0.02);
    let mut eu = vec![0.9f64, -0.6];
    for _ in 0..500 {
        ee.step(&mut eu, &nl);
    }
    feed(eu[0]);
    feed(eu[1]);
    // RK45 + PI controller: full trajectory state including controller
    // memory and counters (the step SEQUENCE is part of the contract).
    let pi = PiController::default();
    let mut st = AdaptiveState::new(0.0, &[1.0, 0.0], 0.1);
    rk45_adaptive(&mut st, &ho_rhs, 10.0, 1e-9, 1e-12, &pi, 100_000);
    feed(st.u[0]);
    feed(st.u[1]);
    feed(st.h);
    feed(st.err_prev);
    feed(st.accepted as f64);
    feed(st.rejected as f64);
    log("time-golden", "info", &format!("{acc:#018x}"));
    assert_eq!(
        acc, GOLDEN_HASH,
        "fs-time bits changed: {acc:#018x} vs {GOLDEN_HASH:#018x} — bump only with semantic \
         justification (golden-evidence policy)"
    );
}
