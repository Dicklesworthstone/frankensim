//! fs-ad — forward-mode automatic differentiation (plan §6.6 regime 1) and
//! the [`Real`] generic-scalar contract.
//!
//! THE CONTRACT THAT MATTERS: physics operators written generic over
//! [`Real`] run unchanged on `f64` (production) and on [`dual::Dual`]
//! (derivative stress-testing) — which is how "stress-test the adjoint
//! against forward duals" becomes a one-liner in the gradient gate
//! (plan §8.7: a solver without a passing gradient check cannot merge).
//!
//! Determinism: the `f64` implementation of `Real` routes elementary
//! functions through fs-math's STRICT det module, so making code generic
//! over `Real` cannot silently reintroduce platform libm — genericity and
//! cross-ISA bit-determinism compose instead of fighting.
//!
//! Reverse mode (FrankenTorch tape bridge) and adjoint-by-construction (IFT,
//! revolve checkpointing) are the fs-ad-adjoint-infra bead's scope.

pub mod dual;

pub use dual::{Dual, Dual64, gradient, jvp, second_directional};

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// The generic scalar: everything a FrankenSim kernel may ask of a number.
/// Implemented by `f64` (strict-mode via fs-math) and by `Dual<T, N>` for
/// any `T: Real` (hence nested duals = higher-order derivatives).
///
/// Design notes:
/// - Deliberately SMALL: only operations with well-defined derivative rules
///   belong here. Comparisons/branching use `value()` (the primal), which is
///   the standard forward-AD convention for control flow — DOCUMENTED
///   caveat: branches are differentiated per-branch (kinks give one-sided
///   derivatives).
/// - `powi` over `powf`: integer powers have exact derivative rules; general
///   `pow` arrives with fs-math extensions.
pub trait Real:
    Copy
    + core::fmt::Debug
    + PartialOrd
    + core::ops::Add<Output = Self>
    + core::ops::Sub<Output = Self>
    + core::ops::Mul<Output = Self>
    + core::ops::Div<Output = Self>
    + core::ops::Neg<Output = Self>
{
    /// Additive identity.
    fn zero() -> Self;
    /// Multiplicative identity.
    fn one() -> Self;
    /// Lift a constant (derivative-free) value.
    fn from_f64(v: f64) -> Self;
    /// The primal value (used for branching/comparison in generic code).
    fn value(self) -> f64;
    /// Fused multiply-add: `self * a + b`.
    #[must_use]
    fn mul_add(self, a: Self, b: Self) -> Self;
    /// Reciprocal.
    #[must_use]
    fn recip(self) -> Self;
    /// Square root.
    #[must_use]
    fn sqrt(self) -> Self;
    /// Absolute value (subgradient 0 at the kink — documented convention).
    #[must_use]
    fn abs(self) -> Self;
    /// e^x (strict-mode on f64).
    #[must_use]
    fn exp(self) -> Self;
    /// Natural log (strict-mode on f64).
    #[must_use]
    fn ln(self) -> Self;
    /// sin (strict-mode on f64; fs-math domain notes apply).
    #[must_use]
    fn sin(self) -> Self;
    /// cos (strict-mode on f64; fs-math domain notes apply).
    #[must_use]
    fn cos(self) -> Self;
    /// tanh (strict-mode on f64).
    #[must_use]
    fn tanh(self) -> Self;
    /// Integer power.
    #[must_use]
    fn powi(self, n: i32) -> Self;
}

impl Real for f64 {
    fn zero() -> Self {
        0.0
    }

    fn one() -> Self {
        1.0
    }

    fn from_f64(v: f64) -> Self {
        v
    }

    fn value(self) -> f64 {
        self
    }

    fn mul_add(self, a: Self, b: Self) -> Self {
        f64::mul_add(self, a, b)
    }

    fn recip(self) -> Self {
        1.0 / self
    }

    fn sqrt(self) -> Self {
        fs_math::det::sqrt(self)
    }

    fn abs(self) -> Self {
        f64::abs(self)
    }

    fn exp(self) -> Self {
        fs_math::det::exp(self)
    }

    fn ln(self) -> Self {
        fs_math::det::ln(self)
    }

    fn sin(self) -> Self {
        fs_math::det::sin(self)
    }

    fn cos(self) -> Self {
        fs_math::det::cos(self)
    }

    fn tanh(self) -> Self {
        fs_math::det::tanh(self)
    }

    fn powi(self, n: i32) -> Self {
        f64::powi(self, n)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lcg(seed: &mut u64) -> f64 {
        *seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((*seed >> 11) as f64) / (1u64 << 53) as f64
    }

    /// A nasty composite exercising every Real operation with chain/product/
    /// quotient rules stacked three deep.
    fn gauntlet<T: Real>(x: T, y: T) -> T {
        let a = (x * y).sin() + (x - y).tanh() * (x * x + T::one()).ln();
        let b = (a.abs() + T::from_f64(0.5)).sqrt().exp();
        let c = (b / (y * y + T::one())).cos();
        c.mul_add(a, b.recip()) + x.powi(3) * y.powi(-2)
    }

    #[test]
    fn dual_gradient_matches_central_differences() {
        let mut seed = 0xAD_F00D_u64;
        let h = 1e-6;
        for _ in 0..500 {
            let x = lcg(&mut seed) * 2.0 + 0.2;
            let y = lcg(&mut seed) * 2.0 + 0.2;
            let (_, grad) = gradient([x, y], |[dx, dy]| gauntlet(dx, dy));
            // Central differences on the SAME strict-mode f64 path.
            let gx = (gauntlet(x + h, y) - gauntlet(x - h, y)) / (2.0 * h);
            let gy = (gauntlet(x, y + h) - gauntlet(x, y - h)) / (2.0 * h);
            for (ad, fd) in [(grad[0], gx), (grad[1], gy)] {
                let scale = ad.abs().max(fd.abs()).max(1.0);
                assert!(
                    (ad - fd).abs() / scale < 5e-9,
                    "dual {ad} vs FD {fd} at ({x},{y})"
                );
            }
        }
        println!(
            "{{\"suite\":\"fs-ad\",\"case\":\"grad-vs-fd\",\"verdict\":\"pass\",\"detail\":\"500 random points, rel<5e-9\"}}"
        );
    }

    #[test]
    fn primal_path_is_bitwise_identical_to_scalar() {
        // Running through Dual must not perturb the VALUE channel: same
        // strict functions, same order — bit-identical primal.
        let mut seed = 7u64;
        for _ in 0..2_000 {
            let x = lcg(&mut seed) * 3.0 + 0.1;
            let y = lcg(&mut seed) * 3.0 + 0.1;
            let scalar = gauntlet(x, y);
            let (dual_val, _) = gradient([x, y], |[dx, dy]| gauntlet(dx, dy));
            assert_eq!(
                scalar.to_bits(),
                dual_val.to_bits(),
                "dual perturbed the primal at ({x},{y})"
            );
        }
    }

    #[test]
    fn known_analytic_derivatives() {
        // f(x) = sin(x²): f' = 2x·cos(x²), f'' = 2cos(x²) − 4x²·sin(x²).
        let f = |v: [Dual<Dual64<1>, 1>; 1]| (v[0] * v[0]).sin();
        for x in [0.3, 1.1, 2.7, -0.9] {
            let (val, d1, d2) = second_directional([x], [1.0], f);
            let z = x * x;
            assert!((val - fs_math::det::sin(z)).abs() < 1e-15);
            let want1 = 2.0 * x * fs_math::det::cos(z);
            let want2 = 2.0 * fs_math::det::cos(z) - 4.0 * z * fs_math::det::sin(z);
            assert!((d1 - want1).abs() < 1e-13, "f' {d1} vs {want1} at {x}");
            assert!((d2 - want2).abs() < 1e-12, "f'' {d2} vs {want2} at {x}");
        }
    }

    #[test]
    fn jvp_is_directional_derivative() {
        let f = |v: [Dual64<1>; 3]| v[0] * v[1] + v[2].exp() * v[0].sin();
        let x = [0.7, -1.2, 0.4];
        let dir = [0.5, 2.0, -1.0];
        let (_, dv) = jvp(x, dir, f);
        // Compare against the full gradient contracted with the direction.
        let (_, g) = gradient(x, |[a, b, c]| a * b + c.exp() * a.sin());
        let want: f64 = g.iter().zip(&dir).map(|(gi, di)| gi * di).sum();
        assert!((dv - want).abs() < 1e-14, "jvp {dv} vs grad·v {want}");
    }

    #[test]
    fn generic_newton_differentiates_through_convergence() {
        // Solve u³ = c by Newton, GENERIC over Real; differentiating through
        // the converged iteration must give d(c^(1/3))/dc = (1/3)c^(-2/3).
        fn cbrt_newton<T: Real>(c: T) -> T {
            let mut u = T::one();
            for _ in 0..60 {
                u = u - (u.powi(3) - c) / (T::from_f64(3.0) * u.powi(2));
            }
            u
        }
        for c in [0.5, 2.0, 8.0, 27.0] {
            let (val, grad) = gradient([c], |[d]| cbrt_newton(d));
            let want_val = c.powf(1.0 / 3.0);
            let want_grad = (1.0 / 3.0) * c.powf(-2.0 / 3.0);
            assert!((val - want_val).abs() < 1e-12, "cbrt {val} vs {want_val}");
            assert!(
                (grad[0] - want_grad).abs() / want_grad < 1e-10,
                "d/dc {} vs {want_grad} at c={c}",
                grad[0]
            );
        }
        println!(
            "{{\"suite\":\"fs-ad\",\"case\":\"generic-solver\",\"verdict\":\"pass\",\"detail\":\"Newton cbrt differentiated through convergence\"}}"
        );
    }

    #[test]
    fn abs_kink_and_sqrt_zero_conventions() {
        let (_, g) = gradient([0.0], |[x]| x.abs());
        assert_eq!(
            g[0].to_bits(),
            0.0f64.to_bits(),
            "abs'(0) = 0 by documented convention"
        );
        let (_, g) = gradient([2.0], |[x]| x.abs());
        assert_eq!(g[0].to_bits(), 1.0f64.to_bits());
        let (_, g) = gradient([-2.0], |[x]| x.abs());
        assert_eq!(g[0].to_bits(), (-1.0f64).to_bits());
        // sqrt at 0: derivative is honestly unbounded (inf), never silently
        // clamped.
        let (_, g) = gradient([0.0], |[x]| x.sqrt());
        assert!(g[0].is_infinite(), "sqrt'(0) must be inf, got {}", g[0]);
    }

    #[test]
    fn multi_lane_gradients_agree_with_single_lane() {
        // Dual<4> computing 4 partials at once must equal 4 × Dual<1> runs.
        let f4 =
            |v: [Dual64<4>; 4]| (v[0] * v[1]).sin() + (v[2] / v[3]).exp() + v[0].tanh() * v[3].ln();
        let x = [0.9, 1.7, 0.3, 2.2];
        let (_, g4) = gradient(x, f4);
        for i in 0..4 {
            let mut dir = [0.0; 4];
            dir[i] = 1.0;
            let (_, gi) = jvp(x, dir, |v: [Dual64<1>; 4]| {
                (v[0] * v[1]).sin() + (v[2] / v[3]).exp() + v[0].tanh() * v[3].ln()
            });
            assert_eq!(
                g4[i].to_bits(),
                gi.to_bits(),
                "lane {i}: packed {} vs single {gi}",
                g4[i]
            );
        }
    }

    #[test]
    fn version_is_stamped() {
        assert!(!VERSION.is_empty());
    }
}
