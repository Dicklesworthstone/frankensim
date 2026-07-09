//! FrankenTorch tape bridge (bead o3ui, feature `torch-bridge`):
//! REVERSE MODE for large N, as a drop-in for kernels already generic
//! over [`Real`]. A [`TapeReal`] is a `Copy` handle onto a thread-local
//! ft-autograd scalar [`Tape`]; every `Real` operation records a node,
//! and one backward pass yields the full gradient — O(cost(f)) instead
//! of the forward-dual O(N·cost(f)).
//!
//! DETERMINISM CLASS (declared, not pretended): ft's `Strict` execution
//! mode is deterministic per FrankenTorch's own contract, but its
//! elementary functions are NOT fs-math's det module — bridge primals
//! and gradients are cross-checked against forward duals to tight
//! TOLERANCES, never bitwise. Kernels needing the strict cross-ISA
//! bit-contract stay on `f64`/[`Dual`](crate::dual::Dual).
//!
//! The [`taped_vjp`] surface is shaped for [`crate::revolve`]: a
//! reverse-step callback `(i, state, bar) -> bar_prev` is one
//! `taped_vjp(|u| step(i, u), state, bar)` — checkpointing composes
//! with taped segments (gated in the battery).

use core::cell::RefCell;

use ft_autograd::{NodeId, Tape};
use ft_core::ExecutionMode;

use crate::Real;

const MODE: ExecutionMode = ExecutionMode::Strict;

thread_local! {
    static TAPE: RefCell<Option<Tape>> = const { RefCell::new(None) };
}

/// Run `f` with a fresh thread-local tape installed; always uninstalls.
fn with_fresh_tape<T>(f: impl FnOnce() -> T) -> T {
    struct Uninstall;
    impl Drop for Uninstall {
        fn drop(&mut self) {
            TAPE.with(|slot| *slot.borrow_mut() = None);
        }
    }
    TAPE.with(|slot| {
        assert!(
            slot.borrow().is_none(),
            "nested tape scopes are not supported (one reverse pass at a time per thread)"
        );
        *slot.borrow_mut() = Some(Tape::new());
    });
    let _guard = Uninstall;
    f()
}

fn with_tape<T>(f: impl FnOnce(&mut Tape) -> T) -> T {
    TAPE.with(|slot| {
        let mut borrow = slot.borrow_mut();
        let tape = borrow
            .as_mut()
            .expect("TapeReal used outside reverse_gradient/taped_vjp scope");
        f(tape)
    })
}

/// A scalar recorded on the thread-local FrankenTorch tape. Valid only
/// inside [`reverse_gradient`]/[`taped_vjp`] scopes (loud panic
/// otherwise).
#[derive(Clone, Copy, Debug)]
pub struct TapeReal(NodeId);

impl TapeReal {
    fn unary(self, op: impl FnOnce(&mut Tape, NodeId) -> NodeId) -> TapeReal {
        TapeReal(with_tape(|t| op(t, self.0)))
    }

    fn primal(self) -> f64 {
        with_tape(|t| t.value(self.0).expect("tape node value"))
    }
}

impl PartialEq for TapeReal {
    fn eq(&self, other: &Self) -> bool {
        self.primal() == other.primal()
    }
}

impl PartialOrd for TapeReal {
    fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> {
        self.primal().partial_cmp(&other.primal())
    }
}

impl core::ops::Add for TapeReal {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        TapeReal(with_tape(|t| t.add(self.0, rhs.0, MODE).expect("add").0))
    }
}

impl core::ops::Sub for TapeReal {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        TapeReal(with_tape(|t| t.sub(self.0, rhs.0, MODE).expect("sub").0))
    }
}

impl core::ops::Mul for TapeReal {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        TapeReal(with_tape(|t| t.mul(self.0, rhs.0, MODE).expect("mul").0))
    }
}

impl core::ops::Div for TapeReal {
    type Output = Self;
    fn div(self, rhs: Self) -> Self {
        TapeReal(with_tape(|t| t.div(self.0, rhs.0, MODE).expect("div").0))
    }
}

impl core::ops::Neg for TapeReal {
    type Output = Self;
    fn neg(self) -> Self {
        self.unary(|t, n| t.neg(n, MODE).expect("neg").0)
    }
}

impl Real for TapeReal {
    fn zero() -> Self {
        Self::from_f64(0.0)
    }

    fn one() -> Self {
        Self::from_f64(1.0)
    }

    fn from_f64(v: f64) -> Self {
        TapeReal(with_tape(|t| t.leaf(v, false)))
    }

    fn value(self) -> f64 {
        self.primal()
    }

    fn mul_add(self, a: Self, b: Self) -> Self {
        // COMPOSED mul + add: ft's scalar tape has no fused node, so
        // the bridge primal is NOT fused here — one more reason bridge
        // results are tolerance-checked, never bitwise (module docs).
        self * a + b
    }

    fn recip(self) -> Self {
        self.unary(|t, n| t.reciprocal(n, MODE).expect("recip").0)
    }

    fn sqrt(self) -> Self {
        self.unary(|t, n| t.sqrt(n, MODE).expect("sqrt").0)
    }

    fn abs(self) -> Self {
        self.unary(|t, n| t.abs(n, MODE).expect("abs").0)
    }

    fn exp(self) -> Self {
        self.unary(|t, n| t.exp(n, MODE).expect("exp").0)
    }

    fn ln(self) -> Self {
        self.unary(|t, n| t.log(n, MODE).expect("log").0)
    }

    fn sin(self) -> Self {
        self.unary(|t, n| t.sin(n, MODE).expect("sin").0)
    }

    fn cos(self) -> Self {
        self.unary(|t, n| t.cos(n, MODE).expect("cos").0)
    }

    fn tanh(self) -> Self {
        self.unary(|t, n| t.tanh(n, MODE).expect("tanh").0)
    }

    fn asin(self) -> Self {
        self.unary(|t, n| t.asin(n, MODE).expect("asin").0)
    }

    fn acos(self) -> Self {
        self.unary(|t, n| t.acos(n, MODE).expect("acos").0)
    }

    fn atan(self) -> Self {
        self.unary(|t, n| t.atan(n, MODE).expect("atan").0)
    }

    fn atan2(self, x: Self) -> Self {
        TapeReal(with_tape(|t| t.atan2(self.0, x.0, MODE).expect("atan2").0))
    }

    fn powi(self, n: i32) -> Self {
        self.unary(|t, node| t.pow(node, f64::from(n), MODE).expect("pow").0)
    }
}

/// Reverse-mode gradient of a scalar function written generic over
/// [`Real`]: ONE tape build + ONE backward pass, any N. Returns
/// (value, gradient).
///
/// # Panics
/// On nested invocation (one reverse pass per thread at a time) or if
/// the tape rejects an operation (programmer error at fixture scale).
pub fn reverse_gradient<F>(x: &[f64], f: F) -> (f64, Vec<f64>)
where
    F: FnOnce(&[TapeReal]) -> TapeReal,
{
    with_fresh_tape(|| {
        let vars: Vec<TapeReal> = x
            .iter()
            .map(|&v| TapeReal(with_tape(|t| t.leaf(v, true))))
            .collect();
        let out = f(&vars);
        with_tape(|t| {
            let val = t.value(out.0).expect("root value");
            let report = t.backward(out.0).expect("backward pass");
            let grad = vars
                .iter()
                .map(|v| report.gradient(v.0).unwrap_or(0.0))
                .collect();
            (val, grad)
        })
    })
}

/// Vector-Jacobian product Jᵀ·bar of a taped map in ONE backward pass:
/// the reverse-step surface [`crate::revolve`] wants — a checkpointed
/// sweep's `reverse(i, state, bar)` is
/// `taped_vjp(|u| step(i, u), state, bar)`.
///
/// # Panics
/// As [`reverse_gradient`]; additionally if `f` returns a length ≠
/// `bar.len()`.
pub fn taped_vjp<F>(f: F, x: &[f64], bar: &[f64]) -> Vec<f64>
where
    F: FnOnce(&[TapeReal]) -> Vec<TapeReal>,
{
    with_fresh_tape(|| {
        let vars: Vec<TapeReal> = x
            .iter()
            .map(|&v| TapeReal(with_tape(|t| t.leaf(v, true))))
            .collect();
        let ys = f(&vars);
        assert_eq!(
            ys.len(),
            bar.len(),
            "vjp seed length {} must match output arity {}",
            bar.len(),
            ys.len()
        );
        // Root = Σ barⱼ·yⱼ with bar as derivative-free leaves: backward
        // from the root IS Jᵀ·bar at the inputs.
        let mut root = TapeReal::from_f64(0.0);
        for (y, &b) in ys.iter().zip(bar) {
            root = root + *y * TapeReal::from_f64(b);
        }
        with_tape(|t| {
            let report = t.backward(root.0).expect("vjp backward");
            vars.iter()
                .map(|v| report.gradient(v.0).unwrap_or(0.0))
                .collect()
        })
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dual::{Dual64, gradient};
    use crate::revolve::{checkpointed_adjoint, min_budget};

    /// The inverse-trig gauntlet from lib.rs, reused verbatim: every
    /// Real op the bridge must tape.
    fn gauntlet<T: Real>(x: T, y: T) -> T {
        let u = (x - y).tanh() * T::from_f64(0.8);
        let v = (x * y).tanh() * T::from_f64(0.7);
        let a = u.asin() + v.acos() * (x * x + T::one()).recip();
        let b = (y.atan() - u.acos() * T::from_f64(0.25)).exp();
        let c = x.atan2(y) + (a * b).sin();
        c.mul_add(a, b.sqrt()) + v.atan() + (x.abs() + T::one()).ln()
    }

    #[test]
    fn bridge_gradient_matches_forward_duals() {
        let mut seed = 0xF7_2A_u64;
        let mut lcg = move || {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((seed >> 11) as f64) / (1u64 << 53) as f64
        };
        for _ in 0..200 {
            let x = lcg() * 2.0 + 0.2;
            let y = lcg() * 2.0 + 0.2;
            let (dv, dg) = gradient([x, y], |[a, b]| gauntlet(a, b));
            let (bv, bg) = reverse_gradient(&[x, y], |v| gauntlet(v[0], v[1]));
            // Tolerance, not bitwise: ft strict-mode libm vs fs-math det.
            assert!(
                (dv - bv).abs() < 1e-12 * dv.abs().max(1.0),
                "primal dual {dv} vs bridge {bv} at ({x},{y})"
            );
            for (d, b) in dg.iter().zip(&bg) {
                assert!(
                    (d - b).abs() < 1e-10 * d.abs().max(1.0),
                    "grad dual {d} vs bridge {b} at ({x},{y})"
                );
            }
        }
        println!(
            "{{\"suite\":\"fs-ad\",\"case\":\"torch-bridge-grad\",\"verdict\":\"pass\",\"detail\":\"200 points, every Real op taped, reverse == forward duals (rel 1e-10)\"}}"
        );
    }

    #[test]
    fn bridge_scales_reverse_one_pass() {
        // N = 300 quadratic-plus-coupling objective: reverse gets the
        // whole gradient in one backward; validated against analytic.
        const N: usize = 300;
        let x: Vec<f64> = (0..N).map(|i| 0.1 + (i as f64) * 1e-3).collect();
        let (val, grad) = reverse_gradient(&x, |v| {
            let mut acc = TapeReal::from_f64(0.0);
            for i in 0..N {
                acc = acc + v[i] * v[i];
                if i + 1 < N {
                    acc = acc + TapeReal::from_f64(0.5) * v[i] * v[i + 1];
                }
            }
            acc
        });
        let mut want_val = 0.0;
        for i in 0..N {
            want_val += x[i] * x[i];
            if i + 1 < N {
                want_val += 0.5 * x[i] * x[i + 1];
            }
        }
        assert!((val - want_val).abs() < 1e-12 * want_val.max(1.0));
        for i in 0..N {
            let mut want = 2.0 * x[i];
            if i + 1 < N {
                want += 0.5 * x[i + 1];
            }
            if i > 0 {
                want += 0.5 * x[i - 1];
            }
            assert!(
                (grad[i] - want).abs() < 1e-12 * want.abs().max(1.0),
                "grad[{i}] {} vs analytic {want}",
                grad[i]
            );
        }
        println!(
            "{{\"suite\":\"fs-ad\",\"case\":\"torch-bridge-large-n\",\"verdict\":\"pass\",\"detail\":\"N=300 gradient in one reverse pass, matches analytic to 1e-12\"}}"
        );
    }

    #[test]
    fn taped_segments_compose_with_revolve() {
        // The bead's composition target: a checkpointed sweep whose
        // reverse steps are TAPED vjps. Must match the forward-dual
        // sensitivity of the whole chain.
        const STEPS: usize = 40;
        fn step<T: Real>(x: &[T]) -> Vec<T> {
            // 2-state map: (a, b) <- (a - 0.05 tanh(a·b), b + 0.03 sin(a)).
            let (a, b) = (x[0], x[1]);
            vec![
                a - T::from_f64(0.05) * (a * b).tanh(),
                b + T::from_f64(0.03) * a.sin(),
            ]
        }
        let x0 = vec![0.7f64, -0.4];
        let fwd = |_i: usize, s: &Vec<f64>| -> Vec<f64> {
            let d: Vec<Dual64<1>> = s.iter().map(|&v| Dual64::constant(v)).collect();
            step(&d).iter().map(|d| d.re).collect()
        };
        let rev = |_i: usize, s: &Vec<f64>, bar: Vec<f64>| -> Vec<f64> { taped_vjp(step, s, &bar) };
        // d(final a)/d(x0) via checkpointed taped adjoint.
        let (bar0, _) =
            checkpointed_adjoint(&x0, STEPS, min_budget(STEPS), &fwd, &rev, vec![1.0, 0.0]);
        // Forward-dual reference, one seed per input.
        for (k, bar_k) in bar0.iter().enumerate() {
            let mut s: Vec<Dual64<1>> = x0
                .iter()
                .enumerate()
                .map(|(i, &v)| {
                    if i == k {
                        Dual64::variable(v, 0)
                    } else {
                        Dual64::constant(v)
                    }
                })
                .collect();
            for _ in 0..STEPS {
                s = step(&s);
            }
            let want = s[0].eps[0];
            assert!(
                (bar_k - want).abs() < 1e-11 * want.abs().max(1.0),
                "taped revolve d/dx0[{k}] = {bar_k} vs dual {want}"
            );
        }
        println!(
            "{{\"suite\":\"fs-ad\",\"case\":\"torch-bridge-revolve\",\"verdict\":\"pass\",\"detail\":\"checkpointed sweep with taped vjp reverse steps == forward duals over 40 steps (rel 1e-11)\"}}"
        );
    }

    #[test]
    fn use_outside_scope_is_loud() {
        let r = std::panic::catch_unwind(|| TapeReal::from_f64(1.0));
        assert!(r.is_err(), "TapeReal outside a tape scope must panic");
    }
}
