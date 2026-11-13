use fs_spectral::service::{EigenBackend, EigenQuery, EigenService, ServiceError, SymmetricOp};
use std::cell::Cell;
struct Diagonal;
impl SymmetricOp for Diagonal {
    fn dim(&self) -> usize {
        24
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        for i in 0..24 {
            y[i] = (1.0 + i as f64) * x[i];
        }
    }
}
#[test]
fn callback_tick_cancels_transactionally_and_replays_the_same_native_state() {
    let query = EigenQuery {
        k: 2,
        largest: true,
        tol: 1e-10,
        steps_per_tick: 1,
    };
    let mut accepted = EigenService::new(EigenBackend::Lobpcg, 24, query).unwrap();
    accepted.tick_cancellable(&Diagonal, &|| false).unwrap();
    let mut expected = accepted.clone();
    let polls = Cell::new(0);
    assert!(matches!(
        accepted.tick_cancellable(&Diagonal, &|| {
            polls.set(polls.get() + 1);
            polls.get() >= 4
        }),
        Err(ServiceError::Cancelled)
    ));
    assert_eq!(accepted.ticks(), expected.ticks());
    let resumed = accepted.tick_cancellable(&Diagonal, &|| false).unwrap();
    let straight = expected.tick_cancellable(&Diagonal, &|| false).unwrap();
    assert_eq!(resumed.ticks, straight.ticks);
    for (a, b) in resumed.pairs.iter().zip(&straight.pairs) {
        assert_eq!(a.value().to_bits(), b.value().to_bits());
        assert_eq!(a.residual().to_bits(), b.residual().to_bits());
        assert_eq!(a.vector(), b.vector());
    }
}
/// `Diagonal` for its first `good` applications, then NaN: an operator whose own
/// work allowance runs out inside a tick writes NaN and keeps its refusal
/// (fs-solid's prestress modes do).
struct NanAfter {
    good: usize,
    calls: Cell<usize>,
}
impl SymmetricOp for NanAfter {
    fn dim(&self) -> usize {
        24
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        let call = self.calls.get();
        self.calls.set(call + 1);
        if call < self.good {
            Diagonal.apply(x, y);
        } else {
            y.fill(f64::NAN);
        }
    }
}
#[test]
fn callback_lobpcg_tick_refuses_a_non_finite_operator_and_rolls_back() {
    let query = EigenQuery {
        k: 2,
        largest: true,
        tol: 1e-10,
        steps_per_tick: 1,
    };
    let mut accepted = EigenService::new(EigenBackend::Lobpcg, 24, query).unwrap();
    accepted.tick_cancellable(&Diagonal, &|| false).unwrap();
    let mut expected = accepted.clone();
    // The current block's two actions are finite; the search space's are not, so
    // the Ritz block that reaches LOBPCG's QR is non-finite.
    let poisoned = NanAfter {
        good: 2,
        calls: Cell::new(0),
    };
    assert!(matches!(
        accepted.tick_cancellable(&poisoned, &|| false),
        Err(ServiceError::NonFiniteOperator)
    ));
    assert!(poisoned.calls.get() > 2);
    assert_eq!(accepted.ticks(), expected.ticks());
    let resumed = accepted.tick_cancellable(&Diagonal, &|| false).unwrap();
    let straight = expected.tick_cancellable(&Diagonal, &|| false).unwrap();
    assert_eq!(resumed.ticks, straight.ticks);
    for (a, b) in resumed.pairs.iter().zip(&straight.pairs) {
        assert_eq!(a.value().to_bits(), b.value().to_bits());
        assert_eq!(a.residual().to_bits(), b.residual().to_bits());
        assert_eq!(a.vector(), b.vector());
    }
}
