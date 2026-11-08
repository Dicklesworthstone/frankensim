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
