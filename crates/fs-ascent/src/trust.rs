//! Trust-region Newton–Krylov with a reusable numerical engine and a Cx adapter.
//! The engine owns accepted state, radius policy, callback ceilings and restart;
//! this facade adds cancellation without changing the numerical operation tree.

mod engine;
pub use engine::{
    PreconditionerError, TrustRegionProgress, TrustRegionReport, TrustRegionRunReport, TrustRegionState,
    hv_fd_of_gradients, trust_region_newton,
};
use fs_exec::Cx;
use crate::stop::StopRule;

impl TrustRegionState {
    /// As [`Self::run`], polling before and after every Hessian-vector and
    /// objective/gradient callback, including inside Steihaug-CG. Individual
    /// callbacks and vector operations are not interruptible. Cancellation
    /// discards the uncommitted trial, preserving the accepted point, radius,
    /// history and iteration count. Spent callbacks remain charged. Resumption
    /// with a fresh context restarts that trial from the accepted point; it
    /// does not replay earlier accepted iterations or retain the Krylov basis.
    pub fn run_cancellable(
        &mut self,
        fg: crate::FnGrad<'_>,
        hv_at: crate::FnHv<'_>,
        rule: &StopRule,
        max_iters: usize,
        cx: &Cx,
    ) -> TrustRegionRunReport {
        self.run_with_pause(fg, hv_at, rule, max_iters, &mut || {
            cx.checkpoint().is_err()
        })
    }

    /// Cancellable form of [`Self::run_bounded`]. Cancellation is checked
    /// before budget attribution and immediately after every attempted
    /// callback. This bounds callback counts, not callback-internal work,
    /// allocation or elapsed time; those remain the model/context's budgets.
    pub fn run_bounded_cancellable(
        &mut self,
        fg: crate::FnGrad<'_>,
        hv_at: crate::FnHv<'_>,
        rule: &StopRule,
        max_iters: usize,
        max_hessian_products: usize,
        cx: &Cx,
    ) -> TrustRegionRunReport {
        self.run_with_limits(fg, hv_at, rule, max_iters, max_hessian_products,
            &mut || cx.checkpoint().is_err())
    }

    /// Cancellable preconditioned Newton--Krylov. Cancellation is polled around
    /// every inverse action as well as the objective and Hessian callbacks.
    /// Both derivative-work ceilings and their spent counts survive a pause.
    /// See [`Self::run_preconditioned_bounded`] for the fixed-SPD contract.
    #[allow(clippy::too_many_arguments)]
    pub fn run_preconditioned_bounded_cancellable(
        &mut self,
        fg: crate::FnGrad<'_>,
        hv_at: crate::FnHv<'_>,
        inverse: crate::FnHv<'_>,
        rule: &StopRule,
        max_iters: usize,
        max_hessian_products: usize,
        max_preconditioner_applications: usize,
        cx: &Cx,
    ) -> Result<TrustRegionRunReport, PreconditionerError> {
        self.run_with_preconditioner(fg, hv_at, Some(inverse), rule, max_iters,
            max_hessian_products, max_preconditioner_applications,
            &mut || cx.checkpoint().is_err())
    }

}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::StopReason;
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
    use std::cell::Cell;

    // G4: cover both the legacy and bounded Cx adapters with cancellation
    // requested from each actual callback, not a substituted polling object.
    #[test]
    fn callback_cancellation_preserves_the_incumbent_and_spent_allowance() {
        for bounded in [false, true] {
            for request_at in 1..=3 {
                let gate = CancelGate::new();
                let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
                let calls = Cell::new(0usize);
                let mut fg = |x: &[f64]| ((x[0] - 2.0).powi(2), vec![2.0 * (x[0] - 2.0)]);
                let mut state = TrustRegionState::new(&[0.0], &mut fg);
                let observe = || {
                    assert!(!gate.is_requested(), "callback ran after cancellation");
                    calls.set(calls.get() + 1);
                    if calls.get() == request_at { gate.request(); }
                };
                pool.scope(|arena| {
                    let cx = Cx::new(&gate, arena,
                        StreamKey { seed: 0, kernel_id: 1, tile: 0, iteration: 0 },
                        Budget::INFINITE, ExecMode::Deterministic);
                    let mut sample = |x: &[f64]| { observe(); fg(x) };
                    let mut product = |_: &[f64], v: &[f64]| { observe(); vec![2.0 * v[0]] };
                    let rule = StopRule::Budget(2);
                    let result = if bounded {
                        state.run_bounded_cancellable(&mut sample, &mut product, &rule, 10, 2, &cx)
                    } else {
                        state.run_cancellable(&mut sample, &mut product, &rule, 10, &cx)
                    };
                    assert_eq!(result.progress, TrustRegionProgress::Paused);
                });
                assert_eq!(calls.get(), request_at);
                let report = state.report();
                assert_eq!(report.x, vec![0.0]);
                assert_eq!(report.f, 4.0);
                assert_eq!(report.iters, 0);
                assert_eq!(state.gradient(), &[-4.0]);
                assert_eq!(state.radius(), 1.0);
                assert_eq!(state.history(), &[4.0]);
                assert_eq!(report.hv_evals, request_at.min(2));
                assert_eq!(report.evals, 1 + usize::from(request_at == 3));
                if request_at == 3 {
                    let result = state.run_bounded(
                        &mut |_| panic!("spent objective allowance reissued"),
                        &mut |_, _| panic!("spent allowance reached Hessian"),
                        &StopRule::Budget(2), 10, 2,
                    );
                    assert_eq!(result.progress, TrustRegionProgress::Stopped(StopReason::Budget));
                    assert_eq!(state.report().iters, 0);
                    assert_eq!(state.radius(), 1.0);
                    assert_eq!(state.history(), &[4.0]);
                }
            }
        }
    }
}
