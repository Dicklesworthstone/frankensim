//! Sparse cancellation must reach the work inside a solver's operator callback.

use fs_solver::{CsrOp, FgmresState, LinearOp, OperatorApplyError, SolverRunProgress};
use fs_sparse::{APPLY_CHECKPOINT_WORK, Coo, precond::IdentityPrecond};
use std::cell::Cell;

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|value| value.to_bits()).collect()
}

struct ObservedCsr {
    operator: CsrOp,
    polls: Cell<usize>,
}

impl LinearOp for ObservedCsr {
    fn n(&self) -> usize {
        self.operator.n()
    }

    fn apply(&self, x: &[f64], output: &mut [f64]) {
        self.operator.apply(x, output);
    }

    fn apply_transpose(&self, x: &[f64], output: &mut [f64]) {
        self.operator.apply_transpose(x, output);
    }

    fn apply_controlled(
        &self,
        x: &[f64],
        output: &mut [f64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<(), OperatorApplyError> {
        self.operator.apply_controlled(x, output, &mut || {
            self.polls.set(self.polls.get() + 1);
            cancelled()
        })
    }
}

#[test]
fn fgmres_cancels_inside_a_wide_sparse_row_and_resumes_bitwise() {
    let n = 2 * APPLY_CHECKPOINT_WORK + 1;
    let mut coo = Coo::new(n, n);
    coo.push(0, 0, 4.0);
    for column in 1..n {
        coo.push(0, column, 0.001);
        coo.push(column, column, 3.0);
    }
    let operator = ObservedCsr {
        operator: CsrOp::general(coo.assemble()),
        polls: Cell::new(0),
    };
    let rhs = vec![1.0; n];
    let initial = FgmresState::new(&rhs, 3);
    let mut expected = initial.clone();
    let expected_report = expected.run(&operator, &IdentityPrecond, &rhs, 1e-12, 4);
    assert!(expected_report.converged_euclidean());
    let mut state = initial.clone();
    // Entry + staged-output initialization + pre-kernel poll have finished.
    // This next checkpoint is in the first row's arithmetic, before it ends.
    let cancel_at = n.div_ceil(APPLY_CHECKPOINT_WORK) + 3;
    let outcome = state.run_controlled(&operator, &IdentityPrecond, &rhs, 1e-12, 4, &mut || {
        operator.polls.get() == cancel_at
    });
    assert_eq!(outcome.progress, SolverRunProgress::Paused);
    assert_eq!(operator.polls.get(), cancel_at);
    assert_eq!(outcome.callbacks.operator, 1);
    assert_eq!(outcome.callbacks.preconditioner, 0);
    assert_eq!(bits(&state.x), bits(&initial.x));
    assert_eq!(state.iters, initial.iters);
    assert_eq!(bits(&state.history), bits(&initial.history));
    assert_eq!(
        state.rel_residual().to_bits(),
        initial.rel_residual().to_bits()
    );

    let resumed = state.run_controlled(&operator, &IdentityPrecond, &rhs, 1e-12, 4, &mut || false);
    assert_eq!(resumed.progress, SolverRunProgress::Complete);
    assert!(resumed.report.converged_euclidean());
    assert_eq!(bits(&state.x), bits(&expected.x));
    assert_eq!(bits(&state.history), bits(&expected.history));
    assert_eq!(state.iters, expected.iters);
}

#[test]
fn controlled_operator_faults_are_not_reported_as_cancellation() {
    struct NoWorkspace;
    impl LinearOp for NoWorkspace {
        fn n(&self) -> usize {
            2
        }
        fn apply(&self, _: &[f64], _: &mut [f64]) {
            panic!("controlled solver skipped the fallible operator hook");
        }
        // Never called (apply panics); the transpose is the same refusal.
        fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
            self.apply(x, y);
        }
        fn apply_controlled(
            &self,
            _: &[f64],
            output: &mut [f64],
            _: &mut dyn FnMut() -> bool,
        ) -> Result<(), OperatorApplyError> {
            output[0] = 100.0;
            Err(OperatorApplyError::Allocation)
        }
    }
    let mut state = FgmresState::new(&[1.0, 2.0], 2);
    let initial = state.clone();
    let outcome = state.run_controlled(
        &NoWorkspace,
        &IdentityPrecond,
        &[1.0, 2.0],
        1e-12,
        2,
        &mut || false,
    );
    assert_eq!(
        outcome.progress,
        SolverRunProgress::OperatorFailed(OperatorApplyError::Allocation)
    );
    assert_eq!(outcome.callbacks.operator, 1);
    assert_eq!(bits(&state.x), bits(&initial.x));
    assert_eq!(state.iters, initial.iters);
    assert_eq!(bits(&state.history), bits(&initial.history));

    let rectangular = CsrOp::general(Coo::new(2, 3).assemble());
    let outcome = state.run_controlled(
        &rectangular,
        &IdentityPrecond,
        &[1.0, 2.0],
        1e-12,
        2,
        &mut || false,
    );
    assert_eq!(
        outcome.progress,
        SolverRunProgress::OperatorFailed(OperatorApplyError::DimensionMismatch {
            buffer: "input",
            expected: 3,
            actual: 2
        },)
    );
    assert_eq!(bits(&state.x), bits(&initial.x));
    assert_eq!(state.iters, initial.iters);
}

#[test]
fn controlled_csr_transpose_is_object_safe_and_transactional() {
    let mut coo = Coo::new(2, 2);
    coo.push(0, 0, 2.0);
    coo.push(0, 1, -3.0);
    coo.push(1, 0, 7.0);
    coo.push(1, 1, 4.0);
    let csr = CsrOp::general(coo.assemble());
    let operator: &dyn LinearOp = &csr;
    let input = [1.5, -2.0];
    let mut expected = [0.0; 2];
    operator.apply_transpose(&input, &mut expected);
    let mut output = [17.0, -0.0];
    let before = bits(&output);
    assert_eq!(
        operator.apply_transpose_controlled(&input, &mut output, &mut || true),
        Err(OperatorApplyError::Interrupted),
    );
    assert_eq!(bits(&output), before);
    operator
        .apply_transpose_controlled(&input, &mut output, &mut || false)
        .unwrap();
    assert_eq!(bits(&output), bits(&expected));
}

#[test]
fn a_predicate_panic_inside_sparse_work_keeps_its_original_unwind() {
    let operator = ObservedCsr {
        operator: CsrOp::symmetric(Coo::new(1, 1).assemble()),
        polls: Cell::new(0),
    };
    let mut state = FgmresState::new(&[1.0], 1);
    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        state.run_controlled(&operator, &IdentityPrecond, &[1.0], 1e-12, 1, &mut || {
            assert_eq!(operator.polls.get(), 0, "predicate failure");
            false
        })
    }));
    assert!(
        result.is_err(),
        "predicate panic became a numerical callback report"
    );
    assert_eq!(state.x, [0.0]);
    assert_eq!(state.iters, 0);
    assert!(state.history.is_empty());
}
