//! fs-solver — the solver stack (plan §8.9). Layer: L3 FLUX.
//!
//! Matrix-free Krylov methods (CG, MINRES, GMRES(m)) bound by the four
//! workspace contract obligations from day one:
//!
//! - RESUMABLE: solver state is plain data — `clone()` is a
//!   checkpoint, and split runs are bitwise-equal to straight runs
//!   (the fs-time `AdaptiveState` pattern, tested the same way).
//! - CANCELLABLE: iteration granularity — every state is complete
//!   between iterations. FGMRES and Newton additionally accept a Cx,
//!   poll numerical callbacks, and discard interrupted staged work.
//! - DETERMINISTIC: all inner products go through the fixed-shape
//!   chunked reduction (fs-tilelang's combiner — shape depends on
//!   length only, never on threads or tiers).
//! - ADJOINT-EQUIPPED: the operator trait carries `apply_transpose`,
//!   and transposed solves run through the same machinery (tested to
//!   converge comparably to primal — the IFT contract's enabler).
//!
//! Error transparency: every solve returns a residual HISTORY and a
//! structured stall diagnosis instead of a timeout mystery.

pub mod block;
pub mod checked_cg;
/// Source-refusing numerical Jacobian actions for the existing nonlinear owner.
pub mod finite_difference;
/// Outward residual evaluation and checked discrete algebraic goal bounds.
pub mod goal;
pub mod krylov;
pub mod least_squares;
pub mod lsqr;
pub mod mixed;
pub mod nonlinear;
pub mod op;
/// Bounded PSD nonnegative quadratic baseline; no adjoint/uniqueness claim.
pub mod orthant;
pub mod pmg;
pub mod stokes;

pub use block::{
    BlockError, BlockOperator, BlockOperator2, BlockOperator3, BlockSchur2,
    RealEquivalentComplexOp, RectLinearOp, SchurSolveSign, SquareBlock, ZeroBlock,
};
pub use checked_cg::{CheckedCgConfig, CheckedCgError, CheckedCgSolution, checked_cg};
pub use krylov::{
    CgState, GmresState, MinresState, PminresState, ResidualClaim, SolveReport, StallDiagnosis,
};
pub use mixed::{CsrF32, MixedReport, mixed_cg_refine};
pub use nonlinear::{
    AdmittedLinearSolver, CancellableSolveReport, DefinitenessEvidence, FgmresState,
    FlexiblePreconditioner, Globalization, GlobalizationDecision, LineSearchConfig,
    LinearSolverKind, LinearSystemFinding, LinearSystemVerifier, LinearVerificationError,
    NewtonError, NewtonIteration, NewtonKrylovConfig, NewtonKrylovState, NewtonReport,
    NewtonStallDiagnosis, NonlinearProblem, NullspaceEvidence, PreconditionerClass,
    SolverAdmissionError, SolverCallback, SolverCallbackCounts, SolverRunProgress,
    SourceCompatibility, SymmetryEvidence, TrustRegionConfig, VerifiedLinearSystem,
    admit_linear_solver, verify_linear_system,
};
pub use op::{CsrOp, LinearOp};
pub use pmg::{MaskedTensorOp, PMultigrid};
pub use stokes::{StokesBlockDiag, StokesOp, StokesSystem};

/// The single L1 spectral service consumed by solver clients. This re-export
/// prevents a second eigensolver/health namespace from growing in fs-solver.
pub use fs_spectral as spectral_service;

/// Deterministic inner product: elementwise products folded through
/// the fixed-shape chunked combiner (shape = f(length) only, never on threads).
#[must_use]
pub fn dot(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len(), "dot length mismatch");
    let prods: Vec<f64> = a.iter().zip(b).map(|(x, y)| x * y).collect();
    fs_tilelang::deterministic_sum(&prods)
}

/// Deterministic 2-norm, with amplitude scaling before the squared reduction.
///
/// Power-of-two scaling retains the fixed 64-entry reduction order and needs
/// no allocated workspace. Squared accumulation cannot overflow; final
/// rescaling still reports an unrepresentable norm as infinity. A nonzero
/// finite vector cannot acquire a zero norm just because its squares underflow.
#[must_use]
pub fn norm2(a: &[f64]) -> f64 {
    let mut maximum = 0.0_f64;
    for &value in a {
        // f64::max ignores NaNs. A callback's invalid output must propagate
        // to solver admission instead of disappearing from its residual norm.
        if value.is_nan() {
            return f64::NAN;
        }
        maximum = maximum.max(value.abs());
    }
    if maximum == 0.0 || maximum.is_infinite() {
        return maximum;
    }
    // Remove the significand to obtain a normal power of two. Subnormal
    // vectors use MIN_POSITIVE: even the smallest entry then scales to 2^-52,
    // whose square is representable. The reciprocal is exact and finite.
    let scale = f64::from_bits(maximum.to_bits() & (0x7ff_u64 << 52)).max(f64::MIN_POSITIVE);
    let inverse_scale = scale.recip();
    let mut total = 0.0;
    for chunk in a.chunks(fs_tilelang::REDUCTION_CHUNK) {
        let mut partial = 0.0;
        for &value in chunk {
            let normalized = value * inverse_scale;
            partial += normalized * normalized;
        }
        total += partial;
    }
    scale * fs_math::det::sqrt(total)
}

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

#[cfg(test)]
mod norm_tests {
    use super::*;

    #[test]
    fn norm_preserves_nonzero_amplitudes_without_overflowing_their_squares() {
        for amplitude in [1e-300, 1e-200, 1e-9, 1.0, 1e200, 1e300] {
            let actual = norm2(&[3.0 * amplitude, -4.0 * amplitude]);
            assert!((actual / amplitude - 5.0).abs() < 2e-15);
        }
        let quantum = f64::from_bits(1);
        assert_eq!(
            norm2(&[3.0 * quantum, 4.0 * quantum]).to_bits(),
            (5.0 * quantum).to_bits()
        );
        assert_eq!(norm2(&[f64::MAX]), f64::MAX);
    }

    #[test]
    fn norm_preserves_zero_and_distinguishes_invalid_and_unrepresentable_values() {
        assert_eq!(norm2(&[]).to_bits(), 0.0_f64.to_bits());
        assert_eq!(norm2(&[0.0, -0.0]).to_bits(), 0.0_f64.to_bits());
        assert_eq!(norm2(&[f64::INFINITY, 1.0]), f64::INFINITY);
        assert_eq!(norm2(&[f64::MAX, f64::MAX]), f64::INFINITY);
        for values in [[f64::NAN, 0.0], [f64::INFINITY, f64::NAN]] {
            assert!(norm2(&values).is_nan());
        }
    }

    #[test]
    fn norm_retains_the_original_chunk_order_at_ordinary_scales() {
        for count in [1, 63, 64, 65, 129, 1025] {
            let values: Vec<_> = (0..count)
                .map(|i| ((i * 17 % 41) as f64 - 20.0) / 13.0)
                .collect();
            assert_eq!(
                norm2(&values).to_bits(),
                fs_math::det::sqrt(dot(&values, &values)).to_bits()
            );
        }
    }
}
