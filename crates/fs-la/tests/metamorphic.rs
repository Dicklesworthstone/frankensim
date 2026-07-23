//! Gauntlet G3 relations for the production dense-operator surface.
//!
//! These checks supplement the crate's existing G0 adjoint identity and
//! numerical-oracle batteries; they do not replace those pinned cases.

use fs_propcheck::Shrink;
use fs_propcheck::metamorphic::{
    RelationCase, RelationObservation, Tolerance, check_relation, unit_rescaling,
};

type Matrix2 = ((f64, f64), (f64, f64));
type Vector2 = (f64, f64);
type GemvCase = (Matrix2, Vector2);

/// A nonidentity coherent vector-unit rescaling, represented exactly as 2^k.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct UnitScaleExponent(i64);

impl UnitScaleExponent {
    fn factor(self) -> f64 {
        fs_math::det::powi(
            2.0,
            i32::try_from(self.0).expect("generated unit exponent fits i32"),
        )
    }
}

impl Shrink for UnitScaleExponent {
    fn shrink_candidates(&self) -> Vec<Self> {
        let exponent = self.0;
        debug_assert_ne!(exponent, 0, "unit rescaling must be nonidentity");
        if exponent.unsigned_abs() == 1 {
            return Vec::new();
        }

        let sign = exponent.signum();
        let mut candidates = vec![sign, exponent / 2, exponent - sign];
        candidates.retain(|candidate| {
            *candidate != 0 && candidate.unsigned_abs() < exponent.unsigned_abs()
        });
        candidates.dedup();
        candidates.into_iter().map(Self).collect()
    }
}

fn finite_component_margin(
    tolerance: Tolerance,
    expected: [f64; 2],
    observed: [f64; 2],
) -> RelationObservation {
    let first = tolerance.evaluate_scalar(expected[0], observed[0]);
    let second = tolerance.evaluate_scalar(expected[1], observed[1]);
    RelationObservation::new(
        first.margin().min(second.margin()),
        "both GEMV components obey the declared scale equivariance",
    )
}

#[test]
fn g3_gemv_is_equivariant_under_coherent_vector_rescaling() {
    let operator = |&(matrix, vector): &GemvCase| {
        // `matrix` is dimensionless. Both vector components share one unit, so
        // one positive factor coherently rescales the input and output units.
        let a = [matrix.0.0, matrix.0.1, matrix.1.0, matrix.1.1];
        let x = [vector.0, vector.1];
        let mut y = [0.0; 2];
        fs_la::gemm_f64(2, 1, 2, 1.0, &a, &x, 0.0, &mut y);
        y
    };
    let relation = unit_rescaling(
        "gemv-vector-scale-equivariance",
        Tolerance::AbsoluteRelative {
            max_abs: 2.0e-12,
            max_relative: 2.0e-12,
        },
        |&(matrix, vector): &GemvCase, exponent: &UnitScaleExponent| {
            let scale = exponent.factor();
            (matrix, (vector.0 * scale, vector.1 * scale))
        },
        |base: &[f64; 2], transformed: &[f64; 2], exponent: &UnitScaleExponent, tolerance| {
            let scale = exponent.factor();
            finite_component_margin(tolerance, [base[0] * scale, base[1] * scale], *transformed)
        },
    );

    check_relation(
        "fs-la::gemm_f64",
        0x2ACE_0001,
        384,
        |stream| {
            let scalar = |stream: &mut fs_propcheck::Stream| stream.f64_in(-16.0, 16.0);
            let raw_exponent = stream.int_in(-4, 3);
            let exponent = if raw_exponent >= 0 {
                raw_exponent + 1
            } else {
                raw_exponent
            };
            RelationCase::new(
                (
                    (
                        (scalar(stream), scalar(stream)),
                        (scalar(stream), scalar(stream)),
                    ),
                    (scalar(stream), scalar(stream)),
                ),
                UnitScaleExponent(exponent),
            )
        },
        &operator,
        &relation,
    );
}
