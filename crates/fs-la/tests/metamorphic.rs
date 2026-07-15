//! Gauntlet G3 relations for the production dense-operator surface.
//!
//! These checks supplement the crate's existing G0 adjoint identity and
//! numerical-oracle batteries; they do not replace those pinned cases.

use fs_propcheck::metamorphic::{
    RelationCase, RelationObservation, Tolerance, check_relation, unit_rescaling,
};

type Matrix2 = ((f64, f64), (f64, f64));
type Vector2 = (f64, f64);
type GemvCase = (Matrix2, Vector2);

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
        |&(matrix, vector): &GemvCase, &scale: &f64| (matrix, (vector.0 * scale, vector.1 * scale)),
        |base: &[f64; 2], transformed: &[f64; 2], &scale: &f64, tolerance| {
            finite_component_margin(tolerance, [base[0] * scale, base[1] * scale], *transformed)
        },
    );

    check_relation(
        "fs-la::gemm_f64",
        0x2ACE_0001,
        384,
        |stream| {
            let scalar = |stream: &mut fs_propcheck::Stream| stream.f64_in(-16.0, 16.0);
            RelationCase::new(
                (
                    (
                        (scalar(stream), scalar(stream)),
                        (scalar(stream), scalar(stream)),
                    ),
                    (scalar(stream), scalar(stream)),
                ),
                stream.f64_in(-4.0, 4.0),
            )
        },
        &operator,
        &relation,
    );
}
