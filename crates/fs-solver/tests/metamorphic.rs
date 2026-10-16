//! Gauntlet G3 relations for the production linear-operator adapter.
//!
//! This supplements the existing nonsymmetric GMRES transpose/LU battery; it
//! deliberately makes only a bounded exact-linear fixture claim.

use fs_propcheck::metamorphic::{
    RelationCase, RelationObservation, Tolerance, adjoint_finite_difference, check_relation,
};
use fs_solver::{CsrOp, LinearOp, dot};

type Vector3 = (f64, f64, f64);
type LinearCase = (Vector3, Vector3);

fn array(vector: Vector3) -> [f64; 3] {
    [vector.0, vector.1, vector.2]
}

fn canonical_zero(value: f64) -> f64 {
    if value == 0.0 { 0.0 } else { value }
}

fn nonsymmetric_fixture() -> CsrOp {
    let mut matrix = fs_sparse::Coo::new(3, 3);
    for (row, column, value) in [
        (0, 0, 2.0),
        (0, 1, -1.0),
        (0, 2, 3.0),
        (1, 1, 4.0),
        (1, 2, 1.0),
        (2, 0, -2.0),
        (2, 1, 5.0),
        (2, 2, 2.0),
    ] {
        matrix.push(row, column, value);
    }
    CsrOp::general(matrix.assemble())
}

#[test]
fn g3_general_csr_action_matches_its_adjoint_directional_difference() {
    let matrix = nonsymmetric_fixture();
    let objective = [1.0, -2.0, 3.0];
    let operator = |&(x, direction): &LinearCase| {
        let mut applied = [0.0; 3];
        matrix.apply(&array(x), &mut applied);
        let mut adjoint = [0.0; 3];
        matrix.apply_transpose(&objective, &mut adjoint);
        (dot(&objective, &applied), dot(&array(direction), &adjoint))
    };
    let relation = adjoint_finite_difference(
        "general-csr-linear-adjoint-difference",
        Tolerance::Exact,
        |&(x, direction): &LinearCase, &step: &i64| {
            let scale = step as f64;
            (
                (
                    x.0 + scale * direction.0,
                    x.1 + scale * direction.1,
                    x.2 + scale * direction.2,
                ),
                direction,
            )
        },
        |&(base_value, base_adjoint): &(f64, f64),
         &(transformed_value, transformed_adjoint): &(f64, f64),
         &step: &i64,
         tolerance| {
            let difference = canonical_zero(transformed_value - base_value);
            let expected = canonical_zero(base_adjoint * step as f64);
            let directional = tolerance.evaluate_scalar(expected, difference);
            let stable_adjoint = tolerance.evaluate_scalar(
                canonical_zero(base_adjoint),
                canonical_zero(transformed_adjoint),
            );
            RelationObservation::new(
                directional.margin().min(stable_adjoint.margin()),
                "finite difference equals the transposed-action directional derivative",
            )
        },
    );

    check_relation(
        "fs-solver::CsrOp::general",
        0x2ACE_0004,
        384,
        |stream| {
            let mut direction = (
                stream.int_in(-8, 8) as f64,
                stream.int_in(-8, 8) as f64,
                stream.int_in(-8, 8) as f64,
            );
            if direction == (0.0, 0.0, 0.0) {
                direction.0 = 1.0;
            }
            RelationCase::new(
                (
                    (
                        stream.int_in(-8, 8) as f64,
                        stream.int_in(-8, 8) as f64,
                        stream.int_in(-8, 8) as f64,
                    ),
                    direction,
                ),
                stream.int_in(-3, 3),
            )
        },
        &operator,
        &relation,
    );
}
