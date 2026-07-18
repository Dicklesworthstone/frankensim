//! G0/G3 conformance baseline for the compiled manifold authority.
//!
//! This fixture deliberately exercises only `fs_opt::Manifold`, the descriptor
//! used by admitted problems and the public retraction boundary. It does not
//! import the separate `fs-ascent` implementation. The tests certify descriptor
//! dimensions, wire/admission preservation, and the currently public
//! retractions. They make no claim yet about a unified projection/transport
//! implementation, canonical quaternion-antipode identity, solver convergence,
//! or `fs-ascent` consumer migration.

#![deny(unsafe_code)]

use fs_opt::{Manifold, OptError, ProblemBuilder, Sense, parse, problem_hash, serialize};
use fs_qty::Dims;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct ExpectedDimensions {
    point: u32,
    parameter: u32,
    tangent: u32,
}

fn assert_dimensions(manifold: Manifold, expected: ExpectedDimensions) {
    assert_eq!(manifold.point_dim(), Some(expected.point));
    assert_eq!(manifold.param_dim(), Some(expected.parameter));
    assert_eq!(manifold.tangent_dim(), Some(expected.tangent));
}

fn rotate_vector(quaternion: &[f64], vector: [f64; 3]) -> [f64; 3] {
    let [w, x, y, z] = quaternion else {
        panic!("SO(3) authority returned non-quaternion point storage");
    };
    let cross = [
        y * vector[2] - z * vector[1],
        z * vector[0] - x * vector[2],
        x * vector[1] - y * vector[0],
    ];
    let cross_twice = [
        y * cross[2] - z * cross[1],
        z * cross[0] - x * cross[2],
        x * cross[1] - y * cross[0],
    ];
    [
        vector[0] + 2.0 * (w * cross[0] + cross_twice[0]),
        vector[1] + 2.0 * (w * cross[1] + cross_twice[1]),
        vector[2] + 2.0 * (w * cross[2] + cross_twice[2]),
    ]
}

fn assert_close(left: f64, right: f64, tolerance: f64) {
    assert!(
        (left - right).abs() <= tolerance,
        "{left:?} differs from {right:?} by more than {tolerance:?}"
    );
}

/// G0: raw descriptors, sealed problems, admission, and the canonical wire
/// round-trip all retain one exact point/parameter/tangent dimension table.
#[test]
fn g0_descriptor_admission_and_wire_paths_share_the_layout_table() {
    let cases = [
        (
            "euclidean",
            Manifold::Rn { dim: 3 },
            ExpectedDimensions {
                point: 3,
                parameter: 3,
                tangent: 3,
            },
        ),
        (
            "sphere",
            Manifold::Sphere { ambient: 4 },
            ExpectedDimensions {
                point: 4,
                parameter: 4,
                tangent: 3,
            },
        ),
        (
            "rotation",
            Manifold::So3,
            ExpectedDimensions {
                point: 4,
                parameter: 3,
                tangent: 3,
            },
        ),
        (
            "frame",
            Manifold::Stiefel { n: 4, p: 2 },
            ExpectedDimensions {
                point: 8,
                parameter: 8,
                tangent: 5,
            },
        ),
        (
            "zero-dimensional-frame",
            Manifold::Stiefel { n: 1, p: 1 },
            ExpectedDimensions {
                point: 1,
                parameter: 1,
                tangent: 0,
            },
        ),
    ];

    let mut builder = ProblemBuilder::new();
    for (name, manifold, expected) in cases {
        assert_dimensions(manifold, expected);
        builder
            .var(name, manifold, Dims::NONE)
            .expect("valid authority-table descriptor");
    }
    let zero = builder
        .konst(0.0, Dims::NONE)
        .expect("finite dimensionless objective");
    builder
        .objective(zero, Sense::Minimize, 1.0)
        .expect("scalar objective");
    let problem = builder.finish();
    let admission = problem.admit().expect("authority fixture admits");

    assert_eq!(problem.vars().len(), cases.len());
    for (variable, (_, expected_manifold, expected)) in problem.vars().iter().zip(cases) {
        assert_eq!(variable.manifold, expected_manifold);
        assert_dimensions(variable.manifold, expected);
    }

    let canonical = serialize(&problem);
    let decoded = parse(&canonical).expect("canonical authority fixture parses");
    assert_eq!(serialize(&decoded), canonical);
    assert_eq!(problem_hash(&decoded), problem_hash(&problem));
    assert_eq!(
        decoded
            .admit()
            .expect("decoded fixture admits")
            .semantic_id(),
        admission.semantic_id()
    );
    for (variable, (_, expected_manifold, expected)) in decoded.vars().iter().zip(cases) {
        assert_eq!(variable.manifold, expected_manifold);
        assert_dimensions(variable.manifold, expected);
    }
}

/// G0/G3: SO(3) points are `(w, x, y, z)` unit quaternions, while public
/// retraction accepts a three-coordinate body-frame Lie increment and composes
/// its exponential on the right.
/// Antipodal input quaternions must still induce the same physical rotation;
/// this does not require or claim a canonical representative yet.
#[test]
fn g0_g3_so3_uses_quaternion_points_and_three_coordinate_increments() {
    let manifold = Manifold::So3;
    assert_dimensions(
        manifold,
        ExpectedDimensions {
            point: 4,
            parameter: 3,
            tangent: 3,
        },
    );

    // A quarter-turn about x followed by a body-frame quarter-turn about y.
    // The noncommuting fixture distinguishes q*exp(step/2) from both a no-op
    // and exp(step/2)*q: the latter would negate the final z component.
    let inverse_sqrt_two = std::f64::consts::FRAC_1_SQRT_2;
    let base = [inverse_sqrt_two, inverse_sqrt_two, 0.0, 0.0];
    let antipode = base.map(|component| -component);
    let omega = [0.0, std::f64::consts::FRAC_PI_2, 0.0];
    let rotated = manifold
        .retract(&base, &omega)
        .expect("three-coordinate SO(3) increment");
    let rotated_antipode = manifold
        .retract(&antipode, &omega)
        .expect("antipodal SO(3) base");

    assert_eq!(rotated.len(), 4);
    assert_close(rotated.iter().map(|value| value * value).sum(), 1.0, 2e-15);
    for (actual, expected) in rotated.iter().zip([0.5, 0.5, 0.5, 0.5]) {
        assert_close(*actual, expected, 2e-15);
    }
    for basis in [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]] {
        let action = rotate_vector(&rotated, basis);
        let antipodal_action = rotate_vector(&rotated_antipode, basis);
        for (left, right) in action.into_iter().zip(antipodal_action) {
            assert_close(left, right, 2e-15);
        }
    }

    for malformed in [&omega[..2], &[omega[0], omega[1], omega[2], 0.0][..]] {
        assert!(matches!(
            manifold.retract(&base, malformed),
            Err(OptError::RetractionLen {
                input: "retraction step",
                expected: 3,
                got,
            }) if got == malformed.len() as u64
        ));
    }
}

/// G0/G3: Stiefel points and QR-retraction parameters use column-major
/// ambient `n*p` storage even though the intrinsic tangent dimension is
/// `n*p - p*(p+1)/2`. Repeated QR landing is bit-stable for one fixed input.
#[test]
fn g0_g3_stiefel_distinguishes_ambient_parameters_from_tangent_dimension() {
    let manifold = Manifold::Stiefel { n: 4, p: 2 };
    assert_dimensions(
        manifold,
        ExpectedDimensions {
            point: 8,
            parameter: 8,
            tangent: 5,
        },
    );

    let base = [
        1.0, 0.0, 0.0, 0.0, // first column
        0.0, 1.0, 0.0, 0.0, // second column
    ];
    let ambient_step = [
        0.0, 0.0, 0.5, 0.0, // tilt first column toward row 2
        0.0, 0.0, 0.0, 0.25, // tilt second column toward row 3
    ];
    let landed = manifold
        .retract(&base, &ambient_step)
        .expect("full ambient QR increment");
    let replay = manifold
        .retract(&base, &ambient_step)
        .expect("deterministic QR replay");
    assert_eq!(
        landed
            .iter()
            .map(|value| value.to_bits())
            .collect::<Vec<_>>(),
        replay
            .iter()
            .map(|value| value.to_bits())
            .collect::<Vec<_>>()
    );
    let expected = [
        2.0 / 5.0_f64.sqrt(),
        0.0,
        1.0 / 5.0_f64.sqrt(),
        0.0,
        0.0,
        4.0 / 17.0_f64.sqrt(),
        0.0,
        1.0 / 17.0_f64.sqrt(),
    ];
    for (actual, expected) in landed.iter().zip(expected) {
        assert_close(*actual, expected, 2e-15);
    }

    for column in 0..2 {
        for against in 0..=column {
            let dot = (0..4)
                .map(|row| landed[column * 4 + row] * landed[against * 4 + row])
                .sum::<f64>();
            let expected = if column == against { 1.0 } else { 0.0 };
            assert_close(dot, expected, 2e-15);
        }
    }

    let intrinsic_sized_step = [0.0; 5];
    assert!(matches!(
        manifold.retract(&base, &intrinsic_sized_step),
        Err(OptError::RetractionLen {
            input: "retraction step",
            expected: 8,
            got: 5,
        })
    ));

    let rank_deficient_step = [
        0.0, 0.0, 0.0, 0.0, // preserve first column
        1.0, -1.0, 0.0, 0.0, // collapse second column onto the first
    ];
    assert!(matches!(
        manifold.retract(&base, &rank_deficient_step),
        Err(OptError::RetractionDomain {
            manifold: "Stiefel",
            ..
        })
    ));
}

/// G0: malformed descriptors are refused as descriptor failures before the
/// public retraction boundary can reinterpret empty storage as a valid point.
#[test]
fn g0_malformed_descriptors_fail_before_storage_dispatch() {
    let malformed = [
        Manifold::Rn { dim: 0 },
        Manifold::Sphere { ambient: 1 },
        Manifold::Stiefel { n: 4, p: 0 },
        Manifold::Stiefel { n: 2, p: 3 },
        Manifold::Stiefel { n: u32::MAX, p: 2 },
    ];

    for manifold in malformed {
        let mut builder = ProblemBuilder::new();
        assert!(matches!(
            builder.var("invalid", manifold, Dims::NONE),
            Err(OptError::ManifoldInvalid { .. })
        ));
        assert!(matches!(
            manifold.retract(&[], &[]),
            Err(OptError::ManifoldInvalid { .. })
        ));
    }
}
