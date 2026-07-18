//! G0/G3 conformance for the authoritative runtime manifold operations.
//!
//! Fixed finite fixtures independently check ambient-gradient pullback,
//! retraction-curve velocity, and tangent transport for Rn, Sphere, SO(3), and
//! Stiefel. The SO(3) fixture distinguishes four-coordinate quaternion points
//! from three-coordinate right/body parameters. The Stiefel fixture compares
//! the differentiated production QR program with an external central
//! difference and treats that differential as a non-isometric transport.
//!
//! This target makes no SO(3) sign-seam/global-chart, Stiefel transport
//! isometry, arbitrary-conditioning, solver-convergence, persistence,
//! cross-ISA, bounded-work/cancellation-latency, or performance claim.

#![deny(unsafe_code)]

use fs_opt::{Manifold, OptError, RetractionCurve};

const SO3_BASE: [f64; 4] = [0.5, 0.5, 0.5, 0.5];
const SO3_DIRECTION: [f64; 3] = [0.7, -0.4, 0.9];
const SO3_AMBIENT_GRADIENT: [f64; 4] = [1.25, -0.5, 2.0, -1.5];

const SPHERE: Manifold = Manifold::Sphere { ambient: 4 };
const SPHERE_BASE: [f64; 4] = [0.5; 4];
const SPHERE_DIRECTION: [f64; 4] = [0.75, -0.25, 0.5, -0.5];
const SPHERE_AMBIENT_GRADIENT: [f64; 4] = [-0.5, 1.25, 0.75, -1.0];

const N: usize = 4;
const P: usize = 2;
const STORAGE: usize = N * P;
const STIEFEL: Manifold = Manifold::Stiefel { n: 4, p: 2 };
const STIEFEL_BASE: [f64; STORAGE] = [
    0.5, 0.5, 0.5, 0.5, // first column
    0.5, -0.5, 0.5, -0.5, // second column
];
const STIEFEL_STEP: [f64; STORAGE] = [
    0.0, 0.0, 0.25, -0.25, // first candidate column increment
    0.125, 0.25, -0.125, 0.0, // second candidate column increment
];
const STIEFEL_TANGENT: [f64; STORAGE] = [
    0.4375, 0.0625, -0.0625, -0.4375, // 3/8 x1 + 1/2 z0
    -0.0625, -0.3125, -0.3125, -0.0625, // -3/8 x0 + 1/4 z1
];
const STIEFEL_AMBIENT_GRADIENT: [f64; STORAGE] = [0.75, -0.5, 0.25, 1.0, -0.25, 0.5, 1.25, -0.75];

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|value| value.to_bits()).collect()
}

fn dot(left: &[f64], right: &[f64]) -> f64 {
    assert_eq!(left.len(), right.len());
    left.iter().zip(right).map(|(a, b)| a * b).sum()
}

fn norm(values: &[f64]) -> f64 {
    dot(values, values).sqrt()
}

fn max_error(left: &[f64], right: &[f64]) -> f64 {
    assert_eq!(left.len(), right.len());
    left.iter()
        .zip(right)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max)
}

fn scaled(values: &[f64], scale: f64) -> Vec<f64> {
    values.iter().map(|value| scale * value).collect()
}

fn affine_step(step: &[f64], variation: &[f64], scale: f64) -> Vec<f64> {
    step.iter()
        .zip(variation)
        .map(|(base, delta)| scale.mul_add(*delta, *base))
        .collect()
}

fn central_curve_velocity(
    manifold: Manifold,
    point: &[f64],
    direction: &[f64],
    alpha: f64,
    h: f64,
) -> Vec<f64> {
    let plus = manifold
        .retract(point, &scaled(direction, alpha + h))
        .expect("positive curve sample");
    let minus = manifold
        .retract(point, &scaled(direction, alpha - h))
        .expect("negative curve sample");
    plus.iter()
        .zip(&minus)
        .map(|(positive, negative)| (positive - negative) / (2.0 * h))
        .collect()
}

fn central_retraction_differential(
    manifold: Manifold,
    point: &[f64],
    step: &[f64],
    variation: &[f64],
    h: f64,
) -> Vec<f64> {
    let plus = manifold
        .retract(point, &affine_step(step, variation, h))
        .expect("positive differential sample");
    let minus = manifold
        .retract(point, &affine_step(step, variation, -h))
        .expect("negative differential sample");
    plus.iter()
        .zip(&minus)
        .map(|(positive, negative)| (positive - negative) / (2.0 * h))
        .collect()
}

fn column(values: &[f64], index: usize) -> &[f64] {
    &values[index * N..(index + 1) * N]
}

fn stiefel_tangent_residual(point: &[f64], tangent: &[f64]) -> f64 {
    let mut residual = 0.0_f64;
    for row in 0..P {
        for column_index in 0..P {
            residual = residual.max(
                (dot(column(point, row), column(tangent, column_index))
                    + dot(column(tangent, row), column(point, column_index)))
                .abs(),
            );
        }
    }
    residual
}

fn independent_stiefel_projection(point: &[f64], ambient: &[f64]) -> Vec<f64> {
    let mut gram = [[0.0; P]; P];
    for (left, row) in gram.iter_mut().enumerate() {
        for (right, value) in row.iter_mut().enumerate() {
            *value = dot(column(point, left), column(ambient, right));
        }
    }
    let mut projected = vec![0.0; STORAGE];
    for output_column in 0..P {
        for row in 0..N {
            let correction: f64 = (0..P)
                .map(|basis| {
                    point[basis * N + row]
                        * 0.5
                        * (gram[basis][output_column] + gram[output_column][basis])
                })
                .sum();
            projected[output_column * N + row] = ambient[output_column * N + row] - correction;
        }
    }
    projected
}

fn independent_so3_parameter_gradient(point: &[f64], ambient: &[f64]) -> [f64; 3] {
    let [w, x, y, z] = [point[0], point[1], point[2], point[3]];
    let [gw, gx, gy, gz] = [ambient[0], ambient[1], ambient[2], ambient[3]];
    [
        0.5 * (-x * gw + w * gx + z * gy - y * gz),
        0.5 * (-y * gw - z * gx + w * gy + x * gz),
        0.5 * (-z * gw + y * gx - x * gy + w * gz),
    ]
}

fn independent_so3_curve_velocity(point: &[f64], direction: &[f64]) -> [f64; 4] {
    let [w, x, y, z] = [point[0], point[1], point[2], point[3]];
    let [a, b, c] = [direction[0], direction[1], direction[2]];
    [
        -0.5 * (x * a + y * b + z * c),
        0.5 * (w * a + y * c - z * b),
        0.5 * (w * b - x * c + z * a),
        0.5 * (w * c + x * b - y * a),
    ]
}

fn assert_curve_landing(
    manifold: Manifold,
    point: &[f64],
    direction: &[f64],
    alpha: f64,
    curve: &RetractionCurve,
) {
    let direct = manifold
        .retract(point, &scaled(direction, alpha))
        .expect("direct public retraction");
    assert_eq!(bits(&curve.point), bits(&direct));
}

fn assert_point_parameter_duality(
    manifold: Manifold,
    point: &[f64],
    direction: &[f64],
    ambient_gradient: &[f64],
) {
    let curve = manifold
        .retract_curve(point, direction, 0.0)
        .expect("zero-parameter retraction curve");
    let parameter_gradient = manifold
        .parameter_gradient(point, ambient_gradient)
        .expect("ambient gradient pullback");
    let parameter_pairing = dot(&parameter_gradient, &curve.velocity);
    let point_pairing = if manifold == Manifold::So3 {
        dot(
            ambient_gradient,
            &independent_so3_curve_velocity(point, direction),
        )
    } else {
        dot(ambient_gradient, &curve.velocity)
    };
    let tolerance = 32.0 * f64::EPSILON * (1.0 + parameter_pairing.abs());
    assert!((parameter_pairing - point_pairing).abs() <= tolerance);
}

#[test]
fn g0_parameter_gradients_respect_point_and_parameter_representations() {
    let rn = Manifold::Rn { dim: 2 };
    assert_eq!(
        bits(
            &rn.parameter_gradient(&[1.0, -2.0], &[3.0, -4.0])
                .expect("Rn parameter gradient")
        ),
        bits(&[3.0, -4.0])
    );

    let sphere = SPHERE
        .parameter_gradient(&SPHERE_BASE, &SPHERE_AMBIENT_GRADIENT)
        .expect("Sphere parameter gradient");
    let sphere_radial = dot(&SPHERE_BASE, &SPHERE_AMBIENT_GRADIENT);
    let expected_sphere: [f64; 4] = core::array::from_fn(|index| {
        SPHERE_AMBIENT_GRADIENT[index] - sphere_radial * SPHERE_BASE[index]
    });
    assert!(max_error(&sphere, &expected_sphere) <= 2.0 * f64::EPSILON);
    assert!(dot(&SPHERE_BASE, &sphere).abs() <= 2.0 * f64::EPSILON);

    let so3 = Manifold::So3
        .parameter_gradient(&SO3_BASE, &SO3_AMBIENT_GRADIENT)
        .expect("SO(3) parameter gradient");
    let expected_so3 = independent_so3_parameter_gradient(&SO3_BASE, &SO3_AMBIENT_GRADIENT);
    assert_eq!(so3.len(), 3);
    assert_eq!(bits(&so3), bits(&expected_so3));
    let antipode = SO3_BASE.map(|value| -value);
    let antipodal_gradient = SO3_AMBIENT_GRADIENT.map(|value| -value);
    assert_eq!(
        bits(
            &Manifold::So3
                .parameter_gradient(&antipode, &antipodal_gradient)
                .expect("antipodally equivariant SO(3) pullback")
        ),
        bits(&so3)
    );

    let stiefel = STIEFEL
        .parameter_gradient(&STIEFEL_BASE, &STIEFEL_AMBIENT_GRADIENT)
        .expect("Stiefel parameter gradient");
    let expected_stiefel = independent_stiefel_projection(&STIEFEL_BASE, &STIEFEL_AMBIENT_GRADIENT);
    assert_eq!(stiefel.len(), STORAGE);
    assert!(max_error(&stiefel, &expected_stiefel) <= 4.0 * f64::EPSILON);
    assert!(stiefel_tangent_residual(&STIEFEL_BASE, &stiefel) <= 4.0e-15);
}

#[test]
fn g0_parameter_gradients_are_dual_to_zero_parameter_curve_velocities() {
    assert_point_parameter_duality(
        Manifold::Rn { dim: 2 },
        &[1.0, -2.0],
        &[0.25, 0.75],
        &[3.0, -4.0],
    );
    assert_point_parameter_duality(
        SPHERE,
        &SPHERE_BASE,
        &SPHERE_DIRECTION,
        &SPHERE_AMBIENT_GRADIENT,
    );
    assert_point_parameter_duality(
        Manifold::So3,
        &SO3_BASE,
        &SO3_DIRECTION,
        &SO3_AMBIENT_GRADIENT,
    );
    assert_point_parameter_duality(
        STIEFEL,
        &STIEFEL_BASE,
        &STIEFEL_TANGENT,
        &STIEFEL_AMBIENT_GRADIENT,
    );
}

#[test]
#[allow(clippy::too_many_lines)]
fn g3_curve_outputs_match_public_retractions_and_independent_velocities() {
    let rn = Manifold::Rn { dim: 2 };
    let rn_point = [1.0, -2.0];
    let rn_direction = [0.25, 0.75];
    let rn_curve = rn
        .retract_curve(&rn_point, &rn_direction, 0.5)
        .expect("Rn curve");
    assert_curve_landing(rn, &rn_point, &rn_direction, 0.5, &rn_curve);
    assert_eq!(bits(&rn_curve.velocity), bits(&rn_direction));

    let sphere_alpha = 0.25;
    let sphere_curve = SPHERE
        .retract_curve(&SPHERE_BASE, &SPHERE_DIRECTION, sphere_alpha)
        .expect("Sphere curve");
    assert_curve_landing(
        SPHERE,
        &SPHERE_BASE,
        &SPHERE_DIRECTION,
        sphere_alpha,
        &sphere_curve,
    );
    let moved: Vec<f64> = SPHERE_BASE
        .iter()
        .zip(SPHERE_DIRECTION)
        .map(|(point, direction)| *point + sphere_alpha * direction)
        .collect();
    let moved_norm = norm(&moved);
    let radial = dot(&sphere_curve.point, &SPHERE_DIRECTION);
    let expected_sphere: Vec<f64> = SPHERE_DIRECTION
        .iter()
        .zip(&sphere_curve.point)
        .map(|(direction, point)| (*direction - radial * *point) / moved_norm)
        .collect();
    assert!(max_error(&sphere_curve.velocity, &expected_sphere) <= 4.0 * f64::EPSILON);
    assert!(dot(&sphere_curve.point, &sphere_curve.velocity).abs() <= 8.0 * f64::EPSILON);
    let sphere_coarse = central_curve_velocity(
        SPHERE,
        &SPHERE_BASE,
        &SPHERE_DIRECTION,
        sphere_alpha,
        1.0 / 1024.0,
    );
    let sphere_fine = central_curve_velocity(
        SPHERE,
        &SPHERE_BASE,
        &SPHERE_DIRECTION,
        sphere_alpha,
        1.0 / 2048.0,
    );
    let sphere_coarse_error = max_error(&sphere_coarse, &sphere_curve.velocity);
    let sphere_fine_error = max_error(&sphere_fine, &sphere_curve.velocity);
    assert!(sphere_fine_error < sphere_coarse_error * 0.27);
    assert!(sphere_fine_error <= 2.0e-7);

    let so3_alpha = 0.375;
    let so3_curve = Manifold::So3
        .retract_curve(&SO3_BASE, &SO3_DIRECTION, so3_alpha)
        .expect("SO(3) curve");
    assert_curve_landing(
        Manifold::So3,
        &SO3_BASE,
        &SO3_DIRECTION,
        so3_alpha,
        &so3_curve,
    );
    assert_eq!(so3_curve.velocity.len(), 3);
    assert_eq!(bits(&so3_curve.velocity), bits(&SO3_DIRECTION));
    let antipodal_curve = Manifold::So3
        .retract_curve(
            &SO3_BASE.map(|coordinate| -coordinate),
            &SO3_DIRECTION,
            so3_alpha,
        )
        .expect("antipodal SO(3) curve");
    assert_eq!(bits(&antipodal_curve.point), bits(&so3_curve.point));
    assert_eq!(bits(&antipodal_curve.velocity), bits(&so3_curve.velocity));
    let expected_so3_point_velocity =
        independent_so3_curve_velocity(&so3_curve.point, &SO3_DIRECTION);
    let so3_coarse = central_curve_velocity(
        Manifold::So3,
        &SO3_BASE,
        &SO3_DIRECTION,
        so3_alpha,
        1.0 / 1024.0,
    );
    let so3_fine = central_curve_velocity(
        Manifold::So3,
        &SO3_BASE,
        &SO3_DIRECTION,
        so3_alpha,
        1.0 / 2048.0,
    );
    let so3_coarse_error = max_error(&so3_coarse, &expected_so3_point_velocity);
    let so3_fine_error = max_error(&so3_fine, &expected_so3_point_velocity);
    assert!(so3_fine_error < so3_coarse_error * 0.27);
    assert!(so3_fine_error <= 2.0e-7);
    assert!(dot(&so3_curve.point, &expected_so3_point_velocity).abs() <= 2.0e-16);

    let stiefel_alpha = 0.3;
    let stiefel_curve = STIEFEL
        .retract_curve(&STIEFEL_BASE, &STIEFEL_TANGENT, stiefel_alpha)
        .expect("Stiefel curve");
    assert_curve_landing(
        STIEFEL,
        &STIEFEL_BASE,
        &STIEFEL_TANGENT,
        stiefel_alpha,
        &stiefel_curve,
    );
    let coarse = central_curve_velocity(
        STIEFEL,
        &STIEFEL_BASE,
        &STIEFEL_TANGENT,
        stiefel_alpha,
        1.0 / 512.0,
    );
    let fine = central_curve_velocity(
        STIEFEL,
        &STIEFEL_BASE,
        &STIEFEL_TANGENT,
        stiefel_alpha,
        1.0 / 1024.0,
    );
    let coarse_error = max_error(&coarse, &stiefel_curve.velocity);
    let fine_error = max_error(&fine, &stiefel_curve.velocity);
    assert!(fine_error < coarse_error * 0.27);
    assert!(fine_error <= 2.0e-7);
    assert!(stiefel_tangent_residual(&stiefel_curve.point, &stiefel_curve.velocity) <= 2.0e-14);
}

#[test]
#[allow(clippy::too_many_lines)]
fn g3_transport_is_tangent_replayable_and_matches_declared_geometry() {
    let rn = Manifold::Rn { dim: 2 };
    let rn_from = [1.0, -2.0];
    let rn_step = [0.5, 0.25];
    let rn_to = rn.retract(&rn_from, &rn_step).expect("Rn destination");
    let rn_vector = [3.0, -4.0];
    assert_eq!(
        bits(
            &rn.transport_parameter(&rn_from, &rn_step, &rn_to, &rn_vector)
                .expect("Rn transport")
        ),
        bits(&rn_vector)
    );

    let sphere = Manifold::Sphere { ambient: 3 };
    let sphere_from = [1.0, 0.0, 0.0];
    let sphere_step = [0.0, 0.5, 0.0];
    let sphere_to = sphere
        .retract(&sphere_from, &sphere_step)
        .expect("Sphere destination");
    let first = [0.0, 1.0, 0.0];
    let second = [0.0, 0.0, 1.0];
    let moved_first = sphere
        .transport_parameter(&sphere_from, &sphere_step, &sphere_to, &first)
        .expect("first Sphere transport");
    let moved_second = sphere
        .transport_parameter(&sphere_from, &sphere_step, &sphere_to, &second)
        .expect("second Sphere transport");
    assert!(dot(&sphere_to, &moved_first).abs() <= 8.0 * f64::EPSILON);
    assert!(dot(&sphere_to, &moved_second).abs() <= 8.0 * f64::EPSILON);
    assert!((dot(&first, &second) - dot(&moved_first, &moved_second)).abs() <= 8.0 * f64::EPSILON);
    assert!((norm(&first) - norm(&moved_first)).abs() <= 8.0 * f64::EPSILON);

    let so3_step = [0.1, -0.2, 0.3];
    let so3_to = Manifold::So3
        .retract(&SO3_BASE, &so3_step)
        .expect("SO(3) destination");
    let so3_vector = [1.0, 2.0, -0.5];
    let so3_transport = Manifold::So3
        .transport_parameter(&SO3_BASE, &so3_step, &so3_to, &so3_vector)
        .expect("SO(3) body-coordinate transport");
    assert_eq!(bits(&so3_transport), bits(&so3_vector));

    let stiefel_step_direction = independent_stiefel_projection(&STIEFEL_BASE, &STIEFEL_STEP);
    assert!(stiefel_tangent_residual(&STIEFEL_BASE, &stiefel_step_direction) <= 2.0e-15);
    let stiefel_step = scaled(&stiefel_step_direction, 0.25);
    let stiefel_to = STIEFEL
        .retract(&STIEFEL_BASE, &stiefel_step)
        .expect("Stiefel destination");
    let transported = STIEFEL
        .transport_parameter(&STIEFEL_BASE, &stiefel_step, &stiefel_to, &STIEFEL_TANGENT)
        .expect("differentiated QR transport");
    let coarse = central_retraction_differential(
        STIEFEL,
        &STIEFEL_BASE,
        &stiefel_step,
        &STIEFEL_TANGENT,
        1.0 / 512.0,
    );
    let fine = central_retraction_differential(
        STIEFEL,
        &STIEFEL_BASE,
        &stiefel_step,
        &STIEFEL_TANGENT,
        1.0 / 1024.0,
    );
    let coarse_error = max_error(&coarse, &transported);
    let fine_error = max_error(&fine, &transported);
    assert!(fine_error < coarse_error * 0.27);
    assert!(fine_error <= 2.0e-7);
    assert!(stiefel_tangent_residual(&stiefel_to, &transported) <= 2.0e-14);
    assert_eq!(
        bits(
            &STIEFEL
                .transport_parameter(&STIEFEL_BASE, &stiefel_step, &stiefel_to, &STIEFEL_TANGENT,)
                .expect("Stiefel transport replay")
        ),
        bits(&transported)
    );
}

#[test]
fn g0_malformed_operation_shapes_and_nonfinite_inputs_fail_through_typed_errors() {
    assert!(matches!(
        Manifold::So3.parameter_gradient(&SO3_BASE[..3], &SO3_AMBIENT_GRADIENT),
        Err(OptError::RetractionLen {
            input: "manifold operation point",
            expected: 4,
            got: 3,
        })
    ));
    assert!(matches!(
        Manifold::So3.parameter_gradient(&SO3_BASE, &SO3_AMBIENT_GRADIENT[..3]),
        Err(OptError::RetractionLen {
            input: "ambient manifold gradient",
            expected: 4,
            got: 3,
        })
    ));
    let mut nonfinite_gradient = SO3_AMBIENT_GRADIENT;
    nonfinite_gradient[2] = f64::NAN;
    assert!(matches!(
        Manifold::So3.parameter_gradient(&SO3_BASE, &nonfinite_gradient),
        Err(OptError::RetractionNonFinite {
            input: "ambient manifold gradient",
            component: 2,
            bits,
        }) if bits == f64::NAN.to_bits()
    ));
    let mut nonfinite_direction = SO3_DIRECTION;
    nonfinite_direction[1] = f64::INFINITY;
    assert!(matches!(
        Manifold::So3.retract_curve(&SO3_BASE, &nonfinite_direction, 1.0),
        Err(OptError::RetractionNonFinite {
            input: "retraction curve direction",
            component: 1,
            bits,
        }) if bits == f64::INFINITY.to_bits()
    ));
    assert!(matches!(
        Manifold::So3.retract_curve(&SO3_BASE, &[0.0; 4], 1.0),
        Err(OptError::RetractionLen {
            input: "retraction curve direction",
            expected: 3,
            got: 4,
        })
    ));
    assert!(matches!(
        Manifold::So3.retract_curve(&SO3_BASE, &SO3_DIRECTION, f64::NAN),
        Err(OptError::RetractionNonFinite {
            input: "retraction curve alpha",
            component: 0,
            bits,
        }) if bits == f64::NAN.to_bits()
    ));
}

#[test]
fn g0_malformed_operation_domains_fail_through_typed_errors() {
    assert!(matches!(
        SPHERE.validate_parameter_tangent(&SPHERE_BASE, &SPHERE_BASE),
        Err(OptError::RetractionDomain {
            manifold: "Sphere",
            what: "parameter vector must belong to the point tangent space",
            location: None,
            ..
        })
    ));
    assert!(matches!(
        SPHERE.retract_curve(&[1.0, 1.0, 0.0, 0.0], &SPHERE_DIRECTION, 0.0),
        Err(OptError::RetractionDomain {
            manifold: "Sphere",
            what: "point must have unit norm before retraction",
            location: None,
            measurement_bits,
        }) if measurement_bits == 2.0_f64.to_bits()
    ));

    let sphere = Manifold::Sphere { ambient: 3 };
    let from = [1.0, 0.0, 0.0];
    let step = [0.0, 0.25, 0.0];
    let incompatible_to = [0.0, 1.0, 0.0];
    assert!(matches!(
        sphere.transport_parameter(&from, &step, &incompatible_to, &[0.0, 1.0, 0.0]),
        Err(OptError::RetractionDomain {
            manifold: "Sphere",
            what: "transport destination must equal the authoritative retraction",
            location: None,
            ..
        })
    ));

    let mut collapse_first_column = [0.0; STORAGE];
    for (increment, coordinate) in collapse_first_column[..N]
        .iter_mut()
        .zip(&STIEFEL_BASE[..N])
    {
        *increment = -*coordinate;
    }
    assert!(matches!(
        STIEFEL.retract_curve(&STIEFEL_BASE, &collapse_first_column, 1.0),
        Err(OptError::RetractionDomain {
            manifold: "Stiefel",
            what: "candidate column is rank-deficient",
            location: Some((0, 0)),
            ..
        })
    ));
}
