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

fn accepted_near_unit_scale() -> f64 {
    // 32,768 ulps above one gives a nonzero norm defect comfortably inside the
    // public 1e-10 manifold-membership envelope.
    f64::from_bits(1.0_f64.to_bits() + 32_768)
}

fn xorshift64(state: &mut u64) -> u64 {
    *state ^= *state << 13;
    *state ^= *state >> 7;
    *state ^= *state << 17;
    *state
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

fn independent_stiefel_sylvester_projection(
    point: &[f64],
    ambient: &[f64],
) -> (Vec<f64>, [[f64; P]; P]) {
    // Independent closed-form solve for the symmetric 2x2 correction S in
    //     M S + S M = X^T G + G^T X,  M = X^T X.
    // This deliberately does not reuse production's Richardson iteration.
    let m00 = dot(column(point, 0), column(point, 0));
    let m01 = dot(column(point, 0), column(point, 1));
    let m11 = dot(column(point, 1), column(point, 1));
    let a00 = 2.0 * dot(column(point, 0), column(ambient, 0));
    let a01 = dot(column(point, 0), column(ambient, 1)) + dot(column(ambient, 0), column(point, 1));
    let a11 = 2.0 * dot(column(point, 1), column(ambient, 1));
    let denominator = m00 + m11 - m01 * m01 * (m00.recip() + m11.recip());
    let s01 = (a01 - m01 * (a00 / (2.0 * m00) + a11 / (2.0 * m11))) / denominator;
    let s00 = (0.5 * a00 - m01 * s01) / m00;
    let s11 = (0.5 * a11 - m01 * s01) / m11;
    let correction = [[s00, s01], [s01, s11]];

    let mut projected = vec![0.0; STORAGE];
    for output_column in 0..P {
        for row in 0..N {
            let normal = point[row] * correction[0][output_column]
                + point[N + row] * correction[1][output_column];
            projected[output_column * N + row] = ambient[output_column * N + row] - normal;
        }
    }
    (projected, correction)
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
    assert!(
        max_error(&sphere, &expected_sphere)
            <= 16.0 * f64::EPSILON * (1.0 + norm(&expected_sphere))
    );
    assert!(dot(&SPHERE_BASE, &sphere).abs() <= 16.0 * f64::EPSILON);

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
    assert!(
        max_error(&stiefel, &expected_stiefel)
            <= 32.0 * f64::EPSILON * (1.0 + norm(&expected_stiefel))
    );
    assert!(stiefel_tangent_residual(&STIEFEL_BASE, &stiefel) <= 8.0e-15);
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
#[allow(clippy::too_many_lines)] // one audit keeps all near-manifold projection obligations together
fn g0_projection_is_total_and_idempotent_on_accepted_near_manifold_points() {
    let scale = accepted_near_unit_scale();
    let sphere = Manifold::Sphere { ambient: 3 };
    let sphere_point = [scale, 0.0, 0.0];
    let sphere_norm_defect = dot(&sphere_point, &sphere_point) - 1.0;
    assert!(
        sphere_norm_defect != 0.0 && sphere_norm_defect.abs() <= 1.0e-10,
        "fixture must be accepted but not bit-perfect: scale={scale:.17e}; norm_defect={sphere_norm_defect:.17e}"
    );
    let sphere_radial = sphere
        .parameter_gradient(&sphere_point, &sphere_point)
        .expect("near-unit Sphere radial gradient must project without refusal");
    let sphere_radial_scale = sphere_radial
        .iter()
        .map(|value| value.abs())
        .fold(0.0_f64, f64::max);
    assert!(
        sphere_radial_scale <= 4.0 * f64::EPSILON,
        "purely radial Sphere gradient must collapse within roundoff without refusal: point={sphere_point:?}; projected={sphere_radial:?}; max={sphere_radial_scale:.17e}"
    );

    let sphere_landed = sphere
        .retract(&[1.0, 0.0, 0.0], &[0.125, 0.75, -0.25])
        .expect("non-axis Sphere retraction");
    let sphere_first = sphere
        .parameter_gradient(&sphere_landed, &[0.75, -1.25, 0.5])
        .expect("Sphere projection at retraction-produced point");
    let sphere_second = sphere
        .parameter_gradient(&sphere_landed, &sphere_first)
        .expect("Sphere projection idempotence");
    let sphere_idempotence_error = max_error(&sphere_first, &sphere_second);
    assert!(
        sphere_idempotence_error <= 2.0e-14 * (1.0 + norm(&sphere_first)),
        "Sphere projection must be idempotent at a retraction-produced point: point={sphere_landed:?}; first={sphere_first:?}; second={sphere_second:?}; error={sphere_idempotence_error:.17e}"
    );
    sphere
        .validate_parameter_tangent(&sphere_landed, &sphere_first)
        .expect("Sphere projection postcondition");

    let stiefel_line = Manifold::Stiefel { n: 3, p: 1 };
    let stiefel_near = [scale, 0.0, 0.0];
    let stiefel_radial = stiefel_line
        .parameter_gradient(&stiefel_near, &stiefel_near)
        .expect("near-orthonormal Stiefel normal gradient must project without refusal");
    let stiefel_radial_scale = stiefel_radial
        .iter()
        .map(|value| value.abs())
        .fold(0.0_f64, f64::max);
    assert!(
        stiefel_radial_scale <= 4.0 * f64::EPSILON,
        "purely normal Stiefel(n,1) gradient must collapse within roundoff without refusal: point={stiefel_near:?}; projected={stiefel_radial:?}; max={stiefel_radial_scale:.17e}"
    );

    let cross_gram = 1.0 / (1_u64 << 35) as f64;
    let mut stiefel_crossed = STIEFEL_BASE;
    for row in 0..N {
        stiefel_crossed[N + row] += cross_gram * STIEFEL_BASE[row];
    }
    let observed_cross_gram = dot(column(&stiefel_crossed, 0), column(&stiefel_crossed, 1));
    assert!(
        observed_cross_gram != 0.0 && observed_cross_gram.abs() <= 1.0e-10,
        "fixture must exercise an admitted nonzero off-diagonal Gram entry: requested={cross_gram:.17e}; observed={observed_cross_gram:.17e}; point={stiefel_crossed:?}"
    );
    let crossed_normal = STIEFEL
        .parameter_gradient(&stiefel_crossed, &stiefel_crossed)
        .expect("generalized projection of a frame-normal gradient");
    let crossed_normal_scale = crossed_normal
        .iter()
        .map(|value| value.abs())
        .fold(0.0_f64, f64::max);
    assert!(
        crossed_normal_scale <= 1.0e-14,
        "generalized projection must remove a normal gradient at an accepted nonorthonormal frame within its declared resolution: point={stiefel_crossed:?}; projected={crossed_normal:?}; max={crossed_normal_scale:.17e}"
    );
    let crossed_first = STIEFEL
        .parameter_gradient(&stiefel_crossed, &STIEFEL_AMBIENT_GRADIENT)
        .expect("generalized projection at an accepted nonorthonormal frame");
    let crossed_second = STIEFEL
        .parameter_gradient(&stiefel_crossed, &crossed_first)
        .expect("generalized near-frame projection idempotence");
    let crossed_idempotence_error = max_error(&crossed_first, &crossed_second);
    let crossed_tangency_residual = stiefel_tangent_residual(&stiefel_crossed, &crossed_first);
    assert!(
        crossed_idempotence_error <= 2.0e-12 * (1.0 + norm(&crossed_first)),
        "generalized projection must be idempotent at an accepted nonorthonormal frame: first={crossed_first:?}; second={crossed_second:?}; error={crossed_idempotence_error:.17e}"
    );
    assert!(
        crossed_tangency_residual <= 6.4e-9 * norm(&crossed_first),
        "generalized projection must satisfy the near-frame tangent equation: point={stiefel_crossed:?}; projected={crossed_first:?}; residual={crossed_tangency_residual:.17e}"
    );

    let (direct_projection, direct_correction) =
        independent_stiefel_sylvester_projection(&stiefel_crossed, &STIEFEL_AMBIENT_GRADIENT);
    let direct_error = max_error(&crossed_first, &direct_projection);
    let direct_tangency = stiefel_tangent_residual(&stiefel_crossed, &direct_projection);
    let direct_gram = [
        [
            dot(column(&stiefel_crossed, 0), column(&stiefel_crossed, 0)),
            dot(column(&stiefel_crossed, 0), column(&stiefel_crossed, 1)),
        ],
        [
            dot(column(&stiefel_crossed, 1), column(&stiefel_crossed, 0)),
            dot(column(&stiefel_crossed, 1), column(&stiefel_crossed, 1)),
        ],
    ];
    let direct_rhs = [
        [
            2.0 * dot(
                column(&stiefel_crossed, 0),
                column(&STIEFEL_AMBIENT_GRADIENT, 0),
            ),
            dot(
                column(&stiefel_crossed, 0),
                column(&STIEFEL_AMBIENT_GRADIENT, 1),
            ) + dot(
                column(&STIEFEL_AMBIENT_GRADIENT, 0),
                column(&stiefel_crossed, 1),
            ),
        ],
        [
            dot(
                column(&stiefel_crossed, 1),
                column(&STIEFEL_AMBIENT_GRADIENT, 0),
            ) + dot(
                column(&STIEFEL_AMBIENT_GRADIENT, 1),
                column(&stiefel_crossed, 0),
            ),
            2.0 * dot(
                column(&stiefel_crossed, 1),
                column(&STIEFEL_AMBIENT_GRADIENT, 1),
            ),
        ],
    ];
    let mut sylvester_residual = 0.0_f64;
    let mut rhs_scale = 0.0_f64;
    for row in 0..P {
        for column in 0..P {
            let lhs: f64 = (0..P)
                .map(|basis| {
                    direct_gram[row][basis] * direct_correction[basis][column]
                        + direct_correction[row][basis] * direct_gram[basis][column]
                })
                .sum();
            sylvester_residual = sylvester_residual.max((lhs - direct_rhs[row][column]).abs());
            rhs_scale = rhs_scale.max(direct_rhs[row][column].abs());
        }
    }
    assert!(
        direct_error <= 2.0e-14 * (1.0 + norm(&direct_projection)),
        "Richardson projection must match the independent closed-form symmetric-Sylvester oracle: point={stiefel_crossed:?}; production={crossed_first:?}; oracle={direct_projection:?}; S={direct_correction:?}; error={direct_error:.17e}"
    );
    assert!(
        sylvester_residual <= 64.0 * f64::EPSILON * (1.0 + rhs_scale),
        "closed-form oracle must independently satisfy M*S + S*M = X^T*G + G^T*X: M={direct_gram:?}; S={direct_correction:?}; rhs={direct_rhs:?}; residual={sylvester_residual:.17e}"
    );
    assert!(
        direct_tangency <= 2.0e-15 * (1.0 + norm(&direct_projection)),
        "independent Sylvester projection must be orthogonal to the frame-normal space: point={stiefel_crossed:?}; oracle={direct_projection:?}; S={direct_correction:?}; residual={direct_tangency:.17e}"
    );

    let stiefel_landed = STIEFEL
        .retract(&STIEFEL_BASE, &scaled(&STIEFEL_STEP, 0.17))
        .expect("nontrivial Stiefel retraction");
    let stiefel_normal = STIEFEL
        .parameter_gradient(&stiefel_landed, &stiefel_landed)
        .expect("Stiefel normal projection at retraction-produced point");
    let normal_scale = stiefel_normal
        .iter()
        .map(|value| value.abs())
        .fold(0.0_f64, f64::max);
    assert!(
        normal_scale <= 1.0e-14,
        "frame-normal Stiefel gradient must collapse within the declared numerical envelope: point={stiefel_landed:?}; projected={stiefel_normal:?}; max={normal_scale:.17e}"
    );

    let stiefel_first = STIEFEL
        .parameter_gradient(&stiefel_landed, &STIEFEL_AMBIENT_GRADIENT)
        .expect("generalized Stiefel projection");
    let stiefel_second = STIEFEL
        .parameter_gradient(&stiefel_landed, &stiefel_first)
        .expect("generalized Stiefel projection idempotence");
    let stiefel_idempotence_error = max_error(&stiefel_first, &stiefel_second);
    assert!(
        stiefel_idempotence_error <= 2.0e-12 * (1.0 + norm(&stiefel_first)),
        "Stiefel projection must be idempotent at a retraction-produced frame: first={stiefel_first:?}; second={stiefel_second:?}; error={stiefel_idempotence_error:.17e}"
    );
    assert!(
        stiefel_tangent_residual(&stiefel_landed, &stiefel_first) <= 6.4e-9,
        "generalized projection must land tangent: point={stiefel_landed:?}; projected={stiefel_first:?}; residual={:.17e}",
        stiefel_tangent_residual(&stiefel_landed, &stiefel_first)
    );
}

#[test]
fn g3_so3_parameter_gradient_matches_an_independent_scalar_objective_difference() {
    let expected = Manifold::So3
        .parameter_gradient(&SO3_BASE, &SO3_AMBIENT_GRADIENT)
        .expect("SO(3) parameter gradient");
    let mut coarse_error = 0.0_f64;
    let mut fine_error = 0.0_f64;

    for parameter in 0..3 {
        let directional_difference = |h: f64| {
            let mut plus_step = [0.0; 3];
            let mut minus_step = [0.0; 3];
            plus_step[parameter] = h;
            minus_step[parameter] = -h;
            let plus = Manifold::So3
                .retract(&SO3_BASE, &plus_step)
                .expect("positive SO(3) scalar-objective sample");
            let minus = Manifold::So3
                .retract(&SO3_BASE, &minus_step)
                .expect("negative SO(3) scalar-objective sample");
            (dot(&SO3_AMBIENT_GRADIENT, &plus) - dot(&SO3_AMBIENT_GRADIENT, &minus)) / (2.0 * h)
        };
        let coarse = directional_difference(1.0 / 1024.0);
        let fine = directional_difference(1.0 / 2048.0);
        coarse_error = coarse_error.max((coarse - expected[parameter]).abs());
        fine_error = fine_error.max((fine - expected[parameter]).abs());
    }

    assert!(
        fine_error < 0.27 * coarse_error && fine_error <= 2.0e-7,
        "independent scalar-objective differences must converge to the declared right/body pullback: expected={expected:?}; coarse_error={coarse_error:.17e}; fine_error={fine_error:.17e}"
    );
}

#[test]
fn g0_extreme_finite_sphere_and_stiefel_gradients_project_without_intermediate_overflow() {
    let third = 0.28_f64.sqrt();
    let point = [0.6, 0.6, third];
    let ambient = [f64::MAX, f64::MAX, -0.5 * f64::MAX];
    let scaled_ambient = [1.0, 1.0, -0.5];
    let norm_sq = dot(&point, &point);
    let normal = dot(&point, &scaled_ambient) / norm_sq;
    let expected_scaled: [f64; 3] =
        core::array::from_fn(|index| scaled_ambient[index] - normal * point[index]);
    let expected: [f64; 3] = core::array::from_fn(|index| expected_scaled[index] * f64::MAX);
    assert!(
        expected.iter().all(|value| value.is_finite()),
        "fixture must have a finite representable projection even though its naive first two dot terms overflow: point={point:?}; scaled_expected={expected_scaled:?}; expected={expected:?}"
    );

    let sphere = Manifold::Sphere { ambient: 3 };
    let sphere_projection = sphere
        .parameter_gradient(&point, &ambient)
        .expect("scale-safe Sphere projection");
    sphere
        .validate_parameter_tangent(&point, &sphere_projection)
        .expect("extreme Sphere projection tangent postcondition");
    let sphere_relative_error = max_error(&sphere_projection, &expected) / f64::MAX;
    assert!(
        sphere_projection.iter().all(|value| value.is_finite()) && sphere_relative_error <= 2.0e-15,
        "Sphere projection must avoid intermediate overflow and match the independently scaled oracle: point={point:?}; ambient={ambient:?}; output={sphere_projection:?}; expected={expected:?}; relative_error={sphere_relative_error:.17e}"
    );

    let stiefel = Manifold::Stiefel { n: 3, p: 1 };
    let stiefel_projection = stiefel
        .parameter_gradient(&point, &ambient)
        .expect("scale-safe Stiefel(n,1) projection");
    stiefel
        .validate_parameter_tangent(&point, &stiefel_projection)
        .expect("extreme Stiefel projection tangent postcondition");
    let stiefel_relative_error = max_error(&stiefel_projection, &expected) / f64::MAX;
    assert!(
        stiefel_projection.iter().all(|value| value.is_finite())
            && stiefel_relative_error <= 2.0e-15,
        "Stiefel(n,1) projection must avoid intermediate overflow and match the independently scaled oracle: point={point:?}; ambient={ambient:?}; output={stiefel_projection:?}; expected={expected:?}; relative_error={stiefel_relative_error:.17e}"
    );

    let nonrepresentable_point = [0.5, 0.75_f64.sqrt()];
    let nonrepresentable = (Manifold::Sphere { ambient: 2 })
        .parameter_gradient(&nonrepresentable_point, &[f64::MAX, -f64::MAX]);
    assert!(
        matches!(
            nonrepresentable,
            Err(OptError::RetractionNonFinite {
                input: "parameter gradient output",
                component: 0,
                bits,
            }) if bits == f64::INFINITY.to_bits()
        ),
        "a genuinely nonrepresentable projected component must refuse at its exact output lane rather than overflow internally or canonicalize to zero: point={nonrepresentable_point:?}; result={nonrepresentable:?}"
    );
}

#[test]
fn g0_compensated_tangency_reduction_accepts_high_dimensional_exact_cancellation() {
    // This 2^20 fixture is the smallest convenient power-of-two version of
    // the max-admission counterexample. It retains the same four-quarter
    // cancellation, but uses only 16 MiB for the two vectors rather than the
    // 256 MiB required by two 2^24-lane f64 payloads.
    const AMBIENT: usize = 1 << 20;
    const QUARTER: usize = AMBIENT / 4;
    let point_coordinate = 1.0 / 1024.0;
    let tiny = 1.0 / (1_u64 << 35) as f64;
    let point = vec![point_coordinate; AMBIENT];
    let mut tangent = vec![1.0; AMBIENT];
    tangent[QUARTER..2 * QUARTER].fill(tiny);
    tangent[2 * QUARTER..3 * QUARTER].fill(-1.0);
    tangent[3 * QUARTER..].fill(-tiny);

    let naive_residual = dot(&point, &tangent).abs();
    assert!(
        naive_residual > 6.4e-9,
        "fixture must independently reproduce the ordered left-fold false rejection: ambient={AMBIENT}; point_coordinate={point_coordinate:.17e}; tiny={tiny:.17e}; naive_residual={naive_residual:.17e}"
    );
    (Manifold::Sphere {
        ambient: AMBIENT as u32,
    })
    .validate_parameter_tangent(&point, &tangent)
    .expect("compensated reduction must accept the exactly cancelling tangent");
    (Manifold::Stiefel {
        n: AMBIENT as u32,
        p: 1,
    })
    .validate_parameter_tangent(&point, &tangent)
    .expect("compensated Gram reduction must accept the exactly cancelling Stiefel tangent");
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
fn g3_sphere_transport_normalizes_accepted_bases_and_fails_closed_near_antipodes() {
    let sphere = Manifold::Sphere { ambient: 3 };
    let scale = accepted_near_unit_scale();
    let from = [scale, 0.0, 0.0];
    let tangent = [0.0, 1.0, 0.0];

    let admitted_gap = 0.02_f64;
    let admitted_target = [-admitted_gap.cos(), admitted_gap.sin(), 0.0];
    let admitted_step: [f64; 3] =
        core::array::from_fn(|index| admitted_target[index] - from[index]);
    let admitted_to = sphere
        .retract(&from, &admitted_step)
        .expect("conditioned near-antipodal destination");
    let transported = sphere
        .transport_parameter(&from, &admitted_step, &admitted_to, &tangent)
        .expect("normalized-base Sphere transport above conditioning floor");
    let norm_error = (norm(&transported) - norm(&tangent)).abs();
    assert!(
        norm_error <= 2.0e-8,
        "accepted near-antipodal transport must preserve norm: from={from:?}; to={admitted_to:?}; input={tangent:?}; output={transported:?}; error={norm_error:.17e}"
    );
    assert!(
        dot(&admitted_to, &transported).abs() <= 6.4e-9 * norm(&transported),
        "accepted near-antipodal transport must land tangent: to={admitted_to:?}; output={transported:?}; residual={:.17e}",
        dot(&admitted_to, &transported).abs()
    );

    let refused_gap = 1.0e-4_f64;
    let refused_target = [-refused_gap.cos(), refused_gap.sin(), 0.0];
    let refused_step: [f64; 3] = core::array::from_fn(|index| refused_target[index] - from[index]);
    let refused_to = sphere
        .retract(&from, &refused_step)
        .expect("ill-conditioned near-antipodal destination remains a valid point");
    let refusal = sphere.transport_parameter(&from, &refused_step, &refused_to, &tangent);
    assert!(
        matches!(
            refusal,
            Err(OptError::RetractionDomain {
                manifold: "Sphere",
                what: "Sphere transport is undefined or ill-conditioned near antipodal points",
                location: None,
                ..
            })
        ),
        "ill-conditioned near-antipodal transport must fail closed: from={from:?}; to={refused_to:?}; result={refusal:?}"
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

#[test]
fn g0_invalid_descriptors_and_expanded_transport_payloads_fail_with_attribution() {
    assert!(matches!(
        (Manifold::Sphere { ambient: 1 }).parameter_gradient(&[1.0], &[1.0]),
        Err(OptError::ManifoldInvalid { .. })
    ));
    assert!(matches!(
        (Manifold::Stiefel { n: 2, p: 3 }).parameter_gradient(&[1.0; 6], &[0.0; 6]),
        Err(OptError::ManifoldInvalid { .. })
    ));

    let sphere = Manifold::Sphere { ambient: 3 };
    let from = [1.0, 0.0, 0.0];
    let step = [0.0, 0.25, 0.0];
    let to = sphere.retract(&from, &step).expect("valid destination");
    assert!(matches!(
        sphere.transport_parameter(&from, &step[..2], &to, &[0.0, 1.0, 0.0]),
        Err(OptError::RetractionLen {
            input: "transport retraction step",
            expected: 3,
            got: 2,
        })
    ));
    assert!(matches!(
        sphere.transport_parameter(&from, &step, &to, &[0.0, 1.0]),
        Err(OptError::RetractionLen {
            input: "manifold tangent parameter",
            expected: 3,
            got: 2,
        })
    ));
    let nonfinite_vector = [0.0, f64::NEG_INFINITY, 0.0];
    assert!(matches!(
        sphere.transport_parameter(&from, &step, &to, &nonfinite_vector),
        Err(OptError::RetractionNonFinite {
            input: "manifold tangent parameter",
            component: 1,
            bits,
        }) if bits == f64::NEG_INFINITY.to_bits()
    ));
    assert!(matches!(
        sphere.transport_parameter(&from, &step, &to[..2], &[0.0, 1.0, 0.0]),
        Err(OptError::RetractionLen {
            input: "manifold operation point",
            expected: 3,
            got: 2,
        })
    ));
}

#[test]
fn g0_seeded_malformed_payloads_reject_deterministically_with_exact_locations() {
    const SEED: u64 = 0x7A22_0711_5EED_0001;
    let mut state = SEED;
    let sphere = Manifold::Sphere { ambient: 3 };
    let point = [1.0, 0.0, 0.0];

    for case in 0..64_u32 {
        let sample = xorshift64(&mut state);
        let lane = (sample as usize) % 3;
        let payload_bits = 0x7ff8_0000_0000_0000_u64 | (sample & 0x0007_ffff_ffff_ffff);
        let payload = f64::from_bits(payload_bits);
        let mut direction = [0.25, -0.5, 0.75];
        direction[lane] = payload;
        let result = sphere.retract_curve(&point, &direction, 0.5);
        assert!(
            matches!(
                result,
                Err(OptError::RetractionNonFinite {
                    input: "retraction curve direction",
                    component,
                    bits,
                }) if component as usize == lane && bits == payload_bits
            ),
            "seeded nonfinite case must retain exact location and payload: seed={SEED:#018x}; case={case}; lane={lane}; bits={payload_bits:#018x}; result={result:?}"
        );

        let excess = 1.0e-4 * (1.0 + ((sample >> 32) & 0xff) as f64);
        let off_manifold = [1.0 + excess, 0.0, 0.0];
        let domain_result = sphere.parameter_gradient(&off_manifold, &[0.0, 1.0, 0.0]);
        assert!(
            matches!(
                domain_result,
                Err(OptError::RetractionDomain {
                    manifold: "Sphere",
                    what: "point must have unit norm before retraction",
                    location: None,
                    ..
                })
            ),
            "seeded off-manifold case must fail through the point-domain rule: seed={SEED:#018x}; case={case}; excess={excess:.17e}; result={domain_result:?}"
        );
    }
}
