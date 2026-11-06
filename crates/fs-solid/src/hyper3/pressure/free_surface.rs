//! Exact positive-part moments in the original simplex barycentric chart.
//!
//! A plane cuts a triangle into at most two wet triangles and a tetrahedron
//! into at most three wet tetrahedra. Positive decompositions avoid subtracting
//! an almost equal dry volume or dividing a potential by a small gradient.

use super::{HyperTetError, finite};

/// The positive simplex formula uses only rational arithmetic. Keeping this
/// small scalar seam lets the native AD and exact interval divided differences
/// differentiate the SAME clipped-volume moment, including very thin wet caps.
pub(super) trait MomentScalar:
    Copy
    + core::ops::Add<Output = Self>
    + core::ops::Sub<Output = Self>
    + core::ops::Mul<Output = Self>
    + core::ops::Div<Output = Self>
    + core::ops::Neg<Output = Self>
{
    fn constant(value: f64) -> Self;
    fn value(self) -> f64;
}

impl<R: fs_ad::Real> MomentScalar for R {
    fn constant(value: f64) -> Self {
        Self::from_f64(value)
    }
    fn value(self) -> f64 {
        fs_ad::Real::value(self)
    }
}

const FULL_BASIS: [[f64; 3]; 3] = [
    [1.0 / 6.0, 1.0 / 12.0, 1.0 / 12.0],
    [1.0 / 12.0, 1.0 / 6.0, 1.0 / 12.0],
    [1.0 / 12.0, 1.0 / 12.0, 1.0 / 6.0],
];

pub(super) struct FaceMoments {
    /// `integral N_i p dA / full_face_area`.
    pub weights_pa: [f64; 3],
    /// Derivative of each weighted moment with respect to the nodal scalar.
    /// For the positive linear moment this is integral_wet N_i N_j / area;
    /// for its square it is integral_wet 2 N_i N_j q / area.
    pub basis: [[f64; 3]; 3],
}

impl FaceMoments {
    pub fn constant(pressure_pa: f64) -> Self {
        Self {
            weights_pa: [pressure_pa / 3.0; 3],
            basis: FULL_BASIS,
        }
    }

    pub fn affine(values: [f64; 3]) -> Result<Self, HyperTetError> {
        let mut weights_pa = [0.0; 3];
        for local in 0..3 {
            weights_pa[local] = finite(
                values[local] / 6.0
                    + values[(local + 1) % 3] / 12.0
                    + values[(local + 2) % 3] / 12.0,
            )?;
        }
        Ok(Self {
            weights_pa,
            basis: FULL_BASIS,
        })
    }

    fn dry() -> Self {
        Self {
            weights_pa: [0.0; 3],
            basis: [[0.0; 3]; 3],
        }
    }
}

/// Fractions from a wet vertex to the plane and from the plane to a dry one.
/// Rescale before adding opposite pressure magnitudes; calculate both fractions
/// directly so a thin wet region does not become `1 - nearly_one`.
fn edge_fractions<R: MomentScalar>(wet: R, dry: R) -> (R, R) {
    // This common scale cancels algebraically, so holding its primal constant
    // also gives the exact derivatives and divided differences of the ratio.
    let scale = R::constant(wet.value().max(-dry.value()));
    let positive = wet / scale;
    let negative = -dry / scale;
    let sum = positive + negative;
    (positive / sum, negative / sum)
}

fn vertex(node: usize) -> [f64; 3] {
    let mut barycentric = [0.0; 3];
    barycentric[node] = 1.0;
    barycentric
}

fn intersection(wet: usize, dry: usize, fractions: (f64, f64)) -> [f64; 3] {
    let mut barycentric = [0.0; 3];
    barycentric[wet] = fractions.1;
    barycentric[dry] = fractions.0;
    barycentric
}

/// Integrate products of two P1 fields on one wet subtriangle. All pressure
/// values and geometric coefficients are nonnegative. TANGENT=false removes
/// the basis-product work from a force-only evaluation.
fn add_triangle<const TANGENT: bool, const SQUARED: bool>(
    result: &mut FaceMoments,
    area_fraction: f64,
    basis: [[f64; 3]; 3],
    pressure: [f64; 3],
) -> Result<(), HyperTetError> {
    if area_fraction == 0.0 {
        return Ok(());
    }
    if SQUARED {
        return add_square_triangle::<TANGENT>(result, area_fraction, basis, pressure);
    }
    let sums: [f64; 3] = std::array::from_fn(|i| basis[0][i] + basis[1][i] + basis[2][i]);
    let scaled_pressure = pressure.map(|value| value / 12.0);
    let pressure_sum = scaled_pressure.into_iter().sum::<f64>();
    for row in 0..3 {
        let diagonal = (0..3)
            .map(|vertex| basis[vertex][row] * scaled_pressure[vertex])
            .sum::<f64>();
        let contribution = finite(area_fraction * finite(sums[row] * pressure_sum + diagonal)?)?;
        result.weights_pa[row] = finite(result.weights_pa[row] + contribution)?;
        if TANGENT {
            for column in 0..3 {
                let diagonal = (0..3)
                    .map(|vertex| basis[vertex][row] * basis[vertex][column])
                    .sum::<f64>();
                let contribution = area_fraction * (sums[row] * sums[column] + diagonal) / 12.0;
                result.basis[row][column] = finite(result.basis[row][column] + contribution)?;
            }
        }
    }
    Ok(())
}

pub(super) fn face_moments<const TANGENT: bool>(
    values: [f64; 3],
) -> Result<FaceMoments, HyperTetError> {
    face_power_moments::<TANGENT, false>(values)
}

/// Exact weighted positive quadratic field and its first derivative.
pub(super) fn face_square_moments<const TANGENT: bool>(
    values: [f64; 3],
) -> Result<FaceMoments, HyperTetError> {
    face_power_moments::<TANGENT, true>(values)
}

/// Exact cubic triangle moment from the Dirichlet(1,1,1) barycentric law.
/// All wet values and weights are nonnegative; no pressure interpolation
/// replaces the quadratic field. Its derivative uses the same cubic rule.
fn add_square_triangle<const TANGENT: bool>(
    result: &mut FaceMoments,
    area_fraction: f64,
    basis: [[f64; 3]; 3],
    values: [f64; 3],
) -> Result<(), HyperTetError> {
    let sum = values.iter().sum::<f64>();
    let square_sum = values.iter().map(|value| value * value).sum::<f64>();
    let sums: [f64; 3] = core::array::from_fn(|i| basis.iter().map(|v| v[i]).sum());
    let mixed: [f64; 3] = core::array::from_fn(|i| (0..3).map(|v| basis[v][i] * values[v]).sum());
    for row in 0..3 {
        let diagonal = (0..3)
            .map(|v| basis[v][row] * values[v] * values[v])
            .sum::<f64>();
        let mean = finite(
            (sums[row] * (sum * sum + square_sum) + 2.0 * mixed[row] * sum + 2.0 * diagonal) / 60.0,
        )?;
        result.weights_pa[row] = finite(result.weights_pa[row] + area_fraction * mean)?;
        if TANGENT {
            for column in 0..3 {
                let product = (0..3)
                    .map(|v| basis[v][row] * basis[v][column])
                    .sum::<f64>();
                let diagonal = (0..3)
                    .map(|v| basis[v][row] * basis[v][column] * values[v])
                    .sum::<f64>();
                let mean = finite(
                    (sums[row] * sums[column] * sum
                        + product * sum
                        + mixed[row] * sums[column]
                        + mixed[column] * sums[row]
                        + 2.0 * diagonal)
                        / 30.0,
                )?;
                result.basis[row][column] =
                    finite(result.basis[row][column] + area_fraction * mean)?;
            }
        }
    }
    Ok(())
}

fn face_power_moments<const TANGENT: bool, const SQUARED: bool>(
    values: [f64; 3],
) -> Result<FaceMoments, HyperTetError> {
    let wet_count = values.iter().filter(|&&value| value > 0.0).count();
    if wet_count == 0 {
        // On an entirely coplanar face, choose the dry-side generalized
        // derivative. Isolated zero-pressure edges/vertices have zero measure.
        return Ok(FaceMoments::dry());
    }
    if values.iter().all(|&value| value >= 0.0) {
        if !SQUARED {
            return FaceMoments::affine(values);
        }
        let mut result = FaceMoments::dry();
        add_square_triangle::<TANGENT>(
            &mut result,
            1.0,
            [vertex(0), vertex(1), vertex(2)],
            values,
        )?;
        return Ok(result);
    }
    let mut result = FaceMoments::dry();
    if wet_count == 1 {
        let wet = values.iter().position(|&value| value > 0.0).unwrap();
        let first = (wet + 1) % 3;
        let second = (wet + 2) % 3;
        let a = edge_fractions(values[wet], values[first]);
        let b = edge_fractions(values[wet], values[second]);
        add_triangle::<TANGENT, SQUARED>(
            &mut result,
            a.0 * b.0,
            [
                vertex(wet),
                intersection(wet, first, a),
                intersection(wet, second, b),
            ],
            [values[wet], 0.0, 0.0],
        )?;
    } else {
        let dry = values.iter().position(|&value| value < 0.0).unwrap();
        let first = (dry + 1) % 3;
        let second = (dry + 2) % 3;
        let a = edge_fractions(values[first], values[dry]);
        let b = edge_fractions(values[second], values[dry]);
        let cut_a = intersection(first, dry, a);
        let cut_b = intersection(second, dry, b);
        // Fan of the wet quadrilateral [first, second, cut_b, cut_a].
        // These area ratios are products, not small determinant differences.
        add_triangle::<TANGENT, SQUARED>(
            &mut result,
            b.0,
            [vertex(first), vertex(second), cut_b],
            [values[first], values[second], 0.0],
        )?;
        add_triangle::<TANGENT, SQUARED>(
            &mut result,
            a.0 * b.1,
            [vertex(first), cut_b, cut_a],
            [values[first], 0.0, 0.0],
        )?;
    }
    Ok(result)
}

/// `integral_tet max(0,p) dV / full_tet_volume`. A clipped simplex is split
/// into positive-volume tetrahedra; the average of an affine field on each
/// tetrahedron is the arithmetic mean of its four vertex values.
pub(super) fn tet_mean(values: [f64; 4]) -> Result<f64, HyperTetError> {
    finite(tet_mean_scalar(values))
}

pub(super) fn tet_mean_scalar<R: MomentScalar>(values: [R; 4]) -> R {
    tet_power_mean_scalar::<R, false>(values)
}

pub(super) fn tet_square_mean(values: [f64; 4]) -> Result<f64, HyperTetError> {
    finite(tet_square_mean_scalar(values))
}

/// Positive quadratic moment on the SAME admitted clipped tetrahedra.
pub(super) fn tet_square_mean_scalar<R: MomentScalar>(values: [R; 4]) -> R {
    tet_power_mean_scalar::<R, true>(values)
}

// mean(q^2) = sum_{i<=j} q_i q_j / 10 on a tetrahedron. Scale
// each nonnegative product before summing; retain the original linear order.
fn full_mean<R: MomentScalar, const SQUARED: bool>(values: [R; 4]) -> R {
    if !SQUARED {
        return values.into_iter().fold(R::constant(0.0), |sum, value| {
            sum + value / R::constant(4.0)
        });
    }
    let mut mean = R::constant(0.0);
    for i in 0..4 {
        for j in i..4 {
            mean = mean + (values[i] / R::constant(10.0)) * values[j];
        }
    }
    mean
}

fn tet_power_mean_scalar<R: MomentScalar, const SQUARED: bool>(values: [R; 4]) -> R {
    if values.iter().all(|&value| value.value() <= 0.0) {
        // At an entirely coplanar simplex select the dry-side derivative.
        return R::constant(0.0);
    }
    let quarter = R::constant(4.0);
    if values.iter().all(|&value| value.value() >= 0.0) {
        return full_mean::<R, SQUARED>(values);
    }
    let mut wet = [0; 4];
    let mut dry = [0; 4];
    let (mut wet_count, mut dry_count) = (0, 0);
    for (index, value) in values.into_iter().enumerate() {
        if value.value() > 0.0 {
            wet[wet_count] = index;
            wet_count += 1;
        } else {
            dry[dry_count] = index;
            dry_count += 1;
        }
    }
    match wet_count {
        0 => R::constant(0.0),
        1 => {
            let p = values[wet[0]];
            let fraction = dry[..3].iter().fold(R::constant(1.0), |product, &index| {
                product * edge_fractions(p, values[index]).0
            });
            let mean = if SQUARED {
                (p / R::constant(10.0)) * p
            } else {
                p / quarter
            };
            mean * fraction
        }
        2 => {
            let (a, b, c, d) = (
                values[wet[0]],
                values[wet[1]],
                values[dry[0]],
                values[dry[1]],
            );
            let ac = edge_fractions(a, c);
            let ad = edge_fractions(a, d);
            let bc = edge_fractions(b, c);
            let bd = edge_fractions(b, d);
            // [A, AC, AD, BD], [A, AC, BC, BD], [A, B, BC, BD].
            let first = ac.0 * ad.0 * bd.1;
            let second = ac.0 * bc.1 * bd.0;
            let third = bc.0 * bd.0;
            if SQUARED {
                let a_mean = (a / R::constant(10.0)) * a;
                let ab_mean = full_mean::<R, true>([a, b, R::constant(0.0), R::constant(0.0)]);
                a_mean * first + a_mean * second + ab_mean * third
            } else {
                (a / quarter) * first + (a / quarter) * second + (a / quarter + b / quarter) * third
            }
        }
        3 => {
            let (a, b, c, d) = (
                values[wet[0]],
                values[wet[1]],
                values[wet[2]],
                values[dry[0]],
            );
            let ad = edge_fractions(a, d);
            let bd = edge_fractions(b, d);
            let cd = edge_fractions(c, d);
            // [A, B, C, CD], [A, B, BD, CD], [A, AD, BD, CD].
            if SQUARED {
                let abc = full_mean::<R, true>([a, b, c, R::constant(0.0)]);
                let ab = full_mean::<R, true>([a, b, R::constant(0.0), R::constant(0.0)]);
                let aa = (a / R::constant(10.0)) * a;
                abc * cd.0 + ab * cd.1 * bd.0 + aa * bd.1 * cd.1 * ad.0
            } else {
                (a / quarter + b / quarter + c / quarter) * cd.0
                    + (a / quarter + b / quarter) * cd.1 * bd.0
                    + (a / quarter) * bd.1 * cd.1 * ad.0
            }
        }
        _ => unreachable!("fully wet tetrahedron handled above"),
    }
}
