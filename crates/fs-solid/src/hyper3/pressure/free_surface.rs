//! Exact positive-part moments in the original simplex barycentric chart.
//!
//! A plane cuts a triangle into at most two wet triangles and a tetrahedron
//! into at most three wet tetrahedra. Positive decompositions avoid subtracting
//! an almost equal dry volume or dividing a potential by a small gradient.

use super::{HyperTetError, finite};

const FULL_BASIS: [[f64; 3]; 3] = [
    [1.0 / 6.0, 1.0 / 12.0, 1.0 / 12.0],
    [1.0 / 12.0, 1.0 / 6.0, 1.0 / 12.0],
    [1.0 / 12.0, 1.0 / 12.0, 1.0 / 6.0],
];

pub(super) struct FaceMoments {
    /// `integral N_i p dA / full_face_area`.
    pub weights_pa: [f64; 3],
    /// `integral_wet N_i N_j dA / full_face_area`.
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
fn edge_fractions(wet: f64, dry: f64) -> (f64, f64) {
    let scale = wet.max(-dry);
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
fn add_triangle<const TANGENT: bool>(
    result: &mut FaceMoments,
    area_fraction: f64,
    basis: [[f64; 3]; 3],
    pressure: [f64; 3],
) -> Result<(), HyperTetError> {
    if area_fraction == 0.0 {
        return Ok(());
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
    let wet_count = values.iter().filter(|&&value| value > 0.0).count();
    if wet_count == 0 {
        // On an entirely coplanar face, choose the dry-side generalized
        // derivative. Isolated zero-pressure edges/vertices have zero measure.
        return Ok(FaceMoments::dry());
    }
    if values.iter().all(|&value| value >= 0.0) {
        return FaceMoments::affine(values);
    }
    let mut result = FaceMoments::dry();
    if wet_count == 1 {
        let wet = values.iter().position(|&value| value > 0.0).unwrap();
        let first = (wet + 1) % 3;
        let second = (wet + 2) % 3;
        let a = edge_fractions(values[wet], values[first]);
        let b = edge_fractions(values[wet], values[second]);
        add_triangle::<TANGENT>(
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
        add_triangle::<TANGENT>(
            &mut result,
            b.0,
            [vertex(first), vertex(second), cut_b],
            [values[first], values[second], 0.0],
        )?;
        add_triangle::<TANGENT>(
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
    if values.iter().all(|&value| value >= 0.0) {
        return finite(values.into_iter().map(|value| value / 4.0).sum());
    }
    let mut wet = [0; 4];
    let mut dry = [0; 4];
    let (mut wet_count, mut dry_count) = (0, 0);
    for (index, value) in values.into_iter().enumerate() {
        if value > 0.0 {
            wet[wet_count] = index;
            wet_count += 1;
        } else {
            dry[dry_count] = index;
            dry_count += 1;
        }
    }
    let mean = match wet_count {
        0 => 0.0,
        1 => {
            let p = values[wet[0]];
            let fraction = dry[..3]
                .iter()
                .map(|&index| edge_fractions(p, values[index]).0)
                .product::<f64>();
            (p / 4.0) * fraction
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
            (a / 4.0) * first + (a / 4.0) * second + (a / 4.0 + b / 4.0) * third
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
            (a / 4.0 + b / 4.0 + c / 4.0) * cd.0
                + (a / 4.0 + b / 4.0) * cd.1 * bd.0
                + (a / 4.0) * bd.1 * cd.1 * ad.0
        }
        _ => unreachable!("fully wet tetrahedron handled above"),
    };
    finite(mean)
}
