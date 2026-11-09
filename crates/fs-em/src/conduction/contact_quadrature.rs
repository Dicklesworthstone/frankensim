//! Positive quadrature for source-sampled electrical contact traces.
//!
//! Edge composite Simpson integrates cubic polynomials exactly. Each triangle
//! subdivision uses positive vertex/midpoint/centroid weights (1/20, 2/15,
//! 9/20), also exact through degree three. These fixed rules are discretizations,
//! not error certificates for an arbitrary nonlinear interface response.

use std::collections::BTreeMap;

#[derive(Clone, Debug, PartialEq)]
pub(crate) enum QuadratureError {
    Invalid(&'static str),
    Cancelled,
}

#[derive(Clone, Debug)]
pub(crate) enum ContactCoefficient {
    Uniform(f64),
    Sampled {
        subdivisions: usize,
        resistance_ohm_m2: Vec<f64>,
    },
}

/// An actual integration point of the admitted source-sampled contact law.
/// Basis coordinates refer to the original paired trace vertices. Surface
/// weights already include planar extrusion or axisymmetric revolution.
#[derive(Clone, Debug, PartialEq)]
pub struct ConductionContactQuadraturePoint<const N: usize> {
    pub basis: [f64; N],
    pub area_weight_m2: f64,
    pub voltage_jump_v: f64,
    pub current_density_a_m2: f64,
    pub resistance_ohm_m2: f64,
    pub power_w: f64,
}

#[derive(Clone, Debug)]
pub(crate) struct RulePoint<const N: usize> {
    pub basis: [f64; N],
    pub weight: f64,
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), QuadratureError> {
    if cancelled() {
        Err(QuadratureError::Cancelled)
    } else {
        Ok(())
    }
}

fn finite(value: f64) -> Result<f64, QuadratureError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(QuadratureError::Invalid(
            "contact quadrature arithmetic is not representable",
        ))
    }
}

fn add(sum: &mut f64, correction: &mut f64, value: f64) -> Result<(), QuadratureError> {
    let adjusted = finite(value - *correction)?;
    let next = finite(*sum + adjusted)?;
    *correction = finite((next - *sum) - adjusted)?;
    *sum = next;
    Ok(())
}

pub(crate) fn edge_point_count(subdivisions: usize) -> Result<usize, QuadratureError> {
    subdivisions
        .checked_mul(2)
        .and_then(|n| n.checked_add(1))
        .filter(|_| subdivisions > 0)
        .ok_or(QuadratureError::Invalid(
            "positive bounded contact subdivisions required",
        ))
}

pub(crate) fn triangle_point_count(subdivisions: usize) -> Result<usize, QuadratureError> {
    subdivisions
        .checked_mul(subdivisions)
        .and_then(|n| n.checked_mul(3))
        .and_then(|n| n.checked_add(subdivisions.checked_mul(3)?))
        .and_then(|n| n.checked_add(1))
        .filter(|_| subdivisions > 0 && subdivisions.checked_mul(6).is_some())
        .ok_or(QuadratureError::Invalid(
            "positive bounded triangular contact subdivisions required",
        ))
}

pub(crate) fn edge_rule(
    subdivisions: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<RulePoint<2>>, QuadratureError> {
    poll(cancelled)?;
    let count = edge_point_count(subdivisions)?;
    let intervals = count - 1;
    let mut points = Vec::with_capacity(count);
    for i in 0..count {
        poll(cancelled)?;
        let t = i as f64 / intervals as f64;
        let factor = if i == 0 || i == intervals {
            1.0
        } else if i % 2 == 1 {
            4.0
        } else {
            2.0
        };
        points.push(RulePoint {
            basis: [1.0 - t, t],
            weight: factor / (3.0 * intervals as f64),
        });
    }
    Ok(points)
}

pub(crate) fn triangle_rule(
    subdivisions: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<RulePoint<3>>, QuadratureError> {
    poll(cancelled)?;
    let count = triangle_point_count(subdivisions)?;
    let denominator = 6 * subdivisions;
    let cell_weight = 1.0 / subdivisions as f64 / subdivisions as f64;
    // Exact integer barycentric keys identify shared quadrature points before
    // converting to f64. Original element order fixes all weight reductions.
    let mut points = BTreeMap::<[usize; 3], f64>::new();
    let node = |i: usize, j: usize| [6 * i, 6 * j, 6 * (subdivisions - i - j)];
    let mut add_triangle = |vertices: [[usize; 3]; 3]| -> Result<(), QuadratureError> {
        for vertex in vertices {
            poll(cancelled)?;
            let sum = points.entry(vertex).or_insert(0.0);
            *sum = finite(*sum + cell_weight / 20.0)?;
        }
        for [i, j] in [[0, 1], [1, 2], [2, 0]] {
            poll(cancelled)?;
            let key = std::array::from_fn(|axis| vertices[i][axis] / 2 + vertices[j][axis] / 2);
            let sum = points.entry(key).or_insert(0.0);
            *sum = finite(*sum + cell_weight * (2.0 / 15.0))?;
        }
        let key = std::array::from_fn(|axis| {
            vertices[0][axis] / 3 + vertices[1][axis] / 3 + vertices[2][axis] / 3
        });
        let sum = points.entry(key).or_insert(0.0);
        *sum = finite(*sum + cell_weight * (9.0 / 20.0))?;
        Ok(())
    };
    for i in 0..subdivisions {
        for j in 0..subdivisions - i {
            add_triangle([node(i, j), node(i + 1, j), node(i, j + 1)])?;
            if i + j + 1 < subdivisions {
                add_triangle([node(i + 1, j), node(i + 1, j + 1), node(i, j + 1)])?;
            }
        }
    }
    if points.len() != count {
        return Err(QuadratureError::Invalid(
            "triangular contact quadrature lost its admitted point count",
        ));
    }
    let mut result = Vec::with_capacity(count);
    for (key, weight) in points {
        poll(cancelled)?;
        result.push(RulePoint {
            basis: key.map(|v| v as f64 / denominator as f64),
            weight,
        });
    }
    Ok(result)
}

#[derive(Clone)]
struct ContactPoint<const N: usize> {
    basis: [f64; N],
    area_m2: f64,
    resistance_ohm_m2: f64,
}

/// One coefficient sample per admitted integration point. The same positive
/// outer products supply assembly, current and irreversible surface work.
#[derive(Clone)]
pub(crate) struct ContactQuadrature<const N: usize> {
    points: Vec<ContactPoint<N>>,
    matrix_s: [[f64; N]; N],
}

impl<const N: usize> ContactQuadrature<N> {
    pub(crate) fn new(
        rule: Vec<RulePoint<N>>,
        area_m2: f64,
        endpoint_measure_factors: [f64; N],
        resistance_ohm_m2: Vec<f64>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, QuadratureError> {
        poll(cancelled)?;
        if N == 0
            || rule.is_empty()
            || rule.len() != resistance_ohm_m2.len()
            || !area_m2.is_finite()
            || area_m2 <= 0.0
            || endpoint_measure_factors
                .iter()
                .any(|v| !v.is_finite() || *v < 0.0)
            || endpoint_measure_factors.iter().all(|v| *v == 0.0)
        {
            return Err(QuadratureError::Invalid(
                "complete positive contact quadrature geometry required",
            ));
        }
        let mut points = Vec::with_capacity(rule.len());
        let mut matrix_s = [[0.0; N]; N];
        let mut correction = [[0.0; N]; N];
        for (point, resistance) in rule.into_iter().zip(resistance_ohm_m2) {
            poll(cancelled)?;
            if !resistance.is_finite()
                || resistance <= 0.0
                || !point.weight.is_finite()
                || point.weight <= 0.0
                || point
                    .basis
                    .iter()
                    .any(|v| !v.is_finite() || !(0.0..=1.0).contains(v))
                || (point.basis.iter().sum::<f64>() - 1.0).abs() > 16.0 * f64::EPSILON
            {
                return Err(QuadratureError::Invalid(
                    "positive source resistance and barycentric contact weights required",
                ));
            }
            let measure = finite(
                point
                    .basis
                    .iter()
                    .zip(endpoint_measure_factors)
                    .map(|(a, b)| a * b)
                    .sum(),
            )?;
            let area = finite(finite(area_m2 * point.weight)? * measure)?;
            // A trace endpoint on the rotation axis has exactly zero swept area.
            // All nonzero quadrature measures must retain representable area.
            if area < 0.0 || (measure > 0.0 && area == 0.0) {
                return Err(QuadratureError::Invalid(
                    "positive contact quadrature area is not representable",
                ));
            }
            let conductance = finite(area / resistance)?;
            if area > 0.0 && conductance <= 0.0 {
                return Err(QuadratureError::Invalid(
                    "positive contact quadrature conductance is not representable",
                ));
            }
            for i in 0..N {
                for j in i..N {
                    let value = finite(conductance * (point.basis[i] * point.basis[j]))?;
                    add(&mut matrix_s[i][j], &mut correction[i][j], value)?;
                    matrix_s[j][i] = matrix_s[i][j];
                }
            }
            points.push(ContactPoint {
                basis: point.basis,
                area_m2: area,
                resistance_ohm_m2: resistance,
            });
        }
        for i in 0..N {
            if matrix_s[i][i] <= 0.0 {
                return Err(QuadratureError::Invalid(
                    "contact quadrature lost a positive trace coordinate",
                ));
            }
        }
        Ok(Self { points, matrix_s })
    }

    pub(crate) fn matrix_entry(&self, i: usize, j: usize) -> f64 {
        self.matrix_s[i][j]
    }

    pub(crate) fn response(
        &self,
        jump: [f64; N],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<ConductionContactQuadraturePoint<N>>, QuadratureError> {
        poll(cancelled)?;
        if jump.iter().any(|v| !v.is_finite()) {
            return Err(QuadratureError::Invalid(
                "finite contact voltage jumps required",
            ));
        }
        let mut result = Vec::with_capacity(self.points.len());
        for point in &self.points {
            poll(cancelled)?;
            let (mut voltage, mut correction) = (0.0, 0.0);
            for (&basis, &value) in point.basis.iter().zip(&jump) {
                add(&mut voltage, &mut correction, finite(basis * value)?)?;
            }
            let density = finite(voltage / point.resistance_ohm_m2)?;
            let current = finite(point.area_m2 * density)?;
            let power = finite(voltage * current)?;
            if point.area_m2 > 0.0 && voltage != 0.0 && power == 0.0 {
                return Err(QuadratureError::Invalid(
                    "nonzero contact loss is not representable",
                ));
            }
            result.push(ConductionContactQuadraturePoint {
                basis: point.basis,
                area_weight_m2: point.area_m2,
                voltage_jump_v: voltage,
                current_density_a_m2: density,
                resistance_ohm_m2: point.resistance_ohm_m2,
                power_w: power,
            });
        }
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn near(a: f64, b: f64) {
        assert!(
            (a - b).abs() <= 2e-13 * a.abs().max(b.abs()).max(1.0),
            "{a} != {b}"
        );
    }

    #[test]
    fn edge_constant_law_retains_exact_planar_and_radial_moments() {
        for n in [1, 2, 7] {
            for factors in [[1.0, 1.0], [0.0, 2.0], [1.5, 0.5]] {
                let q = ContactQuadrature::new(
                    edge_rule(n, &mut || false).unwrap(),
                    2.0,
                    factors,
                    vec![4.0; edge_point_count(n).unwrap()],
                    &mut || false,
                )
                .unwrap();
                near(q.matrix_entry(0, 0), (1.0 + factors[0]) / 12.0);
                near(q.matrix_entry(1, 1), (1.0 + factors[1]) / 12.0);
                near(q.matrix_entry(0, 1), 1.0 / 12.0);
                let samples = q.response([1.0, 3.0], &mut || false).unwrap();
                near(samples.iter().map(|p| p.area_weight_m2).sum(), 2.0);
                near(
                    samples.iter().map(|p| p.power_w).sum(),
                    (factors[0] + factors[1] * 9.0 + 16.0) / 12.0,
                );
            }
        }
    }

    #[test]
    fn triangle_rule_has_positive_weights_and_exact_cubic_moments() {
        for n in [1, 2, 3, 5] {
            let rule = triangle_rule(n, &mut || false).unwrap();
            assert_eq!(rule.len(), triangle_point_count(n).unwrap());
            assert!(rule.iter().all(|p| p.weight > 0.0));
            near(rule.iter().map(|p| p.weight).sum(), 1.0);
            for axis in 0..3 {
                near(
                    rule.iter().map(|p| p.weight * p.basis[axis]).sum(),
                    1.0 / 3.0,
                );
                near(
                    rule.iter().map(|p| p.weight * p.basis[axis].powi(2)).sum(),
                    1.0 / 6.0,
                );
                near(
                    rule.iter().map(|p| p.weight * p.basis[axis].powi(3)).sum(),
                    1.0 / 10.0,
                );
            }
            near(
                rule.iter()
                    .map(|p| p.weight * p.basis.iter().product::<f64>())
                    .sum(),
                1.0 / 60.0,
            );
            let resistance = rule.iter().map(|p| 1.0 / (1.0 + p.basis[0])).collect();
            let q = ContactQuadrature::new(rule, 1.0, [1.0; 3], resistance, &mut || false).unwrap();
            near(q.matrix_entry(0, 0), 4.0 / 15.0);
            near(q.matrix_entry(0, 1), 7.0 / 60.0);
            near(q.matrix_entry(1, 1), 1.0 / 5.0);
        }
    }

    #[test]
    fn variable_contact_heat_refines_without_changing_operator_work_identity() {
        let mut error = f64::INFINITY;
        for n in [1, 2, 4, 8] {
            let rule = edge_rule(n, &mut || false).unwrap();
            let resistance = rule
                .iter()
                .map(|p| 1.0 / (1.0 + p.basis[1].powi(2)))
                .collect();
            let q = ContactQuadrature::new(rule, 1.0, [1.0; 2], resistance, &mut || false).unwrap();
            let values = [1.0, 3.0];
            let samples = q.response(values, &mut || false).unwrap();
            let power: f64 = samples.iter().map(|p| p.power_w).sum();
            let current: f64 = samples
                .iter()
                .map(|p| p.area_weight_m2 * p.current_density_a_m2)
                .sum();
            let mut operator_power = 0.0;
            let mut operator_current = 0.0;
            for i in 0..2 {
                for j in 0..2 {
                    operator_power += values[i] * q.matrix_entry(i, j) * values[j];
                    operator_current += q.matrix_entry(i, j) * values[j];
                }
            }
            near(power, operator_power);
            near(current, operator_current);
            near(current, 17.0 / 6.0);
            let next = (power - 97.0 / 15.0).abs();
            assert!(next < error / 10.0);
            error = next;
        }
    }

    #[test]
    fn edge_reversal_preserves_current_and_heat_with_nonuniform_sources() {
        let rule = edge_rule(3, &mut || false).unwrap();
        let resistance: Vec<_> = rule.iter().map(|p| 1.0 + p.basis[1]).collect();
        let a = ContactQuadrature::new(
            rule.clone(),
            2.0,
            [0.0, 2.0],
            resistance.clone(),
            &mut || false,
        )
        .unwrap();
        let b = ContactQuadrature::new(
            rule,
            2.0,
            [2.0, 0.0],
            resistance.into_iter().rev().collect(),
            &mut || false,
        )
        .unwrap();
        let a = a.response([-1.0, 3.0], &mut || false).unwrap();
        let b = b.response([3.0, -1.0], &mut || false).unwrap();
        near(
            a.iter().map(|p| p.power_w).sum(),
            b.iter().map(|p| p.power_w).sum(),
        );
        near(
            a.iter()
                .map(|p| p.area_weight_m2 * p.current_density_a_m2)
                .sum(),
            b.iter()
                .map(|p| p.area_weight_m2 * p.current_density_a_m2)
                .sum(),
        );
    }

    #[test]
    fn invalid_extents_coefficients_and_cancelled_work_refuse() {
        assert!(edge_point_count(0).is_err());
        assert!(edge_point_count(usize::MAX).is_err());
        assert!(triangle_point_count(0).is_err());
        assert!(triangle_point_count(usize::MAX).is_err());
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(
                ContactQuadrature::new(
                    edge_rule(1, &mut || false).unwrap(),
                    1.0,
                    [1.0; 2],
                    vec![bad; 3],
                    &mut || false
                )
                .is_err()
            );
        }
        let mut polls = 0;
        assert!(matches!(
            triangle_rule(4, &mut || {
                polls += 1;
                polls == 17
            }),
            Err(QuadratureError::Cancelled)
        ));
        let q = ContactQuadrature::new(
            edge_rule(2, &mut || false).unwrap(),
            1.0,
            [1.0; 2],
            vec![2.0; 5],
            &mut || false,
        )
        .unwrap();
        let original = q.response([1.0, 2.0], &mut || false).unwrap();
        assert!(matches!(
            q.response([1.0, 2.0], &mut || true),
            Err(QuadratureError::Cancelled)
        ));
        assert_eq!(q.response([1.0, 2.0], &mut || false).unwrap(), original);
    }
}
