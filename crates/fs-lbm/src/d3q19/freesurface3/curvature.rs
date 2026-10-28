//! Curvature from complete, consistently oriented VOF height columns.
//!
//! For a graph q=h(u,v), with liquid on q<h, outward curvature is
//! -[h_uu(1+h_v²)+h_vv(1+h_u²)-2h_uv h_u h_v]/(1+h_u²+h_v²)^(3/2).
//! Integrating unsmoothed fractions preserves the column volume. A complete
//! 3×3×7 stencil is required; the caller retains its wall/unresolved fallback.
//! Reference: Popinet, JCP 228 (2009), 5838–5866,
//! <https://doi.org/10.1016/j.jcp.2009.04.042> (fixed complete-height rung only).

/// Sum of principal curvatures in inverse cell widths, positive for a drop.
/// `normal` points out of the liquid. No allocation or field mutation occurs.
pub(super) fn height_function(
    normal: [f64; 3],
    mut sample: impl FnMut([i64; 3]) -> Option<f64>,
) -> Option<f64> {
    if normal.iter().any(|n| !n.is_finite()) {
        return None;
    }
    let mut axis = 0;
    for d in 1..3 {
        if normal[d].abs() > normal[axis].abs() {
            axis = d;
        }
    }
    if normal[axis] == 0.0 {
        return None;
    }
    let outward = if normal[axis] > 0.0 { 1 } else { -1 };
    let tangent = match axis {
        0 => [1, 2],
        1 => [0, 2],
        _ => [0, 1],
    };
    let mut heights = [[0.0; 3]; 3];
    const TOLERANCE: f64 = 64.0 * f64::EPSILON;
    for (a, column) in heights.iter_mut().enumerate() {
        for (b, height) in column.iter_mut().enumerate() {
            let mut previous = 1.0;
            for depth in -3..=3 {
                let mut offset = [0; 3];
                offset[axis] = depth * outward;
                offset[tangent[0]] = a as i64 - 1;
                offset[tangent[1]] = b as i64 - 1;
                let fill = sample(offset)?;
                if !fill.is_finite()
                    || !(-TOLERANCE..=1.0 + TOLERANCE).contains(&fill)
                    || fill > previous + TOLERANCE
                    || (depth == -3 && fill < 1.0 - TOLERANCE)
                    || (depth == 3 && fill > TOLERANCE)
                {
                    return None;
                }
                let fill = fill.clamp(0.0, 1.0);
                *height += fill;
                previous = fill;
            }
        }
    }
    let h = heights;
    let hu = (h[2][1] - h[0][1]) * 0.5;
    let hv = (h[1][2] - h[1][0]) * 0.5;
    let huu = h[2][1] + h[0][1] - 2.0 * h[1][1];
    let hvv = h[1][2] + h[1][0] - 2.0 * h[1][1];
    let huv = (h[2][2] + h[0][0] - h[0][2] - h[2][0]) * 0.25;
    let metric = 1.0 + hu * hu + hv * hv;
    let curvature = -(huu * (1.0 + hv * hv) + hvv * (1.0 + hu * hu) - 2.0 * huv * hu * hv)
        / (metric * metric.sqrt());
    curvature.is_finite().then_some(curvature)
}

#[cfg(test)]
mod tests {
    use super::*;

    // Tensor Gauss quadrature independently integrates cell volumes below the
    // prescribed graph. Summing a complete height column telescopes to the
    // integral of the smooth height, including cells cut across several layers.
    fn graph_fill(offset: [i64; 3], axis: usize, height: impl Fn(f64, f64) -> f64) -> f64 {
        let tangent = match axis {
            0 => [1, 2],
            1 => [0, 2],
            _ => [0, 1],
        };
        let nodes = [
            -0.861_136_311_594_052_6,
            -0.339_981_043_584_856_3,
            0.339_981_043_584_856_3,
            0.861_136_311_594_052_6,
        ];
        let weights = [
            0.347_854_845_137_453_8,
            0.652_145_154_862_546_1,
            0.652_145_154_862_546_1,
            0.347_854_845_137_453_8,
        ];
        let mut volume = 0.0;
        for i in 0..4 {
            for j in 0..4 {
                let u = offset[tangent[0]] as f64 + nodes[i] * 0.5;
                let v = offset[tangent[1]] as f64 + nodes[j] * 0.5;
                volume += weights[i]
                    * weights[j]
                    * 0.25
                    * (height(u, v) - offset[axis] as f64 + 0.5).clamp(0.0, 1.0);
            }
        }
        volume
    }

    #[test]
    fn g1_inclined_plane_has_no_capillary_pressure_in_every_axis() {
        for axis in 0..3 {
            for (a, b) in [(0.0, 0.0), (0.25, -0.4), (-0.5, 0.3)] {
                let mut normal = [0.0; 3];
                normal[axis] = 1.0;
                let curvature = height_function(normal, |p| {
                    Some(graph_fill(p, axis, |u, v| 0.2 + a * u + b * v))
                })
                .unwrap();
                assert!(curvature.abs() < 2e-14, "{axis}: {curvature}");
            }
        }
    }

    #[test]
    fn g1_quadratic_graph_retains_both_principal_and_mixed_curvatures() {
        let (a, b, uu, uv, vv) = (0.2_f64, -0.3_f64, -0.1_f64, 0.04_f64, -0.16_f64);
        let graph = |u: f64, v: f64| {
            0.15 + a * u + b * v + 0.5 * uu * u * u + uv * u * v + 0.5 * vv * v * v
        };
        let metric = 1.0 + a * a + b * b;
        let exact = -(uu * (1.0 + b * b) + vv * (1.0 + a * a) - 2.0 * uv * a * b)
            / (metric * metric.sqrt());
        for reverse in [false, true] {
            let sign = if reverse { -1.0 } else { 1.0 };
            let curvature = height_function([-sign * a, -sign * b, sign], |p| {
                let f = graph_fill(p, 2, graph);
                Some(if reverse { 1.0 - f } else { f })
            })
            .unwrap();
            assert!(
                (curvature - sign * exact).abs() < 2e-14,
                "{curvature} vs {exact}"
            );
        }
    }

    #[test]
    fn g1_sphere_laplace_curvature_converges_under_refinement() {
        for (u_fraction, v_fraction) in [(0.0_f64, 0.0_f64), (0.25_f64, 0.2_f64)] {
            let mut previous = f64::INFINITY;
            for radius in [6.0_f64, 12.0_f64, 24.0_f64] {
                let u0 = radius * u_fraction;
                let v0 = radius * v_fraction;
                let z0 = (radius * radius - u0 * u0 - v0 * v0).sqrt();
                let graph = |u: f64, v: f64| {
                    (radius * radius - (u + u0).powi(2) - (v + v0).powi(2)).sqrt() - z0 + 0.2
                };
                let kappa =
                    height_function([u0, v0, z0], |p| Some(graph_fill(p, 2, graph))).unwrap();
                let relative = (kappa / (2.0 / radius) - 1.0).abs();
                assert!(relative < 0.03, "R={radius}: {relative}");
                assert!(
                    relative < previous * 0.4,
                    "R={radius}: {relative} >= {previous}"
                );
                previous = relative;
            }
        }
    }

    #[test]
    fn g0_incomplete_or_multiple_interfaces_do_not_mint_a_height() {
        let normal = [0.0, 0.0, 1.0];
        let flat = |p: [i64; 3]| (0.7 - p[2] as f64).clamp(0.0, 1.0);
        assert!(
            height_function(normal, |p| {
                if p == [1, 1, 3] { None } else { Some(flat(p)) }
            })
            .is_none()
        );
        assert!(height_function(normal, |_| Some(1.0)).is_none());
        assert!(
            height_function(normal, |p| { Some(if p[2] == -1 { 0.2 } else { flat(p) }) }).is_none()
        );
        assert!(height_function([0.0; 3], |p| Some(flat(p))).is_none());
    }
}
