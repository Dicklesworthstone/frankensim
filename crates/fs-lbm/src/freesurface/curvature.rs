//! Resolved-interface height-function curvature in lattice-cell units.
//!
//! Three transverse columns integrate seven unsmoothed VOF fractions each.
//! Pure phases at both ends and one monotone transition are required; the
//! caller retains its wall/contact or unresolved-interface fallback otherwise.
//! The signed liquid-column height gives -h''/(1+h'^2)^(3/2), positive for a
//! convex liquid disk. This is a spatial geometry approximation, not a
//! capillary-wave, wetting-calibration or surface-energy certificate.
//!
//! Reference: Popinet, J. Comput. Phys. 228 (2009), 5838–5866,
//! <https://doi.org/10.1016/j.jcp.2009.04.042>.

/// Estimate the outward liquid curvature from a complete local height stencil.
/// The supplied fill gradient only selects the dominant integration axis.
pub(super) fn height_curvature2(
    gradient: [f64; 2],
    mut fill: impl FnMut(i32, i32) -> Option<f64>,
) -> Option<f64> {
    if gradient.iter().any(|v| !v.is_finite()) {
        return None;
    }
    let axis = usize::from(gradient[1].abs() > gradient[0].abs());
    if gradient[axis].abs() <= 1e-12 {
        return None;
    }
    let mut heights = [0.0; 3];
    let mut common_orientation = 0.0;
    const TOL: f64 = 64.0 * f64::EPSILON;
    for (column, height) in heights.iter_mut().enumerate() {
        let transverse = column as i32 - 1;
        let mut values = [0.0; 7];
        for (depth, value) in values.iter_mut().enumerate() {
            let along = depth as i32 - 3;
            let (dx, dy) = if axis == 0 {
                (along, transverse)
            } else {
                (transverse, along)
            };
            *value = fill(dx, dy)?;
            if !value.is_finite() || !(0.0..=1.0).contains(value) {
                return None;
            }
        }
        let orientation = if values[0] >= 1.0 - TOL && values[6] <= TOL {
            1.0
        } else if values[0] <= TOL && values[6] >= 1.0 - TOL {
            -1.0
        } else {
            return None;
        };
        if column > 0 && orientation != common_orientation {
            return None;
        }
        common_orientation = orientation;
        if values
            .windows(2)
            .any(|pair| orientation * (pair[1] - pair[0]) > TOL)
        {
            return None;
        }
        *height = values.iter().sum();
    }
    let slope = 0.5 * (heights[2] - heights[0]);
    let second = heights[2] - 2.0 * heights[1] + heights[0];
    let metric = 1.0 + slope * slope;
    let curvature = -second / (metric * metric.sqrt());
    curvature.is_finite().then_some(curvature)
}

#[cfg(test)]
mod tests {
    use super::height_curvature2;

    fn circular_cap_fill(radius: f64, x: i32, y: i32) -> f64 {
        // Exact cell integral of the upper circular graph. For these resolved
        // pole stencils its entire graph lies inside the y=0 cell; other rows
        // are pure. This avoids using the estimator as its own geometry oracle.
        let primitive = |t: f64| {
            0.5 * (t * (radius * radius - t * t).sqrt() + radius * radius * (t / radius).asin())
        };
        let center = f64::from(x);
        let height = primitive(center + 0.5) - primitive(center - 0.5) - radius;
        (height - (f64::from(y) - 0.5)).clamp(0.0, 1.0)
    }

    #[test]
    fn planar_interfaces_have_zero_curvature_in_both_axis_orders() {
        for slope in [-0.6, 0.0, 0.35] {
            for axis in 0..2 {
                let gradient = if axis == 0 {
                    [-1.0, slope]
                } else {
                    [slope, -1.0]
                };
                let kappa = height_curvature2(gradient, |x, y| {
                    let (transverse, depth) = if axis == 0 { (y, x) } else { (x, y) };
                    Some(
                        (0.17 + slope * f64::from(transverse) - (f64::from(depth) - 0.5))
                            .clamp(0.0, 1.0),
                    )
                })
                .expect("resolved planar stencil");
                assert!(
                    kappa.abs() < 2e-14,
                    "axis={axis}, slope={slope}, kappa={kappa}"
                );
            }
        }
    }

    #[test]
    fn circular_laplace_curvature_converges_at_second_order() {
        let mut previous_error = f64::INFINITY;
        for radius in [6.0, 12.0, 24.0] {
            let kappa =
                height_curvature2([0.0, -1.0], |x, y| Some(circular_cap_fill(radius, x, y)))
                    .expect("resolved circular cap");
            let relative_error = (radius * kappa - 1.0).abs();
            assert!(kappa > 0.0);
            assert!(relative_error < 0.011);
            assert!(relative_error < 0.27 * previous_error);
            previous_error = relative_error;
        }
    }

    #[test]
    fn curvature_rotates_and_reverses_between_drop_and_bubble() {
        let radius = 12.0;
        let drop =
            height_curvature2([0.0, -1.0], |x, y| Some(circular_cap_fill(radius, x, y))).unwrap();
        let rotated =
            height_curvature2([-1.0, 0.0], |x, y| Some(circular_cap_fill(radius, y, x))).unwrap();
        let bubble = height_curvature2([0.0, 1.0], |x, y| {
            Some(1.0 - circular_cap_fill(radius, x, y))
        })
        .unwrap();
        assert_eq!(drop, rotated);
        assert!((drop + bubble).abs() < 2e-14);
    }

    #[test]
    fn incomplete_wall_and_folded_columns_request_fallback() {
        assert_eq!(height_curvature2([0.0, -1.0], |_, _| Some(0.5)), None);
        assert_eq!(height_curvature2([0.0, 0.0], |_, _| Some(0.0)), None);
        assert_eq!(
            height_curvature2([0.0, -1.0], |x, y| {
                if x == 1 && y == 0 {
                    None
                } else {
                    Some((0.5 - f64::from(y)).clamp(0.0, 1.0))
                }
            }),
            None
        );
        assert_eq!(
            height_curvature2([0.0, -1.0], |_, y| {
                Some([1.0, 1.0, 0.3, 0.7, 0.0, 0.0, 0.0][(y + 3) as usize])
            }),
            None
        );
        assert_eq!(height_curvature2([0.0, -1.0], |_, _| Some(f64::NAN)), None);
    }
}
