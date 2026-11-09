//! Bounded integration of the original positive conductivity on each source segment.
//! This is an embedded quadrature ESTIMATE, never an interval certificate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasConductivityQuadrature {
    pub absolute_tolerance_w_m_k: f64,
    pub relative_tolerance: f64,
    /// Total bisections across all crossed source segments in one mean query.
    pub maximum_subdivisions: usize,
    /// At most 32; fixed stack storage bounds every individual segment.
    pub maximum_depth: u8,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasConductivityMean {
    pub mean_w_m_k: f64,
    /// Sum of weighted embedded estimates, including a floating-point floor.
    /// This is neither a rigorous bound nor a physical heat/energy defect.
    pub estimated_absolute_error_w_m_k: f64,
    pub point_evaluations: usize,
    pub subdivisions: usize,
    pub source_segments: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub enum Error<E> {
    Invalid,
    Arithmetic,
    WorkLimit,
    Resolution,
    Integrand(E),
}

impl GasConductivityQuadrature {
    pub fn validate(self) -> bool {
        self.absolute_tolerance_w_m_k.is_finite()
            && self.absolute_tolerance_w_m_k >= 0.0
            && self.relative_tolerance.is_finite()
            && self.relative_tolerance >= 0.0
            && self.relative_tolerance < 1.0
            && (self.absolute_tolerance_w_m_k > 0.0 || self.relative_tolerance > 0.0)
            && self.maximum_depth <= 32
    }
}

#[derive(Clone, Copy, Default)]
struct Panel {
    left: f64,
    right: f64,
    endpoints: [f64; 2],
    depth: u8,
}

#[derive(Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) {
        let y = value - self.correction;
        let next = self.value + y;
        self.correction = (next - self.value) - y;
        self.value = next;
    }
}

fn sample<E>(
    temperature: f64,
    count: &mut usize,
    f: &mut impl FnMut(f64) -> Result<f64, E>,
) -> Result<f64, Error<E>> {
    *count = count.checked_add(1).ok_or(Error::Arithmetic)?;
    let value = f(temperature).map_err(Error::Integrand)?;
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(Error::Arithmetic)
    }
}

// QUADPACK dqk15 Gauss-7/Kronrod-15 nodes and weights. Integrating the
// normalized mean avoids large absolute potentials and dimensional overflow.
fn panel<E>(
    p: Panel,
    count: &mut usize,
    f: &mut impl FnMut(f64) -> Result<f64, E>,
) -> Result<[f64; 3], Error<E>> {
    const X: [f64; 8] = [
        0.9914553711208126,
        0.9491079123427585,
        0.8648644233597691,
        0.7415311855993945,
        0.5860872354676911,
        0.4058451513773972,
        0.2077849550078985,
        0.0,
    ];
    const K: [f64; 8] = [
        0.02293532201052922,
        0.06309209262997855,
        0.1047900103222502,
        0.1406532597155259,
        0.1690047266392679,
        0.1903505780647854,
        0.2044329400752989,
        0.2094821410847278,
    ];
    const G: [f64; 4] = [
        0.1294849661688697,
        0.2797053914892767,
        0.3818300505051189,
        0.4179591836734694,
    ];
    let width = p.right - p.left;
    let center = sample(p.left + 0.5 * width, count, f)?;
    let mut values = [[0.0; 2]; 7];
    let mut kronrod = Sum::default();
    let mut gauss = Sum::default();
    kronrod.add((0.5 * K[7]) * center);
    gauss.add((0.5 * G[3]) * center);
    for i in 0..7 {
        for side in 0..2 {
            let fraction = 0.5 * (1.0 + if side == 0 { -X[i] } else { X[i] });
            let value = sample(p.left + fraction * width, count, f)?;
            values[i][side] = value;
            kronrod.add((0.5 * K[i]) * value);
            if i % 2 == 1 {
                gauss.add((0.5 * G[i / 2]) * value);
            }
        }
    }
    let mean = kronrod.value;
    let mut variation = Sum::default();
    variation.add((0.5 * K[7]) * (center - mean).abs());
    for i in 0..7 {
        for value in values[i] {
            variation.add((0.5 * K[i]) * (value - mean).abs());
        }
    }
    let raw_error = (mean - gauss.value).abs();
    let mut error = raw_error;
    if variation.value > 0.0 && raw_error > 0.0 {
        let ratio = (200.0 * (raw_error / variation.value)).min(1.0);
        error = variation.value * ratio * ratio.sqrt();
    }
    error = error.max((50.0 * f64::EPSILON) * mean);
    if !mean.is_finite() || mean <= 0.0 || !error.is_finite() {
        return Err(Error::Arithmetic);
    }
    Ok([mean, error, center])
}

/// The caller has prepared strictly increasing finite source knots. Every
/// panel remains within one original segment; no value fit is constructed.
pub fn mean<E>(
    knots: &[f64],
    from: f64,
    to: f64,
    controls: GasConductivityQuadrature,
    mut f: impl FnMut(f64) -> Result<f64, E>,
) -> Result<GasConductivityMean, Error<E>> {
    if !controls.validate()
        || !from.is_finite()
        || from <= 0.0
        || !to.is_finite()
        || to <= 0.0
        || knots.len() < 2
        || from.min(to) < knots[0]
        || from.max(to) > knots[knots.len() - 1]
    {
        return Err(Error::Invalid);
    }
    let mut result = GasConductivityMean {
        mean_w_m_k: 0.0,
        estimated_absolute_error_w_m_k: 0.0,
        point_evaluations: 0,
        subdivisions: 0,
        source_segments: 0,
    };
    if from == to {
        result.mean_w_m_k = sample(from, &mut result.point_evaluations, &mut f)?;
        return Ok(result);
    }
    let lower = from.min(to);
    let upper = from.max(to);
    let width = upper - lower;
    let first = knots.partition_point(|&t| t <= lower).saturating_sub(1);
    let end = knots.partition_point(|&t| t < upper);
    let mut total = Sum::default();
    let mut errors = Sum::default();
    for segment in first..end {
        let left = lower.max(knots[segment]);
        let right = upper.min(knots[segment + 1]);
        result.source_segments += 1;
        let endpoints = [
            sample(left, &mut result.point_evaluations, &mut f)?,
            sample(right, &mut result.point_evaluations, &mut f)?,
        ];
        let mut stack = [Panel::default(); 34];
        stack[0] = Panel {
            left,
            right,
            endpoints,
            depth: 0,
        };
        let mut pending = 1;
        while pending != 0 {
            pending -= 1;
            let p = stack[pending];
            let [value, error, center] = panel(p, &mut result.point_evaluations, &mut f)?;
            // An endpoint layer missed by both interior rules is unresolved.
            // This guard supplements, rather than certifies, the error estimate.
            let endpoint_resolved = value >= 1e-5 * p.endpoints[0].max(p.endpoints[1]);
            let tolerance = controls.absolute_tolerance_w_m_k + controls.relative_tolerance * value;
            if !tolerance.is_finite() {
                return Err(Error::Arithmetic);
            }
            if endpoint_resolved && error <= tolerance {
                let weight = (p.right - p.left) / width;
                let contribution = weight * value;
                if weight <= 0.0 || !contribution.is_finite() || contribution <= 0.0 {
                    return Err(Error::Arithmetic);
                }
                total.add(contribution);
                errors.add(weight * error);
            } else {
                if result.subdivisions >= controls.maximum_subdivisions {
                    return Err(Error::WorkLimit);
                }
                let middle = p.left + 0.5 * (p.right - p.left);
                if p.depth >= controls.maximum_depth
                    || middle <= p.left
                    || middle >= p.right
                    || pending + 2 > stack.len()
                {
                    return Err(Error::Resolution);
                }
                result.subdivisions += 1;
                // Left before right gives a deterministic physical-temperature order.
                stack[pending] = Panel {
                    left: middle,
                    right: p.right,
                    endpoints: [center, p.endpoints[1]],
                    depth: p.depth + 1,
                };
                stack[pending + 1] = Panel {
                    left: p.left,
                    right: middle,
                    endpoints: [p.endpoints[0], center],
                    depth: p.depth + 1,
                };
                pending += 2;
            }
        }
    }
    result.mean_w_m_k = total.value;
    result.estimated_absolute_error_w_m_k = errors.value;
    let tolerance = controls.absolute_tolerance_w_m_k + controls.relative_tolerance * total.value;
    if !total.value.is_finite()
        || total.value <= 0.0
        || !errors.value.is_finite()
        || !tolerance.is_finite()
    {
        return Err(Error::Arithmetic);
    }
    if errors.value > tolerance {
        return Err(Error::Resolution);
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn controls() -> GasConductivityQuadrature {
        GasConductivityQuadrature {
            absolute_tolerance_w_m_k: 1e-12,
            relative_tolerance: 1e-11,
            maximum_subdivisions: 2048,
            maximum_depth: 32,
        }
    }
    #[test]
    fn original_breaks_and_reverse_mean_are_preserved() {
        let f = |t: f64| {
            Ok::<_, ()>(if t < 2.0 {
                2.0 + t
            } else {
                4.0 + 0.5 * (t - 2.0).powi(2)
            })
        };
        let a = mean(&[1.0, 2.0, 5.0], 1.0, 5.0, controls(), f).unwrap();
        let b = mean(&[1.0, 2.0, 5.0], 5.0, 1.0, controls(), f).unwrap();
        assert_eq!(a, b);
        assert!((a.mean_w_m_k - 5.0).abs() < 2e-14);
        assert_eq!(a.source_segments, 2);
        assert_eq!(a.point_evaluations, 34);
        assert!((a.mean_w_m_k - f(3.0).unwrap()).abs() > 0.4);
        let point = mean(&[1.0, 2.0, 5.0], 2.0, 2.0, controls(), f).unwrap();
        assert_eq!(point.mean_w_m_k, 4.0);
        assert_eq!(point.point_evaluations, 1);
        assert_eq!(point.estimated_absolute_error_w_m_k, 0.0);
    }
    #[test]
    fn resolves_sharp_positive_endpoint_layer_against_log_integral() {
        let epsilon = 1e-6;
        let value = mean(&[1.0, 2.0], 1.0, 2.0, controls(), |t| {
            Ok::<_, ()>(1.0 / ((t - 1.0) + epsilon))
        })
        .unwrap();
        let exact = ((1.0 + epsilon) / epsilon).ln();
        assert!(
            (value.mean_w_m_k - exact).abs() < 2e-9,
            "{value:?} != {exact}"
        );
        assert!(value.subdivisions > 0);
        assert!(value.estimated_absolute_error_w_m_k <= 1e-12 + 1e-11 * value.mean_w_m_k);
    }
    #[test]
    fn normalized_means_avoid_overflowing_dimensional_integrals() {
        let value = mean(
            &[1e100, 3e100],
            1e100,
            3e100,
            GasConductivityQuadrature {
                absolute_tolerance_w_m_k: 0.0,
                ..controls()
            },
            |_| Ok::<_, ()>(1e300),
        )
        .unwrap();
        assert!((value.mean_w_m_k / 1e300 - 1.0).abs() < 1e-14);
    }
    #[test]
    fn work_depth_invalid_controls_and_nonpositive_values_refuse() {
        let f = |t: f64| Ok::<_, ()>((10.0 * t).exp());
        let limited = GasConductivityQuadrature {
            maximum_subdivisions: 0,
            ..controls()
        };
        assert_eq!(
            mean(&[1.0, 2.0], 1.0, 2.0, limited, f),
            Err(Error::WorkLimit)
        );
        let shallow = GasConductivityQuadrature {
            maximum_depth: 0,
            ..controls()
        };
        assert_eq!(
            mean(&[1.0, 2.0], 1.0, 2.0, shallow, f),
            Err(Error::Resolution)
        );
        for bad in [f64::NAN, -1.0, f64::INFINITY] {
            assert_eq!(
                mean(
                    &[1.0, 2.0],
                    1.0,
                    2.0,
                    GasConductivityQuadrature {
                        absolute_tolerance_w_m_k: bad,
                        ..controls()
                    },
                    f
                ),
                Err(Error::Invalid)
            );
        }
        assert_eq!(
            mean(&[1.0, 2.0], 1.0, 2.0, controls(), |_| Ok::<_, ()>(0.0)),
            Err(Error::Arithmetic)
        );
    }
    #[test]
    fn caller_checkpoint_error_returns_without_partial_mean() {
        let mut calls = 0;
        let result = mean(&[1.0, 2.0], 1.0, 2.0, controls(), |_| {
            calls += 1;
            if calls == 4 {
                Err("cancelled")
            } else {
                Ok(1.0)
            }
        });
        assert_eq!(result, Err(Error::Integrand("cancelled")));
        assert_eq!(calls, 4);
    }
}
