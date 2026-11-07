//! Piecewise-linear prescribed viscosity, with local (not cumulative) integrals.

/// One authored time/kinematic-viscosity sample. Interpolation is explicitly
/// piecewise linear in physical time; no constitutive law is inferred between
/// samples. Discontinuous changes require separate histories on either side.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ViscosityKnot {
    /// Supported physical time [s].
    pub time_s: f64,
    /// Spatially uniform kinematic viscosity [m²/s].
    pub kinematic_viscosity_m2_s: f64,
}

/// Invalid history support or a time/coefficient that cannot be represented.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ViscosityHistoryError {
    /// Invalid authored time, coefficient or step bound.
    Invalid(&'static str),
    /// Requested window leaves the declared interpolation support.
    OutsideSupport,
    /// No later represented time fits the strict step bound.
    TimeStagnation,
    /// A positive mean/exposure underflows, or an exposure overflows.
    Unrepresentable,
}
impl core::fmt::Display for ViscosityHistoryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "viscosity history: {self:?}")
    }
}
impl std::error::Error for ViscosityHistoryError {}

/// A borrowed, admitted homogeneous viscosity history. No extrapolation,
/// temperature solve, spatial viscosity gradient or stochastic path is implied.
#[derive(Debug, Clone, Copy)]
pub struct ViscosityHistory<'a> {
    knots: &'a [ViscosityKnot],
}

/// One physical interval wholly inside a single linear coefficient segment.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ViscosityInterval {
    /// Exact represented start and end times [s].
    pub time_s: [f64; 2],
    /// Exact interval means of the declared linear law on its first/second
    /// physical half. Keeping these separate captures the order of diffusion
    /// and advection when viscosity evolves during a mechanical step.
    pub half_mean_viscosities_m2_s: [f64; 2],
}

impl<'a> ViscosityHistory<'a> {
    /// Admit at least two finite ordered time samples with finite nonnegative
    /// viscosity. Zero is an explicitly authored inviscid value.
    pub fn piecewise_linear(knots: &'a [ViscosityKnot]) -> Result<Self, ViscosityHistoryError> {
        if knots.len() < 2
            || knots.iter().any(|k| {
                !k.time_s.is_finite()
                    || !k.kinematic_viscosity_m2_s.is_finite()
                    || k.kinematic_viscosity_m2_s < 0.0
            })
            || knots.windows(2).any(|k| {
                let dt = k[1].time_s - k[0].time_s;
                !dt.is_finite() || dt <= 0.0
            })
        {
            return Err(ViscosityHistoryError::Invalid(
                "two or more ordered finite times and nonnegative finite viscosities required",
            ));
        }
        Ok(Self { knots })
    }

    pub(super) const fn from_admitted(knots: &'a [ViscosityKnot]) -> Self {
        Self { knots }
    }

    /// Original authored samples, unchanged by admission.
    #[must_use]
    pub fn knots(self) -> &'a [ViscosityKnot] {
        self.knots
    }

    /// Inclusive first/last supported times [s].
    #[must_use]
    pub fn support_s(self) -> [f64; 2] {
        [
            self.knots[0].time_s,
            self.knots[self.knots.len() - 1].time_s,
        ]
    }

    fn admit_window(self, start: f64, end: f64) -> Result<(), ViscosityHistoryError> {
        let [lower, upper] = self.support_s();
        if !start.is_finite() || !end.is_finite() || end < start || !(end - start).is_finite() {
            return Err(ViscosityHistoryError::Invalid(
                "finite forward window required",
            ));
        }
        if start < lower || end > upper {
            return Err(ViscosityHistoryError::OutsideSupport);
        }
        Ok(())
    }

    /// Original interpolation at one supported time, including both endpoints.
    pub fn viscosity_at(self, time_s: f64) -> Result<f64, ViscosityHistoryError> {
        self.admit_window(time_s, time_s)?;
        let next = self.knots.partition_point(|k| k.time_s <= time_s);
        if next == self.knots.len() {
            return Ok(self.knots[next - 1].kinematic_viscosity_m2_s);
        }
        let left = self.knots[next - 1];
        let right = self.knots[next];
        let a = left.kinematic_viscosity_m2_s;
        let b = right.kinematic_viscosity_m2_s;
        if time_s == left.time_s || a == b {
            return Ok(a);
        }
        // Form the product before a tiny time fraction can underflow. The
        // result lies in the coefficient hull; min only repairs roundoff.
        let distance = if a < b {
            time_s - left.time_s
        } else {
            right.time_s - time_s
        };
        let value = (a.min(b) + product_ratio((b - a).abs(), distance, right.time_s - left.time_s))
            .min(a.max(b));
        if value == 0.0 && a.max(b) > 0.0 {
            return Err(ViscosityHistoryError::Unrepresentable);
        }
        Ok(value)
    }

    /// Next bounded physical step, cut exactly at the next coefficient knot.
    /// Quarter-point coefficient interpolation gives exact half-interval means
    /// without constructing possibly unrepresentable absolute midpoint times.
    pub fn next_interval(
        self,
        start_s: f64,
        end_s: f64,
        maximum_step_s: f64,
    ) -> Result<Option<ViscosityInterval>, ViscosityHistoryError> {
        self.admit_window(start_s, end_s)?;
        if !maximum_step_s.is_finite() || maximum_step_s <= 0.0 {
            return Err(ViscosityHistoryError::Invalid(
                "positive finite maximum step required",
            ));
        }
        if start_s == end_s {
            return Ok(None);
        }
        let next = self.knots.partition_point(|k| k.time_s <= start_s);
        let limit = end_s.min(self.knots[next].time_s);
        let mut end = if maximum_step_s >= limit - start_s {
            limit
        } else {
            start_s + maximum_step_s
        };
        if end - start_s > maximum_step_s {
            end = end.next_down();
        }
        if end <= start_s || end > limit {
            return Err(ViscosityHistoryError::TimeStagnation);
        }
        let a = self.viscosity_at(start_s)?;
        let b = self.viscosity_at(end)?;
        let halves = [mix(a, b, 0.25), mix(a, b, 0.75)];
        if a.max(b) > 0.0 && halves.contains(&0.0) {
            return Err(ViscosityHistoryError::Unrepresentable);
        }
        Ok(Some(ViscosityInterval {
            time_s: [start_s, end],
            half_mean_viscosities_m2_s: halves,
        }))
    }

    /// Mean viscosity over a supported window, integrating each linear segment
    /// analytically. Local intervals are accumulated without subtracting large
    /// cumulative exposures. A zero-length window returns its endpoint value.
    pub fn mean_viscosity(self, start_s: f64, end_s: f64) -> Result<f64, ViscosityHistoryError> {
        self.admit_window(start_s, end_s)?;
        if start_s == end_s {
            return self.viscosity_at(start_s);
        }
        let mut left = start_s;
        let mut mean = 0.0;
        while left < end_s {
            let next = self.knots.partition_point(|k| k.time_s <= left);
            let right = end_s.min(self.knots[next].time_s);
            let a = self.viscosity_at(left)?;
            let b = self.viscosity_at(right)?;
            let local = mix(a, b, 0.5);
            if local == 0.0 && a.max(b) > 0.0 {
                return Err(ViscosityHistoryError::Unrepresentable);
            }
            if left == start_s {
                mean = local;
            } else if mean != local {
                let span = right - start_s;
                let upper = mean.max(local);
                let lower = mean.min(local);
                let weight = if mean < local {
                    right - left
                } else {
                    left - start_s
                };
                mean = (lower + product_ratio(upper - lower, weight, span)).min(upper);
                if mean == 0.0 && upper > 0.0 {
                    return Err(ViscosityHistoryError::Unrepresentable);
                }
            }
            if !mean.is_finite() {
                return Err(ViscosityHistoryError::Unrepresentable);
            }
            left = right;
        }
        Ok(mean)
    }

    /// Integrated kinematic viscosity [m²], the diffusion exposure of the
    /// declared homogeneous law. A physically nonzero underflow refuses;
    /// `mean_viscosity` remains available when the dimensional product is too
    /// small or large for f64 but the native kernel can use separate factors.
    pub fn integrated_viscosity_m2(
        self,
        start_s: f64,
        end_s: f64,
    ) -> Result<f64, ViscosityHistoryError> {
        let mean = self.mean_viscosity(start_s, end_s)?;
        let duration = end_s - start_s;
        let exposure = mean * duration;
        if !exposure.is_finite() || (mean > 0.0 && duration > 0.0 && exposure == 0.0) {
            return Err(ViscosityHistoryError::Unrepresentable);
        }
        Ok(exposure)
    }
}

fn mix(a: f64, b: f64, fraction: f64) -> f64 {
    if a == b || fraction == 0.0 {
        a
    } else if fraction == 1.0 {
        b
    } else if a < b {
        fraction.mul_add(b - a, a)
    } else {
        (1.0 - fraction).mul_add(a - b, b)
    }
}

// Positive product/ratio without first overflowing its numerator or losing a
// tiny time fraction which a large coefficient makes physically representable.
fn product_ratio(a: f64, b: f64, divisor: f64) -> f64 {
    let (large, small) = (a.max(b), a.min(b));
    let ratio = large / divisor;
    if ratio.is_finite() {
        ratio * small
    } else {
        large * (small / divisor)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn knot(time_s: f64, nu: f64) -> ViscosityKnot {
        ViscosityKnot {
            time_s,
            kinematic_viscosity_m2_s: nu,
        }
    }

    #[test]
    fn linear_history_has_distinct_exact_strang_half_exposures() {
        let knots = [knot(0.0, 0.2), knot(2.0, 1.0)];
        let history = ViscosityHistory::piecewise_linear(&knots).unwrap();
        let interval = history.next_interval(0.0, 2.0, 3.0).unwrap().unwrap();
        assert_eq!(interval.time_s, [0.0, 2.0]);
        assert!((interval.half_mean_viscosities_m2_s[0] - 0.4).abs() < 1e-15);
        assert!((interval.half_mean_viscosities_m2_s[1] - 0.8).abs() < 1e-15);
        assert!((history.integrated_viscosity_m2(0.0, 2.0).unwrap() - 1.2).abs() < 1e-15);
    }

    #[test]
    fn coefficient_events_are_exact_and_segmented_exposure_adds() {
        let knots = [
            knot(0.0, 0.1),
            knot(0.25, 0.7),
            knot(0.75, 0.0),
            knot(1.0, 0.3),
        ];
        let history = ViscosityHistory::piecewise_linear(&knots).unwrap();
        let mut time = 0.0;
        let mut ends = Vec::new();
        let mut exposure = 0.0;
        while let Some(interval) = history.next_interval(time, 1.0, 0.3125).unwrap() {
            let [a, b] = interval.time_s;
            exposure += (b - a)
                * mix(
                    interval.half_mean_viscosities_m2_s[0],
                    interval.half_mean_viscosities_m2_s[1],
                    0.5,
                );
            ends.push(b);
            time = b;
        }
        assert_eq!(ends, vec![0.25, 0.5625, 0.75, 1.0]);
        assert!((exposure - (0.1 + 0.175 + 0.0375)).abs() < 2e-16);
        assert!((exposure - history.integrated_viscosity_m2(0.0, 1.0).unwrap()).abs() < 2e-16);
    }

    #[test]
    fn local_integrals_retain_tiny_late_exposure_and_large_clock_halves() {
        let tiny = 2.0_f64.powi(-40);
        let knots = [knot(0.0, 1e15), knot(1.0, 1e-9), knot(1.0 + tiny, 1e-9)];
        let history = ViscosityHistory::piecewise_linear(&knots).unwrap();
        assert_eq!(
            history.integrated_viscosity_m2(1.0, 1.0 + tiny).unwrap(),
            tiny * 1e-9
        );
        let start = 1e16;
        let knots = [knot(start, 1.0), knot(start + 2.0, 5.0)];
        let history = ViscosityHistory::piecewise_linear(&knots).unwrap();
        assert_eq!(start + 1.0, start);
        let interval = history
            .next_interval(start, start + 2.0, 2.0)
            .unwrap()
            .unwrap();
        assert_eq!(interval.half_mean_viscosities_m2_s, [2.0, 4.0]);
    }

    #[test]
    fn autonomous_histories_retain_viscosity_bits_and_zero_intervals() {
        for nu in [0.0, f64::from_bits(1), 0.0013, f64::MAX] {
            let knots = [knot(0.0, nu), knot(0.3, nu), knot(1.0, nu)];
            let history = ViscosityHistory::piecewise_linear(&knots).unwrap();
            assert_eq!(
                history.mean_viscosity(0.1, 0.9).unwrap().to_bits(),
                nu.to_bits()
            );
            assert_eq!(
                history
                    .next_interval(0.1, 0.9, 0.2)
                    .unwrap()
                    .unwrap()
                    .half_mean_viscosities_m2_s,
                [nu; 2]
            );
            assert!(history.next_interval(0.3, 0.3, 1.0).unwrap().is_none());
        }
    }

    #[test]
    fn representable_products_survive_tiny_time_fractions_and_step_caps_are_strict() {
        let knots = [knot(0.0, 0.0), knot(1e300, 1e300)];
        let history = ViscosityHistory::piecewise_linear(&knots).unwrap();
        assert_eq!(history.viscosity_at(1e-300).unwrap(), 1e-300);
        let knots = [knot(1e16, 0.2), knot(1e16 + 8.0, 0.2)];
        let history = ViscosityHistory::piecewise_linear(&knots).unwrap();
        let interval = history
            .next_interval(1e16, 1e16 + 8.0, 3.0)
            .unwrap()
            .unwrap();
        assert_eq!(interval.time_s[1] - interval.time_s[0], 2.0);
        assert_eq!(
            history.next_interval(1e16, 1e16 + 8.0, 1.0),
            Err(ViscosityHistoryError::TimeStagnation)
        );
        let knots = [knot(0.0, 0.0), knot(1.0, f64::from_bits(1))];
        let history = ViscosityHistory::piecewise_linear(&knots).unwrap();
        assert_eq!(
            history.next_interval(0.0, 1.0, 1.0),
            Err(ViscosityHistoryError::Unrepresentable)
        );
    }

    #[test]
    fn invalid_and_unsupported_histories_refuse_without_extrapolation() {
        for knots in [
            [knot(0.0, 1.0), knot(0.0, 2.0)],
            [knot(0.0, -1.0), knot(1.0, 2.0)],
        ] {
            assert!(ViscosityHistory::piecewise_linear(&knots).is_err());
        }
        let knots = [knot(2.0, 1.0), knot(3.0, 2.0)];
        let history = ViscosityHistory::piecewise_linear(&knots).unwrap();
        assert_eq!(
            history.mean_viscosity(1.0, 2.5),
            Err(ViscosityHistoryError::OutsideSupport)
        );
        assert_eq!(
            history.viscosity_at(3.5),
            Err(ViscosityHistoryError::OutsideSupport)
        );
        assert!(history.next_interval(2.0, 3.0, 0.0).is_err());
    }
}
