//! Retained reservoir-pressure forcing for every native transient owner.

use super::ChtError;

/// Piecewise-linear scalar modulation of the authored reservoir pressure
/// differences. One to sixteen `(time_s, factor)` points are retained inline;
/// times must be finite, non-negative and strictly increasing, and factors
/// finite. Negative factors reverse the prescribed head, while zero equalizes
/// the reservoirs within each disconnected fluid component.
///
/// A native march samples the program once at each step's **end** time. Between
/// knots the factor is linear; before the first and after the last knot it is
/// held at that endpoint's value. A one-point program is constant. Resolve
/// changes with the time grid: knots do not split steps or reset BE/BDF2 history.
///
/// The factor acts on `p_authored - p_reference`, retaining each component's
/// original datum. Momentum and thermal backflow share the same projected
/// flux. Existing physical-heating admission rules still apply. It does not
/// prescribe absolute thermodynamic pressure, fan curves or source EOS pressure.
/// Independent opening histories use [`PressureOffsetProgram`]. Source-material
/// queries retain their original pressure context.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PressureProgram {
    points: [(f64, f64); 16],
    len: usize,
}

impl PressureProgram {
    /// Maximum number of inline time/factor knots.
    pub const MAX_POINTS: usize = 16;

    /// The unchanged authored differential pressure profile.
    pub const IDENTITY: Self = Self {
        points: [(0.0, 1.0); 16],
        len: 1,
    };

    /// Retain a bounded pressure program without a heap allocation.
    ///
    /// # Errors
    /// [`ChtError::InvalidInput`] for zero or more than sixteen points,
    /// non-finite entries, negative time, or non-increasing knot times.
    pub fn new(points: &[(f64, f64)]) -> Result<Self, ChtError> {
        Self::from_points(points, "unsteady.pressure_program")
    }

    fn from_points(points: &[(f64, f64)], field: &'static str) -> Result<Self, ChtError> {
        let invalid = |reason: &str| ChtError::InvalidInput {
            field,
            reason: reason.into(),
        };
        if points.is_empty() || points.len() > Self::MAX_POINTS {
            return Err(invalid("requires one to sixteen (time_s, value) points"));
        }
        for (i, &(time, factor)) in points.iter().enumerate() {
            if !time.is_finite() || time < 0.0 || !factor.is_finite() {
                return Err(invalid(
                    "times must be finite and non-negative, values finite",
                ));
            }
            if i > 0 && time <= points[i - 1].0 {
                return Err(invalid("knot times must strictly increase"));
            }
        }
        let mut stored = [(0.0, 0.0); 16];
        stored[..points.len()].copy_from_slice(points);
        Ok(Self {
            points: stored,
            len: points.len(),
        })
    }

    /// Original ordered time/factor knots.
    #[must_use]
    pub fn points(&self) -> &[(f64, f64)] {
        &self.points[..self.len]
    }

    /// Pressure factor at a finite, non-negative time, with endpoint holds.
    ///
    /// Equal-factor segments return that exact factor. Interpolation avoids
    /// overflowing the difference between large, opposite-signed factors.
    /// The native solver separately checks that the sampled factor times every
    /// active differential reservoir pressure remains representable.
    ///
    /// # Errors
    /// [`ChtError::InvalidInput`] for non-finite or negative sample time.
    pub fn sample(&self, time_s: f64) -> Result<f64, ChtError> {
        self.sample_at(time_s, "unsteady.pressure_program.time_s")
    }

    fn sample_at(&self, time_s: f64, field: &'static str) -> Result<f64, ChtError> {
        if !time_s.is_finite() || time_s < 0.0 {
            return Err(ChtError::InvalidInput {
                field,
                reason: "sample time must be finite and non-negative".into(),
            });
        }
        let points = self.points();
        if time_s <= points[0].0 {
            return Ok(points[0].1);
        }
        for pair in points.windows(2) {
            let ((t0, f0), (t1, f1)) = (pair[0], pair[1]);
            if time_s > t1 {
                continue;
            }
            if time_s == t1 || f0 == f1 {
                return Ok(f1);
            }
            let weight = (time_s - t0) / (t1 - t0);
            let factor = if f0.is_sign_negative() == f1.is_sign_negative() {
                f0 + weight * (f1 - f0)
            } else {
                (1.0 - weight) * f0 + weight * f1
            };
            return Ok(factor);
        }
        Ok(points[self.len - 1].1)
    }
}

impl Default for PressureProgram {
    fn default() -> Self {
        Self::IDENTITY
    }
}

/// Independently prescribed pressure offset on one original opening plane, Pa.
///
/// Attach to `UnsteadyConfig::opening_pressure_programs` in `Face3::ALL` order.
/// The offset is added AFTER the global differential factor:
/// `p_face(t) = p_ref0 + factor(t)*(p_authored - p_ref0) + offset_side(t)`.
/// It translates the plane's original spatial profile without scaling it.
/// One to sixteen finite signed-Pa knots use linear interpolation and endpoint
/// holds, sampled once per original step endpoint. A single knot is constant.
///
/// The native solver subtracts the offset of each component's own anchor
/// opening before momentum and pressure work; common pressure changes perform
/// no incompressible work. Reported pressure restores that current reference.
/// Thermal backflow follows the same projected signed flux. This does not
/// change source thermodynamic pressure, interpolate between distinct openings
/// on one plane, or relax the original physical-heating admission rules.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PressureOffsetProgram {
    values: PressureProgram,
}

impl PressureOffsetProgram {
    /// No additional reservoir pressure.
    pub const ZERO: Self = Self {
        values: PressureProgram {
            points: [(0.0, 0.0); 16],
            len: 1,
        },
    };

    /// Retain `(time_s, offset_pa)` knots without a heap allocation.
    ///
    /// # Errors
    /// [`ChtError::InvalidInput`] for zero or more than sixteen points,
    /// non-finite entries, negative time, or non-increasing knot times.
    pub fn new(points: &[(f64, f64)]) -> Result<Self, ChtError> {
        Ok(Self {
            values: PressureProgram::from_points(points, "unsteady.opening_pressure_program")?,
        })
    }

    /// Original ordered `(time_s, offset_pa)` knots.
    #[must_use]
    pub fn points(&self) -> &[(f64, f64)] {
        self.values.points()
    }

    /// Offset in Pa at a finite, non-negative time, with endpoint holds.
    ///
    /// # Errors
    /// [`ChtError::InvalidInput`] for non-finite or negative sample time.
    pub fn sample_pa(&self, time_s: f64) -> Result<f64, ChtError> {
        self.values
            .sample_at(time_s, "unsteady.opening_pressure_program.time_s")
    }
}

impl Default for PressureOffsetProgram {
    fn default() -> Self {
        Self::ZERO
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn interpolation_holds_endpoints_and_preserves_signed_authored_knots() {
        let knots = [(1.0, 2.0), (3.0, 0.0), (5.0, -4.0)];
        let program = PressureProgram::new(&knots).unwrap();
        assert_eq!(program.points(), &knots);
        for (time, expected) in [
            (0.0, 2.0),
            (1.0, 2.0),
            (2.0, 1.0),
            (3.0, 0.0),
            (4.0, -2.0),
            (5.0, -4.0),
            (f64::MAX, -4.0),
        ] {
            assert_eq!(program.sample(time).unwrap(), expected);
        }
        for time in [0.0, 0.5, f64::MAX] {
            assert_eq!(
                PressureProgram::default().sample(time).unwrap().to_bits(),
                1.0_f64.to_bits()
            );
        }
        let constant = PressureProgram::new(&[(1.0, -2.0)]).unwrap();
        assert_eq!(constant.sample(0.0).unwrap(), -2.0);
        assert_eq!(constant.sample(f64::MAX).unwrap(), -2.0);
    }

    #[test]
    fn extreme_finite_interpolation_cannot_create_an_infinite_factor() {
        let max = f64::MAX;
        for values in [
            (max, -max),
            (-max, max),
            (max, 0.5 * max),
            (-max, -0.5 * max),
        ] {
            let program = PressureProgram::new(&[(0.0, values.0), (max, values.1)]).unwrap();
            for fraction in [0.0, 0.25, 0.5, 0.75, 1.0] {
                let factor = program.sample(fraction * max).unwrap();
                assert!(factor.is_finite());
                assert!(factor >= values.0.min(values.1));
                assert!(factor <= values.0.max(values.1));
            }
        }
        let cancel = PressureProgram::new(&[(0.0, max), (2.0, -max)]).unwrap();
        assert_eq!(cancel.sample(1.0).unwrap(), 0.0);
    }

    #[test]
    fn opening_offset_programs_preserve_pa_units_and_independent_history() {
        let left = PressureOffsetProgram::new(&[(0.0, 2.0), (2.0, -2.0)]).unwrap();
        let right = PressureOffsetProgram::new(&[(1.0, -1.0), (3.0, 3.0)]).unwrap();
        assert_eq!(left.points(), &[(0.0, 2.0), (2.0, -2.0)]);
        for (time, left_pa, right_pa) in [
            (0.0, 2.0, -1.0),
            (1.0, 0.0, -1.0),
            (2.0, -2.0, 1.0),
            (4.0, -2.0, 3.0),
        ] {
            assert_eq!(left.sample_pa(time).unwrap(), left_pa);
            assert_eq!(right.sample_pa(time).unwrap(), right_pa);
            assert_eq!(
                PressureOffsetProgram::default().sample_pa(time).unwrap(),
                0.0
            );
        }
        assert!(matches!(
            PressureOffsetProgram::new(&[(0.0, f64::NAN)]),
            Err(ChtError::InvalidInput {
                field: "unsteady.opening_pressure_program",
                ..
            })
        ));
        assert!(matches!(
            left.sample_pa(f64::INFINITY),
            Err(ChtError::InvalidInput {
                field: "unsteady.opening_pressure_program.time_s",
                ..
            })
        ));
    }

    #[test]
    fn malformed_program_and_sample_times_are_refused() {
        for knots in [
            vec![],
            vec![(0.0, 1.0); 17],
            vec![(-1.0, 1.0)],
            vec![(0.0, 1.0), (0.0, 2.0)],
            vec![(1.0, 1.0), (0.0, 2.0)],
            vec![(f64::INFINITY, 1.0)],
            vec![(f64::NAN, 1.0)],
            vec![(0.0, f64::INFINITY)],
            vec![(0.0, f64::NAN)],
        ] {
            assert!(matches!(
                PressureProgram::new(&knots),
                Err(ChtError::InvalidInput {
                    field: "unsteady.pressure_program",
                    ..
                })
            ));
        }
        for time in [-1.0, f64::INFINITY, f64::NEG_INFINITY, f64::NAN] {
            assert!(PressureProgram::IDENTITY.sample(time).is_err());
        }
    }
}
