//! Shared scalar reductions for complete native thermal-trajectory ensembles.
//! These are empirical moments, never confidence or optional-stopping bounds.

use super::{CheckedSum, OpticalHeatingError};
use fs_math::det;
use fs_thermal::enthalpy::DepositionError;

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(super) struct ScalarEstimate {
    pub mean: f64,
    pub standard_error: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub(super) enum StatisticsError {
    Invalid(&'static str),
    Heating(OpticalHeatingError),
}
impl From<OpticalHeatingError> for StatisticsError {
    fn from(error: OpticalHeatingError) -> Self {
        Self::Heating(error)
    }
}

/// Reduce real native coordinates or compensated interval energy increments.
/// A caller-provided paired difference preserves correction terms that would
/// otherwise vanish if two large cumulative accounts were rounded first.
pub(super) fn measure_scalar<E>(
    count: usize,
    base: Option<f64>,
    paired: bool,
    multiplier: f64,
    mut values: impl FnMut(usize) -> Result<(f64, Option<f64>, Option<f64>), E>,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<(ScalarEstimate, f64, f64, f64), E>
where
    E: From<StatisticsError> + From<DepositionError>,
{
    let mut fine = Moments::default();
    let mut coarse = Moments::default();
    let mut difference = Moments::default();
    let origin = base.unwrap_or(0.0);
    let mut magnitude = origin.abs();
    for replica in 0..count {
        checkpoint()?;
        let (a, b, paired_difference) = values(replica)?;
        magnitude = magnitude.max(a.abs());
        fine.observe(finite(a - origin)?)?;
        if let Some(b) = b {
            magnitude = magnitude.max(b.abs());
            coarse.observe(finite(b - origin)?)?;
            difference.observe(match paired_difference {
                Some(value) => value,
                None => finite(a - b)?,
            })?;
        }
    }
    let mut estimate = fine.estimate()?;
    let scale = if base.is_some() {
        estimate.mean.abs().max(if paired {
            coarse.estimate()?.mean.abs()
        } else {
            0.0
        })
    } else {
        0.0
    };
    let width = if paired {
        let delta = difference.estimate()?;
        finite(delta.mean.abs() + multiplier * delta.standard_error)?
    } else {
        finite(multiplier * estimate.standard_error)?
    };
    estimate.mean = finite(origin + estimate.mean)?;
    Ok((estimate, width, scale, magnitude))
}

/// Signed, centered and scaled Welford moments. Centering avoids losing thermal
/// variability against an arbitrary caloric datum; scaling avoids squaring SI
/// magnitudes that would overflow or underflow the second central moment.
#[derive(Clone, Copy, Default)]
struct Moments {
    count: usize,
    anchor: f64,
    scale: f64,
    mean: f64,
    m2: f64,
}
impl Moments {
    fn observe(&mut self, value: f64) -> Result<(), StatisticsError> {
        finite(value)?;
        if self.count == 0 {
            self.anchor = value;
        }
        let shifted = finite(value - self.anchor)?;
        if shifted.abs() > self.scale {
            let factor = self.scale / shifted.abs();
            self.mean *= factor;
            self.m2 *= factor * factor;
            self.scale = shifted.abs();
        }
        let normalized = if self.scale == 0.0 {
            0.0
        } else {
            shifted / self.scale
        };
        self.count += 1;
        let delta = normalized - self.mean;
        self.mean += delta / self.count as f64;
        self.m2 += delta * (normalized - self.mean);
        Ok(())
    }
    fn estimate(self) -> Result<ScalarEstimate, StatisticsError> {
        if self.count < 2 {
            return Err(StatisticsError::Invalid("independent replicas required"));
        }
        let standard_error = finite(
            self.scale * det::sqrt(self.m2.max(0.0) / self.count as f64 / (self.count - 1) as f64),
        )?;
        // A positive variance must not become an exact-zero uncertainty gate.
        // In particular, k*SE can be representable even when SE itself rounds
        // to zero. Refuse that unresolved statistic before applying k; keep the
        // ordinary finite-scale calculation and the caller's allowance intact.
        if self.scale > 0.0 && self.m2 > 0.0 && standard_error == 0.0 {
            return Err(StatisticsError::Invalid(
                "positive ensemble variance has an unrepresentable standard error",
            ));
        }
        Ok(ScalarEstimate {
            mean: finite(self.anchor + self.mean * self.scale)?,
            standard_error,
        })
    }
}

pub(super) fn sum_difference(a: CheckedSum, b: CheckedSum) -> Result<f64, StatisticsError> {
    let mut difference = CheckedSum::default();
    difference.add(a.sum - b.sum)?;
    difference.add(a.correction)?;
    difference.add(-b.correction)?;
    Ok(difference.value()?)
}

fn finite(value: f64) -> Result<f64, StatisticsError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(StatisticsError::Invalid(
            "unrepresentable ensemble statistic",
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Debug)]
    struct Refused;
    impl From<StatisticsError> for Refused {
        fn from(_: StatisticsError) -> Self {
            Self
        }
    }
    impl From<DepositionError> for Refused {
        fn from(_: DepositionError) -> Self {
            Self
        }
    }

    #[test]
    fn paired_heat_statistics_keep_interval_corrections_below_historic_ulp() {
        let mut base = CheckedSum::default();
        base.add(1e25).unwrap();
        let mut fine = [base; 2];
        let mut coarse = [base; 2];
        fine[0].add(-0.125).unwrap();
        fine[1].add(-0.25).unwrap();
        coarse[0].add(-0.0625).unwrap();
        coarse[1].add(-0.125).unwrap();
        assert_eq!(fine[0].value().unwrap(), coarse[0].value().unwrap());
        let (estimate, width, scale, _) = measure_scalar::<Refused>(
            2,
            Some(0.0),
            true,
            3.0,
            |i| {
                Ok((
                    sum_difference(fine[i], base)?,
                    Some(sum_difference(coarse[i], base)?),
                    Some(sum_difference(fine[i], coarse[i])?),
                ))
            },
            &mut || Ok(()),
        )
        .unwrap();
        assert_eq!(estimate.mean, -0.1875);
        assert_eq!(estimate.standard_error, 0.0625);
        assert_eq!(width, 0.1875);
        assert_eq!(scale, 0.1875);
    }

    #[test]
    fn nonzero_subnormal_variance_cannot_report_zero_multiplied_error() {
        let minimum = f64::from_bits(1);
        // These two samples have SE=minimum/2, which rounds to zero, while
        // 3*SE rounds to 2*minimum. Reporting width0 would incorrectly pass
        // a requested allowance of minimum rather than refuse the lost SE.
        let result = measure_scalar::<Refused>(
            2,
            Some(0.0),
            false,
            3.0,
            |index| Ok((if index == 0 { 0.0 } else { minimum }, None, None)),
            &mut || Ok(()),
        );
        assert!(matches!(result, Err(Refused)));
    }
}
