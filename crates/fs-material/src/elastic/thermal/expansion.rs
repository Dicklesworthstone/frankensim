//! Original instantaneous-CTE paths for the existing thermal potential.
//!
//! This caches source derivatives and retained finite increments. The original
//! expansion integration owner still supplies every endpoint free strain and
//! receipt. Continuous CTE knots preserve Psi/S/U; only alpha' and capacity
//! have one-sided limits. Elastic coefficient admission is unchanged.

use super::*;
use crate::state_point::{affine_temperature_slope, pinned_temperature_line};
use fs_matdb::{InterpolationPolicy, MAX_PROPERTY_LINE_KNOTS, PropertyValue};
use fs_math::dd::Dd;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use std::sync::Arc;

#[derive(Clone, Debug, PartialEq)]
enum Path {
    Affine { slope: f64 },
    Curve(Arc<[(f64, f64)]>),
}

#[derive(Clone, Debug, PartialEq)]
pub(super) struct ExpansionPath {
    pub(super) name: &'static str,
    range: [f64; 2],
    reference_temperature_k: f64,
    reference: f64,
    maximum_abs: f64,
    maximum_slope_abs: f64,
    path: Path,
}

impl ExpansionPath {
    pub(super) fn new(
        card: &MaterialCard,
        source: &ResolvedMaterialStatePoint,
        name: &'static str,
        range: [f64; 2],
    ) -> Result<Self, ThermoelasticError> {
        let reference = value(source, name, Dims([0, 0, 0, -1, 0, 0]))?;
        let property = source.property(name).ok_or(ThermoelasticError::Invalid(
            "original expansion coefficient is missing",
        ))?;
        let claim = card
            .claims()
            .claim(property.answer().receipt.selected)
            .ok_or(ThermoelasticError::Invalid(
                "original expansion claim is absent",
            ))?;
        let PropertyValue::Curve {
            abscissa, knots, ..
        } = &claim.value
        else {
            return Self::affine(card, source, name, range, reference);
        };
        if abscissa != "T" {
            // Preserve previously admitted affine fixed-context source slices.
            return Self::affine(card, source, name, range, reference);
        }
        let [lo, hi] = range;
        let reference_temperature_k = temperature(source)?;
        if !lo.is_finite()
            || !hi.is_finite()
            || lo <= 0.0
            || lo >= hi
            || !(lo..=hi).contains(&reference_temperature_k)
            || claim.interpolation != InterpolationPolicy::LinearInside
            || !(2..=MAX_PROPERTY_LINE_KNOTS).contains(&knots.len())
        {
            return Err(ThermoelasticError::Invalid(
                "expansion needs a bounded continuous linear-inside source and positive temperature interval",
            ));
        }
        let context = source.query();
        let absolute = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ));
        if context.axis_quantities().get("T") != Some(&absolute)
            || !property.has_only_fixed_external_conditions(Some("T"), false)
        {
            return Err(ThermoelasticError::Invalid(
                "expansion requires explicit absolute temperature and fixed external source conditions",
            ));
        }
        source.at_query_point(card, &context)?;
        let point = |t| {
            context.clone().with("T", t).map_err(|_| {
                ThermoelasticError::Invalid("expansion source query is unrepresentable")
            })
        };
        let admitted = pinned_temperature_line(
            card,
            &claim.key,
            "T",
            &point(lo)?,
            &point(hi)?,
            property.answer().receipt.selected,
        )?;
        if admitted.first().map(|x| x.0) != Some(lo)
            || admitted.last().map(|x| x.0) != Some(hi)
            || admitted.iter().any(|x| !x.1.is_finite())
        {
            return Err(ThermoelasticError::Invalid(
                "expansion source must cover its complete thermal interval",
            ));
        }
        let maximum_abs = admitted.iter().fold(0.0_f64, |m, x| m.max(x.1.abs()));
        // Derive from the ORIGINAL curve, not rounded clipped endpoints. This
        // retains a nonzero tangent even across a one-ULP temperature interval.
        if let Ok(slope) = affine_temperature_slope(knots, knots[0].0, knots[0].1) {
            if slope == 0.0 && knots.iter().any(|x| x.1 != reference) {
                return Err(ThermoelasticError::Invalid(
                    "nonconstant expansion derivative underflows",
                ));
            }
            return Ok(Self {
                name,
                range,
                reference_temperature_k,
                reference,
                maximum_abs,
                maximum_slope_abs: slope.abs(),
                path: Path::Affine { slope },
            });
        }
        let mut maximum_slope_abs = 0.0_f64;
        for pair in knots.windows(2) {
            if pair[1].0 > lo && pair[0].0 < hi {
                maximum_slope_abs = maximum_slope_abs.max(segment_slope(pair[0], pair[1])?.abs());
            }
        }
        Ok(Self {
            name,
            range,
            reference_temperature_k,
            reference,
            maximum_abs,
            maximum_slope_abs,
            path: Path::Curve(knots.clone().into()),
        })
    }

    fn affine(
        card: &MaterialCard,
        source: &ResolvedMaterialStatePoint,
        name: &'static str,
        range: [f64; 2],
        reference: f64,
    ) -> Result<Self, ThermoelasticError> {
        let (reference_temperature_k, slope, ends) =
            affine_source_path(card, source, name, "T", range, false).map_err(field)?;
        if slope == 0.0 && ends[0] != ends[1] {
            return Err(ThermoelasticError::Invalid(
                "nonconstant expansion derivative underflows",
            ));
        }
        Ok(Self {
            name,
            range,
            reference_temperature_k,
            reference,
            maximum_abs: ends[0].abs().max(ends[1].abs()),
            maximum_slope_abs: slope.abs(),
            path: Path::Affine { slope },
        })
    }

    pub(super) fn maximum_abs(&self) -> f64 {
        self.maximum_abs
    }
    pub(super) fn maximum_slope_abs(&self) -> f64 {
        self.maximum_slope_abs
    }

    /// Actual native alpha, source-side alpha', and classical smoothness.
    /// Interior knots choose the increasing-temperature side; the admitted
    /// upper endpoint uses its interior, decreasing-temperature side.
    pub(super) fn at(
        &self,
        source: &ResolvedMaterialStatePoint,
        t: f64,
    ) -> Result<(f64, f64, bool), ThermoelasticError> {
        let actual = value(source, self.name, Dims([0, 0, 0, -1, 0, 0]))?;
        match &self.path {
            Path::Affine { slope } => {
                let expected =
                    finite(slope.mul_add(t - self.reference_temperature_k, self.reference))?;
                if (actual - expected).abs()
                    > 512.0 * f64::EPSILON * self.maximum_abs.max(actual.abs()).max(expected.abs())
                {
                    return Err(ThermoelasticError::Invalid(
                        "expansion source left its admitted affine path",
                    ));
                }
                Ok((actual, *slope, true))
            }
            Path::Curve(knots) => {
                let right = knots.partition_point(|x| x.0 <= t).saturating_sub(1);
                let mut index = right.min(knots.len() - 2);
                if t == self.range[1] && index > 0 && knots[index].0 == t {
                    index -= 1;
                }
                let slope = segment_slope(knots[index], knots[index + 1])?;
                let smooth = if right > 0 && right + 1 < knots.len() && knots[right].0 == t {
                    segment_slope(knots[right - 1], knots[right])?
                        == segment_slope(knots[right], knots[right + 1])?
                } else {
                    true
                };
                Ok((actual, slope, smooth))
            }
        }
    }

    /// [integral(alpha dT), alpha(T1)-alpha(T0)] from original source segments.
    /// No rounded endpoint free strain or coefficient is subtracted.
    pub(super) fn increments(
        &self,
        t0: f64,
        t1: f64,
        alpha0: f64,
        alpha1: f64,
    ) -> Result<[Dd; 2], ThermoelasticError> {
        let dd = Dd::from_f64;
        let dt = dd(t1) - dd(t0);
        if t0 == t1 {
            return Ok([Dd::ZERO; 2]);
        }
        if let Path::Affine { slope } = &self.path {
            return Ok([dd(0.5) * (dd(alpha0) + dd(alpha1)) * dt, dd(*slope) * dt]);
        }
        let Path::Curve(knots) = &self.path else {
            unreachable!()
        };
        let (mut lo, hi, sign) = if t0 < t1 {
            (t0, t1, 1.0)
        } else {
            (t1, t0, -1.0)
        };
        let mut index = knots
            .partition_point(|x| x.0 <= lo)
            .saturating_sub(1)
            .min(knots.len() - 2);
        let mut free = Dd::ZERO;
        let mut change = Dd::ZERO;
        while lo < hi {
            let (a, b) = (knots[index], knots[index + 1]);
            let end = hi.min(b.0);
            if lo < a.0 || end <= lo {
                return Err(ThermoelasticError::Invalid(
                    "finite expansion interval left original source support",
                ));
            }
            let slope = (dd(b.1) - dd(a.1)) / (dd(b.0) - dd(a.0));
            let at = |t: f64| {
                if t - a.0 <= b.0 - t {
                    dd(a.1) + slope * (dd(t) - dd(a.0))
                } else {
                    dd(b.1) - slope * (dd(b.0) - dd(t))
                }
            };
            let width = dd(end) - dd(lo);
            free = free + dd(0.5) * (at(lo) + at(end)) * width;
            change = change + slope * width;
            lo = end;
            if lo < hi {
                index += 1;
            }
        }
        let result = [dd(sign) * free, dd(sign) * change];
        if result.iter().any(|x| !x.to_f64().is_finite()) {
            return Err(ThermoelasticError::Invalid(
                "retained expansion increment is unrepresentable",
            ));
        }
        Ok(result)
    }
}

fn segment_slope(a: (f64, f64), b: (f64, f64)) -> Result<f64, ThermoelasticError> {
    let slope = finite((b.1 - a.1) / (b.0 - a.0))?;
    if slope == 0.0 && b.1 != a.1 {
        return Err(ThermoelasticError::Invalid(
            "nonconstant expansion derivative underflows",
        ));
    }
    Ok(slope)
}
