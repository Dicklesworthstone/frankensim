//! One bounded finite-body exchange driver for original caloric intervals.
//!
//! State identity, energy-coordinate admission and entropy stay with each owner.
//! This numerical layer stores original interval endpoints and capacities only;
//! its one extensive Q is never a second thermal state.

use super::{FinitePhaseExchangeError, FinitePhaseExchangeControls, Result, positive, unresolved};
use crate::phase::caloric;
use fs_math::{dd::Dd, det};

pub(crate) struct Work<'a> {
    pub(crate) cancelled: &'a mut dyn FnMut() -> bool,
    pub(crate) controls: FinitePhaseExchangeControls,
    pub(crate) charts: usize,
    pub(crate) panels: usize,
}
impl Work<'_> {
    pub(crate) fn poll(&mut self) -> Result<()> {
        if (self.cancelled)() {
            Err(FinitePhaseExchangeError::Cancelled)
        } else {
            Ok(())
        }
    }
    fn panel(&mut self) -> Result<()> {
        self.poll()?;
        if self.panels == self.controls.maximum_quadrature_panels {
            return Err(FinitePhaseExchangeError::WorkLimit);
        }
        self.panels += 1;
        Ok(())
    }
}

#[derive(Clone, Copy)]
pub(crate) struct Segment {
    pub(crate) start_q: f64,
    pub(crate) end_q: f64,
    pub(crate) start_t: f64,
    pub(crate) end_t: f64,
    /// Canonical original phase endpoint; native owners finalize their own e.
    pub(crate) end_energy: Option<f64>,
    /// Original source capacity in traversal order; zero denotes a plateau.
    /// energy_factor converts its energy unit into extensive body joules.
    pub(crate) capacities: [f64; 2],
    pub(crate) energy_factor: f64,
    pub(crate) source_energy_width: f64,
}
#[derive(Clone, Copy)]
struct Point {
    temperature: f64,
    capacity: f64,
}

impl Segment {
    fn point(self, heat: f64) -> Result<Point> {
        if heat < self.start_q || heat > self.end_q || !heat.is_finite() {
            return Err(FinitePhaseExchangeError::OutsideSourceSupport);
        }
        if self.capacities[0] == 0.0 {
            return Ok(Point {
                temperature: self.start_t,
                capacity: 0.0,
            });
        }
        let from_start = heat - self.start_q;
        let from_end = self.end_q - heat;
        let reverse = from_end < from_start;
        let (energy, capacities) = if reverse {
            (from_end, [self.capacities[1], self.capacities[0]])
        } else {
            (from_start, self.capacities)
        };
        // Query the existing inverse from the nearer ORIGINAL source endpoint.
        // Neither an absolute datum nor a rounded total-capacity product replaces
        // the original capacity or its source energy coordinate.
        let local = caloric::temperature_fraction_from_bounds(
            [0.0, self.source_energy_width],
            (self.end_t - self.start_t).abs(),
            capacities,
            energy / self.energy_factor,
        )?;
        let fraction = if reverse { 1.0 - local } else { local };
        Ok(Point {
            temperature: caloric::interpolate(self.start_t, self.end_t, fraction),
            capacity: positive(
                self.energy_factor * caloric::capacity_at_fraction(self.capacities, fraction),
                "local total capacity",
            )?,
        })
    }
    fn slope(self) -> f64 {
        if self.capacities[0] == self.capacities[1] {
            0.0
        } else {
            self.energy_factor
                * ((self.capacities[1] - self.capacities[0]) / (self.end_t - self.start_t)).abs()
        }
    }
}

/// Prepared ORIGINAL capacity intervals and their admitted extensive-heat limit.
/// The separate limit can clip native e support without a replacement T(e) chord.
pub(crate) struct NumericPath {
    pub(crate) initial_temperature: f64,
    pub(crate) segments: Vec<Segment>,
    pub(crate) limit: f64,
}
impl NumericPath {
    fn limit(&self) -> f64 { self.limit }
    fn after(&self, heat: f64) -> Segment {
        self.segments[self.segments.partition_point(|s| s.end_q <= heat)
            .min(self.segments.len() - 1)]
    }
}

pub(crate) struct Pair<'a> {
    hot: &'a NumericPath,
    cold: &'a NumericPath,
    breaks: Vec<f64>,
}

#[derive(Clone, Copy, Default)]
pub(crate) struct Solution {
    pub(crate) heat: f64,
    pub(crate) residual: f64,
    pub(crate) allowance: f64,
    pub(crate) heat_error: f64,
    pub(crate) iterations: usize,
}
#[derive(Clone, Copy)]
struct Measure {
    exposure: f64,
    error: f64,
    complete: bool,
}

impl<'a> Pair<'a> {
    pub(crate) fn new(hot: &'a NumericPath, cold: &'a NumericPath) -> Self {
        let limit = hot.limit().min(cold.limit());
        let mut breaks = vec![0.0, limit];
        breaks.extend(
            hot.segments
                .iter()
                .chain(&cold.segments)
                .map(|s| s.end_q)
                .filter(|q| *q > 0.0 && *q < limit),
        );
        breaks.sort_by(f64::total_cmp);
        breaks.dedup();
        Self { hot, cold, breaks }
    }
    pub(crate) fn solve(&self, exposure: f64, work: &mut Work<'_>) -> Result<Solution> {
        if self.affine() { self.analytic(exposure, work) } else { self.nonlinear(exposure, work) }
    }
    fn driving(&self) -> f64 {
        self.hot.initial_temperature - self.cold.initial_temperature
    }
    fn limit(&self) -> f64 {
        self.breaks[self.breaks.len() - 1]
    }
    fn affine(&self) -> bool {
        self.hot
            .segments
            .iter()
            .chain(&self.cold.segments)
            .all(|s| s.capacities[0] == s.capacities[1])
    }
    fn span(&self, left: f64, right: f64) -> Span {
        Span {
            hot: self.hot.after(left),
            cold: self.cold.after(left),
            left,
            right,
        }
    }

    fn analytic(&self, exposure: f64, work: &mut Work<'_>) -> Result<Solution> {
        let mut used = Dd::ZERO;
        let mut errors = Dd::ZERO;
        for bounds in self.breaks.windows(2) {
            work.poll()?;
            let span = self.span(bounds[0], bounds[1]);
            let (hot, cold) = span.points(span.left)?;
            let gap = hot.temperature - cold.temperature;
            if gap <= 0.0 {
                return Err(unresolved(
                    "analytic trajectory reached an unresolved equilibrium",
                ));
            }
            let remaining = exposure - used.to_f64();
            if !remaining.is_finite() || remaining <= 0.0 {
                return Err(unresolved("remaining analytic exposure"));
            }
            let slope = inverse_capacity(hot.capacity) + inverse_capacity(cold.capacity);
            if !slope.is_finite() {
                return Err(unresolved("analytic gap slope"));
            }
            let attenuation = slope * remaining;
            let loss = if slope == 0.0 {
                0.0
            } else {
                -det::expm1(-attenuation)
            };
            let effective = if slope == 0.0 || attenuation == 0.0 {
                remaining
            } else if attenuation <= 0.5 {
                remaining * (loss / attenuation)
            } else {
                loss / slope
            };
            let heat = gap * effective;
            let width = span.right - span.left;
            let end_gap = span.gap(span.right)?;
            // A final exponential may round onto its equilibrium. It never
            // consumes another interval or crosses that equilibrium.
            if heat < width || (heat == width && end_gap <= 0.0) {
                let total_heat = if heat == width {
                    span.right
                } else {
                    span.left + heat
                };
                positive(total_heat, "represented nonzero pair heat")?;
                let final_gap = if slope == 0.0 {
                    gap
                } else {
                    gap * det::exp(-attenuation)
                };
                let allowance = errors.to_f64() + 64.0 * f64::EPSILON * (exposure + used.to_f64());
                let heat_error = final_gap * allowance
                    + 64.0
                        * f64::EPSILON
                        * (total_heat.abs() + effective * (hot.temperature + cold.temperature));
                if !allowance.is_finite()
                    || !heat_error.is_finite()
                    || allowance > work.controls.exposure_tolerance(exposure)?
                    || heat_error > work.controls.heat_tolerance(total_heat)?
                {
                    return Err(unresolved("analytic pair error budget"));
                }
                return Ok(Solution {
                    heat: total_heat,
                    residual: 0.0,
                    allowance,
                    heat_error,
                    iterations: 0,
                });
            }
            if end_gap <= 0.0 {
                return Err(unresolved("analytic equilibrium rounding"));
            }
            let measure = span.affine_measure()?;
            used = used + Dd::from_f64(measure.exposure);
            errors = errors + Dd::from_f64(measure.error);
            let residual = used.to_f64() - exposure;
            let allowance = errors.to_f64();
            let heat_error = self.driving() * (residual.abs() + allowance);
            if residual.abs() + allowance <= work.controls.exposure_tolerance(exposure)?
                && heat_error <= work.controls.heat_tolerance(span.right)?
            {
                return Ok(Solution {
                    heat: span.right,
                    residual,
                    allowance,
                    heat_error,
                    iterations: 0,
                });
            }
            if residual >= 0.0 {
                return Err(unresolved("analytic boundary exposure rounding"));
            }
        }
        Err(FinitePhaseExchangeError::OutsideSourceSupport)
    }

    fn nonlinear(&self, exposure: f64, work: &mut Work<'_>) -> Result<Solution> {
        let drive = self.driving();
        let limit = self.limit();
        let exposure_tolerance = work.controls.exposure_tolerance(exposure)?;
        if exposure > limit / drive {
            let tolerance =
                exposure_tolerance.min(work.controls.heat_tolerance(limit)? / drive) / 8.0;
            match self.measure(limit, tolerance, exposure, work) {
                Ok(m) if m.complete && exposure > m.exposure + m.error => {
                    return Err(FinitePhaseExchangeError::OutsideSourceSupport);
                }
                Ok(_) | Err(FinitePhaseExchangeError::Unresolved("thermal pinch")) => {}
                Err(error) => return Err(error),
            }
        }
        let (mut lower, mut upper) = (0.0, limit);
        let mut heat = (exposure * drive).min(0.5 * upper);
        positive(heat, "nonzero pair heat")?;
        for iteration in 0..work.controls.maximum_root_iterations {
            work.poll()?;
            let heat_tolerance = work.controls.heat_tolerance(heat)?;
            let tolerance = exposure_tolerance.min(heat_tolerance / drive) / 8.0;
            let measured = match self.measure(heat, tolerance, exposure + exposure_tolerance, work)
            {
                Ok(m) => m,
                Err(FinitePhaseExchangeError::Unresolved("thermal pinch")) => {
                    upper = heat;
                    heat = lower + 0.5 * (upper - lower);
                    continue;
                }
                Err(error) => return Err(error),
            };
            let residual = measured.exposure - exposure;
            // The gap only decreases: dX/dQ=1/gap >= 1/initial_gap.
            let heat_error = (drive * (residual.abs() + measured.error)).min(upper - lower);
            if measured.complete
                && residual.abs() + measured.error <= exposure_tolerance
                && heat_error <= heat_tolerance
            {
                return Ok(Solution {
                    heat,
                    residual,
                    allowance: measured.error,
                    heat_error,
                    iterations: iteration + 1,
                });
            }
            if residual > measured.error {
                upper = heat;
            } else if residual < -measured.error {
                lower = heat;
            } else {
                return Err(unresolved("exposure residual below numerical resolution"));
            }
            let next = lower + 0.5 * (upper - lower);
            if next <= lower || next >= upper {
                return Err(unresolved("pair heat bracket"));
            }
            heat = next;
        }
        Err(FinitePhaseExchangeError::NonConvergence)
    }

    fn measure(
        &self,
        heat: f64,
        tolerance: f64,
        ceiling: f64,
        work: &mut Work<'_>,
    ) -> Result<Measure> {
        let mut sum = Dd::ZERO;
        let mut errors = Dd::ZERO;
        let mut pending = Vec::new();
        for bounds in self.breaks.windows(2) {
            if bounds[0] >= heat {
                break;
            }
            work.poll()?;
            let span = self.span(bounds[0], bounds[1].min(heat));
            let first_gap = span.gap(span.left)?;
            let last_gap = span.gap(span.right)?;
            if last_gap <= 0.0 || first_gap <= 0.0 {
                return Err(unresolved("thermal pinch"));
            }
            // The decreasing gap makes width/first_gap a lower time bound.
            // Incomplete bounds may reject a trial, never accept one.
            let (hot, cold) = span.points(span.left)?;
            let upper_gap = first_gap + 16.0 * f64::EPSILON * (hot.temperature + cold.temperature);
            let low = (span.right - span.left) / upper_gap;
            let lower_bound = sum.to_f64() - errors.to_f64() + low;
            let lower_error = 64.0 * f64::EPSILON * lower_bound.abs();
            if lower_bound - lower_error > ceiling {
                return Ok(Measure {
                    exposure: lower_bound,
                    error: lower_error,
                    complete: false,
                });
            }
            if span.affine() {
                let m = span.affine_measure()?;
                sum = sum + Dd::from_f64(m.exposure);
                errors = errors + Dd::from_f64(m.error);
            } else {
                pending.push((span.left, span.right));
                while let Some((left, right)) = pending.pop() {
                    work.panel()?;
                    let m = span.panel(left, right)?;
                    let allocation = tolerance * ((right - left) / heat);
                    if m.error <= allocation {
                        sum = sum + Dd::from_f64(m.exposure);
                        errors = errors + Dd::from_f64(m.error);
                    } else {
                        let middle = left + 0.5 * (right - left);
                        if middle <= left || middle >= right {
                            return Err(unresolved("pair quadrature precision"));
                        }
                        pending.push((middle, right));
                        pending.push((left, middle));
                    }
                    if sum.to_f64() - errors.to_f64() > ceiling {
                        return Ok(Measure {
                            exposure: sum.to_f64(),
                            error: errors.to_f64(),
                            complete: false,
                        });
                    }
                }
            }
            if sum.to_f64() - errors.to_f64() > ceiling {
                return Ok(Measure {
                    exposure: sum.to_f64(),
                    error: errors.to_f64(),
                    complete: false,
                });
            }
        }
        let exposure = positive(sum.to_f64(), "integrated pair exposure")?;
        let error = errors.to_f64();
        if !error.is_finite() {
            return Err(unresolved("pair exposure allowance"));
        }
        Ok(Measure {
            exposure,
            error,
            complete: true,
        })
    }
}

#[derive(Clone, Copy)]
struct Span {
    hot: Segment,
    cold: Segment,
    left: f64,
    right: f64,
}

impl Span {
    fn points(self, q: f64) -> Result<(Point, Point)> {
        Ok((self.hot.point(q)?, self.cold.point(q)?))
    }
    fn gap(self, q: f64) -> Result<f64> {
        let (hot, cold) = self.points(q)?;
        Ok(hot.temperature - cold.temperature)
    }
    fn affine(self) -> bool {
        self.hot.capacities[0] == self.hot.capacities[1]
            && self.cold.capacities[0] == self.cold.capacities[1]
    }
    fn affine_measure(self) -> Result<Measure> {
        let (ha, ca) = self.points(self.left)?;
        let (hb, cb) = self.points(self.right)?;
        let first = ha.temperature - ca.temperature;
        let last = hb.temperature - cb.temperature;
        if first <= 0.0 || last <= 0.0 {
            return Err(unresolved("thermal pinch"));
        }
        let ratio = (last - first) / first;
        let width = self.right - self.left;
        let exposure = if ratio == 0.0 {
            width / first
        } else if ratio.abs() < 0.5 {
            (width / first) * ratio.ln_1p() / ratio
        } else {
            width * (det::ln(last) - det::ln(first)) / (last - first)
        };
        positive(exposure, "affine pair exposure")?;
        let error = 64.0
            * f64::EPSILON
            * exposure
            * (1.0 + (ha.temperature + hb.temperature + ca.temperature + cb.temperature) / last);
        if !error.is_finite() {
            return Err(unresolved("affine pair exposure allowance"));
        }
        Ok(Measure {
            exposure,
            error,
            complete: true,
        })
    }
    fn panel(self, left: f64, right: f64) -> Result<Measure> {
        let width = right - left;
        let (ha, ca) = self.points(left)?;
        let (hb, cb) = self.points(right)?;
        let (hm, cm) = self.points(left + 0.5 * width)?;
        let (da, db, dm) = (
            ha.temperature - ca.temperature,
            hb.temperature - cb.temperature,
            hm.temperature - cm.temperature,
        );
        let minimum = da.min(db).min(dm);
        if minimum <= 0.0 {
            return Err(unresolved("thermal pinch"));
        }
        let rate = |a: Point, b: Point, slope: f64| {
            if a.capacity == 0.0 {
                (0.0, 0.0)
            } else {
                let minimum = a.capacity.min(b.capacity);
                let r = width / minimum;
                (r, (slope / minimum) * r)
            }
        };
        let (rh, ah) = rate(ha, hb, self.hot.slope());
        let (rc, ac) = rate(ca, cb, self.cold.slope());
        let v1 = (rh + rc) / minimum;
        let v2 = (rh * ah + rc * ac) / minimum;
        let v3 = 3.0 * (rh * ah * ah + rc * ac * ac) / minimum;
        let v4 = 15.0 * (rh * ah * ah * ah + rc * ac * ac * ac) / minimum;
        // Positive original capacities bound all four temperature derivatives.
        // This is Simpson's normalized width/2880 fourth-derivative bound.
        let truncation = (width / minimum) / 2880.0
            * (24.0 * v1.powi(4) + 36.0 * v1 * v1 * v2 + 6.0 * v2 * v2 + 8.0 * v1 * v3 + v4);
        let exposure = (width / 6.0) * (1.0 / da + 4.0 / dm + 1.0 / db);
        let scale = ha.temperature + hb.temperature + ca.temperature + cb.temperature;
        let roundoff = 64.0 * f64::EPSILON * exposure * (1.0 + scale / minimum);
        if !exposure.is_finite() || exposure <= 0.0 || truncation.is_nan() || roundoff.is_nan() {
            return Err(unresolved("pair quadrature arithmetic"));
        }
        Ok(Measure {
            exposure,
            error: truncation + roundoff,
            complete: true,
        })
    }
}

fn inverse_capacity(capacity: f64) -> f64 {
    if capacity == 0.0 {
        0.0
    } else {
        capacity.recip()
    }
}
