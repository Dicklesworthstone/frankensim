//! Direct reciprocal heat exchange between two finite original phase charts.
//!
//! One signed Q debits one fixed mass and credits the other. For held
//! conductance G, X=G*dt satisfies integral_0^Q dq/(T_hot(q)-T_cold(q)).
//! Original sensible Cp intervals and latent plateaus are retained. Affine
//! T(h) paths advance analytically; nonlinear original-Cp paths use bounded
//! source-split quadrature and scalar inversion, with explicit error accounts.
//!
//! This conserves the pair's declared enthalpy account under the original
//! pressure/work convention. Density and volume changes do not supply pV work
//! or establish constant-volume internal-energy conservation. There is no
//! external supplied heat, mass transfer, evolving conductance or spatial model.

use crate::phase::{
    EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState, PhaseStateError, caloric,
};
use fs_math::{dd::Dd, det};

#[cfg(test)]
mod tests;

/// One fixed mass and its original accepted caloric state.
#[derive(Clone, Copy, Debug)]
pub struct FinitePhaseBody<'a> {
    pub curve: &'a EquilibriumEnthalpyPhaseCurve,
    pub initial: EquilibriumPhaseState,
    pub mass_kg: f64,
}

/// Physical error gates and whole-call work bounds.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FinitePhaseExchangeControls {
    pub absolute_heat_tolerance_j: f64,
    pub relative_heat_tolerance: f64,
    pub absolute_exposure_tolerance_j_k: f64,
    pub relative_exposure_tolerance: f64,
    /// Sum of prepared original intervals across both bodies. The final
    /// entropy query visits at most this many intervals again.
    pub maximum_chart_intervals: usize,
    pub maximum_root_iterations: usize,
    /// Total nonlinear Simpson panels across every trial and source probe.
    pub maximum_quadrature_panels: usize,
}

impl Default for FinitePhaseExchangeControls {
    fn default() -> Self {
        Self {
            absolute_heat_tolerance_j: 1e-7,
            relative_heat_tolerance: 1e-10,
            absolute_exposure_tolerance_j_k: 1e-9,
            relative_exposure_tolerance: 1e-10,
            maximum_chart_intervals: 16_384,
            maximum_root_iterations: 96,
            maximum_quadrature_panels: 131_072,
        }
    }
}

/// Atomic proposal in the supplied first/second body order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FinitePhaseExchange {
    pub initial: [EquilibriumPhaseState; 2],
    pub states: [EquilibriumPhaseState; 2],
    pub conductance_time_j_k: f64,
    /// Positive from the first body to the second; the only exchanged heat.
    pub heat_first_to_second_j: f64,
    /// Retained relative heat per mass, negative for the donating body.
    pub specific_heat_j_kg: [f64; 2],
    pub represented_energy_change_j: [f64; 2],
    /// Represented stored energy minus the signed heat supplied to each body.
    /// These are separate from trajectory error and require a caller energy gate.
    pub energy_residual_j: [f64; 2],
    /// Original-chart entropy changes of the represented accepted states.
    /// Latent entropy is included exactly once; no endpoint-flux quadrature
    /// replaces the original caloric integral.
    pub entropy_change_j_k: [f64; 2],
    pub total_entropy_change_j_k: f64,
    /// Includes represented-state energy errors and trajectory heat allowance.
    /// The signed entropy total is never clipped to hide numerical error.
    pub entropy_allowance_j_k: f64,
    /// Required exposure minus requested exposure. On the analytic path the
    /// remaining exposure is used directly, so this residual is zero.
    pub exposure_residual_j_k: f64,
    /// Quadrature truncation plus explicit floating-point allowance. This is
    /// not an interval-arithmetic certificate.
    pub exposure_error_allowance_j_k: f64,
    /// Flow-solution error, excluding the separately exposed stored-h residuals.
    pub heat_error_bound_j: f64,
    pub chart_intervals: usize,
    pub quadrature_panels: usize,
    /// Zero for the analytic affine-T(h)/latent path and stationary states.
    pub root_iterations: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub enum FinitePhaseExchangeError {
    Invalid(&'static str),
    Phase(PhaseStateError),
    OutsideSourceSupport,
    Unresolved(&'static str),
    WorkLimit,
    NonConvergence,
    Cancelled,
}

impl From<PhaseStateError> for FinitePhaseExchangeError {
    fn from(error: PhaseStateError) -> Self {
        if matches!(error, PhaseStateError::EntropyEvaluationCancelled) {
            Self::Cancelled
        } else {
            Self::Phase(error)
        }
    }
}
impl core::fmt::Display for FinitePhaseExchangeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite phase exchange: {self:?}")
    }
}
impl std::error::Error for FinitePhaseExchangeError {}

type Result<T> = core::result::Result<T, FinitePhaseExchangeError>;

/// Integrate isolated fixed-conductance contact on both original caloric paths.
///
/// X=G*dt has units J/K. Each body has its own mass, material identity, pressure
/// context and enthalpy datum. No equilibrium crossing or source extrapolation
/// is admitted. A source boundary with positive temperature gap has finite
/// available exposure; a vanishing gap is approached asymptotically.
///
/// Nonlinear intervals use an analytic Simpson truncation bound and bounded
/// inversion, not a timestep or a replacement h/T chord. Original intervals
/// are prepared once. No histories or accepted source states are mutated.
/// The caller retains spatial lumping, geometry and both energy-residual gates.
pub fn exchange(
    first: FinitePhaseBody<'_>,
    second: FinitePhaseBody<'_>,
    conductance_time_j_k: f64,
    controls: FinitePhaseExchangeControls,
) -> Result<FinitePhaseExchange> {
    exchange_with_checkpoint(first, second, conductance_time_j_k, controls, || false)
}

/// Same exchange with cancellation throughout preparation, quadrature,
/// inversion and the final original entropy queries.
pub fn exchange_with_checkpoint(
    first: FinitePhaseBody<'_>,
    second: FinitePhaseBody<'_>,
    conductance_time_j_k: f64,
    controls: FinitePhaseExchangeControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<FinitePhaseExchange> {
    let mut work = Work {
        cancelled: &mut cancelled,
        controls,
        charts: 0,
        panels: 0,
    };
    work.poll()?;
    controls.validate()?;
    nonnegative(conductance_time_j_k, "conductance-time exposure")?;
    for body in [first, second] {
        positive(body.mass_kg, "fixed body mass")?;
        body.curve.advance_specific_energy(body.initial, 0.0)?;
    }
    let first_hot = first.initial.temperature_k() >= second.initial.temperature_k();
    if conductance_time_j_k == 0.0
        || first.initial.temperature_k() == second.initial.temperature_k()
    {
        return finish(
            [first, second],
            None,
            conductance_time_j_k,
            Solution::default(),
            &mut work,
        );
    }
    let (hot, cold) = if first_hot {
        (first, second)
    } else {
        (second, first)
    };
    let hot = Path::new(hot, -1.0, &mut work)?;
    let cold = Path::new(cold, 1.0, &mut work)?;
    let pair = Pair::new(hot, cold);
    let solution = if pair.affine() {
        pair.analytic(conductance_time_j_k, &mut work)?
    } else {
        pair.nonlinear(conductance_time_j_k, &mut work)?
    };
    finish(
        [first, second],
        Some((&pair, first_hot)),
        conductance_time_j_k,
        solution,
        &mut work,
    )
}

impl FinitePhaseExchangeControls {
    fn validate(self) -> Result<()> {
        positive(self.absolute_heat_tolerance_j, "absolute heat tolerance")?;
        positive(
            self.absolute_exposure_tolerance_j_k,
            "absolute exposure tolerance",
        )?;
        nonnegative(self.relative_heat_tolerance, "relative heat tolerance")?;
        nonnegative(
            self.relative_exposure_tolerance,
            "relative exposure tolerance",
        )?;
        if self.relative_heat_tolerance >= 1.0
            || self.relative_exposure_tolerance >= 1.0
            || self.maximum_chart_intervals == 0
            || self.maximum_root_iterations == 0
            || self.maximum_quadrature_panels == 0
        {
            return Err(FinitePhaseExchangeError::Invalid(
                "finite exchange work controls",
            ));
        }
        Ok(())
    }
    fn heat_tolerance(self, heat: f64) -> Result<f64> {
        positive(
            self.absolute_heat_tolerance_j + self.relative_heat_tolerance * heat.abs(),
            "combined heat tolerance",
        )
    }
    fn exposure_tolerance(self, exposure: f64) -> Result<f64> {
        positive(
            self.absolute_exposure_tolerance_j_k + self.relative_exposure_tolerance * exposure,
            "combined exposure tolerance",
        )
    }
}

struct Work<'a> {
    cancelled: &'a mut dyn FnMut() -> bool,
    controls: FinitePhaseExchangeControls,
    charts: usize,
    panels: usize,
}
impl Work<'_> {
    fn poll(&mut self) -> Result<()> {
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
struct Segment {
    start_q: f64,
    end_q: f64,
    start_t: f64,
    end_t: f64,
    end_h: f64,
    /// Original specific Cp in traversal order; zero denotes a plateau.
    capacities: [f64; 2],
    mass_kg: f64,
    specific_width: f64,
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
        // Neither an absolute h datum nor a rounded m*Cp product replaces Cp.
        let local = caloric::temperature_fraction_from_bounds(
            [0.0, self.specific_width],
            (self.end_t - self.start_t).abs(),
            capacities,
            energy / self.mass_kg,
        )?;
        let fraction = if reverse { 1.0 - local } else { local };
        Ok(Point {
            temperature: caloric::interpolate(self.start_t, self.end_t, fraction),
            capacity: positive(
                self.mass_kg * caloric::capacity_at_fraction(self.capacities, fraction),
                "local total capacity",
            )?,
        })
    }
    fn slope(self) -> f64 {
        if self.capacities[0] == self.capacities[1] {
            0.0
        } else {
            self.mass_kg
                * ((self.capacities[1] - self.capacities[0]) / (self.end_t - self.start_t)).abs()
        }
    }
}

struct Path<'a> {
    body: FinitePhaseBody<'a>,
    direction: f64,
    segments: Vec<Segment>,
}
impl<'a> Path<'a> {
    fn new(body: FinitePhaseBody<'a>, direction: f64, work: &mut Work<'_>) -> Result<Self> {
        let knots = body.curve.knots();
        let initial_h = body.initial.specific_enthalpy_j_kg();
        let upper = if direction > 0.0 {
            knots.partition_point(|k| k.specific_enthalpy_j_kg <= initial_h)
        } else {
            knots.partition_point(|k| k.specific_enthalpy_j_kg < initial_h)
        };
        if upper == 0 || upper == knots.len() {
            return Err(FinitePhaseExchangeError::OutsideSourceSupport);
        }
        let mut index = upper - 1;
        let mut h = initial_h;
        let mut total = Dd::ZERO;
        let mut segments = Vec::new();
        loop {
            work.poll()?;
            if work.charts == work.controls.maximum_chart_intervals {
                return Err(FinitePhaseExchangeError::WorkLimit);
            }
            work.charts += 1;
            let (lower, upper) = (knots[index], knots[index + 1]);
            let (start, end) = if direction > 0.0 {
                (lower, upper)
            } else {
                (upper, lower)
            };
            let capacities = if lower.temperature_k == upper.temperature_k {
                [0.0; 2]
            } else {
                let cp = body
                    .curve
                    .sensible_heat_capacities_on_interval(index)
                    .unwrap_or(
                        [(upper.specific_enthalpy_j_kg - lower.specific_enthalpy_j_kg)
                            / (upper.temperature_k - lower.temperature_k);
                            2],
                    );
                positive(body.mass_kg * cp[0], "lower total capacity")?;
                positive(body.mass_kg * cp[1], "upper total capacity")?;
                if direction > 0.0 { cp } else { [cp[1], cp[0]] }
            };
            let start_q = if segments.is_empty() {
                -body.mass_kg * (start.specific_enthalpy_j_kg - initial_h).abs()
            } else {
                total.to_f64()
            };
            let increment = positive(
                body.mass_kg * (end.specific_enthalpy_j_kg - h).abs(),
                "original source heat interval",
            )?;
            total = total + Dd::from_f64(increment);
            let end_q = total.to_f64();
            if !start_q.is_finite() || !end_q.is_finite() || end_q <= start_q {
                return Err(unresolved("original source heat breakpoint"));
            }
            segments.push(Segment {
                start_q,
                end_q,
                start_t: start.temperature_k,
                end_t: end.temperature_k,
                end_h: end.specific_enthalpy_j_kg,
                capacities,
                mass_kg: body.mass_kg,
                specific_width: (end.specific_enthalpy_j_kg - start.specific_enthalpy_j_kg).abs(),
            });
            h = end.specific_enthalpy_j_kg;
            if direction > 0.0 {
                if index + 2 == knots.len() {
                    break;
                }
                index += 1;
            } else {
                if index == 0 {
                    break;
                }
                index -= 1;
            }
        }
        Ok(Self {
            body,
            direction,
            segments,
        })
    }
    fn limit(&self) -> f64 {
        self.segments[self.segments.len() - 1].end_q
    }
    fn after(&self, heat: f64) -> Segment {
        self.segments[self
            .segments
            .partition_point(|s| s.end_q <= heat)
            .min(self.segments.len() - 1)]
    }
    fn at(&self, heat: f64) -> Result<Point> {
        let index = self.segments.partition_point(|s| s.end_q < heat);
        self.segments
            .get(index)
            .ok_or(FinitePhaseExchangeError::OutsideSourceSupport)?
            .point(heat)
    }
    fn state(&self, heat: f64) -> Result<EquilibriumPhaseState> {
        if heat == 0.0 {
            return Ok(self.body.initial);
        }
        let index = self.segments.partition_point(|s| s.end_q < heat);
        let segment = self
            .segments
            .get(index)
            .ok_or(FinitePhaseExchangeError::OutsideSourceSupport)?;
        if heat == segment.end_q {
            Ok(self.body.curve.state_at_specific_enthalpy(segment.end_h)?)
        } else {
            Ok(self.body.curve.advance_specific_energy(
                self.body.initial,
                self.direction * heat / self.body.mass_kg,
            )?)
        }
    }
}

struct Pair<'a> {
    hot: Path<'a>,
    cold: Path<'a>,
    breaks: Vec<f64>,
}

#[derive(Clone, Copy, Default)]
struct Solution {
    heat: f64,
    residual: f64,
    allowance: f64,
    heat_error: f64,
    iterations: usize,
}
#[derive(Clone, Copy)]
struct Measure {
    exposure: f64,
    error: f64,
    complete: bool,
}

impl<'a> Pair<'a> {
    fn new(hot: Path<'a>, cold: Path<'a>) -> Self {
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
    fn driving(&self) -> f64 {
        self.hot.body.initial.temperature_k() - self.cold.body.initial.temperature_k()
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

fn finish(
    bodies: [FinitePhaseBody<'_>; 2],
    pair: Option<(&Pair<'_>, bool)>,
    exposure: f64,
    solution: Solution,
    work: &mut Work<'_>,
) -> Result<FinitePhaseExchange> {
    work.poll()?;
    let (states, signed_heat) = if let Some((pair, first_hot)) = pair {
        let hot = pair.hot.state(solution.heat)?;
        let cold = pair.cold.state(solution.heat)?;
        if hot.temperature_k() < cold.temperature_k() {
            return Err(unresolved(
                "represented pair state crosses thermal equilibrium",
            ));
        }
        if first_hot {
            ([hot, cold], solution.heat)
        } else {
            ([cold, hot], -solution.heat)
        }
    } else {
        ([bodies[0].initial, bodies[1].initial], 0.0)
    };
    let mut stored = [0.0; 2];
    let mut residuals = [0.0; 2];
    let mut specific = [0.0; 2];
    let mut entropy = [0.0; 2];
    let mut entropy_allowance = 0.0;
    let mut intervals_left = work.controls.maximum_chart_intervals;
    for (side, heat) in [-signed_heat, signed_heat].into_iter().enumerate() {
        work.poll()?;
        let body = bodies[side];
        specific[side] = heat / body.mass_kg;
        stored[side] = body.mass_kg
            * (states[side].specific_enthalpy_j_kg() - body.initial.specific_enthalpy_j_kg());
        residuals[side] = stored[side] - heat;
        if !specific[side].is_finite()
            || !stored[side].is_finite()
            || !residuals[side].is_finite()
            || (heat != 0.0 && (specific[side] == 0.0 || stored[side] == 0.0))
        {
            return Err(unresolved("represented finite-body heat"));
        }
        if heat != 0.0 {
            let mean = body.curve.thermodynamic_mean_temperature(
                body.initial.specific_enthalpy_j_kg(),
                states[side].specific_enthalpy_j_kg(),
                intervals_left,
                || (work.cancelled)(),
            )?;
            intervals_left -= mean.intervals;
            entropy[side] = body.mass_kg * mean.specific_entropy_change_j_kg_k;
        }
        entropy_allowance += 64.0 * f64::EPSILON * entropy[side].abs()
            + (residuals[side].abs() + solution.heat_error) / body.curve.knots()[0].temperature_k;
    }
    let total_entropy = (Dd::from_f64(entropy[0]) + Dd::from_f64(entropy[1])).to_f64();
    if !entropy.iter().all(|s| s.is_finite())
        || !total_entropy.is_finite()
        || !entropy_allowance.is_finite()
        || total_entropy < -entropy_allowance
    {
        return Err(unresolved("original pair entropy balance"));
    }
    work.poll()?;
    Ok(FinitePhaseExchange {
        initial: [bodies[0].initial, bodies[1].initial],
        states,
        conductance_time_j_k: exposure,
        heat_first_to_second_j: signed_heat,
        specific_heat_j_kg: specific,
        represented_energy_change_j: stored,
        energy_residual_j: residuals,
        entropy_change_j_k: entropy,
        total_entropy_change_j_k: total_entropy,
        entropy_allowance_j_k: entropy_allowance,
        exposure_residual_j_k: solution.residual,
        exposure_error_allowance_j_k: solution.allowance,
        heat_error_bound_j: solution.heat_error,
        chart_intervals: work.charts,
        quadrature_panels: work.panels,
        root_iterations: solution.iterations,
    })
}

fn inverse_capacity(capacity: f64) -> f64 {
    if capacity == 0.0 {
        0.0
    } else {
        capacity.recip()
    }
}
fn unresolved(what: &'static str) -> FinitePhaseExchangeError {
    FinitePhaseExchangeError::Unresolved(what)
}
fn positive(value: f64, field: &'static str) -> Result<f64> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(FinitePhaseExchangeError::Invalid(field))
    }
}
fn nonnegative(value: f64, field: &'static str) -> Result<()> {
    if value.is_finite() && value >= 0.0 {
        Ok(())
    } else {
        Err(FinitePhaseExchangeError::Invalid(field))
    }
}
