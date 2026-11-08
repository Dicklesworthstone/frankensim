//! Countercurrent water/steam exchange on the original continuous Gibbs chart.
//!
//! One reciprocal duty Q determines both outlet enthalpies. Along the water
//! direction q, the sensible stream has remaining heat progress |Q|-q. Required
//! UA is the integral of dq / |T_water(q)-T_sensible(|Q|-q)|. Original water
//! saturation endpoints and original Cp knots split every evaluation.
//!
//! Monotone original-temperature ranges admit the WHOLE gap on each accepted
//! panel, with explicit inverse/rounding allowances; endpoint gaps alone are
//! insufficient. The existing adaptive panel traversal is shared with the
//! single-phase exchanger. Its water panels use embedded Simpson error
//! ESTIMATES, not the affine-Cp owner's analytic fourth-derivative bound.
//! No resampled caloric law, latent Cp, new EOS or pressure work is introduced.

pub use super::co_current_water::WaterExchangerState;
use super::co_current_water::{CoCurrentWaterExchangerError, CoCurrentWaterHeatExchanger};
use super::phase::WaterPhaseExchangerConfig;
use super::{
    ExchangerConfig, ExchangerError, ExchangerStream, Measure, SourcePath, integrate_panels, poll,
    positive,
};
use fs_exec::Cx;
use fs_material::phase::water::{WaterError, WaterState};
use fs_math::dd::Dd;
use fs_qty::Temperature;
use fs_thermal::enthalpy::ThermalLaw;
use std::cell::Cell;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CountercurrentWaterExchangerConfig {
    pub endpoint: WaterPhaseExchangerConfig,
    /// Heat-bracket and UA tolerances, outer iterations and panels per UA
    /// query. The water quadrature uses an estimate, not an analytic bound.
    pub numerics: ExchangerConfig,
    /// Shared across preparation and ALL trial duties, including rejected
    /// pinches. Each quadrature water query has at most one paired Cp inverse.
    /// Final source calorimetry additionally uses its native bounded inverse.
    pub maximum_water_source_queries: usize,
}
impl Default for CountercurrentWaterExchangerConfig {
    fn default() -> Self {
        Self {
            endpoint: WaterPhaseExchangerConfig {
                absolute_heat_tolerance_w: 1e-5,
                ..WaterPhaseExchangerConfig::default()
            },
            numerics: ExchangerConfig {
                absolute_heat_tolerance_w: 1e-3,
                relative_heat_tolerance: 1e-8,
                absolute_conductance_tolerance_w_per_k: 1e-6,
                relative_conductance_tolerance: 1e-8,
                max_iterations: 80,
                max_quadrature_panels: 16_384,
            },
            maximum_water_source_queries: 250_000,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct CountercurrentWaterExchangerState {
    /// Original source endpoint/entropy audit. Its SSPRK2 trajectory is None;
    /// this owner solves the countercurrent boundary-value problem by UA.
    pub balance: WaterExchangerState,
    pub required_conductance_w_per_k: f64,
    pub quadrature_error_estimate_w_per_k: f64,
    pub source_inverse_allowance_w_per_k: f64,
    /// Numerical bracket, conditional on the reported quadrature estimate.
    pub heat_bracket_width_w: f64,
    /// Whole-panel monotone-source lower bound, with numerical allowances.
    /// This is not an interval-arithmetic certificate of the full calculation.
    pub minimum_temperature_gap_lower_bound_k: f64,
    pub iterations: usize,
    pub panels: usize,
    pub water_source_queries: usize,
    pub quadrature_enthalpy_iterations: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub enum CountercurrentWaterExchangerError {
    Caloric(ExchangerError),
    Water(WaterError),
    SourceBalance(CoCurrentWaterExchangerError),
    Interrupted,
}
impl From<ExchangerError> for CountercurrentWaterExchangerError {
    fn from(error: ExchangerError) -> Self {
        match error {
            ExchangerError::Interrupted => Self::Interrupted,
            other => Self::Caloric(other),
        }
    }
}
impl From<WaterError> for CountercurrentWaterExchangerError {
    fn from(error: WaterError) -> Self {
        match error {
            WaterError::Cancelled => Self::Interrupted,
            other => Self::Water(other),
        }
    }
}
impl From<CoCurrentWaterExchangerError> for CountercurrentWaterExchangerError {
    fn from(error: CoCurrentWaterExchangerError) -> Self {
        match error {
            CoCurrentWaterExchangerError::Interrupted => Self::Interrupted,
            other => Self::SourceBalance(other),
        }
    }
}
impl core::fmt::Display for CountercurrentWaterExchangerError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "original countercurrent water exchanger: {self:?}")
    }
}
impl std::error::Error for CountercurrentWaterExchangerError {}

/// Continuous countercurrent plug flow, prescribed UA and independent isobars.
/// The retained caloric composition supplies the SAME original endpoint audit;
/// its co-current time integrator is never invoked here.
#[derive(Clone, Debug, PartialEq)]
pub struct CountercurrentWaterHeatExchanger {
    caloric: CoCurrentWaterHeatExchanger,
}
impl CountercurrentWaterHeatExchanger {
    pub fn new(
        sensible: ExchangerStream,
        water_inlet: WaterState,
        water_mass_flow_kg_per_s: f64,
        conductance_w_per_k: f64,
    ) -> Result<Self, CountercurrentWaterExchangerError> {
        Ok(Self {
            caloric: CoCurrentWaterHeatExchanger::new(
                sensible,
                water_inlet,
                water_mass_flow_kg_per_s,
                conductance_w_per_k,
            )?,
        })
    }
    pub fn sensible_stream(&self) -> &ExchangerStream {
        self.caloric.sensible_stream()
    }
    pub fn water_inlet(&self) -> WaterState {
        self.caloric.water_inlet()
    }
    pub fn water_mass_flow_kg_per_s(&self) -> f64 {
        self.caloric.water_mass_flow_kg_per_s()
    }
    pub fn conductance_w_per_k(&self) -> f64 {
        self.caloric.conductance_w_per_k()
    }

    pub fn solve(
        &self,
        cx: &Cx<'_>,
        config: CountercurrentWaterExchangerConfig,
    ) -> Result<CountercurrentWaterExchangerState, CountercurrentWaterExchangerError> {
        self.validate(cx, config)?;
        let driving =
            (self.water_inlet().temperature_k() - self.sensible_stream().inlet_temperature_k).abs();
        let ua = self.conductance_w_per_k();
        if ua == 0.0 || driving == 0.0 {
            // Equal inlet temperatures do not identify UA from zero duty;
            // retain the prescribed value without claiming an inverse rating.
            return self.publish(cx, 0.0, UaMeasure::zero(ua, driving), 0.0, 0, None, config);
        }
        let pair = WaterPair::new(cx, self, config)?;
        let ua_tolerance = config.numerics.ua_tolerance(ua)?;
        let tolerance_at = |q| -> Result<f64, ExchangerError> {
            Ok(ua_tolerance.min(config.numerics.heat_tolerance(q)? / (4.0 * driving)) / 8.0)
        };
        let (mut lower, mut upper) = (0.0, pair.limit);
        if !pair.limit_is_pinch {
            match pair.measure(
                cx,
                upper,
                tolerance_at(upper)?,
                Some(ua + ua_tolerance),
                config,
            ) {
                Ok(measured) if ua > measured.total.ua + measured.total.error + ua_tolerance => {
                    return Err(ExchangerError::OutsideSourceSupport.into());
                }
                Ok(_) | Err(CountercurrentWaterExchangerError::Caloric(ExchangerError::Pinch)) => {}
                Err(error) => return Err(error),
            }
        }
        let mut q = positive((ua * driving).min(0.5 * upper), "nonzero trial heat")?;
        for iteration in 1..=config.numerics.max_iterations {
            poll(cx)?;
            let measured =
                match pair.measure(cx, q, tolerance_at(q)?, Some(ua + ua_tolerance), config) {
                    Ok(value) => value,
                    Err(CountercurrentWaterExchangerError::Caloric(ExchangerError::Pinch)) => {
                        upper = q;
                        q = midpoint(lower, upper)?;
                        continue;
                    }
                    Err(error) => return Err(error),
                };
            if !measured.total.complete {
                upper = q;
                q = midpoint(lower, upper)?;
                continue;
            }
            let residual = measured.total.ua - ua;
            // On the admitted path every gap <= the inlet driving gap, and
            // increasing Q reduces the opposite stream's local gap. Hence
            // dUA/dQ >= 1/driving. The reported numerical uncertainty gives
            // this additional bracket restriction; it remains an estimate.
            let heat_allowance = driving * (residual.abs() + measured.total.error);
            lower = lower.max(q - heat_allowance);
            upper = upper.min(q + heat_allowance);
            if !lower.is_finite() || !upper.is_finite() || upper < lower {
                return Err(ExchangerError::Unresolved("countercurrent heat bracket").into());
            }
            if residual.abs() + measured.total.error <= ua_tolerance
                && upper - lower <= config.numerics.heat_tolerance(q)?
            {
                return self.publish(
                    cx,
                    pair.sign * q,
                    measured,
                    upper - lower,
                    iteration,
                    Some(&pair),
                    config,
                );
            }
            if residual > measured.total.error {
                upper = upper.min(q);
            } else if residual < -measured.total.error {
                lower = lower.max(q);
            } else {
                return Err(ExchangerError::Unresolved(
                    "countercurrent UA below source resolution",
                )
                .into());
            }
            q = midpoint(lower, upper)?;
        }
        Err(ExchangerError::NonConvergence {
            iterations: config.numerics.max_iterations,
        }
        .into())
    }

    /// Rate a specified signed duty and return its required UA. This method
    /// does not use the object's prescribed UA. Whole-path pinch/source gates
    /// and original endpoint/entropy audits are identical to `solve`.
    pub fn rate_duty(
        &self,
        cx: &Cx<'_>,
        heat_rate_water_to_sensible_w: f64,
        config: CountercurrentWaterExchangerConfig,
    ) -> Result<CountercurrentWaterExchangerState, CountercurrentWaterExchangerError> {
        self.validate(cx, config)?;
        let heat = heat_rate_water_to_sensible_w;
        let gap = self.water_inlet().temperature_k() - self.sensible_stream().inlet_temperature_k;
        if !heat.is_finite()
            || (heat != 0.0 && (gap == 0.0 || heat.is_sign_positive() != gap.is_sign_positive()))
        {
            return Err(ExchangerError::Invalid("countercurrent heat direction").into());
        }
        if heat == 0.0 {
            return self.publish(
                cx,
                heat,
                UaMeasure::zero(0.0, gap.abs()),
                0.0,
                0,
                None,
                config,
            );
        }
        let pair = WaterPair::new(cx, self, config)?;
        if heat.abs() > pair.limit {
            return Err(ExchangerError::OutsideSourceSupport.into());
        }
        let tolerance = config.numerics.absolute_conductance_tolerance_w_per_k / 4.0;
        let measured = pair.measure(cx, heat.abs(), tolerance, None, config)?;
        self.publish(cx, heat, measured, 0.0, 0, Some(&pair), config)
    }

    fn validate(
        &self,
        cx: &Cx<'_>,
        config: CountercurrentWaterExchangerConfig,
    ) -> Result<(), CountercurrentWaterExchangerError> {
        poll(cx)?;
        config
            .endpoint
            .validate()
            .map_err(CoCurrentWaterExchangerError::from)?;
        config.numerics.validate()?;
        if config.maximum_water_source_queries == 0 {
            return Err(ExchangerError::Invalid("water query budget").into());
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn publish(
        &self,
        cx: &Cx<'_>,
        heat: f64,
        measured: UaMeasure,
        width: f64,
        iterations: usize,
        pair: Option<&WaterPair<'_>>,
        config: CountercurrentWaterExchangerConfig,
    ) -> Result<CountercurrentWaterExchangerState, CountercurrentWaterExchangerError> {
        let out = match pair {
            Some(pair) => pair.sensible.at(heat.abs())?.temperature,
            None => self.sensible_stream().inlet_temperature_k,
        };
        let balance = self.caloric.finish(cx, heat, out, None, config.endpoint)?;
        poll(cx)?;
        Ok(CountercurrentWaterExchangerState {
            balance,
            required_conductance_w_per_k: measured.total.ua,
            quadrature_error_estimate_w_per_k: measured.quadrature,
            source_inverse_allowance_w_per_k: measured.inverse,
            heat_bracket_width_w: width,
            minimum_temperature_gap_lower_bound_k: measured.minimum_gap,
            iterations,
            panels: pair.map_or(0, |p| p.panels.get()),
            water_source_queries: pair.map_or(0, |p| p.queries.get()),
            quadrature_enthalpy_iterations: pair.map_or(0, |p| p.enthalpy_iterations.get()),
        })
    }
}

fn midpoint(a: f64, b: f64) -> Result<f64, ExchangerError> {
    let mid = a + 0.5 * (b - a);
    if !mid.is_finite() || mid <= a || mid >= b {
        return Err(ExchangerError::Unresolved("countercurrent heat bracket"));
    }
    Ok(mid)
}

#[derive(Clone, Copy)]
struct UaMeasure {
    total: Measure,
    quadrature: f64,
    inverse: f64,
    minimum_gap: f64,
}
impl UaMeasure {
    fn zero(ua: f64, gap: f64) -> Self {
        Self {
            total: Measure {
                ua,
                ..Measure::default()
            },
            quadrature: 0.0,
            inverse: 0.0,
            minimum_gap: gap,
        }
    }
}
#[derive(Clone, Copy)]
struct Sample {
    water: f64,
    sensible: f64,
    water_error: f64,
    sensible_error: f64,
}
impl Sample {
    fn gap(self, sign: f64) -> f64 {
        sign * (self.water - self.sensible)
    }
    fn error(self) -> f64 {
        self.water_error + self.sensible_error
    }
}

struct WaterPair<'a> {
    water: WaterState,
    water_mass: f64,
    sensible: SourcePath<'a>,
    sign: f64,
    limit: f64,
    limit_is_pinch: bool,
    queries: Cell<usize>,
    enthalpy_iterations: Cell<usize>,
    panels: Cell<usize>,
    maximum_queries: usize,
}
impl<'a> WaterPair<'a> {
    fn new(
        cx: &Cx<'_>,
        owner: &'a CountercurrentWaterHeatExchanger,
        config: CountercurrentWaterExchangerConfig,
    ) -> Result<Self, CountercurrentWaterExchangerError> {
        let water = owner.water_inlet();
        let sensible = owner.sensible_stream();
        let sign = (water.temperature_k() - sensible.inlet_temperature_k).signum();
        let source = SourcePath::new(cx, sensible, water.temperature_k(), sign, config.numerics)?;
        let mut pair = Self {
            water,
            water_mass: owner.water_mass_flow_kg_per_s(),
            sensible: source,
            sign,
            limit: 0.0,
            limit_is_pinch: false,
            queries: Cell::new(0),
            enthalpy_iterations: Cell::new(0),
            panels: Cell::new(0),
            maximum_queries: config.maximum_water_source_queries,
        };
        let bounds = water.law().enthalpy_bounds_j_kg();
        let (low, _) = pair.water_state(cx, Dd::from_f64(bounds[0]))?;
        let (high, _) = pair.water_state(cx, Dd::from_f64(bounds[1]))?;
        let target = sensible
            .inlet_temperature_k
            .clamp(low.temperature_k(), high.temperature_k());
        pair.query()?;
        let boundary = if target == water.law().saturation().temperature().value() {
            water
                .law()
                .state_at_quality(if sign > 0.0 { 1.0 } else { 0.0 }, || {
                    cx.checkpoint().is_err()
                })?
        } else {
            water
                .law()
                .state_at_temperature(Temperature::new(target), || cx.checkpoint().is_err())?
        };
        let available = Dd::from_f64(sign)
            * Dd::from_f64(pair.water_mass)
            * (Dd::from_f64(water.specific_enthalpy_j_kg())
                - Dd::from_f64(boundary.specific_enthalpy_j_kg()));
        let mut water_limit = available.to_f64();
        if available.lt(Dd::from_f64(water_limit)) {
            water_limit = water_limit.next_down();
        }
        if water_limit <= 0.0 {
            return Err(ExchangerError::OutsideSourceSupport.into());
        }
        pair.limit = positive(
            water_limit.min(pair.sensible.limit()),
            "countercurrent source heat limit",
        )?;
        pair.limit_is_pinch = (pair.limit == water_limit && target == sensible.inlet_temperature_k)
            || (pair.limit == pair.sensible.limit()
                && pair.sensible.anchors.last().unwrap().temperature == water.temperature_k());
        Ok(pair)
    }
    fn query(&self) -> Result<(), CountercurrentWaterExchangerError> {
        if self.queries.get() >= self.maximum_queries {
            return Err(ExchangerError::WorkLimit.into());
        }
        self.queries.set(self.queries.get() + 1);
        Ok(())
    }
    fn water_state(
        &self,
        cx: &Cx<'_>,
        target: Dd,
    ) -> Result<(WaterState, f64), CountercurrentWaterExchangerError> {
        poll(cx)?;
        self.query()?;
        let bounds = self.water.law().enthalpy_bounds_j_kg();
        if target.lt(Dd::from_f64(bounds[0])) || Dd::from_f64(bounds[1]).lt(target) {
            return Err(ExchangerError::OutsideSourceSupport.into());
        }
        let solved = self
            .water
            .law()
            .state_at_specific_enthalpy(target.to_f64(), || cx.checkpoint().is_err())?;
        self.enthalpy_iterations
            .set(self.enthalpy_iterations.get() + solved.iterations());
        let state = solved.into_state();
        let h_error = (Dd::from_f64(state.specific_enthalpy_j_kg()) - target)
            .abs()
            .to_f64()
            + 64.0 * f64::EPSILON * (state.specific_enthalpy_j_kg().abs() + target.to_f64().abs());
        let error = h_error * self.water.law().maximum_temperature_slope()
            + 64.0 * f64::EPSILON * state.temperature_k();
        Ok((state, error))
    }
    fn point(
        &self,
        cx: &Cx<'_>,
        heat: f64,
        q: f64,
    ) -> Result<Sample, CountercurrentWaterExchangerError> {
        let target = Dd::from_f64(self.water.specific_enthalpy_j_kg())
            - Dd::from_f64(self.sign) * Dd::from_f64(q) / Dd::from_f64(self.water_mass);
        let (water, water_error) = self.water_state(cx, target)?;
        let progress = heat - q;
        let sensible = self.sensible.at(progress)?;
        let i = self
            .sensible
            .anchors
            .partition_point(|a| a.heat <= progress)
            .saturating_sub(1);
        let anchor = self.sensible.anchors[i];
        let dh = self
            .sensible
            .stream
            .caloric
            .specific_enthalpy_change_j_kg(anchor.temperature, sensible.temperature)
            .map_err(ExchangerError::from)?;
        let residual = (Dd::from_f64(anchor.heat)
            + Dd::from_f64(self.sign)
                * Dd::from_f64(self.sensible.stream.mass_flow_kg_per_s)
                * Dd::from_f64(dh)
            - Dd::from_f64(progress))
        .abs()
        .to_f64();
        let sensible_error = (residual + 64.0 * f64::EPSILON * progress)
            * self.sensible.stream.caloric.maximum_temperature_slope()
            / self.sensible.stream.mass_flow_kg_per_s
            + 64.0 * f64::EPSILON * sensible.temperature;
        Ok(Sample {
            water: water.temperature_k(),
            sensible: sensible.temperature,
            water_error,
            sensible_error,
        })
    }
    fn breaks(&self, heat: f64) -> Vec<f64> {
        let mut breaks = vec![0.0, heat];
        for state in [
            self.water.law().saturation().liquid(),
            self.water.law().saturation().vapor(),
        ] {
            let q = self.sign
                * self.water_mass
                * (self.water.specific_enthalpy_j_kg() - state.specific_enthalpy().value());
            if q > 0.0 && q < heat {
                breaks.push(q);
            }
        }
        breaks.extend(
            self.sensible
                .anchors
                .iter()
                .map(|a| heat - a.heat)
                .filter(|q| *q > 0.0 && *q < heat),
        );
        breaks.sort_by(f64::total_cmp);
        breaks.dedup();
        breaks
    }
    fn measure(
        &self,
        cx: &Cx<'_>,
        heat: f64,
        tolerance: f64,
        stop_above: Option<f64>,
        config: CountercurrentWaterExchangerConfig,
    ) -> Result<UaMeasure, CountercurrentWaterExchangerError> {
        let breaks = self.breaks(heat);
        if breaks.len() > config.numerics.max_quadrature_panels {
            return Err(ExchangerError::WorkLimit.into());
        }
        let mut visited = 0;
        let result = (|| {
            let mut lower = Dd::ZERO;
            for span in breaks.windows(2) {
                visited += 1;
                let a = self.point(cx, heat, span[0])?;
                let b = self.point(cx, heat, span[1])?;
                if a.gap(self.sign) + a.error() <= 0.0 || b.gap(self.sign) + b.error() <= 0.0 {
                    return Err(ExchangerError::Pinch.into());
                }
                let maximum_gap =
                    self.sign * (a.water - b.sensible) + a.water_error + b.sensible_error;
                positive(maximum_gap, "countercurrent maximum gap")?;
                lower = lower + Dd::from_f64((span[1] - span[0]) / maximum_gap);
            }
            // This rejection uses original monotone ranges, not the later
            // embedded-Simpson error estimate. Incomplete ratings never pass.
            if stop_above
                .is_some_and(|ceiling| lower.to_f64() * (1.0 - 64.0 * f64::EPSILON) > ceiling)
            {
                return Ok(UaMeasure {
                    total: Measure {
                        ua: lower.to_f64(),
                        error: 64.0 * f64::EPSILON * lower.to_f64(),
                        complete: false,
                    },
                    quadrature: 0.0,
                    inverse: 0.0,
                    minimum_gap: 0.0,
                });
            }
            let (mut quadrature, mut inverse) = (Dd::ZERO, Dd::ZERO);
            let mut minimum = f64::INFINITY;
            let total = integrate_panels(
                cx,
                &breaks,
                heat,
                tolerance,
                config.numerics.max_quadrature_panels,
                &mut visited,
                None,
                |a, b| {
                    let panel = self.panel(cx, heat, a, b)?;
                    if panel.total.error <= tolerance * ((b - a) / heat) {
                        quadrature = quadrature + Dd::from_f64(panel.quadrature);
                        inverse = inverse + Dd::from_f64(panel.inverse);
                        minimum = minimum.min(panel.minimum_gap);
                    }
                    Ok::<_, CountercurrentWaterExchangerError>(panel.total)
                },
            )?;
            Ok(UaMeasure {
                total,
                quadrature: quadrature.to_f64(),
                inverse: inverse.to_f64(),
                minimum_gap: minimum,
            })
        })();
        self.panels.set(self.panels.get() + visited);
        result
    }
    fn panel(
        &self,
        cx: &Cx<'_>,
        heat: f64,
        a: f64,
        b: f64,
    ) -> Result<UaMeasure, CountercurrentWaterExchangerError> {
        let width = b - a;
        let samples = [a, a + 0.25 * width, a + 0.5 * width, a + 0.75 * width, b]
            .map(|q| self.point(cx, heat, q));
        let [a, b, c, d, e] = samples;
        let points = [a?, b?, c?, d?, e?];
        let mut minimum = f64::INFINITY;
        for point in points {
            if point.gap(self.sign) + point.error() <= 0.0 {
                return Err(ExchangerError::Pinch.into());
            }
            if point.gap(self.sign) <= point.error() {
                return Err(ExchangerError::Unresolved(
                    "countercurrent gap below source inverse resolution",
                )
                .into());
            }
        }
        for pair in points.windows(2) {
            minimum = minimum.min(
                self.sign * (pair[1].water - pair[0].sensible)
                    - pair[1].water_error
                    - pair[0].sensible_error,
            );
        }
        // Both oriented stream temperatures decrease with q. These four
        // cross-endpoint ranges cover every point of the panel, including an
        // unsampled interior pinch. An unresolved range requests subdivision.
        if minimum <= 0.0 {
            return Ok(UaMeasure {
                total: Measure {
                    ua: 0.0,
                    error: f64::INFINITY,
                    complete: false,
                },
                quadrature: 0.0,
                inverse: 0.0,
                minimum_gap: minimum,
            });
        }
        let reciprocal = points.map(|p| 1.0 / p.gap(self.sign));
        let coarse = (width / 6.0) * (reciprocal[0] + 4.0 * reciprocal[2] + reciprocal[4]);
        let fine = (width / 12.0)
            * (reciprocal[0]
                + 4.0 * reciprocal[1]
                + 2.0 * reciprocal[2]
                + 4.0 * reciprocal[3]
                + reciprocal[4]);
        let quadrature = (fine - coarse).abs() / 15.0;
        let inverse = (width / 12.0)
            * points
                .iter()
                .zip([1.0, 4.0, 2.0, 4.0, 1.0])
                .map(|(p, w)| w * p.error() / (p.gap(self.sign) * minimum))
                .sum::<f64>()
            + 128.0 * f64::EPSILON * fine;
        positive(fine, "countercurrent integrated panel")?;
        if !quadrature.is_finite() || !inverse.is_finite() {
            return Err(ExchangerError::Unresolved("countercurrent panel arithmetic").into());
        }
        Ok(UaMeasure {
            total: Measure {
                ua: fine,
                error: quadrature + inverse,
                complete: true,
            },
            quadrature,
            inverse,
            minimum_gap: minimum,
        })
    }
}
