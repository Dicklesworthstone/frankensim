//! Steady heat exchange between two finite, independently sourced streams.
//!
//! Both fluids keep their original isobaric Cp(T), source support and entropy
//! integral. There is one signed heat transfer, with equal and opposite stream
//! enthalpy changes. Independent pressure and enthalpy datums are legitimate:
//! the streams exchange heat, never mass. UA is prescribed; there is no storage,
//! phase change, pressure loss, pump work or external heat source in this owner.
//!
//! For heat progress q along the hotter stream, required UA is
//! integral_0^Q dq / (T_hot(q)-T_cold(q)). Counterflow uses cold progress Q-q.
//! Original source knots split this integral. Within each paired affine-Cp
//! interval, squared capacity rates are affine in q; their equality locates
//! every possible interior pinch. Numerical integration has an analytic
//! Simpson truncation bound and an explicit floating-point allowance, not an
//! interval-arithmetic certificate. Both watt and UA gates must pass.
//! Extremely strong exchange that rounds Q onto a thermal pinch is refused;
//! a finite represented inverse UA cannot then be independently audited.

use std::sync::Arc;

use fs_exec::Cx;
use fs_material::thermal::transport::ThermalTransportError;
use fs_material::thermal::transport::isobaric::IsobaricHeatCapacityLaw;
use fs_math::{dd::Dd, det};

mod response;
pub use response::{ExchangerInletResponse, ExchangerResponseConfig};

#[cfg(test)]
mod tests;

/// Relative stream directions through the exchanger.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FlowArrangement {
    CoCurrent,
    Counterflow,
}

/// One independently sourced, prescribed positive mass flow.
#[derive(Clone, Debug, PartialEq)]
pub struct ExchangerStream {
    caloric: Arc<IsobaricHeatCapacityLaw>,
    inlet_temperature_k: f64,
    mass_flow_kg_per_s: f64,
}

impl ExchangerStream {
    pub fn new(
        inlet_temperature_k: f64,
        mass_flow_kg_per_s: f64,
        caloric: Arc<IsobaricHeatCapacityLaw>,
    ) -> Result<Self, ExchangerError> {
        positive(mass_flow_kg_per_s, "mass flow")?;
        let cp = caloric.specific_heat_cp_j_kg_k(inlet_temperature_k)?;
        positive(mass_flow_kg_per_s * cp, "inlet capacity rate")?;
        Ok(Self {
            caloric,
            inlet_temperature_k,
            mass_flow_kg_per_s,
        })
    }

    #[must_use]
    pub fn caloric_law(&self) -> &IsobaricHeatCapacityLaw {
        &self.caloric
    }

    #[must_use]
    pub const fn inlet_temperature_k(&self) -> f64 {
        self.inlet_temperature_k
    }

    #[must_use]
    pub const fn mass_flow_kg_per_s(&self) -> f64 {
        self.mass_flow_kg_per_s
    }
}

/// Explicit physical tolerances and bounded nonlinear/integration work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExchangerConfig {
    pub absolute_heat_tolerance_w: f64,
    pub relative_heat_tolerance: f64,
    pub absolute_conductance_tolerance_w_per_k: f64,
    pub relative_conductance_tolerance: f64,
    pub max_iterations: usize,
    /// Maximum Simpson panels visited per UA evaluation; also bounds each
    /// stream's prepared source-knot count. Total work is bounded by this
    /// limit times max_iterations plus endpoint admission.
    pub max_quadrature_panels: usize,
}

impl Default for ExchangerConfig {
    fn default() -> Self {
        Self {
            absolute_heat_tolerance_w: 1e-7,
            relative_heat_tolerance: 1e-10,
            absolute_conductance_tolerance_w_per_k: 1e-9,
            relative_conductance_tolerance: 1e-10,
            max_iterations: 96,
            max_quadrature_panels: 16_384,
        }
    }
}

/// A fully admitted exchanger result; no partial stream state escapes refusal.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExchangerState {
    pub first_outlet_temperature_k: f64,
    pub second_outlet_temperature_k: f64,
    /// Positive from the first stream to the second; the sole exchanged Q.
    pub heat_rate_first_to_second_w: f64,
    /// Relative heat coordinates, retained independently of rounded absolute T.
    pub first_specific_enthalpy_change_j_kg: f64,
    pub second_specific_enthalpy_change_j_kg: f64,
    /// Native original-source entropy integrals, with the unresolved endpoint
    /// energy residual divided by T_out as a first-order rounding correction.
    pub first_entropy_change_rate_w_per_k: f64,
    pub second_entropy_change_rate_w_per_k: f64,
    pub entropy_generation_w_per_k: f64,
    /// Numerical allowance on the signed entropy balance; it is never clipped.
    pub entropy_roundoff_allowance_w_per_k: f64,
    /// Native source calorimetry at represented outlets minus retained heat.
    /// These expose absolute-temperature rounding instead of hiding it.
    pub first_endpoint_energy_residual_w: f64,
    pub second_endpoint_energy_residual_w: f64,
    /// The endpoint audit admits |residual| <= the configured heat tolerance
    /// plus this explicit allowance for representing absolute temperatures.
    /// This allowance does not relax the independent solved-Q/UA gates.
    pub first_endpoint_energy_roundoff_allowance_w: f64,
    pub second_endpoint_energy_roundoff_allowance_w: f64,
    pub conductance_residual_w_per_k: f64,
    pub conductance_error_allowance_w_per_k: f64,
    pub heat_error_bound_w: f64,
    /// Zero for the exact constant-Cp path or a zero-transfer state.
    pub iterations: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub enum ExchangerError {
    Invalid(&'static str),
    Source(ThermalTransportError),
    /// Finite source support is exhausted before the requested UA is reached.
    OutsideSourceSupport,
    /// A trial crosses or touches an endpoint or interior thermal pinch.
    Pinch,
    Unresolved(&'static str),
    WorkLimit,
    NonConvergence {
        iterations: usize,
    },
    Interrupted,
}

impl From<ThermalTransportError> for ExchangerError {
    fn from(error: ThermalTransportError) -> Self {
        Self::Source(error)
    }
}

impl core::fmt::Display for ExchangerError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source-caloric exchanger: {self:?}")
    }
}
impl std::error::Error for ExchangerError {}

/// Immutable exchanger with independently admitted streams and prescribed UA.
#[derive(Clone, Debug, PartialEq)]
pub struct CaloricHeatExchanger {
    first: ExchangerStream,
    second: ExchangerStream,
    arrangement: FlowArrangement,
    conductance_w_per_k: f64,
}

impl CaloricHeatExchanger {
    pub fn new(
        first: ExchangerStream,
        second: ExchangerStream,
        arrangement: FlowArrangement,
        conductance_w_per_k: f64,
    ) -> Result<Self, ExchangerError> {
        nonnegative(conductance_w_per_k, "conductance")?;
        Ok(Self {
            first,
            second,
            arrangement,
            conductance_w_per_k,
        })
    }

    #[must_use]
    pub const fn streams(&self) -> [&ExchangerStream; 2] {
        [&self.first, &self.second]
    }

    #[must_use]
    pub const fn arrangement(&self) -> FlowArrangement {
        self.arrangement
    }

    #[must_use]
    pub const fn conductance_w_per_k(&self) -> f64 {
        self.conductance_w_per_k
    }

    /// Integrate the original source calorics, then enforce independent source
    /// endpoint energy and entropy audits. A supported but source-limited path
    /// is refused when additional UA would drive either outlet out of support.
    pub fn solve(
        &self,
        cx: &Cx<'_>,
        config: ExchangerConfig,
    ) -> Result<ExchangerState, ExchangerError> {
        poll(cx)?;
        config.validate()?;
        let first_is_hot = self.first.inlet_temperature_k >= self.second.inlet_temperature_k;
        let (hot, cold) = if first_is_hot {
            (&self.first, &self.second)
        } else {
            (&self.second, &self.first)
        };
        let driving = hot.inlet_temperature_k - cold.inlet_temperature_k;
        if self.conductance_w_per_k == 0.0 || driving == 0.0 {
            return self.finish(
                cx,
                0.0,
                self.first.inlet_temperature_k,
                self.second.inlet_temperature_k,
                Measure::default(),
                0.0,
                0,
                config,
            );
        }
        let hot = SourcePath::new(cx, hot, cold.inlet_temperature_k, -1.0, config)?;
        let cold = SourcePath::new(cx, cold, hot.stream.inlet_temperature_k, 1.0, config)?;
        let upper = hot.limit().min(cold.limit());
        positive(upper, "available source-supported heat")?;
        let pair = Pair {
            hot,
            cold,
            arrangement: self.arrangement,
            driving,
        };
        let ua = self.conductance_w_per_k;
        let ua_tolerance = config.ua_tolerance(ua)?;
        let mut lower = 0.0;
        let mut upper_bracket = upper;
        let constant = pair.constant_capacities();
        let mut heat = if let Some((hot_rate, cold_rate)) = constant {
            constant_heat(ua, driving, hot_rate, cold_rate, self.arrangement)?
        } else {
            (ua * driving).min(0.5 * upper)
        };
        if heat > upper {
            return Err(ExchangerError::OutsideSourceSupport);
        }
        positive(heat, "nonzero heat transfer")?;

        // A positive-gap source boundary has finite admissible UA. A thermal
        // pinch instead bounds Q asymptotically and is a valid upper bracket.
        if constant.is_none() && ua > upper / driving {
            let tolerance = (ua_tolerance.min(config.heat_tolerance(upper)? / driving)) / 8.0;
            match pair.measure(cx, upper, tolerance, config, Some(ua)) {
                Ok(measured) if ua > measured.ua + measured.error => {
                    return Err(ExchangerError::OutsideSourceSupport);
                }
                Ok(_) | Err(ExchangerError::Pinch) => {}
                Err(error) => return Err(error),
            }
        }

        for iteration in 0..config.max_iterations {
            poll(cx)?;
            let heat_tolerance = config.heat_tolerance(heat)?;
            let integration_tolerance = ua_tolerance.min(heat_tolerance / driving) / 8.0;
            let measured = match pair.measure(
                cx,
                heat,
                integration_tolerance,
                config,
                Some(ua + ua_tolerance),
            ) {
                Ok(measured) => measured,
                Err(ExchangerError::Pinch) if constant.is_none() => {
                    upper_bracket = heat;
                    heat = lower + 0.5 * (upper_bracket - lower);
                    continue;
                }
                Err(error) => return Err(error),
            };
            let residual = measured.ua - ua;
            // For either arrangement, d(UA_required)/dQ >= 1/driving.
            // This turns the UA residual plus integration allowance into a
            // watt error bound, independently of the bisection bracket width.
            let error = (driving * (residual.abs() + measured.error)).min(upper_bracket - lower);
            if measured.complete
                && residual.abs() + measured.error <= ua_tolerance
                && error <= heat_tolerance
            {
                let hot_out = pair.hot.at(heat)?.temperature;
                let cold_out = pair.cold.at(heat)?.temperature;
                let (signed, first_out, second_out) = if first_is_hot {
                    (heat, hot_out, cold_out)
                } else {
                    (-heat, cold_out, hot_out)
                };
                return self.finish(
                    cx,
                    signed,
                    first_out,
                    second_out,
                    Measure {
                        ua: residual,
                        ..measured
                    },
                    error,
                    if constant.is_some() { 0 } else { iteration + 1 },
                    config,
                );
            }
            if constant.is_some() {
                return Err(ExchangerError::Unresolved("constant-Cp heat/UA gate"));
            }
            if residual > measured.error {
                upper_bracket = heat;
            } else if residual < -measured.error {
                lower = heat;
            } else {
                return Err(ExchangerError::Unresolved(
                    "UA residual below numerical resolution",
                ));
            }
            let next = lower + 0.5 * (upper_bracket - lower);
            if next <= lower || next >= upper_bracket {
                return Err(ExchangerError::Unresolved("heat bracket cannot be refined"));
            }
            heat = next;
        }
        Err(ExchangerError::NonConvergence {
            iterations: config.max_iterations,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn finish(
        &self,
        cx: &Cx<'_>,
        heat: f64,
        first_out: f64,
        second_out: f64,
        measured: Measure,
        heat_error: f64,
        iterations: usize,
        config: ExchangerConfig,
    ) -> Result<ExchangerState, ExchangerError> {
        poll(cx)?;
        let a = endpoint(&self.first, first_out, -heat, config)?;
        let b = endpoint(&self.second, second_out, heat, config)?;
        let entropy = (a.entropy + b.entropy).to_f64();
        let entropy_allowance = a.entropy_allowance + b.entropy_allowance;
        if !entropy.is_finite() || !entropy_allowance.is_finite() || entropy < -entropy_allowance {
            return Err(ExchangerError::Unresolved(
                "original source entropy balance",
            ));
        }
        poll(cx)?;
        Ok(ExchangerState {
            first_outlet_temperature_k: first_out,
            second_outlet_temperature_k: second_out,
            heat_rate_first_to_second_w: heat,
            first_specific_enthalpy_change_j_kg: a.heat_per_mass,
            second_specific_enthalpy_change_j_kg: b.heat_per_mass,
            first_entropy_change_rate_w_per_k: a.entropy.to_f64(),
            second_entropy_change_rate_w_per_k: b.entropy.to_f64(),
            entropy_generation_w_per_k: entropy,
            entropy_roundoff_allowance_w_per_k: entropy_allowance,
            first_endpoint_energy_residual_w: a.energy_residual,
            second_endpoint_energy_residual_w: b.energy_residual,
            first_endpoint_energy_roundoff_allowance_w: a.energy_roundoff,
            second_endpoint_energy_roundoff_allowance_w: b.energy_roundoff,
            conductance_residual_w_per_k: measured.ua,
            conductance_error_allowance_w_per_k: measured.error,
            heat_error_bound_w: heat_error,
            iterations,
        })
    }
}

impl ExchangerConfig {
    fn validate(self) -> Result<(), ExchangerError> {
        positive(self.absolute_heat_tolerance_w, "absolute heat tolerance")?;
        nonnegative(self.relative_heat_tolerance, "relative heat tolerance")?;
        positive(
            self.absolute_conductance_tolerance_w_per_k,
            "absolute UA tolerance",
        )?;
        nonnegative(self.relative_conductance_tolerance, "relative UA tolerance")?;
        if self.relative_heat_tolerance >= 1.0
            || self.relative_conductance_tolerance >= 1.0
            || self.max_iterations == 0
            || self.max_quadrature_panels < 2
        {
            return Err(ExchangerError::Invalid("solver limits"));
        }
        Ok(())
    }
    fn heat_tolerance(self, heat: f64) -> Result<f64, ExchangerError> {
        positive(
            self.absolute_heat_tolerance_w + self.relative_heat_tolerance * heat.abs(),
            "combined heat tolerance",
        )
    }
    fn ua_tolerance(self, ua: f64) -> Result<f64, ExchangerError> {
        positive(
            self.absolute_conductance_tolerance_w_per_k + self.relative_conductance_tolerance * ua,
            "combined UA tolerance",
        )
    }
}

#[derive(Clone, Copy)]
struct Anchor {
    heat: f64,
    temperature: f64,
    cp: f64,
    capacity: f64,
}

/// Heat coordinates at original source knots. Local inversions start from the
/// nearer original endpoint, keeping each query within one source interval.
struct SourcePath<'a> {
    stream: &'a ExchangerStream,
    direction: f64,
    anchors: Vec<Anchor>,
}

#[derive(Clone, Copy)]
struct Point {
    temperature: f64,
    capacity: f64,
}

impl<'a> SourcePath<'a> {
    fn new(
        cx: &Cx<'_>,
        stream: &'a ExchangerStream,
        other_inlet: f64,
        direction: f64,
        config: ExchangerConfig,
    ) -> Result<Self, ExchangerError> {
        let span = stream.caloric.temperature_span_k();
        let boundary = other_inlet.clamp(span[0], span[1]);
        let initial = stream.inlet_temperature_k;
        if boundary == initial {
            return Err(ExchangerError::OutsideSourceSupport);
        }
        let mut temperatures = vec![initial];
        for (temperature, _) in &stream.caloric.source_line().knots {
            poll(cx)?;
            if *temperature > initial.min(boundary) && *temperature < initial.max(boundary) {
                if temperatures.len() >= config.max_quadrature_panels - 1 {
                    return Err(ExchangerError::WorkLimit);
                }
                temperatures.push(*temperature);
            }
        }
        temperatures[1..].sort_by(|a, b| {
            if direction > 0.0 {
                a.total_cmp(b)
            } else {
                b.total_cmp(a)
            }
        });
        temperatures.push(boundary);
        let mut anchors: Vec<Anchor> = Vec::with_capacity(temperatures.len());
        let mut total = Dd::ZERO;
        for temperature in temperatures {
            poll(cx)?;
            let cp = stream.caloric.specific_heat_cp_j_kg_k(temperature)?;
            let capacity = positive(stream.mass_flow_kg_per_s * cp, "source capacity rate")?;
            if let Some(previous) = anchors.last() {
                let increment = direction
                    * stream
                        .caloric
                        .specific_enthalpy_change_j_kg(previous.temperature, temperature)?;
                total = total + Dd::from_f64(stream.mass_flow_kg_per_s) * Dd::from_f64(increment);
                if !total.to_f64().is_finite() || total.to_f64() <= previous.heat {
                    return Err(ExchangerError::Unresolved("source heat breakpoint"));
                }
            }
            anchors.push(Anchor {
                heat: total.to_f64(),
                temperature,
                cp,
                capacity,
            });
        }
        Ok(Self {
            stream,
            direction,
            anchors,
        })
    }

    fn limit(&self) -> f64 {
        self.anchors[self.anchors.len() - 1].heat
    }

    fn at(&self, heat: f64) -> Result<Point, ExchangerError> {
        if !heat.is_finite() || heat < 0.0 || heat > self.limit() {
            return Err(ExchangerError::OutsideSourceSupport);
        }
        let index = self
            .anchors
            .partition_point(|a| a.heat <= heat)
            .saturating_sub(1);
        let a = self.anchors[index];
        if heat == a.heat {
            return Ok(Point {
                temperature: a.temperature,
                capacity: a.capacity,
            });
        }
        let b = self.anchors[index + 1];
        let anchor = if heat - a.heat <= b.heat - heat { a } else { b };
        let change = self.direction * (heat - anchor.heat) / self.stream.mass_flow_kg_per_s;
        let cp = anchor.cp;
        // The native inverse refuses motion smaller than absolute T can
        // represent. Only after that refusal, integrate the half-ULP rounding
        // cell on this original affine-Cp interval. Normal queries use the
        // native inverse alone; weak heat keeps its independent coordinate.
        let temperature = match self
            .stream
            .caloric
            .temperature_after_specific_heat(anchor.temperature, change)
        {
            Ok(temperature) => temperature,
            Err(error) => {
                let neighbor = if change >= 0.0 {
                    anchor.temperature.next_up()
                } else {
                    anchor.temperature.next_down()
                };
                let span = self.stream.caloric.temperature_span_k();
                if neighbor < span[0] || neighbor > span[1] {
                    return Err(error.into());
                }
                let neighbor_cp = self.stream.caloric.specific_heat_cp_j_kg_k(neighbor)?;
                let rounding_heat =
                    (neighbor - anchor.temperature).abs() * (0.375 * cp + 0.125 * neighbor_cp);
                if change.abs() < rounding_heat
                    || (change.abs() == rounding_heat && anchor.temperature.to_bits() & 1 == 0)
                {
                    anchor.temperature
                } else {
                    return Err(error.into());
                }
            }
        };
        let capacity = positive(
            self.stream.mass_flow_kg_per_s
                * self.stream.caloric.specific_heat_cp_j_kg_k(temperature)?,
            "local capacity rate",
        )?;
        Ok(Point {
            temperature,
            capacity,
        })
    }

    fn slope(&self, heat: f64) -> f64 {
        let i = self
            .anchors
            .partition_point(|a| a.heat <= heat)
            .saturating_sub(1)
            .min(self.anchors.len() - 2);
        let a = self.anchors[i];
        let b = self.anchors[i + 1];
        self.stream.mass_flow_kg_per_s * ((b.cp - a.cp) / (b.temperature - a.temperature))
    }

    fn constant_capacity(&self) -> Option<f64> {
        let rate = self.anchors[0].capacity;
        self.anchors
            .iter()
            .all(|a| a.cp == self.anchors[0].cp)
            .then_some(rate)
    }
}

struct Pair<'a> {
    hot: SourcePath<'a>,
    cold: SourcePath<'a>,
    arrangement: FlowArrangement,
    driving: f64,
}

/// Original paired-interval geometry shared by value and response quadrature.
#[derive(Clone, Copy)]
struct PanelGeometry {
    width: f64,
    minimum: f64,
    gaps: [f64; 3],
    derivatives: [f64; 4],
    capacity_min: [f64; 2],
    alpha: [f64; 2],
    capacities: [[f64; 3]; 2],
    temperature_scale: f64,
}

#[derive(Clone, Copy, Debug)]
struct Measure {
    ua: f64,
    error: f64,
    complete: bool,
}

impl Default for Measure {
    fn default() -> Self {
        Self {
            ua: 0.0,
            error: 0.0,
            complete: true,
        }
    }
}

impl Pair<'_> {
    fn constant_capacities(&self) -> Option<(f64, f64)> {
        Some((
            self.hot.constant_capacity()?,
            self.cold.constant_capacity()?,
        ))
    }

    fn cold_progress(&self, heat: f64, progress: f64) -> f64 {
        match self.arrangement {
            FlowArrangement::CoCurrent => progress,
            FlowArrangement::Counterflow => heat - progress,
        }
    }

    fn point(&self, heat: f64, progress: f64) -> Result<(Point, Point), ExchangerError> {
        Ok((
            self.hot.at(progress)?,
            self.cold.at(self.cold_progress(heat, progress))?,
        ))
    }

    fn measure(
        &self,
        cx: &Cx<'_>,
        heat: f64,
        tolerance: f64,
        config: ExchangerConfig,
        stop_above: Option<f64>,
    ) -> Result<Measure, ExchangerError> {
        if let Some((hot, cold)) = self.constant_capacities() {
            let first_gap = match self.arrangement {
                FlowArrangement::CoCurrent => self.driving,
                FlowArrangement::Counterflow => self.driving - heat / cold,
            };
            let last_gap = self.driving
                - heat / hot
                - if self.arrangement == FlowArrangement::CoCurrent {
                    heat / cold
                } else {
                    0.0
                };
            if first_gap <= 0.0 || last_gap <= 0.0 {
                return Err(ExchangerError::Pinch);
            }
            let ratio = (last_gap - first_gap) / first_gap;
            let ua = if ratio == 0.0 {
                heat / first_gap
            } else if ratio.abs() < 0.5 {
                (heat / first_gap) * (ratio.ln_1p() / ratio)
            } else {
                heat * (det::ln(last_gap) - det::ln(first_gap)) / (last_gap - first_gap)
            };
            let scale = self.hot.stream.inlet_temperature_k + self.cold.stream.inlet_temperature_k;
            let error = 64.0 * f64::EPSILON * ua * (1.0 + scale / first_gap.min(last_gap));
            positive(ua, "integrated conductance")?;
            return Ok(Measure {
                ua,
                error,
                complete: true,
            });
        }

        let breaks = self.source_breaks(heat, config)?;
        let mut sum = Dd::ZERO;
        let mut errors = Dd::ZERO;
        let mut visited = 0usize;
        let mut lower_bound = Dd::ZERO;
        // Inspect every original paired interval before refinement: a pinch
        // near the far end must not masquerade as quadrature exhaustion.
        for interval in breaks.windows(2) {
            poll(cx)?;
            self.panel(heat, interval[0], interval[1])?;
            let (h, _) = self.point(heat, interval[0])?;
            let cold_position = if self.arrangement == FlowArrangement::Counterflow {
                interval[1]
            } else {
                interval[0]
            };
            let (_, c) = self.point(heat, cold_position)?;
            let maximum_gap = h.temperature - c.temperature
                + 16.0 * f64::EPSILON * (h.temperature + c.temperature);
            lower_bound = lower_bound + Dd::from_f64((interval[1] - interval[0]) / maximum_gap);
            visited += 1;
        }
        if stop_above
            .is_some_and(|ceiling| lower_bound.to_f64() * (1.0 - 64.0 * f64::EPSILON) > ceiling)
        {
            return Ok(Measure {
                ua: lower_bound.to_f64(),
                error: 64.0 * f64::EPSILON * lower_bound.to_f64(),
                complete: false,
            });
        }
        let mut pending = Vec::new();
        for interval in breaks.windows(2) {
            pending.push((interval[0], interval[1]));
            while let Some((left, right)) = pending.pop() {
                poll(cx)?;
                visited += 1;
                if visited > config.max_quadrature_panels {
                    return Err(ExchangerError::WorkLimit);
                }
                let panel = self.panel(heat, left, right)?;
                let allocation = tolerance * ((right - left) / heat);
                if panel.error <= allocation {
                    sum = sum + Dd::from_f64(panel.ua);
                    errors = errors + Dd::from_f64(panel.error);
                    // A partial integral above the rejection ceiling already
                    // proves that the trial heat is too large. Incomplete
                    // integrals can tighten the bracket but never pass a gate.
                    if stop_above.is_some_and(|ceiling| sum.to_f64() - errors.to_f64() > ceiling) {
                        return Ok(Measure {
                            ua: sum.to_f64(),
                            error: errors.to_f64(),
                            complete: false,
                        });
                    }
                } else {
                    let middle = left + 0.5 * (right - left);
                    if middle <= left || middle >= right {
                        return Err(ExchangerError::Unresolved("UA quadrature precision"));
                    }
                    pending.push((middle, right));
                    pending.push((left, middle));
                }
            }
        }
        let ua = positive(sum.to_f64(), "integrated conductance")?;
        let error = errors.to_f64();
        if !error.is_finite() {
            return Err(ExchangerError::Unresolved("UA quadrature allowance"));
        }
        Ok(Measure {
            ua,
            error,
            complete: true,
        })
    }

    fn source_breaks(
        &self,
        heat: f64,
        config: ExchangerConfig,
    ) -> Result<Vec<f64>, ExchangerError> {
        let mut breaks = vec![0.0, heat];
        breaks.extend(
            self.hot
                .anchors
                .iter()
                .map(|a| a.heat)
                .filter(|&q| q > 0.0 && q < heat),
        );
        breaks.extend(
            self.cold
                .anchors
                .iter()
                .map(|a| self.cold_progress(heat, a.heat))
                .filter(|&q| q > 0.0 && q < heat),
        );
        breaks.sort_by(f64::total_cmp);
        breaks.dedup();
        if breaks.len() > config.max_quadrature_panels {
            return Err(ExchangerError::WorkLimit);
        }
        Ok(breaks)
    }

    #[inline]
    fn panel_geometry(
        &self,
        heat: f64,
        left: f64,
        right: f64,
    ) -> Result<PanelGeometry, ExchangerError> {
        let width = right - left;
        let middle = left + 0.5 * width;
        let (ha, ca) = self.point(heat, left)?;
        let (hb, cb) = self.point(heat, right)?;
        let (hm, cm) = self.point(heat, middle)?;
        let da = ha.temperature - ca.temperature;
        let db = hb.temperature - cb.temperature;
        let dm = hm.temperature - cm.temperature;
        let mut minimum = da.min(db).min(dm);
        if self.arrangement == FlowArrangement::Counterflow {
            // d(T_h-T_c)/dq vanishes iff the positive capacity rates agree.
            // Their squares are affine on this ORIGINAL paired Cp interval.
            let scale = ha
                .capacity
                .max(hb.capacity)
                .max(ca.capacity)
                .max(cb.capacity);
            let fa =
                ((ha.capacity - ca.capacity) / scale) * (ha.capacity / scale + ca.capacity / scale);
            let fb =
                ((hb.capacity - cb.capacity) / scale) * (hb.capacity / scale + cb.capacity / scale);
            if (fa < 0.0 && fb > 0.0) || (fa > 0.0 && fb < 0.0) {
                let critical = left + width * (fa / (fa - fb));
                let (h, c) = self.point(heat, critical)?;
                minimum = minimum.min(h.temperature - c.temperature);
            }
        }
        if minimum <= 0.0 {
            return Err(ExchangerError::Pinch);
        }
        let hot_min = ha.capacity.min(hb.capacity);
        let hot_max = ha.capacity.max(hb.capacity);
        let cold_min = ca.capacity.min(cb.capacity);
        let cold_max = ca.capacity.max(cb.capacity);
        let rh = width / hot_min;
        let rc = width / cold_min;
        let ah = (self.hot.slope(middle).abs() / hot_min) * rh;
        let ac = (self.cold.slope(self.cold_progress(heat, middle)).abs() / cold_min) * rc;
        let d1 = if self.arrangement == FlowArrangement::CoCurrent {
            rh + rc
        } else {
            (-rh + width / cold_max)
                .abs()
                .max((-width / hot_max + rc).abs())
        };
        let v1 = d1 / minimum;
        let v2 = (rh * ah + rc * ac) / minimum;
        let v3 = 3.0 * (rh * ah * ah + rc * ac * ac) / minimum;
        let v4 = 15.0 * (rh * ah * ah * ah + rc * ac * ac * ac) / minimum;
        Ok(PanelGeometry {
            width,
            minimum,
            gaps: [da, dm, db],
            derivatives: [v1, v2, v3, v4],
            capacity_min: [hot_min, cold_min],
            alpha: [ah, ac],
            capacities: [
                [ha.capacity, hm.capacity, hb.capacity],
                [ca.capacity, cm.capacity, cb.capacity],
            ],
            temperature_scale: ha.temperature + hb.temperature + ca.temperature + cb.temperature,
        })
    }

    fn panel(&self, heat: f64, left: f64, right: f64) -> Result<Measure, ExchangerError> {
        let geometry = self.panel_geometry(heat, left, right)?;
        let PanelGeometry {
            width,
            minimum,
            gaps: [da, dm, db],
            derivatives: [v1, v2, v3, v4],
            temperature_scale: scale,
            ..
        } = geometry;
        // Simpson error on the normalized interval u in [0,1]:
        // |E| <= width/2880 * sup |d^4(1/gap)/du^4|.
        let truncation = (width / minimum) / 2880.0
            * (24.0 * v1 * v1 * v1 * v1 + 36.0 * v1 * v1 * v2 + 6.0 * v2 * v2 + 8.0 * v1 * v3 + v4);
        let ua = (width / 6.0) * (1.0 / da + 4.0 / dm + 1.0 / db);
        let roundoff = 64.0 * f64::EPSILON * ua * (1.0 + scale / minimum);
        if !ua.is_finite() || ua <= 0.0 || truncation.is_nan() || roundoff.is_nan() {
            return Err(ExchangerError::Unresolved("UA panel arithmetic"));
        }
        Ok(Measure {
            ua,
            error: truncation + roundoff,
            complete: true,
        })
    }
}

fn constant_heat(
    ua: f64,
    driving: f64,
    hot: f64,
    cold: f64,
    arrangement: FlowArrangement,
) -> Result<f64, ExchangerError> {
    let minimum = hot.min(cold);
    let maximum = hot.max(cold);
    let ratio = minimum / maximum;
    let ntu = ua / minimum;
    let effectiveness = match arrangement {
        FlowArrangement::CoCurrent => -det::expm1(-ntu * (1.0 + ratio)) / (1.0 + ratio),
        FlowArrangement::Counterflow => {
            let difference = (maximum - minimum) / maximum;
            if difference == 0.0 {
                if ntu.is_infinite() {
                    1.0
                } else {
                    ntu / (1.0 + ntu)
                }
            } else {
                let transfer = -det::expm1(-ntu * difference);
                if transfer == 0.0 {
                    ntu / (1.0 + ntu)
                } else {
                    transfer / (difference + ratio * transfer)
                }
            }
        }
    };
    let heat = if ntu == 0.0 {
        ua * driving
    } else {
        (minimum * effectiveness) * driving
    };
    positive(heat, "constant-Cp heat")
}

struct Endpoint {
    heat_per_mass: f64,
    energy_residual: f64,
    energy_roundoff: f64,
    entropy: Dd,
    entropy_allowance: f64,
}

fn endpoint(
    stream: &ExchangerStream,
    outlet: f64,
    heat: f64,
    config: ExchangerConfig,
) -> Result<Endpoint, ExchangerError> {
    let law = &stream.caloric;
    let cp_in = law.specific_heat_cp_j_kg_k(stream.inlet_temperature_k)?;
    let cp_out = law.specific_heat_cp_j_kg_k(outlet)?;
    let native_heat = Dd::from_f64(stream.mass_flow_kg_per_s)
        * Dd::from_f64(law.specific_enthalpy_change_j_kg(stream.inlet_temperature_k, outlet)?);
    let residual = native_heat - Dd::from_f64(heat);
    let energy_residual = residual.to_f64();
    let roundoff = 64.0
        * f64::EPSILON
        * (heat.abs()
            + stream.mass_flow_kg_per_s * cp_in * stream.inlet_temperature_k
            + stream.mass_flow_kg_per_s * cp_out * outlet);
    if !roundoff.is_finite()
        || !energy_residual.is_finite()
        || energy_residual.abs() > config.heat_tolerance(heat)? + roundoff
    {
        return Err(ExchangerError::Unresolved(
            "original source endpoint energy audit",
        ));
    }
    let native_entropy = Dd::from_f64(stream.mass_flow_kg_per_s)
        * Dd::from_f64(
            law.specific_entropy_change_at_temperatures_j_kg_k(stream.inlet_temperature_k, outlet)?,
        );
    let correction = residual / Dd::from_f64(outlet);
    let entropy = native_entropy - correction;
    let entropy_allowance =
        128.0 * f64::EPSILON * (native_entropy.to_f64().abs() + heat.abs() / outlet)
            + correction.to_f64().abs()
                * (energy_residual.abs() / (stream.mass_flow_kg_per_s * cp_out))
                / outlet;
    let heat_per_mass = heat / stream.mass_flow_kg_per_s;
    if !heat_per_mass.is_finite() || (heat != 0.0 && heat_per_mass == 0.0) {
        return Err(ExchangerError::Unresolved("specific heat coordinate"));
    }
    Ok(Endpoint {
        heat_per_mass,
        energy_residual,
        energy_roundoff: roundoff,
        entropy,
        entropy_allowance,
    })
}

fn positive(value: f64, field: &'static str) -> Result<f64, ExchangerError> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(ExchangerError::Invalid(field))
    }
}
fn nonnegative(value: f64, field: &'static str) -> Result<(), ExchangerError> {
    if value.is_finite() && value >= 0.0 {
        Ok(())
    } else {
        Err(ExchangerError::Invalid(field))
    }
}
fn poll(cx: &Cx<'_>) -> Result<(), ExchangerError> {
    cx.checkpoint().map_err(|_| ExchangerError::Interrupted)
}
