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
    EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState, PhaseStateError,
};
use fs_math::dd::Dd;

pub(crate) mod kernel;
use kernel::{NumericPath, Pair, Segment, Solution, Work};

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
    let pair = Pair::new(&hot.numeric, &cold.numeric);
    let solution = pair.solve(conductance_time_j_k, &mut work)?;
    finish(
        [first, second],
        Some((&hot, &cold, first_hot)),
        conductance_time_j_k,
        solution,
        &mut work,
    )
}

impl FinitePhaseExchangeControls {
    pub(crate) fn validate(self) -> Result<()> {
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

struct Path<'a> {
    body: FinitePhaseBody<'a>,
    direction: f64,
    numeric: NumericPath,
}
impl<'a> Path<'a> {
    fn new(body: FinitePhaseBody<'a>, direction: f64, work: &mut Work<'_>) -> Result<Self> {
        Ok(Self { body, direction, numeric: prepare_path(body, direction, None, work)? })
    }
    fn state(&self, heat: f64) -> Result<EquilibriumPhaseState> {
        if heat == 0.0 { return Ok(self.body.initial); }
        let index = self.numeric.segments.partition_point(|s| s.end_q < heat);
        let segment = self.numeric.segments.get(index)
            .ok_or(FinitePhaseExchangeError::OutsideSourceSupport)?;
        if heat == segment.end_q {
            Ok(self.body.curve.state_at_specific_enthalpy(segment.end_energy.ok_or_else(||
                unresolved("original phase endpoint is unavailable"))?)?)
        } else {
            Ok(self.body.curve.advance_specific_energy(
                self.body.initial, self.direction * heat / self.body.mass_kg,
            )?)
        }
    }
}

/// Prepare the SAME original chart; an optional native-h limit retains a
/// homogeneous conduction layer's narrower admitted support. It does not make
/// a surrogate curve or assign another entropy/phase convention.
pub(crate) fn prepare_path(
    body: FinitePhaseBody<'_>,
    direction: f64,
    admitted_end_h: Option<f64>,
    work: &mut Work<'_>,
) -> Result<NumericPath> {
    let knots = body.curve.knots();
    let initial_h = body.initial.specific_enthalpy_j_kg();
    let native_limit = if let Some(limit_h) = admitted_end_h {
        let last_h = if direction > 0.0 {
            knots[knots.len() - 1].specific_enthalpy_j_kg
        } else { knots[0].specific_enthalpy_j_kg };
        if !limit_h.is_finite() || direction * (limit_h - initial_h) <= 0.0
            || direction * (last_h - limit_h) < 0.0
        {
            return Err(FinitePhaseExchangeError::OutsideSourceSupport);
        }
        Some(positive(body.mass_kg * (limit_h - initial_h).abs(), "native phase heat limit")?)
    } else { None };
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
        let (start, end) = if direction > 0.0 { (lower, upper) } else { (upper, lower) };
        let capacities = if lower.temperature_k == upper.temperature_k {
            [0.0; 2]
        } else {
            let cp = body.curve.sensible_heat_capacities_on_interval(index).unwrap_or(
                [(upper.specific_enthalpy_j_kg - lower.specific_enthalpy_j_kg)
                    / (upper.temperature_k - lower.temperature_k); 2],
            );
            positive(body.mass_kg * cp[0], "lower total capacity")?;
            positive(body.mass_kg * cp[1], "upper total capacity")?;
            if direction > 0.0 { cp } else { [cp[1], cp[0]] }
        };
        let start_q = if segments.is_empty() {
            -body.mass_kg * (start.specific_enthalpy_j_kg - initial_h).abs()
        } else { total.to_f64() };
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
            start_q, end_q, start_t: start.temperature_k, end_t: end.temperature_k,
            end_energy: Some(end.specific_enthalpy_j_kg), capacities, energy_factor: body.mass_kg,
            source_energy_width: (end.specific_enthalpy_j_kg - start.specific_enthalpy_j_kg).abs(),
        });
        // A narrower prepared layer clips Q, never the original inverse's
        // endpoint geometry. The existing full-chart path remains unchanged.
        if native_limit.is_some_and(|limit| end_q >= limit) { break; }
        h = end.specific_enthalpy_j_kg;
        if direction > 0.0 {
            if index + 2 == knots.len() { break; }
            index += 1;
        } else {
            if index == 0 { break; }
            index -= 1;
        }
    }
    Ok(NumericPath { initial_temperature: body.initial.temperature_k(),
        limit: native_limit.map_or(total.to_f64(), |limit| total.to_f64().min(limit)), segments })
}

fn finish(
    bodies: [FinitePhaseBody<'_>; 2],
    pair: Option<(&Path<'_>, &Path<'_>, bool)>,
    exposure: f64,
    solution: Solution,
    work: &mut Work<'_>,
) -> Result<FinitePhaseExchange> {
    work.poll()?;
    let (states, signed_heat) = if let Some((hot_path, cold_path, first_hot)) = pair {
        let hot = hot_path.state(solution.heat)?;
        let cold = cold_path.state(solution.heat)?;
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

pub(crate) fn unresolved(what: &'static str) -> FinitePhaseExchangeError {
    FinitePhaseExchangeError::Unresolved(what)
}
pub(crate) fn positive(value: f64, field: &'static str) -> Result<f64> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(FinitePhaseExchangeError::Invalid(field))
    }
}
pub(crate) fn nonnegative(value: f64, field: &'static str) -> Result<()> {
    if value.is_finite() && value >= 0.0 {
        Ok(())
    } else {
        Err(FinitePhaseExchangeError::Invalid(field))
    }
}
