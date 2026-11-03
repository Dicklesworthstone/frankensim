//! Reversible finite-time uptake by independent sorbents in one finite solution.
//!
//! Each source declares `dq_i/dt = k_i (q_eq,i(c) - q_i)`. The equilibrium
//! loading is evaluated by the original [`AqueousLangmuir`] owner, at the
//! current shared concentration, not at a precomputed final equilibrium.
//! The linear-driving-force closure and finite liquid balance follow Brandani,
//! *Adsorption* 27 (2021), 353–368, equations (1), (6), and (7):
//! <https://doi.org/10.1007/s10450-020-00258-9>. A constant lumped rate is a
//! separately sourced approximation; neither the isotherm nor a carbon name
//! supplies it. No pore-diffusion mechanism or calibrated CIP/CIL result follows.
//!
//! Backward Euler gives a convex combination of old and equilibrium loadings.
//! One monotone scalar balance determines the shared dissolved stock, retaining
//! site capacities and allowing desorption. Full versus two-half steps estimate
//! local L1 amount error; the positive two-half solution is accepted without
//! Richardson extrapolation or inventory normalization. Reported estimates are
//! not certified global bounds and exclude source/model uncertainty.
//!
//! Volume, dry masses, temperature, and declared additional chemistry remain
//! fixed. All populations track the same analytical element; there is no aqueous
//! speciation, charge balance, co-adsorbate competition, dissolution, flow,
//! attrition, thermal response, or inference of preg-robbing from organic carbon.
//! A failed or cancelled advance leaves all accepted inventories and time intact.

use super::{
    AdsorbentLoading, AqueousAdsorptionError, AqueousLangmuir, CONCENTRATION_AXIS,
    TEMPERATURE_AXIS, concentration_spec, divide, fractions, multiply, nonnegative, positive,
    temperature_spec,
};
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, MatDbError, QueryPoint};
use fs_qty::{
    Amount, AmountConcentration, Frequency, Mass, QuantitySpec, Temperature, Time, Volume,
};
use std::sync::Arc;

/// Version-one constant linear-driving-force rate source.
pub const LDF_LAW: &str = "aqueous-langmuir-linear-driving-force";
/// Nonnegative relaxation rate in inverse seconds, supplied independently of affinity.
pub const RATE_PARAMETER: &str = "loading-relaxation-rate";
/// Absolute batch clock, never reset by a later advance.
pub const ELAPSED_TIME_AXIS: &str = "elapsed-time";
/// Actual adsorbed analytical-element loading in mol/kg dry sorbent.
pub const LOADING_AXIS: &str = "adsorbent-loading";
/// Actual fixed liquid volume, in cubic metres.
pub const VOLUME_AXIS: &str = "solution-volume";
/// Maximum independently declared populations in one bath.
pub const MAX_POPULATIONS: usize = 128;

/// Original source refusal, invalid layout, or bounded numerical failure.
#[derive(Clone, Debug, PartialEq)]
pub enum AdsorptionKineticsError {
    /// Original equilibrium, source-card, or typed-condition admission failure.
    Source(AqueousAdsorptionError),
    /// Invalid kinetic declaration, finite inventory, or numerical settings.
    Invalid(&'static str),
    /// A physical quantity or accepted inventory balance cannot be represented.
    Arithmetic(&'static str),
    /// A common concentration cannot be solved inside every original source box.
    OutsideConcentrationSupport,
    /// Explicit step, root-iteration, or scalar-work budget was exhausted.
    Budget,
    /// Cancellation was observed before publishing an accepted advance.
    Cancelled,
}
impl core::fmt::Display for AdsorptionKineticsError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "aqueous adsorption kinetics: {self:?}")
    }
}
impl std::error::Error for AdsorptionKineticsError {}
impl From<AqueousAdsorptionError> for AdsorptionKineticsError {
    fn from(e: AqueousAdsorptionError) -> Self {
        Self::Source(e)
    }
}
impl From<MatDbError> for AdsorptionKineticsError {
    fn from(e: MatDbError) -> Self {
        Self::Source(e.into())
    }
}
type Result<T> = core::result::Result<T, AdsorptionKineticsError>;

/// Original equilibrium owner plus an independent source-declared kinetic rate.
#[derive(Clone, Debug, PartialEq)]
pub struct LangmuirKinetics {
    equilibrium: AqueousLangmuir,
    source: ConstitutiveModelCard,
    rate: Frequency,
}
impl LangmuirKinetics {
    /// Bind a constant, nonnegative rate to the original equilibrium card.
    /// Both source boxes apply. The kinetic source must explicitly cover time,
    /// instantaneous concentration/loading, fixed volume, and temperature.
    pub fn from_card(equilibrium: AqueousLangmuir, source: ConstitutiveModelCard) -> Result<Self> {
        source.validate()?;
        require(
            source.law.0 == LDF_LAW
                && source.law_version == 1
                && source.state_schema_version == 0
                && source.initial_state == InitialStatePolicy::ZeroInternalState
                && source.parameters.len() == 1
                && !source.sources.is_empty()
                && !source.provenance.source.trim().is_empty()
                && !source.provenance.license.trim().is_empty(),
            "require a sourced version-one constant LDF card with one rate",
        )?;
        let parameter = source
            .parameters
            .get(RATE_PARAMETER)
            .ok_or(AdsorptionKineticsError::Invalid("missing LDF rate"))?;
        require(
            parameter.dims == Frequency::DIMS,
            "LDF rate must have inverse-second dimensions",
        )?;
        nonnegative(parameter.value, "LDF rate must be finite and nonnegative")?;
        for (axis, &(low, high)) in source.validity.bounds() {
            require(
                low.is_finite()
                    && high.is_finite()
                    && low <= high
                    && source.validity.axis_quantities().contains_key(axis),
                "all kinetic source bounds must be finite and typed",
            )?;
        }
        for (axis, quantity, positive_lower) in owned_axes() {
            let (low, high) =
                source
                    .validity
                    .bound(axis)
                    .ok_or(AdsorptionKineticsError::Invalid(
                        "missing owned kinetic source axis",
                    ))?;
            require(
                source.validity.axis_quantities().get(axis) == Some(&quantity)
                    && low >= 0.0
                    && high > 0.0
                    && (!positive_lower || low > 0.0)
                    && (axis != ELAPSED_TIME_AXIS || low == 0.0),
                "invalid owned kinetic source support",
            )?;
        }
        let rate = Frequency::new(parameter.value);
        Ok(Self {
            equilibrium,
            source,
            rate,
        })
    }
    /// Original Langmuir isotherm and analytical-element identity.
    pub fn equilibrium(&self) -> &AqueousLangmuir {
        &self.equilibrium
    }
    /// Original independent rate evidence, parameters and applicability.
    pub fn source(&self) -> &ConstitutiveModelCard {
        &self.source
    }
    /// Original constant relaxation rate; zero explicitly disables transfer.
    pub const fn rate(&self) -> Frequency {
        self.rate
    }

    fn admit(&self, point: &QueryPoint) -> Result<()> {
        self.equilibrium.admit(point)?;
        for (axis, &(minimum, maximum)) in self.source.validity.bounds() {
            let value = *point
                .axes()
                .get(axis)
                .ok_or(AdsorptionKineticsError::Invalid(
                    "missing original kinetic source condition",
                ))?;
            require(
                point.axis_quantities().get(axis)
                    == self.source.validity.axis_quantities().get(axis),
                "kinetic source condition has wrong quantity",
            )?;
            if value < minimum || value > maximum {
                return Err(AqueousAdsorptionError::OutsideValidity {
                    axis: axis.clone(),
                    value,
                    minimum,
                    maximum,
                }
                .into());
            }
        }
        Ok(())
    }
}

/// One explicitly identified sorbent with its original dry mass and initial loading.
#[derive(Clone, Debug, PartialEq)]
pub struct SorbentPopulation {
    id: String,
    model: LangmuirKinetics,
    mass: Mass,
    initially_adsorbed: Amount,
    conditions: QueryPoint,
    capacity: f64,
}
impl SorbentPopulation {
    /// Declare a positive dry mass and a stock between zero and its original capacity.
    /// Conditions contain additional source axes only; actual owned values are supplied by the bath.
    pub fn new(
        id: impl Into<String>,
        model: LangmuirKinetics,
        mass: Mass,
        initially_adsorbed: Amount,
        conditions: QueryPoint,
    ) -> Result<Self> {
        let id = id.into();
        require(
            !id.trim().is_empty() && id.len() <= 256 && !id.chars().any(char::is_control),
            "population ID needs 1..=256 printable bytes",
        )?;
        positive(mass.value(), "dry sorbent mass must be positive")?;
        nonnegative(
            initially_adsorbed.value(),
            "initial adsorbed amount must be nonnegative",
        )?;
        let capacity = multiply(
            mass.value(),
            model.equilibrium.capacity().value(),
            "total sorbent capacity",
        )?;
        require(
            initially_adsorbed.value() <= capacity,
            "initial adsorbed amount exceeds source capacity",
        )?;
        require(
            owned_axes()
                .iter()
                .all(|(axis, _, _)| !conditions.axes().contains_key(*axis)),
            "caller cannot override actual time, concentration, loading, temperature or volume",
        )?;
        Ok(Self {
            id,
            model,
            mass,
            initially_adsorbed,
            conditions,
            capacity,
        })
    }
    /// Stable caller identity; its text does not imply an adsorption rate.
    pub fn id(&self) -> &str {
        &self.id
    }
    /// Original equilibrium and independent kinetic source owners.
    pub fn model(&self) -> &LangmuirKinetics {
        &self.model
    }
    /// Original fixed dry mass.
    pub const fn mass(&self) -> Mass {
        self.mass
    }
    /// Original adsorbed stock, including any preloading.
    pub const fn initially_adsorbed(&self) -> Amount {
        self.initially_adsorbed
    }
    /// Original additional chemical/context coordinates.
    pub fn conditions(&self) -> &QueryPoint {
        &self.conditions
    }
}

/// Accepted analytical-metal stock and independently evolved vacant sites.
#[derive(Clone, Debug, PartialEq)]
pub struct SorbentState {
    /// Original caller population identity.
    pub id: String,
    /// Current adsorbed analytical-element amount.
    pub adsorbed: Amount,
    /// Current loading per original dry sorbent mass.
    pub loading: AdsorbentLoading,
    /// Unoccupied original site capacity, in mol of the same analytical element.
    pub vacant_capacity: Amount,
    /// Signed change from original adsorbed stock; negative is net desorption.
    pub cumulative_transfer: Amount,
}

/// Explicit local error tolerance and bounded integration work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdsorptionNumerics {
    /// Absolute allowance for the full-vs-two-half L1 amount difference per accepted step.
    pub absolute_amount_error: Amount,
    /// Relative allowance multiplied by the original total metal, in (0, 0.01].
    pub relative_amount_error: f64,
    /// Largest attempted step, independent of the persistent absolute clock.
    pub maximum_step: Time,
    /// Total accepted plus rejected full/two-half attempts, at most 100,000.
    pub max_attempts: usize,
    /// Scalar bisection iterations per implicit step, in 1..=2048.
    pub max_root_iterations: usize,
    /// Shared scalar/population work limit, at most ten million.
    pub max_work: usize,
}
impl Default for AdsorptionNumerics {
    fn default() -> Self {
        Self {
            absolute_amount_error: Amount::new(1e-8),
            relative_amount_error: 1e-6,
            maximum_step: Time::new(1.0),
            max_attempts: 16_384,
            max_root_iterations: 128,
            max_work: 1_000_000,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
struct Layout {
    populations: Vec<SorbentPopulation>,
    initial_dissolved: f64,
    total: f64,
    volume: Volume,
    temperature: Temperature,
    minimum_concentration: f64,
    maximum_concentration: f64,
}
#[derive(Clone)]
struct RawState {
    dissolved: f64,
    adsorbed: Vec<f64>,
    vacant: Vec<f64>,
}

/// Persistent private inventories in one closed, homogeneous finite solution bath.
#[derive(Clone, Debug, PartialEq)]
pub struct AdsorptionBatch {
    layout: Arc<Layout>,
    dissolved: Amount,
    populations: Vec<SorbentState>,
    elapsed: Time,
    accumulated_error_estimate: Amount,
}

/// One accepted advance; reported transfers refer to this complete call.
#[derive(Clone, Debug)]
pub struct AdsorptionAdvance {
    layout: Arc<Layout>,
    /// Accepted clock before this call.
    pub previous_elapsed: Time,
    /// Accepted absolute clock after this call.
    pub elapsed: Time,
    /// Dissolved stock at the beginning of this call.
    pub initial_dissolved: Amount,
    /// Final dissolved stock in the common bath.
    pub dissolved: Amount,
    /// Actual shared final solution concentration.
    pub concentration: AmountConcentration,
    /// Canonical final population accounts.
    pub populations: Vec<SorbentState>,
    /// Signed change in each population during this call, in canonical order.
    pub transfers: Vec<(String, Amount)>,
    /// Actual source points at the beginning and end of this call.
    pub initial_points: Vec<QueryPoint>,
    /// Actual final source points with the persistent clock and loading.
    pub final_points: Vec<QueryPoint>,
    /// Independently summed final-minus-original analytical metal, in mol.
    pub conservation_residual: Amount,
    /// Largest accepted local full-vs-two-half L1 difference in this call.
    pub maximum_local_error: Amount,
    /// Sum of accepted local estimates over the complete retained history.
    /// It is not a certified global error bound.
    pub accumulated_error_estimate: Amount,
    /// Accepted two-half-step advances.
    pub accepted_steps: usize,
    /// Rejected trials whose error exceeded the local tolerance.
    pub rejected_steps: usize,
    /// Work charged to this call, including rejected trials and root evaluations.
    pub work: usize,
}
impl AdsorptionAdvance {
    /// Original sorbents, equilibrium cards, independent rate cards and context.
    pub fn models(&self) -> &[SorbentPopulation] {
        &self.layout.populations
    }
}

impl AdsorptionBatch {
    /// Initialize 1..=128 positive-mass populations tracking the same element.
    pub fn new(
        mut populations: Vec<SorbentPopulation>,
        dissolved: Amount,
        volume: Volume,
        temperature: Temperature,
    ) -> Result<Self> {
        require(
            !populations.is_empty() && populations.len() <= MAX_POPULATIONS,
            "require 1..=128 sorbent populations",
        )?;
        nonnegative(dissolved.value(), "dissolved metal must be nonnegative")?;
        positive(volume.value(), "solution volume must be positive")?;
        positive(temperature.value(), "absolute temperature must be positive")?;
        populations.sort_by(|a, b| a.id.cmp(&b.id));
        require(
            populations.windows(2).all(|p| p[0].id != p[1].id),
            "duplicate sorbent population identity",
        )?;
        require(
            populations.iter().all(|p| {
                p.model.equilibrium.element() == populations[0].model.equilibrium.element()
            }),
            "all sorbents must track the same analytical element",
        )?;
        let total = sum(std::iter::once(dissolved.value())
            .chain(populations.iter().map(|p| p.initially_adsorbed.value())))?;
        let mut minimum_concentration: f64 = 0.0;
        let mut maximum_concentration = f64::MAX;
        for p in &populations {
            for source in [p.model.equilibrium.source(), &p.model.source] {
                let (low, high) = source
                    .validity
                    .bound(CONCENTRATION_AXIS)
                    .expect("admitted source axis");
                minimum_concentration = minimum_concentration.max(low);
                maximum_concentration = maximum_concentration.min(high);
            }
        }
        require(
            minimum_concentration <= maximum_concentration,
            "sorbent concentration source boxes do not intersect",
        )?;
        let raw = RawState {
            dissolved: dissolved.value(),
            adsorbed: populations
                .iter()
                .map(|p| p.initially_adsorbed.value())
                .collect(),
            vacant: populations
                .iter()
                .map(|p| p.capacity - p.initially_adsorbed.value())
                .collect(),
        };
        let layout = Arc::new(Layout {
            populations,
            initial_dissolved: dissolved.value(),
            total,
            volume,
            temperature,
            minimum_concentration,
            maximum_concentration,
        });
        actual_points(&layout, &raw, 0.0)?;
        let populations = states(&layout, &raw)?;
        Ok(Self {
            layout,
            dissolved,
            populations,
            elapsed: Time::new(0.0),
            accumulated_error_estimate: Amount::new(0.0),
        })
    }
    /// Original immutable source-bound populations, in canonical ID order.
    pub fn models(&self) -> &[SorbentPopulation] {
        &self.layout.populations
    }
    /// Current shared dissolved analytical-metal stock.
    pub const fn dissolved(&self) -> Amount {
        self.dissolved
    }
    /// Current accepted adsorbed inventories and site accounts.
    pub fn populations(&self) -> &[SorbentState] {
        &self.populations
    }
    /// Absolute accepted clock, retained across calls.
    pub const fn elapsed(&self) -> Time {
        self.elapsed
    }
    /// Fixed initial total analytical metal in both phases.
    pub fn total_metal(&self) -> Amount {
        Amount::new(self.layout.total)
    }

    /// Advance all populations simultaneously; publication is atomic for the whole call.
    /// Every full, half and accepted end point must satisfy both original cards.
    pub fn advance(
        &mut self,
        duration: Time,
        numerics: AdsorptionNumerics,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AdsorptionAdvance> {
        nonnegative(
            duration.value(),
            "advance duration must be finite and nonnegative",
        )?;
        let tolerance = numerics.tolerance(self.layout.total)?;
        let end = self.elapsed.value() + duration.value();
        require(end.is_finite(), "absolute adsorption clock overflows")?;
        if duration.value() > 0.0 && end == self.elapsed.value() {
            return Err(AdsorptionKineticsError::Arithmetic(
                "positive advance is below clock resolution",
            ));
        }
        let mut work = Work {
            used: 0,
            maximum: numerics.max_work,
            cancelled: &mut cancelled,
        };
        work.tick()?;
        let initial = RawState {
            dissolved: self.dissolved.value(),
            adsorbed: self
                .populations
                .iter()
                .map(|p| p.adsorbed.value())
                .collect(),
            vacant: self
                .populations
                .iter()
                .map(|p| p.vacant_capacity.value())
                .collect(),
        };
        let initial_points = actual_points(&self.layout, &initial, self.elapsed.value())?;
        actual_points(&self.layout, &initial, end)?; // Reject unsupported absolute time before integration.
        let mut raw = initial.clone();
        let mut clock = self.elapsed.value();
        let mut step = numerics.maximum_step.value().min(duration.value());
        let (mut accepted, mut rejected, mut largest_error, mut accumulated) =
            (0, 0, 0.0_f64, self.accumulated_error_estimate.value());
        while clock < end {
            work.tick()?;
            if accepted + rejected >= numerics.max_attempts {
                return Err(AdsorptionKineticsError::Budget);
            }
            step = step.min(end - clock);
            let half = multiply(step, 0.5, "half-step duration")?;
            let middle = clock + half;
            let finish = if step == end - clock {
                end
            } else {
                clock + step
            };
            if middle <= clock || finish <= middle {
                return Err(AdsorptionKineticsError::Arithmetic(
                    "substep clock loses positive duration",
                ));
            }
            let full = implicit(&self.layout, &raw, step, numerics, &mut work)?;
            actual_points(&self.layout, &full, finish)?;
            let first = implicit(&self.layout, &raw, half, numerics, &mut work)?;
            actual_points(&self.layout, &first, middle)?;
            let refined = implicit(&self.layout, &first, half, numerics, &mut work)?;
            actual_points(&self.layout, &refined, finish)?;
            let error = sum(
                std::iter::once((refined.dissolved - full.dissolved).abs()).chain(
                    refined
                        .adsorbed
                        .iter()
                        .zip(&full.adsorbed)
                        .map(|(a, b)| (a - b).abs()),
                ),
            )?;
            if error > tolerance {
                rejected += 1;
                step = half;
                continue;
            }
            raw = refined;
            clock = finish;
            accepted += 1;
            largest_error = largest_error.max(error);
            accumulated = sum([accumulated, error].into_iter())?;
            if error < tolerance * 0.25 {
                step = (step * 2.0).min(numerics.maximum_step.value());
            }
        }
        let final_points = actual_points(&self.layout, &raw, end)?;
        let populations = states(&self.layout, &raw)?;
        let residual = balance(&self.layout, &raw)?;
        let mut answer = AdsorptionAdvance {
            layout: self.layout.clone(),
            previous_elapsed: self.elapsed,
            elapsed: Time::new(end),
            initial_dissolved: self.dissolved,
            dissolved: Amount::new(raw.dissolved),
            concentration: AmountConcentration::new(divide(
                raw.dissolved,
                self.layout.volume.value(),
                "final dissolved concentration",
            )?),
            transfers: populations
                .iter()
                .zip(&self.populations)
                .map(|(after, before)| {
                    (
                        after.id.clone(),
                        Amount::new(after.adsorbed.value() - before.adsorbed.value()),
                    )
                })
                .collect(),
            populations: populations.clone(),
            initial_points,
            final_points,
            conservation_residual: Amount::new(residual),
            maximum_local_error: Amount::new(largest_error),
            accumulated_error_estimate: Amount::new(accumulated),
            accepted_steps: accepted,
            rejected_steps: rejected,
            work: work.used,
        };
        work.tick()?;
        answer.work = work.used;
        self.dissolved = answer.dissolved;
        self.populations = populations;
        self.elapsed = answer.elapsed;
        self.accumulated_error_estimate = answer.accumulated_error_estimate;
        Ok(answer)
    }
}

fn implicit(
    layout: &Layout,
    old: &RawState,
    step: f64,
    numerics: AdsorptionNumerics,
    work: &mut Work<'_, impl FnMut() -> bool>,
) -> Result<RawState> {
    let weights = layout
        .populations
        .iter()
        .map(|p| {
            let a = multiply(p.model.rate.value(), step, "rate-times-step")?;
            Ok((
                1.0 / (1.0 + a),
                divide(a, 1.0 + a, "implicit relaxation fraction")?,
            ))
        })
        .collect::<Result<Vec<_>>>()?;
    if weights.iter().all(|(_, beta)| *beta == 0.0) || layout.total == 0.0 {
        return Ok(old.clone());
    }
    let mobile = sum(std::iter::once(Ok(old.dissolved))
        .chain(
            weights
                .iter()
                .zip(&old.adsorbed)
                .map(|((_, beta), n)| multiply(*beta, *n, "mobilizable metal")),
        )
        .collect::<std::result::Result<Vec<_>, _>>()?
        .into_iter())?;
    if mobile == 0.0 {
        return Ok(old.clone());
    }
    let root_tolerance = multiply(mobile, 64.0 * f64::EPSILON, "scalar balance tolerance")?;
    let mut low = multiply(
        layout.minimum_concentration,
        layout.volume.value(),
        "lower source amount",
    )?;
    let total_concentration = divide(
        layout.total,
        layout.volume.value(),
        "total metal concentration",
    )?;
    let mut high = if total_concentration <= layout.maximum_concentration {
        layout.total
    } else {
        multiply(
            layout.maximum_concentration,
            layout.volume.value(),
            "upper source amount",
        )?
    };
    let mut evaluate = |dissolved: f64| -> Result<(RawState, f64)> {
        work.tick()?;
        let c = divide(
            dissolved,
            layout.volume.value(),
            "implicit shared concentration",
        )?;
        let mut next = RawState {
            dissolved,
            adsorbed: Vec::with_capacity(weights.len()),
            vacant: Vec::with_capacity(weights.len()),
        };
        for (i, (p, &(alpha, beta))) in layout.populations.iter().zip(&weights).enumerate() {
            work.tick()?;
            if beta == 0.0 {
                next.adsorbed.push(old.adsorbed[i]);
                next.vacant.push(old.vacant[i]);
                continue;
            }
            // The original owner evaluates the Langmuir occupied/vacant fractions.
            let (occupied, vacant) = fractions(c, p.model.equilibrium.half_saturation().value())?;
            let target = multiply(p.capacity, occupied, "equilibrium target stock")?;
            let free_target = multiply(p.capacity, vacant, "equilibrium vacant sites")?;
            let adsorbed = interpolate(old.adsorbed[i], target, alpha, beta)?;
            let vacant = interpolate(old.vacant[i], free_target, alpha, beta)?;
            if adsorbed < 0.0
                || adsorbed > p.capacity
                || vacant < 0.0
                || vacant > p.capacity
                || sum([adsorbed, vacant, -p.capacity].into_iter())?.abs()
                    > 128.0 * f64::EPSILON * p.capacity
            {
                return Err(AdsorptionKineticsError::Arithmetic(
                    "implicit stock or site balance is unrepresentable",
                ));
            }
            next.adsorbed.push(adsorbed);
            next.vacant.push(vacant);
        }
        let residual = balance(layout, &next)?;
        Ok((next, residual))
    };
    let (left, left_residual) = evaluate(low)?;
    if left_residual > root_tolerance {
        return Err(AdsorptionKineticsError::OutsideConcentrationSupport);
    }
    if low > 0.0 && left_residual.abs() <= root_tolerance {
        return Ok(left);
    }
    let (right, right_residual) = evaluate(high)?;
    if right_residual < -root_tolerance {
        return Err(AdsorptionKineticsError::OutsideConcentrationSupport);
    }
    if high > 0.0 && right_residual.abs() <= root_tolerance {
        return Ok(right);
    }
    for _ in 0..numerics.max_root_iterations {
        let middle = low + (high - low) * 0.5;
        if middle == low || middle == high {
            return Err(AdsorptionKineticsError::Arithmetic(
                "scalar balance is below inventory resolution",
            ));
        }
        let (next, residual) = evaluate(middle)?;
        if residual.abs() <= root_tolerance && middle > 0.0 {
            return Ok(next);
        }
        if residual > 0.0 {
            high = middle;
        } else {
            low = middle;
        }
    }
    Err(AdsorptionKineticsError::Budget)
}

// Evaluate the same convex combination from its nearer endpoint. In particular,
// equal endpoints remain identical; alpha*C + beta*C can round above C even
// when the exact combination is C. No clipping or mass normalization is used.
fn interpolate(old: f64, target: f64, alpha: f64, beta: f64) -> Result<f64> {
    if old == target {
        return Ok(old);
    }
    let (base, other, weight) = if beta <= alpha {
        (old, target, beta)
    } else {
        (target, old, alpha)
    };
    let difference = other - base;
    let increment = multiply(weight, difference.abs(), "implicit loading increment")?;
    let value = if difference > 0.0 {
        base + increment
    } else {
        base - increment
    };
    if !value.is_finite() || value < old.min(target) || value > old.max(target) {
        return Err(AdsorptionKineticsError::Arithmetic(
            "convex loading interpolation loses range",
        ));
    }
    Ok(value)
}

fn actual_points(layout: &Layout, state: &RawState, elapsed: f64) -> Result<Vec<QueryPoint>> {
    let concentration = divide(
        state.dissolved,
        layout.volume.value(),
        "actual dissolved concentration",
    )?;
    layout
        .populations
        .iter()
        .zip(&state.adsorbed)
        .map(|(p, adsorbed)| {
            let point = p
                .conditions
                .clone()
                .with_quantity(
                    TEMPERATURE_AXIS,
                    temperature_spec(),
                    layout.temperature.value(),
                )?
                .with_quantity(CONCENTRATION_AXIS, concentration_spec(), concentration)?
                .with_quantity(
                    LOADING_AXIS,
                    QuantitySpec::dimensional(AdsorbentLoading::DIMS),
                    divide(*adsorbed, p.mass.value(), "actual adsorbent loading")?,
                )?
                .with_quantity(
                    VOLUME_AXIS,
                    QuantitySpec::dimensional(Volume::DIMS),
                    layout.volume.value(),
                )?
                .with_quantity(
                    ELAPSED_TIME_AXIS,
                    QuantitySpec::dimensional(Time::DIMS),
                    elapsed,
                )?;
            p.model.admit(&point)?;
            Ok(point)
        })
        .collect()
}
fn states(layout: &Layout, raw: &RawState) -> Result<Vec<SorbentState>> {
    layout
        .populations
        .iter()
        .enumerate()
        .map(|(i, p)| {
            Ok(SorbentState {
                id: p.id.clone(),
                adsorbed: Amount::new(raw.adsorbed[i]),
                loading: AdsorbentLoading::new(divide(
                    raw.adsorbed[i],
                    p.mass.value(),
                    "accepted sorbent loading",
                )?),
                vacant_capacity: Amount::new(raw.vacant[i]),
                cumulative_transfer: Amount::new(raw.adsorbed[i] - p.initially_adsorbed.value()),
            })
        })
        .collect()
}
fn balance(layout: &Layout, raw: &RawState) -> Result<f64> {
    sum([raw.dissolved, -layout.initial_dissolved]
        .into_iter()
        .chain(
            raw.adsorbed
                .iter()
                .zip(&layout.populations)
                .flat_map(|(n, p)| [*n, -p.initially_adsorbed.value()]),
        ))
}
fn owned_axes() -> [(&'static str, QuantitySpec, bool); 5] {
    [
        (
            ELAPSED_TIME_AXIS,
            QuantitySpec::dimensional(Time::DIMS),
            false,
        ),
        (CONCENTRATION_AXIS, concentration_spec(), false),
        (TEMPERATURE_AXIS, temperature_spec(), true),
        (
            LOADING_AXIS,
            QuantitySpec::dimensional(AdsorbentLoading::DIMS),
            false,
        ),
        (VOLUME_AXIS, QuantitySpec::dimensional(Volume::DIMS), true),
    ]
}
impl AdsorptionNumerics {
    fn tolerance(self, total: f64) -> Result<f64> {
        require(
            self.absolute_amount_error.value().is_finite()
                && self.absolute_amount_error.value() > 0.0
                && self.relative_amount_error.is_finite()
                && self.relative_amount_error > 0.0
                && self.relative_amount_error <= 0.01
                && self.maximum_step.value().is_finite()
                && self.maximum_step.value() > 0.0
                && (1..=100_000).contains(&self.max_attempts)
                && (1..=2048).contains(&self.max_root_iterations)
                && (1..=10_000_000).contains(&self.max_work),
            "invalid adsorption numerical tolerance or work budget",
        )?;
        sum([
            self.absolute_amount_error.value(),
            multiply(
                self.relative_amount_error,
                total,
                "relative amount tolerance",
            )?,
        ]
        .into_iter())
    }
}
struct Work<'a, F: FnMut() -> bool> {
    used: usize,
    maximum: usize,
    cancelled: &'a mut F,
}
impl<F: FnMut() -> bool> Work<'_, F> {
    fn tick(&mut self) -> Result<()> {
        if (self.cancelled)() {
            return Err(AdsorptionKineticsError::Cancelled);
        }
        if self.used >= self.maximum {
            return Err(AdsorptionKineticsError::Budget);
        }
        self.used += 1;
        Ok(())
    }
}
fn require(condition: bool, message: &'static str) -> Result<()> {
    if condition {
        Ok(())
    } else {
        Err(AdsorptionKineticsError::Invalid(message))
    }
}
fn sum(values: impl Iterator<Item = f64>) -> Result<f64> {
    let (mut sum, mut correction) = (0.0_f64, 0.0_f64);
    for value in values {
        let next = sum + value;
        correction += if sum.abs() >= value.abs() {
            (sum - next) + value
        } else {
            (value - next) + sum
        };
        sum = next;
    }
    let result = sum + correction;
    if result.is_finite() {
        Ok(result)
    } else {
        Err(AdsorptionKineticsError::Arithmetic(
            "amount accumulation overflows",
        ))
    }
}

#[cfg(test)]
#[path = "kinetics/tests.rs"]
mod tests;
