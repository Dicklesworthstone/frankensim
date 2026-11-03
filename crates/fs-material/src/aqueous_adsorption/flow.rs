//! Continuous ideal mixing with retained, finite-capacity sorbent populations.
//!
//! For one analytical element, `dD/dt = Q*(c_in-D/V)-sum(dA_i/dt)` and
//! `dA_i/dt = k_i*(M_i*q_eq(D/V)-A_i)`. All populations share the actual
//! dissolved inventory. The original Langmuir owner supplies `q_eq`; a
//! DISTINCT local rate declaration supplies k. A batch clock or residence
//! correlation cannot supply that local kinetic authority.
//!
//! Backward Euler reduces each step to a monotone scalar balance for D.
//! Adsorbed and vacant stocks use bounded convex updates. Solution outlet
//! amounts are independently integrated as `h*Q*D_new/V`, never supplied by a
//! balance plug. Full-step/two-half comparison uses divisor ONE: this method
//! is first order, and its local estimates are not certified global errors.
//!
//! The Langmuir form follows USGS PHREEQC example 19, equation 35:
//! <https://water.usgs.gov/water-resources/software/PHREEQC/documentation/phreeqc3-html/phreeqc3-81.htm>.
//! The ideal mixed-tank mole balance follows Fogler, chapter 1:
//! <https://www.pearson.com/content/dam/one-dot-com/one-dot-com/hei-samples/science/9780135337554.pdf>.
//! Neither source supplies the caller's kinetic coefficient or ore calibration.
//! Fixed liquid volume, temperature, dry sorbent masses and chemistry are
//! explicit. This is analytical-element partitioning, not charged-species
//! removal, CIL dissolution, selective carbon motion, regeneration chemistry,
//! pore diffusion, attrition or an inferred plant-recovery forecast.

use super::kinetics::SorbentState;
use super::{
    AdsorbentLoading, AqueousAdsorptionError, AqueousLangmuir, CONCENTRATION_AXIS,
    TEMPERATURE_AXIS, concentration_spec, divide, fractions, multiply, nonnegative, positive,
    temperature_spec,
};
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, MatDbError, QueryPoint};
use fs_qty::{
    Amount, AmountConcentration, Dimensionless, Frequency, Mass, QuantitySpec, Temperature, Time,
    Volume, VolumetricFlowRate,
};
use std::sync::Arc;

/// Original local transfer law; deliberately distinct from the finite-batch law.
pub const LOCAL_LDF_LAW: &str = "local-aqueous-langmuir-linear-driving-force";
/// Constant loading-relaxation rate, in inverse seconds.
pub const RATE_PARAMETER: &str = "loading-relaxation-rate";
/// Actual loaded analytical element per kilogram of retained dry sorbent.
pub const LOADING_AXIS: &str = "adsorbent-loading";
/// Actual fixed, full, well-mixed liquid volume.
pub const VOLUME_AXIS: &str = "solution-volume";
/// Maximum independently declared retained populations.
pub const MAX_POPULATIONS: usize = 128;

/// A failed call leaves the complete accepted inventory and history unchanged.
#[derive(Clone, Debug, PartialEq)]
pub enum AdsorptionFlowError {
    /// Original equilibrium/rate source or typed-coordinate refusal.
    Source(AqueousAdsorptionError),
    /// Invalid source declaration, layout, physical input or numerical setting.
    Invalid(&'static str),
    /// Positive physical amounts or numerical balances lost representable range.
    Arithmetic(&'static str),
    /// The implicit common concentration lies outside original source support.
    OutsideConcentrationSupport,
    /// Requested integration could not complete in the authorized work.
    Budget,
    /// Cancellation was observed before atomic publication.
    Cancelled,
}
impl core::fmt::Display for AdsorptionFlowError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "continuous aqueous adsorption: {self:?}")
    }
}
impl std::error::Error for AdsorptionFlowError {}
impl From<AqueousAdsorptionError> for AdsorptionFlowError {
    fn from(e: AqueousAdsorptionError) -> Self {
        Self::Source(e)
    }
}
impl From<MatDbError> for AdsorptionFlowError {
    fn from(e: MatDbError) -> Self {
        Self::Source(e.into())
    }
}
type Result<T> = core::result::Result<T, AdsorptionFlowError>;

/// Original isotherm and independent local kinetic evidence.
#[derive(Clone, Debug, PartialEq)]
pub struct LocalLangmuirTransfer {
    equilibrium: AqueousLangmuir,
    source: ConstitutiveModelCard,
    rate: Frequency,
}
impl LocalLangmuirTransfer {
    /// Bind one finite nonnegative typed rate to the original isotherm.
    /// Both cards retain their own applicability and additional conditions.
    pub fn from_card(equilibrium: AqueousLangmuir, source: ConstitutiveModelCard) -> Result<Self> {
        source.validate()?;
        require(
            source.law.0 == LOCAL_LDF_LAW
                && source.law_version == 1
                && source.state_schema_version == 0
                && source.initial_state == InitialStatePolicy::ZeroInternalState
                && source.parameters.len() == 1
                && !source.sources.is_empty()
                && !source.provenance.source.trim().is_empty()
                && !source.provenance.license.trim().is_empty(),
            "require original version-one local LDF source with one rate and zero hidden state",
        )?;
        let p = source
            .parameters
            .get(RATE_PARAMETER)
            .ok_or(AdsorptionFlowError::Invalid("missing local LDF rate"))?;
        require(
            p.dims == Frequency::DIMS,
            "local LDF rate needs inverse-second dimensions",
        )?;
        nonnegative(p.value, "local LDF rate must be finite and nonnegative")?;
        for (axis, &(lo, hi)) in source.validity.bounds() {
            require(
                lo.is_finite()
                    && hi.is_finite()
                    && lo <= hi
                    && source.validity.axis_quantities().contains_key(axis)
                    && ![
                        "elapsed-time",
                        "residence-time",
                        "mean-residence-time",
                        "batch-volume",
                    ]
                    .contains(&axis.as_str()),
                "local source needs finite typed conditions, not batch age or residence",
            )?;
        }
        for (axis, quantity, positive_lower) in owned_axes() {
            let (lo, hi) = source
                .validity
                .bound(axis)
                .ok_or(AdsorptionFlowError::Invalid(
                    "missing local kinetic source axis",
                ))?;
            require(
                source.validity.axis_quantities().get(axis) == Some(&quantity)
                    && lo >= 0.0
                    && hi > 0.0
                    && (!positive_lower || lo > 0.0),
                "invalid local kinetic source support",
            )?;
        }
        let rate = Frequency::new(p.value);
        Ok(Self {
            equilibrium,
            source,
            rate,
        })
    }
    /// Original Langmuir source owner, including the analytical-element basis.
    pub fn equilibrium(&self) -> &AqueousLangmuir {
        &self.equilibrium
    }
    /// Original local rate card; no fitted or stronger evidence is manufactured.
    pub fn source(&self) -> &ConstitutiveModelCard {
        &self.source
    }
    /// Declared constant inverse-second rate; zero disables adsorption transfer.
    pub const fn rate(&self) -> Frequency {
        self.rate
    }

    fn admit(&self, point: &QueryPoint) -> Result<()> {
        self.equilibrium.admit(point)?;
        for (axis, &(minimum, maximum)) in self.source.validity.bounds() {
            let value = *point.axes().get(axis).ok_or(AdsorptionFlowError::Invalid(
                "missing original local condition",
            ))?;
            require(
                point.axis_quantities().get(axis)
                    == self.source.validity.axis_quantities().get(axis),
                "local source condition has wrong quantity",
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

/// A retained dry sorbent with explicit loading and original independent sources.
#[derive(Clone, Debug, PartialEq)]
pub struct RetainedSorbent {
    id: String,
    model: LocalLangmuirTransfer,
    mass: Mass,
    initial_adsorbed: Amount,
    capacity: f64,
    conditions: QueryPoint,
}
impl RetainedSorbent {
    /// Declare a fixed positive mass, initial adsorbed stock, and extra source conditions.
    pub fn new(
        id: impl Into<String>,
        model: LocalLangmuirTransfer,
        mass: Mass,
        initial_adsorbed: Amount,
        conditions: QueryPoint,
    ) -> Result<Self> {
        let id = id.into();
        require(
            !id.trim().is_empty() && id.len() <= 256 && !id.chars().any(char::is_control),
            "population ID needs 1..=256 printable bytes",
        )?;
        positive(mass.value(), "retained dry sorbent mass must be positive")?;
        nonnegative(
            initial_adsorbed.value(),
            "initial adsorbed stock must be nonnegative",
        )?;
        let capacity = multiply(
            mass.value(),
            model.equilibrium.capacity().value(),
            "sorbent capacity",
        )?;
        require(
            initial_adsorbed.value() <= capacity,
            "initial sorbent stock exceeds capacity",
        )?;
        require(
            conditions.axes().keys().all(|axis| {
                !owned_axes().iter().any(|(owned, _, _)| *owned == axis)
                    && (model.source.validity.bound(axis).is_some()
                        || model.equilibrium.source().validity.bound(axis).is_some())
            }),
            "only additional original source conditions belong in population context",
        )?;
        Ok(Self {
            id,
            model,
            mass,
            initial_adsorbed,
            capacity,
            conditions,
        })
    }
    /// Explicit caller identity, with no implicit carbon performance meaning.
    pub fn id(&self) -> &str {
        &self.id
    }
    /// Original equilibrium and local-rate declarations.
    pub fn model(&self) -> &LocalLangmuirTransfer {
        &self.model
    }
    /// Fixed retained dry mass [kg].
    pub const fn mass(&self) -> Mass {
        self.mass
    }
    /// Declared initial adsorbed analytical-element stock [mol].
    pub const fn initial_adsorbed(&self) -> Amount {
        self.initial_adsorbed
    }
    /// Original additional conditions; actual concentration/loading are not caller substitutes.
    pub fn conditions(&self) -> &QueryPoint {
        &self.conditions
    }
}

/// Numerical tolerances and hard work bounds; no experimental uncertainty claim.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdsorptionFlowBudget {
    /// Absolute allowance per accepted full-vs-two-half L1 stock/outlet comparison.
    pub absolute_amount_error: Amount,
    /// Additional local allowance as a fraction of interval-available metal, in [0, 0.01].
    pub relative_amount_error: Dimensionless,
    /// Largest attempted step.
    pub maximum_step: Time,
    /// Accepted plus rejected trials, at most 100,000.
    pub maximum_attempts: usize,
    /// Scalar iterations per implicit step, at most 2048.
    pub maximum_root_iterations: usize,
    /// Total scalar/population work including rejected trials, at most ten million.
    pub maximum_work: usize,
}
impl Default for AdsorptionFlowBudget {
    fn default() -> Self {
        Self {
            absolute_amount_error: Amount::new(1e-9),
            relative_amount_error: Dimensionless::new(1e-6),
            maximum_step: Time::new(1.0),
            maximum_attempts: 16_384,
            maximum_root_iterations: 128,
            maximum_work: 1_000_000,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
struct Layout {
    populations: Vec<RetainedSorbent>,
    volume: Volume,
    temperature: Temperature,
    initial_total: f64,
    minimum_concentration: f64,
    maximum_concentration: f64,
}
#[derive(Clone, Debug, PartialEq)]
struct RawState {
    dissolved: f64,
    adsorbed: Vec<f64>,
    vacant: Vec<f64>,
    outlet: f64,
}

/// Persistent inventories in one full ideal tank; sorbent masses remain retained.
#[derive(Clone, Debug, PartialEq)]
pub struct AdsorptionCstr {
    layout: Arc<Layout>,
    raw: RawState,
    elapsed: Time,
    cumulative_inlet: Amount,
    cumulative_outlet: Amount,
    accumulated_error: Amount,
}

/// One completed piecewise-constant flow/feed interval with independent accounts.
#[derive(Clone, Debug)]
pub struct AdsorptionFlowAdvance {
    layout: Arc<Layout>,
    /// Accepted absolute time before this call; it is not a source residence time.
    pub previous_elapsed: Time,
    /// Accepted absolute time after this call.
    pub elapsed: Time,
    /// Actual constant liquid flow during the interval, including zero for a stop.
    pub flow: VolumetricFlowRate,
    /// Actual constant incoming analytical-element concentration.
    pub inlet_concentration: AmountConcentration,
    /// Initial liquid stock for this interval.
    pub initial_dissolved: Amount,
    /// Final liquid stock, still in the tank rather than recovered metal.
    pub dissolved: Amount,
    /// Final liquid concentration [mol/m3].
    pub concentration: AmountConcentration,
    /// Actual integrated incoming metal Q * duration * inlet concentration.
    pub inlet_amount: Amount,
    /// Independently integrated solution outlet for THIS interval.
    pub outlet_amount: Amount,
    /// Cumulative incoming metal over the accepted history.
    pub cumulative_inlet: Amount,
    /// Cumulative solution outlet over the accepted history.
    pub cumulative_outlet: Amount,
    /// Canonical final retained populations, with cumulative signed adsorption.
    pub populations: Vec<SorbentState>,
    /// Signed adsorption during THIS interval; negative is desorption.
    pub transfers: Vec<(String, Amount)>,
    /// Original-source conditions at interval start and finish.
    pub initial_points: Vec<QueryPoint>,
    /// Original-source conditions at the accepted final state.
    pub final_points: Vec<QueryPoint>,
    /// Final stocks + interval outlet - initial stocks - interval inlet [mol].
    pub interval_balance_residual: Amount,
    /// Final stocks + all outlet - original stocks - all inlet [mol].
    pub cumulative_balance_residual: Amount,
    /// Largest accepted local full-vs-two-half L1 stock/outlet difference.
    pub maximum_local_error: Amount,
    /// Sum of accepted local estimates over the history, not a global-error certificate.
    pub accumulated_error_estimate: Amount,
    /// Accepted two-half-step updates.
    pub accepted_steps: usize,
    /// Rejected numerical trials.
    pub rejected_steps: usize,
    /// Scalar/population work charged, including unsuccessful trials.
    pub work: usize,
}
impl AdsorptionFlowAdvance {
    /// Original sorbent masses, isotherms, local rates and declared conditions.
    pub fn models(&self) -> &[RetainedSorbent] {
        &self.layout.populations
    }
}

impl AdsorptionCstr {
    /// Admit 1..=128 populations for the same analytical element and complete initial stocks.
    pub fn new(
        mut populations: Vec<RetainedSorbent>,
        dissolved: Amount,
        volume: Volume,
        temperature: Temperature,
    ) -> Result<Self> {
        require(
            !populations.is_empty() && populations.len() <= MAX_POPULATIONS,
            "require 1..=128 retained sorbent populations",
        )?;
        nonnegative(
            dissolved.value(),
            "initial dissolved stock must be nonnegative",
        )?;
        positive(volume.value(), "full solution volume must be positive")?;
        positive(temperature.value(), "absolute temperature must be positive")?;
        populations.sort_by(|a, b| a.id.cmp(&b.id));
        require(
            populations.windows(2).all(|p| p[0].id != p[1].id),
            "duplicate sorbent identity",
        )?;
        require(
            populations.iter().all(|p| {
                p.model.equilibrium.element() == populations[0].model.equilibrium.element()
            }),
            "sorbents must track the same analytical element",
        )?;
        let initial_total = sum(std::iter::once(dissolved.value())
            .chain(populations.iter().map(|p| p.initial_adsorbed.value())))?;
        let (mut minimum_concentration, mut maximum_concentration) = (0.0_f64, f64::MAX);
        for p in &populations {
            for source in [p.model.equilibrium.source(), &p.model.source] {
                let (low, high) = source.validity.bound(CONCENTRATION_AXIS).ok_or(
                    AdsorptionFlowError::Invalid("missing admitted concentration support"),
                )?;
                minimum_concentration = minimum_concentration.max(low);
                maximum_concentration = maximum_concentration.min(high);
            }
        }
        require(
            minimum_concentration <= maximum_concentration,
            "source concentration boxes do not intersect",
        )?;
        let raw = RawState {
            dissolved: dissolved.value(),
            adsorbed: populations
                .iter()
                .map(|p| p.initial_adsorbed.value())
                .collect(),
            vacant: populations
                .iter()
                .map(|p| p.capacity - p.initial_adsorbed.value())
                .collect(),
            outlet: 0.0,
        };
        let layout = Arc::new(Layout {
            populations,
            volume,
            temperature,
            initial_total,
            minimum_concentration,
            maximum_concentration,
        });
        actual_points(&layout, &raw)?;
        Ok(Self {
            layout,
            raw,
            elapsed: Time::new(0.0),
            cumulative_inlet: Amount::new(0.0),
            cumulative_outlet: Amount::new(0.0),
            accumulated_error: Amount::new(0.0),
        })
    }
    /// Original immutable retained populations in canonical order.
    pub fn models(&self) -> &[RetainedSorbent] {
        &self.layout.populations
    }
    /// Actual accepted dissolved inventory [mol].
    pub fn dissolved(&self) -> Amount {
        Amount::new(self.raw.dissolved)
    }
    /// Actual accepted adsorbed/site inventories.
    pub fn populations(&self) -> Result<Vec<SorbentState>> {
        states(&self.layout, &self.raw)
    }
    /// Accepted absolute integration clock.
    pub const fn elapsed(&self) -> Time {
        self.elapsed
    }
    /// Constant active solution volume.
    pub fn volume(&self) -> Volume {
        self.layout.volume
    }
    /// Constant admitted absolute temperature.
    pub fn temperature(&self) -> Temperature {
        self.layout.temperature
    }
    /// Accepted cumulative incoming metal [mol].
    pub const fn cumulative_inlet(&self) -> Amount {
        self.cumulative_inlet
    }
    /// Accepted cumulative solution export [mol].
    pub const fn cumulative_outlet(&self) -> Amount {
        self.cumulative_outlet
    }

    /// Advance one actual flow/feed interval atomically. A stopped flow leaves
    /// sorption active. Carbon remains retained; no carbon export is implied.
    pub fn advance(
        &mut self,
        flow: VolumetricFlowRate,
        inlet_concentration: AmountConcentration,
        duration: Time,
        budget: AdsorptionFlowBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AdsorptionFlowAdvance> {
        nonnegative(flow.value(), "liquid flow must be finite and nonnegative")?;
        nonnegative(
            inlet_concentration.value(),
            "inlet concentration must be finite and nonnegative",
        )?;
        nonnegative(
            duration.value(),
            "interval duration must be finite and nonnegative",
        )?;
        let incoming = multiply(
            multiply(flow.value(), duration.value(), "incoming volume")?,
            inlet_concentration.value(),
            "incoming metal",
        )?;
        let initial_total = total(&self.raw)?;
        let available = sum([initial_total, incoming].into_iter())?;
        let tolerance = budget.tolerance(available)?;
        let end = sum([self.elapsed.value(), duration.value()].into_iter())?;
        require(
            duration.value() == 0.0 || end > self.elapsed.value(),
            "positive interval loses clock resolution",
        )?;
        let mut work = Work {
            used: 0,
            maximum: budget.maximum_work,
            cancelled: &mut cancelled,
        };
        work.tick()?;
        let initial_points = actual_points(&self.layout, &self.raw)?;
        let initial = self.raw.clone();
        let mut raw = initial.clone();
        raw.outlet = 0.0;
        let mut clock = self.elapsed.value();
        let mut step = budget.maximum_step.value().min(duration.value());
        let (mut accepted, mut rejected, mut largest, mut accumulated) =
            (0, 0, 0.0_f64, self.accumulated_error.value());
        while clock < end {
            work.tick()?;
            if accepted + rejected >= budget.maximum_attempts {
                return Err(AdsorptionFlowError::Budget);
            }
            let remaining = end - clock;
            step = step.min(remaining);
            // Avoid leaving a final remainder whose two substep clocks cannot
            // be represented. Split the remaining interval in two NOW; no
            // time is discarded and the requested maximum step is not enlarged.
            if step < remaining {
                let next = clock + step;
                let tail_middle = next.midpoint(end);
                if next < end && (tail_middle <= next || tail_middle >= end) {
                    step = remaining * 0.5;
                }
            }
            let half = multiply(step, 0.5, "half step")?;
            let midpoint = clock + half;
            let finish = if step == end - clock {
                end
            } else {
                clock + step
            };
            require(
                midpoint > clock && finish > midpoint,
                "substep loses clock resolution",
            )?;
            let mut start = raw.clone();
            start.outlet = 0.0;
            let full = implicit(
                &self.layout,
                &start,
                flow,
                inlet_concentration,
                step,
                budget,
                &mut work,
            )?;
            actual_points(&self.layout, &full)?;
            let first = implicit(
                &self.layout,
                &start,
                flow,
                inlet_concentration,
                half,
                budget,
                &mut work,
            )?;
            actual_points(&self.layout, &first)?;
            let fine = implicit(
                &self.layout,
                &first,
                flow,
                inlet_concentration,
                half,
                budget,
                &mut work,
            )?;
            actual_points(&self.layout, &fine)?;
            let error = sum([
                (fine.dissolved - full.dissolved).abs(),
                (fine.outlet - full.outlet).abs(),
            ]
            .into_iter()
            .chain(
                fine.adsorbed
                    .iter()
                    .zip(&full.adsorbed)
                    .map(|(a, b)| (a - b).abs()),
            ))?;
            if error > tolerance {
                rejected += 1;
                step = half;
                continue;
            }
            let outlet = sum([raw.outlet, fine.outlet].into_iter())?;
            raw = fine;
            raw.outlet = outlet;
            clock = finish;
            accepted += 1;
            largest = largest.max(error);
            accumulated = sum([accumulated, error].into_iter())?;
            if error < tolerance * 0.25 {
                step = (step * 2.0).min(budget.maximum_step.value());
            }
        }
        let final_points = actual_points(&self.layout, &raw)?;
        let cumulative_inlet = sum([self.cumulative_inlet.value(), incoming].into_iter())?;
        let cumulative_outlet = sum([self.cumulative_outlet.value(), raw.outlet].into_iter())?;
        let final_total = total(&raw)?;
        let interval_residual =
            sum([final_total, raw.outlet, -initial_total, -incoming].into_iter())?;
        let cumulative_residual = sum([
            final_total,
            cumulative_outlet,
            -self.layout.initial_total,
            -cumulative_inlet,
        ]
        .into_iter())?;
        let populations = states(&self.layout, &raw)?;
        let mut answer = AdsorptionFlowAdvance {
            layout: self.layout.clone(),
            previous_elapsed: self.elapsed,
            elapsed: Time::new(end),
            flow,
            inlet_concentration,
            initial_dissolved: Amount::new(initial.dissolved),
            dissolved: Amount::new(raw.dissolved),
            concentration: AmountConcentration::new(divide(
                raw.dissolved,
                self.layout.volume.value(),
                "final solution concentration",
            )?),
            inlet_amount: Amount::new(incoming),
            outlet_amount: Amount::new(raw.outlet),
            cumulative_inlet: Amount::new(cumulative_inlet),
            cumulative_outlet: Amount::new(cumulative_outlet),
            transfers: self
                .layout
                .populations
                .iter()
                .enumerate()
                .map(|(i, p)| {
                    (
                        p.id.clone(),
                        Amount::new(raw.adsorbed[i] - initial.adsorbed[i]),
                    )
                })
                .collect(),
            populations,
            initial_points,
            final_points,
            interval_balance_residual: Amount::new(interval_residual),
            cumulative_balance_residual: Amount::new(cumulative_residual),
            maximum_local_error: Amount::new(largest),
            accumulated_error_estimate: Amount::new(accumulated),
            accepted_steps: accepted,
            rejected_steps: rejected,
            work: work.used,
        };
        work.tick()?;
        answer.work = work.used;
        raw.outlet = 0.0;
        self.raw = raw;
        self.elapsed = answer.elapsed;
        self.cumulative_inlet = answer.cumulative_inlet;
        self.cumulative_outlet = answer.cumulative_outlet;
        self.accumulated_error = answer.accumulated_error_estimate;
        Ok(answer)
    }
}

fn implicit(
    layout: &Layout,
    old: &RawState,
    flow: VolumetricFlowRate,
    inlet: AmountConcentration,
    step: f64,
    budget: AdsorptionFlowBudget,
    work: &mut Work<'_, impl FnMut() -> bool>,
) -> Result<RawState> {
    let volume = multiply(flow.value(), step, "substep incoming volume")?;
    let exposure = divide(volume, layout.volume.value(), "flow exposure")?;
    let incoming = multiply(volume, inlet.value(), "substep incoming metal")?;
    let weights = layout
        .populations
        .iter()
        .map(|p| {
            let a = multiply(p.model.rate.value(), step, "local rate exposure")?;
            Ok((
                1.0 / (1.0 + a),
                divide(a, 1.0 + a, "implicit relaxation fraction")?,
            ))
        })
        .collect::<Result<Vec<_>>>()?;
    let mut mobile = vec![old.dissolved, incoming];
    for ((_, beta), adsorbed) in weights.iter().zip(&old.adsorbed) {
        mobile.push(multiply(*beta, *adsorbed, "mobilizable sorbed metal")?);
    }
    let mobile = sum(mobile.into_iter())?;
    if mobile == 0.0 {
        return Ok(old.clone());
    }
    let root_tolerance = multiply(mobile, 64.0 * f64::EPSILON, "implicit balance tolerance")?;
    let maximum = divide(mobile, 1.0 + exposure, "maximum implicit dissolved stock")?;
    let mut low = multiply(
        layout.minimum_concentration,
        layout.volume.value(),
        "lower source amount",
    )?;
    let mut high = if divide(
        maximum,
        layout.volume.value(),
        "maximum implicit concentration",
    )? <= layout.maximum_concentration
    {
        maximum
    } else {
        multiply(
            layout.maximum_concentration,
            layout.volume.value(),
            "upper source amount",
        )?
    };
    if low > high {
        return Err(AdsorptionFlowError::OutsideConcentrationSupport);
    }
    let mut evaluate = |dissolved: f64| -> Result<(RawState, f64)> {
        work.tick()?;
        let c = divide(
            dissolved,
            layout.volume.value(),
            "implicit solution concentration",
        )?;
        let outlet = multiply(exposure, dissolved, "integrated solution outlet")?;
        let mut next = RawState {
            dissolved,
            adsorbed: Vec::with_capacity(weights.len()),
            vacant: Vec::with_capacity(weights.len()),
            outlet: sum([old.outlet, outlet].into_iter())?,
        };
        let mut residual_terms = vec![dissolved, outlet, -mobile];
        for (i, (p, &(alpha, beta))) in layout.populations.iter().zip(&weights).enumerate() {
            work.tick()?;
            if beta == 0.0 {
                next.adsorbed.push(old.adsorbed[i]);
                next.vacant.push(old.vacant[i]);
                continue;
            }
            let (occupied, vacant) = fractions(c, p.model.equilibrium.half_saturation().value())?;
            let target = multiply(p.capacity, occupied, "equilibrium loaded stock")?;
            let free_target = multiply(p.capacity, vacant, "equilibrium vacant stock")?;
            residual_terms.push(multiply(beta, target, "implicit target contribution")?);
            let loaded = interpolate(old.adsorbed[i], target, alpha, beta)?;
            let free = interpolate(old.vacant[i], free_target, alpha, beta)?;
            if loaded > p.capacity
                || free > p.capacity
                || sum([loaded, free, -p.capacity].into_iter())?.abs()
                    > 128.0 * f64::EPSILON * p.capacity
            {
                return Err(AdsorptionFlowError::Arithmetic(
                    "implicit site balance lost range",
                ));
            }
            next.adsorbed.push(loaded);
            next.vacant.push(free);
        }
        Ok((next, sum(residual_terms.into_iter())?))
    };
    // Exact stationary stocks remain a candidate; no arbitrary concentration floor.
    if old.dissolved >= low && old.dissolved <= high && old.dissolved > 0.0 {
        let (same, residual) = evaluate(old.dissolved)?;
        if residual.abs() <= root_tolerance {
            return Ok(same);
        }
    }
    let (left, residual) = evaluate(low)?;
    if residual > root_tolerance {
        return Err(AdsorptionFlowError::OutsideConcentrationSupport);
    }
    if low > 0.0 && residual.abs() <= root_tolerance {
        return Ok(left);
    }
    let (right, residual) = evaluate(high)?;
    if residual < -root_tolerance {
        return Err(AdsorptionFlowError::OutsideConcentrationSupport);
    }
    if high > 0.0 && residual.abs() <= root_tolerance {
        return Ok(right);
    }
    for _ in 0..budget.maximum_root_iterations {
        let middle = low.midpoint(high);
        if middle <= low || middle >= high {
            return Err(AdsorptionFlowError::Arithmetic(
                "implicit root is below stock resolution",
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
    Err(AdsorptionFlowError::Budget)
}

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
    let delta = multiply(weight, difference.abs(), "implicit loading change")?;
    let value = if difference > 0.0 {
        base + delta
    } else {
        base - delta
    };
    if !value.is_finite() || value < old.min(target) || value > old.max(target) {
        return Err(AdsorptionFlowError::Arithmetic(
            "convex loading update lost range",
        ));
    }
    Ok(value)
}
fn actual_points(layout: &Layout, state: &RawState) -> Result<Vec<QueryPoint>> {
    let concentration = divide(
        state.dissolved,
        layout.volume.value(),
        "actual solution concentration",
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
                    divide(*adsorbed, p.mass.value(), "actual sorbent loading")?,
                )?
                .with_quantity(
                    VOLUME_AXIS,
                    QuantitySpec::dimensional(Volume::DIMS),
                    layout.volume.value(),
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
                cumulative_transfer: Amount::new(raw.adsorbed[i] - p.initial_adsorbed.value()),
            })
        })
        .collect()
}
fn total(raw: &RawState) -> Result<f64> {
    sum(std::iter::once(raw.dissolved).chain(raw.adsorbed.iter().copied()))
}
fn owned_axes() -> [(&'static str, QuantitySpec, bool); 4] {
    [
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
impl AdsorptionFlowBudget {
    fn tolerance(self, available: f64) -> Result<f64> {
        require(
            self.absolute_amount_error.value().is_finite()
                && self.absolute_amount_error.value() > 0.0
                && self.relative_amount_error.value().is_finite()
                && (0.0..=0.01).contains(&self.relative_amount_error.value())
                && self.maximum_step.value().is_finite()
                && self.maximum_step.value() > 0.0
                && (1..=100_000).contains(&self.maximum_attempts)
                && (1..=2048).contains(&self.maximum_root_iterations)
                && (1..=10_000_000).contains(&self.maximum_work),
            "invalid adsorption flow numerical budget",
        )?;
        sum([
            self.absolute_amount_error.value(),
            multiply(
                self.relative_amount_error.value(),
                available,
                "relative numerical allowance",
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
            return Err(AdsorptionFlowError::Cancelled);
        }
        if self.used >= self.maximum {
            return Err(AdsorptionFlowError::Budget);
        }
        self.used += 1;
        Ok(())
    }
}
fn require(condition: bool, message: &'static str) -> Result<()> {
    if condition {
        Ok(())
    } else {
        Err(AdsorptionFlowError::Invalid(message))
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
    let value = sum + correction;
    if value.is_finite() {
        Ok(value)
    } else {
        Err(AdsorptionFlowError::Arithmetic(
            "amount accumulation overflow",
        ))
    }
}

#[cfg(test)]
#[path = "flow/tests.rs"]
mod tests;
