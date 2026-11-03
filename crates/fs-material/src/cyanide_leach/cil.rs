//! Simultaneous finite-bath gold dissolution and reversible carbon loading.
//!
//! Every stage uses the original cyanidation, acid/base and optional competing-
//! metal owners, followed by each original local Langmuir transfer owner. An
//! explicitly conserved neutral ion pair removes one Au(CN)2- and one declared
//! monovalent counterion together; desorption returns both. Cyanide already in
//! that pair stays bound and never becomes free cyanide through adsorption.
//!
//! Gold-dissolution progress, dissolved complex, counterion, dissolved oxygen,
//! loaded Au and vacant sites have independent integrated stocks. Product carbon and
//! preg-robbing populations share one actual solution but retain separate
//! sources and reports. The declared role supplies no adsorption parameters.
//!
//! Fixed liquid volume, temperature, dry masses and finite initial dissolved O2
//! are explicit. A sufficient initial counterion reserve is required rather
//! than inventing counterion-limited kinetics. Initial reagent doses are kept;
//! new doses, gas exchange, carbon movement and regeneration are not modeled.
//! RK4 full/two-half differences divided by 15 estimate local inventory error;
//! they are not certified bounds or ore/carbon calibration evidence.

use super::*;
use crate::aqueous_adsorption::{
    AdsorbentLoading,
    flow::{AdsorptionFlowError, LocalLangmuirRate, MAX_POPULATIONS, RetainedSorbent},
};

mod ion_pair;
pub use ion_pair::IonPairTransfer;

/// Actual concentration of the explicitly transferred aqueous counterion [mol/m³].
pub const COUNTERION_AXIS: &str = "ion-pair-counterion-concentration";

/// A refused call publishes neither a partial inventory nor a partial clock.
#[derive(Clone, Debug, PartialEq)]
pub enum CilError {
    /// Original gold, acid/base, metal-complexation or exact-reaction refusal.
    Cyanidation(CyanidationError),
    /// Original isotherm or independent local transfer source refusal.
    Adsorption(AdsorptionFlowError),
    /// Missing or inconsistent physical declaration.
    Invalid(&'static str),
    /// Required inventory, clock, or balance cannot be represented.
    Arithmetic(&'static str),
    /// The authorized integration or nested chemical work was exhausted.
    Budget,
    /// Cancellation was observed before atomic publication.
    Cancelled,
}
impl core::fmt::Display for CilError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "carbon in leach: {self:?}")
    }
}
impl std::error::Error for CilError {}
impl From<CyanidationError> for CilError {
    fn from(value: CyanidationError) -> Self {
        match value {
            CyanidationError::Budget => Self::Budget,
            CyanidationError::Cancelled => Self::Cancelled,
            other => Self::Cyanidation(other),
        }
    }
}
impl From<AdsorptionFlowError> for CilError {
    fn from(value: AdsorptionFlowError) -> Self {
        Self::Adsorption(value)
    }
}
type Result<T> = core::result::Result<T, CilError>;

/// Explicit process destination; this label does not supply material properties.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CarbonRole {
    /// Recoverable product carbon.
    ProductCarbon,
    /// Independently declared naturally occurring preg-robbing sorbent.
    PregRobbing,
}

/// Original retained-sorbent owner and explicit destination of its captured gold.
#[derive(Clone, Debug, PartialEq)]
pub struct CarbonPopulation {
    role: CarbonRole,
    sorbent: RetainedSorbent,
}
impl CarbonPopulation {
    /// Retain one original Au isotherm/rate pair and its explicit dry mass.
    pub fn new(role: CarbonRole, sorbent: RetainedSorbent) -> Result<Self> {
        ensure(
            sorbent.model().equilibrium().element().as_str() == "Au",
            "CIL populations must declare analytical Au",
        )?;
        ensure(
            owned_axes()
                .iter()
                .all(|(axis, _)| !sorbent.conditions().axes().contains_key(*axis))
                && !sorbent.conditions().axes().contains_key(COUNTERION_AXIS),
            "CIL owns actual solution coordinates and counterion concentration",
        )?;
        Ok(Self { role, sorbent })
    }
    /// Explicit recovery destination.
    pub const fn role(&self) -> CarbonRole {
        self.role
    }
    /// Original immutable isotherm, local kinetics, mass and extra conditions.
    pub fn sorbent(&self) -> &RetainedSorbent {
        &self.sorbent
    }
}

/// Source-admitted actual state of one retained population.
#[derive(Clone, Debug, PartialEq)]
pub struct CilSorbentState {
    /// Stable caller population identity.
    pub id: String,
    /// Explicit product or preg-robbing destination.
    pub role: CarbonRole,
    /// Loaded Au, retaining its two cyanide ligands and counterion [mol].
    pub adsorbed: Amount,
    /// Independently integrated available sites [mol].
    pub vacant_capacity: Amount,
    /// Actual Au loading per unit dry mass [mol/kg].
    pub loading: AdsorbentLoading,
    /// Actual signed local derivative and both original source owners.
    pub transfer: LocalLangmuirRate,
}

/// Complete accepted common-solution and retained-sorbent state.
#[derive(Clone, Debug)]
pub struct CilState {
    /// Absolute original batch clock.
    pub elapsed: Time,
    /// Actual dissolved complex, free CN, pH and optional competing-metal solution.
    pub solution: CyanidationState,
    /// Canonically ordered, independently retained populations.
    pub populations: Vec<CilSorbentState>,
    /// Independently integrated aqueous ion-pair counterion [mol].
    pub dissolved_counterion: Amount,
}

/// Independently recomputed final-minus-initial analytical balances [mol].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CilConservation {
    /// Accessible, aqueous and all sorbed Au.
    pub gold: Amount,
    /// Uncomplexed, metal-bound, aqueous-gold and sorbed-gold cyanide ligands.
    pub cyanide: Amount,
    /// Aqueous plus sorbed ion-pair counterion.
    pub counterion: Amount,
    /// Dissolved O2 plus the original gold equation's oxygen consumption.
    pub dissolved_oxygen: Amount,
    /// Largest absolute loaded-plus-vacant site residual.
    pub maximum_site_residual: Amount,
}

/// Local inventory-error allowance and hard runtime limits, not source uncertainty.
#[derive(Clone, Copy, Debug)]
pub struct CilNumerics {
    /// Largest trial interval, strictly positive.
    pub maximum_step: Time,
    /// Cumulative allowance for local L1 inventory-error estimates in this call [mol].
    pub absolute_amount_error: Amount,
    /// Additional allowance as a fraction of all initial tracked Au, in [0, 0.01].
    pub relative_amount_error: Dimensionless,
    /// Accepted plus rejected trials, in 1..=65,536.
    pub maximum_attempts: usize,
    /// Shared work including nested speciation and cancellation, in 1..=10,000,000.
    pub maximum_work: usize,
}
impl Default for CilNumerics {
    fn default() -> Self {
        Self {
            maximum_step: Time::new(1.0),
            absolute_amount_error: Amount::new(1e-10),
            relative_amount_error: Dimensionless::new(1e-7),
            maximum_attempts: 8192,
            maximum_work: 2_000_000,
        }
    }
}

#[derive(Clone, Debug)]
struct CilLayout {
    leach: Arc<Layout>,
    pair: IonPairTransfer,
    populations: Vec<CarbonPopulation>,
    capacities: Vec<f64>,
    total_gold: f64,
    total_cyanide: f64,
    total_counterion: f64,
}

/// Persistent simultaneous CIL state; failures leave its complete history intact.
#[derive(Clone, Debug)]
pub struct CarbonInLeachBatch {
    layout: Arc<CilLayout>,
    state: CilState,
    accumulated_error: f64,
}

/// One atomic accepted advance with separate useful and preg-robbing transfers.
#[derive(Clone, Debug)]
pub struct CilAdvance {
    layout: Arc<CilLayout>,
    pub state: CilState,
    pub interval_gold_dissolved: Amount,
    /// Signed net capture; negative denotes release from product carbon.
    pub interval_product_capture: Amount,
    /// Signed net capture; negative denotes release from preg-robbing material.
    pub interval_preg_robbing_capture: Amount,
    pub conservation: CilConservation,
    pub accepted_steps: usize,
    pub rejected_steps: usize,
    pub work: usize,
    pub maximum_local_error: Amount,
    pub accumulated_error_estimate: Amount,
}
impl CilAdvance {
    /// Original source-backed gold and solution-chemistry owner.
    pub fn model(&self) -> &GoldCyanidation {
        &self.layout.leach.model
    }
    /// Original explicitly conserved neutral-pair transfer declaration.
    pub fn ion_pair(&self) -> &IonPairTransfer {
        &self.layout.pair
    }
    /// Original population sources, fixed masses and recovery destinations.
    pub fn models(&self) -> &[CarbonPopulation] {
        &self.layout.populations
    }
    /// Original accepted dissolved reagent inputs retained at attachment.
    pub fn reagent_inputs(&self) -> dosing::ReagentInputs {
        self.layout.leach.reagents
    }
}

impl CarbonInLeachBatch {
    /// Attach explicit sorbents to a fresh finite-oxygen bath; initial doses survive.
    ///
    /// Elapsed time and gold progress must both be zero. Preloaded populations
    /// bring their own already-bound gold, cyanide and counterion inventories.
    /// The aqueous counterion reserve must cover all initially aqueous and
    /// potentially dissolved gold without introducing a new kinetic limitation.
    pub fn from_batch(
        fresh: GoldCyanidationBatch,
        pair: IonPairTransfer,
        mut populations: Vec<CarbonPopulation>,
    ) -> Result<Self> {
        ensure(
            fresh.elapsed().value() == 0.0 && fresh.state.gold_dissolved.value() == 0.0,
            "attach CIL only to a fresh bath with zero clock and gold progress",
        )?;
        ensure(
            populations.len() <= MAX_POPULATIONS,
            "too many CIL populations",
        )?;
        populations.sort_by(|a, b| a.sorbent.id().cmp(b.sorbent.id()));
        ensure(
            populations
                .windows(2)
                .all(|p| p[0].sorbent.id() != p[1].sorbent.id()),
            "CIL population IDs must be unique",
        )?;
        let leach = fresh.layout;
        for species in [pair.counterion(), pair.adsorbed_pair()] {
            ensure(
                leach
                    .model
                    .chemistry
                    .families()
                    .iter()
                    .all(|family| !family.species().contains(species))
                    && leach.competition.as_ref().is_none_or(|metals| {
                        metals
                            .model
                            .families()
                            .iter()
                            .all(|family| !family.species().contains(species))
                    }),
                "ion-pair species cannot also be owned by a solution family",
            )?;
        }
        ReactionNetwork::new(
            leach
                .model
                .chemistry
                .reaction_network()
                .reactions()
                .iter()
                .cloned()
                .chain([leach.model.reaction.clone(), pair.reaction().clone()])
                .collect(),
        )
        .map_err(CyanidationError::from)?;
        let spectators = dosing::spectators(&leach)?;
        ensure(
            !spectators
                .iter()
                .any(|ion| &ion.species == pair.adsorbed_pair()),
            "the adsorbed neutral pair cannot also be a dissolved spectator",
        )?;
        let counterion = spectators
            .iter()
            .find(|ion| &ion.species == pair.counterion())
            .ok_or(CilError::Invalid(
                "declare the aqueous ion-pair counterion explicitly",
            ))?;
        ensure(
            counterion.charge == 1,
            "ion-pair counterion must have charge +1",
        )?;
        let initial_counterion = counterion.amount.value();
        let required = sum([leach.charge.gold_complex.value(), leach.limit].into_iter())?;
        ensure(
            initial_counterion >= required,
            "insufficient counterion reserve for possible gold loading",
        )?;
        let capacities = populations
            .iter()
            .map(|p| {
                product(
                    p.sorbent.mass().value(),
                    p.sorbent.model().equilibrium().capacity().value(),
                )
            })
            .collect::<super::Result<Vec<_>>>()?;
        let adsorbed: Vec<_> = populations
            .iter()
            .map(|p| p.sorbent.initial_adsorbed().value())
            .collect();
        let preloaded = sum(adsorbed.iter().copied())?;
        let total_gold = sum([
            leach.charge.accessible_gold.value(),
            leach.charge.gold_complex.value(),
            preloaded,
        ]
        .into_iter())?;
        let total_cyanide = sum([
            dosing::cyanide_total(&leach)?,
            product(2.0, leach.charge.gold_complex.value())?,
            product(2.0, preloaded)?,
        ]
        .into_iter())?;
        let total_counterion = sum([initial_counterion, preloaded].into_iter())?;
        let stocks = Stocks {
            gold: 0.0,
            dissolved: leach.charge.gold_complex.value(),
            counterion: initial_counterion,
            oxygen: leach.charge.dissolved_oxygen.value(),
            vacant: capacities
                .iter()
                .zip(&adsorbed)
                .map(|(c, a)| c - a)
                .collect(),
            adsorbed,
        };
        let layout = Arc::new(CilLayout {
            leach,
            pair,
            populations,
            capacities,
            total_gold,
            total_cyanide,
            total_counterion,
        });
        let mut never = || false;
        let mut work = Work {
            used: 0,
            maximum: 10_000_000,
            cancelled: &mut never,
            stopped: None,
        };
        let state = evaluate_cil(&layout, 0.0, &stocks, &mut work)?;
        conservation_cil(&layout, &state)?;
        Ok(Self {
            layout,
            state,
            accumulated_error: fresh.accumulated_error,
        })
    }

    /// Complete accepted snapshot, including all original source points.
    pub fn state(&self) -> &CilState {
        &self.state
    }
    /// Original gold and solution-chemistry owner.
    pub fn model(&self) -> &GoldCyanidation {
        &self.layout.leach.model
    }
    /// Explicit conserved neutral-pair transfer owner.
    pub fn ion_pair(&self) -> &IonPairTransfer {
        &self.layout.pair
    }
    /// Original population source owners in canonical ID order.
    pub fn models(&self) -> &[CarbonPopulation] {
        &self.layout.populations
    }
    /// Original initial charge, before any separately retained reagent inputs.
    pub fn initial_charge(&self) -> &CyanidationCharge {
        &self.layout.leach.charge
    }
    /// Accepted reagent doses preceding attachment.
    pub fn reagent_inputs(&self) -> dosing::ReagentInputs {
        self.layout.leach.reagents
    }

    /// Advance all dissolved, sorbed and vacant stocks together, or publish nothing.
    pub fn advance(
        &mut self,
        duration: Time,
        numerics: CilNumerics,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CilAdvance> {
        nonnegative(
            duration.value(),
            "CIL duration must be finite and nonnegative",
        )?;
        // Reuse the original scalar owner's parameter admission and time-weighted
        // allowance; the coupled local estimate below covers every tracked stock.
        let tolerance = CyanidationNumerics {
            maximum_step: numerics.maximum_step,
            absolute_gold_error: numerics.absolute_amount_error,
            relative_gold_error: numerics.relative_amount_error,
            maximum_attempts: numerics.maximum_attempts,
            maximum_work: numerics.maximum_work,
        }
        .tolerance(self.layout.total_gold)?;
        let end = self.state.elapsed.value() + duration.value();
        positive_or_zero(end, "CIL absolute clock overflows")?;
        ensure(
            duration.value() == 0.0 || end > self.state.elapsed.value(),
            "CIL duration is below clock resolution",
        )?;
        let mut work = Work {
            used: 0,
            maximum: numerics.maximum_work,
            cancelled: &mut cancelled,
            stopped: None,
        };
        work.tick()?;
        let (minimum, maximum) = self
            .layout
            .leach
            .model
            .source
            .validity
            .bound(ELAPSED_TIME_AXIS)
            .expect("original time domain");
        if !(minimum..=maximum).contains(&end) {
            return Err(CyanidationError::OutsideValidity {
                axis: ELAPSED_TIME_AXIS.into(),
                value: end,
                minimum,
                maximum,
            }
            .into());
        }
        let mut state = self.state.clone();
        let mut step = numerics.maximum_step.value().min(duration.value());
        let (mut accepted, mut rejected, mut largest, mut accumulated) =
            (0, 0, 0.0_f64, self.accumulated_error);
        while state.elapsed.value() < end {
            work.tick()?;
            if accepted + rejected >= numerics.maximum_attempts {
                return Err(CilError::Budget);
            }
            let remaining = end - state.elapsed.value();
            step = step.min(remaining);
            let finish = if step == remaining {
                end
            } else {
                state.elapsed.value() + step
            };
            let middle = state.elapsed.value() + 0.5 * (finish - state.elapsed.value());
            if middle <= state.elapsed.value() || finish <= middle {
                return Err(CilError::Arithmetic(
                    "positive CIL substep is below clock resolution",
                ));
            }
            let trial = (|| -> Result<_> {
                let full = rk4_cil(&self.layout, &state, finish, &mut work)?;
                let first = rk4_cil(&self.layout, &state, middle, &mut work)?;
                let refined = rk4_cil(&self.layout, &first, finish, &mut work)?;
                Ok((full, refined))
            })();
            let (full, refined) = match trial {
                Err(CilError::Cyanidation(CyanidationError::TrialDepletion)) => {
                    rejected += 1;
                    step = product(step, 0.5)?;
                    continue;
                }
                other => other?,
            };
            let error = Stocks::from_state(&refined).error(&Stocks::from_state(&full))?;
            let local = product(
                tolerance,
                quotient(finish - state.elapsed.value(), duration.value())?,
            )?;
            if error > local {
                rejected += 1;
                step = product(step, 0.5)?;
                continue;
            }
            conservation_cil(&self.layout, &refined)?;
            state = refined;
            accepted += 1;
            largest = largest.max(error);
            accumulated = sum([accumulated, error].into_iter())?;
            if error < local / 32.0 {
                step = (step * 2.0).min(numerics.maximum_step.value());
            }
        }
        let capture = |role| -> Result<Amount> {
            Ok(Amount::new(sum(state
                .populations
                .iter()
                .zip(&self.state.populations)
                .filter(|(p, _)| p.role == role)
                .map(|(after, before)| {
                    after.adsorbed.value() - before.adsorbed.value()
                }))?))
        };
        let answer = CilAdvance {
            layout: self.layout.clone(),
            interval_gold_dissolved: Amount::new(
                state.solution.gold_dissolved.value() - self.state.solution.gold_dissolved.value(),
            ),
            interval_product_capture: capture(CarbonRole::ProductCarbon)?,
            interval_preg_robbing_capture: capture(CarbonRole::PregRobbing)?,
            conservation: conservation_cil(&self.layout, &state)?,
            state: state.clone(),
            accepted_steps: accepted,
            rejected_steps: rejected,
            work: {
                work.tick()?;
                work.used
            },
            maximum_local_error: Amount::new(largest),
            accumulated_error_estimate: Amount::new(accumulated),
        };
        self.state = state;
        self.accumulated_error = accumulated;
        Ok(answer)
    }
}

#[derive(Clone, Debug)]
struct Stocks {
    gold: f64,
    dissolved: f64,
    counterion: f64,
    oxygen: f64,
    adsorbed: Vec<f64>,
    vacant: Vec<f64>,
}
impl Stocks {
    fn from_state(state: &CilState) -> Self {
        Self {
            gold: state.solution.gold_dissolved.value(),
            dissolved: state.solution.gold_complex.value(),
            counterion: state.dissolved_counterion.value(),
            oxygen: state.solution.dissolved_oxygen.value(),
            adsorbed: state
                .populations
                .iter()
                .map(|p| p.adsorbed.value())
                .collect(),
            vacant: state
                .populations
                .iter()
                .map(|p| p.vacant_capacity.value())
                .collect(),
        }
    }
    fn values(&self) -> impl Iterator<Item = f64> + '_ {
        [self.gold, self.dissolved, self.counterion, self.oxygen]
            .into_iter()
            .chain(self.adsorbed.iter().copied())
            .chain(self.vacant.iter().copied())
    }
    fn error(&self, other: &Self) -> Result<f64> {
        Ok(quotient(
            sum(self
                .values()
                .zip(other.values())
                .map(|(a, b)| (a - b).abs()))?,
            15.0,
        )?)
    }
    fn shifted(&self, rate: &Self, time: f64, layout: &CilLayout) -> Result<Self> {
        Ok(Self {
            gold: progress(self.gold, time, rate.gold, layout.leach.limit)?,
            dissolved: signed_stock(self.dissolved, rate.dissolved, time, f64::MAX)?,
            counterion: signed_stock(self.counterion, rate.counterion, time, f64::MAX)?,
            oxygen: signed_stock(self.oxygen, rate.oxygen, time, f64::MAX)?,
            adsorbed: self
                .adsorbed
                .iter()
                .zip(&rate.adsorbed)
                .zip(&layout.capacities)
                .map(|((&a, &r), &cap)| signed_stock(a, r, time, cap))
                .collect::<Result<_>>()?,
            vacant: self
                .vacant
                .iter()
                .zip(&rate.vacant)
                .zip(&layout.capacities)
                .map(|((&v, &r), &cap)| signed_stock(v, r, time, cap))
                .collect::<Result<_>>()?,
        })
    }
    fn average(rates: [&Self; 4]) -> Result<Self> {
        let component = |values: [f64; 4]| -> Result<f64> {
            Ok(sum(values
                .into_iter()
                .zip([6.0, 3.0, 3.0, 6.0])
                .map(|(r, d)| signed_scale(r, 1.0 / d))
                .collect::<Result<Vec<_>>>()?
                .into_iter())?)
        };
        Ok(Self {
            gold: component(rates.map(|r| r.gold))?,
            dissolved: component(rates.map(|r| r.dissolved))?,
            counterion: component(rates.map(|r| r.counterion))?,
            oxygen: component(rates.map(|r| r.oxygen))?,
            adsorbed: (0..rates[0].adsorbed.len())
                .map(|i| component(rates.map(|r| r.adsorbed[i])))
                .collect::<Result<_>>()?,
            vacant: (0..rates[0].vacant.len())
                .map(|i| component(rates.map(|r| r.vacant[i])))
                .collect::<Result<_>>()?,
        })
    }
}
fn signed_scale(value: f64, factor: f64) -> Result<f64> {
    let magnitude = product(value.abs(), factor)?;
    Ok(if value < 0.0 { -magnitude } else { magnitude })
}
fn signed_stock(initial: f64, rate: f64, time: f64, maximum: f64) -> Result<f64> {
    if rate == 0.0 {
        return Ok(initial);
    }
    let increment = signed_scale(rate, time)?;
    let next = sum([initial, increment].into_iter())?;
    if next < 0.0 || next > maximum || (next == 0.0 && initial > 0.0) {
        return Err(CyanidationError::TrialDepletion.into());
    }
    if next == initial {
        return Err(CilError::Arithmetic(
            "nonzero CIL transfer is below inventory resolution",
        ));
    }
    Ok(next)
}
fn derivatives(layout: &CilLayout, state: &CilState) -> Result<Stocks> {
    let adsorbed = state
        .populations
        .iter()
        .zip(&layout.populations)
        .map(|(state, model)| {
            signed_scale(
                state.transfer.loading_rate.value(),
                model.sorbent.mass().value(),
            )
        })
        .collect::<Result<Vec<_>>>()?;
    let uptake = sum(adsorbed.iter().copied())?;
    let gold = state.solution.rate.value();
    Ok(Stocks {
        gold,
        dissolved: sum([gold, -uptake].into_iter())?,
        counterion: -uptake,
        oxygen: -quotient(gold, 4.0)?,
        vacant: adsorbed.iter().map(|r| -r).collect(),
        adsorbed,
    })
}
fn rk4_cil(
    layout: &CilLayout,
    initial: &CilState,
    end: f64,
    work: &mut Work<'_, impl FnMut() -> bool>,
) -> Result<CilState> {
    let duration = end - initial.elapsed.value();
    let half = product(duration, 0.5)?;
    let middle = initial.elapsed.value() + half;
    let original = Stocks::from_state(initial);
    let first = derivatives(layout, initial)?;
    let second_state = evaluate_cil(
        layout,
        middle,
        &original.shifted(&first, half, layout)?,
        work,
    )?;
    let second = derivatives(layout, &second_state)?;
    let third_state = evaluate_cil(
        layout,
        middle,
        &original.shifted(&second, half, layout)?,
        work,
    )?;
    let third = derivatives(layout, &third_state)?;
    let fourth_state = evaluate_cil(
        layout,
        end,
        &original.shifted(&third, duration, layout)?,
        work,
    )?;
    let fourth = derivatives(layout, &fourth_state)?;
    let average = Stocks::average([&first, &second, &third, &fourth])?;
    evaluate_cil(
        layout,
        end,
        &original.shifted(&average, duration, layout)?,
        work,
    )
}
fn evaluate_cil(
    layout: &CilLayout,
    elapsed: f64,
    stocks: &Stocks,
    work: &mut Work<'_, impl FnMut() -> bool>,
) -> Result<CilState> {
    let mut spectators = dosing::spectators(&layout.leach)?;
    let counterion = spectators
        .iter_mut()
        .find(|ion| &ion.species == layout.pair.counterion())
        .expect("explicit counterion retained");
    counterion.amount = Amount::new(stocks.counterion);
    let solution = evaluate_with_partition(
        &layout.leach,
        elapsed,
        stocks.gold,
        Some(stocks.oxygen),
        Some(stocks.dissolved),
        Some(&spectators),
        work,
    )?;
    let charge = &layout.leach.charge;
    let concentration =
        AmountConcentration::new(quotient(stocks.dissolved, charge.solution_volume.value())?);
    let mut populations = Vec::with_capacity(layout.populations.len());
    for (i, population) in layout.populations.iter().enumerate() {
        work.tick()?;
        let sorbent = &population.sorbent;
        let loading = AdsorbentLoading::new(quotient(stocks.adsorbed[i], sorbent.mass().value())?);
        let mut conditions = sorbent.conditions().clone();
        // Reuse every actual leach coordinate. The adsorption owner separately
        // supplies semantic absolute temperature and dimensional liquid volume.
        for (axis, dims) in owned_axes() {
            if matches!(axis, TEMPERATURE_AXIS | VOLUME_AXIS) {
                continue;
            }
            let value = solution.conditions.axes()[axis];
            conditions = conditions
                .with_quantity(axis, QuantitySpec::dimensional(dims), value)
                .map_err(CyanidationError::from)?;
        }
        conditions = conditions
            .with_quantity(
                COUNTERION_AXIS,
                QuantitySpec::dimensional(AmountConcentration::DIMS),
                quotient(stocks.counterion, charge.solution_volume.value())?,
            )
            .map_err(CyanidationError::from)?;
        let transfer = sorbent.model().rate_at(
            concentration,
            loading,
            charge.solution_volume,
            charge.temperature,
            &conditions,
        )?;
        populations.push(CilSorbentState {
            id: sorbent.id().into(),
            role: population.role,
            adsorbed: Amount::new(stocks.adsorbed[i]),
            vacant_capacity: Amount::new(stocks.vacant[i]),
            loading,
            transfer,
        });
    }
    Ok(CilState {
        elapsed: Time::new(elapsed),
        solution,
        populations,
        dissolved_counterion: Amount::new(stocks.counterion),
    })
}
fn conservation_cil(layout: &CilLayout, state: &CilState) -> Result<CilConservation> {
    let solution = &state.solution;
    let adsorbed = sum(state.populations.iter().map(|p| p.adsorbed.value()))?;
    let gold = sum([
        solution.accessible_gold.value(),
        solution.gold_complex.value(),
        adsorbed,
        -layout.total_gold,
    ]
    .into_iter())?;
    let cyanide = sum([
        solution.uncomplexed_cyanide.value(),
        bound_cyanide(solution),
        product(2.0, solution.gold_complex.value())?,
        product(2.0, adsorbed)?,
        -layout.total_cyanide,
    ]
    .into_iter())?;
    let counterion = sum([
        state.dissolved_counterion.value(),
        adsorbed,
        -layout.total_counterion,
    ]
    .into_iter())?;
    let oxygen = sum([
        solution.dissolved_oxygen.value(),
        quotient(solution.gold_dissolved.value(), 4.0)?,
        -layout.leach.charge.dissolved_oxygen.value(),
    ]
    .into_iter())?;
    let mut maximum_site_residual = 0.0_f64;
    for (p, &capacity) in state.populations.iter().zip(&layout.capacities) {
        let residual = sum([p.adsorbed.value(), p.vacant_capacity.value(), -capacity].into_iter())?;
        check_balance(residual, capacity)?;
        maximum_site_residual = maximum_site_residual.max(residual.abs());
    }
    for (residual, scale) in [
        (gold, layout.total_gold),
        (
            cyanide - ligand_balance_residual(solution),
            layout.total_cyanide,
        ),
        (counterion, layout.total_counterion),
        (oxygen, layout.leach.charge.dissolved_oxygen.value()),
    ] {
        check_balance(residual, scale)?;
    }
    Ok(CilConservation {
        gold: Amount::new(gold),
        cyanide: Amount::new(cyanide),
        counterion: Amount::new(counterion),
        dissolved_oxygen: Amount::new(oxygen),
        maximum_site_residual: Amount::new(maximum_site_residual),
    })
}
fn check_balance(residual: f64, scale: f64) -> Result<()> {
    if residual.abs() > 512.0 * f64::EPSILON * scale {
        Err(CilError::Arithmetic(
            "independently integrated CIL inventory balance does not close",
        ))
    } else {
        Ok(())
    }
}
fn ensure(condition: bool, reason: &'static str) -> Result<()> {
    if condition {
        Ok(())
    } else {
        Err(CilError::Invalid(reason))
    }
}

#[cfg(test)]
mod tests;
