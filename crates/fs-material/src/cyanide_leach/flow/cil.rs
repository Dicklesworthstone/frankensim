//! Continuous cyanidation with retained product and preg-robbing sorbents.
//!
//! The existing cyanidation CSTR owns liquid/slurry residence and local leach
//! chemistry. Explicit dry sorbent populations stay in this tank while the
//! conserved Au(CN)2-/counterion pair loads or desorbs through the existing
//! Langmuir owners. A symmetric carbon-transfer / liquid-leach-flow /
//! carbon-transfer split keeps the two original runtime owners separate while
//! full-vs-two-half trials control their coupling error.
//!
//! Carbon does not move between tanks in this owner. Screen losses, carbon
//! residence distributions and counter-current carbon transfer require a
//! separate solid-transport balance.

use super::*;
use crate::{
    aqueous_adsorption::{
        AdsorbentLoading,
        flow::{AdsorptionFlowError, MAX_POPULATIONS},
    },
    cyanide_leach::cil::{
        COUNTERION_AXIS, CarbonPopulation, CarbonRole, CilSorbentState, IonPairTransfer,
    },
};

#[derive(Clone, Debug, PartialEq)]
pub enum FlowCilError {
    Cyanidation(CyanidationError),
    Adsorption(AdsorptionFlowError),
    Invalid(&'static str),
    Arithmetic(&'static str),
    Budget,
    Cancelled,
}
impl core::fmt::Display for FlowCilError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "continuous carbon-in-leach: {self:?}")
    }
}
impl std::error::Error for FlowCilError {}
impl From<CyanidationError> for FlowCilError {
    fn from(value: CyanidationError) -> Self {
        match value {
            CyanidationError::Budget => Self::Budget,
            CyanidationError::Cancelled => Self::Cancelled,
            other => Self::Cyanidation(other),
        }
    }
}
impl From<AdsorptionFlowError> for FlowCilError {
    fn from(value: AdsorptionFlowError) -> Self {
        Self::Adsorption(value)
    }
}
type CilFlowResult<T> = core::result::Result<T, FlowCilError>;

#[derive(Clone, Debug)]
pub struct FlowCilState {
    pub elapsed: Time,
    pub solution: CyanidationFlowState,
    pub populations: Vec<CilSorbentState>,
    pub dissolved_counterion: Amount,
}

#[derive(Clone, Debug, PartialEq)]
pub struct FlowCilConservation {
    pub gold: Amount,
    pub cyanide: Amount,
    pub counterion: Amount,
    pub dissolved_oxygen: Amount,
    pub maximum_site_residual: Amount,
}

#[derive(Clone, Debug)]
pub struct CarbonInLeachCstr {
    tank: CyanidationCstr,
    pair: Arc<IonPairTransfer>,
    populations: Arc<Vec<CarbonPopulation>>,
    capacities: Arc<Vec<f64>>,
    initial_adsorbed: Arc<Vec<f64>>,
    adsorbed: Vec<f64>,
    vacant: Vec<f64>,
    counterion_index: usize,
    state: FlowCilState,
    coupling_error: f64,
}

#[derive(Clone, Debug)]
pub struct FlowCilAdvance {
    pair: Arc<IonPairTransfer>,
    populations: Arc<Vec<CarbonPopulation>>,
    pub previous_elapsed: Time,
    pub elapsed: Time,
    pub flow: VolumetricFlowRate,
    pub feed: CyanidationFeed,
    pub state: FlowCilState,
    pub inlet: CyanidationInventory,
    pub outlet: CyanidationInventory,
    pub cumulative_inlet: CyanidationInventory,
    pub cumulative_outlet: CyanidationInventory,
    pub interval_gold_dissolved: Amount,
    pub cumulative_gold_dissolved: Amount,
    pub interval_product_capture: Amount,
    pub interval_preg_robbing_capture: Amount,
    pub interval_conservation: FlowCilConservation,
    pub cumulative_conservation: FlowCilConservation,
    pub maximum_coupling_error: Amount,
    pub accumulated_coupling_error: Amount,
    pub liquid_error_estimate: Amount,
    pub accepted_steps: usize,
    pub rejected_steps: usize,
}
impl FlowCilAdvance {
    pub fn ion_pair(&self) -> &IonPairTransfer {
        &self.pair
    }
    pub fn models(&self) -> &[CarbonPopulation] {
        self.populations.as_slice()
    }
}

#[derive(Clone)]
struct CoupledTrial {
    tank: CyanidationCstr,
    adsorbed: Vec<f64>,
    vacant: Vec<f64>,
}
impl CoupledTrial {
    fn values(&self) -> impl Iterator<Item = f64> + '_ {
        self.tank
            .raw
            .iter()
            .chain(&self.tank.inlet)
            .chain(&self.tank.outlet)
            .copied()
            .chain(std::iter::once(self.tank.gold_dissolved))
            .chain(self.adsorbed.iter().copied())
            .chain(self.vacant.iter().copied())
    }
    fn error(&self, other: &Self) -> CilFlowResult<f64> {
        Ok(sum(
            self.values()
                .zip(other.values())
                .map(|(a, b)| (a - b).abs()),
        )?)
    }
}

impl CarbonInLeachCstr {
    pub fn from_tank(
        tank: CyanidationCstr,
        pair: IonPairTransfer,
        mut populations: Vec<CarbonPopulation>,
    ) -> CilFlowResult<Self> {
        require_local(
            tank.elapsed().value() == 0.0
                && tank.gold_dissolved == 0.0
                && tank.inlet.iter().all(|n| *n == 0.0)
                && tank.outlet.iter().all(|n| *n == 0.0),
            "attach retained carbon only to a fresh continuous leach tank",
        )?;
        require_local(
            !populations.is_empty() && populations.len() <= MAX_POPULATIONS,
            "require 1..=128 retained CIL populations",
        )?;
        populations.sort_by(|a, b| a.sorbent().id().cmp(b.sorbent().id()));
        require_local(
            populations
                .windows(2)
                .all(|p| p[0].sorbent().id() != p[1].sorbent().id()),
            "continuous CIL population IDs must be unique",
        )?;
        let layout = &tank.layout;
        let model = &layout.model;
        for species in [pair.counterion(), pair.adsorbed_pair()] {
            require_local(
                model.chemistry.families().iter().all(|family| {
                    family.species().iter().all(|candidate| candidate != species)
                }) && layout.competition.as_ref().is_none_or(|metals| {
                    metals.model.families().iter().all(|family| {
                        family.species().iter().all(|candidate| candidate != species)
                    })
                }),
                "ion-pair species cannot also be owned by a solution family",
            )?;
        }
        ReactionNetwork::new(
            layout
                .competition
                .as_ref()
                .map_or(model.chemistry.reaction_network(), |m| m.model.reaction_network())
                .reactions()
                .iter()
                .cloned()
                .chain([model.reaction.clone(), pair.reaction().clone()])
                .collect(),
        )
        .map_err(CyanidationError::from)?;

        let spectator_offset = 4 + layout.charge.other_families.len();
        let spectator_position = layout
            .charge
            .spectators
            .iter()
            .position(|ion| &ion.species == pair.counterion())
            .ok_or(FlowCilError::Invalid(
                "declare the ion-pair counterion as a flowing spectator",
            ))?;
        require_local(
            layout.charge.spectators[spectator_position].charge == 1,
            "ion-pair counterion must have charge +1",
        )?;
        require_local(
            layout
                .charge
                .spectators
                .iter()
                .all(|ion| &ion.species != pair.adsorbed_pair()),
            "adsorbed neutral pair cannot also be a dissolved spectator",
        )?;
        let counterion_index = spectator_offset + spectator_position;
        require_local(
            tank.raw[counterion_index] >= sum([tank.raw[0], tank.raw[3]].into_iter())?,
            "initial counterion must cover all accessible and aqueous-complex Au",
        )?;

        let capacities = populations
            .iter()
            .map(|p| {
                product(
                    p.sorbent().mass().value(),
                    p.sorbent().model().equilibrium().capacity().value(),
                )
            })
            .collect::<core::result::Result<Vec<_>, _>>()?;
        let adsorbed = populations
            .iter()
            .map(|p| p.sorbent().initial_adsorbed().value())
            .collect::<Vec<_>>();
        let vacant = capacities
            .iter()
            .zip(&adsorbed)
            .map(|(capacity, loaded)| capacity - loaded)
            .collect::<Vec<_>>();
        let pair = Arc::new(pair);
        let populations = Arc::new(populations);
        let capacities = Arc::new(capacities);
        let initial_adsorbed = Arc::new(adsorbed.clone());
        let accumulated = sum([self.coupling_error, call_error].into_iter())?;
        let state = admit_state(
            &tank,
            &pair,
            &populations,
            &capacities,
            counterion_index,
            &adsorbed,
            &vacant,
            2_000_000,
            &mut || false,
        )?;
        Ok(Self {
            tank,
            pair,
            populations,
            capacities,
            initial_adsorbed,
            adsorbed,
            vacant,
            counterion_index,
            state,
            coupling_error: 0.0,
        })
    }

    pub fn state(&self) -> &FlowCilState {
        &self.state
    }
    pub fn tank(&self) -> &CyanidationCstr {
        &self.tank
    }
    pub fn ion_pair(&self) -> &IonPairTransfer {
        &self.pair
    }
    pub fn models(&self) -> &[CarbonPopulation] {
        self.populations.as_slice()
    }

    pub fn advance(
        &mut self,
        duration: Time,
        flow: VolumetricFlowRate,
        feed: &CyanidationFeed,
        numerics: FlowNumerics,
        mut cancelled: impl FnMut() -> bool,
    ) -> CilFlowResult<FlowCilAdvance> {
        nonnegative(duration.value(), "continuous CIL duration must be nonnegative")?;
        nonnegative(flow.value(), "continuous CIL flow must be nonnegative")?;
        if cancelled() {
            return Err(FlowCilError::Cancelled);
        }

        let counterion_feed = feed
            .spectators
            .iter()
            .find(|ion| &ion.species == self.pair.counterion())
            .ok_or(FlowCilError::Invalid("missing ion-pair counterion feed"))?;
        require_local(
            counterion_feed.charge == 1
                && counterion_feed.concentration.value()
                    >= sum(
                        [feed.accessible_gold.value(), feed.gold_complex.value()].into_iter(),
                    )?,
            "counterion feed must cover all incoming accessible and complexed Au",
        )?;

        // A zero-duration clone validates the complete feed/source surface with
        // the original CSTR owner and cannot publish into this wrapper.
        let mut preflight = self.tank.clone();
        preflight.advance(
            Time::new(0.0),
            flow,
            feed,
            numerics,
            || cancelled(),
        )?;

        let start = self.tank.elapsed().value();
        let end = add(start, duration.value())?;
        let scale = sum(
            self.tank
                .raw
                .iter()
                .copied()
                .chain(self.capacities.iter().copied()),
        )?;
        let tolerance = CyanidationNumerics {
            maximum_step: numerics.maximum_step,
            absolute_gold_error: numerics.absolute_amount_error,
            relative_gold_error: numerics.relative_amount_error,
            maximum_attempts: numerics.maximum_attempts,
            maximum_work: numerics.maximum_work,
        }
        .tolerance(scale)?;
        let before = self.tank.clone();
        let before_adsorbed = self.adsorbed.clone();
        let mut current = CoupledTrial {
            tank: self.tank.clone(),
            adsorbed: self.adsorbed.clone(),
            vacant: self.vacant.clone(),
        };
        let (mut clock, mut step) = (start, numerics.maximum_step.value().min(duration.value()));
        let (mut accepted, mut rejected) = (0usize, 0usize);
        let (mut largest, mut call_error) = (0.0_f64, 0.0_f64);

        while clock < end {
            if cancelled() {
                return Err(FlowCilError::Cancelled);
            }
            if accepted + rejected >= numerics.maximum_attempts {
                return Err(FlowCilError::Budget);
            }
            let finish = if clock + step >= end
                || end - (clock + step) <= 8.0 * f64::EPSILON * end.abs()
            {
                end
            } else {
                clock + step
            };
            let middle = clock.midpoint(finish);
            let h = finish - clock;
            if h <= 0.0 || middle <= clock || middle >= finish {
                return Err(FlowCilError::Arithmetic(
                    "continuous CIL coupling clock loses resolution",
                ));
            }
            let trial = (|| -> CilFlowResult<_> {
                let full = split(
                    current.clone(),
                    &self.pair,
                    &self.populations,
                    &self.capacities,
                    self.counterion_index,
                    flow,
                    feed,
                    numerics,
                    clock,
                    finish,
                    &mut cancelled,
                )?;
                let first = split(
                    current.clone(),
                    &self.pair,
                    &self.populations,
                    &self.capacities,
                    self.counterion_index,
                    flow,
                    feed,
                    numerics,
                    clock,
                    middle,
                    &mut cancelled,
                )?;
                let refined = split(
                    first,
                    &self.pair,
                    &self.populations,
                    &self.capacities,
                    self.counterion_index,
                    flow,
                    feed,
                    numerics,
                    middle,
                    finish,
                    &mut cancelled,
                )?;
                Ok((full, refined))
            })();
            let (full, refined) = match trial {
                Err(FlowCilError::Cyanidation(CyanidationError::TrialDepletion)) => {
                    rejected += 1;
                    step = h * 0.5;
                    continue;
                }
                other => other?,
            };
            let error = full.error(&refined)?;
            let allowance = product(tolerance, quotient(h, duration.value())?)?;
            let next = sum([call_error, error].into_iter())?;
            if error > allowance || next > tolerance {
                rejected += 1;
                step = h * 0.5;
                continue;
            }
            current = refined;
            call_error = next;
            largest = largest.max(error);
            accepted += 1;
            clock = finish;
            step = if error <= allowance / 8.0 {
                (h * 2.0).min(numerics.maximum_step.value())
            } else {
                h
            };
        }

        let state = admit_state(
            &current.tank,
            &self.pair,
            &self.populations,
            &self.capacities,
            self.counterion_index,
            &current.adsorbed,
            &current.vacant,
            numerics.maximum_work,
            &mut cancelled,
        )?;
        let interval_in = differences(&current.tank.inlet, &before.inlet)?;
        let interval_out = differences(&current.tank.outlet, &before.outlet)?;
        let interval_gold = difference(
            current.tank.gold_dissolved,
            before.gold_dissolved,
            "cumulative gold progress regressed",
        )?;
        let interval_conservation = conservation_cil_flow(
            &current.tank.layout,
            &self.pair,
            self.counterion_index,
            &self.capacities,
            &state,
            &before.raw,
            &before_adsorbed,
            &interval_in,
            &interval_out,
            interval_gold,
        )?;
        let cumulative_conservation = conservation_cil_flow(
            &current.tank.layout,
            &self.pair,
            self.counterion_index,
            &self.capacities,
            &state,
            &current.tank.layout.initial,
            &self.initial_adsorbed,
            &current.tank.inlet,
            &current.tank.outlet,
            current.tank.gold_dissolved,
        )?;
        let capture = |role| -> CilFlowResult<Amount> {
            Ok(Amount::new(sum(
                state
                    .populations
                    .iter()
                    .zip(&self.state.populations)
                    .filter(|(after, _)| after.role == role)
                    .map(|(after, prior)| after.adsorbed.value() - prior.adsorbed.value()),
            )?))
        };
        if cancelled() {
            return Err(FlowCilError::Cancelled);
        }
        let answer = FlowCilAdvance {
            pair: self.pair.clone(),
            populations: self.populations.clone(),
            previous_elapsed: Time::new(start),
            elapsed: Time::new(end),
            flow,
            feed: feed.clone(),
            state: state.clone(),
            inlet: inventory(&current.tank.layout, &interval_in),
            outlet: inventory(&current.tank.layout, &interval_out),
            cumulative_inlet: inventory(&current.tank.layout, &current.tank.inlet),
            cumulative_outlet: inventory(&current.tank.layout, &current.tank.outlet),
            interval_gold_dissolved: Amount::new(interval_gold),
            cumulative_gold_dissolved: Amount::new(current.tank.gold_dissolved),
            interval_product_capture: capture(CarbonRole::ProductCarbon)?,
            interval_preg_robbing_capture: capture(CarbonRole::PregRobbing)?,
            interval_conservation,
            cumulative_conservation,
            maximum_coupling_error: Amount::new(largest),
            accumulated_coupling_error: Amount::new(accumulated),
            liquid_error_estimate: Amount::new(current.tank.accumulated_error),
            accepted_steps: accepted,
            rejected_steps: rejected,
        };
        self.tank = current.tank;
        self.adsorbed = current.adsorbed;
        self.vacant = current.vacant;
        self.state = state;
        self.coupling_error = accumulated;
        Ok(answer)
    }
}

#[allow(clippy::too_many_arguments)]
fn split(
    mut state: CoupledTrial,
    pair: &IonPairTransfer,
    populations: &[CarbonPopulation],
    capacities: &[f64],
    counterion_index: usize,
    flow: VolumetricFlowRate,
    feed: &CyanidationFeed,
    numerics: FlowNumerics,
    start: f64,
    end: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> CilFlowResult<CoupledTrial> {
    let half = product(end - start, 0.5)?;
    adsorption_step(
        &mut state,
        pair,
        populations,
        capacities,
        counterion_index,
        half,
        start,
        numerics.maximum_work,
        cancelled,
    )?;
    state.tank.advance(
        Time::new(end - start),
        flow,
        feed,
        FlowNumerics {
            maximum_step: Time::new(numerics.maximum_step.value().min(end - start)),
            ..numerics
        },
        || cancelled(),
    )?;
    adsorption_step(
        &mut state,
        pair,
        populations,
        capacities,
        counterion_index,
        half,
        end,
        numerics.maximum_work,
        cancelled,
    )?;
    Ok(state)
}

#[allow(clippy::too_many_arguments)]
fn adsorption_step(
    state: &mut CoupledTrial,
    pair: &IonPairTransfer,
    populations: &[CarbonPopulation],
    capacities: &[f64],
    counterion_index: usize,
    duration: f64,
    elapsed: f64,
    maximum_work: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> CilFlowResult<()> {
    if cancelled() {
        return Err(FlowCilError::Cancelled);
    }
    if duration == 0.0 {
        return Ok(());
    }
    let initial = state.clone();
    let first = adsorption_rates(
        &initial,
        pair,
        populations,
        capacities,
        counterion_index,
        elapsed,
        maximum_work,
        cancelled,
    )?;
    let second_state = shifted_adsorption(
        &initial,
        &first,
        pair,
        populations,
        capacities,
        counterion_index,
        0.5 * duration,
        elapsed,
        maximum_work,
        cancelled,
    )?;
    let second = adsorption_rates(
        &second_state,
        pair,
        populations,
        capacities,
        counterion_index,
        elapsed,
        maximum_work,
        cancelled,
    )?;
    let third_state = shifted_adsorption(
        &initial,
        &second,
        pair,
        populations,
        capacities,
        counterion_index,
        0.5 * duration,
        elapsed,
        maximum_work,
        cancelled,
    )?;
    let third = adsorption_rates(
        &third_state,
        pair,
        populations,
        capacities,
        counterion_index,
        elapsed,
        maximum_work,
        cancelled,
    )?;
    let fourth_state = shifted_adsorption(
        &initial,
        &third,
        pair,
        populations,
        capacities,
        counterion_index,
        duration,
        elapsed,
        maximum_work,
        cancelled,
    )?;
    let fourth = adsorption_rates(
        &fourth_state,
        pair,
        populations,
        capacities,
        counterion_index,
        elapsed,
        maximum_work,
        cancelled,
    )?;
    let average = average_rates([&first, &second, &third, &fourth])?;
    *state = shifted_adsorption(
        &initial,
        &average,
        pair,
        populations,
        capacities,
        counterion_index,
        duration,
        elapsed,
        maximum_work,
        cancelled,
    )?;
    Ok(())
}

fn adsorption_rates(
    state: &CoupledTrial,
    pair: &IonPairTransfer,
    populations: &[CarbonPopulation],
    capacities: &[f64],
    counterion_index: usize,
    elapsed: f64,
    maximum_work: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> CilFlowResult<Vec<f64>> {
    let admitted = admit_state(
        &state.tank,
        pair,
        populations,
        capacities,
        counterion_index,
        &state.adsorbed,
        &state.vacant,
        maximum_work,
        cancelled,
    )?;
    admitted
        .populations
        .iter()
        .zip(populations)
        .map(|(actual, population)| -> CilFlowResult<f64> {
            let magnitude = product(
                actual.transfer.loading_rate.value().abs(),
                population.sorbent().mass().value(),
            )?;
            Ok(if actual.transfer.loading_rate.value() < 0.0 {
                -magnitude
            } else {
                magnitude
            })
        })
        .collect()
}

#[allow(clippy::too_many_arguments)]
fn shifted_adsorption(
    initial: &CoupledTrial,
    rates: &[f64],
    pair: &IonPairTransfer,
    populations: &[CarbonPopulation],
    capacities: &[f64],
    counterion_index: usize,
    duration: f64,
    elapsed: f64,
    maximum_work: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> CilFlowResult<CoupledTrial> {
    let total = sum(rates.iter().copied())?;
    let mut next = initial.clone();
    next.tank.raw[3] = signed_shift(
        initial.tank.raw[3],
        -total,
        duration,
        f64::MAX,
    )?;
    next.tank.raw[counterion_index] = signed_shift(
        initial.tank.raw[counterion_index],
        -total,
        duration,
        f64::MAX,
    )?;
    for i in 0..rates.len() {
        next.adsorbed[i] = signed_shift(
            initial.adsorbed[i],
            rates[i],
            duration,
            capacities[i],
        )?;
        next.vacant[i] = signed_shift(
            initial.vacant[i],
            -rates[i],
            duration,
            capacities[i],
        )?;
    }
    let admitted = admit_state(
        &next.tank,
        pair,
        populations,
        capacities,
        counterion_index,
        &next.adsorbed,
        &next.vacant,
        maximum_work,
        cancelled,
    )?;
    next.tank.state = admitted.solution;
    Ok(next)
}

fn signed_shift(initial: f64, rate: f64, duration: f64, maximum: f64) -> CilFlowResult<f64> {
    if rate == 0.0 || duration == 0.0 {
        return Ok(initial);
    }
    let change = product(rate.abs(), duration)?;
    let next = sum(
        [initial, if rate < 0.0 { -change } else { change }].into_iter(),
    )?;
    if next < 0.0 || next > maximum || (initial > 0.0 && next == 0.0) {
        return Err(CyanidationError::TrialDepletion.into());
    }
    if next == initial {
        return Err(FlowCilError::Arithmetic(
            "nonzero carbon transfer is below inventory resolution",
        ));
    }
    Ok(next)
}

fn signed_fraction(value: f64, divisor: f64) -> CilFlowResult<f64> {
    let magnitude = quotient(value.abs(), divisor)?;
    Ok(if value < 0.0 { -magnitude } else { magnitude })
}
fn average_rates(rates: [&Vec<f64>; 4]) -> CilFlowResult<Vec<f64>> {
    (0..rates[0].len())
        .map(|i| {
            Ok(sum(
                [rates[0][i], rates[1][i], rates[2][i], rates[3][i]]
                    .into_iter()
                    .zip([6.0, 3.0, 3.0, 6.0])
                    .map(|(value, divisor)| signed_fraction(value, divisor))
                    .collect::<CilFlowResult<Vec<_>>>()?
                    .into_iter(),
            )?)
        })
        .collect()
}

#[allow(clippy::too_many_arguments)]
fn admit_state(
    tank: &CyanidationCstr,
    pair: &IonPairTransfer,
    populations: &[CarbonPopulation],
    capacities: &[f64],
    counterion_index: usize,
    adsorbed: &[f64],
    vacant: &[f64],
    maximum_work: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> CilFlowResult<FlowCilState> {
    if cancelled() {
        return Err(FlowCilError::Cancelled);
    }
    let mut proxy = || cancelled();
    let mut work = dosing::work(maximum_work, &mut proxy)?;
    let solution = evaluate(
        &tank.layout,
        &tank.raw,
        tank.elapsed().value(),
        &mut work,
    )?;
    let volume = tank.layout.charge.solution_volume;
    let concentration = AmountConcentration::new(quotient(tank.raw[3], volume.value())?);
    let counterion = tank.raw[counterion_index];
    let mut states = Vec::with_capacity(populations.len());
    for (index, population) in populations.iter().enumerate() {
        work.tick()?;
        let sorbent = population.sorbent();
        let loading = AdsorbentLoading::new(quotient(
            adsorbed[index],
            sorbent.mass().value(),
        )?);
        let mut conditions = sorbent.conditions().clone();
        for (axis, dims) in owned_axes() {
            if matches!(axis, TEMPERATURE_AXIS | VOLUME_AXIS) {
                continue;
            }
            let value = if axis == ELAPSED_TIME_AXIS {
                tank.elapsed().value()
            } else {
                *solution
                    .conditions
                    .axes()
                    .get(axis)
                    .ok_or(FlowCilError::Invalid(
                        "missing actual solution coordinate for carbon source",
                    ))?
            };
            conditions = conditions
                .with_quantity(axis, QuantitySpec::dimensional(dims), value)
                .map_err(CyanidationError::from)?;
        }
        conditions = conditions
            .with_quantity(
                COUNTERION_AXIS,
                QuantitySpec::dimensional(AmountConcentration::DIMS),
                quotient(counterion, volume.value())?,
            )
            .map_err(CyanidationError::from)?;
        let transfer = sorbent.model().rate_at(
            concentration,
            loading,
            volume,
            tank.layout.charge.temperature,
            &conditions,
        )?;
        states.push(CilSorbentState {
            id: sorbent.id().into(),
            role: population.role(),
            adsorbed: Amount::new(adsorbed[index]),
            vacant_capacity: Amount::new(vacant[index]),
            loading,
            transfer,
        });
    }
    let state = FlowCilState {
        elapsed: tank.elapsed(),
        solution,
        populations: states,
        dissolved_counterion: Amount::new(counterion),
    };
    site_balance(&state.populations, capacities)?;
    Ok(state)
}

fn site_balance(states: &[CilSorbentState], capacities: &[f64]) -> CilFlowResult<f64> {
    let mut largest = 0.0_f64;
    for (state, &capacity) in states.iter().zip(capacities) {
        let residual = sum(
            [
                state.adsorbed.value(),
                state.vacant_capacity.value(),
                -capacity,
            ]
            .into_iter(),
        )?;
        if residual.abs() > 512.0 * f64::EPSILON * capacity.max(1.0) {
            return Err(FlowCilError::Arithmetic(
                "retained carbon site inventory does not close",
            ));
        }
        largest = largest.max(residual.abs());
    }
    Ok(largest)
}

#[allow(clippy::too_many_arguments)]
fn conservation_cil_flow(
    layout: &FlowLayout,
    pair: &IonPairTransfer,
    counterion_index: usize,
    capacities: &[f64],
    state: &FlowCilState,
    initial_raw: &[f64],
    initial_adsorbed: &[f64],
    incoming: &[f64],
    outgoing: &[f64],
    gold: f64,
) -> CilFlowResult<FlowCilConservation> {
    let final_adsorbed = sum(state.populations.iter().map(|p| p.adsorbed.value()))?;
    let initial_adsorbed = sum(initial_adsorbed.iter().copied())?;
    let adsorption = final_adsorbed - initial_adsorbed;
    let mut adjusted_in = incoming.to_vec();
    let mut adjusted_out = outgoing.to_vec();
    if adsorption >= 0.0 {
        adjusted_out[3] = add(adjusted_out[3], adsorption)?;
        adjusted_out[counterion_index] = add(adjusted_out[counterion_index], adsorption)?;
    } else {
        let released = -adsorption;
        adjusted_in[3] = add(adjusted_in[3], released)?;
        adjusted_in[counterion_index] = add(adjusted_in[counterion_index], released)?;
    }
    let adjusted = conservation(
        layout,
        &state.solution,
        initial_raw,
        &adjusted_in,
        &adjusted_out,
        gold,
    )?;
    let counterion = adjusted
        .components
        .spectators
        .iter()
        .find(|ion| &ion.species == pair.counterion())
        .ok_or(FlowCilError::Invalid(
            "counterion residual missing from component account",
        ))?
        .amount;
    Ok(FlowCilConservation {
        gold: adjusted.gold,
        cyanide: adjusted.cyanide,
        counterion,
        dissolved_oxygen: adjusted.dissolved_oxygen,
        maximum_site_residual: Amount::new(site_balance(
            &state.populations,
            capacities,
        )?),
    })
}

fn difference(after: f64, before: f64, reason: &'static str) -> CilFlowResult<f64> {
    let value = after - before;
    if !value.is_finite() || value < 0.0 {
        Err(FlowCilError::Arithmetic(reason))
    } else {
        Ok(value)
    }
}
fn differences(after: &[f64], before: &[f64]) -> CilFlowResult<Vec<f64>> {
    after
        .iter()
        .zip(before)
        .map(|(&a, &b)| difference(a, b, "cumulative liquid stream regressed"))
        .collect()
}
fn require_local(condition: bool, reason: &'static str) -> CilFlowResult<()> {
    if condition {
        Ok(())
    } else {
        Err(FlowCilError::Invalid(reason))
    }
}

#[cfg(test)]
#[path = "cil/tests.rs"]
mod tests;
