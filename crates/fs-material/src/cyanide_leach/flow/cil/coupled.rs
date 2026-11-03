//! Continuous, aerated gold leaching with explicitly retained carbon populations.
//!
//! The original mixed-tank owner transports every liquid/slurry analytical
//! component and solves local cyanidation with live acid/base and metal-ligand
//! equilibrium. The original gas owner supplies sourced solubility, kLa and
//! independent uptake/stripping. Original Langmuir owners supply reversible
//! loading of one declared neutral gold/counterion pair.
//!
//! All supplied sorbents are retained at fixed dry mass, including a population
//! labeled PregRobbing. That role identifies a destination, not entrainment;
//! this owner does not predict washout of naturally carbonaceous ore solids.
//! The mobile feed and initial tank must supply at least one declared
//! countercation per accessible or aqueous gold atom. The invariant
//! C-D-S >= 0 is preserved by neutral capture, dissolution and equal-residence
//! transport, without manufacturing counterion-limited kinetic authority.
//!
//! Symmetric transport/gas/adsorption/reaction/adsorption/gas/transport splitting
//! evolves mobile stocks, loaded/vacant sites, liquid exports, dissolution and
//! gross gas amounts independently. Full/two-half raw differences control
//! their joint L1 local error; no smooth-order divisor is assumed at gas-flux
//! reversals. Stage source admission is not an enclosure of the unsplit path.
//! Fixed volume, temperature and activity basis remain explicit. No carbon
//! harvesting, regeneration, entrained-carbon transport or plant calibration
//! is inferred.

use super::super::aeration::{AeratedCyanidationCstr, AeratedProblem, Boundary, GasTrial};
use super::super::*;
use crate::aqueous_adsorption::{
    AdsorbentLoading, AqueousAdsorptionError,
    flow::{AdsorptionFlowError, MAX_POPULATIONS},
};
use crate::cyanide_leach::aeration::{
    OXYGEN_PARTIAL_PRESSURE_AXIS, OxygenTransfer, SATURATION_CONCENTRATION_AXIS,
};
use crate::cyanide_leach::cil::{
    COUNTERION_AXIS, CarbonPopulation, CarbonRole, CilConservation, CilError, CilSorbentState,
    IonPairTransfer,
};
use fs_matdb::MaterialAnswer;
use fs_qty::Pressure;

type Result<T> = core::result::Result<T, CilError>;

#[derive(Clone, Debug)]
struct CaptureLayout {
    pair: IonPairTransfer,
    populations: Vec<CarbonPopulation>,
    capacities: Vec<f64>,
    initial_adsorbed: Vec<f64>,
    counterion_component: usize,
    counterion_spectator: usize,
}

/// Actual shared flowing solution, retained carbon and prescribed gas boundary.
#[derive(Clone, Debug)]
pub struct CilFlowState {
    pub elapsed: Time,
    pub solution: CyanidationFlowState,
    pub populations: Vec<CilSorbentState>,
    pub dissolved_counterion: Amount,
    pub oxygen_partial_pressure: Pressure,
    pub oxygen_saturation: Amount,
    /// Original gas-transfer source point at this actual accepted state.
    pub transfer_conditions: QueryPoint,
}

/// Continuous CIL retaining original liquid, gas, chemical and sorbent owners.
#[derive(Clone, Debug)]
pub struct AeratedCilCstr {
    aerated: AeratedCyanidationCstr,
    capture: Arc<CaptureLayout>,
    state: CilFlowState,
}

/// One atomic coupled interval. Liquid streams remain separate from gas and carbon.
#[derive(Clone, Debug)]
pub struct CilFlowAdvance {
    /// Original liquid/source receipt, with balances including explicit carbon/gas exchange.
    pub flow: CyanidationFlowAdvance,
    pub state: CilFlowState,
    /// Signed net uptake by product carbon during this interval.
    pub interval_product_capture: Amount,
    /// Signed net uptake by explicitly retained preg-robbing material.
    pub interval_preg_robbing_capture: Amount,
    pub interval_conservation: CilConservation,
    pub cumulative_conservation: CilConservation,
    pub oxygen_partial_pressure: Pressure,
    pub oxygen_saturation: Amount,
    /// Pinned measured-Hcp receipt at the held temperature, when selected.
    /// The Gibbs-equilibrium route retains `None` without invented measurements.
    pub henry_answer: Option<MaterialAnswer>,
    pub oxygen_in: Amount,
    pub oxygen_out: Amount,
    pub cumulative_oxygen_in: Amount,
    pub cumulative_oxygen_out: Amount,
    pub transfer_conditions: QueryPoint,
    capture: Arc<CaptureLayout>,
    transfer: Arc<OxygenTransfer>,
}
impl CilFlowAdvance {
    pub fn model(&self) -> &LocalGoldCyanidation {
        self.flow.model()
    }
    pub fn transfer(&self) -> &OxygenTransfer {
        &self.transfer
    }
    pub fn ion_pair(&self) -> &IonPairTransfer {
        &self.capture.pair
    }
    pub fn models(&self) -> &[CarbonPopulation] {
        &self.capture.populations
    }
}

impl AeratedCilCstr {
    /// Attach retained carbon and the original gas owner to a fresh or advanced tank.
    /// Pressure/extras explicitly admit every initial gas and adsorption source.
    #[allow(clippy::too_many_arguments)]
    pub fn from_tank(
        tank: CyanidationCstr,
        transfer: OxygenTransfer,
        pair: IonPairTransfer,
        populations: Vec<CarbonPopulation>,
        oxygen_partial_pressure: Pressure,
        transfer_conditions: &QueryPoint,
    ) -> Result<Self> {
        Self::from_aerated(
            AeratedCyanidationCstr::from_tank(tank, transfer)?,
            pair,
            populations,
            oxygen_partial_pressure,
            transfer_conditions,
        )
    }

    /// Retain earlier gas/liquid exports, reaction progress, clock and error history.
    /// Preloaded carbon brings its own complete bound-pair inventory at attachment.
    pub fn from_aerated(
        aerated: AeratedCyanidationCstr,
        pair: IonPairTransfer,
        mut populations: Vec<CarbonPopulation>,
        oxygen_partial_pressure: Pressure,
        transfer_conditions: &QueryPoint,
    ) -> Result<Self> {
        ensure(
            populations.len() <= MAX_POPULATIONS,
            "too many retained CIL populations",
        )?;
        populations.sort_by(|a, b| a.sorbent().id().cmp(b.sorbent().id()));
        ensure(
            populations
                .windows(2)
                .all(|p| p[0].sorbent().id() != p[1].sorbent().id()),
            "retained CIL population IDs must be unique",
        )?;
        let tank = &aerated.tank;
        let layout = &tank.layout;
        let model = &layout.model;
        for species in [pair.counterion(), pair.adsorbed_pair()] {
            ensure(
                model
                    .chemistry
                    .families()
                    .iter()
                    .all(|f| !f.species().contains(species))
                    && layout.competition.as_ref().is_none_or(|m| {
                        m.model
                            .families()
                            .iter()
                            .all(|f| !f.species().contains(species))
                    }),
                "pair species cannot also belong to a dissolved chemical family",
            )?;
        }
        ensure(
            !layout
                .charge
                .spectators
                .iter()
                .any(|s| &s.species == pair.adsorbed_pair()),
            "the bound pair cannot also be a liquid spectator",
        )?;
        let counterion_spectator = layout
            .charge
            .spectators
            .iter()
            .position(|s| &s.species == pair.counterion())
            .ok_or(CilError::Invalid(
                "declare the mobile pair counterion explicitly",
            ))?;
        ensure(
            layout.charge.spectators[counterion_spectator].charge == 1,
            "the pair counterion must retain charge +1",
        )?;
        let aqueous = layout
            .competition
            .as_ref()
            .map_or(model.chemistry.reaction_network(), |m| {
                m.model.reaction_network()
            });
        ReactionNetwork::new(
            aqueous
                .reactions()
                .iter()
                .cloned()
                .chain([
                    model.reaction.clone(),
                    aerated.transfer.reaction()?,
                    pair.reaction().clone(),
                ])
                .collect(),
        )
        .map_err(CyanidationError::from)?;
        let counterion_component = 4 + layout.charge.other_families.len() + counterion_spectator;
        reserve(&tank.raw, counterion_component)?;
        for population in &populations {
            ensure(
                [OXYGEN_PARTIAL_PRESSURE_AXIS, SATURATION_CONCENTRATION_AXIS]
                    .iter()
                    .all(|axis| !population.sorbent().conditions().axes().contains_key(*axis)),
                "CIL owns actual gas pressure and oxygen saturation",
            )?;
        }
        let capacities = populations
            .iter()
            .map(|p| {
                product(
                    p.sorbent().mass().value(),
                    p.sorbent().model().equilibrium().capacity().value(),
                )
            })
            .collect::<super::super::Result<Vec<_>>>()?;
        let initial_adsorbed: Vec<_> = populations
            .iter()
            .map(|p| p.sorbent().initial_adsorbed().value())
            .collect();
        let vacant: Vec<_> = capacities
            .iter()
            .zip(&initial_adsorbed)
            .map(|(c, a)| c - a)
            .collect();
        let capture = Arc::new(CaptureLayout {
            pair,
            populations,
            capacities,
            initial_adsorbed,
            counterion_component,
            counterion_spectator,
        });
        let (boundary, _) = boundary(
            &aerated,
            &capture,
            oxygen_partial_pressure,
            transfer_conditions,
        )?;
        let state = carbon_state(
            &capture,
            &tank.state,
            &capture.initial_adsorbed,
            &vacant,
            &boundary,
        )?;
        Ok(Self {
            aerated,
            capture,
            state,
        })
    }

    pub fn state(&self) -> &CilFlowState {
        &self.state
    }
    pub fn model(&self) -> &LocalGoldCyanidation {
        self.aerated.tank.model()
    }
    pub fn transfer(&self) -> &OxygenTransfer {
        &self.aerated.transfer
    }
    pub fn ion_pair(&self) -> &IonPairTransfer {
        &self.capture.pair
    }
    pub fn models(&self) -> &[CarbonPopulation] {
        &self.capture.populations
    }
    pub fn initial_charge(&self) -> &CyanidationCharge {
        self.aerated.tank.initial_charge()
    }
    pub fn cumulative_inlet(&self) -> CyanidationInventory {
        self.aerated.tank.cumulative_inlet()
    }
    pub fn cumulative_outlet(&self) -> CyanidationInventory {
        self.aerated.tank.cumulative_outlet()
    }
    pub fn cumulative_gold_dissolved(&self) -> Amount {
        self.aerated.tank.cumulative_gold_dissolved()
    }
    pub fn oxygen_in(&self) -> Amount {
        Amount::new(self.aerated.oxygen_in)
    }
    pub fn oxygen_out(&self) -> Amount {
        Amount::new(self.aerated.oxygen_out)
    }

    /// Advance constant feed/flow/pressure, allowing all three to change between calls.
    /// Zero flow preserves retained carbon while dissolution and gas exchange continue.
    #[allow(clippy::too_many_arguments)]
    pub fn advance(
        &mut self,
        duration: Time,
        flow: VolumetricFlowRate,
        feed: &CyanidationFeed,
        oxygen_partial_pressure: Pressure,
        transfer_conditions: &QueryPoint,
        numerics: FlowNumerics,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CilFlowAdvance> {
        self.advance_inner(
            duration,
            flow,
            feed,
            oxygen_partial_pressure,
            transfer_conditions,
            numerics,
            &mut cancelled,
        )
    }

    // One driver instantiation retains the same callback/work contract for all
    // callers, instead of duplicating the whole nested chemical solver per closure.
    #[allow(clippy::too_many_arguments)]
    fn advance_inner(
        &mut self,
        duration: Time,
        flow: VolumetricFlowRate,
        feed: &CyanidationFeed,
        oxygen_partial_pressure: Pressure,
        transfer_conditions: &QueryPoint,
        numerics: FlowNumerics,
        mut cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<CilFlowAdvance> {
        nonnegative(
            duration.value(),
            "continuous CIL duration must be nonnegative",
        )?;
        nonnegative(flow.value(), "continuous CIL flow must be nonnegative")?;
        let tank = &self.aerated.tank;
        let layout = &tank.layout;
        let start = tank.elapsed().value();
        let end = add(start, duration.value())?;
        let (boundary, henry_answer) = boundary(
            &self.aerated,
            &self.capture,
            oxygen_partial_pressure,
            transfer_conditions,
        )?;
        let (minimum, maximum) = self
            .aerated
            .transfer
            .source()
            .validity
            .bound(ELAPSED_TIME_AXIS)
            .expect("original gas clock domain");
        if !(minimum..=maximum).contains(&end) {
            return Err(CyanidationError::OutsideValidity {
                axis: ELAPSED_TIME_AXIS.into(),
                value: end,
                minimum,
                maximum,
            }
            .into());
        }
        let mut work = dosing::work(numerics.maximum_work, &mut cancelled)?;
        work.tick()?;
        let initial_adsorbed: Vec<_> = self
            .state
            .populations
            .iter()
            .map(|p| p.adsorbed.value())
            .collect();
        let initial_vacant: Vec<_> = self
            .state
            .populations
            .iter()
            .map(|p| p.vacant_capacity.value())
            .collect();
        carbon_state(
            &self.capture,
            &tank.state,
            &initial_adsorbed,
            &initial_vacant,
            &boundary,
        )?;
        let concentration = feed_values(layout, feed, &mut work)?;
        reserve(&concentration, self.capture.counterion_component)?;
        let volume = layout.charge.solution_volume.value();
        let exposure = quotient(product(flow.value(), duration.value())?, volume)?;
        let target = concentration
            .iter()
            .map(|c| product(*c, volume))
            .collect::<super::super::Result<Vec<_>>>()?;
        let incoming = target
            .iter()
            .map(|n| product(*n, exposure))
            .collect::<super::super::Result<Vec<_>>>()?;
        let saturation = boundary.saturation;
        let gas_scale = product(
            product(
                self.aerated.transfer.coefficient().value(),
                duration.value(),
            )?,
            tank.raw[2].max(saturation.value()).max(target[2]),
        )?;
        let scale = sum(tank
            .raw
            .iter()
            .chain(&incoming)
            .chain(&initial_adsorbed)
            .copied()
            .chain(std::iter::once(gas_scale)))?;
        let tolerance = CyanidationNumerics {
            maximum_step: numerics.maximum_step,
            absolute_gold_error: numerics.absolute_amount_error,
            relative_gold_error: numerics.relative_amount_error,
            maximum_attempts: numerics.maximum_attempts,
            maximum_work: numerics.maximum_work,
        }
        .tolerance(scale)?;
        let problem = CoupledProblem {
            gas: AeratedProblem {
                liquid: Problem {
                    layout,
                    target: &target,
                    flow,
                },
                boundary,
            },
            capture: &self.capture,
        };
        let mut current = CaptureTrial {
            gas: GasTrial::new(tank.raw.clone()),
            adsorbed: initial_adsorbed.clone(),
            vacant: initial_vacant,
        };
        let (mut clock, mut step) = (start, numerics.maximum_step.value().min(duration.value()));
        let (mut accepted, mut rejected, mut largest, mut error_sum) = (0, 0, 0.0_f64, 0.0_f64);
        while clock < end {
            work.tick()?;
            if accepted + rejected >= numerics.maximum_attempts {
                return Err(CilError::Budget);
            }
            let remaining = end - clock;
            // Balance the final two proposed steps instead of leaving a tiny
            // floating-point remainder after a long sequence of accepted steps.
            // Both halves remain within the configured maximum step and both
            // are integrated under the ordinary coupled error controller.
            let final_pair = remaining > step && remaining - step < step;
            let finish = if remaining <= step {
                end
            } else if final_pair {
                clock.midpoint(end)
            } else {
                clock + step
            };
            let middle = clock.midpoint(finish);
            let h = finish - clock;
            if h <= 0.0 || middle <= clock || middle >= finish {
                return Err(CilError::Arithmetic(
                    "continuous CIL clock loses resolution",
                ));
            }
            let initial = CaptureTrial {
                gas: GasTrial::new(current.gas.liquid.raw.clone()),
                adsorbed: current.adsorbed.clone(),
                vacant: current.vacant.clone(),
            };
            let trial = (|| -> Result<_> {
                let coarse = problem.split(initial.clone(), clock, finish, &mut work)?;
                let first = problem.split(initial, clock, middle, &mut work)?;
                let fine = problem.split(first, middle, finish, &mut work)?;
                Ok((coarse, fine))
            })();
            let (coarse, fine) = match trial {
                Ok(pair) => pair,
                Err(error) if retryable(&error) => {
                    rejected += 1;
                    step = h * 0.5;
                    continue;
                }
                Err(error) => return Err(error),
            };
            let error = sum(coarse
                .values()
                .zip(fine.values())
                .map(|(a, b)| (a - b).abs()))?;
            let allowance = product(tolerance, quotient(h, duration.value())?)?;
            let next_error = sum([error_sum, error].into_iter())?;
            if error > allowance || next_error > tolerance {
                rejected += 1;
                step = h * 0.5;
                continue;
            }
            current.gas.liquid.raw = fine.gas.liquid.raw;
            for (total, increment) in current
                .gas
                .liquid
                .outlet
                .iter_mut()
                .zip(fine.gas.liquid.outlet)
            {
                *total = add(*total, increment)?;
            }
            current.gas.liquid.gold = add(current.gas.liquid.gold, fine.gas.liquid.gold)?;
            current.gas.inward = add(current.gas.inward, fine.gas.inward)?;
            current.gas.outward = add(current.gas.outward, fine.gas.outward)?;
            current.adsorbed = fine.adsorbed;
            current.vacant = fine.vacant;
            error_sum = next_error;
            largest = largest.max(error);
            accepted += 1;
            clock = finish;
            step = if final_pair {
                (end - finish).min(numerics.maximum_step.value())
            } else if error <= allowance / 8.0 {
                (h * 2.0).min(numerics.maximum_step.value())
            } else {
                h
            };
        }
        let solution = evaluate(layout, &current.gas.liquid.raw, end, &mut work)?;
        let state = carbon_state(
            &self.capture,
            &solution,
            &current.adsorbed,
            &current.vacant,
            &problem.gas.boundary,
        )?;
        let cumulative_inlet = tank
            .inlet
            .iter()
            .zip(&incoming)
            .map(|(a, b)| add(*a, *b))
            .collect::<super::super::Result<Vec<_>>>()?;
        let cumulative_outlet = tank
            .outlet
            .iter()
            .zip(&current.gas.liquid.outlet)
            .map(|(a, b)| add(*a, *b))
            .collect::<super::super::Result<Vec<_>>>()?;
        let gold = add(tank.gold_dissolved, current.gas.liquid.gold)?;
        let oxygen_in = add(self.aerated.oxygen_in, current.gas.inward)?;
        let oxygen_out = add(self.aerated.oxygen_out, current.gas.outward)?;
        let accumulated_error = sum([tank.accumulated_error, error_sum].into_iter())?;
        let (interval_liquid, interval_conservation) = accounts(
            layout,
            &self.capture,
            &state,
            &tank.raw,
            &incoming,
            &current.gas.liquid.outlet,
            current.gas.liquid.gold,
            current.gas.inward,
            current.gas.outward,
            &initial_adsorbed,
        )?;
        let (cumulative_liquid, cumulative_conservation) = accounts(
            layout,
            &self.capture,
            &state,
            &layout.initial,
            &cumulative_inlet,
            &cumulative_outlet,
            gold,
            oxygen_in,
            oxygen_out,
            &self.capture.initial_adsorbed,
        )?;
        let capture = |role| -> Result<Amount> {
            Ok(Amount::new(sum(self
                .capture
                .populations
                .iter()
                .zip(&current.adsorbed)
                .zip(&initial_adsorbed)
                .filter(|((p, _), _)| p.role() == role)
                .map(|((_, after), before)| after - before))?))
        };
        let stoichiometric_amounts = self
            .model()
            .reaction()
            .amounts(Amount::new(quotient(gold, 4.0)?))
            .map_err(CyanidationError::from)?;
        let mut answer = CilFlowAdvance {
            flow: CyanidationFlowAdvance {
                layout: layout.clone(),
                previous_elapsed: Time::new(start),
                elapsed: Time::new(end),
                flow,
                feed: feed.clone(),
                state: solution.clone(),
                inlet: inventory(layout, &incoming),
                outlet: inventory(layout, &current.gas.liquid.outlet),
                cumulative_inlet: inventory(layout, &cumulative_inlet),
                cumulative_outlet: inventory(layout, &cumulative_outlet),
                interval_gold_dissolved: Amount::new(current.gas.liquid.gold),
                cumulative_gold_dissolved: Amount::new(gold),
                stoichiometric_amounts,
                interval_conservation: interval_liquid,
                cumulative_conservation: cumulative_liquid,
                maximum_local_error: Amount::new(largest),
                accumulated_error_estimate: Amount::new(accumulated_error),
                accepted_steps: accepted,
                rejected_steps: rejected,
                work: 0,
            },
            interval_product_capture: capture(CarbonRole::ProductCarbon)?,
            interval_preg_robbing_capture: capture(CarbonRole::PregRobbing)?,
            interval_conservation,
            cumulative_conservation,
            oxygen_partial_pressure,
            oxygen_saturation: saturation,
            henry_answer,
            oxygen_in: Amount::new(current.gas.inward),
            oxygen_out: Amount::new(current.gas.outward),
            cumulative_oxygen_in: Amount::new(oxygen_in),
            cumulative_oxygen_out: Amount::new(oxygen_out),
            transfer_conditions: state.transfer_conditions.clone(),
            state: state.clone(),
            capture: self.capture.clone(),
            transfer: self.aerated.transfer.clone(),
        };
        work.tick()?;
        answer.flow.work = work.used;
        let tank = &mut self.aerated.tank;
        tank.raw = current.gas.liquid.raw;
        tank.state = solution;
        tank.inlet = cumulative_inlet;
        tank.outlet = cumulative_outlet;
        tank.gold_dissolved = gold;
        tank.accumulated_error = accumulated_error;
        self.aerated.oxygen_in = oxygen_in;
        self.aerated.oxygen_out = oxygen_out;
        self.state = state;
        Ok(answer)
    }
}

fn boundary<'a>(
    aerated: &'a AeratedCyanidationCstr,
    capture: &CaptureLayout,
    pressure: Pressure,
    extra: &'a QueryPoint,
) -> Result<(Boundary<'a>, Option<MaterialAnswer>)> {
    ensure(
        !aerated
            .tank
            .layout
            .conditions
            .axes()
            .contains_key(COUNTERION_AXIS)
            && !extra.axes().contains_key(COUNTERION_AXIS),
        "CIL owns the actual mobile counterion source coordinate",
    )?;
    aerated
        .transfer
        .validate_extras(&aerated.tank.layout.conditions, extra)?;
    for population in &capture.populations {
        aerated
            .transfer
            .validate_extras(population.sorbent().conditions(), extra)?;
    }
    // CIL owns this coordinate before either original gas source is queried.
    // Henry claims retain only their admitted temperature/reference flags;
    // kLa continues to admit the complete changing liquid point at each stage.
    let chemical_conditions = aerated
        .tank
        .state
        .conditions
        .clone()
        .with_quantity(
            COUNTERION_AXIS,
            QuantitySpec::dimensional(AmountConcentration::DIMS),
            quotient(
                aerated.tank.raw[capture.counterion_component],
                aerated.tank.layout.charge.solution_volume.value(),
            )?,
        )
        .map_err(CyanidationError::from)?;
    let (saturation, henry_answer) = aerated.transfer.equilibrium_saturation(
        pressure,
        aerated.tank.layout.charge.amount_at_unit_activity,
        &chemical_conditions,
        extra,
    )?;
    Ok((
        Boundary {
            transfer: &aerated.transfer,
            pressure,
            saturation,
            extra,
        },
        henry_answer,
    ))
}

fn counterion_extra(boundary: &Boundary<'_>, amount: f64, volume: f64) -> Result<QueryPoint> {
    Ok(boundary
        .extra
        .clone()
        .with_quantity(
            COUNTERION_AXIS,
            QuantitySpec::dimensional(AmountConcentration::DIMS),
            quotient(amount, volume)?,
        )
        .map_err(CyanidationError::from)?)
}

fn carbon_state(
    capture: &CaptureLayout,
    solution: &CyanidationFlowState,
    adsorbed: &[f64],
    vacant: &[f64],
    boundary: &Boundary<'_>,
) -> Result<CilFlowState> {
    let volume = Volume::new(solution.conditions.axes()[VOLUME_AXIS]);
    let temperature = Temperature::new(solution.conditions.axes()[TEMPERATURE_AXIS]);
    let counterion = solution.inventory.spectators[capture.counterion_spectator].amount;
    let gas_extra = counterion_extra(boundary, counterion.value(), volume.value())?;
    let transfer_conditions = Boundary {
        transfer: boundary.transfer,
        pressure: boundary.pressure,
        saturation: boundary.saturation,
        extra: &gas_extra,
    }
    .admit(solution)?;
    let concentration = AmountConcentration::new(quotient(
        solution.inventory.gold_complex.value(),
        volume.value(),
    )?);
    let mut populations = Vec::with_capacity(capture.populations.len());
    for (i, population) in capture.populations.iter().enumerate() {
        let sorbent = population.sorbent();
        let loading = AdsorbentLoading::new(quotient(adsorbed[i], sorbent.mass().value())?);
        let mut conditions = sorbent.conditions().clone();
        for (axis, dims) in owned_axes() {
            if matches!(axis, TEMPERATURE_AXIS | VOLUME_AXIS) {
                continue;
            }
            conditions = conditions
                .with_quantity(
                    axis,
                    QuantitySpec::dimensional(dims),
                    transfer_conditions.axes()[axis],
                )
                .map_err(CyanidationError::from)?;
        }
        for (axis, dims, value) in [
            (
                COUNTERION_AXIS,
                AmountConcentration::DIMS,
                quotient(counterion.value(), volume.value())?,
            ),
            (
                OXYGEN_PARTIAL_PRESSURE_AXIS,
                Pressure::DIMS,
                boundary.pressure.value(),
            ),
            (
                SATURATION_CONCENTRATION_AXIS,
                AmountConcentration::DIMS,
                quotient(boundary.saturation.value(), volume.value())?,
            ),
        ] {
            conditions = conditions
                .with_quantity(axis, QuantitySpec::dimensional(dims), value)
                .map_err(CyanidationError::from)?;
        }
        let transfer =
            sorbent
                .model()
                .rate_at(concentration, loading, volume, temperature, &conditions)?;
        populations.push(CilSorbentState {
            id: sorbent.id().into(),
            role: population.role(),
            adsorbed: Amount::new(adsorbed[i]),
            vacant_capacity: Amount::new(vacant[i]),
            loading,
            transfer,
        });
    }
    Ok(CilFlowState {
        elapsed: solution.elapsed,
        solution: solution.clone(),
        populations,
        dissolved_counterion: counterion,
        oxygen_partial_pressure: boundary.pressure,
        oxygen_saturation: boundary.saturation,
        transfer_conditions,
    })
}

#[derive(Clone)]
struct CaptureTrial {
    gas: GasTrial,
    adsorbed: Vec<f64>,
    vacant: Vec<f64>,
}
impl CaptureTrial {
    fn values(&self) -> impl Iterator<Item = f64> + '_ {
        self.gas
            .values()
            .chain(self.adsorbed.iter().copied())
            .chain(self.vacant.iter().copied())
    }
}
struct CoupledProblem<'a> {
    gas: AeratedProblem<'a>,
    capture: &'a CaptureLayout,
}
impl CoupledProblem<'_> {
    fn admit(
        &self,
        state: &CaptureTrial,
        time: f64,
        work: &mut Work<'_, impl FnMut() -> bool>,
    ) -> Result<()> {
        carbon_state(
            self.capture,
            &evaluate(self.gas.liquid.layout, &state.gas.liquid.raw, time, work)?,
            &state.adsorbed,
            &state.vacant,
            &self.gas.boundary,
        )?;
        Ok(())
    }
    fn split(
        &self,
        mut state: CaptureTrial,
        start: f64,
        end: f64,
        work: &mut Work<'_, impl FnMut() -> bool>,
    ) -> Result<CaptureTrial> {
        let half = product(end - start, 0.5)?;
        {
            let mut admit = |solution: &CyanidationFlowState| {
                carbon_state(
                    self.capture,
                    solution,
                    &state.adsorbed,
                    &state.vacant,
                    &self.gas.boundary,
                )
                .map(|_| ())
            };
            self.gas.liquid.transport_admitted(
                &mut state.gas.liquid,
                half,
                start,
                work,
                &mut admit,
            )?;
        }
        self.gas_substep(&mut state.gas, half, start, work)?;
        self.admit(&state, start, work)?;
        self.adsorb(&mut state, half, start, work)?;
        {
            let mut admit = |solution: &CyanidationFlowState| {
                carbon_state(
                    self.capture,
                    solution,
                    &state.adsorbed,
                    &state.vacant,
                    &self.gas.boundary,
                )
                .map(|_| ())
            };
            self.gas.liquid.reaction_admitted(
                &mut state.gas.liquid,
                start,
                end,
                work,
                &mut admit,
            )?;
        }
        self.adsorb(&mut state, half, end, work)?;
        self.gas_substep(&mut state.gas, half, end, work)?;
        self.admit(&state, end, work)?;
        {
            let mut admit = |solution: &CyanidationFlowState| {
                carbon_state(
                    self.capture,
                    solution,
                    &state.adsorbed,
                    &state.vacant,
                    &self.gas.boundary,
                )
                .map(|_| ())
            };
            self.gas.liquid.transport_admitted(
                &mut state.gas.liquid,
                half,
                end,
                work,
                &mut admit,
            )?;
        }
        Ok(state)
    }

    fn gas_substep(
        &self,
        state: &mut GasTrial,
        h: f64,
        time: f64,
        work: &mut Work<'_, impl FnMut() -> bool>,
    ) -> Result<()> {
        // Neutral gas exchange leaves the liquid counterion unchanged. Bind its
        // actual value to the original source for this exact original substep.
        let liquid = &self.gas.liquid;
        let extra = counterion_extra(
            &self.gas.boundary,
            state.liquid.raw[self.capture.counterion_component],
            liquid.layout.charge.solution_volume.value(),
        )?;
        AeratedProblem {
            liquid: Problem {
                layout: liquid.layout,
                target: liquid.target,
                flow: liquid.flow,
            },
            boundary: Boundary {
                transfer: self.gas.boundary.transfer,
                pressure: self.gas.boundary.pressure,
                saturation: self.gas.boundary.saturation,
                extra: &extra,
            },
        }
        .gas(state, h, time, work)?;
        Ok(())
    }

    fn adsorb(
        &self,
        state: &mut CaptureTrial,
        h: f64,
        time: f64,
        work: &mut Work<'_, impl FnMut() -> bool>,
    ) -> Result<()> {
        let initial = SorptionStocks::from_trial(state, self.capture.counterion_component);
        // The local sorption law is autonomous. The surrounding split endpoint
        // owns the absolute source clock, as for the original exact gas substep.
        let rate = |stocks: &SorptionStocks, work: &mut Work<'_, _>| -> Result<SorptionStocks> {
            work.tick()?;
            let mut raw = state.gas.liquid.raw.clone();
            raw[3] = stocks.dissolved;
            raw[self.capture.counterion_component] = stocks.counterion;
            let solution = evaluate(self.gas.liquid.layout, &raw, time, work)?;
            let actual = carbon_state(
                self.capture,
                &solution,
                &stocks.adsorbed,
                &stocks.vacant,
                &self.gas.boundary,
            )?;
            let adsorbed = actual
                .populations
                .iter()
                .zip(&self.capture.populations)
                .map(|(s, p)| {
                    signed_product(s.transfer.loading_rate.value(), p.sorbent().mass().value())
                })
                .collect::<Result<Vec<_>>>()?;
            let total = sum(adsorbed.iter().copied())?;
            Ok(SorptionStocks {
                dissolved: -total,
                counterion: -total,
                vacant: adsorbed.iter().map(|r| -r).collect(),
                adsorbed,
            })
        };
        let first = rate(&initial, work)?;
        let second = rate(
            &initial.shifted(&first, h * 0.5, &self.capture.capacities)?,
            work,
        )?;
        let third = rate(
            &initial.shifted(&second, h * 0.5, &self.capture.capacities)?,
            work,
        )?;
        let fourth = rate(&initial.shifted(&third, h, &self.capture.capacities)?, work)?;
        let average = SorptionStocks::average([&first, &second, &third, &fourth])?;
        let result = initial.shifted(&average, h, &self.capture.capacities)?;
        state.gas.liquid.raw[3] = result.dissolved;
        state.gas.liquid.raw[self.capture.counterion_component] = result.counterion;
        state.adsorbed = result.adsorbed;
        state.vacant = result.vacant;
        self.admit(state, time, work)
    }
}

struct SorptionStocks {
    dissolved: f64,
    counterion: f64,
    adsorbed: Vec<f64>,
    vacant: Vec<f64>,
}
impl SorptionStocks {
    fn from_trial(state: &CaptureTrial, counterion: usize) -> Self {
        Self {
            dissolved: state.gas.liquid.raw[3],
            counterion: state.gas.liquid.raw[counterion],
            adsorbed: state.adsorbed.clone(),
            vacant: state.vacant.clone(),
        }
    }
    fn shifted(&self, rate: &Self, h: f64, capacities: &[f64]) -> Result<Self> {
        Ok(Self {
            dissolved: displaced(self.dissolved, rate.dissolved, h, f64::MAX)?,
            counterion: displaced(self.counterion, rate.counterion, h, f64::MAX)?,
            adsorbed: self
                .adsorbed
                .iter()
                .zip(&rate.adsorbed)
                .zip(capacities)
                .map(|((&a, &r), &cap)| displaced(a, r, h, cap))
                .collect::<Result<_>>()?,
            vacant: self
                .vacant
                .iter()
                .zip(&rate.vacant)
                .zip(capacities)
                .map(|((&v, &r), &cap)| displaced(v, r, h, cap))
                .collect::<Result<_>>()?,
        })
    }
    fn average(rates: [&Self; 4]) -> Result<Self> {
        let component = |values: [f64; 4]| -> Result<f64> {
            Ok(sum(values
                .into_iter()
                .zip([6.0, 3.0, 3.0, 6.0])
                .map(|(r, d)| signed_product(r, 1.0 / d))
                .collect::<Result<Vec<_>>>()?
                .into_iter())?)
        };
        Ok(Self {
            dissolved: component(rates.map(|r| r.dissolved))?,
            counterion: component(rates.map(|r| r.counterion))?,
            adsorbed: (0..rates[0].adsorbed.len())
                .map(|i| component(rates.map(|r| r.adsorbed[i])))
                .collect::<Result<_>>()?,
            vacant: (0..rates[0].vacant.len())
                .map(|i| component(rates.map(|r| r.vacant[i])))
                .collect::<Result<_>>()?,
        })
    }
}
fn signed_product(rate: f64, time: f64) -> Result<f64> {
    let magnitude = product(rate.abs(), time)?;
    Ok(if rate < 0.0 { -magnitude } else { magnitude })
}
fn displaced(initial: f64, rate: f64, h: f64, maximum: f64) -> Result<f64> {
    if rate == 0.0 {
        return Ok(initial);
    }
    let change = signed_product(rate, h)?;
    let next = sum([initial, change].into_iter())?;
    if next < 0.0 || next > maximum || (next == 0.0 && initial > 0.0) {
        return Err(CyanidationError::TrialDepletion.into());
    }
    if next == initial {
        return Err(CilError::Arithmetic(
            "nonzero carbon transfer is below inventory resolution",
        ));
    }
    Ok(next)
}
fn reserve(values: &[f64], counterion: usize) -> Result<()> {
    ensure(
        values[counterion] >= sum([values[0], values[3]].into_iter())?,
        "mobile counterion must cover accessible and dissolved gold in tank and feed",
    )
}

#[allow(clippy::too_many_arguments)]
fn accounts(
    layout: &FlowLayout,
    capture: &CaptureLayout,
    state: &CilFlowState,
    initial: &[f64],
    incoming: &[f64],
    outgoing: &[f64],
    gold: f64,
    gas_in: f64,
    gas_out: f64,
    initial_adsorbed: &[f64],
) -> Result<(FlowConservation, CilConservation)> {
    // Independently integrated bound stocks supply signed neutral-pair exchange.
    // Published inlet/outlet fields remain liquid-only; neither is a balance plug.
    let net = sum(state
        .populations
        .iter()
        .zip(initial_adsorbed)
        .map(|(p, old)| p.adsorbed.value() - old))?;
    let mut all_in = incoming.to_vec();
    let mut all_out = outgoing.to_vec();
    let account = if net >= 0.0 {
        &mut all_out
    } else {
        &mut all_in
    };
    for i in [3, capture.counterion_component] {
        account[i] = sum([account[i], net.abs()].into_iter())?;
    }
    let liquid = super::super::aeration::open_conservation(
        layout,
        &state.solution,
        initial,
        &all_in,
        &all_out,
        gold,
        gas_in,
        gas_out,
    )?;
    let mut maximum_site_residual = 0.0_f64;
    for (p, capacity) in state.populations.iter().zip(&capture.capacities) {
        let residual = sum([p.adsorbed.value(), p.vacant_capacity.value(), -capacity].into_iter())?;
        require_balance(residual, *capacity)?;
        maximum_site_residual = maximum_site_residual.max(residual.abs());
    }
    let coupled = CilConservation {
        gold: liquid.gold,
        cyanide: liquid.cyanide,
        counterion: liquid.components.spectators[capture.counterion_spectator].amount,
        dissolved_oxygen: liquid.dissolved_oxygen,
        maximum_site_residual: Amount::new(maximum_site_residual),
    };
    Ok((liquid, coupled))
}
fn retryable(error: &CilError) -> bool {
    matches!(
        error,
        CilError::Cyanidation(
            CyanidationError::TrialDepletion | CyanidationError::OutsideValidity { .. }
        ) | CilError::Adsorption(AdsorptionFlowError::Source(
            AqueousAdsorptionError::OutsideValidity { .. }
        ))
    )
}
fn ensure(condition: bool, reason: &'static str) -> Result<()> {
    if condition {
        Ok(())
    } else {
        Err(CilError::Invalid(reason))
    }
}

#[cfg(test)]
#[path = "coupled/tests.rs"]
mod tests;
