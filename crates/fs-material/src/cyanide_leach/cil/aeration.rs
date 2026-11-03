//! Maintained-pressure oxygen exchange during simultaneous carbon-in-leach.
//!
//! This is the direct coupling of the existing CIL owner with the existing
//! prescribed-pressure oxygen-transfer owner. Gold dissolution, shared aqueous
//! chemistry, reversible product/preg-robbing carbon loading and gas exchange
//! advance in one RK4 state. Oxygen supplied by gas is therefore available to
//! the same rate stages that determine adsorption and pH/free-CN chemistry.
//!
//! The gas reservoir is maintained at a caller-prescribed O2 partial pressure;
//! it is not a finite headspace, blower or power model. Carbon remains retained
//! in this batch owner. Continuous slurry flow and inter-tank carbon movement
//! remain separate physical balances.

use super::*;
use crate::cyanide_leach::aeration::{
    OXYGEN_PARTIAL_PRESSURE_AXIS, OxygenTransfer, SATURATION_CONCENTRATION_AXIS,
};
use fs_qty::Pressure;

/// Simultaneous CIL with the original maintained-pressure oxygen-transfer owner.
#[derive(Clone, Debug)]
pub struct AeratedCarbonInLeachBatch {
    batch: CarbonInLeachBatch,
    transfer: Arc<OxygenTransfer>,
    oxygen_in: f64,
    oxygen_out: f64,
    accumulated_error: f64,
}

/// One atomic accepted CIL interval plus independent gas-transfer accounts.
#[derive(Clone, Debug)]
pub struct AeratedCilAdvance {
    /// Existing simultaneous leach/carbon report; its oxygen conservation field
    /// includes the separately reported gas boundary below.
    pub cil: CilAdvance,
    pub oxygen_partial_pressure: Pressure,
    pub oxygen_saturation: Amount,
    pub interval_oxygen_in: Amount,
    pub interval_oxygen_out: Amount,
    pub cumulative_oxygen_in: Amount,
    pub cumulative_oxygen_out: Amount,
    pub net_oxygen_transfer: Amount,
    /// Actual endpoint admitted by the independent transfer source.
    pub transfer_conditions: QueryPoint,
    transfer: Arc<OxygenTransfer>,
}
impl AeratedCilAdvance {
    pub fn transfer(&self) -> &OxygenTransfer {
        &self.transfer
    }
}

impl AeratedCarbonInLeachBatch {
    /// Attach gas exchange to freshly admitted simultaneous CIL.
    ///
    /// The closed-bath oxygen ceiling is replaced by the original Au/CN ceiling;
    /// gas exchange starts at zero. Carbon inventories, preloading, reagent
    /// history, competing-metal chemistry and every original source are retained.
    pub fn from_batch(closed: CarbonInLeachBatch, transfer: OxygenTransfer) -> Result<Self> {
        ensure(
            closed.state.elapsed.value() == 0.0
                && closed.state.solution.gold_dissolved.value() == 0.0,
            "aerated CIL attachment requires a fresh zero-progress bath",
        )?;
        let external_gas = transfer.solubility().gas();
        let chemical = &closed.layout.leach;
        let model = &chemical.model;
        ensure(
            model.chemistry.families().iter().all(|family| {
                family.species().iter().all(|species| species != external_gas)
            }) && chemical
                .charge
                .spectators
                .iter()
                .all(|ion| &ion.species != external_gas)
                && chemical.competition.as_ref().is_none_or(|metals| {
                    metals.model.families().iter().all(|family| {
                        family.species().iter().all(|species| species != external_gas)
                    })
                })
                && closed.layout.pair.counterion() != external_gas
                && closed.layout.pair.adsorbed_pair() != external_gas,
            "external oxygen cannot also belong to the retained liquid or sorbed chemistry",
        )?;
        ensure(
            [OXYGEN_PARTIAL_PRESSURE_AXIS, SATURATION_CONCENTRATION_AXIS]
                .iter()
                .all(|axis| {
                    !chemical.conditions.axes().contains_key(*axis)
                        && model.source.validity.bound(axis).is_none()
                }),
            "gas pressure and saturation belong to the transfer source",
        )?;
        let gas = transfer.solubility().standard();
        let water = model.chemistry.water().standard();
        ensure(
            gas.temperature() == water.temperature()
                && gas.reference_pressure() == water.reference_pressure()
                && gas.reference_axis() == water.reference_axis(),
            "gas and aqueous sources need identical temperature and reference conventions",
        )?;
        let aqueous = chemical
            .competition
            .as_ref()
            .map_or(model.chemistry.reaction_network(), |metals| {
                metals.model.reaction_network()
            });
        ReactionNetwork::new(
            aqueous
                .reactions()
                .iter()
                .cloned()
                .chain([
                    model.reaction.clone(),
                    closed.layout.pair.reaction().clone(),
                    transfer.solubility().reaction().clone(),
                ])
                .collect(),
        )
        .map_err(CyanidationError::from)?;

        let mut layout = (*closed.layout).clone();
        let mut leach = (*layout.leach).clone();
        leach.limit = leach
            .charge
            .accessible_gold
            .value()
            .min(product(0.5, dosing::cyanide_total(&leach)?)?);
        // The closed CIL attachment may previously have admitted a smaller
        // counterion reserve because finite liquid O2 limited attainable Au
        // dissolution. Opening oxygen removes that ceiling, so re-admit the
        // neutral-pair counterion against the expanded Au/CN progress limit.
        let available_counterion = closed.state.dissolved_counterion.value();
        let required_counterion =
            sum([leach.charge.gold_complex.value(), leach.limit].into_iter())?;
        ensure(
            available_counterion >= required_counterion,
            "insufficient counterion reserve for aerated CIL gold loading",
        )?;
        layout.leach = Arc::new(leach);
        let layout = Arc::new(layout);

        let mut never = || false;
        let mut work = Work {
            used: 0,
            maximum: 10_000_000,
            cancelled: &mut never,
            stopped: None,
        };
        let state = evaluate_cil(&layout, 0.0, &Stocks::from_state(&closed.state), &mut work)?;
        let batch = CarbonInLeachBatch {
            layout,
            state,
            accumulated_error: closed.accumulated_error,
        };
        open_conservation(&batch.layout, &batch.state, 0.0, 0.0)?;
        Ok(Self {
            accumulated_error: batch.accumulated_error,
            batch,
            transfer: Arc::new(transfer),
            oxygen_in: 0.0,
            oxygen_out: 0.0,
        })
    }

    pub fn state(&self) -> &CilState {
        self.batch.state()
    }
    pub fn model(&self) -> &GoldCyanidation {
        self.batch.model()
    }
    pub fn transfer(&self) -> &OxygenTransfer {
        &self.transfer
    }
    pub fn ion_pair(&self) -> &IonPairTransfer {
        self.batch.ion_pair()
    }
    pub fn models(&self) -> &[CarbonPopulation] {
        self.batch.models()
    }
    pub fn initial_charge(&self) -> &CyanidationCharge {
        self.batch.initial_charge()
    }
    pub fn reagent_inputs(&self) -> dosing::ReagentInputs {
        self.batch.reagent_inputs()
    }
    pub fn oxygen_in(&self) -> Amount {
        Amount::new(self.oxygen_in)
    }
    pub fn oxygen_out(&self) -> Amount {
        Amount::new(self.oxygen_out)
    }

    /// Add dissolved NaCN or NaOH without resetting gas-transfer history.
    ///
    /// No gas pressure is retained between advances, so the next advance admits
    /// its own pressure and transfer context against the newly re-speciated CIL
    /// state. The gas source and cumulative in/out accounts remain unchanged.
    pub fn dose_reagent(
        &mut self,
        reagent: dosing::DissolvedSodiumReagent,
        amount: Amount,
        maximum_work: usize,
        cancelled: impl FnMut() -> bool,
    ) -> Result<CilDoseReceipt> {
        self.batch
            .dose_reagent(reagent, amount, maximum_work, cancelled)
    }

    /// Advance at fixed prescribed O2 partial pressure. A later call may change
    /// pressure. Source refusal, work exhaustion and cancellation publish no
    /// partial chemistry, carbon inventory, gas account or clock.
    pub fn advance(
        &mut self,
        duration: Time,
        oxygen_partial_pressure: Pressure,
        transfer_conditions: &QueryPoint,
        numerics: CilNumerics,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AeratedCilAdvance> {
        nonnegative(
            duration.value(),
            "aerated CIL duration must be finite and nonnegative",
        )?;
        let start = self.batch.state.elapsed.value();
        let end = sum([start, duration.value()].into_iter())?;
        ensure(
            duration.value() == 0.0 || end > start,
            "aerated CIL duration is below clock resolution",
        )?;
        let saturation = self.transfer.solubility().equilibrium_aqueous_amount(
            oxygen_partial_pressure,
            self.batch.layout.leach.charge.amount_at_unit_activity,
        )?;
        self.transfer
            .validate_extras(&self.batch.layout.leach.conditions, transfer_conditions)?;
        check_clock(&self.batch.layout.leach.model.source, end)?;
        check_clock(self.transfer.source(), end)?;

        let scale = sum(
            [
                self.batch.layout.total_gold,
                sum(self.batch.layout.capacities.iter().copied())?,
                self.batch.state.solution.dissolved_oxygen.value(),
                saturation.value(),
            ]
            .into_iter(),
        )?;
        let tolerance = CyanidationNumerics {
            maximum_step: numerics.maximum_step,
            absolute_gold_error: numerics.absolute_amount_error,
            relative_gold_error: numerics.relative_amount_error,
            maximum_attempts: numerics.maximum_attempts,
            maximum_work: numerics.maximum_work,
        }
        .tolerance(scale)?;

        let boundary = Boundary {
            transfer: &self.transfer,
            pressure: oxygen_partial_pressure,
            saturation,
            extra: transfer_conditions,
        };
        let mut work = Work {
            used: 0,
            maximum: numerics.maximum_work,
            cancelled: &mut cancelled,
            stopped: None,
        };
        work.tick()?;
        boundary.rate(&self.batch.state)?;
        let initial_state = self.batch.state.clone();
        let mut current = Sample {
            state: initial_state.clone(),
            inward: 0.0,
            outward: 0.0,
        };
        let mut step = numerics.maximum_step.value().min(duration.value());
        let (mut accepted, mut rejected) = (0usize, 0usize);
        let mut largest = 0.0_f64;
        let mut accumulated = self.accumulated_error;

        while current.state.elapsed.value() < end {
            work.tick()?;
            if accepted + rejected >= numerics.maximum_attempts {
                return Err(CilError::Budget);
            }
            let remaining = end - current.state.elapsed.value();
            step = step.min(remaining);
            let finish = if step == remaining {
                end
            } else {
                current.state.elapsed.value() + step
            };
            let middle = current.state.elapsed.value() + 0.5 * (finish - current.state.elapsed.value());
            if middle <= current.state.elapsed.value() || finish <= middle {
                return Err(CilError::Arithmetic(
                    "positive aerated CIL substep is below clock resolution",
                ));
            }
            let trial = (|| -> Result<_> {
                let full = rk4(&self.batch.layout, &boundary, &current, finish, &mut work)?;
                let first = rk4(&self.batch.layout, &boundary, &current, middle, &mut work)?;
                let refined = rk4(&self.batch.layout, &boundary, &first, finish, &mut work)?;
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
            let chemistry_error =
                Stocks::from_state(&refined.state).error(&Stocks::from_state(&full.state))?;
            let gas_error = sum(
                [
                    (refined.inward - full.inward).abs(),
                    (refined.outward - full.outward).abs(),
                ]
                .into_iter(),
            )?;
            let error = sum([chemistry_error, gas_error].into_iter())?;
            let local = product(
                tolerance,
                quotient(finish - current.state.elapsed.value(), duration.value())?,
            )?;
            let next_accumulated = sum([accumulated, error].into_iter())?;
            if error > local || next_accumulated > tolerance {
                rejected += 1;
                step = product(step, 0.5)?;
                continue;
            }
            let trial_in = add_positive(self.oxygen_in, refined.inward)?;
            let trial_out = add_positive(self.oxygen_out, refined.outward)?;
            open_conservation(&self.batch.layout, &refined.state, trial_in, trial_out)?;
            current = refined;
            accepted += 1;
            largest = largest.max(error);
            accumulated = next_accumulated;
            if error < local / 32.0 {
                step = (step * 2.0).min(numerics.maximum_step.value());
            }
        }

        let (_, final_point) = boundary.rate(&current.state)?;
        let oxygen_in = add_positive(self.oxygen_in, current.inward)?;
        let oxygen_out = add_positive(self.oxygen_out, current.outward)?;
        let conservation =
            open_conservation(&self.batch.layout, &current.state, oxygen_in, oxygen_out)?;
        let capture = |role| -> Result<Amount> {
            Ok(Amount::new(sum(
                current
                    .state
                    .populations
                    .iter()
                    .zip(&initial_state.populations)
                    .filter(|(p, _)| p.role == role)
                    .map(|(after, before)| after.adsorbed.value() - before.adsorbed.value()),
            )?))
        };
        let mut cil = CilAdvance {
            layout: self.batch.layout.clone(),
            state: current.state.clone(),
            interval_gold_dissolved: Amount::new(
                current.state.solution.gold_dissolved.value()
                    - initial_state.solution.gold_dissolved.value(),
            ),
            interval_product_capture: capture(CarbonRole::ProductCarbon)?,
            interval_preg_robbing_capture: capture(CarbonRole::PregRobbing)?,
            conservation,
            accepted_steps: accepted,
            rejected_steps: rejected,
            work: 0,
            maximum_local_error: Amount::new(largest),
            accumulated_error_estimate: Amount::new(accumulated),
        };
        work.tick()?;
        cil.work = work.used;
        let answer = AeratedCilAdvance {
            cil,
            oxygen_partial_pressure,
            oxygen_saturation: saturation,
            interval_oxygen_in: Amount::new(current.inward),
            interval_oxygen_out: Amount::new(current.outward),
            cumulative_oxygen_in: Amount::new(oxygen_in),
            cumulative_oxygen_out: Amount::new(oxygen_out),
            net_oxygen_transfer: Amount::new(sum([oxygen_in, -oxygen_out].into_iter())?),
            transfer_conditions: final_point,
            transfer: self.transfer.clone(),
        };
        self.batch.state = current.state;
        self.batch.accumulated_error = accumulated;
        self.accumulated_error = accumulated;
        self.oxygen_in = oxygen_in;
        self.oxygen_out = oxygen_out;
        Ok(answer)
    }
}

struct Boundary<'a> {
    transfer: &'a OxygenTransfer,
    pressure: Pressure,
    saturation: Amount,
    extra: &'a QueryPoint,
}
impl Boundary<'_> {
    fn rate(&self, state: &CilState) -> Result<(f64, QueryPoint)> {
        self.transfer
            .rate_at_conditions(
                &state.solution.conditions,
                state.solution.dissolved_oxygen,
                self.pressure,
                self.saturation,
                self.extra,
            )
            .map_err(CilError::from)
    }
}

#[derive(Clone)]
struct Sample {
    state: CilState,
    inward: f64,
    outward: f64,
}
#[derive(Clone)]
struct Rates {
    cil: Stocks,
    inward: f64,
    outward: f64,
}
fn rates(layout: &CilLayout, boundary: &Boundary<'_>, sample: &Sample) -> Result<Rates> {
    let mut cil = derivatives(layout, &sample.state)?;
    let (flux, _) = boundary.rate(&sample.state)?;
    cil.oxygen = sum([cil.oxygen, flux].into_iter())?;
    Ok(Rates {
        cil,
        inward: flux.max(0.0),
        outward: (-flux).max(0.0),
    })
}
fn displaced(
    layout: &CilLayout,
    boundary: &Boundary<'_>,
    initial: &Sample,
    elapsed: f64,
    rate: &Rates,
    work: &mut Work<'_, impl FnMut() -> bool>,
) -> Result<Sample> {
    let dt = elapsed - initial.state.elapsed.value();
    let stocks = Stocks::from_state(&initial.state).shifted(&rate.cil, dt, layout)?;
    let state = evaluate_cil(layout, elapsed, &stocks, work)?;
    let sample = Sample {
        state,
        inward: add_positive(initial.inward, product(dt, rate.inward)?)?,
        outward: add_positive(initial.outward, product(dt, rate.outward)?)?,
    };
    boundary.rate(&sample.state)?;
    Ok(sample)
}
fn average(rates: [&Rates; 4]) -> Result<Rates> {
    let cil = Stocks::average([
        &rates[0].cil,
        &rates[1].cil,
        &rates[2].cil,
        &rates[3].cil,
    ])?;
    let component = |values: [f64; 4]| -> Result<f64> {
        Ok(sum(
            values
                .into_iter()
                .zip([6.0, 3.0, 3.0, 6.0])
                .map(|(r, d)| signed_scale(r, 1.0 / d))
                .collect::<Result<Vec<_>>>()?
                .into_iter(),
        )?)
    };
    Ok(Rates {
        cil,
        inward: component(rates.map(|r| r.inward))?,
        outward: component(rates.map(|r| r.outward))?,
    })
}
fn rk4(
    layout: &CilLayout,
    boundary: &Boundary<'_>,
    initial: &Sample,
    end: f64,
    work: &mut Work<'_, impl FnMut() -> bool>,
) -> Result<Sample> {
    let duration = end - initial.state.elapsed.value();
    let half = product(duration, 0.5)?;
    let middle = initial.state.elapsed.value() + half;
    let first = rates(layout, boundary, initial)?;
    let second_state = displaced(layout, boundary, initial, middle, &first, work)?;
    let second = rates(layout, boundary, &second_state)?;
    let third_state = displaced(layout, boundary, initial, middle, &second, work)?;
    let third = rates(layout, boundary, &third_state)?;
    let fourth_state = displaced(layout, boundary, initial, end, &third, work)?;
    let fourth = rates(layout, boundary, &fourth_state)?;
    displaced(
        layout,
        boundary,
        initial,
        end,
        &average([&first, &second, &third, &fourth])?,
        work,
    )
}
fn open_conservation(
    layout: &CilLayout,
    state: &CilState,
    oxygen_in: f64,
    oxygen_out: f64,
) -> Result<CilConservation> {
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
        oxygen_out,
        -layout.leach.charge.dissolved_oxygen.value(),
        -oxygen_in,
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
    ] {
        check_balance(residual, scale)?;
    }
    let oxygen_scale = sum(
        [
            layout.leach.charge.dissolved_oxygen.value(),
            quotient(solution.gold_dissolved.value(), 4.0)?,
            oxygen_in,
            oxygen_out,
        ]
        .into_iter(),
    )?;
    if oxygen.abs() > 2048.0 * f64::EPSILON * oxygen_scale.max(1.0) {
        return Err(CilError::Arithmetic(
            "open-boundary CIL oxygen balance does not close",
        ));
    }
    Ok(CilConservation {
        gold: Amount::new(gold),
        cyanide: Amount::new(cyanide),
        counterion: Amount::new(counterion),
        dissolved_oxygen: Amount::new(oxygen),
        maximum_site_residual: Amount::new(maximum_site_residual),
    })
}
fn add_positive(initial: f64, increment: f64) -> Result<f64> {
    let next = sum([initial, increment].into_iter())?;
    if increment > 0.0 && next <= initial {
        Err(CilError::Arithmetic(
            "positive oxygen transfer is below numeric resolution",
        ))
    } else {
        Ok(next)
    }
}
fn check_clock(source: &ConstitutiveModelCard, end: f64) -> Result<()> {
    let (minimum, maximum) = source
        .validity
        .bound(ELAPSED_TIME_AXIS)
        .ok_or(CilError::Invalid("source is missing the absolute clock"))?;
    if end < minimum || end > maximum {
        Err(CyanidationError::OutsideValidity {
            axis: ELAPSED_TIME_AXIS.into(),
            value: end,
            minimum,
            maximum,
        }
        .into())
    } else {
        Ok(())
    }
}

#[cfg(test)]
#[path = "aeration/tests.rs"]
mod tests;
