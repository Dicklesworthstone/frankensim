//! Continuous gold leaching with independently sourced oxygen uptake/stripping.
//!
//! The original tank owns liquid streams, local leach kinetics and simultaneous
//! aqueous chemistry. The original OxygenTransfer owns gas solubility, kLa and
//! its source domain. Its maintained reservoir adds kLa*(n_sat-n_O2) to the
//! dissolved-O2 balance; no gas supply capacity or aerator power is inferred.
//! The model reference and clean-water/slurry distinction are retained in
//! crate::cyanide_leach::aeration.
//!
//! Symmetric liquid transport / gas relaxation / chemical reaction splitting
//! retains positive finite stocks. Gross inward/outward gas amounts are evolved
//! independently, never reconstructed from a final balance. Raw full/two-half
//! differences control stocks, liquid exports, Au progress and both gas totals:
//! no smooth-order divisor is assumed across a directional-flux reversal.
//! These are local estimates, not global bounds or source/model uncertainty.

use super::*;
use crate::cyanide_leach::aeration::{
    OXYGEN_PARTIAL_PRESSURE_AXIS, OxygenTransfer, SATURATION_CONCENTRATION_AXIS,
};
use fs_qty::Pressure;

/// A continuous liquid tank with a separately declared maintained gas reservoir.
/// Attachment retains the tank's original charge, elapsed clock and all history.
#[derive(Clone, Debug)]
pub struct AeratedCyanidationCstr {
    tank: CyanidationCstr,
    transfer: Arc<OxygenTransfer>,
    oxygen_in: f64,
    oxygen_out: f64,
}

/// One atomic interval, with liquid streams and gas exchange reported separately.
#[derive(Clone, Debug)]
pub struct AeratedFlowAdvance {
    /// Original liquid/reaction report. Its oxygen/component conservation fields
    /// additionally include the separately published inward/outward gas amounts.
    pub flow: CyanidationFlowAdvance,
    pub oxygen_partial_pressure: Pressure,
    pub oxygen_saturation: Amount,
    /// Pinned measured-Hcp receipt at the actual held temperature, when selected.
    pub henry_answer: Option<fs_matdb::MaterialAnswer>,
    pub oxygen_in: Amount,
    pub oxygen_out: Amount,
    pub cumulative_oxygen_in: Amount,
    pub cumulative_oxygen_out: Amount,
    /// Actual final state admitted by the original transfer source, including clock.
    pub transfer_conditions: QueryPoint,
    transfer: Arc<OxygenTransfer>,
}
impl AeratedFlowAdvance {
    pub fn transfer(&self) -> &OxygenTransfer {
        &self.transfer
    }
}
impl AeratedCyanidationCstr {
    /// Attach to an admitted tank, including one already advanced without gas.
    /// Earlier liquid/reaction history remains part of the cumulative balance;
    /// the newly attached gas boundary starts with zero inward/outward transfer.
    pub fn from_tank(tank: CyanidationCstr, transfer: OxygenTransfer) -> Result<Self> {
        let layout = &tank.layout;
        let model = &layout.model;
        let external_gas = transfer.external_gas();
        require(
            model.chemistry.families().iter().all(|family| {
                family
                    .species()
                    .iter()
                    .all(|species| species.as_str() != external_gas)
            }) && layout
                .charge
                .spectators
                .iter()
                .all(|ion| ion.species.as_str() != external_gas)
                && layout.competition.as_ref().is_none_or(|metals| {
                    metals.model.families().iter().all(|family| {
                        family
                            .species()
                            .iter()
                            .all(|species| species.as_str() != external_gas)
                    })
                }),
            "external oxygen cannot also be an owned liquid family or spectator",
        )?;
        require(
            [OXYGEN_PARTIAL_PRESSURE_AXIS, SATURATION_CONCENTRATION_AXIS]
                .iter()
                .all(|axis| {
                    !layout.conditions.axes().contains_key(*axis)
                        && model.source.validity.bound(axis).is_none()
                }),
            "gas pressure and saturation belong to the transfer source",
        )?;
        let water = model.chemistry.water().standard();
        transfer.admit_chemistry(water)?;
        let aqueous = layout
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
                .chain([model.reaction.clone(), transfer.reaction()?])
                .collect(),
        )?;
        Ok(Self {
            tank,
            transfer: Arc::new(transfer),
            oxygen_in: 0.0,
            oxygen_out: 0.0,
        })
    }
    /// Read-only liquid owner; direct advancement cannot bypass gas accounting.
    pub fn tank(&self) -> &CyanidationCstr {
        &self.tank
    }
    pub fn transfer(&self) -> &OxygenTransfer {
        &self.transfer
    }
    pub fn oxygen_in(&self) -> Amount {
        Amount::new(self.oxygen_in)
    }
    pub fn oxygen_out(&self) -> Amount {
        Amount::new(self.oxygen_out)
    }

    /// Advance at constant feed, flow and prescribed oxygen partial pressure.
    /// All three may change between calls. Temperature, volume and the original
    /// activity basis stay fixed. A zero-flow interval still permits gas exchange
    /// and dissolution; a zero-kLa interval still permits liquid transport.
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
    ) -> Result<AeratedFlowAdvance> {
        nonnegative(
            duration.value(),
            "aerated-flow duration must be nonnegative",
        )?;
        nonnegative(flow.value(), "aerated slurry flow must be nonnegative")?;
        let tank = &self.tank;
        let layout = &tank.layout;
        let start = tank.elapsed().value();
        let end = add(start, duration.value())?;
        let (saturation, henry_answer) = self.transfer.equilibrium_saturation(
            oxygen_partial_pressure,
            layout.charge.amount_at_unit_activity,
            &tank.state.conditions,
            transfer_conditions,
        )?;
        self.transfer
            .validate_extras(&layout.conditions, transfer_conditions)?;
        let boundary = Boundary {
            transfer: &self.transfer,
            pressure: oxygen_partial_pressure,
            saturation,
            extra: transfer_conditions,
        };
        let mut work = dosing::work(numerics.maximum_work, &mut cancelled)?;
        boundary.admit(&tank.state)?;
        let (minimum, maximum) = self
            .transfer
            .source()
            .validity
            .bound(ELAPSED_TIME_AXIS)
            .expect("original transfer clock domain");
        if end < minimum || end > maximum {
            return Err(CyanidationError::OutsideValidity {
                axis: ELAPSED_TIME_AXIS.into(),
                value: end,
                minimum,
                maximum,
            });
        }
        let concentration = feed_values(layout, feed, &mut work)?;
        let volume = layout.charge.solution_volume.value();
        let exposure = quotient(product(flow.value(), duration.value())?, volume)?;
        let target = concentration
            .iter()
            .map(|c| product(*c, volume))
            .collect::<Result<Vec<_>>>()?;
        let incoming = target
            .iter()
            .map(|n| product(*n, exposure))
            .collect::<Result<Vec<_>>>()?;
        let gas_scale = product(
            product(self.transfer.coefficient().value(), duration.value())?,
            tank.raw[2].max(saturation.value()).max(target[2]),
        )?;
        let scale = sum(tank
            .raw
            .iter()
            .chain(&incoming)
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
        let liquid = Problem {
            layout,
            target: &target,
            flow,
        };
        let problem = AeratedProblem { liquid, boundary };
        let mut current = GasTrial::new(tank.raw.clone());
        let (mut clock, mut step) = (start, numerics.maximum_step.value().min(duration.value()));
        let (mut accepted, mut rejected) = (0, 0);
        let (mut error_sum, mut largest) = (0.0_f64, 0.0_f64);
        while clock < end {
            work.tick()?;
            if accepted + rejected >= numerics.maximum_attempts {
                return Err(CyanidationError::Budget);
            }
            let finish =
                if clock + step >= end || end - (clock + step) <= 8.0 * f64::EPSILON * end.abs() {
                    end
                } else {
                    clock + step
                };
            let middle = clock.midpoint(finish);
            let h = finish - clock;
            if h <= 0.0 || middle <= clock || middle >= finish {
                return Err(CyanidationError::Arithmetic(
                    "aerated-flow clock loses resolution",
                ));
            }
            // Reset interval counters for each error trial: large historical
            // gas or liquid totals must not hide a small current increment.
            let initial = GasTrial::new(current.liquid.raw.clone());
            let trial = (|| {
                let coarse = problem.split(initial.clone(), clock, finish, &mut work)?;
                let first = problem.split(initial, clock, middle, &mut work)?;
                let fine = problem.split(first, middle, finish, &mut work)?;
                Ok((coarse, fine))
            })();
            let (coarse, fine) = match trial {
                Ok(pair) => pair,
                Err(
                    CyanidationError::TrialDepletion | CyanidationError::OutsideValidity { .. },
                ) => {
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
            current.liquid.raw = fine.liquid.raw;
            for (total, increment) in current.liquid.outlet.iter_mut().zip(fine.liquid.outlet) {
                *total = add(*total, increment)?;
            }
            current.liquid.gold = add(current.liquid.gold, fine.liquid.gold)?;
            current.inward = add(current.inward, fine.inward)?;
            current.outward = add(current.outward, fine.outward)?;
            error_sum = next_error;
            largest = largest.max(error);
            accepted += 1;
            clock = finish;
            step = if error <= allowance / 8.0 {
                (h * 2.0).min(numerics.maximum_step.value())
            } else {
                h
            };
        }
        let state = evaluate(layout, &current.liquid.raw, end, &mut work)?;
        let actual_transfer_conditions = problem.boundary.admit(&state)?;
        let cumulative_inlet = tank
            .inlet
            .iter()
            .zip(&incoming)
            .map(|(a, b)| add(*a, *b))
            .collect::<Result<Vec<_>>>()?;
        let cumulative_outlet = tank
            .outlet
            .iter()
            .zip(&current.liquid.outlet)
            .map(|(a, b)| add(*a, *b))
            .collect::<Result<Vec<_>>>()?;
        let gold = add(tank.gold_dissolved, current.liquid.gold)?;
        let oxygen_in = add(self.oxygen_in, current.inward)?;
        let oxygen_out = add(self.oxygen_out, current.outward)?;
        let accumulated_error = sum([tank.accumulated_error, error_sum].into_iter())?;
        let interval_conservation = open_conservation(
            layout,
            &state,
            &tank.raw,
            &incoming,
            &current.liquid.outlet,
            current.liquid.gold,
            current.inward,
            current.outward,
        )?;
        let cumulative_conservation = open_conservation(
            layout,
            &state,
            &layout.initial,
            &cumulative_inlet,
            &cumulative_outlet,
            gold,
            oxygen_in,
            oxygen_out,
        )?;
        let stoichiometric_amounts = layout
            .model
            .reaction
            .amounts(Amount::new(quotient(gold, 4.0)?))?;
        let mut answer = AeratedFlowAdvance {
            flow: CyanidationFlowAdvance {
                layout: layout.clone(),
                previous_elapsed: Time::new(start),
                elapsed: Time::new(end),
                flow,
                feed: feed.clone(),
                state: state.clone(),
                inlet: inventory(layout, &incoming),
                outlet: inventory(layout, &current.liquid.outlet),
                cumulative_inlet: inventory(layout, &cumulative_inlet),
                cumulative_outlet: inventory(layout, &cumulative_outlet),
                interval_gold_dissolved: Amount::new(current.liquid.gold),
                cumulative_gold_dissolved: Amount::new(gold),
                stoichiometric_amounts,
                interval_conservation,
                cumulative_conservation,
                maximum_local_error: Amount::new(largest),
                accumulated_error_estimate: Amount::new(accumulated_error),
                accepted_steps: accepted,
                rejected_steps: rejected,
                work: 0,
            },
            oxygen_partial_pressure,
            oxygen_saturation: saturation,
            henry_answer,
            oxygen_in: Amount::new(current.inward),
            oxygen_out: Amount::new(current.outward),
            cumulative_oxygen_in: Amount::new(oxygen_in),
            cumulative_oxygen_out: Amount::new(oxygen_out),
            transfer_conditions: actual_transfer_conditions,
            transfer: self.transfer.clone(),
        };
        work.tick()?;
        answer.flow.work = work.used;
        self.tank.raw = current.liquid.raw;
        self.tank.state = state;
        self.tank.inlet = cumulative_inlet;
        self.tank.outlet = cumulative_outlet;
        self.tank.gold_dissolved = gold;
        self.tank.accumulated_error = accumulated_error;
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
    fn admit(&self, state: &CyanidationFlowState) -> Result<QueryPoint> {
        // The transfer source owns the absolute process clock. It does not turn
        // the local leach card into a residence-time or batch-age correlation.
        let point = state.conditions.clone().with_quantity(
            ELAPSED_TIME_AXIS,
            QuantitySpec::dimensional(Time::DIMS),
            state.elapsed.value(),
        )?;
        self.transfer
            .rate_at_conditions(
                &point,
                state.inventory.dissolved_oxygen,
                self.pressure,
                self.saturation,
                self.extra,
            )
            .map(|(_, point)| point)
    }
}
#[derive(Clone)]
struct GasTrial {
    liquid: Trial,
    inward: f64,
    outward: f64,
}
impl GasTrial {
    fn new(raw: Vec<f64>) -> Self {
        Self {
            liquid: Trial {
                outlet: vec![0.0; raw.len()],
                raw,
                gold: 0.0,
            },
            inward: 0.0,
            outward: 0.0,
        }
    }
    fn values(&self) -> impl Iterator<Item = f64> + '_ {
        self.liquid.values().chain([self.inward, self.outward])
    }
}
struct AeratedProblem<'a> {
    liquid: Problem<'a>,
    boundary: Boundary<'a>,
}
impl AeratedProblem<'_> {
    fn gas(
        &self,
        state: &mut GasTrial,
        h: f64,
        time: f64,
        work: &mut Work<'_, impl FnMut() -> bool>,
    ) -> Result<()> {
        work.tick()?;
        let before = state.liquid.raw[2];
        let saturation = self.boundary.saturation.value();
        let exposure = product(self.boundary.transfer.coefficient().value(), h)?;
        if exposure > 0.0 && before != saturation {
            let e = det::exp(-exposure);
            let u = -det::expm1(-exposure);
            let magnitude = product((saturation - before).abs(), u)?;
            let after = sum([product(before, e)?, product(saturation, u)?].into_iter())?;
            if after == before || (before > 0.0 && saturation == 0.0 && e == 0.0) {
                return Err(CyanidationError::Arithmetic(
                    "finite gas exchange loses inventory resolution",
                ));
            }
            state.liquid.raw[2] = after;
            if saturation > before {
                state.inward = add(state.inward, magnitude)?;
            } else {
                state.outward = add(state.outward, magnitude)?;
            }
        }
        self.boundary.admit(&evaluate(
            self.liquid.layout,
            &state.liquid.raw,
            time,
            work,
        )?)?;
        Ok(())
    }
    fn split(
        &self,
        mut state: GasTrial,
        start: f64,
        end: f64,
        work: &mut Work<'_, impl FnMut() -> bool>,
    ) -> Result<GasTrial> {
        let half = product(end - start, 0.5)?;
        let mut admit = |solution: &CyanidationFlowState| self.boundary.admit(solution).map(|_| ());
        self.liquid
            .transport_admitted(&mut state.liquid, half, start, work, &mut admit)?;
        self.gas(&mut state, half, start, work)?;
        self.liquid
            .reaction_admitted(&mut state.liquid, start, end, work, &mut admit)?;
        self.gas(&mut state, half, end, work)?;
        self.liquid
            .transport_admitted(&mut state.liquid, half, end, work, &mut admit)?;
        Ok(state)
    }
}
#[allow(clippy::too_many_arguments)]
fn open_conservation(
    layout: &FlowLayout,
    state: &CyanidationFlowState,
    initial: &[f64],
    incoming: &[f64],
    outgoing: &[f64],
    gold: f64,
    gas_in: f64,
    gas_out: f64,
) -> Result<FlowConservation> {
    // The original liquid arrays and charge remain untouched. Both independently
    // integrated gas accounts participate explicitly in the combined O2 balance.
    let mut all_in = incoming.to_vec();
    let mut all_out = outgoing.to_vec();
    all_in[2] = sum([all_in[2], gas_in].into_iter())?;
    all_out[2] = sum([all_out[2], gas_out].into_iter())?;
    conservation(layout, state, initial, &all_in, &all_out, gold)
}

#[cfg(test)]
#[path = "aeration/tests.rs"]
mod tests;
