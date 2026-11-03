//! Prescribed-pressure oxygen exchange for continuous retained-carbon CIL.
//!
//! This owner couples three already-declared pieces without inventing new
//! constitutive data: the continuous cyanidation CSTR, retained product/preg-
//! robbing carbon populations, and the maintained-pressure O2 transfer owner.
//! Exact gas relaxation is Strang-split around one continuous-CIL advance;
//! full-vs-two-half trials bound the additional coupling error.
//!
//! The gas reservoir remains prescribed-pressure and infinite. Carbon remains
//! retained in this tank. No blower power, finite headspace, screen loss or
//! inter-tank carbon motion is inferred.

use super::*;
use crate::cyanide_leach::aeration::{
    OXYGEN_PARTIAL_PRESSURE_AXIS, OxygenTransfer, SATURATION_CONCENTRATION_AXIS,
};
use fs_math::det;
use fs_qty::Pressure;

/// Continuous CIL with a separately sourced maintained O2 reservoir.
#[derive(Clone, Debug)]
pub struct AeratedCarbonInLeachCstr {
    cil: CarbonInLeachCstr,
    transfer: Arc<OxygenTransfer>,
    oxygen_in: f64,
    oxygen_out: f64,
    coupling_error: f64,
}

/// One accepted interval with liquid, retained-carbon and gas accounts.
#[derive(Clone, Debug)]
pub struct AeratedFlowCilAdvance {
    pub cil: FlowCilAdvance,
    pub oxygen_partial_pressure: Pressure,
    pub oxygen_saturation: Amount,
    pub interval_oxygen_in: Amount,
    pub interval_oxygen_out: Amount,
    pub cumulative_oxygen_in: Amount,
    pub cumulative_oxygen_out: Amount,
    pub net_oxygen_transfer: Amount,
    pub transfer_conditions: QueryPoint,
    pub maximum_gas_coupling_error: Amount,
    pub accumulated_gas_coupling_error: Amount,
    transfer: Arc<OxygenTransfer>,
}
impl AeratedFlowCilAdvance {
    pub fn transfer(&self) -> &OxygenTransfer {
        &self.transfer
    }
}

impl AeratedCarbonInLeachCstr {
    /// Attach gas exchange to an admitted continuous CIL tank, including one
    /// already advanced without gas. Existing liquid/carbon history is retained;
    /// the gas boundary starts with zero cumulative transfer.
    pub fn from_tank(cil: CarbonInLeachCstr, transfer: OxygenTransfer) -> CilFlowResult<Self> {
        let layout = &cil.tank.layout;
        let model = &layout.model;
        let external_gas = transfer.solubility().gas();
        require_local(
            model.chemistry.families().iter().all(|family| {
                family.species().iter().all(|species| species != external_gas)
            }) && layout
                .charge
                .spectators
                .iter()
                .all(|ion| &ion.species != external_gas)
                && layout.competition.as_ref().is_none_or(|metals| {
                    metals.model.families().iter().all(|family| {
                        family.species().iter().all(|species| species != external_gas)
                    })
                })
                && cil.pair.counterion() != external_gas
                && cil.pair.adsorbed_pair() != external_gas,
            "external oxygen cannot also belong to retained liquid or sorbed chemistry",
        )?;
        require_local(
            [OXYGEN_PARTIAL_PRESSURE_AXIS, SATURATION_CONCENTRATION_AXIS]
                .iter()
                .all(|axis| {
                    !layout.conditions.axes().contains_key(*axis)
                        && model.source.validity.bound(axis).is_none()
                }),
            "gas pressure and saturation belong to the transfer source",
        )?;
        let gas = transfer.solubility().standard();
        let water = model.chemistry.water().standard();
        require_local(
            gas.temperature() == water.temperature()
                && gas.reference_pressure() == water.reference_pressure()
                && gas.reference_axis() == water.reference_axis(),
            "gas and liquid sources need identical reference conditions",
        )?;
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
                .chain([
                    model.reaction.clone(),
                    cil.pair.reaction().clone(),
                    transfer.solubility().reaction().clone(),
                ])
                .collect(),
        )
        .map_err(CyanidationError::from)?;
        Ok(Self {
            cil,
            transfer: Arc::new(transfer),
            oxygen_in: 0.0,
            oxygen_out: 0.0,
            coupling_error: 0.0,
        })
    }

    pub fn tank(&self) -> &CarbonInLeachCstr {
        &self.cil
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

    /// Advance at constant liquid/slurry feed, flow and prescribed O2 pressure.
    /// All may change on the next call. Refusal publishes no partial gas,
    /// liquid, carbon or clock state.
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
    ) -> CilFlowResult<AeratedFlowCilAdvance> {
        nonnegative(duration.value(), "aerated continuous CIL duration must be nonnegative")?;
        nonnegative(flow.value(), "aerated continuous CIL flow must be nonnegative")?;
        if cancelled() {
            return Err(FlowCilError::Cancelled);
        }
        let start = self.cil.tank.elapsed().value();
        let end = add(start, duration.value())?;
        let saturation = self
            .transfer
            .solubility()
            .equilibrium_aqueous_amount(
                oxygen_partial_pressure,
                self.cil.tank.layout.charge.amount_at_unit_activity,
            )
            .map_err(CyanidationError::from)?;
        self.transfer
            .validate_extras(&self.cil.tank.layout.conditions, transfer_conditions)
            .map_err(FlowCilError::from)?;
        check_transfer_clock(self.transfer.source(), end)?;

        let boundary = Boundary {
            transfer: &self.transfer,
            pressure: oxygen_partial_pressure,
            saturation,
            extra: transfer_conditions,
        };
        boundary.admit(&self.cil.state)?;

        // Validate feed/counterion/source declarations atomically without
        // advancing the accepted owner.
        let mut preflight = self.cil.clone();
        preflight.advance(
            Time::new(0.0),
            flow,
            feed,
            numerics,
            || cancelled(),
        )?;

        let gas_scale = product(
            product(self.transfer.coefficient().value(), duration.value())?,
            self.cil.tank.raw[2].max(saturation.value()),
        )?;
        let scale = sum(
            self.cil
                .tank
                .raw
                .iter()
                .copied()
                .chain(self.cil.capacities.iter().copied())
                .chain(std::iter::once(gas_scale)),
        )?;
        let tolerance = CyanidationNumerics {
            maximum_step: numerics.maximum_step,
            absolute_gold_error: numerics.absolute_amount_error,
            relative_gold_error: numerics.relative_amount_error,
            maximum_attempts: numerics.maximum_attempts,
            maximum_work: numerics.maximum_work,
        }
        .tolerance(scale)?;

        let before = self.cil.clone();
        let mut current = AeratedTrial {
            cil: self.cil.clone(),
            inward: 0.0,
            outward: 0.0,
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
                    "aerated continuous CIL clock loses resolution",
                ));
            }
            let trial = (|| -> CilFlowResult<_> {
                let full = split(
                    current.clone(),
                    &boundary,
                    flow,
                    feed,
                    numerics,
                    clock,
                    finish,
                    &mut cancelled,
                )?;
                let first = split(
                    current.clone(),
                    &boundary,
                    flow,
                    feed,
                    numerics,
                    clock,
                    middle,
                    &mut cancelled,
                )?;
                let refined = split(
                    first,
                    &boundary,
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
                Err(FlowCilError::Cyanidation(CyanidationError::TrialDepletion))
                | Err(FlowCilError::Cyanidation(CyanidationError::OutsideValidity { .. })) => {
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

        let gas_accumulated = sum([self.coupling_error, call_error].into_iter())?;
        let interval_in = differences(&current.cil.tank.inlet, &before.tank.inlet)?;
        let interval_out = differences(&current.cil.tank.outlet, &before.tank.outlet)?;
        let interval_gold = difference(
            current.cil.tank.gold_dissolved,
            before.tank.gold_dissolved,
            "cumulative gold progress regressed",
        )?;
        let interval_gas_in = current.inward;
        let interval_gas_out = current.outward;
        let oxygen_in = add(self.oxygen_in, interval_gas_in)?;
        let oxygen_out = add(self.oxygen_out, interval_gas_out)?;

        let state = admit_state(
            &current.cil.tank,
            &current.cil.pair,
            &current.cil.populations,
            &current.cil.capacities,
            current.cil.counterion_index,
            &current.cil.adsorbed,
            &current.cil.vacant,
            numerics.maximum_work,
            &mut cancelled,
        )?;
        let interval_conservation = open_conservation(
            &current.cil,
            &state,
            &before.tank.raw,
            &before.adsorbed,
            &interval_in,
            &interval_out,
            interval_gold,
            interval_gas_in,
            interval_gas_out,
        )?;
        let cumulative_conservation = open_conservation(
            &current.cil,
            &state,
            &current.cil.tank.layout.initial,
            &current.cil.initial_adsorbed,
            &current.cil.tank.inlet,
            &current.cil.tank.outlet,
            current.cil.tank.gold_dissolved,
            oxygen_in,
            oxygen_out,
        )?;
        let capture = |role| -> CilFlowResult<Amount> {
            Ok(Amount::new(sum(
                state
                    .populations
                    .iter()
                    .zip(&before.state.populations)
                    .filter(|(after, _)| after.role == role)
                    .map(|(after, prior)| after.adsorbed.value() - prior.adsorbed.value()),
            )?))
        };
        let actual_transfer_conditions = boundary.admit(&state)?;
        if cancelled() {
            return Err(FlowCilError::Cancelled);
        }

        let cil = FlowCilAdvance {
            pair: current.cil.pair.clone(),
            populations: current.cil.populations.clone(),
            previous_elapsed: Time::new(start),
            elapsed: Time::new(end),
            flow,
            feed: feed.clone(),
            state: state.clone(),
            inlet: inventory(&current.cil.tank.layout, &interval_in),
            outlet: inventory(&current.cil.tank.layout, &interval_out),
            cumulative_inlet: inventory(&current.cil.tank.layout, &current.cil.tank.inlet),
            cumulative_outlet: inventory(&current.cil.tank.layout, &current.cil.tank.outlet),
            interval_gold_dissolved: Amount::new(interval_gold),
            cumulative_gold_dissolved: Amount::new(current.cil.tank.gold_dissolved),
            interval_product_capture: capture(CarbonRole::ProductCarbon)?,
            interval_preg_robbing_capture: capture(CarbonRole::PregRobbing)?,
            interval_conservation,
            cumulative_conservation,
            maximum_coupling_error: Amount::new(largest),
            accumulated_coupling_error: Amount::new(current.cil.coupling_error),
            liquid_error_estimate: Amount::new(current.cil.tank.accumulated_error),
            accepted_steps: accepted,
            rejected_steps: rejected,
        };
        let answer = AeratedFlowCilAdvance {
            cil,
            oxygen_partial_pressure,
            oxygen_saturation: saturation,
            interval_oxygen_in: Amount::new(interval_gas_in),
            interval_oxygen_out: Amount::new(interval_gas_out),
            cumulative_oxygen_in: Amount::new(oxygen_in),
            cumulative_oxygen_out: Amount::new(oxygen_out),
            net_oxygen_transfer: Amount::new(sum([oxygen_in, -oxygen_out].into_iter())?),
            transfer_conditions: actual_transfer_conditions,
            maximum_gas_coupling_error: Amount::new(largest),
            accumulated_gas_coupling_error: Amount::new(gas_accumulated),
            transfer: self.transfer.clone(),
        };

        current.cil.state = state;
        self.cil = current.cil;
        self.oxygen_in = oxygen_in;
        self.oxygen_out = oxygen_out;
        self.coupling_error = gas_accumulated;
        Ok(answer)
    }
}

#[derive(Clone)]
struct AeratedTrial {
    cil: CarbonInLeachCstr,
    inward: f64,
    outward: f64,
}
impl AeratedTrial {
    fn values(&self) -> impl Iterator<Item = f64> + '_ {
        self.cil
            .tank
            .raw
            .iter()
            .chain(&self.cil.tank.inlet)
            .chain(&self.cil.tank.outlet)
            .copied()
            .chain(std::iter::once(self.cil.tank.gold_dissolved))
            .chain(self.cil.adsorbed.iter().copied())
            .chain(self.cil.vacant.iter().copied())
            .chain([self.inward, self.outward])
    }
    fn error(&self, other: &Self) -> CilFlowResult<f64> {
        Ok(sum(
            self.values()
                .zip(other.values())
                .map(|(a, b)| (a - b).abs()),
        )?)
    }
}

struct Boundary<'a> {
    transfer: &'a OxygenTransfer,
    pressure: Pressure,
    saturation: Amount,
    extra: &'a QueryPoint,
}
impl Boundary<'_> {
    fn admit(&self, state: &FlowCilState) -> CilFlowResult<QueryPoint> {
        let point = state
            .solution
            .conditions
            .clone()
            .with_quantity(
                ELAPSED_TIME_AXIS,
                QuantitySpec::dimensional(Time::DIMS),
                state.elapsed.value(),
            )
            .map_err(CyanidationError::from)?;
        self.transfer
            .rate_at_conditions(
                &point,
                state.solution.inventory.dissolved_oxygen,
                self.pressure,
                self.saturation,
                self.extra,
            )
            .map(|(_, point)| point)
            .map_err(FlowCilError::from)
    }
}

#[allow(clippy::too_many_arguments)]
fn split(
    mut trial: AeratedTrial,
    boundary: &Boundary<'_>,
    flow: VolumetricFlowRate,
    feed: &CyanidationFeed,
    numerics: FlowNumerics,
    start: f64,
    end: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> CilFlowResult<AeratedTrial> {
    let half = product(end - start, 0.5)?;
    gas_step(&mut trial, boundary, half, numerics.maximum_work, cancelled)?;
    trial.cil.advance(
        Time::new(end - start),
        flow,
        feed,
        FlowNumerics {
            maximum_step: Time::new(numerics.maximum_step.value().min(end - start)),
            ..numerics
        },
        || cancelled(),
    )?;
    gas_step(&mut trial, boundary, half, numerics.maximum_work, cancelled)?;
    Ok(trial)
}

fn gas_step(
    trial: &mut AeratedTrial,
    boundary: &Boundary<'_>,
    duration: f64,
    maximum_work: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> CilFlowResult<()> {
    if cancelled() {
        return Err(FlowCilError::Cancelled);
    }
    if duration == 0.0 {
        boundary.admit(&trial.cil.state)?;
        return Ok(());
    }
    let before = trial.cil.tank.raw[2];
    let saturation = boundary.saturation.value();
    let exposure = product(boundary.transfer.coefficient().value(), duration)?;
    if exposure > 0.0 && before != saturation {
        let e = det::exp(-exposure);
        let u = -det::expm1(-exposure);
        let magnitude = product((saturation - before).abs(), u)?;
        let after = sum([product(before, e)?, product(saturation, u)?].into_iter())?;
        if after == before || (before > 0.0 && saturation == 0.0 && e == 0.0) {
            return Err(FlowCilError::Arithmetic(
                "finite gas exchange loses oxygen inventory resolution",
            ));
        }
        trial.cil.tank.raw[2] = after;
        if saturation > before {
            trial.inward = add(trial.inward, magnitude)?;
        } else {
            trial.outward = add(trial.outward, magnitude)?;
        }
    }
    let state = admit_state(
        &trial.cil.tank,
        &trial.cil.pair,
        &trial.cil.populations,
        &trial.cil.capacities,
        trial.cil.counterion_index,
        &trial.cil.adsorbed,
        &trial.cil.vacant,
        maximum_work,
        cancelled,
    )?;
    boundary.admit(&state)?;
    trial.cil.tank.state = state.solution.clone();
    trial.cil.state = state;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn open_conservation(
    cil: &CarbonInLeachCstr,
    state: &FlowCilState,
    initial_raw: &[f64],
    initial_adsorbed: &[f64],
    incoming: &[f64],
    outgoing: &[f64],
    gold: f64,
    gas_in: f64,
    gas_out: f64,
) -> CilFlowResult<FlowCilConservation> {
    let mut all_in = incoming.to_vec();
    let mut all_out = outgoing.to_vec();
    all_in[2] = sum([all_in[2], gas_in].into_iter())?;
    all_out[2] = sum([all_out[2], gas_out].into_iter())?;
    conservation_cil_flow(
        &cil.tank.layout,
        &cil.pair,
        cil.counterion_index,
        &cil.capacities,
        state,
        initial_raw,
        initial_adsorbed,
        &all_in,
        &all_out,
        gold,
    )
}

fn check_transfer_clock(source: &ConstitutiveModelCard, end: f64) -> CilFlowResult<()> {
    let (minimum, maximum) = source
        .validity
        .bound(ELAPSED_TIME_AXIS)
        .ok_or(FlowCilError::Invalid("transfer source is missing absolute clock"))?;
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
