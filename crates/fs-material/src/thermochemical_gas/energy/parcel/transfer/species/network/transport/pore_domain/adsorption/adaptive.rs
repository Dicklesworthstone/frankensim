//! Temporal accuracy over the original unsplit adsorbing-pore candidates.
//!
//! Each probe is exactly one existing native interval. A full interval estimates
//! the first-order local error of two half intervals; only the latter can enter
//! accepted state and accounts. No state extrapolation or replacement physics is
//! introduced. This local asymptotic estimate is not a global error certificate.

use super::*;

/// Componentwise absolute scales and a relative allowance on state CHANGES.
/// Large carrier inventories or absolute temperatures cannot hide small changes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdsorbingPoreTimeControls {
    pub gas_species_tolerance_kg: f64,
    pub bound_species_tolerance_mol: f64,
    pub gas_temperature_tolerance_k: f64,
    pub solid_temperature_tolerance_k: f64,
    /// Multiplies the largest coarse/fine change from the last accepted state.
    /// Must lie in [0,1); zero uses only absolute tolerances.
    pub relative_tolerance: f64,
}

/// Accepted whole-duration physical accounts plus bounded accuracy/work results.
#[derive(Clone, Debug, PartialEq)]
pub struct AdsorbingPoreAdvanceReport {
    pub interval_s: [f64; 2],
    pub accepted_intervals: usize,
    pub rejected_intervals: usize,
    /// Counts all native full/half probes, including rejected and failed probes.
    pub evaluations: usize,
    /// Gas species, bound species, gas temperature, solid temperature order.
    pub maximum_accepted_error_ratios: [f64; 4],
    pub suggested_step_s: f64,
    /// Only accepted half steps contribute transfers and energy/entropy.
    /// Original-to-final species/mass admission uses one fixed physical budget.
    /// `pore.attempts` equals evaluations; the final native source report remains.
    pub physical: AdsorbingPoreStep,
}

#[derive(Clone, Debug, PartialEq)]
pub struct AdsorbingPoreAdvanceCandidate {
    domain: AdsorbingPoreDomain,
    report: AdsorbingPoreAdvanceReport,
}
impl AdsorbingPoreAdvanceCandidate {
    pub fn domain(&self) -> &AdsorbingPoreDomain {
        &self.domain
    }
    pub fn report(&self) -> &AdsorbingPoreAdvanceReport {
        &self.report
    }
    pub fn into_parts(self) -> (AdsorbingPoreDomain, AdsorbingPoreAdvanceReport) {
        (self.domain, self.report)
    }
}

/// Every refusal preserves the complete original gas, site, thermal and clock state.
#[derive(Clone, Debug, PartialEq)]
pub enum AdsorbingPoreAdvanceError {
    Invalid(&'static str),
    Physical { time_s: f64, step_s: f64, error: PoreDomainError },
    Accuracy { time_s: f64, step_s: f64, error_ratios: [f64; 4] },
    TimeResolution { time_s: f64, target_s: f64 },
    WorkLimit { time_s: f64, accepted_substeps: usize, evaluations: usize },
    Cancelled,
}
impl core::fmt::Display for AdsorbingPoreAdvanceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive adsorbing pore: {self:?}")
    }
}
impl std::error::Error for AdsorbingPoreAdvanceError {}
use AdsorbingPoreAdvanceError as Error;

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() { Err(Error::Cancelled) } else { Ok(()) }
}
fn physical(error: PoreDomainError, time_s: f64, step_s: f64) -> Error {
    if error == PoreDomainError::Cancelled { Error::Cancelled } else { Error::Physical { time_s, step_s, error } }
}
fn can_retry(error: &Error) -> bool {
    match error {
        Error::Physical { error: PoreDomainError::Adsorption { error, .. }, .. } => super::retryable(error),
        Error::Physical { error: PoreDomainError::Transport { error, .. }, .. } => super::super::retryable(error),
        Error::Accuracy { .. } => true,
        _ => false,
    }
}
fn add(total: &mut Sum, value: f64) -> Result<(), Error> {
    total.add(value).map_err(|_| Error::Invalid("unrepresentable adaptive pore account"))
}
fn ratio(old: f64, coarse: f64, fine: f64, absolute: f64, relative: f64) -> Result<f64, Error> {
    let allowance = absolute + relative * (coarse - old).abs().max((fine - old).abs());
    let ratio = (coarse - fine).abs() / allowance;
    if !positive(allowance) || !ratio.is_finite() {
        Err(Error::Invalid("unrepresentable componentwise pore error allowance"))
    } else {
        Ok(ratio)
    }
}

fn errors(
    old: &AdsorbingPoreDomain,
    coarse: &AdsorbingPoreDomain,
    fine: &AdsorbingPoreDomain,
    controls: AdsorbingPoreTimeControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[f64; 4], Error> {
    let mut result = [0.0_f64; 4];
    // Resolve each original mixture inventory once, not once per species. The
    // retained catalog includes absent species without seeding/dropping traces.
    for (&id, before) in &old.domain.parcels {
        poll(cancelled)?;
        let (a, b) = (&coarse.domain.parcels[&id], &fine.domain.parcels[&id]);
        let get = |gas: &ClosedGasParcel| {
            inventory(gas.model(), gas.mass_kg()).map_err(|error| physical(gas_error(error), old.domain.time_s, 0.0))
        };
        let (m0, ma, mb) = (get(before)?, get(a)?, get(b)?);
        for source in old.catalog.thermodynamics().components() {
            let name = source.species();
            result[0] = result[0].max(ratio(
                m0.get(name).copied().unwrap_or(0.0),
                ma.get(name).copied().unwrap_or(0.0),
                mb.get(name).copied().unwrap_or(0.0),
                controls.gas_species_tolerance_kg,
                controls.relative_tolerance,
            )?);
        }
        result[2] = result[2].max(ratio(
            before.state().state().temperature,
            a.state().state().temperature,
            b.state().state().temperature,
            controls.gas_temperature_tolerance_k,
            controls.relative_tolerance,
        )?);
    }
    for (&cell, before) in &old.solids {
        poll(cancelled)?;
        let (a, b) = (&coarse.solids[&cell], &fine.solids[&cell]);
        for (species, moles) in before.loadings() {
            result[1] = result[1].max(ratio(
                moles,
                a.adsorbed_moles(species).unwrap(),
                b.adsorbed_moles(species).unwrap(),
                controls.bound_species_tolerance_mol,
                controls.relative_tolerance,
            )?);
        }
        result[3] = result[3].max(ratio(
            before.temperature_k(),
            a.temperature_k(),
            b.temperature_k(),
            controls.solid_temperature_tolerance_k,
            controls.relative_tolerance,
        )?);
    }
    Ok(result)
}

// Integrate only the already admitted native reports; all thermodynamics stay
// in that owner. These maps preserve gross direction changes across the window.
#[derive(Default)]
struct Accepted {
    energy: Sum,
    entropy: Sum,
    solid_energy: Sum,
    faces: BTreeMap<(usize, SpeciesId), [Sum; 3]>,
    adsorption: BTreeMap<(usize, SpeciesId), [Sum; 2]>,
    conduction: BTreeMap<u64, ([usize; 2], [Sum; 2])>,
    gas_heat: BTreeMap<u64, (usize, [Sum; 2])>,
}
impl Accepted {
    fn record(&mut self, step: &AdsorbingPoreStep, cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
        add(&mut self.energy, step.pore.total_energy_residual_j)?;
        add(&mut self.entropy, step.pore.entropy_change_j_k)?;
        add(&mut self.solid_energy, step.solid_energy_change_j)?;
        for row in &step.pore.faces {
            poll(cancelled)?;
            let sum = self.faces.entry((row.face, row.species.clone())).or_default();
            for (total, value) in
                sum.iter_mut().zip([row.mass_left_to_right_kg, row.mass_right_to_left_kg, row.enthalpy_left_to_right_j])
            {
                add(total, value)?;
            }
        }
        for row in &step.adsorption {
            poll(cancelled)?;
            let sum = self.adsorption.entry((row.cell, row.species.clone())).or_default();
            add(&mut sum[0], row.gross_adsorption_mol)?;
            add(&mut sum[1], row.gross_desorption_mol)?;
        }
        for row in &step.conduction {
            poll(cancelled)?;
            let (_, sum) = self.conduction.entry(row.id).or_insert_with(|| (row.cells, Default::default()));
            add(&mut sum[0], row.gross_a_to_b_j)?;
            add(&mut sum[1], row.gross_b_to_a_j)?;
        }
        for row in &step.gas_heat {
            poll(cancelled)?;
            let (_, sum) = self.gas_heat.entry(row.id).or_insert_with(|| (row.cell, Default::default()));
            add(&mut sum[0], row.gross_gas_to_solid_j)?;
            add(&mut sum[1], row.gross_solid_to_gas_j)?;
        }
        Ok(())
    }

    fn finish(
        self,
        original: &AdsorbingPoreDomain,
        final_state: &AdsorbingPoreDomain,
        mut last: AdsorbingPoreStep,
        controls: PoreDomainSolve,
        entropy_limit: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AdsorbingPoreStep, Error> {
        let (mass, species, pore_mass, bound_mass) =
            super::inventories(original, &final_state.domain, &final_state.solids, cancelled)
                .map_err(|e| physical(e, final_state.domain.time_s, 0.0))?;
        if mass.abs() > controls.transport.exchange.mass_tolerance_kg
            || species > controls.transport.exchange.mass_tolerance_kg
            || self.energy.value.abs() > controls.transport.exchange.energy_tolerance_j
            || self.entropy.value < -entropy_limit
        {
            return Err(Error::Invalid(
                "complete adaptive pore duration exceeded its single original species/mass/energy/entropy budget",
            ));
        }
        last.pore.total_mass_residual_kg = mass;
        last.pore.maximum_species_mass_residual_kg = species;
        last.pore.total_energy_residual_j = self.energy.value;
        last.pore.entropy_change_j_k = self.entropy.value;
        last.pore.pore_mass_change_kg = pore_mass;
        last.adsorbed_mass_change_kg = bound_mass;
        last.solid_energy_change_j = self.solid_energy.value;
        last.pore.faces = self
            .faces
            .into_iter()
            .map(|((face, species), sum)| PoreFaceTransfer {
                face,
                species,
                mass_left_to_right_kg: sum[0].value,
                mass_right_to_left_kg: sum[1].value,
                enthalpy_left_to_right_j: sum[2].value,
            })
            .collect();
        last.adsorption = self
            .adsorption
            .into_iter()
            .map(|((cell, species), sum)| PoreAdsorptionTransfer {
                cell,
                species,
                gross_adsorption_mol: sum[0].value,
                gross_desorption_mol: sum[1].value,
            })
            .collect();
        last.conduction = self
            .conduction
            .into_iter()
            .map(|(id, (cells, sum))| PoreConductionTransfer {
                id,
                cells,
                gross_a_to_b_j: sum[0].value,
                gross_b_to_a_j: sum[1].value,
            })
            .collect();
        last.gas_heat = self
            .gas_heat
            .into_iter()
            .map(|(id, (cell, sum))| PoreGasHeatTransfer {
                id,
                cell,
                gross_gas_to_solid_j: sum[0].value,
                gross_solid_to_gas_j: sum[1].value,
            })
            .collect();
        Ok(last)
    }
}

#[allow(clippy::too_many_arguments)]
fn probe<C: FnMut() -> bool>(
    state: &AdsorbingPoreDomain,
    dt: f64,
    accepted: usize,
    evaluations: &mut usize,
    maximum: usize,
    advance: &mut impl FnMut(&AdsorbingPoreDomain, f64, &mut C) -> Result<AdsorbingPoreCandidate, PoreDomainError>,
    cancelled: &mut C,
) -> Result<AdsorbingPoreCandidate, Error> {
    poll(cancelled)?;
    if *evaluations == maximum {
        return Err(Error::WorkLimit {
            time_s: state.domain.time_s,
            accepted_substeps: 2 * accepted,
            evaluations: *evaluations,
        });
    }
    *evaluations += 1;
    advance(state, dt, cancelled).map_err(|e| physical(e, state.domain.time_s, dt))
}

impl AdsorbingPoreDomain {
    /// Propose the complete duration to an absolute target with componentwise
    /// gas species/temperature and bound species/solid-temperature accuracy.
    /// Fixed faces, source claims, geometry and local controls remain unchanged.
    ///
    /// `controls.maximum_step_s` bounds comparison spans; its minimum is a retry
    /// floor (a shorter final remainder still needs accuracy). Accepted physical
    /// halves may be smaller. `maximum_substeps` bounds accepted half steps;
    /// `maximum_attempts` bounds ALL native probes, and must be at least three.
    /// Physical failures retain their original typed cause; only existing typed
    /// timestep failures or temporal error may retry. Source/domain errors refuse.
    /// Only accepted fine paths spend inventory. The final candidate checks one
    /// original-to-final physical budget and publishes no prefix on any failure.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate_adaptive_to<C: FnMut() -> bool>(
        &self,
        target_s: f64,
        faces: &[PoreFace<'_>],
        endpoints: &[GasEnergySolve],
        controls: PoreDomainSolve,
        adsorption_controls: &BTreeMap<usize, CompetitiveSolve>,
        links: &[AdsorbentConductionLink<'_, '_>],
        gas_links: &[AdsorbentGasHeatLink<'_>],
        accuracy: AdsorbingPoreTimeControls,
        mut cancelled: C,
    ) -> Result<AdsorbingPoreAdvanceCandidate, Error> {
        poll(&mut cancelled)?;
        let start = self.domain.time_s;
        if !target_s.is_finite()
            || target_s <= start
            || !positive(target_s - start)
            || ![
                controls.minimum_step_s,
                controls.maximum_step_s,
                accuracy.gas_species_tolerance_kg,
                accuracy.bound_species_tolerance_mol,
                accuracy.gas_temperature_tolerance_k,
                accuracy.solid_temperature_tolerance_k,
            ]
            .iter()
            .all(|x| positive(*x))
            || controls.minimum_step_s > controls.maximum_step_s
            || !(2..=MAX_PORE_SUBSTEPS).contains(&controls.maximum_substeps)
            || !(3..=MAX_PORE_ATTEMPTS).contains(&controls.maximum_attempts)
            || !accuracy.relative_tolerance.is_finite()
            || !(0.0..1.0).contains(&accuracy.relative_tolerance)
        {
            return Err(Error::Invalid(
                "forward target, positive temporal scales and bounded comparison/evaluation controls required",
            ));
        }
        // Equal floor/ceiling prevent native retries from hiding a physical
        // refusal behind an inner work limit. Geometry/source admission and all
        // original local/global physical checks still run for this one interval.
        let mut advance = |state: &Self, dt, cancelled: &mut C| {
            state.step_candidate_with_heat(
                dt,
                faces,
                endpoints,
                PoreDomainSolve {
                    maximum_step_s: dt,
                    minimum_step_s: dt,
                    maximum_substeps: 1,
                    maximum_attempts: 1,
                    ..controls
                },
                adsorption_controls,
                links,
                gas_links,
                cancelled,
            )
        };
        let mut state = self.clone();
        let mut ceiling = controls.maximum_step_s;
        let (mut accepted, mut rejected, mut evaluations) = (0, 0, 0);
        let (mut minimum, mut maximum) = (f64::INFINITY, 0.0_f64);
        let mut maxima = [0.0_f64; 4];
        let mut ledger = Accepted::default();
        let mut last = None;
        while state.domain.time_s < target_s {
            poll(&mut cancelled)?;
            let time = state.domain.time_s;
            if 2 * accepted + 2 > controls.maximum_substeps {
                return Err(Error::WorkLimit { time_s: time, accepted_substeps: 2 * accepted, evaluations });
            }
            let end = if ceiling >= target_s - time { target_s } else { time + ceiling };
            let span = end - time;
            let middle = time + 0.5 * span;
            if !(time < middle && middle < end) {
                return Err(Error::TimeResolution { time_s: time, target_s: end });
            }
            let (h1, h2) = (middle - time, end - middle);
            let attempted = (|| {
                let coarse = probe(
                    &state,
                    span,
                    accepted,
                    &mut evaluations,
                    controls.maximum_attempts,
                    &mut advance,
                    &mut cancelled,
                )?;
                let first = probe(
                    &state,
                    h1,
                    accepted,
                    &mut evaluations,
                    controls.maximum_attempts,
                    &mut advance,
                    &mut cancelled,
                )?;
                if first.domain.domain.time_s.to_bits() != middle.to_bits() {
                    return Err(Error::TimeResolution { time_s: time, target_s: middle });
                }
                let second = probe(
                    first.domain(),
                    h2,
                    accepted,
                    &mut evaluations,
                    controls.maximum_attempts,
                    &mut advance,
                    &mut cancelled,
                )?;
                if coarse.domain.domain.time_s.to_bits() != end.to_bits()
                    || second.domain.domain.time_s.to_bits() != end.to_bits()
                {
                    return Err(Error::TimeResolution { time_s: time, target_s: end });
                }
                let mut ratios = errors(&state, coarse.domain(), second.domain(), accuracy, &mut cancelled)?;
                // Fine first-order local-error factor for unequal represented halves.
                let (u, v) = (h1 / span, h2 / span);
                let factor = (u * u + v * v) / (2.0 * u * v);
                for ratio in &mut ratios {
                    *ratio *= factor;
                }
                if ratios.iter().any(|r| *r > 1.0) {
                    return Err(Error::Accuracy { time_s: time, step_s: span, error_ratios: ratios });
                }
                Ok((first, second, ratios))
            })();
            match attempted {
                Ok((first, second, ratios)) => {
                    ledger.record(first.step(), &mut cancelled)?;
                    ledger.record(second.step(), &mut cancelled)?;
                    let (next, step) = second.into_parts();
                    state = next;
                    last = Some(step);
                    accepted += 1;
                    minimum = minimum.min(h1.min(h2));
                    maximum = maximum.max(h1.max(h2));
                    for (max, r) in maxima.iter_mut().zip(ratios) {
                        *max = max.max(r);
                    }
                    ceiling = if ratios.iter().all(|r| *r <= 0.25) {
                        span + span.min(controls.maximum_step_s - span)
                    } else {
                        span
                    };
                    ceiling = ceiling.clamp(controls.minimum_step_s, controls.maximum_step_s);
                }
                Err(error) => {
                    let smaller = (0.5 * span).max(controls.minimum_step_s);
                    if can_retry(&error) && smaller < span {
                        rejected += 1;
                        ceiling = smaller;
                    } else {
                        return Err(error);
                    }
                }
            }
        }
        let entropy_limit = adsorption_controls
            .values()
            .fold(controls.transport.entropy_tolerance_j_k, |limit, solve| limit.min(solve.entropy_tolerance_j_k));
        let mut step = ledger.finish(
            self,
            &state,
            last.expect("positive duration accepted fine steps"),
            controls,
            entropy_limit,
            &mut cancelled,
        )?;
        step.pore.duration_s = target_s - start;
        step.pore.accepted_substeps = 2 * accepted;
        step.pore.attempts = evaluations;
        step.pore.minimum_accepted_step_s = minimum;
        step.pore.maximum_accepted_step_s = maximum;
        poll(&mut cancelled)?;
        Ok(AdsorbingPoreAdvanceCandidate {
            domain: state,
            report: AdsorbingPoreAdvanceReport {
                interval_s: [start, target_s],
                accepted_intervals: accepted,
                rejected_intervals: rejected,
                evaluations,
                maximum_accepted_error_ratios: maxima,
                suggested_step_s: ceiling,
                physical: step,
            },
        })
    }

    /// Publish every gas/solid inventory, temperature and clock together only
    /// after the entire accuracy-controlled duration and final cancellation gate.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_adaptive_to(
        &mut self,
        target_s: f64,
        faces: &[PoreFace<'_>],
        endpoints: &[GasEnergySolve],
        controls: PoreDomainSolve,
        adsorption_controls: &BTreeMap<usize, CompetitiveSolve>,
        links: &[AdsorbentConductionLink<'_, '_>],
        gas_links: &[AdsorbentGasHeatLink<'_>],
        accuracy: AdsorbingPoreTimeControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AdsorbingPoreAdvanceReport, Error> {
        let candidate = self.step_candidate_adaptive_to(
            target_s,
            faces,
            endpoints,
            controls,
            adsorption_controls,
            links,
            gas_links,
            accuracy,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        let (next, report) = candidate.into_parts();
        *self = next;
        Ok(report)
    }
}
