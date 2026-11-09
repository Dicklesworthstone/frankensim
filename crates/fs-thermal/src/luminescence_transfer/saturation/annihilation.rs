//! Triplet-pair quenching and delayed fluorescence with conserved emitters.
//!
//! The original pair law removes two triplets per event, creates an authored
//! mean number of singlets, and returns the remaining emitters to ground.
//! Each positive exponential-midpoint stage freezes lambda = 2*k*N_T/V inside
//! the SAME ground/singlet/triplet propagator as pumping and all original decay.
//! Its temporary T->S and T->G coordinates do not become physical crossing or
//! dark-decay channels: original events and pair heat are reconstructed directly.
//! No pair kick, clipped population or residual-to-heat correction is applied.

use super::super::pump_column::{Sum, cap, kinetic_energy_residual, poll, product3};
use super::super::{
    SingletTripletDecay, SingletTripletError, finite, nonnegative, positive_product,
};
use super::{SaturableSingletTripletPlan, SaturableSingletTripletStep};

pub use super::super::annihilation::TripletAnnihilation;
pub use super::super::recycling::annihilation::AnnihilationLimits;

#[path = "annihilation/integrate.rs"]
pub(crate) mod integrate;
pub(crate) use integrate::{PairHistory, midpoint};

/// One physical candidate plus the controlled nonlinear integration work.
#[derive(Clone, Debug, PartialEq)]
pub struct SaturableAnnihilationStep {
    pub kinetics: SaturableSingletTripletStep,
    /// Zero for the unchanged linear path, one for exact isolated pair kinetics.
    pub substeps: usize,
    /// Includes every discarded grid, predictor and corrector evaluation.
    pub kinetic_evaluations: usize,
    /// Componentwise energy discrepancy, not a rigorous error enclosure.
    pub max_energy_discrepancy_j: f64,
}

impl SaturableSingletTripletPlan {
    /// Couple an optional authored pair law to the original fixed pump hazards.
    /// Limits and event energetics are admitted even at zero duration. Original
    /// lifetimes use the native exponential; a graded true-time endpoint tail
    /// resolves fast delayed fluorescence without a lifetime-sized global grid.
    pub fn advance_with_annihilation(
        &self,
        populations: [f64; 3],
        law: Option<TripletAnnihilation>,
        limits: AnnihilationLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SaturableAnnihilationStep, SingletTripletError> {
        let active = validate(
            &[self.model],
            &[law],
            &[populations],
            self.duration_s,
            limits,
            &mut cancelled,
        )?;
        if !active || self.duration_s == 0.0 {
            return Ok(SaturableAnnihilationStep {
                kinetics: self.advance_with_cancel(populations, &mut cancelled)?,
                substeps: 0,
                kinetic_evaluations: 1,
                max_energy_discrepancy_j: 0.0,
            });
        }
        if self.pump_rates_s == [0.0; 2] && is_stable(self.model) {
            let mut kinetics = self.advance_with_cancel(populations, &mut cancelled)?;
            apply_isolated(
                self.model,
                populations,
                law,
                self.duration_s,
                &mut kinetics,
                self.relative_tolerance,
            )?;
            return Ok(SaturableAnnihilationStep {
                kinetics,
                substeps: 1,
                kinetic_evaluations: 1,
                max_energy_discrepancy_j: 0.0,
            });
        }
        let mut evaluations = 0usize;
        let fastest = fastest_rate(&[self.model]).max(self.pump_rates_s.iter().sum());
        let result = midpoint(
            &[populations],
            &[law],
            self.duration_s,
            fastest,
            limits,
            None,
            |state, rates, duration, cancelled| {
                evaluations =
                    evaluations
                        .checked_add(1)
                        .ok_or(SingletTripletError::Unrepresentable(
                            "finite pair evaluation counter",
                        ))?;
                let frozen = frozen_model(self.model, law, rates[0])?;
                let plan = SaturableSingletTripletPlan::prepare(
                    frozen,
                    self.pump_rates_s,
                    self.pump_energies_j,
                    duration,
                    self.limits,
                    &mut *cancelled,
                )?;
                let mut step = plan.advance_with_cancel(state[0], &mut *cancelled)?;
                restore(
                    self.model,
                    frozen,
                    state[0],
                    law,
                    rates[0],
                    &mut step,
                    self.relative_tolerance,
                )?;
                Ok(Some(step))
            },
            |step| audit(self.model, populations, step, self.relative_tolerance),
            |step| energy_components(self.model, step, self.duration_s),
            &mut cancelled,
        )?;
        Ok(SaturableAnnihilationStep {
            kinetics: result.step,
            substeps: result.substeps,
            kinetic_evaluations: evaluations,
            max_energy_discrepancy_j: result.discrepancy_j,
        })
    }
}

pub(crate) fn validate(
    models: &[SingletTripletDecay],
    laws: &[Option<TripletAnnihilation>],
    initial: &[[f64; 3]],
    duration: f64,
    limits: AnnihilationLimits,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<bool, SingletTripletError> {
    poll(cancelled)?;
    if laws.len() != models.len()
        || initial.len() != models.len()
        || !initial.iter().flatten().copied().all(nonnegative)
        || !nonnegative(duration)
        || !nonnegative(limits.absolute_energy_tolerance_j)
        || !nonnegative(limits.relative_tolerance)
        || limits.relative_tolerance >= 1.0
        || (limits.relative_tolerance == 0.0 && limits.absolute_energy_tolerance_j == 0.0)
        || limits.maximum_substeps == 0
    {
        return Err(SingletTripletError::Invalid(
            "finite annihilation populations, laws or limits",
        ));
    }
    let mut active = false;
    for (model, law) in models.iter().zip(laws) {
        poll(cancelled)?;
        if let Some(law) = law {
            law.advance(model.excitation_energies_j(), [0.0; 2], 0.0)?;
            active |= law.pair_rate_m3_s() > 0.0;
        }
    }
    Ok(active)
}

pub(crate) fn is_stable(model: SingletTripletDecay) -> bool {
    model.total_rates_s == [0.0; 2]
}

pub(crate) fn fastest_rate(models: &[SingletTripletDecay]) -> f64 {
    models
        .iter()
        .flat_map(|m| m.total_rates_s)
        .fold(0.0_f64, f64::max)
}

pub(crate) fn frozen_rates(
    state: &[[f64; 3]],
    laws: &[Option<TripletAnnihilation>],
) -> Result<Vec<f64>, SingletTripletError> {
    state
        .iter()
        .zip(laws)
        .map(|(p, law)| law.map_or(Ok(0.0), |law| law.frozen_removal_rate_s(p[2])))
        .collect()
}

pub(crate) fn frozen_model(
    model: SingletTripletDecay,
    law: Option<TripletAnnihilation>,
    lambda: f64,
) -> Result<SingletTripletDecay, SingletTripletError> {
    if lambda == 0.0 {
        return Ok(model);
    }
    let law = law.ok_or(SingletTripletError::Invalid(
        "missing active finite pair law",
    ))?;
    let singlet_rate = product3(lambda, law.singlet_yield(), 0.5)?;
    let ground_rate = positive_product(
        lambda,
        1.0 - 0.5 * law.singlet_yield(),
        "pair ground return rate",
    )?;
    let mut dark = model.nonradiative_rates_s();
    let mut crossing = model.crossing_rates_s();
    dark[1] = finite(dark[1] + ground_rate, "combined finite triplet removal")?;
    crossing[1] = finite(
        crossing[1] + singlet_rate,
        "combined finite triplet conversion",
    )?;
    SingletTripletDecay::try_new(
        model.excitation_energies_j(),
        model.radiative_rates_s(),
        dark,
        crossing,
    )
}

/// Restore physical channels from the dominant native removal count; its
/// fraction is at least 1/3. Never subtract a large artificial bath heat.
pub(crate) fn restore(
    original: SingletTripletDecay,
    frozen: SingletTripletDecay,
    initial: [f64; 3],
    law: Option<TripletAnnihilation>,
    lambda: f64,
    step: &mut SaturableSingletTripletStep,
    tolerance: f64,
) -> Result<(), SingletTripletError> {
    if lambda > 0.0 {
        let law = law.ok_or(SingletTripletError::Invalid("missing finite pair law"))?;
        let fractions = frozen.fractions[1];
        let events = [
            step.radiative_events[1],
            step.nonradiative_events[1],
            step.crossing_events[1],
        ];
        let mut dominant = 0;
        for i in 1..3 {
            if fractions[i] > fractions[dominant] {
                dominant = i;
            }
        }
        let removals = finite(
            events[dominant] / fractions[dominant],
            "finite triplet removals",
        )?;
        let total_rate = frozen.total_rates_s[1];
        step.annihilation_events = positive_product(
            channel(removals, lambda, total_rate)?,
            0.5,
            "triplet pair events",
        )?;
        (step.annihilation_singlets, step.annihilation_heat_j) =
            law.event_channels(original.excitation_energies_j(), step.annihilation_events)?;
        step.nonradiative_events[1] =
            channel(removals, original.nonradiative_rates_s()[1], total_rate)?;
        step.crossing_events[1] = channel(removals, original.crossing_rates_s()[1], total_rate)?;
        let energies = original.excitation_energies_j();
        step.nonradiative_heat_j[1] = positive_product(
            step.nonradiative_events[1],
            energies[1],
            "original triplet heat",
        )?;
        step.reverse_crossing_heat_absorbed_j = positive_product(
            step.crossing_events[1],
            energies[0] - energies[1],
            "original reverse crossing heat",
        )?;
    }
    audit(original, initial, step, tolerance)
}

fn channel(removals: f64, rate: f64, total: f64) -> Result<f64, SingletTripletError> {
    let fraction = rate / total;
    if rate > 0.0 && fraction == 0.0 {
        return Err(SingletTripletError::Unrepresentable(
            "positive finite pair channel",
        ));
    }
    positive_product(removals, fraction, "physical finite triplet events")
}

pub(crate) fn apply_isolated(
    model: SingletTripletDecay,
    initial: [f64; 3],
    law: Option<TripletAnnihilation>,
    duration: f64,
    step: &mut SaturableSingletTripletStep,
    tolerance: f64,
) -> Result<(), SingletTripletError> {
    if let Some(law) = law {
        let pair = law.advance(
            model.excitation_energies_j(),
            [initial[1], initial[2]],
            duration,
        )?;
        step.annihilation_events = pair.pair_events;
        step.annihilation_singlets = pair.singlet_production;
        step.annihilation_heat_j = pair.nonradiative_heat_j;
        let ground_return = positive_product(
            pair.pair_events,
            2.0 - law.singlet_yield(),
            "pair emitter ground return",
        )?;
        step.populations = [
            finite(initial[0] + ground_return, "pair ground population")?,
            pair.populations[0],
            pair.populations[1],
        ];
        for target in 0..2 {
            step.stored_energy_j[target] = positive_product(
                step.populations[target + 1],
                model.excitation_energies_j()[target],
                "pair excitation store",
            )?;
        }
    }
    audit(model, initial, step, tolerance)
}

pub(crate) fn admit_excited(
    initial: [f64; 3],
    step: &SaturableSingletTripletStep,
    tolerance: f64,
) -> Result<(), SingletTripletError> {
    if step.annihilation_events == 0.0 {
        return super::super::admit_population_balance(
            &[initial[1], initial[2]],
            &step.pump_events,
            &[step.populations[1], step.populations[2]],
            &step.radiative_events,
            &step.nonradiative_events,
            tolerance,
        );
    }
    let mut residual = Sum::default();
    let mut scale = Sum::default();
    for i in 0..2 {
        for term in [
            step.populations[i + 1],
            -initial[i + 1],
            -step.pump_events[i],
            step.radiative_events[i],
            step.nonradiative_events[i],
        ] {
            residual.add(term)?;
            scale.add(term.abs())?;
        }
    }
    let lost = finite(
        step.annihilation_events + (step.annihilation_events - step.annihilation_singlets),
        "pair excitation loss",
    )?;
    residual.add(lost)?;
    scale.add(lost)?;
    let residual = residual.value()?;
    let allowance = tolerance * scale.value()?;
    if residual.abs() > allowance {
        return Err(SingletTripletError::PopulationBalance {
            residual,
            allowance,
        });
    }
    Ok(())
}

pub(crate) fn audit(
    model: SingletTripletDecay,
    initial: [f64; 3],
    step: &mut SaturableSingletTripletStep,
    tolerance: f64,
) -> Result<(), SingletTripletError> {
    admit_excited(initial, step, tolerance)?;
    let mut residual = Sum::default();
    let mut scale = Sum::default();
    for i in 0..3 {
        residual.add(step.populations[i])?;
        residual.add(-initial[i])?;
        scale.add(step.populations[i])?;
        scale.add(initial[i])?;
    }
    step.population_residual = residual.value()?;
    let allowance = tolerance * scale.value()?;
    if step.population_residual.abs() > allowance {
        return Err(SingletTripletError::PopulationBalance {
            residual: step.population_residual,
            allowance,
        });
    }
    step.residual_j = kinetic_energy_residual(model, initial, step, tolerance)?;
    Ok(())
}

pub(crate) fn visit_kinetics(
    step: &mut SaturableSingletTripletStep,
    visitor: &mut impl FnMut(&mut f64) -> Result<(), SingletTripletError>,
) -> Result<(), SingletTripletError> {
    for values in [
        &mut step.pump_events,
        &mut step.absorbed_pump_energy_j,
        &mut step.pump_relaxation_heat_j,
        &mut step.radiative_events,
        &mut step.nonradiative_events,
        &mut step.crossing_events,
        &mut step.radiative_decay_energy_j,
        &mut step.nonradiative_heat_j,
    ] {
        for value in values {
            visitor(value)?;
        }
    }
    for value in [
        &mut step.crossing_heat_j,
        &mut step.reverse_crossing_heat_absorbed_j,
        &mut step.annihilation_events,
        &mut step.annihilation_singlets,
        &mut step.annihilation_heat_j,
    ] {
        visitor(value)?;
    }
    Ok(())
}

impl PairHistory for SaturableSingletTripletStep {
    fn populations(&self) -> Vec<[f64; 3]> {
        vec![self.populations]
    }
    fn visit_integrated(
        &mut self,
        visitor: &mut impl FnMut(&mut f64) -> Result<(), SingletTripletError>,
    ) -> Result<(), SingletTripletError> {
        visit_kinetics(self, visitor)
    }
}

pub(crate) fn energy_components(
    model: SingletTripletDecay,
    step: &SaturableSingletTripletStep,
    duration: f64,
) -> Result<Vec<f64>, SingletTripletError> {
    let mut values = Vec::new();
    for pair in [
        step.stored_energy_j,
        step.absorbed_pump_energy_j,
        step.radiative_decay_energy_j,
        step.nonradiative_heat_j,
        step.pump_relaxation_heat_j,
    ] {
        values.extend(pair);
    }
    values.extend([
        step.crossing_heat_j,
        step.reverse_crossing_heat_absorbed_j,
        step.annihilation_heat_j,
    ]);
    let energies = model.excitation_energies_j();
    values.push(positive_product(
        step.populations[0],
        energies[0],
        "ground convergence weight",
    )?);
    for target in 0..2 {
        values.push(product3(
            step.stored_energy_j[target],
            model.radiative_rates_s()[target],
            duration,
        )?);
        for events in [
            step.pump_events[target],
            step.radiative_events[target],
            step.nonradiative_events[target],
            step.crossing_events[target],
        ] {
            values.push(positive_product(
                events,
                energies[target],
                "event convergence weight",
            )?);
        }
    }
    values.push(positive_product(
        step.annihilation_events,
        energies[1],
        "pair convergence weight",
    )?);
    values.push(positive_product(
        step.annihilation_singlets,
        energies[0],
        "pair singlet convergence weight",
    )?);
    Ok(values)
}

#[cfg(test)]
#[path = "annihilation/oracles.rs"]
mod oracles;
