//! Positive exponential-midpoint evolution of photon recycling and local
//! triplet-triplet annihilation. The authored nonlinear vector field is frozen
//! using lambda = 2*k*Nt/V, then the ORIGINAL linear rates, crossing, pumping
//! and spectral photon return evolve together through the native positive owner.
//! A half-interval predictor supplies the midpoint density for a full corrector;
//! only corrector histories contribute to the returned physical accounts.
//!
//! Algebraically, the frozen triplet removal is split into a singlet-producing
//! rate lambda*phi/2 and a dark rate lambda*(1-phi/2). These are internal reaction
//! coordinates, not an authored reverse-crossing mechanism or thermal bath.
//! Original physical channels and pair events are reconstructed directly from
//! integrated triplet removals; artificial bath heat is never subtracted into
//! a reported physical channel.
//!
//! Grids double from 1 to 2 to 4. Endpoint stores AND integrated physical energy
//! channels must agree componentwise, including energy-weighted event counts.
//! The method is second order for resolved smooth nonlinear kinetics; doubling
//! is a convergence indicator, not a rigorous error enclosure. Fast ORIGINAL
//! lifetimes are handled by the native exponential rather than operator kicks.
//! Fast endpoint components can still exhibit stiff order reduction; a zero
//! absolute allowance demands relative accuracy even in tiny singlet stores.
//! Such a request may exhaust the explicit grid budget and is never relaxed.
//! The final base-grid slot is physically divided as h/2,h/4,...,delta,delta
//! until delta resolves the fastest ORIGINAL lifetime; delta also halves on
//! every base-grid refinement. This short graded tail
//! follows the local fast population at the actual endpoint without forcing
//! the entire grid onto that lifetime. Endpoint radiative power times the whole
//! interval is also compared, so tiny fast stores cannot conceal light errors.
//! Each frozen nonlinear exposure must be <=1, preventing an over-damped
//! predictor from causing false coarse/fine agreement. Exact isolated
//! annihilation and inactive-law linear evolution have analytic fast paths.
//! No extrapolation, normalization or clipping changes a positive candidate.

use super::super::annihilation::TripletAnnihilation;
use super::{
    CompensatedSum, SingletTripletError, SingletTripletRecyclingPlan, SingletTripletRecyclingStep,
    cap, finite, nonnegative, poll, positive_product, sums, total,
};

/// Accuracy of the refinement comparison and the largest permitted uniform grid.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AnnihilationLimits {
    /// Absolute allowance [J] for EACH compared energy component.
    pub absolute_energy_tolerance_j: f64,
    /// Relative allowance times the larger magnitude of that component.
    pub relative_tolerance: f64,
    /// Largest candidate grid in actual nonlinear segments, including the
    /// graded endpoint tail. Base grids double; each trial is independently
    /// bounded and cancellation remains available during all attempted work.
    /// Coupled composition needs two distinct grids; isolated annihilation uses one.
    pub maximum_substeps: usize,
}

impl Default for AnnihilationLimits {
    fn default() -> Self {
        Self {
            absolute_energy_tolerance_j: 0.0,
            relative_tolerance: 1e-8,
            maximum_substeps: 1024,
        }
    }
}

/// Whole-interval candidate with linear and annihilation channels kept distinct.
#[derive(Clone, Debug, PartialEq)]
pub struct AnnihilationRecyclingStep {
    /// Final populations/stores and summed LINEAR event/energy channels.
    /// In particular, nonradiative_heat_j excludes annihilation heat.
    /// Its residual is recomputed for the whole interval OMITTING annihilation
    /// heat; it therefore contains that heat plus the signed arithmetic residual.
    pub linear: SingletTripletRecyclingStep,
    /// Integrated triplet-pair events, in original active-cell order.
    pub annihilation_events: Vec<f64>,
    /// Integrated produced singlets, in original active-cell order.
    pub annihilation_singlets: Vec<f64>,
    /// Annihilation heat [J], separate from linear nonradiative decay.
    pub annihilation_heat_j: Vec<f64>,
    /// Accepted nonlinear substeps; one for exact isolated annihilation, zero for linear evolution.
    pub substeps: usize,
    /// Largest absolute energy-component difference from the preceding grid.
    /// This is not a rigorous bound on the accepted solution's error.
    pub max_energy_discrepancy_j: f64,
    /// Whole-interval external energy residual INCLUDING annihilation heat.
    pub residual_j: f64,
}

impl SingletTripletRecyclingPlan {
    /// Evolve fixed linear coefficients, constant pump rates, and authored local
    /// annihilation laws without mutating the supplied inventories.
    ///
    /// One optional law is required per active cell. Limits are checked even
    /// when duration or every pair rate is zero; those cases use the original
    /// full-duration linear plan exactly. Active laws retain their per-event
    /// energy admission. Exhausted refinement returns a bounded-resource error.
    ///
    /// Midpoint convergence uses the supplied tolerances. Independently, the
    /// whole-interval energy balance uses ONLY the prepared linear owner's
    /// relative tolerance; the integration allowance cannot weaken it. Only NET thermal
    /// exchange enters that balance scale, never gross internal photon traffic.
    pub fn advance_with_annihilation(
        &self,
        populations: &[[f64; 2]],
        pump_rates_s: &[[f64; 2]],
        laws: &[Option<TripletAnnihilation>],
        limits: AnnihilationLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AnnihilationRecyclingStep, SingletTripletError> {
        poll(&mut cancelled)?;
        if !nonnegative(limits.absolute_energy_tolerance_j)
            || !nonnegative(limits.relative_tolerance)
            || limits.relative_tolerance >= 1.0
            || (limits.absolute_energy_tolerance_j == 0.0 && limits.relative_tolerance == 0.0)
            || limits.maximum_substeps == 0
        {
            return Err(SingletTripletError::Invalid(
                "annihilation refinement limits",
            ));
        }
        let cells = self.cell_count();
        if laws.len() != cells
            || populations.len() != cells
            || pump_rates_s.len() != cells
            || !populations
                .iter()
                .chain(pump_rates_s)
                .flatten()
                .copied()
                .all(nonnegative)
        {
            return Err(SingletTripletError::Invalid(
                "annihilation cell laws, populations or pumping",
            ));
        }
        let mut active = false;
        for (law, model) in laws.iter().zip(self.models()) {
            poll(&mut cancelled)?;
            if let Some(law) = law.filter(|law| law.pair_rate_m3_s() > 0.0) {
                law.advance(model.excitation_energies_j(), [0.0; 2], 0.0)?;
                active = true;
            }
        }
        if !active || self.duration_s() == 0.0 {
            let linear = self.advance_with_cancel(populations, pump_rates_s, &mut cancelled)?;
            let residual_j = linear.residual_j;
            return Ok(AnnihilationRecyclingStep {
                linear,
                annihilation_events: vec![0.0; cells],
                annihilation_singlets: vec![0.0; cells],
                annihilation_heat_j: vec![0.0; cells],
                substeps: 0,
                max_energy_discrepancy_j: 0.0,
                residual_j,
            });
        }
        let isolated = pump_rates_s.iter().flatten().all(|&rate| rate == 0.0)
            && self.models().iter().all(|model| {
                model
                    .radiative_rates_s()
                    .into_iter()
                    .chain(model.nonradiative_rates_s())
                    .chain(model.crossing_rates_s())
                    .all(|rate| rate == 0.0)
            });
        if isolated {
            return isolated_step(self, populations, pump_rates_s, laws, &mut cancelled);
        }
        cap("annihilation uniform substeps", 2, limits.maximum_substeps)?;
        let mut substeps = 1usize;
        let mut previous: Option<(usize, Vec<f64>)> = None;
        loop {
            poll(&mut cancelled)?;
            cap(
                "annihilation uniform substeps",
                substeps,
                limits.maximum_substeps,
            )?;
            if let Some(mut candidate) = uniform(
                self,
                populations,
                pump_rates_s,
                laws,
                limits,
                substeps,
                &mut cancelled,
            )? {
                let current = energy_components(self, &candidate, &mut cancelled)?;
                if let Some((old_segments, old_components)) = &previous
                    && candidate.substeps > *old_segments
                    && !(substeps == 2 && *old_segments > 1)
                {
                    // Graded base grids 1 and 2 share their bulk partition;
                    // only the tiny final tail changes. That comparison cannot
                    // establish accuracy of the unresolved early interval.
                    let mut accepted = true;
                    let mut maximum = 0.0_f64;
                    for (index, (&old, &new)) in old_components.iter().zip(&current).enumerate() {
                        if index % 256 == 0 {
                            poll(&mut cancelled)?;
                        }
                        let discrepancy = (new - old).abs();
                        maximum = maximum.max(discrepancy);
                        accepted &= within(
                            discrepancy,
                            old.max(new),
                            limits.absolute_energy_tolerance_j,
                            limits.relative_tolerance,
                        );
                    }
                    if accepted {
                        candidate.max_energy_discrepancy_j = maximum;
                        poll(&mut cancelled)?;
                        return Ok(candidate);
                    }
                }
                previous = Some((candidate.substeps, current));
            } else {
                // No candidate from a grid that violates the nonlinear exposure
                // admission may be used as evidence for convergence.
                previous = None;
            }
            substeps = substeps
                .checked_mul(2)
                .ok_or(SingletTripletError::Unrepresentable(
                    "annihilation substep counter",
                ))?;
        }
    }
}

fn isolated_step(
    plan: &SingletTripletRecyclingPlan,
    initial: &[[f64; 2]],
    pump: &[[f64; 2]],
    laws: &[Option<TripletAnnihilation>],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<AnnihilationRecyclingStep, SingletTripletError> {
    let cells = plan.cell_count();
    let mut result = AnnihilationRecyclingStep {
        linear: plan.advance_with_cancel(initial, pump, &mut *cancelled)?,
        annihilation_events: vec![0.0; cells],
        annihilation_singlets: vec![0.0; cells],
        annihilation_heat_j: vec![0.0; cells],
        substeps: 1,
        max_energy_discrepancy_j: 0.0,
        residual_j: 0.0,
    };
    for (cell, law) in laws.iter().enumerate() {
        poll(cancelled)?;
        if let Some(law) = law.filter(|law| law.pair_rate_m3_s() > 0.0) {
            let energies = plan.models()[cell].excitation_energies_j();
            let step = law.advance(energies, initial[cell], plan.duration_s())?;
            result.linear.populations[cell] = step.populations;
            for state in 0..2 {
                result.linear.stored_energy_j[cell][state] = positive_product(
                    step.populations[state],
                    energies[state],
                    "final annihilation store",
                )?;
            }
            result.annihilation_events[cell] = step.pair_events;
            result.annihilation_singlets[cell] = step.singlet_production;
            result.annihilation_heat_j[cell] = step.nonradiative_heat_j;
        }
    }
    balance(plan, initial, &mut result)?;
    poll(cancelled)?;
    Ok(result)
}

#[allow(clippy::too_many_arguments)]
fn uniform(
    plan: &SingletTripletRecyclingPlan,
    initial: &[[f64; 2]],
    pump: &[[f64; 2]],
    laws: &[Option<TripletAnnihilation>],
    limits: AnnihilationLimits,
    substeps: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Option<AnnihilationRecyclingStep>, SingletTripletError> {
    let durations = segment_durations(plan, substeps, limits.maximum_substeps, cancelled)?;
    let actual_substeps = durations.len();
    let cells = plan.cell_count();
    let mut state = initial.to_vec();
    let mut integrated = Vec::new();
    let mut annihilation: Vec<[CompensatedSum; 3]> = (0..cells)
        .map(|_| {
            [
                CompensatedSum::default(),
                CompensatedSum::default(),
                CompensatedSum::default(),
            ]
        })
        .collect();
    let mut last = None;
    for duration in durations {
        poll(cancelled)?;
        let half = duration * 0.5;
        if half == 0.0 || half * 2.0 != duration {
            return Err(SingletTripletError::Unrepresentable(
                "annihilation predictor duration",
            ));
        }
        let initial_rates = frozen_rates(&state, laws)?;
        if initial_rates.iter().any(|&rate| rate * duration > 1.0) {
            return Ok(None);
        }
        let (_, predictor) =
            frozen_step(plan, &state, pump, laws, &initial_rates, half, cancelled)?;
        let midpoint_rates = frozen_rates(&predictor.populations, laws)?;
        if midpoint_rates.iter().any(|&rate| rate * duration > 1.0) {
            return Ok(None);
        }
        let (frozen, mut corrector) = frozen_step(
            plan,
            &state,
            pump,
            laws,
            &midpoint_rates,
            duration,
            cancelled,
        )?;
        for cell in 0..cells {
            poll(cancelled)?;
            let lambda = midpoint_rates[cell];
            if lambda == 0.0 {
                continue;
            }
            let law = laws[cell].ok_or(SingletTripletError::Invalid(
                "missing active annihilation law",
            ))?;
            let model = plan.models()[cell];
            let frozen_model = frozen.models()[cell];
            let total_rate = frozen_model.total_rates_s[1];
            let fractions = frozen_model.fractions[1];
            let events = [
                corrector.radiative_events[cell][1],
                corrector.nonradiative_events[cell][1],
                corrector.crossing_events[cell][1],
            ];
            let mut dominant = 0;
            for channel in 1..3 {
                if fractions[channel] > fractions[dominant] {
                    dominant = channel;
                }
            }
            // Dominant fraction >=1/3: recover native integrated removals
            // without an ill-conditioned subtraction or a tiny divisor.
            let removals = finite(
                events[dominant] / fractions[dominant],
                "integrated frozen triplet removals",
            )?;
            let pairs = positive_product(
                channel_events(removals, lambda, total_rate)?,
                0.5,
                "integrated annihilation pairs",
            )?;
            let (singlets, heat) = law.event_channels(model.excitation_energies_j(), pairs)?;
            let nonradiative =
                channel_events(removals, model.nonradiative_rates_s()[1], total_rate)?;
            let reverse = channel_events(removals, model.crossing_rates_s()[1], total_rate)?;
            let energies = model.excitation_energies_j();
            corrector.nonradiative_events[cell][1] = nonradiative;
            corrector.nonradiative_heat_j[cell][1] = positive_product(
                nonradiative,
                energies[1],
                "physical triplet nonradiative heat",
            )?;
            corrector.crossing_events[cell][1] = reverse;
            corrector.reverse_crossing_heat_absorbed_j[cell] = positive_product(
                reverse,
                energies[0] - energies[1],
                "physical reverse-crossing heat",
            )?;
            for (sum, value) in annihilation[cell].iter_mut().zip([pairs, singlets, heat]) {
                sum.add(value)?;
            }
        }
        accumulate(&mut integrated, &corrector, cancelled)?;
        state.copy_from_slice(&corrector.populations);
        last = Some(corrector);
    }
    let mut linear = last.ok_or(SingletTripletError::Invalid("empty annihilation grid"))?;
    for (index, (output, sum)) in integrated_mut(&mut linear).zip(&integrated).enumerate() {
        if index % 256 == 0 {
            poll(cancelled)?;
        }
        *output = total(sum, "integrated linear channel")?;
    }
    let mut result = AnnihilationRecyclingStep {
        linear,
        annihilation_events: Vec::with_capacity(cells),
        annihilation_singlets: Vec::with_capacity(cells),
        annihilation_heat_j: Vec::with_capacity(cells),
        substeps: actual_substeps,
        max_energy_discrepancy_j: 0.0,
        residual_j: 0.0,
    };
    for channels in annihilation {
        result
            .annihilation_events
            .push(total(&channels[0], "integrated annihilation pairs")?);
        result
            .annihilation_singlets
            .push(total(&channels[1], "integrated annihilation singlets")?);
        result
            .annihilation_heat_j
            .push(total(&channels[2], "integrated annihilation heat")?);
    }
    balance(plan, initial, &mut result)?;
    poll(cancelled)?;
    Ok(Some(result))
}

fn segment_durations(
    plan: &SingletTripletRecyclingPlan,
    base_substeps: usize,
    maximum: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<f64>, SingletTripletError> {
    let duration = plan.duration_s() / base_substeps as f64;
    if duration == 0.0 || duration * base_substeps as f64 != plan.duration_s() {
        return Err(SingletTripletError::Unrepresentable(
            "annihilation split duration",
        ));
    }
    let fastest = plan
        .models()
        .iter()
        .flat_map(|model| model.total_rates_s)
        .fold(0.0_f64, f64::max);
    let mut tail = Vec::new();
    let mut remaining = duration;
    // Keep the dyadic tail depth tied to the WHOLE interval. Refining the
    // base grid then refines delta as well; no fixed tail error is shared by
    // every successive candidate.
    while (remaining * base_substeps as f64) * fastest > 1.0 {
        poll(cancelled)?;
        let half = remaining * 0.5;
        if half == 0.0 || half * 2.0 != remaining {
            return Err(SingletTripletError::Unrepresentable(
                "annihilation endpoint-tail duration",
            ));
        }
        remaining = half;
        tail.push(remaining);
        let required =
            base_substeps
                .checked_add(tail.len())
                .ok_or(SingletTripletError::Unrepresentable(
                    "annihilation tail work counter",
                ))?;
        cap(
            "annihilation segments including endpoint tail",
            required,
            maximum,
        )?;
    }
    // h/2+h/4+...+delta+delta is exactly the original final slot in real
    // arithmetic. Every halving above is admitted only when binary64 retains it.
    let mut durations = vec![duration; base_substeps - 1];
    durations.extend(tail);
    durations.push(remaining);
    Ok(durations)
}

fn frozen_rates(
    populations: &[[f64; 2]],
    laws: &[Option<TripletAnnihilation>],
) -> Result<Vec<f64>, SingletTripletError> {
    populations
        .iter()
        .zip(laws)
        .map(|(state, law)| law.map_or(Ok(0.0), |law| law.frozen_removal_rate_s(state[1])))
        .collect()
}

fn frozen_step(
    plan: &SingletTripletRecyclingPlan,
    populations: &[[f64; 2]],
    pump: &[[f64; 2]],
    laws: &[Option<TripletAnnihilation>],
    rates: &[f64],
    duration: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(SingletTripletRecyclingPlan, SingletTripletRecyclingStep), SingletTripletError> {
    let mut models = plan.models().to_vec();
    for (cell, &lambda) in rates.iter().enumerate() {
        poll(cancelled)?;
        if lambda == 0.0 {
            continue;
        }
        let law = laws[cell].ok_or(SingletTripletError::Invalid(
            "missing frozen annihilation law",
        ))?;
        let model = models[cell];
        let mut nonradiative = model.nonradiative_rates_s();
        let mut crossing = model.crossing_rates_s();
        let singlet_rate = positive_product(
            positive_product(
                lambda,
                law.singlet_yield(),
                "frozen annihilation singlet rate",
            )?,
            0.5,
            "frozen annihilation singlet rate",
        )?;
        let dark_rate = positive_product(
            lambda,
            1.0 - 0.5 * law.singlet_yield(),
            "frozen annihilation dark rate",
        )?;
        nonradiative[1] = finite(
            nonradiative[1] + dark_rate,
            "combined frozen triplet dark rate",
        )?;
        crossing[1] = finite(
            crossing[1] + singlet_rate,
            "combined frozen triplet transfer rate",
        )?;
        models[cell] = super::SingletTripletDecay::try_new(
            model.excitation_energies_j(),
            model.radiative_rates_s(),
            nonradiative,
            crossing,
        )?;
    }
    let frozen = SingletTripletRecyclingPlan::prepare(
        &models,
        plan.channels(),
        plan.absorption_edges(),
        duration,
        plan.fraction_roundoff,
        plan.limits,
        &mut *cancelled,
    )?;
    let step = frozen.advance_with_cancel(populations, pump, &mut *cancelled)?;
    Ok((frozen, step))
}

fn channel_events(removals: f64, rate: f64, total_rate: f64) -> Result<f64, SingletTripletError> {
    let fraction = rate / total_rate;
    if rate > 0.0 && fraction == 0.0 {
        return Err(SingletTripletError::Unrepresentable(
            "positive original competing fraction",
        ));
    }
    positive_product(removals, fraction, "physical competing events")
}

// The endpoint populations/stores and the residual are intentionally absent.
// The two enumerations share this order so every integrated field is compensated.
fn integrated_values(step: &SingletTripletRecyclingStep) -> impl Iterator<Item = &f64> {
    [
        step.injected_excitations.as_flattened(),
        step.radiative_events.as_flattened(),
        step.nonradiative_events.as_flattened(),
        step.crossing_events.as_flattened(),
        step.pumped_energy_j.as_flattened(),
        step.radiative_decay_energy_j.as_flattened(),
        step.nonradiative_heat_j.as_flattened(),
        step.relaxation_heat_j.as_flattened(),
        step.recycled_excitations.as_flattened(),
        step.recycled_energy_j.as_flattened(),
        &step.crossing_heat_j,
        &step.reverse_crossing_heat_absorbed_j,
        &step.emitted_energy_j,
        &step.escaped_energy_j,
        &step.unresolved_energy_j,
        &step.absorbed_energy_j,
        &step.absorption_heat_j,
    ]
    .into_iter()
    .flatten()
}
fn integrated_mut(step: &mut SingletTripletRecyclingStep) -> impl Iterator<Item = &mut f64> {
    [
        step.injected_excitations.as_flattened_mut(),
        step.radiative_events.as_flattened_mut(),
        step.nonradiative_events.as_flattened_mut(),
        step.crossing_events.as_flattened_mut(),
        step.pumped_energy_j.as_flattened_mut(),
        step.radiative_decay_energy_j.as_flattened_mut(),
        step.nonradiative_heat_j.as_flattened_mut(),
        step.relaxation_heat_j.as_flattened_mut(),
        step.recycled_excitations.as_flattened_mut(),
        step.recycled_energy_j.as_flattened_mut(),
        &mut step.crossing_heat_j,
        &mut step.reverse_crossing_heat_absorbed_j,
        &mut step.emitted_energy_j,
        &mut step.escaped_energy_j,
        &mut step.unresolved_energy_j,
        &mut step.absorbed_energy_j,
        &mut step.absorption_heat_j,
    ]
    .into_iter()
    .flatten()
}
fn accumulate(
    integrated: &mut Vec<CompensatedSum>,
    step: &SingletTripletRecyclingStep,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), SingletTripletError> {
    if integrated.is_empty() {
        *integrated = sums(integrated_values(step).count());
    }
    for (index, (sum, value)) in integrated
        .iter_mut()
        .zip(integrated_values(step))
        .enumerate()
    {
        if index % 256 == 0 {
            poll(cancelled)?;
        }
        sum.add(*value)?;
    }
    Ok(())
}

fn energy_components(
    plan: &SingletTripletRecyclingPlan,
    candidate: &AnnihilationRecyclingStep,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<f64>, SingletTripletError> {
    let step = &candidate.linear;
    let mut values = Vec::with_capacity(integrated_values(step).count() + 7 * plan.cell_count());
    for (cell, model) in plan.models().iter().enumerate() {
        poll(cancelled)?;
        let energies = model.excitation_energies_j();
        for state in 0..2 {
            values.push(step.stored_energy_j[cell][state]);
            let endpoint_radiative_power = positive_product(
                step.stored_energy_j[cell][state],
                model.radiative_rates_s()[state],
                "endpoint radiative power",
            )?;
            values.push(positive_product(
                endpoint_radiative_power,
                plan.duration_s(),
                "endpoint radiative interval energy",
            )?);
            for counts in [
                &step.injected_excitations,
                &step.radiative_events,
                &step.nonradiative_events,
                &step.crossing_events,
                &step.recycled_excitations,
            ] {
                values.push(positive_product(
                    counts[cell][state],
                    energies[state],
                    "energy-weighted linear events",
                )?);
            }
            for channel in [
                &step.pumped_energy_j,
                &step.radiative_decay_energy_j,
                &step.nonradiative_heat_j,
                &step.relaxation_heat_j,
                &step.recycled_energy_j,
            ] {
                values.push(channel[cell][state]);
            }
        }
        values.push(step.crossing_heat_j[cell]);
        values.push(step.reverse_crossing_heat_absorbed_j[cell]);
        let pair_energy = if energies[1] <= f64::MAX * 0.5 {
            positive_product(
                candidate.annihilation_events[cell],
                2.0 * energies[1],
                "energy-weighted annihilation pairs",
            )?
        } else {
            positive_product(
                positive_product(
                    candidate.annihilation_events[cell],
                    energies[1],
                    "energy-weighted annihilation pairs",
                )?,
                2.0,
                "energy-weighted annihilation pairs",
            )?
        };
        values.push(pair_energy);
        values.push(positive_product(
            candidate.annihilation_singlets[cell],
            energies[0],
            "energy-weighted annihilation singlets",
        )?);
        values.push(candidate.annihilation_heat_j[cell]);
    }
    for (index, value) in step
        .emitted_energy_j
        .iter()
        .chain(&step.escaped_energy_j)
        .chain(&step.unresolved_energy_j)
        .chain(&step.absorbed_energy_j)
        .chain(&step.absorption_heat_j)
        .enumerate()
    {
        if index % 256 == 0 {
            poll(cancelled)?;
        }
        values.push(*value);
    }
    Ok(values)
}

fn balance(
    plan: &SingletTripletRecyclingPlan,
    initial: &[[f64; 2]],
    candidate: &mut AnnihilationRecyclingStep,
) -> Result<(), SingletTripletError> {
    let step = &mut candidate.linear;
    let mut external = CompensatedSum::default();
    let mut scale = CompensatedSum::default();
    let mut thermal = CompensatedSum::default();
    for (cell, model) in plan.models().iter().enumerate() {
        let energies = model.excitation_energies_j();
        for state in 0..2 {
            for term in [
                positive_product(
                    initial[cell][state],
                    energies[state],
                    "initial excitation energy",
                )?,
                step.pumped_energy_j[cell][state],
                -step.stored_energy_j[cell][state],
            ] {
                external.add(term)?;
                scale.add(term.abs())?;
            }
            thermal.add(step.nonradiative_heat_j[cell][state])?;
            thermal.add(step.relaxation_heat_j[cell][state])?;
        }
        thermal.add(step.crossing_heat_j[cell])?;
        thermal.add(-step.reverse_crossing_heat_absorbed_j[cell])?;
    }
    for &heat in &step.absorption_heat_j {
        thermal.add(heat)?;
    }
    for &loss in step
        .escaped_energy_j
        .iter()
        .chain(&step.unresolved_energy_j)
    {
        external.add(-loss)?;
        scale.add(loss)?;
    }
    step.residual_j = difference(&external, &thermal)?;
    for &heat in &candidate.annihilation_heat_j {
        thermal.add(heat)?;
    }
    candidate.residual_j = difference(&external, &thermal)?;
    scale.add(total(&thermal, "net coupled thermal exchange")?.abs())?;
    let scale = total(&scale, "coupled external energy scale")?;
    let allowance_j = finite(
        plan.limits.relative_tolerance * scale,
        "coupled balance allowance",
    )?;
    if candidate.residual_j.abs() > allowance_j {
        return Err(SingletTripletError::EnergyBalance {
            residual_j: candidate.residual_j,
            allowance_j,
        });
    }
    Ok(())
}

fn difference(a: &CompensatedSum, b: &CompensatedSum) -> Result<f64, SingletTripletError> {
    let mut sum = CompensatedSum::default();
    for term in [a.value, -b.value, a.correction, -b.correction] {
        sum.add(term)?;
    }
    total(&sum, "coupled energy residual")
}
// Rearranged comparison avoids overflowing abs + rel*scale when an explicitly
// generous finite absolute allowance already admits the finite discrepancy.
fn within(error: f64, scale: f64, absolute: f64, relative: f64) -> bool {
    error <= absolute || error - absolute <= relative * scale
}

#[cfg(test)]
#[path = "annihilation/tests.rs"]
mod tests;
