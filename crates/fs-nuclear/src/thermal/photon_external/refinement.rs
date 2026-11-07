//! Whole-trajectory time refinement of native photon/thermal feedback.
//!
//! Material transport and caloric evolution retain their original owners.
//! Comparisons measure final native state and accepted photon-energy histories;
//! conservation roundoff and spatial transport accuracy are separate questions.

use super::*;

/// Temporal agreement and resources for the complete candidate ladder.
/// All absolute tolerances are positive. Relative scales use changes from the
/// input native state and photon energy in this call, never historical totals
/// or the arbitrary zero of a caloric law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhotonTimeRefinement {
    pub initial_max_dt_s: f64,
    pub minimum_max_dt_s: f64,
    /// At least two: every nonzero advance compares coarse, fine and finer runs.
    pub maximum_halvings: usize,
    /// Total native steps over all candidates, including discarded coarse work.
    /// The original thermal per-advance budget also remains binding.
    pub maximum_total_substeps: usize,
    /// Rebuilds across every candidate; the initial retained response is reused.
    /// The original feedback per-advance preparation limit is also enforced.
    pub maximum_total_transport_preparations: usize,
    /// Original transport work across all candidates, including nested spatial
    /// refinement. Zero is useful for fixed-response or unchanged-temperature
    /// paths: it forbids new solves, not cache use.
    pub maximum_total_transport_work: usize,
    pub specific_energy_absolute_tolerance_j_kg: f64,
    pub temperature_absolute_tolerance_k: f64,
    /// Applies separately to each receiver's credited energy and each incident,
    /// absorbed, escaped and thermal-credit account.
    pub photon_energy_absolute_tolerance_j: f64,
    pub relative_tolerance: f64,
    /// Original thermal conservation allowance for each complete trajectory.
    /// Photon schedule/transport roundoff retains its original owner budget.
    pub conservation_tolerance_j: f64,
}

/// Measured temporal discrepancies; these are not rigorous global error bounds.
/// No extrapolated state, repaired energy account or resampled history is used.
#[derive(Clone, Debug, PartialEq)]
pub struct PhotonTimeConvergence {
    pub thermal: EnergyAudit,
    /// Accepted photon credits per original receiver over this call [J].
    /// Other boundary/conduction heat is not included.
    pub receiver_thermal_credit_j: BTreeMap<CellId, f64>,
    pub halvings: usize,
    pub total_substeps: usize,
    pub total_transport_preparations: usize,
    pub total_transport_work: usize,
    pub accepted_max_dt_s: f64,
    pub largest_accepted_interval_s: f64,
    pub maximum_specific_energy_difference_j_kg: f64,
    pub maximum_temperature_difference_k: f64,
    pub maximum_receiver_energy_difference_j: f64,
    pub maximum_photon_account_difference_j: f64,
    pub normalized_difference: f64,
    /// Every compared coordinate differs by at most 64 binary64 epsilons of
    /// its represented magnitude. Energy scales cover only this call. This
    /// relaxes only the decreasing-difference gate, never the caller tolerance.
    pub roundoff_limited: bool,
    /// log2(previous/final normalized discrepancy), when both actual interval
    /// maxima halved and both discrepancies are resolved above roundoff.
    /// Nonsmooth material laws/events need not exhibit first-order convergence.
    pub observed_order: Option<f64>,
}

#[derive(Clone, Copy, Debug)]
struct Difference {
    specific_energy: f64,
    temperature: f64,
    receiver_energy: f64,
    photon_account: f64,
    normalized: f64,
    roundoff_limited: bool,
}
impl Default for Difference {
    fn default() -> Self {
        Self {
            specific_energy: 0.0,
            temperature: 0.0,
            receiver_energy: 0.0,
            photon_account: 0.0,
            normalized: 0.0,
            roundoff_limited: true,
        }
    }
}
impl Difference {
    fn include(
        &mut self,
        delta: f64,
        change_scale: f64,
        represented_scale: f64,
        absolute: f64,
        relative: f64,
    ) -> Result<(), FissionThermalError> {
        let delta = finite(delta)?;
        let change_scale = finite(change_scale)?;
        let represented_scale = finite(represented_scale)?;
        let allowance = finite(absolute + relative * change_scale)?;
        self.normalized = self.normalized.max(delta / allowance);
        self.roundoff_limited &= delta <= 64.0 * f64::EPSILON * represented_scale;
        Ok(())
    }
}

#[derive(Default)]
struct Work {
    substeps: usize,
    preparations: usize,
    transport: usize,
}

impl ExternalPhotonHeating {
    /// Refine complete original thermal/photon trajectories from one input
    /// state, including native events, source knots, conduction and calorics.
    /// Accept only the finest original trajectory after successive differences
    /// meet the declared tolerances and decrease (or both reach roundoff).
    ///
    /// Every new rung halves the previous largest ACTUAL native interval,
    /// additionally bounded by its requested maximum. Thus a large requested
    /// step already clipped by CFL or events cannot produce three identical
    /// time partitions and masquerade as successful refinement.
    ///
    /// Physical-domain, numerical-accounting, cancellation and resource errors
    /// propagate immediately. Only excessive temporal disagreement requests
    /// another rung. No failed candidate mutates either input owner. Work across
    /// discarded candidates is charged, while trial resource limits never
    /// replace the returned source's original budgets.
    ///
    /// This estimates final state and receiver/account energy agreement for the
    /// prescribed geometry, spectral grid and chosen spatial transport solve.
    /// It does not establish spatial/angular/spectral accuracy or measured-law
    /// validity. A zero-duration call retains native event handling without a
    /// fictitious refinement ladder.
    pub fn advance_refined<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        controls: PhotonTimeRefinement,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(DepositionSimulation<'m>, Self, PhotonTimeConvergence), FissionThermalError> {
        checkpoint()?;
        self.validate(thermal, self.time_s, &mut checkpoint)?;
        validate_controls(controls)?;
        let mut work = Work::default();
        let mut dt_s = controls.initial_max_dt_s;
        let mut previous = run(
            thermal,
            self,
            target_s,
            dt_s,
            controls,
            &mut work,
            &mut checkpoint,
        )?;
        if target_s == thermal.time_s() {
            let report = report(
                &previous,
                Difference::default(),
                0,
                dt_s,
                &work,
                None,
                &mut checkpoint,
            )?;
            return Ok((previous.thermal, previous.source, report));
        }
        let mut previous_difference: Option<(Difference, bool)> = None;
        for halving in 1..=controls.maximum_halvings {
            let next_dt = 0.5 * dt_s.min(previous.largest_interval_s);
            if !positive(next_dt) || next_dt < controls.minimum_max_dt_s {
                return Err(convergence_failure(previous_difference, dt_s));
            }
            dt_s = next_dt;
            let current = run(
                thermal,
                self,
                target_s,
                dt_s,
                controls,
                &mut work,
                &mut checkpoint,
            )?;
            if !positive(current.largest_interval_s)
                || current.largest_interval_s >= previous.largest_interval_s
            {
                return Err(FissionThermalError::Invalid(
                    "photon refinement did not shorten actual native intervals",
                ));
            }
            let ratio = previous.largest_interval_s / current.largest_interval_s;
            let halved = (ratio - 2.0).abs() <= 128.0 * f64::EPSILON;
            let difference = compare(thermal, &previous, &current, controls, &mut checkpoint)?;
            if let Some((prior, prior_halved)) = previous_difference {
                let decreasing =
                    difference.normalized == 0.0 || difference.normalized < prior.normalized;
                let unresolved = prior.roundoff_limited && difference.roundoff_limited;
                if difference.normalized <= 1.0 && (decreasing || unresolved) {
                    let order = if prior_halved
                        && halved
                        && !prior.roundoff_limited
                        && !difference.roundoff_limited
                        && prior.normalized.is_finite()
                        && prior.normalized > 0.0
                        && difference.normalized > 0.0
                    {
                        Some(
                            (fs_math::det::ln(prior.normalized)
                                - fs_math::det::ln(difference.normalized))
                                / fs_math::det::ln(2.0),
                        )
                    } else {
                        None
                    };
                    let report = report(
                        &current,
                        difference,
                        halving,
                        dt_s,
                        &work,
                        order,
                        &mut checkpoint,
                    )?;
                    return Ok((current.thermal, current.source, report));
                }
            }
            previous_difference = Some((difference, halved));
            previous = current;
        }
        Err(convergence_failure(previous_difference, dt_s))
    }
}

fn validate_controls(c: PhotonTimeRefinement) -> Result<(), FissionThermalError> {
    if !positive(c.initial_max_dt_s)
        || !positive(c.minimum_max_dt_s)
        || c.minimum_max_dt_s > c.initial_max_dt_s
        || !(2..=64).contains(&c.maximum_halvings)
        || !positive(c.specific_energy_absolute_tolerance_j_kg)
        || !positive(c.temperature_absolute_tolerance_k)
        || !positive(c.photon_energy_absolute_tolerance_j)
        || !c.relative_tolerance.is_finite()
        || !(0.0..1.0).contains(&c.relative_tolerance)
        || !nonnegative(c.conservation_tolerance_j)
    {
        return Err(FissionThermalError::Invalid(
            "positive photon time/error bounds, two to 64 halvings and finite relative tolerance in [0,1) required",
        ));
    }
    Ok(())
}

fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}

fn convergence_failure(
    previous: Option<(Difference, bool)>,
    timestep_s: f64,
) -> FissionThermalError {
    FissionThermalError::PhotonTimeConvergence {
        normalized_difference: previous.map_or(f64::INFINITY, |p| p.0.normalized),
        timestep_s,
    }
}

// This wrapper only observes the original source's accepted watts. It does not
// reapply a response, resample a material or independently integrate temperature.
#[derive(Clone)]
struct ReceiverObserver {
    source: ExternalPhotonHeating,
    receiver_credit: BTreeMap<CellId, Sum>,
    pending_credit: Option<Vec<(CellId, f64)>>,
    largest_interval_s: f64,
}
impl DepositionObserver for ReceiverObserver {
    type Error = FissionThermalError;

    fn begin(&mut self) -> Result<(), Self::Error> {
        self.source.begin()?;
        if self.pending_credit.is_some() {
            return Err(FissionThermalError::Invalid(
                "uncommitted photon receiver credit",
            ));
        }
        self.receiver_credit
            .values_mut()
            .for_each(|v| *v = Sum::default());
        self.largest_interval_s = 0.0;
        Ok(())
    }

    fn requires_step_interval(&self) -> bool {
        true
    }

    fn next_transition_time_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Option<f64>, Self::Error> {
        self.source.next_transition_time_s(thermal, checkpoint)
    }

    fn step_sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        start_s: f64,
        duration_s: f64,
        end_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        let watts = self
            .source
            .step_sources(thermal, start_s, duration_s, end_s, checkpoint)?;
        let mut credits = Vec::with_capacity(watts.len());
        for (&id, &power) in &watts {
            checkpoint()?;
            credits.push((id, finite(power * duration_s)?));
        }
        self.pending_credit = Some(credits);
        Ok(watts)
    }

    fn event(
        &mut self,
        event: &Event,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<(), Self::Error> {
        self.source.event(event, thermal)
    }

    fn step(
        &mut self,
        start_s: f64,
        duration_s: f64,
        old_temperatures_k: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.source
            .step(start_s, duration_s, old_temperatures_k, thermal, checkpoint)?;
        let credits = self
            .pending_credit
            .take()
            .ok_or(FissionThermalError::Invalid(
                "photon refinement needs the original prepared cell credits",
            ))?;
        for (id, energy) in credits {
            checkpoint()?;
            self.receiver_credit
                .get_mut(&id)
                .ok_or(FissionThermalError::Invalid(
                    "photon refinement receiver coverage changed",
                ))?
                .add(energy)?;
        }
        self.largest_interval_s = self.largest_interval_s.max(duration_s);
        Ok(())
    }
}

struct Candidate<'m> {
    thermal: DepositionSimulation<'m>,
    source: ExternalPhotonHeating,
    thermal_audit: EnergyAudit,
    receiver_credit: BTreeMap<CellId, Sum>,
    largest_interval_s: f64,
}

#[allow(clippy::too_many_arguments)]
fn run<'m>(
    thermal: &DepositionSimulation<'m>,
    source: &ExternalPhotonHeating,
    target_s: f64,
    dt_s: f64,
    controls: PhotonTimeRefinement,
    work: &mut Work,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<Candidate<'m>, FissionThermalError> {
    let remaining_steps = controls
        .maximum_total_substeps
        .checked_sub(work.substeps)
        .ok_or(RadiationError::Budget("photon refinement total substeps"))?;
    let remaining = PhotonFeedbackWork {
        preparations: controls
            .maximum_total_transport_preparations
            .checked_sub(work.preparations)
            .ok_or(RadiationError::Budget(
                "photon refinement total preparations",
            ))?,
        transport_work: controls
            .maximum_total_transport_work
            .checked_sub(work.transport)
            .ok_or(RadiationError::Budget(
                "photon refinement total transport work",
            ))?,
    };
    let mut source = source.clone();
    if let Some(feedback) = &mut source.feedback {
        feedback.limit_trial_work(Some(remaining));
    }
    let observer = ReceiverObserver {
        receiver_credit: source
            .receivers
            .keys()
            .map(|&id| (id, Sum::default()))
            .collect(),
        source,
        pending_credit: None,
        largest_interval_s: 0.0,
    };
    let (thermal, mut observer, thermal_audit) = thermal.advance_observed_limited(
        &observer,
        target_s,
        StepConfig {
            max_dt_s: dt_s,
            energy_tolerance_j: controls.conservation_tolerance_j,
        },
        remaining_steps,
        &mut *checkpoint,
    )?;
    observer.source.validate(&thermal, target_s, checkpoint)?;
    observer.source.admit_audit(observer.source.audit)?;
    work.substeps = work
        .substeps
        .checked_add(thermal_audit.substeps)
        .filter(|&n| n <= controls.maximum_total_substeps)
        .ok_or(RadiationError::Budget("photon refinement total substeps"))?;
    if let Some(used) = observer.source.feedback_work() {
        work.preparations = work
            .preparations
            .checked_add(used.preparations)
            .filter(|&n| n <= controls.maximum_total_transport_preparations)
            .ok_or(RadiationError::Budget(
                "photon refinement total preparations",
            ))?;
        work.transport = work
            .transport
            .checked_add(used.transport_work)
            .filter(|&n| n <= controls.maximum_total_transport_work)
            .ok_or(RadiationError::Budget(
                "photon refinement total transport work",
            ))?;
    }
    if let Some(feedback) = &mut observer.source.feedback {
        feedback.limit_trial_work(None);
    }
    checkpoint()?;
    Ok(Candidate {
        thermal,
        source: observer.source,
        thermal_audit,
        receiver_credit: observer.receiver_credit,
        largest_interval_s: observer.largest_interval_s,
    })
}

fn compare(
    initial: &DepositionSimulation<'_>,
    coarse: &Candidate<'_>,
    fine: &Candidate<'_>,
    controls: PhotonTimeRefinement,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<Difference, FissionThermalError> {
    if coarse.thermal.time_s() != fine.thermal.time_s()
        || coarse.thermal.cells().len() != fine.thermal.cells().len()
        || coarse.receiver_credit.len() != fine.receiver_credit.len()
    {
        return Err(FissionThermalError::Invalid(
            "photon refinement candidates disagree on native endpoint or coverage",
        ));
    }
    let mut difference = Difference::default();
    for (&id, f) in fine.thermal.cells() {
        checkpoint()?;
        let c = coarse
            .thermal
            .cells()
            .get(&id)
            .ok_or(FissionThermalError::Invalid(
                "photon refinement cell mismatch",
            ))?;
        let dh = finite((f.specific_h_j_kg - c.specific_h_j_kg).abs())?;
        // A birth has no accepted input datum; its absolute allowance applies.
        let h_scale = initial.cells().get(&id).map_or(0.0, |base| {
            (f.specific_h_j_kg - base.specific_h_j_kg)
                .abs()
                .max((c.specific_h_j_kg - base.specific_h_j_kg).abs())
        });
        difference.specific_energy = difference.specific_energy.max(dh);
        difference.include(
            dh,
            h_scale,
            f.specific_h_j_kg.abs().max(c.specific_h_j_kg.abs()),
            controls.specific_energy_absolute_tolerance_j_kg,
            controls.relative_tolerance,
        )?;
        let tf = fine.thermal.temperature_k(id)?;
        let tc = coarse.thermal.temperature_k(id)?;
        let dt = finite((tf - tc).abs())?;
        let t_scale = if initial.cells().contains_key(&id) {
            let t0 = initial.temperature_k(id)?;
            (tf - t0).abs().max((tc - t0).abs())
        } else {
            0.0
        };
        difference.temperature = difference.temperature.max(dt);
        difference.include(
            dt,
            t_scale,
            tf.abs().max(tc.abs()),
            controls.temperature_absolute_tolerance_k,
            controls.relative_tolerance,
        )?;
    }
    for (&id, fine_energy) in &fine.receiver_credit {
        checkpoint()?;
        let coarse_energy = coarse
            .receiver_credit
            .get(&id)
            .ok_or(FissionThermalError::Invalid(
                "photon refinement receiver mismatch",
            ))?;
        let change = compare_energy(*coarse_energy, *fine_energy, &mut difference, controls)?;
        difference.receiver_energy = difference.receiver_energy.max(change);
    }
    let c = &coarse.source.account;
    let f = &fine.source.account;
    for (a, b) in [
        (c.incident, f.incident),
        (c.absorbed, f.absorbed),
        (c.escaped, f.escaped),
        (c.credit, f.credit),
    ] {
        checkpoint()?;
        let change = compare_energy(a, b, &mut difference, controls)?;
        difference.photon_account = difference.photon_account.max(change);
    }
    Ok(difference)
}

fn compare_energy(
    coarse: Sum,
    fine: Sum,
    difference: &mut Difference,
    controls: PhotonTimeRefinement,
) -> Result<f64, FissionThermalError> {
    let mut delta = Sum::default();
    delta.add(fine.sum - coarse.sum)?;
    delta.add(fine.correction)?;
    delta.add(-coarse.correction)?;
    let change = delta.value()?.abs();
    // The original observer resets these accounts once per complete candidate.
    let scale = coarse.value()?.abs().max(fine.value()?.abs());
    difference.include(
        change,
        scale,
        scale,
        controls.photon_energy_absolute_tolerance_j,
        controls.relative_tolerance,
    )?;
    Ok(change)
}

#[allow(clippy::too_many_arguments)]
fn report(
    candidate: &Candidate<'_>,
    difference: Difference,
    halvings: usize,
    dt_s: f64,
    work: &Work,
    order: Option<f64>,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<PhotonTimeConvergence, FissionThermalError> {
    let mut receiver_thermal_credit_j = BTreeMap::new();
    for (&id, sum) in &candidate.receiver_credit {
        checkpoint()?;
        receiver_thermal_credit_j.insert(id, sum.value()?);
    }
    checkpoint()?;
    Ok(PhotonTimeConvergence {
        thermal: candidate.thermal_audit,
        receiver_thermal_credit_j,
        halvings,
        total_substeps: work.substeps,
        total_transport_preparations: work.preparations,
        total_transport_work: work.transport,
        accepted_max_dt_s: dt_s,
        largest_accepted_interval_s: candidate.largest_interval_s,
        maximum_specific_energy_difference_j_kg: difference.specific_energy,
        maximum_temperature_difference_k: difference.temperature,
        maximum_receiver_energy_difference_j: difference.receiver_energy,
        maximum_photon_account_difference_j: difference.photon_account,
        normalized_difference: difference.normalized,
        roundoff_limited: difference.roundoff_limited,
        observed_order: order,
    })
}
