//! Adaptive pressure/recession feedback along a prescribed local sliding path.
//!
//! Each cell follows `s_i(lambda) = total_distance_i * lambda`, 0 <= lambda <= 1.
//! The original k/H sources and the selected normal boundary condition stay
//! fixed. Full versus two-half-step comparisons estimate local discretization
//! error. Only the actual two-half-step history is accepted; extrapolating
//! irreversible state could violate positivity, pressure caps, or plastic memory.
//!
//! This controls the beginning-pressure quadrature used by the native wear
//! owner. It does not derive slip from friction, supply a global error bound,
//! resolve debris transport, or turn parameter work into heat. Arbitrary changes
//! of the spatial sliding pattern require successive prescribed path segments.

use fs_blake3::ContentHash;
use fs_exec::Cx;

use crate::normal_patch::{
    FiniteGapContactError, FiniteGapPlasticDriveError, FiniteGapPlasticError,
    FiniteGapPlasticLoadControls, FiniteGapPlasticLoadError,
};

use super::{
    SourcePlasticArchardWear, SourcePlasticContactError, SourcePlasticDriveError,
    SourcePlasticLoadError, SourcePlasticWearAdvance, SourcePlasticWearError,
    SourcePlasticWearState, add_compensated, checkpoint,
};

/// Local error thresholds and finite work limits for one complete sliding path.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WearIntegrationControls {
    /// Positive absolute tolerance for recession, plastic depth, and rigid approach [m].
    pub absolute_depth_tolerance_m: f64,
    /// Positive absolute tolerance for each cell's pressure [Pa].
    pub absolute_pressure_tolerance_pa: f64,
    /// Nonnegative relative tolerance, applied per cell. Depth scales use the
    /// current interval's increments, never cumulative wear from earlier paths.
    pub relative_tolerance: f64,
    /// Maximum accepted native wear substeps. A nonzero interval accepts two.
    pub maximum_accepted_substeps: usize,
    /// Maximum native wear evaluations, including every rejected trial.
    /// A load-controlled evaluation also obeys its own normal-solve budget.
    pub maximum_step_evaluations: usize,
}

/// Refusal never publishes a partial sliding path or changes the predecessor.
#[derive(Clone, Debug, PartialEq)]
pub enum WearIntegrationError {
    /// Invalid path, tolerance, or finite arithmetic.
    Invalid(&'static str),
    /// Original source or native wear/contact law refused the step.
    Step(SourcePlasticWearError),
    /// The requested complete path needs more admitted work than the caller allowed.
    BudgetExceeded {
        /// Accepted fine substeps in the discarded candidate path.
        accepted_substeps: usize,
        /// All wear evaluations already attempted.
        step_evaluations: usize,
        /// Most recent step-dependent physical refusal, if any.
        last_refusal: Option<SourcePlasticWearError>,
    },
    /// Further subdivision cannot retain a positive path/sliding increment.
    ResolutionLimit {
        /// Completed fraction in the discarded candidate path.
        completed_fraction: f64,
        /// Interval fraction that cannot be split or represented.
        interval_fraction: f64,
        /// Most recent step-dependent physical refusal, if any.
        last_refusal: Option<SourcePlasticWearError>,
    },
}
impl From<SourcePlasticWearError> for WearIntegrationError {
    fn from(error: SourcePlasticWearError) -> Self {
        Self::Step(error)
    }
}
impl core::fmt::Display for WearIntegrationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive source Archard wear: {self:?}")
    }
}
impl std::error::Error for WearIntegrationError {}

/// Accounts from accepted physical substeps only, never rejected trial histories.
#[derive(Clone, Debug, PartialEq)]
pub struct WearIntegrationReceipt {
    /// Exact local tolerances and work limits used for this complete path.
    pub controls: WearIntegrationControls,
    /// Actual accepted native wear steps, including both halves of each interval.
    pub accepted_substeps: usize,
    /// Intervals rejected by the local estimate or step-dependent physical limits.
    pub rejected_intervals: usize,
    /// All attempted native wear evaluations, including discarded trials.
    pub step_evaluations: usize,
    /// Largest accepted full-versus-half recession, plastic-depth, or approach difference [m].
    pub maximum_accepted_local_depth_error_m: f64,
    /// Largest accepted full-versus-half cell-pressure difference [Pa].
    pub maximum_accepted_local_pressure_error_pa: f64,
    /// Largest accepted component error divided by its local tolerance, <= 1.
    /// This local estimate is not an enclosure or a global path error guarantee.
    pub maximum_accepted_normalized_error: f64,
    /// Sum of accepted pressure-weighted local unsigned sliding measures [N m].
    pub pressure_weighted_sliding_j: f64,
    /// Sum of accepted original-body A/B removed-volume increments [m^3].
    pub removed_volume_m3: [f64; 2],
    /// Sum of accepted irreversible plastic work, distinct from wear volume [J].
    pub plastic_work_j: f64,
    /// Normal boundary work; zero at held approach [J].
    /// At held device position this is released spring energy, not actuator work.
    pub normal_work_j: f64,
    /// Sum of accepted elastic-storage changes, using the native stable difference [J].
    pub elastic_energy_change_j: f64,
    /// Signed parameter work, distinct from actuator work and heat [J].
    pub parameter_work_j: f64,
    /// Stored elastic energy before the complete path [J].
    pub initial_elastic_energy_j: f64,
    /// Stored elastic energy at its final accepted equilibrium [J].
    pub final_elastic_energy_j: f64,
    /// Floating residual `delta(U) + plastic_work - normal_work - parameter_work` [J].
    pub closure_residual_j: f64,
}

/// A complete sliding path with original source identity and retained histories.
#[derive(Clone, Debug, PartialEq)]
pub struct IntegratedWearAdvance {
    /// Complete original contact, coefficient, hardness, and applicability owner.
    pub source_identity: ContentHash,
    /// Candidate at the complete prescribed sliding distance, with both body inventories.
    pub candidate: SourcePlasticWearState,
    /// Accepted path accounts and the explicitly local numerical estimate.
    pub receipt: WearIntegrationReceipt,
}

impl SourcePlasticArchardWear {
    /// Integrate a fixed proportional local-sliding path at held rigid approach.
    /// Pressure is re-equilibrated after every candidate recession update.
    pub fn integrate_wear(
        &self,
        before: &SourcePlasticWearState,
        total_local_slip_m: &[f64],
        controls: WearIntegrationControls,
        cx: &Cx<'_>,
    ) -> Result<IntegratedWearAdvance, WearIntegrationError> {
        self.integrate_wear_with(before, total_local_slip_m, controls, cx, |state, slip| {
            self.advance_wear(state, slip, cx)
        })
    }

    /// Integrate the same local-sliding path while maintaining a prescribed load.
    ///
    /// The predecessor must already be equilibrated at this load. The native
    /// owner changes approach and separately accounts normal actuator work.
    /// Arbitrary load changes use `advance_normal_at_load` before starting wear.
    pub fn integrate_wear_at_load(
        &self,
        before: &SourcePlasticWearState,
        total_local_slip_m: &[f64],
        target_force_n: f64,
        load_controls: FiniteGapPlasticLoadControls,
        controls: WearIntegrationControls,
        cx: &Cx<'_>,
    ) -> Result<IntegratedWearAdvance, WearIntegrationError> {
        if !target_force_n.is_finite() || target_force_n < 0.0 {
            return Err(WearIntegrationError::Invalid(
                "finite nonnegative maintained normal force",
            ));
        }
        self.contact
            .admit_load_controls(load_controls)
            .map_err(SourcePlasticWearError::from)?;
        self.integrate_wear_with(before, total_local_slip_m, controls, cx, |state, slip| {
            self.advance_wear_at_load(state, slip, target_force_n, load_controls, cx)
                .map(|advance| advance.wear)
        })
    }

    // The boundary condition changes the physical leaf, not path integration.
    pub(super) fn integrate_wear_with(
        &self,
        before: &SourcePlasticWearState,
        total_local_slip_m: &[f64],
        controls: WearIntegrationControls,
        cx: &Cx<'_>,
        mut step: impl FnMut(
            &SourcePlasticWearState,
            &[f64],
        ) -> Result<SourcePlasticWearAdvance, SourcePlasticWearError>,
    ) -> Result<IntegratedWearAdvance, WearIntegrationError> {
        self.admit(before)?;
        checkpoint(cx)?;
        if !controls.absolute_depth_tolerance_m.is_finite()
            || controls.absolute_depth_tolerance_m <= 0.0
            || !controls.absolute_pressure_tolerance_pa.is_finite()
            || controls.absolute_pressure_tolerance_pa <= 0.0
            || !controls.relative_tolerance.is_finite()
            || controls.relative_tolerance < 0.0
            || controls.maximum_accepted_substeps == 0
            || controls.maximum_step_evaluations == 0
        {
            return Err(WearIntegrationError::Invalid(
                "positive finite absolute tolerances and work limits; nonnegative relative tolerance",
            ));
        }
        if total_local_slip_m.len() != before.contact.contact_state().pressure_pa().len()
            || total_local_slip_m
                .iter()
                .any(|value| !value.is_finite() || *value < 0.0)
        {
            return Err(WearIntegrationError::Invalid(
                "finite nonnegative total slip at every original cell",
            ));
        }
        let mut runner = Runner {
            controls,
            evaluations: 0,
            accepted: 0,
            last_refusal: None,
        };
        let mut accounts = Accounts::default();
        let mut rejected = 0;
        let mut maximum_error = LocalError::default();
        if total_local_slip_m.iter().all(|value| *value == 0.0) {
            // Still pass through the original leaf to check source/load ownership
            // and obtain its actual stored energy; zero sliding changes no history.
            let advance = runner
                .evaluate(before, total_local_slip_m, cx, &mut step)?
                .ok_or_else(|| {
                    runner.last_refusal.clone().map_or(
                        WearIntegrationError::Invalid("missing zero-path response"),
                        WearIntegrationError::Step,
                    )
                })?;
            accounts.record(&advance)?;
            runner.accepted = 1;
            return finish(
                self.identity(),
                advance.candidate,
                accounts,
                runner,
                rejected,
                maximum_error,
                cx,
            );
        }

        let mut current = before.clone();
        let mut lambda = 0.0_f64;
        let mut width = 1.0_f64;
        let mut full_slip = vec![0.0; total_local_slip_m.len()];
        let mut first_slip = full_slip.clone();
        let mut second_slip = full_slip.clone();
        while lambda < 1.0 {
            checkpoint(cx)?;
            if runner
                .accepted
                .checked_add(2)
                .is_none_or(|count| count > controls.maximum_accepted_substeps)
            {
                return Err(runner.budget_error());
            }
            let end = if width >= 1.0 - lambda {
                1.0
            } else {
                lambda + width
            };
            let interval = end - lambda;
            let middle = lambda + 0.5 * interval;
            let resolution = || WearIntegrationError::ResolutionLimit {
                completed_fraction: lambda,
                interval_fraction: interval,
                last_refusal: runner.last_refusal.clone(),
            };
            if !(middle > lambda && end > middle) {
                return Err(resolution());
            }
            for (index, distance) in total_local_slip_m.iter().copied().enumerate() {
                if index % 64 == 0 {
                    checkpoint(cx)?;
                }
                // Endpoint differences retain the prescribed total after many
                // unequal accepted intervals, instead of repeatedly adding D*h.
                let start_distance = distance * lambda;
                let middle_distance = distance * middle;
                let end_distance = if end == 1.0 { distance } else { distance * end };
                full_slip[index] = end_distance - start_distance;
                first_slip[index] = middle_distance - start_distance;
                second_slip[index] = end_distance - middle_distance;
                if distance > 0.0 && (first_slip[index] <= 0.0 || second_slip[index] <= 0.0) {
                    return Err(resolution());
                }
            }
            let Some(coarse) = runner.evaluate(&current, &full_slip, cx, &mut step)? else {
                rejected += 1;
                width = 0.5 * interval;
                continue;
            };
            let Some(first) = runner.evaluate(&current, &first_slip, cx, &mut step)? else {
                rejected += 1;
                width = 0.5 * interval;
                continue;
            };
            let Some(second) = runner.evaluate(&first.candidate, &second_slip, cx, &mut step)?
            else {
                rejected += 1;
                width = 0.5 * interval;
                continue;
            };
            let error = local_error(&current, &coarse.candidate, &second.candidate, controls, cx)?;
            if error.normalized > 1.0 {
                rejected += 1;
                runner.last_refusal = None;
                width = 0.5 * interval;
                continue;
            }
            // Commit neither trial here: only the returned whole-path candidate
            // is public. These accepted fine accounts remain local until success.
            accounts.record(&first)?;
            accounts.record(&second)?;
            runner.accepted += 2;
            runner.last_refusal = None;
            maximum_error.depth = maximum_error.depth.max(error.depth);
            maximum_error.pressure = maximum_error.pressure.max(error.pressure);
            maximum_error.normalized = maximum_error.normalized.max(error.normalized);
            current = second.candidate;
            lambda = end;
            width = if error.normalized <= 0.25 {
                2.0 * interval
            } else {
                interval
            };
        }
        finish(
            self.identity(),
            current,
            accounts,
            runner,
            rejected,
            maximum_error,
            cx,
        )
    }
}

struct Runner {
    controls: WearIntegrationControls,
    evaluations: usize,
    accepted: usize,
    last_refusal: Option<SourcePlasticWearError>,
}
impl Runner {
    fn budget_error(&self) -> WearIntegrationError {
        WearIntegrationError::BudgetExceeded {
            accepted_substeps: self.accepted,
            step_evaluations: self.evaluations,
            last_refusal: self.last_refusal.clone(),
        }
    }
    fn evaluate(
        &mut self,
        before: &SourcePlasticWearState,
        slip: &[f64],
        cx: &Cx<'_>,
        step: &mut impl FnMut(
            &SourcePlasticWearState,
            &[f64],
        ) -> Result<SourcePlasticWearAdvance, SourcePlasticWearError>,
    ) -> Result<Option<SourcePlasticWearAdvance>, WearIntegrationError> {
        checkpoint(cx)?;
        if self.evaluations == self.controls.maximum_step_evaluations {
            return Err(self.budget_error());
        }
        self.evaluations += 1;
        match step(before, slip) {
            Ok(advance) => Ok(Some(advance)),
            Err(error) if retryable(&error) => {
                checkpoint(cx)?;
                self.last_refusal = Some(error);
                Ok(None)
            }
            Err(error) => Err(error.into()),
        }
    }
}

// A large explicit geometry increment can breach a physical bound even when
// smaller, pressure-updated increments remain admitted. Source/shape, precision,
// cancellation, and hard resource refusals are not retried as smaller steps.
fn retryable(error: &SourcePlasticWearError) -> bool {
    match error {
        SourcePlasticWearError::Contact(SourcePlasticContactError::OutsideApplicability {
            ..
        })
        | SourcePlasticWearError::Load(SourcePlasticLoadError::Load(
            FiniteGapPlasticLoadError::TargetOutsideApproachBounds { .. },
        ))
        | SourcePlasticWearError::Drive(SourcePlasticDriveError::Drive(
            FiniteGapPlasticDriveError::TargetOutsideApproachBounds { .. },
        )) => true,
        SourcePlasticWearError::Contact(SourcePlasticContactError::Contact(native))
        | SourcePlasticWearError::Load(SourcePlasticLoadError::Contact(
            SourcePlasticContactError::Contact(native),
        ))
        | SourcePlasticWearError::Load(SourcePlasticLoadError::Load(
            FiniteGapPlasticLoadError::Contact(native),
        ))
        | SourcePlasticWearError::Drive(SourcePlasticDriveError::Contact(
            SourcePlasticContactError::Contact(native),
        ))
        | SourcePlasticWearError::Drive(SourcePlasticDriveError::Drive(
            FiniteGapPlasticDriveError::Contact(native),
        )) => matches!(
            native,
            FiniteGapPlasticError::Contact(
                FiniteGapContactError::ActiveSetDidNotConverge { .. }
                    | FiniteGapContactError::ContactReachedDomainBoundary { .. }
            )
        ),
        _ => false,
    }
}

#[derive(Clone, Copy, Default)]
struct LocalError {
    depth: f64,
    pressure: f64,
    normalized: f64,
}
fn local_error(
    before: &SourcePlasticWearState,
    coarse: &SourcePlasticWearState,
    fine: &SourcePlasticWearState,
    controls: WearIntegrationControls,
    cx: &Cx<'_>,
) -> Result<LocalError, WearIntegrationError> {
    checkpoint(cx)?;
    let mut error = LocalError::default();
    for side in 0..2 {
        for (index, ((old, full), half)) in before.recession_m[side]
            .iter()
            .zip(&coarse.recession_m[side])
            .zip(&fine.recession_m[side])
            .enumerate()
        {
            if index % 64 == 0 {
                checkpoint(cx)?;
            }
            depth_error(&mut error, *full - *old, *half - *old, controls)?;
        }
    }
    let initial = before.contact.contact_state();
    let full = coarse.contact.contact_state();
    let half = fine.contact.contact_state();
    depth_error(
        &mut error,
        full.approach_m() - initial.approach_m(),
        half.approach_m() - initial.approach_m(),
        controls,
    )?;
    for (index, ((old, coarse), fine)) in initial
        .plastic_depth_m()
        .iter()
        .zip(full.plastic_depth_m())
        .zip(half.plastic_depth_m())
        .enumerate()
    {
        if index % 64 == 0 {
            checkpoint(cx)?;
        }
        depth_error(&mut error, *coarse - *old, *fine - *old, controls)?;
    }
    for (index, (coarse, fine)) in full
        .pressure_pa()
        .iter()
        .zip(half.pressure_pa())
        .enumerate()
    {
        if index % 64 == 0 {
            checkpoint(cx)?;
        }
        let difference = (*coarse - *fine).abs();
        let tolerance = controls.absolute_pressure_tolerance_pa
            + controls.relative_tolerance * coarse.abs().max(fine.abs());
        error.pressure = error.pressure.max(difference);
        error.normalized = error.normalized.max(normalized(difference, tolerance)?);
    }
    checkpoint(cx)?;
    Ok(error)
}
fn depth_error(
    error: &mut LocalError,
    coarse: f64,
    fine: f64,
    controls: WearIntegrationControls,
) -> Result<(), WearIntegrationError> {
    let difference = (coarse - fine).abs();
    let tolerance = controls.absolute_depth_tolerance_m
        + controls.relative_tolerance * coarse.abs().max(fine.abs());
    error.depth = error.depth.max(difference);
    error.normalized = error.normalized.max(normalized(difference, tolerance)?);
    Ok(())
}
fn normalized(difference: f64, tolerance: f64) -> Result<f64, WearIntegrationError> {
    if !difference.is_finite() || !tolerance.is_finite() || tolerance <= 0.0 {
        return Err(WearIntegrationError::Invalid(
            "finite local error and tolerance",
        ));
    }
    // An overflowing ratio still unequivocally rejects the attempted interval.
    Ok(difference / tolerance)
}

#[derive(Clone, Copy, Default)]
struct Sum {
    value: f64,
    compensation: f64,
}
impl Sum {
    fn add(&mut self, value: f64) -> Result<(), WearIntegrationError> {
        if !value.is_finite() {
            return Err(WearIntegrationError::Invalid(
                "finite accepted work and volume accounts",
            ));
        }
        add_compensated(&mut self.value, &mut self.compensation, value);
        if !self.value.is_finite() {
            return Err(WearIntegrationError::Invalid(
                "finite integrated work and volume accounts",
            ));
        }
        Ok(())
    }
}
#[derive(Default)]
struct Accounts {
    measure: Sum,
    volumes: [Sum; 2],
    plastic: Sum,
    normal: Sum,
    elastic: Sum,
    parameter: Sum,
    initial_energy: Option<f64>,
    final_energy: f64,
}
impl Accounts {
    fn record(&mut self, advance: &SourcePlasticWearAdvance) -> Result<(), WearIntegrationError> {
        let receipt = &advance.receipt;
        self.initial_energy
            .get_or_insert(receipt.previous_elastic_energy_j);
        self.final_energy = receipt.elastic_energy_j;
        self.measure.add(advance.pressure_weighted_sliding_j)?;
        for side in 0..2 {
            self.volumes[side].add(advance.wear[side].increment_m3)?;
        }
        self.plastic.add(receipt.response.plastic_work_j)?;
        self.normal.add(receipt.normal_work_j)?;
        self.elastic.add(receipt.elastic_energy_change_j)?;
        self.parameter.add(receipt.parameter_work_j)?;
        Ok(())
    }
}
fn finish(
    source_identity: ContentHash,
    candidate: SourcePlasticWearState,
    accounts: Accounts,
    runner: Runner,
    rejected: usize,
    error: LocalError,
    cx: &Cx<'_>,
) -> Result<IntegratedWearAdvance, WearIntegrationError> {
    checkpoint(cx)?;
    let initial = accounts
        .initial_energy
        .ok_or(WearIntegrationError::Invalid(
            "missing accepted path account",
        ))?;
    let closure = accounts.elastic.value + accounts.plastic.value
        - accounts.normal.value
        - accounts.parameter.value;
    if !initial.is_finite() || !accounts.final_energy.is_finite() || !closure.is_finite() {
        return Err(WearIntegrationError::Invalid(
            "finite complete path energy account",
        ));
    }
    Ok(IntegratedWearAdvance {
        source_identity,
        candidate,
        receipt: WearIntegrationReceipt {
            controls: runner.controls,
            accepted_substeps: runner.accepted,
            rejected_intervals: rejected,
            step_evaluations: runner.evaluations,
            maximum_accepted_local_depth_error_m: error.depth,
            maximum_accepted_local_pressure_error_pa: error.pressure,
            maximum_accepted_normalized_error: error.normalized,
            pressure_weighted_sliding_j: accounts.measure.value,
            removed_volume_m3: accounts.volumes.map(|volume| volume.value),
            plastic_work_j: accounts.plastic.value,
            normal_work_j: accounts.normal.value,
            elastic_energy_change_j: accounts.elastic.value,
            parameter_work_j: accounts.parameter.value,
            initial_elastic_energy_j: initial,
            final_elastic_energy_j: accounts.final_energy,
            closure_residual_j: closure,
        },
    })
}
