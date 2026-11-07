//! Adaptive pressure and exposed-hardness feedback along actual local slip.
//!
//! Each original cell follows s_i(lambda) = total_slip_i * lambda. Full versus
//! two-half-step comparisons control local error, accepting only the actual
//! fine histories. Original depth integrals and normal boundaries remain with
//! their existing physical owners; irreversible state is never extrapolated.
//! Error estimates are local, not global enclosures or a slip/friction law.

use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_material::hardness::depth::DepthHardnessError;

pub use super::super::wear::integration::WearIntegrationControls;
use super::super::wear::integration::WearIntegrationReceipt;
use super::{
    DepthHardnessWearAdvance, DepthHardnessWearError, SourceDepthHardnessWear, add, checkpoint,
};
use crate::normal_patch::plastic_source::SourcePlasticContactError;
use crate::normal_patch::{
    FiniteGapContactError, FiniteGapPlasticDrive, FiniteGapPlasticDriveError,
    FiniteGapPlasticError, FiniteGapPlasticLoadControls, FiniteGapPlasticLoadError,
};

/// Refusal discards the complete trial path, including privately accepted steps.
#[derive(Clone, Debug, PartialEq)]
pub enum DepthHardnessIntegrationError {
    /// Invalid path, tolerance or finite arithmetic.
    Invalid(&'static str),
    /// Original material, contact or boundary refused the physical step.
    Step(DepthHardnessWearError),
    /// The complete path needs more work than the caller admitted.
    BudgetExceeded {
        /// Fine substeps in the discarded path.
        accepted_substeps: usize,
        /// All attempted native wear evaluations, including rejected trials.
        step_evaluations: usize,
        /// Latest step-dependent physical refusal, if any.
        last_refusal: Option<DepthHardnessWearError>,
    },
    /// Further subdivision cannot retain a positive path or local slip increment.
    ResolutionLimit {
        /// Completed fraction of the discarded candidate path.
        completed_fraction: f64,
        /// Interval that cannot be divided or represented.
        interval_fraction: f64,
        /// Latest step-dependent physical refusal, if any.
        last_refusal: Option<DepthHardnessWearError>,
    },
}
impl From<DepthHardnessWearError> for DepthHardnessIntegrationError {
    fn from(error: DepthHardnessWearError) -> Self {
        Self::Step(error)
    }
}
impl core::fmt::Display for DepthHardnessIntegrationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive exposed-depth hardness wear: {self:?}")
    }
}
impl std::error::Error for DepthHardnessIntegrationError {}

/// Accepted fine-path accounts and the explicitly local error estimate.
#[derive(Clone, Debug, PartialEq)]
pub struct DepthHardnessIntegrationReceipt {
    /// Existing path controls, accepted work/volume and depth/pressure errors.
    /// Both native continuations contribute once to each accepted wear step.
    pub path: WearIntegrationReceipt,
    /// Largest accepted error in either body's exposed H or the native cap [Pa].
    /// Uses the pressure tolerance; its ratio also enters the normalized error.
    pub maximum_accepted_local_hardness_error_pa: f64,
}

/// Complete successor at the entire prescribed proportional sliding path.
#[derive(Clone, Debug, PartialEq)]
pub struct DepthHardnessIntegratedWearAdvance {
    /// Original ordered material profiles, interface, grid and physical bounds.
    pub source_identity: ContentHash,
    /// Complete original-source contact, exposed depths and body inventories.
    pub candidate: SourceDepthHardnessWear,
    /// Accepted physical accounts, without any rejected trial work or volume.
    pub receipt: DepthHardnessIntegrationReceipt,
}

impl SourceDepthHardnessWear {
    /// Adapt a proportional local-sliding path at held approach.
    ///
    /// Absolute pressure tolerance also controls both original body hardnesses
    /// and native caps, even where pressure stays below a cap. Depth relative
    /// scales use this interval's increments, not earlier accumulated history.
    pub fn integrate_wear(
        &self,
        total_local_slip_m: &[f64],
        controls: WearIntegrationControls,
        cx: &Cx<'_>,
    ) -> Result<DepthHardnessIntegratedWearAdvance, DepthHardnessIntegrationError> {
        self.integrate_with(total_local_slip_m, controls, cx, |state, slip| {
            state.advance_wear_integrated(slip, cx)
        })
    }

    /// Adapt actual local slip while retaining the predecessor's applied force.
    /// Every physical leaf controls both removal and cap exposure at that force.
    pub fn integrate_wear_at_load(
        &self,
        total_local_slip_m: &[f64],
        target_force_n: f64,
        load_controls: FiniteGapPlasticLoadControls,
        controls: WearIntegrationControls,
        cx: &Cx<'_>,
    ) -> Result<DepthHardnessIntegratedWearAdvance, DepthHardnessIntegrationError> {
        checkpoint(cx)?;
        if !target_force_n.is_finite() || target_force_n < 0.0 {
            return Err(DepthHardnessIntegrationError::Invalid(
                "finite nonnegative maintained force",
            ));
        }
        super::within(
            "load-control maximum approach [m]",
            load_controls.maximum_approach_m,
            self.sources.envelope.maximum_approach_m,
        )
        .map_err(DepthHardnessWearError::from)?;
        self.integrate_with(total_local_slip_m, controls, cx, |state, slip| {
            state
                .advance_wear_integrated_at_load(slip, target_force_n, load_controls, cx)
                .map(|advance| advance.wear)
        })
    }

    /// Adapt actual local slip at fixed finite spring stiffness and position.
    /// Normal work is released spring storage, separate from parameter work.
    pub fn integrate_wear_with_drive(
        &self,
        total_local_slip_m: &[f64],
        drive: FiniteGapPlasticDrive,
        load_controls: FiniteGapPlasticLoadControls,
        controls: WearIntegrationControls,
        cx: &Cx<'_>,
    ) -> Result<DepthHardnessIntegratedWearAdvance, DepthHardnessIntegrationError> {
        checkpoint(cx)?;
        if !drive.stiffness_n_per_m.is_finite()
            || drive.stiffness_n_per_m <= 0.0
            || !drive.position_m.is_finite()
        {
            return Err(DepthHardnessIntegrationError::Invalid(
                "finite positive spring stiffness and finite position",
            ));
        }
        super::within(
            "load-control maximum approach [m]",
            load_controls.maximum_approach_m,
            self.sources.envelope.maximum_approach_m,
        )
        .map_err(DepthHardnessWearError::from)?;
        self.integrate_with(total_local_slip_m, controls, cx, |state, slip| {
            state
                .advance_wear_integrated_with_drive(slip, drive, load_controls, cx)
                .map(|advance| advance.wear)
        })
    }

    fn integrate_with(
        &self,
        total: &[f64],
        controls: WearIntegrationControls,
        cx: &Cx<'_>,
        mut step: impl FnMut(&Self, &[f64]) -> Result<DepthHardnessWearAdvance, DepthHardnessWearError>,
    ) -> Result<DepthHardnessIntegratedWearAdvance, DepthHardnessIntegrationError> {
        checkpoint(cx)?;
        if !controls.absolute_depth_tolerance_m.is_finite()
            || controls.absolute_depth_tolerance_m <= 0.0
            || !controls.absolute_pressure_tolerance_pa.is_finite()
            || controls.absolute_pressure_tolerance_pa <= 0.0
            || !controls.relative_tolerance.is_finite()
            || controls.relative_tolerance < 0.0
            || controls.maximum_accepted_substeps == 0
            || controls.maximum_step_evaluations == 0
            || total.len() != self.contact.pressure_pa().len()
            || total
                .iter()
                .any(|distance| !distance.is_finite() || *distance < 0.0)
        {
            return Err(DepthHardnessIntegrationError::Invalid(
                "finite nonnegative original-grid slip; positive absolute tolerances and budgets; nonnegative relative tolerance",
            ));
        }
        let mut run = Runner {
            controls,
            evaluations: 0,
            accepted: 0,
            last_refusal: None,
        };
        let mut accounts = Accounts::default();
        let mut rejected = 0;
        let mut maximum = LocalError::default();
        if total.iter().all(|distance| *distance == 0.0) {
            // Zero slip still admits the original requested normal boundary.
            let next = run.evaluate(self, total, cx, &mut step)?.ok_or_else(|| {
                run.last_refusal.clone().map_or(
                    DepthHardnessIntegrationError::Invalid("missing zero-path response"),
                    DepthHardnessIntegrationError::Step,
                )
            })?;
            accounts.record(self, &next, total, cx)?;
            run.accepted = 1;
            return finish(
                self.source_identity(),
                next.candidate,
                accounts,
                run,
                rejected,
                maximum,
                cx,
            );
        }
        let mut current = self.clone();
        let mut lambda = 0.0_f64;
        let mut width = 1.0_f64;
        let mut full = vec![0.0; total.len()];
        let mut left = full.clone();
        let mut right = full.clone();
        while lambda < 1.0 {
            checkpoint(cx)?;
            if run
                .accepted
                .checked_add(2)
                .is_none_or(|n| n > controls.maximum_accepted_substeps)
            {
                return Err(run.budget_error());
            }
            let end = if width >= 1.0 - lambda {
                1.0
            } else {
                lambda + width
            };
            let interval = end - lambda;
            let middle = lambda + 0.5 * interval;
            let resolution = || DepthHardnessIntegrationError::ResolutionLimit {
                completed_fraction: lambda,
                interval_fraction: interval,
                last_refusal: run.last_refusal.clone(),
            };
            if !(lambda < middle && middle < end) {
                return Err(resolution());
            }
            for (i, distance) in total.iter().copied().enumerate() {
                if i % 64 == 0 {
                    checkpoint(cx)?;
                }
                let a = distance * lambda;
                let b = distance * middle;
                let c = if end == 1.0 { distance } else { distance * end };
                full[i] = c - a;
                left[i] = b - a;
                right[i] = c - b;
                if distance > 0.0 && (left[i] <= 0.0 || right[i] <= 0.0) {
                    return Err(resolution());
                }
            }
            let Some(coarse) = run.evaluate(&current, &full, cx, &mut step)? else {
                rejected += 1;
                width = 0.5 * interval;
                continue;
            };
            let Some(first) = run.evaluate(&current, &left, cx, &mut step)? else {
                rejected += 1;
                width = 0.5 * interval;
                continue;
            };
            let Some(second) = run.evaluate(&first.candidate, &right, cx, &mut step)? else {
                rejected += 1;
                width = 0.5 * interval;
                continue;
            };
            let error = local_error(&current, &coarse.candidate, &second.candidate, controls, cx)?;
            if error.normalized > 1.0 {
                rejected += 1;
                run.last_refusal = None;
                width = 0.5 * interval;
                continue;
            }
            accounts.record(&current, &first, &left, cx)?;
            accounts.record(&first.candidate, &second, &right, cx)?;
            run.accepted += 2;
            run.last_refusal = None;
            maximum.depth = maximum.depth.max(error.depth);
            maximum.pressure = maximum.pressure.max(error.pressure);
            maximum.hardness = maximum.hardness.max(error.hardness);
            maximum.normalized = maximum.normalized.max(error.normalized);
            current = second.candidate;
            lambda = end;
            width = if error.normalized <= 0.25 {
                2.0 * interval
            } else {
                interval
            };
        }
        finish(
            self.source_identity(),
            current,
            accounts,
            run,
            rejected,
            maximum,
            cx,
        )
    }
}

struct Runner {
    controls: WearIntegrationControls,
    evaluations: usize,
    accepted: usize,
    last_refusal: Option<DepthHardnessWearError>,
}
impl Runner {
    fn budget_error(&self) -> DepthHardnessIntegrationError {
        DepthHardnessIntegrationError::BudgetExceeded {
            accepted_substeps: self.accepted,
            step_evaluations: self.evaluations,
            last_refusal: self.last_refusal.clone(),
        }
    }
    fn evaluate(
        &mut self,
        before: &SourceDepthHardnessWear,
        slip: &[f64],
        cx: &Cx<'_>,
        step: &mut impl FnMut(
            &SourceDepthHardnessWear,
            &[f64],
        ) -> Result<DepthHardnessWearAdvance, DepthHardnessWearError>,
    ) -> Result<Option<DepthHardnessWearAdvance>, DepthHardnessIntegrationError> {
        checkpoint(cx)?;
        if self.evaluations == self.controls.maximum_step_evaluations {
            return Err(self.budget_error());
        }
        self.evaluations += 1;
        match step(before, slip) {
            Ok(next) => Ok(Some(next)),
            Err(error) => {
                checkpoint(cx)?;
                if retryable(&error, before, slip, cx)? {
                    self.last_refusal = Some(error);
                    Ok(None)
                } else {
                    Err(error.into())
                }
            }
        }
    }
}

// Retry only physical overshoot that pressure-updated smaller steps can cure.
// Never hide cancellation, numerical underflow or a native solver work limit.
fn retryable(
    error: &DepthHardnessWearError,
    before: &SourceDepthHardnessWear,
    slip: &[f64],
    cx: &Cx<'_>,
) -> Result<bool, DepthHardnessIntegrationError> {
    Ok(match error {
        DepthHardnessWearError::Contact(SourcePlasticContactError::OutsideApplicability {
            ..
        })
        | DepthHardnessWearError::Load(FiniteGapPlasticLoadError::TargetOutsideApproachBounds {
            ..
        })
        | DepthHardnessWearError::Drive(
            FiniteGapPlasticDriveError::TargetOutsideApproachBounds { .. },
        ) => true,
        DepthHardnessWearError::Contact(SourcePlasticContactError::Contact(native))
        | DepthHardnessWearError::Load(FiniteGapPlasticLoadError::Contact(native))
        | DepthHardnessWearError::Drive(FiniteGapPlasticDriveError::Contact(native)) => matches!(
            native,
            FiniteGapPlasticError::Contact(
                FiniteGapContactError::ContactReachedDomainBoundary { .. }
            )
        ),
        DepthHardnessWearError::Hardness(DepthHardnessError::Invalid(_)) => {
            // Material errors are not matched by their prose. Diagnose support
            // overshoot against the original integral only after such a refusal.
            for (i, (&pressure, &distance)) in
                before.contact.pressure_pa().iter().zip(slip).enumerate()
            {
                if i % 64 == 0 {
                    checkpoint(cx)?;
                }
                for side in 0..2 {
                    let profile = &before.sources.hardness[side];
                    let depth = before.recession_m[side][i];
                    let measure = (before.sources.k[side] * pressure) * distance;
                    // Preserve the leaf's cell/body evaluation order, so a
                    // different cell's overshoot cannot mask an earlier
                    // arithmetic refusal that smaller steps cannot resolve.
                    if profile.advance_by_integral(depth, measure).is_err() {
                        return Ok(depth < profile.maximum_depth_m()
                            && profile
                                .hardness_integral_between(depth, profile.maximum_depth_m())
                                .is_ok_and(|capacity| measure > capacity));
                    }
                }
            }
            false
        }
        _ => false,
    })
}

#[derive(Clone, Copy, Default)]
struct LocalError {
    depth: f64,
    pressure: f64,
    hardness: f64,
    normalized: f64,
}
fn local_error(
    before: &SourceDepthHardnessWear,
    coarse: &SourceDepthHardnessWear,
    fine: &SourceDepthHardnessWear,
    controls: WearIntegrationControls,
    cx: &Cx<'_>,
) -> Result<LocalError, DepthHardnessIntegrationError> {
    let mut error = LocalError::default();
    let depth_tolerance = controls.absolute_depth_tolerance_m;
    let pressure_tolerance = controls.absolute_pressure_tolerance_pa;
    let relative = controls.relative_tolerance;
    for side in 0..2 {
        let profile = &before.sources.hardness[side];
        for i in 0..before.recession_m[side].len() {
            if i % 64 == 0 {
                checkpoint(cx)?;
            }
            let old = before.recession_m[side][i];
            let full = coarse.recession_m[side][i];
            let half = fine.recession_m[side][i];
            component(
                &mut error.depth,
                &mut error.normalized,
                full - old,
                half - old,
                depth_tolerance,
                relative,
            )?;
            let full_h = profile
                .hardness_pa_at(full)
                .map_err(DepthHardnessWearError::from)?;
            let half_h = profile
                .hardness_pa_at(half)
                .map_err(DepthHardnessWearError::from)?;
            component(
                &mut error.hardness,
                &mut error.normalized,
                full_h,
                half_h,
                pressure_tolerance,
                relative,
            )?;
        }
    }
    component(
        &mut error.depth,
        &mut error.normalized,
        coarse.contact.approach_m() - before.contact.approach_m(),
        fine.contact.approach_m() - before.contact.approach_m(),
        depth_tolerance,
        relative,
    )?;
    for i in 0..before.contact.pressure_pa().len() {
        if i % 64 == 0 {
            checkpoint(cx)?;
        }
        let old = before.contact.plastic_depth_m()[i];
        component(
            &mut error.depth,
            &mut error.normalized,
            coarse.contact.plastic_depth_m()[i] - old,
            fine.contact.plastic_depth_m()[i] - old,
            depth_tolerance,
            relative,
        )?;
        component(
            &mut error.pressure,
            &mut error.normalized,
            coarse.contact.pressure_pa()[i],
            fine.contact.pressure_pa()[i],
            pressure_tolerance,
            relative,
        )?;
        component(
            &mut error.hardness,
            &mut error.normalized,
            coarse.model.config().hardness_pa[i],
            fine.model.config().hardness_pa[i],
            pressure_tolerance,
            relative,
        )?;
    }
    Ok(error)
}
fn component(
    maximum: &mut f64,
    normalized: &mut f64,
    coarse: f64,
    fine: f64,
    absolute: f64,
    relative: f64,
) -> Result<(), DepthHardnessIntegrationError> {
    let difference = (coarse - fine).abs();
    let tolerance = absolute + relative * coarse.abs().max(fine.abs());
    if !difference.is_finite() || !tolerance.is_finite() || tolerance <= 0.0 {
        return Err(DepthHardnessIntegrationError::Invalid(
            "finite local error and tolerance",
        ));
    }
    *maximum = maximum.max(difference);
    *normalized = normalized.max(difference / tolerance);
    Ok(())
}

#[derive(Clone, Copy, Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) -> Result<(), DepthHardnessIntegrationError> {
        if !value.is_finite() {
            return Err(DepthHardnessIntegrationError::Invalid(
                "finite accepted account",
            ));
        }
        add(&mut self.value, &mut self.correction, value);
        if !self.value.is_finite() {
            return Err(DepthHardnessIntegrationError::Invalid(
                "finite integrated account",
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
    fn record(
        &mut self,
        before: &SourceDepthHardnessWear,
        advance: &DepthHardnessWearAdvance,
        slip: &[f64],
        cx: &Cx<'_>,
    ) -> Result<(), DepthHardnessIntegrationError> {
        let grid = before.model.config().grid;
        let area = grid.cell_width_m * grid.cell_depth_m;
        for (i, (&pressure, &distance)) in before.contact.pressure_pa().iter().zip(slip).enumerate()
        {
            if i % 64 == 0 {
                checkpoint(cx)?;
            }
            self.measure.add((pressure * area) * distance)?;
        }
        for side in 0..2 {
            self.volumes[side].add(advance.wear[side].increment_m3)?;
        }
        for receipt in [&advance.recession, &advance.hardness] {
            self.initial_energy
                .get_or_insert(receipt.previous_elastic_energy_j);
            self.final_energy = receipt.elastic_energy_j;
            self.plastic.add(receipt.response.plastic_work_j)?;
            self.normal.add(receipt.normal_work_j)?;
            self.elastic.add(receipt.elastic_energy_change_j)?;
            self.parameter.add(receipt.parameter_work_j)?;
        }
        Ok(())
    }
}
fn finish(
    source_identity: ContentHash,
    candidate: SourceDepthHardnessWear,
    accounts: Accounts,
    run: Runner,
    rejected: usize,
    maximum: LocalError,
    cx: &Cx<'_>,
) -> Result<DepthHardnessIntegratedWearAdvance, DepthHardnessIntegrationError> {
    checkpoint(cx)?;
    let initial = accounts
        .initial_energy
        .ok_or(DepthHardnessIntegrationError::Invalid(
            "missing accepted path account",
        ))?;
    let closure = accounts.elastic.value + accounts.plastic.value
        - accounts.normal.value
        - accounts.parameter.value;
    if !initial.is_finite() || !accounts.final_energy.is_finite() || !closure.is_finite() {
        return Err(DepthHardnessIntegrationError::Invalid(
            "finite complete path energy account",
        ));
    }
    Ok(DepthHardnessIntegratedWearAdvance {
        source_identity,
        candidate,
        receipt: DepthHardnessIntegrationReceipt {
            maximum_accepted_local_hardness_error_pa: maximum.hardness,
            path: WearIntegrationReceipt {
                controls: run.controls,
                accepted_substeps: run.accepted,
                rejected_intervals: rejected,
                step_evaluations: run.evaluations,
                maximum_accepted_local_depth_error_m: maximum.depth,
                maximum_accepted_local_pressure_error_pa: maximum.pressure,
                maximum_accepted_normalized_error: maximum.normalized,
                pressure_weighted_sliding_j: accounts.measure.value,
                removed_volume_m3: accounts.volumes.map(|sum| sum.value),
                plastic_work_j: accounts.plastic.value,
                normal_work_j: accounts.normal.value,
                elastic_energy_change_j: accounts.elastic.value,
                parameter_work_j: accounts.parameter.value,
                initial_elastic_energy_j: initial,
                final_elastic_energy_j: accounts.final_energy,
                closure_residual_j: closure,
            },
        },
    })
}
