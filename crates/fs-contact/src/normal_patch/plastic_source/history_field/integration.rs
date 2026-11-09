//! Adaptive pressure feedback for frozen accepted bulk-hardness fields.
//!
//! Each cell follows s_i(lambda) = total_slip_i * lambda. The original mapped
//! H and k remain fixed; native contact re-equilibrates pressure after removal.
//! Only actual fine histories are accepted. Local full/two-half-step estimates
//! do not certify global error or replace the caller's actual sliding path.

use fs_blake3::ContentHash;
use fs_exec::Cx;

pub use super::super::wear::integration::{WearIntegrationControls, WearIntegrationReceipt};
use super::{
    HistoryHardnessContactError, HistoryHardnessControl, HistoryHardnessWearAdvance,
    SourceHistoryHardnessWear, checkpoint,
};
use crate::normal_patch::{
    FiniteGapContactError, FiniteGapPlasticDriveError, FiniteGapPlasticError,
    FiniteGapPlasticLoadError,
};

/// A refusal returns no partial path or partially advanced physical history.
#[derive(Clone, Debug, PartialEq)]
pub enum HistoryHardnessIntegrationError {
    /// Invalid prescribed path, tolerance or finite arithmetic.
    Invalid(&'static str),
    /// Original source or native physical boundary refused an interval.
    Step(HistoryHardnessContactError),
    /// Finite native work allowance exhausted; all private candidates are discarded.
    BudgetExceeded {
        accepted_substeps: usize,
        step_evaluations: usize,
        last_refusal: Option<HistoryHardnessContactError>,
    },
    /// Further subdivision cannot retain positive local sliding increments.
    ResolutionLimit {
        completed_fraction: f64,
        interval_fraction: f64,
        last_refusal: Option<HistoryHardnessContactError>,
    },
}
impl From<HistoryHardnessContactError> for HistoryHardnessIntegrationError {
    fn from(error: HistoryHardnessContactError) -> Self {
        Self::Step(error)
    }
}
impl core::fmt::Display for HistoryHardnessIntegrationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive native-history hardness wear: {self:?}")
    }
}
impl std::error::Error for HistoryHardnessIntegrationError {}

/// Complete path successor with original mapped hardness and separate body inventories.
#[derive(Clone, Debug, PartialEq)]
pub struct HistoryHardnessIntegratedWearAdvance {
    pub source_identity: ContentHash,
    pub candidate: SourceHistoryHardnessWear,
    /// Accounts contain accepted fine steps only, never discarded trial work.
    pub receipt: WearIntegrationReceipt,
}

impl SourceHistoryHardnessWear {
    /// Integrate a prescribed proportional local-sliding path with pressure feedback.
    ///
    /// The accepted approach, maintained resultant or fixed loading spring must
    /// remain the same throughout this path. Frozen original point fields retain
    /// their source, cell map and accumulated bulk strain. Each body's removal
    /// uses its own mapped hardness; no averaging of hardness or plastic state
    /// is introduced. Refine tolerances to check path convergence.
    ///
    /// Input distances are actual unsigned slip, not released body motion or net
    /// travel. This adds no friction closure, bulk return map, thermal evolution,
    /// heat allocation or global energy theorem. Changing the slip pattern or
    /// accepted bulk field requires a separate physical path segment.
    pub fn integrate_wear(
        &self,
        total_local_slip_m: &[f64],
        control: HistoryHardnessControl,
        controls: WearIntegrationControls,
        cx: &Cx<'_>,
    ) -> Result<HistoryHardnessIntegratedWearAdvance, HistoryHardnessIntegrationError> {
        self.admit_control(control, true, cx)?;
        let total = total_local_slip_m;
        if !controls.absolute_depth_tolerance_m.is_finite()
            || controls.absolute_depth_tolerance_m <= 0.0
            || !controls.absolute_pressure_tolerance_pa.is_finite()
            || controls.absolute_pressure_tolerance_pa <= 0.0
            || !controls.relative_tolerance.is_finite()
            || controls.relative_tolerance < 0.0
            || controls.maximum_accepted_substeps == 0
            || controls.maximum_step_evaluations == 0
            || total.len() != self.contact.pressure_pa().len()
            || total.iter().any(|distance| !distance.is_finite() || *distance < 0.0)
        {
            return Err(HistoryHardnessIntegrationError::Invalid(
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
            // Admit the original boundary even when no material is removed.
            let next = run.evaluate(self, total, control, cx)?.ok_or_else(|| {
                run.last_refusal.clone().map_or(
                    HistoryHardnessIntegrationError::Invalid("missing zero-path response"),
                    HistoryHardnessIntegrationError::Step,
                )
            })?;
            accounts.record(&next)?;
            run.accepted = 1;
            return finish(self.source_identity(), next.candidate, accounts, run, 0, maximum, cx);
        }

        // Borrow the accepted predecessor until the first fine interval succeeds.
        let mut current: Option<Self> = None;
        let mut lambda = 0.0_f64;
        let mut width = 1.0_f64;
        let mut full = vec![0.0; total.len()];
        let mut left = full.clone();
        let mut right = full.clone();
        while lambda < 1.0 {
            checkpoint(cx)?;
            if run.accepted.checked_add(2)
                .is_none_or(|n| n > controls.maximum_accepted_substeps)
            {
                return Err(run.budget_error());
            }
            let end = if width >= 1.0 - lambda { 1.0 } else { lambda + width };
            let interval = end - lambda;
            let middle = lambda + 0.5 * interval;
            let resolution = || HistoryHardnessIntegrationError::ResolutionLimit {
                completed_fraction: lambda,
                interval_fraction: interval,
                last_refusal: run.last_refusal.clone(),
            };
            if !(lambda < middle && middle < end) {
                return Err(resolution());
            }
            for (i, distance) in total.iter().copied().enumerate() {
                if i % 64 == 0 { checkpoint(cx)?; }
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
            let before = current.as_ref().unwrap_or(self);
            let Some(coarse) = run.evaluate(before, &full, control, cx)? else {
                rejected += 1;
                width = 0.5 * interval;
                continue;
            };
            let Some(first) = run.evaluate(before, &left, control, cx)? else {
                rejected += 1;
                width = 0.5 * interval;
                continue;
            };
            let Some(second) = run.evaluate(&first.candidate, &right, control, cx)? else {
                rejected += 1;
                width = 0.5 * interval;
                continue;
            };
            let error = local_error(before, &coarse.candidate, &second.candidate, controls, cx)?;
            if error.normalized > 1.0 {
                rejected += 1;
                run.last_refusal = None;
                width = 0.5 * interval;
                continue;
            }
            accounts.record(&first)?;
            accounts.record(&second)?;
            run.accepted += 2;
            run.last_refusal = None;
            maximum.depth = maximum.depth.max(error.depth);
            maximum.pressure = maximum.pressure.max(error.pressure);
            maximum.normalized = maximum.normalized.max(error.normalized);
            current = Some(second.candidate);
            lambda = end;
            width = if error.normalized <= 0.25 { 2.0 * interval } else { interval };
        }
        finish(
            self.source_identity(),
            current.ok_or(HistoryHardnessIntegrationError::Invalid("missing accepted path"))?,
            accounts, run, rejected, maximum, cx,
        )
    }
}

struct Runner {
    controls: WearIntegrationControls,
    evaluations: usize,
    accepted: usize,
    last_refusal: Option<HistoryHardnessContactError>,
}
impl Runner {
    fn budget_error(&self) -> HistoryHardnessIntegrationError {
        HistoryHardnessIntegrationError::BudgetExceeded {
            accepted_substeps: self.accepted,
            step_evaluations: self.evaluations,
            last_refusal: self.last_refusal.clone(),
        }
    }
    fn evaluate(
        &mut self,
        before: &SourceHistoryHardnessWear,
        slip: &[f64],
        control: HistoryHardnessControl,
        cx: &Cx<'_>,
    ) -> Result<Option<HistoryHardnessWearAdvance>, HistoryHardnessIntegrationError> {
        checkpoint(cx)?;
        if self.evaluations == self.controls.maximum_step_evaluations {
            return Err(self.budget_error());
        }
        self.evaluations += 1;
        match before.advance_wear(slip, control, cx) {
            Ok(next) => Ok(Some(next)),
            Err(error) => {
                checkpoint(cx)?;
                // Only physical overshoots may improve after pressure-updated
                // subdivision. Numerical resolution and solver budgets never retry.
                if retryable(&error) {
                    self.last_refusal = Some(error);
                    Ok(None)
                } else {
                    Err(error.into())
                }
            }
        }
    }
}
fn retryable(error: &HistoryHardnessContactError) -> bool {
    match error {
        HistoryHardnessContactError::Contact(
            super::super::SourcePlasticContactError::OutsideApplicability { .. },
        )
        | HistoryHardnessContactError::Load(
            FiniteGapPlasticLoadError::TargetOutsideApproachBounds { .. },
        )
        | HistoryHardnessContactError::Drive(
            FiniteGapPlasticDriveError::TargetOutsideApproachBounds { .. },
        ) => true,
        HistoryHardnessContactError::Contact(
            super::super::SourcePlasticContactError::Contact(native),
        )
        | HistoryHardnessContactError::Load(FiniteGapPlasticLoadError::Contact(native))
        | HistoryHardnessContactError::Drive(FiniteGapPlasticDriveError::Contact(native)) => matches!(
            native,
            FiniteGapPlasticError::Contact(FiniteGapContactError::ContactReachedDomainBoundary { .. })
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
    before: &SourceHistoryHardnessWear,
    coarse: &SourceHistoryHardnessWear,
    fine: &SourceHistoryHardnessWear,
    controls: WearIntegrationControls,
    cx: &Cx<'_>,
) -> Result<LocalError, HistoryHardnessIntegrationError> {
    let mut error = LocalError::default();
    let depth_tolerance = controls.absolute_depth_tolerance_m;
    let relative = controls.relative_tolerance;
    for i in 0..before.contact.pressure_pa().len() {
        if i % 64 == 0 { checkpoint(cx)?; }
        for side in 0..2 {
            let old = before.recession_m[side][i];
            component(
                &mut error.depth, &mut error.normalized,
                coarse.recession_m[side][i] - old,
                fine.recession_m[side][i] - old,
                depth_tolerance, relative,
            )?;
        }
        let old = before.contact.plastic_depth_m()[i];
        component(
            &mut error.depth, &mut error.normalized,
            coarse.contact.plastic_depth_m()[i] - old,
            fine.contact.plastic_depth_m()[i] - old,
            depth_tolerance, relative,
        )?;
        component(
            &mut error.pressure, &mut error.normalized,
            coarse.contact.pressure_pa()[i], fine.contact.pressure_pa()[i],
            controls.absolute_pressure_tolerance_pa, relative,
        )?;
    }
    component(
        &mut error.depth, &mut error.normalized,
        coarse.contact.approach_m() - before.contact.approach_m(),
        fine.contact.approach_m() - before.contact.approach_m(),
        depth_tolerance, relative,
    )?;
    Ok(error)
}
fn component(
    maximum: &mut f64,
    normalized: &mut f64,
    coarse: f64,
    fine: f64,
    absolute: f64,
    relative: f64,
) -> Result<(), HistoryHardnessIntegrationError> {
    let difference = (coarse - fine).abs();
    let tolerance = absolute + relative * coarse.abs().max(fine.abs());
    if !difference.is_finite() || !tolerance.is_finite() || tolerance <= 0.0 {
        return Err(HistoryHardnessIntegrationError::Invalid("finite local error and tolerance"));
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
    fn add(&mut self, value: f64) -> Result<(), HistoryHardnessIntegrationError> {
        if !value.is_finite() {
            return Err(HistoryHardnessIntegrationError::Invalid("finite accepted account"));
        }
        let corrected = value - self.correction;
        let next = self.value + corrected;
        self.correction = (next - self.value) - corrected;
        self.value = next;
        if !self.value.is_finite() {
            return Err(HistoryHardnessIntegrationError::Invalid("finite integrated account"));
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
    fn record(&mut self, advance: &HistoryHardnessWearAdvance)
        -> Result<(), HistoryHardnessIntegrationError>
    {
        self.measure.add(advance.pressure_weighted_sliding_j)?;
        for side in 0..2 { self.volumes[side].add(advance.wear[side].increment_m3)?; }
        let receipt = &advance.receipt;
        self.initial_energy.get_or_insert(receipt.previous_elastic_energy_j);
        self.final_energy = receipt.elastic_energy_j;
        self.plastic.add(receipt.response.plastic_work_j)?;
        self.normal.add(receipt.normal_work_j)?;
        self.elastic.add(receipt.elastic_energy_change_j)?;
        self.parameter.add(receipt.parameter_work_j)?;
        Ok(())
    }
}
fn finish(
    source_identity: ContentHash,
    candidate: SourceHistoryHardnessWear,
    accounts: Accounts,
    run: Runner,
    rejected: usize,
    maximum: LocalError,
    cx: &Cx<'_>,
) -> Result<HistoryHardnessIntegratedWearAdvance, HistoryHardnessIntegrationError> {
    checkpoint(cx)?;
    let initial = accounts.initial_energy.ok_or(
        HistoryHardnessIntegrationError::Invalid("missing accepted path account"),
    )?;
    let closure = accounts.elastic.value + accounts.plastic.value
        - accounts.normal.value - accounts.parameter.value;
    if !initial.is_finite() || !accounts.final_energy.is_finite() || !closure.is_finite() {
        return Err(HistoryHardnessIntegrationError::Invalid("finite complete path energy account"));
    }
    Ok(HistoryHardnessIntegratedWearAdvance {
        source_identity,
        candidate,
        receipt: WearIntegrationReceipt {
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
    })
}
