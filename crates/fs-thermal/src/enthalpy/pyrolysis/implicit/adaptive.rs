//! Time-resolved decomposition, conduction and finite thermal storage.
//!
//! One original backward-Euler interval is compared with two original half
//! intervals. Only the two fine states and their physical releases are accepted;
//! no mass, composition, heat or enthalpy is extrapolated. An accuracy refusal
//! halves the unpublished comparison. Original nonlinear, source-domain and
//! conservation refusals propagate immediately instead of changing the model.
//!
//! Specific storage, temperature, separate virgin/char inventories and outgoing
//! mass/enthalpy all constrain acceptance. This resolves latent storage even
//! while a heat-only body's temperature stays on a plateau. Estimates are local
//! time discrepancies, not certified global error or material validation.

use super::*;

/// Accuracy controls separate from original algebraic conservation tolerances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdaptiveControls {
    /// Initial full comparison span, capped by `PyrolysisControls::maximum_dt_s`.
    pub initial_dt_s: f64,
    /// Smallest retried full span; a final target remainder may be shorter.
    pub minimum_dt_s: f64,
    /// Whole-call comparisons, including discarded accuracy trials. Each uses
    /// exactly three original bounded implicit intervals.
    pub maximum_comparisons: usize,
    /// Consecutive halvings before an accuracy refusal; zero disables retries.
    pub maximum_reductions: usize,
    pub storage_absolute_tolerance_j_kg: f64,
    pub temperature_absolute_tolerance_k: f64,
    /// Applied separately to virgin, char and outgoing volatile mass.
    pub mass_absolute_tolerance_kg: f64,
    /// Per-cell outgoing energy, relative to the initial cell's storage datum.
    pub volatile_energy_absolute_tolerance_j: f64,
    /// Scales the actual change during a comparison, never a gross formation
    /// enthalpy or an unchanged initial inventory. Must lie in [0, 1).
    pub relative_tolerance: f64,
}
impl Default for AdaptiveControls {
    fn default() -> Self {
        Self {
            initial_dt_s: 0.02,
            minimum_dt_s: 1e-8,
            maximum_comparisons: 10_000,
            maximum_reductions: 32,
            storage_absolute_tolerance_j_kg: 1.0,
            temperature_absolute_tolerance_k: 1e-3,
            mass_absolute_tolerance_kg: 1e-8,
            volatile_energy_absolute_tolerance_j: 1e-4,
            relative_tolerance: 1e-3,
        }
    }
}

/// Maximum raw coarse/fine discrepancies and the corrected fine-error ratio.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct TimeError {
    pub maximum_storage_difference_j_kg: f64,
    pub maximum_temperature_difference_k: f64,
    pub maximum_inventory_difference_kg: f64,
    pub maximum_volatile_mass_difference_kg: f64,
    pub maximum_volatile_energy_difference_j: f64,
    pub normalized: f64,
}

/// All original implicit work, including coarse and accuracy-rejected solves.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct AdaptiveWork {
    pub implicit_intervals: usize,
    pub nonlinear_iterations: usize,
    pub heat_evaluations: usize,
    pub local_cell_solves: usize,
    pub temperature_evaluations: usize,
    pub linear_solves: usize,
    pub line_probes: usize,
}
impl AdaptiveWork {
    fn include(&mut self, part: &ImplicitReport) -> Result<(), DepositionError> {
        self.implicit_intervals =
            add_count(self.implicit_intervals, part.pyrolysis.thermal.substeps)?;
        self.nonlinear_iterations =
            add_count(self.nonlinear_iterations, part.nonlinear_iterations)?;
        self.heat_evaluations = add_count(self.heat_evaluations, part.heat_evaluations)?;
        self.local_cell_solves = add_count(self.local_cell_solves, part.local_cell_solves)?;
        self.temperature_evaluations = add_count(
            self.temperature_evaluations,
            part.pyrolysis.temperature_evaluations,
        )?;
        self.linear_solves = add_count(self.linear_solves, part.linear_solves)?;
        self.line_probes = add_count(self.line_probes, part.line_probes)?;
        Ok(())
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct AdaptiveReport {
    /// Original physical accounts and numerical work of accepted fine intervals.
    /// Coarse/rejected heat, mass and volatile enthalpy never enter this receipt.
    pub physical: ImplicitReport,
    /// Complete numerical effort; `physical` is the accepted subset of this work.
    pub work: AdaptiveWork,
    pub comparisons: usize,
    pub accepted_comparisons: usize,
    pub accuracy_rejections: usize,
    pub maximum_accepted_normalized_error: f64,
    pub last_error: Option<TimeError>,
    /// Suggested next full comparison span; no hidden controller history exists.
    pub suggested_dt_s: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum AdaptiveError {
    Implicit(ImplicitError),
    Invalid(&'static str),
    ComparisonBudget,
    TimeAccuracy {
        time_s: f64,
        span_s: f64,
        error: TimeError,
    },
}
impl From<ImplicitError> for AdaptiveError {
    fn from(error: ImplicitError) -> Self {
        Self::Implicit(error)
    }
}
impl From<DepositionError> for AdaptiveError {
    fn from(error: DepositionError) -> Self {
        Self::Implicit(error.into())
    }
}
impl core::fmt::Display for AdaptiveError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive pyrolysis: {self:?}")
    }
}
impl std::error::Error for AdaptiveError {}

impl<'a> PyrolysisSimulation<'a> {
    /// Resolve the original fixed-geometry implicit model's local time error.
    ///
    /// Each accepted comparison publishes two unmodified native half steps.
    /// The original whole-call fine-substep and conservation ceilings persist
    /// across comparisons, while `maximum_comparisons` bounds ALL attempted
    /// native work. Geometry, contacts, materials and authored boundaries remain
    /// fixed throughout this call. This is not a moving-slab adaptive advance.
    ///
    /// Accuracy retries use only dyadic halving. Accepted spans double only
    /// below one quarter of the allowed first-order local error. A material,
    /// nonlinear, arithmetic or cancellation refusal never triggers a retry.
    /// The original physical state is unchanged on every error; only numerical
    /// workspace is reusable scratch. No artificial residual mass is retained.
    pub fn advance_adaptive_implicit_to(
        &self,
        target_s: f64,
        controls: PyrolysisControls,
        solve: ImplicitSolve,
        accuracy: AdaptiveControls,
        workspace: &mut ImplicitWorkspace,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, AdaptiveReport), AdaptiveError> {
        checkpoint()?;
        if !target_s.is_finite()
            || target_s < self.time_s()
            || !positive(accuracy.initial_dt_s)
            || !positive(accuracy.minimum_dt_s)
            || accuracy.minimum_dt_s > accuracy.initial_dt_s
            || accuracy.minimum_dt_s > controls.maximum_dt_s
            || accuracy.maximum_comparisons == 0
            || !positive(accuracy.storage_absolute_tolerance_j_kg)
            || !positive(accuracy.temperature_absolute_tolerance_k)
            || !positive(accuracy.mass_absolute_tolerance_kg)
            || !positive(accuracy.volatile_energy_absolute_tolerance_j)
            || !accuracy.relative_tolerance.is_finite()
            || !(0.0..1.0).contains(&accuracy.relative_tolerance)
        {
            return Err(AdaptiveError::Invalid(
                "adaptive pyrolysis controls or target",
            ));
        }
        // Original admission validates every algebraic control and the exact
        // workspace axis even for an idle request; no physical step is taken.
        let (mut next, physical) =
            self.advance_implicit_to(self.time_s(), controls, solve, workspace, &mut checkpoint)?;
        let mut report = AdaptiveReport {
            physical,
            work: AdaptiveWork::default(),
            comparisons: 0,
            accepted_comparisons: 0,
            accuracy_rejections: 0,
            maximum_accepted_normalized_error: 0.0,
            last_error: None,
            suggested_dt_s: accuracy.initial_dt_s.min(controls.maximum_dt_s),
        };
        let mut account = Account::default();
        let mut releases: BTreeMap<CellId, [Sum; 3]> = self
            .cells()
            .keys()
            .map(|&id| (id, Default::default()))
            .collect();
        let mut reductions = 0usize;
        while next.time_s() < target_s {
            checkpoint()?;
            if self
                .thermal
                .budget
                .max_substeps
                .saturating_sub(report.physical.pyrolysis.thermal.substeps)
                < 2
            {
                return Err(
                    DepositionError::Budget("adaptive pyrolysis fine half intervals").into(),
                );
            }
            if report.comparisons >= accuracy.maximum_comparisons {
                return Err(AdaptiveError::ComparisonBudget);
            }
            let start = next.time_s();
            let remaining = finite(target_s - start)?;
            let end = if report.suggested_dt_s >= remaining {
                target_s
            } else {
                finite(start + report.suggested_dt_s)?
            };
            let middle = f64::midpoint(start, end);
            if !(start < middle && middle < end && end <= target_s) {
                return Err(AdaptiveError::Invalid("no representable adaptive midpoint"));
            }
            let dt = finite(end - start)?;
            report.comparisons = add_count(report.comparisons, 1)?;
            let mut interval = |state: &Self, stop: f64| {
                let (state, result) = state.advance_implicit_to(
                    stop,
                    PyrolysisControls {
                        maximum_dt_s: stop - state.time_s(),
                        ..controls
                    },
                    solve,
                    workspace,
                    &mut checkpoint,
                )?;
                report.work.include(&result)?;
                Ok::<_, AdaptiveError>((state, result))
            };
            let (coarse, coarse_report) = interval(&next, end)?;
            let (half, first) = interval(&next, middle)?;
            let (fine, second) = interval(&half, end)?;
            // For BE the leading local error scales with dt^2. Account for
            // representationally unequal halves, without extrapolating state.
            let a = (middle - start) / dt;
            let b = (end - middle) / dt;
            let correction = finite((a * a + b * b) / (2.0 * a * b))?;
            let error = time_error(
                &next,
                &coarse,
                &half,
                &fine,
                [
                    &coarse_report.pyrolysis,
                    &first.pyrolysis,
                    &second.pyrolysis,
                ],
                accuracy,
                correction,
                &mut checkpoint,
            )?;
            if error.normalized <= 1.0 {
                for part in [&first, &second] {
                    include_accepted(&mut report.physical, part, &mut account, &mut releases)?;
                }
                check_account(&account, controls)?;
                report.accepted_comparisons = add_count(report.accepted_comparisons, 1)?;
                report.maximum_accepted_normalized_error = report
                    .maximum_accepted_normalized_error
                    .max(error.normalized);
                report.last_error = Some(error);
                report.suggested_dt_s = if error.normalized <= 0.25 {
                    dt.min(0.5 * controls.maximum_dt_s) * 2.0
                } else {
                    dt.min(controls.maximum_dt_s)
                }
                .max(accuracy.minimum_dt_s);
                next = fine;
                reductions = 0;
            } else {
                report.accuracy_rejections = add_count(report.accuracy_rejections, 1)?;
                let cause = AdaptiveError::TimeAccuracy {
                    time_s: start,
                    span_s: dt,
                    error,
                };
                if reductions >= accuracy.maximum_reductions || dt <= accuracy.minimum_dt_s {
                    return Err(cause);
                }
                let wanted = (0.5 * dt).max(accuracy.minimum_dt_s);
                let mut smaller = finite(start + wanted)?;
                if smaller - start < accuracy.minimum_dt_s {
                    smaller = smaller.next_up();
                }
                if !(start < smaller && smaller < end) {
                    return Err(cause);
                }
                report.suggested_dt_s = smaller - start;
                reductions = add_count(reductions, 1)?;
            }
        }
        finish_account(&mut report.physical, &account, &releases)?;
        checkpoint()?;
        Ok((next, report))
    }
}

#[allow(clippy::too_many_arguments)]
fn time_error(
    old: &PyrolysisSimulation<'_>,
    coarse: &PyrolysisSimulation<'_>,
    half: &PyrolysisSimulation<'_>,
    fine: &PyrolysisSimulation<'_>,
    audits: [&PyrolysisAudit; 3],
    controls: AdaptiveControls,
    correction: f64,
    checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
) -> Result<TimeError, DepositionError> {
    let mut error = TimeError::default();
    for (&id, cell) in old.cells() {
        checkpoint()?;
        let mut component = |before: f64, c: f64, f: f64, tolerance: f64| {
            let difference = finite(f - c)?.abs();
            let change = finite(c - before)?.abs().max(finite(f - before)?.abs());
            let allowance = finite(tolerance + controls.relative_tolerance * change)?;
            let ratio = correction * (difference / allowance);
            if ratio.is_nan() || ratio < 0.0 || !positive(allowance) {
                return Err(DepositionError::Arithmetic);
            }
            error.normalized = error.normalized.max(ratio);
            Ok(difference)
        };
        error.maximum_storage_difference_j_kg =
            error.maximum_storage_difference_j_kg.max(component(
                cell.specific_h_j_kg,
                coarse.cells()[&id].specific_h_j_kg,
                fine.cells()[&id].specific_h_j_kg,
                controls.storage_absolute_tolerance_j_kg,
            )?);
        error.maximum_temperature_difference_k =
            error.maximum_temperature_difference_k.max(component(
                old.temperature_k(id)?,
                coarse.temperature_k(id)?,
                fine.temperature_k(id)?,
                controls.temperature_absolute_tolerance_k,
            )?);
        if let Some(initial) = old.inventories.get(&id) {
            let c = coarse.inventories[&id];
            let f = fine.inventories[&id];
            for (before, c, f) in [
                (initial.virgin_mass_kg, c.virgin_mass_kg, f.virgin_mass_kg),
                (initial.char_mass_kg, c.char_mass_kg, f.char_mass_kg),
            ] {
                error.maximum_inventory_difference_kg = error.maximum_inventory_difference_kg.max(
                    component(before, c, f, controls.mass_absolute_tolerance_kg)?,
                );
            }
            let [c, first, second] = audits.map(|audit| audit.releases[&id]);
            error.maximum_volatile_mass_difference_kg =
                error.maximum_volatile_mass_difference_kg.max(component(
                    0.0,
                    c.mass_kg,
                    finite(first.mass_kg + second.mass_kg)?,
                    controls.mass_absolute_tolerance_kg,
                )?);
            let model = old.materials[&cell.material_id];
            // Requery the ORIGINAL immutable endpoint law and subtract the
            // common datum BEFORE multiplication. Rounded integrated absolute
            // release enthalpy cannot recover this difference afterward.
            // Precision remains limited by the original chart's returned h(T).
            let relative = |state: &PyrolysisSimulation<'_>, release: VolatileRelease| {
                if release.mass_kg == 0.0 {
                    return Ok(0.0);
                }
                let gas_h = finite(
                    model
                        .volatile
                        .specific_enthalpy_j_kg(state.inventories[&id].temperature_k)?,
                )?;
                finite(release.mass_kg * finite(gas_h - cell.specific_h_j_kg)?)
            };
            let mut fine_energy = Sum::default();
            fine_energy.add(relative(half, first)?)?;
            fine_energy.add(relative(fine, second)?)?;
            error.maximum_volatile_energy_difference_j =
                error.maximum_volatile_energy_difference_j.max(component(
                    0.0,
                    relative(coarse, c)?,
                    fine_energy.value(),
                    controls.volatile_energy_absolute_tolerance_j,
                )?);
        }
    }
    checkpoint()?;
    Ok(error)
}

fn include_accepted(
    total: &mut ImplicitReport,
    part: &ImplicitReport,
    account: &mut Account,
    releases: &mut BTreeMap<CellId, [Sum; 3]>,
) -> Result<(), DepositionError> {
    let a = &part.pyrolysis;
    account.energy.add(a.thermal.residual_j)?;
    account.absolute_energy.add(a.absolute_energy_residual_j)?;
    account.stored.add(a.thermal.stored_energy_change_j)?;
    account.mass.add(a.mass_residual_kg)?;
    account.external.add(a.thermal.external_heat_j)?;
    account.kinetic_mass.add(a.kinetic_mass_roundoff_kg)?;
    account.kinetic_energy.add(a.kinetic_energy_roundoff_j)?;
    account
        .heat_roundoff
        .add(part.absolute_heat_balance_roundoff_j)?;
    total.pyrolysis.thermal.substeps =
        add_count(total.pyrolysis.thermal.substeps, a.thermal.substeps)?;
    total.pyrolysis.temperature_evaluations = add_count(
        total.pyrolysis.temperature_evaluations,
        a.temperature_evaluations,
    )?;
    total.pyrolysis.maximum_cell_energy_residual_j = total
        .pyrolysis
        .maximum_cell_energy_residual_j
        .max(a.maximum_cell_energy_residual_j);
    total.nonlinear_iterations = add_count(total.nonlinear_iterations, part.nonlinear_iterations)?;
    total.heat_evaluations = add_count(total.heat_evaluations, part.heat_evaluations)?;
    total.local_cell_solves = add_count(total.local_cell_solves, part.local_cell_solves)?;
    total.linear_solves = add_count(total.linear_solves, part.linear_solves)?;
    total.line_probes = add_count(total.line_probes, part.line_probes)?;
    for (&id, part) in &a.releases {
        let out = releases.get_mut(&id).ok_or(DepositionError::Arithmetic)?;
        for (sum, value) in
            out.iter_mut()
                .zip([part.mass_kg, part.enthalpy_j, part.consumed_virgin_mass_kg])
        {
            sum.add(value)?;
        }
    }
    Ok(())
}

fn check_account(account: &Account, controls: PyrolysisControls) -> Result<(), DepositionError> {
    if account.absolute_energy.value() > controls.total_energy_tolerance_j {
        return Err(DepositionError::EnergyBalance(
            account.absolute_energy.value(),
        ));
    }
    if account.absolute_energy.value() + account.kinetic_energy.value()
        > controls.total_energy_tolerance_j
    {
        return Err(DepositionError::Budget(
            "adaptive pyrolysis accumulated reaction-energy roundoff",
        ));
    }
    if account.mass.value().abs() > controls.mass_tolerance_kg {
        return Err(DepositionError::InvalidInput(
            "adaptive pyrolysis solid/volatile mass balance",
        ));
    }
    if account.mass.value().abs() + account.kinetic_mass.value() > controls.mass_tolerance_kg {
        return Err(DepositionError::Budget(
            "adaptive pyrolysis accumulated kinetic inventory roundoff",
        ));
    }
    Ok(())
}

fn finish_account(
    total: &mut ImplicitReport,
    account: &Account,
    releases: &BTreeMap<CellId, [Sum; 3]>,
) -> Result<(), DepositionError> {
    let audit = &mut total.pyrolysis;
    audit.thermal.external_heat_j = account.external.value();
    audit.thermal.stored_energy_change_j = account.stored.value();
    audit.thermal.residual_j = account.energy.value();
    audit.absolute_energy_residual_j = account.absolute_energy.value();
    audit.mass_residual_kg = account.mass.value();
    audit.kinetic_mass_roundoff_kg = account.kinetic_mass.value();
    audit.kinetic_energy_roundoff_j = account.kinetic_energy.value();
    total.absolute_heat_balance_roundoff_j = account.heat_roundoff.value();
    let mut mass = Sum::default();
    let mut energy = Sum::default();
    for (&id, sums) in releases {
        let release = VolatileRelease {
            mass_kg: sums[0].value(),
            enthalpy_j: sums[1].value(),
            consumed_virgin_mass_kg: sums[2].value(),
        };
        mass.add(release.mass_kg)?;
        energy.add(release.enthalpy_j)?;
        *audit
            .releases
            .get_mut(&id)
            .ok_or(DepositionError::Arithmetic)? = release;
    }
    audit.volatile_mass_kg = mass.value();
    audit.volatile_enthalpy_j = energy.value();
    Ok(())
}

fn add_count(a: usize, b: usize) -> Result<usize, DepositionError> {
    a.checked_add(b)
        .ok_or(DepositionError::Budget("adaptive pyrolysis work counter"))
}
