//! Complete diffusion windows with componentwise pressure and amount accuracy.
//!
//! One original backward-Euler step is compared with two original half steps;
//! only the half-step trajectory advances storage and integrated face amounts.
//! The first-order local estimator and change-based scales match the existing
//! finite-mixture membrane driver. This owner also covers prescribed reservoirs
//! and finite molar pressure capacities without requiring a gas caloric model.
//! A source-bound adapter calls its own native step for every probe.
//!
//! Boundary pressures are held fixed for the requested window. Split calls at
//! prescribed boundary changes; no schedule or interpolated forcing is inferred.
//! Temporal estimates are not whole-window error bounds or spatial certificates.
//! Every native probe also passes a stricter algebraic amount audit: one percent
//! of its pressure/storage accuracy and a duration-weighted share of the window
//! balance budget. Original native acceptance controls are never relaxed.

use super::finite::{FiniteMembraneDiffusionState, FiniteMembraneDiffusionWorkspace};
use super::{
    MembraneDiffusionError, MembraneDiffusionSolve, MembraneDiffusionState, MembraneDiffusionStep,
    MembraneDiffusionWorkspace, Sum, nonnegative, positive,
};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MembraneDiffusionAdvanceControls {
    /// Original per-step algebraic amount acceptance. An additional residual
    /// audit keeps algebraic errors below the requested temporal/balance budgets.
    pub step: MembraneDiffusionSolve,
    pub initial_step_s: f64,
    /// Floor for retried comparison intervals. Accepted halves can be smaller;
    /// a shorter final remainder is attempted without relaxing accuracy.
    pub minimum_step_s: f64,
    /// A final one-ULP clock remainder can extend this ceiling by that roundoff.
    pub maximum_step_s: f64,
    /// Original implicit calls, including all rejected full/half probes.
    pub maximum_evaluations: usize,
    pub absolute_pressure_tolerance_pa: f64,
    /// Every sorbed cell, finite reservoir and integrated face is compared.
    pub absolute_amount_tolerance_mol: f64,
    /// Multiplies each coordinate's change, never its absolute pressure datum.
    pub relative_tolerance: f64,
    /// One fixed absolute allowance for every whole-window physical balance.
    pub amount_balance_tolerance_mol: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum MembraneDiffusionAdvanceError<E = MembraneDiffusionError> {
    Invalid(&'static str),
    Native(E),
    Accuracy {
        time_s: f64,
        step_s: f64,
        error_ratio: f64,
    },
    WorkLimit {
        time_s: f64,
        evaluations: usize,
    },
    TimeResolution {
        time_s: f64,
        target_s: f64,
    },
    AmountBalance {
        cell: Option<usize>,
        residual_mol: f64,
        tolerance_mol: f64,
    },
    ReservoirBalance {
        reservoir: usize,
        residual_mol: f64,
        tolerance_mol: f64,
    },
    Cancelled,
}
impl<E: core::fmt::Debug> core::fmt::Display for MembraneDiffusionAdvanceError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive membrane diffusion: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for MembraneDiffusionAdvanceError<E> {}
use MembraneDiffusionAdvanceError as Error;
type Result<T, E = MembraneDiffusionError> = core::result::Result<T, Error<E>>;

/// Separate full/half caches retain native factors for repeated interval sizes.
#[derive(Debug)]
pub struct MembraneDiffusionAdvanceWorkspace<W = MembraneDiffusionWorkspace> {
    pub(crate) full: W,
    pub(crate) half: W,
}

#[derive(Clone, Debug, PartialEq)]
pub struct FiniteMembraneAdvanceBalance {
    pub pressures_before_pa: [f64; 2],
    pub pressures_after_pa: [f64; 2],
    pub amount_changes_mol: [f64; 2],
    pub reservoir_residuals_mol: [f64; 2],
    /// Reservoir plus sorbed changes, computed from original-to-final states.
    pub closed_amount_residual_mol: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct MembraneDiffusionAdvanceReport {
    pub interval_s: [f64; 2],
    pub accepted_intervals: usize,
    pub accepted_substeps: usize,
    pub rejected_intervals: usize,
    pub evaluations: usize,
    pub minimum_accepted_substep_s: f64,
    pub maximum_accepted_substep_s: f64,
    pub maximum_accepted_error_ratio: f64,
    pub suggested_step_s: f64,
    /// Signed left-to-right amounts from accepted physical substeps only.
    pub face_amounts_mol: Vec<f64>,
    /// Sum of absolute transfers, retaining physical direction reversals.
    pub gross_face_amounts_mol: Vec<f64>,
    pub stored_before_mol: f64,
    pub stored_after_mol: f64,
    pub stored_change_mol: f64,
    pub amount_residual_mol: f64,
    pub maximum_cell_amount_residual_mol: f64,
    /// None means prescribed external pressures; it does not imply finite gas.
    pub finite_reservoirs: Option<FiniteMembraneAdvanceBalance>,
    pub last_substep: MembraneDiffusionStep,
}

#[derive(Clone, Debug, PartialEq)]
pub struct MembraneDiffusionAdvanceCandidate<S = MembraneDiffusionState> {
    state: S,
    report: MembraneDiffusionAdvanceReport,
}
impl<S> MembraneDiffusionAdvanceCandidate<S> {
    pub const fn state(&self) -> &S {
        &self.state
    }
    pub const fn report(&self) -> &MembraneDiffusionAdvanceReport {
        &self.report
    }
    pub fn into_parts(self) -> (S, MembraneDiffusionAdvanceReport) {
        (self.state, self.report)
    }
}

#[derive(Clone, Copy)]
pub(crate) struct Reservoirs {
    pub capacities: [f64; 2],
    pub pressures: [f64; 2],
}

// Only original numerical/source-bound membrane owners implement this adapter.
// Cell identities, coefficients, geometry and finite capacities remain fixed.
pub(crate) trait Native: Clone {
    type Workspace;
    type Error;
    fn membrane(&self) -> &MembraneDiffusionState;
    fn reservoirs(&self) -> Option<Reservoirs> {
        None
    }
    fn is_cancelled(error: &Self::Error) -> bool;
    fn native_step(
        &self,
        reservoirs: [f64; 2],
        duration: f64,
        solve: MembraneDiffusionSolve,
        work: &mut Self::Workspace,
        cancelled: &mut impl FnMut() -> bool,
    ) -> core::result::Result<(Self, MembraneDiffusionStep), Self::Error>;
}
impl Native for MembraneDiffusionState {
    type Workspace = MembraneDiffusionWorkspace;
    type Error = MembraneDiffusionError;
    fn membrane(&self) -> &MembraneDiffusionState {
        self
    }
    fn is_cancelled(error: &Self::Error) -> bool {
        *error == MembraneDiffusionError::Cancelled
    }
    fn native_step(
        &self,
        reservoirs: [f64; 2],
        duration: f64,
        solve: MembraneDiffusionSolve,
        work: &mut Self::Workspace,
        cancelled: &mut impl FnMut() -> bool,
    ) -> core::result::Result<(Self, MembraneDiffusionStep), Self::Error> {
        Ok(self
            .step_candidate(reservoirs, duration, solve, work, cancelled)?
            .into_parts())
    }
}
impl Native for FiniteMembraneDiffusionState {
    type Workspace = FiniteMembraneDiffusionWorkspace;
    type Error = MembraneDiffusionError;
    fn membrane(&self) -> &MembraneDiffusionState {
        self.membrane()
    }
    fn reservoirs(&self) -> Option<Reservoirs> {
        Some(Reservoirs {
            capacities: self.capacities().map(|c| c.value_si()),
            pressures: self.reservoir_pressures_pa(),
        })
    }
    fn is_cancelled(error: &Self::Error) -> bool {
        *error == MembraneDiffusionError::Cancelled
    }
    fn native_step(
        &self,
        _: [f64; 2],
        duration: f64,
        solve: MembraneDiffusionSolve,
        work: &mut Self::Workspace,
        cancelled: &mut impl FnMut() -> bool,
    ) -> core::result::Result<(Self, MembraneDiffusionStep), Self::Error> {
        let (state, report) = self
            .step_candidate(duration, solve, work, cancelled)?
            .into_parts();
        Ok((state, report.membrane))
    }
}

fn poll<E>(cancelled: &mut impl FnMut() -> bool) -> Result<(), E> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn finite<E>(value: f64) -> Result<f64, E> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(Error::Invalid(
            "unrepresentable adaptive membrane arithmetic",
        ))
    }
}
fn add<E>(sum: &mut Sum, value: f64) -> Result<(), E> {
    sum.add(value)
        .map_err(|_| Error::Invalid("unrepresentable integrated membrane amount"))
}
fn validate<E>(start: f64, duration: f64, c: MembraneDiffusionAdvanceControls) -> Result<(), E> {
    if !positive(duration)
        || !positive(c.initial_step_s)
        || !positive(c.minimum_step_s)
        || !positive(c.maximum_step_s)
        || c.minimum_step_s > c.initial_step_s
        || c.initial_step_s > c.maximum_step_s
        || !(3..=100_000).contains(&c.maximum_evaluations)
        || !positive(c.absolute_pressure_tolerance_pa)
        || !positive(c.absolute_amount_tolerance_mol)
        || !positive(c.amount_balance_tolerance_mol)
        || !c.relative_tolerance.is_finite()
        || !(0.0..1.0).contains(&c.relative_tolerance)
        || c.step.validate().is_err()
    {
        return Err(Error::Invalid(
            "positive forward duration, bounded work and explicit accuracy/balance controls required",
        ));
    }
    if !nonnegative(start) || !(start + duration).is_finite() || start + duration <= start {
        return Err(Error::TimeResolution {
            time_s: start,
            target_s: start + duration,
        });
    }
    Ok(())
}
fn check_clock<S: Native>(state: &S, target: f64) -> Result<(), S::Error> {
    let time = state.membrane().elapsed_s();
    if time.to_bits() != target.to_bits() {
        Err(Error::TimeResolution {
            time_s: time,
            target_s: target,
        })
    } else {
        Ok(())
    }
}
fn probe<S: Native>(
    state: &S,
    reservoirs: [f64; 2],
    span: f64,
    window: f64,
    controls: MembraneDiffusionAdvanceControls,
    work: &mut S::Workspace,
    calls: &mut usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(S, MembraneDiffusionStep), S::Error> {
    poll(cancelled)?;
    if *calls >= controls.maximum_evaluations {
        return Err(Error::WorkLimit {
            time_s: state.membrane().elapsed_s(),
            evaluations: *calls,
        });
    }
    *calls += 1;
    let result = state
        .native_step(reservoirs, span, controls.step, work, cancelled)
        .map_err(|e| {
            if S::is_cancelled(&e) {
                Error::Cancelled
            } else {
                Error::Native(e)
            }
        })?;
    algebraic_balance(state, &result.0, &result.1, window, controls, cancelled)?;
    Ok(result)
}
// Native solves are direct factorizations: their controls admit or refuse a
// measured residual; they do not tune iteration accuracy. Keep their original
// checks, including relative-only requests, and independently bound each amount
// equation before estimating temporal error. The amount allowance is shared
// across inventories so large solubility/capacity contrasts cannot hide a large
// state error behind a small per-row pressure error. Accepted durations sum to
// the window, so the balance allowance never grows with the number of steps.
fn algebraic_balance<S: Native>(
    old: &S,
    new: &S,
    report: &MembraneDiffusionStep,
    window: f64,
    controls: MembraneDiffusionAdvanceControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), S::Error> {
    let model = &old.membrane().model;
    let inventories = model.cells.len() + if old.reservoirs().is_some() { 2 } else { 0 };
    let global = (0.01 * controls.absolute_amount_tolerance_mol)
        .min((0.01 * controls.amount_balance_tolerance_mol) * (report.duration_s / window));
    let per_inventory = global / inventories as f64;
    if !positive(per_inventory) {
        return Err(Error::Invalid(
            "algebraic share of membrane accuracy/balance budget underflowed",
        ));
    }
    let allowance = |capacity: f64| -> Result<f64, S::Error> {
        let amount = per_inventory.min((0.01 * controls.absolute_pressure_tolerance_pa) * capacity);
        if positive(amount) {
            Ok(amount)
        } else {
            Err(Error::Invalid(
                "algebraic pressure/storage budget underflowed",
            ))
        }
    };
    for (i, cell) in model.cells.iter().enumerate() {
        if i % 128 == 0 {
            poll(cancelled)?;
        }
        let capacity = cell.capacity_mol_m2_pa * model.area_m2;
        let delta = finite(
            capacity
                * (new.membrane().partial_pressures_pa[i] - old.membrane().partial_pressures_pa[i]),
        )?;
        let incoming = face_amount(report, i, model.area_m2)?;
        let outgoing = face_amount(report, i + 1, model.area_m2)?;
        balance(
            Some(i),
            finite(delta - (incoming - outgoing))?,
            allowance(capacity)?,
        )?;
    }
    balance(None, report.amount_residual_mol, global)?;
    if let Some(before) = old.reservoirs() {
        let after = new.reservoirs().expect("fixed finite reservoirs");
        let outward = [
            face_amount(report, 0, model.area_m2)?,
            -face_amount(report, model.faces.len() - 1, model.area_m2)?,
        ];
        let mut total = Sum::default();
        add(&mut total, report.stored_change_mol)?;
        for side in 0..2 {
            let delta =
                finite(before.capacities[side] * (after.pressures[side] - before.pressures[side]))?;
            let residual = finite(delta + outward[side])?;
            let tolerance = allowance(before.capacities[side])?;
            if residual.abs() > tolerance {
                return Err(Error::ReservoirBalance {
                    reservoir: side,
                    residual_mol: residual,
                    tolerance_mol: tolerance,
                });
            }
            add(&mut total, delta)?;
        }
        balance(None, total.value, global)?;
    }
    Ok(())
}
fn ratio<E>(difference: f64, change: f64, absolute: f64, relative: f64) -> Result<f64, E> {
    let allowance = finite(absolute + relative * change)?;
    Ok(finite(difference)?.abs() / allowance)
}
fn coordinate_error<E>(
    old: f64,
    a: f64,
    b: f64,
    capacity: f64,
    c: MembraneDiffusionAdvanceControls,
) -> Result<f64, E> {
    let difference = a - b;
    let change = (a - old).abs().max((b - old).abs());
    let pressure = ratio(
        difference,
        change,
        c.absolute_pressure_tolerance_pa,
        c.relative_tolerance,
    )?;
    let amount = ratio(
        difference * capacity,
        change * capacity,
        c.absolute_amount_tolerance_mol,
        c.relative_tolerance,
    )?;
    Ok(pressure.max(amount))
}
fn face_amount<E>(report: &MembraneDiffusionStep, face: usize, area: f64) -> Result<f64, E> {
    finite((report.face_molar_fluxes_mol_m2_s[face] * area) * report.duration_s)
}
fn distance<S: Native>(
    old: &S,
    a: &S,
    b: &S,
    reports: [&MembraneDiffusionStep; 3],
    c: MembraneDiffusionAdvanceControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, S::Error> {
    let model = &old.membrane().model;
    let mut error = 0.0_f64;
    for i in 0..model.cells.len() {
        if i % 128 == 0 {
            poll(cancelled)?;
        }
        error = error.max(coordinate_error(
            old.membrane().partial_pressures_pa[i],
            a.membrane().partial_pressures_pa[i],
            b.membrane().partial_pressures_pa[i],
            model.cells[i].capacity_mol_m2_pa * model.area_m2,
            c,
        )?);
    }
    if let Some(initial) = old.reservoirs() {
        let (a, b) = (
            a.reservoirs().expect("fixed finite reservoirs"),
            b.reservoirs().expect("fixed finite reservoirs"),
        );
        for side in 0..2 {
            error = error.max(coordinate_error(
                initial.pressures[side],
                a.pressures[side],
                b.pressures[side],
                initial.capacities[side],
                c,
            )?);
        }
    }
    for face in 0..model.faces.len() {
        if face % 128 == 0 {
            poll(cancelled)?;
        }
        let coarse = face_amount(reports[0], face, model.area_m2)?;
        let first = face_amount(reports[1], face, model.area_m2)?;
        let second = face_amount(reports[2], face, model.area_m2)?;
        let fine = finite(first + second)?;
        error = error.max(ratio(
            coarse - fine,
            coarse.abs().max(first.abs() + second.abs()),
            c.absolute_amount_tolerance_mol,
            c.relative_tolerance,
        )?);
    }
    Ok(error)
}

pub(crate) fn advance<S: Native>(
    initial: &S,
    reservoirs: [f64; 2],
    duration: f64,
    controls: MembraneDiffusionAdvanceControls,
    work: &mut MembraneDiffusionAdvanceWorkspace<S::Workspace>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<MembraneDiffusionAdvanceCandidate<S>, S::Error> {
    poll(&mut cancelled)?;
    let start = initial.membrane().elapsed_s();
    validate(start, duration, controls)?;
    if !reservoirs.into_iter().all(nonnegative) {
        return Err(Error::Invalid(
            "nonnegative prescribed reservoir pressures required",
        ));
    }
    let target = start + duration;
    let window = target - start;
    let model = &initial.membrane().model;
    let mut signed: Vec<Sum> = (0..model.faces.len()).map(|_| Sum::default()).collect();
    let mut gross: Vec<Sum> = (0..model.faces.len()).map(|_| Sum::default()).collect();
    let mut state = initial.clone();
    let (mut accepted, mut rejected, mut evaluations) = (0, 0, 0);
    let (mut minimum, mut maximum, mut maximum_error) = (f64::INFINITY, 0.0_f64, 0.0_f64);
    let mut step = controls.initial_step_s;
    let mut last = None;
    while state.membrane().elapsed_s() < target {
        poll(&mut cancelled)?;
        let time = state.membrane().elapsed_s();
        let planned = if step >= target - time {
            target
        } else {
            time + step
        };
        let end = if planned < target && planned.next_up() >= target {
            target
        } else {
            planned
        };
        let span = end - time;
        let middle = time + 0.5 * span;
        if !(time < middle && middle < end) {
            return Err(Error::TimeResolution {
                time_s: time,
                target_s: end,
            });
        }
        let (h1, h2) = (middle - time, end - middle);
        let (coarse, full) = probe(
            &state,
            reservoirs,
            span,
            window,
            controls,
            &mut work.full,
            &mut evaluations,
            &mut cancelled,
        )?;
        let (half, first) = probe(
            &state,
            reservoirs,
            h1,
            window,
            controls,
            &mut work.half,
            &mut evaluations,
            &mut cancelled,
        )?;
        let (fine, second) = probe(
            &half,
            reservoirs,
            h2,
            window,
            controls,
            &mut work.half,
            &mut evaluations,
            &mut cancelled,
        )?;
        check_clock(&half, middle)?;
        check_clock(&coarse, end)?;
        check_clock(&fine, end)?;
        let (u, v) = (h1 / span, h2 / span);
        let factor = (u * u + v * v) / (2.0 * u * v);
        let error = factor
            * distance(
                &state,
                &coarse,
                &fine,
                [&full, &first, &second],
                controls,
                &mut cancelled,
            )?;
        if !error.is_finite() || error > 1.0 {
            let smaller = (0.5 * span).max(controls.minimum_step_s);
            if smaller >= span {
                return Err(Error::Accuracy {
                    time_s: time,
                    step_s: span,
                    error_ratio: error,
                });
            }
            rejected += 1;
            step = smaller;
            continue;
        }
        for report in [&first, &second] {
            for face in 0..model.faces.len() {
                if face % 128 == 0 {
                    poll(&mut cancelled)?;
                }
                let amount = face_amount(report, face, model.area_m2)?;
                add(&mut signed[face], amount)?;
                add(&mut gross[face], amount.abs())?;
            }
        }
        accepted += 1;
        minimum = minimum.min(h1.min(h2));
        maximum = maximum.max(h1.max(h2));
        maximum_error = maximum_error.max(error);
        state = fine;
        last = Some(second);
        step = if error <= 0.25 {
            (2.0 * span).min(controls.maximum_step_s)
        } else {
            span
        };
        step = step.max(controls.minimum_step_s);
    }
    check_clock(&state, target)?;
    let mut change = Sum::default();
    let mut maximum_cell_residual = 0.0_f64;
    for i in 0..model.cells.len() {
        if i % 128 == 0 {
            poll(&mut cancelled)?;
        }
        let delta = finite(
            (model.cells[i].capacity_mol_m2_pa * model.area_m2)
                * (state.membrane().partial_pressures_pa[i]
                    - initial.membrane().partial_pressures_pa[i]),
        )?;
        add(&mut change, delta)?;
        let residual = finite(delta - (signed[i].value - signed[i + 1].value))?;
        balance(Some(i), residual, controls.amount_balance_tolerance_mol)?;
        maximum_cell_residual = maximum_cell_residual.max(residual.abs());
    }
    let last_face = model.faces.len() - 1;
    let residual = finite(change.value - (signed[0].value - signed[last_face].value))?;
    balance(None, residual, controls.amount_balance_tolerance_mol)?;
    let finite_reservoirs = if let Some(before) = initial.reservoirs() {
        let after = state.reservoirs().expect("fixed finite reservoirs");
        let mut changes = [0.0; 2];
        let mut residuals = [0.0; 2];
        let outward = [signed[0].value, -signed[last_face].value];
        for side in 0..2 {
            changes[side] =
                finite(before.capacities[side] * (after.pressures[side] - before.pressures[side]))?;
            residuals[side] = finite(changes[side] + outward[side])?;
            if residuals[side].abs() > controls.amount_balance_tolerance_mol {
                return Err(Error::ReservoirBalance {
                    reservoir: side,
                    residual_mol: residuals[side],
                    tolerance_mol: controls.amount_balance_tolerance_mol,
                });
            }
        }
        let mut total = Sum::default();
        add(&mut total, changes[0])?;
        add(&mut total, change.value)?;
        add(&mut total, changes[1])?;
        balance(None, total.value, controls.amount_balance_tolerance_mol)?;
        Some(FiniteMembraneAdvanceBalance {
            pressures_before_pa: before.pressures,
            pressures_after_pa: after.pressures,
            amount_changes_mol: changes,
            reservoir_residuals_mol: residuals,
            closed_amount_residual_mol: total.value,
        })
    } else {
        None
    };
    let report = MembraneDiffusionAdvanceReport {
        interval_s: [start, target],
        accepted_intervals: accepted,
        accepted_substeps: 2 * accepted,
        rejected_intervals: rejected,
        evaluations,
        minimum_accepted_substep_s: minimum,
        maximum_accepted_substep_s: maximum,
        maximum_accepted_error_ratio: maximum_error,
        suggested_step_s: step,
        face_amounts_mol: signed.into_iter().map(|s| s.value).collect(),
        gross_face_amounts_mol: gross.into_iter().map(|s| s.value).collect(),
        stored_before_mol: initial
            .membrane()
            .stored_amount_mol()
            .map_err(|_| Error::Invalid("unrepresentable original sorbed amount"))?,
        stored_after_mol: state
            .membrane()
            .stored_amount_mol()
            .map_err(|_| Error::Invalid("unrepresentable final sorbed amount"))?,
        stored_change_mol: change.value,
        amount_residual_mol: residual,
        maximum_cell_amount_residual_mol: maximum_cell_residual,
        finite_reservoirs,
        last_substep: last.expect("positive completed window contains native steps"),
    };
    poll(&mut cancelled)?;
    Ok(MembraneDiffusionAdvanceCandidate { state, report })
}
fn balance<E>(cell: Option<usize>, residual: f64, tolerance: f64) -> Result<(), E> {
    if residual.abs() > tolerance {
        Err(Error::AmountBalance {
            cell,
            residual_mol: residual,
            tolerance_mol: tolerance,
        })
    } else {
        Ok(())
    }
}

impl MembraneDiffusionState {
    pub fn adaptive_workspace(
        &self,
    ) -> core::result::Result<MembraneDiffusionAdvanceWorkspace, MembraneDiffusionError> {
        Ok(MembraneDiffusionAdvanceWorkspace {
            full: MembraneDiffusionWorkspace::new(self.model.cell_count())?,
            half: MembraneDiffusionWorkspace::new(self.model.cell_count())?,
        })
    }
    /// Complete a fixed-boundary duration with pressure, stored-amount and
    /// integrated-throughput accuracy estimates. All native source-free model
    /// limits remain. On any refusal this state and its clock are unchanged.
    pub fn advance_candidate(
        &self,
        reservoirs_pa: [f64; 2],
        duration_s: f64,
        controls: MembraneDiffusionAdvanceControls,
        work: &mut MembraneDiffusionAdvanceWorkspace,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MembraneDiffusionAdvanceCandidate> {
        advance(self, reservoirs_pa, duration_s, controls, work, cancelled)
    }
}
impl FiniteMembraneDiffusionState {
    pub fn adaptive_workspace(
        &self,
    ) -> core::result::Result<
        MembraneDiffusionAdvanceWorkspace<FiniteMembraneDiffusionWorkspace>,
        MembraneDiffusionError,
    > {
        let n = self.membrane().model().cell_count();
        Ok(MembraneDiffusionAdvanceWorkspace {
            full: FiniteMembraneDiffusionWorkspace::new(n)?,
            half: FiniteMembraneDiffusionWorkspace::new(n)?,
        })
    }
    /// Evolve both finite pressure-capacity reservoirs and the sorbed film on
    /// the same accepted timeline. No EOS, thermostat or caloric law is inferred.
    pub fn advance_candidate(
        &self,
        duration_s: f64,
        controls: MembraneDiffusionAdvanceControls,
        work: &mut MembraneDiffusionAdvanceWorkspace<FiniteMembraneDiffusionWorkspace>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MembraneDiffusionAdvanceCandidate<Self>> {
        advance(self, [0.0; 2], duration_s, controls, work, cancelled)
    }
}

#[cfg(test)]
mod tests;
