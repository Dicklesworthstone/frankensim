//! Competitive breakthrough windows through the original membrane time driver.
//!
//! Every native probe advances all species simultaneously. Original per-species
//! fields remain separate in pressure, mobile amount, occupied-site amount and
//! integrated-face comparisons. Only accepted half-step histories survive;
//! rejected competition, source failures and late cancellation publish nothing.
//! Coefficients, geometry and temperature are fixed during a window. The local
//! first-order estimate is not a spatial or whole-window error certificate.

use super::super::adaptive::{
    self as driver, FiniteMembraneAdvanceBalance, MembraneDiffusionAdvanceControls,
    MembraneDiffusionAdvanceError as Error, MembraneDiffusionAdvanceReport,
    MembraneDiffusionAdvanceWorkspace, Native, Reservoirs,
};
use super::super::trapping::adaptive::TrappingMembraneAdvanceReport;
use super::*;

type AdvanceResult<T, E = CompetitiveMembraneError> = core::result::Result<T, Error<E>>;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CompetitiveMembraneAdvanceControls {
    pub temporal: MembraneDiffusionAdvanceControls,
    pub maximum_iterations: usize,
    pub maximum_backtracks: usize,
    pub maximum_linear_sweeps: usize,
    pub linear_relative_tolerance: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct CompetitiveMembraneAdvanceReport {
    /// Canonical SpeciesId order. Every row has the SAME accepted timeline and
    /// work counters, but its own original cell/bath/face balances and scales.
    pub species: Vec<(SpeciesId, TrappingMembraneAdvanceReport)>,
    pub last_substep: CompetitiveMembraneStep,
}

#[derive(Clone, Debug, PartialEq)]
pub struct CompetitiveMembraneAdvanceCandidate<S = CompetitiveMembraneDiffusionState> {
    state: S,
    report: CompetitiveMembraneAdvanceReport,
}
impl<S> CompetitiveMembraneAdvanceCandidate<S> {
    pub const fn state(&self) -> &S {
        &self.state
    }
    pub const fn report(&self) -> &CompetitiveMembraneAdvanceReport {
        &self.report
    }
    pub fn into_parts(self) -> (S, CompetitiveMembraneAdvanceReport) {
        (self.state, self.report)
    }
}

/// Source-bound owners can retain their own admission and actual-gas rebuild
/// on EVERY joint probe. The temporal driver never replaces physical owners.
pub(crate) trait CompetitiveNative: Clone {
    type Workspace;
    type Error: From<CompetitiveMembraneError>;
    fn competition(&self) -> &CompetitiveMembraneDiffusionState;
    fn finite_reservoirs(&self) -> Option<([f64; 2], &[[f64; 2]])> {
        None
    }
    fn is_cancelled(error: &Self::Error) -> bool;
    /// Source wrappers can compare additional ORIGINAL inventories (for
    /// example actual EOS gas composition), retaining the joint raw estimator.
    fn additional_distance(
        &self,
        _coarse: &Self,
        _fine: &Self,
        _controls: MembraneDiffusionAdvanceControls,
        _cancelled: &mut impl FnMut() -> bool,
    ) -> AdvanceResult<f64, Self::Error> {
        Ok(0.0)
    }
    fn competitive_step(
        &self,
        pressures: &[[f64; 2]],
        duration: f64,
        solve: CompetitiveMembraneSolve,
        work: &mut Self::Workspace,
        cancelled: &mut impl FnMut() -> bool,
    ) -> core::result::Result<(Self, CompetitiveMembraneStep), Self::Error>;
}
impl CompetitiveNative for CompetitiveMembraneDiffusionState {
    type Workspace = CompetitiveMembraneWorkspace;
    type Error = CompetitiveMembraneError;
    fn competition(&self) -> &Self {
        self
    }
    fn is_cancelled(e: &Self::Error) -> bool {
        *e == CompetitiveMembraneError::Cancelled
    }
    fn competitive_step(
        &self,
        pressures: &[[f64; 2]],
        duration: f64,
        solve: CompetitiveMembraneSolve,
        work: &mut Self::Workspace,
        cancelled: &mut impl FnMut() -> bool,
    ) -> core::result::Result<(Self, CompetitiveMembraneStep), Self::Error> {
        Ok(self
            .step_candidate(pressures, duration, solve, work, cancelled)?
            .into_parts())
    }
}
impl CompetitiveNative for FiniteCompetitiveMembraneDiffusionState {
    type Workspace = CompetitiveMembraneWorkspace;
    type Error = CompetitiveMembraneError;
    fn competition(&self) -> &CompetitiveMembraneDiffusionState {
        self.membrane()
    }
    fn finite_reservoirs(&self) -> Option<([f64; 2], &[[f64; 2]])> {
        Some((
            self.capacities().map(MolarPressureCapacity::value_si),
            self.reservoir_partial_pressures_pa(),
        ))
    }
    fn is_cancelled(e: &Self::Error) -> bool {
        *e == CompetitiveMembraneError::Cancelled
    }
    fn competitive_step(
        &self,
        _: &[[f64; 2]],
        duration: f64,
        solve: CompetitiveMembraneSolve,
        work: &mut Self::Workspace,
        cancelled: &mut impl FnMut() -> bool,
    ) -> core::result::Result<(Self, CompetitiveMembraneStep), Self::Error> {
        Ok(self
            .step_candidate(duration, solve, work, cancelled)?
            .into_parts())
    }
}

// A borrowed view of one ORIGINAL channel, solely for the shared audit and
// distance routines. It cannot independently advance a competing species.
struct View<'a, E> {
    channel: &'a Channel,
    baths: Option<Reservoirs>,
    error: core::marker::PhantomData<E>,
}
impl<E> Clone for View<'_, E> {
    fn clone(&self) -> Self {
        Self {
            channel: self.channel,
            baths: self.baths,
            error: core::marker::PhantomData,
        }
    }
}
fn view<S: CompetitiveNative>(state: &S, species: usize) -> View<'_, S::Error> {
    View {
        channel: &state.competition().channels[species],
        baths: state
            .finite_reservoirs()
            .map(|(capacities, pressures)| Reservoirs {
                capacities,
                pressures: pressures[species],
            }),
        error: core::marker::PhantomData,
    }
}
impl<E: From<CompetitiveMembraneError>> Native for View<'_, E> {
    type Workspace = ();
    type Error = E;
    fn membrane(&self) -> &MembraneDiffusionState {
        &self.channel.mobile
    }
    fn immobile_amounts(&self) -> Option<&[f64]> {
        Some(&self.channel.trapped)
    }
    fn reservoirs(&self) -> Option<Reservoirs> {
        self.baths
    }
    fn is_cancelled(_: &E) -> bool {
        false
    }
    fn native_step(
        &self,
        _: [f64; 2],
        _: f64,
        _: MembraneDiffusionSolve,
        _: &mut (),
        _: &mut impl FnMut() -> bool,
    ) -> core::result::Result<(Self, MembraneDiffusionStep), E> {
        Err(E::from(invalid(
            "a competitive audit view cannot advance independently",
        )))
    }
}

struct Ledger {
    signed: Vec<Sum>,
    gross: Vec<Sum>,
}
impl Ledger {
    fn new(n: usize) -> Self {
        Self {
            signed: (0..n).map(|_| Sum::default()).collect(),
            gross: (0..n).map(|_| Sum::default()).collect(),
        }
    }
    fn add(
        &mut self,
        report: &MembraneDiffusionStep,
        area: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> core::result::Result<(), CompetitiveMembraneError> {
        for (face, &flux) in report.face_molar_fluxes_mol_m2_s.iter().enumerate() {
            if face % 128 == 0 {
                poll(cancelled)?;
            }
            let amount = finite((flux * area) * report.duration_s)?;
            self.signed[face].add(amount)?;
            self.gross[face].add(amount.abs())?;
        }
        Ok(())
    }
}
impl Clone for Ledger {
    fn clone(&self) -> Self {
        let copy = |xs: &[Sum]| {
            xs.iter()
                .map(|s| Sum {
                    value: s.value,
                    correction: s.correction,
                })
                .collect()
        };
        Self {
            signed: copy(&self.signed),
            gross: copy(&self.gross),
        }
    }
}

#[derive(Clone)]
struct Controlled<S> {
    state: S,
    boundaries: Arc<[[f64; 2]]>,
    controls: CompetitiveMembraneAdvanceControls,
    minimum_pressure_capacity: f64,
    ledgers: Vec<Ledger>,
    last: Option<CompetitiveMembraneStep>,
    transport: Vec<MembraneDiffusionStep>,
}
fn transport(
    step: &TrappingMembraneStep,
) -> core::result::Result<MembraneDiffusionStep, CompetitiveMembraneError> {
    Ok(MembraneDiffusionStep {
        duration_s: step.duration_s,
        reservoir_partial_pressures_pa: step.transport.reservoir_partial_pressures_pa,
        face_molar_fluxes_mol_m2_s: step.transport.face_molar_fluxes_mol_m2_s.clone(),
        face_partial_pressures_pa: step.transport.face_partial_pressures_pa.clone(),
        stored_before_mol: finite(step.mobile_before_mol + step.trapped_before_mol)?,
        stored_after_mol: finite(step.mobile_after_mol + step.trapped_after_mol)?,
        stored_change_mol: finite(step.mobile_change_mol + step.trapped_change_mol)?,
        left_reservoir_to_membrane_mol: step.boundary_transfers_mol[0],
        membrane_to_right_reservoir_mol: step.boundary_transfers_mol[1],
        amount_residual_mol: step.membrane_amount_residual_mol,
        maximum_cell_amount_residual_mol: step.maximum_cell_amount_residual_mol,
    })
}
fn audit<E>(residual: f64, tolerance: f64, cell: Option<usize>) -> AdvanceResult<(), E> {
    if !residual.is_finite() || residual.abs() > tolerance {
        Err(Error::AmountBalance {
            cell,
            residual_mol: residual,
            tolerance_mol: tolerance,
        })
    } else {
        Ok(())
    }
}
fn arithmetic<E: From<CompetitiveMembraneError>>(e: MembraneDiffusionError) -> Error<E> {
    Error::Native(E::from(e.into()))
}

impl<S: CompetitiveNative> Native for Controlled<S> {
    type Workspace = S::Workspace;
    type Error = S::Error;
    fn membrane(&self) -> &MembraneDiffusionState {
        &self.state.competition().channels[0].mobile
    }
    fn immobile_amounts(&self) -> Option<&[f64]> {
        Some(&self.state.competition().channels[0].trapped)
    }
    fn reservoirs(&self) -> Option<Reservoirs> {
        view(&self.state, 0).baths
    }
    fn is_cancelled(e: &Self::Error) -> bool {
        S::is_cancelled(e)
    }
    fn accuracy_controls(
        &self,
        requested: MembraneDiffusionSolve,
        mut limit: f64,
    ) -> MembraneDiffusionSolve {
        // All channels have the same inventory count, so the driver's first
        // field supplies the common duration-weighted share. Preserve stricter
        // pressure requirements of every other species' distinct Henry capacity.
        limit = limit.min(
            (0.01 * self.controls.temporal.absolute_pressure_tolerance_pa)
                * self.minimum_pressure_capacity,
        );
        if requested.absolute_amount_tolerance_mol > 0.0 {
            MembraneDiffusionSolve {
                absolute_amount_tolerance_mol: requested
                    .absolute_amount_tolerance_mol
                    .min(0.5 * limit),
                relative_amount_tolerance: 0.0,
            }
        } else {
            requested
        }
    }
    fn native_step(
        &self,
        _: [f64; 2],
        duration: f64,
        balance: MembraneDiffusionSolve,
        work: &mut Self::Workspace,
        cancelled: &mut impl FnMut() -> bool,
    ) -> core::result::Result<(Self, MembraneDiffusionStep), Self::Error> {
        let (state, step) = self.state.competitive_step(
            &self.boundaries,
            duration,
            CompetitiveMembraneSolve {
                balance,
                maximum_iterations: self.controls.maximum_iterations,
                maximum_backtracks: self.controls.maximum_backtracks,
                maximum_linear_sweeps: self.controls.maximum_linear_sweeps,
                linear_relative_tolerance: self.controls.linear_relative_tolerance,
            },
            work,
            cancelled,
        )?;
        let mut ledgers = self.ledgers.clone();
        let mut reports = Vec::with_capacity(step.species.len());
        for (s, (_, row)) in step.species.iter().enumerate() {
            let report = transport(row).map_err(S::Error::from)?;
            ledgers[s]
                .add(
                    &report,
                    state.competition().channels[s].mobile.model.area_m2,
                    cancelled,
                )
                .map_err(S::Error::from)?;
            reports.push(report);
        }
        let first = reports[0].clone();
        Ok((
            Self {
                state,
                ledgers,
                boundaries: self.boundaries.clone(),
                controls: self.controls,
                minimum_pressure_capacity: self.minimum_pressure_capacity,
                last: Some(step),
                transport: reports,
            },
            first,
        ))
    }
    fn additional_algebraic_balance(
        &self,
        new: &Self,
        window: f64,
        controls: MembraneDiffusionAdvanceControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> AdvanceResult<(), Self::Error> {
        for s in 0..self.state.competition().channels.len() {
            if cancelled() {
                return Err(Error::Cancelled);
            }
            let old = view(&self.state, s);
            let next = view(&new.state, s);
            let report = &new.transport[s];
            if s > 0 {
                driver::algebraic_balance(&old, &next, report, window, controls, cancelled)?;
            }
            let (_, share) = driver::algebraic_budget(&old, report.duration_s, window, controls)?;
            let mut limit = share;
            for (i, cell) in old.channel.mobile.model.cells.iter().enumerate() {
                if i % 128 == 0 && cancelled() {
                    return Err(Error::Cancelled);
                }
                limit = limit.min(driver::pressure_allowance::<Self::Error>(
                    share,
                    cell.capacity_mol_m2_pa * old.channel.mobile.model.area_m2,
                    controls,
                )?);
            }
            // The native report includes its independently evaluated kinetic
            // residual. A balanced mobile+trapped sum cannot hide lost uptake.
            audit(report.maximum_cell_amount_residual_mol, limit, None)?;
        }
        Ok(())
    }
    fn additional_distance(
        &self,
        coarse: &Self,
        first_half: &Self,
        fine: &Self,
        controls: MembraneDiffusionAdvanceControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> AdvanceResult<f64, Self::Error> {
        let mut error = 0.0_f64;
        for s in 1..self.state.competition().channels.len() {
            if cancelled() {
                return Err(Error::Cancelled);
            }
            error = error.max(driver::distance(
                &view(&self.state, s),
                &view(&coarse.state, s),
                &view(&fine.state, s),
                [
                    &coarse.transport[s],
                    &first_half.transport[s],
                    &fine.transport[s],
                ],
                controls,
                cancelled,
            )?);
        }
        // Shared vacancies are an additional amount coordinate. Sum occupancy
        // CHANGES to avoid cancellation against a nearly saturated capacity.
        for i in 0..self.state.competition().capacities.len() {
            if i % 128 == 0 && cancelled() {
                return Err(Error::Cancelled);
            }
            let (mut a, mut b, mut difference) = (Sum::default(), Sum::default(), Sum::default());
            for s in 0..self.state.competition().channels.len() {
                let initial = self.state.competition().channels[s].trapped[i];
                let full = coarse.state.competition().channels[s].trapped[i];
                let half = fine.state.competition().channels[s].trapped[i];
                a.add(full - initial).map_err(arithmetic::<Self::Error>)?;
                b.add(half - initial).map_err(arithmetic::<Self::Error>)?;
                difference
                    .add(full - half)
                    .map_err(arithmetic::<Self::Error>)?;
            }
            let allowance = finite(
                controls.absolute_amount_tolerance_mol
                    + controls.relative_tolerance * a.value.abs().max(b.value.abs()),
            )
            .map_err(arithmetic::<Self::Error>)?;
            error = error.max(difference.value.abs() / allowance);
        }
        Ok(error.max(self.state.additional_distance(
            &coarse.state,
            &fine.state,
            controls,
            cancelled,
        )?))
    }
}

fn species_report<S: CompetitiveNative>(
    initial: &S,
    final_state: &S,
    s: usize,
    ledger: &Ledger,
    mut report: MembraneDiffusionAdvanceReport,
    last: MembraneDiffusionStep,
    tolerance: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> AdvanceResult<TrappingMembraneAdvanceReport, S::Error> {
    let old = view(initial, s);
    let new = view(final_state, s);
    let model = &old.channel.mobile.model;
    let (mut mobile, mut trapped, mut before, mut after) = (
        Sum::default(),
        Sum::default(),
        Sum::default(),
        Sum::default(),
    );
    let mut max_residual = 0.0_f64;
    for i in 0..model.cells.len() {
        if i % 128 == 0 && cancelled() {
            return Err(Error::Cancelled);
        }
        let delta_mobile = finite(
            (model.cells[i].capacity_mol_m2_pa * model.area_m2)
                * (new.channel.mobile.partial_pressures_pa[i]
                    - old.channel.mobile.partial_pressures_pa[i]),
        )
        .map_err(arithmetic::<S::Error>)?;
        let delta_trapped = finite(new.channel.trapped[i] - old.channel.trapped[i])
            .map_err(arithmetic::<S::Error>)?;
        mobile.add(delta_mobile).map_err(arithmetic::<S::Error>)?;
        trapped.add(delta_trapped).map_err(arithmetic::<S::Error>)?;
        before
            .add(old.channel.trapped[i])
            .map_err(arithmetic::<S::Error>)?;
        after
            .add(new.channel.trapped[i])
            .map_err(arithmetic::<S::Error>)?;
        let residual = finite(
            (delta_mobile + delta_trapped) - (ledger.signed[i].value - ledger.signed[i + 1].value),
        )
        .map_err(arithmetic::<S::Error>)?;
        audit(residual, tolerance, Some(i))?;
        max_residual = max_residual.max(residual.abs());
    }
    let change = finite(mobile.value + trapped.value).map_err(arithmetic::<S::Error>)?;
    let outward = [
        ledger.signed[0].value,
        -ledger.signed[model.faces.len() - 1].value,
    ];
    let residual = finite(change - (outward[0] + outward[1])).map_err(arithmetic::<S::Error>)?;
    audit(residual, tolerance, None)?;
    let baths = if let Some(old) = old.baths {
        let new = new.baths.expect("fixed competitive finite baths");
        let (mut changes, mut residuals) = ([0.0; 2], [0.0; 2]);
        let mut closed = Sum::default();
        closed.add(change).map_err(arithmetic::<S::Error>)?;
        for side in 0..2 {
            changes[side] =
                finite(old.capacities[side] * (new.pressures[side] - old.pressures[side]))
                    .map_err(arithmetic::<S::Error>)?;
            residuals[side] =
                finite(changes[side] + outward[side]).map_err(arithmetic::<S::Error>)?;
            if residuals[side].abs() > tolerance {
                return Err(Error::ReservoirBalance {
                    reservoir: side,
                    residual_mol: residuals[side],
                    tolerance_mol: tolerance,
                });
            }
            closed.add(changes[side]).map_err(arithmetic::<S::Error>)?;
        }
        audit(closed.value, tolerance, None)?;
        Some(FiniteMembraneAdvanceBalance {
            pressures_before_pa: old.pressures,
            pressures_after_pa: new.pressures,
            amount_changes_mol: changes,
            reservoir_residuals_mol: residuals,
            closed_amount_residual_mol: closed.value,
        })
    } else {
        None
    };
    let mobile_before = old
        .channel
        .mobile
        .stored_amount_mol()
        .map_err(arithmetic::<S::Error>)?;
    let mobile_after = new
        .channel
        .mobile
        .stored_amount_mol()
        .map_err(arithmetic::<S::Error>)?;
    report.face_amounts_mol = ledger.signed.iter().map(|s| s.value).collect();
    report.gross_face_amounts_mol = ledger.gross.iter().map(|s| s.value).collect();
    report.stored_before_mol =
        finite(mobile_before + before.value).map_err(arithmetic::<S::Error>)?;
    report.stored_after_mol = finite(mobile_after + after.value).map_err(arithmetic::<S::Error>)?;
    report.stored_change_mol = change;
    report.amount_residual_mol = residual;
    report.maximum_cell_amount_residual_mol = max_residual;
    report.finite_reservoirs = baths;
    report.last_substep = last;
    Ok(TrappingMembraneAdvanceReport {
        transport: report,
        mobile_before_mol: mobile_before,
        mobile_after_mol: mobile_after,
        mobile_change_mol: mobile.value,
        trapped_before_mol: before.value,
        trapped_after_mol: after.value,
        trapped_change_mol: trapped.value,
    })
}

pub(crate) fn advance<S: CompetitiveNative>(
    initial: &S,
    boundaries: &[[f64; 2]],
    duration: f64,
    controls: CompetitiveMembraneAdvanceControls,
    work: &mut MembraneDiffusionAdvanceWorkspace<S::Workspace>,
    mut cancelled: impl FnMut() -> bool,
) -> AdvanceResult<CompetitiveMembraneAdvanceCandidate<S>, S::Error> {
    if cancelled() {
        return Err(Error::Cancelled);
    }
    let competition = initial.competition();
    if initial.finite_reservoirs().is_none()
        && (boundaries.len() != competition.channels.len()
            || boundaries.iter().flatten().any(|&q| !nonnegative(q)))
    {
        return Err(Error::Invalid(
            "every competitive species needs canonical finite nonnegative boundary pressures",
        ));
    }
    let mut minimum_pressure_capacity = f64::INFINITY;
    for channel in &competition.channels {
        for (i, cell) in channel.mobile.model.cells.iter().enumerate() {
            if i % 128 == 0 && cancelled() {
                return Err(Error::Cancelled);
            }
            minimum_pressure_capacity = minimum_pressure_capacity
                .min(cell.capacity_mol_m2_pa * channel.mobile.model.area_m2);
        }
    }
    let controlled = Controlled {
        state: initial.clone(),
        boundaries: Arc::from(boundaries),
        controls,
        minimum_pressure_capacity,
        ledgers: competition
            .channels
            .iter()
            .map(|c| Ledger::new(c.mobile.model.faces.len()))
            .collect(),
        last: None,
        transport: Vec::new(),
    };
    let (next, timeline) = driver::advance(
        &controlled,
        [0.0; 2],
        duration,
        controls.temporal,
        work,
        &mut cancelled,
    )?
    .into_parts();
    let mut species = Vec::with_capacity(competition.channels.len());
    for (s, channel) in competition.channels.iter().enumerate() {
        let report = species_report(
            initial,
            &next.state,
            s,
            &next.ledgers[s],
            timeline.clone(),
            next.transport[s].clone(),
            controls.temporal.amount_balance_tolerance_mol,
            &mut cancelled,
        )?;
        species.push((channel.species.clone(), report));
    }
    if cancelled() {
        return Err(Error::Cancelled);
    }
    Ok(CompetitiveMembraneAdvanceCandidate {
        state: next.state,
        report: CompetitiveMembraneAdvanceReport {
            species,
            last_substep: next.last.expect("completed positive competitive window"),
        },
    })
}

impl CompetitiveMembraneDiffusionState {
    /// The byte ceiling covers BOTH retained native workspaces together.
    pub fn adaptive_workspace(
        &self,
        maximum_bytes: usize,
    ) -> Result<MembraneDiffusionAdvanceWorkspace<CompetitiveMembraneWorkspace>> {
        Ok(MembraneDiffusionAdvanceWorkspace {
            full: self.workspace(maximum_bytes / 2)?,
            half: self.workspace(maximum_bytes / 2)?,
        })
    }
    pub fn advance_candidate(
        &self,
        partial_pressures_pa: &[[f64; 2]],
        duration_s: f64,
        controls: CompetitiveMembraneAdvanceControls,
        work: &mut MembraneDiffusionAdvanceWorkspace<CompetitiveMembraneWorkspace>,
        cancelled: impl FnMut() -> bool,
    ) -> AdvanceResult<CompetitiveMembraneAdvanceCandidate> {
        advance(
            self,
            partial_pressures_pa,
            duration_s,
            controls,
            work,
            cancelled,
        )
    }
}
impl FiniteCompetitiveMembraneDiffusionState {
    pub fn adaptive_workspace(
        &self,
        maximum_bytes: usize,
    ) -> Result<MembraneDiffusionAdvanceWorkspace<CompetitiveMembraneWorkspace>> {
        Ok(MembraneDiffusionAdvanceWorkspace {
            full: self.workspace(maximum_bytes / 2)?,
            half: self.workspace(maximum_bytes / 2)?,
        })
    }
    pub fn advance_candidate(
        &self,
        duration_s: f64,
        controls: CompetitiveMembraneAdvanceControls,
        work: &mut MembraneDiffusionAdvanceWorkspace<CompetitiveMembraneWorkspace>,
        cancelled: impl FnMut() -> bool,
    ) -> AdvanceResult<CompetitiveMembraneAdvanceCandidate<Self>> {
        advance(self, &[], duration_s, controls, work, cancelled)
    }
}

#[cfg(test)]
mod tests;
