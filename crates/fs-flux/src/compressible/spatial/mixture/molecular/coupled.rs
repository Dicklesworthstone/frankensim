use super::super::{
    BodyForce3d, Error, ForceChange3d, ForceControl3d, ForcingWorkspace3d, StepReport, forcing,
};
use super::*;
pub use crate::compressible::mixture::molecular::MolecularStepControl;

/// Complete accepted volumetric multicomponent molecular/advection interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolecularStepReport3d {
    /// Common physical duration; each molecular half represents half.
    pub elapsed: Time,
    /// Complete rejected trials before the accepted interval.
    pub coupling_retries: usize,
    /// First molecular half interval.
    pub first: MolecularReport,
    /// Central original spatial advection, including supplied body forcing.
    pub advection: StepReport,
    /// Second molecular half interval.
    pub second: MolecularReport,
    /// Original extensive field.
    pub before: Totals,
    /// Final extensive field.
    pub after: Totals,
    /// Complete inward molecular and advective boundary exchange.
    pub boundary_change: Totals,
    /// Body impulse and signed work, separate from boundary exchange.
    pub source_change: ForceChange3d,
    /// `after - before - boundary_change - source_change`.
    pub conservation_residual: Totals,
    /// Physical chemical-plus-kinetic energy residual from conservative differences.
    pub absolute_energy_residual_j: f64,
}

/// Atomic nonreacting multicomponent viscous, thermal and advective flow.
///
/// The prepared molecular operator retains its original source, physical walls
/// and optional SGS configuration. Its requested-duration halves surround the
/// original three-dimensional mixture Euler transaction. No chemistry, frozen
/// composition surrogate or additional viscous-heating source is introduced.
///
/// Shortened or refused phases restart the complete interval from accepted
/// cells. Final cancellation preserves both owners' published state and loads.
/// Symmetric composition permits second-order smooth evolution when the
/// constituent paths retain their order; limiting and first-order retries can
/// reduce it. No entropy, turbulence or experimental-validation claim is added.
pub struct MolecularEuler3d<L: TransportLaw> {
    accepted: MixtureEuler3d,
    trial: MixtureEuler3d,
    molecular: MolecularTransport3d<L>,
    boundary_impulses: Vec<f64>,
    boundary_trial: Vec<f64>,
    maximum_workspace_bytes: usize,
    last_step: Option<MolecularStepReport3d>,
}
impl<L: TransportLaw> MolecularEuler3d<L> {
    /// Both full fields, prepared molecular/SGS buffers and complete boundary
    /// buffers, including incoming excess capacity. Immutable sources, bounded
    /// stack scratch and provider storage are excluded.
    pub fn required_workspace_bytes(
        flow: &MixtureEuler3d,
        molecular: &MolecularTransport3d<L>,
    ) -> Result<usize, DissipationError> {
        flow.workspace_bytes()
            .checked_add(MixtureEuler3d::required_workspace_bytes(
                flow.gas(),
                flow.grid(),
            )?)
            .and_then(|bytes| bytes.checked_add(molecular.workspace_bytes()))
            .and_then(|bytes| {
                bytes.checked_add(
                    flow.boundary_impulses
                        .len()
                        .checked_mul(2 * size_of::<f64>())?,
                )
            })
            .ok_or(DissipationError::Invalid(
                "coupled spatial mixture workspace size overflow",
            ))
    }

    /// Bind the original field and its already-prepared molecular operator.
    /// The stored complete budget also admits externally supplied forcing scratch,
    /// even for an exactly zero source. Both owners remain immutable to callers.
    pub fn new(
        flow: MixtureEuler3d,
        molecular: MolecularTransport3d<L>,
        maximum_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DissipationError> {
        poll(&mut cancelled)?;
        if molecular.gas() != flow.gas()
            || molecular.grid != flow.grid
            || molecular.boundaries != flow.boundaries
        {
            return Err(DissipationError::Invalid(
                "coupled spatial mixture molecular binding differs",
            ));
        }
        let required = Self::required_workspace_bytes(&flow, &molecular)?;
        admit_workspace(required, maximum_workspace_bytes)?;
        let mut states = allocated(flow.states.len(), 0.0)?;
        copy_values(&flow.states, &mut states, &mut cancelled)?;
        let mut trial = MixtureEuler3d::new(
            flow.gas.clone(),
            flow.grid,
            states,
            flow.boundaries,
            MixtureEuler3d::required_workspace_bytes(flow.gas(), flow.grid())?,
            &mut cancelled,
        )?;
        copy_physical(&flow, &mut trial, &mut cancelled)?;
        let mut boundary_impulses = allocated(flow.boundary_impulses.len(), 0.0)?;
        copy_values(
            &flow.boundary_impulses,
            &mut boundary_impulses,
            &mut cancelled,
        )?;
        let boundary_trial = allocated(flow.boundary_impulses.len(), 0.0)?;
        let result = Self {
            accepted: flow,
            trial,
            molecular,
            boundary_impulses,
            boundary_trial,
            maximum_workspace_bytes,
            last_step: None,
        };
        admit_workspace(result.workspace_bytes(), maximum_workspace_bytes)?;
        poll(&mut cancelled)?;
        Ok(result)
    }

    /// Complete accepted field and clock. Its Euler metadata describes the
    /// middle phase; this owner's report and face loads describe the whole split.
    pub const fn flow(&self) -> &MixtureEuler3d {
        &self.accepted
    }
    /// Bound original operator; its published report describes the second half.
    pub const fn molecular(&self) -> &MolecularTransport3d<L> {
        &self.molecular
    }
    /// Actual owned numerical capacities, following state-buffer exchanges.
    pub fn workspace_bytes(&self) -> usize {
        self.accepted.workspace_bytes()
            + self.trial.workspace_bytes()
            + self.molecular.workspace_bytes()
            + (self.boundary_impulses.capacity() + self.boundary_trial.capacity())
                * size_of::<f64>()
    }
    /// Last completely accepted physical interval.
    pub const fn last_step(&self) -> Option<&MolecularStepReport3d> {
        self.last_step.as_ref()
    }
    /// Complete outward exchange through each physical face.
    pub fn boundary_impulses(
        &self,
        side: BoundarySide,
    ) -> impl ExactSizeIterator<Item = BoundaryImpulse<'_>> {
        boundary_view(
            &self.boundary_impulses,
            self.accepted.grid(),
            self.accepted.gas().species_count(),
            side,
        )
    }
    /// Recover both accepted owners and their original wall/source configuration.
    pub fn into_parts(self) -> (MixtureEuler3d, MolecularTransport3d<L>) {
        (self.accepted, self.molecular)
    }

    /// Advance at most the requested duration, restarting both molecular halves
    /// when the central interval shortens. Exhausted whole-step retries return
    /// a `DissipationError::StepFailure` with zero committed elapsed time.
    pub fn step(
        &mut self,
        maximum_duration: Time,
        controls: MolecularStepControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularStepReport3d, DissipationError> {
        self.step_impl(maximum_duration, controls, None, cancelled)
    }

    /// Drive the same original three-dimensional split with held acceleration or
    /// force density. Source work stays separate from molecular and wall energy.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_forcing(
        &mut self,
        maximum_duration: Time,
        controls: MolecularStepControl,
        force: BodyForce3d<'_>,
        force_controls: ForceControl3d,
        workspace: &mut ForcingWorkspace3d,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularStepReport3d, DissipationError> {
        self.step_impl(
            maximum_duration,
            controls,
            Some((force, force_controls, workspace)),
            cancelled,
        )
    }

    fn step_impl(
        &mut self,
        maximum_duration: Time,
        controls: MolecularStepControl,
        mut forcing: Option<forcing::Context<'_, '_>>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularStepReport3d, DissipationError> {
        poll(&mut cancelled)?;
        let forcing_bytes = forcing
            .as_ref()
            .map_or(0, |(_, _, work)| work.workspace_bytes());
        let active = if let Some((force, control, workspace)) = forcing.as_mut() {
            control.admit::<Error>()?;
            workspace.admit(&self.accepted)?;
            force.admit::<Error>(self.accepted.cell_count(), &mut cancelled)?
        } else {
            false
        };
        if !active {
            forcing = None;
        }
        let bytes =
            self.workspace_bytes()
                .checked_add(forcing_bytes)
                .ok_or(DissipationError::Invalid(
                    "coupled spatial mixture workspace size overflow",
                ))?;
        admit_workspace(bytes, self.maximum_workspace_bytes)?;
        controls.molecular.validate()?;
        if !positive(maximum_duration.value())
            || controls.maximum_coupling_retries > 32
            || !positive(controls.advection.cfl)
            || controls.advection.cfl > 0.45
            || controls.advection.maximum_halvings > 32
        {
            return Err(DissipationError::Invalid(
                "coupled spatial mixture duration or controls",
            ));
        }
        let n = self.accepted.gas.species_count();
        let width = n + 4;
        let before = checked_totals(
            &self.accepted.states,
            &self.accepted.grid,
            n,
            &mut cancelled,
        )?;
        let mut rate: f64 = 0.0;
        for row in self.accepted.states.chunks_exact(width) {
            let p = cell(row, n).primitive(&self.accepted.gas, &mut cancelled)?;
            if let Some((_, control, _)) = forcing.as_ref() {
                control.admit_velocity::<Error>(p.velocity, p.sound_speed)?;
            }
            rate = rate.max(forcing::acoustic_rate(p, &self.accepted.grid)?);
        }
        // Evolved source stages and faces remain authoritative; this initial
        // headroom avoids a predictably excessive first molecular interval.
        let mut dt = maximum_duration
            .value()
            .min(0.8 * controls.advection.cfl / rate);
        let mut last = DissipationError::Invalid("no coupled spatial mixture attempt");
        let mut last_dt = dt;
        for retry in 0..=controls.maximum_coupling_retries {
            poll(&mut cancelled)?;
            if !positive(0.5 * dt)
                || !(self.accepted.time + dt).is_finite()
                || self.accepted.time + dt <= self.accepted.time
            {
                return Err(DissipationError::Invalid(
                    "unrepresentable coupled spatial mixture increment",
                ));
            }
            copy_physical(&self.accepted, &mut self.trial, &mut cancelled)?;
            let mut shortened = None;
            let attempt = (|| -> Result<MolecularStepReport3d, DissipationError> {
                let first = self.molecular.advance_unpublished(
                    &mut self.trial,
                    Time::new(0.5 * dt),
                    controls.molecular,
                    self.maximum_workspace_bytes,
                    &mut cancelled,
                )?;
                copy_values(
                    self.molecular.trial_impulses(),
                    &mut self.boundary_trial,
                    &mut cancelled,
                )?;
                // Admit the immutable source/context once at the outer owner,
                // preserving native stage checks inside the original Euler step.
                let advection = if let Some((force, control, workspace)) = forcing.as_mut() {
                    self.trial.step_impl(
                        Time::new(dt),
                        controls.advection,
                        Some((*force, *control, &mut **workspace)),
                        &mut cancelled,
                    )?
                } else {
                    self.trial
                        .step(Time::new(dt), controls.advection, &mut cancelled)?
                };
                if advection.elapsed.value() != dt {
                    shortened = Some(advection.elapsed.value());
                    return Err(DissipationError::Invalid(
                        "advection shortened coupled spatial mixture interval",
                    ));
                }
                add_impulses(
                    &self.trial.boundary_impulses,
                    &mut self.boundary_trial,
                    &mut cancelled,
                )?;
                let second = self.molecular.advance_unpublished(
                    &mut self.trial,
                    Time::new(0.5 * dt),
                    controls.molecular,
                    self.maximum_workspace_bytes,
                    &mut cancelled,
                )?;
                add_impulses(
                    self.molecular.trial_impulses(),
                    &mut self.boundary_trial,
                    &mut cancelled,
                )?;
                if let Some((_, control, _)) = forcing.as_ref() {
                    let final_rate = forcing::admit_field(
                        &self.trial.gas,
                        &self.trial.grid,
                        &self.trial.states,
                        *control,
                        &mut cancelled,
                    )?;
                    stepping::check_cfl::<Error>(final_rate * dt, controls.advection.cfl)?;
                }
                let after =
                    checked_totals(&self.trial.states, &self.trial.grid, n, &mut cancelled)?;
                let outward = sum_rows(&self.boundary_trial, width, 1.0, &mut cancelled)?;
                let boundary_change = Totals::from_row(&outward.map(|value| -value), n);
                let source_change = advection.source_change;
                let mut residual = [0.0; WIDTH_LIMIT];
                for (d, value) in residual[..width].iter_mut().enumerate() {
                    *value =
                        after.component(d) - before.component(d) - boundary_change.component(d);
                }
                for axis in 0..3 {
                    residual[n + axis] -= source_change.impulse[axis];
                }
                residual[n + 3] -= source_change.work;
                if !residual[..width].iter().all(|value| value.is_finite()) {
                    return Err(DissipationError::Invalid(
                        "coupled spatial mixture balance overflow",
                    ));
                }
                let conservation_residual = Totals::from_row(&residual[..width], n);
                let absolute_energy_residual_j =
                    conservation_residual.absolute_energy(&self.trial.gas)?;
                poll(&mut cancelled)?;
                Ok(MolecularStepReport3d {
                    elapsed: Time::new(dt),
                    coupling_retries: retry,
                    first,
                    advection,
                    second,
                    before,
                    after,
                    boundary_change,
                    source_change,
                    conservation_residual,
                    absolute_energy_residual_j,
                })
            })();
            match attempt {
                Ok(report) => {
                    std::mem::swap(&mut self.accepted, &mut self.trial);
                    std::mem::swap(&mut self.boundary_impulses, &mut self.boundary_trial);
                    self.molecular.publish(report.second);
                    self.last_step = Some(report);
                    return Ok(report);
                }
                Err(DissipationError::Cancelled) => return Err(DissipationError::Cancelled),
                Err(error) => {
                    last = error;
                    last_dt = dt;
                    dt = shortened.unwrap_or(0.5 * dt);
                }
            }
        }
        Err(DissipationError::StepFailure {
            elapsed_time_s: 0.0,
            step_s: last_dt,
            last: Box::new(last),
        })
    }
}

fn copy_physical(
    source: &MixtureEuler3d,
    target: &mut MixtureEuler3d,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), DissipationError> {
    copy_values(&source.states, &mut target.states, cancelled)?;
    copy_values(
        &source.boundary_impulses,
        &mut target.boundary_impulses,
        cancelled,
    )?;
    target.time = source.time;
    target.last_step = source.last_step;
    Ok(())
}

fn add_impulses(
    source: &[f64],
    target: &mut [f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), DissipationError> {
    for (source, target) in source.chunks(4096).zip(target.chunks_mut(4096)) {
        poll(cancelled)?;
        for (source, target) in source.iter().zip(target) {
            *target += source;
        }
    }
    Ok(())
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
