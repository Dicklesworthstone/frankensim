use super::super::{
    BodyForce2d, Error, ForceChange2d, ForceControl2d, ForcingWorkspace2d, StepReport, forcing,
};
use super::stabilized::{IntervalFailure, StepWork};
use super::*;
use crate::compressible::StepControl;

/// Shared interval and bounded work for molecular/advection/molecular evolution.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolecularStepControl {
    /// Original conservative advection, reconstruction and CFL controls.
    pub advection: StepControl,
    /// Original requested-duration molecular substep controls.
    pub molecular: DissipationControl,
    /// Maximum whole-interval restarts after a shortened or refused phase.
    pub maximum_coupling_retries: usize,
}
impl Default for MolecularStepControl {
    fn default() -> Self {
        Self {
            advection: StepControl::default(),
            molecular: DissipationControl::default(),
            maximum_coupling_retries: 16,
        }
    }
}

/// One shared interval with RKL2 molecular work bounded across both halves.
/// The molecular maximum step caps physical accuracy, not only linear stability.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StabilizedMolecularStepControl {
    /// Original conservative advection, reconstruction and CFL controls.
    pub advection: StepControl,
    /// Original molecular macrostep controls and a whole-call operator budget.
    pub molecular: StabilizedDissipationControl,
    /// Maximum whole-interval restarts after a shortened or refused phase.
    pub maximum_coupling_retries: usize,
}
impl Default for StabilizedMolecularStepControl {
    fn default() -> Self {
        Self {
            advection: StepControl::default(),
            molecular: StabilizedDissipationControl::default(),
            maximum_coupling_retries: 16,
        }
    }
}
impl StabilizedMolecularStepControl {
    pub(in crate::compressible) const fn split(self) -> MolecularStepControl {
        MolecularStepControl {
            advection: self.advection,
            molecular: self.molecular.transport,
            maximum_coupling_retries: self.maximum_coupling_retries,
        }
    }
}

/// Complete accepted multicomponent viscous, thermal and advective interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolecularStepReport {
    /// Actual common physical duration; each molecular half represents half.
    pub elapsed: Time,
    /// Complete rejected trials before the accepted interval.
    pub coupling_retries: usize,
    /// First molecular half interval.
    pub first: MolecularReport,
    /// Central original advection, including body forcing when supplied.
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
    pub source_change: ForceChange2d,
    /// `after - before - boundary_change - source_change`.
    pub conservation_residual: Totals,
    /// Physical chemical-plus-kinetic energy residual from conservative differences.
    pub absolute_energy_residual_j: f64,
}

/// Complete physical split with the actual cost of its original RKL2 operators.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StabilizedMolecularStepReport {
    /// Original physical balances and reports for the accepted complete split.
    pub transport: MolecularStepReport,
    /// Includes both halves and every discarded complete-interval probe.
    pub operator_evaluations: usize,
    /// Largest stage count in the final accepted pair of molecular halves.
    pub maximum_stages_used: usize,
    /// Largest scaled parabolic number in the final accepted molecular halves.
    pub maximum_scaled_parabolic_number: f64,
}

/// Atomic nonreacting multicomponent molecular/advection evolution.
///
/// Prepared Maxwell--Stefan, Fourier and Newtonian transport surrounds the
/// original mixture Euler transaction. No reaction mechanism or frozen-composition
/// surrogate is introduced. Every shortened interval restarts both molecular
/// halves from the accepted field. Source refusal and cancellation retain all
/// accepted cells, time, reports and boundary transfers.
///
/// Symmetric splitting permits second-order smooth evolution when its constituent
/// paths retain their order; limiting and first-order retries can reduce it.
/// This is not an entropy, invariant-domain or experimental-validation certificate.
pub struct MolecularEuler2d<L: TransportLaw> {
    accepted: MixtureEuler2d,
    trial: MixtureEuler2d,
    molecular: MolecularTransport2d<L>,
    boundary_impulses: Vec<f64>,
    boundary_trial: Vec<f64>,
    maximum_workspace_bytes: usize,
    last_step: Option<MolecularStepReport>,
}
impl<L: TransportLaw> MolecularEuler2d<L> {
    /// Both field workspaces, prepared molecular buffers and complete boundary
    /// buffers, including incoming excess capacity and retained reservoir rows.
    /// Immutable sources, bounded stack scratch and provider storage are excluded.
    pub fn required_workspace_bytes(
        flow: &MixtureEuler2d,
        molecular: &MolecularTransport2d<L>,
    ) -> Result<usize, DissipationError> {
        flow.workspace_bytes()
            .checked_add(MixtureEuler2d::required_workspace_bytes(
                flow.gas(),
                flow.grid(),
            )?)
            .and_then(|bytes| bytes.checked_add(flow.reservoirs.workspace_bytes()))
            .and_then(|bytes| bytes.checked_add(molecular.workspace_bytes()))
            .and_then(|bytes| {
                bytes.checked_add(
                    flow.boundary_impulses
                        .len()
                        .checked_mul(2 * size_of::<f64>())?,
                )
            })
            .ok_or(DissipationError::Invalid(
                "coupled mixture workspace size overflow",
            ))
    }

    /// Complete coupled storage plus any unallocated RKL2 face/boundary scratch.
    /// Add supplied forcing workspace to this allowance before calling `new`.
    pub fn required_stabilized_workspace_bytes(
        flow: &MixtureEuler2d,
        molecular: &MolecularTransport2d<L>,
    ) -> Result<usize, DissipationError> {
        Self::required_workspace_bytes(flow, molecular)?
            .checked_add(molecular.additional_stabilized_workspace_bytes()?)
            .ok_or(DissipationError::Invalid(
                "coupled stabilized mixture workspace size overflow",
            ))
    }

    /// Bind the original field and its already-prepared molecular operator.
    /// The complete budget is also enforced on steps, including supplied forcing
    /// scratch even when its force is exactly zero. No mutable owner access is exposed.
    pub fn new(
        flow: MixtureEuler2d,
        molecular: MolecularTransport2d<L>,
        maximum_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DissipationError> {
        poll(&mut cancelled)?;
        if molecular.gas() != flow.gas()
            || molecular.grid != flow.grid
            || molecular.boundaries != flow.boundaries
        {
            return Err(DissipationError::Invalid(
                "coupled mixture molecular binding differs",
            ));
        }
        flow.reservoirs.admit_molecular(flow.boundaries)?;
        let required = Self::required_workspace_bytes(&flow, &molecular)?;
        admit_workspace(required, maximum_workspace_bytes)?;
        let mut states = allocated(flow.states.len(), 0.0)?;
        copy_values(&flow.states, &mut states, &mut cancelled)?;
        let mut trial = MixtureEuler2d::new(
            flow.gas.clone(),
            flow.grid,
            states,
            flow.boundaries,
            MixtureEuler2d::required_workspace_bytes(flow.gas(), flow.grid())?,
            &mut cancelled,
        )?;
        trial.reservoirs = flow.reservoirs.try_clone(&mut cancelled)?;
        trial.workspace_bytes += trial.reservoirs.workspace_bytes();
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

    /// Complete accepted field and clock; its Euler metadata describes the
    /// middle phase. This owner's report and face loads describe the full split.
    pub const fn flow(&self) -> &MixtureEuler2d {
        &self.accepted
    }
    /// Bound original operator; its published report describes the second half.
    pub const fn molecular(&self) -> &MolecularTransport2d<L> {
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
    pub const fn last_step(&self) -> Option<&MolecularStepReport> {
        self.last_step.as_ref()
    }
    /// Complete outward exchange on every physical face.
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
    /// Recover both accepted owners and their published metadata.
    pub fn into_parts(self) -> (MixtureEuler2d, MolecularTransport2d<L>) {
        (self.accepted, self.molecular)
    }

    /// Advance at most the requested duration with both molecular halves using
    /// the exact accepted central interval. Exhausted whole-step retries return
    /// a `DissipationError::StepFailure` with zero committed elapsed time.
    pub fn step(
        &mut self,
        maximum_duration: Time,
        controls: MolecularStepControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularStepReport, DissipationError> {
        self.step_impl(maximum_duration, controls, None, None, cancelled)
    }

    /// Compose held acceleration or force density with the same original
    /// molecular halves. Body work and species-referenced boundary enthalpy
    /// enter the complete physical-energy balance exactly once.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_forcing(
        &mut self,
        maximum_duration: Time,
        controls: MolecularStepControl,
        force: BodyForce2d<'_>,
        force_controls: ForceControl2d,
        workspace: &mut ForcingWorkspace2d,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularStepReport, DissipationError> {
        self.step_impl(
            maximum_duration,
            controls,
            Some((force, force_controls, workspace)),
            None,
            cancelled,
        )
    }

    /// Compose the original Maxwell--Stefan, enthalpy, heat and stress operators
    /// with RKL2 on both sides of the same conservative Euler transaction.
    /// One evaluation counter covers both halves and all discarded trials;
    /// exhaustion cannot be reset by shortening the coupled interval.
    ///
    /// Every stage retains its original physical/source admission and signed
    /// boundary quadrature. The macrostep cap remains an accuracy choice: the
    /// larger linear stability interval is not a nonlinear positivity or entropy
    /// guarantee. Field, clock, loads and reports publish only after full success.
    pub fn step_stabilized(
        &mut self,
        maximum_duration: Time,
        controls: StabilizedMolecularStepControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<StabilizedMolecularStepReport, DissipationError> {
        let mut work = StepWork::new(controls.molecular);
        let transport = self.step_impl(
            maximum_duration,
            controls.split(),
            None,
            Some(&mut work),
            cancelled,
        )?;
        Ok(StabilizedMolecularStepReport {
            transport,
            operator_evaluations: work.evaluations,
            maximum_stages_used: work.maximum_stages_used,
            maximum_scaled_parabolic_number: work.maximum_scaled_parabolic_number,
        })
    }

    /// Retain conservative body impulse/work and final Mach admission between
    /// the stabilized molecular halves. Even zero-force scratch counts toward
    /// the complete owner's workspace allowance.
    #[allow(clippy::too_many_arguments)]
    pub fn step_stabilized_with_forcing(
        &mut self,
        maximum_duration: Time,
        controls: StabilizedMolecularStepControl,
        force: BodyForce2d<'_>,
        force_controls: ForceControl2d,
        workspace: &mut ForcingWorkspace2d,
        cancelled: impl FnMut() -> bool,
    ) -> Result<StabilizedMolecularStepReport, DissipationError> {
        let mut work = StepWork::new(controls.molecular);
        let transport = self.step_impl(
            maximum_duration,
            controls.split(),
            Some((force, force_controls, workspace)),
            Some(&mut work),
            cancelled,
        )?;
        Ok(StabilizedMolecularStepReport {
            transport,
            operator_evaluations: work.evaluations,
            maximum_stages_used: work.maximum_stages_used,
            maximum_scaled_parabolic_number: work.maximum_scaled_parabolic_number,
        })
    }

    fn molecular_half(
        &mut self,
        duration: Time,
        controls: DissipationControl,
        stabilization: Option<&mut StepWork>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MolecularReport, DissipationError> {
        match stabilization {
            None => self.molecular.advance_unpublished(
                &mut self.trial,
                duration,
                controls,
                self.maximum_workspace_bytes,
                cancelled,
            ),
            Some(work) => {
                let report = self.molecular.advance_stabilized_unpublished(
                    &mut self.trial,
                    duration,
                    work.controls,
                    self.maximum_workspace_bytes,
                    &mut work.evaluations,
                    cancelled,
                )?;
                work.record(
                    report.maximum_stages_used,
                    report.maximum_scaled_parabolic_number,
                );
                Ok(report.transport)
            }
        }
    }

    fn step_impl(
        &mut self,
        maximum_duration: Time,
        controls: MolecularStepControl,
        mut forcing: Option<forcing::Context<'_, '_>>,
        mut stabilization: Option<&mut StepWork>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularStepReport, DissipationError> {
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
                    "coupled mixture workspace size overflow",
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
                "coupled mixture duration or controls",
            ));
        }
        if let Some(work) = stabilization.as_ref() {
            work.controls.validate()?;
            // The original operator checks its own field plus scratch. Reserve
            // the accepted owner, complete loads and supplied forcing separately
            // so lazy RKL2 storage cannot spend their portion of the same budget.
            let outside = bytes
                .checked_sub(self.trial.workspace_bytes())
                .and_then(|value| value.checked_sub(self.molecular.workspace_bytes()))
                .ok_or(DissipationError::Invalid(
                    "coupled stabilized mixture workspace size overflow",
                ))?;
            let allowance = self.maximum_workspace_bytes.checked_sub(outside).ok_or(
                DissipationError::Invalid("coupled stabilized mixture total workspace budget"),
            )?;
            self.molecular
                .reserve_stabilized_workspace(&self.trial, allowance, &mut cancelled)?;
        }
        let n = self.accepted.gas.species_count();
        let width = n + 3;
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
        // Avoid paying for a predictably excessive first molecular interval;
        // actual evolved faces and source stages still own central CFL admission.
        let mut dt = maximum_duration
            .value()
            .min(0.8 * controls.advection.cfl / rate);
        let mut last = DissipationError::Invalid("no coupled mixture attempt");
        let mut last_dt = dt;
        for retry in 0..=controls.maximum_coupling_retries {
            poll(&mut cancelled)?;
            if let Some(work) = stabilization.as_deref_mut() {
                work.reset_extrema();
            }
            if !positive(0.5 * dt)
                || !(self.accepted.time + dt).is_finite()
                || self.accepted.time + dt <= self.accepted.time
            {
                return Err(DissipationError::Invalid(
                    "unrepresentable coupled mixture increment",
                ));
            }
            copy_physical(&self.accepted, &mut self.trial, &mut cancelled)?;
            let mut shortened = None;
            let attempt = (|| -> Result<MolecularStepReport, DissipationError> {
                let first = self.molecular_half(
                    Time::new(0.5 * dt),
                    controls.molecular,
                    stabilization.as_deref_mut(),
                    &mut cancelled,
                )?;
                copy_values(
                    self.molecular.trial_impulses(),
                    &mut self.boundary_trial,
                    &mut cancelled,
                )?;
                // Source, control and scratch geometry are already admitted by
                // this owner and remain immutable across complete-step retries.
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
                        "advection shortened coupled mixture interval",
                    ));
                }
                add_impulses(
                    &self.trial.boundary_impulses,
                    &mut self.boundary_trial,
                    &mut cancelled,
                )?;
                let second = self.molecular_half(
                    Time::new(0.5 * dt),
                    controls.molecular,
                    stabilization.as_deref_mut(),
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
                    cartesian::check_cfl::<Error>(final_rate * dt, controls.advection.cfl)?;
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
                for axis in 0..2 {
                    residual[n + axis] -= source_change.impulse[axis];
                }
                residual[n + 2] -= source_change.work;
                if !residual[..width].iter().all(|value| value.is_finite()) {
                    return Err(DissipationError::Invalid(
                        "coupled mixture balance overflow",
                    ));
                }
                let conservation_residual = Totals::from_row(&residual[..width], n);
                let absolute_energy_residual_j =
                    conservation_residual.absolute_energy(&self.trial.gas)?;
                poll(&mut cancelled)?;
                Ok(MolecularStepReport {
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
                Err(error) if stabilization.is_some() && error.is_terminal() => {
                    return Err(error.step_failure(0.0, dt));
                }
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

fn admit_workspace(bytes: usize, maximum: usize) -> Result<(), DissipationError> {
    if bytes > maximum || bytes > MAX_MIXTURE_WORKSPACE_BYTES {
        return Err(DissipationError::Invalid(
            "coupled mixture total workspace budget",
        ));
    }
    Ok(())
}

fn copy_physical(
    source: &MixtureEuler2d,
    target: &mut MixtureEuler2d,
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
