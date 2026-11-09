//! Original reacting sources composed with conservative moving volumes.

use super::*;
use crate::compressible::AxisBoundaries;
use crate::compressible::spatial::mixture::molecular::{DissipationError, WallTransfer};
use crate::compressible::spatial::{BoundaryCondition, MeshMotion};

/// Complete reacting ALE interval, with directly integrated molecular wall heat
/// and work when the original molecular operator is explicitly supplied.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MovingReactiveStepReport3d {
    /// Original species, all momenta, chemistry and physical-energy balances.
    pub reactive: ReactiveStepReport3d,
    /// Molecular stress work into the gas [J], ordered XMin, XMax, YMin, YMax,
    /// ZMin, ZMax. Pressure work remains in `reactive.transport`.
    pub molecular_wall_work_j: [f64; 6],
    /// Fourier and Navier-partition heat into the gas [J], in the same order.
    /// These terms are evaluated before adding mechanical work.
    pub wall_heat_j: [f64; 6],
}

pub(super) struct MovingTrialFailure {
    pub(super) error: ReactiveError,
    pub(super) shortened: Option<f64>,
}
impl From<ReactiveError> for MovingTrialFailure {
    fn from(error: ReactiveError) -> Self {
        Self {
            error,
            shortened: None,
        }
    }
}
impl From<super::super::Error> for MovingTrialFailure {
    fn from(error: super::super::Error) -> Self {
        Self::from(ReactiveError::from(error))
    }
}
impl From<DissipationError> for MovingTrialFailure {
    fn from(error: DissipationError) -> Self {
        Self::from(ReactiveError::from(error))
    }
}

impl ReactingEuler3d {
    /// Advance C(dt/2), ALE(dt), C(dt/2) on a translating or single-axis
    /// deforming closed or periodic Cartesian domain.
    ///
    /// The original prepared cell chemistry sees each phase's actual volume,
    /// keeps all three momenta fixed and accounts for the original species
    /// energy references. Native moving-face fluxes evolve extensive species,
    /// momentum and energy with the discrete geometric conservation law.
    /// Moving reflecting walls exchange pressure impulse/work and no species.
    ///
    /// Final source temperature and relative acoustic limits are checked after
    /// the last chemistry half. A shortened wave phase restarts the entire
    /// physical interval; cells, geometry, clock, reports and loads publish
    /// together. The selected chemistry time method retains its own accuracy;
    /// symmetric splitting alone is not a second-order or entropy guarantee.
    pub fn step_moving(
        &mut self,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        motion: MeshMotion,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MovingReactiveStepReport3d, ReactiveError> {
        self.step_moving_impl::<ConstantTransport>(
            maximum_duration,
            controls,
            None,
            DissipationControl::default(),
            super::super::MAX_MIXTURE_WORKSPACE_BYTES,
            motion,
            cancelled,
        )
    }

    /// Add native molecular half intervals around C/ALE/C, using the old and
    /// new grid respectively. Original Maxwell--Stefan/enthalpy transport,
    /// Fourier heat and full moving-wall stress/work remain source-owned.
    ///
    /// The explicit byte budget covers both reacting frames, chemistry and
    /// molecular buffers. Shortening or refusal restarts D/C/ALE/C/D together.
    /// No source binding, accepted molecular report or boundary load changes
    /// until final publication. Open ports, retained SGS, tangential finite
    /// wall inertia and simultaneous multi-axis deformation are not composed.
    #[allow(clippy::too_many_arguments)]
    pub fn step_moving_with_dissipation<L: TransportLaw>(
        &mut self,
        operator: &mut MolecularTransport3d<L>,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        dissipation_controls: DissipationControl,
        maximum_workspace_bytes: usize,
        motion: MeshMotion,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MovingReactiveStepReport3d, ReactiveError> {
        self.step_moving_impl(
            maximum_duration,
            controls,
            Some(operator),
            dissipation_controls,
            maximum_workspace_bytes,
            motion,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_moving_impl<L: TransportLaw>(
        &mut self,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        mut operator: Option<&mut MolecularTransport3d<L>>,
        dissipation_controls: DissipationControl,
        maximum_workspace_bytes: usize,
        motion: MeshMotion,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MovingReactiveStepReport3d, ReactiveError> {
        poll(&mut cancelled)?;
        self.admit_moving(
            maximum_duration.value(),
            controls,
            operator.as_deref(),
            dissipation_controls,
            maximum_workspace_bytes,
            motion,
        )?;
        let before = checked_totals(
            &self.accepted.states,
            &self.accepted.grid,
            self.accepted.gas.species_count(),
            &mut cancelled,
        )?;
        let rate = endpoint_rate(
            &self.accepted,
            motion,
            self.reactor.controls(),
            &mut cancelled,
        )?;
        let mut dt = maximum_duration
            .value()
            .min(0.8 * controls.transport.cfl / rate);
        for retry in 0..=controls.maximum_coupling_retries {
            match self.trial_moving(
                dt,
                controls,
                operator.as_deref_mut(),
                dissipation_controls,
                maximum_workspace_bytes,
                motion,
                before,
                retry,
                &mut cancelled,
            ) {
                Ok(report) => {
                    self.publish_moving(report, operator.as_deref_mut());
                    return Ok(report);
                }
                Err(MovingTrialFailure {
                    error: ReactiveError::Cancelled,
                    ..
                }) => {
                    return Err(ReactiveError::Cancelled);
                }
                Err(failure) => {
                    if retry == controls.maximum_coupling_retries {
                        return Err(ReactiveError::RetryLimit {
                            attempts: retry + 1,
                            last: Box::new(failure.error),
                        });
                    }
                    dt = failure.shortened.unwrap_or(0.5 * dt);
                }
            }
        }
        unreachable!("bounded moving reactive retries return on the final attempt")
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn admit_moving<L: TransportLaw>(
        &self,
        maximum_duration: f64,
        controls: ReactiveStepControl,
        operator: Option<&MolecularTransport3d<L>>,
        dissipation_controls: DissipationControl,
        maximum_workspace_bytes: usize,
        motion: MeshMotion,
    ) -> Result<(), ReactiveError> {
        if !positive(maximum_duration)
            || controls.maximum_coupling_retries > 32
            || !positive(controls.transport.cfl)
            || controls.transport.cfl > 0.45
            || controls.transport.maximum_halvings > 32
        {
            return Err(ReactiveError::Invalid(
                "moving reactive duration, CFL or retry controls",
            ));
        }
        let bytes = self
            .workspace_bytes()
            .checked_add(operator.map_or(0, |op| op.workspace_bytes()))
            .ok_or(ReactiveError::Invalid(
                "moving reactive workspace size overflow",
            ))?;
        if bytes > maximum_workspace_bytes || bytes > super::super::MAX_MIXTURE_WORKSPACE_BYTES {
            return Err(ReactiveError::Invalid(
                "moving reactive total workspace budget",
            ));
        }
        self.accepted.admit_motion(motion)?;
        for axis in 0..3 {
            if let AxisBoundaries::Ends { lower, upper } = self.accepted.boundaries.axis(axis)
                && [lower, upper]
                    .into_iter()
                    .any(|bc| bc != BoundaryCondition::Reflecting)
            {
                return Err(ReactiveError::Invalid(
                    "moving reactive volumes require closed or periodic boundaries",
                ));
            }
        }
        if let Some(operator) = operator {
            dissipation_controls.validate()?;
            operator.admit_moving_field(&self.accepted, motion)?;
        }
        Ok(())
    }

    /// Unpublished fixed-duration proposal used by both explicit motion and
    /// finite-body coupling. All source halves restart from accepted storage.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn trial_moving<L: TransportLaw>(
        &mut self,
        dt: f64,
        controls: ReactiveStepControl,
        mut operator: Option<&mut MolecularTransport3d<L>>,
        dissipation_controls: DissipationControl,
        maximum_workspace_bytes: usize,
        motion: MeshMotion,
        before: Totals,
        retry: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MovingReactiveStepReport3d, MovingTrialFailure> {
        poll(cancelled)?;
        if !positive(0.5 * dt)
            || !(self.accepted.time + dt).is_finite()
            || self.accepted.time + dt <= self.accepted.time
        {
            return Err(ReactiveError::Invalid("unrepresentable moving reactive increment").into());
        }
        if motion != MeshMotion::default() {
            motion
                .grid_after(&self.accepted.grid, dt)
                .map_err(|error| {
                    ReactiveError::Invalid(match error {
                        fs_gas::dynamics::Error::Invalid(reason) => reason,
                        _ => "moving reactive geometry",
                    })
                })?;
        }
        copy_physical(&self.accepted, &mut self.trial, cancelled)?;
        let mut chemistry = ChemistryAccumulator::default();
        let mut transfer = WallTransfer::default();
        let first = if let Some(operator) = operator.as_deref_mut() {
            let (report, exchange) = operator.advance_moving_unpublished(
                &mut self.trial,
                Time::new(0.5 * dt),
                dissipation_controls,
                maximum_workspace_bytes,
                motion,
                cancelled,
            )?;
            copy_field(
                operator.trial_impulses(),
                &mut self.boundary_trial,
                cancelled,
            )?;
            transfer = exchange;
            Some(report)
        } else {
            None
        };
        react_cells(
            &mut self.trial,
            &mut self.reactor,
            0.5 * dt,
            &mut chemistry,
            cancelled,
        )?;
        let transport =
            self.trial
                .step_moving(Time::new(dt), controls.transport, motion, &mut *cancelled)?;
        if transport.elapsed.value() != dt {
            return Err(MovingTrialFailure {
                error: ReactiveError::Invalid("ALE shortened moving reactive interval"),
                shortened: Some(transport.elapsed.value()),
            });
        }
        if first.is_some() {
            add_impulses(
                &self.trial.boundary_impulses,
                &mut self.boundary_trial,
                cancelled,
            )?;
        } else {
            copy_field(
                &self.trial.boundary_impulses,
                &mut self.boundary_trial,
                cancelled,
            )?;
        }
        react_cells(
            &mut self.trial,
            &mut self.reactor,
            0.5 * dt,
            &mut chemistry,
            cancelled,
        )?;
        let second = if let Some(operator) = operator.as_deref_mut() {
            let (report, exchange) = operator.advance_moving_unpublished(
                &mut self.trial,
                Time::new(0.5 * dt),
                dissipation_controls,
                maximum_workspace_bytes,
                motion,
                cancelled,
            )?;
            add_impulses(
                operator.trial_impulses(),
                &mut self.boundary_trial,
                cancelled,
            )?;
            transfer.add_scaled(exchange, 1.0, &mut WallTransfer::default())?;
            Some(report)
        } else {
            None
        };
        // Both reacting endpoints and viscously heated endpoints remain in
        // the original chemical interval and the moving-face acoustic budget.
        let rate = endpoint_rate(&self.trial, motion, self.reactor.controls(), cancelled)?;
        super::super::stepping::check_cfl::<super::super::Error>(
            rate * dt,
            controls.transport.cfl,
        )?;
        let n = self.trial.gas.species_count();
        let width = n + 4;
        let after = checked_totals(&self.trial.states, &self.trial.grid, n, cancelled)?;
        let outward = sum_rows(&self.boundary_trial, width, 1.0, cancelled)?;
        let boundary_change = Totals::from_row(&outward.map(|v| -v), n);
        let source_change = transport.source_change;
        let chemistry = chemistry.report;
        let mut delta = [0.0; MAX_SPECIES + 4];
        for (d, value) in delta[..width].iter_mut().enumerate() {
            *value = after.component(d) - before.component(d) - boundary_change.component(d);
        }
        for axis in 0..3 {
            delta[n + axis] -= source_change.impulse[axis];
        }
        delta[n + 3] -= source_change.work;
        let mass_residual_kg = sum(delta[..n].iter().copied());
        // Chemical energy is conservative only after its original species
        // references are restored. Form small conservative differences first.
        let absolute_energy_residual_j = sum(std::iter::once(delta[n + 3])
            .chain((0..n).map(|k| delta[k] * self.trial.gas.energy_references()[k])));
        for (k, value) in delta[..n].iter_mut().enumerate() {
            *value -= chemistry.species_change_kg[k];
        }
        delta[n + 3] -= chemistry.referenced_energy_change_j;
        if !delta[..width].iter().all(|v| v.is_finite())
            || !mass_residual_kg.is_finite()
            || !absolute_energy_residual_j.is_finite()
            || !chemistry.entropy_change_j_k.is_finite()
        {
            return Err(ReactiveError::Invalid("moving reactive balance overflow").into());
        }
        let report = MovingReactiveStepReport3d {
            reactive: ReactiveStepReport3d {
                elapsed: Time::new(dt),
                coupling_retries: retry,
                transport,
                chemistry,
                dissipation: first
                    .zip(second)
                    .map(|(first, second)| SplitDissipationReport3d { first, second }),
                before,
                after,
                boundary_change,
                source_change,
                conservation_residual: Totals::from_row(&delta[..width], n),
                mass_residual_kg,
                absolute_energy_residual_j,
            },
            molecular_wall_work_j: transfer.mechanical,
            wall_heat_j: transfer.heat,
        };
        poll(cancelled)?;
        Ok(report)
    }

    pub(super) fn publish_moving<L: TransportLaw>(
        &mut self,
        report: MovingReactiveStepReport3d,
        operator: Option<&mut MolecularTransport3d<L>>,
    ) {
        std::mem::swap(&mut self.accepted, &mut self.trial);
        std::mem::swap(&mut self.boundary_impulses, &mut self.boundary_trial);
        if let (Some(operator), Some(dissipation)) = (operator, report.reactive.dissipation) {
            operator.publish_moving(self.accepted.grid, dissipation.second);
        }
        self.last_step = Some(report.reactive);
    }
}

fn endpoint_rate(
    flow: &MixtureEuler3d,
    motion: MeshMotion,
    chemistry: IsochoricKineticsSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, ReactiveError> {
    let n = flow.gas.species_count();
    let mut rate: f64 = 0.0;
    for (index, row) in flow.states.chunks_exact(n + 4).enumerate() {
        let p = cell(row, n).primitive(&flow.gas, &mut *cancelled)?;
        if !(chemistry.minimum_temperature_k..=chemistry.maximum_temperature_k)
            .contains(&p.temperature)
        {
            return Err(ReactiveError::Invalid(
                "moving endpoint is outside chemistry interval",
            ));
        }
        let coordinates = flow.grid.coordinates(index);
        let local = (0..3)
            .map(|axis| {
                let speed = if motion == MeshMotion::default() {
                    p.velocity[axis].abs()
                } else {
                    let lo = motion.velocity(&flow.grid, axis, coordinates[axis]);
                    let hi = motion.velocity(&flow.grid, axis, coordinates[axis] + 1);
                    (p.velocity[axis] - lo)
                        .abs()
                        .max((p.velocity[axis] - hi).abs())
                };
                (speed + p.sound_speed) / flow.grid.spacing[axis]
            })
            .sum::<f64>();
        if !local.is_finite() {
            return Err(ReactiveError::Invalid(
                "moving reactive relative CFL overflow",
            ));
        }
        rate = rate.max(local);
    }
    Ok(rate)
}
