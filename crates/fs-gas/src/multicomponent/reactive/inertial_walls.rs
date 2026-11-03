//! Finite-rate chemistry, resolved pressure waves and two finite inertial walls.

use super::{
    ChemistryReport, FlowError, ReactingEuler1d, ReactiveError, ReactiveStepControl,
    ReactiveStepReport, TransportStepReport, coupled_report, poll, positive, react_cells, sum,
};
use crate::dynamics::inertial_walls::{
    InertialWallLoad, InertialWallsControl, InertialWallsStep, validate_controls,
};
use crate::dynamics::piston::{Piston, PistonError};
use crate::multicomponent::{Boundaries, Boundary};
use fs_qty::Time;

/// One atomic chemistry/pressure/two-wall step on actual moving cell volumes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReactiveInertialWallsStep {
    /// Complete source and transport balances, including both chemistry halves.
    pub flow: ReactiveStepReport,
    /// Native simultaneous two-wall solve between the chemistry half steps.
    /// Its gas endpoints and energy balance cover that middle stage only.
    pub walls: InertialWallsStep<TransportStepReport>,
    /// Actual physical chemical-gas energy change [J], including endpoint
    /// species-reference changes rather than only stored gauged energy.
    pub gas_energy_change_j: f64,
    /// Complete gas plus BOTH walls momentum defect [N s], subtracting signed
    /// exterior-pressure and additional-force impulses. No fixed support exists.
    pub combined_momentum_residual_n_s: f64,
    /// Complete chemical-gas plus wall kinetic-energy defect [J], including
    /// exterior-pressure work and subtracting additional-force work.
    pub total_energy_residual_j: f64,
}

impl ReactingEuler1d {
    /// Advance the original chemistry and both finite gas end walls together.
    ///
    /// The composition is C(dt/2), the native simultaneous two-wall ALE step
    /// P(dt), then C(dt/2). The first chemistry half uses the accepted physical
    /// cell volumes; the second uses their actual moved values. The existing
    /// reactor and species-reference correction supply chemical conversion once,
    /// without adding a second heat-release source or an artificial support.
    ///
    /// Both midpoint-speed brackets constrain the acoustic and contraction
    /// timestep. If the native pressure solve shortens that duration, the entire
    /// chemistry/motion composition restarts from accepted gas and BOTH accepted
    /// wall velocities. Source and final-balance refusals also restart the whole
    /// attempt. Existing field, geometry, source and amount buffers are reused.
    ///
    /// Wall order is [left, right], with positive velocities and additional
    /// forces along the increasing axis. The original wall tolerances admit
    /// complete gas-plus-walls momentum and chemical-energy balances as well as
    /// the middle mechanical stage; choose them consistently with the original
    /// per-cell chemistry accuracy. Publication of species, geometry, clock and
    /// both wall velocities follows the final cancellation poll.
    ///
    /// This closed constant-area composition does not advance separately bound
    /// molecular or thermal owners. The original chemistry/splitting and
    /// unequal-caloric-contact limitations remain; symmetric splitting alone is
    /// not a second-order, entropy-stability or combustion-validity theorem.
    ///
    /// # Errors
    /// Invalid inputs, source/geometry/mechanical refusal, exhausted work or
    /// cancellation leave every accepted gas and wall state unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_inertial_walls(
        &mut self,
        walls: &mut [Piston; 2],
        maximum_duration: Time,
        controls: ReactiveStepControl,
        loads: [InertialWallLoad; 2],
        wall_controls: InertialWallsControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveInertialWallsStep, ReactiveError> {
        poll(&mut cancelled)?;
        validate_controls(
            *walls,
            maximum_duration,
            controls.transport,
            loads,
            wall_controls,
        )?;
        if controls.maximum_coupling_retries > 32
            || self.accepted.boundaries
                != (Boundaries::Ends {
                    left: Boundary::Reflecting,
                    right: Boundary::Reflecting,
                })
            || self.accepted.mesh.uniform_cross_section().is_none()
        {
            return Err(ReactiveError::Invalid(
                "reacting inertial walls require bounded retries and a constant-area closed interval",
            ));
        }
        let before = self.accepted.totals();
        let area = self.accepted.mesh.face_areas()[0];
        let bounds = wall_controls.midpoint_velocity_bracket_m_s;
        let low = bounds[0][0].min(bounds[1][0]);
        let high = bounds[0][1].max(bounds[1][1]);
        let mut duration = maximum_duration.value();
        for i in 0..self.accepted.cell_count() {
            let primitive = self.accepted.gas().primitive(
                self.accepted.cell(i).expect("admitted cell index"),
                &mut cancelled,
            )?;
            let relative = (primitive.velocity - low)
                .abs()
                .max((primitive.velocity - high).abs());
            let rate =
                (relative + primitive.sound_speed) * (area / self.accepted.mesh.volumes()[i]);
            if !positive(rate) {
                return Err(ReactiveError::Invalid(
                    "reacting inertial walls relative-wave CFL overflow",
                ));
            }
            duration = duration.min(0.8 * controls.transport.cfl / rate);
        }
        let contraction = bounds[0][1] - bounds[1][0];
        if !contraction.is_finite() {
            return Err(ReactiveError::Invalid(
                "reacting inertial walls contraction bound overflow",
            ));
        }
        if contraction > 0.0 {
            let length = self.accepted.mesh.positions()[self.accepted.cell_count()]
                - self.accepted.mesh.positions()[0];
            duration = duration.min(0.8 * controls.transport.cfl * (length / contraction));
        }
        let original_walls = *walls;
        let mut last = ReactiveError::Invalid("no admissible reacting two-wall increment");
        for retry in 0..=controls.maximum_coupling_retries {
            poll(&mut cancelled)?;
            if !positive(duration)
                || !positive(0.5 * duration)
                || !positive(self.accepted.time + duration - self.accepted.time)
                || !(self.accepted.time + duration).is_finite()
            {
                return Err(ReactiveError::Invalid(
                    "reacting inertial walls cannot advance the clock",
                ));
            }
            self.trial
                .copy_physical_with_motion_from(&self.accepted, &mut cancelled)?;
            let mut trial_walls = original_walls;
            let mut chemistry = ChemistryReport::default();
            let attempt = (|| {
                react_cells(
                    &mut self.trial,
                    &mut self.reactor,
                    &mut self.amounts,
                    0.5 * duration,
                    &mut chemistry,
                    &mut cancelled,
                )?;
                let mechanical = self.trial.step_with_inertial_walls(
                    &mut trial_walls,
                    Time::new(duration),
                    controls.transport,
                    loads,
                    wall_controls,
                    &mut cancelled,
                )?;
                if mechanical.gas.elapsed.value() != duration {
                    return Ok((mechanical.gas.elapsed.value(), None));
                }
                react_cells(
                    &mut self.trial,
                    &mut self.reactor,
                    &mut self.amounts,
                    0.5 * duration,
                    &mut chemistry,
                    &mut cancelled,
                )?;
                let flow = coupled_report(
                    &self.trial.gas,
                    before,
                    self.trial.totals(),
                    mechanical.gas,
                    chemistry,
                    None,
                    retry,
                )?;
                // Restore the actual pressure boundary work to the original
                // endpoint-based chemical energy residual. The species datum
                // is already counted once by coupled_report.
                let gas_energy_change_j =
                    flow.absolute_energy_residual_j + flow.transport.boundary_change.energy;
                let total_energy_residual_j =
                    sum(core::iter::once(gas_energy_change_j).chain((0..2).map(|i| {
                        mechanical.kinetic_energy_change_j[i]
                            + mechanical.exterior_pressure_work_j[i]
                            - mechanical.external_force_work_j[i]
                    })));
                let combined_momentum_residual_n_s = sum(core::iter::once(
                    flow.after.momentum - flow.before.momentum,
                )
                .chain((0..2).map(|i| {
                    let external_impulse = duration
                        * (loads[i].external_force_n
                            + [1.0, -1.0][i] * area * loads[i].exterior_pressure_pa);
                    original_walls[i].mass_kg()
                        * (trial_walls[i].velocity_m_s() - original_walls[i].velocity_m_s())
                        - external_impulse
                })));
                if ![
                    gas_energy_change_j,
                    total_energy_residual_j,
                    combined_momentum_residual_n_s,
                ]
                .iter()
                .all(|x| x.is_finite())
                {
                    return Err(ReactiveError::Invalid(
                        "non-finite complete chemical/two-wall balance",
                    ));
                }
                if total_energy_residual_j.abs() > wall_controls.energy_tolerance_j
                    || combined_momentum_residual_n_s.abs() > wall_controls.momentum_tolerance_n_s
                {
                    return Err(ReactiveError::Piston(PistonError::Balance {
                        momentum_residual_n_s: combined_momentum_residual_n_s.abs(),
                        energy_residual_j: total_energy_residual_j.abs(),
                    }));
                }
                Ok::<_, ReactiveError>((
                    duration,
                    Some(ReactiveInertialWallsStep {
                        flow,
                        walls: mechanical,
                        gas_energy_change_j,
                        combined_momentum_residual_n_s,
                        total_energy_residual_j,
                    }),
                ))
            })();
            match attempt {
                Ok((shorter, None)) => {
                    duration = shorter;
                    last = ReactiveError::Invalid(
                        "inertial walls shortened the complete reacting duration",
                    );
                }
                Ok((_, Some(report))) => {
                    poll(&mut cancelled)?;
                    core::mem::swap(&mut self.accepted, &mut self.trial);
                    *walls = trial_walls;
                    return Ok(report);
                }
                Err(ReactiveError::Cancelled) => return Err(ReactiveError::Cancelled),
                Err(error @ ReactiveError::Flow(FlowError::Allocation))
                | Err(
                    error @ ReactiveError::Piston(PistonError::MixtureGas(FlowError::Allocation)),
                ) => {
                    return Err(error);
                }
                Err(error) => {
                    last = error;
                    duration *= 0.5;
                }
            }
        }
        Err(ReactiveError::RetryLimit {
            attempts: controls.maximum_coupling_retries + 1,
            last: Box::new(last),
        })
    }
}
