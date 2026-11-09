//! Reciprocal finite guided pistons driven by original reacting spatial gas.

use super::super::WIDTH_LIMIT;
use super::super::molecular::MolecularWall3d;
use super::*;
use crate::compressible::spatial::MeshMotion;
use crate::compressible::spatial::inertial_walls::wall_torque;
pub use crate::compressible::spatial::molecular::{
    InertialWallLoad, InertialWallsControl, Piston, PistonError,
};
use crate::compressible::spatial::molecular::{admit_planes, bounded_duration, displacement};
use fs_gas::dynamics::inertial_walls::{WallGasExchange, solve_coupled};

mod body_forcing;

/// Original reacting source or mechanical refusal, with atomic publication.
pub type ReactiveInertialWallsError =
    crate::compressible::spatial::molecular::MolecularInertialWallsError<ReactiveError>;

/// Complete reacting gas, two finite pistons and explicit support/guide loads.
/// The mechanical energy gate uses actual chemical-plus-kinetic gas energy;
/// chemistry's referenced-energy change is also admitted independently.
pub type ReactiveInertialWallsStep =
    crate::compressible::spatial::molecular::MolecularInertialWallsStep<MovingReactiveStepReport3d>;

impl From<ReactiveError> for ReactiveInertialWallsError {
    fn from(error: ReactiveError) -> Self {
        if error == ReactiveError::Cancelled {
            Self::Cancelled
        } else {
            Self::Transport(error)
        }
    }
}

impl ReactingEuler3d {
    /// Advance the original reacting species field and two finite guided pistons.
    ///
    /// Every midpoint-speed proposal replays C(dt/2), ALE(dt), C(dt/2) from
    /// the same accepted cells. The chemistry halves use their actual old/new
    /// cell volumes and keep all three cell momenta fixed. Reaction pressure
    /// drives the original nonlinear pressure impulse and conjugate piston
    /// work; no separate combustion-force or chemical-heating term is added.
    ///
    /// `axis` is x=0, y=1, z=2. The selected minimum/maximum ends reflect;
    /// other directions are periodic or stationary reflecting supports. Both
    /// pistons retain their original finite mass and receive the negative of
    /// the accepted gas-side impulse. Ideal guides carry the reported angular
    /// impulse while retaining each plane's orientation. This is the original
    /// bounded simultaneous mechanical solve, not a new FSI algorithm.
    ///
    /// Source, caloric, relative-wave and chemistry controls apply to every
    /// proposal. Refusals restart the complete interval. The physical energy
    /// gate includes actual species-reference changes; a second gate accounts
    /// for the source's referenced-energy change exactly once. Cells, geometry,
    /// clock, reports, loads and both pistons publish after one final poll.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_inertial_walls(
        &mut self,
        axis: usize,
        walls: &mut [Piston; 2],
        maximum_duration: Time,
        transport_controls: ReactiveStepControl,
        loads: [InertialWallLoad; 2],
        controls: InertialWallsControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveInertialWallsStep, ReactiveInertialWallsError> {
        self.piston_step::<ConstantTransport>(
            None,
            axis,
            walls,
            maximum_duration,
            transport_controls,
            DissipationControl::default(),
            super::super::MAX_MIXTURE_WORKSPACE_BYTES,
            None,
            loads,
            controls,
            cancelled,
        )
    }

    /// Compose original molecular transport, reacting gas and finite pistons.
    ///
    /// Each fixed-duration proposal uses D/2--C/2--ALE--C/2--D/2 on the actual
    /// original and new geometry. Maxwell--Stefan species transport, partial
    /// enthalpies, Newtonian stress, Fourier heat and the source reactor keep
    /// their original owners. Every piston receives pressure plus molecular
    /// traction with conjugate mechanical work; directly integrated Fourier
    /// and Navier-partition heat remains separate external thermal exchange.
    ///
    /// All authored tangential wall velocities must be zero. Fixed transverse
    /// walls may absorb axial viscous momentum; the complete side-support
    /// impulse is explicit in the combined momentum gate. Translation and
    /// orientation guides carry the actual transverse and angular reactions
    /// without mechanical work. Source, molecular and acoustic refusals or a
    /// shorter ALE interval restart the complete chemistry/mechanics trial.
    /// The byte budget covers both prepared reacting and molecular owners.
    /// Neither owner nor either piston publishes before the final checkpoint.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_dissipation_and_inertial_walls<L: TransportLaw>(
        &mut self,
        operator: &mut MolecularTransport3d<L>,
        axis: usize,
        walls: &mut [Piston; 2],
        maximum_duration: Time,
        transport_controls: ReactiveStepControl,
        dissipation_controls: DissipationControl,
        maximum_workspace_bytes: usize,
        loads: [InertialWallLoad; 2],
        controls: InertialWallsControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveInertialWallsStep, ReactiveInertialWallsError> {
        self.piston_step(
            Some(operator),
            axis,
            walls,
            maximum_duration,
            transport_controls,
            dissipation_controls,
            maximum_workspace_bytes,
            None,
            loads,
            controls,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn piston_step<L: TransportLaw>(
        &mut self,
        mut operator: Option<&mut MolecularTransport3d<L>>,
        axis: usize,
        walls: &mut [Piston; 2],
        maximum_duration: Time,
        transport_controls: ReactiveStepControl,
        dissipation_controls: DissipationControl,
        maximum_workspace_bytes: usize,
        mut forcing: Option<forcing::Context<'_, '_>>,
        loads: [InertialWallLoad; 2],
        controls: InertialWallsControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveInertialWallsStep, ReactiveInertialWallsError> {
        poll(&mut cancelled)?;
        controls.validate(
            *walls,
            maximum_duration,
            transport_controls.transport,
            loads,
        )?;
        let remaining_bytes =
            self.prepare_moving_forcing(&mut forcing, maximum_workspace_bytes, &mut cancelled)?;
        self.admit_moving(
            maximum_duration.value(),
            transport_controls,
            operator.as_deref(),
            dissipation_controls,
            remaining_bytes,
            MeshMotion::default(),
        )?;
        let grid = self.accepted.grid;
        let wall_laws = operator
            .as_deref()
            .map_or([MolecularWall3d::SlipAdiabatic; 6], |op| {
                BoundarySide::ALL.map(|side| op.wall(side))
            });
        let sides = admit_planes(grid, self.accepted.boundaries, wall_laws, axis)?;
        let area = grid.face_area(sides[0]) * grid.face_count(sides[0]) as f64;
        let bounds = controls.midpoint_velocity_bracket_m_s;
        let low = bounds[0][0].min(bounds[1][0]);
        let high = bounds[0][1].max(bounds[1][1]);
        let n = self.accepted.gas.species_count();
        let width = n + 4;
        let mut rate = 0.0_f64;
        for row in self.accepted.states.chunks_exact(width) {
            let p = cell(row, n)
                .primitive(&self.accepted.gas, &mut cancelled)
                .map_err(ReactiveError::from)?;
            if let Some((_, control, _)) = &forcing {
                control
                    .admit_velocity::<super::super::Error>(p.velocity, p.sound_speed)
                    .map_err(ReactiveError::from)?;
            }
            let local: f64 = (0..3)
                .map(|a| {
                    let relative = if a == axis {
                        (p.velocity[a] - low)
                            .abs()
                            .max((p.velocity[a] - high).abs())
                    } else {
                        p.velocity[a].abs()
                    };
                    (relative + p.sound_speed) / grid.spacing[a]
                })
                .sum();
            if !positive(local) {
                return Err(PistonError::Invalid("reactive inertial-plane CFL overflow").into());
            }
            rate = rate.max(local);
        }
        let mut dt = bounded_duration(
            maximum_duration.value(),
            transport_controls.transport.cfl,
            rate,
            grid,
            axis,
            bounds,
        )?;
        let before = checked_totals(&self.accepted.states, &grid, n, &mut cancelled)
            .map_err(ReactiveError::from)?;
        let original = *walls;
        let mut gas_trials = 0;
        let mut side_impulse = [0.0; 3];
        let mut inward = [[0.0; WIDTH_LIMIT]; 2];
        for retry in 0..=transport_controls.maximum_coupling_retries {
            poll(&mut cancelled)?;
            let mut shortened = None;
            let result = solve_coupled(
                original,
                Time::new(dt),
                transport_controls.transport,
                loads,
                controls,
                area,
                |speed| -> Result<
                    WallGasExchange<MovingReactiveStepReport3d>,
                    ReactiveInertialWallsError,
                > {
                    gas_trials += 1;
                    let mut motion = MeshMotion::default();
                    motion.lower_velocity[axis] = speed[0];
                    motion.upper_velocity[axis] = speed[1];
                    let gas = self
                        .trial_moving(
                            dt,
                            transport_controls,
                            operator.as_deref_mut(),
                            dissipation_controls,
                            maximum_workspace_bytes,
                            motion,
                            before,
                            retry,
                            forcing
                                .as_mut()
                                .map(|(force, control, work)| (*force, *control, &mut **work)),
                            &mut cancelled,
                        )
                        .map_err(|failure| {
                            shortened = failure.shortened;
                            ReactiveInertialWallsError::from(failure.error)
                        })?;
                    let mut work = [0.0; 2];
                    for end in 0..2 {
                        let faces = grid.boundary_range(sides[end]);
                        let values = faces.start * width..faces.end * width;
                        inward[end] = sum_rows(
                            &self.boundary_trial[values.clone()],
                            width,
                            -1.0,
                            &mut cancelled,
                        )
                        .map_err(ReactiveError::from)?;
                        // The middle ALE owner's face account contains the
                        // pressure work. Chemistry contributes no face flux.
                        let pressure = sum_rows(
                            &self.trial.boundary_impulses[values],
                            width,
                            -1.0,
                            &mut cancelled,
                        )
                        .map_err(ReactiveError::from)?[n + 3];
                        work[end] = pressure + gas.molecular_wall_work_j[2 * axis + end];
                    }
                    side_impulse.fill(0.0);
                    let mut correction = [0.0; 3];
                    for side in BoundarySide::ALL
                        .into_iter()
                        .filter(|side| side.axis() != axis)
                    {
                        let faces = grid.boundary_range(side);
                        let values = sum_rows(
                            &self.boundary_trial[faces.start * width..faces.end * width],
                            width,
                            -1.0,
                            &mut cancelled,
                        )
                        .map_err(ReactiveError::from)?;
                        for a in 0..3 {
                            add(&mut side_impulse[a], &mut correction[a], values[n + a]);
                        }
                    }
                    let after = gas.reactive.after;
                    // Chemical conversion changes the stored reference chart,
                    // but is internal to the physical gas+piston energy. Use
                    // actual extensive species changes before multiplying each
                    // source reference, never large absolute energy totals.
                    let physical_change = sum(std::iter::once(after.energy - before.energy)
                        .chain((0..n).map(|k| {
                            (after.species_masses[k] - before.species_masses[k])
                                * self.accepted.gas.energy_references()[k]
                        }))
                        .chain(gas.wall_heat_j.map(|q| -q))
                        .chain(
                            (gas.reactive.source_change.work != 0.0)
                                .then_some(-gas.reactive.source_change.work),
                        ));
                    let travel = displacement(grid, self.trial.grid, axis, speed)?;
                    Ok(WallGasExchange {
                        gas,
                        momentum_change_n_s: after.momentum[axis]
                            - before.momentum[axis]
                            - side_impulse[axis]
                            - gas.reactive.source_change.impulse[axis],
                        energy_change_j: physical_change,
                        boundary_impulse_n_s: inward.map(|q| q[n + axis]),
                        boundary_work_j: work,
                        displacement_m: travel,
                    })
                },
                terminal,
            )
            .and_then(|mechanics| {
                // Independently close the stored energy equation. Its chemical
                // source appears exactly once; it is not external wall heat.
                let reactive = &mechanics.gas.reactive;
                let referenced = sum([
                    reactive.after.energy - before.energy,
                    -reactive.chemistry.referenced_energy_change_j,
                ]
                .into_iter()
                .chain(mechanics.gas.wall_heat_j.map(|q| -q))
                .chain((reactive.source_change.work != 0.0).then_some(-reactive.source_change.work))
                .chain((0..2).map(|i| {
                    mechanics.kinetic_energy_change_j[i] + mechanics.exterior_pressure_work_j[i]
                        - mechanics.external_force_work_j[i]
                })));
                if !referenced.is_finite() || referenced.abs() > controls.energy_tolerance_j {
                    return Err(PistonError::Balance {
                        momentum_residual_n_s: mechanics.combined_momentum_residual_n_s.abs(),
                        energy_residual_j: referenced.abs(),
                    }
                    .into());
                }
                Ok(mechanics)
            });
            match result {
                Ok(mut mechanics) => {
                    let next_walls = [
                        Piston::new(original[0].mass_kg(), mechanics.velocity_after_m_s[0])?,
                        Piston::new(original[1].mass_kg(), mechanics.velocity_after_m_s[1])?,
                    ];
                    let mut torque = [[0.0; 3]; 2];
                    for end in 0..2 {
                        let faces = grid.boundary_range(sides[end]);
                        torque[end] = wall_torque(
                            &grid,
                            sides[end],
                            self.boundary_trial[faces.start * width..faces.end * width]
                                .chunks_exact(width)
                                .map(|q| [q[n], q[n + 1], q[n + 2]]),
                            &mut cancelled,
                        )
                        .map_err(|error| match error {
                            fs_gas::dynamics::Error::Cancelled => ReactiveError::Cancelled,
                            fs_gas::dynamics::Error::Invalid(reason) => {
                                ReactiveError::Invalid(reason)
                            }
                            _ => ReactiveError::Invalid("reactive piston guide torque overflow"),
                        })?;
                    }
                    mechanics.gas_trials = gas_trials;
                    let report = ReactiveInertialWallsStep {
                        axis,
                        mechanics,
                        side_boundary_impulse_n_s: side_impulse,
                        guide_impulse_n_s: inward.map(|q| {
                            std::array::from_fn(|a| if a == axis { 0.0 } else { q[n + a] })
                        }),
                        gas_angular_impulse_n_m_s: torque,
                        guide_angular_impulse_n_m_s: torque.map(|q| q.map(|v| -v)),
                    };
                    poll(&mut cancelled)?;
                    self.publish_moving(mechanics.gas, operator.as_deref_mut());
                    *walls = next_walls;
                    return Ok(report);
                }
                Err(error) if terminal(&error) => return Err(error),
                Err(error) => {
                    if retry == transport_controls.maximum_coupling_retries {
                        return Err(ReactiveInertialWallsError::StepRetries {
                            attempts: retry + 1,
                            last: Box::new(error),
                        });
                    }
                    dt = shortened.unwrap_or(dt * 0.5);
                }
            }
        }
        unreachable!("bounded reacting piston attempts return on the final interval")
    }
}

fn terminal(error: &ReactiveInertialWallsError) -> bool {
    matches!(error, ReactiveInertialWallsError::Cancelled)
}
