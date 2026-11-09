//! Reciprocal finite piston response to original species, stress and heat.

use super::*;
use crate::compressible::spatial::MeshMotion;
use crate::compressible::spatial::inertial_walls::wall_torque;
use crate::compressible::spatial::mixture::sum;
pub use crate::compressible::spatial::molecular::{
    InertialWallLoad, InertialWallsControl, Piston, PistonError,
};
use crate::compressible::spatial::molecular::{admit_planes, bounded_duration, displacement};
use fs_gas::dynamics::inertial_walls::{WallGasExchange, solve_coupled};

/// Shared atomic coupling error, defaulting to the original mixture source.
pub type MolecularInertialWallsError<E = DissipationError> =
    crate::compressible::spatial::molecular::MolecularInertialWallsError<E>;
/// Shared reciprocal wall receipt, defaulting to the moving mixture report.
pub type MolecularInertialWallsStep<R = MovingMolecularStepReport3d> =
    crate::compressible::spatial::molecular::MolecularInertialWallsStep<R>;

type CouplingError = MolecularInertialWallsError<DissipationError>;

impl<L: TransportLaw> MolecularEuler3d<L> {
    /// Advance the evolving molecular mixture and two finite guided pistons.
    ///
    /// Every midpoint-speed proposal starts from the same accepted mixture and
    /// performs its complete D/2--ALE--D/2 interval on actual old/new volumes.
    /// The original pressure and full molecular stresses transfer conjugate
    /// impulse/work to both finite masses. Source-owned Maxwell--Stefan,
    /// enthalpy transport and Fourier heat remain active throughout each trial.
    ///
    /// `axis` is x=0, y=1, z=2; end order is [minimum, maximum]. The selected
    /// ends reflect; other axes are periodic or fixed reflecting walls. Every
    /// authored tangential wall velocity must be zero. Finite Navier slip and
    /// prescribed wall temperatures retain their original constitutive laws.
    /// Transverse fixed walls can absorb axial viscous momentum; their actual
    /// support impulse is reported separately from the two pistons. Ideal
    /// translation/orientation guides carry transverse and angular impulse.
    ///
    /// Molecular wall heat is integrated directly, before mechanical work is
    /// added, and remains external thermal exchange. Both physical and species-
    /// referenced gas-plus-piston energies must pass the mechanical tolerance.
    /// Formation references multiply actual species-mass differences, avoiding
    /// subtraction of large absolute chemical energies. No caloric correction,
    /// species normalization or independent passive-scalar transport is used.
    ///
    /// Native source, acoustic and parabolic controls bound each proposal.
    /// A refused or shortened gas proposal restarts the complete interval.
    /// Publication of cells, geometry, time, reports, loads and both pistons is
    /// atomic after the final checkpoint. Prepared field buffers are reused.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_inertial_walls(
        &mut self,
        axis: usize,
        walls: &mut [Piston; 2],
        maximum_duration: Time,
        transport_controls: MolecularStepControl,
        loads: [InertialWallLoad; 2],
        controls: InertialWallsControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<
        MolecularInertialWallsStep<MovingMolecularStepReport3d>,
        MolecularInertialWallsError<DissipationError>,
    > {
        poll(&mut cancelled)?;
        controls.validate(
            *walls,
            maximum_duration,
            transport_controls.advection,
            loads,
        )?;
        self.admit_moving(
            maximum_duration.value(),
            transport_controls,
            MeshMotion::default(),
        )?;
        let grid = self.accepted.grid;
        let sides = admit_planes(grid, self.accepted.boundaries, self.molecular.walls, axis)?;
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
                .map_err(DissipationError::from)?;
            let local: f64 = (0..3)
                .map(|a| {
                    let speed = if a == axis {
                        (p.velocity[a] - low)
                            .abs()
                            .max((p.velocity[a] - high).abs())
                    } else {
                        p.velocity[a].abs()
                    };
                    (speed + p.sound_speed) / grid.spacing[a]
                })
                .sum();
            if !positive(local) {
                return Err(
                    PistonError::Invalid("molecular mixture inertial-plane CFL overflow").into(),
                );
            }
            rate = rate.max(local);
        }
        let mut dt = bounded_duration(
            maximum_duration.value(),
            transport_controls.advection.cfl,
            rate,
            grid,
            axis,
            bounds,
        )?;
        let before = checked_totals(&self.accepted.states, &grid, n, &mut cancelled)
            .map_err(DissipationError::from)?;
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
                transport_controls.advection,
                loads,
                controls,
                area,
                |speed| -> Result<WallGasExchange<MovingMolecularStepReport3d>, CouplingError> {
                    gas_trials += 1;
                    let mut motion = MeshMotion::default();
                    motion.lower_velocity[axis] = speed[0];
                    motion.upper_velocity[axis] = speed[1];
                    let gas = self
                        .trial_moving(
                            dt,
                            transport_controls,
                            motion,
                            before,
                            retry,
                            &mut cancelled,
                        )
                        .map_err(|failure| {
                            shortened = failure.shortened;
                            CouplingError::from(failure.error)
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
                        .map_err(DissipationError::from)?;
                        // The central native ALE owner retains only pressure
                        // work; molecular halves retain their own full work.
                        let pressure = sum_rows(
                            &self.trial.boundary_impulses[values],
                            width,
                            -1.0,
                            &mut cancelled,
                        )
                        .map_err(DissipationError::from)?[n + 3];
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
                        .map_err(DissipationError::from)?;
                        for a in 0..3 {
                            let adjusted = values[n + a] - correction[a];
                            let next = side_impulse[a] + adjusted;
                            correction[a] = (next - side_impulse[a]) - adjusted;
                            side_impulse[a] = next;
                        }
                    }
                    let after = gas.transport.after;
                    let physical_change = sum(std::iter::once(after.energy - before.energy)
                        .chain((0..n).map(|k| {
                            (after.species_masses[k] - before.species_masses[k])
                                * self.accepted.gas.energy_references()[k]
                        }))
                        .chain(gas.wall_heat_j.map(|q| -q)));
                    let displacement = displacement(grid, self.trial.grid, axis, speed)?;
                    Ok(WallGasExchange {
                        gas,
                        momentum_change_n_s: after.momentum[axis]
                            - before.momentum[axis]
                            - side_impulse[axis],
                        energy_change_j: physical_change,
                        boundary_impulse_n_s: inward.map(|q| q[n + axis]),
                        boundary_work_j: work,
                        displacement_m: displacement,
                    })
                },
                terminal,
            )
            .and_then(|mechanics| {
                // Neither energy chart may hide a defect in the other.
                let referenced = sum(std::iter::once(
                    mechanics.gas.transport.after.energy - before.energy,
                )
                .chain(mechanics.gas.wall_heat_j.map(|q| -q))
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
                            fs_gas::dynamics::Error::Cancelled => DissipationError::Cancelled,
                            fs_gas::dynamics::Error::Invalid(reason) => {
                                DissipationError::Invalid(reason)
                            }
                            _ => {
                                DissipationError::Invalid("molecular mixture guide torque overflow")
                            }
                        })?;
                    }
                    mechanics.gas_trials = gas_trials;
                    let report = MolecularInertialWallsStep {
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
                    self.publish_moving(mechanics.gas);
                    *walls = next_walls;
                    return Ok(report);
                }
                Err(error) if terminal(&error) => return Err(error),
                Err(error) => {
                    if retry == transport_controls.maximum_coupling_retries {
                        return Err(MolecularInertialWallsError::StepRetries {
                            attempts: retry + 1,
                            last: Box::new(error),
                        });
                    }
                    dt = shortened.unwrap_or(dt * 0.5);
                }
            }
        }
        unreachable!("bounded coupled attempts return on the final interval")
    }
}

fn terminal(error: &CouplingError) -> bool {
    match error {
        MolecularInertialWallsError::Cancelled => true,
        MolecularInertialWallsError::Transport(error) => error.is_terminal(),
        _ => false,
    }
}
