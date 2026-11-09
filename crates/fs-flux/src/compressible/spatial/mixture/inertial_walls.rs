//! Original evolving species driving two finite guided piston planes.
//!
//! The native two-wall mechanical solve consumes this field's actual ALE
//! pressure impulses and work. Every trial retains original local composition,
//! all momenta and the species energy chart, with no uniform-cavity substitute.

use super::super::inertial_walls::wall_torque;
pub use super::super::inertial_walls::{
    InertialWallLoad, InertialWallsControl, Piston, SpatialInertialWallsStep,
};
use super::*;
use crate::compressible::{AxisBoundaries, positive};
use fs_gas::dynamics::inertial_walls::{WallGasExchange, solve_coupled};
pub use fs_gas::dynamics::piston::PistonError;

impl MixtureEuler3d {
    /// Advance every original species and both finite end planes together.
    ///
    /// `axis` selects x=0, y=1 or z=2; wall order is [minimum, maximum]. Both
    /// wall speeds and extra forces use the positive coordinate direction.
    /// Positions and face areas belong to the accepted field grid. Both ends
    /// are reflecting; transverse sides are periodic or stationary slip walls.
    /// Finite masses, exterior pressures and additional axial forces are
    /// explicit. Ideal orientation guides carry the reported pressure torque.
    ///
    /// Every midpoint-speed proposal restarts the original source-owned ALE
    /// operator from the same accepted species field at one fixed duration.
    /// Shared mechanical controls bound Newton, line-search and complete-step
    /// retry work. Pressure impulse and its conjugate wall work update both
    /// finite inertias. The native field, boundary and root scratch are reused;
    /// no mixture source or full field is cloned per trial.
    ///
    /// Species cannot cross any piston. Referenced and physical gas-plus-wall
    /// energies must BOTH satisfy `energy_tolerance_j`; the latter includes
    /// each actual species-mass change paired with its original energy
    /// reference. Those products are summed from differences, without
    /// subtracting large absolute formation-energy totals. The returned
    /// mechanical total-energy residual is in that physical chart; the gas
    /// report retains all species and referenced/physical transport residuals.
    /// No caloric correction or species renormalization is applied.
    ///
    /// Gas cells, geometry, time, wall states, loads and report publish together
    /// only after the final gate. Refusal or cancellation preserves them all.
    /// This sealed, adiabatic inviscid lane excludes reaction, molecular
    /// transport, leakage, rotation, contact stops, prescribed exteriors and
    /// deformable piston surfaces. It inherits the original unequal-caloric
    /// contact limitation and makes no general FSI-stability claim.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_inertial_walls(
        &mut self,
        axis: usize,
        walls: &mut [Piston; 2],
        maximum_duration: Time,
        gas_controls: StepControl,
        loads: [InertialWallLoad; 2],
        controls: InertialWallsControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SpatialInertialWallsStep<StepReport>, PistonError> {
        stepping::poll::<Error>(&mut cancelled)?;
        controls.validate(*walls, maximum_duration, gas_controls, loads)?;
        if axis >= 3 {
            return Err(PistonError::Invalid("inertial wall axis must be x, y or z"));
        }
        let reflecting = AxisBoundaries::Ends {
            lower: BoundaryCondition::Reflecting,
            upper: BoundaryCondition::Reflecting,
        };
        if self.boundaries.axis(axis) != reflecting
            || (0..3).any(|a| {
                self.boundaries.axis(a) != reflecting
                    && self.boundaries.axis(a) != AxisBoundaries::Periodic
            })
        {
            return Err(PistonError::Invalid(
                "inertial mixture planes need closed ends and closed or periodic transverse sides",
            ));
        }
        self.admit_motion(MeshMotion::default())?;
        let grid = self.grid;
        let sides = [BoundarySide::ALL[2 * axis], BoundarySide::ALL[2 * axis + 1]];
        let area = grid.face_area(sides[0]) * grid.face_count(sides[0]) as f64;
        if !positive(area) {
            return Err(PistonError::Invalid("unrepresentable inertial plane area"));
        }
        let bounds = controls.midpoint_velocity_bracket_m_s;
        let low = bounds[0][0].min(bounds[1][0]);
        let high = bounds[0][1].max(bounds[1][1]);
        let n = self.gas.species_count();
        let width = n + 4;
        let mut rate = 0.0_f64;
        for (q, p) in self
            .states
            .chunks_exact(width)
            .zip(&mut self.work.primitive)
        {
            *p = cell(q, n).primitive(&self.gas, &mut cancelled)?;
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
                return Err(PistonError::Invalid(
                    "inertial mixture relative-wave CFL overflow",
                ));
            }
            rate = rate.max(local);
        }
        let mut dt = maximum_duration.value().min(0.9 * gas_controls.cfl / rate);
        let contraction = bounds[0][1] - bounds[1][0];
        if !contraction.is_finite() {
            return Err(PistonError::Invalid(
                "inertial plane contraction bound overflow",
            ));
        }
        if contraction > 0.0 {
            let length = grid.shape[axis] as f64 * grid.spacing[axis];
            dt = dt.min(0.9 * gas_controls.cfl * (length / contraction));
        }
        let original = *walls;
        let before = checked_totals(&self.states, &grid, n, &mut cancelled)?;
        let mut gas_trials = 0;
        let mut first_order_retry = false;
        let mut next_grid = grid;
        for halvings in 0..=gas_controls.maximum_halvings {
            stepping::poll::<Error>(&mut cancelled)?;
            if !positive(dt) || self.time + dt <= self.time || !(self.time + dt).is_finite() {
                return Err(PistonError::Invalid(
                    "inertial mixture planes cannot advance the clock",
                ));
            }
            let reconstruction = if first_order_retry {
                Reconstruction::FirstOrder
            } else {
                gas_controls.reconstruction
            };
            let result = solve_coupled(
                original,
                Time::new(dt),
                gas_controls,
                loads,
                controls,
                area,
                |speed| -> Result<WallGasExchange<StepReport>, PistonError> {
                    gas_trials += 1;
                    let mut motion = MeshMotion::default();
                    motion.lower_velocity[axis] = speed[0];
                    motion.upper_velocity[axis] = speed[1];
                    let candidate_grid = motion.grid_after(&grid, dt)?;
                    let stats = moving::attempt(
                        &self.gas,
                        &grid,
                        &candidate_grid,
                        self.boundaries,
                        &self.reservoirs,
                        motion,
                        &self.states,
                        &mut self.work,
                        dt,
                        gas_controls.cfl,
                        reconstruction,
                        &mut cancelled,
                    )?;
                    let after =
                        checked_totals(&self.work.candidate, &candidate_grid, n, &mut cancelled)?;
                    let gas = moving::reported_step(
                        &self.gas,
                        before,
                        after,
                        &self.work.boundary_trial,
                        dt,
                        stats,
                        halvings,
                        reconstruction != gas_controls.reconstruction,
                        &mut cancelled,
                    )?;
                    let mut inward = [[0.0; WIDTH_LIMIT]; 2];
                    for end in 0..2 {
                        let faces = grid.boundary_range(sides[end]);
                        inward[end] = sum_rows(
                            &self.work.boundary_trial[faces.start * width..faces.end * width],
                            width,
                            -1.0,
                            &mut cancelled,
                        )?;
                    }
                    let positions = |g: &CartesianGrid| {
                        [
                            g.origin[axis],
                            g.origin[axis] + g.shape[axis] as f64 * g.spacing[axis],
                        ]
                    };
                    let start = positions(&grid);
                    let end = positions(&candidate_grid);
                    let displacement = std::array::from_fn::<_, 2, _>(|i| end[i] - start[i]);
                    if (0..2).any(|i| {
                        !displacement[i].is_finite() || (speed[i] != 0.0 && end[i] == start[i])
                    }) {
                        return Err(PistonError::Invalid(
                            "inertial plane travel is not representable in the accepted geometry",
                        ));
                    }
                    // Closed species equations retain formation references as
                    // conserved storage, never a heat-release source.
                    let physical_change = sum(std::iter::once(after.energy - before.energy).chain(
                        (0..n).map(|d| {
                            (after.species_masses[d] - before.species_masses[d])
                                * self.gas.energy_references()[d]
                        }),
                    ));
                    next_grid = candidate_grid;
                    Ok(WallGasExchange {
                        gas,
                        momentum_change_n_s: after.momentum[axis] - before.momentum[axis],
                        energy_change_j: physical_change,
                        boundary_impulse_n_s: inward.map(|q| q[n + axis]),
                        boundary_work_j: inward.map(|q| q[n + 3]),
                        displacement_m: displacement,
                    })
                },
                |error| {
                    matches!(
                        error,
                        PistonError::Cancelled | PistonError::MixtureGas(Error::Allocation)
                    )
                },
            )
            .and_then(|mechanics| {
                let referenced = sum(std::iter::once(
                    mechanics.gas.after.energy - mechanics.gas.before.energy,
                )
                .chain((0..2).map(|i| {
                    mechanics.kinetic_energy_change_j[i] + mechanics.exterior_pressure_work_j[i]
                        - mechanics.external_force_work_j[i]
                })));
                if !referenced.is_finite() || referenced.abs() > controls.energy_tolerance_j {
                    return Err(PistonError::Balance {
                        momentum_residual_n_s: mechanics.combined_momentum_residual_n_s.abs(),
                        energy_residual_j: referenced.abs(),
                    });
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
                            self.work.boundary_trial[faces.start * width..faces.end * width]
                                .chunks_exact(width)
                                .map(|q| [q[n], q[n + 1], q[n + 2]]),
                            &mut cancelled,
                        )?;
                    }
                    mechanics.gas_trials = gas_trials;
                    let report = SpatialInertialWallsStep {
                        axis,
                        mechanics,
                        gas_angular_impulse_n_m_s: torque,
                        guide_angular_impulse_n_m_s: torque.map(|v| v.map(|x| -x)),
                    };
                    stepping::poll::<Error>(&mut cancelled)?;
                    std::mem::swap(&mut self.states, &mut self.work.candidate);
                    self.boundary_impulses
                        .copy_from_slice(&self.work.boundary_trial);
                    self.grid = next_grid;
                    self.time += dt;
                    self.last_step = Some(mechanics.gas);
                    *walls = next_walls;
                    return Ok(report);
                }
                Err(error @ PistonError::Cancelled)
                | Err(error @ PistonError::MixtureGas(Error::Allocation)) => return Err(error),
                Err(error) => {
                    if halvings == gas_controls.maximum_halvings {
                        return Err(PistonError::StepRetries {
                            attempts: halvings + 1,
                            last: Box::new(error),
                        });
                    }
                    if matches!(&error, PistonError::MixtureGas(cause) if !matches!(cause, Error::CflExceeded { .. }))
                    {
                        first_order_retry = true;
                    }
                    dt *= 0.5;
                }
            }
        }
        unreachable!("bounded coupled attempts return on the final interval")
    }
}

#[cfg(test)]
mod tests;
