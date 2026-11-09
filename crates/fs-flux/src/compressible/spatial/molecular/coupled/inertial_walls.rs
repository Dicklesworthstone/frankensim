//! Finite guided planes driven by native pressure, molecular stress and heat.

use super::*;
use crate::compressible::spatial::inertial_walls::wall_torque;
pub use fs_gas::dynamics::inertial_walls::{InertialWallLoad, InertialWallsControl};
use fs_gas::dynamics::inertial_walls::{InertialWallsStep, WallGasExchange, solve_coupled};
pub use fs_gas::dynamics::piston::{Piston, PistonError};

/// A complete transport or mechanical refusal preserves all accepted owners.
#[derive(Clone, Debug, PartialEq)]
pub enum MolecularInertialWallsError<E = MolecularError> {
    /// Original source, caloric, molecular, ALE or geometry refusal.
    Transport(E),
    /// Original finite inertia, nonlinear solve or physical-balance refusal.
    Mechanics(PistonError),
    /// Cancellation takes precedence over retries and publication.
    Cancelled,
    /// No complete coupled interval passed its original gates.
    StepRetries {
        /// Number of complete attempted intervals.
        attempts: usize,
        /// Last concrete refusal.
        last: Box<MolecularInertialWallsError<E>>,
    },
}

impl<E: core::fmt::Debug> core::fmt::Display for MolecularInertialWallsError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "molecular gas/inertial-wall coupling: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for MolecularInertialWallsError<E> {}
impl<E> From<PistonError> for MolecularInertialWallsError<E> {
    fn from(error: PistonError) -> Self {
        if error == PistonError::Cancelled {
            Self::Cancelled
        } else {
            Self::Mechanics(error)
        }
    }
}
impl From<MolecularError> for MolecularInertialWallsError {
    fn from(error: MolecularError) -> Self {
        if error == MolecularError::Cancelled {
            Self::Cancelled
        } else {
            Self::Transport(error)
        }
    }
}
impl From<fs_gas::multicomponent::dissipation::DissipationError>
    for MolecularInertialWallsError<fs_gas::multicomponent::dissipation::DissipationError>
{
    fn from(error: fs_gas::multicomponent::dissipation::DissipationError) -> Self {
        if error == fs_gas::multicomponent::dissipation::DissipationError::Cancelled {
            Self::Cancelled
        } else {
            Self::Transport(error)
        }
    }
}

/// Accepted molecular flow, finite wall response and explicit guide reactions.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolecularInertialWallsStep<R = MovingMolecularStepReport3d> {
    /// Translation axis: x=0, y=1, z=2. End order is [minimum, maximum].
    pub axis: usize,
    /// Actual pressure-plus-viscous wall impulse/work and finite-mass response.
    /// Combined axial momentum subtracts the separately reported side-boundary
    /// support impulse. Total energy subtracts directly integrated wall heat;
    /// its residual is in the original source's physical energy chart.
    pub mechanics: InertialWallsStep<R>,
    /// Inward gas impulse on all sides transverse to `axis` [N s]. These fixed
    /// supports can absorb axial viscous momentum, despite doing no wall work.
    pub side_boundary_impulse_n_s: [f64; 3],
    /// Impulse applied by each ideal translation guide to its piston [N s].
    /// Its axial component is zero. The other components cancel the actual
    /// tangential gas force so that each plane translates only along `axis`.
    pub guide_impulse_n_s: [[f64; 3]; 2],
    /// Actual gas angular impulse about each piston centre [N m s], including
    /// every pressure and viscous face impulse.
    pub gas_angular_impulse_n_m_s: [[f64; 3]; 2],
    /// Opposite orientation-guide impulse [N m s]. These guides do no work.
    pub guide_angular_impulse_n_m_s: [[f64; 3]; 2],
}

impl<L: TransportLaw> MolecularEuler3d<L> {
    /// Advance the complete molecular field and two finite guided pistons.
    ///
    /// Every nonlinear midpoint-speed proposal restarts the same accepted
    /// field and executes the original molecular/ALE/molecular split at one
    /// fixed duration. Native pressure and viscous impulses drive both masses;
    /// their conjugate mechanical work uses those very same moving walls.
    /// Fourier and Navier-partition heat are integrated separately and remain
    /// explicit thermal exchange, never a force or a kinetic-energy correction.
    ///
    /// The selected axis has two reflecting ends. All other directions are
    /// periodic or stationary reflecting walls; their actual viscous support
    /// impulse remains in the total momentum balance. Authored tangential wall
    /// velocities must be zero. No-slip, finite Navier slip, free slip and
    /// prescribed-temperature walls retain their original constitutive laws.
    /// Ideal translation/orientation guides carry the reported transverse and
    /// angular impulses. The thermal reservoir has no evolved finite capacity.
    ///
    /// Original source and parabolic/acoustic bounds apply to every proposal.
    /// Complete failures shorten the physical interval, never an individual
    /// molecular half inside the mechanical root. Cells, geometry, clock,
    /// loads, reports and both pistons publish after one final checkpoint.
    /// Prepared cell/face buffers are reused; no field clone is made per solve.
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
    ) -> Result<MolecularInertialWallsStep, MolecularInertialWallsError> {
        poll(&mut cancelled).map_err(MolecularError::from)?;
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
        let mut rate = 0.0_f64;
        for state in &self.accepted.states {
            let p = state
                .primitive(&self.accepted.gas, &mut cancelled)
                .map_err(MolecularError::from)?;
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
                return Err(PistonError::Invalid("molecular inertial-plane CFL overflow").into());
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
        let before = checked_totals(&self.accepted.states, grid.cell_volume(), &mut cancelled)
            .map_err(MolecularError::from)?;
        let original = *walls;
        let mut gas_trials = 0;
        let mut side_impulse = [0.0; 3];
        let mut inward = [Totals::default(); 2];
        for retry in 0..=transport_controls.maximum_coupling_retries {
            poll(&mut cancelled).map_err(MolecularError::from)?;
            let mut shortened = None;
            let result = solve_coupled(
                original,
                Time::new(dt),
                transport_controls.advection,
                loads,
                controls,
                area,
                |speed| -> Result<
                    WallGasExchange<MovingMolecularStepReport3d>,
                    MolecularInertialWallsError,
                > {
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
                            MolecularInertialWallsError::from(failure.error)
                        })?;
                    let mut work = [0.0; 2];
                    for end in 0..2 {
                        let range = grid.boundary_range(sides[end]);
                        let impulse = sum_arrays(
                            self.boundary_trial[range.clone()].iter().map(|q| q.array()),
                            &mut cancelled,
                        )
                        .map_err(MolecularError::from)?;
                        inward[end] = Totals::from_array(impulse.map(|q| -q));
                        // Molecular diffusion does not alter the central
                        // ALE owner's published pressure-only face account.
                        let pressure = sum_arrays(
                            self.trial.boundary_impulses[range]
                                .iter()
                                .map(|q| [-q.energy]),
                            &mut cancelled,
                        )
                        .map_err(MolecularError::from)?[0];
                        work[end] = pressure + gas.molecular_wall_work_j[2 * axis + end];
                    }
                    side_impulse = sum_arrays(
                        BoundarySide::ALL
                            .into_iter()
                            .filter(|side| side.axis() != axis)
                            .flat_map(|side| self.boundary_trial[grid.boundary_range(side)].iter())
                            .map(|q| q.momentum.map(|v| -v)),
                        &mut cancelled,
                    )
                    .map_err(MolecularError::from)?;
                    let after = gas.transport.after;
                    let energy = sum_arrays(
                        [
                            after.energy - before.energy,
                            (after.mass - before.mass) * self.accepted.gas.energy_reference(),
                        ]
                        .into_iter()
                        .chain(gas.wall_heat_j.map(|q| -q))
                        .map(|q| [q]),
                        &mut cancelled,
                    )
                    .map_err(MolecularError::from)?[0];
                    let displacement = displacement(grid, self.trial.grid, axis, speed)?;
                    Ok(WallGasExchange {
                        gas,
                        momentum_change_n_s: after.momentum[axis]
                            - before.momentum[axis]
                            - side_impulse[axis],
                        energy_change_j: energy,
                        boundary_impulse_n_s: inward.map(|q| q.momentum[axis]),
                        boundary_work_j: work,
                        displacement_m: displacement,
                    })
                },
                terminal,
            )
            .and_then(|mechanics| {
                // The physical source chart and the independently evolved
                // referenced-energy equation must both close with wall heat.
                let referenced = sum_arrays(
                    std::iter::once(mechanics.gas.transport.after.energy - before.energy)
                        .chain(mechanics.gas.wall_heat_j.map(|q| -q))
                        .chain((0..2).map(|i| {
                            mechanics.kinetic_energy_change_j[i]
                                + mechanics.exterior_pressure_work_j[i]
                                - mechanics.external_force_work_j[i]
                        }))
                        .map(|q| [q]),
                    &mut cancelled,
                )
                .map_err(MolecularError::from)?[0];
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
                        torque[end] = wall_torque(
                            &grid,
                            sides[end],
                            self.boundary_trial[grid.boundary_range(sides[end])]
                                .iter()
                                .map(|q| q.momentum),
                            &mut cancelled,
                        )
                        .map_err(MolecularError::from)?;
                    }
                    mechanics.gas_trials = gas_trials;
                    let report = MolecularInertialWallsStep {
                        axis,
                        mechanics,
                        side_boundary_impulse_n_s: side_impulse,
                        guide_impulse_n_s: inward.map(|q| {
                            std::array::from_fn(|a| if a == axis { 0.0 } else { q.momentum[a] })
                        }),
                        gas_angular_impulse_n_m_s: torque,
                        guide_angular_impulse_n_m_s: torque.map(|q| q.map(|v| -v)),
                    };
                    poll(&mut cancelled).map_err(MolecularError::from)?;
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

fn terminal(error: &MolecularInertialWallsError) -> bool {
    match error {
        MolecularInertialWallsError::Cancelled => true,
        MolecularInertialWallsError::Transport(error) => error.is_terminal(),
        _ => false,
    }
}

pub(in crate::compressible::spatial) fn admit_planes(
    grid: CartesianGrid,
    boundaries: Boundaries3d,
    walls: [MolecularWall3d; 6],
    axis: usize,
) -> Result<[BoundarySide; 2], PistonError> {
    if axis >= 3 {
        return Err(PistonError::Invalid("inertial wall axis must be x, y or z"));
    }
    let reflecting = AxisBoundaries::Ends {
        lower: BoundaryCondition::Reflecting,
        upper: BoundaryCondition::Reflecting,
    };
    if boundaries.axis(axis) != reflecting
        || (0..3).any(|a| {
            boundaries.axis(a) != reflecting && boundaries.axis(a) != AxisBoundaries::Periodic
        })
        || walls
            .iter()
            .any(|wall| wall.stress_boundary().0 != [0.0; 3])
    {
        return Err(PistonError::Invalid(
            "molecular inertial planes require closed sides and no prescribed tangential wall motion",
        ));
    }
    let sides = [BoundarySide::ALL[2 * axis], BoundarySide::ALL[2 * axis + 1]];
    if !positive(grid.face_area(sides[0]) * grid.face_count(sides[0]) as f64) {
        return Err(PistonError::Invalid("unrepresentable inertial plane area"));
    }
    Ok(sides)
}

pub(in crate::compressible::spatial) fn bounded_duration(
    maximum: f64,
    cfl: f64,
    rate: f64,
    grid: CartesianGrid,
    axis: usize,
    bounds: [[f64; 2]; 2],
) -> Result<f64, PistonError> {
    let mut dt = maximum.min(0.9 * cfl / rate);
    let contraction = bounds[0][1] - bounds[1][0];
    if !contraction.is_finite() {
        return Err(PistonError::Invalid(
            "inertial plane contraction bound overflow",
        ));
    }
    if contraction > 0.0 {
        let length = grid.shape[axis] as f64 * grid.spacing[axis];
        dt = dt.min(0.9 * cfl * (length / contraction));
    }
    Ok(dt)
}

pub(in crate::compressible::spatial) fn displacement(
    before: CartesianGrid,
    after: CartesianGrid,
    axis: usize,
    speed: [f64; 2],
) -> Result<[f64; 2], PistonError> {
    let positions = |g: CartesianGrid| {
        [
            g.origin[axis],
            g.origin[axis] + g.shape[axis] as f64 * g.spacing[axis],
        ]
    };
    let start = positions(before);
    let end = positions(after);
    let displacement = std::array::from_fn(|i| end[i] - start[i]);
    if (0..2).any(|i| !displacement[i].is_finite() || (speed[i] != 0.0 && end[i] == start[i])) {
        return Err(PistonError::Invalid(
            "inertial plane travel is not representable in the accepted geometry",
        ));
    }
    Ok(displacement)
}
