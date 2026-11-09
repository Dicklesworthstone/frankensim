//! Two finite guided pistons driven by the original volumetric Euler field.
//!
//! The shared gas mechanics solves both midpoint speeds together. Every
//! proposal reuses the native ALE residual and extensive SSP RK2 update at one
//! fixed duration. Its actual pressure impulse and work update the same finite
//! masses; no endpoint-pressure approximation or reduced cavity replaces it.

use super::*;
pub use fs_gas::dynamics::inertial_walls::{
    InertialWallLoad, InertialWallsControl, InertialWallsStep,
};
use fs_gas::dynamics::inertial_walls::{WallGasExchange, solve_coupled};
pub use fs_gas::dynamics::piston::Piston;
use fs_gas::dynamics::piston::PistonError;

/// Original flow/mechanics refusals; no accepted state changes on error.
#[derive(Clone, Debug, PartialEq)]
pub enum InertialWallsError {
    /// Original volumetric gas, source, geometry or CFL refusal.
    Flow(Error),
    /// Original finite-wall control, inertia, root or physical-balance refusal.
    Mechanics(PistonError),
    /// Cancellation takes precedence over retries or publication.
    Cancelled,
    /// Every complete coupled interval was rejected.
    StepRetries {
        /// Number of complete attempts, including the initial interval.
        attempts: usize,
        /// Last concrete flow or mechanics refusal.
        last: Box<InertialWallsError>,
    },
}

impl core::fmt::Display for InertialWallsError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "volumetric gas/inertial-wall coupling: {self:?}")
    }
}
impl std::error::Error for InertialWallsError {}
impl From<Error> for InertialWallsError {
    fn from(error: Error) -> Self {
        if error == Error::Cancelled {
            Self::Cancelled
        } else {
            Self::Flow(error)
        }
    }
}
impl From<PistonError> for InertialWallsError {
    fn from(error: PistonError) -> Self {
        if error == PistonError::Cancelled {
            Self::Cancelled
        } else {
            Self::Mechanics(error)
        }
    }
}

/// Accepted field and finite-piston response on one physical time interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpatialInertialWallsStep<R = StepReport> {
    /// Translation axis: x=0, y=1, z=2. Wall order is [minimum, maximum].
    pub axis: usize,
    /// Original three-dimensional gas report and shared two-wall balances.
    /// Combined momentum here is the component along `axis`; stationary
    /// transverse-wall reactions remain in the original gas boundary account.
    pub mechanics: InertialWallsStep<R>,
    /// Actual gas angular impulse on each piston about its geometric centre
    /// [N m s], from the original accepted per-face pressure impulses.
    pub gas_angular_impulse_n_m_s: [[f64; 3]; 2],
    /// Equal and opposite angular impulse supplied by ideal orientation guides
    /// [N m s]. The guides hold both planes parallel and do no rotational work.
    pub guide_angular_impulse_n_m_s: [[f64; 3]; 2],
}

impl Euler3d {
    /// Advance the original gas and both finite, guided end planes together.
    ///
    /// The selected `axis` has two reflecting ends. Both pistons translate
    /// along its positive coordinate; their positions and common area belong
    /// to this field's accepted grid. The other directions are periodic or
    /// have stationary slip walls. Both finite masses, exterior pressures and
    /// additional centre-applied axial forces are explicit. Orientation guides
    /// carry the reported pressure torque instead of silently discarding it.
    ///
    /// The original two-wall solver selects midpoint velocities using native
    /// integrated pressure impulses. All Newton/Jacobian/backtracking proposals
    /// restart from the same cells at one duration and use the existing field
    /// scratch. A complete refusal halves that duration; no proposal advances
    /// physical time. Each wall's impulse/work and the combined axial momentum
    /// and gas-plus-wall energy must pass the caller's independent tolerances.
    /// No energy correction is applied. On acceptance, gas, geometry, time,
    /// wall states, per-face loads and the gas report publish atomically.
    ///
    /// This opt-in engineering lane retains the original frozen-composition,
    /// inviscid caloric model and source domain. It does not admit gravity,
    /// open transverse sides, retained exterior profiles, leakage, rotation,
    /// contact stops, heat transfer, or deformation of the piston surfaces.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_inertial_walls(
        &mut self,
        axis: usize,
        walls: &mut [Piston; 2],
        maximum_duration: Time,
        gas_controls: StepControl,
        loads: [InertialWallLoad; 2],
        controls: InertialWallsControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<SpatialInertialWallsStep, InertialWallsError> {
        self.step_inertial_walls_impl(
            axis,
            walls,
            maximum_duration,
            gas_controls,
            loads,
            controls,
            None,
            cancelled,
        )
    }

    /// Couple the original moving gas and finite pistons under held body forcing.
    ///
    /// Native source half-kicks surround every nonlinear proposal's ALE solve.
    /// Their actual old/new volumes supply all three body impulses and kinetic
    /// work. Only the pressure-face exchange drives the piston masses; gas-body
    /// impulse and work remain explicit in the gas report and are subtracted
    /// from the combined mechanical acceptance residual. Every forced cell,
    /// reconstructed trace and endpoint obeys the authored physical Mach cap
    /// and mesh-relative wave bound. Zero force follows the original method.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_inertial_walls_and_forcing(
        &mut self,
        axis: usize,
        walls: &mut [Piston; 2],
        maximum_duration: Time,
        gas_controls: StepControl,
        loads: [InertialWallLoad; 2],
        controls: InertialWallsControl,
        force: BodyForce3d<'_>,
        force_controls: ForceControl3d,
        workspace: &mut ForcingWorkspace3d,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SpatialInertialWallsStep, InertialWallsError> {
        force_controls.admit::<Error>()?;
        let active = force.admit::<Error>(self.states.len(), &mut cancelled)?;
        workspace.admit(self.states.len())?;
        if !active {
            return self.step_with_inertial_walls(
                axis,
                walls,
                maximum_duration,
                gas_controls,
                loads,
                controls,
                cancelled,
            );
        }
        self.step_inertial_walls_impl(
            axis,
            walls,
            maximum_duration,
            gas_controls,
            loads,
            controls,
            Some((force, force_controls, workspace)),
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_inertial_walls_impl(
        &mut self,
        axis: usize,
        walls: &mut [Piston; 2],
        maximum_duration: Time,
        gas_controls: StepControl,
        loads: [InertialWallLoad; 2],
        controls: InertialWallsControl,
        mut forcing: Option<forcing::Context<'_, '_>>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SpatialInertialWallsStep, InertialWallsError> {
        poll(&mut cancelled)?;
        controls.validate(*walls, maximum_duration, gas_controls, loads)?;
        if axis >= 3 {
            return Err(PistonError::Invalid("inertial wall axis must be x, y or z").into());
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
                "inertial planes need two reflecting ends and closed or periodic transverse sides",
            )
            .into());
        }
        MeshMotion::default().admit(self)?;
        let grid = self.grid;
        let sides = [BoundarySide::ALL[2 * axis], BoundarySide::ALL[2 * axis + 1]];
        let area = grid.face_area(sides[0]) * grid.face_count(sides[0]) as f64;
        if !positive(area) {
            return Err(PistonError::Invalid("unrepresentable inertial plane area").into());
        }
        let bounds = controls.midpoint_velocity_bracket_m_s;
        let low = bounds[0][0].min(bounds[1][0]);
        let high = bounds[0][1].max(bounds[1][1]);
        let mut rate = 0.0_f64;
        for (p, q) in self.work.primitive.iter_mut().zip(&self.states) {
            *p = q.primitive(&self.gas, &mut cancelled)?;
            if let Some((_, control, _)) = &forcing {
                control.admit_velocity::<Error>(p.velocity, p.sound_speed)?;
            }
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
                    PistonError::Invalid("inertial planes relative-wave CFL overflow").into(),
                );
            }
            rate = rate.max(local);
        }
        let mut dt = maximum_duration.value().min(0.9 * gas_controls.cfl / rate);
        let contraction = bounds[0][1] - bounds[1][0];
        if !contraction.is_finite() {
            return Err(PistonError::Invalid("inertial plane contraction bound overflow").into());
        }
        if contraction > 0.0 {
            let length = grid.shape[axis] as f64 * grid.spacing[axis];
            dt = dt.min(0.9 * gas_controls.cfl * (length / contraction));
        }
        let original = *walls;
        let before = checked_totals(&self.states, grid.cell_volume(), &mut cancelled)?;
        let mut gas_trials = 0;
        let mut first_order_retry = false;
        let mut next_grid = grid;
        for halvings in 0..=gas_controls.maximum_halvings {
            poll(&mut cancelled)?;
            if !positive(dt) || self.time + dt <= self.time || !(self.time + dt).is_finite() {
                return Err(
                    PistonError::Invalid("inertial planes cannot advance the clock").into(),
                );
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
                |speed| -> Result<WallGasExchange<StepReport>, InertialWallsError> {
                    gas_trials += 1;
                    let mut motion = MeshMotion::default();
                    motion.lower_velocity[axis] = speed[0];
                    motion.upper_velocity[axis] = speed[1];
                    let candidate_grid = motion.grid_after(&grid, dt)?;
                    let (stats, source) = match forcing.as_mut() {
                        Some((force, control, workspace)) => forcing::attempt_moving(
                            &self.gas,
                            &grid,
                            &candidate_grid,
                            self.boundaries,
                            motion,
                            &self.states,
                            &mut self.work,
                            workspace,
                            *force,
                            *control,
                            dt,
                            gas_controls.cfl,
                            reconstruction,
                            &mut cancelled,
                        )?,
                        None => (
                            moving::attempt(
                                &self.gas,
                                &grid,
                                &candidate_grid,
                                self.boundaries,
                                motion,
                                &self.states,
                                &mut self.work,
                                dt,
                                gas_controls.cfl,
                                reconstruction,
                                &mut cancelled,
                            )?,
                            ForceChange3d::default(),
                        ),
                    };
                    let after = checked_totals(
                        &self.work.candidate,
                        candidate_grid.cell_volume(),
                        &mut cancelled,
                    )?;
                    let mut gas = moving::reported_step(
                        before,
                        after,
                        &self.work.boundary_trial,
                        dt,
                        stats,
                        halvings,
                        reconstruction != gas_controls.reconstruction,
                        &mut cancelled,
                    )?;
                    if forcing.is_some() {
                        gas = moving::with_source(gas, source)?;
                    }
                    let mut inward = [Totals::default(); 2];
                    for end in 0..2 {
                        let sum = sum_arrays(
                            self.work.boundary_trial[grid.boundary_range(sides[end])]
                                .iter()
                                .map(|q| q.array()),
                            &mut cancelled,
                        )?;
                        inward[end] = Totals::from_array(sum.map(|q| -q));
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
                        )
                        .into());
                    }
                    next_grid = candidate_grid;
                    Ok(WallGasExchange {
                        gas,
                        momentum_change_n_s: after.momentum[axis]
                            - before.momentum[axis]
                            - source.impulse[axis],
                        energy_change_j: if forcing.is_some() {
                            sum_arrays(
                                [
                                    after.energy - before.energy,
                                    (after.mass - before.mass) * self.gas.energy_reference(),
                                    -source.work,
                                ]
                                .into_iter()
                                .map(|x| [x]),
                                &mut cancelled,
                            )?[0]
                        } else {
                            after.energy - before.energy
                        },
                        boundary_impulse_n_s: inward.map(|q| q.momentum[axis]),
                        boundary_work_j: inward.map(|q| q.energy),
                        displacement_m: displacement,
                    })
                },
                |error| {
                    matches!(
                        error,
                        InertialWallsError::Cancelled | InertialWallsError::Flow(Error::Allocation)
                    )
                },
            )
            .and_then(|mechanics| {
                if forcing.is_some() {
                    // Check the independently evolved referenced equation as
                    // well as the physical-energy gate in the shared solve.
                    let residual = sum_arrays(
                        [
                            mechanics.gas.after.energy - before.energy,
                            -mechanics.gas.source_change.work,
                        ]
                        .into_iter()
                        .chain((0..2).map(|end| {
                            mechanics.kinetic_energy_change_j[end]
                                + mechanics.exterior_pressure_work_j[end]
                                - mechanics.external_force_work_j[end]
                        }))
                        .map(|x| [x]),
                        &mut cancelled,
                    )?[0];
                    if !residual.is_finite() || residual.abs() > controls.energy_tolerance_j {
                        return Err(PistonError::Balance {
                            momentum_residual_n_s: mechanics.combined_momentum_residual_n_s.abs(),
                            energy_residual_j: residual.abs(),
                        }
                        .into());
                    }
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
                            self.work.boundary_trial[grid.boundary_range(sides[end])]
                                .iter()
                                .map(|q| q.momentum),
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
                    poll(&mut cancelled)?;
                    std::mem::swap(&mut self.states, &mut self.work.candidate);
                    self.boundary_impulses
                        .copy_from_slice(&self.work.boundary_trial);
                    self.grid = next_grid;
                    self.time += dt;
                    self.last_step = Some(mechanics.gas);
                    *walls = next_walls;
                    return Ok(report);
                }
                Err(error @ InertialWallsError::Cancelled)
                | Err(error @ InertialWallsError::Flow(Error::Allocation)) => return Err(error),
                Err(error) => {
                    if halvings == gas_controls.maximum_halvings {
                        return Err(InertialWallsError::StepRetries {
                            attempts: halvings + 1,
                            last: Box::new(error),
                        });
                    }
                    if matches!(&error, InertialWallsError::Flow(cause) if !matches!(cause, Error::CflExceeded { .. }))
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

pub(super) fn wall_torque(
    grid: &CartesianGrid,
    side: BoundarySide,
    impulses: impl Iterator<Item = [f64; 3]>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[f64; 3], Error> {
    let [a, b] = tangents(side.axis());
    let torque = sum_arrays(
        impulses.enumerate().map(|(face, j)| {
            // Relative coordinates avoid subtracting a large translated origin.
            // Motion is purely normal, so both stage centres have these arms.
            let mut arm = [0.0; 3];
            arm[a] = (face % grid.shape[a]) as f64 * grid.spacing[a] + 0.5 * grid.spacing[a]
                - 0.5 * (grid.shape[a] as f64 * grid.spacing[a]);
            arm[b] = (face / grid.shape[a]) as f64 * grid.spacing[b] + 0.5 * grid.spacing[b]
                - 0.5 * (grid.shape[b] as f64 * grid.spacing[b]);
            [
                arm[1] * j[2] - arm[2] * j[1],
                arm[2] * j[0] - arm[0] * j[2],
                arm[0] * j[1] - arm[1] * j[0],
            ]
        }),
        cancelled,
    )?;
    if torque.into_iter().any(|x| !x.is_finite()) {
        return Err(Error::Invalid("inertial plane angular impulse overflow"));
    }
    Ok(torque)
}

#[cfg(test)]
mod tests;
