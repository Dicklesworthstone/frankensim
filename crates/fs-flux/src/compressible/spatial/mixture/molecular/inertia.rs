//! Reciprocal tangential wall inertia on the original molecular SSP trial.
//!
//! The existing nonlinear owner finds a held midpoint velocity for each native
//! substep. Native outward tangent impulse changes carrier momentum and its
//! midpoint mechanical work changes kinetic storage. Thermal boundary energy
//! remains a separate reservoir exchange; normal impulse belongs to the guide.

use super::*;
use fs_qty::Mass;
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{NewtonKrylovConfig, NewtonKrylovState, SolverRunProgress};
use std::cell::{Cell, RefCell};

/// Finite mass of one uniformly translating, tangential no-slip wall.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InertialWall3d {
    pub side: BoundarySide,
    /// Complete translating mass [kg], not an areal density.
    pub mass: Mass,
}

/// Work and accuracy bounds for the original Newton--Krylov interface solve.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WallInertiaControl3d {
    /// Cumulative native spatial-stage evaluations, including probes and retries.
    pub maximum_field_evaluations: usize,
    /// Cumulative attempted SSP substeps across both halves and outer retries.
    pub maximum_molecular_trials: usize,
    /// At most 32 outer Newton iterations per SSP trial.
    pub maximum_newton_iterations: usize,
    /// Maximum-component numerical Jacobian perturbation [m/s].
    pub difference_step_m_s: f64,
    /// Independent per-component carrier impulse gate [kg m/s].
    pub impulse_tolerance_kg_m_s: f64,
    /// Independent carrier kinetic/mechanical-work gate [J].
    pub work_tolerance_j: f64,
    pub newton: NewtonKrylovConfig,
}
impl Default for WallInertiaControl3d {
    fn default() -> Self {
        Self {
            maximum_field_evaluations: 4096,
            maximum_molecular_trials: 1024,
            maximum_newton_iterations: 16,
            difference_step_m_s: 1e-5,
            impulse_tolerance_kg_m_s: 1e-10,
            work_tolerance_j: 1e-10,
            newton: NewtonKrylovConfig {
                relative_tolerance: 0.0,
                linear_restart: 12,
                max_linear_cycles: 4,
                ..NewtonKrylovConfig::default()
            },
        }
    }
}
impl WallInertiaControl3d {
    /// Conservative numerical-capacity reservation for this bounded solve.
    ///
    /// At most 12 coordinates, restart 12, 8 linear cycles and 32 Newton
    /// iterations bound the native Krylov vectors, difference arrays and
    /// cloned iteration telemetry below this reservation. Allocator metadata
    /// and arbitrary constitutive-provider allocations retain their original
    /// exclusions. The reservation is included in operator workspace_bytes().
    pub const fn solver_workspace_bytes() -> usize {
        64 * 1024
    }

    fn validate(self) -> Result<(), DissipationError> {
        self.newton
            .validate()
            .map_err(|_| invalid("finite-wall Newton controls"))?;
        if self.maximum_field_evaluations == 0
            || self.maximum_molecular_trials == 0
            || !(1..=32).contains(&self.maximum_newton_iterations)
            || self.newton.linear_restart > 12
            || self.newton.max_linear_cycles > 8
            || !positive(self.difference_step_m_s)
            || !positive(self.impulse_tolerance_kg_m_s)
            || !positive(self.work_tolerance_j)
        {
            return Err(invalid("finite-wall work, accuracy or workspace controls"));
        }
        Ok(())
    }
}

/// Mechanical storage and actual native exchange for one finite wall.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WallInertiaExchange3d {
    pub mass: Mass,
    pub velocity_before: [f64; 3],
    pub velocity_after: [f64; 3],
    /// Outward fluid impulse received by the two translating components.
    pub tangential_impulse: [f64; 3],
    /// Molecular normal impulse received by the fixed guide; no normal work.
    pub guide_impulse: f64,
    pub momentum_residual: [f64; 3],
    pub kinetic_energy_change_j: f64,
    /// Sum of each accepted midpoint velocity dotted with its native impulse.
    pub mechanical_work_j: f64,
    /// Native outward energy minus mechanical work, for the thermal reservoir.
    pub thermal_exchange_j: f64,
    /// Kinetic storage change minus received mechanical work.
    pub work_residual_j: f64,
}

/// Complete accepted molecular interval; counters cover the enclosing call.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WallInertiaReport3d {
    entries: [Option<WallInertiaExchange3d>; 6],
    pub field_evaluations: usize,
    pub molecular_trials: usize,
}
impl WallInertiaReport3d {
    pub fn wall(&self, side: BoundarySide) -> Option<&WallInertiaExchange3d> {
        self.entries[side_index(side)].as_ref()
    }
    pub(in crate::compressible::spatial::mixture) fn apply(
        &self,
        walls: &mut [MolecularWall3d; 6],
    ) {
        for (wall, entry) in walls.iter_mut().zip(self.entries) {
            if let (MolecularWall3d::NoSlip { velocity, .. }, Some(entry)) = (wall, entry) {
                *velocity = entry.velocity_after;
            }
        }
    }
}

#[derive(Clone, Copy)]
pub(super) struct WallInertia {
    pub masses: [f64; 6],
    pub controls: WallInertiaControl3d,
}
#[derive(Default)]
pub(in crate::compressible::spatial::mixture) struct WallWork {
    fields: Cell<usize>,
    trials: Cell<usize>,
}
impl WallWork {
    pub(super) fn field(&self, model: Option<WallInertia>) -> Result<(), DissipationError> {
        if let Some(model) = model {
            if self.fields.get() >= model.controls.maximum_field_evaluations {
                return Err(invalid("finite-wall cumulative field-evaluation budget"));
            }
            self.fields.set(self.fields.get() + 1);
        }
        Ok(())
    }
    pub(super) fn trial(&self, model: Option<WallInertia>) -> Result<(), DissipationError> {
        if let Some(model) = model {
            if self.trials.get() >= model.controls.maximum_molecular_trials {
                return Err(invalid("finite-wall cumulative SSP-trial budget"));
            }
            self.trials.set(self.trials.get() + 1);
        }
        Ok(())
    }
}

impl<L: TransportLaw> MolecularTransport3d<L> {
    /// Bind free tangential inertia after configuring the no-slip wall models.
    ///
    /// Initial velocity is the existing wall velocity. Masses remain immutable;
    /// subsequent wall reconfiguration on these sides refuses. Include
    /// WallInertiaControl3d::solver_workspace_bytes() in the original complete
    /// operator budget. Body forcing may drive the fluid; a finite wall receives
    /// only its actual reciprocal fluid shear, with no normal geometry motion.
    pub fn with_wall_inertia(
        mut self,
        walls: &[InertialWall3d],
        controls: WallInertiaControl3d,
    ) -> Result<Self, DissipationError> {
        controls.validate()?;
        if self.inertia.is_some() || walls.is_empty() || walls.len() > 6 {
            return Err(invalid(
                "one nonempty finite-wall binding of at most six sides",
            ));
        }
        let mut masses = [0.0; 6];
        for wall in walls {
            let index = side_index(wall.side);
            if masses[index] != 0.0
                || !positive(wall.mass.value())
                || physical_condition(self.boundaries, wall.side)
                    != Some(BoundaryCondition::Reflecting)
                || !matches!(self.walls[index], MolecularWall3d::NoSlip { .. })
            {
                return Err(invalid(
                    "distinct reflecting no-slip walls with finite positive mass required",
                ));
            }
            let speed = velocity(self.walls[index])
                .into_iter()
                .fold(0.0_f64, f64::hypot);
            if !(wall.mass.value() * speed).is_finite()
                || !(0.5 * wall.mass.value() * speed * speed).is_finite()
            {
                return Err(invalid(
                    "finite-wall initial momentum or kinetic storage overflow",
                ));
            }
            masses[index] = wall.mass.value();
        }
        self.inertia = Some(WallInertia { masses, controls });
        let bytes = self
            .workspace_bytes()
            .checked_add(self.bound_flow_workspace_bytes)
            .ok_or(invalid("finite-wall complete workspace size overflow"))?;
        admit_workspace(bytes, self.maximum_workspace_bytes)?;
        Ok(self)
    }

    pub(super) fn solve_wall_trial(
        &mut self,
        dt: f64,
        controls: DissipationControl,
        drive: Option<(&WallDrive3d, f64)>,
        elapsed: f64,
        initial: &[MolecularWall3d; 6],
        work: &WallWork,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Stage, [MolecularWall3d; 6], [f64; 6]), DissipationError> {
        let model = self.inertia.expect("finite-wall trial");
        let mut coordinates = [(0usize, 0usize); 12];
        let mut count = 0;
        let mut guess = Vec::with_capacity(12);
        for side in BoundarySide::ALL {
            let index = side_index(side);
            if model.masses[index] != 0.0 {
                let velocity = velocity(initial[index]);
                for axis in tangents(side.axis()) {
                    coordinates[count] = (index, axis);
                    count += 1;
                    guess.push(velocity[axis]);
                }
            }
        }
        let owner = RefCell::new(self);
        let cancel = RefCell::new(cancelled);
        let evaluate = |x: &[f64]| {
            poll(&mut || (*cancel.borrow_mut())())?;
            let mut midpoint = *initial;
            for (value, &(index, axis)) in x.iter().zip(&coordinates[..count]) {
                if !value.is_finite() {
                    return Err(invalid("finite-wall midpoint velocity"));
                }
                if let MolecularWall3d::NoSlip { velocity, .. } = &mut midpoint[index] {
                    velocity[axis] = *value;
                }
            }
            let mut owner = owner.borrow_mut();
            let stage = owner.ssp_trial(
                dt,
                controls,
                drive,
                elapsed,
                Some(&midpoint),
                None,
                work,
                &mut || (*cancel.borrow_mut())(),
            )?;
            let impulse = owner.wall_impulses(&model);
            let residual = coordinates[..count]
                .iter()
                .map(|&(index, axis)| {
                    2.0 * model.masses[index]
                        * (velocity(midpoint[index])[axis] - velocity(initial[index])[axis])
                        - impulse[index][axis]
                })
                .collect::<Vec<_>>();
            if residual.iter().any(|v| !v.is_finite()) {
                return Err(invalid("finite-wall midpoint impulse overflow"));
            }
            Ok((stage, midpoint, impulse, residual))
        };
        let problem = FallibleFiniteDifference::new(
            count,
            model.controls.difference_step_m_s,
            |x: &[f64]| evaluate(x).map(|(_, _, _, residual)| residual),
        )
        .map_err(numerical)?;
        let construction = NewtonKrylovState::new(&problem, guess, model.controls.newton);
        if let Some(error) = problem.take_failure() {
            return Err(numerical(error));
        }
        let mut solver = construction.map_err(|_| invalid("finite-wall Newton initialization"))?;
        let result = solver.run_controlled(
            &problem,
            model.controls.maximum_newton_iterations,
            &mut || (*cancel.borrow_mut())() || problem.has_failed(),
        );
        if let Some(error) = problem.take_failure() {
            return Err(numerical(error));
        }
        match result.progress {
            SolverRunProgress::Complete => {}
            SolverRunProgress::Paused => return Err(DissipationError::Cancelled),
            _ => return Err(invalid("finite-wall nonlinear callback interrupted")),
        }
        if !result.report.converged {
            return Err(invalid("finite-wall nonlinear solve did not converge"));
        }
        // The final physical replay, not a differenced/line-search probe, owns
        // the candidate field and face impulses left in the original buffers.
        let (stage, midpoint, impulse, residual) = evaluate(&solver.x)?;
        if residual
            .iter()
            .any(|v| v.abs() > model.controls.impulse_tolerance_kg_m_s)
        {
            return Err(invalid("finite-wall independent impulse gate"));
        }
        let mut endpoint = *initial;
        let mut mechanical_work = [0.0; 6];
        for side in BoundarySide::ALL {
            let index = side_index(side);
            let mass = model.masses[index];
            if mass == 0.0 {
                continue;
            }
            let before = velocity(initial[index]);
            let mid = velocity(midpoint[index]);
            let after =
                std::array::from_fn(|axis| 2.0_f64.mul_add(mid[axis] - before[axis], before[axis]));
            if after.iter().any(|v| !v.is_finite()) {
                return Err(invalid("finite-wall endpoint velocity overflow"));
            }
            // Gate the actually represented endpoint too: rounding the
            // midpoint extrapolation must not hide impulse error in a later
            // substep with opposite sign.
            for axis in tangents(side.axis()) {
                let actual = mass * (after[axis] - before[axis]) - impulse[index][axis];
                if !actual.is_finite() || actual.abs() > model.controls.impulse_tolerance_kg_m_s {
                    return Err(invalid("finite-wall represented endpoint impulse gate"));
                }
            }
            let work = super::super::sum((0..3).map(|axis| mid[axis] * impulse[index][axis]));
            let kinetic = kinetic_change(mass, before, after);
            if !work.is_finite()
                || !kinetic.is_finite()
                || (kinetic - work).abs() > model.controls.work_tolerance_j
            {
                return Err(invalid("finite-wall kinetic and mechanical-work gate"));
            }
            if let MolecularWall3d::NoSlip { velocity, .. } = &mut endpoint[index] {
                *velocity = after;
            }
            mechanical_work[index] = work;
        }
        poll(&mut || (*cancel.borrow_mut())())?;
        Ok((stage, endpoint, mechanical_work))
    }

    fn wall_impulses(&self, model: &WallInertia) -> [[f64; 3]; 6] {
        std::array::from_fn(|index| {
            if model.masses[index] == 0.0 {
                return [0.0; 3];
            }
            let side = BoundarySide::ALL[index];
            std::array::from_fn(|axis| {
                if axis == side.axis() {
                    return 0.0;
                }
                super::super::sum(
                    boundary_view(
                        &self.work.boundary_step,
                        &self.grid,
                        self.gas().species_count(),
                        side,
                    )
                    .map(|q| q.momentum[axis]),
                )
            })
        })
    }

    pub(super) fn inertia_report(
        &self,
        initial: &[MolecularWall3d; 6],
        endpoint: &[MolecularWall3d; 6],
        mechanical_work: [f64; 6],
        work: &WallWork,
    ) -> Result<Option<WallInertiaReport3d>, DissipationError> {
        let Some(model) = self.inertia else {
            return Ok(None);
        };
        let mut entries = [None; 6];
        for side in BoundarySide::ALL {
            let index = side_index(side);
            let mass = model.masses[index];
            if mass == 0.0 {
                continue;
            }
            let mut impulse = std::array::from_fn(|axis| {
                super::super::sum(
                    boundary_view(
                        &self.work.boundary_trial,
                        &self.grid,
                        self.gas().species_count(),
                        side,
                    )
                    .map(|q| q.momentum[axis]),
                )
            });
            let guide_impulse = impulse[side.axis()];
            impulse[side.axis()] = 0.0;
            let outward_energy = super::super::sum(
                boundary_view(
                    &self.work.boundary_trial,
                    &self.grid,
                    self.gas().species_count(),
                    side,
                )
                .map(|q| q.energy),
            );
            let before = velocity(initial[index]);
            let after = velocity(endpoint[index]);
            let residual =
                std::array::from_fn(|axis| mass * (after[axis] - before[axis]) - impulse[axis]);
            let kinetic = kinetic_change(mass, before, after);
            let thermal = outward_energy - mechanical_work[index];
            let energy_residual = kinetic - mechanical_work[index];
            if residual
                .iter()
                .any(|v| !v.is_finite() || v.abs() > model.controls.impulse_tolerance_kg_m_s)
                || !thermal.is_finite()
                || !guide_impulse.is_finite()
                || !energy_residual.is_finite()
                || energy_residual.abs() > model.controls.work_tolerance_j
            {
                return Err(invalid(
                    "finite-wall complete-interval momentum or work gate",
                ));
            }
            entries[index] = Some(WallInertiaExchange3d {
                mass: Mass::new(mass),
                velocity_before: before,
                velocity_after: after,
                tangential_impulse: impulse,
                guide_impulse,
                momentum_residual: residual,
                kinetic_energy_change_j: kinetic,
                mechanical_work_j: mechanical_work[index],
                thermal_exchange_j: thermal,
                work_residual_j: energy_residual,
            });
        }
        Ok(Some(WallInertiaReport3d {
            entries,
            field_evaluations: work.fields.get(),
            molecular_trials: work.trials.get(),
        }))
    }
}

fn velocity(wall: MolecularWall3d) -> [f64; 3] {
    match wall {
        MolecularWall3d::NoSlip { velocity, .. } => velocity,
        _ => [0.0; 3],
    }
}
fn kinetic_change(mass: f64, before: [f64; 3], after: [f64; 3]) -> f64 {
    0.5 * mass
        * super::super::sum(
            (0..3).map(|axis| (after[axis] - before[axis]) * (after[axis] + before[axis])),
        )
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn invalid(reason: &'static str) -> DissipationError {
    DissipationError::Invalid(reason)
}
fn numerical(error: FiniteDifferenceError<DissipationError>) -> DissipationError {
    match error {
        FiniteDifferenceError::Source(error) => error,
        FiniteDifferenceError::Invalid(_) => invalid("finite-wall numerical Jacobian stencil"),
    }
}
