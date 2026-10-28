//! Native rigid shear/spin coupled to the existing affine contact free energy.
//!
//! Three impulse coordinates determine the original bodies' midpoint material
//! point motion. The original retained contact/enthalpy owner resolves physical
//! microslip heat, per-surface reversible entropy exchange, source friction and
//! conduction-derived partition at that same motion. No effective inertia,
//! replacement contact history, additional heat capacity or heat fraction is used.
//!
//! This is a frozen-configuration phase for one supplied compressive elliptical
//! patch. Normal dynamics and subsequent configuration evolution remain with
//! their owners. Endpoint constitutive integration remains first order; an
//! inner Newton residual is not a temporal-accuracy or uniqueness claim.
//! The original thermal cells are the designated physical contact volumes;
//! they need not represent the complete thermal inventory of either rigid body.

use super::resolved::{
    ResolvedCaloricContactCandidate, ResolvedCaloricContactReport, ResolvedCaloricContactSolve,
    ResolvedCaloricContactSource,
};
use super::*;
use fs_mbd::{DynamicsError, MassProperties, RigidBodyState, Vec3};

/// Simultaneous retained contacts sharing original rigid bodies and caloric cells.
pub mod network;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidCaloricContactSolve {
    pub contact: ResolvedCaloricContactSolve,
    pub maximum_motion_trials: usize,
    pub maximum_newton_iterations: usize,
    pub difference_step_n_s: f64,
    /// Third impulse coordinate is torsional impulse divided by this length.
    /// It conditions the solver; it does not replace a physical patch radius.
    pub torsional_coordinate_length_m: f64,
    pub geometry_tolerance_m: f64,
    pub linear_momentum_tolerance_n_s: f64,
    pub angular_momentum_tolerance_nm_s: f64,
    pub energy_tolerance_j: f64,
    pub newton: NewtonKrylovConfig,
}

#[derive(Debug)]
pub enum RigidCaloricContactError<E = Infallible> {
    Invalid(&'static str),
    Dynamics(DynamicsError),
    Contact(CaloricContactError<E>),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    Interrupted(SolverRunProgress),
    MotionBudget,
    Cancelled,
}
impl<E> From<DynamicsError> for RigidCaloricContactError<E> {
    fn from(error: DynamicsError) -> Self {
        Self::Dynamics(error)
    }
}
impl<E> From<CaloricContactError<E>> for RigidCaloricContactError<E> {
    fn from(error: CaloricContactError<E>) -> Self {
        Self::Contact(error)
    }
}
impl<E> From<NewtonError> for RigidCaloricContactError<E> {
    fn from(error: NewtonError) -> Self {
        Self::Newton(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for RigidCaloricContactError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native rigid caloric contact: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for RigidCaloricContactError<E> {}

pub struct RigidCaloricContactReport<S> {
    pub contact: ResolvedCaloricContactReport<S>,
    pub nonlinear: NewtonReport,
    pub motion_trials: usize,
    pub contact_trials: usize,
    pub applied_impulse_on_a_n_s: Vec3,
    pub applied_spin_impulse_on_a_nm_s: Vec3,
    pub midpoint_relative_velocity_mps: [f64; 2],
    pub midpoint_relative_spin_rad_s: f64,
    pub body_momentum_residual_n_s: [Vec3; 2],
    pub body_angular_momentum_residual_nm_s: [Vec3; 2],
    pub linear_momentum_change_n_s: Vec3,
    pub world_angular_momentum_change_nm_s: Vec3,
    pub kinetic_energy_change_j: f64,
    /// Actual native kinetic change plus work entering the original contact.
    pub mechanical_energy_residual_j: f64,
    /// Delta(K + U_contact + H_caloric) + D_numerical - Q_external.
    pub total_energy_residual_j: f64,
}

/// All original owner endpoints; publish the returned parts together.
pub struct RigidCaloricContactCandidate<'a, 't, S> {
    bodies: [RigidBodyState; 2],
    thermal: DepositionSimulation<'t>,
    contact: CaloricPartialSlip<'a>,
    report: RigidCaloricContactReport<S>,
}
impl<'a, 't, S> RigidCaloricContactCandidate<'a, 't, S> {
    pub const fn bodies(&self) -> &[RigidBodyState; 2] {
        &self.bodies
    }
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn contact(&self) -> &CaloricPartialSlip<'a> {
        &self.contact
    }
    pub const fn report(&self) -> &RigidCaloricContactReport<S> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        [RigidBodyState; 2],
        DepositionSimulation<'t>,
        CaloricPartialSlip<'a>,
        RigidCaloricContactReport<S>,
    ) {
        (self.bodies, self.thermal, self.contact, self.report)
    }
}

impl<'a> CaloricPartialSlip<'a> {
    /// Close native body momenta, retained shear/spin and finite thermal feedback.
    ///
    /// Bodies/arms follow the original ordered A/B caloric surfaces. Each actual
    /// body point must coincide with their shared thermal face centre; no point
    /// is moved to repair geometry. Only the contact impulse acts in this phase.
    /// Original poses and masses/inertias remain fixed. Other force/torque kicks
    /// and geometry motion belong to the surrounding integrator.
    ///
    /// The inner owner retains its original two-temperature affine free energy,
    /// source-bound friction, and caloric objects. Source speed axes use actual
    /// midpoint translational slip; a source without spin applicability continues
    /// to refuse torsional motion and admits only the two shear impulse unknowns;
    /// its third initial guess must be zero. The callback provides only OTHER heat and
    /// endpoint conductivity, never a prescribed contact heat split.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn step_rigid_resolved_candidate<'t, S, E>(
        &self,
        bodies: &[RigidBodyState; 2],
        properties: [MassProperties; 2],
        arms_body_m: [Vec3; 2],
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        frame: TangentFrame,
        creepage_reference_speed_mps: f64,
        ownership: &GeneralizedWorkOwnership,
        controls: RigidCaloricContactSolve,
        impulse_guess_n_s: [f64; 3],
        energy_guess_j: &[f64],
        interface_temperature_guess_k: f64,
        resolve: impl Fn(&DepositionSimulation<'t>, f64) -> Result<ResolvedCaloricContactSource<S>, E>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<RigidCaloricContactCandidate<'a, 't, S>, RigidCaloricContactError<E>> {
        let cancel = RefCell::new(cancelled);
        let poll = || {
            if (*cancel.borrow_mut())() {
                Err(RigidCaloricContactError::Cancelled)
            } else {
                Ok(())
            }
        };
        poll()?;
        controls.newton.validate()?;
        let dt = target_s - self.time_s();
        if !dt.is_finite()
            || dt <= 0.0
            || controls.maximum_motion_trials == 0
            || controls.maximum_newton_iterations == 0
            || impulse_guess_n_s.iter().any(|v| !v.is_finite())
            || [
                controls.difference_step_n_s,
                controls.torsional_coordinate_length_m,
                controls.geometry_tolerance_m,
                controls.linear_momentum_tolerance_n_s,
                controls.angular_momentum_tolerance_nm_s,
                controls.energy_tolerance_j,
                creepage_reference_speed_mps,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
        {
            return Err(RigidCaloricContactError::Invalid(
                "positive bounded rigid/contact solve required",
            ));
        }
        self.admit_rigid_geometry(
            bodies,
            properties,
            arms_body_m,
            thermal,
            frame,
            controls.geometry_tolerance_m,
        )?;
        let axes = [
            vector(frame.longitudinal()),
            vector(frame.lateral()),
            vector(frame.normal()),
        ];
        let before = [
            bodies[0].point_kinematics(properties[0], arms_body_m[0])?,
            bodies[1].point_kinematics(properties[1], arms_body_m[1])?,
        ];
        let length = controls.torsional_coordinate_length_m;
        // The original material binding has no torsional speed applicability.
        // Do not finite-difference an inadmissible spin coordinate when its
        // translational branch is used. Actual mechanically induced spin still
        // reaches the original owner and is refused there, never projected out.
        let dimension = if self.material_friction.is_some() {
            2
        } else {
            3
        };
        if dimension == 2 && impulse_guess_n_s[2] != 0.0 {
            return Err(RigidCaloricContactError::Invalid(
                "translational material source cannot own a spin impulse",
            ));
        }
        let trials = Cell::new(0usize);
        let contact_trials = Cell::new(0usize);
        let evaluate = |x: &[f64]| -> Result<
            (
                [RigidBodyState; 2],
                ResolvedCaloricContactCandidate<'a, 't, S>,
                [f64; 3],
                [f64; 2],
                f64,
            ),
            RigidCaloricContactError<E>,
        > {
            poll()?;
            if trials.get() >= controls.maximum_motion_trials {
                return Err(RigidCaloricContactError::MotionBudget);
            }
            trials.set(trials.get() + 1);
            if x.len() != dimension || x.iter().any(|v| !v.is_finite()) {
                return Err(RigidCaloricContactError::Invalid(
                    "finite impulse coordinates matching source applicability required",
                ));
            }
            let impulse = axes[0].scale(x[0]).add(axes[1].scale(x[1]));
            let spin_coordinate = x.get(2).copied().unwrap_or(0.0);
            let couple = axes[2].scale(length * spin_coordinate);
            let mut next = *bodies;
            for i in 0..2 {
                let sign = if i == 0 { 1.0 } else { -1.0 };
                let force = next[i]
                    .apply_impulse_at_body_point(
                        properties[i],
                        arms_body_m[i],
                        impulse.scale(sign),
                    )?
                    .state_after;
                next[i] = RigidBodyState::new(
                    force.pose(),
                    force.linear_momentum_world(),
                    force.angular_momentum_body().add(
                        force
                            .pose()
                            .orientation()
                            .rotate_world_to_body(couple.scale(sign)),
                    ),
                )?;
            }
            let after = [
                next[0].point_kinematics(properties[0], arms_body_m[0])?,
                next[1].point_kinematics(properties[1], arms_body_m[1])?,
            ];
            let relative = before[0]
                .point_velocity_world
                .add(after[0].point_velocity_world)
                .sub(before[1].point_velocity_world)
                .sub(after[1].point_velocity_world)
                .scale(0.5);
            let relative_spin = before[0]
                .angular_velocity_world
                .add(after[0].angular_velocity_world)
                .sub(before[1].angular_velocity_world)
                .sub(after[1].angular_velocity_world)
                .scale(0.5)
                .dot(axes[2]);
            let common_spin = before[0]
                .angular_velocity_world
                .add(after[0].angular_velocity_world)
                .add(before[1].angular_velocity_world)
                .add(after[1].angular_velocity_world)
                .scale(0.25)
                .dot(axes[2]);
            let slip = [relative.dot(axes[0]), relative.dot(axes[1])];
            let contact = self.step_resolved_candidate(
                thermal,
                target_s,
                frame,
                PartialSlipKinematics {
                    creepage: slip.map(|v| v / creepage_reference_speed_mps),
                    rolling_speed_mps: creepage_reference_speed_mps,
                    torsional_spin_rad_per_s: relative_spin,
                    contact_spin_rad_per_s: common_spin,
                    dt_s: dt,
                },
                ownership,
                controls.contact,
                energy_guess_j,
                interface_temperature_guess_k,
                &resolve,
                || (*cancel.borrow_mut())(),
            )?;
            contact_trials.set(
                contact_trials
                    .get()
                    .checked_add(contact.report().caloric.contact_trials)
                    .ok_or(RigidCaloricContactError::MotionBudget)?,
            );
            let law = &contact.report().caloric.contact;
            let residual = [
                x[0] - dt * law.tangent_force_n[0],
                x[1] - dt * law.tangent_force_n[1],
                spin_coordinate - dt * law.torsional_moment_nm / length,
            ];
            if residual.iter().any(|v| !v.is_finite()) {
                return Err(RigidCaloricContactError::Invalid(
                    "unrepresentable contact impulse residual",
                ));
            }
            Ok((next, contact, residual, slip, relative_spin))
        };
        let problem = FallibleFiniteDifference::new(
            dimension,
            controls.difference_step_n_s,
            |x: &[f64]| evaluate(x).map(|v| v.2[..dimension].to_vec()),
        )
        .map_err(rigid_difference_error)?;
        let construction = NewtonKrylovState::new(
            &problem,
            impulse_guess_n_s[..dimension].to_vec(),
            controls.newton,
        );
        if let Some(error) = problem.take_failure() {
            return Err(rigid_difference_error(error));
        }
        let mut solver = construction?;
        let result =
            solver.run_controlled(&problem, controls.maximum_newton_iterations, &mut || {
                (*cancel.borrow_mut())() || problem.has_failed()
            });
        if let Some(error) = problem.take_failure() {
            return Err(rigid_difference_error(error));
        }
        poll()?;
        if result.progress != SolverRunProgress::Complete {
            return Err(RigidCaloricContactError::Interrupted(result.progress));
        }
        if !result.report.converged {
            return Err(RigidCaloricContactError::NotConverged(result.report));
        }
        let (next, contact, _, slip, spin) = evaluate(&solver.x)?;
        let law = &contact.report().caloric.contact;
        let impulse = axes[0].scale(solver.x[0]).add(axes[1].scale(solver.x[1]));
        let couple = axes[2].scale(length * solver.x.get(2).copied().unwrap_or(0.0));
        let expected_impulse = vector(law.wrench.force_n).scale(dt);
        let expected_couple = vector(law.wrench.torque_nm).scale(dt);
        let mut linear = [Vec3::ZERO; 2];
        let mut angular = [Vec3::ZERO; 2];
        let mut dp = Vec3::ZERO;
        let mut dl = Vec3::ZERO;
        let mut kinetic = 0.0;
        for i in 0..2 {
            let sign = if i == 0 { 1.0 } else { -1.0 };
            let actual = next[i]
                .linear_momentum_world()
                .sub(bodies[i].linear_momentum_world());
            let rotation = bodies[i].pose().orientation();
            let actual_angular = rotation.rotate_body_to_world(
                next[i]
                    .angular_momentum_body()
                    .sub(bodies[i].angular_momentum_body()),
            );
            linear[i] = actual.sub(expected_impulse.scale(sign));
            angular[i] = actual_angular.sub(
                rotation
                    .rotate_body_to_world(arms_body_m[i])
                    .cross(expected_impulse)
                    .add(expected_couple)
                    .scale(sign),
            );
            dp = dp.add(actual);
            dl = dl
                .add(actual_angular)
                .add(bodies[i].pose().position_world().cross(actual));
            kinetic +=
                next[i].kinetic_energy(properties[i])? - bodies[i].kinetic_energy(properties[i])?;
        }
        let caloric = &contact.report().caloric;
        let mechanical = kinetic + caloric.contact.generalized_work.work_into_interface_j;
        let total = kinetic
            + caloric.contact_internal_energy_change_j
            + caloric.thermal.stored_energy_change_j
            + caloric.contact.generalized_work.numerical_dissipation_j
            - caloric.external_heat_j;
        if linear
            .into_iter()
            .chain([dp])
            .any(|v| !v.is_finite() || norm(v) > controls.linear_momentum_tolerance_n_s)
            || angular
                .into_iter()
                .chain([dl])
                .any(|v| !v.is_finite() || norm(v) > controls.angular_momentum_tolerance_nm_s)
            || [mechanical, total]
                .iter()
                .any(|v| !v.is_finite() || v.abs() > controls.energy_tolerance_j)
        {
            return Err(RigidCaloricContactError::Invalid(
                "native rigid momentum, work or combined thermal balance failed",
            ));
        }
        poll()?;
        let (thermal, contact, report) = contact.into_parts();
        Ok(RigidCaloricContactCandidate {
            bodies: next,
            thermal,
            contact,
            report: RigidCaloricContactReport {
                contact: report,
                nonlinear: result.report,
                motion_trials: trials.get(),
                contact_trials: contact_trials.get(),
                applied_impulse_on_a_n_s: impulse,
                applied_spin_impulse_on_a_nm_s: couple,
                midpoint_relative_velocity_mps: slip,
                midpoint_relative_spin_rad_s: spin,
                body_momentum_residual_n_s: linear,
                body_angular_momentum_residual_nm_s: angular,
                linear_momentum_change_n_s: dp,
                world_angular_momentum_change_nm_s: dl,
                kinetic_energy_change_j: kinetic,
                mechanical_energy_residual_j: mechanical,
                total_energy_residual_j: total,
            },
        })
    }

    fn admit_rigid_geometry<E>(
        &self,
        bodies: &[RigidBodyState; 2],
        properties: [MassProperties; 2],
        arms: [Vec3; 2],
        thermal: &DepositionSimulation<'_>,
        frame: TangentFrame,
        tolerance: f64,
    ) -> Result<(), RigidCaloricContactError<E>> {
        let route = self
            .resolved_heat_route()
            .ok_or(RigidCaloricContactError::Invalid(
                "resolved thermal interface required",
            ))?;
        let patch = self.checkpoint.normal_patch();
        let area =
            core::f64::consts::PI * patch.semi_axis_longitudinal_m() * patch.semi_axis_lateral_m();
        if !area.is_finite() || (area - route.contact_area_m2()).abs() > 128.0 * f64::EPSILON * area
        {
            return Err(RigidCaloricContactError::Invalid(
                "native ellipse and thermal face area differ",
            ));
        }
        let cells = route.cells();
        let a = thermal
            .cells()
            .get(&cells[0])
            .ok_or(RigidCaloricContactError::Invalid(
                "inactive thermal contact cell",
            ))?
            .geometry;
        let b = thermal
            .cells()
            .get(&cells[1])
            .ok_or(RigidCaloricContactError::Invalid(
                "inactive thermal contact cell",
            ))?
            .geometry;
        let centre = [a.start_m, a.end_m]
            .into_iter()
            .find(|p| *p == b.start_m || *p == b.end_m)
            .ok_or(RigidCaloricContactError::Invalid(
                "thermal cells have no shared face",
            ))?;
        let axis = vector(a.end_m).sub(vector(a.start_m));
        if norm(axis.cross(vector(frame.normal()))) > tolerance {
            return Err(RigidCaloricContactError::Invalid(
                "native normal and thermal cell centreline differ",
            ));
        }
        for i in 0..2 {
            let point = bodies[i]
                .point_kinematics(properties[i], arms[i])?
                .point_world;
            if norm(point.sub(vector(centre))) > tolerance {
                return Err(RigidCaloricContactError::Invalid(
                    "native material points must meet the thermal face",
                ));
            }
        }
        Ok(())
    }
}
fn vector(v: [f64; 3]) -> Vec3 {
    Vec3::new(v[0], v[1], v[2])
}
fn norm(v: Vec3) -> f64 {
    v.x.hypot(v.y).hypot(v.z)
}
fn rigid_difference_error<E>(
    error: FiniteDifferenceError<RigidCaloricContactError<E>>,
) -> RigidCaloricContactError<E> {
    match error {
        FiniteDifferenceError::Source(source) => source,
        FiniteDifferenceError::Invalid(message) => RigidCaloricContactError::Invalid(message),
    }
}
