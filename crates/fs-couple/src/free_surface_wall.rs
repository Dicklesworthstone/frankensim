//! Finite wall inertia coupled to the original material free-surface liquid.
//!
//! The original guided owner translates a flooded wall along one tangent. The
//! native liquid supplies its actual reflected-population impulse and work;
//! the native rigid owner supplies the wall's representable momentum change.
//! Their shared midpoint wall velocity is solved before either state escapes.
//! Its wall cannot spin; the rotor and planar submodules add axial and joint
//! tangential/axial motion. All retain a fixed plane without cover/vacate,
//! normal motion, finite plate edges or deformable geometry.
//!
//! The wall's kinetic work is checked against native liquid boundary work.
//! This is not an energy-conserving liquid discretization or a heat model:
//! original isothermal VOF viscosity, forcing, capillarity and conversion retain
//! their own numerical/physical scope. Timestep and spatial refinement remain
//! necessary, particularly for a light wall. All trials start from the same
//! accepted fluid; trial history never becomes physical elapsed time.

/// Joint tangential translation and axial spin through the same liquid links.
pub mod planar;
/// Fixed-center finite-inertia rotation driven by the original liquid links.
pub mod rotor;

use fs_exec::CancelGate;
use fs_lbm::d3q19::freesurface3::PlanarWall3;
use fs_lbm::material::MaterialLbmError;
use fs_lbm::material::free_surface::{
    MaterialFreeSurface3, MaterialPlanarRigidWallMotion3, MaterialPlanarRigidWallReport3,
    MaterialPlanarWallMotionReport3,
};
use fs_math::dd::Dd;
use fs_mbd::{DynamicsError, MassProperties, MomentumImpulseReceipt, Pose, RigidBodyState, Vec3};

/// Explicit midpoint domain, physical-balance and complete native-trial allowances.
#[derive(Clone, Copy, Debug)]
pub struct GuidedWallControls {
    /// Finite ordered bracket for the signed MIDPOINT wall speed [m/s].
    /// Both bounds must satisfy the original material scenario speed limit.
    /// Native liquid admission also applies at every attempted bracket point.
    pub midpoint_speed_bracket_m_s: [f64; 2],
    /// Counts every native liquid trial, including bracket and rejected trials.
    pub maximum_native_trials: usize,
    /// Sum of complete domain cell counts over all native trials. This is a
    /// conservative grid-work measure, not a count of scalar operations.
    pub maximum_trial_cells: usize,
    /// Original native liquid relative mass gate on each attempted timestep.
    pub maximum_relative_mass_defect: f64,
    /// Absolute scalar-root and represented vector-impulse allowance [kg m/s].
    pub momentum_absolute_kg_m_s: f64,
    /// Full orbital-plus-spin impulse allowance about `balance_origin_m` [kg m²/s].
    pub angular_momentum_absolute_kg_m2_s: f64,
    /// Fixed inertial origin for the complete wall angular-momentum balance [m].
    pub balance_origin_m: Vec3,
    /// Absolute mismatch of native wall and rigid midpoint velocities [m/s].
    pub midpoint_absolute_m_s: f64,
    /// Absolute mechanical boundary-work defect allowance [J].
    pub work_absolute_j: f64,
    /// Relative allowance on increments/work/speeds, with no unit-sized floor
    /// and no scaling by the large numerical normal-pressure impulse; in [0, 1).
    pub relative_tolerance: f64,
}

/// A refused trial never changes either original physical owner.
#[derive(Clone, Debug, PartialEq)]
pub enum GuidedWallError {
    /// Invalid physical/budget data or unrepresentable arithmetic.
    Invalid(&'static str),
    /// Original liquid source, scaling or native population refusal.
    Fluid(MaterialLbmError),
    /// Original rigid state/inertia/impulse refusal.
    Mechanics(DynamicsError),
    /// The next native trial would exceed the explicit solve ceiling.
    TrialBudget {
        /// Native trials already attempted successfully before this refusal.
        completed: usize,
        /// Caller-supplied maximum number of native trials.
        maximum: usize,
    },
    /// The next native trial would exceed its accumulated full-domain work cap.
    WorkBudget {
        /// Accumulated trial cells including the refused next trial.
        required_trial_cells: usize,
        /// Caller-supplied accumulated trial-cell ceiling.
        maximum: usize,
    },
    /// The admitted speed interval does not contain the coupled momentum root.
    Unbracketed {
        /// Native affine prediction of the impulse residual at the lower speed [kg m/s].
        lower_residual_kg_m_s: f64,
        /// Native affine prediction of the impulse residual at the upper speed [kg m/s].
        upper_residual_kg_m_s: f64,
    },
    /// The final native replay did not meet the scalar midpoint balance.
    Unresolved {
        /// Native trials already evaluated.
        native_trials: usize,
        /// Actual final native impulse defect [kg m/s].
        residual_kg_m_s: f64,
    },
    /// Actual native endpoint impulse, shared velocity or work failed admission.
    Balance {
        /// Represented wall increment minus guided and external impulses.
        momentum_residual_kg_m_s: Vec3,
        /// Represented orbital-plus-spin increment minus all applied moments.
        angular_momentum_residual_kg_m2_s: Vec3,
        /// Actual rigid midpoint minus imposed native liquid wall velocity.
        midpoint_residual_m_s: Vec3,
        /// Mechanical boundary-work defect, never assigned to heat.
        work_residual_j: f64,
    },
    /// Cancellation observed before any complete candidate escaped.
    Cancelled,
}
impl core::fmt::Display for GuidedWallError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "guided liquid wall refused: {self:?}")
    }
}
impl std::error::Error for GuidedWallError {}
impl From<MaterialLbmError> for GuidedWallError {
    fn from(error: MaterialLbmError) -> Self {
        Self::Fluid(error)
    }
}
impl From<DynamicsError> for GuidedWallError {
    fn from(error: DynamicsError) -> Self {
        Self::Mechanics(error)
    }
}

/// Complete liquid populations and guided wall state at one accepted endpoint.
/// Clones are independent physical candidates, with shared immutable sources.
#[derive(Clone)]
pub struct GuidedFreeSurfaceWall<'a> {
    liquid: MaterialFreeSurface3<'a>,
    properties: MassProperties,
    body: RigidBodyState,
    wall: PlanarWall3,
    tangent: Vec3,
    domain_cells: usize,
    time_s: f64,
    accepted_steps: u64,
}

/// Both native endpoints plus their actual momentum and boundary-work accounts.
pub struct GuidedWallStep<'a> {
    /// Continue subsequent intervals from this complete accepted candidate.
    pub motion: GuidedFreeSurfaceWall<'a>,
    /// Translational view of the original SI liquid wall-exchange receipt.
    pub wall: MaterialPlanarWallMotionReport3,
    /// Complete native link force/moment receipt at the free half-drift COM.
    /// Zero spin retains the original uniform-velocity population arithmetic.
    pub wrench: MaterialPlanarRigidWallReport3,
    /// Fixed-pose native impulse receipt. `motion.body()` additionally includes
    /// midpoint translation; orientation and angular momentum are unchanged.
    pub momentum: MomentumImpulseReceipt,
    /// Held guide drive integrated over the native physical timestep.
    pub external_impulse_world: Vec3,
    /// Translational guide impulse ON the wall in the native numerical pressure
    /// datum. This closes raw momentum; it is not an absolute physical load.
    pub raw_guide_impulse_world: Vec3,
    /// Guide's pure angular impulse about `wrench.motion.origin_m` [kg m²/s].
    /// This cancels the actual native liquid moment and keeps wall spin zero.
    pub raw_guide_angular_impulse_world: Vec3,
    /// Guide reaction inferred from the native pressure-perturbation load.
    /// This diagnostic is NOT substituted into raw population momentum.
    pub pressure_perturbation_guide_impulse_world: Vec3,
    /// Guide couple from the perturbation traction at the same COM [kg m²/s].
    pub pressure_perturbation_guide_angular_impulse_world: Vec3,
    /// Physical normal/transverse guide impulse when an absolute pressure
    /// reference was explicitly bound to the liquid. No reference is invented.
    pub absolute_guide_impulse_world: Option<Vec3>,
    /// Guide couple including an explicitly bound absolute pressure [kg m²/s].
    pub absolute_guide_angular_impulse_world: Option<Vec3>,
    /// Inertial origin of the following complete angular account [m].
    pub balance_origin_m: Vec3,
    /// Actual wall orbital-plus-spin momentum increment about that origin [kg m²/s].
    pub angular_momentum_change_world: Vec3,
    /// Actual angular increment minus liquid, guide and applied-force moments.
    pub angular_momentum_residual_kg_m2_s: Vec3,
    /// Actual velocity imposed in the accepted native liquid trial.
    pub midpoint_velocity_world: Vec3,
    /// Represented rigid midpoint minus that imposed liquid wall velocity.
    pub midpoint_residual_m_s: Vec3,
    /// Actual wall momentum increment minus guided wall and external impulses.
    pub momentum_residual_kg_m_s: Vec3,
    /// Held external force work through the actual rigid midpoint velocity.
    pub external_work_j: f64,
    /// Ideally zero for a stationary guide; retained floating residual work.
    pub raw_guide_work_j: f64,
    /// Wall kinetic change + native liquid wall work - external/guide work.
    /// Computed from impulse and interface defects to retain small residuals.
    pub mechanical_work_residual_j: f64,
    /// Actual full native trials, including the bracketing endpoints.
    pub native_trials: usize,
    /// Sum of original full-domain cell counts across all native trials.
    pub trial_cells: usize,
}

struct Trial<'a> {
    liquid: MaterialFreeSurface3<'a>,
    wall: MaterialPlanarRigidWallReport3,
    speed: f64,
    residual: f64,
    scale: f64,
}

impl<'a> GuidedFreeSurfaceWall<'a> {
    /// Bind a finite wall to an existing frozen-material liquid trajectory.
    ///
    /// `tangent_world` must be unit length within roundoff and exactly tangent
    /// to `wall`; it is normalized once. The native body must have zero spin
    /// and translate only along that guide, up to arithmetic roundoff. Its pose
    /// tracks material translation along the unchanged planar domain wall.
    /// The original liquid link arms supply the guide's full reaction couple;
    /// orbital and spin momentum remain separate from pressure-datum diagnostics.
    ///
    /// `maximum_cells` bounds the complete domain, including gas and solid
    /// cells, before any trial clone or native plane scan. Time below is elapsed
    /// since this binding; the liquid's original step count is retained.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        liquid: MaterialFreeSurface3<'a>,
        properties: MassProperties,
        body: RigidBodyState,
        wall: PlanarWall3,
        tangent_world: Vec3,
        maximum_cells: usize,
        gate: &CancelGate,
    ) -> Result<Self, GuidedWallError> {
        poll(gate)?;
        let length = norm(tangent_world);
        if !tangent_world.is_finite()
            || (length - 1.0).abs() > 64.0 * f64::EPSILON
            || component(tangent_world, normal_axis(wall)) != 0.0
            || body.angular_momentum_body() != Vec3::ZERO
        {
            return Err(invalid("unit tangent guide and zero body spin required"));
        }
        let tangent = tangent_world.scale(length.recip());
        let velocity = body.center_of_mass_velocity_world(properties)?;
        if norm(velocity.sub(tangent.scale(velocity.dot(tangent))))
            > 128.0 * f64::EPSILON * norm(velocity)
        {
            return Err(invalid("initial wall velocity must lie along its guide"));
        }
        let domain_cells = liquid
            .lattice()
            .domain_cells()
            .into_iter()
            .try_fold(1_usize, |n, x| n.checked_mul(x))
            .ok_or_else(|| invalid("complete liquid domain cell count overflow"))?;
        if domain_cells == 0 || domain_cells > maximum_cells {
            return Err(invalid(
                "liquid domain exceeds the guided-wall cell ceiling",
            ));
        }
        let max_speed = liquid
            .scaling()
            .velocity_si(liquid.scaling().flow().u_lattice)?;
        if norm(velocity) > max_speed {
            return Err(invalid(
                "initial wall speed exceeds the original liquid scenario",
            ));
        }
        poll(gate)?;
        Ok(Self {
            liquid,
            properties,
            body,
            wall,
            tangent,
            domain_cells,
            time_s: 0.0,
            accepted_steps: 0,
        })
    }

    /// Accepted original material source, native populations and liquid history.
    pub const fn liquid(&self) -> &MaterialFreeSurface3<'a> {
        &self.liquid
    }
    /// Accepted native rigid momentum and material translation pose.
    pub const fn body(&self) -> RigidBodyState {
        self.body
    }
    /// Original fixed wall mass and inertia, without an added-mass replacement.
    pub const fn properties(&self) -> MassProperties {
        self.properties
    }
    /// Admitted fixed unit direction of the translational guide.
    pub const fn tangent_world(&self) -> Vec3 {
        self.tangent
    }
    /// Accepted physical time elapsed since this wall/liquid binding [s].
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    /// Accepted physical intervals, excluding all rejected native trials.
    pub const fn accepted_steps(&self) -> u64 {
        self.accepted_steps
    }

    /// Advance one ORIGINAL liquid timestep with a held signed guide force [N].
    ///
    /// Solve 2m(w-v0) = J_wall(w)·t + F dt using the original reflection's affine
    /// response. One pilot supplies its exact scalar resistance; one replay at
    /// the proposed midpoint supplies the accepted populations and full wrench.
    /// An already balanced pilot needs no replay. The explicit speed interval
    /// still bounds admission. No empirical drag/added mass or second fluid
    /// algorithm is supplied; every native failure refuses the whole interval.
    ///
    /// Successful liquid and body candidates escape together only after actual
    /// representable impulse, shared midpoint, work, pose and cancellation
    /// checks. Failed iterations do not advance time or consume physical mass.
    pub fn step_candidate(
        &self,
        external_force_n: f64,
        controls: GuidedWallControls,
        gate: &CancelGate,
    ) -> Result<GuidedWallStep<'a>, GuidedWallError> {
        poll(gate)?;
        self.admit(external_force_n, controls)?;
        let mut trials = 0;
        let mut trial_cells = 0;
        let initial = self
            .body
            .center_of_mass_velocity_world(self.properties)?
            .dot(self.tangent);
        let [lower, upper] = controls.midpoint_speed_bracket_m_s;
        let pilot = self.trial(
            initial.clamp(lower, upper),
            external_force_n,
            controls,
            &mut trials,
            &mut trial_cells,
            gate,
        )?;
        if converged(&pilot, controls)? {
            return self.finish(pilot, external_force_n, controls, trials, trial_cells, gate);
        }
        let axis = normal_axis(self.wall);
        let tangents = [
            component(self.tangent, (axis + 1) % 3),
            component(self.tangent, (axis + 2) % 3),
        ];
        let mut resistance = Dd::ZERO;
        for i in 0..2 {
            for j in 0..2 {
                resistance = resistance
                    + Dd::from_f64(tangents[i])
                        * Dd::from_f64(pilot.wall.fluid_response[i][j])
                        * Dd::from_f64(tangents[j]);
            }
        }
        let resistance = finite(resistance.to_f64())?;
        if resistance < 0.0 {
            return Err(invalid("native guided wall resistance must be nonnegative"));
        }
        let half_slope =
            Dd::from_f64(self.properties.mass()) + Dd::from_f64(0.5) * Dd::from_f64(resistance);
        let speed = finite(
            (Dd::from_f64(pilot.speed)
                - Dd::from_f64(0.5) * Dd::from_f64(pilot.residual) / half_slope)
                .to_f64(),
        )?;
        if speed < lower || speed > upper {
            let residual_at = |value| {
                finite(
                    (Dd::from_f64(pilot.residual)
                        + Dd::from_f64(2.0)
                            * half_slope
                            * (Dd::from_f64(value) - Dd::from_f64(pilot.speed)))
                    .to_f64(),
                )
            };
            return Err(GuidedWallError::Unbracketed {
                lower_residual_kg_m_s: residual_at(lower)?,
                upper_residual_kg_m_s: residual_at(upper)?,
            });
        }
        // The immutable pilot is discarded before allocating the accepted
        // replay. Trials never coexist as additional physical liquid histories.
        drop(pilot);
        let candidate = self.trial(
            speed,
            external_force_n,
            controls,
            &mut trials,
            &mut trial_cells,
            gate,
        )?;
        if !converged(&candidate, controls)? {
            return Err(GuidedWallError::Unresolved {
                native_trials: trials,
                residual_kg_m_s: candidate.residual,
            });
        }
        self.finish(
            candidate,
            external_force_n,
            controls,
            trials,
            trial_cells,
            gate,
        )
    }

    fn admit(&self, force: f64, c: GuidedWallControls) -> Result<(), GuidedWallError> {
        let [lo, hi] = c.midpoint_speed_bracket_m_s;
        let max_speed = self
            .liquid
            .scaling()
            .velocity_si(self.liquid.scaling().flow().u_lattice)?;
        if !force.is_finite()
            || !c.balance_origin_m.is_finite()
            || !lo.is_finite()
            || !hi.is_finite()
            || lo >= hi
            || lo.abs().max(hi.abs()) > max_speed
            || c.maximum_native_trials == 0
            || c.maximum_trial_cells == 0
            || !c.maximum_relative_mass_defect.is_finite()
            || c.maximum_relative_mass_defect < 0.0
            || !c.relative_tolerance.is_finite()
            || c.relative_tolerance < 0.0
            || c.relative_tolerance >= 1.0
            || [
                c.momentum_absolute_kg_m_s,
                c.angular_momentum_absolute_kg_m2_s,
                c.midpoint_absolute_m_s,
                c.work_absolute_j,
            ]
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.0)
        {
            return Err(invalid(
                "finite speed bracket, force, physical tolerances and positive trial budgets required",
            ));
        }
        finite(force * self.liquid.scaling().timestep_s())?;
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn trial(
        &self,
        speed: f64,
        force: f64,
        c: GuidedWallControls,
        trials: &mut usize,
        trial_cells: &mut usize,
        gate: &CancelGate,
    ) -> Result<Trial<'a>, GuidedWallError> {
        poll(gate)?;
        if *trials >= c.maximum_native_trials {
            return Err(GuidedWallError::TrialBudget {
                completed: *trials,
                maximum: c.maximum_native_trials,
            });
        }
        let required = trial_cells
            .checked_add(self.domain_cells)
            .ok_or_else(|| invalid("native trial grid-work overflow"))?;
        if required > c.maximum_trial_cells {
            return Err(GuidedWallError::WorkBudget {
                required_trial_cells: required,
                maximum: c.maximum_trial_cells,
            });
        }
        *trial_cells = required;
        *trials += 1;
        let before_velocity = self.body.center_of_mass_velocity_world(self.properties)?;
        let dt = self.liquid.scaling().timestep_s();
        // A kick at the free half-drift COM followed by the other half drift
        // produces the same midpoint translation used by `finish`. Moments
        // therefore share the actual impulse configuration, even off-center.
        let origin = self
            .body
            .pose()
            .position_world()
            .add(before_velocity.scale(0.5 * dt));
        let (liquid, wall) = self.liquid.planar_rigid_wall_motion_candidate_si(
            MaterialPlanarRigidWallMotion3 {
                wall: self.wall,
                origin_m: array(origin),
                velocity_at_origin_m_s: array(self.tangent.scale(speed)),
                angular_speed_rad_s: 0.0,
            },
            c.maximum_relative_mass_defect,
            gate,
        )?;
        poll(gate)?;
        let initial = before_velocity.dot(self.tangent);
        let change = finite((self.properties.mass() * (speed - initial)) * 2.0)?;
        let impulse = finite(vector(wall.wall_impulse_kg_m_s).dot(self.tangent))?;
        let external = force * dt;
        let residual = finite((change - impulse) - external)?;
        let scale = change.abs().max(impulse.abs()).max(external.abs());
        Ok(Trial {
            liquid,
            wall,
            speed,
            residual,
            scale,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn finish(
        &self,
        trial: Trial<'a>,
        force: f64,
        c: GuidedWallControls,
        native_trials: usize,
        trial_cells: usize,
        gate: &CancelGate,
    ) -> Result<GuidedWallStep<'a>, GuidedWallError> {
        poll(gate)?;
        let dt = self.liquid.scaling().timestep_s();
        let wall_impulse = vector(trial.wall.wall_impulse_kg_m_s);
        let guided = self.tangent.scale(wall_impulse.dot(self.tangent));
        let external = self.tangent.scale(force * dt);
        // Apply only the resultant permitted impulse, avoiding cancellation of
        // a large numerical normal pressure against its guide reaction.
        let momentum = self.body.apply_momentum_impulse(
            self.properties,
            self.tangent
                .scale(wall_impulse.dot(self.tangent) + force * dt),
            Vec3::ZERO,
        )?;
        let before_v = self.body.center_of_mass_velocity_world(self.properties)?;
        let after_v = momentum
            .state_after
            .center_of_mass_velocity_world(self.properties)?;
        let midpoint = before_v.scale(0.5).add(after_v.scale(0.5));
        let imposed = vector(trial.wall.motion.velocity_at_origin_m_s);
        let midpoint_residual = midpoint.sub(imposed);
        let guide = guided.sub(wall_impulse);
        let momentum_residual = momentum.linear_increment_world.sub(guided).sub(external);
        let external_work = finite(midpoint.dot(external))?;
        let guide_work = finite(midpoint.dot(guide))?;
        // Algebraically K1-K0 + W_fluid - W_external - W_guide.
        // Keep the actual impulse defect and velocity mismatch before summing
        // so a large common kinetic baseline cannot hide a small lost impulse.
        let work_residual = finite(
            midpoint.dot(momentum_residual)
                + midpoint_residual.dot(wall_impulse)
                + (imposed.dot(wall_impulse) - trial.wall.wall_work_j),
        )?;
        let impulse_scale = norm(momentum.linear_increment_world)
            .max(norm(guided))
            .max(norm(external));
        let velocity_scale = norm(before_v).max(norm(after_v)).max(norm(imposed));
        let work_scale = momentum
            .work
            .kinetic_energy_change_j
            .abs()
            .max(external_work.abs())
            .max(guide_work.abs())
            .max(trial.wall.wall_work_j.abs());
        let impulse_limit = allowance(
            c.momentum_absolute_kg_m_s,
            c.relative_tolerance,
            impulse_scale,
        )?;
        let midpoint_limit = allowance(
            c.midpoint_absolute_m_s,
            c.relative_tolerance,
            velocity_scale,
        )?;
        let work_limit = allowance(c.work_absolute_j, c.relative_tolerance, work_scale)?;
        let max_speed = self
            .liquid
            .scaling()
            .velocity_si(self.liquid.scaling().flow().u_lattice)?;
        if norm(after_v) > max_speed {
            return Err(invalid(
                "accepted wall endpoint exceeds the original liquid speed bound",
            ));
        }
        let displacement = midpoint.scale(dt);
        let old_pose = self.body.pose();
        let pose = Pose::new(
            old_pose.position_world().add(displacement),
            old_pose.orientation(),
        )?;
        let represented = pose.position_world().sub(old_pose.position_world());
        let displacement_limit = finite(
            dt * allowance(
                c.midpoint_absolute_m_s,
                c.relative_tolerance,
                norm(midpoint),
            )?,
        )?;
        if !displacement.is_finite() || norm(represented.sub(displacement)) > displacement_limit {
            return Err(invalid(
                "guided wall displacement is not representable within its velocity tolerance",
            ));
        }
        let body = RigidBodyState::new(
            pose,
            momentum.state_after.linear_momentum_world(),
            momentum.state_after.angular_momentum_body(),
        )?;
        let guide_from = |load: Vec3| self.tangent.scale(load.dot(self.tangent)).sub(load);
        let perturbation_guide =
            guide_from(vector(trial.wall.pressure_perturbation_wall_force_n).scale(dt));
        let absolute_guide = trial
            .wall
            .absolute_wall_force_n
            .map(|load| guide_from(vector(load).scale(dt)));
        let guide_couple = vector(trial.wall.wall_angular_impulse_kg_m2_s).scale(-1.0);
        let perturbation_couple =
            vector(trial.wall.pressure_perturbation_wall_torque_nm).scale(-dt);
        let absolute_couple = trial
            .wall
            .absolute_wall_torque_nm
            .map(|moment| vector(moment).scale(-dt));
        if !guide.is_finite()
            || !perturbation_guide.is_finite()
            || absolute_guide.is_some_and(|value| !value.is_finite())
            || !guide_couple.is_finite()
            || !perturbation_couple.is_finite()
            || absolute_couple.is_some_and(|value| !value.is_finite())
        {
            return Err(invalid("guide reaction is unrepresentable"));
        }
        // Difference actual stored endpoints before forming the orbital moment.
        // delta(x cross p) = (x0-o) cross delta_p + delta_x cross p1.
        // Spin stays exactly zero. Liquid and guide couples at the free
        // half-drift COM cancel, while their net force retains its orbital arm.
        let difference = |after: Vec3, before: Vec3| {
            core::array::from_fn(|d| Dd::from_f64(array(after)[d]) - Dd::from_f64(array(before)[d]))
        };
        let p1 = body.linear_momentum_world();
        let delta_p = difference(p1, self.body.linear_momentum_world());
        let initial_arm = difference(old_pose.position_world(), c.balance_origin_m);
        let travel = difference(pose.position_world(), old_pose.position_world());
        let initial_orbital = cross_dd(initial_arm, delta_p);
        let transport_orbital = cross_dd(travel, array(p1).map(Dd::from_f64));
        let angular_change: [Dd; 3] =
            core::array::from_fn(|d| initial_orbital[d] + transport_orbital[d]);
        let net_impulse = core::array::from_fn(|d| {
            Dd::from_f64(array(guided)[d]) + Dd::from_f64(array(external)[d])
        });
        let impulse_arm = difference(vector(trial.wall.motion.origin_m), c.balance_origin_m);
        let angular_input = cross_dd(impulse_arm, net_impulse);
        let angular_residual = vector(core::array::from_fn(|d| {
            (angular_change[d] - angular_input[d]).to_f64()
        }));
        let angular_change = vector(angular_change.map(Dd::to_f64));
        let angular_limit = allowance(
            c.angular_momentum_absolute_kg_m2_s,
            c.relative_tolerance,
            norm(angular_change).max(norm(vector(angular_input.map(Dd::to_f64)))),
        )?;
        if !momentum_residual.is_finite()
            || !angular_change.is_finite()
            || !angular_residual.is_finite()
            || !midpoint_residual.is_finite()
            || norm(momentum_residual) > impulse_limit
            || norm(angular_residual) > angular_limit
            || norm(midpoint_residual) > midpoint_limit
            || work_residual.abs() > work_limit
        {
            return Err(GuidedWallError::Balance {
                momentum_residual_kg_m_s: momentum_residual,
                angular_momentum_residual_kg_m2_s: angular_residual,
                midpoint_residual_m_s: midpoint_residual,
                work_residual_j: work_residual,
            });
        }
        let time_s = finite(self.time_s + dt)?;
        if time_s <= self.time_s {
            return Err(invalid("guided wall clock cannot advance"));
        }
        let accepted_steps = self
            .accepted_steps
            .checked_add(1)
            .ok_or_else(|| invalid("guided wall accepted-step counter overflow"))?;
        let motion = Self {
            liquid: trial.liquid,
            properties: self.properties,
            body,
            wall: self.wall,
            tangent: self.tangent,
            domain_cells: self.domain_cells,
            time_s,
            accepted_steps,
        };
        poll(gate)?;
        Ok(GuidedWallStep {
            motion,
            wall: translational_receipt(trial.wall),
            wrench: trial.wall,
            momentum,
            external_impulse_world: external,
            raw_guide_impulse_world: guide,
            raw_guide_angular_impulse_world: guide_couple,
            pressure_perturbation_guide_impulse_world: perturbation_guide,
            pressure_perturbation_guide_angular_impulse_world: perturbation_couple,
            absolute_guide_impulse_world: absolute_guide,
            absolute_guide_angular_impulse_world: absolute_couple,
            balance_origin_m: c.balance_origin_m,
            angular_momentum_change_world: angular_change,
            angular_momentum_residual_kg_m2_s: angular_residual,
            midpoint_velocity_world: imposed,
            midpoint_residual_m_s: midpoint_residual,
            momentum_residual_kg_m_s: momentum_residual,
            external_work_j: external_work,
            raw_guide_work_j: guide_work,
            mechanical_work_residual_j: work_residual,
            native_trials,
            trial_cells,
        })
    }
}

// A view of one native receipt, not another boundary solve or force model.
fn translational_receipt(
    wrench: MaterialPlanarRigidWallReport3,
) -> MaterialPlanarWallMotionReport3 {
    MaterialPlanarWallMotionReport3 {
        wall: wrench.motion.wall,
        velocity_m_s: wrench.motion.velocity_at_origin_m_s,
        fluid_impulse_kg_m_s: wrench.fluid_impulse_kg_m_s,
        wall_impulse_kg_m_s: wrench.wall_impulse_kg_m_s,
        reference_pressure_fluid_impulse_kg_m_s: wrench.reference_pressure_fluid_impulse_kg_m_s,
        pressure_perturbation_wall_force_n: wrench.pressure_perturbation_wall_force_n,
        absolute_wall_force_n: wrench.absolute_wall_force_n,
        fluid_work_j: wrench.fluid_work_j,
        wall_work_j: wrench.wall_work_j,
        wetted_area_m2: wrench.wetted_area_m2,
        mass_before_kg: wrench.mass_before_kg,
        mass_after_kg: wrench.mass_after_kg,
        relative_mass_defect: wrench.relative_mass_defect,
        maximum_speed_m_s: wrench.maximum_speed_m_s,
    }
}

fn cross_dd(a: [Dd; 3], b: [Dd; 3]) -> [Dd; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn converged(trial: &Trial<'_>, controls: GuidedWallControls) -> Result<bool, GuidedWallError> {
    Ok(trial.residual.abs()
        <= allowance(
            controls.momentum_absolute_kg_m_s,
            controls.relative_tolerance,
            trial.scale,
        )?)
}
fn allowance(absolute: f64, relative: f64, scale: f64) -> Result<f64, GuidedWallError> {
    finite(absolute + relative * scale)
}
fn invalid(reason: &'static str) -> GuidedWallError {
    GuidedWallError::Invalid(reason)
}
fn finite(value: f64) -> Result<f64, GuidedWallError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("unrepresentable guided wall arithmetic"))
    }
}
fn poll(gate: &CancelGate) -> Result<(), GuidedWallError> {
    if gate.is_requested() {
        Err(GuidedWallError::Cancelled)
    } else {
        Ok(())
    }
}
fn vector(value: [f64; 3]) -> Vec3 {
    Vec3::new(value[0], value[1], value[2])
}
fn array(value: Vec3) -> [f64; 3] {
    [value.x, value.y, value.z]
}
fn norm(value: Vec3) -> f64 {
    value.x.hypot(value.y).hypot(value.z)
}
fn component(value: Vec3, axis: usize) -> f64 {
    array(value)[axis]
}
fn normal_axis(wall: PlanarWall3) -> usize {
    match wall {
        PlanarWall3::XMin | PlanarWall3::XMax => 0,
        PlanarWall3::YMin | PlanarWall3::YMax => 1,
        PlanarWall3::ZMin | PlanarWall3::ZMax => 2,
    }
}
