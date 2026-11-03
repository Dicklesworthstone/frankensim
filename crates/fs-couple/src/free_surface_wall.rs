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

/// Fixed-center finite-inertia rotation driven by the original liquid links.
pub mod rotor;
/// Joint tangential translation and axial spin through the same liquid links.
pub mod planar;

use fs_exec::CancelGate;
use fs_lbm::d3q19::freesurface3::PlanarWall3;
use fs_lbm::material::MaterialLbmError;
use fs_lbm::material::free_surface::{MaterialFreeSurface3, MaterialPlanarWallMotionReport3};
use fs_mbd::{DynamicsError, MassProperties, MomentumImpulseReceipt, Pose, RigidBodyState, Vec3};

/// Explicit root, physical-balance and complete native-trial allowances.
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
    /// The admitted speed bracket does not enclose the coupled momentum root.
    Unbracketed {
        /// Coupled impulse residual at the lower speed [kg m/s].
        lower_residual_kg_m_s: f64,
        /// Coupled impulse residual at the upper speed [kg m/s].
        upper_residual_kg_m_s: f64,
    },
    /// No strictly interior representable speed remains before root acceptance.
    Unresolved {
        /// Native trials already evaluated.
        native_trials: usize,
        /// Smaller-magnitude endpoint impulse defect [kg m/s].
        residual_kg_m_s: f64,
    },
    /// Actual native endpoint impulse, shared velocity or work failed admission.
    Balance {
        /// Represented wall increment minus guided and external impulses.
        momentum_residual_kg_m_s: Vec3,
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
    /// Unmodified original SI liquid wall-exchange receipt.
    pub wall: MaterialPlanarWallMotionReport3,
    /// Fixed-pose native impulse receipt. `motion.body()` additionally includes
    /// midpoint translation; orientation and angular momentum are unchanged.
    pub momentum: MomentumImpulseReceipt,
    /// Held guide drive integrated over the native physical timestep.
    pub external_impulse_world: Vec3,
    /// Translational guide impulse ON the wall in the native numerical pressure
    /// datum. This closes raw momentum; it is not an absolute physical load.
    pub raw_guide_impulse_world: Vec3,
    /// Guide reaction inferred from the native pressure-perturbation load.
    /// This diagnostic is NOT substituted into raw population momentum.
    pub pressure_perturbation_guide_impulse_world: Vec3,
    /// Physical normal/transverse guide impulse when an absolute pressure
    /// reference was explicitly bound to the liquid. No reference is invented.
    pub absolute_guide_impulse_world: Option<Vec3>,
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
    wall: MaterialPlanarWallMotionReport3,
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
    /// Native link arms are unavailable, so the rotational guide reaction is
    /// unresolved and no total-angular-momentum claim is made.
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
    /// Solve 2m(w-v0) = J_wall(w)·t + F dt. Bracketed secant trials use the
    /// complete native wall response; a midpoint fallback guarantees bracket
    /// contraction. No empirical drag/added mass or force-frame replacement is
    /// supplied. Native trial failures refuse the whole interval. A user may
    /// provide a narrower physical bracket or refine the original liquid dt.
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
        let lower = self.trial(
            controls.midpoint_speed_bracket_m_s[0],
            external_force_n,
            controls,
            &mut trials,
            &mut trial_cells,
            gate,
        )?;
        if converged(&lower, controls)? {
            return self.finish(lower, external_force_n, controls, trials, trial_cells, gate);
        }
        let mut lo = (lower.speed, lower.residual);
        drop(lower);
        let upper = self.trial(
            controls.midpoint_speed_bracket_m_s[1],
            external_force_n,
            controls,
            &mut trials,
            &mut trial_cells,
            gate,
        )?;
        if converged(&upper, controls)? {
            return self.finish(upper, external_force_n, controls, trials, trial_cells, gate);
        }
        let mut hi = (upper.speed, upper.residual);
        drop(upper);
        if lo.1 >= 0.0 || hi.1 <= 0.0 {
            return Err(GuidedWallError::Unbracketed {
                lower_residual_kg_m_s: lo.1,
                upper_residual_kg_m_s: hi.1,
            });
        }
        loop {
            // Scale before division; forming hi.residual-lo.residual can overflow.
            let scale = (-lo.1).max(hi.1);
            let a = -lo.1 / scale;
            let b = hi.1 / scale;
            let fraction = a / (a + b);
            let fraction = if (0.1..=0.9).contains(&fraction) {
                fraction
            } else {
                0.5
            };
            let speed = lo.0 * (1.0 - fraction) + hi.0 * fraction;
            if !speed.is_finite() || speed <= lo.0 || speed >= hi.0 {
                return Err(GuidedWallError::Unresolved {
                    native_trials: trials,
                    residual_kg_m_s: if -lo.1 < hi.1 { lo.1 } else { hi.1 },
                });
            }
            let candidate = self.trial(
                speed,
                external_force_n,
                controls,
                &mut trials,
                &mut trial_cells,
                gate,
            )?;
            if converged(&candidate, controls)? {
                return self.finish(
                    candidate,
                    external_force_n,
                    controls,
                    trials,
                    trial_cells,
                    gate,
                );
            }
            if candidate.residual < 0.0 {
                lo = (speed, candidate.residual);
            } else {
                hi = (speed, candidate.residual);
            }
        }
    }

    fn admit(&self, force: f64, c: GuidedWallControls) -> Result<(), GuidedWallError> {
        let [lo, hi] = c.midpoint_speed_bracket_m_s;
        let max_speed = self
            .liquid
            .scaling()
            .velocity_si(self.liquid.scaling().flow().u_lattice)?;
        if !force.is_finite()
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
        let (liquid, wall) = self.liquid.planar_wall_motion_candidate_si(
            self.wall,
            array(self.tangent.scale(speed)),
            c.maximum_relative_mass_defect,
            gate,
        )?;
        poll(gate)?;
        let initial = self
            .body
            .center_of_mass_velocity_world(self.properties)?
            .dot(self.tangent);
        let change = finite((self.properties.mass() * (speed - initial)) * 2.0)?;
        let impulse = finite(vector(wall.wall_impulse_kg_m_s).dot(self.tangent))?;
        let external = force * self.liquid.scaling().timestep_s();
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
        let imposed = vector(trial.wall.velocity_m_s);
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
        if !momentum_residual.is_finite()
            || !midpoint_residual.is_finite()
            || norm(momentum_residual) > impulse_limit
            || norm(midpoint_residual) > midpoint_limit
            || work_residual.abs() > work_limit
        {
            return Err(GuidedWallError::Balance {
                momentum_residual_kg_m_s: momentum_residual,
                midpoint_residual_m_s: midpoint_residual,
                work_residual_j: work_residual,
            });
        }
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
        if !guide.is_finite()
            || !perturbation_guide.is_finite()
            || absolute_guide.is_some_and(|value| !value.is_finite())
        {
            return Err(invalid("guide reaction is unrepresentable"));
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
            wall: trial.wall,
            momentum,
            external_impulse_world: external,
            raw_guide_impulse_world: guide,
            pressure_perturbation_guide_impulse_world: perturbation_guide,
            absolute_guide_impulse_world: absolute_guide,
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

