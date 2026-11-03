//! A finite axial rotor driven by the original resolved free-surface wall.
//!
//! The COM is fixed and one principal inertia axis is the stationary plane's
//! normal. The guide carries all translation and transverse torque; axial
//! torque changes the actual native body spin and attitude. Two original fluid
//! candidates suffice: a pilot measures the reflection's affine response, then
//! final replay uses the solved common midpoint angular speed. The pilot is
//! discarded and advances no physical time. This is fixed flooded geometry,
//! without tilting, cover/vacate, finite plate edges or caloric heat.

use super::{array, norm, normal_axis, vector};
use fs_exec::CancelGate;
use fs_lbm::d3q19::freesurface3::PlanarWall3;
use fs_lbm::material::MaterialLbmError;
use fs_lbm::material::free_surface::{
    MaterialFreeSurface3, MaterialPlanarRigidWallMotion3, MaterialPlanarRigidWallReport3,
};
use fs_math::dd::Dd;
use fs_mbd::{
    DynamicsError, Gravity, MassProperties, RigidBodyIntegrator, RigidBodyState, StepReceipt, Vec3,
    Wrench,
};

/// Explicit physical and full-grid work bounds for one original liquid tick.
#[derive(Clone, Copy, Debug)]
pub struct PlanarRotorControls {
    /// Sum of complete domain cell counts across the pilot and final replay.
    pub maximum_trial_cells: usize,
    /// Original native mass-conservation allowance on each fluid candidate.
    pub maximum_relative_mass_defect: f64,
    /// Maximum absolute midpoint angular displacement in this step [rad].
    pub maximum_rotation_rad: f64,
    /// Maximum absolute mechanical work exchanged with the wall [J].
    pub maximum_absolute_wall_work_j: f64,
    /// Absolute actual angular-impulse and affine-replay allowance [kg m²/s].
    pub angular_momentum_absolute_kg_m2_s: f64,
    /// Absolute imposed/actual midpoint angular-speed mismatch [rad/s].
    pub midpoint_absolute_rad_s: f64,
    /// Absolute rigid kinetic/boundary/external work defect [J].
    pub work_absolute_j: f64,
    /// Relative increment/work tolerance in [0,1), with no unit-sized floor
    /// or scaling by the large numerical normal-pressure force.
    pub relative_tolerance: f64,
}

/// A refused candidate publishes neither fluid nor native rigid motion.
#[derive(Clone, Debug, PartialEq)]
pub enum PlanarRotorError {
    /// Invalid physical state, unsupported guide, limit or finite arithmetic.
    Invalid(&'static str),
    /// Original material, scaling, local speed or native fluid refusal.
    Fluid(MaterialLbmError),
    /// Original native body/inertia/motion refusal.
    Mechanics(DynamicsError),
    /// Two full-domain native candidates exceed the explicit work bound.
    WorkBudget {
        required_trial_cells: usize,
        maximum: usize,
    },
    /// Actual represented endpoints fail the declared conservation/interface gates.
    Balance {
        angular_momentum_residual_kg_m2_s: Vec3,
        midpoint_residual_rad_s: Vec3,
        response_residual_kg_m2_s: f64,
        work_residual_j: f64,
    },
    /// Cancellation observed before final candidate publication.
    Cancelled,
}

impl core::fmt::Display for PlanarRotorError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "planar liquid rotor refused: {self:?}")
    }
}
impl std::error::Error for PlanarRotorError {}
impl From<MaterialLbmError> for PlanarRotorError {
    fn from(value: MaterialLbmError) -> Self {
        Self::Fluid(value)
    }
}
impl From<DynamicsError> for PlanarRotorError {
    fn from(value: DynamicsError) -> Self {
        Self::Mechanics(value)
    }
}

/// Original fluid and finite native rotor retained together for continuation.
#[derive(Clone)]
pub struct FreeSurfacePlanarRotor<'a> {
    liquid: MaterialFreeSurface3<'a>,
    properties: MassProperties,
    body: RigidBodyState,
    wall: PlanarWall3,
    axis_body: Vec3,
    inertia_kg_m2: f64,
    domain_cells: usize,
    elapsed_time_s: f64,
    accepted_steps: u64,
}

/// One accepted common fluid/rotor tick, with all guide reactions explicit.
pub struct PlanarRotorStep<'a> {
    /// Both original accepted owners; continue by stepping this motion.
    pub motion: FreeSurfacePlanarRotor<'a>,
    /// Final original reflection, source/scaling, actual wrench and work.
    pub wall: MaterialPlanarRigidWallReport3,
    /// Original MBD step, including its actual attitude and momentum endpoints.
    pub mechanics: StepReceipt,
    /// Applied axial actuator impulse about the fixed COM [kg m²/s].
    pub external_angular_impulse_world: Vec3,
    /// Stationary guide impulse cancelling all raw wall translation [kg m/s].
    /// Includes numerical pressure; this is a momentum account, not a physical
    /// absolute support load.
    pub raw_guide_impulse_world: Vec3,
    /// Guide angular impulse cancelling transverse raw wall moments [kg m²/s].
    pub raw_guide_angular_impulse_world: Vec3,
    /// Physical guide force impulse about the declared pressure reference.
    pub pressure_perturbation_guide_impulse_world: Vec3,
    /// Physical transverse guide angular impulse about that reference.
    pub pressure_perturbation_guide_angular_impulse_world: Vec3,
    /// Physical guide impulse with an explicitly supplied absolute reference.
    pub absolute_guide_impulse_world: Option<Vec3>,
    /// Physical transverse guide angular impulse with that same reference.
    pub absolute_guide_angular_impulse_world: Option<Vec3>,
    /// Angular speed imposed on every final fluid link [rad/s].
    pub midpoint_angular_speed_rad_s: f64,
    /// Actual native world midpoint angular velocity minus imposed velocity.
    pub midpoint_residual_rad_s: Vec3,
    /// Actual spatial spin increment minus wall, guide and actuator impulses.
    pub angular_momentum_residual_kg_m2_s: Vec3,
    /// Final native axial wall impulse minus measured affine prediction.
    pub response_residual_kg_m2_s: f64,
    /// Actual rotor kinetic increment from native momenta and midpoint pairing.
    pub kinetic_energy_change_j: f64,
    /// Actuator impulse paired with actual native midpoint angular velocity.
    pub external_work_j: f64,
    /// Actual guide torque paired with the native midpoint; ideally zero.
    pub raw_guide_work_j: f64,
    /// Kinetic change + liquid wall work - actuator and guide work [J].
    /// A numerical defect, never a heat destination or energy correction.
    pub mechanical_work_residual_j: f64,
    /// Exactly two original native candidates: pilot and accepted replay.
    pub native_trials: usize,
    /// Their total complete-domain cell work.
    pub trial_cells: usize,
}

impl<'a> FreeSurfacePlanarRotor<'a> {
    /// Bind a fixed-COM rotor whose normal is a principal axis of its ORIGINAL
    /// inertia. Initial spin must lie on that axis (to arithmetic roundoff)
    /// and initial linear momentum must be exactly zero. No state is projected.
    /// COM may lie off the wall plane; actual link arms retain that geometry.
    /// The full domain is bounded before any native clone or plane scan.
    pub fn new(
        liquid: MaterialFreeSurface3<'a>,
        properties: MassProperties,
        body: RigidBodyState,
        wall: PlanarWall3,
        maximum_cells: usize,
        gate: &CancelGate,
    ) -> Result<Self, PlanarRotorError> {
        poll(gate)?;
        body.kinetic_energy(properties)?;
        if body.linear_momentum_world() != Vec3::ZERO {
            return Err(invalid(
                "fixed rotor COM requires zero initial linear momentum",
            ));
        }
        let normal = axis_vector(normal_axis(wall));
        let axis_body = body.pose().orientation().rotate_world_to_body(normal);
        let unit_momentum = properties.angular_momentum_body_checked(axis_body)?;
        let inertia_kg_m2 = finite(axis_body.dot(unit_momentum))?;
        if inertia_kg_m2 <= 0.0
            || norm(unit_momentum.sub(axis_body.scale(inertia_kg_m2)))
                > 128.0 * f64::EPSILON * inertia_kg_m2
        {
            return Err(invalid(
                "wall normal must be a principal axis of the original rotor inertia",
            ));
        }
        let spin = body.angular_momentum_body();
        if norm(spin.sub(axis_body.scale(axis_body.dot(spin)))) > 128.0 * f64::EPSILON * norm(spin)
        {
            return Err(invalid(
                "initial rotor spin must lie on the fixed wall-normal axis",
            ));
        }
        let domain_cells = liquid
            .lattice()
            .domain_cells()
            .into_iter()
            .try_fold(1_usize, |n, x| n.checked_mul(x))
            .ok_or_else(|| invalid("rotor domain cell count overflowed"))?;
        if domain_cells == 0 || domain_cells > maximum_cells {
            return Err(invalid("liquid domain exceeds the rotor cell ceiling"));
        }
        let result = Self {
            liquid,
            properties,
            body,
            wall,
            axis_body,
            inertia_kg_m2,
            domain_cells,
            elapsed_time_s: 0.0,
            accepted_steps: 0,
        };
        result.admit_endpoint_speed(body)?;
        poll(gate)?;
        Ok(result)
    }

    /// Original accepted material liquid and population history.
    pub const fn liquid(&self) -> &MaterialFreeSurface3<'a> {
        &self.liquid
    }
    /// Original native rotor state, including actual attitude.
    pub const fn body(&self) -> RigidBodyState {
        self.body
    }
    /// Original inertia and mass; no apparent fluid inertia replaces them.
    pub const fn properties(&self) -> MassProperties {
        self.properties
    }
    /// Elapsed physical time since this binding, excluding pilot trials.
    pub const fn elapsed_time_s(&self) -> f64 {
        self.elapsed_time_s
    }
    /// Number of jointly accepted native fluid/body steps.
    pub const fn accepted_steps(&self) -> u64 {
        self.accepted_steps
    }
    /// Actual signed native angular speed about the fixed normal [rad/s].
    pub fn angular_speed_rad_s(&self) -> Result<f64, PlanarRotorError> {
        self.body_speed(self.body)
    }

    fn body_speed(&self, body: RigidBodyState) -> Result<f64, PlanarRotorError> {
        let omega = self
            .properties
            .angular_velocity_body_checked(body.angular_momentum_body())?;
        finite(omega.dot(self.axis_body))
    }

    fn prescribed(&self, angular_speed_rad_s: f64) -> MaterialPlanarRigidWallMotion3 {
        MaterialPlanarRigidWallMotion3 {
            wall: self.wall,
            origin_m: array(self.body.pose().position_world()),
            velocity_at_origin_m_s: [0.0; 3],
            angular_speed_rad_s,
        }
    }

    fn admit_endpoint_speed(&self, body: RigidBodyState) -> Result<(), PlanarRotorError> {
        let bound = self
            .liquid
            .scaling()
            .velocity_si(self.liquid.scaling().flow().u_lattice)?;
        let actual = self.prescribed(self.body_speed(body)?).maximum_speed_m_s(
            self.liquid.lattice().domain_cells(),
            self.liquid.scaling().cell_size_m(),
        )?;
        if actual > bound {
            return Err(invalid(
                "rotor endpoint exceeds the original material local-speed domain",
            ));
        }
        Ok(())
    }

    /// Advance a held axial actuator torque and reciprocal liquid traction.
    ///
    /// The pilot's measured response C solves
    /// `omega_mid = omega0 + (J_wall(omega0)+J_external)/(2 I + C)`.
    /// Final native replay supplies the actual torque to ONE original MBD step.
    /// Its real momentum, attitude, midpoint and kinetic increment must pass
    /// the declared gates. No body/fluid solver is replayed as a correction.
    /// The original fixed liquid timestep remains authoritative; users must
    /// refine it and the grid for trajectory accuracy, especially for light rotors.
    pub fn step_candidate(
        &self,
        external_torque_nm: f64,
        controls: PlanarRotorControls,
        gate: &CancelGate,
    ) -> Result<PlanarRotorStep<'a>, PlanarRotorError> {
        poll(gate)?;
        controls.validate()?;
        finite(external_torque_nm)?;
        let trial_cells = self
            .domain_cells
            .checked_mul(2)
            .ok_or_else(|| invalid("rotor native trial work overflowed"))?;
        if trial_cells > controls.maximum_trial_cells {
            return Err(PlanarRotorError::WorkBudget {
                required_trial_cells: trial_cells,
                maximum: controls.maximum_trial_cells,
            });
        }
        let axis = normal_axis(self.wall);
        let dt = self.liquid.scaling().timestep_s();
        let omega0 = self.body_speed(self.body)?;
        let external = finite(external_torque_nm * dt)?;
        let pilot = {
            let (candidate, report) = self.liquid.planar_rigid_wall_motion_candidate_si(
                self.prescribed(omega0),
                controls.maximum_relative_mass_defect,
                gate,
            )?;
            drop(candidate);
            report
        };
        poll(gate)?;
        let response = pilot.fluid_response[2][2];
        if !response.is_finite() || response < 0.0 {
            return Err(invalid(
                "native rotor angular response must be finite and nonnegative",
            ));
        }
        let denominator = finite(self.inertia_kg_m2 + 0.5 * response)?;
        let increment = (Dd::from_f64(pilot.wall_angular_impulse_kg_m2_s[axis])
            + Dd::from_f64(external))
            * Dd::from_f64(0.5)
            / Dd::from_f64(denominator);
        let imposed = finite((Dd::from_f64(omega0) + increment).to_f64())?;
        if finite(imposed * dt)?.abs() > controls.maximum_rotation_rad {
            return Err(invalid(
                "rotor midpoint rotation exceeds the declared step bound",
            ));
        }
        let (liquid, wall) = self.liquid.planar_rigid_wall_motion_candidate_si(
            self.prescribed(imposed),
            controls.maximum_relative_mass_defect,
            gate,
        )?;
        if wall.fluid_work_j.abs() > controls.maximum_absolute_wall_work_j {
            return Err(invalid(
                "rotor wall work exceeds the declared exchange bound",
            ));
        }
        poll(gate)?;
        let wall_axial = wall.wall_angular_impulse_kg_m2_s[axis];
        let total_impulse = Dd::from_f64(wall_axial) + Dd::from_f64(external);
        let torque = finite((total_impulse / Dd::from_f64(dt)).to_f64())?;
        let mechanics = RigidBodyIntegrator::new(Gravity::ZERO).step(
            self.body,
            self.properties,
            Wrench {
                force_world: Vec3::ZERO,
                torque_body: self.axis_body.scale(torque),
            },
            dt,
        )?;
        let body = mechanics.state_after;
        self.admit_endpoint_speed(body)?;
        if body.pose().position_world() != self.body.pose().position_world()
            || body.linear_momentum_world() != Vec3::ZERO
        {
            return Err(invalid("original rotor step changed the fixed COM"));
        }
        let actual_axis = body
            .pose()
            .orientation()
            .rotate_body_to_world(self.axis_body);
        if norm(actual_axis.sub(axis_vector(axis))) > 256.0 * f64::EPSILON {
            return Err(invalid(
                "native rotor attitude departed from its fixed normal axis",
            ));
        }
        let omega_before = self
            .properties
            .angular_velocity_body_checked(self.body.angular_momentum_body())?;
        let omega_after = self
            .properties
            .angular_velocity_body_checked(body.angular_momentum_body())?;
        let world_before = self
            .body
            .pose()
            .orientation()
            .rotate_body_to_world(omega_before);
        let world_after = body.pose().orientation().rotate_body_to_world(omega_after);
        let midpoint_dd: [Dd; 3] = core::array::from_fn(|d| {
            (Dd::from_f64(array(world_before)[d]) + Dd::from_f64(array(world_after)[d]))
                * Dd::from_f64(0.5)
        });
        let midpoint = vector(midpoint_dd.map(Dd::to_f64));
        if finite(norm(midpoint) * dt)?.abs() > controls.maximum_rotation_rad {
            return Err(invalid(
                "actual native rotor rotation exceeds the declared step bound",
            ));
        }
        let midpoint_residual = midpoint.sub(axis_vector(axis).scale(imposed));
        let before_spin = self
            .body
            .pose()
            .orientation()
            .rotate_body_to_world(self.body.angular_momentum_body());
        let after_spin = body
            .pose()
            .orientation()
            .rotate_body_to_world(body.angular_momentum_body());
        let angular_increment: [Dd; 3] = core::array::from_fn(|d| {
            Dd::from_f64(array(after_spin)[d]) - Dd::from_f64(array(before_spin)[d])
        });
        let angular_residual = vector(core::array::from_fn(|d| {
            (angular_increment[d] - if d == axis { total_impulse } else { Dd::ZERO }).to_f64()
        }));
        let response_residual = finite(
            (Dd::from_f64(wall_axial) - Dd::from_f64(pilot.wall_angular_impulse_kg_m2_s[axis])
                + Dd::from_f64(response) * (Dd::from_f64(imposed) - Dd::from_f64(omega0)))
            .to_f64(),
        )?;
        let mut kinetic = Dd::ZERO;
        for d in 0..3 {
            let average = (Dd::from_f64(array(omega_before)[d])
                + Dd::from_f64(array(omega_after)[d]))
                * Dd::from_f64(0.5);
            let change = Dd::from_f64(array(body.angular_momentum_body())[d])
                - Dd::from_f64(array(self.body.angular_momentum_body())[d]);
            kinetic = kinetic + average * change;
        }
        let raw_guide = vector(wall.wall_impulse_kg_m_s).scale(-1.0);
        let raw_guide_angular = transverse_reaction(wall.wall_angular_impulse_kg_m2_s, axis, 1.0);
        let external_work = midpoint_dd[axis] * Dd::from_f64(external);
        let guide_work = (0..3).fold(Dd::ZERO, |sum, d| {
            sum + midpoint_dd[d] * Dd::from_f64(array(raw_guide_angular)[d])
        });
        let work_residual = finite(
            (kinetic + Dd::from_f64(wall.fluid_work_j) - external_work - guide_work).to_f64(),
        )?;
        let kinetic_energy_change_j = finite(kinetic.to_f64())?;
        let external_work_j = finite(external_work.to_f64())?;
        let raw_guide_work_j = finite(guide_work.to_f64())?;
        let impulse_scale = wall_axial
            .abs()
            .max(external.abs())
            .max(angular_increment[axis].to_f64().abs());
        let speed_scale = imposed.abs().max(norm(midpoint));
        let work_scale = kinetic_energy_change_j
            .abs()
            .max(wall.fluid_work_j.abs())
            .max(external_work_j.abs())
            .max(raw_guide_work_j.abs());
        if !angular_residual.is_finite()
            || !midpoint_residual.is_finite()
            || norm(angular_residual)
                > allowance(
                    controls.angular_momentum_absolute_kg_m2_s,
                    controls.relative_tolerance,
                    impulse_scale,
                )?
            || response_residual.abs()
                > allowance(
                    controls.angular_momentum_absolute_kg_m2_s,
                    controls.relative_tolerance,
                    impulse_scale,
                )?
            || norm(midpoint_residual)
                > allowance(
                    controls.midpoint_absolute_rad_s,
                    controls.relative_tolerance,
                    speed_scale,
                )?
            || work_residual.abs()
                > allowance(
                    controls.work_absolute_j,
                    controls.relative_tolerance,
                    work_scale,
                )?
        {
            return Err(PlanarRotorError::Balance {
                angular_momentum_residual_kg_m2_s: angular_residual,
                midpoint_residual_rad_s: midpoint_residual,
                response_residual_kg_m2_s: response_residual,
                work_residual_j: work_residual,
            });
        }
        let time = finite(self.elapsed_time_s + dt)?;
        if time <= self.elapsed_time_s {
            return Err(invalid(
                "rotor elapsed time cannot represent its next native tick",
            ));
        }
        let accepted_steps = self
            .accepted_steps
            .checked_add(1)
            .ok_or_else(|| invalid("rotor accepted-step count overflowed"))?;
        let pressure_perturbation_guide_impulse_world =
            scaled_reaction(wall.pressure_perturbation_wall_force_n, dt)?;
        let pressure_perturbation_guide_angular_impulse_world = checked_vector(
            transverse_reaction(wall.pressure_perturbation_wall_torque_nm, axis, dt),
        )?;
        let absolute_guide_impulse_world = wall
            .absolute_wall_force_n
            .map(|force| scaled_reaction(force, dt))
            .transpose()?;
        let absolute_guide_angular_impulse_world = wall
            .absolute_wall_torque_nm
            .map(|moment| checked_vector(transverse_reaction(moment, axis, dt)))
            .transpose()?;
        poll(gate)?;
        Ok(PlanarRotorStep {
            motion: Self {
                liquid,
                properties: self.properties,
                body,
                wall: self.wall,
                axis_body: self.axis_body,
                inertia_kg_m2: self.inertia_kg_m2,
                domain_cells: self.domain_cells,
                elapsed_time_s: time,
                accepted_steps,
            },
            wall,
            mechanics,
            external_angular_impulse_world: axis_vector(axis).scale(external),
            raw_guide_impulse_world: raw_guide,
            raw_guide_angular_impulse_world: raw_guide_angular,
            pressure_perturbation_guide_impulse_world,
            pressure_perturbation_guide_angular_impulse_world,
            absolute_guide_impulse_world,
            absolute_guide_angular_impulse_world,
            midpoint_angular_speed_rad_s: imposed,
            midpoint_residual_rad_s: midpoint_residual,
            angular_momentum_residual_kg_m2_s: angular_residual,
            response_residual_kg_m2_s: response_residual,
            kinetic_energy_change_j,
            external_work_j,
            raw_guide_work_j,
            mechanical_work_residual_j: work_residual,
            native_trials: 2,
            trial_cells,
        })
    }
}

impl PlanarRotorControls {
    fn validate(self) -> Result<(), PlanarRotorError> {
        if [
            self.maximum_relative_mass_defect,
            self.maximum_absolute_wall_work_j,
            self.angular_momentum_absolute_kg_m2_s,
            self.midpoint_absolute_rad_s,
            self.work_absolute_j,
        ]
        .into_iter()
        .any(|value| !value.is_finite() || value < 0.0)
            || !self.maximum_rotation_rad.is_finite()
            || self.maximum_rotation_rad <= 0.0
            || !self.relative_tolerance.is_finite()
            || !(0.0..1.0).contains(&self.relative_tolerance)
        {
            return Err(invalid(
                "rotor bounds must be finite with nonnegative allowances and positive rotation limit",
            ));
        }
        Ok(())
    }
}

fn axis_vector(axis: usize) -> Vec3 {
    vector(core::array::from_fn(|d| if d == axis { 1.0 } else { 0.0 }))
}
fn transverse_reaction(moment: [f64; 3], axis: usize, scale: f64) -> Vec3 {
    vector(core::array::from_fn(|d| {
        if d == axis { 0.0 } else { -moment[d] * scale }
    }))
}
fn scaled_reaction(force: [f64; 3], dt: f64) -> Result<Vec3, PlanarRotorError> {
    checked_vector(vector(force).scale(-dt))
}
fn checked_vector(value: Vec3) -> Result<Vec3, PlanarRotorError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("unrepresentable rotor vector"))
    }
}
fn allowance(absolute: f64, relative: f64, scale: f64) -> Result<f64, PlanarRotorError> {
    finite(absolute + relative * scale)
}
fn invalid(reason: &'static str) -> PlanarRotorError {
    PlanarRotorError::Invalid(reason)
}
fn finite(value: f64) -> Result<f64, PlanarRotorError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("unrepresentable rotor arithmetic"))
    }
}
fn poll(gate: &CancelGate) -> Result<(), PlanarRotorError> {
    if gate.is_requested() {
        Err(PlanarRotorError::Cancelled)
    } else {
        Ok(())
    }
}
