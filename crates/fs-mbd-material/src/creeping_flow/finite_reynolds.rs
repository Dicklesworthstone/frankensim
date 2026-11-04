//! Finite-Reynolds drag on the original source-bound, nonrotating sphere.
//!
//! Schiller--Naumann gives F = -6*pi*mu*a*(1 + 0.15*Re^0.687)*(v-U),
//! with Re = 2*a*|v-U|/nu. The resistance form has a regular zero-speed
//! Stokes limit; no divergent drag coefficient or artificial speed floor is
//! evaluated. See van Wachem et al., JFM 1018 (2025), A41, section 3.2,
//! equations (3.5)-(3.6): <https://doi.org/10.1017/jfm.2025.10526>.
//!
//! This empirical isolated-sphere mean-drag correlation is limited to Re <= 1000 here.
//! Ambient flow is the UNDISTURBED uniform reservoir velocity, not a local
//! particle-disturbed CFD sample. No spin, lift, wake interaction, history,
//! added mass, wall, rarefaction or compressibility correction is inferred.
//! Continuum, Newtonian and incompressible applicability remains a caller/source
//! obligation: the shared transport snapshot alone cannot enforce Mach or Knudsen.
//! The fluid/body timescale gate is an approximation screen, not a bound on
//! the neglected transient flow or the empirical correlation error.

use super::*;

const EXPONENT: f64 = 0.687;

/// Explicit empirical and quasi-steady applicability limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteReynoldsLimits {
    /// Diameter-based relative Reynolds ceiling, in (0,1000].
    pub maximum_reynolds: f64,
    /// In (0,0.1]. Bounds (a²/nu)*c_tangent/m at every accepted stage.
    /// Heavy bodies in light fluids are the intended quasi-steady regime.
    pub maximum_diffusion_ratio: f64,
}

/// Shared source, spherical geometry and native mass owner with nonlinear drag.
#[derive(Debug)]
pub struct FiniteReynoldsSphere<'a, F: FluidTransport + ?Sized> {
    // Reuse the original source/geometry preparation, never its Stokes-only
    // regime or stepping methods. Its rotary resistance is not used here.
    sphere: CreepingSphere<'a, F>,
    maximum_reynolds: f64,
}

/// Instantaneous finite-Reynolds drag and buoyancy, excluding body weight.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteReynoldsLoad {
    /// Native force at the COM and zero torque; this lane refuses spin.
    pub wrench: Wrench,
    /// Resistance opposite the body/fluid relative velocity [N].
    pub drag_force_world_n: Vec3,
    /// Fully immersed Archimedes force [N].
    pub buoyancy_force_world_n: Vec3,
    /// Nonnegative mechanical loss in the prescribed fluid reservoir [W].
    pub dissipated_power_w: f64,
    /// Power supplied by maintaining the background flow, F_drag dot U [W].
    pub background_power_w: f64,
    /// Diameter-based Reynolds number at this relative velocity.
    pub reynolds: f64,
    /// Largest eigenvalue of the drag/relative-velocity tangent [kg/s].
    pub tangent_resistance_kg_s: f64,
}

/// Atomic native candidate and separate physical/numerical accounts.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteReynoldsStep {
    /// Original integrator's accepted state and mechanical diagnostics.
    pub dynamics: StepReceipt,
    /// Fluid load at the supplied starting state.
    pub load_before: FiniteReynoldsLoad,
    /// Load supplied to the native integrator by the nonlinear midpoint solve.
    pub load_midpoint: FiniteReynoldsLoad,
    /// Fluid load at the actual native endpoint.
    pub load_after: FiniteReynoldsLoad,
    /// Drag loss deposited in the prescribed fluid reservoir [J].
    pub viscous_heat_j: f64,
    /// Work supplied by maintaining the uniform background flow [J].
    pub background_work_j: f64,
    /// Work supplied by the separately held external COM force [J].
    pub applied_work_j: f64,
    /// Kinetic + body-gravity + displaced-fluid potential change [J].
    pub stored_energy_change_j: f64,
    /// Storage + heat - background work - applied work [J]; never physical heat.
    pub energy_defect_j: f64,
    /// Largest fluid-diffusion/body-relaxation ratio of the three stages.
    pub diffusion_ratio: f64,
    /// Difference between solved and actual native midpoint velocities [m/s].
    pub midpoint_velocity_residual_m_s: f64,
    /// Number of scalar residual/tangent evaluations in the bounded solve.
    pub midpoint_iterations: usize,
}

impl<'a, F: FluidTransport + ?Sized> FiniteReynoldsSphere<'a, F> {
    /// Declare the same single, uniform, fully immersed no-slip sphere used by
    /// the Stokes owner, with an explicit finite-Reynolds empirical extension.
    /// Density, viscosity, displaced mass and body inertia retain their source
    /// owners. No body/medium state evolves during this preparation.
    pub fn new(
        source: &'a F,
        radius_m: f64,
        gravity: Gravity,
        ambient_velocity_world: Vec3,
        limits: FiniteReynoldsLimits,
    ) -> Result<Self, Error> {
        if !small_positive(limits.maximum_reynolds, 1000.0) {
            return Err(Error::Invalid(
                "sphere drag Reynolds ceiling must be in (0,1000]",
            ));
        }
        let sphere = CreepingSphere::new(
            source,
            radius_m,
            gravity,
            ambient_velocity_world,
            CreepingFlowLimits {
                maximum_reynolds: 0.1,
                maximum_diffusion_ratio: limits.maximum_diffusion_ratio,
            },
        )?;
        Ok(Self {
            sphere,
            maximum_reynolds: limits.maximum_reynolds,
        })
    }

    /// Original fluid owner, including its material/mixture source receipts.
    pub const fn source(&self) -> &'a F {
        self.sphere.source()
    }
    /// The single frozen transport snapshot shared by force and regime checks.
    pub const fn transport(&self) -> FluidTransportState {
        self.sphere.transport()
    }
    /// Declared spherical boundary radius [m].
    pub const fn radius_m(&self) -> f64 {
        self.sphere.radius_m()
    }
    /// Same gravity used for body weight and displaced-fluid buoyancy.
    pub const fn gravity(&self) -> Gravity {
        self.sphere.gravity()
    }

    /// Force at the supplied state, without advancing the native body.
    /// Nonzero spin refuses rather than silently omitting rotary drag or lift.
    pub fn load(
        &self,
        body: &PrincipalBody<'_>,
        state: RigidBodyState,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FiniteReynoldsLoad, Error> {
        poll(&mut cancelled)?;
        self.sphere.admit_body(body)?;
        let result = self.load_frame(body.frame(), state)?;
        poll(&mut cancelled)?;
        Ok(result)
    }

    fn load_frame(
        &self,
        frame: &PrincipalFrame,
        state: RigidBodyState,
    ) -> Result<FiniteReynoldsLoad, Error> {
        if state.angular_momentum_body() != Vec3::ZERO {
            return Err(Error::Invalid(
                "finite-Reynolds sphere drag requires zero spin",
            ));
        }
        let relative = vector(
            state
                .center_of_mass_velocity_world(frame.properties())?
                .sub(self.sphere.ambient_velocity_world),
        )?;
        let (reynolds, correction) = self.correction(norm(relative)?)?;
        let drag = vector(relative.scale(-self.sphere.resistance_kg_s * (1.0 + correction)))?;
        let buoyancy = vector(
            self.sphere
                .gravity
                .acceleration_world()
                .scale(-self.sphere.fluid_mass_kg),
        )?;
        let dissipated_power_w = finite(-drag.dot(relative))?;
        let background_power_w = finite(drag.dot(self.sphere.ambient_velocity_world))?;
        if dissipated_power_w < 0.0 || relative != Vec3::ZERO && dissipated_power_w == 0.0 {
            return Err(Error::Invalid(
                "relative motion needs representable positive drag dissipation",
            ));
        }
        Ok(FiniteReynoldsLoad {
            wrench: Wrench {
                force_world: vector(drag.add(buoyancy))?,
                torque_body: Vec3::ZERO,
            },
            drag_force_world_n: drag,
            buoyancy_force_world_n: buoyancy,
            dissipated_power_w,
            background_power_w,
            reynolds,
            tangent_resistance_kg_s: positive(
                self.sphere.resistance_kg_s * (1.0 + (1.0 + EXPONENT) * correction),
            )?,
        })
    }

    fn correction(&self, speed: f64) -> Result<(f64, f64), Error> {
        let reynolds = self
            .sphere
            .transport
            .reynolds(speed, 2.0 * self.sphere.radius_m)?;
        if reynolds > self.maximum_reynolds {
            return Err(Error::Reynolds {
                translation: reynolds,
                rotation: 0.0,
                maximum: self.maximum_reynolds,
            });
        }
        Ok((reynolds, 0.15 * fs_math::det::pow(reynolds, EXPONENT)))
    }

    /// Close nonlinear drag at the original integrator's midpoint velocity.
    ///
    /// Radial symmetry reduces the velocity solve to one monotone scalar:
    /// s*(1 + dt*c(s)/(2*m)) = |v-U + dt*g_buoyant/2|. Safeguarded Newton
    /// uses its analytic tangent, bounded by 64 iterations without allocation.
    /// The native rigid-body integrator remains the only state update owner.
    ///
    /// The existing controls now bound dt*c_tangent/m, at both endpoints and
    /// the midpoint. Keeping that fraction <=1 prevents drag-only sign reversal.
    /// This is a second-order discrete model; time refinement and correlation
    /// validity remain necessary even when its work balance closes to roundoff.
    pub fn step(
        &self,
        body: &PrincipalBody<'_>,
        state: RigidBodyState,
        controls: CreepingStepControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<FiniteReynoldsStep, Error> {
        self.step_with_wrench(body, state, Wrench::ZERO, controls, cancelled)
    }

    /// Solve drag and a held external COM force at the SAME native midpoint.
    /// The world-frame force excludes gravity and fluid loading, which this
    /// owner already supplies. Nonzero applied torque is refused because this
    /// correlation admits nonrotating spheres only. An offset-point force
    /// cannot be treated as a COM force unless its torque is independently zero.
    ///
    /// The radial nonlinear solve is unchanged: its predictor now includes
    /// the imposed acceleration. External drive work is reported separately
    /// and subtracted in the work-energy balance, never deposited as heat.
    /// A rapidly changing drive requires time refinement and an appropriate
    /// unsteady fluid model; passing the quasi-steady gate is not such proof.
    pub fn step_with_wrench(
        &self,
        body: &PrincipalBody<'_>,
        state: RigidBodyState,
        applied: Wrench,
        controls: CreepingStepControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FiniteReynoldsStep, Error> {
        poll(&mut cancelled)?;
        self.sphere.admit_body(body)?;
        vector(applied.force_world)?;
        vector(applied.torque_body)?;
        if applied.torque_body != Vec3::ZERO {
            return Err(Error::Invalid(
                "finite-Reynolds sphere drag requires zero applied torque",
            ));
        }
        if !small_positive(controls.maximum_relaxation_fraction, 1.0)
            || !controls.duration_s.is_finite()
            || controls.duration_s <= 0.0
            || !controls.maximum_energy_defect_j.is_finite()
            || controls.maximum_energy_defect_j <= 0.0
        {
            return Err(Error::Invalid(
                "positive finite step/work controls and relaxation fraction in (0,1] required",
            ));
        }
        let frame = body.frame();
        let properties = body.properties();
        let load_before = self.load_frame(frame, state)?;
        let mut diffusion_ratio =
            self.check_timescales(properties.mass(), load_before, controls)?;
        let dt = controls.duration_s;
        let half_dt = positive(0.5 * dt)?;
        let velocity = state.center_of_mass_velocity_world(properties)?;
        let relative_before = vector(velocity.sub(self.sphere.ambient_velocity_world))?;
        let acceleration = vector(
            self.sphere
                .gravity
                .acceleration_world()
                .scale(1.0 - self.sphere.fluid_mass_kg / properties.mass())
                .add(applied.force_world.scale(1.0 / properties.mass())),
        )?;
        let rhs = vector(relative_before.add(acceleration.scale(half_dt)))?;
        let rhs_speed = norm(rhs)?;
        let factor = positive(half_dt * (self.sphere.resistance_kg_s / properties.mass()))?;
        let (mid_speed, midpoint_iterations) =
            self.solve_midpoint(rhs_speed, factor, &mut cancelled)?;
        let relative_mid = if rhs_speed == 0.0 {
            Vec3::ZERO
        } else {
            vector(rhs.scale(mid_speed / rhs_speed))?
        };
        let midpoint = RigidBodyState::new(
            state.pose(),
            vector(
                relative_mid
                    .add(self.sphere.ambient_velocity_world)
                    .scale(properties.mass()),
            )?,
            Vec3::ZERO,
        )?;
        let load_midpoint = self.load_frame(frame, midpoint)?;
        diffusion_ratio = diffusion_ratio.max(self.check_timescales(
            properties.mass(),
            load_midpoint,
            controls,
        )?);
        let total = Wrench {
            force_world: vector(load_midpoint.wrench.force_world.add(applied.force_world))?,
            torque_body: Vec3::ZERO,
        };
        let dynamics =
            RigidBodyIntegrator::new(self.sphere.gravity).step(state, properties, total, dt)?;
        let load_after = self.load_frame(frame, dynamics.state_after)?;
        diffusion_ratio =
            diffusion_ratio.max(self.check_timescales(properties.mass(), load_after, controls)?);
        let actual_mid = vector(
            velocity.scale(0.5).add(
                dynamics
                    .state_after
                    .center_of_mass_velocity_world(properties)?
                    .scale(0.5),
            ),
        )?;
        let solved_mid = midpoint.center_of_mass_velocity_world(properties)?;
        let midpoint_velocity_residual_m_s = norm(actual_mid.sub(solved_mid))?;
        let viscous_heat_j = finite(dt * load_midpoint.dissipated_power_w)?;
        let background_work_j = finite(dt * load_midpoint.background_power_w)?;
        let applied_power_w = finite(applied.force_world.dot(actual_mid))?;
        let applied_work_j = finite(dt * applied_power_w)?;
        if load_midpoint.dissipated_power_w > 0.0 && viscous_heat_j == 0.0
            || load_midpoint.background_power_w != 0.0 && background_work_j == 0.0
            || applied_power_w != 0.0 && applied_work_j == 0.0
        {
            return Err(Error::Invalid(
                "sphere drag heat or work vanished in conversion to energy",
            ));
        }
        let dx = dynamics
            .state_after
            .pose()
            .position_world()
            .sub(state.pose().position_world());
        let stored_energy_change_j = finite(
            dynamics.diagnostics_after.translational_kinetic_energy
                - dynamics.diagnostics_before.translational_kinetic_energy
                - (properties.mass() - self.sphere.fluid_mass_kg)
                    * self.sphere.gravity.acceleration_world().dot(dx),
        )?;
        let energy_defect_j =
            finite(stored_energy_change_j + viscous_heat_j - background_work_j - applied_work_j)?;
        if energy_defect_j.abs() > controls.maximum_energy_defect_j {
            return Err(Error::EnergyDefect {
                defect_j: energy_defect_j,
                maximum_j: controls.maximum_energy_defect_j,
            });
        }
        poll(&mut cancelled)?;
        Ok(FiniteReynoldsStep {
            dynamics,
            load_before,
            load_midpoint,
            load_after,
            viscous_heat_j,
            background_work_j,
            applied_work_j,
            stored_energy_change_j,
            energy_defect_j,
            diffusion_ratio,
            midpoint_velocity_residual_m_s,
            midpoint_iterations,
        })
    }

    fn check_timescales(
        &self,
        mass: f64,
        load: FiniteReynoldsLoad,
        controls: CreepingStepControls,
    ) -> Result<f64, Error> {
        let relaxation_rate = positive(load.tangent_resistance_kg_s / mass)?;
        let diffusion_ratio = positive(self.sphere.diffusion_time_s * relaxation_rate)?;
        if diffusion_ratio > self.sphere.limits.maximum_diffusion_ratio {
            return Err(Error::Unsteady {
                diffusion_ratio,
                maximum: self.sphere.limits.maximum_diffusion_ratio,
            });
        }
        if finite(controls.duration_s * relaxation_rate)? > controls.maximum_relaxation_fraction {
            return Err(Error::Invalid(
                "step exceeds nonlinear drag relaxation fraction; reduce dt",
            ));
        }
        Ok(diffusion_ratio)
    }

    fn solve_midpoint(
        &self,
        rhs: f64,
        factor: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(f64, usize), Error> {
        if rhs == 0.0 {
            return Ok((0.0, 0));
        }
        // A root outside the law's domain is refused; never evaluate the
        // empirical correction at a nonphysical high-Re predictor.
        let maximum_speed = positive(
            (self.maximum_reynolds * self.sphere.transport.kinematic_viscosity_m2_s())
                / (2.0 * self.sphere.radius_m),
        )?;
        let mut lower = 0.0;
        let mut upper = (rhs / (1.0 + factor)).min(maximum_speed);
        let mut speed = upper;
        for iteration in 1..=64 {
            poll(cancelled)?;
            let (_, correction) = self.correction(speed)?;
            let residual = finite(speed * (1.0 + factor * (1.0 + correction)) - rhs)?;
            if residual.abs() <= 16.0 * f64::EPSILON * rhs {
                return Ok((speed, iteration));
            }
            if iteration == 1 && residual < 0.0 {
                return Err(Error::Invalid(
                    "sphere drag midpoint exceeds the declared Reynolds ceiling",
                ));
            }
            if residual > 0.0 {
                upper = speed;
            } else {
                lower = speed;
            }
            let derivative = positive(1.0 + factor * (1.0 + (1.0 + EXPONENT) * correction))?;
            let newton = speed - residual / derivative;
            speed = if newton > lower && newton < upper {
                newton
            } else {
                lower + 0.5 * (upper - lower)
            };
        }
        Err(Error::Invalid(
            "sphere drag midpoint solve did not converge",
        ))
    }
}
