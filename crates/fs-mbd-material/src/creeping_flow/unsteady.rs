//! Driven unsteady Stokes motion with reciprocal fluid inertia and Basset memory.
//!
//! For the translational channel of a no-slip sphere in quiescent unbounded Newtonian fluid,
//! Z(s) = 6*pi*mu*a*(1 + a*sqrt(s/nu)) + rho*V*s/2. The body's weight and
//! Archimedes force are separate from this force/velocity impedance. See
//! Candelier, Mehaddi, Mehlig & Magnaudet, JFM 954 (2023), A25, equation (1.1):
//! <https://doi.org/10.1017/jfm.2022.1015>.
//!
//! sqrt(s) = integral_0^infinity s/(s+lambda) d(lambda)/(pi*sqrt(lambda)).
//! With x=a*sqrt(lambda/nu), positive log-midpoint quadrature gives branches
//! c_j*s/(s+lambda_j), m_j=c_j/lambda_j. Each branch exchanges equal/opposite
//! force with the body, stores m_j*u_j^2/2 and dissipates c_j*(v-u_j)^2.
//! These are spectral fluid coordinates, not a reconstructed finite-volume
//! velocity field. The two omitted tails and quadrature need refinement.
//!
//! Midpoint elimination of the branch momenta supplies ONE wrench to the
//! existing rigid-body integrator. It remains the only owner of body momentum,
//! position and attitude; no alternate body trajectory is advanced here.
//! The explicit [`rotation`] owner adds the distinct causal rotary impedance
//! and one shared native Cayley step; the original owner remains nonrotating.
//!
//! A uniformly accelerating ambient flow adds its displaced-fluid pressure
//! force m_f*a_U separately from added mass -m_a*d(v-U)/dt. All diffusive
//! branches act on the slip v-U. This is the spatially uniform limit of Maxey
//! & Riley, Physics of Fluids 26 (1983), 883–889, DOI 10.1063/1.864230.

use super::*;
use std::sync::Arc;

/// Rotational fluid memory coupled to the same translating native body.
pub mod rotation;

mod prehistory;
pub use prehistory::SlipPrehistory;

/// Numerical approximation of the continuum history kernel, not material data.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BassetResolution {
    /// Smallest positive dimensionless diffusion wavenumber a*sqrt(lambda/nu).
    pub minimum_wavenumber_radius: f64,
    /// Largest wavenumber; finite and greater than the minimum.
    pub maximum_wavenumber_radius: f64,
    /// Logarithmic midpoint intervals in 1..=4096; bounds all history storage.
    pub intervals: usize,
}

/// One positive fluid-memory branch, fixed by source transport and geometry.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BassetBranch {
    pub resistance_kg_s: f64,
    pub mass_kg: f64,
    pub rate_per_s: f64,
}

/// Numerical acceptance, independent of physical dissipation and model error.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UnsteadyStepControls {
    pub duration_s: f64,
    pub maximum_energy_defect_j: f64,
    pub maximum_momentum_defect_ns: f64,
}

/// One accepted native body step and its separately represented fluid accounts.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UnsteadySphereStep {
    pub dynamics: StepReceipt,
    /// Includes buoyancy, ambient pressure, steady resistance, memory and added mass.
    /// With spin history, its body torque is the work-conjugate Cayley pullback
    /// of the physical world torque in `rotation`, not an instantaneous rotation.
    pub fluid_wrench: Wrench,
    pub history_force_world_n: Vec3,
    /// Reaction to acceleration RELATIVE to the uniform ambient flow.
    pub added_mass_force_world_n: Vec3,
    /// Displaced-fluid pressure force m_fluid*a_ambient, excluding buoyancy.
    pub ambient_pressure_force_world_n: Vec3,
    /// Steady Stokes loss plus all branch dashpot loss, never stored energy.
    pub viscous_heat_j: f64,
    /// Work of the separately supplied centre-of-mass force and optional
    /// body-frame couple, using the accepted native midpoint velocities.
    pub applied_work_j: f64,
    /// Work supplied by the prescribed translating/accelerating far field.
    /// This includes pressure-force work at the body and disturbance-reaction
    /// work at the ambient velocity; it is neither viscous heat nor storage.
    pub background_work_j: f64,
    /// Change of represented added-mass energy `m_added*|v-U|²/2`.
    pub added_mass_energy_change_j: f64,
    /// Sum of all represented translational and optional rotary fluid storage.
    pub memory_energy_change_j: f64,
    /// Added mass and branch momentum RELATIVE to the ambient flow, expressed
    /// in world axes. This is not the total momentum of a finite fluid field.
    pub spectral_fluid_momentum_change_ns: Vec3,
    /// Impulse delivered to the eliminated Stokes exterior by relative motion.
    pub far_field_impulse_ns: Vec3,
    pub momentum_residual_ns: Vec3,
    /// Dry-body kinetic + gravity/buoyancy potential + represented fluid energy.
    pub stored_energy_change_j: f64,
    pub energy_defect_j: f64,
    /// Present for a rotating owner; its heat/storage also enter the totals above.
    pub rotation: Option<rotation::RotationalHistoryStep>,
}

/// One restartable body and its complete accepted viscous-memory history.
/// The immutable source, original material mass and quadrature remain bound.
pub struct UnsteadySphereMotion<'s, 'f, 'b, F: FluidTransport + ?Sized> {
    sphere: &'s CreepingSphere<'f, F>,
    body: &'b PrincipalBody<'b>,
    resolution: BassetResolution,
    branches: Arc<[BassetBranch]>,
    fluid_velocities: Vec<Vec3>,
    rotation: Option<rotation::RotationHistory>,
    ambient_velocity_world: Vec3,
    state: RigidBodyState,
    time_s: f64,
}

impl<F: FluidTransport + ?Sized> Clone for UnsteadySphereMotion<'_, '_, '_, F> {
    fn clone(&self) -> Self {
        Self {
            sphere: self.sphere,
            body: self.body,
            resolution: self.resolution,
            branches: Arc::clone(&self.branches),
            fluid_velocities: self.fluid_velocities.clone(),
            rotation: self.rotation.clone(),
            ambient_velocity_world: self.ambient_velocity_world,
            state: self.state,
            time_s: self.time_s,
        }
    }
}

impl<'f, F: FluidTransport + ?Sized> CreepingSphere<'f, F> {
    /// Start the original sphere comoving with its uniform ambient fluid.
    ///
    /// This retains the original source snapshot, geometry, native gravity and
    /// Reynolds cap. The quasi-steady diffusion/relaxation restriction is replaced
    /// by actual added mass and viscous memory; it is NOT relaxed by density.
    /// Initial translation must be comoving and angular momentum must vanish.
    /// Author the initial native momentum as `mass * ambient_velocity`, with
    /// equality in momentum units. A velocity division can hide a one-ULP
    /// momentum mismatch or create apparent slip from exact co-motion, so
    /// rounded velocity equality alone does not admit the initial state.
    /// All relative memory starts at zero; arbitrary slip and prior history
    /// are not inferred.
    /// Requires no slip: Navier slip has a different unsteady impedance and
    /// cannot inherit these no-slip Basset or rotational-memory branches.
    /// A reservoir selected with `with_rotational_lift` refuses here: the
    /// quasi-steady lift law cannot be silently dropped or superposed on this
    /// history model, including when spin is introduced by a later drive.
    /// No spin, thermal shape/source change, free surface, wall, other particle,
    /// nonlinear convection, compressibility or Basset prehistory is inferred.
    ///
    /// Refine both spectral cutoffs, branch count and physical dt over the
    /// intended observation window. Omitted slow modes affect long times;
    /// omitted fast modes affect the startup layer. The exact steady Stokes
    /// resistance and instantaneous added mass rho*V/2 are never fitted.
    pub fn unsteady_motion<'s, 'b>(
        &'s self,
        body: &'b PrincipalBody<'b>,
        state: RigidBodyState,
        resolution: BassetResolution,
        cancelled: impl FnMut() -> bool,
    ) -> Result<UnsteadySphereMotion<'s, 'f, 'b, F>, Error> {
        self.prepare_unsteady_motion(body, state, resolution, None, false, cancelled)
    }

    /// Start with explicit translational preparation instead of discarding
    /// the history of an already moving particle.
    ///
    /// The current native momentum supplies the endpoint slip; [`SlipPrehistory`]
    /// declares how it was reached. Its exact branch responses and nonzero
    /// kinetic storage enter the SAME subsequent native midpoint, including
    /// when ambient acceleration or an applied force changes at time zero.
    /// No past work is charged to the new accounting window. Zero spin and
    /// all original material, geometry, Reynolds and resolution gates remain.
    pub fn unsteady_motion_with_prehistory<'s, 'b>(
        &'s self,
        body: &'b PrincipalBody<'b>,
        state: RigidBodyState,
        resolution: BassetResolution,
        prehistory: SlipPrehistory,
        cancelled: impl FnMut() -> bool,
    ) -> Result<UnsteadySphereMotion<'s, 'f, 'b, F>, Error> {
        self.prepare_unsteady_motion(
            body,
            state,
            resolution,
            Some(prehistory),
            false,
            cancelled,
        )
    }

    fn prepare_unsteady_motion<'s, 'b>(
        &'s self,
        body: &'b PrincipalBody<'b>,
        state: RigidBodyState,
        resolution: BassetResolution,
        prehistory: Option<SlipPrehistory>,
        allow_spin: bool,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<UnsteadySphereMotion<'s, 'f, 'b, F>, Error> {
        poll(&mut cancelled)?;
        self.require_no_slip()?;
        if self.rotational_lift_enabled() {
            return Err(Error::Invalid(
                "quasi-steady rotational lift cannot be combined with unsteady sphere history",
            ));
        }
        self.admit_body(body)?;
        self.load_frame(body.frame(), state, &mut cancelled)?;
        let lo = resolution.minimum_wavenumber_radius;
        let hi = resolution.maximum_wavenumber_radius;
        let comoving_momentum =
            vector(self.ambient_velocity_world.scale(body.properties().mass()))?;
        if prehistory.is_none() && state.linear_momentum_world() != comoving_momentum {
            return Err(Error::Invalid(
                "unsteady sphere needs canonical comoving momentum or explicit slip prehistory",
            ));
        }
        if !allow_spin && state.angular_momentum_body() != Vec3::ZERO {
            return Err(Error::Invalid(
                "unsteady translational sphere needs zero spin; use the rotating history owner",
            ));
        }
        if !lo.is_finite()
            || lo <= 0.0
            || !hi.is_finite()
            || hi <= lo
            || !(1..=4096).contains(&resolution.intervals)
        {
            return Err(Error::Invalid(
                "unsteady sphere needs zero spin and finite positive ordered spectral cutoffs with 1..=4096 intervals",
            ));
        }
        let slip = vector(
            state
                .linear_momentum_world()
                .sub(comoving_momentum)
                .scale(body.properties().mass().recip()),
        )?;
        let added_mass = positive(0.5 * self.fluid_mass_kg)?;
        let mut fluid_storage = 0.5 * added_mass * slip.norm_squared();
        let width = positive((hi.ln() - lo.ln()) / resolution.intervals as f64)?;
        let mut branches = Vec::with_capacity(resolution.intervals);
        let mut fluid_velocities = Vec::with_capacity(resolution.intervals);
        for i in 0..resolution.intervals {
            poll(&mut cancelled)?;
            let x = positive((lo.ln() + (i as f64 + 0.5) * width).exp())?;
            let resistance =
                positive((2.0 / core::f64::consts::PI) * self.resistance_kg_s * x * width)?;
            // Geometry/density determine storage independently of viscosity.
            let mass = positive(
                12.0 * self.transport.density_kg_m3() * self.radius_m.powi(3) * width / x,
            )?;
            let rate = positive(resistance / mass)?;
            let fraction = match prehistory {
                None => 0.0,
                Some(preparation) => preparation.memory_fraction(rate).ok_or(Error::Invalid(
                    "slip prehistory requires a finite positive representable preparation timescale",
                ))?,
            };
            let velocity = vector(slip.scale(fraction))?;
            fluid_storage += 0.5 * mass * velocity.norm_squared();
            fluid_velocities.push(velocity);
            branches.push(BassetBranch {
                resistance_kg_s: resistance,
                mass_kg: mass,
                rate_per_s: rate,
            });
        }
        // Nonzero prehistory brings initial added/branch storage into the
        // accounting window; its finite spectral reference must be representable.
        finite(fluid_storage)?;
        poll(&mut cancelled)?;
        Ok(UnsteadySphereMotion {
            sphere: self,
            body,
            resolution,
            branches: branches.into(),
            fluid_velocities,
            rotation: None,
            ambient_velocity_world: self.ambient_velocity_world,
            state,
            time_s: 0.0,
        })
    }
}

impl<F: FluidTransport + ?Sized> UnsteadySphereMotion<'_, '_, '_, F> {
    pub const fn state(&self) -> RigidBodyState {
        self.state
    }
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub const fn resolution(&self) -> BassetResolution {
        self.resolution
    }
    pub fn branches(&self) -> &[BassetBranch] {
        &self.branches
    }
    /// Spectral branch velocities RELATIVE to the uniform far field, expressed
    /// in world axes. These are not sampled physical fluid-field velocities.
    pub fn fluid_velocities_world_m_s(&self) -> &[Vec3] {
        &self.fluid_velocities
    }
    /// Accepted uniform far-field velocity, retained across step and restart.
    pub const fn ambient_velocity_world_m_s(&self) -> Vec3 {
        self.ambient_velocity_world
    }
    pub fn added_mass_kg(&self) -> f64 {
        0.5 * self.sphere.fluid_mass_kg
    }

    /// Actual discrete force/velocity impedance [real, imaginary], exp(+i*w*t).
    /// Compare against 6*pi*mu*a*(1+a*sqrt(i*w/nu)) + i*w*rho*V/2.
    pub fn impedance_kg_s(&self, omega_rad_s: f64) -> Result<[f64; 2], Error> {
        positive(omega_rad_s)?;
        let mut z = [
            self.sphere.resistance_kg_s,
            self.added_mass_kg() * omega_rad_s,
        ];
        for branch in self.branches.iter() {
            let scale = omega_rad_s.max(branch.rate_per_s);
            let w = omega_rad_s / scale;
            let rate = branch.rate_per_s / scale;
            let denominator = positive(w * w + rate * rate)?;
            z[0] += branch.resistance_kg_s * w * w / denominator;
            z[1] += branch.resistance_kg_s * w * rate / denominator;
        }
        Ok([finite(z[0])?, finite(z[1])?])
    }

    /// Bound on omitted continuum tails only [kg/s], excluding quadrature/time
    /// errors: 2*c/pi*(x_min + omega*a^2/(nu*x_max)). Not a whole-model bound.
    pub fn omitted_tail_bound_kg_s(&self, omega_rad_s: f64) -> Result<f64, Error> {
        positive(omega_rad_s)?;
        finite(
            (2.0 / core::f64::consts::PI)
                * self.sphere.resistance_kg_s
                * (self.resolution.minimum_wavenumber_radius
                    + omega_rad_s * self.sphere.diffusion_time_s
                        / self.resolution.maximum_wavenumber_radius),
        )
    }

    /// Advance all body/fluid histories atomically with the native midpoint.
    /// Work is O(branches), storage O(branches), independent of elapsed steps.
    /// Each source coefficient is already prepared; there is no per-step source
    /// lookup, dense matrix solve or growing convolution history. Native Reynolds,
    /// finite arithmetic, separate energy/impulse gates and cancellation precede
    /// publication. Midpoint passivity does not certify time or spectral accuracy.
    pub fn step(
        &self,
        controls: UnsteadyStepControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<(Self, UnsteadySphereStep), Error> {
        self.step_with_wrench(Wrench::ZERO, controls, cancelled)
    }

    /// Advance with a world-frame centre-of-mass force held over this step.
    ///
    /// The drive enters the same native midpoint solve as fluid added mass,
    /// steady drag and every existing memory branch. Calling this repeatedly
    /// with time-dependent forces retains all accepted viscous prehistory.
    /// Sample a smooth drive at the temporal midpoint and refine dt to resolve
    /// it; piecewise-constant drive discontinuities belong at step boundaries.
    /// Applied work uses the actual native midpoint velocity and is distinct
    /// from both stored fluid energy and dissipated heat. Momentum admission
    /// includes the applied impulse. The nonrotating translational Basset law
    /// refuses every nonzero applied torque; it does not infer spin drag.
    /// This holds the currently accepted ambient velocity, including after a
    /// previous [`Self::step_with_flow_acceleration`] call.
    pub fn step_with_wrench(
        &self,
        applied: Wrench,
        controls: UnsteadyStepControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<(Self, UnsteadySphereStep), Error> {
        self.advance_with_wrench(applied, Vec3::ZERO, controls, None, cancelled)
    }

    /// Advance in a spatially uniform flow with held world acceleration [m/s²].
    ///
    /// The accepted ambient velocity advances by `a_ambient * dt`; no unrelated
    /// endpoint can reset it or discard history. Sample a smooth acceleration
    /// and applied wrench at the midpoint and refine dt. Changing acceleration
    /// at a boundary is allowed; an ambient velocity jump/impulse is not.
    ///
    /// The original source mass gives the pressure-gradient force `m_fluid*a`.
    /// Added mass and Basset history respond to slip acceleration and slip,
    /// respectively. For a neutrally buoyant, comoving unforced sphere the
    /// density-contrast driving term vanishes, so it follows the ambient flow.
    /// Background source work is explicit, with no artificial dissipation.
    /// Spatial gradients, vorticity, Faxen corrections and convection remain
    /// outside this uniform-flow model. Refusal changes neither accepted body,
    /// ambient velocity nor memory.
    pub fn step_with_flow_acceleration(
        &self,
        applied: Wrench,
        ambient_acceleration_world_m_s2: Vec3,
        controls: UnsteadyStepControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<(Self, UnsteadySphereStep), Error> {
        self.advance_with_wrench(
            applied,
            ambient_acceleration_world_m_s2,
            controls,
            None,
            cancelled,
        )
    }

    fn advance_with_wrench(
        &self,
        applied: Wrench,
        ambient_acceleration: Vec3,
        controls: UnsteadyStepControls,
        rotational_controls: Option<rotation::RotatingStepControls>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(Self, UnsteadySphereStep), Error> {
        poll(&mut cancelled)?;
        vector(applied.force_world)?;
        vector(applied.torque_body)?;
        vector(ambient_acceleration)?;
        if self.rotation.is_none() && applied.torque_body != Vec3::ZERO {
            return Err(Error::Invalid(
                "unsteady translational sphere drag requires zero applied torque",
            ));
        }
        if self.rotation.is_some() != rotational_controls.is_some() {
            return Err(Error::Invalid(
                "rotational history requires its explicit controls",
            ));
        }
        let dt = positive(controls.duration_s)?;
        positive(controls.maximum_energy_defect_j)?;
        positive(controls.maximum_momentum_defect_ns)?;
        let time_s = finite(self.time_s + dt)?;
        if time_s <= self.time_s {
            return Err(Error::Invalid("unsteady clock did not advance"));
        }
        let p = self.body.properties();
        let v0 = self.state.center_of_mass_velocity_world(p)?;
        let ambient_change = vector(ambient_acceleration.scale(dt))?;
        let ambient_after = vector(self.ambient_velocity_world.add(ambient_change))?;
        for (acceleration, before, after) in [
            (
                ambient_acceleration.x,
                self.ambient_velocity_world.x,
                ambient_after.x,
            ),
            (
                ambient_acceleration.y,
                self.ambient_velocity_world.y,
                ambient_after.y,
            ),
            (
                ambient_acceleration.z,
                self.ambient_velocity_world.z,
                ambient_after.z,
            ),
        ] {
            if acceleration != 0.0 && before == after {
                return Err(Error::Invalid("unresolved ambient velocity increment"));
            }
        }
        // Subtract before dividing: canonical co-motion remains exactly zero
        // even when (m*U)*recip(m) does not round-trip to U. This introduces no
        // velocity tolerance, clipping or alternate body trajectory.
        let ambient_momentum_before = vector(self.ambient_velocity_world.scale(p.mass()))?;
        let slip0 = vector(
            self.state
                .linear_momentum_world()
                .sub(ambient_momentum_before)
                .scale(p.mass().recip()),
        )?;
        let half_dt = positive(0.5 * dt)?;
        let added_mass = self.added_mass_kg();
        let moving_mass = positive(p.mass() + added_mass)?;
        let effective_weight = vector(
            self.sphere
                .gravity
                .acceleration_world()
                .scale(p.mass() - self.sphere.fluid_mass_kg),
        )?;
        // Eliminate each branch's native midpoint relation. Retain the summed
        // intercept and slope, so evaluating the solved history force needs no
        // second quadrature traversal or additional workspace.
        let mut memory_resistance = 0.0;
        let mut history_at_start = Vec3::ZERO;
        for (branch, &u0) in self.branches.iter().zip(&self.fluid_velocities) {
            poll(&mut cancelled)?;
            let denominator = positive(1.0 + half_dt * branch.rate_per_s)?;
            let c = positive(branch.resistance_kg_s / denominator)?;
            memory_resistance = finite(memory_resistance + c)?;
            history_at_start = vector(history_at_start.sub(slip0.sub(u0).scale(c)))?;
        }
        let resistance = finite(self.sphere.resistance_kg_s + memory_resistance)?;
        let force = vector(
            self.sphere
                .gravity
                .acceleration_world()
                .sub(ambient_acceleration)
                .scale(p.mass() - self.sphere.fluid_mass_kg)
                .add(applied.force_world)
                .sub(slip0.scale(self.sphere.resistance_kg_s))
                .add(history_at_start),
        )?;
        // Solve the slip acceleration directly. Adding two separately large
        // pressure/added-mass forces first would lose the exact neutral limit.
        let acceleration =
            vector(force.scale(1.0 / positive(moving_mass + half_dt * resistance)?))?;
        let predicted_midpoint = vector(slip0.add(acceleration.scale(half_dt)))?;
        let history_force =
            vector(history_at_start.sub(acceleration.scale(half_dt * memory_resistance)))?;
        let added_force = vector(acceleration.scale(-added_mass))?;
        let buoyancy = vector(
            self.sphere
                .gravity
                .acceleration_world()
                .scale(-self.sphere.fluid_mass_kg),
        )?;
        let ambient_pressure = self.sphere.ambient_pressure_force(ambient_acceleration)?;
        let disturbance = vector(
            history_force
                .add(added_force)
                .sub(predicted_midpoint.scale(self.sphere.resistance_kg_s)),
        )?;
        let mut wrench = Wrench {
            force_world: vector(buoyancy.add(ambient_pressure).add(disturbance))?,
            torque_body: Vec3::ZERO,
        };
        let total_wrench = Wrench {
            force_world: vector(wrench.force_world.add(applied.force_world))?,
            torque_body: Vec3::ZERO,
        };
        let (dynamics, rotation_trial) = match (&self.rotation, rotational_controls) {
            (Some(history), Some(rotational)) => {
                let trial = history.native_candidate(
                    self,
                    total_wrench.force_world,
                    applied.torque_body,
                    rotational,
                    &mut cancelled,
                )?;
                wrench.torque_body = trial.fluid_torque_body_nm;
                (trial.dynamics, Some(trial))
            }
            _ => (
                RigidBodyIntegrator::new(self.sphere.gravity).step(
                    self.state,
                    p,
                    total_wrench,
                    dt,
                )?,
                None,
            ),
        };
        self.sphere.load_frame_in_flow(
            self.body.frame(),
            dynamics.state_after,
            ambient_after,
            &mut cancelled,
        )?;
        let v1 = dynamics.state_after.center_of_mass_velocity_world(p)?;
        let vm = vector(v0.add(v1).scale(0.5))?;
        let ambient_momentum_after = vector(ambient_after.scale(p.mass()))?;
        let slip1 = vector(
            dynamics
                .state_after
                .linear_momentum_world()
                .sub(ambient_momentum_after)
                .scale(p.mass().recip()),
        )?;
        let slip_midpoint = vector(slip0.add(slip1).scale(0.5))?;
        let ambient_midpoint = vector(self.ambient_velocity_world.add(ambient_after).scale(0.5))?;
        let mut fluid_velocities = Vec::with_capacity(self.branches.len());
        let mut memory_energy = 0.0;
        let mut heat = finite(dt * self.sphere.resistance_kg_s * slip_midpoint.norm_squared())?;
        let mut fluid_impulse = vector(slip1.sub(slip0).scale(added_mass))?;
        for (branch, &u0) in self.branches.iter().zip(&self.fluid_velocities) {
            poll(&mut cancelled)?;
            let alpha = finite(half_dt * branch.rate_per_s)?;
            let du = vector(slip_midpoint.sub(u0).scale(2.0 * (alpha / (1.0 + alpha))))?;
            let u1 = vector(u0.add(du))?;
            let um = vector(u0.add(u1).scale(0.5))?;
            memory_energy = finite(memory_energy + branch.mass_kg * du.dot(um))?;
            fluid_impulse = vector(fluid_impulse.add(du.scale(branch.mass_kg)))?;
            heat =
                finite(heat + dt * branch.resistance_kg_s * slip_midpoint.sub(um).norm_squared())?;
            fluid_velocities.push(u1);
        }
        let (rotation, rotation_receipt) =
            match (&self.rotation, rotation_trial, rotational_controls) {
                (Some(history), Some(trial), Some(rotational)) => {
                    let (next, receipt) = history.finish(
                        self,
                        trial,
                        applied.torque_body,
                        rotational,
                        &mut cancelled,
                    )?;
                    memory_energy = finite(memory_energy + receipt.memory_energy_change_j)?;
                    heat = finite(heat + receipt.viscous_heat_j)?;
                    (Some(next), Some(receipt))
                }
                _ => (None, None),
            };
        let added_energy = finite(added_mass * slip1.sub(slip0).dot(slip_midpoint))?;
        let dx = dynamics
            .state_after
            .pose()
            .position_world()
            .sub(self.state.pose().position_world());
        let mut stored = finite(
            dynamics
                .state_after
                .linear_momentum_world()
                .sub(self.state.linear_momentum_world())
                .dot(vm)
                - effective_weight.dot(dx)
                + added_energy
                + memory_energy,
        )?;
        let mut applied_power = finite(applied.force_world.dot(vm))?;
        if let Some(trial) = rotation_trial {
            stored = finite(
                stored
                    + (dynamics.diagnostics_after.rotational_kinetic_energy
                        - dynamics.diagnostics_before.rotational_kinetic_energy),
            )?;
            applied_power =
                finite(applied_power + applied.torque_body.dot(trial.omega_body_midpoint))?;
        }
        let applied_work = finite(dt * applied_power)?;
        if applied_power != 0.0 && applied_work == 0.0 {
            return Err(Error::Invalid(
                "unsteady drive work vanished in conversion to energy",
            ));
        }
        let background_power =
            finite(ambient_pressure.dot(vm) + disturbance.dot(ambient_midpoint))?;
        let background_work = finite(dt * background_power)?;
        if background_power != 0.0 && background_work == 0.0 {
            return Err(Error::Invalid(
                "unsteady background work vanished in conversion to energy",
            ));
        }
        let energy_defect = finite(stored + heat - applied_work - background_work)?;
        if energy_defect.abs() > controls.maximum_energy_defect_j {
            return Err(Error::EnergyDefect {
                defect_j: energy_defect,
                maximum_j: controls.maximum_energy_defect_j,
            });
        }
        let far_field = vector(slip_midpoint.scale(dt * self.sphere.resistance_kg_s))?;
        let momentum_residual = vector(
            dynamics
                .state_after
                .linear_momentum_world()
                .sub(self.state.linear_momentum_world())
                .add(fluid_impulse)
                .add(far_field)
                .sub(
                    vector(
                        effective_weight
                            .add(applied.force_world)
                            .add(ambient_pressure),
                    )?
                    .scale(dt),
                ),
        )?;
        if norm(momentum_residual)? > controls.maximum_momentum_defect_ns {
            return Err(Error::MomentumDefect {
                residual_ns: momentum_residual,
                maximum_ns: controls.maximum_momentum_defect_ns,
            });
        }
        poll(&mut cancelled)?;
        let next = Self {
            sphere: self.sphere,
            body: self.body,
            resolution: self.resolution,
            branches: Arc::clone(&self.branches),
            fluid_velocities,
            rotation,
            ambient_velocity_world: ambient_after,
            state: dynamics.state_after,
            time_s,
        };
        Ok((
            next,
            UnsteadySphereStep {
                dynamics,
                fluid_wrench: wrench,
                history_force_world_n: history_force,
                added_mass_force_world_n: added_force,
                ambient_pressure_force_world_n: ambient_pressure,
                viscous_heat_j: heat,
                applied_work_j: applied_work,
                background_work_j: background_work,
                added_mass_energy_change_j: added_energy,
                memory_energy_change_j: memory_energy,
                spectral_fluid_momentum_change_ns: fluid_impulse,
                far_field_impulse_ns: far_field,
                momentum_residual_ns: momentum_residual,
                stored_energy_change_j: stored,
                energy_defect_j: energy_defect,
                rotation: rotation_receipt,
            },
        ))
    }
}
