//! Unsteady Stokes settling with reciprocal fluid inertia and Basset memory.
//!
//! For a nonrotating no-slip sphere in a quiescent unbounded Newtonian fluid,
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

use super::*;
use std::sync::Arc;

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
    /// Includes buoyancy, steady resistance, memory and added-mass reaction.
    pub fluid_wrench: Wrench,
    pub history_force_world_n: Vec3,
    pub added_mass_force_world_n: Vec3,
    /// Steady Stokes loss plus all branch dashpot loss, never stored energy.
    pub viscous_heat_j: f64,
    pub added_mass_energy_change_j: f64,
    pub memory_energy_change_j: f64,
    /// Added mass and branch generalized momentum, not total field momentum.
    pub spectral_fluid_momentum_change_ns: Vec3,
    /// Impulse delivered to the eliminated stationary Stokes exterior.
    pub far_field_impulse_ns: Vec3,
    pub momentum_residual_ns: Vec3,
    /// Dry-body kinetic + gravity/buoyancy potential + represented fluid energy.
    pub stored_energy_change_j: f64,
    pub energy_defect_j: f64,
}

/// One restartable body and its complete accepted viscous-memory history.
/// The immutable source, original material mass and quadrature remain bound.
pub struct UnsteadySphereMotion<'s, 'f, 'b, F: FluidTransport + ?Sized> {
    sphere: &'s CreepingSphere<'f, F>,
    body: &'b PrincipalBody<'b>,
    resolution: BassetResolution,
    branches: Arc<[BassetBranch]>,
    fluid_velocities: Vec<Vec3>,
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
            state: self.state,
            time_s: self.time_s,
        }
    }
}

impl<'f, F: FluidTransport + ?Sized> CreepingSphere<'f, F> {
    /// Release the original sphere from rest, with quiescent fluid memory.
    ///
    /// This retains the original source snapshot, geometry, native gravity and
    /// Reynolds cap. The quasi-steady diffusion/relaxation restriction is replaced
    /// by actual added mass and viscous memory; it is NOT relaxed by density.
    /// Both initial momenta and the prescribed ambient velocity must be zero.
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
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<UnsteadySphereMotion<'s, 'f, 'b, F>, Error> {
        poll(&mut cancelled)?;
        self.admit_body(body)?;
        self.load_frame(body.frame(), state, &mut cancelled)?;
        let lo = resolution.minimum_wavenumber_radius;
        let hi = resolution.maximum_wavenumber_radius;
        if self.ambient_velocity_world != Vec3::ZERO
            || state.linear_momentum_world() != Vec3::ZERO
            || state.angular_momentum_body() != Vec3::ZERO
            || !lo.is_finite()
            || lo <= 0.0
            || !hi.is_finite()
            || hi <= lo
            || !(1..=4096).contains(&resolution.intervals)
        {
            return Err(Error::Invalid(
                "unsteady sphere needs rest, quiescent fluid and finite positive ordered spectral cutoffs with 1..=4096 intervals",
            ));
        }
        positive(0.5 * self.fluid_mass_kg)?;
        let width = positive((hi.ln() - lo.ln()) / resolution.intervals as f64)?;
        let mut branches = Vec::with_capacity(resolution.intervals);
        for i in 0..resolution.intervals {
            poll(&mut cancelled)?;
            let x = positive((lo.ln() + (i as f64 + 0.5) * width).exp())?;
            let resistance =
                positive((2.0 / core::f64::consts::PI) * self.resistance_kg_s * x * width)?;
            // Geometry/density determine storage independently of viscosity.
            let mass = positive(
                12.0 * self.transport.density_kg_m3() * self.radius_m.powi(3) * width / x,
            )?;
            branches.push(BassetBranch {
                resistance_kg_s: resistance,
                mass_kg: mass,
                rate_per_s: positive(resistance / mass)?,
            });
        }
        poll(&mut cancelled)?;
        Ok(UnsteadySphereMotion {
            sphere: self,
            body,
            resolution,
            branches: branches.into(),
            fluid_velocities: vec![Vec3::ZERO; resolution.intervals],
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
    pub fn fluid_velocities_world_m_s(&self) -> &[Vec3] {
        &self.fluid_velocities
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
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(Self, UnsteadySphereStep), Error> {
        poll(&mut cancelled)?;
        let dt = positive(controls.duration_s)?;
        positive(controls.maximum_energy_defect_j)?;
        positive(controls.maximum_momentum_defect_ns)?;
        let time_s = finite(self.time_s + dt)?;
        if time_s <= self.time_s {
            return Err(Error::Invalid("unsteady clock did not advance"));
        }
        let p = self.body.properties();
        let v0 = self.state.center_of_mass_velocity_world(p)?;
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
            history_at_start = vector(history_at_start.sub(v0.sub(u0).scale(c)))?;
        }
        let resistance = finite(self.sphere.resistance_kg_s + memory_resistance)?;
        let force = vector(
            effective_weight
                .sub(v0.scale(self.sphere.resistance_kg_s))
                .add(history_at_start),
        )?;
        let acceleration =
            vector(force.scale(1.0 / positive(moving_mass + half_dt * resistance)?))?;
        let predicted_midpoint = vector(v0.add(acceleration.scale(half_dt)))?;
        let history_force =
            vector(history_at_start.sub(acceleration.scale(half_dt * memory_resistance)))?;
        let added_force = vector(acceleration.scale(-added_mass))?;
        let buoyancy = vector(
            self.sphere
                .gravity
                .acceleration_world()
                .scale(-self.sphere.fluid_mass_kg),
        )?;
        let wrench = Wrench {
            force_world: vector(
                buoyancy
                    .add(history_force)
                    .add(added_force)
                    .sub(predicted_midpoint.scale(self.sphere.resistance_kg_s)),
            )?,
            torque_body: Vec3::ZERO,
        };
        let dynamics =
            RigidBodyIntegrator::new(self.sphere.gravity).step(self.state, p, wrench, dt)?;
        self.sphere
            .load_frame(self.body.frame(), dynamics.state_after, &mut cancelled)?;
        let v1 = dynamics.state_after.center_of_mass_velocity_world(p)?;
        let vm = vector(v0.add(v1).scale(0.5))?;
        let mut fluid_velocities = Vec::with_capacity(self.branches.len());
        let mut memory_energy = 0.0;
        let mut heat = finite(dt * self.sphere.resistance_kg_s * vm.norm_squared())?;
        let mut fluid_impulse = vector(v1.sub(v0).scale(added_mass))?;
        for (branch, &u0) in self.branches.iter().zip(&self.fluid_velocities) {
            poll(&mut cancelled)?;
            let alpha = finite(half_dt * branch.rate_per_s)?;
            let du = vector(vm.sub(u0).scale(2.0 * (alpha / (1.0 + alpha))))?;
            let u1 = vector(u0.add(du))?;
            let um = vector(u0.add(u1).scale(0.5))?;
            memory_energy = finite(memory_energy + branch.mass_kg * du.dot(um))?;
            fluid_impulse = vector(fluid_impulse.add(du.scale(branch.mass_kg)))?;
            heat = finite(heat + dt * branch.resistance_kg_s * vm.sub(um).norm_squared())?;
            fluid_velocities.push(u1);
        }
        let added_energy = finite(added_mass * v1.sub(v0).dot(vm))?;
        let dx = dynamics
            .state_after
            .pose()
            .position_world()
            .sub(self.state.pose().position_world());
        let stored = finite(
            dynamics.diagnostics_after.translational_kinetic_energy
                - dynamics.diagnostics_before.translational_kinetic_energy
                - effective_weight.dot(dx)
                + added_energy
                + memory_energy,
        )?;
        let energy_defect = finite(stored + heat)?;
        if energy_defect.abs() > controls.maximum_energy_defect_j {
            return Err(Error::EnergyDefect {
                defect_j: energy_defect,
                maximum_j: controls.maximum_energy_defect_j,
            });
        }
        let far_field = vector(vm.scale(dt * self.sphere.resistance_kg_s))?;
        let momentum_residual = vector(
            dynamics
                .state_after
                .linear_momentum_world()
                .sub(self.state.linear_momentum_world())
                .add(fluid_impulse)
                .add(far_field)
                .sub(effective_weight.scale(dt)),
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
                viscous_heat_j: heat,
                added_mass_energy_change_j: added_energy,
                memory_energy_change_j: memory_energy,
                spectral_fluid_momentum_change_ns: fluid_impulse,
                far_field_impulse_ns: far_field,
                momentum_residual_ns: momentum_residual,
                stored_energy_change_j: stored,
                energy_defect_j: energy_defect,
            },
        ))
    }
}
