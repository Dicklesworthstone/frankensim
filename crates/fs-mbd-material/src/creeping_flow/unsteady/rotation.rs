//! Causal spin resistance on the same native translating sphere.
//!
//! For a no-slip sphere in quiescent unbounded Newtonian fluid, with
//! tau=a^2/nu and d=8*pi*mu*a^3, the exact rotary impedance is
//! Z_r(s)=d*(1+tau*s/(3*(1+sqrt(tau*s)))). There is no instantaneous
//! added moment of inertia. See "The Resistance of an Arbitrary Body in
//! Confined Unsteady Stokes Flow", Fluids 9 (2024), 260, equations (10)-(11):
//! <https://www.mdpi.com/2311-5521/9/11/260>.
//!
//! Its positive diffusive measure, using x=a*sqrt(lambda/nu), is
//! d*c(x)=2*d*x^3/(3*pi*(1+x^2))*d(log(x)). Each branch stores
//! J_j*|omega_j|^2/2, with J_j=c_j/lambda_j. These spectral angular
//! velocities live in WORLD coordinates; rotating them with the body would
//! invent fluid precession. Translation and spin share one native Cayley step.
//! Its finite-step world-torque pullback closes angular impulse and power
//! together: Rbar^-1*T = (I-dt*[omega_mid]x/2)*R0^T*T.

use super::*;

/// One source- and geometry-derived rotary-memory mode.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RotationalMemoryBranch {
    pub resistance_kg_m2_s: f64,
    pub inertia_kg_m2: f64,
    pub rate_per_s: f64,
}

/// Numerical controls for translation, fluid spin history and native attitude.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RotatingStepControls {
    pub mechanical: UnsteadyStepControls,
    /// Absolute tolerance on world angular impulse [N m s], including native
    /// solve/roundoff error. The Cayley pullback conserves the discrete exchange.
    pub maximum_angular_momentum_defect_nms: f64,
    /// In (0,0.5]; bound dt*|omega_mid|, the native Cayley rotation coordinate.
    pub maximum_rotation_increment_rad: f64,
    /// In [64*epsilon,1e-8], relative closure of fluid/body torque.
    pub relative_torque_tolerance: f64,
    /// Bounded attitude/load coupling, in 1..=64.
    pub maximum_coupling_iterations: usize,
}

/// Rotary part of the accepted body/fluid work and impulse receipt.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RotationalHistoryStep {
    /// Physical world-frame drag moment. The native body torque is its Cayley
    /// pullback, recorded separately in the parent fluid wrench.
    pub fluid_torque_world_nm: Vec3,
    pub history_torque_world_nm: Vec3,
    pub memory_energy_change_j: f64,
    pub viscous_heat_j: f64,
    pub spectral_fluid_angular_momentum_change_nms: Vec3,
    pub far_field_angular_impulse_nms: Vec3,
    pub angular_momentum_residual_nms: Vec3,
    pub coupling_iterations: usize,
    pub relative_torque_residual: f64,
    pub rotation_increment_rad: f64,
}

#[derive(Clone)]
pub(super) struct RotationHistory {
    branches: Arc<[RotationalMemoryBranch]>,
    velocities: Vec<Vec3>,
}

/// Rest-started translation and rotation with separate fluid-memory channels.
/// The original material rigid body remains the sole owner of pose/momentum.
/// The continuum/nonlinear-convection limits of the parent sphere still apply.
pub struct RotatingSphereMotion<'s, 'f, 'b, F: FluidTransport + ?Sized> {
    inner: UnsteadySphereMotion<'s, 'f, 'b, F>,
}

impl<F: FluidTransport + ?Sized> Clone for RotatingSphereMotion<'_, '_, '_, F> {
    fn clone(&self) -> Self {
        Self {
            inner: self.inner.clone(),
        }
    }
}

impl<'f, F: FluidTransport + ?Sized> CreepingSphere<'f, F> {
    /// Admit a quiescent fluid and stationary sphere, including spin memory.
    ///
    /// Both momenta initially vanish. No prehistory, ambient vorticity, wall,
    /// finite-Reynolds lift/Magnus force, or nonspherical coupling is inferred.
    /// Spectral cutoffs/count and dt need refinement over the observation window.
    pub fn unsteady_rotating_motion<'s, 'b>(
        &'s self,
        body: &'b PrincipalBody<'b>,
        state: RigidBodyState,
        resolution: BassetResolution,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RotatingSphereMotion<'s, 'f, 'b, F>, Error> {
        let mut inner = self.unsteady_motion(body, state, resolution, &mut cancelled)?;
        inner.rotation = Some(RotationHistory::prepare(&inner, &mut cancelled)?);
        poll(&mut cancelled)?;
        Ok(RotatingSphereMotion { inner })
    }
}

impl<F: FluidTransport + ?Sized> RotatingSphereMotion<'_, '_, '_, F> {
    pub const fn state(&self) -> RigidBodyState {
        self.inner.state
    }
    pub const fn time_s(&self) -> f64 {
        self.inner.time_s
    }
    pub const fn resolution(&self) -> BassetResolution {
        self.inner.resolution
    }
    pub fn translational_branches(&self) -> &[BassetBranch] {
        &self.inner.branches
    }
    pub fn fluid_velocities_world_m_s(&self) -> &[Vec3] {
        &self.inner.fluid_velocities
    }
    pub fn rotational_branches(&self) -> &[RotationalMemoryBranch] {
        self.inner
            .rotation
            .as_ref()
            .map_or(&[], |r| r.branches.as_ref())
    }
    pub fn fluid_angular_velocities_world_rad_s(&self) -> &[Vec3] {
        self.inner
            .rotation
            .as_ref()
            .map_or(&[], |r| r.velocities.as_slice())
    }
    pub fn added_mass_kg(&self) -> f64 {
        self.inner.added_mass_kg()
    }
    pub fn translational_impedance_kg_s(&self, omega_rad_s: f64) -> Result<[f64; 2], Error> {
        self.inner.impedance_kg_s(omega_rad_s)
    }

    /// Actual discrete rotary impedance [real, imaginary], exp(+i*omega*t).
    /// Its exact limit is d*(1+i*omega*tau/(3*(1+sqrt(i*omega*tau)))).
    pub fn rotational_impedance_kg_m2_s(&self, omega_rad_s: f64) -> Result<[f64; 2], Error> {
        positive(omega_rad_s)?;
        let mut z = [self.inner.sphere.rotary_resistance_kg_m2_s, 0.0];
        for branch in self.rotational_branches() {
            let scale = omega_rad_s.max(branch.rate_per_s);
            let w = omega_rad_s / scale;
            let rate = branch.rate_per_s / scale;
            let denominator = positive(w * w + rate * rate)?;
            z[0] = finite(z[0] + branch.resistance_kg_m2_s * w * w / denominator)?;
            z[1] = finite(z[1] + branch.resistance_kg_m2_s * w * rate / denominator)?;
        }
        Ok(z)
    }

    /// Bound on omitted continuum tails only; excludes quadrature/time errors.
    /// With p=omega*tau, the bound is 2*d/(3*pi)*(x_min^3/3+p/x_max).
    pub fn rotational_omitted_tail_bound_kg_m2_s(&self, omega_rad_s: f64) -> Result<f64, Error> {
        positive(omega_rad_s)?;
        let r = self.inner.resolution;
        finite(
            (2.0 / (3.0 * core::f64::consts::PI))
                * self.inner.sphere.rotary_resistance_kg_m2_s
                * (r.minimum_wavenumber_radius.powi(3) / 3.0
                    + omega_rad_s * self.inner.sphere.diffusion_time_s
                        / r.maximum_wavenumber_radius),
        )
    }

    pub fn step(
        &self,
        controls: RotatingStepControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<(Self, UnsteadySphereStep), Error> {
        self.step_with_wrench(Wrench::ZERO, controls, cancelled)
    }

    /// Hold an external COM force in world coordinates and free couple in the
    /// native principal-body coordinates over the step. Sample smooth drives
    /// at the midpoint; discontinuities belong at step boundaries. The native
    /// torque closure uses world-fixed fluid history and admits its
    /// actual work and angular impulse before returning any candidate.
    pub fn step_with_wrench(
        &self,
        applied: Wrench,
        controls: RotatingStepControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<(Self, UnsteadySphereStep), Error> {
        let (inner, receipt) = self.inner.advance_with_wrench(
            applied,
            controls.mechanical,
            Some(controls),
            cancelled,
        )?;
        Ok((Self { inner }, receipt))
    }
}

#[derive(Clone, Copy)]
pub(super) struct RotationTrial {
    pub dynamics: StepReceipt,
    pub fluid_torque_body_nm: Vec3,
    pub omega_body_midpoint: Vec3,
    omega_world_midpoint: Vec3,
    fluid_torque_world_nm: Vec3,
    history_torque_world_nm: Vec3,
    iterations: usize,
    residual: f64,
    rotation_increment_rad: f64,
}

impl RotationHistory {
    fn prepare<F: FluidTransport + ?Sized>(
        motion: &UnsteadySphereMotion<'_, '_, '_, F>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        let r = motion.resolution;
        let width = positive(
            (r.maximum_wavenumber_radius.ln() - r.minimum_wavenumber_radius.ln())
                / r.intervals as f64,
        )?;
        // Dividing the two continuum spectral measures gives
        // (d/(3*c))*x^2/(1+x^2). Their diffusion rates remain identical.
        let factor = positive(
            motion.sphere.rotary_resistance_kg_m2_s / (3.0 * motion.sphere.resistance_kg_s),
        )?;
        let mut branches = Vec::with_capacity(r.intervals);
        for (i, translation) in motion.branches.iter().enumerate() {
            poll(cancelled)?;
            let x = positive((r.minimum_wavenumber_radius.ln() + (i as f64 + 0.5) * width).exp())?;
            let ratio = if x <= 1.0 {
                x * x / (1.0 + x * x)
            } else {
                1.0 / (1.0 + (1.0 / x) * (1.0 / x))
            };
            let resistance = positive(translation.resistance_kg_s * factor * ratio)?;
            let inertia = positive(translation.mass_kg * factor * ratio)?;
            branches.push(RotationalMemoryBranch {
                resistance_kg_m2_s: resistance,
                inertia_kg_m2: inertia,
                rate_per_s: positive(resistance / inertia)?,
            });
        }
        poll(cancelled)?;
        Ok(Self {
            branches: branches.into(),
            velocities: vec![Vec3::ZERO; r.intervals],
        })
    }

    pub(super) fn native_candidate<F: FluidTransport + ?Sized>(
        &self,
        motion: &UnsteadySphereMotion<'_, '_, '_, F>,
        total_force_world_n: Vec3,
        applied_torque_body_nm: Vec3,
        controls: RotatingStepControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<RotationTrial, Error> {
        positive(controls.maximum_angular_momentum_defect_nms)?;
        if !controls.maximum_rotation_increment_rad.is_finite()
            || controls.maximum_rotation_increment_rad <= 0.0
            || controls.maximum_rotation_increment_rad > 0.5
            || !controls.relative_torque_tolerance.is_finite()
            || !(64.0 * f64::EPSILON..=1e-8).contains(&controls.relative_torque_tolerance)
            || !(1..=64).contains(&controls.maximum_coupling_iterations)
        {
            return Err(Error::Invalid(
                "invalid rotational history coupling controls",
            ));
        }
        let dt = positive(controls.mechanical.duration_s)?;
        let half_dt = positive(0.5 * dt)?;
        let d = motion.sphere.rotary_resistance_kg_m2_s;
        let mut slope = d;
        let mut intercept = Vec3::ZERO;
        for (branch, &u0) in self.branches.iter().zip(&self.velocities) {
            poll(cancelled)?;
            let c =
                positive(branch.resistance_kg_m2_s / positive(1.0 + half_dt * branch.rate_per_s)?)?;
            slope = finite(slope + c)?;
            intercept = vector(intercept.add(u0.scale(c)))?;
        }
        let p = motion.body.properties();
        let inertia = p.principal_inertia_body();
        let l0 = motion.state.angular_momentum_body();
        let q0 = motion.state.pose().orientation();
        let history_body = vector(q0.rotate_world_to_body(intercept))?;
        let drive = vector(history_body.add(applied_torque_body_nm))?;
        let diagonal = vector(Vec3::new(
            positive(inertia.x + half_dt * slope)?,
            positive(inertia.y + half_dt * slope)?,
            positive(inertia.z + half_dt * slope)?,
        ))?;
        let skew = vector(history_body.scale(half_dt * half_dt))?;
        let mut omega_guess = p.angular_velocity_body_checked(l0)?;
        let integrator = RigidBodyIntegrator::new(motion.sphere.gravity);
        let mut residual = f64::INFINITY;
        for iteration in 1..=controls.maximum_coupling_iterations {
            poll(cancelled)?;
            // Eliminate the viscous stiffness and finite-rotation pullback
            // together. Gyro correction retains the original admitted inertia,
            // including its roundoff-sized departure from spherical symmetry.
            let momentum_guess = vector(Vec3::new(
                inertia.x * omega_guess.x,
                inertia.y * omega_guess.y,
                inertia.z * omega_guess.z,
            ))?;
            let right =
                vector(l0.add(drive.sub(omega_guess.cross(momentum_guess)).scale(half_dt)))?;
            let predicted = diagonal_minus_cross(diagonal, skew, right)?;
            let used = vector(
                history_body
                    .sub(predicted.scale(slope))
                    .sub(predicted.cross(history_body).scale(half_dt)),
            )?;
            let dynamics = integrator.step_cayley(
                motion.state,
                p,
                Wrench {
                    force_world: total_force_world_n,
                    torque_body: vector(used.add(applied_torque_body_nm))?,
                },
                dt,
            )?;
            let omega = p.angular_velocity_body_checked(
                l0.scale(0.5)
                    .add(dynamics.state_after.angular_momentum_body().scale(0.5)),
            )?;
            let increment = finite(dt * norm(omega)?)?;
            if increment > controls.maximum_rotation_increment_rad {
                return Err(Error::Invalid(
                    "unsteady sphere exceeds its declared rotation increment; reduce dt",
                ));
            }
            // R0*omega_mid is the spatial velocity conjugate to the native
            // Cayley coordinate. Rbar=(R0+R1)/2 maps the native torque impulse.
            let world = vector(q0.rotate_body_to_world(omega))?;
            let desired_world = vector(intercept.sub(world.scale(slope)))?;
            let desired = vector(
                history_body
                    .sub(omega.scale(slope))
                    .sub(omega.cross(history_body).scale(half_dt)),
            )?;
            let scale = norm(desired)?
                .max(norm(used)?)
                .max(norm(applied_torque_body_nm)?);
            residual = if scale == 0.0 {
                norm(desired.sub(used))?
            } else {
                finite(norm(desired.sub(used))? / scale)?
            };
            if residual <= controls.relative_torque_tolerance {
                poll(cancelled)?;
                return Ok(RotationTrial {
                    dynamics,
                    fluid_torque_body_nm: used,
                    omega_body_midpoint: omega,
                    omega_world_midpoint: world,
                    fluid_torque_world_nm: desired_world,
                    history_torque_world_nm: vector(desired_world.add(world.scale(d)))?,
                    iterations: iteration,
                    residual,
                    rotation_increment_rad: increment,
                });
            }
            omega_guess = omega;
        }
        Err(Error::RotationalCoupling {
            relative_residual: residual,
            tolerance: controls.relative_torque_tolerance,
            iterations: controls.maximum_coupling_iterations,
        })
    }

    pub(super) fn finish<F: FluidTransport + ?Sized>(
        &self,
        motion: &UnsteadySphereMotion<'_, '_, '_, F>,
        trial: RotationTrial,
        applied_torque_body_nm: Vec3,
        controls: RotatingStepControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Self, RotationalHistoryStep), Error> {
        let dt = controls.mechanical.duration_s;
        let half_dt = 0.5 * dt;
        let omega = trial.omega_world_midpoint;
        let mut velocities = Vec::with_capacity(self.velocities.len());
        let mut energy = 0.0;
        let mut impulse = Vec3::ZERO;
        let mut heat = finite(dt * motion.sphere.rotary_resistance_kg_m2_s * omega.norm_squared())?;
        for (branch, &u0) in self.branches.iter().zip(&self.velocities) {
            poll(cancelled)?;
            let alpha = finite(half_dt * branch.rate_per_s)?;
            let du = vector(omega.sub(u0).scale(2.0 * (alpha / (1.0 + alpha))))?;
            let u1 = vector(u0.add(du))?;
            let um = vector(u0.scale(0.5).add(u1.scale(0.5)))?;
            energy = finite(energy + branch.inertia_kg_m2 * du.dot(um))?;
            impulse = vector(impulse.add(du.scale(branch.inertia_kg_m2)))?;
            heat = finite(heat + dt * branch.resistance_kg_m2_s * omega.sub(um).norm_squared())?;
            velocities.push(u1);
        }
        let before = motion.state;
        let after = trial.dynamics.state_after;
        let l0 = vector(
            before
                .pose()
                .orientation()
                .rotate_body_to_world(before.angular_momentum_body()),
        )?;
        let l1 = vector(
            after
                .pose()
                .orientation()
                .rotate_body_to_world(after.angular_momentum_body()),
        )?;
        let far_field = vector(omega.scale(dt * motion.sphere.rotary_resistance_kg_m2_s))?;
        let applied = vector(
            before
                .pose()
                .orientation()
                .rotate_body_to_world(applied_torque_body_nm)
                .scale(0.5)
                .add(
                    after
                        .pose()
                        .orientation()
                        .rotate_body_to_world(applied_torque_body_nm)
                        .scale(0.5),
                )
                .scale(dt),
        )?;
        let residual = vector(l1.sub(l0).add(impulse).add(far_field).sub(applied))?;
        if norm(residual)? > controls.maximum_angular_momentum_defect_nms {
            return Err(Error::AngularMomentumDefect {
                residual_nms: residual,
                maximum_nms: controls.maximum_angular_momentum_defect_nms,
            });
        }
        poll(cancelled)?;
        Ok((
            Self {
                branches: Arc::clone(&self.branches),
                velocities,
            },
            RotationalHistoryStep {
                fluid_torque_world_nm: trial.fluid_torque_world_nm,
                history_torque_world_nm: trial.history_torque_world_nm,
                memory_energy_change_j: energy,
                viscous_heat_j: heat,
                spectral_fluid_angular_momentum_change_nms: impulse,
                far_field_angular_impulse_nms: far_field,
                angular_momentum_residual_nms: residual,
                coupling_iterations: trial.iterations,
                relative_torque_residual: trial.residual,
                rotation_increment_rad: trial.rotation_increment_rad,
            },
        ))
    }
}

// (diag(a)-[k]x)^-1 b. Symmetric diagonal scaling reduces the solve to
// (I-[z]x)^-1=(I+[z]x+z*z^T)/(1+|z|^2), whose denominator is positive.
fn diagonal_minus_cross(a: Vec3, k: Vec3, b: Vec3) -> Result<Vec3, Error> {
    let roots = Vec3::new(
        fs_math::det::sqrt(positive(a.x)?),
        fs_math::det::sqrt(positive(a.y)?),
        fs_math::det::sqrt(positive(a.z)?),
    );
    let z = vector(Vec3::new(
        k.x / roots.y / roots.z,
        k.y / roots.x / roots.z,
        k.z / roots.x / roots.y,
    ))?;
    let right = vector(Vec3::new(b.x / roots.x, b.y / roots.y, b.z / roots.z))?;
    let inverse = 1.0 / positive(1.0 + z.norm_squared())?;
    let scaled = vector(
        right
            .add(z.cross(right))
            .add(z.scale(z.dot(right)))
            .scale(inverse),
    )?;
    vector(Vec3::new(
        scaled.x / roots.x,
        scaled.y / roots.y,
        scaled.z / roots.z,
    ))
}
