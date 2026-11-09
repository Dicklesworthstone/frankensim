//! Simultaneous affine force and drag; no separated kick or work trapezoid.

use super::*;

mod flow;

impl IsotropicLinearDrag {
    /// Conservative relative-speed bound over one prescribed affine-force path.
    ///
    /// For held bath velocity, let `S0=|v0-U|`, `lambda=c/m`, and let `A` be
    /// the larger norm of the two supplied endpoint accelerations, including
    /// gravity. The positive exponential convolution gives
    /// `|v(t)-U| <= exp(-lambda*t)*S0 + integral_0^t exp(-lambda*u)du*A`.
    /// This scalar comparison is monotone, so its maximum is at zero or at the
    /// final time even when the actual vector velocity has an interior maximum.
    /// It has the continuous zero-drag limit `S0+duration*A`.
    ///
    /// The bound may be conservative for opposing or turning force histories;
    /// subdividing at their actual source endpoints can tighten admission. This
    /// is floating-point model admission, not an interval certificate. Inputs
    /// and fixed-spin assumptions are identical to [`Self::step_affine_force`].
    pub fn affine_force_relative_speed_bound_m_s(
        self,
        state: RigidBodyState,
        wrench_before: Wrench,
        force_after_world: Vec3,
        duration_s: f64,
    ) -> Result<f64, LinearDragError> {
        validate_affine_inputs(state, wrench_before, force_after_world, duration_s)?;
        let mass = dd(self.properties.mass());
        let momentum = components(state.linear_momentum_world()).map(dd);
        let ambient = components(self.ambient_velocity_world).map(dd);
        let gravity = components(self.gravity.acceleration_world()).map(dd);
        let first = components(wrench_before.force_world).map(dd);
        let last = components(force_after_world).map(dd);
        relative_speed_bound(
            core::array::from_fn(|i| (momentum[i] - mass * ambient[i]) / mass),
            core::array::from_fn(|i| (first[i] + mass * gravity[i]) / mass),
            core::array::from_fn(|i| (last[i] + mass * gravity[i]) / mass),
            Moments::new(
                self.translation_resistance_kg_s,
                self.properties.mass(),
                duration_s,
            )?,
        )
    }

    /// Advance a world force prescribed by its two endpoints over this interval.
    ///
    /// The force is `F(s)=(1-s)*wrench_before.force_world+s*force_after_world`,
    /// `s=t/duration_s`. These are authored source endpoints, not samples of an
    /// unknown velocity/configuration-dependent law. Gravity, bath velocity and
    /// resistance remain held. The initial wrench's body torque must be zero;
    /// the same native exact free-spin flow owns attitude and rotary loss.
    ///
    /// Translation, displacement, impulse, applied work and slip heat integrate
    /// the simultaneous affine/exponential path. In particular, mean applied
    /// force times displacement is not its work, and an impulse followed by a
    /// held-drag step is not this solution. Independent positive-square heat
    /// retains reversals and small relaxation numbers without a terminal-speed
    /// subtraction. Zero translation resistance uses the exact polynomial limit.
    ///
    /// Adjacent affine intervals may share their actual force endpoint. This
    /// does not choose source-event times or admit an accelerating fluid model:
    /// fluid pressure-gradient, added-mass and memory forces remain the physical
    /// caller's responsibility. A material caller must admit the whole velocity
    /// path, since an affine drive can produce an interior speed maximum.
    /// Storage and numerical residuals use the represented native endpoint, as
    /// in [`Self::step`]. Fixed stack storage, with no nonlinear solve or allocation.
    pub fn step_affine_force(
        self,
        state: RigidBodyState,
        wrench_before: Wrench,
        force_after_world: Vec3,
        duration_s: f64,
    ) -> Result<LinearDragStep, LinearDragError> {
        validate_affine_inputs(state, wrench_before, force_after_world, duration_s)?;
        if wrench_before.force_world == force_after_world {
            return self.step(state, wrench_before, duration_s);
        }
        let mass = dd(self.properties.mass());
        let dt = dd(duration_s);
        let resistance = dd(self.translation_resistance_kg_s);
        let exposure = relaxation_exposure(
            self.translation_resistance_kg_s,
            self.properties.mass(),
            duration_s,
        )?;
        let before_p = components(state.linear_momentum_world()).map(dd);
        let first = components(wrench_before.force_world).map(dd);
        let last = components(force_after_world).map(dd);
        let weight = components(self.gravity.acceleration_world()).map(|g| mass * dd(g));
        let total_first: [Dd; 3] = core::array::from_fn(|i| first[i] + weight[i]);
        let total_last: [Dd; 3] = core::array::from_fn(|i| last[i] + weight[i]);
        let ambient = components(self.ambient_velocity_world).map(dd);
        let mut path = TranslationPath {
            after_p: [Dd::ZERO; 3],
            displacement: [Dd::ZERO; 3],
            drag_impulse: [Dd::ZERO; 3],
            translation_heat: Dd::ZERO,
            applied_work: Dd::ZERO,
            total_impulse: core::array::from_fn(|i| {
                (total_first[i] * dd(0.5) + total_last[i] * dd(0.5)) * dt
            }),
        };
        let small = (exposure > 0.0 && exposure < 2.0).then(|| CurvatureMoments::new(exposure));
        let stiff = (exposure >= 2.0).then(|| ExponentialMoments::new(exposure));
        for i in 0..3 {
            let initial_velocity = before_p[i] / mass;
            let (mean_velocity, weighted_velocity) = if exposure == 0.0 {
                // No bath arithmetic: zero resistance must remain independent
                // of even a very large authored bath velocity.
                let first_increment = total_first[i] / mass * dt;
                let last_increment = total_last[i] / mass * dt;
                path.after_p[i] = before_p[i] + path.total_impulse[i];
                (
                    initial_velocity + first_increment / dd(3.0) + last_increment / dd(6.0),
                    initial_velocity * dd(0.5)
                        + first_increment * dd(5.0) / dd(24.0)
                        + last_increment / dd(8.0),
                )
            } else if let Some(moments) = &small {
                let x = dd(exposure);
                let relative = (before_p[i] - mass * ambient[i]) / mass;
                let b0 = total_first[i] / mass * dt;
                let b1 = (last[i] - first[i]) / mass * dt;
                let initial_slope = b0 - x * relative;
                let curvature = b1 - x * initial_slope;
                let mean_increment = initial_slope * dd(0.5) + curvature * moments.mean;
                let mean_relative = relative + mean_increment;
                path.after_p[i] = before_p[i]
                    + mass * (initial_slope + curvature * moments.endpoint);
                path.drag_impulse[i] = -resistance * dt * mean_relative;
                path.translation_heat = path.translation_heat
                    + positive_mean_square(
                        resistance * dt,
                        mean_relative,
                        initial_slope + dd(12.0) * curvature * moments.covariance,
                        curvature,
                        moments.residual_variance,
                    )?;
                (
                    initial_velocity + mean_increment,
                    initial_velocity * dd(0.5)
                        + initial_slope / dd(3.0)
                        + curvature * moments.weighted,
                )
            } else {
                let moments = stiff.as_ref().expect("positive exposure selects one basis");
                // r(s)=B+C*s+A*exp(-x*s). Use force/resistance directly so
                // a long interval does not first form acceleration*duration.
                let terminal_first = total_first[i] / resistance;
                let terminal_last = total_last[i] / resistance;
                let slope = terminal_last - terminal_first;
                let offset = terminal_first - slope / dd(exposure);
                let amplitude = (before_p[i] - mass * ambient[i]) / mass - offset;
                let transient_mean = amplitude * moments.mean;
                let transient_weighted =
                    (transient_mean - attenuated(amplitude, exposure)) / dd(exposure);
                let mean_relative = offset + slope * dd(0.5) + transient_mean;
                path.after_p[i] = mass * (ambient[i] + offset + slope)
                    + attenuated(before_p[i] - mass * (ambient[i] + offset), exposure);
                path.drag_impulse[i] = -resistance * dt * mean_relative;
                path.translation_heat = path.translation_heat
                    + positive_mean_square(
                        resistance * dt,
                        mean_relative,
                        slope + dd(12.0) * amplitude * moments.covariance,
                        amplitude,
                        moments.residual_variance,
                    )?;
                (
                    ambient[i] + mean_relative,
                    (ambient[i] + offset) * dd(0.5) + slope / dd(3.0) + transient_weighted,
                )
            };
            path.displacement[i] = dt * mean_velocity;
            // The endpoint-force basis integrates F(s)*v(s) independently.
            // The weighted moment is integral_0^1 s*v(s) ds, not v_end/2.
            path.applied_work = path.applied_work
                + dt * (first[i] * (mean_velocity - weighted_velocity)
                    + last[i] * weighted_velocity);
        }
        self.finish(state, duration_s, path)
    }
}

fn relative_speed_bound(
    initial: [Dd; 3],
    acceleration_before: [Dd; 3],
    acceleration_after: [Dd; 3],
    moments: Moments,
) -> Result<f64, LinearDragError> {
    let norm = |value: [Dd; 3]| -> Result<f64, LinearDragError> {
        let value = vector(value)?;
        scalar(dd(value.x.hypot(value.y).hypot(value.z)))
    };
    let initial = norm(initial)?;
    let acceleration = norm(acceleration_before)?.max(norm(acceleration_after)?);
    let endpoint = scalar(
        attenuated(dd(initial), moments.exposure) + dd(acceleration) * moments.integral_s,
    )?;
    Ok(initial.max(endpoint))
}

fn validate_affine_inputs(
    state: RigidBodyState,
    wrench_before: Wrench,
    force_after_world: Vec3,
    duration_s: f64,
) -> Result<(), LinearDragError> {
    state.validate()?;
    wrench_before.validate()?;
    if !force_after_world.is_finite() {
        return Err(DynamicsError::NonFinite("linear_drag.force_after_world").into());
    }
    if wrench_before.torque_body != Vec3::ZERO {
        return Err(LinearDragError::Invalid(
            "exact free-spin flow does not admit an applied torque",
        ));
    }
    if !duration_s.is_finite() || duration_s <= 0.0 {
        return Err(DynamicsError::InvalidStepDuration.into());
    }
    Ok(())
}

// Orthogonal projection onto {1,s-1/2} leaves a nonnegative residual norm.
// These three independent squares integrate the trajectory's physical heat;
// no stored-energy or source-work remainder enters this calculation.
fn positive_mean_square(
    weight: Dd,
    mean: Dd,
    linear: Dd,
    amplitude: Dd,
    residual_variance: Dd,
) -> Result<Dd, LinearDragError> {
    if scalar(residual_variance)? < 0.0 {
        return Err(LinearDragError::Invalid("negative affine-path variance"));
    }
    // Multiply a small weight before squaring a potentially large velocity.
    Ok((weight * mean) * mean
        + ((weight / dd(12.0)) * linear) * linear
        + ((weight * residual_variance) * amplitude) * amplitude)
}

// For x<2, r(s)=r0+d0*s+k*u(s), u=s²*phi2(x*s). All these moments are
// regular at zero. The variance is ||u-proj_{1,s}u||², with limit 1/720.
struct CurvatureMoments {
    endpoint: Dd,
    mean: Dd,
    weighted: Dd,
    covariance: Dd,
    residual_variance: Dd,
}
impl CurvatureMoments {
    fn new(exposure: f64) -> Self {
        let x = dd(exposure);
        let mut mean = Dd::ZERO;
        let mut covariance = Dd::ZERO;
        let mut residual_variance = Dd::ZERO;
        let mut moment_term = dd(1.0) / dd(24.0);
        let mut variance_term = dd(1.0) / dd(40_320.0);
        let mut power_two = 16.0;
        // At x<=2 the first omitted variance term is below 1.5e-24.
        // Coefficients below are exact integers through this bounded order.
        for order in 0..=32 {
            let n = f64::from(order);
            mean = mean + moment_term * dd(n + 4.0);
            covariance = covariance + moment_term * dd(0.5 * (n + 2.0));
            let numerator = power_two * (n * n * n + 5.0 * n * n + 2.0 * n + 16.0)
                - 4.0 * (n + 5.0) * (n + 10.0);
            residual_variance = residual_variance + variance_term * dd(numerator);
            moment_term = moment_term * -x / dd(n + 5.0);
            variance_term = variance_term * -x / dd(n + 9.0);
            power_two *= 2.0;
        }
        Self {
            endpoint: dd(0.5) - x * mean,
            mean,
            weighted: mean * dd(0.5) + covariance,
            covariance,
            residual_variance,
        }
    }
}

// For x>=2, the decaying exponential is sufficiently independent of {1,s}
// to evaluate its residual variance directly with retained arithmetic.
struct ExponentialMoments {
    mean: Dd,
    covariance: Dd,
    residual_variance: Dd,
}
impl ExponentialMoments {
    fn new(exposure: f64) -> Self {
        let x = dd(exposure);
        let lost = dd(decay::decayed_fraction(exposure));
        let remaining = dd(decay::remaining_fraction(exposure));
        let mean = lost / x;
        let weighted = (mean - remaining) / x;
        let covariance = weighted - mean * dd(0.5);
        let mean_square = mean * (dd(1.0) - lost * dd(0.5));
        Self {
            mean,
            covariance,
            residual_variance: mean_square - mean * mean - dd(12.0) * covariance * covariance,
        }
    }
}
