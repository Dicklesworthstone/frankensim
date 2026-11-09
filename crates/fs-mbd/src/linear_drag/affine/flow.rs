//! Affine bath translation uses the same exact slip moments as affine forcing.

use super::*;

impl IsotropicLinearDrag {
    /// Advance simultaneous affine world force and uniform bath velocity.
    ///
    /// The prepared bath velocity is U0 and `ambient_velocity_after_world` is
    /// U1; both histories interpolate their authored endpoints over this step.
    /// This integrates the moving-damper law `m v'=F+m g-c(v-U(t))`. It does
    /// not infer an accelerating fluid's pressure, added mass or memory forces.
    /// Those belong to the physical caller and its admitted constitutive owner.
    ///
    /// Slip heat, applied-force work and signed background work are independent
    /// path integrals. In particular, background work pairs drag with U(t), not
    /// with either endpoint or an endpoint-mean bath velocity. Free-spin decay,
    /// represented-state storage and numerical residuals retain their original
    /// owners. There is no allocation, quadrature or nonlinear solve.
    ///
    /// For continuation prepare the next bath at the accepted U1. Holding U
    /// reproduces `step_affine_force` exactly; zero drag ignores bath arithmetic.
    pub fn step_affine_force_and_flow(
        self,
        state: RigidBodyState,
        wrench_before: Wrench,
        force_after_world: Vec3,
        ambient_velocity_after_world: Vec3,
        duration_s: f64,
    ) -> Result<LinearDragStep, LinearDragError> {
        validate_affine_inputs(state, wrench_before, force_after_world, duration_s)?;
        validate_flow_endpoint(ambient_velocity_after_world)?;
        if ambient_velocity_after_world == self.ambient_velocity_world
            || self.translation_resistance_kg_s == 0.0
        {
            return self.step_affine_force(state, wrench_before, force_after_world, duration_s);
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
        let ambient_before = components(self.ambient_velocity_world).map(dd);
        let ambient_after = components(ambient_velocity_after_world).map(dd);
        let mut path = TranslationPath {
            after_p: [Dd::ZERO; 3],
            displacement: [Dd::ZERO; 3],
            drag_impulse: [Dd::ZERO; 3],
            translation_heat: Dd::ZERO,
            applied_work: Dd::ZERO,
            total_impulse: core::array::from_fn(|i| {
                (first[i] * dd(0.5) + last[i] * dd(0.5) + weight[i]) * dt
            }),
        };
        let small = (exposure < 2.0).then(|| CurvatureMoments::new(exposure));
        let stiff = (exposure >= 2.0).then(|| ExponentialMoments::new(exposure));
        let mut background_work = Dd::ZERO;
        for i in 0..3 {
            let initial_velocity = before_p[i] / mass;
            let relative = (before_p[i] - mass * ambient_before[i]) / mass;
            let bath_change = ambient_after[i] - ambient_before[i];
            let (mean_velocity, weighted_velocity, mean_slip, weighted_slip) =
                if let Some(moments) = &small {
                    let x = dd(exposure);
                    // The absolute initial slope avoids adding DeltaU back
                    // after a potentially severe relative-frame cancellation.
                    let initial_slope = total_first[i] / mass * dt - x * relative;
                    let relative_slope = initial_slope - bath_change;
                    let curvature = (last[i] - first[i]) / mass * dt - x * relative_slope;
                    let mean_increment = initial_slope * dd(0.5) + curvature * moments.mean;
                    let mean_slip = relative + relative_slope * dd(0.5) + curvature * moments.mean;
                    let weighted_slip = relative * dd(0.5)
                        + relative_slope / dd(3.0)
                        + curvature * moments.weighted;
                    path.after_p[i] =
                        before_p[i] + mass * (initial_slope + curvature * moments.endpoint);
                    path.translation_heat = path.translation_heat
                        + positive_mean_square(
                            resistance * dt,
                            mean_slip,
                            relative_slope + dd(12.0) * curvature * moments.covariance,
                            curvature,
                            moments.residual_variance,
                        )?;
                    (
                        initial_velocity + mean_increment,
                        initial_velocity * dd(0.5)
                            + initial_slope / dd(3.0)
                            + curvature * moments.weighted,
                        mean_slip,
                        weighted_slip,
                    )
                } else {
                    let moments = stiff.as_ref().expect("exposure selects one exact basis");
                    // Relative forcing includes -m*dU/dt. Balance it before
                    // dividing by resistance; exact co-acceleration stays quiet.
                    let terminal_first = (total_first[i] - mass * (bath_change / dt)) / resistance;
                    let slope = (last[i] - first[i]) / resistance;
                    let offset = terminal_first - slope / dd(exposure);
                    let amplitude = relative - offset;
                    let transient_mean = amplitude * moments.mean;
                    let transient_weighted =
                        (transient_mean - attenuated(amplitude, exposure)) / dd(exposure);
                    let mean_slip = offset + slope * dd(0.5) + transient_mean;
                    let weighted_slip = offset * dd(0.5) + slope / dd(3.0) + transient_weighted;
                    path.after_p[i] = mass * (ambient_after[i] + offset + slope)
                        + attenuated(before_p[i] - mass * (ambient_before[i] + offset), exposure);
                    path.translation_heat = path.translation_heat
                        + positive_mean_square(
                            resistance * dt,
                            mean_slip,
                            slope + dd(12.0) * amplitude * moments.covariance,
                            amplitude,
                            moments.residual_variance,
                        )?;
                    (
                        ambient_before[i] * dd(0.5) + ambient_after[i] * dd(0.5) + mean_slip,
                        ambient_before[i] / dd(6.0) + ambient_after[i] / dd(3.0) + weighted_slip,
                        mean_slip,
                        weighted_slip,
                    )
                };
            path.displacement[i] = dt * mean_velocity;
            path.drag_impulse[i] = -resistance * dt * mean_slip;
            path.applied_work = path.applied_work
                + dt * (first[i] * (mean_velocity - weighted_velocity)
                    + last[i] * weighted_velocity);
            background_work = background_work
                - resistance
                    * dt
                    * (ambient_before[i] * (mean_slip - weighted_slip)
                        + ambient_after[i] * weighted_slip);
        }
        self.finish_with_background_work(state, duration_s, path, background_work)
    }

    /// Conservative whole-path speed bound relative to the affine bath.
    ///
    /// Uses the original positive-convolution bound with relative accelerations
    /// `F0/m+g-dU/dt` and `F1/m+g-dU/dt`. Interior slip maxima are covered;
    /// physical fluid validity and source bandwidth remain the caller's job.
    pub fn affine_force_and_flow_relative_speed_bound_m_s(
        self,
        state: RigidBodyState,
        wrench_before: Wrench,
        force_after_world: Vec3,
        ambient_velocity_after_world: Vec3,
        duration_s: f64,
    ) -> Result<f64, LinearDragError> {
        validate_affine_inputs(state, wrench_before, force_after_world, duration_s)?;
        validate_flow_endpoint(ambient_velocity_after_world)?;
        if ambient_velocity_after_world == self.ambient_velocity_world {
            return self.affine_force_relative_speed_bound_m_s(
                state,
                wrench_before,
                force_after_world,
                duration_s,
            );
        }
        let mass = dd(self.properties.mass());
        let dt = dd(duration_s);
        let momentum = components(state.linear_momentum_world()).map(dd);
        let ambient = components(self.ambient_velocity_world).map(dd);
        let after = components(ambient_velocity_after_world).map(dd);
        let gravity = components(self.gravity.acceleration_world()).map(dd);
        let first = components(wrench_before.force_world).map(dd);
        let last = components(force_after_world).map(dd);
        relative_speed_bound(
            core::array::from_fn(|i| (momentum[i] - mass * ambient[i]) / mass),
            core::array::from_fn(|i| {
                (first[i] + mass * gravity[i]) / mass - (after[i] - ambient[i]) / dt
            }),
            core::array::from_fn(|i| {
                (last[i] + mass * gravity[i]) / mass - (after[i] - ambient[i]) / dt
            }),
            Moments::new(
                self.translation_resistance_kg_s,
                self.properties.mass(),
                duration_s,
            )?,
        )
    }
}

fn validate_flow_endpoint(velocity: Vec3) -> Result<(), LinearDragError> {
    if velocity.is_finite() {
        Ok(())
    } else {
        Err(DynamicsError::NonFinite("linear_drag.ambient_velocity_after_world").into())
    }
}
