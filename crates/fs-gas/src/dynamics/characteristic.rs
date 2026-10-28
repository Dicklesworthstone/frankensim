//! Compact rigid-reservoir interfaces for the original frozen caloric source.
//!
//! An outgoing duct wave supplies one compatibility relation. Subsonic inflow
//! accelerates isentropically from the reservoir stagnation state; outflow
//! sees its static backpressure and carries the duct's total enthalpy into the
//! reservoir. The ideal contraction may choke at a sonic inlet. A stationary
//! ghost passed to a two-sided flux would impose an extra acoustic impedance.
//!
//! These are estimated classical-wave boundary closures. No discharge loss,
//! finite neck inertia, multidimensional mixing, or entropy theorem is added.

mod wave;

use super::{Conserved, Error, FaceFlux, FrozenGas, Primitive, poll, positive};
use fs_math::det;
use fs_qty::{Pressure, Temperature};
use wave::{cold_outflow, pressure_wave, sample_outflow};

impl FrozenGas {
    /// One-sided flux at a compact, stationary reservoir attached to a duct.
    ///
    /// Both states use this exact frozen source and energy reference. Reservoir
    /// momentum must be zero. Its pressure/temperature are stagnation conditions
    /// for incoming gas; gas leaving the duct sees its backpressure. The outward
    /// normal is -1 at a left end and +1 at a right end; returned flux is always
    /// oriented along increasing duct coordinate. No state is extrapolated beyond
    /// the declared caloric interval. A supported sonic inlet limits inflow.
    pub fn reservoir_flux(
        &self,
        interior: Conserved,
        reservoir: Conserved,
        outward_normal: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FaceFlux, Error> {
        let inside = self.decode(interior, None, &mut cancelled)?;
        let outside = self.decode(reservoir, None, &mut cancelled)?;
        self.plenum_flux(
            interior,
            inside,
            reservoir,
            outside,
            outward_normal,
            &mut cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn plenum_flux(
        &self,
        interior: Conserved,
        inside: Primitive,
        reservoir: Conserved,
        outside: Primitive,
        outward_normal: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<FaceFlux, Error> {
        poll(cancelled)?;
        if !matches!(outward_normal, -1.0 | 1.0)
            || reservoir.momentum != 0.0
            || !self.admissible_state(interior)
            || !self.admissible_state(reservoir)
        {
            return Err(Error::Invalid(
                "reservoir normal, stationary state or energy reference",
            ));
        }
        let inside = Primitive {
            velocity: outward_normal * inside.velocity,
            ..inside
        };
        if inside.velocity >= inside.sound_speed && outside.pressure <= inside.pressure {
            // Every expansion characteristic departs the interval. Its unused
            // far-field star need not lie in the interior source's caloric window.
            return self.outlet_flux(inside, outward_normal, inside.velocity + inside.sound_speed);
        }
        let stopped = match pressure_wave(self, inside, outside.pressure, cancelled) {
            Ok(wave) => Some(wave),
            Err(error) if hot_trial(&error) && outside.pressure > inside.pressure => None,
            Err(error) if cold_trial(&error) && outside.pressure < inside.pressure => {
                let (port, signal) = cold_outflow(self, inside, cancelled)?;
                return self.outlet_flux(port, outward_normal, signal);
            }
            Err(error) => return Err(error),
        };
        let scale = inside.velocity.abs() + inside.sound_speed + outside.sound_speed;
        let velocity_tolerance = self.domain().absolute_energy_tolerance / scale
            + (self.domain().relative_energy_tolerance + 64.0 * f64::EPSILON) * scale;
        if let Some(stopped) = stopped.filter(|wave| wave.primitive.velocity >= 0.0) {
            // The outgoing shock/fan may lie wholly outside the domain. Sampling
            // its x/t=0 state prevents a supersonic exit being overprescribed.
            let port = sample_outflow(self, inside, stopped, cancelled)?;
            return self.outlet_flux(port, outward_normal, stopped.signal_speed);
        }

        // For inflow, solve v + u_i - f_i(p(v)) = 0 with v >= 0.
        // dh = -v dv and ds=0 give dp/dv=-rho*v, so its derivative is
        // 1 + rho*v*f_i'(p) > 0 on an admitted classical wave. The bracket
        // v=-u_i+f_i(p0) normally encloses the root without a sonic solve.
        let source_drop =
            self.enthalpy_difference(self.domain().minimum_temperature, outside.temperature)?;
        if !positive(source_drop) {
            return Err(Error::CaloricPath(
                "reservoir inflow has no supported cooling interval",
            ));
        }
        let source_speed = det::sqrt(2.0 * source_drop);
        if !positive(source_speed) {
            return Err(Error::Invalid("reservoir inflow speed is unrepresentable"));
        }
        let (mut hi, upper, choked) = if let Some(stopped) = stopped {
            let speed = (-stopped.primitive.velocity).min(source_speed);
            let port = self.inlet_state(outside, speed, cancelled)?;
            if speed < port.sound_speed {
                (speed, port, false)
            } else {
                self.inlet_limit(outside, cancelled)?
            }
        } else {
            // The stopped shock at p0 is only an auxiliary state. Isentropic
            // acceleration lowers the actual port pressure; a valid match may
            // exist even if that unused stopped shock would exceed Tmax.
            self.inlet_limit(outside, cancelled)?
        };
        let upper_wave = pressure_wave(self, inside, upper.pressure, cancelled)?;
        let upper_residual = hi + inside.velocity - upper_wave.jump;
        if upper_residual < -velocity_tolerance {
            if !choked {
                return Err(Error::CaloricPath(
                    "reservoir characteristic leaves supported source path",
                ));
            }
            return self.inlet_flux(
                outside,
                upper,
                outward_normal,
                upper_wave
                    .signal_speed
                    .max(inside.velocity.abs() + inside.sound_speed),
                cancelled,
            );
        }
        if upper_residual.abs() <= velocity_tolerance {
            return self.inlet_flux(
                outside,
                upper,
                outward_normal,
                upper_wave
                    .signal_speed
                    .max(inside.velocity.abs() + inside.sound_speed),
                cancelled,
            );
        }
        let mut lo = 0.0;
        let mut speed = hi;
        for _ in 0..self.domain().maximum_iterations {
            poll(cancelled)?;
            let port = self.inlet_state(outside, speed, cancelled)?;
            let wave = match pressure_wave(self, inside, port.pressure, cancelled) {
                Ok(wave) => wave,
                Err(error) if hot_trial(&error) && port.pressure > inside.pressure => {
                    // Search toward lower pressure without accepting the refused
                    // trial as a physical state. Acceptance still requires an
                    // admitted wave with the actual compatibility residual small.
                    lo = speed;
                    speed = lo.midpoint(hi);
                    continue;
                }
                Err(error) => return Err(error),
            };
            let residual = speed + inside.velocity - wave.jump;
            if residual.abs() <= velocity_tolerance {
                return self.inlet_flux(
                    outside,
                    port,
                    outward_normal,
                    wave.signal_speed
                        .max(inside.velocity.abs() + inside.sound_speed),
                    cancelled,
                );
            }
            if residual > 0.0 {
                hi = speed;
            } else {
                lo = speed;
            }
            let derivative = 1.0 + port.density * speed * wave.derivative;
            if !positive(derivative) {
                return Err(Error::CaloricPath("nonclassical reservoir characteristic"));
            }
            let trial = speed - residual / derivative;
            speed = if trial > lo && trial < hi {
                trial
            } else {
                lo.midpoint(hi)
            };
        }
        Err(Error::TemperatureIterations)
    }

    fn outlet_flux(&self, port: Primitive, normal: f64, signal: f64) -> Result<FaceFlux, Error> {
        let q = self.conserved(port.density, port.velocity, port.pressure)?;
        checked_flux(
            normal * q.momentum,
            q.momentum * port.velocity + port.pressure,
            normal * (q.energy + port.pressure) * port.velocity,
            signal.max(port.velocity.abs() + port.sound_speed),
        )
    }

    /// Sonic limit, or the cold source endpoint if that is reached first.
    /// Temperature bracketing avoids evaluating an unused extreme low-pressure
    /// isentrope before locating a supported sonic throat.
    fn inlet_limit(
        &self,
        reservoir: Primitive,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(f64, Primitive, bool), Error> {
        let mut lo = self.domain().minimum_temperature;
        let mut hi = reservoir.temperature;
        let at_temperature = |t: f64| -> Result<(f64, f64), Error> {
            let kinetic_twice = 2.0 * self.enthalpy_difference(t, reservoir.temperature)?;
            let properties = self
                .source()
                .evaluate_properties(Temperature::new(t), Pressure::new(reservoir.pressure))?;
            let acoustic_squared = (properties.mass_specific_cp().value()
                / properties.mass_specific_cv().value())
                * self.gas_constant()
                * t;
            if !kinetic_twice.is_finite() || !positive(acoustic_squared) {
                return Err(Error::CaloricPath("reservoir sonic source state"));
            }
            Ok((kinetic_twice, acoustic_squared))
        };
        let (kinetic, acoustic) = at_temperature(lo)?;
        if kinetic < acoustic {
            let speed = det::sqrt(kinetic);
            let port = self.inlet_at_temperature(reservoir, lo, speed)?;
            return Ok((speed, port, false));
        }
        for _ in 0..self.domain().maximum_iterations {
            poll(cancelled)?;
            let t = lo.midpoint(hi);
            let (kinetic, acoustic) = at_temperature(t)?;
            let tolerance = self.domain().absolute_energy_tolerance
                + (self.domain().relative_energy_tolerance + 64.0 * f64::EPSILON)
                    * kinetic.max(acoustic);
            let residual = kinetic - acoustic;
            if residual.abs() <= tolerance {
                let speed = det::sqrt(kinetic);
                let port = self.inlet_at_temperature(reservoir, t, speed)?;
                return Ok((speed, port, true));
            }
            if residual > 0.0 {
                lo = t;
            } else {
                hi = t;
            }
        }
        Err(Error::TemperatureIterations)
    }

    fn enthalpy_difference(&self, low: f64, high: f64) -> Result<f64, Error> {
        let value = self
            .source()
            .specific_internal_energy_change(Temperature::new(low), Temperature::new(high))?
            .value()
            + self.gas_constant() * (high - low);
        if !value.is_finite() || value < 0.0 {
            return Err(Error::CaloricPath("reservoir source enthalpy difference"));
        }
        Ok(value)
    }

    /// Static state from reservoir stagnation enthalpy and entropy at speed v.
    fn inlet_state(
        &self,
        reservoir: Primitive,
        speed: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Primitive, Error> {
        let target = 0.5 * speed * speed;
        if !speed.is_finite() || speed < 0.0 || !target.is_finite() {
            return Err(Error::Invalid("reservoir inlet kinetic energy"));
        }
        let mut lo = self.domain().minimum_temperature;
        let mut hi = reservoir.temperature;
        let available = self.enthalpy_difference(lo, hi)?;
        let tolerance = self.domain().absolute_energy_tolerance
            + self.domain().relative_energy_tolerance * target;
        if target > available + 8.0 * f64::EPSILON * available {
            return Err(Error::CaloricPath(
                "reservoir acceleration leaves cold source domain",
            ));
        }
        let cp0 = self
            .source()
            .evaluate_properties(Temperature::new(hi), Pressure::new(reservoir.pressure))?
            .mass_specific_cp()
            .value();
        let mut t = (hi - target / cp0).clamp(lo, hi);
        for _ in 0..self.domain().maximum_iterations {
            poll(cancelled)?;
            let residual = self.enthalpy_difference(t, reservoir.temperature)? - target;
            if residual.abs() <= tolerance {
                return self.inlet_at_temperature(reservoir, t, speed);
            }
            if residual > 0.0 {
                lo = t;
            } else {
                hi = t;
            }
            let cp = self
                .source()
                .evaluate_properties(Temperature::new(t), Pressure::new(reservoir.pressure))?
                .mass_specific_cp()
                .value();
            let trial = t + residual / cp;
            t = if trial > lo && trial < hi {
                trial
            } else {
                lo.midpoint(hi)
            };
        }
        Err(Error::TemperatureIterations)
    }

    fn inlet_at_temperature(
        &self,
        reservoir: Primitive,
        t: f64,
        speed: f64,
    ) -> Result<Primitive, Error> {
        let entropy = self
            .source()
            .specific_entropy_change(
                Temperature::new(reservoir.temperature),
                Pressure::new(reservoir.pressure),
                Temperature::new(t),
                Pressure::new(reservoir.pressure),
            )?
            .value();
        let pressure = reservoir.pressure * det::exp(entropy / self.gas_constant());
        let density = (pressure / self.gas_constant()) / t;
        self.from_rho_u_p(density, -speed, pressure)
    }

    fn inlet_flux(
        &self,
        reservoir: Primitive,
        port: Primitive,
        outward_normal: f64,
        signal_speed: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<FaceFlux, Error> {
        poll(cancelled)?;
        // The upstream reservoir supplies total enthalpy once. Use its original
        // sensible enthalpy directly, retaining the small thermal difference
        // independently of both the formation gauge and the kinetic head.
        let enthalpy =
            self.sensible(reservoir.temperature)? + self.gas_constant() * reservoir.temperature;
        self.conserved(port.density, port.velocity, port.pressure)?;
        let mass = outward_normal * port.density * port.velocity;
        checked_flux(
            mass,
            port.density * port.velocity * port.velocity + port.pressure,
            mass * enthalpy,
            signal_speed.max(port.velocity.abs() + port.sound_speed),
        )
    }
}

fn hot_trial(error: &Error) -> bool {
    matches!(
        error,
        Error::CaloricPath(
            "characteristic shock exceeds the hot source path"
                | "characteristic shock state exceeds source temperature maximum"
        )
    )
}
fn cold_trial(error: &Error) -> bool {
    matches!(
        error,
        Error::CaloricPath(
            "characteristic rarefaction exceeds the cold source path"
                | "characteristic rarefaction state below source temperature minimum"
        )
    )
}

fn checked_flux(
    mass: f64,
    momentum: f64,
    energy: f64,
    signal_speed: f64,
) -> Result<FaceFlux, Error> {
    if !mass.is_finite() || !positive(momentum) || !energy.is_finite() || !positive(signal_speed) {
        return Err(Error::Invalid("non-finite reservoir flux"));
    }
    Ok(FaceFlux {
        mass,
        momentum,
        energy,
        signal_speed,
        used_hll_fallback: false,
    })
}

#[cfg(test)]
mod tests;
