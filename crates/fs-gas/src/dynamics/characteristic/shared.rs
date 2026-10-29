//! Characteristic compact-volume port over two borrowed local compositions.
//!
//! The duct law owns the outgoing pressure wave; the reservoir law owns the
//! stagnation-to-static acceleration and choking condition. Their shared
//! pressure and normal velocity couple across the advected composition contact.
//! Owners form the conservative donor flux after this physical state is admitted.

use super::wave::{cold_outflow, pressure_wave, sample_outflow};
use crate::dynamics::wall::shared::{CaloricLaw, WallError as Error};
use crate::dynamics::{Primitive, positive};
use fs_math::det;

/// Shared converging-nozzle closure for finite compact-volume valves.
pub(crate) mod nozzle;

#[derive(Clone, Copy)]
pub(crate) struct PortState {
    /// Velocity is oriented outward from the duct, independent of its axis.
    pub primitive: Primitive,
    /// The reservoir supplies composition and total enthalpy when true.
    pub inflow: bool,
    pub signal_speed: f64,
}

pub(crate) fn reservoir_state<L: CaloricLaw>(
    duct: &L,
    reservoir: &L,
    inside: Primitive,
    outside: Primitive,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<PortState, Error<L::Error>> {
    poll(cancelled)?;
    if inside.velocity >= inside.sound_speed && outside.pressure <= inside.pressure {
        return port(inside, false, inside.velocity + inside.sound_speed);
    }
    let stopped = match pressure_wave(duct, inside, outside.pressure, cancelled) {
        Ok(wave) => Some(wave),
        Err(error) if hot_trial(&error) && outside.pressure > inside.pressure => None,
        Err(error) if cold_trial(&error) && outside.pressure < inside.pressure => {
            let (state, signal) = cold_outflow(duct, inside, cancelled)?;
            return port(state, false, signal);
        }
        Err(error) => return Err(error),
    };
    let scale = inside.velocity.abs() + inside.sound_speed + outside.sound_speed;
    let velocity_tolerance = duct.domain().absolute_energy_tolerance / scale
        + (duct.domain().relative_energy_tolerance + 64.0 * f64::EPSILON) * scale;
    if let Some(stopped) = stopped.filter(|wave| wave.primitive.velocity >= 0.0) {
        let state = sample_outflow(duct, inside, stopped, cancelled)?;
        return port(state, false, stopped.signal_speed);
    }

    // Inflow compatibility: v + u_i - f_i(p(v)) = 0, v >= 0. The reservoir
    // isentrope gives dp/dv=-rho_res*v while the outgoing wave uses duct Cp/Cv.
    // Thus the derivative remains 1+rho_res*v*f_i'(p), even across composition.
    let source_drop = enthalpy_difference(
        reservoir,
        reservoir.domain().minimum_temperature,
        outside.temperature,
        cancelled,
    )?;
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
        let speed = -stopped.primitive.velocity;
        if speed >= source_speed {
            // This auxiliary compatibility bound exhausts the entire caloric
            // path. Find the supported sonic throat before constructing an
            // unused cold, supersonic conservative state at Tmin.
            inlet_limit(reservoir, outside, cancelled)?
        } else {
            let state = inlet_state(reservoir, outside, speed, cancelled)?;
            if speed < state.sound_speed {
                (speed, state, false)
            } else {
                inlet_limit(reservoir, outside, cancelled)?
            }
        }
    } else {
        // An unsupported auxiliary stopped shock does not disqualify a
        // supported lower-pressure accelerating inlet.
        inlet_limit(reservoir, outside, cancelled)?
    };
    let upper_wave = pressure_wave(duct, inside, upper.pressure, cancelled)?;
    let upper_residual = hi + inside.velocity - upper_wave.jump;
    let signal = upper_wave
        .signal_speed
        .max(inside.velocity.abs() + inside.sound_speed);
    if upper_residual < -velocity_tolerance {
        if !choked {
            return Err(Error::CaloricPath(
                "reservoir characteristic leaves supported source path",
            ));
        }
        return port(upper, true, signal);
    }
    if upper_residual.abs() <= velocity_tolerance {
        return port(upper, true, signal);
    }
    let mut lo = 0.0;
    let mut speed = hi;
    for _ in 0..duct.domain().maximum_iterations {
        poll(cancelled)?;
        let state = inlet_state(reservoir, outside, speed, cancelled)?;
        let wave = match pressure_wave(duct, inside, state.pressure, cancelled) {
            Ok(wave) => wave,
            Err(error) if hot_trial(&error) && state.pressure > inside.pressure => {
                lo = speed;
                speed = lo.midpoint(hi);
                continue;
            }
            Err(error) => return Err(error),
        };
        let residual = speed + inside.velocity - wave.jump;
        if residual.abs() <= velocity_tolerance {
            return port(
                state,
                true,
                wave.signal_speed
                    .max(inside.velocity.abs() + inside.sound_speed),
            );
        }
        if residual > 0.0 {
            hi = speed;
        } else {
            lo = speed;
        }
        let derivative = 1.0 + state.density * speed * wave.derivative;
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
    Err(Error::Iterations)
}

fn inlet_limit<L: CaloricLaw>(
    gas: &L,
    reservoir: Primitive,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(f64, Primitive, bool), Error<L::Error>> {
    let mut lo = gas.domain().minimum_temperature;
    let mut hi = reservoir.temperature;
    let at_temperature = |t: f64, cancelled: &mut _| -> Result<(f64, f64), Error<L::Error>> {
        let kinetic = 2.0 * enthalpy_difference(gas, t, reservoir.temperature, cancelled)?;
        let capacities = gas
            .heat_capacities(t, reservoir.pressure, cancelled)
            .map_err(Error::Source)?;
        let acoustic = (capacities.cp / capacities.cv) * gas.gas_constant() * t;
        if !kinetic.is_finite() || !positive(acoustic) {
            return Err(Error::CaloricPath("reservoir sonic source state"));
        }
        Ok((kinetic, acoustic))
    };
    let (kinetic, acoustic) = at_temperature(lo, cancelled)?;
    if kinetic < acoustic {
        let speed = det::sqrt(kinetic);
        let state = inlet_at_temperature(gas, reservoir, lo, speed, cancelled)?;
        return Ok((speed, state, false));
    }
    for _ in 0..gas.domain().maximum_iterations {
        poll(cancelled)?;
        let t = lo.midpoint(hi);
        let (kinetic, acoustic) = at_temperature(t, cancelled)?;
        let tolerance = gas.domain().absolute_energy_tolerance
            + (gas.domain().relative_energy_tolerance + 64.0 * f64::EPSILON)
                * kinetic.max(acoustic);
        let residual = kinetic - acoustic;
        if residual.abs() <= tolerance {
            let speed = det::sqrt(kinetic);
            let state = inlet_at_temperature(gas, reservoir, t, speed, cancelled)?;
            return Ok((speed, state, true));
        }
        if residual > 0.0 {
            lo = t;
        } else {
            hi = t;
        }
    }
    Err(Error::Iterations)
}

fn enthalpy_difference<L: CaloricLaw>(
    gas: &L,
    low: f64,
    high: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, Error<L::Error>> {
    let value = gas
        .internal_energy_change(low, high, cancelled)
        .map_err(Error::Source)?
        + gas.gas_constant() * (high - low);
    if !value.is_finite() || value < 0.0 {
        return Err(Error::CaloricPath("reservoir source enthalpy difference"));
    }
    Ok(value)
}

fn inlet_state<L: CaloricLaw>(
    gas: &L,
    reservoir: Primitive,
    speed: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Primitive, Error<L::Error>> {
    let target = 0.5 * speed * speed;
    if !speed.is_finite() || speed < 0.0 || !target.is_finite() {
        return Err(Error::Invalid("reservoir inlet kinetic energy"));
    }
    let mut lo = gas.domain().minimum_temperature;
    let mut hi = reservoir.temperature;
    let available = enthalpy_difference(gas, lo, hi, cancelled)?;
    let tolerance =
        gas.domain().absolute_energy_tolerance + gas.domain().relative_energy_tolerance * target;
    if target > available + 8.0 * f64::EPSILON * available {
        return Err(Error::CaloricPath(
            "reservoir acceleration leaves cold source domain",
        ));
    }
    let cp0 = gas
        .heat_capacities(hi, reservoir.pressure, cancelled)
        .map_err(Error::Source)?
        .cp;
    let mut t = (hi - target / cp0).clamp(lo, hi);
    for _ in 0..gas.domain().maximum_iterations {
        poll(cancelled)?;
        let residual = enthalpy_difference(gas, t, reservoir.temperature, cancelled)? - target;
        if residual.abs() <= tolerance {
            return inlet_at_temperature(gas, reservoir, t, speed, cancelled);
        }
        if residual > 0.0 {
            lo = t;
        } else {
            hi = t;
        }
        let cp = gas
            .heat_capacities(t, reservoir.pressure, cancelled)
            .map_err(Error::Source)?
            .cp;
        let trial = t + residual / cp;
        t = if trial > lo && trial < hi {
            trial
        } else {
            lo.midpoint(hi)
        };
    }
    Err(Error::Iterations)
}

fn inlet_at_temperature<L: CaloricLaw>(
    gas: &L,
    reservoir: Primitive,
    t: f64,
    speed: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Primitive, Error<L::Error>> {
    let entropy = gas
        .entropy_change(
            reservoir.temperature,
            reservoir.pressure,
            t,
            reservoir.pressure,
            cancelled,
        )
        .map_err(Error::Source)?;
    let pressure = reservoir.pressure * det::exp(entropy / gas.gas_constant());
    let density = (pressure / gas.gas_constant()) / t;
    gas.primitive(density, -speed, pressure, cancelled)
        .map_err(Error::Source)
}

fn port<E>(primitive: Primitive, inflow: bool, signal: f64) -> Result<PortState, Error<E>> {
    let signal_speed = signal.max(primitive.velocity.abs() + primitive.sound_speed);
    if !positive(primitive.pressure)
        || !positive(primitive.density)
        || !primitive.velocity.is_finite()
        || !positive(signal_speed)
    {
        return Err(Error::Invalid("unrepresentable characteristic port"));
    }
    Ok(PortState {
        primitive,
        inflow,
        signal_speed,
    })
}
fn hot_trial<E>(error: &Error<E>) -> bool {
    matches!(
        error,
        Error::CaloricPath(
            "characteristic shock exceeds the hot source path"
                | "characteristic shock state exceeds source temperature maximum"
        )
    )
}
fn cold_trial<E>(error: &Error<E>) -> bool {
    matches!(
        error,
        Error::CaloricPath(
            "characteristic rarefaction exceeds the cold source path"
                | "characteristic rarefaction state below source temperature minimum"
        )
    )
}
fn poll<E>(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error<E>> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
