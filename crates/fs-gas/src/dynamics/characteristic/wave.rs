//! One-sided thermally perfect pressure waves, in the outward axial coordinate.
//!
//! The retained source supplies the shock Hugoniot and rarefaction entropy and
//! invariant. Classical compressive shocks/expansion fans are admitted here;
//! no gamma-law replacement or unsupported-temperature extension is used.

use super::super::{Error, FrozenGas, Primitive, poll, positive};
use fs_math::det;
use fs_qty::{Pressure, Temperature};

#[derive(Clone, Copy, Debug)]
pub(super) struct PressureWave {
    pub(super) primitive: Primitive,
    /// Signed velocity decrement f(p): positive for compression.
    pub(super) jump: f64,
    pub(super) derivative: f64,
    pub(super) signal_speed: f64,
    /// Shock speed, or the rarefaction head u_i-c_i.
    pub(super) head_speed: f64,
    /// Shock speed, or the rarefaction tail u_star-c_star.
    pub(super) tail_speed: f64,
}

pub(super) fn pressure_wave(
    gas: &FrozenGas,
    inside: Primitive,
    target_pressure: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<PressureWave, Error> {
    poll(cancelled)?;
    if !positive(target_pressure)
        || !positive(inside.density)
        || !positive(inside.pressure)
        || !positive(inside.temperature)
        || !positive(inside.sound_speed)
        || !inside.velocity.is_finite()
        || !(gas.domain().minimum_temperature..=gas.domain().maximum_temperature)
            .contains(&inside.temperature)
    {
        return Err(Error::Invalid("characteristic pressure-wave input"));
    }
    if target_pressure == inside.pressure {
        return checked_wave(PressureWave {
            primitive: inside,
            jump: 0.0,
            derivative: 1.0 / (inside.density * inside.sound_speed),
            signal_speed: inside.velocity.abs() + inside.sound_speed,
            head_speed: inside.velocity - inside.sound_speed,
            tail_speed: inside.velocity - inside.sound_speed,
        });
    }
    if target_pressure > inside.pressure {
        compression(gas, inside, target_pressure, cancelled)
    } else {
        expansion(gas, inside, target_pressure, cancelled)
    }
}

fn compression(
    gas: &FrozenGas,
    inside: Primitive,
    pressure: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<PressureWave, Error> {
    let r = gas.gas_constant();
    let ti = inside.temperature;
    let dp = pressure - inside.pressure;
    let vi = 1.0 / inside.density;
    let alpha = dp / pressure;
    let cp_i = cp(gas, ti, inside.pressure)?;
    let maximum = gas.domain().maximum_temperature;
    if ti >= maximum || !positive(vi) {
        return Err(Error::CaloricPath(
            "characteristic shock exceeds the hot source path",
        ));
    }
    // z = Delta T / Delta p. Factoring out Delta p retains the acoustic
    // weak-shock limit even when adding Delta T to Ti rounds back to Ti.
    // Hugoniot: z*(Cp_secant - .5*alpha*R) = .5*(vi + R*Ti/pb).
    let rhs = 0.5 * (vi + r * ti / pressure);
    let residual_tolerance = (gas.domain().relative_energy_tolerance + 64.0 * f64::EPSILON) * vi;
    let evaluate = |z: f64| -> Result<(f64, f64, f64), Error> {
        let raw_t = ti + dp * z;
        // z's upper bound maps to Tmax; remove only endpoint rounding.
        let t = if raw_t > maximum && raw_t - maximum <= 4.0 * f64::EPSILON * maximum {
            maximum
        } else {
            raw_t
        };
        if !t.is_finite() || t < ti || t > maximum {
            return Err(Error::CaloricPath("characteristic shock temperature path"));
        }
        let delta_t = t - ti;
        let cp_secant = if delta_t == 0.0 {
            cp_i
        } else {
            gas.source()
                .specific_internal_energy_change(Temperature::new(ti), Temperature::new(t))?
                .value()
                / delta_t
                + r
        };
        let derivative = cp(gas, t, pressure)? - 0.5 * alpha * r;
        let residual = z * (cp_secant - 0.5 * alpha * r) - rhs;
        if !positive(cp_secant) || !positive(derivative) || !residual.is_finite() {
            return Err(Error::CaloricPath(
                "characteristic shock caloric derivative",
            ));
        }
        Ok((residual, derivative, t))
    };
    let mut lo = 0.0;
    // Positive Cv implies z < vi/R; this also bounds very weak-shock work.
    let mut hi = ((maximum - ti) / dp).min(vi / r);
    if !positive(hi) || evaluate(hi)?.0 < -residual_tolerance {
        return Err(Error::CaloricPath(
            "characteristic shock state exceeds source temperature maximum",
        ));
    }
    let mut z = (rhs / (cp_i - 0.5 * alpha * r)).clamp(lo, hi);
    for _ in 0..gas.domain().maximum_iterations {
        poll(cancelled)?;
        let (residual, thermal_derivative, t) = evaluate(z)?;
        if residual.abs() <= residual_tolerance {
            let factor = (vi - r * z) / pressure;
            if !positive(factor) {
                return Err(Error::CaloricPath("noncompressive pressure-shock branch"));
            }
            let root = det::sqrt(factor);
            let jump = dp * root;
            let velocity = inside.velocity - jump;
            let density = (pressure / r) / t;
            let primitive = gas.from_rho_u_p(density, velocity, pressure)?;
            let vb = r * t / pressure;
            let dt_dp = 0.5 * (vi + vb * (inside.pressure / pressure)) / thermal_derivative;
            let dv_dp = (r / pressure) * (dt_dp - t / pressure);
            let derivative = (factor - dv_dp) / (2.0 * root);
            let speed = inside.velocity - vi / root;
            let tolerance = velocity_tolerance(gas, inside, primitive);
            if !positive(derivative)
                || !speed.is_finite()
                || speed > inside.velocity - inside.sound_speed + tolerance
                || speed < primitive.velocity - primitive.sound_speed - tolerance
                || speed >= primitive.velocity
            {
                return Err(Error::CaloricPath(
                    "source shock does not satisfy the classical compressive wave ordering",
                ));
            }
            let ds = gas
                .source()
                .specific_entropy_change(
                    Temperature::new(ti),
                    Pressure::new(inside.pressure),
                    Temperature::new(t),
                    Pressure::new(pressure),
                )?
                .value();
            let entropy_tolerance = gas.domain().absolute_energy_tolerance / ti
                + gas.domain().relative_energy_tolerance * dp * vi / ti
                + 128.0 * f64::EPSILON * cp_i.max(cp(gas, t, pressure)?);
            if !ds.is_finite() || ds < -entropy_tolerance {
                return Err(Error::CaloricPath(
                    "characteristic shock has inadmissible entropy production",
                ));
            }
            return checked_wave(PressureWave {
                primitive,
                jump,
                derivative,
                signal_speed: (inside.velocity.abs() + inside.sound_speed)
                    .max(primitive.velocity.abs() + primitive.sound_speed)
                    .max(speed.abs()),
                head_speed: speed,
                tail_speed: speed,
            });
        }
        if residual < 0.0 {
            lo = z;
        } else {
            hi = z;
        }
        let trial = z - residual / thermal_derivative;
        z = if trial > lo && trial < hi {
            trial
        } else {
            lo.midpoint(hi)
        };
    }
    Err(Error::TemperatureIterations)
}

fn expansion(
    gas: &FrozenGas,
    inside: Primitive,
    pressure: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<PressureWave, Error> {
    let ti = inside.temperature;
    let r = gas.gas_constant();
    let cp_i = cp(gas, ti, inside.pressure)?;
    // Obtain log(pi/pb) from the source's cancellation-safe pressure entropy
    // change, rather than forming an overflowing or rounded pressure ratio.
    let log_ratio = gas
        .source()
        .specific_entropy_change(
            Temperature::new(ti),
            Pressure::new(inside.pressure),
            Temperature::new(ti),
            Pressure::new(pressure),
        )?
        .value()
        / r;
    if !positive(log_ratio) || ti <= gas.domain().minimum_temperature {
        return Err(Error::CaloricPath(
            "characteristic rarefaction exceeds the cold source path",
        ));
    }
    // z = (Ti-Tb)/log(pi/pb), using the source Cp/T secant. This preserves
    // the nonzero weak-wave increment if Tb rounds to Ti.
    let evaluate = |z: f64| -> Result<(f64, f64, f64), Error> {
        let raw_t = ti - log_ratio * z;
        let minimum = gas.domain().minimum_temperature;
        let t = if raw_t < minimum && minimum - raw_t <= 4.0 * f64::EPSILON * ti {
            minimum
        } else {
            raw_t
        };
        if !t.is_finite() || t < minimum || t > ti {
            return Err(Error::CaloricPath(
                "characteristic rarefaction temperature path",
            ));
        }
        let delta_t = ti - t;
        let secant = if delta_t == 0.0 {
            cp_i / ti
        } else {
            -gas.source()
                .specific_entropy_change(
                    Temperature::new(ti),
                    Pressure::new(inside.pressure),
                    Temperature::new(t),
                    Pressure::new(inside.pressure),
                )?
                .value()
                / delta_t
        };
        let derivative = cp(gas, t, pressure)? / t;
        let residual = secant * z - r;
        if !positive(secant) || !positive(derivative) || !residual.is_finite() {
            return Err(Error::CaloricPath(
                "characteristic rarefaction caloric derivative",
            ));
        }
        Ok((residual, derivative, t))
    };
    let residual_tolerance = (gas.domain().relative_energy_tolerance + 64.0 * f64::EPSILON) * r;
    let mut lo = 0.0;
    let mut hi = ((ti - gas.domain().minimum_temperature) / log_ratio).min(ti);
    if !positive(hi) || evaluate(hi)?.0 < -residual_tolerance {
        return Err(Error::CaloricPath(
            "characteristic rarefaction state below source temperature minimum",
        ));
    }
    let mut z = (r * ti / cp_i).clamp(lo, hi);
    for _ in 0..gas.domain().maximum_iterations {
        poll(cancelled)?;
        let (residual, derivative, t) = evaluate(z)?;
        if residual.abs() <= residual_tolerance {
            let delta_t = log_ratio * z;
            let tolerance = base_velocity_tolerance(gas, inside);
            let integral = if t == ti {
                delta_t * (cp_i / inside.sound_speed)
            } else {
                let integral = gas.rarefaction_integral(t, ti, 0.25 * tolerance, cancelled)?;
                if integral == 0.0 {
                    // The square-root quadrature coordinates may coincide for
                    // a one-ulp temperature interval. Its source tangent is
                    // the representable weak-wave limit.
                    delta_t * (cp_i / inside.sound_speed)
                } else {
                    integral * (delta_t / (ti - t))
                }
            };
            let primitive =
                gas.from_rho_u_p((pressure / r) / t, inside.velocity + integral, pressure)?;
            let head = inside.velocity - inside.sound_speed;
            let tail = primitive.velocity - primitive.sound_speed;
            if tail < head - velocity_tolerance(gas, inside, primitive) {
                return Err(Error::CaloricPath(
                    "source rarefaction does not satisfy classical fan ordering",
                ));
            }
            return checked_wave(PressureWave {
                primitive,
                jump: -integral,
                derivative: 1.0 / (primitive.density * primitive.sound_speed),
                signal_speed: (inside.velocity.abs() + inside.sound_speed)
                    .max(primitive.velocity.abs() + primitive.sound_speed),
                head_speed: head,
                tail_speed: tail,
            });
        }
        if residual < 0.0 {
            lo = z;
        } else {
            hi = z;
        }
        let trial = z - residual / derivative;
        z = if trial > lo && trial < hi {
            trial
        } else {
            lo.midpoint(hi)
        };
    }
    Err(Error::TemperatureIterations)
}

/// Sample the one-sided pressure wave at the fixed interface x/t=0.
/// A supersonic departing wave leaves the interior trace untouched; when a
/// rarefaction straddles the interface, its sonic state supplies the flux.
pub(super) fn sample_outflow(
    gas: &FrozenGas,
    inside: Primitive,
    wave: PressureWave,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Primitive, Error> {
    poll(cancelled)?;
    if wave.head_speed >= 0.0 {
        return Ok(inside);
    }
    if wave.tail_speed <= 0.0 {
        return Ok(wave.primitive);
    }
    if wave.primitive.pressure >= inside.pressure {
        return Err(Error::CaloricPath("unsupported straddling pressure shock"));
    }
    let mut lo = wave.primitive.temperature;
    let mut hi = inside.temperature;
    let tolerance = velocity_tolerance(gas, inside, wave.primitive);
    for _ in 0..gas.domain().maximum_iterations {
        poll(cancelled)?;
        let t = lo.midpoint(hi);
        let primitive = fan_state(gas, inside, t, 0.25 * tolerance, cancelled)?;
        let residual = primitive.velocity - primitive.sound_speed;
        if residual.abs() <= tolerance {
            return Ok(primitive);
        }
        if residual > 0.0 {
            lo = t;
        } else {
            hi = t;
        }
        if lo >= hi || lo.midpoint(hi) == lo || lo.midpoint(hi) == hi {
            return Err(Error::TemperatureIterations);
        }
    }
    Err(Error::TemperatureIterations)
}

/// Recover a supported sonic outlet when the requested downstream pressure
/// would extend an unused part of the expansion fan below the source window.
/// The admitted cold endpoint must already lie beyond the sonic characteristic;
/// otherwise even the interface state requires unsupported colder gas.
pub(super) fn cold_outflow(
    gas: &FrozenGas,
    inside: Primitive,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Primitive, f64), Error> {
    poll(cancelled)?;
    let minimum = gas.domain().minimum_temperature;
    let tolerance = base_velocity_tolerance(gas, inside);
    let endpoint = fan_state(gas, inside, minimum, 0.25 * tolerance, cancelled)?;
    let head = inside.velocity - inside.sound_speed;
    let tail = endpoint.velocity - endpoint.sound_speed;
    if tail < 0.0 {
        return Err(Error::CaloricPath(
            "sonic outflow state lies below the source temperature minimum",
        ));
    }
    let wave = checked_wave(PressureWave {
        primitive: endpoint,
        jump: inside.velocity - endpoint.velocity,
        derivative: 1.0 / (endpoint.density * endpoint.sound_speed),
        signal_speed: (inside.velocity.abs() + inside.sound_speed)
            .max(endpoint.velocity.abs() + endpoint.sound_speed),
        head_speed: head,
        tail_speed: tail,
    })?;
    let sampled = sample_outflow(gas, inside, wave, cancelled)?;
    Ok((
        sampled,
        wave.signal_speed
            .max(sampled.velocity.abs() + sampled.sound_speed),
    ))
}

fn fan_state(
    gas: &FrozenGas,
    inside: Primitive,
    temperature: f64,
    tolerance: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Primitive, Error> {
    let ds = gas
        .source()
        .specific_entropy_change(
            Temperature::new(inside.temperature),
            Pressure::new(inside.pressure),
            Temperature::new(temperature),
            Pressure::new(inside.pressure),
        )?
        .value();
    let pressure = inside.pressure * det::exp(ds / gas.gas_constant());
    let velocity = inside.velocity
        + gas.rarefaction_integral(temperature, inside.temperature, tolerance, cancelled)?;
    gas.from_rho_u_p(
        (pressure / gas.gas_constant()) / temperature,
        velocity,
        pressure,
    )
}

fn cp(gas: &FrozenGas, temperature: f64, pressure: f64) -> Result<f64, Error> {
    let value = gas
        .source()
        .evaluate_properties(Temperature::new(temperature), Pressure::new(pressure))?
        .mass_specific_cp()
        .value();
    if positive(value) {
        Ok(value)
    } else {
        Err(Error::CaloricPath("characteristic wave has nonpositive Cp"))
    }
}
fn base_velocity_tolerance(gas: &FrozenGas, state: Primitive) -> f64 {
    let speed = state.velocity.abs() + state.sound_speed;
    gas.domain().absolute_energy_tolerance / speed
        + (gas.domain().relative_energy_tolerance + 128.0 * f64::EPSILON) * speed
}
fn velocity_tolerance(gas: &FrozenGas, inside: Primitive, star: Primitive) -> f64 {
    base_velocity_tolerance(gas, inside).max(base_velocity_tolerance(gas, star))
}
fn checked_wave(wave: PressureWave) -> Result<PressureWave, Error> {
    if !wave.jump.is_finite()
        || !positive(wave.derivative)
        || !positive(wave.signal_speed)
        || !wave.head_speed.is_finite()
        || !wave.tail_speed.is_finite()
    {
        return Err(Error::Invalid(
            "unrepresentable characteristic pressure wave",
        ));
    }
    Ok(wave)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dynamics::tests::gas;

    fn relative(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance * expected.abs().max(1.0),
            "actual={actual:.16e}, expected={expected:.16e}"
        );
    }

    #[test]
    fn variable_cp_wave_matches_independent_polynomial_and_quadrature() {
        let slope = 0.15;
        let gas = gas(slope, 0.0);
        let inside = gas.from_rho_u_p(1.0, 0.1, 1.0).unwrap();
        for pressure in [1.7, 3.0, 9.0] {
            let wave = pressure_wave(&gas, inside, pressure, &mut || false).unwrap();
            // Independent caloric law h(T)-h(1)=3.5(T-1)+.075(T²-1).
            let a = 0.5 * slope;
            let b = 3.5 - 0.5 * (pressure - 1.0) / pressure;
            let c = -(3.5 + a + 0.5 * (pressure - 1.0));
            let t = -2.0 * c / (b + det::sqrt(b * b - 4.0 * a * c));
            let jump = det::sqrt((pressure - 1.0) * (1.0 - t / pressure));
            relative(wave.primitive.temperature, t, 3e-12);
            relative(wave.jump, jump, 3e-12);
            let entropy = 3.5 * det::ln(t) + slope * (t - 1.0) - det::ln(pressure);
            assert!(entropy > 0.0);
        }
        let pressure = 0.15;
        let wave = pressure_wave(&gas, inside, pressure, &mut || false).unwrap();
        // Independent entropy equation of this analytic Cp(T) fixture.
        let (mut lo, mut hi) = (0.001_f64, 1.0_f64);
        for _ in 0..80 {
            let t = lo.midpoint(hi);
            let residual = 3.5 * det::ln(t) + slope * (t - 1.0) - det::ln(pressure);
            if residual < 0.0 {
                lo = t;
            } else {
                hi = t;
            }
        }
        let t = lo.midpoint(hi);
        let panels = 2000;
        let dx = (1.0 - t) / panels as f64;
        let integrand = |temperature: f64| {
            let cp = 3.5 + slope * temperature;
            cp / det::sqrt(cp / (cp - 1.0) * temperature)
        };
        let mut sum = integrand(t) + integrand(1.0);
        for i in 1..panels {
            sum += if i % 2 == 0 { 2.0 } else { 4.0 } * integrand(t + i as f64 * dx);
        }
        relative(wave.primitive.temperature, t, 3e-12);
        relative(wave.jump, -dx * sum / 3.0, 3e-10);
    }
}
