//! Incoming acoustic characteristic over a borrowed original caloric source.
//!
//! The reference state is projected onto the interior isentrope. Using a fixed
//! remote temperature instead would spuriously turn an advected entropy contact
//! into an acoustic reflection. No source or workspace allocation is needed.

use super::wave::{pressure_wave, sample_outflow};
use crate::dynamics::wall::shared::{CaloricLaw, WallError as Error, WallSolver};
use crate::dynamics::{Primitive, positive};
use fs_math::det;

pub(crate) fn outlet<L: CaloricLaw>(
    gas: &L,
    inside: Primitive,
    reference_pressure: f64,
    reference_velocity: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Primitive, f64), Error<L::Error>> {
    poll(cancelled)?;
    if !positive(reference_pressure) || !reference_velocity.is_finite() {
        return Err(Error::Invalid("characteristic outlet pressure or velocity"));
    }
    if inside.velocity < 0.0 || reference_velocity < 0.0 {
        return Err(Error::Invalid(
            "characteristic outlet backflow requires an exterior thermodynamic state",
        ));
    }
    let signal = inside.velocity + inside.sound_speed;
    if inside.velocity >= inside.sound_speed
        || (inside.pressure == reference_pressure && inside.velocity == reference_velocity)
    {
        return Ok((inside, signal));
    }
    let scale = signal.max(reference_velocity + inside.sound_speed);
    let tolerance = gas.domain().absolute_energy_tolerance / scale
        + (gas.domain().relative_energy_tolerance + 64.0 * f64::EPSILON) * scale;
    let reference_temperature = isentropic_temperature(gas, inside, reference_pressure, cancelled)?;
    // I(T)=integral Cp(T)/c(T) dT. Only differences are formed: there is no
    // absolute invariant that could depend on a caloric or entropy datum.
    let reference_integral = integral(
        gas,
        inside.temperature,
        reference_temperature,
        0.125 * tolerance,
        cancelled,
    )?;
    // J+ comes from inside; J- comes from the reference on that SAME entropy.
    let residual = (inside.velocity - reference_velocity) + reference_integral;
    if residual.abs() <= tolerance {
        return Ok((inside, signal));
    }
    let velocity = inside.velocity - 0.5 * residual;
    if !velocity.is_finite() || velocity < 0.0 {
        return Err(Error::Invalid(
            "characteristic outlet backflow requires an exterior thermodynamic state",
        ));
    }
    let temperature =
        characteristic_temperature(gas, inside, 0.5 * residual, tolerance, cancelled)?;
    let entropy_change = gas
        .entropy_change(
            inside.temperature,
            inside.pressure,
            temperature,
            inside.pressure,
            cancelled,
        )
        .map_err(Error::Source)?;
    let pressure = inside.pressure * det::exp(entropy_change / gas.gas_constant());
    let density = (pressure / gas.gas_constant()) / temperature;
    if !positive(pressure) || !positive(density) {
        return Err(Error::Invalid(
            "unrepresentable characteristic outlet state",
        ));
    }
    let face = gas
        .primitive(density, velocity, pressure, cancelled)
        .map_err(Error::Source)?;
    if velocity >= face.sound_speed {
        if pressure >= inside.pressure {
            return Err(Error::CaloricPath(
                "nonclassical sonic transition in characteristic outlet",
            ));
        }
        // A subsonic interior can reach sonic discharge along an expansion.
        // Sample that incoming fan at x/t=0 instead of prescribing its remote
        // supersonic tail. The same source owns this already-bounded fan.
        let wave = pressure_wave(gas, inside, pressure, cancelled)?;
        let sonic = sample_outflow(gas, inside, wave, cancelled)?;
        return Ok((sonic, signal.max(wave.signal_speed)));
    }
    Ok((face, signal.max(velocity + face.sound_speed)))
}

fn isentropic_temperature<L: CaloricLaw>(
    gas: &L,
    inside: Primitive,
    pressure: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, Error<L::Error>> {
    if pressure == inside.pressure {
        return Ok(inside.temperature);
    }
    let initial = inside.temperature;
    let domain = gas.domain();
    let (mut low, mut high) = if pressure < inside.pressure {
        (domain.minimum_temperature, initial)
    } else {
        (initial, domain.maximum_temperature)
    };
    let entropy = |temperature: f64, cancelled: &mut _| {
        gas.entropy_change(initial, inside.pressure, temperature, pressure, cancelled)
            .map_err(Error::Source)
    };
    if entropy(low, cancelled)? > 0.0 || entropy(high, cancelled)? < 0.0 {
        return Err(Error::CaloricPath(
            "characteristic reference pressure leaves the source isentrope",
        ));
    }
    let mut temperature = initial;
    for _ in 0..domain.maximum_iterations {
        poll(cancelled)?;
        let cp = gas
            .heat_capacities(temperature, pressure, cancelled)
            .map_err(Error::Source)?
            .cp;
        let residual = entropy(temperature, cancelled)?;
        let tolerance = domain.absolute_energy_tolerance / temperature
            + (domain.relative_energy_tolerance + 32.0 * f64::EPSILON) * cp;
        if residual.abs() <= tolerance {
            return Ok(temperature);
        }
        if residual < 0.0 {
            low = temperature;
        } else {
            high = temperature;
        }
        let candidate = temperature - residual * (temperature / cp);
        temperature = if candidate > low && candidate < high {
            candidate
        } else {
            low + 0.5 * (high - low)
        };
    }
    Err(Error::Iterations)
}

fn characteristic_temperature<L: CaloricLaw>(
    gas: &L,
    inside: Primitive,
    target: f64,
    tolerance: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, Error<L::Error>> {
    let initial = inside.temperature;
    let domain = gas.domain();
    let (mut low, mut high) = if target < 0.0 {
        (domain.minimum_temperature, initial)
    } else {
        (initial, domain.maximum_temperature)
    };
    let edge = if target < 0.0 { low } else { high };
    let bound = integral(gas, initial, edge, 0.125 * tolerance, cancelled)?;
    if (target < 0.0 && target < bound) || (target > 0.0 && target > bound) {
        return Err(Error::CaloricPath(
            "incoming acoustic characteristic leaves the source temperature path",
        ));
    }
    let mut temperature = initial;
    for _ in 0..domain.maximum_iterations {
        poll(cancelled)?;
        let residual = integral(gas, initial, temperature, 0.125 * tolerance, cancelled)? - target;
        if residual.abs() <= 0.5 * tolerance {
            return Ok(temperature);
        }
        if residual < 0.0 {
            low = temperature;
        } else {
            high = temperature;
        }
        let capacities = gas
            .heat_capacities(temperature, inside.pressure, cancelled)
            .map_err(Error::Source)?;
        let sound_speed =
            det::sqrt((capacities.cp / capacities.cv) * gas.gas_constant() * temperature);
        let candidate = temperature - residual * sound_speed / capacities.cp;
        temperature = if candidate > low && candidate < high {
            candidate
        } else {
            low + 0.5 * (high - low)
        };
    }
    Err(Error::Iterations)
}

fn integral<L: CaloricLaw>(
    gas: &L,
    initial: f64,
    final_temperature: f64,
    tolerance: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, Error<L::Error>> {
    let solver = WallSolver::new(gas);
    if final_temperature >= initial {
        solver.rarefaction_integral(initial, final_temperature, tolerance, cancelled)
    } else {
        Ok(-solver.rarefaction_integral(final_temperature, initial, tolerance, cancelled)?)
    }
}

fn poll<E>(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error<E>> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests;
