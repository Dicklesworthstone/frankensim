//! Compact rigid-reservoir interfaces for the original frozen caloric source.
//!
//! The borrowed shared kernel retains one physical characteristic closure for
//! frozen and evolving mixtures. This adapter forms the original source-gauged
//! flux and preserves the public frozen reservoir API.

pub(crate) mod shared;
mod wave;

use super::{Conserved, Error, FaceFlux, FrozenGas, Primitive, poll, positive};
use crate::dynamics::wall::shared::WallError;
#[cfg(test)]
use wave::pressure_wave;

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
        let port =
            shared::reservoir_state(self, self, inside, outside, cancelled).map_err(wave_error)?;
        if port.inflow {
            self.inlet_flux(
                outside,
                port.primitive,
                outward_normal,
                port.signal_speed,
                cancelled,
            )
        } else {
            self.outlet_flux(port.primitive, outward_normal, port.signal_speed)
        }
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

    fn inlet_flux(
        &self,
        reservoir: Primitive,
        port: Primitive,
        outward_normal: f64,
        signal_speed: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<FaceFlux, Error> {
        poll(cancelled)?;
        // Retain the original sensible stagnation enthalpy independently of the
        // formation reference and the kinetic head; supply it exactly once.
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

fn wave_error(error: WallError<Error>) -> Error {
    match error {
        WallError::Source(error) => error,
        WallError::Invalid(reason) => Error::Invalid(reason),
        WallError::CaloricPath(reason) => Error::CaloricPath(reason),
        WallError::Iterations => Error::TemperatureIterations,
        WallError::Cancelled => Error::Cancelled,
    }
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
