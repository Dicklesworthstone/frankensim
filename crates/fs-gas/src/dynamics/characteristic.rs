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

/// The interior side of one classical pressure wave at a duct end.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EndWave {
    /// State between the wave and the contact. Its velocity is along the
    /// increasing duct coordinate, like every other primitive.
    pub primitive: Primitive,
    /// Derivative of that velocity with respect to the end pressure [m/(s Pa)].
    pub velocity_per_pressure: f64,
    /// Largest characteristic or shock speed the wave involves [m/s].
    pub signal_speed: f64,
}

impl FrozenGas {
    /// The state an interior trace reaches when a duct end imposes
    /// `end_pressure`: this source's shock Hugoniot for compression, its
    /// isentropic fan for expansion. This is the interior's half of an end
    /// Riemann problem; a caller-owned port supplies the other relation
    /// between end pressure and end velocity. `outward_normal` is -1 at a left
    /// end and +1 at a right end. No gamma-law replacement is used.
    pub fn end_wave(
        &self,
        interior: Primitive,
        outward_normal: f64,
        end_pressure: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<EndWave, Error> {
        poll(&mut cancelled)?;
        if !matches!(outward_normal, -1.0 | 1.0) {
            return Err(Error::Invalid("end wave outward normal"));
        }
        let inside = Primitive {
            velocity: outward_normal * interior.velocity,
            ..interior
        };
        let wave =
            wave::pressure_wave(self, inside, end_pressure, &mut cancelled).map_err(wave_error)?;
        Ok(EndWave {
            primitive: Primitive {
                velocity: outward_normal * wave.primitive.velocity,
                ..wave.primitive
            },
            velocity_per_pressure: -outward_normal * wave.derivative,
            signal_speed: wave.signal_speed,
        })
    }

    /// Isentropic discharge of stationary stagnation gas `donor` through a
    /// converging throat into `backpressure`, along this source's isentrope.
    /// Returns the throat state, whose velocity is the positive exit speed,
    /// and whether the throat is sonic (choked).
    pub fn nozzle_exit(
        &self,
        donor: Primitive,
        backpressure: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(Primitive, bool), Error> {
        shared::nozzle::state(self, donor, backpressure, &mut cancelled).map_err(wave_error)
    }

    /// Referenced specific stagnation enthalpy [J/kg] of gas at rest at
    /// `temperature`: the energy each kilogram entering from rest carries.
    pub fn stagnation_enthalpy(&self, temperature: f64) -> Result<f64, Error> {
        let enthalpy = self.sensible(temperature)? + self.gas_constant() * temperature;
        if !enthalpy.is_finite() {
            return Err(Error::Invalid("unrepresentable stagnation enthalpy"));
        }
        Ok(enthalpy)
    }

    /// Euler flux through a stationary face that carries `face`, oriented
    /// along the increasing duct coordinate. `signal_speed` is raised to the
    /// face's own `|u| + c` when smaller.
    pub fn face_flux(&self, face: Primitive, signal_speed: f64) -> Result<FaceFlux, Error> {
        let q = self.conserved(face.density, face.velocity, face.pressure)?;
        checked_flux(
            q.momentum,
            q.momentum * face.velocity + face.pressure,
            (q.energy + face.pressure) * face.velocity,
            signal_speed.max(face.velocity.abs() + face.sound_speed),
        )
    }

    /// Flux of gas that entered from rest at `stagnation_temperature` and
    /// crosses a stationary face with density, velocity and pressure `face`.
    /// Energy is the mass flux times the stagnation enthalpy, exactly once.
    pub fn inflow_flux(
        &self,
        face: Primitive,
        stagnation_temperature: f64,
        signal_speed: f64,
    ) -> Result<FaceFlux, Error> {
        let enthalpy = self.stagnation_enthalpy(stagnation_temperature)?;
        let mass = face.density * face.velocity;
        checked_flux(
            mass,
            mass * face.velocity + face.pressure,
            mass * enthalpy,
            signal_speed.max(face.velocity.abs() + face.sound_speed),
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
