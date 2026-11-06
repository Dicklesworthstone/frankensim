//! Compact rigid-reservoir interfaces for the original frozen caloric source.
//!
//! The borrowed shared kernel retains one physical characteristic closure for
//! frozen and evolving mixtures. This adapter forms the original source-gauged
//! flux and preserves the public frozen reservoir API.

mod acoustic;
pub(crate) mod shared;
mod wave;

use super::{Conserved, Error, FaceFlux, FrozenGas, Primitive, poll, positive};
use crate::dynamics::wall::shared::WallError;
#[cfg(test)]
use wave::pressure_wave;

impl FrozenGas {
    /// Acoustic characteristic outlet with reference absolute pressure [Pa]
    /// and signed axial velocity [m/s].
    ///
    /// For subsonic outflow, only the incoming acoustic characteristic is
    /// prescribed. The outgoing characteristic and entropy come from the
    /// interior. Pressure is free to fluctuate: a smooth outgoing isentropic
    /// simple wave and an outgoing constant-pressure entropy contact pass
    /// without reflection by this local closure. The original Cp(T) supplies
    /// the characteristic integral; no constant-gamma replacement is used.
    ///
    /// The reference pressure is continued along the interior isentrope, so
    /// it does not impose a remote temperature or erase an entropy contact.
    /// Both velocities use the increasing duct coordinate; `outward_normal`
    /// is -1 at a left end and +1 at a right end. Supersonic outflow uses the
    /// interior flux and an expansion crossing sonic speed is sampled at its
    /// sonic state. Backflow and unsupported caloric paths refuse. This is a
    /// one-dimensional smooth-wave condition, not a nonreflection guarantee
    /// for shocks or interactions with spatial entropy gradients.
    pub fn characteristic_outlet_flux(
        &self,
        interior: Conserved,
        outward_normal: f64,
        reference_pressure: f64,
        reference_velocity: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FaceFlux, Error> {
        let inside = self.decode(interior, None, &mut cancelled)?;
        let (face, signal) = self.characteristic_outlet_state(
            inside,
            outward_normal,
            reference_pressure,
            reference_velocity,
            &mut cancelled,
        )?;
        self.face_flux(face, signal)
    }

    pub(super) fn characteristic_outlet_state(
        &self,
        inside: Primitive,
        outward_normal: f64,
        reference_pressure: f64,
        reference_velocity: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Primitive, f64), Error> {
        poll(cancelled)?;
        if !matches!(outward_normal, -1.0 | 1.0) {
            return Err(Error::Invalid("characteristic outlet outward normal"));
        }
        let (face, signal) = acoustic::outlet(
            self,
            Primitive {
                velocity: outward_normal * inside.velocity,
                ..inside
            },
            reference_pressure,
            outward_normal * reference_velocity,
            cancelled,
        )
        .map_err(wave_error)?;
        Ok((
            Primitive {
                velocity: outward_normal * face.velocity,
                ..face
            },
            signal,
        ))
    }

    /// Outgoing Euler flux selected by a positive absolute backpressure [Pa].
    ///
    /// The original caloric source supplies compression shocks and isentropic
    /// expansions. A fan crossing the boundary is sampled at its sonic state;
    /// a departing supersonic wave leaves the interior trace unchanged. The
    /// normal is -1 or +1; flux stays oriented along the increasing coordinate.
    /// Backflow is refused because pressure alone cannot supply entering gas
    /// temperature or entropy. A fixed pressure is not a nonreflecting outlet.
    pub fn pressure_outlet_flux(
        &self,
        interior: Conserved,
        outward_normal: f64,
        backpressure: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FaceFlux, Error> {
        let inside = self.decode(interior, None, &mut cancelled)?;
        let (face, signal) =
            self.pressure_outlet_state(inside, outward_normal, backpressure, &mut cancelled)?;
        self.face_flux(face, signal)
    }

    pub(super) fn pressure_outlet_state(
        &self,
        inside: Primitive,
        outward_normal: f64,
        backpressure: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Primitive, f64), Error> {
        poll(cancelled)?;
        if !matches!(outward_normal, -1.0 | 1.0) {
            return Err(Error::Invalid("pressure outlet outward normal"));
        }
        let port = shared::pressure_outlet_state(
            self,
            Primitive {
                velocity: outward_normal * inside.velocity,
                ..inside
            },
            backpressure,
            cancelled,
        )
        .map_err(wave_error)?;
        poll(cancelled)?;
        Ok((
            Primitive {
                velocity: outward_normal * port.primitive.velocity,
                ..port.primitive
            },
            port.signal_speed,
        ))
    }

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
    /// end and +1 at a right end. No gamma-law replacement is used. The supplied
    /// temperature and sound speed must agree with this original gas source;
    /// stale or caller-invented auxiliary thermodynamics refuse.
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
        require_matching_thermodynamics(
            interior,
            self.from_rho_u_p(interior.density, interior.velocity, interior.pressure)?,
        )?;
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
    /// and whether the throat is sonic (choked). Even at equal pressures the
    /// donor must carry temperature and sound speed from this gas source.
    pub fn nozzle_exit(
        &self,
        donor: Primitive,
        backpressure: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(Primitive, bool), Error> {
        poll(&mut cancelled)?;
        require_matching_thermodynamics(
            donor,
            self.from_rho_u_p(donor.density, donor.velocity, donor.pressure)?,
        )?;
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

    /// Static state reached from stagnation temperature [K] at caller-owned
    /// pressure [Pa] and signed axial velocity [m/s].
    ///
    /// The original source solves `h(T0) - h(T) = velocity^2 / 2` with its
    /// temperature-dependent heat capacity and bounded caloric inversion. No
    /// state is extrapolated beyond the declared temperature interval. This
    /// closes the energy of an adiabatic, work-free inlet, including a jet that
    /// mixes to a slower face velocity. The caller supplies the pressure/drop
    /// law and any choking limit; this relation does not impose an isentrope.
    pub fn stagnation_to_static_state(
        &self,
        stagnation_temperature: f64,
        pressure: f64,
        velocity: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Primitive, Error> {
        let temperature = shared::static_temperature(
            self,
            stagnation_temperature,
            pressure,
            velocity.abs(),
            &mut cancelled,
        )
        .map_err(wave_error)?;
        poll(&mut cancelled)?;
        let density = (pressure / self.gas_constant()) / temperature;
        let face = self.state_at(density, velocity, pressure)?.1;
        poll(&mut cancelled)?;
        Ok(face)
    }

    /// Euler flux through a stationary face that carries `face`, oriented
    /// along the increasing duct coordinate. `signal_speed` is raised to the
    /// face's own `|u| + c` when smaller. The face's thermodynamics must agree
    /// with this source, including the sound speed used for the CFL bound.
    pub fn face_flux(&self, face: Primitive, signal_speed: f64) -> Result<FaceFlux, Error> {
        let (q, admitted) = self.state_at(face.density, face.velocity, face.pressure)?;
        require_matching_thermodynamics(face, admitted)?;
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
        require_matching_thermodynamics(
            face,
            self.from_rho_u_p(face.density, face.velocity, face.pressure)?,
        )?;
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

// Public Primitive fields can be authored independently. Only the original
// source may supply their redundant thermodynamic coordinates. Allow arithmetic
// round trips (including conservative decoding), not an energy-solve tolerance
// that could admit a materially different acoustic speed. Keep the caller's
// valid bits, so a zero-strength wave remains the exact original state.
fn require_matching_thermodynamics(supplied: Primitive, source: Primitive) -> Result<(), Error> {
    for (value, expected) in [
        (supplied.pressure, source.pressure),
        (supplied.temperature, source.temperature),
        (supplied.sound_speed, source.sound_speed),
    ] {
        if !positive(value) || ((value - expected) / expected).abs() > 64.0 * f64::EPSILON {
            return Err(Error::Invalid(
                "primitive thermodynamics disagree with the original gas source",
            ));
        }
    }
    Ok(())
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
