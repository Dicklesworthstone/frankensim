//! One-sided thermally perfect Euler wall Riemann problem.
//!
//! The shared borrowed-caloric kernel connects stopped gas to an incoming shock
//! or a retreating isentropic rarefaction. This adapter retains the original
//! frozen source, its error types and its existing wall/characteristic APIs.

pub(crate) mod shared;

use super::{CaloricDomain, Conserved, Error, FaceFlux, FrozenGas, Primitive, poll, positive};
use fs_qty::{Pressure, Temperature};
use shared::{CaloricLaw, HeatCapacities, WallError, WallSolver, WallState};

impl FrozenGas {
    /// Reflect one interior trace against a stationary adiabatic wall.
    /// `outward_normal` is -1 at the left face and +1 at the right face.
    pub(super) fn wall_flux(
        &self,
        inside: Conserved,
        primitive: Primitive,
        outward_normal: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<FaceFlux, Error> {
        self.moving_wall_flux(inside, primitive, outward_normal, 0.0, cancelled)
    }

    /// Flux relative to an impermeable face moving along the normal axis.
    /// The gas momentum/energy remain in the original laboratory frame.
    pub(super) fn moving_wall_flux(
        &self,
        inside: Conserved,
        primitive: Primitive,
        outward_normal: f64,
        wall_velocity: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<FaceFlux, Error> {
        poll(cancelled)?;
        if !self.admissible_state(inside)
            || !matches!(outward_normal, -1.0 | 1.0)
            || !positive(primitive.density)
            || !positive(primitive.pressure)
            || !positive(primitive.temperature)
            || !positive(primitive.sound_speed)
            || !primitive.velocity.is_finite()
        {
            return Err(Error::Invalid("wall normal or interior trace"));
        }
        let wall = if wall_velocity == 0.0 {
            self.stopped_wall(primitive, outward_normal * primitive.velocity, cancelled)?
        } else {
            WallSolver::new(self)
                .moving_wall(primitive, outward_normal, wall_velocity, cancelled)
                .map_err(wall_error)?
        };
        if !positive(wall.density) || !positive(wall.temperature) {
            return Err(Error::Invalid(
                "unrepresentable wall density or temperature",
            ));
        }
        let work = if wall_velocity == 0.0 {
            0.0
        } else {
            wall.pressure * wall_velocity
        };
        if !work.is_finite() {
            return Err(Error::Invalid("moving wall pressure work overflow"));
        }
        poll(cancelled)?;
        Ok(FaceFlux {
            mass: 0.0,
            momentum: wall.pressure,
            energy: work,
            signal_speed: wall.signal_speed,
            used_hll_fallback: false,
        })
    }

    fn stopped_wall(
        &self,
        primitive: Primitive,
        toward_wall: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<WallState, Error> {
        WallSolver::new(self)
            .stopped_wall(primitive, toward_wall, cancelled)
            .map_err(wall_error)
    }
}

impl CaloricLaw for FrozenGas {
    type Error = Error;
    fn domain(&self) -> CaloricDomain {
        FrozenGas::domain(self)
    }
    fn gas_constant(&self) -> f64 {
        FrozenGas::gas_constant(self)
    }
    fn heat_capacities(
        &self,
        t: f64,
        p: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<HeatCapacities, Error> {
        poll(cancelled)?;
        let properties = self
            .source()
            .evaluate_properties(Temperature::new(t), Pressure::new(p))?;
        Ok(HeatCapacities {
            cp: properties.mass_specific_cp().value(),
            cv: properties.mass_specific_cv().value(),
        })
    }
    fn internal_energy_change(
        &self,
        initial: f64,
        final_t: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, Error> {
        poll(cancelled)?;
        Ok(self
            .source()
            .specific_internal_energy_change(Temperature::new(initial), Temperature::new(final_t))?
            .value())
    }
    fn entropy_change(
        &self,
        initial_t: f64,
        initial_p: f64,
        final_t: f64,
        final_p: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, Error> {
        poll(cancelled)?;
        Ok(self
            .source()
            .specific_entropy_change(
                Temperature::new(initial_t),
                Pressure::new(initial_p),
                Temperature::new(final_t),
                Pressure::new(final_p),
            )?
            .value())
    }
    fn sound_speed(
        &self,
        density: f64,
        pressure: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, Error> {
        poll(cancelled)?;
        Ok(self.from_rho_u_p(density, 0.0, pressure)?.sound_speed)
    }
    fn primitive(
        &self,
        density: f64,
        velocity: f64,
        pressure: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Primitive, Error> {
        poll(cancelled)?;
        self.from_rho_u_p(density, velocity, pressure)
    }
    fn next_seam(
        &self,
        start: f64,
        mut end: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, Error> {
        for component in self.source().components() {
            poll(cancelled)?;
            for region in component.regions() {
                let seam = region.temperature_min().value();
                if seam > start && seam < end {
                    end = seam;
                }
            }
        }
        Ok(end)
    }
}

fn wall_error(error: WallError<Error>) -> Error {
    match error {
        WallError::Source(error) => error,
        WallError::Invalid(reason) => Error::Invalid(reason),
        WallError::CaloricPath(reason) => Error::CaloricPath(reason),
        WallError::Iterations => Error::TemperatureIterations,
        WallError::Cancelled => Error::Cancelled,
    }
}

#[cfg(test)]
#[path = "wall_tests.rs"]
mod tests;
