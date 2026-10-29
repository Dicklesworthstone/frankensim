//! Species-bearing compact-reservoir flux with one shared characteristic solve.

use super::wall::LocalCaloric;
use super::{Cell, Error, FaceFlux, MAX_COMPONENTS, MixtureGas, Primitive, poll, positive};
use crate::dynamics::characteristic::shared::reservoir_state;
use crate::dynamics::wall::shared::WallError;

impl MixtureGas {
    /// Flux between a duct trace and a compact stationary gas reservoir.
    ///
    /// Both states use this exact species bank, canonical species order and
    /// per-species energy references, but their compositions may differ. Inflow
    /// carries reservoir species and stagnation enthalpy; outflow carries the
    /// interior species and its resolved total enthalpy. Each side retains its
    /// original temperature-dependent caloric law. Sonic inflow/outflow and
    /// departing supersonic characteristics use the shared physical closure.
    ///
    /// The normal is -1 at a left end, +1 at a right end. Output is oriented
    /// along increasing duct coordinate and has species_count()+2 entries:
    /// all species mass fluxes, momentum flux, then referenced total-energy
    /// flux. Output remains unchanged on refusal or cancellation. Successful
    /// evaluation allocates no heap workspace or replacement gas owner.
    pub fn reservoir_flux(
        &self,
        interior: Cell<'_>,
        reservoir: Cell<'_>,
        outward_normal: f64,
        out: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FaceFlux, Error> {
        let inside = self.primitive(interior, &mut cancelled)?;
        let outside = self.primitive(reservoir, &mut cancelled)?;
        self.plenum_flux(
            interior,
            inside,
            reservoir,
            outside,
            outward_normal,
            out,
            &mut cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn plenum_flux(
        &self,
        interior: Cell<'_>,
        inside: Primitive,
        reservoir: Cell<'_>,
        outside: Primitive,
        outward_normal: f64,
        out: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<FaceFlux, Error> {
        poll(cancelled)?;
        let n = self.species_count();
        if out.len() != n + 2
            || !matches!(outward_normal, -1.0 | 1.0)
            || reservoir.momentum != 0.0
            || !self.admissible_cell(interior)
            || !self.admissible_cell(reservoir)
        {
            return Err(Error::Invalid(
                "mixture reservoir dimension, normal or stationary source state",
            ));
        }
        let duct_law = LocalCaloric::new(self, interior.partial_densities)?;
        let reservoir_law = LocalCaloric::new(self, reservoir.partial_densities)?;
        let inside = Primitive {
            velocity: outward_normal * inside.velocity,
            ..inside
        };
        let port = reservoir_state(&duct_law, &reservoir_law, inside, outside, cancelled)
            .map_err(characteristic_error)?;
        let donor = if port.inflow {
            &reservoir_law
        } else {
            &duct_law
        };
        let state = port.primitive;
        let mut candidate = [0.0; MAX_COMPONENTS];
        let admitted = donor.state_at(
            state.density,
            state.velocity,
            state.pressure,
            &mut candidate[..n + 2],
        )?;
        let velocity = outward_normal * state.velocity;
        let mass = admitted.density * velocity;
        let momentum = candidate[n] * state.velocity + state.pressure;
        let energy = if port.inflow {
            // The stored reservoir energy already subtracts the original
            // per-species references. Its stationary total enthalpy is supplied
            // once, even when the donor composition differs from the duct.
            mass * (reservoir.energy / outside.density + outside.pressure / outside.density)
        } else {
            (candidate[n + 1] + state.pressure) * velocity
        };
        for flux in &mut candidate[..n] {
            *flux *= velocity;
        }
        candidate[n] = momentum;
        candidate[n + 1] = energy;
        if !candidate[..n + 2].iter().all(|x| x.is_finite())
            || !positive(momentum)
            || !positive(port.signal_speed)
        {
            return Err(Error::Invalid(
                "unrepresentable species-bearing reservoir flux",
            ));
        }
        poll(cancelled)?;
        out.copy_from_slice(&candidate[..n + 2]);
        Ok(FaceFlux {
            signal_speed: port.signal_speed,
            used_hll_fallback: false,
        })
    }
}

fn characteristic_error(error: WallError<Error>) -> Error {
    match error {
        WallError::Source(error) => error,
        WallError::Invalid(reason) => Error::Invalid(reason),
        WallError::CaloricPath(reason) => Error::CharacteristicPath(reason),
        WallError::Iterations => Error::TemperatureIterations,
        WallError::Cancelled => Error::Cancelled,
    }
}
