//! Species-bearing acoustic/pressure outlets and compact-reservoir fluxes.

mod acoustic;
mod pressure;

use super::wall::LocalCaloric;
use super::{
    Cell, Error, FaceFlux, MAX_NORMAL_COMPONENTS, MAX_SPATIAL_COMPONENTS, MixtureGas, NormalCell,
    NormalPrimitive, Primitive, SpatialCell, SpatialPrimitive, poll, positive,
};
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

    /// Planar compact-reservoir port with original donor species and energy.
    ///
    /// Reservoir momentum must be exactly zero. Its pressure and temperature
    /// supply inflow stagnation conditions; incoming gas has zero shear.
    /// Outflow retains interior composition and transverse kinetic energy.
    /// Output is species, normal momentum, transverse momentum and full energy.
    /// The one-sided choking/reversal law is that of [`Self::reservoir_flux`].
    pub fn reservoir_normal_flux(
        &self,
        interior: NormalCell<'_>,
        reservoir: Cell<'_>,
        outward_normal: f64,
        out: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FaceFlux, Error> {
        let inside = self.normal_primitive(interior, &mut cancelled)?;
        let outside = self.primitive(reservoir, &mut cancelled)?;
        self.reservoir_components(
            interior.partial_densities,
            inside,
            [interior.transverse_momentum / inside.density],
            reservoir,
            outside,
            outward_normal,
            out,
            &mut cancelled,
        )
    }

    /// Source-admitted reconstruction path to [`Self::reservoir_normal_flux`].
    /// A complete source refusal or cancellation leaves output unchanged.
    pub fn reservoir_normal_flux_from_primitive(
        &self,
        interior: NormalPrimitive<'_>,
        reservoir: Cell<'_>,
        outward_normal: f64,
        out: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FaceFlux, Error> {
        poll(&mut cancelled)?;
        let mut state = [0.0; MAX_NORMAL_COMPONENTS];
        let inside =
            self.normal_conserved_into(interior, &mut state[..self.species_count() + 3])?;
        let outside = self.primitive(reservoir, &mut cancelled)?;
        self.reservoir_components(
            interior.partial_densities,
            inside,
            [interior.transverse_velocity],
            reservoir,
            outside,
            outward_normal,
            out,
            &mut cancelled,
        )
    }

    /// Spatial compact-reservoir port retaining both departing shear velocities.
    ///
    /// Both states use the original canonical species bank and energy gauges.
    /// Incoming gas takes the stationary reservoir composition and stagnation
    /// enthalpy, with exactly zero transverse flux. Departing gas transports
    /// both shear momenta and complete kinetic energy. No finite storage,
    /// exterior velocity trace or empirical discharge coefficient is inferred.
    pub fn reservoir_spatial_flux(
        &self,
        interior: SpatialCell<'_>,
        reservoir: Cell<'_>,
        outward_normal: f64,
        out: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FaceFlux, Error> {
        let inside = self.spatial_primitive(interior, &mut cancelled)?;
        let outside = self.primitive(reservoir, &mut cancelled)?;
        self.reservoir_components(
            interior.partial_densities,
            inside,
            interior.transverse_momentum.map(|m| m / inside.density),
            reservoir,
            outside,
            outward_normal,
            out,
            &mut cancelled,
        )
    }

    /// Source-admitted reconstruction path to [`Self::reservoir_spatial_flux`].
    /// Both shear precision checks and atomic publication remain binding.
    pub fn reservoir_spatial_flux_from_primitive(
        &self,
        interior: SpatialPrimitive<'_>,
        reservoir: Cell<'_>,
        outward_normal: f64,
        out: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FaceFlux, Error> {
        poll(&mut cancelled)?;
        let mut state = [0.0; MAX_SPATIAL_COMPONENTS];
        let inside =
            self.spatial_conserved_into(interior, &mut state[..self.species_count() + 4])?;
        let outside = self.primitive(reservoir, &mut cancelled)?;
        self.reservoir_components(
            interior.partial_densities,
            inside,
            interior.transverse_velocity,
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
        self.reservoir_components(
            interior.partial_densities,
            inside,
            [],
            reservoir,
            outside,
            outward_normal,
            out,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn reservoir_components<const TANGENTS: usize>(
        &self,
        partial_densities: &[f64],
        inside: Primitive,
        transverse_velocity: [f64; TANGENTS],
        reservoir: Cell<'_>,
        outside: Primitive,
        outward_normal: f64,
        out: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<FaceFlux, Error> {
        poll(cancelled)?;
        let n = self.species_count();
        let energy_index = n + TANGENTS + 1;
        if TANGENTS > 2
            || out.len() != energy_index + 1
            || !matches!(outward_normal, -1.0 | 1.0)
            || reservoir.momentum != 0.0
            || !self.admissible_cell(reservoir)
            || !transverse_velocity.iter().all(|v| v.is_finite())
        {
            return Err(Error::Invalid(
                "compact reservoir dimension, normal, stationary state or shear",
            ));
        }
        let duct_law = LocalCaloric::new(self, partial_densities)?;
        let reservoir_law = LocalCaloric::new(self, reservoir.partial_densities)?;
        let inside = Primitive {
            velocity: outward_normal * inside.velocity,
            ..inside
        };
        let port = reservoir_state(&duct_law, &reservoir_law, inside, outside, cancelled)
            .map_err(characteristic_error)?;
        if !port.inflow {
            return self.outgoing_characteristic_flux(
                &duct_law,
                port.primitive,
                transverse_velocity,
                outward_normal,
                port.signal_speed,
                out,
                cancelled,
            );
        }
        let state = port.primitive;
        let mut candidate = [0.0; MAX_SPATIAL_COMPONENTS];
        let admitted = reservoir_law.state_at(
            state.density,
            state.velocity,
            state.pressure,
            &mut candidate[..n + 2],
        )?;
        let velocity = outward_normal * state.velocity;
        let mass = admitted.density * velocity;
        let momentum = candidate[n] * state.velocity + state.pressure;
        // Incoming gas accelerates from rest. Its source-gauged stagnation
        // enthalpy supplies pressure work and normal kinetic energy once;
        // interior transverse velocity cannot create extra reservoir energy.
        let energy =
            mass * (reservoir.energy / outside.density + outside.pressure / outside.density);
        for value in &mut candidate[..n] {
            *value *= velocity;
        }
        candidate[n] = momentum;
        candidate[n + 1..energy_index].fill(0.0);
        candidate[energy_index] = energy;
        let signal_speed = port.signal_speed.max(velocity.abs() + admitted.sound_speed);
        if !candidate[..=energy_index].iter().all(|x| x.is_finite())
            || !positive(momentum)
            || !positive(signal_speed)
        {
            return Err(Error::Invalid(
                "unrepresentable species-bearing reservoir flux",
            ));
        }
        poll(cancelled)?;
        out.copy_from_slice(&candidate[..=energy_index]);
        Ok(FaceFlux {
            signal_speed,
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
