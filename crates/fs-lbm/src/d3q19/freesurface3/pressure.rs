//! Prescribed external pressure on the original gas reconstruction.
//!
//! The liquid lattice represents pressure differences about rho=1. Its
//! reference density is not the density of the external gas. The normal load
//! is p_liquid - p_reference = p_gas - p_reference + sigma*kappa.

use super::*;

#[cfg(test)]
mod tests;

/// One actual missing incoming population at a liquid/gas link.
///
/// The gas cell is the pull source `liquid_cell - E3[incoming_direction]`.
/// It is an active native gas cell; inactive and out-of-domain sources remain
/// walls and are never passed to a pressure callback. Coordinates describe
/// lattice cell centers, not a reconstructed interface intersection or normal.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct GasPressureLink3 {
    /// Liquid collision cell receiving the reconstructed population.
    pub liquid_cell: [usize; 3],
    /// Signed coordinates of the actual neighboring gas source.
    pub gas_cell: [i64; 3],
    /// Incoming native D3Q19 population index, in `1..Q3`.
    pub incoming_direction: usize,
}

#[derive(Clone, Copy)]
pub(crate) enum GasPressureField3<'a> {
    Cell(&'a dyn Fn([usize; 3]) -> Result<f64, FreeSurfaceError3>),
    Link(&'a dyn Fn(GasPressureLink3) -> Result<f64, FreeSurfaceError3>),
}

pub(super) struct GasPressure3<'a> {
    field: GasPressureField3<'a>,
    maximum_density_excursion: f64,
}

impl<'a> GasPressure3<'a> {
    pub(super) fn new(
        maximum_density_excursion: f64,
        field: GasPressureField3<'a>,
    ) -> Result<Self, FreeSurfaceError3> {
        if !maximum_density_excursion.is_finite()
            || maximum_density_excursion <= 0.0
            || maximum_density_excursion >= 1.0
        {
            return Err(FreeSurfaceError3::InvalidGasPressure {
                coordinates: None,
                reason: "density-excursion bound must be finite and in (0,1)",
            });
        }
        Ok(Self {
            field,
            maximum_density_excursion,
        })
    }

    fn density(
        &self,
        capillary_density: f64,
        link: GasPressureLink3,
    ) -> Result<f64, FreeSurfaceError3> {
        let pressure = match self.field {
            GasPressureField3::Cell(field) => field(link.liquid_cell)?,
            GasPressureField3::Link(field) => field(link)?,
        };
        let density = capillary_density + pressure / CS2_3;
        if !pressure.is_finite()
            || !density.is_finite()
            || density <= 0.0
            || (density - 1.0).abs() > self.maximum_density_excursion
        {
            return Err(FreeSurfaceError3::InvalidGasPressure {
                coordinates: Some(link.liquid_cell),
                reason: "gas plus capillary pressure exceeds the declared density-excursion bound",
            });
        }
        Ok(density)
    }
}

// Curvature and velocity are shared by a cell's links. A cell field is sampled
// exactly once; a link field is sampled at every actual gas link. Reusing the
// last density keeps uniform link loads bit-identical to cell loads and avoids
// repeated equilibrium evaluation without a map, allocation or second pass.
#[derive(Default)]
pub(super) struct GasEquilibrium3 {
    capillary_density: Option<f64>,
    loaded_density: Option<f64>,
    equilibrium: Option<[f64; Q3]>,
}

impl GasEquilibrium3 {
    pub(super) fn at_link(
        &mut self,
        reference_density: impl FnOnce() -> f64,
        velocity: [f64; 3],
        boundary: Option<&GasPressure3<'_>>,
        link: GasPressureLink3,
    ) -> Result<&[f64; Q3], FreeSurfaceError3> {
        let per_link =
            boundary.is_some_and(|load| matches!(load.field, GasPressureField3::Link(_)));
        if self.loaded_density.is_none() || per_link {
            let capillary = *self.capillary_density.get_or_insert_with(reference_density);
            let density = match boundary {
                Some(load) => load.density(capillary, link)?,
                None => capillary,
            };
            if self
                .loaded_density
                .is_none_or(|previous| previous.to_bits() != density.to_bits())
            {
                self.equilibrium = Some(equilibrium3(density, velocity));
                self.loaded_density = Some(density);
            }
        }
        Ok(self
            .equilibrium
            .as_ref()
            .expect("reconstructed gas boundary"))
    }
}

impl FreeSurface3 {
    /// Advance with an externally prescribed gas-pressure field in lattice
    /// pressure units, relative to the same reference as liquid density one.
    ///
    /// The callback is sampled once per exposed liquid interface cell, in the
    /// native deterministic order, at its integer cell coordinates. It must
    /// describe one frozen time-level boundary field for this step. Returning
    /// zero reproduces the ordinary atmospheric reconstruction exactly.
    ///
    /// The imposed load adds to the existing sigma*kappa traction. The combined
    /// reconstruction must satisfy |rho_boundary - 1| <= the caller's finite
    /// `maximum_density_excursion` in (0,1). This numerical bound is necessary,
    /// not an accuracy certificate; lattice compressibility is not a liquid or
    /// gas equation of state. Gas inertia, gas evolution and external pressure
    /// work are not computed by this single-phase isothermal owner.
    ///
    /// The field is supplied on EVERY step; it is never retained as a stale
    /// interface map after liquid advances or new sparse tiles activate.
    /// Gas-pressure refusal or observed cancellation precedes population/mass
    /// publication; existing post-publication conversion error semantics remain.
    /// With either conservative momentum mode, the complete pressure/flow/
    /// conversion candidate publishes together, and any refusal preserves the
    /// prior native state.
    pub fn step_with_gas_pressure(
        &mut self,
        gate: &CancelGate,
        maximum_density_excursion: f64,
        pressure: impl Fn([usize; 3]) -> Result<f64, FreeSurfaceError3>,
    ) -> Result<(), FreeSurfaceError3> {
        check_cancel(gate)?;
        let boundary = GasPressure3::new(
            maximum_density_excursion,
            GasPressureField3::Cell(&pressure),
        )?;
        self.step_with_pressure_boundary(gate, &boundary)
    }

    /// Advance with independently prescribed gas pressure on each exposed link.
    ///
    /// Each actual missing incoming gas population receives its own pressure
    /// perturbation in lattice units about the density-one reference. The
    /// callback sees the exact liquid and neighboring gas cell and is sampled
    /// in native cell order, then increasing incoming D3Q19 direction. It must
    /// describe one frozen boundary field; distinct neighbors of the SAME
    /// liquid cell may supply distinct loads. Uniform link pressure reproduces
    /// [`Self::step_with_gas_pressure`] population arithmetic exactly.
    ///
    /// Every link retains the original capillary jump, physical Guo source
    /// frame and combined density-excursion bound. Link directions are lattice
    /// characteristics, not reconstructed geometric surface normals. This is
    /// prescribed normal loading without gas topology, gas dynamics, physical
    /// gas volume or pressure-work evolution. The same atomicity and native
    /// conversion error semantics as [`Self::step_with_gas_pressure`] apply.
    pub fn step_with_gas_link_pressure(
        &mut self,
        gate: &CancelGate,
        maximum_density_excursion: f64,
        pressure: impl Fn(GasPressureLink3) -> Result<f64, FreeSurfaceError3>,
    ) -> Result<(), FreeSurfaceError3> {
        check_cancel(gate)?;
        let boundary = GasPressure3::new(
            maximum_density_excursion,
            GasPressureField3::Link(&pressure),
        )?;
        self.step_with_pressure_boundary(gate, &boundary)
    }

    fn step_with_pressure_boundary(
        &mut self,
        gate: &CancelGate,
        boundary: &GasPressure3<'_>,
    ) -> Result<(), FreeSurfaceError3> {
        if self.redistribution.preserves_momentum() || self.body_acceleration.is_some() {
            let mut candidate = self.clone();
            candidate.step_inner(gate, None, None, Some(boundary))?;
            check_cancel(gate)?;
            *self = candidate;
            return Ok(());
        }
        self.step_inner(gate, None, None, Some(boundary))
    }
}
