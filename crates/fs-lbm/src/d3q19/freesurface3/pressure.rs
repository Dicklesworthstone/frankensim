//! Prescribed external pressure on the original gas reconstruction.
//!
//! The liquid lattice represents pressure differences about rho=1. Its
//! reference density is not the density of the external gas. The normal load
//! is p_liquid - p_reference = p_gas - p_reference + sigma*kappa.

use super::*;

pub(super) struct GasPressure3<'a> {
    field: &'a dyn Fn([usize; 3]) -> Result<f64, FreeSurfaceError3>,
    maximum_density_excursion: f64,
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
    pub fn step_with_gas_pressure(
        &mut self,
        gate: &CancelGate,
        maximum_density_excursion: f64,
        pressure: impl Fn([usize; 3]) -> Result<f64, FreeSurfaceError3>,
    ) -> Result<(), FreeSurfaceError3> {
        check_cancel(gate)?;
        if !maximum_density_excursion.is_finite()
            || maximum_density_excursion <= 0.0
            || maximum_density_excursion >= 1.0
        {
            return Err(FreeSurfaceError3::InvalidGasPressure {
                coordinates: None,
                reason: "density-excursion bound must be finite and in (0,1)",
            });
        }
        self.step_inner(
            gate,
            None,
            None,
            Some(&GasPressure3 {
                field: &pressure,
                maximum_density_excursion,
            }),
        )
    }

    pub(super) fn loaded_reference_density(
        &self,
        key: u64,
        lane: usize,
        gas_pressure: Option<&GasPressure3<'_>>,
    ) -> Result<f64, FreeSurfaceError3> {
        let capillary_density = self.reference_density(key, lane);
        let Some(boundary) = gas_pressure else {
            return Ok(capillary_density);
        };
        let (x, y, z) = Self::global_of(key, lane);
        let coordinates = [x as usize, y as usize, z as usize];
        let pressure = (boundary.field)(coordinates)?;
        let density = capillary_density + pressure / CS2_3;
        if !pressure.is_finite()
            || !density.is_finite()
            || density <= 0.0
            || (density - 1.0).abs() > boundary.maximum_density_excursion
        {
            return Err(FreeSurfaceError3::InvalidGasPressure {
                coordinates: Some(coordinates),
                reason: "gas plus capillary pressure exceeds the declared density-excursion bound",
            });
        }
        Ok(density)
    }
}
