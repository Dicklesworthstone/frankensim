//! Finite-rate water phase compression in the native pressure/volume-flow chain.
//!
//! One original thermal state, authored relaxation time and frequency supply
//! both propagation and momentum impedance. The retained material owns the
//! equilibrium/frozen limits and the relaxation model's applicability; this
//! adapter adds no acoustic fit, independent density, or gas transport data.
//! The existing stable cylinder, cone, profile and compliant-wall owners do
//! all propagation. This is a local small-signal response about the retained
//! water state; the pipe does not evolve its phase inventory or heat balance.

use super::{FluidDuctError, bad, positive};
use crate::compose::{ComposeOptions, PreparedChain};
use crate::{ConeSlices, Duct, SegmentWave, WallModel};
use fs_material::phase::water::acoustics::relaxation::{
    RelaxingRegionalWater, RelaxingWaterResponse,
};

/// A source-bound water response at one positive angular frequency.
///
/// The borrowed material retains the complete original thermal state and both
/// phase-response limits. Its complex response is evaluated once here. Wave
/// and chain methods always use that same frequency, response and original
/// density; they cannot accept an unrelated material or wavenumber.
pub struct RelaxingWaterDuct<'a> {
    material: &'a RelaxingRegionalWater,
    omega_rad_s: f64,
    response: RelaxingWaterResponse,
}

impl<'a> RelaxingWaterDuct<'a> {
    /// Bind an already admitted original material to exp(-i omega t).
    ///
    /// Zero frequency is outside the native traveling-wave duct model even
    /// though the constitutive owner has a finite static compressibility.
    pub fn try_new(
        material: &'a RelaxingRegionalWater,
        omega_rad_s: f64,
    ) -> Result<Self, FluidDuctError> {
        if !positive(omega_rad_s) {
            return Err(bad("relaxing water duct requires positive finite frequency").into());
        }
        Ok(Self {
            material,
            omega_rad_s,
            response: material.response(omega_rad_s)?,
        })
    }

    /// Original material, including its thermal source and authored time scale.
    pub fn material(&self) -> &'a RelaxingRegionalWater {
        self.material
    }

    /// The angular frequency shared by the retained response and all pipe ports.
    pub fn angular_frequency_rad_s(&self) -> f64 {
        self.omega_rad_s
    }

    /// Original complex modulus, compressibility, density response and wave.
    pub fn response(&self) -> &RelaxingWaterResponse {
        &self.response
    }

    /// Original pressure/velocity impedance and its area-scaled volume-flow
    /// counterpart at one physical radius. Density is the thermal state's
    /// exact stored density. Native momentum uses z = rho * omega / k; its
    /// operation order may differ by roundoff from the constitutive diagnostic.
    pub fn segment_wave(&self, radius_m: f64) -> Result<SegmentWave, FluidDuctError> {
        super::segment_wave_with_wavenumber(
            self.material.frozen(),
            radius_m,
            self.omega_rad_s,
            self.response.wavenumber_rad_m(),
        )
    }

    /// Prepare the existing stable pressure/volume-flow chain with the same
    /// bulk phase relaxation at every radius. The frozen endpoint supplies
    /// only the coherent original density and positive real power reference;
    /// the actual complex wave is always the retained finite-rate response.
    ///
    /// Native geometry, block, conditioning, wall and cancellation admission
    /// remain in force. No wall diffusion or radiation load is inferred.
    pub fn duct_chain(
        &self,
        duct: &Duct,
        wall: WallModel<'_>,
        cone_slices: ConeSlices,
        options: ComposeOptions,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedChain, FluidDuctError> {
        super::duct_chain_with_wavenumber(
            duct,
            self.material.frozen(),
            self.omega_rad_s,
            self.response.wavenumber_rad_m(),
            wall,
            cone_slices,
            options,
            cancelled,
        )
    }
}

#[cfg(test)]
mod tests;
