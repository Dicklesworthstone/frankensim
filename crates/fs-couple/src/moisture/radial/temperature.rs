//! Prescribed uniform temperature events on the original curved inventory.

use super::*;

impl MoistureRadial {
    /// Current prescribed absolute temperature, distinct from the immutable
    /// density and stress-free datum held by [`Self::reference`].
    pub fn temperature_kelvin(&self) -> f64 {
        self.field().temperature_kelvin()
    }

    /// Re-query the original pinned diffusivity at every occupied shell content
    /// after an externally prescribed uniform-temperature event. The exact
    /// multiaxis source support is checked; a rectangular validity envelope is
    /// not permission to interpolate outside the original source triangulation.
    ///
    /// Water, each shell's dry mass, geometry, accepted time, source pins and
    /// reference identity survive unchanged. Degenerate D(0)=0 remains admitted
    /// by the native Kirchhoff path. Any unsupported cell or cancellation
    /// refuses the entire event before a replacement is returned.
    ///
    /// Cut diffusion calls at temperature events and provide boundary/reservoir
    /// laws qualified at the new temperature. This does not supply heat balance,
    /// thermal expansion or temperature evolution. Native hygroscopic mechanics
    /// retains its refusal of an unqualified change in reference temperature.
    pub fn at_prescribed_temperature(
        &self,
        temperature_kelvin: f64,
        cx: &Cx<'_>,
    ) -> Result<Self, MoistureError> {
        let transport = binding::at_prescribed_temperature(
            &self.card,
            &self.transport,
            temperature_kelvin,
            &self.moisture,
            true,
            cx,
        )?;
        let mut next = self.clone();
        next.transport = transport;
        poll(cx)?;
        Ok(next)
    }
}
