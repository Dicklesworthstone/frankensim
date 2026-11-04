//! Mechanical impulse on the original thermal liquid trajectory.

use super::MaterialFreeSurfaceEnthalpy3;
use crate::material::free_surface::{
    MaterialFreeSurfaceImpulseControls3, MaterialFreeSurfaceImpulseError,
    MaterialFreeSurfaceImpulseReport3,
};
use fs_exec::CancelGate;

impl MaterialFreeSurfaceEnthalpy3<'_, '_> {
    /// Apply a resolved instantaneous impulse density [N s/m³] to the retained
    /// liquid, using the whole geometric cell volume and original SI scaling.
    ///
    /// The native liquid mass determines each velocity jump. Its canonical
    /// density and force-corrected viscous stress are retained, so pressure,
    /// caloric/source states, phase fractions, rheology samples, wetting history
    /// and held thermal/buoyancy forces remain bound to the same trajectory.
    /// No mass route, accepted step or physical time advances.
    ///
    /// The returned signed mechanical work changes liquid kinetic energy. It
    /// is not deposited as enthalpy or substituted for molecular dissipation.
    /// Subsequent native transport and the selected heating owners account for
    /// evolution of the resulting motion. This is not a pressure projection,
    /// compressible shock or paired solid/liquid collision law.
    ///
    /// Native speed/balance, physical-unit receipt and cancellation admission
    /// all precede publication. A refusal leaves both flow and thermal state
    /// intact without cloning the thermal or material owner.
    pub fn apply_impulse_density(
        &mut self,
        impulse_n_s_m3: &[[f64; 3]],
        controls: MaterialFreeSurfaceImpulseControls3,
        gate: &CancelGate,
    ) -> Result<MaterialFreeSurfaceImpulseReport3, MaterialFreeSurfaceImpulseError> {
        self.flow
            .apply_impulse_density(impulse_n_s_m3, controls, gate)
    }
}
