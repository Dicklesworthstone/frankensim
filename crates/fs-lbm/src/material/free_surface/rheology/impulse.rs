//! Instantaneous momentum on the original source-rheology liquid.

use super::MaterialFreeSurfaceRheology3;
use crate::material::free_surface::{
    MaterialFreeSurfaceImpulseControls3, MaterialFreeSurfaceImpulseError,
    MaterialFreeSurfaceImpulseReport3,
};
use fs_exec::CancelGate;

impl MaterialFreeSurfaceRheology3<'_, '_> {
    /// Apply a resolved instantaneous impulse density [N s/m³] using each
    /// original cell's retained liquid mass and fixed physical units.
    ///
    /// The shared native kernel retains intrinsic density and force-corrected
    /// nonequilibrium stress. Thus the original pressure/rate source samples,
    /// collision coefficients and material identities remain in force without
    /// a second constitutive solve. The next native collision owns relaxation
    /// of the imposed velocity jump. Geometry, fill, mass carry and held forces
    /// do not change, and no collision/transport step or heat is fabricated.
    ///
    /// Actual momentum and signed kinetic work are admitted in SI units before
    /// any populations publish, with the original scenario speed ceiling and
    /// bounded cancellation. Errors preserve the complete accepted flow and
    /// material history. No full owner clone or extra population store is used.
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
