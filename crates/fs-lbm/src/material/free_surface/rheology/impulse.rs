//! Instantaneous momentum on the original source-rheology liquid.

use super::MaterialFreeSurfaceRheology3;
use crate::material::free_surface::{
    MaterialFreeSurfaceImpulseControls3, MaterialFreeSurfaceImpulseError,
    MaterialFreeSurfaceMomentum3, MaterialFreeSurfaceScaling,
    MaterialFreeSurfaceImpulseReport3,
    MaterialFreeSurfaceWrenchImpulse3, MaterialFreeSurfaceWrenchImpulseControls3,
    MaterialFreeSurfaceWrenchImpulseReport3,
};
use fs_exec::CancelGate;

impl MaterialFreeSurfaceRheology3<'_, '_> {
    /// Deposit a total physical wrench on the original retained liquid mass.
    /// The native Guo-frame update preserves intrinsic density and physical
    /// viscous stress, retaining the source's pressure/rate samples and local
    /// collision coefficients. No constitutive query, time or heat is added.
    /// SI linear/angular balances and cancellation precede publication.
    pub fn apply_wrench_impulse(
        &mut self,
        cells: &[usize],
        wrench: MaterialFreeSurfaceWrenchImpulse3,
        controls: MaterialFreeSurfaceWrenchImpulseControls3,
        gate: &CancelGate,
    ) -> Result<MaterialFreeSurfaceWrenchImpulseReport3, MaterialFreeSurfaceImpulseError> {
        self.flow.apply_wrench_impulse(cells, wrench, controls, gate)
    }

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

impl crate::material::free_surface::impulse::sealed::Sealed
    for MaterialFreeSurfaceRheology3<'_, '_>
{
}

impl MaterialFreeSurfaceMomentum3 for MaterialFreeSurfaceRheology3<'_, '_> {
    fn momentum_scaling(&self) -> MaterialFreeSurfaceScaling<'_> {
        *self.flow.scaling()
    }

    fn prepare_native_impulse_stencil(
        &mut self,
        indices: &[usize],
        controls: MaterialFreeSurfaceImpulseControls3,
        gate: &CancelGate,
    ) -> Result<
        crate::d3q19::freesurface3::FreeSurfaceImpulseStencil3<'_>,
        MaterialFreeSurfaceImpulseError,
    > {
        self.flow
            .prepare_native_impulse_stencil(indices, controls, gate)
    }
}
