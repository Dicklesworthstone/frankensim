//! A settled AC view of an original source-temperature material trial.

use super::{Result, SpatialInduction, SpatialInductionMaterialTrial, poll};
use crate::spatial::vector_potential::induction::harmonic::{
    SpatialHarmonicInductionControls, SpatialHarmonicInductionResponse,
};
use fs_material::conductor::{ResolvedConductorState, spatial::SpatialConductorMaterial};
use fs_math::c64::C64;

/// Source-refreshed native operators available only to the harmonic solver.
///
/// Geometry, trace, gauge and electrical source selections are inherited from
/// the original material trial. This view cannot publish a transient endpoint:
/// its temperature refresh has neither advanced A nor paid magnetic parameter
/// work. A thermal consumer must declare its settled-carrier approximation or
/// supply a separate transient/thermodynamic account.
pub struct SpatialHarmonicMaterialTrial<'m, C: SpatialConductorMaterial = ResolvedConductorState> {
    source: SpatialInduction<'m, C>,
}

impl<'m, C: SpatialConductorMaterial> SpatialInductionMaterialTrial<'_, 'm, C> {
    /// Reuse the admitted endpoint operators for settled frequency responses.
    /// No geometry or material matrix is assembled a second time. The original
    /// accepted history stays unchanged and is not exposed through this view.
    pub fn into_harmonic(
        self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialHarmonicMaterialTrial<'m, C>> {
        poll(cancelled)?;
        let mut source = (*self.original).clone();
        source.prepared = self.prepared;
        poll(cancelled)?;
        Ok(SpatialHarmonicMaterialTrial { source })
    }
}

impl<'m, C: SpatialConductorMaterial> SpatialHarmonicMaterialTrial<'m, C> {
    /// The complete original electrical sources at this proposed temperature.
    pub fn conductors(&self) -> &[Option<C>] {
        self.source.conductors()
    }

    /// Resolve the original sparse A–phi equations on the refreshed sources.
    /// All original charge, Ampere, power, motion and skin-depth gates apply.
    /// Peak phasors use `exp(-i omega t)`; the original clock is not advanced.
    pub fn solve_harmonic(
        &self,
        impressed_current_a_m2: &[[C64; 3]],
        velocity_m_s: &[[f64; 3]],
        omega_rad_s: f64,
        controls: SpatialHarmonicInductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialHarmonicInductionResponse<'_, 'm, C>> {
        self.source.solve_harmonic(
            impressed_current_a_m2,
            velocity_m_s,
            omega_rad_s,
            controls,
            cancelled,
        )
    }
}
