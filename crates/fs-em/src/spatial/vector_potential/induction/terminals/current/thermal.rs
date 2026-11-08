//! Source-temperature AC refresh of original physical current electrodes.

use super::{
    ResolvedConductorState, Result, SpatialConductorMaterial, SpatialCurrentInduction,
    SpatialMagneticMaterial, WeightedAssemblyLimits, poll,
};
use crate::spatial::vector_potential::induction::harmonic::{
    SpatialHarmonicInductionControls, SpatialHarmonicInductionResponse,
};
use fs_math::c64::C64;

/// A source-refreshed total-current electrode owner available only to settled AC.
/// Grouped scalar unknowns and the original component return currents remain
/// native; the spatial current distribution is never imposed. This view does
/// not publish a transient material state or supply magnetic parameter work.
pub struct SpatialHarmonicCurrentMaterialTrial<
    'm,
    C: SpatialConductorMaterial = ResolvedConductorState,
> {
    source: SpatialCurrentInduction<'m, C>,
}

impl<'m, C: SpatialConductorMaterial> SpatialCurrentInduction<'m, C> {
    /// Reuse the original source-temperature refresh and grouped charge space.
    /// Accepted A, physical time, electrode locations and source pins stay fixed.
    pub fn prepare_harmonic_materials_trial(
        &self,
        materials: Vec<SpatialMagneticMaterial>,
        conductors: Vec<Option<C>>,
        assembly: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialHarmonicCurrentMaterialTrial<'m, C>> {
        let prepared = self
            .field
            .prepare_materials(materials, conductors, assembly, cancelled)?;
        let mut source = (*self).clone();
        source.field.prepared = prepared;
        poll(cancelled)?;
        Ok(SpatialHarmonicCurrentMaterialTrial { source })
    }
}

impl<'m, C: SpatialConductorMaterial> SpatialHarmonicCurrentMaterialTrial<'m, C> {
    /// Actual source-pinned conductivity at this proposed temperature.
    pub fn conductors(&self) -> &[Option<C>] {
        self.source.conductors()
    }

    /// One native A–phi solve with actual electrode voltage/current reactions.
    /// Peak phasors use `exp(-i omega t)`; currents balance on each conductor.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_harmonic_current(
        &self,
        incoming_current_a: &[C64],
        impressed_current_a_m2: &[[C64; 3]],
        velocity_m_s: &[[f64; 3]],
        omega_rad_s: f64,
        controls: SpatialHarmonicInductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialHarmonicInductionResponse<'_, 'm, C>> {
        self.source.solve_harmonic_current(
            incoming_current_a,
            impressed_current_a_m2,
            velocity_m_s,
            omega_rad_s,
            controls,
            cancelled,
        )
    }
}
