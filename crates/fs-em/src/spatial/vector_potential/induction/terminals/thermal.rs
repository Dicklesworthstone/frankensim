//! Source-temperature AC refresh of original physical voltage electrodes.

use super::{
    ResolvedConductorState, Result, SpatialConductorMaterial, SpatialMagneticMaterial,
    SpatialTerminalInduction, WeightedAssemblyLimits, poll,
};
use crate::spatial::vector_potential::induction::harmonic::{
    SpatialHarmonicInductionControls, SpatialHarmonicInductionResponse,
};
use fs_math::c64::C64;

/// A source-refreshed voltage-electrode owner available only to settled AC.
/// The original physical electrode layout, relative references and charge
/// space remain intact. This view cannot publish a transient endpoint or pay
/// magnetic parameter work; its consumer owns the thermal approximation.
pub struct SpatialHarmonicVoltageMaterialTrial<
    'm,
    C: SpatialConductorMaterial = ResolvedConductorState,
> {
    source: SpatialTerminalInduction<'m, C>,
}

impl<'m, C: SpatialConductorMaterial> SpatialTerminalInduction<'m, C> {
    /// Prepare the original material-temperature refresh without advancing A,
    /// changing electrode traces or publishing a new transient source state.
    pub fn prepare_harmonic_materials_trial(
        &self,
        materials: Vec<SpatialMagneticMaterial>,
        conductors: Vec<Option<C>>,
        assembly: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialHarmonicVoltageMaterialTrial<'m, C>> {
        let prepared = self
            .field
            .prepare_materials(materials, conductors, assembly, cancelled)?;
        let mut source = (*self).clone();
        source.field.prepared = prepared;
        poll(cancelled)?;
        Ok(SpatialHarmonicVoltageMaterialTrial { source })
    }
}

impl<'m, C: SpatialConductorMaterial> SpatialHarmonicVoltageMaterialTrial<'m, C> {
    /// Actual source-pinned conductivity at this proposed temperature.
    pub fn conductors(&self) -> &[Option<C>] {
        self.source.conductors()
    }

    /// One native A–phi solve with original electrode reactions and power.
    /// Peak phasors use `exp(-i omega t)`; no clock or history is advanced.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_harmonic_voltage(
        &self,
        voltage_v: &[C64],
        impressed_current_a_m2: &[[C64; 3]],
        velocity_m_s: &[[f64; 3]],
        omega_rad_s: f64,
        controls: SpatialHarmonicInductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialHarmonicInductionResponse<'_, 'm, C>> {
        self.source.solve_harmonic_voltage(
            voltage_v,
            impressed_current_a_m2,
            velocity_m_s,
            omega_rad_s,
            controls,
            cancelled,
        )
    }
}
