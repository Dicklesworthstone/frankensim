//! Original molecular calorics in full sparse three-dimensional pressure waves.
//!
//! Each conforming material region retains its RelaxationState, exact positive
//! compressibility realization and original density. The frozen beta_inf builds
//! the FE pressure mass once; explicit volume memories supply the relaxed
//! increment once. No equilibrium-mass double count, fitted attenuation or
//! additional bulk viscosity is introduced. This owner resolves neither shear
//! nor thermal diffusion; the existing resolved thermoviscous owner does.

use fs_feec::acoustics::helmholtz::heterogeneous::time::{
    ResistiveBoundary,
    bulk::{BulkInteriorError, BulkMemoryArm, BulkMemoryRegion, BulkTimeLimits, PreparedBulkTime},
};
use fs_feec::acoustics::helmholtz::heterogeneous::{
    HeterogeneousHelmholtzFe, SpatialResistiveBoundary,
};
use fs_feec::acoustics::helmholtz::{AcousticFeError, HelmholtzFe, Medium};
use fs_material::{
    MaterialError,
    relaxation::{RelaxationState, compliance::CompressibilityRelaxation},
};
use fs_math::det;

/// Original frozen molecular source occupying complete original tetrahedra.
#[derive(Clone, Debug)]
pub struct MolecularBulkRegion<'s> {
    pub tetrahedra: Vec<usize>,
    pub source: &'s RelaxationState,
}
#[derive(Debug)]
pub enum MolecularBulkError {
    Invalid(&'static str),
    Cancelled,
    Material(MaterialError),
    Geometry(AcousticFeError),
    Time(BulkInteriorError),
}
impl core::fmt::Display for MolecularBulkError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(v) => write!(f, "molecular pressure interior: {v}"),
            Self::Cancelled => f.write_str("molecular pressure interior cancelled"),
            Self::Material(e) => e.fmt(f),
            Self::Geometry(e) => e.fmt(f),
            Self::Time(e) => e.fmt(f),
        }
    }
}
impl std::error::Error for MolecularBulkError {}
impl From<MaterialError> for MolecularBulkError {
    fn from(e: MaterialError) -> Self {
        Self::Material(e)
    }
}
impl From<AcousticFeError> for MolecularBulkError {
    fn from(e: AcousticFeError) -> Self {
        Self::Geometry(e)
    }
}
impl From<BulkInteriorError> for MolecularBulkError {
    fn from(e: BulkInteriorError) -> Self {
        Self::Time(e)
    }
}

/// Original sources and full frozen-density/compressibility spatial forms.
pub struct MolecularBulkInterior<'f, 'c, 's> {
    sources: Vec<MolecularBulkRegion<'s>>,
    laws: Vec<CompressibilityRelaxation>,
    assignments: Vec<usize>,
    model: HeterogeneousHelmholtzFe<'f, 'c>,
    memory_regions: Vec<BulkMemoryRegion>,
}
impl<'f, 'c, 's> MolecularBulkInterior<'f, 'c, 's> {
    /// Bind a complete disjoint conforming partition, with no default medium.
    ///
    /// The background pressure must agree between sources; stationary pressure
    /// jumps would require an additional equilibrium-force model. Temperature
    /// and composition may differ and remain frozen at their original values.
    /// At most 256 regions are admitted here; per-step resource bounds remain
    /// explicit in the prepared sparse owner. No source is queried during steps.
    pub fn new(
        geometry: &'f HelmholtzFe<'c>,
        regions: &[MolecularBulkRegion<'s>],
        quadrature_order: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, MolecularBulkError> {
        if cancelled() {
            return Err(MolecularBulkError::Cancelled);
        }
        if regions.is_empty() || regions.len() > 256 {
            return Err(MolecularBulkError::Invalid(
                "one through 256 original molecular regions required",
            ));
        }
        let mut assignments = vec![usize::MAX; geometry.tetrahedra().len()];
        let mut laws = Vec::with_capacity(regions.len());
        let mut media = Vec::with_capacity(regions.len());
        let mut memory_regions = Vec::new();
        let pressure = regions[0].source.gas().pressure;
        for (index, region) in regions.iter().enumerate() {
            if cancelled() {
                return Err(MolecularBulkError::Cancelled);
            }
            if region.tetrahedra.is_empty() {
                return Err(MolecularBulkError::Invalid(
                    "every molecular region must occupy original tetrahedra",
                ));
            }
            let gas = region.source.gas();
            gas.validate_coherence()?;
            if gas.pressure != pressure {
                return Err(MolecularBulkError::Invalid(
                    "original molecular regions require a common stationary background pressure",
                ));
            }
            for &element in &region.tetrahedra {
                if element >= assignments.len() || assignments[element] != usize::MAX {
                    return Err(MolecularBulkError::Invalid(
                        "molecular regions must be disjoint original tetrahedra",
                    ));
                }
                assignments[element] = index;
            }
            let law = region.source.compressibility_relaxation()?;
            let c = det::sqrt(law.frozen_compressibility_pa_inv().recip() / gas.density);
            if !c.is_finite() || c <= 0.0 {
                return Err(MolecularBulkError::Invalid(
                    "original frozen molecular sound speed is unrepresentable",
                ));
            }
            media.push(Medium {
                density_kg_m3: gas.density,
                sound_speed_m_s: c,
            });
            if !law.arms().is_empty() {
                memory_regions.push(BulkMemoryRegion {
                    tetrahedra: region.tetrahedra.clone(),
                    arms: law
                        .arms()
                        .iter()
                        .map(|a| BulkMemoryArm {
                            delta_compressibility_pa_inv: a.delta_compressibility_pa_inv,
                            relaxation_time_s: a.relaxation_time_s,
                        })
                        .collect(),
                });
            }
            laws.push(law);
        }
        if assignments.iter().any(|&i| i == usize::MAX) {
            return Err(MolecularBulkError::Invalid(
                "molecular sources must cover every original tetrahedron",
            ));
        }
        let model = HeterogeneousHelmholtzFe::new(
            geometry,
            quadrature_order,
            &mut |element, _| Ok::<_, MolecularBulkError>(media[assignments[element]]),
            cancelled,
        )?;
        if cancelled() {
            return Err(MolecularBulkError::Cancelled);
        }
        Ok(Self {
            sources: regions.to_vec(),
            laws,
            assignments,
            model,
            memory_regions,
        })
    }
    pub const fn model(&self) -> &HeterogeneousHelmholtzFe<'f, 'c> {
        &self.model
    }
    pub fn sources(&self) -> &[MolecularBulkRegion<'s>] {
        &self.sources
    }
    pub fn compressibility_laws(&self) -> &[CompressibilityRelaxation] {
        &self.laws
    }
    pub fn region_for_tetrahedron(&self, element: usize) -> Option<usize> {
        self.assignments.get(element).copied()
    }
    /// Full pressure-field propagation with independent original memory heat.
    /// Midpoint is passive and second order; wave phase and reservoir times
    /// require independent time refinement. It is not an exact/L-stable scheme.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_time<'m>(
        &'m self,
        duration_s: f64,
        boundaries: &[ResistiveBoundary],
        spatial: &[&SpatialResistiveBoundary<'_, '_, '_>],
        pressure_release_faces: &[[u32; 3]],
        limits: BulkTimeLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedBulkTime<'m, 'f, 'c>, MolecularBulkError> {
        Ok(self.model.prepare_bulk_time(
            duration_s,
            &self.memory_regions,
            boundaries,
            spatial,
            pressure_release_faces,
            limits,
            cancelled,
        )?)
    }
}

#[cfg(test)]
mod tests;
