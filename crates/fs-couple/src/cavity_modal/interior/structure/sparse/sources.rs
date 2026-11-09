//! Physical injection through the original sparse full-coordinate midpoint owner.

use super::super::time::sources::SourceDrivenTimeError;
use super::{
    PreparedSparseStructuralInterior, SparseStructuralTimeError, SparseStructuralTimeState,
    SparseStructuralTimeStep,
};
use fs_exec::Cx;
use fs_feec::acoustics::helmholtz::heterogeneous::{
    BoundaryVelocitySource, HeterogeneousError, PreparedPressureSourceLoads,
    PressureSourceLoadLimits, VolumeVelocitySource,
};

pub struct SourceDrivenSparseStructuralInterior<'p, 'a, 'f, 'c> {
    original: &'p PreparedSparseStructuralInterior<'a, 'f, 'c>,
    sources: PreparedPressureSourceLoads<'a, 'f, 'c>,
}

impl<'a, 'f, 'c> PreparedSparseStructuralInterior<'a, 'f, 'c> {
    /// Bind complete original source columns without constructing dense spatial
    /// M/K/C forms. Additional storage is source-count times pressure DOFs,
    /// separately bounded before allocation, plus one load vector per query.
    pub fn with_sources(
        &self,
        boundary_sources: &[&BoundaryVelocitySource<'_, '_, '_>],
        volume_sources: &[&VolumeVelocitySource<'_, '_, '_>],
        limits: PressureSourceLoadLimits,
        cx: &Cx<'_>,
    ) -> Result<SourceDrivenSparseStructuralInterior<'_, 'a, 'f, 'c>, HeterogeneousError> {
        let sources = self.interior.prepare_source_loads(
            boundary_sources,
            volume_sources,
            limits,
            &mut || cx.checkpoint().is_err(),
        )?;
        Ok(SourceDrivenSparseStructuralInterior {
            original: self,
            sources,
        })
    }
}

impl<'p> SourceDrivenSparseStructuralInterior<'p, '_, '_, '_> {
    /// Supply actual midpoint source amplitudes to the same native sparse
    /// solve. Its original full equations, independent physical work and state
    /// binding remain authoritative. The immutable initial state survives every
    /// refusal; source work is included once, without a peak-phasor factor.
    pub fn step(
        &self,
        initial: &SparseStructuralTimeState<'p>,
        structural_force: &[f64],
        boundary_amplitudes_m_s: &[f64],
        volume_amplitudes_m3_s: &[f64],
        cx: &Cx<'_>,
    ) -> Result<SparseStructuralTimeStep<'p>, SourceDrivenTimeError<SparseStructuralTimeError>>
    {
        let mut load = vec![0.0; self.sources.pressure_dofs()];
        self.sources.project_into(
            self.original.interior,
            boundary_amplitudes_m_s,
            volume_amplitudes_m3_s,
            &mut load,
            &mut || cx.checkpoint().is_err(),
        )?;
        self.original
            .step(initial, structural_force, &load, cx)
            .map_err(SourceDrivenTimeError::Evolution)
    }
}
