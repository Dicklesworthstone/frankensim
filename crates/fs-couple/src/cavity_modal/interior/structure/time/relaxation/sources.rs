//! Physical source columns retain every original Maxwell history and work term.

use super::super::sources::{SourceDrivenTimeError, SourceDrivenTimeWorkspace};
use super::{
    PreparedStructuralInteriorRelaxationTime, StructuralInteriorRelaxationTimeError,
    StructuralInteriorRelaxationTimeFrame, StructuralInteriorRelaxationTimeState,
    StructuralInteriorRelaxationTimeStep, StructuralInteriorRelaxationTimeWorkspace,
};
use fs_exec::Cx;
use fs_feec::acoustics::helmholtz::heterogeneous::{
    BoundaryVelocitySource, HeterogeneousError, PreparedPressureSourceLoads,
    PressureSourceLoadLimits, VolumeVelocitySource,
};

pub struct SourceDrivenStructuralRelaxationTime<'p, 'm, 'a, 'f, 'c> {
    original: &'p PreparedStructuralInteriorRelaxationTime<'m, 'a, 'f, 'c>,
    sources: PreparedPressureSourceLoads<'a, 'f, 'c>,
}

impl<'m, 'a, 'f, 'c> PreparedStructuralInteriorRelaxationTime<'m, 'a, 'f, 'c> {
    /// Cache the original physical load columns without replacing any native
    /// Maxwell strain map, relaxation time, history or pressure coordinate.
    pub fn with_sources(
        &self,
        boundary_sources: &[&BoundaryVelocitySource<'_, '_, '_>],
        volume_sources: &[&VolumeVelocitySource<'_, '_, '_>],
        limits: PressureSourceLoadLimits,
        cx: &Cx<'_>,
    ) -> Result<SourceDrivenStructuralRelaxationTime<'_, 'm, 'a, 'f, 'c>, HeterogeneousError> {
        let sources = self.model.interior.prepare_source_loads(
            boundary_sources,
            volume_sources,
            limits,
            &mut || cx.checkpoint().is_err(),
        )?;
        Ok(SourceDrivenStructuralRelaxationTime {
            original: self,
            sources,
        })
    }
}

impl<'m, 'a, 'f, 'c> SourceDrivenStructuralRelaxationTime<'_, 'm, 'a, 'f, 'c> {
    pub fn workspace(
        &self,
    ) -> SourceDrivenTimeWorkspace<StructuralInteriorRelaxationTimeWorkspace<'m, 'a, 'f, 'c>> {
        SourceDrivenTimeWorkspace {
            native: self.original.workspace(),
            load: vec![0.0; self.sources.pressure_dofs()],
        }
    }

    pub fn step(
        &self,
        initial: &StructuralInteriorRelaxationTimeState<'m, 'a, 'f, 'c>,
        structural_force: &[f64],
        boundary_amplitudes_m_s: &[f64],
        volume_amplitudes_m3_s: &[f64],
        cx: &Cx<'_>,
    ) -> Result<
        StructuralInteriorRelaxationTimeStep<'m, 'a, 'f, 'c>,
        SourceDrivenTimeError<StructuralInteriorRelaxationTimeError>,
    > {
        let mut load = vec![0.0; self.sources.pressure_dofs()];
        self.sources.project_into(
            self.original.model.interior,
            boundary_amplitudes_m_s,
            volume_amplitudes_m3_s,
            &mut load,
            &mut || cx.checkpoint().is_err(),
        )?;
        self.original
            .step(initial, structural_force, &load, cx)
            .map_err(SourceDrivenTimeError::Evolution)
    }

    /// Advance the original fields, all material memories and clock atomically.
    /// Cached source projection and native step_into allocate nothing. Held
    /// source work is included once in acoustic_input_work_j; every original
    /// branch retains its independently integrated work, storage and heat.
    pub fn step_into<'w>(
        &self,
        state: &mut StructuralInteriorRelaxationTimeState<'m, 'a, 'f, 'c>,
        structural_force: &[f64],
        boundary_amplitudes_m_s: &[f64],
        volume_amplitudes_m3_s: &[f64],
        workspace: &'w mut SourceDrivenTimeWorkspace<
            StructuralInteriorRelaxationTimeWorkspace<'m, 'a, 'f, 'c>,
        >,
        cx: &Cx<'_>,
    ) -> Result<
        &'w StructuralInteriorRelaxationTimeFrame,
        SourceDrivenTimeError<StructuralInteriorRelaxationTimeError>,
    > {
        self.sources.project_into(
            self.original.model.interior,
            boundary_amplitudes_m_s,
            volume_amplitudes_m3_s,
            &mut workspace.load,
            &mut || cx.checkpoint().is_err(),
        )?;
        self.original
            .step_into(
                state,
                structural_force,
                &workspace.load,
                &mut workspace.native,
                cx,
            )
            .map_err(SourceDrivenTimeError::Evolution)
    }
}
