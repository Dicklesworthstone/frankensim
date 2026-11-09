//! Original distributed injection on the complete structural/acoustic trajectory.

use super::{
    PreparedStructuralInteriorTime, StructuralInteriorTimeError, StructuralInteriorTimeFrame,
    StructuralInteriorTimeState, StructuralInteriorTimeStep, StructuralInteriorTimeWorkspace,
};
use fs_exec::Cx;
use fs_feec::acoustics::helmholtz::heterogeneous::{
    BoundaryVelocitySource, HeterogeneousError, PreparedPressureSourceLoads,
    PressureSourceLoadLimits, VolumeVelocitySource,
};

/// Original load admission and original physical evolution errors stay distinct.
#[derive(Debug)]
pub enum SourceDrivenTimeError<E> {
    Source(HeterogeneousError),
    Evolution(E),
}
impl<E: core::fmt::Display> core::fmt::Display for SourceDrivenTimeError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Source(e) => e.fmt(f),
            Self::Evolution(e) => e.fmt(f),
        }
    }
}
impl<E: std::error::Error + 'static> std::error::Error for SourceDrivenTimeError<E> {}
impl<E> From<HeterogeneousError> for SourceDrivenTimeError<E> {
    fn from(value: HeterogeneousError) -> Self {
        Self::Source(value)
    }
}

/// One original evolution workspace plus one admitted pressure-load scratch.
pub struct SourceDrivenTimeWorkspace<W> {
    pub(super) native: W,
    pub(super) load: Vec<f64>,
}

/// Retains the original whole-system exact propagator; no second stepper exists.
pub struct SourceDrivenStructuralTime<'p, 'm, 'a, 'f, 'c> {
    original: &'p PreparedStructuralInteriorTime<'m, 'a, 'f, 'c>,
    sources: PreparedPressureSourceLoads<'a, 'f, 'c>,
}

impl<'m, 'a, 'f, 'c> PreparedStructuralInteriorTime<'m, 'a, 'f, 'c> {
    /// Bind original physical source columns once. Source limits admit the map
    /// and one load buffer independently of the original complete-state cache.
    pub fn with_sources(
        &self,
        boundary_sources: &[&BoundaryVelocitySource<'_, '_, '_>],
        volume_sources: &[&VolumeVelocitySource<'_, '_, '_>],
        limits: PressureSourceLoadLimits,
        cx: &Cx<'_>,
    ) -> Result<SourceDrivenStructuralTime<'_, 'm, 'a, 'f, 'c>, HeterogeneousError> {
        let sources = self.model.interior.prepare_source_loads(
            boundary_sources,
            volume_sources,
            limits,
            &mut || cx.checkpoint().is_err(),
        )?;
        Ok(SourceDrivenStructuralTime {
            original: self,
            sources,
        })
    }
}

impl<'m, 'a, 'f, 'c> SourceDrivenStructuralTime<'_, 'm, 'a, 'f, 'c> {
    pub fn workspace(
        &self,
    ) -> SourceDrivenTimeWorkspace<StructuralInteriorTimeWorkspace<'m, 'a, 'f, 'c>> {
        SourceDrivenTimeWorkspace {
            native: self.original.workspace(),
            load: vec![0.0; self.sources.pressure_dofs()],
        }
    }

    /// Candidate-returning convenience path. The accepted state is unchanged.
    pub fn step(
        &self,
        initial: &StructuralInteriorTimeState<'m, 'a, 'f, 'c>,
        structural_force: &[f64],
        boundary_amplitudes_m_s: &[f64],
        volume_amplitudes_m3_s: &[f64],
        cx: &Cx<'_>,
    ) -> Result<
        StructuralInteriorTimeStep<'m, 'a, 'f, 'c>,
        SourceDrivenTimeError<StructuralInteriorTimeError>,
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

    /// Exact held-amplitude evolution using cached contractions and no warm
    /// allocation. Source work is already included once in acoustic input as
    /// the actual pressure-trajectory integral dotted with this physical load.
    /// There is no peak-phasor factor, inferred heating or harmonic coefficient
    /// substituted for a causal law. All native fields/clock publish together
    /// only after the original storage, loss, interface and cancellation checks.
    pub fn step_into<'w>(
        &self,
        state: &mut StructuralInteriorTimeState<'m, 'a, 'f, 'c>,
        structural_force: &[f64],
        boundary_amplitudes_m_s: &[f64],
        volume_amplitudes_m3_s: &[f64],
        workspace: &'w mut SourceDrivenTimeWorkspace<
            StructuralInteriorTimeWorkspace<'m, 'a, 'f, 'c>,
        >,
        cx: &Cx<'_>,
    ) -> Result<&'w StructuralInteriorTimeFrame, SourceDrivenTimeError<StructuralInteriorTimeError>>
    {
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
