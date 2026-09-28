//! Cancellable feasible-baseline construction and exact accepted-state recovery.
use super::*;

/// Cooperative boundaries during baseline projection and independent replay.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProjectedSetupStage {
    /// Before validation, cloning the source or allocating numerical buffers.
    Prepare,
    /// Projection or verification of the declared material/fixed-node policy.
    Projection(VolumeProjectionStage),
    /// Independent final mechanics, including bounded batches of CG iterations.
    Evaluation(DesignEvaluationStage),
    /// Feasible baseline is complete but the optimizer has not been returned.
    Publish,
}

impl ProjectedOptimizer {
    /// Exact prescribed nodes retained by this owner, in canonical order.
    #[must_use]
    pub fn fixed_nodes(&self) -> &[(usize, f64)] { &self.fixed }

    /// The admitted numerical area policy; refinement cannot silently relax it.
    #[must_use]
    pub const fn projection_settings(&self) -> VolumeProjectionSettings { self.projection }

    /// Establish a feasible baseline with cancellation during projection and CG.
    ///
    /// The source is borrowed: stopping discards only staged work. The same
    /// fixed-node and area gates as [`Self::new`] apply. Baseline preparation
    /// includes a second, no-change projection check before the independent
    /// solve; its `Prepare` boundary is also delivered to the controller.
    ///
    /// # Errors
    /// Refuses invalid controls, unattainable area, or a failed baseline solve.
    pub fn new_controlled<B>(
        geometry: &GridSdf,
        fixture: Cantilever,
        settings: OptimizeSettings,
        fixed: Vec<(usize, f64)>,
        projection: VolumeProjectionSettings,
        controls: ProjectedSettings,
        mut control: impl FnMut(ProjectedSetupStage) -> ControlFlow<B>,
    ) -> Result<ControlFlow<B, Self>, CutFemError> {
        validate_controls(settings, projection, controls)?;
        if let ControlFlow::Break(reason) = control(ProjectedSetupStage::Prepare) {
            return Ok(ControlFlow::Break(reason));
        }
        if geometry.n() != (1usize << settings.level) {
            return Err(refused("projected descent requires a level-matched input lattice"));
        }
        let _ = OptimizeCheckpoint::new(geometry.clone(), fixture, settings)?;
        let mut trial = geometry.clone();
        if let ControlFlow::Break(reason) = project_material_volume_controlled(
            &mut trial, settings.level, &fixed, projection,
            |stage| control(ProjectedSetupStage::Projection(stage)),
        )? {
            return Ok(ControlFlow::Break(reason));
        }
        let checkpoint = OptimizeCheckpoint::new(trial, fixture, settings)?;
        Self::from_checkpoint_controlled(&checkpoint, fixed, projection, controls, control)
    }

    /// Re-admit an exact accepted checkpoint with cancellation inside its replay.
    ///
    /// Verification may not change a single node bit. The source checkpoint is
    /// borrowed and survives cancellation or refusal; no partial optimizer is
    /// returned. The resulting baseline describes this resumed segment only.
    ///
    /// # Errors
    /// Refuses changed fixed nodes, infeasible area, invalid policy or PDE failure.
    pub fn from_checkpoint_controlled<B>(
        checkpoint: &OptimizeCheckpoint,
        fixed: Vec<(usize, f64)>,
        projection: VolumeProjectionSettings,
        controls: ProjectedSettings,
        mut control: impl FnMut(ProjectedSetupStage) -> ControlFlow<B>,
    ) -> Result<ControlFlow<B, Self>, CutFemError> {
        validate_controls(checkpoint.settings(), projection, controls)?;
        if let ControlFlow::Break(reason) = control(ProjectedSetupStage::Prepare) {
            return Ok(ControlFlow::Break(reason));
        }
        let mut checked = checkpoint.geometry().clone();
        if let ControlFlow::Break(reason) = project_material_volume_controlled(
            &mut checked, checkpoint.settings().level, &fixed, projection,
            |stage| control(ProjectedSetupStage::Projection(stage)),
        )? {
            return Ok(ControlFlow::Break(reason));
        }
        if checked.nodes().iter().zip(checkpoint.geometry().nodes())
            .any(|(left, right)| left.to_bits() != right.to_bits())
        {
            return Err(refused("resumed projected state does not already satisfy fixed nodes and material area"));
        }
        let current = match evaluate_compliance_design_controlled(
            checkpoint.geometry(), checkpoint.fixture(), checkpoint.settings(), controls.poll_iters,
            |stage| control(ProjectedSetupStage::Evaluation(stage)),
        )? {
            ControlFlow::Continue(current) => current,
            ControlFlow::Break(reason) => return Ok(ControlFlow::Break(reason)),
        };
        if (current.volume - projection.target).abs() > projection.tolerance {
            return Err(refused("canonical baseline area disagrees with the projected material constraint"));
        }
        if let ControlFlow::Break(reason) = control(ProjectedSetupStage::Publish) {
            return Ok(ControlFlow::Break(reason));
        }
        Ok(ControlFlow::Continue(Self {
            checkpoint: checkpoint.clone(), fixed, projection, controls, baseline: current, current,
        }))
    }
}
