//! Exact resumable state for single-load level-set compliance optimization.
//!
//! A checkpoint contains only durable numerical state: the exact level-set
//! field, global iteration ordinal, current augmented-Lagrange multiplier, and
//! immutable problem declaration. The displacement solution is deliberately not
//! serialized; each resumed step re-solves the retained geometry canonically
//! before deriving sensitivities. This avoids replaying prior geometry updates
//! while preserving the global hole-nucleation schedule exactly.

use super::*;
use super::engine::{optimize_compliance_segment, optimize_compliance_segment_controlled};
use std::convert::Infallible;
use std::ops::ControlFlow;

/// Durable optimizer state sufficient for exact deterministic continuation.
#[derive(Debug, Clone)]
pub struct OptimizeCheckpoint {
    geometry: GridSdf,
    fixture: Cantilever,
    settings: OptimizeSettings,
    next_iteration: usize,
    ell: f64,
}

impl OptimizeCheckpoint {
    /// Admit a new optimization checkpoint at iteration zero without solving a
    /// PDE or mutating the supplied geometry.
    ///
    /// # Errors
    /// Returns the same typed admission refusals as the fixed-run optimizer.
    pub fn new(
        geometry: GridSdf,
        fixture: Cantilever,
        settings: OptimizeSettings,
    ) -> Result<Self, CutFemError> {
        Self::restore(geometry, fixture, settings, 0, settings.ell0)
    }

    /// Restore an exact checkpoint from durable geometry and scalar state.
    ///
    /// `next_iteration` is the global zero-based ordinal of the next update;
    /// therefore a checkpoint after three completed updates stores `3`. The
    /// multiplier must be the `ell` published by the last completed row (or the
    /// declared `ell0` at iteration zero).
    ///
    /// # Errors
    /// Refuses an ordinal past the declared total iteration count, invalid
    /// multiplier state, or any invalid geometry/problem declaration.
    pub fn restore(
        geometry: GridSdf,
        fixture: Cantilever,
        settings: OptimizeSettings,
        next_iteration: usize,
        ell: f64,
    ) -> Result<Self, CutFemError> {
        if next_iteration > settings.iterations {
            return Err(invalid_input(format!(
                "checkpoint next iteration {next_iteration} exceeds declared total {}",
                settings.iterations
            )));
        }
        if !(ell.is_finite() && ell >= 0.0) {
            return Err(invalid_input(
                "checkpoint augmented-Lagrange multiplier must be finite and nonnegative",
            ));
        }
        let mut admitted = geometry.clone();
        let mut validation = settings;
        validation.iterations = 0;
        optimize_compliance_segment(
            &mut admitted,
            fixture,
            validation,
            next_iteration,
            ell,
        )?;
        Ok(Self { geometry, fixture, settings, next_iteration, ell })
    }

    /// Exact retained level-set geometry.
    #[must_use]
    pub fn geometry(&self) -> &GridSdf {
        &self.geometry
    }

    /// Global ordinal of the next update.
    #[must_use]
    pub const fn next_iteration(&self) -> usize {
        self.next_iteration
    }

    /// Current augmented-Lagrange multiplier.
    #[must_use]
    pub const fn ell(&self) -> f64 {
        self.ell
    }

    /// Immutable fixture declaration.
    #[must_use]
    pub const fn fixture(&self) -> Cantilever {
        self.fixture
    }

    /// Immutable optimizer declaration, including the total iteration target.
    #[must_use]
    pub const fn settings(&self) -> OptimizeSettings {
        self.settings
    }

    /// Number of declared updates not yet completed.
    #[must_use]
    pub fn remaining(&self) -> usize {
        self.settings.iterations.saturating_sub(self.next_iteration)
    }

    /// Whether the declared iteration target has been reached.
    #[must_use]
    pub fn is_complete(&self) -> bool {
        self.next_iteration == self.settings.iterations
    }

    /// Consume the checkpoint and return the retained geometry.
    #[must_use]
    pub fn into_geometry(self) -> GridSdf {
        self.geometry
    }

    /// Advance exactly one globally numbered update transactionally.
    ///
    /// The checkpoint changes only after the complete candidate geometry has
    /// been canonically solved and its one-row report has been produced. A
    /// refusal therefore leaves geometry, iteration ordinal, and multiplier
    /// unchanged.
    ///
    /// # Errors
    /// Propagates the canonical optimizer/CutFEM refusal for the attempted step.
    pub fn advance_one(&mut self) -> Result<Option<OptimizeReport>, CutFemError> {
        match self.advance_one_controlled(usize::MAX, |_| ControlFlow::<Infallible>::Continue(()))? {
            ControlFlow::Continue(report) => Ok(report),
            ControlFlow::Break(never) => match never {},
        }
    }

    /// Advance one update with interruption inside both canonical CG solves.
    ///
    /// `poll_iters` bounds additional CG iterations between callbacks. The
    /// callback sees the named stage and may return its own stop reason. A
    /// `Break` leaves geometry, multiplier and global ordinal BITWISE unchanged,
    /// even after evolution or a completed candidate solve. Retrying performs
    /// this update again, never replays earlier accepted geometry updates.
    /// `Continue(None)` means the declared update count was already complete.
    ///
    /// Assembly, sensitivity smoothing, advection and individual sparse/vector
    /// operations remain non-preemptible; checks bracket those stages. This is
    /// cooperative interruption, not a hard wall-time or memory guarantee.
    ///
    /// # Errors
    /// Refuses a zero polling interval or propagates a numerical refusal. An
    /// interrupted solve is `Ok(Break(reason))`, not a failed/converged field.
    pub fn advance_one_controlled<B>(
        &mut self,
        poll_iters: usize,
        mut control: impl FnMut(CheckpointStage) -> ControlFlow<B>,
    ) -> Result<ControlFlow<B, Option<OptimizeReport>>, CutFemError> {
        if poll_iters == 0 {
            return Err(invalid_input("checkpoint CG poll interval must be positive"));
        }
        if self.is_complete() {
            return Ok(ControlFlow::Continue(None));
        }
        if let ControlFlow::Break(reason) = control(CheckpointStage::Prepare) {
            return Ok(ControlFlow::Break(reason));
        }
        let mut candidate = self.geometry.clone();
        let mut step_settings = self.settings;
        step_settings.iterations = 1;
        let report = match optimize_compliance_segment_controlled(
            &mut candidate,
            self.fixture,
            step_settings,
            self.next_iteration,
            self.ell,
            poll_iters,
            &mut control,
        )? {
            ControlFlow::Continue(report) => report,
            ControlFlow::Break(reason) => return Ok(ControlFlow::Break(reason)),
        };
        if report.rows.len() != 1 || report.ell.len() != 1 {
            return Err(invalid_input(
                "internal checkpoint step did not produce exactly one evaluated row",
            ));
        }
        let next_iteration = self.next_iteration.checked_add(1)
            .ok_or_else(|| invalid_input("checkpoint iteration ordinal overflow"))?;
        let ell = report.ell[0];
        if !(ell.is_finite() && ell >= 0.0) {
            return Err(invalid_input("checkpoint step produced invalid multiplier state"));
        }
        if let ControlFlow::Break(reason) = control(CheckpointStage::Publish) {
            return Ok(ControlFlow::Break(reason));
        }
        self.geometry = candidate;
        self.next_iteration = next_iteration;
        self.ell = ell;
        Ok(ControlFlow::Continue(Some(report)))
    }

    /// Advance at most `max_steps` updates and concatenate their evaluated
    /// evidence in global iteration order.
    ///
    /// # Errors
    /// Stops on the first refused step, retaining all earlier successfully
    /// committed checkpoint updates and leaving the refused step unpublished.
    pub fn advance(&mut self, max_steps: usize) -> Result<OptimizeReport, CutFemError> {
        let mut combined = OptimizeReport::default();
        for _ in 0..max_steps.min(self.remaining()) {
            let Some(step) = self.advance_one()? else { break };
            append(&mut combined, step);
        }
        Ok(combined)
    }
}

fn append(target: &mut OptimizeReport, mut source: OptimizeReport) {
    target.compliance.append(&mut source.compliance);
    target.volume.append(&mut source.volume);
    target.ell.append(&mut source.ell);
    target.audits.append(&mut source.audits);
    target.events.append(&mut source.events);
    target.snapshots.append(&mut source.snapshots);
    target.load_pad_nodes.append(&mut source.load_pad_nodes);
    target.rows.append(&mut source.rows);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::optimize::optimize_compliance;

    fn beam(level: u32) -> GridSdf {
        GridSdf::from_fn(1usize << level, &|_, y| (y - 0.5).abs() - 0.35)
    }

    fn settings() -> OptimizeSettings {
        OptimizeSettings {
            level: 3,
            iterations: 3,
            move_cells: 0.05,
            nucleation_period: 2,
            hole_radius_cells: 1.0,
            ..OptimizeSettings::default()
        }
    }

    #[test]
    fn checkpoint_steps_match_uninterrupted_global_trajectory() {
        let settings = settings();
        let fixture = Cantilever { load: 1.0, band: 0.125 };
        let initial = beam(settings.level);
        let mut uninterrupted_geometry = initial.clone();
        let uninterrupted = optimize_compliance(
            &mut uninterrupted_geometry,
            fixture,
            settings,
        ).expect("uninterrupted trajectory");

        let mut checkpoint = OptimizeCheckpoint::new(initial, fixture, settings)
            .expect("admitted checkpoint");
        let resumed = checkpoint.advance(settings.iterations)
            .expect("checkpoint trajectory");
        assert!(checkpoint.is_complete());
        assert_eq!(checkpoint.geometry().nodes(), uninterrupted_geometry.nodes());
        assert_eq!(resumed.rows, uninterrupted.rows);
        assert_eq!(resumed.snapshots, uninterrupted.snapshots);
        assert_eq!(resumed.ell, uninterrupted.ell);
    }

    #[test]
    fn durable_midpoint_restore_matches_uninterrupted_tail() {
        let settings = settings();
        let fixture = Cantilever { load: 1.0, band: 0.125 };
        let initial = beam(settings.level);
        let mut checkpoint = OptimizeCheckpoint::new(initial.clone(), fixture, settings)
            .expect("admitted checkpoint");
        let first = checkpoint.advance_one().expect("first step").expect("row");
        let retained_geometry = checkpoint.geometry().clone();
        let retained_iteration = checkpoint.next_iteration();
        let retained_ell = checkpoint.ell();

        let mut restored = OptimizeCheckpoint::restore(
            retained_geometry,
            fixture,
            settings,
            retained_iteration,
            retained_ell,
        ).expect("restored checkpoint");
        let tail = restored.advance(usize::MAX).expect("restored tail");

        let mut uninterrupted_geometry = initial;
        let uninterrupted = optimize_compliance(
            &mut uninterrupted_geometry,
            fixture,
            settings,
        ).expect("uninterrupted trajectory");
        let mut rows = first.rows;
        rows.extend(tail.rows);
        assert_eq!(rows, uninterrupted.rows);
        assert_eq!(restored.geometry().nodes(), uninterrupted_geometry.nodes());
        assert_eq!(restored.ell().to_bits(), uninterrupted.ell.last().unwrap().to_bits());
    }

    #[test]
    fn malformed_restore_refuses_without_touching_input_geometry() {
        let settings = settings();
        let geometry = beam(settings.level);
        let original = geometry.nodes().to_vec();
        assert!(OptimizeCheckpoint::restore(
            geometry.clone(),
            Cantilever { load: 1.0, band: 0.125 },
            settings,
            settings.iterations + 1,
            0.0,
        ).is_err());
        assert_eq!(geometry.nodes(), original.as_slice());
        assert!(OptimizeCheckpoint::restore(
            geometry.clone(),
            Cantilever { load: 1.0, band: 0.125 },
            settings,
            0,
            f64::NAN,
        ).is_err());
        assert_eq!(geometry.nodes(), original.as_slice());
    }

    fn geometry_bits(state: &OptimizeCheckpoint) -> Vec<u64> {
        state.geometry().nodes().iter().map(|value| value.to_bits()).collect()
    }

    #[test]
    fn interrupted_solves_and_late_publication_leave_the_accepted_checkpoint_exact() {
        let settings = settings();
        let fixture = Cantilever { load: 1.0, band: 0.125 };
        let mut retained = OptimizeCheckpoint::new(beam(settings.level), fixture, settings)
            .expect("checkpoint");
        retained.advance_one().expect("retain one real accepted update");
        let mut reference = retained.clone();
        let expected = reference.advance_one().expect("reference update").expect("row");
        for target in [
            CheckpointStage::InitialSolve(2),
            CheckpointStage::CandidateSolve(2),
            CheckpointStage::Publish,
        ] {
            let mut state = retained.clone();
            let stopped = state.advance_one_controlled(2, |stage| {
                if stage == target { ControlFlow::Break(stage) }
                else { ControlFlow::Continue(()) }
            }).expect("interruption is not a numerical error");
            assert!(matches!(stopped, ControlFlow::Break(stage) if stage == target));
            assert_eq!(geometry_bits(&state), geometry_bits(&retained));
            assert_eq!(state.next_iteration(), retained.next_iteration());
            assert_eq!(state.ell().to_bits(), retained.ell().to_bits());
            let ControlFlow::Continue(Some(retry)) = state.advance_one_controlled(
                3, |_| ControlFlow::<()>::Continue(()),
            ).expect("retry") else { panic!("retry did not complete") };
            assert_eq!(retry.rows, expected.rows);
            assert_eq!(retry.snapshots, expected.snapshots);
            assert_eq!(geometry_bits(&state), geometry_bits(&reference));
            assert_eq!(state.ell().to_bits(), reference.ell().to_bits());
            assert_eq!(state.next_iteration(), reference.next_iteration());
        }
    }

    #[test]
    fn controlled_batches_keep_the_global_nucleation_schedule_and_trajectory() {
        let settings = settings();
        let fixture = Cantilever { load: 1.0, band: 0.125 };
        let initial = beam(settings.level);
        let mut full_geometry = initial.clone();
        let full = optimize_compliance(&mut full_geometry, fixture, settings).expect("fixed run");
        let mut state = OptimizeCheckpoint::new(initial, fixture, settings).expect("checkpoint");
        let mut rows = Vec::new();
        for batch in [1, 3, 7] {
            let ControlFlow::Continue(Some(report)) = state.advance_one_controlled(
                batch, |_| ControlFlow::<()>::Continue(()),
            ).expect("controlled update") else { panic!("missing update") };
            rows.extend(report.rows);
        }
        assert!(state.is_complete());
        assert_eq!(rows, full.rows);
        assert_eq!(geometry_bits(&state), full_geometry.nodes().iter().map(|v| v.to_bits()).collect::<Vec<_>>());
    }

    #[test]
    fn zero_poll_refuses_and_preflight_stop_allocates_no_trial() {
        let settings = settings();
        let mut state = OptimizeCheckpoint::new(
            beam(settings.level), Cantilever { load: 1.0, band: 0.125 }, settings,
        ).expect("checkpoint");
        let before = geometry_bits(&state);
        assert!(state.advance_one_controlled(0, |_| -> ControlFlow<()> {
            panic!("invalid control must be refused before polling")
        }).is_err());
        let stopped = state.advance_one_controlled(32, |stage| {
            assert_eq!(stage, CheckpointStage::Prepare);
            ControlFlow::Break("cancelled")
        }).expect("preflight stop");
        assert!(matches!(stopped, ControlFlow::Break("cancelled")));
        assert_eq!(geometry_bits(&state), before);
        assert_eq!(state.next_iteration(), 0);
        assert_eq!(state.ell().to_bits(), settings.ell0.to_bits());
    }
}
