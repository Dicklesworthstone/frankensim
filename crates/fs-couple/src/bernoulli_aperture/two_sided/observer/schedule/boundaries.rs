//! Shared physical key programs on the original two-sided mechanical clock.
use super::{BoundarySchedule, GestureCompileError, ScheduledTwoSidedPressure, ValveSide, invalid};
use crate::bernoulli_aperture::network::boundaries::{BoundaryDrive, BoundaryNetworkFrame};
use crate::bernoulli_aperture::performance::boundaries::{BoundaryControl, BoundaryProgram};

impl ScheduledTwoSidedPressure {
    /// Attach one complete program to an already keyed acoustic side. Events
    /// address mechanical samples, including substeps inside one audio frame.
    /// All keys retain their original material/contact and chimney histories.
    /// Zero force releases the actuator without clearing any physical state.
    ///
    /// Initial values and changes on BOTH sides share max_controls with pressure
    /// and source-plate forces. Boundary admission and source-force projection
    /// share max_compile_work, independently of the existing pressure compiler.
    /// Attachment order cannot renew either budget. A keyed side without an
    /// explicit program refuses playback before any physical sample advances.
    pub fn with_boundary_controls(
        mut self,
        side: ValveSide,
        program: BoundaryProgram,
    ) -> Result<Self, GestureCompileError> {
        if self.poisoned || self.source.poisoned
            || self.source.mechanical_samples_rendered() != 0
            || self.boundaries[side.index()].is_some() {
            return Err(invalid("boundary controls require one cold attachment per keyed acoustic side"));
        }
        let (assignments, work) = self.remaining_control_budgets()?;
        let boundaries = self.source.system.boundary_apertures(side)
            .ok_or_else(|| invalid("boundary controls require the side's original physical apertures"))?;
        let samples = self.config.samples.checked_mul(self.source.mechanics_per_output() as u64)
            .ok_or_else(|| invalid("boundary control mechanical window overflow"))?;
        self.boundaries[side.index()] = Some(BoundarySchedule::compile(
            boundaries, program, samples, assignments, work,
        )?);
        Ok(self)
    }

    pub(super) fn remaining_control_budgets(&self) -> Result<(usize, u64), GestureCompileError> {
        let mut assignments = self.config.max_controls.checked_sub(self.controls.len())
            .ok_or_else(|| invalid("pressure assignments exceed the combined storage budget"))?;
        let mut work = self.config.max_compile_work;
        for force in self.plate_forces.iter().flatten() {
            assignments = assignments.checked_sub(force.events.len())
                .ok_or_else(|| invalid("source forces exceed the combined control storage budget"))?;
            work = work.checked_sub(force.work)
                .ok_or_else(|| invalid("source forces exceed the auxiliary compilation work budget"))?;
        }
        for boundary in self.boundaries.iter().flatten() {
            assignments = assignments.checked_sub(boundary.assignment_count())
                .ok_or_else(|| invalid("boundaries exceed the combined control storage budget"))?;
            work = work.checked_sub(boundary.work)
                .ok_or_else(|| invalid("boundaries exceed the auxiliary compilation work budget"))?;
        }
        Ok((assignments, work))
    }

    /// Last accepted held drives, or explicit initial drives before playback.
    #[must_use]
    pub fn held_boundary_drives(&self, side: ValveSide) -> &[BoundaryDrive] {
        self.boundaries[side.index()].as_ref().map_or(&[], BoundarySchedule::held)
    }
    /// Events published with successfully observed mechanical samples only.
    #[must_use]
    pub fn applied_boundary_controls(&self, side: ValveSide) -> &[BoundaryControl] {
        self.boundaries[side.index()].as_ref().map_or(&[], BoundarySchedule::applied)
    }
    /// Remaining events in mechanical (sample, node) order.
    #[must_use]
    pub fn pending_boundary_controls(&self, side: ValveSide) -> &[BoundaryControl] {
        self.boundaries[side.index()].as_ref().map_or(&[], BoundarySchedule::pending)
    }
    /// Last key/graph work account ingested by the original pressure receiver.
    #[must_use]
    pub fn last_boundary_frame(&self, side: ValveSide) -> Option<&BoundaryNetworkFrame> {
        self.source.last_frame()?.boundaries[side.index()].as_ref()
    }
}
