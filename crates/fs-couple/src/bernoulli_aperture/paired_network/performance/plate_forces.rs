//! Physical source-mesh loads for the existing independently moving reed player.
use super::{PairedAperturePerformance, ForceSchedule, GestureCompileError, invalid};
use crate::bernoulli_aperture::dynamic::force::PlateForcePort;
use crate::bernoulli_aperture::performance::force::{ApertureForceEvent, ApertureForceProgram};

impl PairedAperturePerformance {
    /// Bind physical node/patch loads to ONE retained source blade. Total input
    /// force [N] is on that blade, not divided across the pair or interpreted as
    /// pressure. The source plate projects force AND actuator velocity through
    /// its actual mode, so material/geometry substitutions change this coupling.
    /// Call once per desired blade; each may have different sites and events.
    ///
    /// Reuses the scalar valve's original footprint and force-schedule compiler.
    /// Unsorted events retain stable equal-time input order (last write per port
    /// wins). Overlapping ports add. Zero releases only the applied load; neither
    /// blade's material memory, contact energy nor traveling waves are reset.
    ///
    /// Up to 32 ports per blade. max_controls bounds pressure assignments plus
    /// BOTH raw force programs and boundary assignments. max_compile_work
    /// independently bounds pressure sampling and COMBINED auxiliary boundary
    /// admission and force projection/grouping work across blades.
    /// Already-generalized `with_blade_forces` controls are an alternative input
    /// representation and cannot be combined with source-mesh programs.
    ///
    /// # Errors
    /// Invalid blade/site/force/event, missing source plates, duplicate or late
    /// attachment, poisoned renderer, mixed force representations, or exhausted
    /// shared compilation/control capacity. No physical force is clipped.
    pub fn with_plate_forces(mut self, blade: usize, program: ApertureForceProgram)
        -> Result<Self, GestureCompileError>
    {
        if blade>=2 || self.poisoned || self.system.aperture().accepted_steps()!=0
            || self.forces.is_some() || self.plate_forces[blade].is_some() {
            return Err(invalid("plate forces require one attachment per cold blade and no already-generalized force program"));
        }
        let (events, work) = self.remaining_control_budgets()?;
        let plates=self.system.aperture().plates()
            .ok_or_else(||invalid("physical blade footprints require the original source plates"))?;
        self.plate_forces[blade]=Some(ForceSchedule::compile(&plates[blade],program,
            self.config.samples,events,work)?);
        Ok(self)
    }

    fn plate_force_schedule(&self, blade: usize) -> Option<&ForceSchedule> {
        self.plate_forces.get(blade).and_then(Option::as_ref)
    }
    /// Work-conjugate force/velocity ports projected from this original blade.
    /// Empty for an absent blade/program or already-generalized controls.
    #[must_use]
    pub fn plate_force_ports(&self, blade: usize) -> &[PlateForcePort] {
        self.plate_force_schedule(blade).map_or(&[],|f|f.ports.as_slice())
    }
    /// Physical plate-load assignments published with successful source samples.
    #[must_use]
    pub fn applied_plate_force_controls(&self, blade: usize) -> &[ApertureForceEvent] {
        self.plate_force_schedule(blade).map_or(&[],|f|&f.events[..f.applied])
    }
    /// Remaining physical plate-load assignments, in stable execution order.
    #[must_use]
    pub fn pending_plate_force_controls(&self, blade: usize) -> &[ApertureForceEvent] {
        self.plate_force_schedule(blade).map_or(&[],|f|&f.events[f.applied..])
    }
    /// Accepted generalized forces on both opening coordinates [N]. These are
    /// either the supplied generalized values or the actual source-mesh sums;
    /// never a pressure, gain or force silently shared between the two blades.
    #[must_use]
    pub const fn held_generalized_forces_n(&self) -> [f64; 2] { self.held_forces_n }
}
