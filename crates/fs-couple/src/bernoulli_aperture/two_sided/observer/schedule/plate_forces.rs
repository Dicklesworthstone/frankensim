//! Source-mesh actuation through the original force compiler and two-sided solve.
use super::{ForceSchedule, GestureCompileError, ScheduledTwoSidedPressure, invalid};
use crate::bernoulli_aperture::dynamic::force::PlateForcePort;
use crate::bernoulli_aperture::performance::force::{ApertureForceEvent, ApertureForceProgram};
use crate::bernoulli_aperture::two_sided::LoadedAperture;

impl ScheduledTwoSidedPressure {
    /// Schedule the physical force footprint already defined by one original
    /// material modal plate. The source's node weights and every signed modal
    /// force/velocity projection remain unchanged. Each event addresses port 0
    /// of this blade and supplies force in newtons, positive toward opening.
    ///
    /// Events execute before their mechanical sample, including within an audio
    /// frame, with stable last-write semantics for equal-time assignments.
    /// They share the existing pressure/key/force budgets and publish only with
    /// accepted physics. Zero releases the actuator without resetting any state.
    /// These are prescribed loads, not a lip/tongue contact or feedback model.
    ///
    /// # Errors
    /// Non-material or absent blade, absent source force footprint, invalid
    /// event, late/duplicate attachment, held-force conflict, or shared budgets.
    pub fn with_material_plate_force_controls(mut self, blade: usize,
        events: Vec<ApertureForceEvent>) -> Result<Self, GestureCompileError>
    {
        if blade >= 2 || self.poisoned || self.source.poisoned
            || self.source.mechanical_samples_rendered() != 0
            || self.plate_forces[blade].is_some() || self.source.held.blade_forces_n[blade] != 0.0 {
            return Err(invalid("material source forces require one cold attachment per unloaded blade"));
        }
        let LoadedAperture::Modal(source) = self.source.system().valve() else {
            return Err(invalid("source-defined force controls require original material modal plates"));
        };
        if source.material_plates().get(blade).is_none() {
            return Err(invalid("source force names an absent material plate"));
        }
        let mut probe = [0.0; 2];
        probe[blade] = 1.0;
        source.validate_forces(probe).map_err(|error|
            GestureCompileError::Render(crate::render::RenderError::Voice(error)))?;
        let (assignments, work) = self.remaining_control_budgets()?;
        let samples = self.config.samples.checked_mul(self.source.mechanics_per_output() as u64)
            .ok_or_else(|| invalid("material source force mechanical window overflow"))?;
        self.plate_forces[blade] = Some(ForceSchedule::source_port(events, samples, assignments, work)?);
        Ok(self)
    }

    /// Bind physical node/patch loads to one actual source plate before playback.
    /// Blade 0 is the scalar opening coordinate for a single-coordinate source;
    /// 0 and 1 identify the independent plates of a paired source. No source mesh
    /// is inferred for an authored scalar valve. For a symmetric opposed scalar
    /// reduction the original port retains its total-force/equal-split meaning.
    ///
    /// Events use MECHANICAL sample addresses, not output frames. Loads at a
    /// substep inside an audio frame change that solve without losing preceding
    /// flows. The existing force compiler preserves unsorted input, stable last
    /// assignment per port at equal time, and signed sums of overlapping ports.
    /// Zero releases an actuator only; material, contact and both ducts keep state.
    ///
    /// Each bound blade starts with zero prescribed force. Nonzero manual force
    /// on that blade refuses rather than being silently replaced; an unbound
    /// blade's already-held force is retained. Further manual changes are refused
    /// once a source-mesh program is attached. This is prescribed actuation, not
    /// a tongue/lip contact law or an additional source of pressure work.
    ///
    /// max_controls bounds pressure assignments, BOTH raw force programs and
    /// every moving-boundary drive. max_compile_work separately bounds combined
    /// force projection/grouping, typed sampling and boundary compilation;
    /// it remains the pressure compiler's independent work allowance. The same
    /// source-owned force/velocity projection and numerical solve do all physics.
    ///
    /// # Errors
    /// Missing source plate, invalid blade/site/force/event, late or duplicate
    /// attachment, conflicting held load, poisoned source, or combined budgets.
    pub fn with_plate_forces(mut self, blade: usize, program: ApertureForceProgram)
        -> Result<Self, GestureCompileError>
    {
        if blade >= 2 || self.poisoned || self.source.poisoned
            || self.source.mechanical_samples_rendered() != 0
            || self.plate_forces[blade].is_some() || self.source.held.blade_forces_n[blade] != 0.0 {
            return Err(invalid("plate forces require one cold attachment per existing unloaded blade"));
        }
        let plate = match self.source.system().valve() {
            LoadedAperture::Single(a) if blade == 0 => a.plate_reduction(),
            LoadedAperture::Paired(a) => a.plates().map(|p| &p[blade]),
            _ => None,
        }.ok_or_else(|| invalid("physical force footprints require this blade's retained source plate"))?;
        let (events, work) = self.remaining_control_budgets()?;
        let samples = self.config.samples.checked_mul(self.source.mechanics_per_output() as u64)
            .ok_or_else(|| invalid("physical force mechanical window overflow"))?;
        self.plate_forces[blade] = Some(ForceSchedule::compile(plate, program, samples, events, work)?);
        Ok(self)
    }

    fn plate_force_schedule(&self, blade: usize) -> Option<&ForceSchedule> {
        self.plate_forces.get(blade).and_then(Option::as_ref)
    }
    /// Work-conjugate ports on the retained plate, in authored footprint order.
    /// Empty for an unbound blade or a material-modal source-defined footprint;
    /// that footprint remains in the original material plate's reduction options.
    #[must_use]
    pub fn plate_force_ports(&self, blade: usize) -> &[PlateForcePort] {
        self.plate_force_schedule(blade).map_or(&[], |f| f.ports.as_slice())
    }
    /// Raw assignments committed with successfully ingested MECHANICAL samples,
    /// including any accepted prefix inside an unfinished, failed audio callback.
    #[must_use]
    pub fn applied_plate_force_controls(&self, blade: usize) -> &[ApertureForceEvent] {
        self.plate_force_schedule(blade).map_or(&[], |f| &f.events[..f.applied])
    }
    /// Uncommitted assignments in stable mechanical-time execution order.
    #[must_use]
    pub fn pending_plate_force_controls(&self, blade: usize) -> &[ApertureForceEvent] {
        self.plate_force_schedule(blade).map_or(&[], |f| &f.events[f.applied..])
    }
    /// Last accepted input forces [N]. Scalar/paired footprint programs supply
    /// generalized forces; material-modal programs retain physical port forces
    /// for the original basis projection. Failed candidates never publish loads.
    #[must_use]
    pub fn held_generalized_forces_n(&self) -> [f64; 2] {
        std::array::from_fn(|i| self.plate_forces[i].as_ref()
            .map_or(self.source.held.blade_forces_n[i], |f| f.held))
    }
}
