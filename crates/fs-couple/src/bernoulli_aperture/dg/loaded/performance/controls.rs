//! Reuse the original source-mesh force and pressure compilers.
use super::{ForceSchedule, GestureCompileError, LoadedDgPerformance, PressureSchedule, invalid};
use crate::bernoulli_aperture::dynamic::force::PlateForcePort;
use crate::bernoulli_aperture::performance::force::{ApertureForceEvent, ApertureForceProgram};
use crate::bernoulli_aperture::two_sided::LoadedAperture;
use crate::bernoulli_aperture::two_sided::observer::schedule::spatial_pressures::MaterialPressureEvent;
use crate::render::RenderError;

impl LoadedDgPerformance {
    fn admit_force_attachment(&self, blade: usize) -> Result<(), GestureCompileError> {
        if blade >= 2 || self.poisoned || self.system.aperture().accepted_steps() != 0
            || self.plate_forces[blade].is_some()
        {
            return Err(invalid("loaded DG force controls require one cold attachment per blade"));
        }
        Ok(())
    }

    /// Original material-modal physical force footprint, in newtons. The retained
    /// basis projects it onto every admitted mode during the same coupled solve.
    /// Events address port zero on this blade and use mechanical sample numbers.
    /// Equal-time order and release semantics belong to the shared force compiler.
    pub fn with_material_plate_force_controls(
        mut self, blade: usize, events: Vec<ApertureForceEvent>,
    ) -> Result<Self, GestureCompileError> {
        self.admit_force_attachment(blade)?;
        let LoadedAperture::Modal(source) = self.system.aperture() else {
            return Err(invalid("source-defined force controls require original material modal plates"));
        };
        if source.material_plates().get(blade).is_none() {
            return Err(invalid("source force names an absent material plate"));
        }
        let mut probe = [0.0; 2];
        probe[blade] = 1.0;
        source.validate_forces(probe)
            .map_err(|error| GestureCompileError::Render(RenderError::Voice(error)))?;
        let (assignments, work) = self.remaining_control_budgets()?;
        self.plate_forces[blade] = Some(ForceSchedule::source_port(
            events, self.config.samples, assignments, work,
        )?);
        Ok(self)
    }

    /// Force in newtons on the original complete solid nodal footprint. Every
    /// selected mode, including a zero-opening mode, uses the same source basis
    /// as pressure and swept volume. The existing source-port scheduler owns
    /// equal-time assignments, release and transactional accepted-sample state.
    pub fn with_material_solid_force_controls(
        mut self, events: Vec<ApertureForceEvent>,
    ) -> Result<Self, GestureCompileError> {
        self.admit_force_attachment(0)?;
        let LoadedAperture::Modal(source) = self.system.aperture() else {
            return Err(invalid("solid force controls require an original material solid"));
        };
        if source.material_solid().is_none() && source.hereditary_solid().is_none() {
            return Err(invalid("solid force controls require an original material solid"));
        }
        source.validate_forces([1.0, 0.0])
            .map_err(|error| GestureCompileError::Render(RenderError::Voice(error)))?;
        let (assignments, work) = self.remaining_control_budgets()?;
        self.plate_forces[0] = Some(ForceSchedule::source_port(
            events, self.config.samples, assignments, work,
        )?);
        Ok(self)
    }

    /// Project explicit physical node/patch loads from a retained scalar or paired
    /// source plate. Material-modal sources instead use their complete force port
    /// through with_material_plate_force_controls; no modal column is discarded.
    pub fn with_plate_forces(
        mut self, blade: usize, program: ApertureForceProgram,
    ) -> Result<Self, GestureCompileError> {
        self.admit_force_attachment(blade)?;
        let plate = match self.system.aperture() {
            LoadedAperture::Single(a) if blade == 0 => a.plate_reduction(),
            LoadedAperture::Paired(a) => a.plates().map(|p| &p[blade]),
            _ => None,
        }.ok_or_else(|| invalid("physical force footprints require this blade's retained source plate"))?;
        let (assignments, work) = self.remaining_control_budgets()?;
        self.plate_forces[blade] = Some(ForceSchedule::compile(
            plate, program, self.config.samples, assignments, work,
        )?);
        Ok(self)
    }

    /// Signed pressure amplitude [Pa] on one previously bound material pressure
    /// trace. Its source-mesh shape and swept-volume projection stay in the valve.
    /// This prescribed source is distinct from blowing pressure and physical force;
    /// its work is already counted once by the original mechanical owner.
    pub fn with_material_pressure_controls(
        mut self, blade: usize, events: Vec<MaterialPressureEvent>,
    ) -> Result<Self, GestureCompileError> {
        if blade >= 2 || self.poisoned || self.system.aperture().accepted_steps() != 0
            || self.spatial_pressures[blade].is_some()
        {
            return Err(invalid("loaded DG pressure controls require one cold attachment per material port"));
        }
        let mut probe = [0.0; 2];
        probe[blade] = 1.0;
        self.system.aperture().validate_pressures(probe)
            .map_err(|error| GestureCompileError::Render(RenderError::Voice(error)))?;
        let (assignments, work) = self.remaining_control_budgets()?;
        self.spatial_pressures[blade] = Some(PressureSchedule::compile(
            events, self.config.samples, assignments, work,
        )?);
        Ok(self)
    }

    /// Source-bound work-conjugate ports, empty for material-modal source footprints.
    #[must_use]
    pub fn plate_force_ports(&self, blade: usize) -> &[PlateForcePort] {
        self.plate_forces.get(blade).and_then(Option::as_ref)
            .map_or(&[], |f| f.ports.as_slice())
    }
    /// Raw physical force assignments whose source sample has been accepted.
    #[must_use]
    pub fn applied_plate_force_controls(&self, blade: usize) -> &[ApertureForceEvent] {
        self.plate_forces.get(blade).and_then(Option::as_ref)
            .map_or(&[], |f| &f.events[..f.applied])
    }
    /// Raw physical force assignments still pending.
    #[must_use]
    pub fn pending_plate_force_controls(&self, blade: usize) -> &[ApertureForceEvent] {
        self.plate_forces.get(blade).and_then(Option::as_ref)
            .map_or(&[], |f| &f.events[f.applied..])
    }
    /// Last accepted force amplitudes, never ahead of the source histories.
    #[must_use]
    pub fn held_generalized_forces_n(&self) -> [f64; 2] {
        std::array::from_fn(|i| self.plate_forces[i].as_ref().map_or(0.0, |f| f.held))
    }
    /// Raw pressure assignments whose source sample has been accepted.
    #[must_use]
    pub fn applied_material_pressure_controls(&self, blade: usize) -> &[MaterialPressureEvent] {
        self.spatial_pressures.get(blade).and_then(Option::as_ref)
            .map_or(&[], |p| &p.events[..p.next])
    }
    /// Raw source-mesh pressure assignments still pending.
    #[must_use]
    pub fn pending_material_pressure_controls(&self, blade: usize) -> &[MaterialPressureEvent] {
        self.spatial_pressures.get(blade).and_then(Option::as_ref)
            .map_or(&[], |p| &p.events[p.next..])
    }
    /// Last accepted source-mesh pressure amplitudes in pascals.
    #[must_use]
    pub fn held_material_pressures_pa(&self) -> [f64; 2] {
        std::array::from_fn(|i| self.spatial_pressures[i].as_ref().map_or(0.0, |p| p.held))
    }
}
