//! Material history reacts through the original shared body and physical bridges.
use super::*;
use crate::plate_viscosity::relaxation::{
    InitialPlateMemory, ModalPlateRelaxation, PlateBendingRelaxation, PlateRelaxationLimits,
};
use crate::plate_viscosity::relaxation::tensor::{TensorPlateRelaxation,TensorPlateRelaxationLimits};
use fs_material::visco_tensor::source::TensorPronyProperties;

impl MaterialStringBankSource {
    /// All original regional spectra and their full same-basis strain factors.
    /// This image survives gas/torque rebuilding and multibow/PCM ownership.
    #[must_use]
    pub const fn material_chart_relaxation(&self) -> Option<&ModalPlateRelaxation> {
        self.material_chart_relaxation.as_ref()
    }

    /// Complete source tensors and every original normal/shear history.
    pub const fn material_chart_tensor_relaxation(&self) -> Option<&TensorPlateRelaxation> {
        self.material_chart_tensor_relaxation.as_ref()
    }

    pub(in crate::bowed_string::runtime::reciprocal::bank) fn validate_body_memory_frequencies(
        &self, system: &FrictionModalSystem, gate: &CancelGate,
    ) -> Result<(), ReciprocalStringError> {
        poll(gate)?;
        self.validate_source_memory_frequencies(system, gate)?;
        if let Some(image) = &self.material_chart_relaxation {
            // A coupled body receives forcing and relaxation below its bare
            // resonances as well as at the instantaneous coupled frequencies.
            // Validate all source coefficients, even explicitly zero arms.
            image.source().validate_band(0.0, system.linear_frequency_upper_bound_rad_s(),
                1.0 / f64::from(self.options.sample_rate_hz))
                .map_err(ReciprocalStringError::Material)?;
        }
        if let Some(image) = &self.string_bending_relaxation {
            image.validate_band(system.linear_frequency_upper_bound_rad_s(),
                1.0 / f64::from(self.options.sample_rate_hz))?;
        }
        if let Some(image)=&self.material_chart_tensor_relaxation {
            image.validate_band(system.linear_frequency_upper_bound_rad_s(),
                1.0/f64::from(self.options.sample_rate_hz)).map_err(ReciprocalStringError::Material)?;
        }
        Ok(())
    }

    // Re-embed the retained body factors without resolving material or modes
    // again. The common owner orders source strings, body, then added layers.
    pub(in crate::bowed_string::runtime::reciprocal::bank) fn body_memory_branches(
        &self, system: &FrictionModalSystem, gate: &CancelGate,
    ) -> Result<Vec<fs_phs::RelaxationBranch>, ReciprocalStringError> {
        poll(gate)?;
        let mut branches = Vec::new();
        let image=self.material_chart_relaxation.as_ref().map(|image|
            (image.mode_count(),image.coupling_rows_per_s(),image.relaxation_times_s()))
            .or_else(||self.material_chart_tensor_relaxation.as_ref().map(|image|
                (image.mode_count(),image.coupling_rows_per_s(),image.relaxation_times_s())));
        if let Some((modes,rows,times)) = image {
            let body = self.strings.len();
            let component = system.components().get(body)
                .ok_or_else(|| invalid("regional memory lost its original body component"))?;
            if component.modes().len() != modes
                || component.modes().iter().any(|m| m.damping_ratio != 0.0) {
                return Err(invalid("regional body memory requires the unchanged complete elastic body basis"));
            }
            let start: usize = system.components()[..body].iter().map(|c| c.modes().len()).sum();
            let dimension = system.mode_count().checked_mul(2)
                .ok_or_else(|| invalid("body memory state embedding overflow"))?;
            branches.extend(rows.iter().zip(times)
                .map(|(row, &time)| {
                    let mut projection = vec![0.0; dimension];
                    for (i, &value) in row.iter().enumerate() { projection[2*(start+i)] = value; }
                    fs_phs::RelaxationBranch { projection, stiffness: 1.0, relaxation_time_s: time }
                }));
        }
        poll(gate)?;
        Ok(branches)
    }

    pub(in crate::bowed_string::runtime::reciprocal::bank) fn restore_body_tensor_constraints(
        &self,system:FrictionModalSystem,cards:&[BowedStringCard],gate:&CancelGate,
    ) -> Result<FrictionModalSystem,ReciprocalStringError> {
        let Some(image)=&self.material_chart_tensor_relaxation else {return Ok(system);};
        let start:usize=system.components()[..self.strings.len()].iter().map(|c|c.modes().len()).sum();
        let memory=self.source_memory_count(cards,gate)?;
        system.with_plane_stress_material_points(image.bindings(start,memory),gate)
            .map_err(ReciprocalStringError::Friction)
    }
}

impl ReciprocalStringBank {
    /// Resolve independent regional material spectra into the ORIGINAL body's
    /// retained DKT basis and attach their internal strains to the SAME coupled
    /// string/body/cavity/friction owner. Material changes therefore modify both
    /// physical pressure and the bridge reactions back onto every string.
    /// Original string Prony histories remain unchanged; new regional body
    /// histories follow them and precede any explicit layered-string addition.
    ///
    /// Each region's selected excess modulus scales its actual rotated bending
    /// tensor; all modal cross terms survive. The equilibrium elastic pencil,
    /// conservative prestress, nodal bridge footprints, gas and observer remain
    /// in that one basis. This is not a general tensor-Prony plane-stress law.
    ///
    /// The existing dense exact capacity and setup allowance include every real
    /// first-order memory state. Prepared mechanical-only propagation refuses.
    /// Source validity conservatively covers DC through the complete coupled
    /// instantaneous stiffness bound. This does not qualify finite-amplitude
    /// harmonic content or reproduce a measured violin/cello/bass by itself.
    ///
    /// # Errors
    /// Requires one cold attachment to an original chart with no selected scalar
    /// or tensor viscosity or authored body decay. Missing source properties,
    /// unresolved positive forms, excessive state/projection work, source/clock
    /// refusal and cancellation return no partially assembled performer.
    pub fn with_material_plate_relaxation(
        mut self, properties: &[(String, String)], initial: InitialPlateMemory,
        limits: PlateRelaxationLimits, gate: &CancelGate,
    ) -> Result<Self, ReciprocalStringError> {
        poll(gate)?;
        let source_count = self.source.source_memory_count(&self.cards, gate)?;
        if self.system.samples_rendered() != 0 || self.source.material_chart_relaxation.is_some()
            || self.source.material_chart_tensor_relaxation.is_some()
            || self.source.material_chart_viscosity.is_some()
            || self.source.string_bending_relaxation.is_some()
            || self.system.material_relaxation_branches().len() != source_count {
            return Err(invalid("body material memory requires one cold attachment without parallel intrinsic damping"));
        }
        let body = self.source.material_chart.as_ref()
            .ok_or_else(|| invalid("body material memory requires an original regional plate chart"))?;
        let basis = self.source.material_chart_basis.as_ref()
            .ok_or_else(|| invalid("body material memory requires the original full DKT basis"))?;
        if body.reduction.damping_ratio != 0.0 || body.specimen.bending_viscosity().is_some() {
            return Err(invalid("body material memory cannot duplicate authored or selected intrinsic body loss"));
        }
        let law = PlateBendingRelaxation::from_chart(&body.specimen, properties, gate)
            .map_err(ReciprocalStringError::Material)?;
        let image = law.project(&body.specimen, basis, limits, gate)
            .map_err(ReciprocalStringError::Material)?;
        let memory = match initial {
            InitialPlateMemory::Relaxed => None,
            InitialPlateMemory::Unrelaxed => Some(vec![0.0; image.coupling_rows_per_s().len()]),
        };
        self.source.material_chart_relaxation = Some(image);
        let branches = self.source.body_memory_branches(&self.system, gate)?;
        self.system = self.system.with_additional_material_relaxation(branches, memory, gate)
            .map_err(ReciprocalStringError::Friction)?;
        self.source.validate_body_memory_frequencies(&self.system, gate)?;
        self.source.validate_material_bridge_frequencies(&self.system, gate)?;
        poll(gate)?;
        Ok(self)
    }

    /// Attach complete source-derived 3D tensor memory to the original body.
    /// The shared DKT quadrature and plane-stress constraint owners retain all
    /// six histories per branch/point, including uncoupled histories. The exact
    /// original string/body/bridge/friction solve carries their force and heat.
    /// Surface strain is admitted atomically at each physical endpoint, through
    /// static preload, nonlinear contacts, thermal feedback and pressure output.
    /// This fixed-state homogeneous plate model does not infer shell, wet wood,
    /// finite-strain or temperature evolution; original source/capacity limits apply.
    pub fn with_tensor_material_plate_relaxation(self,
        properties:&[TensorPronyProperties], initial:InitialPlateMemory,
        limits:TensorPlateRelaxationLimits, gate:&CancelGate,
    ) -> Result<Self,ReciprocalStringError> {
        self.with_tensor_plate_memory(properties,initial,limits,None,gate)
    }

    /// Share only bit-identically timed material coordinates through an exact
    /// orthogonal change of variables. Every original point history and source
    /// strain remains reconstructible; the original mechanical basis is unchanged.
    /// This admits relaxed or unrelaxed initial histories, not arbitrary
    /// independently authored point histories. No mode or relaxation rate is
    /// removed or merged to satisfy capacity, and all exact-owner caps remain.
    pub fn with_shared_tensor_material_plate_relaxation(self,
        properties:&[TensorPronyProperties], initial:InitialPlateMemory,
        limits:TensorPlateRelaxationLimits, max_factor_terms:usize, gate:&CancelGate,
    ) -> Result<Self,ReciprocalStringError> {
        self.with_tensor_plate_memory(properties,initial,limits,Some(max_factor_terms),gate)
    }

    fn with_tensor_plate_memory(mut self,properties:&[TensorPronyProperties],
        initial:InitialPlateMemory,limits:TensorPlateRelaxationLimits,
        sharing:Option<usize>,gate:&CancelGate,
    )->Result<Self,ReciprocalStringError> {
        poll(gate)?;
        let source_count=self.source.source_memory_count(&self.cards,gate)?;
        if self.system.samples_rendered()!=0 || self.source.material_chart_relaxation.is_some()
            || self.source.material_chart_tensor_relaxation.is_some()
            || self.source.material_chart_viscosity.is_some()
            || self.source.string_bending_relaxation.is_some()
            || self.system.material_relaxation_branches().len()!=source_count {
            return Err(invalid("complete tensor body memory requires one cold binding without duplicate intrinsic damping"));
        }
        let body=self.source.material_chart.as_ref()
            .ok_or_else(||invalid("complete tensor body memory requires its original regional plate chart"))?;
        let basis=self.source.material_chart_basis.as_ref()
            .ok_or_else(||invalid("complete tensor body memory requires its original full DKT basis"))?;
        if body.reduction.damping_ratio!=0.0 || body.specimen.bending_viscosity().is_some() {
            return Err(invalid("complete tensor body memory cannot duplicate authored or selected intrinsic decay"));
        }
        let image=match sharing {
            Some(work)=>TensorPlateRelaxation::from_chart_shared(
                &body.specimen,basis,properties,limits,work,gate),
            None=>TensorPlateRelaxation::from_chart(&body.specimen,basis,properties,limits,gate),
        }.map_err(ReciprocalStringError::Material)?;
        let memory=match initial {
            InitialPlateMemory::Relaxed=>None,
            InitialPlateMemory::Unrelaxed=>Some(vec![0.0;image.coupling_rows_per_s().len()]),
        };
        self.source.material_chart_tensor_relaxation=Some(image);
        let branches=self.source.body_memory_branches(&self.system,gate)?;
        self.system=self.system.with_additional_material_relaxation(branches,memory,gate)
            .map_err(ReciprocalStringError::Friction)?;
        self.system=self.source.restore_body_tensor_constraints(self.system,&self.cards,gate)?;
        self.source.validate_body_memory_frequencies(&self.system,gate)?;
        self.source.validate_material_bridge_frequencies(&self.system,gate)?;
        poll(gate)?;
        Ok(self)
    }
}
