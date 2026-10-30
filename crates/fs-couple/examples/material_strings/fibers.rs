//! Data-only fiber geometry and material selection for the existing renderer.
use super::*;
pub use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::performance::fiber_motion::{
    FiberSupportKnot, FiberSupportProgram,
};
pub use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::ribbon::motion::FiberSupportMotionLimits;
pub use fs_couple::fiber_ribbon::strain::FiberStrainLimits;
use fs_couple::fiber_ribbon::{FiberRibbonGeometry, FiberRibbonResolution, TensionedFiberRibbon};
use fs_material::section::relaxation::{
    LayerRelaxationProperties, MAX_SECTION_RELAXATION_BRANCHES,
};
use fs_material::section::{
    ConcentricElasticSection, ConcentricMaterialLayer, axial::YOUNG_VISCOSITY_PROPERTY,
};
use fs_qty::Time;

/// Source roles for the fiber's isotropic Young relaxation spectrum. The same
/// law acts on axial extension and bending; these are not independent loss fits.
#[derive(Clone, Debug)]
pub struct Memory {
    pub properties: LayerRelaxationProperties,
    pub max_states: usize,
}

impl Memory {
    fn requirements(&self) -> Result<Vec<Requirement>> {
        if self.max_states == 0 || self.properties.branches.len() > MAX_SECTION_RELAXATION_BRANCHES
        {
            return Err(
                "fiber memory needs a positive history budget and bounded source arms".into(),
            );
        }
        let mut requirements = vec![Requirement::try_new(
            &self.properties.equilibrium_modulus,
            Pressure::DIMS,
            Admit::StrictlyPositive,
        )?];
        for (index, (modulus, time)) in self.properties.branches.iter().enumerate() {
            if self.properties.branches[..index].contains(&(modulus.clone(), time.clone())) {
                return Err("fiber memory repeats a source-property pair".into());
            }
            requirements.push(Requirement::try_new(
                modulus,
                Pressure::DIMS,
                Admit::NonNegative,
            )?);
            requirements.push(Requirement::try_new(
                time,
                Time::DIMS,
                Admit::StrictlyPositive,
            )?);
        }
        Ok(requirements)
    }
}

#[derive(Clone, Debug)]
pub struct Input {
    pub pack: String,
    pub radius_m: f64,
    pub length_m: f64,
    pub fiber_count: usize,
    pub tension_n: f64,
    pub station_fraction: f64,
    pub axial_modes: usize,
    pub transverse_modes: usize,
    pub memory: Option<Memory>,
    pub strain: Option<(FiberStrainLimits, NonlinearModalConfig)>,
    pub motion: Option<FiberSupportMotionLimits>,
    pub support_program: Option<FiberSupportProgram>,
    pub normal_program: Option<FiberSupportProgram>,
}
impl Input {
    pub fn mode_count(&self, spec: &input::Spec) -> Result<usize> {
        let (_, cells) = spec
            .bow_ribbon
            .ok_or("bow_fibers requires bow_ribbon width/cells")?;
        if self.axial_modes == 0
            || self.transverse_modes == 0
            || self.fiber_count == 0
            || !(0.0 < self.station_fraction && self.station_fraction < 1.0)
        {
            return Err(
                "bow_fibers requires positive span modes/fiber count and interior fiber station"
                    .into(),
            );
        }
        if let Some(memory) = &self.memory {
            memory.requirements()?;
            if self.strain.is_some() {
                return Err("fiber relaxation requires linear strain; nonlinear material-history transport is not available".into());
            }
        }
        if let Some((limits, config)) = self.strain {
            limits.validate()?;
            config.validate()?;
            if self.axial_modes > fs_nlmodal::longitudinal::MAX_LONGITUDINAL_MODES
                || self.transverse_modes > fs_nlmodal::longitudinal::MAX_TRANSVERSE_MODES
            {
                return Err(
                    "fiber strain requires an explicitly supported axial/normal resolution".into(),
                );
            }
            let string_config = spec
                .geometric_stretching
                .as_ref()
                .map(|(_, config)| config)
                .or_else(|| spec.longitudinal_waves.as_ref().map(|(_, config)| config));
            if string_config.is_some_and(|existing| *existing != config) {
                return Err(
                    "fiber and string strain must use the same joint solve configuration".into(),
                );
            }
        }
        if let Some(program) = &self.support_program {
            program.validate(spec.samples, spec.rate, &CancelGate::new_clock_free())?;
            if self.motion.is_some() {
                return Err("fiber motion and fiber program are exclusive support controls".into());
            }
        }
        if let Some(limits) = self.motion {
            limits.validate()?;
        }
        if let Some(program) = &self.normal_program {
            program.validate(spec.samples, spec.rate, &CancelGate::new_clock_free())?;
            if self.motion.is_none() && self.support_program.is_none() {
                return Err(
                    "normal fiber support requires an explicitly retained axial support frame"
                        .into(),
                );
            }
        }
        // An accelerating stroke must opt into source-derived inertia and
        // contact travel; omission retains the explicit frozen-station rung.
        for track in spec.schedule.tracks() {
            let GestureValue::Bow {
                velocity_m_per_s: initial,
                ..
            } = track.initial
            else {
                return Err("fiber ribbon needs a complete bow drive".into());
            };
            if self.support_program.is_some() && initial != 0.0 {
                return Err("fiber support program requires zero competing gesture speed".into());
            }
            for event in &track.events {
                if self.motion.is_none()
                    && !matches!(event.value, GestureValue::Bow { velocity_m_per_s, .. }
                    if velocity_m_per_s.to_bits() == initial.to_bits())
                {
                    return Err("frozen-station fiber ribbon requires constant support speed; load may vary".into());
                }
            }
        }
        self.axial_modes
            .checked_add(self.transverse_modes)
            .and_then(|n| n.checked_mul(cells))
            .ok_or_else(|| "fiber mode count overflow".into())
    }

    pub fn compile(
        &self,
        spec: &input::Spec,
        catalog: &mut Catalog,
        loader: &mut impl FnMut(&str) -> Result<MaterialPack>,
        gate: &CancelGate,
    ) -> Result<TensionedFiberRibbon> {
        self.mode_count(spec)?;
        poll(gate)?;
        // Profile ownership is available here, after parsing the performance
        // and before resolving any fiber material. Dimension counting alone
        // cannot reject a film whose separate thermal profile is not loaded yet.
        if (self.motion.is_some() || self.support_program.is_some())
            && ((spec.contact_thermal && !spec.dry_contact_thermal)
                || spec.friction.shear_layers().is_some())
            && !spec.retained_string_films
        {
            return Err("moving fibers require dry heat routes or a complete v3 film binding attached to the original string patches".into());
        }
        let mut requirements = vec![
            Requirement::try_new("density", Density::DIMS, Admit::StrictlyPositive)?,
            Requirement::try_new("young_modulus", Pressure::DIMS, Admit::StrictlyPositive)?,
            Requirement::try_new(ETA, DynViscosity::DIMS, Admit::NonNegative)?,
            Requirement::try_new(
                YOUNG_VISCOSITY_PROPERTY,
                DynViscosity::DIMS,
                Admit::NonNegative,
            )?,
        ];
        if let Some(memory) = &self.memory {
            for requirement in memory.requirements()? {
                // The equilibrium selector may explicitly name young_modulus.
                // Keep one requirement per identical source role, while refusing
                // an alias that assigns incompatible physical quantities.
                if let Some(existing) = requirements
                    .iter()
                    .find(|existing| existing.name() == requirement.name())
                {
                    if existing != &requirement {
                        return Err("fiber memory source roles conflict with the original material requirements".into());
                    }
                } else {
                    requirements.push(requirement);
                }
            }
        }
        let material = resolve_properties(&self.pack, catalog, loader, spec, &requirements)?;
        let section = ConcentricElasticSection::try_new(
            0.0,
            vec![ConcentricMaterialLayer {
                outer_radius_m: self.radius_m,
                material,
            }],
        )?;
        let (width_m, cells) = spec
            .bow_ribbon
            .ok_or("bow_fibers requires a physical ribbon")?;
        let memory = self
            .memory
            .as_ref()
            .map(|binding| {
                section
                    .clone()
                    .with_prony_bending(std::slice::from_ref(&binding.properties))
            })
            .transpose()?;
        let ribbon = TensionedFiberRibbon::try_new(
            FiberRibbonGeometry {
                length_m: self.length_m,
                width_m,
                fiber_count: self.fiber_count,
                tension_n: self.tension_n,
            },
            FiberRibbonResolution {
                cells,
                axial_modes: self.axial_modes,
                transverse_modes: self.transverse_modes,
            },
            section
                .clone()
                .with_axial_viscosity(YOUNG_VISCOSITY_PROPERTY)?,
            section.with_kelvin_voigt_bending(ETA)?,
            spec.rate,
            modal_budget(),
            gate,
        )?;
        let ribbon = match (&self.memory, memory) {
            (Some(binding), Some(memory)) => {
                ribbon.with_material_relaxation(memory, binding.max_states, gate)?
            }
            _ => ribbon,
        };
        match self.strain {
            Some((limits, config)) => Ok(ribbon.with_axial_transverse_strain(limits, config)?),
            None => Ok(ribbon),
        }
    }
}
