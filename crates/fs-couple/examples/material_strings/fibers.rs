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
    thermal::ConcentricThermoelasticSection,
};
use fs_couple::material_assembly::concentric::{
    with_concentric_section, with_concentric_thermal_extension,
};
use fs_qty::Time;
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::compliant_ribbon::HairRibbonProperties;

/// Explicit low-inertia axial reduction of parallel fibers. The selected
/// source point remains distinct from the simulation clock: choosing a static
/// modulus is a declared constitutive continuation, not broadband measurement.
/// Density remains in the retained source inputs for the resolved-wave comparison;
/// the spring-dashpot owner does not introduce a fictitious lumped fiber mass.
#[derive(Clone, Debug)]
pub struct QuasistaticInput {
    pub pack: String,
    pub radius_m: f64,
    pub length_m: f64,
    pub fiber_count: usize,
    pub station_fraction: f64,
    pub axes: Vec<(String, f64)>,
    pub young_property: String,
    pub density_property: String,
    /// None is only the explicitly selected lossless reference. A missing
    /// property in the viscosity form is an error, never zero damping.
    pub viscosity_property: Option<String>,
}

impl QuasistaticInput {
    fn requirements(&self) -> Result<Vec<Requirement>> {
        let mut result = vec![
            Requirement::try_new(&self.young_property, Pressure::DIMS, Admit::StrictlyPositive)?,
            Requirement::try_new(&self.density_property, Density::DIMS, Admit::StrictlyPositive)?,
        ];
        if let Some(name) = &self.viscosity_property {
            result.push(Requirement::try_new(name, DynViscosity::DIMS, Admit::NonNegative)?);
        }
        for (index, item) in result.iter().enumerate() {
            if result[..index].iter().any(|previous| previous.name() == item.name()) {
                return Err("quasistatic fiber source roles must name distinct physical properties".into());
            }
        }
        Ok(result)
    }

    fn source_axes(&self, spec: &input::Spec) -> Result<Vec<(String, f64)>> {
        let mut axes = self.axes.clone();
        axes.extend(spec.material_states.get(&self.pack).into_iter().flatten().cloned());
        if axes.len() > fs_couple::material_source::MAX_RECORD_AXES
            || axes.iter().enumerate().any(|(i, (name, value))|
                name.is_empty() || name.len() > 128 || !value.is_finite()
                    || axes[..i].iter().any(|(other, _)| other == name))
        {
            return Err("quasistatic fiber coordinates require at most eight distinct finite source axes".into());
        }
        Ok(axes)
    }

    pub fn preflight(&self, spec: &input::Spec) -> Result<()> {
        let (_, cells) = spec.bow_ribbon.ok_or("bow_compliance requires a finite bow_ribbon")?;
        if spec.bow_fibers.is_some() {
            return Err("choose quasistatic bow_compliance or resolved bow_fibers, not both".into());
        }
        if !self.radius_m.is_finite() || self.radius_m <= 0.0
            || !self.length_m.is_finite() || self.length_m <= 0.0
            || self.fiber_count == 0 || self.fiber_count > (1 << 24)
            || !self.station_fraction.is_finite()
            || !(0.0 < self.station_fraction && self.station_fraction < 1.0)
        {
            return Err("quasistatic fibers require positive geometry/count and an interior support station".into());
        }
        self.requirements()?;
        self.source_axes(spec)?;
        for track in spec.schedule.tracks() {
            let string = fs_scenario::gesture::path::legacy_member(&track.target, "bow")
                .ok_or("quasistatic ribbon requires an original bow track")?;
            if spec.strings.get(string).is_none_or(|wire| wire.bending_modes < 4 * cells) {
                return Err("quasistatic ribbon requires at least four string modes per contact cell".into());
            }
        }
        Ok(())
    }

    pub fn compile(
        &self,
        spec: &input::Spec,
        catalog: &mut Catalog,
        loader: &mut impl FnMut(&str) -> Result<MaterialPack>,
        gate: &CancelGate,
    ) -> Result<HairRibbonProperties> {
        self.preflight(spec)?;
        poll(gate)?;
        load_pack(&self.pack, catalog, loader)?;
        let card = &catalog[&self.pack].card;
        let query = fs_couple::material_source::MaterialRecord {
            pack: self.pack.clone(), state: card.id().process.clone(), axes: self.source_axes(spec)?,
        }.query_point(card)?;
        let source = resolve_material_state_point(card, &query, &self.requirements()?, Pick::SingleClaimOnly)?;
        let value = |name: &str| -> Result<f64> {
            Ok(source.property(name).ok_or("missing resolved fiber property")?.value_si())
        };
        let youngs_modulus_pa = value(&self.young_property)?;
        let density_kg_m3 = value(&self.density_property)?;
        let viscosity = self.viscosity_property.as_ref().map(|name| value(name)).transpose()?.unwrap_or(0.0);
        // Both fixed-end branches act in parallel at x = gamma L. The same
        // geometric factor converts E to S and material viscosity to D; width
        // and numerical cell count never change either total.
        let area = self.fiber_count as f64 * core::f64::consts::PI * self.radius_m.powi(2);
        let factor = area / (self.length_m * self.station_fraction * (1.0 - self.station_fraction));
        let damping = factor * viscosity;
        let stiffness = factor * youngs_modulus_pa;
        let (width_m, _) = spec.bow_ribbon.ok_or("missing ribbon geometry")?;
        if !stiffness.is_finite() || stiffness <= 0.0 || !damping.is_finite()
            || (viscosity > 0.0 && damping == 0.0)
        {
            return Err("source-derived quasistatic fiber response is unrepresentable".into());
        }
        poll(gate)?;
        Ok(HairRibbonProperties {
            youngs_modulus_pa, hair_radius_m: self.radius_m,
            hairs_in_contact: self.fiber_count as f64, free_length_m: self.length_m,
            excitation_fraction: self.station_fraction, density_kg_m3,
            distributed_damping_kg_s_m: damping / width_m, reference_width_m: width_m,
        })
    }
}

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
    /// FixedTension is the TOTAL bundle force. FixedExtension applies the same
    /// support spacing to each original fiber before summing their forces.
    pub prestress: StringPrestress,
    /// Stress-free reference temperature for the explicit thermal-extension form.
    pub thermal_reference_kelvin: Option<f64>,
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
        self.validate_support()?;
        if let Some(memory) = &self.memory {
            memory.requirements()?;
            if !matches!(self.prestress, StringPrestress::FixedTension(_)) {
                return Err("fixed-extension fiber preload requires an initialized hereditary strain history; use an explicit maintained force with fiber memory".into());
            }
            // The joint material-history owner admits the selected nonlinear
            // strain. Its geometric mean adds one state per original arm and
            // cell; the authored budget is checked when the section resolves.
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

    /// Admit the authored boundary condition before any material-pack I/O.
    /// Thermal expansion may make an initially shorter support spacing taut;
    /// the source-bound section, not this syntax check, decides its resultant.
    fn validate_support(&self) -> Result<()> {
        if !self.length_m.is_finite() || self.length_m <= 0.0 {
            return Err("fiber support spacing must be finite and positive".into());
        }
        match self.prestress {
            StringPrestress::FixedTension(tension) => {
                if !tension.is_finite() || tension <= 0.0 || self.thermal_reference_kelvin.is_some() {
                    return Err("maintained fiber tension must be positive and cannot also prescribe thermal extension".into());
                }
            }
            StringPrestress::FixedExtension { stress_free_length_m, linear_strain_limit } => {
                if !stress_free_length_m.is_finite() || stress_free_length_m <= 0.0
                    || !linear_strain_limit.is_finite()
                    || !(0.0 < linear_strain_limit && linear_strain_limit < 1.0)
                {
                    return Err("fiber extension needs a positive stress-free length and a strain limit in (0,1)".into());
                }
                let strain = (self.length_m - stress_free_length_m) / stress_free_length_m;
                if !strain.is_finite() || strain.abs() > linear_strain_limit
                    || (self.thermal_reference_kelvin.is_none() && strain <= 0.0)
                    || self.thermal_reference_kelvin.is_some_and(|t| !t.is_finite() || t <= 0.0)
                {
                    return Err("fiber support extension or reference temperature is inadmissible".into());
                }
            }
            _ => return Err("fiber supports require maintained total tension or explicit fixed extension".into()),
        }
        Ok(())
    }

    // The original section/support owner computes ONE physical fiber's force.
    // Numerical ribbon cells do not enter this boundary condition. The model
    // later partitions the original bundle force, mass and stiffness together.
    fn support_tension_n(&self, section: &ConcentricElasticSection, catalog: &Catalog) -> Result<f64> {
        self.validate_support()?;
        if let StringPrestress::FixedTension(tension) = self.prestress {
            return Ok(tension);
        }
        let input = PrestressedString {
            length_m: self.length_m, tension_n: 0.0, lin_density_kg_m: 1.0,
            axial_stiffness_n: 1.0, width_m: 1.0, n_modes: self.transverse_modes,
            damping_ratio: 0.0, rayleigh: None, bending_stiffness_n_m2: 0.0,
            kelvin_voigt_bending: None, relaxation_bending: None,
            polarization_detune: 0.0, moving_end: false, tension_loss_factor: None, bending_loss_factor: None,
        };
        let fiber = if let Some(reference) = self.thermal_reference_kelvin {
            let StringPrestress::FixedExtension { stress_free_length_m, linear_strain_limit } = self.prestress else {
                return Err("thermal fiber supports require a fixed-extension prescription".into());
            };
            let path = construction::expansion(&self.pack, &section.layers()[0].material, catalog, reference)?;
            let thermal = ConcentricThermoelasticSection::try_new(section.clone(), vec![path])?;
            with_concentric_thermal_extension(input, &thermal, stress_free_length_m, linear_strain_limit)?.0
        } else {
            with_concentric_section(input, section, self.prestress)?
        };
        let total = fiber.tension_n * self.fiber_count as f64;
        if !total.is_finite() || total <= 0.0 {
            return Err("source-derived total fiber tension is unrepresentable".into());
        }
        Ok(total)
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
        let tension_n = self.support_tension_n(&section, catalog)?;
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
                tension_n,
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
