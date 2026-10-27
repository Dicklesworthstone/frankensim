//! Data-only construction over the original specimen, concentric-section and
//! torsion owners. No section formula, support law or material fallback lives here.
use super::*;
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::section::thermal::ConcentricThermoelasticSection;
use fs_material::state_point::{IntegratedIsotropicThermalExpansion, integrate_isotropic_thermal_expansion};
use fs_couple::material_assembly::concentric::with_concentric_thermal_extension;
use fs_couple::string_specimen::with_uniform_circular_thermal_extension;

// Use the exact material card and current coordinates of each physical layer.
// The shared owner integrates alpha(T) across every knot and checks the whole
// path; an endpoint coefficient times delta-T is not a thermal support law.
fn expansion(name: &str, material: &ResolvedMaterialStatePoint, catalog: &Catalog,
    reference: f64) -> Result<IntegratedIsotropicThermalExpansion>
{
    let mut current = QueryPoint::new();
    for (axis, value) in material.query_point() {
        current = match material.axis_quantities().get(axis) {
            Some(quantity) => current.with_quantity(axis.clone(), *quantity, *value)?,
            None => current.with(axis.clone(), *value)?,
        };
    }
    let point = current.clone().with("T", reference)?;
    Ok(integrate_isotropic_thermal_expansion(&catalog[name].card,
        &point, &current, Pick::SingleClaimOnly)?)
}

impl input::Wire {
    pub(super) fn resolve_source(&self, catalog: &mut Catalog,
        loader: &mut impl FnMut(&str) -> Result<NormalizedPack>, spec: &input::Spec,
        gate: &CancelGate) -> Result<MaterialStringSource>
    {
        poll(gate)?;
        let thermal = match self.thermal_reference_kelvin {
            Some(reference) => {
                let StringPrestress::FixedExtension { stress_free_length_m, linear_strain_limit } = self.prestress else {
                    return Err("thermal support requires a fixed-extension prescription".into());
                };
                if !reference.is_finite() || reference <= 0.0 {
                    return Err("thermal support requires positive reference temperature".into());
                }
                Some((reference, stress_free_length_m, linear_strain_limit))
            },
            None => None,
        };
        let twist = self.twist_modes > 0;
        let core = resolve(&self.pack, catalog, loader, spec, false, twist, &self.bending_loss)?;
        let mut source = match &self.section {
            input::WireSection::CloseWound { core_radius_m, pack, wire_radius_m, pitch_m } => {
                if thermal.is_some() || twist || spec.longitudinal_waves.is_some() || self.torsion_memory.is_some() {
                    return Err("close-wound reduction supports frozen core flexure and spring extension; resolved winding thermal/torsional/axial-wave laws are unavailable".into());
                }
                let elastic = ConcentricElasticSection::try_new(0.0, vec![ConcentricMaterialLayer {
                    outer_radius_m: *core_radius_m, material: core,
                }])?;
                let core: fs_material::section::winding::WoundCoreSection = match &self.bending_loss {
                    input::BendingLoss::KelvinVoigt => elastic.with_kelvin_voigt_bending(ETA)?.into(),
                    input::BendingLoss::Prony { properties } => elastic.with_prony_bending(&[
                        fs_material::section::relaxation::LayerBendingRelaxationProperties {
                            equilibrium_modulus: fs_couple::string_specimen::EQUILIBRIUM_YOUNG_MODULUS_PROPERTY.into(),
                            branches: properties.iter().map(|p| (p.modulus.clone(), p.relaxation_time.clone())).collect(),
                        },
                    ])?.into(),
                };
                // The wrap's spring law needs density and G, not an invented
                // flexural viscosity or the core's Young modulus.
                let wrap = resolve_properties(pack, catalog, loader, spec, &[
                    Requirement::try_new("density", Density::DIMS, Admit::StrictlyPositive)?,
                    Requirement::try_new(SHEAR, Pressure::DIMS, Admit::StrictlyPositive)?,
                ])?;
                let section = fs_material::section::winding::CloseWoundSection::try_new(core,
                    fs_material::section::winding::RoundWireWinding {
                        wire_radius_m: *wire_radius_m, pitch_m: *pitch_m,
                    }, wrap, SHEAR)?;
                if (section.outer_radius_m() - self.radius).abs() > 8.0 * f64::EPSILON * self.radius {
                    return Err("winding geometry disagrees with the original finished radius".into());
                }
                MaterialStringSource::Wound(fs_couple::material_assembly::wound::ResolvedWoundStringSpecimen::try_new(
                    template(self), section, self.prestress)?)
            }
            input::WireSection::Homogeneous => {
                // V1 and explicit V2 fixed-tension specimens take precisely the
                // old route. The extension variant binds the original support law.
                let specimen = if let Some((reference, length, strain_limit)) = thermal {
                    let expansion = expansion(&self.pack, &core, catalog, reference)?;
                    with_uniform_circular_thermal_extension(template(self),
                        StringGeometryConstraint::FixedRadius(self.radius), &core, &expansion,
                        length, strain_limit)?
                } else { match self.prestress {
                    StringPrestress::FixedTension(_) =>
                        with_uniform_circular_material_state(template(self), self.radius, &core)?,
                    _ => with_uniform_circular_material_and_constraints(template(self),
                        StringGeometryConstraint::FixedRadius(self.radius), &core, self.prestress)?,
                } };
                let specimen = match &self.bending_loss {
                    input::BendingLoss::KelvinVoigt => specimen.with_kelvin_voigt_bending_loss()?,
                    input::BendingLoss::Prony { properties } => specimen.with_prony_bending_loss(properties)?,
                };
                MaterialStringSource::Homogeneous(specimen)
            }
            input::WireSection::Bonded { core_radius_m, sleeves } => {
                let mut layers = Vec::with_capacity(sleeves.len() + 1);
                layers.push(ConcentricMaterialLayer { outer_radius_m: *core_radius_m, material: core });
                for sleeve in sleeves {
                    poll(gate)?;
                    layers.push(ConcentricMaterialLayer { outer_radius_m: sleeve.outer_radius_m,
                        material: resolve(&sleeve.pack, catalog, loader, spec, false, twist, &self.bending_loss)? });
                }
                let elastic = ConcentricElasticSection::try_new(0.0, layers)?;
                if elastic.outer_radius_m().to_bits() != self.radius.to_bits() {
                    return Err("resolved section disagrees with the original string radius".into());
                }
                let prestress = if let Some((reference, length, strain_limit)) = thermal {
                    let names = std::iter::once(self.pack.as_str())
                        .chain(sleeves.iter().map(|layer| layer.pack.as_str()));
                    let mut paths = Vec::with_capacity(elastic.layers().len());
                    for (layer, name) in elastic.layers().iter().zip(names) {
                        poll(gate)?;
                        paths.push(expansion(name, &layer.material, catalog, reference)?);
                    }
                    let section = ConcentricThermoelasticSection::try_new(elastic.clone(), paths)?;
                    // The shared support owner checks EVERY layer's total,
                    // thermal and elastic strain, and refuses a slack resultant.
                    // Opposing layer stresses must not be hidden by averaging.
                    let (_, response) = with_concentric_thermal_extension(template(self),
                        &section, length, strain_limit)?;
                    if section.expansion_paths().iter().all(|path| path.free_linear_strain() == 0.0) {
                        // Preserve ordinary fixed-extension arithmetic exactly,
                        // but only AFTER validating all explicit thermal sources.
                        self.prestress
                    } else {
                        StringPrestress::FixedTension(response.force_n())
                    }
                } else { self.prestress };
                // Frozen-state lowering: the original reference prescription and
                // complete expansion cards remain in Ready.spec / Ready.catalog.
                // All later mechanical rebuilds use this same resultant, E/rho,
                // layer geometry and material loss. No live heat evolution is implied.
                match &self.bending_loss {
                    input::BendingLoss::KelvinVoigt =>
                        MaterialStringSource::Concentric { template: template(self),
                            section: elastic.with_kelvin_voigt_bending(ETA)?, prestress },
                    input::BendingLoss::Prony { properties } => {
                        // The same declared property roles are resolved in EACH
                        // original layer. Values and time constants remain local;
                        // no viscosity claim is requested or synthesized.
                        let bindings = vec![fs_material::section::relaxation::LayerBendingRelaxationProperties {
                            equilibrium_modulus: fs_couple::string_specimen::EQUILIBRIUM_YOUNG_MODULUS_PROPERTY.into(),
                            branches: properties.iter().map(|p|
                                (p.modulus.clone(), p.relaxation_time.clone())).collect(),
                        }; elastic.layers().len()];
                        MaterialStringSource::ConcentricProny { template: template(self),
                            section: elastic.with_prony_bending(&bindings)?, prestress }
                    }
                }
            }
        };
        // Torsion resolves from these SAME original layer states, not a second
        // query at averaged conditions or a mass copied from the nominal core.
        if twist { source = source.with_torsion(SHEAR, SHEAR_ETA, self.twist_modes)?; }
        if let Some(layers) = &self.torsion_memory {
            source = source.with_torsional_relaxation(layers)?;
        }
        if let Some(directions) = self.polarizations {
            source = source.with_transverse_polarizations(directions)?;
        }
        poll(gate)?;
        Ok(source)
    }
}
