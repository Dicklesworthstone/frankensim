//! Data-only construction over the original specimen, concentric-section and
//! torsion owners. No section formula, support law or material fallback lives here.
use super::*;
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::state_point::integrate_isotropic_thermal_expansion;
use fs_couple::string_specimen::with_uniform_circular_thermal_extension;

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
                if !reference.is_finite() || reference <= 0.0 || !matches!(&self.section, input::WireSection::Homogeneous) {
                    return Err("thermal support requires positive reference temperature and a homogeneous string".into());
                }
                Some((reference, stress_free_length_m, linear_strain_limit))
            },
            None => None,
        };
        let twist = self.twist_modes > 0;
        let core = resolve(&self.pack, catalog, loader, spec, false, twist)?;
        let mut source = match &self.section {
            input::WireSection::Homogeneous => {
                // V1 and explicit V2 fixed-tension specimens take precisely the
                // old route. The extension variant binds the original support law.
                let specimen = if let Some((reference, length, strain_limit)) = thermal {
                    // Reuse the exact card/current state; integrate the selected
                    // alpha(T) path across every knot, never endpoint-alpha*dT.
                    let mut current = QueryPoint::new();
                    for (axis, value) in core.query_point() {
                        current = match core.axis_quantities().get(axis) {
                            Some(quantity) => current.with_quantity(axis.clone(), *quantity, *value)?,
                            None => current.with(axis.clone(), *value)?,
                        };
                    }
                    let point = current.clone().with("T", reference)?;
                    let expansion = integrate_isotropic_thermal_expansion(&catalog[&self.pack].card,
                        &point, &current, Pick::SingleClaimOnly)?;
                    with_uniform_circular_thermal_extension(template(self),
                        StringGeometryConstraint::FixedRadius(self.radius), &core, &expansion,
                        length, strain_limit)?
                } else { match self.prestress {
                    StringPrestress::FixedTension(_) =>
                        with_uniform_circular_material_state(template(self), self.radius, &core)?,
                    _ => with_uniform_circular_material_and_constraints(template(self),
                        StringGeometryConstraint::FixedRadius(self.radius), &core, self.prestress)?,
                } }.with_kelvin_voigt_bending_loss()?;
                MaterialStringSource::Homogeneous(specimen)
            }
            input::WireSection::Bonded { core_radius_m, sleeves } => {
                let mut layers = Vec::with_capacity(sleeves.len() + 1);
                layers.push(ConcentricMaterialLayer { outer_radius_m: *core_radius_m, material: core });
                for sleeve in sleeves {
                    poll(gate)?;
                    layers.push(ConcentricMaterialLayer { outer_radius_m: sleeve.outer_radius_m,
                        material: resolve(&sleeve.pack, catalog, loader, spec, false, twist)? });
                }
                let elastic = ConcentricElasticSection::try_new(0.0, layers)?;
                if elastic.outer_radius_m().to_bits() != self.radius.to_bits() {
                    return Err("resolved section disagrees with the original string radius".into());
                }
                MaterialStringSource::Concentric { template: template(self),
                    section: elastic.with_kelvin_voigt_bending(ETA)?, prestress: self.prestress }
            }
        };
        // Torsion resolves from these SAME original layer states, not a second
        // query at averaged conditions or a mass copied from the nominal core.
        if twist { source = source.with_torsion(SHEAR, SHEAR_ETA, self.twist_modes)?; }
        if let Some(directions) = self.polarizations {
            source = source.with_transverse_polarizations(directions)?;
        }
        poll(gate)?;
        Ok(source)
    }
}
