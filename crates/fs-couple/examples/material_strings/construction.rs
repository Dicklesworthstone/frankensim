//! Data-only construction over the original specimen, concentric-section and
//! torsion owners. No section formula, support law or material fallback lives here.
use super::*;
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};

impl input::Wire {
    pub(super) fn resolve_source(&self, catalog: &mut Catalog,
        loader: &mut impl FnMut(&str) -> Result<NormalizedPack>, spec: &input::Spec,
        gate: &CancelGate) -> Result<MaterialStringSource>
    {
        poll(gate)?;
        let twist = self.twist_modes > 0;
        let core = resolve(&self.pack, catalog, loader, spec, false, twist)?;
        let mut source = match &self.section {
            input::WireSection::Homogeneous => {
                // V1 and explicit V2 fixed-tension specimens take precisely the
                // old route. The extension variant binds the original support law.
                let specimen = match self.prestress {
                    StringPrestress::FixedTension(_) =>
                        with_uniform_circular_material_state(template(self), self.radius, &core)?,
                    _ => with_uniform_circular_material_and_constraints(template(self),
                        StringGeometryConstraint::FixedRadius(self.radius), &core, self.prestress)?,
                }.with_kelvin_voigt_bending_loss()?;
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
        poll(gate)?;
        Ok(source)
    }
}
