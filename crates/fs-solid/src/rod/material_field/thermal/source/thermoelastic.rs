//! Add the original pinned expansion claim to the original mechanical bundle.
//! No claim selection, geometry or material identity changes at this boundary.

use super::*;
use fs_material::section_temperature::{
    ConcentricSectionTemperatureLaw, thermoelastic::rod::ConcentricRodThermoelasticLaw,
};
use fs_material::state_point::INVERSE_TEMPERATURE_DIMS;

impl<'a> RodMaterialSource<'a> {
    /// Prepare complete source-derived axial/bending/shear/torsion potentials.
    /// The thermal owner admits affine elastic moduli and scalar rho/alpha;
    /// unsupported tables remain supported by prescribed mechanics only.
    pub fn thermoelastic_laws(
        &self,
        gate: &CancelGate,
    ) -> Result<Vec<ConcentricRodThermoelasticLaw<'a>>, SpatialRodError> {
        checkpoint(gate)?;
        let mut requirements = self.requirements.clone();
        requirements.push(ScalarPropertyRequirement::try_new(
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
            INVERSE_TEMPERATURE_DIMS,
            ScalarAdmissibility::Finite,
        )?);
        let mut laws = Vec::with_capacity(self.segments.len());
        for (segment, reference) in self.segments.iter().zip(self.reference.sections()) {
            let mut layers = Vec::with_capacity(segment.layers.len());
            for (source, original) in segment.layers.iter().zip(reference.elastic().layers()) {
                checkpoint(gate)?;
                let mut pins = original
                    .material
                    .properties()
                    .iter()
                    .map(|p| {
                        (
                            p.requirement().name().to_owned(),
                            p.answer().receipt.selected,
                        )
                    })
                    .collect::<Vec<_>>();
                let MaterialPropertySelection::PinnedByProperty(expansion) =
                    &source.expansion_selection
                else {
                    return Err(SpatialRodError::Invalid(
                        "rod expansion source must remain pinned",
                    ));
                };
                pins.extend(expansion.iter().cloned());
                let material = resolve_material_state_point(
                    source.card,
                    &source.reference_point,
                    &requirements,
                    MaterialPropertySelection::PinnedByProperty(pins),
                )?;
                layers.push(ConcentricMaterialLayer {
                    outer_radius_m: source.outer_radius_m,
                    material,
                });
            }
            let section = ConcentricElasticSection::try_new(segment.inner_radius_m, layers)?;
            let cards = segment.layers.iter().map(|s| s.card).collect::<Vec<_>>();
            let temperature =
                ConcentricSectionTemperatureLaw::try_new(section, &cards, &vec!["T"; cards.len()])?;
            laws.push(ConcentricRodThermoelasticLaw::try_new(
                temperature,
                &self.shear_property,
                segment.shear_correction,
            )?);
        }
        checkpoint(gate)?;
        Ok(laws)
    }
}
