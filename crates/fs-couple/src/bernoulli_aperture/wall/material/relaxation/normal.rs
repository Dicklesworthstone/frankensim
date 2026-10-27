//! The axial strain behind an exactly condensed normal Kelvin--Voigt law.
//! Its identity is the source tensor, not a fabricated modulus/time property.
use super::*;

#[derive(Clone, Debug, PartialEq)]
pub(super) struct NormalPlaneHistory {
    card: ContentHash,
    axes: [usize; 2],
    viscosity_properties: [String; 3],
    selected: Vec<(String, ClaimId)>,
    elastic_ratio: f64,
    mismatch: f64,
}

impl WallRelaxation {
    /// Called only by the original exact normal-plane condensation, after its
    /// passivity, nonzero-memory and representability checks. Retain all source
    /// roles and selected claims, not just the generated positive coefficients.
    pub(crate) fn with_normal_plane_history(
        mut self,
        material: &ResolvedMaterialStatePoint,
        axes: [usize; 2],
        viscosity_properties: [&str; 3],
        elastic_ratio: f64,
        mismatch: f64,
    ) -> Result<Self, AcousticRealizeError> {
        if self.source_arms.is_some() || self.normal_history.is_some()
            || self.material.terms.len() != 1 || self.material.terms[0].0 <= 0.0
            || !elastic_ratio.is_finite() || !mismatch.is_finite() || mismatch == 0.0
        {
            return Err(invalid("generated wall history needs its one positive original normal-plane mode"));
        }
        let requirements = super::super::LocalMaterialWall::orthotropic_normal_viscosity_requirements(
            axes, viscosity_properties,
        )?;
        let selected = requirements.iter().map(|requirement| {
            let name = requirement.name();
            material.property(name).map(|property|
                (name.to_owned(), property.answer().receipt.selected))
                .ok_or_else(|| invalid("generated wall history is missing an original source claim"))
        }).collect::<Result<Vec<_>, _>>()?;
        self.source_identity = Some(material.identity());
        self.normal_history = Some(NormalPlaneHistory {
            card: material.card_identity(), axes,
            viscosity_properties: viscosity_properties.map(str::to_owned),
            selected, elastic_ratio, mismatch,
        });
        Ok(self)
    }

    /// In the original tensor, c=C_ha/C_aa, d=c-D_ha/D_aa and
    /// y=z/K = radius*(epsilon_a+c*epsilon_h)/d. Preserving physical axial
    /// strain and radial displacement therefore requires
    /// y_new=(d_old/d_new)*y_old+(c_new-c_old)/d_new*w.
    /// Preserving y itself would reset epsilon_a whenever c or d changed.
    ///
    /// None selects the existing source-pair Maxwell transfer. Mixing those
    /// histories with a generated mode, changing source roles/axes/claims, or
    /// losing the induced mode refuses. The native wall owner separately keeps
    /// thickness, geometry and axial restraint fixed and budgets actual work.
    pub(crate) fn normal_history_extension_map(
        &self,
        replacement: &Self,
    ) -> Result<Option<[f64; 2]>, AcousticRealizeError> {
        let (old, new) = match (&self.normal_history, &replacement.normal_history) {
            (None, None) => return Ok(None),
            (Some(old), Some(new)) => (old, new),
            _ => return Err(invalid("generated axial strain cannot become an unrelated Maxwell history")),
        };
        if old.card != new.card || old.axes != new.axes
            || old.viscosity_properties != new.viscosity_properties || old.selected != new.selected
        {
            return Err(invalid("live normal-plane memory requires the same card, frame, property roles and selected claims"));
        }
        let difference = new.elastic_ratio - old.elastic_ratio;
        let map = [old.mismatch / new.mismatch, difference / new.mismatch];
        if map.iter().any(|v| !v.is_finite()) || map[0] == 0.0
            || (difference != 0.0 && map[1] == 0.0)
        {
            return Err(invalid("normal-plane axial-strain transfer is not representable"));
        }
        Ok(Some(map))
    }
}
