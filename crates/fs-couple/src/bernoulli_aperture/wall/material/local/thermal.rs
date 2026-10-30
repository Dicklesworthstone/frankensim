//! Temperature changes of the original fixed-mass material wall source.

use super::{AcousticRealizeError, LocalMaterialWall, WallMaterialLoss, invalid};
use fs_matdb::MaterialCard;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{QuantitySpec, Temperature};

impl LocalMaterialWall {
    /// Re-evaluate the original material wall bundle at one finite body's
    /// temperature. The shared-body owner settles the native rebind's signed
    /// work separately; this operation only resolves the unchanged source law.
    ///
    /// All original claims, coordinate conventions, thickness, restraint and
    /// applicability survive. A fixed geometry and inventory require unchanged
    /// density. Elastic, source Young-viscosity Kelvin--Voigt and source-pair
    /// Maxwell laws retain their original constitutive owner. Native Maxwell
    /// transfer preserves each original arm's elastic extension and viscous
    /// offset; no memory is authored, matched by coefficient or reset here.
    /// Principal-aligned orthotropic elasticity and proportional Kelvin--Voigt
    /// response retain their original axes, restraint and time-property role.
    /// Generated normal-plane viscosity/memory, phase change and thermal
    /// expansion still require their own physical transfer model.
    pub(crate) fn at_temperature_pinned(
        &self,
        card: &MaterialCard,
        temperature_axis: &str,
        initial_temperature_k: f64,
        temperature_k: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AcousticRealizeError> {
        if cancelled() {
            return Err(AcousticRealizeError::Cancelled);
        }
        let quantity = self
            .material
            .axis_quantities()
            .get(temperature_axis)
            .copied();
        let temperature_coordinate = quantity
            == Some(QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )))
            || temperature_axis == "T"
                && quantity.is_none_or(|q| q == QuantitySpec::dimensional(Temperature::DIMS));
        if !initial_temperature_k.is_finite()
            || initial_temperature_k <= 0.0
            || !temperature_k.is_finite()
            || temperature_k <= 0.0
            || !temperature_coordinate
            || self.material.card_identity() != card.content_hash()
            || self
                .material
                .query_point()
                .iter()
                .find(|(name, _)| name == temperature_axis)
                .map(|(_, value)| *value)
                != Some(initial_temperature_k)
        {
            return Err(invalid(
                "thermal wall feedback requires the original material card and accepted body temperature",
            ));
        }
        let point = self
            .material
            .query()
            .with(temperature_axis, temperature_k)
            .map_err(|error| AcousticRealizeError::Nonlinear(error.to_string()))?;
        let material = self
            .material
            .at_query_point(card, &point)
            .map_err(|error| AcousticRealizeError::Nonlinear(error.to_string()))?;
        let next = if let Some(axes) = self.orthotropic_axes {
            let time = self.proportional_viscosity_time_property.as_deref();
            // The separate normal-plane viscosity constructor can also produce
            // Kelvin or Maxwell loss. Its tensor-property roles and induced
            // axial-strain transfer cannot be replaced by proportional damping.
            match (&self.loss, time) {
                (WallMaterialLoss::Elastic, None) | (WallMaterialLoss::KelvinVoigt, Some(_)) => {}
                _ => {
                    return Err(invalid(
                        "thermal orthotropic wall feedback requires original elasticity or proportional viscosity; normal-plane laws need their own source transfer",
                    ));
                }
            }
            Self::new_orthotropic(
                material,
                self.thickness_m,
                self.axial,
                axes,
                time,
                self.band_hz,
            )?
        } else {
            if self.proportional_viscosity_time_property.is_some() {
                return Err(invalid(
                    "thermal wall feedback cannot change a directional viscosity role",
                ));
            }
            let loss = match &self.loss {
                WallMaterialLoss::Elastic => WallMaterialLoss::Elastic,
                // The complete pinned bundle retains the original viscosity.
                WallMaterialLoss::KelvinVoigt => WallMaterialLoss::KelvinVoigt,
                WallMaterialLoss::Maxwell(source) => {
                    WallMaterialLoss::Maxwell(source.at_material_pinned(&self.material, &material)?)
                }
            };
            Self::new(material, self.thickness_m, self.axial, loss, self.band_hz)?
        };
        if next.density != self.density {
            return Err(invalid(
                "fixed-geometry thermal wall feedback cannot change material density or mass",
            ));
        }
        if cancelled() {
            return Err(AcousticRealizeError::Cancelled);
        }
        Ok(next)
    }
}
