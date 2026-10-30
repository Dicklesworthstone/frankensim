//! Temperature changes of the original fixed-mass elastic material source.

use super::{AcousticRealizeError, LocalMaterialWall, WallMaterialLoss, invalid};
use fs_matdb::MaterialCard;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{QuantitySpec, Temperature};

impl LocalMaterialWall {
    /// Re-evaluate the original elastic wall bundle at one finite body's
    /// temperature. The shared-body owner settles the native rebind's signed
    /// work separately; this operation only resolves the unchanged source law.
    ///
    /// All original claims, coordinate conventions, thickness, restraint and
    /// applicability survive. A fixed geometry and inventory require unchanged
    /// density. This bounded adapter has no anisotropic, hereditary, phase or
    /// thermal-expansion transfer model; it cannot substitute one silently.
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
        if self.orthotropic_axes.is_some()
            || self.proportional_viscosity_time_property.is_some()
            || self.loss != WallMaterialLoss::Elastic
            || !initial_temperature_k.is_finite()
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
                "thermal wall feedback requires the original isotropic elastic card and accepted body temperature",
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
        let next = Self::new(
            material,
            self.thickness_m,
            self.axial,
            WallMaterialLoss::Elastic,
            self.band_hz,
        )?;
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
