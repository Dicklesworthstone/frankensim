//! Pinned temperature continuation of original ordered-interface spectra.
use super::*;
use fs_qty::Temperature;

impl ResolvedViscoelasticBond {
    /// Requery every SAME selected coefficient claim at one new temperature.
    /// Geometry, branch identities/counts, model-use bounds and every other
    /// physical query coordinate remain fixed. This does not advance memory.
    pub fn at_temperature(
        &self,
        card: &InterfaceSystemCard,
        temperature_axis: &str,
        temperature_k: f64,
    ) -> Result<Self, MaterialStatePointError> {
        self.require_temperature_axis(temperature_axis)?;
        if !temperature_k.is_finite() || temperature_k <= 0.0 {
            return Err(invalid("bond temperature requires finite positive kelvin"));
        }
        let point = self
            .source
            .query()
            .with(temperature_axis, temperature_k)
            .map_err(|source| MaterialStatePointError::Query {
                property: format!("bond temperature:{temperature_axis}"),
                source,
            })?;
        let resolved = self.source.at_query_point(card, &point)?;
        let next = Self::try_from_resolved(card, resolved, self.law.limits())?;
        self.require_temperature_successor(&next, temperature_axis)?;
        Ok(next)
    }

    /// Establish source authority BEFORE transferring live physical memory.
    /// Canonical indexed modulus/time claims identify branches; equal numeric
    /// coefficients or equal counts do not permit replacement or reordering.
    /// The subsequent numeric transfer preserves those indexed spring strains.
    pub fn require_temperature_successor(
        &self,
        next: &Self,
        temperature_axis: &str,
    ) -> Result<(), MaterialStatePointError> {
        self.require_temperature_axis(temperature_axis)?;
        next.require_temperature_axis(temperature_axis)?;
        if self.source.card_identity() != next.source.card_identity()
            || self.surfaces != next.surfaces
            || self.source.axis_quantities() != next.source.axis_quantities()
            || self.law.thickness_m() != next.law.thickness_m()
            || self.law.limits() != next.law.limits()
            || self.law.normal_spectrum().terms.len() != next.law.normal_spectrum().terms.len()
            || self.law.shear_spectrum().terms.len() != next.law.shear_spectrum().terms.len()
        {
            return Err(invalid(
                "bond temperature successor changed original source, geometry, bounds or branches",
            ));
        }
        let a = self.source.query_point();
        let b = next.source.query_point();
        if a.len() != b.len()
            || a.iter().zip(b).any(|((axis, value), (other, new))| {
                axis != other || (axis != temperature_axis && value != new)
            })
        {
            return Err(invalid(
                "bond temperature successor changed a held query coordinate",
            ));
        }
        let a = self.source.properties();
        let b = next.source.properties();
        if a.len() != b.len()
            || a.iter().zip(b).any(|(old, new)| {
                old.requirement() != new.requirement()
                    || old.answer().receipt.selected != new.answer().receipt.selected
            })
        {
            return Err(invalid(
                "bond temperature successor changed original indexed coefficient claims",
            ));
        }
        Ok(())
    }

    fn require_temperature_axis(&self, axis: &str) -> Result<(), MaterialStatePointError> {
        if self.source.axis_quantities().get(axis)
            != Some(&QuantitySpec::dimensional(Temperature::DIMS))
            || !self
                .source
                .query_point()
                .iter()
                .any(|(name, value)| name == axis && value.is_finite() && *value > 0.0)
        {
            return Err(invalid(
                "bond temperature continuation requires an original kelvin coordinate",
            ));
        }
        Ok(())
    }
}
