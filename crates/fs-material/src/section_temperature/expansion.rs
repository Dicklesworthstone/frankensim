//! Retained thermal expansion paths on the SAME temperature-driven section.
use super::{ConcentricSectionTemperatureLaw, invalid};
use crate::section::{ConcentricElasticSection, thermal::ConcentricThermoelasticSection};
use crate::state_point::{MaterialPropertySelection, MaterialStatePointError,
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY as ALPHA,
    integrate_isotropic_thermal_expansion};

impl ConcentricSectionTemperatureLaw<'_> {
    /// Re-evaluate the complete layer bundles AND integrate their original
    /// source expansion paths from the common stress-free reference temperature.
    /// Alpha must be explicitly present in each original bundle. Its exact
    /// source claim stays pinned along with mechanics, heat and electricity.
    /// Negative expansion is admissible; missing expansion never means zero.
    ///
    /// The existing expansion owner currently uses the canonical `T` axis.
    /// Other temperature axes remain supported by ordinary section evaluation,
    /// but are not silently renamed to enter this thermoelastic rung.
    /// Geometry is fixed; this supplies eigenstrain/stress, not remeshing or a
    /// heat/work balance. See `ConcentricThermoelasticSection` for its regime.
    pub fn thermoelastic_at_temperatures_k(
        &self, temperatures_k: &[f64],
    ) -> Result<ConcentricThermoelasticSection, MaterialStatePointError> {
        self.bind_expansion_paths(self.at_temperatures_k(temperatures_k)?)
    }

    // Heating supplies its already resolved endpoint, so observing strain does
    // not resample material properties at a second, possibly stale temperature.
    pub(super) fn bind_expansion_paths(
        &self, section: ConcentricElasticSection,
    ) -> Result<ConcentricThermoelasticSection, MaterialStatePointError> {
        let mut paths = Vec::with_capacity(self.bindings.len());
        for ((binding, reference), layer) in self.bindings.iter()
            .zip(self.reference.layers()).zip(section.layers())
        {
            if binding.axis != "T" {
                return Err(invalid("section expansion requires the existing canonical T coordinate"));
            }
            let alpha = reference.material.property(ALPHA)
                .ok_or_else(|| invalid("every thermally constrained layer needs explicit expansion data"))?;
            let temperature = layer.material.query_point().iter().find(|(key, _)| key == "T")
                .ok_or_else(|| invalid("thermal section endpoint lacks T"))?.1;
            let point = binding.point.clone().with("T", temperature)
                .map_err(|_| invalid("section expansion endpoint violates its coordinate convention"))?;
            paths.push(integrate_isotropic_thermal_expansion(
                binding.card, &binding.point, &point,
                MaterialPropertySelection::PinnedByProperty(vec![(
                    ALPHA.to_owned(), alpha.answer().receipt.selected,
                )]),
            )?);
        }
        ConcentricThermoelasticSection::try_new(section, paths)
    }
}

#[cfg(test)]
mod tests;
