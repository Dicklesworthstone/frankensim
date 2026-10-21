//! Temperature-driven sections with one retained material assignment per layer.
//!
//! Every property in the original bundle is re-evaluated, including auxiliary
//! electrical, thermal and loss properties. The exact original claims remain
//! pinned; changing temperature cannot silently select another phase, process,
//! pressure, composition or property table. Geometry is fixed on this rung.
//! This is a constitutive evaluation, not a heat integrator or an assertion of
//! constant mass when a source density varies with temperature.

/// Live layer-resolved electrical heating through the shared thermal timeline.
pub mod heating;
/// Source-pinned expansion and stress on the same layer temperatures.
pub mod expansion;

use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{QuantitySpec, Temperature};

use crate::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use crate::state_point::MaterialStatePointError;

#[derive(Clone, Debug)]
struct LayerBinding<'a> {
    card: &'a MaterialCard,
    point: QueryPoint,
    axis: String,
}

/// A reusable data-only temperature law for the existing layered section.
///
/// Source cards are borrowed, not copied or replaced by equivalent materials.
/// The section's layer order and each complete property requirement survive
/// every evaluation, including distinct same-name tensor/test-context keys.
#[derive(Clone, Debug)]
pub struct ConcentricSectionTemperatureLaw<'a> {
    reference: ConcentricElasticSection,
    bindings: Vec<LayerBinding<'a>>,
}

impl<'a> ConcentricSectionTemperatureLaw<'a> {
    /// Bind the exact cards behind an already resolved section. Supply one card
    /// and one absolute-temperature coordinate per layer, in spatial order.
    /// Legacy untyped `T` is kelvin; other names require a temperature descriptor.
    /// Typed temperature differences are refused rather than reinterpreted.
    ///
    /// # Errors
    /// Refuses mismatched card/count/temperature data. Original material
    /// requirements and source pins are preserved, never reconstructed from
    /// chemistry names or from a fixed list of assumed properties.
    pub fn try_new(
        reference: ConcentricElasticSection,
        cards: &[&'a MaterialCard],
        temperature_axes: &[&str],
    ) -> Result<Self, MaterialStatePointError> {
        if cards.len() != reference.layers().len() || temperature_axes.len() != cards.len() {
            return Err(invalid("section temperature law needs one card and axis per layer"));
        }
        let absolute = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature, ValueForm::Static,
        ));
        let mut bindings = Vec::with_capacity(cards.len());
        for ((layer, &card), &axis) in reference.layers().iter().zip(cards).zip(temperature_axes) {
            if card.content_hash() != layer.material.card_identity() {
                return Err(invalid("section temperature card differs from the original layer card"));
            }
            let descriptor_ok = match layer.material.axis_quantities().get(axis) {
                Some(q) => *q == absolute || *q == QuantitySpec::dimensional(Temperature::DIMS),
                None => axis == "T",
            };
            if !descriptor_ok || !layer.material.query_point().iter()
                .any(|(name, value)| name == axis && value.is_finite() && *value > 0.0)
            {
                return Err(invalid("section temperature law requires positive absolute-temperature coordinates"));
            }
            bindings.push(LayerBinding {
                card, point: layer.material.query(), axis: axis.to_owned(),
            });
        }
        Ok(Self { reference, bindings })
    }

    /// Original mechanical section, geometry, source states and auxiliary data.
    #[must_use]
    pub const fn reference(&self) -> &ConcentricElasticSection { &self.reference }

    /// Resolve mechanics AND every retained auxiliary property at the supplied
    /// per-layer absolute temperatures [K]. The whole candidate is returned
    /// only after every material query and section integral succeeds.
    /// Source support is checked at each requested state; this is not an
    /// assertion of continuous thermal-path or phase-stability support.
    pub fn at_temperatures_k(
        &self,
        temperatures_k: &[f64],
    ) -> Result<ConcentricElasticSection, MaterialStatePointError> {
        if temperatures_k.len() != self.bindings.len()
            || temperatures_k.iter().any(|t| !t.is_finite() || *t <= 0.0)
        {
            return Err(invalid("section needs one positive finite temperature per layer"));
        }
        let mut layers = Vec::with_capacity(self.bindings.len());
        for ((binding, geometry), &temperature) in self.bindings.iter()
            .zip(self.reference.layers()).zip(temperatures_k)
        {
            let point = binding.point.clone().with(&binding.axis, temperature)
                .map_err(|_| invalid("section temperature violates its retained convention"))?;
            let material = geometry.material.at_query_point(binding.card, &point)?;
            layers.push(ConcentricMaterialLayer { outer_radius_m: geometry.outer_radius_m, material });
        }
        ConcentricElasticSection::try_new(self.reference.inner_radius_m(), layers)
    }
}

fn invalid(quantity: &'static str) -> MaterialStatePointError {
    MaterialStatePointError::InvalidDerived { quantity }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conductor::{
        ELECTRICAL_RESISTIVITY_DIMS, ELECTRICAL_RESISTIVITY_PROPERTY, OhmicDrive,
    };
    use crate::state_point::{
        DENSITY_PROPERTY, MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS,
        SPECIFIC_HEAT_CAPACITY_PROPERTY, ScalarAdmissibility, ScalarPropertyRequirement,
        THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY, YOUNG_MODULUS_PROPERTY,
        resolve_material_state_point,
    };
    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        ClaimSet, InterpolationPolicy, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue,
        Provenance, UncertaintyModel,
    };
    use fs_qty::{Density, Dims, Pressure};

    pub(super) fn card(scale: f64) -> (MaterialCard, Vec<ScalarPropertyRequirement>) {
        let mut claims = ClaimSet::new();
        let mut requirements = Vec::new();
        for (name, dims, values) in [
            (DENSITY_PROPERTY, Density::DIMS, [8000.0, 8000.0]),
            (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, [2e11, 1e11]),
            (SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS, [400.0, 800.0]),
            (THERMAL_CONDUCTIVITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS, [20.0, 40.0]),
            (ELECTRICAL_RESISTIVITY_PROPERTY, ELECTRICAL_RESISTIVITY_DIMS, [1e-7 * scale, 2e-7 * scale]),
            ("auxiliary", Dims([0; 6]), [2.0, 3.0]),
        ] {
            claims.insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                    knots: vec![(300.0, values[0]), (400.0, values[1])], dims },
                validity: ValidityDomain::unconstrained().with("T", 300.0, 400.0)
                    .with("pressure", 100000.0, 100000.0),
                uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
                observations: vec![], provenance: Provenance {
                    source: "synthetic section temperature test, not measurements".into(),
                    license: "test-only".into(), artifact: None },
            }).unwrap();
            requirements.push(ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive).unwrap());
        }
        (MaterialCard::assemble(MaterialStateId { chemistry: "synthetic".into(),
            phase: "solid".into(), process: "test".into(), revision: 0 }, claims, vec![]).unwrap(), requirements)
    }

    pub(super) fn section(card: &MaterialCard, requirements: &[ScalarPropertyRequirement]) -> ConcentricElasticSection {
        let point = QueryPoint::new().with("T", 300.0).unwrap().with("pressure", 100000.0).unwrap();
        let material = resolve_material_state_point(card, &point, requirements,
            MaterialPropertySelection::SingleClaimOnly).unwrap();
        ConcentricElasticSection::try_new(0.0, vec![
            ConcentricMaterialLayer { outer_radius_m: 0.001, material: material.clone() },
            ConcentricMaterialLayer { outer_radius_m: 0.002, material },
        ]).unwrap()
    }

    #[test]
    fn temperature_updates_mechanics_heat_and_current_from_the_same_pins() {
        let (card, requirements) = card(1.0);
        let reference = section(&card, &requirements);
        let law = ConcentricSectionTemperatureLaw::try_new(reference.clone(), &[&card, &card], &["T", "T"]).unwrap();
        let hot = law.at_temperatures_k(&[400.0, 300.0]).unwrap();
        assert!(hot.axial_stiffness_n() < reference.axial_stiffness_n());
        assert!(hot.bending_stiffness_n_m2() < reference.bending_stiffness_n_m2());
        assert_eq!(hot.line_density_kg_m(), reference.line_density_kg_m());
        let old = reference.thermal_transport().unwrap().with_ohmic_transport(1.0).unwrap().evaluate(OhmicDrive::Voltage(1.0), 10.0).unwrap();
        let new = hot.thermal_transport().unwrap().with_ohmic_transport(1.0).unwrap().evaluate(OhmicDrive::Voltage(1.0), 10.0).unwrap();
        assert!(new.joule_power_w < old.joule_power_w);
        assert!(new.layers[0].heat_capacity_j_k > old.layers[0].heat_capacity_j_k);
        assert_eq!(new.layers[1], old.layers[1]);
        assert!(new.axial_heat_flow_w.abs() > old.axial_heat_flow_w.abs());
        for (a, b) in reference.layers()[0].material.properties().iter().zip(hot.layers()[0].material.properties()) {
            assert_eq!(a.answer().receipt.selected, b.answer().receipt.selected);
        }
        assert_eq!(hot.layers()[0].material.property("auxiliary").unwrap().value_si(), 3.0);
        assert_eq!(hot.layers()[0].material.query_point().iter().find(|(k, _)| k == "pressure").unwrap().1, 100000.0);
        assert_eq!(law.reference(), &reference);
    }

    #[test]
    fn wrong_card_counts_and_unsupported_temperatures_refuse() {
        let (a, requirements) = card(1.0);
        let (b, _) = card(2.0);
        let reference = section(&a, &requirements);
        assert!(ConcentricSectionTemperatureLaw::try_new(reference.clone(), &[&a], &["T"]).is_err());
        assert!(ConcentricSectionTemperatureLaw::try_new(reference.clone(), &[&a, &b], &["T", "T"]).is_err());
        assert!(ConcentricSectionTemperatureLaw::try_new(reference.clone(), &[&a, &a], &["pressure", "T"]).is_err());
        let law = ConcentricSectionTemperatureLaw::try_new(reference.clone(), &[&a, &a], &["T", "T"]).unwrap();
        for temperatures in [vec![300.0], vec![300.0, 401.0], vec![f64::NAN, 300.0], vec![0.0, 300.0]] {
            assert!(law.at_temperatures_k(&temperatures).is_err());
            assert_eq!(law.reference(), &reference);
        }
        let restored = law.at_temperatures_k(&[300.0, 300.0]).unwrap();
        assert_eq!(restored.axial_stiffness_n(), reference.axial_stiffness_n());
        assert_eq!(restored.area_m2(), reference.area_m2());
    }
}
