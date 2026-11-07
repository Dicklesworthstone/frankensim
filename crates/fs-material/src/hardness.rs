//! Original indentation hardness shared by normal contact and wear.
//!
//! Hardness in pascals is a force/area response under a declared indentation
//! test, not a yield stress, elastic modulus, fracture toughness or empirical
//! HRC/Shore reading. Source keys retain the test apparatus, load/dwell history,
//! protocol, specimen and area convention. No scale conversion or material-name
//! inference occurs here.

/// Original hardness evaluated along the native accumulated plastic history.
pub mod plastic_history;

/// Original hardness at material depths exposed by actual surface removal.
pub mod depth;

/// Original hardness at accepted native thermoplastic temperature/history.
pub mod thermal_history;

use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::{MaterialCard, PropertyKey, QueryPoint};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Pressure, QuantitySpec};

use crate::state_point::{
    MaterialStatePointError, ResolvedMaterialStatePoint, ResolvedScalarProperty,
};

/// Static force/area hardness, distinct from stress and fluid pressure.
pub const INDENTATION_HARDNESS_QUANTITY: QuantitySpec = QuantitySpec::semantic(SemanticType::new(
    QuantityKind::IndentationHardness,
    ValueForm::Static,
));

/// Admit one explicitly selected positive indentation-hardness value [Pa].
///
/// Typed inputs require the exact original test context. Existing callers may
/// explicitly select a dimension-only Pa claim as a *declared* hardness
/// convention; that path does not make it a test-qualified measurement or
/// change its source authority. Typed pressure/stress, empirical scales and
/// tensor components cannot stand in for hardness. No property-name inference
/// occurs. Consumers separately admit the source coordinates they can own.
pub fn indentation_hardness_property<'a>(
    resolved: &'a ResolvedMaterialStatePoint,
    key: &PropertyKey,
) -> Result<&'a ResolvedScalarProperty, MaterialStatePointError> {
    let quantity = key.quantity();
    if (quantity != INDENTATION_HARDNESS_QUANTITY
        && quantity != QuantitySpec::dimensional(Pressure::DIMS))
        || key.elastic_component().is_some()
        || key.strain_component().is_some()
        || key.stress_component().is_some()
        || (quantity == INDENTATION_HARDNESS_QUANTITY && key.hardness_test().is_none())
    {
        return Err(invalid(
            "indentation hardness requires its exact static Pa convention and test context",
        ));
    }
    let property = resolved.property_by_key(key).ok_or_else(|| {
        invalid("indentation hardness is missing its exact original property key")
    })?;
    if !(property.value_si().is_finite() && property.value_si() > 0.0) {
        return Err(invalid(
            "indentation hardness must be finite and strictly positive in pascals",
        ));
    }
    Ok(property)
}

/// Typed indentation hardness projected from a complete immutable material
/// bundle. Extra properties and their exact source selections remain intact.
///
/// Refresh evaluates the same claims at a new physical state. It does not
/// evolve plastic strain, wear, temperature or test conditions. A consumer
/// must check [`Self::property`] source axes against the coordinates it owns
/// before using this endpoint in a constitutive path.
#[derive(Clone, Debug, PartialEq)]
pub struct IndentationHardnessStatePoint {
    resolved: ResolvedMaterialStatePoint,
    key: PropertyKey,
    property_index: usize,
    identity: ContentHash,
}

impl IndentationHardnessStatePoint {
    /// Adopt an explicitly typed original hardness claim and its exact test.
    /// Dimension-only estimates are intentionally confined to the borrowed
    /// legacy adapter [`indentation_hardness_property`].
    pub fn try_from_resolved(
        resolved: ResolvedMaterialStatePoint,
        key: PropertyKey,
    ) -> Result<Self, MaterialStatePointError> {
        if key.quantity() != INDENTATION_HARDNESS_QUANTITY {
            return Err(invalid(
                "source indentation hardness requires the pressure-valued hardness kind",
            ));
        }
        let property = indentation_hardness_property(&resolved, &key)?;
        let property_index = resolved
            .properties()
            .iter()
            .position(|candidate| core::ptr::eq(candidate, property))
            .ok_or_else(|| invalid("hardness property is not retained by its source bundle"))?;
        let mut identity = DomainHasher::new("org.frankensim.fs-material.indentation-hardness.v1");
        identity.update(resolved.identity().as_bytes());
        identity.update(&(property_index as u64).to_le_bytes());
        Ok(Self {
            resolved,
            key,
            property_index,
            identity: identity.finalize(),
        })
    }

    /// Complete original bundle, including all selected coupled properties.
    #[must_use]
    pub const fn resolved(&self) -> &ResolvedMaterialStatePoint {
        &self.resolved
    }

    /// Exact hardness key, including apparatus, loading, protocol and specimen.
    #[must_use]
    pub const fn key(&self) -> &PropertyKey {
        &self.key
    }

    /// Original selected value, evidence, usage receipt and source coordinates.
    #[must_use]
    pub fn property(&self) -> &ResolvedScalarProperty {
        &self.resolved.properties()[self.property_index]
    }

    /// Positive force/area hardness in coherent SI [Pa].
    #[must_use]
    pub fn hardness_pa(&self) -> f64 {
        self.property().value_si()
    }

    /// Identity binds the complete source bundle and its selected hardness key.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Refresh all retained properties through their exact original sources.
    /// Changed cards, quantity conventions, test context, missing coverage or
    /// nonpositive hardness refuse without altering this state.
    pub fn at_query_point(
        &self,
        card: &MaterialCard,
        point: &QueryPoint,
    ) -> Result<Self, MaterialStatePointError> {
        Self::try_from_resolved(self.resolved.at_query_point(card, point)?, self.key.clone())
    }
}

fn invalid(quantity: &'static str) -> MaterialStatePointError {
    MaterialStatePointError::InvalidDerived { quantity }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::state_point::{
        MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
        resolve_material_state_point,
    };
    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        ClaimSet, HardnessLoadStep, HardnessTestContext, InterpolationPolicy, MaterialStateId,
        ObservationDataset, PropertyClaim, PropertyValue, Provenance, UncertaintyModel,
    };
    use fs_qty::{Dims, QtyAny, Temperature};

    fn point(temperature: f64) -> QueryPoint {
        QueryPoint::new().with("T", temperature).unwrap()
    }

    fn fixture() -> (MaterialCard, [PropertyKey; 2]) {
        let provenance = Provenance {
            source: "synthetic indentation fixture; no measured material claim".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        };
        let mut claims = ClaimSet::new();
        let observation = claims
            .register_observation(ObservationDataset {
                specimen: "synthetic homogeneous solid".into(),
                method: "synthetic force/projected-area test".into(),
                artifact: fs_blake3::hash_domain("fixture", b"source indentation hardness"),
                caveats: "numerical fixture only".into(),
                provenance: provenance.clone(),
            })
            .unwrap();
        let keys = [10.0, 20.0].map(|load| {
            PropertyKey::with_quantity("H", INDENTATION_HARDNESS_QUANTITY)
                .with_hardness_test(
                    HardnessTestContext::new(
                        "synthetic spherical diamond indenter, radius 1 mm",
                        vec![
                            HardnessLoadStep::new(
                                QtyAny::new(load, Dims([1, 1, -2, 0, 0, 0])),
                                QtyAny::new(5.0, Dims([0, 0, 1, 0, 0, 0])),
                            )
                            .unwrap(),
                        ],
                        "synthetic projected-contact-area protocol revision 1",
                        observation,
                    )
                    .unwrap(),
                )
                .unwrap()
        });
        for (key, factor) in keys.iter().zip([1.0, 2.0]) {
            claims
                .insert_claim(PropertyClaim {
                    key: key.clone(),
                    value: PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: Temperature::DIMS,
                        knots: vec![(280.0, factor * 5e8), (380.0, factor * 2.5e8)],
                        dims: Pressure::DIMS,
                    },
                    validity: ValidityDomain::unconstrained().with("T", 280.0, 380.0),
                    uncertainty: UncertaintyModel::Unstated,
                    interpolation: InterpolationPolicy::LinearInside,
                    observations: vec![observation],
                    provenance: provenance.clone(),
                })
                .unwrap();
        }
        let card = MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic".into(),
                phase: "solid".into(),
                process: "synthetic indentation specimen".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap();
        (card, keys)
    }

    fn resolved(card: &MaterialCard, keys: &[PropertyKey]) -> ResolvedMaterialStatePoint {
        let requirements = keys
            .iter()
            .map(|key| {
                ScalarPropertyRequirement::try_with_key(key, ScalarAdmissibility::Finite).unwrap()
            })
            .collect::<Vec<_>>();
        resolve_material_state_point(
            card,
            &point(280.0),
            &requirements,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap()
    }

    #[test]
    fn g3_source_hardness_refresh_retains_exact_tests_claims_and_full_bundle() {
        let (card, keys) = fixture();
        let bundle = resolved(&card, &keys);
        assert!(
            bundle.property("H").is_none(),
            "test contexts must not collapse by name"
        );
        let source =
            IndentationHardnessStatePoint::try_from_resolved(bundle.clone(), keys[0].clone())
                .unwrap();
        let other =
            IndentationHardnessStatePoint::try_from_resolved(bundle.clone(), keys[1].clone())
                .unwrap();
        assert_eq!(source.hardness_pa(), 5e8);
        assert_eq!(other.hardness_pa(), 1e9);
        assert_ne!(source.identity(), other.identity());
        assert_eq!(source.resolved(), &bundle);
        assert_eq!(
            source.property().requirement().hardness_test(),
            keys[0].hardness_test()
        );
        assert_eq!(
            source.property().source_axes(),
            &[("T".to_string(), Temperature::DIMS)]
        );
        assert!(
            source
                .property()
                .has_only_fixed_external_conditions(Some("T"), false)
        );
        let next = source.at_query_point(&card, &point(330.0)).unwrap();
        assert_eq!(next.hardness_pa(), 3.75e8);
        assert_eq!(next.resolved().properties().len(), 2);
        assert_eq!(
            next.resolved()
                .property_by_key(&keys[1])
                .unwrap()
                .value_si(),
            7.5e8
        );
        assert_eq!(next.key(), source.key());
        assert_eq!(
            next.property().answer().receipt.selected,
            source.property().answer().receipt.selected
        );
        assert_ne!(next.identity(), source.identity());
        assert!(source.at_query_point(&card, &point(400.0)).is_err());
        let replacement = MaterialCard::supersede(&card, card.claims().clone(), vec![]).unwrap();
        assert!(source.at_query_point(&replacement, &point(330.0)).is_err());
        assert_eq!(
            source.hardness_pa(),
            5e8,
            "failed refresh leaves original untouched"
        );
    }

    #[test]
    fn g0_hardness_admission_refuses_wrong_semantics_context_and_nonpositive_values() {
        let (card, keys) = fixture();
        let bundle = resolved(&card, &keys[..1]);
        assert!(indentation_hardness_property(&bundle, &keys[1]).is_err());
        assert!(
            indentation_hardness_property(
                &bundle,
                &PropertyKey::with_quantity("H", INDENTATION_HARDNESS_QUANTITY)
            )
            .is_err()
        );
        for kind in [
            QuantityKind::Pressure,
            QuantityKind::Stress,
            QuantityKind::Hardness(fs_qty::semantic::HardnessScale::Vickers),
        ] {
            let wrong = PropertyKey::with_quantity(
                "H",
                QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static)),
            );
            assert!(indentation_hardness_property(&bundle, &wrong).is_err());
        }
        let raw_key = PropertyKey::new("declared-hardness", Pressure::DIMS);
        for value in [5e8, 0.0, -1.0] {
            let mut raw = card.claims().claims_for("H")[0].1.clone();
            raw.key = raw_key.clone();
            raw.value = PropertyValue::Scalar {
                value,
                dims: Pressure::DIMS,
            };
            raw.interpolation = InterpolationPolicy::ConstantWithinValidity;
            raw.observations.clear();
            let mut claims = ClaimSet::new();
            claims.insert_claim(raw).unwrap();
            let declared = MaterialCard::assemble(card.id().clone(), claims, vec![]).unwrap();
            let bundle = resolved(&declared, &[raw_key.clone()]);
            assert_eq!(
                indentation_hardness_property(&bundle, &raw_key).is_ok(),
                value > 0.0
            );
            assert!(
                IndentationHardnessStatePoint::try_from_resolved(bundle, raw_key.clone()).is_err()
            );
        }
    }
}
