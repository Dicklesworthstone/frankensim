//! Wood moisture: the equilibrium moisture content a declared environment
//! gives, and the clear-wood stiffness and density at that moisture content,
//! in the printed forms of the FPL Wood Handbook (FPL-GTR-282). Every
//! coefficient is a claim of a material card; this module holds only the forms.
//!
//! - Equation (4-5), the Hailwood-Horrobin EMC fitted by Simpson (1973), with
//!   W, K, K1 and K2 quadratic in the Celsius temperature
//!   (pack `wood-emc-hailwood-horrobin-fpl-gtr282`, valid over Table 4-2).
//! - Equation (5-3), P = P12 (P12/Pg)^((12 - M)/(Mp - 12)) below Mp, applied
//!   with the bending modulus's 12 % and green values.
//! - Equation (4-9), shrinkage S_x = S_0 (1 - x/MCfs) below the fiber
//!   saturation point, which with the (1 + M) moisture mass gives density.
//!
//! Moisture contents are fractions (0.12 is 12 %), as the pack compiler
//! stores a `%` claim. The elastic RATIOS of a reference-ratio card are not
//! moisture laws in the source: a caller that scales every modulus by
//! [`WoodMoistureLaw::modulus_factor`] holds them at their 12 % values, an
//! assumption it must state.
use crate::MaterialError;
use crate::state_point::{
    MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, is_fixed_external_coordinate, resolve_material_state_point,
};
use fs_blake3::ContentHash;
use fs_matdb::{EvaluationDecision, MaterialCard, PropertyKey, QueryPoint, SelectionPolicy};
use fs_math::det;
use fs_qty::semantic::{MoistureBasis, QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, QuantitySpec};

/// The FPL tables' reference moisture content, 12 %.
pub const REFERENCE_MOISTURE_CONTENT: f64 = 0.12;
const PER_KELVIN: Dims = Dims([0, 0, 0, -1, 0, 0]);
const PER_KELVIN_SQUARED: Dims = Dims([0, 0, 0, -2, 0, 0]);
const PRESSURE: Dims = Dims([-1, 1, -2, 0, 0, 0]);

fn claim(
    card: &MaterialCard,
    name: &str,
    point: &QueryPoint,
    dims: Dims,
) -> Result<f64, MaterialError> {
    let answer = card
        .claims()
        .query(name, point, SelectionPolicy::SingleClaimOnly)
        .map_err(|e| MaterialError::Parameters {
            what: format!("{name}: {e:?}"),
        })?;
    let sample = &answer.evidence.value;
    if sample.dims != dims || !sample.value.is_finite() {
        return Err(MaterialError::Parameters {
            what: format!("{name}: dims {:?}, value {}", sample.dims.0, sample.value),
        });
    }
    Ok(sample.value)
}

/// Equation (4-5): the equilibrium moisture content (fraction) of wood at
/// `temperature_k` and `relative_humidity` (0..1), from a card carrying the
/// Hailwood-Horrobin coefficients. Outside the card's validity it refuses.
///
/// # Errors
/// [`MaterialError::Parameters`] for a missing or mis-dimensioned
/// coefficient, a state outside the card's validity, or a nonphysical result.
pub fn equilibrium_moisture_content(
    card: &MaterialCard,
    temperature_k: f64,
    relative_humidity: f64,
) -> Result<f64, MaterialError> {
    let point = QueryPoint::new()
        .with("temperature", temperature_k)
        .and_then(|p| p.with("relative_humidity", relative_humidity))
        .map_err(|e| MaterialError::Parameters {
            what: format!("EMC state: {e:?}"),
        })?;
    let quadratic = |name: &str| -> Result<[f64; 3], MaterialError> {
        let stem = format!("wood_emc_hailwood_horrobin_{name}");
        Ok([
            claim(card, &format!("{stem}_constant"), &point, Dims::NONE)?,
            claim(card, &format!("{stem}_per_celsius"), &point, PER_KELVIN)?,
            claim(
                card,
                &format!("{stem}_per_celsius_squared"),
                &point,
                PER_KELVIN_SQUARED,
            )?,
        ])
    };
    let numerator = claim(
        card,
        "wood_emc_hailwood_horrobin_numerator",
        &point,
        Dims::NONE,
    )?;
    // The printed polynomials take the Celsius temperature.
    let t = temperature_k - 273.15;
    let at = |c: [f64; 3]| c[0] + c[1] * t + c[2] * t * t;
    let (w, k, k1, k2) = (
        at(quadratic("w")?),
        at(quadratic("k")?),
        at(quadratic("k1")?),
        at(quadratic("k2")?),
    );
    let kh = k * relative_humidity;
    let percent = numerator / w
        * (kh / (1.0 - kh)
            + (k1 * kh + 2.0 * k1 * k2 * kh * kh) / (1.0 + k1 * kh + k1 * k2 * kh * kh));
    if !(percent.is_finite() && percent >= 0.0 && kh < 1.0) {
        return Err(MaterialError::Parameters {
            what: format!(
                "EMC at {temperature_k} K, RH {relative_humidity} is not physical ({percent} %)"
            ),
        });
    }
    Ok(percent / 100.0)
}

/// One species' moisture laws, from its card's claims at a moisture content.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WoodMoistureLaw {
    /// Bending modulus at 12 % [Pa], the P12 of Equation (5-3).
    pub bending_modulus_12_pa: f64,
    /// Green bending modulus [Pa], the Pg of Equation (5-3).
    pub bending_modulus_green_pa: f64,
    /// Mp, where the green and dry branches of Equation (5-3) meet (fraction).
    pub intersection_moisture_content: f64,
    /// S_0, volumetric shrinkage green to ovendry (fraction of green volume).
    pub volumetric_shrinkage: f64,
    /// MCfs, the fiber saturation point of Equation (4-9) (fraction).
    pub fiber_saturation: f64,
    /// The moisture content the laws were admitted at (fraction).
    pub moisture_content: f64,
}

impl WoodMoistureLaw {
    /// The laws of `card` at `moisture_content`. The 12 % bending modulus is
    /// read at 12 %; the law parameters at `moisture_content`, so a moisture
    /// content outside their stated domain (below Mp, below MCfs) refuses.
    ///
    /// # Errors
    /// [`MaterialError::Parameters`] for a missing claim, a moisture content
    /// outside a claim's validity, or inconsistent parameters.
    pub fn from_card(card: &MaterialCard, moisture_content: f64) -> Result<Self, MaterialError> {
        let at = |m: f64| {
            QueryPoint::new()
                .with("moisture_content", m)
                .map_err(|e| MaterialError::Parameters {
                    what: format!("moisture state: {e:?}"),
                })
        };
        let (reference, here) = (at(REFERENCE_MOISTURE_CONTENT)?, at(moisture_content)?);
        let law = Self {
            bending_modulus_12_pa: claim(
                card,
                "modulus_of_elasticity_bending",
                &reference,
                PRESSURE,
            )?,
            bending_modulus_green_pa: claim(
                card,
                "moisture_law_green_bending_modulus",
                &here,
                PRESSURE,
            )?,
            intersection_moisture_content: claim(
                card,
                "moisture_law_intersection_moisture_content",
                &here,
                Dims::NONE,
            )?,
            volumetric_shrinkage: claim(
                card,
                "volumetric_shrinkage_green_to_ovendry",
                &here,
                Dims::NONE,
            )?,
            fiber_saturation: claim(card, "fiber_saturation_moisture_content", &here, Dims::NONE)?,
            moisture_content,
        };
        let consistent = law.bending_modulus_green_pa > 0.0
            && law.bending_modulus_12_pa > 0.0
            && law.intersection_moisture_content > REFERENCE_MOISTURE_CONTENT
            && (0.0..1.0).contains(&law.volumetric_shrinkage)
            && law.fiber_saturation > REFERENCE_MOISTURE_CONTENT
            && (0.0..law.intersection_moisture_content.min(law.fiber_saturation))
                .contains(&moisture_content);
        if !consistent {
            return Err(MaterialError::Parameters {
                what: format!("inconsistent wood moisture laws {law:?}"),
            });
        }
        Ok(law)
    }

    /// Equation (5-3) as a factor on a 12 % property: (P12/Pg)^((12 - M)/(Mp - 12)).
    #[must_use]
    pub fn modulus_factor(&self) -> f64 {
        let exponent = (REFERENCE_MOISTURE_CONTENT - self.moisture_content)
            / (self.intersection_moisture_content - REFERENCE_MOISTURE_CONTENT);
        det::pow(
            self.bending_modulus_12_pa / self.bending_modulus_green_pa,
            exponent,
        )
    }

    /// Density at M over density at 12 %: moisture mass (1 + M) over the
    /// Equation (4-9) volume, V(M) = V_green (1 - S_0 (1 - M/MCfs)).
    #[must_use]
    pub fn density_factor(&self) -> f64 {
        let volume = |m| remaining_fraction(self.volumetric_shrinkage, m, self.fiber_saturation);
        (1.0 + self.moisture_content) / (1.0 + REFERENCE_MOISTURE_CONTENT)
            * volume(REFERENCE_MOISTURE_CONTENT)
            / volume(self.moisture_content)
    }
}

fn remaining_fraction(shrinkage: f64, moisture: f64, fiber_saturation: f64) -> f64 {
    1.0 - shrinkage * (1.0 - moisture / fiber_saturation)
}

/// Original source roles for the FPL Equation (4-9) dimensional model.
/// All moisture values are water mass divided by ovendry wood mass. An
/// unqualified dimensional coordinate is accepted only with this explicit
/// dry-mass convention; wet-mass moisture and relative humidity are different.
#[derive(Clone, Copy, Debug)]
pub struct WoodMoistureSource<'a> {
    /// Density at exactly 12% dry-basis moisture, in kg/m³.
    pub density: &'a PropertyKey,
    /// Radial green-to-ovendry shrinkage as a fraction of green length.
    pub radial_shrinkage: &'a PropertyKey,
    /// Tangential green-to-ovendry shrinkage as a fraction of green length.
    pub tangential_shrinkage: &'a PropertyKey,
    /// Independent volumetric green-to-ovendry shrinkage of green volume.
    pub volumetric_shrinkage: &'a PropertyKey,
    /// Source-declared fiber saturation; a handbook approximation remains one.
    pub fiber_saturation: &'a PropertyKey,
    /// Original source coordinate name, without aliasing or normalization.
    pub moisture_axis: &'a str,
    /// Complete original query for density at 12% moisture.
    pub reference_point: &'a QueryPoint,
    /// Complete original query for the four constant model parameters.
    pub shrinkage_point: &'a QueryPoint,
}

/// Source-bound moisture dimensions and density relative to 12% moisture.
/// The independent radial, tangential and volume averages do not define one
/// deformation tensor. In particular, no longitudinal change is inferred.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedWoodMoisture {
    density_source: ResolvedMaterialStatePoint,
    shrinkage_source: ResolvedMaterialStatePoint,
    reference_density: f64,
    shrinkages: [f64; 3],
    fiber_saturation: f64,
    moisture_span: [f64; 2],
}

/// A prescribed-moisture state of the retained handbook model. Moisture
/// exchange changes total mass; density times volume divided by `(1 + M)`
/// preserves ovendry mass. This does not predict drying rate or moisture flow.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WoodMoistureState {
    /// Immutable card providing density and all four model coefficients.
    pub material_card_identity: ContentHash,
    /// Prescribed water mass / ovendry wood mass.
    pub moisture_content: f64,
    /// Actual mass, including the prescribed water, per current volume [kg/m³].
    pub density_kg_m3: f64,
    /// Radial length at this moisture divided by radial length at 12%.
    pub radial_length_factor: f64,
    /// Tangential length at this moisture divided by tangential length at 12%.
    pub tangential_length_factor: f64,
    /// Volume at this moisture divided by volume at 12%, using its own source.
    pub volume_factor: f64,
}

impl ResolvedWoodMoisture {
    /// Bind exact source keys and contexts without freezing varying laws or
    /// replacing missing coefficients. Source applicability bounds intersect
    /// the below-fiber-saturation model domain; they are never widened.
    #[allow(clippy::too_many_lines)]
    pub fn from_card(
        card: &MaterialCard,
        source: WoodMoistureSource<'_>,
        selection: MaterialPropertySelection,
    ) -> Result<Self, MaterialError> {
        let dry_mass = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::Moisture(MoistureBasis::DryMass),
            ValueForm::Static,
        ));
        for point in [source.reference_point, source.shrinkage_point] {
            if !point.axes().contains_key(source.moisture_axis)
                || !point
                    .axis_quantities()
                    .get(source.moisture_axis)
                    .is_some_and(|q| *q == QuantitySpec::dimensional(Dims::NONE) || *q == dry_mass)
            {
                return Err(wood_error(
                    "a typed dry-basis moisture source coordinate is required",
                ));
            }
        }
        if source.reference_point.axes()[source.moisture_axis] != REFERENCE_MOISTURE_CONTENT {
            return Err(wood_error(
                "density must be queried at exactly 12% dry-basis moisture",
            ));
        }
        for (axis, value) in source.reference_point.axes() {
            if axis != source.moisture_axis
                && let Some(other) = source.shrinkage_point.axes().get(axis)
                && (value.to_bits() != other.to_bits()
                    || source.reference_point.axis_quantities().get(axis)
                        != source.shrinkage_point.axis_quantities().get(axis))
            {
                return Err(wood_error(
                    "density and shrinkage sources must agree on shared non-moisture conditions",
                ));
            }
        }
        let requirement = |key: &PropertyKey, dims, admissibility| {
            if key.quantity() != QuantitySpec::dimensional(dims) {
                return Err(wood_error(
                    "wood model keys must retain their original SI quantity dimensions",
                ));
            }
            ScalarPropertyRequirement::try_with_key(key, admissibility)
                .map_err(|error| wood_error(error.to_string()))
        };
        let (density_selection, shrinkage_selection) = match selection {
            MaterialPropertySelection::PinnedByProperty(pins) => {
                if pins.len() != 5 {
                    return Err(wood_error(
                        "the source plan must pin exactly five wood parameters",
                    ));
                }
                let (density_pins, shrinkage_pins): (Vec<_>, Vec<_>) =
                    pins.into_iter().partition(|(name, id)| {
                        name == source.density.name()
                            && card
                                .claims()
                                .claim(*id)
                                .is_some_and(|claim| claim.key == *source.density)
                    });
                (
                    MaterialPropertySelection::PinnedByProperty(density_pins),
                    MaterialPropertySelection::PinnedByProperty(shrinkage_pins),
                )
            }
            policy => (policy.clone(), policy),
        };
        let density_source = resolve_material_state_point(
            card,
            source.reference_point,
            &[requirement(
                source.density,
                Density::DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )?],
            density_selection,
        )
        .map_err(|error| wood_error(error.to_string()))?;
        let keys = [
            source.radial_shrinkage,
            source.tangential_shrinkage,
            source.volumetric_shrinkage,
            source.fiber_saturation,
        ];
        let requirements = keys
            .iter()
            .map(|key| requirement(key, Dims::NONE, ScalarAdmissibility::NonNegative))
            .collect::<Result<Vec<_>, _>>()?;
        let shrinkage_source = resolve_material_state_point(
            card,
            source.shrinkage_point,
            &requirements,
            shrinkage_selection,
        )
        .map_err(|error| wood_error(error.to_string()))?;
        for (state, point) in [
            (&density_source, source.reference_point),
            (&shrinkage_source, source.shrinkage_point),
        ] {
            admit_wood_parameters(state, point, source.moisture_axis)?;
        }
        let density = &density_source.properties()[0];
        if density
            .answer()
            .evidence
            .model
            .validity
            .bound(source.moisture_axis)
            != Some((REFERENCE_MOISTURE_CONTENT, REFERENCE_MOISTURE_CONTENT))
        {
            return Err(wood_error(
                "the reference density source must be specific to 12% moisture",
            ));
        }
        let values = keys.map(|key| {
            shrinkage_source
                .property_by_key(key)
                .expect("the exact required key was resolved")
                .value_si()
        });
        let [radial, tangential, volumetric, fiber_saturation] = values;
        if [radial, tangential, volumetric]
            .iter()
            .any(|s| !(0.0..1.0).contains(s))
            || !(REFERENCE_MOISTURE_CONTENT..=1.0).contains(&fiber_saturation)
            || fiber_saturation == REFERENCE_MOISTURE_CONTENT
        {
            return Err(wood_error(
                "shrinkages must be fractions below one and fiber saturation above 12%",
            ));
        }
        let mut moisture_span = [0.0_f64, fiber_saturation];
        for property in shrinkage_source.properties() {
            let (low, high) = property
                .answer()
                .evidence
                .model
                .validity
                .bound(source.moisture_axis)
                .ok_or_else(|| {
                    wood_error("model parameters require an explicit moisture applicability domain")
                })?;
            moisture_span[0] = moisture_span[0].max(low);
            moisture_span[1] = moisture_span[1].min(high);
        }
        if !(moisture_span[0]..=moisture_span[1]).contains(&REFERENCE_MOISTURE_CONTENT) {
            return Err(wood_error(
                "model parameter domains must include the 12% reference",
            ));
        }
        let reference_density = density.value_si();
        Ok(Self {
            density_source,
            shrinkage_source,
            reference_density,
            shrinkages: [radial, tangential, volumetric],
            fiber_saturation,
            moisture_span,
        })
    }

    /// Original reference-density query, card and selected-source receipt.
    #[must_use]
    pub const fn density_source(&self) -> &ResolvedMaterialStatePoint {
        &self.density_source
    }

    /// Original R/T/V/FSP parameter query and all four selected-source receipts.
    #[must_use]
    pub const fn shrinkage_source(&self) -> &ResolvedMaterialStatePoint {
        &self.shrinkage_source
    }

    /// Source-provided dry-basis fiber-saturation approximation.
    #[must_use]
    pub const fn fiber_saturation(&self) -> f64 {
        self.fiber_saturation
    }

    /// Intersection of the parameter domains and `[0, fiber_saturation]`.
    #[must_use]
    pub const fn moisture_span(&self) -> [f64; 2] {
        self.moisture_span
    }

    /// Evaluate Equation (4-9) at prescribed dry-basis moisture. Ratios of
    /// remaining dimensions give factors relative to the 12% reference, not
    /// subtraction of two shrinkage percentages. Volume uses its own data.
    pub fn at_moisture_content(
        &self,
        moisture_content: f64,
    ) -> Result<WoodMoistureState, MaterialError> {
        if !moisture_content.is_finite()
            || !(self.moisture_span[0]..=self.moisture_span[1]).contains(&moisture_content)
        {
            return Err(wood_error(
                "moisture is outside the admitted below-fiber-saturation domain",
            ));
        }
        let [
            radial_length_factor,
            tangential_length_factor,
            volume_factor,
        ] = self.shrinkages.map(|s| {
            remaining_fraction(s, moisture_content, self.fiber_saturation)
                / remaining_fraction(s, REFERENCE_MOISTURE_CONTENT, self.fiber_saturation)
        });
        let density_kg_m3 = self.reference_density
            * ((1.0 + moisture_content) / (1.0 + REFERENCE_MOISTURE_CONTENT))
            / volume_factor;
        if !density_kg_m3.is_finite() || density_kg_m3 <= 0.0 {
            return Err(wood_error(
                "moisture-adjusted density is not finite and positive",
            ));
        }
        Ok(WoodMoistureState {
            material_card_identity: self.density_source.card_identity(),
            moisture_content,
            density_kg_m3,
            radial_length_factor,
            tangential_length_factor,
            volume_factor,
        })
    }
}

fn wood_error(what: impl Into<String>) -> MaterialError {
    MaterialError::Parameters { what: what.into() }
}

fn admit_wood_parameters(
    state: &ResolvedMaterialStatePoint,
    point: &QueryPoint,
    moisture_axis: &str,
) -> Result<(), MaterialError> {
    for property in state.properties() {
        if !property.source_axes().is_empty()
            || property.answer().receipt.decision != EvaluationDecision::ConstantWithinValidity
        {
            return Err(wood_error(
                "wood reference parameters must be static scalar source laws",
            ));
        }
        for (axis, &(low, high)) in property.answer().evidence.model.validity.bounds() {
            if axis == moisture_axis {
                continue;
            }
            let quantity = point.axis_quantities().get(axis).copied();
            let value = point.axes().get(axis).copied();
            let flag = quantity == Some(QuantitySpec::dimensional(Dims::NONE))
                && low == high
                && (low == 0.0 || low == 1.0)
                && value == Some(low);
            if !flag && !is_fixed_external_coordinate(axis, quantity, value) {
                return Err(wood_error(
                    "wood parameters require fixed typed conditions and explicit source applicability flags",
                ));
            }
        }
    }
    Ok(())
}

/// Equation (5-3) on explicit values (moisture contents as fractions), for a
/// property with no card: P12 (P12/Pg)^((0.12 - M)/(Mp - 0.12)).
#[must_use]
pub fn moisture_adjusted_property(
    p12: f64,
    green: f64,
    intersection: f64,
    moisture_content: f64,
) -> f64 {
    p12 * det::pow(
        p12 / green,
        (REFERENCE_MOISTURE_CONTENT - moisture_content)
            / (intersection - REFERENCE_MOISTURE_CONTENT),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        ClaimSet, InterpolationPolicy, MaterialStateId, PropertyClaim, PropertyValue, Provenance,
        UncertaintyModel,
    };

    fn fixture(
        axis: &str,
        values: [f64; 5],
        span: [f64; 2],
    ) -> (MaterialCard, [PropertyKey; 5], QueryPoint) {
        let keys = [
            "rho-source",
            "radial-source",
            "tangential-source",
            "volume-source",
            "fsp-source",
        ]
        .map(|name| {
            PropertyKey::new(
                name,
                if name == "rho-source" {
                    Density::DIMS
                } else {
                    Dims::NONE
                },
            )
        });
        let point = QueryPoint::new()
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), 0.12)
            .unwrap()
            .with_quantity(
                "engineering-reference-model",
                QuantitySpec::dimensional(Dims::NONE),
                1.0,
            )
            .unwrap();
        let mut claims = ClaimSet::new();
        for (index, (key, value)) in keys.iter().zip(values).enumerate() {
            let bounds = if index == 0 { [0.12, 0.12] } else { span };
            claims
                .insert_claim(PropertyClaim {
                    key: key.clone(),
                    value: PropertyValue::Scalar {
                        value,
                        dims: key.dims(),
                    },
                    validity: ValidityDomain::default()
                        .with_quantity(
                            axis,
                            QuantitySpec::dimensional(Dims::NONE),
                            bounds[0],
                            bounds[1],
                        )
                        .with_quantity(
                            "engineering-reference-model",
                            QuantitySpec::dimensional(Dims::NONE),
                            1.0,
                            1.0,
                        ),
                    uncertainty: UncertaintyModel::Unstated,
                    interpolation: InterpolationPolicy::ConstantWithinValidity,
                    observations: Vec::new(),
                    provenance: Provenance {
                        source: "manufactured independent moisture geometry coefficients".into(),
                        license: "CC0-1.0".into(),
                        artifact: None,
                    },
                })
                .unwrap();
        }
        let card = MaterialCard::assemble(
            MaterialStateId {
                chemistry: "manufactured moisture model".into(),
                phase: "solid".into(),
                process: "numerical regression only".into(),
                revision: 0,
            },
            claims,
            Vec::new(),
        )
        .unwrap();
        (card, keys, point)
    }

    fn source<'a>(
        keys: &'a [PropertyKey; 5],
        point: &'a QueryPoint,
        axis: &'a str,
    ) -> WoodMoistureSource<'a> {
        WoodMoistureSource {
            density: &keys[0],
            radial_shrinkage: &keys[1],
            tangential_shrinkage: &keys[2],
            volumetric_shrinkage: &keys[3],
            fiber_saturation: &keys[4],
            moisture_axis: axis,
            reference_point: point,
            shrinkage_point: point,
        }
    }

    #[test]
    fn sourced_moisture_geometry_preserves_dry_mass_and_independent_volume() {
        for axis in ["moisture-content", "moisture_content"] {
            let (card, keys, point) = fixture(axis, [560.0, 0.04, 0.08, 0.14, 0.30], [0.0, 0.30]);
            let law = ResolvedWoodMoisture::from_card(
                &card,
                source(&keys, &point, axis),
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap();
            let reference = law.at_moisture_content(0.12).unwrap();
            assert_eq!(reference.density_kg_m3, 560.0);
            assert_eq!(reference.radial_length_factor, 1.0);
            assert_eq!(reference.tangential_length_factor, 1.0);
            assert_eq!(reference.volume_factor, 1.0);
            let dry = law.at_moisture_content(0.0).unwrap();
            // Independent green dimensions: R12=.976, T12=.952, V12=.916.
            assert!((dry.radial_length_factor - 0.96 / 0.976).abs() < 1e-14);
            assert!((dry.tangential_length_factor - 0.92 / 0.952).abs() < 1e-14);
            assert!((dry.volume_factor - 0.86 / 0.916).abs() < 1e-14);
            assert!((dry.density_kg_m3 - 500.0 * 0.916 / 0.86).abs() < 1e-10);
            assert!(
                (dry.volume_factor - dry.radial_length_factor * dry.tangential_length_factor).abs()
                    > 0.01
            );
            for moisture in [0.0, 0.08, 0.12, 0.20, 0.30] {
                let state = law.at_moisture_content(moisture).unwrap();
                assert!(
                    (state.density_kg_m3 * state.volume_factor / (1.0 + moisture) - 500.0).abs()
                        < 1e-10
                );
                assert_eq!(state.material_card_identity, card.content_hash());
            }
            for state in [law.density_source(), law.shrinkage_source()] {
                assert_eq!(state.axis_quantities(), point.axis_quantities());
                for property in state.properties() {
                    card.claims()
                        .verify_receipt(&property.answer().receipt)
                        .unwrap();
                }
            }
            let pins = [law.density_source(), law.shrinkage_source()]
                .into_iter()
                .flat_map(|state| state.properties())
                .map(|property| {
                    (
                        property.requirement().name().to_owned(),
                        property.answer().receipt.selected,
                    )
                })
                .collect();
            let pinned = ResolvedWoodMoisture::from_card(
                &card,
                source(&keys, &point, axis),
                MaterialPropertySelection::PinnedByProperty(pins),
            )
            .unwrap();
            assert_eq!(
                pinned.at_moisture_content(0.20).unwrap(),
                law.at_moisture_content(0.20).unwrap()
            );
        }
    }

    #[test]
    fn sourced_moisture_refuses_changed_basis_and_unsupported_domains() {
        let axis = "moisture-content";
        let (card, keys, point) = fixture(axis, [560.0, 0.04, 0.08, 0.14, 0.30], [0.05, 0.20]);
        let law = ResolvedWoodMoisture::from_card(
            &card,
            source(&keys, &point, axis),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        assert_eq!(law.moisture_span(), [0.05, 0.20]);
        for moisture in [0.0, 0.049, 0.201, 0.30, f64::NAN, f64::INFINITY] {
            assert!(law.at_moisture_content(moisture).is_err());
        }
        let wet_mass = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::Moisture(MoistureBasis::WetMass),
            ValueForm::Static,
        ));
        for invalid in [
            point.clone().with(axis, 0.08).unwrap(),
            point.clone().with_quantity(axis, wet_mass, 0.12).unwrap(),
            point
                .clone()
                .with("engineering-reference-model", 0.0)
                .unwrap(),
        ] {
            assert!(
                ResolvedWoodMoisture::from_card(
                    &card,
                    source(&keys, &invalid, axis),
                    MaterialPropertySelection::SingleClaimOnly
                )
                .is_err()
            );
        }
        let pressure =
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static));
        let reference = point
            .clone()
            .with_quantity("pressure", pressure, 100_000.0)
            .unwrap();
        for conflicting in [
            reference.clone().with("pressure", 200_000.0).unwrap(),
            reference
                .clone()
                .with_quantity("pressure", QuantitySpec::dimensional(PRESSURE), 100_000.0)
                .unwrap(),
        ] {
            assert!(
                ResolvedWoodMoisture::from_card(
                    &card,
                    WoodMoistureSource {
                        reference_point: &reference,
                        shrinkage_point: &conflicting,
                        ..source(&keys, &point, axis)
                    },
                    MaterialPropertySelection::SingleClaimOnly
                )
                .is_err()
            );
        }
        for values in [
            [560.0, 0.04, 0.08, 1.0, 0.30],
            [560.0, 0.04, 0.08, 0.14, 0.12],
        ] {
            let (bad, keys, point) = fixture(axis, values, [0.0, 0.30]);
            assert!(
                ResolvedWoodMoisture::from_card(
                    &bad,
                    source(&keys, &point, axis),
                    MaterialPropertySelection::SingleClaimOnly
                )
                .is_err()
            );
        }
    }

    /// FPL's own worked example under Equation (5-3): white ash modulus of
    /// rupture at 8 %, P12 = 103,000 kPa, Pg = 66,000 kPa, Mp = 24 %:
    /// "P8 = 103,000 (103,000/66,000)^(4/12) = 119,500 kPa".
    #[test]
    fn equation_5_3_reproduces_the_handbook_example() {
        let p8 = moisture_adjusted_property(103_000.0, 66_000.0, 0.24, 0.08);
        assert!((p8 - 119_500.0).abs() < 50.0, "{p8}");
        assert_eq!(
            moisture_adjusted_property(103_000.0, 66_000.0, 0.24, 0.12),
            103_000.0
        );
    }
}
