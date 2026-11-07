//! Binary liquid kinetics from one original, immutable material source.
//!
//! The ordered species assignment explicitly identifies the source's first
//! species MOLE-fraction coordinate. Temperature and composition evolve;
//! every other typed external coordinate, including pressure, stays at its
//! original reference value. No diffusivity is inferred from viscosity, the
//! regular-solution interaction parameter, or a signed Fick coefficient.
//!
//! Every face evaluation qualifies both endpoints and their arithmetic-mean
//! composition against the same pinned claim. The source is evaluated AT the
//! face, not averaged between endpoint diffusivities. This supplies a positive
//! stage-lagged Maxwell–Stefan mobility to a transport owner; it neither solves
//! transport nor certifies source support between the three actual samples.
//! Queries and full source receipts allocate on each call. Bit-identical
//! samples share one resolution within the call and clone its receipt. The
//! original card is shared; no mutable cache or replacement source is introduced.

use super::{BinaryLiquidError, poll};
use crate::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, is_fixed_external_coordinate,
    resolve_material_state_point,
};
use fs_blake3::ContentHash;
use fs_matdb::{ClaimId, MAX_PROPERTY_USAGE_AXIS_BYTES, MaterialCard, PropertyKey, QueryPoint};
use fs_qty::{
    Dims, QuantitySpec, Temperature,
    semantic::{CompositionBasis, QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{SpeciesId, solution::BINARY_COMPOSITION_SUM_TOLERANCE_V1};
use std::sync::Arc;

/// Binary Maxwell–Stefan diffusivity [m²/s], independent of thermodynamic factors.
pub const BINARY_LIQUID_DIFFUSIVITY_DIMS: Dims = Dims([2, 0, -1, 0, 0, 0]);

/// Immutable binary kinetic source with an explicit, ordered species binding.
///
/// Equality compares the full original card identity, pinned claim, property,
/// reference context, axes and species order. Independently reconstructed
/// identical bindings compare equal; changing a held pressure does not.
#[derive(Clone, Debug)]
pub struct MaterialBinaryDiffusivity {
    species: [SpeciesId; 2],
    material: Arc<MaterialCard>,
    card_identity: ContentHash,
    property: PropertyKey,
    selected: ClaimId,
    reference: QueryPoint,
    temperature_axis: String,
    mole_fraction_axis: String,
}

impl PartialEq for MaterialBinaryDiffusivity {
    fn eq(&self, other: &Self) -> bool {
        self.card_identity == other.card_identity
            && self.species == other.species
            && self.property == other.property
            && self.selected == other.selected
            && self.reference == other.reference
            && self.temperature_axis == other.temperature_axis
            && self.mole_fraction_axis == other.mole_fraction_axis
    }
}

/// Original endpoint and face receipts for one actual transport stage.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedBinaryDiffusivity {
    /// Source evaluated at the arithmetic-mean composition [m²/s]. Zero is an
    /// explicit closed kinetic path; negative or nonfinite values are refused.
    pub diffusivity_m2_s: f64,
    pub left: ResolvedMaterialStatePoint,
    pub face: ResolvedMaterialStatePoint,
    pub right: ResolvedMaterialStatePoint,
}

impl MaterialBinaryDiffusivity {
    /// Select the original kinetic claim once at the actual reference state.
    ///
    /// Absolute temperature and first-species mole fraction must have explicit
    /// static quantity declarations. Other coordinates must be typed external
    /// conditions; stress, strain, rates and fields need their own evolving
    /// physical owner and are not silently frozen by this adapter. Source
    /// validity and interpolation support remain authoritative at every call.
    pub fn new(
        species: [SpeciesId; 2],
        material: Arc<MaterialCard>,
        property: PropertyKey,
        selection: MaterialPropertySelection,
        reference: QueryPoint,
        temperature_axis: &str,
        mole_fraction_axis: &str,
    ) -> Result<Self, BinaryLiquidError> {
        if species[0] == species[1]
            || property.quantity() != QuantitySpec::dimensional(BINARY_LIQUID_DIFFUSIVITY_DIMS)
            || property.hardness_test().is_some()
            || property.elastic_component().is_some()
            || property.strain_component().is_some()
            || property.stress_component().is_some()
        {
            return Err(BinaryLiquidError::Invalid(
                "distinct ordered species and an explicit scalar m2/s property required",
            ));
        }
        if !valid_axis(temperature_axis)
            || !valid_axis(mole_fraction_axis)
            || temperature_axis == mole_fraction_axis
            || reference.axis_quantities().get(temperature_axis)
                != Some(&typed(QuantityKind::AbsoluteTemperature))
            || reference.axis_quantities().get(mole_fraction_axis)
                != Some(&typed(QuantityKind::Composition(
                    CompositionBasis::MoleFraction,
                )))
            || !reference
                .axes()
                .get(temperature_axis)
                .is_some_and(|t| t.is_finite() && *t > 0.0)
            || !reference
                .axes()
                .get(mole_fraction_axis)
                .is_some_and(|x| (0.0..=1.0).contains(x))
        {
            return Err(BinaryLiquidError::Invalid(
                "distinct typed absolute-temperature and first-species mole-fraction axes required",
            ));
        }
        if reference.axes().iter().any(|(axis, &value)| {
            let quantity = reference.axis_quantities().get(axis).copied();
            !valid_axis(axis)
                || quantity.and_then(QuantitySpec::semantic_type).is_none()
                || (axis != temperature_axis
                    && quantity.is_some_and(|q| q.dims() == Temperature::DIMS))
                || !is_fixed_external_coordinate(axis, quantity, Some(value))
        }) {
            return Err(BinaryLiquidError::Invalid(
                "diffusion source requires one owned temperature and typed fixed external conditions",
            ));
        }
        let requirement =
            ScalarPropertyRequirement::try_with_key(&property, ScalarAdmissibility::NonNegative)?;
        let resolved = resolve_material_state_point(
            &material,
            &reference,
            core::slice::from_ref(&requirement),
            selection,
        )?;
        let selected = resolved
            .property_by_key(&property)
            .expect("original resolver retains the exact required property")
            .answer()
            .receipt
            .selected;
        let claim = material
            .claims()
            .claim(selected)
            .expect("selected claim belongs to the original immutable card");
        if claim
            .validity
            .bounds()
            .keys()
            .any(|axis| !reference.axis_quantities().contains_key(axis))
            || claim.value.axes().any(|(axis, dims)| {
                reference
                    .axis_quantities()
                    .get(axis)
                    .is_none_or(|quantity| quantity.dims() != dims)
            })
        {
            return Err(BinaryLiquidError::Invalid(
                "binary diffusivity source has an unowned physical coordinate",
            ));
        }
        Ok(Self {
            species,
            material,
            card_identity: resolved.card_identity(),
            property,
            selected,
            reference,
            temperature_axis: temperature_axis.to_owned(),
            mole_fraction_axis: mole_fraction_axis.to_owned(),
        })
    }

    pub const fn species(&self) -> &[SpeciesId; 2] {
        &self.species
    }
    pub fn material(&self) -> &MaterialCard {
        &self.material
    }
    pub const fn property(&self) -> &PropertyKey {
        &self.property
    }
    pub const fn selected_claim(&self) -> ClaimId {
        self.selected
    }
    pub const fn reference_query(&self) -> &QueryPoint {
        &self.reference
    }

    /// Resolve both actual endpoints and the arithmetic-mean face at one T.
    ///
    /// Fractions follow the bound species order and must lie on the binary
    /// simplex within the shared thermochemical rounding allowance. Coordinates
    /// are not clipped or renormalized. A positive second-species trace erased
    /// by the source's single first-species coordinate refuses explicitly.
    /// Cancellation is checked before each synchronous source resolution and
    /// before publication. Failures return no partial result or mutated owner.
    pub fn resolve_face(
        &self,
        temperature: Temperature,
        fractions: [[f64; 2]; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ResolvedBinaryDiffusivity, BinaryLiquidError> {
        poll(&mut cancelled)?;
        if !temperature.value().is_finite() || temperature.value() <= 0.0 {
            return Err(BinaryLiquidError::Invalid(
                "binary diffusivity requires positive finite absolute temperature",
            ));
        }
        for endpoint in fractions {
            validate_fractions(endpoint)?;
        }
        let mean = [
            (fractions[0][0] + fractions[1][0]) * 0.5,
            (fractions[0][1] + fractions[1][1]) * 0.5,
        ];
        if mean
            .iter()
            .enumerate()
            .any(|(i, &value)| (fractions[0][i] > 0.0 || fractions[1][i] > 0.0) && value == 0.0)
        {
            return Err(BinaryLiquidError::Arithmetic(
                "positive component lost in arithmetic face composition",
            ));
        }
        validate_fractions(mean)?;
        let requirement = ScalarPropertyRequirement::try_with_key(
            &self.property,
            ScalarAdmissibility::NonNegative,
        )?;
        let mut at = |x: f64,
                      existing: Option<&ResolvedMaterialStatePoint>|
         -> Result<ResolvedMaterialStatePoint, BinaryLiquidError> {
            poll(&mut cancelled)?;
            if let Some(existing) = existing {
                return Ok(existing.clone());
            }
            let point = self
                .reference
                .clone()
                .with(&self.temperature_axis, temperature.value())
                .and_then(|point| point.with(&self.mole_fraction_axis, x))
                .map_err(|source| MaterialStatePointError::Query {
                    property: self.property.name().to_owned(),
                    source,
                })?;
            Ok(resolve_material_state_point(
                &self.material,
                &point,
                core::slice::from_ref(&requirement),
                MaterialPropertySelection::PinnedByProperty(vec![(
                    self.property.name().to_owned(),
                    self.selected,
                )]),
            )?)
        };
        let left = at(fractions[0][0], None)?;
        let face = at(
            mean[0],
            (mean[0].to_bits() == fractions[0][0].to_bits()).then_some(&left),
        )?;
        let right = at(
            fractions[1][0],
            if fractions[1][0].to_bits() == fractions[0][0].to_bits() {
                Some(&left)
            } else if fractions[1][0].to_bits() == mean[0].to_bits() {
                Some(&face)
            } else {
                None
            },
        )?;
        let diffusivity_m2_s = face
            .property_by_key(&self.property)
            .expect("original resolver retains the exact required property")
            .value_si();
        poll(&mut cancelled)?;
        Ok(ResolvedBinaryDiffusivity {
            diffusivity_m2_s,
            left,
            face,
            right,
        })
    }
}

fn valid_axis(axis: &str) -> bool {
    !axis.trim().is_empty()
        && axis.len() <= MAX_PROPERTY_USAGE_AXIS_BYTES
        && !axis.chars().any(char::is_control)
}

fn typed(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

fn validate_fractions(fractions: [f64; 2]) -> Result<(), BinaryLiquidError> {
    if fractions
        .iter()
        .any(|x| !x.is_finite() || !(0.0..=1.0).contains(x))
        || (fractions[0] + fractions[1] - 1.0).abs() > BINARY_COMPOSITION_SUM_TOLERANCE_V1
    {
        return Err(BinaryLiquidError::Invalid(
            "binary diffusivity endpoints must be normalized nonnegative mole fractions",
        ));
    }
    if fractions[0] == 1.0 && fractions[1] > 0.0 {
        return Err(BinaryLiquidError::Arithmetic(
            "nonzero component lost in the source's first-species coordinate",
        ));
    }
    Ok(())
}
