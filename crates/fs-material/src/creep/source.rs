//! Immutable material-card binding for the native Norton/Arrhenius creep law.
//!
//! All seven coefficients come from explicit source claims. Their complete
//! temperature line must be supported and constant: the declared Arrhenius law
//! owns the temperature dependence of kinetics, while elastic storage is fixed.
//! This permits temperature changes without silently discarding coefficient
//! work. Temperature-dependent elasticity, primary/tertiary creep, damage and
//! fitted-data qualification require other explicit models.

use super::{CreepError, CreepState, NortonCreep, NortonCreepParameters};
use crate::state_point::{
    MaterialPropertySelection, MaterialStatePointError, POISSON_RATIO_PROPERTY,
    ResolvedMaterialStatePoint, ScalarAdmissibility, ScalarPropertyRequirement,
    YOUNG_MODULUS_PROPERTY, resolve_material_state_point,
};
use crate::{IsotropicElastic, MaterialAdmissibility, SmallStrainLaw, Tangent6, Voigt};
use fs_evidence::ModelCard;
use fs_matdb::{ClaimSelection, MaterialCard, QueryPoint};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};

/// Reference equivalent stress [Pa], never an inferred yield threshold.
pub const CREEP_REFERENCE_STRESS_PROPERTY: &str = "creep_reference_stress";
/// Equivalent strain rate at the reference stress and temperature [1/s].
pub const CREEP_REFERENCE_RATE_PROPERTY: &str = "creep_reference_strain_rate";
/// Dimensionless stress exponent of the secondary-creep power law.
pub const CREEP_STRESS_EXPONENT_PROPERTY: &str = "creep_stress_exponent";
/// Molar activation energy [J/mol], not energy per atom or per mass.
pub const CREEP_ACTIVATION_ENERGY_PROPERTY: &str = "creep_activation_energy";
/// Absolute reference temperature [K].
pub const CREEP_REFERENCE_TEMPERATURE_PROPERTY: &str = "creep_reference_temperature";
pub const CREEP_RATE_DIMS: Dims = Dims([0, 0, -1, 0, 0, 0]);
pub const CREEP_ACTIVATION_ENERGY_DIMS: Dims = Dims([2, 1, -2, 0, 0, -1]);

/// Exact source requirements; no yield stress, hardening or missing-data default.
pub fn norton_requirements() -> Result<Vec<ScalarPropertyRequirement>, MaterialStatePointError> {
    use ScalarAdmissibility::{NonNegative, StrictlyPositive};
    [
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, StrictlyPositive),
        (
            POISSON_RATIO_PROPERTY,
            Dims([0; 6]),
            ScalarAdmissibility::OpenInterval {
                lower: -1.0,
                upper: 0.5,
            },
        ),
        (
            CREEP_REFERENCE_STRESS_PROPERTY,
            Pressure::DIMS,
            StrictlyPositive,
        ),
        (
            CREEP_REFERENCE_RATE_PROPERTY,
            CREEP_RATE_DIMS,
            StrictlyPositive,
        ),
        (
            CREEP_STRESS_EXPONENT_PROPERTY,
            Dims([0; 6]),
            StrictlyPositive,
        ),
        (
            CREEP_ACTIVATION_ENERGY_PROPERTY,
            CREEP_ACTIVATION_ENERGY_DIMS,
            NonNegative,
        ),
        (
            CREEP_REFERENCE_TEMPERATURE_PROPERTY,
            Temperature::DIMS,
            StrictlyPositive,
        ),
    ]
    .into_iter()
    .map(|(name, dims, domain)| ScalarPropertyRequirement::try_new(name, dims, domain))
    .collect()
}

/// One explicit temperature interval at otherwise unchanged source conditions.
#[derive(Clone, Debug)]
pub struct NortonSourceQuery<'a> {
    pub lower: &'a QueryPoint,
    pub upper: &'a QueryPoint,
    /// Explicit absolute-temperature coordinate [K]. Preserve its source schema:
    /// legacy scalar applicability limits may omit an axis descriptor; any
    /// supplied descriptor must denote absolute temperature or kelvin dimensions.
    pub temperature_axis: &'a str,
    pub selection: MaterialPropertySelection,
    /// Maximum admitted tensor strain magnitude, in (0, 0.05].
    pub strain_limit: f64,
}

/// Source admission preserves the original refusal instead of inventing data.
#[derive(Clone, Debug, PartialEq)]
pub enum CreepSourceError {
    Material(MaterialStatePointError),
    Mechanical(CreepError),
    Invalid(&'static str),
    Support { property: String, reason: String },
    Cancelled,
}
impl core::fmt::Display for CreepSourceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Material(e) => write!(f, "creep source: {e}"),
            Self::Mechanical(e) => write!(f, "creep law: {e}"),
            Self::Invalid(s) => write!(f, "creep source: {s}"),
            Self::Support { property, reason } => write!(f, "creep source {property}: {reason}"),
            Self::Cancelled => f.write_str("creep source: cancelled"),
        }
    }
}
impl std::error::Error for CreepSourceError {}
impl From<MaterialStatePointError> for CreepSourceError {
    fn from(value: MaterialStatePointError) -> Self {
        Self::Material(value)
    }
}
impl From<CreepError> for CreepSourceError {
    fn from(value: CreepError) -> Self {
        Self::Mechanical(value)
    }
}

/// Supported, source-pinned family of physical steps. Preparation checks the
/// complete original temperature lines once, outside constitutive hot loops.
#[derive(Clone, Debug)]
pub struct NortonCreepSource<'a> {
    card: &'a MaterialCard,
    reference: ResolvedMaterialStatePoint,
    elastic: IsotropicElastic,
    parameters: NortonCreepParameters,
    temperature_axis: String,
}
impl<'a> NortonCreepSource<'a> {
    /// Admit the complete interval, with unchanged nonthermal coordinates and
    /// original claim identities. Every original knot is checked; matching
    /// endpoints alone cannot conceal an interior modulus or rate change.
    pub fn new(
        card: &'a MaterialCard,
        query: NortonSourceQuery<'_>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, CreepSourceError> {
        poll(&mut cancelled)?;
        let range = temperature_span(query.lower, query.upper, query.temperature_axis)?;
        let reference = resolve_material_state_point(
            card,
            query.lower,
            &norton_requirements()?,
            query.selection,
        )?;
        for property in reference.properties() {
            poll(&mut cancelled)?;
            if !property.has_only_fixed_external_conditions(Some(query.temperature_axis), false) {
                return Err(CreepSourceError::Invalid(
                    "creep coefficients cannot freeze evolving stress, strain, rate or frequency axes",
                ));
            }
            let selected = property.answer().receipt.selected;
            let claim = card
                .claims()
                .claim(selected)
                .ok_or(CreepSourceError::Invalid(
                    "selected creep coefficient is absent from the original card",
                ))?;
            let line = card
                .claims()
                .query_line_typed(
                    &claim.key,
                    query.temperature_axis,
                    query.lower,
                    query.upper,
                    ClaimSelection::Pinned(selected),
                )
                .map_err(|e| CreepSourceError::Support {
                    property: property.requirement().name().into(),
                    reason: e.to_string(),
                })?;
            for (index, (_, answer)) in line.knots.iter().enumerate() {
                if index % 256 == 0 {
                    poll(&mut cancelled)?;
                }
                if answer.evidence.value.value != property.value_si() {
                    return Err(CreepSourceError::Invalid(
                        "Norton reference coefficients and elastic storage must be constant over the complete thermal span",
                    ));
                }
            }
        }
        let value = |name| {
            reference
                .property(name)
                .expect("required source coefficient")
                .value_si()
        };
        let elastic = IsotropicElastic::new(
            value(YOUNG_MODULUS_PROPERTY),
            value(POISSON_RATIO_PROPERTY),
            query.strain_limit,
        )
        .map_err(|_| CreepSourceError::Invalid("invalid sourced isotropic elastic coefficients"))?;
        let parameters = NortonCreepParameters {
            reference_stress_pa: value(CREEP_REFERENCE_STRESS_PROPERTY),
            reference_rate_s: value(CREEP_REFERENCE_RATE_PROPERTY),
            stress_exponent: value(CREEP_STRESS_EXPONENT_PROPERTY),
            activation_energy_j_mol: value(CREEP_ACTIVATION_ENERGY_PROPERTY),
            reference_temperature_k: value(CREEP_REFERENCE_TEMPERATURE_PROPERTY),
            temperature_range_k: range,
        };
        // Unit duration is admission only; every executed step supplies its own dt.
        NortonCreep::new(elastic.clone(), parameters.clone(), range[0], 1.0)?;
        NortonCreep::new(elastic.clone(), parameters.clone(), range[1], 1.0)?;
        Ok(Self {
            card,
            reference,
            elastic,
            parameters,
            temperature_axis: query.temperature_axis.into(),
        })
    }

    /// Complete original source bundle at the lower admitted temperature.
    pub fn reference(&self) -> &ResolvedMaterialStatePoint {
        &self.reference
    }
    pub const fn parameters(&self) -> &NortonCreepParameters {
        &self.parameters
    }

    /// Absolute-temperature coordinate retained in every original source receipt.
    pub fn temperature_axis(&self) -> &str {
        &self.temperature_axis
    }

    /// Bind the exact physical interval and endpoint temperature. This refreshes
    /// receipts with pinned claims and changes kinetics only; it owns no history
    /// and cannot advance a solver clock or mutate an accepted material state.
    pub fn at_temperature(
        &self,
        temperature_k: f64,
        duration_s: f64,
    ) -> Result<ResolvedNortonCreep, CreepSourceError> {
        let law = NortonCreep::new(
            self.elastic.clone(),
            self.parameters.clone(),
            temperature_k,
            duration_s,
        )?;
        let point = self
            .reference
            .query()
            .with(&self.temperature_axis, temperature_k)
            .map_err(|e| CreepSourceError::Support {
                property: self.temperature_axis.clone(),
                reason: e.to_string(),
            })?;
        let material = self.reference.at_query_point(self.card, &point)?;
        Ok(ResolvedNortonCreep { material, law })
    }
}

/// Native creep law together with every original selected source receipt.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedNortonCreep {
    material: ResolvedMaterialStatePoint,
    law: NortonCreep,
}
impl ResolvedNortonCreep {
    pub const fn law(&self) -> &NortonCreep {
        &self.law
    }
    pub fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }
}
impl SmallStrainLaw for ResolvedNortonCreep {
    type State = CreepState;
    fn initial_state(&self) -> CreepState {
        self.law.initial_state()
    }
    fn stress(&self, strain: &Voigt, state: &CreepState) -> Voigt {
        self.law.stress(strain, state)
    }
    fn tangent(&self, strain: &Voigt, state: &CreepState) -> Tangent6 {
        self.law.tangent(strain, state)
    }
    fn update_state(&self, strain: &Voigt, state: &CreepState) -> CreepState {
        self.law.update_state(strain, state)
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        self.law.admissibility()
    }
    fn card(&self) -> ModelCard {
        self.law.card()
    }
}

fn temperature_span(
    lower: &QueryPoint,
    upper: &QueryPoint,
    axis: &str,
) -> Result<[f64; 2], CreepSourceError> {
    let low = lower
        .axes()
        .get(axis)
        .copied()
        .ok_or(CreepSourceError::Invalid("missing lower temperature"))?;
    let high = upper
        .axes()
        .get(axis)
        .copied()
        .ok_or(CreepSourceError::Invalid("missing upper temperature"))?;
    let absolute = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    if !low.is_finite()
        || !high.is_finite()
        || low <= 0.0
        || low >= high
        || lower.axis_quantities() != upper.axis_quantities()
        || lower.axes().len() != upper.axes().len()
        || lower.axes().iter().any(|(name, value)| {
            upper
                .axes()
                .get(name)
                .is_none_or(|other| name != axis && value.to_bits() != other.to_bits())
        })
        || lower
            .axis_quantities()
            .get(axis)
            .is_some_and(|q| *q != absolute && *q != QuantitySpec::dimensional(Temperature::DIMS))
    {
        return Err(CreepSourceError::Invalid(
            "one increasing positive absolute-temperature axis with its original schema and unchanged source context is required",
        ));
    }
    Ok([low, high])
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), CreepSourceError> {
    if cancelled() {
        Err(CreepSourceError::Cancelled)
    } else {
        Ok(())
    }
}
