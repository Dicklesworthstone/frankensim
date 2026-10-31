//! Original torsional strain-life calibration with explicit isotropic elasticity.
//!
//! Four measured/fitted shear coefficients retain one source fit and specimen
//! population. The separately calibrated normal-stress sensitivity retains its
//! own receipt and life support. No axial-to-shear approximation, inferred yield
//! normalization, or default sensitivity is introduced. Engineering shear and
//! the signed/tension-only normal-stress conventions have distinct source keys.
//!
//! This adapter explicitly chooses isotropic elasticity: G = E / (2 * (1 + nu)).
//! Original E, nu and INITIAL yield-stress answers remain attached so a native
//! J2 consumer can match the actual constitutive sources, not just equal numbers.
//! Source power terms fitted against N and 2N are converted without changing the
//! original claims. Predictions stay within the intersection of fitted support.

use std::sync::Arc;

use fs_blake3::ContentHash;
use fs_matdb::{
    ClaimSelection, MatDbError, MaterialAnswer, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, QueryPoint,
};
use fs_math::det;
use fs_qty::semantic::{QuantityKind, ValueForm};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};

use super::{
    FatemiSocieCalibration, FatemiSocieError, FatemiSocieNormalStress, FatemiSocieParameters,
};
use crate::fatigue::FatigueScope;

/// Explicit per-property selection order. Equal dimensions do not make axial
/// and shear fatigue data, engineering and tensor strain, or yield conventions
/// interchangeable.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FatemiSocieProperty {
    YoungsModulus,
    PoissonRatio,
    InitialYieldStress,
    ShearFatigueStrengthCoefficient,
    ShearFatigueStrengthExponent,
    ShearFatigueDuctilityCoefficient,
    ShearFatigueDuctilityExponent,
    NormalStressSensitivity,
}

impl FatemiSocieProperty {
    pub const ALL: [Self; 8] = [
        Self::YoungsModulus,
        Self::PoissonRatio,
        Self::InitialYieldStress,
        Self::ShearFatigueStrengthCoefficient,
        Self::ShearFatigueStrengthExponent,
        Self::ShearFatigueDuctilityCoefficient,
        Self::ShearFatigueDuctilityExponent,
        Self::NormalStressSensitivity,
    ];

    pub const fn index(self) -> usize {
        self as usize
    }

    /// The sensitivity source explicitly names the normal-stress convention
    /// used in its calibration. Switching convention cannot reuse that claim
    /// implicitly. Other property keys do not depend on the convention.
    pub fn property_key(self, normal_stress: FatemiSocieNormalStress) -> PropertyKey {
        let (name, dims) = match self {
            Self::YoungsModulus => ("young_modulus", Pressure::DIMS),
            Self::PoissonRatio => ("poisson_ratio", Dims::NONE),
            Self::InitialYieldStress => ("yield_stress", Pressure::DIMS),
            Self::ShearFatigueStrengthCoefficient => {
                ("shear_fatigue_strength_coefficient", Pressure::DIMS)
            }
            Self::ShearFatigueStrengthExponent => ("shear_fatigue_strength_exponent", Dims::NONE),
            Self::ShearFatigueDuctilityCoefficient => (
                "engineering_shear_fatigue_ductility_coefficient",
                Dims::NONE,
            ),
            Self::ShearFatigueDuctilityExponent => ("shear_fatigue_ductility_exponent", Dims::NONE),
            Self::NormalStressSensitivity => (
                match normal_stress {
                    FatemiSocieNormalStress::Signed => {
                        "fatemi_socie_signed_normal_stress_sensitivity"
                    }
                    FatemiSocieNormalStress::TensionOnly => {
                        "fatemi_socie_tension_only_normal_stress_sensitivity"
                    }
                },
                Dims::NONE,
            ),
        };
        PropertyKey::new(name, dims)
    }
}

/// Independent variable in BOTH original shear strain-life power terms.
/// The output life always counts full cycles; this is the original fit basis.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FatemiSociePowerBasis {
    FullCycles,
    Reversals,
}

impl FatemiSociePowerBasis {
    pub const fn life_axis(self) -> &'static str {
        match self {
            Self::FullCycles => "cycles-to-failure",
            Self::Reversals => "reversals-to-failure",
        }
    }
}

/// Explicit source conventions and applicability labels. The baseline ratio
/// is tau_min/tau_max for the original fully reversed torsional calibration.
/// Labels declare transfer scope; they do not prove component applicability.
#[derive(Debug, Clone, PartialEq)]
pub struct FatemiSocieSourceRequest {
    pub power_basis: FatemiSociePowerBasis,
    pub normal_stress: FatemiSocieNormalStress,
    pub temperature_axis: String,
    pub shear_stress_ratio_axis: String,
    pub surface_and_geometry: String,
    pub environment: String,
    pub load_basis: String,
}

/// Original source evidence shared by the immutable calibration and consumers.
#[derive(Debug, Clone, PartialEq)]
pub struct FatemiSocieCalibrationSource {
    card_identity: ContentHash,
    material_state: MaterialStateId,
    request: FatemiSocieSourceRequest,
    query: QueryPoint,
    selections: [ClaimSelection; 8],
    properties: Vec<(PropertyClaim, MaterialAnswer)>,
}

impl FatemiSocieCalibrationSource {
    pub const fn card_identity(&self) -> ContentHash {
        self.card_identity
    }
    pub fn material_state(&self) -> &MaterialStateId {
        &self.material_state
    }
    pub fn request(&self) -> &FatemiSocieSourceRequest {
        &self.request
    }
    pub fn query(&self) -> &QueryPoint {
        &self.query
    }
    pub fn selections(&self) -> &[ClaimSelection; 8] {
        &self.selections
    }
    pub fn property(&self, property: FatemiSocieProperty) -> &MaterialAnswer {
        &self.properties[property.index()].1
    }
    pub fn claim(&self, property: FatemiSocieProperty) -> &PropertyClaim {
        &self.properties[property.index()].0
    }

    /// Require the original physical sources and complete external state.
    /// Only calibration life and reversed torsion-ratio axes are excluded from
    /// the operating query. The initial yield source is not the evolving yield
    /// surface or a separately fitted cyclic-yield measurement.
    pub fn matches_material_state(
        &self,
        card_identity: ContentHash,
        material: &MaterialStateId,
        point: &QueryPoint,
        youngs: &MaterialAnswer,
        poisson: &MaterialAnswer,
        initial_yield: &MaterialAnswer,
    ) -> bool {
        self.card_identity == card_identity
            && &self.material_state == material
            && [
                (FatemiSocieProperty::YoungsModulus, youngs),
                (FatemiSocieProperty::PoissonRatio, poisson),
                (FatemiSocieProperty::InitialYieldStress, initial_yield),
            ]
            .iter()
            .all(|(property, native)| {
                let original = self.property(*property);
                original.receipt.selected == native.receipt.selected
                    && original.evidence.value == native.evidence.value
            })
            && self.query.axes().iter().all(|(axis, value)| {
                axis == &self.request.shear_stress_ratio_axis
                    || axis == self.request.power_basis.life_axis()
                    || point.axes().get(axis) == Some(value)
                        && point.axis_quantities().get(axis)
                            == self.query.axis_quantities().get(axis)
            })
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum FatemiSocieSourceError {
    Source(MatDbError),
    Model(FatemiSocieError),
    Invalid(&'static str),
}

impl core::fmt::Display for FatemiSocieSourceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source-bound Fatemi-Socie calibration refused: {self:?}")
    }
}
impl std::error::Error for FatemiSocieSourceError {}
impl From<MatDbError> for FatemiSocieSourceError {
    fn from(error: MatDbError) -> Self {
        Self::Source(error)
    }
}
impl From<FatemiSocieError> for FatemiSocieSourceError {
    fn from(error: FatemiSocieError) -> Self {
        Self::Model(error)
    }
}

fn require(condition: bool, reason: &'static str) -> Result<(), FatemiSocieSourceError> {
    if condition {
        Ok(())
    } else {
        Err(FatemiSocieSourceError::Invalid(reason))
    }
}

/// Resolve original engineering-shear strain-life data and explicitly isotropic
/// elasticity. Selections follow `FatemiSocieProperty::ALL` and may individually
/// pin original claims. All parameters support the actual absolute temperature;
/// fixed scalars and source-supported temperature curves are admitted.
///
/// The four torsional coefficients must share fit provenance, observations and
/// life bounds. The independently calibrated k may have different provenance
/// and support; its life interval intersects the torsional fit. Each fitted
/// parameter explicitly declares the reversed torsional baseline R = -1.
/// E, nu and initial yield remain physical-state properties independent of those
/// calibration axes. No life-dependent coefficient is frozen into a prediction.
pub fn resolve_isotropic_fatemi_socie(
    card: &MaterialCard,
    point: &QueryPoint,
    request: &FatemiSocieSourceRequest,
    selections: [ClaimSelection; 8],
) -> Result<FatemiSocieCalibration, FatemiSocieSourceError> {
    let life_axis = request.power_basis.life_axis();
    require(
        !request.temperature_axis.trim().is_empty()
            && !request.shear_stress_ratio_axis.trim().is_empty()
            && request.temperature_axis != request.shear_stress_ratio_axis
            && request.temperature_axis != life_axis
            && request.shear_stress_ratio_axis != life_axis,
        "declare distinct temperature, torsional baseline ratio and life axes",
    )?;
    let temperature = point
        .axis_quantities()
        .get(&request.temperature_axis)
        .ok_or(FatemiSocieSourceError::Invalid(
            "typed absolute temperature is missing",
        ))?;
    require(
        temperature.semantic_type().is_some_and(|q| {
            q.kind() == QuantityKind::AbsoluteTemperature && q.form() == ValueForm::Static
        }),
        "source temperature must be semantic static absolute temperature in kelvin",
    )?;
    let temperature_k =
        *point
            .axes()
            .get(&request.temperature_axis)
            .ok_or(FatemiSocieSourceError::Invalid(
                "source temperature value is missing",
            ))?;
    require(
        temperature_k.is_finite() && temperature_k > 0.0,
        "source temperature must be finite positive kelvin",
    )?;
    let dimensionless = QuantitySpec::dimensional(Dims::NONE);
    require(
        point.axis_quantities().get(life_axis) == Some(&dimensionless)
            && point
                .axis_quantities()
                .get(&request.shear_stress_ratio_axis)
                == Some(&dimensionless)
            && point.axes().get(&request.shear_stress_ratio_axis) == Some(&-1.0),
        "declare typed source life and fully reversed torsion R = -1",
    )?;

    let mut properties = Vec::with_capacity(8);
    for property in FatemiSocieProperty::ALL {
        let key = property.property_key(request.normal_stress);
        let answer = match selections[property.index()] {
            ClaimSelection::Policy(policy) => card.claims().query_typed(&key, point, policy)?,
            ClaimSelection::Pinned(pin) => card.claims().query_pinned_typed(&key, point, pin)?,
        };
        card.claims().verify_receipt(&answer.receipt)?;
        let claim =
            card.claims()
                .claim(answer.receipt.selected)
                .ok_or(FatemiSocieSourceError::Invalid(
                    "selected original source claim is absent",
                ))?;
        require(
            claim
                .validity
                .axis_quantities()
                .get(&request.temperature_axis)
                == Some(temperature)
                && claim
                    .validity
                    .bounds()
                    .contains_key(&request.temperature_axis),
            "every physical and fitted parameter needs an explicit matching source temperature domain",
        )?;
        require(
            match &claim.value {
                PropertyValue::Scalar { .. } => true,
                PropertyValue::Curve {
                    abscissa,
                    abscissa_dims,
                    ..
                } => abscissa == &request.temperature_axis && *abscissa_dims == Temperature::DIMS,
                _ => false,
            },
            "parameters must be fixed scalars or source functions of the declared temperature",
        )?;
        if property.index() >= FatemiSocieProperty::ShearFatigueStrengthCoefficient.index() {
            require(
                claim
                    .validity
                    .axis_quantities()
                    .get(&request.shear_stress_ratio_axis)
                    == Some(&dimensionless)
                    && claim
                        .validity
                        .bounds()
                        .get(&request.shear_stress_ratio_axis)
                        == Some(&(-1.0, -1.0))
                    && claim.validity.axis_quantities().get(life_axis) == Some(&dimensionless)
                    && claim.validity.bounds().contains_key(life_axis),
                "fitted parameters need explicit R = -1 and typed finite life support",
            )?;
        } else {
            require(
                !claim.validity.bounds().contains_key(life_axis)
                    && !claim
                        .validity
                        .bounds()
                        .contains_key(&request.shear_stress_ratio_axis),
                "E, nu and initial yield must be physical-state properties, not fatigue-life fits",
            )?;
        }
        properties.push((claim.clone(), answer));
    }

    let fit = &properties[FatemiSocieProperty::ShearFatigueStrengthCoefficient.index()].0;
    let fit_domain = fit.validity.bounds()[life_axis];
    for (claim, _) in &properties[4..7] {
        require(
            claim.provenance == fit.provenance
                && claim.observations == fit.observations
                && claim.validity.bounds()[life_axis] == fit_domain,
            "the four torsional coefficients must retain one source fit, observation population and life interval",
        )?;
    }
    let sensitivity_domain = properties[FatemiSocieProperty::NormalStressSensitivity.index()]
        .0
        .validity
        .bounds()[life_axis];
    let minimum_source_life = match request.power_basis {
        FatemiSociePowerBasis::FullCycles => 0.5,
        FatemiSociePowerBasis::Reversals => 1.0,
    };
    for original in [fit_domain, sensitivity_domain] {
        require(
            original.0.is_finite()
                && original.0 >= minimum_source_life
                && original.1.is_finite()
                && original.1 > original.0,
            "each original fit needs a finite increasing life interval starting at one reversal or later",
        )?;
    }
    let domain = (
        fit_domain.0.max(sensitivity_domain.0),
        fit_domain.1.min(sensitivity_domain.1),
    );
    require(
        domain.0 < domain.1,
        "torsional fit and normal-stress sensitivity need a common finite life interval",
    )?;
    let value = |property: FatemiSocieProperty| properties[property.index()].1.evidence.value.value;
    let youngs = value(FatemiSocieProperty::YoungsModulus);
    let poisson = value(FatemiSocieProperty::PoissonRatio);
    require(
        youngs.is_finite()
            && youngs > 0.0
            && poisson.is_finite()
            && poisson > -1.0
            && poisson < 0.5,
        "isotropic elasticity needs positive E and -1 < nu < 0.5",
    )?;
    let b = value(FatemiSocieProperty::ShearFatigueStrengthExponent);
    let c = value(FatemiSocieProperty::ShearFatigueDuctilityExponent);
    let convert = |coefficient: f64, exponent: f64| {
        if request.power_basis == FatemiSociePowerBasis::Reversals || coefficient == 0.0 {
            coefficient
        } else {
            coefficient * det::exp(-exponent * core::f64::consts::LN_2)
        }
    };
    let parameters = FatemiSocieParameters {
        shear_modulus_pa: youngs / (2.0 * (1.0 + poisson)),
        shear_fatigue_strength_coefficient_pa: convert(
            value(FatemiSocieProperty::ShearFatigueStrengthCoefficient),
            b,
        ),
        shear_fatigue_strength_exponent: b,
        shear_fatigue_ductility_coefficient: convert(
            value(FatemiSocieProperty::ShearFatigueDuctilityCoefficient),
            c,
        ),
        shear_fatigue_ductility_exponent: c,
        normal_stress_sensitivity: value(FatemiSocieProperty::NormalStressSensitivity),
        initial_yield_stress_pa: value(FatemiSocieProperty::InitialYieldStress),
        normal_stress: request.normal_stress,
    };
    let life_scale = if request.power_basis == FatemiSociePowerBasis::Reversals {
        0.5
    } else {
        1.0
    };
    let scope = FatigueScope {
        source_reference: fit.provenance.source.clone(),
        material: card.id().chemistry.clone(),
        process: card.id().process.clone(),
        surface_and_geometry: request.surface_and_geometry.clone(),
        environment: request.environment.clone(),
        load_basis: request.load_basis.clone(),
        temperature_k,
    };
    let mut calibration = FatemiSocieCalibration::new(
        scope,
        parameters,
        domain.0 * life_scale,
        domain.1 * life_scale,
    )?;
    calibration.source = Some(Arc::new(FatemiSocieCalibrationSource {
        card_identity: card.content_hash(),
        material_state: card.id().clone(),
        request: request.clone(),
        query: point.clone(),
        selections,
        properties,
    }));
    Ok(calibration)
}
