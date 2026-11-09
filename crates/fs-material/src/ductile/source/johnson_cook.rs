//! A complete original Johnson--Cook initiation fit, without refitting data.
//!
//! The versioned model identity fixes `exp(-d3*eta)`, tensile-positive mean
//! stress, equivalent PLASTIC strain rate and the parent's normalized absolute
//! temperature. Original Johnson--Cook coefficients using `exp(+d3*eta)` must
//! be converted explicitly by their author before constructing this model card.
//! All eight scalar parameters belong to one pinned constitutive model; no
//! unrelated curves, coupon elongations or independently selected fits are
//! combined. The original model, calibration hashes and provenance are retained.
//!
//! Unlike an isothermal tabulated locus, this model's original typed temperature,
//! triaxiality and rate bounds remain live at every plastic increment. Other
//! external conditions stay fixed. Its rectangle is explicitly authored model
//! support, never a hull inferred from scattered observations. Runtime response
//! delegates to the same admitted numeric law and analytic derivatives; no
//! source query, allocation, new calibration or thermal-work model is introduced.

use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy};
use fs_qty::Temperature;

use super::*;
use crate::ductile::{JohnsonCookInitiation, JohnsonCookParameters};

/// Version-1 law with the explicit `exp(-d3*eta)` triaxiality convention.
pub const JOHNSON_COOK_INITIATION_LAW_ID: &str = "ductile-johnson-cook-minus-triaxiality";

const PARAMETERS: [(&str, Dims); 8] = [
    ("d1", Dims::NONE),
    ("d2", Dims::NONE),
    ("d3", Dims::NONE),
    ("d4", Dims::NONE),
    ("d5", Dims::NONE),
    ("reference_strain_rate_s_inv", RATE_DIMS),
    ("transition_temperature_k", Temperature::DIMS),
    ("melting_temperature_k", Temperature::DIMS),
];

/// Original material state, whole-fit identity, calibration and fixed context.
#[derive(Clone, Debug, PartialEq)]
pub struct JohnsonCookCalibrationSource {
    card_identity: ContentHash,
    material_state: MaterialStateId,
    model: ConstitutiveModelCard,
    query: QueryPoint,
    temperature_axis: String,
}

impl JohnsonCookCalibrationSource {
    pub const fn card_identity(&self) -> ContentHash {
        self.card_identity
    }

    pub fn material_state(&self) -> &MaterialStateId {
        &self.material_state
    }

    /// Complete immutable original fit, including all source hashes/provenance.
    pub const fn model_card(&self) -> &ConstitutiveModelCard {
        &self.model
    }

    pub fn model_identity(&self) -> ContentHash {
        self.model.content_hash()
    }

    /// Original complete typed observation used for admission, not frozen physics.
    pub const fn query(&self) -> &QueryPoint {
        &self.query
    }

    pub fn temperature_axis(&self) -> &str {
        &self.temperature_axis
    }

    fn is_runtime_axis(&self, axis: &str) -> bool {
        axis == STRESS_TRIAXIALITY_AXIS
            || axis == PLASTIC_STRAIN_RATE_AXIS
            || axis == self.temperature_axis
    }

    fn admits_runtime_coordinate(
        &self,
        axis: &str,
        value: f64,
        quantity: Option<&QuantitySpec>,
    ) -> bool {
        value.is_finite()
            && quantity == self.query.axis_quantities().get(axis)
            && self
                .model
                .validity
                .bounds()
                .get(axis)
                .is_some_and(|&(lo, hi)| value >= lo && value <= hi)
    }

    /// Bind the same material card/state and complete held environment.
    /// Runtime eta, rate and temperature may be absent from a mechanical bundle;
    /// any present coordinate must retain its exact type and original support.
    /// Their actual trial values remain checked by the numeric initiation law.
    pub fn matches_material_state(
        &self,
        card_identity: ContentHash,
        material: &MaterialStateId,
        point: &QueryPoint,
    ) -> bool {
        self.card_identity == card_identity
            && &self.material_state == material
            && self.query.axes().iter().all(|(axis, value)| {
                self.is_runtime_axis(axis)
                    || point.axes().get(axis) == Some(value)
                        && point.axis_quantities().get(axis)
                            == self.query.axis_quantities().get(axis)
            })
            && point.axes().iter().all(|(axis, value)| {
                if self.is_runtime_axis(axis) {
                    self.admits_runtime_coordinate(axis, *value, point.axis_quantities().get(axis))
                } else {
                    self.query.axes().get(axis) == Some(value)
                        && self.query.axis_quantities().get(axis)
                            == point.axis_quantities().get(axis)
                }
            })
    }

    pub fn matches_material(&self, material: &ResolvedMaterialStatePoint) -> bool {
        self.card_identity == material.card_identity()
            && &self.material_state == material.material()
            && self.query.axes().iter().all(|(axis, value)| {
                self.is_runtime_axis(axis)
                    || material.state_coordinate(axis) == Some(*value)
                        && material.axis_quantities().get(axis)
                            == self.query.axis_quantities().get(axis)
            })
            && material.query_point().iter().all(|(axis, value)| {
                if self.is_runtime_axis(axis) {
                    self.admits_runtime_coordinate(
                        axis,
                        *value,
                        material.axis_quantities().get(axis),
                    )
                } else {
                    self.query.axes().get(axis) == Some(value)
                        && self.query.axis_quantities().get(axis)
                            == material.axis_quantities().get(axis)
                }
            })
    }
}

/// Prepared rate/temperature-dependent initiation retaining one original fit.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedJohnsonCookInitiation {
    law: JohnsonCookInitiation,
    source: JohnsonCookCalibrationSource,
}

impl ResolvedJohnsonCookInitiation {
    pub const fn law(&self) -> &JohnsonCookInitiation {
        &self.law
    }

    pub const fn source(&self) -> &JohnsonCookCalibrationSource {
        &self.source
    }
}

impl DuctileInitiationLaw for ResolvedJohnsonCookInitiation {
    fn domain(&self) -> &DuctileDomain {
        self.law.domain()
    }

    fn fracture_strain(&self, eta: f64, rate: f64, temperature: f64) -> Result<f64, DuctileError> {
        self.law.fracture_strain(eta, rate, temperature)
    }
}

impl DifferentiableDuctileInitiationLaw for ResolvedJohnsonCookInitiation {
    fn fracture_strain_gradient(
        &self,
        eta: f64,
        rate: f64,
        temperature: f64,
    ) -> Result<DuctileLocusGradient, DuctileError> {
        self.law.fracture_strain_gradient(eta, rate, temperature)
    }
}

/// Resolve one explicitly pinned whole Johnson--Cook fit from an original card.
///
/// Law/state version 1 requires virgin initiation history (zero accumulated
/// equivalent plastic strain and utilization). `d1` through `d5` are finite
/// dimensionless coefficients; reference rate is positive 1/s; the two
/// temperature parameters are absolute kelvin with `0 < transition < melting`.
/// Extra or missing parameters refuse rather than choose an implicit convention.
///
/// Original validity must explicitly bound typed tensile-positive triaxiality,
/// positive equivalent plastic strain rate, and one named semantic absolute
/// temperature. All three runtime bounds must be finite. The model's entire
/// rectangle is checked by the numeric constructor, including factor positivity.
/// A supplied observation must lie inside it; it never expands that support.
pub fn resolve_johnson_cook_initiation(
    card: &MaterialCard,
    model: ContentHash,
    point: &QueryPoint,
    temperature_axis: &str,
) -> Result<ResolvedJohnsonCookInitiation, DuctileSourceError> {
    require(
        !temperature_axis.trim().is_empty()
            && ![
                STRESS_TRIAXIALITY_AXIS,
                PLASTIC_STRAIN_RATE_AXIS,
                LODE_PARAMETER_AXIS,
            ]
            .contains(&temperature_axis),
        "Johnson-Cook needs one distinct named absolute-temperature coordinate",
    )?;
    let original = card
        .models()
        .iter()
        .find(|candidate| candidate.content_hash() == model)
        .ok_or(DuctileSourceError::Invalid(
            "pinned Johnson-Cook model is absent from the original material card",
        ))?;
    original.validate()?;
    require(
        original.law.0 == JOHNSON_COOK_INITIATION_LAW_ID
            && original.law_version == 1
            && original.state_schema_version == 1
            && original.initial_state == InitialStatePolicy::ZeroInternalState,
        "Johnson-Cook initiation requires the explicit minus-triaxiality law and virgin law/state v1",
    )?;
    require(
        !original.sources.is_empty(),
        "source-bound Johnson-Cook requires original calibration artifact identities",
    )?;
    require(
        original.parameters.len() == PARAMETERS.len()
            && PARAMETERS.iter().all(|&(name, dims)| {
                original
                    .parameters
                    .get(name)
                    .is_some_and(|p| p.dims == dims && p.value.is_finite())
            }),
        "Johnson-Cook requires all eight exact scalar parameter names and SI dimensions",
    )?;
    let temperature = point
        .axis_quantities()
        .get(temperature_axis)
        .copied()
        .ok_or(DuctileSourceError::Invalid(
            "Johnson-Cook source temperature is untyped",
        ))?;
    require(
        temperature.semantic_type().is_some_and(|q| {
            q.kind() == QuantityKind::AbsoluteTemperature && q.form() == ValueForm::Static
        }),
        "Johnson-Cook source temperature must be static absolute kelvin",
    )?;
    for (axis, quantity) in [
        (
            STRESS_TRIAXIALITY_AXIS,
            QuantitySpec::dimensional(Dims::NONE),
        ),
        (
            PLASTIC_STRAIN_RATE_AXIS,
            QuantitySpec::dimensional(RATE_DIMS),
        ),
        (temperature_axis, temperature),
    ] {
        require(
            original.validity.axis_quantities().get(axis) == Some(&quantity)
                && point.axis_quantities().get(axis) == Some(&quantity)
                && original.validity.bounds().contains_key(axis),
            "Johnson-Cook needs original typed triaxiality, plastic-rate and absolute-temperature bounds",
        )?;
    }
    require(
        original
            .validity
            .contains_typed(point.axes(), point.axis_quantities()),
        "Johnson-Cook observation is missing, mistyped or outside original calibration support",
    )?;
    for (axis, &value) in point.axes() {
        if [
            STRESS_TRIAXIALITY_AXIS,
            PLASTIC_STRAIN_RATE_AXIS,
            temperature_axis,
        ]
        .contains(&axis.as_str())
        {
            continue;
        }
        let quantity = point.axis_quantities().get(axis).copied();
        require(
            is_fixed_external_coordinate(axis, quantity, Some(value))
                && axis != "T"
                && quantity.is_none_or(|q| {
                    q.semantic_type()
                        .is_none_or(|t| t.kind() != QuantityKind::AbsoluteTemperature)
                }),
            "Johnson-Cook cannot freeze another evolving coordinate or second temperature axis",
        )?;
    }
    let bounds = |axis: &str| {
        let &(lo, hi) = &original.validity.bounds()[axis];
        [lo, hi]
    };
    let parameter = |name: &str| original.parameters[name].value;
    let law = JohnsonCookInitiation::new(
        JohnsonCookParameters {
            d1: parameter("d1"),
            d2: parameter("d2"),
            d3: parameter("d3"),
            d4: parameter("d4"),
            d5: parameter("d5"),
            reference_strain_rate_s_inv: parameter("reference_strain_rate_s_inv"),
            transition_temperature_k: parameter("transition_temperature_k"),
            melting_temperature_k: parameter("melting_temperature_k"),
        },
        DuctileDomain {
            triaxiality: bounds(STRESS_TRIAXIALITY_AXIS),
            strain_rate_s_inv: bounds(PLASTIC_STRAIN_RATE_AXIS),
            temperature_k: bounds(temperature_axis),
        },
    )?;
    Ok(ResolvedJohnsonCookInitiation {
        law,
        source: JohnsonCookCalibrationSource {
            card_identity: card.content_hash(),
            material_state: card.id().clone(),
            model: original.clone(),
            query: point.clone(),
            temperature_axis: temperature_axis.into(),
        },
    })
}
