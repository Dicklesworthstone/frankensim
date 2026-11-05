//! Original fracture-equivalent plastic strain versus stress-triaxiality data.
//!
//! The exact ordinate is accumulated equivalent plastic strain at ductile
//! initiation. Engineering elongation, total strain at rupture, reduction of
//! area and cyclic fatigue ductility are different quantities and never enter
//! through an alias. Triaxiality is mean tensile stress divided by von Mises
//! stress; compression therefore has negative triaxiality.
//!
//! Construction explicitly authors a piecewise-linear locus between original
//! observations. It does not change the claim's interpolation policy or create
//! measured values between knots. Every source knot and receipt is retained.
//! This one-dimensional reduction is frozen at the selected absolute
//! temperature, and uses only the original positive plastic-strain-rate
//! applicability interval. Both are checked against actual runtime inputs.
//! Lode-angle dependence and other unsupported evolving coordinates refuse.
//! Initiation is an empirical criterion, not stiffness degradation, fracture
//! energy, a separated surface, or experimentally validated component failure.

use fs_blake3::ContentHash;
use fs_matdb::{
    ClaimSelection, MatDbError, MaterialAnswer, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, QueryPoint,
};
use fs_qty::semantic::{QuantityKind, ValueForm};
use fs_qty::{Dims, QuantitySpec};

use super::{DuctileDomain, DuctileError, DuctileInitiationLaw, TabulatedDuctileLocus};
use crate::state_point::{ResolvedMaterialStatePoint, is_fixed_external_coordinate};

/// Local accumulated equivalent plastic strain at ductile initiation.
pub const FRACTURE_PLASTIC_STRAIN_PROPERTY: &str = "ductile-fracture-equivalent-plastic-strain";
/// `eta = mean tensile stress / von Mises stress`, dimensionless.
pub const STRESS_TRIAXIALITY_AXIS: &str = "stress-triaxiality";
/// Actual accumulated equivalent plastic strain rate, in inverse seconds.
/// A crosshead speed or cyclic frequency cannot supply this coordinate.
pub const PLASTIC_STRAIN_RATE_AXIS: &str = "equivalent-plastic-strain-rate";

const RATE_DIMS: Dims = Dims([0, 0, -1, 0, 0, 0]);

/// Fixed calibration temperature and bounded source preparation work.
/// Material, process and specimen provenance remain in the original card.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DuctileSourceRequest {
    /// Source axis carrying semantic static absolute temperature, in kelvin.
    pub temperature_axis: String,
    /// Maximum original knot count; at least two are needed for this locus.
    pub maximum_points: usize,
}

/// Original material identity, condition, claim, uncertainty and query receipts.
#[derive(Debug, Clone, PartialEq)]
pub struct DuctileCalibrationSource {
    card_identity: ContentHash,
    material_state: MaterialStateId,
    claim: PropertyClaim,
    request: DuctileSourceRequest,
    query: QueryPoint,
    selection: ClaimSelection,
    observations: Vec<(f64, MaterialAnswer)>,
}

impl DuctileCalibrationSource {
    pub const fn card_identity(&self) -> ContentHash {
        self.card_identity
    }
    pub fn material_state(&self) -> &MaterialStateId {
        &self.material_state
    }
    pub fn claim(&self) -> &PropertyClaim {
        &self.claim
    }
    pub fn request(&self) -> &DuctileSourceRequest {
        &self.request
    }
    pub fn query(&self) -> &QueryPoint {
        &self.query
    }
    pub const fn selection(&self) -> ClaimSelection {
        self.selection
    }
    /// Original `(triaxiality, answer)` observations, in source order.
    pub fn observations(&self) -> &[(f64, MaterialAnswer)] {
        &self.observations
    }

    /// Match the native mechanical material's original card, named state and
    /// complete fixed context. Only triaxiality and plastic strain rate are
    /// excluded: their actual runtime values are checked by the numeric law.
    /// Temperature remains the original fixed calibration temperature.
    pub fn matches_material_state(
        &self,
        card_identity: ContentHash,
        material: &MaterialStateId,
        point: &QueryPoint,
    ) -> bool {
        self.card_identity == card_identity
            && &self.material_state == material
            && self.query.axes().iter().all(|(axis, value)| {
                axis == STRESS_TRIAXIALITY_AXIS
                    || axis == PLASTIC_STRAIN_RATE_AXIS
                    || point.axes().get(axis) == Some(value)
                        && point.axis_quantities().get(axis)
                            == self.query.axis_quantities().get(axis)
            })
            && point.axes().iter().all(|(axis, value)| {
                axis == STRESS_TRIAXIALITY_AXIS
                    || axis == PLASTIC_STRAIN_RATE_AXIS
                    || self.query.axes().get(axis) == Some(value)
                        && self.query.axis_quantities().get(axis)
                            == point.axis_quantities().get(axis)
            })
    }

    /// Direct binding for the shared bundle used by source-resolved J2 laws.
    pub fn matches_material(&self, material: &ResolvedMaterialStatePoint) -> bool {
        self.matches_material_state(
            material.card_identity(),
            material.material(),
            &material.query(),
        )
    }
}

/// Prepared numeric locus together with immutable original-source evidence.
/// Evaluation neither queries a database nor allocates source objects.
#[derive(Debug, Clone, PartialEq)]
pub struct ResolvedDuctileLocus {
    law: TabulatedDuctileLocus,
    source: DuctileCalibrationSource,
}

impl ResolvedDuctileLocus {
    pub const fn law(&self) -> &TabulatedDuctileLocus {
        &self.law
    }
    pub const fn source(&self) -> &DuctileCalibrationSource {
        &self.source
    }
}

impl DuctileInitiationLaw for ResolvedDuctileLocus {
    fn domain(&self) -> &DuctileDomain {
        self.law.domain()
    }

    fn fracture_strain(
        &self,
        triaxiality: f64,
        strain_rate_s_inv: f64,
        temperature_k: f64,
    ) -> Result<f64, DuctileError> {
        self.law
            .fracture_strain(triaxiality, strain_rate_s_inv, temperature_k)
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum DuctileSourceError {
    Source(MatDbError),
    Model(DuctileError),
    Invalid(&'static str),
}

impl core::fmt::Display for DuctileSourceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source-bound ductile initiation refused: {self:?}")
    }
}
impl std::error::Error for DuctileSourceError {}
impl From<MatDbError> for DuctileSourceError {
    fn from(error: MatDbError) -> Self {
        Self::Source(error)
    }
}
impl From<DuctileError> for DuctileSourceError {
    fn from(error: DuctileError) -> Self {
        Self::Model(error)
    }
}

fn require(accepted: bool, reason: &'static str) -> Result<(), DuctileSourceError> {
    if accepted {
        Ok(())
    } else {
        Err(DuctileSourceError::Invalid(reason))
    }
}

fn query(
    card: &MaterialCard,
    key: &PropertyKey,
    point: &QueryPoint,
    selection: ClaimSelection,
) -> Result<MaterialAnswer, DuctileSourceError> {
    let answer = match selection {
        ClaimSelection::Policy(policy) => card.claims().query_typed(key, point, policy)?,
        ClaimSelection::Pinned(pin) => card.claims().query_pinned_typed(key, point, pin)?,
    };
    card.claims().verify_receipt(&answer.receipt)?;
    Ok(answer)
}

/// Author a bounded piecewise-linear initiation locus from every original knot.
/// The supplied point must select an actual source knot, a semantic absolute
/// temperature and a source-declared positive equivalent-plastic-strain rate.
/// No missing rate bound, temperature, fracture strain or source is inferred.
///
/// The original rate interval is an applicability statement, not a fitted rate
/// response. The temperature is held at the original query value even when the
/// claim admits a wider interval. Different temperatures require separately
/// resolved loci; history transfer and thermal coupling remain caller-owned.
pub fn resolve_ductile_locus(
    card: &MaterialCard,
    point: &QueryPoint,
    request: &DuctileSourceRequest,
    selection: ClaimSelection,
) -> Result<ResolvedDuctileLocus, DuctileSourceError> {
    require(
        request.maximum_points >= 2
            && !request.temperature_axis.trim().is_empty()
            && request.temperature_axis != STRESS_TRIAXIALITY_AXIS
            && request.temperature_axis != PLASTIC_STRAIN_RATE_AXIS,
        "declare a distinct temperature axis and at least two source knots",
    )?;
    let key = PropertyKey::new(FRACTURE_PLASTIC_STRAIN_PROPERTY, Dims::NONE);
    let selected = query(card, &key, point, selection)?;
    let claim = card
        .claims()
        .claim(selected.receipt.selected)
        .ok_or(DuctileSourceError::Invalid("selected ductile source is absent"))?;
    let knots = match &claim.value {
        PropertyValue::Curve {
            abscissa,
            abscissa_dims,
            knots,
            dims,
        } if abscissa == STRESS_TRIAXIALITY_AXIS
            && *abscissa_dims == Dims::NONE
            && *dims == Dims::NONE =>
        {
            knots
        }
        _ => {
            return Err(DuctileSourceError::Invalid(
                "original equivalent-plastic-strain versus tensile-positive triaxiality table required",
            ));
        }
    };
    require(
        knots.len() >= 2 && knots.len() <= request.maximum_points,
        "original ductile table exceeds the declared knot budget or has fewer than two knots",
    )?;
    for (axis, quantity) in [
        (STRESS_TRIAXIALITY_AXIS, QuantitySpec::dimensional(Dims::NONE)),
        (PLASTIC_STRAIN_RATE_AXIS, QuantitySpec::dimensional(RATE_DIMS)),
    ] {
        require(
            point.axis_quantities().get(axis) == Some(&quantity)
                && claim.validity.axis_quantities().get(axis) == Some(&quantity)
                && claim.validity.bounds().contains_key(axis),
            "source triaxiality and equivalent-plastic-strain rate need their exact typed bounds",
        )?;
    }
    let temperature = point
        .axis_quantities()
        .get(&request.temperature_axis)
        .ok_or(DuctileSourceError::Invalid("source temperature is untyped"))?;
    require(
        temperature.semantic_type().is_some_and(|q| {
            q.kind() == QuantityKind::AbsoluteTemperature && q.form() == ValueForm::Static
        }) && claim.validity.axis_quantities().get(&request.temperature_axis)
            == Some(temperature)
            && claim.validity.bounds().contains_key(&request.temperature_axis),
        "source temperature must be declared static absolute kelvin",
    )?;
    let temperature_k = *point
        .axes()
        .get(&request.temperature_axis)
        .ok_or(DuctileSourceError::Invalid("source temperature is absent"))?;
    require(
        point
            .axes()
            .get(STRESS_TRIAXIALITY_AXIS)
            .is_some_and(|eta| knots.iter().any(|knot| knot.0 == *eta)),
        "select an original triaxiality knot, not an interpolated source anchor",
    )?;
    for axis in claim.validity.bounds().keys() {
        if axis == STRESS_TRIAXIALITY_AXIS
            || axis == PLASTIC_STRAIN_RATE_AXIS
            || axis == &request.temperature_axis
        {
            continue;
        }
        require(
            is_fixed_external_coordinate(
                axis,
                point.axis_quantities().get(axis).copied(),
                point.axes().get(axis).copied(),
            ) && point.axis_quantities().get(axis).is_none_or(|quantity| {
                quantity.semantic_type().is_none_or(|q| {
                    q.kind() != QuantityKind::AbsoluteTemperature
                })
            }) && axis != "T",
            "a one-dimensional ductile locus cannot freeze an evolving or second temperature axis",
        )?;
    }
    let &(minimum_rate, maximum_rate) = claim
        .validity
        .bounds()
        .get(PLASTIC_STRAIN_RATE_AXIS)
        .ok_or(DuctileSourceError::Invalid("plastic strain rate bounds are absent"))?;
    let domain = DuctileDomain {
        triaxiality: [knots[0].0, knots[knots.len() - 1].0],
        strain_rate_s_inv: [minimum_rate, maximum_rate],
        temperature_k: [temperature_k, temperature_k],
    };
    let law = TabulatedDuctileLocus::new(knots.clone(), domain)?;
    let mut observations = Vec::with_capacity(knots.len());
    for &(eta, strain) in knots {
        let at = point.clone().with(STRESS_TRIAXIALITY_AXIS, eta)?;
        let answer = query(card, &key, &at, selection)?;
        require(
            answer.receipt.selected == selected.receipt.selected
                && answer.evidence.value.value == strain,
            "every original knot must retain its selected source and exact plastic strain",
        )?;
        observations.push((eta, answer));
    }
    Ok(ResolvedDuctileLocus {
        law,
        source: DuctileCalibrationSource {
            card_identity: card.content_hash(),
            material_state: card.id().clone(),
            claim: claim.clone(),
            request: request.clone(),
            query: point.clone(),
            selection,
            observations,
        },
    })
}
