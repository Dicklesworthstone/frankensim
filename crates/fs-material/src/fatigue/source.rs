//! Original material-card S-N observations with explicit life/stress conventions.
//!
//! `resolve_log_log_sn` authors the same bounded log-log model as `SnCurve`.
//! It does not relabel that model as the source's interpolation policy: every
//! original knot, uncertainty, validity domain, selection and query receipt is
//! retained. No tensile-strength fit, endurance limit, runout interpretation,
//! temperature law or mean-stress transformation is inferred.

use std::sync::Arc;

use fs_blake3::ContentHash;
use fs_matdb::{
    ClaimSelection, MatDbError, MaterialAnswer, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, QueryPoint,
};
use fs_qty::semantic::{QuantityKind, ValueForm};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};

use super::{FatigueError, FatigueScope, SnCurve, SnPoint};

/// The original ordinate is identified by a distinct property key. A generic
/// strength, modulus or maximum stress cannot silently become an S-N amplitude.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SnStressMeasure {
    Amplitude,
    Range,
}
impl SnStressMeasure {
    pub fn property_key(self) -> PropertyKey {
        PropertyKey::new(
            match self {
                Self::Amplitude => "fatigue-stress-amplitude",
                Self::Range => "fatigue-stress-range",
            },
            Pressure::DIMS,
        )
    }
}

/// Two reversals make one cycle. The source axis names keep the conventions
/// distinct despite their identical dimension vectors.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SnLifeMeasure {
    Cycles,
    Reversals,
}
impl SnLifeMeasure {
    pub const fn axis_name(self) -> &'static str {
        match self {
            Self::Cycles => "cycles-to-failure",
            Self::Reversals => "reversals-to-failure",
        }
    }
}

/// The source must explicitly admit the query temperature and fully reversed
/// stress ratio R = sigma_min/sigma_max = -1. Other context coordinates travel
/// in the original typed query. Labels describe applicability; they do not
/// independently establish specimen or component equivalence.
#[derive(Debug, Clone, PartialEq)]
pub struct SnSourceRequest {
    pub stress_measure: SnStressMeasure,
    pub life_measure: SnLifeMeasure,
    pub temperature_axis: String,
    pub stress_ratio_axis: String,
    pub surface_and_geometry: String,
    pub environment: String,
    pub load_basis: String,
    pub maximum_points: usize,
}

/// Retained inside the immutable `SnCurve`, including when it is shared by a
/// live fatigue counter or a structural gauge. Source authority remains exactly
/// that of the original material claim; interpolation is an authored model.
#[derive(Debug, Clone, PartialEq)]
pub struct SnCalibrationSource {
    card_identity: ContentHash,
    material_state: MaterialStateId,
    claim: PropertyClaim,
    request: SnSourceRequest,
    query: QueryPoint,
    selection: ClaimSelection,
    /// Original life coordinate and its original stress answer, in source order.
    observations: Vec<(f64, MaterialAnswer)>,
}
impl SnCalibrationSource {
    pub const fn card_identity(&self) -> ContentHash {
        self.card_identity
    }
    pub fn material_state(&self) -> &MaterialStateId {
        &self.material_state
    }
    pub fn claim(&self) -> &PropertyClaim {
        &self.claim
    }
    pub fn request(&self) -> &SnSourceRequest {
        &self.request
    }
    pub fn query(&self) -> &QueryPoint {
        &self.query
    }
    pub const fn selection(&self) -> ClaimSelection {
        self.selection
    }
    pub fn observations(&self) -> &[(f64, MaterialAnswer)] {
        &self.observations
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum SnSourceError {
    Source(MatDbError),
    Fatigue(FatigueError),
    Invalid(&'static str),
}
impl core::fmt::Display for SnSourceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source-bound S-N calibration refused: {self:?}")
    }
}
impl std::error::Error for SnSourceError {}
impl From<MatDbError> for SnSourceError {
    fn from(error: MatDbError) -> Self {
        Self::Source(error)
    }
}
impl From<FatigueError> for SnSourceError {
    fn from(error: FatigueError) -> Self {
        Self::Fatigue(error)
    }
}

fn require(condition: bool, message: &'static str) -> Result<(), SnSourceError> {
    if condition {
        Ok(())
    } else {
        Err(SnSourceError::Invalid(message))
    }
}

fn query(
    card: &MaterialCard,
    key: &PropertyKey,
    point: &QueryPoint,
    selection: ClaimSelection,
) -> Result<MaterialAnswer, SnSourceError> {
    let answer = match selection {
        ClaimSelection::Policy(policy) => card.claims().query_typed(key, point, policy)?,
        ClaimSelection::Pinned(pin) => card.claims().query_pinned_typed(key, point, pin)?,
    };
    card.claims().verify_receipt(&answer.receipt)?;
    Ok(answer)
}

/// Build a finite-life log-log model from every original S-N knot. This explicit
/// constructor is the model-authoring choice, including for tabulated-only data.
/// Exact source knots and their receipts are never changed or synthesized.
///
/// The supplied query must select the source at an actual knot. All knots must
/// remain admitted under that same selection and frozen external context. Only
/// fully reversed normal-stress calibrations are admitted; R != -1 requires a
/// separately justified fatigue model, not relabeling or a second Goodman fit.
/// Source observations are lives TO FAILURE, not censored runouts. No endurance
/// plateau is inferred from the final measured point.
pub fn resolve_log_log_sn(
    card: &MaterialCard,
    point: &QueryPoint,
    request: &SnSourceRequest,
    selection: ClaimSelection,
) -> Result<SnCurve, SnSourceError> {
    let life_axis = request.life_measure.axis_name();
    require(
        request.maximum_points >= 2
            && !request.temperature_axis.trim().is_empty()
            && !request.stress_ratio_axis.trim().is_empty()
            && request.temperature_axis != request.stress_ratio_axis
            && request.temperature_axis != life_axis
            && request.stress_ratio_axis != life_axis,
        "declare distinct temperature, stress-ratio and life coordinates and at least two knots",
    )?;
    let key = request.stress_measure.property_key();
    let selected = query(card, &key, point, selection)?;
    let claim = card
        .claims()
        .claim(selected.receipt.selected)
        .ok_or(SnSourceError::Invalid("selected S-N source is absent"))?;
    let knots = match &claim.value {
        PropertyValue::Curve {
            abscissa,
            abscissa_dims,
            knots,
            dims,
        } if abscissa == life_axis && *abscissa_dims == Dims::NONE && *dims == Pressure::DIMS => {
            knots
        }
        _ => {
            return Err(SnSourceError::Invalid(
                "an original stress-versus-life table is required",
            ));
        }
    };
    require(
        knots.len() <= request.maximum_points,
        "original S-N table exceeds the declared point budget",
    )?;
    let dimensionless = QuantitySpec::dimensional(Dims::NONE);
    require(
        point.axis_quantities().get(life_axis) == Some(&dimensionless)
            && claim.validity.axis_quantities().get(life_axis) == Some(&dimensionless),
        "life must retain its explicitly typed cycles/reversals source coordinate",
    )?;
    require(
        point.axis_quantities().get(&request.stress_ratio_axis) == Some(&dimensionless)
            && claim
                .validity
                .axis_quantities()
                .get(&request.stress_ratio_axis)
                == Some(&dimensionless)
            && point.axes().get(&request.stress_ratio_axis) == Some(&-1.0)
            && claim.validity.bounds().get(&request.stress_ratio_axis) == Some(&(-1.0, -1.0)),
        "fully reversed calibration needs a source-declared R = -1, not an assumed mean stress",
    )?;
    let temperature_quantity = point
        .axis_quantities()
        .get(&request.temperature_axis)
        .ok_or(SnSourceError::Invalid(
            "temperature must have an explicit quantity descriptor",
        ))?;
    require(
        (*temperature_quantity == QuantitySpec::dimensional(Temperature::DIMS)
            || temperature_quantity.semantic_type().is_some_and(|q| {
                q.kind() == QuantityKind::AbsoluteTemperature && q.form() == ValueForm::Static
            }))
            && claim
                .validity
                .axis_quantities()
                .get(&request.temperature_axis)
                == Some(temperature_quantity)
            && claim
                .validity
                .bounds()
                .contains_key(&request.temperature_axis),
        "source temperature must be absolute kelvin and admitted by the selected claim",
    )?;
    let temperature_k = *point
        .axes()
        .get(&request.temperature_axis)
        .ok_or(SnSourceError::Invalid("source temperature is missing"))?;
    require(
        point
            .axes()
            .get(life_axis)
            .is_some_and(|n| knots.iter().any(|p| p.0 == *n)),
        "select an original life knot; the adapter does not mint an interpolated source anchor",
    )?;

    let mut observations = Vec::with_capacity(knots.len());
    let mut points = Vec::with_capacity(knots.len());
    for &(life, stress) in knots {
        let at = point.clone().with(life_axis, life)?;
        let answer = query(card, &key, &at, selection)?;
        require(
            answer.receipt.selected == selected.receipt.selected
                && answer.evidence.value.value == stress,
            "all original knots must retain the same selected source and exact stress",
        )?;
        let amplitude_pa = match request.stress_measure {
            SnStressMeasure::Amplitude => stress,
            SnStressMeasure::Range => 0.5 * stress,
        };
        let cycles_to_failure = match request.life_measure {
            SnLifeMeasure::Cycles => life,
            SnLifeMeasure::Reversals => 0.5 * life,
        };
        points.push(SnPoint {
            amplitude_pa,
            cycles_to_failure,
        });
        observations.push((life, answer));
    }
    // Original life increases while stress decreases. The shared S-N owner uses
    // increasing amplitude. Reversal preserves every observation; no sorting,
    // fitting, averaging, filtering of runouts, or monotonicity repair occurs.
    points.reverse();
    let id = card.id();
    let scope = FatigueScope {
        source_reference: claim.provenance.source.clone(),
        material: format!("{}; {}", id.chemistry, id.phase),
        process: id.process.clone(),
        surface_and_geometry: request.surface_and_geometry.clone(),
        environment: request.environment.clone(),
        load_basis: request.load_basis.clone(),
        temperature_k,
    };
    let mut curve = SnCurve::new(scope, points, None)?;
    curve.source = Some(Arc::new(SnCalibrationSource {
        card_identity: card.content_hash(),
        material_state: id.clone(),
        claim: claim.clone(),
        request: request.clone(),
        query: point.clone(),
        selection,
        observations,
    }));
    Ok(curve)
}
