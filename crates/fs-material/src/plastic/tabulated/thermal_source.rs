//! Original thermal hardening observations on explicitly interpreted J2 owners.
//!
//! A complete Cartesian grid retains every measured accumulated-plastic-strain
//! knot and every original sample receipt. The explicit receiving model makes
//! the yield stress piecewise linear in strain and affine in absolute
//! temperature. A TabulatedOnly claim stays observational: its interpolation
//! policy is never rewritten. A LinearInside claim must additionally describe
//! the same globally affine plane through its own cells, with complete source
//! support. Triangulated and bilinear interpolation are not interchangeable.
//!
//! Preparation is bounded and outside constitutive iterations. No temperature
//! extrapolation, fitted thermal kink, engineering-strain conversion, stress
//! relabeling, or derivative from a single temperature is admitted.

use super::source::{ACCUMULATED_PLASTIC_STRAIN_AXIS, YieldStressMeasure};
use super::thermal::AffineTabulatedHardening;
use super::{HardeningPoint, TabulatedError, TabulatedHardening};
use crate::plastic::finite::thermal::{FiniteJ2ThermalError, FiniteJ2ThermalTabulatedHardening};
use crate::state_point::{affine_temperature_slope, is_fixed_external_coordinate};
use fs_blake3::ContentHash;
use fs_matdb::{
    ClaimSelection, EnvelopeAnswer, InterpolationPolicy, MatDbError, MaterialAnswer, MaterialCard,
    MaterialStateId, PropertyClaim, PropertyKey, PropertySupportError, PropertyValue, QueryPoint,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};

/// An explicit family interval and preparation budget; no defaults infer support.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermalHardeningCurveRequest {
    pub temperature_axis: String,
    pub reference_temperature_k: f64,
    pub temperature_range_k: [f64; 2],
    pub maximum_samples: usize,
    pub maximum_knots: usize,
}

/// Original source, selected claim, fixed context and every original observation.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermalHardeningCurveSource {
    card_identity: ContentHash,
    material_state: MaterialStateId,
    stress_measure: YieldStressMeasure,
    context: QueryPoint,
    request: ThermalHardeningCurveRequest,
    selection: ClaimSelection,
    selection_answer: MaterialAnswer,
    claim: PropertyClaim,
    observations: Vec<(QueryPoint, MaterialAnswer)>,
    continuous_support: Option<EnvelopeAnswer>,
}

impl ThermalHardeningCurveSource {
    pub const fn card_identity(&self) -> ContentHash {
        self.card_identity
    }
    pub const fn material_state(&self) -> &MaterialStateId {
        &self.material_state
    }
    /// Exact original stress interpretation; no current/reference-volume conversion.
    pub const fn stress_measure(&self) -> YieldStressMeasure {
        self.stress_measure
    }
    /// Original reference-temperature context, excluding accumulated strain.
    pub const fn context(&self) -> &QueryPoint {
        &self.context
    }
    pub const fn request(&self) -> &ThermalHardeningCurveRequest {
        &self.request
    }
    pub const fn selection(&self) -> ClaimSelection {
        self.selection
    }
    /// Original selection at zero plastic strain and the lower family temperature.
    pub const fn selection_answer(&self) -> &MaterialAnswer {
        &self.selection_answer
    }
    pub const fn claim(&self) -> &PropertyClaim {
        &self.claim
    }
    /// Every vertex, in original source order, including outside a narrower family.
    pub fn observations(&self) -> &[(QueryPoint, MaterialAnswer)] {
        &self.observations
    }
    /// Present only when original LinearInside cells cover the complete rectangle.
    pub const fn continuous_support(&self) -> Option<&EnvelopeAnswer> {
        self.continuous_support.as_ref()
    }
}

/// Affine ordinate evolution with its original source evidence and operating span.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedThermalTabulatedHardening {
    hardening: FiniteJ2ThermalTabulatedHardening,
    source: ThermalHardeningCurveSource,
}

impl ResolvedThermalTabulatedHardening {
    pub const fn hardening(&self) -> &FiniteJ2ThermalTabulatedHardening {
        &self.hardening
    }
    pub const fn source(&self) -> &ThermalHardeningCurveSource {
        &self.source
    }
}

/// Original Cauchy yield observations and their shared affine storage preparation.
/// This is a distinct source entry point from the original finite Kirchhoff law.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedThermalCauchyTabulatedHardening {
    hardening: AffineTabulatedHardening,
    source: ThermalHardeningCurveSource,
}

impl ResolvedThermalCauchyTabulatedHardening {
    pub const fn hardening(&self) -> &AffineTabulatedHardening {
        &self.hardening
    }
    pub const fn source(&self) -> &ThermalHardeningCurveSource {
        &self.source
    }
}

/// Selection, original support, native curve, thermal arithmetic or admission refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum ThermalHardeningSourceError {
    Source(MatDbError),
    Support(PropertySupportError),
    Curve(TabulatedError),
    Native(FiniteJ2ThermalError),
    Invalid(&'static str),
}

impl core::fmt::Display for ThermalHardeningSourceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Source(e) => write!(f, "thermal hardening source: {e}"),
            Self::Support(e) => write!(f, "thermal hardening support: {e}"),
            Self::Curve(e) => e.fmt(f),
            Self::Native(e) => e.fmt(f),
            Self::Invalid(e) => write!(f, "thermal hardening source admission: {e}"),
        }
    }
}
impl std::error::Error for ThermalHardeningSourceError {}
impl From<MatDbError> for ThermalHardeningSourceError {
    fn from(e: MatDbError) -> Self {
        Self::Source(e)
    }
}
impl From<PropertySupportError> for ThermalHardeningSourceError {
    fn from(e: PropertySupportError) -> Self {
        Self::Support(e)
    }
}
impl From<TabulatedError> for ThermalHardeningSourceError {
    fn from(e: TabulatedError) -> Self {
        Self::Curve(e)
    }
}
impl From<FiniteJ2ThermalError> for ThermalHardeningSourceError {
    fn from(e: FiniteJ2ThermalError) -> Self {
        Self::Native(e)
    }
}

fn require(accepted: bool, reason: &'static str) -> Result<(), ThermalHardeningSourceError> {
    if accepted {
        Ok(())
    } else {
        Err(ThermalHardeningSourceError::Invalid(reason))
    }
}

fn query(
    card: &MaterialCard,
    key: &PropertyKey,
    point: &QueryPoint,
    selection: ClaimSelection,
) -> Result<MaterialAnswer, ThermalHardeningSourceError> {
    let answer = match selection {
        ClaimSelection::Policy(policy) => card.claims().query_typed(key, point, policy)?,
        ClaimSelection::Pinned(pin) => card.claims().query_pinned_typed(key, point, pin)?,
    };
    card.claims().verify_receipt(&answer.receipt)?;
    Ok(answer)
}

fn same_slope(a: f64, b: f64) -> bool {
    a.is_finite() && b.is_finite() && (a - b).abs() <= 128.0 * f64::EPSILON * a.abs().max(b.abs())
}

fn coordinate_index(coordinates: &[f64], value: f64) -> Result<usize, ThermalHardeningSourceError> {
    // Source coordinates remain untouched in the receipts. Signed zero has one
    // physical grid position, consistent with the native strain admission.
    let value = if value == 0.0 { 0.0 } else { value };
    coordinates
        .binary_search_by(|probe| probe.total_cmp(&value))
        .map_err(|_| ThermalHardeningSourceError::Invalid("missing original grid coordinate"))
}

/// Prepare affine Kirchhoff hardening from an original complete Cartesian grid.
///
/// The context already contains the explicitly typed absolute reference
/// temperature and excludes accumulated strain. Both family endpoints must be
/// measured temperatures. Every intervening original temperature row constrains
/// the derivative; non-affine columns refuse rather than lose thermal entropy.
/// A reference between measured rows is explicit native-model interpolation.
/// The returned curve always spans the full original plastic-strain axis.
pub fn resolve_thermal_tabulated_hardening(
    card: &MaterialCard,
    context: &QueryPoint,
    request: &ThermalHardeningCurveRequest,
    selection: ClaimSelection,
) -> Result<ResolvedThermalTabulatedHardening, ThermalHardeningSourceError> {
    let (reference, slopes, source) = prepare_thermal_hardening(
        card,
        context,
        request,
        selection,
        YieldStressMeasure::Kirchhoff,
    )?;
    Ok(ResolvedThermalTabulatedHardening {
        hardening: FiniteJ2ThermalTabulatedHardening::new(reference, slopes)?,
        source,
    })
}

/// Prepare an original Cauchy stress surface for the small-strain thermal family.
/// The same complete-grid, affine-temperature, original-cell and receipt gates
/// apply. Only the exact Cauchy property can supply this source; a Kirchhoff
/// claim is never accepted as an alias or converted by assuming unit Jacobian.
pub fn resolve_thermal_cauchy_tabulated_hardening(
    card: &MaterialCard,
    context: &QueryPoint,
    request: &ThermalHardeningCurveRequest,
    selection: ClaimSelection,
) -> Result<ResolvedThermalCauchyTabulatedHardening, ThermalHardeningSourceError> {
    let (reference, slopes, source) = prepare_thermal_hardening(
        card,
        context,
        request,
        selection,
        YieldStressMeasure::Cauchy,
    )?;
    Ok(ResolvedThermalCauchyTabulatedHardening {
        hardening: AffineTabulatedHardening::new(reference, slopes)?,
        source,
    })
}

fn prepare_thermal_hardening(
    card: &MaterialCard,
    context: &QueryPoint,
    request: &ThermalHardeningCurveRequest,
    selection: ClaimSelection,
    stress_measure: YieldStressMeasure,
) -> Result<(TabulatedHardening, Vec<f64>, ThermalHardeningCurveSource), ThermalHardeningSourceError>
{
    let [lower_k, upper_k] = request.temperature_range_k;
    let reference_k = request.reference_temperature_k;
    require(
        !request.temperature_axis.trim().is_empty()
            && request.temperature_axis != ACCUMULATED_PLASTIC_STRAIN_AXIS
            && lower_k.is_finite()
            && upper_k.is_finite()
            && reference_k.is_finite()
            && lower_k > 0.0
            && upper_k > lower_k
            && reference_k >= lower_k
            && reference_k <= upper_k
            && request.maximum_samples >= 4
            && request.maximum_knots >= 2,
        "declare a positive temperature span, contained reference, distinct axis and bounded grid",
    )?;
    let temperature_quantity = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let strain_quantity = QuantitySpec::dimensional(Dims::NONE);
    require(
        context.axes().get(&request.temperature_axis) == Some(&reference_k)
            && context.axis_quantities().get(&request.temperature_axis)
                == Some(&temperature_quantity),
        "context must contain the declared semantic absolute reference temperature",
    )?;
    for (axis, &value) in context.axes() {
        if axis == &request.temperature_axis {
            continue;
        }
        let quantity = context.axis_quantities().get(axis).copied();
        require(
            axis != ACCUMULATED_PLASTIC_STRAIN_AXIS
                && axis != "T"
                && quantity != Some(temperature_quantity)
                && is_fixed_external_coordinate(axis, quantity, Some(value)),
            "only one evolving temperature and fixed external conditions are admitted",
        )?;
    }
    let origin = context
        .clone()
        .with(&request.temperature_axis, lower_k)?
        .with_quantity(ACCUMULATED_PLASTIC_STRAIN_AXIS, strain_quantity, 0.0)?;
    let key = PropertyKey::new(stress_measure.property_name(), Pressure::DIMS);
    let selected = query(card, &key, &origin, selection)?;
    let claim = card.claims().claim(selected.receipt.selected).ok_or(
        ThermalHardeningSourceError::Invalid("selected thermal hardening claim is absent"),
    )?;
    let PropertyValue::Surface2D {
        axes,
        axis_dims,
        samples,
        triangles,
        dims,
    } = &claim.value
    else {
        return Err(ThermalHardeningSourceError::Invalid(
            "thermal hardening needs an original strain-temperature surface",
        ));
    };
    let strain_axis = if axes[0] == ACCUMULATED_PLASTIC_STRAIN_AXIS {
        0
    } else {
        1
    };
    let temperature_axis = 1 - strain_axis;
    require(
        axes[strain_axis] == ACCUMULATED_PLASTIC_STRAIN_AXIS
            && axes[temperature_axis] == request.temperature_axis
            && axis_dims[strain_axis] == Dims::NONE
            && axis_dims[temperature_axis] == Temperature::DIMS
            && *dims == Pressure::DIMS
            && claim
                .validity
                .axis_quantities()
                .get(ACCUMULATED_PLASTIC_STRAIN_AXIS)
                == Some(&strain_quantity)
            && claim
                .validity
                .axis_quantities()
                .get(&request.temperature_axis)
                == Some(&temperature_quantity)
            && claim
                .validity
                .bounds()
                .contains_key(ACCUMULATED_PLASTIC_STRAIN_AXIS)
            && claim
                .validity
                .bounds()
                .contains_key(&request.temperature_axis),
        "original source axes need exact accumulated-strain and absolute-temperature meaning and bounds",
    )?;
    require(
        samples.len() >= 4 && samples.len() <= request.maximum_samples,
        "original sample count exceeds the preparation budget or cannot form a thermal grid",
    )?;
    for axis in claim.validity.bounds().keys() {
        if axis == ACCUMULATED_PLASTIC_STRAIN_AXIS || axis == &request.temperature_axis {
            continue;
        }
        require(
            is_fixed_external_coordinate(
                axis,
                context.axis_quantities().get(axis).copied(),
                context.axes().get(axis).copied(),
            ),
            "thermal hardening cannot freeze an evolving response, rate or history coordinate",
        )?;
    }
    let mut strains = Vec::with_capacity(samples.len());
    let mut temperatures = Vec::with_capacity(samples.len());
    for (point, value) in samples {
        require(
            point.iter().all(|coordinate| coordinate.is_finite()) && value.is_finite(),
            "thermal grid observations must be finite",
        )?;
        strains.push(if point[strain_axis] == 0.0 {
            0.0
        } else {
            point[strain_axis]
        });
        temperatures.push(point[temperature_axis]);
    }
    strains.sort_by(f64::total_cmp);
    strains.dedup();
    temperatures.sort_by(f64::total_cmp);
    temperatures.dedup();
    require(
        strains.len() >= 2
            && strains.len() <= request.maximum_knots
            && strains[0] == 0.0
            && temperatures.len() >= 2
            && temperatures[0] > 0.0
            && strains.len().checked_mul(temperatures.len()) == Some(samples.len()),
        "a complete original Cartesian grid with zero-origin strain and at least two temperatures is required",
    )?;
    let lower_row = coordinate_index(&temperatures, lower_k)?;
    let upper_row = coordinate_index(&temperatures, upper_k)?;
    let mut grid = vec![f64::NAN; samples.len()];
    for (point, value) in samples {
        let index = coordinate_index(&temperatures, point[temperature_axis])? * strains.len()
            + coordinate_index(&strains, point[strain_axis])?;
        require(grid[index].is_nan(), "duplicate original grid position")?;
        grid[index] = *value;
    }
    require(
        grid.iter().all(|value| value.is_finite()),
        "missing original Cartesian grid observation",
    )?;
    let row = |index: usize| -> Vec<HardeningPoint> {
        strains
            .iter()
            .enumerate()
            .map(|(column, &plastic_strain)| HardeningPoint {
                plastic_strain,
                yield_stress_pa: grid[index * strains.len() + column],
            })
            .collect()
    };
    // Every segment's endpoint differences are affine in temperature, so these
    // two original rows prove positive nondecreasing resistance throughout.
    TabulatedHardening::new(row(lower_row))?;
    TabulatedHardening::new(row(upper_row))?;
    let mut points = Vec::with_capacity(strains.len());
    let mut slopes = Vec::with_capacity(strains.len());
    for (column, &plastic_strain) in strains.iter().enumerate() {
        let line: Vec<_> = (lower_row..=upper_row)
            .map(|index| (temperatures[index], grid[index * strains.len() + column]))
            .collect();
        let reference_yield = match coordinate_index(&temperatures, reference_k) {
            Ok(index) => grid[index * strains.len() + column],
            Err(_) => {
                let first = line[0];
                let last = line[line.len() - 1];
                ((last.1 - first.1) / (last.0 - first.0)).mul_add(reference_k - first.0, first.1)
            }
        };
        let slope = affine_temperature_slope(&line, reference_k, reference_yield)
            .map_err(ThermalHardeningSourceError::Invalid)?;
        points.push(HardeningPoint {
            plastic_strain,
            yield_stress_pa: reference_yield,
        });
        slopes.push(slope);
    }
    let pin = ClaimSelection::Pinned(selected.receipt.selected);
    let continuous_support = match claim.interpolation {
        InterpolationPolicy::TabulatedOnly => {
            require(
                triangles.is_empty(),
                "the observational grid must not declare competing source interpolation cells",
            )?;
            None
        }
        InterpolationPolicy::LinearInside => {
            // Bilinear grid interpolation normally differs from affine source
            // triangles. This bounded route admits only one shared plane.
            let strain_line: Vec<_> = points
                .iter()
                .map(|point| (point.plastic_strain, point.yield_stress_pa))
                .collect();
            let strain_slope =
                affine_temperature_slope(&strain_line, 0.0, points[0].yield_stress_pa)
                    .map_err(ThermalHardeningSourceError::Invalid)?;
            require(
                slopes.iter().all(|&slope| same_slope(slope, slopes[0])),
                "source triangles cannot be replaced by bilinear hardening interpolation",
            )?;
            let prepared = claim.value.prepare_surface_response()?;
            for triangle in triangles {
                let mut center = [0.0; 2];
                for &vertex in triangle {
                    for axis in 0..2 {
                        center[axis] += samples[vertex as usize].0[axis] / 3.0;
                    }
                }
                let gradient = prepared.gradient_for_axes(center, [true, true])?;
                require(
                    same_slope(gradient[strain_axis], strain_slope)
                        && same_slope(gradient[temperature_axis], slopes[0]),
                    "original source cells must have the same admitted affine strain-temperature derivatives",
                )?;
            }
            let upper = context
                .clone()
                .with(&request.temperature_axis, upper_k)?
                .with_quantity(
                    ACCUMULATED_PLASTIC_STRAIN_AXIS,
                    strain_quantity,
                    strains[strains.len() - 1],
                )?;
            Some(
                card.claims()
                    .query_envelope_typed(&key, &origin, &upper, pin)?,
            )
        }
        _ => {
            return Err(ThermalHardeningSourceError::Invalid(
                "thermal hardening needs a discrete measured grid or compatible original affine cells",
            ));
        }
    };
    let hardening = TabulatedHardening::new(points)?;
    let mut observations = Vec::with_capacity(samples.len());
    for (point, value) in samples {
        let at = context
            .clone()
            .with(&request.temperature_axis, point[temperature_axis])?
            .with_quantity(
                ACCUMULATED_PLASTIC_STRAIN_AXIS,
                strain_quantity,
                point[strain_axis],
            )?;
        let answer = query(card, &key, &at, pin)?;
        require(
            answer.receipt.selected == selected.receipt.selected
                && answer.evidence.value.value == *value,
            "each original observation must retain its pinned claim and exact yield stress",
        )?;
        observations.push((at, answer));
    }
    Ok((
        hardening,
        slopes,
        ThermalHardeningCurveSource {
            card_identity: card.content_hash(),
            material_state: card.id().clone(),
            stress_measure,
            context: context.clone(),
            request: request.clone(),
            selection,
            selection_answer: selected,
            claim: claim.clone(),
            observations,
            continuous_support,
        },
    ))
}
