//! Prepare one original single-phase material source line once for trial fluxes.
use std::sync::Arc;
use super::{ConductivityCurve, DepositionError};
use fs_thermal::enthalpy::conduction::implicit::ImplicitConductionState;
use crate::phase::{EquilibriumEnthalpyPhaseCurve, PhaseStateError, SolidLiquidPhase};
use crate::thermal::source::{ThermalMaterialError, ThermalMaterialPoint};
use crate::state_point::{MaterialPropertySelection, ScalarAdmissibility,
    ScalarPropertyRequirement, THERMAL_CONDUCTIVITY_DIMS, resolve_material_state_point};
use fs_matdb::{ClaimSelection, PropertyKey, PropertyLineSlice, PropertySupportError};
use fs_qty::{QuantitySpec, semantic::{QuantityKind, SemanticType, ValueForm}};

/// Original source/caloric failures are retained rather than replaced by k=const.
#[derive(Clone, Debug, PartialEq)]
pub enum LayerSourceError {
    Invalid(&'static str),
    Phase(PhaseStateError),
    Material(ThermalMaterialError),
    Source(PropertySupportError),
    Thermal(DepositionError),
}
impl From<DepositionError> for LayerSourceError {
    fn from(error: DepositionError) -> Self { Self::Thermal(error) }
}
impl core::fmt::Display for LayerSourceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result { write!(f, "layer conduction: {self:?}") }
}
impl std::error::Error for LayerSourceError {}

/// Source-qualified k(T) on an explicit all-solid or all-liquid branch of the
/// SAME caloric curve. No source query, receipt construction or state hashing
/// occurs during a nonlinear flux trial. The original source line is retained
/// behind an Arc so accepted thermal payloads can share it without copying.
#[derive(Clone, Debug)]
pub struct ResolvedConductionLayer<'a> {
    caloric: &'a EquilibriumEnthalpyPhaseCurve,
    conductivity: ConductivityCurve,
    source: Arc<PropertyLineSlice>,
    phase: SolidLiquidPhase,
    enthalpy_range_j_kg: [f64; 2],
}
impl<'a> ResolvedConductionLayer<'a> {
    /// Resolve continuously supported k at fixed pressure/composition/etc. and
    /// fixed phase fraction. Every original source breakpoint is retained.
    /// Gaps, discrete-only data, changing source identity, hidden frequency
    /// dependence and source spans crossing the caloric phase boundary refuse.
    ///
    /// The source-line owner bounds admission to 4096 knots; it is currently a
    /// synchronous source query. Cancellation is checked around that operation
    /// and every 256 copied knots. Preparation never changes thermal history.
    #[allow(clippy::too_many_arguments)]
    pub fn new(caloric: &'a EquilibriumEnthalpyPhaseCurve,
        point: ThermalMaterialPoint<'_>, key: &PropertyKey, selection: ClaimSelection,
        phase: SolidLiquidPhase, temperature_range_k: [f64; 2],
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>) -> Result<Self, LayerSourceError> {
        checkpoint()?;
        let invalid = LayerSourceError::Invalid;
        let q = key.quantity();
        if q != QuantitySpec::dimensional(THERMAL_CONDUCTIVITY_DIMS)
            && q != QuantitySpec::semantic(SemanticType::new(QuantityKind::ThermalConductivity, ValueForm::Static)) {
            return Err(invalid("a scalar static thermal conductivity key is required"));
        }
        let fraction = match phase {
            SolidLiquidPhase::Solid => 0.0,
            SolidLiquidPhase::Liquid => 1.0,
            SolidLiquidPhase::SolidLiquid => return Err(invalid("a resolved conduction layer needs one pure phase")),
        };
        let knots = caloric.knots();
        // The caloric owner already admits monotone fractions. Find the pure
        // prefix/suffix logarithmically, without rescanning a large h chart.
        let (first, last) = if fraction == 0.0 {
            let end = knots.partition_point(|k| k.liquid_mass_fraction == 0.0);
            (0, end.checked_sub(1).ok_or(invalid("the caloric chart has no requested pure phase"))?)
        } else {
            (knots.partition_point(|k| k.liquid_mass_fraction < 1.0), knots.len() - 1)
        };
        let [lo, hi] = temperature_range_k;
        if first >= last || !lo.is_finite() || !hi.is_finite() || lo >= hi
            || lo < knots[first].temperature_k || hi > knots[last].temperature_k {
            return Err(invalid("conductivity support must remain inside a finite pure-phase temperature span"));
        }
        let enthalpy_range_j_kg = [
            source_enthalpy_bound(caloric, first, last, lo, false)?,
            source_enthalpy_bound(caloric, first, last, hi, true)?,
        ];
        if enthalpy_range_j_kg[0] >= enthalpy_range_j_kg[1]
            || !(enthalpy_range_j_kg[1] - enthalpy_range_j_kg[0]).is_finite() {
            return Err(invalid("conductivity source needs a representable original enthalpy span"));
        }
        let anchor = caloric.state_at_specific_enthalpy(knots[first].specific_enthalpy_j_kg)
            .map_err(LayerSourceError::Phase)?;
        let lower = point.at(anchor).map_err(LayerSourceError::Material)?
            .with(point.temperature_axis, lo)
            .map_err(|e| LayerSourceError::Material(ThermalMaterialError::Coordinate(e.to_string())))?;
        let upper = lower.clone().with(point.temperature_axis, hi)
            .map_err(|e| LayerSourceError::Material(ThermalMaterialError::Coordinate(e.to_string())))?;
        let line = point.material.claims().query_line_typed(key, point.temperature_axis,
            &lower, &upper, selection).map_err(LayerSourceError::Source)?;
        checkpoint()?;
        let first_answer = line.knots.first().ok_or(invalid("source line has no supported knot"))?;
        let requirement = ScalarPropertyRequirement::try_with_key(key, ScalarAdmissibility::StrictlyPositive)
            .map_err(|e| LayerSourceError::Material(ThermalMaterialError::Material(e)))?;
        let bundle = resolve_material_state_point(point.material, &lower, &[requirement],
            MaterialPropertySelection::PinnedByProperty(vec![(key.name().to_owned(), first_answer.1.receipt.selected)]))
            .map_err(|e| LayerSourceError::Material(ThermalMaterialError::Material(e)))?;
        let property = bundle.properties().first().ok_or(invalid("conductivity bundle is empty"))?;
        point.admit_static(property).map_err(LayerSourceError::Material)?;
        let mut values = Vec::with_capacity(line.knots.len());
        for (i, (t, answer)) in line.knots.iter().enumerate() {
            if i % 256 == 0 { checkpoint()?; }
            values.push([*t, answer.evidence.value.value]);
        }
        let conductivity = ConductivityCurve::new(values)?;
        checkpoint()?;
        Ok(Self { caloric, conductivity, source: Arc::new(line), phase,
            enthalpy_range_j_kg })
    }

    pub fn caloric_curve(&self) -> &'a EquilibriumEnthalpyPhaseCurve { self.caloric }
    pub fn conductivity_curve(&self) -> &ConductivityCurve { &self.conductivity }
    pub fn source_line(&self) -> &Arc<PropertyLineSlice> { &self.source }
    pub const fn phase(&self) -> SolidLiquidPhase { self.phase }

    /// Original enthalpy support common to the selected pure phase and k(T) source.
    pub const fn enthalpy_range_j_kg(&self) -> [f64; 2] {
        self.enthalpy_range_j_kg
    }

    /// Original chart temperature, with the prepared phase/source domain enforced.
    /// A melting plateau cannot be admitted merely because its T is unchanged.
    pub fn temperature_at_enthalpy(&self, h: f64) -> Result<f64, LayerSourceError> {
        if !h.is_finite() || !(self.enthalpy_range_j_kg[0]..=self.enthalpy_range_j_kg[1]).contains(&h) {
            return Err(LayerSourceError::Invalid("layer endpoint left its prepared pure-phase branch"));
        }
        // Preparation inverted these source limits on the original chart and
        // checked the round trip. Retain exact source temperatures at the two
        // resulting h bounds, so inverse rounding cannot create extrapolation.
        let bounds = self.conductivity.temperature_range_k();
        let t = if h == self.enthalpy_range_j_kg[0] {
            bounds[0]
        } else if h == self.enthalpy_range_j_kg[1] {
            bounds[1]
        } else {
            self.caloric.temperature_at_specific_enthalpy(h).map_err(LayerSourceError::Phase)?
        };
        self.conductivity.conductivity_w_m_k(t)?;
        Ok(t)
    }

    /// Original caloric inverse and analytic k(T) tangent in the same h datum.
    /// At a restricted upper endpoint use the admitted left Cp, including a
    /// melting-onset endpoint whose unused right interval is a latent plateau.
    /// This remains one pure-phase layer; it supplies no mixed-phase profile.
    pub fn conduction_state_at_enthalpy(
        &self, h: f64,
    ) -> Result<ImplicitConductionState, DepositionError> {
        let temperature = self.temperature_at_enthalpy(h)
            .map_err(|_| DepositionError::MaterialDomain)?;
        let curve = self.caloric;
        let temperature_slope = if h == self.enthalpy_range_j_kg[1] {
            let knots = curve.knots();
            let upper = knots.partition_point(|knot| knot.specific_enthalpy_j_kg < h)
                .min(knots.len() - 1);
            let lower = upper.checked_sub(1).ok_or(DepositionError::MaterialDomain)?;
            if h == knots[upper].specific_enthalpy_j_kg {
                curve.sensible_heat_capacities_on_interval(lower).map_or_else(
                    || (knots[upper].temperature_k - knots[lower].temperature_k)
                        / (knots[upper].specific_enthalpy_j_kg - knots[lower].specific_enthalpy_j_kg),
                    |[_, capacity]| capacity.recip(),
                )
            } else {
                curve.temperature_derivative_at_specific_enthalpy(h)
                    .map_err(|_| DepositionError::MaterialDomain)?
            }
        } else {
            curve.temperature_derivative_at_specific_enthalpy(h)
                .map_err(|_| DepositionError::MaterialDomain)?
        };
        let (conductivity, conductivity_temperature_slope) = self.conductivity
            .conductivity_and_temperature_derivative(temperature)?;
        let conductivity_slope = conductivity_temperature_slope * temperature_slope;
        if !temperature_slope.is_finite() || temperature_slope < 0.0
            || !conductivity_slope.is_finite()
            || (conductivity_temperature_slope != 0.0 && temperature_slope != 0.0
                && conductivity_slope == 0.0) {
            return Err(DepositionError::Arithmetic);
        }
        Ok(ImplicitConductionState {
            temperature_k: temperature,
            temperature_slope_kg_k_j: temperature_slope,
            conductivity_w_m_k: conductivity,
            conductivity_slope_kg_w_m_k_j: conductivity_slope,
        })
    }
}

// Invert each k-source endpoint once, preserving original latent coordinates:
// the lower temperature keeps the first h and the upper keeps the last h on
// an equal-temperature interval of the admitted pure phase.
fn source_enthalpy_bound(
    curve: &EquilibriumEnthalpyPhaseCurve, first: usize, last: usize,
    temperature: f64, upper: bool,
) -> Result<f64, LayerSourceError> {
    let knots = &curve.knots()[first..=last];
    let end = knots.partition_point(|knot| knot.temperature_k <= temperature);
    if end > 0 && knots[end - 1].temperature_k == temperature {
        let index = if upper { end - 1 } else {
            knots.partition_point(|knot| knot.temperature_k < temperature)
        };
        return Ok(knots[index].specific_enthalpy_j_kg);
    }
    let lower = end.checked_sub(1).ok_or(LayerSourceError::Invalid(
        "conductivity temperature has no original pure-phase enthalpy",
    ))?;
    let next = knots.get(end).ok_or(LayerSourceError::Invalid(
        "conductivity temperature has no original pure-phase enthalpy",
    ))?;
    let h = curve.specific_enthalpy_on_sensible_interval(first + lower, temperature)
        .map_err(LayerSourceError::Phase)?;
    if !h.is_finite() || h <= knots[lower].specific_enthalpy_j_kg
        || h >= next.specific_enthalpy_j_kg {
        return Err(LayerSourceError::Invalid(
            "conductivity temperature has no representable interior enthalpy",
        ));
    }
    let original = curve.values_at_specific_enthalpy(h).map_err(LayerSourceError::Phase)?;
    if (original.temperature_k - temperature).abs()
        > 128.0 * f64::EPSILON * temperature.max(original.temperature_k) {
        return Err(LayerSourceError::Invalid(
            "source enthalpy cannot reproduce the original conductivity temperature",
        ));
    }
    Ok(h)
}
