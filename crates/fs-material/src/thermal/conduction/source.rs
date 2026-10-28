//! Prepare one original single-phase material source line once for trial fluxes.
use std::sync::Arc;
use super::{ConductivityCurve, DepositionError};
use crate::phase::{EquilibriumEnthalpyPhaseCurve, PhaseStateError, SolidLiquidPhase};
use crate::thermal::source::{ThermalMaterialError, ThermalMaterialPoint};
use crate::state_point::{MaterialPropertySelection, ScalarAdmissibility,
    ScalarPropertyRequirement, THERMAL_CONDUCTIVITY_DIMS, resolve_material_state_point};
use fs_matdb::{ClaimSelection, PropertyKey, PropertyLineSlice, PropertySupportError};
use fs_qty::{QuantitySpec, semantic::{QuantityKind, SemanticType, ValueForm}};

/// Original source/caloric failures are retained rather than replaced by k=const.
#[derive(Debug)]
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
        point: ThermalMaterialPoint<'a>, key: &PropertyKey, selection: ClaimSelection,
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
            enthalpy_range_j_kg: [knots[first].specific_enthalpy_j_kg, knots[last].specific_enthalpy_j_kg] })
    }

    pub fn caloric_curve(&self) -> &'a EquilibriumEnthalpyPhaseCurve { self.caloric }
    pub fn conductivity_curve(&self) -> &ConductivityCurve { &self.conductivity }
    pub fn source_line(&self) -> &Arc<PropertyLineSlice> { &self.source }
    pub const fn phase(&self) -> SolidLiquidPhase { self.phase }

    /// Original chart temperature, with the prepared phase/source domain enforced.
    /// A melting plateau cannot be admitted merely because its T is unchanged.
    pub fn temperature_at_enthalpy(&self, h: f64) -> Result<f64, LayerSourceError> {
        if !h.is_finite() || !(self.enthalpy_range_j_kg[0]..=self.enthalpy_range_j_kg[1]).contains(&h) {
            return Err(LayerSourceError::Invalid("layer endpoint left its prepared pure-phase branch"));
        }
        let t = self.caloric.temperature_at_specific_enthalpy(h).map_err(LayerSourceError::Phase)?;
        self.conductivity.conductivity_w_m_k(t)?;
        Ok(t)
    }
}
