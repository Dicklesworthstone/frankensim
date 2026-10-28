//! Source-backed sensible-enthalpy transport through material-bound channels.
//!
//! The existing phase owner supplies a strictly liquid, piecewise-linear T(h)
//! chart. Within each caloric interval the existing constant-cp `AirPath` march
//! is exact; this driver stops at source knots before continuing with the next
//! slope. Heat is mdot*(h_out-h_in), never a frozen cp times total delta T.
//!
//! Hydraulic mass flow and each region's convection coefficient remain frozen
//! at their original material query. This changes the caloric transport law,
//! not the momentum solve, film-property iteration, boiling model or source
//! uncertainty. The Cp-to-enthalpy chart's interpolation limits still apply.

pub mod coupled;

use fs_convection::CorrelationRobinBoundary;
use fs_exec::Cx;
use fs_material::phase::{
    EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState,
    PhaseStateError, ResolvedHeatCapacityCurve,
};
use fs_material::state_point::{DENSITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_PROPERTY};
use fs_math::det;
use fs_qty::{MassFlowRate, Temperature};

use crate::conjugate::{AirPath, AirSegment};
use super::{MaterialFlowPath, MaterialHeatExchange};
use super::hydraulic::{MaterialFlowResistance, ResistanceError, ResistanceEvaluation};
use super::super::{FrozenFlowMedium, MediumFlowError};

/// Source mismatch, unsupported thermal trajectory, or original owner refusal.
#[derive(Debug)]
pub enum LiquidEnthalpyError {
    /// Unsupported source, geometry, or unrepresentable arithmetic.
    Invalid(&'static str),
    /// Original material phase/chart refusal.
    Phase(PhaseStateError),
    /// Original channel/correlation refusal.
    Flow(MediumFlowError),
    /// Original pressure-loss owner refused this region.
    Resistance { region: String, source: ResistanceError },
    /// Cancellation publishes no partial march or boundary rows.
    Cancelled,
}
impl core::fmt::Display for LiquidEnthalpyError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "liquid enthalpy flow: {self:?}")
    }
}
impl std::error::Error for LiquidEnthalpyError {}
impl From<PhaseStateError> for LiquidEnthalpyError {
    fn from(value: PhaseStateError) -> Self { Self::Phase(value) }
}
impl From<crate::AirflowError> for LiquidEnthalpyError {
    fn from(value: crate::AirflowError) -> Self { Self::Flow(MediumFlowError::Airflow(value)) }
}

/// One accepted region, retaining its actual caloric states and convection law.
#[derive(Clone, Debug, PartialEq)]
pub struct LiquidEnthalpyRegion<'a> {
    region: String,
    inlet: EquilibriumPhaseState,
    outlet: EquilibriumPhaseState,
    heat_rate_w: f64,
    reference_temperature_k: f64,
    exchange: MaterialHeatExchange<'a>,
    hydraulic_loss: Option<ResistanceEvaluation>,
}
impl<'a> LiquidEnthalpyRegion<'a> {
    /// Original solid trace name.
    #[must_use]
    pub fn region(&self) -> &str { &self.region }
    /// Source-bound entering state.
    #[must_use]
    pub const fn inlet(&self) -> EquilibriumPhaseState { self.inlet }
    /// Source-bound leaving state.
    #[must_use]
    pub const fn outlet(&self) -> EquilibriumPhaseState { self.outlet }
    /// Signed heat received from this wall, W; hydraulic heat is separate.
    #[must_use]
    pub const fn heat_rate_w(&self) -> f64 { self.heat_rate_w }
    /// Original hydraulic evaluation when pressure-loss heating was requested.
    #[must_use]
    pub const fn hydraulic_loss(&self) -> Option<&ResistanceEvaluation> { self.hydraulic_loss.as_ref() }
    /// Pressure work converted to fluid heat in this region, W.
    #[must_use]
    pub fn hydraulic_heat_rate_w(&self) -> f64 {
        self.hydraulic_loss.as_ref().map_or(0.0, |loss| loss.dissipated_power_w)
    }
    /// Effective Robin reference reproducing this region's actual heat, K.
    #[must_use]
    pub const fn reference_temperature_k(&self) -> f64 { self.reference_temperature_k }
    /// Original frozen transport and correlation source.
    #[must_use]
    pub const fn exchange(&self) -> &MaterialHeatExchange<'a> { &self.exchange }
    /// Lower exactly this heat exchange to the existing conduction boundary.
    /// The caller cannot substitute another reference, conductivity or length.
    pub fn robin_boundary(&self) -> Result<CorrelationRobinBoundary, MediumFlowError> {
        self.exchange.robin_boundary(Temperature::new(self.reference_temperature_k))
    }
}

/// A complete conservative sensible-enthalpy march, including the source chart.
#[derive(Clone, Debug, PartialEq)]
pub struct LiquidEnthalpyMarch<'a, 'c> {
    source: &'c ResolvedHeatCapacityCurve,
    regions: Vec<LiquidEnthalpyRegion<'a>>,
    mass_flow_kg_per_s: f64,
    total_heat_rate_w: f64,
    hydraulic_heat_rate_w: f64,
    enthalpy_gain_w: f64,
    balance_residual_w: f64,
}
impl<'a, 'c> LiquidEnthalpyMarch<'a, 'c> {
    /// Cp/density source and every receipt used to build the caloric chart.
    #[must_use]
    pub const fn source(&self) -> &'c ResolvedHeatCapacityCurve { self.source }
    /// Results in exactly the material path's region order.
    #[must_use]
    pub fn regions(&self) -> &[LiquidEnthalpyRegion<'a>] { &self.regions }
    /// One conserved mass-flow rate through every region, kg/s.
    #[must_use]
    pub const fn mass_flow_kg_per_s(&self) -> f64 { self.mass_flow_kg_per_s }
    /// Sum of signed wall heat transfers, W.
    #[must_use]
    pub const fn total_heat_rate_w(&self) -> f64 { self.total_heat_rate_w }
    /// Total original pressure-loss power deposited in the fluid, W.
    #[must_use]
    pub const fn hydraulic_heat_rate_w(&self) -> f64 { self.hydraulic_heat_rate_w }
    /// Actual `mdot*(h_out-h_in)` from the source-backed endpoint states, W.
    /// This independent observation is not reconstructed from heat inputs.
    #[must_use]
    pub const fn enthalpy_gain_w(&self) -> f64 { self.enthalpy_gain_w }
    /// Actual advected enthalpy minus reported wall and hydraulic heat, W.
    /// This exposes representation/reduction error; it is not physical heat.
    #[must_use]
    pub const fn balance_residual_w(&self) -> f64 { self.balance_residual_w }
    /// Source-bound final fluid state.
    #[must_use]
    pub fn outlet(&self) -> EquilibriumPhaseState {
        self.regions.last().expect("admitted nonempty material path").outlet
    }
}

impl<'a> MaterialFlowPath<'a> {
    /// March a source-resolved liquid using its temperature-dependent enthalpy.
    ///
    /// The chart must select the SAME Cp/density claims and fixed nonthermal
    /// source coordinates as the flow medium. The named temperature axis must
    /// be present in both. Solid/mixed charts and unrelated cards refuse. Each
    /// fluid state must stay inside the chart; a hotter/cooler wall may be
    /// outside it as long as the resulting fluid trajectory is supported.
    ///
    /// Within each piece the original `AirPath` exponential exchange is used.
    /// At most one traversal of the chart is needed per isothermal wall region.
    /// Cancellation is checked at source receipts, regions and knot crossings.
    /// The original path, source and any previous march remain unchanged on error.
    pub fn march_liquid_enthalpy<'c>(
        &self, cx: &Cx<'_>, source: &'c ResolvedHeatCapacityCurve,
        temperature_axis: &str, wall_temperatures_k: &[f64],
    ) -> Result<LiquidEnthalpyMarch<'a, 'c>, LiquidEnthalpyError> {
        self.march_liquid_enthalpy_inner(cx, source, temperature_axis, wall_temperatures_k, None)
    }

    /// Include the original `dp*Q` loss as uniformly distributed fluid heat.
    ///
    /// Choosing this reduced model declares equal elevation, equal flow areas,
    /// unchanged bulk kinetic energy, constant source density, and complete local
    /// conversion of each declared pressure loss into fluid heat. It does not
    /// model gas throttling, pump efficiency or heat retained by a porous skeleton.
    /// Resistance and convection remain at the original frozen material query.
    ///
    /// Within each caloric piece, the same exact wall march integrates
    /// `mdot dh/dA = htc*(Twall-T(h)) + P_loss/A_region` by shifting its driving
    /// temperature. Hydraulic heat is never charged to the wall's Robin row.
    ///
    /// # Errors
    /// Original source/correlation/Reynolds refusals, incompatible geometry,
    /// nonconstant density, unsupported caloric trajectory, or cancellation.
    pub fn march_liquid_enthalpy_with_dissipation<'c>(
        &self, cx: &Cx<'_>, source: &'c ResolvedHeatCapacityCurve,
        temperature_axis: &str, wall_temperatures_k: &[f64],
        resistances: &[MaterialFlowResistance],
    ) -> Result<LiquidEnthalpyMarch<'a, 'c>, LiquidEnthalpyError> {
        self.march_liquid_enthalpy_inner(cx, source, temperature_axis,
            wall_temperatures_k, Some(resistances))
    }

    fn march_liquid_enthalpy_inner<'c>(
        &self, cx: &Cx<'_>, source: &'c ResolvedHeatCapacityCurve,
        temperature_axis: &str, wall_temperatures_k: &[f64],
        resistances: Option<&[MaterialFlowResistance]>,
    ) -> Result<LiquidEnthalpyMarch<'a, 'c>, LiquidEnthalpyError> {
        checkpoint(cx)?;
        self.admit_caloric_source(cx, source, temperature_axis)?;
        if let Some(resistances) = resistances {
            let FrozenFlowMedium::Liquid(liquid) = self.medium else { unreachable!("admitted liquid") };
            admit_dissipative_density(source, liquid.density_kg_m3())?;
            let area = self.exchanges[0].geometry.flow_area;
            if resistances.len() != self.exchanges.len()
                || resistances.iter().zip(&self.exchanges).any(|(r, e)|
                    r.flow_area() != e.geometry.flow_area || e.geometry.flow_area != area)
            {
                return Err(LiquidEnthalpyError::Invalid("dissipative liquid heat requires one resistance per region and equal flow areas"));
            }
        }
        if wall_temperatures_k.len() != self.path.segments().len()
            || wall_temperatures_k.iter().any(|&t| !positive(t))
        {
            return Err(LiquidEnthalpyError::Invalid("one positive finite wall temperature is required per region"));
        }
        let curve = source.curve();
        let mut state = state_at_temperature(curve, self.path.inlet_temperature_k())?;
        let initial_enthalpy = state.specific_enthalpy_j_kg();
        let mut regions = Vec::with_capacity(self.exchanges.len());
        let mut total = 0.0;
        let mut hydraulic_total = 0.0;
        for (index, ((segment, exchange), &wall)) in self.path.segments().iter().zip(&self.exchanges).zip(wall_temperatures_k).enumerate() {
            checkpoint(cx)?;
            let hydraulic_loss = if let Some(resistances) = resistances {
                Some(resistances[index].evaluate(&self.medium, MassFlowRate::new(self.mass_flow_kg_per_s))
                    .map_err(|source| LiquidEnthalpyError::Resistance { region: segment.region().to_owned(), source })?)
            } else { None };
            let power = hydraulic_loss.as_ref().map_or(0.0, |loss| loss.dissipated_power_w);
            let outlet = march_region_with_heat(cx, curve, state, self.mass_flow_kg_per_s, segment, wall, power)?;
            let enthalpy_gain = self.mass_flow_kg_per_s
                * (outlet.specific_enthalpy_j_kg() - state.specific_enthalpy_j_kg());
            let heat = enthalpy_gain - power;
            let reference = wall - heat / (segment.htc_w_per_m2_k() * segment.area_m2());
            total += heat;
            hydraulic_total += power;
            if !heat.is_finite() || !total.is_finite() || !hydraulic_total.is_finite()
                || !(total + hydraulic_total).is_finite() || !positive(reference) {
                return Err(LiquidEnthalpyError::Invalid("enthalpy transfer or effective Robin reference is unrepresentable"));
            }
            regions.push(LiquidEnthalpyRegion {
                region: segment.region().to_owned(), inlet: state, outlet,
                heat_rate_w: heat, reference_temperature_k: reference, exchange: exchange.clone(), hydraulic_loss,
            });
            state = outlet;
        }
        checkpoint(cx)?;
        let enthalpy_gain_w = self.mass_flow_kg_per_s * (state.specific_enthalpy_j_kg() - initial_enthalpy);
        let (partial, first_tail) = fs_math::eft::two_sum(enthalpy_gain_w, -total);
        let (rounded, second_tail) = fs_math::eft::two_sum(partial, -hydraulic_total);
        let balance_residual_w = rounded + (first_tail + second_tail);
        if !enthalpy_gain_w.is_finite() || !balance_residual_w.is_finite() {
            return Err(LiquidEnthalpyError::Invalid("advected enthalpy or heat balance is unrepresentable"));
        }
        Ok(LiquidEnthalpyMarch { source, regions, mass_flow_kg_per_s: self.mass_flow_kg_per_s,
            total_heat_rate_w: total, hydraulic_heat_rate_w: hydraulic_total,
            enthalpy_gain_w, balance_residual_w })
    }

    fn admit_caloric_source(
        &self, cx: &Cx<'_>, source: &ResolvedHeatCapacityCurve, axis: &str,
    ) -> Result<(), LiquidEnthalpyError> {
        let FrozenFlowMedium::Liquid(liquid) = self.medium else {
            return Err(LiquidEnthalpyError::Invalid("a solid-liquid caloric chart cannot represent a gas"));
        };
        let material = liquid.material();
        let curve = source.curve();
        if curve.material_card_identity() != material.card_identity()
            || curve.knots().iter().any(|k| k.liquid_mass_fraction != 1.0)
            || curve.knots().windows(2).any(|k| k[1].temperature_k <= k[0].temperature_k)
        {
            return Err(LiquidEnthalpyError::Invalid("caloric chart must be fully liquid and belong to the flow material card"));
        }
        for receipt in source.receipts() {
            checkpoint(cx)?;
            if receipt.property != DENSITY_PROPERTY && receipt.property != SPECIFIC_HEAT_CAPACITY_PROPERTY {
                return Err(LiquidEnthalpyError::Invalid("caloric chart uses different source property keys"));
            }
            let property = material.property(&receipt.property)
                .ok_or(LiquidEnthalpyError::Invalid("flow material lacks the caloric property"))?;
            let original = &property.answer().receipt;
            if receipt.selected != original.selected
                || receipt.axis_quantities != original.axis_quantities
                || !receipt.query_point.iter().any(|(name, _)| name.as_str() == axis)
                || !original.query_point.iter().any(|(name, _)| name.as_str() == axis)
                || !receipt.query_point.iter().filter(|(name, _)| name.as_str() != axis)
                    .eq(original.query_point.iter().filter(|(name, _)| name.as_str() != axis))
            {
                return Err(LiquidEnthalpyError::Invalid("caloric and transport sources must select identical claims and nonthermal query conditions"));
            }
        }
        Ok(())
    }
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), LiquidEnthalpyError> {
    cx.checkpoint().map_err(|_| LiquidEnthalpyError::Cancelled)
}
fn positive(value: f64) -> bool { value.is_finite() && value > 0.0 }

fn admit_dissipative_density(source: &ResolvedHeatCapacityCurve, density: f64) -> Result<(), LiquidEnthalpyError> {
    if source.curve().knots().iter().any(|k| k.bulk_density_kg_m3 != density) {
        return Err(LiquidEnthalpyError::Invalid("hydraulic heating requires constant source density over the caloric span"));
    }
    Ok(())
}

fn march_region_with_heat(
    cx: &Cx<'_>, curve: &EquilibriumEnthalpyPhaseCurve, state: EquilibriumPhaseState,
    mass_flow: f64, segment: &AirSegment, wall: f64, hydraulic_power_w: f64,
) -> Result<EquilibriumPhaseState, LiquidEnthalpyError> {
    if !hydraulic_power_w.is_finite() || hydraulic_power_w < 0.0 {
        return Err(LiquidEnthalpyError::Invalid("hydraulic heat must be finite and nonnegative"));
    }
    let hydraulic_rise = hydraulic_power_w / (segment.htc_w_per_m2_k() * segment.area_m2());
    let driving = wall + hydraulic_rise;
    if hydraulic_power_w != 0.0 && (hydraulic_rise == 0.0 || driving == wall) {
        return Err(LiquidEnthalpyError::Invalid("nonzero hydraulic driving-temperature increment is unrepresentable"));
    }
    if !positive(driving) {
        return Err(LiquidEnthalpyError::Invalid("hydraulic-heating driving temperature is unrepresentable"));
    }
    march_region(cx, curve, state, mass_flow, segment, driving)
}

fn state_at_temperature(
    curve: &EquilibriumEnthalpyPhaseCurve, temperature: f64,
) -> Result<EquilibriumPhaseState, LiquidEnthalpyError> {
    let knots = curve.knots();
    if !temperature.is_finite() || temperature < knots[0].temperature_k
        || temperature > knots[knots.len() - 1].temperature_k
    {
        return Err(LiquidEnthalpyError::Invalid("fluid temperature is outside the caloric source span"));
    }
    let index = knots.partition_point(|k| k.temperature_k < temperature);
    if index < knots.len() && knots[index].temperature_k == temperature {
        return Ok(curve.state_at_specific_enthalpy(knots[index].specific_enthalpy_j_kg)?);
    }
    let lo = knots[index - 1];
    let hi = knots[index];
    let fraction = (temperature - lo.temperature_k) / (hi.temperature_k - lo.temperature_k);
    Ok(curve.state_at_specific_enthalpy(
        lo.specific_enthalpy_j_kg + fraction * (hi.specific_enthalpy_j_kg - lo.specific_enthalpy_j_kg),
    )?)
}

fn march_region(
    cx: &Cx<'_>, curve: &EquilibriumEnthalpyPhaseCurve, mut state: EquilibriumPhaseState,
    mass_flow: f64, segment: &AirSegment, wall: f64,
) -> Result<EquilibriumPhaseState, LiquidEnthalpyError> {
    let knots = curve.knots();
    let mut area = segment.area_m2();
    for _ in 0..knots.len() {
        checkpoint(cx)?;
        let temperature = state.temperature_k();
        if temperature == wall { return Ok(state); }
        let heating = wall > temperature;
        let enthalpy = state.specific_enthalpy_j_kg();
        let upper = if heating {
            knots.partition_point(|k| k.specific_enthalpy_j_kg <= enthalpy)
        } else {
            knots.partition_point(|k| k.specific_enthalpy_j_kg < enthalpy)
        };
        if upper == 0 || upper == knots.len() {
            return Err(LiquidEnthalpyError::Invalid("heat exchange would leave the source-supported liquid span"));
        }
        let lo = knots[upper - 1];
        let hi = knots[upper];
        let cp = (hi.specific_enthalpy_j_kg - lo.specific_enthalpy_j_kg) / (hi.temperature_k - lo.temperature_k);
        let local = AirPath::new(temperature, mass_flow, cp,
            vec![AirSegment::new(segment.region(), area, segment.htc_w_per_m2_k())?])?;
        let step = local.march(&[wall])?;
        let increment = step.total_heat_rate_w / mass_flow;
        let trial_enthalpy = enthalpy + increment;
        // The original exponential owns a nonzero physical transfer even when
        // division underflows or a large caloric datum hides its addition. Do
        // not accept unchanged fluid storage and turn that loss into wall heat.
        if step.total_heat_rate_w != 0.0 && (increment == 0.0 || trial_enthalpy == enthalpy) {
            return Err(LiquidEnthalpyError::Invalid("nonzero specific-enthalpy increment is unrepresentable"));
        }
        if !trial_enthalpy.is_finite() {
            return Err(LiquidEnthalpyError::Invalid("specific-enthalpy advance is unrepresentable"));
        }
        if (lo.specific_enthalpy_j_kg..=hi.specific_enthalpy_j_kg).contains(&trial_enthalpy) {
            return Ok(curve.state_at_specific_enthalpy(trial_enthalpy)?);
        }
        let boundary = if heating { hi } else { lo };
        let consumed = area_to_knot(temperature, wall, boundary, mass_flow, cp, segment.htc_w_per_m2_k())?;
        let remaining = area - consumed;
        if !positive(consumed) || !positive(remaining) || remaining >= area {
            return Err(LiquidEnthalpyError::Invalid("caloric knot crossing is not representable in the remaining exchange area"));
        }
        area = remaining;
        state = curve.state_at_specific_enthalpy(boundary.specific_enthalpy_j_kg)?;
    }
    Err(LiquidEnthalpyError::Invalid("caloric interval traversal exceeded the chart size"))
}

fn area_to_knot(
    temperature: f64, wall: f64, knot: EnthalpyPhaseKnot, mass_flow: f64, cp: f64, htc: f64,
) -> Result<f64, LiquidEnthalpyError> {
    let driving = wall - temperature;
    let fraction = (knot.temperature_k - temperature) / driving;
    if !(fraction > 0.0 && fraction < 1.0) {
        return Err(LiquidEnthalpyError::Invalid("the next caloric knot is not between the fluid and wall temperatures"));
    }
    // Invert the same exponential exchange. A short local series avoids
    // cancellation in -ln(1-x); its omitted term is <= x^6/[6(1-x)].
    let ntu = if fraction < 1e-4 {
        fraction * (1.0 + fraction * (0.5 + fraction * (1.0/3.0 + fraction * (0.25 + fraction/5.0))))
    } else {
        -det::ln((wall - knot.temperature_k) / driving)
    };
    let area = ((mass_flow * cp) / htc) * ntu;
    if !positive(area) {
        return Err(LiquidEnthalpyError::Invalid("area to caloric knot is unrepresentable"));
    }
    Ok(area)
}
