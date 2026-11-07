//! The same sheet specimen in the existing bounded enthalpy/transport owners.
//! Apparent reference density is ONLY mu/t, not a sourced solid-fiber density.
use super::*;
use crate::phase::{EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState,
    HeatCapacityKnot, PhaseStateError, SolidLiquidPhase};
use crate::state_point::{SPECIFIC_HEAT_CAPACITY_PROPERTY as CP,
    SPECIFIC_HEAT_CAPACITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY as K,
    THERMAL_CONDUCTIVITY_DIMS, ResolvedScalarProperty};
use fs_matdb::{ClaimId, ClaimSelection, MatDbError, PropertyKey, PropertySupportError};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

/// A source query, geometry conversion or caloric-chart refusal.
#[derive(Debug)]
pub enum SheetThermalError {
    /// Original common source/constitutive refusal.
    Sheet(SheetError),
    /// Original query-coordinate admission.
    Query(MatDbError),
    /// Original continuous source-line refusal, including its witness context.
    Support(Box<PropertySupportError>),
    /// Original bounded enthalpy-chart construction/evaluation.
    Phase(PhaseStateError),
    /// A missing, ambiguous or unsupported physical representation.
    Invalid(&'static str),
}
impl From<SheetError> for SheetThermalError {fn from(e:SheetError)->Self{Self::Sheet(e)}}
impl From<MaterialStatePointError> for SheetThermalError {fn from(e:MaterialStatePointError)->Self{Self::Sheet(e.into())}}
impl From<MatDbError> for SheetThermalError {fn from(e:MatDbError)->Self{Self::Query(e)}}
impl From<PropertySupportError> for SheetThermalError {fn from(e:PropertySupportError)->Self{Self::Support(Box::new(e))}}
impl From<PhaseStateError> for SheetThermalError {fn from(e:PhaseStateError)->Self{Self::Phase(e)}}
impl core::fmt::Display for SheetThermalError {
    fn fmt(&self,f:&mut core::fmt::Formatter<'_>)->core::fmt::Result {write!(f,"sheet thermal: {self:?}")}
}
impl std::error::Error for SheetThermalError {}

/// Restore ALL original conditioning coordinates and their quantity meanings;
/// change only T. No conditions or replacement material are guessed.
fn point_at(sheet:&ResolvedSheet,temperature:f64)->Result<QueryPoint,SheetThermalError> {
    if !temperature.is_finite() || temperature<=0.0 {return Err(SheetThermalError::Invalid("positive absolute temperature required"));}
    let absolute=QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,ValueForm::Static));
    if sheet.material().axis_quantities().get("T").is_some_and(|&q|
        q!=absolute && q!=QuantitySpec::dimensional(fs_qty::Temperature::DIMS)) {
        return Err(SheetThermalError::Invalid("sheet thermal T must retain its absolute-temperature meaning"));
    }
    Ok(sheet.material().query().with("T",temperature)?)
}
fn sample(sheet:&ResolvedSheet,card:&MaterialCard,claim:ClaimId,key:&'static str,dims:Dims,
    temperature:f64)->Result<ResolvedMaterialStatePoint,SheetThermalError> {
    if sheet.material().card_identity()!=card.content_hash() {
        return Err(SheetThermalError::Invalid("thermal source must be the original sheet material card"));
    }
    let point=point_at(sheet,temperature)?;
    let requirement=ScalarPropertyRequirement::try_new(key,dims,ScalarAdmissibility::StrictlyPositive)?;
    Ok(resolve_material_state_point(card,&point,&[requirement],
        MaterialPropertySelection::PinnedByProperty(vec![(key.into(),claim)]))?)
}

// Existing characterized scalar/T-only records retain their conditioning
// selectors and applicability. Additional response axes need the common typed
// held-environment policy; this owner evolves temperature alone.
fn thermal_source_supported(property:&ResolvedScalarProperty)->bool {
    (property.source_axes().iter().all(|(axis,_)|axis=="T")
        && property.is_constant_at_fixed_temperature())
        || property.has_only_fixed_external_conditions(Some("T"),false)
}

/// A source-retaining single-phase sensible-heat chart for the original sheet.
/// Its curve implements the existing `fs_thermal::enthalpy::ThermalLaw`.
/// No integration/transport algorithm or transient storage is duplicated here.
#[derive(Clone,Debug)]
pub struct SheetCaloric {
    sheet:ResolvedSheet,
    curve:EquilibriumEnthalpyPhaseCurve,
    heat_capacity_samples:Vec<ResolvedMaterialStatePoint>,
}
impl SheetCaloric {
    /// Integrate the explicitly selected Cp source through the ORIGINAL chart
    /// constructor. The canonical source-line owner retains every crossed curve
    /// knot and sampled-surface edge in the selected Kelvin span. Source holes,
    /// source switching and exact-only spans refuse. Other conditioning
    /// coordinates and typed descriptors remain fixed; new environmental
    /// response axes require explicit static quantity meanings.
    ///
    /// The apparent density mu/t is a reference geometry conversion for chart
    /// bookkeeping ONLY. Actual thermal cells must use reference mu*area mass.
    /// The sheet's measured construction density is never replaced with a fiber
    /// density, and absent Cp never borrows a nearby catalog specimen's value.
    /// No latent phase transition, compression, moisture or moving-area work is
    /// inferred. The bound controls the existing Cp-to-h interpolation, not data
    /// uncertainty or time-discretization error.
    pub fn try_new(sheet:ResolvedSheet,card:&MaterialCard,heat_capacity_claim:ClaimId,
        temperature_k:[f64;2],reference_h_j_kg:f64,maximum_chart_error_j_kg:f64)
        ->Result<Self,SheetThermalError> {
        let [lo,hi]=temperature_k;
        if !lo.is_finite() || !hi.is_finite() || lo<=0.0 || hi<=lo {
            return Err(SheetThermalError::Invalid("finite increasing positive Kelvin span required"));
        }
        let lower=sample(&sheet,card,heat_capacity_claim,CP,SPECIFIC_HEAT_CAPACITY_DIMS,lo)?;
        if !thermal_source_supported(lower.property(CP)
            .ok_or(SheetThermalError::Invalid("missing resolved Cp"))?) {
            return Err(SheetThermalError::Invalid("Cp needs a temperature path at supported held conditions, not a frozen foreign-axis law"));
        }
        let line=card.claims().query_line_typed(
            &PropertyKey::new(CP,SPECIFIC_HEAT_CAPACITY_DIMS),"T",
            &lower.query(),&point_at(&sheet,hi)?,ClaimSelection::Pinned(heat_capacity_claim))?;
        let density=sheet.areal_density_kg_m2()/sheet.thickness_m();
        let mut samples=Vec::with_capacity(line.knots.len());let mut knots=Vec::with_capacity(line.knots.len());
        for (temperature,_) in line.knots {
            let material=sample(&sheet,card,heat_capacity_claim,CP,SPECIFIC_HEAT_CAPACITY_DIMS,temperature)?;
            let cp=material.property(CP).ok_or(SheetThermalError::Invalid("missing resolved Cp"))?.value_si();
            knots.push(HeatCapacityKnot {temperature_k:temperature,specific_heat_capacity_j_kg_k:cp,
                bulk_density_kg_m3:density});
            samples.push(material);
        }
        let curve=EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity(card.content_hash(),
            SolidLiquidPhase::Solid,reference_h_j_kg,&knots,maximum_chart_error_j_kg)?;
        Ok(Self {sheet,curve,heat_capacity_samples:samples})
    }
    /// Whole original sheet bundle, independent mechanical channels retained.
    pub const fn sheet(&self)->&ResolvedSheet {&self.sheet}
    /// Existing material enthalpy curve, not a wrapper around another integrator.
    pub const fn curve(&self)->&EquilibriumEnthalpyPhaseCurve {&self.curve}
    /// Actual selected source evaluations, including every interior Cp knot.
    pub fn heat_capacity_samples(&self)->&[ResolvedMaterialStatePoint] {&self.heat_capacity_samples}
    /// Same single-phase inversion already used by resolved bulk caloric charts.
    pub fn state_at_temperature(&self,t:f64)->Result<EquilibriumPhaseState,SheetThermalError> {
        let h=crate::heat_capacity_state::enthalpy_at_temperature(&self.curve,t)?;
        Ok(self.curve.state_at_specific_enthalpy(h)?)
    }
}

/// Resolve the same specimen's positive conductivity at a declared temperature.
/// This works independently of heat capacity: a conductivity-only specimen may
/// support steady Fourier transport without being falsely transient-ready.
/// Retain the returned common source bundle with the consumer's numeric result.
/// Through-thickness textile conductivity is NOT automatically an isotropic
/// three-dimensional conductivity tensor; direction/geometry remain explicit.
/// Typed pressure, moisture and composition remain at the original held state.
pub fn conductivity_at_temperature(sheet:&ResolvedSheet,card:&MaterialCard,claim:ClaimId,t:f64)
    ->Result<ResolvedMaterialStatePoint,SheetThermalError> {
    let resolved=sample(sheet,card,claim,K,THERMAL_CONDUCTIVITY_DIMS,t)?;
    if !thermal_source_supported(resolved.property(K)
        .ok_or(SheetThermalError::Invalid("missing conductivity"))?) {
        return Err(SheetThermalError::Invalid("conductivity needs supported held conditions, not a frozen foreign-axis law"));
    }
    Ok(resolved)
}
