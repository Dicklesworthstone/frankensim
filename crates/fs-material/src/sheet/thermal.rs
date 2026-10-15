//! The same sheet specimen in the existing bounded enthalpy/transport owners.
//! Apparent reference density is ONLY mu/t, not a sourced solid-fiber density.
use super::*;
use crate::phase::{EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState,
    HeatCapacityKnot, PhaseStateError, SolidLiquidPhase};
use crate::state_point::{SPECIFIC_HEAT_CAPACITY_PROPERTY as CP,
    SPECIFIC_HEAT_CAPACITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY as K,
    THERMAL_CONDUCTIVITY_DIMS};
use fs_matdb::{ClaimId, InterpolationPolicy, MatDbError, PropertyValue};

/// A source query, geometry conversion or caloric-chart refusal.
#[derive(Debug)]
pub enum SheetThermalError {
    /// Original common source/constitutive refusal.
    Sheet(SheetError),
    /// Original query-coordinate admission.
    Query(MatDbError),
    /// Original bounded enthalpy-chart construction/evaluation.
    Phase(PhaseStateError),
    /// A missing, ambiguous or unsupported physical representation.
    Invalid(&'static str),
}
impl From<SheetError> for SheetThermalError {fn from(e:SheetError)->Self{Self::Sheet(e)}}
impl From<MaterialStatePointError> for SheetThermalError {fn from(e:MaterialStatePointError)->Self{Self::Sheet(e.into())}}
impl From<MatDbError> for SheetThermalError {fn from(e:MatDbError)->Self{Self::Query(e)}}
impl From<PhaseStateError> for SheetThermalError {fn from(e:PhaseStateError)->Self{Self::Phase(e)}}
impl core::fmt::Display for SheetThermalError {
    fn fmt(&self,f:&mut core::fmt::Formatter<'_>)->core::fmt::Result {write!(f,"sheet thermal: {self:?}")}
}
impl std::error::Error for SheetThermalError {}

/// Restore ALL original conditioning coordinates and their quantity meanings;
/// change only T. No conditions or replacement material are guessed.
fn point_at(sheet:&ResolvedSheet,temperature:f64)->Result<QueryPoint,SheetThermalError> {
    if !temperature.is_finite() || temperature<=0.0 {return Err(SheetThermalError::Invalid("positive absolute temperature required"));}
    let material=sheet.material();let mut point=QueryPoint::new();
    for (name,value) in material.query_point() {
        let value=*value;
        point=if let Some(&quantity)=material.axis_quantities().get(name) {
            point.with_quantity(name.clone(),quantity,value)?
        } else {point.with(name.clone(),value)?};
    }
    Ok(point.with("T",temperature)?)
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
    /// constructor. Include all source breakpoints in the selected Kelvin span;
    /// tabulated-only knots are not permission to interpolate or evolve between
    /// them. Other conditioning coordinates and typed descriptors remain fixed.
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
        let source=card.claims().claim(heat_capacity_claim)
            .ok_or(SheetThermalError::Invalid("missing selected sheet heat-capacity source"))?;
        let mut temperatures=vec![lo];
        match (&source.value,source.interpolation) {
            (PropertyValue::Scalar {..},InterpolationPolicy::ConstantWithinValidity)=>{},
            (PropertyValue::Curve {abscissa,abscissa_dims,knots,..},InterpolationPolicy::LinearInside)
                if abscissa=="T" && *abscissa_dims==Dims([0,0,0,1,0,0]) && knots.len()<=4096=> {
                temperatures.extend(knots.iter().map(|p|p.0).filter(|&t|lo<t && t<hi));
            },
            _=>return Err(SheetThermalError::Invalid("Cp needs an admitted constant or bounded Kelvin-linear law, not exact-only samples or a different axis")),
        }
        temperatures.push(hi);
        let density=sheet.areal_density_kg_m2()/sheet.thickness_m();
        let mut samples=Vec::with_capacity(temperatures.len());let mut knots=Vec::with_capacity(temperatures.len());
        for temperature in temperatures {
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
pub fn conductivity_at_temperature(sheet:&ResolvedSheet,card:&MaterialCard,claim:ClaimId,t:f64)
    ->Result<ResolvedMaterialStatePoint,SheetThermalError> {
    let resolved=sample(sheet,card,claim,K,THERMAL_CONDUCTIVITY_DIMS,t)?;
    if !resolved.property(K).ok_or(SheetThermalError::Invalid("missing conductivity"))?
        .is_constant_at_fixed_temperature() {
        return Err(SheetThermalError::Invalid("conductivity must be constant or temperature-dependent, not a frozen foreign-axis law"));
    }
    Ok(resolved)
}
