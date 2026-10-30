//! Reversible thermal potential of the original isotropic elastic law.
//!
//! At fixed total Mandel strain, e=epsilon-integral(alpha(T) dT)*I and
//! Psi=V e:C(T):e/2. The selected E(T), nu(T) and instantaneous alpha(T)
//! are independently affine. Bulk and shear stiffness remain their exact
//! rational functions of those sources, including first and second thermal
//! derivatives. S=-dPsi/dT, U=Psi+T S and C_epsilon=-T d2Psi/dT2.
//! The temperature derivative is at fixed TOTAL strain, not along a solved
//! equilibrium path. The native elastic response and original source paths
//! remain authoritative; this is not a new mechanical or thermal integrator.

use super::IsotropicElastic;
use crate::{MaterialError, SmallStrainLaw};
use crate::field_storage::{FieldStorageError, admit_source_temperature_refresh,
    thermal::affine_source_path};
use crate::state_point::{ResolvedMaterialStatePoint, MaterialStatePointError,
    YOUNG_MODULUS_PROPERTY, POISSON_RATIO_PROPERTY,
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, MaterialPropertySelection,
    integrate_isotropic_thermal_expansion};
use fs_matdb::MaterialCard;
use fs_qty::{Dims, Pressure, QuantitySpec};

#[derive(Debug)]
pub enum ThermoelasticError {
    Source(MaterialStatePointError),
    Material(MaterialError),
    Invalid(&'static str),
}
impl From<MaterialStatePointError> for ThermoelasticError {
    fn from(error: MaterialStatePointError) -> Self { Self::Source(error) }
}
impl From<MaterialError> for ThermoelasticError {
    fn from(error: MaterialError) -> Self { Self::Material(error) }
}
impl core::fmt::Display for ThermoelasticError {
    fn fmt(&self,f:&mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f,"isotropic thermoelastic source: {self:?}")
    }
}
impl std::error::Error for ThermoelasticError {}
fn field(error: FieldStorageError) -> ThermoelasticError {
    match error { FieldStorageError::Material(error)=>ThermoelasticError::Source(error),
        FieldStorageError::Invalid(message)=>ThermoelasticError::Invalid(message) }
}
fn finite(value:f64) -> Result<f64,ThermoelasticError> {
    if value.is_finite() {Ok(value)}
    else {Err(ThermoelasticError::Invalid("thermoelastic arithmetic is unrepresentable"))}
}
fn value(source:&ResolvedMaterialStatePoint,name:&str,dims:Dims) -> Result<f64,ThermoelasticError> {
    let property=source.property(name).ok_or(ThermoelasticError::Invalid("thermoelastic property is missing or ambiguous"))?;
    if property.requirement().quantity()!=QuantitySpec::dimensional(dims) {
        return Err(ThermoelasticError::Invalid("thermoelastic source requires the original dimension-only static coefficient schema"));
    }
    finite(property.value_si())
}
fn temperature(source:&ResolvedMaterialStatePoint) -> Result<f64,ThermoelasticError> {
    source.query_point().iter().find_map(|(name,value)|(name=="T").then_some(*value))
        .ok_or(ThermoelasticError::Invalid("canonical thermoelastic T is missing"))
}

/// One pinned, differentiable, same-card path. Density and inertia remain with
/// the original reference inventory; no current density is synthesized here.
#[derive(Clone,Debug,PartialEq)]
pub struct AffineIsotropicThermoelastic {
    card:MaterialCard,
    elastic_reference:ResolvedMaterialStatePoint,
    expansion_reference:ResolvedMaterialStatePoint,
    temperature_range_k:[f64;2],
    reference_temperature_k:f64,
    youngs_slope_pa_k:f64,
    maximum_youngs_pa:f64,
    poisson:f64,
    poisson_slope_per_k:f64,
    poisson_endpoints:[f64;2],
    expansion_per_k:f64,
    expansion_slope_per_k2:f64,
    maximum_absolute_expansion_per_k:f64,
    maximum_absolute_mandel_strain:f64,
}

/// Same original source response at one admitted physical temperature.
#[derive(Clone,Debug,PartialEq)]
pub struct ResolvedThermoelastic {
    pub elastic:ResolvedMaterialStatePoint,
    pub expansion:ResolvedMaterialStatePoint,
    material:IsotropicElastic,
    temperature_k:f64,
    youngs_slope_pa_k:f64,
    poisson_slope_per_k:f64,
    expansion_per_k:f64,
    expansion_slope_per_k2:f64,
    free_linear_strain:f64,
    maximum_absolute_mandel_strain:f64,
}

/// Integrated elastic-field potentials, separate from the caloric chart.
#[derive(Clone,Copy,Debug,PartialEq)]
pub struct ThermoelasticState {
    pub elastic_strain_mandel:[f64;6],
    pub stress_mandel_pa:[f64;6],
    /// Local d(sigma)/dT at fixed TOTAL Mandel strain [Pa/K].
    /// This mixed derivative is work-conjugate to -d(entropy)/d(strain).
    pub stress_temperature_derivative_mandel_pa_k:[f64;6],
    pub free_energy_j:f64,
    pub entropy_j_k:f64,
    pub internal_energy_j:f64,
    /// Fixed-total-strain contribution, which can have either sign.
    pub heat_capacity_j_k:f64,
}

impl AffineIsotropicThermoelastic {
    /// Admit the complete source interval, retaining every original claim and
    /// fixed external coordinate. The expansion reference is the original
    /// stress-free state; current elasticity may start at another temperature.
    /// E, nu and instantaneous alpha may each be affine over the whole
    /// interval; E stays positive and nu stays in (-1,1/2). Expansion is
    /// integrated by the original pinned source owner. A slope discontinuity
    /// is refused: this smooth potential does not invent a transition entropy.
    pub fn new(card:&MaterialCard,elastic:&ResolvedMaterialStatePoint,
        expansion_reference:&ResolvedMaterialStatePoint,temperature_range_k:[f64;2],
        maximum_absolute_mandel_strain:f64) -> Result<Self,ThermoelasticError> {
        if !maximum_absolute_mandel_strain.is_finite()
            || maximum_absolute_mandel_strain<=0.0 || maximum_absolute_mandel_strain>=1.0 {
            return Err(ThermoelasticError::Invalid("a finite small-strain component bound in (0,1) is required"));
        }
        let reference=elastic.at_query_point(card,&expansion_reference.query())?;
        admit_source_temperature_refresh(elastic,&reference).map_err(field)?;
        let e=value(&reference,YOUNG_MODULUS_PROPERTY,Pressure::DIMS)?;
        let poisson=value(&reference,POISSON_RATIO_PROPERTY,Dims([0;6]))?;
        let alpha=value(expansion_reference,LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
            Dims([0,0,0,-1,0,0]))?;
        IsotropicElastic::new(e,poisson,maximum_absolute_mandel_strain)?;
        let (reference_temperature_k,youngs_slope_pa_k,ends)=affine_source_path(
            card,&reference,YOUNG_MODULUS_PROPERTY,"T",temperature_range_k,true).map_err(field)?;
        let (_,poisson_slope_per_k,poisson_endpoints)=affine_source_path(card,&reference,POISSON_RATIO_PROPERTY,
            "T",temperature_range_k,false).map_err(field)?;
        let (_,expansion_slope_per_k2,alpha_endpoints)=affine_source_path(card,expansion_reference,
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,"T",temperature_range_k,false).map_err(field)?;
        for (slope,endpoints) in [(youngs_slope_pa_k,ends),
            (poisson_slope_per_k,poisson_endpoints),(expansion_slope_per_k2,alpha_endpoints)] {
            if slope==0.0 && endpoints[0]!=endpoints[1] {
                return Err(ThermoelasticError::Invalid("nonconstant thermoelastic source derivative underflows"));
            }
        }
        let law=Self {card:card.clone(),elastic_reference:reference,
            expansion_reference:expansion_reference.clone(),temperature_range_k,
            reference_temperature_k,youngs_slope_pa_k,maximum_youngs_pa:ends[0].max(ends[1]),
            poisson,poisson_slope_per_k,poisson_endpoints,expansion_per_k:alpha,
            expansion_slope_per_k2,
            maximum_absolute_expansion_per_k:alpha_endpoints[0].abs().max(alpha_endpoints[1].abs()),
            maximum_absolute_mandel_strain};
        for t in temperature_range_k {law.at_temperature(t)?;}
        Ok(law)
    }
    pub fn card(&self)->&MaterialCard {&self.card}
    pub fn temperature_range_k(&self)->[f64;2] {self.temperature_range_k}
    pub fn reference_temperature_k(&self)->f64 {self.reference_temperature_k}
    pub fn maximum_absolute_mandel_strain(&self)->f64 {self.maximum_absolute_mandel_strain}

    pub fn at_temperature(&self,t:f64)->Result<ResolvedThermoelastic,ThermoelasticError> {
        if !t.is_finite() || !(self.temperature_range_k[0]..=self.temperature_range_k[1]).contains(&t) {
            return Err(ThermoelasticError::Invalid("temperature leaves the admitted thermoelastic source interval"));
        }
        let point=self.elastic_reference.query().with("T",t)
            .map_err(|_|ThermoelasticError::Invalid("thermoelastic temperature query is unrepresentable"))?;
        let elastic=self.elastic_reference.at_query_point(&self.card,&point)?;
        let expansion=self.expansion_reference.at_query_point(&self.card,&point)?;
        let e=value(&elastic,YOUNG_MODULUS_PROPERTY,Pressure::DIMS)?;
        let poisson=value(&elastic,POISSON_RATIO_PROPERTY,Dims([0;6]))?;
        let alpha=value(&expansion,LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
            Dims([0,0,0,-1,0,0]))?;
        let delta=t-self.reference_temperature_k;
        for (actual,reference,slope,scale) in [
            (e,value(&self.elastic_reference,YOUNG_MODULUS_PROPERTY,Pressure::DIMS)?,
                self.youngs_slope_pa_k,self.maximum_youngs_pa),
            (poisson,self.poisson,self.poisson_slope_per_k,
                self.poisson_endpoints[0].abs().max(self.poisson_endpoints[1].abs())),
            (alpha,self.expansion_per_k,self.expansion_slope_per_k2,
                self.maximum_absolute_expansion_per_k),
        ] {
            let expected=finite(slope.mul_add(delta,reference))?;
            if (actual-expected).abs()>512.0*f64::EPSILON*scale.max(actual.abs()).max(expected.abs()) {
                return Err(ThermoelasticError::Invalid("thermoelastic sources left their admitted affine paths"));
            }
        }
        let selected=self.expansion_reference.property(LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY)
            .ok_or(ThermoelasticError::Invalid("original expansion claim is missing"))?
            .answer().receipt.selected;
        let integrated=integrate_isotropic_thermal_expansion(&self.card,
            &self.expansion_reference.query(),&point,MaterialPropertySelection::PinnedByProperty(
                vec![(LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY.to_owned(),selected)]))?;
        let material=IsotropicElastic::new(e,poisson,self.maximum_absolute_mandel_strain)?;
        Ok(ResolvedThermoelastic {elastic,expansion,material,temperature_k:t,
            youngs_slope_pa_k:self.youngs_slope_pa_k,poisson_slope_per_k:self.poisson_slope_per_k,
            expansion_per_k:alpha,expansion_slope_per_k2:self.expansion_slope_per_k2,
            free_linear_strain:integrated.free_linear_strain(),
            maximum_absolute_mandel_strain:self.maximum_absolute_mandel_strain})
    }

    /// Match the actual native source fields; a different card, selection,
    /// external condition or same-name coefficient cannot replace live history.
    pub fn at_materials(&self,elastic:&ResolvedMaterialStatePoint,expansion:&ResolvedMaterialStatePoint)
        -> Result<ResolvedThermoelastic,ThermoelasticError> {
        admit_source_temperature_refresh(&self.elastic_reference,elastic).map_err(field)?;
        admit_source_temperature_refresh(&self.expansion_reference,expansion).map_err(field)?;
        let t=temperature(elastic)?;
        if t.to_bits()!=temperature(expansion)?.to_bits() || elastic.query()!=expansion.query() {
            return Err(ThermoelasticError::Invalid("elasticity and expansion must share the actual physical query"));
        }
        let at=self.at_temperature(t)?;
        for (expected,actual) in [(&at.elastic,elastic),(&at.expansion,expansion)] {
            if expected.properties().iter().zip(actual.properties())
                .any(|(a,b)|a.value_si().to_bits()!=b.value_si().to_bits()) {
                return Err(ThermoelasticError::Invalid("native elasticity and caloric potential selected different source values"));
            }
        }
        Ok(at)
    }

    /// Conservative lower bound over the full temperature and TOTAL Mandel
    /// component envelope. Add this to the actual V*rho0*dh/dT before coupling.
    pub fn minimum_heat_capacity_j_k(&self,volume_m3:f64)->Result<f64,ThermoelasticError> {
        if !volume_m3.is_finite() || volume_m3<=0.0 {
            return Err(ThermoelasticError::Invalid("positive finite reference volume is required"));
        }
        // Positive arithmetic is rounded outward. Denominator minima are
        // rounded downward before reciprocation; an unrepresentable bound
        // refuses instead of admitting an unproved positive total capacity.
        let nu_low=self.poisson_endpoints[0].min(self.poisson_endpoints[1]);
        let nu_high=self.poisson_endpoints[0].max(self.poisson_endpoints[1]);
        let d=(1.0-2.0*nu_high).next_down();
        let q=(1.0+nu_low).next_down();
        if d<=0.0 || q<=0.0 {
            return Err(ThermoelasticError::Invalid("thermoelastic capacity denominator is unrepresentable"));
        }
        let inv_d=upper_div(1.0,d)?;
        let inv_q=upper_div(1.0,q)?;
        let e=self.maximum_youngs_pa.next_up();
        let et=self.youngs_slope_pa_k.abs();
        let nut=self.poisson_slope_per_k.abs();
        let bulk=upper_div(upper_mul(e,inv_d)?,3.0)?;
        let bulk_first=upper_add(upper_div(upper_mul(et,inv_d)?,3.0)?,
            upper_mul(2.0,upper_mul(nut,upper_mul(bulk,inv_d)?)?)?)?;
        let bulk_second=upper_mul(4.0,upper_mul(nut,upper_mul(bulk_first,inv_d)?)?)?;
        let shear=upper_div(upper_mul(e,inv_q)?,2.0)?;
        let shear_first=upper_add(upper_div(upper_mul(et,inv_q)?,2.0)?,
            upper_mul(nut,upper_mul(shear,inv_q)?)?)?;
        let shear_second=upper_mul(2.0,upper_mul(nut,upper_mul(shear_first,inv_q)?)?)?;
        let alpha=self.maximum_absolute_expansion_per_k;
        let alpha_first=self.expansion_slope_per_k2.abs();
        let span=finite((self.temperature_range_k[0]-self.reference_temperature_k).abs()
            .max((self.temperature_range_k[1]-self.reference_temperature_k).abs()).next_up())?;
        let strain=self.maximum_absolute_mandel_strain;
        let trace=upper_mul(3.0,upper_add(strain,upper_mul(alpha,span)?)?)?;
        // The limit is on SIX MANDEL components, whose squared norm is the
        // tensor Frobenius norm. Isotropic free strain does not change dev(e).
        let deviator_squared=upper_mul(6.0,upper_mul(strain,strain)?)?;
        let linear=upper_mul(upper_add(upper_mul(6.0,upper_mul(bulk_first,alpha)?)?,
            upper_mul(3.0,upper_mul(bulk,alpha_first)?)?)?,trace)?;
        let expansion=upper_mul(9.0,upper_mul(bulk,upper_mul(alpha,alpha)?)?)?;
        let volumetric=upper_div(upper_mul(bulk_second,upper_mul(trace,trace)?)?,2.0)?;
        let distortional=upper_mul(shear_second,deviator_squared)?;
        let magnitude=upper_mul(upper_mul(self.temperature_range_k[1],volume_m3)?,
            upper_add(upper_add(linear,expansion)?,upper_add(volumetric,distortional)?)?)?;
        // nu(T) produces stiffness curvature even when thermal expansion is zero.
        Ok(-magnitude)
    }
}

// Bounds use only nonnegative operands. Preserve exact zero while rounding
// positive underflow upward; infinity is an explicit admission failure.
fn upper_add(a:f64,b:f64)->Result<f64,ThermoelasticError> {
    if a==0.0 && b==0.0 {Ok(0.0)} else {finite((a+b).next_up())}
}
fn upper_mul(a:f64,b:f64)->Result<f64,ThermoelasticError> {
    if a==0.0 || b==0.0 {Ok(0.0)} else {finite((a*b).next_up())}
}
fn upper_div(a:f64,b:f64)->Result<f64,ThermoelasticError> {
    if a==0.0 {Ok(0.0)} else {finite((a/b).next_up())}
}

impl ResolvedThermoelastic {
    pub fn temperature_k(&self)->f64 {self.temperature_k}
    pub fn free_strain(&self)->Result<f64,ThermoelasticError> {
        Ok(self.free_linear_strain)
    }
    pub fn state(&self,total_strain_mandel:[f64;6],volume_m3:f64)
        -> Result<ThermoelasticState,ThermoelasticError> {
        if !volume_m3.is_finite() || volume_m3<=0.0
            || total_strain_mandel.iter().any(|v|!v.is_finite()
                ||v.abs()>self.maximum_absolute_mandel_strain) {
            return Err(ThermoelasticError::Invalid("strain/volume leaves the admitted thermoelastic envelope"));
        }
        let free=self.free_strain()?;
        let mut elastic_strain_mandel=total_strain_mandel;
        for e in &mut elastic_strain_mandel[..3] {*e=finite(*e-free)?;}
        if elastic_strain_mandel.iter().any(|e|e.abs()>self.maximum_absolute_mandel_strain) {
            return Err(ThermoelasticError::Invalid("elastic strain leaves the admitted small-strain envelope"));
        }
        // The existing material law uses tensor-shear Voigt coordinates.
        // Convert both strain and stress explicitly; the work dot product is Mandel.
        let root_two=fs_math::det::sqrt(2.0);
        let mut voigt=elastic_strain_mandel;
        for e in &mut voigt[3..] {*e/=root_two;}
        let mut stress_mandel_pa=self.material.stress(&voigt,&());
        for s in &mut stress_mandel_pa[3..] {*s*=root_two;}
        let free_energy_j=finite(0.5*volume_m3*elastic_strain_mandel.iter()
            .zip(stress_mandel_pa).map(|(e,s)|e*s).sum::<f64>())?;
        let trace=finite(elastic_strain_mandel[..3].iter().sum())?;
        let mean=trace/3.0;
        let deviator_squared=finite(elastic_strain_mandel.iter().enumerate()
            .map(|(k,e)|if k<3 {(e-mean).powi(2)}else{e*e}).sum())?;
        let e=self.material.youngs;
        let nu=self.material.poisson;
        let d=1.0-2.0*nu; let q=1.0+nu;
        let bulk=finite(e/(3.0*d))?;
        let bulk_first=finite((self.youngs_slope_pa_k+2.0*e*self.poisson_slope_per_k/d)/(3.0*d))?;
        let shear_first=finite((self.youngs_slope_pa_k-e*self.poisson_slope_per_k/q)/(2.0*q))?;
        if bulk_first==0.0 && shear_first==0.0
            && (self.youngs_slope_pa_k!=0.0 || self.poisson_slope_per_k!=0.0) {
            return Err(ThermoelasticError::Invalid("nonconstant thermoelastic stiffness derivatives underflow"));
        }
        let bulk_second=finite(4.0*self.poisson_slope_per_k*bulk_first/d)?;
        let shear_second=finite(-2.0*self.poisson_slope_per_k*shear_first/q)?;
        let alpha=self.expansion_per_k;
        let stress_temperature_derivative_mandel_pa_k=core::array::from_fn(|k| {
            if k<3 {
                bulk_first*trace-3.0*bulk*alpha
                    +2.0*shear_first*(elastic_strain_mandel[k]-mean)
            } else {2.0*shear_first*elastic_strain_mandel[k]}
        });
        let entropy_j_k=finite(volume_m3*(3.0*bulk*alpha*trace
            -0.5*bulk_first*trace*trace-shear_first*deviator_squared))?;
        let internal_energy_j=finite(self.temperature_k.mul_add(entropy_j_k,free_energy_j))?;
        let heat_capacity_j_k=finite(self.temperature_k*volume_m3
            *((6.0*bulk_first*alpha+3.0*bulk*self.expansion_slope_per_k2)*trace
                -9.0*bulk*alpha*alpha-0.5*bulk_second*trace*trace
                -shear_second*deviator_squared))?;
        if stress_mandel_pa.iter().chain(stress_temperature_derivative_mandel_pa_k.iter())
            .any(|v|!v.is_finite()) || free_energy_j<0.0 {
            return Err(ThermoelasticError::Invalid("native thermoelastic response is unrepresentable"));
        }
        Ok(ThermoelasticState {elastic_strain_mandel,stress_mandel_pa,
            stress_temperature_derivative_mandel_pa_k,free_energy_j,
            entropy_j_k,internal_energy_j,heat_capacity_j_k})
    }
}

/// Source-derived two-temperature Hertz compliance and thermal derivatives.
pub mod contact;

/// Isentropic incremental stiffness from this same source thermal potential.
pub mod acoustic;
