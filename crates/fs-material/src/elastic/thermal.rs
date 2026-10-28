//! Reversible thermal potential of the original isotropic elastic law.
//!
//! At fixed total Mandel strain, e=epsilon-alpha*(T-Tr)*I and
//! Psi=V e:C(T):e/2. For affine E(T), fixed nu and constant instantaneous
//! expansion alpha, S=V alpha tr(sigma)-(E'/E)Psi and U=Psi+T S.
//! C_epsilon=T V [2(E'/E)alpha tr(sigma)-3 E alpha^2/(1-2nu)].
//! The temperature derivative is at fixed TOTAL strain, not along a solved
//! equilibrium path. The native elastic response and original source paths
//! remain authoritative; this is not a new mechanical or thermal integrator.

use super::IsotropicElastic;
use crate::{MaterialError, SmallStrainLaw};
use crate::field_storage::{FieldStorageError, admit_source_temperature_refresh,
    thermal::affine_source_path};
use crate::state_point::{ResolvedMaterialStatePoint, MaterialStatePointError,
    YOUNG_MODULUS_PROPERTY, POISSON_RATIO_PROPERTY,
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY};
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
    expansion_per_k:f64,
    maximum_absolute_mandel_strain:f64,
}

/// Same original source response at one admitted physical temperature.
#[derive(Clone,Debug,PartialEq)]
pub struct ResolvedThermoelastic {
    pub elastic:ResolvedMaterialStatePoint,
    pub expansion:ResolvedMaterialStatePoint,
    material:IsotropicElastic,
    temperature_k:f64,
    reference_temperature_k:f64,
    youngs_slope_pa_k:f64,
    expansion_per_k:f64,
    maximum_absolute_mandel_strain:f64,
}

/// Integrated elastic-field potentials, separate from the caloric chart.
#[derive(Clone,Copy,Debug,PartialEq)]
pub struct ThermoelasticState {
    pub elastic_strain_mandel:[f64;6],
    pub stress_mandel_pa:[f64;6],
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
    /// E must be positive affine, nu fixed in (-1,1/2), and alpha constant.
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
        let (_,nu_slope,_)=affine_source_path(card,&reference,POISSON_RATIO_PROPERTY,
            "T",temperature_range_k,false).map_err(field)?;
        let (_,alpha_slope,_)=affine_source_path(card,expansion_reference,
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,"T",temperature_range_k,false).map_err(field)?;
        if nu_slope!=0.0 || alpha_slope!=0.0 {
            return Err(ThermoelasticError::Invalid("this potential requires constant Poisson ratio and instantaneous thermal expansion"));
        }
        let law=Self {card:card.clone(),elastic_reference:reference,
            expansion_reference:expansion_reference.clone(),temperature_range_k,
            reference_temperature_k,youngs_slope_pa_k,maximum_youngs_pa:ends[0].max(ends[1]),
            poisson,expansion_per_k:alpha,maximum_absolute_mandel_strain};
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
        let expected=finite(self.youngs_slope_pa_k.mul_add(t-self.reference_temperature_k,
            value(&self.elastic_reference,YOUNG_MODULUS_PROPERTY,Pressure::DIMS)?))?;
        if (e-expected).abs()>512.0*f64::EPSILON*e.abs().max(expected.abs())
            || value(&elastic,POISSON_RATIO_PROPERTY,Dims([0;6]))?.to_bits()!=self.poisson.to_bits()
            || value(&expansion,LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
                Dims([0,0,0,-1,0,0]))?.to_bits()!=self.expansion_per_k.to_bits() {
            return Err(ThermoelasticError::Invalid("thermoelastic sources left their admitted affine/fixed paths"));
        }
        let material=IsotropicElastic::new(e,self.poisson,self.maximum_absolute_mandel_strain)?;
        Ok(ResolvedThermoelastic {elastic,expansion,material,temperature_k:t,
            reference_temperature_k:self.reference_temperature_k,youngs_slope_pa_k:self.youngs_slope_pa_k,
            expansion_per_k:self.expansion_per_k,maximum_absolute_mandel_strain:self.maximum_absolute_mandel_strain})
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
        if self.expansion_per_k==0.0 {return Ok(0.0);}
        let span=(self.temperature_range_k[0]-self.reference_temperature_k).abs()
            .max((self.temperature_range_k[1]-self.reference_temperature_k).abs());
        let elastic_trace_bound=finite(3.0*(self.maximum_absolute_mandel_strain
            +self.expansion_per_k.abs()*span))?;
        let magnitude=finite(self.temperature_range_k[1]*volume_m3/(1.0-2.0*self.poisson)
            *(2.0*self.youngs_slope_pa_k.abs()*self.expansion_per_k.abs()*elastic_trace_bound
                +3.0*self.maximum_youngs_pa*self.expansion_per_k*self.expansion_per_k))?;
        if magnitude<=0.0 {
            return Err(ThermoelasticError::Invalid("nonzero thermoelastic capacity bound is unrepresentable"));
        }
        finite(-magnitude*(1.0+4096.0*f64::EPSILON))
    }
}

impl ResolvedThermoelastic {
    pub fn temperature_k(&self)->f64 {self.temperature_k}
    pub fn free_strain(&self)->Result<f64,ThermoelasticError> {
        finite(self.expansion_per_k*(self.temperature_k-self.reference_temperature_k))
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
        let trace=finite(stress_mandel_pa[..3].iter().sum())?;
        let ratio=finite(self.youngs_slope_pa_k/self.material.youngs)?;
        if self.youngs_slope_pa_k!=0.0 && ratio==0.0 {
            return Err(ThermoelasticError::Invalid("thermoelastic derivative ratio underflows"));
        }
        let entropy_j_k=finite(volume_m3*self.expansion_per_k*trace-ratio*free_energy_j)?;
        let internal_energy_j=finite(self.temperature_k.mul_add(entropy_j_k,free_energy_j))?;
        let heat_capacity_j_k=finite(self.temperature_k*volume_m3
            *(2.0*ratio*self.expansion_per_k*trace
                -3.0*self.material.youngs*self.expansion_per_k*self.expansion_per_k
                    /(1.0-2.0*self.material.poisson)))?;
        if stress_mandel_pa.iter().any(|v|!v.is_finite()) || free_energy_j<0.0 {
            return Err(ThermoelasticError::Invalid("native thermoelastic response is unrepresentable"));
        }
        Ok(ThermoelasticState {elastic_strain_mandel,stress_mandel_pa,free_energy_j,
            entropy_j_k,internal_energy_j,heat_capacity_j_k})
    }
}

/// Source-derived two-temperature Hertz compliance and thermal derivatives.
pub mod contact;
