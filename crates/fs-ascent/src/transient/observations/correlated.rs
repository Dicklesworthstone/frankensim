//! Fixed shared-error factors for correlated transient sensor fitting.
//!
//! C = diag(sigma_i^2) + U U^T, with U supplied in each reading's units and
//! sigma the INDEPENDENT noise scale. Columns describe independent unit-normal
//! shared sources; their loadings can couple different channels and timestamps.
//! No dense observation-by-observation covariance or inverse is constructed.
//!
//! With r_i=(prediction_i-reading_i)/sigma_i and A_ij=U_ij/sigma_i, solve
//! (I+A^T A)z=A^T r and v=r-Az. The GLS cost is (||v||^2+||z||^2)/2 and its
//! prediction derivative is v_i/sigma_i. This avoids subtracting two nearly
//! equal quadratic costs. Supplied noise is fixed, not fitted; its constant
//! log-determinant is omitted. No posterior or physical-validation claim.

use super::{ObservedFamily, ObservedModel, SensorData, SensorFamily, SensorLoss, SensorModel};
use crate::transient::{TransientFamily, TransientModel};
use fs_la::factor::{Cholesky, FactorError, cholesky};
use fs_time::adaptive::adjoint::{AdjointError, OdeVjp, trajectory::{
    RecordedRk45, ReplayBudget, TrajectoryError,
    samples::{SampleObjective, SampledGradient, joint::JointSampleObjective},
}};
use std::sync::Arc;

/// The low-rank solve stays within fs-la's unblocked Cholesky panel; long
/// observation loops poll cancellation, and the bounded dense solve is polled
/// before and after. This is not a cancellation-latency guarantee.
pub const MAX_SHARED_FACTORS: usize = 32;

#[derive(Debug, Clone, Copy)]
pub struct SharedNoiseLimits {
    /// Conservative owned scalar envelope, excluding existing SensorData and
    /// caller inputs: m*k + 2*k*k + 6*m + 6*k. RK45 memory is separate.
    pub max_components: usize,
    /// Conservative scalar-visit envelope for preparation and one score.
    /// Repeated evaluations remain bounded by the enclosing study/campaign.
    pub max_work: usize,
    /// Required relative stationarity residual of the profiled source solve.
    /// Positive, finite and below one; no tolerance clipping or silent fallback.
    pub relative_tolerance: f64,
}
#[derive(Debug, Clone, PartialEq)]
pub enum SharedNoiseError {
    Invalid(&'static str),
    Limit(&'static str),
    Factor(FactorError),
    NonFinite,
    Unresolved,
    Cancelled,
}
impl std::fmt::Display for SharedNoiseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result { write!(f, "shared sensor noise: {self:?}") }
}
impl std::error::Error for SharedNoiseError {}
impl From<SharedNoiseError> for TrajectoryError {
    fn from(error: SharedNoiseError) -> Self {
        match error {
            SharedNoiseError::Cancelled => Self::Step(AdjointError::Cancelled),
            other => Self::Observation(other.to_string()),
        }
    }
}
fn poll<C: FnMut() -> bool>(cancelled: &mut C) -> Result<(), SharedNoiseError> {
    if cancelled() { Err(SharedNoiseError::Cancelled) } else { Ok(()) }
}
fn zeros(n: usize) -> Result<Vec<f64>, SharedNoiseError> {
    let mut out=Vec::new();out.try_reserve_exact(n).map_err(|_|SharedNoiseError::Limit("allocation"))?;
    out.resize(n,0.0);Ok(out)
}

#[derive(Debug, Clone, PartialEq)]
pub struct SharedNoiseScore {
    pub value: f64,
    pub prediction_bar: Vec<f64>,
    /// Minimizer z in the declared unit-normal source coordinates, not a
    /// posterior covariance, a calibrated source estimate, or an extra fit DOF.
    pub source_offsets: Vec<f64>,
    pub relative_stationarity: f64,
}

/// Prepared once, then shared immutably across parameter trials. Covariance
/// rows are bound to the exact SensorData order, not sorted independently.
#[derive(Debug)]
pub struct SharedNoise {
    data: SensorData,
    factors: usize,
    normalized: Vec<f64>,
    factor: Option<Cholesky>,
    tolerance: f64,
}
impl SharedNoise {
    pub fn new<C: FnMut() -> bool>(
        data: SensorData, factors: usize, loadings: &[f64], limits: SharedNoiseLimits, cancelled: &mut C,
    ) -> Result<Self, SharedNoiseError> {
        let m=data.readings().len();
        let mk=m.checked_mul(factors).ok_or(SharedNoiseError::Limit("factor dimensions"))?;
        if factors>MAX_SHARED_FACTORS || loadings.len()!=mk {
            return Err(SharedNoiseError::Invalid("expected m-by-k loadings with k <= 32"));
        }
        let kk=factors*factors;
        let components=mk.checked_add(2*kk).and_then(|v|m.checked_mul(6).and_then(|w|v.checked_add(w)))
            .and_then(|v|v.checked_add(6*factors)).ok_or(SharedNoiseError::Limit("scalar extent"))?;
        let work=mk.checked_mul(factors+17).and_then(|v|v.checked_add(factors*kk+16*kk+16*factors))
            .and_then(|v|m.checked_mul(16).and_then(|w|v.checked_add(w))).ok_or(SharedNoiseError::Limit("work extent"))?;
        if components>limits.max_components || work>limits.max_work {
            return Err(SharedNoiseError::Limit("shared-noise preparation/score envelope"));
        }
        if !limits.relative_tolerance.is_finite() || limits.relative_tolerance<=0.0 || limits.relative_tolerance>=1.0 {
            return Err(SharedNoiseError::Invalid("stationarity tolerance must lie in (0,1)"));
        }
        poll(cancelled)?;
        let mut normalized=zeros(mk)?;let mut active=false;
        for (i,reading) in data.readings().iter().enumerate() {
            poll(cancelled)?;
            if reading.loss()!=SensorLoss::Quadratic {
                return Err(SharedNoiseError::Invalid("correlated GLS requires quadratic readings; Huber whitening is not implied"));
            }
            for j in 0..factors {
                let value=loadings[i*factors+j];let a=value/reading.sigma();
                if !value.is_finite() || !a.is_finite() || (value!=0.0 && a==0.0) { return Err(SharedNoiseError::NonFinite); }
                normalized[i*factors+j]=a;active|=a!=0.0;
            }
        }
        let factor=if !active {None} else {
            let mut gram=zeros(kk)?;
            for j in 0..factors {
                for l in 0..=j {
                    let mut value=if j==l {1.0} else {0.0};
                    for i in 0..m {
                        if i%256==0 {poll(cancelled)?;}
                        value=normalized[i*factors+j].mul_add(normalized[i*factors+l],value);
                    }
                    if !value.is_finite() {return Err(SharedNoiseError::NonFinite);}
                    gram[j*factors+l]=value;gram[l*factors+j]=value;
                }
            }
            poll(cancelled)?;
            let factor=cholesky(&gram,factors).map_err(SharedNoiseError::Factor)?;
            poll(cancelled)?;Some(factor)
        };
        poll(cancelled)?;
        Ok(Self {data,factors,normalized,factor,tolerance:limits.relative_tolerance})
    }
    pub fn data(&self) -> &SensorData { &self.data }
    pub fn factors(&self) -> usize { self.factors }
    pub fn is_diagonal(&self) -> bool { self.factor.is_none() }

    /// One atomic profiled GLS evaluation. Up to two residual corrections use
    /// the retained factor; an unresolved solve refuses instead of treating
    /// correlated readings as independent. Numerical tolerances are not interval
    /// certificates. Parameter-dependent covariance is outside this model.
    pub fn score<C: FnMut() -> bool>(&self, predictions: &[f64], cancelled: &mut C)
        -> Result<SharedNoiseScore, SharedNoiseError>
    {
        let m=self.data.readings().len();let k=self.factors;
        if predictions.len()!=m {return Err(SharedNoiseError::Invalid("prediction count"));}
        poll(cancelled)?;
        let mut residual=zeros(m)?;let mut z=zeros(k)?;
        for (i,(&value,reading)) in predictions.iter().zip(self.data.readings()).enumerate() {
            poll(cancelled)?;
            let difference=value-reading.value();
            residual[i]=if difference.is_finite() {difference/reading.sigma()}
                else {value/reading.sigma()-reading.value()/reading.sigma()};
            if !value.is_finite() || !residual[i].is_finite() {return Err(SharedNoiseError::NonFinite);}
            for j in 0..k {z[j]=self.normalized[i*k+j].mul_add(residual[i],z[j]);}
        }
        if z.iter().any(|v|!v.is_finite()) {return Err(SharedNoiseError::NonFinite);}
        if let Some(factor)=&self.factor {poll(cancelled)?;factor.solve(&mut z);poll(cancelled)?;}
        let mut v=zeros(m)?;let mut correction=zeros(k)?;let mut relative=0.0f64;
        for attempt in 0..=2 {
            if z.iter().any(|v|!v.is_finite()) {return Err(SharedNoiseError::NonFinite);}
            for i in 0..m {
                poll(cancelled)?;v[i]=residual[i];
                for j in 0..k {v[i]=(-self.normalized[i*k+j]).mul_add(z[j],v[i]);}
            }
            relative=0.0;
            for j in 0..k {
                let mut value=-z[j];let mut magnitude=z[j].abs();
                for i in 0..m {
                    if i%256==0 {poll(cancelled)?;}
                    value=self.normalized[i*k+j].mul_add(v[i],value);
                    magnitude+=self.normalized[i*k+j].abs()*v[i].abs();
                }
                if !value.is_finite() || !magnitude.is_finite() {return Err(SharedNoiseError::NonFinite);}
                correction[j]=value;
                let ratio=if magnitude==0.0 {0.0} else {value.abs()/magnitude};
                relative=relative.max(ratio);
            }
            if relative<=self.tolerance {break;}
            if attempt==2 {return Err(SharedNoiseError::Unresolved);}
            let factor=self.factor.as_ref().ok_or(SharedNoiseError::Unresolved)?;
            poll(cancelled)?;factor.solve(&mut correction);poll(cancelled)?;
            for (a,b) in z.iter_mut().zip(&correction) {*a+=b;}
        }
        let mut value=0.0f64;
        for (i,item) in v.iter().chain(&z).enumerate() {
            if i%256==0 {poll(cancelled)?;}value=(0.5*item).mul_add(*item,value);
        }
        for (i,(seed,reading)) in v.iter_mut().zip(self.data.readings()).enumerate() {
            if i%256==0 {poll(cancelled)?;}*seed/=reading.sigma();
        }
        if !value.is_finite() || v.iter().any(|x|!x.is_finite()) {return Err(SharedNoiseError::NonFinite);}
        poll(cancelled)?;
        Ok(SharedNoiseScore {value,prediction_bar:v,source_offsets:z,relative_stationarity:relative})
    }
}

/// Upgrade an existing observed family without changing its dynamics, sensor
/// maps, optimizer, or campaign interface. Sources couple rows WITHIN this
/// family; separate campaign experiments remain separate covariance blocks.
pub struct CorrelatedFamily<F> { observed: ObservedFamily<F>, noise: Arc<SharedNoise> }
impl<F> CorrelatedFamily<F> {
    pub fn new<C: FnMut() -> bool>(observed: ObservedFamily<F>, factors: usize, loadings: &[f64],
        limits: SharedNoiseLimits, cancelled: &mut C) -> Result<Self, SharedNoiseError>
    {
        let noise=Arc::new(SharedNoise::new(observed.data().clone(),factors,loadings,limits,cancelled)?);
        Ok(Self {observed,noise})
    }
    pub fn data(&self) -> &SensorData {self.noise.data()}
    pub fn noise(&self) -> &SharedNoise {&self.noise}
}
pub struct CorrelatedModel<M> { observed: ObservedModel<M>, noise: Arc<SharedNoise> }
impl<F: SensorFamily> TransientFamily for CorrelatedFamily<F> {
    type Model=CorrelatedModel<F::Model>;
    fn bounds(&self)->&[[f64;2]] {self.observed.bounds()}
    fn sample_times(&self)->&[f64] {self.data().times()}
    fn instantiate(&self,point:&[f64])->Result<Self::Model,String> {
        Ok(CorrelatedModel {observed:self.observed.instantiate(point)?,noise:self.noise.clone()})
    }
}
impl<M: SensorModel> OdeVjp for CorrelatedModel<M> {
    fn dimension(&self)->usize {self.observed.dimension()}
    fn parameter_count(&self)->usize {self.observed.parameter_count()}
    fn rhs(&self,t:f64,x:&[f64],out:&mut[f64]) {self.observed.rhs(t,x,out);}
    fn rhs_vjp(&self,t:f64,x:&[f64],b:&[f64],xb:&mut[f64],pb:&mut[f64])->Result<(),String> {
        self.observed.rhs_vjp(t,x,b,xb,pb)
    }
}
impl<M: SensorModel> SampleObjective for CorrelatedModel<M> {
    fn evaluate(&self,i:usize,t:f64,x:&[f64],xb:&mut[f64],pb:&mut[f64])->Result<f64,String> {
        if !self.noise.is_diagonal() {return Err("correlated observations require the joint trajectory pullback".into());}
        self.observed.evaluate(i,t,x,xb,pb)
    }
}
impl<M: SensorModel> TransientModel for CorrelatedModel<M> {
    fn initial_values(&self)->&[f64] {self.observed.initial_values()}
    fn initial_vjp(&self,b:&[f64],pb:&mut[f64])->Result<(),String> {self.observed.initial_vjp(b,pb)}
    fn observation_pullback<C: FnMut() -> bool>(&self,recording:&RecordedRk45<'_,Self>,budget:ReplayBudget,
        cancelled:&mut C)->Result<SampledGradient,TrajectoryError>
    {
        if self.noise.is_diagonal() {recording.pullback_samples(self,budget,cancelled)}
        else {recording.pullback_joint(self,budget,self.data_len(),cancelled)}
    }
}
impl<M> CorrelatedModel<M> {fn data_len(&self)->usize {self.noise.data().readings().len()}}
impl<M: SensorModel> JointSampleObjective for CorrelatedModel<M> {
    fn observe(&self,i:usize,t:f64,x:&[f64])->Result<f64,String> {
        let reading=self.noise.data().readings().get(i).ok_or("unknown correlated sample")?;
        if t!=reading.time() || x.len()!=self.dimension() {return Err("correlated sample binding mismatch".into());}
        self.observed.model.predict(reading.channel(),t,x)
    }
    fn loss<C: FnMut() -> bool>(&self,predictions:&[f64],seeds:&mut[f64],direct:&mut[f64],cancelled:&mut C)
        ->Result<f64,TrajectoryError>
    {
        if seeds.len()!=self.data_len() || direct.len()!=self.parameter_count() {
            return Err(TrajectoryError::Observation("correlated derivative dimension mismatch".into()));
        }
        let result=self.noise.score(predictions,cancelled)?;
        seeds.copy_from_slice(&result.prediction_bar);direct.fill(0.0);Ok(result.value)
    }
    #[allow(clippy::too_many_arguments)]
    fn observation_vjp(&self,i:usize,t:f64,x:&[f64],seed:f64,xb:&mut[f64],pb:&mut[f64])->Result<(),String> {
        let reading=self.noise.data().readings().get(i).ok_or("unknown correlated sample")?;
        if t!=reading.time() || x.len()!=self.dimension() || xb.len()!=self.dimension() || pb.len()!=self.parameter_count() {
            return Err("correlated derivative binding mismatch".into());
        }
        self.observed.model.prediction_vjp(reading.channel(),t,x,seed,xb,pb)
    }
}

#[cfg(test)]
mod tests;
