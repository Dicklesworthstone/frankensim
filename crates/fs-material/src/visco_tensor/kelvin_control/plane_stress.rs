//! Original Kelvin/Maxwell plane stress with resolved transverse viscosity.
//!
//! Write epsilon=P*u+S*v, with P the relaxed equilibrium lift and S the three
//! transverse selectors. In energy coordinates z=[G_F*v;q], the constraint is
//! z'=-A*z+H*u', A=C*D_FF^-1*C^T+diag(0,1/tau)=F^T*F in the definite case.
//! For singular transverse viscosity, original null traction constrains z;
//! native QR restricts the positive factor to that admissible subspace. Native
//! QR/SVD retains every original Maxwell coordinate. Total physical heat is
//! the restricted positive quadratic form plus the
//! positive residual Kelvin Gram applied to u'. The original tensor reaction
//! supplies independent work; no Kelvin/Maxwell heat partition is inferred.

use super::{TensorKelvinStressControlBudget, checked, zeros};
use crate::plane_stress_viscosity::MaterialPlane;
use crate::visco::tensor::{HarmonicEnergySpectrum, HarmonicWorkSum as Sum, StrainParts};
use crate::visco::{MAX_HARMONIC_MATERIAL_LINES, TensorPronyHarmonicSpectrum};
use crate::visco_tensor::{TensorProny, TensorPronyError, TensorPronyHistory, Voigt,
    apply, dot, finite, invalid, poll, tensor};
use crate::visco_tensor::stress_control::TensorPlaneStressHarmonicBudget;
use fs_fft::harmonic_window::{HarmonicDecayIntervalLimits, HarmonicPhaseInterval,
    PreparedHarmonicDecayInterval, ProductSpectrumError, phase_angular_interval_from_clock_dd};
use fs_la::factor::{Qr, qr};
use fs_math::{c64::C64, dd::Dd};
use std::sync::Arc;

mod prepare;
mod harmonic;
mod semidefinite;
mod ramp;
pub use ramp::{TensorKelvinPlaneStressRampStep, TensorKelvinPlaneStressRampTangent};
pub use harmonic::TensorKelvinPlaneStressHarmonicWindow;
type Result<T> = core::result::Result<T, TensorPronyError>;
const LINES: usize = MAX_HARMONIC_MATERIAL_LINES;
const ENGINEERING: [f64;3] = [1.0,1.0,core::f64::consts::FRAC_1_SQRT_2];
type PlaneParts = [[Dd;2];3];

/// Unweighted physical interval account per original reference volume. Material
/// heat is the combined original Kelvin and Maxwell heat, not a branch split.
#[derive(Clone, Debug, PartialEq)]
pub struct TensorKelvinPlaneStressHarmonicAccount {
    pub start_stored_energy_j_m3: f64,
    pub end_stored_energy_j_m3: f64,
    pub stored_energy_change_j_m3: f64,
    pub material_heat_j_m3: f64,
    pub stress_work_j_m3: f64,
    pub energy_balance_residual_j_m3: f64,
}

/// Separate source-bound endpoint candidate; every tensor uses tensor shears.
#[derive(Clone, Debug, PartialEq)]
pub struct TensorKelvinPlaneStressHarmonicPoint {
    pub history: TensorPronyHistory,
    pub account: TensorKelvinPlaneStressHarmonicAccount,
    pub strain_increment: Voigt,
    pub strain_rate: Voigt,
    /// Conservative whole-interval full-tensor bounds, not an outward certificate.
    pub strain_bounds: [Voigt;2],
}

/// Causal plane stress of the unchanged full anisotropic Kelvin/Maxwell law.
/// The original viscosity may be semidefinite. Besides a resolved SPD transverse
/// block, original independent projected Gram columns admit transverse ranks
/// zero, one and two. Algebraic strain constraints retain every original branch;
/// unresolved/redundant projected ranks refuse rather than receive a loss floor.
/// Source parameters and temperature are fixed, with no fitted branch
/// reduction, initial strain jump or replacement of original history memory.
pub struct PreparedTensorKelvinPlaneStress {
    model: Arc<TensorProny>,
    plane: MaterialPlane,
    inside: [usize;3],
    outside: [usize;3],
    equilibrium: Root3,
    viscosity: Option<Root3>,
    constraint: Option<semidefinite::ConstrainedPlaneStress>,
    map: [Voigt;3],
    direct_gram: [[f64;6];6],
    direct_work: [[f64;6];6],
    rates: Vec<f64>,
    rate_roots: Vec<f64>,
    modes: Vec<f64>,
    forcing: Vec<[f64;3]>,
    strain_modes: Vec<Voigt>,
}

struct Root3 { qr: Qr, scale: f64 }
impl Root3 {
    fn new(gram: &[[f64;6];6], axes: [usize;3]) -> Result<Self> {
        let scale=gram.iter().flat_map(|row|axes.map(|i|row[i].abs())).fold(0.0_f64,f64::max);
        if !scale.is_finite() || scale<=0.0 { return Err(invalid("missing positive transverse material factor")); }
        let rows:Vec<_>=gram.iter().flat_map(|row|axes.map(|i|row[i]/scale)).collect();
        let qr=qr(&rows,6,3);
        if (0..3).any(|i|!qr.r(i,i).is_finite() || qr.r(i,i).abs()<=64.0*f64::EPSILON) {
            return Err(invalid("transverse Kelvin/equilibrium factor is not resolved positive definite"));
        }
        Ok(Self {qr,scale})
    }
    fn solve(&self,rhs:[Dd;3],transpose:bool)->Result<[Dd;3]> {
        let mut x=[Dd::ZERO;3];
        for index in 0..3 {
            let i=if transpose {index} else {2-index};
            let mut value=rhs[i]/Dd::from_f64(self.scale);
            for j in 0..3 {
                if (transpose && j<i) || (!transpose && j>i) {
                    value=value-Dd::from_f64(if transpose {self.qr.r(j,i)} else {self.qr.r(i,j)})*x[j];
                }
            }
            x[i]=number(value/Dd::from_f64(self.qr.r(i,i)))?;
        }
        Ok(x)
    }
    fn apply(&self,value:[Dd;3])->Result<[Dd;3]> {
        let mut out=[Dd::ZERO;3];
        for i in 0..3 {
            let mut sum=Sum::default();
            for j in i..3 { for word in [value[j].hi,value[j].lo] {
                sum.add([self.scale,self.qr.r(i,j),word,1.0]);
            } }
            out[i]=number(sum.retained_value())?;
        }
        Ok(out)
    }
}
fn number(value:Dd)->Result<Dd> {
    if value.hi.is_finite() && value.lo.is_finite() {Ok(value)}
    else {Err(invalid("unrepresentable causal Kelvin plane-stress response"))}
}
fn moment_error(error:ProductSpectrumError)->TensorPronyError {
    if error==ProductSpectrumError::Cancelled {TensorPronyError::Cancelled}
    else {invalid("causal Kelvin plane-stress moment is unresolved")}
}
fn add(sum:&mut Sum,value:Dd,factor:f64) {
    sum.add([value.hi,factor,1.0,1.0]); sum.add([value.lo,factor,1.0,1.0]);
}
fn l1(value:[Dd;2])->f64 {value.iter().map(|v|v.hi.abs()+v.lo.abs()).sum()}

impl PreparedTensorKelvinPlaneStress {
    pub const fn plane(&self)->MaterialPlane {self.plane}
    /// Complete original equilibrium and Maxwell coordinates, including those
    /// participating in algebraic constraints; no branch coordinate is removed.
    pub fn state_coordinates(&self)->usize {3+6*self.model.branches().len()}
    /// Original resolved rank of the transverse Kelvin viscosity block.
    pub fn transverse_viscosity_rank(&self)->usize {self.constraint.as_ref().map_or(3,|c|c.rank)}
}
