//! Physical forcing and ideal reactions about an admitted affine equilibrium.
//!
//! Peak phasors use exp(-i*omega*t). The original native allowed-direction
//! M/C/K pencil is solved by fs-modal. Original joint forces are projected by
//! N^T and motion is lifted by N. Full physical rows then determine ideal
//! reaction increments, whose original A^T map and actual virtual work are
//! checked independently. Neither a multiplier mass nor penalty is introduced.

use super::{ConstrainedQuadraticModel, ConstrainedLinearizationError, product, represented};
use crate::articulated::motion::reconfigure::constraints::{
    ConstraintError, apply_rows, apply_transpose_rows,
};
use fs_la::factor::{Qr,qr};
use fs_math::{c64::C64,dd::Dd};
use fs_modal::quadratic_modes::{PoleOptions,QuadraticModeError};
use fs_modal::quadratic_modes::forced::{
    MechanicalHarmonicBalance,PreparedQuadraticForcedResponse,QuadraticForcedError,
    QuadraticForcedOptions,evaluate_harmonic_balance,
};

#[derive(Clone,Copy,Debug)]
pub struct ConstrainedForcedControls {
    /// Absolute allowed numerical mean power of the ACTUAL original ideal
    /// reactions [W]. Signed work remains visible; nothing is clipped to zero.
    pub maximum_reaction_power_w:f64,
}

#[derive(Clone,Debug,PartialEq)]
pub enum ConstrainedForcedError {
    Original(ConstrainedLinearizationError),
    Forced(QuadraticForcedError),
    Constraint(ConstraintError),
    Invalid(&'static str),
    Budget {quantity:&'static str,requested:u128,maximum:u128},
    Manifold {row:usize,position_amplitude:f64,velocity_amplitude:f64},
    ReactionWork {power_w:f64,maximum_w:f64},
    Cancelled,
}
impl From<ConstrainedLinearizationError> for ConstrainedForcedError {
    fn from(e:ConstrainedLinearizationError)->Self {Self::Original(e)}
}
impl From<ConstraintError> for ConstrainedForcedError {
    fn from(e:ConstraintError)->Self {
        if e==ConstraintError::Cancelled {Self::Cancelled} else {Self::Constraint(e)}
    }
}
impl From<QuadraticForcedError> for ConstrainedForcedError {
    fn from(e:QuadraticForcedError)->Self {
        if e==QuadraticForcedError::Cancelled {Self::Cancelled} else {Self::Forced(e)}
    }
}
impl From<QuadraticModeError> for ConstrainedForcedError {
    fn from(e:QuadraticModeError)->Self {QuadraticForcedError::from(e).into()}
}
impl core::fmt::Display for ConstrainedForcedError {
    fn fmt(&self,f:&mut core::fmt::Formatter<'_>)->core::fmt::Result {
        write!(f,"constrained physical forced response: {self:?}")
    }
}
impl std::error::Error for ConstrainedForcedError {}
type E=ConstrainedForcedError;
fn poll(cancelled:&mut impl FnMut()->bool)->Result<(),E> {
    if cancelled(){Err(E::Cancelled)}else{Ok(())}
}

/// Local perturbations only: add their real time-dependent values to the
/// separately retained equilibrium. Reactions are increments, not a replacement
/// for the mean ideal reactions that balanced the original loaded root.
#[derive(Clone,Debug,PartialEq)]
pub struct ConstrainedForcedResponse {
    pub reduced_displacement:Vec<C64>,
    pub displacement:Vec<C64>,
    pub velocity:Vec<C64>,
    /// Constraint multipliers in original row order, units and scaling.
    pub multipliers:Vec<C64>,
    /// ACTUAL native A^T*lambda in physical joint effort coordinates.
    pub reaction_effort:Vec<C64>,
    pub position_residuals:Vec<C64>,
    pub velocity_residuals:Vec<C64>,
    pub reduced_relative_force_residual:f64,
    pub balance:MechanicalHarmonicBalance,
}

/// Reusable reduced dynamic factor and original row-space reaction factor.
/// No original model/load query is needed for subsequent force phasors.
pub struct PreparedConstrainedForcedResponse<'a> {
    source:&'a ConstrainedQuadraticModel,
    original:Option<PreparedQuadraticForcedResponse<'a>>,
    reaction_factor:Option<Qr>,
    row_scales:Vec<f64>,
    omega:f64,
    maximum_relative_residual:f64,
    controls:ConstrainedForcedControls,
}

impl ConstrainedQuadraticModel {
    /// Prepare physical forcing at a nonnegative angular frequency, including
    /// DC. The native general quadratic owner retains signed K and nonsymmetric
    /// C/K. Its pole-preparation budget independently bounds source admission.
    ///
    /// Forced options additionally admit the FULL physical reaction problem:
    /// 192*n^2+1024*n bytes and 64*n^3+512*n^2+1024*n work units, excluding
    /// source/caller arrays and the separately admitted quadratic preparation.
    /// A fully locked source has no dynamic factor/probes but does recover its
    /// physical support forces. A genuine reduced resonance still refuses.
    pub fn prepare_forced_response(&self,omega:f64,pole_options:PoleOptions,
        forced_options:QuadraticForcedOptions,controls:ConstrainedForcedControls,
        mut cancelled:impl FnMut()->bool)
        ->Result<PreparedConstrainedForcedResponse<'_>,E> {
        poll(&mut cancelled)?;
        if !omega.is_finite() || omega<0.0 || !controls.maximum_reaction_power_w.is_finite()
            || controls.maximum_reaction_power_w<0.0 {
            return Err(E::Invalid("finite frequency and explicit ideal-reaction work allowance required"));
        }
        let n=self.subspace.physical_dimension();
        let size=n as u128;
        for (quantity,requested,maximum) in [
            ("physical forced payload bytes",192*size*size+1024*size,
                forced_options.maximum_payload_bytes as u128),
            ("physical forced work",64*size*size*size+512*size*size+1024*size,
                forced_options.maximum_work as u128),
        ] {
            if requested>maximum {return Err(E::Budget {quantity,requested,maximum});}
        }
        // Shared residual-policy admission also covers empty/fully locked sources.
        evaluate_harmonic_balance(&[],&[],&[],&[],&[],&[],
            forced_options.maximum_relative_residual,&mut cancelled)?;
        let original=self.prepare_modes(pole_options,&mut cancelled)?
            .map(|p|p.into_forced_response(omega,forced_options,&mut cancelled)).transpose()?;
        let m=self.rows.len();
        let mut row_scales=Vec::with_capacity(m);
        let reaction_factor=if m==0 {None} else {
            let mut transpose=vec![0.0;n*m];
            for (j,row) in self.rows.iter().enumerate() {
                poll(&mut cancelled)?;
                let mut row_scale=0.0_f64;
                for (&a,&s) in row.coefficients.iter().zip(&self.coordinate_scales) {
                    row_scale=row_scale.max(represented(product(a,s)?)?.abs());
                }
                if row_scale==0.0 {return Err(E::Invalid("original reaction row is unrepresentable"));}
                row_scales.push(row_scale);
                for i in 0..n {
                    transpose[i*m+j]=represented(product(row.coefficients[i],
                        self.coordinate_scales[i])?/Dd::from_f64(row_scale))?;
                    if row.coefficients[i]!=0.0 && transpose[i*m+j]==0.0 {
                        return Err(E::Invalid("reaction row normalization erased a coefficient"));
                    }
                }
            }
            poll(&mut cancelled)?;
            let factor=qr(&transpose,n,m);
            poll(&mut cancelled)?;
            for i in 0..m {
                let diagonal=factor.r(i,i).abs();
                if !diagonal.is_finite() || diagonal<self.minimum_scaled_qr_diagonal {
                    return Err(E::Invalid("original reaction row QR is unresolved"));
                }
            }
            Some(factor)
        };
        poll(&mut cancelled)?;
        Ok(PreparedConstrainedForcedResponse {source:self,original,reaction_factor,row_scales,
            omega,maximum_relative_residual:forced_options.maximum_relative_residual,controls})
    }
}

impl PreparedConstrainedForcedResponse<'_> {
    pub fn angular_frequency_rad_s(&self)->f64 {self.omega}
    pub fn physical_dimension(&self)->usize {self.source.subspace.physical_dimension()}
    pub fn reduced_dimension(&self)->usize {self.source.dimension()}

    /// Resolve one original generalized-force phasor without changing the source.
    /// Full native mass and allowed K*N/C*N columns recover physical reactions;
    /// they are not inferred from a reduced equilibrium report or load ratio.
    /// Original mean-position residual plus oscillating row amplitude must fit
    /// each original position tolerance; rate amplitude fits its rate tolerance.
    /// This local solution does not assert finite-motion joint-limit reachability
    /// or decay of the source equilibrium's unstable homogeneous modes.
    pub fn solve(&self,applied:&[C64],mut cancelled:impl FnMut()->bool)
        ->Result<ConstrainedForcedResponse,E> {
        poll(&mut cancelled)?;
        let n=self.physical_dimension();
        let r=self.reduced_dimension();
        if applied.len()!=n || applied.iter().any(|z|!z.re.is_finite()||!z.im.is_finite()) {
            return Err(E::Invalid("complete finite original effort phasor required"));
        }
        let reduced_force=project(self.source,applied,&mut cancelled)?;
        let (reduced_displacement,reduced_rate,reduced_residual)=if let Some(owner)=&self.original {
            let solved=owner.solve(&reduced_force,&mut cancelled)?;
            (solved.displacement,solved.velocity,solved.balance.relative_force_residual)
        } else {(Vec::new(),Vec::new(),0.0)};
        let displacement=self.source.expand_mode(&reduced_displacement,&mut cancelled)?;
        let velocity=derivative(&displacement,self.omega,&mut cancelled)?;
        let acceleration=derivative(&velocity,self.omega,&mut cancelled)?;
        let (inertial,velocity_force,stiffness)=if r==0 {
            (vec![C64::ZERO;n],vec![C64::ZERO;n],vec![C64::ZERO;n])
        } else {
            (action(&self.source.original_mass,n,n,&acceleration,&mut cancelled)?,
             action(&self.source.original_velocity_lift,n,r,&reduced_rate,&mut cancelled)?,
             action(&self.source.original_stiffness_lift,n,r,&reduced_displacement,&mut cancelled)?)
        };
        let mut needed=Vec::with_capacity(n);
        for i in 0..n {
            poll(&mut cancelled)?;
            needed.push(C64::new(
                represented(Dd::from_f64(inertial[i].re)+Dd::from_f64(velocity_force[i].re)
                    +Dd::from_f64(stiffness[i].re)-Dd::from_f64(applied[i].re))?,
                represented(Dd::from_f64(inertial[i].im)+Dd::from_f64(velocity_force[i].im)
                    +Dd::from_f64(stiffness[i].im)-Dd::from_f64(applied[i].im))?));
        }
        let multipliers=self.reaction_multipliers(&needed,&mut cancelled)?;
        let real:Vec<_>=multipliers.iter().map(|z|z.re).collect();
        let imaginary:Vec<_>=multipliers.iter().map(|z|z.im).collect();
        let real=apply_transpose_rows(&self.source.rows,&real,n,&mut cancelled)?;
        let imaginary=apply_transpose_rows(&self.source.rows,&imaginary,n,&mut cancelled)?;
        let reaction_effort=real.into_iter().zip(imaginary).map(|(re,im)|C64::new(re,im))
            .collect::<Vec<_>>();
        let position_residuals=rows(self.source,&displacement,&mut cancelled)?;
        let velocity_residuals=rows(self.source,&velocity,&mut cancelled)?;
        for (i,row) in self.source.rows.iter().enumerate() {
            poll(&mut cancelled)?;
            let position=position_residuals[i].abs();
            let velocity=velocity_residuals[i].abs();
            let total_position=self.source.mean_position_residuals[i].abs()+position;
            if !total_position.is_finite() || !velocity.is_finite()
                || total_position>row.position_tolerance || velocity>row.velocity_tolerance {
                return Err(E::Manifold {row:i,position_amplitude:position,velocity_amplitude:velocity});
            }
        }
        let balance=evaluate_harmonic_balance(applied,&reaction_effort,&inertial,
            &velocity_force,&stiffness,&velocity,self.maximum_relative_residual,&mut cancelled)?;
        if balance.average_reaction_power_w.abs()>self.controls.maximum_reaction_power_w {
            return Err(E::ReactionWork {power_w:balance.average_reaction_power_w,
                maximum_w:self.controls.maximum_reaction_power_w});
        }
        poll(&mut cancelled)?;
        Ok(ConstrainedForcedResponse {reduced_displacement,displacement,velocity,multipliers,
            reaction_effort,position_residuals,velocity_residuals,
            reduced_relative_force_residual:reduced_residual,balance})
    }

    fn reaction_multipliers(&self,needed:&[C64],cancelled:&mut impl FnMut()->bool)
        ->Result<Vec<C64>,E> {
        let Some(factor)=&self.reaction_factor else {return Ok(Vec::new());};
        let mut out=vec![C64::ZERO;self.source.rows.len()];
        for imaginary in [false,true] {
            poll(cancelled)?;
            let mut right=needed.iter().zip(&self.source.coordinate_scales).map(|(z,&s)|
                represented(product(if imaginary{z.im}else{z.re},s)?))
                .collect::<Result<Vec<_>,ConstrainedLinearizationError>>()?;
            let amplitude=right.iter().map(|v|v.abs()).fold(0.0_f64,f64::max);
            if amplitude==0.0 {continue;}
            for value in &mut right {
                let old=*value;
                *value=represented(Dd::from_f64(old)/Dd::from_f64(amplitude))?;
                if old!=0.0 && *value==0.0 {return Err(E::Invalid("reaction effort normalization underflow"));}
            }
            poll(cancelled)?;
            let solution=factor.solve_ls(&right);
            poll(cancelled)?;
            for (i,&mu) in solution.iter().enumerate() {
                let value=represented(Dd::from_f64(mu)*Dd::from_f64(amplitude)
                    /Dd::from_f64(self.row_scales[i]))?;
                if imaginary {out[i].im=value;} else {out[i].re=value;}
            }
        }
        Ok(out)
    }
}

fn project(source:&ConstrainedQuadraticModel,force:&[C64],cancelled:&mut impl FnMut()->bool)
    ->Result<Vec<C64>,E> {
    let re:Vec<_>=force.iter().map(|z|z.re).collect();
    let im:Vec<_>=force.iter().map(|z|z.im).collect();
    let re=source.project_effort(&re,&mut *cancelled)?;
    let im=source.project_effort(&im,&mut *cancelled)?;
    Ok(re.into_iter().zip(im).map(|(re,im)|C64::new(re,im)).collect())
}
fn derivative(q:&[C64],omega:f64,cancelled:&mut impl FnMut()->bool)->Result<Vec<C64>,E> {
    let mut out=Vec::with_capacity(q.len());
    for z in q {
        poll(cancelled)?;
        out.push(C64::new(represented(product(omega,z.im)?)?,
            represented(product(-omega,z.re)?)?));
    }
    Ok(out)
}
fn action(matrix:&[f64],n:usize,m:usize,q:&[C64],cancelled:&mut impl FnMut()->bool)
    ->Result<Vec<C64>,E> {
    let mut out=Vec::with_capacity(n);
    for i in 0..n {
        poll(cancelled)?;
        let mut re=Dd::ZERO;let mut im=Dd::ZERO;
        for j in 0..m {
            re=re+product(matrix[i*m+j],q[j].re)?;
            im=im+product(matrix[i*m+j],q[j].im)?;
        }
        out.push(C64::new(represented(re)?,represented(im)?));
    }
    Ok(out)
}
fn rows(source:&ConstrainedQuadraticModel,q:&[C64],cancelled:&mut impl FnMut()->bool)
    ->Result<Vec<C64>,E> {
    let re:Vec<_>=q.iter().map(|z|z.re).collect();
    let im:Vec<_>=q.iter().map(|z|z.im).collect();
    let re=apply_rows(&source.rows,&re,&mut *cancelled)?;
    let im=apply_rows(&source.rows,&im,&mut *cancelled)?;
    Ok(re.into_iter().zip(im).map(|(re,im)|C64::new(re,im)).collect())
}
