//! Complete forced response of the original general quadratic pencil.
//!
//! Peak phasors use exp(-i*omega*t): (K-omega^2*M-i*omega*C) q=F.
//! Existing quadratic admission, binary physical scaling and polynomial
//! assembly are retained. Existing cancellable complex LU supplies the solve;
//! no mode truncation, damping floor or spectral reconstruction is introduced.
//! This is a particular response, not a claim that unstable transients decay.

use super::{PreparedQuadraticModes, QuadraticModeError, finite, scale, scale_complex,
    quadratic_fraction, polynomial::ComplexDd};
use fs_la::eigen_complex::{LuComplex, LuComplexError, lu_complex_with_cancel, with_lu_workers};
use fs_math::{c64::C64, dd::Dd};

#[derive(Clone, Copy, Debug)]
pub struct QuadraticForcedOptions {
    /// Componentwise relative balance of the original separate force actions.
    /// This is a backward residual, not a forward-response bound near resonance.
    pub maximum_relative_residual: f64,
    /// Additional factor/solve payload, beyond original quadratic preparation
    /// and caller source arrays: conservatively 96*n^2+512*n bytes.
    pub maximum_payload_bytes: usize,
    /// Additional factor/solve work score: 16*n^3+128*n^2+512*n.
    pub maximum_work: u64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum QuadraticForcedError {
    Original(QuadraticModeError),
    Linear(LuComplexError),
    Invalid(&'static str),
    Residual { measured: f64, maximum: f64 },
    Cancelled,
}
impl From<QuadraticModeError> for QuadraticForcedError {
    fn from(e: QuadraticModeError) -> Self {
        if e == QuadraticModeError::Cancelled { Self::Cancelled } else { Self::Original(e) }
    }
}
impl From<LuComplexError> for QuadraticForcedError {
    fn from(e: LuComplexError) -> Self {
        if e == LuComplexError::Cancelled { Self::Cancelled } else { Self::Linear(e) }
    }
}
impl core::fmt::Display for QuadraticForcedError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "original quadratic forced response: {self:?}")
    }
}
impl std::error::Error for QuadraticForcedError {}
type E = QuadraticForcedError;
fn poll(cancelled: &mut impl FnMut()->bool) -> Result<(), E> {
    if cancelled() { Err(E::Cancelled) } else { Ok(()) }
}
fn exponent(z: C64) -> Option<i32> {
    [z.re,z.im].into_iter().filter(|v| *v!=0.0).map(|v|quadratic_fraction(v).1).max()
}
fn divide_checked(z: C64, a: f64) -> Result<C64,E> {
    let out=C64::new(z.re/a,z.im/a);
    if !finite(out) || (z.re!=0.0 && out.re==0.0) || (z.im!=0.0 && out.im==0.0) {
        return Err(E::Invalid("force/operator row scaling is not representable"));
    }
    Ok(out)
}

/// Complete physical row balance and signed cycle-mean powers.
/// Velocity-operator work may be active. Skew stiffness may supply/absorb work;
/// neither that work nor ideal support work is silently called heat.
#[derive(Clone, Debug, PartialEq)]
pub struct MechanicalHarmonicBalance {
    /// M*qdd+C*qdot+K*q-F-R in original conjugate effort units.
    pub force_residual: Vec<C64>,
    pub relative_force_residual: f64,
    pub average_input_power_w: f64,
    pub average_reaction_power_w: f64,
    pub average_inertial_power_w: f64,
    pub average_velocity_power_w: f64,
    pub average_stiffness_power_w: f64,
    /// Input+reaction-inertial-velocity-stiffness, independently accumulated.
    pub power_balance_residual_w: f64,
}

/// Shared original force/work account for a complete mechanical harmonic state.
/// The caller supplies independently evaluated force actions and ACTUAL reaction
/// amplitudes. Zero-dimensional fully constrained models are legitimate here.
#[allow(clippy::too_many_arguments)]
pub fn evaluate_harmonic_balance(
    applied: &[C64], reaction: &[C64], inertial: &[C64], velocity: &[C64],
    stiffness: &[C64], rate: &[C64], maximum_relative_residual: f64,
    mut cancelled: impl FnMut()->bool,
) -> Result<MechanicalHarmonicBalance,E> {
    poll(&mut cancelled)?;
    let n=applied.len();
    if n>super::MAX_NONLINEAR_POLE_DIMENSION
        || !maximum_relative_residual.is_finite() || maximum_relative_residual<=0.0
        || maximum_relative_residual>=1.0
        || [applied,reaction,inertial,velocity,stiffness,rate].iter()
            .any(|v|v.len()!=n || v.iter().any(|&z|!finite(z))) {
        return Err(E::Invalid("complete bounded physical forces/rates and residual tolerance required"));
    }
    let mut force_residual=Vec::with_capacity(n);
    let mut relative=0.0_f64;
    for i in 0..n {
        poll(&mut cancelled)?;
        let entries=[inertial[i],velocity[i],stiffness[i],applied[i],reaction[i]];
        let e=entries.iter().filter_map(|&z|exponent(z)).max().unwrap_or(0);
        let mut sum=ComplexDd::ZERO;
        let mut envelope=0.0_f64;
        for (index,&z) in entries.iter().enumerate() {
            let z=scale_complex(z,-e)?;
            envelope+=z.abs();
            let z=ComplexDd::from_complex(z);
            sum=if index<3 {sum.add(z)} else {sum.sub(z)};
        }
        let defect=sum.complex();
        let measured=if envelope==0.0 && defect==C64::ZERO {0.0} else {defect.abs()/envelope};
        if !measured.is_finite() {
            return Err(E::Invalid("original harmonic force residual is not representable"));
        }
        relative=relative.max(measured);
        force_residual.push(scale_complex(defect,e)?);
    }
    if relative>maximum_relative_residual {
        return Err(E::Residual {measured:relative,maximum:maximum_relative_residual});
    }
    let input=mean_power(applied,rate,&mut cancelled)?;
    let reaction=mean_power(reaction,rate,&mut cancelled)?;
    let inertia=mean_power(inertial,rate,&mut cancelled)?;
    let velocity=mean_power(velocity,rate,&mut cancelled)?;
    let stiffness=mean_power(stiffness,rate,&mut cancelled)?;
    let balance=(Dd::from_f64(input)+Dd::from_f64(reaction)-Dd::from_f64(inertia)
        -Dd::from_f64(velocity)-Dd::from_f64(stiffness)).to_f64();
    if !balance.is_finite() {return Err(E::Invalid("physical harmonic power balance is not representable"));}
    poll(&mut cancelled)?;
    Ok(MechanicalHarmonicBalance {force_residual,relative_force_residual:relative,
        average_input_power_w:input,average_reaction_power_w:reaction,
        average_inertial_power_w:inertia,average_velocity_power_w:velocity,
        average_stiffness_power_w:stiffness,power_balance_residual_w:balance})
}

fn mean_power(force:&[C64],rate:&[C64],cancelled:&mut impl FnMut()->bool)->Result<f64,E> {
    let ef=force.iter().filter_map(|&z|exponent(z)).max().unwrap_or(0);
    let ev=rate.iter().filter_map(|&z|exponent(z)).max().unwrap_or(0);
    let mut sum=Dd::ZERO;
    for (&f,&v) in force.iter().zip(rate) {
        poll(cancelled)?;
        let f=scale_complex(f,-ef)?;
        let v=scale_complex(v,-ev)?;
        sum=sum+Dd::from_f64(f.re)*Dd::from_f64(v.re)
            +Dd::from_f64(f.im)*Dd::from_f64(v.im);
    }
    Ok(scale((sum*Dd::from_f64(0.5)).to_f64(),ef+ev)?)
}

/// Source-coordinate displacement/velocity plus the independent full force account.
#[derive(Clone, Debug, PartialEq)]
pub struct QuadraticForcedResponse {
    pub displacement: Vec<C64>,
    pub velocity: Vec<C64>,
    pub balance: MechanicalHarmonicBalance,
}

/// One reusable native factor at one angular frequency. This owns original
/// quadratic preparation while continuing to borrow its authoritative M/C/K.
pub struct PreparedQuadraticForcedResponse<'a> {
    original: PreparedQuadraticModes<'a>,
    omega: f64,
    s: C64,
    row_scales: Vec<f64>,
    factor: LuComplex,
    options: QuadraticForcedOptions,
}

impl<'a> PreparedQuadraticModes<'a> {
    /// Consume this original preparation to bind one physical drive frequency.
    /// DC is static compliance; genuine singular static/undamped resonance
    /// refuses through the original LU. Signed K and nonsymmetric C/K survive.
    pub fn into_forced_response(self, angular_frequency_rad_s:f64,
        options:QuadraticForcedOptions, mut cancelled:impl FnMut()->bool)
        ->Result<PreparedQuadraticForcedResponse<'a>,E> {
        poll(&mut cancelled)?;
        if !angular_frequency_rad_s.is_finite() || angular_frequency_rad_s<0.0
            || !options.maximum_relative_residual.is_finite()
            || options.maximum_relative_residual<=0.0 || options.maximum_relative_residual>=1.0 {
            return Err(E::Invalid("nonnegative frequency and explicit force residual tolerance required"));
        }
        let n=self.n as u128;
        super::budget("forced numerical payload bytes",96*n*n+512*n,options.maximum_payload_bytes as u128)?;
        super::budget("forced arithmetic work",16*n*n*n+128*n*n+512*n,options.maximum_work as u128)?;
        let s=scale_complex(C64::new(0.0,-angular_frequency_rad_s),-self.time_exponent)?;
        // The SAME compensated polynomial owner serves poles and forced states.
        let evaluated=self.evaluate(s,&mut cancelled).map_err(QuadraticModeError::from)?;
        let mut operator=evaluated.matrix;
        let mut row_scales=Vec::with_capacity(self.n);
        for i in 0..self.n {
            poll(&mut cancelled)?;
            let row=&mut operator[i*self.n..(i+1)*self.n];
            let magnitude=row.iter().map(|z|z.abs()).fold(0.0_f64,f64::max);
            if !magnitude.is_finite() || magnitude==0.0 {
                return Err(LuComplexError::Singular {index:i}.into());
            }
            for value in row {*value=divide_checked(*value,magnitude)?;}
            row_scales.push(magnitude);
        }
        // One bounded serial native factor; no hidden worker/payload expansion.
        let factor=with_lu_workers(1,||lu_complex_with_cancel(&operator,self.n,&mut cancelled))?;
        poll(&mut cancelled)?;
        Ok(PreparedQuadraticForcedResponse {original:self,omega:angular_frequency_rad_s,
            s,row_scales,factor,options})
    }
}

impl PreparedQuadraticForcedResponse<'_> {
    pub fn angular_frequency_rad_s(&self)->f64 {self.omega}
    pub fn dimension(&self)->usize {self.original.n}

    /// Reuse the original factor for one complete physical force phasor.
    /// Coordinate scaling is dual: Fhat_i=F_i/(D_i*a^2), q_i=y_i/D_i.
    /// A common force amplitude is deferred until physical recovery; no source
    /// component is silently truncated. Refusals leave this factor reusable.
    pub fn solve(&self, applied:&[C64],mut cancelled:impl FnMut()->bool)
        ->Result<QuadraticForcedResponse,E> {
        poll(&mut cancelled)?;
        let n=self.original.n;
        if applied.len()!=n || applied.iter().any(|&z|!finite(z)) {
            return Err(E::Invalid("complete finite physical force phasor required"));
        }
        let time=self.original.time_exponent;
        let coordinates=&self.original.coordinate_exponents;
        let amplitude=applied.iter().enumerate().filter_map(|(i,&z)|
            exponent(z).map(|e|e-coordinates[i]-2*time)).max().unwrap_or(0);
        let mut rhs=Vec::with_capacity(n);
        for i in 0..n {
            poll(&mut cancelled)?;
            let mapped=scale_complex(applied[i],-coordinates[i]-2*time-amplitude)?;
            rhs.push(divide_checked(mapped,self.row_scales[i])?);
        }
        self.factor.solve_with_cancel(&mut rhs,&mut cancelled)?;
        let mut displacement=Vec::with_capacity(n);
        let mut velocity=Vec::with_capacity(n);
        let mut y=Vec::with_capacity(n);
        for i in 0..n {
            poll(&mut cancelled)?;
            let q=scale_complex(rhs[i],amplitude-coordinates[i])?;
            displacement.push(q);
            // Reconstruct from the actual PUBLISHED q, not hidden solve values.
            let actual=scale_complex(q,coordinates[i]-amplitude)?;
            y.push(actual);
            velocity.push(scale_complex(ComplexDd::from_complex(actual)
                .mul(ComplexDd::from_complex(self.s)).complex(),amplitude-coordinates[i]+time)?);
        }
        let mut actions=[vec![C64::ZERO;n],vec![C64::ZERO;n],vec![C64::ZERO;n]];
        for i in 0..n {
            poll(&mut cancelled)?;
            for order in 0..3 {
                let mut action=ComplexDd::ZERO;
                for j in 0..n {
                    let coefficient=scale(self.original.source[order][i*n+j],
                        -coordinates[i]-coordinates[j]-order as i32*time)?;
                    action=action.add(ComplexDd::from_complex(y[j]).scale(coefficient));
                }
                for _ in order..2 {action=action.mul(ComplexDd::from_complex(self.s));}
                actions[order][i]=scale_complex(action.complex(),coordinates[i]+2*time+amplitude)?;
            }
        }
        let balance=evaluate_harmonic_balance(applied,&vec![C64::ZERO;n],
            &actions[0],&actions[1],&actions[2],&velocity,
            self.options.maximum_relative_residual,&mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(QuadraticForcedResponse {displacement,velocity,balance})
    }
}

#[cfg(test)]
mod tests;
