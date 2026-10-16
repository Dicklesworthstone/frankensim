//! Intrinsic loss in the same moving-boundary coordinates as the strings.
//!
//! With z=q+beta*b, distributed velocity drag is diagonal in z plus residual
//! endpoint-mass drag. Strain-rate loss acts on q_dot=z_dot-beta*b_dot, so it
//! must react on the bridge too. The nonnegative power is
//! sum(a*z_dot²+d*(z_dot-beta*b_dot)²)+endpoint_drag*b_dot².
//! The RT-0425 reduction additionally includes viscous T/L chord traction;
//! bending alone cannot damp an uncurved moving chord.
//!
//! For the scalar Euler--Bernoulli reduction of Chabassier et al. (2014),
//! equations (3.68), (3.71), (3.73), take eta_phi=eta_u and omit rotary inertia:
//! https://www.numdam.org/item/10.1051/m2an/2013136.pdf
//!
//! A symmetric midpoint velocity flow surrounds the conservative bank step.
//! Its small board Schur complement is formed cold with positive terms; the
//! hot path is allocation-free and groups strings sharing one bridge shape.
//! This is a passive second-order flow, not an exact damping exponential.

use super::{StringMode, StringPort, identity, inverse_spd};
use std::ops::Range;

#[derive(Clone, Copy, Debug)]
pub struct ModeLoss {
    /// Distributed absolute-velocity rate, 2 R [1/s].
    pub mass_rate: f64,
    /// Relative strain-rate coefficient [1/s].
    pub strain_rate: f64,
}

#[derive(Clone, Debug)]
pub struct Prepared {
    pub modes: Vec<ModeLoss>,
    /// Per-span residual-mass drag plus any viscous chord traction [kg/s].
    pub endpoint_drag: Vec<f64>,
    duration_s: f64,
    enabled: bool,
    inverse: Vec<f64>,
    cross_response: Vec<f64>,
    schur_inverse: Vec<f64>,
    rhs: Vec<f64>,
    board_midpoint: Vec<f64>,
    board_force: Vec<f64>,
}

impl Prepared {
    #[allow(clippy::too_many_arguments)]
    pub fn new(strings: &[StringPort], modes: &[StringMode], groups: &[Range<usize>],
        losses: Vec<ModeLoss>, endpoint_drag: Vec<f64>, board_count: usize,
        duration_s: f64) -> Result<Self, String> {
        if losses.len()!=modes.len() || endpoint_drag.len()!=strings.len()
            || !duration_s.is_finite() || duration_s<=0.0
            || losses.iter().any(|c|!c.mass_rate.is_finite() || c.mass_rate<0.0
                || !c.strain_rate.is_finite() || c.strain_rate<0.0)
            || endpoint_drag.iter().any(|c|!c.is_finite() || *c<0.0) {
            return Err("invalid physical string damping coefficients".into());
        }
        let enabled=losses.iter().any(|c|c.mass_rate!=0.0 || c.strain_rate!=0.0)
            || endpoint_drag.iter().any(|c|*c!=0.0);
        let half=0.5*duration_s;
        let mut inverse=vec![1.0;modes.len()];
        let mut cross_response=vec![0.0;modes.len()];
        for (k,c) in losses.iter().enumerate() {
            inverse[k]=1.0/(1.0+half*(c.mass_rate+c.strain_rate));
            cross_response[k]=half*c.strain_rate*modes[k].beta*inverse[k];
        }
        let r=board_count;
        let mut schur=identity(r);
        for group in groups {
            let g=&strings[group.start].bridge;
            let mut weight=0.0;
            for si in group.clone() {
                weight+=half*endpoint_drag[si];
                for k in strings[si].modes.clone() {
                    let c=losses[k];
                    // d*beta² - half*d²*beta²/(1+half*(a+d)), without
                    // cancellation between two large nonnegative terms.
                    weight+=half*c.strain_rate*modes[k].beta.powi(2)
                        *(1.0+half*c.mass_rate)*inverse[k];
                }
            }
            for i in 0..r {for j in i..r {
                let value=weight*g[i]*g[j];
                schur[i*r+j]+=value;
                if i!=j {schur[j*r+i]+=value;}
            }}
        }
        if inverse.iter().chain(&cross_response).chain(&schur).any(|x|!x.is_finite())
            || inverse.iter().any(|x|*x<=0.0) {
            return Err("physical string damping preparation overflow".into());
        }
        let schur_inverse=if enabled {inverse_spd(&schur,r)?} else {schur};
        Ok(Self {modes:losses,endpoint_drag,duration_s,enabled,inverse,cross_response,
            schur_inverse,rhs:vec![0.0;r],board_midpoint:vec![0.0;r],board_force:vec![0.0;r]})
    }

    /// Continuous physical dissipation [W], also used to audit harmonic loss.
    pub fn power(&self, strings: &[StringPort], modes: &[StringMode], v: &[f64]) -> f64 {
        let n=modes.len();
        let mut power=0.0;
        for (si,s) in strings.iter().enumerate() {
            let b=s.bridge.iter().zip(&v[n..]).map(|(g,v)|g*v).sum::<f64>();
            power+=self.endpoint_drag[si]*b*b;
            for k in s.modes.clone() {
                let c=self.modes[k];let relative=v[k]-modes[k].beta*b;
                power+=c.mass_rate*v[k]*v[k]+c.strain_rate*relative*relative;
            }
        }
        power
    }

    /// Apply v'=-C v over the prepared duration using implicit midpoint.
    /// Return the actual kinetic-energy decrease in the bank coordinates.
    /// All bridge reactions use the SAME midpoint as their string forces.
    pub fn apply(&mut self, strings: &[StringPort], modes: &[StringMode],
        groups: &[Range<usize>], v: &mut [f64]) -> f64 {
        if !self.enabled {return 0.0;}
        let n=modes.len();let r=self.rhs.len();
        self.rhs.copy_from_slice(&v[n..]);
        for group in groups {
            let g=&strings[group.start].bridge;
            let mut source=0.0;
            for si in group.clone() {for k in strings[si].modes.clone() {
                source+=self.cross_response[k]*v[k];
            }}
            for (rhs,g) in self.rhs.iter_mut().zip(g) {*rhs+=g*source;}
        }
        for i in 0..r {
            self.board_midpoint[i]=(0..r).map(|j|self.schur_inverse[i*r+j]*self.rhs[j]).sum();
        }
        self.board_force.fill(0.0);
        let mut loss=0.0;
        for group in groups {
            let g=&strings[group.start].bridge;
            let b=g.iter().zip(&self.board_midpoint).map(|(g,v)|g*v).sum::<f64>();
            let mut reaction=0.0;
            for si in group.clone() {
                reaction+=self.endpoint_drag[si]*b;
                for k in strings[si].modes.clone() {
                    let c=self.modes[k];let old=v[k];
                    let midpoint=self.inverse[k]*old+self.cross_response[k]*b;
                    let relative=midpoint-modes[k].beta*b;
                    v[k]=old-self.duration_s*(c.mass_rate*midpoint+c.strain_rate*relative);
                    loss+=0.5*(old-v[k])*(old+v[k]);
                    reaction-=c.strain_rate*modes[k].beta*relative;
                }
            }
            for (force,g) in self.board_force.iter_mut().zip(g) {*force+=g*reaction;}
        }
        for j in 0..r {
            let k=n+j;let old=v[k];
            v[k]=old-self.duration_s*self.board_force[j];
            loss+=0.5*(old-v[k])*(old+v[k]);
        }
        loss
    }
}

#[cfg(test)]
#[path = "string_damping_tests.rs"]
mod tests;
