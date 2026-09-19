//! Spatial, retarded-time pressure observation of the actual soundboard modes.
//!
//! The time-domain Rayleigh half-space integral is
//! p(x,t)=rho/(2*pi) integral_S a_n(y,t-|x-y|/c)/|x-y| dS.
//! This is the distributed counterpart of fs-couple's compact plate observer,
//! not a second plate solver or a sampled room response. Geometry, signed mode
//! shapes, phase cancellation, distance and propagation time determine the FIR.
//! See Wood, Rayleigh integral eq. (2):
//! https://euphonics.org/4-3-2-the-rayleigh-integral-and-the-baffled-piston/
//!
//! IMPORTANT MODEL: an infinite rigid baffle around a flat panel, prescribed
//! motion and stationary microphone. This is NOT full piano/room/lid scattering
//! or acoustic backreaction. Both sides of an unbaffled board need exterior BEM.
//! Finite differences give interval-average acceleration; linear fractional
//! delays have less than one-sample wavefront uncertainty. Mesh quadrature,
//! structural bandwidth and output sampling still limit acoustic accuracy.
//! This is not the previous arbitrary volume-velocity gain in Pa/(m^3/s).
use super::board_geometry::SurfaceSample;
use super::linear::Bank;
use fs_bem::helmholtz::Medium;
use fs_math::det;

const MAX_MODES: usize = 32;
const MAX_HISTORY: usize = 32_768;

/// Time-major modal convolution. All allocations and geometry operations are
/// cold; one output step has a bounded multiply-add count and no allocation.
pub struct Microphone {
    kernels: Vec<(usize, [f64; MAX_MODES])>,
    history: Vec<[f64; MAX_MODES]>,
    previous_velocity: [f64; MAX_MODES],
    head: usize,
    modes: usize,
    rate: f64,
    pub position_m: [f64; 3],
    pub delay_samples: (usize, usize),
}
impl Microphone {
    pub fn new(surface: &[SurfaceSample], bank: &Bank, rate: u32,
        position_m: [f64; 3], medium: Medium) -> Result<Self,String> {
        // This is the very same coordinate map used for bridge force/flow,
        // retained by Bank rather than reconstructed by another eigensolve.
        let loaded = surface.iter().map(|p| Ok(SurfaceSample {
            position_m: p.position_m, area_m2: p.area_m2,
            mode_shape: bank.project_board_shape(&p.mode_shape)?,
        })).collect::<Result<Vec<_>,String>>()?;
        Self::from_loaded(&loaded, bank.board_count, rate, position_m, medium)
    }
    fn from_loaded(surface: &[SurfaceSample], modes: usize, rate: u32,
        position_m: [f64; 3], medium: Medium) -> Result<Self,String> {
        if !(1..=MAX_MODES).contains(&modes) || !(8_000..=192_000).contains(&rate)
            || surface.is_empty() || surface.len()>120_000
            || position_m.iter().any(|v|!v.is_finite()) || position_m[2]<0.05
            || !medium.density.is_finite() || medium.density<=0.0
            || !medium.sound_speed.is_finite() || medium.sound_speed<=0.0 {
            return Err("invalid surface/microphone/medium budget; microphone must be >=5 cm above z=0".into());
        }
        let mut bins=std::collections::BTreeMap::<usize,[f64;MAX_MODES]>::new();
        for p in surface {
            if !p.area_m2.is_finite() || p.area_m2<=0.0
                || p.position_m.iter().any(|v|!v.is_finite()) || p.position_m[2]!=0.0
                || p.mode_shape.len()!=modes || p.mode_shape.iter().any(|v|!v.is_finite()) {
                return Err("Rayleigh sample must be a finite flat-panel area with all modal shapes".into());
            }
            let radius=det::sqrt((0..3).map(|i|(position_m[i]-p.position_m[i]).powi(2)).sum());
            let delay=radius*f64::from(rate)/medium.sound_speed;
            if !delay.is_finite() || delay<1.0 || delay>=(MAX_HISTORY-2) as f64 {
                return Err("microphone propagation delay exceeds the causal history budget".into());
            }
            let first=delay.floor() as usize; let fraction=delay-first as f64;
            let gain=medium.density*p.area_m2/(2.0*std::f64::consts::PI*radius);
            for (d,weight) in [(first,1.0-fraction),(first+1,fraction)] {
                let row=bins.entry(d).or_insert([0.0;MAX_MODES]);
                for (i,value) in p.mode_shape.iter().enumerate() {row[i]+=gain*weight*value;}
            }
        }
        if bins.values().flatten().any(|v|!v.is_finite()) {return Err("radiation kernel overflow".into());}
        let first=*bins.first_key_value().ok_or("empty radiation kernel")?.0;
        let last=*bins.last_key_value().ok_or("empty radiation kernel")?.0;
        Ok(Self { kernels: bins.into_iter().collect(), history: vec![[0.0;MAX_MODES];last+1],
            previous_velocity:[0.0;MAX_MODES],head:0,modes,rate:f64::from(rate),
            position_m,delay_samples:(first,last) })
    }
    pub fn multiply_adds_per_sample(&self)->usize {self.kernels.len()*self.modes}

    /// The model starts from rest. Supply the loaded board's END velocity after
    /// each accepted audio sample. Differencing actual velocities includes the
    /// contact/coupling forces and damping; using only -omega^2*q would not.
    /// Input errors do not mutate history, and failed output rolls back its slot.
    pub fn step(&mut self, velocity: &[f64])->Result<f64,String> {
        if velocity.len()!=self.modes || velocity.iter().any(|v|!v.is_finite()) {
            return Err("microphone needs one finite velocity per loaded board mode".into());
        }
        let mut acceleration=[0.0;MAX_MODES];
        for i in 0..self.modes {acceleration[i]=(velocity[i]-self.previous_velocity[i])*self.rate;}
        if acceleration.iter().any(|v|!v.is_finite()) {return Err("surface acceleration overflow".into());}
        let old=self.history[self.head];self.history[self.head]=acceleration;
        let mut pressure=0.0;
        for (delay,weights) in &self.kernels {
            let index=if self.head>=*delay {self.head-*delay}else{self.head+self.history.len()-*delay};
            for (i,weight) in weights[..self.modes].iter().enumerate() {
                pressure+=weight*self.history[index][i];
            }
        }
        if !pressure.is_finite() {self.history[self.head]=old;return Err("radiated pressure overflow".into());}
        self.previous_velocity[..self.modes].copy_from_slice(velocity);
        self.head+=1;if self.head==self.history.len(){self.head=0;}
        Ok(pressure)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn medium()->Medium {Medium{density:1.2,sound_speed:320.0}}
    fn patch(x:f64,shape:f64)->SurfaceSample {
        SurfaceSample{position_m:[x,0.0,0.0],area_m2:0.1,mode_shape:vec![shape]}
    }
    #[test]
    fn known_acceleration_arrives_after_propagation_with_physical_pressure_units() {
        let mut mic=Microphone::from_loaded(&[patch(0.0,1.0)],1,32_000,[0.,0.,1.],medium()).unwrap();
        let gain=1.2*0.1/(2.0*std::f64::consts::PI);
        for n in 0..200 {
            let p=mic.step(&[(n+1) as f64/32_000.0]).unwrap();
            if n<100 {assert_eq!(p,0.0);}else{assert!((p/gain-1.0).abs()<1e-12);}
        }
    }
    #[test]
    fn opposite_surface_motion_cancels_on_axis_but_not_off_axis() {
        let surface=[patch(-0.3,1.0),patch(0.3,-1.0)];
        let mut center=Microphone::from_loaded(&surface,1,48_000,[0.,0.,1.],medium()).unwrap();
        let mut side=Microphone::from_loaded(&surface,1,48_000,[0.8,0.,1.],medium()).unwrap();
        let mut peak:f64=0.0;
        for n in 0..2000 {
            let velocity=[(n as f64*0.07).sin()];
            assert!(center.step(&velocity).unwrap().abs()<1e-12);
            peak=peak.max(side.step(&velocity).unwrap().abs());
        }
        assert!(peak>1.0);
    }
    #[test]
    fn pressure_is_linear_and_invalid_samples_leave_history_untouched() {
        let surface=[patch(0.0,1.0)];
        let mut a=Microphone::from_loaded(&surface,1,48_000,[0.,0.,1.],medium()).unwrap();
        let mut b=Microphone::from_loaded(&surface,1,48_000,[0.,0.,1.],medium()).unwrap();
        for n in 0..600 {
            assert!(a.step(&[f64::NAN]).is_err());
            let v=(n as f64*0.08).sin();
            assert_eq!(a.step(&[v]).unwrap()*2.0,b.step(&[2.0*v]).unwrap());
        }
        assert!(Microphone::from_loaded(&surface,1,48_000,[0.,0.,0.],medium()).is_err());
    }
}
