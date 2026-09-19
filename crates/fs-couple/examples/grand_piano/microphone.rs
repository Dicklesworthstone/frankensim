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
use super::linear::{Bank, MAX_BOARD_MODES};
use fs_bem::helmholtz::Medium;
use fs_math::det;

#[path = "antialias.rs"]
mod antialias;

const MAX_HISTORY: usize = 32_768;

/// Time-major modal convolution. All allocations and geometry operations are
/// cold; one output step has a bounded multiply-add count and no allocation.
pub struct Microphone {
    kernels: Vec<(usize, Vec<f64>)>,
    /// Interleaved frames with exactly `modes` entries, not MAX_BOARD_MODES
    /// padding. Increasing the admitted bandwidth does not multiply memory
    /// for an existing one-mode or four-mode soundboard.
    history: Vec<f64>,
    history_frames: usize,
    previous_velocity: Vec<f64>,
    acceleration: Vec<f64>,
    decimator: antialias::Decimator,
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
    /// Production path: use every mechanics substep, not just its final velocity.
    pub fn new_multirate(surface: &[SurfaceSample], bank: &Bank, rate: u32,
        position_m: [f64; 3], medium: Medium) -> Result<Self,String> {
        if rate == 0 || bank.rate % rate != 0 {
            return Err("microphone requires an integer mechanics/audio rate ratio".into());
        }
        let mut mic = Self::new(surface, bank, rate, position_m, medium)?;
        mic.set_input_ratio((bank.rate / rate) as usize)?;
        Ok(mic)
    }
    fn set_input_ratio(&mut self, ratio: usize) -> Result<(), String> {
        self.decimator = antialias::Decimator::new(ratio, self.modes)?;
        self.acceleration = vec![0.0; ratio * self.modes];
        Ok(())
    }
    pub fn filter_delay_samples(&self) -> f64 { self.decimator.delay_output_frames() }
    fn from_loaded(surface: &[SurfaceSample], modes: usize, rate: u32,
        position_m: [f64; 3], medium: Medium) -> Result<Self,String> {
        if !(1..=MAX_BOARD_MODES).contains(&modes) || !(8_000..=192_000).contains(&rate)
            || surface.is_empty() || surface.len()>120_000
            || position_m.iter().any(|v|!v.is_finite()) || position_m[2]<0.05
            || !medium.density.is_finite() || medium.density<=0.0
            || !medium.sound_speed.is_finite() || medium.sound_speed<=0.0 {
            return Err("invalid surface/microphone/medium budget; microphone must be >=5 cm above z=0".into());
        }
        let mut bins=std::collections::BTreeMap::<usize,Vec<f64>>::new();
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
                let row=bins.entry(d).or_insert_with(||vec![0.0;modes]);
                for (i,value) in p.mode_shape.iter().enumerate() {row[i]+=gain*weight*value;}
            }
        }
        if bins.values().flatten().any(|v|!v.is_finite()) {return Err("radiation kernel overflow".into());}
        let first=*bins.first_key_value().ok_or("empty radiation kernel")?.0;
        let last=*bins.last_key_value().ok_or("empty radiation kernel")?.0;
        let history_frames=last+1;
        let samples=history_frames.checked_mul(modes).ok_or("microphone history size overflow")?;
        Ok(Self { kernels: bins.into_iter().collect(), history: vec![0.0;samples],history_frames,
            previous_velocity:vec![0.0;modes],acceleration:vec![0.0;modes],head:0,modes,rate:f64::from(rate),
            decimator:antialias::Decimator::new(1,modes)?,
            position_m,delay_samples:(first,last) })
    }
    pub fn multiply_adds_per_sample(&self)->usize {
        self.kernels.len()*self.modes + self.decimator.multiplies_per_output_frame()
    }

    /// Legacy one-rate input. Multirate callers must supply the complete trace.
    pub fn step(&mut self, velocity: &[f64])->Result<f64,String> {
        self.step_trace(velocity)
    }

    /// Interleaved END velocities from EVERY accepted mechanics substep.
    /// Difference at mechanics rate, low-pass/decimate acceleration, then apply
    /// spatial propagation at audio rate. The filter adds explicit causal delay.
    /// All candidate values are checked before either delay line is committed.
    pub fn step_trace(&mut self, velocity: &[f64])->Result<f64,String> {
        if velocity.len()!=self.acceleration.len() || velocity.iter().any(|v|!v.is_finite()) {
            return Err("microphone needs the complete finite substep velocity trace".into());
        }
        let ratio=self.decimator.input_frames();
        let mechanics_rate=self.rate*ratio as f64;
        for s in 0..ratio { for i in 0..self.modes {
            let previous=if s==0 {self.previous_velocity[i]} else {velocity[(s-1)*self.modes+i]};
            self.acceleration[s*self.modes+i]=(velocity[s*self.modes+i]-previous)*mechanics_rate;
        }}
        let filtered=self.decimator.preview(&self.acceleration)?;
        let mut pressure=0.0;
        for (delay,weights) in &self.kernels {
            let frame=if self.head>=*delay {self.head-*delay}else{self.head+self.history_frames-*delay};
            let offset=frame*self.modes;
            for (i,weight) in weights.iter().enumerate() {
                pressure+=weight*self.history[offset+i];
            }
        }
        if !pressure.is_finite() {return Err("radiated pressure overflow".into());}
        let offset=self.head*self.modes;
        self.history[offset..offset+self.modes].copy_from_slice(filtered);
        self.previous_velocity.copy_from_slice(&velocity[(ratio-1)*self.modes..]);
        self.decimator.commit();
        self.head+=1;if self.head==self.history_frames{self.head=0;}
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

    #[test]
    fn highest_retained_mode_reaches_pressure_without_a_32_mode_truncation() {
        for modes in [1,33,64,MAX_BOARD_MODES] {
            let mut shape=vec![0.0;modes];shape[modes-1]=1.0;
            let sample=SurfaceSample{position_m:[0.,0.,0.],area_m2:0.1,mode_shape:shape};
            let mut mic=Microphone::from_loaded(&[sample],modes,32_000,[0.,0.,1.],medium()).unwrap();
            assert_eq!(mic.history.len(),mic.history_frames*modes);
            assert!(mic.kernels.iter().all(|(_,row)|row.len()==modes));
            let pointers=(mic.history.as_ptr(),mic.acceleration.as_ptr(),mic.previous_velocity.as_ptr());
            let mut velocity=vec![0.0;modes];
            let gain=1.2*0.1/(2.0*std::f64::consts::PI);
            for n in 0..600 {
                velocity[modes-1]=(n+1) as f64/32_000.0;
                let p=mic.step(&velocity).unwrap();
                if n<100 {assert_eq!(p,0.0);}else{assert!((p/gain-1.0).abs()<1e-11);}
            }
            assert_eq!(pointers,(mic.history.as_ptr(),mic.acceleration.as_ptr(),mic.previous_velocity.as_ptr()));
        }
        let too_many=SurfaceSample{position_m:[0.,0.,0.],area_m2:0.1,mode_shape:vec![1.;MAX_BOARD_MODES+1]};
        assert!(Microphone::from_loaded(&[too_many],MAX_BOARD_MODES+1,32_000,[0.,0.,1.],medium()).is_err());
    }
    #[test]
    fn rejected_acceleration_does_not_advance_any_mode_or_delay() {
        let sample=SurfaceSample{position_m:[0.,0.,0.],area_m2:0.1,mode_shape:vec![1.;64]};
        let mut mic=Microphone::from_loaded(&[sample],64,48_000,[0.,0.,1.],medium()).unwrap();
        for _ in 0..50 {mic.step(&[0.1;64]).unwrap();}
        let history=mic.history.clone();let previous=mic.previous_velocity.clone();let head=mic.head;
        assert!(mic.step(&[f64::MAX;64]).is_err());
        assert_eq!(mic.history,history);assert_eq!(mic.previous_velocity,previous);assert_eq!(mic.head,head);
    }
    #[test]
    fn an_intrasample_impact_is_not_erased_by_endpoint_downsampling() {
        let mut mic=Microphone::from_loaded(&[patch(0.0,1.0)],1,32_000,[0.,0.,1.],medium()).unwrap();
        mic.set_input_ratio(4).unwrap();
        assert_eq!(mic.filter_delay_samples(),44.0);
        let mut peak=0.0_f64;
        for frame in 0..400 {
            // The audio endpoint is zero, but a real within-sample motion occurred.
            let trace=if frame==0 {[0.01,0.0,0.0,0.0]} else {[0.0;4]};
            peak=peak.max(mic.step_trace(&trace).unwrap().abs());
        }
        assert!(peak>1e-5,"lost the substep impact");
    }
    #[test]
    fn multirate_pressure_keeps_dc_units_and_rejected_input_does_not_change_time() {
        for ratio in [1,3,4,8,16] {
            let surface=[patch(0.0,1.0)];
            let mut a=Microphone::from_loaded(&surface,1,32_000,[0.,0.,1.],medium()).unwrap();
            let mut b=Microphone::from_loaded(&surface,1,32_000,[0.,0.,1.],medium()).unwrap();
            a.set_input_ratio(ratio).unwrap();b.set_input_ratio(ratio).unwrap();
            let gain=1.2*0.1/(2.0*std::f64::consts::PI);
            for frame in 0..500 {
                assert!(a.step_trace(&vec![f64::MAX;ratio]).is_err());
                let trace:Vec<f64>=(0..ratio).map(|s|(frame*ratio+s+1) as f64/(32_000*ratio) as f64).collect();
                let p=a.step_trace(&trace).unwrap();assert_eq!(p,b.step_trace(&trace).unwrap());
                if frame>300 {assert!((p/gain-1.0).abs()<1e-10);}
            }
        }
    }

}
