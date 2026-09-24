//! Shared stationary baffled-surface pressure from physical velocity histories.
//! The existing piano receiver and moving-valve outlets use this one owner.
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
use crate::pcm_wav::decimate as antialias;
use fs_math::det;

/// Largest retained input basis (the established observer/decimator limit).
pub const MAX_RAYLEIGH_INPUTS: usize = 128;

/// A positive-area quadrature point on the stationary z=0 baffle.
/// Normal velocity is the dot product of `mode_shape` with the input vector.
/// For a unit volume-flow input on area A, use a shape coefficient of 1/A.
#[derive(Clone, Debug)]
pub struct SurfaceSample {
    /// Position in metres, in the common fixed baffle frame.
    pub position_m: [f64; 3],
    /// Actual quadrature area [m²], not a normalized mixing weight.
    pub area_m2: f64,
    /// Signed normal-velocity coefficients for the complete input basis.
    pub mode_shape: Vec<f64>,
}

/// Homogeneous exterior medium; separate from a solid material state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RayleighMedium {
    /// Density [kg/m³].
    pub density: f64,
    /// Sound speed [m/s].
    pub sound_speed: f64,
}

const MAX_HISTORY: usize = 32_768;

/// Explicit geometry and numerical resolution of an infinite-baffle circular
/// outlet. The center is the origin and outward normal is +z. Its one input is
/// total outward volume flow Q [m³/s], with uniform normal velocity Q/(pi a²).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CircularOutletReceiver {
    /// Stationary receiver relative to the outlet plane [m].
    pub position_m: [f64; 3],
    /// Equal-area radial annuli, 1..=128.
    pub radial_rings: usize,
    /// Uniform azimuthal samples per annulus, 8..=512.
    pub angular_points: usize,
    /// Declared resolved band upper limit [Hz]. At most 0.1 of the mechanical
    /// rate, limiting differencing/interpolation errors. Not a spectral filter
    /// or proof that a nonlinear source has no content above this band.
    pub maximum_frequency_hz: f64,
}


/// Time-major modal convolution. All allocations and geometry operations are
/// cold; one output step has a bounded multiply-add count and no allocation.
pub struct BaffledPressure {
    kernels: Vec<(usize, Vec<f64>)>,
    /// Interleaved frames with exactly `modes` entries, not MAX_RAYLEIGH_INPUTS
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
impl BaffledPressure {
    /// Prepare the existing Rayleigh integral over caller-supplied velocity shapes.
    /// `rate` is the observation clock. All allocations happen at construction.
    /// The surface must lie at z=0; the receiver is at least 5 cm into z>0.
    /// This is prescribed-flow radiation, NOT a radiation-reaction load.
    ///
    /// # Errors
    /// Invalid complete geometry/medium/basis, sample rate or history budget.
    pub fn new(surface: &[SurfaceSample], inputs: usize, rate: u32,
        position_m: [f64; 3], medium: RayleighMedium) -> Result<Self, String> {
        Self::from_loaded(surface, inputs, rate, position_m, medium)
    }

    /// Observe every mechanics-frame velocity before decimating acceleration.
    /// The input and output clocks must have an integer ratio in 1..=16.
    ///
    /// # Errors
    /// Clock mismatch or the same geometry/admission refusals as `new`.
    pub fn new_multirate(surface: &[SurfaceSample], inputs: usize,
        mechanics_rate: u32, output_rate: u32, position_m: [f64; 3],
        medium: RayleighMedium) -> Result<Self, String> {
        if output_rate == 0 || mechanics_rate == 0 || mechanics_rate % output_rate != 0 {
            return Err("microphone requires an integer mechanics/audio rate ratio".into());
        }
        let mut receiver = Self::new(surface, inputs, output_rate, position_m, medium)?;
        receiver.set_input_ratio((mechanics_rate / output_rate) as usize)?;
        Ok(receiver)
    }
    /// Observe a uniform circular outlet's actual volume flow in exterior Pa.
    /// Uses this same surface kernel, not a compact pressure gain. Positive-area
    /// midpoint quadrature retains off-axis travel times and interference.
    /// The existing terminal load remains separate: this is a ONE-WAY prescribed
    /// flow observer, not an automatically matched radiation impedance.
    ///
    /// # Errors
    /// Invalid geometry, medium, count/band budget or cell phase span > pi/4
    /// in the declared band. Spatial and temporal refinement remain necessary
    /// for quantitative accuracy; no continuum error enclosure is claimed.
    pub fn circular_outlet(radius_m: f64, mechanical_rate: u32,
        receiver: CircularOutletReceiver, medium: RayleighMedium) -> Result<Self,String> {
        if !radius_m.is_finite() || radius_m<=0.0
            || !(1..=128).contains(&receiver.radial_rings)
            || !(8..=512).contains(&receiver.angular_points)
            || !receiver.maximum_frequency_hz.is_finite() || receiver.maximum_frequency_hz<=0.0
            || receiver.maximum_frequency_hz>0.1*f64::from(mechanical_rate)
            || !medium.sound_speed.is_finite() || medium.sound_speed<=0.0 {
            return Err("circular receiver requires physical geometry and explicit bounded resolution".into());
        }
        let area=std::f64::consts::PI*radius_m*radius_m;
        let cell=radius_m/det::sqrt(receiver.radial_rings as f64)
            +std::f64::consts::TAU*radius_m/receiver.angular_points as f64;
        let phase=std::f64::consts::TAU*receiver.maximum_frequency_hz*cell/medium.sound_speed;
        if !area.is_finite() || area<=0.0 || !(1.0/area).is_finite()
            || !phase.is_finite() || phase>std::f64::consts::FRAC_PI_4 {
            return Err("outlet area is unrepresentable or quadrature underresolves the declared acoustic band".into());
        }
        let points=receiver.radial_rings*receiver.angular_points;
        let patch_area=area/points as f64;
        let mut surface=Vec::with_capacity(points);
        for ring in 0..receiver.radial_rings {
            let r=radius_m*det::sqrt((ring as f64+0.5)/receiver.radial_rings as f64);
            for angle in 0..receiver.angular_points {
                let theta=std::f64::consts::TAU*(angle as f64+0.5)/receiver.angular_points as f64;
                surface.push(SurfaceSample { position_m:[r*det::cos(theta),r*det::sin(theta),0.0],
                    area_m2:patch_area,mode_shape:vec![1.0/area] });
            }
        }
        Self::new(&surface,1,mechanical_rate,receiver.position_m,medium)
    }

    fn set_input_ratio(&mut self, ratio: usize) -> Result<(), String> {
        self.decimator = antialias::Decimator::new(ratio, self.modes)?;
        self.acceleration = vec![0.0; ratio * self.modes];
        Ok(())
    }
    pub fn filter_delay_samples(&self) -> f64 { self.decimator.delay_output_frames() }
    fn from_loaded(surface: &[SurfaceSample], modes: usize, rate: u32,
        position_m: [f64; 3], medium: RayleighMedium) -> Result<Self,String> {
        if !(1..=MAX_RAYLEIGH_INPUTS).contains(&modes) || !(8_000..=192_000).contains(&rate)
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
    fn medium()->RayleighMedium {RayleighMedium{density:1.2,sound_speed:320.0}}
    fn patch(x:f64,shape:f64)->SurfaceSample {
        SurfaceSample{position_m:[x,0.0,0.0],area_m2:0.1,mode_shape:vec![shape]}
    }
    #[test]
    fn known_acceleration_arrives_after_propagation_with_physical_pressure_units() {
        let mut mic=BaffledPressure::from_loaded(&[patch(0.0,1.0)],1,32_000,[0.,0.,1.],medium()).unwrap();
        let gain=1.2*0.1/(2.0*std::f64::consts::PI);
        for n in 0..200 {
            let p=mic.step(&[(n+1) as f64/32_000.0]).unwrap();
            if n<100 {assert_eq!(p,0.0);}else{assert!((p/gain-1.0).abs()<1e-12);}
        }
    }
    #[test]
    fn opposite_surface_motion_cancels_on_axis_but_not_off_axis() {
        let surface=[patch(-0.3,1.0),patch(0.3,-1.0)];
        let mut center=BaffledPressure::from_loaded(&surface,1,48_000,[0.,0.,1.],medium()).unwrap();
        let mut side=BaffledPressure::from_loaded(&surface,1,48_000,[0.8,0.,1.],medium()).unwrap();
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
        let mut a=BaffledPressure::from_loaded(&surface,1,48_000,[0.,0.,1.],medium()).unwrap();
        let mut b=BaffledPressure::from_loaded(&surface,1,48_000,[0.,0.,1.],medium()).unwrap();
        for n in 0..600 {
            assert!(a.step(&[f64::NAN]).is_err());
            let v=(n as f64*0.08).sin();
            assert_eq!(a.step(&[v]).unwrap()*2.0,b.step(&[2.0*v]).unwrap());
        }
        assert!(BaffledPressure::from_loaded(&surface,1,48_000,[0.,0.,0.],medium()).is_err());
    }

    #[test]
    fn highest_retained_mode_reaches_pressure_without_a_32_mode_truncation() {
        for modes in [1,33,64,MAX_RAYLEIGH_INPUTS] {
            let mut shape=vec![0.0;modes];shape[modes-1]=1.0;
            let sample=SurfaceSample{position_m:[0.,0.,0.],area_m2:0.1,mode_shape:shape};
            let mut mic=BaffledPressure::from_loaded(&[sample],modes,32_000,[0.,0.,1.],medium()).unwrap();
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
        let too_many=SurfaceSample{position_m:[0.,0.,0.],area_m2:0.1,mode_shape:vec![1.;MAX_RAYLEIGH_INPUTS+1]};
        assert!(BaffledPressure::from_loaded(&[too_many],MAX_RAYLEIGH_INPUTS+1,32_000,[0.,0.,1.],medium()).is_err());
    }
    #[test]
    fn rejected_acceleration_does_not_advance_any_mode_or_delay() {
        let sample=SurfaceSample{position_m:[0.,0.,0.],area_m2:0.1,mode_shape:vec![1.;64]};
        let mut mic=BaffledPressure::from_loaded(&[sample],64,48_000,[0.,0.,1.],medium()).unwrap();
        for _ in 0..50 {mic.step(&[0.1;64]).unwrap();}
        let history=mic.history.clone();let previous=mic.previous_velocity.clone();let head=mic.head;
        assert!(mic.step(&[f64::MAX;64]).is_err());
        assert_eq!(mic.history,history);assert_eq!(mic.previous_velocity,previous);assert_eq!(mic.head,head);
    }
    #[test]
    fn an_intrasample_impact_is_not_erased_by_endpoint_downsampling() {
        let mut mic=BaffledPressure::from_loaded(&[patch(0.0,1.0)],1,32_000,[0.,0.,1.],medium()).unwrap();
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
            let mut a=BaffledPressure::from_loaded(&surface,1,32_000,[0.,0.,1.],medium()).unwrap();
            let mut b=BaffledPressure::from_loaded(&surface,1,32_000,[0.,0.,1.],medium()).unwrap();
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

#[cfg(test)]
mod circular_tests {
    use super::*;
    fn config() -> CircularOutletReceiver { CircularOutletReceiver {
        position_m:[0.0,0.0,0.2],radial_rings:32,angular_points:64,maximum_frequency_hz:1000.0,
    } }
    #[test]
    fn disk_flow_matches_the_independent_on_axis_rayleigh_integral() {
        let (radius,rate,rho,c,f)=(0.025,48000,1.2,320.0,400.0);
        let config=config();
        let mut mic=BaffledPressure::circular_outlet(radius,rate,config,RayleighMedium{density:rho,sound_speed:c}).unwrap();
        let area=std::f64::consts::PI*radius*radius;
        let w=std::f64::consts::TAU*f; let q=1e-5;
        let (near,far)=(config.position_m[2]/c,(config.position_m[2].powi(2)+radius*radius).sqrt()/c);
        let mut error=0.0;let mut signal=0.0;
        for n in 0..4800 {
            let t=n as f64/f64::from(rate);
            let p=mic.step(&[q*(w*t).sin()]).unwrap();
            if n>1000 {
                // Differencing returns interval-average acceleration centered half
                // a sample before this source sample. The independent continuous
                // on-axis integral eliminates the surface integration entirely.
                let tc=t-0.5/f64::from(rate);
                let exact=rho*c/area*q*((w*(tc-near)).sin()-(w*(tc-far)).sin());
                error+=(p-exact).powi(2); signal+=exact.powi(2);
            }
        }
        assert!((error/signal).sqrt()<0.002, "relative RMS {}",(error/signal).sqrt());
    }
    #[test]
    fn prescribed_flow_has_causal_distance_dependent_pressure_and_no_dc_sound() {
        let medium=RayleighMedium{density:1.2,sound_speed:320.0};
        let mut close=BaffledPressure::circular_outlet(0.007,48000,config(),medium).unwrap();
        let mut distant=config();distant.position_m[2]=0.4;
        let mut far=BaffledPressure::circular_outlet(0.007,48000,distant,medium).unwrap();
        assert!(far.delay_samples.0>close.delay_samples.1);
        let (mut p1,mut p2)=(0.0_f64,0.0_f64);
        for n in 0..600 {
            let a=close.step(&[1e-5]).unwrap();let b=far.step(&[1e-5]).unwrap();
            if n<close.delay_samples.0 {assert_eq!(a,0.0);}
            p1=p1.max(a.abs());p2=p2.max(b.abs());
            if n>far.delay_samples.1+2 {assert_eq!(a,0.0);assert_eq!(b,0.0);}
        }
        assert!(p1>0.0 && p2>0.0);
    }
    #[test]
    fn coarse_or_unrepresentable_outlets_refuse_instead_of_inventing_a_point_source() {
        let medium=RayleighMedium{density:1.2,sound_speed:320.0};
        let mut coarse=config();coarse.radial_rings=1;coarse.angular_points=8;
        assert!(BaffledPressure::circular_outlet(0.2,48000,coarse,medium).is_err());
        for radius in [0.0,-1.0,f64::NAN,f64::MAX] {
            assert!(BaffledPressure::circular_outlet(radius,48000,config(),medium).is_err());
        }
        let mut high=config();high.maximum_frequency_hz=4801.0;
        assert!(BaffledPressure::circular_outlet(0.007,48000,high,medium).is_err());
    }
}
