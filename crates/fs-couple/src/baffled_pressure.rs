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
//! The optional bandlimited constructors replace linear fractional delays with
//! a windowed sinc and explicitly account for its 32-sample extra latency.
//! Differentiating constructors additionally reconstruct acceleration with the
//! shared bandlimited derivative, retaining high partials without finite-
//! difference sinc attenuation. All numerical delay is reported, never hidden.
use crate::pcm_wav::decimate as antialias;
use fs_math::det;

#[path = "baffled_pressure_bandlimited.rs"]
mod bandlimited;

#[path = "baffled_pressure_path.rs"]
mod path;

/// Continuous-frequency pressure, gradient and particle velocity on this surface basis.
#[path = "baffled_pressure_harmonic.rs"]
pub mod harmonic;

/// Material-derived absorption/dispersion on every original surface path.
#[path = "baffled_pressure_gas.rs"]
pub mod gas;

/// Multiple physical openings through the same retarded-time pressure owner.
#[path = "baffled_pressure_openings.rs"]
pub mod openings;

/// Independently resolved physical tiles through bounded original receivers.
#[path = "baffled_pressure_tiles.rs"]
pub mod tiles;

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

impl From<&fs_material::gas::GasState> for RayleighMedium {
    /// Exact view of a caller-resolved medium, without a second gas evaluation.
    fn from(gas: &fs_material::gas::GasState) -> Self {
        Self {
            density: gas.density,
            sound_speed: gas.sound_speed,
        }
    }
}

const MAX_HISTORY: usize = 32_768;
/// Largest transposed tap ring (f64 entries, 16 MiB). Larger kernels keep the
/// direct form; both forms produce bit-identical pressure.
const MAX_TAP_RING: usize = 1 << 21;

/// One kernel's contribution from one history frame: products in input order,
/// four lanes, pairwise lane sum, then the scalar tail. The direct and the
/// transposed forms both call this, so each tap has one fixed rounding order.
#[inline(always)]
fn kernel_tap(weights: &[f64], frame: &[f64]) -> f64 {
    let mut s=[0.0_f64;4];
    let mut cw=weights.chunks_exact(4);
    let mut ch=frame.chunks_exact(4);
    for (w,h) in (&mut cw).zip(&mut ch) {
        s[0]+=w[0]*h[0];
        s[1]+=w[1]*h[1];
        s[2]+=w[2]*h[2];
        s[3]+=w[3]*h[3];
    }
    let mut tap=(s[0]+s[1])+(s[2]+s[3]);
    for (w,h) in cw.remainder().iter().zip(ch.remainder()) {
        tap+=w*h;
    }
    tap
}

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
    /// Declared resolved band upper limit [Hz]. Legacy finite-difference images
    /// admit at most 0.1 of their mechanical rate; differentiating images admit
    /// at most 0.45 of their output rate. This controls spatial resolution; it
    /// is not a filter or proof that a nonlinear source stays inside the band.
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
    /// Transposed evaluation, a pure layout choice. When frame f is committed,
    /// row f of this ring (parallel to `history`) receives tap_k(f) for every
    /// kernel k. Pressure at head t then sums tap_k(t-d_k) in kernel order.
    /// These are the same products, per-tap order and sum order as the direct
    /// form, so the result is bit-identical. Never-written rows hold +0.0,
    /// which equals a tap of the zero-initialized history.
    taps: Vec<f64>,
    /// Kernel rows packed contiguously in kernel order (transposed form only).
    packed_weights: Vec<f64>,
    /// Per kernel, (history_frames - delay) * kernels + k; reduce mod taps.len().
    tap_offsets: Vec<usize>,
    /// Address of the `history` allocation the ring mirrors; 0 = direct form.
    /// A constructor that replaces the kernels also replaces `history` (a new
    /// allocation while the old one is alive), so a stale ring is never used.
    tap_history: usize,
    previous_velocity: Vec<f64>,
    acceleration: Vec<f64>,
    decimator: antialias::Decimator,
    /// Cold choice: the prepared filter consumes velocity, not finite differences.
    differentiate_velocity: bool,
    /// `differentiate_pressure`: the bandlimited derivative on the summed pressure (one
    /// channel) while `decimator` passes velocity through.
    pressure_derivative: Option<antialias::Decimator>,
    /// Fixed delay introduced only by the optional bandlimited kernel.
    interpolation_delay: usize,
    /// Common causal padding of independently prepared material paths.
    material_path_delay: usize,
    head: usize,
    modes: usize,
    rate: f64,
    pub position_m: [f64; 3],
    /// FIR support, including fractional interpolation taps; not a group delay.
    pub delay_samples: (usize, usize),
}
/// A copy steps bit for bit as the original would from here. It keeps the
/// original's evaluation form: when the transposed ring mirrors the original's
/// history, the copy's ring names the copy's own history (a derived clone
/// would keep the original's address and fall back to the direct form).
impl Clone for BaffledPressure {
    fn clone(&self) -> Self {
        let transposed=self.tap_history==self.history.as_ptr() as usize;
        let mut copy=Self { kernels:self.kernels.clone(),history:self.history.clone(),
            history_frames:self.history_frames,taps:self.taps.clone(),
            packed_weights:self.packed_weights.clone(),tap_offsets:self.tap_offsets.clone(),tap_history:0,
            previous_velocity:self.previous_velocity.clone(),acceleration:self.acceleration.clone(),
            decimator:self.decimator.clone(),differentiate_velocity:self.differentiate_velocity,
            pressure_derivative:self.pressure_derivative.clone(),
            interpolation_delay:self.interpolation_delay,material_path_delay:self.material_path_delay,
            head:self.head,modes:self.modes,rate:self.rate,
            position_m:self.position_m,delay_samples:self.delay_samples };
        if transposed {copy.tap_history=copy.history.as_ptr() as usize;}
        copy
    }
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
    /// Bandlimited acceleration and travel-time reconstruction at one clock.
    /// The original constructors remain available as finite-difference references.
    ///
    /// # Errors
    /// The same geometry, medium, clock and causal-history admission as `new`.
    pub fn new_differentiating(surface: &[SurfaceSample], inputs: usize, rate: u32,
        position_m: [f64; 3], medium: RayleighMedium) -> Result<Self, String> {
        Self::new_multirate_differentiating(surface, inputs, rate, rate, position_m, medium)
    }

    /// Bandlimited SI acceleration from EVERY accepted mechanics-frame velocity,
    /// followed by the existing signed Rayleigh surface/travel-time kernel.
    /// The derivative is fused with the final anti-alias stage while cold. It
    /// never filters mechanical forces or changes the loaded surface basis.
    ///
    /// At 4x oversampling numerical latency is 44 + 8 + 32 = 84 output frames;
    /// at 1x it is 48 + 32 = 80. Physical propagation is additional. There is
    /// no two-sample derivative's implicit half-mechanics-frame time shift.
    /// Reconstruction is tested through 0.45 of the output rate, not at Nyquist.
    /// It cannot restore missing structural modes or certify continuum acoustics.
    ///
    /// # Errors
    /// Invalid complete geometry, medium, clocks or bounded history allocation.
    pub fn new_multirate_differentiating(surface: &[SurfaceSample], inputs: usize,
        mechanics_rate: u32, output_rate: u32, position_m: [f64; 3],
        medium: RayleighMedium) -> Result<Self, String> {
        let mut receiver = Self::new_multirate_bandlimited(surface, inputs, mechanics_rate,
            output_rate, position_m, medium)?;
        receiver.differentiate_velocity = true;
        receiver.set_input_ratio((mechanics_rate / output_rate) as usize)?;
        receiver.prepare_taps();
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
        if receiver.maximum_frequency_hz>0.1*f64::from(mechanical_rate) {
            return Err("circular receiver requires physical geometry and explicit bounded resolution".into());
        }
        let surface = Self::circular_surface(radius_m, receiver, medium)?;
        Self::new(&surface,1,mechanical_rate,receiver.position_m,medium)
    }

    // The same physical quadrature for every circular-outlet observation image.
    // Temporal admission belongs to its selected reconstruction, not the mesh.
    pub(crate) fn circular_surface(radius_m: f64, receiver: CircularOutletReceiver,
        medium: RayleighMedium) -> Result<Vec<SurfaceSample>, String>
    {
        if !radius_m.is_finite() || radius_m<=0.0
            || !(1..=128).contains(&receiver.radial_rings)
            || !(8..=512).contains(&receiver.angular_points)
            || !receiver.maximum_frequency_hz.is_finite() || receiver.maximum_frequency_hz<=0.0
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
        Ok(surface)
    }

    fn set_input_ratio(&mut self, ratio: usize) -> Result<(), String> {
        self.decimator = if self.differentiate_velocity {
            antialias::Decimator::new_differentiating(ratio, self.modes, self.rate as u32)?
        } else { antialias::Decimator::new(ratio, self.modes)? };
        self.acceleration = vec![0.0; ratio * self.modes];
        Ok(())
    }
    /// Anti-alias, derivative and interpolation latency, excluding physical travel.
    pub fn filter_delay_samples(&self) -> f64 {
        self.decimator.delay_output_frames()
            + self.pressure_derivative.as_ref().map_or(0.0, antialias::Decimator::delay_output_frames)
            + self.interpolation_delay as f64
            + self.material_path_delay as f64
    }
    /// Cold, before stepping: a differentiating receiver at one mechanics frame per output
    /// frame takes the bandlimited derivative of the summed pressure instead of every input's
    /// velocity. The derivative and the Rayleigh kernels are both linear and time-invariant, and
    /// both start from rest: D(sum_k w_k . v(t-d_k)) = sum_k w_k . D(v)(t-d_k). The same filter
    /// then runs on one channel instead of `inputs` (the live piano's r32 board: 32 x 48 tap
    /// pairs per sample). Only the rounding changes. `mode_pressures` is refused afterwards: the
    /// history then holds velocities. Returns false, changing nothing, for other receivers:
    /// their derivative is convolved into a decimation stage, which does not commute.
    pub fn differentiate_pressure(&mut self) -> Result<bool, String> {
        if !self.differentiate_velocity || self.decimator.input_frames() != 1 {return Ok(false);}
        let derivative = antialias::Decimator::new_differentiating(1, 1, self.rate as u32)?;
        debug_assert_eq!(derivative.delay_output_frames(), self.decimator.delay_output_frames());
        self.decimator = antialias::Decimator::new(1, self.modes)?;
        self.pressure_derivative = Some(derivative);
        Ok(true)
    }
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
        let mut receiver = Self { kernels: bins.into_iter().collect(), history: vec![0.0;samples],history_frames,
            taps:Vec::new(),packed_weights:Vec::new(),tap_offsets:Vec::new(),tap_history:0,
            previous_velocity:vec![0.0;modes],acceleration:vec![0.0;modes],head:0,modes,rate:f64::from(rate),
            decimator:antialias::Decimator::new(1,modes)?,differentiate_velocity:false,pressure_derivative:None,interpolation_delay:0,
            material_path_delay:0,
            position_m,delay_samples:(first,last) };
        receiver.prepare_taps();
        Ok(receiver)
    }

    /// Cold: build the transposed tap ring for the current kernels and the
    /// still-zero history. Constructors call it after their final kernel
    /// choice. Rings above MAX_TAP_RING keep the direct form.
    fn prepare_taps(&mut self) {
        let (frames,count)=(self.history_frames,self.kernels.len());
        match frames.checked_mul(count) {
            Some(ring) if ring<=MAX_TAP_RING => {
                self.taps=vec![0.0;ring];
                self.packed_weights=self.kernels.iter().flat_map(|(_,row)|row.iter().copied()).collect();
                self.tap_offsets=self.kernels.iter().enumerate()
                    .map(|(k,(delay,_))|(frames-delay)*count+k).collect();
                self.tap_history=self.history.as_ptr() as usize;
            }
            _ => {
                self.taps=Vec::new();self.packed_weights=Vec::new();self.tap_offsets=Vec::new();
                self.tap_history=0;
            }
        }
    }
    pub fn multiply_adds_per_sample(&self)->usize {
        self.kernels.len()*self.modes + self.decimator.multiplies_per_output_frame()
            + self.pressure_derivative.as_ref().map_or(0, antialias::Decimator::multiplies_per_output_frame)
    }

    /// Pressure contribution of each input mode at the most recently emitted
    /// sample. Call after `step_trace`; this only reads the delayed history and
    /// does not advance the receiver. `sum(output)` reconstructs that sample's
    /// pressure, subject to floating-point summation order.
    pub fn mode_pressures(&self, output: &mut [f64]) -> Result<(), String> {
        if output.len() != self.modes {
            return Err("modal pressure output must have one entry per input mode".into());
        }
        if self.pressure_derivative.is_some() {
            return Err("modal pressures need the per-input derivative; differentiate_pressure is on".into());
        }
        output.fill(0.0);
        let emitted_head = if self.head == 0 {self.history_frames - 1} else {self.head - 1};
        for (delay, weights) in &self.kernels {
            let frame = if emitted_head >= *delay {emitted_head - *delay}
                else {emitted_head + self.history_frames - *delay};
            let offset = frame * self.modes;
            for (i, weight) in weights.iter().enumerate() {
                output[i] += weight * self.history[offset + i];
            }
        }
        Ok(())
    }

    /// Legacy one-rate input. Multirate callers must supply the complete trace.
    pub fn step(&mut self, velocity: &[f64])->Result<f64,String> {
        self.step_trace(velocity)
    }

    /// Interleaved END velocities from EVERY accepted mechanics substep.
    /// The differentiating image passes velocity to its prepared derivative;
    /// legacy images first difference at mechanics rate. Both anti-alias before
    /// spatial propagation at audio rate and retain explicit causal delay.
    /// All candidate values are checked before either delay line is committed.
    pub fn step_trace(&mut self, velocity: &[f64])->Result<f64,String> {
        let pressure = self.preview_trace(velocity)?;
        self.accept_trace(velocity);
        Ok(pressure)
    }

    // A distributed surface may need several original bounded input bases.
    // Preview changes only pending filter buffers; all histories wait until
    // every piece and their physical pressure sum have been admitted.
    fn preview_trace(&mut self, velocity: &[f64])->Result<f64,String> {
        if velocity.len()!=self.acceleration.len() || velocity.iter().any(|v|!v.is_finite()) {
            return Err("microphone needs the complete finite substep velocity trace".into());
        }
        let ratio=self.decimator.input_frames();
        if self.differentiate_velocity {
            self.decimator.preview(velocity)?;
        } else {
            let mechanics_rate=self.rate*ratio as f64;
            for s in 0..ratio { for i in 0..self.modes {
                let previous=if s==0 {self.previous_velocity[i]} else {velocity[(s-1)*self.modes+i]};
                self.acceleration[s*self.modes+i]=(velocity[s*self.modes+i]-previous)*mechanics_rate;
            }}
            self.decimator.preview(&self.acceleration)?;
        }
        let mut pressure=0.0;
        let transposed=self.tap_history==self.history.as_ptr() as usize;
        if !transposed {
            for (delay,weights) in &self.kernels {
                let frame=if self.head>=*delay {self.head-*delay}else{self.head+self.history_frames-*delay};
                let offset=frame*self.modes;
                pressure+=kernel_tap(weights,&self.history[offset..offset+self.modes]);
            }
        } else {
            // tap_k(head-d_k), already formed when that frame was committed.
            let ring=self.taps.len();
            let base=self.head*self.kernels.len();
            for &offset in &self.tap_offsets {
                let mut index=base+offset;
                if index>=ring {index-=ring;}
                pressure+=self.taps[index];
            }
        }
        if !pressure.is_finite() {return Err("radiated pressure overflow".into());}
        let pressure=match &mut self.pressure_derivative {
            Some(derivative)=>derivative.preview(&[pressure])?[0],
            None=>pressure,
        };
        Ok(pressure)
    }

    fn accept_trace(&mut self, velocity: &[f64]) {
        let ratio = self.decimator.input_frames();
        let filtered = self.decimator.pending_output();
        let transposed = self.tap_history == self.history.as_ptr() as usize;
        let offset=self.head*self.modes;
        self.history[offset..offset+self.modes].copy_from_slice(filtered);
        if transposed {
            let count=self.kernels.len();
            let row=&mut self.taps[self.head*count..(self.head+1)*count];
            for (tap,weights) in row.iter_mut().zip(self.packed_weights.chunks_exact(self.modes)) {
                *tap=kernel_tap(weights,filtered);
            }
        }
        self.previous_velocity.copy_from_slice(&velocity[(ratio-1)*self.modes..]);
        self.decimator.commit();
        if let Some(derivative)=&mut self.pressure_derivative {derivative.commit();}
        self.head+=1;if self.head==self.history_frames{self.head=0;}
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn g3_resolved_rayleigh_medium_preserves_density_and_sound_speed_bits() {
        use fs_material::gas::{GasSpec, GasState};
        let mut gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 308.15, 90000.0).unwrap();
        gas.density = f64::from_bits(gas.density.to_bits() + 1);
        gas.sound_speed = f64::from_bits(gas.sound_speed.to_bits() + 1);
        let medium = super::RayleighMedium::from(&gas);
        assert_eq!(medium.density.to_bits(), gas.density.to_bits());
        assert_eq!(medium.sound_speed.to_bits(), gas.sound_speed.to_bits());
    }

    use super::*;
    fn medium()->RayleighMedium {RayleighMedium{density:1.2,sound_speed:320.0}}
    fn patch(x:f64,shape:f64)->SurfaceSample {
        SurfaceSample{position_m:[x,0.0,0.0],area_m2:0.1,mode_shape:vec![shape]}
    }
    #[test]
    fn modal_pressures_reconstruct_delayed_output_without_advancing_time() {
        let surface = [SurfaceSample {position_m:[0.0,0.0,0.0],area_m2:0.1,
            mode_shape:vec![1.0,-0.5]},
            SurfaceSample {position_m:[0.25,0.0,0.0],area_m2:0.1,
            mode_shape:vec![-0.25,1.0]}];
        let mut mic=BaffledPressure::new(&surface,2,32_000,[0.0,0.0,1.0],medium()).unwrap();
        let mut components=[0.0;2];
        let mut heard = false;
        assert!(mic.mode_pressures(&mut [0.0;1]).is_err());
        for sample in 0..200 {
            let velocity=[if sample<40 {0.0} else {0.2},
                if sample<55 {0.0} else {-0.1}];
            let pressure=mic.step_trace(&velocity).unwrap();
            mic.mode_pressures(&mut components).unwrap();
            let first=components;
            mic.mode_pressures(&mut components).unwrap();
            assert_eq!(components,first);
            assert!((components.iter().sum::<f64>()-pressure).abs()<1e-12);
            heard |= components.iter().any(|x| x.abs() > 0.0);
        }
        assert!(heard);
    }
    /// The bandlimited derivative of the summed pressure is the same filter as the derivative
    /// of every input (both linear, time-invariant and from rest): same delay, fewer
    /// multiplies, the same pressure up to rounding. Modal pressures are refused, and a
    /// multirate receiver (derivative inside a decimation stage) is left as it was.
    #[test]
    fn the_summed_pressure_derivative_equals_the_per_input_derivative() {
        let surface = [SurfaceSample {position_m:[0.0,0.0,0.0],area_m2:0.1,
            mode_shape:vec![1.0,-0.5,0.25]},
            SurfaceSample {position_m:[0.25,0.1,0.0],area_m2:0.1,
            mode_shape:vec![-0.25,1.0,0.5]}];
        let make=||BaffledPressure::new_multirate_differentiating(&surface,3,48_000,48_000,
            [0.1,0.0,1.0],medium()).unwrap();
        let (mut each,mut summed)=(make(),make());
        assert!(summed.differentiate_pressure().unwrap());
        assert_eq!(summed.filter_delay_samples(),each.filter_delay_samples());
        assert!(summed.multiply_adds_per_sample()<each.multiply_adds_per_sample());
        assert!(summed.mode_pressures(&mut [0.0;3]).is_err());
        let (mut error,mut signal)=(0.0_f64,0.0_f64);
        for n in 0..4000 {
            let t=n as f64/48_000.0;
            let velocity=[(t*2.0*std::f64::consts::PI*220.0).sin()*(-t*3.0).exp(),
                if n<300 {0.0} else {(t*2.0*std::f64::consts::PI*1730.0).cos()*0.2},
                (n%97) as f64*1e-3];
            let (a,b)=(each.step_trace(&velocity).unwrap(),summed.step_trace(&velocity).unwrap());
            error+=(a-b).powi(2);
            signal+=a*a;
        }
        assert!(signal>0.0 && error<=1e-24*signal,"{error} {signal}");
        // A copy keeps stepping the summed form bit for bit.
        let mut copy=summed.clone();
        for n in 0..100 {
            let velocity=[n as f64*1e-4,0.0,-1e-3];
            assert_eq!(copy.step_trace(&velocity).unwrap().to_bits(),summed.step_trace(&velocity).unwrap().to_bits());
        }
        let mut multirate=BaffledPressure::new_multirate_differentiating(&surface,3,96_000,48_000,
            [0.1,0.0,1.0],medium()).unwrap();
        let multiplies=multirate.multiply_adds_per_sample();
        assert!(!multirate.differentiate_pressure().unwrap());
        assert_eq!(multirate.multiply_adds_per_sample(),multiplies);
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
    #[test]
    fn transposed_taps_match_the_direct_kernel_sum_bit_for_bit() {
        for (modes,ratio,differentiating) in [(5,1,false),(32,1,true),(37,4,false),(3,3,true)] {
            let surface:Vec<SurfaceSample>=(0..60).map(|p| SurfaceSample{
                position_m:[0.031*p as f64-0.9,0.017*((p*7)%23) as f64-0.2,0.0],
                area_m2:0.01+0.0003*p as f64,
                mode_shape:(0..modes).map(|i|((p*modes+i) as f64*0.71).sin()).collect()}).collect();
            let build=|| {
                let mic=if differentiating {
                    BaffledPressure::new_multirate_differentiating(&surface,modes,48_000*ratio as u32,48_000,[0.3,0.4,1.0],medium())
                } else {
                    BaffledPressure::new_multirate(&surface,modes,48_000*ratio as u32,48_000,[0.3,0.4,1.0],medium())
                };
                mic.unwrap()
            };
            let mut transposed=build();
            let mut direct=build();
            assert_eq!(transposed.tap_history,transposed.history.as_ptr() as usize);
            assert_eq!(transposed.taps.len(),transposed.history_frames*transposed.kernels.len());
            assert!(transposed.kernels.len()>8);
            direct.taps=Vec::new();direct.packed_weights=Vec::new();direct.tap_offsets=Vec::new();
            direct.tap_history=0;
            // Replacing kernels after construction (the plain bandlimited
            // constructor) must fall back to the direct form, never a stale ring.
            let replaced=BaffledPressure::new_multirate_bandlimited(&surface,modes,48_000*ratio as u32,48_000,
                [0.3,0.4,1.0],medium()).unwrap();
            assert_ne!(replaced.tap_history,replaced.history.as_ptr() as usize);
            let (mut a,mut b)=(vec![0.0;modes],vec![0.0;modes]);
            for n in 0..3*transposed.history_frames {
                let trace:Vec<f64>=(0..ratio*modes).map(|i|((n*ratio*modes+i) as f64*0.37).sin()*1e-3).collect();
                let p=transposed.step_trace(&trace).unwrap();
                assert_eq!(p.to_bits(),direct.step_trace(&trace).unwrap().to_bits(),"modes {modes} ratio {ratio} n {n}");
                transposed.mode_pressures(&mut a).unwrap();direct.mode_pressures(&mut b).unwrap();
                assert_eq!(a,b);
            }
            // A copy taken mid-stream keeps its form (the transposed ring names
            // the copy's own history) and steps bit for bit as the original.
            for mut original in [transposed,direct] {
                let mut copy=original.clone();
                let form=|m:&BaffledPressure| m.tap_history==m.history.as_ptr() as usize;
                assert_eq!(form(&copy),form(&original));
                assert_eq!(form(&copy),!copy.taps.is_empty());
                for n in 0..2*copy.history_frames {
                    let trace:Vec<f64>=(0..ratio*modes).map(|i|((n*ratio*modes+i) as f64*0.53).cos()*1e-3).collect();
                    assert_eq!(copy.step_trace(&trace).unwrap().to_bits(),original.step_trace(&trace).unwrap().to_bits());
                }
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
            assert!(BaffledPressure::circular_outlet(radius,48000,coarse,medium).is_err());
        }
        let mut high=config();high.maximum_frequency_hz=4801.0;
        assert!(BaffledPressure::circular_outlet(0.007,48000,high,medium).is_err());
    }
}

#[cfg(test)]
#[path = "baffled_pressure_differentiating_tests.rs"]
mod differentiating_tests;

#[cfg(test)]
#[path = "baffled_pressure_qualification.rs"]
mod qualification;
