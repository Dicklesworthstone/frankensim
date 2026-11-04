//! Terminal radiation from the same regional TMM response as the inlet port.
//!
//! `p+ = (Zin + Zc) Qin / 2`, so the terminal flow transfer is
//! `Qout/p+ = 2 (Qout/Qin)/(Zin + Zc)`. Differentiate that transfer on the
//! sample grid before realization: its DC is exactly zero, with no guessed
//! steady flow through the viscothermal duct. The existing Rayleigh owner
//! supplies the spatial kernel; neither constitutive laws nor a radiation
//! integral or reed stepper are reproduced here.

use super::RegionalDuct;
use crate::driving_point::DrivingPointError;
use crate::pcm_wav::baffled::{BaffledPressure, CircularOutletReceiver, RayleighMedium};
use fs_duct::{
    DuctError, DuctResponse, HoleState, LossModel, OpenHoleFlow, Segment, Termination, WallModel,
};
use fs_fft::{C64 as FftC64, Fft};
use fs_material::gas::GasState;
use fs_math::{c64::C64, det};
use fs_phs::WallPin;
use fs_vfit::discretize::DelayedFilter;

/// Original transfer or shared-observer refusal, without dropping its cause.
#[derive(Clone, Debug, PartialEq)]
pub enum OutletError {
    /// Existing TMM, Fourier or FIR admission/runtime error.
    Transfer(DrivingPointError),
    /// Original diagnostic from the shared Rayleigh owner.
    Receiver(String),
}
impl core::fmt::Display for OutletError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self { Self::Transfer(e) => write!(f, "{e}"), Self::Receiver(e) => write!(f, "outlet receiver: {e}") }
    }
}
impl std::error::Error for OutletError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self { Self::Transfer(e) => Some(e), Self::Receiver(_) => None }
    }
}
impl From<DrivingPointError> for OutletError {
    fn from(error: DrivingPointError) -> Self { Self::Transfer(error) }
}

/// Largest outlet grid: the characteristic line's own ladder ceiling
/// ([`crate::driving_point::LineLadder`]'s default `max_taps`), since the reed adapter passes the
/// line's converged grid. A low-loss real bore converges above 8192 taps.
pub const MAX_OUTLET_FFT_SAMPLES: usize = 65_536;

/// Bound the explicit bulk delay plus sampled residual, before allocation.
pub const MAX_OUTLET_TRANSFER_SAMPLES: usize = 2 * MAX_OUTLET_FFT_SAMPLES;

/// Physical binding and finite-history choices of an outlet observation.
#[derive(Clone, Copy, Debug)]
pub struct OutletReceiverInfo {
    /// Homogeneous exterior used by BOTH the radiation load and receiver.
    /// Defaults to the internal terminal medium unless explicitly bound.
    pub gas: GasState,
    /// Actual final axial radius [m], never an independent observer gain.
    pub radius_m: f64,
    /// Stationary receiver and quadrature in the outlet's outward baffle frame.
    pub receiver: CircularOutletReceiver,
    /// Clock of both the source wave and pressure samples [Hz].
    pub sample_rate_hz: u32,
    /// Power-of-two TMM sampling grid; not a continuum accuracy certificate.
    pub fft_samples: usize,
    /// Whole samples of axial travel factored out before sampling the residual.
    pub axial_delay_samples: usize,
    /// Support of the original Rayleigh quadrature's physical travel delays.
    pub exterior_delay_samples: (usize, usize),
}

/// Baffled terminal pressure driven by inlet outgoing characteristic pressure.
///
/// Both stages are the existing `DelayedFilter`. Construction is cold and
/// bounded; a sample allocates nothing. This observes the flanged terminal
/// ONLY: side-hole and structural radiation and rooms are not added.
/// The homogeneous exterior may be bound independently of the internal gas.
/// No reaction force is added to the source. Its feedback must
/// use the same duct and flanged load. The sampled finite-history image is not
/// an exact broadband time-domain solution or a certified radiation match.
#[derive(Clone, Debug)]
pub struct OutletPressure {
    acceleration: DelayedFilter,
    radiation: DelayedFilter,
    info: OutletReceiverInfo,
    poisoned: bool,
}

impl OutletPressure {
    /// Original media, geometry, clocks and delay support.
    #[must_use]
    pub const fn info(&self) -> OutletReceiverInfo { self.info }

    /// Advance one source-clock sample from outgoing pressure [Pa] to receiver Pa.
    ///
    /// # Errors
    /// Nonfinite input refuses before state changes. A later arithmetic failure
    /// poisons this complete cascade, including failures after stage one commits.
    pub fn step(&mut self, outgoing_pressure_pa: f64) -> Result<f64, OutletError> {
        if self.poisoned { return Err(invalid("outlet observer has partial state after a failed sample").into()); }
        if !outgoing_pressure_pa.is_finite() { return Err(invalid("outgoing pressure must be finite").into()); }
        let result = self.acceleration.push(outgoing_pressure_pa)
            .and_then(|acceleration| self.radiation.push(acceleration));
        match result {
            Ok(pressure) => Ok(pressure),
            Err(error) => { self.poisoned = true; Err(DrivingPointError::Discrete(error).into()) }
        }
    }
}

impl RegionalDuct {
    /// Observe the actual flanged outlet in its bound homogeneous exterior gas.
    ///
    /// `fft_samples` must be 256..=[`MAX_OUTLET_FFT_SAMPLES`] and a power of two. The reed adapter
    /// obtains it from its already-admitted characteristic line, so feedback
    /// and output use the same sampling grid without a second size heuristic.
    /// Axial whole-sample delay is removed before the inverse transform and
    /// restored as exact zero taps: no circular wrap can emit before that delay.
    /// The fractional remainder and late echoes retain a finite-period sampling
    /// approximation. Increase the grid for convergence; no error enclosure is
    /// implied. The derivative is the same backward difference as the receiver.
    ///
    /// # Errors
    /// Invalid sampling/history/receiver budgets or the existing TMM, material,
    /// loss, wall, radiation or Fourier/filter owner's refusal. No fallback.
    pub fn baffled_outlet(&self, sample_rate_hz: u32, fft_samples: usize,
        receiver: CircularOutletReceiver, wall: Option<&WallPin>)
        -> Result<OutletPressure, OutletError> {
        self.baffled_outlet_with_wall(sample_rate_hz, fft_samples, receiver, WallModel::from_pin(wall))
    }

    /// The same physical outlet transfer with a complete local material wall.
    pub fn baffled_outlet_with_wall(&self, sample_rate_hz: u32, fft_samples: usize,
        receiver: CircularOutletReceiver, wall: WallModel<'_>)
        -> Result<OutletPressure, OutletError> {
        let speed = self.exterior_gas().sound_speed;
        self.baffled_outlet_delayed_at(sample_rate_hz, fft_samples, receiver, wall, speed)
    }

    /// The same transfer with the receiver's propagation delays referenced to
    /// `delay_speed_m_s` instead of the exterior's equilibrium speed `c_0`.
    ///
    /// A dispersive observer law applied after the receiver supplies its own
    /// low-frequency excess delay `r (1/c_0 - 1/c_inf)`; its receipt names the
    /// speed that makes the total low-frequency delay `r / c_0` exactly once.
    /// The Rayleigh kernel uses the speed only in its delays.
    ///
    /// # Errors
    /// As [`Self::baffled_outlet_with_wall`], and a non-finite or non-positive
    /// speed.
    pub fn baffled_outlet_delayed_at(&self, sample_rate_hz: u32, fft_samples: usize,
        receiver: CircularOutletReceiver, wall: WallModel<'_>, delay_speed_m_s: f64)
        -> Result<OutletPressure, OutletError> {
        let terminal = self.regions.last().expect("admitted nonempty regional duct");
        let radius_m = terminal.duct.segments.last().expect("admitted axial terminal").outlet_radius();
        let gas = *self.exterior_gas();
        // Geometry/medium/phase/history admission belongs to the shared receiver.
        let mic = BaffledPressure::circular_outlet(radius_m, sample_rate_hz, receiver,
            RayleighMedium { density: gas.density, sound_speed: delay_speed_m_s })
            .map_err(OutletError::Receiver)?;
        let exterior_delay_samples = mic.delay_samples;
        let (kernel, axial_delay_samples) = acceleration_kernel(sample_rate_hz, fft_samples,
            self.one_way_time_s(), self.inlet_impedance(), |omega| {
                self.response_with_wall(omega, LossModel::Bessel, Termination::FlangedOpen, wall)
            })?;
        let dt = 1.0 / f64::from(sample_rate_hz);
        let acceleration = DelayedFilter::from_impulse_response(dt, kernel)
            .map_err(DrivingPointError::Discrete)?;
        let radiation = rayleigh_acceleration_filter(mic, sample_rate_hz)?;
        Ok(OutletPressure { acceleration, radiation, poisoned: false,
            info: OutletReceiverInfo { gas, radius_m, receiver, sample_rate_hz,
                fft_samples, axial_delay_samples, exterior_delay_samples } })
    }

    /// [`Self::baffled_outlet_delayed_at`] plus one baffled circular-disk observer
    /// per open or vented tone hole in the chain, evaluated in a single frequency
    /// sweep and positioned at `[receiver.x + (L - x_hole), receiver.y, receiver.z]`.
    pub fn baffled_outlet_and_open_holes_delayed_at(
        &self,
        sample_rate_hz: u32,
        fft_samples: usize,
        receiver: CircularOutletReceiver,
        wall: WallModel<'_>,
        delay_speed_m_s: f64,
    ) -> Result<(OutletPressure, Vec<OutletPressure>), OutletError> {
        let sites = self.open_hole_sites();
        if sites.is_empty() {
            let bell = self.baffled_outlet_delayed_at(
                sample_rate_hz,
                fft_samples,
                receiver,
                wall,
                delay_speed_m_s,
            )?;
            return Ok((bell, Vec::new()));
        }
        let terminal = self.regions.last().expect("admitted nonempty regional duct");
        let radius_m = terminal
            .duct
            .segments
            .last()
            .expect("admitted axial terminal")
            .outlet_radius();
        let gas = *self.exterior_gas();
        let medium = RayleighMedium {
            density: gas.density,
            sound_speed: delay_speed_m_s,
        };
        let bell_mic = BaffledPressure::circular_outlet(radius_m, sample_rate_hz, receiver, medium)
            .map_err(OutletError::Receiver)?;
        let bell_exterior_delay = bell_mic.delay_samples;
        let total_length_m = self.total_axial_length_m();
        let mut hole_mics = Vec::with_capacity(sites.len());
        for site in &sites {
            let mut hole_receiver = receiver;
            hole_receiver.position_m[0] =
                receiver.position_m[0] + (total_length_m - site.axial_position_m);
            let mic = BaffledPressure::circular_outlet(
                site.effective_radius_m,
                sample_rate_hz,
                hole_receiver,
                medium,
            )
            .map_err(OutletError::Receiver)?;
            hole_mics.push((mic, hole_receiver));
        }
        let hole_times: Vec<f64> = sites.iter().map(|s| s.one_way_time_s).collect();
        let ((bell_kernel, bell_delay), hole_kernels) = acceleration_kernels_with_holes(
            sample_rate_hz,
            fft_samples,
            self.one_way_time_s(),
            &hole_times,
            self.inlet_impedance(),
            |omega| {
                self.response_with_wall_and_holes(
                    omega,
                    LossModel::Bessel,
                    Termination::FlangedOpen,
                    wall,
                )
            },
        )?;
        let dt = 1.0 / f64::from(sample_rate_hz);
        let bell_acc = DelayedFilter::from_impulse_response(dt, bell_kernel)
            .map_err(DrivingPointError::Discrete)?;
        let bell_rad = rayleigh_acceleration_filter(bell_mic, sample_rate_hz)?;
        let bell = OutletPressure {
            acceleration: bell_acc,
            radiation: bell_rad,
            poisoned: false,
            info: OutletReceiverInfo {
                gas,
                radius_m,
                receiver,
                sample_rate_hz,
                fft_samples,
                axial_delay_samples: bell_delay,
                exterior_delay_samples: bell_exterior_delay,
            },
        };
        let mut holes = Vec::with_capacity(sites.len());
        for ((site, (mic, hole_receiver)), (kernel, axial_delay_samples)) in
            sites.iter().zip(hole_mics).zip(hole_kernels)
        {
            let exterior_delay_samples = mic.delay_samples;
            let acceleration = DelayedFilter::from_impulse_response(dt, kernel)
                .map_err(DrivingPointError::Discrete)?;
            let radiation = rayleigh_acceleration_filter(mic, sample_rate_hz)?;
            holes.push(OutletPressure {
                acceleration,
                radiation,
                poisoned: false,
                info: OutletReceiverInfo {
                    gas,
                    radius_m: site.effective_radius_m,
                    receiver: hole_receiver,
                    sample_rate_hz,
                    fft_samples,
                    axial_delay_samples,
                    exterior_delay_samples,
                },
            });
        }
        Ok((bell, holes))
    }

    fn total_axial_length_m(&self) -> f64 {
        self.regions
            .iter()
            .flat_map(|region| &region.duct.segments)
            .map(Segment::axial_length)
            .sum()
    }

    fn open_hole_sites(&self) -> Vec<OpenHoleSite> {
        if self.continuous.is_some() {
            return Vec::new();
        }
        let mut sites = Vec::new();
        let mut axial_x = 0.0_f64;
        let mut axial_time_s = 0.0_f64;
        for region in &self.regions {
            for segment in &region.duct.segments {
                if let Segment::ToneHole {
                    hole_radius,
                    chimney_height,
                    bore_radius,
                    state,
                } = *segment
                {
                    let sigma = match state {
                        HoleState::Open => 1.0,
                        HoleState::Closed => 0.0,
                        HoleState::Vent(s) => s.clamp(0.0, 1.0),
                    };
                    if sigma > 0.0 {
                        let inner = fs_phs::side_hole_inner_length(hole_radius, bore_radius);
                        sites.push(OpenHoleSite {
                            axial_position_m: axial_x,
                            effective_radius_m: hole_radius * sigma.sqrt(),
                            one_way_time_s: axial_time_s
                                + (chimney_height + inner) / region.gas.sound_speed,
                        });
                    }
                } else {
                    let length = segment.axial_length();
                    axial_x += length;
                    axial_time_s += length / region.gas.sound_speed;
                }
            }
        }
        sites
    }
}

struct OpenHoleSite {
    axial_position_m: f64,
    effective_radius_m: f64,
    one_way_time_s: f64,
}

fn invalid(what: &'static str) -> DrivingPointError { DrivingPointError::Invalid { what } }

// Sample the derivative DIRECTLY. Sampling flow with an invented zero-frequency
// value would add a spurious finite-length step/echo after differentiation.
fn acceleration_kernel(sample_rate_hz: u32, n: usize, one_way_time_s: f64, zc: f64,
    mut response: impl FnMut(f64) -> Result<DuctResponse, DuctError>)
    -> Result<(Vec<f64>, usize), DrivingPointError> {
    if sample_rate_hz == 0 || !(256..=MAX_OUTLET_FFT_SAMPLES).contains(&n) || !n.is_power_of_two()
        || !zc.is_finite() || zc <= 0.0 || !one_way_time_s.is_finite() || one_way_time_s <= 0.0 {
        return Err(invalid("outlet transfer requires physical clock, impedance, travel time and a power-of-two grid in 256..=65536"));
    }
    let rate = f64::from(sample_rate_hz);
    let delay = (one_way_time_s * rate).floor();
    if !delay.is_finite() || delay < 1.0 || delay > (MAX_OUTLET_TRANSFER_SAMPLES - n) as f64 {
        return Err(invalid("axial delay and response exceed the outlet transfer history budget"));
    }
    let delay = delay as usize;
    let fft = Fft::new(n);
    let mut spectrum = vec![FftC64::new(0.0, 0.0); n];
    for k in 1..=n / 2 {
        let theta = core::f64::consts::TAU * k as f64 / n as f64;
        let r = response(theta * rate).map_err(DrivingPointError::Duct)?;
        let denominator = C64::new(r.impedance.re + zc, r.impedance.im);
        let flow = C64::new(r.u_mouth_over_u_in.re, r.u_mouth_over_u_in.im)
            .scale(2.0) * denominator.recip();
        // Acoustic e^{-iwt}: conjugate of the digital backward difference.
        let derivative = C64::new(rate * (1.0 - det::cos(theta)), -rate * det::sin(theta));
        let acceleration = flow * derivative;
        let phase = theta * delay as f64;
        let residual = C64::new(acceleration.re, -acceleration.im)
            * C64::new(det::cos(phase), det::sin(phase));
        if !residual.re.is_finite() || !residual.im.is_finite() {
            return Err(invalid("nonfinite terminal acceleration transfer"));
        }
        // A real sampled sequence has a real Nyquist bin. Interior frequencies
        // retain their full complex transfer; DC stays zero by differentiation.
        spectrum[k] = FftC64::new(residual.re, if k == n / 2 { 0.0 } else { residual.im });
        if k != n / 2 { spectrum[n - k] = FftC64::new(residual.re, -residual.im); }
    }
    let mut scratch = vec![FftC64::new(0.0, 0.0); n];
    fft.inverse(&mut spectrum, &mut scratch);
    let mut kernel = vec![0.0; delay];
    kernel.extend(spectrum.iter().map(|value| value.re));
    Ok((kernel, delay))
}

#[allow(clippy::type_complexity)]
fn acceleration_kernels_with_holes(
    sample_rate_hz: u32,
    n: usize,
    bell_time_s: f64,
    hole_times_s: &[f64],
    zc: f64,
    mut response: impl FnMut(f64) -> Result<(DuctResponse, Vec<OpenHoleFlow>), DuctError>,
) -> Result<((Vec<f64>, usize), Vec<(Vec<f64>, usize)>), DrivingPointError> {
    if sample_rate_hz == 0
        || !(256..=MAX_OUTLET_FFT_SAMPLES).contains(&n)
        || !n.is_power_of_two()
        || !zc.is_finite()
        || zc <= 0.0
        || !bell_time_s.is_finite()
        || bell_time_s <= 0.0
    {
        return Err(invalid(
            "outlet transfer requires physical clock, impedance, travel time and a power-of-two grid in 256..=65536",
        ));
    }
    let rate = f64::from(sample_rate_hz);
    let max_delay = (MAX_OUTLET_TRANSFER_SAMPLES - n) as f64;
    let check_delay = |time_s: f64| -> Result<usize, DrivingPointError> {
        let delay = (time_s * rate).floor();
        if !delay.is_finite() || delay < 1.0 || delay > max_delay {
            return Err(invalid(
                "axial delay and response exceed the outlet transfer history budget",
            ));
        }
        Ok(delay as usize)
    };
    let bell_delay = check_delay(bell_time_s)?;
    let mut hole_delays = Vec::with_capacity(hole_times_s.len());
    for &time_s in hole_times_s {
        hole_delays.push(check_delay(time_s)?);
    }
    let fft = Fft::new(n);
    let mut bell_spectrum = vec![FftC64::new(0.0, 0.0); n];
    let mut hole_spectra = vec![vec![FftC64::new(0.0, 0.0); n]; hole_times_s.len()];
    let store_bin = |spectrum: &mut [FftC64], k: usize, theta: f64, delay: usize, acc: C64| -> Result<(), DrivingPointError> {
        let phase = theta * delay as f64;
        let residual =
            C64::new(acc.re, -acc.im) * C64::new(det::cos(phase), det::sin(phase));
        if !residual.re.is_finite() || !residual.im.is_finite() {
            return Err(invalid("nonfinite terminal acceleration transfer"));
        }
        spectrum[k] = FftC64::new(residual.re, if k == n / 2 { 0.0 } else { residual.im });
        if k != n / 2 {
            spectrum[n - k] = FftC64::new(residual.re, -residual.im);
        }
        Ok(())
    };
    for k in 1..=n / 2 {
        let theta = core::f64::consts::TAU * k as f64 / n as f64;
        let (r, hole_flows) = response(theta * rate).map_err(DrivingPointError::Duct)?;
        if hole_flows.len() != hole_times_s.len() {
            return Err(invalid("open-hole transfer count mismatch across frequencies"));
        }
        let denominator = C64::new(r.impedance.re + zc, r.impedance.im);
        let two_over_denom = denominator.recip().scale(2.0);
        let derivative = C64::new(rate * (1.0 - det::cos(theta)), -rate * det::sin(theta));
        let bell_flow = C64::new(r.u_mouth_over_u_in.re, r.u_mouth_over_u_in.im) * two_over_denom;
        store_bin(&mut bell_spectrum, k, theta, bell_delay, bell_flow * derivative)?;
        for ((spectrum, &delay), hole) in hole_spectra
            .iter_mut()
            .zip(&hole_delays)
            .zip(&hole_flows)
        {
            let hole_flow =
                C64::new(hole.u_hole_over_u_in.re, hole.u_hole_over_u_in.im) * two_over_denom;
            store_bin(spectrum, k, theta, delay, hole_flow * derivative)?;
        }
    }
    let mut scratch = vec![FftC64::new(0.0, 0.0); n];
    fft.inverse(&mut bell_spectrum, &mut scratch);
    let mut bell_kernel = vec![0.0; bell_delay];
    bell_kernel.extend(bell_spectrum.iter().map(|value| value.re));
    let mut hole_kernels = Vec::with_capacity(hole_spectra.len());
    for (mut spectrum, delay) in hole_spectra.into_iter().zip(hole_delays) {
        fft.inverse(&mut spectrum, &mut scratch);
        let mut kernel = vec![0.0; delay];
        kernel.extend(spectrum.iter().map(|value| value.re));
        hole_kernels.push((kernel, delay));
    }
    Ok(((bell_kernel, bell_delay), hole_kernels))
}

// Identify the existing linear Rayleigh operator using a unit volume-flow step.
// Its backward difference is one pulse of amplitude `rate`, so division by that
// known excitation yields pressure per unit volume acceleration. This is a cold
// image of the shared owner, NOT another transcription of its integral/geometry.
fn rayleigh_acceleration_filter(mut mic: BaffledPressure, rate: u32)
    -> Result<DelayedFilter, OutletError> {
    let n = (mic.delay_samples.1 + 1).max(4);
    let mut kernel = Vec::with_capacity(n);
    for _ in 0..n {
        kernel.push(mic.step(&[1.0]).map_err(OutletError::Receiver)?
            / f64::from(rate));
    }
    DelayedFilter::from_impulse_response(1.0 / f64::from(rate), kernel)
        .map_err(|error| DrivingPointError::Discrete(error).into())
}

#[cfg(test)]
mod tests;
