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
use fs_duct::{DuctError, DuctResponse, LossModel, Termination, WallModel};
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
        let terminal = self.regions.last().expect("admitted nonempty regional duct");
        let radius_m = terminal.duct.segments.last().expect("admitted axial terminal").outlet_radius();
        let gas = *self.exterior_gas();
        // Geometry/medium/phase/history admission belongs to the shared receiver.
        let mic = BaffledPressure::circular_outlet(radius_m, sample_rate_hz, receiver,
            RayleighMedium { density: gas.density, sound_speed: gas.sound_speed })
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
