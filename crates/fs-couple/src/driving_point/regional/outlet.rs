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
    /// The outlet grid would replay part of the transfer one window late above
    /// [`OUTLET_REPLAY_GAIN_BOUND`] re its direct path. The regional observers
    /// double their grid first and refuse only past [`MAX_OUTLET_FFT_SAMPLES`].
    Replay {
        /// Worst |replayed| / |direct| amplitude ratio over the gated band.
        gain: f64,
        /// Frequency of that worst bin [Hz].
        frequency_hz: f64,
        /// Pre-arrival latency the observer already applies [samples].
        latency_samples: usize,
        /// Outlet grid [samples]: the largest tried by a regional observer,
        /// otherwise the requested one; a longer grid separates the replay further.
        fft_samples: usize,
    },
}
impl core::fmt::Display for OutletError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Transfer(e) => write!(f, "{e}"),
            Self::Receiver(e) => write!(f, "outlet receiver: {e}"),
            Self::Replay {
                gain,
                frequency_hz,
                latency_samples,
                fft_samples,
            } => write!(
                f,
                "outlet transfer replays {:.1} dB re its direct path at {frequency_hz:.1} Hz one window late \
                 (bound {:.0} dB, pre-arrival latency {latency_samples} samples, outlet grid \
                 {fft_samples})",
                20.0 * gain.log10(),
                20.0 * OUTLET_REPLAY_GAIN_BOUND.log10()
            ),
        }
    }
}
impl std::error::Error for OutletError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Transfer(e) => Some(e),
            Self::Receiver(_) | Self::Replay { .. } => None,
        }
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

/// Largest admitted replay: |replayed| / |direct| amplitude at every gated bin
/// (-60 dB). A numerical tolerance set by the verifier (mail 14776), not physics.
pub const OUTLET_REPLAY_GAIN_BOUND: f64 = 1.0e-3;

/// Common observer latency on an `n`-point outlet grid: `n / 4` samples.
///
/// A sampled transfer spreads its arrival over taps before the whole-sample
/// delay (fractional arrival, finite band). A raw inverse transform puts them
/// at the END of the window, and the observer then replays the sound one
/// window late, at the direct path's level near the fundamental (clarinet G3:
/// +0.1 dB at 16384 taps, i.e. a 0.341 s hold after every breath-off). Keeping
/// the last `n / 4` taps before the arrival instead delays every observer on
/// the grid by the same amount, so bell and side openings stay aligned; it is
/// the quarter-window bound of the shared observation admission (4cd00d850).
#[must_use]
pub const fn outlet_pre_arrival_latency(n: usize) -> usize {
    n / 4
}

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
    /// Power-of-two TMM sampling grid of this observer's kernel: the requested
    /// grid or, for a regional observer, the smallest doubling of it whose
    /// replay the gate admits. Not a continuum accuracy certificate.
    pub fft_samples: usize,
    /// Whole samples of axial travel factored out before sampling the residual.
    pub axial_delay_samples: usize,
    /// Support of the original Rayleigh quadrature's physical travel delays.
    pub exterior_delay_samples: (usize, usize),
    /// Output latency shared by every observer on this grid
    /// ([`outlet_pre_arrival_latency`]); the physical arrival is emitted this
    /// many samples after `axial + exterior` and nothing precedes `axial + exterior`.
    pub pre_arrival_latency_samples: usize,
    /// Worst |replayed| / |direct| amplitude ratio this kernel admitted, and its
    /// frequency [Hz] (0 when no bin lies in the gated band).
    pub replay_gain: (f64, f64),
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
    /// Shared cold realization of a physical terminal observer. The callback
    /// supplies the original owner's `Qout/p+` [m^3/(Pa s)] with exp(-i omega t),
    /// using the same terminal load and inlet wave reference as its feedback.
    /// This keeps the derivative, delay factoring and Rayleigh radiation in one
    /// place for equilibrium and molecular constitutive owners.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn from_flow_response(
        gas: GasState,
        radius_m: f64,
        sample_rate_hz: u32,
        fft_samples: usize,
        one_way_time_s: f64,
        receiver: CircularOutletReceiver,
        delay_speed_m_s: f64,
        flow: impl FnMut(f64) -> Result<C64, DrivingPointError>,
    ) -> Result<Self, OutletError> {
        let mic = BaffledPressure::circular_outlet(radius_m, sample_rate_hz, receiver,
            RayleighMedium { density: gas.density, sound_speed: delay_speed_m_s })
            .map_err(OutletError::Receiver)?;
        let exterior_delay_samples = mic.delay_samples;
        let (residual, axial_delay_samples) =
            acceleration_residual(sample_rate_hz, fft_samples, one_way_time_s, flow)?;
        let latency = outlet_pre_arrival_latency(fft_samples);
        let band = replay_band(one_way_time_s, receiver);
        let replay_gain = admit_replay(&residual, latency, sample_rate_hz, band)?;
        let kernel = realize_with_latency(axial_delay_samples, &residual, latency);
        let dt = 1.0 / f64::from(sample_rate_hz);
        let acceleration = DelayedFilter::from_impulse_response(dt, kernel)
            .map_err(DrivingPointError::Discrete)?;
        let radiation = rayleigh_acceleration_filter(mic, sample_rate_hz)?;
        Ok(Self {
            acceleration,
            radiation,
            poisoned: false,
            info: OutletReceiverInfo {
                gas,
                radius_m,
                receiver,
                sample_rate_hz,
                fft_samples,
                axial_delay_samples,
                exterior_delay_samples,
                pre_arrival_latency_samples: latency,
                replay_gain,
            },
        })
    }

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

/// Time-domain WDF 2-port nonlinear orifice shunt at the first open tone hole.
///
/// The linear characteristic line uses the zero-vortex passive TMM (`R_v = 0`).
/// At the open-hole port (`Port 2`, numerical WDF reference resistance `R_2 > 0`),
/// the incident unthrottled volume flow from `Port 1` is `u_0(t) = (s21 * p^+)(t)`
/// and the Thévenin flow including past port reflections is
/// `u_th(t) = u_0(t) + 2 h_22(t) / (R_2 (1 + r_refl))`, with instantaneous port
/// impedance `Z_eff = R_2 (1 + r_refl) / (1 - r_refl)` (both independent of the
/// numerical choice `R_2`).
/// The orifice solves the quasi-steady sharp-edged separation law
/// `Delta p_v(U_h) = K |U_eff| U_eff` with `U_eff = max(|U_h| - U_crit, 0)` and
/// `U_crit = 0` for a sharp-edged orifice (`r_c / b -> 0`, Atig et al. 2004
/// *JSV* 269:1119-1130; Dalmont et al. 2002 *JASA* 112(5):2292-2305 Eq. 18),
/// updates the port scattering state with `a_2(t) = (Delta p_v(t) - h_22(t)) / (1 + r_refl)`,
/// and returns the pressure perturbation `delta p^-(t) = (s12 * Delta p_v)(t)` to the reed.
#[derive(Clone, Debug)]
pub struct RegisterHoleShuntLine {
    s12: fs_vfit::scattering::FirScatteringPort,
    s21: fs_vfit::scattering::FirScatteringPort,
    s22: fs_vfit::scattering::FirScatteringPort,
    r2_pa_s_m3: f64,
    u_crit_m3_s: f64,
    k_quad_pa_s2_m6: f64,
}

impl RegisterHoleShuntLine {
    /// Incoming pressure perturbation at Port 1 (the reed) from prior accepted
    /// orifice pressure drops `Delta p_v`.
    #[must_use]
    pub fn s12_history(&self) -> f64 {
        self.s12.response().history
    }

    /// Return a copy of this shunt line with `k_quad_pa_s2_m6 = 0.0` (the linear
    /// vortex-free limit) for small-amplitude / passive-limit verification.
    #[must_use]
    pub fn without_vortex_loss(&self) -> Self {
        let mut copy = self.clone();
        copy.k_quad_pa_s2_m6 = 0.0;
        copy
    }

    /// Prime sample `-1` with the initial outgoing reed wave `p_plus_prev` and `a_2 = 0`.
    ///
    /// # Errors
    /// Nonfinite input or discrete filter failure.
    pub fn prime(&mut self, p_plus_prev: f64) -> Result<(), DrivingPointError> {
        self.s21
            .push(p_plus_prev)
            .map_err(|_| invalid("register-hole s21 priming failed"))?;
        self.s22
            .push(0.0)
            .map_err(|_| invalid("register-hole s22 priming failed"))?;
        self.s12
            .push(0.0)
            .map_err(|_| invalid("register-hole s12 priming failed"))?;
        Ok(())
    }

    /// Solve the memoryless open-hole orifice at the current sample and
    /// advance `s21`, `s22`, and `s12` with the accepted `p_plus`, `a_2`, and `Delta p_v`.
    ///
    /// # Errors
    /// Nonfinite state or quadratic-orifice solver failure.
    pub fn step(&mut self, p_plus: f64) -> Result<f64, DrivingPointError> {
        let s22_resp = self.s22.response();
        let r2_refl = s22_resp.direct;
        let u_0 = self.s21.response().history;
        let h22 = s22_resp.history;
        if !r2_refl.is_finite()
            || r2_refl <= -1.0
            || r2_refl >= 1.0
            || !u_0.is_finite()
            || !h22.is_finite()
        {
            return Err(invalid(
                "register-hole shunt port left the finite passive domain",
            ));
        }
        let one_minus_r = 1.0 - r2_refl;
        let one_plus_r = 1.0 + r2_refl;
        let z_eff = self.r2_pa_s_m3 * one_plus_r / one_minus_r;
        let u_th = u_0 + 2.0 * h22 / (self.r2_pa_s_m3 * one_plus_r);
        let u_excess = u_th.abs() - self.u_crit_m3_s;
        let (u_h, dp_v, a2) = if u_excess <= 0.0 || self.k_quad_pa_s2_m6 == 0.0 {
            let a2 = -h22 / one_plus_r;
            if !a2.is_finite() {
                return Err(invalid("register-hole reflected wave is not finite"));
            }
            (u_th, 0.0, a2)
        } else {
            let disc = (self.k_quad_pa_s2_m6 / z_eff) * u_excess;
            if !disc.is_finite() || disc < 0.0 {
                return Err(invalid(
                    "register-hole orifice quadratic discriminant is not finite",
                ));
            }
            let u_eff = u_excess / (0.5 + det::sqrt(0.25 + disc));
            let u_h = u_th.signum() * (self.u_crit_m3_s + u_eff);
            let dp_v = u_th.signum() * self.k_quad_pa_s2_m6 * u_eff * u_eff;
            let a2 = (dp_v - h22) / one_plus_r;
            if !a2.is_finite() || !dp_v.is_finite() {
                return Err(invalid("register-hole reflected wave is not finite"));
            }
            (u_h, dp_v, a2)
        };
        self.s21
            .push(p_plus)
            .map_err(|_| invalid("register-hole s21 step failed"))?;
        self.s22
            .push(a2)
            .map_err(|_| invalid("register-hole s22 step failed"))?;
        self.s12
            .push(dp_v)
            .map_err(|_| invalid("register-hole s12 step failed"))?;
        Ok(u_h)
    }
}

impl RegionalDuct {
    /// Observe the actual flanged outlet in its bound homogeneous exterior gas.
    ///
    /// `fft_samples` must be 256..=[`MAX_OUTLET_FFT_SAMPLES`] and a power of two. The reed adapter
    /// obtains it from its already-admitted characteristic line, so feedback
    /// and output start from the same sampling grid without a second size heuristic.
    /// Axial whole-sample delay is removed before the inverse transform and
    /// restored as exact zero taps: no circular wrap can emit before that delay.
    /// The window's last [`outlet_pre_arrival_latency`] taps (the transfer's
    /// pre-arrival spread) are emitted right after those zeros, so every
    /// observer on the grid is that many samples late. A grid on which part of
    /// the response would still replay one window late above
    /// [`OUTLET_REPLAY_GAIN_BOUND`] doubles for the observer's kernel only (the
    /// grid used is in [`OutletReceiverInfo::fft_samples`]); past
    /// [`MAX_OUTLET_FFT_SAMPLES`] it refuses.
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
        let zc = self.inlet_impedance();
        if !zc.is_finite() || zc <= 0.0 {
            return Err(invalid("outlet transfer requires finite positive inlet impedance").into());
        }
        // The response does not depend on the grid (the residual's taps agree
        // across grids); a grid shorter than it puts its later part, or its
        // pre-arrival spread, in the replayed taps. Example: the 0.5 m, 1.5 mm
        // cylinder-and-cone test duct gates at -44.4 / -49.5 / -72.1 dB on
        // 256 / 512 / 1024 points with the same tap profile, so it is observed
        // on 1024. Only the observer's kernel moves; the feedback line stays on
        // the requested grid. Past the largest grid the refusal stands.
        let mut grid = fft_samples;
        loop {
            let observed = OutletPressure::from_flow_response(
                gas,
                radius_m,
                sample_rate_hz,
                grid,
                self.one_way_time_s(),
                receiver,
                delay_speed_m_s,
                |omega| {
                    let response = self.response_with_wall(
                        omega,
                        LossModel::Bessel,
                        Termination::FlangedOpen,
                        wall,
                    )?;
                    let denominator = C64::new(response.impedance.re + zc, response.impedance.im);
                    Ok(
                        C64::new(response.u_mouth_over_u_in.re, response.u_mouth_over_u_in.im)
                            .scale(2.0)
                            * denominator.recip(),
                    )
                },
            );
            match observed {
                Err(OutletError::Replay { .. }) if grid < MAX_OUTLET_FFT_SAMPLES => grid *= 2,
                observed => return observed,
            }
        }
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
        self.baffled_outlet_and_open_holes_and_shunt_delayed_at(
            sample_rate_hz,
            fft_samples,
            receiver,
            wall,
            delay_speed_m_s,
        )
        .map(|(bell, holes, _)| (bell, holes))
    }

    /// [`Self::baffled_outlet_and_open_holes_delayed_at`] plus an optional
    /// time-domain WDF 2-port nonlinear orifice shunt line at the first open
    /// tone hole using its intrinsic acoustic wave impedance `R_2 = Z_{0,port}`.
    pub fn baffled_outlet_and_open_holes_and_shunt_delayed_at(
        &self,
        sample_rate_hz: u32,
        fft_samples: usize,
        receiver: CircularOutletReceiver,
        wall: WallModel<'_>,
        delay_speed_m_s: f64,
    ) -> Result<(OutletPressure, Vec<OutletPressure>, Option<RegisterHoleShuntLine>), OutletError> {
        self.baffled_outlet_and_open_holes_and_shunt_delayed_at_with_port_scale(
            sample_rate_hz,
            fft_samples,
            receiver,
            wall,
            delay_speed_m_s,
            1.0,
        )
    }

    /// Build the time-domain WDF 2-port nonlinear orifice shunt line at the
    /// first open tone hole (if any) on the `fft_samples` grid without building
    /// baffled spatial Rayleigh receivers (used by the legacy observer path so
    /// in-loop reed/bore state stays bitwise identical to the baffled path).
    pub fn register_hole_shunt_line(
        &self,
        sample_rate_hz: u32,
        fft_samples: usize,
        wall: WallModel<'_>,
    ) -> Result<Option<RegisterHoleShuntLine>, OutletError> {
        let sites = self.open_hole_sites();
        if sites.is_empty() {
            return Ok(None);
        }
        let hole_times: Vec<f64> = sites.iter().map(|s| s.one_way_time_s).collect();
        let first_bore_time_s = sites[0].bore_one_way_time_s;
        let (_, _, shunt_line) = acceleration_kernels_with_holes(
            sample_rate_hz,
            fft_samples,
            self.one_way_time_s(),
            &hole_times,
            first_bore_time_s,
            self.inlet_impedance(),
            1.0,
            |omega| {
                self.response_with_wall_and_holes(
                    omega,
                    LossModel::Bessel,
                    Termination::FlangedOpen,
                    wall,
                )
            },
        )?;
        Ok(shunt_line)
    }

    /// [`Self::baffled_outlet_and_open_holes_and_shunt_delayed_at`] with the
    /// numerical WDF port reference resistance scaled by `port_resistance_scale`
    /// (`R_2 = port_resistance_scale * Z_{0,port}`).
    ///
    /// Because `R_2` is purely a numerical WDF port choice, scaling `R_2` leaves
    /// the solved orifice flow `U_h(t)`, pressure drop `Delta p_v(t)`, and
    /// reflected reed wave `delta p^-(t)` invariant to solver tolerance.
    pub fn baffled_outlet_and_open_holes_and_shunt_delayed_at_with_port_scale(
        &self,
        sample_rate_hz: u32,
        fft_samples: usize,
        receiver: CircularOutletReceiver,
        wall: WallModel<'_>,
        delay_speed_m_s: f64,
        port_resistance_scale: f64,
    ) -> Result<(OutletPressure, Vec<OutletPressure>, Option<RegisterHoleShuntLine>), OutletError> {
        let sites = self.open_hole_sites();
        if sites.is_empty() {
            let bell = self.baffled_outlet_delayed_at(
                sample_rate_hz,
                fft_samples,
                receiver,
                wall,
                delay_speed_m_s,
            )?;
            return Ok((bell, Vec::new(), None));
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
        let first_bore_time_s = sites[0].bore_one_way_time_s;
        let kernels_on = |n: usize| {
            acceleration_kernels_with_holes(
                sample_rate_hz,
                n,
                self.one_way_time_s(),
                &hole_times,
                first_bore_time_s,
                self.inlet_impedance(),
                port_resistance_scale,
                |omega| {
                    self.response_with_wall_and_holes(
                        omega,
                        LossModel::Bessel,
                        Termination::FlangedOpen,
                        wall,
                    )
                },
            )
        };
        // The shunt line is in the loop: it stays on the requested grid.
        let ((mut bell_residual, mut bell_delay), mut hole_residuals, shunt_line) =
            kernels_on(fft_samples)?;
        // One grid, one latency and one band (the whole bore's) for every mouth,
        // so the summed observers keep their relative timing and are gated alike.
        // A grid that does not hold some mouth's response to the bound doubles
        // for every mouth together (observer kernels only), as in
        // `baffled_outlet_delayed_at`.
        let band = replay_band(self.one_way_time_s(), receiver);
        let mut grid = fft_samples;
        let (bell_replay, hole_replays) = loop {
            let latency = outlet_pre_arrival_latency(grid);
            let gated =
                admit_replay(&bell_residual, latency, sample_rate_hz, band).and_then(|bell| {
                    hole_residuals
                        .iter()
                        .map(|(residual, _)| admit_replay(residual, latency, sample_rate_hz, band))
                        .collect::<Result<Vec<_>, _>>()
                        .map(|holes| (bell, holes))
                });
            match gated {
                Ok(replays) => break replays,
                Err(_) if grid < MAX_OUTLET_FFT_SAMPLES => {
                    grid *= 2;
                    let ((bell, delay), holes, _) = kernels_on(grid)?;
                    bell_residual = bell;
                    bell_delay = delay;
                    hole_residuals = holes;
                }
                Err(refusal) => return Err(refusal),
            }
        };
        let latency = outlet_pre_arrival_latency(grid);
        let mut hole_kernels = Vec::with_capacity(hole_residuals.len());
        for ((residual, delay), replay) in hole_residuals.into_iter().zip(hole_replays) {
            hole_kernels.push((
                realize_with_latency(delay, &residual, latency),
                delay,
                replay,
            ));
        }
        let bell_kernel = realize_with_latency(bell_delay, &bell_residual, latency);
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
                fft_samples: grid,
                axial_delay_samples: bell_delay,
                exterior_delay_samples: bell_exterior_delay,
                pre_arrival_latency_samples: latency,
                replay_gain: bell_replay,
            },
        };
        let mut holes = Vec::with_capacity(sites.len());
        for ((site, (mic, hole_receiver)), (kernel, axial_delay_samples, replay_gain)) in
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
                    fft_samples: grid,
                    axial_delay_samples,
                    exterior_delay_samples,
                    pre_arrival_latency_samples: latency,
                    replay_gain,
                },
            });
        }
        Ok((bell, holes, shunt_line))
    }

    fn total_axial_length_m(&self) -> f64 {
        self.regions
            .iter()
            .flat_map(|region| &region.duct.segments)
            .map(Segment::axial_length)
            .sum()
    }

    fn open_hole_sites(&self) -> Vec<OpenHoleSite> {
        if let Some(field) = &self.continuous {
            return field.open_hole_sites().into_iter().map(|site| OpenHoleSite {
                axial_position_m: site.axial_position_m,
                effective_radius_m: site.effective_radius_m,
                bore_one_way_time_s: site.bore_one_way_time_s,
                one_way_time_s: site.one_way_time_s,
            }).collect();
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
                            bore_one_way_time_s: axial_time_s,
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
    bore_one_way_time_s: f64,
    one_way_time_s: f64,
}

fn invalid(what: &'static str) -> DrivingPointError { DrivingPointError::Invalid { what } }

/// Gated band [Hz]: from the quarter-wave frequency of the whole bore's
/// one-way travel time (no tone the bore sounds lies lower; nearer DC the
/// kernel's exact zero DC forces |replayed| = |direct|) to the receiver's
/// declared band top.
fn replay_band(bore_one_way_time_s: f64, receiver: CircularOutletReceiver) -> (f64, f64) {
    (0.25 / bore_one_way_time_s, receiver.maximum_frequency_hz)
}

/// Worst |replayed| / |direct| amplitude ratio over the bins of `band_hz`, and its
/// frequency, for a circular residual whose whole-sample arrival is tap 0 and
/// whose last `latency` taps will be emitted before that arrival.
///
/// Replayed = taps `[n/2, n - latency)`: emitted at least half a window after
/// the arrival, past the offsets this grid resolves, i.e. the one-window replay
/// of pre-arrival spread (or an unresolved tail). Direct = every other tap.
/// Both are exact at the grid's bins.
pub(crate) fn replay_gain(
    residual: &[f64],
    latency: usize,
    sample_rate_hz: u32,
    band_hz: (f64, f64),
) -> (f64, f64) {
    let n = residual.len();
    let replayed = |j: usize| j >= n / 2 && j < n - latency.min(n / 2);
    let mut direct: Vec<FftC64> = residual
        .iter()
        .enumerate()
        .map(|(j, &v)| FftC64::new(if replayed(j) { 0.0 } else { v }, 0.0))
        .collect();
    let mut late: Vec<FftC64> = residual
        .iter()
        .enumerate()
        .map(|(j, &v)| FftC64::new(if replayed(j) { v } else { 0.0 }, 0.0))
        .collect();
    let fft = Fft::new(n);
    let mut scratch = vec![FftC64::new(0.0, 0.0); n];
    fft.forward(&mut direct, &mut scratch);
    fft.forward(&mut late, &mut scratch);
    let bin_hz = f64::from(sample_rate_hz) / n as f64;
    let mut worst = (0.0, 0.0);
    for k in 1..=n / 2 {
        let frequency = k as f64 * bin_hz;
        if frequency < band_hz.0 || frequency > band_hz.1 {
            continue;
        }
        let d = det::hypot(direct[k].re, direct[k].im);
        let w = det::hypot(late[k].re, late[k].im);
        let gain = if d > 0.0 {
            w / d
        } else if w > 0.0 {
            f64::INFINITY
        } else {
            0.0
        };
        if gain > worst.0 {
            worst = (gain, frequency);
        }
    }
    worst
}

fn admit_replay(
    residual: &[f64],
    latency: usize,
    sample_rate_hz: u32,
    band_hz: (f64, f64),
) -> Result<(f64, f64), OutletError> {
    let (gain, frequency_hz) = replay_gain(residual, latency, sample_rate_hz, band_hz);
    if gain.is_nan() || gain > OUTLET_REPLAY_GAIN_BOUND {
        return Err(OutletError::Replay {
            gain,
            frequency_hz,
            latency_samples: latency,
            fft_samples: residual.len(),
        });
    }
    Ok((gain, frequency_hz))
}

/// `delay` exact zeros, then the window rotated so its last `latency` taps
/// precede the arrival. The length stays `delay + n`.
fn realize_with_latency(delay: usize, residual: &[f64], latency: usize) -> Vec<f64> {
    let split = residual.len() - latency;
    let mut kernel = vec![0.0; delay];
    kernel.extend_from_slice(&residual[split..]);
    kernel.extend_from_slice(&residual[..split]);
    kernel
}

// Sample the derivative DIRECTLY. Sampling flow with an invented zero-frequency
// value would add a spurious finite-length step/echo after differentiation.
// Returns the circular residual (arrival at tap 0) and the removed whole delay.
fn acceleration_residual(
    sample_rate_hz: u32,
    n: usize,
    one_way_time_s: f64,
    mut flow_response: impl FnMut(f64) -> Result<C64, DrivingPointError>,
) -> Result<(Vec<f64>, usize), DrivingPointError> {
    if sample_rate_hz == 0 || !(256..=MAX_OUTLET_FFT_SAMPLES).contains(&n) || !n.is_power_of_two()
        || !one_way_time_s.is_finite() || one_way_time_s <= 0.0 {
        return Err(invalid("outlet transfer requires physical clock, travel time and a power-of-two grid in 256..=65536"));
    }
    let rate = f64::from(sample_rate_hz);
    let delay = (one_way_time_s * rate).floor();
    if !delay.is_finite() || delay < 1.0 || delay > (MAX_OUTLET_TRANSFER_SAMPLES - n) as f64 {
        return Err(invalid("axial delay and response exceed the outlet transfer history budget"));
    }
    let delay = delay as usize;
    let mut bins = Vec::with_capacity(n / 2);
    for k in 1..=n / 2 {
        let theta = core::f64::consts::TAU * k as f64 / n as f64;
        let flow = flow_response(theta * rate)?;
        // Acoustic e^{-iwt}: conjugate of the digital backward difference.
        let derivative = C64::new(rate * (1.0 - det::cos(theta)), -rate * det::sin(theta));
        bins.push(residual_bin(theta, delay, flow * derivative)?);
    }
    Ok((circular_residual(n, bins), delay))
}

/// One positive-frequency bin of the circular residual: the conjugated
/// acceleration transfer with the whole-sample `delay` removed.
fn residual_bin(theta: f64, delay: usize, acceleration: C64) -> Result<C64, DrivingPointError> {
    let phase = theta * delay as f64;
    let residual =
        C64::new(acceleration.re, -acceleration.im) * C64::new(det::cos(phase), det::sin(phase));
    if !residual.re.is_finite() || !residual.im.is_finite() {
        return Err(invalid("nonfinite terminal acceleration transfer"));
    }
    Ok(residual)
}

/// The real circular residual from its bins k = 1..=n/2. DC stays zero by
/// differentiation; interior bins keep their full complex transfer except for
/// the smooth Nyquist rotation of [`rotate_to_real_nyquist`].
fn circular_residual(n: usize, mut bins: Vec<C64>) -> Vec<f64> {
    rotate_to_real_nyquist(&mut bins);
    let mut spectrum = vec![FftC64::new(0.0, 0.0); n];
    for (index, residual) in bins.into_iter().enumerate() {
        let k = index + 1;
        spectrum[k] = FftC64::new(residual.re, if k == n / 2 { 0.0 } else { residual.im });
        if k != n / 2 {
            spectrum[n - k] = FftC64::new(residual.re, -residual.im);
        }
    }
    let mut scratch = vec![FftC64::new(0.0, 0.0); n];
    Fft::new(n).inverse(&mut spectrum, &mut scratch);
    spectrum.iter().map(|value| value.re).collect()
}

/// Fraction of the Nyquist frequency where the Nyquist phase rotation starts.
pub const OUTLET_NYQUIST_ROTATION_START: f64 = 0.8;

/// A real sampled kernel needs a real Nyquist bin. Dropping the imaginary part
/// of that one bin is a phase jump (up to pi/2) in the periodic spectrum, and
/// the jump rings as 1/|t| on BOTH sides of the arrival. Its far pre-arrival
/// half then wraps into the replay taps at any grid: on the G3 site card the
/// bell replayed at -39.9 dB at L = n/4 and gained only ~6 dB per doubling of L,
/// whatever n. This applies the same rotation (to the nearer real half-axis)
/// smoothly over [`OUTLET_NYQUIST_ROTATION_START`] * Nyquist .. Nyquist with a
/// raised-cosine weight. Bins below the start are untouched, |transfer| is kept
/// at every bin, and the start is a fraction of the sample rate only, never a
/// receiver's admission band. A real Nyquist bin (any whole-sample delay of a
/// real transfer) is left exactly as it was.
fn rotate_to_real_nyquist(bins: &mut [C64]) {
    let Some(&nyquist) = bins.last() else {
        return;
    };
    let phase = det::atan2(nyquist.im, nyquist.re);
    let correction = if phase.abs() <= core::f64::consts::FRAC_PI_2 {
        -phase
    } else if phase > 0.0 {
        core::f64::consts::PI - phase
    } else {
        -core::f64::consts::PI - phase
    };
    if correction == 0.0 {
        return;
    }
    let half = bins.len() as f64;
    for (index, bin) in bins.iter_mut().enumerate() {
        let fraction = (index + 1) as f64 / half;
        if fraction <= OUTLET_NYQUIST_ROTATION_START {
            continue;
        }
        let x = (fraction - OUTLET_NYQUIST_ROTATION_START) / (1.0 - OUTLET_NYQUIST_ROTATION_START);
        let angle = correction * (0.5 - 0.5 * det::cos(core::f64::consts::PI * x));
        *bin = *bin * C64::new(det::cos(angle), det::sin(angle));
    }
}

#[allow(clippy::type_complexity, clippy::too_many_arguments)]
fn acceleration_kernels_with_holes(
    sample_rate_hz: u32,
    n: usize,
    bell_time_s: f64,
    hole_times_s: &[f64],
    first_bore_time_s: f64,
    zc: f64,
    port_resistance_scale: f64,
    mut response: impl FnMut(f64) -> Result<(DuctResponse, Vec<OpenHoleFlow>), DuctError>,
) -> Result<
    (
        (Vec<f64>, usize),
        Vec<(Vec<f64>, usize)>,
        Option<RegisterHoleShuntLine>,
    ),
    DrivingPointError,
> {
    if sample_rate_hz == 0
        || !(256..=MAX_OUTLET_FFT_SAMPLES).contains(&n)
        || !n.is_power_of_two()
        || !zc.is_finite()
        || zc <= 0.0
        || !port_resistance_scale.is_finite()
        || port_resistance_scale <= 0.0
        || !bell_time_s.is_finite()
        || !first_bore_time_s.is_finite()
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
    let mut bell_bins = Vec::with_capacity(n / 2);
    let mut hole_bins = vec![Vec::with_capacity(n / 2); hole_times_s.len()];
    let mut shunt_spec: Option<(f64, f64, Vec<C64>, Vec<C64>)> = None;
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
        bell_bins.push(residual_bin(theta, bell_delay, bell_flow * derivative)?);
        for ((bins, &delay), hole) in hole_bins.iter_mut().zip(&hole_delays).zip(&hole_flows) {
            let hole_flow =
                C64::new(hole.u_hole_over_u_in.re, hole.u_hole_over_u_in.im) * two_over_denom;
            bins.push(residual_bin(theta, delay, hole_flow * derivative)?);
        }
        if let Some(hole0) = hole_flows.first()
            && hole0.shunt_port_reference_pa_s_m3 > 0.0
            && hole0.orifice_quadratic_coeff_pa_s2_m6 >= 0.0
        {
            let (_, _, s12_pos, s0_pos) = shunt_spec.get_or_insert_with(|| {
                (
                    hole0.shunt_port_reference_pa_s_m3,
                    hole0.orifice_quadratic_coeff_pa_s2_m6,
                    vec![C64::ZERO; n / 2 + 1],
                    vec![C64::ZERO; n / 2 + 1],
                )
            });
            let zc_over_denom = two_over_denom.scale(0.5 * zc);
            s12_pos[k] = hole0.u_chimney_in_over_u_in * zc_over_denom;
            s0_pos[k] = hole0.shunt_port_reflection;
        }
    }
    // Circular residuals (arrival at tap 0); the caller gates and realizes them.
    let bell_kernel = circular_residual(n, bell_bins);
    let hole_kernels: Vec<(Vec<f64>, usize)> = hole_bins
        .into_iter()
        .zip(hole_delays)
        .map(|(bins, delay)| (circular_residual(n, bins), delay))
        .collect();
    let shunt_line = if let Some((z0_port, k_quad, mut s12_pos, mut s0_pos)) = shunt_spec {
        if let Ok((r0, holes0)) = response(1.0e-3)
            && let Some(h0) = holes0.first()
        {
            let denom0 = C64::new(r0.impedance.re + zc, r0.impedance.im);
            let zc_over_denom0 = denom0.recip().scale(zc);
            let s12_0 = (h0.u_chimney_in_over_u_in * zc_over_denom0).re;
            if s12_0.is_finite() {
                s12_pos[0] = C64::new(s12_0, 0.0);
            }
            if h0.shunt_port_reflection.re.is_finite() {
                s0_pos[0] = C64::new(h0.shunt_port_reflection.re.clamp(-1.0, 1.0), 0.0);
            }
        }
        let mut s12_ir = crate::driving_point::realize_reflectance_ir_from_positive_spectrum(
            &s12_pos, n, 0.0, false,
        );
        let reg_delay = ((first_bore_time_s * rate).floor().max(1.0) as usize).min(s12_ir.len());
        for slot in &mut s12_ir[..reg_delay] {
            *slot = 0.0;
        }
        let scale_21 = 2.0 / zc;
        let s21_ir: Vec<f64> = s12_ir.iter().map(|&h| h * scale_21).collect();
        let mut s0_ir = crate::driving_point::realize_reflectance_ir_from_positive_spectrum(
            &s0_pos, n, 0.0, false,
        );
        for slot in &mut s0_ir[n / 2 + 1..] {
            *slot = 0.0;
        }
        let fft_os = Fft::new(2 * n);
        let mut s22_buf = vec![FftC64::new(0.0, 0.0); 2 * n];
        let mut scratch_os = vec![FftC64::new(0.0, 0.0); 2 * n];
        let s0_direct = s0_ir[0].clamp(-0.999, 0.999);
        s0_ir[0] = s0_direct;
        for (dst, &src) in s22_buf[1..=n / 2].iter_mut().zip(&s0_ir[1..=n / 2]) {
            *dst = FftC64::new(src, 0.0);
        }
        fft_os.forward(&mut s22_buf, &mut scratch_os);
        let max_tail_mag = s22_buf
            .iter()
            .map(|c| det::hypot(c.re, c.im))
            .fold(0.0_f64, f64::max);
        let max_allowed = (1.0 - 1.0e-3) - s0_direct.abs();
        if max_tail_mag > max_allowed && max_tail_mag > 0.0 {
            let damp = max_allowed / max_tail_mag;
            for slot in &mut s0_ir[1..=n / 2] {
                *slot *= damp;
            }
        }
        let r2 = z0_port * port_resistance_scale;
        let s22_ir = if (port_resistance_scale - 1.0).abs() <= f64::EPSILON {
            s0_ir
        } else {
            let mu = (z0_port - r2) / (z0_port + r2);
            s22_buf.fill(FftC64::new(0.0, 0.0));
            for (dst, &src) in s22_buf.iter_mut().zip(&s0_ir) {
                *dst = FftC64::new(src, 0.0);
            }
            fft_os.forward(&mut s22_buf, &mut scratch_os);
            for c in &mut s22_buf {
                let s0 = C64::new(c.re, c.im);
                let num = C64::new(s0.re + mu, s0.im);
                let den = C64::new(1.0 + mu * s0.re, mu * s0.im);
                let w = num * den.recip();
                *c = FftC64::new(w.re, w.im);
            }
            fft_os.inverse(&mut s22_buf, &mut scratch_os);
            s22_buf[..n].iter().map(|c| c.re).collect()
        };
        let dt = 1.0 / rate;
        let s12 = fs_vfit::scattering::FirScatteringPort::new(dt, s12_ir)
            .map_err(|_| invalid("register-hole s12 scattering port failed"))?;
        let s21 = fs_vfit::scattering::FirScatteringPort::new(dt, s21_ir)
            .map_err(|_| invalid("register-hole s21 scattering port failed"))?;
        let s22 = fs_vfit::scattering::FirScatteringPort::new(dt, s22_ir)
            .map_err(|_| invalid("register-hole s22 scattering port failed"))?;
        if !k_quad.is_finite() || k_quad < 0.0 || !r2.is_finite() || r2 <= 0.0 {
            return Err(invalid("register-hole quadratic loss coefficient is invalid"));
        }
        Some(RegisterHoleShuntLine {
            s12,
            s21,
            s22,
            r2_pa_s_m3: r2,
            u_crit_m3_s: 0.0,
            k_quad_pa_s2_m6: k_quad,
        })
    } else {
        None
    };
    Ok(((bell_kernel, bell_delay), hole_kernels, shunt_line))
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
