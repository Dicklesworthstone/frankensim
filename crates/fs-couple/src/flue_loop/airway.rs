//! The player's airway ahead of the windway (opt-in, [`FlueVoice::with_upstream_airway`]): the duct from the windway's
//! inlet back to the lungs, as [`crate::upstream_airway`] realizes it, in place of the ideal blowing source.
//!
//! **The source.** With an airway the gesture's blowing pressure is the pressure at the airway's far end (the lungs).
//! Each sample the windway's inlet, the foot, sees the airway's port `p = P_oc - Z Q` ([`UpstreamAirway::port`]), with
//! `Q = u h W` the windway's volume flux at its bulk speed `u` over the exit `h W`.
//!
//! **The foot pressure.** The windway's own quasi-steady law closes the sample. With no windway law the jet takes the
//! inviscid `P_f = rho u^2 / 2`; with the fully developed slit (`windway slit`,
//! [`crate::bernoulli_aperture::viscous::ViscousSlit`]) `P_f = a u + rho u^2 / 2` with `a = 12 mu L / h^2`. Against the
//! port, `P_oc = (a + Z h W) u + rho u^2 / 2`, a quadratic in `u` whose one non-negative root is taken in the
//! rationalized form `u = 2 P_oc / (a' + sqrt(a'^2 + 2 rho P_oc))`, `a' = a + Z h W`, and the foot pressure is the law's
//! own `P_f = (a + rho u / 2) u`, with no cancellation. With no slit term this is `sqrt(P_f) = (-k + sqrt(k^2 + 4
//! P_oc)) / 2`, `k = Z h W sqrt(2 / rho)`. With the developing slit (`windway slit-developing`,
//! [`DevelopingSlit`], Shah 1978) the foot pressure is the plenum pressure `P_s(u)` that law holds at `u`, and `P_oc =
//! P_s(u) + Z h W u` has one root, since the right side rises strictly with `u`. It is bisected to adjacent
//! floating-point values as the law's own held state is ([`DevelopingSlit::against_series`]); with `Z = 0` it is that
//! held state's speed, bit for bit. Shah's validity is judged twice: when the windway is attached, at the held
//! pressure it is given (the file route gives it the gesture's held lung pressure), and again when the airway is
//! attached, at the steady foot pressure of that lung pressure, `P_f = P_lung - R_dc Q`
//! ([`crate::upstream_airway::AirwayChain::dc_resistance`]), which the windway's held state then holds. A card must
//! pass both, so a lung pressure past one of Shah's upper bounds refuses by that bound's name even when its foot,
//! lower by the airway's drop, is inside it. The attack and release are not covered, as without an airway.
//!
//! **What reads it.** Every per-sample reader of the blowing pressure reads the foot pressure instead: the jet's bulk
//! speed (its flux), the card's speed, the breath noise's foot pressure and the blowing work. After the sample is
//! accepted, [`UpstreamAirway::commit`] takes `Q` and sends its wave down the airway. Without an airway nothing here
//! runs: the foot pressure is the gesture's, the ideal port, bit for bit.
//!
//! **A one-way windway.** The jet's law runs one way, as the voice without an airway refuses a negative blowing
//! pressure. On a sample whose open-circuit pressure is not positive the windway passes no flow: the jet's speed and
//! flux and the foot pressure it reads are zero, and the airway is committed a zero flow, so it sees a closed end that
//! sample. Air drawn back through the windway has no law here. Those samples are counted and the most negative
//! open-circuit pressure is kept ([`FlueVoice::upstream_airway_closures`]), so a reverse drive is reported, never
//! hidden. From rest, until the lung pressure's wave reaches the windway, the open-circuit pressure is the source
//! FIR's taps ahead of the airway's travel time, zero up to its realization and of either sign.
//!
//! **The arrival.** The lung pressure's wave reaches the windway after the airway's travel time `tau = L f_s / c`,
//! rarely a whole number of samples. Before it, the band-limited source FIR carries the sinc's pre-ringing of that
//! sub-sample arrival, alternating in sign each sample, so samples before `tau` can close the windway; nothing in the
//! realization lies at negative time. The closures are kept before and from the arrival apart: before it the most
//! negative open-circuit pressure is judged against [`sampled_step_undershoot`] at the airway's own `tau` times the
//! closed-end jump the chain itself realizes ([`closed_end_jump_pa`]: the source FIR's step height at its first maximum
//! after the arrival times the drive of the largest lung pressure there, twice that pressure only where the step height
//! is one), from it every sample against the no-reverse-drive bound. A
//! later fall of the lung pressure rings the same way about its own arrival, which then lies after `tau`: a render
//! with a release names those closures.
//!
//! **The jet at each sample's speed.** A no-flow sample has no jet: nothing convects and the jet reads no retarded
//! mouth velocity ([`FlueAirwayFrame::jet_transit_s`] is infinite). An open sample's jet convects by the transit law of
//! [`FlueVoice::transit_s`] at THAT sample's speed, never at the speed the sample before left, which after a no-flow
//! sample is zero. No delay is clamped: at or above the jet's minimum speed the delay lies inside its retained history.
//!
//! **A two-way windway (opt-in, [`FlueVoice::with_two_way_windway`]; MB 18587).** By default the windway runs one way
//! (above). With the two-way windway a negative open-circuit pressure draws air BACK through the windway by the same
//! law, signed: `P_oc - Z A u = a u + rho u |u| / 2`, whose one root is `u = sign(P_oc) 2 |P_oc| / (a' + sqrt(a'^2 + 2
//! rho |P_oc|))` and whose foot pressure is `(a + rho |u| / 2) u`, odd in `P_oc` bit for bit and zero at zero; on a
//! positive open-circuit pressure it is the one-way law's own arithmetic, so those samples keep their bits. The
//! developing slit's closure is taken odd the same way: `u = sign(P_oc) U(|P_oc|)` and foot pressure `sign(P_oc)
//! P_s(|u|)`, the forward law on the magnitude. The airway
//! is committed the reverse flux, so it no longer sees a closed end. The reverse coefficients are the forward ones, a
//! DECLARED STAND-IN: the measured windway tapers (opening to exit area about 1.7, Bouterse 2010) and its exit carries
//! the block's chamfer, so drawn back it enters at the chamfered exit (an entrance loss nobody prints) and leaves into
//! the mouth at the wider opening. The quadratic coefficient lies between `rho / 2` (the exit-height head, taken) and
//! `rho / (2 x 1.7^2)` (outflow at the opening's area). A reverse sample has no jet, as a no-flow sample has none:
//! the jet reads zero speed and zero foot pressure, so nothing convects, the breath noise reads zero and the blowing
//! work takes nothing from the reverse flux. Limits: the air drawn back comes from the window side at ambient pressure
//! (the airway and the pipe still exchange no waves); the airway stays filled with the source gas, mixing has no model;
//! the reverse coefficient is the stand-in above.
//!
//! **Refused by name.** The rows that hold the jet on the HELD blowing pressure, which with an airway is no longer the
//! foot's: the unsteady windway (its channel is driven by the gesture's pressure, not by the port) and the card
//! centreline speed (held at attachment). Also the separated-jet slit law, which no flue row selects. The developing
//! windway is admitted when it is attached before the airway, which judges its validity at the steady foot pressure;
//! attached after, it is refused by name ([`DEVELOPING_WINDWAY_AFTER_AIRWAY`]).
//!
//! **Limits, not hidden.** The windway sees the airway through its quasi-steady law alone, and the jet's speed reads
//! the foot pressure and nothing else, so the bore's acoustic pressure at the mouth never reaches the airway: the
//! airway and the pipe exchange no waves. With a held lung pressure the airway settles to `P_f = P_lung - R_dc Q`
//! (both FIRs carry exact DC bins), and between holds it shapes the attack and the release. Every tube is filled with
//! the voice's own source gas. The forced jet's reference pressure (the Reynolds number of its integration) stays the
//! gesture's held pressure, and with the developing slit so does the exit profile the forced jet takes from it (its
//! core over bulk speed), both fixed when the jet is attached, before the airway. Everything [`crate::upstream_airway`]
//! states as its limits applies.

use super::{FlueError, FlueVoice};
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::developing::DevelopingSlit;
use crate::bernoulli_aperture::viscous::ViscousSlit;
use crate::driving_point::{DrivingPointError, LineLadder};
use crate::upstream_airway::{AirwayChain, UpstreamAirway, UpstreamPort};

/// The airway judges the developing windway's validity at the steady foot pressure, so the windway comes first.
pub const DEVELOPING_WINDWAY_AFTER_AIRWAY: &str =
    "the developing windway must be attached before the upstream airway, which judges its validity at the steady foot pressure";
/// The unsteady windway's channel is driven by the gesture's pressure.
pub const AIRWAY_REFUSES_UNSTEADY_WINDWAY: &str =
    "an upstream airway refuses the unsteady windway: its channel is driven by the gesture's pressure, not by the airway's port";
/// The card centreline speed is held at attachment on the held blowing pressure.
pub const AIRWAY_REFUSES_CARD_CENTRELINE: &str =
    "an upstream airway refuses the card centreline speed: it is held on the held blowing pressure, not on the airway's foot pressure";
/// No flue row selects the separated-jet slit law; the foot closure is written for the fully developed one.
pub const AIRWAY_REFUSES_SEPARATED_SLIT: &str =
    "an upstream airway takes no windway law or the fully developed slit, not the separated-jet slit";
/// The accepted samples on which the windway passed no flow (module doc, a one-way windway), in all and split at the
/// source FIR's arrival (module doc, the arrival). Sample indices count the voice's own clock from its first accepted
/// sample.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct FlueAirwayClosures {
    /// Accepted samples whose open-circuit pressure was not positive.
    pub samples: u64,
    /// Runs of consecutive such samples.
    pub runs: u64,
    /// The first such sample's index; `None` when there was none.
    pub first_sample: Option<u64>,
    /// The last such sample's index; `None` when there was none.
    pub last_sample: Option<u64>,
    /// The most negative open-circuit pressure on those samples [Pa]; zero when none was below zero.
    pub least_open_circuit_pa: f64,
    /// The source FIR's arrival `tau` [samples]: the airway's travel time from the lungs to the windway, `L f_s / c`,
    /// as the airway realizes it. A sample `n < tau` lies before the arrival.
    pub arrival_samples: f64,
    /// The largest lung pressure on the accepted samples before the arrival [Pa]: the drive of the closed-end jump
    /// ([`closed_end_jump_pa`]) the pre-arrival ringing is judged against ([`sampled_step_undershoot`]).
    pub lung_before_arrival_pa: f64,
    /// The most negative open-circuit pressure on the no-flow samples before the arrival [Pa]; zero when none.
    pub least_before_arrival_pa: f64,
    /// The no-flow samples from the arrival on (`n >= tau`).
    pub samples_from_arrival: u64,
    /// The first of them; `None` when there was none.
    pub first_sample_from_arrival: Option<u64>,
    /// Their most negative open-circuit pressure [Pa]; zero when none was below zero.
    pub least_from_arrival_pa: f64,
    /// The most negative windway flux committed to the airway [m^3/s]: zero with the one-way windway, the largest flow
    /// drawn back with the two-way windway ([`FlueVoice::with_two_way_windway`]), on whose samples counted here a
    /// negative open-circuit pressure passed reverse flow rather than none.
    pub least_flow_m3_s: f64,
}

impl FlueAirwayClosures {
    /// Count the accepted sample `index` at the lung pressure `lung_pa` if its open-circuit pressure `open_circuit_pa`
    /// passed no flow, before or from the arrival.
    fn note(&mut self, index: u64, lung_pa: f64, open_circuit_pa: f64) {
        // A sample index is far below 2^53, so the cast is exact.
        #[allow(clippy::cast_precision_loss)]
        let before = (index as f64) < self.arrival_samples;
        if before {
            self.lung_before_arrival_pa = self.lung_before_arrival_pa.max(lung_pa);
        }
        if open_circuit_pa > 0.0 {
            return;
        }
        if self.last_sample.is_none_or(|last| last + 1 != index) {
            self.runs += 1;
        }
        self.samples += 1;
        self.first_sample.get_or_insert(index);
        self.last_sample = Some(index);
        self.least_open_circuit_pa = self.least_open_circuit_pa.min(open_circuit_pa);
        if before {
            self.least_before_arrival_pa = self.least_before_arrival_pa.min(open_circuit_pa);
        } else {
            self.samples_from_arrival += 1;
            self.first_sample_from_arrival.get_or_insert(index);
            self.least_from_arrival_pa = self.least_from_arrival_pa.min(open_circuit_pa);
        }
    }
}

/// The deepest undershoot `B` before its arrival of the ideal lossless band-limited unit step delayed by
/// `arrival_samples` (`tau`) and sampled from 0: `B = -min over n < tau of sum_{k=0..=n} sinc(k - tau)`, zero when no
/// partial sum is negative (MB 18430).
///
/// A causal chain has nothing before its arrival. A band-limited realization of a delay that is not a whole number of
/// samples cannot be zero on `[0, tau)`: its samples there are the sinc's tails, alternating in sign and rising toward
/// the arrival, so a step's running sum dips below zero just before it. A realization with loss and the Nyquist rolloff
/// rings no more than this lossless brickwall one at the same `tau`, so `B` times the closed-end jump bounds the
/// open-circuit pressure before the arrival. A whole-sample delay gives zero. Over every fraction, `B` is 0.1338 at
/// `tau` near 20 samples, 0.1381 near 50, 0.1404 near 200 and 0.1411 in the limit of a long delay (all at a fraction
/// near 0.42). Under about 10 samples the step's missing negative-time part raises it (0.217 near 1.4).
#[must_use]
pub fn sampled_step_undershoot(arrival_samples: f64) -> f64 {
    let sinc = |x: f64| if x == 0.0 { 1.0 } else { (core::f64::consts::PI * x).sin() / (core::f64::consts::PI * x) };
    let (mut sum, mut least, mut n) = (0.0_f64, 0.0_f64, 0.0_f64);
    while n < arrival_samples {
        sum += sinc(n - arrival_samples);
        least = least.min(sum);
        n += 1.0;
    }
    0.0 - least
}

/// The chain's own step height at the windway (MB 18441): the running sum of `airway`'s realized source FIR at its
/// first maximum at or after the arrival `arrival_samples`, its first overshoot included; zero when the FIR ends first.
/// A chain that narrows toward the windway raises it above one.
#[must_use]
pub fn source_step_height(airway: &UpstreamAirway, arrival_samples: f64) -> f64 {
    let (mut sum, mut step, mut n) = (0.0_f64, 0.0_f64, 0.0_f64);
    let mut rising = false;
    for &g in airway.source_taps() {
        sum += g;
        if n >= arrival_samples {
            if rising && sum < step {
                break;
            }
            rising = true;
            step = sum;
        }
        n += 1.0;
    }
    step
}

/// The closed-end jump [Pa] a lung pressure `lung_pa` makes at the windway on `airway`: the chain's own step height
/// ([`source_step_height`]) times the drive `2 lung_pa / (1 - r0)` the port's open-circuit pressure carries
/// ([`UpstreamAirway::port`]). It is `2 lung_pa` only where the step height is one.
#[must_use]
pub fn closed_end_jump_pa(airway: &UpstreamAirway, arrival_samples: f64, lung_pa: f64) -> f64 {
    2.0 * lung_pa / (1.0 - airway.receipt().same_sample_reflection) * source_step_height(airway, arrival_samples)
}

/// One accepted sample at the windway's inlet with an airway attached.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FlueAirwayFrame {
    /// The gesture's blowing pressure, held at the airway's far end (the lungs) [Pa].
    pub lung_pa: f64,
    /// The Thevenin source the windway saw this sample.
    pub port: UpstreamPort,
    /// The windway's solved foot pressure [Pa]: the jet's, or negative on a reverse sample of the two-way windway,
    /// where the jet reads zero.
    pub foot_pa: f64,
    /// The windway's volume flux `u h W` committed to the airway [m^3/s], negative when drawn back.
    pub flow_m3_s: f64,
    /// The pressure the airway returned on commit [Pa]; the foot pressure up to rounding.
    pub valve_pa: f64,
    /// The jet's convective transit this sample used [s]: the law of [`FlueVoice::transit_s`] at this sample's own
    /// speed, infinite on a no-flow sample, which has no jet.
    pub jet_transit_s: f64,
}

/// The realized airway, the last accepted frame and the samples it passed no flow.
pub(super) struct FlueUpstream {
    pub(super) airway: UpstreamAirway,
    pub(super) frame: Option<FlueAirwayFrame>,
    pub(super) closures: FlueAirwayClosures,
    /// The two-way windway ([`FlueVoice::with_two_way_windway`]); `false` runs it one way.
    pub(super) two_way: bool,
}

/// This sample's port and its closure, before the sample is accepted: the windway's signed foot pressure and bulk
/// speed, negative only on a reverse sample of the two-way windway.
#[derive(Debug, Clone, Copy)]
pub(super) struct AirwayFoot {
    pub(super) port: UpstreamPort,
    pub(super) foot_pa: f64,
    pub(super) speed_m_s: f64,
}

impl AirwayFoot {
    /// The foot pressure the jet reads: the windway's, or zero on a reverse sample, which has no jet.
    pub(super) fn jet_foot_pa(&self) -> f64 {
        if self.foot_pa > 0.0 { self.foot_pa } else { 0.0 }
    }

    /// The jet's bulk speed: the windway's, or zero on a reverse sample.
    pub(super) fn jet_speed_m_s(&self) -> f64 {
        if self.speed_m_s > 0.0 { self.speed_m_s } else { 0.0 }
    }
}

pub(super) fn airway_error(error: DrivingPointError) -> FlueError {
    match error {
        DrivingPointError::Invalid { what } => FlueError::Invalid { what },
        other => FlueError::Physics(AcousticRealizeError::Nonlinear(format!("upstream airway refused: {other}"))),
    }
}

/// The slit's linear coefficient `a = 12 mu L / h^2` [Pa s/m], as the slit law forms it; zero with no windway law.
fn slit_linear(law: Option<ViscousSlit>, opening_m: f64) -> Result<f64, FlueError> {
    match law {
        None => Ok(0.0),
        Some(law) if law.is_separated_jet() => Err(FlueError::Invalid { what: AIRWAY_REFUSES_SEPARATED_SLIT }),
        Some(law) => Ok((12.0 * (law.dynamic_viscosity_pa_s() * law.length_m()) / opening_m) / opening_m),
    }
}

/// The windway's closure against `port` (module doc): the bulk speed `u >= 0` with `P_oc - Z A u = a u + rho u^2 / 2`,
/// returned as `(foot pressure, u)`, the foot pressure the law's own `(a + rho u / 2) u`. A non-positive open-circuit
/// pressure passes no flow, `(0, 0)` (a one-way windway).
///
/// # Errors
/// A non-finite open-circuit pressure, a non-finite or negative coefficient, or a non-finite result.
pub fn foot_pressure(port: UpstreamPort, linear_pa_s_m: f64, exit_area_m2: f64, density: f64)
    -> Result<(f64, f64), FlueError>
{
    closure(port, linear_pa_s_m, exit_area_m2, density, false)
}

/// The two-way windway's closure against `port` (module doc): the signed bulk speed `u` with `P_oc - Z A u = a u + rho
/// u |u| / 2`, returned as `(foot pressure, u)`, the foot pressure `(a + rho |u| / 2) u`. Odd in the open-circuit
/// pressure bit for bit, `(0, 0)` at zero, and on a positive one the same bits as [`foot_pressure`].
///
/// # Errors
/// As [`foot_pressure`].
pub fn foot_pressure_two_way(port: UpstreamPort, linear_pa_s_m: f64, exit_area_m2: f64, density: f64)
    -> Result<(f64, f64), FlueError>
{
    closure(port, linear_pa_s_m, exit_area_m2, density, true)
}

/// Both closures: one way (`two_way` false) a non-positive open-circuit pressure passes no flow; two way it draws
/// the flow back by the same law on its magnitude.
fn closure(port: UpstreamPort, linear_pa_s_m: f64, exit_area_m2: f64, density: f64, two_way: bool)
    -> Result<(f64, f64), FlueError>
{
    let p = port.open_circuit_pa;
    if !p.is_finite() {
        return Err(FlueError::Invalid { what: "the airway's open-circuit pressure is not finite" });
    }
    let a = linear_pa_s_m + port.impedance * exit_area_m2;
    if !(a.is_finite() && a >= 0.0 && density.is_finite() && density > 0.0) {
        return Err(FlueError::Invalid { what: "the windway's closure against the airway needs finite non-negative coefficients" });
    }
    if p == 0.0 || (p < 0.0 && !two_way) {
        return Ok((0.0, 0.0));
    }
    // On a positive pressure `drive` is `p` itself, so the forward arithmetic is the one-way law's, bit for bit.
    let drive = p.abs();
    let magnitude = 2.0 * drive / (a + (a * a + 2.0 * density * drive).sqrt());
    let speed = if p < 0.0 { -magnitude } else { magnitude };
    let foot = (linear_pa_s_m + 0.5 * density * magnitude) * speed;
    if !(speed.is_finite() && foot.is_finite()) {
        return Err(FlueError::Invalid { what: "the windway's closure against the airway left the finite set" });
    }
    Ok((foot, speed))
}

/// The developing slit's closure against `port` (module doc): the bulk speed `u` with `P_oc - Z A u = P_s(u)`, Shah's
/// plenum pressure at `u` ([`DevelopingSlit::against_series`]), returned as `(foot pressure, u)`, the foot pressure the
/// law's own `P_s(u)`. One way a non-positive open-circuit pressure passes no flow, `(0, 0)`; two ways the forward law
/// is taken on its magnitude and signed, odd bit for bit.
///
/// # Errors
/// A non-finite open-circuit pressure, a non-finite or negative impedance, or a non-finite result.
fn developing_closure(law: &DevelopingSlit, height_m: f64, port: UpstreamPort, exit_area_m2: f64, two_way: bool)
    -> Result<(f64, f64), FlueError>
{
    let p = port.open_circuit_pa;
    if !p.is_finite() {
        return Err(FlueError::Invalid { what: "the airway's open-circuit pressure is not finite" });
    }
    let series = port.impedance * exit_area_m2;
    if !(series.is_finite() && series >= 0.0) {
        return Err(FlueError::Invalid { what: "the windway's closure against the airway needs finite non-negative coefficients" });
    }
    if p == 0.0 || (p < 0.0 && !two_way) {
        return Ok((0.0, 0.0));
    }
    let (foot, speed) = law.against_series(height_m, p.abs(), series);
    let (foot, speed) = if p < 0.0 { (-foot, -speed) } else { (foot, speed) };
    if !(speed.is_finite() && foot.is_finite()) {
        return Err(FlueError::Invalid { what: "the windway's closure against the airway left the finite set" });
    }
    Ok((foot, speed))
}

impl FlueVoice {
    /// Put the player's airway ahead of the windway (module doc) before any accepted sample: `chain` from the windway's
    /// inlet to the lungs, realized in the voice's own source gas on its own clock with `ladder`. From then on the
    /// blowing pressure is the lung pressure and the jet reads the foot pressure each sample. An attached developing
    /// windway is judged again at the steady foot pressure of its held pressure, now the lungs' (module doc), and holds
    /// that state from then on.
    ///
    /// # Errors
    /// Late or repeated attachment, a poisoned voice, the unsteady windway or the card centreline speed (by name,
    /// [`AIRWAY_REFUSES_UNSTEADY_WINDWAY`], [`AIRWAY_REFUSES_CARD_CENTRELINE`]), the separated-jet slit, a clock that is
    /// not a whole number of hertz, a developing windway whose steady foot pressure lies outside Shah's validity (by the
    /// bound's own name), or the airway's own realization refusal.
    pub fn with_upstream_airway(mut self, chain: &AirwayChain, ladder: &LineLadder) -> Result<Self, FlueError> {
        if self.poisoned || self.samples != 0 || self.upstream_airway.is_some() {
            return Err(FlueError::Invalid { what: "an upstream airway requires one cold attachment" });
        }
        if let Some((law, held)) = self.developing_windway {
            // The held pressure is the lungs' now; the windway's steady foot is below it by R_dc Q.
            let (width, height) = (self.geometry.flue_width_m, self.geometry.flue_height_m);
            let series = chain.dc_resistance(&self.gas) * width * height;
            let (foot_pa, _) = law.against_series(height, held.pressure_pa, series);
            let foot = law.held(width, height, foot_pa).map_err(FlueError::Physics)?;
            foot.validity().map_err(|what| FlueError::Invalid { what })?;
            self.developing_windway = Some((law, foot));
        }
        if self.unsteady_windway.is_some() {
            return Err(FlueError::Invalid { what: AIRWAY_REFUSES_UNSTEADY_WINDWAY });
        }
        if self.card_centreline.is_some() {
            return Err(FlueError::Invalid { what: AIRWAY_REFUSES_CARD_CENTRELINE });
        }
        slit_linear(self.viscous_windway, self.geometry.flue_height_m)?;
        let rate = (1.0 / self.dt).round();
        if !(rate >= 1.0 && rate <= f64::from(u32::MAX)) || 1.0 / rate != self.dt {
            return Err(FlueError::Invalid { what: "an upstream airway needs a flue clock of a whole number of hertz" });
        }
        // The guard above bounds `rate` to the u32 range.
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let airway = UpstreamAirway::realize(chain, &self.gas, rate as u32, ladder).map_err(airway_error)?;
        // The source FIR's arrival, as the realization computes its transit.
        let arrival_samples = chain.path_length_m() / self.gas.sound_speed * rate;
        let closures = FlueAirwayClosures { arrival_samples, ..FlueAirwayClosures::default() };
        self.upstream_airway = Some(Box::new(FlueUpstream { airway, frame: None, closures, two_way: false }));
        Ok(self)
    }

    /// Run the attached airway's windway two ways (module doc, MB 18587): a negative open-circuit pressure draws air
    /// back by the windway's own law, signed, with the forward coefficients (a declared stand-in). Opt-in; without it
    /// the windway runs one way, bit for bit.
    ///
    /// # Errors
    /// No airway attached, a voice that has accepted a sample, a poisoned voice, or a repeated call.
    pub fn with_two_way_windway(mut self) -> Result<Self, FlueError> {
        let cold = !self.poisoned && self.samples == 0;
        let Some(up) = self.upstream_airway.as_deref_mut().filter(|up| cold && !up.two_way) else {
            return Err(FlueError::Invalid { what: "a two-way windway requires one cold call after the upstream airway" });
        };
        up.two_way = true;
        Ok(self)
    }

    /// Whether the attached airway's windway runs two ways; `false` without an airway.
    #[must_use]
    pub fn upstream_airway_two_way(&self) -> bool {
        self.upstream_airway.as_deref().is_some_and(|up| up.two_way)
    }

    /// The realized airway ahead of the windway, if any.
    #[must_use]
    pub fn upstream_airway(&self) -> Option<&UpstreamAirway> {
        self.upstream_airway.as_deref().map(|up| &up.airway)
    }

    /// The last accepted sample at the windway's inlet with an airway attached.
    #[must_use]
    pub fn upstream_airway_frame(&self) -> Option<&FlueAirwayFrame> {
        self.upstream_airway.as_deref().and_then(|up| up.frame.as_ref())
    }

    /// The accepted samples on which the windway passed no flow, and the most negative open-circuit pressure among
    /// them (module doc, a one-way windway), if an airway is attached.
    #[must_use]
    pub fn upstream_airway_closures(&self) -> Option<FlueAirwayClosures> {
        self.upstream_airway.as_deref().map(|up| up.closures)
    }

    /// This sample's port at the lung pressure `lung_pa` and the windway's closure against it; `None` without an airway.
    pub(super) fn airway_foot(&mut self, lung_pa: f64) -> Result<Option<AirwayFoot>, FlueError> {
        // Without an airway nothing is read or computed.
        if self.upstream_airway.is_none() {
            return Ok(None);
        }
        let linear = slit_linear(self.viscous_windway, self.geometry.flue_height_m)?;
        let height = self.geometry.flue_height_m;
        let area = height * self.geometry.flue_width_m;
        let rho = self.rho;
        let developing = self.developing_windway.map(|(law, _)| law);
        let Some(up) = self.upstream_airway.as_deref_mut() else { return Ok(None) };
        let port = up.airway.port(lung_pa).map_err(airway_error)?;
        let (foot_pa, speed_m_s) = match &developing {
            Some(law) => developing_closure(law, height, port, area, up.two_way)?,
            None => closure(port, linear, area, rho, up.two_way)?,
        };
        Ok(Some(AirwayFoot { port, foot_pa, speed_m_s }))
    }

    /// Close an accepted sample: commit the windway's flux at `foot`'s speed and keep the frame, with the jet's transit
    /// `jet_transit_s` [s] this sample used.
    pub(super) fn airway_commit(&mut self, lung_pa: f64, foot: AirwayFoot, jet_transit_s: f64) -> Result<(), FlueError> {
        let flow_m3_s = foot.speed_m_s * self.geometry.flue_height_m * self.geometry.flue_width_m;
        // This sample's index: the voice counts it once it is accepted, after this commit.
        let index = self.samples;
        let Some(up) = self.upstream_airway.as_deref_mut() else { return Ok(()) };
        let valve_pa = up.airway.commit(flow_m3_s).map_err(airway_error)?;
        up.frame = Some(FlueAirwayFrame { lung_pa, port: foot.port, foot_pa: foot.foot_pa, flow_m3_s, valve_pa, jet_transit_s });
        up.closures.note(index, lung_pa, foot.port.open_circuit_pa);
        up.closures.least_flow_m3_s = up.closures.least_flow_m3_s.min(flow_m3_s);
        // The jet's stored speed: zero on a reverse sample, which has no jet.
        self.jet_speed_m_s = foot.jet_speed_m_s();
        Ok(())
    }
}

#[cfg(test)]
#[allow(clippy::float_cmp)]
mod tests {
    use super::*;
    use crate::bernoulli_aperture::network::{NetworkNode, TubeNetworkSpec, TubeSection};
    use crate::flue_loop::{FlueControl, FlueGeometry, network};
    use crate::upstream_airway::AirwaySection;
    use fs_duct::LossModel;
    use fs_exec::CancelGate;
    use fs_material::gas::{GasSpec, GasState};

    // Plumbing fixtures on a 48 kHz clock, not instrument specimens: a straight pipe, a two-section duct (its first
    // section's Z_s near the windway's own dP/dQ, so the loop rings down fast) ahead of eight generations of
    // Weibel's tree, rigid.
    const RATE: u32 = 48_000;
    const HOLD_PA: f64 = 300.0;
    // The held pressure of the existing developing-windway and centreline tests, inside Shah's validity here.
    const HELD_PA: f64 = 100.0;

    fn gas() -> GasState {
        GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
    }

    fn geometry() -> FlueGeometry {
        FlueGeometry { flue_width_m: 0.010, flue_height_m: 0.001, cut_up_m: 0.004, labium_offset_m: 0.0,
            jet_angle_rad: 0.0, bore_area_m2: core::f64::consts::PI * 0.006 * 0.006, mouth_loss_rel: 0.02 }
    }

    fn voice() -> FlueVoice {
        let gas = gas();
        let source = network::FlueNetworkSpec {
            graph: TubeNetworkSpec {
                nodes: vec![NetworkNode::Inlet, NetworkNode::Termination { reflection: -0.7 }],
                sections: vec![TubeSection { nodes: [0, 1], length_m: 0.3, radius_m: 0.006, max_length_error_m: 0.005 }],
                sound_speed_m_s: gas.sound_speed, max_wave_memory_bytes: 1_000_000,
            },
            gases: vec![gas], viscothermal: vec![], minimum_jet_speed_m_s: 1.0,
            maximum_jet_mach: 0.3, max_jet_history_samples: 4096,
        };
        FlueVoice::from_network(Some(&fs_aeroac::jetcard::mint_tonal_interim_card().unwrap()),
            90.0, 1e-6, geometry(), source, RATE, 3, 0, &CancelGate::new_clock_free()).unwrap()
    }

    fn chain() -> AirwayChain {
        let duct = AirwayChain::new(vec![
            AirwaySection::from_area(0.02, 2.0e-4, 1, None).unwrap(),
            AirwaySection::from_area(0.10, 3.0e-4, 1, None).unwrap(),
        ], LossModel::Bessel).unwrap();
        duct.ahead_of(&AirwayChain::weibel_saab_2021_b(8, None).unwrap()).unwrap()
    }

    fn with_airway(v: FlueVoice) -> FlueVoice {
        v.with_upstream_airway(&chain(), &LineLadder::default()).unwrap()
    }

    fn refusal<T>(result: Result<T, FlueError>) -> &'static str {
        match result {
            Err(FlueError::Invalid { what }) => what,
            Err(other) => panic!("refused, not by name: {other:?}"),
            Ok(_) => panic!("admitted"),
        }
    }

    // One control and one sample at a time, so every sample carries its own lung pressure.
    fn step(v: &mut FlueVoice, lung_pa: f64) -> FlueAirwayFrame {
        v.apply(FlueControl::SetBlowingPressure(lung_pa)).unwrap();
        let mut out = [0.0];
        v.step_block(&mut out).unwrap();
        assert!(out[0].is_finite());
        *v.upstream_airway_frame().expect("a frame with an airway")
    }

    /// The windway law the voice holds, evaluated by its own owner: the bulk speed at a foot pressure.
    fn law_speed(v: &FlueVoice, foot_pa: f64) -> f64 {
        let g = v.geometry;
        if let Some((law, _)) = v.developing_windway() {
            return law.held(g.flue_width_m, g.flue_height_m, foot_pa).unwrap().bulk_speed_m_s;
        }
        match v.viscous_windway() {
            None => (2.0 * foot_pa / v.rho).sqrt(),
            Some(law) => law.evaluate(g.flue_width_m, g.flue_height_m, v.rho, foot_pa).unwrap().flow_m3_s
                / (g.flue_width_m * g.flue_height_m),
        }
    }

    fn closures_json(c: FlueAirwayClosures) -> String {
        let index = |i: Option<u64>| i.map_or_else(|| "null".to_owned(), |i| i.to_string());
        format!("{{\"samples\":{},\"runs\":{},\"first_sample\":{},\"last_sample\":{},\"least_open_circuit_pa\":{:e},\
            \"arrival_samples\":{:e},\"lung_before_arrival_pa\":{:e},\"least_before_arrival_pa\":{:e},\
            \"samples_from_arrival\":{},\"first_sample_from_arrival\":{},\"least_from_arrival_pa\":{:e},\"least_flow_m3_s\":{:e}}}",
            c.samples, c.runs, index(c.first_sample), index(c.last_sample), c.least_open_circuit_pa,
            c.arrival_samples, c.lung_before_arrival_pa, c.least_before_arrival_pa,
            c.samples_from_arrival, index(c.first_sample_from_arrival), c.least_from_arrival_pa, c.least_flow_m3_s)
    }

    /// The windway law the voice holds, signed, evaluated by its own owner: the bulk speed at a foot pressure of either
    /// sign (the slit law's `evaluate` takes a signed drop; the inviscid head is `sign(P) sqrt(2 |P| / rho)`).
    fn signed_law_speed(v: &FlueVoice, foot_pa: f64) -> f64 {
        let g = v.geometry;
        match v.viscous_windway() {
            None => (2.0 * foot_pa.abs() / v.rho).sqrt().copysign(foot_pa),
            Some(law) => law.evaluate(g.flue_width_m, g.flue_height_m, v.rho, foot_pa).unwrap().flow_m3_s
                / (g.flue_width_m * g.flue_height_m),
        }
    }

    /// The signed foot pressure by bisection on `P_oc - Z A u(P_f) - P_f`, decreasing in `P_f`, over `[P_oc, 0]` for a
    /// negative open-circuit pressure.
    fn direct_reverse_foot(v: &FlueVoice, port: UpstreamPort) -> f64 {
        let area = v.geometry.flue_width_m * v.geometry.flue_height_m;
        let residual = |foot: f64| port.open_circuit_pa - port.impedance * area * signed_law_speed(v, foot) - foot;
        let (mut lo, mut hi) = (port.open_circuit_pa, 0.0);
        assert!(residual(lo) >= 0.0 && residual(hi) <= 0.0);
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if mid <= lo || mid >= hi { break; }
            if residual(mid) >= 0.0 { lo = mid; } else { hi = mid; }
        }
        0.5 * (lo + hi)
    }

    fn two_way(v: FlueVoice) -> FlueVoice {
        v.with_two_way_windway().unwrap()
    }

    /// MB 18430/18441 on `v`'s closures, held pressure `held_pa`: before the arrival the least open-circuit pressure is
    /// at least `-B(tau)` times the chain's own closed-end jump, from the arrival on at least `-1e-6` of the held
    /// pressure. `B` at the airway's own `tau` ([`sampled_step_undershoot`]) is the lossless brickwall's, at most 0.1412
    /// for any `tau` of 20 samples or more. The jump is the chain's realized step height times the drive
    /// ([`closed_end_jump_pa`]); twice the lung pressure is printed beside it, being the jump only where the step is one.
    /// Returns the line's fields.
    fn judge(v: &FlueVoice, held_pa: f64, what: &str) -> String {
        let c = v.upstream_airway_closures().expect("an airway");
        let airway = v.upstream_airway().expect("an airway");
        let b = sampled_step_undershoot(c.arrival_samples);
        let step = source_step_height(airway, c.arrival_samples);
        let jump = closed_end_jump_pa(airway, c.arrival_samples, c.lung_before_arrival_pa);
        let fields = format!("\"undershoot_b\":{b:e},\"step_height\":{step:e},\"jump_chain_pa\":{jump:e},\"jump_2_lung_pa\":{:e},\"no_flow\":{}",
            2.0 * c.lung_before_arrival_pa, closures_json(c));
        assert!(c.arrival_samples >= 20.0 && b <= 0.1412, "{what}: B {b:e} at tau {:e}", c.arrival_samples);
        assert!(c.least_before_arrival_pa >= -b * jump, "{what}: before the arrival, past B of the chain's jump: {fields}");
        assert!(c.least_from_arrival_pa >= -1.0e-6 * held_pa, "{what}: a reverse drive from the arrival on: {fields}");
        fields
    }

    // The closures' line fields without judging them, for a print ahead of the asserts.
    fn judged_fields(v: &FlueVoice) -> String {
        let c = v.upstream_airway_closures().expect("an airway");
        let airway = v.upstream_airway().expect("an airway");
        format!("\"undershoot_b\":{:e},\"step_height\":{:e},\"jump_chain_pa\":{:e},\"jump_2_lung_pa\":{:e},\"no_flow\":{}",
            sampled_step_undershoot(c.arrival_samples), source_step_height(airway, c.arrival_samples),
            closed_end_jump_pa(airway, c.arrival_samples, c.lung_before_arrival_pa), 2.0 * c.lung_before_arrival_pa,
            closures_json(c))
    }

    /// A realized FIR against its physical arrival `arrival` [samples]: the largest tap, then the energy fraction and
    /// the most negative running sum of the taps before the arrival, and the same over the last quarter of the
    /// circular buffer, where a wrapped precursor would sit. A causal realization has neither.
    fn ir_json(taps: &[f64], arrival: f64) -> String {
        let total: f64 = taps.iter().map(|h| h * h).sum();
        let part = |taps: &[f64]| {
            let (mut sum, mut least, mut at) = (0.0_f64, 0.0_f64, 0_usize);
            for (k, h) in taps.iter().enumerate() {
                sum += h;
                if sum < least { least = sum; at = k; }
            }
            (taps.iter().map(|h| h * h).sum::<f64>() / total, least, at)
        };
        let before = (arrival.floor() as usize).min(taps.len());
        let tail = 3 * (taps.len() / 4);
        let (pre_energy, pre_least, pre_at) = part(&taps[..before]);
        let (tail_energy, tail_least, tail_at) = part(&taps[tail..]);
        let peak = (0..taps.len()).fold(0, |best, k| if taps[k].abs() > taps[best].abs() { k } else { best });
        format!("{{\"taps\":{},\"arrival\":{arrival:e},\"peak_sample\":{peak},\"pre_arrival\":{{\"energy\":{pre_energy:e},\
            \"least_cumsum\":{pre_least:e},\"at\":{pre_at}}},\"last_quarter\":{{\"from\":{tail},\"energy\":{tail_energy:e},\
            \"least_cumsum\":{tail_least:e},\"at\":{}}}}}", taps.len(), tail + tail_at)
    }

    /// The windways the closure tests run: no law, the fully developed slit, and Shah's developing slit held at the
    /// tests' 300 Pa (inside its validity on this slit: Chen's shift 0.009 of the drop, `x+` 0.0035, `Re_delta*` 199).
    const WINDWAYS: [&str; 3] = ["none", "slit", "developing"];

    fn windway_voice(windway: &str) -> FlueVoice {
        match windway {
            "slit" => voice().with_viscous_windway(0.06).unwrap(),
            "developing" => voice().with_developing_windway(0.015, HOLD_PA).unwrap(),
            _ => voice(),
        }
    }

    /// The foot pressure by bisection on `P_oc - Z A u(P_f) - P_f`, decreasing in `P_f` over `[0, P_oc]`.
    fn direct_foot(v: &FlueVoice, port: UpstreamPort) -> f64 {
        let area = v.geometry.flue_width_m * v.geometry.flue_height_m;
        let residual = |foot: f64| port.open_circuit_pa - port.impedance * area * law_speed(v, foot) - foot;
        let (mut lo, mut hi) = (0.0, port.open_circuit_pa);
        assert!(residual(lo) >= 0.0 && residual(hi) <= 0.0);
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if mid <= lo || mid >= hi { break; }
            if residual(mid) >= 0.0 { lo = mid; } else { hi = mid; }
        }
        0.5 * (lo + hi)
    }

    /// The closed form against a direct solve on unsteady samples: a 10 ms rise from rest and 20 ms held, one sample
    /// per step, with no windway law, the fully developed slit and the developing slit. On every sample the foot pressure matches the
    /// bisection root of the same port to 1e-11 of the open-circuit pressure, on every sample that passes flow the
    /// airway's own returned pressure matches it to 1e-9 (a no-flow sample sees a closed end, so the airway holds its
    /// open-circuit pressure there while the foot is zero), and the sample is counted unsteady where the foot pressure
    /// moved more than 1e-6 of the held pressure from the sample before. Hundreds of them are. Samples that passed no
    /// flow (a non-positive open-circuit pressure, the one-way windway) are printed and judged (`judge`): before the
    /// arrival within the band-limited step's own undershoot of the chain's jump, from it none below -1e-6 of the held
    /// pressure, no reverse drive. Every windway prints before any is asserted.
    #[test]
    fn ua_flue_001_the_closed_form_foot_pressure_is_the_direct_solve_on_unsteady_samples() {
        let mut results = Vec::new();
        for windway in WINDWAYS {
            let mut v = with_airway(windway_voice(windway));
            let lungs = (0..480).map(|n| HOLD_PA * f64::from(n) / 480.0).chain(core::iter::repeat_n(HOLD_PA, 960));
            let (mut worst, mut worst_valve, mut unsteady, mut previous) = (0.0_f64, 0.0_f64, 0usize, 0.0_f64);
            for lung in lungs {
                let frame = step(&mut v, lung);
                let scale = frame.port.open_circuit_pa.max(1.0);
                let direct = if frame.port.open_circuit_pa > 0.0 { direct_foot(&v, frame.port) } else { 0.0 };
                worst = worst.max((frame.foot_pa - direct).abs() / scale);
                if frame.port.open_circuit_pa > 0.0 {
                    worst_valve = worst_valve.max((frame.valve_pa - frame.foot_pa).abs() / scale);
                }
                if (frame.foot_pa - previous).abs() > 1.0e-6 * HOLD_PA { unsteady += 1; }
                previous = frame.foot_pa;
                assert!(frame.foot_pa >= 0.0 && frame.foot_pa <= frame.port.open_circuit_pa.max(0.0));
                let flux_speed = frame.flow_m3_s / (v.geometry.flue_width_m * v.geometry.flue_height_m);
                assert!((v.jet_speed_m_s() - flux_speed).abs() <= 1e-14 * flux_speed, "the committed flux is the jet's");
            }
            let closed = v.upstream_airway_closures().expect("an airway");
            println!("{{\"suite\":\"fs-couple\",\"case\":\"ua-flue-001\",\"windway\":\"{windway}\",\"worst_direct\":{worst:e},\"worst_valve\":{worst_valve:e},\"unsteady_samples\":{unsteady},{}}}",
                judged_fields(&v));
            results.push((windway, v, closed, worst, worst_valve, unsteady));
        }
        for (windway, v, closed, worst, worst_valve, unsteady) in results {
            // The first sample's lung pressure is zero, so it passes no flow.
            assert_eq!(closed.first_sample, Some(0), "windway {windway}");
            judge(&v, HOLD_PA, &format!("windway {windway}"));
            assert!(worst < 1.0e-11, "windway {windway}: closed form against the direct solve {worst:e}");
            assert!(worst_valve < 1.0e-9, "windway {windway}: the airway's pressure against the foot {worst_valve:e}");
            assert!(unsteady >= 100, "windway {windway}: only {unsteady} unsteady samples");
        }
    }

    /// A held lung pressure settles to the Poiseuille drop: after three lengths of both FIRs held at 300 Pa, `P_lung -
    /// P_f = R_dc Q` to 1e-6 of the drop, with no windway law, the slit and the developing slit, and the foot pressure
    /// moves less than 1e-9 of the lung pressure over the last length. The developing slit's held state, judged at
    /// attachment, is the settled foot's to 1e-6. The step from rest also prints where its no-flow samples sit
    /// against the realized airway (MB 18339): both FIRs against their physical arrivals (`ir_json`; the source
    /// transfer arrives after one transit, the reflectance after two), then, for each sample up to two transits and
    /// 32 more, the running sum of the source taps (the step response), the open-circuit pressure and whether the
    /// sample passed no flow. Every windway prints before any is asserted.
    #[test]
    fn ua_flue_002_a_held_lung_pressure_settles_to_the_poiseuille_drop() {
        let gas = gas();
        let r_dc = chain().dc_resistance(&gas);
        let transit = chain().path_length_m() / gas.sound_speed * f64::from(RATE);
        let onset = (2.0 * transit).ceil() as usize + 32;
        let mut results = Vec::new();
        for windway in WINDWAYS {
            let mut v = with_airway(windway_voice(windway));
            let airway = v.upstream_airway().unwrap();
            let length = airway.reflectance_taps().len() + airway.source_taps().len() + 8;
            let (source, reflect) = (airway.source_taps().to_vec(), airway.reflectance_taps().to_vec());
            let mut frame = step(&mut v, HOLD_PA);
            let mut onset_pa = vec![frame.port.open_circuit_pa];
            for _ in 1..2 * length {
                frame = step(&mut v, HOLD_PA);
                if onset_pa.len() < onset { onset_pa.push(frame.port.open_circuit_pa); }
            }
            let mut running = 0.0_f64;
            let rows: Vec<String> = onset_pa.iter().zip(&source).enumerate().map(|(n, (p, g))| {
                running += g;
                format!("[{n},{running:e},{p:e},{}]", u8::from(*p <= 0.0))
            }).collect();
            println!("{{\"suite\":\"fs-couple\",\"case\":\"ua-flue-002-onset\",\"windway\":\"{windway}\",\"transit\":{transit:e},\"source_ir\":{},\"reflectance_ir\":{},\"step_n_cumsum_poc_noflow\":[{}]}}",
                ir_json(&source, transit), ir_json(&reflect, 2.0 * transit), rows.join(","));
            let before = frame.foot_pa;
            for _ in 0..length { frame = step(&mut v, HOLD_PA); }
            let drop = HOLD_PA - frame.foot_pa;
            let miss = (drop - r_dc * frame.flow_m3_s).abs() / (r_dc * frame.flow_m3_s);
            let drift = (frame.foot_pa - before).abs() / HOLD_PA;
            // The developing slit's held state against the settled foot; zero for the other windways.
            let held_miss = v.developing_windway().map_or(0.0, |(_, held)| (held.pressure_pa - frame.foot_pa).abs() / frame.foot_pa);
            let closed = v.upstream_airway_closures().expect("an airway");
            println!("{{\"suite\":\"fs-couple\",\"case\":\"ua-flue-002\",\"windway\":\"{windway}\",\"r_dc\":{r_dc:e},\"flow\":{:e},\"drop_pa\":{drop:e},\"miss\":{miss:e},\"drift\":{drift:e},\"held_miss\":{held_miss:e},\"length\":{length},{}}}",
                frame.flow_m3_s, judged_fields(&v));
            results.push((windway, v, closed, drop, miss, drift, held_miss));
        }
        for (windway, v, closed, drop, miss, drift, held_miss) in results {
            // The record's arrival is the source FIR's own transit.
            assert!((closed.arrival_samples - transit).abs() <= 1e-12 * transit, "windway {windway}: arrival {:e}", closed.arrival_samples);
            judge(&v, HOLD_PA, &format!("windway {windway}"));
            assert!(drop > 0.0 && miss < 1.0e-6, "windway {windway}: P_lung - P_f against R_dc Q: {miss:e}");
            assert!(drift < 1.0e-9, "windway {windway}: the foot pressure still moves: {drift:e}");
            assert!(held_miss < 1.0e-6, "windway {windway}: the held state against the settled foot: {held_miss:e}");
        }
    }

    /// The two rows that hold the jet on the held blowing pressure refuse the airway by name, in either order. The
    /// developing windway is admitted before the airway, which judges it at the steady foot, and refused by name after.
    #[test]
    fn ua_flue_003_the_held_pressure_rows_refuse_by_name_and_the_developing_windway_comes_first() {
        let ladder = LineLadder::default();
        assert!(voice().with_developing_windway(0.015, HELD_PA).unwrap().with_upstream_airway(&chain(), &ladder).is_ok());
        assert_eq!(refusal(with_airway(voice()).with_developing_windway(0.015, HELD_PA)), DEVELOPING_WINDWAY_AFTER_AIRWAY);
        assert_eq!(refusal(voice().with_unsteady_windway(0.015, 8, 8).unwrap().with_upstream_airway(&chain(), &ladder)),
            AIRWAY_REFUSES_UNSTEADY_WINDWAY);
        assert_eq!(refusal(with_airway(voice()).with_unsteady_windway(0.015, 8, 8)), AIRWAY_REFUSES_UNSTEADY_WINDWAY);
        assert_eq!(refusal(voice().with_card_centreline_speed(HELD_PA).unwrap().with_upstream_airway(&chain(), &ladder)),
            AIRWAY_REFUSES_CARD_CENTRELINE);
        assert_eq!(refusal(with_airway(voice()).with_card_centreline_speed(HELD_PA)), AIRWAY_REFUSES_CARD_CENTRELINE);
        // And the fully developed slit is admitted either side of the airway's closure.
        assert!(voice().with_viscous_windway(0.06).unwrap().with_upstream_airway(&chain(), &ladder).is_ok());
    }

    /// Without an airway the voice keeps no frame, counts no closures and takes the gesture's pressure; with one, a
    /// non-positive open-circuit pressure passes no flow, a positive one does; a non-finite one refuses.
    #[test]
    fn ua_flue_004_no_airway_keeps_no_frame_and_the_windway_runs_one_way() {
        let mut v = voice();
        v.apply(FlueControl::SetBlowingPressure(HOLD_PA)).unwrap();
        let mut out = [0.0; 64];
        v.step_block(&mut out).unwrap();
        assert!(v.upstream_airway().is_none() && v.upstream_airway_frame().is_none());
        assert_eq!(v.upstream_airway_closures(), None);
        assert_eq!(v.jet_speed_m_s(), (2.0 * HOLD_PA / v.rho).sqrt());
        for p in [-1.0e-3, -0.0, 0.0] {
            let port = UpstreamPort { open_circuit_pa: p, impedance: 1.0e6 };
            assert_eq!(foot_pressure(port, 0.0, 1.0e-5, 1.2).unwrap(), (0.0, 0.0), "{p}");
        }
        let (foot, speed) = foot_pressure(UpstreamPort { open_circuit_pa: 1.0e-3, impedance: 1.0e6 }, 0.0, 1.0e-5, 1.2).unwrap();
        assert!(foot > 0.0 && speed > 0.0 && foot < 1.0e-3);
        assert!(foot_pressure(UpstreamPort { open_circuit_pa: f64::NAN, impedance: 1.0e6 }, 0.0, 1.0e-5, 1.2).is_err());
        // A fresh airway counts nothing and holds only its arrival.
        let c = with_airway(voice()).upstream_airway_closures().expect("an airway");
        assert!(c.arrival_samples > 0.0, "arrival {:e}", c.arrival_samples);
        assert_eq!(c, FlueAirwayClosures { arrival_samples: c.arrival_samples, ..FlueAirwayClosures::default() });
    }

    /// The closures count samples, runs of consecutive samples, the first and last index and the most negative
    /// open-circuit pressure; a positive one is not counted. They split at the arrival: an index below it counts before,
    /// with the largest lung pressure there, any other from it. From rest at zero lung pressure every sample closes, on
    /// the voice's own clock from index 0.
    #[test]
    fn ua_flue_005_closures_count_runs_on_the_voice_clock() {
        let mut c = FlueAirwayClosures { arrival_samples: 2.5, ..FlueAirwayClosures::default() };
        for (index, lung, p) in [(0, 0.0, 0.0), (1, 10.0, -1.0e-3), (2, 20.0, 5.0), (3, 30.0, -2.0e-3), (7, 70.0, -0.0), (8, 80.0, 1.0e-9)] {
            c.note(index, lung, p);
        }
        assert_eq!((c.samples, c.runs, c.first_sample, c.last_sample), (4, 3, Some(0), Some(7)));
        assert_eq!(c.least_open_circuit_pa, -2.0e-3);
        // Indexes 0 to 2 lie before the arrival at 2.5 (index 2 passed flow but carries its lung pressure); 3 and 7 from it.
        assert_eq!((c.lung_before_arrival_pa, c.least_before_arrival_pa), (20.0, -1.0e-3));
        assert_eq!((c.samples_from_arrival, c.first_sample_from_arrival, c.least_from_arrival_pa), (2, Some(3), -2.0e-3));
        let mut v = with_airway(voice());
        for _ in 0..3 {
            let frame = step(&mut v, 0.0);
            assert_eq!((frame.foot_pa, frame.flow_m3_s), (0.0, 0.0));
        }
        let c = v.upstream_airway_closures().unwrap();
        assert_eq!((c.samples, c.runs, c.first_sample, c.last_sample, c.least_open_circuit_pa),
            (3, 1, Some(0), Some(2), 0.0));
        assert_eq!((c.lung_before_arrival_pa, c.least_before_arrival_pa, c.samples_from_arrival), (0.0, 0.0, 0));
    }

    // One held sample: its output's bits and its frame.
    fn held(v: &mut FlueVoice) -> (u64, FlueAirwayFrame) {
        v.apply(FlueControl::SetBlowingPressure(HOLD_PA)).unwrap();
        let mut out = [0.0];
        v.step_block(&mut out).unwrap();
        assert!(out[0].is_finite());
        (out[0].to_bits(), *v.upstream_airway_frame().expect("a frame with an airway"))
    }

    /// The fixture's jet is the island chain, whose convection reads the retarded mouth velocity, and its card speed is
    /// the bulk speed: the transit at speed `u` is `W / (c_p/U u)`.
    fn island_transit(v: &FlueVoice, u: f64) -> f64 {
        assert!(v.forced.is_none() && v.nolle.is_none() && v.auvray.is_none() && v.card_centreline.is_none());
        v.geometry.cut_up_m / (v.island.convection_ratio * u)
    }

    /// An open sample after a no-flow sample convects at its OWN speed (MB 18359/18363). A no-flow sample commits the
    /// jet speed zero. Read from that stored zero, the next open sample's transit was unbounded, past the jet's retained
    /// history, and once the jet had passed its minimum speed the voice refused ("jet convection exceeds its finite
    /// retained history"). Planted, with no windway law and with the slit: 480 held samples (the jet passes its minimum
    /// speed), then the stored speed a no-flow sample leaves (zero), then 64 held samples. A: each renders. B: each
    /// one's transit is `W / (c_p/U u)` of its own committed speed to 1e-12 and inside the retained history, and a twin
    /// voice that kept its true stored speed gives the same output bits and frames, so no reader of the stored speed is
    /// left in an airway sample.
    #[test]
    fn ua_flue_006_an_open_sample_after_a_no_flow_sample_convects_at_its_own_speed() {
        for slit in [false, true] {
            let fresh = || if slit { with_airway(voice().with_viscous_windway(0.06).unwrap()) } else { with_airway(voice()) };
            let (mut planted, mut twin) = (fresh(), fresh());
            for _ in 0..480 {
                assert_eq!(held(&mut planted), held(&mut twin), "slit {slit}: the twins differ before the plant");
            }
            assert!(planted.reached_minimum_speed, "slit {slit}: the jet never passed its minimum speed");
            planted.jet_speed_m_s = 0.0;
            let (dt, ring) = (planted.dt, (planted.history.len() - 1) as f64);
            let area = planted.geometry.flue_width_m * planted.geometry.flue_height_m;
            for n in 0..64 {
                let (bits, frame) = held(&mut planted);
                let u = frame.flow_m3_s / area;
                assert!(u >= planted.minimum_jet_speed_m_s(), "slit {slit}, sample {n}: speed {u:e}");
                let want = island_transit(&planted, u);
                assert!((frame.jet_transit_s - want).abs() <= 1e-12 * want,
                    "slit {slit}, sample {n}: transit {:e} against W/(c_p/U u) {want:e}", frame.jet_transit_s);
                assert!(frame.jet_transit_s / dt <= ring, "slit {slit}, sample {n}: past the retained history");
                assert_eq!((bits, frame), held(&mut twin), "slit {slit}, sample {n}: the stored speed reached the sample");
            }
        }
    }

    /// A no-flow onset renders (MB 18363): the held step from rest of `ua_flue_002`, whose first samples pass no flow
    /// on this realization. C: every sample renders, every no-flow sample has no jet (no flux, an infinite transit), and
    /// every open sample convects at its own speed to 1e-12. C needs its no-flow samples: with none to test it fails.
    #[test]
    fn ua_flue_007_a_no_flow_onset_renders() {
        for slit in [false, true] {
            let mut v = if slit { with_airway(voice().with_viscous_windway(0.06).unwrap()) } else { with_airway(voice()) };
            let area = v.geometry.flue_width_m * v.geometry.flue_height_m;
            let mut no_flow = 0usize;
            for n in 0..2048 {
                let (_, frame) = held(&mut v);
                if frame.port.open_circuit_pa <= 0.0 {
                    no_flow += 1;
                    assert_eq!((frame.foot_pa, frame.flow_m3_s, frame.jet_transit_s), (0.0, 0.0, f64::INFINITY),
                        "slit {slit}, sample {n}: a no-flow sample with a jet");
                } else {
                    let u = frame.flow_m3_s / area;
                    let want = island_transit(&v, u);
                    assert!(u > 0.0 && (frame.jet_transit_s - want).abs() <= 1e-12 * want,
                        "slit {slit}, sample {n}: transit {:e} against W/(c_p/U u) {want:e}", frame.jet_transit_s);
                }
            }
            println!("{{\"suite\":\"fs-couple\",\"case\":\"ua-flue-007\",\"slit\":{slit},\"no_flow_samples\":{no_flow},{}}}",
                judged_fields(&v));
            assert!(no_flow > 0, "slit {slit}: no no-flow sample to test");
        }
    }

    /// The undershoot `B` (MB 18430), planted: zero at a whole-sample delay; at the fixture's own arrival
    /// 49.57553507152286 samples 0.1284407524 to 1e-9; its largest over 200 fractions near 50 samples inside
    /// (0.1380, 0.1382); at most 0.1412 near 1000 samples. And the chain's own step height on the fixture is the
    /// running sum of its source taps at the first maximum from the arrival on, at least one (the chain narrows).
    #[test]
    fn ua_flue_008_the_band_limited_step_undershoot() {
        for tau in [20.0, 50.0, 200.0] {
            assert!(sampled_step_undershoot(tau) < 1e-12, "tau {tau}: {:e}", sampled_step_undershoot(tau));
        }
        let at = sampled_step_undershoot(49.575_535_071_522_86);
        assert!((at - 0.128_440_752_4).abs() < 1e-9, "the fixture's B {at:e}");
        let near_50 = (0..200).map(|k| sampled_step_undershoot(50.0 + f64::from(k) / 200.0)).fold(0.0_f64, f64::max);
        assert!(near_50 > 0.1380 && near_50 < 0.1382, "the largest B near 50 samples {near_50:e}");
        let near_1000 = (0..50).map(|k| sampled_step_undershoot(1000.0 + f64::from(k) / 50.0)).fold(0.0_f64, f64::max);
        assert!(near_1000 <= 0.1412, "the largest B near 1000 samples {near_1000:e}");
        let v = with_airway(voice());
        let (airway, tau) = (v.upstream_airway().unwrap(), v.upstream_airway_closures().unwrap().arrival_samples);
        assert!((tau - 49.575_535_071_522_86).abs() < 1e-9, "the fixture's arrival {tau:e}");
        let step = source_step_height(airway, tau);
        // The step height is a running sum of the taps at an index from the arrival on, and no later sum before the
        // first decrease is larger.
        let sums: Vec<f64> = airway.source_taps().iter().scan(0.0_f64, |s, g| { *s += g; Some(*s) }).collect();
        let from = tau.ceil() as usize;
        let peak = (from..sums.len()).find(|&k| k + 1 == sums.len() || sums[k + 1] < sums[k]).unwrap();
        println!("{{\"suite\":\"fs-couple\",\"case\":\"ua-flue-008\",\"tau\":{tau:e},\"b\":{at:e},\"b_max_near_50\":{near_50:e},\"b_max_near_1000\":{near_1000:e},\"step_height\":{step:e},\"step_sample\":{peak}}}");
        assert_eq!(step, sums[peak], "the first maximum from the arrival");
        assert!(step >= 1.0, "the fixture's chain narrows toward the windway: step {step:e}");
    }

    /// The two-way windway keeps the one-way bits wherever the windway passes forward flow (MB 18587, test a). No held
    /// step from rest on this realization keeps its open-circuit pressure positive (before the arrival the source FIR's
    /// band-limited pre-ringing goes negative, `ua_flue_002`), so the test is in two parts. A, the law: on 4000 positive
    /// open-circuit pressures from 1e-9 to 1e6 Pa, at several impedances and linear coefficients, the two-way closure
    /// returns the one-way closure's bits. B, the voice: on the fixture's held step, with no windway law and with the
    /// slit, a two-way and a one-way voice give the same output bits and frames on every sample before the first
    /// negative open-circuit pressure, and from that sample on the two-way voice has drawn air back (a negative flux).
    #[test]
    fn ua_flue_009_the_two_way_windway_keeps_the_forward_bits() {
        for (impedance, linear) in [(1.0e3, 0.0), (2.55e7, 0.0), (1.0e6, 6.0e3), (2.55e7, 6.0e3)] {
            for k in 0..1000 {
                let p = 1.0e-9 * 10f64.powf(f64::from(k) * 15.0 / 999.0);
                let port = UpstreamPort { open_circuit_pa: p, impedance };
                let one = foot_pressure(port, linear, 1.0e-5, 1.2).unwrap();
                let two = foot_pressure_two_way(port, linear, 1.0e-5, 1.2).unwrap();
                assert_eq!((one.0.to_bits(), one.1.to_bits()), (two.0.to_bits(), two.1.to_bits()), "P_oc {p:e}, Z {impedance:e}");
            }
        }
        for slit in [false, true] {
            let fresh = || if slit { with_airway(voice().with_viscous_windway(0.06).unwrap()) } else { with_airway(voice()) };
            let (mut one, mut two) = (fresh(), two_way(fresh()));
            assert!(two.upstream_airway_two_way() && !one.upstream_airway_two_way());
            let mut first_reverse = None;
            for n in 0..1024 {
                let (bits_one, frame_one) = held(&mut one);
                let (bits_two, frame_two) = held(&mut two);
                if first_reverse.is_none() && frame_two.port.open_circuit_pa < 0.0 {
                    first_reverse = Some(n);
                }
                match first_reverse {
                    None => assert_eq!((bits_one, frame_one), (bits_two, frame_two), "slit {slit}, sample {n}: the forward bits moved"),
                    Some(at) if at == n => assert!(frame_two.flow_m3_s < 0.0 && frame_one.flow_m3_s == 0.0,
                        "slit {slit}, sample {n}: the first negative open-circuit pressure drew no air back"),
                    Some(_) => {}
                }
            }
            let at = first_reverse.expect("a negative open-circuit pressure on the held step");
            println!("{{\"suite\":\"fs-couple\",\"case\":\"ua-flue-009\",\"slit\":{slit},\"same_bits_through\":{at},\"two_way\":{}}}",
                closures_json(two.upstream_airway_closures().expect("an airway")));
        }
    }

    /// The two-way closed form against a direct solve on its reverse samples (MB 18587, test b): the fixture's held step
    /// and the ramp of `ua_flue_001`, with no windway law and with the slit. On every sample whose open-circuit pressure
    /// is negative the foot pressure matches the bisection root of the signed law (the slit's own `evaluate` on a signed
    /// drop) to 1e-11 of that pressure, the airway's returned pressure matches the foot to 1e-9, the flux is negative,
    /// and the jet has none: zero stored speed and an infinite transit. Without reverse samples it fails.
    #[test]
    fn ua_flue_010_the_reverse_closed_form_is_the_direct_solve() {
        for slit in [false, true] {
            let mut results = Vec::new();
            for ramp in [false, true] {
                let mut v = two_way(if slit { with_airway(voice().with_viscous_windway(0.06).unwrap()) } else { with_airway(voice()) });
                let lungs: Vec<f64> = if ramp {
                    (0..480).map(|n| HOLD_PA * f64::from(n) / 480.0).chain(core::iter::repeat_n(HOLD_PA, 960)).collect()
                } else { vec![HOLD_PA; 1440] };
                let (mut worst, mut worst_valve, mut reverse) = (0.0_f64, 0.0_f64, 0usize);
                for lung in lungs {
                    let frame = step(&mut v, lung);
                    if frame.port.open_circuit_pa >= 0.0 { continue; }
                    reverse += 1;
                    let scale = -frame.port.open_circuit_pa;
                    worst = worst.max((frame.foot_pa - direct_reverse_foot(&v, frame.port)).abs() / scale);
                    worst_valve = worst_valve.max((frame.valve_pa - frame.foot_pa).abs() / scale);
                    assert!(frame.flow_m3_s < 0.0 && frame.foot_pa < 0.0 && frame.foot_pa >= frame.port.open_circuit_pa);
                    assert_eq!((v.jet_speed_m_s(), frame.jet_transit_s), (0.0, f64::INFINITY), "a reverse sample with a jet");
                }
                println!("{{\"suite\":\"fs-couple\",\"case\":\"ua-flue-010\",\"slit\":{slit},\"ramp\":{ramp},\"reverse_samples\":{reverse},\"worst_direct\":{worst:e},\"worst_valve\":{worst_valve:e},\"closures\":{}}}",
                    closures_json(v.upstream_airway_closures().expect("an airway")));
                results.push((ramp, reverse, worst, worst_valve));
            }
            for (ramp, reverse, worst, worst_valve) in results {
                assert!(reverse > 0, "slit {slit}, ramp {ramp}: no reverse sample to test");
                assert!(worst < 1.0e-11, "slit {slit}, ramp {ramp}: closed form against the direct solve {worst:e}");
                assert!(worst_valve < 1.0e-9, "slit {slit}, ramp {ramp}: the airway's pressure against the foot {worst_valve:e}");
            }
        }
    }

    /// The two-way law is odd and passes through zero (MB 18587, test c and its addition): on pressures from 1e-300 to
    /// 1e6 Pa at several impedances and linear coefficients, `u(-P) = -u(P)` and `P_f(-P) = -P_f(P)` bit for bit; at
    /// `P = 0` and `P = -0` both are exactly zero, finite, and the smallest pressures give finite speeds that shrink toward
    /// zero, so the closed form has no kink and no NaN at zero.
    #[test]
    fn ua_flue_011_the_two_way_law_is_odd_and_zero_at_zero() {
        for (impedance, linear) in [(1.0e3, 0.0), (2.55e7, 0.0), (1.0e6, 6.0e3), (2.55e7, 6.0e3)] {
            let law = |p: f64| foot_pressure_two_way(UpstreamPort { open_circuit_pa: p, impedance }, linear, 1.0e-5, 1.2).unwrap();
            for zero in [0.0, -0.0] {
                let (foot, speed) = law(zero);
                assert!(foot == 0.0 && speed == 0.0 && foot.is_finite() && speed.is_finite(), "Z {impedance:e}: at {zero}");
            }
            let mut previous = 0.0_f64;
            for k in 0..=1000 {
                let p = 1.0e-300 * 10f64.powf(f64::from(k) * 306.0 / 1000.0);
                let ((foot, speed), (foot_back, speed_back)) = (law(p), law(-p));
                assert_eq!((foot_back.to_bits(), speed_back.to_bits()), ((-foot).to_bits(), (-speed).to_bits()), "Z {impedance:e}: P {p:e}");
                assert!(speed.is_finite() && speed > 0.0 && speed >= previous, "Z {impedance:e}: P {p:e} speed {speed:e}");
                previous = speed;
            }
            assert!(law(1.0e-300).1 < 1.0e-290, "Z {impedance:e}: the speed does not shrink to zero");
        }
    }

    /// The two-way windway is opt-in on an attached airway before any sample, once; a voice without it runs one way.
    #[test]
    fn ua_flue_012_the_two_way_windway_is_one_cold_opt_in() {
        const WHAT: &str = "a two-way windway requires one cold call after the upstream airway";
        assert_eq!(refusal(voice().with_two_way_windway()), WHAT);
        assert_eq!(refusal(two_way(with_airway(voice())).with_two_way_windway()), WHAT);
        let mut v = with_airway(voice());
        step(&mut v, HOLD_PA);
        assert_eq!(refusal(v.with_two_way_windway()), WHAT);
        assert!(!voice().upstream_airway_two_way() && !with_airway(voice()).upstream_airway_two_way());
        assert!(two_way(with_airway(voice())).upstream_airway_two_way());
    }

    /// The developing slit's closure at zero impedance is the law's own held speed, bit for bit, and its foot the
    /// law's plenum pressure there. Against the chain's DC resistance, two ways it keeps the one-way bits on a positive
    /// open-circuit pressure and is odd bit for bit; one way a non-positive one passes no flow; a non-finite one refuses.
    #[test]
    fn ua_flue_013_the_developing_closure_at_zero_impedance_is_the_held_law_bit_for_bit() {
        let (g, gas) = (geometry(), gas());
        let law = DevelopingSlit::from_gas(&gas, 0.015).unwrap();
        let (h, area) = (g.flue_height_m, g.flue_width_m * g.flue_height_m);
        let r_dc = chain().dc_resistance(&gas);
        for p in [1.0e-3, 1.0, 11.0, HELD_PA, HOLD_PA, 2000.0] {
            let held = law.held(g.flue_width_m, h, p).unwrap();
            let (foot, speed) = developing_closure(&law, h, UpstreamPort { open_circuit_pa: p, impedance: 0.0 }, area, false).unwrap();
            assert_eq!(speed.to_bits(), held.bulk_speed_m_s.to_bits(), "{p} Pa: the held speed");
            assert_eq!(foot.to_bits(), law.against_series(h, p, 0.0).0.to_bits(), "{p} Pa: the plenum pressure");
            let port = UpstreamPort { open_circuit_pa: p, impedance: r_dc };
            let forward = developing_closure(&law, h, port, area, true).unwrap();
            assert_eq!(forward, developing_closure(&law, h, port, area, false).unwrap(), "{p} Pa: the one-way bits");
            assert!(forward.1 < speed && forward.0 < p, "{p} Pa: the airway lowers the foot");
            let back = UpstreamPort { open_circuit_pa: -p, impedance: r_dc };
            let reverse = developing_closure(&law, h, back, area, true).unwrap();
            assert_eq!(((-reverse.0).to_bits(), (-reverse.1).to_bits()), (forward.0.to_bits(), forward.1.to_bits()), "{p} Pa: odd");
            assert_eq!(developing_closure(&law, h, back, area, false).unwrap(), (0.0, 0.0), "{p} Pa: one way");
        }
        assert_eq!(developing_closure(&law, h, UpstreamPort { open_circuit_pa: 0.0, impedance: r_dc }, area, true).unwrap(), (0.0, 0.0));
        assert!(developing_closure(&law, h, UpstreamPort { open_circuit_pa: f64::NAN, impedance: r_dc }, area, true).is_err());
        assert!(developing_closure(&law, h, UpstreamPort { open_circuit_pa: 1.0, impedance: -1.0 }, area, true).is_err());
    }

    /// With an airway the developing windway is judged at the steady foot pressure, not at the lungs. A lung pressure
    /// above the lowest held pressure the slit admits by half the airway's steady drop there has its steady foot
    /// below that bound: admitted without the airway, refused with it by the bound's own name. Well inside the
    /// validity, the windway's held state becomes the steady foot's.
    #[test]
    fn ua_flue_014_the_developing_windway_is_judged_at_the_steady_foot() {
        let (g, gas) = (geometry(), gas());
        let law = DevelopingSlit::from_gas(&gas, 0.015).unwrap();
        let held = |p: f64| law.held(g.flue_width_m, g.flue_height_m, p).unwrap();
        // The lowest admitted held pressure on this slit (Chen's low-Reynolds bound), bisected between 1 and 50 Pa.
        let (mut lo, mut hi) = (1.0_f64, 50.0_f64);
        assert!(held(lo).validity().is_err() && held(hi).validity().is_ok());
        loop {
            let mid = 0.5 * (lo + hi);
            if mid <= lo || mid >= hi { break; }
            if held(mid).validity().is_ok() { hi = mid; } else { lo = mid; }
        }
        let series = chain().dc_resistance(&gas) * g.flue_width_m * g.flue_height_m;
        let lung = hi + 0.5 * series * held(hi).bulk_speed_m_s;
        let (foot, _) = law.against_series(g.flue_height_m, lung, series);
        assert!(foot < hi && held(lung).validity().is_ok(), "lung {lung}, foot {foot}, bound {hi}");
        let bound = held(foot).validity().unwrap_err();
        let ladder = LineLadder::default();
        let at_lung = voice().with_developing_windway(0.015, lung).unwrap();
        assert_eq!(refusal(at_lung.with_upstream_airway(&chain(), &ladder)), bound);
        let v = voice().with_developing_windway(0.015, HOLD_PA).unwrap().with_upstream_airway(&chain(), &ladder).unwrap();
        let (steady, _) = law.against_series(g.flue_height_m, HOLD_PA, series);
        assert!(steady < HOLD_PA);
        assert_eq!(v.developing_windway().unwrap().1, &held(steady));
    }
}
