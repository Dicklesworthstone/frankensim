//! The forced Bickley jet's causal response ([`fs_aeroac::forced_jet`]) as the island's labium
//! displacement (rulings 7376 and 7518; beads pexfp.135/.138/.321). At constant speed the
//! labium displacement is `eta_n = sum_k g_k v_(n-k)` [m], where `v` is the exit transverse velocity: the mouth
//! velocity `q / S_m` plus the onset seed. `g_k` [s] are the hat weights of one integrated step
//! response, and `g_0` multiplies the current mouth velocity, which is still unknown inside the
//! same-sample solve. It is opt-in by name and replaces the island's band chain and authored
//! displacement gain. The tanh cut, mouth loss and window are separate. It is not the default:
//! of the four gates of the switch it passes fg_004 but lands high on the three register gates
//! (fj_001, fg_003, gas_direction_amount; bead pexfp.135.1), and the default switches only once
//! all four pass unchanged.
//!
//! The exit profile is the Bickley jet with the exit's momentum flux and centre speed:
//! - A slit windway's Poiseuille exit gives `b = 0.4 h` and `U0 = 1.5 U`.
//! - No windway (a lip jet of uniform speed) gives `b = 0.75 h` and `U0 = U`.
//! - A developing slit ([`JetExit::developing`]) gives `U0` its inviscid core speed and `b` the
//!   momentum flux of the core and its two wall layers, between those two.
//!
//! This is Verge's matching as Segoufin states it (thesis, Univ. Paris 6, 2000, sec. 5.2
//! p. 69): the centre speed from Bernoulli and `b` from the momentum conserved between the
//! exit and the jet.
//!
//! The station is `cut-up / b`. The step response is integrated once, at the Reynolds number
//! `U0 b / nu` of the reference (held) pressure. At variable speed its clock is
//! `tau(t) = integral U0(t) / b dt`: each accepted interval retains its own advance.
//! A new speed never retroactively rescales old input intervals. This is the explicit
//! frozen-response extension `dz/dt = (U0(t)/b) A z + B v`, with the SAME fixed-Re
//! operator A and physical transverse velocity v, not a new variable-Re flow solve.
//! Velocity is linear within each mechanical interval and speed is held on that interval.
//! At zero speed the jet output/forcing is disabled and this response's clock freezes;
//! molecular decay of a stopped jet is outside this fixed-Re image. A trajectory with
//! constant positive speed from rest retains the original cached hat arithmetic exactly.
//! Any positive speed is admitted (ruling 8132, bead pexfp.135.2): input older than the first
//! accepted interval is the rest state, and a jet whose retained intervals cannot reach the
//! window freezes this clock as at rest while its voice takes the declared state.
//!
//! Not claimed:
//! - The base jet's own spreading over the cut-up.
//! - The labium's hydrodynamic back-reaction.
//! - The Reynolds number of attack and release, and a developing exit's profile there: it is
//!   held at the reference pressure.
//! - The exit's shape beyond its centre speed and momentum flux. A developing exit is flatter
//!   than a Bickley jet: de la Cuadra (thesis, Stanford, 2005) fits Nolle profiles with
//!   n = 2-4 1 mm past a 1 mm flue.
//! - Memory beyond the window, where the response is held.

use std::sync::{Mutex, PoisonError};

use fs_aeroac::AeroacError;
use fs_aeroac::forced_jet::{
    FRONT_MARGIN, ForcedJetResolution, JetStepResponse, bickley_step_response,
};

use crate::bernoulli_aperture::developing::DevelopingSlitState;

/// Bickley half-width over flue height of a Poiseuille slit exit (its momentum flux at 1.5 x
/// the bulk speed).
pub const POISEUILLE_HALF_WIDTH_OVER_HEIGHT: f64 = 0.4;
/// Centre over bulk speed of a Poiseuille slit exit.
pub const POISEUILLE_CENTRE_OVER_BULK: f64 = 1.5;
/// Bickley half-width over flue height of a uniform exit (its momentum flux at the bulk speed).
pub const UNIFORM_HALF_WIDTH_OVER_HEIGHT: f64 = 0.75;
/// Shape factor `delta*/theta` of a developing slit's wall layers: the parabolic layer, which is
/// the Poiseuille profile once the two layers meet. Van Zon's linear layer (3, Segoufin 2000
/// app. A) would move the tin whistle's `b` by at most 3.4 %.
pub const DEVELOPING_LAYER_SHAPE_FACTOR: f64 = 2.5;

/// The windway exit the Bickley jet is matched to: its centre speed over the bulk speed, and the
/// half-width over flue height that carries the exit's momentum flux at that centre speed.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct JetExit {
    /// Centre over bulk speed `U0/U`.
    pub centre_over_bulk: f64,
    /// Bickley half-width over flue height `b/h`.
    pub half_width_over_height: f64,
}

impl JetExit {
    /// A lip jet of uniform speed (no windway).
    pub const UNIFORM: Self = Self {
        centre_over_bulk: 1.0,
        half_width_over_height: UNIFORM_HALF_WIDTH_OVER_HEIGHT,
    };
    /// A fully developed slit's Poiseuille exit.
    pub const POISEUILLE: Self = Self {
        centre_over_bulk: POISEUILLE_CENTRE_OVER_BULK,
        half_width_over_height: POISEUILLE_HALF_WIDTH_OVER_HEIGHT,
    };

    /// A developing slit's exit at its held state. The centre leaves at the inviscid core speed
    /// `u_c`, and each wall layer displaces `delta* = (h/2)(1 - U/u_c)` of Shah's bulk speed `U`.
    /// The momentum flux `rho u_c^2 (h - 2 delta* - 2 theta)`, with `theta = delta*/H`, gives
    /// `b = 0.75 (h - 2 delta* (1 + 1/H))`. An exit whose layers have met is the Poiseuille exit.
    #[must_use]
    pub fn developing(state: &DevelopingSlitState, flue_height_m: f64) -> Self {
        if state.developed {
            return Self::POISEUILLE;
        }
        let layers = 2.0 * state.displacement_thickness_m
            * (1.0 + 1.0 / DEVELOPING_LAYER_SHAPE_FACTOR) / flue_height_m;
        Self {
            centre_over_bulk: state.core_speed_m_s / state.bulk_speed_m_s,
            half_width_over_height: UNIFORM_HALF_WIDTH_OVER_HEIGHT * (1.0 - layers),
        }
    }
}
/// Half-height of the resolved vortical band [b]; the exact exterior closes it.
pub const BAND_HALF_HEIGHT: f64 = 8.0;
/// Centreline Courant number of the integration.
pub const COURANT: f64 = 0.8;
/// Most worker threads for the integration (the result does not depend on the count).
pub const MAX_WORKERS: usize = 8;
/// The production row, `jet_response forced-bickley 0.1 90`, on which register physics is graded.
pub const PRODUCTION_SPEC: ForcedJetSpec = ForcedJetSpec { spacing: 0.1, window: 90.0 };
/// The coarse row numerical fixtures declare: the same physics on a coarse grid and a shorter
/// window. On the whistle D5 card it renders within 0.3 cents and 0.3 dB of [`PRODUCTION_SPEC`].
pub const COARSE_SPEC: ForcedJetSpec = ForcedJetSpec { spacing: 0.25, window: 40.0 };
/// Most integrated responses one process retains for reuse.
pub const RETAINED_RESPONSES: usize = 32;

// Integrated responses by their exact inputs. The integration is deterministic and its bits do
// not depend on the worker count, so a voice built again from the same flue, gas and row reuses
// the identical response instead of integrating it again.
static RESPONSES: Mutex<Vec<([u64; 5], JetStepResponse)>> = Mutex::new(Vec::new());

fn integrated(station: f64, resolution: &ForcedJetResolution)
    -> Result<JetStepResponse, AeroacError>
{
    let key = [station.to_bits(), resolution.spacing.to_bits(), resolution.window.to_bits(),
        resolution.reynolds.to_bits(), resolution.streamwise_cells as u64];
    let retained = || RESPONSES.lock().unwrap_or_else(PoisonError::into_inner);
    if let Some((_, step)) = retained().iter().find(|(k, _)| *k == key) {
        return Ok(step.clone());
    }
    let step = bickley_step_response(station, resolution)?;
    let mut cache = retained();
    if !cache.iter().any(|(k, _)| *k == key) {
        if cache.len() == RETAINED_RESPONSES {
            cache.remove(0);
        }
        cache.push((key, step.clone()));
    }
    Ok(step)
}

/// Card row `jet_response forced-bickley SPACING WINDOW`: the grid spacing [b] and the recorded
/// window [b/U0] of the step response.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ForcedJetSpec {
    /// Grid spacing [b].
    pub spacing: f64,
    /// Recorded window [b/U0].
    pub window: f64,
}

/// One accepted physical-time velocity interval and its convective-clock advance.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
struct JetInterval {
    start_velocity_m_s: f64,
    end_velocity_m_s: f64,
    advance: f64,
}

/// The integrated response, its jet scales and the retained exit-velocity intervals.
#[derive(Debug, Clone)]
pub struct ForcedJetKernel {
    spec: ForcedJetSpec,
    step: JetStepResponse,
    half_width_m: f64,
    centre_over_bulk: f64,
    reynolds: f64,
    dt: f64,
    ring: Vec<JetInterval>,
    head: usize,
    used: usize,
    last_velocity_m_s: f64,
    uniform_bulk_bits: Option<u64>,
    uniform_intervals: usize,
    nonuniform_history: bool,
    prepared_bulk_bits: u64,
    prepared_advance: f64,
    // The response's weight count at the cached speed; `weights` holds at most the ring's
    // length of its leading weights.
    weight_count: usize,
    weights: Vec<f64>,
    weights_bulk_bits: u64,
    current: (f64, f64),
}

fn refusal(error: AeroacError) -> &'static str {
    match error {
        AeroacError::NonFinite { what } | AeroacError::InvalidParameter { what } => what,
        AeroacError::NotConverged { .. } => "forced jet integration did not converge",
    }
}

impl ForcedJetKernel {
    /// Integrate the step response for this flue. `exit` is the windway exit the jet is matched
    /// to; `reference_bulk_m_s` is the bulk speed at the held pressure, which sets the Reynolds
    /// number; `minimum_bulk_m_s` sizes the retained input history, which may not exceed
    /// `max_taps` samples.
    ///
    /// # Errors
    /// Non-positive or non-finite inputs, a resolution the integrator refuses, or more retained
    /// samples at the slowest jet than `max_taps`.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        spec: ForcedJetSpec,
        exit: JetExit,
        flue_height_m: f64,
        cut_up_m: f64,
        kinematic_viscosity_m2_s: f64,
        reference_bulk_m_s: f64,
        minimum_bulk_m_s: f64,
        max_taps: usize,
        sample_period_s: f64,
    ) -> Result<Self, &'static str> {
        let (width_ratio, centre_over_bulk) = (exit.half_width_over_height, exit.centre_over_bulk);
        let half_width_m = width_ratio * flue_height_m;
        let inputs = [
            spec.spacing,
            spec.window,
            centre_over_bulk,
            half_width_m,
            cut_up_m,
            kinematic_viscosity_m2_s,
            reference_bulk_m_s,
            minimum_bulk_m_s,
            sample_period_s,
        ];
        if !inputs.iter().all(|v| v.is_finite() && *v > 0.0) {
            return Err("forced jet needs positive finite geometry, gas, speeds and window");
        }
        let station = cut_up_m / half_width_m;
        let reynolds =
            centre_over_bulk * reference_bulk_m_s * half_width_m / kinematic_viscosity_m2_s;
        let reach = (spec.window + FRONT_MARGIN).max(station + FRONT_MARGIN);
        let cells = ((reach / spec.spacing).ceil() as usize)
            .max(32)
            .next_power_of_two();
        let workers = std::thread::available_parallelism()
            .map_or(1, std::num::NonZeroUsize::get)
            .min(MAX_WORKERS);
        let resolution = ForcedJetResolution {
            spacing: spec.spacing,
            streamwise_cells: cells,
            half_height: BAND_HALF_HEIGHT,
            window: spec.window,
            courant: COURANT,
            reynolds,
            workers,
        };
        let step = integrated(station, &resolution).map_err(refusal)?;
        let slowest = step
            .hat_weights(
                centre_over_bulk * minimum_bulk_m_s,
                half_width_m,
                sample_period_s,
            )
            .map_err(refusal)?
            .len();
        if slowest > max_taps {
            return Err("forced jet history at the slowest admitted jet exceeds the declared jet history budget");
        }
        Ok(Self {
            spec,
            step,
            half_width_m,
            centre_over_bulk,
            reynolds,
            dt: sample_period_s,
            ring: vec![JetInterval::default(); slowest],
            head: 0,
            used: 0,
            last_velocity_m_s: 0.0,
            uniform_bulk_bits: None,
            uniform_intervals: 0,
            nonuniform_history: false,
            prepared_bulk_bits: 0,
            prepared_advance: 0.0,
            weight_count: 0,
            weights: Vec::new(),
            weights_bulk_bits: 0,
            current: (0.0, 0.0),
        })
    }

    /// The card row this response was integrated from.
    #[must_use]
    pub const fn spec(&self) -> ForcedJetSpec {
        self.spec
    }

    /// The integrated step response.
    #[must_use]
    pub const fn step_response(&self) -> &JetStepResponse {
        &self.step
    }

    /// Bickley half-width of the exit [m].
    #[must_use]
    pub const fn half_width_m(&self) -> f64 {
        self.half_width_m
    }

    /// Centre over bulk speed of the exit.
    #[must_use]
    pub const fn centre_over_bulk(&self) -> f64 {
        self.centre_over_bulk
    }

    /// `U0 b / nu` the response was integrated at.
    #[must_use]
    pub const fn reynolds(&self) -> f64 {
        self.reynolds
    }

    /// Held displacement per exit velocity [s] at bulk speed `bulk_m_s`: `b s(window) / U0`.
    #[must_use]
    pub fn held_gain_s(&self, bulk_m_s: f64) -> f64 {
        self.half_width_m / (self.centre_over_bulk * bulk_m_s) * self.step.held_gain()
    }

    /// Prepare the known displacement [m] and current-velocity coefficient [s].
    /// Accepted intervals keep their individual convective ages when speed changes.
    /// Repeated trials do not advance the clock or publish any velocity.
    ///
    /// A still jet emits no displacement and does not advance this fixed-Re clock.
    ///
    /// Any positive speed is admitted (ruling 8132). Intervals older than every accepted one
    /// precede the jet's first motion and read the rest state, exactly zero, so the response is
    /// exact whenever the retained intervals reach either its window or past every accepted
    /// interval: unbounded history at zero memory. Otherwise (a jet slower than the declared
    /// minimum for longer than the ring retains, or one too slow to advance this clock) this
    /// returns `None`, prepares no displacement and freezes the clock as at rest; the voice then
    /// takes its declared state below its minimum speed, or refuses.
    ///
    /// # Errors
    /// Nonfinite/negative speed or nonfinite response. A refusal leaves the accepted history
    /// unchanged.
    pub fn prepare(&mut self, bulk_m_s: f64) -> Result<Option<(f64, f64)>, &'static str> {
        if !bulk_m_s.is_finite() || bulk_m_s < 0.0 {
            return Err("forced jet speed must be finite and nonnegative");
        }
        self.prepared_bulk_bits = 0;
        self.prepared_advance = 0.0;
        self.current = (0.0, 0.0);
        if bulk_m_s == 0.0 {
            return Ok(Some(self.current));
        }
        let bits = bulk_m_s.to_bits();
        let advance = self.centre_over_bulk * bulk_m_s * self.dt / self.half_width_m;
        if !advance.is_finite() {
            return Err("forced jet convective advance must be finite and positive");
        }
        if advance <= 0.0 {
            return Ok(None);
        }
        let n = self.ring.len();
        let current = if !self.nonuniform_history
            && self.uniform_bulk_bits.is_none_or(|old| old == bits)
        {
            // Keep the original summation and interpolation for a uniformly moving jet.
            if bits != self.weights_bulk_bits || self.weights.is_empty() {
                // No read reaches past the ring, so at most its length is weighted.
                (self.weight_count, self.weights) = self.step
                    .leading_hat_weights(self.centre_over_bulk * bulk_m_s, self.half_width_m,
                        self.dt, n)
                    .map_err(refusal)?;
                self.weights_bulk_bits = bits;
            }
            let taps = if self.weight_count <= n {
                self.weight_count
            } else if self.used < n {
                self.weight_count.min(self.used + 1)
            } else {
                return Ok(None);
            };
            let mut known = 0.0;
            for (k, g) in self.weights.iter().enumerate().take(taps).skip(1) {
                known += g * self.ring[(self.head + n - k) % n].end_velocity_m_s;
            }
            (known, self.weights[0])
        } else {
            let (current_weight, previous_weight) = self.step
                .interval_hat_weights(0.0, advance, self.dt).map_err(refusal)?;
            let mut known = previous_weight * self.last_velocity_m_s;
            let mut age = advance;
            for k in 0..self.used {
                if age >= self.step.window() { break; }
                let interval = self.ring[(self.head + n - 1 - k) % n];
                let (newer, older) = self.step
                    .interval_hat_weights(age, interval.advance, self.dt).map_err(refusal)?;
                known += newer * interval.end_velocity_m_s + older * interval.start_velocity_m_s;
                age += interval.advance;
            }
            if age < self.step.window() && self.used == n {
                return Ok(None);
            }
            (known, current_weight)
        };
        if !current.0.is_finite() || !current.1.is_finite() {
            return Err("forced jet response left the finite displacement domain");
        }
        self.prepared_bulk_bits = bits;
        self.prepared_advance = advance;
        self.current = current;
        Ok(Some(current))
    }

    /// The prepared sample's known displacement [m] and current-velocity weight [s]:
    /// `eta = known + g_0 v` for this sample's exit velocity `v`.
    #[must_use]
    pub const fn current(&self) -> (f64, f64) {
        self.current
    }

    /// Commit the prepared interval and its accepted exit velocity [m/s].
    /// The coupled voice validates velocity before this publication. At rest, or after a
    /// [`Self::prepare`] that returned `None`, no velocity is injected; old intervals retain
    /// their convective ages.
    pub fn push(&mut self, velocity_m_s: f64) {
        if self.prepared_advance == 0.0 {
            self.nonuniform_history |= self.used != 0;
            self.uniform_intervals = 0;
            self.last_velocity_m_s = 0.0;
            return;
        }
        let same_speed = self.uniform_bulk_bits == Some(self.prepared_bulk_bits);
        self.nonuniform_history |= self.uniform_bulk_bits.is_some() && !same_speed;
        self.uniform_intervals = if same_speed {
            (self.uniform_intervals + 1).min(self.ring.len())
        } else { 1 };
        self.uniform_bulk_bits = Some(self.prepared_bulk_bits);
        // Once the complete bounded history is uniform again, the original
        // cached convolution is valid without reinterpreting any old interval.
        if self.uniform_intervals == self.ring.len() {
            self.nonuniform_history = false;
        }
        self.ring[self.head] = JetInterval {
            start_velocity_m_s: self.last_velocity_m_s,
            end_velocity_m_s: velocity_m_s,
            advance: self.prepared_advance,
        };
        self.head = (self.head + 1) % self.ring.len();
        self.used = (self.used + 1).min(self.ring.len());
        self.last_velocity_m_s = velocity_m_s;
    }
}

#[cfg(test)]
mod tests;
