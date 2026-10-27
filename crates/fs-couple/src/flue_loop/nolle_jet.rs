//! Opt-in hosted jet transfer of Nolle (1998) as published by Price, Johnston and McKinnon,
//! J. Acoust. Soc. Am. 138, 3282 (2015), arXiv:1502.02170, Eqs. 4-10 and Fig. 5.
//!
//! The island's default receptivity chain (two low-pass poles, a high-pass, a leaky integral and an
//! authored gain) is NOT touched; a description selects this transfer by name (`jet_transfer
//! nolle-price-2015 ...`). The flue jet then answers the mouth flow `Q` (into the pipe) with the
//! published linear law, in Price's variables
//!
//! ```text
//! q_j  = J(w) (alpha_p q + alpha_j q_j)            (Eq. 4;  q = -Q, the flow into the head)
//! J(w) = -g (U0 / W~) (1 / (i w)) (1 - exp((mu - i w / c_p) W~))   (Eq. 10)
//! Z_h  = r + i w M + i w M_d alpha_p J / (1 - alpha_j J)            (Eq. 5, M_d = 0.88 rho / H)
//! ```
//!
//! so that with `Xi(t) = integral of (alpha_p Q - alpha_j dQ_j) / (H W~) dt` the fluctuating jet flow
//! is `dQ_j = g U0 H (Xi - E * Xi)` with `E(w) = exp(mu(w) W~) exp(-i w W~ / c_p)`. The island turns
//! `dQ_j` into the labium displacement `eta = dQ_j / (u H)` (`u` = bulk speed, `U0 = 1.5 u` for the
//! Poiseuille windway), cuts the jet with its own tanh law and drives the mouth with the source
//! `M_d dQ_j/dt`.
//!
//! PROVENANCE (no number here is fitted to audio):
//! - `mu b (Str_b)` and the theory `c_p / U0 (Str_b)` curves are Price's Fig. 5 (Bickley jet,
//!   Mattingly and Criminale; Nolle's tables), READ BY EYE from the published figure. The theory curve
//!   is kept only for the closed-form gate check: with it Price's measured gain bands are missed by
//!   +30.6 / +7.5 / +15.4 / +28.2 % (his own model edges 0.0705 / 0.2871 / 0.4964 / 0.7690 against the
//!   measured 0.054 / 0.267 / 0.43 / 0.60), so the host uses a MEASURED phase velocity instead.
//! - `c_p / U_j` = 0.3 for a long formation channel and 0.5 for a short one is de la Cuadra's
//!   measurement (Stanford thesis 2005, Sec. 2.2.6, Fig. 2.14a, CO2 jets at Re about 200, constant over
//!   Str_h 0.05-1.1). It was adopted AFTER the theory table failed Price's band, so Price's band is
//!   consistency, not an independent test. The class is chosen per jet geometry, a priori, from the
//!   thesis' own classification (lips: "a very short channel", Sec. 3.1).
//! - `g = 0.145` is Price's fitted jet gain, measured with the pipe flow BLOCKED; his derivation assumes
//!   1 for an unblocked jet and he could not measure it. The playing-condition gain is UNKNOWN.
//! - `alpha_p = 2/pi`, `alpha_j = 0.38`, `M_d = 0.88 rho / H` are Verge's two-dimensional potential-flow
//!   constants as used by Price; `b = 2h/5` is Price's momentum-conserving Bickley width.
//!
//! NUMERICS: `E` is realised as one delay-line tap at `W~/c_p` plus the tap bump of the quadrature
//! of `exp(mu W~) - 1` (inverse transform of the tabulated growth). The bump is symmetric about the
//! carrier and has mass before the jet leaves the flue; that acausal part (about 1e-4 of the energy) is
//! lumped at the first causal tap and the tail beyond the retained span at the last tap, which keeps the
//! low-frequency gain. The realisation error against the closed form is measured at construction and
//! reported; construction refuses above 15 %. The leaky integral is advanced exactly for held
//! input, with leak `0.01 U0 / W~` (below Price's first band edge, St 0.054). A read-only
//! affine next-step response supports simultaneous jet/mouth coupling without advancing history.
//! Gauntlet tier: G0/G1/G2 (held-input integration, closed-form edges,
//! time-domain lock-in against the closed form). NO CLAIM: linear small-signal law only (the labium cut
//! supplies the saturation), Bickley-profile growth without spreading or Reynolds dependence.

use core::f64::consts::{PI, TAU};

/// Price Fig. 5: Strouhal number `Str_b = w b / U0` of the `mu b` samples (read by eye).
pub const MU_B_STROUHAL: [f64; 16] = [
    0.0, 0.05, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.33,
];
/// Price Fig. 5: dimensionless spatial growth rate `mu b` of the sinuous mode (read by eye).
pub const MU_B: [f64; 16] = [
    0.04, 0.14, 0.19, 0.265, 0.275, 0.275, 0.25, 0.225, 0.19, 0.165, 0.135, 0.105, 0.08, 0.055,
    0.03, 0.0,
];
/// Price Fig. 5: Strouhal number of the theory phase-velocity samples (read by eye).
pub const CP_STROUHAL: [f64; 12] = [0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.33];
/// Price Fig. 5: theory phase velocity over the centre speed (read by eye).
pub const CP_OVER_U0: [f64; 12] = [
    0.295, 0.345, 0.37, 0.39, 0.42, 0.45, 0.475, 0.50, 0.55, 0.60, 0.64, 0.665,
];
/// Price Eq. 4: pipe-flow coupling of the disturbance (Verge, 2/pi).
pub const ALPHA_P: f64 = 2.0 / PI;
/// Price Eq. 4: direct jet-source coupling of the disturbance (Verge).
pub const ALPHA_J: f64 = 0.38;
/// Price Sec. III: drive inertance `M_d = 0.88 rho / H` (Verge's 2-D potential flow).
pub const DRIVE_INERTANCE_COEFF: f64 = 0.88;
/// Price Sec. V: jet gain fitted to the blocked-pipe reflectometer data.
pub const PRICE_JET_GAIN: f64 = 0.145;
/// de la Cuadra 2005 Fig. 2.14a: `c_p / U_j` of a long (10 mm) formation channel.
pub const LONG_CHANNEL_CP_OVER_U0: f64 = 0.3;
/// de la Cuadra 2005 Fig. 2.14a: `c_p / U_j` of a short (1 mm) formation channel.
pub const SHORT_CHANNEL_CP_OVER_U0: f64 = 0.5;
/// Price Sec. III: Bickley width `b = 2h/5` from momentum conservation.
pub const BICKLEY_WIDTH_OVER_HEIGHT: f64 = 0.4;
/// Poiseuille centre speed over bulk speed (`Q_j = 2/3 h U0 H`, Price Sec. III).
pub const POISEUILLE_CENTRE_OVER_BULK: f64 = 1.5;

const LEAK_REL: f64 = 0.01;
const KERNEL_STEP_HAT: f64 = 0.02;
const KERNEL_TABLE_SPAN_HAT: f64 = 12.0;
const KERNEL_LATE_SPAN_HAT: f64 = 8.0;
const KERNEL_QUADRATURE_POINTS: usize = 4000;
const TAP_SPACING_HAT: f64 = 0.05;
const MAX_REALISATION_ERROR: f64 = 0.15;
const MAX_RING_SAMPLES: usize = 1 << 23;

fn interpolate(x: f64, xs: &[f64], ys: &[f64]) -> f64 {
    if x <= xs[0] {
        return ys[0];
    }
    for i in 1..xs.len() {
        if x <= xs[i] {
            return ys[i - 1] + (ys[i] - ys[i - 1]) * (x - xs[i - 1]) / (xs[i] - xs[i - 1]);
        }
    }
    ys[ys.len() - 1]
}

/// Price Fig. 5 growth rate `mu b` at `Str_b = w b / U0`.
#[must_use]
pub fn mu_b(strouhal_b: f64) -> f64 {
    interpolate(strouhal_b.abs(), &MU_B_STROUHAL, &MU_B)
}

/// Price Fig. 5 theory phase velocity `c_p / U0` at `Str_b`.
#[must_use]
pub fn cp_over_u0_theory(strouhal_b: f64) -> f64 {
    interpolate(strouhal_b.abs(), &CP_STROUHAL, &CP_OVER_U0)
}

/// Phase-velocity law of the closed-form jet transfer.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum PhaseVelocity {
    /// Price's Fig. 5 theory curve (fails the measured gain band; kept for the gate check).
    TheoryTable,
    /// A constant `c_p / U0` (the measured channel-class value).
    Constant(f64),
}

impl PhaseVelocity {
    fn over_u0(self, strouhal_b: f64) -> f64 {
        match self {
            PhaseVelocity::TheoryTable => cp_over_u0_theory(strouhal_b),
            PhaseVelocity::Constant(c) => c,
        }
    }
}

/// Price Eq. 10 as `(re, im)`, with `omega_hat = w W~ / U0` and `w_over_b = W~ / b`
/// (so `U0 / W~ = 1`): `J = g/omega_hat * (-G sin(phi) + i (1 - G cos(phi)))`,
/// `G = exp(mu W~)`, `phi = omega_hat / (c_p / U0)`.
#[must_use]
pub fn jet_transfer_closed_form(
    omega_hat: f64,
    w_over_b: f64,
    phase: PhaseVelocity,
    gain: f64,
) -> (f64, f64) {
    let strouhal_b = omega_hat / w_over_b;
    let growth = (mu_b(strouhal_b) * w_over_b).exp();
    let angle = omega_hat / phase.over_u0(strouhal_b);
    (
        -gain * growth * angle.sin() / omega_hat,
        gain * (1.0 - growth * angle.cos()) / omega_hat,
    )
}

/// Strouhal numbers `St = f W~ / U0` in `[st_lo, st_hi]` where `Im J` changes sign (the edges of the
/// amplifier's gain bands: `|S_h|^2 > 1` where `Im J > 0`, Price Sec. III).
#[must_use]
pub fn gain_band_edges(
    w_over_b: f64,
    phase: PhaseVelocity,
    st_lo: f64,
    st_hi: f64,
    steps: usize,
) -> Vec<f64> {
    let mut edges = Vec::new();
    let mut previous: Option<bool> = None;
    for i in 0..=steps {
        let st = st_lo + (st_hi - st_lo) * (i as f64) / (steps as f64);
        let (_, im) = jet_transfer_closed_form(TAU * st, w_over_b, phase, 1.0);
        let gain = im > 0.0;
        if let Some(before) = previous {
            if before != gain {
                edges.push(st);
            }
        }
        previous = Some(gain);
    }
    edges
}

/// de la Cuadra's channel classes for the measured phase velocity.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ChannelClass {
    /// A long formation channel (windway): `c_p = 0.3 U_j`.
    Long,
    /// A short channel such as the player's lips: `c_p = 0.5 U_j`.
    Short,
}

impl ChannelClass {
    /// The measured `c_p / U_j` of the class.
    #[must_use]
    pub fn phase_velocity_over_u0(self) -> f64 {
        match self {
            ChannelClass::Long => LONG_CHANNEL_CP_OVER_U0,
            ChannelClass::Short => SHORT_CHANNEL_CP_OVER_U0,
        }
    }
    /// The description-file token of the class.
    #[must_use]
    pub fn token(self) -> &'static str {
        match self {
            ChannelClass::Long => "long-channel",
            ChannelClass::Short => "short-channel",
        }
    }
}

/// The published constants of the hosted transfer (fixed; a description chooses only the class and
/// the chamfer depth, never a gain).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NolleSpec {
    /// Channel class of the jet.
    pub channel: ChannelClass,
    /// Measured `c_p / U0` of the class.
    pub phase_velocity_over_u0: f64,
    /// Jet gain `g` (Price: 0.145, blocked pipe).
    pub jet_gain: f64,
    /// Pipe-flow coupling `alpha_p`.
    pub alpha_p: f64,
    /// Jet-source coupling `alpha_j`.
    pub alpha_j: f64,
    /// Drive inertance coefficient (`M_d = coeff * rho / H`).
    pub drive_inertance_coeff: f64,
    /// Chamfer depth added to the cut-up for the travel distance `W~ = W + d` [m].
    pub chamfer_depth_m: f64,
    /// Integrator leak over `U0 / W~` (numerical regularisation).
    pub leak_rel: f64,
}

impl NolleSpec {
    /// Price's constants with the measured phase velocity of the channel class.
    ///
    /// # Errors
    /// A negative or non-finite chamfer depth.
    pub fn price_2015(
        channel: ChannelClass,
        chamfer_depth_m: f64,
    ) -> Result<NolleSpec, &'static str> {
        if !(chamfer_depth_m.is_finite() && chamfer_depth_m >= 0.0) {
            return Err("chamfer depth must be finite and non-negative");
        }
        Ok(NolleSpec {
            channel,
            phase_velocity_over_u0: channel.phase_velocity_over_u0(),
            jet_gain: PRICE_JET_GAIN,
            alpha_p: ALPHA_P,
            alpha_j: ALPHA_J,
            drive_inertance_coeff: DRIVE_INERTANCE_COEFF,
            chamfer_depth_m,
            leak_rel: LEAK_REL,
        })
    }
}

/// Excess growth kernel `g_s(t') = (1/pi) int_0^Omega (exp(mu W~) - 1) cos(w t') dw` on the grid
/// `t' = i KERNEL_STEP_HAT` (units `W~ / U0`), `Omega = 1.33 W~ / b` (`mu b` vanishes beyond).
fn excess_kernel(w_over_b: f64) -> Vec<f64> {
    let omega_max = MU_B_STROUHAL[MU_B_STROUHAL.len() - 1] * w_over_b;
    let dw = omega_max / KERNEL_QUADRATURE_POINTS as f64;
    let mut omegas = Vec::with_capacity(KERNEL_QUADRATURE_POINTS);
    let mut excess = Vec::with_capacity(KERNEL_QUADRATURE_POINTS);
    for k in 0..KERNEL_QUADRATURE_POINTS {
        let w = (k as f64 + 0.5) * dw;
        omegas.push(w);
        excess.push((mu_b(w / w_over_b) * w_over_b).exp() - 1.0);
    }
    let count = (KERNEL_TABLE_SPAN_HAT / KERNEL_STEP_HAT).floor() as usize + 1;
    let mut kernel = Vec::with_capacity(count);
    for i in 0..count {
        let t = i as f64 * KERNEL_STEP_HAT;
        let mut sum = 0.0;
        for (w, e) in omegas.iter().zip(excess.iter()) {
            sum += e * (w * t).cos();
        }
        kernel.push(sum * dw / PI);
    }
    kernel
}

/// The hosted jet state of one voice.
#[derive(Clone, Debug)]
pub struct NolleJet {
    spec: NolleSpec,
    centre_over_bulk: f64,
    w_tilde_m: f64,
    w_over_b: f64,
    tau_hat: f64,
    kernel: Vec<f64>,
    early_index: usize,
    late_index: usize,
    early_lump: f64,
    late_lump: f64,
    realisation_error: f64,
    dt: f64,
    xi: f64,
    ring: Vec<f64>,
    head: usize,
}

/// Read-only next-step condensation of the retained jet history. The two flows
/// are held over the mechanical interval and are unknowns of the host's solve.
/// In particular, the jet's direct hydrodynamic feedback is not yesterday's flow.
#[derive(Clone, Copy, Debug)]
pub(crate) struct NolleStep {
    xi_history_m: f64,
    xi_mouth_s_m2: f64,
    xi_jet_s_m2: f64,
    eta_history_m: f64,
    pub(crate) eta_mouth_s_m2: f64,
    pub(crate) eta_jet_s_m2: f64,
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct NolleCandidate {
    xi_m: f64,
    pub(crate) eta_m: f64,
}

impl NolleStep {
    pub(crate) fn candidate(
        self,
        mouth_flow_m3_s: f64,
        jet_flow_m3_s: f64,
    ) -> Result<NolleCandidate, &'static str> {
        let xi_m = self.xi_history_m
            + self.xi_mouth_s_m2 * mouth_flow_m3_s
            + self.xi_jet_s_m2 * jet_flow_m3_s;
        let eta_m = self.eta_history_m
            + self.eta_mouth_s_m2 * mouth_flow_m3_s
            + self.eta_jet_s_m2 * jet_flow_m3_s;
        if ![mouth_flow_m3_s, jet_flow_m3_s, xi_m, eta_m]
            .iter()
            .all(|v| v.is_finite())
        {
            return Err("nolle jet trial left the finite domain");
        }
        Ok(NolleCandidate { xi_m, eta_m })
    }
}

impl NolleJet {
    /// Build the transfer for one flue geometry.
    ///
    /// `centre_over_bulk` is `U0 / u` (1.5 for a Poiseuille windway, 1 for the ideal Bernoulli jet);
    /// `min_centre_speed_m_s` sizes the history (slower jets are refused by `step`).
    ///
    /// # Errors
    /// Non-finite or non-positive geometry, a history over 2^23 samples, or a tap realisation that
    /// misses the closed form by more than 15 % over St 0.03-0.7.
    pub fn new(
        spec: NolleSpec,
        centre_over_bulk: f64,
        flue_height_m: f64,
        cut_up_m: f64,
        dt: f64,
        min_centre_speed_m_s: f64,
    ) -> Result<NolleJet, &'static str> {
        let positive = |x: f64| x.is_finite() && x > 0.0;
        if !(positive(centre_over_bulk)
            && positive(flue_height_m)
            && positive(cut_up_m)
            && positive(dt)
            && positive(min_centre_speed_m_s)
            && positive(spec.phase_velocity_over_u0)
            && positive(spec.jet_gain)
            && positive(spec.leak_rel)
            && spec.alpha_p.is_finite()
            && spec.alpha_j.is_finite()
            && positive(spec.drive_inertance_coeff)
            && spec.chamfer_depth_m.is_finite()
            && spec.chamfer_depth_m >= 0.0)
        {
            return Err("nolle jet needs finite positive geometry, speeds and constants");
        }
        let w_tilde_m = cut_up_m + spec.chamfer_depth_m;
        let w_over_b = w_tilde_m / (BICKLEY_WIDTH_OVER_HEIGHT * flue_height_m);
        let tau_hat = 1.0 / spec.phase_velocity_over_u0;
        let kernel = excess_kernel(w_over_b);
        let last = kernel.len() - 1;
        let early_index = ((tau_hat / KERNEL_STEP_HAT).floor() as usize).min(last);
        let late_index = ((KERNEL_LATE_SPAN_HAT / KERNEL_STEP_HAT).floor() as usize).min(last);
        let early_lump: f64 = kernel[early_index + 1..].iter().sum::<f64>() * KERNEL_STEP_HAT;
        let late_lump: f64 = kernel[late_index + 1..].iter().sum::<f64>() * KERNEL_STEP_HAT;
        let mut jet = NolleJet {
            spec,
            centre_over_bulk,
            w_tilde_m,
            w_over_b,
            tau_hat,
            kernel,
            early_index,
            late_index,
            early_lump,
            late_lump,
            realisation_error: 0.0,
            dt,
            xi: 0.0,
            ring: Vec::new(),
            head: 0,
        };
        jet.realisation_error = jet.measure_realisation_error();
        if !(jet.realisation_error <= MAX_REALISATION_ERROR) {
            return Err(
                "tap realisation of the jet transfer misses the closed form by more than 15 %",
            );
        }
        let late_hat = late_index as f64 * KERNEL_STEP_HAT;
        let max_delay = (tau_hat + late_hat) * w_tilde_m / (min_centre_speed_m_s * dt);
        if !max_delay.is_finite() || max_delay + 4.0 > MAX_RING_SAMPLES as f64 {
            return Err("jet history exceeds the retained-sample limit");
        }
        jet.ring = vec![0.0; max_delay.ceil() as usize + 4];
        Ok(jet)
    }

    /// Frequency response of the tap set (units of the carrier removed) against `exp(mu W~) - 1`:
    /// the worst `|E_taps - E| / |E|` over St 0.03-0.7.
    fn measure_realisation_error(&self) -> f64 {
        let early_hat = self.early_index as f64 * KERNEL_STEP_HAT;
        let late_hat = self.late_index as f64 * KERNEL_STEP_HAT;
        let mut worst = 0.0f64;
        for j in 0..24 {
            let st = 0.03 + (0.7 - 0.03) * (j as f64) / 23.0;
            let w = TAU * st;
            let full = (mu_b(w / self.w_over_b) * self.w_over_b).exp() - 1.0;
            let mut re = self.kernel[0] * KERNEL_STEP_HAT;
            let mut im = 0.0;
            for (i, g) in self.kernel.iter().enumerate().skip(1) {
                let t = i as f64 * KERNEL_STEP_HAT;
                let weight = g * KERNEL_STEP_HAT;
                if i <= self.late_index {
                    re += weight * (w * t).cos();
                    im -= weight * (w * t).sin();
                }
                if i <= self.early_index {
                    re += weight * (w * t).cos();
                    im += weight * (w * t).sin();
                }
            }
            re += self.early_lump * (w * early_hat).cos() + self.late_lump * (w * late_hat).cos();
            im += self.early_lump * (w * early_hat).sin() - self.late_lump * (w * late_hat).sin();
            let error = ((re - full).powi(2) + im * im).sqrt() / (1.0 + full).abs();
            worst = worst.max(error);
        }
        worst
    }

    /// Worst relative error of the tap realisation against the closed form (St 0.03-0.7).
    #[must_use]
    pub fn realisation_error(&self) -> f64 {
        self.realisation_error
    }

    /// Travel distance over Bickley width, `W~ / b`.
    #[must_use]
    pub fn w_over_b(&self) -> f64 {
        self.w_over_b
    }

    /// Jet transit `W~ / c_p` at a centre speed [s] (infinite at rest).
    #[must_use]
    pub fn transit_s(&self, centre_speed_m_s: f64) -> f64 {
        self.w_tilde_m / (self.spec.phase_velocity_over_u0 * centre_speed_m_s)
    }

    /// Centre speed over bulk speed, `U0 / u`, used by this transfer.
    #[must_use]
    pub fn centre_over_bulk(&self) -> f64 {
        self.centre_over_bulk
    }

    /// The constants in force.
    #[must_use]
    pub fn spec(&self) -> NolleSpec {
        self.spec
    }

    /// Accepted leaky integral of the transverse disturbance [m]. This is
    /// physical source memory, not the labium displacement after convection.
    #[must_use]
    pub fn receptivity_displacement_m(&self) -> f64 {
        self.xi
    }

    fn push(&mut self, value: f64) {
        let len = self.ring.len();
        self.ring[self.head] = value;
        self.head = (self.head + 1) % len;
    }

    // Newest entry (delay 0) is the one just pushed. Requires 0 <= delay <= len - 2.
    fn at(&self, delay_samples: f64) -> f64 {
        let len = self.ring.len();
        let whole = delay_samples.floor();
        let fraction = delay_samples - whole;
        let k = whole as usize;
        let sample = |back: usize| self.ring[(self.head + len - 1 - back) % len];
        (1.0 - fraction) * sample(k) + fraction * sample(k + 1)
    }

    // The ring still ends at n-1. Separate the unknown xi_n from the actual
    // retarded interpolation; a sub-sample transit must not acquire a full delay.
    fn pending_at(&self, delay_samples: f64) -> (f64, f64) {
        if delay_samples < 1.0 {
            (1.0 - delay_samples, delay_samples * self.at(0.0))
        } else {
            (0.0, self.at(delay_samples - 1.0))
        }
    }

    fn kernel_at(&self, t_hat: f64) -> f64 {
        let position = t_hat.abs() / KERNEL_STEP_HAT;
        let i = position.floor() as usize;
        if i + 1 >= self.kernel.len() {
            return 0.0;
        }
        let fraction = position - i as f64;
        (1.0 - fraction) * self.kernel[i] + fraction * self.kernel[i + 1]
    }

    /// Condense one held-input step without changing source memory. Every
    /// delayed tap is affine in the same interval's endpoint displacement.
    pub(crate) fn prepare_step(
        &self,
        centre_speed_m_s: f64,
        flue_width_m: f64,
        seed_velocity_m_s: f64,
    ) -> Result<NolleStep, &'static str> {
        if !(centre_speed_m_s.is_finite()
            && centre_speed_m_s >= 0.0
            && flue_width_m.is_finite()
            && flue_width_m > 0.0
            && seed_velocity_m_s.is_finite())
        {
            return Err("nolle jet step needs finite inputs and a positive flue width");
        }
        if centre_speed_m_s == 0.0 {
            return Ok(NolleStep {
                xi_history_m: 0.0,
                xi_mouth_s_m2: 0.0,
                xi_jet_s_m2: 0.0,
                eta_history_m: 0.0,
                eta_mouth_s_m2: 0.0,
                eta_jet_s_m2: 0.0,
            });
        }
        let dt = self.dt;
        let w = self.w_tilde_m;
        let dt_hat = centre_speed_m_s * dt / w;
        let leak = self.spec.leak_rel * centre_speed_m_s / w;
        // xi' + leak*xi = disturbance. exp_m1 retains the source integral when
        // leak*dt is small; this removes Euler's decay error and stability limit.
        let decay = (-leak * dt).exp();
        let input_time = -(-leak * dt).exp_m1() / leak;
        let xi_history_m = decay * self.xi + input_time * seed_velocity_m_s;
        let xi_mouth_s_m2 = input_time * self.spec.alpha_p / (flue_width_m * w);
        let xi_jet_s_m2 = -input_time * self.spec.alpha_j / (flue_width_m * w);
        let carrier_delay = self.tau_hat / dt_hat;
        let spacing_hat = dt_hat.max(TAP_SPACING_HAT);
        let early_hat = self.early_index as f64 * KERNEL_STEP_HAT;
        let late_hat = self.late_index as f64 * KERNEL_STEP_HAT;
        let latest = carrier_delay + late_hat / dt_hat;
        if !dt_hat.is_finite()
            || dt_hat <= 0.0
            || !latest.is_finite()
            || !(latest + 2.0 <= self.ring.len() as f64)
        {
            return Err("jet transit exceeds the retained history");
        }
        let (carrier_direct, carrier_history) = self.pending_at(carrier_delay);
        let mut direct = (1.0 + self.kernel[0] * spacing_hat) * carrier_direct;
        let mut history = (1.0 + self.kernel[0] * spacing_hat) * carrier_history;
        let mut accumulate = |weight: f64, delay: f64| {
            let (current, past) = self.pending_at(delay);
            direct += weight * current;
            history += weight * past;
        };
        let late_taps = (late_hat / spacing_hat).floor() as usize;
        for k in 1..=late_taps {
            let offset_hat = k as f64 * spacing_hat;
            accumulate(
                self.kernel_at(offset_hat) * spacing_hat,
                carrier_delay + offset_hat / dt_hat,
            );
        }
        let early_taps = (early_hat / spacing_hat).floor() as usize;
        for k in 1..=early_taps {
            let offset_hat = k as f64 * spacing_hat;
            accumulate(
                self.kernel_at(offset_hat) * spacing_hat,
                (carrier_delay - offset_hat / dt_hat).max(0.0),
            );
        }
        accumulate(
            self.early_lump,
            (carrier_delay - early_hat / dt_hat).max(0.0),
        );
        accumulate(self.late_lump, carrier_delay + late_hat / dt_hat);
        let gain = self.spec.jet_gain * self.centre_over_bulk;
        let response = NolleStep {
            xi_history_m,
            xi_mouth_s_m2,
            xi_jet_s_m2,
            eta_history_m: gain * ((1.0 - direct) * xi_history_m - history),
            eta_mouth_s_m2: gain * (1.0 - direct) * xi_mouth_s_m2,
            eta_jet_s_m2: gain * (1.0 - direct) * xi_jet_s_m2,
        };
        if ![
            decay,
            input_time,
            xi_history_m,
            xi_mouth_s_m2,
            xi_jet_s_m2,
            response.eta_history_m,
            response.eta_mouth_s_m2,
            response.eta_jet_s_m2,
        ]
        .iter()
        .all(|v| v.is_finite())
            || input_time <= 0.0
        {
            return Err("nolle jet left the finite domain");
        }
        Ok(response)
    }

    // Called only after the complete physical host has accepted this step.
    pub(crate) fn commit(&mut self, candidate: NolleCandidate) {
        self.xi = candidate.xi_m;
        self.push(candidate.xi_m);
    }

    /// Advance one sample with explicitly held mouth and fluctuating jet flows,
    /// returning labium displacement `eta` [m], positive towards the pipe.
    /// A coupled host instead uses the affine preview to solve these flows at
    /// the same time level. No source memory is changed on a refused step.
    /// At zero centre speed there is no jet and no response.
    pub fn step(
        &mut self,
        centre_speed_m_s: f64,
        mouth_flow_m3_s: f64,
        jet_flow_m3_s: f64,
        flue_width_m: f64,
        seed_velocity_m_s: f64,
    ) -> Result<f64, &'static str> {
        let candidate = self
            .prepare_step(centre_speed_m_s, flue_width_m, seed_velocity_m_s)?
            .candidate(mouth_flow_m3_s, jet_flow_m3_s)?;
        self.commit(candidate);
        Ok(candidate.eta_m)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const PRICE_STATED_EDGES: [f64; 4] = [0.0705, 0.2871, 0.4964, 0.7690];
    const PRICE_MEASURED_EDGES: [f64; 4] = [0.054, 0.267, 0.43, 0.60];
    // Price head: W~ = 4.83 mm + 0.7 mm chamfer, h = 1.09 mm, b = 2h/5, W~/b = 12.68.
    const PRICE_W_OVER_B: f64 = 5.53e-3 / (0.4 * 1.09e-3);

    fn first_four(edges: Vec<f64>) -> [f64; 4] {
        assert!(edges.len() >= 4, "{edges:?}");
        [edges[0], edges[1], edges[2], edges[3]]
    }

    #[test]
    fn theory_table_reproduces_prices_stated_model_edges() {
        // G2: the by-eye digitisation of Fig. 5 reproduces the four edges Price prints for his model.
        let edges = first_four(gain_band_edges(
            PRICE_W_OVER_B,
            PhaseVelocity::TheoryTable,
            0.002,
            1.0,
            8000,
        ));
        for (got, stated) in edges.iter().zip(PRICE_STATED_EDGES.iter()) {
            assert!(
                (got / stated - 1.0).abs() < 0.04,
                "{edges:?} vs {PRICE_STATED_EDGES:?}"
            );
        }
    }

    #[test]
    fn theory_table_fails_the_25_percent_gate_against_the_measured_bands() {
        // The pre-registered gate (every edge within 25 % of the measured St bands) FAILS for the
        // theory phase velocity: this test pins that finding. Do not "fix" it by editing the tables.
        let edges = first_four(gain_band_edges(
            PRICE_W_OVER_B,
            PhaseVelocity::TheoryTable,
            0.002,
            1.0,
            8000,
        ));
        let misses = edges
            .iter()
            .zip(PRICE_MEASURED_EDGES.iter())
            .filter(|(got, measured)| (*got / *measured - 1.0).abs() > 0.25)
            .count();
        assert!(misses >= 2, "{edges:?}");
    }

    #[test]
    fn measured_long_channel_phase_velocity_is_consistent_with_the_measured_bands() {
        // CONSISTENCY ONLY: c_p = 0.3 U0 was adopted after the theory table failed (not independent).
        let edges = first_four(gain_band_edges(
            PRICE_W_OVER_B,
            PhaseVelocity::Constant(LONG_CHANNEL_CP_OVER_U0),
            0.002,
            1.0,
            8000,
        ));
        for (got, measured) in edges.iter().zip(PRICE_MEASURED_EDGES.iter()) {
            assert!((got / measured - 1.0).abs() < 0.25, "{edges:?}");
        }
    }

    fn lock_in(
        jet: &mut NolleJet,
        speed: f64,
        width: f64,
        flow_amp: f64,
        period: usize,
    ) -> (f64, f64) {
        let dt = jet.dt;
        let omega = TAU / (period as f64 * dt);
        let settle = 30 * period;
        let measure = 10 * period;
        let (mut a, mut b) = (0.0, 0.0);
        for n in 1..=settle + measure {
            let flow_prev = flow_amp * (omega * (n as f64 - 1.0) * dt).cos();
            let eta = jet.step(speed, flow_prev, 0.0, width, 0.0).expect("step");
            if n > settle {
                a += eta * (omega * n as f64 * dt).cos();
                b += eta * (omega * n as f64 * dt).sin();
            }
        }
        (2.0 * a / measure as f64, 2.0 * b / measure as f64)
    }

    #[test]
    fn step_matches_the_closed_form_in_the_time_domain() {
        // G2: lock-in of the sample recursion against eta = -alpha_p J Q / (u H) (alpha_j = 0),
        // across the first gain band and above, long-channel windway geometry (Price head).
        let mut spec = NolleSpec::price_2015(ChannelClass::Long, 0.7e-3).unwrap();
        spec.alpha_j = 0.0;
        let dt = 1.0 / 288_000.0;
        let (h, cut_up, width) = (1.09e-3, 4.83e-3, 14.7e-3);
        let u_bulk = 8.7;
        let u0 = POISEUILLE_CENTRE_OVER_BULK * u_bulk;
        for period in [1200_usize, 600, 300] {
            let mut jet =
                NolleJet::new(spec, POISEUILLE_CENTRE_OVER_BULK, h, cut_up, dt, 0.5 * u0).unwrap();
            let q0 = 1.0e-6;
            let (a, b) = lock_in(&mut jet, u0, width, q0, period);
            let omega = TAU / (period as f64 * dt);
            let omega_hat = omega * (cut_up + 0.7e-3) / u0;
            let (jr, ji) = jet_transfer_closed_form(
                omega_hat,
                jet.w_over_b(),
                PhaseVelocity::Constant(LONG_CHANNEL_CP_OVER_U0),
                PRICE_JET_GAIN,
            );
            let scale = -ALPHA_P * q0 / (u_bulk * width);
            let (er, ei) = (scale * jr, scale * ji);
            // eta = A cos + B sin  <=>  phasor A - iB.
            let (mr, mi) = (a, -b);
            let error = ((mr - er).powi(2) + (mi - ei).powi(2)).sqrt() / (er * er + ei * ei).sqrt();
            assert!(
                error < 0.15,
                "period {period}: got ({mr:e}, {mi:e}) expected ({er:e}, {ei:e}) error {error}"
            );
        }
    }

    #[test]
    fn construction_reports_a_small_realisation_error_for_both_channel_classes() {
        let dt = 1.0 / 288_000.0;
        let long = NolleJet::new(
            NolleSpec::price_2015(ChannelClass::Long, 0.0).unwrap(),
            POISEUILLE_CENTRE_OVER_BULK,
            1.48e-3,
            5.3e-3,
            dt,
            0.3,
        )
        .unwrap();
        let short = NolleJet::new(
            NolleSpec::price_2015(ChannelClass::Short, 0.0).unwrap(),
            1.0,
            1.0e-3,
            4.0e-3,
            dt,
            1.0,
        )
        .unwrap();
        assert!(long.realisation_error() < MAX_REALISATION_ERROR);
        assert!(short.realisation_error() < MAX_REALISATION_ERROR);
    }

    #[test]
    fn rest_gives_no_response_and_bad_input_is_refused() {
        let dt = 1.0 / 288_000.0;
        let mut jet = NolleJet::new(
            NolleSpec::price_2015(ChannelClass::Long, 0.0).unwrap(),
            1.0,
            1.0e-3,
            5.0e-3,
            dt,
            1.0,
        )
        .unwrap();
        assert_eq!(jet.step(0.0, 1.0e-6, 0.0, 8.0e-3, 0.0).unwrap(), 0.0);
        assert!(jet.step(f64::NAN, 0.0, 0.0, 8.0e-3, 0.0).is_err());
        assert!(jet.step(10.0, 0.0, 0.0, 0.0, 0.0).is_err());
        assert!(NolleSpec::price_2015(ChannelClass::Long, -1.0).is_err());
    }

    #[test]
    fn g1_held_receptivity_is_exact_on_three_source_clocks() {
        // Independent forced first-order ODE solution, including both physical
        // feedback flows and the seed. No gain, material or pole is retuned.
        let (speed, width, mouth, direct_jet, seed) = (15.0, 0.01, 2e-6, -3e-7, 1e-5);
        let duration = 0.001;
        for rate in [48_000, 96_000, 192_000] {
            let spec = NolleSpec::price_2015(ChannelClass::Long, 0.0007).unwrap();
            let mut jet =
                NolleJet::new(spec, 1.5, 0.00109, 0.00483, 1.0 / f64::from(rate), 1.0).unwrap();
            jet.xi = 2e-7;
            let leak = spec.leak_rel * speed / jet.w_tilde_m;
            let equilibrium = ((spec.alpha_p * mouth - spec.alpha_j * direct_jet)
                / (width * jet.w_tilde_m)
                + seed)
                / leak;
            let expected = equilibrium + (jet.xi - equilibrium) * (-leak * duration).exp();
            for _ in 0..rate / 1000 {
                jet.step(speed, mouth, direct_jet, width, seed).unwrap();
            }
            assert!(
                (jet.receptivity_displacement_m() - expected).abs()
                    < 512.0 * f64::EPSILON * expected.abs(),
                "rate {rate}"
            );
        }
    }

    #[test]
    fn g1_pending_retarded_trace_preserves_affine_time_and_current_sample() {
        let mut jet = NolleJet::new(
            NolleSpec::price_2015(ChannelClass::Short, 0.0).unwrap(),
            1.0,
            0.001,
            0.004,
            1.0 / 48_000.0,
            1.0,
        )
        .unwrap();
        for head in [0, 3, jet.ring.len() - 1] {
            jet.head = head;
            for back in 0..jet.ring.len() {
                let slot = (head + jet.ring.len() - 1 - back) % jet.ring.len();
                jet.ring[slot] = 2.0 - 0.7 * (back + 1) as f64;
            }
            for delay in [0.0, 0.25, 0.75, 1.0, 1.25, 3.7] {
                let (current, history) = jet.pending_at(delay);
                assert!((current * 2.0 + history - (2.0 - 0.7 * delay)).abs() < 2e-15);
            }
        }
    }

    #[test]
    fn g0_refused_steps_and_affine_trials_preserve_the_complete_source_history() {
        let mut jet = NolleJet::new(
            NolleSpec::price_2015(ChannelClass::Long, 0.0).unwrap(),
            1.0,
            0.001,
            0.004,
            1.0 / 48_000.0,
            2.0,
        )
        .unwrap();
        for _ in 0..7 {
            jet.step(20.0, 1e-6, 2e-7, 0.01, 0.0).unwrap();
        }
        let saved = jet.clone();
        let step = jet.prepare_step(20.0, 0.01, 0.0).unwrap();
        for q in [-1e-5, 0.0, 1e-5] {
            step.candidate(q, -q).unwrap();
        }
        assert!(jet.step(0.01, 1e-6, 0.0, 0.01, 0.0).is_err());
        assert!(jet.step(20.0, f64::NAN, 0.0, 0.01, 0.0).is_err());
        assert_eq!(jet.xi, saved.xi);
        assert_eq!(jet.head, saved.head);
        assert_eq!(jet.ring, saved.ring);
    }
}
