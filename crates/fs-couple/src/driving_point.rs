//! Characteristic-port composition: a viscothermal TMM driving-point
//! becomes a [`fs_vfit::discretize::DelayedFilter`].
//!
//! This file is the coupling-layer filling of that nameless port. A
//! bore, a muffler, an HVAC run, and a pulse tube are the same object.
/// Static regional media composed at the existing TMM pressure/flow ports.
pub mod regional;
/// Explicitly terminated generic acoustic-network frequency-domain ports.
pub mod network;
/// Material molecular reservoirs in the same characteristic-port realization.
pub mod molecular;
/// Mouth volume flow of a duct from the wave incident on its inlet.
pub mod outlet;
mod baked;

use fs_duct::compose::input_reflection_wall;
use fs_duct::{
    Duct, DuctError, DuctResponse, HoleState, LossModel, MouthLoad, Segment, Termination,
    input_impedance_wall, segment_wave,
};
use fs_fft::{C64 as FftC64, Fft};
use fs_material::gas::GasState;
use fs_math::c64::C64;
use fs_phs::WallPin;
use fs_vfit::discretize::{DelayedFilter, DiscretizeError, RealizeError, reflectance};

/// Passivity tolerance on the 4x-oversampled reflectance grid (`|R(e^{iωΔt})| <= 1 + 1e-12`).
pub const PASSIVITY_TOL: f64 = 1.0e-12;

/// The static reflection is supplied separately from positive-frequency bins.
/// Tables have an explicit DC limit; they never extrapolate or splice a fit.
#[derive(Clone, Copy)]
pub(crate) enum LineDc {
    Analytic(Termination),
    Explicit(f64),
}
impl From<Termination> for LineDc {
    fn from(value: Termination) -> Self { Self::Analytic(value) }
}

/// Convergence ladder specification for [`realize_line`] and [`realize_impedance_line`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LineLadder {
    /// Relative tail energy ceiling in `[3N/4 .. N)` (`default = 1.0e-10`).
    pub tail_tol: f64,
    /// Maximum resonance-peak shift between adjacent rungs (`default = 0.05` cents).
    pub peak_tol_cents: f64,
    /// Hard upper bound on DFT length (`default = 65_536`, at most `524_288`).
    pub max_taps: usize,
}

/// Largest ceiling a caller may declare for a line ladder; the default stays lower.
const MAX_LINE_TAPS: usize = 524_288;

impl Default for LineLadder {
    fn default() -> Self {
        Self {
            tail_tol: 1.0e-10,
            peak_tol_cents: 0.05,
            max_taps: 65_536,
        }
    }
}

/// Per-rung diagnostic recorded during [`realize_line`] convergence.
#[derive(Debug, Clone, PartialEq)]
pub struct LineRung {
    /// Realized DFT length at this rung.
    pub n_fft: usize,
    /// Relative tail energy in `[3N/4 .. N)` over total tap energy.
    pub tail_energy: f64,
    /// Worst shift in cents across the first three resonances vs the prior rung (`+∞` at rung 0).
    pub peak_shift_cents: f64,
}

/// Explicit time-domain pole or delay state peeled from or accompanying a line realization.
#[derive(Debug, Clone, PartialEq)]
pub enum PoleState {
    /// `ω = 0` integrator for closed volumes (`p_comp[n+1] = p_comp[n] + gain * dt * u[n]`),
    /// peeled before the DFT so `Z(0)` is not zeroed.
    Compliance {
        /// Static cavity stiffness [Pa/m³]: `ρ c² / V` for lossless gas,
        /// `ρ c² / (γ V)` for the Bessel isothermal low-frequency limit.
        gain: f64,
    },
    /// An analytically established pole stepped by exact ZOH
    /// `z_{n+1} = exp((-sigma + i*omega)*dt) * z_n + ...` so ring-down length is independent of `n_fft`.
    Resonance {
        /// Angular frequency `ω` [rad/s].
        omega: f64,
        /// Decay rate `σ` [1/s].
        sigma: f64,
        /// Complex modal residue `r` under `e^{-iωt}`.
        residue: C64,
    },
    /// Bulk round-trip propagation delay extracted from geometry.
    BulkDelay {
        /// Integer delay samples.
        samples: usize,
        /// Fractional Thiran allpass delay in `[0, 1)`.
        frac: f64,
    },
}

/// Receipt of a converged driving-point line realization.
#[derive(Debug, Clone, PartialEq)]
pub struct LineRealization {
    /// Converged DFT tap count.
    pub n_fft: usize,
    /// Relative tail energy in `[3N/4 .. N)` at the converged rung.
    pub tail_energy: f64,
    /// Peak shift in cents at the converged rung vs the preceding rung.
    pub peak_shift_cents: f64,
    /// Explicit pole and bulk-delay states attached to the realization.
    pub pole_states: Vec<PoleState>,
    /// Rungs of the accepted tap construction. An untapered fallback records
    /// its own two measured rungs, never a mixed tapered/untapered comparison.
    pub rungs_tried: Vec<LineRung>,
}

/// Typed refusal from a driving-point characteristic line.
#[derive(Debug, Clone, PartialEq)]
pub enum DrivingPointError {
    /// Geometry or delay is not realizable.
    Invalid {
        /// Which check failed.
        what: &'static str,
    },
    /// Convergence ladder hit `max_taps` before `tail_energy` and `peak_shift_cents` converged.
    LadderUnconverged {
        /// Ceiling tap count reached.
        max_taps: usize,
        /// Measured relative tail energy in `[3N/4 .. N)` at the ceiling rung.
        tail_energy: f64,
        /// Required tail energy tolerance.
        tail_tol: f64,
        /// Measured resonance peak shift [cents] at the ceiling rung.
        peak_shift_cents: f64,
        /// Required peak shift tolerance [cents].
        peak_tol_cents: f64,
    },
    /// A transfer ladder hit `max_taps` before `tail_energy` and `tap_change` converged. A transfer has no peak-shift
    /// gate: it gates on the change of its shared taps.
    TransferLadderUnconverged {
        /// Ceiling tap count reached.
        max_taps: usize,
        /// Measured relative tail energy in `[3N/4 .. N)` at the ceiling rung.
        tail_energy: f64,
        /// Required tail energy tolerance.
        tail_tol: f64,
        /// Measured RMS change of the taps shared with the rung below, over the ceiling rung's RMS (`+∞` at rung 0).
        tap_change: f64,
        /// Required tap change tolerance: the square root of `tail_tol`.
        tap_tol: f64,
    },
    /// TMM refusal.
    Duct(DuctError),
    /// Original constitutive-state refusal from a material-backed line.
    Material(fs_material::MaterialError),
    /// Vector-fit / bilinear refusal.
    Realize(RealizeError),
    /// Discrete runtime refusal.
    Discrete(DiscretizeError),
}

impl core::fmt::Display for DrivingPointError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid { what } => write!(f, "FS-COUPLE-DRIVE: {what}"),
            Self::LadderUnconverged {
                max_taps,
                tail_energy,
                tail_tol,
                peak_shift_cents,
                peak_tol_cents,
            } => write!(
                f,
                "FS-COUPLE-DRIVE-UNCONVERGED: ladder reached max_taps={max_taps} with \
                 tail_energy={tail_energy:.3e} (tol={tail_tol:.3e}) and \
                 peak_shift_cents={peak_shift_cents:.4} (tol={peak_tol_cents:.4})"
            ),
            Self::TransferLadderUnconverged {
                max_taps,
                tail_energy,
                tail_tol,
                tap_change,
                tap_tol,
            } => write!(
                f,
                "FS-COUPLE-DRIVE-UNCONVERGED: ladder reached max_taps={max_taps} with \
                 tail_energy={tail_energy:.3e} (tol={tail_tol:.3e}) and \
                 tap_change={tap_change:.3e} (tol={tap_tol:.3e})"
            ),
            Self::Duct(e) => write!(f, "{e}"),
            Self::Material(e) => write!(f, "{e}"),
            Self::Realize(e) => write!(f, "{e}"),
            Self::Discrete(e) => write!(f, "{e}"),
        }
    }
}

impl std::error::Error for DrivingPointError {}

impl From<DuctError> for DrivingPointError {
    fn from(error: DuctError) -> Self { Self::Duct(error) }
}

impl From<fs_duct::molecular::MolecularDuctError> for DrivingPointError {
    fn from(error: fs_duct::molecular::MolecularDuctError) -> Self {
        match error {
            fs_duct::molecular::MolecularDuctError::Duct(error) => Self::Duct(error),
            fs_duct::molecular::MolecularDuctError::Material(error) => Self::Material(error),
        }
    }
}

fn duct_axial_length(physics: &Duct) -> f64 {
    physics
        .segments
        .iter()
        .map(Segment::axial_length)
        .sum()
}

/// Sealed-volume geometric compliance gain `g_comp = ρ c² / V` [Pa / m³] when the
/// duct has a rigid closed termination (`Termination::Closed`), rigid walls
/// (`wall.is_none()`), and no open or vented tone holes. Returns `None` for any
/// open, vented, or yielding-wall duct.
#[must_use]
pub fn sealed_volume_compliance_gain(
    physics: &Duct,
    gas: &GasState,
    termination: Termination,
    wall: Option<&WallPin>,
) -> Option<f64> {
    if !matches!(termination, Termination::Closed) || wall.is_some() {
        return None;
    }
    let mut volume = 0.0_f64;
    for seg in &physics.segments {
        match *seg {
            // The optional closed-volume shortcut has no profile-volume
            // quadrature. Full frequency-domain propagation remains available.
            Segment::Profile { .. } => return None,
            Segment::Cylinder { radius, length } => {
                volume += core::f64::consts::PI * radius * radius * length;
            }
            Segment::Cone {
                inlet_radius: r1,
                outlet_radius: r2,
                length,
            } => {
                volume += (core::f64::consts::PI * length / 3.0) * (r1 * r1 + r1 * r2 + r2 * r2);
            }
            Segment::ToneHole {
                hole_radius,
                chimney_height,
                state,
                ..
            } => match state {
                HoleState::Closed => {
                    volume += core::f64::consts::PI * hole_radius * hole_radius * chimney_height;
                }
                HoleState::Open => return None,
                HoleState::Vent(s) => {
                    if s > 0.0 {
                        return None;
                    }
                    volume += core::f64::consts::PI * hole_radius * hole_radius * chimney_height;
                }
            },
        }
    }
    if volume > 0.0 && volume.is_finite() {
        let bulk_modulus = gas.density * gas.sound_speed * gas.sound_speed;
        let gain = bulk_modulus / volume;
        gain.is_finite().then_some(gain)
    } else {
        None
    }
}

fn bulk_delay_pole_state(geo_delay: f64) -> PoleState {
    let samples = geo_delay.floor().max(0.0) as usize;
    let frac = (geo_delay - samples as f64).clamp(0.0, 1.0 - f64::EPSILON);
    PoleState::BulkDelay { samples, frac }
}

fn validate_ladder(ladder: &LineLadder) -> Result<usize, DrivingPointError> {
    if !(ladder.tail_tol > 0.0
        && ladder.tail_tol.is_finite()
        && ladder.peak_tol_cents > 0.0
        && ladder.peak_tol_cents.is_finite())
    {
        return Err(DrivingPointError::Invalid {
            what: "line ladder tolerances must be positive and finite",
        });
    }
    if !(256..=MAX_LINE_TAPS).contains(&ladder.max_taps) || !ladder.max_taps.is_power_of_two() {
        return Err(DrivingPointError::Invalid {
            what: "line ladder max_taps must be a power of two in [256, 524288]",
        });
    }
    Ok(ladder.max_taps)
}

/// Evaluate `|Z(ω)/Zc|² = |1 + R(ω)|² / |1 - R(ω)|²` from the continuous DTFT
/// `R(ω) = Σ_{n=0}^{N-1} h[n] e^{+i ω n Δt}` (`e^{-iωt}` convention).
fn dtft_norm_impedance_sq(ir: &[f64], dt: f64, omega: f64) -> f64 {
    let theta = omega * dt;
    let (sin_d, cos_d) = theta.sin_cos();
    let mut cos_n = 1.0_f64;
    let mut sin_n = 0.0_f64;
    let mut re = 0.0_f64;
    let mut im = 0.0_f64;
    for (idx, &h) in ir.iter().enumerate() {
        if idx != 0 && (idx & 63) == 0 {
            let (s, c) = (theta * idx as f64).sin_cos();
            cos_n = c;
            sin_n = s;
        }
        re += h * cos_n;
        im += h * sin_n;
        let next_cos = cos_n * cos_d - sin_n * sin_d;
        let next_sin = sin_n * cos_d + cos_n * sin_d;
        cos_n = next_cos;
        sin_n = next_sin;
    }
    let num = (1.0 + re) * (1.0 + re) + im * im;
    let den = ((1.0 - re) * (1.0 - re) + im * im).max(1.0e-30);
    num / den
}

/// Existing phase-tapered minimum-phase C^3 Nyquist rolloff on a positive
/// frequency grid. Computing it once on a ladder's ceiling and subsampling
/// even bins keeps the numerical transfer function fixed between rungs.
/// This window is a realization treatment, not a passivity certificate.
pub(crate) fn reflectance_rolloff(n_fft: usize) -> Vec<C64> {
    let half = n_fft / 2;
    let fft = Fft::new(n_fft);
    let mut scratch = vec![FftC64::new(0.0, 0.0); n_fft];

    // Phase-tapered minimum-phase C^3 Nyquist rolloff W(θ) = exp(-α(θ) - i η(θ) H[α](θ))
    // with α(u) = 36 t^4 on u >= 0.97 and a C^3 septic smoothstep η(u) on u ∈ [0.50, 0.95].
    // For u <= 0.50, W(θ) = 1 + 0i identically (zero low-frequency phase distortion);
    // for u >= 0.95, η(u) = 1 so W(θ) matches the causal minimum-phase rolloff near Nyquist.
    let mut alpha = vec![0.0_f64; n_fft];
    let u0 = 0.97_f64;
    for k in 1..=half {
        let u = k as f64 / half as f64;
        if u > u0 {
            let t = (u - u0) / (1.0 - u0);
            let a = 36.0 * t * t * t * t;
            alpha[k] = a;
            if k != half {
                alpha[n_fft - k] = a;
            }
        }
    }
    let mut cep = alpha
        .iter()
        .map(|&a| FftC64::new(-a, 0.0))
        .collect::<Vec<_>>();
    fft.inverse(&mut cep, &mut scratch);
    cep[0].im = 0.0;
    for item in cep.iter_mut().take(half).skip(1) {
        item.re *= 2.0;
        item.im *= 2.0;
    }
    if n_fft.is_multiple_of(2) {
        cep[half].im = 0.0;
    }
    for item in cep.iter_mut().skip(half + 1) {
        item.re = 0.0;
        item.im = 0.0;
    }
    fft.forward(&mut cep, &mut scratch);

    let mut rolloff = vec![C64::ONE; half + 1];
    for k in 1..=half {
        let u = k as f64 / half as f64;
        let eta = if u <= 0.50 {
            0.0
        } else if u >= 0.95 {
            1.0
        } else {
            let s = (u - 0.50) / (0.95 - 0.50);
            let s4 = s * s * s * s;
            s4 * (35.0 - 84.0 * s + 70.0 * s * s - 20.0 * s * s * s)
        };
        let mag = cep[k].re.exp();
        let (s, c) = (eta * cep[k].im).sin_cos();
        rolloff[k] = C64::new(mag * c, mag * s);
    }
    rolloff
}

/// Realize a causal reflectance impulse response from positive-frequency samples
/// `r_pos[0..=n_fft/2]` (`e^{-iωt}` convention) using a phase-tapered minimum-phase
/// $C^3$ Nyquist rolloff `W(θ) = exp(-α(θ) - i η(θ) H[α](θ))` above `u0 = 0.97`.
///
/// Existing callers retain the same rung-local rolloff and arithmetic. The
/// complete-matrix ladder explicitly supplies a fixed reference through
/// [`realize_reflectance_ir_with_rolloff`] instead.
pub(crate) fn realize_reflectance_ir_from_positive_spectrum(
    r_pos: &[C64],
    n_fft: usize,
    _geo_delay: f64,
    fold_precursor: bool,
) -> Vec<f64> {
    let rolloff = reflectance_rolloff(n_fft);
    realize_reflectance_ir_with_rolloff(r_pos, &rolloff, n_fft, fold_precursor)
}

/// The same inverse transform and precursor treatment with an explicit shared
/// rolloff. Callers supply one coefficient per positive-frequency sample and
/// remain responsible for final causal convergence and passivity admission.
pub(crate) fn realize_reflectance_ir_with_rolloff(
    r_pos: &[C64],
    rolloff: &[C64],
    n_fft: usize,
    fold_precursor: bool,
) -> Vec<f64> {
    let half = n_fft / 2;
    debug_assert_eq!(r_pos.len(), half + 1);
    debug_assert_eq!(rolloff.len(), half + 1);
    let fft = Fft::new(n_fft);
    let mut scratch = vec![FftC64::new(0.0, 0.0); n_fft];
    let mut buf = vec![FftC64::new(0.0, 0.0); n_fft];
    buf[0] = FftC64::new(r_pos[0].re, 0.0);
    for k in 1..=half {
        let rac = r_pos[k];
        // Acoustic e^{-iωt} → DFT e^{+iωt} so the IR is causal (DEF-01).
        let r_dft = C64::new(rac.re, -rac.im);
        let rw = r_dft * rolloff[k];
        buf[k] = FftC64::new(rw.re, rw.im);
        if k != half {
            buf[n_fft - k] = FftC64::new(rw.re, -rw.im);
        } else {
            buf[half] = FftC64::new(rw.re, 0.0);
        }
    }
    fft.inverse(&mut buf, &mut scratch);
    let mut ir: Vec<f64> = buf.iter().map(|c| c.re).collect();
    if fold_precursor {
        // Fold the localized precursor skirt (from tapering η on [0.50, 0.95])
        // at the end of the circular buffer back into its causal mirror.
        let pre_k = (n_fft / 16).min(192);
        for m in 1..=pre_k {
            ir[m] += ir[n_fft - m];
            ir[n_fft - m] = 0.0;
        }
    }
    ir
}

/// Locate the first three impedance resonance peaks `[Hz]` and maximum 4x-oversampled
/// reflectance magnitude `max |R(e^{iωΔt})|` from a realized impulse response `ir`.
pub(crate) fn analyze_reflectance_ir(
    ir: &[f64],
    dt: f64,
    geo_delay: f64,
) -> (f64, Vec<f64>) {
    let n_fft = ir.len();
    let n_os = 4 * n_fft;
    let fft_os = Fft::new(n_os);
    let mut os_buf = vec![FftC64::new(0.0, 0.0); n_os];
    for (dst, &src) in os_buf[..n_fft].iter_mut().zip(ir) {
        *dst = FftC64::new(src, 0.0);
    }
    let mut os_scratch = vec![FftC64::new(0.0, 0.0); n_os];
    fft_os.forward(&mut os_buf, &mut os_scratch);

    let half_os = n_os / 2;
    let mut max_r_abs = 0.0_f64;
    let mut norm_z_sq = vec![0.0_f64; half_os + 1];
    for (m, c) in os_buf[..=half_os].iter().enumerate() {
        let r_abs = (c.re * c.re + c.im * c.im).sqrt();
        max_r_abs = max_r_abs.max(r_abs);
        // Under e^{-iωt}, R(ω_m) = conj(os_buf[m]); |1 ± R|² is invariant to conjugation.
        let num = (1.0 + c.re) * (1.0 + c.re) + c.im * c.im;
        let den = ((1.0 - c.re) * (1.0 - c.re) + c.im * c.im).max(1.0e-30);
        norm_z_sq[m] = num / den;
    }

    let mode_spacing = n_os as f64 / geo_delay.max(2.0);
    let w = ((0.20 * mode_spacing).round() as usize).clamp(2, 64);
    let m_min = ((0.30 * mode_spacing).floor() as usize)
        .max(w + 1)
        .min(half_os / 2);
    let m_max = ((0.90 * half_os as f64) as usize)
        .max(m_min + w)
        .min(half_os.saturating_sub(w + 1));
    let d_omega = core::f64::consts::TAU / (n_os as f64 * dt);

    let collect_peaks = |min_prominence_ratio: f64| -> Vec<f64> {
        let mut peaks = Vec::with_capacity(3);
        let mut m = m_min;
        while m <= m_max && peaks.len() < 3 {
            let g_cur = norm_z_sq[m];
            if g_cur > norm_z_sq[m - 1] && g_cur >= norm_z_sq[m + 1] {
                let is_local_max = norm_z_sq[m - w..=m + w].iter().all(|&v| g_cur >= v);
                let flank = norm_z_sq[m - w].min(norm_z_sq[m + w]);
                if is_local_max && g_cur > min_prominence_ratio * flank {
                    // Golden-section search on continuous DTFT inside [ω_{m-1}, ω_{m+1}].
                    let mut a = (m - 1) as f64 * d_omega;
                    let mut b = (m + 1) as f64 * d_omega;
                    let inv_phi = 0.5 * (5.0_f64.sqrt() - 1.0);
                    let mut c = b - inv_phi * (b - a);
                    let mut d = a + inv_phi * (b - a);
                    let mut fc = dtft_norm_impedance_sq(ir, dt, c);
                    let mut fd = dtft_norm_impedance_sq(ir, dt, d);
                    for _ in 0..16 {
                        if fc >= fd {
                            b = d;
                            d = c;
                            fd = fc;
                            c = b - inv_phi * (b - a);
                            fc = dtft_norm_impedance_sq(ir, dt, c);
                        } else {
                            a = c;
                            c = d;
                            fc = fd;
                            d = a + inv_phi * (b - a);
                            fd = dtft_norm_impedance_sq(ir, dt, d);
                        }
                    }
                    let omega_peak = 0.5 * (a + b);
                    peaks.push(omega_peak / core::f64::consts::TAU);
                    m += w.max(3);
                    continue;
                }
            }
            m += 1;
        }
        peaks
    };

    let mut peaks = collect_peaks(1.01);
    if peaks.is_empty() {
        peaks = collect_peaks(1.0);
    }
    (max_r_abs, peaks)
}

// A Nyquist window can make an already causal lossless response nonpassive
// between its sample bins. Try the original samples without that alteration,
// but only after the ordinary construction has failed its passivity gate.
// Both rungs must independently fit the delay, tail and unchanged power bound;
// their first three peaks must converge. No source is queried again, no gain
// is repaired, and no negative-time precursor is folded into positive time.
fn untapered_reflectance_pair(
    r_pos: &[C64],
    n_fft: usize,
    dt: f64,
    geo_delay: f64,
    ladder: &LineLadder,
) -> Option<(Vec<f64>, Vec<LineRung>)> {
    if n_fft < 512 || n_fft as f64 * 0.5 <= geo_delay + 2.0
        || r_pos.iter().any(|z| !z.re.is_finite() || !z.im.is_finite())
        || r_pos[0].im != 0.0 || r_pos.last()?.im.abs() > PASSIVITY_TOL
    {
        return None;
    }
    let mut rungs = Vec::with_capacity(2);
    let mut previous_peaks: Vec<f64> = Vec::new();
    for stride in [2, 1] {
        let n = n_fft / stride;
        let samples: Vec<_> = r_pos.iter().step_by(stride).copied().collect();
        let ir = realize_reflectance_ir_with_rolloff(
            &samples, &vec![C64::ONE; samples.len()], n, false,
        );
        if ir.iter().any(|h| !h.is_finite()) {
            return None;
        }
        let total: f64 = ir.iter().map(|h| h * h).sum();
        if !total.is_finite() {
            return None;
        }
        let tail: f64 = ir[3 * (n / 4)..].iter().map(|h| h * h).sum();
        let tail_energy = if total > 0.0 { tail / total } else { 0.0 };
        let (peak_reflection, peaks) = analyze_reflectance_ir(&ir, dt, geo_delay);
        if !peak_reflection.is_finite() || peak_reflection > 1.0 + PASSIVITY_TOL
            || !tail_energy.is_finite() || tail_energy >= ladder.tail_tol
        {
            return None;
        }
        let peak_shift_cents = if rungs.is_empty() {
            f64::INFINITY
        } else if previous_peaks.is_empty() && peaks.is_empty() {
            0.0
        } else if previous_peaks.len() != peaks.len() {
            f64::INFINITY
        } else {
            previous_peaks.iter().zip(&peaks).take(3)
                .map(|(&a, &b)| (1200.0 * (b / a).log2()).abs())
                .fold(0.0_f64, f64::max)
        };
        rungs.push(LineRung { n_fft: n, tail_energy, peak_shift_cents });
        if stride == 1 {
            return (peak_shift_cents < ladder.peak_tol_cents).then_some((ir, rungs));
        }
        previous_peaks = peaks;
    }
    None
}

/// Realize a causal characteristic line via the convergence ladder [`LineLadder`],
/// returning both the [`DelayedFilter`] and the [`LineRealization`] receipt.
///
/// Starts at `max(ceil(4 * delay_samples), 256).next_power_of_two()` and doubles
/// `n_fft` until both `tail_energy < ladder.tail_tol` and the first three resonance
/// peaks shift by `< ladder.peak_tol_cents` cents vs the preceding rung, or refuses
/// with [`DrivingPointError::LadderUnconverged`] at `ladder.max_taps`.
///
/// # Errors
/// Empty geometry, non-positive sample rate or port impedance, TMM refusal,
/// non-passive converged reflectance, or ladder non-convergence at `ladder.max_taps`.
pub fn realize_line<'a>(
    physics: &Duct,
    gas: &GasState,
    termination: impl Into<MouthLoad<'a>>,
    sample_rate_hz: u32,
    zc: f64,
    wall: Option<&WallPin>,
    ladder: &LineLadder,
) -> Result<(DelayedFilter, LineRealization), DrivingPointError> {
    let termination = match termination.into() {
        MouthLoad::Analytic(t) => t,
        MouthLoad::Tabulated(table) => return baked::characteristic(
            physics, gas, table, sample_rate_hz, usize::MAX, zc, wall, None, Some(ladder)),
        MouthLoad::Matched => return Err(baked::matched_refusal()),
    };
    let length = duct_axial_length(physics);
    if !(length > 0.0 && length.is_finite()) {
        return Err(DrivingPointError::Invalid {
            what: "line needs positive axial length",
        });
    }
    if sample_rate_hz == 0 || !zc.is_finite() || zc <= 0.0 {
        return Err(DrivingPointError::Invalid {
            what: "positive sample rate and finite positive port impedance required",
        });
    }
    let dt = 1.0 / f64::from(sample_rate_hz);
    let geo_delay = 2.0 * length / gas.sound_speed / dt;
    let mut pole_states = Vec::new();
    if let Some(gain) = sealed_volume_compliance_gain(physics, gas, termination, wall) {
        pole_states.push(PoleState::Compliance { gain });
    }
    pole_states.push(bulk_delay_pole_state(geo_delay));
    realize_characteristic_response(
        geo_delay,
        termination,
        sample_rate_hz,
        zc,
        pole_states,
        ladder,
        |omega| duct_reflectance(physics, gas, omega, LossModel::Bessel, termination, wall, zc),
    )
}

/// `R` referenced to `zc` from an input impedance, as every line realizes it.
pub(crate) fn impedance_reflectance(response: &DuctResponse, zc: f64) -> C64 {
    reflectance(C64::new(response.impedance.re, response.impedance.im), zc)
}

/// Uniform-gas duct reflectance referenced to `zc` (DEF-01). Where the stable
/// chain represents Z, R comes from Z exactly as before. At a lossless
/// resonance `|1 - r_in|` vanishes: Z is unrepresentable (fs-duct refuses
/// above its condition bound) while R stays finite, so that bin takes R from
/// the same chain's input reflection without forming Z. Every other refusal,
/// including the chain's own feedback bound, still propagates.
pub(crate) fn duct_reflectance(
    duct: &Duct,
    gas: &GasState,
    omega: f64,
    loss: LossModel,
    termination: Termination,
    wall: Option<&WallPin>,
    zc: f64,
) -> Result<C64, DuctError> {
    match input_impedance_wall(duct, gas, omega, loss, termination, wall) {
        Ok(response) => Ok(impedance_reflectance(&response, zc)),
        Err(DuctError::Singular) => {
            input_reflection_wall(duct, gas, omega, loss, termination, wall)?.referenced_to(zc)
        }
        Err(e) => Err(e),
    }
}

/// `reflectance_at(omega)` is R referenced to `zc`, never an impedance: a
/// characteristic line must stay realizable at lossless Z poles.
pub(crate) fn realize_characteristic_response<E>(
    geo_delay: f64,
    termination: impl Into<LineDc>,
    sample_rate_hz: u32,
    zc: f64,
    pole_states: Vec<PoleState>,
    ladder: &LineLadder,
    mut reflectance_at: impl FnMut(f64) -> Result<C64, E>,
) -> Result<(DelayedFilter, LineRealization), DrivingPointError>
where
    DrivingPointError: From<E>,
{
    if sample_rate_hz == 0 || !zc.is_finite() || zc <= 0.0 {
        return Err(DrivingPointError::Invalid {
            what: "positive sample rate and finite positive port impedance required",
        });
    }
    if !(geo_delay >= 2.0 && geo_delay.is_finite()) {
        return Err(DrivingPointError::Invalid {
            what: "round-trip delay does not fit the realized history",
        });
    }
    let max_taps = validate_ladder(ladder)?;
    let dt = 1.0 / f64::from(sample_rate_hz);
    let r_dc = match termination.into() {
        LineDc::Explicit(value) => {
            if !value.is_finite() || value.abs() > 1.0 + PASSIVITY_TOL {
                return Err(DrivingPointError::Invalid { what: "static reflection must be finite and passive" });
            }
            C64::new(value, 0.0)
        }
        LineDc::Analytic(Termination::Closed) => C64::new(1.0, 0.0),
        LineDc::Analytic(Termination::IdealOpen
        | Termination::UnflangedOpen
        | Termination::UnflangedThinWall { .. }
        | Termination::FlangedOpen) => {
            if let Ok(r0) = reflectance_at(1.0e-3) {
                if r0.re.is_finite() {
                    C64::new(r0.re.clamp(-1.0, 1.0), 0.0)
                } else {
                    C64::new(-1.0, 0.0)
                }
            } else {
                C64::new(-1.0, 0.0)
            }
        }
    };

    let demand_4x = (4.0 * geo_delay).ceil();
    if !demand_4x.is_finite() {
        return Err(DrivingPointError::Invalid {
            what: "round-trip delay is not finite",
        });
    }
    // Start at the 4-round-trip power-of-two rung, capped to max_taps / 2 when
    // max_taps >= 512 so at least two rungs are evaluated to measure peak_shift_cents.
    let base_n0 = (demand_4x.min(max_taps as f64) as usize)
        .max(256)
        .next_power_of_two()
        .min(max_taps);
    let n0 = if max_taps >= 512 {
        base_n0.min(max_taps / 2).max(256)
    } else {
        base_n0
    };

    let mut n_fft = n0;
    let mut r_pos: Vec<C64> = Vec::new();
    let mut prev_peaks_hz: Vec<f64> = Vec::new();
    let mut rungs_tried: Vec<LineRung> = Vec::new();

    loop {
        let half = n_fft / 2;
        let mut next_r_pos = vec![C64::ZERO; half + 1];
        next_r_pos[0] = r_dc;
        if r_pos.is_empty() {
            for (k, slot) in next_r_pos.iter_mut().enumerate().skip(1) {
                let omega = core::f64::consts::TAU * k as f64 / (n_fft as f64 * dt);
                *slot = reflectance_at(omega).map_err(DrivingPointError::from)?;
            }
        } else {
            // Reuse even bins from the previous rung (bit-identical frequencies).
            for (k_prev, &val) in r_pos.iter().enumerate() {
                next_r_pos[2 * k_prev] = val;
            }
            for k in (1..half).step_by(2) {
                let omega = core::f64::consts::TAU * k as f64 / (n_fft as f64 * dt);
                next_r_pos[k] = reflectance_at(omega).map_err(DrivingPointError::from)?;
            }
        }
        r_pos = next_r_pos;

        let ir = realize_reflectance_ir_from_positive_spectrum(&r_pos, n_fft, geo_delay, true);

        let total_energy: f64 = ir.iter().map(|&h| h * h).sum();
        let tail_start = 3 * (n_fft / 4);
        let tail_sum: f64 = ir[tail_start..].iter().map(|&h| h * h).sum();
        let tail_energy = if total_energy > 0.0 {
            tail_sum / total_energy
        } else {
            0.0
        };

        let (max_r_abs, peaks_hz) = analyze_reflectance_ir(&ir, dt, geo_delay);
        let peak_shift_cents = if rungs_tried.is_empty() {
            f64::INFINITY
        } else if prev_peaks_hz.is_empty() && peaks_hz.is_empty() {
            0.0
        } else if prev_peaks_hz.is_empty() || peaks_hz.is_empty() {
            f64::INFINITY
        } else {
            let take = prev_peaks_hz.len().min(peaks_hz.len()).min(3);
            prev_peaks_hz[..take]
                .iter()
                .zip(&peaks_hz[..take])
                .map(|(&f_prev, &f_cur)| (1200.0 * (f_cur / f_prev).log2()).abs())
                .fold(0.0_f64, f64::max)
        };

        rungs_tried.push(LineRung {
            n_fft,
            tail_energy,
            peak_shift_cents,
        });

        let fits_delay = (n_fft as f64) > geo_delay + 2.0;
        let converged = fits_delay
            && tail_energy < ladder.tail_tol
            && peak_shift_cents < ladder.peak_tol_cents;

        if converged {
            if max_r_abs > 1.0 + PASSIVITY_TOL {
                if let Some((raw, raw_rungs)) =
                    untapered_reflectance_pair(&r_pos, n_fft, dt, geo_delay, ladder)
                {
                    let last = raw_rungs.last().expect("two untapered rungs admitted");
                    let realization = LineRealization {
                        n_fft, tail_energy: last.tail_energy,
                        peak_shift_cents: last.peak_shift_cents,
                        pole_states, rungs_tried: raw_rungs,
                    };
                    let filter = DelayedFilter::from_impulse_response(dt, raw)
                        .map_err(DrivingPointError::Discrete)?;
                    return Ok((filter, realization));
                }
                return Err(DrivingPointError::Invalid {
                    what: "converged characteristic line violates 4x-oversampled passivity",
                });
            }
            let filter =
                DelayedFilter::from_impulse_response(dt, ir).map_err(DrivingPointError::Discrete)?;
            let realization = LineRealization {
                n_fft,
                tail_energy,
                peak_shift_cents,
                pole_states,
                rungs_tried,
            };
            return Ok((filter, realization));
        }

        if n_fft >= max_taps {
            return Err(DrivingPointError::LadderUnconverged {
                max_taps,
                tail_energy,
                tail_tol: ladder.tail_tol,
                peak_shift_cents,
                peak_tol_cents: ladder.peak_tol_cents,
            });
        }

        prev_peaks_hz = peaks_hz;
        n_fft *= 2;
    }
}

/// Build a causal characteristic line from a TMM duct.
///
/// Samples `R(ω) = (Z−Zc)/(Z+Zc)` on the DFT grid of the run and
/// inverse-transforms it, using [`LineLadder::default`] to select a
/// converged power-of-two FFT length where wrap-around tail energy and
/// resonance peak shifts have vanished.
///
/// # Errors
/// Empty geometry, a delay that does not fit `n` samples, TMM, or DFT.
pub fn characteristic_line<'a>(
    physics: &Duct,
    gas: &GasState,
    termination: impl Into<MouthLoad<'a>>,
    sample_rate_hz: u32,
    n: usize,
    zc: f64,
    wall: Option<&WallPin>,
) -> Result<DelayedFilter, DrivingPointError> {
    characteristic_line_dense(physics, gas, termination, sample_rate_hz, n, zc, wall, None)
}

/// [`characteristic_line`] with an optional explicit realization length.
///
/// `n_fft_override: Some(v)` pins a specific DFT length (power of two,
/// clamped to `[256, 8192]` before rounding per `u6pph`); `None` runs the
/// convergence-driven [`LineLadder::default`] up to `65_536` taps.
///
/// # Errors
/// As [`characteristic_line`].
#[allow(clippy::too_many_arguments)]
pub fn characteristic_line_dense<'a>(
    physics: &Duct,
    gas: &GasState,
    termination: impl Into<MouthLoad<'a>>,
    sample_rate_hz: u32,
    n: usize,
    zc: f64,
    wall: Option<&WallPin>,
    n_fft_override: Option<usize>,
) -> Result<DelayedFilter, DrivingPointError> {
    characteristic_line_realized(
        physics,
        gas,
        termination,
        sample_rate_hz,
        n,
        zc,
        wall,
        n_fft_override,
    )
    .map(|(line, _)| line)
}

/// [`characteristic_line_dense`] returning both the [`DelayedFilter`] and the
/// [`LineRealization`] receipt.
///
/// # Errors
/// As [`characteristic_line`].
#[allow(clippy::too_many_arguments)]
pub fn characteristic_line_realized<'a>(
    physics: &Duct,
    gas: &GasState,
    termination: impl Into<MouthLoad<'a>>,
    sample_rate_hz: u32,
    n: usize,
    zc: f64,
    wall: Option<&WallPin>,
    n_fft_override: Option<usize>,
) -> Result<(DelayedFilter, LineRealization), DrivingPointError> {
    characteristic_line_realized_with_ladder(
        physics,
        gas,
        termination,
        sample_rate_hz,
        n,
        zc,
        wall,
        n_fft_override,
        None,
    )
}

/// [`characteristic_line_realized`] with an optional caller-supplied [`LineLadder`].
///
/// # Errors
/// As [`characteristic_line`].
#[allow(clippy::too_many_arguments)]
pub fn characteristic_line_realized_with_ladder<'a>(
    physics: &Duct,
    gas: &GasState,
    termination: impl Into<MouthLoad<'a>>,
    sample_rate_hz: u32,
    n: usize,
    zc: f64,
    wall: Option<&WallPin>,
    n_fft_override: Option<usize>,
    ladder: Option<&LineLadder>,
) -> Result<(DelayedFilter, LineRealization), DrivingPointError> {
    let termination = match termination.into() {
        MouthLoad::Analytic(t) => t,
        MouthLoad::Tabulated(table) => return baked::characteristic(
            physics, gas, table, sample_rate_hz, n, zc, wall, n_fft_override, ladder),
        MouthLoad::Matched => return Err(baked::matched_refusal()),
    };
    let length = duct_axial_length(physics);
    if !(length > 0.0) {
        return Err(DrivingPointError::Invalid {
            what: "line needs positive axial length",
        });
    }
    let dt = 1.0 / f64::from(sample_rate_hz);
    let geo_delay = 2.0 * length / gas.sound_speed / dt;
    let mut pole_states = Vec::new();
    if let Some(gain) = sealed_volume_compliance_gain(physics, gas, termination, wall) {
        pole_states.push(PoleState::Compliance { gain });
    }
    pole_states.push(bulk_delay_pole_state(geo_delay));
    if n_fft_override.is_none()
        && let Some(ladder) = ladder
    {
        if !(geo_delay >= 2.0 && geo_delay < n as f64 - 2.0) {
            return Err(DrivingPointError::Invalid {
                what: "round-trip delay does not fit the realized history",
            });
        }
        return realize_characteristic_response(
            geo_delay,
            termination,
            sample_rate_hz,
            zc,
            pole_states,
            ladder,
            |omega| duct_reflectance(physics, gas, omega, LossModel::Bessel, termination, wall, zc),
        );
    }
    characteristic_line_response_with_poles(
        geo_delay,
        termination,
        sample_rate_hz,
        n,
        zc,
        n_fft_override,
        pole_states,
        |omega| duct_reflectance(physics, gas, omega, LossModel::Bessel, termination, wall, zc),
    )
}

/// [`characteristic_line_realized_with_ladder`] for a bore filled with `gas` that radiates into a different `exterior`
/// (a tract of breath into room air, say). Only a tabulated mouth load is taken, checked against the `exterior` it was
/// computed in. Every bore quantity (the delay, `zc`, the propagation and its loss) stays `gas`'s, and the load is
/// absolute, so nothing is converted at the mouth. An analytic termination would be evaluated in the bore's gas, so it
/// is not offered here. With `exterior` equal to `gas` this is the tabulated branch of
/// [`characteristic_line_realized_with_ladder`] exactly.
///
/// # Errors
/// As [`characteristic_line`], or a table whose mouth or exterior gas does not match.
#[allow(clippy::too_many_arguments)]
pub fn characteristic_line_realized_with_ladder_exterior(
    physics: &Duct,
    gas: &GasState,
    exterior: &GasState,
    table: &fs_duct::TabulatedLoad,
    sample_rate_hz: u32,
    n: usize,
    zc: f64,
    wall: Option<&WallPin>,
    n_fft_override: Option<usize>,
    ladder: Option<&LineLadder>,
) -> Result<(DelayedFilter, LineRealization), DrivingPointError> {
    baked::characteristic_in(physics, gas, exterior, table, sample_rate_hz, n, zc, wall, n_fft_override, ladder)
}

// Both uniform and regional ports use the identical sampling, Fourier sign,
// finite history and DelayedFilter admission; only the physical response differs.
fn characteristic_line_response(
    geo_delay: f64,
    termination: Termination,
    sample_rate_hz: u32,
    n: usize,
    zc: f64,
    n_fft_override: Option<usize>,
    reflectance_at: impl FnMut(f64) -> Result<C64, DuctError>,
) -> Result<DelayedFilter, DrivingPointError> {
    let pole_states = vec![bulk_delay_pole_state(geo_delay)];
    characteristic_line_response_with_poles(
        geo_delay,
        termination,
        sample_rate_hz,
        n,
        zc,
        n_fft_override,
        pole_states,
        reflectance_at,
    )
    .map(|(line, _)| line)
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn characteristic_line_response_with_poles(
    geo_delay: f64,
    termination: impl Into<LineDc>,
    sample_rate_hz: u32,
    n: usize,
    zc: f64,
    n_fft_override: Option<usize>,
    pole_states: Vec<PoleState>,
    reflectance_at: impl FnMut(f64) -> Result<C64, DuctError>,
) -> Result<(DelayedFilter, LineRealization), DrivingPointError> {
    characteristic_line_response_with_poles_ladder(geo_delay, termination, sample_rate_hz, n, zc, n_fft_override, pole_states, None, reflectance_at)
}

/// [`characteristic_line_response_with_poles`] with the convergence ladder a caller declares (`None` = [`LineLadder::default`], the bits of the function above).
/// A fixed `n_fft_override` does not use the ladder.
#[allow(clippy::too_many_arguments)]
pub(crate) fn characteristic_line_response_with_poles_ladder(
    geo_delay: f64,
    termination: impl Into<LineDc>,
    sample_rate_hz: u32,
    n: usize,
    zc: f64,
    n_fft_override: Option<usize>,
    pole_states: Vec<PoleState>,
    ladder: Option<&LineLadder>,
    mut reflectance_at: impl FnMut(f64) -> Result<C64, DuctError>,
) -> Result<(DelayedFilter, LineRealization), DrivingPointError> {
    let termination = termination.into();
    if sample_rate_hz == 0 || !zc.is_finite() || zc <= 0.0 {
        return Err(DrivingPointError::Invalid {
            what: "positive sample rate and finite positive port impedance required",
        });
    }
    let dt = 1.0 / f64::from(sample_rate_hz);
    if !(geo_delay >= 2.0 && geo_delay < n as f64 - 2.0) {
        return Err(DrivingPointError::Invalid {
            what: "round-trip delay does not fit the realized history",
        });
    }
    let Some(override_v) = n_fft_override else {
        return realize_characteristic_response(
            geo_delay,
            termination,
            sample_rate_hz,
            zc,
            pole_states,
            ladder.unwrap_or(&LineLadder::default()),
            reflectance_at,
        );
    };
    // Clamp INSIDE the documented budget BEFORE rounding: rounding
    // first overflows/saturates on oversized requests (usize::MAX-style
    // inputs panic before the cap ever applies — bead frankensim-u6pph).
    let n_fft = override_v.clamp(256, 8192).next_power_of_two();
    let fft = Fft::new(n_fft);
    let mut buf = vec![FftC64::new(0.0, 0.0); n_fft];
    let enforce_passivity = matches!(termination, LineDc::Explicit(_));
    buf[0] = match termination {
        LineDc::Explicit(value) => {
            if !value.is_finite() || value.abs() > 1.0 + PASSIVITY_TOL {
                return Err(DrivingPointError::Invalid { what: "static reflection must be finite and passive" });
            }
            FftC64::new(value, 0.0)
        }
        LineDc::Analytic(Termination::Closed) => FftC64::new(1.0, 0.0),
        LineDc::Analytic(Termination::IdealOpen
        | Termination::UnflangedOpen
        | Termination::UnflangedThinWall { .. }
        | Termination::FlangedOpen) => FftC64::new(-1.0, 0.0),
    };
    for k in 1..=n_fft / 2 {
        let omega = core::f64::consts::TAU * k as f64 / (n_fft as f64 * dt);
        let rac = reflectance_at(omega).map_err(DrivingPointError::Duct)?;
        // Acoustic e^{-iωt} → DFT e^{+iωt} so the IR is causal.
        buf[k] = FftC64::new(rac.re, -rac.im);
        if k != n_fft / 2 {
            buf[n_fft - k] = FftC64::new(rac.re, rac.im);
        }
    }
    let mut scratch = vec![FftC64::new(0.0, 0.0); n_fft];
    fft.inverse(&mut buf, &mut scratch);
    let ir: Vec<f64> = buf.iter().map(|c| c.re).collect();
    if enforce_passivity && analyze_reflectance_ir(&ir, dt, geo_delay).0 > 1.0 + PASSIVITY_TOL {
        return Err(DrivingPointError::Invalid { what: "tabulated characteristic line violates 4x-oversampled passivity" });
    }
    let total_energy: f64 = ir.iter().map(|&h| h * h).sum();
    let tail_energy = if total_energy > 0.0 {
        ir[3 * (n_fft / 4)..].iter().map(|&h| h * h).sum::<f64>() / total_energy
    } else {
        0.0
    };
    let line = DelayedFilter::from_impulse_response(dt, ir).map_err(DrivingPointError::Discrete)?;
    let realization = LineRealization {
        n_fft,
        tail_energy,
        peak_shift_cents: 0.0,
        pole_states,
        rungs_tried: vec![LineRung {
            n_fft,
            tail_energy,
            peak_shift_cents: 0.0,
        }],
    };
    Ok((line, realization))
}

/// Exact lossless normal modes below Nyquist for a uniform rigid cylinder.
/// Real-axis viscothermal wavenumbers do not determine complex poles or
/// residues: those responses stay in the original frequency-domain owner.
fn extract_lossless_cylinder_poles(
    physics: &Duct,
    gas: &GasState,
    loss: LossModel,
    termination: Termination,
    sample_rate_hz: u32,
    wall: Option<&WallPin>,
) -> Result<Vec<(f64, f64, C64)>, DrivingPointError> {
    if !matches!(loss, LossModel::Lossless)
        || wall.is_some()
        || physics.segments.len() != 1
        || !matches!(termination, Termination::Closed | Termination::IdealOpen)
    {
        return Ok(Vec::new());
    }
    let Segment::Cylinder { radius, length } = physics.segments[0] else {
        return Ok(Vec::new());
    };
    if !(radius > 0.0 && length > 0.0 && radius.is_finite() && length.is_finite()) {
        return Ok(Vec::new());
    }
    segment_wave(gas, radius, 1.0, LossModel::Lossless).map_err(DrivingPointError::Duct)?;
    let area = core::f64::consts::PI * radius * radius;
    let volume = area * length;
    let residue_re = gas.density * gas.sound_speed * gas.sound_speed / volume;
    let omega_nyq = core::f64::consts::PI * f64::from(sample_rate_hz);
    let mut poles = Vec::new();
    let mut m = 0usize;
    loop {
        let mode_half_waves = match termination {
            Termination::Closed => {
                m += 1;
                m as f64
            }
            Termination::IdealOpen => {
                let hw = m as f64 + 0.5;
                m += 1;
                hw
            }
            _ => break,
        };
        let omega_m = mode_half_waves * core::f64::consts::PI * gas.sound_speed / length;
        if omega_m >= 0.999 * omega_nyq {
            break;
        }
        poles.push((omega_m, 0.0, C64::from_re(residue_re)));
    }
    Ok(poles)
}

/// Realize a causal driving-point impedance line with explicit `PoleState::Compliance`
/// and `PoleState::Resonance` decomposition and [`LineRealization`] receipt.
///
/// Exact resonance states are available for the analytic lossless cylinder.
/// Dispersive losses remain in the original TMM response. Its residual FIR
/// doubles until measured tail energy and the first three response peaks agree
/// between adjacent rungs. Every even frequency bin is reused. High-Q decay or
/// Nyquist truncation that cannot fit the declared history refuses; no lossless
/// frequency, fitted damping, taper, or fabricated convergence is substituted.
/// `WideTube` has no admitted static limit and must use `Bessel` for this
/// full-band time-domain realization.
///
/// # Errors
/// Invalid geometry/clock/ladder, TMM refusal, unresolved FIR, or discrete filter refusal.
pub fn realize_impedance_line(
    physics: &Duct,
    gas: &GasState,
    loss: LossModel,
    termination: Termination,
    sample_rate_hz: u32,
    wall: Option<&WallPin>,
    ladder: &LineLadder,
) -> Result<(DelayedFilter, LineRealization), DrivingPointError> {
    let length = duct_axial_length(physics);
    if !(length > 0.0 && length.is_finite()) {
        return Err(DrivingPointError::Invalid {
            what: "line needs positive axial length",
        });
    }
    if sample_rate_hz == 0 {
        return Err(DrivingPointError::Invalid {
            what: "positive sample rate required",
        });
    }
    if matches!(loss, LossModel::WideTube) {
        return Err(DrivingPointError::Invalid {
            what: "full-band impedance needs a DC-admissible loss law; use Bessel instead of WideTube",
        });
    }
    let max_taps = validate_ladder(ladder)?;
    let dt = 1.0 / f64::from(sample_rate_hz);
    let geo_delay = 2.0 * length / gas.sound_speed / dt;

    let compliance_gain = impedance_compliance_gain(physics, gas, loss, termination, wall);
    let resonance_poles = extract_lossless_cylinder_poles(
        physics,
        gas,
        loss,
        termination,
        sample_rate_hz,
        wall,
    )?;

    let mut pole_states = Vec::new();
    if compliance_gain > 0.0 {
        pole_states.push(PoleState::Compliance {
            gain: compliance_gain,
        });
    }
    for &(omega, sigma, residue) in &resonance_poles {
        pole_states.push(PoleState::Resonance {
            omega,
            sigma,
            residue,
        });
    }
    pole_states.push(bulk_delay_pole_state(geo_delay));

    let n_fft = ((4.0 * geo_delay).ceil() as usize)
        .clamp(256, max_taps.min(8192))
        .next_power_of_two();

    // Purely lossless single-cylinder closed/ideal-open cavity: Mittag-Leffler
    // modal decomposition retains all poles below Nyquist with their exact
    // held-input, interval-mean response. There is no finite-FIR wrap-around
    // echo. This is a band-limited modal truncation, not a claim that omitted
    // above-Nyquist modes contribute zero to the continuous impedance.
    if matches!(loss, LossModel::Lossless) && !resonance_poles.is_empty() {
        let ir = vec![0.0_f64; n_fft];
        let line = DelayedFilter::from_impulse_response(dt, ir)
            .and_then(|f| f.with_exact_poles(compliance_gain, &resonance_poles))
            .map_err(DrivingPointError::Discrete)?;
        let realization = LineRealization {
            n_fft,
            tail_energy: 0.0,
            peak_shift_cents: 0.0,
            pole_states,
            rungs_tried: vec![LineRung {
                n_fft,
                tail_energy: 0.0,
                peak_shift_cents: 0.0,
            }],
        };
        return Ok((line, realization));
    }

    let response_at = |omega| {
        input_impedance_wall(physics, gas, omega, loss, termination, wall)
            .map(|response| impedance_observable(response, termination))
    };
    // Preserve the finite real part of the original quasi-static response;
    // the singular imaginary compliance is carried by its explicit state.
    let dc = response_at(1.0e-3).map_err(DrivingPointError::Duct)?.re;
    realize_impedance_response(
        geo_delay, dt, dc, compliance_gain, pole_states, ladder, response_at,
    )
}

fn impedance_compliance_gain(
    physics: &Duct, gas: &GasState, loss: LossModel,
    termination: Termination, wall: Option<&WallPin>,
) -> f64 {
    let adiabatic = sealed_volume_compliance_gain(physics, gas, termination, wall)
        .unwrap_or(0.0);
    // The original Bessel thermal function tends to one at DC: the gas bulk
    // modulus is rho*c^2/gamma, not its high-frequency adiabatic value.
    if matches!(loss, LossModel::Bessel) { adiabatic / gas.gamma } else { adiabatic }
}

fn impedance_observable(response: DuctResponse, termination: Termination) -> C64 {
    match termination {
        Termination::Closed | Termination::IdealOpen => response.impedance,
        Termination::UnflangedOpen
        | Termination::UnflangedThinWall { .. }
        | Termination::FlangedOpen => response.p_mouth_over_u_in,
    }
}

/// The complete playback response, including its retained static storage.
fn impedance_ir_response(ir: &[f64], dt: f64, compliance_gain: f64, omega: f64) -> C64 {
    let theta = omega * dt;
    let (sin_step, cos_step) = fs_math::det::sin_cos(theta);
    let mut phase = C64::ONE;
    let step = C64::new(cos_step, sin_step);
    let mut response = C64::ZERO;
    for (index, &tap) in ir.iter().enumerate() {
        if index != 0 && index.is_multiple_of(64) {
            let (sine, cosine) = fs_math::det::sin_cos(theta * index as f64);
            phase = C64::new(cosine, sine);
        }
        response = response + phase.scale(tap);
        phase = phase * step;
    }
    if compliance_gain > 0.0 {
        let (sine, cosine) = fs_math::det::sin_cos(0.5 * theta);
        response.im += 0.5 * compliance_gain * dt * cosine / sine;
    }
    response
}

// Fourfold zero padding locates peaks of the ACTUAL FIR + compliance response.
// Refinement reads its continuous DTFT, not the samples that authored the FIR.
// This is a finite sampled peak check, not a proof of spectrum completeness.
fn impedance_ir_peaks(
    ir: &[f64], dt: f64, compliance_gain: f64,
) -> Result<Vec<f64>, DrivingPointError> {
    let count = 4 * ir.len();
    let fft = Fft::new(count);
    let mut data = vec![FftC64::new(0.0, 0.0); count];
    let scale = ir.iter().map(|value| value.abs()).fold(0.0_f64, f64::max)
        .max(compliance_gain * dt);
    if scale == 0.0 {
        return Ok(Vec::new());
    }
    if !scale.is_finite() {
        return Err(DrivingPointError::Invalid { what: "unrepresentable impedance response scale" });
    }
    for (entry, &value) in data.iter_mut().zip(ir) {
        entry.re = value / scale;
    }
    let mut scratch = vec![FftC64::new(0.0, 0.0); count];
    fft.forward(&mut data, &mut scratch);
    let spacing = core::f64::consts::TAU / (count as f64 * dt);
    let mut amplitude = vec![0.0; count / 2 + 1];
    for index in 1..=count / 2 {
        let mut value = C64::new(data[index].re, -data[index].im);
        if compliance_gain > 0.0 {
            let (sine, cosine) = fs_math::det::sin_cos(core::f64::consts::PI * index as f64 / count as f64);
            value.im += 0.5 * (compliance_gain * dt / scale) * cosine / sine;
        }
        amplitude[index] = value.abs();
        if !amplitude[index].is_finite() {
            return Err(DrivingPointError::Invalid { what: "nonfinite impedance peak response" });
        }
    }
    let mut peaks = Vec::with_capacity(3);
    for index in 2..count / 2 {
        if amplitude[index] <= amplitude[index - 1]
            || amplitude[index] <= amplitude[index + 1]
        {
            continue;
        }
        let mut left = (index - 1) as f64 * spacing;
        let mut right = (index + 1) as f64 * spacing;
        let ratio = 0.5 * (5.0_f64.sqrt() - 1.0);
        let mut a = right - ratio * (right - left);
        let mut b = left + ratio * (right - left);
        let mut fa = impedance_ir_response(ir, dt, compliance_gain, a).abs();
        let mut fb = impedance_ir_response(ir, dt, compliance_gain, b).abs();
        for _ in 0..32 {
            if fa >= fb {
                right = b;
                b = a;
                fb = fa;
                a = right - ratio * (right - left);
                fa = impedance_ir_response(ir, dt, compliance_gain, a).abs();
            } else {
                left = a;
                a = b;
                fa = fb;
                b = left + ratio * (right - left);
                fb = impedance_ir_response(ir, dt, compliance_gain, b).abs();
            }
        }
        if !fa.is_finite() || !fb.is_finite() {
            return Err(DrivingPointError::Invalid { what: "nonfinite refined impedance response" });
        }
        peaks.push(f64::midpoint(left, right) / core::f64::consts::TAU);
        if peaks.len() == 3 { break; }
    }
    Ok(peaks)
}

#[allow(clippy::too_many_arguments)]
fn realize_impedance_response(
    geo_delay: f64, dt: f64, dc_impedance: f64, compliance_gain: f64,
    pole_states: Vec<PoleState>, ladder: &LineLadder,
    mut response_at: impl FnMut(f64) -> Result<C64, DuctError>,
) -> Result<(DelayedFilter, LineRealization), DrivingPointError> {
    let max_taps = validate_ladder(ladder)?;
    if !geo_delay.is_finite() || geo_delay <= 0.0 || !dc_impedance.is_finite() {
        return Err(DrivingPointError::Invalid { what: "finite positive delay and finite static response required" });
    }
    let initial = ((4.0 * geo_delay).ceil().min(max_taps as f64) as usize)
        .max(256).next_power_of_two().min(max_taps);
    let mut count = initial.min((max_taps / 2).max(256));
    let mut samples: Vec<C64> = Vec::new();
    let mut previous_peaks: Vec<f64> = Vec::new();
    let mut rungs = Vec::new();
    loop {
        let mut next = vec![C64::ZERO; count / 2 + 1];
        next[0] = C64::from_re(dc_impedance);
        if samples.is_empty() {
            for (index, value) in next.iter_mut().enumerate().skip(1) {
                *value = response_at(core::f64::consts::TAU * index as f64 / (count as f64 * dt))?;
            }
        } else {
            for (index, &value) in samples.iter().enumerate() { next[2 * index] = value; }
            for index in (1..count / 2).step_by(2) {
                next[index] = response_at(core::f64::consts::TAU * index as f64 / (count as f64 * dt))?;
            }
        }
        samples = next;
        let ir = impedance_ir_from_spectrum(dt, &samples, compliance_gain, &[])?;
        let scale = ir.iter().map(|value| value.abs()).fold(0.0_f64, f64::max);
        let tail_energy = if scale > 0.0 {
            let total: f64 = ir.iter().map(|value| (value / scale).powi(2)).sum();
            let tail: f64 = ir[3 * (count / 4)..].iter().map(|value| (value / scale).powi(2)).sum();
            tail / total
        } else { 0.0 };
        let peaks = impedance_ir_peaks(&ir, dt, compliance_gain)?;
        let peak_shift_cents = if rungs.is_empty() || peaks.len() != previous_peaks.len() {
            f64::INFINITY
        } else {
            previous_peaks.iter().zip(&peaks)
                .map(|(&a, &b)| (1200.0 * (b / a).log2()).abs())
                .fold(0.0_f64, f64::max)
        };
        rungs.push(LineRung { n_fft: count, tail_energy, peak_shift_cents });
        if count as f64 > geo_delay + 2.0
            && tail_energy < ladder.tail_tol && peak_shift_cents < ladder.peak_tol_cents
        {
            let line = DelayedFilter::from_impulse_response(dt, ir)
                .and_then(|line| line.with_exact_poles(compliance_gain, &[]))
                .map_err(DrivingPointError::Discrete)?;
            return Ok((line, LineRealization {
                n_fft: count, tail_energy, peak_shift_cents, pole_states, rungs_tried: rungs,
            }));
        }
        if count == max_taps {
            return Err(DrivingPointError::LadderUnconverged {
                max_taps, tail_energy, tail_tol: ladder.tail_tol,
                peak_shift_cents, peak_tol_cents: ladder.peak_tol_cents,
            });
        }
        previous_peaks = peaks;
        count *= 2;
    }
}

#[allow(clippy::too_many_arguments)]
fn build_impedance_filter(
    physics: &Duct,
    gas: &GasState,
    loss: LossModel,
    termination: Termination,
    dt: f64,
    n_fft: usize,
    wall: Option<&WallPin>,
    compliance_gain: f64,
    resonance_poles: &[(f64, f64, C64)],
) -> Result<DelayedFilter, DrivingPointError> {
    let response_at = |omega| {
        input_impedance_wall(physics, gas, omega, loss, termination, wall)
            .map(|response| impedance_observable(response, termination))
    };
    let dc = response_at(1.0e-3).map_err(DrivingPointError::Duct)?.re;
    build_impedance_response_filter(dt, n_fft, dc, compliance_gain, resonance_poles, response_at)
}

fn build_impedance_response_filter(
    dt: f64, n_fft: usize, dc_impedance: f64, compliance_gain: f64,
    resonance_poles: &[(f64, f64, C64)],
    mut response_at: impl FnMut(f64) -> Result<C64, DuctError>,
) -> Result<DelayedFilter, DrivingPointError> {
    let mut samples = vec![C64::ZERO; n_fft / 2 + 1];
    samples[0] = C64::from_re(dc_impedance);
    for (index, sample) in samples.iter_mut().enumerate().skip(1) {
        let omega = core::f64::consts::TAU * index as f64 / (n_fft as f64 * dt);
        *sample = response_at(omega).map_err(DrivingPointError::Duct)?;
    }
    let ir = impedance_ir_from_spectrum(dt, &samples, compliance_gain, resonance_poles)?;
    DelayedFilter::from_impulse_response(dt, ir)
        .and_then(|f| f.with_exact_poles(compliance_gain, resonance_poles))
        .map_err(DrivingPointError::Discrete)
}

fn impedance_ir_from_spectrum(
    dt: f64, samples: &[C64], compliance_gain: f64,
    resonance_poles: &[(f64, f64, C64)],
) -> Result<Vec<f64>, DrivingPointError> {
    if samples.iter().any(|value| !value.re.is_finite() || !value.im.is_finite()) {
        return Err(DrivingPointError::Invalid { what: "nonfinite original impedance sample" });
    }
    let n_fft = 2 * (samples.len() - 1);
    // Prepare the SAME coefficients that the committing playback owner uses.
    // The four zero taps merely satisfy its FIR-port shape; no output is sampled
    // and no history is advanced while evaluating the response below.
    let pole_response = DelayedFilter::from_impulse_response(dt, vec![0.0; 4])
        .and_then(|f| f.with_exact_poles(0.0, resonance_poles))
        .map_err(DrivingPointError::Discrete)?;
    let fft = Fft::new(n_fft);
    let mut buf = vec![FftC64::new(0.0, 0.0); n_fft];
    let pole_dc = pole_response.resonance_response_e_minus_iwt(0.0)
        .map_err(DrivingPointError::Discrete)?;
    buf[0] = FftC64::new(samples[0].re - pole_dc.re, 0.0);
    for k in 1..=n_fft / 2 {
        let omega = core::f64::consts::TAU * k as f64 / (n_fft as f64 * dt);
        let mut z = samples[k];
        if compliance_gain > 0.0 && k < n_fft / 2 {
            // Peel the discrete trapezoidal compliance integrator:
            // Z_comp(ω_k) = +i (g_comp * dt / 2) cot(ω_k dt / 2) under e^{-iωt}.
            let half_theta = 0.5 * omega * dt;
            let cot = half_theta.cos() / half_theta.sin();
            z.im -= 0.5 * compliance_gain * dt * cot;
        }
        z = z - pole_response.resonance_response_e_minus_iwt(omega)
            .map_err(DrivingPointError::Discrete)?;
        // Acoustic e^{-iωt} → DFT e^{+iωt} so the IR is causal,
        // matching characteristic_line_response (DEF-01).
        buf[k] = FftC64::new(z.re, -z.im);
        if k != n_fft / 2 {
            buf[n_fft - k] = FftC64::new(z.re, z.im);
        }
    }
    let mut scratch = vec![FftC64::new(0.0, 0.0); n_fft];
    fft.inverse(&mut buf, &mut scratch);
    let ir: Vec<f64> = buf.iter().map(|c| c.re).collect();
    if ir.iter().any(|value| !value.is_finite()) {
        return Err(DrivingPointError::Invalid { what: "nonfinite impedance impulse response" });
    }
    Ok(ir)
}

/// Causal driving-point impedance as the same FIR + pole-state port object.
///
/// Samples `Z(ω)` (closed) or the mouth transfer `p/u` (open) on the
/// DFT grid of the run, peeling the Bessel isothermal `PoleState::Compliance`
/// (`g_comp = ρ c² / (γ V)`)
/// before the IFFT for sealed cavities so a closed volume retains its static
/// compliance without `1/(iω)` wrap-around aliasing. Like [`characteristic_line`],
/// TMM `Z` (`e^{-iωt}`) is conjugated onto the DFT grid (`e^{+iωt}`) so the
/// realized impulse response is causal (DEF-01).
///
/// # Errors
/// Empty geometry, TMM, or DFT.
pub fn impedance_line<'a>(
    physics: &Duct,
    gas: &GasState,
    termination: impl Into<MouthLoad<'a>>,
    sample_rate_hz: u32,
    n: usize,
    wall: Option<&WallPin>,
) -> Result<DelayedFilter, DrivingPointError> {
    let termination = match termination.into() {
        MouthLoad::Analytic(t) => t,
        MouthLoad::Tabulated(table) => return baked::impedance(physics, gas, table, sample_rate_hz, n, wall),
        MouthLoad::Matched => return Err(baked::matched_refusal()),
    };
    if physics.segments.is_empty() {
        return Err(DrivingPointError::Invalid {
            what: "line needs positive axial length",
        });
    }
    if sample_rate_hz == 0 {
        return Err(DrivingPointError::Invalid {
            what: "positive sample rate required",
        });
    }
    let dt = 1.0 / f64::from(sample_rate_hz);
    let n_fft = n.clamp(256, 8192).next_power_of_two();
    let compliance_gain = impedance_compliance_gain(
        physics, gas, LossModel::Bessel, termination, wall,
    );
    build_impedance_filter(
        physics,
        gas,
        LossModel::Bessel,
        termination,
        dt,
        n_fft,
        wall,
        compliance_gain,
        &[],
    )
}

#[cfg(test)]
mod convention_tests;
#[cfg(test)]
mod ladder_tests;

