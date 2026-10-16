//! Fluctuation strength [vacil] — implemented from the PUBLISHED
//! model equations of Osses Vecchi, García León & Kohlrausch,
//! "Modelling the sensation of fluctuation strength", ICA 2016 paper
//! ICA2016-113 (openly hosted; PDF archived in the session record).
//! NO reference code was ported: the only open implementations are
//! CC BY-NC (SQAT) or unlicensed (fluctuation-strength-TUe), so this
//! module implements the paper's equations directly — methods are
//! not copyrightable, code is — and validates against the PUBLISHED
//! Fastl-Zwicker values reproduced in the paper's Table 1.
//!
//! Model (the paper's Eq. 1, adapted from the Daniel-Weber/Aures
//! roughness chain this crate already ports in [`crate::roughness`]):
//!
//! `FS = C_FS * sum_{i=1..47} (m*_i)^1.7 * |k_{i-2} k_i|^1.7 * g(z_i)`
//!
//! with the excitation-pattern stages of Eqs. 2-9: Terhardt slopes
//! (S1 = 27 dB/Bark, S2 = 24 + 230/f - 0.2 L), the Zwicker-Terhardt
//! Bark formula, 47 channels at 0.5 i Bark, generalised modulation
//! depth `m* = rms(h_BP)/h_0` with a 3:1 COMPRESSION above 0.7
//! (instead of roughness's clamp at 1), the modulation band-pass
//! `H(fmod)` band-pass (the paper: 4th-order LPF at 12 Hz x
//! 2nd-order HPF at 3.1 Hz; THIS realisation: see below), the
//! normalised cross covariance between channels i and i±2, and
//! `g(z)` falling linearly 1 -> 0.5 over 15..23.5 Bark.
//!
//! Stated deviations / choices the paper leaves open (all disclosed;
//! the gates in the conformance tests measure their net effect
//! against the published values):
//! - Frame: 131072 samples (2.73 s at 48 kHz, power of two for
//!   fs-fft; the paper uses 2-s frames at exactly 0.5 Hz resolution).
//!   Single-frame API v1 (no 90%-overlap frame train); the paper's
//!   50-ms raised-cosine gates are applied to the frame.
//! - `H(fmod)` filter TYPE is unstated; Butterworth-style magnitude
//!   responses are used, applied ZERO-PHASE in the frequency domain
//!   (identical filtering for every channel, so the cross covariance
//!   is unaffected by the phase choice), and the corner parameters
//!   are FITTED to the published AM-tone curve per the paper's own
//!   methodology (see `F_HP`/`F_LP` docs).
//! - `C_FS` is re-derived on the paper's own reference stimulus
//!   (1 kHz, 60 dB, m = 1, fmod = 4 Hz := 1.00 vacil) because the
//!   published 0.2490 is specific to their unstated filter; the
//!   paper itself fixes the constant this way.
//! - The absolute-threshold gate reuses the Daniel-Weber
//!   threshold-in-quiet table ([`crate::dw_tables`]) — the chain the
//!   model is adapted from.
//! - Edge channels (i = 1, 2 and 46, 47) use the single available
//!   cross-covariance factor (the paper's Eq. 1 leaves the edges
//!   undefined; this mirrors the roughness reference's convention).
//!
//! THE VALUES ARE NEVER A SUBSTITUTE FOR HUMAN LISTENING (crate law).

use crate::dw_tables::{A0_BARK, A0_DB, LTQ_R_BARK, LTQ_R_DB};
use crate::{PsychoError, tables_interp};
use fs_fft::{C64, Fft};
use fs_math::det;

/// Analysis frame length (samples): 2.73 s at 48 kHz, `Delta f`
/// = 0.366 Hz (the paper: 2 s / 0.5 Hz; power-of-two deviation
/// stated in the module doc).
pub const FS_FRAME: usize = 131_072;

/// Number of critical-band channels.
const N_CHANNEL: usize = 47;
/// Compression knee and ratio for the generalised modulation depth.
const M_KNEE: f64 = 0.7;
const M_RATIO: f64 = 3.0;
/// Model exponents (paper section 3.1).
const P_M: f64 = 1.7;
const P_K: f64 = 1.7;
/// H(fmod) band-pass parameters — FITTED to the published AM-tone
/// curve per the paper's own methodology ("the bandpass filter H was
/// fitted by using 1-kHz AM tones"): the paper's stated 3.1-12 Hz
/// corners belong to their unstated IIR realisation; THIS
/// realisation (pure Butterworth-style magnitudes, zero phase) fits
/// the same published data with a 1st-order high-pass at 1.7 Hz and
/// a 3rd-order low-pass at 13 Hz (executed: the literal 3.1 Hz
/// 2nd-order high-pass read FS(1 Hz) = 0.05 vs the published 0.39).
const F_HP: f64 = 1.7;
const F_LP: f64 = 13.0;
/// Calibration constant: derived on the paper's reference stimulus
/// (1 kHz, 60 dB SPL, m = 1, fmod = 4 Hz := 1.00 vacil) by running
/// THIS implementation once; the non-ignored battery asserts the
/// reference reads 1.0 against this pinned constant (a REGRESSION
/// pin — re-derivation itself lives in the ignored calibration
/// probe). The paper's own 0.2490 belongs to their unstated filter
/// realisation.
pub const C_FS: f64 = 0.12753996479900015;

/// Zwicker-Terhardt critical-band rate [Bark] (paper Eq. 3; the
/// canonical 7.6e-4 coefficient — the ICA PDF's OCR renders it
/// ambiguously, settled by the z(1000 Hz) ~ 8.5 Bark anchor and
/// cross-checked against the independent Daniel-Weber BARK_FREQS
/// table by a conformance test).
#[must_use]
pub fn bark_z(f_hz: f64) -> f64 {
    13.0 * det::atan(7.6e-4 * f_hz) + 3.5 * det::atan((f_hz / 7500.0) * (f_hz / 7500.0))
}

/// |H(fmod)|: fitted band-pass magnitude — 3rd-order low-pass at
/// 13 Hz times 1st-order high-pass at 1.7 Hz (BOTH the corners and
/// the ORDERS differ from the paper's stated 4th/2nd @ 12/3.1 Hz
/// realisation; fitted to the published AM curve per the paper's
/// own methodology — see `F_HP`/`F_LP`). Zero-phase application.
fn h_mod(f: f64) -> f64 {
    if f <= 0.0 {
        return 0.0;
    }
    // 3rd-order low-pass magnitude, 1st-order high-pass magnitude.
    let lp = 1.0 / det::sqrt(1.0 + det::pow(f / F_LP, 6.0));
    let r = f / F_HP;
    let hp = r / det::sqrt(1.0 + r * r);
    lp * hp
}

/// g(z): unity below 15 Bark, falling linearly to 0.5 at 23.5 Bark
/// (paper section 3.1).
fn g_weight(z: f64) -> f64 {
    if z < 15.0 {
        1.0
    } else {
        1.0 - 0.5 * (z - 15.0) / (23.5 - 15.0)
    }
}

/// The paper's 3:1 modulation-depth compression above the 0.7 knee
/// (section 2.1.3; worked example pinned by test: 0.85 -> 0.75).
/// Public because the tone battery never drives m* past the knee
/// (review-executed: a clamp mutant passed every tone test) — the
/// stage is pinned directly.
#[must_use]
pub fn compress_modulation_depth(m: f64) -> f64 {
    if m > M_KNEE {
        M_KNEE + (m - M_KNEE) / M_RATIO
    } else {
        m
    }
}

/// Pearson normalised cross covariance (paper Eq. 9).
fn cross_covariance(a: &[f64], b: &[f64]) -> f64 {
    #[allow(clippy::cast_precision_loss)]
    let n = a.len() as f64;
    let ma = a.iter().sum::<f64>() / n;
    let mb = b.iter().sum::<f64>() / n;
    let mut num = 0.0;
    let mut da = 0.0;
    let mut db = 0.0;
    for (&x, &y) in a.iter().zip(b) {
        num += (x - ma) * (y - mb);
        da += (x - ma) * (x - ma);
        db += (y - mb) * (y - mb);
    }
    let denom = det::sqrt(da * db);
    if denom > 0.0 { num / denom } else { 0.0 }
}

/// Fluctuation strength of one analysis frame (sound pressure in
/// Pa), [vacil]. Uses the first [`FS_FRAME`] samples.
///
/// # Errors
/// [`PsychoError::DegenerateSignal`] on a short frame;
/// [`PsychoError::NonFinite`] on bad samples or rate.
#[allow(clippy::too_many_lines)] // one published pipeline, told in the paper's stage order
pub fn fluctuation_strength_frame(pcm_pa: &[f64], sample_rate: f64) -> Result<f64, PsychoError> {
    if pcm_pa.len() < FS_FRAME {
        return Err(PsychoError::DegenerateSignal {
            what: "fluctuation strength needs a full analysis frame",
        });
    }
    if !sample_rate.is_finite() || sample_rate <= 0.0 {
        return Err(PsychoError::NonFinite {
            what: "sample rate",
        });
    }
    if pcm_pa.iter().take(FS_FRAME).any(|v| !v.is_finite()) {
        return Err(PsychoError::NonFinite { what: "pcm sample" });
    }
    let n = FS_FRAME;
    #[allow(clippy::cast_precision_loss)]
    let nf = n as f64;
    let fs = sample_rate;
    // --- 50-ms raised-cosine on/off gates (paper section 2.1). ---
    #[allow(clippy::cast_sign_loss, clippy::cast_possible_truncation)]
    let ramp = (0.05 * fs) as usize;
    let fft = Fft::new(n);
    let mut scratch = vec![C64::new(0.0, 0.0); n];
    let mut buf: Vec<C64> = (0..n)
        .map(|i| {
            #[allow(clippy::cast_precision_loss)]
            let g = if i < ramp {
                0.5 - 0.5 * det::cos(core::f64::consts::PI * i as f64 / ramp as f64)
            } else if i >= n - ramp {
                0.5 - 0.5 * det::cos(core::f64::consts::PI * (n - 1 - i) as f64 / ramp as f64)
            } else {
                1.0
            };
            C64::new(pcm_pa[i] * g, 0.0)
        })
        .collect();
    fft.forward(&mut buf, &mut scratch);
    let half = n / 2;
    // --- Component levels with a0 transmission, thresholded. ---
    // Amplitude of component k: 2/N * |X_k| (one-sided real signal).
    let mut level = vec![f64::NEG_INFINITY; half];
    let mut zk = vec![0.0f64; half];
    let mut audible: Vec<usize> = Vec::new();
    for k in 1..half {
        #[allow(clippy::cast_precision_loss)]
        let f = k as f64 * fs / nf;
        let z = bark_z(f);
        zk[k] = z;
        let amp = 2.0 * det::sqrt(buf[k].norm_sq()) / nf;
        if amp <= 0.0 {
            continue;
        }
        let a0_db = tables_interp(z, &A0_BARK, &A0_DB);
        let l = 20.0 * det::ln(amp / 2.0e-5) / det::ln(10.0) + a0_db;
        level[k] = l;
        let thr = tables_interp(z, &LTQ_R_BARK, &LTQ_R_DB);
        if l > thr {
            audible.push(k);
        }
    }
    if audible.is_empty() {
        return Ok(0.0);
    }
    // --- Per-channel excitation, envelope, modulation depth. ---
    let mut h_bp: Vec<Vec<f64>> = Vec::with_capacity(N_CHANNEL);
    let mut m_star = [0.0f64; N_CHANNEL];
    let mut exc = vec![C64::new(0.0, 0.0); n];
    let mut env = vec![C64::new(0.0, 0.0); n];
    #[allow(clippy::needless_range_loop)] // ch indexes m_star AND h_bp push order
    for ch in 0..N_CHANNEL {
        #[allow(clippy::cast_precision_loss)]
        let zi = f64::midpoint(ch as f64, 1.0);
        // Excitation spectrum at observation point zi (Eqs. 4-5).
        for e in &mut exc {
            *e = C64::new(0.0, 0.0);
        }
        for &k in &audible {
            let lk = level[k];
            let dz = zi - zk[k];
            let l_contrib = if dz > 0.0 {
                // Component below the observation point: upper slope.
                #[allow(clippy::cast_precision_loss)]
                let f = k as f64 * fs / nf;
                let s2 = 24.0 + 230.0 / f - 0.2 * lk;
                lk - s2 * dz
            } else {
                // Component above: lower slope, 27 dB/Bark.
                lk + 27.0 * dz
            };
            if l_contrib <= -80.0 {
                continue; // negligible tail (numerical economy)
            }
            let amp = 2.0e-5 * det::pow(10.0, l_contrib / 20.0);
            // Keep the component's phase (paper: "keeping the same
            // phase of the component at k").
            let mag = det::sqrt(buf[k].norm_sq());
            if mag > 0.0 {
                let scale = amp * nf / 2.0 / mag;
                exc[k] = C64::new(buf[k].re * scale, buf[k].im * scale);
            }
        }
        // Hermitian mirror so the IFFT is real.
        for k in 1..half {
            exc[n - k] = C64::new(exc[k].re, -exc[k].im);
        }
        let mut tmp = exc.clone();
        fft.inverse(&mut tmp, &mut scratch);
        // Full-wave rectified envelope and its DC value (Eq. 6).
        for (e, t) in env.iter_mut().zip(&tmp) {
            *e = C64::new(t.re.abs(), 0.0);
        }
        let h0 = env.iter().map(|c| c.re).sum::<f64>() / nf;
        // Weighted envelope h_BP (Eq. 7): H(fmod) in the frequency
        // domain, zero-phase.
        fft.forward(&mut env, &mut scratch);
        for (k, e) in env.iter_mut().enumerate() {
            let kk = if k <= n / 2 { k } else { n - k };
            #[allow(clippy::cast_precision_loss)]
            let fmod = kk as f64 * fs / nf;
            let w = h_mod(fmod);
            *e = C64::new(e.re * w, e.im * w);
        }
        fft.inverse(&mut env, &mut scratch);
        let series: Vec<f64> = env.iter().map(|c| c.re).collect();
        let rms = det::sqrt(series.iter().map(|v| v * v).sum::<f64>() / nf);
        // Generalised modulation depth with 3:1 compression above
        // the 0.7 knee (paper section 2.1.3).
        let m = if h0 > 0.0 { rms / h0 } else { 0.0 };
        m_star[ch] = compress_modulation_depth(m);
        h_bp.push(series);
    }
    // --- Cross covariances between channels i and i+2 (Eq. 9). ---
    let mut k_fwd = [0.0f64; N_CHANNEL];
    for i in 0..N_CHANNEL - 2 {
        k_fwd[i] = cross_covariance(&h_bp[i], &h_bp[i + 2]);
    }
    // --- Total FS (Eq. 1) with edge-channel single-factor rule. ---
    let mut total = 0.0f64;
    for ch in 0..N_CHANNEL {
        #[allow(clippy::cast_precision_loss)]
        let zi = f64::midpoint(ch as f64, 1.0);
        let k_factor = if ch < 2 {
            k_fwd[ch].abs()
        } else if ch >= N_CHANNEL - 2 {
            k_fwd[ch - 2].abs()
        } else {
            (k_fwd[ch - 2] * k_fwd[ch]).abs()
        };
        total += det::pow(m_star[ch], P_M) * det::pow(k_factor, P_K) * g_weight(zi);
    }
    Ok(C_FS * total)
}
