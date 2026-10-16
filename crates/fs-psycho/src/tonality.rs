//! ECMA-74 annex D tonality: tone-to-noise ratio (TNR) and prominence
//! ratio (PR) with prominence classification, ported
//! statement-for-statement from the Apache-2.0 MoSQITo reference
//! implementation (github.com/Eomys/MoSQITo, commit d990c33:
//! `sq_metrics/tonality/{tone_to_noise_ecma,prominence_ratio_ecma}`,
//! Bray screening with 1/24-octave smoothed-spectrum criteria; totals
//! per ECMA TR/108 over the prominent tones).
//!
//! Scope facts, stated:
//! - Single stationary block, Pa input, power-of-two length (fs-fft;
//!   the reference FFTs the whole signal at its native length — the
//!   reference was re-run at the fixture lengths to mint every pin).
//! - The spectrum front end is the reference's `comp_spectrum`:
//!   sum-normalized Hann window, one-sided spectrum scaled by the
//!   empirical 1.42, frequency axis labeled `(k+1) fs/n` (the
//!   reference's own one-bin offset, kept), dB re 2e-5 Pa.
//! - Detection band is the reference's (89.1, 11200) Hz.
//! - Disclosed deviations from the reference, all on paths its own
//!   fixtures never exercise: (a) the reference paints the smoothed
//!   spectrum into a numpy `empty` array and can leave the LAST bin
//!   uninitialized memory — here the tail is painted with the last
//!   band's value (deterministic completion, head and tail); (b) a
//!   zero spectrum magnitude reads -140 dB, matching the reference's
//!   amp2db 2e-12 substitution; (c) hearing-threshold coefficients
//!   (unreachable through the detection band) clamp to the nearest
//!   band where the reference raises NameError; (d) an empty
//!   1/24-octave smoothing band REFUSES instead of the reference's
//!   corrupted bookkeeping; (e) a screening walk past the band start
//!   clamps where numpy index-wraps.
//!
//! THE VALUES ARE NEVER A SUBSTITUTE FOR HUMAN LISTENING (crate law).

use crate::PsychoError;
use fs_fft::{C64, Fft};
use fs_math::det;

/// One detected tonal component.
#[derive(Debug, Clone)]
pub struct Tone {
    /// Tone frequency [Hz] (the spectral-line label of the peak).
    pub frequency_hz: f64,
    /// TNR or PR value [dB].
    pub ratio_db: f64,
    /// ECMA-74 prominence criterion.
    pub prominent: bool,
}

/// A tonality result: the detected tones plus the total over the
/// PROMINENT tones (ECMA TR/108 energy sum; 0.0 when no tone is
/// prominent — the reference's own convention, kept).
#[derive(Debug, Clone)]
pub struct Tonality {
    /// Detected tones (all with positive ratio, prominent or not).
    pub tones: Vec<Tone>,
    /// Total level over prominent tones [dB]; 0.0 if none.
    pub total_db: f64,
}

/// Detection band bounds [Hz] (reference constants).
const F_LOW: f64 = 89.1;
const F_HIGH: f64 = 11_200.0;

#[allow(clippy::cast_precision_loss)]
fn usize_f64(x: usize) -> f64 {
    x as f64
}

/// numpy-style sign.
fn sign(x: f64) -> f64 {
    if x > 0.0 {
        1.0
    } else if x < 0.0 {
        -1.0
    } else {
        0.0
    }
}

fn db_to_pow(db: f64) -> f64 {
    det::pow(10.0, db / 10.0)
}

fn pow_to_db(p: f64) -> f64 {
    10.0 * det::ln(p) / det::ln(10.0)
}

/// First index minimizing |xs[i] - target| (numpy argmin: first
/// occurrence wins).
fn argmin_abs(xs: &[f64], target: f64) -> usize {
    let mut best = 0usize;
    let mut best_v = f64::INFINITY;
    for (i, &x) in xs.iter().enumerate() {
        let d = (x - target).abs();
        if d < best_v {
            best_v = d;
            best = i;
        }
    }
    best
}

/// The reference's `comp_spectrum` (db=True): sum-normalized Hann,
/// one-sided x1.42, axis (k+1)fs/n, dB re 2e-5.
fn spectrum_db(pcm_pa: &[f64], sample_rate: f64) -> (Vec<f64>, Vec<f64>) {
    let n = pcm_pa.len();
    let mut win = vec![0.0f64; n];
    let mut wsum = 0.0;
    for (i, w) in win.iter_mut().enumerate() {
        *w = 0.5 - 0.5 * det::cos(2.0 * core::f64::consts::PI * usize_f64(i) / usize_f64(n - 1));
        wsum += *w;
    }
    let fft = Fft::new(n);
    let mut scratch = vec![C64::new(0.0, 0.0); n];
    let mut buf: Vec<C64> = pcm_pa
        .iter()
        .zip(&win)
        .map(|(x, w)| C64::new(x * w / wsum, 0.0))
        .collect();
    fft.forward(&mut buf, &mut scratch);
    let half = n / 2;
    let db: Vec<f64> = buf[..half]
        .iter()
        .map(|c| {
            let m = det::sqrt(c.norm_sq()) * 1.42;
            if m > 0.0 {
                20.0 * det::ln(m / 2.0e-5) / det::ln(10.0)
            } else {
                // The reference's amp2db substitutes 2e-12 for a
                // zero amplitude: exactly -140 dB re 2e-5.
                -140.0
            }
        })
        .collect();
    let freq: Vec<f64> = (0..half)
        .map(|k| usize_f64(k + 1) * sample_rate / usize_f64(n))
        .collect();
    (db, freq)
}

/// The reference's `_getFrequencies(90, 11200, 24, G=10, fr=1000)`
/// 1/24-octave band edges, including its loop shape (the final band
/// is appended before the `f2 <= fend` check fails) and the
/// first/last edge clamps applied by `_spectrum_smoothing`.
fn third_of_24_octave_bands() -> Vec<(f64, f64)> {
    let g = det::pow(10.0, 3.0 / 10.0);
    let fr = 1000.0;
    let b = 24.0;
    let mut bands: Vec<(f64, f64)> = Vec::new();
    let mut x = -1000.0f64;
    let mut f2 = 0.0f64;
    while f2 <= F_HIGH {
        // b = 24 is EVEN: the reference's even-b midband formula
        // (2x - 59)/(2b) applies (porting the odd-b branch shifted
        // every band edge and rewired the candidate set — executed).
        let fm = det::pow(g, (2.0 * x - 59.0) / (2.0 * b)) * fr;
        let f1 = det::pow(g, -1.0 / (2.0 * b)) * fm;
        f2 = det::pow(g, 1.0 / (2.0 * b)) * fm;
        if f2 >= 90.0 {
            bands.push((f1, f2));
        }
        x += 1.0;
    }
    if let Some(last) = bands.last_mut() {
        last.1 = F_HIGH;
    }
    if let Some(first) = bands.first_mut() {
        first.0 = 90.0;
    }
    bands
}

/// The reference's `_spectrum_smoothing` for one segment: per-band
/// energy mean painted piecewise-constant onto the frequency axis.
/// Faithful details: masked bin 0 is EXCLUDED from every band's mean
/// (the reference's `bin_index > stop - nperseg` filter). Honest
/// divergences from the reference, both disclosed: an EMPTY
/// 1/24-octave band REFUSES (`DegenerateSignal`) instead of entering
/// the reference's corrupted empty-band bookkeeping (a
/// double-decrement loop over a mutated array — review-executed at
/// 96 kHz / 8192 samples, where the reference and any completion of
/// it genuinely diverge), and the bins BEFORE the first painted band
/// and AFTER the last one — numpy `empty` uninitialized memory in
/// the reference — carry the first/last band's value
/// (deterministic completion).
fn smooth_spectrum(freqs: &[f64], spec_db: &[f64]) -> Result<Vec<f64>, PsychoError> {
    let bands = third_of_24_octave_bands();
    let mut vals: Vec<(f64, f64, f64)> = Vec::new(); // (f1, f2, value)
    for &(f1, f2) in &bands {
        let mut sum = 0.0f64;
        let mut count = 0usize;
        for (i, (f, s)) in freqs.iter().zip(spec_db).enumerate() {
            if i > 0 && *f >= f1 && *f <= f2 {
                sum += db_to_pow(*s);
                count += 1;
            }
        }
        if count == 0 {
            return Err(PsychoError::DegenerateSignal {
                what: "frequency resolution too coarse for 1/24-octave smoothing (empty band)",
            });
        }
        vals.push((f1, f2, pow_to_db(sum / usize_f64(count))));
    }
    let mut smooth = vec![0.0f64; freqs.len()];
    for &(f1, f2, v) in &vals {
        let low = argmin_abs(freqs, f1);
        let high = argmin_abs(freqs, f2);
        for s in &mut smooth[low..high] {
            *s = v;
        }
    }
    if let Some(&(f1, _, v)) = vals.first() {
        let low = argmin_abs(freqs, f1);
        for s in &mut smooth[..low] {
            *s = v;
        }
    }
    if let Some(&(_, f2, v)) = vals.last() {
        let high = argmin_abs(freqs, f2);
        for s in &mut smooth[high..] {
            *s = v;
        }
    }
    Ok(smooth)
}

/// ECMA-74 annex D.7.1 lower threshold of hearing (reference `_LTH`).
fn lth(f: f64) -> f64 {
    // Outside 20 Hz..22.05 kHz the reference raises NameError; the
    // detection band makes those unreachable — clamp to the nearest
    // band (disclosed).
    let (fmean, fstd, a1, a2, a3, a4, a5) = if f < 305.0 {
        (
            167.5, 87.3212, 1.415_532, -2.451_068, 1.498_869, -6.983_224, 8.621_226,
        )
    } else if f < 2230.0 {
        (
            1157.5, 488.582, 0.397_994, -0.891_839, -0.815_138, -1.221_319, -7.600_754,
        )
    } else if f < 14_000.0 {
        (
            7250.0,
            3033.25,
            1.584_978,
            -2.766_599,
            -6.906_191_2,
            10.138_553,
            -3.149_339,
        )
    } else {
        (
            16_990.0,
            4049.0,
            -5.775_593,
            -9.200_034,
            26.591_15,
            52.167_12,
            15.615_520_48,
        )
    };
    let ff = (f - fmean) / fstd;
    a1 * ff * ff * ff * ff + a2 * ff * ff * ff + a3 * ff * ff + a4 * ff + a5
}

/// ECMA-74 annex D.8 critical band centered on `f0` (reference
/// `_critical_band`).
fn critical_band(f0: f64) -> (f64, f64) {
    let delta_fc = 25.0 + 75.0 * det::pow(1.0 + 1.4 * det::pow(f0 / 1000.0, 2.0), 0.69);
    if f0 < 500.0 {
        (f0 - delta_fc / 2.0, f0 + delta_fc / 2.0)
    } else {
        let f1 = -delta_fc / 2.0 + det::sqrt(delta_fc * delta_fc + 4.0 * f0 * f0) / 2.0;
        (f1, f1 + delta_fc)
    }
}

/// ECMA-74 annex D.10 contiguous lower critical band.
fn lower_critical_band(f0: f64) -> (f64, f64) {
    let (f2, _) = critical_band(f0);
    let (c0, c1, c2) = if f0 < 171.4 {
        (20.0, 0.0, 0.0)
    } else if f0 <= 1600.0 {
        (-149.5, 1.001, -6.9e-5)
    } else {
        (6.8, 0.806, -8.2e-6)
    };
    (c0 + c1 * f0 + c2 * f0 * f0, f2)
}

/// ECMA-74 annex D.10 contiguous upper critical band.
fn upper_critical_band(f0: f64) -> (f64, f64) {
    let (_, f1) = critical_band(f0);
    let (c0, c1, c2) = if f0 <= 1600.0 {
        (149.5, 1.035, 7.7e-5)
    } else {
        (3.3, 1.215, 2.16e-5)
    };
    (f1, c0 + c1 * f0 + c2 * f0 * f0)
}

/// Peak-level correction (reference `_peak_level`), including its
/// literal `Ltemp - abs(spec[temp]) > 0` walk condition (the abs of a
/// dB value — a reference quirk, replicated).
fn peak_level(spec: &[f64], peak_index: usize) -> f64 {
    let li = spec[peak_index];
    let mut l = spec[peak_index];
    let m = spec.len();
    // Right walk.
    let mut temp = peak_index + 1;
    if temp != m {
        let mut ltemp = li;
        while ltemp - spec[temp].abs() > 0.0 {
            if li - spec[temp] < 10.0 {
                ltemp = spec[temp];
                l = pow_to_db(db_to_pow(l) + db_to_pow(spec[temp]));
                temp += 1;
                if temp == m {
                    temp -= 1;
                    ltemp = -1.0;
                }
            } else {
                ltemp = -1.0;
            }
        }
    }
    // Left walk.
    if peak_index > 0 {
        let mut temp = peak_index - 1;
        let mut ltemp = li;
        while ltemp - spec[temp].abs() > 0.0 {
            if li - spec[temp] < 10.0 {
                ltemp = spec[temp];
                l = pow_to_db(db_to_pow(l) + db_to_pow(spec[temp]));
                if temp == 0 {
                    ltemp = -1.0;
                } else {
                    temp -= 1;
                }
            } else {
                ltemp = -1.0;
            }
        }
    }
    l
}

/// Reference `_find_highest_tone`: the two highest candidates within
/// the critical band centered on `index_list`'s tone `ind`, deleting
/// the lower ones (recursive re-centering exactly as the reference).
fn find_highest_tone(
    freqs: &[f64],
    spec: &[f64],
    index: &mut Vec<usize>,
    ind: usize,
) -> (usize, Option<usize>) {
    let f = freqs[ind];
    let (f1, f2) = critical_band(f);
    let low_limit_idx = argmin_abs(freqs, f1);
    let high_limit_idx = argmin_abs(freqs, f2);

    let multiple_idx: Vec<usize> = index
        .iter()
        .copied()
        .filter(|&i| i > low_limit_idx && i < high_limit_idx)
        .collect();

    if multiple_idx.len() > 1 {
        // Descending sort by level. numpy's default argsort is an
        // UNSTABLE introsort; Rust's sort_by is stable — divergence
        // needs exactly-tied dB levels inside one critical band
        // (practically unreachable; disclosed).
        let mut order: Vec<usize> = (0..multiple_idx.len()).collect();
        order.sort_by(|&a, &b| {
            spec[multiple_idx[b]]
                .partial_cmp(&spec[multiple_idx[a]])
                .expect("finite dB")
        });
        let ind_p = multiple_idx[order[0]];
        let ind_s = multiple_idx[order[1]];
        for &s in &order[2..] {
            index.retain(|&i| i != multiple_idx[s]);
        }
        if ind_p == ind {
            (ind_p, Some(ind_s))
        } else {
            find_highest_tone(freqs, spec, index, ind_p)
        }
    } else {
        (ind, None)
    }
}

/// Reference `_screening_for_tones` (the "smoothed" Bray method,
/// single segment): local maxima 6 dB over the 1/24-octave smoothed
/// spectrum and 10 dB over the hearing threshold, then the tonal
/// width check against half..full critical bandwidth.
fn screening_for_tones(freqs: &[f64], spec_db: &[f64]) -> Result<Vec<usize>, PsychoError> {
    let m = spec_db.len();
    let smooth = smooth_spectrum(freqs, spec_db)?;
    // Criteria 1-3.
    let mut index: Vec<usize> = Vec::new();
    for j in 1..m - 1 {
        if sign(spec_db[j + 1] - spec_db[j]) - sign(spec_db[j] - spec_db[j - 1]) < 0.0
            && spec_db[j] > smooth[j] + 6.0
            && spec_db[j] > lth(freqs[j]) + 10.0
        {
            index.push(j);
        }
    }
    // Width check.
    let mut tones: Vec<usize> = Vec::new();
    while let Some(&first) = index.first() {
        let mut peak_index = first;
        // `low_limit` is signed: it starts at the ORIGINAL candidate
        // but the left walk starts from the right-walk-MIGRATED peak
        // and can take more steps than `low_limit`'s starting value
        // (review-executed panic: usize underflow at coarse
        // resolutions). The reference lets it go negative and numpy
        // WRAPS `freqs[low_limit]` into a garbage width; here the
        // width read clamps at the band edge instead — a disclosed
        // deviation on a reference-buggy path.
        let mut low_limit = peak_index.cast_signed();
        let mut high_limit = peak_index;
        // Right walk.
        let mut temp = peak_index + 1;
        while temp < m && spec_db[temp] > smooth[temp] + 6.0 && temp + 1 < m {
            if spec_db[temp] > spec_db[peak_index] {
                peak_index = temp;
            }
            high_limit += 1;
            temp += 1;
        }
        // Left walk (the reference's wrap-read at -1 is dead: the
        // bound test fails regardless of the wrapped comparison).
        let mut temp = peak_index.cast_signed() - 1;
        while temp >= 0 && spec_db[temp.cast_unsigned()] > smooth[temp.cast_unsigned()] + 6.0 {
            if spec_db[temp.cast_unsigned()] > spec_db[peak_index] {
                peak_index = temp.cast_unsigned();
            }
            low_limit -= 1;
            temp -= 1;
        }
        let (f1, f2) = critical_band(freqs[peak_index]);
        let cb_width = f2 - f1;
        let t_width = freqs[high_limit] - freqs[low_limit.max(0).cast_unsigned()];
        if t_width < cb_width {
            tones.push(peak_index);
        }
        index.retain(|&i| i > high_limit);
    }
    Ok(tones)
}

/// Validate + build the masked detection-band view.
fn detection_band(pcm_pa: &[f64], sample_rate: f64) -> Result<(Vec<f64>, Vec<f64>), PsychoError> {
    let n = pcm_pa.len();
    if !n.is_power_of_two() || n < 4096 {
        return Err(PsychoError::Shape {
            what: "tonality needs a power-of-two block of at least 4096 samples",
        });
    }
    if !sample_rate.is_finite() || sample_rate <= 0.0 {
        return Err(PsychoError::NonFinite {
            what: "sample rate",
        });
    }
    if pcm_pa.iter().any(|v| !v.is_finite()) {
        return Err(PsychoError::NonFinite { what: "pcm sample" });
    }
    let (db, freq) = spectrum_db(pcm_pa, sample_rate);
    let lo = freq.partition_point(|&f| f <= F_LOW);
    let hi = freq.partition_point(|&f| f < F_HIGH);
    if hi.saturating_sub(lo) < 64 {
        return Err(PsychoError::DegenerateSignal {
            what: "detection band (89.1..11200 Hz) has too few spectral lines",
        });
    }
    Ok((freq[lo..hi].to_vec(), db[lo..hi].to_vec()))
}

/// Shared: record a tone with the method's prominence rule and
/// accumulate.
struct ToneAccum {
    tones: Vec<Tone>,
    prominent_pow: f64,
}

impl ToneAccum {
    fn push(&mut self, f: f64, ratio: f64, prominent: bool) {
        if prominent {
            self.prominent_pow += db_to_pow(ratio);
        }
        self.tones.push(Tone {
            frequency_hz: f,
            ratio_db: ratio,
            prominent,
        });
    }

    fn finish(self) -> Tonality {
        let total_db = if self.prominent_pow > 0.0 {
            pow_to_db(self.prominent_pow)
        } else {
            0.0
        };
        Tonality {
            tones: self.tones,
            total_db,
        }
    }
}

/// Tone-to-noise ratio per ECMA-74 annex D (reference
/// `_tnr_main_calc`, single segment), with T-TNR over the prominent
/// tones per ECMA TR/108.
///
/// # Errors
/// [`PsychoError::Shape`] on a non-power-of-two or too-short block;
/// [`PsychoError::NonFinite`] on bad samples/rate;
/// [`PsychoError::DegenerateSignal`] when the detection band is
/// unresolvable at this rate/length.
#[allow(clippy::too_many_lines)] // one reference routine, told in its order
pub fn tone_to_noise_ecma(pcm_pa: &[f64], sample_rate: f64) -> Result<Tonality, PsychoError> {
    let (fr, spec) = detection_band(pcm_pa, sample_rate)?;
    let mut peaks = screening_for_tones(&fr, &spec)?;
    let mut acc = ToneAccum {
        tones: Vec::new(),
        prominent_pow: 0.0,
    };
    // The reference indexes the UNMASKED axis here with masked
    // indices (its own bug); on a uniform axis the DIFFERENCES it
    // takes are identical, so the masked axis is used with no
    // numerical effect.
    let df = fr[1] - fr[0];
    while !peaks.is_empty() {
        let ind = peaks[0];
        let (ind_p, ind_s) = if peaks.len() > 1 {
            find_highest_tone(&fr, &spec, &mut peaks, ind)
        } else {
            (ind, None)
        };
        let (lt, low_limit_idx, high_limit_idx, delta_fc, delta_ft) = if let Some(ind_s) = ind_s {
            let fp = fr[ind_p];
            let fs_ = fr[ind_s];
            // Proximity criterion.
            let delta_f = 21.0
                * det::pow(
                    10.0,
                    det::pow(1.2 * det::ln(fp / 212.0).abs() / det::ln(10.0), 1.8),
                );
            if (fs_ - fp).abs() < delta_f {
                let lp = peak_level(&spec, ind_p);
                let ls = peak_level(&spec, ind_s);
                let lt = pow_to_db(db_to_pow(lp) + db_to_pow(ls));
                let (f1, f2) = critical_band(fp);
                peaks.retain(|&i| i != ind_s);
                (
                    lt,
                    argmin_abs(&fr, f1),
                    argmin_abs(&fr, f2),
                    f2 - f1,
                    2.0 * df,
                )
            } else {
                let (f1, f2) = critical_band(fr[ind_p]);
                (
                    spec[ind_p],
                    argmin_abs(&fr, f1),
                    argmin_abs(&fr, f2),
                    f2 - f1,
                    df,
                )
            }
        } else {
            let (f1, f2) = critical_band(fr[ind_p]);
            (
                peak_level(&spec, ind_p),
                argmin_abs(&fr, f1),
                argmin_abs(&fr, f2),
                f2 - f1,
                df,
            )
        };
        let spec_sum: f64 = spec[low_limit_idx..high_limit_idx]
            .iter()
            .map(|&s| db_to_pow(s))
            .sum();
        let ltot = pow_to_db(spec_sum);
        let delta_ftot = fr[high_limit_idx] - fr[low_limit_idx];
        let ln_noise = pow_to_db(db_to_pow(ltot) - db_to_pow(lt))
            + pow_to_db(delta_fc / (delta_ftot - delta_ft));
        let f = fr[ind_p];
        let delta_t = lt - ln_noise;
        if delta_t > 0.0 {
            let prominent = if f < 1000.0 {
                delta_t >= 8.0 + 8.33 * det::ln(1000.0 / f) / det::ln(10.0)
            } else {
                delta_t >= 8.0
            };
            acc.push(f, delta_t, prominent);
        }
        peaks.retain(|&i| i != ind_p);
    }
    Ok(acc.finish())
}

/// Prominence ratio per ECMA-74 annex D (reference `_pr_main_calc`,
/// single segment), with T-PR over the prominent tones.
///
/// # Errors
/// Same refusal set as [`tone_to_noise_ecma`].
pub fn prominence_ratio_ecma(pcm_pa: &[f64], sample_rate: f64) -> Result<Tonality, PsychoError> {
    let (fr, spec) = detection_band(pcm_pa, sample_rate)?;
    let mut peaks = screening_for_tones(&fr, &spec)?;
    let mut acc = ToneAccum {
        tones: Vec::new(),
        prominent_pow: 0.0,
    };
    let band_level = |lo: usize, hi: usize| -> f64 {
        let s: f64 = spec[lo..hi].iter().map(|&v| db_to_pow(v)).sum();
        if s == 0.0 { 0.0 } else { pow_to_db(s) }
    };
    while !peaks.is_empty() {
        let ind = peaks[0];
        let ind = if peaks.len() > 1 {
            find_highest_tone(&fr, &spec, &mut peaks, ind).0
        } else {
            ind
        };
        let ft = fr[ind];
        // Middle critical band.
        let (f1, f2) = critical_band(ft);
        let low_limit_idx = argmin_abs(&fr, f1);
        let high_limit_idx = argmin_abs(&fr, f2);
        let lm = band_level(low_limit_idx, high_limit_idx);
        // Lower contiguous band.
        let (f1, f2) = lower_critical_band(ft);
        let ll = band_level(argmin_abs(&fr, f1), argmin_abs(&fr, f2));
        let delta_f = f2 - f1;
        // Upper contiguous band.
        let (f1, f2) = upper_critical_band(ft);
        let lu = band_level(argmin_abs(&fr, f1), argmin_abs(&fr, f2));
        let delta = if ft <= 171.4 {
            pow_to_db(db_to_pow(lm))
                - pow_to_db(f64::midpoint(
                    (100.0 / delta_f) * db_to_pow(ll),
                    db_to_pow(lu),
                ))
        } else {
            pow_to_db(db_to_pow(lm)) - pow_to_db(f64::midpoint(db_to_pow(ll), db_to_pow(lu)))
        };
        if delta > 0.0 {
            let prominent = if ft <= 1000.0 {
                delta >= 9.0 + 10.0 * det::ln(1000.0 / ft) / det::ln(10.0)
            } else {
                delta >= 9.0
            };
            acc.push(ft, delta, prominent);
        }
        // Remove every candidate inside the middle critical band.
        peaks.retain(|&i| i < low_limit_idx || i > high_limit_idx);
    }
    Ok(acc.finish())
}
