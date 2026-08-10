//! Fitted noise-source tables — the product deliverable of bead
//! 9ok02: dipole source SPECTRAL SHAPES versus blowing velocity,
//! extracted from jet-labium runs, exported with provenance and the
//! scope statement, and consumed by a demo synthesizer.
//!
//! Everything here is a SHAPE/SCALING authority
//! ([`crate::SCOPE_STATEMENT`]): spectra are stored as RELATIVE dB
//! over log-spaced Strouhal bands, normalized per entry; the
//! velocity scaling of total radiated dipole power is reported as a
//! MEASURED exponent, never asserted equal to a textbook value.

use crate::jetlab::{JetLabiumConfig, run_jet_labium};
use crate::{AeroacError, SCOPE_STATEMENT};
use fs_math::det;

/// Log-spaced Strouhal bands per table entry.
pub const N_BANDS: usize = 16;
/// Band range in Strouhal (f * delta / u).
const ST_LO: f64 = 0.02;
const ST_HI: f64 = 4.0;

/// One table row: the jet condition and its dipole spectral shape.
#[derive(Debug, Clone)]
pub struct NoiseEntry {
    /// Jet peak velocity [lu/step] (the blowing-pressure proxy).
    pub u_jet: f64,
    /// Jet Reynolds number.
    pub reynolds: f64,
    /// Max Mach diagnostic of the run.
    pub mach_max_lattice: f64,
    /// Plate-vs-fringe flux imbalance of the run.
    pub flux_imbalance: f64,
    /// Transverse-force RMS [lattice units] WITHIN the table's
    /// Strouhal band range (Parseval over the banded periodogram) —
    /// the oscillatory dipole strength, excluding out-of-band drift
    /// (executed: total RMS was drift-dominated and read a flat
    /// velocity scaling).
    pub force_rms: f64,
    /// Relative band power DENSITY levels [dB] (band power divided
    /// by the band's bin count — density, so the shape is
    /// record-length-independent; a band-SUM convention broke the
    /// synth round trip by the log-band width ratio, executed),
    /// normalized so the strongest band is 0 dB; bands are
    /// log-spaced in Strouhal over [`ST_LO`], [`ST_HI`].
    pub band_db: [f64; N_BANDS],
}

/// The fitted table: entries over a velocity sweep plus the measured
/// power-law exponent, provenance, and the scope law.
#[derive(Debug, Clone)]
pub struct NoiseTable {
    /// Sweep entries, ascending in `u_jet`.
    pub entries: Vec<NoiseEntry>,
    /// Least-squares slope of `ln(force_rms^2)` vs `ln(u_jet)` — the
    /// MEASURED dipole-power scaling exponent (reported, not
    /// prescribed).
    pub power_exponent: f64,
    /// Geometry provenance (the config of the first entry, which
    /// fixes everything but `u_jet`).
    pub geometry: JetLabiumConfig,
    /// The honest-scope statement (the marketing-mutation guard
    /// asserts its presence in every export).
    pub scope: &'static str,
}

/// Band center Strouhal values (log-spaced).
#[must_use]
pub fn band_centers() -> [f64; N_BANDS] {
    let mut c = [0.0; N_BANDS];
    let ratio = det::ln(ST_HI / ST_LO);
    for (i, v) in c.iter_mut().enumerate() {
        #[allow(clippy::cast_precision_loss)]
        let t = (i as f64 + 0.5) / N_BANDS as f64;
        *v = ST_LO * det::exp(ratio * t);
    }
    c
}

/// Band index for a Strouhal value (edge-based log spacing) — public
/// so consumers analyze with the SAME folding the synth applies.
#[must_use]
pub fn band_of(st: f64) -> Option<usize> {
    if st < ST_LO || st >= ST_HI {
        return None;
    }
    let t = det::ln(st / ST_LO) / det::ln(ST_HI / ST_LO);
    #[allow(
        clippy::cast_precision_loss,
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss
    )]
    Some(((t * N_BANDS as f64) as usize).min(N_BANDS - 1))
}

/// Build the table by sweeping `u_values` over the fixed `geometry`
/// (its `u_jet` field is overwritten per entry).
///
/// # Errors
/// Propagates run refusals; refuses fewer than 2 velocities,
/// non-ascending sweeps, or a run whose diagnostics leave the
/// low-Mach regime (`mach > 0.25`) — a table row from an
/// out-of-regime run would be a fabricated authority.
pub fn fit_noise_table(
    geometry: &JetLabiumConfig,
    u_values: &[f64],
) -> Result<NoiseTable, AeroacError> {
    if u_values.len() < 2 || u_values.windows(2).any(|w| w[1] <= w[0]) {
        return Err(AeroacError::InvalidParameter {
            what: "velocity sweep must be ascending with at least 2 points",
        });
    }
    let mut entries = Vec::with_capacity(u_values.len());
    for &u in u_values {
        let mut cfg = geometry.clone();
        cfg.u_jet = u;
        let run = run_jet_labium(&cfg)?;
        let d = &run.diagnostics;
        if d.mach_max_lattice > 0.25 {
            return Err(AeroacError::InvalidParameter {
                what: "sweep point left the low-Mach regime (diagnostic gate)",
            });
        }
        let n = run.force_series.len();
        #[allow(clippy::cast_precision_loss)]
        let nf = n as f64;
        let mean = run.force_series.iter().map(|f| f[1]).sum::<f64>() / nf;
        // Band power via a Hann-windowed periodogram folded into the
        // Strouhal bands.
        let fft = fs_fft::Fft::new(n);
        let mut buf: Vec<fs_fft::C64> = run
            .force_series
            .iter()
            .enumerate()
            .map(|(i, f)| {
                #[allow(clippy::cast_precision_loss)]
                let w = 0.5 - 0.5 * det::cos(2.0 * core::f64::consts::PI * i as f64 / (nf - 1.0));
                fs_fft::C64::new((f[1] - mean) * w, 0.0)
            })
            .collect();
        let mut scratch = vec![fs_fft::C64::new(0.0, 0.0); n];
        fft.forward(&mut buf, &mut scratch);
        let delta = 2.0 * cfg.slot_half;
        let mut band_pow = [0.0f64; N_BANDS];
        let mut band_bins = [0usize; N_BANDS];
        for (k, c) in buf[..n / 2].iter().enumerate().skip(1) {
            #[allow(clippy::cast_precision_loss)]
            let st = (k as f64 / nf) * delta / u;
            if let Some(b) = band_of(st) {
                band_pow[b] += c.norm_sq();
                band_bins[b] += 1;
            }
        }
        // Density: divide each band by its bin count (empty bands
        // stay zero and read the -30 dB floor below).
        for (p, &m) in band_pow.iter_mut().zip(&band_bins) {
            if m > 0 {
                #[allow(clippy::cast_precision_loss)]
                let mf = m as f64;
                *p /= mf;
            }
        }
        let peak = band_pow.iter().copied().fold(f64::MIN, f64::max);
        if !peak.is_finite() || peak <= 0.0 {
            return Err(AeroacError::InvalidParameter {
                what: "no band energy in the Strouhal range (degenerate run)",
            });
        }
        // Band-limited RMS via Parseval (relative units; the Hann
        // window's coherent-gain factor is common to every entry so
        // the SCALING claim is unaffected).
        let force_rms = det::sqrt(band_pow.iter().sum::<f64>() / (nf * nf));
        let mut band_db = [0.0f64; N_BANDS];
        for (db, p) in band_db.iter_mut().zip(&band_pow) {
            *db = 10.0 * det::ln((p / peak).max(1e-30)) / det::ln(10.0);
        }
        entries.push(NoiseEntry {
            u_jet: u,
            reynolds: d.reynolds,
            mach_max_lattice: d.mach_max_lattice,
            flux_imbalance: (d.flux_plate_plane - d.flux_fringe_plane).abs()
                / d.flux_plate_plane.abs(),
            force_rms,
            band_db,
        });
    }
    // Measured power exponent: slope of ln(rms^2) vs ln(u).
    #[allow(clippy::cast_precision_loss)]
    let np = entries.len() as f64;
    let xs: Vec<f64> = entries.iter().map(|e| det::ln(e.u_jet)).collect();
    let ys: Vec<f64> = entries.iter().map(|e| 2.0 * det::ln(e.force_rms)).collect();
    let sx: f64 = xs.iter().sum();
    let sy: f64 = ys.iter().sum();
    let sxx: f64 = xs.iter().map(|x| x * x).sum();
    let sxy: f64 = xs.iter().zip(&ys).map(|(x, y)| x * y).sum();
    let power_exponent = (np * sxy - sx * sy) / (np * sxx - sx * sx);
    Ok(NoiseTable {
        entries,
        power_exponent,
        geometry: geometry.clone(),
        scope: SCOPE_STATEMENT,
    })
}

impl NoiseTable {
    /// Serialize to a JSON export (hand-rolled, no schema claims):
    /// entries, measured exponent, geometry provenance, and the
    /// scope statement — the marketing-mutation guard asserts the
    /// latter's presence.
    #[must_use]
    #[allow(clippy::format_push_string)] // export builder, clarity over micro-alloc
    pub fn to_json(&self) -> String {
        let mut s = String::from("{\"kind\":\"fs-aeroac-noise-table\",\"entries\":[");
        for (i, e) in self.entries.iter().enumerate() {
            if i > 0 {
                s.push(',');
            }
            s.push_str(&format!(
                "{{\"u_jet\":{:.6},\"reynolds\":{:.1},\"mach_max\":{:.4},\"flux_imbalance\":{:.6},\"force_rms\":{:e},\"band_db\":{:?}}}",
                e.u_jet, e.reynolds, e.mach_max_lattice, e.flux_imbalance, e.force_rms, e.band_db
            ));
        }
        s.push_str(&format!(
            "],\"power_exponent_measured\":{:.3},\"geometry\":{{\"nx\":{},\"ny\":{},\"slot_half\":{},\"edge_distance\":{},\"nozzle_thickness\":{}}},\"scope\":\"{}\"}}",
            self.power_exponent,
            self.geometry.nx,
            self.geometry.ny,
            self.geometry.slot_half,
            self.geometry.edge_distance,
            self.geometry.nozzle_thickness,
            self.scope
        ));
        s
    }

    /// The demo-synth consumer: shape deterministic broadband noise
    /// to the table entry nearest `u_jet`, returning `n` samples
    /// whose smoothed spectrum follows the entry's band shape
    /// (RELATIVE shape only — the output is dimensionless and the
    /// scope law applies).
    ///
    /// # Errors
    /// [`AeroacError::InvalidParameter`] on an empty table,
    /// non-power-of-two `n`, or `u_jet` outside the sweep range
    /// (extrapolating a fitted table is a silent claim — refused).
    pub fn synthesize(&self, u_jet: f64, n: usize, seed: u32) -> Result<Vec<f64>, AeroacError> {
        if self.entries.is_empty() {
            return Err(AeroacError::InvalidParameter {
                what: "empty table",
            });
        }
        if !n.is_power_of_two() || n < 256 {
            return Err(AeroacError::InvalidParameter {
                what: "synthesis length must be a power of two >= 256",
            });
        }
        let lo = self.entries.first().expect("nonempty").u_jet;
        let hi = self.entries.last().expect("nonempty").u_jet;
        if !(lo..=hi).contains(&u_jet) {
            return Err(AeroacError::InvalidParameter {
                what: "u_jet outside the fitted sweep (no extrapolation)",
            });
        }
        let entry = self
            .entries
            .iter()
            .min_by(|a, b| {
                (a.u_jet - u_jet)
                    .abs()
                    .partial_cmp(&(b.u_jet - u_jet).abs())
                    .expect("finite")
            })
            .expect("nonempty");
        // White noise from a splitmix64-style hash (an LCG's
        // lattice-line correlations are NOT spectrally white —
        // executed: band-level errors of ~7 dB in the round-trip),
        // shaped in the frequency domain by the band gains.
        let mut noise: Vec<fs_fft::C64> = (0..n)
            .map(|i| {
                let mut z = (u64::from(seed) << 32)
                    .wrapping_add(i as u64)
                    .wrapping_add(0x9E37_79B9_7F4A_7C15);
                z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
                z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
                z ^= z >> 31;
                #[allow(clippy::cast_precision_loss)]
                fs_fft::C64::new((z >> 11) as f64 / 9_007_199_254_740_992.0 - 0.5, 0.0)
            })
            .collect();
        let fft = fs_fft::Fft::new(n);
        let mut scratch = vec![fs_fft::C64::new(0.0, 0.0); n];
        fft.forward(&mut noise, &mut scratch);
        let delta = 2.0 * self.geometry.slot_half;
        #[allow(clippy::cast_precision_loss)]
        let nf = n as f64;
        for (k, c) in noise.iter_mut().enumerate() {
            let kk = if k <= n / 2 { k } else { n - k };
            #[allow(clippy::cast_precision_loss)]
            let st = (kk as f64 / nf) * delta / entry.u_jet;
            let gain = band_of(st).map_or(0.0, |b| det::pow(10.0, entry.band_db[b] / 20.0));
            *c = fs_fft::C64::new(c.re * gain, c.im * gain);
        }
        fft.inverse(&mut noise, &mut scratch);
        #[allow(clippy::cast_precision_loss)]
        Ok(noise.iter().map(|c| c.re * nf).collect())
    }
}
