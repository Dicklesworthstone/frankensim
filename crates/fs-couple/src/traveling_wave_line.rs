//! TOMBSTONE-3ez8g.13.5 — DO NOT WIRE NEW CONSUMERS.
//!
//! This one-pole `TravelingWaveLine` is DEAD CODE, superseded by the
//! exact-FIR characteristic line (`DelayedFilter` + `driving_point`).
//! It has ZERO in-tree consumers (swept at disposition time) and the
//! crate CONTRACT negates it twice: "There is no one-pole
//! TravelingWaveLine fallback on those paths" and "The one-pole
//! TravelingWaveLine is not that path". The program NEVER list bans
//! resurrecting it; this header makes the ban physical (music bead
//! `frankensim-music-v8-root-3ez8g.13.5`, tombstone option — deletion
//! requires explicit written owner permission and remains available).
//!
//! Every public item is `#[deprecated]` and `#[doc(hidden)]`: wiring
//! it in is LOUD. If you believe you need this module, you need the
//! exact-FIR characteristic line instead.
//!
//! Original description (retained for the record): a 1-D
//! traveling-wave line — fractional delay of `2L/c` with a reflection
//! gain from a termination and a mild viscothermal one-pole; the
//! frequency-domain law lives in fs-duct.

use fs_duct::{Duct, DuctError, LossModel, Termination, input_impedance, segment_wave};
use fs_material::gas::GasState;
use fs_math::det;

/// Typed refusal from a traveling-wave line.
#[deprecated(note = "TOMBSTONE-3ez8g.13.5: dead code superseded by the exact-FIR \
            characteristic line (DelayedFilter + driving_point); the CONTRACT \
            forbids one-pole fallbacks on the characteristic paths")]
#[doc(hidden)]
#[derive(Debug, Clone, PartialEq)]
pub enum TravelingWaveError {
    /// Geometry or delay is not realizable.
    Invalid {
        /// Which check failed.
        what: &'static str,
    },
    /// TMM refusal.
    Duct(DuctError),
}

#[allow(deprecated)]
impl core::fmt::Display for TravelingWaveError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid { what } => write!(f, "FS-COUPLE-LINE: {what}"),
            Self::Duct(e) => write!(f, "{e}"),
        }
    }
}

#[allow(deprecated)]
impl std::error::Error for TravelingWaveError {}

/// Fractional-delay traveling-wave line.
#[deprecated(note = "TOMBSTONE-3ez8g.13.5: dead code superseded by the exact-FIR \
            characteristic line (DelayedFilter + driving_point); the CONTRACT \
            forbids one-pole fallbacks on the characteristic paths")]
#[doc(hidden)]
#[derive(Debug, Clone)]
pub struct TravelingWaveLine {
    buf: Vec<f64>,
    write: usize,
    delay_int: usize,
    frac: f64,
    gain: f64,
    pole: f64,
    y: f64,
}

#[allow(deprecated)]
impl TravelingWaveLine {
    /// Build a round-trip line from a TMM duct.
    ///
    /// Gain is `±|R|` at the quarter-wave from Bessel `Z_in`.
    /// The one-pole is mild HF leak from `Im k` at `f0` vs `3 f0`.
    ///
    /// # Errors
    /// Empty geometry, a delay that does not fit `n` samples, or TMM.
    pub fn from_duct(
        physics: &Duct,
        gas: &GasState,
        termination: Termination,
        sample_rate_hz: u32,
        n: usize,
        zc: f64,
    ) -> Result<Self, TravelingWaveError> {
        let length: f64 = physics
            .segments
            .iter()
            .map(|s| match *s {
                fs_duct::Segment::Cylinder { length, .. }
                | fs_duct::Segment::Cone { length, .. } => length,
                fs_duct::Segment::ToneHole { .. } => 0.0,
            })
            .sum();
        if !(length > 0.0) {
            return Err(TravelingWaveError::Invalid {
                what: "line needs positive axial length",
            });
        }
        let dt = 1.0 / f64::from(sample_rate_hz);
        let delay_f = 2.0 * length / gas.sound_speed / dt;
        if !(delay_f >= 2.0 && delay_f < n as f64 - 2.0) {
            return Err(TravelingWaveError::Invalid {
                what: "round-trip delay does not fit the realized history",
            });
        }
        let delay_int = delay_f.floor() as usize;
        let frac = delay_f - delay_int as f64;
        let omega0 = core::f64::consts::PI * gas.sound_speed / (2.0 * length);
        let z = input_impedance(physics, gas, omega0, LossModel::Bessel, termination)
            .map_err(TravelingWaveError::Duct)?
            .impedance;
        let mag = reflection_magnitude(z.re, z.im, zc).clamp(0.2, 0.99);
        let sign = match termination {
            Termination::Closed => 1.0,
            _ => -1.0,
        };
        let radius = physics
            .segments
            .first()
            .map_or(0.007, fs_duct::Segment::outlet_radius);
        let k1 = segment_wave(gas, radius, omega0, LossModel::Bessel)
            .map_err(TravelingWaveError::Duct)?;
        let k3 = segment_wave(gas, radius, 3.0 * omega0, LossModel::Bessel)
            .map_err(TravelingWaveError::Duct)?;
        let att1 = det::exp(-2.0 * k1.wavenumber.im * length);
        let att3 = det::exp(-2.0 * k3.wavenumber.im * length);
        let ratio = (att3 / att1.max(1.0e-9)).clamp(0.4, 1.0);
        let pole = (1.0 - ratio).clamp(0.0, 0.35);
        Ok(Self {
            buf: vec![0.0; delay_int + 4],
            write: 0,
            delay_int,
            frac,
            gain: sign * mag,
            pole,
            y: 0.0,
        })
    }

    /// Inject an outgoing wave and return the new incoming wave.
    pub fn push(&mut self, outgoing: f64) -> f64 {
        self.buf[self.write] = outgoing;
        let n = self.buf.len();
        let i1 = (self.write + n - self.delay_int) % n;
        let i0 = (i1 + n - 1) % n;
        let delayed = (1.0 - self.frac) * self.buf[i1] + self.frac * self.buf[i0];
        self.y = self.gain * delayed * (1.0 - self.pole) + self.pole * self.y;
        self.write = (self.write + 1) % n;
        self.y
    }

    /// Incoming wave from the previous push.
    #[must_use]
    pub fn incoming(&self) -> f64 {
        self.y
    }
}

fn reflection_magnitude(zr: f64, zi: f64, zc: f64) -> f64 {
    let dr = zr + zc;
    let di = zi;
    let den = dr * dr + di * di;
    if den < 1.0e-30 {
        return 0.0;
    }
    let nr = zr - zc;
    let re = (nr * dr + zi * di) / den;
    let im = (zi * dr - nr * di) / den;
    det::sqrt(re * re + im * im)
}
