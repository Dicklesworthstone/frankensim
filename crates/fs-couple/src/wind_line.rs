//! Wind articulation runtime (music bead
//! `frankensim-music-v8-root-3ez8g.6.3`): per-fingering plane-wave
//! characteristic lines with the carry-history + crossfade switch (the
//! same D17 lift the brass MM bank ships), a register-vent lane, and
//! the char -> VFIT hop with a REPLAY-PRIME lift — plus the machinery
//! that makes "when to hop" a MEASURED policy instead of doctrine.
//!
//! THE HOP: a settled note may leave the exact-FIR characteristic line
//! for the cheaper VFIT hold image (a rational fit of the same
//! fingering's reflectance). The D17 lift is REPLAY PRIMING: the
//! incoming IIR line is fed the outgoing-wave history the FIR line
//! already carries (outputs discarded), so its internal state matches
//! the recent signal before the crossfade begins. Hop-readiness is
//! GATED: hopping to an image whose registry gate is not green is
//! structurally refused ([`WindLineError::ImageNotGated`]) — an
//! ungated image cannot be entered mid-performance, by construction.
//!
//! Policy honesty: the settle detector's parameters and the measured
//! click-vs-hop-timing table are DATA (the committed
//! `data/claims/wind-hop-policy.tsv` artifact), consumable by the
//! selector; the fixture that mints them keeps its deliberately-early
//! hop as the falsifier showing the click the policy exists to avoid.

use fs_duct::{Duct, FingeringTable, LossModel, Termination, impedance_sweep};
use fs_material::gas::GasState;
use fs_math::c64::C64;
use fs_vfit::discretize::{
    DelayedFilter, DiscretizeError, modulate_delay, realize_tabulated, reflectance,
};
use fs_vfit::vf::FitOptions;

use crate::driving_point::{DrivingPointError, characteristic_line};

/// Typed refusals from the wind line runtime.
#[derive(Debug)]
pub enum WindLineError {
    /// A structural parameter is unusable.
    Invalid {
        /// Diagnosis.
        what: &'static str,
    },
    /// The line realization refused.
    Realize(DrivingPointError),
    /// A runtime sample refused.
    Discrete(DiscretizeError),
    /// A hop targeted an image whose registry gate is not green — an
    /// ungated image cannot be entered mid-performance.
    ImageNotGated {
        /// The refused image id.
        image: &'static str,
    },
}

impl core::fmt::Display for WindLineError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            WindLineError::Invalid { what } => write!(f, "FS-COUPLE-WINDLINE: {what}"),
            WindLineError::Realize(e) => write!(f, "FS-COUPLE-WINDLINE realize: {e:?}"),
            WindLineError::Discrete(e) => write!(f, "FS-COUPLE-WINDLINE line: {e:?}"),
            WindLineError::ImageNotGated { image } => write!(
                f,
                "FS-COUPLE-WINDLINE-GATE: image {image:?} is not gate-green; an ungated \
                 image cannot be entered mid-performance"
            ),
        }
    }
}

impl core::error::Error for WindLineError {}

impl From<DrivingPointError> for WindLineError {
    fn from(e: DrivingPointError) -> Self {
        WindLineError::Realize(e)
    }
}

impl From<DiscretizeError> for WindLineError {
    fn from(e: DiscretizeError) -> Self {
        WindLineError::Discrete(e)
    }
}

/// Which reflection image is currently playing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ActiveImage {
    /// The exact-FIR characteristic line (the phrase image).
    CharLine,
    /// The rational VFIT hold image (the settled-note image).
    VfitHold,
}

/// One switch/hop record.
#[derive(Debug, Clone)]
pub struct WindLiftRecord {
    /// What happened.
    pub action: &'static str,
    /// Samples carried / replayed into the incoming state.
    pub carried: usize,
    /// Crossfade length [samples].
    pub fade_samples: usize,
}

/// Same-sample reflection boundary, before the next outgoing wave is committed.
/// `p_minus = direct_reflection * p_plus + history_pressure_pa`; both terms
/// use the NEXT crossfade weight. The history contains no current input.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WindBoundary {
    /// Instantaneous pressure reflectance (dimensionless).
    pub direct_reflection: f64,
    /// Delayed and filtered reflection already determined by wave history [Pa].
    pub history_pressure_pa: f64,
}

/// Per-fingering plane-wave line bank with the char->VFIT hop.
///
/// Internally each line stores the strictly causal remainder one sample ahead:
/// R(z) = r0 + z^-1 H_tail(z). The existing filter owner evaluates H_tail;
/// no filter equation or wave propagation law is duplicated here. A junction
/// uses `boundary_response`, then commits its solved outgoing wave with `push`.
/// The original `incoming` accessor remains a LAST-sample observation.
pub struct WindLineBank {
    fir_lines: Vec<DelayedFilter>,
    fir_taps: Vec<Vec<f64>>,
    vfit_lines: Vec<Option<DelayedFilter>>,
    // Immutable zero-state images: every re-entry primes a fresh state.
    vfit_templates: Vec<Option<DelayedFilter>>,
    labels: Vec<String>,
    active: usize,
    image: ActiveImage,
    vfit_gate_green: bool,
    dt: f64,
    zc: f64,
    fade_old: Option<DelayedFilter>,
    fade_old_direct: f64,
    fade_remaining: usize,
    fade_total: usize,
    history: std::collections::VecDeque<f64>,
    history_capacity: usize,
    last_incoming: f64,
    // A downstream failure after another filter advanced cannot be resumed.
    poisoned: bool,
    lift_log: Vec<WindLiftRecord>,
}

impl WindLineBank {
    /// Realize the exact-FIR characteristic line per fingering of the
    /// typed table. `vfit_gate_green` states whether the VFIT hold
    /// image's registry gate is green (the hop refuses otherwise —
    /// scope-4's structural impossibility).
    ///
    /// # Errors
    /// [`WindLineError`] on realization refusal.
    pub fn new(
        table: &FingeringTable,
        gas: &GasState,
        termination: Termination,
        sample_rate_hz: u32,
        zc: f64,
        vfit_gate_green: bool,
    ) -> Result<WindLineBank, WindLineError> {
        if sample_rate_hz == 0 || !(zc > 0.0 && zc.is_finite()) {
            return Err(WindLineError::Invalid {
                what: "sample rate and zc must be positive",
            });
        }
        let dt = 1.0 / f64::from(sample_rate_hz);
        let mut fir_lines = Vec::new();
        let mut fir_taps = Vec::new();
        let mut labels = Vec::new();
        for label in table.labels() {
            let duct = table.duct(label).map_err(|_| WindLineError::Invalid {
                what: "fingering table refused a label it listed",
            })?;
            let line =
                characteristic_line(&duct, gas, termination, sample_rate_hz, 4096, zc, None)?;
            let taps = probe_taps(&line)?;
            fir_lines.push(history_line(dt, &taps, &[])?);
            fir_taps.push(taps);
            labels.push(label.to_string());
        }
        let n = fir_lines.len();
        if n == 0 {
            return Err(WindLineError::Invalid { what: "wind bank requires at least one fingering" });
        }
        let history_capacity = fir_taps.iter().map(Vec::len).max().unwrap_or(0);
        Ok(WindLineBank {
            fir_lines,
            fir_taps,
            vfit_lines: vec![None; n],
            vfit_templates: vec![None; n],
            labels,
            active: 0,
            image: ActiveImage::CharLine,
            vfit_gate_green,
            dt,
            zc,
            fade_old: None,
            fade_old_direct: 0.0,
            fade_remaining: 0,
            fade_total: 0,
            history: std::collections::VecDeque::with_capacity(history_capacity),
            history_capacity,
            last_incoming: 0.0,
            poisoned: false,
            lift_log: Vec::new(),
        })
    }

    /// The active image.
    #[must_use]
    pub fn image(&self) -> ActiveImage {
        self.image
    }

    /// The active fingering's label (from the [`FingeringTable`]).
    #[must_use]
    pub fn active_label(&self) -> &str {
        &self.labels[self.active]
    }

    /// The lift log.
    #[must_use]
    pub fn lift_log(&self) -> &[WindLiftRecord] {
        &self.lift_log
    }

    fn check_live(&self) -> Result<(), WindLineError> {
        if self.poisoned {
            Err(WindLineError::Invalid { what: "wind line has a partial physical step; reconstruct it" })
        } else { Ok(()) }
    }

    fn active_line(&self) -> &DelayedFilter {
        match self.image {
            ActiveImage::CharLine => &self.fir_lines[self.active],
            ActiveImage::VfitHold => self.vfit_lines[self.active]
                .as_ref().expect("hold image active implies realized"),
        }
    }

    fn active_direct(&self) -> f64 {
        match self.image {
            ActiveImage::CharLine => self.fir_taps[self.active][0],
            ActiveImage::VfitHold => 0.0, // explicitly delayed by construction
        }
    }

    /// Reflection for the NEXT junction solve, without advancing any state.
    ///
    /// The FIR's direct tap is not delayed or discarded. Fitted hold images
    /// retain their original physical transit time and residual filter. During
    /// a lift BOTH the direct term and history use the very same weight that
    /// `push` will commit. This is an affine wave boundary, not an impedance fit
    /// or a time-varying scattering-energy certificate.
    ///
    /// # Errors
    /// A previous partial step or an unrepresentable mixed boundary.
    pub fn boundary_response(&self) -> Result<WindBoundary, WindLineError> {
        self.check_live()?;
        let mut direct = self.active_direct();
        let mut history = self.active_line().incoming();
        if let Some(old) = self.fade_old.as_ref() {
            let w = 1.0 - (self.fade_remaining - 1) as f64 / self.fade_total as f64;
            direct = (1.0 - w) * self.fade_old_direct + w * direct;
            history = (1.0 - w) * old.incoming() + w * history;
        }
        if !direct.is_finite() || !history.is_finite() {
            return Err(WindLineError::Invalid { what: "wind boundary is not finite" });
        }
        Ok(WindBoundary { direct_reflection: direct, history_pressure_pa: history })
    }

    /// Commit an outgoing wave and return its SAME-sample incoming reflection.
    ///
    /// Junctions solve against `boundary_response` first; `incoming` alone is
    /// only a diagnostic of the previous step. No allocation occurs here.
    /// Invalid input leaves state unchanged. If a filter fails after another
    /// participant advanced, the bank is poisoned rather than silently resumed.
    ///
    /// # Errors
    /// A nonfinite input, reflected wave, filter output, or partial prior step.
    pub fn push(&mut self, outgoing: f64) -> Result<f64, WindLineError> {
        let boundary = self.boundary_response()?;
        let incoming = boundary.direct_reflection * outgoing + boundary.history_pressure_pa;
        if !outgoing.is_finite() || !incoming.is_finite() {
            return Err(WindLineError::Invalid { what: "wind wave must be finite" });
        }
        let advanced = match self.image {
            ActiveImage::CharLine => self.fir_lines[self.active].push(outgoing),
            ActiveImage::VfitHold => self.vfit_lines[self.active]
                .as_mut().expect("hold image active implies realized").push(outgoing),
        };
        if let Err(error) = advanced {
            self.poisoned = true;
            return Err(error.into());
        }
        if let Some(old) = self.fade_old.as_mut() {
            if let Err(error) = old.push(outgoing) {
                self.poisoned = true;
                return Err(error.into());
            }
            self.fade_remaining -= 1;
            if self.fade_remaining == 0 { self.fade_old = None; }
        }
        // Retain the longest COMPLETE FIR horizon, including during hold.
        // Pop before push so a full ring never transiently reallocates.
        if self.history.len() == self.history_capacity { self.history.pop_front(); }
        self.history.push_back(outgoing);
        self.last_incoming = incoming;
        Ok(incoming)
    }

    /// Re-read the last COMPLETED incoming wave without advancing.
    /// Use `boundary_response` for the next physical junction solve.
    #[must_use]
    pub fn incoming(&self) -> f64 { self.last_incoming }

    fn begin_fade(&mut self, old: DelayedFilter, direct: f64, samples: usize) {
        self.fade_old = if samples > 0 { Some(old) } else { None };
        self.fade_old_direct = direct;
        self.fade_remaining = samples;
        self.fade_total = samples;
    }

    /// Switch fingering on the char image (carry-history + crossfade —
    /// the MM bank's lift, plane-wave edition). Switching while on the
    /// hold image first hops back to char (a fingering change is a
    /// phrase event).
    ///
    /// # Errors
    /// [`WindLineError`] on an unknown fingering.
    pub fn switch_fingering(
        &mut self,
        index: usize,
        fade_samples: usize,
    ) -> Result<(), WindLineError> {
        self.check_live()?;
        if index >= self.fir_lines.len() {
            return Err(WindLineError::Invalid { what: "unknown fingering index" });
        }
        let old = self.active_line().clone();
        let direct = self.active_direct();
        // The active char line is stale while a hold image plays. The shared
        // outgoing history is the physical source for EVERY fingering lift.
        let history: Vec<f64> = self.history.iter().copied().collect();
        let lifted = history_line(self.dt, &self.fir_taps[index], &history)?;
        self.lift_log.push(WindLiftRecord {
            action: "switch-fingering/carry+fade", carried: history.len(), fade_samples,
        });
        self.begin_fade(old, direct, fade_samples);
        self.fir_lines[index] = lifted;
        self.active = index;
        self.image = ActiveImage::CharLine;
        Ok(())
    }

    /// Hop the ACTIVE fingering to the VFIT hold image: realize (once)
    /// the rational fit of the same fingering's reflectance, REPLAY the
    /// recent outgoing history into it (outputs discarded — the D17
    /// replay-prime lift), then crossfade. REFUSES when the hold
    /// image's gate is not green (scope 4: structurally impossible).
    ///
    /// # Errors
    /// [`WindLineError::ImageNotGated`] / realization refusals.
    pub fn hop_to_hold(
        &mut self,
        duct: &Duct,
        gas: &GasState,
        termination: Termination,
        fade_samples: usize,
    ) -> Result<(), WindLineError> {
        self.check_live()?;
        if !self.vfit_gate_green {
            return Err(WindLineError::ImageNotGated { image: "vfit-hold" });
        }
        if self.image == ActiveImage::VfitHold { return Ok(()); }
        if self.vfit_templates[self.active].is_none() {
            self.vfit_templates[self.active] =
                Some(realize_hold_line(duct, gas, termination, self.dt, self.zc)?);
        }
        // A finite replay window is still an approximate IIR-state lift, not
        // the infinite prehistory. Never replay it on top of a previous visit.
        let mut hold = self.vfit_templates[self.active].as_ref().expect("realized").clone();
        let carried = self.history.len();
        for &x in &self.history { hold.push(x)?; }
        let old = self.active_line().clone();
        let direct = self.active_direct();
        self.lift_log.push(WindLiftRecord {
            action: "hop-to-hold/replay-prime+fade", carried, fade_samples,
        });
        self.begin_fade(old, direct, fade_samples);
        self.vfit_lines[self.active] = Some(hold);
        self.image = ActiveImage::VfitHold;
        Ok(())
    }

    /// Hop back to the char image on gesture resumption (replay-prime
    /// the FIR line's ring buffer from the recent history via the
    /// carry-history constructor).
    ///
    /// # Errors
    /// [`WindLineError`] refusals.
    pub fn hop_back_to_char(&mut self, fade_samples: usize) -> Result<(), WindLineError> {
        self.check_live()?;
        if self.image == ActiveImage::CharLine { return Ok(()); }
        let history: Vec<f64> = self.history.iter().copied().collect();
        let lifted = history_line(self.dt, &self.fir_taps[self.active], &history)?;
        let old = self.active_line().clone();
        let direct = self.active_direct();
        self.lift_log.push(WindLiftRecord {
            action: "hop-back-to-char/replay-prime+fade", carried: history.len(), fade_samples,
        });
        self.begin_fade(old, direct, fade_samples);
        self.fir_lines[self.active] = lifted;
        self.image = ActiveImage::CharLine;
        Ok(())
    }
}

/// Build the strictly causal remainder using the original filter owner.
/// Its last output is the history term needed by the NEXT boundary solve.
/// Priming the newest sample is necessary: the carry-history constructor
/// restores the ring, but deliberately leaves `incoming()` at zero.
pub(crate) fn history_line(dt: f64, taps: &[f64], history: &[f64])
    -> Result<DelayedFilter, WindLineError>
{
    if taps.len() < 4 || taps.iter().any(|x| !x.is_finite()) {
        return Err(WindLineError::Invalid { what: "wind reflection needs at least four finite taps" });
    }
    let mut tail = taps[1..].to_vec();
    tail.resize(tail.len().max(4), 0.0);
    if let Some((&newest, earlier)) = history.split_last() {
        let mut line = DelayedFilter::from_impulse_response_with_history(dt, tail, earlier)?;
        line.push(newest)?;
        Ok(line)
    } else {
        Ok(DelayedFilter::from_impulse_response(dt, tail)?)
    }
}

/// Probe a zero-state line's FIR taps by impulse; no caller state changes.
fn probe_taps(line: &DelayedFilter) -> Result<Vec<f64>, WindLineError> {
    let mut probe = line.clone();
    (0..4096usize).map(|k| probe.push(if k == 0 { 1.0 } else { 0.0 })
        .map_err(WindLineError::Discrete)).collect()
}

/// Realize the VFIT hold image for one fingering: sweep the TMM
/// reflectance, conjugate onto the fs-vfit axis (the pinned trap), fit,
/// and realize with the geometric round-trip delay peeled into the
/// filter's own delay slot.
fn realize_hold_line(
    duct: &Duct,
    gas: &GasState,
    termination: Termination,
    dt: f64,
    zc: f64,
) -> Result<DelayedFilter, WindLineError> {
    let sweep = impedance_sweep(
        duct,
        gas,
        core::f64::consts::TAU * 120.0,
        core::f64::consts::TAU * 3000.0,
        160,
        LossModel::Bessel,
        termination,
    )
    .map_err(|_| WindLineError::Invalid {
        what: "hold-image sweep refused",
    })?;
    let omega: Vec<f64> = sweep.iter().map(|r| r.omega).collect();
    // Acoustic e^{-iwt} -> fs-vfit s = +iw: CONJUGATE (the ~50%-error
    // trap, pinned in the clarinet casebook).
    let h: Vec<C64> = sweep
        .iter()
        .map(|r| reflectance(r.impedance, zc).conj())
        .collect();
    let length: f64 =
        duct.segments
            .iter()
            .map(fs_duct::Segment::axial_length)
            .sum();
    let delay_samples = 2.0 * length / gas.sound_speed / dt;
    let mut opts = FitOptions::new(12);
    opts.iterations = 12;
    opts.fit_e = false;
    let prewarp = core::f64::consts::TAU * 600.0;
    // Preserve the ORIGINAL delay peeling and residual fit. Only the delay
    // buffer loses one sample: the bank's pending-history boundary contributes
    // that sample back. Passing D-1 to from_tabulated would instead change the
    // peeled fit and silently introduce the extra delay again.
    if !delay_samples.is_finite() || delay_samples < 3.0 {
        return Err(WindLineError::Invalid {
            what: "hold image needs at least three samples of physical round-trip transit",
        });
    }
    let peeled = modulate_delay(&omega, &h, delay_samples * dt, 1.0)?;
    let filter = realize_tabulated(&omega, &peeled, dt, &opts, prewarp)
        .map_err(|_| WindLineError::Invalid { what: "hold-image rational realization refused" })?;
    let line = DelayedFilter::new(delay_samples - 1.0, filter)?;
    let grid: Vec<f64> = (1..=64)
        .map(|k| core::f64::consts::TAU * 3000.0 * f64::from(k) / 64.0)
        .collect();
    // DEF-03: admit the hold image only when its fit is already
    // scattering-passive on this grid; a global 1/peak gain would rescale the
    // whole reflection. Measured on the wind_line tests' fingerings
    // (2026-10-19): the former clip never fired (gain 1.0).
    if !(line.scattering_peak(&grid) <= 1.0) {
        return Err(WindLineError::Invalid {
            what: "VFIT hold image is not scattering-passive (|R| > 1 on its grid); no gain clip",
        });
    }
    Ok(line)
}

#[cfg(test)]
mod wind_articulation_tests {
    pub(crate) mod helpers {
        pub(crate) use super::{Loop, air20, fingerings};
    }
    use super::*;
    use crate::reed_bore::{blowing_envelope, solve_reed_wave};
    use fs_duct::{Fingering, HoleState, Segment};
    use fs_material::gas::{GasSpec, GasState};
    use fs_scenario::BeatingReed;

    const RATE: u32 = 24_000;

    pub(crate) fn air20() -> GasState {
        GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("air")
    }

    fn verdict(case: &str, pass: bool, detail: &str) {
        println!(
            "{{\"suite\":\"fs-couple\",\"case\":\"{case}\",\"verdict\":\"{}\",\"detail\":\"{detail}\"}}",
            if pass { "pass" } else { "fail" }
        );
        assert!(pass, "case {case}: {detail}");
    }

    /// Clarinet-class bore with one register vent at `vent_frac` of the
    /// closed length and one tone hole near the foot.
    fn wind_template(vent_frac: f64) -> Duct {
        let bore = 0.0075f64;
        let l = 0.45f64;
        Duct {
            segments: vec![
                Segment::Cylinder {
                    radius: bore,
                    length: l * vent_frac,
                },
                Segment::ToneHole {
                    hole_radius: 0.0018,
                    chimney_height: 0.006,
                    bore_radius: bore,
                    state: HoleState::Closed,
                },
                Segment::Cylinder {
                    radius: bore,
                    length: l * (0.82 - vent_frac),
                },
                Segment::ToneHole {
                    hole_radius: 0.004,
                    chimney_height: 0.003,
                    bore_radius: bore,
                    state: HoleState::Closed,
                },
                Segment::Cylinder {
                    radius: bore,
                    length: l * 0.18,
                },
            ],
        }
    }

    pub(crate) fn fingerings(vent_frac: f64) -> FingeringTable {
        FingeringTable::try_new(
            wind_template(vent_frac),
            vec![
                Fingering {
                    label: "low".to_string(),
                    holes: vec![HoleState::Closed, HoleState::Closed],
                },
                Fingering {
                    label: "vented".to_string(),
                    holes: vec![HoleState::Open, HoleState::Closed],
                },
                Fingering {
                    label: "short".to_string(),
                    holes: vec![HoleState::Closed, HoleState::Open],
                },
            ],
        )
        .expect("table")
    }

    fn reed() -> BeatingReed {
        BeatingReed {
            rest_opening_m: 4.0e-4,
            width_m: 0.012,
            closing_pressure_pa: 2600.0,
            blowing_pressure_pa: 2000.0,
            attack_s: 0.01,
            mass_kg: 0.0,
            stiffness_n_m: 0.0,
            damping_ratio: 0.35,
        }
    }

    /// G3: before/after hold-sweep peaks. For this uniform cylinder the
    /// old interpolation selected WideTube over the WHOLE declared band,
    /// so that surviving benchmark arm is the exact old-law comparator.
    /// First uncaptured native measurement: ts1, aeb0c36ee98b, release,
    /// 2026-10-16: max |shift|=0.013889075618217473 cents. The declared
    /// regression bound is five times that measured shift. Receipt retained
    /// in brightowl-pexfp76/hold-first-native-measurement.log (bead .76).
    #[test]
    fn hold_image_bessel_shift_recorded() {
        use fs_duct::{Segment, impedance_peaks, segment_wave};
        const FIRST_NATIVE_MAX_CENTS: f64 = 1.3889075618217473e-2;
        let tolerance = 5.0 * FIRST_NATIVE_MAX_CENTS;
        let start = std::time::Instant::now();
        let gas = air20();
        let duct = Duct {
            segments: vec![Segment::Cylinder {
                radius: 0.0075,
                length: 0.33,
            }],
        };
        assert!(
            segment_wave(&gas, 0.0075, core::f64::consts::TAU * 120.0, LossModel::WideTube)
                .unwrap()
                .shear_number >= 10.0
        );
        let sweep = |loss| {
            impedance_sweep(
                &duct, &gas,
                core::f64::consts::TAU * 120.0, core::f64::consts::TAU * 3000.0,
                32000, loss, Termination::FlangedOpen,
            ).unwrap()
        };
        let old = sweep(LossModel::WideTube);
        let new = sweep(LossModel::Bessel);
        let peaks = |s: &[fs_duct::DuctResponse]| -> Vec<f64> {
            impedance_peaks(s).into_iter().take(4).map(|i| {
                let left = fs_math::det::ln(s[i - 1].impedance.abs());
                let mid = fs_math::det::ln(s[i].impedance.abs());
                let right = fs_math::det::ln(s[i + 1].impedance.abs());
                let offset = 0.5 * (left - right) / (left - 2.0 * mid + right);
                s[i].omega + offset * (s[i + 1].omega - s[i].omega)
            }).collect()
        };
        let old_peaks = peaks(&old);
        let new_peaks = peaks(&new);
        assert_eq!(old_peaks.len(), 4);
        assert_eq!(new_peaks.len(), 4);
        let mut max_abs_cents = 0.0_f64;
        for (index, (&before, &after)) in old_peaks.iter().zip(&new_peaks).enumerate() {
            let cents = 1200.0 * fs_math::det::ln(after / before) / core::f64::consts::LN_2;
            assert!(cents.is_finite());
            max_abs_cents = max_abs_cents.max(cents.abs());
            println!("{{\"case\":\"hold_image_bessel_shift_recorded\",\"inputs\":\"authored uniform-cylinder hold fixture, old rv>=10 branch=WideTube\",\"rung\":{index},\"axis\":\"peak\",\"observed\":{cents:e},\"old_hz\":{},\"new_hz\":{},\"oracle\":\"historical wide branch\",\"tol\":{tolerance:e},\"verdict\":\"pass\",\"elapsed_ms\":{}}}", before/core::f64::consts::TAU, after/core::f64::consts::TAU, start.elapsed().as_millis());
        }
        assert!(max_abs_cents > 0.0 && max_abs_cents < tolerance,
            "hold sweep shifted {max_abs_cents:e} cents; measured band {tolerance:e}");
        println!("{{\"case\":\"hold-image-max-absolute-cents\",\"inputs\":\"authored uniform-cylinder hold fixture\",\"observed\":{max_abs_cents:e},\"oracle\":\"historical wide branch\",\"tol\":{tolerance:e},\"verdict\":\"pass\",\"elapsed_ms\":{}}}", start.elapsed().as_millis());
    }

    pub(crate) struct Loop {
        pub(crate) bank: WindLineBank,
        zc: f64,
        rho: f64,
        p_plus_prev: f64,
        t: f64,
    }

    impl Loop {
        pub(crate) fn new(table: &FingeringTable, gas: &GasState, vfit_green: bool) -> Loop {
            let bore = 0.0075f64;
            let zc = gas.characteristic_impedance / (core::f64::consts::PI * bore * bore);
            let bank =
                WindLineBank::new(table, gas, Termination::FlangedOpen, RATE, zc, vfit_green)
                    .expect("bank");
            Loop {
                bank,
                zc,
                rho: gas.density,
                p_plus_prev: 5.0,
                t: 0.0,
            }
        }

        pub(crate) fn block(&mut self, out: &mut [f64]) {
            self.block_boundary_counterfactual(out, false);
        }

        // Test-only counterfactual for the pre-40aa49e66 junction timing.
        // Ordinary fixtures always take the correct same-sample branch.
        fn block_boundary_counterfactual(&mut self, out: &mut [f64], lagged: bool) {
            let dt = 1.0 / f64::from(RATE);
            for slot in out.iter_mut() {
                let p_m = blowing_envelope(reed(), self.t);
                let (direct,history) = if lagged {
                    (0.0,self.bank.incoming())
                } else {
                    let boundary=self.bank.boundary_response().expect("wave boundary");
                    (boundary.direct_reflection,boundary.history_pressure_pa)
                };
                let p_plus = solve_reed_wave(
                    reed(),
                    self.rho,
                    self.zc,
                    direct,
                    history,
                    p_m,
                    self.p_plus_prev,
                    0.0,
                )
                .expect("reed solve");
                let p_minus = self.bank.push(p_plus).expect("push");
                self.p_plus_prev = p_plus;
                *slot = p_plus + if lagged {history} else {p_minus};
                self.t += dt;
            }
        }
    }

    /// Fundamental estimator by NORMALIZED AUTOCORRELATION: the
    /// smallest lag that is a local maximum within 15% of the best
    /// in-range peak. THREE estimator traps were executed en route
    /// (the tuning record): a global FFT peak returns the dominant
    /// HARMONIC (the .6.4 lesson, re-caught here by the TMM peak map —
    /// a lock at a |Z|-weak peak is a misread, not physics); even-factor
    /// HPS (f,2f,3f) is STRUCTURALLY BLIND on a closed-pipe spectrum
    /// (the |X(2f)| factor lands on the ABSENT even harmonic and zeroes
    /// the true score — it picked f0/2 where the spectrum probe shows
    /// no line); and ODD-factor products (f,3f,5f) cannot disambiguate
    /// f0 from 3 f0 because the odd series is SELF-SIMILAR under 3x.
    /// The period domain has neither ambiguity: the smallest strong
    /// autocorrelation lag IS the period.
    fn dominant_hz(block: &[f64]) -> f64 {
        let mean = block.iter().sum::<f64>() / block.len() as f64;
        let x: Vec<f64> = block.iter().map(|v| v - mean).collect();
        let e0: f64 = x.iter().map(|v| v * v).sum();
        if e0 <= 0.0 {
            return 0.0;
        }
        let lag_min = (f64::from(RATE) / 900.0) as usize;
        let lag_max = ((f64::from(RATE) / 60.0) as usize).min(x.len() / 2);
        let mut ac = vec![0.0f64; lag_max + 1];
        for (lag, slot) in ac.iter_mut().enumerate().take(lag_max + 1).skip(lag_min) {
            let mut acc = 0.0f64;
            for i in 0..x.len() - lag {
                acc += x[i] * x[i + lag];
            }
            *slot = acc / e0;
        }
        let best = ac[lag_min..=lag_max].iter().fold(0.0f64, |m, &v| m.max(v));
        if best <= 0.0 {
            return 0.0;
        }
        for lag in lag_min + 1..lag_max {
            if ac[lag] > 0.85 * best && ac[lag] >= ac[lag - 1] && ac[lag] >= ac[lag + 1] {
                // Parabolic refine on the autocorr peak.
                let (ya, yb, yc) = (ac[lag - 1], ac[lag], ac[lag + 1]);
                let den = ya - 2.0 * yb + yc;
                let shift = if den.abs() > 1e-15 {
                    0.5 * (ya - yc) / den
                } else {
                    0.0
                };
                return f64::from(RATE) / (lag as f64 + shift);
            }
        }
        0.0
    }

    fn block_rms(signal: &[f64], b: usize, len: usize) -> f64 {
        let seg = &signal[b * len..(b + 1) * len];
        let mean = seg.iter().sum::<f64>() / len as f64;
        (seg.iter().map(|x| (x - mean) * (x - mean)).sum::<f64>() / len as f64).sqrt()
    }

    #[test]
    fn wa_001_sigma_transient_and_slur() {
        // Sustained low note, then the fingering slur (low -> short) with
        // carry+fade; the envelope holds through the change and the lock
        // lands on the short fingering's pitch.
        let gas = air20();
        let table = fingerings(0.30);
        let mut lp = Loop::new(&table, &gas, false);
        let len = 600usize;
        let mut signal = Vec::new();
        let mut block = vec![0.0f64; len];
        for b in 0..48usize {
            if b == 24 {
                lp.bank.switch_fingering(2, 96).expect("slur");
            }
            lp.block(&mut block);
            signal.extend_from_slice(&block);
        }
        let f_low = dominant_hz(&signal[10 * len..18 * len]);
        let f_short = dominant_hz(&signal[40 * len..48 * len]);
        let pre_rms = block_rms(&signal, 22, len);
        let through_min = (24..28)
            .map(|b| block_rms(&signal, b, len))
            .fold(f64::INFINITY, f64::min);
        let up_cents = 1200.0 * (f_short / f_low).log2();
        // The short fingering opens the big foot hole: the pitch RISES.
        let pass =
            up_cents > 150.0 && up_cents < 600.0 && through_min > 0.4 * pre_rms && f_low > 150.0;
        verdict(
            "wa-001-sigma-slur",
            pass,
            &format!(
                "slur low {f_low:.1} -> short {f_short:.1} Hz (+{up_cents:.0} cents); \
                 envelope through the change {through_min:.0} vs pre {pre_rms:.0}"
            ),
        );
    }

    #[test]
    fn wa_002_register_vent_flips_the_twelfth() {
        // Opening the register vent mid-note flips the lock toward the
        // TWELFTH (3x, closed-pipe odd series) — emergent, and POSITION
        // MATTERS: the same vent at the wrong station fails to flip
        // (the falsifier).
        let gas = air20();
        let run = |vent_frac: f64| -> (f64, f64) {
            let table = fingerings(vent_frac);
            let mut lp = Loop::new(&table, &gas, false);
            let len = 600usize;
            let mut signal = Vec::new();
            let mut block = vec![0.0f64; len];
            for b in 0..56usize {
                if b == 24 {
                    lp.bank.switch_fingering(1, 96).expect("vent");
                }
                lp.block(&mut block);
                signal.extend_from_slice(&block);
            }
            let f_before = dominant_hz(&signal[12 * len..20 * len]);
            let f_after = dominant_hz(&signal[44 * len..52 * len]);
            (f_before, f_after)
        };
        // MEASURED at d54edeb33 (ts1, release) after the Bessel-Zc and
        // flanged-mouth fixes. A sweep of this 1.8 mm / 6 mm vent flips to
        // the twelfth (x3.07..3.31) for stations 0.12..0.24 L. From 0.30 L
        // on, the lock stays on the vented first peak (x1.61 at 0.30, x1.14
        // at 0.70). The TMM authority agrees: only near-reed vents pull the
        // vented first peak (~300 Hz) below the second. That is the classic
        // register-hole rule near the third mode's pressure node. The old
        // "inverted rule at 0.70 L" was measured on the unconjugated Zc.
        // Both landed locks are cross-checked against the TMM authority's
        // peaks of their own vented geometry: the lock follows the
        // impedance peaks, emergent.
        let (f0_good, f1_good) = run(0.18);
        let ratio_good = f1_good / f0_good;
        let (f0_bad, f1_bad) = run(0.70);
        let ratio_bad = f1_bad / f0_bad;
        let tmm_peak_near = |vent_frac: f64, f: f64| -> f64 {
            use fs_duct::{impedance_peaks, impedance_sweep};
            let table = fingerings(vent_frac);
            let duct = table.duct("vented").expect("duct");
            let sweep = impedance_sweep(
                &duct,
                &gas,
                core::f64::consts::TAU * 80.0,
                core::f64::consts::TAU * 1200.0,
                6000,
                LossModel::Bessel,
                Termination::FlangedOpen,
            )
            .expect("sweep");
            impedance_peaks(&sweep)
                .iter()
                .map(|&i| sweep[i].omega / core::f64::consts::TAU)
                .min_by(|a, b| (a - f).abs().partial_cmp(&(b - f).abs()).expect("finite"))
                .expect("peaks")
        };
        let cents = |a: f64, b: f64| (1200.0 * (b / a).log2()).abs();
        let good_peak = tmm_peak_near(0.18, f1_good);
        let bad_peak = tmm_peak_near(0.70, f1_bad);
        let good_on_peak = cents(good_peak, f1_good);
        let bad_on_peak = cents(bad_peak, f1_bad);
        let pass = ratio_good > 2.8
            && ratio_good < 3.6
            && ratio_bad < 2.0
            && good_on_peak < 60.0
            && bad_on_peak < 60.0;
        verdict(
            "wa-002-register-vent",
            pass,
            &format!(
                "vent at 0.18L: {f0_good:.1} -> {f1_good:.1} Hz (x{ratio_good:.2}, the \
                 twelfth-class flip; {good_on_peak:.0} cents from the TMM peak \
                 {good_peak:.1}); falsifier at 0.70L: {f0_bad:.1} -> {f1_bad:.1} \
                 (x{ratio_bad:.2}, stays low; {bad_on_peak:.0} cents from its own peak \
                 {bad_peak:.1}) — position matters AND the lock follows the authority"
            ),
        );
    }

    /// G3 diagnostic: isolate boundary timing while keeping the exact-ZK
    /// owner, current filter/history implementation and all physical data fixed.
    /// This is not an alternate acceptance image or a claimed restoration.
    #[test]
    fn wa_002_register_vent_sample_convention_counterfactual() {
        let start=std::time::Instant::now();
        let gas=air20();
        let table=fingerings(0.70);
        let mut same_sample_signal:Option<Vec<f64>>=None;
        for lagged in [false,true] {
            let mut lp=Loop::new(&table,&gas,false);
            let len=600_usize;
            let mut signal=Vec::new();
            let mut block=vec![0.0;len];
            for b in 0..56 {
                if b==24 { lp.bank.switch_fingering(1,96).expect("vent"); }
                lp.block_boundary_counterfactual(&mut block,lagged);
                assert!(block.iter().all(|v|v.is_finite()));
                signal.extend_from_slice(&block);
            }
            let before=dominant_hz(&signal[12*len..20*len]);
            let after=dominant_hz(&signal[44*len..52*len]);
            if before<=0.0 || after<=0.0 {
                println!("{{\"case\":\"wa002-boundary-timing-counterfactual\",\"inputs\":\"unchanged vent .70L,1.8mm/6mm chimney, exact ZK and24kHz current FIR/history\",\"lagged\":{lagged},\"observed\":[{before},{after}],\"oracle\":\"same-sample versus previous-sample boundary\",\"tol\":\"diagnostic only; acceptance ratio2.8..3.6 unchanged\",\"verdict\":\"refused\",\"elapsed_us\":{}}}",start.elapsed().as_micros());
            } else {
                let ratio=after/before;
                println!("{{\"case\":\"wa002-boundary-timing-counterfactual\",\"inputs\":\"unchanged vent .70L,1.8mm/6mm chimney, exact ZK and24kHz current FIR/history\",\"lagged\":{lagged},\"observed\":[{before},{after}],\"ratio\":{ratio},\"oracle\":\"same-sample versus previous-sample boundary\",\"tol\":\"diagnostic only; acceptance ratio2.8..3.6 unchanged\",\"verdict\":\"pass\",\"elapsed_us\":{}}}",start.elapsed().as_micros());
            }
            if let Some(first)=&same_sample_signal {
                assert!(first.iter().zip(&signal).any(|(a,b)|a.to_bits()!=b.to_bits()),
                    "boundary-timing counterfactual is not live");
            } else { same_sample_signal=Some(signal); }
        }
    }

    #[test]
    #[ignore = "probe: raw spectrum of the settled unvented note"]
    fn zz_probe_spectrum() {
        use fs_fft::{C64 as FftC64, Fft};
        let gas = air20();
        let table = fingerings(0.30);
        let mut lp = Loop::new(&table, &gas, false);
        let len = 600usize;
        let mut signal = Vec::new();
        let mut block = vec![0.0f64; len];
        for _ in 0..24usize {
            lp.block(&mut block);
            signal.extend_from_slice(&block);
        }
        let seg = &signal[12 * len..12 * len + 4096];
        let n = 4096usize;
        let mean = seg.iter().sum::<f64>() / seg.len() as f64;
        let mut buf: Vec<FftC64> = (0..n).map(|k| FftC64::new(seg[k] - mean, 0.0)).collect();
        let mut scratch = vec![FftC64::new(0.0, 0.0); n];
        Fft::new(n).forward(&mut buf, &mut scratch);
        let df = f64::from(RATE) / n as f64;
        let mags: Vec<f64> = buf[..n / 2]
            .iter()
            .map(|c| (c.re * c.re + c.im * c.im).sqrt())
            .collect();
        let peak = mags.iter().fold(0.0f64, |m, &v| m.max(v));
        for k in ((40.0 / df) as usize)..((1000.0 / df) as usize) {
            if mags[k] > 0.05 * peak && mags[k] > mags[k - 1] && mags[k] > mags[k + 1] {
                println!("line at {:.1} Hz, rel {:.3}", k as f64 * df, mags[k] / peak);
            }
        }
    }

    #[test]
    #[ignore = "probe: TMM peak map for the vented geometries"]
    fn zz_probe_vent_peaks() {
        use fs_duct::{impedance_peaks, impedance_sweep};
        let gas = air20();
        for vent_frac in [0.15f64, 0.2, 0.25, 0.30, 0.4, 0.5, 0.6, 0.70] {
            let table = fingerings(vent_frac);
            let duct = table.duct("vented").expect("duct");
            let sweep = impedance_sweep(
                &duct,
                &gas,
                core::f64::consts::TAU * 80.0,
                core::f64::consts::TAU * 1200.0,
                6000,
                LossModel::Bessel,
                Termination::FlangedOpen,
            )
            .expect("sweep");
            let peaks = impedance_peaks(&sweep);
            let table_hz: Vec<f64> = peaks
                .iter()
                .take(4)
                .map(|&i| (sweep[i].omega / core::f64::consts::TAU * 10.0).round() / 10.0)
                .collect();
            let mags: Vec<f64> = peaks
                .iter()
                .take(4)
                .map(|&i| (sweep[i].impedance.abs() / 1.0e6).round())
                .collect();
            println!("vent at {vent_frac:.2}L: peaks {table_hz:?} Hz, |Z| {mags:?} MPa.s/m3");
        }
    }

    #[test]
    fn wa_003_gate_refusal_is_structural() {
        // A bank whose VFIT image is NOT gate-green cannot hop, by
        // construction.
        let gas = air20();
        let table = fingerings(0.30);
        let mut lp = Loop::new(&table, &gas, false);
        let duct = table.duct("low").expect("duct");
        let refused = lp
            .bank
            .hop_to_hold(&duct, &gas, Termination::FlangedOpen, 96);
        let pass = matches!(refused, Err(WindLineError::ImageNotGated { .. }));
        verdict(
            "wa-003-gate-refusal",
            pass,
            &format!("hop to an ungated image refused structurally: {pass}"),
        );
    }
}

#[cfg(test)]
mod wind_hop_tests {
    use super::wind_articulation_tests::helpers::*;
    use fs_duct::Termination;

    fn verdict(case: &str, pass: bool, detail: &str) {
        println!(
            "{{\"suite\":\"fs-couple\",\"case\":\"{case}\",\"verdict\":\"{}\",\"detail\":\"{detail}\"}}",
            if pass { "pass" } else { "fail" }
        );
        assert!(pass, "case {case}: {detail}");
    }

    /// Settle detector: relative windowed-RMS drift below `eps` for `m`
    /// consecutive blocks. THE PARAMETERS ARE THE POLICY DATA.
    const SETTLE_EPS: f64 = 0.05;
    const SETTLE_M: usize = 4;

    fn run_with_hop(hop_block: Option<usize>) -> (Vec<f64>, Vec<bool>, usize) {
        let gas = air20();
        let table = fingerings(0.30);
        // The VFIT hold image's registry gate IS green (the .6.1 wind
        // gates review) — the hop is legitimately enterable.
        let mut lp = Loop::new(&table, &gas, true);
        let duct = table.duct("low").expect("duct");
        let len = 600usize;
        let mut signal = Vec::new();
        let mut block = vec![0.0f64; len];
        let mut rms_prev = 0.0f64;
        let mut consec = 0usize;
        let mut settled = Vec::new();
        let mut hopped_at = usize::MAX;
        for b in 0..40usize {
            if Some(b) == hop_block {
                lp.bank
                    .hop_to_hold(&duct, &gas, Termination::FlangedOpen, 96)
                    .expect("hop");
                hopped_at = b;
            }
            lp.block(&mut block);
            signal.extend_from_slice(&block);
            let mean = block.iter().sum::<f64>() / len as f64;
            let rms =
                (block.iter().map(|x| (x - mean) * (x - mean)).sum::<f64>() / len as f64).sqrt();
            if rms_prev > 0.0 && ((rms - rms_prev) / rms_prev).abs() < SETTLE_EPS {
                consec += 1;
            } else {
                consec = 0;
            }
            settled.push(consec >= SETTLE_M);
            rms_prev = rms;
        }
        (signal, settled, hopped_at)
    }

    fn click_metric(signal: &[f64], hop_block: usize, len: usize) -> f64 {
        // Worst sample-step in the 3-block window from the hop,
        // normalized by the pre-hop RMS.
        let pre = &signal[(hop_block - 2) * len..hop_block * len];
        let mean = pre.iter().sum::<f64>() / pre.len() as f64;
        let pre_rms =
            (pre.iter().map(|x| (x - mean) * (x - mean)).sum::<f64>() / pre.len() as f64).sqrt();
        let seg = &signal[hop_block * len..(hop_block + 3) * len];
        seg.windows(2)
            .map(|w| (w[1] - w[0]).abs())
            .fold(0.0f64, f64::max)
            / pre_rms.max(1e-300)
    }

    #[test]
    fn wa_004_hop_policy_is_measured_data() {
        // The control run finds when the note settles; the hop is then
        // measured EARLY (during the attack — the falsifier showing the
        // click the policy exists to prevent) and SETTLED.
        let (_, settled, _) = run_with_hop(None);
        let first_settled = settled
            .iter()
            .position(|&s| s)
            .expect("the note must settle");
        let early_block = 2usize;
        assert!(
            !settled[early_block],
            "the early hop must be before the detector settles"
        );
        let settled_block = (first_settled + 2).max(early_block + 1);
        let len = 600usize;
        let (sig_early, _, _) = run_with_hop(Some(early_block));
        let (sig_settled, _, _) = run_with_hop(Some(settled_block));
        let click_early = click_metric(&sig_early, early_block, len);
        let click_settled = click_metric(&sig_settled, settled_block, len);
        // Sustain after the settled hop: the hold image keeps the note.
        let tail = &sig_settled[36 * len..40 * len];
        let tail_mean = tail.iter().sum::<f64>() / tail.len() as f64;
        let tail_rms = (tail
            .iter()
            .map(|x| (x - tail_mean) * (x - tail_mean))
            .sum::<f64>()
            / tail.len() as f64)
            .sqrt();
        // THE POLICY ARTIFACT (also committed as
        // data/claims/wind-hop-policy.tsv by the mint below).
        println!(
            "{{\"suite\":\"fs-couple\",\"case\":\"wa-004-hop-policy\",\
             \"settle_eps\":{SETTLE_EPS},\"settle_m\":{SETTLE_M},\
             \"first_settled_block\":{first_settled},\
             \"click_early\":{click_early:.3},\"click_settled\":{click_settled:.3},\
             \"tail_rms\":{tail_rms:.0}}}"
        );
        let pass = click_settled < 0.8 * click_early && tail_rms > 500.0;
        verdict(
            "wa-004-hop-policy",
            pass,
            &format!(
                "settled hop click {click_settled:.3} < 0.8x early-hop falsifier \
                 {click_early:.3}; note survives on the hold image (tail rms {tail_rms:.0})"
            ),
        );
    }

    #[test]
    #[ignore = "minting run: writes data/claims/wind-hop-policy.tsv"]
    fn mint_hop_policy_artifact() {
        let (_, settled, _) = run_with_hop(None);
        let first_settled = settled.iter().position(|&s| s).expect("settles");
        let len = 600usize;
        let mut rows = Vec::new();
        for hop_block in [2usize, 5, 8, 12, 16, 24, 32] {
            let (sig, _, _) = run_with_hop(Some(hop_block));
            let click = click_metric(&sig, hop_block, len);
            rows.push((hop_block, settled[hop_block.min(settled.len() - 1)], click));
        }
        let mut out = String::new();
        out.push_str("# frankensim-wind-hop-policy-v1\n");
        let _ = core::fmt::Write::write_fmt(
            &mut out,
            format_args!(
                "# settle detector: relative windowed-RMS drift < {SETTLE_EPS} for {SETTLE_M} \
             consecutive 25ms blocks; first settled block (control run): {first_settled}\n"
            ),
        );
        out.push_str("# minted by wind_hop_tests::mint_hop_policy_artifact (fs-couple)\n");
        out.push_str("hop_block\tsettled_at_hop\tclick_metric\n");
        for (b, s, c) in &rows {
            let _ = core::fmt::Write::write_fmt(&mut out, format_args!("{b}\t{s}\t{c:.4}\n"));
        }
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(std::path::Path::parent)
            .expect("root")
            .to_path_buf();
        std::fs::write(root.join("data/claims/wind-hop-policy.tsv"), out).expect("write");
        println!("minted wind-hop-policy.tsv: {rows:?}");
    }

    #[test]
    fn committed_hop_policy_is_consistent() {
        // The committed policy artifact must show the measured pattern.
        // MEASURED NUANCE: the hop one block BEFORE formal settle
        // (b=5, click 0.538) is already nearly as clean as settled —
        // the detector is conservatively LATE, which is the right side
        // to err on. The policy-relevant discrimination is therefore
        // against the TRUE ATTACK region (blocks <= 3, click 0.864):
        // the worst settled click must sit below the attack clicks.
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(std::path::Path::parent)
            .expect("root")
            .to_path_buf();
        let text = std::fs::read_to_string(root.join("data/claims/wind-hop-policy.tsv"))
            .expect("committed hop policy (mint test)");
        assert!(text.starts_with("# frankensim-wind-hop-policy-v1"));
        let mut settled_clicks = Vec::new();
        let mut attack_clicks = Vec::new();
        for line in text
            .lines()
            .filter(|l| !l.starts_with('#') && !l.starts_with("hop_block"))
        {
            let cols: Vec<&str> = line.split('\t').collect();
            let block: usize = cols[0].parse().expect("block");
            let settled: bool = cols[1].parse().expect("settled");
            let click: f64 = cols[2].parse().expect("click");
            if settled {
                settled_clicks.push(click);
            } else if block <= 3 {
                attack_clicks.push(click);
            }
        }
        assert!(!settled_clicks.is_empty() && !attack_clicks.is_empty());
        let worst_settled = settled_clicks.iter().fold(0.0f64, |m, &v| m.max(v));
        let best_attack = attack_clicks.iter().fold(f64::INFINITY, |m, &v| m.min(v));
        assert!(
            worst_settled < best_attack,
            "policy inversion: settled {worst_settled} vs attack {best_attack}"
        );
        println!(
            "{{\"suite\":\"fs-couple\",\"case\":\"wind-hop-policy-artifact\",\"verdict\":\
             \"pass\",\"worst_settled\":{worst_settled:.3},\"best_attack\":{best_attack:.3}}}"
        );
    }
}


#[cfg(test)]
mod same_sample_boundary_tests {
    use super::*;
    use fs_vfit::discretize::{Biquad, DigitalFilter};

    fn bank(taps: Vec<Vec<f64>>) -> WindLineBank {
        let dt = 1.0 / 48_000.0;
        let n = taps.len();
        let capacity = taps.iter().map(Vec::len).max().unwrap();
        WindLineBank {
            fir_lines: taps.iter().map(|h| history_line(dt, h, &[]).unwrap()).collect(),
            fir_taps: taps, vfit_lines: vec![None; n], vfit_templates: vec![None; n],
            labels: (0..n).map(|i| i.to_string()).collect(), active: 0,
            image: ActiveImage::CharLine, vfit_gate_green: true, dt, zc: 1.0e6,
            fade_old: None, fade_old_direct: 0.0, fade_remaining: 0, fade_total: 0,
            history: std::collections::VecDeque::with_capacity(capacity),
            history_capacity: capacity, last_incoming: 0.0, poisoned: false,
            lift_log: Vec::new(),
        }
    }
    fn filter(dt: f64) -> DigitalFilter {
        DigitalFilter { sections: vec![Biquad { b: [0.2, 0.1, 0.0], a: [-0.3, 0.0] }],
            direct: 0.1, t_s: dt, prewarp: 0.0 }
    }
    fn close(a: f64, b: f64) {
        assert!((a-b).abs() <= 2.0e-12 * a.abs().max(b.abs()).max(1.0), "{a} != {b}");
    }

    #[test]
    fn same_sample_boundary_matches_direct_convolution_and_ring_wrap() {
        let taps = vec![0.15, -0.2, 0.1, 0.05, -0.03, 0.02];
        let mut bank = bank(vec![taps.clone()]);
        let mut signal = Vec::new();
        for n in 0..200 {
            let input = f64::from(n % 17) - 8.0;
            signal.push(input);
            let expected: f64 = taps.iter().enumerate().take(signal.len())
                .map(|(k, h)| h * signal[signal.len()-1-k]).sum();
            let boundary = bank.boundary_response().unwrap();
            assert_eq!(boundary.direct_reflection, taps[0]);
            close(boundary.direct_reflection*input + boundary.history_pressure_pa, expected);
            close(bank.push(input).unwrap(), expected);
            close(bank.incoming(), expected);
        }
    }

    #[test]
    fn fingering_lift_keeps_more_than_512_samples_and_next_fade_weight() {
        let mut a = vec![0.0; 1024];
        let mut b = a.clone();
        a[0] = 0.1; a[700] = 0.25;
        b[0] = -0.2; b[700] = -0.3;
        let mut bank = bank(vec![a.clone(), b.clone()]);
        let mut history = Vec::new();
        for n in 0..800 {
            let x = f64::from(n % 19) - 9.0;
            history.push(x); bank.push(x).unwrap();
        }
        bank.switch_fingering(1, 4).unwrap();
        assert_eq!(bank.lift_log().last().unwrap().carried, 800);
        for step in 1..=6 {
            let x = f64::from(step);
            history.push(x);
            let w = (f64::from(step)/4.0).min(1.0);
            let old = a[0]*x + a[700]*history[history.len()-701];
            let new = b[0]*x + b[700]*history[history.len()-701];
            let boundary = bank.boundary_response().unwrap();
            close(boundary.direct_reflection, (1.0-w)*a[0]+w*b[0]);
            close(bank.push(x).unwrap(), (1.0-w)*old+w*new);
        }
    }

    #[test]
    fn pending_hold_preserves_integer_and_fractional_physical_transit() {
        let dt = 1.0/48_000.0;
        for delay in [3.0, 3.25, 7.75, 12.0] {
            let mut original = DelayedFilter::new(delay, filter(dt)).unwrap();
            let mut pending = DelayedFilter::new(delay-1.0, filter(dt)).unwrap();
            for n in 0..200 {
                let x = f64::from(n % 11)-5.0;
                let predicted = pending.incoming();
                close(original.push(x).unwrap(), predicted);
                pending.push(x).unwrap();
            }
        }
    }

    #[test]
    fn hold_reentry_replays_a_fresh_template_and_uses_actual_recent_waves() {
        let taps = vec![0.1, 0.2, 0.0, -0.1, 0.05, 0.0, 0.0, 0.0];
        let mut bank = bank(vec![taps.clone(), taps.clone()]);
        // The template is already realized; the supplied dummy duct is never
        // swept. Test the state lift, not TMM fitting or source calibration.
        let template = DelayedFilter::new(2.5, filter(bank.dt)).unwrap();
        bank.vfit_templates[0] = Some(template.clone());
        let gas = super::wind_articulation_tests::air20();
        let duct = Duct { segments: Vec::new() };
        for x in [1.0, 2.0, -1.0, 4.0] { bank.push(x).unwrap(); }
        bank.hop_to_hold(&duct, &gas, Termination::FlangedOpen, 4).unwrap();
        for x in [3.0, -2.0, 0.5] { bank.push(x).unwrap(); }
        bank.hop_back_to_char(0).unwrap();
        assert!(bank.fade_old.is_none());
        for x in [5.0, -3.0, 1.0] { bank.push(x).unwrap(); }
        let mut expected = template;
        for &x in &bank.history { expected.push(x).unwrap(); }
        bank.hop_to_hold(&duct, &gas, Termination::FlangedOpen, 0).unwrap();
        assert!(bank.fade_old.is_none());
        close(bank.boundary_response().unwrap().history_pressure_pa, expected.incoming());
        // Switch DIRECTLY from hold: the dormant char ring must not be used.
        for x in [13.0, -7.0, 11.0, 2.0] { bank.push(x).unwrap(); }
        let history: Vec<f64> = bank.history.iter().copied().collect();
        bank.switch_fingering(1, 0).unwrap();
        let x = 3.0;
        let expected = taps[0]*x + taps.iter().enumerate().skip(1)
            .map(|(k,h)| h*history[history.len()-k]).sum::<f64>();
        close(bank.push(x).unwrap(), expected);
    }

    #[test]
    fn nonfinite_input_does_not_mutate_wave_or_lift_history() {
        let mut bank = bank(vec![vec![0.1, 0.2, 0.1, 0.0]]);
        bank.push(2.0).unwrap();
        let before = bank.boundary_response().unwrap();
        let history = bank.history.clone();
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(bank.push(bad).is_err());
            assert_eq!(bank.boundary_response().unwrap(), before);
            assert_eq!(bank.history, history);
        }
        assert_eq!(bank.push(3.0).unwrap(), 0.1*3.0+0.2*2.0);
    }
}
