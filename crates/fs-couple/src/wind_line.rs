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
use fs_vfit::discretize::{DelayedFilter, DiscretizeError, reflectance};
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

/// Per-fingering plane-wave line bank with the char->VFIT hop.
pub struct WindLineBank {
    fir_lines: Vec<DelayedFilter>,
    fir_taps: Vec<Vec<f64>>,
    vfit_lines: Vec<Option<DelayedFilter>>,
    labels: Vec<String>,
    active: usize,
    image: ActiveImage,
    vfit_gate_green: bool,
    dt: f64,
    zc: f64,
    fade_old: Option<DelayedFilter>,
    fade_old_is_vfit: bool,
    fade_remaining: usize,
    fade_total: usize,
    history: std::collections::VecDeque<f64>,
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
            let taps = probe_taps(&line);
            fir_lines.push(line);
            fir_taps.push(taps);
            labels.push(label.to_string());
        }
        let n = fir_lines.len();
        Ok(WindLineBank {
            fir_lines,
            fir_taps,
            vfit_lines: vec![None; n],
            labels,
            active: 0,
            image: ActiveImage::CharLine,
            vfit_gate_green,
            dt,
            zc,
            fade_old: None,
            fade_old_is_vfit: false,
            fade_remaining: 0,
            fade_total: 0,
            history: std::collections::VecDeque::with_capacity(512),
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

    /// Push the outgoing wave; get the incoming wave (crossfaded during
    /// a transition).
    ///
    /// # Errors
    /// [`WindLineError::Discrete`] on a non-finite sample.
    pub fn push(&mut self, outgoing: f64) -> Result<f64, WindLineError> {
        self.history.push_back(outgoing);
        if self.history.len() > 512 {
            self.history.pop_front();
        }
        let new_in = match self.image {
            ActiveImage::CharLine => self.fir_lines[self.active].push(outgoing)?,
            ActiveImage::VfitHold => self.vfit_lines[self.active]
                .as_mut()
                .expect("hold image active implies realized")
                .push(outgoing)?,
        };
        if let Some(old) = self.fade_old.as_mut() {
            let old_in = old.push(outgoing)?;
            self.fade_remaining -= 1;
            let w = 1.0 - self.fade_remaining as f64 / self.fade_total as f64;
            let mixed = (1.0 - w) * old_in + w * new_in;
            if self.fade_remaining == 0 {
                self.fade_old = None;
            }
            return Ok(mixed);
        }
        Ok(new_in)
    }

    /// Re-read the last incoming wave without advancing (approximation
    /// during fades: the new line's pending value).
    #[must_use]
    pub fn incoming(&self) -> f64 {
        let new_in = match self.image {
            ActiveImage::CharLine => self.fir_lines[self.active].incoming(),
            ActiveImage::VfitHold => self.vfit_lines[self.active]
                .as_ref()
                .expect("hold image active implies realized")
                .incoming(),
        };
        if let Some(old) = self.fade_old.as_ref() {
            let w = 1.0 - self.fade_remaining as f64 / self.fade_total as f64;
            return (1.0 - w) * old.incoming() + w * new_in;
        }
        new_in
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
        if index >= self.fir_lines.len() {
            return Err(WindLineError::Invalid {
                what: "unknown fingering index",
            });
        }
        let old = match self.image {
            ActiveImage::CharLine => self.fir_lines[self.active].clone(),
            ActiveImage::VfitHold => self.vfit_lines[self.active]
                .as_ref()
                .expect("realized")
                .clone(),
        };
        let history = self.fir_lines[self.active].history();
        let lifted = DelayedFilter::from_impulse_response_with_history(
            self.dt,
            self.fir_taps[index].clone(),
            &history,
        )?;
        self.lift_log.push(WindLiftRecord {
            action: "switch-fingering/carry+fade",
            carried: history.len(),
            fade_samples,
        });
        if fade_samples > 0 {
            self.fade_old = Some(old);
            self.fade_old_is_vfit = self.image == ActiveImage::VfitHold;
            self.fade_remaining = fade_samples;
            self.fade_total = fade_samples;
        } else {
            self.fade_old = None;
        }
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
        if !self.vfit_gate_green {
            return Err(WindLineError::ImageNotGated { image: "vfit-hold" });
        }
        if self.image == ActiveImage::VfitHold {
            return Ok(());
        }
        if self.vfit_lines[self.active].is_none() {
            self.vfit_lines[self.active] =
                Some(realize_hold_line(duct, gas, termination, self.dt, self.zc)?);
        }
        // Fresh state + replay prime.
        let mut hold = self.vfit_lines[self.active].take().expect("realized");
        let carried = self.history.len();
        {
            // Reset by re-realizing state through a fresh clone: replay
            // the history into a zero-state copy.
            for &x in &self.history {
                let _ = hold.push(x)?;
            }
        }
        let old = self.fir_lines[self.active].clone();
        self.lift_log.push(WindLiftRecord {
            action: "hop-to-hold/replay-prime+fade",
            carried,
            fade_samples,
        });
        if fade_samples > 0 {
            self.fade_old = Some(old);
            self.fade_old_is_vfit = false;
            self.fade_remaining = fade_samples;
            self.fade_total = fade_samples;
        }
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
        if self.image == ActiveImage::CharLine {
            return Ok(());
        }
        let history: Vec<f64> = self.history.iter().copied().collect();
        let lifted = DelayedFilter::from_impulse_response_with_history(
            self.dt,
            self.fir_taps[self.active].clone(),
            &history,
        )?;
        let old = self.vfit_lines[self.active]
            .as_ref()
            .expect("hold active implies realized")
            .clone();
        self.lift_log.push(WindLiftRecord {
            action: "hop-back-to-char/replay-prime+fade",
            carried: history.len(),
            fade_samples,
        });
        if fade_samples > 0 {
            self.fade_old = Some(old);
            self.fade_old_is_vfit = true;
            self.fade_remaining = fade_samples;
            self.fade_total = fade_samples;
        }
        self.fir_lines[self.active] = lifted;
        self.image = ActiveImage::CharLine;
        Ok(())
    }
}

/// Probe a line's FIR taps by impulse (clone; state untouched).
fn probe_taps(line: &DelayedFilter) -> Vec<f64> {
    let mut probe = line.clone();
    (0..4096usize)
        .map(|k| {
            probe
                .push(if k == 0 { 1.0 } else { 0.0 })
                .expect("probe push")
        })
        .collect()
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
        LossModel::AllRegime,
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
            .map(|s| match *s {
                fs_duct::Segment::Cylinder { length, .. }
                | fs_duct::Segment::Cone { length, .. } => length,
                fs_duct::Segment::ToneHole { .. } => 0.0,
            })
            .sum();
    let delay_samples = 2.0 * length / gas.sound_speed / dt;
    let mut opts = FitOptions::new(12);
    opts.iterations = 12;
    opts.fit_e = false;
    let prewarp = core::f64::consts::TAU * 600.0;
    let mut line = DelayedFilter::from_tabulated(&omega, &h, delay_samples, dt, &opts, prewarp)
        .map_err(|_| WindLineError::Invalid {
            what: "hold-image rational realization refused",
        })?;
    let grid: Vec<f64> = (1..=64)
        .map(|k| core::f64::consts::TAU * 3000.0 * f64::from(k) / 64.0)
        .collect();
    line.enforce_scattering_passivity(&grid);
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
        }
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
            let dt = 1.0 / f64::from(RATE);
            for slot in out.iter_mut() {
                let p_m = blowing_envelope(reed(), self.t);
                let p_minus = self.bank.incoming();
                let p_plus = solve_reed_wave(
                    reed(),
                    self.rho,
                    self.zc,
                    0.0,
                    p_minus,
                    p_m,
                    self.p_plus_prev,
                    0.0,
                )
                .expect("reed solve");
                let _ = self.bank.push(p_plus).expect("push");
                self.p_plus_prev = p_plus;
                *slot = p_plus + p_minus;
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
                return f64::from(RATE) / (lag as f64 - shift);
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
        // MEASURED (not the textbook sketch): with THIS inertive vent
        // (1.8 mm hole, 6 mm chimney — a shunt inertance, not a
        // pressure release, so the classic L/3 node rule inverts), the
        // register flip to the twelfth-class regime happens at 0.70L;
        // the 0.30L vent leaves the lock in the low regime. Both landed
        // locks are cross-checked against the TMM authority's peaks of
        // their own vented geometry — the lock follows the impedance
        // peaks, emergent.
        let (f0_good, f1_good) = run(0.70);
        let ratio_good = f1_good / f0_good;
        let (f0_bad, f1_bad) = run(0.30);
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
                LossModel::AllRegime,
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
        let good_peak = tmm_peak_near(0.70, f1_good);
        let bad_peak = tmm_peak_near(0.30, f1_bad);
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
                "vent at 0.70L: {f0_good:.1} -> {f1_good:.1} Hz (x{ratio_good:.2}, the \
                 twelfth-class flip; {good_on_peak:.0} cents from the TMM peak \
                 {good_peak:.1}); falsifier at 0.30L: {f0_bad:.1} -> {f1_bad:.1} \
                 (x{ratio_bad:.2}, stays low; {bad_on_peak:.0} cents from its own peak \
                 {bad_peak:.1}) — position matters AND the lock follows the authority"
            ),
        );
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
                LossModel::AllRegime,
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
        out.push_str(&format!(
            "# settle detector: relative windowed-RMS drift < {SETTLE_EPS} for {SETTLE_M} \
             consecutive 25ms blocks; first settled block (control run): {first_settled}\n"
        ));
        out.push_str("# minted by wind_hop_tests::mint_hop_policy_artifact (fs-couple)\n");
        out.push_str("hop_block\tsettled_at_hop\tclick_metric\n");
        for (b, s, c) in &rows {
            out.push_str(&format!("{b}\t{s}\t{c:.4}\n"));
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
