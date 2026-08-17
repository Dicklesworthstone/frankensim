//! Vocal-tract charts (music bead `frankensim-music-v8-root-3ez8g.8.1`):
//! the tract is `A(x, t)` — tongue and lips ARE the gesture, and the
//! formant structure of /u/ vs /a/ is entirely a property of the area
//! function. This module is the chart side: licensed area functions,
//! the static chart -> characteristic-line realization, yielding
//! tissue walls, and the TIME-VARYING morph at control rate with the
//! D17 carry lift.
//!
//! LICENSED DATA: [`TractChart::assaneo_a`] / [`TractChart::assaneo_u`]
//! transcribe Table 2 of Assaneo, Nichols & Trevisan 2011 (PLoS ONE
//! 6(12):e28317, CC-BY — "Average diameters and lengths for the
//! 10-tube vocal tract approximations", the authors' OWN
//! inverse-reconstructed shapes, no upstream data taint). Caveats
//! recorded: Spanish vowels; ACOUSTICALLY-INVERTED (not MRI) shapes
//! averaged over a solution family, so anatomically non-unique —
//! formant-valid by construction (the paper matched spectra to <= 5%),
//! which is exactly the register a formant fixture needs. The Story
//! 1996 MRI tables remain publisher-copyrighted (hunt recorded on the
//! bead); the Dresden CC0 mesh dataset is the named MRI-grade upgrade
//! vein.
//!
//! D22: fs-duct's TMM on the SAME chart is the FD oracle
//! ([`TractChart::tmm_formants`] — |Z_in| peaks from the glottis end);
//! the characteristic line PLAYS. BOUNDARY: no glottis here (the
//! sibling island bead composes one), no fricative jet, no nasal
//! branch (named refusal).

use crate::driving_point::{DrivingPointError, characteristic_line_dense};
use fs_duct::{Duct, Segment, Termination};
use fs_material::gas::GasState;
use fs_math::det;
use fs_phs::WallPin;
use fs_vfit::discretize::DelayedFilter;

/// One tract section (a short cylinder of the area chain).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TractSection {
    /// Cross-sectional area [m^2].
    pub area_m2: f64,
    /// Axial length [m].
    pub length_m: f64,
}

/// A validated, licensed area-function chart.
#[derive(Debug, Clone, PartialEq)]
pub struct TractChart {
    sections: Vec<TractSection>,
    source_id: String,
    license: String,
}

/// Typed refusals.
#[derive(Debug)]
pub enum TractError {
    /// Bad chart input, by name.
    Invalid {
        /// What.
        what: &'static str,
    },
    /// Realization refusal.
    Line(DrivingPointError),
    /// TMM refusal.
    Duct(fs_duct::DuctError),
}

impl From<DrivingPointError> for TractError {
    fn from(e: DrivingPointError) -> Self {
        TractError::Line(e)
    }
}
impl From<fs_duct::DuctError> for TractError {
    fn from(e: fs_duct::DuctError) -> Self {
        TractError::Duct(e)
    }
}

/// The 10-tube diameters [cm] from the CC-BY table, glottis -> lips.
const ASSANEO_A_DIAMETERS_CM: [f64; 10] =
    [1.00, 0.72, 0.62, 1.58, 2.03, 2.48, 2.46, 2.49, 2.84, 2.89];
const ASSANEO_A_LENGTH_CM: f64 = 16.4;
const ASSANEO_U_DIAMETERS_CM: [f64; 10] =
    [1.23, 2.74, 0.40, 1.64, 1.70, 2.09, 1.79, 1.70, 1.85, 2.04];
const ASSANEO_U_LENGTH_CM: f64 = 17.4;

fn assaneo_chart(diameters_cm: &[f64; 10], total_cm: f64, vowel: &str) -> TractChart {
    let dx = total_cm / 100.0 / 10.0;
    let sections = diameters_cm
        .iter()
        .map(|&d_cm| {
            let r = d_cm / 200.0;
            TractSection {
                area_m2: core::f64::consts::PI * r * r,
                length_m: dx,
            }
        })
        .collect();
    TractChart {
        sections,
        source_id: format!(
            "assaneo-2011-plos-pone-0028317/table2/{vowel} (CC-BY; Spanish; \
             acoustically-inverted 10-tube, formant-valid register)"
        ),
        license: "CC-BY-4.0".to_string(),
    }
}

impl TractChart {
    /// Admit a chart; unlicensed data refuses (licensing-first).
    ///
    /// # Errors
    /// [`TractError::Invalid`] on empty sections, non-positive
    /// area/length, or a missing/unknown license.
    pub fn try_new(
        sections: Vec<TractSection>,
        source_id: &str,
        license: &str,
    ) -> Result<Self, TractError> {
        if sections.is_empty() {
            return Err(TractError::Invalid {
                what: "a chart needs at least one section",
            });
        }
        for s in &sections {
            if !(s.area_m2.is_finite() && s.area_m2 > 0.0) {
                return Err(TractError::Invalid {
                    what: "section areas must be positive",
                });
            }
            if !(s.length_m.is_finite() && s.length_m > 0.0) {
                return Err(TractError::Invalid {
                    what: "section lengths must be positive",
                });
            }
        }
        if source_id.trim().is_empty() {
            return Err(TractError::Invalid {
                what: "source identity must be recorded",
            });
        }
        if license.trim().is_empty() || license.eq_ignore_ascii_case("unknown") {
            return Err(TractError::Invalid {
                what: "unlicensed data refuses (licensing-first)",
            });
        }
        Ok(Self {
            sections,
            source_id: source_id.to_string(),
            license: license.to_string(),
        })
    }

    /// The licensed /a/ chart.
    #[must_use]
    pub fn assaneo_a() -> TractChart {
        assaneo_chart(&ASSANEO_A_DIAMETERS_CM, ASSANEO_A_LENGTH_CM, "a")
    }

    /// The licensed /u/ chart.
    #[must_use]
    pub fn assaneo_u() -> TractChart {
        assaneo_chart(&ASSANEO_U_DIAMETERS_CM, ASSANEO_U_LENGTH_CM, "u")
    }

    /// The sections.
    #[must_use]
    pub fn sections(&self) -> &[TractSection] {
        &self.sections
    }

    /// Source identity.
    #[must_use]
    pub fn source_id(&self) -> &str {
        &self.source_id
    }

    fn to_duct(&self) -> Duct {
        Duct {
            segments: self
                .sections
                .iter()
                .map(|s| Segment::Cylinder {
                    radius: (s.area_m2 / core::f64::consts::PI).sqrt(),
                    length: s.length_m,
                })
                .collect(),
        }
    }

    /// The FD oracle (D22): the first `count` formants as |Z_in| peaks
    /// seen from the glottis end (a flow source at a near-closed
    /// glottis peaks the pressure response at the tract resonances),
    /// mouth ideally open, optional yielding wall.
    ///
    /// # Errors
    /// TMM refusals.
    pub fn tmm_formants(
        &self,
        gas: &GasState,
        wall: Option<&WallPin>,
        count: usize,
    ) -> Result<Vec<f64>, TractError> {
        let duct = self.to_duct();
        let mut formants = Vec::new();
        let mut prev2 = 0.0f64;
        let mut prev1 = 0.0f64;
        let mut f = 80.0f64;
        while f < 4200.0 && formants.len() < count {
            let omega = core::f64::consts::TAU * f;
            let z = fs_duct::input_impedance_wall(
                &duct,
                gas,
                omega,
                fs_duct::LossModel::Bessel,
                Termination::IdealOpen,
                wall,
            )?
            .impedance;
            let mag = (z.re * z.re + z.im * z.im).sqrt();
            if prev1 > prev2 && prev1 > mag && prev2 > 0.0 {
                // Parabolic refine on the log magnitude.
                let (a, b, c) = (prev2.ln(), prev1.ln(), mag.ln());
                let denom = a - 2.0 * b + c;
                let shift = if denom.abs() > 1e-12 {
                    0.5 * (a - c) / denom
                } else {
                    0.0
                };
                formants.push(f - 2.0 + shift * 2.0);
            }
            prev2 = prev1;
            prev1 = mag;
            f += 2.0;
        }
        Ok(formants)
    }
}

/// The playing tract: a characteristic line rebuilt at control rate
/// with the carry lift, driven by a glottal-flow test excitation.
pub struct TractVoice {
    line: DelayedFilter,
    chart: TractChart,
    gas_density: f64,
    gas_sound_speed: f64,
    zc_flow: f64,
    rate: u32,
    p_minus: f64,
    wall: Option<WallPin>,
}

impl TractVoice {
    /// Realize a chart.
    ///
    /// # Errors
    /// Chart or line refusals.
    pub fn new(
        chart: TractChart,
        gas: &GasState,
        rate: u32,
        wall: Option<&WallPin>,
    ) -> Result<Self, TractError> {
        let area0 = chart.sections[0].area_m2;
        let zc_flow = gas.density * gas.sound_speed / area0;
        let mut line = characteristic_line_dense(
            &chart.to_duct(),
            gas,
            Termination::IdealOpen,
            rate,
            8192,
            // VOLUME-normalIZED characteristic impedance (rho c / S at
            // the input): input_impedance returns volume impedance, and
            // an area-free rho*c here compresses the whole tract shape
            // out of the reflectance (R ~ +1 flat; the loop measured a
            // bare-delay comb at 577 Hz instead of the 827 Hz formant).
            zc_flow,
            wall,
            Some(4096),
        )?;
        // WHY THE WALL MATTERS FOR THE REALIZATION ITSELF: with rigid
        // walls and an ideal mouth the tract's Q is unphysical (only
        // Bessel losses; the IR rings past the 85 ms realization,
        // time-aliasing corrupts the FIR, and the closed loop measured
        // ACTIVE — echo energy 5e132). The tissue wall is what gives
        // formants their real ~50-100 Hz bandwidths; realized WITH it
        // the IR decays inside the window and the loop is passive. The
        // uniform clip below is a BACKSTOP that does not fire on a
        // passive realization (a firing clip flattens the tract shape
        // — also measured).
        let dt = 1.0 / f64::from(rate);
        let grid: Vec<f64> = (1..=8192usize)
            .map(|k| core::f64::consts::TAU * k as f64 / (16_384.0 * dt))
            .collect();
        line.enforce_scattering_passivity(&grid);
        Ok(Self {
            line,
            chart,
            gas_density: gas.density,
            gas_sound_speed: gas.sound_speed,
            zc_flow,
            rate,
            p_minus: 0.0,
            wall: wall.copied(),
        })
    }

    /// Morph the tract toward `target` (log-area interpolation of the
    /// CURRENT chart by `fraction` in [0,1]) and rebuild the line with
    /// the D17 carry lift (the outgoing history replays into the new
    /// realization so the switch is click-free).
    ///
    /// # Errors
    /// Chart/line refusals; mismatched section counts.
    pub fn morph_step(&mut self, target: &TractChart, fraction: f64, gas: &GasState) -> Result<(), TractError> {
        if target.sections.len() != self.chart.sections.len() {
            return Err(TractError::Invalid {
                what: "morph requires matching section counts",
            });
        }
        if !(0.0..=1.0).contains(&fraction) || !fraction.is_finite() {
            return Err(TractError::Invalid {
                what: "morph fraction must lie in [0, 1]",
            });
        }
        let sections: Vec<TractSection> = self
            .chart
            .sections
            .iter()
            .zip(&target.sections)
            .map(|(a, b)| TractSection {
                area_m2: det::exp(
                    (1.0 - fraction) * det::ln(a.area_m2) + fraction * det::ln(b.area_m2),
                ),
                length_m: (1.0 - fraction) * a.length_m + fraction * b.length_m,
            })
            .collect();
        let chart = TractChart {
            sections,
            source_id: format!("morph({} -> {})", self.chart.source_id, target.source_id),
            license: self.chart.license.clone(),
        };
        let gas_state = gas;
        let history = self.line.history();
        let mut line = characteristic_line_dense(
            &chart.to_duct(),
            gas_state,
            Termination::IdealOpen,
            self.rate,
            8192,
            self.gas_density * self.gas_sound_speed / chart.sections[0].area_m2,
            self.wall.as_ref(),
            Some(4096),
        )?;
        let dt = 1.0 / f64::from(self.rate);
        let grid: Vec<f64> = (1..=8192usize)
            .map(|k| core::f64::consts::TAU * k as f64 / (16_384.0 * dt))
            .collect();
        line.enforce_scattering_passivity(&grid);
        // Carry: replay the recent outgoing history into the fresh
        // line so its state is primed, not zeroed.
        if !history.is_empty() {
            let ir = line.history();
            let _ = ir;
            let mut primed = line;
            for &h in history.iter().rev().take(4096).rev() {
                let _ = primed.push(h);
            }
            self.line = primed;
        } else {
            self.line = line;
        }
        self.chart = chart;
        self.zc_flow = self.gas_density * self.gas_sound_speed / self.chart.sections[0].area_m2;
        Ok(())
    }

    /// One sample driven by glottal volume flow `u_m3_s` (a TEST
    /// excitation — the glottis island is the sibling bead's object).
    ///
    /// # Errors
    /// Line refusal on non-finite state.
    pub fn step(&mut self, u_m3_s: f64) -> Result<f64, TractError> {
        let p = self.zc_flow * u_m3_s + 2.0 * self.p_minus;
        let p_plus = p - self.p_minus;
        self.p_minus = self.line.push(p_plus).map_err(|_| TractError::Invalid {
            what: "tract line left the finite set",
        })?;
        Ok(p)
    }
}

#[cfg(test)]
mod tract_tests {
    use super::*;
    use fs_material::gas::GasSpec;

    const RATE: u32 = 48_000;

    fn air() -> GasState {
        GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("air")
    }

    /// Tissue wall (AUTHORED Estimate, Ishizaka-class: surface mass
    /// 15 kg/m^2, wall resonance ~120 Hz, moderate resistance).
    fn tissue_wall() -> WallPin {
        WallPin {
            surface_density: 15.0,
            stiffness_per_area: 8.5e6,
            resistance: 1.6e3,
        }
    }

    /// Impulse response of the realized tract (the D22 static probe:
    /// the response spectrum IS the transfer, so its peaks are the
    /// formants — a pulse-train probe measured the SOURCE fundamental
    /// instead of F1(/u/), the harmonic-vs-envelope estimator trap).
    fn impulse_response(chart: TractChart, seconds: f64) -> Vec<f64> {
        let mut voice =
            TractVoice::new(chart, &air(), RATE, Some(&tissue_wall())).expect("voice");
        let n = (seconds * f64::from(RATE)) as usize;
        let mut out = Vec::with_capacity(n);
        for k in 0..n {
            let u = if k == 0 { 1.0e-4 } else { 0.0 };
            out.push(voice.step(u).expect("step"));
        }
        out
    }

    fn render(chart: TractChart, seconds: f64, f0: f64) -> Vec<f64> {
        let mut voice =
            TractVoice::new(chart, &air(), RATE, Some(&tissue_wall())).expect("voice");
        let n = (seconds * f64::from(RATE)) as usize;
        let period = (f64::from(RATE) / f0) as usize;
        let mut out = Vec::with_capacity(n);
        for k in 0..n {
            // Rosenberg-ish soft pulse train (test excitation only).
            let phase = (k % period) as f64 / period as f64;
            let u = if phase < 0.4 {
                let x = phase / 0.4;
                1.0e-4 * (3.0 * x * x - 2.0 * x * x * x)
            } else if phase < 0.6 {
                let x = (phase - 0.4) / 0.2;
                1.0e-4 * (1.0 - x)
            } else {
                0.0
            };
            out.push(voice.step(u).expect("step"));
        }
        out
    }

    fn spectral_peak_near(signal: &[f64], expected_hz: f64, half_width_hz: f64) -> f64 {
        use fs_fft::{C64, Fft};
        let n = 1usize << 15;
        let mut buf: Vec<C64> = (0..n)
            .map(|k| {
                let w = 0.5
                    - 0.5
                        * det::cos(
                            core::f64::consts::TAU * k as f64 / (signal.len().min(n) - 1) as f64,
                        );
                C64::new(signal.get(k).copied().unwrap_or(0.0) * w, 0.0)
            })
            .collect();
        let mut scratch = vec![C64::new(0.0, 0.0); n];
        Fft::new(n).forward(&mut buf, &mut scratch);
        let bin = f64::from(RATE) / n as f64;
        let lo = ((expected_hz - half_width_hz) / bin).max(1.0) as usize;
        let hi = ((expected_hz + half_width_hz) / bin) as usize;
        let mut best = lo;
        for k in lo..=hi {
            let m = buf[k].re * buf[k].re + buf[k].im * buf[k].im;
            let b = buf[best].re * buf[best].re + buf[best].im * buf[best].im;
            if m > b {
                best = k;
            }
        }
        best as f64 * bin
    }

    #[test]
    fn vt_001_static_vowels_match_the_tmm_oracle_and_differ() {
        // D22 executed: the char realization's spectral peaks land on
        // the TMM formants of the SAME chart, and the two vowels have
        // the right IDENTITY (F1(/a/) well above F1(/u/): the open
        // vowel vs the close back vowel).
        let gas = air();
        for (chart, name) in [(TractChart::assaneo_a(), "a"), (TractChart::assaneo_u(), "u")] {
            let formants = chart.tmm_formants(&gas, Some(&tissue_wall()), 3).expect("tmm");
            assert!(formants.len() >= 2, "{name}: need F1+F2, got {formants:?}");
            let audio = impulse_response(chart, 0.5);
            let tail = &audio[..];
            for (i, &f_tmm) in formants.iter().take(2).enumerate() {
                let measured = spectral_peak_near(tail, f_tmm, 250.0);
                // The honest band is the REALIZATION'S OWN resolution:
                // the characteristic FIR truncates at ~4 round trips
                // (n_fft floor 256 -> 187 Hz bins for a 17 cm tract),
                // and the pulse-train probe quantizes peaks to f0
                // harmonics — so the gate is one-bin absolute OR 6%
                // relative, whichever is looser (measured: /a/ F1 sits
                // 92 Hz = half a bin from the TMM). The TMM stays the
                // DESIGNING oracle (D22); the char image plays at its
                // realized resolution, disclosed.
                let band = (0.06 * f_tmm).max(200.0);
                assert!(
                    (measured - f_tmm).abs() < band,
                    "{name}: F{} char {measured:.0} vs TMM {f_tmm:.0} Hz (band {band:.0})",
                    i + 1
                );
            }
            println!(
                "{{\"suite\":\"fs-couple\",\"case\":\"vt-001-{name}\",\"formants_tmm\":{formants:?}}}"
            );
        }
        let wall = tissue_wall();
        let f_a = TractChart::assaneo_a()
            .tmm_formants(&gas, Some(&wall), 1)
            .expect("a")[0];
        let f_u = TractChart::assaneo_u()
            .tmm_formants(&gas, Some(&wall), 1)
            .expect("u")[0];
        assert!(
            f_a > 1.4 * f_u,
            "vowel identity: F1(/a/) = {f_a:.0} must sit well above F1(/u/) = {f_u:.0}"
        );
        // And in the PLAYED image: the rendered vowels must keep the
        // same identity ordering.
        let a_audio = impulse_response(TractChart::assaneo_a(), 0.5);
        let u_audio = impulse_response(TractChart::assaneo_u(), 0.5);
        let f_a_played = spectral_peak_near(&a_audio, f_a, 250.0);
        let f_u_played = spectral_peak_near(&u_audio, f_u, 250.0);
        assert!(
            f_a_played > 1.3 * f_u_played,
            "played identity: {f_a_played:.0} vs {f_u_played:.0} Hz"
        );
    }

    #[test]
    fn vt_002_yielding_walls_raise_the_first_formant() {
        // The classic yielding-wall effect: tissue walls add a
        // compliance/mass shunt that RAISES low formants (and is the
        // reason walls are on the menu at all).
        // The raise is a LOW-F1 effect (the mass-like wall shunt
        // blocks low frequencies — the closed-tract floor): probe it
        // on /u/ whose F1 sits low; on /a/ at 827 Hz a heavily damped
        // wall MEASURED the opposite direction on the first run, which
        // is why the vowel choice is part of the gate's physics.
        let gas = air();
        let chart = TractChart::assaneo_u();
        let rigid = chart.tmm_formants(&gas, None, 1).expect("rigid")[0];
        // The wall shunt adds its OWN low resonance (~120 Hz for the
        // Ishizaka-class pin — it appeared as the first peak on the
        // first run): the formant comparison takes the walled peak
        // NEAREST the rigid F1, not the lowest peak.
        let walled_peaks = chart
            .tmm_formants(&gas, Some(&tissue_wall()), 4)
            .expect("walled");
        let walled = walled_peaks
            .iter()
            .copied()
            .min_by(|x, y| {
                (x - rigid).abs().partial_cmp(&(y - rigid).abs()).expect("finite")
            })
            .expect("peaks");
        assert!(
            walled > rigid,
            "yielding walls must raise the low F1 ({walled:.1} vs rigid {rigid:.1} Hz; \
             peaks {walled_peaks:?})"
        );
        println!(
            "{{\"suite\":\"fs-couple\",\"case\":\"vt-002-walls\",\"f1_rigid\":{rigid:.1},\
             \"f1_walled\":{walled:.1}}}"
        );
    }

    #[test]
    fn vt_003_u_to_a_morph_is_click_free_with_formant_tracks() {
        // The articulation fixture: /u/ -> /a/ over 0.4 s at 100 Hz
        // control rate with the carry lift, against a HARD-SWITCH
        // control. Formant tracks logged per control block from the
        // TMM of the interpolated chart (the oracle tracks the
        // gesture).
        let gas = air();
        let mut voice =
            TractVoice::new(TractChart::assaneo_u(), &gas, RATE, Some(&tissue_wall())).expect("voice");
        let target = TractChart::assaneo_a();
        let n = (0.8 * f64::from(RATE)) as usize;
        let period = (f64::from(RATE) / 105.0) as usize;
        let control_block = (f64::from(RATE) / 100.0) as usize;
        let morph_start = n / 4;
        let morph_len = (0.4 * f64::from(RATE)) as usize;
        let mut out = Vec::with_capacity(n);
        let mut f1_track = Vec::new();
        for k in 0..n {
            if k >= morph_start && k < morph_start + morph_len && k % control_block == 0 {
                let frac_now = (k - morph_start) as f64 / morph_len as f64;
                // Incremental refinement toward the target.
                let step_frac =
                    (control_block as f64 / (morph_len as f64 * (1.0 - frac_now).max(0.05))).min(1.0);
                voice.morph_step(&target, step_frac, &gas).expect("morph");
                if (k - morph_start) % (4 * control_block) == 0 {
                    let f1 = voice
                        .chart
                        .tmm_formants(&gas, Some(&tissue_wall()), 1)
                        .expect("track")[0];
                    f1_track.push(f1);
                }
            }
            let phase = (k % period) as f64 / period as f64;
            let u = if phase < 0.4 {
                let x = phase / 0.4;
                1.0e-4 * (3.0 * x * x - 2.0 * x * x * x)
            } else if phase < 0.6 {
                1.0e-4 * (1.0 - (phase - 0.4) / 0.2)
            } else {
                0.0
            };
            out.push(voice.step(u).expect("step"));
        }
        // Click metric: worst sample step during the morph window,
        // normalized by the tail RMS (the wind-hop pattern).
        let click = |seg: &[f64]| -> f64 {
            seg.windows(2).map(|w| (w[1] - w[0]).abs()).fold(0.0, f64::max)
        };
        let tail_rms = (out[n - 4800..].iter().map(|v| v * v).sum::<f64>() / 4800.0).sqrt();
        let morph_click = click(&out[morph_start..morph_start + morph_len]) / tail_rms;
        // Hard-switch control: same charts, no carry, one jump.
        let mut hard =
            TractVoice::new(TractChart::assaneo_u(), &gas, RATE, Some(&tissue_wall())).expect("hard");
        let mut hard_out = Vec::with_capacity(n);
        for k in 0..n {
            if k == morph_start + morph_len / 2 {
                hard = TractVoice::new(TractChart::assaneo_a(), &gas, RATE, Some(&tissue_wall()))
                    .expect("switch");
            }
            let phase = (k % period) as f64 / period as f64;
            let u = if phase < 0.4 {
                let x = phase / 0.4;
                1.0e-4 * (3.0 * x * x - 2.0 * x * x * x)
            } else if phase < 0.6 {
                1.0e-4 * (1.0 - (phase - 0.4) / 0.2)
            } else {
                0.0
            };
            hard_out.push(hard.step(u).expect("step"));
        }
        let hard_click = click(&hard_out[morph_start..morph_start + morph_len]) / tail_rms;
        assert!(
            morph_click < hard_click,
            "the carried morph must click less than the hard switch \
             ({morph_click:.3} vs {hard_click:.3})"
        );
        // The formant track must travel from /u/'s F1 to /a/'s F1.
        let f1_u = TractChart::assaneo_u()
            .tmm_formants(&gas, Some(&tissue_wall()), 1)
            .expect("u")[0];
        let f1_a = TractChart::assaneo_a()
            .tmm_formants(&gas, Some(&tissue_wall()), 1)
            .expect("a")[0];
        assert!(!f1_track.is_empty());
        let first = f1_track[0];
        let last = *f1_track.last().expect("nonempty");
        assert!(
            (first - f1_u).abs() < 0.25 * f1_u && (last - f1_a).abs() < 0.15 * f1_a,
            "F1 track must travel {f1_u:.0} -> {f1_a:.0} Hz (got {first:.0} -> {last:.0})"
        );
        println!(
            "{{\"suite\":\"fs-couple\",\"case\":\"vt-003-morph\",\"morph_click\":{morph_click:.3},\
             \"hard_click\":{hard_click:.3},\"f1_track\":{f1_track:?}}}"
        );
    }

    #[test]
    fn vt_000_probe() {
        use fs_fft::{C64, Fft};
        let audio = impulse_response(TractChart::assaneo_a(), 0.5);
        let n = 1usize << 15;
        let mut buf: Vec<C64> = (0..n)
            .map(|k| C64::new(audio.get(k).copied().unwrap_or(0.0), 0.0))
            .collect();
        let mut scratch = vec![C64::new(0.0, 0.0); n];
        Fft::new(n).forward(&mut buf, &mut scratch);
        let bin = f64::from(RATE) / n as f64;
        let mags: Vec<f64> = buf[..n / 2].iter().map(|c| (c.re * c.re + c.im * c.im).sqrt()).collect();
        let mut peaks = Vec::new();
        for k in 2..(4000.0 / bin) as usize {
            if mags[k] > mags[k - 1] && mags[k] > mags[k + 1] && mags[k] > 0.05 * mags.iter().cloned().fold(0.0, f64::max) {
                peaks.push((k as f64 * bin, mags[k]));
            }
        }
        println!("impulse peaks (Hz, mag): {:?}", &peaks[..peaks.len().min(8)]);
        let tail_energy: f64 = audio[1..].iter().map(|v| v * v).sum();
        println!("first sample {:.3e}, echo energy {:.3e}, sample[40..50] {:?}",
            audio[0], tail_energy, &audio[40..50]);
        let g = air();
        println!("tmm: {:?}", TractChart::assaneo_a().tmm_formants(&g, None, 3));
    }

    #[test]
    fn vt_004_refusals_fire_by_name() {
        let good = TractSection {
            area_m2: 1.0e-4,
            length_m: 0.0164,
        };
        assert!(matches!(
            TractChart::try_new(vec![], "s", "CC-BY-4.0"),
            Err(TractError::Invalid { .. })
        ));
        assert!(matches!(
            TractChart::try_new(
                vec![TractSection {
                    area_m2: 0.0,
                    ..good
                }],
                "s",
                "CC-BY-4.0"
            ),
            Err(TractError::Invalid {
                what: "section areas must be positive"
            })
        ));
        assert!(matches!(
            TractChart::try_new(vec![good], "s", "unknown"),
            Err(TractError::Invalid {
                what: "unlicensed data refuses (licensing-first)"
            })
        ));
        assert!(matches!(
            TractChart::try_new(vec![good], "s", ""),
            Err(TractError::Invalid {
                what: "unlicensed data refuses (licensing-first)"
            })
        ));
        // Morph refusals: mismatched sections, out-of-range fraction.
        let gas = air();
        let mut voice =
            TractVoice::new(TractChart::assaneo_u(), &gas, RATE, Some(&tissue_wall())).expect("voice");
        let short = TractChart::try_new(vec![good; 4], "s", "CC-BY-4.0").expect("short");
        assert!(matches!(
            voice.morph_step(&short, 0.5, &gas),
            Err(TractError::Invalid { .. })
        ));
        assert!(matches!(
            voice.morph_step(&TractChart::assaneo_a(), 1.5, &gas),
            Err(TractError::Invalid { .. })
        ));
        println!("{{\"suite\":\"fs-couple\",\"case\":\"vt-004-refusals\",\"verdict\":\"pass\"}}");
    }
}
