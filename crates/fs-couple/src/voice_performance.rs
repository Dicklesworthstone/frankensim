//! A glottal source driving a tract area function, observed at a microphone.
//!
//! The file declares physical inputs only: the lumped fold card, an optional
//! tract wall, the mouth termination, the area function as ordered sections,
//! the gas, the subglottal pressure the player applies, and where the
//! microphone stands. Nothing here names a vowel, a register or a singer.
//!
//! Composition, all existing owners: [`GlottalIsland`] (fold pHS, Bernoulli
//! aperture, fold collision) against the tract's characteristic line; the
//! mouth volume flow of the SAME transfer-matrix duct
//! ([`crate::driving_point::outlet`], driven by the island's incident wave);
//! and the baffled circular-opening receiver, whose pressure is the retarded
//! derivative of that flow. No source waveform, formant filter or output gain
//! exists, and the supraglottal pressure is never the output.
//!
//! Limits, not hidden: infinite-baffle radiation (no head diffraction), a
//! plane-wave series tract (no side branches or nasal coupling), a lumped
//! two-mass or one-mass fold, and a rigid or single-pin yielding wall.
//!
//! ```text
//! frankensim-voice-performance-v1
//! # comment lines and blank lines are ignored
//! audio <output_rate_hz> <samples> <mechanics_per_output>
//! full_scale_pa <pa>
//! ambient moist_air <temperature_k> <pressure_pa> <relative_humidity>
//! ambient mixture <temperature_k> <pressure_pa> <mole|mass> <count> <species fraction>...
//! ambient thermochemical <temperature_k> <pressure_pa> <mole|mass> <count> <id fraction>...
//! fold <two_mass|one_dof> <m_lower_kg> <m_upper_kg> <k_lower_n_m> <k_upper_n_m>
//!      <k_couple_n_m> <damping_ratio> <rest_gap_m> <fold_length_m> <thickness_m>
//!      [<thickness_upper_m>]      (the second viscous section; default 0.05 cm)
//! tract_wall <surface_density_kg_m2> <stiffness_n_m3> <resistance_pa_s_m>   (optional)
//! termination <flanged|unflanged|ideal_open>
//! tract_source <identifier without spaces>
//! tract_licence <free text to end of line>
//! section <area_m2> <length_m>                       (repeat, glottis first)
//! microphone <x_m> <y_m> <z_m> <radial_rings> <angular_points> <max_frequency_hz>
//! subglottal <time_s> <pressure_pa>                  (repeat, strictly increasing time)
//! ```

use crate::acoustic_realize::AcousticRealizeError;
use crate::driving_point::LineLadder;
use crate::driving_point::outlet::mouth_flow_line;
use crate::glottis::{FoldCard, FoldCardSource, GlottalIsland};
use crate::pcm_wav::baffled::{BaffledPressure, CircularOutletReceiver, RayleighMedium};
use crate::pcm_wav::observation::PressureRenderer;
use crate::render::RenderError;
use crate::tract::{TractChart, TractSection};
use fs_duct::Termination;
use fs_material::gas::GasState;
use fs_material::gas_species::input::{GasRecord, ThermochemicalSources, read_gas_record};
use fs_phs::WallPin;
use fs_vfit::discretize::DelayedFilter;

/// First line of every voice performance file.
pub const VOICE_PERFORMANCE_SCHEMA: &str = "frankensim-voice-performance-v1";
/// Largest accepted file.
pub const MAX_VOICE_PERFORMANCE_BYTES: usize = 1 << 20;
const MAX_SECTIONS: usize = 1024;
const MAX_SCHEDULE_POINTS: usize = 4096;
const MAX_SECONDS: u64 = 120;
const MAX_BLOCK: usize = 65_536;

/// A parsed, validated, not yet realized performance.
#[derive(Debug, Clone)]
pub struct VoiceDescription {
    /// Output (microphone) clock [Hz].
    pub output_rate_hz: u32,
    /// Output samples to render.
    pub samples: u64,
    /// Integer physical substeps per output sample, 1..=16.
    pub mechanics_per_output: u32,
    /// Pascals mapped to PCM full scale by the command output.
    pub full_scale_pa: f64,
    /// One immutable gas record fills the tract and the exterior.
    gas: GasRecord,
    /// `true` for the two-mass fold, `false` for the one-mass image.
    pub two_mass: bool,
    /// Lumped fold parameters.
    pub fold: FoldCard,
    /// Optional single-pin yielding tract wall.
    pub wall: Option<WallPin>,
    /// Mouth radiation load.
    pub termination: Termination,
    /// Area function, glottis first.
    pub sections: Vec<TractSection>,
    /// Identifier of the area-function source.
    pub tract_source: String,
    /// Licence or terms the data were admitted under.
    pub tract_licence: String,
    /// Microphone in the mouth-plane frame (origin at the mouth centre,
    /// +z out of the baffle) and its spatial resolution.
    pub microphone: CircularOutletReceiver,
    /// Subglottal pressure breakpoints `(time_s, pressure_pa)`.
    pub subglottal: Vec<(f64, f64)>,
}

fn field(rest: &[&str], i: usize, what: &str) -> Result<f64, String> {
    let token = rest.get(i).ok_or_else(|| format!("missing {what}"))?;
    let value: f64 = token
        .parse()
        .map_err(|_| format!("{what} is not a number: {token:?}"))?;
    if value.is_finite() {
        Ok(value)
    } else {
        Err(format!("{what} must be finite"))
    }
}

fn whole(rest: &[&str], i: usize, what: &str) -> Result<u64, String> {
    let token = rest.get(i).ok_or_else(|| format!("missing {what}"))?;
    token
        .parse::<u64>()
        .map_err(|_| format!("{what} is not a non-negative integer: {token:?}"))
}

fn arity(rest: &[&str], n: usize, key: &str) -> Result<(), String> {
    if rest.len() == n {
        Ok(())
    } else {
        Err(format!("`{key}` takes {n} fields, got {}", rest.len()))
    }
}

fn once<T>(slot: &mut Option<T>, value: T, key: &str) -> Result<(), String> {
    if slot.replace(value).is_some() {
        Err(format!("`{key}` may only appear once"))
    } else {
        Ok(())
    }
}

impl VoiceDescription {
    /// Parse and validate a performance file.
    ///
    /// # Errors
    /// Oversize or non-UTF-8 input, an unknown or repeated key, a malformed or
    /// non-finite field, or any non-physical value, each named.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, String> {
        Self::from_bytes_with_gas_sources(bytes, None)
    }

    /// Admit legacy air, named mixtures or source-backed NASA-9 gas records.
    /// All tract, fold-flow and radiation consumers borrow the same state;
    /// source receipts are retained without changing its coefficients.
    /// The medium is stationary; this does not model species transport.
    ///
    /// # Errors
    /// Ordinary file refusals plus missing, ambiguous or out-of-domain gas
    /// sources. NASA-9 records never fall back to fixed heat capacities.
    pub fn from_bytes_with_gas_sources(
        bytes: &[u8],
        sources: Option<&ThermochemicalSources<'_>>,
    ) -> Result<Self, String> {
        if bytes.len() > MAX_VOICE_PERFORMANCE_BYTES {
            return Err("voice performance file exceeds its size limit".into());
        }
        let text = std::str::from_utf8(bytes).map_err(|_| "voice performance must be UTF-8")?;
        let mut lines = text
            .lines()
            .map(str::trim)
            .filter(|l| !l.is_empty() && !l.starts_with('#'));
        if lines.next() != Some(VOICE_PERFORMANCE_SCHEMA) {
            return Err(format!("first line must be {VOICE_PERFORMANCE_SCHEMA}"));
        }
        let mut audio = None;
        let mut full_scale = None;
        let mut gas = None;
        let mut fold = None;
        let mut wall = None;
        let mut termination = None;
        let mut source = None;
        let mut licence = None;
        let mut microphone = None;
        let mut sections: Vec<TractSection> = Vec::new();
        let mut subglottal: Vec<(f64, f64)> = Vec::new();
        for line in lines {
            let mut tokens = line.split_whitespace();
            let key = tokens.next().unwrap_or_default();
            let rest: Vec<&str> = tokens.collect();
            match key {
                "audio" => {
                    arity(&rest, 3, key)?;
                    let rate = whole(&rest, 0, "audio rate")?;
                    let samples = whole(&rest, 1, "audio samples")?;
                    let ratio = whole(&rest, 2, "mechanics per output")?;
                    if !(8_000..=192_000).contains(&rate) || !(1..=16).contains(&ratio) {
                        return Err("audio needs a rate in 8000..=192000 and a 1..=16 ratio".into());
                    }
                    if samples == 0 || samples > MAX_SECONDS * rate {
                        return Err(format!("audio samples must be in 1..={MAX_SECONDS} s"));
                    }
                    once(&mut audio, (rate as u32, samples, ratio as u32), key)?;
                }
                "full_scale_pa" => {
                    arity(&rest, 1, key)?;
                    let pa = field(&rest, 0, "full scale")?;
                    if pa <= 0.0 {
                        return Err("full_scale_pa must be positive".into());
                    }
                    once(&mut full_scale, pa, key)?;
                }
                "ambient" => {
                    let record_text = rest.join(" ");
                    let mut fields = record_text.split_ascii_whitespace();
                    if fields.clone().next() == Some("moist_air") {
                        fields.next();
                    }
                    let record = read_gas_record(&mut fields, None, sources)
                        .map_err(|e| format!("ambient gas refused: {e}"))?;
                    if fields.next().is_some() {
                        return Err("ambient gas contains extra fields".into());
                    }
                    once(&mut gas, record, key)?;
                }
                "fold" => {
                    // Ten fields; an optional eleventh sets the upper-section thickness of the
                    // viscous channel (default: Ishizaka-Flanagan's 0.05 cm, `FoldCard`).
                    if rest.len() != 11 {
                        arity(&rest, 10, key)?;
                    }
                    let two_mass = match rest[0] {
                        "two_mass" => true,
                        "one_dof" => false,
                        _ => return Err("fold model must be `two_mass` or `one_dof`".into()),
                    };
                    let card = FoldCard {
                        mass_lower_kg: field(&rest, 1, "lower mass")?,
                        mass_upper_kg: field(&rest, 2, "upper mass")?,
                        stiffness_lower_n_m: field(&rest, 3, "lower stiffness")?,
                        stiffness_upper_n_m: field(&rest, 4, "upper stiffness")?,
                        coupling_n_m: field(&rest, 5, "coupling stiffness")?,
                        damping_ratio: field(&rest, 6, "damping ratio")?,
                        rest_gap_m: field(&rest, 7, "rest gap")?,
                        fold_length_m: field(&rest, 8, "fold length")?,
                        thickness_m: field(&rest, 9, "thickness")?,
                        thickness_upper_m: if rest.len() == 11 {
                            field(&rest, 10, "upper thickness")?
                        } else {
                            FoldCard::two_mass_standard().thickness_upper_m
                        },
                        // Per-field provenance lives in the file's comments.
                        source: FoldCardSource::AuthoredEstimate,
                    };
                    once(&mut fold, (two_mass, card), key)?;
                }
                "tract_wall" => {
                    arity(&rest, 3, key)?;
                    let pin = WallPin {
                        surface_density: field(&rest, 0, "wall surface density")?,
                        stiffness_per_area: field(&rest, 1, "wall stiffness")?,
                        resistance: field(&rest, 2, "wall resistance")?,
                    };
                    once(&mut wall, pin, key)?;
                }
                "termination" => {
                    arity(&rest, 1, key)?;
                    let load = match rest[0] {
                        "flanged" => Termination::FlangedOpen,
                        "unflanged" => Termination::UnflangedOpen,
                        "ideal_open" => Termination::IdealOpen,
                        _ => return Err("termination must be flanged, unflanged or ideal_open".into()),
                    };
                    once(&mut termination, load, key)?;
                }
                "tract_source" => {
                    arity(&rest, 1, key)?;
                    once(&mut source, rest[0].to_string(), key)?;
                }
                "tract_licence" => {
                    let text = rest.join(" ");
                    if text.is_empty() {
                        return Err("tract_licence needs text".into());
                    }
                    once(&mut licence, text, key)?;
                }
                "section" => {
                    arity(&rest, 2, key)?;
                    if sections.len() >= MAX_SECTIONS {
                        return Err("too many tract sections".into());
                    }
                    sections.push(TractSection {
                        area_m2: field(&rest, 0, "section area")?,
                        length_m: field(&rest, 1, "section length")?,
                    });
                }
                "microphone" => {
                    arity(&rest, 6, key)?;
                    let receiver = CircularOutletReceiver {
                        position_m: [
                            field(&rest, 0, "microphone x")?,
                            field(&rest, 1, "microphone y")?,
                            field(&rest, 2, "microphone z")?,
                        ],
                        radial_rings: whole(&rest, 3, "radial rings")? as usize,
                        angular_points: whole(&rest, 4, "angular points")? as usize,
                        maximum_frequency_hz: field(&rest, 5, "resolved band")?,
                    };
                    if receiver.position_m[2] <= 0.0 {
                        return Err("the microphone must stand in front of the baffle (z > 0)".into());
                    }
                    once(&mut microphone, receiver, key)?;
                }
                "subglottal" => {
                    arity(&rest, 2, key)?;
                    if subglottal.len() >= MAX_SCHEDULE_POINTS {
                        return Err("too many subglottal pressure points".into());
                    }
                    let point = (field(&rest, 0, "time")?, field(&rest, 1, "pressure")?);
                    if point.0 < 0.0 || point.1 < 0.0 {
                        return Err("subglottal time and pressure must be non-negative".into());
                    }
                    if subglottal.last().is_some_and(|last| point.0 <= last.0) {
                        return Err("subglottal times must be strictly increasing".into());
                    }
                    subglottal.push(point);
                }
                other => return Err(format!("unknown key {other:?}")),
            }
        }
        let missing = |what: &str| format!("voice performance needs `{what}`");
        let (output_rate_hz, samples, mechanics_per_output) = audio.ok_or_else(|| missing("audio"))?;
        let (two_mass, fold) = fold.ok_or_else(|| missing("fold"))?;
        if sections.is_empty() {
            return Err(missing("section"));
        }
        if subglottal.is_empty() {
            return Err(missing("subglottal"));
        }
        Ok(Self {
            output_rate_hz,
            samples,
            mechanics_per_output,
            full_scale_pa: full_scale.ok_or_else(|| missing("full_scale_pa"))?,
            gas: gas.ok_or_else(|| missing("ambient"))?,
            two_mass,
            fold,
            wall,
            termination: termination.ok_or_else(|| missing("termination"))?,
            sections,
            tract_source: source.ok_or_else(|| missing("tract_source"))?,
            tract_licence: licence.ok_or_else(|| missing("tract_licence"))?,
            microphone: microphone.ok_or_else(|| missing("microphone"))?,
            subglottal,
        })
    }

    /// The exact shared state filling the tract and the exterior.
    #[must_use]
    pub fn ambient_gas(&self) -> &GasState {
        self.gas.state()
    }

    /// Original source receipts remain attached to the numerical gas view.
    #[must_use]
    pub const fn gas_record(&self) -> &GasRecord {
        &self.gas
    }

    /// Realize the island, the tract's mouth-flow line and the microphone.
    ///
    /// # Errors
    /// Any owner's refusal: the chart licence check, the fold card, a line
    /// that does not converge or is not passive, or an under-resolved band.
    pub fn build(&self) -> Result<VoicePerformance, String> {
        let gas = self.ambient_gas();
        let chart = TractChart::try_new(
            self.sections.clone(),
            &self.tract_source,
            &self.tract_licence,
        )
        .map_err(|e| format!("tract chart refused: {e:?}"))?;
        let duct = chart.to_duct();
        let mechanics_rate = self.output_rate_hz * self.mechanics_per_output;
        let island = GlottalIsland::new_with_wall(
            self.fold,
            self.two_mass,
            &duct,
            gas,
            self.termination,
            mechanics_rate,
            self.wall.as_ref(),
        )
        .map_err(|e| format!("glottal island refused: {e:?}"))?;
        let (mouth, receipt) = mouth_flow_line(
            &duct,
            gas,
            self.termination,
            mechanics_rate,
            island.port_impedance(),
            self.wall.as_ref(),
            &LineLadder::default(),
        )
        .map_err(|e| format!("mouth-flow line refused: {e}"))?;
        let mouth_area = self.sections.last().map_or(0.0, |s| s.area_m2);
        let mouth_radius_m = (mouth_area / core::f64::consts::PI).sqrt();
        let medium = RayleighMedium {
            density: gas.density,
            sound_speed: gas.sound_speed,
        };
        let receiver = BaffledPressure::circular_outlet_multirate_differentiating(
            mouth_radius_m,
            mechanics_rate,
            self.output_rate_hz,
            self.microphone,
            medium,
        )?;
        let passive_formants_hz = chart
            .tmm_formants_terminated(gas, self.wall.as_ref(), self.termination, 4)
            .map_err(|e| format!("{e:?}"));
        let info = VoiceInfo {
            output_rate_hz: self.output_rate_hz,
            mechanics_rate_hz: mechanics_rate,
            samples: self.samples,
            full_scale_pa: self.full_scale_pa,
            two_mass: self.two_mass,
            sections: self.sections.len(),
            tract_length_m: self.sections.iter().map(|s| s.length_m).sum(),
            glottal_radius_m: (self.sections[0].area_m2 / core::f64::consts::PI).sqrt(),
            mouth_radius_m,
            sound_speed_m_s: gas.sound_speed,
            has_wall: self.wall.is_some(),
            mouth_flow_taps: receipt.n_fft,
            mouth_flow_rungs: receipt.rungs.len(),
            passive_formants_hz,
        };
        Ok(VoicePerformance {
            island,
            mouth,
            receiver,
            schedule: Schedule {
                points: self.subglottal.clone(),
                cursor: 0,
            },
            trace: vec![0.0; self.mechanics_per_output as usize],
            ratio: self.mechanics_per_output as usize,
            mechanics_rate,
            info,
            steps: 0,
            completed: 0,
            poisoned: false,
            peak_flow_m3_s: 0.0,
            closed_steps: 0,
        })
    }
}

/// What a realized performance reports about itself.
#[derive(Debug, Clone, PartialEq)]
pub struct VoiceInfo {
    /// Output clock [Hz].
    pub output_rate_hz: u32,
    /// Physical clock [Hz].
    pub mechanics_rate_hz: u32,
    /// Output samples declared.
    pub samples: u64,
    /// Pascals at PCM full scale.
    pub full_scale_pa: f64,
    /// Two-mass (`true`) or one-mass fold.
    pub two_mass: bool,
    /// Tract section count.
    pub sections: usize,
    /// Sum of section lengths [m].
    pub tract_length_m: f64,
    /// Radius of the glottal-end section [m].
    pub glottal_radius_m: f64,
    /// Radius of the mouth-end section, also the radiating opening [m].
    pub mouth_radius_m: f64,
    /// Speed of sound in the declared gas [m/s].
    pub sound_speed_m_s: f64,
    /// Whether a yielding wall was declared.
    pub has_wall: bool,
    /// Converged length of the mouth-flow line.
    pub mouth_flow_taps: usize,
    /// Rungs the mouth-flow ladder traversed.
    pub mouth_flow_rungs: usize,
    /// Passive resonances of exactly this tract, mouth load and wall: the
    /// first four `|Z_in|` peaks seen from a closed glottis
    /// ([`TractChart::tmm_formants_terminated`]), or why it refused. The
    /// glottis opens and closes, so this is the rigid-end limit that the
    /// formants of the radiated sound are compared with, not a prediction of
    /// them.
    pub passive_formants_hz: Result<Vec<f64>, String>,
}

/// Source-side summaries over the physical substeps rendered so far.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VoiceStats {
    /// Largest glottal volume flow [m^3/s].
    pub peak_glottal_flow_m3_s: f64,
    /// Fraction of physical substeps with the folds in contact.
    pub closed_fraction: f64,
}

struct Schedule {
    points: Vec<(f64, f64)>,
    cursor: usize,
}

impl Schedule {
    fn at(&mut self, t: f64) -> f64 {
        let last = self.points.len() - 1;
        if t <= self.points[0].0 {
            return self.points[0].1;
        }
        while self.cursor < last && t >= self.points[self.cursor + 1].0 {
            self.cursor += 1;
        }
        if self.cursor >= last {
            return self.points[last].1;
        }
        let (t0, y0) = self.points[self.cursor];
        let (t1, y1) = self.points[self.cursor + 1];
        y0 + (y1 - y0) * (t - t0) / (t1 - t0)
    }
}

/// A realized performance: the fold, the tract and the microphone, advanced
/// together on one physical clock.
pub struct VoicePerformance {
    island: GlottalIsland,
    mouth: DelayedFilter,
    receiver: BaffledPressure,
    schedule: Schedule,
    trace: Vec<f64>,
    ratio: usize,
    mechanics_rate: u32,
    info: VoiceInfo,
    steps: u64,
    completed: u64,
    poisoned: bool,
    peak_flow_m3_s: f64,
    closed_steps: u64,
}

fn refusal(what: String) -> RenderError {
    RenderError::Voice(AcousticRealizeError::Nonlinear(what))
}

impl VoicePerformance {
    /// The realized description.
    #[must_use]
    pub const fn info(&self) -> &VoiceInfo {
        &self.info
    }

    /// Existing microphone filter latency, in output samples. Physical
    /// source-to-listener travel remains in the pressure and is not included.
    #[must_use]
    pub fn filter_delay_samples(&self) -> f64 {
        self.receiver.filter_delay_samples()
    }

    /// Source-side summaries so far.
    #[must_use]
    pub fn stats(&self) -> VoiceStats {
        VoiceStats {
            peak_glottal_flow_m3_s: self.peak_flow_m3_s,
            closed_fraction: if self.steps == 0 {
                0.0
            } else {
                self.closed_steps as f64 / self.steps as f64
            },
        }
    }

    fn sample(&mut self) -> Result<f64, String> {
        for j in 0..self.ratio {
            let t = self.steps as f64 / f64::from(self.mechanics_rate);
            let drive = self.schedule.at(t);
            let frame = self
                .island
                .step(drive)
                .map_err(|e| format!("glottal island refused: {e:?}"))?;
            let flow = self
                .mouth
                .push(self.island.incident_wave_pa())
                .map_err(|e| format!("mouth-flow line refused: {e:?}"))?;
            self.trace[j] = flow;
            self.peak_flow_m3_s = self.peak_flow_m3_s.max(frame.flow_m3_s);
            if frame.gap_m <= 0.0 {
                self.closed_steps += 1;
            }
            self.steps += 1;
        }
        let pressure = self.receiver.step_trace(&self.trace)?;
        self.completed += 1;
        Ok(pressure)
    }
}

impl PressureRenderer for VoicePerformance {
    fn samples_rendered(&self) -> u64 {
        self.completed
    }

    fn max_block_len(&self) -> usize {
        MAX_BLOCK
    }

    fn validate_sample_rate(&self, sample_rate_hz: u32) -> Result<(), RenderError> {
        if self.poisoned {
            return Err(RenderError::Poisoned);
        }
        if sample_rate_hz != self.info.output_rate_hz {
            return Err(RenderError::Sizing {
                what: "voice performance output clock differs from the requested clock",
            });
        }
        Ok(())
    }

    fn validate_sample_count(&self, samples: u64) -> Result<(), RenderError> {
        if self.poisoned {
            return Err(RenderError::Poisoned);
        }
        if samples > self.info.samples - self.completed {
            return Err(RenderError::Sizing {
                what: "request runs past the samples the voice performance declares",
            });
        }
        Ok(())
    }

    fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        if self.poisoned {
            return Err(RenderError::Poisoned);
        }
        if output.is_empty() {
            return Err(RenderError::EmptyBlock);
        }
        if output.len() > MAX_BLOCK {
            return Err(RenderError::Sizing {
                what: "voice callback exceeds its admitted block length",
            });
        }
        self.validate_sample_count(output.len() as u64)?;
        for slot in output.iter_mut() {
            match self.sample() {
                Ok(pressure) => *slot = pressure,
                Err(what) => {
                    // A failed physical callback leaves a prefix advanced.
                    self.poisoned = true;
                    return Err(refusal(what));
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
#[path = "render/schedule/reed/field_input/tests/support.rs"]
mod gas_test_support;

#[cfg(test)]
mod tests {
    use super::*;

    // Test fixture, not data: a short uniform-ish tract and fold values scaled
    // from the two-mass pack so the build exercises every owner.
    fn fixture(extra: &str) -> String {
        let mut text = String::from(
            "frankensim-voice-performance-v1\n\
             audio 48000 2400 2\n\
             full_scale_pa 2.0\n\
             ambient moist_air 310.15 101325 1.0\n\
             fold two_mass 1.25e-4 2.5e-5 128 12.8 40 0.1 3.57e-4 1.4e-2 2.5e-3\n\
             termination flanged\n\
             tract_source test-fixture\n\
             tract_licence test fixture, not data\n\
             microphone 0 0 0.3 4 16 4000\n\
             subglottal 0 0\n\
             subglottal 0.02 900\n",
        );
        for area in [0.5e-4, 0.6e-4, 1.5e-4, 2.5e-4, 3.0e-4, 3.0e-4, 2.5e-4, 2.0e-4] {
            text.push_str(&format!("section {area:e} 0.02\n"));
        }
        text.push_str(extra);
        text
    }

    fn with_gas(record: &str) -> String {
        fixture("").replace("moist_air 310.15 101325 1.0", record)
    }

    fn pressure(description: &VoiceDescription) -> Vec<f64> {
        let mut performance = description.build().expect("gas reaches all physics owners");
        assert_eq!(
            performance.info().sound_speed_m_s,
            description.ambient_gas().sound_speed
        );
        let mut result = vec![0.0; 1200];
        for block in result.chunks_mut(37) {
            performance.block(block).expect("physical pressure");
        }
        assert!(result.iter().all(|value| value.is_finite()));
        assert!(result.iter().any(|value| value.abs() > 1e-12));
        result
    }

    #[test]
    fn g3_shared_mixture_composition_and_pressure_change_the_radiated_voice() {
        let legacy = VoiceDescription::from_bytes(fixture("").as_bytes()).unwrap();
        assert_eq!(
            *legacy.ambient_gas(),
            GasState::try_new_moist_air(310.15, 101325.0, 1.0).unwrap()
        );
        let records = [
            "mixture 310.15 101325 mole 2 N2 0.8 O2 0.2",
            "mixture 310.15 101325 mole 2 N2 0.8 CO2 0.2",
            "mixture 310.15 80000 mole 2 N2 0.8 O2 0.2",
        ];
        let mut rendered = Vec::new();
        for record in records {
            let description = VoiceDescription::from_bytes(with_gas(record).as_bytes()).unwrap();
            let expected =
                read_gas_record(&mut record.split_ascii_whitespace(), None, None).unwrap();
            assert_eq!(description.gas_record(), &expected);
            rendered.push(pressure(&description));
        }
        assert_ne!(
            rendered[0], rendered[1],
            "composition must reach the audible loop"
        );
        assert_ne!(
            rendered[0], rendered[2],
            "pressure must reach flow and radiation"
        );
    }

    #[test]
    fn g3_nasa9_source_swap_changes_voice_pressure_and_retains_the_original_packs() {
        use fs_material::gas_species::input::ThermochemicalPhase;
        let record = "thermochemical 310.15 101325 mole 1 specimen-gas 1";
        let text = with_gas(record);
        assert!(
            VoiceDescription::from_bytes(text.as_bytes())
                .unwrap_err()
                .contains("explicit source registry")
        );
        let mut rendered = Vec::new();
        for cp_r in [3.5, 4.5] {
            let packs = gas_test_support::packs("specimen-gas", 0.028, cp_r, 200.0);
            let bindings = [packs.binding()];
            let sources = ThermochemicalSources {
                species: &bindings,
                water_species: None,
                phase: ThermochemicalPhase::FrozenGas,
            };
            let description =
                VoiceDescription::from_bytes_with_gas_sources(text.as_bytes(), Some(&sources))
                    .unwrap();
            let expected =
                read_gas_record(&mut record.split_ascii_whitespace(), None, Some(&sources))
                    .unwrap();
            assert_eq!(description.gas_record(), &expected);
            let GasRecord::Thermochemical(gas) = description.gas_record() else {
                panic!("retained source");
            };
            assert_eq!(gas.source_packs()[0].nasa9_pack, packs.nasa9.content_hash());
            rendered.push(pressure(&description));
        }
        assert_ne!(
            rendered[0], rendered[1],
            "data-only heat-capacity swap must reach sound"
        );
    }

    #[test]
    fn gas_refusals_keep_source_admission_strict() {
        for record in [
            "mixture 310.15 101325 mole 2 N2 0.7 O2 0.2",
            "mixture 310.15 101325 mole 1 unavailable 1",
            "mixture 310.15 101325 mole 1 N2 1 extra",
            "moist_air 310.15 101325 1 extra",
            "thermochemical 310.15 101325 mole 1 specimen-gas 1",
        ] {
            assert!(
                VoiceDescription::from_bytes(with_gas(record).as_bytes()).is_err(),
                "{record}"
            );
        }
    }

    #[test]
    fn malformed_and_non_physical_files_are_refused_by_name() {
        for (extra, needle) in [
            ("bogus 1\n", "unknown key"),
            ("audio 48000 10 1\n", "only appear once"),
            ("section 1e-4\n", "takes 2 fields"),
            ("subglottal 0.01 5\n", "strictly increasing"),
            ("microphone 0 0 -1 4 16 1000\n", "in front of the baffle"),
        ] {
            let error = VoiceDescription::from_bytes(fixture(extra).as_bytes()).expect_err(extra);
            assert!(error.contains(needle), "{extra:?} -> {error}");
        }
        let no_licence = fixture("").replace("tract_licence test fixture, not data\n", "");
        let error = VoiceDescription::from_bytes(no_licence.as_bytes()).expect_err("licence");
        assert!(error.contains("tract_licence"), "{error}");
        assert!(VoiceDescription::from_bytes(b"not a voice file\n").is_err());
    }

    #[test]
    fn unlicensed_data_is_refused_at_build() {
        let text = fixture("").replace("test fixture, not data", "unknown");
        let description = VoiceDescription::from_bytes(text.as_bytes()).expect("parses");
        let error = description.build().err().expect("licensing-first refusal");
        assert!(error.contains("tract chart refused"), "{error}");
    }

    #[test]
    fn passive_formants_follow_the_declared_mouth_load() {
        // Eight equal sections are a uniform tube of length L and radius a,
        // closed at the glottis. Ideal-open mouth: F1 = c / (4 L). A flanged
        // mouth lengthens the tube by its end correction, about 0.8216 a.
        let gas = GasState::try_new_moist_air(310.15, 101_325.0, 1.0).expect("gas");
        let area = 2.5e-4;
        let sections = vec![TractSection { area_m2: area, length_m: 0.02 }; 8];
        let chart = TractChart::try_new(sections, "uniform-test-fixture", "test fixture, not data")
            .expect("chart");
        let length = 0.16;
        let radius = (area / core::f64::consts::PI).sqrt();
        let ideal = chart
            .tmm_formants_terminated(&gas, None, Termination::IdealOpen, 1)
            .expect("ideal open")[0];
        let flanged = chart
            .tmm_formants_terminated(&gas, None, Termination::FlangedOpen, 1)
            .expect("flanged")[0];
        let ideal_expected = gas.sound_speed / (4.0 * length);
        let flanged_expected = gas.sound_speed / (4.0 * (length + 0.8216 * radius));
        assert!((ideal / ideal_expected - 1.0).abs() < 0.02, "ideal F1 {ideal} vs {ideal_expected}");
        assert!(
            (flanged / flanged_expected - 1.0).abs() < 0.02,
            "flanged F1 {flanged} vs {flanged_expected}"
        );
        assert!(flanged < ideal, "an end correction must lower the resonance");
        assert_eq!(
            chart.tmm_formants(&gas, None, 1).expect("legacy entry point"),
            vec![ideal],
            "the original oracle is the ideal-open case, unchanged"
        );
    }

    #[test]
    fn subglottal_schedule_interpolates_and_holds() {
        let mut schedule = Schedule {
            points: vec![(0.1, 10.0), (0.2, 30.0)],
            cursor: 0,
        };
        assert_eq!(schedule.at(0.0), 10.0);
        assert!((schedule.at(0.15) - 20.0).abs() < 1.0e-9);
        assert_eq!(schedule.at(0.5), 30.0);
    }

    #[test]
    fn rendering_is_independent_of_callback_partition() {
        let description = VoiceDescription::from_bytes(fixture("").as_bytes()).expect("parses");
        let mut whole = description.build().expect("builds");
        let mut split = description.build().expect("builds");
        let mut a = vec![0.0; 600];
        whole.block(&mut a).expect("whole block");
        let mut b = vec![0.0; 600];
        for chunk in b.chunks_mut(37) {
            split.block(chunk).expect("split block");
        }
        assert!(a.iter().all(|v| v.is_finite()));
        assert!(
            a.iter().zip(&b).all(|(x, y)| x.to_bits() == y.to_bits()),
            "block partition changed the pressure"
        );
        assert_eq!(whole.samples_rendered(), 600);
        let past = split.validate_sample_count(2400);
        assert!(past.is_err(), "a request past the declared end must refuse");
    }
}
