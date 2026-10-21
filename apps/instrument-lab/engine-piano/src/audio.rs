//! Output-rate adapter for the EXISTING live piano, not another integrator.
//! The whole string/board/hammer/felt engine and its Rayleigh observer run at
//! the selected mechanical clock. Every resulting f64 pressure sample reaches
//! the shared anti-alias decimator before conversion to the host's Float32.
//! There is no interpolation of finished audio, endpoint sampling, gain fit,
//! retuning, new constitutive law or new filter design here.
use crate::{Config, Description, Ledger, LivePiano};
use fs_couple::pcm_wav::decimate::Decimator;
use fs_material::gas::GasState;

/// One owned live engine and its persistent output-rate observation history.
/// Factors 1, 2 and 4 are explicit; total mechanical rate remains <=192 kHz.
/// At factor one, successful pressure and Float32 samples retain the original
/// arithmetic. Higher factors consume proportionally more mechanical work.
/// The engine's voice-age/retirement counters remain mechanical-step counters,
/// as in a direct high-rate LivePiano; this adapter does not retime that policy.
#[derive(Clone)]
pub struct PianoAudio {
    piano: LivePiano,
    output_rate: u32,
    substeps: usize,
    decimator: Decimator,
    pressure: [f64; 4],
    frames: u64,
    peak_pa: f64,
    clipped_frames: u64,
    /// A partial physical/observer failure permanently prevents continuation.
    failure: Option<String>,
    last_error: Option<String>,
}

impl PianoAudio {
    /// `cfg.sample_rate` is the OUTPUT rate; the partial band remains in Hz
    /// and cannot exceed 45% of that output rate. All physical source data are
    /// passed unchanged to LivePiano, whose time-dependent maps use the higher
    /// clock. The original room-air convenience constructor owns the medium.
    pub fn new(desc: &Description, cfg: Config, substeps: usize) -> Result<Self, String> {
        Self::build(desc, cfg, substeps, None, None)
    }

    /// Same preparation with an explicitly admitted physical acoustic gas.
    pub fn new_with_gas(desc: &Description, cfg: Config, substeps: usize, gas: GasState)
        -> Result<Self, String> {
        Self::build(desc, cfg, substeps, Some(gas), None)
    }

    /// Canonical source-bound strings follow the same complete mechanical
    /// clock and pressure decimator at 1×, 2× or 4×. No material path bypasses
    /// the original observer or reconstructs coefficients after preparation.
    pub fn new_with_materials(desc: &Description, cfg: Config, substeps: usize, gas: GasState,
        materials: std::sync::Arc<crate::preparation::materials::StringMaterials>)
        -> Result<Self, String> {
        Self::build(desc, cfg, substeps, Some(gas), Some(materials))
    }

    fn build(desc: &Description, mut cfg: Config, substeps: usize, gas: Option<GasState>,
        materials: Option<std::sync::Arc<crate::preparation::materials::StringMaterials>>)
        -> Result<Self, String> {
        let output_rate = cfg.sample_rate;
        if !matches!(substeps, 1 | 2 | 4) || !(8_000..=192_000).contains(&output_rate) {
            return Err("piano audio needs 1, 2 or 4 substeps and an output rate in 8..192 kHz".into());
        }
        let mechanical_rate = output_rate.checked_mul(substeps as u32)
            .filter(|rate| *rate <= 192_000)
            .ok_or("piano mechanical rate exceeds the existing 192 kHz engine limit")?;
        if !cfg.partial_band_hz.is_finite() || cfg.partial_band_hz <= 0.0
            || cfg.partial_band_hz > 0.45 * f64::from(output_rate) {
            return Err("piano partial band must stay within 45% of the OUTPUT rate".into());
        }
        // Clock/filter admission precedes any physical preparation.
        let decimator = Decimator::new(substeps, 1).map_err(str::to_owned)?;
        cfg.sample_rate = mechanical_rate;
        let piano = match (gas, materials) {
            (Some(gas), Some(materials)) => LivePiano::new_with_materials(desc, cfg, gas, materials)?,
            (Some(gas), None) => LivePiano::new_with_gas(desc, cfg, gas)?,
            (None, None) => LivePiano::new(desc, cfg)?,
            (None, Some(_)) => return Err("material playback requires an admitted acoustic gas".into()),
        };
        Ok(Self { piano, output_rate, substeps, decimator, pressure: [0.0; 4],
            frames: 0, peak_pa: 0.0, clipped_frames: 0, failure: None, last_error: None })
    }

    /// Immutable diagnostics only: callers cannot advance the source while
    /// bypassing the decimator or replace the clock of a running performance.
    pub fn piano(&self) -> &LivePiano { &self.piano }
    pub fn output_rate(&self) -> u32 { self.output_rate }
    pub fn mechanical_rate(&self) -> u32 { self.piano.config().sample_rate }
    pub fn substeps(&self) -> usize { self.substeps }
    /// Successfully emitted pressure frames, not mechanical samples.
    pub fn frames(&self) -> u64 { self.frames }
    pub fn peak_pa(&self) -> f64 { self.peak_pa }
    pub fn clipped_frames(&self) -> u64 { self.clipped_frames }
    pub fn is_poisoned(&self) -> bool { self.failure.is_some() }
    pub fn decimation_delay_frames(&self) -> f64 { self.decimator.delay_output_frames() }

    /// Propagation FIR support and numerical filter delay in OUTPUT frames.
    /// The first value follows the original receiver's support convention;
    /// it is not a measured microphone group delay. Decimation latency is
    /// included once, with no trimming, pre-roll or future-sample compensation.
    pub fn microphone_delay_frames(&self) -> (f64, f64) {
        let (travel, filter) = self.piano.microphone_delay_samples();
        let ratio = self.substeps as f64;
        (travel as f64 / ratio, filter / ratio + self.decimation_delay_frames())
    }

    pub fn multiplies_per_output_frame(&self) -> usize {
        self.substeps * self.piano.microphone_multiplies()
            + self.decimator.multiplies_per_output_frame()
    }

    fn healthy(&self) -> Result<(), String> {
        match &self.failure {
            Some(error) => Err(format!("piano audio cannot continue after partial failure: {error}")),
            None => Ok(()),
        }
    }
    fn fail(&mut self, error: String) -> String {
        self.failure = Some(error.clone());
        self.last_error = Some(error.clone());
        error
    }
    /// Reading the diagnostic does NOT unpoison the source/observer pair.
    pub fn take_error(&mut self) -> Option<String> { self.last_error.take() }

    pub fn note_on(&mut self, midi: u8, velocity_m_s: f64) -> Result<(), String> {
        self.healthy()?;
        self.piano.note_on(midi, velocity_m_s)
    }
    pub fn silent_key_down(&mut self, midi: u8) -> Result<(), String> {
        self.healthy()?;
        self.piano.silent_key_down(midi)
    }
    pub fn note_off(&mut self, midi: u8) { self.piano.note_off(midi); }
    pub fn all_notes_off(&mut self) { self.piano.all_notes_off(); }
    pub fn set_sustain(&mut self, value: f64) -> Result<(), String> {
        self.healthy()?;
        self.piano.set_sustain(value)
    }
    pub fn set_sostenuto(&mut self, on: bool) -> Result<(), String> {
        self.healthy()?;
        self.piano.set_sostenuto(on);
        Ok(())
    }
    pub fn set_una_corda(&mut self, on: bool) -> Result<(), String> {
        self.healthy()?;
        self.piano.set_una_corda(on);
        Ok(())
    }
    pub fn set_quality(&mut self, course_modes: usize, band_hz: f64, total_modes: usize)
        -> Result<(), String> {
        self.healthy()?;
        if !band_hz.is_finite() || band_hz <= 0.0 || band_hz > 0.45 * f64::from(self.output_rate) {
            return Err("piano partial band must stay within 45% of the OUTPUT rate".into());
        }
        self.piano.set_quality(course_modes, band_hz, total_modes)
    }
    pub fn snapshot(&mut self, out: &mut Vec<f32>) { self.piano.snapshot(out); }
    pub fn close_ledger(&mut self) -> Ledger { self.piano.close_ledger() }

    /// Emit one pressure sample [Pa] after ALL substeps of the original
    /// instrument and high-rate receiver. Events queued before this call act
    /// at its first mechanical substep, not after the output-rate group.
    /// Failure after physical progress poisons the adapter: a missing acoustic
    /// packet may not be silently retried against already-advanced mechanics.
    pub fn step_pressure(&mut self) -> Result<f64, String> {
        self.healthy()?;
        let next = self.frames.checked_add(1).ok_or("piano output clock exhausted")?;
        if self.piano.counters.samples.checked_add(self.substeps as u64).is_none() {
            return Err("piano mechanical clock exhausted".into());
        }
        for index in 0..self.substeps {
            let value = self.piano.step();
            if let Some(error) = self.piano.take_error() {
                return Err(self.fail(error));
            }
            if !value.is_finite() {
                return Err(self.fail("nonfinite high-rate piano pressure".into()));
            }
            self.pressure[index] = value;
        }
        let pressure = match self.decimator.preview(&self.pressure[..self.substeps]) {
            Ok(output) => output[0],
            Err(error) => return Err(self.fail(error.into())),
        };
        self.decimator.commit();
        self.frames = next;
        self.peak_pa = self.peak_pa.max(pressure.abs());
        if pressure.abs() >= self.piano.config().full_scale_pa { self.clipped_frames += 1; }
        Ok(pressure)
    }

    /// Float output uses the original fixed Pa/full-scale conversion once,
    /// after anti-alias filtering. On error the unproduced suffix is silent,
    /// the successful prefix remains, and the original diagnostic is returned.
    /// Counts refer to pressure frames; a final Float32 overflow can therefore
    /// follow an accepted pressure frame. No peak normalization or clipping.
    pub fn render(&mut self, out: &mut [f32]) -> Result<(), String> {
        let scale = 1.0 / self.piano.config().full_scale_pa;
        for index in 0..out.len() {
            let sample = match self.step_pressure() {
                Ok(pressure) => (pressure * scale) as f32,
                Err(error) => { out[index..].fill(0.0); return Err(error); }
            };
            if !sample.is_finite() {
                out[index..].fill(0.0);
                return Err(self.fail("piano pressure exceeds finite Float32 output".into()));
            }
            out[index] = sample;
        }
        self.healthy()
    }
}
