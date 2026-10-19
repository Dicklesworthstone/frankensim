//! Browser surface (wasm-bindgen). Typed refusals become JS exceptions at
//! construction; hot-path failures poison continuation and are readable via lastError.
use crate::{Config, Description, LivePiano, PianoAudio};
use crate::preparation::{self, PreparationReport, Tuning};
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(js_namespace = console, js_name = error)]
    fn console_error(message: &str);
}

/// A panic aborts the wasm instance; say why on the console first.
#[wasm_bindgen(start)]
pub fn install_panic_report() {
    std::panic::set_hook(Box::new(|info| console_error(&format!("fs-live-piano panic: {info}"))));
}

#[wasm_bindgen]
pub struct PianoEngine {
    piano: Playback,
    snap: Vec<f32>,
    preparation: Option<PreparationReport>,
}

// Keep the legacy single-rate path, including its recoverable note refusals.
// The opt-in multirate adapter has a stricter partial-packet failure policy;
// routing old constructors through it would change default overload behavior.
#[derive(Clone)]
enum Playback { Single(LivePiano), Multirate(PianoAudio) }
impl Playback {
    fn new(desc: &Description, cfg: Config, substeps: usize) -> Result<Self, String> {
        if substeps == 1 { LivePiano::new(desc, cfg).map(Self::Single) }
        else { PianoAudio::new(desc, cfg, substeps).map(Self::Multirate) }
    }
    fn piano(&self) -> &LivePiano {
        match self { Self::Single(p) => p, Self::Multirate(a) => a.piano() }
    }
    fn info(&self) -> String {
        match self { Self::Single(p) => crate::describe(p), Self::Multirate(a) => crate::describe_audio(a) }
    }
    fn note_on(&mut self, midi: u8, speed: f64) -> Result<(), String> {
        match self { Self::Single(p) => p.note_on(midi, speed), Self::Multirate(a) => a.note_on(midi, speed) }
    }
    fn silent_key_down(&mut self, midi: u8) -> Result<(), String> {
        match self { Self::Single(p) => p.silent_key_down(midi), Self::Multirate(a) => a.silent_key_down(midi) }
    }
    fn note_off(&mut self, midi: u8) {
        match self { Self::Single(p) => p.note_off(midi), Self::Multirate(a) => a.note_off(midi) }
    }
    fn all_notes_off(&mut self) {
        match self { Self::Single(p) => p.all_notes_off(), Self::Multirate(a) => a.all_notes_off() }
    }
    fn set_sustain(&mut self, value: f64) -> Result<(), String> {
        match self { Self::Single(p) => p.set_sustain(value), Self::Multirate(a) => a.set_sustain(value) }
    }
    fn set_sostenuto(&mut self, on: bool) -> Result<(), String> {
        match self { Self::Single(p) => { p.set_sostenuto(on); Ok(()) } Self::Multirate(a) => a.set_sostenuto(on) }
    }
    fn set_una_corda(&mut self, on: bool) -> Result<(), String> {
        match self { Self::Single(p) => { p.set_una_corda(on); Ok(()) } Self::Multirate(a) => a.set_una_corda(on) }
    }
    fn set_quality(&mut self, modes: usize, band: f64, total: usize) -> Result<(), String> {
        match self { Self::Single(p) => p.set_quality(modes, band, total),
            Self::Multirate(a) => a.set_quality(modes, band, total) }
    }
    fn render(&mut self, out: &mut [f32]) {
        match self { Self::Single(p) => p.render(out), Self::Multirate(a) => {
            let _ = a.render(out); // retained diagnostic; failed suffix is silent
        } }
    }
    fn snapshot(&mut self, out: &mut Vec<f32>) {
        match self { Self::Single(p) => p.snapshot(out), Self::Multirate(a) => a.snapshot(out) }
    }
    fn take_error(&mut self) -> Option<String> {
        match self { Self::Single(p) => p.take_error(), Self::Multirate(a) => a.take_error() }
    }
    fn is_poisoned(&self) -> bool { matches!(self, Self::Multirate(a) if a.is_poisoned()) }
    fn audio_stats(&self) -> (u64, f64, u64) {
        match self { Self::Single(p) => (p.counters.samples, p.counters.peak_pa, p.counters.clipped_samples),
            Self::Multirate(a) => (a.frames(), a.peak_pa(), a.clipped_frames()) }
    }
}

fn config(sample_rate: u32, course_mode_budget: usize, max_dynamic_modes: usize,
    full_scale_pa: f64) -> Config {
    let mut cfg = Config::live(sample_rate);
    // The browser plays the published RT-0425 string losses and hammer law, K_H on each string
    // with its R_H (Config field docs). The reference admits that law with a supplied scale too:
    // both felt laws key K_H by key number, and the published one only stops dividing it by the unison.
    cfg.rt0425_string_damping = true;
    cfg.rt0425_hammer_stiffness = true;
    cfg.rt0425_hammer_dissipation = true;
    // The silent strings move: every partial to 600 Hz of every string not sounding
    // (671 partials; frankensim-ie2zx). The board and mode cap are timed with them.
    cfg.silent_bank_hz = 600.0;
    if course_mode_budget > 0 { cfg.course_mode_budget = course_mode_budget; }
    if max_dynamic_modes > 0 { cfg.max_dynamic_modes = max_dynamic_modes.min(1 << 24); }
    if full_scale_pa > 0.0 { cfg.full_scale_pa = full_scale_pa; }
    cfg
}

#[wasm_bindgen]
impl PianoEngine {
    /// Original single-rate constructor: embedded board or its exact low prefix.
    /// Pass 0 for a numeric quality argument to take the live default.
    #[wasm_bindgen(constructor)]
    pub fn new(sample_rate: u32, board_coordinates: usize, course_mode_budget: usize, max_dynamic_modes: usize,
        full_scale_pa: f64) -> Result<PianoEngine, JsError> {
        Self::with_oversampling(sample_rate, board_coordinates, course_mode_budget,
            max_dynamic_modes, full_scale_pa, 1)
    }

    /// Run the original mechanics and Rayleigh observer at 1x, 2x or 4x the
    /// host rate, then anti-alias every pressure substep before host output.
    /// The output partial band and physical source values are unchanged.
    /// This is not interpolation of a low-rate piano recording.
    #[wasm_bindgen(js_name = withOversampling)]
    pub fn with_oversampling(sample_rate: u32, board_coordinates: usize, course_mode_budget: usize,
        max_dynamic_modes: usize, full_scale_pa: f64, substeps: usize) -> Result<PianoEngine, JsError> {
        let cfg = config(sample_rate, course_mode_budget, max_dynamic_modes, full_scale_pa);
        let mut desc = Description::embedded().map_err(|e| JsError::new(&e))?;
        if board_coordinates > 0 && board_coordinates < desc.board.len() {
            desc = desc.truncated(board_coordinates).map_err(|e| JsError::new(&e))?;
        }
        let piano = Playback::new(&desc, cfg, substeps).map_err(|e| JsError::new(&e))?;
        Ok(Self { piano, snap: Vec::with_capacity(4096), preparation: None })
    }

    /// Original single-rate external-board constructor, unchanged arguments.
    #[wasm_bindgen(js_name = withBoard)]
    pub fn with_board(sample_rate: u32, board: &[u8], course_mode_budget: usize, max_dynamic_modes: usize,
        full_scale_pa: f64) -> Result<PianoEngine, JsError> {
        Self::with_board_oversampling(sample_rate, board, course_mode_budget,
            max_dynamic_modes, full_scale_pa, 1)
    }

    /// The same higher-rate path with a supplied physical board reduction.
    /// Board coordinates, damping and surface fields are not resynthesized.
    #[wasm_bindgen(js_name = withBoardOversampling)]
    pub fn with_board_oversampling(sample_rate: u32, board: &[u8], course_mode_budget: usize,
        max_dynamic_modes: usize, full_scale_pa: f64, substeps: usize) -> Result<PianoEngine, JsError> {
        let cfg = config(sample_rate, course_mode_budget, max_dynamic_modes, full_scale_pa);
        let desc = Description::from_board(board).map_err(|e| JsError::new(&e))?;
        let piano = Playback::new(&desc, cfg, substeps).map_err(|e| JsError::new(&e))?;
        Ok(Self { piano, snap: Vec::with_capacity(4096), preparation: None })
    }

    /// Prepare physical courses BEFORE all live string/bridge/hammer lowering.
    /// Empty board bytes select the embedded board and its optional exact prefix;
    /// a supplied board uses its complete basis and requires coordinates = 0.
    /// Empty scale text preserves existing courses. Empty tuning preserves their
    /// tensions; otherwise the packet is A4, bass[4], treble[4], band (0 = none).
    /// Imported CSVs are not retuned unless that packet is explicitly supplied.
    #[wasm_bindgen(js_name = withPreparation)]
    pub fn with_preparation(sample_rate: u32, board: &[u8], board_coordinates: usize,
        course_mode_budget: usize, max_dynamic_modes: usize, full_scale_pa: f64,
        substeps: usize, scale_csv: &str, tuning: &[f64]) -> Result<PianoEngine, JsError> {
        let policy = Tuning::from_values(tuning).map_err(|e| JsError::new(&e))?;
        if !board.is_empty() && board_coordinates != 0 {
            return Err(JsError::new("supplied board preparation requires board_coordinates = 0"));
        }
        let mut desc = if board.is_empty() { Description::embedded() }
            else { Description::from_board(board) }.map_err(|e| JsError::new(&e))?;
        if board_coordinates > 0 && board_coordinates < desc.board.len() {
            desc = desc.truncated(board_coordinates).map_err(|e| JsError::new(&e))?;
        }
        let report = preparation::prepare(&mut desc,
            if scale_csv.is_empty() { None } else { Some(scale_csv) }, policy)
            .map_err(|e| JsError::new(&e))?;
        let mut cfg = config(sample_rate, course_mode_budget, max_dynamic_modes, full_scale_pa);
        // RT-0425's per-key losses describe the Model D's own strings. An imported scale keeps the
        // authored law, as the reference refuses --rt0425-string-damping with a supplied scale;
        // retuning the built-in scale (tension only, same strings) keeps the published losses.
        cfg.rt0425_string_damping = scale_csv.is_empty();
        let piano = Playback::new(&desc, cfg, substeps).map_err(|e| JsError::new(&e))?;
        Ok(Self { piano, snap: Vec::with_capacity(4096), preparation: Some(report) })
    }

    /// Copy the physical engine AND pressure-filter history; a timed copy and
    /// an unplayed copy retain the same selected clock and quality policy.
    pub fn copy(&self) -> PianoEngine {
        PianoEngine { piano: self.piano.clone(), snap: Vec::with_capacity(self.snap.capacity()),
            preparation: self.preparation.clone() }
    }

    /// JSON includes output and mechanical rates and total output-frame latency.
    pub fn info(&self) -> String {
        let mut info = self.piano.info();
        if let Some(report) = &self.preparation {
            // Both existing description owners produce one complete JSON object.
            let _ = info.pop();
            info.push_str(",\"scalePreparation\":");
            info.push_str(&report.json());
            info.push('}');
        }
        info
    }

    /// Strike: POST-ESCAPEMENT hammer velocity in m/s. False when refused.
    #[wasm_bindgen(js_name = noteOn)]
    pub fn note_on(&mut self, midi: u8, hammer_velocity_m_s: f64) -> bool {
        self.piano.note_on(midi, hammer_velocity_m_s).is_ok()
    }
    /// Hold a key down without striking it: its dampers lift and its strings, attached
    /// from rest, ring along with what the board drives into them. False when refused.
    #[wasm_bindgen(js_name = silentKeyDown)]
    pub fn silent_key_down(&mut self, midi: u8) -> bool { self.piano.silent_key_down(midi).is_ok() }
    #[wasm_bindgen(js_name = noteOff)]
    pub fn note_off(&mut self, midi: u8) { self.piano.note_off(midi); }
    /// Sustain pedal travel in [0,1]; damper drag scales with (1-travel)^2.
    pub fn sustain(&mut self, value: f64) -> bool { self.piano.set_sustain(value).is_ok() }
    /// Sostenuto pedal: pressing it latches the dampers of the keys held down then
    /// (and of all keys while the sustain pedal is fully down); releasing it drops them.
    pub fn sostenuto(&mut self, on: bool) -> bool { self.piano.set_sostenuto(on).is_ok() }
    /// Una corda pedal: strikes that start afterwards miss the last string of each course.
    #[wasm_bindgen(js_name = unaCorda)]
    pub fn una_corda(&mut self, on: bool) -> bool { self.piano.set_una_corda(on).is_ok() }
    #[wasm_bindgen(js_name = allNotesOff)]
    pub fn all_notes_off(&mut self) { self.piano.all_notes_off(); }
    /// Limits for keys that start sounding afterwards, bounded by OUTPUT rate.
    #[wasm_bindgen(js_name = setQuality)]
    pub fn set_quality(&mut self, course_mode_budget: usize, partial_band_hz: f64, max_dynamic_modes: usize) -> bool {
        self.piano.set_quality(course_mode_budget, partial_band_hz, max_dynamic_modes).is_ok()
    }

    /// `frames` host samples of pressure / full scale (mono Float32Array).
    /// A failed suffix is silent and its diagnostic stays in lastError;
    /// reading the diagnostic cannot resume a partially advanced performance.
    pub fn render(&mut self, frames: usize) -> Vec<f32> {
        let mut out = vec![0.0_f32; frames.min(8192)];
        self.piano.render(&mut out);
        out
    }

    /// Flat physical-source display snapshot, not a second audio rendering.
    pub fn snapshot(&mut self) -> Vec<f32> {
        let mut snap = std::mem::take(&mut self.snap);
        self.piano.snapshot(&mut snap);
        let copy = snap.clone();
        self.snap = snap;
        copy
    }

    /// `[outputFrames, soundingKeys, dynamicModes, dynamicStrings, outputPeakPa, clippedOutputFrames,`
    /// ` forcedRetirements, retirements, contactFailures, rebuilds, noteOns, refusedNoteOns]`.
    /// The array layout is unchanged; acoustic quantities use the host clock.
    pub fn stats(&self) -> Vec<f64> {
        let p = self.piano.piano();
        let c = p.counters;
        let (frames, peak, clipped) = self.piano.audio_stats();
        vec![frames as f64, p.sounding_keys() as f64, p.dynamic_modes() as f64,
            p.dynamic_strings() as f64, peak, clipped as f64,
            c.forced_retirements as f64, c.retirements as f64, c.contact_failures as f64,
            c.rebuilds as f64, c.note_ons as f64, c.refused_note_ons as f64]
    }
    #[wasm_bindgen(js_name = lastError)]
    pub fn last_error(&mut self) -> Option<String> { self.piano.take_error() }
    #[wasm_bindgen(js_name = isPoisoned)]
    pub fn is_poisoned(&self) -> bool { self.piano.is_poisoned() }
    /// Soundboard quadrature points, interleaved x,y in metres (for drawing).
    #[wasm_bindgen(js_name = boardSurfaceXY)]
    pub fn board_surface_xy(&self) -> Vec<f32> { self.piano.piano().surface_xy().to_vec() }
    /// Bare soundboard mode shapes at those points, row-major [point][mode].
    #[wasm_bindgen(js_name = boardSurfaceShapes)]
    pub fn board_surface_shapes(&self) -> Vec<f32> { self.piano.piano().surface_shapes().to_vec() }
    #[wasm_bindgen(js_name = boardModeHz)]
    pub fn board_mode_hz(&self) -> Vec<f64> { self.piano.piano().board_frequencies_hz().to_vec() }
    /// Per key: MIDI, length [m], strike fraction, unison count, first partial [Hz].
    #[wasm_bindgen(js_name = keyGeometry)]
    pub fn key_geometry(&self) -> Vec<f64> {
        self.piano.piano().courses().iter().flat_map(|c| [f64::from(c.midi), c.length_m, c.strike_fraction,
            c.unison as f64, c.partial_hz(1, c.tension_n)]).collect()
    }
    /// Stiff-string partial frequencies [Hz] this model retains for a key.
    #[wasm_bindgen(js_name = partialsHz)]
    pub fn partials_hz(&self, midi: u8) -> Vec<f64> { self.piano.piano().partials_hz(midi) }
}
