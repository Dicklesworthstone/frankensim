//! Browser surface (wasm-bindgen). Typed refusals become JS exceptions at
//! construction; hot-path refusals are counted and readable via `lastError`.
use crate::kernel::{Config, Description, LivePiano};
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
    piano: LivePiano,
    snap: Vec<f32>,
}

#[wasm_bindgen]
impl PianoEngine {
    /// Build the live piano from the embedded Model D scale and soundboard:
    /// the 88-key string tables (through the reference bank), the hammers and
    /// the Rayleigh receiver. Cold: a few hundred milliseconds.
    /// `board_coordinates`: 0 = the complete embedded board (32 coordinates);
    /// 1..=16 keeps only that many exact low FE modes (cheaper, duller).
    /// Pass 0 for any other numeric argument to take the live default.
    #[wasm_bindgen(constructor)]
    pub fn new(sample_rate: u32, board_coordinates: usize, course_mode_budget: usize, max_dynamic_modes: usize,
        full_scale_pa: f64) -> Result<PianoEngine, JsError> {
        let mut cfg = Config::live(sample_rate);
        if course_mode_budget > 0 { cfg.course_mode_budget = course_mode_budget; }
        if max_dynamic_modes > 0 { cfg.max_dynamic_modes = max_dynamic_modes.min(1 << 24); }
        if full_scale_pa > 0.0 { cfg.full_scale_pa = full_scale_pa; }
        let mut desc = Description::embedded().map_err(|e| JsError::new(&e))?;
        if board_coordinates > 0 && board_coordinates < desc.board.len() {
            desc = desc.truncated(board_coordinates).map_err(|e| JsError::new(&e))?;
        }
        let piano = LivePiano::new(&desc, cfg).map_err(|e| JsError::new(&e))?;
        Ok(Self { piano, snap: Vec::with_capacity(4096) })
    }

    /// The same live piano on a soundboard file written by `examples/gen_board.rs`
    /// (the offline reference's own bridge-informed reduction with more
    /// coordinates than the embedded board, e.g. `model-d-board-2400-r64.bin`),
    /// fetched by the page. Same scale, equations and limits as `new`.
    #[wasm_bindgen(js_name = withBoard)]
    pub fn with_board(sample_rate: u32, board: &[u8], course_mode_budget: usize, max_dynamic_modes: usize,
        full_scale_pa: f64) -> Result<PianoEngine, JsError> {
        let mut cfg = Config::live(sample_rate);
        if course_mode_budget > 0 { cfg.course_mode_budget = course_mode_budget; }
        if max_dynamic_modes > 0 { cfg.max_dynamic_modes = max_dynamic_modes.min(1 << 24); }
        if full_scale_pa > 0.0 { cfg.full_scale_pa = full_scale_pa; }
        let desc = Description::from_board(board).map_err(|e| JsError::new(&e))?;
        let piano = LivePiano::new(&desc, cfg).map_err(|e| JsError::new(&e))?;
        Ok(Self { piano, snap: Vec::with_capacity(4096) })
    }

    /// A second engine in exactly this engine's state: it renders bit for bit
    /// what this one would (so a fresh engine can be timed and its copy played
    /// without building the piano twice). Costs a copy of the engine's memory.
    pub fn copy(&self) -> PianoEngine {
        PianoEngine { piano: self.piano.clone(), snap: Vec::with_capacity(self.snap.capacity()) }
    }

    /// JSON statement of exactly what this model includes (see `describe`).
    pub fn info(&self) -> String { crate::describe(&self.piano) }

    /// Strike: POST-ESCAPEMENT hammer velocity in m/s. False when refused.
    #[wasm_bindgen(js_name = noteOn)]
    pub fn note_on(&mut self, midi: u8, hammer_velocity_m_s: f64) -> bool {
        self.piano.note_on(midi, hammer_velocity_m_s).is_ok()
    }
    #[wasm_bindgen(js_name = noteOff)]
    pub fn note_off(&mut self, midi: u8) { self.piano.note_off(midi); }
    /// Sustain pedal travel in [0,1]; damper drag scales with (1-travel)^2.
    pub fn sustain(&mut self, value: f64) -> bool { self.piano.set_sustain(value).is_ok() }
    #[wasm_bindgen(js_name = allNotesOff)]
    pub fn all_notes_off(&mut self) { self.piano.all_notes_off(); }
    /// Limits for keys that start sounding afterwards.
    #[wasm_bindgen(js_name = setQuality)]
    pub fn set_quality(&mut self, course_mode_budget: usize, partial_band_hz: f64, max_dynamic_modes: usize) -> bool {
        self.piano.set_quality(course_mode_budget, partial_band_hz, max_dynamic_modes).is_ok()
    }

    /// `frames` samples of pressure / full scale (mono Float32Array).
    pub fn render(&mut self, frames: usize) -> Vec<f32> {
        let mut out = vec![0.0_f32; frames.min(8192)];
        self.piano.render(&mut out);
        out
    }

    /// Flat display snapshot; layout documented at `LivePiano::snapshot`.
    pub fn snapshot(&mut self) -> Vec<f32> {
        let mut snap = std::mem::take(&mut self.snap);
        self.piano.snapshot(&mut snap);
        let copy = snap.clone();
        self.snap = snap;
        copy
    }

    /// `[samples, soundingKeys, dynamicModes, dynamicStrings, peakPa, clippedSamples,`
    /// ` forcedRetirements, retirements, contactFailures, rebuilds, noteOns, refusedNoteOns]`.
    pub fn stats(&self) -> Vec<f64> {
        let c = self.piano.counters;
        vec![c.samples as f64, self.piano.sounding_keys() as f64, self.piano.dynamic_modes() as f64,
            self.piano.dynamic_strings() as f64, c.peak_pa, c.clipped_samples as f64,
            c.forced_retirements as f64, c.retirements as f64, c.contact_failures as f64,
            c.rebuilds as f64, c.note_ons as f64, c.refused_note_ons as f64]
    }
    #[wasm_bindgen(js_name = lastError)]
    pub fn last_error(&mut self) -> Option<String> { self.piano.take_error() }
    /// Soundboard quadrature points, interleaved x,y in metres (for drawing).
    #[wasm_bindgen(js_name = boardSurfaceXY)]
    pub fn board_surface_xy(&self) -> Vec<f32> { self.piano.surface_xy().to_vec() }
    /// Bare soundboard mode shapes at those points, row-major [point][mode].
    #[wasm_bindgen(js_name = boardSurfaceShapes)]
    pub fn board_surface_shapes(&self) -> Vec<f32> { self.piano.surface_shapes().to_vec() }
    #[wasm_bindgen(js_name = boardModeHz)]
    pub fn board_mode_hz(&self) -> Vec<f64> { self.piano.board_frequencies_hz().to_vec() }
    /// Per key of the scale, five values: MIDI number, speaking length [m],
    /// strike station as a fraction of that length, unison strings, and the
    /// first stiff-string partial [Hz] from T, EI, mu, L (for drawing/labels).
    #[wasm_bindgen(js_name = keyGeometry)]
    pub fn key_geometry(&self) -> Vec<f64> {
        self.piano.courses().iter().flat_map(|c| [f64::from(c.midi), c.length_m, c.strike_fraction,
            c.unison as f64, c.partial_hz(1, c.tension_n)]).collect()
    }
    /// Stiff-string partial frequencies [Hz] this model retains for a key.
    #[wasm_bindgen(js_name = partialsHz)]
    pub fn partials_hz(&self, midi: u8) -> Vec<f64> { self.piano.partials_hz(midi) }
}
