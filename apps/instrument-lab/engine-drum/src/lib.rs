//! fs-live-drum: the live, in-browser two-head drum.
//!
//! This crate is the REDUCED REAL-TIME MODEL shown on the FrankenSim demo
//! site, explicitly separate from the offline reference render. It contains no
//! physics of its own: it ports the minimal `drum-modal` construction of
//! `crates/fs-couple/examples/percussion` onto the public library API and adds
//! only the live-session bookkeeping (stick lifecycle, energy ledger, typed
//! refusals, a flat browser boundary).
//!
//! | Stage | Owner |
//! |---|---|
//! | Circular film mesh, DKT bending + membrane prestress pencil | `fs-plate` `TensionedDisk` |
//! | Eigenmodes inside the declared window | `fs-modal` `slice_window` |
//! | Exact-ZOH modal bodies + sealed-air spring + tip contact | `fs-couple` `LinearImpactSystem` |
//! | Hertz (optionally Hunt-Crossley) contact potential | `fs-dcontact` `Obstacle` |
//! | Pressure at the microphones, anti-aliased to 48 kHz | `fs-couple` `BaffledPressure` + `Decimator` |
//!
//! See `live::PROVENANCE` and `live::NO_CLAIM` for the exact scope.

pub mod engine;
pub mod instrument;
pub mod json;
pub mod live;
pub mod params;

pub use engine::{BenchReport, Ledger, Snapshot, StrikeOutcome, Voice};
pub use instrument::Instrument;
pub use live::LiveDrumCore;
pub use params::{Atmosphere, DrumParams, Refusal};

/// Which fused multiply-add executes FrankenSim's `f64::mul_add` chains in
/// this build. Results are bit-identical on every path; only speed differs.
#[must_use]
pub fn fma_path() -> &'static str {
    #[cfg(all(target_arch = "wasm32", feature = "hardware-fma"))]
    {
        match fs_live_fma::active() {
            fs_live_fma::FmaPath::Hardware => "hardware (wasm relaxed-simd madd, probe-verified fused)",
            fs_live_fma::FmaPath::Software => "software (libm; this host's relaxed madd does not fuse)",
        }
    }
    #[cfg(all(target_arch = "wasm32", not(feature = "hardware-fma")))]
    {
        "software (compiler-builtins)"
    }
    #[cfg(not(target_arch = "wasm32"))]
    {
        "native"
    }
}

#[cfg(target_arch = "wasm32")]
mod wasm {
    use crate::engine::StrikeOutcome;
    use crate::live::{LiveDrumCore, flatten_snapshot};
    use crate::params::{DrumParams, Refusal};
    use wasm_bindgen::prelude::*;

    /// A typed refusal as a JS `Error` with `name = "FsRefusal"` and a
    /// machine-readable `code` property.
    fn refuse(refusal: &Refusal) -> JsValue {
        let error = js_sys::Error::new(&refusal.message);
        error.set_name("FsRefusal");
        let _ = js_sys::Reflect::set(&error, &JsValue::from_str("code"), &JsValue::from_str(refusal.code));
        error.into()
    }

    /// The live drum. All inputs are SI; see `info()` for the built instrument.
    #[wasm_bindgen]
    pub struct LiveDrum {
        core: LiveDrumCore,
    }

    #[wasm_bindgen]
    impl LiveDrum {
        /// Build from the reference defaults overlaid with a flat JSON object
        /// (`"{}"` builds the offline reference's 14 x 6.5 inch drum).
        #[wasm_bindgen(constructor)]
        pub fn new(params_json: &str) -> Result<LiveDrum, JsValue> {
            LiveDrumCore::new(params_json).map(|core| Self { core }).map_err(|r| refuse(&r))
        }

        /// The same, with the atmosphere resolved from `gas_record` through
        /// its source packet: the live kit's `withGasSources` framing. A
        /// `thermochemical ...` record needs it.
        #[wasm_bindgen(js_name = withGasSources)]
        pub fn with_gas_sources(params_json: &str, gas_record: &str, bytes: &[u8], lengths: &[u32],
            water_species: &str, phase: &str) -> Result<LiveDrum, JsValue> {
            LiveDrumCore::with_gas_sources(params_json, gas_record, bytes, lengths, water_species, phase)
                .map(|core| Self { core }).map_err(|r| refuse(&r))
        }

        /// Rebuild from the current parameters overlaid with a flat JSON
        /// object. Returns true when the heads were re-meshed and re-solved.
        /// The new instrument starts at rest and silent. On refusal the old
        /// instrument is untouched.
        pub fn rebuild(&mut self, params_json: &str) -> Result<bool, JsValue> {
            self.core.rebuild(params_json).map_err(|r| refuse(&r))
        }

        /// Same instrument, at rest and silent. Clears a halted voice.
        pub fn reset(&mut self) -> Result<(), JsValue> {
            self.core.reset().map_err(|r| refuse(&r))
        }

        /// Launch the stick at a head-centred station [m] with a downward tip
        /// speed [m/s]. Returns "queued" or "deferred" (tip still in contact:
        /// the launch waits for separation).
        pub fn strike(&mut self, x_m: f64, y_m: f64, speed_m_s: f64) -> Result<String, JsValue> {
            match self.core.strike(x_m, y_m, speed_m_s) {
                Ok(StrikeOutcome::Queued) => Ok("queued".into()),
                Ok(StrikeOutcome::DeferredUntilSeparation) => Ok("deferred".into()),
                Err(r) => Err(refuse(&r)),
            }
        }

        /// Render `frames` frames at 48 kHz. Planar float32:
        /// `[batter-side microphone | resonant-side microphone]`, each sample
        /// = pressure [Pa] / `full_scale_pa`. Never throws: after a mechanics
        /// refusal the voice halts, the output is exact silence, and
        /// `take_refusal()` reports the cause once.
        pub fn render(&mut self, frames: usize) -> Vec<f32> {
            self.core.render_f32(frames.min(16_384))
        }

        /// JSON `{code,message}` of the refusal that halted the voice, once.
        pub fn take_refusal(&mut self) -> Option<String> {
            self.core.take_refusal().map(|r| r.to_json())
        }

        /// Flat snapshot: see `info().snapshot_layout` for the scalar header;
        /// modal q then v (batter, resonant) follow.
        pub fn snapshot(&mut self) -> Vec<f64> {
            let snapshot = self.core.snapshot();
            let inst = self.core.instrument();
            flatten_snapshot(&snapshot, inst.heads[0].mode_count(), inst.heads[1].mode_count())
        }

        /// Transverse displacement of every mesh node [m], positive DOWNWARD
        /// for both heads: `[batter nodes | resonant nodes]`. Computed from
        /// the same accepted state as `snapshot()`; read-only.
        pub fn displacements(&self) -> Vec<f32> {
            self.core.node_displacements()
        }

        /// Static description of the built instrument (JSON).
        pub fn info(&self) -> String {
            self.core.info_json()
        }

        /// Mesh node coordinates `[x0,y0,x1,y1,...]` [m], shared by both heads.
        pub fn mesh_nodes(&self) -> Vec<f32> {
            self.core.instrument().nodes.iter().flat_map(|p| [p[0] as f32, p[1] as f32]).collect()
        }

        /// Mesh triangles `[a0,b0,c0,...]`, counter-clockwise seen from above.
        pub fn mesh_triangles(&self) -> Vec<u32> {
            self.core.instrument().tris.iter().flatten().copied().collect()
        }

        /// Mode shapes of head 0 (batter) or 1 (resonant): node-major
        /// `[node * modes + mode]`, transverse displacement per unit modal
        /// coordinate [1/sqrt(kg)], positive downward.
        pub fn mode_shapes(&self, head: usize) -> Result<Vec<f32>, JsValue> {
            let inst = self.core.instrument();
            let head = inst.heads.get(head)
                .ok_or_else(|| refuse(&Refusal::new("bad_parameter", "head index must be 0 (batter) or 1 (resonant)")))?;
            Ok(head.node_shapes.iter().map(|v| *v as f32).collect())
        }

        /// In-vacuo eigenfrequencies [Hz] of head 0 (batter) or 1 (resonant).
        pub fn frequencies_hz(&self, head: usize) -> Result<Vec<f64>, JsValue> {
            let inst = self.core.instrument();
            inst.heads.get(head).map(|h| h.frequencies_hz())
                .ok_or_else(|| refuse(&Refusal::new("bad_parameter", "head index must be 0 (batter) or 1 (resonant)")))
        }

        /// Render a scratch voice of this instrument at `substeps` for
        /// `frames` frames with a stroke every `strike_interval` frames (0 =
        /// one stroke). The running voice is untouched. Time this call from
        /// JS. Returns `[frames, mechanical_steps, steps_with_stick,
        /// steps_in_contact, strikes, peak_pa, balance_residual_j]`.
        pub fn bench(
            &self, substeps: usize, frames: usize, strike_interval: usize, x_m: f64, y_m: f64, speed_m_s: f64,
        ) -> Result<Vec<f64>, JsValue> {
            let r = self.core
                .bench(substeps, frames.min(480_000), strike_interval, x_m, y_m, speed_m_s)
                .map_err(|r| refuse(&r))?;
            Ok(vec![
                r.frames as f64, r.mechanical_steps as f64, r.steps_with_stick as f64,
                r.steps_in_contact as f64, r.strikes as f64, r.peak_pa, r.ledger.balance_residual_j,
            ])
        }
    }

    /// Bit-compare the hardware and software FMA on `count` adversarial
    /// operand triples. Returns the number of mismatches (0 expected whenever
    /// `fma` reports hardware), or -1 in the baseline build.
    #[wasm_bindgen]
    pub fn fma_self_test(count: u32) -> i32 {
        #[cfg(feature = "hardware-fma")]
        {
            let mut state = 0x9E37_79B9_7F4A_7C15_u64;
            let mut next = move || {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                state
            };
            let mut mismatches = 0;
            for i in 0..count {
                // Random significands with exponents confined to +-200 so the
                // products stay finite, plus near-cancelling addends: the
                // cases where one rounding and two roundings differ.
                let make = |bits: u64, spread: u64| {
                    let exponent = 1023 - spread / 2 + (bits >> 52) % spread;
                    f64::from_bits((bits & 0x800F_FFFF_FFFF_FFFF) | (exponent << 52))
                };
                let a = make(next(), 400);
                let b = make(next(), 400);
                let c = match i % 4 {
                    0 => make(next(), 800),
                    1 => -(a * b),
                    2 => -(a * b) * (1.0 + f64::EPSILON * ((next() % 7) as f64 - 3.0)),
                    _ => make(next(), 40) * a * b,
                };
                let (hard, soft) = fs_live_fma::both(a, b, c);
                if hard.to_bits() != soft.to_bits() && !(hard.is_nan() && soft.is_nan()) {
                    mismatches += 1;
                }
            }
            mismatches
        }
        #[cfg(not(feature = "hardware-fma"))]
        {
            let _ = count;
            -1
        }
    }

    /// The reference defaults as the flat JSON object the constructor accepts.
    #[wasm_bindgen]
    pub fn default_params() -> String {
        DrumParams::default().to_json()
    }

    /// Provenance and no-claim text, `{"provenance":..,"no_claim":..}`.
    #[wasm_bindgen]
    pub fn scope() -> String {
        format!(
            "{{\"provenance\":{},\"no_claim\":{}}}",
            crate::json::string(crate::live::PROVENANCE),
            crate::json::string(crate::live::NO_CLAIM)
        )
    }
}
