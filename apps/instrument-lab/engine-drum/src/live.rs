//! The platform-neutral live drum: one built instrument plus one running
//! voice. The wasm-bindgen surface in lib.rs is a thin shell over this type,
//! so every behaviour here is exercised natively by `cargo test` and the
//! `bench` example.

use crate::engine::{BenchReport, Snapshot, StrikeOutcome, Voice};
use crate::instrument::{HERTZ_ALPHA, Instrument, LAUNCH_GAP_M};
use crate::json;
use crate::params::{DrumParams, OUTPUT_RATE_HZ, Refusal};

/// What the site may say about this engine, verbatim.
pub const PROVENANCE: &str = "Reduced real-time model. FrankenSim's own code compiled to WebAssembly: \
circular DKT plate + membrane-prestress finite-element heads (fs-plate), inertia-certified shift-invert \
eigenmodes (fs-modal), the prepared exact zero-order-hold linear impact image with a simultaneously \
solved sealed-air volume spring and unilateral tip contact (fs-couple LinearImpactSystem, the offline \
`drum-modal` image), an energy-consistent Hertz contact potential (fs-dcontact), and a time-domain \
Rayleigh half-space pressure receiver with the library decimator (fs-couple BaffledPressure). \
One fs-material gas state supplies the cavity bulk modulus, receiver density and sound speed. \
Frequencies follow the declared heads and gas spring; pressure follows their driven motion and the \
declared observation geometry and atmosphere, with fixed full-scale pressure.";

/// What this engine does NOT establish, verbatim.
pub const NO_CLAIM: &str = "Not the offline reference and not a calibrated instrument. Omitted here: \
boundary-element radiation of the closed drum (replaced by two independent infinite-baffle half-spaces, \
one per head), radiation loading on the heads, geometric head stretching (pitch glide), snare wires, \
shell flexibility (rigid shell and hoops), internal air modes (uniform-pressure spring only), vents, \
stick flexure and hand mechanics. Head material, tension, damping ratio, stick profile and Hertz tip \
are editable estimates, not measurements of a specific drum. Only the declared head eigen-window is \
retained; partials above it are absent, not approximated. No samples, noise, normalization, automatic \
gain or pitch control is applied anywhere. The atmosphere is a prescribed, uniform, nonreacting \
ideal gas shared by the cavity and exterior, with the material owner's composition and phase limits; \
no gas evolution, thermal exchange, gas damping or atmospheric absorption is modeled here.";

/// One built instrument and its running voice.
pub struct LiveDrumCore {
    inst: Instrument,
    voice: Voice,
    batter: Vec<f64>,
    resonant: Vec<f64>,
    refusal: Option<Refusal>,
    builds: u64,
}

impl LiveDrumCore {
    /// Build from the reference defaults overlaid with a flat JSON object.
    ///
    /// # Errors
    /// Any typed refusal of parameters, geometry, eigen-solve or owners.
    pub fn new(params_json: &str) -> Result<Self, Refusal> {
        let params = DrumParams::default().overlay_json(params_json)?;
        let inst = Instrument::build(params)?;
        let voice = Voice::new(&inst, inst.params.substeps)?;
        Ok(Self { inst, voice, batter: Vec::new(), resonant: Vec::new(), refusal: None, builds: 1 })
    }

    /// Rebuild from the CURRENT parameters overlaid with a flat JSON object.
    /// Returns true when the heads were re-meshed and re-solved, false when
    /// only numerics/observation changed (the eigen-solve is reused).
    ///
    /// The old instrument is stopped, not faded: the new one starts at rest
    /// and silent. On refusal the old instrument keeps playing untouched.
    ///
    /// # Errors
    /// Any typed refusal; nothing is replaced on failure.
    pub fn rebuild(&mut self, params_json: &str) -> Result<bool, Refusal> {
        let params = self.inst.params.overlay_json(params_json)?;
        params.validate()?;
        if params.same_structure(&self.inst.params) {
            let previous = self.inst.params.clone();
            self.inst.retarget(params);
            match Voice::new(&self.inst, self.inst.params.substeps) {
                Ok(voice) => {
                    self.voice = voice;
                    self.refusal = None;
                    self.builds += 1;
                    Ok(false)
                }
                Err(refusal) => {
                    self.inst.retarget(previous);
                    Err(refusal)
                }
            }
        } else {
            let inst = Instrument::build(params)?;
            let voice = Voice::new(&inst, inst.params.substeps)?;
            self.inst = inst;
            self.voice = voice;
            self.refusal = None;
            self.builds += 1;
            Ok(true)
        }
    }

    /// Same instrument, at rest and silent (clears a halted voice).
    ///
    /// # Errors
    /// Owner refusals while recompiling the rest state.
    pub fn reset(&mut self) -> Result<(), Refusal> {
        self.voice = Voice::new(&self.inst, self.inst.params.substeps)?;
        self.refusal = None;
        Ok(())
    }

    /// The built instrument.
    #[must_use]
    pub const fn instrument(&self) -> &Instrument {
        &self.inst
    }
    /// The running voice.
    #[must_use]
    pub const fn voice(&self) -> &Voice {
        &self.voice
    }

    /// Request a stroke at a head-centred station [m] with a tip speed [m/s].
    ///
    /// # Errors
    /// `strike_outside_head`, `strike_speed`, `halted`.
    pub fn strike(&mut self, x_m: f64, y_m: f64, speed_m_s: f64) -> Result<StrikeOutcome, Refusal> {
        self.voice.strike(&self.inst, [x_m, y_m], speed_m_s)
    }

    /// Render `frames` frames of physical pressure [Pa] at both microphones.
    /// After a refusal the voice is halted and the output is exact silence;
    /// `take_refusal` reports it once.
    pub fn render_pa(&mut self, frames: usize) -> (&[f64], &[f64]) {
        self.batter.resize(frames, 0.0);
        self.resonant.resize(frames, 0.0);
        if let Err(refusal) = self.voice.render_pa(&self.inst, &mut self.batter, &mut self.resonant) {
            self.refusal = Some(refusal);
        }
        (&self.batter, &self.resonant)
    }

    /// Planar float audio `[batter mic (frames) | resonant mic (frames)]`,
    /// each sample = pressure / full_scale_pa. No limiter and no clamp: the
    /// overs are counted in the snapshot, not hidden.
    pub fn render_f32(&mut self, frames: usize) -> Vec<f32> {
        let scale = 1.0 / self.inst.params.full_scale_pa;
        let (batter, resonant) = self.render_pa(frames);
        let mut out = Vec::with_capacity(2 * frames);
        out.extend(batter.iter().map(|p| (p * scale) as f32));
        out.extend(resonant.iter().map(|p| (p * scale) as f32));
        out
    }

    /// The refusal that halted the voice, reported once.
    pub fn take_refusal(&mut self) -> Option<Refusal> {
        self.refusal.take()
    }

    /// Accepted state for drawing.
    pub fn snapshot(&mut self) -> Snapshot {
        self.voice.snapshot(&self.inst)
    }

    /// Transverse displacement of every mesh node [m], positive DOWNWARD for
    /// both heads: `[batter nodes | resonant nodes]`, from the accepted modal
    /// state. Read-only: does not reset the snapshot's pressure peaks.
    #[must_use]
    pub fn node_displacements(&self) -> Vec<f32> {
        let (all_q, _) = self.voice.modal_state();
        let nodes = self.inst.nodes.len();
        let mut out = Vec::with_capacity(2 * nodes);
        let mut offset = 0;
        for head in &self.inst.heads {
            let count = head.mode_count();
            let q = &all_q[offset..offset + count];
            for node in 0..nodes {
                let shapes = &head.node_shapes[node * count..(node + 1) * count];
                out.push(shapes.iter().zip(q).map(|(s, q)| s * q).sum::<f64>() as f32);
            }
            offset += count;
        }
        out
    }

    /// Benchmark a scratch voice of the SAME instrument at another mechanical
    /// rate; the running voice is untouched. The caller measures wall time.
    ///
    /// # Errors
    /// Any refusal while building or rendering the scratch voice.
    pub fn bench(
        &self, substeps: usize, frames: usize, strike_interval: usize, x_m: f64, y_m: f64, speed_m_s: f64,
    ) -> Result<BenchReport, Refusal> {
        Voice::bench(&self.inst, substeps, frames, strike_interval, [x_m, y_m], speed_m_s)
    }

    /// Everything static about the built instrument, as JSON.
    #[must_use]
    pub fn info_json(&self) -> String {
        let inst = &self.inst;
        let p = &inst.params;
        let head = |index: usize| {
            let h = &inst.heads[index];
            format!(
                concat!(
                    "{{\"name\":{},\"modes\":{},\"frequencies_hz\":{},\"film_mass_kg\":{},",
                    "\"modal_areas_m2_per_sqrt_kg\":{},\"damping_per_s\":{},\"max_eigen_residual\":{}}}"
                ),
                json::string(if index == 0 { "batter" } else { "resonant" }),
                h.mode_count(),
                json::numbers(&h.frequencies_hz()),
                json::number(h.film.mass_kg),
                json::numbers(&h.areas),
                json::numbers(&h.damping_per_s),
                json::number(h.max_residual),
            )
        };
        let filter = self.voice.filter_delay_frames();
        let travel = self.voice.travel_delay_frames();
        let gas = p.atmosphere.state();
        format!(
            concat!(
                "{{\"model\":\"reduced real-time model\",\"engine\":\"fs-live-drum {}\",\"fma\":{},",
                "\"sample_rate_hz\":{},\"substeps\":{},\"mechanical_rate_hz\":{},",
                "\"params\":{},",
                "\"atmosphere\":{{\"temperature_k\":{},\"pressure_pa\":{},\"density_kg_m3\":{},\"sound_speed_m_s\":{},\"gamma\":{}}},",
                "\"heads\":[{},{}],",
                "\"coupled_frequencies_hz\":{},",
                "\"mesh\":{{\"nodes\":{},\"triangles\":{},\"radial_intervals\":{},\"azimuths\":{}}},",
                "\"stick\":{{\"full_mass_kg\":{},\"effective_mass_kg\":{},\"hertz_stiffness_n_m15\":{},",
                "\"hertz_exponent\":{},\"contact_loss_s_m\":{},\"launch_gap_m\":{},\"catch_gap_m\":{}}},",
                "\"cavity\":{{\"volume_m3\":{},\"bulk_modulus_pa\":{}}},",
                "\"radiation\":{{\"model\":\"time-domain Rayleigh half-space integral, one infinite rigid baffle per head\",",
                "\"full_scale_pa\":{},\"filter_delay_frames\":{},",
                "\"travel_delay_frames\":[[{},{}],[{},{}]],\"multiply_adds_per_frame\":{}}},",
                "\"snapshot_layout\":{{\"header\":{},\"fields\":{}}},",
                "\"builds\":{},\"provenance\":{},\"no_claim\":{}}}"
            ),
            env!("CARGO_PKG_VERSION"), json::string(crate::fma_path()),
            OUTPUT_RATE_HZ, self.voice.substeps(), self.voice.mechanical_rate_hz(),
            p.to_json(),
            gas.temperature, gas.pressure, gas.density, gas.sound_speed, gas.gamma,
            head(0), head(1),
            json::numbers(&inst.coupled_frequencies_hz),
            inst.nodes.len(), inst.tris.len(), p.radial_intervals, p.azimuths,
            json::number(inst.stick.full_mass_kg), json::number(inst.stick.effective_mass_kg),
            json::number(inst.stick.hertz_stiffness), json::number(HERTZ_ALPHA),
            json::number(p.contact_loss_s_m), json::number(LAUNCH_GAP_M), json::number(p.catch_gap_m),
            json::number(p.volume_m3()), json::number(p.bulk_modulus_pa()),
            json::number(p.full_scale_pa), json::numbers(&filter),
            travel[0].0, travel[0].1, travel[1].0, travel[1].1,
            self.voice.radiation_multiply_adds(),
            SNAPSHOT_HEADER, snapshot_fields_json(),
            self.builds, json::string(PROVENANCE), json::string(NO_CLAIM),
        )
    }
}

/// Fixed scalar slots at the front of the flat snapshot array.
pub const SNAPSHOT_HEADER: usize = 24;
const SNAPSHOT_FIELDS: [&str; SNAPSHOT_HEADER] = [
    "time_s", "mechanical_steps", "stored_j", "dissipated_j", "supplied_j", "recovered_j",
    "balance_residual_j", "max_step_residual_j", "stick_present", "stick_in_contact", "strike_pending",
    "stick_x_m", "stick_y_m", "stick_position_m", "stick_velocity_m_s", "penetration_m",
    "contact_force_n", "cavity_pressure_pa", "strikes_launched", "peak_pa_batter", "peak_pa_resonant",
    "clipped_frames", "batter_modes", "resonant_modes",
];
fn snapshot_fields_json() -> String {
    let mut out = String::from("[");
    for (i, name) in SNAPSHOT_FIELDS.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        out.push_str(&json::string(name));
    }
    out.push(']');
    out
}

/// Flat `Float64Array` image of a snapshot: `SNAPSHOT_HEADER` scalars (see
/// `snapshot_layout.fields` in the info JSON), then q (batter, resonant),
/// then v (batter, resonant).
#[must_use]
pub fn flatten_snapshot(s: &Snapshot, batter_modes: usize, resonant_modes: usize) -> Vec<f64> {
    let flag = |b: bool| if b { 1.0 } else { 0.0 };
    let mut out = Vec::with_capacity(SNAPSHOT_HEADER + s.q.len() + s.v.len());
    out.extend_from_slice(&[
        s.time_s, s.mechanical_steps as f64, s.ledger.stored_j, s.ledger.dissipated_j,
        s.ledger.supplied_j, s.ledger.recovered_j, s.ledger.balance_residual_j,
        s.ledger.max_step_residual_j, flag(s.stick_present), flag(s.stick_in_contact),
        flag(s.strike_pending), s.stick_xy_m[0], s.stick_xy_m[1], s.stick_position_m,
        s.stick_velocity_m_s, s.penetration_m, s.contact_force_n, s.cavity_pressure_pa,
        s.strikes_launched as f64, s.peak_pa[0], s.peak_pa[1], s.clipped_frames as f64,
        batter_modes as f64, resonant_modes as f64,
    ]);
    out.extend_from_slice(&s.q);
    out.extend_from_slice(&s.v);
    out
}
