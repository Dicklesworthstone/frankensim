//! Real-time play: one kit piece on a grid `scale` times coarser than the offline reference,
//! struck anywhere on its playing surface, stepped in blocks at the 48 kHz output rate, and not
//! stepped at all while it is silent. Everything else (heads, shells, air, contacts, microphone)
//! is the same solver as the offline renders.

use crate::kit::{self, CymbalSpec, DrumSpec, Piece, BRONZE, OUT_RATE};
use crate::materials::KitMaterials;
use crate::*;
use fs_material::gas_species::input::{GasRecord, ThermochemicalPhase};
use std::sync::Arc;

/// Every piece the live kit knows, in the order a drummer names them.
pub const PIECES: [&str; 10] = ["kick", "snare", "snare-off", "tom-12", "tom-14", "hihat", "hihat-closed", "hihat-open", "crash", "ride"];

#[derive(Clone)]
struct HiHatPedal {
    rod: usize,
    drive: usize,
    transition_end: f64,
}

/// A clone holds every bit of the piece's state, so it steps exactly as the original would from
/// there on. The browser kit keeps a clone of each grid it builds and times, and plays the chosen
/// grid's clone, so that grid is not settled a second time.
#[derive(Clone)]
pub struct Live {
    pub ins: Instrument,
    pub name: String,
    /// The struck body (batter head, top hi-hat shell, cymbal) and its radius [m].
    pub target: usize,
    pub radius: f64,
    /// (point-mass body, contact) slots shared by successive, non-overlapping impacts.
    strikers: Vec<(usize, usize)>,
    oversample: usize,
    decimator: StreamingDecimator,
    /// Striker face radius (a felt beater) or 0 for a stick tip.
    foot: f64,
    asleep: bool,
    quiet_samples: usize,
    ke_peak: f64,
    sleep_frames: usize,
    sleep_peak: f64,
    sleep_force: f64,
    /// One movable pair of cymbals. Fixed open/closed presets retain their original geometry.
    pedal: Option<HiHatPedal>,
    /// Kinetic energy of every elastic body after the last block [J].
    pub kinetic: f64,
    /// Largest contact force of the last block [N].
    pub force: f64,
    pub strokes: u64,
    // Complete immutable source authority for this trajectory, including
    // canonical species/thermodynamic/transport receipts when supplied.
    atmosphere_source: Option<Arc<GasRecord>>,
}

/// Stays awake this long after the motion has died away.
const QUIET_S: f64 = 0.3;
const SLEEP_FRAMES: usize = 128;
const HIHAT_PREPARE_MIN_S: f64 = 0.5;
const HIHAT_PREPARE_EXTRA_S: f64 = 8.0;

fn gas_record(record: &str) -> Result<GasState, String> {
    let mut fields = record.split_ascii_whitespace();
    let gas = fs_material::gas_species::input::read_gas_state(&mut fields, None)
        .map_err(|e| e.to_string())?;
    if fields.next().is_some() {
        return Err("unexpected fields after the gas record".into());
    }
    Ok(gas)
}

impl Live {
    pub fn new(piece: &str, scale: f64) -> Result<Live, String> {
        Self::with_gas(piece, scale, standard_gas())
    }

    /// One shared hi-hat with a 2 mm rest opening and an explicitly prescribed initial
    /// downward pedal force [N]. A nonzero initial load is settled before playback;
    /// subsequent pedal commands move the same cymbals without resetting their motion.
    /// Preparation renders at least 0.5 s, then allows up to 8 s more for the
    /// ordinary quiet condition. An unsettled or nonfinite initial state refuses.
    pub fn with_hihat_pedal(scale: f64, initial_pedal_n: f64) -> Result<Live, String> {
        if !initial_pedal_n.is_finite() || initial_pedal_n < 0.0 {
            return Err("hi-hat pedal force must be finite and nonnegative [N]".into());
        }
        let mut live = Self::new("hihat", scale)?;
        if initial_pedal_n > 0.0 {
            live.set_hihat_pedal(initial_pedal_n, 0.12)?;
            live.settle_initial_hihat_load()?;
        }
        Ok(live)
    }

    fn settle_initial_hihat_load(&mut self) -> Result<(), String> {
        let minimum = (HIHAT_PREPARE_MIN_S * OUT_RATE).round() as usize;
        let limit = minimum + (HIHAT_PREPARE_EXTRA_S * OUT_RATE).round() as usize;
        let mut rendered = 0;
        let mut scratch = [0.0; SLEEP_FRAMES];
        // Ordinary rendering preserves all physical and causal output state.
        // Only its unchanged quiet gate may stop the motion and clear its filter.
        while rendered < minimum || !self.asleep() {
            if rendered == limit {
                return Err(format!("hi-hat initial load did not become quiescent within {} s",
                    HIHAT_PREPARE_MIN_S + HIHAT_PREPARE_EXTRA_S));
            }
            let boundary = if rendered < minimum { minimum } else { limit };
            let count = scratch.len().min(boundary - rendered);
            self.render(&mut scratch[..count], 1.0);
            if scratch[..count].iter().any(|sample| !sample.is_finite()) || !self.kinetic.is_finite() {
                return Err("hi-hat initial load produced nonfinite state during preparation".into());
            }
            rendered += count;
        }
        Ok(())
    }

    /// A declarative temperature/pressure/humidity or species-mixture record,
    /// using the same parser as native material-driven acoustic scenarios.
    /// This constructs a new physical trajectory; it does not alter a ringing
    /// instrument or erase its stored energy by changing coefficients in place.
    pub fn with_gas_record(piece: &str, scale: f64, record: &str) -> Result<Live, String> {
        Self::with_gas(piece, scale, gas_record(record)?)
    }

    /// A complete source record reaches every original gas-dependent operator
    /// before settling. Copies retain the same immutable source authority.
    pub fn with_resolved_atmosphere(piece: &str, scale: f64, record: GasRecord,
        materials: Option<KitMaterials>) -> Result<Live, String> {
        let mut live = Self::build(piece, scale, *record.state(), materials.as_ref())?;
        live.atmosphere_source = Some(Arc::new(record));
        Ok(live)
    }

    /// Shared canonical source-packet admission for native and browser callers.
    /// Empty `water_species` explicitly means none; phase must be `frozen` or
    /// `liquid-water`. No species, coefficients or phase policy are defaulted.
    pub fn with_gas_source_packet(piece: &str, scale: f64, record: &str,
        bytes: &[u8], lengths: &[u32], water_species: &str, phase: &str,
        materials: Option<KitMaterials>) -> Result<Live, String> {
        let phase = match phase {
            "frozen" => ThermochemicalPhase::FrozenGas,
            "liquid-water" => ThermochemicalPhase::LiquidWaterLimit,
            _ => return Err("gas source phase must be frozen or liquid-water".into()),
        };
        let source = fs_material::gas_species::input::packs::read_source_packet(
            record, bytes, lengths,
            if water_species.is_empty() { None } else { Some(water_species) }, phase,
        ).map_err(|error| error.to_string())?;
        Self::with_resolved_atmosphere(piece, scale, source, materials)
    }

    /// Original thermodynamic models, transport claims and pack identities.
    pub fn atmosphere_source(&self) -> Option<&GasRecord> { self.atmosphere_source.as_deref() }

    /// Number of canonical species source triples retained by this trajectory.
    pub fn gas_source_species(&self) -> usize {
        match self.atmosphere_source() {
            Some(GasRecord::Thermochemical(source)) => source.source_packs().len(),
            _ => 0,
        }
    }

    /// Material-card data and gas state reach the same live solver before
    /// settling; cloning or overlapping strokes retain the complete binding.
    pub fn with_materials_record(piece: &str, scale: f64, record: &str, materials: KitMaterials) -> Result<Live, String> {
        Self::with_materials(piece, scale, gas_record(record)?, materials)
    }

    pub fn with_materials(piece: &str, scale: f64, gas: GasState, materials: KitMaterials) -> Result<Live, String> {
        Self::build(piece, scale, gas, Some(&materials))
    }

    /// Build every gas-dependent loading and observation operator from one
    /// coherent immutable snapshot before settling or playing the instrument.
    pub fn with_gas(piece: &str, scale: f64, gas: GasState) -> Result<Live, String> {
        Self::build(piece, scale, gas, None)
    }

    fn build(piece: &str, scale: f64, gas: GasState, materials: Option<&KitMaterials>) -> Result<Live, String> {
        gas.validate_coherence().map_err(|e| e.to_string())?;
        if !(scale.is_finite() && (1.0..=12.0).contains(&scale)) {
            return Err(format!("grid scale {scale} outside [1, 12]"));
        }
        let drum = |mut d: DrumSpec| -> Result<(Piece, f64, f64), String> {
            d.grid *= scale;
            // Keep the physical preset's internal clock; output remains 48 kHz.
            let foot = if materials.is_some() { 0.0 } else { d.striker.2 };
            let piece = match materials {
                Some(m) => kit::drum_with_materials(&d, &[], gas, m)?,
                None => kit::drum_with_gas(&d, &[], gas)?,
            };
            Ok((piece, d.radius, foot))
        };
        let cym = |mut c: CymbalSpec| -> Result<(Piece, f64, f64), String> {
            c.grid *= scale;
            let piece = match materials {
                Some(m) => kit::cymbal_with_materials(&c, &[], gas, m)?,
                None => kit::cymbal_with_gas(&c, &[], gas)?,
            };
            Ok((piece, c.radius, 0.0))
        };
        let rim = kit::HIHAT_RIM_LAW;
        let hat = kit::hihat_spec().0.radius;
        let hihat = |gap, pedal| match materials {
            Some(m) => kit::hihat_grid_mu_with_materials(gap, pedal, &[], rim, 0.8, scale, kit::RIM_FRICTION, gas, m),
            None => kit::hihat_grid_mu_with_gas(gap, pedal, &[], rim, 0.8, scale, kit::RIM_FRICTION, gas),
        };
        let (mut p, radius, foot) = match piece {
            "kick" => drum(kit::kick())?,
            "snare" => drum(kit::snare())?,
            "snare-off" => drum(kit::snare_wires_off())?,
            "tom-12" => drum(kit::tom(12.0))?,
            "tom-14" => drum(kit::tom(14.0))?,
            "hihat" => (hihat(0.002, 0.0)?, hat, 0.0),
            "hihat-closed" => (hihat(0.0, 80.0)?, hat, 0.0),
            "hihat-open" => (hihat(0.002, 0.0)?, hat, 0.0),
            "crash" => cym(kit::crash())?,
            "ride" => cym(kit::ride())?,
            other => return Err(format!("unknown piece {other}")),
        };
        let (ratio, who) = p.instrument.stability();
        if ratio >= 0.95 {
            return Err(format!("{piece}: leapfrog unstable at scale {scale} ({ratio:.3}, {who})"));
        }
        let ins = &mut p.instrument;
        let stick = ins.bodies.len() - 1; // every kit piece adds its striker last
        let target = 0; // batter head, top hi-hat shell or the cymbal itself
        let contact = if ins.contacts.last().is_some_and(|c| c.b.iter().any(|&(b, _, _)| b == stick)) {
            ins.contacts.len() - 1 // a drum's striker contact (inactive until the first stroke)
        } else {
            // hi-hats and cymbals: one stick contact, moved to wherever the next stroke lands
            let tip = if matches!(piece, "hihat" | "hihat-open" | "hihat-closed") {
                kit::HIHAT_STICK_TIP_RADIUS_M
            } else if piece == "ride" { 0.004 } else if piece == "crash" { 0.008 } else { 0.005 };
            let law = materials.map_or_else(|| kit::hertz_stick(tip, BRONZE.0, BRONZE.1), KitMaterials::impact_law);
            let mut c = Contact::new(vec![(target, 0, 1.0)], vec![(stick, 0, 1.0)], 0.0, law);
            c.active = false;
            ins.contacts.push(c);
            ins.contacts.len() - 1
        };
        // Closed rim contacts frequently reach a rounded Newton fixed point.
        // Open-rim chatter keeps its existing solve policy.
        if piece == "hihat-closed" {
            for rim in &mut ins.contacts[..contact] {
                rim.fixed_point_exit = true;
            }
        }

        // Snare-wire preload retains its existing path. Loaded hi-hats settle
        // below through Live::render, including the microphone/decimator history.
        let settle = if piece == "hihat-closed" { 0 } else {
            (p.settle * p.instrument.rate).round() as usize
        };
        let mut out = [0.0; 1];
        for _ in 0..settle {
            p.instrument.step(&mut out);
        }
        let pedal = if piece == "hihat" {
            let rod = p.instrument.bodies.iter().position(|b| b.name == "rod")
                .ok_or_else(|| "shared hi-hat is missing its pull rod".to_string())?;
            let drive = p.instrument.drives.len();
            p.instrument.drives.push(Drive { body: rod, node: 0, knots: vec![(0.0, 0.0)] });
            Some(HiHatPedal { rod, drive, transition_end: p.instrument.time() })
        } else { None };
        let mut live = Live {
            ins: p.instrument, name: piece.into(), target, radius, strikers: vec![(stick, contact)], foot,
            oversample: p.oversample, decimator: StreamingDecimator::new(p.oversample),
            asleep: false, quiet_samples: 0, ke_peak: 0.0, kinetic: 0.0, force: 0.0, strokes: 0,
            sleep_frames: 0, sleep_peak: 0.0, sleep_force: 0.0,
            pedal,
            atmosphere_source: None,
        };
        if piece == "hihat-closed" {
            live.settle_initial_hihat_load()?;
        }
        Ok(live)
    }

    /// Prescribe downward force on the shared hi-hat's pull rod [N], ramped from its
    /// current force over `transition_seconds`. This is an external mechanical input,
    /// not a stick stroke or a PCM envelope. Zero force releases the return spring.
    /// Fixed `hihat-open`/`hihat-closed` presets refuse because they have different rest gaps.
    pub fn set_hihat_pedal(&mut self, force_n: f64, transition_seconds: f64) -> Result<(), String> {
        if !force_n.is_finite() || force_n < 0.0 {
            return Err("hi-hat pedal force must be finite and nonnegative [N]".into());
        }
        if !transition_seconds.is_finite() || transition_seconds < 0.0 {
            return Err("hi-hat pedal transition must be finite and nonnegative [s]".into());
        }
        let now = self.ins.time();
        let end = now + transition_seconds;
        if !end.is_finite() || (transition_seconds > 0.0 && end <= now) {
            return Err("hi-hat pedal transition is outside the integration clock's range".into());
        }
        let pedal = self.pedal.as_mut().ok_or_else(||
            "pedal control requires the shared hihat instrument".to_string())?;
        let drive = self.ins.drives.get_mut(pedal.drive).filter(|d| d.body == pedal.rod && d.node == 0)
            .ok_or_else(|| "shared hi-hat pedal drive was modified externally".to_string())?;
        let current = drive.at(now);
        if pedal.transition_end <= now && current == -force_n {
            return Ok(());
        }
        drive.knots = if transition_seconds == 0.0 {
            vec![(now, -force_n)]
        } else {
            vec![(now, current), (end, -force_n)]
        };
        pedal.transition_end = end;
        self.asleep = false;
        self.quiet_samples = 0;
        Ok(())
    }

    /// Strike at (x, y) [m] from the centre of the playing surface with the stick (or beater)
    /// moving at `speed` [m/s]; it lands on the next step. Overlapping impacts get separate
    /// strikers on the same surface, including calls made before any samples are rendered.
    pub fn strike(&mut self, x: f64, y: f64, speed: f64) {
        let speed = speed.clamp(0.05, 12.0);
        let r = (x * x + y * y).sqrt();
        let lim = 0.97 * self.radius;
        let (x, y) = if r > lim { (x * lim / r, y * lim / r) } else { (x, y) };
        let fp = self.ins.bodies[self.target].footprint(x, y, self.foot);
        let t = self.ins.time();
        self.ins.strikes.retain(|s| !s.done);
        let (stick, contact) = self.available_striker();
        let c = &mut self.ins.contacts[contact];
        c.a = fp.into_iter().map(|(k, w)| (self.target, k, w)).collect();
        self.ins.strikes.push(Strike { time: t, striker: stick, contact, speed, clear: 2e-4, done: false });
        self.asleep = false;
        self.quiet_samples = 0;
        self.strokes += 1;
    }

    fn available_striker(&mut self) -> (usize, usize) {
        for &(stick, contact) in &self.strikers {
            if self.ins.strikes.iter().any(|s| !s.done && s.striker == stick) {
                continue;
            }
            let c = &self.ins.contacts[contact];
            if !c.active {
                return (stick, contact);
            }
            if c.force != 0.0 {
                continue;
            }
            // Zero force alone also describes a newly launched stick still approaching the
            // surface. Reuse only after separation with no closing relative motion. Read the
            // bodies' current state: later contacts may have changed a cached penetration.
            let mut now = -c.gap;
            let mut previous = -c.gap;
            for (side, sign) in [(&c.a, 1.0), (&c.b, -1.0)] {
                for &(body, node, weight) in side {
                    let b = &self.ins.bodies[body];
                    now += sign * weight * b.w[node];
                    previous += sign * weight * b.wp[node];
                }
            }
            if now <= 0.0 && now <= previous {
                return (stick, contact);
            }
        }

        // Only a non-radiating point mass and its contact are added. Surface indices and
        // microphone taps remain valid; re-preparing the instrument would erase its delay
        // buffers. The contact queue already grows its per-body marks when bodies are added.
        let (template_stick, template_contact) = self.strikers[0];
        let mass = self.ins.bodies[template_stick].mass[0];
        let stick = self.ins.add(Body::point("striker", mass));
        let template = &self.ins.contacts[template_contact];
        let mut c = Contact::new(Vec::new(), vec![(stick, 0, 1.0)], template.gap, template.law);
        c.fixed_point_exit = template.fixed_point_exit;
        c.active = false;
        self.ins.contacts.push(c);
        let contact = self.ins.contacts.len() - 1;
        self.strikers.push((stick, contact));
        (stick, contact)
    }

    pub fn asleep(&self) -> bool {
        self.asleep
    }

    /// Read-only access to the actual instrument, including its source
    /// operators and last solved physical state. Does not advance playback.
    pub fn instrument(&self) -> &Instrument { &self.ins }

    /// Output-filter latency, separate from acoustic travel and device/FIFO latency.
    pub fn filter_latency_frames(&self) -> usize { if self.oversample == 1 { 0 } else { 16 } }

    /// Fill `out` with the next samples at 48 kHz: microphone pressure [Pa] times `gain`.
    /// While asleep nothing is stepped and the block is exact silence.
    pub fn render(&mut self, out: &mut [f32], gain: f64) {
        if out.is_empty() { return; }
        if self.asleep {
            out.fill(0.0);
            self.kinetic = 0.0;
            self.force = 0.0;
            return;
        }
        let mut buf = [0.0f64; 1];
        let mut force = 0.0f64;
        let mut remaining = out;
        while !remaining.is_empty() && !self.asleep {
            let frames = (SLEEP_FRAMES - self.sleep_frames).min(remaining.len());
            let (block,tail) = remaining.split_at_mut(frames);
            remaining = tail;
            let mut peak = 0.0_f64;
            let mut contact_force = 0.0_f64;
            for o in block {
                let mut sample = 0.0;
                for sub in 0..self.oversample {
                    self.ins.step(&mut buf);
                    for &(_,contact) in &self.strikers {
                        contact_force = contact_force.max(self.ins.contacts[contact].force);
                    }
                    if self.oversample == 1 { sample = buf[0]; }
                    else {
                        self.decimator.push(buf[0]);
                        // Offline decimation samples n*factor; the causal image delays it 16 frames.
                        if sub == 0 { sample = self.decimator.value(); }
                    }
                }
                *o = (sample * gain) as f32;
                peak = peak.max(sample.abs());
            }
            force = force.max(contact_force);
            self.sleep_frames += frames;
            self.sleep_peak = self.sleep_peak.max(peak);
            self.sleep_force = self.sleep_force.max(contact_force);
            if self.sleep_frames == SLEEP_FRAMES {
                self.kinetic = self.kinetic_energy();
                self.ke_peak = self.ke_peak.max(self.kinetic);
                let pedal_moving = self.pedal.as_ref().is_some_and(|p| self.ins.time() < p.transition_end);
                if !pedal_moving && self.kinetic <= 1e-7*self.ke_peak && self.sleep_peak < 3e-4 && self.sleep_force == 0.0 {
                    self.quiet_samples += SLEEP_FRAMES;
                    if self.quiet_samples as f64 >= QUIET_S*OUT_RATE {
                        self.asleep = true;
                        self.ke_peak = 0.0;
                        self.decimator.clear();
                    }
                } else { self.quiet_samples = 0; }
                self.sleep_frames = 0;
                self.sleep_peak = 0.0;
                self.sleep_force = 0.0;
            }
        }
        if self.asleep { remaining.fill(0.0); self.kinetic = 0.0; }
        else if self.sleep_frames != 0 { self.kinetic = self.kinetic_energy(); }
        self.force = force;
    }

    fn kinetic_energy(&self) -> f64 {
        let k = self.ins.k();
        let mut ke = 0.0;
        for (index, b) in self.ins.bodies.iter().enumerate() {
            if b.n < 2 && !self.pedal.as_ref().is_some_and(|p| p.rod == index) {
                continue;
            }
            for ((m, w), wp) in b.mass.iter().zip(&b.w).zip(&b.wp) {
                let v = (w - wp) / k;
                ke += 0.5 * m * v * v;
            }
        }
        ke
    }

    /// Node positions of the struck surface (x, y pairs, metres from its centre).
    pub fn surface_xy(&self) -> Vec<f32> {
        self.ins.bodies[self.target].xy.iter().flat_map(|p| [p[0] as f32, p[1] as f32]).collect()
    }

    /// Grid spacing of the struck surface [m].
    pub fn surface_h(&self) -> f64 {
        self.ins.bodies[self.target].h
    }

    /// Transverse velocity of every node of the struck surface [m/s].
    pub fn surface_velocity(&self) -> Vec<f32> {
        let b = &self.ins.bodies[self.target];
        let ik = 1.0 / self.ins.k();
        b.w.iter().zip(&b.wp).map(|(a, c)| ((a - c) * ik) as f32).collect()
    }

    /// Transverse displacement of every node of the struck surface [m].
    pub fn surface_displacement(&self) -> Vec<f32> {
        self.ins.bodies[self.target].w.iter().map(|&a| a as f32).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_live_hihat_uses_the_same_area_scaled_rim_law_as_offline() {
        let coefficients = |law| {
            let Law::Power { k, alpha, mu } = law;
            (k, alpha, mu)
        };
        assert_eq!(coefficients(kit::HIHAT_RIM_LAW), (2.0e9, 1.5, 0.08));
        for scale in [4.0, 7.0] {
            let offline = kit::hihat_scaled(0.002, 0.0, &[(0.0, 2.5)], scale).unwrap();
            let expected = (2.0e9 * scale * scale, 1.5, 0.08);
            assert!(!offline.instrument.contacts.is_empty());
            let strike = offline.instrument.contacts.len() - 1;
            let expected_striker = coefficients(offline.instrument.contacts[strike].law);
            assert_eq!(expected_striker, coefficients(kit::hertz_stick(0.004, BRONZE.0, BRONZE.1)),
                "the source hi-hat uses a 4 mm tip independently of live/native agreement");
            let legacy_striker = coefficients(kit::hertz_stick(0.005, BRONZE.0, BRONZE.1));
            assert_ne!(expected_striker, legacy_striker, "the old generic tip must not satisfy native hi-hat parity");
            for contact in &offline.instrument.contacts[..strike] {
                assert_eq!(coefficients(contact.law), expected);
            }
            for name in ["hihat", "hihat-open", "hihat-closed"] {
                let live = Live::new(name, scale).unwrap();
                let rims = &live.ins.contacts[..live.strikers[0].1];
                assert_eq!(rims.len(), strike, "{name}@{scale}");
                for contact in rims {
                    assert_eq!(coefficients(contact.law), expected, "{name}@{scale}");
                }
                assert_eq!(coefficients(live.ins.contacts[live.strikers[0].1].law), expected_striker,
                    "native/live striker contact mismatch for {name}@{scale}");
            }
        }
    }

    #[test]
    fn uncontrolled_shared_hihat_matches_fixed_open_pressure_and_clock() {
        let mut shared = Live::with_hihat_pedal(7.0, 0.0).unwrap();
        let mut fixed = Live::new("hihat-open", 7.0).unwrap();
        assert_eq!(shared.ins.rate, fixed.ins.rate);
        assert_eq!(shared.ins.step_index, fixed.ins.step_index);
        for live in [&mut shared, &mut fixed] {
            live.strike(0.85 * live.radius, 0.0, 2.5);
        }
        let (mut a, mut b) = (vec![0.0; 8192], vec![0.0; 8192]);
        shared.render(&mut a, 1.0);
        fixed.render(&mut b, 1.0);
        assert!(a.iter().any(|v| *v != 0.0));
        assert_eq!(a.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
            b.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
            "an unused pedal must not change contact policy, clock or pressure");
    }

    #[test]
    fn shared_hihat_pedal_changes_force_without_resetting_ringing_state() {
        let mut live = Live::with_hihat_pedal(7.0, 0.0).unwrap();
        live.strike(0.85 * live.radius, 0.0, 2.5);
        live.render(&mut [0.0; 1024], 1.0);
        let state: Vec<_> = live.ins.bodies.iter().map(|b| (b.w.clone(), b.wp.clone())).collect();
        let rings: Vec<_> = live.ins.mics.iter().map(|m| m.ring.clone()).collect();
        let (step, strokes, bodies, contacts) = (live.ins.step_index, live.strokes,
            live.ins.bodies.len(), live.ins.contacts.len());
        let now = live.ins.time();
        live.set_hihat_pedal(80.0, 0.04).unwrap();
        let pedal = live.pedal.as_ref().unwrap();
        let drive = &live.ins.drives[pedal.drive];
        assert_eq!(drive.at(now), 0.0);
        assert!((drive.at(now + 0.02) + 40.0).abs() < 1e-10);
        assert_eq!(drive.at(now + 0.04), -80.0);
        assert_eq!(live.ins.step_index, step);
        assert_eq!(live.strokes, strokes, "pedal motion must not invent a stick stroke");
        assert_eq!(live.ins.bodies.len(), bodies);
        assert_eq!(live.ins.contacts.len(), contacts);
        assert_eq!(live.ins.bodies.iter().map(|b| (b.w.clone(), b.wp.clone())).collect::<Vec<_>>(), state);
        assert_eq!(live.ins.mics.iter().map(|m| m.ring.clone()).collect::<Vec<_>>(), rings);

        let mut control = live.clone();
        control.set_hihat_pedal(0.0, 0.0).unwrap();
        let mut copy = live.clone();
        let (mut a, mut b, mut unclosed) = (vec![0.0; 4096], vec![0.0; 4096], vec![0.0; 4096]);
        live.render(&mut a, 1.0);
        let mut at = 0;
        for count in [1, 127, 511, 2048, 1409] {
            copy.render(&mut b[at..at + count], 1.0);
            at += count;
        }
        control.render(&mut unclosed, 1.0);
        assert!(a.iter().all(|x| x.is_finite()));
        assert_eq!(a, b, "copies must retain force ramps across arbitrary output blocks");
        assert_ne!(a, unclosed, "a closure must affect the same ringing cymbals");
        let rod = live.pedal.as_ref().unwrap().rod;
        assert_ne!(live.ins.bodies[rod].w, control.ins.bodies[rod].w);
    }

    #[test]
    fn shared_hihat_initial_load_is_settled_on_one_rest_geometry() {
        let open = Live::with_hihat_pedal(7.0, 0.0).unwrap();
        let mut closed = Live::with_hihat_pedal(7.0, 80.0).unwrap();
        assert_eq!(open.ins.bodies.len(), closed.ins.bodies.len());
        assert_eq!(open.ins.contacts.len(), closed.ins.contacts.len());
        for (a, b) in open.ins.bodies.iter().zip(&closed.ins.bodies) {
            assert_eq!(a.origin, b.origin, "the pedal must move bodies, not change their rest gap");
            assert_eq!(a.rest_z, b.rest_z);
        }
        assert!(closed.ins.contacts.iter().any(|c| c.active && c.force > 0.0),
            "a closed initial load must actually reach the rims");
        assert_eq!(closed.strokes, 0);
        let pedal = closed.pedal.as_ref().unwrap();
        assert_eq!(closed.ins.drives[pedal.drive].at(closed.ins.time()), -80.0);
        assert!(closed.ins.time() >= 0.5);
        assert!(closed.asleep(), "held initial state must pass the ordinary quiet gate");

        let mut manual = Live::new("hihat", 7.0).unwrap();
        manual.set_hihat_pedal(80.0, 0.12).unwrap();
        let minimum = (0.5 * OUT_RATE).round() as usize;
        let limit = ((0.5 + 8.0) * OUT_RATE).round() as usize;
        let mut rendered = 0;
        let mut scratch = [0.0; 128];
        while rendered < minimum || !manual.asleep() {
            assert!(rendered < limit, "manual physical preparation exhausted its budget");
            let boundary = if rendered < minimum { minimum } else { limit };
            let count = scratch.len().min(boundary - rendered);
            manual.render(&mut scratch[..count], 1.0);
            assert!(scratch[..count].iter().all(|sample| sample.is_finite()));
            rendered += count;
        }
        assert_eq!(manual.ins.step_index, closed.ins.step_index);
        assert_eq!(manual.ins.mics.iter().map(|mic| mic.ring.clone()).collect::<Vec<_>>(),
            closed.ins.mics.iter().map(|mic| mic.ring.clone()).collect::<Vec<_>>());
        // Compare a fresh physical stroke, not only the now-quiet zero output.
        for live in [&mut manual, &mut closed] { live.strike(0.85 * live.radius, 0.0, 2.5); }
        let (mut a, mut b) = (vec![0.0; 1024], vec![0.0; 1024]);
        manual.render(&mut a, 1.0);
        closed.render(&mut b, 1.0);
        assert!(a.iter().any(|sample| *sample != 0.0));
        assert_eq!(a.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
            b.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
            "initial settling must retain the same causal output history as manual playback");
    }

    #[test]
    fn loaded_hihat_initial_states_are_quiet_with_their_force_and_contact_preserved() {
        for mut live in [Live::new("hihat-closed", 7.0).unwrap(),
            Live::with_hihat_pedal(7.0, 80.0).unwrap()] {
            assert!(live.asleep());
            assert_eq!(live.strokes, 0);
            assert!(live.ins.drives.iter().any(|drive| drive.at(live.ins.time()) == -80.0));
            assert!(live.ins.contacts.iter().any(|contact| contact.active && contact.force > 0.0),
                "quiescence must retain an actual loaded rim contact");
            let state = live.ins.bodies.iter().map(|body| (body.w.clone(), body.wp.clone())).collect::<Vec<_>>();
            let step = live.ins.step_index;
            let mut output = [1.0; 1024];
            live.render(&mut output, 1.0);
            assert!(output.iter().all(|sample| *sample == 0.0));
            assert_eq!(live.ins.step_index, step);
            assert_eq!(live.ins.bodies.iter().map(|body| (body.w.clone(), body.wp.clone())).collect::<Vec<_>>(), state);
            live.strike(0.85 * live.radius, 0.0, 2.5);
            assert!(!live.asleep());
            live.render(&mut output, 1.0);
            assert!(output.iter().all(|sample| sample.is_finite()));
            assert!(output.iter().any(|sample| *sample != 0.0));
        }
        let mut unfinished = Live::new("hihat", 7.0).unwrap();
        unfinished.set_hihat_pedal(80.0, 20.0).unwrap();
        let error = unfinished.settle_initial_hihat_load().unwrap_err();
        assert!(error.contains("did not become quiescent within 8.5 s"));
        assert!(!unfinished.asleep(), "an unfinished external drive cannot be treated as quiet");
        assert!((unfinished.ins.time() - 8.5).abs() < 1e-12);
    }

    #[test]
    fn shared_hihat_pedal_interruptions_are_force_continuous_and_held_commands_are_noops() {
        let mut live = Live::with_hihat_pedal(7.0, 0.0).unwrap();
        live.set_hihat_pedal(80.0, 0.01).unwrap();
        live.render(&mut [0.0; 240], 1.0);
        let now = live.ins.time();
        let pedal = live.pedal.as_ref().unwrap();
        let before = live.ins.drives[pedal.drive].at(now);
        assert!(before < 0.0 && before > -80.0, "the control must interrupt a real ramp");
        live.set_hihat_pedal(0.0, 0.02).unwrap();
        let pedal = live.pedal.as_ref().unwrap();
        assert_eq!(live.ins.drives[pedal.drive].at(now), before);
        assert_eq!(live.ins.drives[pedal.drive].at(now + 0.02), 0.0);
        live.render(&mut [0.0; 1024], 1.0);
        let pedal = live.pedal.as_ref().unwrap();
        let held = live.ins.drives[pedal.drive].knots.clone();
        let mut copy = live.clone();
        live.set_hihat_pedal(0.0, 0.5).unwrap();
        assert_eq!(live.ins.drives[live.pedal.as_ref().unwrap().drive].knots, held);
        let (mut a, mut b) = ([0.0; 256], [0.0; 256]);
        live.render(&mut a, 1.0);
        copy.render(&mut b, 1.0);
        assert_eq!(a, b, "repeating a held force must not restart its trajectory");
    }

    #[test]
    fn shared_hihat_pedal_wakes_sleep_and_retains_unfinished_ramps() {
        let mut live = Live::with_hihat_pedal(7.0, 0.0).unwrap();
        // An exactly resting state is sufficient for the native sleep decision.
        live.asleep = true;
        live.quiet_samples = 48_000;
        live.set_hihat_pedal(80.0, 1.0).unwrap();
        assert!(!live.asleep());
        assert_eq!(live.quiet_samples, 0);
        // Interrupt before the clock advances: this now has a future ramp endpoint
        // but identically zero force and motion. Only the pending-program guard
        // can keep the engine awake through the ordinary 300 ms quiet threshold.
        live.set_hihat_pedal(0.0, 1.0).unwrap();
        let step = live.ins.step_index;
        let mut out = vec![0.0; 16_384];
        live.render(&mut out, 1.0);
        assert!(out.iter().all(|v| *v == 0.0));
        assert!(live.ins.step_index > step);
        assert!(!live.asleep(), "an unfinished external force program cannot be skipped as silence");
        assert_eq!(live.strokes, 0);
        live.set_hihat_pedal(0.0, 0.0).unwrap();
        let pedal = live.pedal.as_ref().unwrap();
        assert_eq!(live.ins.drives[pedal.drive].knots, vec![(live.ins.time(), -0.0)]);
        assert!(live.kinetic_energy().is_finite());
    }

    #[test]
    fn shared_hihat_pedal_refuses_invalid_or_wrong_instrument_commands_atomically() {
        for initial in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(Live::with_hihat_pedal(7.0, initial).is_err());
        }
        let mut live = Live::with_hihat_pedal(7.0, 0.0).unwrap();
        let pedal = live.pedal.as_ref().unwrap();
        let original = live.ins.drives[pedal.drive].knots.clone();
        for (force, ramp) in [(-1.0, 0.0), (f64::NAN, 0.0), (f64::INFINITY, 0.0),
            (80.0, -1.0), (80.0, f64::NAN), (80.0, f64::INFINITY)] {
            assert!(live.set_hihat_pedal(force, ramp).is_err());
            assert_eq!(live.ins.drives[live.pedal.as_ref().unwrap().drive].knots, original);
        }
        for name in ["kick", "hihat-open", "hihat-closed"] {
            assert!(Live::new(name, 7.0).unwrap().set_hihat_pedal(80.0, 0.01).is_err(),
                "{name} must not accept control for a different physical rest geometry");
        }
    }

    /// G0/G3: two MIDI strikes before the next sample must not become just the last strike.
    #[test]
    fn simultaneous_impacts_survive_in_both_call_orders() {
        for (name, scale) in [("kick", 8.0), ("snare-off", 8.0), ("hihat-closed", 7.0)] {
            let template = Live::new(name, scale).unwrap();
            let a = (0.7 * template.radius, 0.1 * template.radius, 1.2);
            let b = (0.3 * template.radius, 0.0, 3.0);
            let render = |hits: &[(f64, f64, f64)]| {
                let mut p = template.clone();
                for &(x, y, speed) in hits {
                    p.strike(x, y, speed);
                }
                let mut out = vec![0.0f32; 4096];
                p.render(&mut out, 1.0);
                assert!(out.iter().all(|x| x.is_finite()), "{name}: non-finite impact");
                assert!(out.iter().any(|x| *x != 0.0), "{name}: inaudible control");
                out
            };
            let (single_a, single_b) = (render(&[a]), render(&[b]));
            assert_ne!(single_a, single_b, "{name}: controls must distinguish the impacts");
            assert_ne!(render(&[a, b]), single_b, "{name}: first impact was overwritten");
            assert_ne!(render(&[b, a]), single_a, "{name}: reversed first impact was overwritten");
        }
    }

    #[test]
    fn a_striker_still_approaching_its_first_contact_is_not_relaunched() {
        let mut p = Live::new("kick", 8.0).unwrap();
        p.strike(0.7 * p.radius, 0.1 * p.radius, 1.2);
        let (stick, contact) = p.strikers[0];
        p.render(&mut [0.0f32; 1], 1.0);
        assert_eq!(p.ins.contacts[contact].force, 0.0, "the first stick is still in flight");
        let old_position = p.ins.contacts[contact].a.clone();
        p.strike(0.3 * p.radius, 0.0, 3.0);
        assert_ne!(p.ins.strikes.last().unwrap().striker, stick);
        assert_eq!(p.ins.contacts[contact].a, old_position);
        assert_eq!(p.strikers.len(), 2);
        let mut copy = p.clone();
        let (mut a, mut b) = (vec![0.0f32; 4096], vec![0.0f32; 4096]);
        p.render(&mut a, 1.0);
        copy.render(&mut b, 1.0);
        assert!(a.iter().all(|x| x.is_finite()));
        assert_eq!(a, b, "an overlapping-impact clone must retain every striker");
    }

    #[test]
    fn adding_an_overlapping_striker_preserves_the_ringing_surface_and_microphones() {
        let mut p = Live::new("kick", 8.0).unwrap();
        p.strike(0.3 * p.radius, 0.0, 3.0);
        p.render(&mut [0.0f32; 256], 1.0);
        let rings = p.ins.mics.iter().map(|m| m.ring.clone()).collect::<Vec<_>>();
        assert!(rings.iter().flatten().any(|x| *x != 0.0), "the control must retain delayed sound");
        let surface = p.ins.bodies[p.target].w.clone();
        let history = p.ins.bodies[p.target].wp.clone();
        let step = p.ins.step_index;
        let bodies = p.ins.bodies.len();
        for _ in 0..3 {
            p.strike(0.6 * p.radius, 0.0, 2.0);
        }
        assert!(p.ins.bodies.len() > bodies);
        assert_eq!(p.ins.step_index, step, "striking cannot retime the score");
        assert_eq!(p.ins.bodies[p.target].w, surface);
        assert_eq!(p.ins.bodies[p.target].wp, history);
        assert_eq!(p.ins.mics.iter().map(|m| m.ring.clone()).collect::<Vec<_>>(), rings);
        assert!(p.ins.bodies[bodies..].iter().all(|b| b.n == 1 && !b.radiates));
        let mut out = [0.0f32; 512];
        p.render(&mut out, 1.0);
        assert!(out.iter().all(|x| x.is_finite()));
    }

    #[test]
    fn released_slots_are_reused_without_stealing_a_pending_impact() {
        let mut p = Live::new("kick", 8.0).unwrap();
        for _ in 0..2 {
            p.strike(0.3 * p.radius, 0.0, 3.0);
        }
        p.render(&mut [0.0f32; 1], 1.0);
        // Put both point masses clear of the surface and moving away: no contact load or
        // pending launch remains, so subsequent impacts can reuse these physical slots.
        for &(stick, contact) in &p.strikers {
            p.ins.bodies[stick].w[0] = 1.0;
            p.ins.bodies[stick].wp[0] = 0.99;
            p.ins.contacts[contact].force = 0.0;
        }
        let bodies = p.ins.bodies.len();
        p.strike(0.2 * p.radius, 0.0, 1.0);
        p.strike(0.5 * p.radius, 0.0, 2.0);
        assert_eq!(p.ins.bodies.len(), bodies, "free slots should not accumulate");
        let pending = p.ins.strikes.iter().filter(|s| !s.done).collect::<Vec<_>>();
        assert_eq!(pending.len(), 2);
        assert_ne!(pending[0].striker, pending[1].striker);
        assert_ne!(pending[0].contact, pending[1].contact);
    }

    #[test]
    fn every_piece_builds_at_the_live_grid_strikes_and_falls_asleep() {
        for name in PIECES {
            let mut p = Live::new(name, 5.0).unwrap_or_else(|e| panic!("{name}: {e}"));
            p.strike(0.3 * p.radius, 0.0, 3.0);
            let mut blk = vec![0.0f32; 480];
            let mut peak = 0.0f32;
            for _ in 0..40 {
                p.render(&mut blk, 1.0);
                peak = blk.iter().fold(peak, |m, v| m.max(v.abs()));
                assert!(blk.iter().all(|v| v.is_finite()), "{name}: non-finite output");
            }
            assert!(peak > 1e-3, "{name}: stroke inaudible (peak {peak} Pa)");
            // long enough for the kick and the snare to die away and stop being stepped
            if matches!(name, "kick" | "snare" | "snare-off") {
                for _ in 0..800 {
                    p.render(&mut blk, 1.0);
                }
                assert!(p.asleep(), "{name}: still awake after 8.4 s (kinetic {:.3e} J)", p.kinetic);
                p.strike(0.0, 0.0, 2.0);
                assert!(!p.asleep());
            }
        }
    }

    /// `Live` exists only after its constructor has stepped the whole settle, so no stroke can reach a
    /// snare whose wires are still coming on, in the piece or in a clone of it.
    #[test]
    fn no_stroke_can_land_before_the_snare_wires_are_on() {
        let live = Live::new("snare", 8.0).unwrap();
        assert!(!live.ins.load_is_coming_on());
        assert!(!live.clone().ins.load_is_coming_on());
    }

    /// A clone, taken settled or mid-ring, renders bit for bit what the original goes on to render.
    #[test]
    fn a_clone_renders_what_the_original_would() {
        for name in PIECES {
            let mut a = Live::new(name, 8.0).unwrap_or_else(|e| panic!("{name}: {e}"));
            let mut b = a.clone();
            let (mut x, mut y) = (vec![0.0f32; 128], vec![0.0f32; 128]);
            for k in 0..240 {
                if k % 60 == 0 {
                    // the second pair is cloned mid-ring, after two strokes
                    if k == 120 {
                        b = a.clone();
                    }
                    a.strike(0.3 * a.radius, 0.0, 3.0);
                    b.strike(0.3 * b.radius, 0.0, 3.0);
                }
                a.render(&mut x, 1.0);
                b.render(&mut y, 1.0);
                let same = x.iter().zip(&y).all(|(p, q)| p.to_bits() == q.to_bits());
                assert!(same && a.kinetic.to_bits() == b.kinetic.to_bits() && a.force.to_bits() == b.force.to_bits(), "{name}: block {k} differs");
            }
        }
    }
}
