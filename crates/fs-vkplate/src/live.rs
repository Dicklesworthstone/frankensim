//! Real-time play: one kit piece on a grid `scale` times coarser than the offline reference,
//! struck anywhere on its playing surface, stepped in blocks at the 48 kHz output rate, and not
//! stepped at all while it is silent. Everything else (heads, shells, air, contacts, microphone)
//! is the same solver as the offline renders.

use crate::kit::{self, CymbalSpec, DrumSpec, Piece, BRONZE, OUT_RATE};
use crate::*;

/// Every piece the live kit knows, in the order a drummer names them.
pub const PIECES: [&str; 9] = ["kick", "snare", "snare-off", "tom-12", "tom-14", "hihat-closed", "hihat-open", "crash", "ride"];

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
    /// Kinetic energy of every elastic body after the last block [J].
    pub kinetic: f64,
    /// Largest contact force of the last block [N].
    pub force: f64,
    pub strokes: u64,
}

/// Stays awake this long after the motion has died away.
const QUIET_S: f64 = 0.3;
const SLEEP_FRAMES: usize = 128;

impl Live {
    pub fn new(piece: &str, scale: f64) -> Result<Live, String> {
        if !(scale.is_finite() && (1.0..=12.0).contains(&scale)) {
            return Err(format!("grid scale {scale} outside [1, 12]"));
        }
        let drum = |mut d: DrumSpec| -> Result<(Piece, f64, f64), String> {
            d.grid *= scale;
            // Keep the physical preset's internal clock; output remains 48 kHz.
            let foot = d.striker.2;
            Ok((kit::drum(&d, &[])?, d.radius, foot))
        };
        let cym = |mut c: CymbalSpec| -> Result<(Piece, f64, f64), String> {
            c.grid *= scale;
            Ok((kit::cymbal(&c, &[])?, c.radius, 0.0))
        };
        let rim = Law::Power { k: 2.0e9, alpha: 1.5, mu: 3.0 };
        let hat = kit::hihat_spec().0.radius;
        let (mut p, radius, foot) = match piece {
            "kick" => drum(kit::kick())?,
            "snare" => drum(kit::snare())?,
            "snare-off" => drum(kit::snare_wires_off())?,
            "tom-12" => drum(kit::tom(12.0))?,
            "tom-14" => drum(kit::tom(14.0))?,
            "hihat-closed" => (kit::hihat_grid(0.0, 80.0, &[], rim, 0.8, scale)?, hat, 0.0),
            "hihat-open" => (kit::hihat_grid(0.002, 0.0, &[], rim, 0.8, scale)?, hat, 0.0),
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
            let tip = if piece == "ride" { 0.004 } else if piece == "crash" { 0.008 } else { 0.005 };
            let mut c = Contact::new(vec![(target, 0, 1.0)], vec![(stick, 0, 1.0)], 0.0, kit::hertz_stick(tip, BRONZE.0, BRONZE.1));
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

        // settle (snare-wire preload, hi-hat pedal) before anyone can play
        let settle = (p.settle * p.instrument.rate).round() as usize;
        let mut out = [0.0; 1];
        for _ in 0..settle {
            p.instrument.step(&mut out);
        }
        Ok(Live {
            ins: p.instrument, name: piece.into(), target, radius, strikers: vec![(stick, contact)], foot,
            oversample: p.oversample, decimator: StreamingDecimator::new(p.oversample),
            asleep: false, quiet_samples: 0, ke_peak: 0.0, kinetic: 0.0, force: 0.0, strokes: 0,
            sleep_frames: 0, sleep_peak: 0.0, sleep_force: 0.0,
        })
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
                if self.kinetic <= 1e-7*self.ke_peak && self.sleep_peak < 3e-4 && self.sleep_force == 0.0 {
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
        for b in &self.ins.bodies {
            if b.n < 2 {
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
