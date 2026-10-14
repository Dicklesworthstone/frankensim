//! Real-time play: one kit piece on a grid `scale` times coarser than the offline reference,
//! struck anywhere on its playing surface, stepped in blocks at the 48 kHz output rate, and not
//! stepped at all while it is silent. Everything else (heads, shells, air, contacts, microphone)
//! is the same solver as the offline renders.

use crate::kit::{self, CymbalSpec, DrumSpec, Piece, BRONZE, OUT_RATE};
use crate::*;

/// Every piece the live kit knows, in the order a drummer names them.
pub const PIECES: [&str; 9] = ["kick", "snare", "snare-off", "tom-12", "tom-14", "hihat-closed", "hihat-open", "crash", "ride"];

pub struct Live {
    pub ins: Instrument,
    pub name: String,
    /// The struck body (batter head, top hi-hat shell, cymbal) and its radius [m].
    pub target: usize,
    pub radius: f64,
    stick: usize,
    contact: usize,
    /// Striker face radius (a felt beater) or 0 for a stick tip.
    foot: f64,
    asleep: bool,
    quiet_samples: usize,
    ke_peak: f64,
    /// Kinetic energy of every elastic body after the last block [J].
    pub kinetic: f64,
    /// Largest contact force of the last block [N].
    pub force: f64,
    pub strokes: u64,
}

/// Stays awake this long after the motion has died away.
const QUIET_S: f64 = 0.3;

impl Live {
    pub fn new(piece: &str, scale: f64) -> Result<Live, String> {
        if !(scale.is_finite() && (1.0..=12.0).contains(&scale)) {
            return Err(format!("grid scale {scale} outside [1, 12]"));
        }
        let drum = |mut d: DrumSpec| -> Result<(Piece, f64, f64), String> {
            d.grid *= scale;
            d.oversample = 1;
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
        if p.oversample != 1 {
            return Err(format!("{piece}: grid scale {scale} still needs {}x oversampling; use a coarser grid", p.oversample));
        }
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
        // settle (snare-wire preload, hi-hat pedal) before anyone can play
        let settle = (p.settle * OUT_RATE).round() as usize;
        let mut out = [0.0; 1];
        for _ in 0..settle {
            p.instrument.step(&mut out);
        }
        Ok(Live {
            ins: p.instrument, name: piece.into(), target, radius, stick, contact, foot,
            asleep: false, quiet_samples: 0, ke_peak: 0.0, kinetic: 0.0, force: 0.0, strokes: 0,
        })
    }

    /// Strike at (x, y) [m] from the centre of the playing surface with the stick (or beater)
    /// moving at `speed` [m/s]; it lands on the next step.
    pub fn strike(&mut self, x: f64, y: f64, speed: f64) {
        let speed = speed.clamp(0.05, 12.0);
        let r = (x * x + y * y).sqrt();
        let lim = 0.97 * self.radius;
        let (x, y) = if r > lim { (x * lim / r, y * lim / r) } else { (x, y) };
        let fp = self.ins.bodies[self.target].footprint(x, y, self.foot);
        let t = self.ins.time();
        let c = &mut self.ins.contacts[self.contact];
        c.a = fp.into_iter().map(|(k, w)| (self.target, k, w)).collect();
        self.ins.strikes.retain(|s| !s.done);
        self.ins.strikes.push(Strike { time: t, striker: self.stick, contact: self.contact, speed, clear: 2e-4, done: false });
        self.asleep = false;
        self.quiet_samples = 0;
        self.strokes += 1;
    }

    pub fn asleep(&self) -> bool {
        self.asleep
    }

    /// Fill `out` with the next samples at 48 kHz: microphone pressure [Pa] times `gain`.
    /// While asleep nothing is stepped and the block is exact silence.
    pub fn render(&mut self, out: &mut [f32], gain: f64) {
        if self.asleep {
            out.fill(0.0);
            self.kinetic = 0.0;
            self.force = 0.0;
            return;
        }
        let mut buf = [0.0f64; 1];
        let mut peak = 0.0f64;
        let mut force = 0.0f64;
        for o in out.iter_mut() {
            self.ins.step(&mut buf);
            *o = (buf[0] * gain) as f32;
            peak = peak.max(buf[0].abs());
            force = force.max(self.ins.contacts[self.contact].force);
        }
        let k = self.ins.k();
        let mut ke = 0.0;
        for (i, b) in self.ins.bodies.iter().enumerate() {
            if i == self.stick || b.n < 2 {
                continue;
            }
            for ((m, w), wp) in b.mass.iter().zip(&b.w).zip(&b.wp) {
                let v = (w - wp) / k;
                ke += 0.5 * m * v * v;
            }
        }
        self.kinetic = ke;
        self.force = force;
        self.ke_peak = self.ke_peak.max(ke);
        // asleep once the motion is 70 dB below its loudest and the microphone hears nothing
        if ke <= 1e-7 * self.ke_peak && peak < 3e-4 && force == 0.0 {
            self.quiet_samples += out.len();
            if self.quiet_samples as f64 >= QUIET_S * OUT_RATE {
                self.asleep = true;
                self.ke_peak = 0.0;
            }
        } else {
            self.quiet_samples = 0;
        }
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
}
