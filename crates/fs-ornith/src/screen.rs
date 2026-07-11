//! Stage 2: SCREEN WIDE. The BEM panel L/D (with a DOCUMENTED drag
//! proxy — the inviscid panel method has no drag of its own: profile
//! drag cd₀(t) plus induced drag cl²/(π·AR·e), the classic smoke
//! model) plus the fs-vpm FLAPPING metric (wake vortices shed by the
//! fs-bem unsteady wake, advected as vortex particles; thrust proxy =
//! streamwise drift of the vorticity field per unit time). Generations
//! are E-RACED through fs-race: dominated candidates die at a fraction
//! of their solve budget with anytime validity, and the payoff is
//! MEASURED (evaluations used vs the fixed-N equivalent), ledgered in
//! the report.

use crate::param::OrnithCandidate;
use fs_bem::panel2d::solve;
use fs_bem::wake2d::WakeSim;
use fs_vpm::{VortexParticle, advect};

/// Screen-stage panel count (smoke fixture).
pub const PANELS: usize = 64;
/// Effective aspect ratio of the smoke wing (drag proxy).
const ASPECT: f64 = 6.0;
/// Oswald efficiency (drag proxy).
const OSWALD: f64 = 0.9;

/// The screening L/D of one candidate (higher is better).
#[must_use]
pub fn lift_to_drag(c: &OrnithCandidate) -> f64 {
    let sol = solve(&c.section(PANELS), c.alpha)
        .expect("ornith candidate must satisfy the bounded panel-solve contract");
    let cl = sol.cl;
    // Profile drag grows with thickness; induced drag with cl².
    let cd0 = 0.006 + 0.06 * c.thickness * c.thickness / 0.0144;
    let cd = cd0 + cl * cl / (std::f64::consts::PI * ASPECT * OSWALD);
    cl / cd
}

/// The flapping-gait maneuver proxy: shed the unsteady BEM wake, hand
/// the vortices to fs-vpm as particles, advect, and measure the
/// streamwise vorticity drift (a thrust-like signature; higher =
/// stronger gait authority).
#[must_use]
pub fn flap_metric(c: &OrnithCandidate) -> f64 {
    let foil = c.section(PANELS);
    let mut sim = WakeSim::new(&foil, c.alpha, 0.08, 0.05)
        .expect("ornith candidate must satisfy the bounded wake contract");
    for _ in 0..24 {
        sim.step()
            .expect("short bounded wake fixture must remain finite");
    }
    // Lift the wake into vortex particles, modulate by the gait.
    let particles: Vec<VortexParticle> = sim
        .wake()
        .iter()
        .enumerate()
        .map(|(k, w)| {
            let phase = c.flap_freq * k as f64;
            let gain = c.flap_amp.mul_add(fs_math::det::sin(phase), 1.0);
            VortexParticle::new(w.pos, w.gamma * gain)
        })
        .collect();
    let x0: f64 =
        particles.iter().map(|p| p.pos[0].abs()).sum::<f64>() / particles.len().max(1) as f64;
    let moved = advect_steps(&particles, 0.1, 12, 0.05);
    let x1: f64 = moved.iter().map(|p| p.pos[0].abs()).sum::<f64>() / moved.len().max(1) as f64;
    (x1 - x0) / 1.2
}

/// Advect helper mirroring fs-vpm's `simulate` with a fixed step count.
fn advect_steps(
    particles: &[VortexParticle],
    dt: f64,
    steps: usize,
    core: f64,
) -> Vec<VortexParticle> {
    let mut state = particles.to_vec();
    for _ in 0..steps {
        state = advect(&state, dt, core);
    }
    state
}

/// The e-raced generation report — the P7 evidence row.
#[derive(Debug, Clone)]
pub struct ScreenReport {
    /// Winner index into the generation.
    pub winner: usize,
    /// Candidates eliminated early.
    pub eliminated: usize,
    /// Race evaluations actually spent.
    pub evaluations_used: u64,
    /// What a fixed-N tournament would have spent.
    pub fixed_n_equivalent: u64,
    /// Deterministic screening scores (−L/D, the race loss).
    pub losses: Vec<f64>,
}

/// Race one generation on the (noisy) screening loss. Losses are
/// normalized onto a declared analytical support per the fs-eproc
/// PairwiseRace contract (the vessel flagship's measured lesson: unscaled
/// ~1e-3 gaps starve the betting e-process).
///
/// # Errors
/// Propagates a structured [`fs_race::RaceError`] if the analytically
/// bounded normalization or an observation violates the e-process input
/// contract.
///
/// # Panics
/// If fewer than two candidates are supplied.
pub fn screen_generation(
    candidates: &[OrnithCandidate],
    seed: u64,
) -> Result<ScreenReport, fs_race::RaceError> {
    assert!(candidates.len() >= 2, "a generation needs candidates");
    let base: Vec<f64> = candidates.iter().map(|c| -lift_to_drag(c)).collect();
    let spread = base.iter().fold(f64::NEG_INFINITY, |m, &v| m.max(v))
        - base.iter().fold(f64::INFINITY, |m, &v| m.min(v));
    let scale = 1.5 / spread.max(1e-9);
    let kills = fs_exec::KillRegistry::new();
    for candidate in 0..candidates.len() {
        let _ = kills.register(candidate as u64);
    }
    let mut loss = |i: usize, t: u64| {
        let mut h = (i as u64) << 32 ^ t ^ seed;
        h ^= h << 13;
        h ^= h >> 7;
        h ^= h << 17;
        #[allow(clippy::cast_precision_loss)]
        let jitter = ((h >> 11) as f64 / (1u64 << 53) as f64 - 0.5) * 0.02;
        (base[i] - base.iter().fold(f64::INFINITY, |m, &v| m.min(v))).mul_add(scale, jitter)
    };
    let out = fs_race::race_field(
        &mut loss,
        candidates.len(),
        // Normalized base lies in [0, 1.5]; jitter has total width 0.02.
        fs_race::RaceSettings::new(fs_race::LossSpan::new(1.52).expect("positive constant")),
        &kills,
    )?;
    Ok(ScreenReport {
        winner: out.winner,
        eliminated: out.eliminated.len(),
        evaluations_used: out.evaluations_used,
        fixed_n_equivalent: out.fixed_n_equivalent,
        losses: base,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn advect_steps_matches_manual() {
        let p = vec![
            VortexParticle::new([0.0, 0.1], 1.0),
            VortexParticle::new([0.0, -0.1], -1.0),
        ];
        let a = advect_steps(&p, 0.05, 3, 0.05);
        let mut b = p;
        for _ in 0..3 {
            b = advect(&b, 0.05, 0.05);
        }
        assert_eq!(a[0].pos[0].to_bits(), b[0].pos[0].to_bits());
    }
}
