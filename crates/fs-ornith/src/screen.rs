//! Stage 2: SCREEN WIDE. The BEM panel L/D (with a DOCUMENTED drag
//! proxy — the inviscid panel method has no drag of its own: profile
//! drag cd₀(t) plus induced drag cl²/(π·AR·e), the classic smoke
//! model) plus the fs-vpm FLAPPING metric (wake vortices shed by the
//! fs-bem unsteady wake, advected as vortex particles; the metric measures
//! circulation-weighted streamwise rearrangement, not thrust). Generations
//! are E-RACED through fs-race: dominated candidates die at a fraction
//! of their solve budget with anytime validity, and the payoff is
//! MEASURED (evaluations used vs the fixed-N equivalent), ledgered in
//! the report.

use crate::param::OrnithCandidate;
use fs_bem::panel2d::solve;
use fs_bem::wake2d::{WakeSim, WakeVortex};
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

/// A screening measure of wake rearrangement under the prescribed gait gain.
/// The BEM wake is handed to VPM with its physical circulation sign retained.
/// The result is `sum |Gamma_i| |x_i(T)-x_i(0)| / (T sum |Gamma_i|)` in
/// freestream-speed units, or zero for a wake with no circulation.
///
/// This measures the amount of streamwise motion of the represented vorticity,
/// independently of the coordinate origin or zero-strength tracer particles.
/// It is neither signed centroid drift nor a thrust/force prediction. The gain
/// modulation remains a screening ansatz, not resolved flapping kinematics.
#[must_use]
pub fn flap_metric(c: &OrnithCandidate) -> f64 {
    let foil = c.section(PANELS);
    let mut sim = WakeSim::new(&foil, c.alpha, 0.08, 0.05)
        .expect("ornith candidate must satisfy the bounded wake contract");
    for _ in 0..24 {
        sim.step()
            .expect("short bounded wake fixture must remain finite");
    }
    let particles = gait_particles(sim.wake(), c.flap_amp, c.flap_freq);
    let moved = advect_steps(&particles, 0.1, 12, 0.05);
    rearrangement_rate(&particles, &moved, 1.2)
}

fn gait_particles(wake: &[WakeVortex], amplitude: f64, frequency: f64) -> Vec<VortexParticle> {
    wake.iter()
        .enumerate()
        .map(|(k, w)| {
            let phase = frequency * k as f64;
            let gain = amplitude.mul_add(fs_math::det::sin(phase), 1.0);
            // BEM is clockwise-positive; VPM is counterclockwise-positive.
            // Convert once, before the original VPM transport kernels.
            VortexParticle::new(w.pos, -w.gamma * gain)
        })
        .collect()
}

fn rearrangement_rate(initial: &[VortexParticle], moved: &[VortexParticle], elapsed: f64) -> f64 {
    debug_assert_eq!(initial.len(), moved.len());
    let strength: f64 = initial.iter().map(|p| p.circulation.abs()).sum();
    if strength == 0.0 {
        return 0.0;
    }
    initial
        .iter()
        .zip(moved)
        .map(|(p, q)| p.circulation.abs() * (q.pos[0] - p.pos[0]).abs())
        .sum::<f64>()
        / (elapsed * strength)
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
    fn clockwise_bem_dipole_keeps_its_physical_transport_direction() {
        let wake = [
            WakeVortex {
                pos: [0.0, 0.5],
                gamma: 1.0,
            },
            WakeVortex {
                pos: [0.0, -0.5],
                gamma: -1.0,
            },
        ];
        let particles = gait_particles(&wake, 0.0, 1.0);
        let moved = advect_steps(&particles, 0.05, 4, 0.05);
        // Clockwise positive above counterclockwise negative translates left.
        // Its separation is constant, so this is the exact finite-core speed.
        let speed = -1.0 / (std::f64::consts::TAU * (1.0 + 0.05 * 0.05));
        for (p, q) in particles.iter().zip(&moved) {
            assert!((q.pos[0] - p.pos[0] - 0.2 * speed).abs() < 2e-15);
            assert!((q.pos[1] - p.pos[1]).abs() < 2e-15);
        }
        assert!((rearrangement_rate(&particles, &moved, 0.2) + speed).abs() < 2e-15);
    }

    #[test]
    fn actual_bem_screen_handoff_preserves_the_gained_induced_field() {
        let candidate = OrnithCandidate::from_genes(&[0.5, 0.5, 0.5, 0.5, 0.5]);
        let mut sim =
            WakeSim::new(&candidate.section(PANELS), candidate.alpha, 0.08, 0.05).unwrap();
        for _ in 0..24 {
            sim.step().unwrap();
        }
        let particles = gait_particles(sim.wake(), candidate.flap_amp, candidate.flap_freq);
        assert!(!particles.is_empty());
        for point in [[0.75, 0.1], [1.5, -0.2], [3.0, 0.4]] {
            let mut expected = [0.0; 2];
            for (k, w) in sim.wake().iter().enumerate() {
                let gain = candidate
                    .flap_amp
                    .mul_add(fs_math::det::sin(candidate.flap_freq * k as f64), 1.0);
                let dx = point[0] - w.pos[0];
                let dy = point[1] - w.pos[1];
                let denominator = std::f64::consts::TAU * (dx * dx + dy * dy + 0.05 * 0.05);
                // Independently evaluate BEM's clockwise-positive velocity.
                expected[0] += w.gamma * gain * dy / denominator;
                expected[1] -= w.gamma * gain * dx / denominator;
            }
            let actual = fs_vpm::induced_velocity(&particles, point, 0.05);
            assert!(expected[0].abs() + expected[1].abs() > 1e-5);
            for axis in 0..2 {
                assert!((actual[axis] - expected[axis]).abs() < 2e-14);
            }
        }
        let moved = advect_steps(&particles, 0.1, 12, 0.05);
        assert_eq!(
            flap_metric(&candidate),
            rearrangement_rate(&particles, &moved, 1.2)
        );
    }

    #[test]
    fn rearrangement_uses_circulation_and_displacements_not_world_positions() {
        let initial = [
            VortexParticle::new([-2.0, 0.0], 2.0),
            VortexParticle::new([1.0, 0.0], -1.0),
            VortexParticle::new([10.0, 0.0], 0.0),
        ];
        let moved = [
            VortexParticle::new([-1.0, 0.0], 2.0),
            VortexParticle::new([-1.0, 0.0], -1.0),
            VortexParticle::new([100.0, 0.0], 0.0),
        ];
        let expected = 2.0 / 3.0; // (2*1 + 1*2) / (2 seconds * 3 strength).
        assert_eq!(rearrangement_rate(&initial, &moved, 2.0), expected);
        let translate = |values: [VortexParticle; 3]| {
            values.map(|p| VortexParticle::new([p.pos[0] + 32.0, p.pos[1] - 5.0], p.circulation))
        };
        assert_eq!(
            rearrangement_rate(&translate(initial), &translate(moved), 2.0),
            expected
        );
        assert_eq!(
            rearrangement_rate(&initial[..2], &moved[..2], 2.0),
            expected
        );
        assert_eq!(rearrangement_rate(&initial[2..], &moved[2..], 2.0), 0.0);
        assert_eq!(rearrangement_rate(&[], &[], 2.0), 0.0);
    }
}
