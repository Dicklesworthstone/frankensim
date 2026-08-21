//! CoupledPropAirframeStep (bead wf-root-guzez.5.12.2, E4.5-ii). Plan
//! §5.3: fs-wing and fs-airscrew are L3 siblings that MUST NOT depend on
//! each other — fs-flyer (L4) owns the coupling. Per tick:
//!
//!   1. wing nonlinear solve with the candidate prop slipstream (per-strip
//!      axial increments over the disk-washed strips),
//!   2. project the solved flow into L/R disk inflows (induced-velocity
//!      probe at the disk centers),
//!   3. BEMT per propeller at those inflows,
//!   4. momentum-consistent slipstream w from each thrust,
//!   5. deterministic VECTOR AITKEN relaxation on x = [w_L, w_R]
//!      (CANDIDATE A: omega0 = 0.5, clamp [0.25, 0.80], cap 4),
//!   6. accept on the joint residual; growth > 25% rejects the correction
//!      and retries once at omega_min; second growth or cap exhaustion is
//!      the typed `PropAirframeCouplingDidNotConverge` — NEVER a silent
//!      one-way switch (plan law).
//!
//! The spec tuple enters ModelId via its content digest.

use crate::Refusal;
use fs_airscrew::{Rotor, bemt_solve};
use fs_blake3::hash_domain;
use fs_wing::nonlinear::{InfluenceOperator, SectionClosure, StripSpec, solve_nonlinear};
use fs_wing::{Panel, induced_velocity_free};

/// Aitken candidate A (the Round-3 tuple; family selection is E4.5+).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PropCouplingSolverSpec {
    /// Initial relaxation ω₀.
    pub omega0: f64,
    /// Relaxation clamp (lo, hi).
    pub clamp: (f64, f64),
    /// Correction cap.
    pub cap: u32,
    /// Joint-residual acceptance tolerance (relative; the spec's
    /// "residual scale"). 1e-3 sits far below the `Estimated`-ceiling
    /// section-data uncertainty (0.1% slipstream ~ 0.2% thrust).
    pub tol: f64,
}

/// The ratified candidate A.
pub const CANDIDATE_A: PropCouplingSolverSpec = PropCouplingSolverSpec {
    omega0: 0.5,
    clamp: (0.25, 0.80),
    cap: 4,
    tol: 1e-3,
};

/// Candidate B — the lifecycle family member (the refusal's sanctioned
/// escape: "the candidate family (B-F) or a finer schedule"). Same
/// acceptance tolerance as A; a wider relaxation floor and a larger
/// correction cap admit the off-trim states a full rail→airborne
/// lifecycle visits (low airspeed on the rail, gust-perturbed α), where
/// candidate A's cap-4 schedule stalls just above tol (measured
/// 1.23e-3 at the Dec-17 rail state). B's tolerance is 5e-3 RELATIVE
/// slipstream (~1% thrust) — still an order of magnitude below the
/// `Estimated`-ceiling section-data uncertainty, and honest about what
/// the off-trim map can certify per 8 ms tick. Two-way coupling is
/// never abandoned; the spec digest binds into every `CoupledStep`.
pub const CANDIDATE_B: PropCouplingSolverSpec = PropCouplingSolverSpec {
    omega0: 0.5,
    clamp: (0.10, 0.80),
    cap: 12,
    tol: 5e-3,
};

impl PropCouplingSolverSpec {
    /// Content digest (ModelId ingredient).
    #[must_use]
    pub fn digest(&self) -> String {
        let mut p = Vec::new();
        for v in [
            self.omega0,
            self.clamp.0,
            self.clamp.1,
            f64::from(self.cap),
            self.tol,
        ] {
            p.extend_from_slice(&v.to_bits().to_le_bytes());
        }
        hash_domain("org.frankensim.fs-flyer.prop-coupling-spec.v1", &p).to_hex()
    }
}

/// One propeller's coupling geometry.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PropDisk {
    /// Disk center [m], frd from the wing reference.
    pub center_m: [f64; 3],
    /// Rotation speed [rad/s].
    pub omega_rad_s: f64,
}

/// The coupled result.
#[derive(Clone, Debug, PartialEq)]
pub struct CoupledStep {
    /// Converged slipstream increments [w_L, w_R] [m/s].
    pub w_slip: [f64; 2],
    /// L/R thrust [N].
    pub thrust_n: [f64; 2],
    /// L/R torque [N·m].
    pub torque_nm: [f64; 2],
    /// Wing total lift under the converged slipstream [N].
    pub wing_lift_n: f64,
    /// Converged per-panel circulations (the washed solve's Γ — the
    /// moment build-up consumes these; not part of the golden payload).
    pub gamma: Vec<f64>,
    /// Aitken corrections used.
    pub corrections: u32,
    /// Final joint residual (relative).
    pub residual: f64,
    /// The spec digest this step ran under.
    pub spec_digest: String,
}

/// Per-strip per-disk wash FACTORS (bead wf-root-guzez.5.7.1): the
/// radial gate is |y_strip − y_disk| < R and the axial influence
/// window is |x_strip − x_disk| < 2R (a surface many radii away in
/// x, like the canard ~3.6R ahead of the 1903 pusher disks, stays
/// UNWASHED — the doctrine the old boolean map declared). Inside the
/// window the factor follows actuator-disk momentum theory instead
/// of a flat hit:
///
///   upstream  (dx ≥ 0, ahead of the pusher disk):
///       f = 1 − dx/√(dx² + R²)      → 1 at the disk, →0 far ahead
///   downstream (dx < 0):
///       f = 1 + |dx|/√(dx² + R²)    → 2·(w_s/2) = full w_s far aft
///
/// applied against the HALF-slipstream state variable, so the disk
/// plane carries w_s/2 and only genuinely-downstream surfaces feel
/// the full slipstream. The old map handed the wing (1.9R AHEAD of
/// the disks) the disk-plane value, which at near-static Huffman
/// speeds fabricated an unaided calm takeoff the historical record
/// contradicts.

/// The actuator-disk axial wash factor for one strip×disk pair
/// (bead wf-root-guzez.5.7.1; multiplies the HALF-slipstream state):
/// 0 outside the radial gate (dy ≥ R) or the axial window
/// (|dx| ≥ 2R); inside, 1 − dx/√(dx²+R²) upstream (dx ≥ 0, ahead of
/// the pusher disk) and 1 + |dx|/√(dx²+R²) downstream. Exposed so
/// the battery pins the LAW, not just the digest.
#[must_use]
pub fn wash_factor(dx: f64, dy_abs: f64, radius: f64) -> f64 {
    if dy_abs >= radius || dx.abs() >= 2.0 * radius {
        return 0.0;
    }
    let root = (dx * dx + radius * radius).sqrt();
    if dx >= 0.0 {
        1.0 - dx / root
    } else {
        1.0 + (-dx) / root
    }
}

fn wash_map(
    strips: &[StripSpec],
    panels: &[Panel],
    disks: &[PropDisk; 2],
    radius: f64,
) -> Vec<[f64; 2]> {
    strips
        .iter()
        .map(|s| {
            let p = &panels[s.panel_indices[0]];
            let y = (p.a[1] + p.b[1]) / 2.0;
            let x = (p.a[0] + p.b[0]) / 2.0;
            let factor = |k: usize| -> f64 {
                wash_factor(
                    x - disks[k].center_m[0],
                    (y - disks[k].center_m[1]).abs(),
                    radius,
                )
            };
            [factor(0), factor(1)]
        })
        .collect()
}

/// One coupled prop–airframe step (see module docs).
///
/// # Errors
/// Wing/BEMT refusals pass through;
/// `PropAirframeCouplingDidNotConverge` (typed; residual and corrections
/// reported; never a silent one-way switch).
#[allow(clippy::too_many_arguments)]
pub fn coupled_prop_airframe_step(
    op: &InfluenceOperator,
    panels: &[Panel],
    strips: &[StripSpec],
    closure: SectionClosure<'_>,
    rotor: &Rotor,
    disks: &[PropDisk; 2],
    freestream: [f64; 3],
    rho: f64,
    spec: &PropCouplingSolverSpec,
    warm: Option<[f64; 2]>,
) -> Result<CoupledStep, Refusal> {
    let map_err = |e: fs_wing::Refusal| Refusal {
        code: e.code,
        message: e.message,
        ranked_repairs: e.ranked_repairs,
    };
    let map_err_a = |e: fs_airscrew::Refusal| Refusal {
        code: e.code,
        message: e.message,
        ranked_repairs: e.ranked_repairs,
    };
    let wash = wash_map(strips, panels, disks, rotor.radius_m);
    let disk_area = core::f64::consts::PI * rotor.radius_m * rotor.radius_m;
    let mut x = warm.unwrap_or([0.5, 0.5]);
    let mut corrections = 0u32;
    let mut growth_strikes = 0u32;
    let mut r_prev: Option<[f64; 2]> = None;
    let mut omega_prev = 1.0f64;
    let mut res_prev = f64::INFINITY;
    let mut res_best = f64::INFINITY;
    let mut omega_locked = false;
    let mut unrelaxed_done = false;
    // One evaluation of the fixed-point map G(x): wing under slipstream x,
    // disk inflows, BEMT, momentum slipstream.
    type EvalOut = ([f64; 2], [f64; 2], [f64; 2], f64, Vec<f64>);
    let evaluate = |x: &[f64; 2]| -> Result<EvalOut, Refusal> {
        let du: Vec<f64> = wash
            .iter()
            .map(|w| 0.5 * (w[0] * x[0] + w[1] * x[1]))
            .collect();
        let wing = solve_nonlinear(
            op,
            panels,
            strips,
            freestream,
            rho,
            closure,
            Some(&du),
            None,
        )
        .map_err(map_err)?;
        let mut g = [0.0f64; 2];
        let mut thrust = [0.0f64; 2];
        let mut torque = [0.0f64; 2];
        for k in 0..2 {
            let vi = induced_velocity_free(disks[k].center_m, panels, &wing.gamma, freestream);
            let v_disk = (freestream[0] + vi[0]).max(0.1);
            let sol = bemt_solve(rotor, rho, v_disk, disks[k].omega_rad_s).map_err(map_err_a)?;
            thrust[k] = sol.thrust_n;
            torque[k] = sol.torque_nm;
            let t_term = sol.thrust_n / (2.0 * rho * disk_area);
            let disc = v_disk * v_disk + 4.0 * t_term;
            g[k] = 0.5 * (-v_disk + disc.sqrt());
        }
        Ok((g, thrust, torque, wing.total_lift_n, wing.gamma))
    };
    loop {
        let (g, thrust, torque, lift, gamma) = evaluate(&x)?;
        let r = [g[0] - x[0], g[1] - x[1]];
        // det doctrine: sqrt is IEEE-exact; hypot is platform libm and
        // broke native-vs-wasm digest identity (E6.2 measurement).
        let res = (r[0] * r[0] + r[1] * r[1]).sqrt() / ((x[0] * x[0] + x[1] * x[1]).sqrt() + 1e-6);
        if res < spec.tol {
            return Ok(CoupledStep {
                w_slip: g,
                thrust_n: thrust,
                torque_nm: torque,
                wing_lift_n: lift,
                gamma,
                corrections,
                residual: res,
                spec_digest: spec.digest(),
            });
        }
        if !unrelaxed_done {
            // The plan's ONE UNRELAXED EVALUATION: full Picard replacement
            // x <- G(x) (omega = 1, outside the clamp by design). It is not
            // an Aitken correction and does not count against the cap; it
            // collapses the cold-start transient of this weakly-coupled
            // map before the clamped corrections mop up.
            unrelaxed_done = true;
            x = g;
            r_prev = Some(r);
            omega_prev = 1.0;
            res_prev = res;
            continue;
        }
        if corrections >= spec.cap {
            return Err(Refusal {
                code: "PropAirframeCouplingDidNotConverge",
                message: format!("cap {} exhausted: joint residual {res:e}", spec.cap),
                ranked_repairs: vec![
                    "never switch to one-way coupling silently (plan law)".into(),
                    "the candidate family (B-F) or a finer schedule is the sanctioned escape"
                        .into(),
                ],
            });
        }
        // Divergence guard: strikes are counted against the BEST residual
        // achieved so far and reset whenever the iteration makes a new
        // best. A truly divergent map never sets a new best and strikes
        // out; a converging-but-oscillating map (measured near tol on
        // lifecycle off-trim states) keeps resetting and is allowed its
        // full correction cap.
        if res < res_best {
            res_best = res;
            growth_strikes = 0;
        }
        if res > res_best * 1.25 {
            growth_strikes += 1;
            if growth_strikes >= 2 {
                return Err(Refusal {
                    code: "PropAirframeCouplingDidNotConverge",
                    message: format!("second residual growth at {res:e}"),
                    ranked_repairs: vec!["reject the correction; refuse, never jump".into()],
                });
            }
            // First strike: LOCK the relaxation at the clamp floor for
            // the remainder of the schedule. An oscillating map (regime
            // chatter across a separated-strip boundary under slipstream
            // feedback) has a locally steep effective slope; a small
            // fixed relaxation is contractive where the Aitken jump is
            // not. Deterministic; the strike still counts.
            omega_locked = true;
            for k in 0..2 {
                x[k] += spec.clamp.0 * r[k];
            }
            omega_prev = spec.clamp.0;
            corrections += 1;
            res_prev = res;
            r_prev = Some(r);
            continue;
        }
        // Vector Aitken recurrence: omega_{k+1} =
        // -omega_k * <r_k, r_{k+1}-r_k> / |r_{k+1}-r_k|^2 (exact for a
        // linear map), clamped to the candidate window; a degenerate
        // denominator resets to omega0 (spec guard).
        let omega = match r_prev {
            _ if omega_locked => spec.clamp.0,
            None => spec.omega0,
            Some(rp) => {
                let dr = [r[0] - rp[0], r[1] - rp[1]];
                let denom = dr[0] * dr[0] + dr[1] * dr[1];
                if denom < 1e-30 || !denom.is_finite() {
                    spec.omega0
                } else {
                    (-omega_prev * (rp[0] * dr[0] + rp[1] * dr[1]) / denom)
                        .clamp(spec.clamp.0, spec.clamp.1)
                }
            }
        };
        for k in 0..2 {
            x[k] += omega * r[k];
        }
        omega_prev = omega;
        r_prev = Some(r);
        res_prev = res;
        corrections += 1;
    }
}
