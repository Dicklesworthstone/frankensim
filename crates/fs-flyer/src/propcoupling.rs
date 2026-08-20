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

/// Strips washed by a disk: |y_strip − y_disk| < R AND within the axial
/// influence window |x_strip − x_disk| < 2R (an actuator disk's velocity
/// perturbation decays over ~R upstream and the full slipstream lives
/// downstream — a surface many radii away in x, like the canard 5 m
/// ahead of the 1903 pusher disks, is NOT washed).
fn wash_map(
    strips: &[StripSpec],
    panels: &[Panel],
    disks: &[PropDisk; 2],
    radius: f64,
) -> Vec<[bool; 2]> {
    strips
        .iter()
        .map(|s| {
            let p = &panels[s.panel_indices[0]];
            let y = (p.a[1] + p.b[1]) / 2.0;
            let x = (p.a[0] + p.b[0]) / 2.0;
            let hit = |k: usize| -> bool {
                (y - disks[k].center_m[1]).abs() < radius
                    && (x - disks[k].center_m[0]).abs() < 2.0 * radius
            };
            [hit(0), hit(1)]
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
    let mut unrelaxed_done = false;
    // One evaluation of the fixed-point map G(x): wing under slipstream x,
    // disk inflows, BEMT, momentum slipstream.
    type EvalOut = ([f64; 2], [f64; 2], [f64; 2], f64, Vec<f64>);
    let evaluate = |x: &[f64; 2]| -> Result<EvalOut, Refusal> {
        let du: Vec<f64> = wash
            .iter()
            .map(|w| 0.5 * (if w[0] { x[0] } else { 0.0 } + if w[1] { x[1] } else { 0.0 }))
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
        let res = (r[0] * r[0] + r[1] * r[1]).sqrt() / (x[0].hypot(x[1]) + 1e-6);
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
        if res > res_prev * 1.25 {
            growth_strikes += 1;
            if growth_strikes >= 2 {
                return Err(Refusal {
                    code: "PropAirframeCouplingDidNotConverge",
                    message: format!("second residual growth at {res:e}"),
                    ranked_repairs: vec!["reject the correction; refuse, never jump".into()],
                });
            }
            // Retry once from the same state at omega_min.
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
