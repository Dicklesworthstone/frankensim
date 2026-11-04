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
//!   6. accept on the joint residual; growth > 25% against the correction's
//!      accepted predecessor rejects the entire trial and retries once from
//!      that predecessor at omega_min; second growth or cap exhaustion is
//!      the typed `PropAirframeCouplingDidNotConverge` — NEVER a silent
//!      one-way switch (plan law).
//! Each wing evaluation continues from the accepted predecessor's circulation.
//! A rejected outer trial does not replace that anchor. The first evaluation
//! still selects its cold branch; this is within-solve, not across-tick memory.
//!
//! The spec tuple enters ModelId via its content digest.

use crate::Refusal;
use fs_airscrew::{BemtPower, Rotor, bemt_solve};
use fs_blake3::hash_domain;
use fs_wing::nonlinear::continuation::{self, ContinuationControls, WingContinuation};
use fs_wing::nonlinear::{
    InfluenceOperator, SectionClosure, SectionFrame, StripSpec, solve_nonlinear,
    solve_nonlinear_with_incident, solve_nonlinear_with_section_frames,
};
use fs_wing::{Panel, induced_velocity_free};

#[path = "propcoupling/spatial.rs"]
pub mod spatial;

#[path = "propcoupling/swirl.rs"]
pub mod swirl;

#[path = "propcoupling/disk_wash.rs"]
pub(crate) mod disk_wash;

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
        // v3 binds accepted-circulation continuation within the outer solve.
        // v2 established the common disk-induced wash convention; v1 had an
        // additional, erroneous factor 1/2 in the legacy loading path.
        hash_domain("org.frankensim.fs-flyer.prop-coupling-spec.v3", &p).to_hex()
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

/// Incident relative-air velocities before prop-induced wash, all sampled at
/// one physical tick in body FRD coordinates. The rotor model consumes the
/// separate disk-center axial components, not disk harmonics or crossflow.
#[derive(Clone, Copy, Debug)]
pub struct IncidentInflow<'a> {
    /// One local velocity per wing/canard strip.
    pub strip_velocity_mps: &'a [[f64; 3]],
    /// Separate left/right disk-center velocities.
    pub disk_velocity_mps: [[f64; 3]; 2],
}

/// The coupled result.
#[derive(Clone, Debug, PartialEq)]
pub struct CoupledStep<P = BemtPower> {
    /// Converged disk-induced axial velocities [vi_L, vi_R] [m/s].
    pub w_slip: [f64; 2],
    /// L/R thrust [N].
    pub thrust_n: [f64; 2],
    /// L/R torque [N·m].
    pub torque_nm: [f64; 2],
    /// L/R signed BEMT power at the accepted disk inflows and shaft speeds.
    /// Thrust power uses local disk incident air, including the wing field;
    /// it is not the complete aircraft/wake energy balance. Induced/profile
    /// terms retain their original annular meanings and are not added as heat.
    /// The default is the original axial `BemtPower`; the explicitly selected
    /// swirl coupling retains its separate angular wake power in `P`.
    pub power: [P; 2],
    /// Wing total lift under the converged slipstream [N].
    pub wing_lift_n: f64,
    /// Converged per-panel circulations (the washed solve's Γ — the
    /// moment build-up consumes these; not part of the golden payload).
    pub gamma: Vec<f64>,
    /// Accepted effective section incidences [rad], including twist and the
    /// wing solver's trailing-sheet induction. Paired with `gamma`.
    pub strip_alpha_eff_rad: Vec<f64>,
    /// Accepted local strip velocities including prop wash. Both compact
    /// uniform and explicit incident inputs return `Some`; force assembly
    /// must use these with `gamma`, not the next `w_slip` predictor.
    pub strip_velocity_mps: Option<Vec<[f64; 3]>>,
    /// Correction attempts used, including a rejected trial and its retry.
    pub corrections: u32,
    /// Final relative residual: joint L2 in legacy mode; maximum per-disk
    /// relative residual in distributed mode (1e-6 m/s denominator safeguard).
    pub residual: f64,
    /// The spec digest this step ran under.
    pub spec_digest: String,
}

/// Per-strip per-disk wash FACTORS (bead wf-root-guzez.5.7.1): the
/// radial gate is |y_strip − y_disk| < R and the axial influence
/// window is |x_strip − x_disk| < 2R (a surface many radii away in
/// x, like the canard ~3.6R ahead of the 1903 pusher disks, stays
/// UNWASHED — the doctrine the old boolean map declared). Inside the
/// window the factor multiplies the DISK-INDUCED velocity `vi`:
///
///   upstream  (dx ≥ 0, ahead of the pusher disk):
///       f = 1 − dx/√(dx² + R²)      → 1 at the disk, →0 far ahead
///   downstream (dx < 0):
///       f = 1 + |dx|/√(dx² + R²)
///
/// At the disk the increment is `vi`; the untruncated axial law tends to
/// `2 vi` far downstream. That factor of two is already in this spatial
/// law: multiplying `vi` by another 1/2 would halve the induction again.
/// The retained finite 2R window is a model truncation, not a far-wake
/// calculation. This spanwise-overlap model also omits wake contraction,
/// swirl, radial loading and vertical disk overlap.

/// The actuator-disk axial wash factor for one strip×disk pair
/// (bead wf-root-guzez.5.7.1; multiplies the DISK-INDUCED velocity `vi`):
/// 0 outside the radial gate (dy ≥ R) or the axial window
/// (|dx| ≥ 2R); inside, 1 − dx/√(dx²+R²) upstream (dx ≥ 0, ahead of
/// the pusher disk) and 1 + |dx|/√(dx²+R²) downstream. Exposed so
/// the battery pins the LAW, not just the digest.
#[must_use]
pub fn wash_factor(dx: f64, dy_abs: f64, radius: f64) -> f64 {
    if dy_abs >= radius {
        return 0.0;
    }
    disk_wash::axial_factor(dx, radius)
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

/// `T = rho A (V + vi) (2 vi)`: the disk increment is `vi`, while the
/// ideal far-wake increment is `2 vi`. `wash_factor` owns the spatial
/// factor, so both loading modes use this same sum without another 1/2.
/// See NASA Glenn, "Propeller Thrust", and Deters/Ananda/Selig (2015),
/// AIAA 2015-2265, equations (3)--(7).
#[inline]
fn strip_wash(factors: [f64; 2], induction: [f64; 2]) -> f64 {
    // Statement-split to preserve the deterministic no-FMA convention.
    let left = factors[0] * induction[0];
    let right = factors[1] * induction[1];
    left + right
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
    coupled_step_impl(
        op, panels, strips, closure, rotor, disks, freestream, rho, spec, warm, None, None, None,
    )
}

/// Two-way propeller/wing solve with nonuniform incident velocities. The
/// converged momentum state is disk-induced velocity `vi`, satisfying
/// `T = 2 rho A vi (V_disk + vi)`; strip wash is `wash_factor * vi`.
/// Local section incidence, circulation, and returned force velocities all
/// receive this same wash. Both loading modes use the disk-induced velocity
/// convention and return the accepted washed force velocities.
/// Each disk must meet its own relative convergence tolerance, so a strongly
/// loaded disk cannot hide an unconverged weak disk. The momentum root uses
/// its rationalized form to retain small nonzero induction at weak loading.
///
/// # Errors
/// As [`coupled_prop_airframe_step`]; `disk-inflow-invalid` for non-finite
/// disk samples or reversed net axial disk inflow; `disk-momentum-invalid`
/// for a non-real/non-finite momentum root; local strip refusals.
#[allow(clippy::too_many_arguments)]
pub fn coupled_prop_airframe_step_with_incident(
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
    incident: IncidentInflow<'_>,
) -> Result<CoupledStep, Refusal> {
    coupled_step_impl(
        op,
        panels,
        strips,
        closure,
        rotor,
        disks,
        freestream,
        rho,
        spec,
        warm,
        Some(incident),
        None,
        None,
    )
}

#[allow(clippy::too_many_arguments)]
fn coupled_step_impl(
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
    incident: Option<IncidentInflow<'_>>,
    spatial: Option<spatial::DiskAxes>,
    section_frames: Option<&[SectionFrame]>,
) -> Result<CoupledStep, Refusal> {
    let map_err_a = |e: fs_airscrew::Refusal| Refusal {
        code: e.code,
        message: e.message,
        ranked_repairs: e.ranked_repairs,
    };
    coupled_step_with_rotors(
        op,
        panels,
        strips,
        closure,
        rotor,
        disks,
        freestream,
        rho,
        spec,
        warm,
        incident,
        spatial,
        section_frames,
        None,
        |_, axial, omega, _applied_wash| {
            let sol = bemt_solve(rotor, rho, axial, omega).map_err(map_err_a)?;
            Ok(RotorLoads {
                thrust_n: sol.thrust_n,
                torque_nm: sol.torque_nm,
                power: sol.power,
            })
        },
    )
}

struct RotorLoads<P> {
    thrust_n: f64,
    torque_nm: f64,
    power: P,
}

/// The same two-way iteration and acceptance for explicitly selected rotor
/// laws. Only the actual rotor evaluation and its power payload vary.
// One caller-owned branch remains fixed through every nonlinear trial. Only
// the wing whose residual/loads are accepted may replace this local output.
struct WingBranch<'a> {
    domain: ContinuationControls,
    previous: Option<&'a WingContinuation>,
    accepted: Option<WingContinuation>,
}

#[allow(clippy::too_many_arguments)]
fn coupled_step_with_rotors<P: Copy>(
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
    incident: Option<IncidentInflow<'_>>,
    spatial: Option<spatial::DiskAxes>,
    section_frames: Option<&[SectionFrame]>,
    mut branch: Option<&mut WingBranch<'_>>,
    mut rotor_solve: impl FnMut(usize, f64, f64, [f64; 2]) -> Result<RotorLoads<P>, Refusal>,
) -> Result<CoupledStep<P>, Refusal> {
    let map_err = |e: fs_wing::Refusal| Refusal {
        code: e.code,
        message: e.message,
        ranked_repairs: e.ranked_repairs,
    };
    // wash_map reads the front panel before the nonlinear wing admission.
    // Admit those indices here so malformed public input cannot panic.
    if strips.is_empty()
        || strips.iter().any(|strip| {
            strip.panel_indices.is_empty() || strip.panel_indices.iter().any(|&i| i >= panels.len())
        })
    {
        return Err(Refusal {
            code: "strips-invalid",
            message: "empty strip or panel index outside the supplied layout".into(),
            ranked_repairs: vec!["supply nonempty strips referencing existing panels".into()],
        });
    }
    if incident.is_some_and(|flow| {
        !flow
            .disk_velocity_mps
            .iter()
            .flatten()
            .all(|v| v.is_finite())
    }) {
        return Err(Refusal {
            code: "disk-inflow-invalid",
            message: "non-finite incident disk velocity".into(),
            ranked_repairs: vec!["sample both disks from the admitted physical atmosphere".into()],
        });
    }
    if section_frames.is_some_and(|frames| {
        frames.len() != strips.len()
            || incident.is_none_or(|flow| flow.strip_velocity_mps.len() != strips.len())
    }) {
        return Err(Refusal {
            code: "section-frame-invalid",
            message: "provide one section frame and incident velocity per strip".into(),
            ranked_repairs: vec![
                "retain the admitted frames with the supplied strip layout".into(),
            ],
        });
    }
    let wash = match spatial {
        Some(axes) => axes.wash_map(strips, panels, disks, rotor.radius_m)?,
        None => wash_map(strips, panels, disks, rotor.radius_m),
    };
    // Exchange the same candidate wake between disks as between disks and
    // strips. Self-induction stays exclusively in the native BEMT/momentum law.
    let disk_axes = match spatial {
        Some(axes) => axes,
        None => spatial::DiskAxes::new([[1.0, 0.0, 0.0]; 2])?,
    };
    let disk_wash = disk_axes.disk_wash_map(disks, rotor.radius_m)?;
    let disk_area = core::f64::consts::PI * rotor.radius_m * rotor.radius_m;
    let mut x = warm.unwrap_or([0.5, 0.5]);
    let mut x_prev = x;
    let mut corrections = 0u32;
    let mut retry_used = false;
    let mut r_prev: Option<[f64; 2]> = None;
    let mut omega_prev = 1.0f64;
    let mut res_prev = f64::INFINITY;
    let mut unrelaxed_done = false;
    // This anchor belongs to the evaluated predecessor of the next proposal.
    // Keep it through rejected trials, just like x_prev and r_prev. Moving the
    // already-owned gamma vector adds no copy of the accepted wing state.
    let mut accepted_gamma: Option<Vec<f64>> = None;
    // One evaluation of the fixed-point map G(x): wing under slipstream x,
    // disk inflows, BEMT, momentum slipstream.
    type EvalOut<P> = (
        [f64; 2],
        [f64; 2],
        [f64; 2],
        [P; 2],
        fs_wing::nonlinear::NonlinearReport,
        Option<WingContinuation>,
    );
    let mut evaluate = |x: &[f64; 2], warm_gamma: Option<&[f64]>| -> Result<EvalOut<P>, Refusal> {
        let du: Vec<f64> = if spatial.is_none() {
            wash.iter()
                .map(|&factors| strip_wash(factors, *x))
                .collect()
        } else {
            Vec::new()
        };
        let (wing, wing_anchor) = if let Some(branch) = branch.as_ref() {
            let velocity: Vec<_> = wash
                .iter()
                .enumerate()
                .map(|(s, &factors)| {
                    let velocity = incident.map_or(freestream, |flow| flow.strip_velocity_mps[s]);
                    match spatial {
                        Some(axes) => axes.velocity(velocity, factors, *x),
                        None => [
                            velocity[0] + strip_wash(factors, *x),
                            velocity[1],
                            velocity[2],
                        ],
                    }
                })
                .collect();
            let solved = if let Some(frames) = section_frames {
                continuation::solve_with_section_frames(
                    op,
                    panels,
                    strips,
                    freestream,
                    rho,
                    closure,
                    frames,
                    &velocity,
                    branch.previous,
                    branch.domain,
                )
            } else {
                continuation::solve_with_incident(
                    op,
                    panels,
                    strips,
                    freestream,
                    rho,
                    closure,
                    &velocity,
                    None,
                    branch.previous,
                    branch.domain,
                )
            }
            .map_err(map_err)?;
            (solved.solution, Some(solved.continuation))
        } else {
            let wing = if let Some(frames) = section_frames {
                let flow = incident.expect("framed coupling requires incident samples");
                let velocity: Vec<_> = flow
                    .strip_velocity_mps
                    .iter()
                    .zip(&wash)
                    .map(|(&v, &factors)| {
                        if let Some(axes) = spatial {
                            axes.velocity(v, factors, *x)
                        } else {
                            [v[0] + strip_wash(factors, *x), v[1], v[2]]
                        }
                    })
                    .collect();
                // Compose actual spatial prop wash before the section projection.
                // The same accepted vectors are returned with gamma below.
                solve_nonlinear_with_section_frames(
                    op, panels, strips, freestream, rho, closure, frames, &velocity, warm_gamma,
                )
            } else {
                match (incident, spatial) {
                    (Some(flow), Some(axes)) => {
                        let velocity: Vec<_> = flow
                            .strip_velocity_mps
                            .iter()
                            .zip(&wash)
                            .map(|(&v, &factors)| axes.velocity(v, factors, *x))
                            .collect();
                        solve_nonlinear_with_incident(
                            op, panels, strips, freestream, rho, closure, &velocity, None,
                            warm_gamma,
                        )
                    }
                    (Some(flow), None) => solve_nonlinear_with_incident(
                        op,
                        panels,
                        strips,
                        freestream,
                        rho,
                        closure,
                        flow.strip_velocity_mps,
                        Some(&du),
                        warm_gamma,
                    ),
                    (None, _) => solve_nonlinear(
                        op,
                        panels,
                        strips,
                        freestream,
                        rho,
                        closure,
                        Some(&du),
                        warm_gamma,
                    ),
                }
            }
            .map_err(map_err)?;
            (wing, None)
        };
        let mut g = [0.0f64; 2];
        let mut thrust = [0.0f64; 2];
        let mut torque = [0.0f64; 2];
        let mut power = [None; 2];
        for k in 0..2 {
            let vi = induced_velocity_free(disks[k].center_m, panels, &wing.gamma, freestream);
            let mut disk_incident = incident.map_or(freestream, |flow| flow.disk_velocity_mps[k]);
            let disk_factors = disk_wash[k].map(|wash| wash.mean);
            if disk_factors != [0.0; 2] {
                disk_incident = match spatial {
                    Some(axes) => axes.velocity(disk_incident, disk_factors, *x),
                    None => {
                        disk_incident[0] += strip_wash(disk_factors, *x);
                        disk_incident
                    }
                };
            }
            let axial = disk_incident[0];
            let v_disk = spatial.map_or(axial + vi[0], |axes| axes.axial(k, disk_incident, vi));
            if !v_disk.is_finite() || v_disk < 0.0 {
                return Err(Refusal {
                    code: "disk-inflow-invalid",
                    message: format!("disk {k} has reversed/non-finite axial inflow {v_disk}"),
                    ranked_repairs: vec![
                        "the axial BEMT model does not admit reversed inflow".into(),
                    ],
                });
            }
            // The consumer may retain the actually evaluated wash, distinct from G(x).
            let sol = rotor_solve(k, v_disk, disks[k].omega_rad_s, *x)?;
            thrust[k] = sol.thrust_n;
            torque[k] = sol.torque_nm;
            power[k] = Some(sol.power);
            let t_term = sol.thrust_n / (2.0 * rho * disk_area);
            g[k] = momentum_disk_root(v_disk, t_term)?;
        }
        Ok((
            g,
            thrust,
            torque,
            power.map(|p| p.expect("both rotors evaluated")),
            wing,
            wing_anchor,
        ))
    };
    loop {
        let (g, thrust, torque, power, wing, wing_anchor) =
            evaluate(&x, accepted_gamma.as_deref())?;
        let r = [g[0] - x[0], g[1] - x[1]];
        // det doctrine: sqrt is IEEE-exact; hypot is platform libm and
        // broke native-vs-wasm digest identity (E6.2 measurement).
        let r02 = r[0] * r[0];
        let r12 = r[1] * r[1];
        let x02 = x[0] * x[0];
        let x12 = x[1] * x[1];
        let res = if incident.is_some() {
            (r[0].abs() / (x[0].abs() + 1e-6)).max(r[1].abs() / (x[1].abs() + 1e-6))
        } else {
            (r02 + r12).sqrt() / ((x02 + x12).sqrt() + 1e-6)
        };
        if res < spec.tol {
            if let Some(branch) = branch.as_deref_mut() {
                branch.accepted = wing_anchor;
            }
            return Ok(CoupledStep {
                w_slip: g,
                thrust_n: thrust,
                torque_nm: torque,
                power,
                wing_lift_n: wing.total_lift_n,
                gamma: wing.gamma,
                strip_alpha_eff_rad: wing.alpha_eff_rad,
                strip_velocity_mps: Some(
                    wash.iter()
                        .enumerate()
                        .map(|(s, factors)| {
                            let velocity =
                                incident.map_or(freestream, |flow| flow.strip_velocity_mps[s]);
                            if let Some(axes) = spatial {
                                return axes.velocity(velocity, *factors, x);
                            }
                            let mut velocity = velocity;
                            // `gamma` belongs to the evaluated x, not the next
                            // predictor g returned as the warm-start state.
                            velocity[0] += strip_wash(*factors, x);
                            velocity
                        })
                        .collect(),
                ),
                corrections,
                residual: res,
                spec_digest: {
                    let base = spatial.map_or_else(
                        || {
                            let digest = spec.digest();
                            if incident.is_some() {
                                digest
                            } else {
                                // Compact input now uses the same local-flow
                                // physics and physical disk admission. Preserve
                                // the already-correct explicit-flow identity.
                                hash_domain(
                                    "org.frankensim.fs-flyer.uniform-prop-coupling.v1",
                                    digest.as_bytes(),
                                )
                                .to_hex()
                            }
                        },
                        |axes| axes.digest(spec),
                    );
                    match section_frames {
                        None => base,
                        Some(frames) => {
                            let mut identity = base.into_bytes();
                            identity.extend_from_slice(&(frames.len() as u64).to_le_bytes());
                            for frame in frames {
                                for value in frame
                                    .chord_direction()
                                    .into_iter()
                                    .chain(frame.span_direction())
                                {
                                    identity.extend_from_slice(&value.to_bits().to_le_bytes());
                                }
                            }
                            hash_domain(
                                "org.frankensim.fs-flyer.section-frame-coupling.v1",
                                &identity,
                            )
                            .to_hex()
                        }
                    }
                },
            });
        }
        if !unrelaxed_done {
            // The plan's ONE UNRELAXED EVALUATION: full Picard replacement
            // x <- G(x) (omega = 1, outside the clamp by design). It is not
            // an Aitken correction and does not count against the cap; it
            // collapses the cold-start transient of this weakly-coupled
            // map before the clamped corrections mop up.
            unrelaxed_done = true;
            accepted_gamma = Some(wing.gamma);
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
        // Judge the correction against the iterate that produced it. An
        // improving floor-relaxed iterate can remain above an earlier best;
        // that is not a second growth event. The unrelaxed initialization
        // precedes this correction schedule and keeps its existing budget.
        if corrections > 0 && res > res_prev * 1.25 {
            if retry_used {
                return Err(Refusal {
                    code: "PropAirframeCouplingDidNotConverge",
                    message: format!("second residual growth at {res:e}"),
                    ranked_repairs: vec!["reject the correction; refuse, never jump".into()],
                });
            }
            // Reject the entire trial: retry FROM its accepted predecessor
            // with that predecessor's residual, never from the grown state.
            // Keep the same predecessor/history until this retry is judged.
            // Both the rejected proposal and the retry consume the existing
            // cap. After the one retry, keep the declared relaxation floor;
            // any subsequent growth still refuses instead of accepting a
            // one-way or unconverged force state.
            retry_used = true;
            let rp = r_prev.expect("a correction has an evaluated predecessor");
            for k in 0..2 {
                x[k] = x_prev[k] + spec.clamp.0 * rp[k];
            }
            omega_prev = spec.clamp.0;
            corrections += 1;
            continue;
        }
        // Vector Aitken recurrence: omega_{k+1} =
        // -omega_k * <r_k, r_{k+1}-r_k> / |r_{k+1}-r_k|^2 (exact for a
        // linear map), clamped to the candidate window; a degenerate
        // denominator resets to omega0 (spec guard).
        let omega = match r_prev {
            _ if retry_used => spec.clamp.0,
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
        accepted_gamma = Some(wing.gamma);
        x_prev = x;
        for k in 0..2 {
            x[k] += omega * r[k];
        }
        omega_prev = omega;
        r_prev = Some(r);
        res_prev = res;
        corrections += 1;
    }
}

/// Solve vi(V+vi)=t without subtracting nearly equal positive numbers.
fn momentum_disk_root(axial_mps: f64, t: f64) -> Result<f64, Refusal> {
    let disc = axial_mps * axial_mps + 4.0 * t;
    if !axial_mps.is_finite()
        || axial_mps < 0.0
        || !t.is_finite()
        || !disc.is_finite()
        || disc < 0.0
    {
        return Err(Refusal {
            code: "disk-momentum-invalid",
            message: format!("non-real/non-finite disk momentum state V={axial_mps}, t={t}"),
            ranked_repairs: vec!["use the admitted real axial actuator-disk branch".into()],
        });
    }
    let denominator = axial_mps + disc.sqrt();
    Ok(if denominator == 0.0 {
        0.0
    } else {
        2.0 * t / denominator
    })
}

#[cfg(test)]
mod tests {
    use super::momentum_disk_root;

    #[test]
    fn momentum_disk_root_retains_weak_loading_and_refuses_nonreal_states() {
        let weak = momentum_disk_root(1.0e6, 1.0e-12).unwrap();
        assert!(
            weak > 0.0,
            "subtractive root would round this induction to zero"
        );
        assert!((weak / 1.0e-18 - 1.0).abs() < 1.0e-15);
        for (v, t) in [(0.0, 4.0), (10.0, 20.0), (10.0, -20.0), (1.0e6, 1.0e-12)] {
            let vi = momentum_disk_root(v, t).unwrap();
            assert!((vi * (v + vi) / t - 1.0).abs() < 1.0e-15);
        }
        assert_eq!(momentum_disk_root(0.0, 0.0).unwrap(), 0.0);
        assert_eq!(momentum_disk_root(10.0, 0.0).unwrap(), 0.0);
        for (v, t) in [(0.0, -1.0), (10.0, -26.0), (f64::INFINITY, 1.0)] {
            assert_eq!(
                momentum_disk_root(v, t).unwrap_err().code,
                "disk-momentum-invalid"
            );
        }
    }
}
