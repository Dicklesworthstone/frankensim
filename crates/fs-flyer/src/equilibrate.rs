//! HeldOnRailEquilibrated closure (beads wf-root-guzez.5.17.1/.2,
//! E4.6d). Plan Round-5: the run begins from a CLOSED initialization —
//! the OU atmosphere drawn stationary at the canonical −3840-tick
//! anchor and advanced to tick 0, the aircraft at its fixed-throttle
//! trim (aero + rotor torque balance = the pinned observable envelope),
//! and the canard mechanism SETTLED onto the trim deflection by a
//! deterministic preroll. The resulting `Tick0State` is MODE-COMPLETE
//! and its digest is frozen BEFORE RunIntentId — structurally: no
//! intent field exists on the type, so intent cannot leak into the
//! initial-state identity.
//!
//! Branch semantics (leaf ii): alternate starts WITHIN a branch
//! converge to the canonical envelope inside DECLARED numerical closure
//! bands (never bitwise — the canonical start defines the digest);
//! consuming a tick-0 state under the wrong branch is the typed
//! `PrelaunchBranchMismatch` refusal.

use crate::Refusal;
use crate::aircraft::{OpenLoopDesign, TrimResult};
use crate::canardmech::{CanardMechanism, MechState};
use fs_atmo::ou::{OuMode, STATIONARY_ANCHOR_TICK, StationaryOuPath};
use fs_blake3::hash_domain;

/// Preroll budget cap [ticks].
pub const MAX_PREROLL_TICKS: u32 = 28_800; // 4 minutes at 120 Hz

/// Mechanism settle tolerance [rad] on |δ − δ_trim|.
pub const SETTLE_TOL_RAD: f64 = 1.0e-4;

/// Mechanism settle tolerance [rad/s] on |δ̇|.
pub const SETTLE_RATE_TOL: f64 = 1.0e-3;

/// The initialization branch (enters the digest).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PrelaunchBranch {
    /// The historical Dec-17 closure: held on the rail, equilibrated.
    HeldOnRailEquilibrated,
    /// Counterfactual in-air start (calibration/testing family).
    FreeAirEquilibrated,
}

impl PrelaunchBranch {
    fn discriminant(self) -> u8 {
        match self {
            PrelaunchBranch::HeldOnRailEquilibrated => 0,
            PrelaunchBranch::FreeAirEquilibrated => 1,
        }
    }
}

/// The equilibration specification.
#[derive(Clone, Debug, PartialEq)]
pub struct EquilibrationSpec {
    /// Branch identity.
    pub branch: PrelaunchBranch,
    /// Study seed (OU stream identity).
    pub seed: u64,
    /// OU gust modes (sigma, tau pairs realized upstream).
    pub ou_modes: Vec<OuMode>,
    /// Requested prehistory anchor [tick] — MUST equal the canonical
    /// anchor (canonical-prehistory invariance made executable).
    pub anchor_tick: i64,
    /// Mechanism preroll budget [ticks at 120 Hz].
    pub preroll_ticks: u32,
    /// Air density [kg/m³].
    pub rho_kg_m3: f64,
    /// Trim start point (V, α, δc, ω).
    pub trim_start: [f64; 4],
}

/// The mode-complete tick-0 state. NO intent field exists — the digest
/// is frozen before any RunIntentId can touch it (plan law).
#[derive(Clone, Debug, PartialEq)]
pub struct Tick0State {
    /// Branch this state was equilibrated under.
    pub branch: PrelaunchBranch,
    /// The trim (pinned observable envelope).
    pub trim: TrimResult,
    /// The settled mechanism state.
    pub mech: MechState,
    /// OU amplitudes at tick 0.
    pub ou_amplitudes: Vec<f64>,
    /// Preroll ticks actually used to settle.
    pub settle_ticks: u32,
    /// The frozen mode-complete digest.
    pub digest: String,
}

/// Equilibrate to tick 0.
///
/// # Errors
/// `prelaunch-anchor-invalid` (non-canonical anchor — refused at ±1
/// tick around the canonical value); `prelaunch-spec-invalid`
/// (preroll cap at cap AND cap+1); `prelaunch-not-settled` (the
/// mechanism failed to reach the trim deflection inside the budget —
/// typed, never a silent acceptance); trim/OU refusals pass through.
pub fn equilibrate(
    design: &OpenLoopDesign,
    mech: &CanardMechanism,
    spec: &EquilibrationSpec,
) -> Result<Tick0State, Refusal> {
    if spec.anchor_tick != STATIONARY_ANCHOR_TICK {
        return Err(Refusal {
            code: "prelaunch-anchor-invalid",
            message: format!(
                "anchor {} is not the canonical {STATIONARY_ANCHOR_TICK} — prehistory-window \
                 choice would change the realization (canonical-prehistory law)",
                spec.anchor_tick
            ),
            ranked_repairs: vec!["use the canonical anchor".into()],
        });
    }
    if spec.preroll_ticks == 0 || spec.preroll_ticks > MAX_PREROLL_TICKS {
        return Err(Refusal {
            code: "prelaunch-spec-invalid",
            message: format!(
                "preroll {} ticks outside [1, {MAX_PREROLL_TICKS}]",
                spec.preroll_ticks
            ),
            ranked_repairs: vec!["8 s (960 ticks) is the settle-budget class".into()],
        });
    }
    // 1. OU prehistory: stationary at the canonical anchor, advanced to 0.
    let mut ou = StationaryOuPath::stationary_at_anchor(spec.seed, spec.ou_modes.clone())
        .map_err(map_atmo)?;
    ou.advance_to(0).map_err(map_atmo)?;
    // 2. Aero/rotor equilibration: the fixed-throttle trim.
    let trim = design.trim(spec.rho_kg_m3, spec.trim_start)?;
    // 3. Mechanism preroll from rest onto the trim deflection.
    let (mech_state, settle_ticks) = settle_mech(
        mech,
        trim.delta_canard_rad,
        MechState {
            delta_rad: 0.0,
            rate_rad_s: 0.0,
        },
        spec.preroll_ticks,
    )?;
    let digest = tick0_digest(spec.branch, &trim, &mech_state, ou.amplitudes(), spec.seed);
    Ok(Tick0State {
        branch: spec.branch,
        trim,
        mech: mech_state,
        ou_amplitudes: ou.amplitudes().to_vec(),
        settle_ticks,
        digest,
    })
}

/// Settle the mechanism onto a target deflection with a deterministic
/// proprioceptive hold force (shared by canonical and alternate starts).
pub(crate) fn settle_mech(
    mech: &CanardMechanism,
    target_rad: f64,
    start: MechState,
    budget_ticks: u32,
) -> Result<(MechState, u32), Refusal> {
    let dt = 1.0 / 120.0;
    let mut st = start;
    for tick in 0..budget_ticks {
        if (st.delta_rad - target_rad).abs() < SETTLE_TOL_RAD
            && st.rate_rad_s.abs() < SETTLE_RATE_TOL
        {
            return Ok((st, tick));
        }
        let force = ((3000.0 * (target_rad - st.delta_rad) - 180.0 * st.rate_rad_s)
            / mech.lever_gain_nm_per_n)
            .clamp(-220.0, 220.0);
        st = mech.step(st, 0.0, force, dt)?.0;
    }
    Err(Refusal {
        code: "prelaunch-not-settled",
        message: format!(
            "mechanism at δ {} rad (target {target_rad}) rate {} after {budget_ticks} ticks",
            st.delta_rad, st.rate_rad_s
        ),
        ranked_repairs: vec![
            "raise the preroll budget".into(),
            "an unsettled prelaunch is a refusal, never a silent start".into(),
        ],
    })
}

/// The mode-complete tick-0 digest (frozen BEFORE RunIntentId; the
/// input set is exactly the state, never intent).
fn tick0_digest(
    branch: PrelaunchBranch,
    trim: &TrimResult,
    mech: &MechState,
    ou_amplitudes: &[f64],
    seed: u64,
) -> String {
    let mut p = Vec::new();
    p.push(branch.discriminant());
    p.extend_from_slice(&seed.to_le_bytes());
    for v in [
        trim.v_mps,
        trim.alpha_rad,
        trim.delta_canard_rad,
        trim.omega_prop_rad_s,
        mech.delta_rad,
        mech.rate_rad_s,
    ] {
        p.extend_from_slice(&v.to_bits().to_le_bytes());
    }
    for a in ou_amplitudes {
        p.extend_from_slice(&a.to_bits().to_le_bytes());
    }
    hash_domain("org.frankensim.fs-flyer.tick0.v1", &p).to_hex()
}

/// Consume a tick-0 state under a declared branch: the typed
/// `PrelaunchBranchMismatch` refusal (leaf ii, executed both ways).
///
/// # Errors
/// `PrelaunchBranchMismatch`.
pub fn admit_for_branch(state: &Tick0State, expected: PrelaunchBranch) -> Result<(), Refusal> {
    if state.branch != expected {
        return Err(Refusal {
            code: "PrelaunchBranchMismatch",
            message: format!(
                "tick-0 state was equilibrated under {:?} but is being consumed under {expected:?}",
                state.branch
            ),
            ranked_repairs: vec![
                "re-equilibrate under the intended branch; cross-branch reuse silently \
                 changes the physics"
                    .into(),
            ],
        });
    }
    Ok(())
}

fn map_atmo(e: fs_atmo::Refusal) -> Refusal {
    Refusal {
        code: e.code,
        message: e.message,
        ranked_repairs: e.ranked_repairs,
    }
}

/// Battery-visible settle entry (leaf ii alternate starts).
///
/// # Errors
/// `prelaunch-not-settled`.
pub fn settle_for_test(
    mech: &CanardMechanism,
    target_rad: f64,
    start: MechState,
    budget_ticks: u32,
) -> Result<(MechState, u32), Refusal> {
    settle_mech(mech, target_rad, start, budget_ticks)
}
