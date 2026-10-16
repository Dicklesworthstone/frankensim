//! Partitioned integrator + time-scale certificate (bead
//! wf-root-guzez.4.4, E3.2b). Plan §5.1.2 "Partitioned integrator":
//! exact discrete transitions for linear aero-memory states, implicit
//! midpoint for stiff hinge/control states, the energy-consistent joint
//! effective-mass solve (E3.2a), and the Lie-group rigid update (E3.2-i
//! spine) — composed in ONE DECLARED ORDER per tick:
//!
//!   1. aero-memory exact transitions (inputs held over the tick)
//!   2. added-mass assembly + JOINT (6+nc) solve  → ν̇, coupling force
//!   3. stiff hinge state: implicit midpoint (unconditionally stable)
//!   4. rigid state: the spine step (Verlet + Strang Lie rotation)
//!
//! The battery verifies the composition EXPLICITLY (a miscomposed twin
//! differs; the declared order is pinned by golden) — never assumed.
//!
//! COUPLING CONTRACT (battery-measured): the caller samples cross-coupling
//! loads/torques AND the memory inputs at the TICK MIDPOINT (a predictor
//! evaluation; the held-input exponential integrator is the midpoint rule
//! then). With
//! midpoint-sampled couplings the composition is order 2; tick-START
//! sampling degrades it to order ~1 — the battery measured exactly that
//! (1.22) before the contract was written, so the clause is enforced by
//! evidence, not convention.
//!
//! Time-scale certificate: every state declares its class — ExactTransition
//! and ImplicitMidpoint are admissible at ANY stiffness (that is their
//! point); an ExplicitResolved state REFUSES when dt/τ exceeds the cap.
//! Certification runs at admission, before any stepping.
//!
//! No-claims: event-localized ground contact is E3.4; augmented-loop
//! gain/phase margins on the REAL aero assembly are V-05b (E4.6a). This
//! module's pole-convergence receipts cover the hinge subsystem's exact
//! discrete map.

use crate::Refusal;
use crate::spine::{Loads, RigidBody, SixDofState, step};
use fs_math::det;

/// dt/τ cap for explicit states (refusals at cap AND cap+1 ulp).
pub const MAX_EXPLICIT_RATIO: f64 = 0.2;
/// State-count cap for a certificate.
pub const MAX_CERT_STATES: usize = 64;

fn refuse(code: &'static str, message: String, repair: &str) -> Refusal {
    Refusal {
        code,
        message,
        ranked_repairs: vec![repair.into()],
    }
}

/// Integrator class of one state (the ladder the certificate records).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum IntegratorClass {
    /// Linear state advanced by its exact discrete transition.
    ExactTransition,
    /// Stiff state advanced by implicit midpoint (A-stable).
    ImplicitMidpoint,
    /// State advanced explicitly — must be resolved (dt/τ ≤ cap).
    ExplicitResolved,
}

/// One row of the time-scale certificate.
#[derive(Clone, Debug, PartialEq)]
pub struct TimeScaleEntry {
    /// State name.
    pub name: &'static str,
    /// Characteristic time [s] (1/ω_n for oscillatory states).
    pub tau_s: f64,
    /// dt/τ ratio at the certified tick.
    pub ratio: f64,
    /// Integrator class.
    pub class: IntegratorClass,
}

/// The admission-time certificate: every state, its ratio, its class.
#[derive(Clone, Debug, PartialEq)]
pub struct TimeScaleCertificate {
    /// Certified timestep [s].
    pub dt_s: f64,
    /// Per-state rows.
    pub entries: Vec<TimeScaleEntry>,
}

/// Certify a state set at `dt_s`.
///
/// # Errors
/// `certificate-states-exceeded` (cap AND cap+1);
/// `timescale-invalid` (non-finite or non-positive τ/dt);
/// `stiffness-unresolved` naming the state and its ratio when an
/// EXPLICIT state exceeds [`MAX_EXPLICIT_RATIO`].
pub fn certify(
    dt_s: f64,
    states: &[(&'static str, f64, IntegratorClass)],
) -> Result<TimeScaleCertificate, Refusal> {
    if states.len() > MAX_CERT_STATES {
        return Err(refuse(
            "certificate-states-exceeded",
            format!("{} states exceed {MAX_CERT_STATES}", states.len()),
            "aggregate substates",
        ));
    }
    if !dt_s.is_finite() || dt_s <= 0.0 {
        return Err(refuse(
            "timescale-invalid",
            format!("dt {dt_s}"),
            "positive finite dt",
        ));
    }
    let mut entries = Vec::with_capacity(states.len());
    for (name, tau_s, class) in states {
        if !tau_s.is_finite() || *tau_s <= 0.0 {
            return Err(refuse(
                "timescale-invalid",
                format!("state {name}: tau {tau_s}"),
                "positive finite characteristic time",
            ));
        }
        let ratio = dt_s / tau_s;
        if *class == IntegratorClass::ExplicitResolved && ratio > MAX_EXPLICIT_RATIO {
            return Err(refuse(
                "stiffness-unresolved",
                format!(
                    "explicit state {name}: dt/tau = {ratio:.4} exceeds {MAX_EXPLICIT_RATIO} — \
                     the state is stiff at this tick"
                ),
                "move the state to ImplicitMidpoint or ExactTransition; never under-resolve",
            ));
        }
        entries.push(TimeScaleEntry {
            name,
            tau_s: *tau_s,
            ratio,
            class: *class,
        });
    }
    Ok(TimeScaleCertificate { dt_s, entries })
}

/// A linear aero-memory state ȧ = (u − a)/τ with its EXACT one-tick
/// transition for input held over the tick:
/// a⁺ = e^(−dt/τ)·a + (1 − e^(−dt/τ))·u.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MemoryState {
    /// Correlation/relaxation time [s].
    pub tau_s: f64,
    /// Current value.
    pub value: f64,
}

impl MemoryState {
    /// Exact transition (see type docs).
    #[must_use]
    pub fn advanced(self, input: f64, dt_s: f64) -> MemoryState {
        let decay = det::exp(-dt_s / self.tau_s);
        MemoryState {
            tau_s: self.tau_s,
            value: decay * self.value + (1.0 - decay) * input,
        }
    }
}

/// Stiff hinge parameters: I_eff·q̈ + c·q̇ + k·q = Q.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HingeParams {
    /// Effective inertia [kg·m²] (structural + added-mass cc block).
    pub inertia: f64,
    /// Damping [N·m·s/rad].
    pub damping: f64,
    /// Stiffness [N·m/rad].
    pub stiffness: f64,
}

/// One implicit-midpoint step of the hinge (Q held at the midpoint):
/// with q_m = (q+q⁺)/2, v_m = (v+v⁺)/2:
///   q⁺ = q + dt·v_m,  v⁺ = v + dt·(Q − c·v_m − k·q_m)/I.
/// Linear 2×2 solved in closed form — A-stable, symplectic, order 2; on
/// the undamped oscillator it conserves H = ½Iv² + ½kq² EXACTLY (a
/// Cayley map), which the battery checks to machine precision.
///
/// # Errors
/// `hinge-params-invalid` (non-finite or non-positive inertia; negative
/// damping/stiffness).
pub fn hinge_implicit_midpoint(
    p: &HingeParams,
    q: f64,
    v: f64,
    torque: f64,
    dt_s: f64,
) -> Result<(f64, f64), Refusal> {
    if !(p.inertia.is_finite() && p.damping.is_finite() && p.stiffness.is_finite())
        || p.inertia <= 0.0
        || p.damping < 0.0
        || p.stiffness < 0.0
    {
        return Err(refuse(
            "hinge-params-invalid",
            format!("{p:?}"),
            "positive inertia; non-negative damping/stiffness",
        ));
    }
    let h = 0.5 * dt_s;
    // Unknowns (dq, dv) with q⁺ = q + dq, v⁺ = v + dv:
    //   dq = dt·(v + dv/2)
    //   I·dv = dt·(Q − c·(v + dv/2) − k·(q + dq/2))
    // Substitute dq:
    //   (I + c·h + k·h·dt/2·... ) — solve the 2×2 directly.
    let a11 = 1.0;
    let a12 = -h;
    let a21 = p.stiffness * h / p.inertia;
    let a22 = 1.0 + (p.damping * h) / p.inertia;
    let b1 = dt_s * v;
    let b2 = dt_s * (torque - p.damping * v - p.stiffness * q) / p.inertia;
    let det_m = a11 * a22 - a12 * a21;
    let dq = (b1 * a22 - a12 * b2) / det_m;
    let dv = (a11 * b2 - b1 * a21) / det_m;
    Ok((q + dq, v + dv))
}

/// The partitioned per-tick state.
#[derive(Clone, Debug, PartialEq)]
pub struct PartitionedState {
    /// Rigid 6-DOF state.
    pub rigid: SixDofState,
    /// Hinge coordinate [rad] and rate [rad/s].
    pub hinge_q: f64,
    /// Hinge rate.
    pub hinge_v: f64,
    /// Aero-memory states.
    pub memory: Vec<MemoryState>,
}

/// One partitioned tick in the DECLARED composition order (module docs).
/// `memory_inputs` are the per-state inputs held over the tick;
/// `hinge_torque` and `rigid_loads` are the non-acceleration loads.
///
/// # Errors
/// Refusals from the stages ([`hinge_implicit_midpoint`], spine
/// admission); `memory-input-mismatch`.
pub fn partitioned_step(
    body: &RigidBody,
    hinge: &HingeParams,
    state: &PartitionedState,
    memory_inputs: &[f64],
    hinge_torque: f64,
    rigid_loads: Loads,
    t_s: f64,
    dt_s: f64,
) -> Result<PartitionedState, Refusal> {
    if memory_inputs.len() != state.memory.len() {
        return Err(refuse(
            "memory-input-mismatch",
            format!(
                "{} inputs vs {} states",
                memory_inputs.len(),
                state.memory.len()
            ),
            "one input per memory state",
        ));
    }
    // Stage 1: exact memory transitions.
    let memory: Vec<MemoryState> = state
        .memory
        .iter()
        .zip(memory_inputs)
        .map(|(m, u)| m.advanced(*u, dt_s))
        .collect();
    // Stage 2 is the caller-side joint solve (E3.2a) that produced the
    // effective hinge inertia and the coupled loads — its outputs enter
    // here as `hinge` (I_eff) and the load arguments.
    // Stage 3: stiff hinge, implicit midpoint.
    let (hinge_q, hinge_v) =
        hinge_implicit_midpoint(hinge, state.hinge_q, state.hinge_v, hinge_torque, dt_s)?;
    // Stage 4: rigid spine step with the tick's loads.
    let rigid = step(body, &state.rigid, t_s, dt_s, |_, _| rigid_loads)?;
    Ok(PartitionedState {
        rigid,
        hinge_q,
        hinge_v,
        memory,
    })
}
