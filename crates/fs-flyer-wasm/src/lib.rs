//! fs-flyer-wasm — the browser (WASM) surface for the Wright Flyer
//! simulation. Layer: L6. Bead frankensim-wf-root-guzez.1.3 (E0.3).
//!
//! This crate is the ONLY boundary between the Wright Flyer presentation
//! plane (three.js) and the simulation plane (FrankenSim physics). It follows
//! the proven `fs-wasm` pattern: its own decoupled workspace, `cdylib` +
//! `rlib`, wasm-bindgen confined to wasm32, and the asupersync canonical
//! browser profile pinned so the fs-exec cone builds on wasm32.
//!
//! v0 scope (the E0.3 scaffold): the HELLO KERNEL — a deterministic free
//! rigid-body spin integrated by `fs_time::lie::rigid_body_step` — plus the
//! two contracts every later API entry inherits:
//!
//! - **Typed-refusal JS envelope.** Every fallible entry returns a JSON
//!   envelope `{"ok": ...}` or `{"refusal": {"code", "message",
//!   "ranked_repairs"}}`. Nothing is silently clamped and nothing traps.
//! - **Determinism digest.** The state trajectory is content-hashed with
//!   `fs-blake3` under a versioned domain, so bit-identical replays are a
//!   checkable claim from the very first kernel, native AND wasm.
//!
//! No-claims: the hello kernel is a torque-free rigid body — it makes no
//! aerodynamic, historical, or Flyer-specific claim. It exists to prove the
//! toolchain, the envelope, and the determinism discipline end to end.

use fs_blake3::hash_domain;
use fs_time::lie::rigid_body_step;

pub mod archive;
pub mod engine;
pub mod fieldlease;
pub mod ring;

/// A SPREAD of pinned build-ups, bit-dumped (E6.2 lane bisection —
/// the Newton walks these input classes).
#[must_use]
pub fn buildup_spread_json() -> String {
    use fs_flyer::aircraft::wright_openloop_v1;
    let design = wright_openloop_v1();
    let probes = [
        [13.05, 0.06, 0.1, 45.0],
        [13.0, 0.062, 0.1, 45.0],
        [13.0, 0.06, 0.102, 45.0],
        [13.0, 0.06, 0.1, 45.25],
        [14.3, 0.041, 0.137, 50.4],
        [11.2, 0.18, 0.05, 60.0],
    ];
    let mut out = String::from("{");
    for (i, p) in probes.iter().enumerate() {
        match design.force_buildup(p[0], p[1], p[2], p[3], 0.0, 1.294) {
            Ok(b) => {
                out.push_str(&format!(
                    "\"p{i}\":[\"{}\",\"{}\",\"{}\",\"{}\",{}],",
                    b.force_n[0].to_bits(),
                    b.force_n[2].to_bits(),
                    b.moment_y_nm.to_bits(),
                    b.torque_imbalance_nm.to_bits(),
                    b.coupled.corrections
                ));
            }
            Err(e) => out.push_str(&format!("\"p{i}\":\"{}\",", e.code)),
        }
    }
    out.pop();
    out.push('}');
    out
}

/// Pinned-input BEMT bit probe (guzez.7.2.1 lane bisection, level 3):
/// `bemt_solve` at the canonical rotor/rho/omega across a v_disk sweep
/// bracketing the diverging p0 spread state. If thrust bits differ at
/// PINNED inputs the divergence lives inside fs-airscrew; if not, it
/// lives upstream in the disk-inflow projection.
#[must_use]
pub fn bemt_probe_json() -> String {
    use fs_airscrew::bemt_solve;
    use fs_flyer::aircraft::wright_rotor_v1;
    let rotor = wright_rotor_v1();
    let sweep = [12.9, 13.0, 13.05, 13.1, 13.2, 14.05];
    let mut out = String::from("{");
    for (i, v) in sweep.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        match bemt_solve(&rotor, 1.294, *v, 45.0) {
            Ok(s) => out.push_str(&format!(
                "\"v{i}\":[\"{}\",\"{}\"]",
                s.thrust_n.to_bits(),
                s.torque_nm.to_bits()
            )),
            Err(e) => out.push_str(&format!("\"v{i}\":\"{}\"", e.code)),
        }
    }
    out.push('}');
    out
}

/// Bit-exact trim iterate trace at the canonical start (guzez.7.2.1
/// lane bisection): the first differing iterate across lanes localizes
/// where cross-lane divergence enters the trim search.
#[must_use]
pub fn trim_trace_json() -> String {
    use fs_flyer::aircraft::wright_openloop_v1;
    let design = wright_openloop_v1();
    let mut trace = Vec::new();
    let outcome = design.trim_traced(1.294, [13.0, 0.06, 0.1, 45.0], &mut trace);
    let mut out = String::from("{\"iters\":[");
    for (i, x) in trace.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        out.push_str(&format!(
            "[\"{}\",\"{}\",\"{}\",\"{}\"]",
            x[0], x[1], x[2], x[3]
        ));
    }
    match outcome {
        Ok(t) => out.push_str(&format!(
            "],\"alpha\":\"{}\",\"v\":\"{}\",\"dc\":\"{}\",\"omega\":\"{}\"}}",
            t.alpha_rad.to_bits(),
            t.v_mps.to_bits(),
            t.delta_canard_rad.to_bits(),
            t.omega_prop_rad_s.to_bits(),
        )),
        Err(e) => out.push_str(&format!("],\"refusal\":\"{}\"}}", e.code)),
    }
    out
}

/// One pinned force build-up, bit-dumped (E6.2 lane bisection).
#[must_use]
pub fn buildup_probe_json() -> String {
    use fs_flyer::aircraft::wright_openloop_v1;
    let design = wright_openloop_v1();
    match design.force_buildup(13.0, 0.06, 0.1, 45.0, 0.0, 1.294) {
        Ok(b) => format!(
            "{{\"fx\":\"{}\",\"fz\":\"{}\",\"my\":\"{}\",\"lift\":\"{}\",\"thrust\":\"{}\",\"wslip\":\"{}\",\"corr\":{},\"res\":\"{}\"}}",
            b.force_n[0].to_bits(),
            b.force_n[2].to_bits(),
            b.moment_y_nm.to_bits(),
            b.lift_n.to_bits(),
            b.thrust_n[0].to_bits(),
            b.coupled.w_slip[0].to_bits(),
            b.coupled.corrections,
            b.coupled.residual.to_bits(),
        ),
        Err(e) => format!("{{\"refusal\":\"{}\"}}", e.code),
    }
}

/// Determinism probe payload (shared by native tests and the wasm
/// export): bit patterns of the det kernel + fma on this platform.
#[must_use]
pub fn det_probe_json() -> String {
    use fs_math::det;
    let fma = 0.1f64.mul_add(0.2, 0.3);
    format!(
        "{{\"exp\":\"{}\",\"sin\":\"{}\",\"ln\":\"{}\",\"atan2\":\"{}\",\"cos\":\"{}\",\"fma\":\"{}\"}}",
        det::exp(0.1).to_bits(),
        det::sin(0.3).to_bits(),
        det::ln(2.5).to_bits(),
        det::atan2(1.0, 2.0).to_bits(),
        det::cos(0.7).to_bits(),
        fma.to_bits(),
    )
}

/// Versioned identity domain for hello-kernel trajectory digests.
pub const HELLO_DIGEST_DOMAIN: &str = "org.frankensim.fs-flyer-wasm.hello-trajectory.v1";

/// Hard cap on requested integration steps (typed refusal above; the cap is
/// asserted at cap AND cap+1 by the battery per workspace law).
pub const MAX_HELLO_STEPS: u32 = 1_000_000;

/// Admitted timestep domain [s] (open interval refusals outside).
pub const MIN_DT_S: f64 = 1.0e-9;
pub const MAX_DT_S: f64 = 1.0;

/// A typed refusal crossing the JS boundary.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Refusal {
    /// Stable machine-readable code.
    pub code: &'static str,
    /// Human-readable diagnosis.
    pub message: String,
    /// Ranked repairs, most likely fix first.
    pub ranked_repairs: Vec<String>,
}

/// Final state of a hello-kernel run.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HelloState {
    /// Unit quaternion (w, x, y, z), body→world.
    pub quaternion: [f64; 4],
    /// Body-frame angular velocity [rad/s].
    pub omega_body: [f64; 3],
    /// Steps actually integrated.
    pub steps: u32,
}

/// Deterministic free rigid-body spin: principal inertia `inertia_kg_m2`,
/// initial unit quaternion `q0`, body angular velocity `omega0`, fixed
/// timestep `dt_s`, `steps` steps of the CG2 Lie-group integrator.
///
/// # Errors
/// Typed refusals for non-finite inputs, non-positive inertia, a non-unit
/// quaternion, an out-of-domain timestep, or a step count above
/// [`MAX_HELLO_STEPS`].
pub fn hello_spin(
    inertia_kg_m2: [f64; 3],
    q0: [f64; 4],
    omega0: [f64; 3],
    dt_s: f64,
    steps: u32,
) -> Result<HelloState, Refusal> {
    let finite = inertia_kg_m2.iter().all(|v| v.is_finite())
        && q0.iter().all(|v| v.is_finite())
        && omega0.iter().all(|v| v.is_finite())
        && dt_s.is_finite();
    if !finite {
        return Err(Refusal {
            code: "non-finite-input",
            message: "every hello-kernel input must be finite".into(),
            ranked_repairs: vec![
                "replace NaN/inf fields with finite values".into(),
                "check JS-side unit conversions for overflow".into(),
            ],
        });
    }
    if inertia_kg_m2.iter().any(|&i| i <= 0.0) {
        return Err(Refusal {
            code: "non-positive-inertia",
            message: format!("principal inertia must be positive, got {inertia_kg_m2:?}"),
            ranked_repairs: vec!["supply a physically valid principal inertia triple".into()],
        });
    }
    let norm2: f64 = q0.iter().map(|v| v * v).sum();
    if (norm2 - 1.0).abs() > 1.0e-9 {
        return Err(Refusal {
            code: "non-unit-quaternion",
            message: format!("q0 must be unit-norm (|q|² − 1 = {:.3e})", norm2 - 1.0),
            ranked_repairs: vec![
                "normalize the quaternion before the call".into(),
                "pass identity [1,0,0,0] if no initial attitude is intended".into(),
            ],
        });
    }
    if !(MIN_DT_S..=MAX_DT_S).contains(&dt_s) {
        return Err(Refusal {
            code: "timestep-outside-domain",
            message: format!("dt_s = {dt_s:e} outside admitted [{MIN_DT_S:e}, {MAX_DT_S:e}]"),
            ranked_repairs: vec!["use a fixed timestep inside the admitted domain".into()],
        });
    }
    if steps > MAX_HELLO_STEPS {
        return Err(Refusal {
            code: "step-budget-exceeded",
            message: format!("steps = {steps} exceeds the admitted cap {MAX_HELLO_STEPS}"),
            ranked_repairs: vec![
                format!("request at most {MAX_HELLO_STEPS} steps per call"),
                "advance long runs in bounded chunks from JS".into(),
            ],
        });
    }
    let (mut q, mut w) = (q0, omega0);
    for _ in 0..steps {
        (q, w) = rigid_body_step(q, w, inertia_kg_m2, dt_s);
    }
    Ok(HelloState {
        quaternion: q,
        omega_body: w,
        steps,
    })
}

/// Content digest of a hello-kernel trajectory: every intermediate state's
/// exact bits, hashed under [`HELLO_DIGEST_DOMAIN`]. Bit-identical runs —
/// including native-vs-wasm under the det:: doctrine — produce identical
/// hex digests; this is the seed of the six-lane golden program (E6.2).
///
/// # Errors
/// The same typed refusals as [`hello_spin`].
pub fn hello_digest(
    inertia_kg_m2: [f64; 3],
    q0: [f64; 4],
    omega0: [f64; 3],
    dt_s: f64,
    steps: u32,
) -> Result<String, Refusal> {
    // Admission is shared with hello_spin by running step zero's checks.
    hello_spin(inertia_kg_m2, q0, omega0, dt_s, 0)?;
    let mut payload = Vec::with_capacity(8 * (7 * steps as usize + 11));
    for value in inertia_kg_m2
        .iter()
        .chain(q0.iter())
        .chain(omega0.iter())
        .chain([dt_s].iter())
    {
        payload.extend_from_slice(&value.to_bits().to_le_bytes());
    }
    payload.extend_from_slice(&steps.to_le_bytes());
    let (mut q, mut w) = (q0, omega0);
    for _ in 0..steps {
        (q, w) = rigid_body_step(q, w, inertia_kg_m2, dt_s);
        for value in q.iter().chain(w.iter()) {
            payload.extend_from_slice(&value.to_bits().to_le_bytes());
        }
    }
    Ok(hash_domain(HELLO_DIGEST_DOMAIN, &payload).to_hex())
}

// ---------------------------------------------------------------------------
// JSON envelope (manual, dependency-free: the boundary speaks documented SI
// doubles; floats are serialized exactly via shortest-roundtrip formatting).
// ---------------------------------------------------------------------------

pub(crate) fn json_escape(value: &str) -> String {
    let mut out = String::with_capacity(value.len());
    for ch in value.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if c.is_control() => {
                use core::fmt::Write as _;
                let _ = write!(out, "\\u{:04x}", c as u32);
            }
            c => out.push(c),
        }
    }
    out
}

/// Render a refusal as the JS envelope.
#[must_use]
pub fn refusal_envelope(refusal: &Refusal) -> String {
    let repairs = refusal
        .ranked_repairs
        .iter()
        .map(|r| format!("\"{}\"", json_escape(r)))
        .collect::<Vec<_>>()
        .join(",");
    format!(
        "{{\"refusal\":{{\"code\":\"{}\",\"message\":\"{}\",\"ranked_repairs\":[{}]}}}}",
        json_escape(refusal.code),
        json_escape(&refusal.message),
        repairs
    )
}

/// Render a hello-kernel success as the JS envelope.
#[must_use]
pub fn hello_envelope(state: &HelloState) -> String {
    format!(
        "{{\"ok\":{{\"quaternion\":[{},{},{},{}],\"omega_body\":[{},{},{}],\"steps\":{}}}}}",
        state.quaternion[0],
        state.quaternion[1],
        state.quaternion[2],
        state.quaternion[3],
        state.omega_body[0],
        state.omega_body[1],
        state.omega_body[2],
        state.steps
    )
}

/// One CG2 step with an external body torque [N·m].
/// Split: apply `I⁻¹ τ Δt`, then the torque-free `rigid_body_step`.
/// This is the aero kernel; `hello_spin` remains torque-free.
///
/// # Errors
/// Same admissions as [`hello_spin`], plus non-finite torque.
pub fn aero_step(
    inertia_kg_m2: [f64; 3],
    q0: [f64; 4],
    omega0: [f64; 3],
    torque_n_m: [f64; 3],
    dt_s: f64,
    steps: u32,
) -> Result<HelloState, Refusal> {
    hello_spin(inertia_kg_m2, q0, omega0, dt_s, 0)?;
    if !torque_n_m.iter().all(|v| v.is_finite()) {
        return Err(Refusal {
            code: "non-finite-torque",
            message: "aero torque must be finite".into(),
            ranked_repairs: vec!["replace NaN/inf torque with a finite body wrench".into()],
        });
    }
    let (mut q, mut w) = (q0, omega0);
    for _ in 0..steps {
        w[0] += torque_n_m[0] / inertia_kg_m2[0] * dt_s;
        w[1] += torque_n_m[1] / inertia_kg_m2[1] * dt_s;
        w[2] += torque_n_m[2] / inertia_kg_m2[2] * dt_s;
        (q, w) = rigid_body_step(q, w, inertia_kg_m2, dt_s);
    }
    Ok(HelloState {
        quaternion: q,
        omega_body: w,
        steps,
    })
}

#[must_use]
pub fn aero_step_json(
    inertia_kg_m2: [f64; 3],
    q0: [f64; 4],
    omega0: [f64; 3],
    torque_n_m: [f64; 3],
    dt_s: f64,
    steps: u32,
) -> String {
    match aero_step(inertia_kg_m2, q0, omega0, torque_n_m, dt_s, steps) {
        Ok(state) => hello_envelope(&state),
        Err(refusal) => refusal_envelope(&refusal),
    }
}

/// Envelope-producing entry shared by native tests and the wasm boundary.
#[must_use]
pub fn hello_spin_json(
    inertia_kg_m2: [f64; 3],
    q0: [f64; 4],
    omega0: [f64; 3],
    dt_s: f64,
    steps: u32,
) -> String {
    match hello_spin(inertia_kg_m2, q0, omega0, dt_s, steps) {
        Ok(state) => hello_envelope(&state),
        Err(refusal) => refusal_envelope(&refusal),
    }
}

// ---------------------------------------------------------------------------
// wasm32-only JS boundary.
// ---------------------------------------------------------------------------
#[cfg(target_arch = "wasm32")]
mod js {
    use core::cell::RefCell;
    use wasm_bindgen::prelude::wasm_bindgen;

    thread_local! {
        static ENGINE: RefCell<super::engine::EngineSlot> =
            RefCell::new(super::engine::EngineSlot::default());
    }

    /// Initialize the Wright Flyer lifecycle engine (E5.1). Replaces
    /// any prior run in this worker. mode: 0=fixed, 1=historical
    /// (member selects the registered pilot), 2=human. Returns the
    /// init envelope (run_intent_id, tick0_digest, trim) or a typed
    /// refusal envelope.
    #[wasm_bindgen]
    #[must_use]
    #[allow(clippy::too_many_arguments)]
    pub fn flyer_engine_init(
        seed: u64,
        rho_kg_m3: f64,
        headwind_mps: f64,
        mode: u32,
        member: u32,
        rail_length_m: f64,
        max_ticks: u64,
        assist: bool,
        catapult: bool,
    ) -> String {
        ENGINE.with(|e| {
            e.borrow_mut().init(
                seed,
                rho_kg_m3,
                headwind_mps,
                mode,
                member,
                rail_length_m,
                max_ticks,
                assist,
                catapult,
            )
        })
    }

    /// One 120 Hz engine step. `has_input` gates whether (lever, warp)
    /// is a ControlInput (Human mode requires it every tick).
    #[wasm_bindgen]
    #[must_use]
    pub fn flyer_engine_step(has_input: bool, lever_force_n: f64, warp_cmd_rad: f64) -> String {
        ENGINE.with(|e| e.borrow_mut().step(has_input, lever_force_n, warp_cmd_rad))
    }

    /// The chained lifecycle digest envelope.
    #[wasm_bindgen]
    #[must_use]
    pub fn flyer_engine_digest() -> String {
        ENGINE.with(|e| e.borrow().digest())
    }

    /// Determinism probe (E6.2 six-lane diagnostics; doc-hidden class):
    /// bit patterns of the det:: kernel + fma on this platform.
    #[wasm_bindgen]
    #[must_use]
    pub fn flyer_det_probe() -> String {
        super::det_probe_json()
    }

    /// Pinned force-buildup bit probe (E6.2 lane bisection).
    #[wasm_bindgen]
    #[must_use]
    pub fn flyer_buildup_probe() -> String {
        super::buildup_probe_json()
    }

    /// The startup determinism self-test (per-lane golden; the app
    /// shows the failure badge on a refusal envelope).
    #[wasm_bindgen]
    #[must_use]
    pub fn flyer_selftest() -> String {
        super::engine::selftest()
    }

    /// Trim iterate bit-trace (guzez.7.2.1 lane bisection).
    #[wasm_bindgen]
    #[must_use]
    pub fn flyer_trim_trace() -> String {
        crate::trim_trace_json()
    }

    /// Pinned-input BEMT bit probe (guzez.7.2.1 lane bisection).
    #[wasm_bindgen]
    #[must_use]
    pub fn flyer_bemt_probe() -> String {
        crate::bemt_probe_json()
    }

    /// Spread probe (E6.2 lane bisection).
    #[wasm_bindgen]
    #[must_use]
    pub fn flyer_buildup_spread() -> String {
        super::buildup_spread_json()
    }

    /// E7.1-ii field-lease self-test: ring -> lease -> §5.5 sample ->
    /// bounded JSON (or a typed refusal envelope).
    #[wasm_bindgen]
    #[must_use]
    pub fn flyer_field_selftest() -> String {
        super::fieldlease::field_selftest_json()
    }

    /// Deterministic free rigid-body spin; returns the typed JSON envelope.
    #[wasm_bindgen]
    #[must_use]
    pub fn flyer_hello_spin(
        ixx: f64,
        iyy: f64,
        izz: f64,
        qw: f64,
        qx: f64,
        qy: f64,
        qz: f64,
        wx: f64,
        wy: f64,
        wz: f64,
        dt_s: f64,
        steps: u32,
    ) -> String {
        super::hello_spin_json([ixx, iyy, izz], [qw, qx, qy, qz], [wx, wy, wz], dt_s, steps)
    }

    /// CG2 attitude step with a body-frame torque [N·m].
    #[wasm_bindgen]
    #[must_use]
    pub fn flyer_aero_step(
        ixx: f64,
        iyy: f64,
        izz: f64,
        qw: f64,
        qx: f64,
        qy: f64,
        qz: f64,
        wx: f64,
        wy: f64,
        wz: f64,
        tx: f64,
        ty: f64,
        tz: f64,
        dt_s: f64,
        steps: u32,
    ) -> String {
        super::aero_step_json(
            [ixx, iyy, izz],
            [qw, qx, qy, qz],
            [wx, wy, wz],
            [tx, ty, tz],
            dt_s,
            steps,
        )
    }

    /// Trajectory content digest (hex) or the refusal envelope.
    #[wasm_bindgen]
    #[must_use]
    pub fn flyer_hello_digest(
        ixx: f64,
        iyy: f64,
        izz: f64,
        qw: f64,
        qx: f64,
        qy: f64,
        qz: f64,
        wx: f64,
        wy: f64,
        wz: f64,
        dt_s: f64,
        steps: u32,
    ) -> String {
        match super::hello_digest([ixx, iyy, izz], [qw, qx, qy, qz], [wx, wy, wz], dt_s, steps) {
            Ok(hex) => format!("{{\"ok\":\"{hex}\"}}"),
            Err(refusal) => super::refusal_envelope(&refusal),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const INERTIA: [f64; 3] = [0.9, 1.1, 1.7];
    const IDENTITY_Q: [f64; 4] = [1.0, 0.0, 0.0, 0.0];
    const OMEGA: [f64; 3] = [0.3, -1.2, 2.1];
    const DT: f64 = 1.0 / 120.0;

    #[test]
    fn aero_step_yaws_when_torque_is_applied() {
        let free = hello_spin(INERTIA, IDENTITY_Q, [0.0, 0.0, 0.0], DT, 120).unwrap();
        let yawed = aero_step(
            INERTIA,
            IDENTITY_Q,
            [0.0, 0.0, 0.0],
            [0.0, 4.0, 0.0],
            DT,
            120,
        )
        .unwrap();
        assert!(
            (yawed.omega_body[1] - free.omega_body[1]).abs() > 1.0e-6,
            "yaw torque must change ω_y"
        );
        let n2: f64 = yawed.quaternion.iter().map(|v| v * v).sum();
        assert!((n2 - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn hello_spin_is_deterministic_and_norm_preserving() {
        let a = hello_spin(INERTIA, IDENTITY_Q, OMEGA, DT, 4_800).unwrap();
        let b = hello_spin(INERTIA, IDENTITY_Q, OMEGA, DT, 4_800).unwrap();
        // Bit identity across two runs in one process.
        for i in 0..4 {
            assert_eq!(a.quaternion[i].to_bits(), b.quaternion[i].to_bits());
        }
        for i in 0..3 {
            assert_eq!(a.omega_body[i].to_bits(), b.omega_body[i].to_bits());
        }
        // Lie-group update preserves the unit norm by construction.
        let norm2: f64 = a.quaternion.iter().map(|v| v * v).sum();
        assert!((norm2 - 1.0).abs() < 1.0e-12, "norm² drift {norm2:e}");
    }

    #[test]
    fn hello_digest_is_stable_and_input_sensitive() {
        let base = hello_digest(INERTIA, IDENTITY_Q, OMEGA, DT, 480).unwrap();
        let again = hello_digest(INERTIA, IDENTITY_Q, OMEGA, DT, 480).unwrap();
        assert_eq!(base, again, "identical runs must digest identically");
        // GOLDEN: pinned trajectory digest for the canonical hello scenario.
        // Native aarch64/x86 and wasm-in-node must all reproduce this exact
        // value (the E6.2 six-lane program consumes it). Golden-bump protocol
        // applies to any change.
        assert_eq!(
            base, "d1ee8e689eebfe8ca4330cab2d256968534962ceaf034dacf788fed1838323f9",
            "canonical hello digest moved — determinism regression or an \
             intentional kernel change requiring the golden-bump protocol"
        );
        let mut omega = OMEGA;
        omega[2] += 1.0e-12;
        let perturbed = hello_digest(INERTIA, IDENTITY_Q, omega, DT, 480).unwrap();
        assert_ne!(base, perturbed, "a 1e-12 input change must move the digest");
    }

    #[test]
    fn refusals_are_typed_ranked_and_cap_exact() {
        // Cap AND cap+1 (workspace law).
        assert!(hello_spin(INERTIA, IDENTITY_Q, OMEGA, DT, MAX_HELLO_STEPS - 1).is_ok());
        let at_cap = hello_spin(INERTIA, IDENTITY_Q, OMEGA, DT, MAX_HELLO_STEPS);
        assert!(at_cap.is_ok(), "the cap itself is admitted");
        let over = hello_spin(INERTIA, IDENTITY_Q, OMEGA, DT, MAX_HELLO_STEPS + 1).unwrap_err();
        assert_eq!(over.code, "step-budget-exceeded");
        assert!(!over.ranked_repairs.is_empty());

        let nan = hello_spin([f64::NAN, 1.0, 1.0], IDENTITY_Q, OMEGA, DT, 1).unwrap_err();
        assert_eq!(nan.code, "non-finite-input");
        let neg = hello_spin([-1.0, 1.0, 1.0], IDENTITY_Q, OMEGA, DT, 1).unwrap_err();
        assert_eq!(neg.code, "non-positive-inertia");
        let skew = hello_spin(INERTIA, [1.0, 0.5, 0.0, 0.0], OMEGA, DT, 1).unwrap_err();
        assert_eq!(skew.code, "non-unit-quaternion");
        let dt = hello_spin(INERTIA, IDENTITY_Q, OMEGA, 2.0, 1).unwrap_err();
        assert_eq!(dt.code, "timestep-outside-domain");
    }

    #[test]
    fn json_envelope_is_wellformed_for_ok_and_refusal() {
        let ok = hello_spin_json(INERTIA, IDENTITY_Q, OMEGA, DT, 10);
        assert!(ok.starts_with("{\"ok\":{"), "{ok}");
        assert!(ok.contains("\"steps\":10"));
        let bad = hello_spin_json(INERTIA, IDENTITY_Q, OMEGA, -1.0, 10);
        assert!(bad.starts_with("{\"refusal\":{"), "{bad}");
        assert!(bad.contains("\"code\":\"timestep-outside-domain\""));
        assert!(bad.contains("ranked_repairs"));
        // Hostile message content must escape cleanly.
        let hostile = refusal_envelope(&Refusal {
            code: "x",
            message: "quote\" backslash\\ newline\n".into(),
            ranked_repairs: vec!["tab\there".into()],
        });
        assert!(hostile.contains("quote\\\" backslash\\\\ newline\\n"));
        assert!(hostile.contains("tab\\there"));
    }
}
