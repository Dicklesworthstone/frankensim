//! E3.2b battery (bead wf-root-guzez.4.4) — V-05a's integrator clauses,
//! executed: exact memory transitions vs the closed form; implicit-midpoint
//! energy conservation (machine precision) + stability where the explicit
//! twin DIVERGES; augmented-pole convergence measured from the discrete
//! one-step map (order ~2 on both damping and frequency); local order of
//! the full partitioned composition (Richardson); the composition order
//! verified EXPLICITLY (miscomposed twin differs; declared order pinned by
//! golden); time-scale certificate caps at cap AND the next float.
//! Repro: cargo test -p fs-flyer --test partitioned_battery

use fs_flyer::partitioned::{
    HingeParams, IntegratorClass, MAX_EXPLICIT_RATIO, MemoryState, PartitionedState, certify,
    hinge_implicit_midpoint, partitioned_step,
};
use fs_flyer::spine::{Loads, RigidBody, SixDofState};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-flyer-e32b\",\"case\":\"{case}\",{payload}}}");
}

const DT: f64 = 1.0 / 120.0;

fn body() -> RigidBody {
    RigidBody {
        mass_kg: 340.17,
        inertia_kgm2: [1787.0, 367.4, 1820.9],
    }
}

fn rest() -> SixDofState {
    SixDofState {
        pos_m: [0.0; 3],
        vel_mps: [0.0; 3],
        quat: [1.0, 0.0, 0.0, 0.0],
        omega_body: [0.0; 3],
    }
}

#[test]
fn memory_exact_transition_matches_closed_form() {
    // a(t) = u + (a0 − u)·e^(−t/τ), sampled at tick multiples: the exact
    // transition composed n times must equal the closed form at each n.
    let (a0, u, tau) = (2.0, 0.5, 0.31);
    let mut m = MemoryState {
        tau_s: tau,
        value: a0,
    };
    for n in 1..=240 {
        m = m.advanced(u, DT);
        let exact = u + (a0 - u) * (-(f64::from(n)) * DT / tau).exp();
        assert!(
            (m.value - exact).abs() < 1e-13,
            "memory diverged from the closed form at tick {n}: {} vs {exact}",
            m.value
        );
    }
    jlog("memory-exact", &format!("\"final\":{}", m.value));
}

#[test]
fn implicit_midpoint_conserves_energy_and_survives_stiffness() {
    // Undamped oscillator at a STIFF frequency (ω·dt ≈ 8.3 — far beyond
    // explicit stability): implicit midpoint conserves H exactly (Cayley).
    let p = HingeParams {
        inertia: 0.05,
        stiffness: 50_000.0,
        damping: 0.0,
    };
    let omega = (p.stiffness / p.inertia).sqrt();
    assert!(
        omega * DT > 2.0,
        "fixture must be genuinely stiff (ω·dt = {})",
        omega * DT
    );
    let (mut q, mut v) = (0.02, 0.0);
    let h0 = 0.5 * p.inertia * v * v + 0.5 * p.stiffness * q * q;
    for _ in 0..10_000 {
        (q, v) = hinge_implicit_midpoint(&p, q, v, 0.0, DT).unwrap();
    }
    let h1 = 0.5 * p.inertia * v * v + 0.5 * p.stiffness * q * q;
    assert!(
        ((h1 - h0) / h0).abs() < 1e-9,
        "implicit midpoint must conserve the oscillator energy: drift {}",
        (h1 - h0) / h0
    );
    // FALSIFIER TWIN: explicit Euler on the same fixture DIVERGES (this is
    // what the certificate exists to prevent).
    let (mut qe, mut ve) = (0.02, 0.0);
    for _ in 0..200 {
        let a = (-p.stiffness * qe - p.damping * ve) / p.inertia;
        qe += DT * ve;
        ve += DT * a;
    }
    let he = 0.5 * p.inertia * ve * ve + 0.5 * p.stiffness * qe * qe;
    assert!(
        he > 100.0 * h0,
        "the explicit twin should have exploded (he = {he:e})"
    );
    jlog(
        "implicit-stiff",
        &format!(
            "\"energy_drift\":{:e},\"explicit_blowup\":{:e}",
            (h1 - h0) / h0,
            he / h0
        ),
    );
}

#[test]
fn augmented_pole_convergence_measured_from_the_discrete_map() {
    // A damped hinge: continuous pole σ_c ± i·ω_dc with σ_c = −c/2I,
    // ω_dc = √(k/I − σ_c²). The DISCRETE one-step map (linear ⇒ built from
    // two basis steps) has eigenvalues λ_d; the recovered pole is
    // ln|λ_d|/dt and arg(λ_d)/dt. Both errors must shrink at order ~2 —
    // pole/phase convergence MEASURED, not assumed (V-05a clause).
    let p = HingeParams {
        inertia: 0.4,
        stiffness: 250.0,
        damping: 3.0,
    };
    let sigma_c = -p.damping / (2.0 * p.inertia);
    let omega_c = (p.stiffness / p.inertia - sigma_c * sigma_c).sqrt();
    let pole_err = |dt: f64| -> (f64, f64) {
        // Columns of the one-step map from basis states.
        let (q1, v1) = hinge_implicit_midpoint(&p, 1.0, 0.0, 0.0, dt).unwrap();
        let (q2, v2) = hinge_implicit_midpoint(&p, 0.0, 1.0, 0.0, dt).unwrap();
        // Eigenvalues of [[q1, q2], [v1, v2]] via trace/determinant.
        let tr = q1 + v2;
        let det = q1 * v2 - q2 * v1;
        let disc = tr * tr - 4.0 * det; // negative: complex pair
        assert!(disc < 0.0, "fixture must stay underdamped at dt = {dt}");
        let re = tr / 2.0;
        let im = (-disc).sqrt() / 2.0;
        let mag = (re * re + im * im).sqrt();
        let sigma_d = mag.ln() / dt;
        let omega_d = im.atan2(re) / dt;
        ((sigma_d - sigma_c).abs(), (omega_d - omega_c).abs())
    };
    let (s1, w1) = pole_err(1.0 / 60.0);
    let (s2, w2) = pole_err(1.0 / 120.0);
    let (s3, w3) = pole_err(1.0 / 240.0);
    let order_sigma = (s1 / s2).log2().min((s2 / s3).log2());
    let order_omega = (w1 / w2).log2().min((w2 / w3).log2());
    assert!(
        order_sigma > 1.6 || s3 < 1e-10,
        "damping-pole convergence order {order_sigma:.2} (errors {s1:e},{s2:e},{s3:e})"
    );
    assert!(
        order_omega > 1.6,
        "frequency-pole convergence order {order_omega:.2} (errors {w1:e},{w2:e},{w3:e})"
    );
    jlog(
        "pole-convergence",
        &format!("\"order_sigma\":{order_sigma:.3},\"order_omega\":{order_omega:.3}"),
    );
}

#[test]
fn certificate_caps_at_cap_and_next_float() {
    // Explicit state exactly AT the cap admits; the next float above
    // refuses naming the state and ratio.
    let tau = 1.0;
    let dt_at = MAX_EXPLICIT_RATIO * tau;
    assert!(
        certify(
            dt_at,
            &[("wake-node", tau, IntegratorClass::ExplicitResolved)]
        )
        .is_ok()
    );
    let dt_above = f64::from_bits(dt_at.to_bits() + 1);
    let refusal = certify(
        dt_above,
        &[("wake-node", tau, IntegratorClass::ExplicitResolved)],
    )
    .unwrap_err();
    assert_eq!(refusal.code, "stiffness-unresolved");
    assert!(refusal.message.contains("wake-node"));
    // Implicit/exact classes admit at ANY stiffness — that is their point.
    let stiff = certify(
        DT,
        &[
            ("hinge", 1e-4, IntegratorClass::ImplicitMidpoint),
            ("aero-memory", 1e-5, IntegratorClass::ExactTransition),
        ],
    )
    .unwrap();
    assert_eq!(stiff.entries.len(), 2);
    assert!(
        stiff.entries[0].ratio > 80.0,
        "genuinely stiff, still certified"
    );
    // Parameter gates.
    assert_eq!(
        certify(DT, &[("x", 0.0, IntegratorClass::ExactTransition)])
            .unwrap_err()
            .code,
        "timescale-invalid"
    );
    jlog(
        "certificate",
        "\"cap\":\"admitted at cap, refused at next float\"",
    );
}

#[test]
fn composition_order_verified_explicitly() {
    // The declared order is memory → (solve) → hinge → rigid. Build ONE
    // tick where the hinge torque depends on the UPDATED memory value (as
    // the real aero will): the declared step uses stage-1's output; a
    // MISCOMPOSED twin that feeds the STALE memory value must differ.
    // This proves the order is load-bearing — verified, not assumed.
    let hinge = HingeParams {
        inertia: 2.5,
        stiffness: 40.0,
        damping: 1.2,
    };
    let state = PartitionedState {
        rigid: rest(),
        hinge_q: 0.05,
        hinge_v: -0.1,
        memory: vec![MemoryState {
            tau_s: 0.2,
            value: 0.0,
        }],
    };
    let input = 1.0;
    // Declared: torque from the ADVANCED memory (stage 1 before stage 3).
    let advanced_mem = state.memory[0].advanced(input, DT).value;
    let declared = partitioned_step(
        &body(),
        &hinge,
        &state,
        &[input],
        30.0 * advanced_mem,
        Loads {
            force_n: [0.0, 0.0, 3336.0],
            moment_nm: [0.0; 3],
        },
        0.0,
        DT,
    )
    .unwrap();
    // Miscomposed: torque from the STALE memory value.
    let miscomposed = partitioned_step(
        &body(),
        &hinge,
        &state,
        &[input],
        30.0 * state.memory[0].value,
        Loads {
            force_n: [0.0, 0.0, 3336.0],
            moment_nm: [0.0; 3],
        },
        0.0,
        DT,
    )
    .unwrap();
    assert!(
        (declared.hinge_v - miscomposed.hinge_v).abs() > 1e-6,
        "composition order must be load-bearing (Δv = {})",
        (declared.hinge_v - miscomposed.hinge_v).abs()
    );
    // And the declared step is bitwise deterministic.
    let again = partitioned_step(
        &body(),
        &hinge,
        &state,
        &[input],
        30.0 * advanced_mem,
        Loads {
            force_n: [0.0, 0.0, 3336.0],
            moment_nm: [0.0; 3],
        },
        0.0,
        DT,
    )
    .unwrap();
    assert_eq!(declared.hinge_v.to_bits(), again.hinge_v.to_bits());
    assert_eq!(
        declared.memory[0].value.to_bits(),
        again.memory[0].value.to_bits()
    );
    jlog(
        "composition",
        &format!(
            "\"delta_v\":{:e}",
            (declared.hinge_v - miscomposed.hinge_v).abs()
        ),
    );
}

#[test]
fn partitioned_local_order_is_two() {
    // Richardson on the coupled fixture: hinge torque tracks the memory
    // state; the rigid heave feels the hinge via a coupling force. The
    // COUPLING CONTRACT applies: cross terms are sampled at the tick
    // MIDPOINT via a cheap predictor (tick-start sampling measured 1.22 —
    // that run is why the contract exists). Observed order must be ~2.
    let hinge = HingeParams {
        inertia: 2.5,
        stiffness: 40.0,
        damping: 1.2,
    };
    let run = |n_per_tick: u32| -> (f64, f64) {
        let dt = DT / f64::from(n_per_tick);
        let mut s = PartitionedState {
            rigid: rest(),
            hinge_q: 0.05,
            hinge_v: 0.0,
            memory: vec![MemoryState {
                tau_s: 0.2,
                value: 0.0,
            }],
        };
        let steps = 60 * n_per_tick;
        let m = body().mass_kg;
        for k in 0..steps {
            let t = f64::from(k) * dt;
            // Midpoint predictor for every cross-coupling ingredient.
            let vz0 = s.rigid.vel_mps[2];
            let az0 = (3336.0 - 400.0 * s.hinge_q - 30.0 * vz0) / m;
            let vz_mid = vz0 + 0.5 * dt * az0;
            let q_mid = s.hinge_q + 0.5 * dt * s.hinge_v;
            let input_mid = (2.0 * (t + 0.5 * dt)).sin();
            let mem_mid = s.memory[0].advanced(input_mid, 0.5 * dt).value;
            let torque = 30.0 * mem_mid - 5.0 * vz_mid;
            let loads = Loads {
                force_n: [0.0, 0.0, 3336.0 - 400.0 * q_mid - 30.0 * vz_mid],
                moment_nm: [0.0; 3],
            };
            // The contract covers memory INPUTS too: midpoint-sampled.
            s = partitioned_step(&body(), &hinge, &s, &[input_mid], torque, loads, t, dt).unwrap();
        }
        (s.hinge_q, s.rigid.vel_mps[2])
    };
    let reference = run(8);
    let err = |n: u32| -> f64 {
        let got = run(n);
        (got.0 - reference.0).abs().max((got.1 - reference.1).abs())
    };
    let (e1, e2) = (err(1), err(2));
    let order = (e1 / e2).log2();
    assert!(
        (1.6..=2.6).contains(&order),
        "partitioned local order {order:.2} not ~2 (errors {e1:e}, {e2:e})"
    );
    jlog(
        "local-order",
        &format!("\"order\":{order:.3},\"e1\":{e1:e},\"e2\":{e2:e}"),
    );
}

#[test]
fn partitioned_golden_digest() {
    // 120 declared-order ticks of the coupled fixture (measure-then-pin).
    let hinge = HingeParams {
        inertia: 2.5,
        stiffness: 40.0,
        damping: 1.2,
    };
    let mut s = PartitionedState {
        rigid: rest(),
        hinge_q: 0.05,
        hinge_v: 0.0,
        memory: vec![MemoryState {
            tau_s: 0.2,
            value: 0.0,
        }],
    };
    for k in 0..120u32 {
        let t = f64::from(k) * DT;
        let input = (2.0 * t).sin();
        let mem_next = s.memory[0].advanced(input, DT).value;
        let torque = 30.0 * mem_next;
        let loads = Loads {
            force_n: [0.0, 0.0, 3336.0],
            moment_nm: [1.0, 0.5, -0.2],
        };
        s = partitioned_step(&body(), &hinge, &s, &[input], torque, loads, t, DT).unwrap();
    }
    let mut payload = Vec::new();
    for v in [s.hinge_q, s.hinge_v, s.memory[0].value] {
        payload.extend_from_slice(&v.to_bits().to_le_bytes());
    }
    for v in s
        .rigid
        .pos_m
        .iter()
        .chain(s.rigid.vel_mps.iter())
        .chain(s.rigid.quat.iter())
    {
        payload.extend_from_slice(&v.to_bits().to_le_bytes());
    }
    let digest =
        fs_blake3::hash_domain("org.frankensim.fs-flyer.e32b-golden.v1", &payload).to_hex();
    jlog("golden", &format!("\"digest\":\"{digest}\""));
    assert_eq!(
        digest, "d5284f5bccac0695a667348578bd72b46e199d0111e258a109abec24cbde546b",
        "partitioned golden moved — determinism regression or an intentional \
         composition change requiring the golden-bump protocol"
    );
}
