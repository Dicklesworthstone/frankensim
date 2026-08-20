//! V-15 core battery (bead wf-root-guzez.5.12.2, E4.5-ii): convergence at
//! the rail state within the candidate-A cap, TWO-WAY LIVENESS
//! (anti-vacuity: the coupling must move BOTH sides), bitwise determinism
//! + warm-start, the typed nonconvergence refusal, spec-digest identity,
//! golden. Repro: cargo test -p fs-flyer --test propcoupling_battery

use fs_airscrew::{BladeStation, Rotor, bemt_solve};
use fs_flyer::propcoupling::{
    CANDIDATE_A, PropCouplingSolverSpec, PropDisk, coupled_prop_airframe_step,
};
use fs_wing::nonlinear::{InfluenceOperator, StripSpec};
use fs_wing::{SurfaceId, flat_surface};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-flyer-v15\",\"case\":\"{case}\",{payload}}}");
}

const RHO: f64 = 1.294;
const OMEGA_350: f64 = 350.0 / 60.0 * 2.0 * std::f64::consts::PI;

fn freestream() -> [f64; 3] {
    [13.86 * 0.9976, 0.0, 13.86 * 0.0699] // alpha ~ 4 deg
}

fn rotor_1903() -> Rotor {
    let deg = std::f64::consts::PI / 180.0;
    Rotor {
        radius_m: 1.2954,
        n_blades: 2,
        camber_ratio: 0.04,
        stations: vec![
            BladeStation { r_over_r: 0.30, chord_m: 0.13, beta_rad: 40.0 * deg },
            BladeStation { r_over_r: 0.45, chord_m: 0.17, beta_rad: 30.0 * deg },
            BladeStation { r_over_r: 0.60, chord_m: 0.20, beta_rad: 23.0 * deg },
            BladeStation { r_over_r: 0.75, chord_m: 0.21, beta_rad: 18.5 * deg },
            BladeStation { r_over_r: 0.88, chord_m: 0.20, beta_rad: 15.5 * deg },
            BladeStation { r_over_r: 0.96, chord_m: 0.16, beta_rad: 14.0 * deg },
        ],
    }
}

fn disks() -> [PropDisk; 2] {
    // Pusher props behind the wing (downstream = +x, the wake side),
    // straddling the centerline.
    [
        PropDisk { center_m: [3.0, -1.7, -0.9], omega_rad_s: OMEGA_350 },
        PropDisk { center_m: [3.0, 1.7, -0.9], omega_rad_s: OMEGA_350 },
    ]
}

fn wing_setup() -> (Vec<fs_wing::Panel>, Vec<StripSpec>) {
    let p = flat_surface(SurfaceId::WingLower, 12.29, 1.981, 0.0, 0.0, 8, 2).unwrap();
    let strips = (0..8)
        .map(|s| StripSpec { panel_indices: vec![s, 8 + s], chord_m: 1.981, twist_rad: 0.0 })
        .collect();
    (p, strips)
}

fn camber_closure(_s: usize, alpha: f64) -> (f64, fs_wing::nonlinear::StripRegime) {
    (2.0 * std::f64::consts::PI * (alpha + 0.1), fs_wing::nonlinear::StripRegime::Attached)
}

#[test]
fn converges_within_candidate_a_and_is_two_way_live() {
    let (p, strips) = wing_setup();
    let fs_v = freestream();
    let op = InfluenceOperator::build(&p, fs_v, RHO).unwrap();
    let r = coupled_prop_airframe_step(
        &op, &p, &strips, &camber_closure, &rotor_1903(), &disks(), fs_v, RHO, &CANDIDATE_A, None,
    )
    .unwrap();
    assert!(r.corrections <= 4, "candidate-A cap respected ({})", r.corrections);
    assert!(r.residual < CANDIDATE_A.tol, "accepted above spec tol: {}", r.residual);
    // ANTI-VACUITY, both arms:
    // (a) wing -> prop: the coupled disk inflow differs from freestream, so
    //     coupled thrust differs from the uncoupled BEMT at V_infinity.
    let uncoupled = bemt_solve(&rotor_1903(), RHO, fs_v[0], OMEGA_350).unwrap().thrust_n;
    assert!(
        (r.thrust_n[0] - uncoupled).abs() > 0.05,
        "wing->prop arm dead: {} vs uncoupled {uncoupled}",
        r.thrust_n[0]
    );
    // (b) prop -> wing: the converged slipstream changes the wing lift vs
    //     the unwashed solve.
    let unwashed = fs_wing::nonlinear::solve_nonlinear(
        &op, &p, &strips, fs_v, RHO, &camber_closure, None, None,
    )
    .unwrap()
    .total_lift_n;
    assert!(
        (r.wing_lift_n - unwashed).abs() > 1.0,
        "prop->wing arm dead: {} vs unwashed {unwashed}",
        r.wing_lift_n
    );
    // Left-right symmetry of the symmetric layout.
    assert!((r.thrust_n[0] - r.thrust_n[1]).abs() < 1e-6 * r.thrust_n[0]);
    jlog(
        "coupled",
        &format!(
            "\"thrust\":{},\"uncoupled\":{uncoupled},\"lift\":{},\"unwashed\":{unwashed},\"corr\":{}",
            r.thrust_n[0], r.wing_lift_n, r.corrections
        ),
    );
}

#[test]
fn determinism_and_warm_start() {
    let (p, strips) = wing_setup();
    let fs_v = freestream();
    let op = InfluenceOperator::build(&p, fs_v, RHO).unwrap();
    let a = coupled_prop_airframe_step(
        &op, &p, &strips, &camber_closure, &rotor_1903(), &disks(), fs_v, RHO, &CANDIDATE_A, None,
    )
    .unwrap();
    let b = coupled_prop_airframe_step(
        &op, &p, &strips, &camber_closure, &rotor_1903(), &disks(), fs_v, RHO, &CANDIDATE_A, None,
    )
    .unwrap();
    assert_eq!(a.w_slip[0].to_bits(), b.w_slip[0].to_bits(), "bitwise repeat");
    let warm = coupled_prop_airframe_step(
        &op, &p, &strips, &camber_closure, &rotor_1903(), &disks(), fs_v, RHO, &CANDIDATE_A,
        Some(a.w_slip),
    )
    .unwrap();
    assert!(warm.corrections <= a.corrections, "warm start must not cost more");
    jlog("determinism", &format!("\"cold\":{},\"warm\":{}", a.corrections, warm.corrections));
}

#[test]
fn nonconvergence_is_typed_and_spec_identity_is_sensitive() {
    let (p, strips) = wing_setup();
    let fs_v = freestream();
    let op = InfluenceOperator::build(&p, fs_v, RHO).unwrap();
    // A zero-cap spec cannot converge from a cold start: typed refusal.
    // cap 0 AND an unreachable tolerance: the free unrelaxed evaluation
    // alone must not satisfy it, so the cap gate is what fires.
    let strangled = PropCouplingSolverSpec { cap: 0, tol: 1e-15, ..CANDIDATE_A };
    let err = coupled_prop_airframe_step(
        &op, &p, &strips, &camber_closure, &rotor_1903(), &disks(), fs_v, RHO, &strangled, None,
    )
    .unwrap_err();
    assert_eq!(err.code, "PropAirframeCouplingDidNotConverge");
    assert!(err.ranked_repairs[0].contains("never switch"), "the no-one-way law travels");
    // Spec digests: candidate A is stable; any tuple change moves it.
    assert_eq!(CANDIDATE_A.digest(), CANDIDATE_A.digest());
    assert_ne!(CANDIDATE_A.digest(), strangled.digest());
    let clamped = PropCouplingSolverSpec { clamp: (0.10, 1.00), ..CANDIDATE_A };
    assert_ne!(CANDIDATE_A.digest(), clamped.digest());
    jlog("refusal+identity", &format!("\"specA\":\"{}\"", &CANDIDATE_A.digest()[..12]));
}

#[test]
fn coupling_golden_digest() {
    let (p, strips) = wing_setup();
    let fs_v = freestream();
    let op = InfluenceOperator::build(&p, fs_v, RHO).unwrap();
    let r = coupled_prop_airframe_step(
        &op, &p, &strips, &camber_closure, &rotor_1903(), &disks(), fs_v, RHO, &CANDIDATE_A, None,
    )
    .unwrap();
    let mut payload = Vec::new();
    for v in [r.w_slip[0], r.w_slip[1], r.thrust_n[0], r.thrust_n[1], r.wing_lift_n] {
        payload.extend_from_slice(&v.to_bits().to_le_bytes());
    }
    let digest =
        fs_blake3::hash_domain("org.frankensim.fs-flyer.v15-golden.v1", &payload).to_hex();
    jlog("golden", &format!("\"digest\":\"{digest}\",\"thrust\":{}", r.thrust_n[0]));
    assert_eq!(
        digest, "707c22f06100b48a34d3a8efc7624568e2afc6a4cda91a06144d53220f2d3328",
        "coupling golden moved — determinism regression or an intentional \
         scheme change requiring the golden-bump protocol"
    );
}
