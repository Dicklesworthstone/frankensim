//! V-15 core battery (bead wf-root-guzez.5.12.2, E4.5-ii): convergence at
//! the rail state within the candidate-A cap, TWO-WAY LIVENESS
//! (anti-vacuity: the coupling must move BOTH sides), bitwise determinism
//! + warm-start, the typed nonconvergence refusal, spec-digest identity,
//! golden. Repro: cargo test -p fs-flyer --test propcoupling_battery

use fs_airscrew::{BladeStation, Rotor, bemt_solve};
use fs_flyer::propcoupling::{
    CANDIDATE_A, CANDIDATE_B, IncidentInflow, PropCouplingSolverSpec, PropDisk,
    coupled_prop_airframe_step, coupled_prop_airframe_step_with_incident, wash_factor,
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
            BladeStation {
                r_over_r: 0.30,
                chord_m: 0.13,
                beta_rad: 40.0 * deg,
            },
            BladeStation {
                r_over_r: 0.45,
                chord_m: 0.17,
                beta_rad: 30.0 * deg,
            },
            BladeStation {
                r_over_r: 0.60,
                chord_m: 0.20,
                beta_rad: 23.0 * deg,
            },
            BladeStation {
                r_over_r: 0.75,
                chord_m: 0.21,
                beta_rad: 18.5 * deg,
            },
            BladeStation {
                r_over_r: 0.88,
                chord_m: 0.20,
                beta_rad: 15.5 * deg,
            },
            BladeStation {
                r_over_r: 0.96,
                chord_m: 0.16,
                beta_rad: 14.0 * deg,
            },
        ],
    }
}

fn disks() -> [PropDisk; 2] {
    // Pusher props BEHIND the wing. Frame is FRD (frame-conventions-v1):
    // +x is the nose, the wake trails toward −x — so aft is NEGATIVE x,
    // just behind the trailing edge (matching the E4.6a-ii design).
    // (The first pin of this fixture had the disks at +3.0, i.e. ahead
    // of the wing; golden re-pinned with this geometry correction and
    // the axial wash window.)
    [
        PropDisk {
            center_m: [-2.5, -1.7, -0.9],
            omega_rad_s: OMEGA_350,
        },
        PropDisk {
            center_m: [-2.5, 1.7, -0.9],
            omega_rad_s: OMEGA_350,
        },
    ]
}

fn wing_setup() -> (Vec<fs_wing::Panel>, Vec<StripSpec>) {
    let p = flat_surface(SurfaceId::WingLower, 12.29, 1.981, 0.0, 0.0, 8, 2).unwrap();
    let strips = (0..8)
        .map(|s| StripSpec {
            panel_indices: vec![s, 8 + s],
            chord_m: 1.981,
            twist_rad: 0.0,
        })
        .collect();
    (p, strips)
}

fn camber_closure(_s: usize, alpha: f64) -> (f64, fs_wing::nonlinear::StripRegime) {
    (
        2.0 * std::f64::consts::PI * (alpha + 0.1),
        fs_wing::nonlinear::StripRegime::Attached,
    )
}

#[test]
fn ordinary_wash_reconstructs_disk_induction_without_a_second_half_factor() {
    let (panels, strips) = wing_setup();
    let velocity = freestream();
    let op = InfluenceOperator::build(&panels, velocity, RHO).unwrap();
    let rotor = rotor_1903();
    let disks = disks();
    let spec = PropCouplingSolverSpec {
        tol: 1.0e-11,
        ..CANDIDATE_B
    };
    let coupled = coupled_prop_airframe_step(
        &op,
        &panels,
        &strips,
        &camber_closure,
        &rotor,
        &disks,
        velocity,
        RHO,
        &spec,
        None,
    )
    .unwrap();
    let area = core::f64::consts::PI * rotor.radius_m * rotor.radius_m;
    for k in 0..2 {
        let induced =
            fs_wing::induced_velocity_free(disks[k].center_m, &panels, &coupled.gamma, velocity);
        let disk_inflow = (velocity[0] + induced[0]).max(0.1);
        let vi = coupled.w_slip[k];
        // Independent mass-flux × far-wake velocity increment. The returned
        // momentum root is already vi, which is half the far-wake increment.
        let mass_flow = RHO * area * (disk_inflow + vi);
        let thrust = mass_flow * (2.0 * vi);
        assert!((thrust - coupled.thrust_n[k]).abs() < 1.0e-9);
    }
    let wash: Vec<f64> = strips
        .iter()
        .map(|strip| {
            let p = panels[strip.panel_indices[0]];
            let x = 0.5 * (p.a[0] + p.b[0]);
            let y = 0.5 * (p.a[1] + p.b[1]);
            let left = wash_factor(
                x - disks[0].center_m[0],
                (y - disks[0].center_m[1]).abs(),
                rotor.radius_m,
            ) * coupled.w_slip[0];
            let right = wash_factor(
                x - disks[1].center_m[0],
                (y - disks[1].center_m[1]).abs(),
                rotor.radius_m,
            ) * coupled.w_slip[1];
            left + right
        })
        .collect();
    let solve = |du: &[f64]| {
        fs_wing::nonlinear::solve_nonlinear(
            &op,
            &panels,
            &strips,
            velocity,
            RHO,
            &camber_closure,
            Some(du),
            None,
        )
        .unwrap()
    };
    let physical = solve(&wash);
    assert!((physical.total_lift_n - coupled.wing_lift_n).abs() < 1.0e-7);
    for (expected, actual) in physical.gamma.iter().zip(&coupled.gamma) {
        assert!((expected - actual).abs() < 1.0e-8);
    }
    let half_wash: Vec<_> = wash.iter().map(|du| 0.5 * du).collect();
    let double_halved = solve(&half_wash);
    assert!(
        (double_halved.total_lift_n - coupled.wing_lift_n).abs() > 0.2,
        "the native fixture must detect an extra factor 1/2 in momentum transfer"
    );
}

#[test]
fn distributed_incident_field_drives_separate_rotors_and_consistent_wash() {
    let (panels, strips) = wing_setup();
    let velocity = freestream();
    let op = InfluenceOperator::build(&panels, velocity, RHO).unwrap();
    let rotor = rotor_1903();
    let disks = disks();
    let local = vec![velocity; strips.len()];
    let disk_flow = [[10.0, 0.0, 0.0], [16.0, 0.0, 0.0]];
    let result = coupled_prop_airframe_step_with_incident(
        &op,
        &panels,
        &strips,
        &camber_closure,
        &rotor,
        &disks,
        velocity,
        RHO,
        &CANDIDATE_B,
        None,
        IncidentInflow {
            strip_velocity_mps: &local,
            disk_velocity_mps: disk_flow,
        },
    )
    .unwrap();
    assert!(result.thrust_n[0] > result.thrust_n[1]);
    assert_ne!(result.torque_nm[0].to_bits(), result.torque_nm[1].to_bits());
    let area = core::f64::consts::PI * rotor.radius_m * rotor.radius_m;
    for k in 0..2 {
        let induced =
            fs_wing::induced_velocity_free(disks[k].center_m, &panels, &result.gamma, velocity);
        let axial = disk_flow[k][0] + induced[0];
        let vi = result.w_slip[k];
        let momentum_thrust = 2.0 * RHO * area * vi * (axial + vi);
        assert!((momentum_thrust - result.thrust_n[k]).abs() < 1.0e-9);
    }
    let accepted = result.strip_velocity_mps.as_ref().unwrap();
    let mut washed = 0;
    let mut disk_checked = [false; 2];
    for (s, strip) in strips.iter().enumerate() {
        let p = panels[strip.panel_indices[0]];
        let x = 0.5 * (p.a[0] + p.b[0]);
        let y = 0.5 * (p.a[1] + p.b[1]);
        let factors = disks.map(|disk| {
            wash_factor(
                x - disk.center_m[0],
                (y - disk.center_m[1]).abs(),
                rotor.radius_m,
            )
        });
        let expected = factors[0] * result.w_slip[0] + factors[1] * result.w_slip[1];
        let actual = accepted[s][0] - velocity[0];
        // Gamma and force velocity belong to evaluated x; w_slip is G(x).
        // Bound their difference by each disk's declared relative residual.
        let bound = CANDIDATE_B.tol / (1.0 - CANDIDATE_B.tol)
            * (factors[0].abs() * (result.w_slip[0].abs() + 1.0e-6)
                + factors[1].abs() * (result.w_slip[1].abs() + 1.0e-6));
        assert!((actual - expected).abs() <= bound + 1.0e-12);
        for k in 0..2 {
            if factors[k] > 0.0 && factors[1 - k] == 0.0 {
                let evaluated = actual / factors[k];
                let relative = (result.w_slip[k] - evaluated).abs() / (evaluated.abs() + 1.0e-6);
                assert!(relative < CANDIDATE_B.tol + 1.0e-11);
                disk_checked[k] = true;
            }
        }
        if expected > 0.0 {
            washed += 1;
        }
    }
    assert!(
        washed > 0,
        "fixture must include strips inside real disk wash"
    );
    assert!(disk_checked.into_iter().all(|checked| checked));
}

#[test]
fn malformed_strip_map_refuses_before_sampling_wash() {
    let (panels, strips) = wing_setup();
    let velocity = freestream();
    let op = InfluenceOperator::build(&panels, velocity, RHO).unwrap();
    let mut empty_strip = strips.clone();
    empty_strip[0].panel_indices.clear();
    let mut invalid_index = strips.clone();
    invalid_index[0].panel_indices[0] = panels.len();
    for malformed in [Vec::new(), empty_strip, invalid_index] {
        let local = vec![velocity; malformed.len()];
        let legacy = coupled_prop_airframe_step(
            &op,
            &panels,
            &malformed,
            &camber_closure,
            &rotor_1903(),
            &disks(),
            velocity,
            RHO,
            &CANDIDATE_B,
            None,
        )
        .unwrap_err();
        let distributed = coupled_prop_airframe_step_with_incident(
            &op,
            &panels,
            &malformed,
            &camber_closure,
            &rotor_1903(),
            &disks(),
            velocity,
            RHO,
            &CANDIDATE_B,
            None,
            IncidentInflow {
                strip_velocity_mps: &local,
                disk_velocity_mps: [velocity; 2],
            },
        )
        .unwrap_err();
        assert_eq!(legacy.code, "strips-invalid");
        assert_eq!(distributed.code, "strips-invalid");
    }
}

#[test]
fn distributed_disk_reversal_refuses_without_uniform_fallback() {
    let (panels, strips) = wing_setup();
    let velocity = freestream();
    let op = InfluenceOperator::build(&panels, velocity, RHO).unwrap();
    let local = vec![velocity; strips.len()];
    for axial in [-100.0, f64::NAN] {
        let error = coupled_prop_airframe_step_with_incident(
            &op,
            &panels,
            &strips,
            &camber_closure,
            &rotor_1903(),
            &disks(),
            velocity,
            RHO,
            &CANDIDATE_B,
            None,
            IncidentInflow {
                strip_velocity_mps: &local,
                disk_velocity_mps: [[axial, 0.0, 0.0]; 2],
            },
        )
        .unwrap_err();
        assert_eq!(error.code, "disk-inflow-invalid");
    }
}

#[test]
fn converges_within_candidate_a_and_is_two_way_live() {
    let (p, strips) = wing_setup();
    let fs_v = freestream();
    let op = InfluenceOperator::build(&p, fs_v, RHO).unwrap();
    let r = coupled_prop_airframe_step(
        &op,
        &p,
        &strips,
        &camber_closure,
        &rotor_1903(),
        &disks(),
        fs_v,
        RHO,
        &CANDIDATE_A,
        None,
    )
    .unwrap();
    assert!(
        r.corrections <= 4,
        "candidate-A cap respected ({})",
        r.corrections
    );
    assert!(
        r.residual < CANDIDATE_A.tol,
        "accepted above spec tol: {}",
        r.residual
    );
    // ANTI-VACUITY, both arms:
    // (a) wing -> prop: the coupled disk inflow differs from freestream, so
    //     coupled thrust differs from the uncoupled BEMT at V_infinity.
    let uncoupled = bemt_solve(&rotor_1903(), RHO, fs_v[0], OMEGA_350)
        .unwrap()
        .thrust_n;
    assert!(
        (r.thrust_n[0] - uncoupled).abs() > 0.05,
        "wing->prop arm dead: {} vs uncoupled {uncoupled}",
        r.thrust_n[0]
    );
    // (b) prop -> wing: the converged slipstream changes the wing lift vs
    //     the unwashed solve.
    let unwashed = fs_wing::nonlinear::solve_nonlinear(
        &op,
        &p,
        &strips,
        fs_v,
        RHO,
        &camber_closure,
        None,
        None,
    )
    .unwrap()
    .total_lift_n;
    // Re-measured 2026-08-21 under the 5.7.1 upstream-wash decay: the
    // wing (ahead of the disks) now feels the DECAYED factor, so the
    // live-arm delta is 0.587 N at this fixture (was ~a few N under
    // the flat map). 0.2 keeps ~3x headroom while still proving the
    // arm carries real coupling, not float noise.
    assert!(
        (r.wing_lift_n - unwashed).abs() > 0.2,
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
        &op,
        &p,
        &strips,
        &camber_closure,
        &rotor_1903(),
        &disks(),
        fs_v,
        RHO,
        &CANDIDATE_A,
        None,
    )
    .unwrap();
    let b = coupled_prop_airframe_step(
        &op,
        &p,
        &strips,
        &camber_closure,
        &rotor_1903(),
        &disks(),
        fs_v,
        RHO,
        &CANDIDATE_A,
        None,
    )
    .unwrap();
    assert_eq!(
        a.w_slip[0].to_bits(),
        b.w_slip[0].to_bits(),
        "bitwise repeat"
    );
    let warm = coupled_prop_airframe_step(
        &op,
        &p,
        &strips,
        &camber_closure,
        &rotor_1903(),
        &disks(),
        fs_v,
        RHO,
        &CANDIDATE_A,
        Some(a.w_slip),
    )
    .unwrap();
    assert!(
        warm.corrections <= a.corrections,
        "warm start must not cost more"
    );
    jlog(
        "determinism",
        &format!("\"cold\":{},\"warm\":{}", a.corrections, warm.corrections),
    );
}

#[test]
fn nonconvergence_is_typed_and_spec_identity_is_sensitive() {
    let (p, strips) = wing_setup();
    let fs_v = freestream();
    let op = InfluenceOperator::build(&p, fs_v, RHO).unwrap();
    // A zero-cap spec cannot converge from a cold start: typed refusal.
    // cap 0 AND an unreachable tolerance: the free unrelaxed evaluation
    // alone must not satisfy it, so the cap gate is what fires.
    let strangled = PropCouplingSolverSpec {
        cap: 0,
        tol: 1e-15,
        ..CANDIDATE_A
    };
    let err = coupled_prop_airframe_step(
        &op,
        &p,
        &strips,
        &camber_closure,
        &rotor_1903(),
        &disks(),
        fs_v,
        RHO,
        &strangled,
        None,
    )
    .unwrap_err();
    assert_eq!(err.code, "PropAirframeCouplingDidNotConverge");
    assert!(
        err.ranked_repairs[0].contains("never switch"),
        "the no-one-way law travels"
    );
    // Spec digests: candidate A is stable; any tuple change moves it.
    assert_eq!(CANDIDATE_A.digest(), CANDIDATE_A.digest());
    assert_ne!(CANDIDATE_A.digest(), strangled.digest());
    let clamped = PropCouplingSolverSpec {
        clamp: (0.10, 1.00),
        ..CANDIDATE_A
    };
    assert_ne!(CANDIDATE_A.digest(), clamped.digest());
    jlog(
        "refusal+identity",
        &format!("\"specA\":\"{}\"", &CANDIDATE_A.digest()[..12]),
    );
}

#[test]
fn grown_native_correction_rolls_back_and_keeps_the_existing_budget() {
    // Full-strength disk wash at an admitted off-trim native state. Its
    // separated-strip map still chatters between branches; a correct retry
    // must not relabel that unresolved physical solve as converged.
    let mut design = fs_flyer::aircraft::wright_openloop_v1();
    for (cap, expected_residual) in [
        (5, 0.086_412_799_529_274_37),
        (12, 0.086_011_921_941_914_73),
    ] {
        design.coupling = PropCouplingSolverSpec { cap, ..CANDIDATE_B };
        let error = design
            .force_buildup_warm(
                10.0,
                -0.54,
                0.116_761_945_162_323_48,
                49.437_682_658_920_99,
                -2.5,
                RHO,
                Some([2.0, 2.0]),
            )
            .unwrap_err();
        assert_eq!(error.code, "PropAirframeCouplingDidNotConverge");
        let prefix = format!("cap {cap} exhausted: joint residual ");
        let residual: f64 = error
            .message
            .strip_prefix(&prefix)
            .expect("a recovering correction must not consume a false growth strike")
            .parse()
            .unwrap();
        // Correction 4 grows to residual 0.20955. Retrying from its
        // predecessor gives 0.08641 at cap 5; the old algorithm instead
        // advanced from the rejected point, leaving 0.18505 and then
        // reporting a false second growth before the cap was exhausted.
        assert!((residual - expected_residual).abs() < 1.0e-8);
        assert!(residual > CANDIDATE_B.tol, "unresolved loads stay refused");
        assert!(error.ranked_repairs[0].contains("never switch"));
        jlog(
            "rollback-budget",
            &format!(
                "\"cap\":{cap},\"residual\":{residual},\"tol\":{}",
                CANDIDATE_B.tol
            ),
        );
    }
}

#[test]
fn correct_disk_wash_resolves_the_previous_member3_stopping_state() {
    // Exact former member-3 midpoint at tick 1443. Keep the native plant,
    // rotor, controls, warm state, candidate B and its tolerance unchanged.
    let mut design = fs_flyer::aircraft::wright_openloop_v1();
    design.coupling = CANDIDATE_B;
    let result = design
        .force_buildup_warm(
            12.109_920_004_298_539,
            -0.527_706_488_358_26,
            0.116_761_945_162_323_48,
            49.437_682_658_920_99,
            -2.9,
            RHO,
            Some([1.660_416_864_469_997_5, 1.660_416_865_101_871_6]),
        )
        .unwrap();
    assert!(result.coupled.residual < CANDIDATE_B.tol);
    assert_eq!(result.coupled.corrections, 0);
    assert!((result.thrust_n[0] - result.thrust_n[1]).abs() < 1.0e-8);
    assert!(result.force_n.iter().all(|f| f.is_finite()));
    jlog(
        "previous-stopping-state",
        &format!(
            "\"residual\":{},\"thrust_left\":{},\"thrust_right\":{}",
            result.coupled.residual, result.thrust_n[0], result.thrust_n[1]
        ),
    );
}

#[test]
fn coupling_golden_digest() {
    let (p, strips) = wing_setup();
    let fs_v = freestream();
    let op = InfluenceOperator::build(&p, fs_v, RHO).unwrap();
    let r = coupled_prop_airframe_step(
        &op,
        &p,
        &strips,
        &camber_closure,
        &rotor_1903(),
        &disks(),
        fs_v,
        RHO,
        &CANDIDATE_A,
        None,
    )
    .unwrap();
    let mut payload = Vec::new();
    for v in [
        r.w_slip[0],
        r.w_slip[1],
        r.thrust_n[0],
        r.thrust_n[1],
        r.wing_lift_n,
    ] {
        payload.extend_from_slice(&v.to_bits().to_le_bytes());
    }
    let digest = fs_blake3::hash_domain("org.frankensim.fs-flyer.v15-golden.v1", &payload).to_hex();
    jlog(
        "golden",
        &format!("\"digest\":\"{digest}\",\"thrust\":{}", r.thrust_n[0]),
    );
    // Golden bump 2026-10-28: both loading paths now apply the disk-induced
    // velocity once, using the same rationalized momentum root. The original
    // 350 RPM geometry, signed annular solve and all reciprocal-coupling,
    // warm-start and liveness gates are unchanged.
    assert_eq!(
        digest, "281389a06cab98db04c3f27a581a1b41ff10f2581dcebc9fc81d66a629c60210",
        "coupling golden moved — determinism regression or an intentional \
         scheme change requiring the golden-bump protocol"
    );
}

#[test]
fn wash_factor_law_upstream_decay_and_downstream_growth() {
    // Bead wf-root-guzez.5.7.1: the per-pair factor LAW, pinned
    // against the actuator-disk momentum forms (per-point oracles).
    use fs_flyer::propcoupling::wash_factor;
    let r: f64 = 1.2954;
    // At the disk plane: exactly 1 (the half-slipstream disk value).
    assert!((wash_factor(0.0, 0.0, r) - 1.0).abs() < 1e-15);
    // Upstream 1.9R (the wing's station): 1 - 1.9/sqrt(1.9^2+1).
    let dx: f64 = 1.9 * r;
    let expect = 1.0 - dx / (dx * dx + r * r).sqrt();
    assert!((wash_factor(dx, 0.0, r) - expect).abs() < 1e-15);
    assert!(
        expect < 0.12,
        "the wing's upstream factor is SMALL: {expect}"
    );
    // Downstream mirror grows toward the full slipstream (factor 2).
    let down = wash_factor(-dx, 0.0, r);
    assert!((down - (1.0 + dx / (dx * dx + r * r).sqrt())).abs() < 1e-15);
    assert!(down > 1.88);
    // Monotone decay upstream (per-point over the window).
    let mut prev = wash_factor(0.0, 0.0, r);
    for k in 1..20 {
        let f = wash_factor(k as f64 * 0.1 * r, 0.0, r);
        assert!(f < prev, "monotone at {k}");
        prev = f;
    }
    // Gates: AT the axial window edge and beyond -> 0; radial gate.
    assert_eq!(wash_factor(2.0 * r, 0.0, r), 0.0);
    assert_eq!(wash_factor(-2.0 * r, 0.0, r), 0.0);
    assert_eq!(
        wash_factor(3.65 * r, 0.0, r),
        0.0,
        "the canard stays unwashed"
    );
    assert_eq!(wash_factor(0.5, r, r), 0.0, "radial gate AT R");
    assert!(
        wash_factor(0.5, 0.99 * r, r) > 0.0,
        "inside the radial gate"
    );
    jlog(
        "wash-factor-law",
        &format!("\"upstream_1p9R\":{expect},\"downstream_1p9R\":{down}"),
    );
}
