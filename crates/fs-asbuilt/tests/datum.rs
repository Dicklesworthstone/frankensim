//! Battery for datum-priority (3-2-1) registration. Naming: `g0_` analytic
//! oracles, `g3_` structural hierarchy invariances, `g4_` cancellation,
//! `g5_` determinism. The e2e case drives a scan-like noisy dataset through
//! the global, datum, and calibrated paths and logs the forensic table the
//! bead demands.
//!
//! The datum A set deliberately carries SIX targets: with only four corner
//! targets a plane fit has hat-diagonal 3/4 and spreads a single-corner
//! deviation into equal-magnitude residuals at every corner, which would
//! falsify the "datum path concentrates the deviation" property this battery
//! locks. Six targets keep the deviated corner's leverage low enough that
//! concentration genuinely holds.

#![allow(clippy::needless_range_loop)] // Fixed 3x3 indices mirror the source modules.
#![allow(clippy::float_cmp)] // G5 bitwise-replay assertions compare exact bits on purpose.

use fs_asbuilt::datum::{
    DatumDegeneracy, DatumError, DatumLabel, DatumSystem, FitSide, register3_datum,
};
use fs_asbuilt::rigid3::{
    Covariance3, CrossFiducialModel3, Fiducial3, MetrologyModel3, Point3,
    estimate_calibrated_rigid3, register3, register3_similarity,
};
use fs_asbuilt::uncertainty::HuberPolicy;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use std::fmt::Write as _;

fn with_cx<R>(cancelled: bool, mode: ExecMode, budget: Budget, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    if cancelled {
        gate.request();
    }
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let clock = fs_exec::VirtualClock::new();
    let result = pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x0A5B_0118,
                kernel_id: 4,
                tile: 0,
                iteration: 0,
            },
            budget,
            mode,
        )
        .with_time_source(&clock);
        f(&cx)
    });
    let stats = pool.stats();
    assert!(
        stats.quiescent(),
        "Cx arena must be quiescent after scope: {}",
        stats.to_json()
    );
    result
}

fn with_default_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    with_cx(false, ExecMode::Deterministic, Budget::INFINITE, f)
}

fn p3(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z).expect("finite test point")
}

fn yaw_matrix(angle: f64) -> [[f64; 3]; 3] {
    let (s, c) = angle.sin_cos();
    [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]]
}

fn apply(rotation: &[[f64; 3]; 3], translation: [f64; 3], p: Point3) -> Point3 {
    p3(
        rotation[0][0] * p.x() + rotation[0][1] * p.y() + rotation[0][2] * p.z() + translation[0],
        rotation[1][0] * p.x() + rotation[1][1] * p.y() + rotation[1][2] * p.z() + translation[1],
        rotation[2][0] * p.x() + rotation[2][1] * p.y() + rotation[2][2] * p.z() + translation[2],
    )
}

/// Drawing-style block fixture. Datum A: six bottom-face targets — four
/// corners plus two edge midpoints (0..=5). Datum B: two targets on the y=0
/// side face at equal height (6, 7). Datum C: one corner target (8). Three
/// free fiducials ride along (9..=11).
fn block_design() -> Vec<Point3> {
    vec![
        p3(0.0, 0.0, 0.0),
        p3(4.0, 0.0, 0.0),
        p3(4.0, 3.0, 0.0),
        p3(0.0, 3.0, 0.0),
        p3(2.0, 0.0, 0.0),
        p3(2.0, 3.0, 0.0),
        p3(0.5, 0.0, 0.5),
        p3(3.5, 0.0, 0.5),
        p3(0.0, 0.5, 0.8),
        p3(2.0, 1.5, 1.0),
        p3(1.0, 2.0, 0.3),
        p3(3.0, 2.5, 0.9),
    ]
}

fn block_system() -> DatumSystem {
    DatumSystem::new(vec![0, 1, 2, 3, 4, 5], vec![6, 7], 8).expect("valid block datum system")
}

const TRUTH_YAW: f64 = core::f64::consts::FRAC_PI_6; // 30 degrees
const TRUTH_TRANSLATION: [f64; 3] = [5.0, -2.0, 0.75];

fn exact_fiducials() -> Vec<Fiducial3> {
    let rotation = yaw_matrix(TRUTH_YAW);
    block_design()
        .iter()
        .map(|p| Fiducial3::new(*p, apply(&rotation, TRUTH_TRANSLATION, *p)))
        .collect()
}

/// Deterministic test-local RNG (xorshift64* + Box-Muller); production
/// randomness lives in fs-rand.
struct TestRng(u64);

impl TestRng {
    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    fn uniform(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / 9_007_199_254_740_992.0 + f64::EPSILON
    }

    fn normal(&mut self) -> f64 {
        let u1 = self.uniform();
        let u2 = self.uniform();
        (-2.0 * u1.ln()).sqrt() * (core::f64::consts::TAU * u2).cos()
    }
}

#[test]
fn g0_hand_worked_block_recovers_the_seeded_pose_exactly() {
    let fiducials = exact_fiducials();
    let system = block_system();
    let datum = with_default_cx(|cx| register3_datum(&fiducials, &system, cx))
        .expect("consistent datums recover");
    let expected = yaw_matrix(TRUTH_YAW);
    for row in 0..3 {
        for column in 0..3 {
            assert!(
                (datum.rotation()[row][column] - expected[row][column]).abs() < 1e-10,
                "rotation[{row}][{column}] drifted"
            );
        }
    }
    for axis in 0..3 {
        assert!((datum.translation()[axis] - TRUTH_TRANSLATION[axis]).abs() < 1e-9);
    }
    for residual in datum.a_out_of_plane() {
        assert!(residual.abs() < 1e-9);
    }
    for residual in datum.b_off_line() {
        assert!(residual.abs() < 1e-9);
    }
    assert!(datum.c_along_line().abs() < 1e-9);
    for norm in datum.residual_norms() {
        assert!(*norm < 1e-9);
    }
    // Noise-free consistent data: datum and global agree.
    assert!(datum.comparison().rotation_delta_rad() < 1e-9);
    for axis in 0..3 {
        assert!(datum.comparison().translation_delta()[axis].abs() < 1e-8);
    }
}

#[test]
fn g3_datum_b_out_of_plane_information_is_structurally_discarded() {
    let baseline = exact_fiducials();
    let system = block_system();
    let reference =
        with_default_cx(|cx| register3_datum(&baseline, &system, cx)).expect("reference");

    // Push BOTH B targets along the measured A normal (z, for a yaw): the
    // hierarchy discards B's out-of-plane component, so the datum pose must
    // not move at all, while the global best fit must.
    let mut perturbed = baseline.clone();
    for index in [6usize, 7] {
        let measured = perturbed[index].measured();
        perturbed[index] = Fiducial3::new(
            perturbed[index].design(),
            p3(measured.x(), measured.y(), measured.z() + 0.2),
        );
    }
    let shifted =
        with_default_cx(|cx| register3_datum(&perturbed, &system, cx)).expect("perturbed");
    for row in 0..3 {
        for column in 0..3 {
            assert!(
                (shifted.rotation()[row][column] - reference.rotation()[row][column]).abs() < 1e-12,
                "datum rotation moved under a B out-of-plane perturbation"
            );
        }
    }
    for axis in 0..3 {
        assert!(
            (shifted.translation()[axis] - reference.translation()[axis]).abs() < 1e-12,
            "datum translation moved under a B out-of-plane perturbation"
        );
    }
    // The same perturbation must move the global fit; the delta diagnostic
    // now shows the disagreement.
    let global_shift: f64 = (0..3)
        .map(|axis| {
            (shifted.comparison().global().translation()[axis]
                - reference.comparison().global().translation()[axis])
                .abs()
        })
        .sum();
    assert!(
        global_shift > 1e-3,
        "global fit should absorb the B perturbation, moved {global_shift:.3e}"
    );
    // The discarded out-of-plane deviation stays visible in the full
    // residual norms, though not in the off-line component B constrains.
    for index in [6usize, 7] {
        assert!(shifted.residual_norms()[index] > 0.19);
    }
}

#[test]
fn g3_datum_c_transverse_information_is_structurally_discarded() {
    let baseline = exact_fiducials();
    let system = block_system();
    let reference =
        with_default_cx(|cx| register3_datum(&baseline, &system, cx)).expect("reference");

    // Move C along the plane normal and the in-plane transverse direction:
    // both components are constrained by A and B, not C, so the datum pose
    // must not move.
    let mut perturbed = baseline.clone();
    let frame = reference.frame();
    let offset = [
        0.3 * frame.plane_normal()[0] + 0.2 * frame.in_plane_transverse()[0],
        0.3 * frame.plane_normal()[1] + 0.2 * frame.in_plane_transverse()[1],
        0.3 * frame.plane_normal()[2] + 0.2 * frame.in_plane_transverse()[2],
    ];
    let measured = perturbed[8].measured();
    perturbed[8] = Fiducial3::new(
        perturbed[8].design(),
        p3(
            measured.x() + offset[0],
            measured.y() + offset[1],
            measured.z() + offset[2],
        ),
    );
    let shifted =
        with_default_cx(|cx| register3_datum(&perturbed, &system, cx)).expect("perturbed");
    for axis in 0..3 {
        assert!(
            (shifted.translation()[axis] - reference.translation()[axis]).abs() < 1e-12,
            "datum translation moved under a C transverse perturbation"
        );
    }
    assert!(shifted.c_along_line().abs() < 1e-9);
    assert!(shifted.residual_norms()[8] > 0.3);
}

#[test]
fn seeded_local_deviation_is_exposed_by_the_datum_path() {
    let mut fiducials = exact_fiducials();
    let system = block_system();
    // Bend one A corner out of plane: a form deviation on the primary datum.
    let deviated = 2usize;
    let deviation = 0.2f64;
    let measured = fiducials[deviated].measured();
    fiducials[deviated] = Fiducial3::new(
        fiducials[deviated].design(),
        p3(measured.x(), measured.y(), measured.z() + deviation),
    );

    let datum = with_default_cx(|cx| register3_datum(&fiducials, &system, cx)).expect("fits");

    // The deviated corner dominates the A out-of-plane residuals. With the
    // six-target A set its plane-fit leverage is ~0.58, so roughly 0.42 of
    // the seeded deviation must surface at the deviated target itself.
    let deviated_a_ordinal = 2usize; // position of fiducial 2 within the A set
    let exposed = datum.a_out_of_plane()[deviated_a_ordinal];
    assert!(
        exposed > 0.3 * deviation,
        "datum path must expose the seeded {deviation} deviation, saw {exposed:.4}"
    );
    for (ordinal, residual) in datum.a_out_of_plane().iter().enumerate() {
        if ordinal != deviated_a_ordinal {
            assert!(
                residual.abs() < exposed,
                "A residual {ordinal} ({residual:.4}) must stay below the deviated target"
            );
        }
    }
    // The two frames must measurably disagree at the deviated target: the
    // datum plane fit (six A targets, leverage ~0.58 at this corner) absorbs
    // more of the deviation locally than the twelve-point global fit, so the
    // delta's SIGN is geometry-dependent — the published diagnostic is the
    // disagreement itself, not a one-sided inequality.
    let disagreement = datum.comparison().residual_norm_deltas()[deviated];
    assert!(
        disagreement.abs() > 0.05 * deviation,
        "the frames must disagree at the deviated target, saw {disagreement:.4}"
    );
    assert!(
        datum.comparison().rotation_delta_rad() > 1e-4,
        "the two poses must visibly disagree"
    );
}

#[test]
fn datum_system_and_registration_refusals_are_typed() {
    assert!(matches!(
        DatumSystem::new(vec![0, 1], vec![6, 7], 8),
        Err(DatumError::System { .. })
    ));
    assert!(matches!(
        DatumSystem::new(vec![0, 1, 2], vec![6], 8),
        Err(DatumError::System { .. })
    ));
    assert!(matches!(
        DatumSystem::new(vec![0, 1, 2], vec![2, 6], 8),
        Err(DatumError::DuplicateIndex { index: 2 })
    ));
    assert!(matches!(
        DatumSystem::new(vec![0, 1, 2], vec![3, 6], 6),
        Err(DatumError::DuplicateIndex { index: 6 })
    ));

    let fiducials = exact_fiducials();
    let out_of_range = DatumSystem::new(vec![0, 1, 2], vec![6, 7], 99).expect("shape valid");
    assert!(matches!(
        with_default_cx(|cx| register3_datum(&fiducials, &out_of_range, cx)),
        Err(DatumError::IndexOutOfRange { index: 99, .. })
    ));

    // (0,0,0), (4,0,0), (2,0,0) are collinear: no plane is defined.
    let collinear_a = DatumSystem::new(vec![0, 1, 4], vec![6, 7], 8).expect("shape valid");
    assert!(matches!(
        with_default_cx(|cx| register3_datum(&fiducials, &collinear_a, cx)),
        Err(DatumError::DegenerateDatum {
            datum: DatumLabel::A,
            side: FitSide::Design,
            diagnosis: DatumDegeneracy::Collinear
        })
    ));

    // B chord parallel to the A normal constrains no in-plane rotation.
    let mut stacked = exact_fiducials();
    for (index, design) in [(6usize, p3(0.5, 0.0, 0.2)), (7usize, p3(0.5, 0.0, 0.9))] {
        stacked[index] = Fiducial3::new(
            design,
            apply(&yaw_matrix(TRUTH_YAW), TRUTH_TRANSLATION, design),
        );
    }
    let system = block_system();
    assert!(matches!(
        with_default_cx(|cx| register3_datum(&stacked, &system, cx)),
        Err(DatumError::BParallelToA {
            side: FitSide::Design
        })
    ));

    // A measured plane exactly perpendicular to the design plane cannot be
    // sign-oriented; the within-90-degrees precondition refuses.
    let quarter = [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]]; // 90 deg about x
    let perpendicular: Vec<Fiducial3> = block_design()
        .iter()
        .map(|p| Fiducial3::new(*p, apply(&quarter, [0.0, 0.0, 0.0], *p)))
        .collect();
    assert!(matches!(
        with_default_cx(|cx| register3_datum(&perpendicular, &system, cx)),
        Err(DatumError::OrientationAmbiguous {
            datum: DatumLabel::A
        })
    ));
}

#[test]
#[allow(clippy::too_many_lines)]
fn e2e_scan_like_dataset_logs_the_full_diagnostic_table() {
    let rotation = yaw_matrix(TRUTH_YAW);
    let sigma = 5e-4f64;
    let mut rng = TestRng(0xE2E_0001);
    let mut fiducials: Vec<Fiducial3> = block_design()
        .iter()
        .map(|p| {
            let q = apply(&rotation, TRUTH_TRANSLATION, *p);
            Fiducial3::new(
                *p,
                p3(
                    q.x() + sigma * rng.normal(),
                    q.y() + sigma * rng.normal(),
                    q.z() + sigma * rng.normal(),
                ),
            )
        })
        .collect();
    // Seed a genuine local deviation on an A corner, well above noise.
    let deviated = 2usize;
    let deviation = 0.05f64;
    let measured = fiducials[deviated].measured();
    fiducials[deviated] = Fiducial3::new(
        fiducials[deviated].design(),
        p3(measured.x(), measured.y(), measured.z() + deviation),
    );

    let system = block_system();
    let covariance = Covariance3::new(sigma * sigma, 0.0, 0.0, sigma * sigma, 0.0, sigma * sigma)
        .expect("iso covariance");
    let model = MetrologyModel3::new(
        vec![covariance; fiducials.len()],
        CrossFiducialModel3::Independent,
        HuberPolicy::new(3.0, 6).expect("huber"),
        "cmm-calibration-2026-07/block-fixture@rev3",
    )
    .expect("model");

    let (global, datum, calibrated, similarity) = with_default_cx(|cx| {
        (
            register3(&fiducials, cx).expect("global fit"),
            register3_datum(&fiducials, &system, cx).expect("datum fit"),
            estimate_calibrated_rigid3(&fiducials, &model, cx).expect("calibrated fit"),
            register3_similarity(&fiducials, 0.01, cx).expect("similarity fit"),
        )
    });

    // The datum path concentrates the deviation on the seeded target.
    let a_ordinal = 2usize;
    let exposed = datum.a_out_of_plane()[a_ordinal];
    assert!(
        exposed > 0.3 * deviation,
        "seeded {deviation} deviation exposed, saw {exposed:.4}"
    );
    for (ordinal, residual) in datum.a_out_of_plane().iter().enumerate() {
        if ordinal != a_ordinal {
            assert!(residual.abs() < exposed);
        }
    }
    let disagreement = datum.comparison().residual_norm_deltas()[deviated];
    assert!(
        disagreement.abs() > 0.05 * deviation,
        "the frames must disagree at the deviated target, saw {disagreement:.4}"
    );
    // No unit error was seeded: the similarity scale sits at one.
    assert!(similarity.scale().suspicion().is_none());
    // Calibrated covariance is published and positive on the diagonal.
    for axis in 0..6 {
        assert!(calibrated.covariance()[axis][axis] > 0.0);
    }

    // Forensic log: poses, covariance diagonal, per-datum residuals, and the
    // datum-versus-global delta, as the acceptance lane demands.
    let mut log = String::new();
    let _ = writeln!(
        log,
        "global: angle={:.6} rad translation={:?} rms={:.3e}",
        global.rotation_angle_rad(),
        global.translation(),
        global.residual_rms()
    );
    let _ = writeln!(
        log,
        "datum: angle={:.6} rad translation={:?}",
        datum.rotation_angle_rad(),
        datum.translation()
    );
    let _ = writeln!(
        log,
        "datum residuals: A={:?} B={:?} C={:.3e}",
        datum.a_out_of_plane(),
        datum.b_off_line(),
        datum.c_along_line()
    );
    let _ = writeln!(
        log,
        "datum-vs-global: rotation_delta={:.3e} rad translation_delta={:?} residual_deltas={:?}",
        datum.comparison().rotation_delta_rad(),
        datum.comparison().translation_delta(),
        datum.comparison().residual_norm_deltas()
    );
    let _ = writeln!(
        log,
        "calibrated: dof={} identity={} covariance_diag={:?}",
        calibrated.degrees_of_freedom(),
        calibrated.model_identity(),
        (0..6)
            .map(|axis| calibrated.covariance()[axis][axis])
            .collect::<Vec<_>>()
    );
    let _ = writeln!(
        log,
        "similarity: scale={:.9} se={:.3e} suspicion={:?}",
        similarity.scale().estimate(),
        similarity.scale().standard_error(),
        similarity.scale().suspicion()
    );
    for marker in [
        "global: angle=",
        "datum: angle=",
        "datum residuals: A=",
        "datum-vs-global: rotation_delta=",
        "calibrated: dof=30",
        "similarity: scale=",
    ] {
        assert!(log.contains(marker), "forensic log lost {marker:?}\n{log}");
    }
    println!("{log}");
}

#[test]
fn g5_datum_replay_is_bitwise() {
    let fiducials = exact_fiducials();
    let system = block_system();
    let (first, second) = with_default_cx(|cx| {
        (
            register3_datum(&fiducials, &system, cx).expect("first"),
            register3_datum(&fiducials, &system, cx).expect("second"),
        )
    });
    assert_eq!(first.rotation(), second.rotation());
    assert_eq!(first.translation(), second.translation());
    assert_eq!(first.residual_norms(), second.residual_norms());
    assert_eq!(
        first.comparison().rotation_delta_rad(),
        second.comparison().rotation_delta_rad()
    );
}

#[test]
fn g4_cancellation_is_structured_and_publishes_nothing() {
    let fiducials = exact_fiducials();
    let system = block_system();
    assert!(matches!(
        with_cx(true, ExecMode::Deterministic, Budget::INFINITE, |cx| {
            register3_datum(&fiducials, &system, cx)
        }),
        Err(DatumError::Cancelled {
            phase: "datum.a-design"
        })
    ));
}
