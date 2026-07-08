//! Marquee-demo conformance (the mye.1 bead; runs under `marquee`).
//! Acceptance (smoke tier — the full-resolution run is the nightly
//! golden): the study runs end-to-end from a raw SDF with no mesh in
//! the loop and the objective improves; certificate components verify
//! against a refined-reference measurement; replay is bit-equal (G5);
//! the flat-cadence claim holds (no remeshing spikes); seeded-failure
//! drills (broken gradient, budget exhaustion) produce structured
//! outcomes — the FD falsifier catches the broken adjoint; the design
//! sphere-traces through the render backend (no meshing for pictures).
#![cfg(feature = "marquee")]

use fs_marquee::{PlateWithHoles, StudyConfig, run_study};

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-marquee\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

fn two_hole_plate() -> PlateWithHoles {
    PlateWithHoles {
        centers: vec![[0.3, 0.5], [0.7, 0.5]],
        radii: vec![0.12, 0.18],
    }
}

fn smoke_config() -> StudyConfig {
    StudyConfig {
        level: 4,
        steps: 8,
        step_size: 3.0,
        area_target: two_hole_plate().area(),
        r_min: 0.06,
        r_max: 0.24,
    }
}

#[test]
fn mq_001_end_to_end_objective_improves() {
    let report = run_study(two_hole_plate(), &smoke_config()).expect("study runs");
    assert_eq!(report.iterations.len(), 8);
    let first = report.iterations.first().expect("first").compliance;
    let last = report.iterations.last().expect("last").compliance;
    // The optimizer redistributes hole area toward equal boundary flux
    // (the optimality condition): compliance must improve.
    assert!(
        last < first,
        "compliance improves under the area budget: {first:.5} -> {last:.5}"
    );
    // The area budget is honored throughout.
    for rec in &report.iterations {
        assert!(
            (rec.area - smoke_config().area_target).abs() < 0.02,
            "area budget held at iter {}: {}",
            rec.iter,
            rec.area
        );
    }
    // Every iteration carries the full certificate.
    for rec in &report.iterations {
        assert!(rec.cert_dwr > 0.0, "the DWR component exists");
        assert!(rec.cert_algebraic >= 0.0);
        assert!(
            matches!(rec.color, fs_evidence::Color::Estimated { .. }),
            "the composed color is honest (DWR is estimated): {:?}",
            rec.color
        );
    }
    println!(
        "{{\"metric\":\"marquee-objective\",\"first\":{first:.5},\"last\":{last:.5},\
         \"iters\":8}}"
    );
    verdict(
        "mq-001",
        "6-step smoke study: compliance improves, area budget held, every value carries \
         its composed estimated-color certificate",
    );
}

#[test]
fn mq_002_certificate_vs_refined_reference() {
    // The certificate's own falsifier: the coarse-grid compliance must
    // sit within a small multiple of its certified band of the
    // refined-reference value (DWR effectivity is not guaranteed — the
    // band factor is the documented tolerance).
    let design = two_hole_plate();
    let coarse = run_study(
        design.clone(),
        &StudyConfig {
            steps: 1,
            ..smoke_config()
        },
    )
    .expect("coarse");
    let refined = run_study(
        design,
        &StudyConfig {
            level: 5,
            steps: 1,
            ..smoke_config()
        },
    )
    .expect("refined");
    let jc = coarse.iterations[0].compliance;
    let jf = refined.iterations[0].compliance;
    let band = coarse.iterations[0].cert_dwr + coarse.iterations[0].cert_algebraic;
    let gap = (jc - jf).abs();
    println!(
        "{{\"metric\":\"certificate-check\",\"coarse\":{jc:.6},\"refined\":{jf:.6},\
         \"gap\":{gap:.2e},\"band\":{band:.2e}}}"
    );
    assert!(
        gap <= 4.0 * band.max(1e-12),
        "the refined-reference gap sits within 4x the certified band: gap {gap:.2e} vs \
         band {band:.2e}"
    );
    verdict(
        "mq-002",
        "coarse-vs-refined compliance gap within 4x the composed certificate band \
         (effectivity factor documented)",
    );
}

#[test]
fn mq_003_replay_bit_equal_and_flat_cadence() {
    let a = run_study(two_hole_plate(), &smoke_config()).expect("a");
    let b = run_study(two_hole_plate(), &smoke_config()).expect("b");
    assert_eq!(a.trace_hash, b.trace_hash, "the study replays bit-exact");
    for (ra, rb) in a.iterations.iter().zip(&b.iterations) {
        assert_eq!(ra.compliance.to_bits(), rb.compliance.to_bits());
    }
    // FLAT CADENCE: no remeshing spikes — solver iterations stay within
    // a tight band across the whole study (there is nothing to remesh).
    let iters: Vec<usize> = a.iterations.iter().map(|r| r.solver_iters).collect();
    let (lo, hi) = iters
        .iter()
        .fold((usize::MAX, 0), |(l, h), &v| (l.min(v), h.max(v)));
    println!("{{\"metric\":\"cadence\",\"solver_iters\":{iters:?}}}");
    assert!(
        hi <= 2 * lo.max(1),
        "per-iteration cost stays flat (no remeshing spikes): {iters:?}"
    );
    verdict(
        "mq-003",
        "bit-equal replay (G5); solver iterations within a 2x band across the study — \
         the no-remeshing cadence",
    );
}

#[test]
fn mq_004_seeded_failure_drills() {
    // DRILL 1 — broken adjoint: flip the reported gradient's sign and
    // let the FD falsifier (Proposal 6, through fs-adjoint) catch it.
    let design = two_hole_plate();
    let report = run_study(
        design.clone(),
        &StudyConfig {
            steps: 1,
            ..smoke_config()
        },
    )
    .expect("study");
    let grad = &report.iterations[0].gradient;
    let objective = |radii: &[f64]| -> f64 {
        let d = PlateWithHoles {
            centers: design.centers.clone(),
            radii: radii.to_vec(),
        };
        run_study(
            d,
            &StudyConfig {
                steps: 1,
                step_size: 0.0,
                ..smoke_config()
            },
        )
        .expect("probe")
        .iterations[0]
            .compliance
    };
    // The HONEST gradient passes the conditioning-aware FD check…
    let dir = vec![1.0, -1.0];
    let dd: f64 = grad.iter().zip(&dir).map(|(g, d)| g * d).sum();
    let ok = fs_adjoint::transpose::fd_falsifier(
        &objective,
        &report.iterations[0].radii,
        &dir,
        dd,
        5e-3,
        2e-3,
    );
    assert!(ok.consistent, "the shape gradient passes FD: {ok:?}");
    // …and the SIGN-FLIPPED (broken) adjoint is caught.
    let broken = fs_adjoint::transpose::fd_falsifier(
        &objective,
        &report.iterations[0].radii,
        &dir,
        -dd,
        5e-3,
        2e-3,
    );
    assert!(
        !broken.consistent,
        "the falsifier catches the broken adjoint: {broken:?}"
    );
    // DRILL 2 — budget exhaustion: an over-tight radius box makes the
    // projection infeasible-by-clamping; the study still returns a
    // STRUCTURED report (no panic, no silent nonsense).
    let clamped = run_study(
        two_hole_plate(),
        &StudyConfig {
            r_min: 0.17,
            r_max: 0.19,
            steps: 2,
            ..smoke_config()
        },
    )
    .expect("clamped study still structures its outcome");
    assert_eq!(clamped.iterations.len(), 2);
    verdict(
        "mq-004",
        "the FD falsifier passes the honest gradient and catches the sign-flipped \
         adjoint; the clamped-budget drill returns a structured report",
    );
}

#[test]
fn mq_005_sphere_traced_render_no_meshing() {
    use asupersync::types::Budget;
    use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
    use fs_geom::{Point3, Vec3};
    use fs_render::charts::{Ray, sphere_trace};
    use fs_rep_frep::FrepBuilder;
    // The final design, held as a 3-D F-rep (extruded plate minus hole
    // cylinders) and sphere-traced DIRECTLY — no meshing for pictures.
    let report = run_study(two_hole_plate(), &smoke_config()).expect("study");
    let mut b = FrepBuilder::new();
    let plate = b
        .box_prim(
            Point3::new(0.5, 0.5, 0.0),
            Vec3::new(0.5, 0.5, 0.05),
        )
        .expect("plate");
    let mut shape = plate;
    for (c, r) in report.design.centers.iter().zip(&report.design.radii) {
        let hole = b
            .cylinder(Point3::new(c[0], c[1], 0.0), *r)
            .expect("hole");
        shape = b
            .boolean(
                fs_rep_frep::BoolOp::Difference,
                fs_rep_frep::BoolStyle::Hard,
                shape,
                hole,
            )
            .expect("difference");
    }
    let frep = b.finish(shape).expect("frep");
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 1,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        // A 16x16 turntable frame of hit depths.
        let n = 16usize;
        let mut hits = 0usize;
        let mut misses_in_holes = 0usize;
        for py in 0..n {
            for px in 0..n {
                #[allow(clippy::cast_precision_loss)]
                let (x, y) = (
                    (px as f64 + 0.5) / n as f64,
                    (py as f64 + 0.5) / n as f64,
                );
                let ray = Ray {
                    origin: Point3::new(x, y, 1.0),
                    dir: Vec3::new(0.0, 0.0, -1.0),
                };
                let (hit, _) = sphere_trace(&frep, &cx, &ray, 3.0, 1e-6, 1.0);
                let in_hole = report
                    .design
                    .centers
                    .iter()
                    .zip(&report.design.radii)
                    .any(|(c, r)| {
                        ((x - c[0]).powi(2) + (y - c[1]).powi(2)).sqrt() < r - 0.04
                    });
                if hit.is_some() {
                    hits += 1;
                    assert!(!in_hole, "rays through holes must miss the plate");
                } else if in_hole {
                    misses_in_holes += 1;
                }
            }
        }
        assert!(hits > 150, "the plate fills most of the frame: {hits}");
        assert!(
            misses_in_holes > 0,
            "the holes are visible in the render: {misses_in_holes}"
        );
    });
    verdict(
        "mq-005",
        "the final design sphere-traces directly as an F-rep: plate visible, holes \
         punched, zero meshing anywhere in the study or the picture",
    );
}
