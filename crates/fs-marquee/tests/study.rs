//! Marquee STUDY conformance (the mye.1 bead; runs under `marquee`).
//! Acceptance (smoke tier — the full-resolution run is the nightly
//! golden): the study runs end-to-end from a raw SDF with no mesh in
//! the loop and the objective improves; certificate components verify
//! against a refined-reference measurement; replay is bit-equal (G5);
//! the flat-cadence claim holds (no remeshing spikes); seeded-failure
//! drills (broken gradient, budget exhaustion) produce structured
//! outcomes — the FD falsifier catches the broken adjoint; the design
//! sphere-traces through the render backend (no meshing for pictures).
#![cfg(feature = "marquee")]

use fs_marquee::study::{AREA_PROJECTION_TOLERANCE, PlateWithHoles, StudyConfig, run_study};

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-marquee/study\",\"case\":\"{case}\",\"verdict\":\"pass\",\
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
        step_size: 1.0,
        area_target: two_hole_plate().area(),
        r_min: 0.08,
        r_max: 0.20,
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
        "compliance improves under the area budget: {first:.6} -> {last:.6}"
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
        assert!(
            rec.cert_algebraic.is_finite() && rec.cert_algebraic >= 0.0,
            "the algebraic term comes from an admitted recomputed Euclidean residual"
        );
        assert!(
            matches!(rec.color, fs_evidence::Color::Estimated { .. }),
            "the composed color is honest (DWR is estimated): {:?}",
            rec.color
        );
    }
    println!(
        "{{\"metric\":\"marquee-objective\",\"first\":{first:.6},\"last\":{last:.6},\
         \"iters\":8}}"
    );
    verdict(
        "mq-001",
        "8-step smoke study: compliance improves, area budget held, every value carries \
         its composed estimated-color certificate",
    );
}

#[test]
fn mq_001b_bound_active_projection_keeps_the_entire_trajectory_feasible() {
    let config = StudyConfig {
        r_min: 0.14,
        r_max: 0.16,
        steps: 3,
        ..smoke_config()
    };
    let report = run_study(two_hole_plate(), &config).expect("bound-active study runs");
    let target = config.area_target;

    assert_eq!(report.iterations.len(), 3);
    assert!(
        report.iterations[0]
            .radii
            .iter()
            .any(|radius| (*radius - config.r_max).abs() < 1e-12),
        "the fixture activates a radius bound"
    );
    for record in &report.iterations {
        assert!(
            (record.area - target).abs() <= AREA_PROJECTION_TOLERANCE,
            "iteration {} honors the declared material-area equality: {} vs {target}",
            record.iter,
            record.area
        );
        assert!(
            record
                .radii
                .iter()
                .all(|radius| *radius >= config.r_min && *radius <= config.r_max),
            "iteration {} stays within the radius box: {:?}",
            record.iter,
            record.radii
        );
    }
    let final_record = report.iterations.last().expect("nonempty trajectory");
    assert_eq!(report.design.radii, final_record.accepted_radii);
    assert!(
        (report.design.area() - target).abs() <= AREA_PROJECTION_TOLERANCE,
        "final accepted design honors the declared material-area equality"
    );
    verdict(
        "mq-001b",
        "bound-active three-step trajectory preserves every radius bound and material area \
         within AREA_PROJECTION_TOLERANCE",
    );
}

#[test]
fn mq_001c_armijo_records_bounded_acceptance_for_an_oversized_proposal() {
    let config = StudyConfig {
        step_size: 64.0,
        steps: 3,
        ..smoke_config()
    };
    let report = run_study(two_hole_plate(), &config).expect("Armijo study runs");

    assert_eq!(report.iterations.len(), 3);
    for record in &report.iterations {
        assert!(
            record.backtracks <= 8,
            "the retry budget is bounded: {:?}",
            record
        );
        if record.accepted_step == 0.0 {
            assert_eq!(record.backtracks, 8, "only exhaustion retains the design");
        } else {
            let expected = config.step_size * 0.5_f64.powi(record.backtracks as i32);
            assert_eq!(
                record.accepted_step, expected,
                "the recorded step is the deterministically accepted trial"
            );
        }
    }
    for pair in report.iterations.windows(2) {
        assert_eq!(
            pair[0].accepted_radii, pair[1].radii,
            "each accepted successor is the next recorded design"
        );
        assert_eq!(
            pair[0].accepted_compliance.to_bits(),
            pair[1].compliance.to_bits(),
            "each accepted successor objective is retained in the next record"
        );
        assert!(
            pair[1].compliance <= pair[0].compliance + AREA_PROJECTION_TOLERANCE,
            "accepted Armijo trajectory does not increase compliance: {} -> {}",
            pair[0].compliance,
            pair[1].compliance
        );
    }
    let final_record = report.iterations.last().expect("nonempty trajectory");
    assert_eq!(report.design.radii, final_record.accepted_radii);
    assert!(
        final_record.accepted_compliance <= final_record.compliance + AREA_PROJECTION_TOLERANCE,
        "the final design carries an Armijo-accepted successor objective"
    );
    verdict(
        "mq-001c",
        "oversized proposal records a bounded Armijo acceptance decision instead of an \
         unconditional radius update",
    );
}

#[test]
fn mq_001d_refuses_extreme_projection_ratios() {
    for (r_min, r_max) in [(1e-300, 1.0), (1e-320, 1e308)] {
        let config = StudyConfig {
            r_min,
            r_max,
            area_target: 0.9,
            ..smoke_config()
        };
        assert!(
            std::panic::catch_unwind(|| run_study(two_hole_plate(), &config)).is_err(),
            "ratio {r_max} / {r_min} must be refused before projection"
        );
    }
}

#[test]
fn mq_001e_refuses_overlapping_or_boundary_clipped_holes() {
    let overlap = PlateWithHoles {
        centers: vec![[0.45, 0.5], [0.55, 0.5]],
        radii: vec![0.1, 0.1],
    };
    assert!(
        std::panic::catch_unwind(|| run_study(overlap, &smoke_config())).is_err(),
        "overlapping disks must be refused before summed-area optimization"
    );

    let boundary_clipped = PlateWithHoles {
        centers: vec![[0.1, 0.5], [0.7, 0.5]],
        radii: vec![0.1, 0.1],
    };
    assert!(
        std::panic::catch_unwind(|| run_study(boundary_clipped, &smoke_config())).is_err(),
        "boundary-clipped disks must be refused before summed-area optimization"
    );
}

#[test]
fn mq_001f_refuses_a_candidate_that_grows_into_overlap() {
    let initially_valid = PlateWithHoles {
        centers: vec![[0.35, 0.5], [0.65, 0.5]],
        radii: vec![0.1, 0.1],
    };
    let config = StudyConfig {
        r_min: 0.05,
        r_max: 0.2,
        area_target: 1.0 - std::f64::consts::PI * (0.2_f64.powi(2) + 0.2_f64.powi(2)),
        ..smoke_config()
    };
    assert!(
        std::panic::catch_unwind(|| run_study(initially_valid, &config)).is_err(),
        "post-projection overlap must be refused before the first solve"
    );
}

#[test]
fn mq_001h_invalid_armijo_trial_backtracks_without_panicking_the_study() {
    let design = PlateWithHoles {
        centers: vec![[0.35, 0.5], [0.65, 0.5]],
        radii: vec![0.05, 0.24],
    };
    let config = StudyConfig {
        r_min: 0.05,
        r_max: 0.25,
        step_size: f64::MAX,
        steps: 1,
        area_target: design.area(),
        ..smoke_config()
    };
    let report = run_study(design, &config).expect("invalid trial is rejected, not fatal");
    let record = report.iterations.first().expect("one iteration");

    assert!(record.backtracks > 0, "the oversized trial was rejected");
    assert!(
        record.accepted_radii[0] + record.accepted_radii[1] < 0.3,
        "the retained design remains disjoint after trial rejection"
    );
    assert_eq!(report.design.radii, record.accepted_radii);
}

#[test]
fn mq_001g_iteration_jsonl_rows_are_complete_and_deterministic() {
    let report = run_study(
        two_hole_plate(),
        &StudyConfig {
            steps: 1,
            ..smoke_config()
        },
    )
    .expect("study runs");
    let record = report.iterations.first().expect("one row");
    let row = record.jsonl_row();

    assert_eq!(
        row,
        record.jsonl_row(),
        "JSONL serialization is deterministic"
    );
    assert!(row.ends_with('\n'), "a JSONL row has a newline delimiter");
    for field in [
        "\"iter\":",
        "\"compliance\":",
        "\"volume\":",
        "\"gradient_norm\":",
        "\"dwr_estimate\":",
        "\"cut_cell_count\":",
        "\"accepted_step\":",
        "\"backtracks\":",
        "\"accepted_compliance\":",
        "\"color_rank\":",
        "\"color_payload\":",
    ] {
        assert!(row.contains(field), "JSONL row retains {field}");
    }
    assert!(record.gradient_norm.is_finite() && record.gradient_norm >= 0.0);
    assert_eq!(record.volume.to_bits(), record.area.to_bits());
    assert_eq!(record.dwr_estimate.to_bits(), record.cert_dwr.to_bits());
    assert!(
        record.cut_cell_count > 0,
        "the CutFEM solve retained cut cells"
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
    let short = StudyConfig {
        steps: 4,
        ..smoke_config()
    };
    let a = run_study(two_hole_plate(), &short).expect("a");
    let b = run_study(two_hole_plate(), &short).expect("b");
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
            r_min: 0.14,
            r_max: 0.16,
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
    let report = run_study(
        two_hole_plate(),
        &StudyConfig {
            steps: 2,
            ..smoke_config()
        },
    )
    .expect("study");
    let mut b = FrepBuilder::new();
    let plate = b
        .box_prim(Point3::new(0.5, 0.5, 0.0), Vec3::new(0.5, 0.5, 0.05))
        .expect("plate");
    let mut shape = plate;
    for (c, r) in report.design.centers.iter().zip(&report.design.radii) {
        let hole = b.cylinder(Point3::new(c[0], c[1], 0.0), *r).expect("hole");
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
                let (x, y) = ((px as f64 + 0.5) / n as f64, (py as f64 + 0.5) / n as f64);
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
                    .any(|(c, r)| ((x - c[0]).powi(2) + (y - c[1]).powi(2)).sqrt() < r - 0.04);
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
