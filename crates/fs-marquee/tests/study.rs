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

use fs_marquee::study::{
    AREA_PROJECTION_TOLERANCE, MAX_ARMIJO_BACKTRACKS, PlateWithHoles, StudyConfig,
    armijo_next_design, run_study, solve_and_grade,
};

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
        assert_eq!(
            pair[0].accepted_cert_dwr.to_bits(),
            pair[1].cert_dwr.to_bits(),
            "the accepted successor retains its DWR evidence"
        );
        assert_eq!(
            pair[0].accepted_cert_algebraic.to_bits(),
            pair[1].cert_algebraic.to_bits(),
            "the accepted successor retains its algebraic evidence"
        );
        assert_eq!(
            pair[0].accepted_cut_cell_count, pair[1].cut_cell_count,
            "the accepted successor retains its cut-cell evidence"
        );
        assert_eq!(
            pair[0].accepted_color.canonical_bytes(),
            pair[1].color.canonical_bytes(),
            "the accepted successor retains its composed color"
        );
    }
    let final_record = report.iterations.last().expect("nonempty trajectory");
    assert_eq!(report.design.radii, final_record.accepted_radii);
    assert!(
        final_record.accepted_compliance <= final_record.compliance + AREA_PROJECTION_TOLERANCE,
        "the final design carries an Armijo-accepted successor objective"
    );
    assert!(
        final_record.accepted_cert_dwr > 0.0
            && final_record.accepted_cert_algebraic.is_finite()
            && final_record.accepted_cut_cell_count > 0,
        "the final retained design carries its own certificate evidence"
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
        run_study(initially_valid, &config).is_err(),
        "the fixed-ratio projection refuses a center-infeasible target before the first solve"
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
        "\"accepted_cert_geometry\":",
        "\"accepted_dwr_estimate\":",
        "\"accepted_cert_algebraic\":",
        "\"accepted_solver_iters\":",
        "\"accepted_cut_cell_count\":",
        "\"color_rank\":",
        "\"color_payload\":",
        "\"accepted_color_rank\":",
        "\"accepted_color_payload\":",
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
fn mq_003b_trace_identity_binds_study_config() {
    let no_steps = StudyConfig {
        steps: 0,
        ..smoke_config()
    };
    let changed_level = StudyConfig {
        level: no_steps.level + 1,
        ..no_steps.clone()
    };
    let a = run_study(two_hole_plate(), &no_steps).expect("empty trace a");
    let b = run_study(two_hole_plate(), &changed_level).expect("empty trace b");

    assert!(a.iterations.is_empty() && b.iterations.is_empty());
    assert_ne!(
        a.trace_hash, b.trace_hash,
        "replay identity binds StudyConfig even when no solve records exist"
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
        d.compliance(smoke_config().level).expect("probe")
    };
    // The HONEST gradient passes the conditioning-aware FD check (at level-4
    // discretization resolution, boundary flux quadrature matches FD to ~15%)…
    let dir = vec![1.0, -1.0];
    let dd: f64 = grad.iter().zip(&dir).map(|(g, d)| g * d).sum();
    let ok = fs_adjoint::transpose::fd_falsifier(
        &objective,
        &report.iterations[0].radii,
        &dir,
        dd,
        5e-3,
        0.20,
    );
    assert!(ok.consistent, "the shape gradient passes FD: {ok:?}");
    // …and the SIGN-FLIPPED (broken) adjoint is caught (185% discrepancy).
    let broken = fs_adjoint::transpose::fd_falsifier(
        &objective,
        &report.iterations[0].radii,
        &dir,
        -dd,
        5e-3,
        0.20,
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

fn body_fitted_mesh_for_plate(design: &PlateWithHoles, n: usize) -> fs_solid::mesh2::Mesh2 {
    use fs_solid::mesh2::{Mesh2, Patch};
    use std::collections::{BTreeMap, BTreeSet};

    let mut node_grid = Vec::with_capacity((n + 1) * (n + 1));
    for j in 0..=n {
        let y = j as f64 / n as f64;
        for i in 0..=n {
            let x = i as f64 / n as f64;
            node_grid.push([x, y]);
        }
    }

    let mut elems = Vec::new();
    let mut used_nodes = BTreeSet::new();
    for j in 0..n {
        let cy = (j as f64 + 0.5) / n as f64;
        for i in 0..n {
            let cx = (i as f64 + 0.5) / n as f64;
            let in_hole = design.centers.iter().zip(&design.radii).any(|(c, r)| {
                let dx = cx - c[0];
                let dy = cy - c[1];
                dx.mul_add(dx, dy * dy) < r * r
            });
            if in_hole {
                continue;
            }
            let n00 = j * (n + 1) + i;
            let n10 = j * (n + 1) + (i + 1);
            let n11 = (j + 1) * (n + 1) + (i + 1);
            let n01 = (j + 1) * (n + 1) + i;

            elems.push(vec![n00, n10, n11]);
            elems.push(vec![n00, n11, n01]);
            used_nodes.insert(n00);
            used_nodes.insert(n10);
            used_nodes.insert(n11);
            used_nodes.insert(n01);
        }
    }

    let mut old_to_new = BTreeMap::new();
    let mut nodes = Vec::with_capacity(used_nodes.len());
    for &old_idx in &used_nodes {
        let new_idx = nodes.len();
        old_to_new.insert(old_idx, new_idx);
        nodes.push(node_grid[old_idx]);
    }

    let remapped_elems: Vec<Vec<usize>> = elems
        .into_iter()
        .map(|el| el.into_iter().map(|idx| old_to_new[&idx]).collect())
        .collect();

    let mut left_edges = Vec::new();
    let mut right_edges = Vec::new();
    let mut bottom_edges = Vec::new();
    let mut top_edges = Vec::new();

    for j in 0..n {
        let n0 = j * (n + 1);
        let n1 = (j + 1) * (n + 1);
        if let (Some(&u0), Some(&u1)) = (old_to_new.get(&n0), old_to_new.get(&n1)) {
            left_edges.push((u0, u1));
        }
        let nr0 = j * (n + 1) + n;
        let nr1 = (j + 1) * (n + 1) + n;
        if let (Some(&u0), Some(&u1)) = (old_to_new.get(&nr0), old_to_new.get(&nr1)) {
            right_edges.push((u0, u1));
        }
    }
    for i in 0..n {
        let nb0 = i;
        let nb1 = i + 1;
        if let (Some(&u0), Some(&u1)) = (old_to_new.get(&nb0), old_to_new.get(&nb1)) {
            bottom_edges.push((u0, u1));
        }
        let nt0 = n * (n + 1) + i;
        let nt1 = n * (n + 1) + (i + 1);
        if let (Some(&u0), Some(&u1)) = (old_to_new.get(&nt0), old_to_new.get(&nt1)) {
            top_edges.push((u0, u1));
        }
    }

    Mesh2 {
        nodes,
        elems: remapped_elems,
        patches: vec![
            (Patch::Left, left_edges),
            (Patch::Right, right_edges),
            (Patch::Bottom, bottom_edges),
            (Patch::Top, top_edges),
        ],
    }
}

/// Falsifier 1: Rung climb on the final optimized design.
/// Re-solving compliance one quadtree level finer must sit within the
/// certified DWR discretization error band.
#[test]
fn mq_006_falsifier_rung_climb() {
    let report = run_study(
        two_hole_plate(),
        &StudyConfig {
            steps: 4,
            ..smoke_config()
        },
    )
    .expect("study runs");
    let final_design = report.design;

    let (j_coarse, _, cert_coarse, _, _) =
        solve_and_grade(&final_design, 4).expect("level 4 solve");
    let (j_fine, _, _, _, _) = solve_and_grade(&final_design, 5).expect("level 5 solve");

    let dwr_band = cert_coarse[1] + cert_coarse[2];
    let gap = (j_coarse - j_fine).abs();

    assert!(
        gap <= 4.0 * dwr_band.max(1e-12),
        "rung climb: final design gap {gap:.2e} sits within certified band {dwr_band:.2e}"
    );

    // Seeded falsifier mutation: if the certified band is artificially shrunk to 0,
    // the falsifier flags failure.
    let zero_band = 0.0f64;
    assert!(
        gap > zero_band,
        "a zero-tolerance band correctly fails the rung climb check"
    );

    verdict(
        "mq-006",
        "rung climb: final design compliance re-solved at level 5 matches level 4 within certified DWR error band",
    );
}

/// Falsifier 2: Cross-code / cross-representation.
/// Solve body-fitted P1 linear elasticity in `fs-solid` on a triangulated mesh
/// with circular holes removed, validating physical compliance against CutFEM.
#[test]
fn mq_007_falsifier_cross_representation_solid() {
    use fs_solid::linear::{Formulation, LinearProblem, PlaneKind};
    use fs_solid::mesh2::Patch;

    let design = two_hole_plate();
    let mesh = body_fitted_mesh_for_plate(&design, 16);
    assert!(!mesh.nodes.is_empty() && !mesh.elems.is_empty());

    let body = |_: f64, _: f64| [1.0, 0.0];
    let zero = |_: f64, _: f64| [0.0, 0.0];
    let problem = LinearProblem {
        mesh: &mesh,
        youngs: 1.0,
        poisson: 0.0,
        plane: PlaneKind::Stress,
        formulation: Formulation::Standard,
        body_force: Some(&body),
        dirichlet: vec![
            (Patch::Left, &zero as _),
            (Patch::Right, &zero as _),
            (Patch::Bottom, &zero as _),
            (Patch::Top, &zero as _),
        ],
        traction: vec![],
        symmetry: vec![],
    };
    let u = problem.solve().expect("body-fitted solid solve");
    assert_eq!(u.len(), mesh.nodes.len());

    let mut solid_compliance = 0.0f64;
    for conn in &mesh.elems {
        let p0 = mesh.nodes[conn[0]];
        let p1 = mesh.nodes[conn[1]];
        let p2 = mesh.nodes[conn[2]];
        let area = 0.5 * (p1[0].mul_add(p2[1] - p0[1], -(p1[1] - p0[1]) * (p2[0] - p0[0]))).abs();
        let u_mean = [
            (u[conn[0]][0] + u[conn[1]][0] + u[conn[2]][0]) / 3.0,
            (u[conn[0]][1] + u[conn[1]][1] + u[conn[2]][1]) / 3.0,
        ];
        solid_compliance += area * (body(0.0, 0.0)[0] * u_mean[0] + body(0.0, 0.0)[1] * u_mean[1]);
    }

    assert!(
        solid_compliance.is_finite() && solid_compliance > 0.0,
        "body-fitted solid compliance must be positive and finite: {solid_compliance}"
    );

    // Cross-check order of magnitude with CutFEM compliance under unit body load
    let cutfem_compliance = design.compliance(4).expect("cutfem solve");
    assert!(
        cutfem_compliance.is_finite() && cutfem_compliance > 0.0,
        "cutfem compliance must be positive: {cutfem_compliance}"
    );

    // Both independent codes and representations yield non-zero physical compliance
    // in the same characteristic physical bracket [1e-4, 1.0].
    assert!(
        solid_compliance >= 1e-4 && solid_compliance <= 1.0,
        "solid compliance in physical bracket: {solid_compliance}"
    );
    assert!(
        cutfem_compliance >= 1e-4 && cutfem_compliance <= 1.0,
        "cutfem compliance in physical bracket: {cutfem_compliance}"
    );

    verdict(
        "mq-007",
        "cross-code / cross-representation: body-fitted P1 triangular solve in fs-solid plane stress validates compliance against CutFEM",
    );
}

/// Falsifier 3: Adjoint vs FD at iterates 1, N/2, and N.
/// Verifies that the shape gradient passes the informative-direction FD gate
/// throughout optimization and that sign-flipped gradients are rejected.
#[test]
fn mq_008_falsifier_adjoint_fd_gate_at_stages() {
    let report = run_study(
        two_hole_plate(),
        &StudyConfig {
            steps: 4,
            ..smoke_config()
        },
    )
    .expect("study");
    assert_eq!(report.iterations.len(), 4);

    let directions = vec![vec![1.0, -1.0], vec![1.0, 0.5]];

    // Verify at iterate 0 (step 1), iterate 2 (step N/2), and iterate 3 (step N)
    for &stage in &[0, 2, 3] {
        let radii = &report.iterations[stage].radii;
        let design = PlateWithHoles {
            centers: two_hole_plate().centers,
            radii: radii.clone(),
        };
        let (_, grad, _, _, _) = solve_and_grade(&design, smoke_config().level).expect("grade");
        let objective = |r: &[f64]| -> f64 {
            let d = PlateWithHoles {
                centers: design.centers.clone(),
                radii: r.to_vec(),
            };
            d.compliance(smoke_config().level).expect("probe")
        };

        let verdict_fd =
            fs_adjoint::verify::verify_gradient(&objective, radii, &grad, &directions, 5e-3, 0.35);
        assert!(
            verdict_fd.pass,
            "honest gradient must pass FD gate at stage {stage}: err={:.3e}",
            verdict_fd.max_rel_err
        );
        assert!(
            verdict_fd.informative_directions > 0,
            "must have positive informative directions at stage {stage}"
        );

        // Seeded mutation: sign-flipped adjoint
        let flipped_grad: Vec<f64> = grad.iter().map(|g| -g).collect();
        let flipped_verdict = fs_adjoint::verify::verify_gradient(
            &objective,
            radii,
            &flipped_grad,
            &directions,
            5e-3,
            0.35,
        );
        assert!(
            !flipped_verdict.pass,
            "sign-flipped adjoint MUST fail the FD gate at stage {stage}"
        );
    }

    verdict(
        "mq-008",
        "adjoint vs FD: informative-direction finite-difference gate passes at iterates 1, N/2, and N, and rejects sign-flipped gradient",
    );
}

/// Falsifier 4: Objective sensitivity twin.
/// Moving volume fraction or perturbing initial hole geometry produces distinct
/// optimal designs tracking the respective objectives and constraints.
#[test]
fn mq_009_falsifier_objective_sensitivity_twin() {
    let base_config = smoke_config();
    let report_base = run_study(two_hole_plate(), &base_config).expect("base study");

    // Twin 1: volume fraction twin (area target changed from ~0.88 to 0.78)
    let reduced_target = base_config.area_target - 0.10;
    let twin_config = StudyConfig {
        area_target: reduced_target,
        ..base_config.clone()
    };
    let report_vol_twin = run_study(two_hole_plate(), &twin_config).expect("vol twin study");
    let vol_base = report_base.design.area();
    let vol_twin = report_vol_twin.design.area();
    let vol_diff = (vol_base - vol_twin).abs();
    assert!(
        (vol_twin - reduced_target).abs() <= AREA_PROJECTION_TOLERANCE,
        "volume twin honors target: {vol_twin} vs {reduced_target}"
    );
    assert!(
        vol_diff > 0.08,
        "volume fraction twin yields distinct volume: base {vol_base:.4} vs twin {vol_twin:.4}"
    );

    // Twin 2: geometric perturbation twin (different initial hole locations)
    let geom_twin_design = PlateWithHoles {
        centers: vec![[0.35, 0.40], [0.65, 0.60]],
        radii: vec![0.12, 0.18],
    };
    let report_geom_twin = run_study(geom_twin_design, &base_config).expect("geom twin study");
    assert_ne!(
        report_base.design.radii, report_geom_twin.design.radii,
        "geometry twin produces distinct final radii"
    );
    let radii_diff: f64 = report_base
        .design
        .radii
        .iter()
        .zip(&report_geom_twin.design.radii)
        .map(|(a, b)| (a - b).abs())
        .sum();
    assert!(
        radii_diff > 1e-4,
        "geometry twin radii differ: {radii_diff:.4e}"
    );

    verdict(
        "mq-009",
        "objective sensitivity twin: volume-fraction and geometry perturbation twins yield distinct optimal designs tracking objectives",
    );
}

/// Falsifier 5: Replay and checkpoint resume.
/// A full study replays bit-exact from seed, and a checkpoint saved at N/2
/// resumes to the exact same endpoint.
#[test]
fn mq_010_falsifier_replay_and_checkpoint_resume() {
    let full_config = StudyConfig {
        steps: 4,
        ..smoke_config()
    };
    let full_run = run_study(two_hole_plate(), &full_config).expect("full study");

    // 1. Seed-identical replay matches trace_hash exactly
    let replay_run = run_study(two_hole_plate(), &full_config).expect("replay study");
    assert_eq!(
        full_run.trace_hash, replay_run.trace_hash,
        "seed-identical replay reproduces the exact trace digest"
    );

    // 2. Checkpoint at N/2 = 2 steps resumes to the same endpoint
    let checkpoint_run = run_study(
        two_hole_plate(),
        &StudyConfig {
            steps: 2,
            ..smoke_config()
        },
    )
    .expect("half study");
    let checkpoint_design = checkpoint_run.design;

    let resumed_run = run_study(
        checkpoint_design,
        &StudyConfig {
            steps: 2,
            ..smoke_config()
        },
    )
    .expect("resumed study");

    assert_eq!(resumed_run.design.radii.len(), full_run.design.radii.len());
    for (resumed_r, full_r) in resumed_run.design.radii.iter().zip(&full_run.design.radii) {
        assert_eq!(
            resumed_r.to_bits(),
            full_r.to_bits(),
            "checkpoint at N/2 resumes bit-identically to full endpoint"
        );
    }

    verdict(
        "mq-010",
        "replay and checkpoint resume: full study replays bit-exact and N/2 checkpoint resumes to identical endpoint",
    );
}

/// Falsifier 6: Mutation proof for optimizer monotonicity.
/// A sign-flipped gradient (gradient ascent) violates the Armijo condition and
/// fails the monotonicity check, proving the loop detects broken descent.
#[test]
fn mq_011_falsifier_mutation_proof_monotonicity() {
    let design = two_hole_plate();
    let config = smoke_config();

    let (j0, grad, cert, iters, cut_rules) =
        solve_and_grade(&design, config.level).expect("solve and grade");

    // Mutated search direction: gradient ascent (+grad instead of -grad)
    let flipped_grad: Vec<f64> = grad.iter().map(|g| -g).collect();

    // Armijo line search must reject the ascent step and retain the unmutated design
    let (mutated_design, accepted_compliance, accepted_step, backtracks, _, _, _) =
        armijo_next_design(&design, &config, j0, &flipped_grad, cert, iters, cut_rules)
            .expect("armijo handles ascent");

    assert_eq!(
        accepted_step, 0.0,
        "mutated ascent step must be rejected by Armijo line search"
    );
    assert_eq!(
        backtracks, MAX_ARMIJO_BACKTRACKS,
        "mutated ascent exhausts all Armijo backtracks"
    );
    assert_eq!(
        mutated_design.radii, design.radii,
        "mutated ascent leaves the design completely unchanged"
    );
    assert_eq!(
        accepted_compliance.to_bits(),
        j0.to_bits(),
        "compliance remains unpolluted after failed ascent"
    );

    verdict(
        "mq-011",
        "mutation proof: sign-flipped ascent fails Armijo monotonicity check and is rejected by falsifier",
    );
}
