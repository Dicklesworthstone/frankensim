//! Exercise the actual command, including its resolved refined-request replay.
use super::*;

fn input(adaptive: bool) -> J {
    let mut root = J::parse(BASE).unwrap();
    put(member(&mut root, "objective"), "gradient", J::Bool(false));
    nonlinear(&mut root);
    radiation(&mut root);
    // Explicit study inputs, not solver-driven enlargement of contact budgets.
    let sides = member(contact(&mut root), "nonmatching");
    put(sides, "max_pair_tests", num(1_000_000.0));
    put(sides, "max_overlap_triangles", num(100_000.0));
    let mut policy = J::parse(r#"{"max_refinements":6,"consecutive_passes":2,"temperature_tolerance_k":1,"max_vertices":20000,"max_tetrahedra":100000}"#).unwrap();
    if adaptive {
        put(&mut policy, "strategy", J::Str("goal-recovery".into()));
        put(&mut policy, "marking_fraction", num(0.5));
    }
    put(&mut root, "mesh_convergence", policy);
    root
}

fn declaration(root: &J) -> &J {
    &root.path(&["solid", "contacts"]).unwrap().as_array().unwrap()[0]
}

#[test]
fn uniform_and_local_radiating_studies_preserve_the_problem_and_replay_their_field() {
    for adaptive in [false, true] {
        let base = input(adaptive);
        let result = run(&base);
        let study = result.get("mesh_convergence").unwrap();
        let history = study.get("history").unwrap().as_array().unwrap();
        assert!(history.len() >= if adaptive { 4 } else { 3 });
        assert_eq!(history.last().unwrap().str_field("arrived_by"), Some("uniform"));
        if adaptive {
            assert_eq!(study.get("global_confirmation"), Some(&J::Bool(true)));
            assert_eq!(history[1].str_field("arrived_by"), Some("marked-edge-stars"));
            assert!(n(study, "total_adjoint_sweeps") > 0.0);
        } else {
            assert_eq!(n(study, "total_adjoint_sweeps"), 0.0);
            assert_eq!(n(&history[1], "tetrahedra"), 336.0);
            assert_eq!(n(&history[2], "tetrahedra"), 2688.0);
        }
        for row in history { near(n(row, "source_w"), 1.0, 1e-7); }
        let resolved = study.get("resolved_request").unwrap();
        assert!(resolved.get("mesh_convergence").is_none());
        assert!(resolved.path(&["solid", "component_power"]).is_none());
        assert!(resolved.path(&["solid", "nodal_source_w_m3"]).is_some());
        assert_eq!(resolved.path(&["solid", "materials"]), base.path(&["solid", "materials"]));
        assert_eq!(resolved.get("radiation"), base.get("radiation"));
        let before = declaration(&base);
        let after = declaration(resolved);
        assert!(after.get("face_pairs").is_none());
        for key in ["name", "source", "resistance_m2_k_w", "side_a_material", "side_b_material"] {
            assert_eq!(before.get(key), after.get(key));
        }
        for key in ["plane_tolerance_m", "coverage_relative_tolerance", "max_pair_tests", "max_overlap_triangles"] {
            assert_eq!(before.get("nonmatching").unwrap().get(key), after.get("nonmatching").unwrap().get(key));
        }
        let flux = &result.get("contacts").unwrap().as_array().unwrap()[0];
        near(n(flux, "area_m2"), 0.01, 1e-12);
        near(n(flux, "conductance_w_k"), 1.0, 1e-10);
        assert_eq!(flux.str_field("discretization"), Some("planar-common-refinement-P1"));
        let rad = result.get("radiation").unwrap();
        near(n(rad, "radiative_out_w") + n(rad, "convective_out_w"), 1.0, 1e-7);
        assert!(n(rad, "radiative_out_w") > 0.001);
        assert_eq!(rad.get("adjoint"), Some(&J::Null));
        assert_eq!(result.get("contact_sensitivities"), Some(&J::Null));
        let replay = run(resolved);
        for key in ["solid_temperatures_k", "objective", "contacts", "walls"] {
            assert_eq!(result.get(key), replay.get(key), "refined replay changed {key}");
        }
        let mut dry = resolved.clone();
        remove(&mut dry, "radiation");
        assert!(peak(&run(&dry)) > peak(&result) + 0.001);
    }
}

#[test]
fn independent_contact_face_order_does_not_change_adaptation() {
    let base = input(true);
    let original = run(&base);
    let mut reversed = base;
    for key in ["side_a_faces", "side_b_faces"] {
        rows(member(member(contact(&mut reversed), "nonmatching"), key)).reverse();
    }
    let reordered = run(&reversed);
    assert_eq!(original.get("solid_temperatures_k"), reordered.get("solid_temperatures_k"));
    assert_eq!(original.path(&["mesh_convergence", "history"]), reordered.path(&["mesh_convergence", "history"]));
}

#[test]
fn nonmatching_adaptation_cannot_skip_global_confirmation_or_raise_geometry_budgets() {
    let mut too_short = input(true);
    put(member(&mut too_short, "mesh_convergence"), "max_refinements", num(2.0));
    put(member(&mut too_short, "mesh_convergence"), "temperature_tolerance_k", num(100.0));
    let failed = output(&too_short);
    assert_eq!(failed.status.code(), Some(6));
    assert!(failed.stdout.is_empty());
    assert!(String::from_utf8_lossy(&failed.stderr).contains("global confirmation"));
    let mut pair_cap = input(false);
    put(member(contact(&mut pair_cap), "nonmatching"), "max_pair_tests", num(10_000.0));
    let failed = output(&pair_cap);
    assert_eq!(failed.status.code(), Some(6));
    assert!(failed.stdout.is_empty());
    let message = String::from_utf8_lossy(&failed.stderr);
    assert!(message.contains("24976") && message.contains("max_pair_tests=10000"));
}

#[test]
fn overlap_and_marker_failures_remain_refusals_not_partial_mesh_results() {
    let base = input(true);
    let mut overlap = base.clone();
    put(member(contact(&mut overlap), "nonmatching"), "max_overlap_triangles", num(1.0));
    let mut dual = base.clone();
    put(member(&mut dual, "budgets"), "derivative_iterations", num(1.0));
    let mut wall = base;
    put(member(&mut wall, "budgets"), "wall_seconds", num(1e-12));
    for failed in [overlap, dual, wall] {
        let result = output(&failed);
        assert!(!result.status.success());
        assert!(result.stdout.is_empty());
    }
}
