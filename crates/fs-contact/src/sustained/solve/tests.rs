use super::*;

fn controls() -> SolveControls {
    SolveControls {
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1.0e-11,
            relative_tolerance: 0.0,
            ..NewtonKrylovConfig::default()
        },
        difference_step_n: 1.0e-3,
        force_residual_tolerance_n: 1.0e-9,
        ..SolveControls::default()
    }
}

fn scales() -> NaturalMapScales {
    NaturalMapScales {
        normal_n_per_m: 25.0,
        tangent_n_per_m: 25.0,
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}, tolerance {tolerance:.3e}"
    );
}

/// Independent constant-force midpoint mechanics, not a replacement body owner.
/// The 2 kg body begins on a plane under 20 N gravity. Each response starts at
/// the same accepted state; no endpoint or solver force is used as a new start.
fn particle_path(force: [f64; 3], tangent_load: f64, tangent_speed: f64) -> [f64; 3] {
    let h = 0.4;
    let half_mobility = h * h / (2.0 * 2.0);
    [
        half_mobility * (force[0] - 20.0),
        h * tangent_speed + half_mobility * (tangent_load + force[1]),
        half_mobility * force[2],
    ]
}

#[test]
fn cold_loaded_rest_slide_and_breakaway_recover_force_impulse_and_work() {
    // (authored tangent N, initial m/s, reaction N, endpoint m/s, slip m, loss J)
    for (load, before, force, after, slip, loss) in [
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        (4.0, 0.0, -4.0, 0.0, 0.0, 0.0),
        (0.0, 3.0, -5.0, 2.0, 1.0, 5.0),
        (6.0, 0.0, -5.0, 0.2, 0.04, 0.2),
    ] {
        let mut calls = Vec::new();
        let mut response = |forces: &[[f64; 3]], stop: &mut dyn FnMut() -> bool| {
            assert!(!stop());
            calls.push(forces.to_vec());
            Ok::<_, &'static str>(vec![particle_path(forces[0], load, before)])
        };
        let report = solve(
            &[0.25],
            &[scales()],
            &[[0.0; 3]],
            controls(),
            &mut response,
            &mut || false,
        )
        .expect("constant-force mechanics must solve from zero reactions");
        assert!(report.nonlinear.converged);
        assert!(report.maximum_residual_n <= 1.0e-9);
        close(report.forces_n[0][0], 20.0, 1.0e-8);
        close(report.forces_n[0][1], force, 1.0e-8);
        close(report.forces_n[0][2], 0.0, 1.0e-8);
        close(report.paths_m[0][0], 0.0, 1.0e-9);
        close(report.paths_m[0][1], slip, 1.0e-9);
        let endpoint = before + 0.4 * (load + report.forces_n[0][1]) / 2.0;
        close(endpoint, after, 1.0e-9);
        close(2.0 * (endpoint - before), 0.4 * (load + force), 1.0e-9);
        close(-report.forces_n[0][1] * report.paths_m[0][1], loss, 1.0e-8);
        let kinetic_change = endpoint * endpoint - before * before;
        close(kinetic_change, load * report.paths_m[0][1] - loss, 1.0e-8);
        assert_eq!(report.response_evaluations, calls.len());
        assert_eq!(calls.last(), Some(&report.forces_n));
        assert!(calls.len() >= 3);
    }
}

#[test]
fn off_diagonal_inertia_couples_normal_reaction_and_sliding_friction() {
    // M = [[2,1],[1,2]] kg. Its inverse couples normal and tangent reactions.
    // Independent force balance gives fn=160/9, ft=-40/9, not 20 and -5 N.
    let mut response = |forces: &[[f64; 3]], _: &mut dyn FnMut() -> bool| {
        let f = forces[0];
        let an = (2.0 * (f[0] - 20.0) - f[1]) / 3.0;
        let at = (-(f[0] - 20.0) + 2.0 * f[1]) / 3.0;
        Ok::<_, &'static str>(vec![[0.08 * an, 1.2 + 0.08 * at, 0.04 * f[2]]])
    };
    let report = solve(
        &[0.25],
        &[scales()],
        &[[0.0; 3]],
        controls(),
        &mut response,
        &mut || false,
    )
    .unwrap();
    close(report.forces_n[0][0], 160.0 / 9.0, 1.0e-8);
    close(report.forces_n[0][1], -40.0 / 9.0, 1.0e-8);
    close(report.paths_m[0][0], 0.0, 1.0e-9);
    close(report.paths_m[0][1], 46.0 / 45.0, 1.0e-9);
    close(
        -report.forces_n[0][1] * report.paths_m[0][1],
        368.0 / 81.0,
        1.0e-8,
    );
}

#[test]
fn simultaneous_contacts_share_the_same_mechanical_response() {
    // Three 1 kg particles with B=[[-1,1,0],[0,-1,1]], authored normal
    // forces [6,-2,-4] N and h=.2 s. Both gaps require the coupled solution
    // [6,4] N. Solving contacts independently cannot balance the middle body.
    let mut response = |forces: &[[f64; 3]], _: &mut dyn FnMut() -> bool| {
        let acceleration = [
            6.0 - forces[0][0],
            -2.0 + forces[0][0] - forces[1][0],
            -4.0 + forces[1][0],
        ];
        Ok::<_, &'static str>(vec![
            [0.02 * (acceleration[1] - acceleration[0]), 0.0, 0.0],
            [0.02 * (acceleration[2] - acceleration[1]), 0.0, 0.0],
        ])
    };
    let report = solve(
        &[0.0, 0.0],
        &[scales(); 2],
        &[[0.0; 3]; 2],
        controls(),
        &mut response,
        &mut || false,
    )
    .unwrap();
    close(report.forces_n[0][0], 6.0, 1.0e-8);
    close(report.forces_n[1][0], 4.0, 1.0e-8);
    for path in &report.paths_m {
        close(path[0], 0.0, 1.0e-9);
    }
    let body_contact_impulses = [
        -0.2 * report.forces_n[0][0],
        0.2 * (report.forces_n[0][0] - report.forces_n[1][0]),
        0.2 * report.forces_n[1][0],
    ];
    close(body_contact_impulses.iter().sum(), 0.0, 1.0e-12);
}

#[test]
fn compliant_midpoint_response_changes_slip_and_preserves_storage_balance() {
    // A 2 kg tangent mass with a 100 N/m spring begins at q=0, v=3 m/s.
    // Under held Coulomb ft=-5 N, midpoint q1=(2 m v0/h + ft)/(2 m/h²+k/2).
    let mut response = |forces: &[[f64; 3]], _: &mut dyn FnMut() -> bool| {
        let f = forces[0];
        Ok::<_, &'static str>(vec![[
            0.0025 * (f[0] - 20.0),
            (120.0 + f[1]) / 450.0,
            f[2] / 400.0,
        ]])
    };
    let report = solve(
        &[0.25],
        &[scales()],
        &[[0.0; 3]],
        controls(),
        &mut response,
        &mut || false,
    )
    .unwrap();
    close(report.forces_n[0][0], 20.0, 1.0e-8);
    close(report.forces_n[0][1], -5.0, 1.0e-8);
    let q = report.paths_m[0][1];
    let v = 2.0 * q / 0.1 - 3.0;
    close(q, 23.0 / 90.0, 1.0e-9);
    close(v, 19.0 / 9.0, 1.0e-8);
    close(v * v + 50.0 * q * q - 9.0, -23.0 / 18.0, 1.0e-8);
    close(-report.forces_n[0][1] * q, 23.0 / 18.0, 1.0e-8);
}

#[test]
fn source_refusal_unwind_and_cancellation_stop_further_mechanical_calls() {
    #[derive(Debug, PartialEq)]
    struct SourceFailure(u64);
    let mut calls = 0;
    let error = solve(
        &[0.25],
        &[scales()],
        &[[0.0; 3]],
        controls(),
        &mut |_, _| {
            calls += 1;
            Err::<Vec<[f64; 3]>, _>(SourceFailure(71))
        },
        &mut || false,
    )
    .unwrap_err();
    assert!(matches!(error, SolveError::Response(SourceFailure(71))));
    assert_eq!(calls, 1);

    let mut calls = 0;
    let error = solve(
        &[0.25],
        &[scales()],
        &[[0.0; 3]],
        controls(),
        &mut |_, _| {
            calls += 1;
            panic!("original response refusal");
            #[allow(unreachable_code)]
            Ok::<Vec<[f64; 3]>, SourceFailure>(vec![])
        },
        &mut || false,
    )
    .unwrap_err();
    assert!(matches!(error, SolveError::ResponsePanicked));
    assert_eq!(calls, 1);

    let mut calls = 0;
    let error = solve(
        &[0.25],
        &[scales()],
        &[[0.0; 3]],
        controls(),
        &mut |_, _| {
            calls += 1;
            Ok::<_, SourceFailure>(vec![[0.0; 3]])
        },
        &mut || true,
    )
    .unwrap_err();
    assert!(matches!(error, SolveError::Cancelled));
    assert_eq!(calls, 0);
}

#[test]
fn source_shape_nonfinite_and_callback_cancellation_are_not_newton_success() {
    let error = solve(
        &[0.25],
        &[scales()],
        &[[0.0; 3]],
        controls(),
        &mut |_, _| Ok::<_, &'static str>(vec![]),
        &mut || false,
    )
    .unwrap_err();
    assert!(matches!(
        error,
        SolveError::ResponseShape {
            expected: 1,
            actual: 0
        }
    ));
    let error = solve(
        &[0.25],
        &[scales()],
        &[[0.0; 3]],
        controls(),
        &mut |_, _| Ok::<_, &'static str>(vec![[f64::NAN, 0.0, 0.0]]),
        &mut || false,
    )
    .unwrap_err();
    assert!(matches!(error, SolveError::NonFinite(_)));

    let cancel = Cell::new(false);
    let mut calls = 0;
    let error = solve(
        &[0.25],
        &[scales()],
        &[[0.0; 3]],
        controls(),
        &mut |_, stop| {
            calls += 1;
            cancel.set(true);
            assert!(stop());
            Ok::<_, &'static str>(vec![[0.0; 3]])
        },
        &mut || cancel.get(),
    )
    .unwrap_err();
    assert!(matches!(error, SolveError::Cancelled));
    assert_eq!(calls, 1);
}

#[test]
fn final_fresh_response_is_budgeted_and_must_independently_pass_the_law() {
    let budget = SolveControls {
        maximum_response_evaluations: 1,
        ..controls()
    };
    let mut calls = 0;
    let error = solve(
        &[0.25],
        &[scales()],
        &[[20.0, 0.0, 0.0]],
        budget,
        &mut |f, _| {
            calls += 1;
            Ok::<_, &'static str>(vec![particle_path(f[0], 0.0, 0.0)])
        },
        &mut || false,
    )
    .unwrap_err();
    assert!(matches!(
        error,
        SolveError::ResponseBudget {
            evaluations: 1,
            maximum: 1
        }
    ));
    assert_eq!(calls, 1);

    // Even an initially converged force is freshly replayed. A body source
    // changing its result cannot leave the cached zero residual authoritative.
    let mut calls = 0;
    let error = solve(
        &[0.25],
        &[scales()],
        &[[20.0, 0.0, 0.0]],
        controls(),
        &mut |f, _| {
            calls += 1;
            let mut path = particle_path(f[0], 0.0, 0.0);
            if calls == 2 {
                path[0] = 0.1;
            }
            Ok::<_, &'static str>(vec![path])
        },
        &mut || false,
    )
    .unwrap_err();
    assert!(matches!(error, SolveError::Residual { .. }));
    assert_eq!(calls, 2);
}

#[test]
fn admission_refuses_unbounded_workspace_before_any_mechanical_evaluation() {
    let mut calls = 0;
    let tiny = SolveControls {
        maximum_workspace_entries: 1,
        ..controls()
    };
    let error = solve(
        &[0.25],
        &[scales()],
        &[[0.0; 3]],
        tiny,
        &mut |_, _| {
            calls += 1;
            Ok::<_, &'static str>(vec![[0.0; 3]])
        },
        &mut || false,
    )
    .unwrap_err();
    assert!(matches!(error, SolveError::WorkspaceBudget { .. }));
    assert_eq!(calls, 0);
}

#[test]
fn strict_normal_feasibility_reevaluates_force_instead_of_clamping_the_gap() {
    let initial = 20.0_f64.next_down();
    let mut calls = Vec::new();
    let report = solve(
        &[0.0],
        &[scales()],
        &[[initial, 0.0, 0.0]],
        controls(),
        &mut |forces, _| {
            calls.push(forces.to_vec());
            Ok::<_, &'static str>(vec![particle_path(forces[0], 0.0, 0.0)])
        },
        &mut || false,
    )
    .unwrap();
    assert!(report.feasibility_refinements > 0);
    assert!(report.paths_m[0][0] >= 0.0);
    assert!(report.forces_n[0][0] >= 20.0);
    assert_eq!(
        report.paths_m[0],
        particle_path(report.forces_n[0], 0.0, 0.0)
    );
    assert_eq!(calls.last(), Some(&report.forces_n));
    assert_eq!(report.response_evaluations, calls.len());
    assert!(report.maximum_residual_n <= controls().force_residual_tolerance_n);
}

#[test]
fn feasibility_polishing_keeps_the_original_response_budget_and_source_error() {
    let initial = [[20.0_f64.next_down(), 0.0, 0.0]];
    let mut limited = controls();
    // Initial Newton residual and first final replay consume both calls.
    limited.maximum_response_evaluations = 2;
    let mut calls = 0;
    let failed = solve(
        &[0.0],
        &[scales()],
        &initial,
        limited,
        &mut |forces, _| {
            calls += 1;
            Ok::<_, &'static str>(vec![particle_path(forces[0], 0.0, 0.0)])
        },
        &mut || false,
    );
    assert!(matches!(
        failed,
        Err(SolveError::ResponseBudget {
            evaluations: 2,
            maximum: 2
        })
    ));
    assert_eq!(calls, 2);

    let mut calls = 0;
    let failed = solve(
        &[0.0],
        &[scales()],
        &initial,
        controls(),
        &mut |forces, _| {
            calls += 1;
            if calls == 3 {
                return Err("original force correction refused");
            }
            Ok(vec![particle_path(forces[0], 0.0, 0.0)])
        },
        &mut || false,
    );
    assert!(matches!(
        failed,
        Err(SolveError::Response("original force correction refused"))
    ));
    assert_eq!(calls, 3);
}

#[test]
fn an_unresponsive_penetrating_source_is_never_turned_into_an_admissible_gap() {
    let mut calls = 0;
    let failed = solve(
        &[0.0],
        &[scales()],
        &[[1.0, 0.0, 0.0]],
        controls(),
        &mut |_, _| {
            calls += 1;
            Ok::<_, &'static str>(vec![[-1.0e-15, 0.0, 0.0]])
        },
        &mut || false,
    );
    assert!(matches!(
        failed,
        Err(SolveError::Feasibility {
            gap_m,
            refinements: 0,
            ..
        }) if gap_m == -1.0e-15
    ));
    // An actual force probe demonstrates that this source cannot open.
    assert_eq!(calls, 3);
}

#[test]
fn near_zero_tensile_trial_is_replaced_only_by_a_fresh_zero_force_response() {
    let mut last = None;
    let report = solve(
        &[0.0],
        &[scales()],
        &[[-1.0e-12, 0.0, 0.0]],
        controls(),
        &mut |forces, _| {
            last = Some(forces[0]);
            Ok::<_, &'static str>(vec![[1.0, 0.0, 0.0]])
        },
        &mut || false,
    )
    .unwrap();
    assert_eq!(report.forces_n[0], [0.0; 3]);
    assert_eq!(report.paths_m[0], [1.0, 0.0, 0.0]);
    assert_eq!(last, Some(report.forces_n[0]));
    assert_eq!(report.feasibility_refinements, 1);
}
