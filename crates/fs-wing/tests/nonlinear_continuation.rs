//! Independent equation, force, continuation-domain and work-budget checks.
//! Repro: cargo test -p fs-wing --test nonlinear_continuation

use fs_wing::nonlinear::continuation::{
    ContinuationControls, MAX_JACOBIAN_PRODUCTS, MAX_NEWTON_ITERATIONS, MAX_RESIDUAL_EVALUATIONS,
    solve_with_incident,
};
use fs_wing::nonlinear::{
    GAMMA_TOL, InfluenceOperator, StripRegime, StripSpec, solve_nonlinear_with_incident,
};
use fs_wing::{Panel, SurfaceId, flat_surface, induced_velocity_free};
use std::cell::Cell;

const RHO: f64 = 1.225;
const FREE: [f64; 3] = [10.0, 0.4, 0.5];

fn controls() -> ContinuationControls {
    ContinuationControls {
        maximum_incidence_rad: 0.6,
        maximum_incidence_change_rad: 0.08,
        maximum_row_condition: 8.0,
        minimum_trailing_clearance_chords: 0.02,
        maximum_newton_iterations: 32,
        maximum_residual_evaluations: 600,
        maximum_jacobian_products: 4096,
        section_difference_rad: 1e-6,
    }
}

fn layout() -> (Vec<Panel>, Vec<StripSpec>) {
    let panels = flat_surface(SurfaceId::WingLower, 6.0, 1.3, 0.0, 0.0, 4, 2).unwrap();
    let strips = (0..4)
        .map(|s| StripSpec {
            panel_indices: vec![s, s + 4],
            chord_m: 1.3,
            twist_rad: 0.005 * (s as f64 - 1.5),
        })
        .collect();
    (panels, strips)
}

fn section(s: usize, alpha: f64) -> (f64, StripRegime) {
    (
        core::f64::consts::TAU * (alpha + 0.015 * (s as f64 + 1.0)),
        StripRegime::Attached,
    )
}

fn subtract(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    core::array::from_fn(|i| a[i] - b[i])
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.into_iter().zip(b).map(|(x, y)| x * y).sum()
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Public complete horseshoes minus an independent finite-bound formula.
/// This does not reuse the nonlinear owner's cached trailing coefficients.
fn trailing_z(point: [f64; 3], panels: &[Panel], gamma: &[f64]) -> f64 {
    let speed = dot(FREE, FREE).sqrt();
    let stream = FREE.map(|x| x / speed);
    let mut z = induced_velocity_free(point, panels, gamma, stream)[2];
    for (panel, g) in panels.iter().zip(gamma) {
        let r1 = subtract(point, panel.a);
        let r2 = subtract(point, panel.b);
        let r0 = subtract(panel.b, panel.a);
        let product = cross(r1, r2);
        let c2 = dot(product, product);
        let n1 = dot(r1, r1).sqrt();
        let n2 = dot(r2, r2).sqrt();
        if c2 >= 1e-8 && n1 >= 1e-8 && n2 >= 1e-8 {
            z -= g * product[2] * (dot(r0, r1) / n1 - dot(r0, r2) / n2)
                / (4.0 * core::f64::consts::PI * c2);
        }
    }
    z
}

/// A separate small dense elimination oracle, not the production Krylov owner.
fn eliminate(mut matrix: Vec<Vec<f64>>, mut rhs: Vec<f64>) -> Vec<f64> {
    let n = rhs.len();
    for k in 0..n {
        let pivot = (k..n)
            .max_by(|&a, &b| matrix[a][k].abs().total_cmp(&matrix[b][k].abs()))
            .unwrap();
        matrix.swap(k, pivot);
        rhs.swap(k, pivot);
        assert!(matrix[k][k].abs() > 1e-8);
        for i in k + 1..n {
            let factor = matrix[i][k] / matrix[k][k];
            for j in k + 1..n {
                matrix[i][j] -= factor * matrix[k][j];
            }
            rhs[i] -= factor * rhs[k];
        }
    }
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let upper: f64 = (i + 1..n).map(|j| matrix[i][j] * x[j]).sum();
        x[i] = (rhs[i] - upper) / matrix[i][i];
    }
    x
}

#[test]
fn nonuniform_section_solution_matches_an_independent_dense_system_and_actual_force() {
    let (panels, strips) = layout();
    let op = InfluenceOperator::build(&panels, FREE, RHO).unwrap();
    let incident: Vec<_> = (0..4)
        .map(|s| [9.5 + 0.3 * s as f64, 0.2 * s as f64, 0.25 + 0.12 * s as f64])
        .collect();
    let wash = [0.3, 0.7, 0.4, 0.9];
    let local: Vec<_> = incident
        .iter()
        .zip(wash)
        .map(|(v, du)| [v[0] + du, v[1], v[2]])
        .collect();
    let solved = solve_with_incident(
        &op,
        &panels,
        &strips,
        FREE,
        RHO,
        &section,
        &incident,
        Some(&wash),
        None,
        controls(),
    )
    .unwrap();
    let mut matrix = vec![vec![0.0; 4]; 4];
    let mut rhs = vec![0.0; 4];
    let mut ratios = vec![vec![]; 4];
    for (t, strip) in strips.iter().enumerate() {
        let total: f64 = strip
            .panel_indices
            .iter()
            .map(|&j| op.linear().gamma[j])
            .sum();
        ratios[t] = strip
            .panel_indices
            .iter()
            .map(|&j| op.linear().gamma[j] / total)
            .collect();
    }
    for (s, strip) in strips.iter().enumerate() {
        let front = panels[strip.panel_indices[0]];
        let point = core::array::from_fn(|i| 0.5 * (front.a[i] + front.b[i]));
        let velocity = local[s][0].hypot(local[s][2]);
        let alpha = fs_math::det::atan2(local[s][2], local[s][0]) + strip.twist_rad;
        rhs[s] = -0.5 * velocity * strip.chord_m * section(s, alpha).0;
        for (t, other) in strips.iter().enumerate() {
            let mut unit = vec![0.0; panels.len()];
            for (k, &j) in other.panel_indices.iter().enumerate() {
                unit[j] = ratios[t][k];
            }
            matrix[s][t] = f64::from(s == t)
                + 0.5 * strip.chord_m * core::f64::consts::TAU * trailing_z(point, &panels, &unit);
        }
    }
    let exact = eliminate(matrix, rhs);
    let mut force = 0.0;
    for (s, strip) in strips.iter().enumerate() {
        for (k, &j) in strip.panel_indices.iter().enumerate() {
            let gamma = exact[s] * ratios[s][k];
            assert!((solved.solution.gamma[j] - gamma).abs() < 2e-10);
            let span = subtract(panels[j].b, panels[j].a);
            force -= RHO * gamma * cross(local[s], span)[2];
        }
    }
    assert!((force - solved.solution.total_lift_n).abs() < 5e-9);
    assert!(solved.solution.residual < GAMMA_TOL);
    assert_eq!(solved.solution.gamma, solved.continuation.gamma());
    let direct = solve_with_incident(
        &op,
        &panels,
        &strips,
        FREE,
        RHO,
        &section,
        &local,
        None,
        None,
        controls(),
    )
    .unwrap();
    assert_eq!(
        solved, direct,
        "wash enters the exact same section velocity"
    );
    let warm = solve_with_incident(
        &op,
        &panels,
        &strips,
        FREE,
        RHO,
        &section,
        &local,
        None,
        Some(&solved.continuation),
        controls(),
    )
    .unwrap();
    assert_eq!(warm.solution.iterations, 0);
    for (a, b) in warm.solution.gamma.iter().zip(&solved.solution.gamma) {
        assert!((a - b).abs() < 2e-14);
    }
}

#[test]
fn smooth_attached_to_blended_transition_continues_without_rejecting_regime_labels() {
    let panels = flat_surface(SurfaceId::WingLower, 8.0, 1.0, 0.0, 0.0, 1, 2).unwrap();
    let mut strips = vec![StripSpec {
        panel_indices: vec![0, 1],
        chord_m: 1.0,
        twist_rad: 0.0,
    }];
    let op = InfluenceOperator::build(&panels, FREE, RHO).unwrap();
    let closure = |_: usize, alpha: f64| {
        let attached = core::f64::consts::TAU * alpha;
        if alpha.abs() <= 0.30 {
            (attached, StripRegime::Attached)
        } else {
            let t = ((alpha.abs() - 0.30) / 0.15).min(1.0);
            let blend = t * t * (3.0 - 2.0 * t);
            let separated = 1.98 * alpha.sin() * alpha.cos();
            (
                attached * (1.0 - blend) + separated * blend,
                if t < 1.0 {
                    StripRegime::Blended
                } else {
                    StripRegime::Separated
                },
            )
        }
    };
    let mut accepted = solve_with_incident(
        &op,
        &panels,
        &strips,
        FREE,
        RHO,
        &closure,
        &[FREE],
        None,
        None,
        controls(),
    )
    .unwrap();
    let initial_branch = accepted.solution.branch_id.clone();
    let mut crossed = false;
    for tick in 1..=19 {
        strips[0].twist_rad = 0.02 * tick as f64;
        let next = solve_with_incident(
            &op,
            &panels,
            &strips,
            FREE,
            RHO,
            &closure,
            &[FREE],
            None,
            Some(&accepted.continuation),
            controls(),
        )
        .unwrap();
        assert!(next.continuation.maximum_incidence_change_rad() < 0.04);
        assert!(next.solution.residual < GAMMA_TOL);
        if next.solution.regimes[0] == StripRegime::Blended {
            crossed = true;
            assert_ne!(next.solution.branch_id, initial_branch);
        }
        accepted = next;
    }
    assert!(
        crossed,
        "the real smooth blend boundary must actually be crossed"
    );
}

#[test]
fn geometry_cancellation_filament_and_incidence_domain_refusals_preserve_the_anchor() {
    let (panels, strips) = layout();
    let op = InfluenceOperator::build(&panels, FREE, RHO).unwrap();
    let incident = vec![FREE; strips.len()];
    let accepted = solve_with_incident(
        &op,
        &panels,
        &strips,
        FREE,
        RHO,
        &section,
        &incident,
        None,
        None,
        controls(),
    )
    .unwrap();
    let saved = accepted.continuation.clone();
    let mut changed = panels.clone();
    changed[0].ctrl[2] += 0.01;
    let changed_op = InfluenceOperator::build(&changed, FREE, RHO).unwrap();
    let error = solve_with_incident(
        &changed_op,
        &changed,
        &strips,
        FREE,
        RHO,
        &section,
        &incident,
        None,
        Some(&accepted.continuation),
        controls(),
    )
    .unwrap_err();
    assert_eq!(error.code, "nonlinear-continuation-stale");
    let mut changed_ownership = strips.clone();
    changed_ownership[0].panel_indices.reverse();
    assert_eq!(
        solve_with_incident(
            &op,
            &panels,
            &changed_ownership,
            FREE,
            RHO,
            &section,
            &incident,
            None,
            Some(&accepted.continuation),
            controls(),
        )
        .unwrap_err()
        .code,
        "nonlinear-continuation-stale"
    );
    let restrictive = ContinuationControls {
        minimum_trailing_clearance_chords: 2.0,
        ..controls()
    };
    assert_eq!(
        solve_with_incident(
            &op,
            &panels,
            &strips,
            FREE,
            RHO,
            &section,
            &incident,
            None,
            Some(&accepted.continuation),
            restrictive,
        )
        .unwrap_err()
        .code,
        "nonlinear-filament-proximity"
    );
    let mut jumped = incident.clone();
    jumped[0][2] += 3.0;
    assert_eq!(
        solve_with_incident(
            &op,
            &panels,
            &strips,
            FREE,
            RHO,
            &section,
            &jumped,
            None,
            Some(&accepted.continuation),
            controls(),
        )
        .unwrap_err()
        .code,
        "nonlinear-continuation-outside-domain"
    );
    for vertical in [0.0, 1e-18] {
        let freestream = [10.0, 0.0, vertical];
        let singular = InfluenceOperator::build(&panels, freestream, RHO).unwrap();
        assert_eq!(
            solve_with_incident(
                &singular,
                &panels,
                &strips,
                freestream,
                RHO,
                &section,
                &vec![freestream; strips.len()],
                None,
                None,
                controls(),
            )
            .unwrap_err()
            .code,
            "nonlinear-row-split-ill-conditioned"
        );
    }
    assert_eq!(accepted.continuation, saved);
}

#[test]
fn solver_work_and_section_source_failures_are_bounded_and_never_accepted() {
    let (panels, strips) = layout();
    let op = InfluenceOperator::build(&panels, FREE, RHO).unwrap();
    let incident = vec![FREE; strips.len()];
    for limited in [
        ContinuationControls {
            maximum_residual_evaluations: 1,
            ..controls()
        },
        ContinuationControls {
            maximum_jacobian_products: 1,
            ..controls()
        },
    ] {
        let calls = Cell::new(0);
        let source = |s, alpha| {
            calls.set(calls.get() + 1);
            section(s, alpha)
        };
        let error = solve_with_incident(
            &op, &panels, &strips, FREE, RHO, &source, &incident, None, None, limited,
        )
        .unwrap_err();
        assert_eq!(error.code, "nonlinear-continuation-budget");
        assert!(
            calls.get()
                <= strips.len()
                    * (limited.maximum_residual_evaluations
                        + 3 * limited.maximum_jacobian_products)
        );
    }
    assert_eq!(
        solve_with_incident(
            &op,
            &panels,
            &strips,
            FREE,
            RHO,
            &|_, _| (f64::NAN, StripRegime::Attached),
            &incident,
            None,
            None,
            controls(),
        )
        .unwrap_err()
        .code,
        "nonlinear-did-not-converge"
    );
    let tiny = ContinuationControls {
        section_difference_rad: f64::MIN_POSITIVE,
        ..controls()
    };
    assert_eq!(
        solve_with_incident(
            &op, &panels, &strips, FREE, RHO, &section, &incident, None, None, tiny,
        )
        .unwrap_err()
        .code,
        "nonlinear-continuation-invalid"
    );
    for rho in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            solve_with_incident(
                &op,
                &panels,
                &strips,
                FREE,
                rho,
                &section,
                &incident,
                None,
                None,
                controls(),
            )
            .unwrap_err()
            .code,
            "nonlinear-continuation-invalid"
        );
    }
    for malformed in [
        vec![],
        vec![0.0; panels.len() + 1],
        vec![f64::NAN; panels.len()],
    ] {
        assert_eq!(
            solve_nonlinear_with_incident(
                &op,
                &panels,
                &strips,
                FREE,
                RHO,
                &section,
                &incident,
                None,
                Some(&malformed),
            )
            .unwrap_err()
            .code,
            "strips-invalid"
        );
    }
    let caps = ContinuationControls {
        maximum_newton_iterations: MAX_NEWTON_ITERATIONS,
        maximum_residual_evaluations: MAX_RESIDUAL_EVALUATIONS,
        maximum_jacobian_products: MAX_JACOBIAN_PRODUCTS,
        ..controls()
    };
    caps.admit().unwrap();
    for invalid in [
        ContinuationControls {
            maximum_newton_iterations: 0,
            ..caps
        },
        ContinuationControls {
            maximum_newton_iterations: MAX_NEWTON_ITERATIONS + 1,
            ..caps
        },
        ContinuationControls {
            maximum_residual_evaluations: 0,
            ..caps
        },
        ContinuationControls {
            maximum_residual_evaluations: MAX_RESIDUAL_EVALUATIONS + 1,
            ..caps
        },
        ContinuationControls {
            maximum_jacobian_products: 0,
            ..caps
        },
        ContinuationControls {
            maximum_jacobian_products: MAX_JACOBIAN_PRODUCTS + 1,
            ..caps
        },
        ContinuationControls {
            maximum_incidence_rad: core::f64::consts::FRAC_PI_2,
            ..caps
        },
        ContinuationControls {
            maximum_incidence_change_rad: 0.0,
            ..caps
        },
        ContinuationControls {
            maximum_row_condition: 0.9,
            ..caps
        },
        ContinuationControls {
            minimum_trailing_clearance_chords: 0.0,
            ..caps
        },
        ContinuationControls {
            section_difference_rad: 0.0,
            ..caps
        },
    ] {
        assert_eq!(
            invalid.admit().unwrap_err().code,
            "nonlinear-continuation-invalid"
        );
    }
}
