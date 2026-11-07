//! Native boundary/Kelvin closure, accepted-geometry coupling and frame tests.

use fs_bem::BemError;
use fs_bem::panel2d::{Airfoil2d, PreparedPanel2d, naca4_symmetric, solve};
use fs_bem::wake2d::{PanelWakeTransportOptions, WakeSim};

fn circulation(wake: &WakeSim) -> f64 {
    let mut sum = wake.history().last().map_or(0.0, |row| row.bound);
    let mut correction = 0.0;
    for vortex in wake.wake() {
        let next = sum + vortex.gamma;
        correction += if sum.abs() >= vortex.gamma.abs() {
            (sum - next) + vortex.gamma
        } else {
            (vortex.gamma - next) + sum
        };
        sum = next;
    }
    sum + correction
}

fn assert_boundary(wake: &WakeSim, panels: &PreparedPanel2d, alpha: f64, core: f64) {
    let particles: Vec<_> = wake
        .wake()
        .iter()
        .map(|vortex| vortex.as_particle())
        .collect();
    let freestream = [fs_math::det::cos(alpha), fs_math::det::sin(alpha)];
    let incident: Vec<_> = panels
        .midpoints()
        .iter()
        .map(|&point| {
            let u = fs_vpm::induced_velocity(&particles, point, core);
            [freestream[0] + u[0], freestream[1] + u[1]]
        })
        .collect();
    // Independent ordinary solve from the entire accepted wake. This does
    // not use the coupled newborn response or its stored midpoint geometry.
    let reference = panels.solve_incident(&incident).unwrap();
    let actual = wake.panel_state().unwrap();
    let tolerance = 5e-11 * (1.0 + reference.bound_circulation.abs());
    assert!((reference.bound_circulation - actual.bound_circulation).abs() < tolerance);
    for (&a, &b) in actual.sources.iter().zip(&reference.sources) {
        assert!((a - b).abs() < tolerance, "source {a} != {b}");
    }
    assert!(actual.max_normal_residual < tolerance);
    assert!(actual.kutta_residual.abs() < tolerance);
    assert!(circulation(wake).abs() < tolerance);
    assert_eq!(
        wake.history().last().unwrap().bound,
        actual.bound_circulation
    );
}

#[test]
fn coupled_wake_preserves_kelvin_and_satisfies_the_native_endpoint_boundary() {
    let foil = naca4_symmetric(0.12, 32).unwrap();
    let panels = PreparedPanel2d::new(&foil).unwrap();
    let (alpha, dt, core) = (0.1, 0.02, 0.03);
    let mut wake = WakeSim::new_panel_coupled(&foil, alpha, dt, core).unwrap();
    assert!(wake.panel_state().is_none());
    let steady = solve(&foil, alpha).unwrap();
    let perimeter: f64 = foil
        .nodes()
        .iter()
        .zip(foil.nodes().iter().cycle().skip(1))
        .map(|(a, b)| (b[0] - a[0]).hypot(b[1] - a[1]))
        .sum();
    assert!((wake.steady_circulation() - steady.gamma * perimeter).abs() < 2e-12);
    for step in 1..=12 {
        let strengths: Vec<_> = wake.wake().iter().map(|v| v.gamma.to_bits()).collect();
        wake.step().unwrap();
        assert_eq!(wake.history().len(), step);
        for (vortex, &bits) in wake.wake().iter().zip(&strengths) {
            assert_eq!(vortex.gamma.to_bits(), bits);
        }
        assert_boundary(&wake, &panels, alpha, core);
    }
    assert!(
        wake.trace_json(1)
            .unwrap()
            .contains("inviscid-panel-kelvin")
    );
}

#[test]
fn prescribed_incidence_reuses_the_body_and_sheds_its_circulation_change() {
    let foil = naca4_symmetric(0.12, 32).unwrap();
    let panels = PreparedPanel2d::new(&foil).unwrap();
    let core = 0.03;
    let mut wake = WakeSim::new_panel_coupled(&foil, 0.04, 0.02, core).unwrap();
    let mut values = Vec::new();
    for alpha in [0.04, 0.08, 0.12, 0.06, 0.0, -0.04] {
        let old: Vec<_> = wake.wake().iter().map(|v| v.gamma.to_bits()).collect();
        wake.step_panel_coupled_at_angle(alpha).unwrap();
        for (vortex, &gamma) in wake.wake().iter().zip(&old) {
            assert_eq!(vortex.gamma.to_bits(), gamma);
        }
        assert_boundary(&wake, &panels, alpha, core);
        values.push(wake.history().last().unwrap().bound);
    }
    assert!(values.iter().any(|x| (x - values[0]).abs() > 1e-4));
    // The ordinary step retains the last accepted incidence.
    wake.step().unwrap();
    assert_boundary(&wake, &panels, -0.04, core);
}

#[test]
fn normalized_panel_wake_is_covariant_under_rigid_geometry_and_unit_changes() {
    let foil = naca4_symmetric(0.12, 24).unwrap();
    let (alpha, dt, core) = (0.08, 0.025, 0.025);
    let mut reference = WakeSim::new_panel_coupled(&foil, alpha, dt, core).unwrap();
    for _ in 0..4 {
        reference.step().unwrap();
    }
    for scale in [0.01, 1.0, 100.0] {
        let angle = 0.37_f64;
        let (s, c) = angle.sin_cos();
        let nodes = foil
            .nodes()
            .iter()
            .map(|p| {
                [
                    scale * (c * p[0] - s * p[1] + 2.0),
                    scale * (s * p[0] + c * p[1] - 1.5),
                ]
            })
            .collect();
        let changed = Airfoil2d::new(nodes).unwrap();
        let mut wake = WakeSim::new_panel_coupled(&changed, alpha + angle, dt, core).unwrap();
        for _ in 0..4 {
            wake.step().unwrap();
        }
        assert_eq!(reference.wake().len(), wake.wake().len());
        for (a, b) in reference.wake().iter().zip(wake.wake()) {
            assert!((a.gamma - b.gamma).abs() < 3e-11);
            assert!((a.pos[0] - b.pos[0]).hypot(a.pos[1] - b.pos[1]) < 3e-11);
        }
        assert!(
            (reference.panel_state().unwrap().bound_circulation
                - wake.panel_state().unwrap().bound_circulation)
                .abs()
                < 3e-11
        );
    }
}

#[test]
fn refused_incidence_preserves_accepted_wake_panel_state_and_history() {
    let foil = naca4_symmetric(0.12, 24).unwrap();
    let mut wake = WakeSim::new_panel_coupled(&foil, 0.08, 0.025, 0.03).unwrap();
    wake.step().unwrap();
    let before = wake.clone();
    for alpha in [f64::NAN, f64::INFINITY, std::f64::consts::PI] {
        assert!(wake.step_panel_coupled_at_angle(alpha).is_err());
        assert_eq!(wake.trace_json(1).unwrap(), before.trace_json(1).unwrap());
        assert_eq!(
            wake.panel_state().unwrap().sources,
            before.panel_state().unwrap().sources
        );
        for (a, b) in wake.wake().iter().zip(before.wake()) {
            assert_eq!(a.pos, b.pos);
            assert_eq!(a.gamma, b.gamma);
        }
    }
    let mut legacy = WakeSim::new(&foil, 0.08, 0.025, 0.03).unwrap();
    assert!(legacy.step_panel_coupled_at_angle(0.1).is_err());
    assert!(legacy.history().is_empty());
    legacy.step().unwrap();
    assert!(legacy.panel_state().is_none());
    assert!(legacy.trace_json(1).unwrap().contains("inviscid-screening"));
}

fn transport(substeps: usize) -> PanelWakeTransportOptions {
    PanelWakeTransportOptions {
        substeps,
        max_coupled_solves: 2 * substeps + 1,
        max_displacement_fraction: None,
    }
}

fn assert_same_native_state(a: &WakeSim, b: &WakeSim) {
    assert_eq!(a.trace_json(1).unwrap(), b.trace_json(1).unwrap());
    assert_eq!(
        a.steady_circulation().to_bits(),
        b.steady_circulation().to_bits()
    );
    assert_eq!(a.wake().len(), b.wake().len());
    for (a, b) in a.wake().iter().zip(b.wake()) {
        assert_eq!(a.pos.map(f64::to_bits), b.pos.map(f64::to_bits));
        assert_eq!(a.gamma.to_bits(), b.gamma.to_bits());
    }
    let a = a.panel_state().unwrap();
    let b = b.panel_state().unwrap();
    assert_eq!(a.sources, b.sources);
    assert_eq!(a.tangential_velocity, b.tangential_velocity);
    assert_eq!(a.bound_circulation.to_bits(), b.bound_circulation.to_bits());
    assert_eq!(
        a.max_normal_residual.to_bits(),
        b.max_normal_residual.to_bits()
    );
    assert_eq!(a.kutta_residual.to_bits(), b.kutta_residual.to_bits());
}

#[test]
fn transport_refinement_keeps_the_release_schedule_and_converges_at_second_order() {
    let foil = naca4_symmetric(0.12, 24).unwrap();
    let panels = PreparedPanel2d::new(&foil).unwrap();
    let (alpha, dt, core) = (0.08, 0.02, 0.04);
    let mut start = WakeSim::new_panel_coupled(&foil, alpha, dt, core).unwrap();
    for _ in 0..3 {
        start.step().unwrap();
    }
    let previous_strengths: Vec<_> = start.wake().iter().map(|v| v.gamma.to_bits()).collect();
    let mut reference = start.clone();
    reference
        .step_panel_coupled_with_transport(0.12, transport(128), || false)
        .unwrap();
    let mut errors = Vec::new();
    for substeps in [1, 2, 4, 8] {
        let mut result = start.clone();
        let report = result
            .step_panel_coupled_with_transport(0.12, transport(substeps), || false)
            .unwrap();
        assert_eq!(report.substeps, substeps);
        assert_eq!(report.coupled_solves, 2 * substeps + 1);
        assert_eq!(result.history().len(), start.history().len() + 1);
        assert_eq!(result.wake().len(), start.wake().len() + 1);
        assert_eq!(result.history().last().unwrap().t, 4.0 * dt);
        for (particle, old) in result.wake().iter().zip(&previous_strengths) {
            assert_eq!(particle.gamma.to_bits(), *old);
        }
        assert_boundary(&result, &panels, 0.12, core);
        let error = result
            .wake()
            .iter()
            .zip(reference.wake())
            .map(|(a, b)| (a.pos[0] - b.pos[0]).hypot(a.pos[1] - b.pos[1]))
            .fold(0.0_f64, f64::max);
        errors.push(error);
        if substeps == 1 {
            let mut original = start.clone();
            original.step_panel_coupled_at_angle(0.12).unwrap();
            assert_same_native_state(&original, &result);
        }
    }
    assert!(
        errors[0] > 1e-9,
        "fixture must resolve transport error: {errors:?}"
    );
    for pair in errors.windows(2) {
        assert!(
            pair[1] < 0.4 * pair[0],
            "midpoint refinement must converge: {errors:?}"
        );
    }
}

#[test]
fn core_and_newborn_body_clearance_guard_can_be_resolved_with_substeps() {
    let foil = naca4_symmetric(0.12, 24).unwrap();
    let panels = PreparedPanel2d::new(&foil).unwrap();
    let mut start = WakeSim::new_panel_coupled(&foil, 0.08, 0.02, 0.03).unwrap();
    start.step().unwrap();
    let mut guarded = start.clone();
    let mut options = transport(1);
    options.max_displacement_fraction = Some(0.5);
    assert!(matches!(
        guarded.step_panel_coupled_with_transport(0.08, options, || false),
        Err(BemError::InvalidScalar {
            name: "panel wake stage displacement fraction",
            ..
        })
    ));
    assert_same_native_state(&guarded, &start);
    options.substeps = 32;
    options.max_coupled_solves = 65;
    let report = guarded
        .step_panel_coupled_with_transport(0.08, options, || false)
        .unwrap();
    assert!(report.max_displacement_fraction > 0.0);
    assert!(report.max_displacement_fraction <= 0.5);
    assert_eq!(guarded.history().len(), start.history().len() + 1);
    assert_eq!(guarded.wake().len(), start.wake().len() + 1);
    assert_boundary(&guarded, &panels, 0.08, 0.03);
}

#[test]
fn transport_budget_and_cancellation_refuse_the_whole_interval() {
    let foil = naca4_symmetric(0.12, 24).unwrap();
    let mut start = WakeSim::new_panel_coupled(&foil, 0.08, 0.02, 0.03).unwrap();
    start.step().unwrap();
    let mut result = start.clone();
    let mut options = transport(4);
    options.max_coupled_solves -= 1;
    assert!(matches!(
        result.step_panel_coupled_with_transport(0.12, options, || false),
        Err(BemError::WorkEnvelopeExceeded {
            operation: "panel wake coupled solves",
            requested: 9,
            max: 8,
        })
    ));
    assert_same_native_state(&result, &start);
    let mut polls = 0usize;
    result
        .step_panel_coupled_with_transport(0.12, transport(4), || {
            polls += 1;
            false
        })
        .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut candidate = start.clone();
        let mut visits = 0usize;
        assert!(matches!(
            candidate.step_panel_coupled_with_transport(0.12, transport(4), || {
                visits += 1;
                visits == stop
            }),
            Err(BemError::Cancelled {
                operation: "panel wake transport"
            })
        ));
        assert_eq!(visits, stop);
        assert_same_native_state(&candidate, &start);
    }
}

#[test]
fn invalid_transport_settings_leave_the_native_state_unchanged() {
    let foil = naca4_symmetric(0.12, 24).unwrap();
    let mut start = WakeSim::new_panel_coupled(&foil, 0.08, 0.02, 0.03).unwrap();
    start.step().unwrap();
    let mut options = vec![
        PanelWakeTransportOptions {
            substeps: 0,
            ..transport(1)
        },
        PanelWakeTransportOptions {
            substeps: usize::MAX,
            max_coupled_solves: usize::MAX,
            ..transport(1)
        },
    ];
    options.extend([f64::NAN, f64::INFINITY, 0.0, -1.0, 1.01].map(|fraction| {
        PanelWakeTransportOptions {
            max_displacement_fraction: Some(fraction),
            ..transport(1)
        }
    }));
    for options in options {
        let mut result = start.clone();
        assert!(
            result
                .step_panel_coupled_with_transport(0.12, options, || false)
                .is_err()
        );
        assert_same_native_state(&result, &start);
    }
}
