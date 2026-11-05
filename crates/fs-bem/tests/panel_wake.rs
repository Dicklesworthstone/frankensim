//! Native boundary/Kelvin closure, accepted-geometry coupling and frame tests.

use fs_bem::panel2d::{Airfoil2d, PreparedPanel2d, naca4_symmetric, solve};
use fs_bem::wake2d::WakeSim;

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
