use super::*;
use fs_contact::capillary_bridge::load_line::{
    BridgeLoadLineConfig, BridgeLoadLineError, PlanarBridgeLoadLine,
};

fn solve_config() -> BridgeLoadLineConfig {
    BridgeLoadLineConfig {
        max_iterations: 32,
        max_evaluations: 160,
        max_backtracks: 16,
        gap_absolute_tolerance_m: 2e-16,
        gap_relative_tolerance: 1e-11,
        maximum_gap_increment_m: 0.01,
        minimum_compatibility_slope: 1e-4,
    }
}

#[test]
fn exact_mean_force_derivative_tracks_attached_ruptured_and_reformed_paths() {
    let m = model(PI / 3.0, 0.05, true);
    let gate = CancelGate::new();
    for (before, gap) in [
        (
            m.initialize(0.6 * m.rupture_gap_m(), true).unwrap(),
            0.8 * m.rupture_gap_m(),
        ),
        (
            m.initialize(0.9 * m.rupture_gap_m(), true).unwrap(),
            1.2 * m.rupture_gap_m(),
        ),
        (
            m.initialize(1.5 * m.formation_gap_m(), false).unwrap(),
            0.8 * m.formation_gap_m(),
        ),
    ] {
        let step = m.advance(before, gap, 10.0, &gate).unwrap();
        let dh = 1e-5 * gap;
        let lo = m.advance(before, gap - dh, 10.0, &gate).unwrap();
        let hi = m.advance(before, gap + dh, 10.0, &gate).unwrap();
        near(
            (hi.mean_normal_force_n - lo.mean_normal_force_n) / (2.0 * dh),
            step.mean_force_gap_derivative_n_per_m,
            3e-8,
            1e-4,
        );
        assert!(step.mean_force_gap_derivative_n_per_m > 0.0);
    }
    let old = m.initialize(0.7 * m.rupture_gap_m(), true).unwrap();
    let rest = m.advance(old, old.gap_m(), 1.0, &gate).unwrap();
    near(
        rest.mean_force_gap_derivative_n_per_m,
        0.5 * rest.profile.force_gap_derivative_n_per_m,
        2e-14,
        0.0,
    );
    let touching = m.initialize(m.formation_gap_m(), false).unwrap();
    let coalescence = m.advance(touching, touching.gap_m(), 1.0, &gate).unwrap();
    assert_eq!(coalescence.event, BridgeEvent::Reformed);
    assert!(coalescence.mean_force_gap_derivative_n_per_m.is_finite());
}

#[test]
fn original_two_body_motion_preserves_momentum_and_surface_work_through_rupture() {
    let m = model(PI / 3.0, 0.05, true);
    let gate = CancelGate::new();
    let (ma, mb, dt) = (2.0, 3.0, 0.001);
    let initial_gap = 0.5 * (m.formation_gap_m() + m.rupture_gap_m());
    let mut state = m.initialize(initial_gap, true).unwrap();
    let mut xa = 0.0;
    let mut xb = initial_gap;
    let mut va = -0.0006;
    let mut vb = 0.0004;
    let momentum = ma * va + mb * vb;
    let kinetic = |va: f64, vb: f64| 0.5 * ma * va * va + 0.5 * mb * vb * vb;
    let initial = kinetic(va, vb) + m.evaluate(state).unwrap().surface_free_energy_j;
    let solver =
        PlanarBridgeLoadLine::new(&m, 0.5 * dt * dt * (1.0 / ma + 1.0 / mb), solve_config())
            .unwrap();
    let mut released = 0.0;
    let mut ruptures = 0;
    for _ in 0..700 {
        let old_kinetic = kinetic(va, vb);
        let free = xb - xa + dt * (vb - va);
        let trial = solver.trial(state, free, dt, &gate).unwrap();
        let force = trial.bridge.mean_normal_force_n;
        xa += dt * va - 0.5 * dt * dt * force / ma;
        xb += dt * vb + 0.5 * dt * dt * force / mb;
        va -= dt * force / ma;
        vb += dt * force / mb;
        let actual = m.advance(state, xb - xa, dt, &gate).unwrap();
        near(actual.mean_normal_force_n, force, 1e-9, 1e-7);
        near(xb - xa, trial.bridge.state().gap_m(), 2e-12, 1e-6);
        near(ma * va + mb * vb, momentum, 2e-13, 1e-3);
        near(
            kinetic(va, vb) - old_kinetic + actual.work_on_liquid_j,
            0.0,
            1e-10,
            1e-8,
        );
        released += actual.irreversible_release_j;
        ruptures += usize::from(actual.event == BridgeEvent::Ruptured);
        near(
            kinetic(va, vb) + actual.profile.surface_free_energy_j + released,
            initial,
            1e-10,
            1e-6,
        );
        assert!(trial.evaluations <= solve_config().max_evaluations);
        state = actual.state();
    }
    assert_eq!(ruptures, 1);
    assert!(!state.attached());
    assert!(released > 0.0);
    assert_eq!(
        m.evaluate(state).unwrap().detached_cap_volumes_m3,
        Some([0.5 * VOLUME; 2])
    );
}

#[test]
fn prescribed_and_compliant_formation_use_actual_path_force_and_common_source() {
    let m = model(PI / 3.0, 0.05, true);
    let gate = CancelGate::new();
    let old = m.initialize(1.1 * m.formation_gap_m(), false).unwrap();
    let endpoint = 0.9 * m.formation_gap_m();
    let expected = m.advance(old, endpoint, 10.0, &gate).unwrap();
    for compliance in [0.0, 0.02] {
        let free = endpoint - compliance * expected.mean_normal_force_n;
        let owner = PlanarBridgeLoadLine::new(&m, compliance, solve_config()).unwrap();
        assert!(std::ptr::eq(owner.source(), &m));
        let result = owner.trial(old, free, 10.0, &gate).unwrap();
        assert_eq!(result.bridge.event, BridgeEvent::Reformed);
        near(result.bridge.state().gap_m(), endpoint, 2e-11, 1e-6);
        near(
            result.bridge.mean_normal_force_n,
            expected.mean_normal_force_n,
            2e-10,
            1e-8,
        );
        assert!(result.bridge.irreversible_release_j > 0.0);
        near(
            result.bridge.surface_energy_change_j + result.bridge.irreversible_release_j,
            result.bridge.work_on_liquid_j,
            3e-12,
            1e-8,
        );
    }
}

#[test]
fn unresolved_attraction_source_limits_and_cancellation_never_publish_a_successor() {
    let m = model(PI / 3.0, 0.05, true);
    let gate = CancelGate::new();
    let old = m.initialize(0.5 * m.rupture_gap_m(), true).unwrap();
    let response = m.advance(old, old.gap_m(), 1.0, &gate).unwrap();
    let unstable = PlanarBridgeLoadLine::new(
        &m,
        2.0 / response.mean_force_gap_derivative_n_per_m,
        solve_config(),
    )
    .unwrap();
    assert!(matches!(
        unstable.trial(old, old.gap_m(), 1.0, &gate),
        Err(BridgeLoadLineError::CompatibilitySlope { .. })
    ));
    let mut c = solve_config();
    c.max_evaluations = 1;
    let limited = PlanarBridgeLoadLine::new(&m, 1e-5, c).unwrap();
    assert_eq!(
        limited.trial(old, old.gap_m() + 1e-5, 1.0, &gate),
        Err(BridgeLoadLineError::Budget)
    );
    let owner = PlanarBridgeLoadLine::new(&m, 1e-5, solve_config()).unwrap();
    let first = owner.trial(old, old.gap_m() + 1e-5, 1.0, &gate).unwrap();
    let cancelled = CancelGate::new();
    cancelled.request();
    assert_eq!(
        owner.trial(old, old.gap_m(), 1.0, &cancelled),
        Err(BridgeLoadLineError::Cancelled)
    );
    assert!(owner.trial(old, -1.0, 1.0, &gate).is_err());
    let foreign = model(PI / 4.0, 0.05, true)
        .initialize(old.gap_m(), true)
        .unwrap();
    assert!(matches!(
        owner.trial(foreign, old.gap_m(), 1.0, &gate),
        Err(BridgeLoadLineError::Bridge(BridgeError::StateMismatch))
    ));
    assert_eq!(
        owner.trial(old, old.gap_m() + 1e-5, 1.0, &gate).unwrap(),
        first
    );
    assert_eq!(old.gap_m(), 0.5 * m.rupture_gap_m());
}
