//! Independent attached-ring dynamics and atomic self-induced shedding checks.
use fs_vpm::farfield3d::WakeCoreEvolutionMode;
use fs_vpm::filament3d::{FilamentWake, SelfInducedShed, WakeRateCertificate};

const CHORD: f64 = 1.0;
const SPAN: f64 = 1.4;
const CORE: f64 = 0.2;

fn ring(gamma: f64) -> FilamentWake {
    let mut wake = FilamentWake::new(
        WakeRateCertificate {
            shed_hz: 2.0,
            n_stations: 1,
            max_rows: 128,
        },
        vec![[0.0, -SPAN / 2.0, 0.0], [0.0, SPAN / 2.0, 0.0]],
    )
    .unwrap();
    wake.shed(&[gamma], [0.0; 3]).unwrap();
    wake.rows[0].nodes[0][0] = -CHORD;
    wake.rows[0].nodes[1][0] = -CHORD;
    wake.rows[0].core2_m2 = CORE * CORE;
    wake
}

fn config(substeps: usize) -> SelfInducedShed {
    SelfInducedShed {
        ambient_velocity_mps: [0.0; 3],
        emitted_core_radius_m: CORE,
        core_mode: WakeCoreEvolutionMode::Frozen,
        substeps,
        max_displacement_fraction: 1.0,
        max_work: 10_000_000,
    }
}

// At either free corner, the two incident straight segments induce zero;
// independently integrating the two opposite edges gives this angular rate.
// Rotation about the fixed span preserves chord, span, and that rate.
fn angular_rate(gamma: f64, core2: f64) -> f64 {
    gamma * SPAN / (4.0 * std::f64::consts::PI * (CHORD * CHORD + SPAN * SPAN + core2).sqrt())
        * (1.0 / (CHORD * CHORD + core2) + 1.0 / (SPAN * SPAN + core2))
}

fn expected_position(angle: f64, y: f64) -> [f64; 3] {
    [-CHORD * angle.cos(), y, CHORD * angle.sin()]
}

fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    (a[0] - b[0]).hypot(a[1] - b[1]).hypot(a[2] - b[2])
}

fn assert_close(a: [f64; 3], b: [f64; 3], tolerance: f64) {
    assert!(
        distance(a, b) < tolerance,
        "{a:?} != {b:?}; error {}",
        distance(a, b)
    );
}

#[test]
fn anchored_rectangle_matches_analytic_motion_and_support_impulse() {
    let mut wake = ring(4.0);
    let line = wake.line_nodes.clone();
    let initial_impulse = wake.hydrodynamic_impulse_per_density().unwrap();
    let angle = angular_rate(4.0, CORE * CORE) / wake.cert.shed_hz;
    wake.shed_self_induced(&[7.5], config(64), || false)
        .unwrap();
    for p in &wake.rows[0].nodes {
        assert_close(*p, expected_position(angle, p[1]), 2e-7);
    }
    assert_eq!(wake.line_nodes, line);
    assert_eq!(wake.rows[0].gamma, [4.0]);
    assert_eq!(wake.rows[1].gamma, [7.5]);
    assert_eq!(wake.rows[0].age_ticks, 1.0);
    assert_eq!(wake.rows[1].age_ticks, 0.0);
    assert_eq!(wake.rows[1].nodes, line);
    let impulse = wake.hydrodynamic_impulse_per_density().unwrap();
    assert_close(
        impulse,
        [
            4.0 * SPAN * CHORD * angle.sin(),
            0.0,
            4.0 * SPAN * CHORD * angle.cos(),
        ],
        2e-6,
    );
    assert!(
        distance(impulse, initial_impulse) > 0.5,
        "fixed support exchanges vector impulse"
    );
    for row in 0..=wake.rows.len() {
        for edge in 0..=wake.cert.n_stations {
            assert_eq!(wake.node_circulation_residual(row, edge), Some(0.0));
        }
    }
}

#[test]
fn substeps_show_second_order_convergence_to_the_independent_ring_solution() {
    let angle = angular_rate(4.0, CORE * CORE) / 2.0;
    let exact = expected_position(angle, -SPAN / 2.0);
    let errors: Vec<_> = [1, 2, 4, 8]
        .into_iter()
        .map(|substeps| {
            let mut wake = ring(4.0);
            wake.shed_self_induced(&[4.0], config(substeps), || false)
                .unwrap();
            distance(wake.rows[0].nodes[0], exact)
        })
        .collect();
    for pair in errors.windows(2) {
        let ratio = pair[0] / pair[1];
        assert!(
            (3.7..4.3).contains(&ratio),
            "RK2 refinement {errors:?}, ratio {ratio}"
        );
    }
}

#[test]
fn circulation_reversal_reverses_roll_up() {
    let mut positive = ring(4.0);
    let mut negative = ring(-4.0);
    positive
        .shed_self_induced(&[4.0], config(16), || false)
        .unwrap();
    negative
        .shed_self_induced(&[-4.0], config(16), || false)
        .unwrap();
    for (a, b) in positive.rows[0].nodes.iter().zip(&negative.rows[0].nodes) {
        assert_close(*a, [b[0], b[1], -b[2]], 2e-14);
        assert!(a[2] > 0.0 && b[2] < 0.0);
    }
}

#[test]
fn zero_circulation_has_exact_uniform_convection_and_one_endpoint_emission() {
    let mut wake = ring(0.0);
    wake.rows[0].age_ticks = 2.25;
    let original = wake.clone();
    let mut settings = config(3);
    settings.ambient_velocity_mps = [0.25, -0.5, 0.125];
    wake.shed_self_induced(&[3.0], settings, || false).unwrap();
    for (before, after) in original.rows[0].nodes.iter().zip(&wake.rows[0].nodes) {
        assert_close(
            *after,
            std::array::from_fn(|c| {
                before[c] + settings.ambient_velocity_mps[c] / original.cert.shed_hz
            }),
            1e-14,
        );
    }
    assert_eq!(wake.rows.len(), original.rows.len() + 1);
    assert_eq!(wake.ticks, original.ticks + 1);
    assert_eq!(wake.rows[0].age_ticks, 3.25);
    assert_eq!(wake.rows[1].gamma, [3.0]);
    assert_eq!(wake.rows[1].core2_m2, CORE * CORE);
    assert_eq!(wake.line_nodes, original.line_nodes);
}

#[test]
fn spreading_uses_elapsed_stage_age_without_double_aging_or_core_reset() {
    let mut wake = ring(4.0);
    wake.rows[0].age_ticks = 2.0;
    let mut settings = config(64);
    settings.core_mode = WakeCoreEvolutionMode::CoreSpreading { nu_eff_m2ps: 0.01 };
    let dt = 1.0 / wake.cert.shed_hz;
    // Independent Simpson integration of the analytic time-dependent angular
    // rate; core²(t) = base² + 4ν(age * shed_dt + t).
    let n = 2048;
    let rate = |t: f64| angular_rate(4.0, CORE * CORE + 0.04 * (2.0 * dt + t));
    let mut integral = rate(0.0) + rate(dt);
    for i in 1..n {
        integral += if i % 2 == 0 { 2.0 } else { 4.0 } * rate(dt * i as f64 / n as f64);
    }
    let angle = integral * dt / (3.0 * n as f64);
    wake.shed_self_induced(&[4.0], settings, || false).unwrap();
    assert_close(
        wake.rows[0].nodes[0],
        expected_position(angle, -SPAN / 2.0),
        2e-7,
    );
    assert_eq!(wake.rows[0].core2_m2, CORE * CORE);
    assert_eq!(wake.rows[0].age_ticks, 3.0);
    assert_eq!(wake.rows[1].core2_m2, CORE * CORE);
    assert_eq!(wake.rows[1].age_ticks, 0.0);
}

#[test]
fn newly_emitted_zero_length_edges_release_under_the_positive_core_guard() {
    let mut wake = FilamentWake::new(
        WakeRateCertificate {
            shed_hz: 100.0,
            n_stations: 1,
            max_rows: 16,
        },
        vec![[0.0, -0.5, 0.0], [0.0, 0.5, 0.0]],
    )
    .unwrap();
    let mut settings = config(4);
    settings.ambient_velocity_mps = [-1.0, 0.0, 0.0];
    settings.emitted_core_radius_m = 0.05;
    settings.max_displacement_fraction = 0.5;
    for _ in 0..3 {
        wake.shed_self_induced(&[1.0], settings, || false).unwrap();
    }
    assert_eq!(wake.rows.len(), 3);
    let old = wake.rows[0].nodes[0];
    let young = wake.rows[1].nodes[0];
    assert!(old[0] < young[0] && young[0] < 0.0);
    assert!(
        wake.rows[0].nodes[0][2] > 0.0,
        "wake actually induces motion"
    );
    assert_eq!(wake.rows[2].nodes, wake.line_nodes);
    for row in 0..=wake.rows.len() {
        for edge in 0..=1 {
            assert_eq!(wake.node_circulation_residual(row, edge), Some(0.0));
        }
    }
}

#[test]
fn invalid_core_interface_work_and_displacement_refuse_atomically() {
    let original = ring(4.0);
    let mut invalid_core = original.clone();
    invalid_core.rows[0].core2_m2 = 0.0;
    let mut attached = original.clone();
    attached.downstream_gamma = Some(vec![4.0]);
    let mut tight = config(2);
    tight.max_work = 64; // Exact boundary: (2 * 2 * (2 + 1) + 1) * (4 + 1).
    let mut too_large = config(2);
    too_large.max_displacement_fraction = 1e-8;
    let mut overflow = config(usize::MAX);
    overflow.max_work = u64::MAX;
    let cases = [
        (invalid_core, config(2), "wake-advection-core-invalid"),
        (attached, config(2), "wake-advection-interface-attached"),
        (original.clone(), tight, "wake-advection-work-exhausted"),
        (original.clone(), too_large, "wake-advection-displacement"),
        (original.clone(), overflow, "wake-advection-work-exhausted"),
    ];
    for (mut wake, settings, code) in cases {
        let before = wake.clone();
        assert_eq!(
            wake.shed_self_induced(&[4.0], settings, || false)
                .unwrap_err()
                .code,
            code
        );
        assert_eq!(wake, before);
    }
    let mut at_cap = original;
    tight.max_work = 65;
    at_cap.shed_self_induced(&[4.0], tight, || false).unwrap();
    let mut extreme = ring(f64::MAX);
    let before = extreme.clone();
    assert_eq!(
        extreme
            .shed_self_induced(&[-f64::MAX], config(1), || false)
            .unwrap_err()
            .code,
        "wake-advection-circulation-invalid"
    );
    assert_eq!(extreme, before);
}

#[test]
fn cancellation_interrupts_inside_a_source_target_sweep_without_publication() {
    let mut wake = FilamentWake::new(
        WakeRateCertificate {
            shed_hz: 100.0,
            n_stations: 3,
            max_rows: 128,
        },
        (0..4).map(|i| [0.0, i as f64, 0.0]).collect(),
    )
    .unwrap();
    for i in 0..100 {
        let g = 0.25 + i as f64 / 128.0;
        wake.shed(&[g, 2.0 * g, 0.5 * g], [-0.01, 0.0, 0.0])
            .unwrap();
    }
    for row in &mut wake.rows {
        row.core2_m2 = 0.01;
    }
    let before = wake.clone();
    let mut polls = 0;
    let error = wake
        .shed_self_induced(&[1.0, 2.0, 0.5], config(1), || {
            polls += 1;
            polls == 4 // start, substep, first target, then its 256-source checkpoint
        })
        .unwrap_err();
    assert_eq!(error.code, "wake-advection-cancelled");
    assert_eq!(polls, 4);
    assert_eq!(wake, before);
}

#[test]
fn dynamics_are_repeatable_and_covariant_under_a_proper_rotation_and_translation() {
    let original = ring(4.0);
    let transform = |p: [f64; 3]| [p[2] + 2.0, p[0] - 3.0, p[1] + 0.5];
    let mut rotated = original.clone();
    for p in &mut rotated.line_nodes {
        *p = transform(*p);
    }
    for row in &mut rotated.rows {
        for p in &mut row.nodes {
            *p = transform(*p);
        }
    }
    let mut reference = original.clone();
    let mut repeated = original;
    for wake in [&mut reference, &mut repeated, &mut rotated] {
        wake.shed_self_induced(&[4.0], config(16), || false)
            .unwrap();
    }
    assert_eq!(reference, repeated);
    for (row, turned) in reference.rows.iter().zip(&rotated.rows) {
        for (a, b) in row.nodes.iter().zip(&turned.nodes) {
            assert_close(transform(*a), *b, 2e-14);
        }
    }
}
