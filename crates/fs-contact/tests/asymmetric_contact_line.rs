//! G1/G3/G4: unlike-surface geometry, work, source and dynamic wetting.
#[path = "asymmetric_contact_line/common.rs"]
mod common;
#[path = "contact_line_bridge/common.rs"]
mod symmetric;
use common::*;
use fs_contact::capillary_bridge::asymmetric_contact_line::*;
use fs_exec::CancelGate;
use fs_material::liquid::surface::adhesion::EquilibriumAngleOrigin;

#[test]
fn g1_independent_circular_volume_plate_resultants_and_energy_gradients() {
    let m = source([0.8, 1.1], [4.0, 13.0], [0.7, 1.0], [0.9, 1.2]);
    let initial = m.initialize(GAP, [0.85, 1.2]).unwrap();
    let p = m.evaluate(initial).unwrap();
    let radius = p.meniscus_radius_m;
    let xc = p.wetted_spans_m[0] / 2.0 + radius * p.contact_angles_radians[0].sin();
    let yc = radius * p.contact_angles_radians[0].cos();
    let n = 8192;
    let dy = GAP / f64::from(n);
    let mut integral = 0.0;
    for i in 0..=n {
        let y = f64::from(i) * dy;
        let width = 2.0 * (xc - (radius * radius - (y - yc) * (y - yc)).sqrt());
        let weight = if i == 0 || i == n {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        integral += weight * width;
    }
    close(WIDTH * integral * dy / 3.0, VOLUME, 2e-23);
    close(p.neck_height_m, yc, 2e-18);
    close(p.neck_width_m, 2.0 * (xc - radius), 2e-18);
    let forces = p
        .wetted_spans_m
        .map(|b| p.liquid_minus_gas_pressure_pa * b * WIDTH);
    let lower = -(forces[0] - 2.0 * SIGMA * WIDTH * p.contact_angles_radians[0].sin());
    let upper = forces[1] - 2.0 * SIGMA * WIDTH * p.contact_angles_radians[1].sin();
    close(lower + upper, 0.0, 2e-18);
    close(upper, p.normal_force_n, 2e-18);
    let epsilon = 1e-8;
    for axis in 0..3 {
        let (mut hp, mut hm) = (GAP, GAP);
        let (mut bp, mut bm) = (p.wetted_spans_m, p.wetted_spans_m);
        if axis == 0 {
            hp += epsilon;
            hm -= epsilon;
        } else {
            bp[axis - 1] += epsilon;
            bm[axis - 1] -= epsilon;
        }
        let ep = m
            .evaluate(
                m.initialize_spans(hp, bp, &CancelGate::new_clock_free())
                    .unwrap(),
            )
            .unwrap()
            .surface_free_energy_j;
        let em = m
            .evaluate(
                m.initialize_spans(hm, bm, &CancelGate::new_clock_free())
                    .unwrap(),
            )
            .unwrap()
            .surface_free_energy_j;
        let reference = if axis == 0 {
            -p.normal_force_n
        } else {
            p.footprint_energy_gradients_n[axis - 1]
        };
        close((ep - em) / (2.0 * epsilon), reference, 2e-11);
    }
}

#[test]
fn g3_symmetric_sources_reduce_to_the_original_complete_discrete_wetting_law() {
    let m = source([1.0; 2], [10.0; 2], [0.95; 2], [1.05; 2]);
    let original = symmetric::source(1.0, 10.0, 0.95, 1.05);
    let mut a = m.initialize(GAP, [1.12; 2]).unwrap();
    let mut b = original.initialize(GAP, 1.12).unwrap();
    for i in 0..30 {
        let h = GAP + 1e-6 * (f64::from(i) / 10.0).sin();
        let sa = m
            .advance(a, h, 0.001, &CancelGate::new_clock_free())
            .unwrap();
        let sb = original
            .advance(b, h, 0.001, &CancelGate::new_clock_free())
            .unwrap();
        for span in sa.state().wetted_spans_m() {
            close(span, sb.state().wetted_span_m(), 2e-14);
        }
        close(sa.mean_normal_force_n, sb.mean_normal_force_n, 3e-13);
        close(
            sa.surface_energy_change_j,
            sb.surface_energy_change_j,
            2e-18,
        );
        close(
            sa.contact_line_dissipation_j,
            sb.contact_line_dissipation_j,
            2e-18,
        );
        close(
            sa.mean_force_gap_derivative_n_per_m,
            sb.mean_force_gap_derivative_n_per_m,
            1e-7,
        );
        a = sa.state();
        b = sb.state();
    }
}

#[test]
fn g1_independent_sources_drive_one_surface_while_the_other_stays_pinned() {
    let m = source([0.9, 1.1], [4.0, 15.0], [0.8, 1.09], [1.0, 1.11]);
    let initial = m.initialize(GAP, [0.94, 1.2]).unwrap();
    let s = m
        .advance(initial, GAP + 1e-6, 0.001, &CancelGate::new_clock_free())
        .unwrap();
    assert_eq!(
        s.regimes,
        [ContactLineRegime::Pinned, ContactLineRegime::Advancing]
    );
    assert_eq!(s.state().wetted_spans_m()[0], initial.wetted_spans_m()[0]);
    assert!(s.state().wetted_spans_m()[1] > initial.wetted_spans_m()[1]);
    assert_eq!(s.contact_line_dissipation_by_surface_j[0], 0.0);
    assert!(
        s.viscous_line_dissipation_by_surface_j[1] > 0.0
            && s.pinning_dissipation_by_surface_j[1] > 0.0
    );
    close(
        s.surface_energy_change_j + s.contact_line_dissipation_j,
        s.work_on_liquid_j,
        5e-19,
    );
    let fast = source([0.9, 1.1], [4.0, 3.0], [0.8, 1.09], [1.0, 1.11]);
    assert_eq!(
        m.wetting()[0].resolved().identity(),
        fast.wetting()[0].resolved().identity()
    );
    assert_ne!(
        m.wetting()[1].resolved().identity(),
        fast.wetting()[1].resolved().identity()
    );
    let f = fast
        .advance(
            fast.initialize(GAP, [0.94, 1.2]).unwrap(),
            GAP + 1e-6,
            0.001,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert!(f.state().wetted_spans_m()[1] > s.state().wetted_spans_m()[1]);
    let recede = m
        .advance(
            m.initialize(GAP, [0.78, 1.1]).unwrap(),
            GAP,
            0.001,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(
        recede.regimes,
        [ContactLineRegime::Receding, ContactLineRegime::Pinned]
    );
    assert!(recede.contact_line_dissipation_by_surface_j[0] > 0.0);
    assert_eq!(recede.contact_line_dissipation_by_surface_j[1], 0.0);
}

#[test]
fn g1_implicit_normal_tangent_matches_actual_coupled_motion_on_each_active_set() {
    let m = source([0.9, 1.1], [4.0, 13.0], [0.88, 1.08], [0.92, 1.12]);
    for angles in [[0.9, 1.1], [1.05, 1.24], [0.75, 0.98], [0.9, 1.24]] {
        let state = m.initialize(GAP, angles).unwrap();
        let h = GAP + 1e-6;
        let eps = 1e-9;
        let s = m
            .advance(state, h, 0.002, &CancelGate::new_clock_free())
            .unwrap();
        let plus = m
            .advance(state, h + eps, 0.002, &CancelGate::new_clock_free())
            .unwrap();
        let minus = m
            .advance(state, h - eps, 0.002, &CancelGate::new_clock_free())
            .unwrap();
        let derivative = (plus.mean_normal_force_n - minus.mean_normal_force_n) / (2.0 * eps);
        close(
            s.mean_force_gap_derivative_n_per_m,
            derivative,
            4e-6 * derivative.abs().max(1.0),
        );
    }
}

#[test]
fn g3_unlike_surface_relaxation_converges_to_independent_two_coordinate_ode() {
    let young = [0.9, 1.1];
    let zeta = [4.0, 13.0];
    let m = source(young, zeta, young, young);
    let initial = m.initialize(GAP, [1.02, 1.2]).unwrap();
    let rhs = |b: [f64; 2]| {
        let angles = angles_at(GAP, b);
        core::array::from_fn::<_, 2, _>(|i| {
            2.0 * SIGMA / zeta[i] * (young[i].cos() - angles[i].cos())
        })
    };
    let end = 0.02;
    let dt = end / 4096.0;
    let mut reference = initial.wetted_spans_m();
    for _ in 0..4096 {
        let k1 = rhs(reference);
        let k2 = rhs(core::array::from_fn(|i| reference[i] + 0.5 * dt * k1[i]));
        let k3 = rhs(core::array::from_fn(|i| reference[i] + 0.5 * dt * k2[i]));
        let k4 = rhs(core::array::from_fn(|i| reference[i] + dt * k3[i]));
        reference = core::array::from_fn(|i| {
            reference[i] + dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0
        });
    }
    let errors = [10, 20, 40].map(|n| {
        let mut state = initial;
        let mut loss = 0.0;
        for _ in 0..n {
            let s = m
                .advance(
                    state,
                    GAP,
                    end / f64::from(n),
                    &CancelGate::new_clock_free(),
                )
                .unwrap();
            loss += s.contact_line_dissipation_j;
            state = s.state();
        }
        close(
            m.evaluate(state).unwrap().surface_free_energy_j + loss,
            m.evaluate(initial).unwrap().surface_free_energy_j,
            5e-18,
        );
        let actual = state.wetted_spans_m();
        (actual[0] - reference[0]).hypot(actual[1] - reference[1])
    });
    assert!(
        errors[0] / errors[1] > 3.5 && errors[1] / errors[2] > 3.5,
        "{errors:?}"
    );
}

#[test]
fn g3_plate_exchange_preserves_force_energy_and_each_original_loss_owner() {
    let (liquid, wetting) = sources([0.8, 1.1], [4.0, 13.0], [0.75, 1.05], [0.85, 1.15]);
    let m = SourceAsymmetricContactLineBridge::new(
        liquid.clone(),
        wetting.clone(),
        VOLUME,
        WIDTH,
        limits(),
        config(),
    )
    .unwrap();
    let other = SourceAsymmetricContactLineBridge::new(
        liquid,
        [wetting[1].clone(), wetting[0].clone()],
        VOLUME,
        WIDTH,
        limits(),
        config(),
    )
    .unwrap();
    assert_eq!(m.wetting()[0], other.wetting()[1]);
    assert_eq!(m.wetting()[1], other.wetting()[0]);
    let a = m.initialize(GAP, [0.9, 1.2]).unwrap();
    let b = other.initialize(GAP, [1.2, 0.9]).unwrap();
    let sa = m
        .advance(a, GAP + 1e-6, 0.003, &CancelGate::new_clock_free())
        .unwrap();
    let sb = other
        .advance(b, GAP + 1e-6, 0.003, &CancelGate::new_clock_free())
        .unwrap();
    close(sa.mean_normal_force_n, sb.mean_normal_force_n, 3e-17);
    close(
        sa.profile.surface_free_energy_j,
        sb.profile.surface_free_energy_j,
        2e-21,
    );
    for i in 0..2 {
        close(
            sa.state().wetted_spans_m()[i],
            sb.state().wetted_spans_m()[1 - i],
            2e-17,
        );
        close(
            sa.contact_line_dissipation_by_surface_j[i],
            sb.contact_line_dissipation_by_surface_j[1 - i],
            2e-21,
        );
    }
    assert_eq!(other.evaluate(a), Err(ContactLineError::StateMismatch));
}

#[test]
fn g0_g4_original_source_identity_domains_accuracy_and_budget_refuse_atomically() {
    let (liquid, wetting) = sources_with_origins(
        [0.8, 1.1],
        [4.0, 13.0],
        [0.7, 1.0],
        [0.9, 1.2],
        [
            EquilibriumAngleOrigin::SourcedAngle,
            EquilibriumAngleOrigin::YoungDupre,
        ],
    );
    let m =
        SourceAsymmetricContactLineBridge::new(liquid, wetting, VOLUME, WIDTH, limits(), config())
            .unwrap();
    assert_eq!(
        m.wetting()[1].angle_origin(),
        EquilibriumAngleOrigin::YoungDupre
    );
    let state = m.initialize(GAP, [0.85, 1.2]).unwrap();
    let saved = m.evaluate(state).unwrap();
    let gate = CancelGate::new_clock_free();
    gate.request();
    assert_eq!(
        m.advance(state, GAP, 0.001, &gate),
        Err(ContactLineError::Cancelled)
    );
    assert!(m.initialize(GAP, [0.1, 1.2]).is_err());
    assert!(
        m.initialize_spans(GAP, [1e-7, 1e-3], &CancelGate::new_clock_free())
            .is_err()
    );
    assert!(
        m.advance(state, GAP + 1e-3, 0.001, &CancelGate::new_clock_free())
            .is_err()
    );
    let (l, w) = sources([0.8, 1.1], [4.0, 13.0], [0.7, 1.0], [0.9, 1.2]);
    let mut cfg = config();
    cfg.contact_line.maximum_evaluations = 1;
    let limited =
        SourceAsymmetricContactLineBridge::new(l, w, VOLUME, WIDTH, limits(), cfg).unwrap();
    let limited_state = limited.initialize(GAP, [0.85, 1.2]).unwrap();
    assert_eq!(
        limited.advance(limited_state, GAP, 0.001, &CancelGate::new_clock_free()),
        Err(ContactLineError::Budget)
    );
    let (l, w) = sources([0.8, 1.1], [4.0, 13.0], [0.7, 1.0], [0.9, 1.2]);
    let mut cfg = config();
    cfg.contact_line.force_absolute_tolerance_n = 1e-30;
    cfg.contact_line.force_relative_tolerance = 1e-16;
    let strict =
        SourceAsymmetricContactLineBridge::new(l, w, VOLUME, WIDTH, limits(), cfg).unwrap();
    assert!(
        strict
            .advance(
                strict.initialize(GAP, [0.85, 1.2]).unwrap(),
                GAP,
                0.001,
                &CancelGate::new_clock_free()
            )
            .is_err()
    );
    assert_eq!(m.evaluate(state).unwrap(), saved);
    assert_eq!(
        m.advance(state, GAP, 0.001, &CancelGate::new_clock_free()),
        m.advance(state, GAP, 0.001, &CancelGate::new_clock_free())
    );
}
