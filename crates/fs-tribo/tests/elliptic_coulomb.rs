//! G1/G3: directional static capacity, associated slip, body recoil and transport.
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::tangential::elliptic::EllipticCoulomb;
use fs_tribo::tangential::matrix::{MatrixTangentialConfig, advance_with_matrix_compliance};
use fs_tribo::tangential::multiple::{
    MultipleTangentialConfig, PreparedTangentialContacts, TangentialContactInput,
};
use fs_tribo::tangential::transport::{FramedTangentialState, TangentialTransportConfig};
use fs_tribo::{
    FrictionLaw, FrictionRegime, TangentialLaw, TangentialSolveConfig, TangentialStepInput,
};

fn input(free: [f64; 2], compliance: f64) -> TangentialStepInput {
    TangentialStepInput {
        normal_n: 10.0,
        slip_increment_m: free,
        dt_s: 0.01,
        contact_temperature_k: Some(300.0),
        nominal_pressure_pa: Some(1e5),
        compliance_m_per_n: compliance,
        solve: TangentialSolveConfig {
            max_iterations: 128,
            maximum_force_n: 100.0,
            force_absolute_tolerance_n: 1e-12,
            force_relative_tolerance: 1e-12,
        },
    }
}

fn law() -> EllipticCoulomb {
    EllipticCoulomb::new([0.8, 0.4], [0.6, 0.2], [1.0, 0.0]).unwrap()
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}; tolerance {tolerance:e}"
    );
}

fn close2(a: [f64; 2], b: [f64; 2], tolerance: f64) {
    close(a[0], b[0], tolerance);
    close(a[1], b[1], tolerance);
}

fn rotate(v: [f64; 2], angle: f64) -> [f64; 2] {
    let (s, c) = angle.sin_cos();
    [c * v[0] - s * v[1], s * v[0] + c * v[1]]
}

fn matrix_controls() -> MatrixTangentialConfig {
    MatrixTangentialConfig {
        maximum_evaluations: 512,
        displacement_tolerance_m: 1e-11,
        energy_tolerance_j: 1e-10,
        matrix_relative_tolerance: 1e-12,
    }
}

#[test]
fn static_ellipse_and_weaker_sliding_have_distinct_directional_thresholds() {
    let law = law();
    let state = law.initial_state();
    // This reaction exceeds the kinetic major-axis limit but lies inside static.
    let sticking = law.advance(&state, &input([0.07, 0.0], 0.01)).unwrap();
    assert_eq!(sticking.regime, FrictionRegime::Sticking);
    close2(sticking.traction_n, [7.0, 0.0], 2e-15);
    assert_eq!(sticking.slip_increment_m, [0.0; 2]);
    assert_eq!(sticking.dissipated_j, 0.0);
    let sliding = law.advance(&state, &input([0.081, 0.0], 0.01)).unwrap();
    assert_eq!(sliding.regime, FrictionRegime::Sliding);
    assert_eq!(sliding.iterations, 0);
    close2(sliding.traction_n, [6.0, 0.0], 1e-11);
    close2(sliding.slip_increment_m, [0.021, 0.0], 1e-13);
    close(sliding.dissipated_j, 0.126, 3e-13);
    // Both components individually fit their static bounds, but not the ellipse.
    let corner = law.advance(&state, &input([0.07, 0.03], 0.01)).unwrap();
    assert_eq!(corner.regime, FrictionRegime::Sliding);
    // Rotating the surface by 90 degrees changes the physical breakaway load.
    let rotated = EllipticCoulomb::new([0.8, 0.4], [0.6, 0.2], [0.0, 1.0]).unwrap();
    let result = rotated
        .advance(&rotated.initial_state(), &input([0.07, 0.0], 0.01))
        .unwrap();
    assert_eq!(result.regime, FrictionRegime::Sliding);
    close2(result.traction_n, [2.0, 0.0], 1e-11);
}

#[test]
fn prescribed_sliding_maximizes_dissipation_on_the_entire_ellipse() {
    let law = law();
    for free in [[0.02, 0.03], [-0.015, 0.008], [0.03, -0.06]] {
        let step = law
            .advance(&law.initial_state(), &input(free, 0.0))
            .unwrap();
        let support = (36.0 * free[0] * free[0] + 4.0 * free[1] * free[1]).sqrt();
        close2(
            step.traction_n,
            [36.0 * free[0] / support, 4.0 * free[1] / support],
            1e-14,
        );
        close(step.dissipated_j, support, 1e-14);
        close(
            (step.traction_n[0] / 6.0).powi(2) + (step.traction_n[1] / 2.0).powi(2),
            1.0,
            1e-14,
        );
        // Independent angular sweep: no other feasible traction dissipates more.
        for j in 0..4096 {
            let theta = std::f64::consts::TAU * j as f64 / 4096.0;
            let other_work = 6.0 * theta.cos() * free[0] + 2.0 * theta.sin() * free[1];
            assert!(other_work <= step.dissipated_j + 1e-14);
        }
        assert!((step.traction_n[0] * free[1] - step.traction_n[1] * free[0]).abs() > 0.01);
        let reversed = law
            .advance(&law.initial_state(), &input(free.map(|x| -x), 0.0))
            .unwrap();
        close2(reversed.traction_n, step.traction_n.map(|x| -x), 1e-14);
        close(reversed.dissipated_j, step.dissipated_j, 1e-14);
    }
}

// Independent minimization over ellipse angle. Production instead solves its
// monotone multiplier; this oracle never uses its projection equation.
fn angular_projection(free: [f64; 2], compliance: [[f64; 2]; 2]) -> [f64; 2] {
    let energy = |theta: f64| {
        let f = [6.0 * theta.cos(), 2.0 * theta.sin()];
        0.5 * (f[0] * (compliance[0][0] * f[0] + compliance[0][1] * f[1])
            + f[1] * (compliance[1][0] * f[0] + compliance[1][1] * f[1]))
            - f[0] * free[0]
            - f[1] * free[1]
    };
    let spacing = std::f64::consts::TAU / 8192.0;
    let best = (0..8192)
        .min_by(|&a, &b| energy(a as f64 * spacing).total_cmp(&energy(b as f64 * spacing)))
        .unwrap();
    let mut lo = (best as f64 - 1.0) * spacing;
    let mut hi = (best as f64 + 1.0) * spacing;
    for _ in 0..80 {
        let a = lo + (hi - lo) / 3.0;
        let b = hi - (hi - lo) / 3.0;
        if energy(a) < energy(b) {
            hi = b;
        } else {
            lo = a;
        }
    }
    let theta = 0.5 * (lo + hi);
    [6.0 * theta.cos(), 2.0 * theta.sin()]
}

#[test]
fn finite_body_recoil_matches_independent_constrained_energy_minimization() {
    let law = law();
    for free in [[0.15, 0.04], [-0.13, 0.07], [0.09, -0.12]] {
        for compliance in [0.001, 0.01, 0.02] {
            let step = law
                .advance(&law.initial_state(), &input(free, compliance))
                .unwrap();
            assert_eq!(step.regime, FrictionRegime::Sliding);
            let expected = angular_projection(free, [[compliance, 0.0], [0.0, compliance]]);
            close2(step.traction_n, expected, 3e-7);
            for i in 0..2 {
                close(
                    step.slip_increment_m[i] + compliance * step.traction_n[i],
                    free[i],
                    3e-16,
                );
            }
            let f = step.traction_n;
            let du = step.slip_increment_m;
            let work = (36.0 * du[0] * du[0] + 4.0 * du[1] * du[1]).sqrt();
            close(step.dissipated_j, work, 1e-11);
            close(du[0] * f[1] / 4.0 - du[1] * f[0] / 36.0, 0.0, 1e-12);
        }
    }
}

#[test]
fn isotropic_and_frictionless_limits_match_the_original_coulomb_owner() {
    for (s, k) in [(0.6, 0.3), (0.3, 0.3), (0.0, 0.0)] {
        let elliptic = EllipticCoulomb::new([s; 2], [k; 2], [3.0, 4.0]).unwrap();
        let original = FrictionLaw::Coulomb {
            static_mu: s,
            kinetic_mu: k,
        };
        for free in [[0.0, 0.0], [0.02, -0.03], [-0.17, 0.03]] {
            for compliance in [0.0, 0.005, 0.02] {
                let x = input(free, compliance);
                let a = elliptic.advance(&elliptic.initial_state(), &x).unwrap();
                let b = original.advance(&original.initial_state(), &x).unwrap();
                assert_eq!(a.regime, b.regime);
                assert_eq!(a.iterations, 0);
                close2(a.traction_n, b.traction_n, 1e-11);
                close2(a.slip_increment_m, b.slip_increment_m, 3e-13);
                close(a.dissipated_j, b.dissipated_j, 2e-12);
            }
        }
    }
    let one_axis = EllipticCoulomb::new([0.0, 0.4], [0.0, 0.2], [1.0, 0.0]).unwrap();
    let step = one_axis
        .advance(&one_axis.initial_state(), &input([0.1, 0.01], 0.01))
        .unwrap();
    assert_eq!(step.regime, FrictionRegime::Sliding);
    close2(step.traction_n, [0.0, 1.0], 1e-15);
    close2(step.slip_increment_m, [0.1, 0.0], 1e-15);
    assert_eq!(step.dissipated_j, 0.0);
    let step = one_axis
        .advance(&one_axis.initial_state(), &input([0.1, 0.05], 0.01))
        .unwrap();
    close2(step.traction_n, [0.0, 2.0], 1e-11);
    close(step.dissipated_j, 0.06, 1e-12);

    // Free motion on the frictionless axis and a static transverse reaction
    // have exactly zero constitutive work in every texture orientation. World
    // dot-product cancellation may have either sign; it is numerical closure
    // error, not a thermal source or a reason to reject the same physics.
    for angle in [0.073, 0.292, 0.365, 0.584, 0.803, 0.949] {
        let texture =
            EllipticCoulomb::new([0.0, 0.4], [0.0, 0.2], rotate([1.0, 0.0], angle)).unwrap();
        let step = texture
            .advance(
                &texture.initial_state(),
                &input(rotate([0.123, 0.007], angle), 0.01),
            )
            .unwrap();
        close2(step.traction_n, rotate([0.0, 0.7], angle), 3e-15);
        close2(step.slip_increment_m, rotate([0.123, 0.0], angle), 1e-16);
        assert_eq!(step.dissipated_j, 0.0);
        let actual_work = step.traction_n[0] * step.slip_increment_m[0]
            + step.traction_n[1] * step.slip_increment_m[1];
        assert_eq!(step.closure_residual_j, -actual_work);
        assert!(step.closure_residual_j.abs() < 2e-17);
    }
}

#[test]
fn original_matrix_and_multiple_contact_consumers_accept_directional_laws() {
    let law = law();
    let state = law.initial_state();
    let c = [[0.008, 0.003], [0.003, 0.016]];
    let x = input([0.16, 0.09], 0.0);
    let step =
        advance_with_matrix_compliance(&law, &state, &x, c, matrix_controls(), || false).unwrap();
    close2(
        step.contact.traction_n,
        angular_projection(x.slip_increment_m, c),
        3e-7,
    );
    assert!(step.body_work_residual_j.abs() < 1e-10);
    let config = MultipleTangentialConfig {
        maximum_contacts: 2,
        maximum_matrix_terms: 16,
        maximum_sweeps: 256,
        maximum_law_evaluations: 65536,
        local: matrix_controls(),
        total_energy_tolerance_j: 1e-10,
    };
    let response = vec![
        0.008, 0.003, 0.001, 0.0, 0.003, 0.016, 0.0, 0.002, 0.001, 0.0, 0.008, 0.003, 0.0, 0.002,
        0.003, 0.016,
    ];
    let prepared = PreparedTangentialContacts::new(2, response.clone(), config, || false).unwrap();
    let sites = [
        TangentialContactInput {
            law: &law,
            state: &state,
            input: x,
        },
        TangentialContactInput {
            law: &law,
            state: &state,
            input: input([-0.15, 0.06], 0.0),
        },
    ];
    let all = prepared.advance(&sites, || false).unwrap();
    let forces: Vec<_> = all.contacts.iter().flat_map(|s| s.traction_n).collect();
    for (i, contact) in all.contacts.iter().enumerate() {
        let mut actual = [0.0; 2];
        for a in 0..2 {
            let row = 2 * i + a;
            actual[a] = sites[i].input.slip_increment_m[a]
                - (0..4)
                    .map(|j| response[4 * row + j] * forces[j])
                    .sum::<f64>();
        }
        // Direct independent support formula at the full coupled motion.
        let work = (36.0 * actual[0] * actual[0] + 4.0 * actual[1] * actual[1]).sqrt();
        close2(
            contact.traction_n,
            [36.0 * actual[0] / work, 4.0 * actual[1] / work],
            5e-9,
        );
        close(contact.dissipated_j, work, 1e-10);
    }
}

#[test]
fn texture_and_slip_transport_preserve_world_response_without_changing_history() {
    let law = law();
    let loaded = law
        .advance(&law.initial_state(), &input([0.02, 0.03], 0.0))
        .unwrap()
        .state;
    let frame = TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap();
    let before = FramedTangentialState::new(&law, loaded.clone(), frame).unwrap();
    let angle = 0.73;
    let changed = before
        .transport(
            &law,
            frame.rotated(angle).unwrap(),
            0.0,
            TangentialTransportConfig {
                absolute_energy_tolerance_j: 0.0,
                relative_energy_tolerance: 0.0,
            },
            || false,
        )
        .unwrap();
    assert_eq!(before.state(), &loaded);
    assert_eq!(changed.stored_energy_residual_j, 0.0);
    close2(
        changed.checkpoint.state().principal_axis(),
        rotate([1.0, 0.0], -angle),
        1e-14,
    );
    close2(
        changed.checkpoint.state().slip_m(),
        rotate(loaded.slip_m(), -angle),
        1e-14,
    );
    let original = law.advance(&loaded, &input([0.13, -0.03], 0.01)).unwrap();
    let rebased = law
        .advance(
            changed.checkpoint.state(),
            &input(rotate([0.13, -0.03], -angle), 0.01),
        )
        .unwrap();
    close2(
        rotate(rebased.traction_n, angle),
        original.traction_n,
        1e-10,
    );
    close(rebased.dissipated_j, original.dissipated_j, 1e-12);
}

#[test]
fn invalid_coefficients_foreign_history_and_exhausted_budgets_refuse_atomically() {
    assert!(EllipticCoulomb::new([0.3, 0.2], [0.4, 0.1], [1.0, 0.0]).is_err());
    assert!(EllipticCoulomb::new([0.3, f64::NAN], [0.2, 0.1], [1.0, 0.0]).is_err());
    assert!(EllipticCoulomb::new([0.3, 0.2], [0.2, 0.1], [0.0, 0.0]).is_err());
    let law = law();
    let state = law.initial_state();
    let saved = state.clone();
    let foreign = EllipticCoulomb::new([0.8, 0.4], [0.5, 0.2], [1.0, 0.0]).unwrap();
    assert!(foreign.advance(&state, &input([0.1, 0.1], 0.01)).is_err());
    let mut x = input([0.15, 0.04], 0.01);
    x.solve.max_iterations = 1;
    assert!(law.advance(&state, &x).is_err());
    x = input([0.07, 0.0], 0.01);
    x.solve.maximum_force_n = 6.5;
    assert!(law.advance(&state, &x).is_err()); // feasible static cannot become weaker sliding
    assert_eq!(state, saved);
    let retry = law.advance(&state, &input([0.15, 0.04], 0.01)).unwrap();
    assert_eq!(
        retry,
        law.advance(&saved, &input([0.15, 0.04], 0.01)).unwrap()
    );
}

#[test]
fn explicit_coefficient_continuation_preserves_slip_texture_and_zero_storage() {
    let law = law();
    let before = law
        .advance(&law.initial_state(), &input([0.02, 0.03], 0.0))
        .unwrap()
        .state;
    let (same, unchanged) = law
        .continue_coefficients(&before, [0.8, 0.4], [0.6, 0.2])
        .unwrap();
    assert_eq!(same, law);
    assert_eq!(unchanged, before);
    let (updated, history) = law
        .continue_coefficients(&before, [0.4, 0.8], [0.1, 0.5])
        .unwrap();
    assert_eq!(history.slip_m(), before.slip_m());
    assert_eq!(history.principal_axis(), before.principal_axis());
    assert_eq!(
        history.contact_temperature_k(),
        before.contact_temperature_k()
    );
    assert_eq!(updated.stored_energy_j(&history), 0.0);
    assert!(law.advance(&history, &input([0.1, 0.1], 0.0)).is_err());
    let next = updated
        .advance(&history, &input([0.02, 0.03], 0.0))
        .unwrap();
    let support = (0.02_f64.powi(2) + (5.0_f64 * 0.03).powi(2)).sqrt();
    close2(
        next.traction_n,
        [0.02 / support, 25.0 * 0.03 / support],
        1e-14,
    );
    close(next.dissipated_j, support, 1e-14);
    assert!(
        law.continue_coefficients(&before, [0.1, 0.2], [0.2, 0.1])
            .is_err()
    );
    assert_eq!(before.slip_m(), [0.02, 0.03]);
}
