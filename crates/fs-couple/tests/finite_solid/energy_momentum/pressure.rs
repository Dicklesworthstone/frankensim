//! Uniform held pressure uses an objective algorithmic potential, with its
//! original endpoint-volume work kept separate from material storage.
use super::*;
use fs_solid::hyper3::pressure::ClosedTetPressureSurface;
use std::sync::Arc;

fn surface(cx: &Cx<'_>) -> Arc<ClosedTetPressureSurface> {
    Arc::new(
        ClosedTetPressureSurface::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap(),
    )
}
fn volume(body: &FiniteSolidBody<'_>) -> f64 {
    let u = body.displacement_m();
    let f: [f64; 9] =
        core::array::from_fn(|ij| NODES[ij % 3 + 1][ij / 3] + u[ij % 3 + 1][ij / 3] - u[0][ij / 3]);
    (f[0] * (f[4] * f[8] - f[5] * f[7]) - f[1] * (f[3] * f[8] - f[5] * f[6])
        + f[2] * (f[3] * f[7] - f[4] * f[6]))
        / 6.0
}
fn evolve<'a>(
    initial: &FiniteSolidBody<'a>,
    p: f64,
    h: f64,
    count: usize,
    cx: &Cx<'_>,
) -> FiniteSolidBody<'a> {
    let reference = invariants(initial);
    let initial_total = reference.0 + p * volume(initial);
    let mut body = initial.clone();
    for _ in 0..count {
        let before_volume = volume(&body);
        let (next, report) = body
            .step_energy_momentum_candidate_with_pressure(h, &[], [0.0; 3], p, controls(1e-10), cx)
            .unwrap();
        let measured = invariants(&next);
        let next_volume = volume(&next);
        near(measured.0 + p * next_volume, initial_total, 1e-9);
        near(
            report.pressure_work_j,
            -p * (next_volume - before_volume),
            1e-12,
        );
        near(report.source_work_j, report.pressure_work_j, 1e-14);
        near(report.integration_defect_j, 0.0, 1e-10);
        assert!(report.support_reactions_n.is_empty());
        for i in 0..3 {
            near(measured.1[i], reference.1[i], 1e-10);
            near(measured.2[i], reference.2[i], 1e-10);
        }
        body = next;
    }
    body
}
fn state_distance(a: &FiniteSolidBody<'_>, b: &FiniteSolidBody<'_>) -> f64 {
    a.native_state()
        .q
        .iter()
        .zip(&b.native_state().q)
        .chain(a.native_state().v.iter().zip(&b.native_state().v))
        .map(|(x, y)| (x - y).powi(2))
        .sum::<f64>()
        .sqrt()
}

#[test]
fn pressure_tumbling_closes_work_and_momenta_and_refines_the_actual_motion() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let (u, v) = initial_motion();
        let initial = FiniteSolidBody::new_with_pressure(
            problem(&material, &[], &[]),
            surface(cx),
            3.0,
            &u,
            &v,
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let pressurized = evolve(&initial, 3.0, 0.02, 125, cx);
        let unloaded = evolve(&initial, 0.0, 0.02, 125, cx);
        assert!(state_distance(&pressurized, &unloaded) > 0.1);
        // Energy closure is not a trajectory-accuracy certificate. Resolve the
        // pressure/material dynamics before checking the asymptotic step order.
        let fine = evolve(&initial, 3.0, 0.00125, 160, cx);
        let medium = evolve(&initial, 3.0, 0.005, 40, cx);
        let coarse = evolve(&initial, 3.0, 0.01, 20, cx);
        let ratio = state_distance(&coarse, &fine) / state_distance(&medium, &fine);
        assert!((3.0..6.0).contains(&ratio), "refinement ratio {ratio}");
        assert_eq!(initial.native_state().steps, 0);
    });
}

fn breathing_oracle(s0: f64, v0: f64, h: f64, p: f64) -> (f64, f64) {
    // On F=s I, symmetry and exact potential work give this scalar divided
    // difference of 20*(s^3-1)^2+p*s^3, without production stress or integration.
    let residual = |s1: f64| {
        let stress =
            (20.0 * (s1.powi(3) + s0.powi(3) - 2.0) + p) * (s1 * s1 + s1 * s0 + s0 * s0) / 3.0;
        2.0 * 0.1 * (s1 - s0 - h * v0) / (h * h) + VOLUME * stress
    };
    let (mut low, mut high) = (0.7, 1.4);
    assert!(residual(low) < 0.0 && residual(high) > 0.0);
    for _ in 0..100 {
        let mid = 0.5 * (low + high);
        if residual(mid) > 0.0 {
            high = mid;
        } else {
            low = mid;
        }
    }
    let s1 = 0.5 * (low + high);
    (s1, 2.0 * (s1 - s0) / h - v0)
}

#[test]
fn changing_signed_pressure_matches_independent_breathing_work_and_support_reactions() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let fixed: Vec<_> = (0..12)
            .filter(|dof| ![3, 7, 11].contains(dof))
            .map(|dof| (dof, 0.0))
            .collect();
        let (mut u, mut v) = ([[0.0; 3]; 4], [[0.0; 3]; 4]);
        for i in 0..3 {
            u[i + 1][i] = 0.1;
            v[i + 1][i] = 0.05;
        }
        let mut body = FiniteSolidBody::new_with_pressure(
            problem(&material, &fixed, &[]),
            surface(cx),
            3.0,
            &u,
            &v,
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let (mut stretch, mut velocity) = (1.1, 0.05);
        for p in [3.0, 3.0, -1.5, -1.5] {
            let expected = breathing_oracle(stretch, velocity, 0.02, p);
            let (next, report) = body
                .step_energy_momentum_candidate_with_pressure(
                    0.02,
                    &[],
                    [0.0; 3],
                    p,
                    controls(1e-10),
                    cx,
                )
                .unwrap();
            for i in 0..3 {
                near(next.native_state().q[i], expected.0 - 1.0, 2e-9);
                near(next.native_state().v[i], expected.1, 2e-8);
            }
            near(
                report.pressure_work_j,
                -p * VOLUME * (expected.0.powi(3) - stretch.powi(3)),
                2e-10,
            );
            near(
                report.strain_after_j,
                VOLUME * 20.0 * (expected.0.powi(3) - 1.0).powi(2),
                2e-9,
            );
            near(report.kinetic_after_j, 0.15 * expected.1.powi(2), 2e-9);
            near(report.integration_defect_j, 0.0, 1e-10);
            let reaction = report
                .support_reactions_n
                .iter()
                .find(|(dof, _)| *dof == 0)
                .unwrap()
                .1;
            near(reaction, 0.15 * (expected.1 - velocity) / 0.02, 2e-8);
            (stretch, velocity) = expected;
            body = next;
        }
    });
}

#[test]
fn exact_zero_pressure_and_pressure_refusals_preserve_native_state() {
    with_cx(|gate, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let (u, v) = initial_motion();
        let p = problem(&material, &[], &[]);
        let plain = FiniteSolidBody::new(p, &u, &v, 0.0, Default::default(), cx).unwrap();
        let bound = FiniteSolidBody::new_with_pressure(
            p,
            surface(cx),
            0.0,
            &u,
            &v,
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        assert_eq!(plain.native_state(), bound.native_state());
        let original = plain
            .step_energy_momentum_candidate(0.02, &[], [0.0; 3], controls(1e-10), cx)
            .unwrap();
        let zero = bound
            .step_energy_momentum_candidate_with_pressure(
                0.02,
                &[],
                [0.0; 3],
                0.0,
                controls(1e-10),
                cx,
            )
            .unwrap();
        assert_eq!(original.0.native_state(), zero.0.native_state());
        assert_eq!(
            original.1.integration_defect_j.to_bits(),
            zero.1.integration_defect_j.to_bits()
        );
        assert_eq!(zero.1.pressure_work_j, 0.0);
        let before = bound.native_state().clone();
        assert!(
            plain
                .step_energy_momentum_candidate_with_pressure(
                    0.02,
                    &[],
                    [0.0; 3],
                    3.0,
                    controls(1e-10),
                    cx
                )
                .is_err()
        );
        assert!(
            bound
                .step_energy_momentum_candidate_with_pressure(
                    0.02,
                    &[],
                    [0.0; 3],
                    f64::NAN,
                    controls(1e-10),
                    cx
                )
                .is_err()
        );
        gate.request();
        assert!(
            bound
                .step_energy_momentum_candidate_with_pressure(
                    0.02,
                    &[],
                    [0.0; 3],
                    3.0,
                    controls(1e-10),
                    cx
                )
                .is_err()
        );
        assert_eq!(bound.native_state(), &before);
    });
}
