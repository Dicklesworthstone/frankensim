//! Independent P1 mass, energy and momentum checks for nonlinear solid motion.
use super::*;

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn initial_motion() -> ([[f64; 3]; 4], [[f64; 3]; 4]) {
    let f = [1.3, 0.18, -0.08, 0.06, 0.87, 0.12, 0.04, -0.09, 1.08];
    let x: [[f64; 3]; 4] =
        NODES.map(|point| core::array::from_fn(|i| (0..3).map(|j| f[3 * i + j] * point[j]).sum()));
    let center: [f64; 3] = core::array::from_fn(|i| x.iter().map(|p| p[i]).sum::<f64>() / 4.0);
    let u = core::array::from_fn(|a| core::array::from_fn(|i| x[a][i] - NODES[a][i]));
    let v = x.map(|p| {
        let r = core::array::from_fn(|i| p[i] - center[i]);
        let spin = cross([0.9, -0.6, 0.7], r);
        core::array::from_fn(|i| [0.1, -0.08, 0.03][i] + spin[i] + [0.17, -0.11, 0.09][i] * r[i])
    });
    (u, v)
}

// Exact integrals for this one-tet, total-mass-one P1 fixture, independent
// of the production sparse mass, Piola stress, and invariant observations.
fn invariants(body: &FiniteSolidBody<'_>) -> (f64, [f64; 3], [f64; 3]) {
    let u = body.displacement_m();
    let v = body.velocity_m_s();
    let x: [[f64; 3]; 4] =
        core::array::from_fn(|a| core::array::from_fn(|i| NODES[a][i] + u[a][i]));
    let sum_x = core::array::from_fn(|i| x.iter().map(|p| p[i]).sum());
    let sum_v: [f64; 3] = core::array::from_fn(|i| v.iter().map(|p| p[i]).sum());
    let product = cross(sum_x, sum_v);
    let angular = core::array::from_fn(|i| {
        (product[i] + (0..4).map(|a| cross(x[a], v[a])[i]).sum::<f64>()) / 20.0
    });
    let kinetic = (sum_v.iter().map(|a| a * a).sum::<f64>()
        + v.iter().flatten().map(|a| a * a).sum::<f64>())
        / 40.0;
    let f: [f64; 9] = core::array::from_fn(|ij| x[ij % 3 + 1][ij / 3] - x[0][ij / 3]);
    let j = f[0] * (f[4] * f[8] - f[5] * f[7]) - f[1] * (f[3] * f[8] - f[5] * f[6])
        + f[2] * (f[3] * f[7] - f[4] * f[6]);
    // Original source-card Ogden alpha=2 reduces to this closed expression.
    let potential = VOLUME
        * (7.5 * (f.iter().map(|a| a * a).sum::<f64>() / j.powf(2.0 / 3.0) - 3.0)
            + 20.0 * (j - 1.0).powi(2));
    (kinetic + potential, sum_v.map(|a| a / 4.0), angular)
}

#[test]
fn freely_tumbling_finite_strain_preserves_energy_and_both_momenta() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let (u, v) = initial_motion();
        let initial = FiniteSolidBody::new(
            problem(&material, &[], &[]),
            &u,
            &v,
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let reference = invariants(&initial);
        let mut ordinary = initial.clone();
        let mut compatible = initial.clone();
        let mut ordinary_drift = 0.0_f64;
        let mut compatible_drift = 0.0_f64;
        for _ in 0..125 {
            ordinary = ordinary
                .step_candidate(0.02, &[], [0.0; 3], controls(1e-2), cx)
                .unwrap()
                .0;
            let (next, report) = compatible
                .step_energy_momentum_candidate(0.02, &[], [0.0; 3], controls(1e-10), cx)
                .unwrap();
            near(report.integration_defect_j, 0.0, 1e-10);
            assert!(report.support_reactions_n.is_empty());
            compatible = next;
            let measured = invariants(&compatible);
            ordinary_drift = ordinary_drift.max((invariants(&ordinary).0 - reference.0).abs());
            compatible_drift = compatible_drift.max((measured.0 - reference.0).abs());
            for i in 0..3 {
                near(measured.1[i], reference.1[i], 1e-10);
                near(measured.2[i], reference.2[i], 1e-10);
            }
        }
        assert!(
            ordinary_drift > 1e-6,
            "fixture must exercise nonlinear midpoint's energy defect"
        );
        assert!(compatible_drift < 1e-9, "energy drift {compatible_drift:e}");
        assert!(compatible_drift < 1e-4 * ordinary_drift);
        assert_eq!(initial.native_state().steps, 0);
        assert_eq!(compatible.native_state().steps, 125);
    });
}

#[test]
fn fixed_supports_and_held_loads_keep_the_same_potential_work_identity() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let fixed: Vec<_> = (0..12)
            .filter(|dof| *dof != 3)
            .map(|dof| (dof, 0.0))
            .collect();
        let force = [(3, 0.2)];
        let mut u = [[0.0; 3]; 4];
        u[1][0] = 0.3;
        let mut v = [[0.0; 3]; 4];
        v[1][0] = 0.1;
        let mut body = FiniteSolidBody::new(
            problem(&material, &fixed, &force),
            &u,
            &v,
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        for _ in 0..8 {
            let before = body.native_state().clone();
            let (next, report) = body
                .step_energy_momentum_candidate(0.02, &force, [0.0; 3], controls(1e-10), cx)
                .unwrap();
            let work = 0.2 * (next.native_state().q[0] - before.q[0]);
            let energy_change = strain_energy(1.0 + next.native_state().q[0], 1.0)
                - strain_energy(1.0 + before.q[0], 1.0)
                + 0.05 * (next.native_state().v[0].powi(2) - before.v[0].powi(2));
            near(report.source_work_j, work, 1e-13);
            near(energy_change, work, 1e-10);
            near(report.integration_defect_j, 0.0, 1e-10);
            assert_eq!(report.support_reactions_n.len(), 11);
            assert_eq!(body.native_state(), &before);
            body = next;
        }
    });
}

#[test]
fn refused_or_cancelled_algorithmic_steps_leave_the_accepted_state_intact() {
    with_cx(|gate, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let (u, v) = initial_motion();
        let body = FiniteSolidBody::new(
            problem(&material, &[], &[]),
            &u,
            &v,
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let before = body.native_state().clone();
        let mut stalled = controls(1e-10);
        stalled.solve.max_newton_iterations = 0;
        assert!(
            body.step_energy_momentum_candidate(0.02, &[], [0.0; 3], stalled, cx)
                .is_err()
        );
        assert_eq!(body.native_state(), &before);
        gate.request();
        assert!(
            body.step_energy_momentum_candidate(0.02, &[], [0.0; 3], controls(1e-10), cx)
                .is_err()
        );
        assert_eq!(body.native_state(), &before);
    });
}
