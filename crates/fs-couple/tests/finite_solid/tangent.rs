//! Optional retained local tangents follow the same native physical trajectory.
use super::*;

#[test]
fn local_midpoint_tangents_match_the_independent_finite_stretch_oracle() {
    with_cx(|_, cx| {
        let card = source_card(6.0, 1.0, false, false);
        let material = bind(&card);
        let fixed: Vec<_> = (0..12).filter(|d| *d != 3).map(|d| (d, 0.0)).collect();
        let forces = [(3, 0.2)];
        let mut u = [[0.0; 3]; 4];
        let mut v = [[0.0; 3]; 4];
        u[1][0] = 0.3;
        v[1][0] = 0.1;
        let initial = FiniteSolidBody::new(
            problem(&material, &fixed, &forces),
            &u,
            &v,
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        assert_eq!(
            initial.tangent_storage(),
            FiniteSolidTangentStorage::Assembled
        );
        let storage = FiniteSolidTangentStorage::LocalHessians {
            maximum_tangent_bytes: 81 * 8,
        };
        let mut local = initial.clone().with_tangent_storage(storage);
        let mut assembled = initial.clone();
        let (mut lambda, mut velocity) = (1.3, 0.1);
        for _ in 0..4 {
            let expected = midpoint_oracle(lambda, velocity, 6.0, 1.0, 0.01, 0.2);
            let (next, report) = local
                .step_candidate(0.01, &forces, [0.0; 3], controls(1e-5), cx)
                .unwrap();
            let (sparse, _) = assembled
                .step_candidate(0.01, &forces, [0.0; 3], controls(1e-5), cx)
                .unwrap();
            near(next.native_state().q[0], expected.0 - 1.0, 2e-9);
            near(next.native_state().v[0], expected.1, 2e-8);
            near(next.native_state().q[0], sparse.native_state().q[0], 2e-10);
            near(report.strain_after_j, strain_energy(expected.0, 1.0), 1e-9);
            assert_eq!(next.tangent_storage(), storage);
            (lambda, velocity) = expected;
            local = next;
            assembled = sparse;
        }
        // A rejected cache allowance cannot advance the original accepted body.
        let bounded =
            initial
                .clone()
                .with_tangent_storage(FiniteSolidTangentStorage::LocalHessians {
                    maximum_tangent_bytes: 81 * 8 - 1,
                });
        let before = bounded.native_state().clone();
        assert!(
            bounded
                .step_candidate(0.01, &forces, [0.0; 3], controls(1e-5), cx)
                .is_err()
        );
        assert_eq!(bounded.native_state(), &before);

        // The separate algorithmic endpoint Jacobian retains its own storage;
        // selecting local ordinary tangents must not change its chain factor.
        let (sparse, _) = initial
            .step_energy_momentum_candidate(0.01, &forces, [0.0; 3], controls(1e-6), cx)
            .unwrap();
        let (same, _) = bounded
            .step_energy_momentum_candidate(0.01, &forces, [0.0; 3], controls(1e-6), cx)
            .unwrap();
        assert_eq!(same.native_state(), sparse.native_state());
    });
}
