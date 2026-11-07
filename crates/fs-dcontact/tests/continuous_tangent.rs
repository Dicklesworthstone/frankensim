//! Continuous contact derivatives retain the original inner law and its tail.
use fs_dcontact::{ContactStorage, Obstacle};
use fs_phs::{QuadraticStorage, Storage};

fn contact(alpha: f64) -> ContactStorage {
    // Non-unit kinetic rows and a coupled internal coordinate: an assumed
    // unit-mass/position-only inner Hessian would give the wrong result.
    let q = vec![
        2.0,
        0.0,
        0.0,
        0.0,
        0.25,
        0.0,
        0.5,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        3.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        1.0 / 3.0,
        0.0,
        0.25,
        0.0,
        0.0,
        0.0,
        7.0,
    ];
    ContactStorage::new(
        Box::new(QuadraticStorage::new(q, 5).unwrap()),
        2,
        vec![
            Obstacle::new(
                vec![1.0, -1.0],
                1,
                2,
                vec![0.25],
                vec![0.7],
                3.0,
                alpha,
                "manufactured internal contact".into(),
            )
            .unwrap(),
        ],
    )
    .unwrap()
}

fn inner_action(d: [f64; 5]) -> [f64; 5] {
    [
        2.0 * d[0] + 0.25 * d[4],
        0.5 * d[1],
        3.0 * d[2],
        d[3] / 3.0,
        0.25 * d[0] + 7.0 * d[4],
    ]
}

fn close(actual: [f64; 5], expected: [f64; 5]) {
    for (a, e) in actual.into_iter().zip(expected) {
        assert!((a - e).abs() < 2e-13 * (1.0 + e.abs()), "{a:e} != {e:e}");
    }
}

#[test]
fn boxed_continuous_tangent_preserves_contact_kinetics_and_memory() {
    let d = [0.3, -0.2, -0.4, 0.7, 0.12];
    for alpha in [1.0, 1.5, 2.0, 3.7] {
        let storage: Box<dyn Storage> = Box::new(contact(alpha));
        let tangent = storage.continuous_tangent().unwrap();
        for penetration in [-0.4_f64, 0.8] {
            let x = [penetration + 0.25, 0.4, 0.0, -0.6, 0.7];
            let mut expected = inner_action(d);
            if penetration > 0.0 {
                let df = 0.7 * 3.0 * alpha * penetration.powf(alpha - 1.0) * (d[0] - d[2]);
                expected[0] += df;
                expected[2] -= df;
            }
            let mut actual = [f64::NAN; 5];
            assert!(tangent.apply(&x, &d, &mut actual));
            close(actual, expected);
        }
    }
}

#[test]
fn continuous_contact_kink_is_not_the_inactive_newton_representative() {
    let storage = contact(1.0);
    let x = [0.25, 0.4, 0.0, -0.6, 0.7];
    let d = [0.3, -0.2, -0.4, 0.7, 0.12];
    let mut out = [0.0; 5];
    assert!(
        !storage
            .continuous_tangent()
            .unwrap()
            .apply(&x, &d, &mut out)
    );
    // The existing Newton-side helper retains its documented inactive branch.
    assert!(storage.hessian_vector_with(&x, &d, &mut out, |x, d, out| {
        storage
            .inner_storage()
            .continuous_tangent()
            .unwrap()
            .apply(x, d, out)
    }));
    close(out, inner_action(d));
    let tangential = [0.3, -0.2, 0.3, 0.7, 0.12];
    assert!(
        storage
            .continuous_tangent()
            .unwrap()
            .apply(&x, &tangential, &mut out)
    );
    close(out, inner_action(tangential));
    // alpha>1 has zero contact stiffness at entry, including Hertz alpha=1.5.
    for alpha in [1.5, 2.0] {
        assert!(
            contact(alpha)
                .continuous_tangent()
                .unwrap()
                .apply(&x, &d, &mut out)
        );
        close(out, inner_action(d));
    }
}

#[test]
fn contact_does_not_advertise_a_missing_inner_derivative() {
    struct Undeclared;
    impl Storage for Undeclared {
        fn hamiltonian(&self, x: &[f64]) -> f64 {
            0.5 * x.iter().map(|v| v * v).sum::<f64>()
        }
        fn gradient(&self, x: &[f64], out: &mut [f64]) {
            out.copy_from_slice(x);
        }
    }
    let storage = ContactStorage::new(Box::new(Undeclared), 1, vec![]).unwrap();
    assert!(storage.continuous_tangent().is_none());
}
