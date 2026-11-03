//! Internal contact must exchange momentum through its original collocation.
//! Energy conservation alone does not detect a force on the center of mass.

use fs_dcontact::{ContactStorage, Obstacle};
use fs_phs::{ConstitutiveDiscreteGradient, Storage, discrete_gradient};

const MASS: [f64; 2] = [2.0, 3.0];
const KINETIC_MATRIX: [f64; 16] = [
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.5,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    1.0 / 3.0,
];

struct Kinetic;

fn kinetic_energy(x: &[f64]) -> f64 {
    0.5 * x[1] * x[1] / MASS[0] + 0.5 * x[3] * x[3] / MASS[1]
}

fn kinetic_gradient(x: &[f64], out: &mut [f64]) {
    out.fill(0.0);
    out[1] = x[1] / MASS[0];
    out[3] = x[3] / MASS[1];
}

impl Storage for Kinetic {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        kinetic_energy(x)
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        kinetic_gradient(x, out);
    }
    fn exact_quadratic_matrix(&self) -> Option<&[f64]> {
        Some(&KINETIC_MATRIX)
    }
}

fn contact(alpha: f64, inner: Box<dyn Storage>) -> ContactStorage {
    ContactStorage::new(
        inner,
        2,
        vec![
            Obstacle::new(
                vec![1.0, -1.0],
                1,
                2,
                vec![0.0],
                vec![1.0],
                6.0,
                alpha,
                "synthetic two-body internal contact".into(),
            )
            .unwrap(),
        ],
    )
    .unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "{actual:e} != {expected:e}"
    );
}

fn work(effort: &[f64], a: &[f64], b: &[f64]) -> f64 {
    effort
        .iter()
        .zip(a)
        .zip(b)
        .map(|((&e, &x), &y)| e * (y - x))
        .sum()
}

#[test]
fn internal_contact_preserves_total_momentum_and_kinetic_effort() {
    for alpha in [1.0, 1.5, 2.0, 3.7] {
        let storage = contact(alpha, Box::new(Kinetic));
        assert!(storage.discrete_gradient_rule().is_some());
        for (left, right) in [(1.0, 1.2), (-0.4, 0.7), (0.7, -0.4), (-0.4, -0.2)] {
            // Unequal displacement increments give a moving center of mass.
            let a = [0.5 * left, 1.0, -0.5 * left, 2.0];
            let b = [0.25 + 0.5 * right, 0.9, 0.25 - 0.5 * right, 2.1];
            let effort = discrete_gradient(&storage, &a, &b);
            assert_eq!(
                effort[0] + effort[2],
                0.0,
                "a closed pair must have zero total internal force"
            );
            assert_eq!(effort[1], f64::midpoint(a[1], b[1]) / MASS[0]);
            assert_eq!(effort[3], f64::midpoint(a[3], b[3]) / MASS[1]);
            close(
                work(&effort, &a, &b),
                storage.hamiltonian(&b) - storage.hamiltonian(&a),
                2e-13,
            );
            let reverse = discrete_gradient(&storage, &b, &a);
            assert_eq!(
                effort, reverse,
                "the contact rule is symmetric in its endpoints"
            );
        }
    }
}

#[test]
fn contact_impulse_is_independent_of_galilean_boost_and_position_origin() {
    let storage = contact(2.0, Box::new(Kinetic));
    let a = [0.5, 1.0, -0.5, 2.0];
    let b = [0.75, 0.9, -0.375, 2.1];
    let dt = 0.125;
    let original = discrete_gradient(&storage, &a, &b);
    for origin in [-128.0, 0.0, 256.0] {
        for boost in [-8.0, 0.0, 16.0] {
            let mut boosted_a = a;
            let mut boosted_b = b;
            for i in 0..2 {
                boosted_a[2 * i] += origin;
                boosted_b[2 * i] += origin + dt * boost;
                boosted_a[2 * i + 1] += MASS[i] * boost;
                boosted_b[2 * i + 1] += MASS[i] * boost;
            }
            let effort = discrete_gradient(&storage, &boosted_a, &boosted_b);
            assert_eq!(dt * effort[0], dt * original[0]);
            assert_eq!(dt * effort[2], dt * original[2]);
            for i in [1, 3] {
                close(effort[i], original[i] + boost, 3e-15);
            }
        }
    }
}

#[test]
fn cubic_contact_matches_independent_polynomial_force_and_endpoint_tangent() {
    let storage = contact(2.0, Box::new(Kinetic));
    let a = [0.5, 1.0, -0.5, 2.0];
    for delta in [0.0, 1e-14, 1e-8, 0.03, 0.7, -0.3] {
        let b = [a[0] + delta, 0.9, a[2], 2.1];
        let left = a[0] - a[2];
        let right = b[0] - b[2];
        let effort = discrete_gradient(&storage, &a, &b);
        // V(p)=2*p^3: exact secant and its derivative, with no energy subtraction.
        let force = 2.0 * (left * left + left * right + right * right);
        close(effort[0], force, 3e-14);
        close(effort[2], -force, 3e-14);
        let direction = [0.3, -0.2, -0.4, 0.7];
        let mut actual = [0.0; 4];
        assert!(
            storage
                .discrete_gradient_rule()
                .unwrap()
                .end_jacobian_action(&a, &b, &direction, &mut actual)
        );
        let response = 2.0 * (left + 2.0 * right) * (direction[0] - direction[2]);
        close(actual[0], response, 3e-14);
        close(actual[2], -response, 3e-14);
        close(actual[1], 0.5 * direction[1] / MASS[0], 1e-15);
        close(actual[3], 0.5 * direction[3] / MASS[1], 1e-15);
    }
}

#[test]
fn endpoint_tangent_matches_the_actual_rule_on_contact_and_crossing_branches() {
    let direction = [0.3, -0.2, -0.4, 0.7];
    for alpha in [1.0, 1.5, 2.5, 8.0] {
        let storage = contact(alpha, Box::new(Kinetic));
        for (left, right) in [
            (0.7, 0.7),
            (0.7, 0.2),
            (-0.3, 0.2),
            (0.3, -0.2),
            (-0.4, -0.2),
        ] {
            let a = [left, 0.4, 0.0, -0.5];
            let b = [right, 0.3, 0.0, -0.4];
            let mut actual = [0.0; 4];
            assert!(
                storage
                    .discrete_gradient_rule()
                    .unwrap()
                    .end_jacobian_action(&a, &b, &direction, &mut actual)
            );
            for h in [1e-5, 5e-6] {
                let plus = std::array::from_fn::<_, 4, _>(|i| b[i] + h * direction[i]);
                let minus = std::array::from_fn::<_, 4, _>(|i| b[i] - h * direction[i]);
                let ep = discrete_gradient(&storage, &a, &plus);
                let em = discrete_gradient(&storage, &a, &minus);
                for i in 0..4 {
                    close(actual[i], (ep[i] - em[i]) / (2.0 * h), 2e-8);
                }
            }
        }
    }
    // At coincident zero penetration choose the declared inactive Newton side.
    let storage = contact(1.0, Box::new(Kinetic));
    let mut tangent = [0.0; 4];
    assert!(
        storage
            .discrete_gradient_rule()
            .unwrap()
            .end_jacobian_action(&[0.0; 4], &[0.0; 4], &[1.0, 0.0, 0.0, 0.0], &mut tangent)
    );
    assert_eq!(tangent, [0.0; 4]);
}

struct OwnedTail;
impl Storage for OwnedTail {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        kinetic_energy(x) + 0.25 * x[4].powi(4)
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        kinetic_gradient(x, out);
        out[4] = x[4].powi(3);
    }
    fn discrete_gradient_rule(&self) -> Option<&dyn ConstitutiveDiscreteGradient> {
        Some(self)
    }
    fn coordinate_scale(&self, coordinate: usize) -> Option<f64> {
        (coordinate == 4).then_some(0.125)
    }
}
impl ConstitutiveDiscreteGradient for OwnedTail {
    fn evaluate(&self, a: &[f64], b: &[f64], midpoint: &[f64], out: &mut [f64]) {
        kinetic_gradient(midpoint, out);
        out[4] = 0.25 * (a[4].powi(3) + a[4].powi(2) * b[4] + a[4] * b[4].powi(2) + b[4].powi(3));
    }
    fn end_jacobian_action(
        &self,
        a: &[f64],
        b: &[f64],
        direction: &[f64],
        out: &mut [f64],
    ) -> bool {
        kinetic_gradient(direction, out);
        out[1] *= 0.5;
        out[3] *= 0.5;
        out[4] = 0.25 * (a[4].powi(2) + 2.0 * a[4] * b[4] + 3.0 * b[4].powi(2)) * direction[4];
        true
    }
}

#[test]
fn distributed_contacts_preserve_the_owned_inner_rule_and_internal_state() {
    let obstacles = vec![
        Obstacle::new(
            vec![1.0, -1.0, -0.5, 0.5, 2.0, -2.0],
            3,
            2,
            vec![0.1, -0.7, 0.4],
            vec![0.3, 0.8, 0.2],
            9.0,
            1.5,
            "synthetic mixed contact orientations and measures".into(),
        )
        .unwrap(),
    ];
    let storage = ContactStorage::new(Box::new(OwnedTail), 2, obstacles).unwrap();
    assert_eq!(storage.coordinate_scale(4), Some(0.125));
    let a = [0.4, 0.5, -0.2, 0.7, 1.0];
    let b = [0.6, 0.4, -0.1, 0.8, 1.2];
    let effort = discrete_gradient(&storage, &a, &b);
    let inner = discrete_gradient(&OwnedTail, &a, &b);
    assert_eq!(effort[0] + effort[2], 0.0);
    for i in [1, 3, 4] {
        assert_eq!(effort[i], inner[i]);
    }
    close(
        work(&effort, &a, &b),
        storage.hamiltonian(&b) - storage.hamiltonian(&a),
        2e-13,
    );
    let direction = [0.3, 0.1, -0.2, -0.4, 0.7];
    let (mut actual, mut expected) = ([0.0; 5], [0.0; 5]);
    assert!(
        storage
            .discrete_gradient_rule()
            .unwrap()
            .end_jacobian_action(&a, &b, &direction, &mut actual)
    );
    assert!(OwnedTail.end_jacobian_action(&a, &b, &direction, &mut expected));
    for i in [1, 3, 4] {
        assert_eq!(actual[i], expected[i]);
    }
    close(actual[0] + actual[2], 0.0, 2e-14);
}

#[test]
fn undeclared_nonlinear_inner_storage_keeps_its_existing_gonzalez_path() {
    struct Undeclared;
    impl Storage for Undeclared {
        fn hamiltonian(&self, x: &[f64]) -> f64 {
            kinetic_energy(x) + 0.25 * x[1].powi(4)
        }
        fn gradient(&self, x: &[f64], out: &mut [f64]) {
            kinetic_gradient(x, out);
            out[1] += x[1].powi(3);
        }
    }
    let storage = contact(2.0, Box::new(Undeclared));
    assert!(storage.discrete_gradient_rule().is_none());
    let a = [0.5, 1.0, -0.5, 2.0];
    let b = [0.8, 0.9, -0.4, 2.1];
    let effort = discrete_gradient(&storage, &a, &b);
    close(
        work(&effort, &a, &b),
        storage.hamiltonian(&b) - storage.hamiltonian(&a),
        2e-13,
    );
}
