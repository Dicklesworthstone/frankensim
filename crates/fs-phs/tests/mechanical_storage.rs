//! Original stiffness admission is separate from strict physical heat factors.
use fs_phs::mechanical::LinearMechanicalSystem;

fn identity(n: usize) -> Vec<f64> {
    let mut result = vec![0.0; n * n];
    for i in 0..n {
        result[i * n + i] = 1.0;
    }
    result
}
fn tetrahedron(center: f64, edge: f64) -> Vec<f64> {
    let mut k = vec![0.0; 16];
    k[0] = center;
    for i in 1..4 {
        k[i] = -edge;
        k[i * 4] = -edge;
        k[i * 4 + i] = edge;
    }
    k
}

#[test]
fn original_quadrature_gauge_roundoff_admits_without_changing_stiffness() {
    // Coefficient images from the ORIGINAL quad1d::gauss_legendre,
    // simplex::duffy_quadrature and HeterogeneousHelmholtzFe::new arithmetic:
    // unit-coordinate P1 tetrahedron, rho=.2, orders 3 and 4. Center and
    // off-diagonal sums accumulate separately; their exact represented
    // center grounds are -6.66e-16 and -2.998e-15 despite the positive
    // gradient-square source. Do not replace these with closed-integral
    // constants, rounded-up diagonals or an artificially positive gauge.
    for (center, edge) in [
        (2.499999999999998, 0.8333333333333328),
        (2.4999999999999964, 0.8333333333333331),
    ] {
        let k = tetrahedron(center, edge);
        let mut mass = identity(4);
        // Unrelated kinetic units must not set the elastic admission scale.
        mass[0] = 1e-150;
        mass[5] = 1e150;
        let system = LinearMechanicalSystem::new(&mass, &[0.0; 16], &k, &[], 4, 0, vec![]).unwrap();
        for column in 0..4 {
            let mut q = [0.0; 4];
            q[column] = 1.0;
            let state = system.state(&q, &[0.0; 4], &[]).unwrap();
            let effort = system.system().effort(&state);
            for row in 0..4 {
                assert_eq!(
                    effort[row].to_bits(),
                    k[row * 4 + column].to_bits(),
                    "original K changed at ({row},{column})"
                );
            }
        }
        assert_eq!(system.system().state_dim(), 8);
    }
}

#[test]
fn original_negative_directions_and_zero_diagonal_coupling_cannot_hide_in_other_units() {
    let mass = identity(2);
    for k in [
        [1.0, 0.0, 0.0, -1e-300],
        [1e300, 0.0, 0.0, -f64::from_bits(1)],
        [0.0, 1e-300, 1e-300, 1e300],
        // Positive diagonals but a negative direction resolved after the
        // diagonal congruence: normalized [[1,2],[2,1]]. Raw global K scale
        // would hide the weak negative eigenvalue behind the stiff axis.
        [1e-240, 2.0, 2.0, 1e240],
        [1.0, 1.00000001, 1.00000001, 1.0],
    ] {
        assert!(
            LinearMechanicalSystem::new(&mass, &[0.0; 4], &k, &[], 2, 0, vec![]).is_err(),
            "unresolved physical source accepted: {k:?}"
        );
    }
    // Raw symmetry's global scale must not conceal a material asymmetry
    // between coordinates with very different physical units.
    let asymmetric = [1e240, 1e-6, 1e-5, 1e-240];
    assert!(LinearMechanicalSystem::new(&mass, &[0.0; 4], &asymmetric, &[], 2, 0, vec![]).is_err());
    assert!(LinearMechanicalSystem::new(&[1e-300], &[0.0], &[-1e-300], &[], 1, 0, vec![]).is_err());
    // Finite normalization alone is insufficient: these normalized
    // off-diagonals are 1.7e308 and overflow an unchecked eigensolve. The
    // original K has two resolved negative eigenvalues, masked in aggregate
    // Q by the huge inverse mass. A principal-block bound refuses first.
    let mass = [1e-308, 0.0, 0.0, 0.0, 1e-308, 0.0, 0.0, 0.0, 1e-308];
    let stiffness = [1e-308, 1.7, 1.7, 1.7, 1e-308, 1.7, 1.7, 1.7, 1e-308];
    assert!(LinearMechanicalSystem::new(&mass, &[0.0; 9], &stiffness, &[], 3, 0, vec![]).is_err());
    // The shared storage gate itself must also fail closed when the
    // eigenvalue calculation becomes nonfinite on invalid finite input.
    let overflow = vec![
        1.0, 1.7e308, 1.7e308, 1.7e308, 1.0, 1.7e308, 1.7e308, 1.7e308, 1.0,
    ];
    assert!(fs_phs::QuadraticStorage::new(overflow, 3).is_err());
}

#[test]
fn normalized_storage_admission_does_not_relax_positive_heat_semantics() {
    let matrix = tetrahedron(2.499999999999998, 0.8333333333333328);
    let mass = identity(4);
    assert!(LinearMechanicalSystem::new(&mass, &[0.0; 16], &matrix, &[], 4, 0, vec![]).is_ok());
    // Exactly the same represented negative ground still refuses when it
    // would be interpreted as physical resistance and positive heat.
    assert!(LinearMechanicalSystem::new(&mass, &matrix, &[0.0; 16], &[], 4, 0, vec![]).is_err());
    let weak_loss = [1.0, 0.0, 0.0, -f64::from_bits(1)];
    assert!(
        LinearMechanicalSystem::new(&identity(2), &weak_loss, &[0.0; 4], &[], 2, 0, vec![])
            .is_err()
    );
}

#[test]
fn valid_scaled_storage_and_exact_free_coordinates_are_retained() {
    let mass = identity(2);
    for k in [
        [1e-300, 1.0, 1.0, 1e300],
        [1e300, 1e-300, 1e-300, 1e-300],
        [1e-300, 0.0, 0.0, 0.0],
    ] {
        let system = LinearMechanicalSystem::new(&mass, &[0.0; 4], &k, &[], 2, 0, vec![]).unwrap();
        for column in 0..2 {
            let mut q = [0.0; 2];
            q[column] = 1.0;
            let state = system.state(&q, &[0.0; 2], &[]).unwrap();
            let effort = system.system().effort(&state);
            for row in 0..2 {
                assert_eq!(effort[row], k[row * 2 + column]);
            }
        }
    }
}
