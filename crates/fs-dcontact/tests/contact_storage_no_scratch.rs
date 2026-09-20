use fs_dcontact::{ContactStorage, Obstacle};
use fs_phs::Storage;

struct Zero;
impl Storage for Zero {
    fn hamiltonian(&self, _: &[f64]) -> f64 { 0.0 }
    fn gradient(&self, _: &[f64], out: &mut [f64]) { out.fill(0.0); }
}

fn storage() -> ContactStorage {
    ContactStorage::new(Box::new(Zero), 2, vec![Obstacle::new(
        vec![1.0, -0.5, -1.0, 2.0, 0.0, 1.0], 3, 2,
        vec![0.001, 0.0, 0.01], vec![1.0, 0.4, 0.8], 2e5, 1.5,
        "synthetic distributed contact derivative fixture".into()).unwrap()]).unwrap()
}

#[test]
fn scalar_contact_evaluation_preserves_the_energy_derivative_and_tail() {
    let s = storage();
    // Two mechanical pairs plus one unrelated memory coordinate.
    for x in [[0.005, 0.1, 0.003, -0.2, 7.0], [-0.004, 0.0, 0.002, 0.0, 9.0]] {
        let mut gradient = [99.0; 5]; s.gradient(&x, &mut gradient);
        for i in [0, 2] {
            let mut plus = x; let mut minus = x;
            plus[i] += 1e-8; minus[i] -= 1e-8;
            let fd = (s.hamiltonian(&plus) - s.hamiltonian(&minus)) / 2e-8;
            assert!((fd - gradient[i]).abs() < 1e-5 * (1.0 + fd.abs()));
        }
        assert_eq!(gradient[1], 0.0); assert_eq!(gradient[3], 0.0); assert_eq!(gradient[4], 0.0);
        assert!(s.hamiltonian(&x) >= 0.0);
    }
}

#[test]
fn scalar_contact_does_not_hide_nan_as_an_open_gap() {
    let s = storage(); let x = [f64::NAN, 0.0, 0.0, 0.0];
    let mut gradient = [0.0; 4]; s.gradient(&x, &mut gradient);
    assert!(s.hamiltonian(&x).is_nan());
    assert!(gradient[0].is_nan());
}
