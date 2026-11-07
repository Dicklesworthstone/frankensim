use super::*;
use fs_dcontact::ContactStorage;
use fs_phs::Storage;

struct Zero;
impl Storage for Zero {
    fn hamiltonian(&self, _: &[f64]) -> f64 {
        0.0
    }
    fn gradient(&self, _: &[f64], out: &mut [f64]) {
        out.fill(0.0);
    }
}
fn close(a: f64, b: f64, coordinate_roundoff: f64) {
    assert!(
        (a - b).abs() <= 3e-12 * a.abs().max(b.abs()) + coordinate_roundoff,
        "{a:e} != {b:e}"
    );
}
fn law(gap: f64, alpha: f64) -> Obstacle {
    Obstacle::new(
        vec![-1.0],
        1,
        1,
        vec![gap],
        vec![0.7],
        5000.0,
        alpha,
        "original synthetic affine point".into(),
    )
    .unwrap()
    .with_internal_loss(0.2)
    .unwrap()
}

#[test]
fn fixed_affine_displacement_preserves_native_potential_secant_and_physical_loss() {
    for alpha in [1.0, 1.5, 2.0] {
        for offset in [-3e-4, 0.0, 3e-4] {
            let original = law(1e-4, alpha);
            let lowered = lower(&original, offset).unwrap();
            let physical = ContactStorage::new(Box::new(Zero), 1, vec![original.clone()]).unwrap();
            let effective = ContactStorage::new(Box::new(Zero), 1, vec![lowered.clone()]).unwrap();
            assert_eq!(lowered.collocation(), original.collocation());
            assert_eq!(lowered.weights(), original.weights());
            assert_eq!(lowered.stiffness(), original.stiffness());
            assert_eq!(lowered.alpha(), original.alpha());
            assert_eq!(lowered.internal_loss(), original.internal_loss());
            assert_eq!(lowered.provenance(), original.provenance());
            for (before, after) in [(-2e-4, 6e-4), (6e-4, -2e-4), (4e-4, 4.1e-4), (4e-4, 4e-4)] {
                let a = OpeningContactStep::new(&original, -(before + offset))
                    .unwrap()
                    .coefficients(-(after + offset))
                    .unwrap();
                let b = OpeningContactStep::new(&lowered, -before)
                    .unwrap()
                    .coefficients(-after)
                    .unwrap();
                // The two algebraic coordinate substitutions round differently
                // exactly at touching. Bound that displacement roundoff through
                // the original force slope; do not compare tiny rounded contact
                // artifacts by their own relative magnitude.
                let dx = 8.0
                    * f64::EPSILON
                    * (before.abs() + after.abs() + offset.abs() + original.gaps()[0].abs());
                let penetration = (before.max(after) + offset - original.gaps()[0]).max(0.0) + dx;
                let force_roundoff = original.weights()[0]
                    * original.stiffness()
                    * alpha
                    * penetration.powf(alpha - 1.0)
                    * dx;
                close(a.0, b.0, force_roundoff);
                close(a.1, b.1, original.internal_loss() * force_roundoff);
                for x in [before, after] {
                    close(
                        physical.hamiltonian(&[-(x + offset), 0.0]),
                        effective.hamiltonian(&[-x, 0.0]),
                        original.weights()[0] * original.stiffness() * penetration.powf(alpha) * dx,
                    );
                }
                let change =
                    effective.hamiltonian(&[-after, 0.0]) - effective.hamiltonian(&[-before, 0.0]);
                let dt = 1e-3;
                let delta = before - after;
                let velocity = delta / dt;
                let loss_force = (b.1 * velocity).min(b.0);
                let loss = loss_force * delta;
                let work = (b.0 - loss_force) * delta;
                assert!(loss >= 0.0);
                assert!(
                    (change + work + loss).abs()
                        <= 3e-12 * (change.abs() + work.abs() + loss).max(1e-25)
                );
            }
        }
    }
    // Positive offset closes a zero-coordinate contact. Negative effective
    // clearance retains its initial potential; it is not clamped away.
    let original = law(1e-4, 1.0);
    let lowered = lower(&original, 3e-4).unwrap();
    assert!(lowered.gaps()[0] < 0.0);
    let storage = ContactStorage::new(Box::new(Zero), 1, vec![lowered]).unwrap();
    close(
        storage.hamiltonian(&[0.0, 0.0]),
        0.5 * 0.7 * 5000.0 * (2e-4_f64).powi(2),
        0.0,
    );
}

#[test]
fn zero_shift_retains_the_original_law_and_clearance_bits() {
    for gap in [0.0, -0.0, 1e-4, -1e-4] {
        for offset in [0.0, -0.0] {
            let original = law(gap, 1.5);
            let lowered = lower(&original, offset).unwrap();
            assert_eq!(lowered.gaps()[0].to_bits(), original.gaps()[0].to_bits());
            assert_eq!(lowered.provenance(), original.provenance());
            assert_eq!(
                lowered.internal_loss().to_bits(),
                original.internal_loss().to_bits()
            );
        }
    }
}

#[test]
fn malformed_point_laws_and_unrepresentable_fixed_clearance_refuse() {
    let original = law(1e-4, 1.5);
    for offset in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(lower(&original, offset).is_err());
    }
    assert!(lower(&law(f64::MAX, 1.0), -f64::MAX).is_err());
    let raw = Obstacle::from_raw_parts(
        vec![-1.0],
        1,
        vec![],
        vec![1.0],
        1.0,
        1.0,
        "invalid shape".into(),
    );
    assert!(lower(&raw, 1e-4).is_err());
    let distributed = Obstacle::new(
        vec![-1.0; 2],
        2,
        1,
        vec![0.0; 2],
        vec![0.5; 2],
        1.0,
        1.0,
        "two points are not one stop".into(),
    )
    .unwrap();
    assert!(lower(&distributed, 1e-4).is_err());
}
