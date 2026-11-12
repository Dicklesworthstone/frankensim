//! Fixed symmetric mass forms with independently known small energy changes.

use super::{LinearOp, kinetic_change};

struct Mass([[f64; 2]; 2]);
impl LinearOp for Mass {
    fn n(&self) -> usize {
        2
    }

    fn apply(&self, x: &[f64], y: &mut [f64]) {
        for (row, value) in self.0.iter().zip(y) {
            *value = row[0] * x[0] + row[1] * x[1];
        }
    }

    // Every instance in this file is a symmetric matrix, so its transpose is itself.
    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
        self.apply(x, y);
    }
}

#[test]
fn coupled_mass_retains_countermotion_beside_large_common_velocity() {
    let mass = Mass([[2.0, 1.0], [1.0, 2.0]]);
    let common = 2.0_f64.powi(54);
    let before = [common, common];
    let after = [common + 4.0, common - 4.0];
    // The initial vector is in the common mode; its mass-inner product with
    // [4,-4] vanishes. The counter-mode has eigenvalue one, so delta K=16 J.
    assert_eq!(
        kinetic_change(&mass, &before, &after, &mut || false).unwrap(),
        16.0
    );
    assert_eq!(
        kinetic_change(&mass, &after, &before, &mut || false).unwrap(),
        -16.0
    );
}

#[test]
fn exchanged_motion_retains_small_storage_after_large_transfers_cancel() {
    let mass = Mass([[1.0, 0.0], [0.0, 1.0]]);
    let large = 2.0_f64.powi(54);
    let before = [large, 1.0];
    let after = [2.0, large];
    // The large inventory moves between coordinates unchanged. Only the
    // independent small speeds differ: (2²-1²)/2 = 1.5 J. The second actual
    // difference large-1 needs its low word before the original mass action.
    assert_eq!(
        kinetic_change(&mass, &before, &after, &mut || false).unwrap(),
        1.5
    );
    assert_eq!(
        kinetic_change(&mass, &after, &before, &mut || false).unwrap(),
        -1.5
    );
}
