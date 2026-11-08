//! Retained contractions for the original spatial material and work account.
//!
//! These are algebraic evaluations, not a second motion or constitutive model.
//! The mass action remains with the native solid operator. Cell storage uses
//! the same fixed C^E and epsilon^S as the material response: in (S,E)
//! coordinates W = S.C^E.S/2 + E.epsilon^S.E/2. The reciprocal e terms still
//! enter both original field equations and cancel in this storage expression.
//! A fixed thermal strain shifts the elastic coordinate to S-S_thermal. Its
//! endpoint difference remains the unshifted total-strain increment.

use fs_math::dd::Dd;

/// Separate leading products from low endpoint words until opposing large
/// coordinate transfers have cancelled. Fixed-size cell contractions allocate
/// nothing and do not replace the owner's original tensors or source work.
#[derive(Default)]
pub(super) struct Products {
    leading: Dd,
    cross: Dd,
    trailing: Dd,
}

impl Products {
    pub(super) fn add(&mut self, left: Dd, right: Dd) {
        self.add_scaled(left, 1.0, right);
    }

    fn add_scaled(&mut self, left: Dd, coefficient: f64, right: Dd) {
        let high = Dd::from_f64(left.hi) * Dd::from_f64(coefficient);
        let low = Dd::from_f64(left.lo) * Dd::from_f64(coefficient);
        self.leading = self.leading + high * Dd::from_f64(right.hi);
        self.cross = self.cross + high * Dd::from_f64(right.lo) + low * Dd::from_f64(right.hi);
        self.trailing = self.trailing + low * Dd::from_f64(right.lo);
    }

    pub(super) fn value(self) -> Dd {
        (self.leading + self.cross) + self.trailing
    }
}

/// Polarization of one original symmetric material block at actual endpoints.
/// Symmetry, positivity and finite-value admission remain with that material.
pub(super) fn quadratic_change<const N: usize>(
    tensor: &[[f64; N]; N],
    before: &[f64; N],
    after: &[f64; N],
) -> Dd {
    let mut products = Products::default();
    for row in 0..N {
        let mean = Dd::from_f64(0.5 * before[row]) + Dd::from_f64(0.5 * after[row]);
        for column in 0..N {
            let delta = Dd::from_f64(after[column]) - Dd::from_f64(before[column]);
            products.add_scaled(mean, tensor[row][column], delta);
        }
    }
    products.value()
}

/// Change of (S-shift).C.(S-shift)/2, retaining the thermal subtraction in the
/// mean and the original endpoint increment separately. Rounding both shifted
/// endpoints first can erase a small physical increment under thermal preload.
pub(super) fn quadratic_change_shifted<const N: usize>(
    tensor: &[[f64; N]; N],
    before: &[f64; N],
    after: &[f64; N],
    shift: &[f64; N],
) -> Dd {
    let mut products = Products::default();
    for row in 0..N {
        let mean = Dd::from_f64(0.5 * before[row]) + Dd::from_f64(0.5 * after[row])
            - Dd::from_f64(shift[row]);
        for column in 0..N {
            let delta = Dd::from_f64(after[column]) - Dd::from_f64(before[column]);
            products.add_scaled(mean, tensor[row][column], delta);
        }
    }
    products.value()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn thermal_reference_keeps_a_strain_increment_below_its_rounded_offset() {
        let shift = 2.0_f64.powi(-15);
        let delta = 2.0_f64.powi(-70);
        // (delta-shift)^2 - shift^2 = -2^-84 + 2^-140 exactly.
        let expected = Dd::from_pair(-2.0_f64.powi(-84), 2.0_f64.powi(-140));
        assert_eq!(
            quadratic_change_shifted(&[[2.0]], &[0.0], &[delta], &[shift]),
            expected
        );
        assert_eq!(
            quadratic_change_shifted(&[[2.0]], &[delta], &[0.0], &[shift]),
            -expected
        );
        assert_eq!(
            quadratic_change(&[[2.0]], &[-shift], &[delta - shift]),
            Dd::ZERO
        );
    }

    #[test]
    fn coupled_material_storage_retains_counteracting_coordinate_changes() {
        let tensor = [[2.0, 1.0], [1.0, 2.0]];
        let datum = 2.0_f64.powi(54);
        let before = [datum, datum];
        let after = [datum + 4.0, datum - 4.0];
        assert_eq!(
            quadratic_change(&tensor, &before, &after),
            Dd::from_f64(16.0)
        );
        assert_eq!(
            quadratic_change(&tensor, &after, &before),
            Dd::from_f64(-16.0)
        );
    }

    #[test]
    fn original_work_keeps_endpoint_difference_and_mean_low_words() {
        let datum = 2.0_f64.powi(55);
        let mut work = Products::default();
        // Two independently supplied channels exchange opposite large amounts.
        // Their exact endpoint work differs by one joule.
        work.add(Dd::ONE, Dd::from_f64(datum) - Dd::from_f64(1.0));
        work.add(-Dd::ONE, Dd::from_f64(datum));
        assert_eq!(work.value(), Dd::from_f64(-1.0));
        let before = [datum, 2.0];
        let after = [1.0, datum];
        let change = quadratic_change(&[[1.0, 0.0], [0.0, 1.0]], &before, &after);
        assert_eq!(change, Dd::from_f64(-1.5));
    }

    #[test]
    fn retained_mechanical_storage_reaches_a_small_heat_balance() {
        let datum = 2.0_f64.powi(54);
        let kinetic = Dd::from_f64(datum) - Dd::ONE;
        let elastic = -Dd::from_f64(datum);
        let supplied = Dd::ZERO;
        let physical_heat = Dd::ONE;
        assert_eq!((kinetic + elastic) + physical_heat - supplied, Dd::ZERO);
        // Rounding each diagnostic before the thermal handoff loses that joule.
        assert_eq!(kinetic.to_f64() + elastic.to_f64() + 1.0, 1.0);
    }
}
