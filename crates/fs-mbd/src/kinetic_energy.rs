//! Shared kinetic-storage arithmetic for native rigid-body and articulated owners.
//!
//! Inputs are already admitted by the original state/mass boundaries. Keeping
//! the quadratic form here lets single-body, pair and coupled consumers use
//! one allocation-free calculation without depending on their work laws.

use fs_math::dd::Dd;

// Reduce product precision layers before combining. Opposite large transfers
// can cancel while a smaller cross or trailing product remains physical.
#[derive(Clone, Copy, Default)]
pub(super) struct MomentumWork {
    leading: Dd,
    cross: Dd,
    trailing: Dd,
}

impl MomentumWork {
    pub(super) fn add_compensated(&mut self, mean: [Dd; 3], momentum: [Dd; 3]) {
        for (mean, delta) in mean.into_iter().zip(momentum) {
            self.leading = self.leading + Dd::from_f64(mean.hi) * Dd::from_f64(delta.hi);
            self.cross = self.cross
                + Dd::from_f64(mean.hi) * Dd::from_f64(delta.lo)
                + Dd::from_f64(mean.lo) * Dd::from_f64(delta.hi);
            self.trailing = self.trailing + Dd::from_f64(mean.lo) * Dd::from_f64(delta.lo);
        }
    }

    pub(super) fn add_fixed_inertia(
        &mut self,
        mass: f64,
        principal_moments: [f64; 3],
        principal_to_body: [f64; 4],
        // Each endpoint holds [world linear momentum, body angular momentum].
        before: [[f64; 3]; 2],
        after: [[f64; 3]; 2],
    ) {
        let midpoint = |old: [f64; 3], new: [f64; 3]| {
            core::array::from_fn(|i| Dd::from_f64(0.5 * old[i]) + Dd::from_f64(0.5 * new[i]))
        };
        let difference = |old: [f64; 3], new: [f64; 3]| {
            core::array::from_fn(|i| Dd::from_f64(new[i]) - Dd::from_f64(old[i]))
        };
        self.add_compensated(
            midpoint(before[0], after[0]).map(|value| value / Dd::from_f64(mass)),
            difference(before[0], after[0]),
        );
        let mut mean = midpoint(before[1], after[1]);
        let mut delta = difference(before[1], after[1]);
        if principal_to_body != [1.0, 0.0, 0.0, 0.0]
            && !(principal_moments[0] == principal_moments[1]
                && principal_moments[1] == principal_moments[2])
        {
            // Projective normalization evaluates the proper rotation denoted
            // by the stored quaternion; it changes no pose or inertia value.
            let [w, x, y, z] = principal_to_body.map(Dd::from_f64);
            let axis = [-x, -y, -z];
            let twice_inverse_norm = Dd::from_f64(2.0) / (w * w + x * x + y * y + z * z);
            let cross = |a: [Dd; 3], b: [Dd; 3]| {
                [
                    a[1] * b[2] - a[2] * b[1],
                    a[2] * b[0] - a[0] * b[2],
                    a[0] * b[1] - a[1] * b[0],
                ]
            };
            let principal = |value: [Dd; 3]| {
                let t = cross(axis, value).map(|value| value * twice_inverse_norm);
                let correction = cross(axis, t);
                core::array::from_fn(|i| value[i] + w * t[i] + correction[i])
            };
            mean = principal(mean);
            delta = principal(delta);
        }
        self.add_compensated(
            core::array::from_fn(|i| mean[i] / Dd::from_f64(principal_moments[i])),
            delta,
        );
    }

    pub(super) fn total(self) -> Dd {
        self.leading + self.cross + self.trailing
    }

    pub(super) fn difference(self, other: Self) -> Dd {
        (self.leading - other.leading)
            + (self.cross - other.cross)
            + (self.trailing - other.trailing)
    }
}

/// Accumulate actual storage change for fixed mass/COM and changing full inertia.
///
/// Twists use `[angular, linear]` about the same material-frame origin. Their
/// world poses may differ: each isotropic COM-translation norm and each
/// material-frame rotational form can be evaluated in its own endpoint frame.
/// Inputs are admitted by the original articulated model/state boundaries.
/// This measures only storage; no source, parameter or reaction law enters it.
pub(super) fn accumulate_spatial_change(
    change: &mut MomentumWork,
    mass: f64,
    center: [f64; 3],
    inertia_before: [f64; 9],
    before: [f64; 6],
    inertia_after: [f64; 9],
    after: [f64; 6],
) {
    if inertia_before == inertia_after && before == after {
        return;
    }
    let half = Dd::from_f64(0.5);
    let center = center.map(Dd::from_f64);
    let old_angular = [before[0], before[1], before[2]].map(Dd::from_f64);
    let new_angular = [after[0], after[1], after[2]].map(Dd::from_f64);
    let center_velocity = |twist: [f64; 6], angular: [Dd; 3]| -> [Dd; 3] {
        let linear = [twist[3], twist[4], twist[5]].map(Dd::from_f64);
        let offset = [
            angular[1] * center[2] - angular[2] * center[1],
            angular[2] * center[0] - angular[0] * center[2],
            angular[0] * center[1] - angular[1] * center[0],
        ];
        core::array::from_fn(|i| linear[i] + offset[i])
    };
    let old_center = center_velocity(before, old_angular);
    let new_center = center_velocity(after, new_angular);
    let mass = Dd::from_f64(mass);
    change.add_compensated(
        core::array::from_fn(|i| half * old_center[i] + half * new_center[i]),
        core::array::from_fn(|i| mass * (new_center[i] - old_center[i])),
    );

    let old_tensor = inertia_before.map(Dd::from_f64);
    let tensor_action = |tensor: [Dd; 9], vector: [Dd; 3]| {
        core::array::from_fn(|row| {
            tensor[3 * row] * vector[0]
                + tensor[3 * row + 1] * vector[1]
                + tensor[3 * row + 2] * vector[2]
        })
    };
    let mean_angular = core::array::from_fn(|i| half * old_angular[i] + half * new_angular[i]);
    let delta_angular = core::array::from_fn(|i| new_angular[i] - old_angular[i]);
    // Exact variable-metric polarization:
    // ΔKrot = mean(ω)·I0·Δω + 1/2 ω1·(I1-I0)·ω1.
    // Form Δω and ΔI BEFORE products with a dominant unchanged spin.
    change.add_compensated(mean_angular, tensor_action(old_tensor, delta_angular));
    if inertia_before != inertia_after {
        let new_tensor = inertia_after.map(Dd::from_f64);
        change.add_compensated(
            new_angular.map(|value| half * value),
            tensor_action(
                core::array::from_fn(|i| new_tensor[i] - old_tensor[i]),
                new_angular,
            ),
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const IDENTITY: [f64; 4] = [1.0, 0.0, 0.0, 0.0];

    #[test]
    fn opposite_transfers_keep_nondyadic_translation_and_spin_storage() {
        let p = 2.0_f64.powi(54);
        for inertia in [2.0, 3.0, 7.0] {
            let mut sum = MomentumWork::default();
            for impulse in [4.0, -4.0] {
                sum.add_fixed_inertia(
                    inertia,
                    [inertia; 3],
                    IDENTITY,
                    [[p, 0.0, 0.0], [0.0, p, 0.0]],
                    [[p + impulse, 0.0, 0.0], [0.0, p + impulse, 0.0]],
                );
            }
            // Each pair leaves J²/m in translation and J²/I in rotation.
            let expected = 32.0 / inertia;
            assert!((sum.total().to_f64() - expected).abs() <= 8.0 * f64::EPSILON * expected);
        }
    }

    #[test]
    fn rotated_inertia_retains_small_storage_beside_unchanged_spin() {
        let half_angle = core::f64::consts::FRAC_PI_8;
        let axes = [half_angle.cos(), 0.0, 0.0, half_angle.sin()];
        let mut sum = MomentumWork::default();
        sum.add_fixed_inertia(
            3.0,
            [2.0, 3.0, 4.0],
            axes,
            [[0.0; 3], [1.0, -2.0, 1e12]],
            [[0.0; 3], [-2.0, 4.0, 1e12]],
        );
        // Independent inverse tensor: xx=yy=5/12, xy=1/12, zz=1/4.
        assert!((sum.total().to_f64() - 21.0 / 8.0).abs() < 2e-14);
    }

    #[test]
    fn native_velocity_rounding_defect_remains_independent_of_storage() {
        let p = 2.0_f64.powi(54);
        let mut stored = MomentumWork::default();
        stored.add_fixed_inertia(
            3.0,
            [3.0; 3],
            IDENTITY,
            [[p, 0.0, 0.0], [0.0; 3]],
            [[p + 4.0, 0.0, 0.0], [0.0; 3]],
        );
        let mut work = MomentumWork::default();
        // These are the exact represented native COM velocities for m=3.
        let mean = Dd::from_f64(0.5 * 6_004_799_503_160_661.0)
            + Dd::from_f64(0.5 * 6_004_799_503_160_662.0);
        work.add_compensated(
            [mean, Dd::ZERO, Dd::ZERO],
            [Dd::from_f64(4.0), Dd::ZERO, Dd::ZERO],
        );
        assert_eq!(stored.total().to_f64(), 24_019_198_012_642_648.0);
        assert_eq!(stored.total().to_f64(), work.total().to_f64());
        assert_eq!(stored.difference(work).to_f64(), 2.0);
    }

    #[test]
    fn unchanged_storage_and_reversed_endpoints_have_exact_signed_changes() {
        let before = [[4.0, 0.0, 1e12], [2.0, 0.0, 1e12]];
        let after = [[2.0, 0.0, 1e12], [0.0, 0.0, 1e12]];
        for (old, new, expected) in [
            (before, after, -4.0),
            (after, before, 4.0),
            (before, before, 0.0),
        ] {
            let mut sum = MomentumWork::default();
            sum.add_fixed_inertia(2.0, [2.0, 3.0, 4.0], IDENTITY, old, new);
            assert_eq!(sum.total().to_f64(), expected);
        }
    }
    const OLD_TENSOR: [f64; 9] = [2.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 4.0];

    #[test]
    fn dominant_unchanged_spin_cannot_hide_signed_parameter_work() {
        let next_tensor = [4.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 4.0];
        let old = [1.0, 0.0, 1e12, 0.0, 0.0, 1e12];
        let new = [0.5, 0.0, 1e12, 0.0, 0.0, 1e12];
        for (i0, v0, i1, v1, expected) in [
            (OLD_TENSOR, old, next_tensor, new, -0.5),
            (next_tensor, new, OLD_TENSOR, old, 0.5),
            (OLD_TENSOR, old, OLD_TENSOR, old, 0.0),
        ] {
            let mut change = MomentumWork::default();
            accumulate_spatial_change(&mut change, 2.0, [0.0; 3], i0, v0, i1, v1);
            assert_eq!(change.total().to_f64(), expected);
        }
    }

    #[test]
    fn changing_nondiagonal_tensor_keeps_the_original_cross_inertia() {
        let next_tensor = [4.0, 1.0, 0.0, 1.0, 3.0, 0.0, 0.0, 0.0, 4.0];
        let mut change = MomentumWork::default();
        accumulate_spatial_change(
            &mut change,
            2.0,
            [0.0; 3],
            OLD_TENSOR,
            [1.0, -1.0, 1e12, 0.0, 0.0, 1e12],
            next_tensor,
            [9.0 / 11.0, -14.0 / 11.0, 1e12, 0.0, 0.0, 1e12],
        );
        // pxy=(2,-3) and Ixy^-1=[[3,-1],[-1,4]]/11 give
        // Kxy=30/11 J after, from 5/2 J before: gain 5/22 J.
        assert!((change.total().to_f64() - 5.0 / 22.0).abs() < 2e-14);
    }

    #[test]
    fn offset_com_uses_centroidal_motion_and_is_covariant_under_axis_relabeling() {
        for (center, tensor, before, after) in [
            (
                [1.0, 0.0, 0.0],
                OLD_TENSOR,
                [1e12, 0.0, 2.0, 0.0, 3.0, 1e12],
                [1e12, 0.0, 1.0, 0.0, 2.0, 1e12],
            ),
            (
                [0.0, 1.0, 0.0],
                [4.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 3.0],
                [2.0, 1e12, 0.0, 1e12, 0.0, 3.0],
                [1.0, 1e12, 0.0, 1e12, 0.0, 2.0],
            ),
        ] {
            let mut change = MomentumWork::default();
            accumulate_spatial_change(&mut change, 2.0, center, tensor, before, tensor, after);
            // Transverse COM speed changes 5 -> 3 m/s, spin rate 2 -> 1:
            // translation loses 16 J and I=4 rotation loses 6 J.
            assert_eq!(change.total().to_f64(), -22.0);
        }
    }

    #[test]
    fn changing_link_transfers_reduce_before_rounding_absolute_storage() {
        let p = 2.0_f64.powi(54);
        let diagonal = |value| [value, 0.0, 0.0, 0.0, value, 0.0, 0.0, 0.0, value];
        let mut change = MomentumWork::default();
        for (next_inertia, delta_rate) in [(4.0, 4.0), (2.0, -8.0)] {
            accumulate_spatial_change(
                &mut change,
                2.0,
                [0.0; 3],
                diagonal(3.0),
                [p, 0.0, 0.0, 0.0, 0.0, 0.0],
                diagonal(next_inertia),
                [p + delta_rate, 0.0, 0.0, 0.0, 0.0, 0.0],
            );
        }
        // 1/2[4(p+4)^2 + 2(p-8)^2 - 6p^2] = 96 J.
        assert_eq!(change.total().to_f64(), 96.0);
    }
}
