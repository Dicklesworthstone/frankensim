//! Actual rigid impulse transport and fixed-body-inertia storage increments.
//!
//! Fluid, film and authored loads own their impulses and work. This leaf only
//! transports those original spin impulses and measures the represented native
//! endpoints; it never derives a physical transfer from an energy residual.

use fs_ga::{GaError, So3, Vec3};
use fs_math::dd::Dd;

/// Apply the original world COM couples directly to stored body momentum.
/// Recovering them from rounded angular velocities loses resolved components
/// and needlessly applies the full inertia and its inverse a second time.
pub(super) fn kick_spin(
    orientation: So3,
    momentum: Vec3,
    impulses: &[[f64; 3]],
) -> Result<Vec3, GaError> {
    let mut total = [Dd::ZERO; 3];
    for impulse in impulses {
        for i in 0..3 {
            total[i] = total[i] + Dd::from_f64(impulse[i]);
        }
    }
    let [x, y, z] = total.map(Dd::to_f64);
    Ok(momentum + orientation.inverse().rotate(Vec3::new(x, y, z))?)
}

/// Product layers remain separate until opposing component/half-step changes
/// cancel. Storage uses the fixed material frame, even when the pose changes.
#[derive(Clone, Copy, Default)]
pub(super) struct KineticChange {
    leading: Dd,
    cross: Dd,
    trailing: Dd,
}

impl KineticChange {
    pub(super) fn add_translation(&mut self, mass: f64, before: [f64; 3], after: [f64; 3]) {
        self.add(
            mean(before, after),
            difference(before, after).map(|x| Dd::from_f64(mass) * x),
        );
    }

    pub(super) fn add_spin(
        &mut self,
        moments: [f64; 3],
        principal_to_body: [f64; 4],
        before: [f64; 3],
        after: [f64; 3],
    ) {
        let mut mean = mean(before, after);
        let mut change = difference(before, after);
        if principal_to_body != [1.0, 0.0, 0.0, 0.0]
            && !(moments[0] == moments[1] && moments[1] == moments[2])
        {
            // Evaluate the proper rotation represented by the cached frame,
            // including its admitted quaternion norm error, in retained form.
            let [w, x, y, z] = principal_to_body.map(Dd::from_f64);
            let axis = [-x, -y, -z];
            let twice_inverse_norm = Dd::from_f64(2.0) / (w * w + x * x + y * y + z * z);
            let principal = |value: [Dd; 3]| {
                let t = cross(axis, value).map(|x| twice_inverse_norm * x);
                let correction = cross(axis, t);
                core::array::from_fn(|i| value[i] + w * t[i] + correction[i])
            };
            mean = principal(mean);
            change = principal(change);
        }
        self.add(
            core::array::from_fn(|i| mean[i] / Dd::from_f64(moments[i])),
            change,
        );
    }

    fn add(&mut self, mean: [Dd; 3], change: [Dd; 3]) {
        for (m, d) in mean.into_iter().zip(change) {
            self.leading = self.leading + Dd::from_f64(m.hi) * Dd::from_f64(d.hi);
            self.cross = self.cross
                + Dd::from_f64(m.hi) * Dd::from_f64(d.lo)
                + Dd::from_f64(m.lo) * Dd::from_f64(d.hi);
            self.trailing = self.trailing + Dd::from_f64(m.lo) * Dd::from_f64(d.lo);
        }
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

fn mean(before: [f64; 3], after: [f64; 3]) -> [Dd; 3] {
    core::array::from_fn(|i| Dd::from_f64(0.5 * before[i]) + Dd::from_f64(0.5 * after[i]))
}

fn difference(before: [f64; 3], after: [f64; 3]) -> [Dd; 3] {
    core::array::from_fn(|i| Dd::from_f64(after[i]) - Dd::from_f64(before[i]))
}

fn cross(a: [Dd; 3], b: [Dd; 3]) -> [Dd; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_ga::So3Tangent;

    #[test]
    fn direct_couple_preserves_spin_erased_by_world_rate_subtraction() {
        let frame = So3::exp(So3Tangent::new(Vec3::new(
            0.0,
            0.0,
            core::f64::consts::FRAC_PI_4,
        )))
        .unwrap();
        let before = Vec3::new(0.03125, 0.0, 0.0);
        let impulse = 2.0_f64.powi(-62);
        let world = frame.rotate(Vec3::new(0.0, impulse, 0.0)).unwrap();
        let old_rate = frame.rotate(before).unwrap();
        assert_eq!(old_rate + world, old_rate);
        let after = kick_spin(frame, before, &[[world.x, world.y, world.z]]).unwrap();
        assert_eq!(after.x, before.x);
        assert!((after.y / impulse - 1.0).abs() < 1e-15);
        let mut energy = KineticChange::default();
        energy.add_spin(
            [1.0; 3],
            [1.0, 0.0, 0.0, 0.0],
            [before.x, before.y, before.z],
            [after.x, after.y, after.z],
        );
        assert!((energy.total().to_f64() / (0.5 * impulse * impulse) - 1.0).abs() < 2e-15);
    }

    #[test]
    fn transverse_kinetic_storage_survives_large_common_translation() {
        let mass = 2.0_f64.powi(80);
        let before = [2.0_f64.powi(-5), 0.0, 0.0];
        let after = [before[0], 2.0_f64.powi(-32), 0.0];
        let absolute = |v: [f64; 3]| 0.5 * mass * v.into_iter().map(|x| x * x).sum::<f64>();
        assert_eq!(absolute(before), absolute(after));
        let mut energy = KineticChange::default();
        energy.add_translation(mass, before, after);
        assert_eq!(energy.total().to_f64(), 32768.0);
    }

    #[test]
    fn full_principal_frame_retains_spin_work_beside_unchanged_rotation() {
        let before = [1.0, 0.0, 2.0_f64.powi(54)];
        let after = [1.0, 2.0, before[2]];
        let mut forward = KineticChange::default();
        let mut reverse = KineticChange::default();
        let frame = [0.6, 0.0, 0.0, 0.8];
        forward.add_spin([2.0, 4.0, 8.0], frame, before, after);
        reverse.add_spin([2.0, 4.0, 8.0], frame, after, before);
        // R has cos=-7/25, sin=24/25. The independent full inverse
        // tensor gives delta_K=2*(A_yx+A_yy)=2*(-42/625+1201/2500).
        assert!((forward.total().to_f64() - 0.8264).abs() < 1e-15);
        assert_eq!(forward.total(), -reverse.total());
    }

    #[test]
    fn free_rotation_defect_remains_separate_from_couple_work() {
        let background = 2.0_f64.powi(54);
        let before = [background, 1.0, 0.0];
        let midpoint = [background, 1.5, 0.0];
        let kicked = [background, 2.0, 0.0];
        let after = [background, 0.5, 0.0];
        let frame = [1.0, 0.0, 0.0, 0.0];
        let mut total = KineticChange::default();
        total.add_spin([1.0; 3], frame, before, after);
        let mut free = KineticChange::default();
        free.add_spin([1.0; 3], frame, before, midpoint);
        free.add_spin([1.0; 3], frame, kicked, after);
        assert_eq!(total.total().to_f64(), -0.375);
        assert_eq!(free.total().to_f64(), -1.25);
        assert_eq!(total.difference(free).to_f64(), 0.875);
    }
}
