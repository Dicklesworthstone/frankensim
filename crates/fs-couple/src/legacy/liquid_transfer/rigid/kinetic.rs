//! Original rigid momentum-form storage through a change of representation.
//!
//! A fixed world reference velocity removes common bulk motion before the
//! kinetic comparison. Its mass and momentum terms are restored explicitly,
//! so neither a real source mismatch nor a native readback defect disappears.
//! Principal-frame spin uses the actual admitted moments and momenta. No
//! geometry, spin, dissipation, heating or release work is reconstructed here.

use fs_math::dd::Dd;

#[derive(Clone, Copy, Debug)]
pub(super) struct RigidStorage {
    pub(super) mass: f64,
    pub(super) momentum_world: [f64; 3],
    pub(super) principal_moments: [f64; 3],
    pub(super) spin_principal: [f64; 3],
}

impl RigidStorage {
    fn relative(self, reference: [f64; 3]) -> Option<([Dd; 3], [Dd; 3])> {
        if !self.mass.is_finite()
            || self.mass <= 0.0
            || self
                .principal_moments
                .iter()
                .any(|x| !x.is_finite() || *x <= 0.0)
            || reference
                .iter()
                .chain(&self.momentum_world)
                .chain(&self.spin_principal)
                .any(|x| !x.is_finite())
        {
            return None;
        }
        let mass = Dd::from_f64(self.mass);
        let momentum = core::array::from_fn(|axis| {
            Dd::from_f64(self.momentum_world[axis]) - mass * Dd::from_f64(reference[axis])
        });
        let kinetic = core::array::from_fn(|axis| {
            let spin = Dd::from_f64(self.spin_principal[axis]);
            Dd::from_f64(0.5)
                * (momentum[axis] * momentum[axis] / mass
                    + spin * spin / Dd::from_f64(self.principal_moments[axis]))
        });
        Some((momentum, kinetic))
    }
}

/// Integrated P1 storage minus the original rigid body's independent storage.
/// P1 relative kinetic entries use world axes; principal spin entries need not
/// share those axes because only their scalar sums are compared.
pub(super) fn p1_change(
    rigid: RigidStorage,
    reference: [f64; 3],
    mass: Dd,
    relative_momentum: [Dd; 3],
    relative_kinetic: [Dd; 3],
) -> Option<Dd> {
    let (original_momentum, original_kinetic) = rigid.relative(reference)?;
    if !Dd::ZERO.lt(mass) {
        return None;
    }
    let delta_mass = mass - Dd::from_f64(rigid.mass);
    let mut relative = Dd::ZERO;
    let mut transport = Dd::ZERO;
    for axis in 0..3 {
        let u = Dd::from_f64(reference[axis]);
        relative = relative + (relative_kinetic[axis] - original_kinetic[axis]);
        transport = transport
            + u * (relative_momentum[axis] - original_momentum[axis])
            + Dd::from_f64(0.5) * u * u * delta_mass;
    }
    let change = relative + transport;
    (change.hi.is_finite() && change.lo.is_finite()).then_some(change)
}

/// Actual native child rigid inventories minus their original whole body.
/// The source is not replaced by a fitted inertia or summed rounded energies.
pub(super) fn partition_change(
    original: RigidStorage,
    parts: &[RigidStorage],
    reference: [f64; 3],
) -> Option<Dd> {
    let mut mass = Dd::ZERO;
    let mut momentum = [Dd::ZERO; 3];
    let mut kinetic = [Dd::ZERO; 3];
    for &part in parts {
        let (part_momentum, part_kinetic) = part.relative(reference)?;
        mass = mass + Dd::from_f64(part.mass);
        for axis in 0..3 {
            momentum[axis] = momentum[axis] + part_momentum[axis];
            kinetic[axis] = kinetic[axis] + part_kinetic[axis];
        }
    }
    p1_change(original, reference, mass, momentum, kinetic)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn original_spin_does_not_disappear_beneath_common_translation() {
        let u = 2.0_f64.powi(54);
        let original = RigidStorage {
            mass: 2.0,
            momentum_world: [2.0 * u, 0.0, 0.0],
            principal_moments: [1.0; 3],
            spin_principal: [0.0, 0.0, 2.0],
        };
        assert_eq!(
            p1_change(
                original,
                [u, 0.0, 0.0],
                Dd::from_f64(2.0),
                [Dd::ZERO; 3],
                [Dd::ZERO; 3]
            ),
            Some(Dd::from_f64(-2.0))
        );
        // A P1 field's world-axis velocity energy need not align with the
        // body's principal spin axis. The actual scalar storage is identical.
        assert_eq!(
            p1_change(
                original,
                [u, 0.0, 0.0],
                Dd::from_f64(2.0),
                [Dd::ZERO; 3],
                [2.0, 0.0, 0.0].map(Dd::from_f64)
            ),
            Some(Dd::ZERO)
        );
    }

    #[test]
    fn actual_p1_mass_and_momentum_errors_keep_their_transport_energy() {
        let original = RigidStorage {
            mass: 2.0,
            momentum_world: [10.0, -2.0, 0.0],
            principal_moments: [1.0; 3],
            spin_principal: [0.0; 3],
        };
        // Native P1 values denote mass 3 at velocity (6,-3,2), U=(4,-2,0).
        // Independent K difference is 3/2*(36+9+4) - (25+1) = 47.5 J.
        assert_eq!(
            p1_change(
                original,
                [4.0, -2.0, 0.0],
                Dd::from_f64(3.0),
                [6.0, -3.0, 6.0].map(Dd::from_f64),
                [6.0, 1.5, 6.0].map(Dd::from_f64)
            ),
            Some(Dd::from_f64(47.5))
        );
    }

    #[test]
    fn partition_keeps_spin_and_orbital_energy_under_a_common_boost() {
        let u = 2.0_f64.powi(54);
        // Two unit-mass, unit-inertia children at x=+-1, rotating at omega_z=1.
        // Total I_z=4: 1 J intrinsic spin plus 1 J orbital kinetic storage.
        let original = RigidStorage {
            mass: 2.0,
            momentum_world: [2.0 * u, 0.0, 0.0],
            principal_moments: [2.0, 4.0, 4.0],
            spin_principal: [0.0, 0.0, 4.0],
        };
        let child = |velocity| RigidStorage {
            mass: 1.0,
            momentum_world: [u, velocity, 0.0],
            principal_moments: [1.0; 3],
            spin_principal: [0.0, 0.0, 1.0],
        };
        let mut children = [child(1.0), child(-1.0)];
        assert_eq!(
            partition_change(original, &children, [u, 0.0, 0.0]),
            Some(Dd::ZERO)
        );
        children[1].spin_principal[2] = 2.0;
        assert_eq!(
            partition_change(original, &children, [u, 0.0, 0.0]),
            Some(Dd::from_f64(1.5))
        );
    }

    #[test]
    fn invalid_original_inertia_and_empty_partitions_refuse() {
        let original = RigidStorage {
            mass: 2.0,
            momentum_world: [0.0; 3],
            principal_moments: [1.0; 3],
            spin_principal: [0.0; 3],
        };
        assert!(partition_change(original, &[], [0.0; 3]).is_none());
        let invalid = RigidStorage {
            principal_moments: [0.0; 3],
            ..original
        };
        assert!(partition_change(invalid, &[original], [0.0; 3]).is_none());
        assert!(partition_change(original, &[original], [f64::INFINITY; 3]).is_none());
    }
}
