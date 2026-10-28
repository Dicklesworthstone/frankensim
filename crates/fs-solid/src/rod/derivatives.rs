//! Local energy differentiation for the common rod Newton owner.
//!
//! Each element is differenced BEFORE assembly, so energy in unrelated elements
//! cannot erase a small force. The coupling graph is structural, not determined
//! by the current stress: slack links must retain their possible tangent blocks.
//! The body/right SO(3) residual Jacobian is not forcibly symmetrized.

use super::{Rod, SolidError};
use fs_ga::{So3Tangent, Vec3};

/// Internal protocol shared by native and source-bound rod energies.
/// Arbitrary callbacks retain the complete (dense) dependency graph.
pub(super) trait RodEnergy<E> {
    fn value(&self, rod: &Rod) -> Result<f64, E>;

    fn derivative(
        &self, base: &Rod, plus: &Rod, minus: &Rod,
        node: usize, component: usize, step: f64,
    ) -> Result<f64, E> {
        let _ = (base, node, component);
        Ok((self.value(plus)? - self.value(minus)?) / (2.0 * step))
    }

    fn coupled(&self, _row_node: usize, _column_node: usize) -> bool { true }
}

impl<E, F: Fn(&Rod) -> Result<f64, E>> RodEnergy<E> for F {
    fn value(&self, rod: &Rod) -> Result<f64, E> { self(rod) }
}

pub(super) struct UniformEnergy;

fn segment_energy(rod: &Rod, segment: usize) -> Result<f64, SolidError> {
    let (g, k) = rod.strains(segment)?;
    Ok(0.5 * rod.l0 * rod.section.twice_energy_per_length(g, k))
}

impl RodEnergy<SolidError> for UniformEnergy {
    fn value(&self, rod: &Rod) -> Result<f64, SolidError> { rod.energy() }

    fn derivative(
        &self, base: &Rod, plus: &Rod, minus: &Rod,
        node: usize, _component: usize, step: f64,
    ) -> Result<f64, SolidError> {
        let mut force = 0.0;
        for segment in node.saturating_sub(1)..(node + 1).min(base.positions.len() - 1) {
            force += (segment_energy(plus, segment)? - segment_energy(minus, segment)?)
                / (2.0 * step);
        }
        Ok(force)
    }

    fn coupled(&self, row_node: usize, column_node: usize) -> bool {
        row_node.abs_diff(column_node) <= 1
    }
}

/// Change exactly one coordinate. Unrelated frames are not re-canonicalized.
/// Callers restore this node from the base before reusing a probe.
pub(super) fn perturb(
    rod: &mut Rod, node: usize, component: usize, increment: f64,
) -> Result<(), SolidError> {
    if component < 3 {
        let old = rod.positions[node][component];
        let next = old + increment;
        if !next.is_finite() || (increment != 0.0 && next == old) {
            return Err(SolidError::InvalidInput {
                what: "rod translation derivative probe is not representable at this position".into(),
            });
        }
        rod.positions[node][component] = next;
    } else {
        let mut angle = [0.0; 3];
        angle[component - 3] = increment;
        rod.frames[node] = rod.frames[node]
            .body_plus(So3Tangent::new(Vec3::new(angle[0], angle[1], angle[2])))
            .map_err(super::group_error)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rod::{RodSection, TipLoad};
    use std::cell::Cell;

    fn rod(segments: usize) -> Rod {
        Rod::straight(1.0, segments, RodSection { ea: 40.0, ga: 20.0, gj: 3.0, ei: 2.0 })
    }

    #[test]
    fn local_element_residual_matches_full_energy_at_finite_rotation() {
        let mut r = rod(4);
        for node in 1..r.positions.len() {
            r.positions[node][1] = 0.01 * (node * node) as f64;
            perturb(&mut r, node, 5, 0.03 * node as f64).unwrap();
        }
        let local = r.residual(&TipLoad::default(), 1.0, &UniformEnergy, &|| Ok(()), 1.0).unwrap();
        let full = r.residual(&TipLoad::default(), 1.0, &Rod::energy, &|| Ok(()), 1.0).unwrap();
        for (a, b) in local.iter().zip(full) {
            assert!((a - b).abs() < 2e-8 * a.abs().max(1.0), "{a:e} != {b:e}");
        }
    }

    #[test]
    fn unrelated_large_energy_cannot_erase_local_internal_force() {
        let mut r = rod(4);
        r.positions[1][0] = 1e8;
        r.positions[4][1] = 1e-3;
        assert!(r.energy().unwrap() > 1e17);
        let local = r.residual(&TipLoad::default(), 1.0, &UniformEnergy, &|| Ok(()), 1.0).unwrap();
        let expected = r.section.ga * 1e-3 / r.l0;
        assert!((local[19] - expected).abs() < 1e-10);
        let full = r.residual(&TipLoad::default(), 1.0, &Rod::energy, &|| Ok(()), 1.0).unwrap();
        assert_eq!(full[19], 0.0);
    }

    struct Counted<'a> { calls: &'a Cell<usize> }
    impl RodEnergy<SolidError> for Counted<'_> {
        fn value(&self, rod: &Rod) -> Result<f64, SolidError> { rod.energy() }
        fn derivative(&self, base: &Rod, plus: &Rod, minus: &Rod,
            node: usize, component: usize, step: f64) -> Result<f64, SolidError> {
            self.calls.set(self.calls.get() + 1);
            UniformEnergy.derivative(base, plus, minus, node, component, step)
        }
        fn coupled(&self, row: usize, column: usize) -> bool { UniformEnergy.coupled(row, column) }
    }

    #[test]
    fn one_tangent_column_visits_only_incident_nodes() {
        let r = rod(32);
        let calls = Cell::new(0);
        let model = Counted { calls: &calls };
        let coordinates: Vec<_> = (0..r.ndof()).collect();
        let values = r.residual_on(&TipLoad::default(), 0.0, &model, &|| Ok(()),
            1.0, &coordinates, Some(6 * 15)).unwrap();
        assert_eq!(calls.get(), 18);
        for (coordinate, value) in values.iter().enumerate() {
            if (coordinate / 6 + 1).abs_diff(16) > 1 { assert_eq!(*value, 0.0); }
        }
    }

    #[test]
    fn arbitrary_energy_callbacks_do_not_assume_nearest_neighbor_coupling() {
        let r = rod(4);
        let energy = |rod: &Rod| -> Result<f64, SolidError> {
            Ok(0.5 * (rod.positions[4][1] - rod.positions[1][1]).powi(2))
        };
        assert!(energy.coupled(1, 4));
        let mut plus = r.clone();
        let mut minus = r.clone();
        perturb(&mut plus, 4, 1, 3e-4).unwrap();
        perturb(&mut minus, 4, 1, -3e-4).unwrap();
        let free: Vec<_> = (0..r.ndof()).collect();
        let rp = plus.residual_on(&TipLoad::default(), 0.0, &energy, &|| Ok(()),
            1.0, &free, Some(19)).unwrap();
        let rm = minus.residual_on(&TipLoad::default(), 0.0, &energy, &|| Ok(()),
            1.0, &free, Some(19)).unwrap();
        assert!(((rp[1] - rm[1]) / 6e-4 + 1.0).abs() < 1e-10);
    }
}
