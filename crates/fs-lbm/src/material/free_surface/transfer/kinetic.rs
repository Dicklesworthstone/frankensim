//! P1-to-cell kinetic storage with the transport velocity removed before squaring.
//!
//! Each disconnected piece retains one fixed reference velocity U. Source and
//! target independently accumulate M, P_rel = integral rho (v-U), and
//! K_rel = integral rho |v-U|^2/2. Their physical change is
//! delta K = delta K_rel + U dot delta P_rel + |U|^2 delta M/2.
//! Thus an unchanged bulk translation cannot erase unresolved deformation or
//! spin, while actual mass/momentum readback errors keep their physical work.
//! This is storage arithmetic, not a change to the remap, a correction source,
//! or a conversion of numerical projection error into heat.

use fs_math::dd::Dd;

/// Retained original P1 mass, relative momentum and component kinetic storage.
///
/// Values refer to the declared fixed world velocity `reference_velocity_m_s`.
/// They support composition with another physical representation without
/// subtracting rounded absolute energies. This view supplies storage, not a
/// force, heat source, exact-arithmetic certificate or new physical frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LiquidKineticFrame {
    pub reference_velocity_m_s: [f64; 3],
    pub mass_kg: Dd,
    /// Integral rho*(v-U) [kg m/s], in the original world axes.
    pub relative_momentum_kg_m_s: [Dd; 3],
    /// Integral rho*(v_axis-U_axis)^2/2 [J], one value per world axis.
    pub relative_kinetic_energy_j: [Dd; 3],
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct Inventory {
    reference: [f64; 3],
    mass: Dd,
    momentum: [Dd; 3],
    kinetic: [Dd; 3],
}

impl Inventory {
    pub(super) const fn new(reference: [f64; 3]) -> Self {
        Self {
            reference,
            mass: Dd::ZERO,
            momentum: [Dd::ZERO; 3],
            kinetic: [Dd::ZERO; 3],
        }
    }

    pub(super) const fn empty_like(self) -> Self {
        Self::new(self.reference)
    }

    pub(super) fn in_frame(self, reference: [f64; 3]) -> Option<LiquidKineticFrame> {
        if !self.finite() || !reference.into_iter().all(f64::is_finite) {
            return None;
        }
        let mut momentum = self.momentum;
        let mut kinetic = self.kinetic;
        for axis in 0..3 {
            let shift = Dd::from_f64(self.reference[axis]) - Dd::from_f64(reference[axis]);
            momentum[axis] = momentum[axis] + self.mass * shift;
            kinetic[axis] = kinetic[axis]
                + shift * self.momentum[axis]
                + Dd::from_f64(0.5) * self.mass * shift * shift;
        }
        momentum
            .iter()
            .chain(&kinetic)
            .all(|value| value.hi.is_finite() && value.lo.is_finite())
            .then_some(LiquidKineticFrame {
                reference_velocity_m_s: reference,
                mass_kg: self.mass,
                relative_momentum_kg_m_s: momentum,
                relative_kinetic_energy_j: kinetic,
            })
    }

    /// Merge independently integrated pieces at one common velocity reference.
    pub(super) fn add_inventory(&mut self, other: Self) -> Option<()> {
        let other = other.in_frame(self.reference)?;
        self.mass = self.mass + other.mass_kg;
        for axis in 0..3 {
            self.momentum[axis] = self.momentum[axis] + other.relative_momentum_kg_m_s[axis];
            self.kinetic[axis] = self.kinetic[axis] + other.relative_kinetic_energy_j[axis];
        }
        self.finite().then_some(())
    }

    /// Exact P1 simplex mass form: integral N_i N_j = V (1+delta_ij)/20.
    /// `mass` is the original geometry owner's rho*V, without mass lumping.
    pub(super) fn add_tetrahedron(&mut self, mass: f64, velocity: [[f64; 3]; 4]) {
        let mass = Dd::from_f64(mass);
        self.mass = self.mass + mass;
        for axis in 0..3 {
            let mut sum = Dd::ZERO;
            let mut squares = Dd::ZERO;
            for velocity in velocity {
                let relative = Dd::from_f64(velocity[axis]) - Dd::from_f64(self.reference[axis]);
                sum = sum + relative;
                squares = squares + relative * relative;
            }
            self.momentum[axis] = self.momentum[axis] + mass * sum * Dd::from_f64(0.25);
            self.kinetic[axis] =
                self.kinetic[axis] + mass * (sum * sum + squares) / Dd::from_f64(40.0);
        }
    }

    /// Exact cubic P1 density/velocity mass form on the original reference tet.
    /// Density weights are normalized in retained precision; mass is the SAME
    /// independently integrated original inventory. Positive square terms keep
    /// subcell motion beside a large common translation without mass lumping.
    pub(super) fn add_density_tetrahedron(
        &mut self,
        mass: f64,
        density: [f64; 4],
        velocity: [[f64; 3]; 4],
    ) {
        if density.iter().all(|rho| *rho == density[0]) {
            self.add_tetrahedron(mass, velocity);
            return;
        }
        let scale = density.into_iter().fold(0.0_f64, f64::max);
        let weights = density.map(|rho| Dd::from_f64(rho) / Dd::from_f64(scale));
        let weight_sum = weights.into_iter().fold(Dd::ZERO, |a, b| a + b);
        let mass = Dd::from_f64(mass);
        self.mass = self.mass + mass;
        for axis in 0..3 {
            let relative =
                velocity.map(|v| Dd::from_f64(v[axis]) - Dd::from_f64(self.reference[axis]));
            let sum = relative.into_iter().fold(Dd::ZERO, |a, b| a + b);
            let squares = relative.into_iter().fold(Dd::ZERO, |a, b| a + b * b);
            let mut first = Dd::ZERO;
            let mut second = Dd::ZERO;
            for i in 0..4 {
                let shifted_sum = sum + relative[i];
                first = first + weights[i] * shifted_sum;
                second = second
                    + weights[i]
                        * (shifted_sum * shifted_sum + squares + relative[i] * relative[i]);
            }
            self.momentum[axis] =
                self.momentum[axis] + mass * first / (Dd::from_f64(5.0) * weight_sum);
            self.kinetic[axis] =
                self.kinetic[axis] + mass * second / (Dd::from_f64(60.0) * weight_sum);
        }
    }

    /// Actual original liquid mass and force-corrected cell velocity.
    pub(super) fn add_cell(&mut self, mass: f64, velocity: [f64; 3]) {
        let mass = Dd::from_f64(mass);
        self.mass = self.mass + mass;
        for axis in 0..3 {
            let relative = Dd::from_f64(velocity[axis]) - Dd::from_f64(self.reference[axis]);
            self.momentum[axis] = self.momentum[axis] + mass * relative;
            self.kinetic[axis] =
                self.kinetic[axis] + Dd::from_f64(0.5) * mass * relative * relative;
        }
    }

    pub(super) fn change_to(self, target: Self) -> Option<Dd> {
        if self.reference != target.reference || !self.finite() || !target.finite() {
            return None;
        }
        let delta_mass = target.mass - self.mass;
        // Keep axes and the three physical terms separate until each endpoint
        // subtraction. Their common transport inventories never enter the sum.
        let mut relative = Dd::ZERO;
        let mut momentum = Dd::ZERO;
        let mut mass = Dd::ZERO;
        for axis in 0..3 {
            let reference = Dd::from_f64(self.reference[axis]);
            relative = relative + (target.kinetic[axis] - self.kinetic[axis]);
            momentum = momentum + reference * (target.momentum[axis] - self.momentum[axis]);
            mass = mass + Dd::from_f64(0.5) * reference * reference * delta_mass;
        }
        let change = relative + momentum + mass;
        (change.hi.is_finite() && change.lo.is_finite()).then_some(change)
    }

    fn finite(self) -> bool {
        self.reference.into_iter().all(f64::is_finite)
            && Dd::ZERO.lt(self.mass)
            && core::iter::once(self.mass)
                .chain(self.momentum)
                .chain(self.kinetic)
                .all(|value| value.hi.is_finite() && value.lo.is_finite())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn p1_variance_survives_unchanged_bulk_transport() {
        // A mass-40 tetrahedron with nodal transverse values (1,-1,0,0)
        // has zero transverse momentum and K_rel = 2. One cell at the mean
        // loses those 2 J even beside an unchanged 2^54 transport speed.
        let speed = 2.0_f64.powi(54);
        let mut source = Inventory::new([speed, 0.0, 0.0]);
        source.add_tetrahedron(
            40.0,
            [
                [speed, 1.0, 0.0],
                [speed, -1.0, 0.0],
                [speed, 0.0, 0.0],
                [speed, 0.0, 0.0],
            ],
        );
        let mut target = source.empty_like();
        target.add_cell(40.0, [speed, 0.0, 0.0]);
        assert_eq!(source.change_to(target), Some(Dd::from_f64(-2.0)));
        assert_eq!(target.change_to(source), Some(Dd::from_f64(2.0)));
        let total = 20.0 * speed * speed;
        assert_eq!(total + 2.0, total);
    }

    #[test]
    fn actual_mass_and_momentum_changes_are_not_discarded_by_the_reference() {
        let mut source = Inventory::new([4.0, -2.0, 0.0]);
        source.add_cell(2.0, [5.0, -1.0, 0.0]);
        let mut target = source.empty_like();
        target.add_cell(3.0, [6.0, -3.0, 2.0]);
        // Independent endpoint energies: (3/2)*(36+9+4) - (2/2)*(25+1).
        assert_eq!(source.change_to(target), Some(Dd::from_f64(47.5)));
        let mut alternate = Inventory::new([-7.0, 8.0, -3.0]);
        alternate.add_cell(2.0, [5.0, -1.0, 0.0]);
        let mut alternate_target = alternate.empty_like();
        alternate_target.add_cell(3.0, [6.0, -3.0, 2.0]);
        assert_eq!(
            alternate.change_to(alternate_target),
            source.change_to(target)
        );
    }

    #[test]
    fn p1_cross_terms_and_actual_target_spin_are_retained() {
        let mut source = Inventory::new([0.0; 3]);
        source.add_tetrahedron(
            40.0,
            [
                [1.0, 2.0, 0.0],
                [3.0, 0.0, 1.0],
                [-1.0, 2.0, 3.0],
                [1.0, 0.0, -2.0],
            ],
        );
        let mut target = source.empty_like();
        target.add_cell(20.0, [2.0, 1.0, 0.5]);
        target.add_cell(20.0, [0.0, 1.0, 0.5]);
        // P1: ((4^2+4^2+2^2) + (5+10+14+5)) = 70 J.
        // Target: 10*(4+1+1/4) + 10*(0+1+1/4) = 65 J.
        assert_eq!(source.change_to(target), Some(Dd::from_f64(-5.0)));
    }

    #[test]
    fn detached_pieces_retain_opposite_errors_before_shared_budgeting() {
        let reference = [2.0_f64.powi(54), 0.0, 0.0];
        let mut a = Inventory::new(reference);
        a.add_cell(2.0, [reference[0], 1.0, 0.0]);
        let mut b = a.empty_like();
        b.add_cell(2.0, [reference[0], 0.0, 0.0]);
        assert_eq!(a.change_to(b).unwrap() + b.change_to(a).unwrap(), Dd::ZERO);
        let absolute = a.change_to(b).unwrap().abs() + b.change_to(a).unwrap().abs();
        assert_eq!(absolute, Dd::from_f64(2.0));
    }

    #[test]
    fn missing_or_nonfinite_inventories_cannot_supply_a_change() {
        let blank = Inventory::new([0.0; 3]);
        let mut source = blank;
        source.add_cell(1.0, [0.0; 3]);
        assert!(source.change_to(blank).is_none());
        let mut incompatible = Inventory::new([1.0; 3]);
        incompatible.add_cell(1.0, [0.0; 3]);
        assert!(source.change_to(incompatible).is_none());
        let mut nonfinite = blank;
        nonfinite.add_cell(1.0, [f64::INFINITY, 0.0, 0.0]);
        assert!(source.change_to(nonfinite).is_none());
    }

    #[test]
    fn reframed_source_keeps_relative_storage_and_transport_defects() {
        let mut source = Inventory::new([4.0, -2.0, 0.0]);
        source.add_tetrahedron(40.0, [[5.0, -1.0, 0.0]; 4]);
        let at_mean = source.in_frame([5.0, -1.0, 0.0]).unwrap();
        assert_eq!(at_mean.mass_kg, Dd::from_f64(40.0));
        assert_eq!(at_mean.relative_momentum_kg_m_s, [Dd::ZERO; 3]);
        assert_eq!(at_mean.relative_kinetic_energy_j, [Dd::ZERO; 3]);
        let at_rest = source.in_frame([0.0; 3]).unwrap();
        assert_eq!(
            at_rest.relative_momentum_kg_m_s,
            [200.0, -40.0, 0.0].map(Dd::from_f64)
        );
        assert_eq!(
            at_rest.relative_kinetic_energy_j,
            [500.0, 20.0, 0.0].map(Dd::from_f64)
        );
        assert!(source.in_frame([f64::NAN, 0.0, 0.0]).is_none());
    }

    #[test]
    fn merged_components_do_not_lose_their_own_reference_velocities() {
        let mean = [2.0_f64.powi(40), 0.0, 0.0];
        let mut a = Inventory::new([mean[0] - 1.0, 0.0, 0.0]);
        a.add_cell(2.0, [mean[0] - 1.0, 3.0, 0.0]);
        let mut b = Inventory::new([mean[0] + 1.0, 0.0, 0.0]);
        b.add_cell(2.0, [mean[0] + 1.0, -3.0, 0.0]);
        let mut total = Inventory::new(mean);
        total.add_inventory(a).unwrap();
        total.add_inventory(b).unwrap();
        let frame = total.in_frame(mean).unwrap();
        assert_eq!(frame.mass_kg, Dd::from_f64(4.0));
        assert_eq!(frame.relative_momentum_kg_m_s, [Dd::ZERO; 3]);
        assert_eq!(
            frame.relative_kinetic_energy_j,
            [2.0, 18.0, 0.0].map(Dd::from_f64)
        );
    }
    #[test]
    fn graded_density_matches_independent_cubic_barycentric_integrals() {
        let rho = [3.0, 7.0, 14.0, 24.0];
        let velocity = [
            [1.0, -2.0, 0.5],
            [-3.0, 1.0, 2.0],
            [2.0, 4.0, -1.0],
            [0.0, -1.0, 3.0],
        ];
        let mut source = Inventory::new([0.0; 3]);
        source.add_density_tetrahedron(2.0, rho, velocity); // V=1/6, sum(rho)=48.
        let frame = source.in_frame([0.0; 3]).unwrap();
        let factorial = |n: usize| (1..=n).product::<usize>() as f64;
        let integral = |indices: &[usize]| {
            let mut counts = [0usize; 4];
            for &index in indices {
                counts[index] += 1;
            }
            counts.into_iter().map(factorial).product::<f64>() / factorial(indices.len() + 3)
        };
        assert_eq!(frame.mass_kg, Dd::from_f64(2.0));
        for axis in 0..3 {
            let mut momentum = 0.0;
            let mut kinetic = 0.0;
            for a in 0..4 {
                for b in 0..4 {
                    momentum += rho[a] * velocity[b][axis] * integral(&[a, b]);
                    for c in 0..4 {
                        kinetic += 0.5
                            * rho[a]
                            * velocity[b][axis]
                            * velocity[c][axis]
                            * integral(&[a, b, c]);
                    }
                }
            }
            assert!((frame.relative_momentum_kg_m_s[axis].to_f64() - momentum).abs() < 2e-14);
            assert!((frame.relative_kinetic_energy_j[axis].to_f64() - kinetic).abs() < 2e-14);
        }
    }

    #[test]
    fn uniform_density_retains_the_original_p1_inventory_bits() {
        let velocity = [
            [5.0, -2.0, 1.0],
            [-3.0, 0.0, 2.0],
            [4.0, 7.0, -1.0],
            [1.0, -5.0, 3.0],
        ];
        let mut original = Inventory::new([4.0, -2.0, 0.5]);
        let mut weighted = original;
        original.add_tetrahedron(37.0, velocity);
        weighted.add_density_tetrahedron(37.0, [7913.0; 4], velocity);
        assert_eq!(weighted, original);
    }

    #[test]
    fn graded_density_variance_survives_bulk_transport_permutation_and_mass_scaling() {
        let speed = 2.0_f64.powi(54);
        let reference = [speed, 0.0, 0.0];
        let rho = [1.0, 2.0, 3.0, 4.0];
        let velocity = [
            [speed, 1.0, 0.0],
            [speed, -1.0, 0.0],
            [speed, 0.0, 0.0],
            [speed, 0.0, 0.0],
        ];
        let mut original = Inventory::new(reference);
        original.add_density_tetrahedron(40.0, rho, velocity);
        let source = original.in_frame(reference).unwrap();
        let expected_k = 26.0 / 15.0;
        assert!((source.relative_momentum_kg_m_s[1].to_f64() + 0.8).abs() < 1e-15);
        assert!((source.relative_kinetic_energy_j[1].to_f64() - expected_k).abs() < 1e-15);
        let mut target = original.empty_like();
        target.add_cell(40.0, [speed, -0.02, 0.0]);
        let change = original.change_to(target).unwrap().to_f64();
        assert!((change - (0.008 - expected_k)).abs() < 1e-14);
        assert_eq!(20.0 * speed * speed + expected_k, 20.0 * speed * speed);
        let permutation = [3, 1, 0, 2];
        let mut reordered = Inventory::new(reference);
        reordered.add_density_tetrahedron(
            40.0,
            permutation.map(|i| rho[i]),
            permutation.map(|i| velocity[i]),
        );
        assert!((reordered.change_to(target).unwrap().to_f64() - change).abs() < 1e-14);
        let mut scaled = Inventory::new(reference);
        scaled.add_density_tetrahedron(320.0, rho.map(|value| 8.0 * value), velocity);
        let frame = scaled.in_frame(reference).unwrap();
        assert_eq!(frame.mass_kg, Dd::from_f64(8.0) * source.mass_kg);
        for axis in 0..3 {
            assert_eq!(
                frame.relative_momentum_kg_m_s[axis],
                Dd::from_f64(8.0) * source.relative_momentum_kg_m_s[axis]
            );
            assert_eq!(
                frame.relative_kinetic_energy_j[axis],
                Dd::from_f64(8.0) * source.relative_kinetic_energy_j[axis]
            );
        }
    }
}
