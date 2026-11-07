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

#[derive(Clone, Copy, Debug)]
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
}
