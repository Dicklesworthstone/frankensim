//! Instantaneous spring-parameter work at the retained lip coordinates.
//!
//! A tension change does not move the lip or change its momentum. Only its
//! spring energy changes; changing the damping coefficient does no instantaneous
//! work. Keep this actuator transfer separate from pneumatic work and heat.

use fs_phs::PortHamiltonian;

/// Evaluate the spring increment separately from the unchanged kinetic energy.
/// Subtracting the two complete Hamiltonians can erase small but real work.
pub(super) fn spring_work_j(
    previous: &PortHamiltonian,
    next: &PortHamiltonian,
    state: &[f64],
    previous_stiffness: f64,
    next_stiffness: f64,
) -> Result<f64, &'static str> {
    if state.len() != 2
        || state.iter().any(|x| !x.is_finite())
        || ![previous_stiffness, next_stiffness]
            .iter()
            .all(|k| k.is_finite() && *k >= 0.0)
        || !previous.hamiltonian(state).is_finite()
        || !next.hamiltonian(state).is_finite()
    {
        return Err("lip tension requires finite retained and candidate storage");
    }
    let change = next_stiffness - previous_stiffness;
    if change == 0.0 || state[0] == 0.0 {
        return Ok(0.0);
    }
    // sqrt(|delta k|) * q avoids first squaring a very small/large coordinate
    // or halving a subnormal stiffness increment. The result is in sqrt(J).
    let scaled = fs_math::det::sqrt(change.abs()) * state[0];
    let work = ((0.5 * scaled) * scaled).copysign(change);
    if !work.is_finite() {
        return Err("lip tension actuator work is not representable");
    }
    // A long silent decay may leave a coordinate whose energy is below the
    // floating-point range. Preserve that original trajectory and represented
    // zero storage instead of introducing a control refusal at inaudible scale.
    Ok(if work == 0.0 { 0.0 } else { work })
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
struct ModeSum {
    sum: f64,
    correction: f64,
}

impl ModeSum {
    fn add(&mut self, work: f64) -> Result<(), &'static str> {
        // Rendering another chunk or repeating an unchanged control must not
        // change the accumulator's rounding history.
        if work == 0.0 {
            return Ok(());
        }
        let next = self.sum + work;
        let correction = if self.sum.abs() >= work.abs() {
            (self.sum - next) + work
        } else {
            (work - next) + self.sum
        };
        self.correction += correction;
        self.sum = next;
        if !work.is_finite() || !self.value().is_finite() {
            return Err("accumulated lip tension work is not representable");
        }
        Ok(())
    }

    fn value(self) -> f64 {
        self.sum + self.correction
    }
}

/// Two physical modes, accumulated independently in accepted-control order.
/// Pending transfers enter the next successful nonempty block exactly once.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(super) struct TensionWork {
    total: [ModeSum; 2],
    pending: [ModeSum; 2],
}

impl TensionWork {
    /// Stage the complete account before the caller publishes a new island.
    pub(super) fn with_increment(self, work_j: [f64; 2]) -> Result<Self, &'static str> {
        let mut next = self;
        for (i, work) in work_j.into_iter().enumerate() {
            next.total[i].add(work)?;
            next.pending[i].add(work)?;
        }
        Ok(next)
    }

    pub(super) fn total_j(self) -> [f64; 2] {
        self.total.map(ModeSum::value)
    }

    pub(super) fn finish_block(&mut self) -> [f64; 2] {
        let work = self.pending.map(ModeSum::value);
        self.pending = [ModeSum::default(); 2];
        work
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_phs::mass_spring_damper;

    #[test]
    fn work_survives_unchanged_kinetic_storage_and_reverses_sign() {
        let low = mass_spring_damper(1.0, 256.0, 0.0).unwrap();
        let high = mass_spring_damper(1.0, 1280.0, 2.0).unwrap();
        let state = [1.0 / 1024.0, 1e10];
        assert_eq!(low.hamiltonian(&state), high.hamiltonian(&state));
        let expected = 1.0 / 2048.0;
        assert_eq!(
            spring_work_j(&low, &high, &state, 256.0, 1280.0),
            Ok(expected)
        );
        assert_eq!(
            spring_work_j(&high, &low, &state, 1280.0, 256.0),
            Ok(-expected)
        );
        assert_eq!(
            spring_work_j(&low, &high, &[0.0, 1e10], 256.0, 1280.0),
            Ok(0.0)
        );
        assert_eq!(spring_work_j(&low, &low, &state, 256.0, 256.0), Ok(0.0));
    }

    #[test]
    fn modes_and_staged_changes_keep_small_signed_work_and_publish_once() {
        let mut work = TensionWork::default();
        for change in [[1e16, -1e16], [1.0, -2.0], [-1e16, 1e16]] {
            work = work.with_increment(change).unwrap();
        }
        assert_eq!(work.total_j(), [1.0, -2.0]);
        assert_eq!(work.finish_block(), [1.0, -2.0]);
        assert_eq!(work.finish_block(), [0.0, 0.0]);
        assert_eq!(work.total_j(), [1.0, -2.0]);
        assert_eq!(work.with_increment([0.0; 2]).unwrap(), work);
        let before = work;
        assert!(work.with_increment([f64::MAX, f64::INFINITY]).is_err());
        assert_eq!(work, before);
        let large = TensionWork::default()
            .with_increment([f64::MAX, 0.0])
            .unwrap();
        assert!(large.with_increment([f64::MAX, 0.0]).is_err());
    }

    #[test]
    fn nonfinite_storage_refuses_but_subnormal_decay_does_not_block_controls() {
        let low = mass_spring_damper(1.0, 1.0, 0.0).unwrap();
        let high = mass_spring_damper(1.0, 4.0, 0.0).unwrap();
        for state in [[f64::NAN, 0.0], [1e154, 0.0], [0.0, f64::INFINITY]] {
            assert!(spring_work_j(&low, &high, &state, 1.0, 4.0).is_err());
        }
        assert_eq!(
            spring_work_j(&low, &high, &[1e-200, 0.0], 1.0, 4.0),
            Ok(0.0)
        );
    }
}
