//! Physical velocity jumps on the shared three-dimensional lattice owners.
//!
//! This uses the same hydrodynamic impulse as the D2Q9 and tiled D3Q19
//! owners. Higher kinetic moments are not translated or thermalized. The
//! caller supplies any reciprocal body reaction and the signed kinetic work.

use super::{Grid3, sources::ForcedLattice3};
use crate::impulse::{ImpulseControls, ImpulseError, ImpulseReport, ImpulseSum, poll};
use fs_exec::CancelGate;

mod sealed {
    pub trait NativeLattice {}
    impl NativeLattice for super::super::D3q19 {}
    impl NativeLattice for crate::d3q27::D3q27 {}
}

/// Native quadratures with the rest population at zero and the shared
/// second-order Hermite/Guo moment identities. Sealed because an arbitrary
/// `Lattice3` implementation need not satisfy these fourth-order identities.
pub trait ImpulseLattice3<const Q: usize>: ForcedLattice3<Q> + sealed::NativeLattice {}
impl ImpulseLattice3<19> for super::D3q19 {}
impl ImpulseLattice3<27> for crate::d3q27::D3q27 {}

impl<L: ImpulseLattice3<Q>, const Q: usize> Grid3<L, Q> {
    /// Apply one time-integrated impulse density per original x-fastest cell.
    ///
    /// Density and the force-corrected nonequilibrium second moment survive;
    /// physical velocity changes by `J/rho`. The original uniform force,
    /// collision law and completed-step count are retained. Static solid
    /// entries must be zero. No collision, streaming, time or heat is added.
    /// The tighter of the original grid Mach bound and `maximum_speed` applies.
    ///
    /// Reuses the grid's existing next-population scratch. Refusal or
    /// cancellation never publishes a prefix, and zero impulse retains bits.
    pub fn apply_impulse_density(
        &mut self,
        impulses: &[[f64; 3]],
        controls: ImpulseControls,
        gate: &CancelGate,
    ) -> Result<ImpulseReport<3>, ImpulseError> {
        self.apply_impulse_density_with_forces(impulses, controls, None, gate)
    }

    pub(super) fn apply_impulse_density_with_forces(
        &mut self,
        impulses: &[[f64; 3]],
        mut controls: ImpulseControls,
        forces: Option<&[[f64; 3]]>,
        gate: &CancelGate,
    ) -> Result<ImpulseReport<3>, ImpulseError> {
        poll(gate)?;
        let count = self.f.len();
        controls.admit(count)?;
        if impulses.len() != count || forces.is_some_and(|field| field.len() != count) {
            return Err(ImpulseError::Invalid(
                "impulse and original force fields must cover the grid",
            ));
        }
        controls.maximum_speed = controls
            .maximum_speed
            .min(self.budget.maximum_mach * L::CS2.sqrt());
        let directions = L::E.map(|direction| direction.map(f64::from));
        let mut sum = ImpulseSum::new();
        for (cell, (&original, &impulse)) in self.f.iter().zip(impulses).enumerate() {
            if cell % 64 == 0 {
                poll(gate)?;
            }
            if self.solid[cell] && impulse != [0.0; 3] {
                return Err(ImpulseError::Invalid(
                    "nonzero impulse requires a fully occupied fluid cell",
                ));
            }
            let mut next = original;
            sum.kick(
                &mut next,
                &directions,
                &L::W,
                [0.0; 3],
                forces.map_or(self.force, |field| field[cell]),
                impulse,
                controls,
                cell,
            )?;
            self.next[cell] = next;
        }
        let report = sum.finish(controls)?;
        poll(gate)?;
        if report.affected_cells != 0 {
            core::mem::swap(&mut self.f, &mut self.next);
        }
        Ok(report)
    }
}

impl<L: ImpulseLattice3<Q>, const Q: usize> super::sources::ForcedGrid3<L, Q> {
    /// Instantaneous `J/rho` velocity jump in the accepted spatial force frame.
    /// The force history remains attached to the populations; subsequent
    /// native steps deliver only their ordinary trapezoidal source impulse.
    /// See [`Grid3::apply_impulse_density`] for stress, budgets and atomicity.
    pub fn apply_impulse_density(
        &mut self,
        impulses: &[[f64; 3]],
        controls: ImpulseControls,
        gate: &CancelGate,
    ) -> Result<ImpulseReport<3>, ImpulseError> {
        self.grid
            .apply_impulse_density_with_forces(impulses, controls, Some(&self.forces), gate)
    }
}
