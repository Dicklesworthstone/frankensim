//! Physical velocity jumps on the shared three-dimensional lattice owners.
//!
//! This uses the same hydrodynamic impulse as the D2Q9 and tiled D3Q19
//! owners. Higher kinetic moments are not translated or thermalized. The
//! caller supplies any reciprocal body reaction and the signed kinetic work.

use super::{Grid3, sources::ForcedLattice3};
use crate::impulse::{
    ImpulseControls, ImpulseError, ImpulseReport, ImpulseSum, WrenchImpulse3,
    WrenchImpulseControls3, WrenchImpulseReport3, poll,
    wrench::{WrenchCell3, WrenchPlan3},
};
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

    /// Apply a total force-and-couple impulse on selected native fluid cells.
    ///
    /// Indices are strictly increasing and refer to the original x-fast cell
    /// order. Actual population masses determine the least mass-weighted
    /// velocity increment that realizes the requested linear and angular
    /// impulse. Cell centers are (x+1/2,y+1/2,z+1/2), in lattice lengths.
    /// Collinear stencils can carry transverse couples, but an unresolved
    /// axial couple refuses instead of manufacturing a rotational degree of
    /// freedom. Nonperiodic topology is required: this API does not invent an
    /// unwrapped angular-momentum chart on a periodic domain.
    ///
    /// The shared Guo-frame impulse retains density and physical viscous
    /// stress. Both actual momentum and torque are checked before publication,
    /// together with the native speed and signed kinetic-work bounds. There
    /// is no collision, pressure projection, time advance or heat assignment.
    /// A coupled body's owner must receive the opposite wrench.
    ///
    /// Uses existing population scratch, one compact stencil vector, and the
    /// original fixed-size 3-by-3 rank solve. No full-grid impulse array or new
    /// flow iteration/collision is introduced. Any refusal/cancellation preserves all accepted fields.
    pub fn apply_wrench_impulse(
        &mut self,
        cells: &[usize],
        wrench: WrenchImpulse3,
        controls: WrenchImpulseControls3,
        gate: &CancelGate,
    ) -> Result<WrenchImpulseReport3, ImpulseError> {
        self.apply_wrench_impulse_with_forces(cells, wrench, controls, None, gate)
    }

    pub(super) fn apply_impulse_density_with_forces(
        &mut self,
        impulses: &[[f64; 3]],
        controls: ImpulseControls,
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
        self.apply_impulse_with_receipt(
            controls,
            forces,
            gate,
            |cell| impulses[cell],
            |report, _, _| Ok(report),
        )
    }

    pub(super) fn apply_wrench_impulse_with_forces(
        &mut self,
        indices: &[usize],
        wrench: WrenchImpulse3,
        controls: WrenchImpulseControls3,
        forces: Option<&[[f64; 3]]>,
        gate: &CancelGate,
    ) -> Result<WrenchImpulseReport3, ImpulseError> {
        poll(gate)?;
        let count = self.f.len();
        controls.impulse.admit(count)?;
        if self.periodic.into_iter().any(|value| value)
            || indices.is_empty()
            || indices.len() > count
            || forces.is_some_and(|field| field.len() != count)
        {
            return Err(ImpulseError::Invalid(
                "wrench requires nonperiodic topology and a nonempty bounded fluid stencil",
            ));
        }
        let mut cells = Vec::new();
        cells
            .try_reserve_exact(indices.len())
            .map_err(|_| ImpulseError::Invalid("wrench stencil allocation refused"))?;
        let [nx, ny, _] = self.dims;
        for (i, &index) in indices.iter().enumerate() {
            if i % 64 == 0 {
                poll(gate)?;
            }
            if index >= count || self.solid[index] {
                return Err(ImpulseError::Invalid(
                    "wrench stencil contains a nonfluid cell",
                ));
            }
            cells.push(WrenchCell3 {
                index,
                position: [
                    (index % nx) as f64 + 0.5,
                    ((index / nx) % ny) as f64 + 0.5,
                    (index / nx / ny) as f64 + 0.5,
                ],
                mass: self.f[index].iter().sum(),
            });
        }
        let plan = WrenchPlan3::prepare(&cells, wrench, controls, gate)?;
        let uniform_force = self.force;
        let directions = L::E.map(|direction| direction.map(f64::from));
        let mut next_cell = 0;
        self.apply_impulse_with_receipt(
            controls.impulse,
            forces,
            gate,
            |index| {
                if next_cell < cells.len() && cells[next_cell].index == index {
                    let impulse = plan.impulse(cells[next_cell]);
                    next_cell += 1;
                    impulse
                } else {
                    [0.0; 3]
                }
            },
            |report, before, after| {
                plan.finish(&cells, report, controls, gate, |index| {
                    let force = forces.map_or(uniform_force, |field| field[index]);
                    let observe = |populations: &[f64; Q]| {
                        core::array::from_fn::<_, 3, _>(|a| {
                            populations
                                .iter()
                                .zip(&directions)
                                .map(|(f, e)| f * e[a])
                                .sum::<f64>()
                                + 0.5 * force[a]
                        })
                    };
                    let old = observe(&before[index]);
                    let new = observe(&after[index]);
                    core::array::from_fn(|a| new[a] - old[a])
                })
            },
        )
    }

    // Both callers admit shape/resource controls before entering this common
    // population owner. The dependent receipt runs before the original swap.
    fn apply_impulse_with_receipt<R>(
        &mut self,
        mut controls: ImpulseControls,
        forces: Option<&[[f64; 3]]>,
        gate: &CancelGate,
        mut impulse_at: impl FnMut(usize) -> [f64; 3],
        admit: impl FnOnce(ImpulseReport<3>, &[[f64; Q]], &[[f64; Q]]) -> Result<R, ImpulseError>,
    ) -> Result<R, ImpulseError> {
        controls.maximum_speed = controls
            .maximum_speed
            .min(self.budget.maximum_mach * L::CS2.sqrt());
        let directions = L::E.map(|direction| direction.map(f64::from));
        let mut sum = ImpulseSum::new();
        for (cell, &original) in self.f.iter().enumerate() {
            if cell % 64 == 0 {
                poll(gate)?;
            }
            let impulse = impulse_at(cell);
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
        let changed = report.affected_cells != 0;
        let report = admit(report, &self.f, &self.next)?;
        poll(gate)?;
        if changed {
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

    /// The same localized wrench in the retained, spatially varying Guo frame.
    /// No force field is detached or rebound during the instantaneous impulse.
    pub fn apply_wrench_impulse(
        &mut self,
        cells: &[usize],
        wrench: WrenchImpulse3,
        controls: WrenchImpulseControls3,
        gate: &CancelGate,
    ) -> Result<WrenchImpulseReport3, ImpulseError> {
        self.grid.apply_wrench_impulse_with_forces(
            cells,
            wrench,
            controls,
            Some(&self.forces),
            gate,
        )
    }
}

#[cfg(test)]
mod wrench_tests;
