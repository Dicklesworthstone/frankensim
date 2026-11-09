//! Spatial force epochs on the original native D3Q19/D3Q27 open grid.
//!
//! Accepted momentum is `sum(f e) + F_n/2`. Collision applies `F_n`; the
//! streamed donor, reservoir target and final velocity use `F_(n+1)` at their
//! OWN cells. Thus the body-force impulse is `(sum F_n + sum F_(n+1))/2`, in
//! addition to the actual wall/reservoir transport. This is the existing
//! native source-epoch convention, with no extra per-step force reframe.
//!
//! The original uniform-force owner is private after attachment: exposing it
//! alone would lose the field required to interpret its stored populations.
//! Original collision choices, pressure profiles and open/wall link ownership
//! remain unchanged. Force and population epochs publish together only after
//! the complete native candidate and reservoir account have passed admission.

use super::{
    BoundarySpec3, CharacteristicOutlet3, EquilibriumLattice3, Face3, GridError3,
    OpenBoundaryGrid3, OpenBoundaryMomentumTransport3, initial, macro_cell,
};
use crate::d3q19::{
    BoundaryProfileError3, InletFlux3, PressureProfileReport3, VelocityProfileReport3,
};
use crate::lattice3::sources::admit_forces;
use fs_exec::{CancelGate, Cx};

mod acceleration;
/// Held implicit porous resistance through the same native open faces.
pub mod resistance;

/// An open grid whose accepted force field is inseparable from its populations.
/// Force values are densities in lattice units, not accelerations. Solid grid
/// entries must be zero; all fields retain the original x-fast cell order.
pub struct ForcedOpenBoundaryGrid3<L: EquilibriumLattice3<Q>, const Q: usize> {
    pub(super) inner: OpenBoundaryGrid3<L, Q>,
    pub(super) forces: Vec<[f64; 3]>,
    // Reused only by density-dependent endpoints. The authored acceleration
    // stays immutable while neighboring reservoir donors still need it.
    pending_forces: Vec<[f64; 3]>,
}

/// Source-coordinate change at a fixed physical velocity. The native source
/// owns the quadrature-specific second and higher moments; replacing a state
/// by equilibrium would destroy its nonequilibrium stress.
fn reframe_cell<L: EquilibriumLattice3<Q>, const Q: usize>(
    populations: [f64; Q],
    collision: L::Collision,
    previous: [f64; 3],
    next: [f64; 3],
    cell: usize,
    maximum_mach: f64,
) -> Result<[f64; Q], GridError3<L::Error>> {
    let (_, velocity) = macro_cell::<L, Q>(&populations, previous, cell, maximum_mach)?;
    initial::admit_velocity::<L, Q>(collision, velocity, cell, maximum_mach)?;
    if previous == next {
        return Ok(populations);
    }
    let delta: [f64; 3] = core::array::from_fn(|axis| next[axis] - previous[axis]);
    if delta.iter().any(|v| !v.is_finite()) {
        return Err(GridError3::MaterialField {
            cell,
            reason: "force-frame difference is not representable",
        });
    }
    let mut candidate = populations;
    initial::attach_force::<L, Q>(&mut candidate, collision, velocity, delta);
    let (_, represented) = macro_cell::<L, Q>(&candidate, next, cell, maximum_mach)?;
    initial::admit_velocity::<L, Q>(collision, represented, cell, maximum_mach)?;
    Ok(candidate)
}

impl<L: EquilibriumLattice3<Q>, const Q: usize> ForcedOpenBoundaryGrid3<L, Q> {
    /// Attach an explicit spatial force field without an instantaneous impulse.
    ///
    /// Density, physical velocity and force-corrected nonequilibrium stress are
    /// preserved by changing the stored half-source frame with the native
    /// source kernel. The input uniform owner may already have evolved. Its
    /// accepted time, boundary targets and pressure profiles are retained.
    /// Identical force cells and solids retain their population bits.
    ///
    /// The original background collision selects the initialization source
    /// family; subsequent optional collision fields keep their native meaning.
    /// The supplied force vector becomes the sole accepted source owner and
    /// adds no adapter allocation beyond the existing population scratch.
    pub fn try_new_preserving_velocity(
        inner: OpenBoundaryGrid3<L, Q>,
        cx: &Cx<'_>,
        forces: Vec<[f64; 3]>,
    ) -> Result<Self, GridError3<L::Error>> {
        Self::try_new_preserving_velocity_with_poll(inner, forces, &mut || {
            cx.checkpoint().map_err(|_| GridError3::Cancelled)
        })
    }

    pub(super) fn try_new_preserving_velocity_with_poll(
        mut inner: OpenBoundaryGrid3<L, Q>,
        forces: Vec<[f64; 3]>,
        poll: &mut impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<Self, GridError3<L::Error>> {
        admit_forces(&inner.grid, &forces, None, poll)?;
        let grid = &mut inner.grid;
        let mut changed = false;
        for (cell, &force) in forces.iter().enumerate() {
            if cell % 64 == 0 {
                poll()?;
            }
            grid.next[cell] = if grid.solid[cell] {
                grid.f[cell]
            } else {
                changed |= force != grid.force;
                reframe_cell::<L, Q>(
                    grid.f[cell],
                    grid.collision,
                    grid.force,
                    force,
                    cell,
                    grid.budget.maximum_mach,
                )?
            };
        }
        poll()?;
        if changed {
            core::mem::swap(&mut grid.f, &mut grid.next);
        }
        Ok(Self {
            inner,
            forces,
            pending_forces: Vec::new(),
        })
    }

    pub fn dimensions(&self) -> [usize; 3] {
        self.inner.dimensions()
    }
    /// Accepted populations in the same epoch as [`Self::forces`].
    pub fn populations(&self) -> &[[f64; Q]] {
        self.inner.populations()
    }
    pub fn completed_steps(&self) -> usize {
        self.inner.completed_steps()
    }
    pub fn boundaries(&self) -> BoundarySpec3 {
        self.inner.boundaries()
    }
    /// Actual accepted force density at each cell. There is no uniform-force
    /// getter that could silently omit the spatial field.
    pub fn forces(&self) -> &[[f64; 3]] {
        &self.forces
    }
    /// Physical velocity includes this cell's accepted half-force correction.
    pub fn macroscopic(&self, cell: usize) -> Result<(f64, [f64; 3]), GridError3<L::Error>> {
        let populations = self
            .inner
            .grid
            .f
            .get(cell)
            .ok_or(GridError3::Input("cell index"))?;
        macro_cell::<L, Q>(
            populations,
            self.forces[cell],
            cell,
            self.inner.grid.budget.maximum_mach,
        )
    }

    /// Update original reservoir targets without changing topology or time.
    pub fn rebind_open_boundaries(
        &mut self,
        boundaries: BoundarySpec3,
        gate: &CancelGate,
    ) -> Result<(), GridError3<L::Error>> {
        self.inner.rebind_open_boundaries(boundaries, gate)
    }

    /// Retain the original spatial-pressure API and its atomic source queries.
    pub fn set_pressure_profile(
        &mut self,
        face: Face3,
        gate: &CancelGate,
        sample: impl FnMut([usize; 3]) -> Result<f64, BoundaryProfileError3>,
    ) -> Result<PressureProfileReport3, BoundaryProfileError3> {
        self.inner.set_pressure_profile(face, gate, sample)
    }
    pub fn clear_pressure_profile(&mut self, face: Face3) {
        self.inner.clear_pressure_profile(face);
    }
    pub fn pressure_profile(&self, face: Face3) -> Option<&[f64]> {
        self.inner.pressure_profile(face)
    }

    /// Attach the native acoustic outlet without rebinding accepted forces.
    pub fn set_characteristic_outlet(
        &mut self,
        face: Face3,
        control: CharacteristicOutlet3,
        gate: &CancelGate,
    ) -> Result<(), GridError3<L::Error>> {
        self.inner.set_characteristic_outlet(face, control, gate)
    }
    pub fn clear_characteristic_outlet(&mut self, face: Face3) {
        self.inner.clear_characteristic_outlet(face);
    }
    pub fn characteristic_outlet(&self, face: Face3) -> Option<CharacteristicOutlet3> {
        self.inner.characteristic_outlet(face)
    }

    /// Hold a resolved velocity trace on the same force-aware native owner.
    /// Sampling and optional normal-volume control do not advance the force
    /// or fluid epochs; reconstruction uses each endpoint's actual force.
    pub fn set_velocity_profile(
        &mut self,
        face: Face3,
        flux: InletFlux3,
        gate: &CancelGate,
        sample: impl FnMut([usize; 3]) -> Result<[f64; 3], BoundaryProfileError3>,
    ) -> Result<VelocityProfileReport3, BoundaryProfileError3> {
        self.inner.set_velocity_profile(face, flux, gate, sample)
    }

    /// Also enforce an inclusive speed bound on the corrected final trace.
    pub fn set_velocity_profile_with_speed_limit(
        &mut self,
        face: Face3,
        flux: InletFlux3,
        maximum_speed_lu: f64,
        gate: &CancelGate,
        sample: impl FnMut([usize; 3]) -> Result<[f64; 3], BoundaryProfileError3>,
    ) -> Result<VelocityProfileReport3, BoundaryProfileError3> {
        self.inner
            .set_velocity_profile_with_speed_limit(face, flux, maximum_speed_lu, gate, sample)
    }

    /// Restore the held uniform fallback without changing physical state.
    pub fn clear_velocity_profile(&mut self, face: Face3) {
        self.inner.clear_velocity_profile(face);
    }

    /// The accepted vector plane, independent of the force-history buffers.
    #[must_use]
    pub fn velocity_profile(&self, face: Face3) -> Option<&[[f64; 3]]> {
        self.inner.velocity_profile(face)
    }

    /// Advance with an authored force field at the NEXT accepted time.
    ///
    /// Collision uses the retained opening field; streamed-state admission
    /// and reservoir reconstruction use the endpoint field. No artificial
    /// half kick is added by swapping source arrays. The returned receipt is
    /// actual reservoir transport; the trapezoidal body-force impulse remains
    /// a separate contribution. Failure preserves both accepted epochs.
    pub fn advance_serial(
        &mut self,
        cx: &Cx<'_>,
        collisions: Option<&[L::Collision]>,
        mut forces_at_end: Vec<[f64; 3]>,
        origin: [f64; 3],
    ) -> Result<OpenBoundaryMomentumTransport3, GridError3<L::Error>> {
        let mut poll = || cx.checkpoint().map_err(|_| GridError3::Cancelled);
        let (receipt, step) =
            self.stage_serial_with_poll(&forces_at_end, collisions, origin, &mut poll)?;
        self.publish_with_force_swap(step, &mut forces_at_end);
        Ok(receipt)
    }

    /// Hold the accepted spatial force for one original native interval.
    pub fn step_serial(
        &mut self,
        cx: &Cx<'_>,
        collisions: Option<&[L::Collision]>,
        origin: [f64; 3],
    ) -> Result<OpenBoundaryMomentumTransport3, GridError3<L::Error>> {
        let mut poll = || cx.checkpoint().map_err(|_| GridError3::Cancelled);
        let (receipt, step) = self.inner.stage_serial_with_forces(
            self.inner.plan,
            collisions,
            Some(&self.forces),
            Some(&self.forces),
            origin,
            &mut poll,
        )?;
        core::mem::swap(&mut self.inner.grid.f, &mut self.inner.grid.next);
        self.inner.grid.completed_steps = step;
        self.inner.publish_wall_exchange();
        Ok(receipt)
    }

    // Thermal coupling uses this same production core and inspects the actual
    // native candidate before its own final source/balance gates publish.
    pub(super) fn stage_serial_with_poll(
        &mut self,
        forces_at_end: &[[f64; 3]],
        collisions: Option<&[L::Collision]>,
        origin: [f64; 3],
        poll: &mut impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<(OpenBoundaryMomentumTransport3, usize), GridError3<L::Error>> {
        self.inner.stage_serial_with_forces(
            self.inner.plan,
            collisions,
            Some(&self.forces),
            Some(forces_at_end),
            origin,
            poll,
        )
    }

    // The caller has admitted THIS endpoint slice in the staged interval and
    // completed every additional gate. Old accepted forces become reusable
    // caller workspace; publication allocates nothing and cannot refuse.
    pub(super) fn publish_with_force_swap(
        &mut self,
        step: usize,
        forces_at_end: &mut Vec<[f64; 3]>,
    ) {
        core::mem::swap(&mut self.inner.grid.f, &mut self.inner.grid.next);
        core::mem::swap(&mut self.forces, forces_at_end);
        self.inner.grid.completed_steps = step;
        self.inner.publish_wall_exchange();
    }

    /// Last accepted native exterior-wall impulse and work, in the same
    /// force-centered transaction as the population and spatial force epochs.
    pub fn last_wall_exchange(&self) -> Option<crate::lattice3::moving_walls::MovingWallMomentum3> {
        self.inner.last_wall_exchange()
    }
}

#[cfg(test)]
mod tests;
