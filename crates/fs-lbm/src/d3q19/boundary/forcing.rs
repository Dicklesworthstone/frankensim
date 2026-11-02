//! Atomic uniform flow preparation and held-force changes in the original Guo frame.

use super::*;

/// Refusal before a uniform flow or force-frame proposal is published.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BoundaryForcingError3 {
    /// Cancellation at a lattice-tile boundary or before publication.
    Cancelled,
    /// Invalid physical data or unsupported boundary/collision combination.
    Invalid(&'static str),
    /// Nonpositive density or nonfinite force, velocity or proposed populations.
    Unphysical {
        /// Original cell coordinates.
        cell: [usize; 3],
    },
}
impl core::fmt::Display for BoundaryForcingError3 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "D3Q19 flow preparation: {self:?}")
    }
}
impl std::error::Error for BoundaryForcingError3 {}

impl BoundaryGrid3 {
    /// Held uniform volume-force density, excluding any per-cell contribution.
    pub const fn uniform_force_density(&self) -> [f64; 3] {
        self.force
    }

    /// Replace the uniform part of a held volume-force density, in lattice units.
    ///
    /// The original Guo source transform preserves physical mass, momentum and
    /// force-corrected stress at each fluid cell. This changes the force held for
    /// subsequent steps; it advances no time and supplies no impulse or work.
    /// Existing per-cell forces and their selected spatial shape remain intact.
    /// A discontinuous pressure-gradient change therefore creates no artificial
    /// half-force impulse. Smooth-history trapezoidal forcing remains the
    /// separately declared `set_force_field_scale` convention.
    ///
    /// Uses existing collision scratch with no new population or force-field
    /// allocation. Cancellation or a late-cell refusal leaves accepted fields
    /// and forcing unchanged. An identical force preserves population bits.
    pub fn rebind_uniform_force(
        &mut self,
        force: [f64; 3],
        gate: &CancelGate,
    ) -> Result<(), BoundaryForcingError3> {
        poll(gate)?;
        if !force.iter().all(|value| value.is_finite()) {
            return Err(BoundaryForcingError3::Invalid(
                "uniform force must be finite",
            ));
        }
        if force != [0.0; 3]
            && (!self.collision_model.supports_body_force()
                || Face3::ALL
                    .into_iter()
                    .any(|face| self.boundaries.face(face).is_open()))
        {
            return Err(BoundaryForcingError3::Invalid(
                "nonzero uniform force requires a forcing-compatible closed or periodic domain",
            ));
        }
        if self.force == force {
            return Ok(());
        }
        for tile in 0..self.f[0].len() {
            poll(gate)?;
            for lane in 0..TILE_CELLS {
                let mut populations = core::array::from_fn(|q| self.f[q][tile].0[lane]);
                if self.solid[tile] & (1u64 << lane) == 0 {
                    let previous = cell_force(self.force, self.force_field.as_deref(), tile, lane);
                    let next = cell_force(force, self.force_field.as_deref(), tile, lane);
                    let cell = self.coords(tile, lane);
                    if !previous.iter().chain(&next).all(|value| value.is_finite()) {
                        return Err(BoundaryForcingError3::Unphysical { cell });
                    }
                    if previous != next {
                        // The original moment observer asserts positive density;
                        // turn a late invalid cell into this atomic typed refusal.
                        admit_populations(&populations, cell)?;
                        let (rho, momentum) = self.raw_moments(cell[0], cell[1], cell[2]);
                        let velocity = core::array::from_fn(|axis| {
                            (momentum[axis] + 0.5 * previous[axis]) / rho
                        });
                        if !velocity.iter().all(|value| value.is_finite()) {
                            return Err(BoundaryForcingError3::Unphysical { cell });
                        }
                        crate::thermal::forcing::reframe3(
                            &mut populations,
                            velocity,
                            previous,
                            next,
                        );
                        admit_populations(&populations, cell)?;
                    }
                }
                for (proposed, value) in self.post.iter_mut().zip(populations) {
                    proposed[tile].0[lane] = value;
                }
            }
        }
        poll(gate)?;
        core::mem::swap(&mut self.f, &mut self.post);
        self.force = force;
        Ok(())
    }

    /// Prepare uniform PHYSICAL density/velocity using the original equilibrium.
    ///
    /// Existing uniform and per-cell forces select each cell's Guo half-source
    /// frame. Thus the requested velocity is the force-corrected physical
    /// velocity even when a force is already held. Solid cells are unchanged.
    ///
    /// This deliberately replaces the fluid initial state and its stress; it
    /// is not a momentum-conserving operation on an evolved flow. No body state,
    /// boundary topology, force, viscosity or clock changes. Material adapters
    /// restrict this initializer to their initial accepted time. Wall mismatch
    /// describes a startup transient, not an already solved steady flow.
    /// The caller supplies physical Mach/regime admission for its scenario.
    pub fn initialize_uniform_equilibrium(
        &mut self,
        density: f64,
        velocity: [f64; 3],
        gate: &CancelGate,
    ) -> Result<(), BoundaryForcingError3> {
        poll(gate)?;
        if !density.is_finite() || density <= 0.0 || !velocity.iter().all(|value| value.is_finite())
        {
            return Err(BoundaryForcingError3::Invalid(
                "initial density must be positive and velocity finite",
            ));
        }
        let equilibrium = equilibrium3(density, velocity);
        if !equilibrium.iter().all(|value| value.is_finite()) {
            return Err(BoundaryForcingError3::Invalid(
                "unrepresentable initial equilibrium",
            ));
        }
        for tile in 0..self.f[0].len() {
            poll(gate)?;
            for lane in 0..TILE_CELLS {
                let mut populations = core::array::from_fn(|q| self.f[q][tile].0[lane]);
                if self.solid[tile] & (1u64 << lane) == 0 {
                    populations = equilibrium;
                    let force = cell_force(self.force, self.force_field.as_deref(), tile, lane);
                    if force != [0.0; 3] {
                        crate::thermal::forcing::reframe3(
                            &mut populations,
                            velocity,
                            [0.0; 3],
                            force,
                        );
                    }
                    admit_populations(&populations, self.coords(tile, lane))?;
                }
                for (proposed, value) in self.post.iter_mut().zip(populations) {
                    proposed[tile].0[lane] = value;
                }
            }
        }
        poll(gate)?;
        core::mem::swap(&mut self.f, &mut self.post);
        Ok(())
    }
}

fn admit_populations(
    populations: &[f64; Q3],
    cell: [usize; 3],
) -> Result<(), BoundaryForcingError3> {
    let density = populations.iter().sum::<f64>();
    if populations.iter().all(|value| value.is_finite()) && density.is_finite() && density > 0.0 {
        Ok(())
    } else {
        Err(BoundaryForcingError3::Unphysical { cell })
    }
}
fn poll(gate: &CancelGate) -> Result<(), BoundaryForcingError3> {
    if gate.is_requested() {
        Err(BoundaryForcingError3::Cancelled)
    } else {
        Ok(())
    }
}
