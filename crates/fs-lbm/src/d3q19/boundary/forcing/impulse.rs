//! Impulsive volume loading with the original tile-major publication scratch.

use super::super::{BoundaryGrid3, E3, TILE_CELLS, W3, cell_force};
use crate::impulse::{ImpulseControls, ImpulseError, ImpulseReport, ImpulseSum, poll};
use fs_exec::CancelGate;

impl BoundaryGrid3 {
    /// Apply instantaneous impulse density to the original fluid populations.
    ///
    /// Inputs cover every original cell in x-fastest order, in lattice momentum
    /// units; solid-cell entries must be zero. The kernel retains physical
    /// nonequilibrium stress and the held Guo force while changing `rho*u` by
    /// the prescribed impulse. It advances no time, collision or streaming.
    /// Later native evolution owns pressure propagation and wall reactions.
    ///
    /// Uses existing collision scratch with no population/force allocation.
    /// A late refusal, exceeded readback budget or cancellation leaves all
    /// accepted populations and fields intact. Zero inputs preserve all bits.
    /// A published nonzero impulse locks the original solid topology.
    ///
    /// Immersed partial-saturation cells are full kinetic cells here. A coupled
    /// owner must account for its own fictitious interior-fluid momentum; this
    /// operation alone does not apply an equal/opposite impulse to a rigid body.
    pub fn apply_impulse_density(
        &mut self,
        impulses: &[[f64; 3]],
        controls: ImpulseControls,
        gate: &CancelGate,
    ) -> Result<ImpulseReport<3>, ImpulseError> {
        poll(gate)?;
        let count = self
            .nx
            .checked_mul(self.ny)
            .and_then(|n| n.checked_mul(self.nz))
            .ok_or(ImpulseError::Invalid("invalid grid size"))?;
        controls.admit(count)?;
        if impulses.len() != count {
            return Err(ImpulseError::Invalid(
                "impulse field must cover the original grid",
            ));
        }
        let directions = E3.map(|(x, y, z)| [f64::from(x), f64::from(y), f64::from(z)]);
        let mut sum = ImpulseSum::new();
        for tile in 0..self.f[0].len() {
            poll(gate)?;
            for lane in 0..TILE_CELLS {
                let [x, y, z] = self.coords(tile, lane);
                let index = (z * self.ny + y) * self.nx + x;
                let impulse = impulses[index];
                if self.solid[tile] & (1u64 << lane) != 0 && impulse != [0.0; 3] {
                    return Err(ImpulseError::Invalid(
                        "solid cells cannot receive a fluid impulse",
                    ));
                }
                let mut next = core::array::from_fn(|q| self.f[q][tile].0[lane]);
                let force = cell_force(self.force, self.force_field.as_deref(), tile, lane);
                sum.kick(
                    &mut next,
                    &directions,
                    &W3,
                    [0.0; 3],
                    force,
                    impulse,
                    controls,
                    index,
                )?;
                for (field, value) in self.post.iter_mut().zip(next) {
                    field[tile].0[lane] = value;
                }
            }
        }
        let report = sum.finish(controls)?;
        poll(gate)?;
        if report.affected_cells != 0 {
            core::mem::swap(&mut self.f, &mut self.post);
            self.topology_locked = true;
        }
        Ok(report)
    }
}
