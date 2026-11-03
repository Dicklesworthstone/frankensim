//! Reservoir exchange of the original streamed and regularized open faces.
//!
//! Count the temporary open-link reflection in the existing collision walk,
//! then the actual reconstruction differences in the native face order. This
//! owns all momentum supplied by that numerical boundary operator, including
//! replacement of tangential populations on face interiors. It is not a
//! residual inferred from a whole-domain conservation check.

use super::super::{BoundaryStepError3, Tile};
use super::{
    BoundaryGrid3, CancelGate, E3, OPP3, PartialSaturationError3, Q3, Result, TILE_CELLS,
    checkpoint, invalid,
};
use fs_math::dd::Dd;

/// Mass and momentum delivered TO the fluid by prescribed velocity/pressure
/// reservoirs in one native D3Q19 tick. The opposite belongs to the external
/// reservoir, not to an immersed body or wall. All quantities use lattice
/// units; angular momentum uses the same fixed origin as the body receipt.
///
/// This is the measured exchange of the regularized population operator,
/// including its non-advective reconstruction corrections. It is not a
/// pressure-work, internal-energy or complete thermodynamic balance.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct OpenBoundaryMomentumTransport3 {
    /// Actual population mass introduced at open faces.
    pub mass_change: f64,
    /// Actual population linear momentum introduced at open faces. With held
    /// Guo forcing this is also the change in physical momentum.
    pub impulse: [f64; 3],
    /// Moment of the actual open-face momentum increments about the origin.
    pub angular_impulse: [f64; 3],
    /// Open incoming links counted once, excluding wall-owned corner links.
    pub open_links: usize,
    /// Native face/cell reconstruction events in the original canonical order.
    /// The boundary owner admits open faces along at most one grid axis.
    pub reconstructed_cells: usize,
}

// Opposing reservoir pressure impulses can dwarf the net transfer. Retain
// their low bits until the shared receipt is rounded, as for body reactions.
#[derive(Default)]
pub(super) struct OpenTransport {
    mass: Dd,
    impulse: [Dd; 3],
    angular: [Dd; 3],
    links: usize,
    cells: usize,
}

impl OpenTransport {
    pub(super) fn add_links(
        &mut self,
        post: &[f64; Q3],
        mask: u32,
        cell: [usize; 3],
        origin: [f64; 3],
    ) {
        if mask == 0 {
            return;
        }
        let arm = core::array::from_fn(|i| Dd::from_f64(cell[i] as f64 + 0.5 - origin[i]));
        for q in 1..Q3 {
            if mask & (1u32 << q) != 0 {
                // The original pull operator temporarily reflects post[opp].
                // Removal of outgoing -c and insertion of incoming +c give
                // 2*c*post[opp], with zero mass change. A half-link moment
                // equals the center moment because the offset is parallel c.
                self.add_momentum(q, Dd::from_f64(post[OPP3[q]]) * Dd::from_f64(2.0), arm);
                self.links += 1;
            }
        }
    }

    fn add_momentum(&mut self, q: usize, population: Dd, arm: [Dd; 3]) {
        let c = [E3[q].0, E3[q].1, E3[q].2];
        let impulse: [Dd; 3] = core::array::from_fn(|i| population * Dd::from_f64(f64::from(c[i])));
        for i in 0..3 {
            let (j, k) = ((i + 1) % 3, (i + 2) % 3);
            self.impulse[i] = self.impulse[i] + impulse[i];
            self.angular[i] = self.angular[i] + arm[j] * impulse[k] - arm[k] * impulse[j];
        }
    }

    fn add_reconstruction(
        &mut self,
        before: [f64; Q3],
        after: [f64; Q3],
        cell: [usize; 3],
        origin: [f64; 3],
    ) {
        let arm = core::array::from_fn(|i| Dd::from_f64(cell[i] as f64 + 0.5 - origin[i]));
        for q in 0..Q3 {
            let difference = Dd::from_f64(after[q]) - Dd::from_f64(before[q]);
            self.mass = self.mass + difference;
            self.add_momentum(q, difference, arm);
        }
        self.cells += 1;
    }

    fn finish(self) -> Result<OpenBoundaryMomentumTransport3> {
        let receipt = OpenBoundaryMomentumTransport3 {
            mass_change: self.mass.to_f64(),
            impulse: self.impulse.map(Dd::to_f64),
            angular_impulse: self.angular.map(Dd::to_f64),
            open_links: self.links,
            reconstructed_cells: self.cells,
        };
        if [receipt.mass_change]
            .into_iter()
            .chain(receipt.impulse)
            .chain(receipt.angular_impulse)
            .all(f64::is_finite)
        {
            Ok(receipt)
        } else {
            Err(invalid("open-boundary momentum account overflowed"))
        }
    }
}

impl BoundaryGrid3 {
    pub(super) fn reconstruct_immersed_open_faces(
        &mut self,
        next: &mut [Vec<Tile>; Q3],
        mut transport: OpenTransport,
        origin: [f64; 3],
        gate: &CancelGate,
    ) -> Result<OpenBoundaryMomentumTransport3> {
        checkpoint(gate)?;
        if transport.links == 0 {
            return transport.finish();
        }
        if self.open_cells.is_none() {
            self.open_cells = Some(self.compile_open_cells());
        }
        let cells = self.open_cells.take().unwrap_or_default();
        let result = (|| {
            for (index, &(cell, interior, face)) in cells.iter().enumerate() {
                if index.is_multiple_of(TILE_CELLS) {
                    checkpoint(gate)?;
                }
                let values = self
                    .open_boundary_populations(next, cell, interior, face, None)
                    .map_err(|source| match source {
                        BoundaryStepError3::Unphysical { cell } => PartialSaturationError3::Cell {
                            cell,
                            reason: "native open-face reconstruction is unphysical",
                        },
                        _ => invalid("native open-face reconstruction refused"),
                    })?;
                let (tile, lane) = self.addr(cell[0], cell[1], cell[2]);
                let before = next.each_ref().map(|field| field[tile].0[lane]);
                transport.add_reconstruction(before, values, cell, origin);
                for (field, value) in next.iter_mut().zip(values) {
                    field[tile].0[lane] = value;
                }
            }
            checkpoint(gate)?;
            transport.finish()
        })();
        self.open_cells = Some(cells);
        result
    }
}
