//! Mechanical loads from the SAME current and magnetic stage as induction.
//!
//! P1 Lorentz density `J cross B` is integrated against each mechanical P1
//! shape function. The consistent tetrahedral mass weights preserve net force,
//! torque and virtual work for every P1 velocity, including rigid motion.
//! Averaging the density or assigning `volume/4` to its vertex samples would
//! generally preserve force but corrupt torque and mechanical work.
//!
//! These are conduction-current volume loads, not magnetization or interface
//! loads. Do not add an overlapping enclosing Maxwell-stress load to the same
//! body. Geometry and motion remain those admitted by the induction candidate;
//! this transfer does not itself update motion or close mechanical feedback.

pub(crate) mod load;

pub use load::{SpatialLorentzLoadError, SpatialLorentzLoadLimits, SpatialLorentzLoads};

use super::{SpatialInductionCandidate, SpatialInductionCell};
use fs_material::conductor::spatial::SpatialConductorMaterial;
use load::{CellLoad, integrate};

impl<C: SpatialConductorMaterial> SpatialInductionCandidate<'_, C> {
    /// Consistent nodal loads on selected original conducting tetrahedra.
    ///
    /// Uses the midpoint `J cross B` from this candidate's accepted field
    /// update and its original mesh. Cell IDs may arrive in any order but must
    /// be distinct and conducting. Selection is canonicalized before summing.
    /// Nodal forces remain in original mesh coordinates and vertex order, so
    /// disjoint selections can supply different mechanical bodies. Nodes not
    /// incident on the selection receive zero force.
    ///
    /// Pairing the returned nodal forces with the prescribed P1 velocity,
    /// then multiplying by this step's duration, recovers exactly the selected
    /// cells' Lorentz work (up to floating-point reduction). A mechanical
    /// receiver must use the same time stage and trajectory for that account.
    /// Cancellation or refusal produces no partial load and changes no state.
    pub fn lorentz_loads(
        &self,
        selected_cells: &[usize],
        limits: SpatialLorentzLoadLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialLorentzLoads, SpatialLorentzLoadError> {
        let problem = self.state().problem();
        recover(
            problem.positions_m(),
            &problem.mesh().tets,
            &self.step().cells,
            selected_cells,
            limits,
            cancelled,
        )
    }
}

impl<C: SpatialConductorMaterial> super::nonlinear::SpatialNonlinearInductionCandidate<'_, C> {
    /// Consistent conduction-current loads from this nonlinear candidate.
    ///
    /// Selection and frame contracts match
    /// [`SpatialInductionCandidate::lorentz_loads`]. The retained nonlinear
    /// update uses its accepted ENDPOINT current and magnetic field, so these
    /// endpoint loads pair with that update's prescribed velocity and duration.
    /// Saturation does not change the P1 mechanical integration rule.
    pub fn lorentz_loads(
        &self,
        selected_cells: &[usize],
        limits: SpatialLorentzLoadLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialLorentzLoads, SpatialLorentzLoadError> {
        let problem = self.state().problem();
        recover(
            problem.positions_m(),
            &problem.mesh().tets,
            &self.step().cells,
            selected_cells,
            limits,
            cancelled,
        )
    }
}

fn recover(
    positions_m: &[[f64; 3]],
    tets: &[[u32; 4]],
    cells: &[SpatialInductionCell],
    selected_cells: &[usize],
    limits: SpatialLorentzLoadLimits,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<SpatialLorentzLoads, SpatialLorentzLoadError> {
    use SpatialLorentzLoadError::{Cancelled, Invalid};
    if cancelled() {
        return Err(Cancelled);
    }
    if selected_cells.is_empty()
        || selected_cells.len() > limits.max_cells
        || positions_m.len() > limits.max_vertices
    {
        return Err(Invalid("empty or over-budget Lorentz load selection"));
    }
    // Induction recovery emits its source-bound cell records in increasing
    // original cell order. Binary lookup avoids scanning every conductor for
    // each body; only the admitted selection is allocated and canonicalized.
    let mut selected = Vec::with_capacity(selected_cells.len());
    for &cell in selected_cells {
        if cancelled() {
            return Err(Cancelled);
        }
        let record = cells
            .binary_search_by_key(&cell, |record| record.cell)
            .ok()
            .and_then(|index| cells.get(index))
            .ok_or(Invalid(
                "Lorentz selection requires original conducting cells",
            ))?;
        let vertices = *tets
            .get(cell)
            .ok_or(Invalid("Lorentz cell is outside its original mesh"))?;
        selected.push((
            cell,
            CellLoad {
                vertices,
                volume_m3: record.volume_m3,
                force_density_n_m3: &record.lorentz_force_density_n_m3,
            },
        ));
    }
    selected.sort_unstable_by_key(|entry| entry.0);
    if selected.windows(2).any(|pair| pair[0].0 == pair[1].0) {
        return Err(Invalid("Lorentz load selection contains a repeated cell"));
    }
    integrate(
        positions_m,
        selected.into_iter().map(|entry| entry.1),
        limits,
        cancelled,
    )
}
