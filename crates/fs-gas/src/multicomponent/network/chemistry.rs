//! Prepared fixed-volume chemical pressure feedback from the original reactor.
use super::Plenum;
use crate::multicomponent::MAX_COMPONENTS;
use crate::multicomponent::reactive::{CellChemistry, ChemistryReport, ReactiveError};
use fs_qty::{Time, Volume};

impl Plenum {
    /// Advance a stationary finite inventory with the original cell chemistry.
    ///
    /// The prepared source owns the exact species bank, rate laws, fixed-U,V
    /// solve and formation-energy convention. Reaction changes the stored
    /// sensible-energy gauge by the opposite species-reference change; no
    /// separate heat-release term is added. Volume and zero momentum persist.
    ///
    /// The already allocated output and returned source receipt are candidates.
    /// The composing owner publishes them together. Every output write follows
    /// the last cancellation poll; successful warm calls allocate no storage.
    /// `cell` is a caller-supplied diagnostic index, not a physics selector.
    pub fn react_candidate_into(
        &self,
        chemistry: &mut CellChemistry,
        duration: Time,
        cell: usize,
        output: &mut Self,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ChemistryReport, ReactiveError> {
        if cancelled() {
            return Err(ReactiveError::Cancelled);
        }
        let n = chemistry.gas().species_count();
        if self.state.partial_densities.len() != n
            || output.state.partial_densities.len() != n
            || self.state.momentum != 0.0
        {
            return Err(ReactiveError::Invalid(
                "prepared stationary plenum chemistry species shape",
            ));
        }
        let mut row = [0.0; MAX_COMPONENTS];
        row[..n].copy_from_slice(&self.state.partial_densities);
        row[n + 1] = self.state.energy;
        let mut report = ChemistryReport::default();
        chemistry.advance_axial(
            &mut row[..n + 2],
            Volume::new(self.volume_m3),
            duration,
            cell,
            &mut report,
            &mut cancelled,
        )?;
        if row[n] != 0.0 || !row[n + 1].is_finite() {
            return Err(ReactiveError::Invalid(
                "unrepresentable stationary plenum chemistry candidate",
            ));
        }
        if cancelled() {
            return Err(ReactiveError::Cancelled);
        }
        output.volume_m3 = self.volume_m3;
        output.state.partial_densities.copy_from_slice(&row[..n]);
        output.state.momentum = 0.0;
        output.state.energy = row[n + 1];
        Ok(report)
    }
}

#[cfg(test)]
mod tests;
