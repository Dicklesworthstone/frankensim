//! Immutable material heat storage in the original network wall transaction.

use super::*;

impl<'law> NetworkWallOperator<'law> {
    /// Bind original material calorics to selected finite walls before evolution.
    ///
    /// Each `Some` consumes a complete fixed-amount storage binding. Original
    /// contact temperatures or h/u and conductances remain authoritative. A
    /// temperature contact's capacity becomes the source's refreshed tangent;
    /// an energy contact retains latent storage directly. `None` retains the
    /// authored wall path. Construct the contact on its same source first.
    /// Distinct contacts still represent distinct wall masses even if their
    /// immutable source tables share an Arc. A shared physical wall needs its
    /// own spatial thermal owner instead of duplicating its storage here.
    ///
    /// Ordinary, reacting and dissipating network compositions reuse this same
    /// binding and original candidate wall temperatures. A failed later cell,
    /// reaction half, source domain or cancellation publishes no accepted wall.
    pub fn with_material_heat_storage(
        mut self,
        branches: Vec<Vec<Option<MaterialWallCaloric<'law>>>>,
        plenums: Vec<Option<MaterialWallCaloric<'law>>>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, WallExchangeError> {
        poll(&mut cancelled)?;
        if self.last_report.is_some()
            || self.trial_active
            || self.caloric.iter().any(Option::is_some)
            || branches.len() != self.branches.len()
            || plenums.len() != self.plenum_volumes.len()
            || !branches
                .iter()
                .flatten()
                .chain(&plenums)
                .any(Option::is_some)
        {
            return Err(WallExchangeError::Invalid(
                "material wall binding state or topology",
            ));
        }
        for (input, shape) in branches.iter().zip(&self.branches) {
            if input.len() != shape.mesh.cell_count() {
                return Err(WallExchangeError::Invalid(
                    "material wall branch cell count",
                ));
            }
        }
        for (index, law) in branches.into_iter().flatten().chain(plenums).enumerate() {
            poll(&mut cancelled)?;
            if let Some(law) = &law {
                self.contacts[index] = law.bind_contact(self.contacts[index])?;
            }
            self.caloric[index] = law;
        }
        poll(&mut cancelled)?;
        Ok(self)
    }

    /// Original material owners in this branch's original cell order.
    pub fn branch_heat_storage(
        &self,
        branch: usize,
    ) -> Option<&[Option<MaterialWallCaloric<'law>>]> {
        self.branches
            .get(branch)
            .map(|shape| &self.caloric[shape.start..shape.start + shape.mesh.cell_count()])
    }

    /// Original material owners in the original finite-plenum order.
    pub fn plenum_heat_storage(&self) -> &[Option<MaterialWallCaloric<'law>>] {
        &self.caloric[self.plenum_start..]
    }
}
