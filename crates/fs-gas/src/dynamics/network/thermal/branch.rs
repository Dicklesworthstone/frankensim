//! Explicit geometry and original wall sources for resolved branch cells.

use super::*;

/// One resolved cell's original stationary wall loss and thermal storage.
/// A contact is one physical wall amount; sharing one wall mass among cells
/// requires a separate thermal owner rather than duplicate energy inventories.
#[derive(Clone, Debug, PartialEq)]
pub struct BranchWallContact<'law> {
    /// Original linear/Darcy drag and adiabatic/reservoir/finite thermal contact.
    pub contact: WallContact,
    /// Original material h/u storage, including nonlinear or latent inventory.
    pub heat_storage: Option<MaterialWallCaloric<'law>>,
}

#[derive(Clone)]
pub(super) struct BranchBinding {
    pub(super) mesh: Mesh1d,
    pub(super) start: usize,
}

impl<'law> NetworkWallOperator<'law> {
    /// Bind selected branch cells before the first source interval.
    ///
    /// Supply one entry per retained branch. None leaves that branch adiabatic
    /// and without wall drag. Some contains exactly one distinct contact per
    /// original cell, with complete conductance h*A [W/K] and wall amount.
    /// The original mesh is retained for exact geometry admission on every use.
    /// No source coefficient or finite thermal mass is inferred from the mesh.
    pub fn with_branch_contacts(
        mut self,
        network: &GasNetwork1d,
        contacts: Vec<Option<Vec<BranchWallContact<'law>>>>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(&mut cancelled)?;
        self.validate_binding(network)?;
        if self.last_report.is_some()
            || self.branches.iter().any(Option::is_some)
            || contacts.len() != network.branches.len()
        {
            return Err(Error::Invalid("branch wall binding state or topology"));
        }
        let mut count = self.plenum_count;
        for (input, branch) in contacts.iter().zip(&network.branches) {
            poll(&mut cancelled)?;
            if let Some(input) = input {
                if input.len() != branch.states.len() {
                    return Err(Error::Invalid("branch wall contact cell count"));
                }
                count = count.checked_add(input.len()).ok_or(Error::Allocation)?;
            }
        }
        let extra = count - self.plenum_count;
        self.volumes
            .try_reserve_exact(extra)
            .map_err(|_| Error::Allocation)?;
        self.contacts
            .try_reserve_exact(extra)
            .map_err(|_| Error::Allocation)?;
        self.candidate_contacts
            .try_reserve_exact(extra)
            .map_err(|_| Error::Allocation)?;
        self.heat_storage
            .try_reserve_exact(extra)
            .map_err(|_| Error::Allocation)?;
        for (b, input) in contacts.into_iter().enumerate() {
            let Some(input) = input else {
                continue;
            };
            poll(&mut cancelled)?;
            let branch = &network.branches[b];
            let start = self.contacts.len();
            for (i, input) in input.into_iter().enumerate() {
                poll(&mut cancelled)?;
                let contact = admit_contact(input.contact, input.heat_storage.as_ref())?;
                self.volumes.push(branch.mesh.volumes()[i]);
                self.contacts.push(contact);
                self.candidate_contacts.push(contact);
                self.heat_storage.push(input.heat_storage);
            }
            self.branches[b] = Some(BranchBinding {
                mesh: branch.mesh.clone(),
                start,
            });
        }
        poll(&mut cancelled)?;
        Ok(self)
    }

    /// Accepted contacts along a bound branch, in its original cell order.
    pub fn branch_contacts(&self, branch: usize) -> Option<&[WallContact]> {
        self.branches
            .get(branch)?
            .as_ref()
            .map(|binding| &self.contacts[binding.start..binding.start + binding.mesh.cell_count()])
    }

    /// The original finite-storage owners along that same branch.
    pub fn branch_heat_storage(
        &self,
        branch: usize,
    ) -> Option<&[Option<MaterialWallCaloric<'law>>]> {
        self.branches.get(branch)?.as_ref().map(|binding| {
            &self.heat_storage[binding.start..binding.start + binding.mesh.cell_count()]
        })
    }
}
