//! Pathline history owned by the complete reacting physical interval.
//!
//! The reactor retains chemistry, while its borrowed molecular operator retains
//! the shared SGS histories. Only the D/C/A/C/D transaction may prepare and
//! commit them; molecular substeps and chemical source solves never age memory.

use super::super::reactive::{ReactingEuler3d, ReactiveError};
use super::*;
use sgs::lagrangian::{LagrangianFilter, LagrangianReports, LagrangianTrial};

impl ReactingEuler3d {
    /// Complete current reactor/operator capacity plus fresh momentum history.
    /// Add any unallocated RKL2 buffers and externally supplied forcing scratch.
    pub fn required_lagrangian_dynamic_sgs_workspace_bytes<L: TransportLaw>(
        &self,
        operator: &MolecularTransport3d<L>,
    ) -> Result<usize, ReactiveError> {
        self.required_lagrangian_reacting_bytes(
            operator,
            SubgridTransport::required_lagrangian_bytes(self.flow())?,
        )
    }

    /// Complete current reactor/operator capacity plus independent stress,
    /// temperature and common-species histories. Their additional payload is
    /// 448 bytes/cell, independent of species count. Add unallocated RKL2
    /// buffers and externally supplied forcing scratch.
    pub fn required_lagrangian_dynamic_species_sgs_workspace_bytes<L: TransportLaw>(
        &self,
        operator: &MolecularTransport3d<L>,
    ) -> Result<usize, ReactiveError> {
        self.required_lagrangian_reacting_bytes(
            operator,
            SubgridTransport::required_lagrangian_species_bytes(self.flow())?,
        )
    }

    /// Attach fresh periodic Lagrangian stress with explicit Prandtl/Schmidt
    /// heat and species transport to this reactor's molecular operator.
    ///
    /// The complete reacting split samples its accepted START state once and
    /// holds the coefficient through both molecular halves, both original
    /// chemistry solves, SSP/RKL stages and all whole-interval retries.
    /// Only complete physical acceptance commits the history and its report.
    /// The next sample uses the previous actual accepted duration.
    ///
    /// Attachment is atomic and requires the same original gas, grid and
    /// boundaries and an operator without another SGS closure. The supplied
    /// cap covers this reactor, the molecular operator and all added storage.
    /// Molecular advancement outside a prepared physical interval refuses.
    pub fn attach_lagrangian_dynamic_sgs<L: TransportLaw>(
        &self,
        operator: &mut MolecularTransport3d<L>,
        config: sgs::MixtureLagrangianSubgridConfig,
        maximum_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), ReactiveError> {
        poll(&mut cancelled)?;
        config.validate()?;
        let required = self.required_lagrangian_dynamic_sgs_workspace_bytes(operator)?;
        admit_workspace(required, maximum_workspace_bytes)?;
        let subgrid = SubgridTransport::new_lagrangian(self.flow(), config)?;
        self.attach_reacting_lagrangian(operator, subgrid, maximum_workspace_bytes, &mut cancelled)
    }

    /// Attach independent physical-time memories for stress, temperature and
    /// common species mixing to the original reacting molecular operator.
    ///
    /// All mass fractions enter the shared species Germano objective without a
    /// transpose or per-species history. Chemistry remains the original
    /// source-owned finite-rate solve; no turbulent reaction-rate or flame
    /// wrinkling model is inferred. Current stage composition, strain, density,
    /// source Cp/Cv, Maxwell--Stefan flux and referenced partial enthalpies stay
    /// active in both conservative molecular halves. No extra heating source
    /// is introduced.
    ///
    /// All three histories and reports commit with the complete D/C/A/C/D step.
    /// Authored memory times, first-order pathline/coefficient integration and
    /// fully periodic Cartesian exploratory admission are unchanged. Attachment
    /// is atomic, starts fresh history, and retains the original source and gas.
    pub fn attach_lagrangian_dynamic_species_sgs<L: TransportLaw>(
        &self,
        operator: &mut MolecularTransport3d<L>,
        config: sgs::MixtureLagrangianSpeciesSubgridConfig,
        maximum_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), ReactiveError> {
        poll(&mut cancelled)?;
        config.validate()?;
        let required = self.required_lagrangian_dynamic_species_sgs_workspace_bytes(operator)?;
        admit_workspace(required, maximum_workspace_bytes)?;
        let subgrid = SubgridTransport::new_lagrangian_species(self.flow(), config)?;
        self.attach_reacting_lagrangian(operator, subgrid, maximum_workspace_bytes, &mut cancelled)
    }

    fn required_lagrangian_reacting_bytes<L: TransportLaw>(
        &self,
        operator: &MolecularTransport3d<L>,
        additional: usize,
    ) -> Result<usize, ReactiveError> {
        operator.admit_lagrangian_flow(self.flow())?;
        if operator.sgs.is_some() {
            return Err(ReactiveError::Invalid(
                "reacting Lagrangian attachment requires a molecular operator without SGS",
            ));
        }
        self.workspace_bytes()
            .checked_add(operator.workspace_bytes())
            .and_then(|bytes| bytes.checked_add(additional))
            .ok_or(ReactiveError::Invalid(
                "reacting Lagrangian workspace size overflow",
            ))
    }

    fn attach_reacting_lagrangian<L: TransportLaw>(
        &self,
        operator: &mut MolecularTransport3d<L>,
        subgrid: SubgridTransport,
        maximum_workspace_bytes: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), ReactiveError> {
        let bytes = self.required_lagrangian_reacting_bytes(operator, subgrid.workspace_bytes())?;
        admit_workspace(bytes, maximum_workspace_bytes)?;
        poll(cancelled)?;
        // No fallible operation follows publication. Existing physical fields,
        // molecular reports and loads were never changed during admission.
        operator.sgs = Some(subgrid);
        Ok(())
    }
}

impl<L: TransportLaw> MolecularTransport3d<L> {
    /// Fixed-ratio pathline controls, or None for independent/other closures.
    pub fn lagrangian_dynamic_subgrid_config(&self) -> Option<sgs::MixtureLagrangianSubgridConfig> {
        self.sgs
            .as_ref()
            .and_then(SubgridTransport::lagrangian_config)
    }

    /// Independent heat/common-species pathline controls, if selected.
    pub fn lagrangian_dynamic_species_subgrid_config(
        &self,
    ) -> Option<sgs::MixtureLagrangianSpeciesSubgridConfig> {
        self.sgs
            .as_ref()
            .and_then(SubgridTransport::lagrangian_species_config)
    }

    /// Last committed physical START sample, distinct from the flow endpoint.
    /// None until the complete composing interval accepts.
    pub fn lagrangian_subgrid_sample_time(&self) -> Option<Time> {
        self.sgs
            .as_ref()
            .and_then(SubgridTransport::lagrangian_sample_time)
    }

    /// Last committed accepted-start momentum fit.
    pub fn lagrangian_subgrid_report(&self) -> Option<sgs::DynamicSgsReport> {
        self.sgs
            .as_ref()
            .and_then(SubgridTransport::lagrangian_report)
    }

    /// Last committed temperature fit, or None for fixed-Prandtl transport.
    pub fn lagrangian_subgrid_heat_report(&self) -> Option<sgs::DynamicScalarReport> {
        self.sgs
            .as_ref()
            .and_then(SubgridTransport::lagrangian_heat_report)
    }

    /// Last committed common-species fit, or None for fixed-Schmidt transport.
    pub fn lagrangian_subgrid_species_report(&self) -> Option<sgs::DynamicScalarReport> {
        self.sgs
            .as_ref()
            .and_then(SubgridTransport::lagrangian_species_report)
    }

    fn admit_lagrangian_flow(&self, flow: &MixtureEuler3d) -> Result<(), DissipationError> {
        if self.gas() != flow.gas() || self.grid != flow.grid || self.boundaries != flow.boundaries
        {
            return Err(DissipationError::Invalid(
                "reacting Lagrangian molecular field ownership differs",
            ));
        }
        Ok(())
    }

    pub(in crate::compressible::spatial::mixture) fn take_lagrangian_history(
        &mut self,
    ) -> Option<LagrangianFilter> {
        self.sgs
            .as_mut()
            .and_then(SubgridTransport::take_lagrangian_filter)
    }

    pub(in crate::compressible::spatial::mixture) fn prepare_lagrangian_history<'a>(
        &mut self,
        filter: &'a mut LagrangianFilter,
        flow: &MixtureEuler3d,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<LagrangianTrial<'a>, DissipationError> {
        // The externally borrowed operator can come from a different field.
        // Validate ownership before indexing any prepared history/staging bank.
        self.admit_lagrangian_flow(flow)?;
        self.sgs
            .as_mut()
            .ok_or(DissipationError::Invalid(
                "missing reacting Lagrangian owner",
            ))?
            .prepare_lagrangian(filter, flow, cancelled)
    }

    pub(in crate::compressible::spatial::mixture) fn restore_lagrangian_history(
        &mut self,
        filter: LagrangianFilter,
        committed: Option<LagrangianReports>,
    ) {
        self.sgs
            .as_mut()
            .expect("a reacting step retains its original Lagrangian owner")
            .restore_lagrangian_filter(filter, committed);
    }
}
