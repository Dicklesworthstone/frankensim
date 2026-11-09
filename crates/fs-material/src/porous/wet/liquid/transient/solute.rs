//! Passive neutral solutes carried by the original porous liquid interval.
//!
//! This is dilute Fick transport in the stationary solid frame:
//! J_k=q*c_k-(theta*D_k/tau)*grad(c_k), theta=phi*S_l.
//! Each solute has its explicit aqueous diffusivity; tau is the transport divisor
//! in this formula. Geometric or acoustic tortuosity conventions require an
//! explicit source-supported conversion. The same caller
//! ordering is required in every constituent array. The caller qualifies D/tau
//! and the passive dilute approximation at the actual aqueous state. Solutes do not alter water
//! pressure, saturation, viscosity, density or retention. Activity/cross-friction,
//! ionic migration and reaction need their own laws. Longitudinal mechanical
//! dispersion has a separate original-source interval adapter.
//! Explicit linear and original-source Langmuir extensions add reversible
//! immobile storage; sorption is never inferred from aqueous D or tortuosity.

use super::{LiquidSaturationCandidate, LiquidSaturationColumn};

mod dispersion;
mod langmuir;
pub use dispersion::{
    DISPERSIVITY_PARAMETER, DispersionCellState, DispersiveSoluteStep, LIQUID_SATURATION_AXIS,
    LONGITUDINAL_DISPERSION_LAW, LongitudinalDispersivity, PORE_SPEED_AXIS, POROSITY_AXIS,
    SoluteDispersionCell, SoluteDispersionControls, SoluteDispersionError, SoluteDispersionReport,
    TRANSPORT_SCALE_AXIS,
};
mod sorption;
mod transport;
pub use langmuir::populations::{
    LangmuirPopulationState, LangmuirSolutePopulation, MultiLangmuirSoluteControls,
    MultiLangmuirSoluteError, MultiLangmuirSoluteStep,
};
pub use langmuir::{
    LangmuirSoluteControls, LangmuirSoluteError, LangmuirSoluteLaw, LangmuirSoluteSorption,
    LangmuirSoluteStep,
};
pub use sorption::{LinearSoluteSorption, SorbingDiluteSoluteReport, SorbingDiluteSoluteStep};
pub use transport::{
    DiluteSoluteBoundary, DiluteSoluteControls, DiluteSoluteError, DiluteSoluteReport,
    DiluteSoluteStep,
};

impl LiquidSaturationCandidate {
    /// Transport with multiple independent original-source sorbent populations
    /// in each cell/constituent. Every population retains its own occupied and
    /// vacant stocks, dry mass and equilibrium or local LDF source. Sparse
    /// bindings are grouped once; their analytic storage changes/tangents enter
    /// one original aqueous solve. No population is updated by operator splitting.
    ///
    /// Distinct IDs denote disjoint site inventories, not competing species on
    /// shared sites. Common concentration support is required; equilibrium
    /// populations must already be equilibrated. See the one-population entry
    /// for physical-area, source, passive-flow and first-order-time assumptions.
    #[allow(clippy::too_many_arguments)]
    pub fn transport_multi_langmuir_dilute_solutes<'a, const N: usize>(
        &self,
        original: &LiquidSaturationColumn,
        initial_aqueous_amount_mol_m2: &[[f64; N]],
        populations: &'a [LangmuirSolutePopulation<'a>],
        transverse_area_m2: f64,
        aqueous_diffusivity_m2_s: &[[f64; N]],
        transport_tortuosity: &[f64],
        boundaries: [DiluteSoluteBoundary<N>; 2],
        controls: MultiLangmuirSoluteControls<N>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MultiLangmuirSoluteStep<N>, MultiLangmuirSoluteError> {
        let mut storage = langmuir::populations::MultiLangmuirStorage::new(
            original.cells.len(),
            self.step.duration_s,
            populations,
            transverse_area_m2,
            controls,
            &mut cancelled,
        )?;
        let result = self.transport_solute_storage(
            original,
            initial_aqueous_amount_mol_m2,
            aqueous_diffusivity_m2_s,
            transport_tortuosity,
            boundaries,
            controls.langmuir.transport,
            &mut storage,
            None,
            &mut cancelled,
        );
        let aqueous = storage.recover(result)?;
        storage.finish(initial_aqueous_amount_mol_m2, aqueous, &mut cancelled)
    }

    /// Prepare source-qualified longitudinal mixing on this exact accepted
    /// water interval. Each adjacent cell uses its OWN endpoint liquid volume
    /// fraction and the SAME face volume/duration. Opposite faces never cancel.
    /// The original hydraulic interval freezes its chosen face flux; this does
    /// not resolve within-interval changes or reversals. Refine water intervals
    /// to resolve those histories. No grid-dependent dispersivity is inferred.
    pub fn with_solute_dispersion<'a>(
        &'a self,
        original: &'a LiquidSaturationColumn,
        cells: &'a [SoluteDispersionCell<'a>],
        controls: SoluteDispersionControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DispersiveSoluteInterval<'a>, SoluteDispersionError> {
        transport::poll(&mut cancelled)?;
        self.admit_solute_interval(original)?;
        if cells.len() != original.cells.len() {
            return Err(SoluteDispersionError::Invalid(
                "one dispersion source binding per original hydraulic cell required",
            ));
        }
        let prepared = dispersion::prepare(
            [original.elapsed_s, self.state.elapsed_s],
            self.step.duration_s,
            cells,
            controls,
            |i| {
                let cell = &original.cells[i];
                let saturation = self.state.saturations[i];
                DispersionCellState {
                    porosity: cell.porosity,
                    liquid_saturation: saturation,
                    liquid_volume_fraction: (cell.capacity_m() * saturation) / cell.thickness_m,
                }
            },
            |face| transport::integrated_flux(self.step.duration_s, self.step.face_flux_m_s[face]),
            &mut cancelled,
        )?;
        Ok(DispersiveSoluteInterval {
            hydraulic: self,
            original,
            prepared,
        })
    }

    /// Transport passive analytical-element solutes with original source-bound
    /// Langmuir storage, optionally at the independently sourced LOCAL LDF rate.
    /// Each cell/constituent slot has independent sites; there is no competition
    /// or charged-species migration. `None` denotes inert, zero-site storage.
    ///
    /// Actual cell liquid volume is the hydraulic water capacity per area times
    /// `transverse_area_m2`; it is never a freely chosen source coordinate.
    /// Retained dry mass and occupied/vacant stocks use the same bulk area basis.
    /// Equilibrium mode requires an already equilibrated original loading.
    /// Finite-rate mode permits disequilibrium and independent desorption.
    ///
    /// Bounded damped Newton reuses the original tridiagonal aqueous stencil.
    /// Every nonlinear trial checks the original constitutive sources at its
    /// actual concentration, loading, volume and prescribed temperature. A
    /// source refusal, exhausted work or unconverged residual publishes nothing.
    /// Newton globalization and first-order endpoint coupling do not establish
    /// global convergence or temporal accuracy. The numerical payload budget
    /// excludes transient original source-owner query/receipt allocations.
    #[allow(clippy::too_many_arguments)]
    pub fn transport_langmuir_dilute_solutes<'a, const N: usize>(
        &self,
        original: &LiquidSaturationColumn,
        initial_aqueous_amount_mol_m2: &[[f64; N]],
        initial_sorbed_amount_mol_m2: &'a [[f64; N]],
        initial_vacant_amount_mol_m2: &'a [[f64; N]],
        sorption: &'a [[Option<LangmuirSoluteSorption<'a>>; N]],
        transverse_area_m2: f64,
        aqueous_diffusivity_m2_s: &[[f64; N]],
        transport_tortuosity: &[f64],
        boundaries: [DiluteSoluteBoundary<N>; 2],
        controls: LangmuirSoluteControls<N>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<LangmuirSoluteStep<N>, LangmuirSoluteError> {
        let mut storage = langmuir::LangmuirStorage::new(
            original.cells.len(),
            self.step.duration_s,
            initial_sorbed_amount_mol_m2,
            initial_vacant_amount_mol_m2,
            sorption,
            transverse_area_m2,
            controls,
            &mut cancelled,
        )?;
        let result = self.transport_solute_storage(
            original,
            initial_aqueous_amount_mol_m2,
            aqueous_diffusivity_m2_s,
            transport_tortuosity,
            boundaries,
            controls.transport,
            &mut storage,
            None,
            &mut cancelled,
        );
        let aqueous = storage.recover(result)?;
        storage.finish(initial_aqueous_amount_mol_m2, aqueous, &mut cancelled)
    }

    /// Transport dilute solutes with explicit reversible linear immobile storage.
    ///
    /// Aqueous and sorbed amounts are both mol per bulk transverse m². Every
    /// law declares B [m] so b_eq=B*c, independent of the changing water volume.
    /// The finite-rate option applies db/dt=k*(B*c-b) implicitly in the SAME
    /// transport solve; equilibrium mode requires an already equilibrated
    /// original bound pool within the absolute amount tolerance.
    ///
    /// The eliminated bound storage changes only the system diagonal and RHS.
    /// Darcy advection and theta*D/tau diffusion retain the original water
    /// geometry. Total constituent equations/balances and the represented bound
    /// exchange are checked independently. The whole proposal is first-order
    /// in time; coefficients are caller-qualified and frozen over the interval.
    /// No nonlinear sites, cross-species competition, sorption heat, drying
    /// deposition or feedback into hydraulics is inferred.
    #[allow(clippy::too_many_arguments)]
    pub fn transport_sorbing_dilute_solutes<const N: usize>(
        &self,
        original: &LiquidSaturationColumn,
        initial_aqueous_amount_mol_m2: &[[f64; N]],
        initial_sorbed_amount_mol_m2: &[[f64; N]],
        sorption_laws: &[[LinearSoluteSorption; N]],
        aqueous_diffusivity_m2_s: &[[f64; N]],
        transport_tortuosity: &[f64],
        boundaries: [DiluteSoluteBoundary<N>; 2],
        controls: DiluteSoluteControls<N>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SorbingDiluteSoluteStep<N>, DiluteSoluteError> {
        let mut storage = sorption::SorptionStorage::new(
            original.cells.len(),
            self.step.duration_s,
            initial_sorbed_amount_mol_m2,
            sorption_laws,
            controls,
            &mut cancelled,
        )?;
        let aqueous = self.transport_solute_storage(
            original,
            initial_aqueous_amount_mol_m2,
            aqueous_diffusivity_m2_s,
            transport_tortuosity,
            boundaries,
            controls,
            &mut storage,
            None,
            &mut cancelled,
        )?;
        storage.finish(initial_aqueous_amount_mol_m2, aqueous, &mut cancelled)
    }

    /// Advance passive dissolved inventories over this SAME accepted water step.
    ///
    /// The candidate binds its exact original saturation/time and native cells/
    /// calibrations. Both original and final water volume must be strictly
    /// positive in every cell; deposition at a drying endpoint requires another
    /// phase/storage owner. Input amounts are mol per bulk transverse m², whereas
    /// boundary and returned concentrations are mol per liquid m³.
    ///
    /// Backward Euler uses the accepted hydraulic face integrals and the actual
    /// final water capacity. This preserves uniform concentration under changing
    /// water storage up to the separately admitted hydraulic residual. Endpoint
    /// theta*D/tau half-cell resistances add in series. An exact zero D seals
    /// molecular exchange only; the original water may still advect solute.
    ///
    /// Donor-cell hydraulic transport and diffusion share one tridiagonal
    /// M-matrix per solute. Native adjacent-pivot LU costs O(cells*species);
    /// positive storage and nonnegative sources give a positive implicit map.
    /// No concentration clipping, limiter, invented water flux or diffusion CFL
    /// is used. Physical cell equations and global constituent balances are
    /// independently checked after solving. Failure/cancellation changes neither
    /// this hydraulic candidate, the original column nor the initial solutes.
    ///
    /// This is one first-order endpoint approximation. Refine physical water
    /// intervals for temporal accuracy; an adaptive water window's cumulative
    /// face report cannot substitute for its individual accepted intervals.
    #[allow(clippy::too_many_arguments)]
    pub fn transport_dilute_solutes<const N: usize>(
        &self,
        original: &LiquidSaturationColumn,
        initial_amount_mol_m2: &[[f64; N]],
        aqueous_diffusivity_m2_s: &[[f64; N]],
        transport_tortuosity: &[f64],
        boundaries: [DiluteSoluteBoundary<N>; 2],
        controls: DiluteSoluteControls<N>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<DiluteSoluteStep<N>, DiluteSoluteError> {
        self.transport_solute_storage(
            original,
            initial_amount_mol_m2,
            aqueous_diffusivity_m2_s,
            transport_tortuosity,
            boundaries,
            controls,
            &mut (),
            None,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn transport_solute_storage<const N: usize>(
        &self,
        original: &LiquidSaturationColumn,
        initial_amount_mol_m2: &[[f64; N]],
        aqueous_diffusivity_m2_s: &[[f64; N]],
        transport_tortuosity: &[f64],
        boundaries: [DiluteSoluteBoundary<N>; 2],
        controls: DiluteSoluteControls<N>,
        storage: &mut impl transport::ImmobileStorage<N>,
        mechanical_bulk_m2_s: Option<&[[f64; 2]]>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DiluteSoluteStep<N>, DiluteSoluteError> {
        if cancelled() {
            return Err(DiluteSoluteError::Cancelled);
        }
        let n = original.cells.len();
        transport::admit_extent(n, controls)?;
        self.admit_solute_interval(original)?;
        let mut before = Vec::with_capacity(n);
        let mut after = Vec::with_capacity(n);
        let mut widths = Vec::with_capacity(n);
        for i in 0..n {
            if i % 128 == 0 && cancelled() {
                return Err(DiluteSoluteError::Cancelled);
            }
            let cell = &original.cells[i];
            before.push(cell.capacity_m() * original.saturations[i]);
            after.push(cell.capacity_m() * self.state.saturations[i]);
            widths.push(cell.thickness_m);
        }
        let mut faces = Vec::with_capacity(n + 1);
        for (i, &flux) in self.step.face_flux_m_s.iter().enumerate() {
            if i % 128 == 0 && cancelled() {
                return Err(DiluteSoluteError::Cancelled);
            }
            let volume = self.step.duration_s * flux;
            if !volume.is_finite() || (flux != 0.0 && volume == 0.0) {
                return Err(DiluteSoluteError::Invalid(
                    "original integrated water flux is unrepresentable",
                ));
            }
            faces.push(volume);
        }
        transport::transport_with_dispersion(
            [original.elapsed_s, self.state.elapsed_s],
            self.step.duration_s,
            &before,
            &after,
            &widths,
            &faces,
            initial_amount_mol_m2,
            aqueous_diffusivity_m2_s,
            transport_tortuosity,
            boundaries,
            controls,
            storage,
            mechanical_bulk_m2_s,
            cancelled,
        )
    }

    fn admit_solute_interval(
        &self,
        original: &LiquidSaturationColumn,
    ) -> Result<(), DiluteSoluteError> {
        if self.original_elapsed_s != original.elapsed_s
            || self.original_saturations != original.saturations
            || self.state.cells != original.cells
            || self.state.elapsed_s != original.elapsed_s + self.step.duration_s
        {
            return Err(DiluteSoluteError::Invalid(
                "solute transport requires the exact original hydraulic interval and material cells",
            ));
        }
        Ok(())
    }
}

/// One prepared physical interval. Each consuming method reuses the original
/// transport/storage owner; preparation arrays move into the result without a
/// deep clone. Failure leaves the borrowed hydraulic and constituent states
/// unchanged. All passive species share the declared mechanical calibration.
pub struct DispersiveSoluteInterval<'a> {
    hydraulic: &'a LiquidSaturationCandidate,
    original: &'a LiquidSaturationColumn,
    prepared: dispersion::PreparedDispersion<'a>,
}
impl<'a> DispersiveSoluteInterval<'a> {
    /// Multiple independently sourced sorbent populations with the prepared
    /// physical molecular/longitudinal dispersion and original aqueous flow.
    #[allow(clippy::too_many_arguments)]
    pub fn transport_multi_langmuir_dilute_solutes<'b, const N: usize>(
        self,
        initial_aqueous_amount_mol_m2: &[[f64; N]],
        populations: &'b [LangmuirSolutePopulation<'b>],
        transverse_area_m2: f64,
        aqueous_diffusivity_m2_s: &[[f64; N]],
        transport_tortuosity: &[f64],
        boundaries: [DiluteSoluteBoundary<N>; 2],
        mut controls: MultiLangmuirSoluteControls<N>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DispersiveSoluteStep<'a, MultiLangmuirSoluteStep<N>>, SoluteDispersionError> {
        controls.langmuir.transport = self
            .prepared
            .transport_controls(controls.langmuir.transport)?;
        let mut storage = langmuir::populations::MultiLangmuirStorage::new(
            self.original.cells.len(),
            self.hydraulic.step.duration_s,
            populations,
            transverse_area_m2,
            controls,
            &mut cancelled,
        )
        .map_err(|error| self.prepared.multi_langmuir_error(error))?;
        let result = self.hydraulic.transport_solute_storage(
            self.original,
            initial_aqueous_amount_mol_m2,
            aqueous_diffusivity_m2_s,
            transport_tortuosity,
            boundaries,
            controls.langmuir.transport,
            &mut storage,
            Some(&self.prepared.report.mechanical_bulk_m2_s),
            &mut cancelled,
        );
        let mut aqueous = storage
            .recover(result)
            .map_err(|error| self.prepared.multi_langmuir_error(error))?;
        aqueous.report.workspace_bytes += self.prepared.report.preparation_storage_bytes;
        let result = storage
            .finish(initial_aqueous_amount_mol_m2, aqueous, &mut cancelled)
            .map_err(|error| self.prepared.multi_langmuir_error(error))?;
        Ok(self.prepared.finish(result))
    }

    pub fn dispersion(&self) -> &SoluteDispersionReport {
        &self.prepared.report
    }

    #[allow(clippy::too_many_arguments)]
    pub fn transport_dilute_solutes<const N: usize>(
        self,
        initial_amount_mol_m2: &[[f64; N]],
        aqueous_diffusivity_m2_s: &[[f64; N]],
        transport_tortuosity: &[f64],
        boundaries: [DiluteSoluteBoundary<N>; 2],
        controls: DiluteSoluteControls<N>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<DispersiveSoluteStep<'a, DiluteSoluteStep<N>>, SoluteDispersionError> {
        let controls = self.prepared.transport_controls(controls)?;
        let mut step = self
            .hydraulic
            .transport_solute_storage(
                self.original,
                initial_amount_mol_m2,
                aqueous_diffusivity_m2_s,
                transport_tortuosity,
                boundaries,
                controls,
                &mut (),
                Some(&self.prepared.report.mechanical_bulk_m2_s),
                cancelled,
            )
            .map_err(|error| self.prepared.transport_error(error))?;
        step.report.workspace_bytes += self.prepared.report.preparation_storage_bytes;
        Ok(self.prepared.finish(step))
    }

    #[allow(clippy::too_many_arguments)]
    pub fn transport_sorbing_dilute_solutes<const N: usize>(
        self,
        initial_aqueous_amount_mol_m2: &[[f64; N]],
        initial_sorbed_amount_mol_m2: &[[f64; N]],
        sorption_laws: &[[LinearSoluteSorption; N]],
        aqueous_diffusivity_m2_s: &[[f64; N]],
        transport_tortuosity: &[f64],
        boundaries: [DiluteSoluteBoundary<N>; 2],
        controls: DiluteSoluteControls<N>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DispersiveSoluteStep<'a, SorbingDiluteSoluteStep<N>>, SoluteDispersionError> {
        let controls = self.prepared.transport_controls(controls)?;
        let mut storage = sorption::SorptionStorage::new(
            self.original.cells.len(),
            self.hydraulic.step.duration_s,
            initial_sorbed_amount_mol_m2,
            sorption_laws,
            controls,
            &mut cancelled,
        )
        .map_err(|error| self.prepared.transport_error(error))?;
        let mut aqueous = self
            .hydraulic
            .transport_solute_storage(
                self.original,
                initial_aqueous_amount_mol_m2,
                aqueous_diffusivity_m2_s,
                transport_tortuosity,
                boundaries,
                controls,
                &mut storage,
                Some(&self.prepared.report.mechanical_bulk_m2_s),
                &mut cancelled,
            )
            .map_err(|error| self.prepared.transport_error(error))?;
        aqueous.report.workspace_bytes += self.prepared.report.preparation_storage_bytes;
        let result = storage
            .finish(initial_aqueous_amount_mol_m2, aqueous, &mut cancelled)
            .map_err(|error| self.prepared.transport_error(error))?;
        Ok(self.prepared.finish(result))
    }

    #[allow(clippy::too_many_arguments)]
    pub fn transport_langmuir_dilute_solutes<'b, const N: usize>(
        self,
        initial_aqueous_amount_mol_m2: &[[f64; N]],
        initial_sorbed_amount_mol_m2: &'b [[f64; N]],
        initial_vacant_amount_mol_m2: &'b [[f64; N]],
        sorption: &'b [[Option<LangmuirSoluteSorption<'b>>; N]],
        transverse_area_m2: f64,
        aqueous_diffusivity_m2_s: &[[f64; N]],
        transport_tortuosity: &[f64],
        boundaries: [DiluteSoluteBoundary<N>; 2],
        mut controls: LangmuirSoluteControls<N>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DispersiveSoluteStep<'a, LangmuirSoluteStep<N>>, SoluteDispersionError> {
        controls.transport = self.prepared.transport_controls(controls.transport)?;
        let mut storage = langmuir::LangmuirStorage::new(
            self.original.cells.len(),
            self.hydraulic.step.duration_s,
            initial_sorbed_amount_mol_m2,
            initial_vacant_amount_mol_m2,
            sorption,
            transverse_area_m2,
            controls,
            &mut cancelled,
        )
        .map_err(|error| self.prepared.langmuir_error(error))?;
        let result = self.hydraulic.transport_solute_storage(
            self.original,
            initial_aqueous_amount_mol_m2,
            aqueous_diffusivity_m2_s,
            transport_tortuosity,
            boundaries,
            controls.transport,
            &mut storage,
            Some(&self.prepared.report.mechanical_bulk_m2_s),
            &mut cancelled,
        );
        let mut aqueous = storage
            .recover(result)
            .map_err(|error| self.prepared.langmuir_error(error))?;
        aqueous.report.workspace_bytes += self.prepared.report.preparation_storage_bytes;
        let result = storage
            .finish(initial_aqueous_amount_mol_m2, aqueous, &mut cancelled)
            .map_err(|error| self.prepared.langmuir_error(error))?;
        Ok(self.prepared.finish(result))
    }
}

#[cfg(test)]
mod adapter_tests;
