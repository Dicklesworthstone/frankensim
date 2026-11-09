//! Internal stage quadrature on the original joule residual and caloric chart.
use super::*;
use crate::enthalpy::entropy::{AdvectiveEntropyContribution, EntropyAudit};

/// The physical endpoint clock and the implicit quadrature weight are distinct.
/// A retained earlier stage is native heat, never an external source or a
/// replacement for the actual initial enthalpy. This ingress is crate-internal.
#[derive(Clone, Copy)]
pub(in crate::enthalpy) struct StageQuadrature<'a> {
    pub(in crate::enthalpy) implicit_dt_s: f64,
    pub(in crate::enthalpy) earlier: Option<&'a StageHeat>,
    pub(in crate::enthalpy) retain_heat: bool,
}

/// Integrated native rates at one physical stage, captured only after an
/// entropy query. Private fields prevent callers substituting external heat
/// for internal conduction or relabeling storage as a rate integral.
pub(in crate::enthalpy) struct StageHeat {
    ids: Vec<CellId>,
    pub(super) heat_j: Vec<f64>,
    pub(super) external_heat_j: f64,
    reservoir_entropy_j_k: f64,
    conduction_generation_j_k: f64,
    source_entropy_j_k: f64,
    has_unclassified_sources: bool,
    advection: Option<AdvectiveEntropyContribution>,
}

impl StageQuadrature<'_> {
    pub(super) fn validate(&self, ids: &[CellId]) -> Result<(), DepositionError> {
        if !positive(self.implicit_dt_s)
            || (self.retain_heat && self.earlier.is_some())
            || self.earlier.is_some_and(|earlier| earlier.ids != ids)
        {
            return Err(DepositionError::InvalidInput(
                "implicit stage needs positive quadrature and the original cell ordering",
            ));
        }
        Ok(())
    }

    pub(super) fn add_earlier_entropy(
        &self,
        audit: &mut EntropyAudit,
    ) -> Result<(), DepositionError> {
        if let Some(earlier) = self.earlier {
            audit.reservoir_entropy_input_j_k =
                finite(audit.reservoir_entropy_input_j_k + earlier.reservoir_entropy_j_k)?;
            audit.conduction_generation_j_k =
                finite(audit.conduction_generation_j_k + earlier.conduction_generation_j_k)?;
            audit.unclassified_source_at_cell_j_k =
                finite(audit.unclassified_source_at_cell_j_k + earlier.source_entropy_j_k)?;
            audit.has_unclassified_sources |= earlier.has_unclassified_sources;
            if let Some(advection) = earlier.advection {
                audit.advection.get_or_insert_with(Default::default).accumulate(advection)?;
            }
            // Storage is still the exact caloric integral from the ORIGINAL
            // physical initial state. The first stage is not an extra step.
            audit.refresh_temporal_defect()?;
        }
        Ok(())
    }
}

impl StageHeat {
    pub(super) fn capture(
        ids: Vec<CellId>,
        heat_j: Vec<f64>,
        external_heat_j: f64,
        entropy: EntropyAudit,
    ) -> Self {
        Self {
            ids,
            heat_j,
            external_heat_j,
            reservoir_entropy_j_k: entropy.reservoir_entropy_input_j_k,
            conduction_generation_j_k: entropy.conduction_generation_j_k,
            source_entropy_j_k: entropy.unclassified_source_at_cell_j_k,
            has_unclassified_sources: entropy.has_unclassified_sources,
            advection: entropy.advection,
        }
    }

    /// Reweight the same native stage rates; do not derive them from a rounded
    /// storage increment or transfer its algebraic defect into later forcing.
    pub(in crate::enthalpy) fn scaled(
        mut self,
        factor: f64,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DepositionError> {
        checkpoint()?;
        if !positive(factor) {
            return Err(DepositionError::InvalidInput(
                "implicit stage heat needs a positive finite weight",
            ));
        }
        for (i, heat) in self.heat_j.iter_mut().enumerate() {
            if i % 256 == 0 {
                checkpoint()?;
            }
            *heat = finite(*heat * factor)?;
        }
        self.external_heat_j = finite(self.external_heat_j * factor)?;
        self.reservoir_entropy_j_k = finite(self.reservoir_entropy_j_k * factor)?;
        self.conduction_generation_j_k = finite(self.conduction_generation_j_k * factor)?;
        self.source_entropy_j_k = finite(self.source_entropy_j_k * factor)?;
        self.advection = self.advection.map(|value| value.scaled(factor)).transpose()?;
        checkpoint()?;
        Ok(self)
    }
}

impl<'m, T> BackwardEulerTrial<'m, T> {
    /// Preserve independent source ownership when signed powers cancel.
    /// Numerical heat and entropy values stay unchanged. Entropy-free Newton
    /// trials remain entropy-free; retained native stage heat carries the same
    /// classification into later quadrature even if its endpoint sources stop.
    pub(in crate::enthalpy) fn preserve_unclassified_sources(&mut self, active: bool) {
        if let Some(audit) = &mut self.entropy {
            audit.has_unclassified_sources |= active;
        }
        if let Some(heat) = &mut self.stage_heat {
            heat.has_unclassified_sources |= active;
        }
    }

    /// Classify only the current endpoint's admitted upwind flow. Earlier
    /// native stages already carry their original classified contributions.
    /// Actual stored energy, rates, residuals and acceptance stay unchanged.
    pub(in crate::enthalpy) fn classify_endpoint_advection(
        &mut self,
        advection: AdvectiveEntropyContribution,
        source_entropy_j_k: f64,
        has_unclassified_sources: bool,
        quadrature: Option<StageQuadrature<'_>>,
    ) -> Result<(), DepositionError> {
        let audit = self.entropy.as_mut().ok_or(DepositionError::InvalidInput(
            "advective classification requires the native entropy account",
        ))?;
        let earlier = quadrature.and_then(|stage| stage.earlier);
        // Rebuild the NONadvective source contribution from its separately
        // sampled owner, not by subtracting two large combined entropy totals.
        audit.unclassified_source_at_cell_j_k = finite(
            source_entropy_j_k + earlier.map_or(0.0, |heat| heat.source_entropy_j_k),
        )?;
        audit.has_unclassified_sources = has_unclassified_sources
            || earlier.is_some_and(|heat| heat.has_unclassified_sources);
        audit.advection.get_or_insert_with(Default::default).accumulate(advection)?;
        audit.refresh_temporal_defect()?;
        if let Some(heat) = &mut self.stage_heat {
            // Retained first-stage heat has no earlier quadrature contribution.
            heat.source_entropy_j_k = source_entropy_j_k;
            heat.has_unclassified_sources = has_unclassified_sources;
            heat.advection = Some(advection);
        }
        Ok(())
    }

    pub(in crate::enthalpy) fn accept_with_stage_heat(
        self,
        cell_tolerance_j: f64,
        total_tolerance_j: f64,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(DepositionSimulation<'m>, T, EnergyAudit, StageHeat), DepositionError> {
        self.check_acceptance(cell_tolerance_j, total_tolerance_j, checkpoint)?;
        let heat = self.stage_heat.ok_or(DepositionError::InvalidInput(
            "implicit stage did not retain its native heat rates",
        ))?;
        Ok((self.next, self.payload, self.audit, heat))
    }
}
