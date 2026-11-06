//! Simultaneous endpoint heat transport and selected liquid/vapor exchanges.
//!
//! The original backward-Euler owner evaluates every cell, contact and thermal
//! boundary at one endpoint. For each evaporating cell its storage equation is
//! `m0*(h1-h0) + dm*L(T1) = dt*P1`. The retained mass `m0` is a coordinate for
//! this reduced open-system balance, not the proposed endpoint inventory.
//! Acceptance installs `m1 = m0-dm` and accounts for the full vapor enthalpy
//! `dm*(h1+L1)` once. There is no sensible-only intermediate temperature.
//!
//! Low-level residuals and bounded native nonlinear solves share one physical
//! owner. Endpoint convergence does not establish time accuracy; the adaptive
//! march controls physical h/T/mass error separately. The film remains prescribed
//! and its receiving mass/energy account belongs to the caller.
//! Geometry and exposure stay fixed over the step;
//! endpoint events remain pending until all coupled owners accept the interval.

use super::ThermalEvaporationError;
use crate::enthalpy::backward_euler::BoundaryExchangeStepCoefficients;
use crate::enthalpy::{CellId, DepositionError, DepositionSimulation};
use crate::evaporation::{EvaporationLaw, GasFilm, LiquidSurfaceState, SurfaceStep};

/// Several phase surfaces sharing the same original heat graph and acceptance.
#[path = "implicit/network.rs"]
pub mod network;
pub use network::{
    EvaporatingCellBinding, EvaporatingCellSource, EvaporatingNetworkBackwardEulerAudit,
    EvaporatingNetworkBackwardEulerStep, EvaporatingNetworkBackwardEulerTrial,
    EvaporatingNetworkVaporEnthalpyAudit, EvaporatingNetworkVaporEnthalpyTrial,
    VaporEnthalpyStepCoefficients,
};

/// Bounded automatic nonlinear solves of this original endpoint residual.
#[path = "implicit/solve.rs"]
pub mod solve;

#[path = "implicit/timeline.rs"]
pub mod timeline;
pub use solve::{
    EvaporationImplicitControls, EvaporationImplicitError, EvaporationImplicitReport,
    EvaporationNetworkImplicitControls, EvaporationNetworkImplicitReport,
};
pub use timeline::{EvaporationImplicitTimelineControls, EvaporationImplicitTimelineReport};

/// Conservative multi-surface advancement with physical h/T/mass time accuracy.
#[path = "implicit/adaptive.rs"]
pub mod adaptive;
pub use adaptive::{
    AdaptiveEvaporationControls, AdaptiveEvaporationError, AdaptiveEvaporationReport,
    EvaporationSurfaceTotals, EvaporationTimeError,
};

/// Physical open-system account for one interval, before endpoint events.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EvaporatingBackwardEulerAudit {
    /// Original sources and prescribed thermal baths, positive inward [J].
    /// Internal contacts/transfers and phase enthalpy are excluded.
    pub external_heat_j: f64,
    /// Actual liquid/network storage change, including the changed mass [J].
    pub stored_energy_change_j: f64,
    /// Final network mass plus vapor transfer minus initial network mass [kg].
    pub mass_balance_residual_kg: f64,
    /// Storage change plus full vapor enthalpy minus external heat [J].
    /// Evaluated after analytically cancelling the common carried-liquid datum;
    /// summing rounded report fields can lose accuracy at large datum offsets.
    pub energy_residual_j: f64,
    /// Local phase receipt. Its heat includes the cell's signed contact heat;
    /// `external_heat_j` above excludes those conservative internal exchanges.
    pub surface: SurfaceStep,
}

impl EvaporatingBackwardEulerAudit {
    pub(super) fn from_network(audit: &EvaporatingNetworkBackwardEulerAudit, id: CellId) -> Self {
        Self {
            external_heat_j: audit.external_heat_j,
            stored_energy_change_j: audit.stored_energy_change_j,
            mass_balance_residual_kg: audit.mass_balance_residual_kg,
            energy_residual_j: audit.energy_residual_j,
            surface: audit.surfaces[&id],
        }
    }
}

/// An event-free step on the original thermal topology and caloric laws.
/// Unknowns are `m0*(h1-h0)` in `cell_ids()` order, including every active cell.
pub struct EvaporatingBackwardEulerStep<'s, 'm, 'l, L: EvaporationLaw + ?Sized> {
    inner: EvaporatingNetworkBackwardEulerStep<'s, 'm, 'l, L>,
    id: CellId,
}

impl<'s, 'm, 'l, L: EvaporationLaw + ?Sized> EvaporatingBackwardEulerStep<'s, 'm, 'l, L> {
    /// Prepare one positive-duration interval after processing start-time events.
    /// The caller explicitly admits the pure-liquid caloric interval. A missing
    /// liquid cell, initially boiling source or unsupported original law refuses.
    /// Events strictly inside the interval refuse; endpoint events stay pending.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        initial: &'s DepositionSimulation<'m>,
        id: CellId,
        law: &'l L,
        film: GasFilm,
        target_s: f64,
        pure_liquid_enthalpy_range_j_kg: [f64; 2],
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, ThermalEvaporationError> {
        let binding = EvaporatingCellBinding {
            id,
            law,
            film,
            pure_liquid_enthalpy_range_j_kg,
        };
        let inner =
            EvaporatingNetworkBackwardEulerStep::new(initial, &[binding], target_s, checkpoint)?;
        Ok(Self { inner, id })
    }

    /// Stable ordering shared by energy coordinates and local joule residuals.
    pub fn cell_ids(&self) -> &[CellId] {
        self.inner.cell_ids()
    }

    pub const fn duration_s(&self) -> f64 {
        self.inner.duration_s()
    }

    pub const fn target_s(&self) -> f64 {
        self.inner.target_s()
    }

    /// Use the original owner's representability projection for solver iterates.
    /// It does not discard any heat from the residual or weaken acceptance.
    pub fn resolvable_increments_j(&self, energy_increments_j: &[f64]) -> Vec<f64> {
        self.inner.resolvable_increments_j(energy_increments_j)
    }

    /// Evaluate one tentative endpoint, with no accepted history changes.
    ///
    /// `resolve` receives original endpoint calorics/topology with the INITIAL
    /// mass retained for BE coordinates, plus the explicit proposed liquid mass
    /// and enthalpy. Mass-dependent sources must use that proposed liquid state.
    /// Return only independently owned heat: the adapter adds phase transfer.
    /// Source, link, bath and conservative-transfer coefficients retain the
    /// meanings of the original `BoundaryExchangeStepCoefficients` ingress.
    ///
    /// Each call performs one interface-property query and two original thermal
    /// rate evaluations: one independent physical heat account and one native
    /// residual including latent exchange. Work is bounded by active topology;
    /// the caller must bound the number of nonlinear evaluations.
    pub fn evaluate<T, E>(
        &self,
        energy_increments_j: &[f64],
        resolve: impl FnOnce(
            &DepositionSimulation<'m>,
            LiquidSurfaceState,
        ) -> Result<BoundaryExchangeStepCoefficients<T>, E>,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<EvaporatingBackwardEulerTrial<'m, T>, E>
    where
        E: From<ThermalEvaporationError> + From<DepositionError>,
    {
        let inner = self.inner.evaluate(
            energy_increments_j,
            |thermal, liquids| resolve(thermal, liquids[&self.id]),
            checkpoint,
        )?;
        let audit = EvaporatingBackwardEulerAudit::from_network(inner.audit(), self.id);
        Ok(EvaporatingBackwardEulerTrial { inner, audit })
    }
}

/// A complete tentative endpoint. No mass or history is published by evaluation.
pub struct EvaporatingBackwardEulerTrial<'m, T> {
    inner: EvaporatingNetworkBackwardEulerTrial<'m, T>,
    audit: EvaporatingBackwardEulerAudit,
}

impl<'m, T> EvaporatingBackwardEulerTrial<'m, T> {
    /// Original local energy residuals, including represented latent transfer.
    pub fn residual_j(&self) -> &[f64] {
        self.inner.residual_j()
    }

    /// Includes the independent physical surface balance, so cancellation in
    /// mixed source arithmetic cannot hide a local defect.
    pub fn maximum_cell_residual_j(&self) -> f64 {
        self.inner.maximum_cell_residual_j()
    }

    pub const fn audit(&self) -> &EvaporatingBackwardEulerAudit {
        &self.audit
    }

    pub fn payload(&self) -> &T {
        self.inner.payload()
    }

    /// Accept only after native local/total energy gates, the independent open
    /// energy account, represented-rate and whole-network mass gates all pass.
    /// Tolerances are finite nonnegative absolute joules and kilograms. Events
    /// at the target remain pending; reference geometry and authored inputs are
    /// unchanged. Algebraic acceptance does not establish timestep accuracy.
    pub fn accept(
        self,
        cell_tolerance_j: f64,
        total_tolerance_j: f64,
        mass_tolerance_kg: f64,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(DepositionSimulation<'m>, T, EvaporatingBackwardEulerAudit), ThermalEvaporationError>
    {
        let (next, payload, _) = self.inner.accept(
            cell_tolerance_j,
            total_tolerance_j,
            mass_tolerance_kg,
            checkpoint,
        )?;
        Ok((next, payload, self.audit))
    }
}

#[cfg(test)]
#[path = "implicit/tests.rs"]
mod tests;
