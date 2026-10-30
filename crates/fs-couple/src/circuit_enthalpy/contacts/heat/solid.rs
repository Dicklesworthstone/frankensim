//! Source-consistent, closed-energy admission for the original loaded 3D solid.
//! The native candidate and observer histories remain the sole state owners.

mod dynamics;
use super::{HeatNetworkError, MaterialHeatObserver, PreparedHeatNetwork};
use crate::solid_enthalpy::plastic::{
    ThermoplasticCandidate, ThermoplasticConfig, ThermoplasticError, ThermoplasticSolid,
};
use fs_exec::Cx;
use fs_material::phase::PhaseStateError;
use fs_material::thermal::source::ThermalMaterialError;
use fs_solid::plastic3::IncrementSettings;
use fs_thermal::enthalpy::{DepositionObserver, DepositionSimulation, StepConfig};
use fs_thermal::observers::{CoupledObservers, ObserverCouplingError};

/// Original lower-level failures, including which observer refused the interval.
#[derive(Debug)]
pub enum SolidHeatError<E> {
    Heat(HeatNetworkError),
    Binding(ThermoplasticError),
    Phase(PhaseStateError),
    State(ThermalMaterialError),
    Coupled(ThermoplasticError<ObserverCouplingError<HeatNetworkError, E>>),
    Invalid(&'static str),
    Cancelled,
}
impl<E: core::fmt::Debug> core::fmt::Display for SolidHeatError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material heat/3D solid: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for SolidHeatError<E> {}

/// Closed solid/caloric account; original component reports remain unchanged.
/// These are numerical balances, not empirical validation or time-error bounds.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SolidHeatAudit {
    /// Signed integrated internal radiative roundoff from actual source epochs.
    pub closed_radiation_net_j: f64,
    /// Authored external heat, excluding plastic/reversible heat AND this net.
    pub external_heat_j: f64,
    /// Native caloric residual with internal radiation excluded from supply.
    pub thermal_residual_j: f64,
    /// Native total work/storage residual with the same correction, not 2*tol.
    pub combined_residual_j: f64,
}

pub(crate) fn closed_account(external: f64, thermal: f64, combined: f64, radiation: f64,
    tolerance: f64) -> Result<SolidHeatAudit, &'static str> {
    if !tolerance.is_finite() || tolerance <= 0.0
        || [external, thermal, combined, radiation].iter().any(|x| !x.is_finite()) {
        return Err("finite work/heat accounts and a positive energy tolerance required");
    }
    let audit = SolidHeatAudit {
        closed_radiation_net_j: radiation,
        external_heat_j: external - radiation,
        thermal_residual_j: thermal + radiation,
        combined_residual_j: combined + radiation,
    };
    if !audit.external_heat_j.is_finite() || !audit.thermal_residual_j.is_finite()
        || !audit.combined_residual_j.is_finite()
        || audit.thermal_residual_j.abs() > tolerance
        || audit.combined_residual_j.abs() > tolerance {
        return Err("internal radiative error cannot masquerade as external solid energy");
    }
    Ok(audit)
}

impl<'m> ThermoplasticSolid<'m> {
    /// Advance the original loaded tetrahedra with a complete material heat network.
    ///
    /// Conductivity, contact resistance and emissivity use the shared network;
    /// plastic history, compatibility and physical heating use `advance` unchanged.
    /// Require identical structural/bulk cards, caloric objects and fixed source
    /// conditions BEFORE solving. Merely sharing a card is insufficient when,
    /// for example, its pressure or processing coordinates differ.
    ///
    /// The original report includes generated radiation in its thermal source
    /// account. The returned additional audit removes only its integrated NET
    /// roundoff from external supply and checks BOTH corrected residuals against
    /// the SAME original joule tolerance. No state or source is advanced again.
    /// Failure returns neither candidate history nor an admitted closed audit.
    ///
    /// `inner` retains independent observer history; its heat/work semantics
    /// remain the caller's contract. Electrical networks refuse. Radiation uses
    /// explicit thermal substeps, requiring timestep refinement; geometry and
    /// the native small-strain/single-solid restrictions are unchanged.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn advance_heat_network_candidate<'t, 'p, 'n, 'a, 'scope, O: DepositionObserver>(
        &self, thermal: &DepositionSimulation<'t>,
        prepared: &'p PreparedHeatNetwork<'n, 'a>, inner: &O,
        target_s: f64, prescribed_m: &[(usize, f64)], nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings, thermal_config: StepConfig, config: ThermoplasticConfig,
        cx: &'p Cx<'scope>,
    ) -> Result<(ThermoplasticCandidate<'m, 't,
        CoupledObservers<MaterialHeatObserver<'p, 'n, 'a, 'scope>, O>>, SolidHeatAudit),
        SolidHeatError<O::Error>> {
        cx.checkpoint().map_err(|_| SolidHeatError::Cancelled)?;
        let observer = prepared.observe(cx, inner.clone(), config.energy_tolerance_j)
            .map_err(SolidHeatError::Heat)?;
        self.admit_heat_network(thermal, prepared, cx)?;
        let candidate = self.advance(thermal, &observer, target_s, prescribed_m, nodal_forces_n,
            settings, thermal_config, config, cx).map_err(SolidHeatError::Coupled)?;
        let report = candidate.report();
        let audit = closed_account(report.external_heat_j, report.thermal_residual_j,
            report.combined_residual_j, candidate.observer().first.closed_radiation_net_j(),
            config.energy_tolerance_j).map_err(SolidHeatError::Invalid)?;
        cx.checkpoint().map_err(|_| SolidHeatError::Cancelled)?;
        Ok((candidate, audit))
    }

    // Both unpowered and finite-circuit solids must describe the same specimen.
    // This admission performs no physical advance and does not discard circuit work.
    pub(crate) fn admit_heat_network<E>(
        &self, thermal: &DepositionSimulation<'_>, prepared: &PreparedHeatNetwork<'_, '_>,
        cx: &Cx<'_>,
    ) -> Result<(), SolidHeatError<E>> {
        cx.checkpoint().map_err(|_| SolidHeatError::Cancelled)?;
        // Reuse the native material-to-cell lowering only for identity/context;
        // conductivity selections and evaluations stay with the prepared network.
        let roles = vec![None; self.materials().len()];
        let expected = self.thermal_material_cells(&roles).map_err(SolidHeatError::Binding)?;
        for cell in &expected {
            cx.checkpoint().map_err(|_| SolidHeatError::Cancelled)?;
            let bound = prepared.network.electrical.bulk.cells.iter().find(|x| x.cell == cell.cell)
                .ok_or(SolidHeatError::Invalid("structural cell lacks a heat-network binding"))?;
            if cell.material_id != bound.material_id || !std::ptr::eq(cell.curve, bound.curve)
                || cell.point.material.content_hash() != bound.point.material.content_hash()
                || cell.point.temperature_axis != bound.point.temperature_axis
                || cell.point.liquid_fraction_axis != bound.point.liquid_fraction_axis {
                return Err(SolidHeatError::Invalid("heat and solid must retain the same cell material and caloric owner"));
            }
            let state = thermal.cells().get(&cell.cell)
                .ok_or(SolidHeatError::Invalid("inactive structural thermal cell"))?;
            let phase = cell.curve.state_at_specific_enthalpy(state.specific_h_j_kg)
                .map_err(SolidHeatError::Phase)?;
            let structural = cell.point.at(phase).map_err(SolidHeatError::State)?;
            let transport = bound.point.at(phase).map_err(SolidHeatError::State)?;
            if structural != transport {
                return Err(SolidHeatError::Invalid("heat and solid material conditions disagree"));
            }
        }
        cx.checkpoint().map_err(|_| SolidHeatError::Cancelled)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn internal_radiation_cannot_hide_an_over_budget_closed_balance() {
        // Each old gate would pass separately; their closed balance must refuse.
        assert!(closed_account(8e-9, 5e-9, 5e-9, 8e-9, 1e-8).is_err());
        assert!(closed_account(-8e-9, -5e-9, -5e-9, -8e-9, 1e-8).is_err());
        let actual = closed_account(2.0 + 1e-9, -1e-9, -1e-9, 1e-9, 1e-8).unwrap();
        assert!((actual.external_heat_j - 2.0).abs() < 1e-15);
        assert_eq!(actual.thermal_residual_j, 0.0);
        assert_eq!(actual.combined_residual_j, 0.0);
    }
}
