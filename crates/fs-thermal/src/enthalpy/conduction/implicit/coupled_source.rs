//! Sparse cross-cell heat-generation feedback on the original enthalpy graph.

use super::{ImplicitConductionChain, ImplicitConductionError, LinearTopology, finite};
use crate::enthalpy::{CellId, DepositionError, DepositionSimulation, ThermalLaw};
use std::collections::BTreeSet;

/// Signed total cell powers and the complete sparse source tangent.
///
/// Powers follow the attachment's cell order, in watts INTO each original
/// cell. Derivatives follow the immutable `power_dependencies` order and are
/// dP_row/dh_column [W kg/J], holding physical time and other coordinates fixed.
/// The inherited h coordinate is each original law's declared specific stored
/// energy. Missing pattern entries assert zero derivatives; no numerical
/// differentiation, temperature surrogate or reference conversion is inferred.
#[derive(Clone, Debug, PartialEq)]
pub struct ImplicitCoupledHeatSourceState {
    pub powers_w: Vec<f64>,
    pub power_enthalpy_derivatives_w_kg_j: Vec<f64>,
}

impl ImplicitCoupledHeatSourceState {
    pub(super) fn admit(
        &self,
        cells: usize,
        derivatives: usize,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DepositionError> {
        checkpoint()?;
        if self.powers_w.len() != cells
            || self.power_enthalpy_derivatives_w_kg_j.len() != derivatives
        {
            return Err(DepositionError::InvalidInput(
                "coupled heat-source response must match its fixed cells and sparse pattern",
            ));
        }
        for (i, value) in self
            .powers_w
            .iter()
            .chain(&self.power_enthalpy_derivatives_w_kg_j)
            .enumerate()
        {
            if i % 256 == 0 {
                checkpoint()?;
            }
            if !value.is_finite() {
                return Err(DepositionError::MaterialDomain);
            }
        }
        checkpoint()?;
        Ok(())
    }
}

/// Immutable source/feedback operator on a declared set of original calorics.
///
/// A domain owner supplies powers and their full cross-cell derivatives from
/// one simultaneous trial state, for example after its stationary electrical
/// solve. Thermal consumes only the exported heat; it does not solve a circuit,
/// infer reaction inventory, or invent a corresponding upstream work debit.
/// Source identity, calibration, joint validity domains and energy ownership
/// remain with that owner. Conductive/contact/internal transfers retain their
/// existing owners rather than being disguised as generated heat.
///
/// Caloric identities and the ordered sparse pattern must remain immutable.
/// Every potentially nonzero derivative must be declared, including one that
/// currently vanishes. Every source query must enforce its actual state/time
/// domain and may neither advance physical history nor publish a partial
/// coupled state. The native thermal owner adds each returned power once.
pub trait ImplicitCoupledHeatSourceLaw {
    fn caloric_laws(&self) -> &[&dyn ThermalLaw];

    /// (Power row, enthalpy column), both indexing the attachment's cell order.
    /// Entries must be unique and in range. Order is retained by every query.
    fn power_dependencies(&self) -> &[(usize, usize)];

    fn state_at_specific_enthalpies_and_time(
        &self,
        specific_h_j_kg: &[f64],
        time_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitCoupledHeatSourceState, DepositionError>;

    /// A finite, strictly future slope change of a continuous prescribed input.
    /// Native fixed/adaptive time owners stop exactly here. Discontinuous
    /// powers still require native events or the interval-mean observer owner.
    fn next_transition_time_s(&self, _after_s: f64) -> Option<f64> {
        None
    }
}

#[derive(Clone, Copy)]
pub(super) struct SourceDependency {
    pub row: usize,
    pub column: usize,
    column_mass_kg: f64,
}

impl SourceDependency {
    /// The unknown is a cell energy increment: dh_column/dE_column=1/m_column.
    /// In R_row=Delta E_row-dt*P_row, the row's mass does not enter this term.
    pub(super) fn energy_jacobian_entry(
        self,
        dt_s: f64,
        power_enthalpy_derivative: f64,
    ) -> Result<f64, DepositionError> {
        finite((-dt_s * power_enthalpy_derivative) / self.column_mass_kg)
    }
}

fn bind_dependencies(
    cells: &[(usize, f64)],
    pattern: &[(usize, usize)],
    maximum_derivatives: usize,
    chain: bool,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<Vec<SourceDependency>, DepositionError> {
    checkpoint()?;
    if pattern.len() > maximum_derivatives {
        return Err(DepositionError::Budget("coupled heat-source derivatives"));
    }
    let mut seen = BTreeSet::new();
    let mut mapped = Vec::new();
    mapped
        .try_reserve_exact(pattern.len())
        .map_err(|_| DepositionError::Budget("coupled heat-source derivative storage"))?;
    for &(row, column) in pattern {
        checkpoint()?;
        if row >= cells.len() || column >= cells.len() || !seen.insert((row, column)) {
            return Err(DepositionError::InvalidInput(
                "coupled heat-source pattern needs unique in-range dependencies",
            ));
        }
        let global_row = cells[row].0;
        let (global_column, column_mass_kg) = cells[column];
        if chain && global_row.abs_diff(global_column) > 1 {
            return Err(DepositionError::InvalidInput(
                "non-neighbor heat-source feedback requires the implicit graph solver",
            ));
        }
        mapped.push(SourceDependency {
            row: global_row,
            column: global_column,
            column_mass_kg,
        });
    }
    checkpoint()?;
    Ok(mapped)
}

pub(super) struct CoupledHeatSource<'a> {
    pub law: &'a dyn ImplicitCoupledHeatSourceLaw,
    pub cells: Vec<CellId>,
    pub dependencies: Vec<SourceDependency>,
}

impl CoupledHeatSource<'_> {
    pub(super) fn endpoint(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitCoupledHeatSourceState, DepositionError> {
        let mut h = Vec::new();
        h.try_reserve_exact(self.cells.len())
            .map_err(|_| DepositionError::Budget("coupled heat-source state storage"))?;
        for id in &self.cells {
            checkpoint()?;
            h.push(thermal.cells[id].specific_h_j_kg);
        }
        let state =
            self.law
                .state_at_specific_enthalpies_and_time(&h, thermal.time_s, checkpoint)?;
        state.admit(self.cells.len(), self.dependencies.len(), checkpoint)?;
        Ok(state)
    }
}

impl<'a> ImplicitConductionChain<'a> {
    /// Attach one bounded coupled source owner to selected original cells.
    ///
    /// This retains every cell's mass/storage coordinate and the source's
    /// complete fixed sparse dependencies. Graphs admit nonlocal couplings;
    /// chains admit only diagonal and neighboring dependencies before any
    /// source solve, even if a forbidden derivative would presently be zero.
    /// `maximum_derivatives` bounds additional Jacobian storage and work.
    ///
    /// Powers are additive with independently authored, observed and local
    /// endpoint sources; callers must not credit the same physical loss twice.
    /// The native external-heat account integrates their sum using the original
    /// BE/SDIRK quadrature. Source entropy remains unclassified. This attachment
    /// does not change source state, material inventory or accepted history.
    pub fn with_coupled_heat_source(
        mut self,
        thermal: &DepositionSimulation<'_>,
        cells: &[CellId],
        law: &'a dyn ImplicitCoupledHeatSourceLaw,
        maximum_derivatives: usize,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, ImplicitConductionError> {
        checkpoint()?;
        if self.coupled_heat_source.is_some()
            || cells.is_empty()
            || cells.len() > self.ids.len()
            || law.caloric_laws().len() != cells.len()
        {
            return Err(ImplicitConductionError::Invalid(
                "one nonempty coupled heat source with matching original calorics is required",
            ));
        }
        self.admit(thermal, &mut checkpoint)?;
        let calorics = law.caloric_laws();
        let mut seen = BTreeSet::new();
        let mut coordinates = Vec::with_capacity(cells.len());
        for (&id, &caloric) in cells.iter().zip(calorics) {
            checkpoint()?;
            let cell = thermal
                .cells
                .get(&id)
                .ok_or(ImplicitConductionError::Invalid(
                    "coupled heat source requires active original cells",
                ))?;
            if !seen.insert(id) || !thermal.uses_material_law(cell.material_id, caloric) {
                return Err(ImplicitConductionError::Invalid(
                    "coupled heat source requires unique cells and their original caloric objects",
                ));
            }
            coordinates.push((self.index[&id], cell.mass_kg));
        }
        let dependencies = bind_dependencies(
            &coordinates,
            law.power_dependencies(),
            maximum_derivatives,
            matches!(self.topology, LinearTopology::Chain),
            &mut checkpoint,
        )?;
        self.coupled_heat_source = Some(CoupledHeatSource {
            law,
            cells: cells.to_vec(),
            dependencies,
        });
        self.endpoint(thermal, &mut checkpoint)?;
        checkpoint()?;
        Ok(self)
    }
}
