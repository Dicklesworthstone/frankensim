//! Reuse the original network's transport selections in composed residuals.
//! Store claim identities, never temperatures, coefficients or physical history.
use super::*;

/// Immutable transport-source choices tied to one borrowed material network.
///
/// Construction resolves the supplied initial state once. Every later resolve
/// evaluates those SAME claims at the actual proposed enthalpy/phase state.
/// No mutable cache or alternate material evaluator is introduced. Keep this
/// binding across callbacks (or accepted steps) when source continuity is
/// required; explicitly prepare again to permit a new policy selection.
/// Caloric-law, card, coordinate, topology and geometry checks remain with the
/// original network. Source-domain refusal never falls back to another claim.
pub struct PreparedMaterialNetwork<'n, 'a> {
    network: &'n MaterialNetwork<'a>,
    bases: BTreeMap<usize, ElectricalPropertyBasis>,
    pins: BTreeMap<CellId, MaterialPropertySelection>,
}

impl PreparedMaterialNetwork<'_, '_> {
    /// Fresh complete material snapshots and coefficients at `thermal`.
    /// The canonical network performs all evaluation and admission. This call
    /// advances no history and is reusable by circuit, radiation, contact and
    /// gas callbacks which retain their own original source results.
    pub fn resolve(
        &self, thermal: &DepositionSimulation<'_>, cancelled: impl FnMut() -> bool,
    ) -> Result<CircuitThermalCoefficients<MaterialNetworkSample>, MaterialNetworkError> {
        self.network.resolve_bound(thermal, &self.bases, Some(&self.pins), cancelled)
    }
}

impl<'a> MaterialNetwork<'a> {
    /// Prepare the legacy explicit resistivity basis for every branch.
    /// This selects sources, not a nominal material or a new thermal state.
    pub fn prepare<'n>(
        &'n self, thermal: &DepositionSimulation<'_>, mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedMaterialNetwork<'n, 'a>, MaterialNetworkError> {
        if cancelled() { return Err(DepositionError::Cancelled.into()); }
        if self.resistors.len() > MAX_MATERIAL_CIRCUIT_BRANCHES {
            return Err(MaterialNetworkError::Invalid("material resistor bindings exceed circuit capacity"));
        }
        let bases = self.resistors.iter()
            .map(|r| (r.destination.branch, ElectricalPropertyBasis::Resistivity)).collect::<Vec<_>>();
        self.prepare_with_bases(thermal, &bases, cancelled)
    }

    /// Prepare exactly one conductivity/resistivity choice per original branch.
    ///
    /// Only the initial policy selection is retained. Values, complete material
    /// snapshots and transport coefficients are re-evaluated by `resolve_bound`
    /// at every trial. Original user-supplied pins remain the selected claims.
    /// Failed preparation or cancellation returns no partially prepared object.
    pub fn prepare_with_bases<'n>(
        &'n self, thermal: &DepositionSimulation<'_>,
        electrical_bases: &[(usize, ElectricalPropertyBasis)], mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedMaterialNetwork<'n, 'a>, MaterialNetworkError> {
        if cancelled() { return Err(DepositionError::Cancelled.into()); }
        let bases = self.basis_plan(electrical_bases)?;
        let initial = self.resolve_bound(thermal, &bases, None, &mut cancelled)?;
        let selection = |material: &ResolvedMaterialStatePoint| {
            MaterialPropertySelection::PinnedByProperty(material.properties().iter()
                .map(|p| (p.requirement().name().to_owned(), p.answer().receipt.selected)).collect())
        };
        let mut pins = BTreeMap::new();
        for sample in &initial.source.cells {
            if cancelled() { return Err(DepositionError::Cancelled.into()); }
            if let Some(material) = &sample.conductivity {
                pins.insert(sample.cell, selection(material));
            }
        }
        // resolve_bound retains original binding order, including unsorted
        // branch input. Electrical-only cells have no conductivity snapshot.
        for (binding, resistor) in self.resistors.iter().zip(&initial.resistors) {
            if cancelled() { return Err(DepositionError::Cancelled.into()); }
            pins.insert(binding.destination.cell, selection(resistor.conductor.material()));
        }
        if cancelled() { return Err(DepositionError::Cancelled.into()); }
        Ok(PreparedMaterialNetwork { network: self, bases, pins })
    }
}
