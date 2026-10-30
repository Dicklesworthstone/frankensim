//! Reuse the original network's transport selections in composed residuals.
//! Store claim identities, never temperatures, coefficients or physical history.
use super::*;

mod extension;
mod expansion;
pub use expansion::{
    ConductorFreeExpansion, ExpandingBoundary, ExpandingBoundarySample,
    ExpandingMaterialNetwork, ExpandingNetworkSample,
};

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
    pub(super) network: &'n MaterialNetwork<'a>,
    pub(super) bases: BTreeMap<usize, ElectricalPropertyBasis>,
    pub(super) pins: BTreeMap<CellId, MaterialPropertySelection>,
}

impl PreparedMaterialNetwork<'_, '_> {
    /// Fresh complete material snapshots and coefficients at `thermal`.
    /// The canonical network performs all evaluation and admission. This call
    /// advances no history and is reusable by circuit, radiation, contact and
    /// gas callbacks which retain their own original source results.
    pub fn resolve(
        &self, thermal: &DepositionSimulation<'_>, cancelled: impl FnMut() -> bool,
    ) -> Result<CircuitThermalCoefficients<MaterialNetworkSample>, MaterialNetworkError> {
        self.network.resolve_bound(
            thermal,
            &self.bases,
            Some(&self.pins),
            None,
            cancelled,
        )
    }

    /// Resolve this original material assembly with additional physical links.
    ///
    /// The supplied map and the authored Fourier bindings must partition all
    /// native links exactly. Added links require zero authored conductance and
    /// finite nonnegative G; omission, overlap and inactive IDs refuse. The
    /// caller retains the original source samples producing these conductances.
    /// Native link geometry, reciprocal heat transfer and integration remain
    /// authoritative. No additional source heat is generated.
    pub fn resolve_with_link_conductances(
        &self,
        thermal: &DepositionSimulation<'_>,
        additional_links: &BTreeMap<u64, f64>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<CircuitThermalCoefficients<MaterialNetworkSample>, MaterialNetworkError> {
        self.network.resolve_bound_with_links(
            thermal, &self.bases, Some(&self.pins), None, Some(additional_links), cancelled,
        )
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
        self.prepare_bound(thermal, bases, None, cancelled)
    }

    /// Select original material claims while another physical owner supplies
    /// the remaining native links. Their coefficients must be re-resolved and
    /// supplied to `resolve_with_link_conductances` at every trial. Ordinary
    /// `resolve` still refuses incomplete Fourier coverage; it never silently
    /// turns these links into insulation or caches their preparation values.
    pub fn prepare_with_link_conductances<'n>(
        &'n self,
        thermal: &DepositionSimulation<'_>,
        electrical_bases: &[(usize, ElectricalPropertyBasis)],
        additional_links: &BTreeMap<u64, f64>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedMaterialNetwork<'n, 'a>, MaterialNetworkError> {
        if cancelled() { return Err(DepositionError::Cancelled.into()); }
        let bases = self.basis_plan(electrical_bases)?;
        self.prepare_bound_with_links(thermal, bases, None, Some(additional_links), cancelled)
    }

    // Share the original pin-selection owner with explicitly coupled electrical
    // axes. Ordinary preparation/resolution still passes no additional axes.
    pub(super) fn prepare_bound<'n>(
        &'n self,
        thermal: &DepositionSimulation<'_>,
        bases: BTreeMap<usize, ElectricalPropertyBasis>,
        electrical_axes: Option<&ElectricalAxes>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedMaterialNetwork<'n, 'a>, MaterialNetworkError> {
        self.prepare_bound_with_links(thermal, bases, electrical_axes, None, cancelled)
    }

    fn prepare_bound_with_links<'n>(
        &'n self,
        thermal: &DepositionSimulation<'_>,
        bases: BTreeMap<usize, ElectricalPropertyBasis>,
        electrical_axes: Option<&ElectricalAxes>,
        additional_links: Option<&BTreeMap<u64, f64>>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedMaterialNetwork<'n, 'a>, MaterialNetworkError> {
        let initial = self.resolve_bound_with_links(
            thermal, &bases, None, electrical_axes, additional_links, &mut cancelled,
        )?;
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

/// Prepared material transport on the original thermal observer timeline.
/// Source choices remain with the borrowed network. Only the most recent
/// accepted endpoint sample is retained; no electrical state is advanced.
#[derive(Clone)]
pub struct MaterialTransportObserver<'p, 'n, 'a> {
    prepared: &'p PreparedMaterialNetwork<'n, 'a>,
    sample: Option<MaterialNetworkSample>,
}

impl MaterialTransportObserver<'_, '_, '_> {
    /// Original material snapshots at the most recent accepted step or event.
    /// None means this observer has not yet observed an accepted thermal state.
    #[must_use]
    pub fn sample(&self) -> Option<&MaterialNetworkSample> {
        self.sample.as_ref()
    }

    fn resolve_checked(
        &self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<CircuitThermalCoefficients<MaterialNetworkSample>, MaterialNetworkError> {
        // The material evaluator accepts a cancellation predicate. Preserve
        // any richer thermal checkpoint failure instead of erasing it.
        let mut failure = None;
        let result = self.prepared.resolve(thermal, || match check() {
            Ok(()) => false,
            Err(error) => {
                failure = Some(error);
                true
            }
        });
        if let Some(error) = failure {
            return Err(error.into());
        }
        result
    }
}

impl<'n, 'a> PreparedMaterialNetwork<'n, 'a> {
    /// Drive sourced conduction together with another physical observer.
    ///
    /// This adapter is for networks with no resistor branches; use the circuit
    /// owner when electrical work and Joule heat are present. The original
    /// resolver evaluates pinned k(T, phase) on every thermal source epoch and
    /// endpoint. The existing observer pair adds sources and rejects competing
    /// conductance owners, even when they happen to provide equal coefficients.
    /// Its first member exposes material samples; its second is the caller's
    /// physical history. Both are cloned and accepted by the shared timeline.
    pub fn observe<'p, O>(
        &'p self,
        inner: O,
    ) -> Result<
        fs_thermal::observers::CoupledObservers<MaterialTransportObserver<'p, 'n, 'a>, O>,
        MaterialNetworkError,
    > {
        if !self.network.resistors.is_empty() {
            return Err(MaterialNetworkError::Invalid(
                "thermal-only observation requires no resistor branches; use the circuit owner",
            ));
        }
        Ok(fs_thermal::observers::CoupledObservers::new(
            MaterialTransportObserver {
                prepared: self,
                sample: None,
            },
            inner,
        ))
    }
}

impl fs_thermal::enthalpy::DepositionObserver for MaterialTransportObserver<'_, '_, '_> {
    type Error = MaterialNetworkError;

    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        Ok(self.resolve_checked(thermal, check)?.sources_w)
    }

    fn conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<u64, f64>, Self::Error> {
        Ok(self.resolve_checked(thermal, check)?.conductances_w_k)
    }

    fn event(
        &mut self,
        _: &fs_thermal::enthalpy::Event,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<(), Self::Error> {
        self.sample = Some(self.prepared.resolve(thermal, || false)?.source);
        Ok(())
    }

    fn step(
        &mut self,
        _: f64,
        _: f64,
        _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.sample = Some(self.resolve_checked(thermal, check)?.source);
        Ok(())
    }
}
