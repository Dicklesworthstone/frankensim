//! Mechanically owned planar heat paths on the original complete assembly.
use super::*;
use std::collections::BTreeSet;

impl PreparedHeatNetwork<'_, '_> {
    /// Resolve pinned bulk and interface laws on current scalar planar paths.
    /// The original series owner uses both current half-distances and current
    /// area, then adds contact resistance once with that same area. Reference
    /// masses, caloric charts, thermal link topology and fixed source coordinates
    /// do not change. Mechanical correspondence belongs to the supplying owner.
    ///
    /// Every path incident to an affected cell must have a current metric. This
    /// route refuses directional conductivity, integrated layer geometry and
    /// untransformed exterior/finite-coolant coverage on those cells. A disjoint
    /// original heat component retains all its sources and fixed geometry.
    pub fn resolve_heat_with_contact_metrics(
        &self, cx: &Cx<'_>, thermal: &DepositionSimulation<'_>,
        pressures: &[ThermalContactPressure<'_>], metrics: &[ThermalLinkMetric],
    ) -> Result<ContactThermalCoefficients<HeatNetworkSample>, HeatNetworkError> {
        self.resolve_with_optional_link_metrics(cx, thermal, None, pressures, metrics)
    }

    /// The same current paths with disjoint additional physical link owners.
    /// Added conductances cannot silently retain obsolete geometry on affected
    /// cells; their complete components must be disjoint from the moving paths.
    pub fn resolve_with_link_conductances_and_contact_metrics(
        &self, cx: &Cx<'_>, thermal: &DepositionSimulation<'_>,
        additional_links: &BTreeMap<u64, f64>, pressures: &[ThermalContactPressure<'_>],
        metrics: &[ThermalLinkMetric],
    ) -> Result<ContactThermalCoefficients<HeatNetworkSample>, HeatNetworkError> {
        self.resolve_with_optional_link_metrics(cx, thermal, Some(additional_links), pressures, metrics)
    }

    pub(super) fn admit_contact_metrics(
        &self, cx: &Cx<'_>, thermal: &DepositionSimulation<'_>,
        additional_links: Option<&BTreeMap<u64, f64>>, metrics: &[ThermalLinkMetric],
    ) -> Result<(), HeatNetworkError> {
        check(cx)?;
        if metrics.is_empty() { return Ok(()); }
        if metrics.len() > thermal.links().len() {
            return Err(HeatNetworkError::Invalid("current thermal metrics exceed active links"));
        }
        let mut selected = BTreeSet::new();
        let mut cells = BTreeSet::new();
        for metric in metrics {
            check(cx)?;
            metric.validate()?;
            let link = thermal.links().get(&metric.link)
                .ok_or(HeatNetworkError::Invalid("current thermal metric names an inactive link"))?;
            if !selected.insert(metric.link) {
                return Err(HeatNetworkError::Invalid("duplicate current thermal path metric"));
            }
            cells.extend([link.a, link.b]);
        }
        if let Some(integrated) = &self.integrated_conduction {
            for binding in integrated.model.bindings() {
                check(cx)?;
                if selected.contains(&binding.link) {
                    return Err(HeatNetworkError::Invalid(
                        "current planar metrics cannot replace integrated layer geometry"));
                }
            }
        }
        for source in &self.radiating_cells {
            check(cx)?;
            if cells.contains(&source.cell) {
                return Err(HeatNetworkError::Invalid(
                    "deformed radiating cells require current surface and enclosure geometry"));
            }
        }
        if let Some(additional) = additional_links {
            for id in additional.keys() {
                check(cx)?;
                let link = thermal.links().get(id)
                    .ok_or(HeatNetworkError::Invalid("additional thermal link is inactive"))?;
                if cells.contains(&link.a) || cells.contains(&link.b) {
                    return Err(HeatNetworkError::Invalid(
                        "deformed thermal cells require current additional-link geometry"));
                }
            }
        }
        check(cx)?;
        Ok(())
    }
}
