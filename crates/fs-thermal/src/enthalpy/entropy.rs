//! Thermal entropy from the same caloric states and conductive interfaces.
//!
//! Birth carries a new material cell into the system. Its incoming entropy is
//! excluded on BOTH sides of this relative account; subsequent thermal aging
//! is integrated from its actual birth enthalpy. No absolute entropy datum is
//! needed. This account does not supply chemical, advective or mechanical entropy.

use super::*;

/// Entropy over accepted thermal intervals, in J/K. This is an accounting
/// identity, not a spatial/time accuracy certificate or a global second-law proof.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct EntropyAudit {
    /// Actual sum of m*integral(dh/T), excluding entropy delivered at cell birth.
    pub thermal_storage_change_j_k: f64,
    /// Signed input Q/T_reservoir at each actual ambient/support temperature.
    pub reservoir_entropy_input_j_k: f64,
    /// Nonnegative production in finite-temperature conductive links and boundaries,
    /// evaluated on the SAME start/end quadrature as the accepted heat flux.
    pub conduction_generation_j_k: f64,
    /// Applied source energy divided by its receiving cell's quadrature temperature.
    /// This is a thermal-storage contribution, NOT an assigned boundary entropy.
    /// The source owner must classify heat, work conversion and its own entropy.
    pub unclassified_source_at_cell_j_k: f64,
    /// Whether ANY active source was nonzero, even if signed contributions cancel.
    pub has_unclassified_sources: bool,
    /// Storage minus reservoir input, conductive production and source-at-cell
    /// contribution. Includes time quadrature, algebraic residual and roundoff;
    /// never relabel this as physical heat or material entropy production.
    pub temporal_defect_j_k: f64,
}

impl EntropyAudit {
    fn accumulate(&mut self, other: Self) -> Result<(), DepositionError> {
        self.thermal_storage_change_j_k =
            finite(self.thermal_storage_change_j_k + other.thermal_storage_change_j_k)?;
        self.reservoir_entropy_input_j_k =
            finite(self.reservoir_entropy_input_j_k + other.reservoir_entropy_input_j_k)?;
        self.conduction_generation_j_k =
            finite(self.conduction_generation_j_k + other.conduction_generation_j_k)?;
        self.unclassified_source_at_cell_j_k =
            finite(self.unclassified_source_at_cell_j_k + other.unclassified_source_at_cell_j_k)?;
        self.has_unclassified_sources |= other.has_unclassified_sources;
        self.temporal_defect_j_k = finite(
            self.thermal_storage_change_j_k
                - self.reservoir_entropy_input_j_k
                - self.conduction_generation_j_k
                - self.unclassified_source_at_cell_j_k,
        )?;
        Ok(())
    }
}

// Fixed-order compensated sums retain small cell/interface terms without
// changing the original enthalpy update or its deterministic reduction order.
#[derive(Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) -> Result<(), DepositionError> {
        let value = finite(value)?;
        let next = finite(self.value + value)?;
        let lost = if self.value.abs() >= value.abs() {
            (self.value - next) + value
        } else {
            (value - next) + self.value
        };
        self.correction = finite(self.correction + lost)?;
        self.value = next;
        Ok(())
    }
    fn total(self) -> Result<f64, DepositionError> {
        finite(self.value + self.correction)
    }
}

impl DepositionSimulation<'_> {
    /// Advance the original explicit thermal owner and integrate entropy on its
    /// actual CFL/event grid. Every active law must provide its exact dh/T integral.
    /// A missing law, domain failure or cancellation publishes no candidate.
    /// Existing `advance_to` performs no entropy evaluation and keeps its behavior.
    ///
    /// Conductive production uses starting temperatures; the signed temporal
    /// defect exposes the difference from exact state entropy. Source watts are
    /// deliberately unclassified: their entropy cannot be inferred from watts
    /// alone. Cost remains O(cells + links) storage and sparse work per substep.
    pub fn advance_to_with_entropy(
        &self,
        target_s: f64,
        config: StepConfig,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, EnergyAudit, EntropyAudit), DepositionError> {
        let (next, observer, energy) =
            self.advance_observed(&EntropyObserver::default(), target_s, config, checkpoint)?;
        Ok((next, energy, observer.audit))
    }

    // `self` is the candidate AFTER changing enthalpy, but with the SAME
    // coefficients/geometry that supplied heat. Temperatures are supplied by
    // the original rate owner at its quadrature point, not re-inferred here.
    pub(super) fn entropy_step_audit(
        &self,
        initial_h: &BTreeMap<CellId, f64>,
        temperatures: &BTreeMap<CellId, f64>,
        duration_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<EntropyAudit, DepositionError> {
        checkpoint()?;
        if initial_h.len() != self.cells.len()
            || temperatures.len() != self.cells.len()
            || !positive(duration_s)
        {
            return Err(DepositionError::InvalidInput(
                "entropy interval must retain its active cells",
            ));
        }
        let mut storage = Sum::default();
        let mut boundary = Sum::default();
        let mut generation = Sum::default();
        let mut source = Sum::default();
        let mut has_source = false;
        let mut covered = BTreeMap::new();
        for (i, (&id, cell)) in self.cells.iter().enumerate() {
            if i % 256 == 0 {
                checkpoint()?;
            }
            let old_h = *initial_h.get(&id).ok_or(DepositionError::Arithmetic)?;
            let ds = self.materials[&cell.material_id]
                .specific_entropy_change_j_kg_k(old_h, cell.specific_h_j_kg)?
                .ok_or(DepositionError::InvalidInput(
                    "thermal law does not provide entropy change",
                ))?;
            storage.add(finite(cell.mass_kg * finite(ds)?)?)?;
            covered.insert(id, cell.boundary.reservoir.map_or(0.0, |r| r.area_m2));
        }
        for (i, link) in self.links.values().enumerate() {
            if i % 256 == 0 {
                checkpoint()?;
            }
            let ta = temperatures[&link.a];
            let tb = temperatures[&link.b];
            let power = finite(link.conductance_w_k * (tb - ta))?;
            // q*(1/Ta - 1/Tb), evaluated without subtracting reciprocals.
            generation.add(finite((power / ta) * ((tb - ta) / tb))?)?;
            let areas = self.contact_areas(link);
            for (id, area) in [(link.a, areas[0]), (link.b, areas[1])] {
                let value = covered.get_mut(&id).ok_or(DepositionError::Arithmetic)?;
                *value = finite(*value + area)?;
            }
        }
        for (i, (&id, cell)) in self.cells.iter().enumerate() {
            if i % 256 == 0 {
                checkpoint()?;
            }
            let t = temperatures[&id];
            let area = cell.geometry.surface_area_m2 - covered[&id];
            if !positive(t) || !nonnegative(area) {
                return Err(DepositionError::InvalidInput(
                    "invalid entropy temperature or exposed area",
                ));
            }
            let ambient = cell.boundary.ambient_k;
            let g = finite(area * cell.boundary.htc_w_m2_k)?;
            let q = finite(g * (ambient - t))?;
            boundary.add(finite(q / ambient)?)?;
            generation.add(finite((q / t) * ((ambient - t) / ambient))?)?;
            if let Some(reservoir) = cell.boundary.reservoir {
                let tr = reservoir.temperature_k;
                let q = finite(reservoir.conductance_w_k * (tr - t))?;
                boundary.add(finite(q / tr)?)?;
                generation.add(finite((q / t) * ((tr - t) / tr))?)?;
            }
            source.add(finite(cell.boundary.source_w / t)?)?;
            has_source |= cell.boundary.source_w != 0.0;
        }
        let mut audit = EntropyAudit {
            thermal_storage_change_j_k: storage.total()?,
            reservoir_entropy_input_j_k: finite(duration_s * boundary.total()?)?,
            conduction_generation_j_k: finite(duration_s * generation.total()?)?,
            unclassified_source_at_cell_j_k: finite(duration_s * source.total()?)?,
            has_unclassified_sources: has_source,
            temporal_defect_j_k: 0.0,
        };
        audit.temporal_defect_j_k = finite(
            audit.thermal_storage_change_j_k
                - audit.reservoir_entropy_input_j_k
                - audit.conduction_generation_j_k
                - audit.unclassified_source_at_cell_j_k,
        )?;
        checkpoint()?;
        Ok(audit)
    }
}

#[derive(Clone, Default)]
struct EntropyObserver {
    initial_h: BTreeMap<CellId, f64>,
    audit: EntropyAudit,
}
impl DepositionObserver for EntropyObserver {
    type Error = DepositionError;
    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.initial_h.clear();
        for (i, (&id, cell)) in thermal.cells.iter().enumerate() {
            if i % 256 == 0 {
                checkpoint()?;
            }
            self.initial_h.insert(id, cell.specific_h_j_kg);
        }
        Ok(BTreeMap::new())
    }
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Ok(())
    }
    fn step(
        &mut self,
        _: f64,
        duration_s: f64,
        temperatures: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.audit.accumulate(thermal.entropy_step_audit(
            &self.initial_h,
            temperatures,
            duration_s,
            checkpoint,
        )?)
    }
}

#[cfg(test)]
#[path = "entropy/tests.rs"]
mod tests;
