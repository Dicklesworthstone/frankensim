//! Continue thermal transport after explicit cell/contact topology growth.
//! Existing cells keep their material source, not just their last coefficient.
use super::*;

impl PreparedMaterialNetwork<'_, '_> {
    /// Original branch-to-cell heat routes from this admitted material network.
    /// Coupled electrical consumers use these instead of an independently
    /// authored destination list that could heat a different material. The
    /// original circuit owner still checks complete active branch coverage.
    #[must_use]
    pub fn resistor_destinations(&self) -> Vec<ResistorCell> {
        self.network.resistors.iter().map(|r| r.destination).collect()
    }

    /// Bind a complete successor thermal network while retaining old source choices.
    ///
    /// Use after accepted births/contact events. Every old cell must remain, with
    /// the same caloric object, material card, material ID and fixed conditions.
    /// Existing conductivity keys and selected claims remain authoritative even
    /// when another claim becomes eligible at the new temperature. A caller
    /// intentionally changing an old source must use ordinary `prepare` instead.
    ///
    /// New cells, and previously unrequested conductivity on old isolated cells,
    /// use the successor's explicit selection policy. Links/half-cell distances
    /// come from that complete successor description; no contact is inferred.
    /// Electrical branches refuse: this does not transfer circuit topology/state.
    /// No temperature, coefficient, event cursor or physical history is stored or
    /// changed here. Failure returns no partially prepared successor.
    pub fn extend_thermal<'n, 'a>(
        &self, network: &'n MaterialNetwork<'a>, thermal: &DepositionSimulation<'_>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedMaterialNetwork<'n, 'a>, MaterialNetworkError> {
        let invalid = MaterialNetworkError::Invalid;
        if cancelled() { return Err(DepositionError::Cancelled.into()); }
        if !self.network.resistors.is_empty() || !network.resistors.is_empty() {
            return Err(invalid("thermal extension cannot discard or introduce circuit branches"));
        }
        if network.cells.len() != thermal.cells().len() || network.links.len() != thermal.links().len()
            || self.network.cells.len() > network.cells.len() {
            return Err(invalid("extension needs complete successor cells/links and every old cell"));
        }
        let mut remaining = BTreeMap::new();
        for cell in self.network.cells {
            if cancelled() { return Err(DepositionError::Cancelled.into()); }
            remaining.insert(cell.cell, cell);
        }
        let mut staged = Vec::with_capacity(network.cells.len());
        let mut seen = std::collections::BTreeSet::new();
        for cell in network.cells {
            if cancelled() { return Err(DepositionError::Cancelled.into()); }
            if !seen.insert(cell.cell) {
                return Err(invalid("duplicate successor material cell"));
            }
            let mut conductivity = cell.conductivity.clone();
            if let Some(old) = remaining.remove(&cell.cell) {
                let a = old.point;
                let b = cell.point;
                if cell.material_id != old.material_id || !std::ptr::eq(cell.curve, old.curve)
                    || a.material.content_hash() != b.material.content_hash()
                    || a.temperature_axis != b.temperature_axis
                    || a.liquid_fraction_axis != b.liquid_fraction_axis
                    || a.context.axis_quantities() != b.context.axis_quantities()
                    || a.context.axes().len() != b.context.axes().len()
                    || a.context.axes().iter().any(|(key, value)| b.context.axes().get(key)
                        .is_none_or(|other| key.as_str() != a.temperature_axis
                            && Some(key.as_str()) != a.liquid_fraction_axis
                            && value.to_bits() != other.to_bits())) {
                    return Err(invalid("existing transport cells must retain their original material and fixed context"));
                }
                if let Some((key, _)) = &old.conductivity {
                    let (next_key, _) = conductivity.as_ref()
                        .ok_or(invalid("existing conductivity cannot be removed during extension"))?;
                    if key.name() != next_key.name() || key.quantity() != next_key.quantity() {
                        return Err(invalid("existing conductivity role cannot change during extension"));
                    }
                    conductivity = Some((*next_key, self.pins.get(&cell.cell)
                        .ok_or(invalid("existing conductivity source selection is absent"))?.clone()));
                }
            }
            staged.push(MaterialCell { cell: cell.cell, material_id: cell.material_id,
                curve: cell.curve, point: cell.point, conductivity });
        }
        if !remaining.is_empty() {
            return Err(invalid("an existing transport cell is absent from the successor"));
        }
        // The original resolver admits the combined old pins/new policies and
        // computes a fresh complete state. No substitute conductivity evaluator.
        let stage = MaterialNetwork { cells: &staged, resistors: &[], links: network.links };
        let prepared = stage.prepare(thermal, &mut cancelled)?;
        let pins = prepared.pins;
        let bases = prepared.bases;
        if cancelled() { return Err(DepositionError::Cancelled.into()); }
        // Keep only owned claim choices, never a borrow of the temporary stage.
        Ok(PreparedMaterialNetwork { network, bases, pins })
    }
}
