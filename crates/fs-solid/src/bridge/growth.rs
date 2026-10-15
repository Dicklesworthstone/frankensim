//! Growing clamped overhangs, reusing the existing section and beam evolution.
//!
//! New cells continue the accepted tip tangent with zero curvature and axial
//! strain. Their reference lengths and incoming masses are explicitly supplied;
//! existing cells, branch histories and geometry are untouched at birth. The
//! next positive-time solve applies the enlarged load field to the whole road.
//! Birth itself is not a new static equilibrium or an impact calculation.
//!
//! This bounded lane starts from a born, single-layer, slender clamped seed.
//! Cross-section and layer identity continue unchanged; each new cell may have
//! its own immutable material and birth temperature. Perfect end bonding is an
//! explicit assumption, not a weld predicted from temperature or proximity.
//! Small-rotation geometry neglects second-order drawing length. No advancing
//! nozzle constraint, drawing tension/slack, right-end anchoring or free jet.

use super::{Bridge, Cell, Config, Error, Support, finite, poll, positive};
use crate::deposition::{self, LayeredStrip};
use fs_material::phase::ResolvedHeatCapacityCurve;
use fs_material::thermovisco::ThermoMaxwell;
use fs_material::thermovisco::source::ResolvedThermoMaxwell;

/// A newly delivered cell, in order along the advancing reference x direction.
#[derive(Debug, Clone, Copy)]
pub struct CellBirth<'a> {
    pub id: u64,
    /// Unstrained projected length in the small-rotation model, m.
    pub length_m: f64,
    /// Delivered reference mass, kg; never inferred from changing density.
    pub mass_kg: f64,
    pub temperature_k: f64,
    pub law: &'a ThermoMaxwell,
}

impl<'a> CellBirth<'a> {
    /// Construct one geometric birth from a single source-compatible material
    /// bundle. Reference mass is `rho_ref * area * length`; the mechanical law
    /// comes from the same resolved source. This deliberately differs from the
    /// ordinary constructor fields, which continue to support process schedules
    /// where delivered mass is independently measured and need not equal a dense
    /// equilibrium-material volume.
    ///
    /// The subsequent bridge append still validates section geometry, law domain,
    /// work budgets and the actual birth temperature. No porosity, bead swelling
    /// or process-specific mass correction is inferred here.
    pub fn from_material_bundle(id: u64, length_m: f64, area_m2: f64,
        temperature_k: f64, source: &'a ResolvedThermoMaxwell,
        caloric: &ResolvedHeatCapacityCurve) -> Result<Self, Error> {
        positive(length_m)?;
        positive(area_m2)?;
        positive(temperature_k)?;
        let density = source.reference_density_kg_m3(caloric)
            .map_err(|_| Error::Invalid("growth law and reference density must share one admitted material source"))?;
        let mass_kg = finite(finite(density * area_m2)? * length_m)?;
        positive(mass_kg)?;
        Ok(Self { id, length_m, mass_kg, temperature_k, law: source.law() })
    }
}

/// Birth accounts only. They are not a force balance or a thermal-energy audit.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BirthReport {
    pub time_s: f64,
    pub added_cells: usize,
    pub incoming_mass_kg: f64,
    pub incoming_reference_volume_m3: f64,
    /// Actual reference-line tip before and after birth, world [x,z] in metres.
    pub previous_tip_m: [f64; 2],
    pub tip_m: [f64; 2],
}

impl<'a> Bridge<'a> {
    /// Append material at this accepted time; no thermal or mechanical aging.
    /// Every newborn is stress-free on the straight continuation of the CURRENT
    /// tip. Old stress, relaxation variables, positions and birth data are copied
    /// unchanged. The supplied mass contributes weight only after its birth.
    ///
    /// Exactly one common rectangular layer must already span the clamped seed.
    /// Width/depth and layer ID are retained, but material data may differ along
    /// x. Multi-layer extension and moving the right pin of a span are refused.
    /// An empty batch is a checked no-op. Clones are bounded by the bridge budget;
    /// validation polls per existing/new cell and immediately before publication.
    pub fn append_cantilever_controlled(&self, births: &[CellBirth<'a>], config: Config,
        cancelled: &mut impl FnMut() -> bool) -> Result<(Self, BirthReport), Error> {
        poll(cancelled)?;
        if self.support != Support::CantileverLeft {
            return Err(Error::Invalid("growth requires a clamped seed with a free advancing tip"));
        }
        positive(config.maximum_deflection_m)?;
        positive(config.section.max_rotation_rad)?;
        if config.section.max_rotation_rad > 0.1 {
            return Err(Error::Regime("growth requires the same small-rotation ceiling"));
        }
        let first = &self.cells[0];
        if first.section.layers().len() != 1 {
            return Err(Error::Invalid("growth requires a born single-layer seed"));
        }
        let (&layer_id, layer) = first.section.layers().iter().next().expect("one layer");
        let geometry = layer.geometry;
        let mut ids = std::collections::BTreeSet::new();
        for cell in &self.cells {
            poll(cancelled)?;
            if cell.section.layers().len() != 1
                || cell.section.layers().get(&layer_id).is_none_or(|l| l.geometry != geometry)
                || cell.section.time_s() != self.time_s {
                return Err(Error::Invalid("growth needs one common section at the accepted time"));
            }
            ids.insert(cell.id);
        }
        if births.len() > self.budget.max_cells.saturating_sub(self.cells.len()) {
            return Err(Error::Budget);
        }
        let last = self.cells.last().expect("nonempty bridge");
        let h = last.length_m();
        let mut w = finite(last.height(h))?;
        let slope = finite(last.slope + last.section.curvature_per_m()*h)?;
        let u = finite(last.u_m + last.section.axial_strain()*h)?;
        if slope.abs() > config.section.max_rotation_rad {
            return Err(Error::Regime("advancing tip rotation exceeds the declared range"));
        }
        let previous_tip_m = self.deformed_position_m(self.length_m, 0.0)?;
        let mut next = self.clone();
        let mut mass = 0.0;
        let mut volume = 0.0;
        for birth in births {
            poll(cancelled)?;
            if !ids.insert(birth.id) { return Err(Error::Invalid("reused growth cell ID")); }
            positive(birth.length_m)?; positive(birth.mass_kg)?;
            let end = finite(next.length_m+birth.length_m)?;
            if end <= next.length_m { return Err(Error::Invalid("unrepresentable growth interval")); }
            let branches = birth.law.branch_count().checked_mul(2).ok_or(Error::Budget)?;
            next.branch_states = next.branch_states.checked_add(branches).ok_or(Error::Budget)?;
            if next.branch_states > next.budget.max_branch_states { return Err(Error::Budget); }
            let build = || -> Result<LayeredStrip<'a>, deposition::Error> {
                let mut empty = LayeredStrip::new(birth.length_m, deposition::Budget {
                    max_layers: self.budget.max_layers_per_cell,
                    max_branch_states: self.budget.max_branch_states,
                })?;
                // The existing empty-section operation advances only its clock.
                // Deposit AFTER this: no newborn material is aged from t=0.
                if self.time_s > 0.0 {
                    empty = empty.advance_to(self.time_s, &std::collections::BTreeMap::new(),
                        deposition::Control::HELD, config.section, || Ok(()))?.0;
                }
                empty.deposit(layer_id, geometry, birth.law, birth.temperature_k)
            };
            let section = build().map_err(|error| Error::Section { cell: birth.id, error })?;
            let line_mass = finite(birth.mass_kg/birth.length_m)?;
            positive(line_mass)?;
            mass = finite(mass+birth.mass_kg)?;
            volume = finite(volume+geometry.area_m2*birth.length_m)?;
            let cell = Cell { id: birth.id, start_m: next.length_m, section,
                mass_per_length_kg_m: line_mass, w_m: w, slope, u_m: u };
            let end_height = finite(cell.height(birth.length_m))?;
            if w.abs().max(end_height.abs()) > config.maximum_deflection_m {
                return Err(Error::Regime("new material leaves the declared deflection range"));
            }
            next.cells.push(cell);
            next.length_m = end;
            w = end_height;
        }
        finite(next.mass_kg())?;
        let tip_m = [finite(next.length_m+u)?, w];
        poll(cancelled)?;
        Ok((next, BirthReport { time_s: self.time_s, added_cells: births.len(),
            incoming_mass_kg: mass, incoming_reference_volume_m3: volume,
            previous_tip_m, tip_m }))
    }

    /// The identical birth operation inside the existing executor scope.
    pub fn append_cantilever(&self, cx: &fs_exec::Cx<'_>, births: &[CellBirth<'a>],
        config: Config) -> Result<(Self, BirthReport), Error> {
        self.append_cantilever_controlled(births, config, &mut || cx.checkpoint().is_err())
    }
}

#[cfg(test)]
mod tests;
