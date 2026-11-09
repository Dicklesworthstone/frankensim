//! State-dependent conductivity on the existing fixed pyrolysis heat graph.
//!
//! Original constituent k(h) sources are evaluated at the same common
//! temperature as reaction and storage, using each ORIGINAL constituent h(T).
//! Explicit Series/Parallel arrangements use actual constituent volumes.
//! Source densities describe those constituent volumes only: disappearing
//! volatile mass does not invent a pore-gas conductivity, porosity law or
//! shrinking geometry. A porous-char source must itself admit that effective
//! material interpretation. All bound face metrics remain fixed.

use super::slab::{ConductivityMixing, mixed_conductivity};
use super::*;
use crate::enthalpy::conduction::{ConductiveLink, ConductivityLaw};
use crate::enthalpy::{ThermalLink, series_contact_conductance_w_k};

/// Explicit source selection on the original material's stored-energy chart.
#[derive(Clone, Copy)]
pub enum PyrolysisConductivity<'a> {
    /// Original k(e) of a heat-only body, including variation along latent h.
    Thermal(&'a dyn ConductivityLaw),
    /// Original virgin and char sources and an explicitly chosen arrangement.
    /// Each source receives its own h(T), never the bulk mixture enthalpy.
    Reactive {
        virgin: &'a dyn ConductivityLaw,
        /// Required when the original material retains char.
        char: Option<&'a dyn ConductivityLaw>,
        /// Constant source-owned [virgin, char] densities [kg/m³]. These only
        /// determine mixture volume fractions; reference geometry stays fixed.
        /// Char density may be zero only when the material produces no char.
        densities_kg_m3: [f64; 2],
        mixing: ConductivityMixing,
        /// Applicability of the constituent densities and arrangement.
        source: &'a str,
    },
}

#[derive(Clone)]
pub(super) struct Conduction<'a> {
    laws: BTreeMap<MaterialId, PyrolysisConductivity<'a>>,
    pub(super) links: BTreeMap<u64, ConductiveLink>,
    /// One current source query per participating cell, shared by its faces.
    conductivity: BTreeMap<CellId, f64>,
}

#[derive(Clone, Copy, Default)]
pub(super) struct ConductivityProbe {
    pub conductivity_w_m_k: f64,
    pub energy_difference_j: f64,
}

impl<'a> PyrolysisSimulation<'a> {
    /// Attach a bounded subset of fixed normal conduction/contact geometries.
    ///
    /// Existing `ConductiveLink` metrics and original series/contact arithmetic
    /// determine G. Unbound links and authored boundary coefficients retain
    /// their existing owner. Bound links cannot also be manually overwritten.
    /// The material sources must be immutable/repeatable across trials.
    /// Receding-slab adoption refuses this fixed-geometry attachment.
    pub fn with_conductivity(
        mut self,
        laws: &[(MaterialId, PyrolysisConductivity<'a>)],
        links: &[ConductiveLink],
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DepositionError> {
        checkpoint()?;
        if self.conduction.is_some()
            || links.is_empty()
            || links.len() > self.thermal.links.len()
            || laws.len() > self.thermal.materials.len()
        {
            return Err(DepositionError::InvalidInput(
                "bounded unique pyrolysis conductivity attachment required",
            ));
        }
        let mut binding = Conduction {
            laws: BTreeMap::new(),
            links: BTreeMap::new(),
            conductivity: BTreeMap::new(),
        };
        for &(id, law) in laws {
            checkpoint()?;
            if !self.thermal.materials.contains_key(&id) || binding.laws.insert(id, law).is_some() {
                return Err(DepositionError::InvalidInput(
                    "unknown or repeated pyrolysis conductivity material",
                ));
            }
            match (self.materials.get(&id), law) {
                (None, PyrolysisConductivity::Thermal(_)) => {}
                (
                    Some(model),
                    PyrolysisConductivity::Reactive {
                        char,
                        densities_kg_m3,
                        source,
                        ..
                    },
                ) if positive(densities_kg_m3[0])
                    && nonnegative(densities_kg_m3[1])
                    && !source.trim().is_empty()
                    && (model.char_mass_yield == 0.0
                        || (char.is_some() && positive(densities_kg_m3[1]))) => {}
                _ => {
                    return Err(DepositionError::InvalidInput(
                        "conductivity must match original reactive or heat-only material",
                    ));
                }
            }
        }
        for &geometry in links {
            checkpoint()?;
            let link =
                self.thermal
                    .links
                    .get(&geometry.link_id)
                    .ok_or(DepositionError::InvalidInput(
                        "inactive pyrolysis conductivity link",
                    ))?;
            if binding.links.insert(link.id, geometry).is_some() {
                return Err(DepositionError::InvalidInput(
                    "repeated pyrolysis conductivity link",
                ));
            }
            // Original series owner validates distances/contact even when k=0.
            conductance(link, geometry, [0.0; 2])?;
            for id in [link.a, link.b] {
                if !binding.laws.contains_key(&self.cells()[&id].material_id) {
                    return Err(DepositionError::InvalidInput(
                        "missing linked pyrolysis conductivity source",
                    ));
                }
                binding.conductivity.insert(id, 0.0);
            }
        }
        self.conduction = Some(binding);
        self.refresh_conductivity(&mut checkpoint)?;
        checkpoint()?;
        Ok(self)
    }

    /// Current resolved scalar k [W/(m K)] on an explicitly attached cell.
    /// None leaves conductivity with the original caller-supplied link owner.
    pub fn conductivity_w_m_k(&self, id: CellId) -> Option<f64> {
        self.conduction.as_ref()?.conductivity.get(&id).copied()
    }

    pub(super) fn refresh_conductivity(
        &mut self,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DepositionError> {
        let Some(binding) = &self.conduction else {
            return Ok(());
        };
        let mut conductivity = BTreeMap::new();
        for &id in binding.conductivity.keys() {
            checkpoint()?;
            conductivity.insert(id, binding.cell_conductivity(self, id, None)?);
        }
        let mut values = BTreeMap::new();
        for (&id, &geometry) in &binding.links {
            checkpoint()?;
            let link = &self.thermal.links[&id];
            values.insert(
                id,
                conductance(
                    link,
                    geometry,
                    [conductivity[&link.a], conductivity[&link.b]],
                )?,
            );
        }
        // Retain the original coefficient admission and paired heat owner.
        self.thermal.set_step_conductances(&values, checkpoint)?;
        self.conduction.as_mut().expect("bound above").conductivity = conductivity;
        Ok(())
    }
}

impl Conduction<'_> {
    pub(super) fn cell_conductivity(
        &self,
        state: &PyrolysisSimulation<'_>,
        id: CellId,
        probe: Option<&CellStep>,
    ) -> Result<f64, DepositionError> {
        let cell = &state.cells()[&id];
        let value = match self.laws[&cell.material_id] {
            PyrolysisConductivity::Thermal(law) => {
                law.conductivity_w_m_k(probe.map_or(cell.specific_h_j_kg, |p| p.specific_h_j_kg))?
            }
            PyrolysisConductivity::Reactive {
                virgin,
                char,
                densities_kg_m3,
                mixing,
                ..
            } => {
                let inventory = probe.map_or(state.inventories[&id], |p| p.inventory);
                let model = state.materials[&cell.material_id];
                let temperature = inventory.temperature_k;
                let kv = if inventory.virgin_mass_kg == 0.0 {
                    0.0
                } else {
                    virgin.conductivity_w_m_k(original_h(model.virgin, temperature)?)?
                };
                let kc = if inventory.char_mass_kg == 0.0 {
                    0.0
                } else {
                    char.ok_or(DepositionError::MaterialDomain)?
                        .conductivity_w_m_k(original_h(
                            model.char_caloric.ok_or(DepositionError::MaterialDomain)?,
                            temperature,
                        )?)?
                };
                mixed_conductivity(
                    [
                        inventory.virgin_mass_kg / densities_kg_m3[0],
                        if inventory.char_mass_kg == 0.0 {
                            0.0
                        } else {
                            inventory.char_mass_kg / densities_kg_m3[1]
                        },
                    ],
                    [kv, kc],
                    mixing,
                )?
            }
        };
        if !nonnegative(value) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(value)
    }

    pub(super) fn probe(
        &self,
        state: &PyrolysisSimulation<'_>,
        id: CellId,
        step: &CellStep,
        delta: f64,
    ) -> Result<ConductivityProbe, DepositionError> {
        Ok(ConductivityProbe {
            conductivity_w_m_k: if self.conductivity.contains_key(&id) {
                self.cell_conductivity(state, id, Some(step))?
            } else {
                0.0
            },
            energy_difference_j: delta,
        })
    }

    /// Complete dG/d(absorbed energy) contribution to the same paired flux.
    /// Original one-cell reaction probes supply T, h AND component inventory;
    /// no full graph/reactive solve is repeated for each matrix column.
    pub(super) fn add_tangent(
        &self,
        state: &PyrolysisSimulation<'_>,
        ids: &[CellId],
        steps: &[CellStep],
        probes: &[ConductivityProbe],
        dt: f64,
        matrix: &mut [f64],
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DepositionError> {
        let n = ids.len();
        for (&id, &geometry) in &self.links {
            checkpoint()?;
            let link = &state.thermal.links[&id];
            let a = ids
                .binary_search(&link.a)
                .map_err(|_| DepositionError::Arithmetic)?;
            let b = ids
                .binary_search(&link.b)
                .map_err(|_| DepositionError::Arithmetic)?;
            let temperature_difference =
                finite(steps[a].inventory.temperature_k - steps[b].inventory.temperature_k)?;
            let k = [self.conductivity[&link.a], self.conductivity[&link.b]];
            for (side, column) in [a, b].into_iter().enumerate() {
                let mut trial = k;
                trial[side] = probes[column].conductivity_w_m_k;
                let g = conductance(link, geometry, trial)?;
                let derivative =
                    finite((g - link.conductance_w_k) / probes[column].energy_difference_j)?;
                let change = finite(dt * finite(derivative * temperature_difference)?)?;
                matrix[a * n + column] = finite(matrix[a * n + column] + change)?;
                matrix[b * n + column] = finite(matrix[b * n + column] - change)?;
            }
        }
        Ok(())
    }
}

fn conductance(
    link: &ThermalLink,
    geometry: ConductiveLink,
    k: [f64; 2],
) -> Result<f64, DepositionError> {
    series_contact_conductance_w_k(
        link.area_m2,
        geometry.half_distances_m,
        k,
        geometry.contact_resistance_m2_k_w,
    )
}
