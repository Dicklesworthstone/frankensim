//! Source refresh at conserved native charge and flux. The change of the
//! dielectric Helmholtz potential is signed parameter work, not resistor heat.
//! A thermal consumer must additionally account for the source's field entropy.

use super::*;
use crate::material_circuit::MaterialResistor;
use fs_material::field_storage::admit_source_temperature_refresh;

/// Atomic parameter proposal; no electrical time or charge history is advanced.
pub struct DebyeRebindCandidate {
    circuit: DebyeCircuit,
    parameter_work_j: f64,
}
impl DebyeRebindCandidate {
    pub const fn circuit(&self) -> &DebyeCircuit {
        &self.circuit
    }
    /// Change of the original electric free energy at fixed individual charges.
    pub const fn parameter_work_j(&self) -> f64 {
        self.parameter_work_j
    }
    pub fn into_parts(self) -> (DebyeCircuit, f64) {
        (self.circuit, self.parameter_work_j)
    }
}

impl DebyeCircuit {
    /// Original graph extent, excluding the prepared internal polarization arms.
    pub fn original_branch_count(&self) -> usize {
        self.graph().branches.len() - 2 * self.source().poles().len()
    }
    /// Original capacitor branch represented by this dielectric material.
    pub fn dielectric_branch(&self) -> usize {
        self.prepared.branch
    }
    /// Explicit source records for original resistors after material binding.
    /// Polarization friction is owned by the dielectric, not by a conductor.
    pub fn material_resistors(&self) -> &[MaterialResistor] {
        &self.prepared.material_resistors
    }

    /// Refresh the complete original dielectric source at fixed geometry,
    /// individual instantaneous/polarization charges, and original fluxes.
    /// Every retained material claim and nonthermal coordinate must agree.
    /// Other resistor values and their source records remain unchanged.
    /// Incompatible ideal clamps refuse without projecting the charge history.
    pub fn rebind_source_candidate(
        &self,
        source: ResolvedDebyeCapacitor,
        inputs: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<DebyeRebindCandidate, MaterialCircuitError> {
        self.rebind(source, None, inputs, cancelled)
    }

    /// The same refresh with a complete original conductor binding for every
    /// non-polarization resistor. Native branch indices are preserved. Once a
    /// branch has a material owner, subsequent refreshes retain that owner's
    /// complete source bundle, geometry and temperature-only context.
    pub fn rebind_source_and_resistors_candidate(
        &self,
        source: ResolvedDebyeCapacitor,
        resistors: Vec<MaterialResistor>,
        inputs: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<DebyeRebindCandidate, MaterialCircuitError> {
        self.rebind(source, Some(resistors), inputs, cancelled)
    }

    fn rebind(
        &self,
        source: ResolvedDebyeCapacitor,
        resistors: Option<Vec<MaterialResistor>>,
        inputs: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DebyeRebindCandidate, MaterialCircuitError> {
        poll(&mut cancelled)?;
        let polarization = self.polarization_charges_c();
        let parameter_work_j = self.source().fixed_charge_free_energy_change_j(
            &source,
            self.instantaneous_charge_c(),
            &polarization,
        )?;
        let original_count = self.original_branch_count();
        let mut graph = CircuitGraph {
            node_count: self.graph().node_count - self.source().poles().len(),
            branches: self.graph().branches[..original_count].to_vec(),
            transformers: self.graph().transformers.clone(),
        };
        let complete_binding_requested = resistors.is_some();
        let mut records = resistors.unwrap_or_else(|| self.material_resistors().to_vec());
        records.sort_by_key(|r| r.branch);
        if complete_binding_requested || !records.is_empty() {
            let count = graph
                .branches
                .iter()
                .filter(|(_, _, b)| matches!(b, Branch::Resistor { .. }))
                .count();
            if count != records.len() || records.windows(2).any(|r| r[0].branch == r[1].branch) {
                return Err(invalid(
                    "every original resistor needs one unique complete material binding",
                ));
            }
            for record in &records {
                poll(&mut cancelled)?;
                let (_, _, branch) = graph
                    .branches
                    .get_mut(record.branch)
                    .ok_or_else(|| invalid("material resistor must name an original branch"))?;
                if !matches!(branch, Branch::Resistor { .. }) {
                    return Err(invalid(
                        "material resistor cannot replace storage or polarization",
                    ));
                }
                if let Some(before) = self
                    .material_resistors()
                    .iter()
                    .find(|r| r.branch == record.branch)
                {
                    if before.conductor.length_m() != record.conductor.length_m()
                        || before.conductor.area_m2() != record.conductor.area_m2()
                    {
                        return Err(invalid(
                            "dielectric circuit conductor geometry must remain fixed",
                        ));
                    }
                    admit_source_temperature_refresh(
                        before.conductor.material(),
                        record.conductor.material(),
                    )?;
                }
                *branch = Branch::Resistor {
                    ohms: record.conductor.resistance_ohm(),
                };
            }
        } else if !self.material_resistors().is_empty() {
            return Err(invalid(
                "dielectric circuit cannot discard original conductor source bindings",
            ));
        }
        let fluxes = self
            .dae()
            .flux_index
            .iter()
            .map(|&i| self.state[i])
            .collect();
        let original_capacitors = self.dae().charge_index.len() - polarization.len();
        let charges = self.dae().charge_index[..original_capacitors]
            .iter()
            .map(|&i| self.state[i])
            .collect();
        let mut rebound = Self::new(
            graph,
            self.dielectric_branch(),
            source,
            fluxes,
            charges,
            polarization,
            inputs,
            &mut cancelled,
        )?;
        // A new prepared owner has no aliases until this candidate is published.
        Arc::get_mut(&mut rebound.prepared)
            .ok_or_else(|| invalid("unexpected shared dielectric preparation"))?
            .material_resistors = records;
        if rebound.dae().flux_index != self.dae().flux_index
            || rebound.dae().charge_index != self.dae().charge_index
            || self
                .dae()
                .flux_index
                .iter()
                .chain(&self.dae().charge_index)
                .any(|&i| rebound.state[i].to_bits() != self.state[i].to_bits())
        {
            return Err(invalid(
                "dielectric refresh changed original conserved coordinates",
            ));
        }
        poll(&mut cancelled)?;
        Ok(DebyeRebindCandidate {
            circuit: rebound,
            parameter_work_j,
        })
    }
}
