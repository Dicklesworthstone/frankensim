//! Material-derived field storage on the existing Kirchhoff circuit runtime.
//!
//! Charge and flux linkage survive material/geometry updates. The change in
//! their Hamiltonian is explicit parameter work, never unearned Joule heat or
//! a reset of electrical state. The caller must couple that work to the body,
//! heat reservoir or external control that actually changed the coefficient.

use fs_material::field_storage::{ResolvedCapacitor, ResolvedInductor};
use fs_phs::circuit::{Branch, CircuitGraph, MutualInductance, assemble_circuit_with_mutual_inductance};

use super::{MaterialCircuit, MaterialCircuitError, MaterialResistor, poll};

/// The original field-storage source and geometry for one graph branch.
#[derive(Clone, Debug, PartialEq)]
pub enum MaterialStorageElement {
    /// Uniform dielectric-gap reduction, C = epsilon A / d.
    Capacitor(ResolvedCapacitor),
    /// Uniform magnetic-path reduction, L = mu N^2 A / l.
    Inductor(ResolvedInductor),
}

/// Bind a material storage element to an ORIGINAL circuit branch index.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialStorage {
    /// Full graph index, not the capacitor-only or inductor-only position.
    pub branch: usize,
    /// Source material and geometry; a numeric graph coefficient is replaced.
    pub element: MaterialStorageElement,
}

/// One branch's contribution to parameter work at fixed conserved coordinates.
/// With mutual coupling, magnetic contributions are `lambda_i * i_i / 2`.
/// An individual contribution may be negative; only their sum is the physical
/// magnetic storage. Mutual energy must never be added a second time.
#[derive(Clone, Debug, PartialEq)]
pub struct StorageParameterWork {
    /// Original graph branch.
    pub branch: usize,
    /// Electric energy or signed magnetic energy contribution before [J].
    pub stored_before_j: f64,
    /// Same charge/total flux evaluated with the complete new storage law [J].
    pub stored_after_j: f64,
    /// Positive means supplied TO the electrical system; negative means returned.
    pub supplied_work_j: f64,
}

/// An uncommitted source/material update with explicitly required parameter work.
/// A proposal never changes the original circuit; the caller publishes this
/// together with the physical owner that pays/receives the reported work.
pub struct MaterialStorageCandidate {
    circuit: MaterialCircuit,
    work: Vec<StorageParameterWork>,
    parameter_work_j: f64,
}
impl MaterialStorageCandidate {
    /// New material coefficients and unchanged conserved charge/flux coordinates.
    pub const fn circuit(&self) -> &MaterialCircuit {
        &self.circuit
    }
    /// Ordered per-branch field-energy changes.
    pub fn work(&self) -> &[StorageParameterWork] {
        &self.work
    }
    /// External parameter work required by this update, with its physical sign.
    pub const fn parameter_work_j(&self) -> f64 {
        self.parameter_work_j
    }
    /// Move the candidate and its work account into a joint physical transaction.
    pub fn into_parts(self) -> (MaterialCircuit, Vec<StorageParameterWork>, f64) {
        (self.circuit, self.work, self.parameter_work_j)
    }
}

fn invalid(message: &'static str) -> MaterialCircuitError {
    MaterialCircuitError::Invalid(message)
}

pub(crate) fn lower(
    graph: &CircuitGraph,
    sources: &mut [MaterialStorage],
) -> Result<CircuitGraph, MaterialCircuitError> {
    if graph.node_count > super::MAX_MATERIAL_CIRCUIT_NODES
        || graph.branches.len() > super::MAX_MATERIAL_CIRCUIT_BRANCHES
    {
        return Err(invalid(
            "material storage graph exceeds the circuit capacity",
        ));
    }
    let count = graph
        .branches
        .iter()
        .filter(|(_, _, branch)| {
            matches!(branch, Branch::Capacitor { .. } | Branch::Inductor { .. })
        })
        .count();
    if sources.len() != count {
        return Err(invalid(
            "every capacitor and inductor requires exactly one material storage binding",
        ));
    }
    sources.sort_by_key(|source| source.branch);
    if sources
        .windows(2)
        .any(|pair| pair[0].branch == pair[1].branch)
    {
        return Err(invalid("duplicate material storage branch"));
    }
    let mut result = graph.clone();
    for source in sources {
        let branch = result
            .branches
            .get_mut(source.branch)
            .ok_or_else(|| invalid("unknown material storage branch"))?;
        branch.2 = match (&source.element, branch.2) {
            (MaterialStorageElement::Capacitor(capacitor), Branch::Capacitor { .. }) => {
                Branch::Capacitor {
                    farads: capacitor.capacitance_f(),
                }
            }
            (MaterialStorageElement::Inductor(inductor), Branch::Inductor { .. }) => {
                Branch::Inductor {
                    henries: inductor.inductance_h(),
                }
            }
            _ => {
                return Err(invalid(
                    "material storage type does not match its circuit branch",
                ));
            }
        };
    }
    Ok(result)
}

// Shared fixed-coordinate account for the native lumped and spatial-conductor
// circuit paths. Both use the original material element's energy evaluator.
pub(crate) fn parameter_work(
    graph: &CircuitGraph,
    original: &[MaterialStorage],
    proposed: &[MaterialStorage],
    fluxes_wb: &[f64],
    charges_c: &[f64],
    mut cancelled: impl FnMut() -> bool,
) -> Result<(Vec<StorageParameterWork>, f64), MaterialCircuitError> {
    poll(&mut cancelled)?;
    if original.len() != proposed.len() {
        return Err(invalid(
            "complete original and proposed field storage required",
        ));
    }
    let mut charges = charges_c.iter();
    let mut fluxes = fluxes_wb.iter();
    let mut source = 0;
    let mut work = Vec::with_capacity(original.len());
    let mut total = 0.0;
    for (index, &(_, _, branch)) in graph.branches.iter().enumerate() {
        poll(&mut cancelled)?;
        let coordinate = match branch {
            Branch::Capacitor { .. } => *charges
                .next()
                .ok_or_else(|| invalid("missing capacitor charge"))?,
            Branch::Inductor { .. } => *fluxes
                .next()
                .ok_or_else(|| invalid("missing inductor flux linkage"))?,
            _ => continue,
        };
        let before = original
            .get(source)
            .ok_or_else(|| invalid("missing original field storage"))?;
        let after = proposed
            .get(source)
            .ok_or_else(|| invalid("missing proposed field storage"))?;
        if before.branch != index || after.branch != index {
            return Err(invalid("field storage must retain original branch order"));
        }
        let evaluate = |element: &MaterialStorageElement| match element {
            MaterialStorageElement::Capacitor(c) => c.stored_energy_j(coordinate),
            MaterialStorageElement::Inductor(l) => l.stored_energy_j(coordinate),
        };
        let before = evaluate(&before.element)?;
        let after = evaluate(&after.element)?;
        let supplied = after - before;
        total += supplied;
        if !supplied.is_finite() || !total.is_finite() {
            return Err(invalid(
                "material storage parameter work is unrepresentable",
            ));
        }
        work.push(StorageParameterWork {
            branch: index,
            stored_before_j: before,
            stored_after_j: after,
            supplied_work_j: supplied,
        });
        source += 1;
    }
    if source != original.len() || charges.next().is_some() || fluxes.next().is_some() {
        return Err(invalid("field coordinates and source bindings disagree"));
    }
    poll(&mut cancelled)?;
    Ok((work, total))
}

impl MaterialCircuit {
    /// Construct the ordinary circuit with ALL R, C and L values supplied by
    /// their material owners. Numeric graph values are placeholders only.
    /// No alternate circuit solver or stored-state representation is introduced.
    ///
    /// Initial source/constraint consistency is checked by the existing DAE.
    /// Uniform-field approximations and source validity belong to the material
    /// reductions. C/L remain fixed during each ordinary `step_candidate`.
    pub fn new_with_storage(
        graph: CircuitGraph,
        resistors: Vec<MaterialResistor>,
        storage: Vec<MaterialStorage>,
        fluxes_wb: Vec<f64>,
        charges_c: Vec<f64>,
        initial_inputs: &[f64],
    ) -> Result<Self, MaterialCircuitError> {
        Self::new_with_storage_and_mutual_inductance(graph, resistors, storage,
            Vec::new(), fluxes_wb, charges_c, initial_inputs)
    }

    /// Source-derived diagonal C/L and explicit signed mutual coupling on the
    /// same native descriptor. Each inductance source owns its self L; authored
    /// M supplies the shared-field cross terms without another energy store.
    /// The full matrix must remain positive definite. No common-core geometry
    /// or coupling coefficient is inferred from the individual material cards.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_storage_and_mutual_inductance(
        graph: CircuitGraph,
        resistors: Vec<MaterialResistor>,
        mut storage: Vec<MaterialStorage>,
        mutual: Vec<MutualInductance>,
        fluxes_wb: Vec<f64>,
        charges_c: Vec<f64>,
        initial_inputs: &[f64],
    ) -> Result<Self, MaterialCircuitError> {
        let graph = lower(&graph, &mut storage)?;
        let mut circuit = Self::new_with_mutual_inductance(graph, resistors, mutual,
            fluxes_wb, charges_c, initial_inputs)?;
        circuit.storage = storage;
        Ok(circuit)
    }

    /// Source-bound storage, retained unchanged through ordinary circuit steps.
    /// Empty on the original constructor's explicitly authored C/L path.
    pub fn storage(&self) -> &[MaterialStorage] {
        &self.storage
    }

    /// Rebind source-derived C/L at fixed charge and flux linkage, without time
    /// advance. This is an impulse-free parameter change with explicit work:
    /// W = sum(q^2/2 * delta(1/C)) + lambda^T delta(L^-1) lambda / 2.
    /// Mutual coefficients remain fixed, while native admission rechecks the
    /// entire new L. Branch work includes the resulting shared-field terms.
    ///
    /// The existing descriptor must admit the SAME coordinates under the new
    /// algebraic constraints. A capacitor clamped by an ideal voltage source,
    /// for example, cannot generally change C with both q and voltage fixed;
    /// that update refuses instead of inventing an unaccounted current impulse.
    ///
    /// This does not solve mechanical force, heat, hysteresis or continuous
    /// parameter evolution. Positive work must be supplied by the changing
    /// material/geometry's physical owner, and negative work returned to it.
    /// Any source, budget, constraint or cancellation refusal leaves `self`
    /// unchanged. Every branch remains on the same circuit topology.
    pub fn rebind_storage_candidate(
        &self,
        storage: Vec<MaterialStorage>,
        inputs: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialStorageCandidate, MaterialCircuitError> {
        poll(&mut cancelled)?;
        let storage_count = self
            .graph
            .branches
            .iter()
            .filter(|(_, _, branch)| {
                matches!(branch, Branch::Capacitor { .. } | Branch::Inductor { .. })
            })
            .count();
        if self.storage.len() != storage_count {
            return Err(invalid(
                "live storage changes require source-bound original capacitors and inductors",
            ));
        }
        // Storage changes must retain the accepted contact drops and sources.
        let mut storage = storage;
        let graph = lower(&self.graph, &mut storage)?;
        let mut candidate = Self::new_with_contacts_and_mutual_inductance(
            graph,
            self.resistors.clone(),
            self.contacts.clone(),
            self.mutual.clone(),
            self.fluxes_wb.clone(),
            self.charges_c.clone(),
            inputs,
        )?;
        candidate.storage = storage;
        let (work, total) = if self.mutual.is_empty() { parameter_work(
            &self.graph,
            &self.storage,
            &candidate.storage,
            &self.fluxes_wb,
            &self.charges_c,
            &mut cancelled,
        )? } else { coupled_parameter_work(self, &candidate, inputs, &mut cancelled)? };
        Ok(MaterialStorageCandidate {
            circuit: candidate,
            work,
            parameter_work_j: total,
        })
    }
}

// Partition the original owner's full quadratic energy, not lambda_i^2/2L_ii.
// The latter drops all induced-current changes when any shared self L changes.
fn coupled_parameter_work(
    before: &MaterialCircuit,
    after: &MaterialCircuit,
    inputs: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<StorageParameterWork>, f64), MaterialCircuitError> {
    let partition = |circuit: &MaterialCircuit| -> Result<Vec<f64>, MaterialCircuitError> {
        let dae = assemble_circuit_with_mutual_inductance(&circuit.graph, &circuit.mutual)?;
        let state = super::extended(&dae, &circuit.fluxes_wb, &circuit.charges_c, inputs)?;
        let mut efforts = vec![0.0; state.len()];
        dae.system.storage().gradient(&state, &mut efforts);
        let mut fluxes = dae.flux_index.iter();
        let mut charges = dae.charge_index.iter();
        let mut values = Vec::with_capacity(circuit.storage.len());
        for &(_, _, branch) in &circuit.graph.branches {
            let index = match branch {
                Branch::Inductor { .. } => fluxes.next(),
                Branch::Capacitor { .. } => charges.next(),
                _ => continue,
            }.ok_or_else(|| invalid("missing coupled field coordinate"))?;
            let energy = 0.5 * state[*index] * efforts[*index];
            if !energy.is_finite() {
                return Err(invalid("coupled field energy contribution is unrepresentable"));
            }
            values.push(energy);
        }
        Ok(values)
    };
    poll(cancelled)?;
    let original = partition(before)?;
    poll(cancelled)?;
    let proposed = partition(after)?;
    let mut work = Vec::with_capacity(original.len());
    let mut total = 0.0;
    for ((source, stored_before_j), stored_after_j) in before.storage.iter().zip(original).zip(proposed) {
        poll(cancelled)?;
        let supplied_work_j = stored_after_j - stored_before_j;
        total += supplied_work_j;
        if !supplied_work_j.is_finite() || !total.is_finite() {
            return Err(invalid("coupled storage parameter work is unrepresentable"));
        }
        work.push(StorageParameterWork { branch: source.branch, stored_before_j,
            stored_after_j, supplied_work_j });
    }
    poll(cancelled)?;
    Ok((work, total))
}
