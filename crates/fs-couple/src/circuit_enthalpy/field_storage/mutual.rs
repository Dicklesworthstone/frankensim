//! Shared-temperature thermodynamics of the native signed inductance matrix.
//!
//! L(T) has source-derived affine diagonal and optionally mutual entries.
//! At fixed total linkage, i=L^-1 lambda, psi=lambda.i/2, s=i.L'.i/2,
//! and C_field=-T (L'i).L^-1.(L'i). A connected winding group belongs to
//! ONE caloric cell; independent-temperature cross derivatives need another law.

use super::*;
use crate::material_circuit::storage::refresh_mutual;
use fs_material::field_storage::linear_magnetic_permeability;
use fs_phs::circuit::assemble_circuit_with_mutual_inductance;

pub(super) struct MutualThermal {
    groups: Vec<Vec<usize>>,
    linked: Vec<bool>,
    links: Vec<MutualInductance>,
    couplings: Vec<MaterialMagneticCoupling>,
    // One constant symmetric dL/dT block per connected group, in binding order.
    // Prepared once; endpoint evaluations perform only two Q*A products.
    derivatives: Vec<Vec<f64>>,
}

impl MutualThermal {
    pub fn new(
        circuit: &FieldCircuitView<'_>,
        bindings: &[FieldThermalBinding],
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<Self, Error> {
        let count = bindings.len();
        let branches = circuit.graph.branches.len();
        let pairs = branches.saturating_mul(branches.saturating_sub(1)) / 2;
        if circuit.mutual.len() > pairs || circuit.magnetic_couplings.len() > pairs {
            return Err(invalid(
                "mutual thermal links exceed the native branch-pair budget",
            ));
        }
        let mut result = Self {
            groups: Vec::new(),
            linked: vec![false; count],
            links: circuit.mutual.to_vec(),
            couplings: circuit.magnetic_couplings.to_vec(),
            derivatives: Vec::new(),
        };
        result.check_links(circuit)?;
        if !circuit.mutual.iter().any(|link| link.henries != 0.0) && result.couplings.is_empty() {
            return Ok(result);
        }
        crate::material_circuit::admit_graph(circuit.graph)?;
        let mut slots = BTreeMap::new();
        for (slot, binding) in bindings.iter().enumerate() {
            if slots.insert(binding.branch, slot).is_some() {
                return Err(invalid(
                    "unique field branches with finite coordinate bounds required",
                ));
            }
        }
        let mut component: Vec<_> = (0..count).collect();
        // A source-governed M can cross zero. It is still a coupled field and
        // contributes entropy there; never infer group membership from M alone.
        let active_pairs = circuit
            .mutual
            .iter()
            .filter(|link| link.henries != 0.0)
            .map(|link| link.branches)
            .chain(result.couplings.iter().map(|coupling| coupling.branches));
        for (first, second) in active_pairs {
            check(cancel)?;
            let (Some(&a), Some(&b)) = (slots.get(&first), slots.get(&second)) else {
                return Err(invalid(
                    "mutual thermal links require original field bindings",
                ));
            };
            if [a, b].iter().any(|&slot| {
                !matches!(
                    bindings[slot].law.original(),
                    ThermalFieldElement::Inductor(_)
                )
            }) {
                return Err(invalid(
                    "mutual thermal links must connect magnetic storage",
                ));
            }
            result.linked[a] = true;
            result.linked[b] = true;
            let old = component[a].max(component[b]);
            let new = component[a].min(component[b]);
            for label in &mut component {
                if *label == old {
                    *label = new;
                }
            }
        }
        let mut groups = BTreeMap::<usize, Vec<usize>>::new();
        for (slot, label) in component.into_iter().enumerate() {
            if result.linked[slot] {
                groups.entry(label).or_default().push(slot);
            }
        }
        for group in groups.into_values() {
            let owner = &bindings[group[0]];
            if group.iter().any(|&slot| {
                bindings[slot].cell != owner.cell
                    || bindings[slot].law.temperature_range_k() != owner.law.temperature_range_k()
            }) {
                return Err(invalid(
                    "coupled magnetic fields require one shared caloric cell and temperature envelope",
                ));
            }
            let width = group.len();
            let mut derivative = vec![0.0; width * width];
            for (row, &slot) in group.iter().enumerate() {
                derivative[row * width + row] = bindings[slot].law.coefficient_slope_per_k();
            }
            for coupling in &result.couplings {
                let a = slots[&coupling.branches.0];
                let b = slots[&coupling.branches.1];
                let Some(row) = group.iter().position(|&slot| slot == a) else {
                    continue;
                };
                let column = group
                    .iter()
                    .position(|&slot| slot == b)
                    .ok_or_else(|| invalid("source mutual pair lost its shared thermal group"))?;
                let slope = |slot: usize| -> Result<f64, Error> {
                    let law = &bindings[slot].law;
                    let mu =
                        linear_magnetic_permeability(law.original().material()).map_err(field)?;
                    // L_ii=mu*G_ii at fixed geometry. The source owner already
                    // admitted mu(T) as affine on the ENTIRE temperature path.
                    let ratio =
                        finite(law.coefficient_slope_per_k() / law.original().coefficient())?;
                    let value = finite(finite(coupling.geometry_m * mu)? * ratio)?;
                    if coupling.geometry_m != 0.0
                        && law.coefficient_slope_per_k() != 0.0
                        && value == 0.0
                    {
                        return Err(invalid(
                            "nonzero source mutual temperature derivative underflows",
                        ));
                    }
                    Ok(value)
                };
                let first = slope(a)?;
                let second = slope(b)?;
                if (first - second).abs() > 256.0 * f64::EPSILON * first.abs().max(second.abs()) {
                    return Err(invalid(
                        "source mutual pair requires one permeability temperature law",
                    ));
                }
                derivative[row * width + column] = first;
                derivative[column * width + row] = first;
            }
            result.derivatives.push(derivative);
            result.groups.push(group);
        }
        Ok(result)
    }

    pub fn contains(&self, slot: usize) -> bool {
        self.linked[slot]
    }

    // Retain the native scaled SPD/inverse admission rather than implementing
    // another inverse or inferring admissibility from pairwise coupling bounds.
    pub fn inverse(
        &self,
        graph: &CircuitGraph,
        links: &[MutualInductance],
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<Vec<f64>, Error> {
        if self.groups.is_empty() {
            return Ok(Vec::new());
        }
        check(cancel)?;
        let dae = assemble_circuit_with_mutual_inductance(graph, links)
            .map_err(|e| CircuitEnthalpyError::Circuit(e.into()))?;
        let stride = dae.system.state_dim();
        let matrix = dae
            .system
            .storage()
            .exact_quadratic_matrix()
            .ok_or_else(|| {
                invalid("native mutual magnetic storage requires its exact quadratic law")
            })?;
        let mut slots = vec![None; graph.branches.len()];
        let mut fluxes = dae.flux_index.iter();
        for (branch, (_, _, kind)) in graph.branches.iter().enumerate() {
            if matches!(kind, Branch::Inductor { .. }) {
                slots[branch] = fluxes.next().copied();
            }
        }
        // The group branch indices are supplied by apply/capacity below. Keep
        // only the magnetic block, in original graph inductor order.
        let indices: Vec<_> = slots.into_iter().flatten().collect();
        let mut inverse = Vec::with_capacity(indices.len() * indices.len());
        for &i in &indices {
            check(cancel)?;
            for &j in &indices {
                inverse.push(matrix[i * stride + j]);
            }
        }
        Ok(inverse)
    }

    fn magnetic_slots(
        graph: &CircuitGraph,
        bindings: &[FieldThermalBinding],
    ) -> Vec<Option<usize>> {
        let mut by_branch = BTreeMap::new();
        for (slot, (branch, _)) in graph
            .branches
            .iter()
            .enumerate()
            .filter(|(_, (_, _, kind))| matches!(kind, Branch::Inductor { .. }))
            .enumerate()
        {
            by_branch.insert(branch, slot);
        }
        bindings
            .iter()
            .map(|binding| by_branch.get(&binding.branch).copied())
            .collect()
    }

    /// Override only coupled branch reports. Each report is a signed partition;
    /// its group's sum is the physical free energy, entropy, U and capacity.
    #[allow(clippy::too_many_arguments)]
    pub fn apply(
        &self,
        graph: &CircuitGraph,
        bindings: &[FieldThermalBinding],
        inverse: &[f64],
        coordinates: &BTreeMap<usize, f64>,
        temperatures: &[f64],
        states: &mut [FieldThermalState],
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<(), Error> {
        if self.groups.is_empty() {
            return Ok(());
        }
        let slots = Self::magnetic_slots(graph, bindings);
        let n = graph
            .branches
            .iter()
            .filter(|(_, _, kind)| matches!(kind, Branch::Inductor { .. }))
            .count();
        for (group, derivative) in self.groups.iter().zip(&self.derivatives) {
            check(cancel)?;
            let temperature = temperatures[group[0]];
            if group.iter().any(|&slot| temperatures[slot] != temperature) {
                return Err(invalid(
                    "a shared magnetic group must use one actual temperature",
                ));
            }
            let mut currents = vec![0.0; group.len()];
            for (i, &slot) in group.iter().enumerate() {
                let row = slots[slot]
                    .ok_or_else(|| invalid("coupled field lost its native flux slot"))?;
                for &other in group {
                    let column = slots[other]
                        .ok_or_else(|| invalid("coupled field lost its native flux slot"))?;
                    currents[i] = finite(
                        inverse[row * n + column]
                            .mul_add(coordinates[&bindings[other].branch], currents[i]),
                    )?;
                }
            }
            // z = Q L' i = -di/dT; w = Q L' z = (d²i/dT²)/2.
            // Derivatives of EACH signed psi_i=lambda_i*i_i/2 then remain
            // consistent with its reported entropy, U and capacity partition.
            let action = |values: &[f64]| -> Result<Vec<f64>, Error> {
                let mut weighted = vec![0.0; group.len()];
                for (row, value) in weighted.iter_mut().enumerate() {
                    check(cancel)?;
                    for (column, &coordinate) in values.iter().enumerate() {
                        *value = finite(
                            derivative[row * group.len() + column].mul_add(coordinate, *value),
                        )?;
                    }
                }
                let mut result = vec![0.0; group.len()];
                for (i, &slot) in group.iter().enumerate() {
                    check(cancel)?;
                    let row = slots[slot].unwrap();
                    for (j, &other) in group.iter().enumerate() {
                        result[i] = finite(
                            inverse[row * n + slots[other].unwrap()]
                                .mul_add(weighted[j], result[i]),
                        )?;
                    }
                }
                Ok(result)
            };
            let z = action(&currents)?;
            let w = action(&z)?;
            for (i, &slot) in group.iter().enumerate() {
                check(cancel)?;
                let linkage = coordinates[&bindings[slot].branch];
                let free_energy_j = finite(0.5 * linkage * currents[i])?;
                let entropy_j_k = finite(0.5 * linkage * z[i])?;
                if linkage != 0.0 && z[i] != 0.0 && entropy_j_k == 0.0 {
                    return Err(invalid("nonzero shared magnetic entropy underflows"));
                }
                states[slot] = FieldThermalState {
                    free_energy_j,
                    entropy_j_k,
                    internal_energy_j: finite(temperature.mul_add(entropy_j_k, free_energy_j))?,
                    heat_capacity_j_k: finite((-temperature * linkage) * w[i])?,
                };
            }
        }
        Ok(())
    }

    pub fn check_links(&self, circuit: &FieldCircuitView<'_>) -> Result<(), Error> {
        if circuit.magnetic_couplings != self.couplings {
            return Err(invalid("thermal trial changed its source magnetic overlap"));
        }
        let expected = refresh_mutual(&self.links, &self.couplings, circuit.storage)?;
        if circuit.mutual != expected {
            return Err(invalid(
                "thermal trial mutual inductance disagrees with its source model",
            ));
        }
        Ok(())
    }

    /// L(T) is affine. Positive endpoint matrices admit its whole interval, and
    /// inverse operator convexity bounds ||L(T)^-1|| by the two endpoint bounds.
    /// C >= -Tmax ||L'||^2 ||lambda_bound||^2 ||L^-1||^3. This is conservative
    /// analytic model admission with numerical padding, not an interval proof.
    pub fn capacity_floors(
        &self,
        graph: &CircuitGraph,
        bindings: &[FieldThermalBinding],
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<Vec<(CellId, f64)>, Error> {
        if self.groups.is_empty() {
            return Ok(Vec::new());
        }
        let slots = Self::magnetic_slots(graph, bindings);
        let n = graph
            .branches
            .iter()
            .filter(|(_, _, kind)| matches!(kind, Branch::Inductor { .. }))
            .count();
        let mut norms = vec![0.0_f64; self.groups.len()];
        for endpoint in 0..2 {
            let mut at = graph.clone();
            let mut storage = Vec::new();
            for group in &self.groups {
                for &slot in group {
                    check(cancel)?;
                    let binding = &bindings[slot];
                    let sample = binding
                        .law
                        .at_temperature(binding.law.temperature_range_k()[endpoint])
                        .map_err(field)?;
                    at.branches[binding.branch].2 = Branch::Inductor {
                        henries: sample.element().coefficient(),
                    };
                    storage.push(MaterialStorage {
                        branch: binding.branch,
                        element: element(sample.element()),
                    });
                }
            }
            let links = refresh_mutual(&self.links, &self.couplings, &storage)?;
            let inverse = self.inverse(&at, &links, cancel)?;
            for (index, group) in self.groups.iter().enumerate() {
                let mut trace = 0.0;
                let mut row_norm = 0.0_f64;
                for &slot in group {
                    let row = slots[slot]
                        .ok_or_else(|| invalid("coupled capacity lost its native flux slot"))?;
                    trace = finite(trace + inverse[row * n + row])?;
                    let mut sum = 0.0;
                    for &other in group {
                        sum = finite(sum + inverse[row * n + slots[other].unwrap()].abs())?;
                    }
                    row_norm = row_norm.max(sum);
                }
                // Native dimensionless inverse checks allow 1e-10 per entry;
                // include that residual scale as well as reduction roundoff.
                let bound = finite(
                    trace.min(row_norm) * (1.0 + 4096.0 * f64::EPSILON * n as f64)
                        / (1.0 - 1e-10 * n as f64),
                )?;
                norms[index] = norms[index].max(bound);
            }
        }
        let mut result = Vec::with_capacity(self.groups.len());
        for ((group, derivative), norm) in self.groups.iter().zip(&self.derivatives).zip(norms) {
            let mut linkage_norm = 0.0_f64;
            let mut slope = 0.0_f64;
            for (row, &slot) in group.iter().enumerate() {
                linkage_norm = linkage_norm.hypot(bindings[slot].maximum_abs_coordinate);
                let mut row_norm = 0.0;
                for &entry in &derivative[row * group.len()..(row + 1) * group.len()] {
                    row_norm = finite(row_norm + entry.abs())?;
                }
                // For symmetric A, max absolute row sum bounds ||A||_2,
                // including every temperature-dependent mutual cross term.
                slope = slope.max(row_norm);
            }
            let bound = if slope == 0.0 || linkage_norm == 0.0 {
                0.0
            } else {
                // Form squared ||L'||||Q||||lambda|| before the final Q factor
                // to avoid cubing inverse inductance in isolation.
                let response = finite(finite(slope * norm)? * linkage_norm)?;
                let magnitude = finite(
                    finite(response * response)?
                        * norm
                        * bindings[group[0]].law.temperature_range_k()[1],
                )?;
                if magnitude == 0.0 {
                    return Err(invalid("nonzero shared magnetic capacity bound underflows"));
                }
                -magnitude
            };
            result.push((bindings[group[0]].cell, bound));
        }
        Ok(result)
    }
}
