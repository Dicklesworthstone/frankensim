use super::super::{
    DiffusionMaterial, MAX_OSMOTIC_SPECIES, OsmoticCompartment, OsmoticNetworkError,
    SolutionDiffusionControls, SolutionDiffusionLink, SolutionDiffusionStepRefusal,
    SolutionDiffusionWorkspace, compensated_add, diffusion_source, finite, implicit, nonnegative,
    poll, positive, source,
};
use super::{ElectroneutralDiffusion, ElectroneutralDiffusionControls};
use fs_material::osmosis::diffusion::electrical::IonicDiffusionResponse;
use fs_solver::{
    LinearOp,
    checked_cg::{CheckedCgConfig, CheckedCgError, checked_cg},
};
use fs_sparse::precond::Precond;

struct Face<'a> {
    endpoints: [usize; 2],
    area_m2: f64,
    weight: f64,
    response: IonicDiffusionResponse<'a>,
    workspace: &'a mut SolutionDiffusionWorkspace,
}

pub(super) fn workspace_bytes(
    nodes: usize,
    links: &[SolutionDiffusionLink],
    controls: ElectroneutralDiffusionControls,
) -> Result<usize, OsmoticNetworkError> {
    let faces = links.iter().filter(|link| link.area_m2 != 0.0).count();
    let mut bytes = faces.checked_mul(core::mem::size_of::<Face<'_>>());
    for link in links
        .iter()
        .filter(|link| link.area_m2 != 0.0 && link.material_sources().is_some())
    {
        let n = link.material().species().len() + 1;
        bytes = bytes.and_then(|b| b.checked_add(n.checked_mul(n)?.checked_mul(8)?));
    }
    // checked_cg owns original/correction iterates, defect vectors, native CG
    // state and deterministic-dot scratch. This deliberately generous payload
    // estimate also covers its temporary vector overlaps and bounded histories.
    bytes
        .and_then(|b| b.checked_add(nodes.checked_mul(40 * core::mem::size_of::<f64>())?))
        .and_then(|b| b.checked_add((controls.maximum_iterations_per_solve + 1).checked_mul(32)?))
        .ok_or(OsmoticNetworkError::Invalid(
            "electroneutral transient storage overflow",
        ))
}

impl ElectroneutralDiffusion {
    fn face_evaluation(&mut self, evaluations: &mut usize) -> Result<(), OsmoticNetworkError> {
        if self.report.face_evaluations >= self.controls.maximum_face_evaluations {
            return Err(OsmoticNetworkError::WorkBudget);
        }
        self.report.face_evaluations += 1;
        *evaluations = evaluations
            .checked_add(1)
            .ok_or(OsmoticNetworkError::WorkBudget)?;
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn sweep(
        &mut self,
        cells: &mut [OsmoticCompartment],
        links: &mut [SolutionDiffusionLink],
        dt: f64,
        controls: SolutionDiffusionControls,
        evaluations: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), OsmoticNetworkError> {
        poll(cancelled)?;
        let n = cells[0].solution.solute_moles().len();
        let mut faces = Vec::new();
        faces
            .try_reserve_exact(links.iter().filter(|link| link.area_m2 != 0.0).count())
            .map_err(|_| OsmoticNetworkError::Invalid("ionic face response allocation failed"))?;
        for link in links.iter_mut().filter(|link| link.area_m2 != 0.0) {
            poll(cancelled)?;
            self.face_evaluation(evaluations)?;
            let ionic = link
                .electrical
                .as_deref()
                .ok_or(OsmoticNetworkError::Invalid(
                    "global potential requires the original ionic assignment on every face",
                ))?;
            let (left, right) = (&cells[link.left].solution, &cells[link.right].solution);
            let mut chemical_flux = [0.0; MAX_OSMOTIC_SPECIES];
            // The response borrows only the immutable original owner/policy.
            // Its matching mutable workspace remains available for final replay.
            let response = match &link.material {
                DiffusionMaterial::Constant(material) => material.current_response_into(
                    left,
                    right,
                    link.distance_m,
                    controls.constitutive,
                    ionic,
                    &mut link.workspace,
                    &mut chemical_flux[..n],
                    &mut *cancelled,
                ),
                DiffusionMaterial::Source(material) => material.current_response_into(
                    left,
                    right,
                    link.distance_m,
                    controls.constitutive,
                    ionic,
                    &mut link.workspace,
                    &mut chemical_flux[..n],
                    &mut *cancelled,
                ),
            }
            .map_err(diffusion_source)?;
            let conductance = response.electric_conductance_mol_m2_s();
            let weight = finite(link.area_m2 * conductance)?;
            if conductance > 0.0 && !positive(weight) {
                return Err(OsmoticNetworkError::Invalid(
                    "ionic edge conductance is unrepresentable",
                ));
            }
            faces.push(Face {
                endpoints: [link.left, link.right],
                area_m2: link.area_m2,
                weight,
                response,
                workspace: &mut link.workspace,
            });
        }
        self.solve_potential(&mut faces, cancelled)?;
        self.amounts.fill(0.0);
        self.amount_correction.fill(0.0);
        self.withdrawn.fill(0.0);
        self.withdrawal_correction.fill(0.0);
        self.charge_rate.fill(0.0);
        self.charge_correction.fill(0.0);
        let mut maximum_edge_rate: f64 = 0.0;
        let mut maximum_potential: f64 = 0.0;
        for face in &mut faces {
            poll(cancelled)?;
            self.face_evaluation(evaluations)?;
            let [left, right] = face.endpoints;
            let potential = finite(self.potentials[right] - self.potentials[left])?;
            let mut flux = [0.0; MAX_OSMOTIC_SPECIES];
            let report = face
                .response
                .flux_at_potential_into(potential, face.workspace, &mut flux[..n], &mut *cancelled)
                .map_err(diffusion_source)?;
            let current = finite(face.area_m2 * report.charge_flux_mol_m2_s)?;
            if report.charge_flux_mol_m2_s != 0.0 && current == 0.0 {
                return Err(OsmoticNetworkError::Invalid(
                    "ionic edge current is unrepresentable",
                ));
            }
            compensated_add(
                &mut self.charge_rate[left],
                &mut self.charge_correction[left],
                current,
            );
            compensated_add(
                &mut self.charge_rate[right],
                &mut self.charge_correction[right],
                -current,
            );
            maximum_edge_rate = maximum_edge_rate.max(current.abs());
            maximum_potential = maximum_potential.max(potential.abs());
            for (i, &rate) in flux[..n].iter().enumerate() {
                let transfer = finite(implicit::integrated_flux(dt, face.area_m2, rate))?;
                if rate != 0.0 && transfer == 0.0 {
                    return Err(OsmoticNetworkError::StepUnderflow);
                }
                let (a, b) = (left * n + i, right * n + i);
                compensated_add(
                    &mut self.amounts[a],
                    &mut self.amount_correction[a],
                    -transfer,
                );
                compensated_add(
                    &mut self.amounts[b],
                    &mut self.amount_correction[b],
                    transfer,
                );
                let donor = if transfer >= 0.0 { a } else { b };
                compensated_add(
                    &mut self.withdrawn[donor],
                    &mut self.withdrawal_correction[donor],
                    transfer.abs(),
                );
            }
        }
        // Free all law/workspace borrows before qualifying the proposed state.
        // Sourced snapshots cannot accidentally survive into another stage.
        drop(faces);
        let mut maximum_node_rate: f64 = 0.0;
        for (compartment, &rate) in self.charge_rate.iter().enumerate() {
            if !rate.is_finite() || rate.abs() > self.controls.maximum_node_charge_rate_mol_s {
                return Err(OsmoticNetworkError::ElectroneutralChargeRate {
                    compartment,
                    residual_mol_s: rate,
                    maximum_mol_s: self.controls.maximum_node_charge_rate_mol_s,
                });
            }
            maximum_node_rate = maximum_node_rate.max(rate.abs());
        }
        self.admit_and_publish(cells, links, n, controls, cancelled)?;
        self.report.maximum_node_charge_rate_mol_s = self
            .report
            .maximum_node_charge_rate_mol_s
            .max(maximum_node_rate);
        self.report.maximum_edge_charge_rate_mol_s = self
            .report
            .maximum_edge_charge_rate_mol_s
            .max(maximum_edge_rate);
        self.report.maximum_dimensionless_potential_difference = self
            .report
            .maximum_dimensionless_potential_difference
            .max(maximum_potential);
        Ok(())
    }

    fn solve_potential(
        &mut self,
        faces: &mut [Face<'_>],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), OsmoticNetworkError> {
        self.rhs.fill(0.0);
        self.diagonal.fill(0.0);
        self.charge_correction.fill(0.0);
        self.potentials.fill(0.0);
        for (i, parent) in self.parents.iter_mut().enumerate() {
            *parent = i;
        }
        for face in faces.iter() {
            poll(cancelled)?;
            let [left, right] = face.endpoints;
            let chemical = finite(face.area_m2 * face.response.chemical_charge_flux_mol_m2_s())?;
            if face.response.chemical_charge_flux_mol_m2_s() != 0.0 && chemical == 0.0 {
                return Err(OsmoticNetworkError::Invalid(
                    "chemical edge current is unrepresentable",
                ));
            }
            compensated_add(
                &mut self.rhs[left],
                &mut self.charge_correction[left],
                -chemical,
            );
            compensated_add(
                &mut self.rhs[right],
                &mut self.charge_correction[right],
                chemical,
            );
            if face.weight == 0.0 {
                if chemical != 0.0 {
                    return Err(OsmoticNetworkError::Invalid(
                        "current without ionic conductance",
                    ));
                }
                continue;
            }
            self.diagonal[left] = finite(self.diagonal[left] + face.weight)?;
            self.diagonal[right] = finite(self.diagonal[right] + face.weight)?;
            let (a, b) = (
                root(&mut self.parents, left),
                root(&mut self.parents, right),
            );
            self.parents[a.max(b)] = a.min(b);
        }
        self.report.potential_solves = self
            .report
            .potential_solves
            .checked_add(1)
            .ok_or(OsmoticNetworkError::WorkBudget)?;
        // Independent components may have entirely different conductance
        // scales. Normalize each block separately; an unrelated strong brine
        // must not underflow a dilute component's otherwise well-posed operator.
        self.charge_rate.fill(0.0);
        for node in 0..self.parents.len() {
            let component = root(&mut self.parents, node);
            self.parents[node] = component;
            self.charge_rate[component] = self.charge_rate[component].max(self.diagonal[node]);
        }
        for face in faces.iter_mut() {
            let original = face.weight;
            if original == 0.0 {
                continue;
            }
            face.weight /= self.charge_rate[self.parents[face.endpoints[0]]];
            if original != 0.0 && !positive(face.weight) {
                return Err(OsmoticNetworkError::Invalid(
                    "scaled ionic conductance is unrepresentable",
                ));
            }
        }
        for node in 0..self.rhs.len() {
            if self.parents[node] == node {
                self.rhs[node] = 0.0;
                self.diagonal[node] = 1.0;
            } else {
                let scale = self.charge_rate[self.parents[node]];
                self.rhs[node] = finite(self.rhs[node] / scale)?;
                self.diagonal[node] /= scale;
                if !positive(self.diagonal[node]) {
                    return Err(OsmoticNetworkError::Invalid(
                        "scaled ionic diagonal is unrepresentable",
                    ));
                }
            }
        }
        if self.rhs.iter().all(|&value| value == 0.0) {
            return Ok(());
        }
        let remaining = self
            .controls
            .maximum_potential_iterations
            .saturating_sub(self.report.potential_iterations);
        if remaining == 0 {
            return Err(OsmoticNetworkError::WorkBudget);
        }
        let graph = Graph {
            faces,
            parents: &self.parents,
        };
        let diagonal = Jacobi(&self.diagonal);
        let mut previous = 0;
        let solved = checked_cg(
            &graph,
            &diagonal,
            &self.rhs,
            CheckedCgConfig {
                tolerance: self.controls.potential_relative_tolerance,
                max_iterations: remaining.min(self.controls.maximum_iterations_per_solve),
                max_corrections: 2,
            },
            |iterations| {
                self.report.potential_iterations += iterations - previous;
                previous = iterations;
                poll(cancelled)
            },
        )
        .map_err(|error| match error {
            CheckedCgError::Interrupted(error) => error,
            CheckedCgError::InvalidInput(what) => OsmoticNetworkError::Invalid(what),
        })?;
        if !solved.report.converged_euclidean() {
            if self.report.potential_iterations >= self.controls.maximum_potential_iterations {
                return Err(OsmoticNetworkError::WorkBudget);
            }
            return Err(OsmoticNetworkError::ElectroneutralPotentialConvergence {
                relative_residual: solved.report.rel_residual,
            });
        }
        self.potentials.copy_from_slice(&solved.x);
        Ok(())
    }

    fn admit_and_publish(
        &mut self,
        cells: &mut [OsmoticCompartment],
        links: &[SolutionDiffusionLink],
        n: usize,
        controls: SolutionDiffusionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), OsmoticNetworkError> {
        let ionic = links[0]
            .electrical
            .as_deref()
            .ok_or(OsmoticNetworkError::Invalid(
                "missing original ionic policy",
            ))?;
        let mut balance = [0.0; MAX_OSMOTIC_SPECIES];
        let mut corrections = [0.0; MAX_OSMOTIC_SPECIES];
        let mut energy = 0.0;
        let mut energy_correction = 0.0;
        let mut proposed_change = false;
        let mut actual_change = false;
        for (node, cell) in cells.iter().enumerate() {
            poll(cancelled)?;
            let old = &cell.solution;
            let volume = old.solvent_volume_m3();
            let next = &mut self.amounts[node * n..(node + 1) * n];
            for i in 0..n {
                let amount = old.solute_moles()[i];
                let withdrawal = finite(self.withdrawn[node * n + i])?;
                let fraction = if withdrawal == 0.0 {
                    0.0
                } else {
                    withdrawal / amount
                };
                if !fraction.is_finite() || fraction > controls.maximum_species_fraction {
                    return Err(OsmoticNetworkError::DiffusionStep(
                        SolutionDiffusionStepRefusal::Withdrawal {
                            species: i,
                            fraction,
                            maximum: controls.maximum_species_fraction,
                        },
                    ));
                }
                proposed_change |= next[i] != 0.0;
                next[i] = finite(amount + next[i])?;
                if !nonnegative(next[i]) {
                    return Err(OsmoticNetworkError::Invalid(
                        "nonnegative ionic endpoint is unrepresentable",
                    ));
                }
                actual_change |= next[i] != amount;
                compensated_add(&mut balance[i], &mut corrections[i], next[i] - amount);
            }
            // Actual rounded inventories, including isolated nodes, are checked;
            // no counterion repair or background-relative charge allowance.
            ionic.admit_inventories(next).map_err(diffusion_source)?;
            old.osmotic_pressure_at(volume, next).map_err(source)?;
            let law = old.nonideal_model().ok_or(OsmoticNetworkError::Invalid(
                "missing original ionic free energy",
            ))?;
            let mut before = [0.0; MAX_OSMOTIC_SPECIES];
            let mut after = [0.0; MAX_OSMOTIC_SPECIES];
            old.concentrations_into(&mut before[..n]).map_err(source)?;
            for i in 0..n {
                after[i] = next[i] / volume;
            }
            let change = law
                .reference_cancelled_free_energy_change_j_m3(&before[..n], &after[..n])
                .map_err(fs_material::osmosis::OsmosisError::from)
                .map_err(source)?;
            compensated_add(
                &mut energy,
                &mut energy_correction,
                finite(volume * change)?,
            );
        }
        if proposed_change && !actual_change {
            return Err(OsmoticNetworkError::StepUnderflow);
        }
        for (species, &residual_mol) in balance[..n].iter().enumerate() {
            if !residual_mol.is_finite()
                || residual_mol.abs() > controls.amount_balance_tolerance_mol
            {
                return Err(OsmoticNetworkError::DiffusionBalance {
                    species,
                    residual_mol,
                });
            }
        }
        finite(energy)?;
        if energy > controls.free_energy_tolerance_j {
            return Err(OsmoticNetworkError::DiffusionStep(
                SolutionDiffusionStepRefusal::FreeEnergy {
                    change_j: energy,
                    maximum_j: controls.free_energy_tolerance_j,
                },
            ));
        }
        for link in links {
            poll(cancelled)?;
            link.qualify_sources(
                [
                    cells[link.left].solution.solvent_volume_m3(),
                    cells[link.right].solution.solvent_volume_m3(),
                ],
                [
                    &self.amounts[link.left * n..(link.left + 1) * n],
                    &self.amounts[link.right * n..(link.right + 1) * n],
                ],
                cancelled,
            )?;
        }
        // The enclosing owner publishes only after the entire requested window.
        for (node, cell) in cells.iter_mut().enumerate() {
            poll(cancelled)?;
            cell.solution
                .replace_inventories(
                    cell.solution.solvent_volume_m3(),
                    &self.amounts[node * n..(node + 1) * n],
                )
                .map_err(source)?;
        }
        Ok(())
    }
}

fn root(parents: &mut [usize], mut node: usize) -> usize {
    while parents[node] != node {
        parents[node] = parents[parents[node]];
        node = parents[node];
    }
    node
}

/// Gauge rows are identity; their columns are removed from all other rows.
/// Positive edge weights and one gauge per conducting component imply SPD.
struct Graph<'a, 'b> {
    faces: &'a [Face<'b>],
    parents: &'a [usize],
}
impl LinearOp for Graph<'_, '_> {
    fn n(&self) -> usize {
        self.parents.len()
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        y.fill(0.0);
        for face in self.faces {
            let [left, right] = face.endpoints;
            let a = if self.parents[left] == left {
                0.0
            } else {
                x[left]
            };
            let b = if self.parents[right] == right {
                0.0
            } else {
                x[right]
            };
            let value = face.weight * (a - b);
            if self.parents[left] != left {
                y[left] += value;
            }
            if self.parents[right] != right {
                y[right] -= value;
            }
        }
        for (node, &parent) in self.parents.iter().enumerate() {
            if parent == node {
                y[node] = x[node];
            }
        }
    }
}
struct Jacobi<'a>(&'a [f64]);
impl Precond for Jacobi<'_> {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        for ((zi, &ri), &di) in z.iter_mut().zip(r).zip(self.0) {
            *zi = ri / di;
        }
    }
}
