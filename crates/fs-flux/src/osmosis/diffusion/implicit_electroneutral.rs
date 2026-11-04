//! Globally coupled backward Euler on the original electroneutral network.
//!
//! The shared ionic coordinate chart preserves each node's original charge and
//! every molecular component's species inventories during nonlinear probes.
//! Every residual re-evaluates the original material faces, global potential
//! and combined-force fluxes. Native Newton/finite differences own nonlinear
//! iteration; the existing electrical owner owns all physical transport.

use super::ionic_coordinates::IonicCoordinates;
use super::{
    ElectroneutralDiffusion, ImplicitMulticomponentDiffusionControls, MAX_COMPARTMENTS,
    MAX_OSMOTIC_SPECIES, OsmoticCompartment, OsmoticNetwork, OsmoticNetworkError,
    SolutionDiffusionControls, SolutionDiffusionLink, finite, poll, positive, source,
};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{NewtonKrylovState, NewtonStallDiagnosis, SolverRunProgress};
use std::cell::RefCell;

/// Original physical amount accuracy and native nonlinear controls, plus a
/// bound on the complete coupled network's independent coordinates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImplicitElectroneutralDiffusionControls {
    /// The original multicomponent nonlinear policy. Its molecular evaluation
    /// budget counts actual face preparations/replays across the whole window,
    /// including every Jacobian stencil, line search and discarded time probe.
    pub nonlinear: ImplicitMulticomponentDiffusionControls,
    /// Maximum independent amount coordinates across all molecular components.
    /// This is not a limit on the original species equations checked at replay.
    pub maximum_unknowns: usize,
}

impl Default for ImplicitElectroneutralDiffusionControls {
    fn default() -> Self {
        Self {
            nonlinear: ImplicitMulticomponentDiffusionControls::default(),
            maximum_unknowns: 4096,
        }
    }
}

impl ImplicitElectroneutralDiffusionControls {
    fn admit(self) -> Result<(), OsmoticNetworkError> {
        self.nonlinear.admit()?;
        if self.maximum_unknowns == 0
            || self.maximum_unknowns > MAX_COMPARTMENTS * MAX_OSMOTIC_SPECIES
        {
            return Err(OsmoticNetworkError::Invalid(
                "global implicit ionic coordinates need a bounded positive unknown count",
            ));
        }
        Ok(())
    }

    fn solver_workspace_bytes(self, dimension: usize) -> Result<usize, OsmoticNetworkError> {
        // The shared native policy caps restart, cycles and outer attempts at
        // the original 64-coordinate owner's limits. Its conservative payload
        // estimate therefore scales in blocks of 64 unknowns; no dense global
        // Jacobian or second numerical solver is retained here.
        dimension
            .max(1)
            .div_ceil(MAX_OSMOTIC_SPECIES)
            .checked_mul(ImplicitMulticomponentDiffusionControls::solver_workspace_bytes())
            .ok_or(OsmoticNetworkError::Invalid(
                "global ionic native solver storage overflow",
            ))
    }
}

#[derive(Clone, Debug)]
pub(crate) struct ImplicitElectroneutralDiffusion {
    controls: ImplicitElectroneutralDiffusionControls,
    trial: Vec<OsmoticCompartment>,
    volumes: Vec<f64>,
    original_amounts: Vec<f64>,
    candidate_amounts: Vec<f64>,
    concentration_caps: Vec<f64>,
    component_ids: Vec<usize>,
    evaluations: usize,
}

impl OsmoticNetwork {
    /// Select simultaneous backward-Euler ionic transport on an already
    /// prepared global electrical network, including cycles and independent
    /// salts. Every nonlinear trial remains charge- and species-constrained.
    ///
    /// The original full/two-half timeline still controls temporal error.
    /// Implicit stability does not waive original gross withdrawal, rounded
    /// inventory charge, species balance, material support or global free-energy
    /// gates. All original node/species equations are replayed before admission,
    /// including eliminated anchors/counterions and exactly absent spectators.
    /// Source, nonlinear and cancellation failures publish no physical state.
    ///
    /// Storage includes reusable original solution buffers, the bounded ionic
    /// chart and a dimension-scaled allowance for native Newton/FD/FGMRES work.
    /// This is a numerical payload estimate, not an allocator or RSS bound;
    /// original material query/receipt allocations keep their existing scope.
    pub fn with_implicit_electroneutral_diffusion(
        mut self,
        controls: ImplicitElectroneutralDiffusionControls,
        maximum_workspace_bytes: usize,
    ) -> Result<Self, OsmoticNetworkError> {
        controls.admit()?;
        if self.electroneutral_diffusion.is_none()
            || self.implicit_diffusion.is_some()
            || self.implicit_multicomponent.is_some()
            || self.implicit_salt.is_some()
        {
            return Err(OsmoticNetworkError::Invalid(
                "global ionic backward Euler requires the original global electrical closure",
            ));
        }
        if self.implicit_electroneutral.is_some() {
            return Err(OsmoticNetworkError::Invalid(
                "global ionic backward Euler is already selected",
            ));
        }
        let count = self.species.len();
        let entries = self
            .cells
            .len()
            .checked_mul(count)
            .ok_or(OsmoticNetworkError::Invalid(
                "global ionic inventory storage overflow",
            ))?;
        let buffers = entries
            .checked_mul(3 * core::mem::size_of::<f64>())
            .and_then(|bytes| {
                bytes.checked_add(self.cells.len().checked_mul(
                    core::mem::size_of::<f64>()
                        + core::mem::size_of::<usize>()
                        + core::mem::size_of::<OsmoticCompartment>(),
                )?)
            });
        let solutions = self.cells.iter().try_fold(0usize, |bytes, cell| {
            bytes.checked_add(cell.solution.heap_bytes())
        });
        let base_required = buffers
            .and_then(|bytes| bytes.checked_add(solutions?))
            .and_then(|bytes| {
                bytes.checked_add(core::mem::size_of::<ImplicitElectroneutralDiffusion>())
            })
            .and_then(|bytes| bytes.checked_add(self.workspace_bytes))
            .ok_or(OsmoticNetworkError::Invalid(
                "global ionic workspace byte count overflow",
            ))?;
        if base_required > maximum_workspace_bytes {
            return Err(OsmoticNetworkError::Workspace {
                required_bytes: base_required,
                budget_bytes: maximum_workspace_bytes,
            });
        }
        let component_ids = molecular_components(self.cells.len(), &self.diffusion);
        let components = component_ids
            .iter()
            .enumerate()
            .filter(|&(i, root)| i == *root)
            .count();
        // Water exchange can release an originally saturated species. Reserve
        // the topology's admitted coordinate bound, not just its initial rank.
        let admitted_dimension =
            ((self.cells.len() - components) * count).min(controls.maximum_unknowns);
        let chart_bytes =
            IonicCoordinates::payload_bound(self.cells.len(), count, admitted_dimension)?;
        let required = base_required
            .checked_add(chart_bytes)
            .and_then(|bytes| {
                bytes.checked_add(controls.solver_workspace_bytes(admitted_dimension).ok()?)
            })
            .ok_or(OsmoticNetworkError::Invalid(
                "global ionic workspace byte count overflow",
            ))?;
        if required > maximum_workspace_bytes {
            return Err(OsmoticNetworkError::Workspace {
                required_bytes: required,
                budget_bytes: maximum_workspace_bytes,
            });
        }
        let volumes: Vec<_> = self
            .cells
            .iter()
            .map(|cell| cell.solution.solvent_volume_m3())
            .collect();
        let original_amounts: Vec<_> = self
            .cells
            .iter()
            .flat_map(|cell| cell.solution.solute_moles().iter().copied())
            .collect();
        let mut concentration_caps = Vec::with_capacity(entries);
        for cell in &self.cells {
            let law = cell
                .solution
                .nonideal_model()
                .ok_or(OsmoticNetworkError::Invalid(
                    "global ionic backward Euler requires original nonideal thermodynamics",
                ))?;
            concentration_caps.extend_from_slice(law.maximum_concentrations_mol_m3());
        }
        let ionic = self.diffusion[0]
            .electrical
            .as_deref()
            .ok_or(OsmoticNetworkError::Invalid(
                "global ionic backward Euler requires original valences",
            ))?;
        IonicCoordinates::new(
            &volumes,
            &original_amounts,
            &concentration_caps,
            ionic,
            &component_ids,
            controls.maximum_unknowns,
            &mut || false,
        )?;
        self.implicit_electroneutral = Some(ImplicitElectroneutralDiffusion {
            controls,
            trial: self.cells.clone(),
            volumes,
            candidate_amounts: original_amounts.clone(),
            original_amounts,
            concentration_caps,
            component_ids,
            evaluations: 0,
        });
        self.workspace_bytes = required;
        Ok(self)
    }

    pub fn implicit_electroneutral_diffusion_controls(
        &self,
    ) -> Option<ImplicitElectroneutralDiffusionControls> {
        self.implicit_electroneutral
            .as_ref()
            .map(|state| state.controls)
    }
}

struct Candidate {
    scaled: Vec<f64>,
    refusal: Option<OsmoticNetworkError>,
    changed: bool,
    nonzero_transfer: bool,
}

impl ImplicitElectroneutralDiffusion {
    pub(crate) fn begin_window(&mut self) {
        self.evaluations = 0;
    }

    #[allow(clippy::too_many_arguments)]
    fn evaluate(
        &mut self,
        chart: &IonicCoordinates,
        x: &[f64],
        electrical: &mut ElectroneutralDiffusion,
        links: &mut [SolutionDiffusionLink],
        dt: f64,
        controls: SolutionDiffusionControls,
        evaluations: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Candidate, OsmoticNetworkError> {
        poll(cancelled)?;
        chart.decode_into(x, &mut self.candidate_amounts, cancelled)?;
        let n = self.original_amounts.len() / self.trial.len();
        for (node, cell) in self.trial.iter_mut().enumerate() {
            poll(cancelled)?;
            cell.solution
                .replace_inventories(
                    self.volumes[node],
                    &self.candidate_amounts[node * n..(node + 1) * n],
                )
                .map_err(source)?;
        }
        // A global physical evaluation prepares and replays each active face.
        // Reserve its complete known face work before starting source/CG work;
        // still charge the actual spent prefix when cancellation/refusal occurs.
        let requested = links
            .iter()
            .filter(|link| link.area_m2 != 0.0)
            .count()
            .checked_mul(2)
            .ok_or(OsmoticNetworkError::WorkBudget)?;
        if requested
            > self
                .controls
                .nonlinear
                .maximum_evaluations
                .saturating_sub(self.evaluations)
        {
            return Err(OsmoticNetworkError::WorkBudget);
        }
        let before = *evaluations;
        let physical =
            electrical.evaluate_transfers(&self.trial, links, dt, controls, evaluations, cancelled);
        self.evaluations = self
            .evaluations
            .checked_add(
                evaluations
                    .checked_sub(before)
                    .ok_or(OsmoticNetworkError::WorkBudget)?,
            )
            .ok_or(OsmoticNetworkError::WorkBudget)?;
        physical?;
        let transfers = electrical.transfers_mol();
        let mut refusal = None;
        let mut changed = false;
        let mut nonzero_transfer = false;
        for (i, &transfer) in transfers.iter().enumerate() {
            poll(cancelled)?;
            let change = finite(self.candidate_amounts[i] - self.original_amounts[i])?;
            let residual = finite(change - transfer)?;
            let maximum = self.controls.nonlinear.tolerance(change)?;
            if refusal.is_none() && residual.abs() > maximum {
                refusal = Some(OsmoticNetworkError::ImplicitDiffusionConvergence {
                    residual_mol: residual,
                    maximum_mol: maximum,
                });
            }
            changed |= change != 0.0;
            nonzero_transfer |= transfer != 0.0;
        }
        let scaled = chart
            .independent_indices()
            .iter()
            .map(|&i| {
                finite(
                    (self.candidate_amounts[i] - self.original_amounts[i] - transfers[i])
                        / self.controls.nonlinear.absolute_residual_tolerance_mol,
                )
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Candidate {
            scaled,
            refusal,
            changed,
            nonzero_transfer,
        })
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn sweep(
        &mut self,
        electrical: &mut ElectroneutralDiffusion,
        cells: &mut [OsmoticCompartment],
        links: &mut [SolutionDiffusionLink],
        dt: f64,
        controls: SolutionDiffusionControls,
        evaluations: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), OsmoticNetworkError> {
        poll(cancelled)?;
        if !positive(dt) || cells.len() != self.trial.len() {
            return Err(OsmoticNetworkError::Invalid(
                "global ionic backward Euler needs original nodes and a finite positive duration",
            ));
        }
        let n = cells[0].solution.solute_moles().len();
        for (node, cell) in cells.iter().enumerate() {
            self.volumes[node] = cell.solution.solvent_volume_m3();
            self.original_amounts[node * n..(node + 1) * n]
                .copy_from_slice(cell.solution.solute_moles());
        }
        let ionic = links[0]
            .electrical
            .as_deref()
            .ok_or(OsmoticNetworkError::Invalid(
                "global ionic backward Euler lost its original valences",
            ))?;
        let chart = IonicCoordinates::new(
            &self.volumes,
            &self.original_amounts,
            &self.concentration_caps,
            ionic,
            &self.component_ids,
            self.controls.maximum_unknowns,
            cancelled,
        )?;
        if chart.dimension() == 0 {
            let candidate = self.evaluate(
                &chart,
                &[],
                electrical,
                links,
                dt,
                controls,
                evaluations,
                cancelled,
            )?;
            if let Some(refusal) = candidate.refusal {
                return Err(refusal);
            }
            if candidate.nonzero_transfer && !candidate.changed {
                return Err(OsmoticNetworkError::StepUnderflow);
            }
            return electrical.publish_candidate(
                cells,
                links,
                &self.candidate_amounts,
                controls,
                cancelled,
            );
        }
        let numerical = self.controls.nonlinear;
        let owner = RefCell::new((self, electrical, links, evaluations));
        let cancel = RefCell::new(cancelled);
        let evaluate = |x: &[f64]| {
            let mut owner = owner.borrow_mut();
            let (state, electrical, links, evaluations) = &mut *owner;
            state.evaluate(
                &chart,
                x,
                electrical,
                links,
                dt,
                controls,
                evaluations,
                &mut || (*cancel.borrow_mut())(),
            )
        };
        let problem = FallibleFiniteDifference::new(
            chart.dimension(),
            numerical.difference_step,
            |x: &[f64]| evaluate(x).map(|candidate| candidate.scaled),
        )
        .map_err(numerical_error)?;
        let construction =
            NewtonKrylovState::new(&problem, chart.initial_guess().to_vec(), numerical.newton);
        if let Some(error) = problem.take_failure() {
            return Err(numerical_error(error));
        }
        let mut solver = construction.map_err(|_| {
            OsmoticNetworkError::Invalid("native global ionic Newton initialization refused")
        })?;
        let mut terminal = None;
        for attempt in 0..=numerical.maximum_iterations {
            // This fresh checkpoint replay is also the exact source/current
            // snapshot whose gross withdrawals the final physical gate checks.
            let candidate = evaluate(&solver.x)?;
            if candidate.refusal.is_none() {
                if candidate.nonzero_transfer && !candidate.changed {
                    return Err(OsmoticNetworkError::StepUnderflow);
                }
                poll(&mut || (*cancel.borrow_mut())())?;
                let mut owner = owner.borrow_mut();
                let (state, electrical, links, _) = &mut *owner;
                return electrical.publish_candidate(
                    cells,
                    links,
                    &state.candidate_amounts,
                    controls,
                    &mut || (*cancel.borrow_mut())(),
                );
            }
            if terminal == Some(None) {
                return Err(candidate
                    .refusal
                    .expect("failed original global ion equations"));
            }
            if terminal.is_some() || attempt == numerical.maximum_iterations {
                return Err(OsmoticNetworkError::ImplicitMulticomponentConvergence {
                    scaled_residual: solver.residual_norm(),
                    diagnosis: terminal
                        .flatten()
                        .or(Some(NewtonStallDiagnosis::BudgetExhausted)),
                });
            }
            let result = solver.run_controlled(&problem, 1, &mut || {
                (*cancel.borrow_mut())() || problem.has_failed()
            });
            if let Some(error) = problem.take_failure() {
                return Err(numerical_error(error));
            }
            match result.progress {
                SolverRunProgress::Complete => {}
                SolverRunProgress::Paused => return Err(OsmoticNetworkError::Cancelled),
                _ => {
                    return Err(OsmoticNetworkError::Invalid(
                        "native global ionic nonlinear callback interrupted",
                    ));
                }
            }
            if result.report.converged
                || result.report.diagnosis != Some(NewtonStallDiagnosis::BudgetExhausted)
            {
                terminal = Some(result.report.diagnosis);
            }
        }
        unreachable!("bounded global Newton checkpoint loop returns at its final attempt")
    }
}

fn numerical_error(error: FiniteDifferenceError<OsmoticNetworkError>) -> OsmoticNetworkError {
    match error {
        FiniteDifferenceError::Source(error) => error,
        FiniteDifferenceError::Invalid(_) => OsmoticNetworkError::Invalid(
            "native global ionic numerical Jacobian refused its finite stencil",
        ),
    }
}

fn molecular_components(nodes: usize, links: &[SolutionDiffusionLink]) -> Vec<usize> {
    fn root(parents: &mut [usize], mut node: usize) -> usize {
        while parents[node] != node {
            parents[node] = parents[parents[node]];
            node = parents[node];
        }
        node
    }
    let mut parents: Vec<_> = (0..nodes).collect();
    for link in links.iter().filter(|link| link.area_m2 != 0.0) {
        let (a, b) = (
            root(&mut parents, link.left),
            root(&mut parents, link.right),
        );
        parents[a.max(b)] = a.min(b);
    }
    for node in 0..nodes {
        parents[node] = root(&mut parents, node);
    }
    parents
}
