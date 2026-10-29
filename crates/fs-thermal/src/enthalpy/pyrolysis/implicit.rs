//! Backward-Euler conduction and decomposition on the original thermal graph.
//!
//! Cell absorbed energies are the nonlinear unknowns. The original `solve_cell`
//! maps each trial energy to its thermally coupled reaction endpoint, including
//! virgin/char storage and outgoing formation enthalpy. The SAME original
//! heat-rate operator evaluates that endpoint; the remaining equation is
//! `absorbed_i - dt * heat_i(T_endpoint) = 0`. No second heat-of-reaction or
//! alternate conductive graph is introduced.
//!
//! A prepared original dense LU solves the bounded Newton system. Its Jacobian
//! uses the original link/Robin conductances and a one-cell caloric derivative;
//! it does not rerun the entire graph for every matrix column. Matrix/vector
//! storage is reused. The original heat owner still allocates its bounded rate
//! maps, and candidate/result storage is allocated per advance.
//!
//! The explicit conduction CFL no longer limits accepted intervals. This is
//! nevertheless first-order time integration, not an accuracy, unconditional
//! nonlinear-convergence, or entropy theorem. Original material/source domains,
//! finite-rate kinetics and physical depletion refusals remain in force.

use super::*;
use fs_la::prepared_lu::{LuWorkspace, LuWorkspaceError};

/// Bounded small dense solve; larger graphs need a separate sparse solver.
pub const MAX_IMPLICIT_CELLS: usize = 128;

/// Nonlinear work controls, independent of the original physical/time budgets.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImplicitSolve {
    /// Endpoint residual evaluations per interval, 1 through 64.
    pub maximum_iterations: usize,
    /// Backtracking trials per Newton update, 1 through 32.
    pub maximum_backtracks: usize,
    /// Minimum one-cell energy probe [J] for `dT / d(absorbed energy)`.
    /// The probe is enlarged to resolve the requested local algebraic accuracy.
    pub energy_difference_j: f64,
}
impl Default for ImplicitSolve {
    fn default() -> Self {
        Self {
            maximum_iterations: 24,
            maximum_backtracks: 24,
            energy_difference_j: 1e-4,
        }
    }
}

/// Numerical scratch only; no physical history or material binding is retained.
pub struct ImplicitWorkspace {
    ids: Vec<CellId>,
    energy: Vec<f64>,
    probe: Vec<f64>,
    residual: Vec<f64>,
    slope: Vec<f64>,
    conductance: Vec<f64>,
    matrix: Vec<f64>,
    rhs: Vec<f64>,
    delta: Vec<f64>,
    steps: Vec<CellStep>,
    lu: LuWorkspace,
}

/// Original conservation account plus actual nonlinear work.
#[derive(Clone, Debug, PartialEq)]
pub struct ImplicitReport {
    /// Original solid/volatile and heat-storage accounts, over accepted steps.
    pub pyrolysis: PyrolysisAudit,
    /// Unperturbed endpoint residual evaluations across all accepted intervals.
    pub nonlinear_iterations: usize,
    /// Original graph heat-rate evaluations, including discarded line probes.
    pub heat_evaluations: usize,
    /// Calls to the original local reaction/temperature solver, including probes.
    pub local_cell_solves: usize,
    /// Original prepared linear solves.
    pub linear_solves: usize,
    /// Backtracking probes evaluated, including rejected trials.
    pub line_probes: usize,
    /// Sum of absolute differences between the original graph's summed cell
    /// input and its independent external input [J]. This floating-point
    /// transfer discrepancy is included in the original energy budget.
    pub absolute_heat_balance_roundoff_j: f64,
}

/// An error publishes neither cells, inventory, releases, nor physical time.
#[derive(Clone, Debug, PartialEq)]
pub enum ImplicitError {
    /// Original physical, source, representation, or work-budget refusal.
    Thermal(DepositionError),
    /// Original prepared linear solve refusal.
    Linear(LuWorkspaceError),
    /// Bounded Newton/backtracking failed to close the actual endpoint residual.
    Iterations,
    /// Original caller cancellation, preserved as a direct path.
    Cancelled,
}
impl core::fmt::Display for ImplicitError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "implicit pyrolysis: {self:?}")
    }
}
impl std::error::Error for ImplicitError {}
impl From<DepositionError> for ImplicitError {
    fn from(error: DepositionError) -> Self {
        if error == DepositionError::Cancelled {
            Self::Cancelled
        } else {
            Self::Thermal(error)
        }
    }
}
impl From<LuWorkspaceError> for ImplicitError {
    fn from(error: LuWorkspaceError) -> Self {
        Self::Linear(error)
    }
}

#[derive(Default)]
struct Work {
    iterations: usize,
    heat: usize,
    cells: usize,
    temperature_evaluations: usize,
    linear: usize,
    line_probes: usize,
}

#[derive(Default)]
struct Account {
    energy: Sum,
    absolute_energy: Sum,
    stored: Sum,
    mass: Sum,
    external: Sum,
    kinetic_mass: Sum,
    kinetic_energy: Sum,
    heat_roundoff: Sum,
}

impl PyrolysisSimulation<'_> {
    /// Allocate numerical storage once for these stable cell-slot identities.
    /// Reuse across source/boundary/conductance or material-state changes is safe:
    /// every physical coefficient and guess is rebuilt from the actual state.
    pub fn implicit_workspace(&self) -> Result<ImplicitWorkspace, ImplicitError> {
        let n = self.cells().len();
        if !(1..=MAX_IMPLICIT_CELLS).contains(&n) {
            return Err(DepositionError::Budget("implicit pyrolysis cell count").into());
        }
        Ok(ImplicitWorkspace {
            ids: self.cells().keys().copied().collect(),
            energy: vec![0.0; n],
            probe: vec![0.0; n],
            residual: vec![0.0; n],
            slope: vec![0.0; n],
            conductance: vec![0.0; n * n],
            matrix: vec![0.0; n * n],
            rhs: vec![0.0; n],
            delta: vec![0.0; n],
            steps: Vec::with_capacity(n),
            lu: LuWorkspace::new(n)?,
        })
    }
}

impl<'a> PyrolysisSimulation<'a> {
    /// Advance implicit conduction and original thermally coupled decomposition.
    ///
    /// `maximum_dt_s` is an accuracy input, not a conduction-CFL ceiling. No
    /// automatic halving hides nonlinear/material refusal. Every accepted
    /// interval closes local and aggregate physical energy and kinetic-roundoff
    /// budgets; cancellation/refusal leaves `self` unchanged. Geometry and
    /// authored coefficients are frozen over each call. A moving-slab owner
    /// must rezone and resolve its coefficients between these advances.
    pub fn advance_implicit_to(
        &self,
        target_s: f64,
        controls: PyrolysisControls,
        solve: ImplicitSolve,
        workspace: &mut ImplicitWorkspace,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, ImplicitReport), ImplicitError> {
        checkpoint()?;
        if !target_s.is_finite()
            || target_s < self.time_s()
            || !positive(controls.maximum_dt_s)
            || !(4..=256).contains(&controls.maximum_temperature_iterations)
            || !nonnegative(controls.cell_energy_tolerance_j)
            || !nonnegative(controls.total_energy_tolerance_j)
            || !nonnegative(controls.mass_tolerance_kg)
            || !(1..=64).contains(&solve.maximum_iterations)
            || !(1..=32).contains(&solve.maximum_backtracks)
            || !positive(solve.energy_difference_j)
            || workspace.ids.len() != self.cells().len()
            || !workspace.ids.iter().eq(self.cells().keys())
        {
            return Err(
                DepositionError::InvalidInput("implicit pyrolysis controls or cell axis").into(),
            );
        }
        prepare_conductance(self, workspace, &mut checkpoint)?;
        let mut next = self.clone();
        let mut audit = PyrolysisAudit {
            thermal: empty_thermal_audit(),
            releases: self
                .cells()
                .keys()
                .map(|&id| (id, VolatileRelease::default()))
                .collect(),
            volatile_mass_kg: 0.0,
            volatile_enthalpy_j: 0.0,
            mass_residual_kg: 0.0,
            maximum_cell_energy_residual_j: 0.0,
            absolute_energy_residual_j: 0.0,
            kinetic_mass_roundoff_kg: 0.0,
            kinetic_energy_roundoff_j: 0.0,
            temperature_evaluations: 0,
        };
        let mut account = Account::default();
        let mut work = Work::default();
        while next.time_s() < target_s {
            checkpoint()?;
            if audit.thermal.substeps >= self.thermal.budget.max_substeps {
                return Err(DepositionError::Budget("implicit pyrolysis intervals").into());
            }
            let initial = next.clone();
            let remaining = target_s - initial.time_s();
            let requested = controls.maximum_dt_s.min(remaining);
            let endpoint = if requested == remaining {
                target_s
            } else {
                initial.time_s() + requested
            };
            let dt = endpoint - initial.time_s();
            if !positive(dt) || !endpoint.is_finite() {
                return Err(DepositionError::Arithmetic.into());
            }
            initialize(&initial, dt, workspace, &mut checkpoint)?;
            let external_w = nonlinear_step(
                &initial,
                &mut next,
                dt,
                controls,
                solve,
                workspace,
                &mut work,
                &mut checkpoint,
            )?;
            accept_interval(
                &initial,
                &next,
                dt,
                external_w,
                controls,
                workspace,
                &mut audit,
                &mut account,
                &mut checkpoint,
            )?;
            next.thermal.time_s = endpoint;
            audit.thermal.substeps += 1;
        }
        audit.thermal.external_heat_j = account.external.value();
        audit.thermal.stored_energy_change_j = account.stored.value();
        audit.thermal.residual_j = account.energy.value();
        audit.absolute_energy_residual_j = account.absolute_energy.value();
        audit.mass_residual_kg = account.mass.value();
        audit.kinetic_mass_roundoff_kg = account.kinetic_mass.value();
        audit.kinetic_energy_roundoff_j = account.kinetic_energy.value();
        audit.temperature_evaluations = work.temperature_evaluations;
        let mut emitted_mass = Sum::default();
        let mut emitted_enthalpy = Sum::default();
        for release in audit.releases.values() {
            emitted_mass.add(release.mass_kg)?;
            emitted_enthalpy.add(release.enthalpy_j)?;
        }
        audit.volatile_mass_kg = emitted_mass.value();
        audit.volatile_enthalpy_j = emitted_enthalpy.value();
        if audit.absolute_energy_residual_j > controls.total_energy_tolerance_j {
            return Err(DepositionError::EnergyBalance(audit.absolute_energy_residual_j).into());
        }
        if audit.absolute_energy_residual_j + audit.kinetic_energy_roundoff_j
            > controls.total_energy_tolerance_j
        {
            return Err(DepositionError::Budget(
                "implicit pyrolysis accumulated reaction-energy roundoff",
            )
            .into());
        }
        if audit.mass_residual_kg.abs() > controls.mass_tolerance_kg {
            return Err(DepositionError::InvalidInput(
                "implicit pyrolysis solid/volatile mass balance",
            )
            .into());
        }
        if audit.mass_residual_kg.abs() + audit.kinetic_mass_roundoff_kg
            > controls.mass_tolerance_kg
        {
            return Err(DepositionError::Budget(
                "implicit pyrolysis accumulated kinetic inventory roundoff",
            )
            .into());
        }
        checkpoint()?;
        Ok((
            next,
            ImplicitReport {
                pyrolysis: audit,
                nonlinear_iterations: work.iterations,
                heat_evaluations: work.heat,
                local_cell_solves: work.cells,
                linear_solves: work.linear,
                line_probes: work.line_probes,
                absolute_heat_balance_roundoff_j: account.heat_roundoff.value(),
            },
        ))
    }
}

fn prepare_conductance(
    initial: &PyrolysisSimulation<'_>,
    w: &mut ImplicitWorkspace,
    checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
) -> Result<(), ImplicitError> {
    let n = w.ids.len();
    w.conductance.fill(0.0);
    for (i, &id) in w.ids.iter().enumerate() {
        checkpoint()?;
        let cell = &initial.cells()[&id];
        // Same coverage owner as heat_rates_with_calorics; no independently
        // reconstructed exposed area or invented endpoint coefficient.
        let ambient = initial.thermal.exposed_area_m2(id)? * cell.boundary.htc_w_m2_k;
        w.conductance[i * n + i] =
            finite(ambient + cell.boundary.reservoir.map_or(0.0, |r| r.conductance_w_k))?;
    }
    for link in initial.thermal.links().values() {
        checkpoint()?;
        let i = w
            .ids
            .binary_search(&link.a)
            .map_err(|_| DepositionError::Arithmetic)?;
        let j = w
            .ids
            .binary_search(&link.b)
            .map_err(|_| DepositionError::Arithmetic)?;
        let g = link.conductance_w_k;
        for k in [i * n + i, j * n + j] {
            w.conductance[k] = finite(w.conductance[k] + g)?;
        }
        for k in [i * n + j, j * n + i] {
            w.conductance[k] = finite(w.conductance[k] - g)?;
        }
    }
    Ok(())
}

fn initialize(
    initial: &PyrolysisSimulation<'_>,
    dt: f64,
    w: &mut ImplicitWorkspace,
    checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
) -> Result<(), ImplicitError> {
    // An isothermal reaction-energy guess lies at the accepted temperature,
    // even when zero absorbed energy would cool below the material domain.
    // This only seeds Newton; the original solve_cell determines every trial.
    for (i, &id) in w.ids.iter().enumerate() {
        checkpoint()?;
        let inventory = initial.inventories[&id];
        let model = initial.materials[&initial.cells()[&id].material_id];
        let h = model.enthalpies(inventory.temperature_k)?;
        let consumed = model
            .kinetics
            .isothermal_step(inventory.virgin_mass_kg, inventory.temperature_k, dt)?
            .consumed_mass_kg;
        w.energy[i] = finite(
            consumed
                * (model.char_mass_yield * (h[1] - h[0])
                    + (1.0 - model.char_mass_yield) * (h[2] - h[0])),
        )?;
    }
    Ok(())
}

fn local_solve(
    initial: &PyrolysisSimulation<'_>,
    id: CellId,
    dt: f64,
    input: f64,
    controls: PyrolysisControls,
    work: &mut Work,
    checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
) -> Result<CellStep, DepositionError> {
    work.cells += 1;
    solve_cell(
        initial.inventories[&id],
        initial.materials[&initial.cells()[&id].material_id],
        dt,
        input,
        PyrolysisControls {
            cell_energy_tolerance_j: 0.25 * controls.cell_energy_tolerance_j,
            ..controls
        },
        &mut || {
            checkpoint()?;
            // Original solve_cell polls immediately before each scalar energy
            // evaluation, including failed domain/bracket trials.
            work.temperature_evaluations += 1;
            Ok(())
        },
    )
}

fn evaluate(
    initial: &PyrolysisSimulation<'_>,
    candidate: &mut PyrolysisSimulation<'_>,
    dt: f64,
    controls: PyrolysisControls,
    ids: &[CellId],
    energy: &[f64],
    residual: &mut [f64],
    steps: &mut Vec<CellStep>,
    work: &mut Work,
    checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
) -> Result<(f64, f64), ImplicitError> {
    steps.clear();
    for (i, &id) in ids.iter().enumerate() {
        let step = local_solve(initial, id, dt, energy[i], controls, work, checkpoint)?;
        let cell = candidate
            .thermal
            .cells
            .get_mut(&id)
            .expect("admitted cell axis");
        cell.mass_kg = finite(step.inventory.virgin_mass_kg + step.inventory.char_mass_kg)?;
        cell.specific_h_j_kg = step.specific_h_j_kg;
        candidate.inventories.insert(id, step.inventory);
        steps.push(step);
    }
    let rates = candidate.thermal.heat_rates_with_calorics(
        f64::MAX,
        |id, cell| {
            let inventory = candidate.inventories[&id];
            Ok((
                inventory.temperature_k,
                candidate.materials[&cell.material_id].mixture_slope(inventory)?,
            ))
        },
        checkpoint,
    )?;
    work.heat += 1;
    let mut maximum = 0.0_f64;
    for (i, &id) in ids.iter().enumerate() {
        residual[i] = finite(energy[i] - finite(dt * rates.watts[&id])?)?;
        maximum = maximum.max(residual[i].abs());
    }
    Ok((maximum, rates.external_w))
}

fn derivative(
    initial: &PyrolysisSimulation<'_>,
    id: CellId,
    dt: f64,
    input: f64,
    temperature: f64,
    controls: PyrolysisControls,
    solve: ImplicitSolve,
    work: &mut Work,
    checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
) -> Result<f64, ImplicitError> {
    let model = initial.materials[&initial.cells()[&id].material_id];
    let middle = 0.5 * model.temperature_range_k[0] + 0.5 * model.temperature_range_k[1];
    let direction = if temperature > middle { -1.0 } else { 1.0 };
    let mut delta = solve
        .energy_difference_j
        .max(256.0 * controls.cell_energy_tolerance_j)
        .max(1e-7 * input.abs());
    for _ in 0..24 {
        for sign in [direction, -direction] {
            checkpoint()?;
            let value = finite(input + sign * delta)?;
            let actual = value - input;
            if actual == 0.0 {
                continue;
            }
            match local_solve(initial, id, dt, value, controls, work, checkpoint) {
                Ok(probe) => {
                    let slope = finite((probe.inventory.temperature_k - temperature) / actual)?;
                    if positive(slope) {
                        return Ok(slope);
                    }
                }
                Err(DepositionError::MaterialDomain | DepositionError::EnergyBalance(_)) => {}
                Err(error) => return Err(error.into()),
            }
        }
        delta *= 0.5;
    }
    Err(DepositionError::Arithmetic.into())
}

fn nonlinear_step(
    initial: &PyrolysisSimulation<'_>,
    candidate: &mut PyrolysisSimulation<'_>,
    dt: f64,
    controls: PyrolysisControls,
    solve: ImplicitSolve,
    w: &mut ImplicitWorkspace,
    work: &mut Work,
    checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
) -> Result<f64, ImplicitError> {
    let n = w.ids.len();
    let threshold = 0.25 * controls.cell_energy_tolerance_j;
    for _ in 0..solve.maximum_iterations {
        checkpoint()?;
        let (norm, external_w) = evaluate(
            initial,
            candidate,
            dt,
            controls,
            &w.ids,
            &w.energy,
            &mut w.residual,
            &mut w.steps,
            work,
            checkpoint,
        )?;
        work.iterations += 1;
        if norm <= threshold {
            return Ok(external_w);
        }
        for (j, &id) in w.ids.iter().enumerate() {
            w.slope[j] = derivative(
                initial,
                id,
                dt,
                w.energy[j],
                w.steps[j].inventory.temperature_k,
                controls,
                solve,
                work,
                checkpoint,
            )?;
            w.rhs[j] = -w.residual[j];
        }
        for i in 0..n {
            checkpoint()?;
            for j in 0..n {
                w.matrix[i * n + j] = finite(
                    dt * w.conductance[i * n + j] * w.slope[j] + if i == j { 1.0 } else { 0.0 },
                )?;
            }
        }
        w.lu.solve_into(&w.matrix, &w.rhs, &mut w.delta)?;
        work.linear += 1;
        let mut scale = 1.0;
        let mut accepted = false;
        for _ in 0..solve.maximum_backtracks {
            checkpoint()?;
            for i in 0..n {
                w.probe[i] = finite(w.energy[i] + scale * w.delta[i])?;
            }
            work.line_probes += 1;
            match evaluate(
                initial,
                candidate,
                dt,
                controls,
                &w.ids,
                &w.probe,
                &mut w.residual,
                &mut w.steps,
                work,
                checkpoint,
            ) {
                Ok((probe_norm, _))
                    if probe_norm <= threshold || probe_norm < norm * (1.0 - 1e-4 * scale) =>
                {
                    w.energy.copy_from_slice(&w.probe);
                    accepted = true;
                    break;
                }
                Ok(_)
                | Err(ImplicitError::Thermal(
                    DepositionError::MaterialDomain | DepositionError::EnergyBalance(_),
                )) => {}
                Err(error) => return Err(error),
            }
            scale *= 0.5;
        }
        if !accepted {
            return Err(ImplicitError::Iterations);
        }
    }
    Err(ImplicitError::Iterations)
}

fn accept_interval(
    initial: &PyrolysisSimulation<'_>,
    candidate: &PyrolysisSimulation<'_>,
    dt: f64,
    external_w: f64,
    controls: PyrolysisControls,
    w: &ImplicitWorkspace,
    audit: &mut PyrolysisAudit,
    account: &mut Account,
    checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
) -> Result<(), ImplicitError> {
    let mut graph_input = Sum::default();
    for (i, &id) in w.ids.iter().enumerate() {
        checkpoint()?;
        let old = &initial.cells()[&id];
        let new = &candidate.cells()[&id];
        let step = &w.steps[i];
        let relative_storage = finite(new.mass_kg * (new.specific_h_j_kg - old.specific_h_j_kg))?;
        let actual_mass_loss = old.mass_kg - new.mass_kg;
        // Same stable cell-relative physical accounting as the original
        // explicit-conduction method, including numerical mass reconciliation.
        let local = finite(
            relative_storage
                + step.release.mass_kg * (step.volatile_h_j_kg - old.specific_h_j_kg)
                + (step.release.mass_kg - actual_mass_loss) * old.specific_h_j_kg
                - w.energy[i]
                + w.residual[i],
        )?;
        if local.abs() > controls.cell_energy_tolerance_j {
            return Err(DepositionError::EnergyBalance(local).into());
        }
        if local.abs() + step.kinetic_energy_roundoff_j > controls.cell_energy_tolerance_j {
            return Err(DepositionError::Budget(
                "implicit pyrolysis local reaction-energy roundoff",
            )
            .into());
        }
        account.energy.add(local)?;
        account.absolute_energy.add(local.abs())?;
        account
            .mass
            .add(finite(step.release.mass_kg - actual_mass_loss)?)?;
        account.stored.add(finite(
            relative_storage - actual_mass_loss * old.specific_h_j_kg,
        )?)?;
        account.kinetic_mass.add(step.kinetic_mass_roundoff_kg)?;
        account.kinetic_energy.add(step.kinetic_energy_roundoff_j)?;
        graph_input.add(finite(w.energy[i] - w.residual[i])?)?;
        audit.maximum_cell_energy_residual_j =
            audit.maximum_cell_energy_residual_j.max(local.abs());
        let release = audit
            .releases
            .get_mut(&id)
            .expect("admitted source cell axis");
        release.mass_kg = finite(release.mass_kg + step.release.mass_kg)?;
        release.enthalpy_j = finite(release.enthalpy_j + step.release.enthalpy_j)?;
        release.consumed_virgin_mass_kg =
            finite(release.consumed_virgin_mass_kg + step.release.consumed_virgin_mass_kg)?;
    }
    let external_j = finite(dt * external_w)?;
    let graph_roundoff = finite(graph_input.value() - external_j)?;
    account.external.add(external_j)?;
    account.energy.add(graph_roundoff)?;
    account.absolute_energy.add(graph_roundoff.abs())?;
    account.heat_roundoff.add(graph_roundoff.abs())?;
    Ok(())
}

#[cfg(test)]
#[path = "implicit/tests.rs"]
mod tests;
