//! Simultaneous conservative backward Euler for a declared reaction mechanism.
//!
//! Solve every extent equation xi_r - dt V rate_r(n0 + S xi, U0, V) = 0
//! with fs-solver's Newton--Krylov/globalization owner. A trial reconstructs
//! ONE original gas state; all channels then query its actual T, pressure and
//! composition. Formation energy remains in the original total-U chart.
//!
//! The finite-difference Jacobian is an approximation used only to find the
//! root. Physical acceptance independently checks every molar residual, exact
//! stoichiometric inventory, source support, mass, U and entropy. Neither a
//! small Newton step nor a converged Krylov iteration accepts chemistry.
//! This is first order, with no uniqueness or unconditional convergence claim.
//! References: Cantera's ideal-gas control-volume equations and SUNDIALS CVODE
//! mathematics, nonlinear equation F = y - h f(y) - y0. Neither is a runtime
//! dependency, and this implementation does not claim CVODE's order/adaptivity.

use super::super::*;
use fs_solver::{NewtonKrylovConfig, NewtonKrylovState, NewtonReport, NonlinearProblem};
use fs_thermochem::ReactionId;
use std::cell::{Cell, RefCell};
use std::collections::BTreeMap;

/// Shared nonlinear policy, bound AFTER combining all source-bound channels.
#[derive(Clone, Debug, PartialEq)]
pub struct ImplicitMechanismSolve {
    /// One finite signed [lower, upper] extent range [mol] per reaction ID.
    /// All ranges contain zero. Unknown/missing IDs refuse; no order guessing.
    pub extent_brackets_mol: BTreeMap<ReactionId, [f64; 2]>,
    /// Independent maximum absolute endpoint-equation residual per channel [mol].
    pub extent_tolerance_mol: f64,
    /// Numerical interior seed fraction in (0, 1/4]. Only reachable initially
    /// absent species are seeded, by conservative reaction extents. This is an
    /// initial Newton guess, NOT an accepted burn or a species floor.
    pub interior_fraction: f64,
    /// Relative centered/one-sided Jacobian stencil on normalized extents.
    /// Must be in [32 epsilon, 0.01]; source support is required at each query.
    pub jacobian_relative_step: f64,
    /// Native Newton outer attempts, 1..=128.
    pub maximum_iterations: usize,
    /// All residual evaluations, including stencils and rejected trials,
    /// 1..=65536. Each original caloric inversion retains its separate cap.
    pub maximum_evaluations: usize,
    /// Native solver policy. Each extent and its residual is divided by that
    /// channel's largest absolute bracket endpoint. Final physical acceptance
    /// uses extent_tolerance_mol even if these tolerances stop Newton earlier.
    /// Restart length and cycles are each capped at 32 for this bounded lane.
    pub newton: NewtonKrylovConfig,
}

impl ImplicitMechanismSolve {
    fn validate(&self, law: &GasReaction) -> Result<(), GasReactionError> {
        if self.extent_brackets_mol.len() != law.channels().count()
            || !positive(self.extent_tolerance_mol)
            || !positive(self.interior_fraction)
            || self.interior_fraction > 0.25
            || !self.jacobian_relative_step.is_finite()
            || !(32.0 * f64::EPSILON..=0.01).contains(&self.jacobian_relative_step)
            || !(1..=128).contains(&self.maximum_iterations)
            || !(1..=65536).contains(&self.maximum_evaluations)
            || !(1..=32).contains(&self.newton.linear_restart)
            || !(1..=32).contains(&self.newton.max_linear_cycles)
        {
            return Err(invalid(
                "invalid coupled reaction bounds, tolerances or work budget",
            ));
        }
        self.newton
            .validate()
            .map_err(|_| invalid("invalid native reaction Newton policy"))?;
        for (kinetics, _) in law.channels() {
            let id = kinetics.reverse_rate_closure().equilibrium().reaction();
            let &[lo, hi] = self.extent_brackets_mol.get(id).ok_or_else(|| {
                invalid("coupled extent bounds must name every original reaction")
            })?;
            if !lo.is_finite() || !hi.is_finite() || lo >= hi || lo > 0.0 || hi < 0.0 {
                return Err(invalid("each coupled reaction bracket must contain zero"));
            }
        }
        Ok(())
    }
}

/// A distinct solved equation for each original stoichiometric column.
#[derive(Clone, Debug, PartialEq)]
pub struct ImplicitChannelResidual {
    pub reaction_id: ReactionId,
    pub extent_bracket_mol: [f64; 2],
    pub residual_mol: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ImplicitMechanismReport {
    pub duration_s: f64,
    /// Canonical channel order, including channels on invariant zero faces.
    pub equations: Vec<ImplicitChannelResidual>,
    pub maximum_residual_mol: f64,
    pub evaluations: usize,
    /// Original caloric inversion work at the accepted endpoint; trial work is
    /// separately bounded by evaluations and each GasEnergySolve policy.
    pub endpoint_energy_iterations: usize,
    /// None only when the original zero-extent state already closes all equations.
    pub nonlinear: Option<NewtonReport>,
}

impl GasReaction {
    /// Solve all reaction extents together at one conservative endpoint. This
    /// selector reaches existing parcel, piston, transport and thermal/circuit
    /// consumers through the same original parcel_candidate call. Optional
    /// dyadic subcycling refines this method's actual trajectories.
    pub fn with_implicit_mechanism(
        mut self,
        policy: ImplicitMechanismSolve,
    ) -> Result<Self, GasReactionError> {
        policy.validate(&self)?;
        if self.implicit.is_some() || self.implicit_mechanism.is_some() {
            return Err(invalid("one reaction integration method may be bound"));
        }
        self.implicit_mechanism = Some(policy);
        Ok(self)
    }

    #[must_use]
    pub fn implicit_mechanism(&self) -> Option<&ImplicitMechanismSolve> {
        self.implicit_mechanism.as_ref()
    }
}

// Negative trial inventories/range violations ask the native globalization to
// shorten its step. Missing sources, arithmetic, work exhaustion and cancellation
// are fatal and latched; line search cannot turn them into a successful result.
enum TrialError {
    OutsideDomain,
    Fatal(GasReactionError),
}
impl From<GasReactionError> for TrialError {
    fn from(error: GasReactionError) -> Self {
        Self::Fatal(error)
    }
}

struct Inventory {
    changes: Vec<f64>,
    amounts: Vec<f64>,
    consumed: f64,
}
struct Evaluation {
    gas: ClosedGasParcel,
    energy_iterations: usize,
    progress: mechanism::Progress,
    residual_mol: Vec<f64>,
}
impl Evaluation {
    fn maximum_residual(&self) -> f64 {
        self.residual_mol
            .iter()
            .map(|v| v.abs())
            .fold(0.0, f64::max)
    }
}
#[derive(Clone)]
struct Jacobian {
    point: Vec<f64>,
    entries: Vec<f64>,
}

struct Problem<'a, C> {
    law: &'a GasReaction,
    gas: &'a ClosedGasParcel,
    initial: Vec<f64>,
    columns: Vec<Vec<f64>>,
    // Structurally unreachable channels stay at zero. Their original rates
    // and residuals are STILL evaluated and admitted at the returned endpoint.
    active: Vec<usize>,
    scales: Vec<f64>,
    bounds: Vec<[f64; 2]>,
    duration_s: f64,
    heat_j: f64,
    controls: GasReactionSolve,
    policy: &'a ImplicitMechanismSolve,
    cancelled: RefCell<&'a mut C>,
    failure: RefCell<Option<GasReactionError>>,
    evaluations: Cell<usize>,
    jacobian: RefCell<Option<Jacobian>>,
}

impl<C: FnMut() -> bool> Problem<'_, C> {
    fn stop(&self) -> bool {
        if self.failure.borrow().is_some() {
            return true;
        }
        if (self.cancelled.borrow_mut())() {
            *self.failure.borrow_mut() = Some(GasReactionError::Cancelled);
            true
        } else {
            false
        }
    }
    fn poll(&self) -> Result<(), GasReactionError> {
        if self.stop() {
            Err(self
                .failure
                .borrow()
                .clone()
                .expect("stopped with a retained reason"))
        } else {
            Ok(())
        }
    }
    fn latch(&self, error: GasReactionError) {
        let mut first = self.failure.borrow_mut();
        if first.is_none() {
            *first = Some(error);
        }
    }
    fn extents(&self, x: &[f64]) -> Result<Vec<f64>, TrialError> {
        self.poll()?;
        if x.len() != self.active.len() {
            return Err(invalid("coupled extent dimension changed").into());
        }
        let mut extents = vec![0.0; self.columns.len()];
        for (&j, &value) in self.active.iter().zip(x) {
            let extent = value * self.scales[j];
            if !extent.is_finite() {
                return Err(invalid("coupled extent is unrepresentable").into());
            }
            if extent < self.bounds[j][0] || extent > self.bounds[j][1] {
                return Err(TrialError::OutsideDomain);
            }
            extents[j] = extent;
        }
        Ok(extents)
    }
    fn inventory(&self, extents: &[f64]) -> Result<Inventory, TrialError> {
        let mut changes = vec![mechanism::Sum::default(); self.initial.len()];
        let mut production = changes.clone();
        let mut withdrawal = changes.clone();
        for (column, &extent) in self.columns.iter().zip(extents) {
            self.poll()?;
            for (i, &coefficient) in column.iter().enumerate() {
                let change = coefficient * extent;
                changes[i].add(change)?;
                if change > 0.0 {
                    production[i].add(change)?;
                }
                if change < 0.0 {
                    withdrawal[i].add(-change)?;
                }
            }
        }
        let mut amounts = Vec::with_capacity(self.initial.len());
        let mut consumed = 0.0_f64;
        for (i, &amount) in self.initial.iter().enumerate() {
            self.poll()?;
            let next = amount + changes[i].value;
            if !next.is_finite() || (changes[i].value != 0.0 && next == amount) {
                return Err(
                    invalid("coupled reaction lost a represented species increment").into(),
                );
            }
            if next < 0.0 {
                return Err(TrialError::OutsideDomain);
            }
            if withdrawal[i].value > 0.0 {
                let available = amount + production[i].value;
                let fraction = withdrawal[i].value / available;
                if !positive(available) || !fraction.is_finite() {
                    return Err(
                        invalid("coupled reaction donor accounting is unrepresentable").into(),
                    );
                }
                consumed = consumed.max(fraction);
            }
            amounts.push(next);
        }
        Ok(Inventory {
            changes: changes.into_iter().map(|s| s.value).collect(),
            amounts,
            consumed,
        })
    }
    fn evaluate(&self, x: &[f64]) -> Result<Evaluation, TrialError> {
        self.poll()?;
        let spent = self.evaluations.get();
        if spent == self.policy.maximum_evaluations {
            return Err(GasReactionError::ImplicitEvaluationLimit {
                evaluations: spent,
                maximum: self.policy.maximum_evaluations,
            }
            .into());
        }
        self.evaluations.set(spent + 1);
        let extents = self.extents(x)?;
        let inventory = self.inventory(&extents)?;
        let (gas, energy_iterations) = self.law.reconstruct_inventory_with_heat(
            self.gas,
            &inventory.changes,
            &inventory.amounts,
            self.heat_j,
            self.controls,
            &mut || self.stop(),
        )?;
        let actual = self.law.inventory(&gas)?;
        let mut reports = Vec::with_capacity(extents.len());
        let mut residual_mol = Vec::with_capacity(extents.len());
        for ((kinetics, source), &extent) in self.law.channels().zip(&extents) {
            self.poll()?;
            let (rate, rate_material) =
                self.law
                    .evaluate_channel(kinetics, source, &gas, &actual, &mut || self.stop())?;
            let value = rate
                .net_progress_rate()
                .ok_or_else(|| {
                    invalid("coupled implicit rate is log-only; no finite endpoint equation")
                })?
                .value();
            let predicted = (value * gas.volume_m3()) * self.duration_s;
            let residual = extent - predicted;
            if !predicted.is_finite() || (value != 0.0 && predicted == 0.0) || !residual.is_finite()
            {
                return Err(
                    invalid("coupled reaction endpoint equation is unrepresentable").into(),
                );
            }
            reports.push(mechanism::ReactionChannelStep {
                reaction_id: kinetics
                    .reverse_rate_closure()
                    .equilibrium()
                    .reaction()
                    .clone(),
                rate,
                rate_material,
                extent_mol: extent,
            });
            residual_mol.push(residual);
        }
        let mut reports = reports.into_iter();
        Ok(Evaluation {
            gas,
            energy_iterations,
            progress: mechanism::Progress {
                primary: reports.next().expect("admitted nonempty mechanism"),
                additional: reports.collect(),
                changes: inventory.changes,
                expected: inventory.amounts,
                consumed: inventory.consumed,
            },
            residual_mol,
        })
    }
    fn normalized_residual(&self, value: &Evaluation) -> Result<Vec<f64>, GasReactionError> {
        let result: Vec<_> = self
            .active
            .iter()
            .map(|&j| value.residual_mol[j] / self.scales[j])
            .collect();
        if result.iter().any(|x| !x.is_finite()) {
            return Err(invalid(
                "normalized coupled reaction residual is unrepresentable",
            ));
        }
        Ok(result)
    }
    // Numerical interior construction only. Reachability prevents a disconnected
    // absent species from appearing merely to help Newton.
    fn interior_seed(&self) -> Result<Vec<f64>, GasReactionError> {
        let mut extents = vec![0.0; self.columns.len()];
        let mut amounts = self.initial.clone();
        for _ in 0..self.active.len() {
            let mut changed = false;
            for &j in &self.active {
                self.poll()?;
                let column = &self.columns[j];
                if column
                    .iter()
                    .zip(&amounts)
                    .all(|(&nu, &n)| nu == 0.0 || n > 0.0)
                {
                    continue;
                }
                let forward = column
                    .iter()
                    .zip(&amounts)
                    .all(|(&nu, &n)| nu >= 0.0 || n > 0.0);
                let reverse = column
                    .iter()
                    .zip(&amounts)
                    .all(|(&nu, &n)| nu <= 0.0 || n > 0.0);
                let sign = if forward && self.bounds[j][1] > extents[j] {
                    1.0
                } else if reverse && self.bounds[j][0] < extents[j] {
                    -1.0
                } else {
                    continue;
                };
                let stock = column
                    .iter()
                    .zip(&amounts)
                    .filter_map(|(&nu, &n)| {
                        if sign * nu < 0.0 {
                            Some(n / nu.abs())
                        } else {
                            None
                        }
                    })
                    .fold(f64::INFINITY, f64::min);
                let room = if sign > 0.0 {
                    self.bounds[j][1] - extents[j]
                } else {
                    extents[j] - self.bounds[j][0]
                };
                let delta = (self.policy.interior_fraction * stock).min(0.5 * room);
                if !positive(delta) {
                    return Err(invalid("no represented conservative interior seed"));
                }
                extents[j] += sign * delta;
                amounts = self.inventory(&extents).map_err(trial_error)?.amounts;
                changed = true;
            }
            if self.active.iter().all(|&j| {
                self.columns[j]
                    .iter()
                    .zip(&amounts)
                    .all(|(&nu, &n)| nu == 0.0 || n > 0.0)
            }) {
                return Ok(self
                    .active
                    .iter()
                    .map(|&j| extents[j] / self.scales[j])
                    .collect());
            }
            if !changed {
                break;
            }
        }
        Err(invalid(
            "declared extent bounds cannot seed the reachable reaction interior",
        ))
    }
    fn step_room(&self, x: &[f64], amounts: &[f64], column_index: usize, sign: f64) -> f64 {
        let j = self.active[column_index];
        let mut room = if sign > 0.0 {
            self.bounds[j][1] / self.scales[j] - x[column_index]
        } else {
            x[column_index] - self.bounds[j][0] / self.scales[j]
        };
        for (&nu, &n) in self.columns[j].iter().zip(amounts) {
            let derivative = sign * nu * self.scales[j];
            if derivative < 0.0 {
                room = room.min(n / -derivative);
            }
        }
        room
    }
    fn matrix(&self, x: &[f64]) -> Result<Jacobian, TrialError> {
        self.poll()?;
        if let Some(cache) = self.jacobian.borrow().as_ref() {
            if cache.point == x {
                return Ok(cache.clone());
            }
        }
        let base = self.evaluate(x)?;
        let base_residual = self.normalized_residual(&base)?;
        let amounts = &base.progress.expected;
        let n = x.len();
        let mut entries = vec![0.0; n * n];
        for j in 0..n {
            self.poll()?;
            let requested = self.policy.jacobian_relative_step * x[j].abs().max(1.0);
            let plus = self.step_room(x, amounts, j, 1.0);
            let minus = self.step_room(x, amounts, j, -1.0);
            let centered = plus > 0.0 && minus > 0.0;
            let sign = if plus > 0.0 { 1.0 } else { -1.0 };
            let room = if centered {
                plus.min(minus)
            } else if sign > 0.0 {
                plus
            } else {
                minus
            };
            let h = requested.min(0.5 * room);
            if !positive(h) {
                return Err(invalid("reaction Jacobian has no feasible extent stencil").into());
            }
            let mut high = x.to_vec();
            high[j] += sign * h;
            if high[j] == x[j] {
                return Err(invalid("reaction Jacobian increment rounded away").into());
            }
            let high_value = self.evaluate(&high)?;
            let high_residual = self.normalized_residual(&high_value)?;
            let (low, low_residual) = if centered {
                let mut low = x.to_vec();
                low[j] -= sign * h;
                if low[j] == x[j] {
                    return Err(invalid("reaction Jacobian increment rounded away").into());
                }
                let value = self.evaluate(&low)?;
                let residual = self.normalized_residual(&value)?;
                (low[j], residual)
            } else {
                (x[j], base_residual.clone())
            };
            let denominator = high[j] - low;
            for i in 0..n {
                let derivative = (high_residual[i] - low_residual[i]) / denominator;
                if !derivative.is_finite() {
                    return Err(invalid("reaction Jacobian is unrepresentable").into());
                }
                entries[i * n + j] = derivative;
            }
        }
        let result = Jacobian {
            point: x.to_vec(),
            entries,
        };
        *self.jacobian.borrow_mut() = Some(result.clone());
        Ok(result)
    }
    fn callback_failure(&self, error: TrialError) {
        if let TrialError::Fatal(error) = error {
            self.latch(error);
        }
    }
}

impl<C: FnMut() -> bool> NonlinearProblem for Problem<'_, C> {
    fn dimension(&self) -> usize {
        self.active.len()
    }
    fn residual(&self, x: &[f64], residual: &mut [f64]) {
        residual.fill(f64::NAN);
        match self
            .evaluate(x)
            .and_then(|value| self.normalized_residual(&value).map_err(Into::into))
        {
            Ok(value) => residual.copy_from_slice(&value),
            Err(error) => self.callback_failure(error),
        }
    }
    fn jacobian_apply(&self, x: &[f64], direction: &[f64], output: &mut [f64]) {
        output.fill(f64::NAN);
        match self.matrix(x) {
            Ok(matrix) => {
                let n = self.dimension();
                for (i, target) in output.iter_mut().enumerate() {
                    if let Err(error) = self.poll() {
                        self.latch(error);
                        return;
                    }
                    *target = matrix.entries[i * n..(i + 1) * n]
                        .iter()
                        .zip(direction)
                        .map(|(a, b)| a * b)
                        .sum();
                }
            }
            Err(error) => self.callback_failure(error),
        }
    }
    fn preconditioner_apply(
        &self,
        x: &[f64],
        _: usize,
        _: usize,
        residual: &[f64],
        output: &mut [f64],
    ) {
        output.fill(f64::NAN);
        match self.matrix(x) {
            Ok(matrix) => {
                let n = self.dimension();
                for i in 0..n {
                    if let Err(error) = self.poll() {
                        self.latch(error);
                        return;
                    }
                    let diagonal = matrix.entries[i * n + i];
                    output[i] = if diagonal.abs() >= f64::EPSILON {
                        residual[i] / diagonal
                    } else {
                        residual[i]
                    };
                }
            }
            Err(error) => self.callback_failure(error),
        }
    }
}

fn trial_error(error: TrialError) -> GasReactionError {
    match error {
        TrialError::Fatal(error) => error,
        TrialError::OutsideDomain => {
            invalid("coupled reaction endpoint lies outside its physical extent domain")
        }
    }
}

pub(in crate::cavity_phs::piston::reactive) fn candidate(
    law: &GasReaction,
    gas: &ClosedGasParcel,
    duration_s: f64,
    controls: GasReactionSolve,
    policy: &ImplicitMechanismSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<GasReactionCandidate, GasReactionError> {
    candidate_with_heat(law, gas, duration_s, 0.0, controls, policy, cancelled)
}

pub(in crate::cavity_phs::piston::reactive) fn candidate_with_heat(
    law: &GasReaction,
    gas: &ClosedGasParcel,
    duration_s: f64,
    heat_j: f64,
    controls: GasReactionSolve,
    policy: &ImplicitMechanismSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<GasReactionCandidate, GasReactionError> {
    poll(cancelled)?;
    policy.validate(law)?;
    let initial = law.inventory(gas)?;
    let components = law.catalog.thermodynamics().components();
    let mut columns = Vec::with_capacity(law.channels().count());
    let mut bounds = Vec::with_capacity(law.channels().count());
    let mut scales = Vec::with_capacity(law.channels().count());
    for (kinetics, _) in law.channels() {
        poll(cancelled)?;
        let equilibrium = kinetics.reverse_rate_closure().equilibrium();
        let mut column = vec![0.0; initial.len()];
        for term in equilibrium.terms() {
            let i = components
                .binary_search_by(|s| s.species().cmp(term.species()))
                .map_err(|_| invalid("bound reaction species disappeared"))?;
            column[i] = term.coefficient() as f64;
        }
        let bracket = policy.extent_brackets_mol[equilibrium.reaction()];
        columns.push(column);
        bounds.push(bracket);
        scales.push(bracket[0].abs().max(bracket[1].abs()));
    }
    // Mass-action zero faces are invariant unless a supported forward/reverse
    // reaction can supply the absent species. This is graph reachability over
    // exact declared stoichiometry, not concentration clipping or a kinetic law.
    let mut reachable: Vec<_> = initial.iter().map(|n| *n > 0.0).collect();
    let mut active = vec![false; columns.len()];
    for _ in 0..columns.len() {
        let mut changed = false;
        for (j, column) in columns.iter().enumerate() {
            poll(cancelled)?;
            let forward = column
                .iter()
                .zip(&reachable)
                .all(|(&nu, &yes)| nu >= 0.0 || yes);
            let reverse = column
                .iter()
                .zip(&reachable)
                .all(|(&nu, &yes)| nu <= 0.0 || yes);
            if forward || reverse {
                if !active[j] {
                    active[j] = true;
                    changed = true;
                }
                for (&nu, yes) in column.iter().zip(&mut reachable) {
                    if nu != 0.0 && !*yes {
                        *yes = true;
                        changed = true;
                    }
                }
            }
        }
        if !changed {
            break;
        }
    }
    let problem = Problem {
        law,
        gas,
        initial,
        columns,
        active: active
            .iter()
            .enumerate()
            .filter_map(|(i, yes)| yes.then_some(i))
            .collect(),
        scales,
        bounds,
        duration_s,
        heat_j,
        controls,
        policy,
        cancelled: RefCell::new(cancelled),
        failure: RefCell::new(None),
        evaluations: Cell::new(0),
        jacobian: RefCell::new(None),
    };
    let original = problem
        .evaluate(&vec![0.0; problem.dimension()])
        .map_err(trial_error)?;
    let (accepted, nonlinear) = if original.maximum_residual() <= policy.extent_tolerance_mol {
        (original, None)
    } else {
        let seed = problem.interior_seed()?;
        let construction = NewtonKrylovState::new(&problem, seed, policy.newton);
        if let Some(error) = problem.failure.borrow().clone() {
            return Err(error);
        }
        let mut state =
            construction.map_err(|_| invalid("native reaction Newton construction refused"))?;
        let run = state.run_controlled(&problem, policy.maximum_iterations, &mut || problem.stop());
        if let Some(error) = problem.failure.borrow().clone() {
            return Err(error);
        }
        let endpoint = problem.evaluate(&state.x).map_err(trial_error)?;
        if !run.report.converged || endpoint.maximum_residual() > policy.extent_tolerance_mol {
            return Err(GasReactionError::ImplicitNotConverged {
                iterations: run.report.iterations,
                residual_mol: endpoint.maximum_residual(),
            });
        }
        (endpoint, Some(run.report))
    };
    // Gross withdrawal is limited against original stock PLUS simultaneous
    // production, which permits reaction intermediates while bounding turnover.
    if accepted.progress.consumed > controls.maximum_consumed_fraction {
        let mut production = vec![0.0; problem.initial.len()];
        let mut withdrawal = production.clone();
        let extents =
            std::iter::once(&accepted.progress.primary).chain(&accepted.progress.additional);
        for (column, report) in problem.columns.iter().zip(extents) {
            for (i, &nu) in column.iter().enumerate() {
                let change = nu * report.extent_mol;
                if change > 0.0 {
                    production[i] += change;
                } else {
                    withdrawal[i] -= change;
                }
            }
        }
        for (i, (&withdrawn, &produced)) in withdrawal.iter().zip(&production).enumerate() {
            if withdrawn > 0.0 {
                let fraction = withdrawn / (problem.initial[i] + produced);
                if fraction > controls.maximum_consumed_fraction {
                    return Err(GasReactionError::ConsumptionLimit {
                        species: components[i].species().clone(),
                        fraction,
                        maximum: controls.maximum_consumed_fraction,
                    });
                }
            }
        }
        return Err(invalid(
            "coupled reaction turnover exceeded the admitted limit",
        ));
    }
    let maximum_residual_mol = accepted.maximum_residual();
    let endpoint_energy_iterations = accepted.energy_iterations;
    let equations = law
        .channels()
        .enumerate()
        .map(|(j, (kinetics, _))| ImplicitChannelResidual {
            reaction_id: kinetics
                .reverse_rate_closure()
                .equilibrium()
                .reaction()
                .clone(),
            extent_bracket_mol: problem.bounds[j],
            residual_mol: accepted.residual_mol[j],
        })
        .collect();
    let mut result = law.admit_progress_with_heat(
        gas,
        &problem.initial,
        accepted.gas,
        accepted.progress,
        heat_j,
        controls,
        &mut || problem.stop(),
    )?;
    result.step.implicit_mechanism = Some(ImplicitMechanismReport {
        duration_s,
        equations,
        maximum_residual_mol,
        evaluations: problem.evaluations.get(),
        endpoint_energy_iterations,
        nonlinear,
    });
    problem.poll()?;
    Ok(result)
}
