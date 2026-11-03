//! Fixed-clock source-tensor dynamics with exact equal-time Maxwell memory.
//!
//! Every bit-identical relaxation time retains one complete displacement-memory
//! vector y: y' = v-y/tau. Each ORIGINAL material factor observes z=A*y, so its
//! spring force, storage and heat are unchanged. No weak direction, source time,
//! mass cross term or quadrature point is discarded. This realization admits
//! explicitly relaxed or instantaneous initial histories, not arbitrary imported
//! point memories. Nonzero thermal eigenstrain requires the existing full-state
//! owner and is refused here.
//!
//! The original fs-phs flow and positive source-square kernels are prepared once.
//! Warm calls use bounded contractions and source-domain checks, independently
//! integrate held or affine-force work and grouped physical heat, and publish atomically.
//! Material properties, geometry, force ports and the positive clock stay fixed.

use super::{Inertia, InitialMemory, TensorMaterialError};
use fs_exec::CancelGate;
use fs_material::visco_tensor::ProjectedTensorProny;
use fs_phs::linear_exact::prepared::PreparedQuadraticFlow;
use fs_phs::linear_exact::{DenseExponentialBudget, PreparedSourceLoss};
use fs_phs::{PhsError, PreparedStepError};

mod affine;
mod prepare;
pub use affine::PreparedTensorAffineStep;

/// Explicit whole-owner resource and physical numerical allowances. No defaults.
#[derive(Clone, Copy, Debug)]
pub struct PreparedTensorControls {
    pub maximum_coordinates: usize,
    pub maximum_material_points: usize,
    /// Complete propagation coordinates, including force records for affine loading.
    /// Accepted physical state remains displacement, momentum and material memory.
    pub maximum_state_coordinates: usize,
    pub maximum_heat_groups: usize,
    /// Includes retained source data, ALL channel caches and simultaneous scratch.
    pub maximum_preparation_bytes: usize,
    /// Includes ALL native flow/loss preparations and source-row compression.
    pub maximum_preparation_multiply_terms: u64,
    /// Complete warm contractions, source strain and original storage checks.
    pub maximum_step_multiply_terms: u64,
    /// Original mass-inverse and source-factor reconstruction allowance.
    pub maximum_relative_operator_error: f64,
    /// Absolute interval AND cumulative physical work-balance allowance [J].
    pub maximum_work_error_j: f64,
    pub maximum_steps: u64,
}

#[derive(Debug)]
pub enum PreparedTensorError {
    Invalid(&'static str),
    Budget { required: u128, maximum: u128 },
    Operator { defect: f64, allowance: f64 },
    Energy { defect_j: f64, allowance_j: f64 },
    Original(TensorMaterialError),
    Propagation(PhsError),
    Cancelled,
}
impl core::fmt::Display for PreparedTensorError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "prepared tensor dynamics refused: {self:?}")
    }
}
impl core::error::Error for PreparedTensorError {}
impl From<TensorMaterialError> for PreparedTensorError {
    fn from(error: TensorMaterialError) -> Self {
        match error {
            TensorMaterialError::Cancelled => Self::Cancelled,
            other => Self::Original(other),
        }
    }
}
impl From<PhsError> for PreparedTensorError {
    fn from(error: PhsError) -> Self {
        Self::Propagation(error)
    }
}
impl From<PreparedStepError> for PreparedTensorError {
    fn from(error: PreparedStepError) -> Self {
        match error {
            PreparedStepError::Cancelled => Self::Cancelled,
            PreparedStepError::Solver(error) => Self::Propagation(error),
        }
    }
}
type Result<T> = core::result::Result<T, PreparedTensorError>;
fn invalid(message: &'static str) -> PreparedTensorError {
    PreparedTensorError::Invalid(message)
}
fn poll(gate: &CancelGate) -> Result<()> {
    if gate.is_requested() {
        Err(PreparedTensorError::Cancelled)
    } else {
        Ok(())
    }
}
fn zeros(n: usize) -> Result<Vec<f64>> {
    let mut result = Vec::new();
    result
        .try_reserve_exact(n)
        .map_err(|_| invalid("allocation refused"))?;
    result.resize(n, 0.0);
    Ok(result)
}
fn finite(values: &[f64]) -> bool {
    values.iter().all(|x| x.is_finite())
}
fn add(sum: &mut f64, correction: &mut f64, value: f64) {
    let increment = value - *correction;
    let next = *sum + increment;
    *correction = (next - *sum) - increment;
    *sum = next;
}
fn bound(required: u128, maximum: u128) -> Result<()> {
    if required > maximum {
        Err(PreparedTensorError::Budget { required, maximum })
    } else {
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, Default)]
struct EnergyParts {
    equilibrium: f64,
    kinetic: f64,
    memory: f64,
}
impl EnergyParts {
    fn total(self) -> f64 {
        self.equilibrium + self.kinetic + self.memory
    }
}

/// Original physical storage, applied work and nonnegative material heat [J].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PreparedTensorEnergy {
    pub equilibrium_before_j: f64,
    pub equilibrium_after_j: f64,
    pub kinetic_before_j: f64,
    pub kinetic_after_j: f64,
    pub memory_before_j: f64,
    pub memory_after_j: f64,
    pub supplied_work_j: f64,
    pub dissipated_j: f64,
    pub cumulative_work_j: f64,
    pub cumulative_dissipated_j: f64,
    /// W - D - Delta H, using a stable original-coordinate energy increment.
    pub balance_error_j: f64,
    pub cumulative_balance_error_j: f64,
}

#[derive(Clone, Debug)]
pub struct PreparedTensorStep {
    pub step: u64,
    pub time_s: f64,
    /// Delta u / h, work-conjugate to the supplied held force [m/s].
    pub mean_velocity_mps: Vec<f64>,
    /// Original Kelvin plus Maxwell heat, in the declared channel order [J].
    pub heat_by_group_j: Vec<f64>,
    pub energy: PreparedTensorEnergy,
}

/// Immutable original source, consistent inertia and fixed-clock exact caches.
/// A heat group may collect several points (for example one physical element).
/// It changes reporting/cost only; every original positive source row survives.
pub struct PreparedTensorBody {
    inertia: Inertia,
    points: Vec<ProjectedTensorProny>,
    heat_groups: Vec<usize>,
    times: Vec<f64>,
    point_times: Vec<Vec<usize>>,
    equilibrium: Vec<f64>,
    flow: PreparedQuadraticFlow,
    heat: Vec<PreparedSourceLoss>,
    heat_scales: Vec<f64>,
    affine: Option<affine::AffinePrepared>,
    states: usize,
    dt_s: f64,
    controls: PreparedTensorControls,
}
impl PreparedTensorBody {
    pub fn coordinates(&self) -> usize {
        self.inertia.dimension()
    }
    pub fn state_dimension(&self) -> usize {
        self.states
    }
    /// Prepared physical plus source-record coordinates. Affine loading adds
    /// exactly one force coordinate per physical displacement; these carry no
    /// mechanical storage and never enter accepted constitutive history.
    pub fn propagation_state_dimension(&self) -> usize {
        self.states
            + if self.affine.is_some() {
                self.coordinates()
            } else {
                0
            }
    }
    /// Increasing exact positive source times; adjacent f64 values stay distinct.
    pub fn relaxation_times_s(&self) -> &[f64] {
        &self.times
    }
    pub fn material_points(&self) -> &[ProjectedTensorProny] {
        &self.points
    }
    pub fn heat_groups(&self) -> &[usize] {
        &self.heat_groups
    }
    pub fn heat_group_count(&self) -> usize {
        self.heat.len()
    }
    pub fn duration_s(&self) -> f64 {
        self.dt_s
    }
    pub fn mass_matrix_kg(&self) -> &[f64] {
        match &self.inertia {
            Inertia::Consistent { mass, .. } => mass,
            Inertia::Diagonal(mass) => mass, // This constructor always retains consistent M.
        }
    }

    /// Explicit initial prehistory. Counters start here and retain all initial
    /// equilibrium, kinetic and Maxwell spring storage. Kelvin prehistory is
    /// not integrated or assigned a fictitious finite jump work.
    pub fn history<'a>(
        &'a self,
        displacement_m: &[f64],
        velocity_mps: &[f64],
        memory: InitialMemory,
        gate: &CancelGate,
    ) -> Result<PreparedTensorHistory<'a>> {
        poll(gate)?;
        let n = self.coordinates();
        if displacement_m.len() != n
            || velocity_mps.len() != n
            || !finite(displacement_m)
            || !finite(velocity_mps)
        {
            return Err(invalid(
                "complete finite initial displacement and velocity required",
            ));
        }
        let mut state = zeros(self.states)?;
        state[..n].copy_from_slice(displacement_m);
        let momentum = self.inertia.momentum(velocity_mps);
        if !finite(&momentum) {
            return Err(invalid("initial momentum overflow"));
        }
        state[n..2 * n].copy_from_slice(&momentum);
        if memory == InitialMemory::Instantaneous {
            for y in state[2 * n..].chunks_exact_mut(n) {
                y.copy_from_slice(displacement_m);
            }
        }
        let velocity = self.inertia.velocity(&momentum);
        let scale = velocity_mps.iter().fold(0.0_f64, |s, x| s.max(x.abs()));
        let defect = velocity
            .iter()
            .zip(velocity_mps)
            .fold(0.0_f64, |s, (a, b)| s.max((a - b).abs()));
        if !finite(&velocity) || defect > self.controls.maximum_relative_operator_error * scale {
            return Err(PreparedTensorError::Operator {
                defect,
                allowance: self.controls.maximum_relative_operator_error * scale,
            });
        }
        self.admit_displacement(displacement_m, gate)?;
        let energy = self.energy(&state, &velocity, gate)?;
        let scratch = Scratch::new(self.states, n, self.affine.is_some())?;
        poll(gate)?;
        Ok(PreparedTensorHistory {
            body: self,
            state,
            velocity,
            scratch,
            steps: 0,
            energy,
            initial_energy_j: energy.total(),
            work_j: 0.0,
            work_correction_j: 0.0,
            heat_j: 0.0,
            heat_correction_j: 0.0,
        })
    }

    fn admit_displacement(&self, displacement: &[f64], gate: &CancelGate) -> Result<()> {
        for point in &self.points {
            poll(gate)?;
            point
                .strain(displacement)
                .map_err(TensorMaterialError::Material)?;
        }
        Ok(())
    }
    fn energy(&self, state: &[f64], velocity: &[f64], gate: &CancelGate) -> Result<EnergyParts> {
        let n = self.coordinates();
        let equilibrium = quadratic(&self.equilibrium, &state[..n], gate)?;
        let kinetic = quadratic(self.mass_matrix_kg(), velocity, gate)?;
        let (mut memory, mut correction) = (0.0, 0.0);
        for (point, times) in self.points.iter().zip(&self.point_times) {
            poll(gate)?;
            for (branch, &time) in point.branches().iter().zip(times) {
                let y = &state[(2 + time) * n..(3 + time) * n];
                for row in branch.rows() {
                    let z = row.iter().zip(y).fold(0.0, |s, (a, b)| a.mul_add(*b, s));
                    add(&mut memory, &mut correction, 0.5 * z * z);
                }
            }
        }
        let result = EnergyParts {
            equilibrium,
            kinetic,
            memory,
        };
        if !result.total().is_finite()
            || equilibrium < -self.controls.maximum_work_error_j
            || kinetic < 0.0
            || memory < 0.0
        {
            return Err(invalid(
                "original physical storage is invalid or unrepresentable",
            ));
        }
        Ok(result)
    }
    fn energy_change(
        &self,
        before: &[f64],
        velocity_before: &[f64],
        after: &[f64],
        velocity_after: &[f64],
        gate: &CancelGate,
    ) -> Result<f64> {
        let n = self.coordinates();
        let mut change = quadratic_change(&self.equilibrium, &before[..n], &after[..n], gate)?;
        let mut correction = 0.0;
        add(
            &mut change,
            &mut correction,
            quadratic_change(self.mass_matrix_kg(), velocity_before, velocity_after, gate)?,
        );
        for (point, times) in self.points.iter().zip(&self.point_times) {
            poll(gate)?;
            for (branch, &time) in point.branches().iter().zip(times) {
                let old = &before[(2 + time) * n..(3 + time) * n];
                let new = &after[(2 + time) * n..(3 + time) * n];
                for row in branch.rows() {
                    let mut mean = 0.0;
                    let mut increment = 0.0;
                    for j in 0..n {
                        mean = row[j].mul_add(0.5 * old[j] + 0.5 * new[j], mean);
                        increment = row[j].mul_add(new[j] - old[j], increment);
                    }
                    add(&mut change, &mut correction, mean * increment);
                }
            }
        }
        if !change.is_finite() {
            return Err(invalid("original storage change overflow"));
        }
        Ok(change)
    }
}

#[derive(Clone)]
struct Scratch {
    augmented: Vec<f64>,
    balanced: Vec<f64>,
    endpoint: Vec<f64>,
    heat_state: Vec<f64>,
    heat_coordinates: Vec<f64>,
    source_state: Vec<f64>,
    physical_candidate: Vec<f64>,
    source_increment: Vec<f64>,
    integrals: Vec<f64>,
}
impl Scratch {
    fn new(physical_states: usize, ports: usize, affine: bool) -> Result<Self> {
        let states = physical_states + if affine { ports } else { 0 };
        Ok(Self {
            augmented: zeros(states + ports + 1)?,
            balanced: zeros(states)?,
            endpoint: zeros(states)?,
            heat_state: zeros(states)?,
            heat_coordinates: zeros(2 * states)?,
            source_state: if affine { zeros(states)? } else { Vec::new() },
            physical_candidate: if affine {
                zeros(physical_states)?
            } else {
                Vec::new()
            },
            source_increment: if affine { zeros(ports)? } else { Vec::new() },
            integrals: if affine { zeros(1)? } else { Vec::new() },
        })
    }
}

/// Accepted motion borrows immutable caches. Clones fork only O(state) history
/// and scratch; they never duplicate dense prepared operators or material laws.
#[derive(Clone)]
pub struct PreparedTensorHistory<'a> {
    body: &'a PreparedTensorBody,
    state: Vec<f64>,
    velocity: Vec<f64>,
    scratch: Scratch,
    steps: u64,
    energy: EnergyParts,
    initial_energy_j: f64,
    work_j: f64,
    work_correction_j: f64,
    heat_j: f64,
    heat_correction_j: f64,
}
impl PreparedTensorHistory<'_> {
    pub fn body(&self) -> &PreparedTensorBody {
        self.body
    }
    /// Complete [u,p,y_tau...] state; y_tau is SPRING displacement memory.
    pub fn state(&self) -> &[f64] {
        &self.state
    }
    pub fn displacement_m(&self) -> &[f64] {
        &self.state[..self.body.coordinates()]
    }
    pub fn velocity_mps(&self) -> &[f64] {
        &self.velocity
    }
    /// One complete block per relaxation_times_s(); original spring z=A*y.
    pub fn memory_displacement_m(&self) -> &[f64] {
        &self.state[2 * self.body.coordinates()..]
    }
    pub fn time_s(&self) -> f64 {
        self.steps as f64 * self.body.dt_s
    }
    pub fn steps(&self) -> u64 {
        self.steps
    }
    pub fn stored_energy_j(&self) -> f64 {
        self.energy.total()
    }
    pub fn supplied_work_j(&self) -> f64 {
        self.work_j
    }
    pub fn dissipated_heat_j(&self) -> f64 {
        self.heat_j
    }
    pub fn balance_residual_j(&self) -> f64 {
        self.energy.total() - self.initial_energy_j + self.heat_j - self.work_j
    }

    /// Exact held-force interval of this fixed source. Finite force jumps are
    /// permitted and apply no impulse. Source strain limits are checked at the
    /// candidate endpoint; this is not a within-interval strain/band certificate.
    /// Any refusal leaves physical state, clock and every work/heat counter intact.
    pub fn advance(&mut self, force_n: &[f64], gate: &CancelGate) -> Result<PreparedTensorStep> {
        if self.body.affine.is_some() {
            return self
                .advance_affine(force_n, force_n, gate)
                .map(|step| step.motion);
        }
        poll(gate)?;
        let body = self.body;
        let n = body.coordinates();
        if force_n.len() != n || !finite(force_n) {
            return Err(invalid("complete finite held force required"));
        }
        let step = self
            .steps
            .checked_add(1)
            .filter(|s| *s <= body.controls.maximum_steps)
            .ok_or_else(|| invalid("accepted step budget exhausted"))?;
        let time_s = step as f64 * body.dt_s;
        if !time_s.is_finite() || time_s <= self.time_s() {
            return Err(invalid("representable later fixed clock required"));
        }
        body.flow.evaluate_into(
            &self.state,
            force_n,
            &mut self.scratch.augmented,
            &mut self.scratch.balanced,
            &mut self.scratch.endpoint,
            &mut [],
        )?;
        poll(gate)?;
        body.admit_displacement(&self.scratch.endpoint[..n], gate)?;
        let velocity = body.inertia.velocity(&self.scratch.endpoint[n..2 * n]);
        if !finite(&velocity) {
            return Err(invalid("candidate velocity overflow"));
        }
        let after = body.energy(&self.scratch.endpoint, &velocity, gate)?;
        let change = body.energy_change(
            &self.state,
            &self.velocity,
            &self.scratch.endpoint,
            &velocity,
            gate,
        )?;
        let mut mean_velocity_mps = zeros(n)?;
        let (mut supplied_work_j, mut correction) = (0.0, 0.0);
        for i in 0..n {
            let du = self.scratch.endpoint[i] - self.state[i];
            add(&mut supplied_work_j, &mut correction, force_n[i] * du);
            mean_velocity_mps[i] = du / body.dt_s;
        }
        for i in 0..body.states {
            self.scratch.heat_state[i] = body.heat_scales[i] * self.state[i];
        }
        let mut heat_by_group_j = zeros(body.heat.len())?;
        let (mut dissipated_j, mut heat_roundoff) = (0.0, 0.0);
        for (out, cache) in heat_by_group_j.iter_mut().zip(&body.heat) {
            poll(gate)?;
            *out = cache.evaluate_into_controlled(
                &self.scratch.heat_state,
                force_n,
                &mut self.scratch.heat_coordinates,
                || gate.is_requested(),
            )?;
            if !out.is_finite() || *out < 0.0 {
                return Err(invalid(
                    "original source heat must be finite and nonnegative",
                ));
            }
            add(&mut dissipated_j, &mut heat_roundoff, *out);
        }
        let (mut work_j, mut work_correction_j) = (self.work_j, self.work_correction_j);
        let (mut heat_j, mut heat_correction_j) = (self.heat_j, self.heat_correction_j);
        add(&mut work_j, &mut work_correction_j, supplied_work_j);
        add(&mut heat_j, &mut heat_correction_j, dissipated_j);
        if !finite(&[
            supplied_work_j,
            dissipated_j,
            work_j,
            work_correction_j,
            heat_j,
            heat_correction_j,
        ]) || !finite(&mean_velocity_mps)
            || heat_j < self.heat_j
        {
            return Err(invalid("physical work/heat ledger overflow"));
        }
        let energy = PreparedTensorEnergy {
            equilibrium_before_j: self.energy.equilibrium,
            equilibrium_after_j: after.equilibrium,
            kinetic_before_j: self.energy.kinetic,
            kinetic_after_j: after.kinetic,
            memory_before_j: self.energy.memory,
            memory_after_j: after.memory,
            supplied_work_j,
            dissipated_j,
            cumulative_work_j: work_j,
            cumulative_dissipated_j: heat_j,
            balance_error_j: supplied_work_j - dissipated_j - change,
            cumulative_balance_error_j: work_j - heat_j - (after.total() - self.initial_energy_j),
        };
        for defect_j in [energy.balance_error_j, energy.cumulative_balance_error_j] {
            if !defect_j.is_finite() || defect_j.abs() > body.controls.maximum_work_error_j {
                return Err(PreparedTensorError::Energy {
                    defect_j,
                    allowance_j: body.controls.maximum_work_error_j,
                });
            }
        }
        poll(gate)?;
        core::mem::swap(&mut self.state, &mut self.scratch.endpoint);
        self.velocity = velocity;
        self.energy = after;
        self.steps = step;
        self.work_j = work_j;
        self.work_correction_j = work_correction_j;
        self.heat_j = heat_j;
        self.heat_correction_j = heat_correction_j;
        Ok(PreparedTensorStep {
            step,
            time_s,
            mean_velocity_mps,
            heat_by_group_j,
            energy,
        })
    }
}

fn quadratic(matrix: &[f64], x: &[f64], gate: &CancelGate) -> Result<f64> {
    let n = x.len();
    let (mut total, mut correction) = (0.0, 0.0);
    for i in 0..n {
        poll(gate)?;
        let row = (0..n).fold(0.0, |s, j| matrix[i * n + j].mul_add(x[j], s));
        add(&mut total, &mut correction, 0.5 * x[i] * row);
    }
    if !total.is_finite() {
        return Err(invalid("original quadratic energy overflow"));
    }
    Ok(total)
}
fn quadratic_change(
    matrix: &[f64],
    before: &[f64],
    after: &[f64],
    gate: &CancelGate,
) -> Result<f64> {
    let n = before.len();
    let (mut total, mut correction) = (0.0, 0.0);
    for i in 0..n {
        poll(gate)?;
        let row = (0..n).fold(0.0, |s, j| {
            matrix[i * n + j].mul_add(after[j] - before[j], s)
        });
        add(
            &mut total,
            &mut correction,
            (0.5 * before[i] + 0.5 * after[i]) * row,
        );
    }
    if !total.is_finite() {
        return Err(invalid("original quadratic energy change overflow"));
    }
    Ok(total)
}
