//! Original-coordinate DKT plate dynamics with Maxwell bending memory.
//!
//! `M v_dot + K_eq u + sum(K_j q_j) = f`, `q_j_dot = v-q_j/tau_j`.
//! Each interval uses midpoint momentum and the exact constitutive solution
//! along its affine displacement path. Branch heat is integrated in physical
//! curvature Gram channels; it is not a work/energy remainder. The sparse
//! effective matrix is factored once for the selected fixed time step.

use crate::{
    NoPlateLoad, PlateDynamics, PlateDynamicsBudget, PlateDynamicsError, PlateDynamicsParameters,
    PlateModel,
};
use fs_material::plane_stress_viscosity::factor_passive_resultant;
use fs_material::visco::PronyRamp;
use fs_plate::viscosity::{PlateBendingViscosity, PlateViscosityBudget, ViscousPlateSection};
use fs_plate::{PlateError, PlateMesh, dkt_curvature_quadrature};
use fs_sparse::{Coo, Csr, DirectOrdering, LdltError, LdltFactor, LdltOptions, SymbolicLdlt};

/// Admission, bounded sparse solve, or cancellation refusal. Failed advances
/// leave all accepted displacement, velocity, material history and ledgers intact.
#[derive(Debug)]
pub enum PlateRelaxationError {
    Invalid(&'static str),
    Budget(&'static str),
    Plate(PlateError),
    Dynamics(PlateDynamicsError),
    Factor(LdltError),
    SolveResidual { relative: f64, allowed: f64 },
    Cancelled,
}
impl std::fmt::Display for PlateRelaxationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "plate relaxation failed: {self:?}")
    }
}
impl std::error::Error for PlateRelaxationError {}

type Error = PlateRelaxationError;
fn poll<C: FnMut() -> bool>(cancelled: &mut C) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(values: &[f64]) -> bool {
    values.iter().all(|x| x.is_finite())
}
fn row(matrix: &Csr, i: usize, values: &[f64]) -> f64 {
    let (columns, coefficients) = matrix.row(i);
    columns
        .iter()
        .zip(coefficients)
        .map(|(&j, &a)| a * values[j])
        .sum()
}

#[derive(Debug)]
struct CurvatureChannels {
    dofs: [Option<usize>; 9],
    curvature: [[[f64; 9]; 3]; 3],
    factor: [[f64; 6]; 3],
    weight: f64,
}
impl CurvatureChannels {
    fn apply(&self, values: &[f64]) -> [[f64; 6]; 3] {
        let local: [f64; 9] = self.dofs.map(|dof| dof.map_or(0.0, |i| values[i]));
        self.project(local)
    }
    fn increment(&self, initial: &[f64], endpoint: &[f64]) -> [[f64; 6]; 3] {
        let local = self
            .dofs
            .map(|dof| dof.map_or(0.0, |i| endpoint[i] - initial[i]));
        self.project(local)
    }
    fn project(&self, local: [f64; 9]) -> [[f64; 6]; 3] {
        self.curvature.map(|b| {
            let strain: [f64; 3] = b.map(|r| r.iter().zip(local).map(|(a, x)| a * x).sum());
            core::array::from_fn(|j| (0..3).map(|i| self.factor[i][j] * strain[i]).sum())
        })
    }
}

/// One spatial excess bending stiffness [N m] and its positive relaxation time.
/// Different times remain distinct; full directional and cross-coordinate
/// stiffness survives. The supplied plate K must contain EQUILIBRIUM stiffness:
/// adding the branch again to K would count its instantaneous stiffness twice.
#[derive(Debug)]
pub struct PlateRelaxationBranch<'a> {
    // Reuse the admitted original DKT assembly, whose numeric Gram operation
    // is identical. Only this wrapper assigns stiffness rather than viscosity
    // units; no viscosity source or fitted damping ratio is manufactured.
    spatial: PlateBendingViscosity<'a>,
    channels: Vec<CurvatureChannels>,
    relaxation_time_s: f64,
}
impl<'a> PlateRelaxationBranch<'a> {
    /// Assemble one branch on the same original mesh/support map as `model`.
    /// Each triangle requires an explicit passive engineering-curvature tensor
    /// in mesh axes; zero tensors are allowed. Geometry identity remains the
    /// caller's responsibility, as for the original spatial plate assembly.
    /// `budget` also bounds 81 curvature and 18 factor scalars per triangle.
    pub fn assemble<C: FnMut() -> bool>(
        mesh: &PlateMesh,
        model: &'a PlateModel,
        excess_bending_n_m: Vec<[f64; 9]>,
        relaxation_time_s: f64,
        budget: PlateViscosityBudget,
        cancelled: &mut C,
    ) -> Result<Self, Error> {
        poll(cancelled)?;
        if !relaxation_time_s.is_finite() || relaxation_time_s <= 0.0 {
            return Err(Error::Invalid(
                "finite positive Maxwell relaxation time required",
            ));
        }
        if mesh.tris.len() > budget.max_elements || excess_bending_n_m.len() != mesh.tris.len() {
            return Err(Error::Invalid(
                "bounded complete per-triangle Maxwell resultants required",
            ));
        }
        let sections = excess_bending_n_m
            .iter()
            .copied()
            .map(ViscousPlateSection::from_resultant)
            .collect::<Result<Vec<_>, _>>()
            .map_err(Error::Plate)?;
        let spatial = PlateBendingViscosity::assemble(mesh, model, sections, budget, cancelled)
            .map_err(Error::Plate)?;
        let mut channels = Vec::with_capacity(mesh.tris.len());
        for (element, (tri, tensor)) in mesh.tris.iter().zip(excess_bending_n_m).enumerate() {
            poll(cancelled)?;
            let (curvature, area2) = dkt_curvature_quadrature(
                &tri.map(|i| mesh.nodes[i].0),
                &tri.map(|i| mesh.nodes[i].1),
                element,
            )
            .map_err(Error::Plate)?;
            let factor = factor_passive_resultant(&tensor)
                .map_err(|_| Error::Invalid("passive Maxwell bending tensor required"))?;
            channels.push(CurvatureChannels {
                dofs: core::array::from_fn(|i| model.dof_map[3 * tri[i / 3] + i % 3]),
                curvature,
                factor,
                weight: area2 / 6.0,
            });
        }
        poll(cancelled)?;
        Ok(Self {
            spatial,
            channels,
            relaxation_time_s,
        })
    }
    pub fn excess_stiffness(&self) -> &Csr {
        self.spatial.matrix()
    }
    pub const fn relaxation_time_s(&self) -> f64 {
        self.relaxation_time_s
    }
    fn storage(&self, values: &[f64]) -> f64 {
        self.channels
            .iter()
            .map(|element| {
                0.5 * element.weight
                    * element
                        .apply(values)
                        .iter()
                        .flatten()
                        .map(|q| q * q)
                        .sum::<f64>()
            })
            .sum()
    }
}

/// Bounds checked before effective-matrix or state/workspace allocation.
#[derive(Debug, Clone, Copy)]
pub struct PlateRelaxationBudget {
    pub dynamics: PlateDynamicsBudget,
    pub max_branches: usize,
    /// Combined M/K/branch entries BEFORE duplicate effective entries merge.
    pub max_effective_entries: usize,
    /// Seven spatial vectors and two branch-history arrays: `(7+2*b)*n`.
    pub max_state_scalars: usize,
}

/// Prepared fixed-step sparse structural/constitutive dynamics. No dense
/// spatial matrix, modal diagonalization, or repeated factorization is used.
#[derive(Debug)]
pub struct RelaxingPlate<'a> {
    model: &'a PlateModel,
    branches: Vec<PlateRelaxationBranch<'a>>,
    ramps: Vec<PronyRamp>,
    effective: Csr,
    factor: LdltFactor,
    step_s: f64,
    solve_relative_tolerance: f64,
}
impl<'a> RelaxingPlate<'a> {
    pub fn prepare<C: FnMut() -> bool>(
        model: &'a PlateModel,
        branches: Vec<PlateRelaxationBranch<'a>>,
        step_s: f64,
        solve_relative_tolerance: f64,
        budget: PlateRelaxationBudget,
        factor_options: &LdltOptions,
        cancelled: &mut C,
    ) -> Result<Self, Error> {
        poll(cancelled)?;
        if !step_s.is_finite()
            || step_s <= 0.0
            || !solve_relative_tolerance.is_finite()
            || solve_relative_tolerance <= 0.0
            || solve_relative_tolerance >= 1.0
        {
            return Err(Error::Invalid(
                "positive finite time step and solve tolerance below one required",
            ));
        }
        let n = model.free;
        if branches.len() > budget.max_branches
            || branches
                .len()
                .checked_mul(2)
                .and_then(|b| b.checked_add(7))
                .and_then(|v| v.checked_mul(n))
                .is_none_or(|v| v > budget.max_state_scalars)
        {
            return Err(Error::Budget("Maxwell branch/state scalar cap exceeded"));
        }
        let mut entries = model
            .m
            .nnz()
            .checked_add(model.k.nnz())
            .ok_or(Error::Budget("effective operator extent overflow"))?;
        for branch in &branches {
            if !core::ptr::eq(branch.spatial.model(), model) {
                return Err(Error::Invalid(
                    "Maxwell branch belongs to a different original plate",
                ));
            }
            entries = entries
                .checked_add(branch.excess_stiffness().nnz())
                .ok_or(Error::Budget("effective operator extent overflow"))?;
        }
        if entries > budget.max_effective_entries {
            return Err(Error::Budget("effective operator entry cap exceeded"));
        }
        // The established original-coordinate admission keeps full consistent
        // inertia and its factor-byte cap, without adding Rayleigh damping.
        PlateDynamics::new_with_mass_options(
            model,
            PlateDynamicsParameters {
                stiffness_scale: 1.0,
                mass_scale: 1.0,
                mass_damping_per_s: 0.0,
                stiffness_damping_s: 0.0,
            },
            &NoPlateLoad,
            budget.dynamics,
            factor_options,
            cancelled,
        )
        .map_err(Error::Dynamics)?;
        let ramps: Vec<_> = branches
            .iter()
            .map(|b| PronyRamp::new(step_s / b.relaxation_time_s))
            .collect();
        let mut effective = Coo::new(n, n);
        for i in 0..n {
            poll(cancelled)?;
            for (matrix, weight) in [(&model.m, 1.0), (&model.k, 0.25 * step_s * step_s)] {
                let (columns, values) = matrix.row(i);
                for (&j, &value) in columns.iter().zip(values) {
                    let mirror = matrix.get(j, i);
                    if (value - mirror).abs() > 64.0 * f64::EPSILON * value.abs().max(mirror.abs())
                    {
                        return Err(Error::Invalid("symmetric physical operators required"));
                    }
                    effective.push(i, j, weight * value);
                }
            }
            for (branch, ramp) in branches.iter().zip(&ramps) {
                let weight = 0.5 * step_s * step_s * ramp.apply(0.0, 1.0).1;
                let (columns, values) = branch.excess_stiffness().row(i);
                for (&j, &value) in columns.iter().zip(values) {
                    effective.push(i, j, weight * value);
                }
            }
        }
        let effective = effective.assemble();
        if (0..n).any(|i| !finite(effective.row(i).1)) {
            return Err(Error::Invalid("effective Maxwell operator overflows"));
        }
        poll(cancelled)?;
        let symbolic =
            SymbolicLdlt::analyze(&effective, DirectOrdering::Amd).map_err(Error::Factor)?;
        poll(cancelled)?;
        let factor = symbolic
            .factor(&effective, factor_options)
            .map_err(Error::Factor)?;
        if factor.inertia().negative != 0 {
            return Err(Error::Invalid(
                "time step gives a nonpositive effective structural operator",
            ));
        }
        poll(cancelled)?;
        Ok(Self {
            model,
            branches,
            ramps,
            effective,
            factor,
            step_s,
            solve_relative_tolerance,
        })
    }
    pub const fn step_s(&self) -> f64 {
        self.step_s
    }
    pub fn branches(&self) -> &[PlateRelaxationBranch<'a>] {
        &self.branches
    }

    /// Initialize each branch's spring displacement/slope coordinates,
    /// flattened in branch then original reduced-coordinate order. Physical
    /// branch curvatures are projected from these coordinates. Zero histories
    /// mean a relaxed initial shape; copies of `displacement` mean an
    /// instantaneous change from the original zero shape. No prehistory is inferred.
    pub fn state<'s>(
        &'s self,
        displacement: &[f64],
        velocity: &[f64],
        branch_strains: &[f64],
    ) -> Result<PlateRelaxationState<'s, 'a>, Error> {
        let n = self.model.free;
        if displacement.len() != n
            || velocity.len() != n
            || branch_strains.len() != n * self.branches.len()
            || !finite(displacement)
            || !finite(velocity)
            || !finite(branch_strains)
        {
            return Err(Error::Invalid(
                "finite compatible original-coordinate state required",
            ));
        }
        let state = PlateRelaxationState {
            owner: self,
            displacement: displacement.to_vec(),
            velocity: velocity.to_vec(),
            history: branch_strains.to_vec(),
            next_displacement: vec![0.0; n],
            next_velocity: vec![0.0; n],
            next_history: vec![0.0; branch_strains.len()],
            rhs: vec![0.0; n],
            mean_velocity: vec![0.0; n],
            increment: vec![0.0; n],
            time_s: 0.0,
            dissipated_j: 0.0,
            supplied_work_j: 0.0,
        };
        if !state.mechanical_energy_j().is_finite() {
            return Err(Error::Invalid("initial physical energy overflows"));
        }
        Ok(state)
    }
}

/// Accepted step's physical work and energy increments [J]. Storage changes
/// use the actual rounded endpoint states. Heat and source work integrate the
/// intended affine path, so state rounding appears in the separate residual;
/// it is never assigned to material heat or hidden in a storage increment.
#[derive(Debug, Clone, Copy)]
pub struct PlateRelaxationStep {
    pub supplied_work_j: f64,
    pub dissipated_j: f64,
    pub kinetic_change_j: f64,
    pub stored_change_j: f64,
    pub balance_residual_j: f64,
    pub relative_solve_residual: f64,
}

/// A resumable material/structural state tied to its immutable prepared owner.
/// Workspace is allocated once here; successful or refused steps allocate none.
#[derive(Debug)]
pub struct PlateRelaxationState<'s, 'a> {
    owner: &'s RelaxingPlate<'a>,
    displacement: Vec<f64>,
    velocity: Vec<f64>,
    history: Vec<f64>,
    next_displacement: Vec<f64>,
    next_velocity: Vec<f64>,
    next_history: Vec<f64>,
    rhs: Vec<f64>,
    mean_velocity: Vec<f64>,
    increment: Vec<f64>,
    time_s: f64,
    dissipated_j: f64,
    supplied_work_j: f64,
}
impl PlateRelaxationState<'_, '_> {
    pub fn displacement(&self) -> &[f64] {
        &self.displacement
    }
    pub fn velocity(&self) -> &[f64] {
        &self.velocity
    }
    /// Branch spring coordinates in displacement/slope units, in the same
    /// flattened order admitted by [`RelaxingPlate::state`].
    pub fn branch_strains(&self) -> &[f64] {
        &self.history
    }
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub const fn dissipated_j(&self) -> f64 {
        self.dissipated_j
    }
    pub const fn supplied_work_j(&self) -> f64 {
        self.supplied_work_j
    }
    pub fn mechanical_energy_j(&self) -> f64 {
        let model = self.owner.model;
        let mut energy: f64 = (0..model.free)
            .map(|i| {
                0.5 * (self.velocity[i] * row(&model.m, i, &self.velocity)
                    + self.displacement[i] * row(&model.k, i, &self.displacement))
            })
            .sum();
        for (branch, q) in self
            .owner
            .branches
            .iter()
            .zip(self.history.chunks_exact(model.free))
        {
            energy += branch.storage(q);
        }
        energy
    }

    /// Advance under the supplied interval-mean force [N, N m] in the original
    /// reduced coordinates. Its work is exactly `mean_force dot delta_u` along
    /// this affine numerical displacement path. Supports and material parameters
    /// are fixed. Cancellation surrounds the bounded synchronous sparse solve.
    pub fn advance<C: FnMut() -> bool>(
        &mut self,
        mean_force: &[f64],
        cancelled: &mut C,
    ) -> Result<PlateRelaxationStep, Error> {
        poll(cancelled)?;
        let owner = self.owner;
        let model = owner.model;
        let n = model.free;
        let dt = owner.step_s;
        let next_time = self.time_s + dt;
        if mean_force.len() != n
            || !finite(mean_force)
            || !next_time.is_finite()
            || next_time <= self.time_s
        {
            return Err(Error::Invalid(
                "finite force and representable advancing time required",
            ));
        }
        for i in 0..n {
            poll(cancelled)?;
            let mut resistance = row(&model.k, i, &self.displacement) - mean_force[i];
            for ((branch, ramp), q) in owner
                .branches
                .iter()
                .zip(&owner.ramps)
                .zip(self.history.chunks_exact(n))
            {
                resistance += ramp.apply(1.0, 0.0).1 * row(branch.excess_stiffness(), i, q);
            }
            self.rhs[i] = row(&model.m, i, &self.velocity) - 0.5 * dt * resistance;
        }
        if !finite(&self.rhs) {
            return Err(Error::Invalid("Maxwell momentum forcing overflows"));
        }
        poll(cancelled)?;
        owner.factor.solve_into(&self.rhs, &mut self.mean_velocity);
        poll(cancelled)?;
        if !finite(&self.mean_velocity) {
            return Err(Error::Invalid("Maxwell velocity solve overflows"));
        }
        let mut relative = 0.0_f64;
        for i in 0..n {
            poll(cancelled)?;
            let (columns, values) = owner.effective.row(i);
            let magnitude = self.rhs[i].abs()
                + columns
                    .iter()
                    .zip(values)
                    .map(|(&j, &a)| (a * self.mean_velocity[j]).abs())
                    .sum::<f64>();
            let residual = (row(&owner.effective, i, &self.mean_velocity) - self.rhs[i]).abs();
            if !magnitude.is_finite() || !residual.is_finite() {
                return Err(Error::Invalid("Maxwell solve residual overflows"));
            }
            relative = relative.max(if magnitude == 0.0 {
                residual
            } else {
                residual / magnitude
            });
            self.increment[i] = dt * self.mean_velocity[i];
            self.next_displacement[i] = self.displacement[i] + self.increment[i];
            self.next_velocity[i] =
                self.mean_velocity[i] + (self.mean_velocity[i] - self.velocity[i]);
        }
        if relative > owner.solve_relative_tolerance {
            return Err(Error::SolveResidual {
                relative,
                allowed: owner.solve_relative_tolerance,
            });
        }
        let mut heat = 0.0;
        let mut stored = 0.0;
        for (j, (branch, ramp)) in owner.branches.iter().zip(&owner.ramps).enumerate() {
            let q = &self.history[j * n..(j + 1) * n];
            for i in 0..n {
                self.next_history[j * n + i] = ramp.apply(q[i], self.increment[i]).0;
            }
            for element in &branch.channels {
                poll(cancelled)?;
                let initial = element.apply(q);
                let increment = element.apply(&self.increment);
                let accepted_increment =
                    element.increment(q, &self.next_history[j * n..(j + 1) * n]);
                for ((&q0, &dq), &change) in initial
                    .iter()
                    .flatten()
                    .zip(increment.iter().flatten())
                    .zip(accepted_increment.iter().flatten())
                {
                    heat += element.weight * ramp.apply(q0, dq).2;
                    stored += element.weight * change * (q0 + 0.5 * change);
                }
            }
        }
        let mut work = 0.0;
        let mut kinetic = 0.0;
        for i in 0..n {
            poll(cancelled)?;
            work += mean_force[i] * self.increment[i];
            let (columns, values) = model.k.row(i);
            let elastic_change: f64 = columns
                .iter()
                .zip(values)
                .map(|(&j, &k)| k * (self.next_displacement[j] - self.displacement[j]))
                .sum();
            stored += (self.next_displacement[i] - self.displacement[i])
                * (row(&model.k, i, &self.displacement) + 0.5 * elastic_change);
            let (columns, values) = model.m.row(i);
            let momentum_change: f64 = columns
                .iter()
                .zip(values)
                .map(|(&j, &m)| m * (self.next_velocity[j] - self.velocity[j]))
                .sum();
            kinetic += (self.next_velocity[i] - self.velocity[i])
                * (row(&model.m, i, &self.velocity) + 0.5 * momentum_change);
        }
        let balance = kinetic + stored + heat - work;
        let total_heat = self.dissipated_j + heat;
        let total_work = self.supplied_work_j + work;
        if !finite(&self.next_displacement)
            || !finite(&self.next_velocity)
            || !finite(&self.next_history)
            || !finite(&[kinetic, stored, heat, work, balance, total_heat, total_work])
            || heat < 0.0
        {
            return Err(Error::Invalid("physical Maxwell state or energy overflows"));
        }
        poll(cancelled)?;
        core::mem::swap(&mut self.displacement, &mut self.next_displacement);
        core::mem::swap(&mut self.velocity, &mut self.next_velocity);
        core::mem::swap(&mut self.history, &mut self.next_history);
        self.time_s = next_time;
        self.dissipated_j = total_heat;
        self.supplied_work_j = total_work;
        Ok(PlateRelaxationStep {
            supplied_work_j: work,
            dissipated_j: heat,
            kinetic_change_j: kinetic,
            stored_change_j: stored,
            balance_residual_j: balance,
            relative_solve_residual: relative,
        })
    }
}
