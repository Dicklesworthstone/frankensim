//! Simultaneous, unregularized finite-patch tangential contact.
//!
//! Positive normal patches and geometry are supplied by the normal owner. The
//! original two-region return map supplies friction, storage and physical heat;
//! this module solves its coupled endpoint velocities, without smoothing away
//! stick/slip events. It is not a normal-complementarity or impact stepper.
//!
//! A caller supplies B = J M^(-1/2), giving passive mobility W = B B^T, including
//! rank-deficient shared-body contacts. Each contact uses [vx, vy, ell*omega]
//! and conjugate [Fx, Fy, M/ell], so every solver coordinate has units m/s.
//! Backward Euler solves v = v_free + dt W f(v). Its mechanical integration
//! loss is reported separately from the law's physical heat and elastic loss.

use fs_exec::CancelGate;
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
};
use fs_tribo::partial_slip::{
    GeneralizedWorkOwnership, NormalPatchView, PartialSlipCheckpoint, PartialSlipError,
    PartialSlipInterface, PartialSlipKinematics, PartialSlipLaw, PartialSlipState, PartialSlipStep,
    TangentFrame, TangentialWrench,
};
use std::cell::Cell;
use std::collections::{BTreeMap, BTreeSet};

/// Identity of this bounded nonsmooth tangential solver.
pub const NONSMOOTH_TANGENTIAL_MODEL_ID: &str = "fs-contact/coupled-partial-slip-backward-euler-v1";

/// One currently loaded contact. The normal owner retains its own solve gate.
pub struct NonsmoothContact<'a> {
    pub law: &'a PartialSlipLaw,
    pub patch: &'a NormalPatchView,
    pub interface: &'a PartialSlipInterface,
    pub frame: TangentFrame,
    /// Physical common frame spin, distinct from relative torsional slip.
    pub common_spin_rad_per_s: f64,
    /// Positive coordinate scale ell [m], with no change to physical torsion.
    pub torsional_length_m: f64,
    /// Positive reference speed used only to express physical slip as creepage.
    pub creepage_reference_speed_mps: f64,
    pub work: GeneralizedWorkOwnership,
}

/// Factor B = J M^(-1/2), in input contact order; three rows per contact.
/// The caller owns the actual inertia and contact Jacobian. Zero columns are
/// allowed for prescribed relative motion. Finite entries have units kg^-1/2.
pub struct TangentialMobility<'a> {
    pub mechanical_dofs: usize,
    pub factor: &'a [f64],
}

/// Numerical budgets and independent physical acceptance tolerances.
#[derive(Debug, Clone, Copy)]
pub struct NonsmoothControls {
    pub maximum_evaluations: usize,
    pub maximum_newton_iterations: usize,
    pub maximum_mechanical_dofs: usize,
    pub difference_step_mps: f64,
    pub velocity_tolerance_mps: f64,
    pub cone_tolerance: f64,
    pub energy_tolerance_j: f64,
    pub newton: NewtonKrylovConfig,
}

/// Complete accepted history, including the contact frame and consumed work.
/// Clone this value for checkpoint/restart; trials never mutate it.
#[derive(Debug, Clone, PartialEq)]
pub struct NonsmoothTangentialState {
    version: u64,
    time_s: f64,
    contacts: BTreeMap<String, PartialSlipCheckpoint>,
    committed_work: BTreeSet<GeneralizedWorkOwnership>,
    maximum_contacts: usize,
    maximum_work_keys: usize,
}

#[derive(Debug)]
pub enum NonsmoothTangentialError {
    Invalid(&'static str),
    Law(PartialSlipError),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    Cancelled,
    EvaluationBudget,
    Capacity,
    DuplicateWork,
    StateMismatch,
    /// A loaded history cannot silently disappear when its patch is omitted.
    MissingActiveContact,
    /// Independent force, cone, flow, or energy checks rejected the candidate.
    Audit(NonsmoothAudit),
}
impl From<PartialSlipError> for NonsmoothTangentialError {
    fn from(error: PartialSlipError) -> Self {
        Self::Law(error)
    }
}
impl From<NewtonError> for NonsmoothTangentialError {
    fn from(error: NewtonError) -> Self {
        Self::Newton(error)
    }
}
impl core::fmt::Display for NonsmoothTangentialError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "nonsmooth tangential contact: {self:?}")
    }
}
impl std::error::Error for NonsmoothTangentialError {}

/// Independent endpoint, regional KKT, and whole-island energy checks.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct NonsmoothAudit {
    pub velocity_residual_mps: f64,
    pub cone_violation: f64,
    /// Regional plastic flow residual and cone complementarity product [J].
    pub complementarity_residual_j: f64,
    /// Work transferred to free mechanical motion and inertial recoil [J].
    /// For dynamic bodies this is their kinetic change from the free endpoint;
    /// prescribed motion instead owns the corresponding signed drive work.
    pub mechanical_work_j: f64,
    pub storage_change_j: f64,
    /// Explicit signed work from changing coefficients at retained old strain.
    /// The coupled coefficient driver must balance this port; it is not heat.
    pub coefficient_work_j: f64,
    pub physical_heat_j: f64,
    pub elastic_numerical_loss_j: f64,
    pub mechanical_numerical_loss_j: f64,
    pub energy_residual_j: f64,
}

/// A solved contact in canonical patch-id order.
#[derive(Debug, Clone, PartialEq)]
pub struct NonsmoothContactResult {
    pub patch_id: String,
    /// [vx, vy, ell*omega], m/s, in this contact's declared frame.
    pub endpoint_velocity_mps: [f64; 3],
    /// [dt Fx, dt Fy, dt M/ell], N s.
    pub impulse_ns: [f64; 3],
    pub action: TangentialWrench,
    pub reaction: TangentialWrench,
    pub law: PartialSlipStep,
}

/// Immutable proposal. Inspect the physical solution, then commit all history.
/// Dropping it rolls back all contact work and state.
#[derive(Debug, Clone)]
pub struct NonsmoothTangentialTrial {
    prior: NonsmoothTangentialState,
    next: NonsmoothTangentialState,
    contacts: Vec<NonsmoothContactResult>,
    audit: NonsmoothAudit,
    nonlinear: Option<NewtonReport>,
    evaluations: usize,
    mass_weighted_velocity_increment: Vec<f64>,
}
impl NonsmoothTangentialTrial {
    pub fn contacts(&self) -> &[NonsmoothContactResult] {
        &self.contacts
    }
    pub const fn audit(&self) -> NonsmoothAudit {
        self.audit
    }
    pub fn nonlinear(&self) -> Option<&NewtonReport> {
        self.nonlinear.as_ref()
    }
    pub const fn evaluations(&self) -> usize {
        self.evaluations
    }
    /// B^T impulse. The mechanical owner applies M^(-1/2) to obtain its actual
    /// generalized-velocity increment and accepts that state with this trial.
    pub fn mass_weighted_velocity_increment(&self) -> &[f64] {
        &self.mass_weighted_velocity_increment
    }
}

impl NonsmoothTangentialState {
    pub fn new(
        maximum_contacts: usize,
        maximum_work_keys: usize,
    ) -> Result<Self, NonsmoothTangentialError> {
        if maximum_contacts == 0 || maximum_work_keys == 0 {
            return Err(NonsmoothTangentialError::Invalid(
                "positive state budgets required",
            ));
        }
        Ok(Self {
            version: 0,
            time_s: 0.0,
            contacts: BTreeMap::new(),
            committed_work: BTreeSet::new(),
            maximum_contacts,
            maximum_work_keys,
        })
    }
    pub const fn version(&self) -> u64 {
        self.version
    }
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub fn contact_checkpoint(&self, patch_id: &str) -> Option<&PartialSlipCheckpoint> {
        self.contacts.get(patch_id)
    }

    /// Solve all loaded contacts simultaneously from this exact accepted state.
    /// Input permutation only permutes B's rows and free velocities; reduction
    /// order is canonical by patch id. Every nonlinear probe reuses old history.
    ///
    /// Explicit coefficient changes use the original law's continuation and
    /// retain regional memory. The signed parameter-work port must be balanced
    /// by the coupled coefficient driver; it is never silently added as heat.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn prepare(
        &self,
        contacts: &[NonsmoothContact<'_>],
        free_velocity_mps: &[[f64; 3]],
        mobility: TangentialMobility<'_>,
        dt_s: f64,
        controls: NonsmoothControls,
        gate: &CancelGate,
    ) -> Result<NonsmoothTangentialTrial, NonsmoothTangentialError> {
        poll(gate)?;
        controls.newton.validate()?;
        let count = contacts.len();
        let dimension = count
            .checked_mul(3)
            .ok_or(NonsmoothTangentialError::Capacity)?;
        let factor_len = dimension
            .checked_mul(mobility.mechanical_dofs)
            .ok_or(NonsmoothTangentialError::Capacity)?;
        if count > self.maximum_contacts
            || mobility.mechanical_dofs > controls.maximum_mechanical_dofs
            || self
                .committed_work
                .len()
                .checked_add(count)
                .is_none_or(|n| n > self.maximum_work_keys)
        {
            return Err(NonsmoothTangentialError::Capacity);
        }
        let target = self.time_s + dt_s;
        if !dt_s.is_finite()
            || dt_s <= 0.0
            || !target.is_finite()
            || target <= self.time_s
            || free_velocity_mps.len() != count
            || mobility.factor.len() != factor_len
            || free_velocity_mps
                .iter()
                .flatten()
                .chain(mobility.factor)
                .any(|v| !v.is_finite())
            || controls.maximum_evaluations == 0
            || controls.maximum_newton_iterations == 0
            || !positive(controls.difference_step_mps)
            || !positive(controls.velocity_tolerance_mps)
            || !positive(controls.cone_tolerance)
            || !positive(controls.energy_tolerance_j)
        {
            return Err(NonsmoothTangentialError::Invalid(
                "finite complete coordinates and positive solve controls required",
            ));
        }
        let mut order: Vec<_> = (0..count).collect();
        order.sort_by(|&a, &b| {
            contacts[a]
                .patch
                .patch_id()
                .cmp(contacts[b].patch.patch_id())
        });
        let mut accepted = Vec::with_capacity(count);
        let mut free = Vec::with_capacity(dimension);
        let mut factor = Vec::with_capacity(factor_len);
        for (slot, &i) in order.iter().enumerate() {
            poll(gate)?;
            let contact = &contacts[i];
            if (slot > 0 && contacts[order[slot - 1]].patch.patch_id() == contact.patch.patch_id())
                || contact.work.patch_id() != contact.patch.patch_id()
                || (slot > 0 && contact.work.interval_id() != contacts[order[0]].work.interval_id())
                || !contact.common_spin_rad_per_s.is_finite()
                || !positive(contact.torsional_length_m)
                || !positive(contact.creepage_reference_speed_mps)
            {
                return Err(NonsmoothTangentialError::Invalid(
                    "distinct contact identities, one interval, and finite physical scales required",
                ));
            }
            if self.committed_work.contains(&contact.work) {
                return Err(NonsmoothTangentialError::DuplicateWork);
            }
            accepted.push(match self.contacts.get(contact.patch.patch_id()) {
                Some(checkpoint) => {
                    checkpoint.law().restore_checkpoint(
                        contact.patch,
                        contact.interface,
                        checkpoint,
                    )?;
                    checkpoint.clone()
                }
                None => PartialSlipCheckpoint::new(
                    contact.patch.clone(),
                    contact.interface.clone(),
                    contact.law.clone(),
                    PartialSlipState::zero(),
                )?,
            });
            free.extend_from_slice(&free_velocity_mps[i]);
            let start = 3 * i * mobility.mechanical_dofs;
            factor.extend_from_slice(&mobility.factor[start..start + 3 * mobility.mechanical_dofs]);
        }
        if self
            .contacts
            .keys()
            .any(|key| !order.iter().any(|&i| contacts[i].patch.patch_id() == key))
        {
            return Err(NonsmoothTangentialError::MissingActiveContact);
        }
        let evaluations = Cell::new(0usize);
        let evaluate = |velocity: &[f64]| -> Result<
            (Vec<f64>, Vec<PartialSlipStep>, Vec<f64>),
            NonsmoothTangentialError,
        > {
            poll(gate)?;
            if evaluations.get() >= controls.maximum_evaluations {
                return Err(NonsmoothTangentialError::EvaluationBudget);
            }
            evaluations.set(evaluations.get() + 1);
            let mut steps = Vec::with_capacity(count);
            let mut impulse = Vec::with_capacity(dimension);
            for (slot, &i) in order.iter().enumerate() {
                poll(gate)?;
                let contact = &contacts[i];
                let v = &velocity[3 * slot..3 * slot + 3];
                let step = contact.law.advance_with_coefficients(
                    contact.patch,
                    contact.interface,
                    contact.frame,
                    PartialSlipKinematics {
                        creepage: [
                            v[0] / contact.creepage_reference_speed_mps,
                            v[1] / contact.creepage_reference_speed_mps,
                        ],
                        rolling_speed_mps: contact.creepage_reference_speed_mps,
                        torsional_spin_rad_per_s: v[2] / contact.torsional_length_m,
                        contact_spin_rad_per_s: contact.common_spin_rad_per_s,
                        dt_s,
                    },
                    &contact.work,
                    &accepted[slot],
                )?;
                impulse.extend_from_slice(&[
                    dt_s * step.tangent_force_n[0],
                    dt_s * step.tangent_force_n[1],
                    dt_s * step.torsional_moment_nm / contact.torsional_length_m,
                ]);
                steps.push(step);
            }
            let weighted = transpose_product(&factor, mobility.mechanical_dofs, &impulse, gate)?;
            let delta = product(
                &factor,
                mobility.mechanical_dofs,
                dimension,
                &weighted,
                gate,
            )?;
            let residual: Vec<_> = velocity
                .iter()
                .zip(&free)
                .zip(&delta)
                .map(|((&v, &v0), &dv)| v - v0 - dv)
                .collect();
            if impulse
                .iter()
                .chain(&residual)
                .chain(&weighted)
                .any(|v| !v.is_finite())
            {
                return Err(NonsmoothTangentialError::Invalid(
                    "unrepresentable contact response",
                ));
            }
            Ok((residual, steps, weighted))
        };
        let (velocity, nonlinear) = if dimension == 0 {
            (Vec::new(), None)
        } else {
            let problem = FallibleFiniteDifference::new(
                dimension,
                controls.difference_step_mps,
                |x: &[f64]| evaluate(x).map(|(residual, _, _)| residual),
            )
            .map_err(numerical)?;
            let construction = NewtonKrylovState::new(&problem, free.clone(), controls.newton);
            if let Some(error) = problem.take_failure() {
                return Err(numerical(error));
            }
            let mut solver = construction?;
            let result =
                solver.run_controlled(&problem, controls.maximum_newton_iterations, &mut || {
                    gate.is_requested() || problem.has_failed()
                });
            if let Some(error) = problem.take_failure() {
                return Err(numerical(error));
            }
            poll(gate)?;
            if result.progress != SolverRunProgress::Complete || !result.report.converged {
                return Err(NonsmoothTangentialError::NotConverged(result.report));
            }
            (solver.x, Some(result.report))
        };
        let (residual, steps, weighted) = evaluate(&velocity)?;
        let mut audit = NonsmoothAudit {
            velocity_residual_mps: residual.iter().copied().map(f64::abs).fold(0.0, f64::max),
            mechanical_numerical_loss_j: 0.5 * weighted.iter().map(|x| x * x).sum::<f64>(),
            ..NonsmoothAudit::default()
        };
        let mut solved = Vec::with_capacity(count);
        let mut next = self.clone();
        for (slot, (step, &i)) in steps.into_iter().zip(&order).enumerate() {
            let contact = &contacts[i];
            let v = [
                velocity[3 * slot],
                velocity[3 * slot + 1],
                velocity[3 * slot + 2],
            ];
            let impulse = [
                dt_s * step.tangent_force_n[0],
                dt_s * step.tangent_force_n[1],
                dt_s * step.torsional_moment_nm / contact.torsional_length_m,
            ];
            audit.mechanical_work_j += (0..3).map(|j| impulse[j] * free[3 * slot + j]).sum::<f64>();
            audit.storage_change_j += step.storage_change_j;
            audit.coefficient_work_j += step.generalized_work.coefficient_work_j;
            audit.physical_heat_j += step.dissipation.heat_j;
            audit.elastic_numerical_loss_j += step.generalized_work.numerical_dissipation_j;
            let (cone, flow) = regional_kkt(contact.law, &step);
            audit.cone_violation = audit.cone_violation.max(cone);
            audit.complementarity_residual_j = audit.complementarity_residual_j.max(flow);
            next.contacts
                .insert(contact.patch.patch_id().to_owned(), step.checkpoint.clone());
            next.committed_work.insert(contact.work.clone());
            solved.push(NonsmoothContactResult {
                patch_id: contact.patch.patch_id().to_owned(),
                endpoint_velocity_mps: v,
                impulse_ns: impulse,
                action: step.wrench,
                reaction: TangentialWrench {
                    force_n: step.wrench.force_n.map(|x| -x),
                    torque_nm: step.wrench.torque_nm.map(|x| -x),
                },
                law: step,
            });
        }
        audit.mechanical_work_j += audit.mechanical_numerical_loss_j;
        audit.energy_residual_j = audit.mechanical_work_j
            + audit.storage_change_j
            + audit.physical_heat_j
            + audit.elastic_numerical_loss_j
            + audit.mechanical_numerical_loss_j
            - audit.coefficient_work_j;
        if [
            audit.velocity_residual_mps,
            audit.cone_violation,
            audit.complementarity_residual_j,
            audit.mechanical_work_j,
            audit.storage_change_j,
            audit.coefficient_work_j,
            audit.physical_heat_j,
            audit.elastic_numerical_loss_j,
            audit.mechanical_numerical_loss_j,
            audit.energy_residual_j,
        ]
        .iter()
        .any(|x| !x.is_finite())
            || audit.velocity_residual_mps > controls.velocity_tolerance_mps
            || audit.cone_violation > controls.cone_tolerance
            || audit.complementarity_residual_j > controls.energy_tolerance_j
            || audit.energy_residual_j.abs() > controls.energy_tolerance_j
        {
            return Err(NonsmoothTangentialError::Audit(audit));
        }
        next.version = self
            .version
            .checked_add(1)
            .ok_or(NonsmoothTangentialError::Capacity)?;
        next.time_s = target;
        poll(gate)?;
        Ok(NonsmoothTangentialTrial {
            prior: self.clone(),
            next,
            contacts: solved,
            audit,
            nonlinear,
            evaluations: evaluations.get(),
            mass_weighted_velocity_increment: weighted,
        })
    }

    /// Publish every solved contact and its work once. The caller must accept
    /// the corresponding mechanical velocity update in the same transaction.
    pub fn commit(
        &mut self,
        trial: NonsmoothTangentialTrial,
    ) -> Result<(), NonsmoothTangentialError> {
        if trial.prior != *self {
            return Err(NonsmoothTangentialError::StateMismatch);
        }
        *self = trial.next;
        Ok(())
    }

    /// Remove a separated contact only after its retained elastic energy has
    /// been explicitly exported to the caller's mechanical/unloading account.
    /// This export is not heat. Recontact then starts with empty shear/twist.
    /// Work keys survive release, preventing duplicate interval consumption.
    pub fn release_contact(
        &mut self,
        work: GeneralizedWorkOwnership,
        exported_storage_j: f64,
        energy_tolerance_j: f64,
    ) -> Result<f64, NonsmoothTangentialError> {
        if !exported_storage_j.is_finite()
            || exported_storage_j < 0.0
            || !positive(energy_tolerance_j)
        {
            return Err(NonsmoothTangentialError::Invalid(
                "finite nonnegative release energy required",
            ));
        }
        if self.committed_work.contains(&work) {
            return Err(NonsmoothTangentialError::DuplicateWork);
        }
        if self.committed_work.len() >= self.maximum_work_keys {
            return Err(NonsmoothTangentialError::Capacity);
        }
        let checkpoint =
            self.contacts
                .get(work.patch_id())
                .ok_or(NonsmoothTangentialError::Invalid(
                    "released contact is not retained",
                ))?;
        let state = checkpoint.state();
        let p = checkpoint.law().parameters();
        let mut stored = 0.0;
        for (region, weight) in [
            1.0 - p.partial_slip_hardening_fraction,
            p.partial_slip_hardening_fraction,
        ]
        .into_iter()
        .enumerate()
        {
            stored += 0.5
                * weight
                * (p.tangential_stiffness_n_per_m
                    * state.region_elastic_displacement_m()[region]
                        .iter()
                        .map(|v| v * v)
                        .sum::<f64>()
                    + p.torsional_stiffness_nm_per_rad
                        * state.region_elastic_twist_rad()[region].powi(2));
        }
        if !stored.is_finite() || (stored - exported_storage_j).abs() > energy_tolerance_j {
            return Err(NonsmoothTangentialError::Invalid(
                "release must balance all retained elastic storage",
            ));
        }
        let version = self
            .version
            .checked_add(1)
            .ok_or(NonsmoothTangentialError::Capacity)?;
        self.contacts.remove(work.patch_id());
        self.committed_work.insert(work);
        self.version = version;
        Ok(stored)
    }
}

fn numerical(error: FiniteDifferenceError<NonsmoothTangentialError>) -> NonsmoothTangentialError {
    match error {
        FiniteDifferenceError::Source(source) => source,
        FiniteDifferenceError::Invalid(message) => NonsmoothTangentialError::Invalid(message),
    }
}
fn poll(gate: &CancelGate) -> Result<(), NonsmoothTangentialError> {
    if gate.is_requested() {
        Err(NonsmoothTangentialError::Cancelled)
    } else {
        Ok(())
    }
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn transpose_product(
    factor: &[f64],
    columns: usize,
    x: &[f64],
    gate: &CancelGate,
) -> Result<Vec<f64>, NonsmoothTangentialError> {
    let mut out = vec![0.0; columns];
    for (row, &value) in x.iter().enumerate() {
        poll(gate)?;
        for col in 0..columns {
            out[col] += factor[row * columns + col] * value;
        }
    }
    Ok(out)
}
fn product(
    factor: &[f64],
    columns: usize,
    rows: usize,
    x: &[f64],
    gate: &CancelGate,
) -> Result<Vec<f64>, NonsmoothTangentialError> {
    (0..rows)
        .map(|row| {
            poll(gate)?;
            Ok((0..columns)
                .map(|col| factor[row * columns + col] * x[col])
                .sum())
        })
        .collect()
}

/// Check the original regional return map's cone and associated-flow KKT.
/// q=effort/capacity is dimensionless; z=plastic_increment*capacity is J.
/// For a slipping region z=lambda*q, lambda>=0, and lambda*(1-|q|)=0.
fn regional_kkt(law: &PartialSlipLaw, step: &PartialSlipStep) -> (f64, f64) {
    let p = law.parameters();
    let hard = p.partial_slip_hardening_fraction;
    let fractions = [
        p.partial_slip_onset_fraction * (1.0 - hard),
        1.0 - p.partial_slip_onset_fraction * (1.0 - hard),
    ];
    let mut cone = 0.0_f64;
    let mut flow = 0.0_f64;
    for (region, weight) in [1.0 - hard, hard].into_iter().enumerate() {
        let kinetic = step.next_state.sliding_regions()[region];
        let capacity = if kinetic {
            [
                step.capacity.kinetic_force_n,
                step.capacity.kinetic_force_n,
                step.capacity.kinetic_torque_nm,
            ]
        } else {
            [
                step.capacity.static_force_n,
                step.capacity.static_force_n,
                step.capacity.static_torque_nm,
            ]
        }
        .map(|v| v * fractions[region]);
        let strain = [
            step.next_state.region_elastic_displacement_m()[region][0],
            step.next_state.region_elastic_displacement_m()[region][1],
            step.next_state.region_elastic_twist_rad()[region],
        ];
        let stiffness = [
            p.tangential_stiffness_n_per_m * weight,
            p.tangential_stiffness_n_per_m * weight,
            p.torsional_stiffness_nm_per_rad * weight,
        ];
        let mut q = [0.0; 3];
        let mut z = [0.0; 3];
        for j in 0..3 {
            q[j] = stiffness[j] * strain[j] / capacity[j];
            z[j] = step.region_plastic_increment[region][j] * capacity[j];
        }
        let norm = q[0].hypot(q[1]).hypot(q[2]);
        let denominator = q.iter().map(|x| x * x).sum::<f64>();
        let lambda = if denominator > 0.0 {
            q.iter().zip(z).map(|(a, b)| a * b).sum::<f64>() / denominator
        } else {
            0.0
        };
        cone = cone.max((norm - 1.0).max(0.0));
        flow = flow
            .max((-lambda).max(0.0))
            .max((lambda * (1.0 - norm)).abs());
        for j in 0..3 {
            flow = flow.max((z[j] - lambda * q[j]).abs());
        }
    }
    (cone, flow)
}
