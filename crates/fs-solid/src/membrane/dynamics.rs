//! Nonlinear elastic sheet motion using the shared physical-momentum Verlet map.
//!
//! The original membrane energy supplies BOTH force kicks, including geometric
//! stiffening and every admitted nonlinear material law. Only free coordinates
//! are integrated; stationary supports and the reference mass are fixed.
//! The original entry uses positive lumping. An explicitly bounded consistent
//! image retains inter-node inertia, factors each scalar component once, and
//! includes the prescribed/free mass coupling in the reported support impulse.
//! No equilibrium tangent, authored spring, numerical damping or mass rescaling
//! replaces the material. Work arrays are retained across steps.
//!
//! This explicit, second-order lane needs a resolved time step. Its signed
//! energy/work defect is integration error, not dissipated heat. Symplecticity
//! does not establish nonlinear stability or mesh/time convergence. Moving
//! supports, contact, and material history are not included. Held signed
//! pressure uses current-area forces at both kicks and exact work over the
//! resolved straight drift. Nonuniform material profiles or a free perimeter
//! can be nonconservative; no symplectic or path-independent pressure-work claim
//! applies in those cases.

use super::{
    MembraneError, MembranePatch, MembranePressure, finite_vectors, inertia::MembraneMassImage,
    poll,
};
use fs_solver::op::LinearOp;
use fs_time::{MomentumVerletError, verlet_momentum_step_with_displacement};
use std::cell::{Cell, RefCell};

mod inertia;
pub use inertia::{
    ConsistentMassBudget, MAX_CONSISTENT_AXIS_NODES, SparseConsistentMassBudget,
    SparseMassStatistics,
};
use inertia::{ConsistentPreparation, InverseMass};

fn kinetic(momentum: &[f64], inverse: &InverseMass, velocity: &[f64]) -> f64 {
    match inverse {
        InverseMass::Lumped(diagonal) => momentum
            .iter()
            .zip(diagonal)
            .map(|(p, inv)| 0.5 * p * (inv * p))
            .sum(),
        InverseMass::Consistent { .. } => momentum
            .iter()
            .zip(velocity)
            .map(|(p, v)| 0.5 * p * v)
            .sum(),
    }
}

// Polarize the original quadratic form: do not subtract large absolute
// kinetic energies to recover a small, resolved mechanical transfer.
fn kinetic_change(
    before: &[f64],
    after: &[f64],
    inverse: &InverseMass,
    free: &[usize],
    old_velocity: &[[f64; 3]],
    new_velocity: &[f64],
) -> f64 {
    match inverse {
        InverseMass::Lumped(diagonal) => before
            .iter()
            .zip(after)
            .zip(diagonal)
            .map(|((old, new), inv)| (new - old) * (inv * (0.5 * old + 0.5 * new)))
            .sum(),
        InverseMass::Consistent { .. } => before
            .iter()
            .zip(after)
            .zip(free)
            .zip(new_velocity)
            .map(|(((old, new), &dof), v)| {
                (new - old) * (0.5 * old_velocity[dof / 3][dof % 3] + 0.5 * v)
            })
            .sum(),
    }
}

/// Physical endpoint and signed work balance of one accepted step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MembraneStep {
    /// Accepted endpoint time [s]. Negative steps can replay conservative motion.
    pub time_s: f64,
    /// Original material strain energy [J], including any static prestretch.
    pub strain_energy_j: f64,
    /// Kinetic energy from the selected original reference mass image [J].
    pub kinetic_energy_j: f64,
    /// Work of the held load over the resolved kinematic drift [J], even when
    /// adding that displacement to the absolute position rounds it away.
    pub dead_load_work_j: f64,
    /// Work of held current-area pressure over the resolved straight drift [J].
    /// It is separate from dead-load work and from numerical integration error.
    pub pressure_work_j: f64,
    /// Change in total mechanical energy minus external work [J]. This is a
    /// signed numerical defect; no physical dissipation is introduced.
    pub energy_defect_j: f64,
    /// Net support impulse [N s], using the same two endpoint force kicks and
    /// the consistent image's prescribed/free inertial coupling when selected.
    /// Supports do no work because their displacement is held fixed.
    pub support_impulse_n_s: [f64; 3],
}

/// Resumable nonlinear motion bound to one immutable source/mesh owner.
///
/// Cloning is a complete checkpoint. The physical state is private and changes
/// only after both material evaluations, all finite checks and the final
/// cancellation poll succeed. Rejected trials may change scratch only.
#[derive(Clone, Debug)]
pub struct MembraneMotion<'a> {
    patch: &'a MembranePatch,
    fixed: Vec<bool>,
    free: Vec<usize>,
    inverse_mass: InverseMass,
    displacement_m: Vec<[f64; 3]>,
    velocity_m_s: Vec<[f64; 3]>,
    q: Vec<f64>,
    momentum: Vec<f64>,
    time_s: f64,
    strain_energy_j: f64,
    kinetic_energy_j: f64,
    trial_q: Vec<f64>,
    trial_momentum: Vec<f64>,
    // Only the consistent image needs a retained inverse-action output. The
    // original diagonal image keeps its allocation and arithmetic path.
    trial_velocity: Vec<f64>,
    trial_displacement: Vec<[f64; 3]>,
    gradient: Vec<[f64; 3]>,
    trial_gradient: Vec<[f64; 3]>,
    // Reused for endpoint pressure forces, then for the resolved full drift.
    pressure_scratch: Vec<[f64; 3]>,
    scratch: Vec<f64>,
}

impl<'a> MembraneMotion<'a> {
    /// Prepare reference inertia and stationary supports once. Prescribed
    /// entries are `(3*node+axis, displacement_m)`. Initial support displacement
    /// must match exactly and support velocity must be zero: admission cannot
    /// silently move a boundary or inject an initial impulse.
    pub fn new(
        patch: &'a MembranePatch,
        prescribed: &[(usize, f64)],
        displacement_m: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        time_s: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MembraneError> {
        Self::prepare(
            patch,
            prescribed,
            displacement_m,
            velocity_m_s,
            time_s,
            None,
            cancelled,
        )
    }

    /// Prepare exact P1 density-weighted consistent reference inertia for the
    /// SAME nonlinear momentum Verlet map. All density/velocity cross moments
    /// are retained before the stationary supports are eliminated.
    ///
    /// Each Cartesian component factors its scalar free-node mass block once;
    /// equal support patterns share the immutable Cholesky factor. The caller's
    /// explicit axis-node cap and `MAX_CONSISTENT_AXIS_NODES` bound that dense
    /// preparation and each later triangular solve. Warm steps allocate no mass
    /// matrix or solve workspace. This is a bounded dense inverse, not a sparse
    /// large-mesh solver or an automatic explicit-stability certificate.
    /// Initial velocity is converted to physical momentum with the original
    /// mass, then observed through this same inverse; the recovered velocity
    /// can differ from the supplied values by factorization roundoff.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_consistent_inertia(
        patch: &'a MembranePatch,
        prescribed: &[(usize, f64)],
        displacement_m: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        time_s: f64,
        budget: ConsistentMassBudget,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MembraneError> {
        Self::prepare(
            patch,
            prescribed,
            displacement_m,
            velocity_m_s,
            time_s,
            Some(ConsistentPreparation::Dense(budget)),
            cancelled,
        )
    }

    /// Prepare the SAME consistent reference inertia through the existing
    /// sparse symmetric direct owner. Exact element mass moments are assembled
    /// into one scalar sparse block per distinct Cartesian support pattern;
    /// no dense free-node mass matrix or 512-node dense leaf is required.
    ///
    /// Reference row sums scale the solver coordinates. The caller chooses the
    /// existing `fs_sparse::direct::DirectOrdering`; nested dissection receives
    /// original reference node coordinates. Each distinct pattern is analyzed
    /// and factored once. Numeric factors and prepared solve permutations are
    /// immutable and shared by equal support patterns and cloned checkpoints.
    /// Warm inverse actions allocate nothing.
    ///
    /// Node and assembly caps are checked before sparse construction. Retained
    /// symbolic bytes, factor fill and estimated work are admitted after the
    /// original synchronous analysis; numeric payload is gated by its original
    /// factor owner. These are payload/work limits, not an RSS or symbolic
    /// construction peak-memory guarantee. Cancellation is checked around those
    /// synchronous leaves and while assembling; it cannot interrupt a leaf.
    /// Accepted physical momentum, support impulses and energy use the original
    /// consistent matrix/inverse and the same nonlinear Verlet map.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_sparse_consistent_inertia(
        patch: &'a MembranePatch,
        prescribed: &[(usize, f64)],
        displacement_m: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        time_s: f64,
        budget: SparseConsistentMassBudget,
        ordering: fs_sparse::direct::DirectOrdering,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MembraneError> {
        Self::prepare(
            patch,
            prescribed,
            displacement_m,
            velocity_m_s,
            time_s,
            Some(ConsistentPreparation::Sparse(budget, ordering)),
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn prepare(
        patch: &'a MembranePatch,
        prescribed: &[(usize, f64)],
        displacement_m: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        time_s: f64,
        consistent: Option<ConsistentPreparation>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MembraneError> {
        poll(&mut cancelled)?;
        let nodes = patch.mesh.nodes.len();
        if !finite_vectors(displacement_m, nodes)
            || !finite_vectors(velocity_m_s, nodes)
            || !time_s.is_finite()
        {
            return Err(MembraneError::Invalid(
                "finite membrane motion state required",
            ));
        }
        let mut fixed = vec![false; 3 * nodes];
        for &(dof, value) in prescribed {
            poll(&mut cancelled)?;
            if dof >= fixed.len() || fixed[dof] || !value.is_finite() {
                return Err(MembraneError::Invalid(
                    "unique finite in-range stationary supports required",
                ));
            }
            if displacement_m[dof / 3][dof % 3] != value || velocity_m_s[dof / 3][dof % 3] != 0.0 {
                return Err(MembraneError::Invalid(
                    "initial motion must satisfy stationary supports exactly",
                ));
            }
            fixed[dof] = true;
        }
        let free: Vec<_> = if consistent.is_some() {
            // Axis-major coordinates permit three direct scalar-block solves
            // without gathering strided values or forming a dense 3N matrix.
            let mut free = Vec::with_capacity(fixed.iter().filter(|&&held| !held).count());
            for axis in 0..3 {
                for node in 0..nodes {
                    if matches!(consistent, Some(ConsistentPreparation::Sparse(..))) {
                        poll(&mut cancelled)?;
                    }
                    let dof = 3 * node + axis;
                    if !fixed[dof] {
                        free.push(dof);
                    }
                }
            }
            free
        } else {
            (0..fixed.len()).filter(|&dof| !fixed[dof]).collect()
        };
        if free.is_empty() {
            return Err(MembraneError::Invalid(
                "motion requires a free displacement coordinate",
            ));
        }
        let scratch_size = free
            .len()
            .checked_mul(4)
            .ok_or(MembraneError::Budget("motion workspace overflow"))?;
        let inverse_mass = if let Some(budget) = consistent {
            InverseMass::consistent(patch, &free, &fixed, budget, &mut cancelled)?
        } else {
            InverseMass::Lumped(
                free.iter()
                    .map(|&dof| 1.0 / patch.nodal_mass_kg[dof / 3])
                    .collect(),
            )
        };
        let q: Vec<_> = free
            .iter()
            .map(|&dof| displacement_m[dof / 3][dof % 3])
            .collect();
        let (momentum, mut trial_velocity): (Vec<_>, Vec<_>) = if consistent.is_some() {
            let full = patch.reference_mass_apply(
                MembraneMassImage::Consistent,
                velocity_m_s,
                &mut cancelled,
            )?;
            (
                free.iter().map(|&dof| full[dof / 3][dof % 3]).collect(),
                free.iter()
                    .map(|&dof| velocity_m_s[dof / 3][dof % 3])
                    .collect(),
            )
        } else {
            (
                free.iter()
                    .map(|&dof| patch.nodal_mass_kg[dof / 3] * velocity_m_s[dof / 3][dof % 3])
                    .collect(),
                Vec::new(),
            )
        };
        if consistent.is_some() {
            inverse_mass.apply(&momentum, &mut trial_velocity);
        }
        let kinetic_energy_j = kinetic(&momentum, &inverse_mass, &trial_velocity);
        if matches!(&inverse_mass, InverseMass::Lumped(diagonal) if diagonal.iter().any(|x| !x.is_finite() || *x <= 0.0))
            || momentum.iter().any(|x| !x.is_finite())
            || !kinetic_energy_j.is_finite()
            || kinetic_energy_j < 0.0
            || trial_velocity.iter().any(|v| !v.is_finite())
        {
            return Err(MembraneError::Invalid(
                "unrepresentable reference inertia or momentum",
            ));
        }
        let mut observed_velocity = velocity_m_s.to_vec();
        if consistent.is_some() {
            for (&dof, &velocity) in free.iter().zip(&trial_velocity) {
                observed_velocity[dof / 3][dof % 3] = velocity;
            }
        }
        let mut gradient = vec![[0.0; 3]; nodes];
        let strain_energy_j = patch.evaluate_into(displacement_m, &mut gradient, &mut cancelled)?;
        if !(strain_energy_j + kinetic_energy_j).is_finite() {
            return Err(MembraneError::Invalid(
                "unrepresentable initial mechanical energy",
            ));
        }
        poll(&mut cancelled)?;
        Ok(Self {
            patch,
            fixed,
            free,
            inverse_mass,
            displacement_m: displacement_m.to_vec(),
            velocity_m_s: observed_velocity,
            trial_q: q.clone(),
            trial_momentum: momentum.clone(),
            trial_velocity,
            q,
            momentum,
            time_s,
            strain_energy_j,
            kinetic_energy_j,
            trial_displacement: displacement_m.to_vec(),
            gradient,
            trial_gradient: vec![[0.0; 3]; nodes],
            pressure_scratch: vec![[0.0; 3]; nodes],
            scratch: vec![0.0; scratch_size],
        })
    }

    /// Original immutable material, reference mesh and reference mass owner.
    pub const fn patch(&self) -> &MembranePatch {
        self.patch
    }
    /// Immutable reference mass image used by this motion and its observations.
    pub fn inertia_image(&self) -> MembraneMassImage {
        self.inverse_mass.image()
    }
    /// Sparse scalar storage and warm-work accounting, when that inverse was
    /// explicitly selected. Shared factors are counted once in stored payload
    /// and once per active component in the inverse-action work count.
    pub fn sparse_inertia_statistics(&self) -> Option<SparseMassStatistics> {
        self.inverse_mass.sparse_statistics()
    }
    /// Accepted full-space displacement [m], including exact fixed values.
    pub fn displacement_m(&self) -> &[[f64; 3]] {
        &self.displacement_m
    }
    /// Accepted full-space velocity [m/s]; constrained components remain zero.
    pub fn velocity_m_s(&self) -> &[[f64; 3]] {
        &self.velocity_m_s
    }
    /// Accepted time [s].
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    /// Original material stored energy [J].
    pub const fn strain_energy_j(&self) -> f64 {
        self.strain_energy_j
    }
    /// Kinetic energy from reference inertia [J].
    pub const fn kinetic_energy_j(&self) -> f64 {
        self.kinetic_energy_j
    }

    /// Advance through the original `fs-time` kick--drift--kick map. The dead
    /// load is held constant for this step, including both endpoint kicks.
    /// Loads may change between calls; their work is accounted per step. A
    /// finite nonzero `dt_s` is required; negative time reverses the same map.
    ///
    /// Assembly and integration reuse prepared O(nodes) buffers and apply
    /// element forces directly without constructing a global stiffness matrix.
    /// The immutable patch permits reuse of the accepted initial force; only
    /// the drifted endpoint requires a new constitutive evaluation each step.
    /// Original constitutive owners retain their own evaluation costs. No
    /// automatic stability limit, exact-energy or real-time throughput claim
    /// is inferred from a successful step; refine time against the returned
    /// signed energy defect and the physical observable of interest.
    pub fn step(
        &mut self,
        dt_s: f64,
        dead_load_n: &[[f64; 3]],
        cancelled: impl FnMut() -> bool,
    ) -> Result<MembraneStep, MembraneError> {
        self.step_with_pressure(dt_s, dead_load_n, 0.0, cancelled)
    }

    /// Advance the same original momentum map with held material pressure [Pa].
    /// Positive pressure acts along the current reference-oriented normal,
    /// initially +z. It is a pressure difference, not an absolute gas pressure.
    ///
    /// Both endpoint kicks assemble the original current-area traction,
    /// including loads at prescribed nodes in the support impulse. Work is
    /// integrated exactly along the resolved straight drift, retaining changing
    /// face areas and increments below the absolute-position ULP. Samples in
    /// `MembranePressure::TriangleLinear` are held in each original triangle's
    /// material coordinates; a spatial fluid field must be resampled by its owner.
    /// Nonuniform profiles or moving perimeter nodes can give path-dependent
    /// work. The positional-force Verlet map remains second order and reversible,
    /// but is not claimed symplectic in those cases. Only uniform pressure on a
    /// fixed perimeter recovers the native conservative pressure-volume load.
    ///
    /// Pressure may change between calls. No gas inventory, acoustic stiffness,
    /// fluid inertia, or pressure-dependent material state is inferred. Pressure
    /// assembly and work reuse prepared O(nodes) storage; a rejected step leaves
    /// the complete accepted mechanical state unchanged. A scalar zero pressure
    /// follows the original unloaded arithmetic path.
    pub fn step_with_pressure<'p>(
        &mut self,
        dt_s: f64,
        dead_load_n: &[[f64; 3]],
        pressure: impl Into<MembranePressure<'p>>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MembraneStep, MembraneError> {
        poll(&mut cancelled)?;
        let pressure_load = pressure.into();
        pressure_load.validate(self.patch.elements.len(), &mut cancelled)?;
        let has_pressure = !pressure_load.is_zero();
        let next_time = self.time_s + dt_s;
        if !dt_s.is_finite()
            || dt_s == 0.0
            || !next_time.is_finite()
            || next_time == self.time_s
            || !finite_vectors(dead_load_n, self.displacement_m.len())
        {
            return Err(MembraneError::Invalid(
                "finite nonzero representable time advance, pressure and nodal dead loads required",
            ));
        }
        self.trial_q.copy_from_slice(&self.q);
        self.trial_momentum.copy_from_slice(&self.momentum);
        self.trial_displacement
            .copy_from_slice(&self.displacement_m);
        let control = RefCell::new(&mut cancelled);
        let stopped = RefCell::new(None);
        let energy = Cell::new(self.strain_energy_j);
        let support_sum = Cell::new([0.0; 3]);
        let first_kick = Cell::new(true);
        let buffers = RefCell::new((
            &mut self.trial_displacement,
            &mut self.trial_gradient,
            &mut self.pressure_scratch,
        ));
        let force = |q: &[f64], output: &mut [f64]| {
            output.fill(f64::NAN);
            if stopped.borrow().is_some() {
                return;
            }
            let mut buffers = buffers.borrow_mut();
            let (u, trial_gradient, pressure) = &mut *buffers;
            for (&dof, &value) in self.free.iter().zip(q) {
                u[dof / 3][dof % 3] = value;
            }
            let gradient = if first_kick.replace(false) {
                &self.gradient
            } else {
                match self
                    .patch
                    .evaluate_into(u, trial_gradient, &mut || (control.borrow_mut())())
                {
                    Ok(value) => {
                        energy.set(value);
                    }
                    Err(error) => {
                        *stopped.borrow_mut() = Some(error);
                        return;
                    }
                }
                trial_gradient
            };
            if has_pressure {
                if let Err(error) =
                    self.patch
                        .pressure_load_into_n(u, pressure_load, pressure, &mut || {
                            (control.borrow_mut())()
                        })
                {
                    *stopped.borrow_mut() = Some(error);
                    return;
                }
            }
            for (out, &dof) in output.iter_mut().zip(&self.free) {
                let mut applied = dead_load_n[dof / 3][dof % 3];
                if has_pressure {
                    applied += pressure[dof / 3][dof % 3];
                }
                *out = applied - gradient[dof / 3][dof % 3];
            }
            let mut reactions = support_sum.get();
            for (dof, &held) in self.fixed.iter().enumerate() {
                if held {
                    let mut applied = dead_load_n[dof / 3][dof % 3];
                    if has_pressure {
                        applied += pressure[dof / 3][dof % 3];
                    }
                    reactions[dof % 3] += gradient[dof / 3][dof % 3] - applied;
                }
            }
            support_sum.set(reactions);
        };
        let result = verlet_momentum_step_with_displacement(
            &mut self.trial_q,
            &mut self.trial_momentum,
            dt_s,
            &self.inverse_mass,
            &force,
            &mut self.scratch,
            &mut || (control.borrow_mut())(),
        );
        if let Some(error) = stopped.into_inner() {
            return Err(error);
        }
        let displacement = result.map_err(|error| match error {
            MomentumVerletError::Cancelled => MembraneError::Cancelled,
            MomentumVerletError::InvalidInput(what) | MomentumVerletError::NonFinite(what) => {
                MembraneError::Invalid(what)
            }
        })?;
        if matches!(self.inverse_mass, InverseMass::Consistent { .. }) {
            self.inverse_mass
                .apply(&self.trial_momentum, &mut self.trial_velocity);
        }
        let kinetic_energy_j = kinetic(
            &self.trial_momentum,
            &self.inverse_mass,
            &self.trial_velocity,
        );
        let strain_energy_j = energy.get();
        let dead_load_work_j: f64 = self
            .free
            .iter()
            .zip(displacement)
            .map(|(&dof, &increment)| dead_load_n[dof / 3][dof % 3] * increment)
            .sum();
        let pressure_work_j = if !has_pressure {
            0.0
        } else {
            self.pressure_scratch.fill([0.0; 3]);
            for (&dof, &increment) in self.free.iter().zip(displacement) {
                self.pressure_scratch[dof / 3][dof % 3] = increment;
            }
            self.patch.pressure_work_from_increment_j(
                &self.displacement_m,
                &self.pressure_scratch,
                pressure_load,
                &mut cancelled,
            )?
        };
        let energy_defect_j = kinetic_change(
            &self.momentum,
            &self.trial_momentum,
            &self.inverse_mass,
            &self.free,
            &self.velocity_m_s,
            &self.trial_velocity,
        ) + (strain_energy_j - self.strain_energy_j)
            - dead_load_work_j
            - pressure_work_j;
        let mut support_impulse_n_s = support_sum.get().map(|f| (0.5 * dt_s) * f);
        if self.inverse_mass.has_support_coupling() {
            // The drift observation has already supplied the held-load work,
            // so its integration workspace is now reusable. Recover delta_v
            // from accepted delta_p rather than subtracting large velocities.
            let n = self.free.len();
            let (delta_p, work) = self.scratch.split_at_mut(n);
            for ((delta, new), old) in delta_p
                .iter_mut()
                .zip(&self.trial_momentum)
                .zip(&self.momentum)
            {
                *delta = new - old;
            }
            self.inverse_mass.apply(delta_p, &mut work[..n]);
            let inertia = self.inverse_mass.support_impulse(&work[..n]);
            for axis in 0..3 {
                support_impulse_n_s[axis] += inertia[axis];
            }
        }
        let invalid_velocity = match &self.inverse_mass {
            InverseMass::Lumped(diagonal) => self
                .trial_momentum
                .iter()
                .zip(diagonal)
                .any(|(p, inv)| !(p * inv).is_finite()),
            InverseMass::Consistent { .. } => self.trial_velocity.iter().any(|v| !v.is_finite()),
        };
        if [
            kinetic_energy_j,
            strain_energy_j,
            kinetic_energy_j + strain_energy_j,
            dead_load_work_j,
            pressure_work_j,
            energy_defect_j,
        ]
        .iter()
        .any(|x| !x.is_finite())
            || support_impulse_n_s.iter().any(|x| !x.is_finite())
            || invalid_velocity
            || kinetic_energy_j < 0.0
        {
            return Err(MembraneError::Invalid(
                "unrepresentable membrane motion or work balance",
            ));
        }
        poll(&mut cancelled)?;
        // Publication follows the final poll: none of these infallible copies
        // can leave half an accepted state after a returned refusal.
        self.q.copy_from_slice(&self.trial_q);
        self.momentum.copy_from_slice(&self.trial_momentum);
        self.displacement_m
            .copy_from_slice(&self.trial_displacement);
        self.gradient.copy_from_slice(&self.trial_gradient);
        match &self.inverse_mass {
            InverseMass::Lumped(diagonal) => {
                for (i, &dof) in self.free.iter().enumerate() {
                    self.velocity_m_s[dof / 3][dof % 3] = self.momentum[i] * diagonal[i];
                }
            }
            InverseMass::Consistent { .. } => {
                for (&dof, &velocity) in self.free.iter().zip(&self.trial_velocity) {
                    self.velocity_m_s[dof / 3][dof % 3] = velocity;
                }
            }
        }
        self.time_s = next_time;
        self.strain_energy_j = strain_energy_j;
        self.kinetic_energy_j = kinetic_energy_j;
        Ok(MembraneStep {
            time_s: next_time,
            strain_energy_j,
            kinetic_energy_j,
            dead_load_work_j,
            pressure_work_j,
            energy_defect_j,
            support_impulse_n_s,
        })
    }
}

/// Finite-amplitude physical loss composed with the original nonlinear motion.
pub mod viscosity;
