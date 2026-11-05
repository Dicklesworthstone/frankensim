//! Nonlinear elastic sheet motion using the shared physical-momentum Verlet map.
//!
//! The original membrane energy supplies BOTH force kicks, including geometric
//! stiffening and every admitted nonlinear material law. Only free coordinates
//! are integrated; stationary supports and the reference lumped mass are fixed.
//! No equilibrium tangent, authored spring, numerical damping or mass rescaling
//! replaces the material. Work arrays are retained across steps.
//!
//! This explicit, second-order lane needs a resolved time step. Its signed
//! energy/work defect is integration error, not dissipated heat. Symplecticity
//! does not establish nonlinear stability or mesh/time convergence. Moving
//! supports, follower pressure, contact, and material history are not included.

use super::{MembraneError, MembranePatch, finite_vectors, poll};
use fs_solver::op::LinearOp;
use fs_time::{MomentumVerletError, verlet_momentum_step_with_displacement};
use std::cell::{Cell, RefCell};

#[derive(Clone, Debug)]
struct InverseMass(Vec<f64>);

impl LinearOp for InverseMass {
    fn n(&self) -> usize {
        self.0.len()
    }

    fn apply(&self, momentum: &[f64], velocity: &mut [f64]) {
        for ((out, p), inverse) in velocity.iter_mut().zip(momentum).zip(&self.0) {
            *out = p * inverse;
        }
    }
}

fn kinetic(momentum: &[f64], inverse: &InverseMass) -> f64 {
    momentum
        .iter()
        .zip(&inverse.0)
        .map(|(p, inv)| 0.5 * p * (inv * p))
        .sum()
}

// Polarize the original quadratic form: do not subtract large absolute
// kinetic energies to recover a small, resolved mechanical transfer.
fn kinetic_change(before: &[f64], after: &[f64], inverse: &InverseMass) -> f64 {
    before
        .iter()
        .zip(after)
        .zip(&inverse.0)
        .map(|((old, new), inv)| (new - old) * (inv * (0.5 * old + 0.5 * new)))
        .sum()
}

/// Physical endpoint and signed work balance of one accepted step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MembraneStep {
    /// Accepted endpoint time [s]. Negative steps can replay conservative motion.
    pub time_s: f64,
    /// Original material strain energy [J], including any static prestretch.
    pub strain_energy_j: f64,
    /// Kinetic energy from the original reference lumped mass [J].
    pub kinetic_energy_j: f64,
    /// Work of the held load over the resolved kinematic drift [J], even when
    /// adding that displacement to the absolute position rounds it away.
    pub dead_load_work_j: f64,
    /// Change in total mechanical energy minus external work [J]. This is a
    /// signed numerical defect; no physical dissipation is introduced.
    pub energy_defect_j: f64,
    /// Net support impulse [N s], using the same two endpoint force kicks.
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
    trial_displacement: Vec<[f64; 3]>,
    gradient: Vec<[f64; 3]>,
    trial_gradient: Vec<[f64; 3]>,
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
        let free: Vec<_> = (0..fixed.len()).filter(|&dof| !fixed[dof]).collect();
        if free.is_empty() {
            return Err(MembraneError::Invalid(
                "motion requires a free displacement coordinate",
            ));
        }
        let scratch_size = free
            .len()
            .checked_mul(4)
            .ok_or(MembraneError::Budget("motion workspace overflow"))?;
        let inverse_mass = InverseMass(
            free.iter()
                .map(|&dof| 1.0 / patch.nodal_mass_kg[dof / 3])
                .collect(),
        );
        let q: Vec<_> = free
            .iter()
            .map(|&dof| displacement_m[dof / 3][dof % 3])
            .collect();
        let momentum: Vec<_> = free
            .iter()
            .map(|&dof| patch.nodal_mass_kg[dof / 3] * velocity_m_s[dof / 3][dof % 3])
            .collect();
        let kinetic_energy_j = kinetic(&momentum, &inverse_mass);
        if inverse_mass.0.iter().any(|x| !x.is_finite() || *x <= 0.0)
            || momentum.iter().any(|x| !x.is_finite())
            || !kinetic_energy_j.is_finite()
        {
            return Err(MembraneError::Invalid(
                "unrepresentable reference inertia or momentum",
            ));
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
            velocity_m_s: velocity_m_s.to_vec(),
            trial_q: q.clone(),
            trial_momentum: momentum.clone(),
            q,
            momentum,
            time_s,
            strain_energy_j,
            kinetic_energy_j,
            trial_displacement: displacement_m.to_vec(),
            gradient,
            trial_gradient: vec![[0.0; 3]; nodes],
            scratch: vec![0.0; scratch_size],
        })
    }

    /// Original immutable material, reference mesh and reference mass owner.
    pub const fn patch(&self) -> &MembranePatch {
        self.patch
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
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MembraneStep, MembraneError> {
        poll(&mut cancelled)?;
        let next_time = self.time_s + dt_s;
        if !dt_s.is_finite()
            || dt_s == 0.0
            || !next_time.is_finite()
            || next_time == self.time_s
            || !finite_vectors(dead_load_n, self.displacement_m.len())
        {
            return Err(MembraneError::Invalid(
                "finite nonzero representable time advance and finite nodal dead loads required",
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
        let buffers = RefCell::new((&mut self.trial_displacement, &mut self.trial_gradient));
        let force = |q: &[f64], output: &mut [f64]| {
            output.fill(f64::NAN);
            if stopped.borrow().is_some() {
                return;
            }
            let mut buffers = buffers.borrow_mut();
            let (u, trial_gradient) = &mut *buffers;
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
            for (out, &dof) in output.iter_mut().zip(&self.free) {
                *out = dead_load_n[dof / 3][dof % 3] - gradient[dof / 3][dof % 3];
            }
            let mut reactions = support_sum.get();
            for (dof, &held) in self.fixed.iter().enumerate() {
                if held {
                    reactions[dof % 3] +=
                        gradient[dof / 3][dof % 3] - dead_load_n[dof / 3][dof % 3];
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
        let kinetic_energy_j = kinetic(&self.trial_momentum, &self.inverse_mass);
        let strain_energy_j = energy.get();
        let dead_load_work_j: f64 = self
            .free
            .iter()
            .zip(displacement)
            .map(|(&dof, &increment)| dead_load_n[dof / 3][dof % 3] * increment)
            .sum();
        let energy_defect_j = kinetic_change(&self.momentum, &self.trial_momentum, &self.inverse_mass)
            + (strain_energy_j - self.strain_energy_j)
            - dead_load_work_j;
        let support_impulse_n_s = support_sum.get().map(|f| (0.5 * dt_s) * f);
        if [
            kinetic_energy_j,
            strain_energy_j,
            kinetic_energy_j + strain_energy_j,
            dead_load_work_j,
            energy_defect_j,
        ]
        .iter()
        .any(|x| !x.is_finite())
            || support_impulse_n_s.iter().any(|x| !x.is_finite())
            || self
                .trial_momentum
                .iter()
                .zip(&self.inverse_mass.0)
                .any(|(p, inv)| !(p * inv).is_finite())
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
        for (i, &dof) in self.free.iter().enumerate() {
            self.velocity_m_s[dof / 3][dof % 3] = self.momentum[i] * self.inverse_mass.0[i];
        }
        self.time_s = next_time;
        self.strain_energy_j = strain_energy_j;
        self.kinetic_energy_j = kinetic_energy_j;
        Ok(MembraneStep {
            time_s: next_time,
            strain_energy_j,
            kinetic_energy_j,
            dead_load_work_j,
            energy_defect_j,
            support_impulse_n_s,
        })
    }
}
