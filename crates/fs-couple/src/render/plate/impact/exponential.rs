//! Exact coupled linear propagation inside the original string/contact model.
//!
//! The retained EI/tension modes, reciprocal volume springs and full viscous
//! matrix form one quadratic block. fs-phs::PreparedNonlinearZoh propagates that
//! block without midpoint frequency warping and solves the original stretching
//! and elastic-contact remainder implicitly. This is a different numerical
//! discretization of the SAME Hamiltonian, not a new string law or fitted filter.
//!
//! Nonlinear/contact accuracy still requires time and modal refinement. Only
//! the quadratic held-input limit is exact. Nonlinear drag,
//! shell/support potentials and nonlinear mapped-strain memory are not admitted
//! by this entry point. Cold linear Maxwell and passive radiation attachments
//! retain their complete original storage and interconnection. The felt entry
//! retains the original WoolFelt potential,
//! accepted compaction history and complete Kelvin memory in this same solve.
//! Nonlinear/contact transients and conditioning still need time refinement.

use super::{
    BodyPotential, ImpactBody, ImpactConfig, ImpactError, ImpactFrame, ImpactSystem,
    MAX_IMPACT_MODES, MechanicalStorage, VolumeSpring, audio::ImpactSource, damping::ViscousDamper,
    felt::FeltPad, invalid, string::StringPotential,
};
use fs_dcontact::Obstacle;

mod memory;
use fs_exec::CancelGate;
use fs_material::fiber::WoolFeltState;
use fs_math::dd::Dd;
use fs_phs::linear_exact::{LinearSplitBudget, PreparedNonlinearZoh, QuadraticPart};
use fs_phs::{PreparedStepError, QuadraticStorage, Storage};
use std::{ops::Deref, rc::Rc};

fn owner(error: fs_phs::PhsError) -> ImpactError {
    ImpactError::Owner(error.to_string())
}
fn step_error(error: PreparedStepError) -> ImpactError {
    match error {
        PreparedStepError::Cancelled => ImpactError::Cancelled,
        PreparedStepError::Solver(error) => ImpactError::PreparedSolve(error),
    }
}

struct StringRemainder {
    start: usize,
    source: StringPotential,
    quadratic: QuadraticStorage,
}

// This wrapper declares a numerical split, while H and grad H continue to call
// the original mechanical owner. It exists only in the opt-in preparation.
pub(super) struct MechanicalSplit {
    source: Rc<MechanicalStorage>,
    quadratic: QuadraticStorage,
    strings: Vec<StringRemainder>,
}
impl MechanicalSplit {
    pub(super) fn new(source: Rc<MechanicalStorage>) -> Result<Self, ImpactError> {
        let dim = 2 * source.modes
            + source
                .pads
                .iter()
                .map(|pad| pad.spec.creep.len())
                .sum::<usize>();
        let mut matrix = vec![0.0; dim * dim];
        for mode in 0..source.modes {
            matrix[(2 * mode + 1) * dim + 2 * mode + 1] = 1.0;
        }
        // Original Kelvin coordinates are sqrt(k)*deformation and store z²/2.
        // Their original k/c resistance remains in the SAME complete operator.
        for coordinate in 2 * source.modes..dim {
            matrix[coordinate * dim + coordinate] = 1.0;
        }
        let mut strings = Vec::new();
        let mut start = 0;
        for body in &source.bodies {
            let (omega, nonlinear) = match body {
                BodyPotential::Linear(omega) => (omega.as_slice(), None),
                BodyPotential::String(string) => (string.omegas(), Some(string)),
                _ => {
                    return Err(invalid(
                        "exponential impact admits original fixed-end strings and linear bodies only",
                    ));
                }
            };
            let count = omega.len();
            for (mode, frequency) in omega.iter().enumerate() {
                let coordinate = 2 * (start + mode);
                matrix[coordinate * dim + coordinate] = frequency * frequency;
            }
            if let Some(string) = nonlinear.filter(|string| !string.is_quadratic()) {
                let mut local = vec![0.0; count * count];
                for (mode, frequency) in omega.iter().enumerate() {
                    local[mode * count + mode] = frequency * frequency;
                }
                strings.push(StringRemainder {
                    start,
                    source: string.clone(),
                    quadratic: QuadraticStorage::new(local, count).map_err(owner)?,
                });
            }
            start += count;
        }
        // Preserve every cross term of each original signed volume coordinate.
        for volume in &source.volumes {
            let stiffness = volume.bulk_modulus_pa / volume.volume_m3;
            for i in 0..source.modes {
                for j in 0..=i {
                    let row = 2 * i;
                    let column = 2 * j;
                    matrix[row * dim + column] += stiffness * volume.areas[i] * volume.areas[j];
                    matrix[column * dim + row] = matrix[row * dim + column];
                }
            }
        }
        Ok(Self {
            source,
            strings,
            quadratic: QuadraticStorage::new(matrix, dim).map_err(owner)?,
        })
    }
}
impl Storage for MechanicalSplit {
    fn hamiltonian(&self, state: &[f64]) -> f64 {
        self.source.hamiltonian(state)
    }
    fn gradient(&self, state: &[f64], out: &mut [f64]) {
        self.source.gradient(state, out);
    }
    fn quadratic_part(&self) -> Option<QuadraticPart> {
        Some(QuadraticPart {
            matrix: self.quadratic.matrix().to_vec(),
            exact: self.strings.is_empty() && self.source.pads.is_empty(),
        })
    }
    fn midpoint_energy_remainder(&self, a: &[f64], b: &[f64], mid: &[f64]) -> Option<f64> {
        let dim = 2 * self.source.modes
            + self
                .source
                .pads
                .iter()
                .map(|pad| pad.spec.creep.len())
                .sum::<usize>();
        if a.len() != dim || b.len() != dim || mid.len() != dim {
            return None;
        }
        let mut result = Dd::from_f64(self.quadratic.midpoint_energy_remainder(a, b, mid)?);
        let (mut qa, mut qb, mut qm) = (
            [0.0; MAX_IMPACT_MODES],
            [0.0; MAX_IMPACT_MODES],
            [0.0; MAX_IMPACT_MODES],
        );
        for string in &self.strings {
            let count = string.source.mode_count();
            for mode in 0..count {
                let coordinate = 2 * (string.start + mode);
                qa[mode] = a[coordinate];
                qb[mode] = b[coordinate];
                qm[mode] = mid[coordinate];
            }
            // Only tiny stable remainders are subtracted. No difference of
            // large total energies is used to recover weak stretching work.
            let complete = string.source.midpoint_energy_remainder(
                &qa[..count],
                &qb[..count],
                &qm[..count],
            )?;
            let linear = string.quadratic.midpoint_energy_remainder(
                &qa[..count],
                &qb[..count],
                &qm[..count],
            )?;
            result = result + Dd::from_f64(complete) - Dd::from_f64(linear);
        }
        let histories = self.source.histories.borrow();
        for (pad, history) in self.source.pads.iter().zip(histories.iter()) {
            result = result
                + Dd::from_f64(pad.midpoint_remainder(a, b, mid, self.source.modes, history));
        }
        Some(result.to_f64())
    }
}

/// The original impact state and acceptance gates with an exponential linear block.
///
/// Original state order, strings, signed contacts, volume energies and full
/// damping survive preparation. The dense owner's complete-state capacity is
/// enforced before assembly; no modes are discarded to fit it. Accepted state
/// is read-only through Deref, so a coarser exponential clock cannot accidentally
/// advance through the midpoint owner. Cold preparation allocates; warm storage
/// callbacks and the existing exponential owner reuse bounded scratch.
pub struct ExponentialImpactSystem {
    inner: ImpactSystem,
    flow: PreparedNonlinearZoh,
    candidate: Vec<f64>,
    output: Vec<f64>,
    histories: Vec<WoolFeltState>,
}
impl ExponentialImpactSystem {
    /// Prepare held-force evolution of original linear/string bodies with
    /// elastic contacts and reciprocal fluid-volume springs.
    ///
    /// Every body frequency and the complete coupled quadratic spectrum must
    /// fit the original 0.9*pi sampled-band guard. The midpoint one-cent bound
    /// is inapplicable to the exact linear block. Nonlinear harmonics and contact
    /// transients still require a sufficiently refined physical clock.
    pub fn new(
        bodies: Vec<ImpactBody>,
        contacts: Vec<Obstacle>,
        volumes: Vec<VolumeSpring>,
        config: ImpactConfig,
        gate: &CancelGate,
    ) -> Result<Self, ImpactError> {
        Self::new_with_dampers(bodies, contacts, volumes, Vec::new(), config, gate)
    }

    /// Retain complete spatial damping, including reciprocal cross-body terms.
    /// No damping ratio, diagonal approximation or lagged drag is substituted.
    pub fn new_with_dampers(
        bodies: Vec<ImpactBody>,
        contacts: Vec<Obstacle>,
        volumes: Vec<VolumeSpring>,
        dampers: Vec<ViscousDamper>,
        config: ImpactConfig,
        gate: &CancelGate,
    ) -> Result<Self, ImpactError> {
        Self::new_with_felt(bodies, contacts, Vec::new(), volumes, dampers, config, gate)
    }

    /// Prepare the original nonlinear felt and all of its Kelvin arms together
    /// with linear/string bodies, reciprocal volumes and spatial damping.
    ///
    /// Kelvin recovery and original body modes use the existing exact linear
    /// block. The original nonlinear felt effort stays in its implicit
    /// discrete-gradient equation; irreversible conditioning commits only after
    /// the same mechanical energy and validity gates. No held-stress substep,
    /// lagged contact force or extra heat account is introduced.
    ///
    /// The complete state, including every Kelvin arm, must fit the dense
    /// owner's unchanged capacity. The material and contact time discretization
    /// still require refinement; exact free recovery does not certify impact
    /// duration, conditioning transitions or nonlinear harmonics.
    pub fn new_with_felt(
        bodies: Vec<ImpactBody>,
        contacts: Vec<Obstacle>,
        pads: Vec<FeltPad>,
        volumes: Vec<VolumeSpring>,
        dampers: Vec<ViscousDamper>,
        config: ImpactConfig,
        gate: &CancelGate,
    ) -> Result<Self, ImpactError> {
        if gate.is_requested() {
            return Err(ImpactError::Cancelled);
        }
        let modes = bodies
            .iter()
            .try_fold(0usize, |n, b| n.checked_add(b.potential.count()))
            .ok_or_else(|| invalid("exponential impact coordinate count overflows"))?;
        let memory = pads
            .iter()
            .try_fold(0usize, |n, pad| n.checked_add(pad.creep.len()))
            .ok_or_else(|| invalid("exponential felt memory count overflows"))?;
        let states = modes
            .checked_mul(2)
            .and_then(|n| n.checked_add(memory))
            .ok_or_else(|| invalid("exponential impact state count overflows"))?;
        if !PreparedNonlinearZoh::fits_reference_capacity(states, modes) {
            return Err(invalid(
                "complete exponential impact state exceeds the original dense owner capacity",
            ));
        }
        if bodies.iter().any(|body| {
            !matches!(
                &body.potential,
                BodyPotential::Linear(_) | BodyPotential::String(_)
            )
        }) || contacts
            .iter()
            .any(|contact| contact.internal_loss() != 0.0)
        {
            return Err(invalid(
                "exponential impact requires fixed-end strings/linear bodies and elastic contacts",
            ));
        }
        let inner = ImpactSystem::assemble_with_dampers(
            bodies, contacts, pads, volumes, dampers, config, false,
        )?;
        Self::prepare_inner(inner, gate)
    }

    // Reprepare the COMPLETE original composition after cold attachment.
    // Both the native acceptance owner and exact propagator share that same
    // storage allocation; no contact-only reconstruction can omit a memory.
    fn prepare_inner(mut inner: ImpactSystem, gate: &CancelGate) -> Result<Self, ImpactError> {
        if gate.is_requested() {
            return Err(ImpactError::Cancelled);
        }
        let (states, modes) = (inner.x.len(), inner.modes);
        if states != inner.system.state_dim()
            || modes != inner.system.port_dim()
            || !PreparedNonlinearZoh::fits_reference_capacity(states, modes)
        {
            return Err(invalid(
                "complete exponential impact state exceeds its original owner dimensions or capacity",
            ));
        }
        let propagation = inner.system.share_storage();
        let flow = PreparedNonlinearZoh::new_owned(
            propagation,
            inner.config.dt_s,
            LinearSplitBudget {
                max_modes: states,
                max_omega_dt: 0.9 * core::f64::consts::PI,
            },
        )
        .map_err(owner)?;
        if gate.is_requested() {
            return Err(ImpactError::Cancelled);
        }
        let histories = inner.histories.borrow().clone();
        Ok(Self {
            inner,
            flow,
            candidate: vec![0.0; states],
            output: vec![0.0; modes],
            histories,
        })
    }

    /// Retain the owner's nonlinear solve limits without changing physical data.
    pub fn set_iteration_limit(&mut self, limit: usize) -> Result<(), ImpactError> {
        self.flow.set_max_iterations(limit).map_err(owner)
    }

    /// Advance one original held generalized-force sample transactionally.
    /// Refusal leaves accepted mechanics, energy, source data and time unchanged.
    pub fn step(
        &mut self,
        external: &[f64],
        gate: &CancelGate,
    ) -> Result<ImpactFrame, ImpactError> {
        if gate.is_requested() {
            return Err(ImpactError::Cancelled);
        }
        if self.inner.sample >= self.inner.config.max_steps {
            return Err(ImpactError::Budget);
        }
        if external.len() != self.inner.modes
            || external.iter().any(|force| {
                !force.is_finite() || force.abs() > self.inner.config.maximum_generalized_force
            })
        {
            return Err(invalid(
                "external generalized force shape or ceiling failed",
            ));
        }
        let before = self.inner.stored_energy_j();
        let record = self
            .flow
            .step_into_controlled(
                &self.inner.x,
                external,
                &mut self.candidate,
                &mut self.output,
                || gate.is_requested(),
            )
            .map_err(step_error)?;
        self.inner.accept_step(
            &mut self.candidate,
            &mut self.histories,
            before,
            record,
            gate,
        )
    }
}
impl Deref for ExponentialImpactSystem {
    type Target = ImpactSystem;
    fn deref(&self) -> &ImpactSystem {
        &self.inner
    }
}
impl ImpactSource for ExponentialImpactSystem {
    type Frame = ImpactFrame;
    fn state(&self) -> &[f64] {
        self.inner.state()
    }
    fn mode_count(&self) -> usize {
        self.inner.modes
    }
    fn samples_rendered(&self) -> u64 {
        self.inner.sample
    }
    fn sample_period_s(&self) -> f64 {
        self.inner.config.dt_s
    }
    fn remaining_steps(&self) -> u64 {
        self.inner.config.max_steps - self.inner.sample
    }
    fn maximum_generalized_force(&self) -> f64 {
        self.inner.config.maximum_generalized_force
    }
    fn advance(&mut self, forces: &[f64], gate: &CancelGate) -> Result<Self::Frame, ImpactError> {
        self.step(forces, gate)
    }
}

#[cfg(test)]
mod tests;
