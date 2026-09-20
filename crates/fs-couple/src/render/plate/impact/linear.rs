//! Prepared linear-body image of the existing impact description.
//!
//! This is a compiler into the existing modal/contact owners, NOT a time
//! integrator. Independent linear bodies use their exact held-force transition;
//! bilateral reactions use the once-factored coupling system; only the normal
//! contact reactions are nonlinear. A two-head drum therefore does not put
//! every head coordinate into the reference host's dense FD Newton solve.
//!
//! The admitted Hamiltonian, damping and contact coefficients are unchanged.
//! Finite-step coupling is not the exact exponential of the coupled continuum,
//! nor bit-identical to the Gonzalez reference. Compare by time refinement.
//! Nonlinear shells, felt, and more than two bodies in a single port are not
//! silently linearized. This image has no pressure observer or inferred loss.
//! The underlying contact path still allocates; no hard-real-time claim follows.

/// Geometry-derived tensioned filaments with reciprocal distributed contact.
pub mod wire;

use super::{BodyPotential, ImpactBody, ImpactError, VolumeSpring};
use crate::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use crate::render::schedule::force::coupled::{
    CoupledModalSystem, ModalAttachment, ModalConnection, ModalCouplingConfig,
    ModalCouplingError,
    contact::{ContactModalSystem, ModalContact, ModalContactConfig},
    contact::multiple::{MultiContactConfig, MultiContactModalSystem, MAX_NORMAL_CONTACTS},
};
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;
use fs_math::c64::C64;

/// One volume spring and its explicitly chosen reference area.
/// x = (sum areas_i*q_i)/A and k = bulk*A^2/V, so k*x^2/2 is EXACTLY
/// the original small-signal volume energy. A changes units/conditioning, not
/// the model. It is not an acoustic piston area or a material parameter.
#[derive(Debug, Clone)]
pub struct VolumeConnection {
    /// The unchanged fluid volume, bulk modulus and signed surface integrals.
    pub spring: VolumeSpring,
    /// Positive displacement-coordinate reference area [m^2].
    pub reference_area_m2: f64,
}

/// Explicit host limits plus the numerical owners' unchanged admission limits.
#[derive(Debug, Clone, Copy)]
pub struct LinearImpactConfig {
    /// Shared fixed mechanics clock [Hz], not the eventual WAV output rate.
    pub sample_rate_hz: u32,
    /// Lifetime accepted-step budget, in 1..=2^53.
    pub max_steps: u64,
    /// Magnitude ceiling on each external generalized force [N/sqrt(kg)].
    pub maximum_generalized_force: f64,
    /// Individual component state/energy limits. Pressure transfers are zero.
    pub component: ModalAcousticTimeBudget,
    /// Coupling setup, frequency, reaction and whole-system energy limits.
    pub coupling: ModalCouplingConfig,
    /// Per-point contact force, penetration and nonlinear residual limits.
    pub contact: ModalContactConfig,
    /// Joint contact-set work limits; used when more than one point is supplied.
    pub multiple: MultiContactConfig,
}

/// One accepted step, without pretending its contact residual is a pHS residual.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct LinearImpactFrame {
    /// One-based accepted mechanical sample.
    pub sample: u64,
    /// Sample endpoint [s].
    pub time_s: f64,
    /// Body, volume and contact storage [J].
    pub stored_energy_j: f64,
    /// Body/connection/contact dissipation [J].
    pub dissipated_energy_j: f64,
    /// Work from external generalized forces only [J].
    pub supplied_work_j: f64,
    /// Uncorrected discrete energy-balance residual [J].
    pub balance_residual_j: f64,
    /// Largest absolute normal-force equation residual [N]; zero without contact.
    pub maximum_contact_residual_n: f64,
    /// Scalar root iterations (one contact) or joint sweeps (several contacts).
    pub nonlinear_iterations: usize,
}

enum Prepared {
    Bilateral(CoupledModalSystem),
    Single(ContactModalSystem),
    Multiple(MultiContactModalSystem),
}
impl Prepared {
    fn components(&self) -> &[ModalAcousticTimeModel] {
        match self {
            Self::Bilateral(s) => s.components(),
            Self::Single(s) => s.components(),
            Self::Multiple(s) => s.components(),
        }
    }
    fn energy(&self) -> Result<f64, ModalCouplingError> {
        match self {
            Self::Bilateral(s) => s.total_energy_j(),
            Self::Single(s) => s.total_energy_j(),
            Self::Multiple(s) => s.total_energy_j(),
        }
    }
}

/// Stateful all-linear impact image, retaining the input body/mode order.
/// All mappings and linear coupling factors are built once. The flat state is
/// refreshed only after an owner accepts a complete step; consumers cannot
/// change a component behind the host's clock or acoustic coordinate maps.
pub struct LinearImpactSystem {
    prepared: Prepared,
    state: Vec<f64>,
    config: LinearImpactConfig,
    frame: LinearImpactFrame,
    contacts: usize,
}

fn invalid(what: &'static str) -> ImpactError { ImpactError::Invalid(what) }
fn owner(error: ModalCouplingError) -> ImpactError {
    match error {
        ModalCouplingError::Cancelled => ImpactError::Cancelled,
        other => ImpactError::Owner(other.to_string()),
    }
}

/// Express a complete signed column as the original owner's two attachments.
/// Negating the second attachment preserves the ORIGINAL column, including a
/// single body's arbitrary modal signs. There is no point-wise normalization.
fn attachments(
    column: &[f64], counts: &[usize],
) -> Result<(ModalAttachment, ModalAttachment), ImpactError> {
    if column.iter().any(|v| !v.is_finite()) {
        return Err(invalid("linear impact port contains nonfinite participation"));
    }
    let mut active = Vec::with_capacity(2);
    let mut offset = 0;
    for (component, &n) in counts.iter().enumerate() {
        let row = &column[offset..offset+n];
        if row.iter().any(|b| *b != 0.0) {
            if active.len() == 2 {
                return Err(invalid("this linear image admits at most two bodies in each volume/contact port"));
            }
            active.push(ModalAttachment { component, shapes: row.to_vec() });
        }
        offset += n;
    }
    match active.len() {
        0 => Ok((ModalAttachment { component: 0, shapes: vec![0.0; counts[0]] },
                 ModalAttachment { component: 0, shapes: vec![0.0; counts[0]] })),
        1 => {
            let left = active.remove(0);
            let right = ModalAttachment { component: left.component, shapes: vec![0.0; left.shapes.len()] };
            Ok((left, right))
        }
        _ => {
            let left = active.remove(0);
            let mut right = active.remove(0);
            for b in &mut right.shapes { *b = -*b; }
            Ok((left, right))
        }
    }
}

impl LinearImpactSystem {
    /// Compile the same physical bodies/obstacles/volumes into existing owners.
    /// No body eigenbasis is changed. Volume coordinate scaling cancels between
    /// stiffness and attachments. Distributed obstacle rows become jointly
    /// solved contact points, with original weights, gaps and constitutive data.
    ///
    /// # Errors
    /// Refuses nonlinear bodies, damping on a free coordinate, unsupported
    /// multi-body ports, malformed data and every downstream owner budget.
    /// Felt is deliberately absent from this signature: no history is discarded.
    pub fn new(
        bodies: Vec<ImpactBody>, contacts: Vec<Obstacle>, volumes: Vec<VolumeConnection>,
        config: LinearImpactConfig, gate: &CancelGate,
    ) -> Result<Self, ImpactError> {
        if gate.is_requested() { return Err(ImpactError::Cancelled); }
        let counts: Vec<_> = bodies.iter().map(|b| b.potential.count()).collect();
        let count = counts.iter().try_fold(0usize, |n, &m| n.checked_add(m))
            .ok_or_else(|| invalid("linear impact mode count overflow"))?;
        if count == 0 || counts.contains(&0) || count > config.coupling.max_modes
            || count > crate::modal_acoustic_time::MAX_TIME_DOMAIN_ACOUSTIC_MODES
            || config.sample_rate_hz == 0 || config.max_steps == 0 || config.max_steps > (1u64 << 53)
            || !config.maximum_generalized_force.is_finite() || config.maximum_generalized_force <= 0.0
            || volumes.len() > config.coupling.max_connections || contacts.len() > 32
        {
            return Err(invalid("linear impact needs bounded bodies, ports, positive clock/force and exact step budget"));
        }
        let point_count = contacts.iter().try_fold(0usize, |n, c| n.checked_add(c.n_points()))
            .ok_or_else(|| invalid("linear impact contact count overflow"))?;
        if point_count > MAX_NORMAL_CONTACTS || (point_count > 1 && point_count > config.multiple.max_contacts) {
            return Err(invalid("linear impact distributed contacts exceed the original contact budget"));
        }
        let mut models = Vec::with_capacity(bodies.len());
        for body in bodies {
            if gate.is_requested() { return Err(ImpactError::Cancelled); }
            let BodyPotential::Linear(omega) = body.potential else {
                return Err(invalid("nonlinear shell storage cannot be compiled into the linear impact image"));
            };
            if body.initial.len() != omega.len() || body.damping_per_s.len() != omega.len()
                || omega.iter().any(|w| !w.is_finite() || *w < 0.0)
                || body.damping_per_s.iter().any(|d| !d.is_finite() || *d < 0.0)
            {
                return Err(invalid("linear impact body state/damping dimensions or coefficients are invalid"));
            }
            let mut modes = Vec::with_capacity(omega.len());
            for (&w, &d) in omega.iter().zip(&body.damping_per_s) {
                if w == 0.0 && d != 0.0 {
                    return Err(invalid("free-coordinate drag needs an explicit supported image; it is not discarded"));
                }
                let zeta = if w == 0.0 { 0.0 } else { 0.5 * (d / w) };
                if d > 0.0 && zeta == 0.0 {
                    return Err(invalid("linear impact modal damping underflow"));
                }
                modes.push(ModalAcousticMode { angular_frequency_rad_s: w, damping_ratio: zeta,
                    pressure_per_modal_velocity: C64::ZERO });
            }
            let mut model = ModalAcousticTimeModel::try_new_with_free_coordinates(
                config.sample_rate_hz, modes, config.component,
            ).map_err(|e| ImpactError::Owner(e.to_string()))?;
            model.restore_states(&body.initial).map_err(|e| ImpactError::Owner(e.to_string()))?;
            models.push(model);
        }
        let mut links = Vec::with_capacity(volumes.len());
        for volume in volumes {
            let v = volume.spring;
            let area = volume.reference_area_m2;
            if v.areas.len() != count || v.areas.iter().any(|a| !a.is_finite())
                || [area, v.bulk_modulus_pa, v.volume_m3].iter().any(|v| !v.is_finite() || *v <= 0.0)
            {
                return Err(invalid("linear impact volume requires finite positive SI data and every signed modal area"));
            }
            let stiffness = (v.bulk_modulus_pa * area) * (area / v.volume_m3);
            if !stiffness.is_finite() || stiffness <= 0.0 {
                return Err(invalid("volume reference-area stiffness is not representable"));
            }
            let column: Vec<_> = v.areas.iter().map(|b| b / area).collect();
            if column.iter().zip(&v.areas).any(|(&c,&a)| a != 0.0 && c == 0.0) {
                return Err(invalid("volume reference-area participation underflow"));
            }
            let (left, right) = attachments(&column, &counts)?;
            links.push(ModalConnection { left, right, stiffness_n_m: stiffness,
                damping_n_s_m: 0.0, rest_extension_m: 0.0 });
        }
        let network = CoupledModalSystem::new(models, links, config.coupling, gate).map_err(owner)?;
        let mut points = Vec::with_capacity(point_count);
        for ob in contacts {
            // Re-admit raw-parts escape values before reading individual rows.
            let ob = Obstacle::new(ob.collocation().to_vec(), ob.n_points(), count,
                ob.gaps().to_vec(), ob.weights().to_vec(), ob.stiffness(), ob.alpha(),
                ob.provenance().to_string()).and_then(|clean| clean.with_internal_loss(ob.internal_loss()))
                .map_err(|e| ImpactError::Owner(e.to_string()))?;
            for point in 0..ob.n_points() {
                let column = &ob.collocation()[point*count..(point+1)*count];
                if column.iter().all(|b| *b == 0.0) {
                    return Err(invalid("contact requires a moving attachment in this condensed image"));
                }
                let (left, right) = attachments(column, &counts)?;
                let law = Obstacle::new(vec![-1.0], 1, 1, vec![ob.gaps()[point]],
                    vec![ob.weights()[point]], ob.stiffness(), ob.alpha(), ob.provenance().to_string())
                    .and_then(|law| law.with_internal_loss(ob.internal_loss()))
                    .map_err(|e| ImpactError::Owner(e.to_string()))?;
                points.push((ModalContact { left, right, law }, config.contact));
            }
        }
        let prepared = match points.len() {
            0 => Prepared::Bilateral(network),
            1 => {
                let (contact, limits) = points.remove(0);
                Prepared::Single(ContactModalSystem::new(network, contact, limits, gate).map_err(owner)?)
            }
            _ => Prepared::Multiple(MultiContactModalSystem::new(network, points, config.multiple, gate).map_err(owner)?),
        };
        let initial_energy = prepared.energy().map_err(owner)?;
        let mut result = Self { prepared, state: vec![0.0; 2*count], config, contacts: point_count,
            frame: LinearImpactFrame { stored_energy_j: initial_energy, ..LinearImpactFrame::default() } };
        result.refresh_state();
        if gate.is_requested() { return Err(ImpactError::Cancelled); }
        Ok(result)
    }

    fn refresh_state(&mut self) {
        let mut i = 0;
        for model in self.prepared.components() {
            for s in model.states() {
                self.state[2*i] = s.displacement_m_sqrt_kg;
                self.state[2*i+1] = s.velocity_m_sqrt_kg_per_s;
                i += 1;
            }
        }
    }
    /// Accepted interleaved q/v values in the ORIGINAL body/mode order.
    #[must_use]
    pub fn state(&self) -> &[f64] { &self.state }
    /// Number of mechanical coordinates, with no acoustic or creep state appended.
    #[must_use]
    pub fn mode_count(&self) -> usize { self.state.len()/2 }
    /// Jointly solved point count after expanding distributed obstacles.
    #[must_use]
    pub const fn contact_count(&self) -> usize { self.contacts }
    /// Actual fixed mechanical period [s].
    #[must_use]
    pub fn sample_period_s(&self) -> f64 { f64::from(self.config.sample_rate_hz).recip() }
    /// Accepted step count; refusals consume no physical time.
    #[must_use]
    pub const fn samples_rendered(&self) -> u64 { self.frame.sample }
    /// Last accepted diagnostic, or initial energy at sample zero.
    #[must_use]
    pub const fn frame(&self) -> &LinearImpactFrame { &self.frame }
    /// Absolute per-coordinate external force ceiling [N/sqrt(kg)].
    #[must_use]
    pub const fn maximum_generalized_force(&self) -> f64 { self.config.maximum_generalized_force }
    /// Remaining lifetime mechanical step allowance.
    #[must_use]
    pub fn remaining_steps(&self) -> u64 { self.config.max_steps-self.frame.sample }
    /// Extend the total budget without altering any physical/filter history.
    pub fn extend_step_budget(&mut self, total: u64) -> Result<(), ImpactError> {
        if total <= self.config.max_steps || total > (1u64 << 53) {
            return Err(invalid("linear impact budget extension must increase an exact-step horizon"));
        }
        self.config.max_steps = total;
        Ok(())
    }

    /// One joint owner step. All external force admission precedes mutation;
    /// cancellation/refusal retains the complete old state and report.
    pub fn step(&mut self, external: &[f64], gate: &CancelGate) -> Result<LinearImpactFrame, ImpactError> {
        if gate.is_requested() { return Err(ImpactError::Cancelled); }
        if self.frame.sample >= self.config.max_steps { return Err(ImpactError::Budget); }
        if external.len() != self.mode_count() || external.iter()
            .any(|f| !f.is_finite() || f.abs() > self.config.maximum_generalized_force) {
            return Err(invalid("linear impact external force exceeds its finite vector budget"));
        }
        let mut next = match &mut self.prepared {
            Prepared::Bilateral(system) => {
                let f = system.step_under_gate(external, gate).map_err(owner)?;
                LinearImpactFrame { sample: f.sample, stored_energy_j: f.modal_energy_j+f.connection_energy_j,
                    dissipated_energy_j: f.component_dissipation_j+f.connection_dissipation_j,
                    supplied_work_j: f.external_work_j, balance_residual_j: f.energy_residual_j,
                    ..LinearImpactFrame::default() }
            }
            Prepared::Single(system) => {
                let f = system.step_under_gate(external, gate).map_err(owner)?;
                LinearImpactFrame { sample: f.sample, stored_energy_j: f.network_energy_j+f.contact_energy_j,
                    dissipated_energy_j: f.network_dissipation_j+f.contact_dissipation_j,
                    supplied_work_j: f.external_work_j, balance_residual_j: f.energy_residual_j,
                    maximum_contact_residual_n: f.constitutive_residual_n.abs(), nonlinear_iterations: f.iterations,
                    ..LinearImpactFrame::default() }
            }
            Prepared::Multiple(system) => {
                let f = system.step_under_gate(external, gate).map_err(owner)?;
                LinearImpactFrame { sample: f.sample, stored_energy_j: f.network_energy_j+f.contact_energy_j,
                    dissipated_energy_j: f.network_dissipation_j+f.contact_dissipation_j+f.friction_dissipation_j,
                    supplied_work_j: f.external_work_j, balance_residual_j: f.energy_residual_j,
                    maximum_contact_residual_n: f.contacts.iter().map(|p|p.constitutive_residual_n.abs()).fold(0.0,f64::max),
                    nonlinear_iterations: f.sweeps, ..LinearImpactFrame::default() }
            }
        };
        // The owner has published. Only infallible bookkeeping follows.
        next.time_s = next.sample as f64 * self.sample_period_s();
        self.refresh_state();
        self.frame = next;
        Ok(next)
    }
}
