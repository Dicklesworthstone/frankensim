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

use super::{BodyPotential, ImpactBody, ImpactError, VolumeSpring, damping::{self, ViscousDamper}};
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
    /// Refuses nonlinear bodies, unsupported multi-body ports, malformed data
    /// and every downstream owner budget. Nonzero free-coordinate drag consumes
    /// one bilateral connection per coordinate; it is never silently dropped.
    /// Felt is deliberately absent from this signature: no history is discarded.
    pub fn new(
        bodies: Vec<ImpactBody>, contacts: Vec<Obstacle>, volumes: Vec<VolumeConnection>,
        config: LinearImpactConfig, gate: &CancelGate,
    ) -> Result<Self, ImpactError> {
        Self::new_with_dampers(bodies, contacts, volumes, Vec::new(), config, gate)
    }

    /// Compile spatial drag into the existing simultaneous bilateral port solve.
    /// A port spans at most two original body components. Diagonal drag on a
    /// zero-frequency coordinate is also lowered to a grounded viscous link.
    /// Its reaction is solved simultaneously with all volume/contact reactions,
    /// not split or evaluated on the previous velocity. This retains the free
    /// inertia without inventing a small resonance. The drag discretization is
    /// second-order and passive, not the exact damped free-body exponential.
    pub fn new_with_dampers(
        bodies: Vec<ImpactBody>, contacts: Vec<Obstacle>, volumes: Vec<VolumeConnection>,
        dampers: Vec<ViscousDamper>, config: LinearImpactConfig, gate: &CancelGate,
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
        damping::validate(&dampers, count)?;
        let active_dampers = dampers.iter().filter(|d| d.damping_n_s_m != 0.0).count();
        if volumes.len().checked_add(active_dampers).is_none_or(|n| n > config.coupling.max_connections) {
            return Err(invalid("volumes and viscous ports exceed the existing connection budget"));
        }
        let point_count = contacts.iter().try_fold(0usize, |n, c| n.checked_add(c.n_points()))
            .ok_or_else(|| invalid("linear impact contact count overflow"))?;
        if point_count > MAX_NORMAL_CONTACTS || (point_count > 1 && point_count > config.multiple.max_contacts) {
            return Err(invalid("linear impact distributed contacts exceed the original contact budget"));
        }
        let mut models = Vec::with_capacity(bodies.len());
        let mut free_drag = Vec::new();
        // The checked sum above bounds this subtraction. Derived diagonal
        // links share the caller's budget with volumes and spatial attachments.
        let remaining_links = config.coupling.max_connections - volumes.len() - active_dampers;
        for (component, body) in bodies.into_iter().enumerate() {
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
            for (mode, (&w, &d)) in omega.iter().zip(&body.damping_per_s).enumerate() {
                if w == 0.0 && d != 0.0 {
                    if free_drag.len() >= remaining_links {
                        return Err(invalid("free-coordinate drag exceeds the remaining bilateral connection budget"));
                    }
                    free_drag.push((component, mode, d));
                }
                let zeta = if w == 0.0 { 0.0 } else { 0.5 * (d / w) };
                if w > 0.0 && d > 0.0 && zeta == 0.0 {
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
        let mut links = Vec::with_capacity(volumes.len()+active_dampers+free_drag.len());
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
        for damper in dampers {
            if damper.damping_n_s_m == 0.0 { continue; }
            let (left, right) = attachments(&damper.weights, &counts)?;
            links.push(ModalConnection { left, right, stiffness_n_m: 0.0,
                damping_n_s_m: damper.damping_n_s_m, rest_extension_m: 0.0 });
        }
        for (component, mode, drag) in free_drag {
            // q,p are already mass-normalized and H_kin=p^2/2. Choosing the
            // port's REFERENCE mass as 1 kg gives b=1/sqrt(1 kg), c=drag*1 kg,
            // so c*b*b^T is exactly the supplied diagonal resistance [1/s].
            // This changes units only: no physical mass or coordinate is added.
            // The same b maps velocity and reaction, dissipating drag*p^2.
            let mut shapes = vec![0.0; counts[component]];
            shapes[mode] = 1.0;
            links.push(ModalConnection {
                left: ModalAttachment { component, shapes },
                right: ModalAttachment { component, shapes: vec![0.0; counts[component]] },
                stiffness_n_m: 0.0, damping_n_s_m: drag, rest_extension_m: 0.0,
            });
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

#[cfg(test)]
mod free_drag_tests {
    use super::*;
    use crate::modal_acoustic_time::ModalAcousticState;

    fn config(rate: u32) -> LinearImpactConfig {
        LinearImpactConfig {
            sample_rate_hz: rate, max_steps: 2_000, maximum_generalized_force: 10_000.0,
            component: ModalAcousticTimeBudget::audible_reference(),
            coupling: ModalCouplingConfig {
                max_modes: 64, max_connections: 8, max_setup_terms: 100_000,
                nyquist_guard_fraction: 0.9, maximum_total_energy_j: 10.0,
                maximum_abs_pressure_pa: 1e6, maximum_abs_connection_force_n: 1e5,
                solve_relative_tolerance: 1e-10, energy_absolute_tolerance_j: 1e-10,
                energy_relative_tolerance: 1e-8,
            },
            contact: ModalContactConfig {
                max_iterations: 100, maximum_force_n: 1e4, maximum_penetration_m: 0.01,
                force_absolute_tolerance_n: 1e-10, force_relative_tolerance: 1e-9,
            },
            multiple: MultiContactConfig { max_contacts: 32, max_sweeps: 100, max_setup_terms: 100_000 },
        }
    }

    #[test]
    fn inertial_drag_preserves_mass_force_work_and_discrete_decay() {
        let gate = CancelGate::new_clock_free(); let dt = 1.0/20_000.0;
        for drag in [0.0, 100.0, 10_000.0, 100_000.0] {
            let (mut body, tip) = ImpactBody::free_mass(0.04, 0.0, 0.03).unwrap();
            body.damping_per_s[0] = drag;
            let mut system = LinearImpactSystem::new(vec![body], vec![], vec![], config(20_000), &gate).unwrap();
            let initial = system.frame().stored_energy_j;
            let (mut q, mut p, mut work, mut loss) = (0.0, 0.006, 0.0, 0.0);
            for step in 0..100 {
                let force = if step < 20 { 0.008*tip } else { 0.0 };
                let next = ((1.0-0.5*drag*dt)*p+dt*force)/(1.0+0.5*drag*dt);
                let displacement = 0.5*dt*(p+next);
                let expected_loss = dt*drag*(0.5*(p+next)).powi(2);
                q += displacement; p = next;
                let frame = system.step(&[force], &gate).unwrap();
                work += frame.supplied_work_j; loss += frame.dissipated_energy_j;
                assert!((system.state()[0]-q).abs() < 1e-12);
                assert!((system.state()[1]-p).abs() < 1e-11);
                assert!((frame.supplied_work_j-force*displacement).abs() < 1e-12);
                assert!((frame.dissipated_energy_j-expected_loss).abs() < 1e-12);
                assert!(frame.dissipated_energy_j >= 0.0);
                assert!((frame.stored_energy_j+loss-initial-work).abs() < 1e-10);
            }
            if drag > 0.0 { assert!(loss > 0.0); }
        }
    }

    #[test]
    fn inertial_drag_refines_toward_continuous_damped_motion_without_artificial_stiffness() {
        let run = |rate: u32| {
            let gate = CancelGate::new_clock_free();
            let (mut body, _) = ImpactBody::free_mass(1.0, 0.0, 0.01).unwrap();
            body.damping_per_s[0] = 400.0;
            let mut s = LinearImpactSystem::new(vec![body], vec![], vec![], config(rate), &gate).unwrap();
            for _ in 0..rate/100 { s.step(&[0.0], &gate).unwrap(); }
            let exact_p = 0.01*(-4.0_f64).exp();
            let exact_q = (0.01-exact_p)/400.0;
            (s.state()[1]-exact_p).abs()+400.0*(s.state()[0]-exact_q).abs()
        };
        let coarse = run(20_000); let fine = run(40_000);
        assert!(coarse > 1e-10 && fine < 0.3*coarse, "{coarse:e} -> {fine:e}");
    }

    #[test]
    fn mixed_body_drag_and_contacts_match_explicit_ports_and_retry_together() {
        let parts = || vec![
            ImpactBody { potential: BodyPotential::Linear(vec![0.0,800.0]),
                damping_per_s: vec![250.0,0.0], initial: vec![
                    ModalAcousticState { displacement_m_sqrt_kg: -1e-5, velocity_m_sqrt_kg_per_s: 0.2 },
                    ModalAcousticState::default()] },
            ImpactBody { potential: BodyPotential::Linear(vec![0.0]), damping_per_s: vec![75.0],
                initial: vec![ModalAcousticState { displacement_m_sqrt_kg: -1e-5,
                    velocity_m_sqrt_kg_per_s: 0.1 }] },
        ];
        let contacts = || vec![fs_dcontact::Obstacle::new(
            vec![1.0,-1.0,0.0, 0.0,-1.0,1.0], 2, 3, vec![0.0;2], vec![1.0;2],
            1e6, 1.5, "two inertial bodies contacting the same modal receiver".into()).unwrap()];
        let gate = CancelGate::new_clock_free();
        let mut a = LinearImpactSystem::new(parts(), contacts(), vec![], config(20_000), &gate).unwrap();
        let mut explicit = parts(); explicit[0].damping_per_s[0] = 0.0; explicit[1].damping_per_s[0] = 0.0;
        let mut b = LinearImpactSystem::new_with_dampers(explicit, contacts(), vec![], vec![
            ViscousDamper { weights: vec![1.0,0.0,0.0], damping_n_s_m: 250.0 },
            ViscousDamper { weights: vec![0.0,0.0,1.0], damping_n_s_m: 75.0 },
        ], config(20_000), &gate).unwrap();
        assert_eq!(a.contact_count(), 2);
        for _ in 0..80 {
            assert_eq!(a.step(&[0.0;3], &gate).unwrap(), b.step(&[0.0;3], &gate).unwrap());
            assert_eq!(a.state(), b.state());
        }
        assert!(a.state()[3].abs() > 1e-7, "real contact must excite the shared receiver");
        let before = a.state().to_vec(); let frame = *a.frame();
        let stopped = CancelGate::new_clock_free(); stopped.request();
        assert!(matches!(a.step(&[0.0;3], &stopped), Err(ImpactError::Cancelled)));
        assert!(a.step(&[0.0,f64::NAN,0.0], &gate).is_err());
        assert_eq!(a.state(), before); assert_eq!(*a.frame(), frame);
        assert_eq!(a.step(&[0.0;3], &gate).unwrap(), b.step(&[0.0;3], &gate).unwrap());
        assert_eq!(a.state(), b.state());
    }

    #[test]
    fn derived_drag_consumes_connection_budget_without_displacing_spatial_ports() {
        let gate = CancelGate::new_clock_free();
        let (mut body, _) = ImpactBody::free_mass(1.0, 0.0, 0.01).unwrap();
        body.damping_per_s[0] = 20.0;
        let pad = ViscousDamper { weights: vec![1.0], damping_n_s_m: 10.0 };
        let mut cfg = config(20_000); cfg.coupling.max_connections = 1;
        assert!(LinearImpactSystem::new_with_dampers(vec![body.clone()], vec![], vec![],
            vec![pad.clone()], cfg, &gate).is_err());
        cfg.coupling.max_connections = 2;
        assert!(LinearImpactSystem::new_with_dampers(vec![body.clone()], vec![], vec![],
            vec![pad], cfg, &gate).is_ok());
        cfg.coupling.max_connections = 0;
        assert!(LinearImpactSystem::new(vec![body.clone()], vec![], vec![], cfg, &gate).is_err());
        body.damping_per_s[0] = 0.0;
        assert!(LinearImpactSystem::new(vec![body.clone()], vec![], vec![], cfg, &gate).is_ok());
        for drag in [-1.0, f64::NAN, f64::INFINITY] {
            body.damping_per_s[0] = drag;
            assert!(LinearImpactSystem::new(vec![body.clone()], vec![], vec![], config(20_000), &gate).is_err());
        }
    }
}
