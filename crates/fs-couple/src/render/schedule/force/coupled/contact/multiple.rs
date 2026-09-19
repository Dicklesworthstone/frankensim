//! Simultaneous compliant normal contacts on one retained modal network.
//!
//! All contacts see x = x_free - S R, with the FULL cross-contact displacement
//! compliance S_ij = b_i^T D_network b_j. A deterministic nonlinear coordinate
//! solve reuses the single-contact bracket solver. Coordinates are iterated to
//! joint convergence, not advanced as separate collisions. Only then does the
//! existing network stage one physical sample with the sum of all reactions.
//! Every law is checked again against those actual candidate endpoints.
//!
//! Nonconvergence is a refusal: no ordering-independent/global convergence or
//! rigid complementarity claim is made. Duplicated/redundant attachment maps
//! are legal when their compliant laws can be resolved within the given budget.
//! Contact laws, scalar roots and time integration retain their existing owners.
//! The contact owner allocates during evaluation. This is not a hard-real-time,
//! rigid-impact, contact-discovery or alias-free sound implementation.
//! Optional 1-D regularized Coulomb friction is attached with `with_friction`;
//! its full mixed compliance and work join the SAME sample transaction.

use super::*;

/// Tangential friction composed with these simultaneous normal contacts.
pub mod friction;
/// Two-direction set-valued Coulomb graph on the same compliant contact network.
pub mod coulomb;

/// Aggregate contact-set work limits, in addition to each contact's own limits.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MultiContactConfig {
    /// Maximum admitted contacts, in 1..=32. At least one must be supplied.
    pub max_contacts: usize,
    /// Maximum joint nonlinear sweeps per physical sample, in 1..=128.
    /// Each coordinate uses its original ModalContactConfig::max_iterations.
    pub max_sweeps: usize,
    /// Setup term budget: p*(n*(k+2)+(k+1)^2)+n*p^2, for n total modes,
    /// p contacts and k bilateral connections. Checked before shape allocation.
    pub max_setup_terms: usize,
}

/// One contact evaluated at the same accepted network endpoints as every other.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ContactPointFrame {
    /// Applied compressive reaction [N], never attractive.
    pub normal_force_n: f64,
    /// Applied reaction minus the law evaluated on actual endpoint motion [N].
    pub constitutive_residual_n: f64,
    /// The contact's original absolute-plus-relative force tolerance [N].
    pub force_tolerance_n: f64,
    /// Positive initial penetration [m].
    pub penetration_before_m: f64,
    /// Positive endpoint penetration [m].
    pub penetration_after_m: f64,
    /// Endpoint contact potential [J].
    pub stored_energy_j: f64,
    /// Nonadhesive contact loss in this sample [J].
    pub dissipation_j: f64,
}

/// Last accepted whole-network sample, including ALL contact potentials/losses.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MultiContactFrame {
    /// One-based complete sample ordinal.
    pub sample: u64,
    /// Sum of original component pressure observations [Pa].
    pub observer_pressure_pa: f64,
    /// Oscillator plus bilateral-spring storage [J].
    pub network_energy_j: f64,
    /// Sum of all contact potentials [J].
    pub contact_energy_j: f64,
    /// Authored external work only; solved contact forces are internal [J].
    pub external_work_j: f64,
    /// Component and bilateral-dashpot losses [J].
    pub network_dissipation_j: f64,
    /// Sum of all normal contact losses [J].
    pub contact_dissipation_j: f64,
    /// Tangential dissipation of the selected friction model [J].
    /// Coulomb reports its signed actual-work defect separately.
    pub friction_dissipation_j: f64,
    /// Whole-system storage change plus losses minus external work [J].
    pub energy_residual_j: f64,
    /// Original network's absolute-plus-relative energy allowance [J].
    pub energy_tolerance_j: f64,
    /// Joint nonlinear sweeps actually performed.
    pub sweeps: usize,
    /// Diagnostics in contact construction order; none describe partial trials.
    pub contacts: Vec<ContactPointFrame>,
    /// Optional tangential diagnostics in normal-contact order. Empty unless
    /// friction was attached; None entries are explicitly frictionless.
    pub friction: Vec<Option<friction::TangentialContactFrame>>,
    /// Two-direction Coulomb diagnostics. Empty for the original 1-D model.
    pub coulomb_friction: Vec<Option<coulomb::CoulombContactFrame>>,
}

struct ContactPoint {
    contact: ModalContact,
    config: ModalContactConfig,
    column: Vec<f64>,
    storage: ContactStorage,
}
impl ContactPoint {
    fn energy(&self, x: f64) -> Result<f64, ModalCouplingError> {
        let value = finite(self.storage.hamiltonian(&[-x, 0.0]))?;
        if value < 0.0 { return Err(invalid("contact potential is negative")); }
        Ok(value)
    }
    fn penetration(&self, x: f64) -> Result<f64, ModalCouplingError> {
        let value = finite(x - self.contact.law.gaps()[0])?.max(0.0);
        limit("contact penetration", value, self.config.maximum_penetration_m)?;
        Ok(value)
    }
}

/// Several normal contacts sharing actual component motion and bilateral loads.
/// Failed trials leave the accepted components, sample clock and report intact.
/// All contact geometry is fixed; contact activity emerges from the supplied gaps.
pub struct MultiContactModalSystem {
    network: CoupledModalSystem,
    points: Vec<ContactPoint>,
    config: MultiContactConfig,
    compliance: Vec<f64>,
    forces: Vec<f64>,
    free_q: Vec<f64>,
    old_x: Vec<f64>,
    free_x: Vec<f64>,
    reactions: Vec<f64>,
    staged_points: Vec<ContactPointFrame>,
    frame: MultiContactFrame,
    friction: Option<friction::TangentialSet>,
}
impl MultiContactModalSystem {
    /// Condense the network's displacement response at every contact pair.
    /// Each obstacle and per-contact budget goes through the single-contact
    /// admission owner. Initial contact storage is included in the total energy.
    /// No contact-loaded static equilibrium or hidden initial impulse is chosen.
    pub fn new(
        network: CoupledModalSystem,
        contacts: Vec<(ModalContact, ModalContactConfig)>,
        config: MultiContactConfig,
        gate: &CancelGate,
    ) -> Result<Self, ModalCouplingError> {
        poll(Some(gate))?;
        if network.samples_rendered() != 0 || !(1..=32).contains(&config.max_contacts)
            || contacts.is_empty() || contacts.len() > config.max_contacts
            || !(1..=128).contains(&config.max_sweeps) {
            return Err(invalid("multiple contacts require a sample-zero network and bounded nonempty contact/sweep counts"));
        }
        let n = network.mode_count();
        let p = contacts.len();
        let k = network.columns.len();
        let terms = n.checked_mul(k+2).and_then(|v| v.checked_add((k+1)*(k+1)))
            .and_then(|v| v.checked_mul(p))
            .and_then(|v| n.checked_mul(p*p).and_then(|w| v.checked_add(w)))
            .ok_or_else(|| invalid("multi-contact setup work overflow"))?;
        if terms > config.max_setup_terms {
            return Err(invalid("multi-contact setup exceeds max_setup_terms"));
        }
        let mut points = Vec::with_capacity(p);
        for (contact, c) in contacts {
            poll(Some(gate))?;
            let column = contact_column(&network, &contact, c)?;
            let storage = ContactStorage::new(Box::new(ZeroStorage), 1, vec![contact.law.clone()])
                .map_err(ModalCouplingError::ContactLaw)?;
            let x = extension(&network.models, &column, 0.0)?;
            SlitContactStep::new(&contact.law, -x).map_err(ModalCouplingError::ContactLaw)?;
            let point = ContactPoint { contact, config: c, column, storage };
            point.penetration(x)?;
            points.push(point);
        }
        let mut compliance = vec![0.0; p*p];
        for j in 0..p {
            let response = network_response(&network, &points[j].column, gate)?;
            for i in 0..p {
                poll(Some(gate))?;
                compliance[i*p+j] = dot(&points[i].column, &response)?;
            }
            if compliance[j*p+j] <= 0.0 {
                return Err(invalid("each contact attachment requires positive representable network compliance"));
            }
        }
        let system = Self {
            network, points, config, compliance, forces: vec![0.0; n], free_q: vec![0.0; n],
            old_x: vec![0.0; p], free_x: vec![0.0; p], reactions: vec![0.0; p],
            staged_points: vec![ContactPointFrame::default(); p], friction: None,
            frame: MultiContactFrame { contacts: vec![ContactPointFrame::default(); p], ..MultiContactFrame::default() },
        };
        limit("total initial energy including contacts", system.total_energy_j()?, system.network.config.maximum_total_energy_j)?;
        poll(Some(gate))?;
        Ok(system)
    }

    /// Accepted components, never trial states.
    #[must_use]
    pub fn components(&self) -> &[ModalAcousticTimeModel] { self.network.components() }
    /// Flattened external-force size, in component then mode order.
    #[must_use]
    pub fn mode_count(&self) -> usize { self.network.mode_count() }
    /// Number of supplied contacts, including currently separated ones.
    #[must_use]
    pub fn contact_count(&self) -> usize { self.points.len() }
    /// Supplied obstacle and source label; no material inference is performed.
    #[must_use]
    pub fn contact_law(&self, index: usize) -> Option<&Obstacle> { self.points.get(index).map(|p| &p.contact.law) }
    /// Shared mechanical period [s].
    #[must_use]
    pub fn sample_period_s(&self) -> f64 { self.network.sample_period_s() }
    /// Number of accepted samples; no progress on a failed joint solve.
    #[must_use]
    pub fn samples_rendered(&self) -> u64 { self.network.samples_rendered() }
    /// Last complete report, absent before the first successful sample.
    #[must_use]
    pub fn last_frame(&self) -> Option<&MultiContactFrame> { (self.frame.sample != 0).then_some(&self.frame) }
    /// Accepted mechanical and contact-potential storage [J].
    pub fn total_energy_j(&self) -> Result<f64, ModalCouplingError> {
        let mut energy = self.network.total_energy_j()?;
        for point in &self.points {
            energy = finite(energy + point.energy(extension(&self.network.models, &point.column, 0.0)?)?)?;
        }
        Ok(energy)
    }
    /// Advance all bodies and contacts together or publish nothing.
    pub fn step(&mut self, external: &[f64]) -> Result<&MultiContactFrame, ModalCouplingError> {
        self.step_inner(external, None)
    }
    /// Same transaction, with cancellation checks throughout the nonlinear solve.
    pub fn step_under_gate(&mut self, external: &[f64], gate: &CancelGate)
        -> Result<&MultiContactFrame, ModalCouplingError>
    { self.step_inner(external, Some(gate)) }

    fn step_inner(&mut self, external: &[f64], gate: Option<&CancelGate>)
        -> Result<&MultiContactFrame, ModalCouplingError>
    {
        poll(gate)?;
        let before = self.total_energy_j()?;
        self.network.prepare_forces(external, gate, false)?;
        let mut at = 0;
        for model in &self.network.models {
            poll(gate)?;
            for (&mode, &state) in model.modes().iter().zip(model.states()) {
                self.free_q[at] = finite(advance_exact_zoh(mode, state, self.network.forces[at], self.network.dt)
                    .displacement_m_sqrt_kg)?;
                at += 1;
            }
        }
        let p = self.points.len();
        let mut laws = Vec::with_capacity(p);
        for i in 0..p {
            poll(gate)?;
            self.old_x[i] = dot(&self.points[i].column, &self.network.old_q)?;
            self.free_x[i] = dot(&self.points[i].column, &self.free_q)?;
            laws.push(SlitContactStep::new(&self.points[i].contact.law, -self.old_x[i])
                .map_err(ModalCouplingError::ContactLaw)?);
        }
        // Scratch guesses restart deterministically after every failed/cancelled
        // attempt. No speculative contact memory leaks into the accepted state.
        self.reactions.fill(0.0);
        if let Some(friction) = &mut self.friction {
            friction.prepare(&self.network.old_q, &self.free_q, gate)?;
        }
        let mut converged = false;
        let mut sweeps = 0;
        let mut worst = (0.0_f64, 1.0_f64);
        for sweep in 1..=self.config.max_sweeps {
            sweeps = sweep;
            for i in 0..p {
                poll(gate)?;
                let mut free = self.free_x[i];
                for j in 0..p {
                    if i != j { free = finite(free - self.compliance[i*p+j]*self.reactions[j])?; }
                }
                if let Some(friction) = &self.friction { free = friction.normal_endpoint(i, free)?; }
                let mut local = self.points[i].config;
                // Seek a tighter coordinate root to leave room for later
                // coordinates. Acceptance still uses the user's ORIGINAL caps.
                if local.force_absolute_tolerance_n*0.125 > 0.0 { local.force_absolute_tolerance_n *= 0.125; }
                if local.force_relative_tolerance*0.125 > 0.0 { local.force_relative_tolerance *= 0.125; }
                self.reactions[i] = match solve_contact(&laws[i], self.old_x[i], free,
                    self.compliance[i*p+i], self.network.dt, local, gate) {
                    Ok((r, _)) => r,
                    // An early coordinate can exceed its bound before another
                    // contact relieves it. Keep a bounded TRIAL, not a clipped
                    // accepted force; the joint residual below must still pass.
                    Err(ModalCouplingError::Budget { what: "contact normal force", .. }) => local.maximum_force_n,
                    Err(error) => return Err(error),
                };
                if let Some(friction) = &mut self.friction {
                    friction.solve_coordinate(i, &self.reactions, self.network.dt, local, gate)?;
                }
            }
            worst = (0.0, 1.0);
            for i in 0..p {
                poll(gate)?;
                let mut x = self.free_x[i];
                for j in 0..p { x = finite(x - self.compliance[i*p+j]*self.reactions[j])?; }
                if let Some(friction) = &self.friction { x = friction.normal_endpoint(i, x)?; }
                let (expected, _) = law_force(&laws[i], self.old_x[i], x, self.network.dt)?;
                let residual = finite(self.reactions[i] - expected)?;
                let tolerance = force_tolerance(self.reactions[i], expected, self.points[i].config)?;
                if residual.abs()/tolerance > worst.0.abs()/worst.1 { worst = (residual, tolerance); }
            }
            if let Some(friction) = &self.friction {
                friction.residuals(&self.reactions, &self.points, self.network.dt, &mut worst, gate)?;
            }
            if worst.0.abs() <= worst.1 { converged = true; break; }
        }
        if !converged {
            return Err(ModalCouplingError::ContactSolve {
                residual_n: worst.0, tolerance_n: worst.1, iterations: sweeps,
            });
        }
        self.forces.copy_from_slice(external);
        for i in 0..p {
            if self.reactions[i] != 0.0 {
                for (f, b) in self.forces.iter_mut().zip(&self.points[i].column) {
                    *f = finite(*f - b*self.reactions[i])?;
                }
            }
        }
        if let Some(friction) = &self.friction { friction.add_forces(&mut self.forces, gate)?; }
        self.network.stage_inner(&self.forces, gate)?;
        let mut contact_energy = 0.0;
        let mut contact_loss = 0.0;
        for i in 0..p {
            poll(gate)?;
            let point = &self.points[i];
            let x1 = extension(&self.network.candidates, &point.column, 0.0)?;
            let penetration = point.penetration(x1)?;
            let (expected, elastic) = law_force(&laws[i], self.old_x[i], x1, self.network.dt)?;
            let tolerance = force_tolerance(self.reactions[i], expected, point.config)?;
            let residual = finite(self.reactions[i] - expected)?;
            if residual.abs() > tolerance {
                return Err(ModalCouplingError::ContactSolve { residual_n: residual, tolerance_n: tolerance, iterations: sweeps });
            }
            let energy = point.energy(x1)?;
            let loss = finite((expected - elastic)*finite(x1-self.old_x[i])?)?;
            if loss < 0.0 { return Err(invalid("contact loss must not create energy")); }
            contact_energy = finite(contact_energy + energy)?;
            contact_loss = finite(contact_loss + loss)?;
            self.staged_points[i] = ContactPointFrame {
                normal_force_n: self.reactions[i], constitutive_residual_n: residual,
                force_tolerance_n: tolerance, penetration_before_m: point.penetration(self.old_x[i])?,
                penetration_after_m: penetration, stored_energy_j: energy, dissipation_j: loss,
            };
        }
        let friction_loss = if let Some(friction) = &mut self.friction {
            friction.stage(&self.network, &self.reactions, &self.points, sweeps, gate)?
        } else { 0.0 };
        let frame = &self.network.staged_frame;
        let network_energy = finite(frame.modal_energy_j + frame.connection_energy_j)?;
        let network_loss = finite(frame.component_dissipation_j + frame.connection_dissipation_j)?;
        let external_work = dot(external, &self.network.free_delta)?;
        let after = finite(network_energy + contact_energy)?;
        limit("total energy including contacts", after, self.network.config.maximum_total_energy_j)?;
        let balance = finite((after-before) + network_loss + contact_loss - external_work)?;
        let residual = if self.friction.is_some() { finite(balance + friction_loss)? } else { balance };
        let scale = before.max(after).max(external_work.abs()).max(network_loss.abs()).max(contact_loss).max(friction_loss);
        let tolerance = finite(self.network.config.energy_absolute_tolerance_j
            + self.network.config.energy_relative_tolerance*scale)?;
        if residual.abs() > tolerance {
            return Err(ModalCouplingError::EnergyBalance { residual_j: residual, tolerance_j: tolerance });
        }
        poll(gate)?;
        self.frame.sample = frame.sample;
        self.frame.observer_pressure_pa = frame.observer_pressure_pa;
        self.frame.network_energy_j = network_energy;
        self.frame.contact_energy_j = contact_energy;
        self.frame.external_work_j = external_work;
        self.frame.network_dissipation_j = network_loss;
        self.frame.contact_dissipation_j = contact_loss;
        self.frame.friction_dissipation_j = friction_loss;
        self.frame.energy_residual_j = residual;
        self.frame.energy_tolerance_j = tolerance;
        self.frame.sweeps = sweeps;
        std::mem::swap(&mut self.frame.contacts, &mut self.staged_points);
        if let Some(friction) = &mut self.friction { friction.publish(&mut self.frame.friction, &mut self.frame.coulomb_friction); }
        self.network.publish_staged();
        Ok(&self.frame)
    }
}
