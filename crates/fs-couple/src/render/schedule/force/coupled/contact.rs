//! A two-body unilateral contact on the existing spring-connected modal network.
//!
//! Relative closure is x = left^T q_left - right^T q_right. The supplied
//! fs-dcontact obstacle has one unit opening coordinate, Phi=[-1], so its
//! opening is -x and penetration is x-gap. Its existing discrete potential
//! gradient and nonadhesive Hunt-Crossley rule supply a NONNEGATIVE normal
//! reaction R, applied as -B R to both bodies, not a prescribed force history.
//!
//! For this fixed linear network, x1 = x_free - S R, where S is the network's
//! exact held-force displacement response INCLUDING bilateral connections.
//! Only this scalar normal reaction is nonlinear. A bounded bracket solve
//! drives the existing contact law; the original modal stepper then stages the
//! actual states. Constitutive residual, penetration and total energy are
//! rechecked on those states before any component is published.
//!
//! One normal contact is supported. This is not a rigid impact, friction law,
//! free rigid-body integrator or a multi-contact complementarity solver. The
//! contact owner allocates during discrete-gradient evaluation; no hard-real-
//! time or no-allocation claim is made. Contact can generate high frequencies:
//! the inherited LINEAR Nyquist screen is not an anti-aliasing guarantee.

use fs_dcontact::{ContactStorage, Obstacle};
use fs_exec::CancelGate;
use fs_phs::Storage;

use super::{
    CoupledModalSystem, ModalAttachment, ModalCouplingError, check_solve, dot,
    extension, finite, invalid, limit, poll,
};
use crate::modal_acoustic_time::{ModalAcousticState, ModalAcousticTimeModel, advance_exact_zoh};
use crate::unilateral_contact::SlitContactStep;

/// One physical attachment pair and an explicitly sourced/authored contact law.
#[derive(Clone, Debug)]
pub struct ModalContact {
    /// Motion on this side increases closure; contact pushes against it.
    pub left: ModalAttachment,
    /// Motion on this side decreases closure; receives the opposite reaction.
    pub right: ModalAttachment,
    /// One-point fs-dcontact obstacle with collocation [-1]. Its gap, weight,
    /// stiffness, exponent, loss coefficient and provenance are used verbatim.
    /// Both Obstacle::new and Obstacle::from_receipt are supported.
    pub law: Obstacle,
}

/// Explicit per-step nonlinear work and physical limits, not clipping targets.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ModalContactConfig {
    /// Maximum bisection refinements; in 1..=128. At most two bracket endpoint
    /// trials precede them and one final actual-state constitutive check follows.
    pub max_iterations: usize,
    /// Upper normal reaction bound [N], also used for root bracketing.
    pub maximum_force_n: f64,
    /// Maximum permitted positive endpoint penetration [m].
    pub maximum_penetration_m: f64,
    /// Absolute constitutive equation residual tolerance [N].
    pub force_absolute_tolerance_n: f64,
    /// Relative residual tolerance, in (0,1). Scaled by the larger of applied
    /// and constitutive reactions, not by the much larger force ceiling.
    pub force_relative_tolerance: f64,
}

/// Last accepted whole-system step, including contact storage and loss.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ContactModalFrame {
    /// One-based completed sample clock.
    pub sample: u64,
    /// Sum of the original pressure observations [Pa].
    pub observer_pressure_pa: f64,
    /// Modal oscillator plus bilateral-spring energy [J].
    pub network_energy_j: f64,
    /// fs-dcontact potential at the accepted endpoint [J].
    pub contact_energy_j: f64,
    /// Work by authored external modal forces, EXCLUDING contact reactions [J].
    pub external_work_j: f64,
    /// Original component plus bilateral-dashpot loss [J].
    pub network_dissipation_j: f64,
    /// Nonadhesive contact loss (constitutive R - conservative R)*delta_x [J].
    pub contact_dissipation_j: f64,
    /// Whole-system storage change + all losses - authored external work [J].
    pub energy_residual_j: f64,
    /// Actual whole-system energy tolerance [J].
    pub energy_tolerance_j: f64,
    /// Applied compressive reaction [N], never attractive.
    pub normal_force_n: f64,
    /// Applied reaction minus law evaluated at ACTUAL accepted endpoints [N].
    pub constitutive_residual_n: f64,
    /// Absolute-plus-relative force tolerance used at those endpoints [N].
    pub force_tolerance_n: f64,
    /// Positive penetration at start and end of this step [m].
    pub penetration_before_m: f64,
    /// Positive penetration at the endpoint [m]; zero means separated/touching.
    pub penetration_after_m: f64,
    /// Bisection refinements actually used (zero for a bracket-endpoint root).
    pub iterations: usize,
}

/// A linear network plus one implicit compliant contact. Accepted states stay
/// transactional across the extra constitutive, penetration and energy gates.
/// Construction begins at network sample zero with caller-supplied vibration;
/// it does NOT silently solve a contact-loaded static equilibrium.
pub struct ContactModalSystem {
    network: CoupledModalSystem,
    contact: ModalContact,
    config: ModalContactConfig,
    column: Vec<f64>,
    compliance_m_per_n: f64,
    forces: Vec<f64>,
    storage: ContactStorage,
    last: Option<ContactModalFrame>,
}

struct ZeroStorage;
impl Storage for ZeroStorage {
    fn hamiltonian(&self, _: &[f64]) -> f64 { 0.0 }
    fn gradient(&self, _: &[f64], out: &mut [f64]) { out.fill(0.0); }
}

impl ContactModalSystem {
    /// Attach an explicitly declared contact to a sample-zero modal network.
    /// Initial vibration and any penetration energy are retained and budgeted.
    /// Malformed raw-parts obstacles are checked before any storage can index
    /// their arrays. Setup adds O(modes*connections + connections^2) work under
    /// the network's already admitted mode/connection ceilings.
    pub fn new(
        network: CoupledModalSystem,
        contact: ModalContact,
        config: ModalContactConfig,
        gate: &CancelGate,
    ) -> Result<Self, ModalCouplingError> {
        poll(Some(gate))?;
        if network.samples_rendered() != 0 {
            return Err(invalid("contact admission requires a sample-zero network; no mid-run energy insertion"));
        }
        let column = contact_column(&network, &contact, config)?;
        let compliance_m_per_n = effective_compliance(&network, &column, gate)?;
        let storage = ContactStorage::new(Box::new(ZeroStorage), 1, vec![contact.law.clone()])
            .map_err(ModalCouplingError::ContactLaw)?;
        let system = Self { forces: vec![0.0; network.mode_count()], network, contact,
            config, column, compliance_m_per_n, storage, last: None };
        let x = extension(&system.network.models, &system.column, 0.0)?;
        system.check_penetration(x)?;
        limit("total initial energy including contact", system.total_energy_j()?, system.network.config.maximum_total_energy_j)?;
        // Validates the shared discrete-law evaluation at the actual starting state.
        SlitContactStep::new(&system.contact.law, -x).map_err(ModalCouplingError::ContactLaw)?;
        poll(Some(gate))?;
        Ok(system)
    }

    /// Accepted original modal components; no mutable bypass of the contact.
    #[must_use]
    pub fn components(&self) -> &[ModalAcousticTimeModel] { self.network.components() }
    /// Flattened external-force count, in component then mode order.
    #[must_use]
    pub fn mode_count(&self) -> usize { self.network.mode_count() }
    /// Shared mechanical period [s].
    #[must_use]
    pub fn sample_period_s(&self) -> f64 { self.network.sample_period_s() }
    /// Complete accepted samples; cancelled/refused trials never advance it.
    #[must_use]
    pub fn samples_rendered(&self) -> u64 { self.network.samples_rendered() }
    /// Supplied obstacle, with its unmodified provenance and physical coefficients.
    #[must_use]
    pub fn contact_law(&self) -> &Obstacle { &self.contact.law }
    /// Last complete contact-inclusive diagnostic, absent before the first sample.
    #[must_use]
    pub const fn last_frame(&self) -> Option<&ContactModalFrame> { self.last.as_ref() }
    /// Accepted network plus contact-potential energy [J].
    pub fn total_energy_j(&self) -> Result<f64, ModalCouplingError> {
        let x = extension(&self.network.models, &self.column, 0.0)?;
        finite(self.network.total_energy_j()? + self.contact_energy(x)?)
    }
    /// Advance one held-force sample; neither state nor last report changes on refusal.
    pub fn step(&mut self, external: &[f64]) -> Result<&ContactModalFrame, ModalCouplingError> {
        self.step_inner(external, None)
    }
    /// Same transaction, polling cancellation throughout root evaluation and staging.
    pub fn step_under_gate(&mut self, external: &[f64], gate: &CancelGate)
        -> Result<&ContactModalFrame, ModalCouplingError>
    {
        self.step_inner(external, Some(gate))
    }

    fn step_inner(&mut self, external: &[f64], gate: Option<&CancelGate>)
        -> Result<&ContactModalFrame, ModalCouplingError>
    {
        poll(gate)?;
        let before = self.total_energy_j()?;
        // This is a scratch prediction, not a component step. Final physical
        // caps are checked after contact has acted, not on unrestrained motion.
        self.network.prepare_forces(external, gate, false)?;
        let x0 = dot(&self.column, &self.network.old_q)?;
        let mut x_free = 0.0;
        let mut index = 0;
        for model in &self.network.models {
            poll(gate)?;
            for (&mode, &state) in model.modes().iter().zip(model.states()) {
                let q = advance_exact_zoh(mode, state, self.network.forces[index], self.network.dt)
                    .displacement_m_sqrt_kg;
                x_free = finite(x_free + self.column[index]*q)?;
                index += 1;
            }
        }
        let law = SlitContactStep::new(&self.contact.law, -x0).map_err(ModalCouplingError::ContactLaw)?;
        let (reaction, iterations) = solve_contact(&law, x0, x_free, self.compliance_m_per_n,
            self.network.dt, self.config, gate)?;
        self.forces.copy_from_slice(external);
        if reaction != 0.0 {
            for (force, b) in self.forces.iter_mut().zip(&self.column) { *force = finite(*force-b*reaction)?; }
        }
        self.network.stage_inner(&self.forces, gate)?;
        let x1 = extension(&self.network.candidates, &self.column, 0.0)?;
        let penetration_after_m = self.check_penetration(x1)?;
        let (expected, elastic) = law_force(&law, x0, x1, self.network.dt)?;
        let force_tolerance_n = force_tolerance(reaction, expected, self.config)?;
        let constitutive_residual_n = finite(reaction-expected)?;
        if constitutive_residual_n.abs() > force_tolerance_n {
            return Err(ModalCouplingError::ContactSolve { residual_n: constitutive_residual_n,
                tolerance_n: force_tolerance_n, iterations });
        }
        let contact_energy_j = self.contact_energy(x1)?;
        let delta = finite(x1-x0)?;
        let contact_dissipation_j = finite((expected-elastic)*delta)?;
        if contact_dissipation_j < 0.0 { return Err(invalid("contact loss must not create energy")); }
        let frame = &self.network.staged_frame;
        let network_energy_j = finite(frame.modal_energy_j + frame.connection_energy_j)?;
        let network_dissipation_j = finite(frame.component_dissipation_j + frame.connection_dissipation_j)?;
        let external_work_j = dot(external, &self.network.free_delta)?;
        let after = finite(network_energy_j + contact_energy_j)?;
        limit("total energy including contact", after, self.network.config.maximum_total_energy_j)?;
        let residual = finite((after-before) + network_dissipation_j + contact_dissipation_j - external_work_j)?;
        let scale = before.max(after).max(external_work_j.abs()).max(network_dissipation_j.abs()).max(contact_dissipation_j);
        let tolerance = finite(self.network.config.energy_absolute_tolerance_j
            + self.network.config.energy_relative_tolerance * scale)?;
        if residual.abs() > tolerance {
            return Err(ModalCouplingError::EnergyBalance { residual_j: residual, tolerance_j: tolerance });
        }
        let accepted = ContactModalFrame {
            sample: frame.sample, observer_pressure_pa: frame.observer_pressure_pa,
            network_energy_j, contact_energy_j, external_work_j, network_dissipation_j,
            contact_dissipation_j, energy_residual_j: residual, energy_tolerance_j: tolerance,
            normal_force_n: reaction, constitutive_residual_n, force_tolerance_n,
            penetration_before_m: (x0-self.contact.law.gaps()[0]).max(0.0),
            penetration_after_m, iterations,
        };
        poll(gate)?;
        self.network.publish_staged();
        Ok(self.last.insert(accepted))
    }

    fn contact_energy(&self, closure: f64) -> Result<f64, ModalCouplingError> {
        let value = finite(self.storage.hamiltonian(&[-closure, 0.0]))?;
        if value < 0.0 { return Err(invalid("contact potential is negative")); }
        Ok(value)
    }
    fn check_penetration(&self, closure: f64) -> Result<f64, ModalCouplingError> {
        let penetration = finite(closure-self.contact.law.gaps()[0])?.max(0.0);
        limit("contact penetration", penetration, self.config.maximum_penetration_m)?;
        Ok(penetration)
    }
}

// Shared single/multiple-contact admission: obstacle, limits and signed basis.
pub(super) fn contact_column(network: &CoupledModalSystem, contact: &ModalContact, config: ModalContactConfig)
    -> Result<Vec<f64>, ModalCouplingError>
{
    if !(1..=128).contains(&config.max_iterations)
        || [config.maximum_force_n, config.maximum_penetration_m,
            config.force_absolute_tolerance_n, config.force_relative_tolerance]
            .iter().any(|x| !x.is_finite() || *x <= 0.0)
        || config.force_relative_tolerance >= 1.0 {
        return Err(invalid("contact requires explicit positive finite force, penetration, iteration and residual budgets"));
    }
    let law = &contact.law;
    if law.n_points() != 1 || law.collocation() != [-1.0]
        || law.gaps().len() != 1 || law.weights().len() != 1
        || !law.gaps()[0].is_finite() || !law.weights()[0].is_finite() || law.weights()[0] < 0.0
        || !law.stiffness().is_finite() || law.stiffness() < 0.0
        || !law.alpha().is_finite() || law.alpha() < 1.0
        || !law.internal_loss().is_finite() || law.internal_loss() < 0.0
        || law.provenance().trim().is_empty() {
        return Err(invalid("contact requires a finite provenance-labelled one-point unit-opening obstacle"));
    }
    let mut column = vec![0.0; network.mode_count()];
    for (attachment, sign) in [(&contact.left, 1.0), (&contact.right, -1.0)] {
        let model = network.models.get(attachment.component)
            .ok_or_else(|| invalid("contact attachment names an unknown component"))?;
        if attachment.shapes.len() != model.modes().len() || attachment.shapes.iter().any(|x| !x.is_finite()) {
            return Err(invalid("contact shapes must match the finite mass-normalized component basis"));
        }
        for (k, b) in attachment.shapes.iter().enumerate() {
            let index = network.offsets[attachment.component] + k;
            column[index] = finite(column[index] + sign*b)?;
        }
    }
    Ok(column)
}

// Condense the ALREADY ADMITTED bilateral network, not a different integrator.
// S = b^T[D - D B sqrt(h) A^-1 sqrt(h) B^T D]b.
fn effective_compliance(network: &CoupledModalSystem, column: &[f64], gate: &CancelGate)
    -> Result<f64, ModalCouplingError>
{
    let response = network_response(network, column, gate)?;
    let s = dot(column, &response)?;
    if s <= 0.0 { return Err(invalid("contact attachment needs positive representable network compliance")); }
    Ok(s)
}

// Full displacement response of the existing bilateral network to a unit
// attachment load. Cross-contact compliance is b_i^T response(b_j).
fn network_response(network: &CoupledModalSystem, column: &[f64], gate: &CancelGate)
    -> Result<Vec<f64>, ModalCouplingError>
{
    let mut d = Vec::with_capacity(column.len());
    for model in &network.models {
        poll(Some(gate))?;
        for &mode in model.modes() {
            d.push(finite(advance_exact_zoh(mode, ModalAcousticState::default(), 1.0, network.dt)
                .displacement_m_sqrt_kg)?);
        }
    }
    let mut response: Vec<f64> = d.iter().zip(column).map(|(d,b)| finite(d*b)).collect::<Result<_,_>>()?;
    let mut rhs = Vec::with_capacity(network.columns.len());
    for (j, b) in network.columns.iter().enumerate() { rhs.push(finite(network.roots[j]*dot(b,&response)?)?); }
    let mut solution = rhs.clone();
    network.factor.solve(&mut solution);
    check_solve(&network.matrix, &solution, &rhs, network.config.solve_relative_tolerance)?;
    for (j, b) in network.columns.iter().enumerate() {
        poll(Some(gate))?;
        for k in 0..response.len() {
            response[k] = finite(response[k] - d[k]*b[k]*network.roots[j]*solution[j])?;
        }
    }
    Ok(response)
}

fn law_force(law: &SlitContactStep, x0: f64, x1: f64, dt: f64)
    -> Result<(f64,f64), ModalCouplingError>
{
    let (elastic,damping) = law.coefficients(-x1).map_err(ModalCouplingError::ContactLaw)?;
    // Opening velocity is minus closure velocity. Nonadhesive unloading uses
    // exactly the existing slit law; it is not an attractive spring at exit.
    let velocity = finite((x1-x0)/dt)?;
    let force = finite(elastic + damping*velocity)?.max(0.0);
    Ok((force,elastic))
}
pub(super) fn force_tolerance(applied: f64, expected: f64, c: ModalContactConfig) -> Result<f64,ModalCouplingError> {
    finite(c.force_absolute_tolerance_n + c.force_relative_tolerance*applied.abs().max(expected.abs()))
}
fn solve_contact(law: &SlitContactStep, x0: f64, free: f64, compliance: f64, dt: f64,
    config: ModalContactConfig, gate: Option<&CancelGate>) -> Result<(f64,usize),ModalCouplingError>
{
    let evaluate = |reaction: f64| -> Result<(f64,f64),ModalCouplingError> {
        poll(gate)?;
        let x = finite(free - compliance*reaction)?;
        let (expected,_) = law_force(law,x0,x,dt)?;
        Ok((finite(reaction-expected)?,force_tolerance(reaction,expected,config)?))
    };
    solve_reaction(evaluate, config)
}

// One bracket owner for dynamic discrete forces and stationary potentials.
pub(super) fn solve_reaction(
    evaluate: impl Fn(f64) -> Result<(f64, f64), ModalCouplingError>,
    config: ModalContactConfig,
) -> Result<(f64, usize), ModalCouplingError> {
    let (at_zero,_) = evaluate(0.0)?;
    if at_zero == 0.0 { return Ok((0.0,0)); }
    let mut lo = 0.0;
    // The convex potential secant and nonadhesive loss are monotone in closure.
    // Thus R(0) bounds the root; the physical force ceiling is never increased.
    let mut hi = (-at_zero).min(config.maximum_force_n);
    let (at_hi,tol_hi) = evaluate(hi)?;
    if at_hi.abs() <= tol_hi { return Ok((hi,0)); }
    if at_hi < 0.0 {
        return Err(ModalCouplingError::Budget { what: "contact normal force", value: hi-at_hi,
            limit: config.maximum_force_n });
    }
    let mut residual = at_hi;
    let mut tolerance = tol_hi;
    let mut used = 0;
    for iteration in 1..=config.max_iterations {
        used = iteration;
        let mid = f64::midpoint(lo,hi);
        let (r,t) = evaluate(mid)?;
        residual=r; tolerance=t;
        if r.abs() <= t { return Ok((mid,iteration)); }
        if mid == lo || mid == hi { break; }
        if r < 0.0 { lo=mid; } else { hi=mid; }
    }
    Err(ModalCouplingError::ContactSolve { residual_n: residual, tolerance_n: tolerance, iterations: used })
}

/// Simultaneous normal contacts sharing the same mechanical network.
pub mod multiple;
