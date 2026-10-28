//! Independent compliant obstacles and prescribed-load friction on ONE network.
//!
//! Normal closure x=N^T q applies -N R; a driven friction port applies +B F.
//! Both mixed responses are retained, including every bilateral body reaction.
//! The original normal potential/secant, normal root, friction root and modal
//! transition own the physics. This module supplies their joint coordinates,
//! not a second contact law or time integrator. By default these reactions are
//! independent of bow loads; explicit `normal_load` bindings use a reaction as
//! the corresponding tangential law's normal load in this SAME joint solve.

use super::*;
use super::super::super::{
    ModalContact, ModalContactConfig, ZeroStorage, contact_column, extension,
    force_tolerance, law_force, solve_contact,
};
use fs_contact::normal_patch::source::SourcePointContact;
use fs_dcontact::{ContactStorage, OpeningContactStep};
use fs_phs::Storage;

/// Prescribed normal support translation with work-conjugate reactions.
pub mod motion;
/// Prescribed attachment motion with a discrete product-rule work balance.
pub mod geometry;

/// Admission for optional normal contacts in the existing joint friction solve.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FrictionObstacleConfig {
    /// At most 32 normal contacts, with at most 64 normal + friction ports.
    pub max_contacts: usize,
    /// Conservative combined setup bound: s*(n*(k+2)+(k+1)^2)+n*s^2,
    /// where s is normal + friction ports, n modes and k bilateral links.
    /// This includes existing friction setup; it is not an incremental allowance.
    pub max_setup_terms: usize,
}

/// Actual normal law and energy on the same endpoint as every friction port.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct FrictionObstacleFrame {
    /// Applied nonnegative reaction [N], with modal load -N R. A changing
    /// attachment uses the mean of its two endpoint columns for this load.
    pub normal_force_n: f64,
    /// Original unilateral potential at the accepted endpoint [J].
    pub stored_energy_j: f64,
    /// Original nonadhesive contact loss over the complete step [J].
    pub dissipation_j: f64,
    /// Positive endpoint penetration [m]; zero means separated or touching.
    pub penetration_m: f64,
    /// Applied minus constitutive reaction at the ACTUAL endpoint [N].
    pub constitutive_residual_n: f64,
    /// This contact's absolute-plus-relative residual allowance [N].
    pub force_tolerance_n: f64,
    /// Accepted support translation from its authored reference [m]. Positive
    /// opens the gap; the original source law and its clearance are unchanged.
    pub support_offset_m: f64,
    /// Work supplied by changing the attachment map, R delta_N^T mean_q [J].
    /// Zero for a contact that only translates along its fixed normal.
    pub attachment_work_j: f64,
    /// Total prescribed geometry work: attachment_work - R delta_offset [J].
    pub actuator_work_j: f64,
    /// This contact's share of the omitted modes' static storage, R_i (G R)_i / 2 [J].
    /// Zero unless a residual compliance is bound.
    pub residual_energy_j: f64,
    /// This contact's share of the declared numerical loss, dR_i (G dR)_i / 2 [J], where dR is
    /// the change of the held reactions over the step. Zero unless a residual compliance is bound.
    pub residual_loss_j: f64,
}

impl MultiFrictionModalSystem {
    /// Bind independent compliant stops/obstacles before any physical sample.
    /// All initial vibration and contact potential are retained and budgeted;
    /// this is not a nonlinear static-equilibrium solve. Each original obstacle
    /// may carry its existing material contact receipt. No new source authority
    /// or friction law is inferred from its coefficients or normal reaction.
    ///
    /// A zero right-hand shape denotes a fixed support. A nonzero right-hand
    /// shape receives the reciprocal reaction on its existing body. Supports
    /// retain their spatial stations unless `with_obstacle_geometry_motion` is
    /// admitted. `with_obstacle_motion` additionally admits normal support
    /// translation; both changes account for their actual actuator work.
    /// Every normal and friction law is checked together before publication.
    /// Bounded coordinate sweeps can refuse; convergence for arbitrary strong
    /// coupling/nonmonotone friction and rigid-contact limits is not promised.
    ///
    /// Empty input is an exact no-op. A nonempty set binds once, at sample zero.
    /// Bound geometric storage retains these same contact trials, including
    /// initial/final potential, loss and support work in the joint DG balance.
    /// Released-friction preload solves these contacts at equilibrium only when
    /// no additional nonlinear storage is bound.
    /// Existing linear bandwidth checks are not a nonlinear anti-aliasing proof.
    pub fn with_obstacles(mut self, contacts: Vec<(ModalContact, ModalContactConfig)>,
        limits: FrictionObstacleConfig, gate: &CancelGate) -> Result<Self, ModalFrictionError>
    {
        poll(Some(gate))?;
        if contacts.is_empty() { return Ok(self); }
        if self.samples_rendered() != 0 || self.obstacles.is_some() {
            return Err(invalid("friction obstacles bind once before the first physical sample").into());
        }
        let set = ObstacleSet::new(&self, contacts, limits, gate)?;
        let energy = finite(self.core.total_energy_j()? + set.energy(self.components())?)?;
        limit("initial energy including friction obstacles", energy,
            self.core.network.config.maximum_total_energy_j)?;
        poll(Some(gate))?;
        self.obstacles = Some(set);
        Ok(self)
    }

    /// Retain source-bound normal laws and their full applicability envelopes.
    ///
    /// Bind selected obstacle indices before playback. The original obstacle must
    /// be the same unit-point Hertz/Hunt--Crossley law as its source; coefficients
    /// cannot override material resolution. The source's bulk/interface state and
    /// local relative curvature remain frozen. Moving geometry changes the
    /// contact coordinates and work, not that declared local contact model.
    ///
    /// Initial motion, both candidate indentations, the ACTUAL secant loss rate
    /// and sample duration must be admitted. Mathematical root iterates do not
    /// publish state; the completed mechanical candidate is checked before any
    /// normal, friction, thermal, geometry or clock publication. Static preload
    /// retains the existing penetration ceiling, bounded by the same envelope.
    pub fn with_obstacle_sources(mut self, sources: Vec<(usize, SourcePointContact)>,
        gate: &CancelGate) -> Result<Self, ModalFrictionError>
    {
        poll(Some(gate))?;
        if sources.is_empty() { return Ok(self); }
        if self.samples_rendered() != 0 {
            return Err(invalid("normal material sources bind before the first physical sample").into());
        }
        let set = self.obstacles.as_ref()
            .ok_or_else(|| invalid("normal material sources require bound obstacles"))?;
        let velocities: Vec<_> = self.components().iter().flat_map(|m|
            m.states().iter().map(|s| s.velocity_m_sqrt_kg_per_s)).collect();
        for (entry, (index, source)) in sources.iter().enumerate() {
            poll(Some(gate))?;
            let point = set.points.get(*index)
                .ok_or_else(|| invalid("normal material source names an absent obstacle"))?;
            if point.source.is_some() || sources[..entry].iter().any(|(prior, _)| prior == index) {
                return Err(invalid("each normal material source binds once").into());
            }
            let law = &point.contact.law;
            if law.stiffness().to_bits() != source.stiffness_n_per_m_3_2().to_bits()
                || law.alpha().to_bits() != 1.5_f64.to_bits()
                || law.weights() != [1.0]
                || law.internal_loss().to_bits() != source.dissipation_s_per_m().to_bits()
                || point.config.maximum_penetration_m > source.envelope().maximum_indentation_m {
                return Err(invalid("normal obstacle law and penetration ceiling must retain the material source").into());
            }
            let offset = set.motion.as_ref().map_or(0.0, |m| m.accepted[*index]);
            let x = finite(extension(self.components(), &point.column, 0.0)? - offset)?;
            source.validate_kinematics(point.penetration(x)?, dot(&point.column, &velocities)?,
                self.sample_period_s()).map_err(ModalCouplingError::ContactMaterial)?;
        }
        poll(Some(gate))?;
        let set = self.obstacles.as_mut().expect("validated obstacle set");
        for (index, source) in sources { set.points[index].source = Some(source); }
        Ok(self)
    }

    /// Complete immutable source receipts for a material-bound normal contact.
    #[must_use]
    pub fn obstacle_source(&self, index: usize) -> Option<&SourcePointContact> {
        self.obstacles.as_ref()?.points.get(index)?.source.as_ref()
    }

    /// Bind the static compliance of modes the network does NOT retain, among the
    /// bound obstacles, once and before the first sample (pexfp.299).
    ///
    /// `compliance` is row-major q x q [m/N], symmetric and positive semidefinite:
    /// entry (i, j) is obstacle i's closure under a unit reaction of obstacle j
    /// through the omitted modes. Each closure then reads x = N^T q - G R. The
    /// omitted modes respond quasi-statically to each step's held reactions, so
    /// the joint solve adds G to the dynamic compliance and the ledger stores
    /// R G R / 2. Over a step whose reactions change by dR, dR G dR / 2 >= 0
    /// leaves the mechanics: that is reported as declared NUMERICAL loss, apart
    /// from every contact's constitutive loss. The bound producer owns G.
    ///
    /// Two obstacles whose rows of the complete joint compliance are exactly
    /// opposed (a pinching pair at one station) act through one net force. A
    /// coordinate sweep between them converges at rate (kA/(1+kA))^2 per sweep,
    /// which G makes slow, so the pair is solved together on its net force.
    ///
    /// Translating or moving supports, released static preload and bound
    /// nonlinear storage refuse: their residual is not wired.
    pub fn with_obstacle_residual_compliance(mut self, compliance: Vec<f64>, gate: &CancelGate)
        -> Result<Self, ModalFrictionError>
    {
        poll(Some(gate))?;
        if self.samples_rendered() != 0 || self.core.has_nonlinear_potentials() {
            return Err(invalid("obstacle residual compliance binds once, before playback and without nonlinear storage").into());
        }
        let set = self.obstacles.as_mut()
            .ok_or_else(|| invalid("obstacle residual compliance requires bound obstacles"))?;
        if set.residual.is_some() || set.motion.is_some() || set.geometry.is_some() {
            return Err(invalid("obstacle residual compliance binds once, on fixed supports").into());
        }
        let q = set.points.len();
        if compliance.len() != q*q || compliance.iter().any(|g| !g.is_finite()) {
            return Err(invalid("obstacle residual compliance needs one finite entry per obstacle pair").into());
        }
        let mut largest = 0.0_f64;
        for i in 0..q {
            if compliance[i*q+i] < 0.0 {
                return Err(invalid("obstacle residual compliance must be positive semidefinite").into());
            }
            largest = largest.max(compliance[i*q+i]);
            for j in 0..i {
                if compliance[i*q+j].to_bits() != compliance[j*q+i].to_bits() {
                    return Err(invalid("obstacle residual compliance must be exactly symmetric").into());
                }
            }
        }
        let (eigenvalues, _) = fs_la::eigen::jacobi_eigh(&compliance, q);
        if eigenvalues.iter().any(|&l| !(l >= -1e-12*largest)) {
            return Err(invalid("obstacle residual compliance must be positive semidefinite").into());
        }
        for (a, g) in set.nn.iter_mut().zip(&compliance) { *a = finite(*a+g)?; }
        let mut paired = vec![false; q];
        let mut blocks = Vec::with_capacity(q);
        for i in 0..q {
            if paired[i] { continue; }
            let partner = (i+1..q).find(|&j| !paired[j] && opposed(&set.nn, q, i, j));
            paired[i] = true;
            match partner {
                Some(j) => { paired[j] = true; blocks.push(Block::Pair(i, j)); },
                None => blocks.push(Block::Single(i)),
            }
        }
        set.residual = Some(Residual { compliance, blocks });
        poll(Some(gate))?;
        Ok(self)
    }

    /// Number of independent normal obstacles, not the prescribed bow loads.
    #[must_use]
    pub fn obstacle_count(&self) -> usize { self.obstacles.as_ref().map_or(0, |o| o.points.len()) }

    /// Bound residual compliance among the obstacles, row-major [m/N], if any.
    #[must_use]
    pub fn obstacle_residual_compliance(&self) -> Option<&[f64]> {
        self.obstacles.as_ref()?.residual.as_ref().map(|r| r.compliance.as_slice())
    }

    /// Obstacle pairs solved together on their net force, in input indices.
    #[must_use]
    pub fn obstacle_pairs(&self) -> Vec<(usize, usize)> {
        self.obstacles.as_ref().and_then(|o| o.residual.as_ref()).map_or_else(Vec::new, |r|
            r.blocks.iter().filter_map(|b| match *b { Block::Pair(i, j) => Some((i, j)), Block::Single(_) => None })
                .collect())
    }

    /// Original geometry, law and its optional source receipt, in input order.
    #[must_use]
    pub fn obstacle(&self, index: usize) -> Option<&ModalContact> {
        self.obstacles.as_ref()?.points.get(index).map(|p| &p.contact)
    }

    /// Last jointly accepted normals; a rejected trial cannot change this slice.
    #[must_use]
    pub fn last_obstacle_frames(&self) -> Option<&[FrictionObstacleFrame]> {
        if !self.accepted { return None; }
        Some(&self.obstacles.as_ref()?.accepted)
    }
}

struct NormalPoint {
    contact: ModalContact,
    source: Option<SourcePointContact>,
    config: ModalContactConfig,
    column: Vec<f64>,
    storage: ContactStorage,
    law: OpeningContactStep,
    old: f64,
    free: f64,
    iterations: usize,
}
impl NormalPoint {
    fn energy(&self, x: f64) -> Result<f64, ModalCouplingError> {
        let energy = finite(self.storage.hamiltonian(&[-x, 0.0]))?;
        if energy < 0.0 { return Err(invalid("obstacle potential must be nonnegative")); }
        Ok(energy)
    }
    fn penetration(&self, x: f64) -> Result<f64, ModalCouplingError> {
        let penetration = finite(x-self.contact.law.gaps()[0])?.max(0.0);
        limit("friction obstacle penetration", penetration, self.config.maximum_penetration_m)?;
        Ok(penetration)
    }
    fn check_source(&self, x: f64, dt: f64) -> Result<(), ModalCouplingError> {
        if let Some(source) = &self.source {
            // Closure includes the source gap, support translation and moving
            // participation map. This is exactly law_force's dissipative rate,
            // including separated-to-contact and contact-to-separated steps.
            let rate = finite(finite(x-self.old)?/dt)?;
            for endpoint in [self.old, x] {
                source.validate_kinematics(self.penetration(endpoint)?, rate, dt)
                    .map_err(ModalCouplingError::ContactMaterial)?;
            }
        }
        Ok(())
    }
    // `share` of the declared tolerance bounds the residual: all of it at the
    // actual endpoint, half of it on a prediction solved with headroom.
    fn check(&self, reaction: f64, x: f64, dt: f64, share: f64) -> Result<(f64, f64, f64), ModalCouplingError> {
        let (expected, elastic) = law_force(&self.law, self.old, x, dt)?;
        let tolerance = force_tolerance(reaction, expected, self.config)?;
        let residual = finite(reaction-expected)?;
        if residual.abs() > share*tolerance {
            return Err(ModalCouplingError::ContactSolve { residual_n: residual,
                tolerance_n: share*tolerance, iterations: self.iterations });
        }
        // Exactly the original nonadhesive constitutive loss, not a new damper
        // or a residual disguised as physical dissipation.
        let loss = finite((expected-elastic)*finite(x-self.old)?)?;
        if loss < 0.0 { return Err(invalid("obstacle loss must not create energy")); }
        Ok((loss, residual, tolerance))
    }
}

pub(super) struct ObstacleSet {
    points: Vec<NormalPoint>,
    limits: FrictionObstacleConfig,
    // Both mixed directions are evaluated; do not assume numerical symmetry.
    nn: Vec<f64>, nt: Vec<f64>, tn: Vec<f64>,
    reactions: Vec<f64>, friction_ports: usize,
    staged: Vec<FrictionObstacleFrame>, accepted: Vec<FrictionObstacleFrame>,
    motion: Option<motion::MotionState>,
    geometry: Option<geometry::GeometryState>,
    absolute_work_j: f64,
    // The prediction and the propagated endpoint are rounded along different
    // paths. A step whose actual endpoint missed is solved once more with
    // headroom: tighter roots and a prediction within half the tolerance.
    headroom: bool,
    endpoint_missed: bool,
    residual: Option<Residual>,
}

// Omitted-mode static compliance among the contacts, already added to nn.
struct Residual {
    compliance: Vec<f64>,
    blocks: Vec<Block>,
}
#[derive(Clone, Copy)]
enum Block { Single(usize), Pair(usize, usize) }

// Exactly opposed 2x2 block of the complete joint compliance A: x_i = f_i - a(R_i - R_j)
// and x_j = f_j + a(R_i - R_j). Bitwise, so a pair is never formed from rounding.
fn opposed(a: &[f64], q: usize, i: usize, j: usize) -> bool {
    let d = a[i*q+i];
    d > 0.0 && a[j*q+j].to_bits() == d.to_bits()
        && a[i*q+j].to_bits() == (-d).to_bits() && a[j*q+i].to_bits() == (-d).to_bits()
}

impl ObstacleSet {
    pub(super) const fn has_residual(&self) -> bool { self.residual.is_some() }
    // Closure lag of contact i through the omitted modes under reactions `r`: (G r)_i.
    fn lag(&self, i: usize, reactions: &[f64]) -> Result<f64, ModalCouplingError> {
        match &self.residual {
            Some(residual) => { let q = self.points.len(); dot(&residual.compliance[i*q..(i+1)*q], reactions) },
            None => Ok(0.0),
        }
    }
    fn accepted_reactions(&self) -> Vec<f64> { self.accepted.iter().map(|f| f.normal_force_n).collect() }
    // Both contacts of an exactly opposed pair on their net force D = R_i - R_j.
    // h(D) = D - L_i(f_i - a D) + L_j(f_j + a D) rises with slope at least one, so
    // [min(0, -h(0)), max(0, -h(0))] brackets the root, as in the scalar owner.
    // Accepts once both laws hold at the closures their OWN reactions imply.
    fn solve_pair(&self, i: usize, j: usize, f_i: f64, f_j: f64, dt: f64, gate: Option<&CancelGate>)
        -> Result<(f64, f64, usize), ModalCouplingError>
    {
        let q = self.points.len();
        let a = self.nn[i*q+i];
        let (pi, pj) = (&self.points[i], &self.points[j]);
        let share = if self.headroom { 0.125 } else { 0.5 };
        let reactions = |d: f64| -> Result<(f64, f64), ModalCouplingError> {
            poll(gate)?;
            let (ri, _) = law_force(&pi.law, pi.old, finite(f_i-a*d)?, dt)?;
            let (rj, _) = law_force(&pj.law, pj.old, finite(f_j+a*d)?, dt)?;
            Ok((ri, rj))
        };
        let settled = |ri: f64, rj: f64| -> Result<bool, ModalCouplingError> {
            let net = finite(ri-rj)?;
            let (ei, _) = law_force(&pi.law, pi.old, finite(f_i-a*net)?, dt)?;
            let (ej, _) = law_force(&pj.law, pj.old, finite(f_j+a*net)?, dt)?;
            Ok((ri-ei).abs() <= share*force_tolerance(ri, ei, pi.config)?
                && (rj-ej).abs() <= share*force_tolerance(rj, ej, pj.config)?)
        };
        let (r0i, r0j) = reactions(0.0)?;
        let h0 = finite(r0j-r0i)?;
        if settled(r0i, r0j)? { return Ok((r0i, r0j, 0)); }
        let (mut lo, mut hi) = if h0 < 0.0 { (0.0, -h0) } else { (-h0, 0.0) };
        let iterations = pi.config.max_iterations.max(pj.config.max_iterations);
        let (mut ri, mut rj) = (r0i, r0j);
        for iteration in 1..=iterations {
            let mid = f64::midpoint(lo, hi);
            (ri, rj) = reactions(mid)?;
            if settled(ri, rj)? { return Ok((ri, rj, iteration)); }
            if mid == lo || mid == hi { break; }
            if finite(mid-(ri-rj))? < 0.0 { lo = mid; } else { hi = mid; }
        }
        let net = finite(ri-rj)?;
        let (expected, _) = law_force(&pi.law, pi.old, finite(f_i-a*net)?, dt)?;
        Err(ModalCouplingError::ContactSolve { residual_n: finite(ri-expected)?,
            tolerance_n: share*force_tolerance(ri, expected, pi.config)?, iterations })
    }
    // Only the joint friction owner reads these current trial reactions. An
    // accepted public report is deliberately not used as a one-sample delay.
    pub(super) fn normal_reaction(&self, index: usize) -> Option<f64> {
        self.reactions.get(index).copied()
    }
    // Whether the last stage refused at an actual endpoint; reading clears it.
    pub(super) fn take_endpoint_miss(&mut self) -> bool { std::mem::take(&mut self.endpoint_missed) }
    pub(super) fn set_headroom(&mut self, headroom: bool) { self.headroom = headroom; }
    fn new(system: &MultiFrictionModalSystem, contacts: Vec<(ModalContact, ModalContactConfig)>,
        limits: FrictionObstacleConfig, gate: &CancelGate) -> Result<Self, ModalCouplingError>
    {
        let q = contacts.len(); let p = system.points.len(); let s = p+q;
        if !(1..=32).contains(&limits.max_contacts) || q > limits.max_contacts || s > 64 {
            return Err(invalid("normal/friction port count exceeds the declared joint dimensions"));
        }
        let network = &system.core.network; let n = network.mode_count(); let k = network.columns.len();
        let terms = n.checked_mul(k+2).and_then(|v| v.checked_add((k+1)*(k+1)))
            .and_then(|v| v.checked_mul(s))
            .and_then(|v| n.checked_mul(s*s).and_then(|w| v.checked_add(w)))
            .ok_or_else(|| invalid("normal/friction setup work overflows"))?;
        if terms > limits.max_setup_terms { return Err(invalid("normal/friction setup exceeds max_setup_terms")); }
        let mut points = Vec::with_capacity(q);
        for (contact, config) in contacts {
            poll(Some(gate))?;
            let column = contact_column(network, &contact, config)?;
            let old = extension(network.components(), &column, 0.0)?;
            let law = OpeningContactStep::new(&contact.law, -old).map_err(ModalCouplingError::ContactLaw)?;
            let storage = ContactStorage::new(Box::new(ZeroStorage), 1, vec![contact.law.clone()])
                .map_err(ModalCouplingError::ContactLaw)?;
            let point = NormalPoint { contact, source: None, config, column, storage, law, old, free: old, iterations: 0 };
            point.penetration(old)?; point.energy(old)?;
            points.push(point);
        }
        let mut nn = vec![0.0; q*q]; let mut nt = vec![0.0; q*p]; let mut tn = vec![0.0; p*q];
        for j in 0..q {
            let response = network_response(network, &points[j].column, gate)?;
            for i in 0..q { poll(Some(gate))?; nn[i*q+j] = dot(&points[i].column, &response)?; }
            for i in 0..p { poll(Some(gate))?; tn[i*q+j] = dot(&system.points[i].column, &response)?; }
            if nn[j*q+j] <= 0.0 { return Err(invalid("normal obstacle needs positive network compliance")); }
        }
        for j in 0..p {
            let response = network_response(network, &system.points[j].column, gate)?;
            for i in 0..q { poll(Some(gate))?; nt[i*p+j] = dot(&points[i].column, &response)?; }
        }
        Ok(Self { points, limits, nn, nt, tn, reactions: vec![0.0; q], friction_ports: p,
            staged: vec![FrictionObstacleFrame::default(); q], accepted: vec![FrictionObstacleFrame::default(); q],
            motion: None, geometry: None, absolute_work_j: 0.0, headroom: false, endpoint_missed: false,
            residual: None })
    }
    pub(super) fn energy(&self, models: &[ModalAcousticTimeModel]) -> Result<f64, ModalCouplingError> {
        let mut total = 0.0;
        let accepted = if self.residual.is_some() { self.accepted_reactions() } else { Vec::new() };
        for (i, point) in self.points.iter().enumerate() {
            let mut x = extension(models, &point.column, 0.0)?;
            if let Some(motion) = &self.motion { x = finite(x-motion.accepted[i])?; }
            if self.residual.is_some() {
                // The omitted modes still hold the accepted reactions' static deflection.
                let lag = self.lag(i, &accepted)?;
                x = finite(x-lag)?;
                total = finite(total+0.5*accepted[i]*lag)?;
            }
            total = finite(total+point.energy(x)?)?;
        }
        Ok(total)
    }
    // Reuse the original stationary potential and inverse stiffness. Dynamic
    // secants/compliances are not static constitutive laws. No new root solver,
    // artificial settling interval or discarded startup vibration is introduced.
    pub(super) fn initialize_equilibrium(&self,
        network: &mut super::super::super::super::CoupledModalSystem,
        external: &[f64], max_sweeps: usize, gate: &CancelGate,
    ) -> Result<f64, ModalCouplingError> {
        poll(Some(gate))?;
        if self.residual.is_some() {
            // The equilibrium owner condenses the retained network only.
            return Err(invalid("released static preload does not include the obstacles' residual compliance"));
        }
        let contacts: Vec<_> = self.points.iter()
            .map(|point| (point.contact.clone(), point.config)).collect();
        let config = super::super::super::multiple::MultiContactConfig {
            max_contacts: self.limits.max_contacts, max_sweeps,
            max_setup_terms: self.limits.max_setup_terms,
        };
        // The original owner publishes only after checking all contacts, total
        // stored energy and actual force balance. Dynamic prepare() reconstructs
        // its disposable secants from the newly accepted state on the next step.
        // No dynamic report or sample is fabricated by this initial condition.
        network.initialize_contact_equilibrium(external, &contacts, config, gate)
    }
    pub(super) fn prepare(&mut self, old_q: &[f64], free_delta: &[f64], offsets: Option<&[f64]>, gate: Option<&CancelGate>)
        -> Result<(), ModalCouplingError>
    {
        // Closure lag through omitted modes, held from the accepted reactions.
        let lags = match &self.residual {
            Some(_) => {
                let accepted = self.accepted_reactions();
                (0..self.points.len()).map(|i| self.lag(i, &accepted)).collect::<Result<Vec<_>, _>>()?
            },
            None => Vec::new(),
        };
        self.reactions.fill(0.0);
        if let Some(motion) = &mut self.motion {
            motion.candidate.copy_from_slice(offsets.unwrap_or(&motion.accepted));
        }
        for (i, point) in self.points.iter_mut().enumerate() {
            poll(gate)?;
            point.old = dot(&point.column, old_q)?;
            if !lags.is_empty() {
                // The free endpoint holds no reaction yet, so no lag; the solve adds G R.
                point.free = finite(point.old+dot(&point.column, free_delta)?)?;
                point.old = finite(point.old-lags[i])?;
                point.law = OpeningContactStep::new(&point.contact.law, -point.old).map_err(ModalCouplingError::ContactLaw)?;
                point.iterations = 0;
                continue;
            }
            if let Some(geometry) = &mut self.geometry {
                if geometry.moving[i] {
                    point.free = finite(dot(&geometry.end_columns[i], old_q)?
                        + dot(&geometry.end_columns[i], free_delta)?)?;
                    geometry.start_transport[i] = dot(&geometry.change_columns[i], old_q)?;
                } else {
                    point.free = finite(point.old+dot(&point.column, free_delta)?)?;
                }
            } else {
                point.free = finite(point.old+dot(&point.column, free_delta)?)?;
            }
            if let Some(motion) = &self.motion {
                point.old = finite(point.old-motion.accepted[i])?;
                point.free = finite(point.free-motion.candidate[i])?;
            }
            point.law = OpeningContactStep::new(&point.contact.law, -point.old).map_err(ModalCouplingError::ContactLaw)?;
            point.iterations = 0;
        }
        Ok(())
    }
    pub(super) fn sweep(&mut self, friction: &[f64], dt: f64, gate: Option<&CancelGate>)
        -> Result<(), ModalCouplingError>
    {
        let q = self.points.len(); let p = self.friction_ports;
        if let Some(residual) = &self.residual {
            for &block in &residual.blocks {
                poll(gate)?;
                let (i, j) = match block { Block::Pair(i, j) => (i, j), Block::Single(i) => (i, i) };
                // Free closure under friction and every contact OUTSIDE this block.
                let mut frees = [0.0; 2];
                for (slot, k) in [i, j].into_iter().enumerate() {
                    let mut free = finite(self.points[k].free+dot(&self.nt[k*p..(k+1)*p], friction)?)?;
                    for m in 0..q { if m != i && m != j { free = finite(free-self.nn[k*q+m]*self.reactions[m])?; } }
                    frees[slot] = free;
                }
                if i == j {
                    let point = &mut self.points[i];
                    let mut config = point.config;
                    if self.headroom {
                        if config.force_absolute_tolerance_n*0.125 > 0.0 { config.force_absolute_tolerance_n *= 0.125; }
                        if config.force_relative_tolerance*0.125 > 0.0 { config.force_relative_tolerance *= 0.125; }
                    }
                    let (reaction, iterations) = solve_contact(&point.law, point.old, frees[0], self.nn[i*q+i], dt, config, gate)?;
                    self.reactions[i] = reaction; point.iterations += iterations;
                } else {
                    let (ri, rj, iterations) = self.solve_pair(i, j, frees[0], frees[1], dt, gate)?;
                    for (k, r) in [(i, ri), (j, rj)] {
                        let limit = self.points[k].config.maximum_force_n;
                        if r > limit {
                            return Err(ModalCouplingError::Budget { what: "contact normal force", value: r, limit });
                        }
                        self.reactions[k] = r; self.points[k].iterations += iterations;
                    }
                }
            }
            return Ok(());
        }
        for i in 0..q {
            poll(gate)?;
            let mut free = finite(self.points[i].free+dot(&self.nt[i*p..(i+1)*p], friction)?)?;
            for j in 0..q { if j != i { free = finite(free-self.nn[i*q+j]*self.reactions[j])?; } }
            let point = &mut self.points[i];
            let mut config = point.config;
            if self.headroom {
                // As in the plain multi-contact owner: a tighter coordinate root
                // leaves room, and acceptance keeps the ORIGINAL caps.
                if config.force_absolute_tolerance_n*0.125 > 0.0 { config.force_absolute_tolerance_n *= 0.125; }
                if config.force_relative_tolerance*0.125 > 0.0 { config.force_relative_tolerance *= 0.125; }
            }
            let (reaction, iterations) = solve_contact(&point.law, point.old, free, self.nn[i*q+i], dt, config, gate)?;
            self.reactions[i] = reaction; point.iterations += iterations;
        }
        Ok(())
    }
    pub(super) fn slip_shift(&self, port: usize) -> Result<f64, ModalCouplingError> {
        let q = self.points.len(); dot(&self.tn[port*q..(port+1)*q], &self.reactions)
    }
    pub(super) fn check_prediction(&self, friction: &[f64], dt: f64, gate: Option<&CancelGate>)
        -> Result<(), ModalCouplingError>
    {
        let q = self.points.len(); let p = self.friction_ports;
        for (i, point) in self.points.iter().enumerate() {
            poll(gate)?;
            let x = finite(finite(point.free+dot(&self.nt[i*p..(i+1)*p], friction)?)?
                - dot(&self.nn[i*q..(i+1)*q], &self.reactions)?)?;
            point.check(self.reactions[i], x, dt, if self.headroom { 0.5 } else { 1.0 })?;
        }
        Ok(())
    }
    pub(super) fn add_forces(&self, forces: &mut [f64], gate: Option<&CancelGate>) -> Result<(), ModalCouplingError> {
        for (i, (point, &reaction)) in self.points.iter().zip(&self.reactions).enumerate() {
            poll(gate)?;
            let column = self.geometry.as_ref().filter(|g| g.moving[i])
                .map_or(&point.column, |g| &g.force_columns[i]);
            if reaction != 0.0 { for (f, b) in forces.iter_mut().zip(column) { *f = finite(*f-b*reaction)?; } }
        }
        Ok(())
    }
    pub(super) fn stage(&mut self, models: &[ModalAcousticTimeModel], dt: f64, gate: Option<&CancelGate>)
        -> Result<(f64, f64, f64, f64), ModalCouplingError>
    {
        let (mut energy, mut loss) = (0.0, 0.0);
        let (mut work, mut absolute_work) = (0.0, 0.0);
        // Change of the held reactions over this step, for the declared numerical loss.
        let change = match &self.residual {
            Some(_) => self.reactions.iter().zip(&self.accepted).map(|(r, f)| finite(r-f.normal_force_n))
                .collect::<Result<Vec<_>, _>>()?,
            None => Vec::new(),
        };
        for (i, point) in self.points.iter().enumerate() {
            poll(gate)?;
            let (residual_energy_j, residual_loss_j, lag) = if self.residual.is_some() {
                let lag = self.lag(i, &self.reactions)?;
                (finite(0.5*self.reactions[i]*lag)?, finite(0.5*change[i]*self.lag(i, &change)?)?, lag)
            } else { (0.0, 0.0, 0.0) };
            let (mut x, attachment_work_j) = if let Some(geometry) = &self.geometry {
                if geometry.moving[i] {
                    let transport = finite(0.5 * finite(geometry.start_transport[i]
                        + extension(models, &geometry.change_columns[i], 0.0)?)?)?;
                    (extension(models, &geometry.end_columns[i], 0.0)?, finite(self.reactions[i]*transport)?)
                } else { (extension(models, &point.column, 0.0)?, 0.0) }
            } else { (extension(models, &point.column, 0.0)?, 0.0) };
            if self.residual.is_some() { x = finite(x-lag)?; }
            let (support_offset_m, translation_work_j) = if let Some(motion) = &self.motion {
                x = finite(x-motion.candidate[i])?;
                let distance = finite(motion.candidate[i]-motion.accepted[i])?;
                (motion.candidate[i], finite(-self.reactions[i]*distance)?)
            } else { (0.0, 0.0) };
            let actuator_work_j = if attachment_work_j == 0.0 { translation_work_j }
                else { finite(translation_work_j+attachment_work_j)? };
            let penetration_m = point.penetration(x)?;
            let (dissipation_j, constitutive_residual_n, force_tolerance_n) = match point.check(self.reactions[i], x, dt, 1.0) {
                Ok(checked) => checked,
                Err(error) => {
                    self.endpoint_missed = matches!(error, ModalCouplingError::ContactSolve { .. });
                    return Err(error);
                },
            };
            let stored_energy_j = point.energy(x)?;
            energy = finite(energy+stored_energy_j)?; loss = finite(loss+dissipation_j)?;
            work = finite(work+actuator_work_j)?;
            absolute_work = finite(absolute_work+translation_work_j.abs()+attachment_work_j.abs())?;
            self.staged[i] = FrictionObstacleFrame { normal_force_n: self.reactions[i], stored_energy_j,
                dissipation_j, penetration_m, constitutive_residual_n, force_tolerance_n,
                support_offset_m, attachment_work_j, actuator_work_j, residual_energy_j, residual_loss_j };
        }
        self.absolute_work_j = absolute_work;
        Ok((energy, loss, work, absolute_work))
    }
    // Omitted-mode storage and declared numerical loss of the staged step [J].
    pub(super) fn staged_residual_j(&self) -> Result<(f64, f64), ModalCouplingError> {
        let mut sums = (0.0, 0.0);
        for frame in &self.staged {
            sums = (finite(sums.0+frame.residual_energy_j)?, finite(sums.1+frame.residual_loss_j)?);
        }
        Ok(sums)
    }
    // Run only on the selected COMPLETE mechanical/thermal candidate, never
    // inside the coordinate roots, outer strain iteration or thermal Jacobian.
    pub(super) fn validate_sources(&self, models: &[ModalAcousticTimeModel], dt: f64,
        gate: Option<&CancelGate>) -> Result<(), ModalCouplingError> {
        for (i, point) in self.points.iter().enumerate() {
            if point.source.is_none() { continue; }
            poll(gate)?;
            let column = self.geometry.as_ref().filter(|g| g.moving[i])
                .map_or(&point.column, |g| &g.end_columns[i]);
            let mut x = extension(models, column, 0.0)?;
            if self.residual.is_some() { x = finite(x-self.lag(i, &self.reactions)?)?; }
            if let Some(motion) = &self.motion { x = finite(x-motion.candidate[i])?; }
            point.check_source(x, dt)?;
        }
        Ok(())
    }
    pub(super) fn publish(&mut self) {
        self.accepted.copy_from_slice(&self.staged);
        if let Some(motion) = &mut self.motion { motion.accepted.copy_from_slice(&motion.candidate); }
        if let Some(geometry) = &self.geometry {
            for (point, endpoint) in self.points.iter_mut().zip(&geometry.end_columns) {
                point.column.copy_from_slice(endpoint);
            }
        }
    }
    // Scale from individual physical supplies, not their possibly cancelling
    // net sum. Read only staged reports from this SAME joint contact solve.
    pub(super) fn staged_absolute_work_j(&self) -> Result<f64, ModalCouplingError> {
        finite(self.absolute_work_j)
    }
}

#[cfg(test)]
mod tests;
