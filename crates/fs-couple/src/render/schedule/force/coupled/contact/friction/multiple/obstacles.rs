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

    /// Number of independent normal obstacles, not the prescribed bow loads.
    #[must_use]
    pub fn obstacle_count(&self) -> usize { self.obstacles.as_ref().map_or(0, |o| o.points.len()) }

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
}
impl ObstacleSet {
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
            let point = NormalPoint { contact, config, column, storage, law, old, free: old, iterations: 0 };
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
            motion: None, geometry: None, absolute_work_j: 0.0, headroom: false, endpoint_missed: false })
    }
    pub(super) fn energy(&self, models: &[ModalAcousticTimeModel]) -> Result<f64, ModalCouplingError> {
        let mut total = 0.0;
        for (i, point) in self.points.iter().enumerate() {
            let mut x = extension(models, &point.column, 0.0)?;
            if let Some(motion) = &self.motion { x = finite(x-motion.accepted[i])?; }
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
        self.reactions.fill(0.0);
        if let Some(motion) = &mut self.motion {
            motion.candidate.copy_from_slice(offsets.unwrap_or(&motion.accepted));
        }
        for (i, point) in self.points.iter_mut().enumerate() {
            poll(gate)?;
            point.old = dot(&point.column, old_q)?;
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
        for (i, point) in self.points.iter().enumerate() {
            poll(gate)?;
            let (mut x, attachment_work_j) = if let Some(geometry) = &self.geometry {
                if geometry.moving[i] {
                    let transport = finite(0.5 * finite(geometry.start_transport[i]
                        + extension(models, &geometry.change_columns[i], 0.0)?)?)?;
                    (extension(models, &geometry.end_columns[i], 0.0)?, finite(self.reactions[i]*transport)?)
                } else { (extension(models, &point.column, 0.0)?, 0.0) }
            } else { (extension(models, &point.column, 0.0)?, 0.0) };
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
                support_offset_m, attachment_work_j, actuator_work_j };
        }
        self.absolute_work_j = absolute_work;
        Ok((energy, loss, work, absolute_work))
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
