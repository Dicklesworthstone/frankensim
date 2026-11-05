//! Simultaneous source-owned adhesion on one complete native modal network.
//!
//! Every approach obeys `delta = delta_free - S F`, with the full original
//! reciprocal compliance `S_ij = B_i^T D B_j`. Deterministic nonlinear sweeps
//! reuse the original source-normal scalar solve, including its branch and
//! signed-work rules. They iterate one joint physical step; they never publish
//! successive independent collisions. Redundant work columns are permitted.
//!
//! Each accepted source is reevaluated on the actual native candidate motion.
//! Original modal, connection, relaxation and radiation states remain in the
//! network, and all contact potentials share its energy account. Scratch and
//! reports are allocated at setup; the contact coupling does not allocate while
//! stepping. The original source evaluators retain their own work budgets.
//!
//! A conservative coupled negative-curvature screen requires
//! `||sqrt(D-) S sqrt(D-)||_infinity < 1`, where `D-` contains the negative
//! parts of the original path slopes. Positive source stiffness is ignored by
//! this sufficient screen, so some resolvable steps require refinement. It is
//! a local floating-point stability check, not a global uniqueness certificate.
//! Attached JKR/Maugis folds, nonconvergence and out-of-domain paths refuse;
//! no snap, contact discovery, evolving frame or shear law is invented.

/// Atomic original-source temperature changes on the shared native network.
pub mod temperature;

use super::{
    BoundNormalContact, CancelGate, CoupledModalSystem, ModalAcousticTimeModel, NormalLawSource,
    NormalLawSourceError, SourceNormalModalConfig, SourceNormalModalError, admit_contact,
    check_contact_work, check_coordinate, check_native_energy, dot, extension, finite,
    force_tolerance, invalid, limit, mean_force, normal_velocity, observation, poll, solve,
};
use crate::render::schedule::force::coupled::contact::multiple::{
    MAX_NORMAL_CONTACTS, MultiContactConfig,
};
use fs_math::{dd::Dd, det};

/// One original source on the common accepted native endpoints.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SourceNormalContactFrame {
    /// Signed path effort actually applied through this contact's complete row.
    pub normal_force_n: f64,
    pub endpoint_normal_force_n: f64,
    pub approach_before_m: f64,
    pub approach_after_m: f64,
    /// Original signed absolute potential; attraction is not dissipative heat.
    pub contact_energy_j: f64,
    pub source_work_j: f64,
    pub constitutive_residual_n: f64,
    pub force_tolerance_n: f64,
    pub coordinate_residual_m: f64,
    /// Original potential change minus its own path work.
    pub contact_work_residual_j: f64,
    pub contact_work_tolerance_j: f64,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct MultipleSourceNormalModalFrame {
    pub sample: u64,
    pub observer_pressure_pa: f64,
    pub network_energy_j: f64,
    pub contact_energy_j: f64,
    pub external_work_j: f64,
    /// Work of the represented sum of contact forces on actual native motion.
    pub contact_work_on_network_j: f64,
    pub network_dissipation_j: f64,
    pub energy_residual_j: f64,
    pub energy_tolerance_j: f64,
    /// Sum of source work plus the represented native contact work.
    pub contact_work_residual_j: f64,
    pub contact_work_tolerance_j: f64,
    pub sweeps: usize,
    pub scalar_iterations: usize,
    pub contacts: Vec<SourceNormalContactFrame>,
}

struct Point {
    contact: BoundNormalContact,
    config: SourceNormalModalConfig,
}
impl Point {
    fn approach(&self, models: &[ModalAcousticTimeModel]) -> Result<f64, SourceNormalModalError> {
        let value = extension(models, self.contact.work_column(), self.contact.gap_m())?;
        self.check_approach(value)?;
        Ok(value)
    }
    fn check_approach(&self, value: f64) -> Result<(), SourceNormalModalError> {
        if !value.is_finite()
            || value < self.config.minimum_approach_m
            || value > self.config.maximum_approach_m
        {
            return Err(SourceNormalModalError::Unresolved(
                "simultaneous source approach leaves its admitted domain",
            ));
        }
        Ok(())
    }
    fn kinematics(&self, approach: f64, rate: f64, dt: f64) -> Result<(), SourceNormalModalError> {
        if let NormalLawSource::HertzPoint(source) = self.contact.source() {
            source
                .validate_kinematics(approach.max(0.0), rate, dt)
                .map_err(NormalLawSourceError::Point)?;
        }
        Ok(())
    }
}

/// Several signed conservative contact sources sharing one retained network.
/// Failed proposals change only scratch, never accepted motion or observations.
pub struct MultipleSourceNormalModalSystem {
    network: CoupledModalSystem,
    points: Vec<Point>,
    config: MultiContactConfig,
    compliance: Vec<f64>,
    forces: Vec<f64>,
    before: Vec<f64>,
    free_increment: Vec<f64>,
    reactions: Vec<f64>,
    solved_approaches: Vec<f64>,
    trial_approaches: Vec<f64>,
    negative_slopes: Vec<f64>,
    frame: MultipleSourceNormalModalFrame,
    staged: MultipleSourceNormalModalFrame,
}

impl MultipleSourceNormalModalSystem {
    /// Bind each original law and its per-contact domain/accuracy limits.
    /// `MultiContactConfig` supplies the same contact-count, joint-sweep and
    /// setup-term ceilings as the existing nonadhesive native contact owner.
    /// Mixed Hertz, finite-gap, DMT and explicitly attached JKR/Maugis sources
    /// are supported. A dissipative point law cannot enter this conservative
    /// path; its loss must be advanced by the original timed-loss owner.
    pub fn new(
        network: CoupledModalSystem,
        contacts: Vec<(BoundNormalContact, SourceNormalModalConfig)>,
        config: MultiContactConfig,
        gate: &CancelGate,
    ) -> Result<Self, SourceNormalModalError> {
        poll(Some(gate))?;
        let p = contacts.len();
        let n = network.mode_count();
        let k = network.columns.len();
        if !(1..=MAX_NORMAL_CONTACTS).contains(&config.max_contacts)
            || p == 0
            || p > config.max_contacts
            || !(1..=128).contains(&config.max_sweeps)
        {
            return Err(
                invalid("bounded nonempty source contact set and sweep budget required").into(),
            );
        }
        let terms = k
            .checked_add(1)
            .and_then(|v| v.checked_mul(v))
            .and_then(|square| n.checked_mul(k + 2).and_then(|v| v.checked_add(square)))
            .and_then(|v| v.checked_mul(p))
            .and_then(|v| n.checked_mul(p * p).and_then(|w| v.checked_add(w)))
            .ok_or_else(|| invalid("simultaneous source setup work overflow"))?;
        if terms > config.max_setup_terms {
            return Err(invalid("simultaneous source setup exceeds max_setup_terms").into());
        }
        let mut points = Vec::with_capacity(p);
        for (contact, limits) in contacts {
            poll(Some(gate))?;
            admit_contact(&network, &contact, limits)?;
            let point = Point {
                contact,
                config: limits,
            };
            let approach = point.approach(network.components())?;
            let (_, force) = observation(point.contact.source(), approach)?;
            limit(
                "initial simultaneous source force",
                force.abs(),
                limits.maximum_force_n,
            )?;
            point.kinematics(
                approach,
                normal_velocity(network.components(), point.contact.work_column())?,
                network.dt,
            )?;
            points.push(point);
        }
        let mut compliance = vec![0.0; p * p];
        for j in 0..p {
            poll(Some(gate))?;
            let response = network.linear_response(points[j].contact.work_column(), gate)?;
            for i in 0..p {
                compliance[i * p + j] = dot(points[i].contact.work_column(), &response)?;
            }
            if compliance[j * p + j] <= 0.0 {
                return Err(invalid("each source contact needs positive native compliance").into());
            }
        }
        let empty_frame = || MultipleSourceNormalModalFrame {
            contacts: vec![SourceNormalContactFrame::default(); p],
            ..MultipleSourceNormalModalFrame::default()
        };
        let owner = Self {
            network,
            points,
            config,
            compliance,
            forces: vec![0.0; n],
            before: vec![0.0; p],
            free_increment: vec![0.0; p],
            reactions: vec![0.0; p],
            solved_approaches: vec![0.0; p],
            trial_approaches: vec![0.0; p],
            negative_slopes: vec![0.0; p],
            frame: empty_frame(),
            staged: empty_frame(),
        };
        limit(
            "initial simultaneous source energy magnitude",
            owner.total_energy_j()?.abs(),
            owner.network.config.maximum_total_energy_j,
        )?;
        poll(Some(gate))?;
        Ok(owner)
    }

    pub const fn network(&self) -> &CoupledModalSystem {
        &self.network
    }
    pub fn components(&self) -> &[ModalAcousticTimeModel] {
        self.network.components()
    }
    pub fn mode_count(&self) -> usize {
        self.network.mode_count()
    }
    pub fn samples_rendered(&self) -> u64 {
        self.network.samples_rendered()
    }
    pub fn sample_period_s(&self) -> f64 {
        self.network.sample_period_s()
    }
    pub fn contact_count(&self) -> usize {
        self.points.len()
    }
    pub fn contact(&self, index: usize) -> Option<&BoundNormalContact> {
        self.points.get(index).map(|point| &point.contact)
    }
    pub const fn config(&self) -> MultiContactConfig {
        self.config
    }
    /// Original complete response matrix in contact-major row order [m/N].
    pub fn compliance_m_per_n(&self) -> &[f64] {
        &self.compliance
    }
    pub fn last_frame(&self) -> Option<&MultipleSourceNormalModalFrame> {
        (self.frame.sample != 0).then_some(&self.frame)
    }
    pub fn total_energy_j(&self) -> Result<f64, SourceNormalModalError> {
        let mut energy = Dd::from_f64(self.network.total_energy_j()?);
        for point in &self.points {
            energy = energy
                + Dd::from_f64(
                    observation(point.contact.source(), point.approach(self.components())?)?.0,
                );
        }
        Ok(finite(energy.to_f64())?)
    }

    pub fn step(
        &mut self,
        external: &[f64],
    ) -> Result<&MultipleSourceNormalModalFrame, SourceNormalModalError> {
        self.step_under_gate(external, &CancelGate::new_clock_free())
    }
    pub fn step_under_gate(
        &mut self,
        external: &[f64],
        gate: &CancelGate,
    ) -> Result<&MultipleSourceNormalModalFrame, SourceNormalModalError> {
        self.stage(external, gate)?;
        self.publish();
        Ok(&self.frame)
    }
    /// A receiving-state refusal or late cancellation preserves all accepted
    /// native histories and the previous frame. Install the returned receiving
    /// candidate only after success, by an infallible move.
    pub fn step_with_admission<R, E>(
        &mut self,
        external: &[f64],
        gate: &CancelGate,
        admit: impl FnOnce(&MultipleSourceNormalModalFrame, &[ModalAcousticTimeModel]) -> Result<R, E>,
    ) -> Result<R, E>
    where
        E: From<SourceNormalModalError>,
    {
        self.stage(external, gate).map_err(E::from)?;
        let candidate = admit(&self.staged, &self.network.candidates)?;
        poll(Some(gate))
            .map_err(SourceNormalModalError::from)
            .map_err(E::from)?;
        self.publish();
        Ok(candidate)
    }
    fn publish(&mut self) {
        self.network.publish_staged();
        core::mem::swap(&mut self.frame, &mut self.staged);
    }

    fn stage(&mut self, external: &[f64], gate: &CancelGate) -> Result<(), SourceNormalModalError> {
        poll(Some(gate))?;
        let network_before = self.network.total_energy_j()?;
        for (i, point) in self.points.iter().enumerate() {
            self.before[i] = point.approach(self.network.components())?;
            point.kinematics(
                self.before[i],
                normal_velocity(self.network.components(), point.contact.work_column())?,
                self.network.dt,
            )?;
            self.reactions[i] = observation(point.contact.source(), self.before[i])?.1;
            self.trial_approaches[i] = self.before[i];
        }
        self.coupled_stability(gate)?;
        self.network.prepare_prediction(external, Some(gate))?;
        for (i, point) in self.points.iter().enumerate() {
            self.free_increment[i] = dot(point.contact.work_column(), &self.network.free_delta)?;
        }
        let (sweeps, iterations) = self.solve(gate)?;
        self.forces.copy_from_slice(external);
        for (point, &reaction) in self.points.iter().zip(&self.reactions) {
            for (force, &weight) in self.forces.iter_mut().zip(point.contact.work_column()) {
                *force = finite(*force - weight * reaction)?;
            }
        }
        self.network.stage_inner(&self.forces, Some(gate))?;
        let mut stored = Dd::ZERO;
        let mut storage_change = Dd::ZERO;
        let mut source_work = Dd::ZERO;
        let mut contact_tolerance = Dd::ZERO;
        let mut contact_scale = 0.0_f64;
        for (i, point) in self.points.iter().enumerate() {
            poll(Some(gate))?;
            let before = self.before[i];
            let after = point.approach(&self.network.candidates)?;
            let increment = finite(after - before)?;
            let coordinate_residual =
                check_coordinate(after, self.trial_approaches[i], point.config)?;
            check_coordinate(
                increment,
                dot(point.contact.work_column(), &self.network.free_delta)?,
                point.config,
            )?;
            let (old_energy, _) = observation(point.contact.source(), before)?;
            let (energy, endpoint) = observation(point.contact.source(), after)?;
            let expected = mean_force(point.contact.source(), before, after)?;
            let force = self.reactions[i];
            let residual = finite(force - expected)?;
            let force_tolerance = force_tolerance(force, expected, point.config)?;
            if residual.abs() > force_tolerance {
                return Err(super::ModalCouplingError::ContactSolve {
                    residual_n: residual,
                    tolerance_n: force_tolerance,
                    iterations,
                }
                .into());
            }
            for value in [force, endpoint] {
                limit(
                    "simultaneous sourced normal force",
                    value.abs(),
                    point.config.maximum_force_n,
                )?;
            }
            point.kinematics(
                after,
                normal_velocity(&self.network.candidates, point.contact.work_column())?,
                self.network.dt,
            )?;
            point.kinematics(
                before.max(after),
                finite(increment / self.network.dt)?,
                self.network.dt,
            )?;
            let work = finite(expected * increment)?;
            let change = finite(energy - old_energy)?;
            let scale = old_energy.abs().max(energy.abs()).max(work.abs());
            let tolerance = finite(
                point.config.contact_energy_absolute_tolerance_j
                    + point.config.contact_energy_relative_tolerance * scale,
            )?;
            check_contact_work(finite(change - work)?, tolerance)?;
            stored = stored + Dd::from_f64(energy);
            storage_change = storage_change + Dd::from_f64(change);
            source_work = source_work + Dd::from_f64(work);
            contact_tolerance = contact_tolerance + Dd::from_f64(tolerance);
            contact_scale = contact_scale.max(scale);
            self.staged.contacts[i] = SourceNormalContactFrame {
                normal_force_n: force,
                endpoint_normal_force_n: endpoint,
                approach_before_m: before,
                approach_after_m: after,
                contact_energy_j: energy,
                source_work_j: work,
                constitutive_residual_n: residual,
                force_tolerance_n: force_tolerance,
                coordinate_residual_m: coordinate_residual,
                contact_work_residual_j: change - work,
                contact_work_tolerance_j: tolerance,
            };
            self.trial_approaches[i] = after;
        }
        self.coupled_stability(gate)?;
        // Independently measure the actual represented force correction. This
        // catches a source contribution erased while summing native forces.
        let mut actual_work = Dd::ZERO;
        for ((held, input), delta) in self
            .forces
            .iter()
            .zip(external)
            .zip(&self.network.free_delta)
        {
            actual_work = actual_work + Dd::from_f64(finite((held - input) * delta)?);
        }
        let work_tolerance = finite(contact_tolerance.to_f64())?;
        let work_residual = finite((actual_work + source_work).to_f64())?;
        check_contact_work(work_residual, work_tolerance)?;
        check_contact_work(
            finite((actual_work + storage_change).to_f64())?,
            work_tolerance,
        )?;
        let frame = &self.network.staged_frame;
        let network_after = finite(frame.total_energy_j())?;
        let network_loss = finite(frame.total_dissipation_j())?;
        let external_work = dot(external, &self.network.free_delta)?;
        let contact_energy = finite(stored.to_f64())?;
        limit(
            "simultaneous native/source energy magnitude",
            finite(network_after + contact_energy)?.abs(),
            self.network.config.maximum_total_energy_j,
        )?;
        let residual = finite(
            (Dd::from_f64(network_after) - Dd::from_f64(network_before)
                + storage_change
                + Dd::from_f64(network_loss)
                - Dd::from_f64(external_work))
            .to_f64(),
        )?;
        let scale = network_before
            .abs()
            .max(network_after.abs())
            .max(contact_scale)
            .max(external_work.abs())
            .max(network_loss.abs());
        let tolerance = finite(
            self.network.config.energy_absolute_tolerance_j
                + self.network.config.energy_relative_tolerance * scale,
        )?;
        check_native_energy(residual, tolerance)?;
        check_native_energy(
            finite(
                (Dd::from_f64(frame.external_work_j) - Dd::from_f64(external_work) - actual_work)
                    .to_f64(),
            )?,
            tolerance,
        )?;
        self.staged.sample = frame.sample;
        self.staged.observer_pressure_pa = frame.observer_pressure_pa;
        self.staged.network_energy_j = network_after;
        self.staged.contact_energy_j = contact_energy;
        self.staged.external_work_j = external_work;
        self.staged.contact_work_on_network_j = finite(actual_work.to_f64())?;
        self.staged.network_dissipation_j = network_loss;
        self.staged.energy_residual_j = residual;
        self.staged.energy_tolerance_j = tolerance;
        self.staged.contact_work_residual_j = work_residual;
        self.staged.contact_work_tolerance_j = work_tolerance;
        self.staged.sweeps = sweeps;
        self.staged.scalar_iterations = iterations;
        self.network.validate_staged_material_points(Some(gate))?;
        poll(Some(gate))?;
        Ok(())
    }

    fn solve(&mut self, gate: &CancelGate) -> Result<(usize, usize), SourceNormalModalError> {
        let p = self.points.len();
        let mut iterations = 0;
        for sweep in 1..=self.config.max_sweeps {
            for (i, point) in self.points.iter().enumerate() {
                poll(Some(gate))?;
                let mut increment = Dd::from_f64(self.free_increment[i]);
                for j in 0..p {
                    if j != i {
                        increment = increment
                            - Dd::from_f64(finite(self.compliance[i * p + j] * self.reactions[j])?);
                    }
                }
                let solution = solve::solve(
                    point.contact.source(),
                    self.before[i],
                    finite(increment.to_f64())?,
                    self.compliance[i * p + i],
                    point.config,
                    gate,
                )?;
                self.reactions[i] = solution.force_n;
                self.solved_approaches[i] = solution.approach_m;
                iterations += solution.iterations;
            }
            let mut converged = true;
            for (i, point) in self.points.iter().enumerate() {
                let mut increment = Dd::from_f64(self.free_increment[i]);
                for j in 0..p {
                    increment = increment
                        - Dd::from_f64(finite(self.compliance[i * p + j] * self.reactions[j])?);
                }
                let to = finite(self.before[i] + increment.to_f64())?;
                point.check_approach(to)?;
                self.trial_approaches[i] = to;
                let expected = mean_force(point.contact.source(), self.before[i], to)?;
                let force_error = finite(self.reactions[i] - expected)?.abs();
                let force_tolerance = force_tolerance(self.reactions[i], expected, point.config)?;
                let scale = self.before[i]
                    .abs()
                    .max(to.abs())
                    .max(self.free_increment[i].abs())
                    .max((self.compliance[i * p + i] * self.reactions[i]).abs());
                let tolerance = (point.config.coordinate_absolute_tolerance_m
                    + point.config.coordinate_relative_tolerance * scale)
                    .min(32.0 * f64::EPSILON * scale);
                converged &= force_error <= force_tolerance
                    && (to - self.solved_approaches[i]).abs() <= tolerance;
            }
            self.coupled_stability(gate)?;
            if converged {
                return Ok((sweep, iterations));
            }
        }
        Err(SourceNormalModalError::Unresolved(
            "simultaneous source nonlinear sweep budget",
        ))
    }

    fn coupled_stability(&mut self, gate: &CancelGate) -> Result<(), SourceNormalModalError> {
        for (i, point) in self.points.iter().enumerate() {
            poll(Some(gate))?;
            self.negative_slopes[i] = negative_slope_root(
                point.contact.source(),
                point.config,
                self.before[i],
                self.trial_approaches[i],
            )?;
        }
        check_coupled_stability(&self.compliance, &self.negative_slopes, gate)
    }
}

// Shared by mechanical paths and source replacement. In particular, a refresh
// must check the final complete set, including untouched adhesive neighbours.
fn negative_slope_root(
    source: &NormalLawSource,
    config: SourceNormalModalConfig,
    before: f64,
    after: f64,
) -> Result<f64, SourceNormalModalError> {
    let slope = source.mean_force_end_derivative(before, after)?.ok_or(
        SourceNormalModalError::Unresolved("original source path derivative unavailable"),
    )?;
    let mut negative = (-finite(slope)?).max(0.0);
    // Across its finite tail DMT has a known global lower slope. Keep
    // that original bound even if the current point is compressive.
    if let NormalLawSource::DmtSphere(law) = source
        && config.minimum_approach_m < 0.0
        && config.maximum_approach_m > -law.interaction_range_m()
    {
        let tail = law
            .evaluate(0.0)
            .map_err(NormalLawSourceError::Dmt)?
            .one_sided_tangents_n_per_m[0];
        negative = negative.max(-0.5 * tail);
    }
    Ok(det::sqrt(negative))
}

fn check_coupled_stability(
    compliance: &[f64],
    negative_slopes: &[f64],
    gate: &CancelGate,
) -> Result<(), SourceNormalModalError> {
    let p = negative_slopes.len();
    for i in 0..p {
        poll(Some(gate))?;
        let mut row = 0.0;
        for j in 0..p {
            row = finite(
                row + negative_slopes[i] * compliance[i * p + j].abs() * negative_slopes[j],
            )?;
        }
        if row >= 1.0 - 64.0 * f64::EPSILON {
            return Err(SourceNormalModalError::Unresolved(
                "coupled adhesive negative curvature needs a smaller native step",
            ));
        }
    }
    Ok(())
}
