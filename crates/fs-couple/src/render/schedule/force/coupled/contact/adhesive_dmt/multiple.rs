//! Simultaneous original DMT normal and two-direction junction shear contacts.
//!
//! Every contact sees the complete held-force response `x = x_free - S F`,
//! with `S_ij = B_i^T D_native B_j` in contact/normal/shear order. Bounded
//! deterministic block sweeps reuse each original three-coordinate source
//! solve. They iterate ONE physical sample; they do not apply serial impacts.
//! The original native network advances once with the sum of all reactions.
//!
//! Original normal storage and junction loss are independently reevaluated on
//! actual candidate-minus-accepted native displacement. Each contact and the
//! represented total force correction pass small-work gates before all native
//! histories, contact occupancies and observations publish together. Adapter
//! scratch and report storage are allocated at construction and reused.
//!
//! Sources stay fixed during mechanical steps; [`temperature`] applies atomic
//! original-source changes at held coordinates. Axes and work columns stay fixed.
//! Straight coordinate paths and the original static/kinetic branch policy
//! retain their limitations: this is
//! not impact/event resolution, a uniqueness or global convergence theorem,
//! contact discovery, heat deposition, torsion or an implicit temperature solve.
//! Redundant columns BETWEEN contacts are allowed; every local three-axis
//! response must pass its original source's mechanical admission.
//! The collective curvature screen is sufficient and conservative on smooth
//! kinetic branches. It does not resolve a static-strength drop or prove
//! uniqueness or convergence of the bounded block iteration.

use super::{
    CoupledModalSystem, DmtAdhesivePath, DmtCoupledConfig, DmtCoupledError, DmtCoupledLoadLine,
    DmtCoupledRegime, DmtCoupledState, ModalAcousticTimeModel, ModalCouplingError,
    SourceAdhesiveShear, SourceDmtSphereContact, finite, invalid, limit, network_response, poll,
    tangential_velocity,
};
use fs_exec::CancelGate;
use fs_matdb::MaterialStateId;
use fs_math::{dd::Dd, det};

pub mod temperature;

/// Hard bound in addition to the caller's quadratic setup and evaluation caps.
pub const MAX_ADHESIVE_DMT_CONTACTS: usize = 128;

/// One original source and its complete fixed native work columns.
pub struct AdhesiveDmtModalContact<'a> {
    pub normal: &'a SourceDmtSphereContact,
    pub shear: &'a SourceAdhesiveShear,
    /// Actual receiving material order, checked against the original interface.
    pub surface_materials: [&'a MaterialStateId; 2],
    /// Physical coordinate datum at the supplied sample-zero native state.
    pub initial_position_m: [f64; 3],
    /// Normal, shear 1, shear 2; include every original native mode.
    pub columns: [Vec<f64>; 3],
    /// Original physical bounds and per-block nonlinear solve budget. Repeated
    /// block solves also consume the aggregate evaluation budget below.
    pub config: DmtCoupledConfig,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MultipleAdhesiveDmtModalConfig {
    pub maximum_contacts: usize,
    /// In 1..=128. Exhaustion refuses the complete physical sample.
    pub maximum_sweeps: usize,
    /// Bounds response construction and the complete 3N by 3N storage/work.
    pub maximum_setup_terms: usize,
    /// Original source-path probes across ALL contacts and sweeps, including
    /// per-step initial path probes and actual-motion replay. Constructor
    /// initialization belongs to cold setup. Each block's full remaining
    /// ceiling is reserved before entering it; unused reservations are returned.
    pub maximum_source_evaluations: usize,
    /// Enforced independently for each contact and summed for spatial transfer.
    /// Unrelated native storage or external work never enlarges these gates.
    pub contact_work_absolute_tolerance_j: f64,
    pub contact_work_relative_tolerance: f64,
}

#[derive(Debug)]
pub enum MultipleAdhesiveDmtModalError {
    Native(ModalCouplingError),
    Contact {
        contact: usize,
        source: DmtCoupledError,
    },
    /// Original individual sources passed, but their shared response did not.
    CoupledStability,
    EvaluationBudget,
    NotConverged {
        sweeps: usize,
    },
    Coordinate {
        contact: usize,
        axis: usize,
        residual_m: f64,
        tolerance_m: f64,
    },
    Work {
        contact: Option<usize>,
        residual_j: f64,
        tolerance_j: f64,
    },
}
impl From<ModalCouplingError> for MultipleAdhesiveDmtModalError {
    fn from(error: ModalCouplingError) -> Self {
        Self::Native(error)
    }
}
impl core::fmt::Display for MultipleAdhesiveDmtModalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "simultaneous native DMT adhesion refused: {self:?}")
    }
}
impl std::error::Error for MultipleAdhesiveDmtModalError {}

/// One original law on the SAME actual native endpoints as every other law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdhesiveDmtContactFrame {
    pub position_before_m: [f64; 3],
    pub position_after_m: [f64; 3],
    /// Positive efforts resist the corresponding coordinate motion.
    pub mean_force_n: [f64; 3],
    pub path: DmtAdhesivePath,
    pub regime: DmtCoupledRegime,
    pub friction_dissipation_j: f64,
    pub prescribed_work_j: f64,
    pub constitutive_residual_n: [f64; 3],
    /// Actual native position minus the complete condensed response.
    pub compatibility_residual_m: [f64; 3],
    pub compatibility_tolerance_m: f64,
    /// Applied local work minus original normal path work and physical shear loss.
    pub work_residual_j: f64,
    pub work_tolerance_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct MultipleAdhesiveDmtModalFrame {
    pub sample: u64,
    pub observer_pressure_pa: f64,
    pub network_energy_j: f64,
    pub surface_energy_j: f64,
    pub external_work_j: f64,
    pub prescribed_work_j: f64,
    pub contact_work_on_network_j: f64,
    pub network_dissipation_j: f64,
    pub friction_dissipation_j: f64,
    pub energy_residual_j: f64,
    pub energy_tolerance_j: f64,
    pub contact_work_residual_j: f64,
    pub contact_work_tolerance_j: f64,
    pub sweeps: usize,
    pub iterations: usize,
    pub evaluations: usize,
    /// Original construction order; prescribed rates use this same order.
    pub contacts: Vec<AdhesiveDmtContactFrame>,
}

struct Point<'a> {
    load_line: DmtCoupledLoadLine<'a>,
    columns: [Vec<f64>; 3],
    state: DmtCoupledState,
}

pub struct MultipleAdhesiveDmtModalSystem<'a> {
    network: CoupledModalSystem,
    points: Vec<Point<'a>>,
    config: MultipleAdhesiveDmtModalConfig,
    compliance: Vec<f64>,
    previous: Vec<DmtCoupledState>,
    candidate: Vec<DmtCoupledState>,
    before_energy: Vec<f64>,
    free: Vec<[f64; 3]>,
    reactions: Vec<[f64; 3]>,
    position: Vec<[f64; 3]>,
    trials: Vec<Option<fs_contact::normal_patch::dmt::coupled::DmtCoupledStep>>,
    forces: Vec<f64>,
    delta: Vec<f64>,
    frame: MultipleAdhesiveDmtModalFrame,
    staged: MultipleAdhesiveDmtModalFrame,
}

impl<'a> MultipleAdhesiveDmtModalSystem<'a> {
    pub fn new(
        network: CoupledModalSystem,
        contacts: Vec<AdhesiveDmtModalContact<'a>>,
        config: MultipleAdhesiveDmtModalConfig,
        gate: &CancelGate,
    ) -> Result<Self, MultipleAdhesiveDmtModalError> {
        poll(Some(gate))?;
        let p = contacts.len();
        let n = network.mode_count();
        let k = network.columns.len();
        if network.samples_rendered() != 0
            || !(1..=MAX_ADHESIVE_DMT_CONTACTS).contains(&config.maximum_contacts)
            || p == 0
            || p > config.maximum_contacts
            || !(1..=128).contains(&config.maximum_sweeps)
            || config.maximum_source_evaluations == 0
            || !config.contact_work_absolute_tolerance_j.is_finite()
            || config.contact_work_absolute_tolerance_j < 0.0
            || !config.contact_work_relative_tolerance.is_finite()
            || !(0.0..1.0).contains(&config.contact_work_relative_tolerance)
        {
            return Err(invalid(
                "bounded sample-zero DMT contact set and finite work controls required",
            )
            .into());
        }
        let rows = 3 * p;
        let terms = k
            .checked_add(1)
            .and_then(|x| x.checked_mul(x))
            .and_then(|x| n.checked_mul(k + 2).and_then(|v| x.checked_add(v)))
            .and_then(|x| x.checked_mul(rows))
            .and_then(|x| {
                rows.checked_mul(rows)
                    .and_then(|v| v.checked_mul(n + 1))
                    .and_then(|v| x.checked_add(v))
            })
            .ok_or_else(|| invalid("simultaneous DMT setup work overflow"))?;
        if terms > config.maximum_setup_terms {
            return Err(invalid("simultaneous DMT setup exceeds its term budget").into());
        }
        for contact in &contacts {
            poll(Some(gate))?;
            let materials = contact
                .normal
                .normal_state()
                .elastic()
                .interface()
                .surface_materials();
            if contact.surface_materials != [&materials[0], &materials[1]]
                || contact.config.maximum_evaluations < 2
                || contact
                    .columns
                    .iter()
                    .any(|v| v.len() != n || v.iter().any(|x| !x.is_finite()))
            {
                return Err(invalid("original ordered DMT surfaces, complete columns and a replay reservation required").into());
            }
        }
        let mut compliance = vec![0.0; rows * rows];
        for j in 0..rows {
            poll(Some(gate))?;
            let response = network_response(&network, &contacts[j / 3].columns[j % 3], gate)?;
            for i in 0..rows {
                compliance[i * rows + j] = dot(&contacts[i / 3].columns[i % 3], &response)?;
            }
        }
        let mut points = Vec::with_capacity(p);
        let mut states = Vec::with_capacity(p);
        let mut frames = Vec::with_capacity(p);
        for (index, contact) in contacts.into_iter().enumerate() {
            poll(Some(gate))?;
            let local = core::array::from_fn(|i| {
                core::array::from_fn(|j| compliance[(3 * index + i) * rows + 3 * index + j])
            });
            let mut controls = contact.config;
            controls.maximum_evaluations -= 1;
            let load_line = DmtCoupledLoadLine::new(contact.normal, contact.shear, local, controls)
                .map_err(|source| source_error(index, source))?;
            let velocity = tangential_velocity(network.components(), &contact.columns, [0.0; 3])?;
            let state = load_line
                .initialize_with_velocity(contact.initial_position_m, velocity)
                .map_err(|source| source_error(index, source))?;
            let path = load_line
                .path(contact.initial_position_m[0], contact.initial_position_m[0])
                .map_err(|source| source_error(index, source))?;
            frames.push(AdhesiveDmtContactFrame {
                position_before_m: contact.initial_position_m,
                position_after_m: contact.initial_position_m,
                mean_force_n: [0.0; 3],
                path,
                regime: state.regime(),
                friction_dissipation_j: 0.0,
                prescribed_work_j: 0.0,
                constitutive_residual_n: [0.0; 3],
                compatibility_residual_m: [0.0; 3],
                compatibility_tolerance_m: 0.0,
                work_residual_j: 0.0,
                work_tolerance_j: 0.0,
            });
            states.push(state);
            points.push(Point {
                load_line,
                columns: contact.columns,
                state,
            });
        }
        coupled_stability(
            points
                .iter()
                .map(|p| (p.load_line.normal(), p.load_line.shear())),
            &compliance,
            gate,
        )?;
        let frame = MultipleAdhesiveDmtModalFrame {
            sample: 0,
            observer_pressure_pa: 0.0,
            network_energy_j: 0.0,
            surface_energy_j: 0.0,
            external_work_j: 0.0,
            prescribed_work_j: 0.0,
            contact_work_on_network_j: 0.0,
            network_dissipation_j: 0.0,
            friction_dissipation_j: 0.0,
            energy_residual_j: 0.0,
            energy_tolerance_j: 0.0,
            contact_work_residual_j: 0.0,
            contact_work_tolerance_j: 0.0,
            sweeps: 0,
            iterations: 0,
            evaluations: 0,
            contacts: frames,
        };
        let owner = Self {
            network,
            points,
            config,
            compliance,
            previous: states.clone(),
            candidate: states,
            before_energy: vec![0.0; p],
            free: vec![[0.0; 3]; p],
            reactions: vec![[0.0; 3]; p],
            position: vec![[0.0; 3]; p],
            trials: vec![None; p],
            forces: vec![0.0; n],
            delta: vec![0.0; n],
            staged: frame.clone(),
            frame,
        };
        limit(
            "simultaneous native/DMT initial energy magnitude",
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
    pub fn contact_count(&self) -> usize {
        self.points.len()
    }
    pub fn samples_rendered(&self) -> u64 {
        self.network.samples_rendered()
    }
    pub fn sample_period_s(&self) -> f64 {
        self.network.sample_period_s()
    }
    pub fn state(&self, contact: usize) -> Option<DmtCoupledState> {
        self.points.get(contact).map(|p| p.state)
    }
    pub fn normal(&self, contact: usize) -> Option<&'a SourceDmtSphereContact> {
        self.points.get(contact).map(|p| p.load_line.normal())
    }
    pub fn shear(&self, contact: usize) -> Option<&'a SourceAdhesiveShear> {
        self.points.get(contact).map(|p| p.load_line.shear())
    }
    pub fn work_columns(&self, contact: usize) -> Option<&[Vec<f64>; 3]> {
        self.points.get(contact).map(|p| &p.columns)
    }
    /// Complete row-major response, contact-major then normal/shear 1/shear 2.
    pub fn compliance_m_per_n(&self) -> &[f64] {
        &self.compliance
    }
    pub const fn config(&self) -> MultipleAdhesiveDmtModalConfig {
        self.config
    }
    pub fn last_frame(&self) -> Option<&MultipleAdhesiveDmtModalFrame> {
        (self.frame.sample != 0).then_some(&self.frame)
    }
    pub fn total_energy_j(&self) -> Result<f64, MultipleAdhesiveDmtModalError> {
        let mut energy = Dd::from_f64(self.network.total_energy_j()?);
        for (i, point) in self.points.iter().enumerate() {
            let normal = point
                .load_line
                .normal()
                .evaluate(point.state.position_m()[0])
                .map_err(|e| source_error(i, e.into()))?;
            energy = energy + Dd::from_f64(normal.reversible_energy_j);
        }
        Ok(finite(energy.to_f64())?)
    }

    /// One rate triple per contact in construction order. Zero rates represent
    /// internal two-body contacts; prescribed motion supplies its own signed work.
    pub fn step(
        &mut self,
        external: &[f64],
        prescribed_rate_m_s: &[[f64; 3]],
    ) -> Result<&MultipleAdhesiveDmtModalFrame, MultipleAdhesiveDmtModalError> {
        self.step_under_gate(external, prescribed_rate_m_s, &CancelGate::new_clock_free())
    }
    pub fn step_under_gate(
        &mut self,
        external: &[f64],
        prescribed_rate_m_s: &[[f64; 3]],
        gate: &CancelGate,
    ) -> Result<&MultipleAdhesiveDmtModalFrame, MultipleAdhesiveDmtModalError> {
        self.stage(external, prescribed_rate_m_s, gate)?;
        self.publish();
        Ok(&self.frame)
    }
    /// A downstream owner stages its candidate before any accepted state changes.
    /// Install the returned candidate infallibly after success. Receiver refusal
    /// or late cancellation preserves all native/contact histories and frames.
    pub fn step_with_admission<R, E>(
        &mut self,
        external: &[f64],
        prescribed_rate_m_s: &[[f64; 3]],
        gate: &CancelGate,
        admit: impl FnOnce(&MultipleAdhesiveDmtModalFrame, &[ModalAcousticTimeModel]) -> Result<R, E>,
    ) -> Result<R, E>
    where
        E: From<MultipleAdhesiveDmtModalError>,
    {
        self.stage(external, prescribed_rate_m_s, gate)
            .map_err(E::from)?;
        let received = admit(&self.staged, self.network.staged_components())?;
        poll(Some(gate))
            .map_err(MultipleAdhesiveDmtModalError::from)
            .map_err(E::from)?;
        self.publish();
        Ok(received)
    }
    fn publish(&mut self) {
        self.network.publish_staged();
        for (point, state) in self.points.iter_mut().zip(&self.candidate) {
            point.state = *state;
        }
        core::mem::swap(&mut self.frame, &mut self.staged);
    }

    fn stage(
        &mut self,
        external: &[f64],
        rates: &[[f64; 3]],
        gate: &CancelGate,
    ) -> Result<(), MultipleAdhesiveDmtModalError> {
        poll(Some(gate))?;
        let p = self.points.len();
        let dt = self.network.dt;
        if rates.len() != p || rates.iter().flatten().any(|v| !v.is_finite()) {
            return Err(invalid(
                "one finite prescribed rate triple per original DMT contact required",
            )
            .into());
        }
        let native_before = self.network.total_energy_j()?;
        self.network.prepare_prediction(external, Some(gate))?;
        let mut evaluations = 0;
        for (i, point) in self.points.iter().enumerate() {
            poll(Some(gate))?;
            self.reserve(evaluations, 1, p)?;
            let old = point.state.position_m();
            self.previous[i] = point
                .load_line
                .with_endpoint_velocity(
                    point.state,
                    tangential_velocity(self.network.components(), &point.columns, rates[i])?,
                )
                .map_err(|source| source_error(i, source))?;
            let path = point
                .load_line
                .path(old[0], old[0])
                .map_err(|source| source_error(i, source))?;
            evaluations += 1;
            self.before_energy[i] = path.contact.reversible_energy_j;
            self.reactions[i] = [path.mean_normal_force_n, 0.0, 0.0];
            for axis in 0..3 {
                self.free[i][axis] = finite(
                    old[axis]
                        + finite(rates[i][axis] * dt)?
                        + dot(&point.columns[axis], &self.network.free_delta)?,
                )?;
            }
        }
        let (sweeps, iterations) = self.solve(&mut evaluations, gate)?;
        for mode in 0..self.forces.len() {
            if mode % 256 == 0 {
                poll(Some(gate))?;
            }
            let mut force = Dd::from_f64(external[mode]);
            for (point, reaction) in self.points.iter().zip(&self.reactions) {
                for axis in 0..3 {
                    force = force - product(point.columns[axis][mode], reaction[axis]);
                }
            }
            self.forces[mode] = finite(force.to_f64())?;
        }
        self.network.stage_inner(&self.forces, Some(gate))?;
        for ((delta, end), start) in self
            .delta
            .iter_mut()
            .zip(
                self.network
                    .staged_components()
                    .iter()
                    .flat_map(|m| m.states()),
            )
            .zip(self.network.components().iter().flat_map(|m| m.states()))
        {
            *delta = finite(end.displacement_m_sqrt_kg - start.displacement_m_sqrt_kg)?;
        }
        let mut stored = Dd::ZERO;
        let mut storage_change = Dd::ZERO;
        let mut path_work = Dd::ZERO;
        let mut friction = Dd::ZERO;
        let mut prescribed_work = Dd::ZERO;
        let mut work_tolerance = Dd::ZERO;
        for (i, point) in self.points.iter().enumerate() {
            poll(Some(gate))?;
            self.reserve(evaluations, 1, p - i - 1)?;
            let old = self.previous[i].position_m();
            let mut actual_position = [0.0; 3];
            let mut compatibility = [0.0; 3];
            let solved = self.trials[i]
                .as_ref()
                .expect("a converged sweep visited each contact");
            for axis in 0..3 {
                actual_position[axis] = finite(
                    old[axis]
                        + finite(rates[i][axis] * dt)?
                        + dot(&point.columns[axis], &self.delta)?,
                )?;
                compatibility[axis] = finite(actual_position[axis] - self.position[i][axis])?;
                if compatibility[axis].abs() > solved.compatibility_tolerance_m {
                    return Err(MultipleAdhesiveDmtModalError::Coordinate {
                        contact: i,
                        axis,
                        residual_m: compatibility[axis],
                        tolerance_m: solved.compatibility_tolerance_m,
                    });
                }
            }
            let actual = point
                .load_line
                .advance_actual(
                    self.previous[i],
                    actual_position,
                    self.reactions[i],
                    solved.regime,
                    dt,
                    gate,
                )
                .map_err(|source| source_error(i, source))?;
            evaluations += 1;
            self.candidate[i] = point
                .load_line
                .with_endpoint_velocity(
                    actual.state(),
                    tangential_velocity(
                        self.network.staged_components(),
                        &point.columns,
                        rates[i],
                    )?,
                )
                .map_err(|source| source_error(i, source))?;
            let after = actual.path.contact.reversible_energy_j;
            let change = Dd::from_f64(after) - Dd::from_f64(self.before_energy[i]);
            let mut applied = Dd::ZERO;
            let mut prescribed = Dd::ZERO;
            let mut constitutive = [0.0; 3];
            for axis in 0..3 {
                applied = applied
                    + product(
                        self.reactions[i][axis],
                        finite(actual_position[axis] - old[axis])?,
                    );
                prescribed =
                    prescribed + product(self.reactions[i][axis], finite(rates[i][axis] * dt)?);
                constitutive[axis] = finite(self.reactions[i][axis] - actual.mean_force_n[axis])?;
            }
            let normal_work = Dd::from_f64(actual.path.normal_work_j);
            let loss = Dd::from_f64(actual.friction_dissipation_j);
            let scale = self.before_energy[i]
                .abs()
                .max(after.abs())
                .max(actual.path.normal_work_j.abs())
                .max(actual.friction_dissipation_j.abs())
                .max(finite(applied.to_f64())?.abs());
            let tolerance = finite(
                self.config.contact_work_absolute_tolerance_j
                    + self.config.contact_work_relative_tolerance * scale,
            )?;
            check_work(Some(i), finite((change - normal_work).to_f64())?, tolerance)?;
            let residual = finite((applied - normal_work - loss).to_f64())?;
            check_work(Some(i), residual, tolerance)?;
            stored = stored + Dd::from_f64(after);
            storage_change = storage_change + change;
            path_work = path_work + normal_work;
            friction = friction + loss;
            prescribed_work = prescribed_work + prescribed;
            work_tolerance = work_tolerance + Dd::from_f64(tolerance);
            self.staged.contacts[i] = AdhesiveDmtContactFrame {
                position_before_m: old,
                position_after_m: actual_position,
                mean_force_n: self.reactions[i],
                path: actual.path,
                regime: actual.regime,
                friction_dissipation_j: actual.friction_dissipation_j,
                prescribed_work_j: finite(prescribed.to_f64())?,
                constitutive_residual_n: constitutive,
                compatibility_residual_m: compatibility,
                compatibility_tolerance_m: solved.compatibility_tolerance_m,
                work_residual_j: residual,
                work_tolerance_j: tolerance,
            };
        }
        let mut contact_work = Dd::ZERO;
        let mut external_work = Dd::ZERO;
        for ((held, input), delta) in self.forces.iter().zip(external).zip(&self.delta) {
            contact_work =
                contact_work + (Dd::from_f64(*held) - Dd::from_f64(*input)) * Dd::from_f64(*delta);
            external_work = external_work + product(*input, *delta);
        }
        let tolerance = finite(work_tolerance.to_f64())?;
        let contact_residual =
            finite((contact_work + path_work + friction - prescribed_work).to_f64())?;
        check_work(None, contact_residual, tolerance)?;
        check_work(
            None,
            finite((contact_work + storage_change + friction - prescribed_work).to_f64())?,
            tolerance,
        )?;
        let native = &self.network.staged_frame;
        let native_after = finite(native.total_energy_j())?;
        let native_loss = finite(native.total_dissipation_j())?;
        let after = finite((Dd::from_f64(native_after) + stored).to_f64())?;
        limit(
            "simultaneous native/DMT energy magnitude",
            after.abs(),
            self.network.config.maximum_total_energy_j,
        )?;
        let residual = finite(
            (Dd::from_f64(native_after) - Dd::from_f64(native_before)
                + storage_change
                + Dd::from_f64(native_loss)
                + friction
                - external_work
                - prescribed_work)
                .to_f64(),
        )?;
        let scale = native_before
            .abs()
            .max(native_after.abs())
            .max(finite(stored.to_f64())?.abs())
            .max(finite(external_work.to_f64())?.abs())
            .max(finite(prescribed_work.to_f64())?.abs())
            .max(native_loss.abs())
            .max(finite(friction.to_f64())?.abs());
        let energy_tolerance = finite(
            self.network.config.energy_absolute_tolerance_j
                + self.network.config.energy_relative_tolerance * scale,
        )?;
        check_native(residual, energy_tolerance)?;
        check_native(
            finite((Dd::from_f64(native.external_work_j) - external_work - contact_work).to_f64())?,
            energy_tolerance,
        )?;
        self.staged.sample = native.sample;
        self.staged.observer_pressure_pa = native.observer_pressure_pa;
        self.staged.network_energy_j = native_after;
        self.staged.surface_energy_j = finite(stored.to_f64())?;
        self.staged.external_work_j = finite(external_work.to_f64())?;
        self.staged.prescribed_work_j = finite(prescribed_work.to_f64())?;
        self.staged.contact_work_on_network_j = finite(contact_work.to_f64())?;
        self.staged.network_dissipation_j = native_loss;
        self.staged.friction_dissipation_j = finite(friction.to_f64())?;
        self.staged.energy_residual_j = residual;
        self.staged.energy_tolerance_j = energy_tolerance;
        self.staged.contact_work_residual_j = contact_residual;
        self.staged.contact_work_tolerance_j = tolerance;
        self.staged.sweeps = sweeps;
        self.staged.iterations = iterations;
        self.staged.evaluations = evaluations;
        self.network.validate_staged_material_points(Some(gate))?;
        poll(Some(gate))?;
        Ok(())
    }

    fn reserve(
        &self,
        used: usize,
        requested: usize,
        final_replay: usize,
    ) -> Result<(), MultipleAdhesiveDmtModalError> {
        if used
            .checked_add(requested)
            .and_then(|n| n.checked_add(final_replay))
            .is_none_or(|n| n > self.config.maximum_source_evaluations)
        {
            return Err(MultipleAdhesiveDmtModalError::EvaluationBudget);
        }
        Ok(())
    }
    fn solve(
        &mut self,
        evaluations: &mut usize,
        gate: &CancelGate,
    ) -> Result<(usize, usize), MultipleAdhesiveDmtModalError> {
        let p = self.points.len();
        let rows = 3 * p;
        let mut iterations = 0usize;
        for sweep in 1..=self.config.maximum_sweeps {
            for i in 0..p {
                poll(Some(gate))?;
                self.reserve(
                    *evaluations,
                    self.points[i].load_line.config().maximum_evaluations,
                    p,
                )?;
                let mut free = [0.0; 3];
                for axis in 0..3 {
                    let mut coordinate = Dd::from_f64(self.free[i][axis]);
                    for j in 0..p {
                        if j != i {
                            for other in 0..3 {
                                coordinate = coordinate
                                    - product(
                                        self.compliance[(3 * i + axis) * rows + 3 * j + other],
                                        self.reactions[j][other],
                                    );
                            }
                        }
                    }
                    free[axis] = finite(coordinate.to_f64())?;
                }
                let solved = self.points[i]
                    .load_line
                    .trial(self.previous[i], free, self.network.dt, gate)
                    .map_err(|source| source_error(i, source))?;
                *evaluations += solved.evaluations;
                iterations = iterations
                    .checked_add(solved.iterations)
                    .ok_or(MultipleAdhesiveDmtModalError::EvaluationBudget)?;
                self.reactions[i] = solved.mean_force_n;
                self.trials[i] = Some(solved);
            }
            let mut converged = true;
            for i in 0..p {
                poll(Some(gate))?;
                let solved = self.trials[i]
                    .as_ref()
                    .expect("the sweep visited every contact");
                for axis in 0..3 {
                    let mut coordinate = Dd::from_f64(self.free[i][axis]);
                    let mut scale = self.free[i][axis]
                        .abs()
                        .max(self.previous[i].position_m()[axis].abs());
                    for j in 0..p {
                        for other in 0..3 {
                            let response = product(
                                self.compliance[(3 * i + axis) * rows + 3 * j + other],
                                self.reactions[j][other],
                            );
                            scale = scale.max(finite(response.to_f64())?.abs());
                            coordinate = coordinate - response;
                        }
                    }
                    self.position[i][axis] = finite(coordinate.to_f64())?;
                    // The original local root owns its own residual tolerance.
                    // Only subsequent OTHER-contact updates must settle to
                    // roundoff before one joint native propagation is admitted.
                    let drift = finite(
                        self.position[i][axis] - solved.state().position_m()[axis]
                            + solved.compatibility_residual_m[axis],
                    )?;
                    let tolerance = solved
                        .compatibility_tolerance_m
                        .min(64.0 * f64::EPSILON * scale);
                    converged &= drift.abs() <= tolerance;
                }
            }
            if converged {
                return Ok((sweep, iterations));
            }
        }
        Err(MultipleAdhesiveDmtModalError::NotConverged {
            sweeps: self.config.maximum_sweeps,
        })
    }
}

fn source_error(contact: usize, source: DmtCoupledError) -> MultipleAdhesiveDmtModalError {
    MultipleAdhesiveDmtModalError::Contact { contact, source }
}
fn product(a: f64, b: f64) -> Dd {
    Dd::from_f64(a) * Dd::from_f64(b)
}
fn dot(a: &[f64], b: &[f64]) -> Result<f64, ModalCouplingError> {
    let mut result = Dd::ZERO;
    for (&a, &b) in a.iter().zip(b) {
        result = result + product(a, b);
    }
    finite(result.to_f64())
}
fn check_work(
    contact: Option<usize>,
    residual_j: f64,
    tolerance_j: f64,
) -> Result<(), MultipleAdhesiveDmtModalError> {
    if residual_j.abs() > tolerance_j {
        return Err(MultipleAdhesiveDmtModalError::Work {
            contact,
            residual_j,
            tolerance_j,
        });
    }
    Ok(())
}
fn check_native(residual_j: f64, tolerance_j: f64) -> Result<(), MultipleAdhesiveDmtModalError> {
    if residual_j.abs() > tolerance_j {
        return Err(ModalCouplingError::EnergyBalance {
            residual_j,
            tolerance_j,
        }
        .into());
    }
    Ok(())
}

fn coupled_stability<'a>(
    sources: impl ExactSizeIterator<Item = (&'a SourceDmtSphereContact, &'a SourceAdhesiveShear)>,
    compliance: &[f64],
    gate: &CancelGate,
) -> Result<(), MultipleAdhesiveDmtModalError> {
    let rows = 3 * sources.len();
    let mut roots = Vec::with_capacity(sources.len());
    for (i, (normal, shear)) in sources.enumerate() {
        poll(Some(gate))?;
        let tail = -normal
            .evaluate(0.0)
            .map_err(|e| source_error(i, e.into()))?
            .one_sided_tangents_n_per_m[0];
        // Original DMT path: Nbar' >= -tail/2 and Abar' <= pi*R/2.
        // The shear-normal symmetric cross term is at most tau_k*pi*R/4.
        // Therefore sym(J_i) >= -eta_i I, eta=(a+hypot(a,2c))/2.
        // This intentionally conservative bound also covers arbitrary slip
        // directions; retaining only normal curvature would miss this coupling.
        let a = finite(0.5 * tail)?;
        let c = finite(
            0.25 * core::f64::consts::PI
                * normal.effective_radius_m()
                * shear.law().kinetic_shear_strength_pa(),
        )?;
        let eta = finite(0.5 * a + 0.5 * a.hypot(2.0 * c))?;
        roots.push(det::sqrt(eta));
    }
    let maximum = compliance.iter().fold(0.0_f64, |a, b| a.max(b.abs()));
    for i in 0..rows {
        poll(Some(gate))?;
        let mut row = Dd::ZERO;
        for j in 0..rows {
            if (compliance[i * rows + j] - compliance[j * rows + i]).abs()
                > 64.0 * f64::EPSILON * maximum
            {
                return Err(MultipleAdhesiveDmtModalError::CoupledStability);
            }
            row = row + product(roots[i / 3] * compliance[i * rows + j].abs(), roots[j / 3]);
        }
        // Native compliance is positive semidefinite. This infinity-norm bound
        // on sqrt(eta) S sqrt(eta) screens collective negative curvature without
        // inverting a possibly redundant contact map or altering any source.
        if finite(row.to_f64())? >= 1.0 - 64.0 * f64::EPSILON {
            return Err(MultipleAdhesiveDmtModalError::CoupledStability);
        }
    }
    Ok(())
}
