//! Simultaneous sourced DMT attraction and real-junction shear on native bodies.
//!
//! Three complete work columns map every retained modal coordinate to signed
//! approach and two local shear coordinates. Their transposes apply the same
//! resisting efforts. The original network supplies the entire reciprocal
//! response, including normal/shear coupling, and advances its actual states.
//! Normal surface work and friction loss are then reevaluated on that motion.
//! Only a complete mechanical + signed surface + friction work account may
//! publish the bodies, occupied friction branch and last frame together.
//! Actual endpoint velocity retains kinetic occupancy even when a step has
//! zero net slip; that step does not certify an event-resolved physical stop.
//!
//! The source-bound DMT and constant-junction model retain their scope: fixed
//! source conditions held over each step, fixed local work columns and straight
//! coordinate paths. Explicit prescribed-temperature transfer retains native
//! history and separately charges the contact potential's parameter work.
//! Gap attraction has no solid shear area. There is no inferred partial-slip
//! stiffness, JKR mixed-mode law, evolving contact frame or heat destination.

use super::{
    CoupledModalSystem, ModalCouplingError, dot, finite, invalid, limit, network_response, poll,
};
use crate::modal_acoustic_time::ModalAcousticTimeModel;
use fs_contact::interface_binding::adhesive_shear::SourceAdhesiveShear;
use fs_contact::normal_patch::adhesive::temperature::SourceDmtThermalRefresh;
use fs_contact::normal_patch::dmt::coupled::{
    DmtAdhesivePath, DmtCoupledConfig, DmtCoupledError, DmtCoupledLoadLine, DmtCoupledRegime,
    DmtCoupledState,
};
use fs_contact::normal_patch::dmt::{DmtError, SourceDmtSphereContact};
use fs_exec::CancelGate;
use fs_matdb::MaterialStateId;
use fs_math::dd::Dd;

/// Failed steps publish neither native motion nor contact history.
#[derive(Debug)]
pub enum AdhesiveDmtModalError {
    Network(ModalCouplingError),
    Contact(DmtCoupledError),
    CoordinateMismatch {
        component: usize,
        residual_m: f64,
        tolerance_m: f64,
    },
    /// Contact work is checked independently of unrelated native stored energy.
    ContactWork {
        residual_j: f64,
        tolerance_j: f64,
    },
}
impl core::fmt::Display for AdhesiveDmtModalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native adhesive DMT trial refused: {self:?}")
    }
}
impl std::error::Error for AdhesiveDmtModalError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Network(error) => Some(error),
            Self::Contact(error) => Some(error),
            Self::CoordinateMismatch { .. } | Self::ContactWork { .. } => None,
        }
    }
}
impl From<ModalCouplingError> for AdhesiveDmtModalError {
    fn from(error: ModalCouplingError) -> Self {
        Self::Network(error)
    }
}
impl From<DmtCoupledError> for AdhesiveDmtModalError {
    fn from(error: DmtCoupledError) -> Self {
        Self::Contact(error)
    }
}
impl From<DmtError> for AdhesiveDmtModalError {
    fn from(error: DmtError) -> Self {
        Self::Contact(error.into())
    }
}

/// Accepted actual-body path and separate mechanical/surface/friction account.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdhesiveDmtModalFrame {
    pub sample: u64,
    pub observer_pressure_pa: f64,
    pub network_energy_j: f64,
    /// Original signed DMT potential; attraction is not dissipative heat.
    pub surface_energy_j: f64,
    pub external_work_j: f64,
    /// Applied resisting effort dotted with prescribed coordinate increment.
    pub prescribed_work_j: f64,
    pub network_dissipation_j: f64,
    /// Sourced kinetic junction strength times integrated real area and slip.
    /// Its eventual thermal recipient remains an explicit downstream decision.
    pub friction_dissipation_j: f64,
    /// Work of the contact correction actually represented in native forcing.
    pub contact_work_on_network_j: f64,
    /// Original normal-potential change plus junction loss and native contact
    /// work, minus work supplied through prescribed contact motion.
    pub contact_work_residual_j: f64,
    /// Original force-accuracy budget times actual travel, plus contact-sized
    /// arithmetic allowance. Unrelated native energy never enlarges this gate.
    pub contact_work_tolerance_j: f64,
    /// Complete native plus contact storage and dissipation balance.
    pub energy_residual_j: f64,
    pub energy_tolerance_j: f64,
    /// Held efforts actually applied to the complete original modal network.
    pub mean_force_n: [f64; 3],
    /// Independent original constitutive evaluation on actual staged motion.
    pub path: DmtAdhesivePath,
    /// Path regime; a detaching successor may already be Released.
    pub regime: DmtCoupledRegime,
    pub position_m: [f64; 3],
    pub compatibility_residual_m: [f64; 3],
    pub compatibility_tolerance_m: f64,
    /// Applied minus actual constitutive effort, separately from compatibility.
    pub constitutive_residual_n: [f64; 3],
    pub iterations: usize,
    /// Includes the reserved actual-motion constitutive evaluation.
    pub evaluations: usize,
}

/// Instantaneous prescribed contact-source change at unchanged native state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdhesiveDmtSourceRefreshReport {
    pub sample: u64,
    pub energy_before_j: f64,
    pub energy_after_j: f64,
    pub surface_energy_before_j: f64,
    pub surface_energy_after_j: f64,
    /// Signed external coefficient work, distinct from slip work and heat.
    pub parameter_work_j: f64,
    /// Contact-potential change minus parameter work, computed without
    /// subtracting the possibly much larger native energy observations.
    pub energy_residual_j: f64,
    pub energy_tolerance_j: f64,
}

/// One source-pinned junction and unchanged native modal component network.
pub struct AdhesiveDmtModalSystem<'a> {
    network: CoupledModalSystem,
    load_line: DmtCoupledLoadLine<'a>,
    state: DmtCoupledState,
    columns: [Vec<f64>; 3],
    forces: Vec<f64>,
    config: DmtCoupledConfig,
    last: Option<AdhesiveDmtModalFrame>,
}

impl<'a> AdhesiveDmtModalSystem<'a> {
    /// Attach the explicitly supplied physical approach/shear coordinates to
    /// the current sample-zero geometry, including nonzero initial modal q.
    /// Columns are complete component-then-mode vectors: `B_i^T delta_q`
    /// changes coordinate i and `-B_i*F_i` is its reciprocal modal reaction.
    /// Ordered receiving material identities must match the original source.
    /// Actual initial tangential body velocity selects the occupied branch;
    /// moving contact therefore cannot start with an artificial static capture.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        network: CoupledModalSystem,
        normal: &'a SourceDmtSphereContact,
        shear: &'a SourceAdhesiveShear,
        surface_materials: [&MaterialStateId; 2],
        initial_position_m: [f64; 3],
        columns: [Vec<f64>; 3],
        config: DmtCoupledConfig,
        gate: &CancelGate,
    ) -> Result<Self, AdhesiveDmtModalError> {
        poll(Some(gate))?;
        if network.samples_rendered() != 0 {
            return Err(invalid("DMT junction attachment requires a sample-zero network").into());
        }
        let materials = normal
            .normal_state()
            .elastic()
            .interface()
            .surface_materials();
        if surface_materials[0] != &materials[0] || surface_materials[1] != &materials[1] {
            return Err(invalid("DMT receiving surfaces differ from the ordered source").into());
        }
        if config.maximum_evaluations < 2
            || columns.iter().any(|column| {
                column.len() != network.mode_count()
                    || column.iter().any(|value| !value.is_finite())
            })
        {
            return Err(invalid(
                "DMT work columns must retain all modes and reserve an actual-state evaluation",
            )
            .into());
        }
        let mut compliance = [[0.0; 3]; 3];
        for j in 0..3 {
            let response = network_response(&network, &columns[j], gate)?;
            for i in 0..3 {
                compliance[i][j] = dot(&columns[i], &response)?;
            }
        }
        let mut solve_config = config;
        solve_config.maximum_evaluations -= 1;
        let load_line = DmtCoupledLoadLine::new(normal, shear, compliance, solve_config)?;
        let initial_velocity = tangential_velocity(network.components(), &columns, [0.0; 3])?;
        let state = load_line.initialize_with_velocity(initial_position_m, initial_velocity)?;
        let system = Self {
            forces: vec![0.0; network.mode_count()],
            network,
            load_line,
            state,
            columns,
            config,
            last: None,
        };
        limit(
            "native/DMT initial energy magnitude",
            system.total_energy_j()?.abs(),
            system.network.config.maximum_total_energy_j,
        )?;
        poll(Some(gate))?;
        Ok(system)
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
    pub const fn state(&self) -> DmtCoupledState {
        self.state
    }
    pub const fn normal(&self) -> &'a SourceDmtSphereContact {
        self.load_line.normal()
    }
    pub const fn shear(&self) -> &'a SourceAdhesiveShear {
        self.load_line.shear()
    }
    pub const fn work_columns(&self) -> &[Vec<f64>; 3] {
        &self.columns
    }
    pub const fn compliance_m_per_n(&self) -> [[f64; 3]; 3] {
        self.load_line.compliance_m_per_n()
    }
    pub const fn config(&self) -> DmtCoupledConfig {
        self.config
    }
    pub const fn last_frame(&self) -> Option<&AdhesiveDmtModalFrame> {
        self.last.as_ref()
    }
    pub fn total_energy_j(&self) -> Result<f64, AdhesiveDmtModalError> {
        Ok(finite(
            self.network.total_energy_j()?
                + self
                    .normal()
                    .evaluate(self.state.position_m()[0])?
                    .reversible_energy_j,
        )?)
    }

    /// Apply a prescribed temperature successor to the same actual native
    /// contact coordinates, without advancing a sample or resetting a branch.
    /// The original normal source owns independent bulk-contact and interface
    /// temperatures; the original selected shear claims follow the interface.
    ///
    /// The complete native motions, work columns, cross-axis compliance and
    /// last historical frame remain unchanged. This operation changes contact
    /// material parameters only; it does not change the modal network's bulk
    /// stiffness, supply thermal transport/expansion, or infer surface entropy.
    /// Surface-potential change is independently charged as signed parameter
    /// work. Changing shear strength creates no heat at fixed actual position.
    /// Any source, stability, energy or cancellation refusal retains this owner.
    pub fn apply_source_refresh(
        &mut self,
        refreshed: &'a SourceDmtThermalRefresh,
        next_shear: &'a SourceAdhesiveShear,
        gate: &CancelGate,
    ) -> Result<AdhesiveDmtSourceRefreshReport, AdhesiveDmtModalError> {
        self.apply_source_refresh_with_admission(refreshed, next_shear, gate, |report, _| {
            Ok(*report)
        })
    }

    /// Stage the receiving source/history owner before publishing either local
    /// source replacement. The callback returns an infallibly installable
    /// candidate; refusal or late cancellation retains every accepted owner.
    pub fn apply_source_refresh_with_admission<R, E>(
        &mut self,
        refreshed: &'a SourceDmtThermalRefresh,
        next_shear: &'a SourceAdhesiveShear,
        gate: &CancelGate,
        admit: impl FnOnce(&AdhesiveDmtSourceRefreshReport, DmtCoupledState) -> Result<R, E>,
    ) -> Result<R, E>
    where
        E: From<AdhesiveDmtModalError>,
    {
        let (next, state, report) = self
            .prepare_source_refresh(refreshed, next_shear, gate)
            .map_err(E::from)?;
        let receiving = admit(&report, state)?;
        poll(Some(gate))
            .map_err(AdhesiveDmtModalError::from)
            .map_err(E::from)?;
        self.load_line = next;
        self.state = state;
        Ok(receiving)
    }

    fn prepare_source_refresh(
        &self,
        refreshed: &'a SourceDmtThermalRefresh,
        next_shear: &'a SourceAdhesiveShear,
        gate: &CancelGate,
    ) -> Result<
        (
            DmtCoupledLoadLine<'a>,
            DmtCoupledState,
            AdhesiveDmtSourceRefreshReport,
        ),
        AdhesiveDmtModalError,
    > {
        poll(Some(gate))?;
        let before = self.total_energy_j()?;
        let surface_before = self
            .normal()
            .evaluate(self.state.position_m()[0])?
            .reversible_energy_j;
        let (next, state) = self
            .load_line
            .with_temperature_refresh(self.state, refreshed, next_shear, gate)?;
        let surface_after = next
            .normal()
            .evaluate(state.position_m()[0])?
            .reversible_energy_j;
        let after = finite(self.network.total_energy_j()? + surface_after)?;
        let work = finite(refreshed.parameter_work_j())?;
        limit(
            "native/DMT refreshed energy magnitude",
            after.abs(),
            self.network.config.maximum_total_energy_j,
        )?;
        let change = Dd::from_f64(surface_after) - Dd::from_f64(surface_before);
        let residual = finite((change - Dd::from_f64(work)).to_f64())?;
        // The sealed source has no mechanical solve error at held geometry.
        // Only original potential/parameter arithmetic enters this allowance.
        let tolerance = finite(
            128.0
                * f64::EPSILON
                * surface_before
                    .abs()
                    .max(surface_after.abs())
                    .max(work.abs()),
        )?;
        check_contact_work(residual, tolerance)?;
        check_contact_work(
            finite((change - Dd::from_f64(refreshed.potential_change_j())).to_f64())?,
            tolerance,
        )?;
        let report = AdhesiveDmtSourceRefreshReport {
            sample: self.samples_rendered(),
            energy_before_j: before,
            energy_after_j: after,
            surface_energy_before_j: surface_before,
            surface_energy_after_j: surface_after,
            parameter_work_j: work,
            energy_residual_j: residual,
            energy_tolerance_j: tolerance,
        };
        poll(Some(gate))?;
        Ok((next, state, report))
    }

    /// Advance held modal forces and a prescribed local coordinate rate [m/s].
    /// Pass zero rates for a wholly internal two-body contact. All physical
    /// timestep, source, force and native work gates remain active.
    pub fn step(
        &mut self,
        external: &[f64],
        prescribed_rate_m_s: [f64; 3],
    ) -> Result<&AdhesiveDmtModalFrame, AdhesiveDmtModalError> {
        self.step_under_gate(external, prescribed_rate_m_s, &CancelGate::new_clock_free())
    }

    pub fn step_under_gate(
        &mut self,
        external: &[f64],
        prescribed_rate_m_s: [f64; 3],
        gate: &CancelGate,
    ) -> Result<&AdhesiveDmtModalFrame, AdhesiveDmtModalError> {
        let (state, accepted) = self.prepare_step(external, prescribed_rate_m_s, gate)?;
        self.network.publish_staged();
        self.state = state;
        Ok(self.last.insert(accepted))
    }

    /// Stage another physical owner's candidate against the original contact
    /// account and actual native endpoint. Install that candidate only after
    /// success. No native state, source branch, clock or frame publishes on
    /// receiver refusal or cancellation, including a late callback request.
    pub fn step_with_admission<R, E>(
        &mut self,
        external: &[f64],
        prescribed_rate_m_s: [f64; 3],
        gate: &CancelGate,
        admit: impl FnOnce(&AdhesiveDmtModalFrame, &[ModalAcousticTimeModel]) -> Result<R, E>,
    ) -> Result<R, E>
    where
        E: From<AdhesiveDmtModalError>,
    {
        let (state, accepted) = self
            .prepare_step(external, prescribed_rate_m_s, gate)
            .map_err(E::from)?;
        let receiving = admit(&accepted, self.network.staged_components())?;
        poll(Some(gate))
            .map_err(AdhesiveDmtModalError::from)
            .map_err(E::from)?;
        self.network.publish_staged();
        self.state = state;
        self.last = Some(accepted);
        Ok(receiving)
    }

    fn prepare_step(
        &mut self,
        external: &[f64],
        prescribed_rate_m_s: [f64; 3],
        gate: &CancelGate,
    ) -> Result<(DmtCoupledState, AdhesiveDmtModalFrame), AdhesiveDmtModalError> {
        poll(Some(gate))?;
        let before = self.total_energy_j()?;
        let old = self.state.position_m();
        let surface_before = self.normal().evaluate(old[0])?.reversible_energy_j;
        let dt = self.network.dt;
        let mut prescribed = [0.0; 3];
        for i in 0..3 {
            prescribed[i] = finite(prescribed_rate_m_s[i] * dt)?;
        }
        // A prescribed support rate participates in actual relative velocity.
        // Its change cannot turn a moving junction into a fresh static contact.
        let previous = self.load_line.with_endpoint_velocity(
            self.state,
            tangential_velocity(
                self.network.components(),
                &self.columns,
                prescribed_rate_m_s,
            )?,
        )?;
        self.network.prepare_prediction(external, Some(gate))?;
        let mut free = [0.0; 3];
        for i in 0..3 {
            free[i] =
                finite(old[i] + prescribed[i] + dot(&self.columns[i], &self.network.free_delta)?)?;
        }
        let solved = self.load_line.trial(previous, free, dt, gate)?;
        for (j, (force, input)) in self.forces.iter_mut().zip(external).enumerate() {
            if j % 256 == 0 {
                poll(Some(gate))?;
            }
            let mut sum = Dd::from_f64(*input);
            for i in 0..3 {
                sum = sum - Dd::from_f64(self.columns[i][j]) * Dd::from_f64(solved.mean_force_n[i]);
            }
            *force = finite(sum.to_f64())?;
        }
        self.network.stage_inner(&self.forces, Some(gate))?;
        // Reconstruct material motion from accepted and staged native states.
        // This also keeps the original scratch's actual-increment semantics.
        for (j, ((change, end), start)) in self
            .network
            .free_delta
            .iter_mut()
            .zip(self.network.candidates.iter().flat_map(|m| m.states()))
            .zip(self.network.models.iter().flat_map(|m| m.states()))
            .enumerate()
        {
            if j % 256 == 0 {
                poll(Some(gate))?;
            }
            *change = finite(end.displacement_m_sqrt_kg - start.displacement_m_sqrt_kg)?;
        }
        let mut actual_position = [0.0; 3];
        let mut compatibility = [0.0; 3];
        for i in 0..3 {
            actual_position[i] =
                finite(old[i] + prescribed[i] + dot(&self.columns[i], &self.network.free_delta)?)?;
            compatibility[i] = finite(actual_position[i] - solved.state().position_m()[i])?;
            if compatibility[i].abs() > solved.compatibility_tolerance_m {
                return Err(AdhesiveDmtModalError::CoordinateMismatch {
                    component: i,
                    residual_m: compatibility[i],
                    tolerance_m: solved.compatibility_tolerance_m,
                });
            }
        }
        let actual = self.load_line.advance_actual(
            previous,
            actual_position,
            solved.mean_force_n,
            solved.regime,
            dt,
            gate,
        )?;
        let frame = &self.network.staged_frame;
        let network_energy = finite(frame.total_energy_j())?;
        let surface_energy = actual.path.contact.reversible_energy_j;
        let after = finite(network_energy + surface_energy)?;
        limit(
            "native/DMT energy magnitude",
            after.abs(),
            self.network.config.maximum_total_energy_j,
        )?;
        let network_loss = finite(frame.total_dissipation_j())?;
        let friction = actual.friction_dissipation_j;
        let loss = finite(network_loss + friction)?;
        let external_work = dot(external, &self.network.free_delta)?;
        let prescribed_work = dot(&solved.mean_force_n, &prescribed)?;
        let mut represented_work = Dd::ZERO;
        for (j, ((held, input), change)) in self
            .forces
            .iter()
            .zip(external)
            .zip(&self.network.free_delta)
            .enumerate()
        {
            if j % 256 == 0 {
                poll(Some(gate))?;
            }
            represented_work = represented_work
                + (Dd::from_f64(*held) - Dd::from_f64(*input)) * Dd::from_f64(*change);
        }
        let contact_work = finite(represented_work.to_f64())?;
        let mut transpose_work = Dd::ZERO;
        for i in 0..3 {
            transpose_work = transpose_work
                + Dd::from_f64(solved.mean_force_n[i])
                    * Dd::from_f64(dot(&self.columns[i], &self.network.free_delta)?);
        }
        let storage_change = Dd::from_f64(surface_energy) - Dd::from_f64(surface_before);
        let work_scale = surface_before
            .abs()
            .max(surface_energy.abs())
            .max(contact_work.abs())
            .max(prescribed_work.abs())
            .max(actual.path.normal_work_j.abs())
            .max(friction.abs());
        let contact_tolerance = contact_work_tolerance(
            self.config,
            old,
            actual_position,
            core::array::from_fn(|i| {
                solved.mean_force_n[i]
                    .abs()
                    .max(actual.mean_force_n[i].abs())
            }),
            solved.regime,
            work_scale,
        )?;
        check_contact_work(
            finite((represented_work + transpose_work).to_f64())?,
            contact_tolerance,
        )?;
        check_contact_work(
            finite((storage_change - Dd::from_f64(actual.path.normal_work_j)).to_f64())?,
            contact_tolerance,
        )?;
        let contact_residual = finite(
            (storage_change + Dd::from_f64(friction) + represented_work
                - Dd::from_f64(prescribed_work))
            .to_f64(),
        )?;
        check_contact_work(contact_residual, contact_tolerance)?;
        let residual = finite((after - before) + loss - external_work - prescribed_work)?;
        let scale = before
            .abs()
            .max(after.abs())
            .max(external_work.abs())
            .max(prescribed_work.abs())
            .max(loss.abs());
        let tolerance = finite(
            self.network.config.energy_absolute_tolerance_j
                + self.network.config.energy_relative_tolerance * scale,
        )?;
        if residual.abs() > tolerance {
            return Err(ModalCouplingError::EnergyBalance {
                residual_j: residual,
                tolerance_j: tolerance,
            }
            .into());
        }
        let mut constitutive = [0.0; 3];
        for i in 0..3 {
            constitutive[i] = finite(solved.mean_force_n[i] - actual.mean_force_n[i])?;
        }
        let accepted = AdhesiveDmtModalFrame {
            sample: frame.sample,
            observer_pressure_pa: frame.observer_pressure_pa,
            network_energy_j: network_energy,
            surface_energy_j: surface_energy,
            external_work_j: external_work,
            prescribed_work_j: prescribed_work,
            network_dissipation_j: network_loss,
            friction_dissipation_j: friction,
            contact_work_on_network_j: contact_work,
            contact_work_residual_j: contact_residual,
            contact_work_tolerance_j: contact_tolerance,
            energy_residual_j: residual,
            energy_tolerance_j: tolerance,
            mean_force_n: solved.mean_force_n,
            path: actual.path,
            regime: actual.regime,
            position_m: actual_position,
            compatibility_residual_m: compatibility,
            compatibility_tolerance_m: solved.compatibility_tolerance_m,
            constitutive_residual_n: constitutive,
            iterations: solved.iterations,
            evaluations: solved.evaluations + 1,
        };
        let state = self.load_line.with_endpoint_velocity(
            actual.state(),
            tangential_velocity(
                self.network.staged_components(),
                &self.columns,
                prescribed_rate_m_s,
            )?,
        )?;
        self.network.validate_staged_material_points(Some(gate))?;
        poll(Some(gate))?;
        Ok((state, accepted))
    }
}

// A constitutive force error performs at most its admitted magnitude times
// actual coordinate travel. A sticking path additionally admits its original
// small coordinate residual, whose maximum work is shear magnitude times that
// residual. Add only contact-sized roundoff, never a fraction of unrelated
// native storage or background work. All terms have units of J.
fn contact_work_tolerance(
    config: DmtCoupledConfig,
    from: [f64; 3],
    to: [f64; 3],
    force: [f64; 3],
    regime: DmtCoupledRegime,
    scale_j: f64,
) -> Result<f64, ModalCouplingError> {
    let mut tolerance = Dd::from_f64(128.0 * f64::EPSILON * scale_j);
    for i in 0..3 {
        let force_error =
            config.force_absolute_tolerance_n + config.force_relative_tolerance * force[i].abs();
        tolerance =
            tolerance + Dd::from_f64(force_error) * Dd::from_f64(finite(to[i] - from[i])?.abs());
    }
    if regime == DmtCoupledRegime::Sticking {
        let coordinate_error = config.coordinate_absolute_tolerance_m
            + config.coordinate_relative_tolerance * finite(to[0] - from[0])?.abs();
        tolerance =
            tolerance + Dd::from_f64(coordinate_error) * Dd::from_f64(force[1].hypot(force[2]));
    }
    finite(tolerance.to_f64())
}

fn check_contact_work(residual_j: f64, tolerance_j: f64) -> Result<(), AdhesiveDmtModalError> {
    if residual_j.abs() > tolerance_j {
        Err(AdhesiveDmtModalError::ContactWork {
            residual_j,
            tolerance_j,
        })
    } else {
        Ok(())
    }
}

fn tangential_velocity(
    components: &[ModalAcousticTimeModel],
    columns: &[Vec<f64>; 3],
    prescribed_rate_m_s: [f64; 3],
) -> Result<[f64; 2], ModalCouplingError> {
    let mut velocity = [prescribed_rate_m_s[1], prescribed_rate_m_s[2]];
    for (index, modal) in components
        .iter()
        .flat_map(|model| model.states())
        .enumerate()
    {
        for tangent in 0..2 {
            velocity[tangent] = finite(
                velocity[tangent] + columns[tangent + 1][index] * modal.velocity_m_sqrt_kg_per_s,
            )?;
        }
    }
    Ok(velocity)
}
