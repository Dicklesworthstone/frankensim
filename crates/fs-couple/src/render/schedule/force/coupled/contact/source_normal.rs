//! Signed sourced normal work on the complete original modal network.
//!
//! The existing native held-force response supplies the contact compliance;
//! its original stage retains every modal, connection, material-memory and
//! radiation state and loss. Only one scalar contact equation is condensed.
//! The same source-owned path effort and endpoint derivative used by
//! `BoundNormalStorage` act through the unchanged physical work column.
//!
//! Sources and geometry are frozen. Contact work uses the straight approach
//! path between native endpoints, as in the existing DMT network coupling.
//! Attached JKR stays on its explicitly selected branch; its fold/domain and
//! unresolved mechanical branch refuse. This is a finite-step conservative
//! coupling, without within-sample attachment or pull-off event localization.

use super::{CoupledModalSystem, ModalCouplingError, dot, extension, finite, invalid, limit, poll};
use crate::modal_acoustic_time::ModalAcousticTimeModel;
use crate::unilateral_contact::binding::BoundNormalContact;
use fs_contact::normal_patch::{NormalLawSource, NormalLawSourceError};
use fs_exec::CancelGate;
use fs_math::dd::Dd;

mod solve;

/// Explicit physical support and accuracy controls; bounds never clip a state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SourceNormalModalConfig {
    /// Total safeguarded scalar refinements, including JKR stable-root isolation.
    pub maximum_iterations: usize,
    /// Maximum absolute held reaction and endpoint constitutive force [N].
    pub maximum_force_n: f64,
    /// Finite caller-selected approach interval inside the original source [m].
    pub minimum_approach_m: f64,
    pub maximum_approach_m: f64,
    pub force_absolute_tolerance_n: f64,
    pub force_relative_tolerance: f64,
    pub coordinate_absolute_tolerance_m: f64,
    pub coordinate_relative_tolerance: f64,
    /// Independent contact-work gate; unrelated native storage cannot enlarge it.
    pub contact_energy_absolute_tolerance_j: f64,
    pub contact_energy_relative_tolerance: f64,
}

#[derive(Debug)]
pub enum SourceNormalModalError {
    Network(ModalCouplingError),
    Source(NormalLawSourceError),
    Unresolved(&'static str),
    CoordinateMismatch { residual_m: f64, tolerance_m: f64 },
    ContactWork { residual_j: f64, tolerance_j: f64 },
}
impl core::fmt::Display for SourceNormalModalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native source-normal contact refused: {self:?}")
    }
}
impl std::error::Error for SourceNormalModalError {}
impl From<ModalCouplingError> for SourceNormalModalError {
    fn from(error: ModalCouplingError) -> Self {
        Self::Network(error)
    }
}
impl From<NormalLawSourceError> for SourceNormalModalError {
    fn from(error: NormalLawSourceError) -> Self {
        Self::Source(error)
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SourceNormalModalFrame {
    pub sample: u64,
    pub observer_pressure_pa: f64,
    /// All original native stores, including admitted material/radiation memory.
    pub network_energy_j: f64,
    /// Original absolute contact potential, including signed adhesive storage.
    pub contact_energy_j: f64,
    /// Authored native force work, excluding the solved contact reaction.
    pub external_work_j: f64,
    /// Actual represented contact force applied to each native displacement.
    pub contact_work_on_network_j: f64,
    pub network_dissipation_j: f64,
    pub energy_residual_j: f64,
    pub energy_tolerance_j: f64,
    pub contact_work_residual_j: f64,
    pub contact_work_tolerance_j: f64,
    /// Signed held path effort; a negative value attracts the native bodies.
    pub normal_force_n: f64,
    pub endpoint_normal_force_n: f64,
    pub approach_before_m: f64,
    pub approach_after_m: f64,
    pub constitutive_residual_n: f64,
    pub force_tolerance_n: f64,
    pub coordinate_residual_m: f64,
    pub iterations: usize,
}

/// One source contact, with the continuing native network as sole state owner.
pub struct SourceNormalModalSystem {
    network: CoupledModalSystem,
    contact: BoundNormalContact,
    config: SourceNormalModalConfig,
    compliance_m_per_n: f64,
    forces: Vec<f64>,
    last: Option<SourceNormalModalFrame>,
}

impl SourceNormalModalSystem {
    /// Attach an original sourced potential and full signed mechanical row.
    /// Native modes, all bilateral connections, rest loads, initial material
    /// histories and admitted radiation stay in their existing owners.
    /// Initial vibration and contact storage are retained; no equilibrium is
    /// imposed. Current support is conservative Hertz, DMT, attached JKR and
    /// retained finite-gap curves. Hunt--Crossley loss and a Maugis derivative
    /// require their own admitted paths and are not silently omitted.
    pub fn new(
        network: CoupledModalSystem,
        contact: BoundNormalContact,
        config: SourceNormalModalConfig,
        gate: &CancelGate,
    ) -> Result<Self, SourceNormalModalError> {
        poll(Some(gate))?;
        network.require_time_domain()?;
        if network.samples_rendered() != 0 || contact.work_column().len() != network.mode_count() {
            return Err(invalid(
                "source contact needs a sample-zero network and complete native row",
            )
            .into());
        }
        if !(1..=128).contains(&config.maximum_iterations)
            || [
                config.maximum_force_n,
                config.force_absolute_tolerance_n,
                config.coordinate_absolute_tolerance_m,
                config.contact_energy_absolute_tolerance_j,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
            || [
                config.force_relative_tolerance,
                config.coordinate_relative_tolerance,
                config.contact_energy_relative_tolerance,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v < 0.0 || *v >= 1.0)
            || !config.minimum_approach_m.is_finite()
            || !config.maximum_approach_m.is_finite()
            || config.minimum_approach_m >= config.maximum_approach_m
        {
            return Err(invalid(
                "finite source-normal domain and positive explicit accuracy/work bounds required",
            )
            .into());
        }
        match contact.source() {
            NormalLawSource::MaugisAttached(_) => {
                return Err(invalid(
                    "native source-normal solve requires the original analytic path derivative",
                )
                .into());
            }
            NormalLawSource::HertzPoint(source) if source.dissipation_s_per_m() != 0.0 => {
                return Err(invalid(
                    "dissipative point contact requires its original timed loss path",
                )
                .into());
            }
            _ => {}
        }
        // Admit the entire explicit scalar bracket through the original owner.
        observation(contact.source(), config.minimum_approach_m)?;
        observation(contact.source(), config.maximum_approach_m)?;
        let response = network.linear_response(contact.work_column(), gate)?;
        let compliance = dot(contact.work_column(), &response)?;
        if compliance <= 0.0 {
            return Err(
                invalid("source contact requires positive original native compliance").into(),
            );
        }
        let owner = Self {
            forces: vec![0.0; network.mode_count()],
            network,
            contact,
            config,
            compliance_m_per_n: compliance,
            last: None,
        };
        let approach = owner.approach(owner.components())?;
        owner.check_approach(approach)?;
        let (_, force) = observation(owner.contact.source(), approach)?;
        limit(
            "initial sourced normal force",
            force.abs(),
            config.maximum_force_n,
        )?;
        owner.point_kinematics(
            approach,
            normal_velocity(owner.components(), owner.contact.work_column())?,
        )?;
        limit(
            "initial native/source energy magnitude",
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
    pub const fn contact(&self) -> &BoundNormalContact {
        &self.contact
    }
    pub const fn config(&self) -> SourceNormalModalConfig {
        self.config
    }
    pub const fn compliance_m_per_n(&self) -> f64 {
        self.compliance_m_per_n
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
    pub const fn last_frame(&self) -> Option<&SourceNormalModalFrame> {
        self.last.as_ref()
    }

    pub fn total_energy_j(&self) -> Result<f64, SourceNormalModalError> {
        Ok(finite(
            self.network.total_energy_j()?
                + observation(self.contact.source(), self.approach(self.components())?)?.0,
        )?)
    }

    pub fn step(
        &mut self,
        external: &[f64],
    ) -> Result<&SourceNormalModalFrame, SourceNormalModalError> {
        self.step_under_gate(external, &CancelGate::new_clock_free())
    }
    pub fn step_under_gate(
        &mut self,
        external: &[f64],
        gate: &CancelGate,
    ) -> Result<&SourceNormalModalFrame, SourceNormalModalError> {
        let accepted = self.stage(external, gate)?;
        self.network.publish_staged();
        Ok(self.last.insert(accepted))
    }

    /// Admit receiving state before any original native body, memory, clock or
    /// port publishes. The callback returns a candidate; installing it after
    /// success must be infallible. A late cancellation preserves all owners.
    pub fn step_with_admission<R, E>(
        &mut self,
        external: &[f64],
        gate: &CancelGate,
        admit: impl FnOnce(&SourceNormalModalFrame, &[ModalAcousticTimeModel]) -> Result<R, E>,
    ) -> Result<R, E>
    where
        E: From<SourceNormalModalError>,
    {
        let accepted = self.stage(external, gate).map_err(E::from)?;
        let receiving = admit(&accepted, &self.network.candidates)?;
        poll(Some(gate))
            .map_err(SourceNormalModalError::from)
            .map_err(E::from)?;
        self.network.publish_staged();
        self.last = Some(accepted);
        Ok(receiving)
    }

    fn stage(
        &mut self,
        external: &[f64],
        gate: &CancelGate,
    ) -> Result<SourceNormalModalFrame, SourceNormalModalError> {
        poll(Some(gate))?;
        let network_before = self.network.total_energy_j()?;
        let before = self.approach(self.components())?;
        self.check_approach(before)?;
        let (energy_before, _) = observation(self.contact.source(), before)?;
        self.point_kinematics(
            before,
            normal_velocity(self.components(), self.contact.work_column())?,
        )?;
        self.network.prepare_prediction(external, Some(gate))?;
        let free_increment = dot(self.contact.work_column(), &self.network.free_delta)?;
        let solution = solve::solve(
            self.contact.source(),
            before,
            free_increment,
            self.compliance_m_per_n,
            self.config,
            gate,
        )?;
        for ((held, input), weight) in self
            .forces
            .iter_mut()
            .zip(external)
            .zip(self.contact.work_column())
        {
            *held = finite(input - weight * solution.force_n)?;
        }
        self.network.stage_inner(&self.forces, Some(gate))?;
        let after = self.approach(&self.network.candidates)?;
        self.check_approach(after)?;
        let increment = finite(after - before)?;
        let native_increment = dot(self.contact.work_column(), &self.network.free_delta)?;
        let coordinate_residual = check_coordinate(after, solution.approach_m, self.config)?;
        check_coordinate(increment, native_increment, self.config)?;
        let (energy_after, endpoint_force) = observation(self.contact.source(), after)?;
        limit(
            "native sourced normal force",
            solution.force_n.abs(),
            self.config.maximum_force_n,
        )?;
        limit(
            "native sourced endpoint force",
            endpoint_force.abs(),
            self.config.maximum_force_n,
        )?;
        self.point_kinematics(
            after,
            normal_velocity(&self.network.candidates, self.contact.work_column())?,
        )?;
        self.point_kinematics(before.max(after), finite(increment / self.network.dt)?)?;
        let expected = mean_force(self.contact.source(), before, after)?;
        let constitutive_residual = finite(solution.force_n - expected)?;
        let force_tolerance = force_tolerance(solution.force_n, expected, self.config)?;
        if constitutive_residual.abs() > force_tolerance {
            return Err(ModalCouplingError::ContactSolve {
                residual_n: constitutive_residual,
                tolerance_n: force_tolerance,
                iterations: solution.iterations,
            }
            .into());
        }
        solve::positive_slope(
            self.contact.source(),
            before,
            after,
            self.compliance_m_per_n,
        )?;

        // Measure the represented force actually sent to each original mode.
        // If a large authored force erased a small contact correction, this
        // independent small-work gate refuses instead of crediting fictitious work.
        let mut actual_work = Dd::ZERO;
        for ((held, input), delta) in self
            .forces
            .iter()
            .zip(external)
            .zip(&self.network.free_delta)
        {
            actual_work = actual_work + Dd::from_f64(finite((held - input) * delta)?);
        }
        let contact_work = finite(actual_work.to_f64())?;
        let source_work = finite(expected * increment)?;
        let storage_change = finite(energy_after - energy_before)?;
        let contact_scale = energy_before
            .abs()
            .max(energy_after.abs())
            .max(contact_work.abs())
            .max(source_work.abs());
        let contact_tolerance = finite(
            self.config.contact_energy_absolute_tolerance_j
                + self.config.contact_energy_relative_tolerance * contact_scale,
        )?;
        check_contact_work(finite(contact_work + source_work)?, contact_tolerance)?;
        check_contact_work(finite(storage_change - source_work)?, contact_tolerance)?;
        let contact_residual = finite(storage_change + contact_work)?;
        check_contact_work(contact_residual, contact_tolerance)?;

        let frame = &self.network.staged_frame;
        let network_after = finite(frame.total_energy_j())?;
        let network_loss = finite(frame.total_dissipation_j())?;
        let external_work = dot(external, &self.network.free_delta)?;
        limit(
            "native/source energy magnitude",
            finite(network_after + energy_after)?.abs(),
            self.network.config.maximum_total_energy_j,
        )?;
        let residual = finite(
            (Dd::from_f64(network_after) - Dd::from_f64(network_before)
                + Dd::from_f64(storage_change)
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
        let accepted = SourceNormalModalFrame {
            sample: frame.sample,
            observer_pressure_pa: frame.observer_pressure_pa,
            network_energy_j: network_after,
            contact_energy_j: energy_after,
            external_work_j: external_work,
            contact_work_on_network_j: contact_work,
            network_dissipation_j: network_loss,
            energy_residual_j: residual,
            energy_tolerance_j: tolerance,
            contact_work_residual_j: contact_residual,
            contact_work_tolerance_j: contact_tolerance,
            normal_force_n: solution.force_n,
            endpoint_normal_force_n: endpoint_force,
            approach_before_m: before,
            approach_after_m: after,
            constitutive_residual_n: constitutive_residual,
            force_tolerance_n: force_tolerance,
            coordinate_residual_m: coordinate_residual,
            iterations: solution.iterations,
        };
        self.network.validate_staged_material_points(Some(gate))?;
        poll(Some(gate))?;
        Ok(accepted)
    }

    fn approach(&self, models: &[ModalAcousticTimeModel]) -> Result<f64, SourceNormalModalError> {
        Ok(extension(
            models,
            self.contact.work_column(),
            self.contact.gap_m(),
        )?)
    }
    fn check_approach(&self, approach: f64) -> Result<(), SourceNormalModalError> {
        if approach < self.config.minimum_approach_m || approach > self.config.maximum_approach_m {
            return Err(SourceNormalModalError::Unresolved(
                "native approach leaves the admitted source domain",
            ));
        }
        Ok(())
    }
    fn point_kinematics(&self, approach: f64, rate: f64) -> Result<(), SourceNormalModalError> {
        if let NormalLawSource::HertzPoint(source) = self.contact.source() {
            source
                .validate_kinematics(approach.max(0.0), rate, self.network.dt)
                .map_err(NormalLawSourceError::Point)?;
        }
        Ok(())
    }
}

fn observation(
    source: &NormalLawSource,
    approach: f64,
) -> Result<(f64, f64), SourceNormalModalError> {
    Ok(
        source.energy_and_force(if source.supports_signed_approach() {
            approach
        } else {
            approach.max(0.0)
        })?,
    )
}
fn mean_force(source: &NormalLawSource, from: f64, to: f64) -> Result<f64, SourceNormalModalError> {
    source
        .mean_force_between(from, to)?
        .ok_or(SourceNormalModalError::Unresolved(
            "original source path work unavailable",
        ))
}
fn force_tolerance(
    a: f64,
    b: f64,
    config: SourceNormalModalConfig,
) -> Result<f64, SourceNormalModalError> {
    Ok(finite(
        config.force_absolute_tolerance_n + config.force_relative_tolerance * a.abs().max(b.abs()),
    )?)
}
fn check_coordinate(
    a: f64,
    b: f64,
    config: SourceNormalModalConfig,
) -> Result<f64, SourceNormalModalError> {
    let residual = finite(a - b)?;
    let tolerance = finite(
        config.coordinate_absolute_tolerance_m
            + config.coordinate_relative_tolerance * a.abs().max(b.abs()),
    )?;
    if residual.abs() > tolerance {
        return Err(SourceNormalModalError::CoordinateMismatch {
            residual_m: residual,
            tolerance_m: tolerance,
        });
    }
    Ok(residual)
}
fn normal_velocity(
    models: &[ModalAcousticTimeModel],
    column: &[f64],
) -> Result<f64, ModalCouplingError> {
    models
        .iter()
        .flat_map(|model| model.states())
        .zip(column)
        .try_fold(0.0, |value, (state, weight)| {
            finite(value + state.velocity_m_sqrt_kg_per_s * weight)
        })
}
fn check_contact_work(residual_j: f64, tolerance_j: f64) -> Result<(), SourceNormalModalError> {
    if residual_j.abs() > tolerance_j {
        return Err(SourceNormalModalError::ContactWork {
            residual_j,
            tolerance_j,
        });
    }
    Ok(())
}
fn check_native_energy(residual_j: f64, tolerance_j: f64) -> Result<(), SourceNormalModalError> {
    if residual_j.abs() > tolerance_j {
        return Err(ModalCouplingError::EnergyBalance {
            residual_j,
            tolerance_j,
        }
        .into());
    }
    Ok(())
}
