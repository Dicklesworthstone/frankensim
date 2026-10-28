//! Simultaneous rolling loss, retained microslip and finite caloric feedback.
//!
//! The original rolling law is evaluated inside the SAME thermal endpoint as
//! the original partial-slip law. Its physical work heats the original cells
//! once; those temperatures immediately feed the retained friction/stiffness
//! model. This does not infer a rolling coefficient, effective temperature,
//! contact frequency, geometric speed mapping or measured calibration.
//!
//! The normal patch is held fixed. The caller supplies separate power-conjugate
//! rolling/contour motion, with a coordinate distinct from all microslip
//! coordinates. A selected actual surface supplies the rolling card's scalar
//! temperature; averaging unequal temperatures would require a different law.
//! Cards retain their original constant coefficients and applicability domains.
//! Endpoint time accuracy remains first order.

use super::*;
use fs_tribo::partial_slip::NormalPatchAuthority;
use fs_tribo::rolling_loss::{
    RollingKinematics, RollingLossChannel, RollingLossCheckpoint, RollingLossError, RollingLossLaw,
    RollingLossState, RollingLossStep, RollingPatchReceipt, RollingWorkOwnership,
};
use fs_tribo::{InterfaceMedium, InterfaceSystemRef};

/// The actual temperature coordinate required by the declared rolling card.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RollingTemperatureSurface {
    SurfaceA,
    SurfaceB,
}
impl RollingTemperatureSurface {
    const fn index(self) -> usize {
        match self {
            Self::SurfaceA => 0,
            Self::SurfaceB => 1,
        }
    }
}

/// Signed native generalized speeds, not inferred from microslip or spin.
/// The outer motion owner must supply their physical geometric mapping.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RollingMotion {
    pub contour_speed_mps: f64,
    pub rolling_rate_rad_s: f64,
    pub excitation_frequency_hz: f64,
}

/// Failures of the extra physical source evaluated inside the existing solve.
#[derive(Debug)]
pub enum RollingEndpointError<E = Infallible> {
    Native(RollingLossError),
    Invalid(&'static str),
    Thermal(DepositionError),
    Source(E),
}

/// A refusal never changes either original mechanical history or caloric state.
#[derive(Debug)]
pub enum CaloricRollingError<E = Infallible> {
    Invalid(&'static str),
    Native(RollingLossError),
    Contact(CaloricContactError<RollingEndpointError<E>>),
}

/// The original rolling result at the accepted nonlinear thermal endpoint.
#[derive(Debug)]
pub struct RollingThermalEndpoint<S> {
    pub rolling: RollingLossStep,
    pub temperature_k: f64,
    pub deposited_heat_j: [f64; 2],
    /// Heat from the caller's other sources/boundaries before adding rolling.
    pub external_heat_j: f64,
    pub source: S,
}

/// Nested contact accounting and the independently admitted combined account.
pub struct CaloricRollingReport<S> {
    /// Its thermal subsystem receives rolling heat as a source, so this INNER
    /// external_heat_j includes rolling work. Use the outer external_heat_j
    /// below for the complete rolling + microslip + finite-store balance.
    pub partial_slip: CaloricContactReport<RollingThermalEndpoint<S>>,
    pub total_mechanical_work_j: f64,
    pub external_heat_j: f64,
    /// Delta H + Delta U_contact + D_numerical - W_microslip - W_rolling - Q_other.
    pub total_energy_residual_j: f64,
}
impl<S> CaloricRollingReport<S> {
    pub const fn rolling(&self) -> &RollingLossStep {
        &self.partial_slip.source.rolling
    }
    pub const fn source(&self) -> &S {
        &self.partial_slip.source.source
    }
}

/// Both original histories and their shared thermal clock, accepted together.
#[derive(Clone)]
pub struct CaloricRollingContact<'a> {
    partial_slip: CaloricPartialSlip<'a>,
    patch: RollingPatchReceipt,
    interface: InterfaceSystemRef,
    law: RollingLossLaw,
    state: RollingLossState,
    checkpoint: Option<RollingLossCheckpoint>,
    ownership: RollingWorkOwnership,
    temperature_surface: RollingTemperatureSurface,
    partition: HeatPartition,
}

pub struct CaloricRollingCandidate<'a, 't, S> {
    thermal: DepositionSimulation<'t>,
    contact: CaloricRollingContact<'a>,
    report: CaloricRollingReport<S>,
}
impl<'a, 't, S> CaloricRollingCandidate<'a, 't, S> {
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn contact(&self) -> &CaloricRollingContact<'a> {
        &self.contact
    }
    pub const fn report(&self) -> &CaloricRollingReport<S> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        CaloricRollingContact<'a>,
        CaloricRollingReport<S>,
    ) {
        (self.thermal, self.contact, self.report)
    }
}

impl<'a> CaloricRollingContact<'a> {
    /// Attach an original rolling channel without resetting either history.
    /// A nonzero rolling history requires its original native checkpoint.
    /// The two native patches must describe the same held elliptical footprint.
    /// Heat allocation is explicit data, independent of microslip allocation;
    /// no loss is exported outside these two finite destinations.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        partial_slip: CaloricPartialSlip<'a>,
        patch: RollingPatchReceipt,
        interface: InterfaceSystemRef,
        law: RollingLossLaw,
        state: RollingLossState,
        checkpoint: Option<RollingLossCheckpoint>,
        ownership: RollingWorkOwnership,
        temperature_surface: RollingTemperatureSurface,
        partition: HeatPartition,
    ) -> Result<Self, CaloricRollingError> {
        let normal = partial_slip.checkpoint.normal_patch();
        let original = partial_slip.checkpoint.interface();
        let area = core::f64::consts::PI
            * normal.semi_axis_longitudinal_m()
            * normal.semi_axis_lateral_m();
        let same_authority = |a, b| {
            matches!(
                (a, b),
                (
                    NormalPatchAuthority::CallerDeclared,
                    InputAuthority::CallerDeclared
                ) | (
                    NormalPatchAuthority::SyntheticFixture,
                    InputAuthority::SyntheticFixture
                ) | (NormalPatchAuthority::Estimated, InputAuthority::Estimated)
            )
        };
        if patch.patch_id() != normal.patch_id()
            || patch.normal_model_id() != normal.normal_model_id()
            || patch.source_id() != normal.source_id()
            || !same_authority(normal.authority(), patch.authority())
            || patch.normal_load_n() != normal.normal_load_n()
            || !area.is_finite()
            || (patch.contact_area_m2() - area).abs() > 32.0 * f64::EPSILON * area
            || interface.ordered_system_id() != original.ordered_interface_id()
            || interface.history_id() != original.history_id()
            || interface.provenance().source_id() != original.source_id()
            || !same_authority(original.authority(), interface.provenance().authority())
            || interface.medium() != InterfaceMedium::Dry
            || partition.other() != 0.0
            || ownership.patch_id() != patch.patch_id()
        {
            return Err(CaloricRollingError::Invalid(
                "rolling needs the same normal patch, dry ordered interface and finite destinations",
            ));
        }
        let (channel, domain) = match &law {
            RollingLossLaw::CoulombContour(c) => {
                (RollingLossChannel::ContourDeformation, c.applicability())
            }
            RollingLossLaw::ViscousContour(c) => {
                (RollingLossChannel::ContourDeformation, c.applicability())
            }
            RollingLossLaw::HystereticRollingMoment(c) => {
                (RollingLossChannel::RollingHysteresis, c.applicability())
            }
        };
        let temperature = partial_slip.temperature_k[temperature_surface.index()];
        if ownership.channel() != channel
            || temperature < domain.temperature_kelvin().minimum()
            || temperature > domain.temperature_kelvin().maximum()
        {
            return Err(CaloricRollingError::Invalid(
                "rolling channel or original surface temperature is outside its card",
            ));
        }
        if let Some(checkpoint) = &checkpoint {
            let restored = law
                .restore_checkpoint(&patch, &interface, &ownership, checkpoint)
                .map_err(CaloricRollingError::Native)?;
            if restored != state {
                return Err(CaloricRollingError::Invalid(
                    "rolling history differs from its original checkpoint",
                ));
            }
        } else if state != RollingLossState::zero() {
            return Err(CaloricRollingError::Invalid(
                "nonzero rolling history requires its original checkpoint",
            ));
        }
        Ok(Self {
            partial_slip,
            patch,
            interface,
            law,
            state,
            checkpoint,
            ownership,
            temperature_surface,
            partition,
        })
    }

    pub const fn partial_slip(&self) -> &CaloricPartialSlip<'a> {
        &self.partial_slip
    }
    pub const fn rolling_state(&self) -> &RollingLossState {
        &self.state
    }
    pub const fn rolling_checkpoint(&self) -> Option<&RollingLossCheckpoint> {
        self.checkpoint.as_ref()
    }
    pub const fn law(&self) -> &RollingLossLaw {
        &self.law
    }

    /// One common endpoint solves both losses and finite thermal storage.
    /// The rolling work owner adopts the microslip interval identity and keeps
    /// its distinct physical coordinate. Every trial reuses the accepted native
    /// histories, never a previous Newton trial. Commit all returned owners.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn step_candidate<'t, S, E>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        frame: TangentFrame,
        kinematics: PartialSlipKinematics,
        ownership: &GeneralizedWorkOwnership,
        rolling_motion: RollingMotion,
        controls: CaloricContactSolve,
        guess_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'t>) -> Result<StepCoefficients<S>, E>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<CaloricRollingCandidate<'a, 't, S>, CaloricRollingError<E>> {
        if ownership.patch_id() != self.patch.patch_id()
            || [
                ownership.longitudinal_coordinate_id(),
                ownership.lateral_coordinate_id(),
                ownership.torsional_coordinate_id(),
            ]
            .contains(&self.ownership.generalized_coordinate_id())
        {
            return Err(CaloricRollingError::Invalid(
                "rolling work must own a distinct generalized coordinate on the same patch",
            ));
        }
        let mut rolling_owner = self.ownership.clone();
        rolling_owner
            .retarget(
                self.patch.patch_id(),
                ownership.interval_id(),
                self.ownership.generalized_coordinate_id(),
                self.ownership.channel(),
            )
            .map_err(CaloricRollingError::Native)?;
        if let Some(checkpoint) = &self.checkpoint {
            let restored = self
                .law
                .restore_checkpoint(&self.patch, &self.interface, &rolling_owner, checkpoint)
                .map_err(CaloricRollingError::Native)?;
            if restored != self.state {
                return Err(CaloricRollingError::Invalid("stale native rolling history"));
            }
        }
        let cancel = RefCell::new(cancelled);
        let poll = || {
            if (*cancel.borrow_mut())() {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        };
        let dt = target_s - self.partial_slip.time_s;
        let proposal = self
            .partial_slip
            .step_candidate(
                thermal,
                target_s,
                frame,
                kinematics,
                ownership,
                controls,
                guess_j,
                |candidate| {
                    let temperature = candidate
                        .temperature_k(self.partial_slip.cells[self.temperature_surface.index()])
                        .map_err(RollingEndpointError::Thermal)?;
                    let motion = RollingKinematics::new(
                        rolling_motion.contour_speed_mps,
                        rolling_motion.rolling_rate_rad_s,
                        kinematics.torsional_spin_rad_per_s,
                        temperature,
                        rolling_motion.excitation_frequency_hz,
                        dt,
                    )
                    .map_err(RollingEndpointError::Native)?;
                    let rolling = self
                        .law
                        .advance(
                            &self.patch,
                            &self.interface,
                            motion,
                            &rolling_owner,
                            &self.state,
                        )
                        .map_err(RollingEndpointError::Native)?;
                    let work = -(rolling.wrench.contour_force_n * rolling_motion.contour_speed_mps
                        + rolling.wrench.rolling_moment_nm * rolling_motion.rolling_rate_rad_s)
                        * dt;
                    let heat = rolling.dissipation.total_heat_j;
                    if !work.is_finite()
                        || !heat.is_finite()
                        || work < 0.0
                        || heat < 0.0
                        || rolling.stored_energy_j != 0.0
                        || rolling.storage_change_j != 0.0
                        || (work - heat).abs() > controls.total_energy_tolerance_j
                        || (work - rolling.generalized_work.work_into_interface_j).abs()
                            > controls.total_energy_tolerance_j
                    {
                        return Err(RollingEndpointError::Invalid(
                            "native rolling work and irreversible heat do not close",
                        ));
                    }
                    let StepCoefficients {
                        mut sources_w,
                        conductances_w_k,
                        payload,
                    } = resolve(candidate).map_err(RollingEndpointError::Source)?;
                    let external_heat_j = dt
                        * candidate
                            .external_step_power_w(&sources_w, &conductances_w_k, poll)
                            .map_err(RollingEndpointError::Thermal)?;
                    let deposited_heat_j = [
                        self.partition.surface_a() * heat,
                        self.partition.surface_b() * heat,
                    ];
                    if !external_heat_j.is_finite()
                        || deposited_heat_j.iter().any(|q| !q.is_finite() || *q < 0.0)
                        || (deposited_heat_j[0] + deposited_heat_j[1] - heat).abs()
                            > controls.total_energy_tolerance_j
                    {
                        return Err(RollingEndpointError::Invalid(
                            "invalid rolling heat allocation",
                        ));
                    }
                    for (id, q) in self.partial_slip.cells.into_iter().zip(deposited_heat_j) {
                        let p = q / dt;
                        let source = sources_w.entry(id).or_insert(0.0);
                        let total = *source + p;
                        if !p.is_finite() || !total.is_finite() || (q != 0.0 && p == 0.0) {
                            return Err(RollingEndpointError::Invalid(
                                "unrepresentable rolling power",
                            ));
                        }
                        *source = total;
                    }
                    Ok(StepCoefficients {
                        sources_w,
                        conductances_w_k,
                        payload: RollingThermalEndpoint {
                            rolling,
                            temperature_k: temperature,
                            deposited_heat_j,
                            external_heat_j,
                            source: payload,
                        },
                    })
                },
                || (*cancel.borrow_mut())(),
            )
            .map_err(CaloricRollingError::Contact)?;
        let (thermal, partial_slip, report) = proposal.into_parts();
        let rolling = &report.source.rolling;
        let external_heat_j = report.source.external_heat_j;
        let total_mechanical_work_j = report.contact.generalized_work.work_into_interface_j
            + rolling.generalized_work.work_into_interface_j;
        let total_energy_residual_j = report.thermal.stored_energy_change_j
            + report.contact_internal_energy_change_j
            + report.contact.generalized_work.numerical_dissipation_j
            - total_mechanical_work_j
            - external_heat_j;
        if !total_mechanical_work_j.is_finite()
            || !total_energy_residual_j.is_finite()
            || total_energy_residual_j.abs() > controls.total_energy_tolerance_j
        {
            return Err(CaloricRollingError::Invalid(
                "rolling, retained contact and finite thermal energy do not balance",
            ));
        }
        poll().map_err(|e| CaloricRollingError::Contact(CaloricContactError::from(e)))?;
        let mut contact = self.clone();
        contact.partial_slip = partial_slip;
        contact.state = rolling.next_state.clone();
        contact.checkpoint = Some(rolling.checkpoint.clone());
        contact.ownership = rolling_owner;
        Ok(CaloricRollingCandidate {
            thermal,
            contact,
            report: CaloricRollingReport {
                partial_slip: report,
                total_mechanical_work_j,
                external_heat_j,
                total_energy_residual_j,
            },
        })
    }
}
