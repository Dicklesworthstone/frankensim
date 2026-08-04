//! Finite-patch rolling and contour-deformation loss candidates.
//!
//! This module keeps three rival rungs separate: a Leine-style Coulomb contour
//! force, a viscous contour force, and a finite-patch hysteretic rolling moment.
//! They retain their source cards and are never selected or blended from an
//! apparatus response. Tangential creepage, torsional microslip, stick, and
//! gross slide are intentionally outside this module.

use core::fmt;

use crate::partial_slip::GeneralizedWorkOwnership;
use crate::{ApplicabilityRange, InputAuthority, InterfaceMedium, InterfaceSystemRef};

const EPSILON: f64 = 128.0 * f64::EPSILON;

/// Generic published-form identity retained by the Coulomb contour rung.
///
/// This identifies the Leine-style contour-force form used here; it is not a
/// claim that the form dominates any particular material system.
pub const LEINE_STYLE_CONTOUR_LAW_ID: &str = "Leine-style Coulomb contour-force rung";

/// Refusal surface for rolling-loss inputs and replay receipts.
#[derive(Debug, Clone, PartialEq)]
pub enum RollingLossError {
    /// A required caller identity was blank.
    MissingIdentity { field: &'static str },
    /// A value was nonfinite or outside this module's physical input domain.
    InvalidInput { field: &'static str },
    /// The declared interface is not dry.
    NotDryInterface { medium: InterfaceMedium },
    /// Temperature or excitation frequency is outside a source card's range.
    OutsideApplicability {
        /// Named quantity that failed the card range.
        field: &'static str,
        /// Submitted value in the documented SI unit.
        value: f64,
        /// Inclusive card lower bound.
        minimum: f64,
        /// Inclusive card upper bound.
        maximum: f64,
    },
    /// A rolling work key does not bind the patch or required loss channel.
    WorkOwnershipMismatch { field: &'static str },
    /// A rolling work interval overlaps a partial-slip interval.
    WorkOwnershipOverlap,
    /// A checkpoint does not bind the supplied receipt, card, or law.
    CheckpointMismatch { field: &'static str },
}

impl fmt::Display for RollingLossError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingIdentity { field } => write!(formatter, "missing identity: {field}"),
            Self::InvalidInput { field } => {
                write!(formatter, "invalid rolling-loss input: {field}")
            }
            Self::NotDryInterface { medium } => write!(
                formatter,
                "rolling loss requires dry medium, got {medium:?}"
            ),
            Self::OutsideApplicability {
                field,
                value,
                minimum,
                maximum,
            } => write!(
                formatter,
                "{field}={value} lies outside declared applicability [{minimum}, {maximum}]"
            ),
            Self::WorkOwnershipMismatch { field } => {
                write!(formatter, "rolling work ownership mismatch: {field}")
            }
            Self::WorkOwnershipOverlap => {
                write!(formatter, "rolling work overlaps the partial-slip interval")
            }
            Self::CheckpointMismatch { field } => write!(formatter, "checkpoint mismatch: {field}"),
        }
    }
}

impl std::error::Error for RollingLossError {}

/// Whether curvature is retained as two principal curvatures or as an explicit
/// caller approximation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PatchCurvature {
    /// Signed principal curvatures [1/m] in a caller-declared patch frame.
    Principal {
        /// First signed principal curvature [1/m].
        first_per_m: f64,
        /// Second signed principal curvature [1/m].
        second_per_m: f64,
    },
    /// A positive equivalent radius [m] substituted for principal curvatures.
    EquivalentRadiusApproximation {
        /// Equivalent radius [m].
        radius_m: f64,
        /// Caller ceiling for the approximation; this module does not upgrade it.
        authority: InputAuthority,
    },
}

impl PatchCurvature {
    fn validate(self) -> Result<(), RollingLossError> {
        match self {
            Self::Principal {
                first_per_m,
                second_per_m,
            } => {
                finite(first_per_m, "first_principal_curvature_per_m")?;
                finite(second_per_m, "second_principal_curvature_per_m")
            }
            Self::EquivalentRadiusApproximation { radius_m, .. } => {
                positive(radius_m, "equivalent_radius_m")
            }
        }
    }

    fn approximation_authority(self) -> Option<InputAuthority> {
        match self {
            Self::Principal { .. } => None,
            Self::EquivalentRadiusApproximation { authority, .. } => Some(authority),
        }
    }
}

/// Immutable finite normal-patch receipt used by the loss rungs.
///
/// It retains the caller's source and authority ceiling. It is not an
/// independent normal-contact admission or a pressure-field solution.
#[derive(Debug, Clone, PartialEq)]
pub struct RollingPatchReceipt {
    patch_id: String,
    normal_model_id: String,
    source_id: String,
    authority: InputAuthority,
    normal_load_n: f64,
    contact_area_m2: f64,
    curvature: PatchCurvature,
}

impl RollingPatchReceipt {
    /// Builds a named finite-patch receipt. A zero load may retain zero area;
    /// a positive load requires a positive contact area.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        patch_id: impl Into<String>,
        normal_model_id: impl Into<String>,
        source_id: impl Into<String>,
        authority: InputAuthority,
        normal_load_n: f64,
        contact_area_m2: f64,
        curvature: PatchCurvature,
    ) -> Result<Self, RollingLossError> {
        let value = Self {
            patch_id: patch_id.into(),
            normal_model_id: normal_model_id.into(),
            source_id: source_id.into(),
            authority,
            normal_load_n,
            contact_area_m2,
            curvature,
        };
        value.validate()?;
        Ok(value)
    }

    /// Caller patch identity.
    #[must_use]
    pub fn patch_id(&self) -> &str {
        &self.patch_id
    }

    /// Caller normal-response model identity.
    #[must_use]
    pub fn normal_model_id(&self) -> &str {
        &self.normal_model_id
    }

    /// Caller source identity for the patch receipt.
    #[must_use]
    pub fn source_id(&self) -> &str {
        &self.source_id
    }

    /// Caller-declared authority ceiling.
    #[must_use]
    pub const fn authority(&self) -> InputAuthority {
        self.authority
    }

    /// Compressive normal load [N].
    #[must_use]
    pub const fn normal_load_n(&self) -> f64 {
        self.normal_load_n
    }

    /// Retained finite contact area [m²].
    #[must_use]
    pub const fn contact_area_m2(&self) -> f64 {
        self.contact_area_m2
    }

    /// Principal-curvature receipt or an explicit equivalent-radius approximation.
    #[must_use]
    pub const fn curvature(&self) -> PatchCurvature {
        self.curvature
    }

    fn validate(&self) -> Result<(), RollingLossError> {
        nonblank(&self.patch_id, "patch_id")?;
        nonblank(&self.normal_model_id, "normal_model_id")?;
        nonblank(&self.source_id, "patch_source_id")?;
        nonnegative(self.normal_load_n, "normal_load_n")?;
        nonnegative(self.contact_area_m2, "contact_area_m2")?;
        if self.normal_load_n > 0.0 {
            positive(self.contact_area_m2, "contact_area_m2")?;
        }
        self.curvature.validate()
    }
}

/// Closed temperature/frequency validity ranges on an immutable source card.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RollingLossApplicability {
    temperature_kelvin: ApplicabilityRange,
    excitation_frequency_hz: ApplicabilityRange,
}

impl RollingLossApplicability {
    /// Builds a dry-loss validity domain. Absolute temperature needs a positive
    /// lower bound; excitation frequency may include zero.
    pub fn new(
        temperature_kelvin: ApplicabilityRange,
        excitation_frequency_hz: ApplicabilityRange,
    ) -> Result<Self, RollingLossError> {
        if temperature_kelvin.minimum() == 0.0 {
            return Err(RollingLossError::InvalidInput {
                field: "temperature_kelvin.minimum",
            });
        }
        Ok(Self {
            temperature_kelvin,
            excitation_frequency_hz,
        })
    }

    /// Declared absolute-temperature interval [K].
    #[must_use]
    pub const fn temperature_kelvin(&self) -> ApplicabilityRange {
        self.temperature_kelvin
    }

    /// Declared excitation-frequency interval [Hz].
    #[must_use]
    pub const fn excitation_frequency_hz(&self) -> ApplicabilityRange {
        self.excitation_frequency_hz
    }

    fn require_contains(self, kinematics: RollingKinematics) -> Result<(), RollingLossError> {
        positive(kinematics.temperature_kelvin, "temperature_kelvin")?;
        require_range(
            self.temperature_kelvin,
            kinematics.temperature_kelvin,
            "temperature_kelvin",
        )?;
        require_range(
            self.excitation_frequency_hz,
            kinematics.excitation_frequency_hz,
            "excitation_frequency_hz",
        )
    }
}

/// Immutable source card for the Leine-style Coulomb contour-force rung.
#[derive(Debug, Clone, PartialEq)]
pub struct CoulombContourCard {
    model_id: String,
    source_id: String,
    authority: InputAuthority,
    contour_coefficient: f64,
    applicability: RollingLossApplicability,
}

impl CoulombContourCard {
    /// Creates `F_c = -sign(v_c) c_c N`, with dimensionless `c_c`.
    pub fn new(
        model_id: impl Into<String>,
        source_id: impl Into<String>,
        authority: InputAuthority,
        contour_coefficient: f64,
        applicability: RollingLossApplicability,
    ) -> Result<Self, RollingLossError> {
        let value = Self {
            model_id: model_id.into(),
            source_id: source_id.into(),
            authority,
            contour_coefficient,
            applicability,
        };
        value.validate()?;
        Ok(value)
    }

    /// Stable source-model identity.
    #[must_use]
    pub fn model_id(&self) -> &str {
        &self.model_id
    }

    /// Caller card source identity.
    #[must_use]
    pub fn source_id(&self) -> &str {
        &self.source_id
    }

    /// Caller authority ceiling.
    #[must_use]
    pub const fn authority(&self) -> InputAuthority {
        self.authority
    }

    /// Dimensionless contour coefficient.
    #[must_use]
    pub const fn contour_coefficient(&self) -> f64 {
        self.contour_coefficient
    }

    /// Card validity ranges.
    #[must_use]
    pub const fn applicability(&self) -> RollingLossApplicability {
        self.applicability
    }

    fn validate(&self) -> Result<(), RollingLossError> {
        nonblank(&self.model_id, "coulomb_model_id")?;
        nonblank(&self.source_id, "coulomb_source_id")?;
        nonnegative(self.contour_coefficient, "contour_coefficient")
    }
}

/// Immutable source card for the viscous contour-force alternative.
#[derive(Debug, Clone, PartialEq)]
pub struct ViscousContourCard {
    model_id: String,
    source_id: String,
    authority: InputAuthority,
    viscous_coefficient_n_s_per_m: f64,
    applicability: RollingLossApplicability,
}

impl ViscousContourCard {
    /// Creates `F_v = -c_v v_c`, where `c_v` has units N s/m.
    pub fn new(
        model_id: impl Into<String>,
        source_id: impl Into<String>,
        authority: InputAuthority,
        viscous_coefficient_n_s_per_m: f64,
        applicability: RollingLossApplicability,
    ) -> Result<Self, RollingLossError> {
        let value = Self {
            model_id: model_id.into(),
            source_id: source_id.into(),
            authority,
            viscous_coefficient_n_s_per_m,
            applicability,
        };
        value.validate()?;
        Ok(value)
    }

    /// Stable source-model identity.
    #[must_use]
    pub fn model_id(&self) -> &str {
        &self.model_id
    }

    /// Caller card source identity.
    #[must_use]
    pub fn source_id(&self) -> &str {
        &self.source_id
    }

    /// Caller authority ceiling.
    #[must_use]
    pub const fn authority(&self) -> InputAuthority {
        self.authority
    }

    /// Viscous contour coefficient [N s/m].
    #[must_use]
    pub const fn viscous_coefficient_n_s_per_m(&self) -> f64 {
        self.viscous_coefficient_n_s_per_m
    }

    /// Card validity ranges.
    #[must_use]
    pub const fn applicability(&self) -> RollingLossApplicability {
        self.applicability
    }

    fn validate(&self) -> Result<(), RollingLossError> {
        nonblank(&self.model_id, "viscous_model_id")?;
        nonblank(&self.source_id, "viscous_source_id")?;
        nonnegative(
            self.viscous_coefficient_n_s_per_m,
            "viscous_coefficient_n_s_per_m",
        )
    }
}

/// Immutable finite-patch source card for the hysteretic rolling-moment candidate.
#[derive(Debug, Clone, PartialEq)]
pub struct HystereticRollingCard {
    model_id: String,
    source_id: String,
    authority: InputAuthority,
    loss_length_m: f64,
    loss_factor: f64,
    applicability: RollingLossApplicability,
}

impl HystereticRollingCard {
    /// Creates `M_h = -sign(omega_r) N l_h eta_h`, where `l_h` is a declared
    /// loss length [m] and `0 <= eta_h <= 1` is a dimensionless loss factor.
    pub fn new(
        model_id: impl Into<String>,
        source_id: impl Into<String>,
        authority: InputAuthority,
        loss_length_m: f64,
        loss_factor: f64,
        applicability: RollingLossApplicability,
    ) -> Result<Self, RollingLossError> {
        let value = Self {
            model_id: model_id.into(),
            source_id: source_id.into(),
            authority,
            loss_length_m,
            loss_factor,
            applicability,
        };
        value.validate()?;
        Ok(value)
    }

    /// Stable source-model identity.
    #[must_use]
    pub fn model_id(&self) -> &str {
        &self.model_id
    }

    /// Caller card source identity.
    #[must_use]
    pub fn source_id(&self) -> &str {
        &self.source_id
    }

    /// Caller authority ceiling.
    #[must_use]
    pub const fn authority(&self) -> InputAuthority {
        self.authority
    }

    /// Declared hysteretic loss length [m].
    #[must_use]
    pub const fn loss_length_m(&self) -> f64 {
        self.loss_length_m
    }

    /// Declared dimensionless loss factor.
    #[must_use]
    pub const fn loss_factor(&self) -> f64 {
        self.loss_factor
    }

    /// Card validity ranges.
    #[must_use]
    pub const fn applicability(&self) -> RollingLossApplicability {
        self.applicability
    }

    fn validate(&self) -> Result<(), RollingLossError> {
        nonblank(&self.model_id, "hysteretic_model_id")?;
        nonblank(&self.source_id, "hysteretic_source_id")?;
        nonnegative(self.loss_length_m, "loss_length_m")?;
        nonnegative(self.loss_factor, "loss_factor")?;
        if self.loss_factor > 1.0 {
            return Err(RollingLossError::InvalidInput {
                field: "loss_factor",
            });
        }
        Ok(())
    }
}

/// Separately parameterized rolling/contour loss rungs.
///
/// This enum deliberately has no blending or preference operation. Consumers
/// must preserve model disagreement instead of combining candidate outputs.
#[derive(Debug, Clone, PartialEq)]
pub enum RollingLossLaw {
    /// Leine-style Coulomb contour force opposing material contour speed.
    CoulombContour(CoulombContourCard),
    /// Linear viscous contour-force alternative.
    ViscousContour(ViscousContourCard),
    /// Finite-patch hysteretic rolling-moment candidate.
    HystereticRollingMoment(HystereticRollingCard),
}

impl RollingLossLaw {
    /// Card model identity for this distinct candidate rung.
    #[must_use]
    pub fn model_id(&self) -> &str {
        match self {
            Self::CoulombContour(card) => card.model_id(),
            Self::ViscousContour(card) => card.model_id(),
            Self::HystereticRollingMoment(card) => card.model_id(),
        }
    }

    /// Caller source identity for this distinct candidate rung.
    #[must_use]
    pub fn source_id(&self) -> &str {
        match self {
            Self::CoulombContour(card) => card.source_id(),
            Self::ViscousContour(card) => card.source_id(),
            Self::HystereticRollingMoment(card) => card.source_id(),
        }
    }

    /// Caller authority ceiling.
    #[must_use]
    pub const fn authority(&self) -> InputAuthority {
        match self {
            Self::CoulombContour(card) => card.authority(),
            Self::ViscousContour(card) => card.authority(),
            Self::HystereticRollingMoment(card) => card.authority(),
        }
    }

    fn applicability(&self) -> RollingLossApplicability {
        match self {
            Self::CoulombContour(card) => card.applicability(),
            Self::ViscousContour(card) => card.applicability(),
            Self::HystereticRollingMoment(card) => card.applicability(),
        }
    }

    fn required_channel(&self) -> RollingLossChannel {
        match self {
            Self::CoulombContour(_) | Self::ViscousContour(_) => {
                RollingLossChannel::ContourDeformation
            }
            Self::HystereticRollingMoment(_) => RollingLossChannel::RollingHysteresis,
        }
    }

    fn validate(&self) -> Result<(), RollingLossError> {
        match self {
            Self::CoulombContour(card) => card.validate(),
            Self::ViscousContour(card) => card.validate(),
            Self::HystereticRollingMoment(card) => card.validate(),
        }
    }

    /// Evaluates one accepted loss interval. The returned state/checkpoint is a
    /// candidate that the caller commits only after its outer solver accepts it.
    pub fn advance(
        &self,
        patch: &RollingPatchReceipt,
        interface: &InterfaceSystemRef,
        kinematics: RollingKinematics,
        ownership: &RollingWorkOwnership,
        state: &RollingLossState,
    ) -> Result<RollingLossStep, RollingLossError> {
        patch.validate()?;
        validate_interface(interface)?;
        self.validate()?;
        self.applicability().require_contains(kinematics)?;
        ownership.validate_for(patch.patch_id(), self.required_channel())?;
        state.validate()?;

        let (state_kind, contour_force_n, rolling_moment_nm) = match self {
            Self::CoulombContour(card) => {
                let magnitude = checked_mul(
                    card.contour_coefficient,
                    patch.normal_load_n,
                    "coulomb_contour_force_n",
                )?;
                (
                    if magnitude == 0.0 || kinematics.contour_speed_mps == 0.0 {
                        RollingLossStateKind::Quiescent
                    } else {
                        RollingLossStateKind::ContourDeformation
                    },
                    -kinematics.contour_speed_mps.signum() * magnitude,
                    0.0,
                )
            }
            Self::ViscousContour(card) => {
                let force = -checked_mul(
                    card.viscous_coefficient_n_s_per_m,
                    kinematics.contour_speed_mps,
                    "viscous_contour_force_n",
                )?;
                (
                    if force == 0.0 {
                        RollingLossStateKind::Quiescent
                    } else {
                        RollingLossStateKind::ContourDeformation
                    },
                    force,
                    0.0,
                )
            }
            Self::HystereticRollingMoment(card) => {
                let magnitude = checked_mul(
                    checked_mul(
                        patch.normal_load_n,
                        card.loss_length_m,
                        "hysteretic_rolling_moment_nm",
                    )?,
                    card.loss_factor,
                    "hysteretic_rolling_moment_nm",
                )?;
                (
                    if magnitude == 0.0 || kinematics.rolling_rate_rad_s == 0.0 {
                        RollingLossStateKind::Quiescent
                    } else {
                        RollingLossStateKind::RollingHysteresis
                    },
                    0.0,
                    -kinematics.rolling_rate_rad_s.signum() * magnitude,
                )
            }
        };
        let endpoint_body_power_w = checked_add(
            checked_mul(
                contour_force_n,
                kinematics.contour_speed_mps,
                "contour_power_w",
            )?,
            checked_mul(
                rolling_moment_nm,
                kinematics.rolling_rate_rad_s,
                "rolling_power_w",
            )?,
            "endpoint_body_power_w",
        )?;
        if endpoint_body_power_w > EPSILON {
            return Err(RollingLossError::InvalidInput {
                field: "non_passive_endpoint_power_w",
            });
        }
        let heat_j = checked_mul(
            -endpoint_body_power_w,
            kinematics.interval_s,
            "rolling_loss_heat_j",
        )?;
        nonnegative(heat_j, "rolling_loss_heat_j")?;
        let next_state = RollingLossState {
            cumulative_heat_j: checked_add(state.cumulative_heat_j, heat_j, "cumulative_heat_j")?,
            accepted_steps: state.accepted_steps.checked_add(1).ok_or(
                RollingLossError::InvalidInput {
                    field: "accepted_steps",
                },
            )?,
        };
        next_state.validate()?;
        let reconstructed_body_power_w = if kinematics.interval_s == 0.0 {
            0.0
        } else {
            -heat_j / kinematics.interval_s
        };
        finite(reconstructed_body_power_w, "reconstructed_body_power_w")?;
        let checkpoint = RollingLossCheckpoint {
            patch: patch.clone(),
            interface: interface.clone(),
            law: self.clone(),
            ownership: ownership.clone(),
            state: next_state.clone(),
        };
        Ok(RollingLossStep {
            state: state_kind,
            wrench: RollingLossWrench {
                contour_force_n,
                rolling_moment_nm,
                spin_moment_nm: 0.0,
            },
            stored_energy_j: 0.0,
            storage_change_j: 0.0,
            dissipation: RollingLossDissipation {
                contour_deformation_heat_j: if matches!(
                    self,
                    Self::CoulombContour(_) | Self::ViscousContour(_)
                ) {
                    heat_j
                } else {
                    0.0
                },
                rolling_hysteresis_heat_j: if matches!(self, Self::HystereticRollingMoment(_)) {
                    heat_j
                } else {
                    0.0
                },
                total_heat_j: heat_j,
            },
            generalized_work: RollingGeneralizedWork {
                ownership: ownership.clone(),
                work_into_interface_j: heat_j,
                endpoint_body_power_w,
                reconstructed_body_power_w,
            },
            next_state,
            checkpoint,
            applicability: self.applicability_kind(),
            uncertainty: RollingLossUncertainty::ModelFormNoNumericalCertificate {
                curvature_approximation_authority: patch.curvature.approximation_authority(),
            },
            patch_authority: patch.authority,
            interface_authority: interface.provenance().authority(),
            card_authority: self.authority(),
        })
    }

    /// Restores a candidate state only when stable receipt lineage and source
    /// card identity match the checkpoint.
    ///
    /// Normal load, contact area, and curvature are resolved for the current
    /// step by [`Self::advance`], and an ownership interval names that current
    /// work interval. They therefore must not freeze a continuous rolling
    /// trajectory at its first accepted step. Patch lineage (patch/model/source/
    /// authority), the complete interface, law card, and work coordinate/channel
    /// remain bound so a changed material, law, interface, or work channel still
    /// refuses deterministically.
    pub fn restore_checkpoint(
        &self,
        patch: &RollingPatchReceipt,
        interface: &InterfaceSystemRef,
        ownership: &RollingWorkOwnership,
        checkpoint: &RollingLossCheckpoint,
    ) -> Result<RollingLossState, RollingLossError> {
        patch.validate()?;
        validate_interface(interface)?;
        self.validate()?;
        ownership.validate_for(patch.patch_id(), self.required_channel())?;
        checkpoint.validate()?;
        if checkpoint.law != *self {
            return Err(RollingLossError::CheckpointMismatch { field: "law" });
        }
        if checkpoint.patch.patch_id != patch.patch_id
            || checkpoint.patch.normal_model_id != patch.normal_model_id
            || checkpoint.patch.source_id != patch.source_id
            || checkpoint.patch.authority != patch.authority
        {
            return Err(RollingLossError::CheckpointMismatch { field: "patch" });
        }
        if checkpoint.interface != *interface {
            return Err(RollingLossError::CheckpointMismatch { field: "interface" });
        }
        if checkpoint.ownership.patch_id != ownership.patch_id
            || checkpoint.ownership.generalized_coordinate_id != ownership.generalized_coordinate_id
            || checkpoint.ownership.channel != ownership.channel
        {
            return Err(RollingLossError::CheckpointMismatch { field: "ownership" });
        }
        Ok(checkpoint.state.clone())
    }

    fn applicability_kind(&self) -> RollingLossApplicabilityKind {
        match self {
            Self::CoulombContour(_) => RollingLossApplicabilityKind::LeineStyleCoulombContour,
            Self::ViscousContour(_) => RollingLossApplicabilityKind::ViscousContour,
            Self::HystereticRollingMoment(_) => {
                RollingLossApplicabilityKind::FinitePatchHystereticRollingMoment
            }
        }
    }
}

/// Rolling and spin kinematics over one caller-owned accepted interval.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RollingKinematics {
    contour_speed_mps: f64,
    rolling_rate_rad_s: f64,
    spin_rate_rad_s: f64,
    temperature_kelvin: f64,
    excitation_frequency_hz: f64,
    interval_s: f64,
}

impl RollingKinematics {
    /// Constructs finite kinematics. The source-card validity check happens at
    /// law evaluation; this constructor does not infer frequency from a rate.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        contour_speed_mps: f64,
        rolling_rate_rad_s: f64,
        spin_rate_rad_s: f64,
        temperature_kelvin: f64,
        excitation_frequency_hz: f64,
        interval_s: f64,
    ) -> Result<Self, RollingLossError> {
        finite(contour_speed_mps, "contour_speed_mps")?;
        finite(rolling_rate_rad_s, "rolling_rate_rad_s")?;
        finite(spin_rate_rad_s, "spin_rate_rad_s")?;
        positive(temperature_kelvin, "temperature_kelvin")?;
        nonnegative(excitation_frequency_hz, "excitation_frequency_hz")?;
        nonnegative(interval_s, "interval_s")?;
        Ok(Self {
            contour_speed_mps,
            rolling_rate_rad_s,
            spin_rate_rad_s,
            temperature_kelvin,
            excitation_frequency_hz,
            interval_s,
        })
    }

    /// Material contact-point contour speed [m/s].
    #[must_use]
    pub const fn contour_speed_mps(&self) -> f64 {
        self.contour_speed_mps
    }

    /// Rolling rate [rad/s].
    #[must_use]
    pub const fn rolling_rate_rad_s(&self) -> f64 {
        self.rolling_rate_rad_s
    }

    /// Spin rate [rad/s], retained but not consumed by these non-microslip rungs.
    #[must_use]
    pub const fn spin_rate_rad_s(&self) -> f64 {
        self.spin_rate_rad_s
    }

    /// Absolute temperature [K].
    #[must_use]
    pub const fn temperature_kelvin(&self) -> f64 {
        self.temperature_kelvin
    }

    /// Caller-declared excitation frequency [Hz].
    #[must_use]
    pub const fn excitation_frequency_hz(&self) -> f64 {
        self.excitation_frequency_hz
    }

    /// Accepted interval duration [s].
    #[must_use]
    pub const fn interval_s(&self) -> f64 {
        self.interval_s
    }
}

/// Distinct deformation-loss work channels.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RollingLossChannel {
    /// Work conjugate to material contact-point contour speed.
    ContourDeformation,
    /// Work conjugate to rolling rate for the hysteretic candidate.
    RollingHysteresis,
}

/// Generalized-work ownership for this module.
///
/// A caller must use a dedicated channel and may conservatively check an
/// interval against partial-slip ownership before composing both responses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RollingWorkOwnership {
    patch_id: String,
    interval_id: String,
    generalized_coordinate_id: String,
    channel: RollingLossChannel,
}

impl RollingWorkOwnership {
    /// Creates an explicit patch/interval/coordinate work key.
    pub fn new(
        patch_id: impl Into<String>,
        interval_id: impl Into<String>,
        generalized_coordinate_id: impl Into<String>,
        channel: RollingLossChannel,
    ) -> Result<Self, RollingLossError> {
        let value = Self {
            patch_id: patch_id.into(),
            interval_id: interval_id.into(),
            generalized_coordinate_id: generalized_coordinate_id.into(),
            channel,
        };
        value.validate_for(&value.patch_id, channel)?;
        Ok(value)
    }

    /// Patch identity owning the work interval.
    #[must_use]
    pub fn patch_id(&self) -> &str {
        &self.patch_id
    }

    /// Caller interval identity.
    #[must_use]
    pub fn interval_id(&self) -> &str {
        &self.interval_id
    }

    /// Dedicated generalized-coordinate identity.
    #[must_use]
    pub fn generalized_coordinate_id(&self) -> &str {
        &self.generalized_coordinate_id
    }

    /// Dedicated deformation-loss channel.
    #[must_use]
    pub const fn channel(&self) -> RollingLossChannel {
        self.channel
    }

    /// Conservatively refuses reusing a partial-slip patch/interval pair.
    ///
    /// `GeneralizedWorkOwnership` intentionally keeps its coordinate labels
    /// private, so this leaf refuses the complete shared interval rather than
    /// guessing whether a hidden coordinate is disjoint.
    pub fn require_disjoint_from_partial_slip(
        &self,
        partial_slip: &GeneralizedWorkOwnership,
    ) -> Result<(), RollingLossError> {
        if self.patch_id == partial_slip.patch_id()
            && self.interval_id == partial_slip.interval_id()
        {
            return Err(RollingLossError::WorkOwnershipOverlap);
        }
        Ok(())
    }

    fn validate_for(
        &self,
        patch_id: &str,
        required_channel: RollingLossChannel,
    ) -> Result<(), RollingLossError> {
        nonblank(&self.patch_id, "work_patch_id")?;
        nonblank(&self.interval_id, "work_interval_id")?;
        nonblank(&self.generalized_coordinate_id, "generalized_coordinate_id")?;
        if self.patch_id != patch_id {
            return Err(RollingLossError::WorkOwnershipMismatch { field: "patch_id" });
        }
        if self.channel != required_channel {
            return Err(RollingLossError::WorkOwnershipMismatch { field: "channel" });
        }
        Ok(())
    }
}

/// Internal state for deterministic, caller-committed replay.
#[derive(Debug, Clone, PartialEq)]
pub struct RollingLossState {
    cumulative_heat_j: f64,
    accepted_steps: u64,
}

impl RollingLossState {
    /// Empty reversible state: these rungs retain no reversible energy.
    #[must_use]
    pub const fn zero() -> Self {
        Self {
            cumulative_heat_j: 0.0,
            accepted_steps: 0,
        }
    }

    /// Total heat accumulated across accepted intervals [J].
    #[must_use]
    pub const fn cumulative_heat_j(&self) -> f64 {
        self.cumulative_heat_j
    }

    /// Number of accepted intervals.
    #[must_use]
    pub const fn accepted_steps(&self) -> u64 {
        self.accepted_steps
    }

    fn validate(&self) -> Result<(), RollingLossError> {
        nonnegative(self.cumulative_heat_j, "cumulative_heat_j")
    }
}

impl Default for RollingLossState {
    fn default() -> Self {
        Self::zero()
    }
}

/// Candidate branch selected by a distinct loss law.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RollingLossStateKind {
    /// Zero load, rate, or card factor produces no loss.
    Quiescent,
    /// A contour-force rung dissipates through contact-point contour speed.
    ContourDeformation,
    /// A finite-patch rolling moment dissipates through rolling rate.
    RollingHysteresis,
}

/// Scalar generalized wrench components on the body.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RollingLossWrench {
    /// Force conjugate to contour speed [N].
    pub contour_force_n: f64,
    /// Moment conjugate to rolling rate [N m].
    pub rolling_moment_nm: f64,
    /// Moment conjugate to spin [N m]; zero for these separately scoped rungs.
    pub spin_moment_nm: f64,
}

/// Passive heat channels [J].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RollingLossDissipation {
    /// Heat from a contour-deformation force.
    pub contour_deformation_heat_j: f64,
    /// Heat from a finite-patch rolling-moment candidate.
    pub rolling_hysteresis_heat_j: f64,
    /// Sum of declared heat channels.
    pub total_heat_j: f64,
}

/// Exact scalar work accounting for an accepted interval.
#[derive(Debug, Clone, PartialEq)]
pub struct RollingGeneralizedWork {
    /// Dedicated generalized-work owner.
    pub ownership: RollingWorkOwnership,
    /// Positive passive work entering the interface [J].
    pub work_into_interface_j: f64,
    /// Wrench power on the body; passive loss is non-positive [W].
    pub endpoint_body_power_w: f64,
    /// `-work_into_interface_j / interval_s`, or zero for zero duration [W].
    pub reconstructed_body_power_w: f64,
}

/// Model scope retained with every candidate response.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RollingLossApplicabilityKind {
    /// Leine-style Coulomb contour-force form.
    LeineStyleCoulombContour,
    /// Linear viscous contour alternative.
    ViscousContour,
    /// Finite-patch hysteretic rolling-moment candidate.
    FinitePatchHystereticRollingMoment,
}

/// Explicit no-certificate uncertainty statement.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RollingLossUncertainty {
    /// The model form has no numerical error certificate; an optional curvature
    /// approximation authority is retained without being upgraded.
    ModelFormNoNumericalCertificate {
        /// Authority for an equivalent-radius approximation, if one was used.
        curvature_approximation_authority: Option<InputAuthority>,
    },
}

/// Lineage-bound rolling-loss replay candidate.
///
/// The checkpoint retains the full accepted-step patch and work receipt for
/// auditability, while restoration deliberately admits current normal geometry
/// and a fresh interval key with stable contact and work-channel lineage.
#[derive(Debug, Clone, PartialEq)]
pub struct RollingLossCheckpoint {
    patch: RollingPatchReceipt,
    interface: InterfaceSystemRef,
    law: RollingLossLaw,
    ownership: RollingWorkOwnership,
    state: RollingLossState,
}

impl RollingLossCheckpoint {
    fn validate(&self) -> Result<(), RollingLossError> {
        self.patch.validate()?;
        validate_interface(&self.interface)?;
        self.law.validate()?;
        self.ownership
            .validate_for(self.patch.patch_id(), self.law.required_channel())?;
        self.state.validate()
    }
}

/// One candidate rolling-loss response.
#[derive(Debug, Clone, PartialEq)]
pub struct RollingLossStep {
    /// Explicit candidate branch.
    pub state: RollingLossStateKind,
    /// Generalized wrench components.
    pub wrench: RollingLossWrench,
    /// Always zero for this irreversible-loss slice.
    pub stored_energy_j: f64,
    /// Always zero for this irreversible-loss slice.
    pub storage_change_j: f64,
    /// Passive heat decomposition.
    pub dissipation: RollingLossDissipation,
    /// Independent work/heat accounting with ownership keys.
    pub generalized_work: RollingGeneralizedWork,
    /// Caller-committed candidate internal state.
    pub next_state: RollingLossState,
    /// Identity-bound candidate restart receipt.
    pub checkpoint: RollingLossCheckpoint,
    /// Distinct model scope.
    pub applicability: RollingLossApplicabilityKind,
    /// Explicit no-certificate uncertainty status.
    pub uncertainty: RollingLossUncertainty,
    /// Caller normal-patch authority ceiling.
    pub patch_authority: InputAuthority,
    /// Ordered interface authority ceiling.
    pub interface_authority: InputAuthority,
    /// Immutable law-card authority ceiling.
    pub card_authority: InputAuthority,
}

fn validate_interface(interface: &InterfaceSystemRef) -> Result<(), RollingLossError> {
    if interface.medium() != InterfaceMedium::Dry {
        return Err(RollingLossError::NotDryInterface {
            medium: interface.medium(),
        });
    }
    nonblank(interface.ordered_system_id(), "ordered_system_id")?;
    nonblank(interface.history_id(), "history_id")?;
    nonblank(interface.provenance().source_id(), "interface_source_id")
}

fn require_range(
    range: ApplicabilityRange,
    value: f64,
    field: &'static str,
) -> Result<(), RollingLossError> {
    nonnegative(value, field)?;
    if value < range.minimum() || value > range.maximum() {
        return Err(RollingLossError::OutsideApplicability {
            field,
            value,
            minimum: range.minimum(),
            maximum: range.maximum(),
        });
    }
    Ok(())
}

fn nonblank(value: &str, field: &'static str) -> Result<(), RollingLossError> {
    if value.trim().is_empty() {
        return Err(RollingLossError::MissingIdentity { field });
    }
    Ok(())
}

fn finite(value: f64, field: &'static str) -> Result<(), RollingLossError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(RollingLossError::InvalidInput { field })
    }
}

fn nonnegative(value: f64, field: &'static str) -> Result<(), RollingLossError> {
    finite(value, field)?;
    if value < 0.0 {
        return Err(RollingLossError::InvalidInput { field });
    }
    Ok(())
}

fn positive(value: f64, field: &'static str) -> Result<(), RollingLossError> {
    finite(value, field)?;
    if value <= 0.0 {
        return Err(RollingLossError::InvalidInput { field });
    }
    Ok(())
}

fn checked_add(a: f64, b: f64, field: &'static str) -> Result<f64, RollingLossError> {
    let value = a + b;
    finite(value, field)?;
    Ok(value)
}

fn checked_mul(a: f64, b: f64, field: &'static str) -> Result<f64, RollingLossError> {
    let value = a * b;
    finite(value, field)?;
    Ok(value)
}
