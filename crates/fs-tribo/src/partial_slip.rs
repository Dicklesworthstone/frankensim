//! Solver-independent finite-patch tangential partial-slip rung.
//!
//! This module intentionally consumes a neutral [`NormalPatchView`] rather
//! than a contact-solver type.  It therefore cannot determine normal contact,
//! evolve patch geometry, or admit a material card.  The caller supplies all
//! such data with named identities and an authority ceiling.
//!
//! The law is a bounded Cattaneo--Mindlin-style *return-mapping rung*, not a
//! resolved pressure/traction field.  It represents a reversible tangential
//! spring core plus a plastic microslip remainder.  `PartialSlip` means that
//! remainder is nonzero in this lumped law; it is **not** a claimed resolved
//! slipping-area fraction.  Rolling-deformation loss is deliberately absent:
//! its energy channel is always zero and must be composed separately.

use core::fmt;

use fs_math::det;

const EPSILON: f64 = 128.0 * f64::EPSILON;

/// Stable identity for this bounded constitutive rung.
pub const PARTIAL_SLIP_MODEL_ID: &str = "fs-tribo/cattaneo-mindlin-style-return-map-v1";

/// Caller-declared ceiling on a normal-patch receipt.
///
/// This is deliberately not an admission or validation claim.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NormalPatchAuthority {
    /// A caller named a source but this module did not independently admit it.
    CallerDeclared,
    /// A synthetic fixture used only for numerical tests or exercises.
    SyntheticFixture,
    /// An explicitly estimated patch input.
    Estimated,
}

/// Neutral, solver-independent normal-patch input in SI units.
///
/// `pressure_second_moment_m2` is the normalized pressure second moment
/// `integral(r^2 p dA) / normal_load`; its square root supplies the torsional
/// lever radius in this lumped rung.  It is not a resolved pressure field.
#[derive(Debug, Clone, PartialEq)]
pub struct NormalPatchView {
    patch_id: String,
    normal_model_id: String,
    source_id: String,
    authority: NormalPatchAuthority,
    normal_load_n: f64,
    semi_axis_longitudinal_m: f64,
    semi_axis_lateral_m: f64,
    pressure_second_moment_m2: f64,
}

impl NormalPatchView {
    /// Constructs a fully named finite, positively loaded normal-patch view.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        patch_id: impl Into<String>,
        normal_model_id: impl Into<String>,
        source_id: impl Into<String>,
        authority: NormalPatchAuthority,
        normal_load_n: f64,
        semi_axis_longitudinal_m: f64,
        semi_axis_lateral_m: f64,
        pressure_second_moment_m2: f64,
    ) -> Result<Self, PartialSlipError> {
        let value = Self {
            patch_id: patch_id.into(),
            normal_model_id: normal_model_id.into(),
            source_id: source_id.into(),
            authority,
            normal_load_n,
            semi_axis_longitudinal_m,
            semi_axis_lateral_m,
            pressure_second_moment_m2,
        };
        value.validate()?;
        Ok(value)
    }

    /// Stable caller patch identity.
    #[must_use]
    pub fn patch_id(&self) -> &str {
        &self.patch_id
    }

    /// Identity of the normal-response model that produced this view.
    #[must_use]
    pub fn normal_model_id(&self) -> &str {
        &self.normal_model_id
    }

    /// Caller source identity for the normal-patch data.
    #[must_use]
    pub fn source_id(&self) -> &str {
        &self.source_id
    }

    /// Authority ceiling retained from the caller.
    #[must_use]
    pub const fn authority(&self) -> NormalPatchAuthority {
        self.authority
    }

    /// Compressive normal load in N.
    #[must_use]
    pub const fn normal_load_n(&self) -> f64 {
        self.normal_load_n
    }

    /// Longitudinal patch semi-axis in m.
    #[must_use]
    pub const fn semi_axis_longitudinal_m(&self) -> f64 {
        self.semi_axis_longitudinal_m
    }

    /// Lateral patch semi-axis in m.
    #[must_use]
    pub const fn semi_axis_lateral_m(&self) -> f64 {
        self.semi_axis_lateral_m
    }

    /// Normalized pressure second moment in m².
    #[must_use]
    pub const fn pressure_second_moment_m2(&self) -> f64 {
        self.pressure_second_moment_m2
    }

    fn validate(&self) -> Result<(), PartialSlipError> {
        nonblank(&self.patch_id, "patch_id")?;
        nonblank(&self.normal_model_id, "normal_model_id")?;
        nonblank(&self.source_id, "normal_patch_source_id")?;
        positive(self.normal_load_n, "normal_load_n")?;
        positive(self.semi_axis_longitudinal_m, "semi_axis_longitudinal_m")?;
        positive(self.semi_axis_lateral_m, "semi_axis_lateral_m")?;
        positive(self.pressure_second_moment_m2, "pressure_second_moment_m2")?;
        Ok(())
    }
}

/// Ordered dry interface and history identity consumed by this law.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PartialSlipInterface {
    ordered_interface_id: String,
    history_id: String,
    source_id: String,
    authority: NormalPatchAuthority,
}

impl PartialSlipInterface {
    /// Creates an ordered, named dry-interface identity.
    pub fn new(
        ordered_interface_id: impl Into<String>,
        history_id: impl Into<String>,
        source_id: impl Into<String>,
        authority: NormalPatchAuthority,
    ) -> Result<Self, PartialSlipError> {
        let value = Self {
            ordered_interface_id: ordered_interface_id.into(),
            history_id: history_id.into(),
            source_id: source_id.into(),
            authority,
        };
        value.validate()?;
        Ok(value)
    }

    /// Ordered interface identity; reversing surfaces requires a new value.
    #[must_use]
    pub fn ordered_interface_id(&self) -> &str {
        &self.ordered_interface_id
    }

    /// History owner identity.
    #[must_use]
    pub fn history_id(&self) -> &str {
        &self.history_id
    }

    /// Caller source identity.
    #[must_use]
    pub fn source_id(&self) -> &str {
        &self.source_id
    }

    /// Caller-declared authority ceiling.
    #[must_use]
    pub const fn authority(&self) -> NormalPatchAuthority {
        self.authority
    }

    fn validate(&self) -> Result<(), PartialSlipError> {
        nonblank(&self.ordered_interface_id, "ordered_interface_id")?;
        nonblank(&self.history_id, "history_id")?;
        nonblank(&self.source_id, "interface_source_id")?;
        Ok(())
    }
}

/// A right-handed orthonormal frame for a contact tangent plane.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TangentFrame {
    longitudinal: [f64; 3],
    lateral: [f64; 3],
    normal: [f64; 3],
}

impl TangentFrame {
    /// Normalizes `normal` and projects/normalizes `longitudinal_hint` into its plane.
    pub fn new(normal: [f64; 3], longitudinal_hint: [f64; 3]) -> Result<Self, PartialSlipError> {
        finite_vec3(normal, "normal")?;
        finite_vec3(longitudinal_hint, "longitudinal_hint")?;
        let normal = normalize(normal, "normal")?;
        let normal_part = dot(normal, longitudinal_hint, "longitudinal_normal_component")?;
        let longitudinal = sub(
            longitudinal_hint,
            scale3(normal, normal_part, "longitudinal_projection")?,
            "longitudinal_projection",
        )?;
        let longitudinal = normalize(longitudinal, "longitudinal_hint")?;
        let lateral = normalize(cross(normal, longitudinal, "lateral")?, "lateral")?;
        Ok(Self {
            longitudinal,
            lateral,
            normal,
        })
    }

    /// Unit longitudinal tangent direction.
    #[must_use]
    pub const fn longitudinal(&self) -> [f64; 3] {
        self.longitudinal
    }

    /// Unit lateral tangent direction.
    #[must_use]
    pub const fn lateral(&self) -> [f64; 3] {
        self.lateral
    }

    /// Unit contact normal, completing a right-handed frame.
    #[must_use]
    pub const fn normal(&self) -> [f64; 3] {
        self.normal
    }

    /// Rotates tangent axes by `angle_rad`, preserving the physical normal.
    pub fn rotated(self, angle_rad: f64) -> Result<Self, PartialSlipError> {
        finite(angle_rad, "tangent_rotation_rad")?;
        let c = det::cos(angle_rad);
        let s = det::sin(angle_rad);
        finite(c, "tangent_rotation_cos")?;
        finite(s, "tangent_rotation_sin")?;
        let longitudinal = add(
            scale3(self.longitudinal, c, "rotated_longitudinal")?,
            scale3(self.lateral, s, "rotated_longitudinal")?,
            "rotated_longitudinal",
        )?;
        Self::new(self.normal, longitudinal)
    }

    fn to_world(self, tangent: [f64; 2]) -> Result<[f64; 3], PartialSlipError> {
        finite_vec2(tangent, "tangent_components")?;
        add(
            scale3(self.longitudinal, tangent[0], "tangent_world")?,
            scale3(self.lateral, tangent[1], "tangent_world")?,
            "tangent_world",
        )
    }
}

/// Declared longitudinal/lateral creepage and torsional spin for one interval.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PartialSlipKinematics {
    /// Dimensionless longitudinal/lateral creepage in the supplied tangent frame.
    pub creepage: [f64; 2],
    /// Non-negative rolling-speed magnitude in m/s used to convert creepage to slip velocity.
    pub rolling_speed_mps: f64,
    /// Relative spin about the contact normal in rad/s.
    pub torsional_spin_rad_per_s: f64,
    /// Positive interval duration in s.
    pub dt_s: f64,
}

impl PartialSlipKinematics {
    fn validate(self) -> Result<(), PartialSlipError> {
        finite_vec2(self.creepage, "creepage")?;
        nonnegative(self.rolling_speed_mps, "rolling_speed_mps")?;
        finite(self.torsional_spin_rad_per_s, "torsional_spin_rad_per_s")?;
        positive(self.dt_s, "dt_s")?;
        if self.rolling_speed_mps == 0.0 && (self.creepage[0] != 0.0 || self.creepage[1] != 0.0) {
            return Err(PartialSlipError::InvalidInput {
                field: "nonzero creepage requires rolling_speed_mps",
            });
        }
        Ok(())
    }

    fn slip_velocity_tangent_mps(self) -> Result<[f64; 2], PartialSlipError> {
        scale2(self.creepage, self.rolling_speed_mps, "creepage_velocity")
    }
}

/// Caller-declared coefficient set for the lumped finite-patch law.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PartialSlipParameters {
    /// Static Coulomb coefficient, dimensionless.
    pub static_mu: f64,
    /// Kinetic Coulomb coefficient, dimensionless and no greater than `static_mu`.
    pub kinetic_mu: f64,
    /// Reversible translational stiffness in N/m.
    pub tangential_stiffness_n_per_m: f64,
    /// Reversible torsional stiffness in N m/rad.
    pub torsional_stiffness_nm_per_rad: f64,
    /// Dimensionless torsional capacity factor in `(0, 1]` multiplying `mu N sqrt(moment)`.
    pub torsional_capacity_factor: f64,
    /// Dimensionless onset of this lumped partial-slip rung in `(0, 1)`.
    pub partial_slip_onset_fraction: f64,
    /// Dimensionless retained hardening fraction in `[0, 1)` during partial slip.
    pub partial_slip_hardening_fraction: f64,
}

impl PartialSlipParameters {
    fn validate(self) -> Result<(), PartialSlipError> {
        positive(self.static_mu, "static_mu")?;
        positive(self.kinetic_mu, "kinetic_mu")?;
        if self.kinetic_mu > self.static_mu {
            return Err(PartialSlipError::InvalidInput {
                field: "kinetic_mu > static_mu",
            });
        }
        positive(
            self.tangential_stiffness_n_per_m,
            "tangential_stiffness_n_per_m",
        )?;
        positive(
            self.torsional_stiffness_nm_per_rad,
            "torsional_stiffness_nm_per_rad",
        )?;
        positive(self.torsional_capacity_factor, "torsional_capacity_factor")?;
        if self.torsional_capacity_factor > 1.0 {
            return Err(PartialSlipError::InvalidInput {
                field: "torsional_capacity_factor > 1",
            });
        }
        positive(
            self.partial_slip_onset_fraction,
            "partial_slip_onset_fraction",
        )?;
        if self.partial_slip_onset_fraction >= 1.0 {
            return Err(PartialSlipError::InvalidInput {
                field: "partial_slip_onset_fraction >= 1",
            });
        }
        nonnegative(
            self.partial_slip_hardening_fraction,
            "partial_slip_hardening_fraction",
        )?;
        if self.partial_slip_hardening_fraction >= 1.0 {
            return Err(PartialSlipError::InvalidInput {
                field: "partial_slip_hardening_fraction >= 1",
            });
        }
        Ok(())
    }
}

/// Explicit model/rung identity and coefficients.
#[derive(Debug, Clone, PartialEq)]
pub struct PartialSlipLaw {
    model_id: String,
    source_id: String,
    parameters: PartialSlipParameters,
}

impl PartialSlipLaw {
    /// Creates a named generic Cattaneo--Mindlin-style return-mapping law.
    pub fn new(
        model_id: impl Into<String>,
        source_id: impl Into<String>,
        parameters: PartialSlipParameters,
    ) -> Result<Self, PartialSlipError> {
        let value = Self {
            model_id: model_id.into(),
            source_id: source_id.into(),
            parameters,
        };
        value.validate()?;
        Ok(value)
    }

    /// Stable caller model identity.
    #[must_use]
    pub fn model_id(&self) -> &str {
        &self.model_id
    }

    /// Caller source identity for coefficients.
    #[must_use]
    pub fn source_id(&self) -> &str {
        &self.source_id
    }

    /// Checked caller coefficients.
    #[must_use]
    pub const fn parameters(&self) -> PartialSlipParameters {
        self.parameters
    }

    /// Advances one interval without mutating the supplied state on refusal.
    #[allow(clippy::too_many_arguments)]
    pub fn advance(
        &self,
        patch: &NormalPatchView,
        interface: &PartialSlipInterface,
        frame: TangentFrame,
        kinematics: PartialSlipKinematics,
        ownership: &GeneralizedWorkOwnership,
        state: &PartialSlipState,
    ) -> Result<PartialSlipStep, PartialSlipError> {
        self.validate()?;
        patch.validate()?;
        interface.validate()?;
        kinematics.validate()?;
        ownership.validate_for_patch(patch.patch_id())?;
        state.validate()?;

        let velocity = kinematics.slip_velocity_tangent_mps()?;
        let delta_displacement = scale2(velocity, kinematics.dt_s, "tangential_displacement")?;
        let delta_twist = checked_mul(
            kinematics.torsional_spin_rad_per_s,
            kinematics.dt_s,
            "torsional_displacement",
        )?;
        let trial_displacement = add2(
            state.elastic_displacement_m,
            delta_displacement,
            "trial_displacement",
        )?;
        let trial_twist = checked_add(state.elastic_twist_rad, delta_twist, "trial_twist")?;
        let trial_force = scale2(
            trial_displacement,
            -self.parameters.tangential_stiffness_n_per_m,
            "trial_force",
        )?;
        let trial_torque = checked_mul(
            -self.parameters.torsional_stiffness_nm_per_rad,
            trial_twist,
            "trial_torque",
        )?;
        let capacity = CoulombCapacity::from_patch(patch, self.parameters)?;
        let trial_utilization = capacity.utilization(trial_force, trial_torque)?;

        let (state_kind, retained_utilization, slip_fraction) = if trial_utilization == 0.0 {
            (PartialSlipStateKind::Sticking, 0.0, 0.0)
        } else if trial_utilization <= self.parameters.partial_slip_onset_fraction {
            (PartialSlipStateKind::Sticking, trial_utilization, 0.0)
        } else if trial_utilization < 1.0 {
            let retained = self.parameters.partial_slip_onset_fraction
                + (trial_utilization - self.parameters.partial_slip_onset_fraction)
                    * self.parameters.partial_slip_hardening_fraction;
            let fraction = (trial_utilization - self.parameters.partial_slip_onset_fraction)
                / (1.0 - self.parameters.partial_slip_onset_fraction);
            (PartialSlipStateKind::PartialSlip, retained, fraction)
        } else {
            (
                PartialSlipStateKind::GrossSlide,
                self.parameters.kinetic_mu / self.parameters.static_mu,
                1.0,
            )
        };
        finite(retained_utilization, "retained_utilization")?;
        let retained_scale = if trial_utilization == 0.0 {
            1.0
        } else {
            retained_utilization / trial_utilization
        };
        finite(retained_scale, "retained_scale")?;
        let next_displacement = scale2(trial_displacement, retained_scale, "next_displacement")?;
        let next_twist = checked_mul(trial_twist, retained_scale, "next_twist")?;
        let force_tangent = scale2(
            next_displacement,
            -self.parameters.tangential_stiffness_n_per_m,
            "tangential_force",
        )?;
        let torsional_moment_nm = checked_mul(
            -self.parameters.torsional_stiffness_nm_per_rad,
            next_twist,
            "torsional_moment",
        )?;
        let output_utilization = capacity.utilization(force_tangent, torsional_moment_nm)?;
        if output_utilization > 1.0 + EPSILON {
            return Err(PartialSlipError::InvalidDerived {
                field: "output_coulomb_utilization",
            });
        }

        let plastic_displacement = sub2(
            trial_displacement,
            next_displacement,
            "plastic_displacement",
        )?;
        let plastic_twist = checked_add(-next_twist, trial_twist, "plastic_twist")?;
        let dissipation_j = nonnegative_sum(
            &[
                -dot2(
                    force_tangent,
                    plastic_displacement,
                    "tangential_plastic_work",
                )?,
                -checked_mul(torsional_moment_nm, plastic_twist, "torsional_plastic_work")?,
            ],
            "microslip_dissipation_j",
        )?;
        let old_storage = stored_energy(
            self.parameters,
            state.elastic_displacement_m,
            state.elastic_twist_rad,
        )?;
        let new_storage = stored_energy(self.parameters, next_displacement, next_twist)?;
        let storage_change_j = checked_add(new_storage, -old_storage, "storage_change_j")?;
        let work_into_interface_j = checked_add(storage_change_j, dissipation_j, "work_closure_j")?;
        let reconstructed_body_power_w = -work_into_interface_j / kinematics.dt_s;
        finite(reconstructed_body_power_w, "reconstructed_body_power_w")?;
        let wrench = TangentialWrench {
            force_n: frame.to_world(force_tangent)?,
            torque_nm: scale3(frame.normal, torsional_moment_nm, "torsional_wrench")?,
        };
        let endpoint_body_power_w = checked_add(
            dot(
                wrench.force_n,
                frame.to_world(velocity)?,
                "endpoint_force_power",
            )?,
            checked_mul(
                torsional_moment_nm,
                kinematics.torsional_spin_rad_per_s,
                "endpoint_torsion_power",
            )?,
            "endpoint_body_power",
        )?;
        let next_state = PartialSlipState {
            elastic_displacement_m: next_displacement,
            elastic_twist_rad: next_twist,
            accepted_steps: state.accepted_steps.checked_add(1).ok_or(
                PartialSlipError::InvalidDerived {
                    field: "accepted_steps overflow",
                },
            )?,
        };
        next_state.validate()?;
        let checkpoint = PartialSlipCheckpoint::new(
            patch.clone(),
            interface.clone(),
            self.clone(),
            next_state.clone(),
        )?;
        Ok(PartialSlipStep {
            state: state_kind,
            slip_partition: SlipPartition {
                microslip_fraction: slip_fraction,
                rolling_deformation_fraction: 0.0,
            },
            wrench,
            tangent_force_n: force_tangent,
            torsional_moment_nm,
            capacity,
            stored_energy_j: new_storage,
            storage_change_j,
            dissipation: DissipationAndHeat {
                tangential_and_torsional_microslip_j: dissipation_j,
                rolling_deformation_loss_j: 0.0,
                heat_j: dissipation_j,
            },
            generalized_work: GeneralizedWork {
                ownership: ownership.clone(),
                work_into_interface_j,
                reconstructed_body_power_w,
                endpoint_body_power_w,
            },
            next_state,
            checkpoint,
            applicability: PartialSlipApplicability::CattaneoMindlinStyleLumpedReturnMap,
            normal_patch_authority: patch.authority,
            interface_authority: interface.authority,
        })
    }

    /// Restores state only when the constitutive law and contact lineage match.
    ///
    /// In particular, a matching model name alone is insufficient: source identity,
    /// caller authority, every coefficient, and ordered interface provenance bind
    /// a checkpoint.  The normal patch's identity, normal-model identity, source,
    /// and authority bind its lineage. Its load, semi-axes, and pressure moment are
    /// *current step data*: [`Self::advance`] recomputes the Coulomb capacity from
    /// them, so a normal solver may evolve a continuous contact without erasing
    /// the tangential elastic state.
    ///
    /// The checkpoint still retains the complete normal-patch receipt from the
    /// preceding accepted step for auditability. This method does not claim that
    /// equal lineage values establish external physical validity; it merely
    /// refuses an ambiguous constitutive or provenance restart.
    pub fn restore_checkpoint(
        &self,
        patch: &NormalPatchView,
        interface: &PartialSlipInterface,
        checkpoint: &PartialSlipCheckpoint,
    ) -> Result<PartialSlipState, PartialSlipError> {
        self.validate()?;
        patch.validate()?;
        interface.validate()?;
        checkpoint.validate()?;
        if checkpoint.law.model_id != self.model_id {
            return Err(PartialSlipError::CheckpointIdentityMismatch {
                field: "law_model_id",
            });
        }
        if checkpoint.law.source_id != self.source_id {
            return Err(PartialSlipError::CheckpointIdentityMismatch {
                field: "law_source_id",
            });
        }
        if checkpoint.law.parameters != self.parameters {
            return Err(PartialSlipError::CheckpointIdentityMismatch {
                field: "law_parameters",
            });
        }
        if checkpoint.patch.patch_id != patch.patch_id {
            return Err(PartialSlipError::CheckpointIdentityMismatch { field: "patch_id" });
        }
        if checkpoint.patch.normal_model_id != patch.normal_model_id {
            return Err(PartialSlipError::CheckpointIdentityMismatch {
                field: "normal_model_id",
            });
        }
        if checkpoint.patch.source_id != patch.source_id {
            return Err(PartialSlipError::CheckpointIdentityMismatch {
                field: "normal_patch_source_id",
            });
        }
        if checkpoint.patch.authority != patch.authority {
            return Err(PartialSlipError::CheckpointIdentityMismatch {
                field: "normal_patch_authority",
            });
        }
        if checkpoint.interface.ordered_interface_id != interface.ordered_interface_id {
            return Err(PartialSlipError::CheckpointIdentityMismatch {
                field: "ordered_interface_id",
            });
        }
        if checkpoint.interface.history_id != interface.history_id {
            return Err(PartialSlipError::CheckpointIdentityMismatch {
                field: "history_id",
            });
        }
        if checkpoint.interface.source_id != interface.source_id {
            return Err(PartialSlipError::CheckpointIdentityMismatch {
                field: "interface_source_id",
            });
        }
        if checkpoint.interface.authority != interface.authority {
            return Err(PartialSlipError::CheckpointIdentityMismatch {
                field: "interface_authority",
            });
        }
        Ok(checkpoint.state.clone())
    }

    fn validate(&self) -> Result<(), PartialSlipError> {
        nonblank(&self.model_id, "law_model_id")?;
        nonblank(&self.source_id, "law_source_id")?;
        self.parameters.validate()
    }
}

/// Explicit contact state emitted by the tangential law.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PartialSlipStateKind {
    /// Entire lumped spring increment is reversible.
    Sticking,
    /// A nonzero lumped microslip remainder was dissipated while capacity remains below static limit.
    PartialSlip,
    /// The combined force/torsion resultant is projected to kinetic Coulomb capacity.
    GrossSlide,
}

/// Reversible internal state suitable for checkpointing and deterministic replay.
#[derive(Debug, Clone, PartialEq)]
pub struct PartialSlipState {
    elastic_displacement_m: [f64; 2],
    elastic_twist_rad: f64,
    accepted_steps: u64,
}

impl PartialSlipState {
    /// Empty history with no stored tangential energy.
    #[must_use]
    pub const fn zero() -> Self {
        Self {
            elastic_displacement_m: [0.0; 2],
            elastic_twist_rad: 0.0,
            accepted_steps: 0,
        }
    }

    /// Reconstructs validated reversible state from a checkpoint decoder.
    pub fn from_checkpoint(
        elastic_displacement_m: [f64; 2],
        elastic_twist_rad: f64,
        accepted_steps: u64,
    ) -> Result<Self, PartialSlipError> {
        let value = Self {
            elastic_displacement_m,
            elastic_twist_rad,
            accepted_steps,
        };
        value.validate()?;
        Ok(value)
    }

    /// Elastic tangent displacement in m.
    #[must_use]
    pub const fn elastic_displacement_m(&self) -> [f64; 2] {
        self.elastic_displacement_m
    }

    /// Elastic torsional displacement in rad.
    #[must_use]
    pub const fn elastic_twist_rad(&self) -> f64 {
        self.elastic_twist_rad
    }

    /// Number of accepted interval updates.
    #[must_use]
    pub const fn accepted_steps(&self) -> u64 {
        self.accepted_steps
    }

    fn validate(&self) -> Result<(), PartialSlipError> {
        finite_vec2(self.elastic_displacement_m, "elastic_displacement_m")?;
        finite(self.elastic_twist_rad, "elastic_twist_rad")
    }
}

/// Lineage-bound snapshot of reversible partial-slip history.
///
/// The snapshot retains the complete inputs from its accepted step instead of a
/// lossy hash, so a decoder can retain a transparent receipt. Numeric normal
/// load and geometry are deliberately not restore identity: they are resolved
/// again by the current normal step. It does not validate the physical truth of
/// caller-declared inputs.
#[derive(Debug, Clone, PartialEq)]
pub struct PartialSlipCheckpoint {
    patch: NormalPatchView,
    interface: PartialSlipInterface,
    law: PartialSlipLaw,
    state: PartialSlipState,
}

impl PartialSlipCheckpoint {
    /// Reconstructs a replay-bound checkpoint from validated decoded inputs.
    pub fn new(
        patch: NormalPatchView,
        interface: PartialSlipInterface,
        law: PartialSlipLaw,
        state: PartialSlipState,
    ) -> Result<Self, PartialSlipError> {
        let value = Self {
            patch,
            interface,
            law,
            state,
        };
        value.validate()?;
        Ok(value)
    }

    /// Patch identity bound into this checkpoint.
    #[must_use]
    pub fn patch_id(&self) -> &str {
        self.patch.patch_id()
    }

    /// History identity bound into this checkpoint.
    #[must_use]
    pub fn history_id(&self) -> &str {
        self.interface.history_id()
    }

    /// Law identity bound into this checkpoint.
    #[must_use]
    pub fn law_model_id(&self) -> &str {
        self.law.model_id()
    }

    /// Complete neutral normal-patch receipt bound into this checkpoint.
    #[must_use]
    pub fn normal_patch(&self) -> &NormalPatchView {
        &self.patch
    }

    /// Complete ordered-interface receipt bound into this checkpoint.
    #[must_use]
    pub fn interface(&self) -> &PartialSlipInterface {
        &self.interface
    }

    /// Complete law/source/coefficient receipt bound into this checkpoint.
    #[must_use]
    pub fn law(&self) -> &PartialSlipLaw {
        &self.law
    }

    /// Stored reversible state.
    #[must_use]
    pub fn state(&self) -> &PartialSlipState {
        &self.state
    }

    fn validate(&self) -> Result<(), PartialSlipError> {
        self.patch.validate()?;
        self.interface.validate()?;
        self.law.validate()?;
        self.state.validate()
    }
}

/// Tangential force and normal-axis torque acting on the body with the declared relative motion.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TangentialWrench {
    /// World-frame force in N.
    pub force_n: [f64; 3],
    /// World-frame torque in N m.
    pub torque_nm: [f64; 3],
}

/// Coupled Coulomb force and torsional capacity in SI units.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CoulombCapacity {
    /// Static tangential-force limit in N.
    pub static_force_n: f64,
    /// Kinetic tangential-force limit in N.
    pub kinetic_force_n: f64,
    /// Static torsional-moment limit in N m.
    pub static_torque_nm: f64,
    /// Kinetic torsional-moment limit in N m.
    pub kinetic_torque_nm: f64,
}

impl CoulombCapacity {
    fn from_patch(
        patch: &NormalPatchView,
        parameters: PartialSlipParameters,
    ) -> Result<Self, PartialSlipError> {
        let radius = patch.pressure_second_moment_m2.sqrt();
        positive(radius, "torsional_pressure_radius_m")?;
        let static_force_n = checked_mul(
            parameters.static_mu,
            patch.normal_load_n,
            "static_force_capacity",
        )?;
        let kinetic_force_n = checked_mul(
            parameters.kinetic_mu,
            patch.normal_load_n,
            "kinetic_force_capacity",
        )?;
        let static_torque_nm = checked_mul(
            checked_mul(static_force_n, radius, "static_torque_capacity")?,
            parameters.torsional_capacity_factor,
            "static_torque_capacity",
        )?;
        let kinetic_torque_nm = checked_mul(
            checked_mul(kinetic_force_n, radius, "kinetic_torque_capacity")?,
            parameters.torsional_capacity_factor,
            "kinetic_torque_capacity",
        )?;
        positive(static_force_n, "static_force_capacity")?;
        positive(static_torque_nm, "static_torque_capacity")?;
        positive(kinetic_force_n, "kinetic_force_capacity")?;
        positive(kinetic_torque_nm, "kinetic_torque_capacity")?;
        Ok(Self {
            static_force_n,
            kinetic_force_n,
            static_torque_nm,
            kinetic_torque_nm,
        })
    }

    fn utilization(self, force: [f64; 2], torque_nm: f64) -> Result<f64, PartialSlipError> {
        let force_ratio = norm2(force, "force_utilization")? / self.static_force_n;
        let torque_ratio = torque_nm.abs() / self.static_torque_nm;
        let value = det::hypot(force_ratio, torque_ratio);
        finite(value, "coulomb_utilization")?;
        Ok(value)
    }
}

/// Generalized-coordinate ownership labels used to prevent work double counting.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct GeneralizedWorkOwnership {
    patch_id: String,
    interval_id: String,
    longitudinal_coordinate_id: String,
    lateral_coordinate_id: String,
    torsional_coordinate_id: String,
}

impl GeneralizedWorkOwnership {
    /// Creates ownership labels for exactly one patch and interval.
    pub fn new(
        patch_id: impl Into<String>,
        interval_id: impl Into<String>,
        longitudinal_coordinate_id: impl Into<String>,
        lateral_coordinate_id: impl Into<String>,
        torsional_coordinate_id: impl Into<String>,
    ) -> Result<Self, PartialSlipError> {
        let value = Self {
            patch_id: patch_id.into(),
            interval_id: interval_id.into(),
            longitudinal_coordinate_id: longitudinal_coordinate_id.into(),
            lateral_coordinate_id: lateral_coordinate_id.into(),
            torsional_coordinate_id: torsional_coordinate_id.into(),
        };
        value.validate_for_patch(&value.patch_id)?;
        Ok(value)
    }

    /// Retargets this owner to another caller interval while reusing its
    /// existing string storage.
    ///
    /// All prospective fields are validated before any field is changed, so a
    /// refusal leaves the previous owner intact. This is equivalent to
    /// constructing a fresh owner from the same five strings, but retains each
    /// field's capacity; no field allocates when its replacement bytes fit.
    pub fn retarget(
        &mut self,
        patch_id: &str,
        interval_id: &str,
        longitudinal_coordinate_id: &str,
        lateral_coordinate_id: &str,
        torsional_coordinate_id: &str,
    ) -> Result<(), PartialSlipError> {
        nonblank(patch_id, "work_patch_id")?;
        nonblank(interval_id, "work_interval_id")?;
        nonblank(longitudinal_coordinate_id, "longitudinal_coordinate_id")?;
        nonblank(lateral_coordinate_id, "lateral_coordinate_id")?;
        nonblank(torsional_coordinate_id, "torsional_coordinate_id")?;

        self.patch_id.clear();
        self.patch_id.push_str(patch_id);
        self.interval_id.clear();
        self.interval_id.push_str(interval_id);
        self.longitudinal_coordinate_id.clear();
        self.longitudinal_coordinate_id
            .push_str(longitudinal_coordinate_id);
        self.lateral_coordinate_id.clear();
        self.lateral_coordinate_id.push_str(lateral_coordinate_id);
        self.torsional_coordinate_id.clear();
        self.torsional_coordinate_id
            .push_str(torsional_coordinate_id);
        Ok(())
    }

    /// Patch identity which owns this work interval.
    #[must_use]
    pub fn patch_id(&self) -> &str {
        &self.patch_id
    }

    /// Caller interval identity.
    #[must_use]
    pub fn interval_id(&self) -> &str {
        &self.interval_id
    }

    /// Longitudinal generalized-coordinate identity.
    #[must_use]
    pub fn longitudinal_coordinate_id(&self) -> &str {
        &self.longitudinal_coordinate_id
    }

    /// Lateral generalized-coordinate identity.
    #[must_use]
    pub fn lateral_coordinate_id(&self) -> &str {
        &self.lateral_coordinate_id
    }

    /// Torsional generalized-coordinate identity.
    #[must_use]
    pub fn torsional_coordinate_id(&self) -> &str {
        &self.torsional_coordinate_id
    }

    fn validate_for_patch(&self, patch_id: &str) -> Result<(), PartialSlipError> {
        nonblank(&self.patch_id, "work_patch_id")?;
        nonblank(&self.interval_id, "work_interval_id")?;
        nonblank(
            &self.longitudinal_coordinate_id,
            "longitudinal_coordinate_id",
        )?;
        nonblank(&self.lateral_coordinate_id, "lateral_coordinate_id")?;
        nonblank(&self.torsional_coordinate_id, "torsional_coordinate_id")?;
        if self.patch_id != patch_id {
            return Err(PartialSlipError::WorkOwnershipMismatch);
        }
        Ok(())
    }
}

/// Partition reported by the lumped law.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SlipPartition {
    /// Model partition only; not a resolved physical slipping-area fraction.
    pub microslip_fraction: f64,
    /// Always zero: rolling-deformation loss belongs to another named mechanism.
    pub rolling_deformation_fraction: f64,
}

/// Dissipated energy channels in J.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DissipationAndHeat {
    /// Combined tangential-creepage and torsional-microslip heat in J.
    pub tangential_and_torsional_microslip_j: f64,
    /// Always zero; this law never models rolling-deformation loss.
    pub rolling_deformation_loss_j: f64,
    /// Heat generated by this dry microslip rung in J.
    pub heat_j: f64,
}

/// Work accounting whose sign is defined as positive into the interface.
#[derive(Debug, Clone, PartialEq)]
pub struct GeneralizedWork {
    /// The unique patch/interval/generalized-coordinate owner.
    pub ownership: GeneralizedWorkOwnership,
    /// Exact discrete work closure: storage change plus heat, in J.
    pub work_into_interface_j: f64,
    /// `-work_into_interface_j / dt_s`, W; negative means body energy is dissipated/stored.
    pub reconstructed_body_power_w: f64,
    /// Endpoint wrench power `F·v + M·ω`, W, retained separately from discrete closure.
    pub endpoint_body_power_w: f64,
}

/// Applicability identity for this generic rung.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PartialSlipApplicability {
    /// Dry elastic finite patch; caller supplied load, moment, coefficients, and frame.
    CattaneoMindlinStyleLumpedReturnMap,
}

/// One accepted constitutive interval.
#[derive(Debug, Clone, PartialEq)]
pub struct PartialSlipStep {
    /// Explicit branch selected for the interval.
    pub state: PartialSlipStateKind,
    /// Lumped microslip partition; it never reports rolling-deformation loss.
    pub slip_partition: SlipPartition,
    /// World-frame tangential wrench on the body.
    pub wrench: TangentialWrench,
    /// Tangential force coordinates in the declared tangent frame, N.
    pub tangent_force_n: [f64; 2],
    /// Normal-axis torsional moment, N m.
    pub torsional_moment_nm: f64,
    /// Coupled Coulomb capacity used for this interval.
    pub capacity: CoulombCapacity,
    /// Reversible stored energy after the update, J.
    pub stored_energy_j: f64,
    /// Change in reversible stored energy, J.
    pub storage_change_j: f64,
    /// Passive heat/dissipation ledger.
    pub dissipation: DissipationAndHeat,
    /// Exact discrete energy closure and ownership keys.
    pub generalized_work: GeneralizedWork,
    /// Candidate internal state; the caller commits it only after solver acceptance.
    pub next_state: PartialSlipState,
    /// Identity-bound restart candidate.
    pub checkpoint: PartialSlipCheckpoint,
    /// Model scope retained with the response.
    pub applicability: PartialSlipApplicability,
    /// Normal-patch authority ceiling.
    pub normal_patch_authority: NormalPatchAuthority,
    /// Interface authority ceiling.
    pub interface_authority: NormalPatchAuthority,
}

/// Total refusal surface for the partial-slip module.
#[derive(Debug, Clone, PartialEq)]
pub enum PartialSlipError {
    /// Required identity was blank.
    MissingIdentity { field: &'static str },
    /// Input was non-finite, non-positive where required, or outside this rung's domain.
    InvalidInput { field: &'static str },
    /// A finite input produced a non-finite candidate.
    InvalidDerived { field: &'static str },
    /// Work keys named a different patch than the normal-patch view.
    WorkOwnershipMismatch,
    /// A checkpoint identity does not match current caller-owned identities.
    CheckpointIdentityMismatch { field: &'static str },
}

impl fmt::Display for PartialSlipError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingIdentity { field } => write!(f, "nonblank identity required: {field}"),
            Self::InvalidInput { field } => write!(f, "invalid partial-slip input: {field}"),
            Self::InvalidDerived { field } => {
                write!(f, "invalid partial-slip derived value: {field}")
            }
            Self::WorkOwnershipMismatch => write!(f, "generalized work ownership patch mismatch"),
            Self::CheckpointIdentityMismatch { field } => {
                write!(f, "partial-slip checkpoint identity mismatch: {field}")
            }
        }
    }
}

impl std::error::Error for PartialSlipError {}

fn stored_energy(
    parameters: PartialSlipParameters,
    displacement: [f64; 2],
    twist: f64,
) -> Result<f64, PartialSlipError> {
    let translation = checked_mul(
        0.5 * parameters.tangential_stiffness_n_per_m,
        dot2(displacement, displacement, "stored_translation")?,
        "stored_translation",
    )?;
    let torsion = checked_mul(
        0.5 * parameters.torsional_stiffness_nm_per_rad,
        checked_mul(twist, twist, "stored_torsion")?,
        "stored_torsion",
    )?;
    nonnegative_sum(&[translation, torsion], "stored_energy")
}

fn nonblank(value: &str, field: &'static str) -> Result<(), PartialSlipError> {
    if value.trim().is_empty() {
        Err(PartialSlipError::MissingIdentity { field })
    } else {
        Ok(())
    }
}

fn finite(value: f64, field: &'static str) -> Result<(), PartialSlipError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(PartialSlipError::InvalidDerived { field })
    }
}

fn positive(value: f64, field: &'static str) -> Result<(), PartialSlipError> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(PartialSlipError::InvalidInput { field })
    }
}

fn nonnegative(value: f64, field: &'static str) -> Result<(), PartialSlipError> {
    if value.is_finite() && value >= 0.0 {
        Ok(())
    } else {
        Err(PartialSlipError::InvalidInput { field })
    }
}

fn finite_vec2(value: [f64; 2], field: &'static str) -> Result<(), PartialSlipError> {
    finite(value[0], field)?;
    finite(value[1], field)
}

fn finite_vec3(value: [f64; 3], field: &'static str) -> Result<(), PartialSlipError> {
    finite(value[0], field)?;
    finite(value[1], field)?;
    finite(value[2], field)
}

fn checked_add(left: f64, right: f64, field: &'static str) -> Result<f64, PartialSlipError> {
    let value = left + right;
    finite(value, field)?;
    Ok(value)
}

fn checked_mul(left: f64, right: f64, field: &'static str) -> Result<f64, PartialSlipError> {
    let value = left * right;
    finite(value, field)?;
    Ok(value)
}

fn scale2(value: [f64; 2], scale: f64, field: &'static str) -> Result<[f64; 2], PartialSlipError> {
    Ok([
        checked_mul(value[0], scale, field)?,
        checked_mul(value[1], scale, field)?,
    ])
}

fn add2(
    left: [f64; 2],
    right: [f64; 2],
    field: &'static str,
) -> Result<[f64; 2], PartialSlipError> {
    Ok([
        checked_add(left[0], right[0], field)?,
        checked_add(left[1], right[1], field)?,
    ])
}

fn sub2(
    left: [f64; 2],
    right: [f64; 2],
    field: &'static str,
) -> Result<[f64; 2], PartialSlipError> {
    Ok([
        checked_add(left[0], -right[0], field)?,
        checked_add(left[1], -right[1], field)?,
    ])
}

fn dot2(left: [f64; 2], right: [f64; 2], field: &'static str) -> Result<f64, PartialSlipError> {
    checked_add(
        checked_mul(left[0], right[0], field)?,
        checked_mul(left[1], right[1], field)?,
        field,
    )
}

fn norm2(value: [f64; 2], field: &'static str) -> Result<f64, PartialSlipError> {
    let result = det::hypot(value[0], value[1]);
    finite(result, field)?;
    Ok(result)
}

fn scale3(value: [f64; 3], scale: f64, field: &'static str) -> Result<[f64; 3], PartialSlipError> {
    Ok([
        checked_mul(value[0], scale, field)?,
        checked_mul(value[1], scale, field)?,
        checked_mul(value[2], scale, field)?,
    ])
}

fn add(left: [f64; 3], right: [f64; 3], field: &'static str) -> Result<[f64; 3], PartialSlipError> {
    Ok([
        checked_add(left[0], right[0], field)?,
        checked_add(left[1], right[1], field)?,
        checked_add(left[2], right[2], field)?,
    ])
}

fn sub(left: [f64; 3], right: [f64; 3], field: &'static str) -> Result<[f64; 3], PartialSlipError> {
    Ok([
        checked_add(left[0], -right[0], field)?,
        checked_add(left[1], -right[1], field)?,
        checked_add(left[2], -right[2], field)?,
    ])
}

fn dot(left: [f64; 3], right: [f64; 3], field: &'static str) -> Result<f64, PartialSlipError> {
    checked_add(
        checked_add(
            checked_mul(left[0], right[0], field)?,
            checked_mul(left[1], right[1], field)?,
            field,
        )?,
        checked_mul(left[2], right[2], field)?,
        field,
    )
}

fn cross(
    left: [f64; 3],
    right: [f64; 3],
    field: &'static str,
) -> Result<[f64; 3], PartialSlipError> {
    Ok([
        checked_add(
            checked_mul(left[1], right[2], field)?,
            -checked_mul(left[2], right[1], field)?,
            field,
        )?,
        checked_add(
            checked_mul(left[2], right[0], field)?,
            -checked_mul(left[0], right[2], field)?,
            field,
        )?,
        checked_add(
            checked_mul(left[0], right[1], field)?,
            -checked_mul(left[1], right[0], field)?,
            field,
        )?,
    ])
}

fn normalize(value: [f64; 3], field: &'static str) -> Result<[f64; 3], PartialSlipError> {
    let magnitude = det::hypot(det::hypot(value[0], value[1]), value[2]);
    positive(magnitude, field)?;
    scale3(value, 1.0 / magnitude, field)
}

fn nonnegative_sum(values: &[f64], field: &'static str) -> Result<f64, PartialSlipError> {
    let mut total = 0.0;
    for &value in values {
        nonnegative(value, field)?;
        total = checked_add(total, value, field)?;
    }
    Ok(total)
}

#[cfg(test)]
mod tests {
    use super::*;

    const TOLERANCE: f64 = 2.0e-12;

    fn patch() -> NormalPatchView {
        NormalPatchView::new(
            "patch-a",
            "normal-card-v1",
            "synthetic/normal-patch",
            NormalPatchAuthority::SyntheticFixture,
            100.0,
            0.02,
            0.01,
            1.0e-4,
        )
        .expect("synthetic normal patch")
    }

    fn interface() -> PartialSlipInterface {
        PartialSlipInterface::new(
            "body-a->support-b",
            "history-a",
            "synthetic/dry-interface",
            NormalPatchAuthority::SyntheticFixture,
        )
        .expect("synthetic interface")
    }

    fn law() -> PartialSlipLaw {
        PartialSlipLaw::new(
            PARTIAL_SLIP_MODEL_ID,
            "synthetic/partial-slip-coefficients",
            PartialSlipParameters {
                static_mu: 0.8,
                kinetic_mu: 0.4,
                tangential_stiffness_n_per_m: 10_000.0,
                torsional_stiffness_nm_per_rad: 100.0,
                torsional_capacity_factor: 0.5,
                partial_slip_onset_fraction: 0.5,
                partial_slip_hardening_fraction: 0.4,
            },
        )
        .expect("synthetic law")
    }

    fn frame() -> TangentFrame {
        TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).expect("frame")
    }

    fn ownership() -> GeneralizedWorkOwnership {
        GeneralizedWorkOwnership::new("patch-a", "interval-0", "qx", "qy", "qspin")
            .expect("ownership")
    }

    fn kinematics(creepage: [f64; 2], spin: f64) -> PartialSlipKinematics {
        PartialSlipKinematics {
            creepage,
            rolling_speed_mps: 1.0,
            torsional_spin_rad_per_s: spin,
            dt_s: 0.001,
        }
    }

    fn step(
        frame: TangentFrame,
        state: &PartialSlipState,
        creepage: [f64; 2],
        spin: f64,
    ) -> PartialSlipStep {
        law()
            .advance(
                &patch(),
                &interface(),
                frame,
                kinematics(creepage, spin),
                &ownership(),
                state,
            )
            .expect("accepted partial-slip step")
    }

    fn close(left: f64, right: f64) {
        assert!(
            (left - right).abs() <= TOLERANCE * left.abs().max(right.abs()).max(1.0),
            "left={left:.17e}, right={right:.17e}"
        );
    }

    fn close_vec(left: [f64; 3], right: [f64; 3]) {
        for index in 0..3 {
            close(left[index], right[index]);
        }
    }

    fn world(frame: TangentFrame, tangent: [f64; 2]) -> [f64; 3] {
        [
            frame.longitudinal()[0] * tangent[0] + frame.lateral()[0] * tangent[1],
            frame.longitudinal()[1] * tangent[0] + frame.lateral()[1] * tangent[1],
            frame.longitudinal()[2] * tangent[0] + frame.lateral()[2] * tangent[1],
        ]
    }

    /// G0: zero relative rate selects stick and carries no invented heat.
    #[test]
    fn zero_rate_is_sticking_with_exact_zero_channels() {
        let response = step(frame(), &PartialSlipState::zero(), [0.0; 2], 0.0);
        assert_eq!(response.state, PartialSlipStateKind::Sticking);
        assert_eq!(response.tangent_force_n, [0.0; 2]);
        assert_eq!(response.torsional_moment_nm, 0.0);
        assert_eq!(response.dissipation.heat_j, 0.0);
        assert_eq!(response.dissipation.rolling_deformation_loss_j, 0.0);
        assert_eq!(response.generalized_work.work_into_interface_j, 0.0);
    }

    /// G1 limiting cases: stick, the declared partial rung, and kinetic gross slide.
    #[test]
    fn force_capacity_branches_are_explicit_and_bounded() {
        let state = PartialSlipState::zero();
        let sticking = step(frame(), &state, [1.0, 0.0], 0.0);
        assert_eq!(sticking.state, PartialSlipStateKind::Sticking);
        close(sticking.tangent_force_n[0], -10.0);

        let partial = step(frame(), &state, [6.0, 0.0], 0.0);
        assert_eq!(partial.state, PartialSlipStateKind::PartialSlip);
        assert!(partial.slip_partition.microslip_fraction > 0.0);
        assert!(partial.slip_partition.microslip_fraction < 1.0);
        assert!(partial.dissipation.heat_j > 0.0);
        assert!(partial.tangent_force_n[0].abs() < partial.capacity.static_force_n);

        let gross = step(frame(), &state, [20.0, 0.0], 0.0);
        assert_eq!(gross.state, PartialSlipStateKind::GrossSlide);
        close(
            gross.tangent_force_n[0].abs(),
            gross.capacity.kinetic_force_n,
        );
        assert_eq!(gross.slip_partition.microslip_fraction, 1.0);
        assert_eq!(gross.dissipation.rolling_deformation_loss_j, 0.0);
    }

    /// G1 pure torsion uses the finite-patch pressure moment without inventing translation.
    #[test]
    fn pure_torsion_is_coulomb_limited_and_separate_from_rolling_loss() {
        let response = step(frame(), &PartialSlipState::zero(), [0.0; 2], 10.0);
        assert_eq!(response.state, PartialSlipStateKind::GrossSlide);
        assert_eq!(response.tangent_force_n, [0.0; 2]);
        close(
            response.torsional_moment_nm.abs(),
            response.capacity.kinetic_torque_nm,
        );
        close_vec(response.wrench.force_n, [0.0; 3]);
        close_vec(
            response.wrench.torque_nm,
            [0.0, 0.0, response.torsional_moment_nm],
        );
        assert_eq!(response.dissipation.rolling_deformation_loss_j, 0.0);
    }

    /// G0 exact discrete accounting: stored-energy change plus heat owns all interval work.
    #[test]
    fn work_heat_and_storage_close_exactly() {
        let response = step(frame(), &PartialSlipState::zero(), [6.0, 2.0], 2.0);
        assert!(response.dissipation.heat_j >= 0.0);
        close(
            response.generalized_work.work_into_interface_j,
            response.storage_change_j + response.dissipation.heat_j,
        );
        close(
            response.generalized_work.reconstructed_body_power_w,
            -response.generalized_work.work_into_interface_j / 0.001,
        );
        assert_eq!(
            response.dissipation.heat_j,
            response.dissipation.tangential_and_torsional_microslip_j
        );
    }

    /// G3 reversal is deterministic and changes the gross-slide force direction.
    #[test]
    fn reversal_preserves_history_then_reverses_resistance() {
        let first = step(frame(), &PartialSlipState::zero(), [20.0, 0.0], 0.0);
        let reversed = step(frame(), &first.next_state, [-20.0, 0.0], 0.0);
        assert_eq!(first.state, PartialSlipStateKind::GrossSlide);
        assert_eq!(reversed.state, PartialSlipStateKind::GrossSlide);
        assert!(first.tangent_force_n[0] < 0.0);
        assert!(reversed.tangent_force_n[0] > 0.0);
        assert!(reversed.dissipation.heat_j >= 0.0);
        assert_eq!(reversed.next_state.accepted_steps(), 2);
    }

    /// G3: passive SO(2) tangent-frame re-expression preserves physical wrench and power.
    #[test]
    fn tangent_frame_rotation_preserves_wrench_and_reconstructed_power() {
        let base = frame();
        let angle = 0.37;
        let rotated = base.rotated(angle).expect("rotation");
        let creepage = [6.0, 2.0];
        // det-ok: the INDEPENDENT oracle arm stays on the platform
        // sequence so the comparison is cross-implementation, not
        // self-referential (tolerance absorbs the ULP-class gap).
        let c = angle.cos(); // det-ok: independent oracle arm (platform)
        let s = angle.sin(); // det-ok: independent oracle arm (platform)
        let rotated_creepage = [
            c * creepage[0] + s * creepage[1],
            -s * creepage[0] + c * creepage[1],
        ];
        let base_response = step(base, &PartialSlipState::zero(), creepage, 1.5);
        let rotated_response = step(rotated, &PartialSlipState::zero(), rotated_creepage, 1.5);
        close_vec(
            base_response.wrench.force_n,
            rotated_response.wrench.force_n,
        );
        close_vec(
            base_response.wrench.torque_nm,
            rotated_response.wrench.torque_nm,
        );
        close(
            base_response.generalized_work.reconstructed_body_power_w,
            rotated_response.generalized_work.reconstructed_body_power_w,
        );
        let base_velocity = world(base, creepage);
        let rotated_velocity = world(rotated, rotated_creepage);
        close_vec(base_velocity, rotated_velocity);
        close(
            dot(base_response.wrench.force_n, base_velocity, "test_power").expect("power")
                + base_response.torsional_moment_nm * 1.5,
            base_response.generalized_work.endpoint_body_power_w,
        );
    }

    /// G3 restart: exact caller inputs restore the deterministic reversible history.
    #[test]
    fn checkpoint_restart_is_exact_with_full_receipt() {
        let response = step(frame(), &PartialSlipState::zero(), [6.0, 0.0], 1.0);
        let restored = law()
            .restore_checkpoint(&patch(), &interface(), &response.checkpoint)
            .expect("matching checkpoint");
        assert_eq!(restored, response.next_state);

        let wrong_patch = NormalPatchView::new(
            "other-patch",
            "normal-card-v1",
            "synthetic/normal-patch",
            NormalPatchAuthority::SyntheticFixture,
            100.0,
            0.02,
            0.01,
            1.0e-4,
        )
        .expect("valid mutation");
        let wrong = PartialSlipCheckpoint::new(wrong_patch, interface(), law(), restored)
            .expect("syntactically valid checkpoint mutation");
        assert_eq!(
            law().restore_checkpoint(&patch(), &interface(), &wrong),
            Err(PartialSlipError::CheckpointIdentityMismatch { field: "patch_id" })
        );
    }

    /// G3: a normal solve may evolve its resolved load and patch geometry while
    /// retaining the same constitutive lineage and tangential history.
    #[test]
    fn checkpoint_admits_evolved_normal_step_data_but_refuses_lineage_change() {
        let first = step(frame(), &PartialSlipState::zero(), [6.0, 0.0], 1.0);
        let evolved_patch = NormalPatchView::new(
            "patch-a",
            "normal-card-v1",
            "synthetic/normal-patch",
            NormalPatchAuthority::SyntheticFixture,
            125.0,
            0.024,
            0.012,
            1.5e-4,
        )
        .expect("same-lineage normal update");
        let restored = law()
            .restore_checkpoint(&evolved_patch, &interface(), &first.checkpoint)
            .expect("evolved normal step data must restore tangential history");
        assert_eq!(restored, first.next_state);

        let second = law()
            .advance(
                &evolved_patch,
                &interface(),
                frame(),
                kinematics([6.0, 0.0], 1.0),
                &ownership(),
                &restored,
            )
            .expect("evolved normal step must advance");
        assert_eq!(second.next_state.accepted_steps(), 2);
        close(second.capacity.static_force_n, 100.0);
        assert_eq!(second.checkpoint.normal_patch(), &evolved_patch);

        let changed_patch_lineage = NormalPatchView::new(
            "other-patch",
            "normal-card-v1",
            "synthetic/normal-patch",
            NormalPatchAuthority::SyntheticFixture,
            125.0,
            0.024,
            0.012,
            1.5e-4,
        )
        .expect("valid lineage mutation");
        assert_eq!(
            law().restore_checkpoint(&changed_patch_lineage, &interface(), &first.checkpoint,),
            Err(PartialSlipError::CheckpointIdentityMismatch { field: "patch_id" })
        );
    }

    /// G3: replay refuses every law/source/authority receipt mutation, not only display IDs.
    #[test]
    fn checkpoint_rejects_coefficients_and_admitted_provenance_mutations() {
        let response = step(frame(), &PartialSlipState::zero(), [6.0, 0.0], 1.0);
        let state = response.next_state;
        let expect_refusal = |current_law: PartialSlipLaw,
                              current_patch: NormalPatchView,
                              current_interface: PartialSlipInterface,
                              field| {
            assert_eq!(
                current_law.restore_checkpoint(
                    &current_patch,
                    &current_interface,
                    &response.checkpoint,
                ),
                Err(PartialSlipError::CheckpointIdentityMismatch { field })
            );
        };

        expect_refusal(
            PartialSlipLaw::new(
                PARTIAL_SLIP_MODEL_ID,
                "different-law-source",
                law().parameters(),
            )
            .expect("valid source mutation"),
            patch(),
            interface(),
            "law_source_id",
        );
        for parameters in [
            PartialSlipParameters {
                static_mu: 0.81,
                ..law().parameters()
            },
            PartialSlipParameters {
                kinetic_mu: 0.41,
                ..law().parameters()
            },
            PartialSlipParameters {
                tangential_stiffness_n_per_m: 10_001.0,
                ..law().parameters()
            },
            PartialSlipParameters {
                torsional_stiffness_nm_per_rad: 101.0,
                ..law().parameters()
            },
            PartialSlipParameters {
                torsional_capacity_factor: 0.6,
                ..law().parameters()
            },
            PartialSlipParameters {
                partial_slip_onset_fraction: 0.6,
                ..law().parameters()
            },
            PartialSlipParameters {
                partial_slip_hardening_fraction: 0.5,
                ..law().parameters()
            },
        ] {
            expect_refusal(
                PartialSlipLaw::new(PARTIAL_SLIP_MODEL_ID, law().source_id(), parameters)
                    .expect("valid coefficient mutation"),
                patch(),
                interface(),
                "law_parameters",
            );
        }

        for (current_patch, field) in [
            (
                NormalPatchView::new(
                    "patch-a",
                    "other-normal-model",
                    "synthetic/normal-patch",
                    NormalPatchAuthority::SyntheticFixture,
                    100.0,
                    0.02,
                    0.01,
                    1.0e-4,
                )
                .expect("valid model mutation"),
                "normal_model_id",
            ),
            (
                NormalPatchView::new(
                    "patch-a",
                    "normal-card-v1",
                    "other-normal-source",
                    NormalPatchAuthority::SyntheticFixture,
                    100.0,
                    0.02,
                    0.01,
                    1.0e-4,
                )
                .expect("valid source mutation"),
                "normal_patch_source_id",
            ),
            (
                NormalPatchView::new(
                    "patch-a",
                    "normal-card-v1",
                    "synthetic/normal-patch",
                    NormalPatchAuthority::Estimated,
                    100.0,
                    0.02,
                    0.01,
                    1.0e-4,
                )
                .expect("valid authority mutation"),
                "normal_patch_authority",
            ),
        ] {
            expect_refusal(law(), current_patch, interface(), field);
        }

        for (current_interface, field) in [
            (
                PartialSlipInterface::new(
                    "other-body->support-b",
                    "history-a",
                    "synthetic/dry-interface",
                    NormalPatchAuthority::SyntheticFixture,
                )
                .expect("valid ordered-interface mutation"),
                "ordered_interface_id",
            ),
            (
                PartialSlipInterface::new(
                    "body-a->support-b",
                    "history-a",
                    "other-interface-source",
                    NormalPatchAuthority::SyntheticFixture,
                )
                .expect("valid interface source mutation"),
                "interface_source_id",
            ),
            (
                PartialSlipInterface::new(
                    "body-a->support-b",
                    "history-a",
                    "synthetic/dry-interface",
                    NormalPatchAuthority::Estimated,
                )
                .expect("valid interface authority mutation"),
                "interface_authority",
            ),
        ] {
            expect_refusal(law(), patch(), current_interface, field);
        }

        let replacement = PartialSlipCheckpoint::new(patch(), interface(), law(), state)
            .expect("replacement checkpoint");
        assert_eq!(replacement, response.checkpoint);
    }

    /// G0 hostile admission: unknown/invalid normal, rate, ownership, and maximum-plus-one state refuse.
    #[test]
    fn invalid_inputs_refuse_without_a_candidate_state() {
        assert!(
            NormalPatchView::new(
                "patch",
                "normal",
                "source",
                NormalPatchAuthority::CallerDeclared,
                0.0,
                0.01,
                0.01,
                1.0e-4
            )
            .is_err()
        );
        assert!(
            PartialSlipLaw::new(
                "law",
                "source",
                PartialSlipParameters {
                    kinetic_mu: 0.9,
                    ..law().parameters()
                }
            )
            .is_err()
        );
        let bad_rate = PartialSlipKinematics {
            creepage: [1.0, 0.0],
            rolling_speed_mps: 0.0,
            torsional_spin_rad_per_s: 0.0,
            dt_s: 0.001,
        };
        assert!(
            law()
                .advance(
                    &patch(),
                    &interface(),
                    frame(),
                    bad_rate,
                    &ownership(),
                    &PartialSlipState::zero()
                )
                .is_err()
        );
        let other_patch = GeneralizedWorkOwnership::new("other", "i", "x", "y", "z").expect("keys");
        assert_eq!(
            law().advance(
                &patch(),
                &interface(),
                frame(),
                kinematics([0.0; 2], 0.0),
                &other_patch,
                &PartialSlipState::zero(),
            ),
            Err(PartialSlipError::WorkOwnershipMismatch)
        );
        let exhausted = PartialSlipState::from_checkpoint([0.0; 2], 0.0, u64::MAX).expect("state");
        assert_eq!(
            law().advance(
                &patch(),
                &interface(),
                frame(),
                kinematics([0.0; 2], 0.0),
                &ownership(),
                &exhausted,
            ),
            Err(PartialSlipError::InvalidDerived {
                field: "accepted_steps overflow"
            })
        );
    }
}
