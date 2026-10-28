//! Solver-independent finite-patch tangential partial-slip rung.
//!
//! This module intentionally consumes a neutral [`NormalPatchView`] rather
//! than a contact-solver type.  It therefore cannot determine normal contact,
//! evolve patch geometry, or admit a material card.  The caller supplies all
//! such data with named identities and an authority ceiling.
//!
//! The law is a two-region Jenkins/Iwan *return-mapping rung*, not a resolved
//! Cattaneo--Mindlin pressure/traction field. Each parallel region retains its
//! own elastic shear and slip branch. `PartialSlip` means one region yields;
//! it is **not** a resolved slipping-area fraction. Backward-Euler work uses
//! the actual endpoint wrench. Its numerical dissipation is reported separately
//! from friction heat. Rolling-deformation loss belongs to another law.
//! Retained shear follows the contact normal's minimum rotation and the supplied
//! common material spin; tangent-gauge changes only re-express its coordinates.

// Explicit parameter-work continuation; ordinary restart remains coefficient-exact.
mod coefficients;

use core::fmt;

use fs_math::det;

const EPSILON: f64 = 128.0 * f64::EPSILON;

/// Stable identity for this bounded constitutive rung.
pub const PARTIAL_SLIP_MODEL_ID: &str = "fs-tribo/two-region-iwan-frame-transport-v3";

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
// Bead frankensim-30cz7: `_id` suffixes are deliberate identity vocabulary;
// renaming would break caller-facing field semantics for a lint preference.
#[allow(clippy::struct_field_names)]
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

    // SO(2) coordinates of the old frame after the minimum normal rotation,
    // followed by the caller's physical common spin about the current normal.
    // Gauge rotation is not material rotation. Construct one isometry for both
    // regions; projecting each shear separately would shorten stored springs.
    fn history_rotation_from(
        self,
        previous: Self,
        common_spin_angle_rad: f64,
    ) -> Result<[f64; 2], PartialSlipError> {
        if self == previous && common_spin_angle_rad == 0.0 {
            return Ok([1.0, 0.0]);
        }
        let mut rotation = if self == previous {
            [1.0, 0.0]
        } else {
            let normal_cosine =
                dot(previous.normal, self.normal, "contact_normal_rotation")?.clamp(-1.0, 1.0);
            let denominator = 1.0 + normal_cosine;
            if denominator <= EPSILON {
                // Antipodal normals do not determine a unique shortest rotation.
                // Require resolved intermediate geometry instead of selecting an
                // arbitrary physical direction from the caller's tangent gauge.
                return Err(PartialSlipError::InvalidInput {
                    field: "antipodal contact-normal transport",
                });
            }
            let projection = dot(
                previous.longitudinal,
                self.normal,
                "transported_longitudinal_projection",
            )? / denominator;
            let transported_longitudinal = sub(
                previous.longitudinal,
                scale3(
                    add(previous.normal, self.normal, "contact_normal_sum")?,
                    projection,
                    "transported_longitudinal",
                )?,
                "transported_longitudinal",
            )?;
            [
                dot(
                    self.longitudinal,
                    transported_longitudinal,
                    "history_rotation_cosine",
                )?,
                dot(
                    self.lateral,
                    transported_longitudinal,
                    "history_rotation_sine",
                )?,
            ]
        };
        if common_spin_angle_rad != 0.0 {
            let cosine = det::cos(common_spin_angle_rad);
            let sine = det::sin(common_spin_angle_rad);
            rotation = [
                cosine * rotation[0] - sine * rotation[1],
                sine * rotation[0] + cosine * rotation[1],
            ];
        }
        let magnitude = norm2(rotation, "history_rotation_norm")?;
        positive(magnitude, "history_rotation_norm")?;
        scale2(rotation, 1.0 / magnitude, "history_rotation")
    }
}

/// Declared relative motion and physical contact-frame spin for one interval.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PartialSlipKinematics {
    /// Dimensionless longitudinal/lateral creepage in the supplied tangent frame.
    pub creepage: [f64; 2],
    /// Non-negative rolling-speed magnitude in m/s used to convert creepage to slip velocity.
    pub rolling_speed_mps: f64,
    /// Relative spin about the contact normal in rad/s.
    pub torsional_spin_rad_per_s: f64,
    /// Common material-frame spin about the current normal in rad/s.
    /// For two rigid surfaces this is `(omega_a + omega_b) dot normal / 2`.
    /// It rotates retained shear without adding slip, twist, work, or heat.
    /// Supply it separately from the arbitrary tangent gauge and relative spin.
    pub contact_spin_rad_per_s: f64,
    /// Positive interval duration in s.
    pub dt_s: f64,
}

impl PartialSlipKinematics {
    fn validate(self) -> Result<(), PartialSlipError> {
        finite_vec2(self.creepage, "creepage")?;
        nonnegative(self.rolling_speed_mps, "rolling_speed_mps")?;
        finite(self.torsional_spin_rad_per_s, "torsional_spin_rad_per_s")?;
        finite(self.contact_spin_rad_per_s, "contact_spin_rad_per_s")?;
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
    /// First region's yield displacement divided by `mu_s N / k`, in `(0, 1)`.
    pub partial_slip_onset_fraction: f64,
    /// Stiffness fraction of the harder region, in `(0, 1)`. The two
    /// regions' capacities sum to the full Coulomb capacity; after the
    /// weaker region yields, this is the remaining elastic slope fraction.
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
        positive(
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
    /// Creates a named two-region Jenkins/Iwan return-mapping law.
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
    #[allow(clippy::too_many_lines)]
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
        let common_spin_angle = checked_mul(
            kinematics.contact_spin_rad_per_s,
            kinematics.dt_s,
            "contact_spin_angle",
        )?;
        let history_rotation = match state.reference_frame {
            Some(previous) => frame.history_rotation_from(previous, common_spin_angle)?,
            None => [1.0, 0.0],
        };
        let capacity = CoulombCapacity::from_patch(patch, self.parameters)?;
        let hard = self.parameters.partial_slip_hardening_fraction;
        let weights = [1.0 - hard, hard];
        // k_weak / k = 1-h, C_weak / C = a(1-h), so first yield
        // occurs at a C/k. The remaining region retains tangent h k.
        let weak_capacity = self.parameters.partial_slip_onset_fraction * weights[0];
        let capacity_fractions = [weak_capacity, 1.0 - weak_capacity];
        let increment = [delta_displacement[0], delta_displacement[1], delta_twist];
        let mut next_state = state.clone();
        next_state.reference_frame = Some(frame);
        let mut force_tangent = [0.0; 2];
        let mut torsional_moment_nm = 0.0;
        let mut old_storage = 0.0;
        let mut new_storage = 0.0;
        let mut dissipation_j = 0.0;
        let mut numerical_dissipation_j = 0.0;
        let mut slip_fraction = 0.0;
        let mut yielding_regions = 0;
        for region in 0..2 {
            let stored_old = [
                state.region_elastic_displacement_m[region][0],
                state.region_elastic_displacement_m[region][1],
                state.region_elastic_twist_rad[region],
            ];
            let old = [
                history_rotation[0] * stored_old[0] - history_rotation[1] * stored_old[1],
                history_rotation[1] * stored_old[0] + history_rotation[0] * stored_old[1],
                stored_old[2],
            ];
            let stiffness = [
                self.parameters.tangential_stiffness_n_per_m * weights[region],
                self.parameters.tangential_stiffness_n_per_m * weights[region],
                self.parameters.torsional_stiffness_nm_per_rad * weights[region],
            ];
            let limits = [
                capacity.static_force_n * capacity_fractions[region],
                capacity.static_force_n * capacity_fractions[region],
                capacity.static_torque_nm * capacity_fractions[region],
            ];
            let result = advance_region(
                old,
                increment,
                stiffness,
                limits,
                self.parameters.kinetic_mu / self.parameters.static_mu,
                state.sliding_regions[region],
            )?;
            next_state.region_elastic_displacement_m[region] =
                [result.elastic[0], result.elastic[1]];
            next_state.region_elastic_twist_rad[region] = result.elastic[2];
            next_state.sliding_regions[region] = result.kinetic_branch;
            for i in 0..2 {
                force_tangent[i] =
                    checked_add(force_tangent[i], -result.effort[i], "tangential_force")?;
            }
            torsional_moment_nm =
                checked_add(torsional_moment_nm, -result.effort[2], "torsional_moment")?;
            for i in 0..3 {
                let delta_elastic = result.elastic[i] - old[i];
                old_storage = checked_add(
                    old_storage,
                    checked_mul(
                        0.5 * stiffness[i] * stored_old[i],
                        stored_old[i],
                        "old_storage",
                    )?,
                    "old_storage",
                )?;
                new_storage = checked_add(
                    new_storage,
                    checked_mul(
                        0.5 * stiffness[i] * result.elastic[i],
                        result.elastic[i],
                        "new_storage",
                    )?,
                    "new_storage",
                )?;
                dissipation_j = checked_add(
                    dissipation_j,
                    checked_mul(result.effort[i], result.plastic[i], "plastic_work")?,
                    "plastic_work",
                )?;
                // The backward-Euler spring defect is not frictional heat.
                numerical_dissipation_j = checked_add(
                    numerical_dissipation_j,
                    checked_mul(
                        0.5 * stiffness[i] * delta_elastic,
                        delta_elastic,
                        "numerical_dissipation",
                    )?,
                    "numerical_dissipation",
                )?;
            }
            if result.yielded {
                yielding_regions += 1;
                slip_fraction += capacity_fractions[region];
            }
        }
        let state_kind = match yielding_regions {
            0 => PartialSlipStateKind::Sticking,
            1 => PartialSlipStateKind::PartialSlip,
            _ => PartialSlipStateKind::GrossSlide,
        };
        let output_utilization = capacity.utilization(force_tangent, torsional_moment_nm)?;
        if output_utilization > 1.0 + EPSILON {
            return Err(PartialSlipError::InvalidDerived {
                field: "output_coulomb_utilization",
            });
        }

        nonnegative(dissipation_j, "microslip_dissipation_j")?;
        nonnegative(numerical_dissipation_j, "numerical_dissipation_j")?;
        let storage_change_j = checked_add(new_storage, -old_storage, "storage_change_j")?;
        let work_into_interface_j = -checked_add(
            dot2(force_tangent, delta_displacement, "tangential_work")?,
            checked_mul(torsional_moment_nm, delta_twist, "torsional_work")?,
            "interface_work",
        )?;
        let closure_residual_j = checked_add(
            checked_add(
                checked_add(storage_change_j, dissipation_j, "physical_energy_change")?,
                numerical_dissipation_j,
                "accounted_energy_change",
            )?,
            -work_into_interface_j,
            "energy_closure_residual",
        )?;
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
        next_state.accepted_steps =
            state
                .accepted_steps
                .checked_add(1)
                .ok_or(PartialSlipError::InvalidDerived {
                    field: "accepted_steps overflow",
                })?;
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
                coefficient_work_j: 0.0,
                reconstructed_body_power_w,
                endpoint_body_power_w,
                numerical_dissipation_j,
                closure_residual_j,
            },
            next_state,
            checkpoint,
            applicability: PartialSlipApplicability::TwoRegionIwanReturnMap,
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
    /// Both regions remain elastic over this increment.
    Sticking,
    /// One region slips while the other remains elastic.
    PartialSlip,
    /// Both regions slip at their respective kinetic Coulomb capacities.
    GrossSlide,
}

/// Reversible internal state suitable for checkpointing and deterministic replay.
#[derive(Debug, Clone, PartialEq)]
pub struct PartialSlipState {
    region_elastic_displacement_m: [[f64; 2]; 2],
    region_elastic_twist_rad: [f64; 2],
    sliding_regions: [bool; 2],
    reference_frame: Option<TangentFrame>,
    accepted_steps: u64,
}

impl PartialSlipState {
    /// Empty history with no stored tangential energy.
    #[must_use]
    pub const fn zero() -> Self {
        Self {
            region_elastic_displacement_m: [[0.0; 2]; 2],
            region_elastic_twist_rad: [0.0; 2],
            sliding_regions: [false; 2],
            reference_frame: None,
            accepted_steps: 0,
        }
    }

    /// Reconstructs both regions' complete history, ordered weaker then harder.
    /// A single aggregate displacement cannot reconstruct partial-slip memory.
    /// The accepted reference frame is required for every nonempty history.
    pub fn from_region_checkpoint(
        region_elastic_displacement_m: [[f64; 2]; 2],
        region_elastic_twist_rad: [f64; 2],
        sliding_regions: [bool; 2],
        reference_frame: Option<TangentFrame>,
        accepted_steps: u64,
    ) -> Result<Self, PartialSlipError> {
        let value = Self {
            region_elastic_displacement_m,
            region_elastic_twist_rad,
            sliding_regions,
            reference_frame,
            accepted_steps,
        };
        value.validate()?;
        Ok(value)
    }

    /// Elastic tangent displacement of the harder region in m.
    /// Use [`Self::region_elastic_displacement_m`] for the complete checkpoint.
    #[must_use]
    pub const fn elastic_displacement_m(&self) -> [f64; 2] {
        self.region_elastic_displacement_m[1]
    }

    /// Elastic torsional displacement of the harder region in rad.
    #[must_use]
    pub const fn elastic_twist_rad(&self) -> f64 {
        self.region_elastic_twist_rad[1]
    }

    /// Elastic tangent displacement of both regions, weaker then harder [m].
    #[must_use]
    pub const fn region_elastic_displacement_m(&self) -> [[f64; 2]; 2] {
        self.region_elastic_displacement_m
    }

    /// Elastic twist of both regions, weaker then harder [rad].
    #[must_use]
    pub const fn region_elastic_twist_rad(&self) -> [f64; 2] {
        self.region_elastic_twist_rad
    }

    /// Retained kinetic branches, including a stationary pause during sliding.
    #[must_use]
    pub const fn sliding_regions(&self) -> [bool; 2] {
        self.sliding_regions
    }

    /// Last accepted frame in which both regional shear coordinates are stored.
    /// Only an empty, never-advanced state may omit this frame.
    #[must_use]
    pub const fn reference_frame(&self) -> Option<TangentFrame> {
        self.reference_frame
    }

    /// Number of accepted interval updates.
    #[must_use]
    pub const fn accepted_steps(&self) -> u64 {
        self.accepted_steps
    }

    fn validate(&self) -> Result<(), PartialSlipError> {
        for region in 0..2 {
            finite_vec2(
                self.region_elastic_displacement_m[region],
                "elastic_displacement_m",
            )?;
            finite(self.region_elastic_twist_rad[region], "elastic_twist_rad")?;
        }
        if self.reference_frame.is_none()
            && (self.accepted_steps != 0
                || self.sliding_regions.iter().any(|&sliding| sliding)
                || self
                    .region_elastic_twist_rad
                    .iter()
                    .any(|&twist| twist != 0.0)
                || self
                    .region_elastic_displacement_m
                    .iter()
                    .flatten()
                    .any(|&shear| shear != 0.0))
        {
            return Err(PartialSlipError::InvalidInput {
                field: "missing reference frame for retained history",
            });
        }
        Ok(())
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
// Bead frankensim-30cz7: `_id` suffixes are deliberate identity vocabulary;
// renaming would break caller-facing field semantics for a lint preference.
#[allow(clippy::struct_field_names)]
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
    /// Actual endpoint-wrench work into the interface, `-F·du - M dtheta` [J].
    /// Equals storage change + heat + numerical dissipation minus explicit
    /// coefficient work, up to roundoff.
    /// This is backward-Euler interval work, not an independently resolved
    /// continuous-time work integral.
    pub work_into_interface_j: f64,
    /// Signed work supplied by an explicit coefficient change at the OLD
    /// elastic history [J]. Zero for [`PartialSlipLaw::advance`]. This is an
    /// independent parameter-work port, not friction heat or wrench work.
    /// Temperature-dependent stiffness additionally requires the mechanical
    /// entropy/internal-energy account; this scalar alone is not caloric heat.
    pub coefficient_work_j: f64,
    /// `-work_into_interface_j / dt_s`, W; negative means body energy is dissipated/stored.
    pub reconstructed_body_power_w: f64,
    /// Endpoint wrench power `F·v + M·ω`, W, retained separately from discrete closure.
    pub endpoint_body_power_w: f64,
    /// Backward-Euler elastic integration loss [J], never exported as heat.
    /// This vanishes under refinement on smooth fixed-branch histories;
    /// a static-to-kinetic jump needs event resolution for that claim.
    pub numerical_dissipation_j: f64,
    /// Storage change + physical heat + numerical loss - endpoint-wrench work
    /// - coefficient work [J].
    pub closure_residual_j: f64,
}

/// Applicability identity for this generic rung.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PartialSlipApplicability {
    /// Two parallel elastic/Coulomb regions; no resolved Hertz traction field.
    TwoRegionIwanReturnMap,
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
    /// Independently measured discrete work, numerical defect, and ownership keys.
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
    MissingIdentity {
        /// The identity field that was blank.
        field: &'static str,
    },
    /// Input was non-finite, non-positive where required, or outside this rung's domain.
    InvalidInput {
        /// The refused input field name.
        field: &'static str,
    },
    /// A finite input produced a non-finite candidate.
    InvalidDerived {
        /// The derived-value field that was refused.
        field: &'static str,
    },
    /// Work keys named a different patch than the normal-patch view.
    WorkOwnershipMismatch,
    /// A checkpoint identity does not match current caller-owned identities.
    CheckpointIdentityMismatch {
        /// The checkpoint identity field that failed to bind.
        field: &'static str,
    },
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

struct RegionStep {
    elastic: [f64; 3],
    plastic: [f64; 3],
    effort: [f64; 3],
    kinetic_branch: bool,
    yielded: bool,
}

// Backward-Euler elastic/Coulomb update of one Jenkins region. The return is
// the maximum-dissipation projection in the elastic compliance metric, not a
// radial force clipping (which is incorrect when torsion and shear have
// different stiffness/capacity ratios).
fn advance_region(
    old: [f64; 3],
    increment: [f64; 3],
    stiffness: [f64; 3],
    limits: [f64; 3],
    kinetic_ratio: f64,
    was_sliding: bool,
) -> Result<RegionStep, PartialSlipError> {
    let mut trial = [0.0; 3];
    let mut effort = [0.0; 3];
    let mut normalized = [0.0; 3];
    let mut loading_work = 0.0;
    for i in 0..3 {
        positive(stiffness[i], "region_stiffness")?;
        positive(limits[i], "region_capacity")?;
        trial[i] = checked_add(old[i], increment[i], "region_trial_displacement")?;
        effort[i] = checked_mul(stiffness[i], trial[i], "region_trial_effort")?;
        normalized[i] = effort[i] / limits[i];
        finite(normalized[i], "region_trial_utilization")?;
        loading_work = checked_add(
            loading_work,
            checked_mul(stiffness[i] * old[i], increment[i], "region_loading_work")?,
            "region_loading_work",
        )?;
    }
    // A pause preserves the kinetic branch. Reversal unloads the springs
    // elastically before the opposite static capacity can be mobilized.
    let continuing = was_sliding && loading_work >= 0.0;
    let bound = if continuing { kinetic_ratio } else { 1.0 };
    let utilization = norm3(normalized)?;
    if utilization <= bound * (1.0 + EPSILON) {
        return Ok(RegionStep {
            elastic: trial,
            plastic: [0.0; 3],
            effort,
            kinetic_branch: continuing,
            yielded: false,
        });
    }
    effort = return_to_ellipse(effort, stiffness, limits, kinetic_ratio)?;
    let mut elastic = [0.0; 3];
    let mut plastic = [0.0; 3];
    for i in 0..3 {
        elastic[i] = effort[i] / stiffness[i];
        finite(elastic[i], "region_elastic_displacement")?;
        plastic[i] = checked_add(trial[i], -elastic[i], "region_plastic_displacement")?;
    }
    Ok(RegionStep {
        elastic,
        plastic,
        effort,
        kinetic_branch: true,
        yielded: true,
    })
}

fn norm3(value: [f64; 3]) -> Result<f64, PartialSlipError> {
    let norm = det::hypot(det::hypot(value[0], value[1]), value[2]);
    finite(norm, "generalized_norm")?;
    Ok(norm)
}

fn return_to_ellipse(
    trial: [f64; 3],
    stiffness: [f64; 3],
    limits: [f64; 3],
    kinetic_ratio: f64,
) -> Result<[f64; 3], PartialSlipError> {
    let mut normalized = [0.0; 3];
    for i in 0..3 {
        let capacity = limits[i] * kinetic_ratio;
        positive(capacity, "region_kinetic_capacity")?;
        normalized[i] = trial[i] / capacity;
        finite(normalized[i], "region_kinetic_utilization")?;
    }
    let norm = norm3(normalized)?;
    // The common pure-shear and pure-spin paths require no iteration.
    if trial[2] == 0.0 || (trial[0] == 0.0 && trial[1] == 0.0) {
        return Ok([trial[0] / norm, trial[1] / norm, trial[2] / norm]);
    }
    let mut rates = [0.0; 3];
    for i in 0..3 {
        rates[i] = stiffness[i] / limits[i] / limits[i];
        positive(rates[i], "region_projection_rate")?;
    }
    let scale = rates.iter().copied().fold(0.0_f64, f64::max);
    for rate in &mut rates {
        *rate /= scale;
        positive(*rate, "normalized_projection_rate")?;
    }
    // The KKT multiplier is the unique nonnegative root of
    // ||q_trial_i / (C_i (1 + lambda k_i/C_i^2))|| = 1.
    // Scale lambda so all rates <= 1, then use bracketed Newton steps.
    let mut lo = 0.0;
    let mut hi = (norm - 1.0) / rates.iter().copied().fold(1.0_f64, f64::min);
    positive(hi, "region_projection_bracket")?;
    let mut multiplier = norm - 1.0;
    for _ in 0..80 {
        let mut point = [0.0; 3];
        let mut denominators = [0.0; 3];
        for i in 0..3 {
            denominators[i] = 1.0 + multiplier * rates[i];
            finite(denominators[i], "region_projection_denominator")?;
            point[i] = normalized[i] / denominators[i];
        }
        let magnitude = norm3(point)?;
        if (magnitude - 1.0).abs() <= EPSILON {
            return Ok([
                trial[0] / denominators[0],
                trial[1] / denominators[1],
                trial[2] / denominators[2],
            ]);
        }
        if magnitude > 1.0 {
            lo = multiplier;
        } else {
            hi = multiplier;
        }
        let mut derivative = 0.0;
        for i in 0..3 {
            let direction = point[i] / magnitude;
            derivative -= magnitude * direction * direction * rates[i] / denominators[i];
        }
        let proposed = multiplier - (magnitude - 1.0) / derivative;
        multiplier = if proposed.is_finite() && proposed > lo && proposed < hi {
            proposed
        } else {
            f64::midpoint(lo, hi)
        };
    }
    Err(PartialSlipError::InvalidDerived {
        field: "region_projection_not_converged",
    })
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
            contact_spin_rad_per_s: 0.0,
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
        // The producer legitimately emits signed zeros in these channels
        // (e.g. a Coulomb force of `-capacity * signum(0.0)`); the original
        // contract is value equality (where -0.0 == 0.0), not bit equality.
        #[allow(clippy::float_cmp)]
        assert_eq!(response.tangent_force_n, [0.0; 2]);
        #[allow(clippy::float_cmp)]
        assert_eq!(response.torsional_moment_nm, 0.0);
        #[allow(clippy::float_cmp)]
        assert_eq!(response.dissipation.heat_j, 0.0);
        #[allow(clippy::float_cmp)]
        assert_eq!(response.dissipation.rolling_deformation_loss_j, 0.0);
        #[allow(clippy::float_cmp)]
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
        assert_eq!(
            gross.slip_partition.microslip_fraction.to_bits(),
            1.0f64.to_bits()
        );
        assert_eq!(
            gross.dissipation.rolling_deformation_loss_j.to_bits(),
            0.0f64.to_bits()
        );
    }

    /// G1 pure torsion uses the finite-patch pressure moment without inventing translation.
    #[test]
    fn pure_torsion_is_coulomb_limited_and_separate_from_rolling_loss() {
        let response = step(frame(), &PartialSlipState::zero(), [0.0; 2], 10.0);
        assert_eq!(response.state, PartialSlipStateKind::GrossSlide);
        #[allow(clippy::float_cmp)]
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
        #[allow(clippy::float_cmp)]
        assert_eq!(response.dissipation.rolling_deformation_loss_j, 0.0);
    }

    /// G0: independently computed wrench work closes with heat and numerical loss.
    #[test]
    fn work_heat_and_storage_close_exactly() {
        let response = step(frame(), &PartialSlipState::zero(), [6.0, 2.0], 2.0);
        assert!(response.dissipation.heat_j >= 0.0);
        close(
            response.generalized_work.work_into_interface_j,
            response.storage_change_j
                + response.dissipation.heat_j
                + response.generalized_work.numerical_dissipation_j,
        );
        close(
            response.generalized_work.work_into_interface_j,
            -0.001
                * (response.tangent_force_n[0] * 6.0
                    + response.tangent_force_n[1] * 2.0
                    + response.torsional_moment_nm * 2.0),
        );
        close(
            response.generalized_work.reconstructed_body_power_w,
            -response.generalized_work.work_into_interface_j / 0.001,
        );
        close(response.generalized_work.closure_residual_j, 0.0);
        close(
            response.generalized_work.endpoint_body_power_w,
            response.generalized_work.reconstructed_body_power_w,
        );
        assert_eq!(
            response.dissipation.heat_j.to_bits(),
            response
                .dissipation
                .tangential_and_torsional_microslip_j
                .to_bits()
        );
    }

    /// G1: purely elastic loading has numerical loss, but no frictional heat.
    #[test]
    fn elastic_loading_and_unloading_measure_actual_endpoint_work() {
        let loaded = step(frame(), &PartialSlipState::zero(), [1.0, 0.0], 0.0);
        close(loaded.stored_energy_j, 0.005);
        close(loaded.dissipation.heat_j, 0.0);
        close(loaded.generalized_work.work_into_interface_j, 0.01);
        close(loaded.generalized_work.numerical_dissipation_j, 0.005);
        close(loaded.generalized_work.closure_residual_j, 0.0);

        let unloaded = step(frame(), &loaded.next_state, [-0.5, 0.0], 0.0);
        close(unloaded.storage_change_j, -0.00375);
        close(unloaded.dissipation.heat_j, 0.0);
        close(unloaded.generalized_work.work_into_interface_j, -0.0025);
        close(unloaded.generalized_work.numerical_dissipation_j, 0.00125);
        close(unloaded.generalized_work.closure_residual_j, 0.0);
        assert!(unloaded.generalized_work.reconstructed_body_power_w > 0.0);
    }

    /// G1/G3: closed elastic cycles expose first-order algorithmic damping.
    #[test]
    fn elastic_cycle_numerical_loss_refines_without_becoming_heat() {
        fn cycle(steps_per_leg: usize) -> (f64, f64) {
            let mut state = PartialSlipState::zero();
            let mut work = 0.0;
            let mut numerical = 0.0;
            let dt = 0.001 / steps_per_leg as f64;
            for direction in [1.0, -1.0] {
                for _ in 0..steps_per_leg {
                    let response = law()
                        .advance(
                            &patch(),
                            &interface(),
                            frame(),
                            PartialSlipKinematics {
                                creepage: [direction, 0.0],
                                rolling_speed_mps: 1.0,
                                torsional_spin_rad_per_s: 0.0,
                                contact_spin_rad_per_s: 0.0,
                                dt_s: dt,
                            },
                            &ownership(),
                            &state,
                        )
                        .expect("elastic cycle");
                    assert_eq!(response.state, PartialSlipStateKind::Sticking);
                    close(response.dissipation.heat_j, 0.0);
                    work += response.generalized_work.work_into_interface_j;
                    numerical += response.generalized_work.numerical_dissipation_j;
                    state = response.next_state;
                }
            }
            for regional_displacement in state.region_elastic_displacement_m() {
                close(regional_displacement[0], 0.0);
            }
            close(work, numerical);
            (work, numerical)
        }
        let coarse = cycle(4);
        let medium = cycle(8);
        let fine = cycle(16);
        close(coarse.1, 0.01 / 4.0);
        close(coarse.1, 2.0 * medium.1);
        close(medium.1, 2.0 * fine.1);
    }

    /// G0/G3: combined translation/spin reversals never hide or heat numerical loss.
    #[test]
    fn mixed_reversal_history_closes_against_world_wrench_work() {
        let rotated = frame().rotated(0.43).expect("rotated frame");
        let mut state = PartialSlipState::zero();
        let mut total_work = 0.0;
        let mut total_heat = 0.0;
        let mut total_numerical = 0.0;
        let mut final_storage = 0.0;
        for (creepage, spin) in [
            ([1.0, -0.5], 0.1),
            ([6.0, 2.0], 2.0),
            ([-20.0, 0.5], -10.0),
            ([0.0, 0.0], 0.0),
            ([10.0, -10.0], 5.0),
            ([-0.25, 0.5], -0.2),
        ] {
            let response = step(rotated, &state, creepage, spin);
            let velocity = world(rotated, creepage);
            let body_power = (0..3)
                .map(|i| {
                    response.wrench.force_n[i] * velocity[i]
                        + response.wrench.torque_nm[i] * rotated.normal()[i] * spin
                })
                .sum::<f64>();
            let work = -body_power * 0.001;
            close(response.generalized_work.work_into_interface_j, work);
            close(response.generalized_work.closure_residual_j, 0.0);
            assert!(response.generalized_work.numerical_dissipation_j >= 0.0);
            assert!(response.dissipation.heat_j >= 0.0);
            total_work += work;
            total_heat += response.dissipation.heat_j;
            total_numerical += response.generalized_work.numerical_dissipation_j;
            final_storage = response.stored_energy_j;
            state = response.next_state;
        }
        close(total_work, final_storage + total_heat + total_numerical);
        assert!(total_numerical > 0.0);
    }

    /// G1: steady saturated sliding transfers real frictional work entirely to heat.
    #[test]
    fn steady_gross_sliding_has_heat_without_numerical_storage_loss() {
        let first = step(frame(), &PartialSlipState::zero(), [20.0, 0.0], 0.0);
        let steady = step(frame(), &first.next_state, [20.0, 0.0], 0.0);
        close(steady.storage_change_j, 0.0);
        close(steady.generalized_work.numerical_dissipation_j, 0.0);
        close(steady.dissipation.heat_j, 40.0 * 0.02);
        close(
            steady.generalized_work.work_into_interface_j,
            steady.dissipation.heat_j,
        );
        close(steady.generalized_work.closure_residual_j, 0.0);
    }

    /// G1: rate-independent partial slip cannot relax or heat a held contact.
    #[test]
    fn loaded_partial_contact_holds_without_relaxation_or_heat() {
        let first = step(frame(), &PartialSlipState::zero(), [6.0, 0.0], 0.0);
        assert_eq!(first.state, PartialSlipStateKind::PartialSlip);
        let mut state = first.next_state.clone();
        for _ in 0..32 {
            let held = step(frame(), &state, [0.0; 2], 0.0);
            assert_eq!(
                held.next_state.region_elastic_displacement_m(),
                first.next_state.region_elastic_displacement_m()
            );
            assert_eq!(
                held.next_state.region_elastic_twist_rad(),
                first.next_state.region_elastic_twist_rad()
            );
            close(held.tangent_force_n[0], first.tangent_force_n[0]);
            close(held.stored_energy_j, first.stored_energy_j);
            assert_eq!(held.dissipation.heat_j, 0.0);
            assert_eq!(held.generalized_work.work_into_interface_j, 0.0);
            assert_eq!(held.generalized_work.numerical_dissipation_j, 0.0);
            state = held.next_state;
        }
        // Resuming the same direction must retain sliding in the weak region;
        // a pause cannot spuriously restore its static force ramp.
        let resumed = step(frame(), &state, [0.01, 0.0], 0.0);
        close(resumed.tangent_force_n[0] - first.tangent_force_n[0], -0.04);
        assert_eq!(resumed.state, PartialSlipStateKind::PartialSlip);
    }

    /// G1: reversal first unloads both retained springs elastically.
    #[test]
    fn partial_slip_reversal_recovers_storage_before_opposite_sliding() {
        let first = step(frame(), &PartialSlipState::zero(), [6.0, 0.0], 0.0);
        let reverse = step(frame(), &first.next_state, [-0.1, 0.0], 0.0);
        assert_eq!(reverse.state, PartialSlipStateKind::Sticking);
        close(reverse.tangent_force_n[0] - first.tangent_force_n[0], 1.0);
        assert_eq!(reverse.dissipation.heat_j, 0.0);
        assert!(reverse.storage_change_j < 0.0);
        assert!(reverse.generalized_work.work_into_interface_j < 0.0);
        close(reverse.generalized_work.closure_residual_j, 0.0);
    }

    /// G1: the elliptic return satisfies maximum-dissipation normality with
    /// unequal translational/torsional compliance, rather than radial clipping.
    #[test]
    fn mixed_spin_and_translation_return_in_the_compliance_metric() {
        let response = step(frame(), &PartialSlipState::zero(), [20.0, 10.0], 10.0);
        assert_eq!(response.state, PartialSlipStateKind::GrossSlide);
        let weights = [0.6, 0.4];
        let fractions = [0.3, 0.7];
        for region in 0..2 {
            let e = response.next_state.region_elastic_displacement_m()[region];
            let twist = response.next_state.region_elastic_twist_rad()[region];
            let force = 10_000.0 * weights[region] * e[0];
            let torque = 100.0 * weights[region] * twist;
            let force_cap = response.capacity.kinetic_force_n * fractions[region];
            let torque_cap = response.capacity.kinetic_torque_nm * fractions[region];
            close(
                (0.02 - e[0]) * force_cap * force_cap / force,
                (0.01 - twist) * torque_cap * torque_cap / torque,
            );
        }
        assert!(response.dissipation.heat_j > 0.0);
        close(response.generalized_work.closure_residual_j, 0.0);
    }

    /// G1/G3: monotone dry sliding has analytic slider heat; the numerical
    /// spring loss decreases under time refinement and never becomes heat.
    #[test]
    fn slider_heat_and_force_are_path_consistent_under_time_refinement() {
        let law = PartialSlipLaw::new(
            PARTIAL_SLIP_MODEL_ID,
            "synthetic/equal-coulomb",
            PartialSlipParameters {
                kinetic_mu: 0.8,
                ..law().parameters()
            },
        )
        .unwrap();
        let mut numerical = Vec::new();
        for steps in [32, 64, 128] {
            let mut state = PartialSlipState::zero();
            let mut heat = 0.0;
            let mut work = 0.0;
            let mut loss = 0.0;
            let mut energy = 0.0;
            for _ in 0..steps {
                let response = law
                    .advance(
                        &patch(),
                        &interface(),
                        frame(),
                        PartialSlipKinematics {
                            creepage: [0.02, 0.0],
                            rolling_speed_mps: 1.0,
                            torsional_spin_rad_per_s: 0.0,
                            contact_spin_rad_per_s: 0.0,
                            dt_s: 1.0 / f64::from(steps),
                        },
                        &ownership(),
                        &state,
                    )
                    .unwrap();
                heat += response.dissipation.heat_j;
                work += response.generalized_work.work_into_interface_j;
                loss += response.generalized_work.numerical_dissipation_j;
                energy = response.stored_energy_j;
                state = response.next_state;
            }
            // Parallel sliders: (k,C) = (6000,24), (4000,56).
            close(
                heat,
                24.0 * (0.02 - 24.0 / 6000.0) + 56.0 * (0.02 - 56.0 / 4000.0),
            );
            close(energy, 0.5 * (24.0 * 24.0 / 6000.0 + 56.0 * 56.0 / 4000.0));
            close(work, energy + heat + loss);
            numerical.push(loss);
        }
        assert!(numerical[1] < 0.6 * numerical[0]);
        assert!(numerical[2] < 0.6 * numerical[1]);
    }

    /// G5: both regional strains and kinetic branches survive serialization.
    #[test]
    fn regional_checkpoint_preserves_partial_slip_memory() {
        let initial_frame = frame().rotated(0.7).unwrap();
        let first = step(initial_frame, &PartialSlipState::zero(), [6.0, 0.0], 0.0);
        let state = &first.next_state;
        let restored = PartialSlipState::from_region_checkpoint(
            state.region_elastic_displacement_m(),
            state.region_elastic_twist_rad(),
            state.sliding_regions(),
            state.reference_frame(),
            state.accepted_steps(),
        )
        .unwrap();
        assert_eq!(&restored, state);
        assert_eq!(
            step(frame(), &restored, [0.1, 0.2], 0.1),
            step(frame(), state, [0.1, 0.2], 0.1)
        );
    }

    /// Warm-state gauge changes cannot rotate physical preload or change work.
    #[test]
    fn changing_tangent_gauge_preserves_loaded_history_and_following_work() {
        let first = step(frame(), &PartialSlipState::zero(), [6.0, 2.0], 1.5);
        let mut fixed_state = first.next_state.clone();
        let mut changing_state = first.next_state;
        for (creepage, spin, angle) in [
            ([0.0, 0.0], 0.0, core::f64::consts::FRAC_PI_2),
            ([0.5, -0.2], 0.1, 0.37),
            ([-8.0, 1.0], -0.4, -0.81),
            ([20.0, 4.0], 1.2, 1.14),
        ] {
            let rotated_frame = frame().rotated(angle).unwrap();
            let cosine = angle.cos(); // det-ok: independent coordinate oracle.
            let sine = angle.sin(); // det-ok: independent coordinate oracle.
            let rotated_creepage = [
                cosine * creepage[0] + sine * creepage[1],
                -sine * creepage[0] + cosine * creepage[1],
            ];
            let fixed = step(frame(), &fixed_state, creepage, spin);
            let changing = step(rotated_frame, &changing_state, rotated_creepage, spin);
            close_vec(changing.wrench.force_n, fixed.wrench.force_n);
            close_vec(changing.wrench.torque_nm, fixed.wrench.torque_nm);
            close(changing.stored_energy_j, fixed.stored_energy_j);
            close(changing.storage_change_j, fixed.storage_change_j);
            close(changing.dissipation.heat_j, fixed.dissipation.heat_j);
            close(
                changing.generalized_work.numerical_dissipation_j,
                fixed.generalized_work.numerical_dissipation_j,
            );
            close(
                changing.generalized_work.work_into_interface_j,
                fixed.generalized_work.work_into_interface_j,
            );
            close(changing.generalized_work.closure_residual_j, 0.0);
            fixed_state = fixed.next_state;
            changing_state = changing.next_state;
        }
    }

    /// Normal tilt transports shear isometrically, without projection damping.
    #[test]
    fn changing_normal_transports_preload_without_heat_or_storage_loss() {
        let first = step(frame(), &PartialSlipState::zero(), [6.0, 2.0], 1.5);
        let tilted = TangentFrame::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]).unwrap();
        let held = step(tilted, &first.next_state, [0.0; 2], 0.0);
        // Minimum rotation from +z to +x is a right-angle rotation about +y.
        close_vec(
            held.wrench.force_n,
            [0.0, first.wrench.force_n[1], -first.wrench.force_n[0]],
        );
        close_vec(held.wrench.torque_nm, [first.torsional_moment_nm, 0.0, 0.0]);
        close(held.stored_energy_j, first.stored_energy_j);
        close(held.storage_change_j, 0.0);
        close(held.dissipation.heat_j, 0.0);
        close(held.generalized_work.numerical_dissipation_j, 0.0);
        close(held.generalized_work.work_into_interface_j, 0.0);
        let returned = step(frame(), &held.next_state, [0.0; 2], 0.0);
        close_vec(returned.wrench.force_n, first.wrench.force_n);
        close_vec(returned.wrench.torque_nm, first.wrench.torque_nm);
        close(returned.stored_energy_j, first.stored_energy_j);
    }

    /// Common material spin rotates shear; only relative spin loads torsion.
    #[test]
    fn common_contact_spin_rotates_preload_without_loading_torsion() {
        let first = step(frame(), &PartialSlipState::zero(), [6.0, 2.0], 1.5);
        let held = law()
            .advance(
                &patch(),
                &interface(),
                frame(),
                PartialSlipKinematics {
                    contact_spin_rad_per_s: core::f64::consts::FRAC_PI_2 / 0.001,
                    ..kinematics([0.0; 2], 0.0)
                },
                &ownership(),
                &first.next_state,
            )
            .unwrap();
        close_vec(
            held.wrench.force_n,
            [-first.wrench.force_n[1], first.wrench.force_n[0], 0.0],
        );
        close_vec(held.wrench.torque_nm, first.wrench.torque_nm);
        close(held.stored_energy_j, first.stored_energy_j);
        close(held.storage_change_j, 0.0);
        close(held.dissipation.heat_j, 0.0);
        close(held.generalized_work.numerical_dissipation_j, 0.0);
        close(held.generalized_work.work_into_interface_j, 0.0);
    }

    /// Transport commutes with a rigid rotation of the complete physical setup.
    #[test]
    fn warm_history_transport_is_covariant_under_three_dimensional_rotation() {
        let rotate = |value: [f64; 3]| [value[2], value[0], value[1]];
        let rotate_frame = |value: TangentFrame| {
            TangentFrame::new(rotate(value.normal()), rotate(value.longitudinal())).unwrap()
        };
        let first = step(frame(), &PartialSlipState::zero(), [6.0, 2.0], 1.5);
        let first_rotated = step(
            rotate_frame(frame()),
            &PartialSlipState::zero(),
            [6.0, 2.0],
            1.5,
        );
        let tilted = TangentFrame::new([0.3, 0.4, 0.5], [0.6, -0.2, 0.1]).unwrap();
        let input = PartialSlipKinematics {
            contact_spin_rad_per_s: 30.0,
            ..kinematics([0.2, -0.1], -0.2)
        };
        let next = law()
            .advance(
                &patch(),
                &interface(),
                tilted,
                input,
                &ownership(),
                &first.next_state,
            )
            .unwrap();
        let rotated = law()
            .advance(
                &patch(),
                &interface(),
                rotate_frame(tilted),
                input,
                &ownership(),
                &first_rotated.next_state,
            )
            .unwrap();
        close_vec(rotated.wrench.force_n, rotate(next.wrench.force_n));
        close_vec(rotated.wrench.torque_nm, rotate(next.wrench.torque_nm));
        close(rotated.stored_energy_j, next.stored_energy_j);
        close(rotated.dissipation.heat_j, next.dissipation.heat_j);
        close(
            rotated.generalized_work.work_into_interface_j,
            next.generalized_work.work_into_interface_j,
        );
    }

    #[test]
    fn surface_order_reversal_preserves_transported_physical_wrench() {
        let reverse_frame = |value: TangentFrame| {
            TangentFrame::new(
                value.normal().map(|component| -component),
                value.longitudinal(),
            )
            .unwrap()
        };
        let first = step(frame(), &PartialSlipState::zero(), [6.0, 2.0], 1.5);
        let first_reversed = step(
            reverse_frame(frame()),
            &PartialSlipState::zero(),
            [-6.0, 2.0],
            1.5,
        );
        let tilted = TangentFrame::new([0.3, 0.4, 0.5], [0.6, -0.2, 0.1]).unwrap();
        let input = PartialSlipKinematics {
            contact_spin_rad_per_s: 30.0,
            ..kinematics([0.2, -0.1], -0.2)
        };
        let next = law()
            .advance(
                &patch(),
                &interface(),
                tilted,
                input,
                &ownership(),
                &first.next_state,
            )
            .unwrap();
        // Reversing surfaces negates both the normal and the relative world
        // velocity. The common-spin scalar negates; relative normal spin does not.
        let reversed_input = PartialSlipKinematics {
            creepage: [-input.creepage[0], input.creepage[1]],
            contact_spin_rad_per_s: -input.contact_spin_rad_per_s,
            ..input
        };
        let reversed = law()
            .advance(
                &patch(),
                &interface(),
                reverse_frame(tilted),
                reversed_input,
                &ownership(),
                &first_reversed.next_state,
            )
            .unwrap();
        close_vec(
            reversed.wrench.force_n,
            next.wrench.force_n.map(|value| -value),
        );
        close_vec(
            reversed.wrench.torque_nm,
            next.wrench.torque_nm.map(|value| -value),
        );
        close(reversed.stored_energy_j, next.stored_energy_j);
        close(reversed.dissipation.heat_j, next.dissipation.heat_j);
        close(
            reversed.generalized_work.work_into_interface_j,
            next.generalized_work.work_into_interface_j,
        );
    }

    #[test]
    fn nonempty_checkpoint_requires_frame_and_antipodal_transport_refuses() {
        let first = step(frame(), &PartialSlipState::zero(), [6.0, 2.0], 1.5);
        let state = &first.next_state;
        assert!(
            PartialSlipState::from_region_checkpoint(
                state.region_elastic_displacement_m(),
                state.region_elastic_twist_rad(),
                state.sliding_regions(),
                None,
                state.accepted_steps(),
            )
            .is_err()
        );
        let reversed_normal = TangentFrame::new([0.0, 0.0, -1.0], [1.0, 0.0, 0.0]).unwrap();
        assert_eq!(
            law().advance(
                &patch(),
                &interface(),
                reversed_normal,
                kinematics([0.0; 2], 0.0),
                &ownership(),
                state,
            ),
            Err(PartialSlipError::InvalidInput {
                field: "antipodal contact-normal transport",
            }),
        );
        assert_eq!(state.reference_frame(), Some(frame()));
    }

    /// G3 reversal is deterministic and changes the gross-slide force direction.
    #[test]
    fn reversal_preserves_history_then_reverses_resistance() {
        let first = step(frame(), &PartialSlipState::zero(), [20.0, 0.0], 0.0);
        let reversed = step(frame(), &first.next_state, [-30.0, 0.0], 0.0);
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
    // Bead frankensim-30cz7: this G3 refusal battery intentionally enumerates
    // every mutation class inline; extraction would obscure the case table.
    #[allow(clippy::too_many_lines)]
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
            contact_spin_rad_per_s: 0.0,
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
        let exhausted = PartialSlipState::from_region_checkpoint(
            [[0.0; 2]; 2],
            [0.0; 2],
            [false; 2],
            Some(frame()),
            u64::MAX,
        )
        .expect("state");
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
