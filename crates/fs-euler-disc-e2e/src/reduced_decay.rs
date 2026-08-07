//! Small-angle, late-stage Euler-disc reduced decay.
//!
//! This numerical reference uses the declared relations
//! `Omega^2 sin(theta) = 4 g / R` and `E approximately equals 3 m g R theta / 2`.
//! It stops at a caller-declared positive validity cutoff before `theta = 0`.
//! It is neither a full rigid-body/contact solve nor a fit to a video or Mould
//! outcome, thin-gap solution, or resolved CFD. The generic dry channel calls
//! `fs-tribo` directly. A separately named source-bound rolling-power channel
//! preserves a published closure without inheriting the generic contour-speed
//! convention. It and the Bildsten boundary-layer channel are energy-only and
//! supply no force/wrench.
//!
//! Exponent recovery and discrete `P * dt` energy closure test only numerical
//! self-consistency of these stated equations. They are neither independent
//! emergence nor experimental or physical validation.

use core::fmt;

use fs_exec::Cx;
use fs_rep_frep::SquatDiscEdgeTreatment;
use fs_tribo::{
    ConstantContourForce, InputAuthority, InterfaceMedium, InterfaceSystemRef, ResistanceInput,
    ResistanceLaw,
};

use crate::specimen::{DiscProfileSpec, ResolvedDiscProfile};

/// Standard gravitational acceleration [m/s^2].
pub const STANDARD_GRAVITY_M_PER_S2: f64 = 9.806_65;
/// Retention bound for deterministic one-shot integration.
pub const MAX_REDUCED_DECAY_STEPS: u32 = 200_000;
/// Declared small-angle applicability ceiling [rad].
///
/// This is an input refusal boundary for the reduced reference, not a
/// certification that every state below it is physically resolved.
pub const MAX_SMALL_ANGLE_THETA_RAD: f64 = 0.2;
/// Coefficient in the published Bildsten energy-only power law.
///
/// With a caller multiplier of one, the channel evaluates
/// `4 sqrt(mu rho) g^(5/4) R^(11/4) theta^(-5/4)` [W].
pub const BILDSTEN_PUBLISHED_POWER_COEFFICIENT: f64 = 4.0;
/// Stable model identity retained in every numerical-reference run.
pub const REDUCED_DECAY_MODEL_ID: &str = "euler-disc-small-angle-late-stage-v1";
/// Fixed bisection count for the encoded-model channel-crossover diagnostic.
pub const CHANNEL_CROSSOVER_BISECTION_STEPS: u32 = 64;
/// Source identity for the Thorne et al. steel-on-glass benchmark declaration.
pub const THORNE_2026_SOURCE_ID: &str = "arxiv:2603.14520v1";
/// Reported diameter of the Table S1 steel disc [m].
pub const THORNE_2026_STEEL_DISC_DIAMETER_M: f64 = 0.075;
/// Reported axial thickness of the Table S1 steel disc [m].
pub const THORNE_2026_STEEL_DISC_THICKNESS_M: f64 = 0.013;
/// Reported mass of the Table S1 steel disc [kg].
pub const THORNE_2026_STEEL_DISC_MASS_KG: f64 = 0.445;
/// Reported circular outer-edge fillet radius [m].
pub const THORNE_2026_STEEL_DISC_FILLET_RADIUS_M: f64 = 0.001_6;
/// Ambient density used for the Figure 3 analytical comparison [kg/m^3].
pub const THORNE_2026_AMBIENT_AIR_DENSITY_KG_PER_M3: f64 = 1.18;
/// Partial-vacuum density used for the Figure 3 analytical comparison [kg/m^3].
pub const THORNE_2026_VACUUM_AIR_DENSITY_KG_PER_M3: f64 = 0.118;
/// Fitted rolling coefficient used for the Figure 3 analytical comparison.
pub const THORNE_2026_FITTED_ROLLING_COEFFICIENT: f64 = 1.0e-4;
/// Standard-air dynamic viscosity pinned by this benchmark [Pa s].
///
/// This is a declared analytical-model input, not a reported measurement of
/// the experiment's chamber air.
pub const THORNE_2026_DECLARED_AIR_VISCOSITY_PA_S: f64 = 1.8e-5;
/// Rendering-oriented initial inclination [rad]. This is a FrankenSim choice,
/// not a digitized initial point from the paper.
pub const THORNE_2026_RENDER_INITIAL_THETA_RAD: f64 = 0.12;
/// Positive analytical validity cutoff [rad]. No `theta = 0` or loss-of-contact
/// event is asserted at this value.
pub const THORNE_2026_RENDER_VALIDITY_CUTOFF_THETA_RAD: f64 = 0.003;
/// Fixed integration step used by the rendering benchmark [s].
pub const THORNE_2026_RENDER_TIMESTEP_S: f64 = 1.0e-4;
/// Bounded retained-step budget used by the rendering benchmark.
pub const THORNE_2026_RENDER_MAXIMUM_STEPS: u32 = 100_000;

/// Refusal from the bounded late-stage reference model.
#[derive(Debug, Clone, PartialEq)]
pub enum ReducedDecayError {
    /// A scalar input was non-finite or outside its declared domain.
    InvalidInput { field: &'static str },
    /// A required source identity was blank.
    MissingIdentity { field: &'static str },
    /// The initial state is not above the positive validity cutoff.
    InitialStateOutsideValidity {
        theta_rad: f64,
        cutoff_theta_rad: f64,
    },
    /// No active dissipation channel was declared.
    NoActiveChannel,
    /// The real dry-law boundary refused the supplied interface/input.
    DryLawRefusal { detail: String },
    /// The exact filleted benchmark profile could not be resolved.
    ProfileRefusal { detail: String },
    /// A checked arithmetic result was not finite.
    NonFiniteDerived { field: &'static str },
    /// Refinement would exceed the retained step bound.
    RefinementStepBudgetOverflow,
    /// Coarse and fine integrations reached different terminal classes.
    RefinementTerminalMismatch {
        /// Requested-step terminal class.
        coarse: ReducedDecayTerminal,
        /// Half-step terminal class.
        fine: ReducedDecayTerminal,
    },
    /// Matching refinement terminals did not reach the declared validity cutoff.
    RefinementIncompleteTerminal {
        /// Shared terminal class.
        terminal: ReducedDecayTerminal,
    },
    /// A caller supplied a malformed run with no retained samples.
    MissingSample,
}

impl fmt::Display for ReducedDecayError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{self:?}")
    }
}

impl std::error::Error for ReducedDecayError {}

fn finite_positive(value: f64, field: &'static str) -> Result<(), ReducedDecayError> {
    if !(value.is_finite() && value > 0.0) {
        return Err(ReducedDecayError::InvalidInput { field });
    }
    Ok(())
}

fn nonblank(value: &str, field: &'static str) -> Result<(), ReducedDecayError> {
    if value.trim().is_empty() {
        return Err(ReducedDecayError::MissingIdentity { field });
    }
    Ok(())
}

/// Caller-declared dry contour channel.  Its force is supplied to the generic
/// `fs-tribo::ConstantContourForce` law; it is not calibrated here.
#[derive(Debug, Clone, PartialEq)]
pub struct DryContourChannel {
    /// Ordered dry interface and history identity.
    pub interface: InterfaceSystemRef,
    /// Caller-declared normal force [N] retained by `fs-tribo` input checking.
    pub normal_force_n: f64,
    /// Caller-declared contour resistance magnitude [N].
    pub contour_force_n: f64,
}

impl DryContourChannel {
    fn validate(&self) -> Result<(), ReducedDecayError> {
        finite_positive(self.normal_force_n, "dry_contour.normal_force_n")?;
        finite_positive(self.contour_force_n, "dry_contour.contour_force_n")
    }
}

/// Direct source-bound rolling-power closure.
///
/// This channel evaluates `Phi_roll = mu m g R Omega` exactly. It is kept
/// separate from [`DryContourChannel`] because routing the same coefficient
/// through [`ConstantContourForce`] would multiply by the generic contour
/// speed `R Omega cos(theta)` and silently change the declared equation. Like
/// the boundary-layer channel, this is energy-only and supplies no wrench.
#[derive(Debug, Clone, PartialEq)]
pub struct PublishedRollingPowerChannel {
    /// Literature or declaration identity for this exact closure.
    pub source_id: String,
    /// Dimensionless fitted rolling coefficient `mu`.
    pub coefficient_mu: f64,
}

impl PublishedRollingPowerChannel {
    fn validate(&self) -> Result<(), ReducedDecayError> {
        nonblank(&self.source_id, "published_rolling.source_id")?;
        finite_positive(self.coefficient_mu, "published_rolling.coefficient_mu")
    }

    fn power_w(
        &self,
        mass_kg: f64,
        gravity_m_per_s2: f64,
        radius_m: f64,
        omega_rad_s: f64,
    ) -> Result<f64, ReducedDecayError> {
        let power = self.coefficient_mu * mass_kg * gravity_m_per_s2 * radius_m * omega_rad_s;
        if !(power.is_finite() && power > 0.0) {
            return Err(ReducedDecayError::NonFiniteDerived {
                field: "published_rolling.power_w",
            });
        }
        Ok(power)
    }
}

/// Source-bound geometry and mass declaration for a literature benchmark.
#[derive(Debug, Clone, PartialEq)]
pub struct LiteratureDiscSpecimen {
    /// Exact literature-version identity.
    pub source_id: String,
    /// Reported disc diameter [m].
    pub diameter_m: f64,
    /// Reported disc thickness [m].
    pub thickness_m: f64,
    /// Reported disc mass [kg].
    pub mass_kg: f64,
    /// Reported outer-edge circular fillet radius [m].
    pub outer_fillet_radius_m: f64,
}

impl LiteratureDiscSpecimen {
    fn validate(&self) -> Result<(), ReducedDecayError> {
        nonblank(&self.source_id, "literature_specimen.source_id")?;
        finite_positive(self.diameter_m, "literature_specimen.diameter_m")?;
        finite_positive(self.thickness_m, "literature_specimen.thickness_m")?;
        finite_positive(self.mass_kg, "literature_specimen.mass_kg")?;
        finite_positive(
            self.outer_fillet_radius_m,
            "literature_specimen.outer_fillet_radius_m",
        )?;
        if 2.0 * self.outer_fillet_radius_m >= self.thickness_m
            || self.outer_fillet_radius_m >= 0.5 * self.diameter_m
        {
            return Err(ReducedDecayError::InvalidInput {
                field: "literature_specimen.fillet_geometry",
            });
        }
        Ok(())
    }

    /// Exact squat-cylinder line/arc profile named by the declaration.
    #[must_use]
    pub fn profile_spec(&self) -> DiscProfileSpec {
        DiscProfileSpec::SolidCylinder {
            outer_radius_m: 0.5 * self.diameter_m,
            thickness_m: self.thickness_m,
            edge_treatment: SquatDiscEdgeTreatment::CircularFillet {
                radius: self.outer_fillet_radius_m,
            },
        }
    }

    /// Resolve the reported shape using the homogeneous density required to
    /// reproduce its reported mass.
    ///
    /// That inferred density is a geometry/mass bookkeeping value, not a
    /// material-characterization claim about the machined steel.
    pub fn resolve_mass_matched_profile(
        &self,
        cx: &Cx<'_>,
    ) -> Result<ResolvedDiscProfile, ReducedDecayError> {
        self.validate()?;
        let unit_density = self.profile_spec().resolve(1.0, cx).map_err(|error| {
            ReducedDecayError::ProfileRefusal {
                detail: error.to_string(),
            }
        })?;
        let density_kg_per_m3 = self.mass_kg / unit_density.mass_properties.volume;
        self.profile_spec()
            .resolve(density_kg_per_m3, cx)
            .map_err(|error| ReducedDecayError::ProfileRefusal {
                detail: error.to_string(),
            })
    }
}

/// Bildsten-style rotating-disc boundary-layer energy closure.
///
/// Its direct energy-only law is `C * 4 sqrt(mu rho) g^(5/4) R^(11/4)
/// theta^(-5/4)`, where `C = 1` is the published coefficient convention.
/// Under the declared late-stage relation, this is proportional to
/// `theta^(-5/4)`, yielding the reference inclination exponent `4/9` while
/// this closure remains applicable.
#[derive(Debug, Clone, PartialEq)]
pub struct BildstenBoundaryLayerChannel {
    /// Named correlation/source identity; no authority is upgraded here.
    pub source_id: String,
    /// Gas density [kg/m^3].
    pub density_kg_per_m3: f64,
    /// Dynamic viscosity [Pa s].
    pub dynamic_viscosity_pa_s: f64,
    /// Caller-declared positive multiplier of the published law. A multiplier
    /// of one means the cited coefficient is used exactly; it is not fitted to
    /// an Euler-disc outcome.
    pub dimensionless_prefactor: f64,
}

impl BildstenBoundaryLayerChannel {
    fn validate(&self) -> Result<(), ReducedDecayError> {
        nonblank(&self.source_id, "bildsten.source_id")?;
        finite_positive(self.density_kg_per_m3, "bildsten.density_kg_per_m3")?;
        finite_positive(
            self.dynamic_viscosity_pa_s,
            "bildsten.dynamic_viscosity_pa_s",
        )?;
        finite_positive(
            self.dimensionless_prefactor,
            "bildsten.dimensionless_prefactor",
        )
    }

    fn power_w(
        &self,
        radius_m: f64,
        gravity_m_per_s2: f64,
        theta_rad: f64,
    ) -> Result<f64, ReducedDecayError> {
        let power = self.dimensionless_prefactor
            * BILDSTEN_PUBLISHED_POWER_COEFFICIENT
            * (self.dynamic_viscosity_pa_s * self.density_kg_per_m3).sqrt()
            * gravity_m_per_s2.powf(1.25)
            * radius_m.powf(2.75)
            * theta_rad.powf(-1.25);
        if !(power.is_finite() && power >= 0.0) {
            return Err(ReducedDecayError::NonFiniteDerived {
                field: "bildsten.power_w",
            });
        }
        Ok(power)
    }
}

/// Bounded input for the late-stage reference integration.
#[derive(Debug, Clone, PartialEq)]
pub struct ReducedDecayInput {
    /// Disc mass [kg].
    pub mass_kg: f64,
    /// Disc radius [m].
    pub radius_m: f64,
    /// Declared gravitational acceleration [m/s^2].
    pub gravity_m_per_s2: f64,
    /// Initial inclination [rad], strictly above the cutoff and at most the
    /// declared small-angle ceiling.
    pub initial_theta_rad: f64,
    /// Positive terminal validity cutoff [rad].
    pub validity_cutoff_theta_rad: f64,
    /// Fixed requested integration step [s].
    pub timestep_s: f64,
    /// Hard maximum retained integration steps.
    pub maximum_steps: u32,
    /// Nonblank identity for the small-angle energy/precession oracle. This is
    /// distinct from a dissipation-law source and upgrades no authority.
    pub small_angle_oracle_source_id: String,
    /// Optional real dry contour channel.
    pub dry_contour: Option<DryContourChannel>,
    /// Optional energy-only boundary-layer channel.
    pub bildsten_boundary_layer: Option<BildstenBoundaryLayerChannel>,
}

impl ReducedDecayInput {
    /// A named numerical-reference setup with both channels active.  It is not
    /// a material card, physical prediction, or outcome fit.
    pub fn nominal_reference() -> Result<Self, ReducedDecayError> {
        Ok(Self {
            mass_kg: 0.12,
            radius_m: 0.038,
            gravity_m_per_s2: STANDARD_GRAVITY_M_PER_S2,
            initial_theta_rad: 0.05,
            validity_cutoff_theta_rad: 0.001,
            timestep_s: 1.0e-5,
            maximum_steps: 100_000,
            small_angle_oracle_source_id: "analytic/euler-disc-small-angle-oracle-v1".to_owned(),
            dry_contour: Some(DryContourChannel {
                interface: InterfaceSystemRef::new(
                    "reduced-decay/disc->support",
                    "reduced-decay/reference-history",
                    "caller-declared/numerical-reference",
                    InputAuthority::CallerDeclared,
                    InterfaceMedium::Dry,
                )
                .map_err(|error| ReducedDecayError::DryLawRefusal {
                    detail: error.to_string(),
                })?,
                normal_force_n: 0.12 * STANDARD_GRAVITY_M_PER_S2,
                contour_force_n: 0.002,
            }),
            bildsten_boundary_layer: Some(BildstenBoundaryLayerChannel {
                source_id: "doi:10.1103/PhysRevE.66.056309".to_owned(),
                density_kg_per_m3: 1.2,
                dynamic_viscosity_pa_s: 1.8e-5,
                dimensionless_prefactor: 1.0,
            }),
        })
    }

    fn validate(&self) -> Result<(), ReducedDecayError> {
        finite_positive(self.mass_kg, "mass_kg")?;
        finite_positive(self.radius_m, "radius_m")?;
        finite_positive(self.gravity_m_per_s2, "gravity_m_per_s2")?;
        finite_positive(self.initial_theta_rad, "initial_theta_rad")?;
        finite_positive(self.validity_cutoff_theta_rad, "validity_cutoff_theta_rad")?;
        finite_positive(self.timestep_s, "timestep_s")?;
        nonblank(
            &self.small_angle_oracle_source_id,
            "small_angle_oracle_source_id",
        )?;
        if self.maximum_steps == 0 || self.maximum_steps > MAX_REDUCED_DECAY_STEPS {
            return Err(ReducedDecayError::InvalidInput {
                field: "maximum_steps",
            });
        }
        if self.initial_theta_rad <= self.validity_cutoff_theta_rad {
            return Err(ReducedDecayError::InitialStateOutsideValidity {
                theta_rad: self.initial_theta_rad,
                cutoff_theta_rad: self.validity_cutoff_theta_rad,
            });
        }
        if self.initial_theta_rad > MAX_SMALL_ANGLE_THETA_RAD {
            return Err(ReducedDecayError::InvalidInput {
                field: "initial_theta_rad_small_angle",
            });
        }
        if let Some(channel) = &self.dry_contour {
            channel.validate()?;
        }
        if let Some(channel) = &self.bildsten_boundary_layer {
            channel.validate()?;
        }
        if self.dry_contour.is_none() && self.bildsten_boundary_layer.is_none() {
            return Err(ReducedDecayError::NoActiveChannel);
        }
        Ok(())
    }

    fn energy_slope_j_per_rad(&self) -> f64 {
        1.5 * self.mass_kg * self.gravity_m_per_s2 * self.radius_m
    }

    /// Late-stage precession rate from `Omega^2 sin(theta) = 4g/R` [rad/s].
    pub fn omega_rad_s(&self, theta_rad: f64) -> Result<f64, ReducedDecayError> {
        finite_positive(theta_rad, "theta_rad")?;
        let value = (4.0 * self.gravity_m_per_s2 / (self.radius_m * theta_rad.sin())).sqrt();
        if !(value.is_finite() && value > 0.0) {
            return Err(ReducedDecayError::NonFiniteDerived {
                field: "omega_rad_s",
            });
        }
        Ok(value)
    }

    /// Approximate late-stage energy `3 m g R theta / 2` [J].
    pub fn energy_j(&self, theta_rad: f64) -> Result<f64, ReducedDecayError> {
        finite_positive(theta_rad, "theta_rad")?;
        let value = self.energy_slope_j_per_rad() * theta_rad;
        if !(value.is_finite() && value > 0.0) {
            return Err(ReducedDecayError::NonFiniteDerived { field: "energy_j" });
        }
        Ok(value)
    }
}

/// Literature-calibrated analytical benchmark for the Thorne et al. 2026
/// 445 g filleted steel disc on glass.
///
/// Geometry and mass are source-bound, `mu = 1e-4` is the paper's fitted
/// rolling coefficient, and the Bildsten multiplier is exactly one. The
/// initial angle, cutoff, timestep, and step ceiling are explicit FrankenSim
/// rendering controls rather than digitized raw observations. This benchmark
/// does not claim access to the paper's raw trajectory corpus or validation of
/// the boundary-layer prefactor by a full fluid-structure solve.
#[derive(Debug, Clone, PartialEq)]
pub struct Thorne2026SteelGlassBenchmark {
    /// Reported filleted steel specimen.
    pub specimen: LiteratureDiscSpecimen,
    /// Small-angle model and rendering integration controls.
    pub decay: ReducedDecayInput,
    /// Direct source-bound rolling-power channel.
    pub published_rolling: PublishedRollingPowerChannel,
}

impl Thorne2026SteelGlassBenchmark {
    /// Construct the ambient-pressure Figure 3 analytical benchmark.
    pub fn ambient() -> Result<Self, ReducedDecayError> {
        Self::with_air_density(THORNE_2026_AMBIENT_AIR_DENSITY_KG_PER_M3)
    }

    /// Construct the paper's 0.1-atmosphere analytical comparison.
    pub fn partial_vacuum() -> Result<Self, ReducedDecayError> {
        Self::with_air_density(THORNE_2026_VACUUM_AIR_DENSITY_KG_PER_M3)
    }

    /// Construct the same source-bound specimen and model at a declared air
    /// density. This supports the paper's ambient/vacuum ablation without
    /// relabeling the supplied density as a measured chamber record.
    pub fn with_air_density(density_kg_per_m3: f64) -> Result<Self, ReducedDecayError> {
        finite_positive(density_kg_per_m3, "thorne_2026.density_kg_per_m3")?;
        let specimen = LiteratureDiscSpecimen {
            source_id: THORNE_2026_SOURCE_ID.to_owned(),
            diameter_m: THORNE_2026_STEEL_DISC_DIAMETER_M,
            thickness_m: THORNE_2026_STEEL_DISC_THICKNESS_M,
            mass_kg: THORNE_2026_STEEL_DISC_MASS_KG,
            outer_fillet_radius_m: THORNE_2026_STEEL_DISC_FILLET_RADIUS_M,
        };
        let decay = ReducedDecayInput {
            mass_kg: THORNE_2026_STEEL_DISC_MASS_KG,
            radius_m: 0.5 * THORNE_2026_STEEL_DISC_DIAMETER_M,
            gravity_m_per_s2: STANDARD_GRAVITY_M_PER_S2,
            initial_theta_rad: THORNE_2026_RENDER_INITIAL_THETA_RAD,
            validity_cutoff_theta_rad: THORNE_2026_RENDER_VALIDITY_CUTOFF_THETA_RAD,
            timestep_s: THORNE_2026_RENDER_TIMESTEP_S,
            maximum_steps: THORNE_2026_RENDER_MAXIMUM_STEPS,
            small_angle_oracle_source_id: THORNE_2026_SOURCE_ID.to_owned(),
            dry_contour: None,
            bildsten_boundary_layer: Some(BildstenBoundaryLayerChannel {
                source_id: THORNE_2026_SOURCE_ID.to_owned(),
                density_kg_per_m3,
                dynamic_viscosity_pa_s: THORNE_2026_DECLARED_AIR_VISCOSITY_PA_S,
                dimensionless_prefactor: 1.0,
            }),
        };
        let benchmark = Self {
            specimen,
            decay,
            published_rolling: PublishedRollingPowerChannel {
                source_id: THORNE_2026_SOURCE_ID.to_owned(),
                coefficient_mu: THORNE_2026_FITTED_ROLLING_COEFFICIENT,
            },
        };
        benchmark.validate()?;
        Ok(benchmark)
    }

    /// Resolve the exact declared filleted profile at its reported mass.
    pub fn resolve_specimen(&self, cx: &Cx<'_>) -> Result<ResolvedDiscProfile, ReducedDecayError> {
        self.validate()?;
        self.specimen.resolve_mass_matched_profile(cx)
    }

    fn validate(&self) -> Result<(), ReducedDecayError> {
        self.specimen.validate()?;
        self.decay.validate()?;
        self.published_rolling.validate()?;
        if self.specimen.source_id != THORNE_2026_SOURCE_ID
            || self.decay.small_angle_oracle_source_id != THORNE_2026_SOURCE_ID
            || self.published_rolling.source_id != THORNE_2026_SOURCE_ID
        {
            return Err(ReducedDecayError::InvalidInput {
                field: "thorne_2026.source_binding",
            });
        }
        for (field, actual, expected) in [
            (
                "thorne_2026.specimen.diameter_m",
                self.specimen.diameter_m,
                THORNE_2026_STEEL_DISC_DIAMETER_M,
            ),
            (
                "thorne_2026.specimen.thickness_m",
                self.specimen.thickness_m,
                THORNE_2026_STEEL_DISC_THICKNESS_M,
            ),
            (
                "thorne_2026.specimen.mass_kg",
                self.specimen.mass_kg,
                THORNE_2026_STEEL_DISC_MASS_KG,
            ),
            (
                "thorne_2026.specimen.outer_fillet_radius_m",
                self.specimen.outer_fillet_radius_m,
                THORNE_2026_STEEL_DISC_FILLET_RADIUS_M,
            ),
            (
                "thorne_2026.decay.mass_kg",
                self.decay.mass_kg,
                THORNE_2026_STEEL_DISC_MASS_KG,
            ),
            (
                "thorne_2026.decay.radius_m",
                self.decay.radius_m,
                0.5 * THORNE_2026_STEEL_DISC_DIAMETER_M,
            ),
            (
                "thorne_2026.rolling.coefficient_mu",
                self.published_rolling.coefficient_mu,
                THORNE_2026_FITTED_ROLLING_COEFFICIENT,
            ),
        ] {
            if actual.to_bits() != expected.to_bits() {
                return Err(ReducedDecayError::InvalidInput { field });
            }
        }
        let bildsten =
            self.decay
                .bildsten_boundary_layer
                .as_ref()
                .ok_or(ReducedDecayError::InvalidInput {
                    field: "thorne_2026.bildsten_boundary_layer",
                })?;
        if bildsten.source_id != THORNE_2026_SOURCE_ID
            || bildsten.dynamic_viscosity_pa_s.to_bits()
                != THORNE_2026_DECLARED_AIR_VISCOSITY_PA_S.to_bits()
            || bildsten.dimensionless_prefactor.to_bits() != 1.0_f64.to_bits()
            || self.decay.dry_contour.is_some()
        {
            return Err(ReducedDecayError::InvalidInput {
                field: "thorne_2026.channel_binding",
            });
        }
        Ok(())
    }
}

/// Separately retained powers at one state.  `bildsten_boundary_layer_w` is an
/// energy-only closure, and no `fs-flux` exterior-wrench channel is inserted
/// without that distinct generic receipt/dependency being cleanly registered.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ChannelPowers {
    /// Dry contour dissipation [W].
    pub dry_contour_w: f64,
    /// Direct source-bound published rolling dissipation [W].
    pub published_rolling_w: f64,
    /// Bildsten boundary-layer energy-only dissipation [W].
    pub bildsten_boundary_layer_w: f64,
}

impl ChannelPowers {
    /// Sum [W].
    #[must_use]
    pub const fn total_w(self) -> f64 {
        self.dry_contour_w + self.published_rolling_w + self.bildsten_boundary_layer_w
    }
}

/// Per-channel accumulated work [J].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ChannelWork {
    /// Dry contour work [J].
    pub dry_contour_j: f64,
    /// Direct source-bound published rolling work [J].
    pub published_rolling_j: f64,
    /// Bildsten boundary-layer energy-only work [J].
    pub bildsten_boundary_layer_j: f64,
}

impl ChannelWork {
    /// Sum [J].
    #[must_use]
    pub const fn total_j(self) -> f64 {
        self.dry_contour_j + self.published_rolling_j + self.bildsten_boundary_layer_j
    }
}

/// Exact scalar inputs retained with every admitted reduced-decay run.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ReducedDecayParameters {
    /// Disc mass [kg].
    pub mass_kg: f64,
    /// Disc radius [m].
    pub radius_m: f64,
    /// Gravitational acceleration [m/s^2].
    pub gravity_m_per_s2: f64,
    /// Initial inclination [rad].
    pub initial_theta_rad: f64,
    /// Positive terminal validity cutoff [rad].
    pub validity_cutoff_theta_rad: f64,
    /// Declared fixed integration step [s].
    pub timestep_s: f64,
    /// Declared retained-step ceiling.
    pub maximum_steps: u32,
}

impl From<&ReducedDecayInput> for ReducedDecayParameters {
    fn from(input: &ReducedDecayInput) -> Self {
        Self {
            mass_kg: input.mass_kg,
            radius_m: input.radius_m,
            gravity_m_per_s2: input.gravity_m_per_s2,
            initial_theta_rad: input.initial_theta_rad,
            validity_cutoff_theta_rad: input.validity_cutoff_theta_rad,
            timestep_s: input.timestep_s,
            maximum_steps: input.maximum_steps,
        }
    }
}

/// One deterministic retained point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ReducedDecaySample {
    /// Elapsed integration time [s].
    pub time_s: f64,
    /// Inclination [rad].
    pub theta_rad: f64,
    /// Derived precession rate [rad/s].
    pub omega_rad_s: f64,
    /// Reduced energy [J].
    pub energy_j: f64,
    /// Active channel powers [W].
    pub powers: ChannelPowers,
    /// Cumulative channel work [J].
    pub work: ChannelWork,
}

/// Why the integration stopped.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReducedDecayTerminal {
    /// The positive validity cutoff was reached; no theta-zero claim is made.
    ValidityCutoff,
    /// The declared bounded work budget was exhausted before the cutoff.
    StepBudgetExhausted,
}

/// Complete numerical-reference output.
#[derive(Debug, Clone, PartialEq)]
pub struct ReducedDecayRun {
    /// Model, source, interface, authority, and no-validation declarations.
    pub provenance: ReducedDecayProvenance,
    /// Exact scalar model and integration inputs.
    pub parameters: ReducedDecayParameters,
    /// Retained deterministic state samples, including initial/final states.
    pub samples: Vec<ReducedDecaySample>,
    /// Structured terminal condition.
    pub terminal: ReducedDecayTerminal,
    /// `initial_energy - final_energy - channel_work` [J].
    pub energy_closure_residual_j: f64,
}

/// Provenance retained with a reduced numerical-reference run.
///
/// These declarations bind the model to caller inputs and explicitly state
/// that this run is not physical validation or an admitted experiment.
#[derive(Debug, Clone, PartialEq)]
pub struct ReducedDecayProvenance {
    /// Stable reduced-model identity.
    pub model_id: &'static str,
    /// Identity of the small-angle energy/precession oracle.
    pub small_angle_oracle_source_id: String,
    /// Model authority ceiling.
    pub model_authority: &'static str,
    /// Physical-validation disposition.
    pub physical_validation: &'static str,
    /// Cancellation capability deliberately absent from this scalar reference.
    pub cancellation_capability: &'static str,
    /// Dry interface system identity, when the channel is active.
    pub dry_interface_system_id: Option<String>,
    /// Dry caller-source identity, when the channel is active.
    pub dry_source_id: Option<String>,
    /// Dry caller authority, when the channel is active.
    pub dry_authority: Option<InputAuthority>,
    /// Bildsten closure source identity, when the channel is active.
    pub bildsten_source_id: Option<String>,
    /// Authority ceiling of the Bildsten multiplier.
    pub bildsten_multiplier_authority: &'static str,
    /// Source identity for a direct published rolling closure, when active.
    pub published_rolling_source_id: Option<String>,
    /// Exact direct rolling coefficient, when active.
    pub published_rolling_coefficient_mu: Option<f64>,
    /// Authority ceiling of the direct rolling coefficient.
    pub published_rolling_coefficient_authority: &'static str,
    /// Exact boundary-layer gas density [kg/m^3], when active.
    pub bildsten_density_kg_per_m3: Option<f64>,
    /// Exact boundary-layer dynamic viscosity [Pa s], when active.
    pub bildsten_dynamic_viscosity_pa_s: Option<f64>,
    /// Exact boundary-layer dimensionless multiplier, when active.
    pub bildsten_dimensionless_prefactor: Option<f64>,
    /// Source-bound specimen declaration, when the run is a literature benchmark.
    pub literature_specimen: Option<LiteratureDiscSpecimen>,
}

impl ReducedDecayRun {
    /// Final retained sample, refusing externally constructed malformed runs.
    pub fn final_sample(&self) -> Result<&ReducedDecaySample, ReducedDecayError> {
        self.samples.last().ok_or(ReducedDecayError::MissingSample)
    }
}

/// Coarse/fine deterministic comparison.  This is numerical refinement
/// evidence only, not physical validation.
#[derive(Debug, Clone, PartialEq)]
pub struct RefinementEvidence {
    /// Requested-step run.
    pub coarse: ReducedDecayRun,
    /// Half-step run with doubled bounded step budget.
    pub fine: ReducedDecayRun,
    /// Final elapsed-time difference [s].
    pub terminal_time_difference_s: f64,
    /// Final cumulative-work difference [J].
    pub total_work_difference_j: f64,
}

/// Why an encoded-model crossover cannot be compared.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChannelCrossoverNotComparable {
    /// No dry contour channel was declared.
    MissingDryContour,
    /// No Bildsten energy-only channel was declared.
    MissingBildstenBoundaryLayer,
}

/// Deterministic dry/Bildsten channel-crossover diagnostic.
///
/// This is derived solely from the encoded reduced power functions. It is not
/// independent mechanism evidence, model selection, or physical validation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ChannelCrossoverDiagnostic {
    /// The channel pair was not both declared.
    NotComparable {
        /// Typed reason no pairwise comparison exists.
        reason: ChannelCrossoverNotComparable,
    },
    /// Both channels were declared but their power difference did not change
    /// sign on the closed initial-to-cutoff inclination interval.
    NoneWithinInterval,
    /// The encoded power functions cross at this inclination [rad].
    AtInclination {
        /// Bisection-derived inclination [rad].
        theta_rad: f64,
    },
}

#[derive(Clone, Copy)]
enum ReducedRunAuthority<'a> {
    NumericalReference,
    Thorne2026(&'a LiteratureDiscSpecimen),
}

/// Computes the deterministic encoded-model dry/Bildsten crossover, if any.
pub fn channel_crossover_diagnostic(
    input: &ReducedDecayInput,
) -> Result<ChannelCrossoverDiagnostic, ReducedDecayError> {
    input.validate()?;
    if input.dry_contour.is_none() {
        return Ok(ChannelCrossoverDiagnostic::NotComparable {
            reason: ChannelCrossoverNotComparable::MissingDryContour,
        });
    }
    if input.bildsten_boundary_layer.is_none() {
        return Ok(ChannelCrossoverDiagnostic::NotComparable {
            reason: ChannelCrossoverNotComparable::MissingBildstenBoundaryLayer,
        });
    }
    crossover_from_difference(input, |theta_rad| {
        channel_power_difference_w(input, theta_rad)
    })
}

/// Computes the direct published-rolling/Bildsten crossover encoded by the
/// Thorne et al. benchmark.
///
/// This diagnostic compares only the two declared analytical power functions;
/// agreement near the paper's reported crossover is not raw-trajectory or
/// full fluid-structure validation.
pub fn thorne_2026_channel_crossover_diagnostic(
    benchmark: &Thorne2026SteelGlassBenchmark,
) -> Result<ChannelCrossoverDiagnostic, ReducedDecayError> {
    benchmark.validate()?;
    crossover_from_difference(&benchmark.decay, |theta_rad| {
        published_channel_power_difference_w(benchmark, theta_rad)
    })
}

fn crossover_from_difference(
    input: &ReducedDecayInput,
    mut difference: impl FnMut(f64) -> Result<f64, ReducedDecayError>,
) -> Result<ChannelCrossoverDiagnostic, ReducedDecayError> {
    let mut lower_theta_rad = input.validity_cutoff_theta_rad;
    let mut upper_theta_rad = input.initial_theta_rad;
    let mut lower_difference_w = difference(lower_theta_rad)?;
    let upper_difference_w = difference(upper_theta_rad)?;
    if lower_difference_w == 0.0 {
        return Ok(ChannelCrossoverDiagnostic::AtInclination {
            theta_rad: lower_theta_rad,
        });
    }
    if upper_difference_w == 0.0 {
        return Ok(ChannelCrossoverDiagnostic::AtInclination {
            theta_rad: upper_theta_rad,
        });
    }
    if lower_difference_w.is_sign_negative() == upper_difference_w.is_sign_negative() {
        return Ok(ChannelCrossoverDiagnostic::NoneWithinInterval);
    }
    for _ in 0..CHANNEL_CROSSOVER_BISECTION_STEPS {
        let midpoint_theta_rad = 0.5 * (lower_theta_rad + upper_theta_rad);
        let midpoint_difference_w = difference(midpoint_theta_rad)?;
        if midpoint_difference_w == 0.0 {
            return Ok(ChannelCrossoverDiagnostic::AtInclination {
                theta_rad: midpoint_theta_rad,
            });
        }
        if midpoint_difference_w.is_sign_negative() == lower_difference_w.is_sign_negative() {
            lower_theta_rad = midpoint_theta_rad;
            lower_difference_w = midpoint_difference_w;
        } else {
            upper_theta_rad = midpoint_theta_rad;
        }
    }
    Ok(ChannelCrossoverDiagnostic::AtInclination {
        theta_rad: 0.5 * (lower_theta_rad + upper_theta_rad),
    })
}

/// Produces the stable, single-record numerical-reference runner output.
///
/// The record carries computed channel accounting and refinement differences;
/// it does not imply experimental validation or a resolved flow field.
#[must_use]
pub fn structured_runner_output(
    run: &ReducedDecayRun,
    refinement: &RefinementEvidence,
) -> Result<String, ReducedDecayError> {
    let final_sample = run.final_sample()?;
    Ok(format!(
        "schema=reduced-decay-v1 model_id={} small_angle_oracle_source_id={} model_authority={} physical_validation={} cancellation_capability={} dry_interface_system_id={} dry_source_id={} dry_authority={:?} published_rolling_source_id={} published_rolling_coefficient_authority={} bildsten_source_id={} bildsten_multiplier_authority={} terminal={:?} time_s={:.12e} theta_rad={:.12e} energy_j={:.12e} dry_work_j={:.12e} published_rolling_work_j={:.12e} bildsten_work_j={:.12e} closure_residual_j={:.12e}\nrefinement_terminal_time_difference_s={:.12e} refinement_total_work_difference_j={:.12e} evidence_scope=numerical-self-consistency-only",
        run.provenance.model_id,
        run.provenance.small_angle_oracle_source_id,
        run.provenance.model_authority,
        run.provenance.physical_validation,
        run.provenance.cancellation_capability,
        run.provenance
            .dry_interface_system_id
            .as_deref()
            .unwrap_or("none"),
        run.provenance.dry_source_id.as_deref().unwrap_or("none"),
        run.provenance.dry_authority,
        run.provenance
            .published_rolling_source_id
            .as_deref()
            .unwrap_or("none"),
        run.provenance.published_rolling_coefficient_authority,
        run.provenance
            .bildsten_source_id
            .as_deref()
            .unwrap_or("none"),
        run.provenance.bildsten_multiplier_authority,
        run.terminal,
        final_sample.time_s,
        final_sample.theta_rad,
        final_sample.energy_j,
        final_sample.work.dry_contour_j,
        final_sample.work.published_rolling_j,
        final_sample.work.bildsten_boundary_layer_j,
        run.energy_closure_residual_j,
        refinement.terminal_time_difference_s,
        refinement.total_work_difference_j,
    ))
}

/// Runs the fixed-step reference with exact per-step energy decrement.
pub fn run_reduced_decay(input: &ReducedDecayInput) -> Result<ReducedDecayRun, ReducedDecayError> {
    run_reduced_decay_internal(input, None, ReducedRunAuthority::NumericalReference)
}

/// Runs the source-bound Thorne et al. analytical benchmark.
pub fn run_thorne_2026_steel_glass_benchmark(
    benchmark: &Thorne2026SteelGlassBenchmark,
) -> Result<ReducedDecayRun, ReducedDecayError> {
    benchmark.validate()?;
    run_reduced_decay_internal(
        &benchmark.decay,
        Some(&benchmark.published_rolling),
        ReducedRunAuthority::Thorne2026(&benchmark.specimen),
    )
}

fn run_reduced_decay_internal(
    input: &ReducedDecayInput,
    published_rolling: Option<&PublishedRollingPowerChannel>,
    authority: ReducedRunAuthority<'_>,
) -> Result<ReducedDecayRun, ReducedDecayError> {
    input.validate()?;
    if let Some(channel) = published_rolling {
        channel.validate()?;
    }
    let provenance = provenance(input, published_rolling, authority);
    let parameters = ReducedDecayParameters::from(input);
    let initial_energy_j = input.energy_j(input.initial_theta_rad)?;
    let mut theta_rad = input.initial_theta_rad;
    let mut time_s = 0.0;
    let mut work = ChannelWork {
        dry_contour_j: 0.0,
        published_rolling_j: 0.0,
        bildsten_boundary_layer_j: 0.0,
    };
    let capacity = usize::try_from(input.maximum_steps)
        .map_err(|_| ReducedDecayError::InvalidInput {
            field: "maximum_steps",
        })?
        .saturating_add(1);
    let mut samples = Vec::with_capacity(capacity);
    let initial_powers = powers_at(input, published_rolling, theta_rad)?;
    samples.push(sample(input, time_s, theta_rad, initial_powers, work)?);

    for _ in 0..input.maximum_steps {
        let powers = powers_at(input, published_rolling, theta_rad)?;
        let total_power_w = powers.total_w();
        if !(total_power_w.is_finite() && total_power_w > 0.0) {
            return Err(ReducedDecayError::InvalidInput {
                field: "active_channel_power_w",
            });
        }
        let time_to_cutoff_s = (theta_rad - input.validity_cutoff_theta_rad)
            * input.energy_slope_j_per_rad()
            / total_power_w;
        let dt_s = input.timestep_s.min(time_to_cutoff_s);
        if !(dt_s.is_finite() && dt_s > 0.0) {
            return Err(ReducedDecayError::NonFiniteDerived {
                field: "integration_step_s",
            });
        }
        work.dry_contour_j += powers.dry_contour_w * dt_s;
        work.published_rolling_j += powers.published_rolling_w * dt_s;
        work.bildsten_boundary_layer_j += powers.bildsten_boundary_layer_w * dt_s;
        time_s += dt_s;
        theta_rad -= total_power_w * dt_s / input.energy_slope_j_per_rad();
        if time_to_cutoff_s <= input.timestep_s {
            theta_rad = input.validity_cutoff_theta_rad;
        }
        let next_powers = powers_at(input, published_rolling, theta_rad)?;
        samples.push(sample(input, time_s, theta_rad, next_powers, work)?);
        if theta_rad <= input.validity_cutoff_theta_rad {
            let final_energy_j = input.energy_j(theta_rad)?;
            return Ok(ReducedDecayRun {
                provenance,
                parameters,
                samples,
                terminal: ReducedDecayTerminal::ValidityCutoff,
                energy_closure_residual_j: initial_energy_j - final_energy_j - work.total_j(),
            });
        }
    }
    let final_energy_j = input.energy_j(theta_rad)?;
    Ok(ReducedDecayRun {
        provenance,
        parameters,
        samples,
        terminal: ReducedDecayTerminal::StepBudgetExhausted,
        energy_closure_residual_j: initial_energy_j - final_energy_j - work.total_j(),
    })
}

/// Runs the requested and half-step numerical references.
pub fn refinement_evidence(
    input: &ReducedDecayInput,
) -> Result<RefinementEvidence, ReducedDecayError> {
    refinement_evidence_internal(input, None, ReducedRunAuthority::NumericalReference)
}

/// Runs requested-step and half-step versions of the complete source-bound
/// Thorne benchmark, retaining both its published rolling and Bildsten power
/// channels in each run.
///
/// This is timestep-refinement evidence for the encoded analytical model. It
/// is not raw-trajectory agreement or full fluid-structure validation.
pub fn thorne_2026_refinement_evidence(
    benchmark: &Thorne2026SteelGlassBenchmark,
) -> Result<RefinementEvidence, ReducedDecayError> {
    benchmark.validate()?;
    refinement_evidence_internal(
        &benchmark.decay,
        Some(&benchmark.published_rolling),
        ReducedRunAuthority::Thorne2026(&benchmark.specimen),
    )
}

fn refinement_evidence_internal(
    input: &ReducedDecayInput,
    published_rolling: Option<&PublishedRollingPowerChannel>,
    authority: ReducedRunAuthority<'_>,
) -> Result<RefinementEvidence, ReducedDecayError> {
    input.validate()?;
    if let Some(channel) = published_rolling {
        channel.validate()?;
    }
    let fine_steps = input
        .maximum_steps
        .checked_mul(2)
        .filter(|steps| *steps <= MAX_REDUCED_DECAY_STEPS)
        .ok_or(ReducedDecayError::RefinementStepBudgetOverflow)?;
    let coarse = run_reduced_decay_internal(input, published_rolling, authority)?;
    let mut fine_input = input.clone();
    fine_input.timestep_s *= 0.5;
    fine_input.maximum_steps = fine_steps;
    let fine = run_reduced_decay_internal(&fine_input, published_rolling, authority)?;
    if coarse.terminal != fine.terminal {
        return Err(ReducedDecayError::RefinementTerminalMismatch {
            coarse: coarse.terminal,
            fine: fine.terminal,
        });
    }
    if coarse.terminal != ReducedDecayTerminal::ValidityCutoff {
        return Err(ReducedDecayError::RefinementIncompleteTerminal {
            terminal: coarse.terminal,
        });
    }
    let coarse_terminal = coarse.final_sample()?;
    let fine_terminal = fine.final_sample()?;
    Ok(RefinementEvidence {
        terminal_time_difference_s: (coarse_terminal.time_s - fine_terminal.time_s).abs(),
        total_work_difference_j: (coarse_terminal.work.total_j() - fine_terminal.work.total_j())
            .abs(),
        coarse,
        fine,
    })
}

fn powers_at(
    input: &ReducedDecayInput,
    published_rolling: Option<&PublishedRollingPowerChannel>,
    theta_rad: f64,
) -> Result<ChannelPowers, ReducedDecayError> {
    let omega_rad_s = input.omega_rad_s(theta_rad)?;
    let dry_contour_w = if let Some(channel) = &input.dry_contour {
        let contour_speed_m_per_s = input.radius_m * theta_rad.cos() * omega_rad_s;
        let response = ConstantContourForce {
            force_n: channel.contour_force_n,
        }
        .evaluate(
            &channel.interface,
            ResistanceInput {
                normal_force_n: channel.normal_force_n,
                angular_speed_rad_s: omega_rad_s,
                contour_speed_mps: contour_speed_m_per_s,
            },
        )
        .map_err(|error| ReducedDecayError::DryLawRefusal {
            detail: error.to_string(),
        })?;
        response.dissipated_power_w()
    } else {
        0.0
    };
    let published_rolling_w = if let Some(channel) = published_rolling {
        channel.power_w(
            input.mass_kg,
            input.gravity_m_per_s2,
            input.radius_m,
            omega_rad_s,
        )?
    } else {
        0.0
    };
    let bildsten_boundary_layer_w = if let Some(channel) = &input.bildsten_boundary_layer {
        channel.power_w(input.radius_m, input.gravity_m_per_s2, theta_rad)?
    } else {
        0.0
    };
    Ok(ChannelPowers {
        dry_contour_w,
        published_rolling_w,
        bildsten_boundary_layer_w,
    })
}

fn channel_power_difference_w(
    input: &ReducedDecayInput,
    theta_rad: f64,
) -> Result<f64, ReducedDecayError> {
    let powers = powers_at(input, None, theta_rad)?;
    let difference_w = powers.dry_contour_w - powers.bildsten_boundary_layer_w;
    if !difference_w.is_finite() {
        return Err(ReducedDecayError::NonFiniteDerived {
            field: "channel_crossover_difference_w",
        });
    }
    Ok(difference_w)
}

fn published_channel_power_difference_w(
    benchmark: &Thorne2026SteelGlassBenchmark,
    theta_rad: f64,
) -> Result<f64, ReducedDecayError> {
    let powers = powers_at(
        &benchmark.decay,
        Some(&benchmark.published_rolling),
        theta_rad,
    )?;
    let difference_w = powers.published_rolling_w - powers.bildsten_boundary_layer_w;
    if !difference_w.is_finite() {
        return Err(ReducedDecayError::NonFiniteDerived {
            field: "thorne_2026.channel_crossover_difference_w",
        });
    }
    Ok(difference_w)
}

fn provenance(
    input: &ReducedDecayInput,
    published_rolling: Option<&PublishedRollingPowerChannel>,
    authority: ReducedRunAuthority<'_>,
) -> ReducedDecayProvenance {
    let dry = input.dry_contour.as_ref();
    let (model_authority, physical_validation, bildsten_multiplier_authority, literature_specimen) =
        match authority {
            ReducedRunAuthority::NumericalReference => (
                "numerical-reference-only",
                "not-claimed",
                "caller-declared",
                None,
            ),
            ReducedRunAuthority::Thorne2026(specimen) => (
                "literature-calibrated-analytical",
                "no-raw-trajectory-or-full-fsi-validation-claimed",
                "source-bound-published-coefficient",
                Some(specimen.clone()),
            ),
        };
    ReducedDecayProvenance {
        model_id: REDUCED_DECAY_MODEL_ID,
        small_angle_oracle_source_id: input.small_angle_oracle_source_id.clone(),
        model_authority,
        physical_validation,
        cancellation_capability: "not-implemented",
        dry_interface_system_id: dry
            .map(|channel| channel.interface.ordered_system_id().to_owned()),
        dry_source_id: dry.map(|channel| channel.interface.provenance().source_id().to_owned()),
        dry_authority: dry.map(|channel| channel.interface.provenance().authority()),
        bildsten_source_id: input
            .bildsten_boundary_layer
            .as_ref()
            .map(|channel| channel.source_id.clone()),
        bildsten_multiplier_authority,
        published_rolling_source_id: published_rolling.map(|channel| channel.source_id.clone()),
        published_rolling_coefficient_mu: published_rolling.map(|channel| channel.coefficient_mu),
        published_rolling_coefficient_authority: if published_rolling.is_some() {
            "source-bound-fitted-coefficient"
        } else {
            "not-present"
        },
        bildsten_density_kg_per_m3: input
            .bildsten_boundary_layer
            .as_ref()
            .map(|channel| channel.density_kg_per_m3),
        bildsten_dynamic_viscosity_pa_s: input
            .bildsten_boundary_layer
            .as_ref()
            .map(|channel| channel.dynamic_viscosity_pa_s),
        bildsten_dimensionless_prefactor: input
            .bildsten_boundary_layer
            .as_ref()
            .map(|channel| channel.dimensionless_prefactor),
        literature_specimen,
    }
}

fn sample(
    input: &ReducedDecayInput,
    time_s: f64,
    theta_rad: f64,
    powers: ChannelPowers,
    work: ChannelWork,
) -> Result<ReducedDecaySample, ReducedDecayError> {
    Ok(ReducedDecaySample {
        time_s,
        theta_rad,
        omega_rad_s: input.omega_rad_s(theta_rad)?,
        energy_j: input.energy_j(theta_rad)?,
        powers,
        work,
    })
}
