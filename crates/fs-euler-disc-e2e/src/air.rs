//! Explicit sector/strip adapter from a tilted rigid disc to `fs_flux::gas_film`.
//!
//! Each azimuthal sector is an independent radial strip of constant equivalent
//! width.  This preserves the sampled area and first radial moment, but makes
//! **no cross-sector-flow, resolved rim, molecular, contact, or outcome claim**.

use core::fmt;

use fs_flux::{
    ContactExclusionMask, GasFilmApplicability, GasFilmBoundaryTopology, GasFilmBudget,
    GasFilmCheckpoint, GasFilmError, GasFilmGrid1d, GasFilmIdentity, GasFilmInput,
    GasFilmInputAuthority, GasFilmReceipt, GasFilmUncertainty, IsothermalIdealGas, MovingWallInput,
    RoughnessPolicy, SlipPolicy, isothermal_compressible_reynolds_model_id,
    solve_isothermal_gas_film_1d,
};

/// Stable identity for this declared independent-sector approximation.
pub const TILTED_DISC_GAS_FILM_ADAPTER_ID: &str = "euler-disc-tilted-sector-gas-film-v1";
const MAX_SECTORS: usize = 128;
const MAX_RADIAL_CELLS: usize = 128;

/// Small explicit world-frame vector used by this adapter.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AirVec3 {
    /// World x component.
    pub x: f64,
    /// World y component.
    pub y: f64,
    /// World z component.
    pub z: f64,
}

impl AirVec3 {
    /// The zero vector.
    pub const ZERO: Self = Self::new(0.0, 0.0, 0.0);

    /// Construct a vector.
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    fn add(self, other: Self) -> Self {
        Self::new(self.x + other.x, self.y + other.y, self.z + other.z)
    }

    fn scale(self, scalar: f64) -> Self {
        Self::new(self.x * scalar, self.y * scalar, self.z * scalar)
    }

    fn dot(self, other: Self) -> f64 {
        self.x
            .mul_add(other.x, self.y.mul_add(other.y, self.z * other.z))
    }

    fn cross(self, other: Self) -> Self {
        Self::new(
            self.y.mul_add(other.z, -(self.z * other.y)),
            self.z.mul_add(other.x, -(self.x * other.z)),
            self.x.mul_add(other.y, -(self.y * other.x)),
        )
    }

    fn norm(self) -> Option<f64> {
        if !self.x.is_finite() || !self.y.is_finite() || !self.z.is_finite() {
            return None;
        }
        let scale = self.x.abs().max(self.y.abs()).max(self.z.abs());
        if scale == 0.0 {
            return Some(0.0);
        }
        let reduced = self.scale(1.0 / scale);
        let norm = scale * reduced.dot(reduced).sqrt();
        norm.is_finite().then_some(norm)
    }

    fn unit(self, field: &'static str) -> Result<Self, AirFilmError> {
        let norm = self.norm().ok_or(AirFilmError::InvalidInput { field })?;
        if norm == 0.0 {
            return Err(AirFilmError::InvalidInput { field });
        }
        let unit = self.scale(1.0 / norm);
        if !unit.x.is_finite() || !unit.y.is_finite() || !unit.z.is_finite() {
            return Err(AirFilmError::NonFiniteDerived { field });
        }
        Ok(unit)
    }
}

/// Typed refusal from the Euler-disc gas-film adapter.
#[derive(Debug, Clone, PartialEq)]
pub enum AirFilmError {
    /// A public scalar/vector/count is invalid.
    InvalidInput { field: &'static str },
    /// A derived value overflowed or became non-finite.
    NonFiniteDerived { field: &'static str },
    /// A requested regime or topology is intentionally unavailable.
    Unavailable { reason: &'static str },
    /// Contact leaves a non-prefix radial gas domain that `fs_flux` cannot represent.
    ContactTopologyUnavailable,
    /// Restart state was not made by the same identity/discretization.
    CheckpointMismatch { field: &'static str },
    /// The delegated generic film operator refused the sector.
    GasFilmRefusal { detail: GasFilmError },
}

impl fmt::Display for AirFilmError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for AirFilmError {}

/// Provenance that remains caller-declared/synthetic and never becomes an admission.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AirFilmIdentity {
    /// Stable caller case identity.
    pub case_id: String,
    /// Adapter identity; must equal [`TILTED_DISC_GAS_FILM_ADAPTER_ID`].
    pub adapter_model_id: String,
    /// World-frame identity.
    pub frame_id: String,
    /// Named prescribed-base source.
    pub base_motion_id: String,
    /// Caller canonical gas-species identity.
    pub gas_species_id: String,
    /// Caller canonical EOS identity.
    pub eos_id: String,
    /// Caller canonical viscosity source identity.
    pub viscosity_source_id: String,
    /// Caller canonical isothermal-model identity.
    pub thermal_model_id: String,
    /// Caller-owned identity for frozen gas/boundary/applicability inputs.
    pub configuration_id: String,
    /// Deterministic replay seed.
    pub deterministic_seed: u64,
    /// Authority carried into every sector request.
    pub authority: GasFilmInputAuthority,
}

/// Rigid disc pose and prescribed rigid-body rates in one world frame.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TiltedDiscKinematics {
    /// Disc centre/COM [m].
    pub center_world_m: AirVec3,
    /// Unit normal pointing away from the base-facing disc surface.
    pub normal_away_from_base_world: AirVec3,
    /// COM translational velocity [m s^-1].
    pub center_velocity_world_m_per_s: AirVec3,
    /// Angular velocity [rad s^-1].
    pub angular_velocity_world_rad_per_s: AirVec3,
}

/// Prescribed horizontal plane base.  Flexible feedback is intentionally absent.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PrescribedPlaneBase {
    /// Plane elevation in the world z direction [m].
    pub height_m: f64,
    /// Uniform vertical plane velocity [m s^-1].
    pub vertical_velocity_m_per_s: f64,
}

/// Explicit sector/radial-grid approximation controls.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AirFilmDiscretization {
    /// Independent azimuthal sector count.
    pub azimuthal_sectors: usize,
    /// Uniform radial cells per sector.
    pub radial_cells: usize,
}

/// Static gas/contact handoff threshold.  It is not a numerical gap floor.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ContactExclusion {
    /// A sampled gap at or below this physical handoff value is excluded [m].
    pub handoff_gap_m: f64,
}

/// Full prescribed-base one-step request.
#[derive(Debug, Clone, PartialEq)]
pub struct TiltedDiscAirFilmInput {
    /// Identities, seed, and authority.
    pub identity: AirFilmIdentity,
    /// Disc radius [m].
    pub disc_radius_m: f64,
    /// Positive disc half thickness [m]. The adapter samples only the flat
    /// base-facing face; chamfers and rim geometry are unavailable.
    pub disc_half_thickness_m: f64,
    /// Rigid disc pose/rates.
    pub disc: TiltedDiscKinematics,
    /// Prescribed plane base.
    pub base: PrescribedPlaneBase,
    /// Independent sector/radial grid.
    pub discretization: AirFilmDiscretization,
    /// Explicit contact handoff, never a floor.
    pub contact_exclusion: ContactExclusion,
    /// Gas properties passed to every sector.
    pub gas: IsothermalIdealGas,
    /// Boundary connectivity for each independent radial strip.
    pub boundary: GasFilmBoundaryTopology,
    /// Continuum/no-slip gate.
    pub slip_policy: SlipPolicy,
    /// Smooth-wall gate.
    pub roughness_policy: RoughnessPolicy,
    /// Continuum/slope/Mach envelope.
    pub applicability: GasFilmApplicability,
    /// Caller-declared input uncertainty, retained but not propagated.
    pub uncertainty: GasFilmUncertainty,
    /// Fresh-state pressure [Pa absolute].
    pub initial_absolute_pressure_pa: f64,
    /// Ambient/reference pressure used only for the body pressure wrench [Pa absolute].
    pub gauge_reference_absolute_pressure_pa: f64,
    /// Step duration [s].
    pub timestep_s: f64,
    /// Bounded nonlinear iteration controls.
    pub budget: GasFilmBudget,
}

/// One resolved radial/azimuthal sample used for geometry and wrench quadrature.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AirFilmSample {
    /// Sector index.
    pub sector_index: usize,
    /// Radial-cell index.
    pub radial_index: usize,
    /// Centre-to-sample vector [m].
    pub lever_arm_world_m: AirVec3,
    /// Exact sampled separation to the prescribed plane [m].
    pub gap_m: f64,
    /// Exact material-point normal gap rate [m s^-1].
    pub gap_rate_m_per_s: f64,
    /// Relative radial wall speed [m s^-1].
    pub radial_relative_velocity_m_per_s: f64,
    /// Relative circumferential wall speed [m s^-1].
    pub circumferential_relative_velocity_m_per_s: f64,
    /// Whether the sample was handed to contact rather than gas.
    pub excluded_by_contact: bool,
}

/// Pressure/shear resultant about the disc COM.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AirFilmWrench {
    /// Net gas-on-body force using `p - gauge_reference` [N].
    pub force_world_n: AirVec3,
    /// Net gas moment about the disc COM [N m].
    pub moment_about_com_world_n_m: AirVec3,
}

/// Aggregated conservative accounting across independent sectors.
#[derive(Debug, Clone, PartialEq)]
pub struct AirFilmReceipt {
    /// Gas force/moment on the body about COM.
    pub wrench: AirFilmWrench,
    /// Gas mass [kg].
    pub gas_mass_kg: f64,
    /// Storage rate [kg s^-1].
    pub storage_rate_kg_per_s: f64,
    /// Signed outward boundary/vent mass flux [kg s^-1].
    pub outward_mass_flux_kg_per_s: f64,
    /// Signed outward boundary/vent enthalpy transport [W].
    pub outward_enthalpy_flux_w: f64,
    /// Conservative mass closure residual [kg s^-1].
    pub mass_closure_residual_kg_per_s: f64,
    /// Mechanical work rate from both walls into gas [W].
    pub wall_power_to_gas_w: f64,
    /// Explicit isothermal viscous heat receipt [W]; not a total energy closure.
    pub viscous_heat_w: f64,
    /// Circumferential no-cross-sector Couette contribution to wall power [W].
    pub circumferential_couette_power_w: f64,
    /// Circumferential no-cross-sector Couette contribution to heat [W].
    pub circumferential_couette_heat_w: f64,
    /// Equivalent rectangular strip width per sector [m].
    pub equal_area_strip_width_m: f64,
    /// `surrogate - exact` integral of radial lever over area [m^3]. It is
    /// exactly negative one quarter of the exact polar-disc first moment.
    pub signed_first_radial_moment_discrepancy_m3: f64,
    /// Largest delegated sector residual [kg m^-2 s^-1].
    pub max_sector_mass_residual_kg_m2_s: f64,
    /// Declared reference used for the body pressure wrench [Pa absolute].
    pub gauge_reference_absolute_pressure_pa: f64,
    /// Retained authority, never promoted.
    pub input_authority: GasFilmInputAuthority,
}

/// Restartable state, bound to this adapter identity and sector count.
#[derive(Debug, Clone, PartialEq)]
pub struct AirFilmCheckpoint {
    /// Adapter identity.
    pub adapter_model_id: String,
    /// Caller case identity.
    pub case_id: String,
    /// Caller-owned frozen configuration identity.
    pub configuration_id: String,
    /// Sector count bound into this state.
    pub azimuthal_sectors: usize,
    /// Ordered sector checkpoints.
    pub sectors: Vec<GasFilmCheckpoint>,
}

/// Accepted adapter result. It does not validate an Euler outcome or air dominance.
#[derive(Debug, Clone, PartialEq)]
pub struct AirFilmStep {
    /// Exact geometry/rate samples used by the approximation.
    pub samples: Vec<AirFilmSample>,
    /// Aggregated body wrench and conservative receipts.
    pub receipt: AirFilmReceipt,
    /// Identity-bound restart state.
    pub checkpoint: AirFilmCheckpoint,
}

fn finite(value: f64, field: &'static str) -> Result<(), AirFilmError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(AirFilmError::InvalidInput { field })
    }
}

fn positive(value: f64, field: &'static str) -> Result<(), AirFilmError> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(AirFilmError::InvalidInput { field })
    }
}

fn checked(value: f64, field: &'static str) -> Result<f64, AirFilmError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(AirFilmError::NonFiniteDerived { field })
    }
}

fn checked_vec(value: AirVec3, field: &'static str) -> Result<AirVec3, AirFilmError> {
    if value.norm().is_some() {
        Ok(value)
    } else {
        Err(AirFilmError::NonFiniteDerived { field })
    }
}

fn canonical(value: &str, field: &'static str) -> Result<(), AirFilmError> {
    if value.is_empty()
        || !value.bytes().all(|byte| {
            byte.is_ascii_alphanumeric() || matches!(byte, b'-' | b'_' | b'.' | b':' | b'/')
        })
    {
        return Err(AirFilmError::InvalidInput { field });
    }
    Ok(())
}

fn sector_basis(normal: AirVec3, angle_rad: f64) -> Result<(AirVec3, AirVec3), AirFilmError> {
    let anchor = if normal.z.abs() < 0.9 {
        AirVec3::new(0.0, 0.0, 1.0)
    } else {
        AirVec3::new(1.0, 0.0, 0.0)
    };
    let u = anchor
        .cross(normal)
        .unit("disc.normal_away_from_base_world")?;
    let v = normal.cross(u).unit("disc.normal_away_from_base_world")?;
    Ok((
        u.scale(angle_rad.cos()).add(v.scale(angle_rad.sin())),
        normal
            .cross(u.scale(angle_rad.cos()).add(v.scale(angle_rad.sin())))
            .unit("sector_tangent")?,
    ))
}

impl TiltedDiscAirFilmInput {
    fn validate(&self) -> Result<AirVec3, AirFilmError> {
        for (field, value) in [
            ("identity.case_id", self.identity.case_id.as_str()),
            (
                "identity.adapter_model_id",
                self.identity.adapter_model_id.as_str(),
            ),
            ("identity.frame_id", self.identity.frame_id.as_str()),
            (
                "identity.base_motion_id",
                self.identity.base_motion_id.as_str(),
            ),
            (
                "identity.gas_species_id",
                self.identity.gas_species_id.as_str(),
            ),
            ("identity.eos_id", self.identity.eos_id.as_str()),
            (
                "identity.viscosity_source_id",
                self.identity.viscosity_source_id.as_str(),
            ),
            (
                "identity.thermal_model_id",
                self.identity.thermal_model_id.as_str(),
            ),
            (
                "identity.configuration_id",
                self.identity.configuration_id.as_str(),
            ),
        ] {
            canonical(value, field)?;
        }
        if self.identity.adapter_model_id != TILTED_DISC_GAS_FILM_ADAPTER_ID {
            return Err(AirFilmError::Unavailable {
                reason: "unrecognized-air-film-adapter",
            });
        }
        positive(self.disc_radius_m, "disc_radius_m")?;
        positive(self.disc_half_thickness_m, "disc_half_thickness_m")?;
        positive(
            self.contact_exclusion.handoff_gap_m,
            "contact_exclusion.handoff_gap_m",
        )?;
        positive(
            self.gauge_reference_absolute_pressure_pa,
            "gauge_reference_absolute_pressure_pa",
        )?;
        finite(self.base.height_m, "base.height_m")?;
        finite(
            self.base.vertical_velocity_m_per_s,
            "base.vertical_velocity_m_per_s",
        )?;
        let normal = self
            .disc
            .normal_away_from_base_world
            .unit("disc.normal_away_from_base_world")?;
        if normal.z <= 0.0 {
            return Err(AirFilmError::Unavailable {
                reason: "disc-normal-must-point-away-from-horizontal-base",
            });
        }
        for (field, value) in [
            ("disc.center_world_m", self.disc.center_world_m),
            (
                "disc.center_velocity_world_m_per_s",
                self.disc.center_velocity_world_m_per_s,
            ),
            (
                "disc.angular_velocity_world_rad_per_s",
                self.disc.angular_velocity_world_rad_per_s,
            ),
        ] {
            if value.norm().is_none() {
                return Err(AirFilmError::InvalidInput { field });
            }
        }
        if !(2..=MAX_SECTORS).contains(&self.discretization.azimuthal_sectors) {
            return Err(AirFilmError::InvalidInput {
                field: "discretization.azimuthal_sectors",
            });
        }
        if !(2..=MAX_RADIAL_CELLS).contains(&self.discretization.radial_cells) {
            return Err(AirFilmError::InvalidInput {
                field: "discretization.radial_cells",
            });
        }
        Ok(normal)
    }
}

/// Evaluate the exact planar-disc separation and material-point velocity at a polar sample.
pub fn sample_tilted_disc_gap(
    input: &TiltedDiscAirFilmInput,
    sector_index: usize,
    radial_index: usize,
) -> Result<AirFilmSample, AirFilmError> {
    let normal = input.validate()?;
    if sector_index >= input.discretization.azimuthal_sectors
        || radial_index >= input.discretization.radial_cells
    {
        return Err(AirFilmError::InvalidInput {
            field: "sample_index",
        });
    }
    let sectors = input.discretization.azimuthal_sectors as f64;
    let radial_cells = input.discretization.radial_cells as f64;
    let angle = 2.0 * core::f64::consts::PI * (sector_index as f64 + 0.5) / sectors;
    let (radial, circumferential) = sector_basis(normal, angle)?;
    let radius = input.disc_radius_m * (radial_index as f64 + 0.5) / radial_cells;
    let lever = radial
        .scale(radius)
        .add(normal.scale(-input.disc_half_thickness_m));
    let point_velocity = input
        .disc
        .center_velocity_world_m_per_s
        .add(input.disc.angular_velocity_world_rad_per_s.cross(lever));
    let gap = checked(
        input.disc.center_world_m.z + lever.z - input.base.height_m,
        "sample_gap_m",
    )?;
    let gap_rate = checked(
        point_velocity.z - input.base.vertical_velocity_m_per_s,
        "sample_gap_rate_m_per_s",
    )?;
    let radial_speed = checked(point_velocity.dot(radial), "sample_radial_speed_m_per_s")?;
    let circumferential_speed = checked(
        point_velocity.dot(circumferential),
        "sample_circumferential_speed_m_per_s",
    )?;
    Ok(AirFilmSample {
        sector_index,
        radial_index,
        lever_arm_world_m: lever,
        gap_m: gap,
        gap_rate_m_per_s: gap_rate,
        radial_relative_velocity_m_per_s: radial_speed,
        circumferential_relative_velocity_m_per_s: circumferential_speed,
        excluded_by_contact: gap <= input.contact_exclusion.handoff_gap_m,
    })
}

/// Advance independent radial strips and aggregate their body wrench/receipts.
pub fn solve_tilted_disc_air_film(
    input: &TiltedDiscAirFilmInput,
    checkpoint: Option<&AirFilmCheckpoint>,
) -> Result<AirFilmStep, AirFilmError> {
    let normal = input.validate()?;
    if let Some(state) = checkpoint {
        if state.adapter_model_id != input.identity.adapter_model_id {
            return Err(AirFilmError::CheckpointMismatch {
                field: "adapter_model_id",
            });
        }
        if state.case_id != input.identity.case_id {
            return Err(AirFilmError::CheckpointMismatch { field: "case_id" });
        }
        if state.configuration_id != input.identity.configuration_id {
            return Err(AirFilmError::CheckpointMismatch {
                field: "configuration_id",
            });
        }
        if state.azimuthal_sectors != input.discretization.azimuthal_sectors
            || state.sectors.len() != input.discretization.azimuthal_sectors
        {
            return Err(AirFilmError::CheckpointMismatch {
                field: "sector_count",
            });
        }
    }
    let sector_angle = 2.0 * core::f64::consts::PI / input.discretization.azimuthal_sectors as f64;
    // Equal-area rectangular surrogate; it deliberately does not preserve the
    // polar-disc first radial moment.
    let strip_width = 0.5 * input.disc_radius_m * sector_angle;
    let first_moment_discrepancy = checked(
        core::f64::consts::PI * input.disc_radius_m.powi(3) / 2.0
            - 2.0 * core::f64::consts::PI * input.disc_radius_m.powi(3) / 3.0,
        "equal_area_surrogate_first_radial_moment_discrepancy",
    )?;
    let mut samples = Vec::with_capacity(
        input.discretization.azimuthal_sectors * input.discretization.radial_cells,
    );
    let mut force = AirVec3::ZERO;
    let mut moment = AirVec3::ZERO;
    let mut gas_mass = 0.0;
    let mut storage = 0.0;
    let mut mass_flux = 0.0;
    let mut enthalpy = 0.0;
    let mut closure = 0.0;
    let mut wall_power = 0.0;
    let mut heat = 0.0;
    let mut circumferential_power = 0.0;
    let mut max_residual = 0.0;
    let mut sector_checkpoints = Vec::with_capacity(input.discretization.azimuthal_sectors);
    for sector in 0..input.discretization.azimuthal_sectors {
        let first = samples.len();
        for radial in 0..input.discretization.radial_cells {
            samples.push(sample_tilted_disc_gap(input, sector, radial)?);
        }
        let sector_samples = &samples[first..];
        let mut excluded = sector_samples
            .iter()
            .map(|sample| sample.excluded_by_contact)
            .collect::<Vec<_>>();
        if excluded.iter().all(|value| *value) {
            return Err(AirFilmError::Unavailable {
                reason: "contact-excludes-entire-sector",
            });
        }
        let first_excluded = excluded
            .iter()
            .position(|value| *value)
            .unwrap_or(excluded.len());
        if excluded[first_excluded..].iter().any(|value| !*value) {
            return Err(AirFilmError::ContactTopologyUnavailable);
        }
        let active = first_excluded;
        if active < 2 {
            return Err(AirFilmError::Unavailable {
                reason: "fewer-than-two-active-radial-cells",
            });
        }
        let gaps = sector_samples
            .iter()
            .map(|sample| sample.gap_m)
            .collect::<Vec<_>>();
        if gaps
            .iter()
            .take(active)
            .any(|gap| !gap.is_finite() || *gap <= input.contact_exclusion.handoff_gap_m)
        {
            return Err(AirFilmError::Unavailable {
                reason: "nonpositive-or-contact-overlap-gap",
            });
        }
        let mean = |values: &[f64], field: &'static str| -> Result<f64, AirFilmError> {
            checked(
                values
                    .iter()
                    .take(active)
                    .try_fold(0.0, |sum, value| checked(sum + value, field))?
                    / active as f64,
                field,
            )
        };
        let rates = sector_samples
            .iter()
            .map(|sample| sample.gap_rate_m_per_s)
            .collect::<Vec<_>>();
        let speeds = sector_samples
            .iter()
            .map(|sample| sample.radial_relative_velocity_m_per_s)
            .collect::<Vec<_>>();
        let radial_angle = 2.0 * core::f64::consts::PI * (sector as f64 + 0.5)
            / input.discretization.azimuthal_sectors as f64;
        let (radial_direction, circumferential_direction) = sector_basis(normal, radial_angle)?;
        let gas_input = GasFilmInput {
            identity: GasFilmIdentity {
                case_id: format!("{}:sector:{sector}", input.identity.case_id),
                model_id: isothermal_compressible_reynolds_model_id().to_owned(),
                gas_species_id: input.identity.gas_species_id.clone(),
                eos_id: input.identity.eos_id.clone(),
                viscosity_source_id: input.identity.viscosity_source_id.clone(),
                thermal_model_id: input.identity.thermal_model_id.clone(),
                frame_id: input.identity.frame_id.clone(),
                deterministic_seed: input
                    .identity
                    .deterministic_seed
                    .wrapping_add(sector as u64),
                authority: input.identity.authority,
            },
            gas: input.gas,
            grid: GasFilmGrid1d {
                length_m: input.disc_radius_m,
                gap_m: gaps,
                contact_exclusion: ContactExclusionMask {
                    excluded: core::mem::take(&mut excluded),
                },
            },
            boundary: input.boundary,
            slip_policy: input.slip_policy.clone(),
            roughness_policy: input.roughness_policy.clone(),
            applicability: input.applicability,
            uncertainty: input.uncertainty,
            wall_motion: MovingWallInput {
                lower_tangential_velocity_m_per_s: 0.0,
                upper_tangential_velocity_m_per_s: mean(&speeds, "sector_mean_radial_speed")?,
                gap_rate_m_per_s: mean(&rates, "sector_mean_gap_rate")?,
            },
            initial_absolute_pressure_pa: input.initial_absolute_pressure_pa,
            timestep_s: input.timestep_s,
            budget: input.budget,
        };
        let step = solve_isothermal_gas_film_1d(
            &gas_input,
            checkpoint.map(|state| &state.sectors[sector]),
        )
        .map_err(|detail| AirFilmError::GasFilmRefusal { detail })?;
        accumulate_receipt(
            &mut gas_mass,
            &mut storage,
            &mut mass_flux,
            &mut enthalpy,
            &mut closure,
            &mut wall_power,
            &mut heat,
            &mut max_residual,
            &step.receipt,
            strip_width,
        )?;
        for (radial, pressure) in step.absolute_pressure_pa.iter().enumerate().take(active) {
            let pressure = pressure.ok_or(AirFilmError::NonFiniteDerived {
                field: "active_sector_pressure",
            })?;
            let area = strip_width * input.disc_radius_m / input.discretization.radial_cells as f64;
            let gauge_pressure = checked(
                pressure - input.gauge_reference_absolute_pressure_pa,
                "gauge_pressure_for_body_wrench",
            )?;
            let pressure_force = normal.scale(checked(gauge_pressure * area, "pressure_force")?);
            let shear =
                step.receipt.upper_wall_shear_pa[radial].ok_or(AirFilmError::NonFiniteDerived {
                    field: "active_sector_shear",
                })?;
            let shear_force =
                radial_direction.scale(checked(-shear * area, "shear_reaction_force")?);
            let circumferential_shear = checked(
                input.gas.dynamic_viscosity_pa_s
                    * sector_samples[radial].circumferential_relative_velocity_m_per_s
                    / sector_samples[radial].gap_m,
                "circumferential_wall_on_gas_shear",
            )?;
            let circumferential_force = circumferential_direction.scale(checked(
                -circumferential_shear * area,
                "circumferential_gas_on_body_force",
            )?);
            let local_power = checked(
                circumferential_shear
                    * sector_samples[radial].circumferential_relative_velocity_m_per_s
                    * area,
                "circumferential_couette_power",
            )?;
            circumferential_power = checked(
                circumferential_power + local_power,
                "circumferential_couette_power",
            )?;
            let local_force = pressure_force.add(shear_force).add(circumferential_force);
            force = force.add(local_force);
            moment = moment.add(sector_samples[radial].lever_arm_world_m.cross(local_force));
        }
        sector_checkpoints.push(step.checkpoint);
    }
    wall_power = checked(
        wall_power + circumferential_power,
        "total_wall_power_to_gas",
    )?;
    heat = checked(heat + circumferential_power, "total_viscous_heat")?;
    force = checked_vec(force, "net_air_force_world_n")?;
    moment = checked_vec(moment, "net_air_moment_about_com_world_n_m")?;
    Ok(AirFilmStep {
        samples,
        receipt: AirFilmReceipt {
            wrench: AirFilmWrench {
                force_world_n: force,
                moment_about_com_world_n_m: moment,
            },
            gas_mass_kg: gas_mass,
            storage_rate_kg_per_s: storage,
            outward_mass_flux_kg_per_s: mass_flux,
            outward_enthalpy_flux_w: enthalpy,
            mass_closure_residual_kg_per_s: closure,
            wall_power_to_gas_w: wall_power,
            viscous_heat_w: heat,
            circumferential_couette_power_w: circumferential_power,
            circumferential_couette_heat_w: circumferential_power,
            equal_area_strip_width_m: strip_width,
            signed_first_radial_moment_discrepancy_m3: first_moment_discrepancy,
            max_sector_mass_residual_kg_m2_s: max_residual,
            gauge_reference_absolute_pressure_pa: input.gauge_reference_absolute_pressure_pa,
            input_authority: input.identity.authority,
        },
        checkpoint: AirFilmCheckpoint {
            adapter_model_id: input.identity.adapter_model_id.clone(),
            case_id: input.identity.case_id.clone(),
            configuration_id: input.identity.configuration_id.clone(),
            azimuthal_sectors: input.discretization.azimuthal_sectors,
            sectors: sector_checkpoints,
        },
    })
}

#[allow(clippy::too_many_arguments)]
fn accumulate_receipt(
    gas_mass: &mut f64,
    storage: &mut f64,
    mass_flux: &mut f64,
    enthalpy: &mut f64,
    closure: &mut f64,
    wall_power: &mut f64,
    heat: &mut f64,
    max_residual: &mut f64,
    receipt: &GasFilmReceipt,
    width: f64,
) -> Result<(), AirFilmError> {
    *gas_mass = checked(*gas_mass + receipt.gas_mass_kg_per_m * width, "gas_mass_kg")?;
    *storage = checked(
        *storage + receipt.storage_rate_kg_per_m_s * width,
        "storage_rate",
    )?;
    *mass_flux = checked(
        *mass_flux
            + width
                * (receipt.left_boundary_outward_mass_flux_kg_per_m_s
                    + receipt.right_boundary_outward_mass_flux_kg_per_m_s
                    + receipt.vent_outward_mass_flux_kg_per_m_s),
        "outward_mass_flux",
    )?;
    *enthalpy = checked(
        *enthalpy
            + width
                * (receipt.left_boundary_outward_enthalpy_flux_w_per_m
                    + receipt.right_boundary_outward_enthalpy_flux_w_per_m
                    + receipt.vent_outward_enthalpy_flux_w_per_m),
        "outward_enthalpy_flux",
    )?;
    *closure = checked(
        *closure + receipt.mass_closure_residual_kg_per_m_s * width,
        "mass_closure",
    )?;
    *wall_power = checked(
        *wall_power + receipt.wall_power_to_gas_w_per_m * width,
        "wall_power",
    )?;
    *heat = checked(*heat + receipt.viscous_heat_w_per_m * width, "viscous_heat")?;
    *max_residual = max_residual.max(receipt.max_mass_residual_kg_m2_s);
    Ok(())
}
