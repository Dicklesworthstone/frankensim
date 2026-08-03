//! Bounded one-dimensional isothermal compressible gas-film Reynolds solver.
//!
//! This is deliberately a small, generic continuum foundation: a uniform
//! Cartesian finite-volume line of unit out-of-plane width.  It solves the
//! isothermal ideal-gas Reynolds mass balance, not a liquid cavitation or EHL
//! model.  Every pressure is absolute [Pa], every gap is explicit [m], and an
//! excluded contact cell is absent from the gas state rather than receiving an
//! anonymous lower-gap floor.

use core::fmt;
use core::fmt::Write;

const MODEL_ID: &str = "fs-flux/isothermal-compressible-reynolds-1d-v1";
const MAX_GEOMETRIC_CELLS: usize = 4_096;
const RESTART_GAP_EVOLUTION_ABSOLUTE_TOLERANCE_M: f64 = 1.0e-15;
const RESTART_GAP_EVOLUTION_RELATIVE_TOLERANCE: f64 = 1.0e-12;

/// Typed refusal or bounded-solve failure for the gas-film foundation.
#[derive(Debug, Clone, PartialEq)]
pub enum GasFilmError {
    /// A named scalar, vector entry, or count is not admitted.
    InvalidInput { field: &'static str },
    /// A semantic identity is empty or contains non-canonical bytes.
    InvalidIdentity { field: &'static str },
    /// A requested model lies outside this deliberately narrow continuum scope.
    Unavailable { reason: &'static str },
    /// The static mask would split or change the retained one-dimensional topology.
    TopologyChangeUnavailable,
    /// A checkpoint is not compatible with the next requested step.
    CheckpointMismatch { field: &'static str },
    /// A derived quantity overflowed or became non-finite.
    NonFiniteDerived { field: &'static str },
    /// The deterministic Picard/Jacobi budget ended before the residual gate.
    IterationBudgetExceeded {
        /// Number of completed iterations.
        iterations: u32,
        /// Largest retained cell mass residual [kg m^-2 s^-1].
        max_mass_residual_kg_m2_s: f64,
    },
}

impl fmt::Display for GasFilmError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for GasFilmError {}

/// Caller-declared authority retained in results; this crate does not mint it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GasFilmInputAuthority {
    /// Inputs are caller declarations, not a materials/EOS admission.
    CallerDeclared,
    /// Inputs are a declared synthetic fixture.
    SyntheticFixture,
}

/// Canonical identities and five-explicit metadata for one film request.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GasFilmIdentity {
    /// Stable caller case identifier.
    pub case_id: String,
    /// Model identity, normally [`MODEL_ID`].
    pub model_id: String,
    /// Gas species/card identity.
    pub gas_species_id: String,
    /// Isothermal ideal-gas EOS identity.
    pub eos_id: String,
    /// Dynamic-viscosity source identity.
    pub viscosity_source_id: String,
    /// Thermal-model identity; this implementation accepts isothermal only.
    pub thermal_model_id: String,
    /// Declared Cartesian frame for the positive line coordinate and wall speeds.
    pub frame_id: String,
    /// Deterministic caller seed retained for replay even though this solver samples no RNG.
    pub deterministic_seed: u64,
    /// Authority of the supplied gas inputs.
    pub authority: GasFilmInputAuthority,
}

/// Isothermal ideal-gas properties in coherent SI units.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IsothermalIdealGas {
    /// Specific gas constant [J kg^-1 K^-1].
    pub specific_gas_constant_j_kg_k: f64,
    /// Uniform absolute temperature [K].
    pub temperature_k: f64,
    /// Dynamic viscosity [Pa s].
    pub dynamic_viscosity_pa_s: f64,
    /// Caller-declared density at the supplied initial absolute pressure [kg m^-3].
    /// It is checked against `p / (R T)`; it is not an independent EOS.
    pub declared_density_kg_m3: f64,
    /// Caller-declared specific enthalpy [J kg^-1] used only to transport a
    /// signed boundary enthalpy receipt. It is not a thermochemical model.
    pub declared_specific_enthalpy_j_kg: f64,
}

/// Boundary topology admitted by the 1-D line model.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GasFilmBoundaryTopology {
    /// Both ends have zero gas mass flux.
    Sealed,
    /// Left and right exterior pressure reservoirs [Pa absolute].
    Open {
        /// Pressure at the left exterior face [Pa absolute].
        left_absolute_pressure_pa: f64,
        /// Pressure at the right exterior face [Pa absolute].
        right_absolute_pressure_pa: f64,
    },
    /// Sealed outer ends with one explicit internal pressure vent.
    Vented {
        /// Active-cell index of the vent.
        cell_index: usize,
        /// Vent reservoir pressure [Pa absolute].
        absolute_pressure_pa: f64,
    },
}

/// Explicit static contact handoff mask. `true` means no gas-film cell exists.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ContactExclusionMask {
    /// One flag per geometric cell.
    pub excluded: Vec<bool>,
}

/// Uniform structured line geometry with an exact per-cell gap [m].
#[derive(Debug, Clone, PartialEq)]
pub struct GasFilmGrid1d {
    /// Total line length [m].
    pub length_m: f64,
    /// Current positive gap at each geometric cell [m].
    pub gap_m: Vec<f64>,
    /// Static gas/contact partition.
    pub contact_exclusion: ContactExclusionMask,
}

/// Relative wall motion used by the Reynolds mass flux and Couette receipt.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MovingWallInput {
    /// Lower-wall line velocity [m s^-1].
    pub lower_tangential_velocity_m_per_s: f64,
    /// Upper-wall line velocity [m s^-1].
    pub upper_tangential_velocity_m_per_s: f64,
    /// Uniform `dh/dt` [m s^-1]; negative values squeeze the gas.
    pub gap_rate_m_per_s: f64,
}

/// Explicit continuum-domain gates.  There is no slip correction in this slice.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GasFilmApplicability {
    /// Caller-supplied molecular mean free path [m] for the declared state.
    pub mean_free_path_m: f64,
    /// Maximum admitted cell Knudsen number `lambda / h`.
    pub maximum_knudsen_number: f64,
    /// Maximum admitted adjacent-cell slope `abs(dh/dx)`.
    pub maximum_gap_slope: f64,
    /// Maximum admitted line Mach number based on caller-declared sound speed.
    pub speed_of_sound_m_per_s: f64,
    /// Maximum admitted wall Mach number.
    pub maximum_mach_number: f64,
}

/// Caller-declared input uncertainty bounds retained without propagation or certification.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GasFilmUncertainty {
    /// Relative viscosity uncertainty bound [1].
    pub viscosity_relative_bound: f64,
    /// Relative gap uncertainty bound [1].
    pub gap_relative_bound: f64,
    /// Relative boundary/initial pressure uncertainty bound [1].
    pub pressure_relative_bound: f64,
}

/// This foundation only admits no-slip continuum flow.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SlipPolicy {
    /// Named no-slip continuum policy.
    NoSlipContinuum { source_id: String },
    /// Any rarefied/slip correction requires a separately admitted model.
    RarefiedSlipRequested { source_id: String },
}

/// Explicit wall-roughness closure for the smooth-gap Reynolds model.
#[derive(Debug, Clone, PartialEq)]
pub enum RoughnessPolicy {
    /// A caller-declared upper bound lies below every active gap.
    ResolvedSmooth {
        /// Identity of the supplied roughness source.
        source_id: String,
        /// Maximum asperity height [m].
        maximum_roughness_m: f64,
    },
    /// Roughness requires a model not provided by this foundation.
    Unresolved { source_id: String },
}

/// Bounded nonlinear iteration and numerical acceptance controls.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GasFilmBudget {
    /// Maximum deterministic sweeps; zero is refused.
    pub maximum_iterations: u32,
    /// Maximum cell mass-balance residual [kg m^-2 s^-1].
    pub mass_residual_tolerance_kg_m2_s: f64,
    /// Damped Jacobi factor in `(0, 1]`.
    pub relaxation: f64,
}

/// One complete one-step film request.
#[derive(Debug, Clone, PartialEq)]
pub struct GasFilmInput {
    /// Semantic identities and replay metadata.
    pub identity: GasFilmIdentity,
    /// Gas EOS and viscosity values.
    pub gas: IsothermalIdealGas,
    /// Grid and explicit contact exclusion.
    pub grid: GasFilmGrid1d,
    /// Open, sealed, or explicit-vent boundary topology.
    pub boundary: GasFilmBoundaryTopology,
    /// No-slip/rarefaction policy.
    pub slip_policy: SlipPolicy,
    /// Smooth-wall applicability or an explicit refusal boundary.
    pub roughness_policy: RoughnessPolicy,
    /// Continuum/thin-film admission envelope.
    pub applicability: GasFilmApplicability,
    /// Caller-declared uncertainties; no uncertainty propagation is claimed.
    pub uncertainty: GasFilmUncertainty,
    /// Moving-wall Couette and squeeze data.
    pub wall_motion: MovingWallInput,
    /// Initial absolute pressure for a fresh state [Pa].
    pub initial_absolute_pressure_pa: f64,
    /// Absolute pressure used only as the named gauge-load reference [Pa].
    ///
    /// This is intentionally distinct from the initial condition, so a restart
    /// or non-equilibrium initial field cannot silently redefine a reported
    /// gauge load.
    pub gauge_reference_absolute_pressure_pa: f64,
    /// Fixed physical step size [s].
    pub timestep_s: f64,
    /// Bounded numerical controls.
    pub budget: GasFilmBudget,
}

/// Restartable state after a fully converged accepted step.
#[derive(Debug, Clone, PartialEq)]
pub struct GasFilmCheckpoint {
    /// Step count beginning at one.
    pub step_index: u64,
    /// Model identity bound into the accepted state.
    pub model_id: String,
    /// Caller case identity bound into the accepted state.
    pub case_id: String,
    /// Input authority retained without promotion across restarts.
    pub input_authority: GasFilmInputAuthority,
    /// Caller-declared input uncertainty retained without propagation.
    pub input_uncertainty: GasFilmUncertainty,
    /// Deterministic invariant-configuration fingerprint required for restart.
    pub configuration_fingerprint: String,
    /// Full accepted request fingerprint retained for provenance only.
    ///
    /// It includes dynamic gap and wall-motion values, so it is deliberately
    /// not compared to a subsequent physically evolving request.
    pub accepted_request_fingerprint: String,
    /// Pressure at active cells only [Pa absolute].
    pub active_absolute_pressure_pa: Vec<f64>,
    /// Exact active-cell gaps at the accepted step [m].
    pub active_gap_m: Vec<f64>,
    /// Number of active prefix cells bound into this checkpoint.
    pub active_cells: usize,
}

/// Pressure/shear/load/work and conservative accounting receipt.
#[derive(Debug, Clone, PartialEq)]
pub struct GasFilmReceipt {
    /// Model identity retained from the request.
    pub model_id: String,
    /// Input authority retained without promotion.
    pub input_authority: GasFilmInputAuthority,
    /// Number of deterministic sweeps used.
    pub iterations: u32,
    /// Maximum accepted cell residual [kg m^-2 s^-1].
    pub max_mass_residual_kg_m2_s: f64,
    /// Gas mass per unit out-of-plane width at the end [kg m^-1].
    pub gas_mass_kg_per_m: f64,
    /// Storage rate `(m_new-m_old)/dt` [kg m^-1 s^-1].
    pub storage_rate_kg_per_m_s: f64,
    /// Net outward left-boundary mass flux [kg m^-1 s^-1].
    pub left_boundary_outward_mass_flux_kg_per_m_s: f64,
    /// Net outward right-boundary mass flux [kg m^-1 s^-1].
    pub right_boundary_outward_mass_flux_kg_per_m_s: f64,
    /// Net outward vent mass flux [kg m^-1 s^-1].
    pub vent_outward_mass_flux_kg_per_m_s: f64,
    /// Signed left-boundary enthalpy transport [W m^-1], from caller-declared h.
    pub left_boundary_outward_enthalpy_flux_w_per_m: f64,
    /// Signed right-boundary enthalpy transport [W m^-1], from caller-declared h.
    pub right_boundary_outward_enthalpy_flux_w_per_m: f64,
    /// Signed vent enthalpy transport [W m^-1], from caller-declared h.
    pub vent_outward_enthalpy_flux_w_per_m: f64,
    /// Conservative global storage-plus-boundary residual [kg m^-1 s^-1].
    pub mass_closure_residual_kg_per_m_s: f64,
    /// Integral of absolute pressure over active gas area [N m^-1].
    pub absolute_pressure_load_n_per_m: f64,
    /// Integral of gauge pressure relative to the declared gauge reference [N m^-1].
    pub gauge_pressure_load_n_per_m: f64,
    /// Per-cell normal traction of gas on the upper wall [Pa]; positive is compressive.
    pub gas_on_upper_wall_normal_traction_pa: Vec<Option<f64>>,
    /// Per-cell normal traction of upper wall on gas [Pa]; exactly the negative of gas-on-wall.
    pub upper_wall_on_gas_normal_traction_pa: Vec<Option<f64>>,
    /// Per-cell tangential traction of gas on the upper wall [Pa].
    pub gas_on_upper_wall_tangential_traction_pa: Vec<Option<f64>>,
    /// Per-cell tangential traction of upper wall on gas [Pa]; exactly the negative of gas-on-wall.
    pub upper_wall_on_gas_tangential_traction_pa: Vec<Option<f64>>,
    /// Legacy-named alias for upper-wall-on-gas tangential traction [Pa].
    /// Excluded cells are `None`; new callers should use the explicit field above.
    pub upper_wall_shear_pa: Vec<Option<f64>>,
    /// Normal moving-gap mechanical power into gas, `-integral(p dh/dt dx)` [W m^-1].
    pub normal_gap_power_to_gas_w_per_m: f64,
    /// Relative-wall mechanical power into gas per unit width [W m^-1].
    pub wall_power_to_gas_w_per_m: f64,
    /// Isothermal viscous heat receipt, equal to non-negative relative-wall power [W m^-1].
    pub viscous_heat_w_per_m: f64,
    /// Maximum admitted Knudsen number.
    pub maximum_knudsen_number: f64,
    /// Maximum adjacent gap slope.
    pub maximum_gap_slope: f64,
}

/// Fully accepted pressure field, receipt, and restart checkpoint.
#[derive(Debug, Clone, PartialEq)]
pub struct GasFilmStep {
    /// One value for every geometric cell; excluded contact cells carry no gas pressure.
    pub absolute_pressure_pa: Vec<Option<f64>>,
    /// Conservative accounting and mechanical receipt.
    pub receipt: GasFilmReceipt,
    /// State that may seed an identical-topology next step.
    pub checkpoint: GasFilmCheckpoint,
}

fn canonical_identity(value: &str, field: &'static str) -> Result<(), GasFilmError> {
    if value.is_empty()
        || !value.bytes().all(|byte| {
            byte.is_ascii_alphanumeric() || matches!(byte, b'-' | b'_' | b'.' | b':' | b'/')
        })
    {
        return Err(GasFilmError::InvalidIdentity { field });
    }
    Ok(())
}

fn positive(value: f64, field: &'static str) -> Result<(), GasFilmError> {
    if !(value.is_finite() && value > 0.0) {
        return Err(GasFilmError::InvalidInput { field });
    }
    Ok(())
}

fn nonnegative(value: f64, field: &'static str) -> Result<(), GasFilmError> {
    if !(value.is_finite() && value >= 0.0) {
        return Err(GasFilmError::InvalidInput { field });
    }
    Ok(())
}

fn finite(value: f64, field: &'static str) -> Result<(), GasFilmError> {
    if !value.is_finite() {
        return Err(GasFilmError::InvalidInput { field });
    }
    Ok(())
}

fn derived(value: f64, field: &'static str) -> Result<f64, GasFilmError> {
    if !value.is_finite() {
        return Err(GasFilmError::NonFiniteDerived { field });
    }
    Ok(value)
}

fn active_prefix(mask: &[bool]) -> Result<usize, GasFilmError> {
    let active = mask.iter().take_while(|excluded| !**excluded).count();
    if active < 2 {
        return Err(GasFilmError::Unavailable {
            reason: "fewer-than-two-active-gas-cells",
        });
    }
    if mask[active..].iter().any(|excluded| !*excluded) {
        return Err(GasFilmError::TopologyChangeUnavailable);
    }
    Ok(active)
}

fn face_flux(
    left_pressure_pa: f64,
    right_pressure_pa: f64,
    left_gap_m: f64,
    right_gap_m: f64,
    distance_m: f64,
    gas: IsothermalIdealGas,
    mean_wall_velocity_m_per_s: f64,
) -> Result<f64, GasFilmError> {
    let pressure = derived(
        0.5 * (left_pressure_pa + right_pressure_pa),
        "face_pressure_pa",
    )?;
    let gap = derived(0.5 * (left_gap_m + right_gap_m), "face_gap_m")?;
    positive(pressure, "face_pressure_pa")?;
    positive(gap, "face_gap_m")?;
    let rt = derived(
        gas.specific_gas_constant_j_kg_k * gas.temperature_k,
        "specific_gas_constant_times_temperature",
    )?;
    let conductance = derived(
        gap.powi(3) * pressure / (12.0 * gas.dynamic_viscosity_pa_s * rt),
        "reynolds_mass_conductance",
    )?;
    let pressure_flux = derived(
        -conductance * (right_pressure_pa - left_pressure_pa) / distance_m,
        "pressure_mass_flux",
    )?;
    let couette_flux = derived(
        pressure * gap * mean_wall_velocity_m_per_s / rt,
        "couette_mass_flux",
    )?;
    let flux = derived(pressure_flux + couette_flux, "face_mass_flux")?;
    Ok(flux)
}

impl GasFilmInput {
    fn validate(&self) -> Result<(usize, f64, f64, f64, f64), GasFilmError> {
        for (field, value) in [
            ("case_id", self.identity.case_id.as_str()),
            ("model_id", self.identity.model_id.as_str()),
            ("gas_species_id", self.identity.gas_species_id.as_str()),
            ("eos_id", self.identity.eos_id.as_str()),
            (
                "viscosity_source_id",
                self.identity.viscosity_source_id.as_str(),
            ),
            ("thermal_model_id", self.identity.thermal_model_id.as_str()),
            ("frame_id", self.identity.frame_id.as_str()),
        ] {
            canonical_identity(value, field)?;
        }
        if self.identity.model_id != MODEL_ID {
            return Err(GasFilmError::Unavailable {
                reason: "unrecognized-gas-film-model-identity",
            });
        }
        match &self.slip_policy {
            SlipPolicy::NoSlipContinuum { source_id } => {
                canonical_identity(source_id, "slip_source_id")?
            }
            SlipPolicy::RarefiedSlipRequested { .. } => {
                return Err(GasFilmError::Unavailable {
                    reason: "slip-or-rarefied-model-not-implemented",
                });
            }
        }
        match &self.roughness_policy {
            RoughnessPolicy::ResolvedSmooth {
                source_id,
                maximum_roughness_m,
            } => {
                canonical_identity(source_id, "roughness_source_id")?;
                nonnegative(*maximum_roughness_m, "roughness.maximum_roughness_m")?;
            }
            RoughnessPolicy::Unresolved { source_id } => {
                canonical_identity(source_id, "roughness_source_id")?;
                return Err(GasFilmError::Unavailable {
                    reason: "roughness-model-not-admitted",
                });
            }
        }
        positive(self.grid.length_m, "grid.length_m")?;
        if self.grid.gap_m.len() != self.grid.contact_exclusion.excluded.len() {
            return Err(GasFilmError::InvalidInput {
                field: "contact_exclusion.length",
            });
        }
        let active = active_prefix(&self.grid.contact_exclusion.excluded)?;
        let cells = self.grid.gap_m.len();
        if cells < 2 {
            return Err(GasFilmError::InvalidInput {
                field: "grid.cells",
            });
        }
        if cells > MAX_GEOMETRIC_CELLS {
            return Err(GasFilmError::Unavailable {
                reason: "gas-film-grid-exceeds-bounded-cell-cap",
            });
        }
        let spacing = derived(self.grid.length_m / cells as f64, "grid.spacing_m")?;
        positive(spacing, "grid.spacing_m")?;
        positive(
            self.gas.specific_gas_constant_j_kg_k,
            "gas.specific_gas_constant_j_kg_k",
        )?;
        positive(self.gas.temperature_k, "gas.temperature_k")?;
        positive(
            self.gas.dynamic_viscosity_pa_s,
            "gas.dynamic_viscosity_pa_s",
        )?;
        positive(
            self.gas.declared_density_kg_m3,
            "gas.declared_density_kg_m3",
        )?;
        positive(
            self.gas.declared_specific_enthalpy_j_kg,
            "gas.declared_specific_enthalpy_j_kg",
        )?;
        positive(
            self.initial_absolute_pressure_pa,
            "initial_absolute_pressure_pa",
        )?;
        positive(
            self.gauge_reference_absolute_pressure_pa,
            "gauge_reference_absolute_pressure_pa",
        )?;
        positive(self.timestep_s, "timestep_s")?;
        positive(
            self.applicability.mean_free_path_m,
            "applicability.mean_free_path_m",
        )?;
        positive(
            self.applicability.maximum_knudsen_number,
            "applicability.maximum_knudsen_number",
        )?;
        positive(
            self.applicability.maximum_gap_slope,
            "applicability.maximum_gap_slope",
        )?;
        positive(
            self.applicability.speed_of_sound_m_per_s,
            "applicability.speed_of_sound_m_per_s",
        )?;
        positive(
            self.applicability.maximum_mach_number,
            "applicability.maximum_mach_number",
        )?;
        nonnegative(
            self.uncertainty.viscosity_relative_bound,
            "uncertainty.viscosity_relative_bound",
        )?;
        nonnegative(
            self.uncertainty.gap_relative_bound,
            "uncertainty.gap_relative_bound",
        )?;
        nonnegative(
            self.uncertainty.pressure_relative_bound,
            "uncertainty.pressure_relative_bound",
        )?;
        if self.budget.maximum_iterations == 0 {
            return Err(GasFilmError::InvalidInput {
                field: "budget.maximum_iterations",
            });
        }
        positive(
            self.budget.mass_residual_tolerance_kg_m2_s,
            "budget.mass_residual_tolerance_kg_m2_s",
        )?;
        if !(self.budget.relaxation.is_finite()
            && self.budget.relaxation > 0.0
            && self.budget.relaxation <= 1.0)
        {
            return Err(GasFilmError::InvalidInput {
                field: "budget.relaxation",
            });
        }
        finite(
            self.wall_motion.lower_tangential_velocity_m_per_s,
            "wall_motion.lower_tangential_velocity_m_per_s",
        )?;
        finite(
            self.wall_motion.upper_tangential_velocity_m_per_s,
            "wall_motion.upper_tangential_velocity_m_per_s",
        )?;
        finite(
            self.wall_motion.gap_rate_m_per_s,
            "wall_motion.gap_rate_m_per_s",
        )?;
        let wall_speed = self
            .wall_motion
            .lower_tangential_velocity_m_per_s
            .abs()
            .max(self.wall_motion.upper_tangential_velocity_m_per_s.abs());
        let mach = derived(
            wall_speed / self.applicability.speed_of_sound_m_per_s,
            "wall_mach_number",
        )?;
        if mach > self.applicability.maximum_mach_number {
            return Err(GasFilmError::Unavailable {
                reason: "wall-mach-outside-declared-envelope",
            });
        }
        let rt = derived(
            self.gas.specific_gas_constant_j_kg_k * self.gas.temperature_k,
            "specific_gas_constant_times_temperature",
        )?;
        let eos_density = derived(
            self.initial_absolute_pressure_pa / rt,
            "initial_eos_density_kg_m3",
        )?;
        let density_scale = eos_density.max(self.gas.declared_density_kg_m3).max(1.0);
        if (eos_density - self.gas.declared_density_kg_m3).abs() > 1.0e-10 * density_scale {
            return Err(GasFilmError::InvalidInput {
                field: "gas.declared_density_kg_m3",
            });
        }
        let mut max_knudsen = 0.0_f64;
        let mut max_slope = 0.0_f64;
        for (index, gap) in self.grid.gap_m.iter().take(active).copied().enumerate() {
            positive(gap, "grid.gap_m")?;
            if let RoughnessPolicy::ResolvedSmooth {
                maximum_roughness_m,
                ..
            } = &self.roughness_policy
            {
                if *maximum_roughness_m >= gap {
                    return Err(GasFilmError::Unavailable {
                        reason: "roughness-not-small-relative-to-gap",
                    });
                }
            }
            let knudsen = derived(self.applicability.mean_free_path_m / gap, "knudsen_number")?;
            max_knudsen = max_knudsen.max(knudsen);
            if index > 0 {
                let slope = derived(
                    (gap - self.grid.gap_m[index - 1]).abs() / spacing,
                    "gap_slope",
                )?;
                max_slope = max_slope.max(slope);
            }
        }
        if max_knudsen > self.applicability.maximum_knudsen_number {
            return Err(GasFilmError::Unavailable {
                reason: "rarefied-knudsen-outside-continuum-envelope",
            });
        }
        if max_slope > self.applicability.maximum_gap_slope {
            return Err(GasFilmError::Unavailable {
                reason: "large-gap-slope-outside-reynolds-envelope",
            });
        }
        match self.boundary {
            GasFilmBoundaryTopology::Sealed => {}
            GasFilmBoundaryTopology::Open {
                left_absolute_pressure_pa,
                right_absolute_pressure_pa,
            } => {
                positive(
                    left_absolute_pressure_pa,
                    "boundary.left_absolute_pressure_pa",
                )?;
                positive(
                    right_absolute_pressure_pa,
                    "boundary.right_absolute_pressure_pa",
                )?;
            }
            GasFilmBoundaryTopology::Vented {
                cell_index,
                absolute_pressure_pa,
            } => {
                if cell_index >= active {
                    return Err(GasFilmError::InvalidInput {
                        field: "boundary.vent_cell_index",
                    });
                }
                positive(absolute_pressure_pa, "boundary.vent_absolute_pressure_pa")?;
            }
        }
        Ok((active, spacing, rt, max_knudsen, max_slope))
    }
}

fn initial_state(
    input: &GasFilmInput,
    active: usize,
    checkpoint: Option<&GasFilmCheckpoint>,
) -> Result<(Vec<f64>, Vec<f64>, u64), GasFilmError> {
    match checkpoint {
        None => {
            let old_gap = input
                .grid
                .gap_m
                .iter()
                .take(active)
                .map(|gap| gap - input.wall_motion.gap_rate_m_per_s * input.timestep_s)
                .collect::<Vec<_>>();
            for gap in &old_gap {
                positive(*gap, "derived.previous_gap_m")?;
            }
            Ok((vec![input.initial_absolute_pressure_pa; active], old_gap, 0))
        }
        Some(state) => {
            if state.model_id != input.identity.model_id {
                return Err(GasFilmError::CheckpointMismatch { field: "model_id" });
            }
            if state.case_id != input.identity.case_id {
                return Err(GasFilmError::CheckpointMismatch { field: "case_id" });
            }
            if state.input_authority != input.identity.authority {
                return Err(GasFilmError::CheckpointMismatch {
                    field: "input_authority",
                });
            }
            if state.configuration_fingerprint != invariant_configuration_fingerprint(input) {
                return Err(GasFilmError::CheckpointMismatch {
                    field: "configuration_fingerprint",
                });
            }
            if state.active_cells != active
                || state.active_absolute_pressure_pa.len() != active
                || state.active_gap_m.len() != active
            {
                return Err(GasFilmError::CheckpointMismatch {
                    field: "active_cells_or_state_length",
                });
            }
            for old in &state.active_gap_m {
                positive(*old, "checkpoint.active_gap_m")?;
            }
            for pressure in &state.active_absolute_pressure_pa {
                positive(*pressure, "checkpoint.active_absolute_pressure_pa")?;
            }
            for (new_gap, old_gap) in input
                .grid
                .gap_m
                .iter()
                .take(active)
                .zip(&state.active_gap_m)
            {
                let expected = derived(
                    *old_gap + input.wall_motion.gap_rate_m_per_s * input.timestep_s,
                    "checkpoint_expected_active_gap_m",
                )?;
                let tolerance = RESTART_GAP_EVOLUTION_ABSOLUTE_TOLERANCE_M
                    + RESTART_GAP_EVOLUTION_RELATIVE_TOLERANCE * expected.abs().max(new_gap.abs());
                if (*new_gap - expected).abs() > tolerance {
                    return Err(GasFilmError::CheckpointMismatch {
                        field: "active_gap_evolution",
                    });
                }
            }
            Ok((
                state.active_absolute_pressure_pa.clone(),
                state.active_gap_m.clone(),
                state.step_index,
            ))
        }
    }
}

/// Fingerprint of fields that cannot change across a restart sequence.
///
/// Current gaps and wall velocities/rate are deliberately excluded.  The
/// latter are per-step data, with gap evolution checked geometrically against
/// the accepted state before iteration.
fn invariant_configuration_fingerprint(input: &GasFilmInput) -> String {
    fn text(output: &mut String, label: &str, value: &str) {
        let _ = write!(output, "{label}={value};");
    }
    fn number(output: &mut String, label: &str, value: f64) {
        let _ = write!(output, "{label}={:016x};", value.to_bits());
    }
    let mut output = String::new();
    text(&mut output, "case", &input.identity.case_id);
    text(&mut output, "model", &input.identity.model_id);
    text(&mut output, "species", &input.identity.gas_species_id);
    text(&mut output, "eos", &input.identity.eos_id);
    text(
        &mut output,
        "viscosity_source",
        &input.identity.viscosity_source_id,
    );
    text(&mut output, "thermal", &input.identity.thermal_model_id);
    text(&mut output, "frame", &input.identity.frame_id);
    let _ = write!(
        output,
        "seed={};authority={:?};",
        input.identity.deterministic_seed, input.identity.authority
    );
    number(&mut output, "R", input.gas.specific_gas_constant_j_kg_k);
    number(&mut output, "T", input.gas.temperature_k);
    number(&mut output, "mu", input.gas.dynamic_viscosity_pa_s);
    number(
        &mut output,
        "rho_declared",
        input.gas.declared_density_kg_m3,
    );
    number(
        &mut output,
        "h_declared",
        input.gas.declared_specific_enthalpy_j_kg,
    );
    number(&mut output, "length", input.grid.length_m);
    let _ = write!(output, "cells={};", input.grid.gap_m.len());
    for excluded in &input.grid.contact_exclusion.excluded {
        let _ = write!(output, "mask={excluded};");
    }
    match input.boundary {
        GasFilmBoundaryTopology::Sealed => text(&mut output, "boundary_kind", "sealed"),
        GasFilmBoundaryTopology::Open { .. } => text(&mut output, "boundary_kind", "open"),
        GasFilmBoundaryTopology::Vented { .. } => text(&mut output, "boundary_kind", "vented"),
    }
    match &input.slip_policy {
        SlipPolicy::NoSlipContinuum { source_id } => text(&mut output, "slip", source_id),
        SlipPolicy::RarefiedSlipRequested { source_id } => {
            text(&mut output, "rarefied_slip", source_id)
        }
    }
    match &input.roughness_policy {
        RoughnessPolicy::ResolvedSmooth {
            source_id,
            maximum_roughness_m,
        } => {
            text(&mut output, "roughness", source_id);
            number(&mut output, "roughness_max", *maximum_roughness_m);
        }
        RoughnessPolicy::Unresolved { source_id } => {
            text(&mut output, "roughness_unresolved", source_id)
        }
    }
    number(
        &mut output,
        "mean_free_path",
        input.applicability.mean_free_path_m,
    );
    number(
        &mut output,
        "knudsen_max",
        input.applicability.maximum_knudsen_number,
    );
    number(
        &mut output,
        "slope_max",
        input.applicability.maximum_gap_slope,
    );
    number(
        &mut output,
        "sound_speed",
        input.applicability.speed_of_sound_m_per_s,
    );
    number(
        &mut output,
        "mach_max",
        input.applicability.maximum_mach_number,
    );
    number(
        &mut output,
        "uncertainty_mu",
        input.uncertainty.viscosity_relative_bound,
    );
    number(
        &mut output,
        "uncertainty_gap",
        input.uncertainty.gap_relative_bound,
    );
    number(
        &mut output,
        "uncertainty_pressure",
        input.uncertainty.pressure_relative_bound,
    );
    number(
        &mut output,
        "gauge_reference",
        input.gauge_reference_absolute_pressure_pa,
    );
    number(&mut output, "timestep", input.timestep_s);
    let _ = write!(output, "iterations={};", input.budget.maximum_iterations);
    number(
        &mut output,
        "residual_tolerance",
        input.budget.mass_residual_tolerance_kg_m2_s,
    );
    number(&mut output, "relaxation", input.budget.relaxation);
    output
}

/// Deterministic, complete restart identity for every semantic input.
///
/// Floats are represented by IEEE-754 bits after validation, preserving exact
/// caller values without locale/formatting dependence.  This is an identity
/// guard, not a cryptographic digest or physical-validation certificate.
fn accepted_request_fingerprint(input: &GasFilmInput) -> String {
    fn text(output: &mut String, label: &str, value: &str) {
        let _ = write!(output, "{label}={value};");
    }
    fn number(output: &mut String, label: &str, value: f64) {
        let _ = write!(output, "{label}={:016x};", value.to_bits());
    }
    fn count(output: &mut String, label: &str, value: usize) {
        let _ = write!(output, "{label}={value};");
    }
    let mut output = String::new();
    text(&mut output, "case", &input.identity.case_id);
    text(&mut output, "model", &input.identity.model_id);
    text(&mut output, "species", &input.identity.gas_species_id);
    text(&mut output, "eos", &input.identity.eos_id);
    text(
        &mut output,
        "viscosity_source",
        &input.identity.viscosity_source_id,
    );
    text(&mut output, "thermal", &input.identity.thermal_model_id);
    text(&mut output, "frame", &input.identity.frame_id);
    let _ = write!(
        output,
        "seed={};authority={:?};",
        input.identity.deterministic_seed, input.identity.authority
    );
    number(&mut output, "R", input.gas.specific_gas_constant_j_kg_k);
    number(&mut output, "T", input.gas.temperature_k);
    number(&mut output, "mu", input.gas.dynamic_viscosity_pa_s);
    number(
        &mut output,
        "rho_declared",
        input.gas.declared_density_kg_m3,
    );
    number(
        &mut output,
        "h_declared",
        input.gas.declared_specific_enthalpy_j_kg,
    );
    number(&mut output, "length", input.grid.length_m);
    count(&mut output, "gap_count", input.grid.gap_m.len());
    for (index, gap) in input.grid.gap_m.iter().copied().enumerate() {
        number(&mut output, &format!("gap_{index}"), gap);
    }
    count(
        &mut output,
        "mask_count",
        input.grid.contact_exclusion.excluded.len(),
    );
    for (index, excluded) in input
        .grid
        .contact_exclusion
        .excluded
        .iter()
        .copied()
        .enumerate()
    {
        let _ = write!(output, "mask_{index}={excluded};");
    }
    match input.boundary {
        GasFilmBoundaryTopology::Sealed => text(&mut output, "boundary", "sealed"),
        GasFilmBoundaryTopology::Open {
            left_absolute_pressure_pa,
            right_absolute_pressure_pa,
        } => {
            text(&mut output, "boundary", "open");
            number(&mut output, "left_reservoir", left_absolute_pressure_pa);
            number(&mut output, "right_reservoir", right_absolute_pressure_pa);
        }
        GasFilmBoundaryTopology::Vented {
            cell_index,
            absolute_pressure_pa,
        } => {
            text(&mut output, "boundary", "vented");
            count(&mut output, "vent_cell", cell_index);
            number(&mut output, "vent_reservoir", absolute_pressure_pa);
        }
    }
    match &input.slip_policy {
        SlipPolicy::NoSlipContinuum { source_id } => text(&mut output, "slip", source_id),
        SlipPolicy::RarefiedSlipRequested { source_id } => {
            text(&mut output, "rarefied_slip", source_id)
        }
    }
    match &input.roughness_policy {
        RoughnessPolicy::ResolvedSmooth {
            source_id,
            maximum_roughness_m,
        } => {
            text(&mut output, "roughness", source_id);
            number(&mut output, "roughness_max", *maximum_roughness_m);
        }
        RoughnessPolicy::Unresolved { source_id } => {
            text(&mut output, "roughness_unresolved", source_id)
        }
    }
    number(
        &mut output,
        "mean_free_path",
        input.applicability.mean_free_path_m,
    );
    number(
        &mut output,
        "knudsen_max",
        input.applicability.maximum_knudsen_number,
    );
    number(
        &mut output,
        "slope_max",
        input.applicability.maximum_gap_slope,
    );
    number(
        &mut output,
        "sound_speed",
        input.applicability.speed_of_sound_m_per_s,
    );
    number(
        &mut output,
        "mach_max",
        input.applicability.maximum_mach_number,
    );
    number(
        &mut output,
        "uncertainty_mu",
        input.uncertainty.viscosity_relative_bound,
    );
    number(
        &mut output,
        "uncertainty_gap",
        input.uncertainty.gap_relative_bound,
    );
    number(
        &mut output,
        "uncertainty_pressure",
        input.uncertainty.pressure_relative_bound,
    );
    number(
        &mut output,
        "wall_lower",
        input.wall_motion.lower_tangential_velocity_m_per_s,
    );
    number(
        &mut output,
        "wall_upper",
        input.wall_motion.upper_tangential_velocity_m_per_s,
    );
    number(&mut output, "gap_rate", input.wall_motion.gap_rate_m_per_s);
    number(
        &mut output,
        "initial_pressure",
        input.initial_absolute_pressure_pa,
    );
    number(
        &mut output,
        "gauge_reference",
        input.gauge_reference_absolute_pressure_pa,
    );
    number(&mut output, "timestep", input.timestep_s);
    let _ = write!(output, "iterations={};", input.budget.maximum_iterations);
    number(
        &mut output,
        "residual_tolerance",
        input.budget.mass_residual_tolerance_kg_m2_s,
    );
    number(&mut output, "relaxation", input.budget.relaxation);
    output
}

fn boundary_fluxes(
    pressure: &[f64],
    gap: &[f64],
    input: &GasFilmInput,
    spacing: f64,
) -> Result<(f64, f64), GasFilmError> {
    let velocity = 0.5
        * (input.wall_motion.lower_tangential_velocity_m_per_s
            + input.wall_motion.upper_tangential_velocity_m_per_s);
    match input.boundary {
        GasFilmBoundaryTopology::Sealed | GasFilmBoundaryTopology::Vented { .. } => Ok((0.0, 0.0)),
        GasFilmBoundaryTopology::Open {
            left_absolute_pressure_pa,
            right_absolute_pressure_pa,
        } => Ok((
            face_flux(
                left_absolute_pressure_pa,
                pressure[0],
                gap[0],
                gap[0],
                0.5 * spacing,
                input.gas,
                velocity,
            )?,
            face_flux(
                pressure[pressure.len() - 1],
                right_absolute_pressure_pa,
                gap[gap.len() - 1],
                gap[gap.len() - 1],
                0.5 * spacing,
                input.gas,
                velocity,
            )?,
        )),
    }
}

fn face_fluxes(
    pressure: &[f64],
    gap: &[f64],
    input: &GasFilmInput,
    spacing: f64,
) -> Result<Vec<f64>, GasFilmError> {
    let mut flux = vec![0.0; pressure.len() + 1];
    let (left, right) = boundary_fluxes(pressure, gap, input, spacing)?;
    flux[0] = left;
    flux[pressure.len()] = right;
    let velocity = 0.5
        * (input.wall_motion.lower_tangential_velocity_m_per_s
            + input.wall_motion.upper_tangential_velocity_m_per_s);
    for face in 1..pressure.len() {
        flux[face] = face_flux(
            pressure[face - 1],
            pressure[face],
            gap[face - 1],
            gap[face],
            spacing,
            input.gas,
            velocity,
        )?;
    }
    Ok(flux)
}

/// Advance one bounded deterministic isothermal gas-film step.
///
/// A successful result is the only publication path.  Iteration-budget failure
/// returns an error and intentionally yields no partial checkpoint/receipt.
pub fn solve_isothermal_gas_film_1d(
    input: &GasFilmInput,
    checkpoint: Option<&GasFilmCheckpoint>,
) -> Result<GasFilmStep, GasFilmError> {
    let (active, spacing, rt, max_knudsen, max_slope) = input.validate()?;
    let (old_pressure, old_gap, old_step_index) = initial_state(input, active, checkpoint)?;
    let gap = &input.grid.gap_m[..active];
    let mut pressure = old_pressure.clone();
    if let GasFilmBoundaryTopology::Vented {
        cell_index,
        absolute_pressure_pa,
    } = input.boundary
    {
        pressure[cell_index] = absolute_pressure_pa;
    }
    let mass_old = old_pressure
        .iter()
        .zip(&old_gap)
        .try_fold(0.0, |total, (pressure, gap)| {
            derived(
                total + pressure * gap * spacing / rt,
                "initial_gas_mass_kg_per_m",
            )
        })?;
    let mut iterations = 0;
    let mut max_residual = 0.0_f64;
    while iterations < input.budget.maximum_iterations {
        let flux = face_fluxes(&pressure, gap, input, spacing)?;
        let mut next = pressure.clone();
        max_residual = 0.0;
        for cell in 0..active {
            if let GasFilmBoundaryTopology::Vented {
                cell_index,
                absolute_pressure_pa,
            } = input.boundary
            {
                if cell == cell_index {
                    next[cell] = absolute_pressure_pa;
                    continue;
                }
            }
            let storage_coefficient =
                derived(gap[cell] / (rt * input.timestep_s), "storage_coefficient")?;
            let rhs = derived(
                old_pressure[cell] * old_gap[cell] / (rt * input.timestep_s),
                "previous_mass_storage",
            )?;
            let residual = derived(
                storage_coefficient * pressure[cell] - rhs
                    + (flux[cell + 1] - flux[cell]) / spacing,
                "cell_mass_residual",
            )?;
            max_residual = max_residual.max(residual.abs());
            let mut diagonal = storage_coefficient;
            let mut neighbors = Vec::with_capacity(2);
            if cell > 0 {
                neighbors.push((pressure[cell - 1], gap[cell - 1], spacing));
            } else if let GasFilmBoundaryTopology::Open {
                left_absolute_pressure_pa,
                ..
            } = input.boundary
            {
                neighbors.push((left_absolute_pressure_pa, gap[cell], 0.5 * spacing));
            }
            if cell + 1 < active {
                neighbors.push((pressure[cell + 1], gap[cell + 1], spacing));
            } else if let GasFilmBoundaryTopology::Open {
                right_absolute_pressure_pa,
                ..
            } = input.boundary
            {
                neighbors.push((right_absolute_pressure_pa, gap[cell], 0.5 * spacing));
            }
            for (neighbor_pressure, neighbor_gap, distance) in neighbors {
                let p_face = 0.5 * (pressure[cell] + neighbor_pressure);
                let h_face = 0.5 * (gap[cell] + neighbor_gap);
                let k = derived(
                    h_face.powi(3) * p_face / (12.0 * input.gas.dynamic_viscosity_pa_s * rt),
                    "linearized_reynolds_conductance",
                )?;
                diagonal = derived(diagonal + k / (distance * spacing), "cell_diagonal")?;
            }
            positive(diagonal, "cell_diagonal")?;
            let proposal = derived(
                pressure[cell] - residual / diagonal,
                "pressure_iteration_proposal",
            )?;
            let updated = derived(
                pressure[cell] + input.budget.relaxation * (proposal - pressure[cell]),
                "pressure_iteration_update",
            )?;
            positive(updated, "absolute_pressure_pa")?;
            next[cell] = updated;
        }
        iterations += 1;
        pressure = next;
        if max_residual <= input.budget.mass_residual_tolerance_kg_m2_s {
            break;
        }
    }
    let flux = face_fluxes(&pressure, gap, input, spacing)?;
    let mut cell_residual = vec![0.0; active];
    max_residual = 0.0;
    for cell in 0..active {
        let storage = derived(
            (pressure[cell] * gap[cell] - old_pressure[cell] * old_gap[cell])
                / (rt * input.timestep_s),
            "final_storage_residual",
        )?;
        let residual = derived(
            storage + (flux[cell + 1] - flux[cell]) / spacing,
            "final_cell_mass_residual",
        )?;
        cell_residual[cell] = residual;
        if !matches!(input.boundary, GasFilmBoundaryTopology::Vented { cell_index, .. } if cell == cell_index)
        {
            max_residual = max_residual.max(residual.abs());
        }
    }
    if max_residual > input.budget.mass_residual_tolerance_kg_m2_s {
        return Err(GasFilmError::IterationBudgetExceeded {
            iterations,
            max_mass_residual_kg_m2_s: max_residual,
        });
    }
    let mass_new = pressure
        .iter()
        .zip(gap)
        .try_fold(0.0, |total, (pressure, gap)| {
            derived(total + pressure * gap * spacing / rt, "gas_mass_kg_per_m")
        })?;
    let storage_rate = derived(
        (mass_new - mass_old) / input.timestep_s,
        "storage_rate_kg_per_m_s",
    )?;
    let left_outward = derived(-flux[0], "left_boundary_outward_mass_flux")?;
    let right_outward = flux[active];
    let vent_outward = match input.boundary {
        GasFilmBoundaryTopology::Vented { cell_index, .. } => derived(
            -cell_residual[cell_index] * spacing,
            "vent_outward_mass_flux",
        )?,
        _ => 0.0,
    };
    let left_enthalpy = derived(
        left_outward * input.gas.declared_specific_enthalpy_j_kg,
        "left_boundary_outward_enthalpy_flux",
    )?;
    let right_enthalpy = derived(
        right_outward * input.gas.declared_specific_enthalpy_j_kg,
        "right_boundary_outward_enthalpy_flux",
    )?;
    let vent_enthalpy = derived(
        vent_outward * input.gas.declared_specific_enthalpy_j_kg,
        "vent_outward_enthalpy_flux",
    )?;
    let closure = derived(
        storage_rate + left_outward + right_outward + vent_outward,
        "mass_closure_residual_kg_per_m_s",
    )?;
    let mut absolute_load = 0.0;
    let mut gauge_load = 0.0;
    let mut gas_on_upper_normal = Vec::with_capacity(input.grid.gap_m.len());
    let mut upper_on_gas_normal = Vec::with_capacity(input.grid.gap_m.len());
    let mut gas_on_upper_tangential = Vec::with_capacity(input.grid.gap_m.len());
    let mut upper_on_gas_tangential = Vec::with_capacity(input.grid.gap_m.len());
    let mut upper_shear = Vec::with_capacity(input.grid.gap_m.len());
    let mut viscous_power = 0.0;
    let mut normal_gap_power = 0.0;
    let relative_velocity = input.wall_motion.upper_tangential_velocity_m_per_s
        - input.wall_motion.lower_tangential_velocity_m_per_s;
    for (index, gap) in input.grid.gap_m.iter().enumerate() {
        if index >= active {
            gas_on_upper_normal.push(None);
            upper_on_gas_normal.push(None);
            gas_on_upper_tangential.push(None);
            upper_on_gas_tangential.push(None);
            upper_shear.push(None);
            continue;
        }
        let pressure_value = pressure[index];
        absolute_load = derived(
            absolute_load + pressure_value * spacing,
            "absolute_pressure_load_n_per_m",
        )?;
        gauge_load = derived(
            gauge_load + (pressure_value - input.gauge_reference_absolute_pressure_pa) * spacing,
            "gauge_pressure_load_n_per_m",
        )?;
        normal_gap_power = derived(
            normal_gap_power - pressure_value * input.wall_motion.gap_rate_m_per_s * spacing,
            "normal_gap_power_to_gas_w_per_m",
        )?;
        let shear = derived(
            input.gas.dynamic_viscosity_pa_s * relative_velocity / gap,
            "upper_wall_shear_pa",
        )?;
        gas_on_upper_normal.push(Some(pressure_value));
        upper_on_gas_normal.push(Some(-pressure_value));
        gas_on_upper_tangential.push(Some(-shear));
        upper_on_gas_tangential.push(Some(shear));
        upper_shear.push(Some(shear));
        viscous_power = derived(
            viscous_power + shear * relative_velocity * spacing,
            "viscous_heat_w_per_m",
        )?;
    }
    if viscous_power < -1.0e-18 {
        return Err(GasFilmError::NonFiniteDerived {
            field: "negative_viscous_dissipation",
        });
    }
    let wall_power = derived(
        viscous_power + normal_gap_power,
        "total_wall_power_to_gas_w_per_m",
    )?;
    let checkpoint_pressure = pressure.clone();
    let active_field = pressure.into_iter().map(Some).collect::<Vec<_>>();
    let mut field = active_field.clone();
    field.resize(input.grid.gap_m.len(), None);
    Ok(GasFilmStep {
        absolute_pressure_pa: field,
        receipt: GasFilmReceipt {
            model_id: input.identity.model_id.clone(),
            input_authority: input.identity.authority,
            input_uncertainty: input.uncertainty,
            iterations,
            max_mass_residual_kg_m2_s: max_residual,
            gas_mass_kg_per_m: mass_new,
            storage_rate_kg_per_m_s: storage_rate,
            left_boundary_outward_mass_flux_kg_per_m_s: left_outward,
            right_boundary_outward_mass_flux_kg_per_m_s: right_outward,
            vent_outward_mass_flux_kg_per_m_s: vent_outward,
            left_boundary_outward_enthalpy_flux_w_per_m: left_enthalpy,
            right_boundary_outward_enthalpy_flux_w_per_m: right_enthalpy,
            vent_outward_enthalpy_flux_w_per_m: vent_enthalpy,
            mass_closure_residual_kg_per_m_s: closure,
            absolute_pressure_load_n_per_m: absolute_load,
            gauge_pressure_load_n_per_m: gauge_load,
            gas_on_upper_wall_normal_traction_pa: gas_on_upper_normal,
            upper_wall_on_gas_normal_traction_pa: upper_on_gas_normal,
            gas_on_upper_wall_tangential_traction_pa: gas_on_upper_tangential,
            upper_wall_on_gas_tangential_traction_pa: upper_on_gas_tangential,
            upper_wall_shear_pa: upper_shear,
            normal_gap_power_to_gas_w_per_m: normal_gap_power,
            wall_power_to_gas_w_per_m: wall_power,
            viscous_heat_w_per_m: viscous_power,
            maximum_knudsen_number: max_knudsen,
            maximum_gap_slope: max_slope,
        },
        checkpoint: GasFilmCheckpoint {
            step_index: old_step_index
                .checked_add(1)
                .ok_or(GasFilmError::NonFiniteDerived {
                    field: "checkpoint_step_index",
                })?,
            model_id: input.identity.model_id.clone(),
            case_id: input.identity.case_id.clone(),
            input_authority: input.identity.authority,
            input_uncertainty: input.uncertainty,
            configuration_fingerprint: invariant_configuration_fingerprint(input),
            accepted_request_fingerprint: accepted_request_fingerprint(input),
            active_absolute_pressure_pa: checkpoint_pressure,
            active_gap_m: gap.to_vec(),
            active_cells: active,
        },
    })
}

/// Canonical model identity for callers constructing a [`GasFilmIdentity`].
pub const fn isothermal_compressible_reynolds_model_id() -> &'static str {
    MODEL_ID
}
