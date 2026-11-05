//! Bounded one-dimensional isothermal compressible gas-film Reynolds solver.
//!
//! This is deliberately a small, generic continuum foundation: a uniform
//! Cartesian finite-volume line of unit out-of-plane width, or a radial
//! axisymmetric disk/annulus with exact annular storage. It solves the
//! isothermal ideal-gas Reynolds mass balance, not a liquid cavitation or EHL
//! model.  Every pressure is absolute [Pa], every gap is explicit [m], and an
//! excluded contact cell is absent from the gas state rather than receiving an
//! anonymous lower-gap floor.

use core::fmt;
use core::fmt::Write;
use fs_la::tridiagonal::{TridiagonalError, TridiagonalWorkspace};

mod axisymmetric;
mod capillary;
pub use axisymmetric::{
    AxisymmetricCapillaryVent, AxisymmetricCapillaryVentReport, AxisymmetricGasFilmReceipt,
    AxisymmetricGasFilmReservoir, AxisymmetricGasFilmReservoirReport, AxisymmetricGasFilmStep,
    solve_isothermal_axisymmetric_gas_film_1d,
};
use capillary::{CapillaryConnection, CapillaryExchange};
pub use capillary::{CapillaryVent, CapillaryVentReport, GasFilmReservoir, GasFilmReservoirReport};

const MODEL_ID: &str = "fs-flux/isothermal-compressible-reynolds-1d-v1";
const MAX_GEOMETRIC_CELLS: usize = 4_096;
const MAX_NEWTON_BACKTRACKS: usize = 32;
const NEWTON_SUFFICIENT_DECREASE: f64 = 1.0e-4;
const NEWTON_POSITIVE_FRACTION: f64 = 0.99;
// Classical first-order Maxwell slip is not a transition/molecular-flow model.
// See Gavasane et al., doi:10.1260/1759-3093.2.2-3.109.
const MAX_MAXWELL_SLIP_KNUDSEN_NUMBER: f64 = 0.1;
const RESTART_GAP_EVOLUTION_ABSOLUTE_TOLERANCE_M: f64 = 1.0e-15;
const RESTART_GAP_EVOLUTION_RELATIVE_TOLERANCE: f64 = 1.0e-12;

/// Typed refusal or bounded-solve failure for the gas-film foundation.
#[derive(Debug, Clone, PartialEq)]
pub enum GasFilmError {
    /// A named scalar, vector entry, or count is not admitted.
    InvalidInput {
        /// Name of invalid input field.
        field: &'static str,
    },
    /// A semantic identity is empty or contains non-canonical bytes.
    InvalidIdentity {
        /// Name of invalid identity field.
        field: &'static str,
    },
    /// A requested model lies outside this deliberately narrow continuum scope.
    Unavailable {
        /// Reason for unavailability.
        reason: &'static str,
    },
    /// The static mask would split or change the retained one-dimensional topology.
    TopologyChangeUnavailable,
    /// A checkpoint is not compatible with the next requested step.
    CheckpointMismatch {
        /// Mismatched field name.
        field: &'static str,
    },
    /// A derived quantity overflowed or became non-finite.
    NonFiniteDerived {
        /// Non-finite derived field name.
        field: &'static str,
    },
    /// The deterministic Newton iteration budget ended before the residual gate.
    IterationBudgetExceeded {
        /// Number of completed iterations.
        iterations: u32,
        /// Largest retained cell mass residual [kg m^-2 s^-1].
        max_mass_residual_kg_m2_s: f64,
    },
    /// The native pivoted tridiagonal pressure Jacobian could not be solved.
    LinearSolveFailure {
        /// Native linear-algebra refusal; no partial gas-film state is published.
        source: TridiagonalError,
    },
    /// No positive-pressure residual-decreasing trial was found in the bounded
    /// Newton line search. No pressure floor or unconverged state is substituted.
    NewtonStepRejected {
        /// Number of attempted Newton corrections.
        iterations: u32,
        /// Last accepted nonlinear mass residual [kg m^-2 s^-1].
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
    /// Declared frame for the positive line coordinate and wall speeds, or the
    /// symmetry axis and positive radial direction for an axisymmetric film.
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
    /// Independently sealed or pressure-open ends. `None` means zero mass flux;
    /// `Some(p)` supplies an absolute reservoir pressure [Pa]. In radial
    /// geometry, left is the inner end and right is the outer end.
    MixedEnds {
        /// Left/inner pressure bath, or a sealed end.
        left_absolute_pressure_pa: Option<f64>,
        /// Right/outer pressure bath, or a sealed end.
        right_absolute_pressure_pa: Option<f64>,
    },
    /// Sealed outer ends with one ideal, zero-impedance pressure vent.
    /// The selected cell is held at the bath pressure with unlimited exchange.
    Vented {
        /// Active-cell index of the vent.
        cell_index: usize,
        /// Vent reservoir pressure [Pa absolute].
        absolute_pressure_pa: f64,
    },
    /// Sealed outer ends with finite-conductance capillaries to a pressure bath.
    /// The cell pressure is solved from storage, film flux and capillary flow;
    /// it is not pinned to the exterior pressure.
    CapillaryVented {
        /// Active-cell index of the capillary connection.
        cell_index: usize,
        /// Isothermal bath pressure [Pa absolute], allowed to vary by step.
        absolute_pressure_pa: f64,
        /// Explicit stationary no-slip capillary geometry.
        capillary: CapillaryVent,
    },
    /// Sealed outer ends exchanging gas with a finite rigid isothermal reservoir.
    /// Reservoir pressure evolves implicitly with the film through the capillary.
    CapillaryReservoir {
        /// Active-cell index of the capillary connection.
        cell_index: usize,
        /// Explicit stationary no-slip capillary geometry.
        capillary: CapillaryVent,
        /// Finite gas inventory initialized from volume and absolute pressure.
        reservoir: GasFilmReservoir,
    },
    /// Sealed radial ends with a total-count capillary bundle to a pressure bath.
    /// Admitted only by the axisymmetric entry point; exchange is homogenized
    /// over the selected annular cell rather than pinned at the symmetry axis.
    AxisymmetricCapillaryVented {
        /// Active radial-cell index of the connection.
        cell_index: usize,
        /// Isothermal bath pressure [Pa absolute], allowed to vary by step.
        absolute_pressure_pa: f64,
        /// Explicit total-count geometry of the stationary no-slip tubes.
        capillary: AxisymmetricCapillaryVent,
    },
    /// Sealed radial ends with a finite total-volume isothermal reservoir.
    /// Admitted only by the axisymmetric entry point. Film and reservoir
    /// exchange the same total mass through the implicit capillary connection.
    AxisymmetricCapillaryReservoir {
        /// Active radial-cell index of the connection.
        cell_index: usize,
        /// Explicit total-count geometry of the stationary no-slip tubes.
        capillary: AxisymmetricCapillaryVent,
        /// Finite total gas volume and initial absolute pressure.
        reservoir: AxisymmetricGasFilmReservoir,
    },
}

impl GasFilmBoundaryTopology {
    fn is_axisymmetric_capillary(self) -> bool {
        matches!(
            self,
            Self::AxisymmetricCapillaryVented { .. } | Self::AxisymmetricCapillaryReservoir { .. }
        )
    }

    fn end_pressures(self) -> [Option<f64>; 2] {
        match self {
            Self::Open {
                left_absolute_pressure_pa,
                right_absolute_pressure_pa,
            } => [
                Some(left_absolute_pressure_pa),
                Some(right_absolute_pressure_pa),
            ],
            Self::MixedEnds {
                left_absolute_pressure_pa,
                right_absolute_pressure_pa,
            } => [left_absolute_pressure_pa, right_absolute_pressure_pa],
            _ => [None, None],
        }
    }
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
    /// Total line length [m], or radial span for the axisymmetric entry point.
    pub length_m: f64,
    /// Current positive gap at each geometric cell [m].
    pub gap_m: Vec<f64>,
    /// Static gas/contact partition.
    pub contact_exclusion: ContactExclusionMask,
}

/// Wall motion used by Reynolds mass flux and signed mechanical work.
#[derive(Debug, Clone, PartialEq)]
pub struct MovingWallInput {
    /// Lower-wall line velocity [m s^-1].
    pub lower_tangential_velocity_m_per_s: f64,
    /// Upper-wall line velocity [m s^-1].
    pub upper_tangential_velocity_m_per_s: f64,
    /// Uniform fallback `dh/dt` [m s^-1]; negative values squeeze the gas.
    ///
    /// This value is used when `gap_rate_profile_m_per_s` is `None`.
    pub gap_rate_m_per_s: f64,
    /// Optional exact per-cell `dh/dt` [m s^-1].
    ///
    /// A profile is required for changing wedge slopes, such as a tilting
    /// rigid surface. When present it has exactly one value per geometric cell
    /// and overrides `gap_rate_m_per_s` in storage/restart/work accounting.
    pub gap_rate_profile_m_per_s: Option<Vec<f64>>,
}

/// Explicit continuum/slip-flow admission gates for the selected wall policy.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GasFilmApplicability {
    /// Caller-supplied molecular mean free path [m] at the declared temperature
    /// and `GasFilmInput::initial_absolute_pressure_pa`.
    /// The isothermal gas uses `lambda(p) = lambda_ref * p_ref / p`.
    pub mean_free_path_m: f64,
    /// Maximum admitted cell or reservoir/vent Knudsen number `lambda(p) / h`.
    pub maximum_knudsen_number: f64,
    /// Maximum admitted adjacent-cell slope `abs(dh/dx)`.
    pub maximum_gap_slope: f64,
    /// Caller-declared sound speed [m s^-1] for the isothermal gas state.
    pub speed_of_sound_m_per_s: f64,
    /// Maximum admitted tangential gas Mach number, including wall velocities
    /// and the interior peak of each Couette--Poiseuille velocity profile.
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

/// Explicit smooth-wall boundary condition for the isothermal gas film.
#[derive(Debug, Clone, PartialEq)]
pub enum SlipPolicy {
    /// Named no-slip continuum policy.
    NoSlipContinuum {
        /// Source identity of the policy.
        source_id: String,
    },
    /// Classical first-order Maxwell velocity slip on both smooth walls.
    ///
    /// Both walls use the same declared tangential momentum accommodation
    /// coefficient `sigma` in `(0, 1]`, with slip length
    /// `b(p) = (2-sigma)/sigma * lambda_ref * p_ref/p`. The caller's Knudsen
    /// cap still applies, and local `Kn > 0.1` is always refused. Temperature
    /// is uniform: there is no thermal creep, temperature jump, or fitted
    /// higher-order/Knudsen-layer correction. This is not a molecular model.
    MaxwellFirstOrder {
        /// Identity of the caller's gas/surface accommodation declaration.
        source_id: String,
        /// Shared wall tangential momentum accommodation coefficient `(0, 1]`.
        tangential_momentum_accommodation: f64,
    },
    /// An unspecified rarefied/slip model is not an explicit Maxwell admission.
    RarefiedSlipRequested {
        /// Source identity of requested slip model.
        source_id: String,
    },
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
    Unresolved {
        /// Source identity of unresolved roughness.
        source_id: String,
    },
}

/// Bounded nonlinear iteration and numerical acceptance controls.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GasFilmBudget {
    /// Maximum deterministic Newton corrections; zero is refused.
    pub maximum_iterations: u32,
    /// Maximum cell mass-balance residual [kg m^-2 s^-1].
    pub mass_residual_tolerance_kg_m2_s: f64,
    /// Maximum Newton step fraction in `(0, 1]`. Bounded backtracking may reduce
    /// it further to preserve positive pressure and decrease the mass residual.
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
    /// Initial absolute pressure for a fresh state [Pa], also the fixed
    /// reference pressure of the declared molecular mean free path.
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
    /// Accepted finite reservoir mass [kg m^-1], present only for
    /// Cartesian `CapillaryReservoir`. Absent for every radial checkpoint.
    /// Volume and gas identity are fixed by the fingerprint.
    pub capillary_reservoir_mass_kg_per_m: Option<f64>,
    /// Accepted total reservoir mass [kg], present only for
    /// `AxisymmetricCapillaryReservoir`. Absent for every Cartesian checkpoint;
    /// the total volume, tube count and radial geometry bind its identity.
    pub capillary_reservoir_mass_kg: Option<f64>,
}

/// Pressure/shear/load/work and conservative accounting receipt.
#[derive(Debug, Clone, PartialEq)]
pub struct GasFilmReceipt {
    /// Model identity retained from the request.
    pub model_id: String,
    /// Input authority retained without promotion.
    pub input_authority: GasFilmInputAuthority,
    /// Number of deterministic Newton corrections used; zero at initial equilibrium.
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
    /// Additional mechanics of a finite capillary connection, when present.
    /// Its dissipation is outside the film and is not included in
    /// `viscous_heat_w_per_m` or counted again as boundary enthalpy.
    pub capillary_vent: Option<CapillaryVentReport>,
    /// Finite connected gas storage, when present. The capillary's signed
    /// exchange is internal to the combined film/reservoir system.
    pub finite_reservoir: Option<GasFilmReservoirReport>,
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
    /// Per-cell Couette-plus-Poiseuille tangential traction of gas on the upper wall [Pa].
    pub gas_on_upper_wall_tangential_traction_pa: Vec<Option<f64>>,
    /// Per-cell tangential traction of upper wall on gas [Pa]; exactly the negative of gas-on-wall.
    pub upper_wall_on_gas_tangential_traction_pa: Vec<Option<f64>>,
    /// Per-cell Couette-plus-Poiseuille tangential traction of gas on the lower wall [Pa].
    pub gas_on_lower_wall_tangential_traction_pa: Vec<Option<f64>>,
    /// Per-cell tangential traction of lower wall on gas [Pa]; exactly the negative of gas-on-wall.
    pub lower_wall_on_gas_tangential_traction_pa: Vec<Option<f64>>,
    /// Legacy-named alias for upper-wall-on-gas tangential traction [Pa].
    /// Excluded cells are `None`; new callers should use the explicit field above.
    pub upper_wall_shear_pa: Vec<Option<f64>>,
    /// Normal moving-gap mechanical power into gas, `-integral(p dh/dt dx)` [W m^-1].
    pub normal_gap_power_to_gas_w_per_m: f64,
    /// Signed work by both translating walls plus normal moving-gap work [W m^-1].
    ///
    /// Pressure-driven gas may do work on a moving wall, so this can be negative.
    pub wall_power_to_gas_w_per_m: f64,
    /// Pressure-gradient force work, `-integral(q dp/dx dx)` [W m^-1].
    ///
    /// Here `q` is the accepted mass flux divided by the face ideal-gas density.
    /// This is not boundary enthalpy transport or a total-energy closure.
    pub pressure_gradient_power_to_gas_w_per_m: f64,
    /// Nonnegative Couette-plus-Poiseuille dissipation [W m^-1], including
    /// gas/wall slip friction when the explicit Maxwell policy is selected.
    ///
    /// Equals tangential wall work plus pressure-gradient force work up to
    /// roundoff; normal gap work is separate. No temperature evolution is solved.
    pub viscous_heat_w_per_m: f64,
    /// Part of `viscous_heat_w_per_m` dissipated at the slipping gas/wall
    /// interfaces [W m^-1]; exactly zero for no-slip. The remainder is bulk
    /// viscous dissipation. No thermal partition between wall and gas is assumed.
    pub wall_slip_dissipation_w_per_m: f64,
    /// Maximum pressure-scaled Knudsen number in the film and directly imposed
    /// pressure boundaries. A finite capillary has its separate radius-based value.
    pub maximum_knudsen_number: f64,
    /// Maximum tangential gas Mach number across the accepted face profiles,
    /// including their wall values and any interior velocity extremum.
    pub maximum_gas_mach_number: f64,
    /// Maximum adjacent gap slope.
    pub maximum_gap_slope: f64,
    /// Caller-declared input uncertainty retained without propagation or promotion.
    pub input_uncertainty: GasFilmUncertainty,
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
    if mask[active..].iter().any(|excluded| !*excluded) {
        return Err(GasFilmError::TopologyChangeUnavailable);
    }
    if active < 2 {
        return Err(GasFilmError::Unavailable {
            reason: "fewer-than-two-active-gas-cells",
        });
    }
    Ok(active)
}

#[derive(Clone, Copy)]
struct FilmFaceMetric {
    // Unit width for a Cartesian strip, circular circumference for a radial film.
    width: f64,
    // Radial faces use the integrated logarithmic resistance, not delta-r.
    distance_m: f64,
}

// Geometry is compiled once. The pressure/slip law, mass balance and mechanical
// quadrature below consume the same measures; there is no second radial solver.
enum FilmMetrics {
    Cartesian {
        spacing_m: f64,
    },
    Axisymmetric {
        inner_radius_m: f64,
        cell_area_m2: Vec<f64>,
        faces: Vec<FilmFaceMetric>,
    },
}

impl FilmMetrics {
    fn cell_measure(&self, cell: usize) -> f64 {
        match self {
            Self::Cartesian { spacing_m } => *spacing_m,
            Self::Axisymmetric { cell_area_m2, .. } => cell_area_m2[cell],
        }
    }

    fn face(&self, face: usize, active: usize) -> FilmFaceMetric {
        match self {
            Self::Cartesian { spacing_m } => FilmFaceMetric {
                width: 1.0,
                distance_m: if face == 0 || face == active {
                    0.5 * spacing_m
                } else {
                    *spacing_m
                },
            },
            Self::Axisymmetric { faces, .. } => faces[face],
        }
    }

    fn bind_fingerprint(&self, mut fingerprint: String) -> String {
        if let Self::Axisymmetric { inner_radius_m, .. } = self {
            let _ = write!(
                fingerprint,
                "geometry=axisymmetric-log-resistance-v1;inner_radius={:016x};",
                inner_radius_m.to_bits()
            );
        }
        fingerprint
    }
}

// Compile the wall law once per step. For Maxwell slip, b*p is constant at
// fixed temperature, so the extra MASS conductance is independent of pressure.
// No-slip keeps the original arithmetic and requires no added allocation.
struct FilmTransport {
    viscosity_pa_s: f64,
    rt: f64,
    mean_wall_velocity_m_per_s: f64,
    slip_pressure_length_pa_m: f64,
}

impl FilmTransport {
    fn new(input: &GasFilmInput, rt: f64) -> Result<Self, GasFilmError> {
        let slip_pressure_length_pa_m = match &input.slip_policy {
            SlipPolicy::NoSlipContinuum { .. } => 0.0,
            SlipPolicy::MaxwellFirstOrder {
                tangential_momentum_accommodation: sigma,
                ..
            } => derived(
                (2.0 - sigma) / sigma
                    * input.applicability.mean_free_path_m
                    * input.initial_absolute_pressure_pa,
                "maxwell_slip_pressure_length_pa_m",
            )?,
            SlipPolicy::RarefiedSlipRequested { .. } => {
                return Err(GasFilmError::Unavailable {
                    reason: "slip-or-rarefied-model-not-implemented",
                });
            }
        };
        Ok(Self {
            viscosity_pa_s: input.gas.dynamic_viscosity_pa_s,
            rt,
            mean_wall_velocity_m_per_s: 0.5
                * (input.wall_motion.lower_tangential_velocity_m_per_s
                    + input.wall_motion.upper_tangential_velocity_m_per_s),
            slip_pressure_length_pa_m,
        })
    }

    fn slip_length_m(&self, pressure_pa: f64) -> Result<f64, GasFilmError> {
        derived(
            self.slip_pressure_length_pa_m / pressure_pa,
            "wall_slip_length_m",
        )
    }

    fn mass_conductance(&self, pressure_pa: f64, gap_m: f64) -> Result<f64, GasFilmError> {
        let no_slip = derived(
            gap_m.powi(3) * pressure_pa / (12.0 * self.viscosity_pa_s * self.rt),
            "reynolds_mass_conductance",
        )?;
        if self.slip_pressure_length_pa_m == 0.0 {
            return Ok(no_slip);
        }
        // rho*h^3*(1+6*b/h)/(12*mu), with b = (b*p)/p.
        let slip =
            gap_m.powi(2) * self.slip_pressure_length_pa_m / (2.0 * self.viscosity_pa_s * self.rt);
        derived(no_slip + slip, "maxwell_reynolds_mass_conductance")
    }
}

fn face_flux(
    left_pressure_pa: f64,
    right_pressure_pa: f64,
    left_gap_m: f64,
    right_gap_m: f64,
    distance_m: f64,
    transport: &FilmTransport,
) -> Result<f64, GasFilmError> {
    let pressure = derived(
        0.5 * (left_pressure_pa + right_pressure_pa),
        "face_pressure_pa",
    )?;
    let gap = derived(0.5 * (left_gap_m + right_gap_m), "face_gap_m")?;
    positive(pressure, "face_pressure_pa")?;
    positive(gap, "face_gap_m")?;
    let conductance = transport.mass_conductance(pressure, gap)?;
    let pressure_flux = derived(
        -conductance * (right_pressure_pa - left_pressure_pa) / distance_m,
        "pressure_mass_flux",
    )?;
    let couette_flux = derived(
        pressure * gap * transport.mean_wall_velocity_m_per_s / transport.rt,
        "couette_mass_flux",
    )?;
    let flux = derived(pressure_flux + couette_flux, "face_mass_flux")?;
    Ok(flux)
}

impl GasFilmInput {
    fn gap_rate_m_per_s(&self, cell: usize) -> f64 {
        self.wall_motion
            .gap_rate_profile_m_per_s
            .as_ref()
            .map_or(self.wall_motion.gap_rate_m_per_s, |profile| profile[cell])
    }

    fn admit_knudsen(&self, pressure_pa: f64, gap_m: f64) -> Result<f64, GasFilmError> {
        let knudsen = derived(
            (self.applicability.mean_free_path_m / gap_m)
                * (self.initial_absolute_pressure_pa / pressure_pa),
            "pressure_scaled_knudsen_number",
        )?;
        if matches!(self.slip_policy, SlipPolicy::MaxwellFirstOrder { .. })
            && knudsen > MAX_MAXWELL_SLIP_KNUDSEN_NUMBER
        {
            return Err(GasFilmError::Unavailable {
                reason: "maxwell-slip-knudsen-outside-first-order-envelope",
            });
        }
        if knudsen > self.applicability.maximum_knudsen_number {
            return Err(GasFilmError::Unavailable {
                reason: "rarefied-knudsen-outside-continuum-envelope",
            });
        }
        Ok(knudsen)
    }

    fn validate(&self) -> Result<(usize, f64, f64, f64), GasFilmError> {
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
                canonical_identity(source_id, "slip_source_id")?;
            }
            SlipPolicy::MaxwellFirstOrder {
                source_id,
                tangential_momentum_accommodation: sigma,
            } => {
                canonical_identity(source_id, "slip_source_id")?;
                if !(sigma.is_finite() && *sigma > 0.0 && *sigma <= 1.0) {
                    return Err(GasFilmError::InvalidInput {
                        field: "slip_policy.tangential_momentum_accommodation",
                    });
                }
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
        if let Some(profile) = &self.wall_motion.gap_rate_profile_m_per_s {
            if profile.len() != self.grid.gap_m.len() {
                return Err(GasFilmError::InvalidInput {
                    field: "wall_motion.gap_rate_profile_m_per_s.length",
                });
            }
            for rate in profile {
                finite(*rate, "wall_motion.gap_rate_profile_m_per_s")?;
            }
        }
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
        let mut max_slope = 0.0_f64;
        for (index, gap) in self.grid.gap_m.iter().take(active).copied().enumerate() {
            positive(gap, "grid.gap_m")?;
            if let RoughnessPolicy::ResolvedSmooth {
                maximum_roughness_m,
                ..
            } = &self.roughness_policy
                && *maximum_roughness_m >= gap
            {
                return Err(GasFilmError::Unavailable {
                    reason: "roughness-not-small-relative-to-gap",
                });
            }
            if index > 0 {
                let slope = derived(
                    (gap - self.grid.gap_m[index - 1]).abs() / spacing,
                    "gap_slope",
                )?;
                max_slope = max_slope.max(slope);
            }
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
            GasFilmBoundaryTopology::MixedEnds {
                left_absolute_pressure_pa,
                right_absolute_pressure_pa,
            } => {
                if let Some(pressure) = left_absolute_pressure_pa {
                    positive(pressure, "boundary.left_absolute_pressure_pa")?;
                }
                if let Some(pressure) = right_absolute_pressure_pa {
                    positive(pressure, "boundary.right_absolute_pressure_pa")?;
                }
            }
            GasFilmBoundaryTopology::Vented {
                cell_index,
                absolute_pressure_pa,
            }
            | GasFilmBoundaryTopology::CapillaryVented {
                cell_index,
                absolute_pressure_pa,
                ..
            }
            | GasFilmBoundaryTopology::AxisymmetricCapillaryVented {
                cell_index,
                absolute_pressure_pa,
                ..
            } => {
                if cell_index >= active {
                    return Err(GasFilmError::InvalidInput {
                        field: "boundary.vent_cell_index",
                    });
                }
                positive(absolute_pressure_pa, "boundary.vent_absolute_pressure_pa")?;
            }
            GasFilmBoundaryTopology::CapillaryReservoir {
                cell_index,
                reservoir,
                ..
            } => {
                if cell_index >= active {
                    return Err(GasFilmError::InvalidInput {
                        field: "boundary.vent_cell_index",
                    });
                }
                positive(
                    reservoir.volume_per_width_m2,
                    "reservoir.volume_per_width_m2",
                )?;
                positive(
                    reservoir.initial_absolute_pressure_pa,
                    "reservoir.initial_absolute_pressure_pa",
                )?;
            }
            GasFilmBoundaryTopology::AxisymmetricCapillaryReservoir {
                cell_index,
                reservoir,
                ..
            } => {
                if cell_index >= active {
                    return Err(GasFilmError::InvalidInput {
                        field: "boundary.vent_cell_index",
                    });
                }
                positive(reservoir.volume_m3, "reservoir.volume_m3")?;
                positive(
                    reservoir.initial_absolute_pressure_pa,
                    "reservoir.initial_absolute_pressure_pa",
                )?;
            }
        }
        if let GasFilmBoundaryTopology::CapillaryVented { capillary, .. }
        | GasFilmBoundaryTopology::CapillaryReservoir { capillary, .. } = self.boundary
        {
            capillary.validate()?;
        }
        if let GasFilmBoundaryTopology::AxisymmetricCapillaryVented { capillary, .. }
        | GasFilmBoundaryTopology::AxisymmetricCapillaryReservoir { capillary, .. } =
            self.boundary
        {
            capillary.validate()?;
        }
        Ok((active, spacing, rt, max_slope))
    }
}

fn initial_state(
    input: &GasFilmInput,
    active: usize,
    checkpoint: Option<&GasFilmCheckpoint>,
    configuration_fingerprint: &str,
) -> Result<(Vec<f64>, Vec<f64>, u64), GasFilmError> {
    match checkpoint {
        None => {
            let old_gap = input
                .grid
                .gap_m
                .iter()
                .take(active)
                .enumerate()
                .map(|(index, gap)| gap - input.gap_rate_m_per_s(index) * input.timestep_s)
                .collect::<Vec<_>>();
            for gap in &old_gap {
                positive(*gap, "derived.previous_gap_m")?;
            }
            Ok((vec![input.initial_absolute_pressure_pa; active], old_gap, 0))
        }
        Some(state) => {
            if matches!(
                input.boundary,
                GasFilmBoundaryTopology::CapillaryReservoir { .. }
            ) != state.capillary_reservoir_mass_kg_per_m.is_some()
            {
                return Err(GasFilmError::CheckpointMismatch {
                    field: "capillary_reservoir_mass_kg_per_m",
                });
            }
            if matches!(
                input.boundary,
                GasFilmBoundaryTopology::AxisymmetricCapillaryReservoir { .. }
            ) != state.capillary_reservoir_mass_kg.is_some()
            {
                return Err(GasFilmError::CheckpointMismatch {
                    field: "capillary_reservoir_mass_kg",
                });
            }
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
            if state.configuration_fingerprint != configuration_fingerprint {
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
            for (index, (new_gap, old_gap)) in input
                .grid
                .gap_m
                .iter()
                .take(active)
                .zip(&state.active_gap_m)
                .enumerate()
            {
                let expected = derived(
                    *old_gap + input.gap_rate_m_per_s(index) * input.timestep_s,
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
        GasFilmBoundaryTopology::MixedEnds {
            left_absolute_pressure_pa,
            right_absolute_pressure_pa,
        } => {
            text(&mut output, "boundary_kind", "mixed-ends");
            let _ = write!(
                output,
                "left_open={};right_open={};",
                left_absolute_pressure_pa.is_some(),
                right_absolute_pressure_pa.is_some()
            );
        }
        GasFilmBoundaryTopology::Vented { .. } => text(&mut output, "boundary_kind", "vented"),
        GasFilmBoundaryTopology::CapillaryVented {
            cell_index,
            capillary,
            ..
        }
        | GasFilmBoundaryTopology::CapillaryReservoir {
            cell_index,
            capillary,
            ..
        } => {
            if let GasFilmBoundaryTopology::CapillaryReservoir { reservoir, .. } = input.boundary {
                text(
                    &mut output,
                    "boundary_kind",
                    "finite-isothermal-capillary-reservoir-v1",
                );
                number(
                    &mut output,
                    "reservoir_volume_per_width",
                    reservoir.volume_per_width_m2,
                );
                number(
                    &mut output,
                    "reservoir_initial_pressure",
                    reservoir.initial_absolute_pressure_pa,
                );
            } else {
                text(&mut output, "boundary_kind", "no-slip-capillary-vent-v1");
            }
            let _ = write!(output, "vent_cell={cell_index};");
            number(&mut output, "capillary_radius", capillary.radius_m);
            number(&mut output, "capillary_length", capillary.length_m);
            number(
                &mut output,
                "capillaries_per_m",
                capillary.capillaries_per_m,
            );
        }
        GasFilmBoundaryTopology::AxisymmetricCapillaryVented {
            cell_index,
            capillary,
            ..
        }
        | GasFilmBoundaryTopology::AxisymmetricCapillaryReservoir {
            cell_index,
            capillary,
            ..
        } => {
            if let GasFilmBoundaryTopology::AxisymmetricCapillaryReservoir { reservoir, .. } =
                input.boundary
            {
                text(
                    &mut output,
                    "boundary_kind",
                    "radial-capillary-reservoir-v1",
                );
                number(&mut output, "reservoir_volume_m3", reservoir.volume_m3);
                number(
                    &mut output,
                    "reservoir_initial_pressure",
                    reservoir.initial_absolute_pressure_pa,
                );
            } else {
                text(&mut output, "boundary_kind", "radial-capillary-vent-v1");
            }
            let _ = write!(output, "vent_cell={cell_index};");
            number(&mut output, "capillary_radius", capillary.radius_m);
            number(&mut output, "capillary_length", capillary.length_m);
            let _ = write!(output, "capillary_count={};", capillary.capillary_count);
        }
    }
    match &input.slip_policy {
        SlipPolicy::NoSlipContinuum { source_id } => text(&mut output, "slip", source_id),
        SlipPolicy::MaxwellFirstOrder {
            source_id,
            tangential_momentum_accommodation,
        } => {
            text(&mut output, "maxwell_first_order", source_id);
            number(
                &mut output,
                "maxwell_tmac",
                *tangential_momentum_accommodation,
            );
        }
        SlipPolicy::RarefiedSlipRequested { source_id } => {
            text(&mut output, "rarefied_slip", source_id);
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
            text(&mut output, "roughness_unresolved", source_id);
        }
    }
    number(
        &mut output,
        "mean_free_path",
        input.applicability.mean_free_path_m,
    );
    number(
        &mut output,
        "mean_free_path_reference_pressure",
        input.initial_absolute_pressure_pa,
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
        GasFilmBoundaryTopology::MixedEnds {
            left_absolute_pressure_pa,
            right_absolute_pressure_pa,
        } => {
            text(&mut output, "boundary", "mixed-ends");
            for (name, pressure) in [
                ("left_reservoir", left_absolute_pressure_pa),
                ("right_reservoir", right_absolute_pressure_pa),
            ] {
                match pressure {
                    Some(pressure) => number(&mut output, name, pressure),
                    None => text(&mut output, name, "sealed"),
                }
            }
        }
        GasFilmBoundaryTopology::Vented {
            cell_index,
            absolute_pressure_pa,
        } => {
            text(&mut output, "boundary", "vented");
            count(&mut output, "vent_cell", cell_index);
            number(&mut output, "vent_reservoir", absolute_pressure_pa);
        }
        GasFilmBoundaryTopology::CapillaryVented {
            cell_index,
            absolute_pressure_pa,
            capillary,
        } => {
            text(&mut output, "boundary", "no-slip-capillary-vent-v1");
            count(&mut output, "vent_cell", cell_index);
            number(&mut output, "vent_reservoir", absolute_pressure_pa);
            number(&mut output, "capillary_radius", capillary.radius_m);
            number(&mut output, "capillary_length", capillary.length_m);
            number(
                &mut output,
                "capillaries_per_m",
                capillary.capillaries_per_m,
            );
        }
        GasFilmBoundaryTopology::CapillaryReservoir {
            cell_index,
            capillary,
            reservoir,
        } => {
            text(
                &mut output,
                "boundary",
                "finite-isothermal-capillary-reservoir-v1",
            );
            count(&mut output, "vent_cell", cell_index);
            number(
                &mut output,
                "reservoir_volume_per_width",
                reservoir.volume_per_width_m2,
            );
            number(
                &mut output,
                "reservoir_initial_pressure",
                reservoir.initial_absolute_pressure_pa,
            );
            number(&mut output, "capillary_radius", capillary.radius_m);
            number(&mut output, "capillary_length", capillary.length_m);
            number(
                &mut output,
                "capillaries_per_m",
                capillary.capillaries_per_m,
            );
        }
        GasFilmBoundaryTopology::AxisymmetricCapillaryVented {
            cell_index,
            absolute_pressure_pa,
            ..
        }
        | GasFilmBoundaryTopology::AxisymmetricCapillaryReservoir {
            cell_index,
            reservoir:
                AxisymmetricGasFilmReservoir {
                    initial_absolute_pressure_pa: absolute_pressure_pa,
                    ..
                },
            ..
        } => {
            // Includes total tube count/volume as well as all other invariants.
            // The bath pressure is dynamic; finite-reservoir initial pressure
            // is already invariant, but is also retained here unambiguously.
            text(
                &mut output,
                "configuration",
                &invariant_configuration_fingerprint(input),
            );
            count(&mut output, "vent_cell", cell_index);
            number(&mut output, "connection_pressure", absolute_pressure_pa);
        }
    }
    match &input.slip_policy {
        SlipPolicy::NoSlipContinuum { source_id } => text(&mut output, "slip", source_id),
        SlipPolicy::MaxwellFirstOrder {
            source_id,
            tangential_momentum_accommodation,
        } => {
            text(&mut output, "maxwell_first_order", source_id);
            number(
                &mut output,
                "maxwell_tmac",
                *tangential_momentum_accommodation,
            );
        }
        SlipPolicy::RarefiedSlipRequested { source_id } => {
            text(&mut output, "rarefied_slip", source_id);
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
            text(&mut output, "roughness_unresolved", source_id);
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
    match &input.wall_motion.gap_rate_profile_m_per_s {
        None => text(&mut output, "gap_rate_profile", "none"),
        Some(profile) => {
            text(&mut output, "gap_rate_profile", "per-cell");
            for value in profile {
                number(&mut output, "gap_rate_cell", *value);
            }
        }
    }
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
    metrics: &FilmMetrics,
    transport: &FilmTransport,
) -> Result<(f64, f64), GasFilmError> {
    let [left, right] = input.boundary.end_pressures();
    let active = pressure.len();
    let left_metric = metrics.face(0, active);
    let right_metric = metrics.face(active, active);
    Ok((
        if let Some(reservoir) = left {
            face_flux(
                reservoir,
                pressure[0],
                gap[0],
                gap[0],
                left_metric.distance_m,
                transport,
            )? * left_metric.width
        } else {
            0.0
        },
        if let Some(reservoir) = right {
            face_flux(
                pressure[active - 1],
                reservoir,
                gap[active - 1],
                gap[active - 1],
                right_metric.distance_m,
                transport,
            )? * right_metric.width
        } else {
            0.0
        },
    ))
}

fn face_fluxes(
    pressure: &[f64],
    gap: &[f64],
    input: &GasFilmInput,
    metrics: &FilmMetrics,
    transport: &FilmTransport,
    flux: &mut [f64],
) -> Result<(), GasFilmError> {
    let (left, right) = boundary_fluxes(pressure, gap, input, metrics, transport)?;
    flux[0] = derived(left, "left_face_mass_rate")?;
    flux[pressure.len()] = derived(right, "right_face_mass_rate")?;
    for face in 1..pressure.len() {
        let metric = metrics.face(face, pressure.len());
        flux[face] = derived(
            face_flux(
                pressure[face - 1],
                pressure[face],
                gap[face - 1],
                gap[face],
                metric.distance_m,
                transport,
            )? * metric.width,
            "interior_face_mass_rate",
        )?;
    }
    Ok(())
}

// One pressure system owns both the nonlinear residual and its exact Jacobian.
// Geometry and the compiled wall law are shared with the final receipts.
struct PressureSystem<'a> {
    input: &'a GasFilmInput,
    metrics: &'a FilmMetrics,
    transport: &'a FilmTransport,
    capillary: Option<&'a CapillaryConnection>,
    gap: &'a [f64],
    old_pressure: &'a [f64],
    old_gap: &'a [f64],
}

impl PressureSystem<'_> {
    fn residual(
        &self,
        pressure: &[f64],
        flux: &mut [f64],
        residual: &mut [f64],
    ) -> Result<(f64, Option<CapillaryExchange>), GasFilmError> {
        face_fluxes(
            pressure,
            self.gap,
            self.input,
            self.metrics,
            self.transport,
            flux,
        )?;
        let exchange = self
            .capillary
            .map(|vent| vent.exchange(pressure[vent.cell]))
            .transpose()?;
        let mut maximum = 0.0_f64;
        for cell in 0..pressure.len() {
            let measure = self.metrics.cell_measure(cell);
            let storage = derived(
                (pressure[cell] * self.gap[cell] - self.old_pressure[cell] * self.old_gap[cell])
                    / (self.transport.rt * self.input.timestep_s),
                "cell_storage_residual",
            )?;
            let mut value = derived(
                storage + (flux[cell + 1] - flux[cell]) / measure,
                "cell_mass_residual",
            )?;
            if self.capillary.is_some_and(|vent| cell == vent.cell) {
                value = derived(
                    value + exchange.expect("capillary was compiled").mass_flow / measure,
                    "vent_cell_mass_residual",
                )?;
            }
            residual[cell] = value;
            if !matches!(self.input.boundary, GasFilmBoundaryTopology::Vented { cell_index, .. } if cell == cell_index)
            {
                maximum = maximum.max(value.abs());
            }
        }
        if let Some(vent) = self.capillary {
            let value = vent.reservoir_residual(exchange.expect("capillary was compiled"))?;
            maximum = maximum.max(value.abs() / self.metrics.cell_measure(vent.cell));
        }
        Ok((maximum, exchange))
    }

    fn jacobian(
        &self,
        pressure: &[f64],
        exchange: Option<CapillaryExchange>,
        lower: &mut [f64],
        diagonal: &mut [f64],
        upper: &mut [f64],
    ) -> Result<(), GasFilmError> {
        let active = pressure.len();
        lower.fill(0.0);
        upper.fill(0.0);
        for (cell, value) in diagonal.iter_mut().enumerate() {
            *value = derived(
                self.gap[cell] / (self.transport.rt * self.input.timestep_s),
                "pressure_storage_derivative",
            )?;
        }
        let ends = self.input.boundary.end_pressures();
        for face in 0..=active {
            let (left_pressure, right_pressure, gap) = if face == 0 {
                let Some(bath) = ends[0] else { continue };
                (bath, pressure[0], self.gap[0])
            } else if face == active {
                let Some(bath) = ends[1] else { continue };
                (pressure[active - 1], bath, self.gap[active - 1])
            } else {
                (
                    pressure[face - 1],
                    pressure[face],
                    0.5 * (self.gap[face - 1] + self.gap[face]),
                )
            };
            let metric = self.metrics.face(face, active);
            // K(p)=a*p+c, including Maxwell b*p=constant. Thus the pressure
            // flux derivatives are +K(p_left)/d and -K(p_right)/d, while
            // Couette advection contributes h*U_mean/(2*R*T) to both.
            let advection =
                gap * self.transport.mean_wall_velocity_m_per_s / (2.0 * self.transport.rt);
            let left = derived(
                metric.width
                    * (self.transport.mass_conductance(left_pressure, gap)? / metric.distance_m
                        + advection),
                "left_pressure_flux_derivative",
            )?;
            let right = derived(
                metric.width
                    * (-self.transport.mass_conductance(right_pressure, gap)? / metric.distance_m
                        + advection),
                "right_pressure_flux_derivative",
            )?;
            if face > 0 {
                let measure = self.metrics.cell_measure(face - 1);
                diagonal[face - 1] = derived(
                    diagonal[face - 1] + left / measure,
                    "pressure_jacobian_diagonal",
                )?;
                if face < active {
                    upper[face - 1] = derived(right / measure, "pressure_jacobian_upper")?;
                }
            }
            if face < active {
                let measure = self.metrics.cell_measure(face);
                diagonal[face] = derived(
                    diagonal[face] - right / measure,
                    "pressure_jacobian_diagonal",
                )?;
                if face > 0 {
                    lower[face - 1] = derived(-left / measure, "pressure_jacobian_lower")?;
                }
            }
        }
        if let Some(vent) = self.capillary {
            diagonal[vent.cell] = derived(
                diagonal[vent.cell]
                    + exchange.expect("capillary was compiled").derivative
                        / self.metrics.cell_measure(vent.cell),
                "capillary_pressure_derivative",
            )?;
        }
        if let GasFilmBoundaryTopology::Vented { cell_index, .. } = self.input.boundary {
            // The pressure pin has zero Newton correction. Neighbor equations
            // retain their coupling to this row; no storage equation is invented.
            diagonal[cell_index] = 1.0;
            if cell_index > 0 {
                lower[cell_index - 1] = 0.0;
            }
            if cell_index + 1 < active {
                upper[cell_index] = 0.0;
            }
        }
        Ok(())
    }
}

// Integrate pressure-driven mechanics over the same face intervals as the
// Reynolds flux. Radial resistance uses a logarithmic dual measure, so pressure
// force work and dissipation close for the very same discrete mass flux.
// Square gradients BEFORE averaging so opposing/sharp gradients retain heat.
fn pressure_driven_mechanics(
    pressure: &[f64],
    gap: &[f64],
    flux: &[f64],
    input: &GasFilmInput,
    metrics: &FilmMetrics,
    transport: &FilmTransport,
) -> Result<(Vec<f64>, f64, f64, f64, f64), GasFilmError> {
    let mean_velocity = transport.mean_wall_velocity_m_per_s;
    let lower_velocity = input.wall_motion.lower_tangential_velocity_m_per_s;
    let upper_velocity = input.wall_motion.upper_tangential_velocity_m_per_s;
    let relative_velocity = derived(upper_velocity - lower_velocity, "relative_wall_velocity")?;
    let mut maximum_speed = 0.0_f64;
    let viscosity = input.gas.dynamic_viscosity_pa_s;
    let active = pressure.len();
    let mut face_shear = Vec::with_capacity(active + 1);
    let mut heat = 0.0;
    let mut pressure_work = 0.0;
    let mut wall_slip_heat = 0.0;
    let end_pressures = input.boundary.end_pressures();
    for (face, mass_flux) in flux.iter().copied().enumerate() {
        let metric = metrics.face(face, active);
        let distance = metric.distance_m;
        let (face_gap, face_pressure, gradient, distance) = if face > 0 && face < active {
            (
                0.5 * (gap[face - 1] + gap[face]),
                0.5 * (pressure[face - 1] + pressure[face]),
                (pressure[face] - pressure[face - 1]) / distance,
                distance,
            )
        } else {
            let cell = if face == 0 { 0 } else { active - 1 };
            match end_pressures[usize::from(face != 0)] {
                Some(reservoir) => {
                    let gradient = if face == 0 {
                        (pressure[cell] - reservoir) / distance
                    } else {
                        (reservoir - pressure[cell]) / distance
                    };
                    (
                        gap[cell],
                        0.5 * (pressure[cell] + reservoir),
                        gradient,
                        distance,
                    )
                }
                None => {
                    // A sealed face has q = h U_mean - h^3 p_x/(12 mu) = 0.
                    // Moving walls therefore require backpressure, not p_x = 0.
                    let gradient = 12.0 * viscosity * mean_velocity / gap[cell].powi(2);
                    (gap[cell], pressure[cell], gradient, distance)
                }
            }
        };
        let face_gap = derived(face_gap, "mechanical_face_gap_m")?;
        let face_pressure = derived(face_pressure, "mechanical_face_pressure_pa")?;
        let slip_length = transport.slip_length_m(face_pressure)?;
        let slip_ratio = derived(slip_length / face_gap, "wall_slip_length_ratio")?;
        let gradient =
            if (face == 0 || face == active) && end_pressures[usize::from(face != 0)].is_none() {
                // Zero flux uses the same slip-enhanced conductance as an open face.
                gradient / (1.0 + 6.0 * slip_ratio)
            } else {
                gradient
            };
        let gradient = derived(gradient, "face_pressure_gradient_pa_per_m")?;
        // Symmetric Maxwell slip gives u(s) = u_gas_lower + V*s + A*s*(s-1),
        // with V = (U_upper-U_lower)/(1+2*b/h) and s = y/h. The gas-side wall
        // values differ from the solid velocities when b > 0.
        // The endpoints and one stationary point give the exact profile peak.
        // A mean-flow or mid-gap check can miss an off-centre gas-speed maximum.
        let pressure_velocity = derived(
            face_gap.powi(2) * gradient / (2.0 * viscosity),
            "poiseuille_velocity_coefficient",
        )?;
        let (gas_lower_velocity, gas_upper_velocity, gas_relative_velocity) = if slip_length == 0.0
        {
            (lower_velocity, upper_velocity, relative_velocity)
        } else {
            let relative = derived(
                relative_velocity / (1.0 + 2.0 * slip_ratio),
                "couette_gas_velocity_difference",
            )?;
            let lower = derived(
                lower_velocity + slip_ratio * (relative - pressure_velocity),
                "lower_wall_gas_velocity",
            )?;
            let upper = derived(
                upper_velocity - slip_ratio * (relative + pressure_velocity),
                "upper_wall_gas_velocity",
            )?;
            (lower, upper, relative)
        };
        maximum_speed = maximum_speed
            .max(gas_lower_velocity.abs())
            .max(gas_upper_velocity.abs());
        if gas_relative_velocity.abs() < pressure_velocity.abs() {
            let stationary_point = 0.5 * (1.0 - gas_relative_velocity / pressure_velocity);
            let velocity = derived(
                gas_lower_velocity
                    + stationary_point
                        * (gas_relative_velocity + pressure_velocity * (stationary_point - 1.0)),
                "interior_gas_velocity",
            )?;
            maximum_speed = maximum_speed.max(velocity.abs());
        }
        let mach = derived(
            maximum_speed / input.applicability.speed_of_sound_m_per_s,
            "gas_mach_number",
        )?;
        if mach > input.applicability.maximum_mach_number {
            return Err(GasFilmError::Unavailable {
                reason: "gas-mach-outside-declared-envelope",
            });
        }
        let shear = derived(0.5 * face_gap * gradient, "pressure_wall_shear_pa")?;
        let conductance = derived(
            face_gap.powi(3) / (12.0 * viscosity),
            "reynolds_volume_conductance",
        )?;
        let slip_conductance = derived(
            face_gap.powi(2) * slip_length / (2.0 * viscosity),
            "maxwell_slip_volume_conductance",
        )?;
        heat = derived(
            heat + (conductance + slip_conductance) * gradient * gradient * distance * metric.width,
            "poiseuille_viscous_heat_w_per_m",
        )?;
        wall_slip_heat = derived(
            wall_slip_heat + slip_conductance * gradient * gradient * distance * metric.width,
            "pressure_wall_slip_dissipation_w_per_m",
        )?;
        let volume_flux = derived(
            mass_flux * transport.rt / face_pressure,
            "face_volume_flux_m2_per_s",
        )?;
        pressure_work = derived(
            pressure_work - volume_flux * gradient * distance,
            "pressure_gradient_power_to_gas_w_per_m",
        )?;
        face_shear.push(shear);
    }
    Ok((
        face_shear,
        heat,
        pressure_work,
        maximum_speed / input.applicability.speed_of_sound_m_per_s,
        wall_slip_heat,
    ))
}

/// Advance one bounded deterministic isothermal gas-film step.
///
/// A successful result is the only publication path.  Iteration-budget failure
/// returns an error and intentionally yields no partial checkpoint/receipt.
pub fn solve_isothermal_gas_film_1d(
    input: &GasFilmInput,
    checkpoint: Option<&GasFilmCheckpoint>,
) -> Result<GasFilmStep, GasFilmError> {
    solve_film(input, checkpoint, None)
}

// The internal integrals use the selected geometry measure. Only the Cartesian
// entry point publishes this per-width receipt directly; the radial entry point
// converts it into its distinct total-unit receipt before publication.
fn solve_film(
    input: &GasFilmInput,
    checkpoint: Option<&GasFilmCheckpoint>,
    radial_inner_radius_m: Option<f64>,
) -> Result<GasFilmStep, GasFilmError> {
    let (active, spacing, rt, max_slope) = input.validate()?;
    let metrics = match radial_inner_radius_m {
        Some(radius) => axisymmetric::compile_metrics(input, active, spacing, radius)?,
        None => {
            if input.boundary.is_axisymmetric_capillary() {
                return Err(GasFilmError::Unavailable {
                    reason: "cartesian-total-unit-capillary-connection",
                });
            }
            FilmMetrics::Cartesian { spacing_m: spacing }
        }
    };
    let configuration_fingerprint =
        metrics.bind_fingerprint(invariant_configuration_fingerprint(input));
    let transport = FilmTransport::new(input, rt)?;
    let (old_pressure, old_gap, old_step_index) =
        initial_state(input, active, checkpoint, &configuration_fingerprint)?;
    let capillary = CapillaryConnection::from_input(input, rt, checkpoint)?;
    let gap = &input.grid.gap_m[..active];
    // Both the starting state and the imposed reservoirs must belong to the
    // continuum model. A low-pressure boundary cannot be hidden by averaging
    // it with a denser interior cell, even for a very small timestep.
    for (pressure, gap) in old_pressure.iter().zip(&old_gap) {
        input.admit_knudsen(*pressure, *gap)?;
    }
    let mut max_knudsen = match input.boundary {
        GasFilmBoundaryTopology::Vented {
            cell_index,
            absolute_pressure_pa,
        } => input.admit_knudsen(absolute_pressure_pa, gap[cell_index])?,
        _ => 0.0,
    };
    let end_pressures = input.boundary.end_pressures();
    for (pressure, end_gap) in end_pressures.into_iter().zip([gap[0], gap[active - 1]]) {
        if let Some(pressure) = pressure {
            max_knudsen = max_knudsen.max(input.admit_knudsen(pressure, end_gap)?);
        }
    }
    if let Some(vent) = &capillary {
        // A capillary bath is outside the film: its Knudsen scale is the
        // tube radius, not the local film gap. Admit both tube endpoints and
        // the initial flow before the implicit pressure solve can hide them.
        vent.admit_initial(old_pressure[vent.cell], input)?;
    }
    let mut pressure = old_pressure.clone();
    if let GasFilmBoundaryTopology::Vented {
        cell_index,
        absolute_pressure_pa,
    } = input.boundary
    {
        pressure[cell_index] = absolute_pressure_pa;
    }
    let mass_old = old_pressure.iter().zip(&old_gap).enumerate().try_fold(
        0.0,
        |total, (cell, (pressure, gap))| {
            derived(
                total + pressure * gap * metrics.cell_measure(cell) / rt,
                "initial_gas_mass_kg_per_m",
            )
        },
    )?;
    let system = PressureSystem {
        input,
        metrics: &metrics,
        transport: &transport,
        capillary: capillary.as_ref(),
        gap,
        old_pressure: &old_pressure,
        old_gap: &old_gap,
    };
    // All nonlinear, line-search and pivoted linear-solve storage is retained
    // through this step. No cell/face vector is allocated inside the iteration.
    let mut flux = vec![0.0; active + 1];
    let mut cell_residual = vec![0.0; active];
    let mut next = pressure.clone();
    let mut correction = vec![0.0; active];
    let mut lower = vec![0.0; active - 1];
    let mut diagonal = vec![0.0; active];
    let mut upper = vec![0.0; active - 1];
    let mut linear = TridiagonalWorkspace::new(active)
        .map_err(|source| GasFilmError::LinearSolveFailure { source })?;
    let (mut max_residual, mut capillary_exchange) =
        system.residual(&pressure, &mut flux, &mut cell_residual)?;
    let mut iterations = 0;
    while max_residual > input.budget.mass_residual_tolerance_kg_m2_s
        && iterations < input.budget.maximum_iterations
    {
        system.jacobian(
            &pressure,
            capillary_exchange,
            &mut lower,
            &mut diagonal,
            &mut upper,
        )?;
        for (update, residual) in correction.iter_mut().zip(&cell_residual) {
            *update = -residual;
        }
        if let GasFilmBoundaryTopology::Vented { cell_index, .. } = input.boundary {
            correction[cell_index] = 0.0;
        }
        linear
            .factor(&lower, &diagonal, &upper, || false)
            .and_then(|()| linear.solve_factored_in_place(&mut correction, || false))
            .map_err(|source| GasFilmError::LinearSolveFailure { source })?;
        if let GasFilmBoundaryTopology::Vented { cell_index, .. } = input.boundary {
            correction[cell_index] = 0.0;
        }
        iterations += 1;
        let mut fraction = input.budget.relaxation;
        for (value, update) in pressure.iter().zip(&correction) {
            if *update < 0.0 {
                fraction = fraction.min(NEWTON_POSITIVE_FRACTION * (*value / -*update));
            }
        }
        let mut accepted = false;
        for _ in 0..MAX_NEWTON_BACKTRACKS {
            for ((trial, value), update) in next.iter_mut().zip(&pressure).zip(&correction) {
                *trial = value + fraction * update;
            }
            if fraction > 0.0 && next.iter().all(|value| value.is_finite() && *value > 0.0) {
                // A nonfinite trial can be shortened; it is never an accepted
                // state. Every trial uses the original mass/reservoir residual.
                if let Ok((trial_residual, trial_exchange)) =
                    system.residual(&next, &mut flux, &mut cell_residual)
                    && (trial_residual <= input.budget.mass_residual_tolerance_kg_m2_s
                        || (trial_residual < max_residual
                            && trial_residual
                                <= (1.0 - NEWTON_SUFFICIENT_DECREASE * fraction) * max_residual))
                {
                    std::mem::swap(&mut pressure, &mut next);
                    max_residual = trial_residual;
                    capillary_exchange = trial_exchange;
                    accepted = true;
                    break;
                }
            }
            fraction *= 0.5;
        }
        if !accepted {
            return Err(GasFilmError::NewtonStepRejected {
                iterations,
                max_mass_residual_kg_m2_s: max_residual,
            });
        }
    }
    if max_residual > input.budget.mass_residual_tolerance_kg_m2_s {
        return Err(GasFilmError::IterationBudgetExceeded {
            iterations,
            max_mass_residual_kg_m2_s: max_residual,
        });
    }
    for (pressure, gap) in pressure.iter().zip(gap) {
        max_knudsen = max_knudsen.max(input.admit_knudsen(*pressure, *gap)?);
    }
    let capillary_report = capillary
        .as_ref()
        .map(|vent| vent.observe(pressure[vent.cell], input))
        .transpose()?;
    let mass_new =
        pressure
            .iter()
            .zip(gap)
            .enumerate()
            .try_fold(0.0, |total, (cell, (pressure, gap))| {
                derived(
                    total + pressure * gap * metrics.cell_measure(cell) / rt,
                    "gas_mass",
                )
            })?;
    let storage_rate = derived(
        (mass_new - mass_old) / input.timestep_s,
        "storage_rate_kg_per_m_s",
    )?;
    let left_outward = derived(-flux[0], "left_boundary_outward_mass_flux")?;
    let right_outward = flux[active];
    let vent_outward = match input.boundary {
        GasFilmBoundaryTopology::Vented { cell_index, .. } => derived(
            -cell_residual[cell_index] * metrics.cell_measure(cell_index),
            "vent_outward_mass_flux",
        )?,
        GasFilmBoundaryTopology::CapillaryVented { .. }
        | GasFilmBoundaryTopology::CapillaryReservoir { .. }
        | GasFilmBoundaryTopology::AxisymmetricCapillaryVented { .. }
        | GasFilmBoundaryTopology::AxisymmetricCapillaryReservoir { .. } => {
            capillary_exchange
                .expect("capillary was compiled")
                .mass_flow
        }
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
    let reservoir_report = match (&capillary, capillary_exchange) {
        (Some(vent), Some(exchange)) => vent.reservoir_report(exchange, input, storage_rate)?,
        _ => None,
    };
    let mut absolute_load = 0.0;
    let mut gauge_load = 0.0;
    let mut gas_on_upper_normal = Vec::with_capacity(input.grid.gap_m.len());
    let mut upper_on_gas_normal = Vec::with_capacity(input.grid.gap_m.len());
    let mut gas_on_upper_tangential = Vec::with_capacity(input.grid.gap_m.len());
    let mut upper_on_gas_tangential = Vec::with_capacity(input.grid.gap_m.len());
    let mut gas_on_lower_tangential = Vec::with_capacity(input.grid.gap_m.len());
    let mut lower_on_gas_tangential = Vec::with_capacity(input.grid.gap_m.len());
    let mut upper_shear = Vec::with_capacity(input.grid.gap_m.len());
    let (
        pressure_shear,
        mut viscous_power,
        pressure_gradient_power,
        max_gas_mach,
        mut wall_slip_power,
    ) = pressure_driven_mechanics(&pressure, gap, &flux, input, &metrics, &transport)?;
    let mut tangential_wall_power = 0.0;
    let mut normal_gap_power = 0.0;
    let relative_velocity = input.wall_motion.upper_tangential_velocity_m_per_s
        - input.wall_motion.lower_tangential_velocity_m_per_s;
    for (index, gap) in input.grid.gap_m.iter().enumerate() {
        if index >= active {
            gas_on_upper_normal.push(None);
            upper_on_gas_normal.push(None);
            gas_on_upper_tangential.push(None);
            upper_on_gas_tangential.push(None);
            gas_on_lower_tangential.push(None);
            lower_on_gas_tangential.push(None);
            upper_shear.push(None);
            continue;
        }
        let measure = metrics.cell_measure(index);
        let pressure_value = pressure[index];
        absolute_load = derived(
            absolute_load + pressure_value * measure,
            "absolute_pressure_load_n_per_m",
        )?;
        gauge_load = derived(
            gauge_load + (pressure_value - input.gauge_reference_absolute_pressure_pa) * measure,
            "gauge_pressure_load_n_per_m",
        )?;
        normal_gap_power = derived(
            normal_gap_power - pressure_value * input.gap_rate_m_per_s(index) * measure,
            "normal_gap_power_to_gas_w_per_m",
        )?;
        let slip_length = transport.slip_length_m(pressure_value)?;
        let couette_gap = derived(gap + 2.0 * slip_length, "couette_effective_gap_m")?;
        let couette_shear = derived(
            input.gas.dynamic_viscosity_pa_s * relative_velocity / couette_gap,
            "couette_wall_shear_pa",
        )?;
        // A local face-average traction. Radial work/heat is integrated directly
        // on the dual faces; axisymmetric tangential wall motion is not admitted.
        let pressure_shear = 0.5 * (pressure_shear[index] + pressure_shear[index + 1]);
        let upper_shear_value = derived(couette_shear + pressure_shear, "upper_wall_shear_pa")?;
        let lower_shear_value = derived(-couette_shear + pressure_shear, "lower_wall_shear_pa")?;
        gas_on_upper_normal.push(Some(pressure_value));
        upper_on_gas_normal.push(Some(-pressure_value));
        gas_on_upper_tangential.push(Some(-upper_shear_value));
        upper_on_gas_tangential.push(Some(upper_shear_value));
        gas_on_lower_tangential.push(Some(-lower_shear_value));
        lower_on_gas_tangential.push(Some(lower_shear_value));
        upper_shear.push(Some(upper_shear_value));
        viscous_power = derived(
            viscous_power + couette_shear * relative_velocity * measure,
            "viscous_heat_w_per_m",
        )?;
        wall_slip_power = derived(
            wall_slip_power
                + (2.0 * slip_length / couette_gap) * couette_shear * relative_velocity * measure,
            "couette_wall_slip_dissipation_w_per_m",
        )?;
        tangential_wall_power = derived(
            tangential_wall_power
                + (upper_shear_value * input.wall_motion.upper_tangential_velocity_m_per_s
                    + lower_shear_value * input.wall_motion.lower_tangential_velocity_m_per_s)
                    * measure,
            "tangential_wall_power_to_gas_w_per_m",
        )?;
    }
    if viscous_power < -1.0e-18 {
        return Err(GasFilmError::NonFiniteDerived {
            field: "negative_viscous_dissipation",
        });
    }
    let wall_power = derived(
        tangential_wall_power + normal_gap_power,
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
            capillary_vent: capillary_report,
            finite_reservoir: reservoir_report,
            mass_closure_residual_kg_per_m_s: closure,
            absolute_pressure_load_n_per_m: absolute_load,
            gauge_pressure_load_n_per_m: gauge_load,
            gas_on_upper_wall_normal_traction_pa: gas_on_upper_normal,
            upper_wall_on_gas_normal_traction_pa: upper_on_gas_normal,
            gas_on_upper_wall_tangential_traction_pa: gas_on_upper_tangential,
            upper_wall_on_gas_tangential_traction_pa: upper_on_gas_tangential,
            gas_on_lower_wall_tangential_traction_pa: gas_on_lower_tangential,
            lower_wall_on_gas_tangential_traction_pa: lower_on_gas_tangential,
            upper_wall_shear_pa: upper_shear,
            normal_gap_power_to_gas_w_per_m: normal_gap_power,
            wall_power_to_gas_w_per_m: wall_power,
            pressure_gradient_power_to_gas_w_per_m: pressure_gradient_power,
            viscous_heat_w_per_m: viscous_power,
            wall_slip_dissipation_w_per_m: wall_slip_power,
            maximum_knudsen_number: max_knudsen,
            maximum_gas_mach_number: max_gas_mach,
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
            configuration_fingerprint,
            accepted_request_fingerprint: metrics
                .bind_fingerprint(accepted_request_fingerprint(input)),
            active_absolute_pressure_pa: checkpoint_pressure,
            active_gap_m: gap.to_vec(),
            active_cells: active,
            capillary_reservoir_mass_kg_per_m: capillary_exchange
                .filter(|_| radial_inner_radius_m.is_none())
                .and_then(|exchange| exchange.reservoir_mass),
            capillary_reservoir_mass_kg: capillary_exchange
                .filter(|_| radial_inner_radius_m.is_some())
                .and_then(|exchange| exchange.reservoir_mass),
        },
    })
}

/// Canonical model identity for callers constructing a [`GasFilmIdentity`].
#[must_use]
pub const fn isothermal_compressible_reynolds_model_id() -> &'static str {
    MODEL_ID
}
