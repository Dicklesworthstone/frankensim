//! L3 (`FLUX`) stirred multiphase liquid/vapor/oil/emulsion vessel
//! [`fs_ports::component::Component`] (`MultiphaseVesselComponent`).
//!
//! Unifies stockpots, saucepans, sauté skillets, pressure cookers, and deep-fryer
//! vats behind one conservative `Component` contract on the shared
//! `janaf-298.15K-formation` enthalpy datum with zero raw `libm` calls (100%
//! [`fs_math::det`]).
//!
//! # Port Layout (All Storage Ports Use `PortRole::FlowIn`, Positive INTO Vessel)
//! 0. `wall_heat`: scalar `Carrier::Heat` (receives `Q_wall` [W], returns wall-contact fluid temp [K])
//! 1. `immersed_food_heat`: `PortSpec::heat_field` width `M >= 1` (receives per-batch `Q_food` [W], returns bath temp [K])
//! 2. `liquid_stream`: `PortSpec::identified_stream` on `kitchen_liquid_v1` (`janaf-298.15K-formation`)
//! 3. `vapor_stream`: `PortSpec::identified_stream` on `kitchen_gas_v1` (`janaf-298.15K-formation`)
//! 4. `stir_rpm`: width-1 `SignalIn` (impeller/whisk/spoon speed [RPM])
//! 5. `p_vessel_pa`: width-1 `SignalIn` (explicit headspace pressure `P_atm(z)` or regulated gauge pressure [Pa])
//! 6. `probes`: width-12 `SignalOut`:
//!    `[T_liq_K, T_coldzone_K, m_liq_kg, x_etoh, x_solids, visc_pa_s, ffa_pct, tpm_pct, T_smoke_K, foam_h_m, fond_adhered_kg, safety_or_boilover_flag]`.

use core::fmt;
use core::num::NonZeroUsize;

use fs_exec::solver::snapshot_v2::StateIdentityCharterV2;
use fs_exec::{Cancelled, Cx};
use fs_matdb::ConstitutiveModelCard;
use fs_math::det::{exp, expm1, ln, pow, sqrt, tanh};
use fs_ports::StableId;
use fs_ports::component::{
    Class, Component, ContractViolation, ElementalMatrix, MolarMass, PortLayout, PortRole,
    PortSpec, Proposal, Receipt, STREAM_CARRIED, STREAM_ENTHALPY_RATE, STREAM_MASS_RATE,
    STREAM_PRESSURE, STREAM_SPECIFIC_ENTHALPY, SpeciesId, StreamQuantity, kitchen_gas_v1,
    kitchen_liquid_v1,
};
use fs_qty::chemistry::ElementId;
use fs_qty::{Dims, Pressure};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;
use fs_thermochem::caloric::proximate::{
    H2O_LIQUID_FORMATION_ENTHALPY_298K_J_KG, JANAF_REFERENCE_TEMPERATURE_K,
};
use fs_thermochem::solution::VleSolverOptionsV1;
use fs_thermochem::solution::nrtl::NrtlModelV1;
use fs_thermochem::water::if97::If97SaturationPair;

/// Fixed width of the `probes` `SignalOut` port on [`MultiphaseVesselComponent`].
pub const VESSEL_PROBE_COUNT: usize = 12;
/// Width of the `liquid_stream` port (`2` header + `8` `kitchen_liquid_v1` carried species).
pub const LIQUID_STREAM_WIDTH: usize = 10;
/// Width of the `vapor_stream` port (`2` header + `9` `kitchen_gas_v1` carried species).
pub const VAPOR_STREAM_WIDTH: usize = 11;

/// Probe index: bulk active liquid temperature [K].
pub const PROBE_T_LIQ_K: usize = 0;
/// Probe index: sub-tube cold-zone temperature [K] (equals `T_liq_K` when disabled).
pub const PROBE_T_COLDZONE_K: usize = 1;
/// Probe index: free liquid inventory mass [kg] (active + cold-zone liquid, excluding adhered fond and vapor).
pub const PROBE_M_LIQ_KG: usize = 2;
/// Probe index: ethanol mole fraction `x_EtOH` in volatile aqueous-ethanol phase.
pub const PROBE_X_ETOH: usize = 3;
/// Probe index: non-volatile dissolved + suspended solids mass fraction in active liquid.
pub const PROBE_X_SOLIDS: usize = 4;
/// Probe index: effective liquid dynamic viscosity [Pa s].
pub const PROBE_VISC_PA_S: usize = 5;
/// Probe index: free fatty acid percentage `100 * m_ffa / m_lipid_total` [%].
pub const PROBE_FFA_PCT: usize = 6;
/// Probe index: total polar materials percentage `100 * tpm_fraction` [%].
pub const PROBE_TPM_PCT: usize = 7;
/// Probe index: degraded cooking-oil smoke point [K] (`0.0` for non-lipid vessels).
pub const PROBE_T_SMOKE_K: usize = 8;
/// Probe index: Bikerman foam column height [m].
pub const PROBE_FOAM_H_M: usize = 9;
/// Probe index: undissolved fond mass adhered to pan floor [kg].
pub const PROBE_FOND_ADHERED_KG: usize = 10;
/// Probe index: safety / boil-over bitmask (`1.0` boil-over, `2.0` smoke, `4.0` flash, `8.0` autoignition, `16.0` oil discard, `32.0` scorch risk).
pub const PROBE_SAFETY_OR_BOILOVER_FLAG: usize = 11;

/// Exact `fs-matdb` law id for kitchen cooking-oil and butter cards.
pub const COOKING_OIL_LAW_ID_V1: &str = "kitchen-cooking-oil";
/// Supported semantic law version for [`COOKING_OIL_LAW_ID_V1`].
pub const COOKING_OIL_LAW_VERSION_V1: u32 = 1;

// Durable snapshot codec tag and scalar field count.
const VESSEL_STATE_CODEC_TAG: u64 = 0x4653_4D50_5653_0001;
const VESSEL_STATE_SCALAR_COUNT: usize = 31;
const AUX_GAS_COUNT: usize = 7;

// Exact IUPAC atomic weights [kg/mol] guaranteeing
// M_TRIGLYCERIDE + 3 * M_H2O == 3 * M_FFA + M_GLYCEROL to machine precision.
const ATOMIC_MASS_C_KG_MOL: f64 = 0.012_011;
const ATOMIC_MASS_H_KG_MOL: f64 = 0.001_007_94;
const ATOMIC_MASS_O_KG_MOL: f64 = 0.015_999_4;

/// Exact molar mass of water `H2O` [kg/mol] used in hydrolysis stoichiometry.
pub const MOLAR_MASS_H2O_KG_MOL: f64 = 2.0 * ATOMIC_MASS_H_KG_MOL + ATOMIC_MASS_O_KG_MOL;
/// Exact molar mass of triolein (`C57H104O6`, representative frying triglyceride) [kg/mol].
pub const MOLAR_MASS_TRIGLYCERIDE_KG_MOL: f64 =
    57.0 * ATOMIC_MASS_C_KG_MOL + 104.0 * ATOMIC_MASS_H_KG_MOL + 6.0 * ATOMIC_MASS_O_KG_MOL;
/// Exact molar mass of oleic free fatty acid (`C18H34O2`) [kg/mol].
pub const MOLAR_MASS_FFA_KG_MOL: f64 =
    18.0 * ATOMIC_MASS_C_KG_MOL + 34.0 * ATOMIC_MASS_H_KG_MOL + 2.0 * ATOMIC_MASS_O_KG_MOL;
/// Exact molar mass of glycerol (`C3H8O3`) [kg/mol].
pub const MOLAR_MASS_GLYCEROL_KG_MOL: f64 =
    3.0 * ATOMIC_MASS_C_KG_MOL + 8.0 * ATOMIC_MASS_H_KG_MOL + 3.0 * ATOMIC_MASS_O_KG_MOL;
/// Molar mass of ethanol `C2H5OH` [kg/mol].
pub const MOLAR_MASS_ETOH_KG_MOL: f64 = 0.046_068_44;
/// Molar mass of sodium chloride `NaCl` [kg/mol].
pub const MOLAR_MASS_NACL_KG_MOL: f64 = 0.058_442_8;
/// Molar mass of sucrose `C12H22O11` [kg/mol].
pub const MOLAR_MASS_SUCROSE_KG_MOL: f64 = 0.342_296_5;
// Effective molar mass for dissolved Maillard/peptide fond solutes [kg/mol].
const MOLAR_MASS_FOND_EFF_KG_MOL: f64 = 0.250;

// Unified janaf-298.15K-formation specific enthalpies [J/kg] and heat capacities [J/(kg K)].
const H_F0_H2O_LIQ_J_KG: f64 = H2O_LIQUID_FORMATION_ENTHALPY_298K_J_KG;
// CODATA/NIST-JANAF H2O(g) formation enthalpy (-241,826 J/mol / 0.018015268 kg/mol).
const H_F0_H2O_VAP_J_KG: f64 = -241_826.0 / 0.018_015_268;
const CP_H2O_LIQ_J_KG_K: f64 = 4_200.0;
const CP_H2O_VAP_J_KG_K: f64 = 1_875.0;

// NIST-JANAF ethanol liquid (-277.6 kJ/mol) and vapor (-235.3 kJ/mol).
const H_F0_ETOH_LIQ_J_KG: f64 = -277_600.0 / MOLAR_MASS_ETOH_KG_MOL;
const H_F0_ETOH_VAP_J_KG: f64 = -235_300.0 / MOLAR_MASS_ETOH_KG_MOL;
const CP_ETOH_LIQ_J_KG_K: f64 = 2_440.0;
const CP_ETOH_VAP_J_KG_K: f64 = 1_420.0;

// NIST-JANAF liquid lipid & hydrolysis products (endothermic Delta H_rxn ~ +78 kJ/mol).
const H_F0_TRIGLYCERIDE_J_KG: f64 = -2_450_000.0;
const H_F0_FFA_J_KG: f64 = -2_690_000.0;
const H_F0_GLYCEROL_J_KG: f64 = -7_260_000.0;
const CP_TRIG_DEFAULT_J_KG_K: f64 = 1_985.0;
const DCP_DT_TRIG_DEFAULT_J_KG_K2: f64 = 2.25;
const CP_FFA_J_KG_K: f64 = 2_100.0;
const CP_GLYCEROL_J_KG_K: f64 = 2_400.0;

// Solutes & internal non-volatile solids.
const H_F0_SALT_J_KG: f64 = -6_967_000.0;
const CP_SALT_J_KG_K: f64 = 880.0;
const H_F0_SUGAR_J_KG: f64 = -6_491_000.0;
const CP_SUGAR_J_KG_K: f64 = 1_244.0;
const H_F0_FOND_J_KG: f64 = -5_000_000.0;
const CP_FOND_J_KG_K: f64 = 1_650.0;
const H_F0_MILK_SOLIDS_J_KG: f64 = -5_200_000.0;
const CP_MILK_SOLIDS_J_KG_K: f64 = 1_750.0;
const H_F0_STARCH_J_KG: f64 = -5_800_000.0;
const CP_STARCH_J_KG_K: f64 = 1_550.0;
const H_F0_PROTEIN_J_KG: f64 = -4_800_000.0;
const CP_PROTEIN_J_KG_K: f64 = 2_000.0;

// Auxiliary headspace gases in kitchen_gas_v1 order after [H2O, C2H5OH]:
// [CO2, O2, N2, CH4, C3H8, CO, grease_vapor].
const AUX_GAS_H_F0_J_KG: [f64; AUX_GAS_COUNT] = [
    -8_941_490.0, // CO2
    0.0,          // O2
    0.0,          // N2
    -4_662_000.0, // CH4
    -2_354_000.0, // C3H8
    -3_946_000.0, // CO
    -2_100_000.0, // grease_vapor
];
const AUX_GAS_CP_J_KG_K: [f64; AUX_GAS_COUNT] = [
    846.0,   // CO2
    918.0,   // O2
    1_040.0, // N2
    2_226.0, // CH4
    1_670.0, // C3H8
    1_040.0, // CO
    1_850.0, // grease_vapor
];
const AUX_GAS_MOLAR_MASS_KG_MOL: [f64; AUX_GAS_COUNT] = [
    0.044_009_5, // CO2
    0.031_998_8, // O2
    0.028_013_4, // N2
    0.016_042_5, // CH4
    0.044_095_6, // C3H8
    0.028_010_1, // CO
    0.282_466_8, // grease_vapor
];

const TEMP_MIN_VALID_K: f64 = 200.0;
const TEMP_MAX_VALID_K: f64 = 1_500.0;

// SI dimensions for `kitchen-cooking-oil` matdb card validation.
const DIMLESS: Dims = Dims([0, 0, 0, 0, 0, 0]);
const TEMP_DIMS: Dims = Dims([0, 0, 0, 1, 0, 0]);
const DENSITY_DIMS: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const DRHO_DT_DIMS: Dims = Dims([-3, 1, 0, -1, 0, 0]);
const CP_DIMS: Dims = Dims([2, 0, -2, -1, 0, 0]);
const DCP_DT_DIMS: Dims = Dims([2, 0, -2, -2, 0, 0]);
const K_DIMS: Dims = Dims([1, 1, -3, -1, 0, 0]);
const DK_DT_DIMS: Dims = Dims([1, 1, -3, -2, 0, 0]);
const VISC_DIMS: Dims = Dims([-1, 1, -1, 0, 0, 0]);
const MOLAR_ENERGY_DIMS: Dims = Dims([2, 1, -2, 0, 0, -1]);
const SURF_TENS_DIMS: Dims = Dims([0, 1, -2, 0, 0, 0]);
const DSIGMA_DT_DIMS: Dims = Dims([0, 1, -2, -1, 0, 0]);
const RATE_DIMS: Dims = Dims([0, 0, -1, 0, 0, 0]);

/// Typed refusal enum for [`MultiphaseVesselComponent`] configuration, state,
/// card loading, and step execution.
#[derive(Clone, Debug, PartialEq)]
pub enum VesselErrorV1 {
    /// Invalid vessel geometry, port, or physical parameter.
    InvalidConfig {
        /// Parameter or field name that failed validation.
        parameter: &'static str,
        /// Stable human-readable explanation.
        reason: &'static str,
    },
    /// Invalid scalar runtime input or state value.
    InvalidScalar {
        /// Offending physical quantity name.
        quantity: &'static str,
        /// Rejected value.
        value: f64,
        /// Stable explanation.
        reason: &'static str,
    },
    /// Vessel headspace pressure `p_vessel_pa` is non-finite, non-positive, or out of range.
    InvalidPressure {
        /// Rejected pressure [Pa].
        pressure_pa: f64,
    },
    /// `fs-matdb` cooking-oil or butter card failed validation.
    InvalidCard {
        /// Card or parameter context.
        context: &'static str,
        /// Stable explanation.
        reason: &'static str,
    },
    /// Stream port total mass flow disagreed with the sum of its carried species rates.
    InconsistentStreamMassFlow {
        /// Stream port identifier (`"liquid_stream"` or `"vapor_stream"`).
        port: &'static str,
        /// Declared `STREAM_MASS_RATE` [kg/s].
        declared_total_kg_s: f64,
        /// Sum of carried species mass flow rates [kg/s].
        species_sum_kg_s: f64,
    },
    /// Outgoing stream flux attempted to withdraw mass from a depleted species inventory.
    DepletedInventory {
        /// Species identifier that would become negative.
        species: &'static str,
    },
    /// Nonzero thermal energy entered a completely empty zero-mass vessel with zero dry heat capacity.
    EmptyVesselThermalSingularity,
    /// Underlying `fs-thermochem` IAPWS-IF97 or NRTL VLE evaluation failed.
    Thermochem(&'static str),
    /// `fs-ports` layout or basis contract violation.
    PortContract(ContractViolation),
    /// Temperature inversion or sub-solver failed to produce a valid physical temperature.
    SolverConvergence {
        /// Stage label.
        stage: &'static str,
        /// Residual magnitude.
        residual: f64,
        /// Target tolerance.
        tolerance: f64,
    },
    /// Cooperative cancellation was triggered via [`fs_exec::Cx`].
    Cancelled(Cancelled),
}

impl fmt::Display for VesselErrorV1 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidConfig { parameter, reason } => {
                write!(f, "invalid vessel config `{parameter}`: {reason}")
            }
            Self::InvalidScalar {
                quantity,
                value,
                reason,
            } => write!(f, "invalid vessel scalar `{quantity}` = {value}: {reason}"),
            Self::InvalidPressure { pressure_pa } => {
                write!(f, "invalid vessel pressure {pressure_pa} Pa")
            }
            Self::InvalidCard { context, reason } => {
                write!(f, "invalid vessel oil/butter card `{context}`: {reason}")
            }
            Self::InconsistentStreamMassFlow {
                port,
                declared_total_kg_s,
                species_sum_kg_s,
            } => write!(
                f,
                "inconsistent stream mass flow on `{port}`: total={declared_total_kg_s} kg/s vs species sum={species_sum_kg_s} kg/s"
            ),
            Self::DepletedInventory { species } => {
                write!(f, "depleted vessel inventory for species `{species}`")
            }
            Self::EmptyVesselThermalSingularity => write!(
                f,
                "nonzero heat/enthalpy entered a completely empty vessel with zero dry heat capacity"
            ),
            Self::Thermochem(msg) => write!(f, "vessel thermochem lookup failed: {msg}"),
            Self::PortContract(err) => write!(f, "vessel port contract error: {err}"),
            Self::SolverConvergence {
                stage,
                residual,
                tolerance,
            } => write!(
                f,
                "vessel solver stage `{stage}` failed: residual={residual} > tol={tolerance}"
            ),
            Self::Cancelled(c) => write!(f, "vessel step cancelled: {c:?}"),
        }
    }
}

impl core::error::Error for VesselErrorV1 {}

impl From<Cancelled> for VesselErrorV1 {
    fn from(value: Cancelled) -> Self {
        Self::Cancelled(value)
    }
}

impl From<ContractViolation> for VesselErrorV1 {
    fn from(value: ContractViolation) -> Self {
        Self::PortContract(value)
    }
}

/// Thermophysical, thermal-safety, emulsion, and aging-kinetics properties for a
/// cooking oil, animal fat, ghee, or butter medium (`kitchen-cooking-oil` v1).
#[derive(Clone, Debug, PartialEq)]
pub struct OilCardV1 {
    /// Triglyceride lipid mass fraction in fresh medium (`in (0, 1]`).
    pub lipid_mass_fraction: f64,
    /// Water mass fraction in fresh medium (`in [0, 0.5]`, ~0.175 for unsalted butter).
    pub water_mass_fraction: f64,
    /// Non-fat milk solids mass fraction (`in [0, 0.2]`, ~0.014 for butter).
    pub milk_solids_mass_fraction: f64,
    /// Salt (`NaCl`) mass fraction (`in [0, 0.1]`, ~0.016 for salted butter).
    pub salt_nacl_mass_fraction: f64,
    /// Reference temperature `T_ref` [K] (`298.15 K`).
    pub t_ref_k: f64,
    /// Liquid density at `T_ref` [kg/m^3].
    pub rho_ref_kg_m3: f64,
    /// Density temperature slope `drho/dT` [kg/(m^3 K)] (typically negative).
    pub drho_dt_kg_m3_k: f64,
    /// Specific heat at `T_ref` [J/(kg K)].
    pub cp_ref_j_kg_k: f64,
    /// Specific heat temperature slope `dc_p/dT` [J/(kg K^2)] (`>= 0`).
    pub dcp_dt_j_kg_k2: f64,
    /// Thermal conductivity at `T_ref` [W/(m K)].
    pub k_ref_w_m_k: f64,
    /// Thermal conductivity temperature slope `dk/dT` [W/(m K^2)].
    pub dk_dt_w_m_k2: f64,
    /// Arrhenius viscosity pre-exponential `mu_0` [Pa s].
    pub mu_0_pa_s: f64,
    /// Arrhenius viscosity activation energy `E_visc` [J/mol].
    pub e_visc_j_mol: f64,
    /// Surface tension at `T_ref` [N/m].
    pub sigma_ref_n_m: f64,
    /// Surface tension temperature slope `dsigma/dT` [N/(m K)].
    pub dsigma_dt_n_m_k: f64,
    /// Fresh smoke point `T_smoke,0` [K].
    pub t_smoke_fresh_k: f64,
    /// Open-cup flash point `T_flash` [K] (~588 K / 315 C for canola).
    pub t_flash_k: f64,
    /// Autoignition temperature `T_autoignition` [K] (~651 K / 378 C for canola).
    pub t_autoignition_k: f64,
    /// Hydrolysis Arrhenius pre-exponential `k_hyd` [1/s].
    pub k_hyd_s_inv: f64,
    /// Hydrolysis activation energy `Ea_hyd` [J/mol].
    pub ea_hyd_j_mol: f64,
    /// Thermo-oxidative TPM Arrhenius pre-exponential `k_ox` [1/s].
    pub k_ox_s_inv: f64,
    /// Oxidation activation energy `Ea_ox` [J/mol].
    pub ea_ox_j_mol: f64,
    /// Thermal polymerization Arrhenius pre-exponential `k_poly` [1/s].
    pub k_poly_s_inv: f64,
    /// Polymerization activation energy `Ea_poly` [J/mol].
    pub ea_poly_j_mol: f64,
    /// Positive smoke-point depression per unit FFA mass fraction `beta_FFA` [K].
    pub beta_smoke_ffa_k: f64,
    /// Positive smoke-point depression per unit TPM mass fraction `beta_TPM` [K].
    pub beta_smoke_tpm_k: f64,
    /// Viscosity TPM exponent `gamma_TPM` in `mu(T, TPM) = mu_0(T) * exp(gamma_TPM * TPM)`.
    pub gamma_tpm: f64,
    /// Fresh oil FFA mass fraction (`0.0005` = 0.05% for refined canola).
    pub ffa_fresh: f64,
    /// Commercial discard threshold for TPM mass fraction (`0.25` = 25%).
    pub tpm_discard_limit: f64,
    /// Commercial discard threshold for FFA mass fraction (`0.02` = 2.0%).
    pub ffa_discard_limit: f64,
    /// Effective steam-oil bubble contact residence time `tau_hyd` [s] coupling
    /// `k_hyd(T) * tau_hyd * (m_dot_steam_food / m_oil)` into `[FFA]` formation.
    pub hydrolysis_contact_time_s: f64,
    /// Effective oil-air surface oxygenation mass-area factor [kg/m^2] coupling
    /// `k_ox(T) * factor * (A_air / m_oil)` into `[TPM]` formation.
    pub oxidation_surface_factor_kg_m2: f64,
}

impl OilCardV1 {
    /// Refined canola oil matching `data/matdb/kitchen-v1/oils/oil-canola-refined`.
    #[must_use]
    pub const fn canola_refined() -> Self {
        Self {
            lipid_mass_fraction: 0.999,
            water_mass_fraction: 0.001,
            milk_solids_mass_fraction: 0.0,
            salt_nacl_mass_fraction: 0.0,
            t_ref_k: JANAF_REFERENCE_TEMPERATURE_K,
            rho_ref_kg_m3: 914.0,
            drho_dt_kg_m3_k: -0.66,
            cp_ref_j_kg_k: 1_985.0,
            dcp_dt_j_kg_k2: 2.25,
            k_ref_w_m_k: 0.168,
            dk_dt_w_m_k2: -0.000_21,
            mu_0_pa_s: 2.4e-6,
            e_visc_j_mol: 24_900.0,
            sigma_ref_n_m: 0.032,
            dsigma_dt_n_m_k: -6.7e-5,
            t_smoke_fresh_k: 477.15,
            t_flash_k: 588.15,
            t_autoignition_k: 651.15,
            k_hyd_s_inv: 110_000.0,
            ea_hyd_j_mol: 84_000.0,
            k_ox_s_inv: 360_000.0,
            ea_ox_j_mol: 88_000.0,
            k_poly_s_inv: 200_000.0,
            ea_poly_j_mol: 94_000.0,
            beta_smoke_ffa_k: 1_250.0,
            beta_smoke_tpm_k: 120.0,
            gamma_tpm: 2.6,
            ffa_fresh: 0.0005,
            tpm_discard_limit: 0.25,
            ffa_discard_limit: 0.02,
            hydrolysis_contact_time_s: 50.0,
            oxidation_surface_factor_kg_m2: 7.5,
        }
    }

    /// Refined peanut oil matching `data/matdb/kitchen-v1/oils/oil-peanut-refined`.
    #[must_use]
    pub const fn peanut_refined() -> Self {
        Self {
            lipid_mass_fraction: 0.999,
            water_mass_fraction: 0.001,
            milk_solids_mass_fraction: 0.0,
            salt_nacl_mass_fraction: 0.0,
            t_ref_k: JANAF_REFERENCE_TEMPERATURE_K,
            rho_ref_kg_m3: 912.0,
            drho_dt_kg_m3_k: -0.65,
            cp_ref_j_kg_k: 1_975.0,
            dcp_dt_j_kg_k2: 2.20,
            k_ref_w_m_k: 0.169,
            dk_dt_w_m_k2: -0.000_21,
            mu_0_pa_s: 2.2e-6,
            e_visc_j_mol: 25_400.0,
            sigma_ref_n_m: 0.0322,
            dsigma_dt_n_m_k: -6.6e-5,
            t_smoke_fresh_k: 505.15,
            t_flash_k: 598.15,
            t_autoignition_k: 665.15,
            k_hyd_s_inv: 95_000.0,
            ea_hyd_j_mol: 85_000.0,
            k_ox_s_inv: 280_000.0,
            ea_ox_j_mol: 89_500.0,
            k_poly_s_inv: 170_000.0,
            ea_poly_j_mol: 95_000.0,
            beta_smoke_ffa_k: 1_250.0,
            beta_smoke_tpm_k: 120.0,
            gamma_tpm: 2.5,
            ffa_fresh: 0.0005,
            tpm_discard_limit: 0.25,
            ffa_discard_limit: 0.02,
            hydrolysis_contact_time_s: 50.0,
            oxidation_surface_factor_kg_m2: 7.5,
        }
    }

    /// Extra-virgin olive oil matching `data/matdb/kitchen-v1/oils/oil-olive-extra-virgin`.
    #[must_use]
    pub const fn extra_virgin_olive() -> Self {
        Self {
            lipid_mass_fraction: 0.998,
            water_mass_fraction: 0.002,
            milk_solids_mass_fraction: 0.0,
            salt_nacl_mass_fraction: 0.0,
            t_ref_k: JANAF_REFERENCE_TEMPERATURE_K,
            rho_ref_kg_m3: 911.0,
            drho_dt_kg_m3_k: -0.65,
            cp_ref_j_kg_k: 1_970.0,
            dcp_dt_j_kg_k2: 2.15,
            k_ref_w_m_k: 0.169,
            dk_dt_w_m_k2: -0.000_20,
            mu_0_pa_s: 2.0e-6,
            e_visc_j_mol: 26_100.0,
            sigma_ref_n_m: 0.0324,
            dsigma_dt_n_m_k: -6.6e-5,
            t_smoke_fresh_k: 433.15,
            t_flash_k: 573.15,
            t_autoignition_k: 635.15,
            k_hyd_s_inv: 140_000.0,
            ea_hyd_j_mol: 83_000.0,
            k_ox_s_inv: 260_000.0,
            ea_ox_j_mol: 89_000.0,
            k_poly_s_inv: 150_000.0,
            ea_poly_j_mol: 95_000.0,
            beta_smoke_ffa_k: 1_150.0,
            beta_smoke_tpm_k: 110.0,
            gamma_tpm: 2.4,
            ffa_fresh: 0.005,
            tpm_discard_limit: 0.25,
            ffa_discard_limit: 0.02,
            hydrolysis_contact_time_s: 50.0,
            oxidation_surface_factor_kg_m2: 7.5,
        }
    }

    /// Unsalted sweet-cream butter matching `data/matdb/kitchen-v1/oils/butter-unsalted`
    /// (`81.1%` butterfat, `17.5%` water emulsion, `1.4%` non-fat milk solids).
    #[must_use]
    pub const fn butter_unsalted() -> Self {
        Self {
            lipid_mass_fraction: 0.811,
            water_mass_fraction: 0.175,
            milk_solids_mass_fraction: 0.014,
            salt_nacl_mass_fraction: 0.0,
            t_ref_k: JANAF_REFERENCE_TEMPERATURE_K,
            rho_ref_kg_m3: 935.0,
            drho_dt_kg_m3_k: -0.52,
            cp_ref_j_kg_k: 2_020.0,
            dcp_dt_j_kg_k2: 1.45,
            k_ref_w_m_k: 0.214,
            dk_dt_w_m_k2: -0.000_25,
            mu_0_pa_s: 1.8e-6,
            e_visc_j_mol: 26_500.0,
            sigma_ref_n_m: 0.0325,
            dsigma_dt_n_m_k: -6.8e-5,
            t_smoke_fresh_k: 450.15,
            t_flash_k: 578.15,
            t_autoignition_k: 643.15,
            k_hyd_s_inv: 180_000.0,
            ea_hyd_j_mol: 82_000.0,
            k_ox_s_inv: 420_000.0,
            ea_ox_j_mol: 86_000.0,
            k_poly_s_inv: 150_000.0,
            ea_poly_j_mol: 92_000.0,
            beta_smoke_ffa_k: 1_200.0,
            beta_smoke_tpm_k: 110.0,
            gamma_tpm: 2.4,
            ffa_fresh: 0.0025,
            tpm_discard_limit: 0.25,
            ffa_discard_limit: 0.02,
            hydrolysis_contact_time_s: 50.0,
            oxidation_surface_factor_kg_m2: 7.5,
        }
    }

    /// Salted sweet-cream butter matching `data/matdb/kitchen-v1/oils/butter-salted`
    /// (`81.0%` butterfat, `16.0%` water, `1.4%` milk solids, `1.6%` NaCl).
    #[must_use]
    pub const fn butter_salted() -> Self {
        Self {
            lipid_mass_fraction: 0.810,
            water_mass_fraction: 0.160,
            milk_solids_mass_fraction: 0.014,
            salt_nacl_mass_fraction: 0.016,
            t_ref_k: JANAF_REFERENCE_TEMPERATURE_K,
            rho_ref_kg_m3: 942.0,
            drho_dt_kg_m3_k: -0.52,
            cp_ref_j_kg_k: 2_010.0,
            dcp_dt_j_kg_k2: 1.44,
            k_ref_w_m_k: 0.212,
            dk_dt_w_m_k2: -0.000_25,
            mu_0_pa_s: 1.9e-6,
            e_visc_j_mol: 26_600.0,
            sigma_ref_n_m: 0.0328,
            dsigma_dt_n_m_k: -6.8e-5,
            t_smoke_fresh_k: 450.15,
            t_flash_k: 578.15,
            t_autoignition_k: 643.15,
            k_hyd_s_inv: 170_000.0,
            ea_hyd_j_mol: 82_000.0,
            k_ox_s_inv: 450_000.0,
            ea_ox_j_mol: 85_500.0,
            k_poly_s_inv: 160_000.0,
            ea_poly_j_mol: 92_000.0,
            beta_smoke_ffa_k: 1_200.0,
            beta_smoke_tpm_k: 110.0,
            gamma_tpm: 2.4,
            ffa_fresh: 0.0025,
            tpm_discard_limit: 0.25,
            ffa_discard_limit: 0.02,
            hydrolysis_contact_time_s: 50.0,
            oxidation_surface_factor_kg_m2: 7.5,
        }
    }

    /// Clarified butter / ghee matching `data/matdb/kitchen-v1/oils/clarified-butter-ghee`.
    #[must_use]
    pub const fn clarified_butter_ghee() -> Self {
        Self {
            lipid_mass_fraction: 0.996,
            water_mass_fraction: 0.002,
            milk_solids_mass_fraction: 0.002,
            salt_nacl_mass_fraction: 0.0,
            t_ref_k: JANAF_REFERENCE_TEMPERATURE_K,
            rho_ref_kg_m3: 911.0,
            drho_dt_kg_m3_k: -0.66,
            cp_ref_j_kg_k: 2_020.0,
            dcp_dt_j_kg_k2: 2.10,
            k_ref_w_m_k: 0.172,
            dk_dt_w_m_k2: -0.000_24,
            mu_0_pa_s: 2.1e-6,
            e_visc_j_mol: 25_800.0,
            sigma_ref_n_m: 0.0315,
            dsigma_dt_n_m_k: -6.5e-5,
            t_smoke_fresh_k: 525.15,
            t_flash_k: 595.15,
            t_autoignition_k: 661.15,
            k_hyd_s_inv: 60_000.0,
            ea_hyd_j_mol: 86_000.0,
            k_ox_s_inv: 240_000.0,
            ea_ox_j_mol: 90_000.0,
            k_poly_s_inv: 120_000.0,
            ea_poly_j_mol: 96_000.0,
            beta_smoke_ffa_k: 1_200.0,
            beta_smoke_tpm_k: 115.0,
            gamma_tpm: 2.5,
            ffa_fresh: 0.001,
            tpm_discard_limit: 0.25,
            ffa_discard_limit: 0.02,
            hydrolysis_contact_time_s: 50.0,
            oxidation_surface_factor_kg_m2: 7.5,
        }
    }

    /// Rendered pork lard / beef tallow matching `data/matdb/kitchen-v1/oils/fat-pork-lard-tallow`.
    #[must_use]
    pub const fn pork_lard_tallow() -> Self {
        Self {
            lipid_mass_fraction: 0.997,
            water_mass_fraction: 0.003,
            milk_solids_mass_fraction: 0.0,
            salt_nacl_mass_fraction: 0.0,
            t_ref_k: JANAF_REFERENCE_TEMPERATURE_K,
            rho_ref_kg_m3: 915.0,
            drho_dt_kg_m3_k: -0.67,
            cp_ref_j_kg_k: 2_005.0,
            dcp_dt_j_kg_k2: 2.05,
            k_ref_w_m_k: 0.174,
            dk_dt_w_m_k2: -0.000_23,
            mu_0_pa_s: 1.9e-6,
            e_visc_j_mol: 26_800.0,
            sigma_ref_n_m: 0.0319,
            dsigma_dt_n_m_k: -6.5e-5,
            t_smoke_fresh_k: 463.15,
            t_flash_k: 588.15,
            t_autoignition_k: 653.15,
            k_hyd_s_inv: 100_000.0,
            ea_hyd_j_mol: 85_000.0,
            k_ox_s_inv: 250_000.0,
            ea_ox_j_mol: 90_000.0,
            k_poly_s_inv: 140_000.0,
            ea_poly_j_mol: 96_000.0,
            beta_smoke_ffa_k: 1_200.0,
            beta_smoke_tpm_k: 115.0,
            gamma_tpm: 2.4,
            ffa_fresh: 0.002,
            tpm_discard_limit: 0.25,
            ffa_discard_limit: 0.02,
            hydrolysis_contact_time_s: 50.0,
            oxidation_surface_factor_kg_m2: 7.5,
        }
    }

    /// Load a validated cooking-oil / butter card from an immutable `fs-matdb`
    /// [`ConstitutiveModelCard`] with law `"kitchen-cooking-oil"` v1.
    pub fn from_cooking_oil_card(card: &ConstitutiveModelCard) -> Result<Self, VesselErrorV1> {
        card.validate().map_err(|_| VesselErrorV1::InvalidCard {
            context: "kitchen-cooking-oil",
            reason: "fs-matdb card failed validation",
        })?;
        if card.law.0 != COOKING_OIL_LAW_ID_V1 {
            return Err(VesselErrorV1::InvalidCard {
                context: "law_id",
                reason: "expected law id `kitchen-cooking-oil`",
            });
        }
        if card.law_version != COOKING_OIL_LAW_VERSION_V1 {
            return Err(VesselErrorV1::InvalidCard {
                context: "law_version",
                reason: "expected law version 1",
            });
        }
        let get = |name: &'static str, expected_dims: Dims| -> Result<f64, VesselErrorV1> {
            let param = card
                .parameters
                .get(name)
                .ok_or(VesselErrorV1::InvalidCard {
                    context: name,
                    reason: "missing required oil card parameter",
                })?;
            if param.dims != expected_dims {
                return Err(VesselErrorV1::InvalidCard {
                    context: name,
                    reason: "oil card parameter SI dimensions mismatch",
                });
            }
            if !param.value.is_finite() {
                return Err(VesselErrorV1::InvalidCard {
                    context: name,
                    reason: "oil card parameter must be finite",
                });
            }
            Ok(param.value)
        };

        let beta_ffa_slope = get("beta_ffa", TEMP_DIMS)?;
        let beta_tpm_slope = get("beta_tpm", TEMP_DIMS)?;
        if beta_ffa_slope > 0.0 || beta_tpm_slope > 0.0 {
            return Err(VesselErrorV1::InvalidCard {
                context: "beta_ffa/beta_tpm",
                reason: "kitchen-cooking-oil smoke-point slopes must be non-positive (dT_smoke/dC <= 0)",
            });
        }

        let parsed = Self {
            lipid_mass_fraction: get("lipid_mass_fraction", DIMLESS)?,
            water_mass_fraction: get("water_mass_fraction", DIMLESS)?,
            milk_solids_mass_fraction: get("milk_solids_mass_fraction", DIMLESS)?,
            salt_nacl_mass_fraction: get("salt_nacl_mass_fraction", DIMLESS)?,
            t_ref_k: get("t_ref", TEMP_DIMS)?,
            rho_ref_kg_m3: get("rho_ref", DENSITY_DIMS)?,
            drho_dt_kg_m3_k: get("drho_dt", DRHO_DT_DIMS)?,
            cp_ref_j_kg_k: get("cp_ref", CP_DIMS)?,
            dcp_dt_j_kg_k2: get("dcp_dt", DCP_DT_DIMS)?,
            k_ref_w_m_k: get("k_ref", K_DIMS)?,
            dk_dt_w_m_k2: get("dk_dt", DK_DT_DIMS)?,
            mu_0_pa_s: get("mu_0", VISC_DIMS)?,
            e_visc_j_mol: get("e_visc", MOLAR_ENERGY_DIMS)?,
            sigma_ref_n_m: get("sigma_ref", SURF_TENS_DIMS)?,
            dsigma_dt_n_m_k: get("dsigma_dt", DSIGMA_DT_DIMS)?,
            t_smoke_fresh_k: get("t_smoke_fresh", TEMP_DIMS)?,
            t_flash_k: get("t_flash", TEMP_DIMS)?,
            t_autoignition_k: get("t_autoignition", TEMP_DIMS)?,
            k_hyd_s_inv: get("k_hyd", RATE_DIMS)?,
            ea_hyd_j_mol: get("ea_hyd", MOLAR_ENERGY_DIMS)?,
            k_ox_s_inv: get("k_ox", RATE_DIMS)?,
            ea_ox_j_mol: get("ea_ox", MOLAR_ENERGY_DIMS)?,
            k_poly_s_inv: get("k_poly", RATE_DIMS)?,
            ea_poly_j_mol: get("ea_poly", MOLAR_ENERGY_DIMS)?,
            beta_smoke_ffa_k: -beta_ffa_slope,
            beta_smoke_tpm_k: -beta_tpm_slope,
            gamma_tpm: get("gamma_tpm", DIMLESS)?,
            ffa_fresh: get("ffa_fresh", DIMLESS)?,
            tpm_discard_limit: get("tpm_discard_limit", DIMLESS)?,
            ffa_discard_limit: 0.02,
            hydrolysis_contact_time_s: 50.0,
            oxidation_surface_factor_kg_m2: 7.5,
        };
        parsed.validate()?;
        Ok(parsed)
    }

    /// Validate all physical and safety bounds on the oil/butter card.
    pub fn validate(&self) -> Result<(), VesselErrorV1> {
        let fracs = [
            self.lipid_mass_fraction,
            self.water_mass_fraction,
            self.milk_solids_mass_fraction,
            self.salt_nacl_mass_fraction,
        ];
        if !fracs
            .iter()
            .all(|v| v.is_finite() && (0.0..=1.0).contains(v))
        {
            return Err(VesselErrorV1::InvalidCard {
                context: "composition",
                reason: "oil/butter mass fractions must be finite and in [0, 1]",
            });
        }
        let comp_sum: f64 = fracs.iter().sum();
        if self.lipid_mass_fraction <= 0.0 || (comp_sum - 1.0).abs() > 1e-3 {
            return Err(VesselErrorV1::InvalidCard {
                context: "composition_sum",
                reason: "oil/butter lipid fraction must be positive and mass fractions must sum to 1",
            });
        }
        let strict_pos = [
            self.t_ref_k,
            self.rho_ref_kg_m3,
            self.cp_ref_j_kg_k,
            self.k_ref_w_m_k,
            self.mu_0_pa_s,
            self.e_visc_j_mol,
            self.sigma_ref_n_m,
            self.t_smoke_fresh_k,
            self.t_flash_k,
            self.t_autoignition_k,
            self.ea_hyd_j_mol,
            self.ea_ox_j_mol,
            self.ea_poly_j_mol,
            self.tpm_discard_limit,
            self.ffa_discard_limit,
        ];
        if !strict_pos.iter().all(|v| v.is_finite() && *v > 0.0) {
            return Err(VesselErrorV1::InvalidCard {
                context: "positive_scalars",
                reason: "thermophysical, activation, and safety temperatures must be finite positive",
            });
        }
        let non_neg = [
            self.dcp_dt_j_kg_k2,
            self.k_hyd_s_inv,
            self.k_ox_s_inv,
            self.k_poly_s_inv,
            self.beta_smoke_ffa_k,
            self.beta_smoke_tpm_k,
            self.gamma_tpm,
            self.ffa_fresh,
            self.hydrolysis_contact_time_s,
            self.oxidation_surface_factor_kg_m2,
        ];
        if !non_neg.iter().all(|v| v.is_finite() && *v >= 0.0) {
            return Err(VesselErrorV1::InvalidCard {
                context: "non_negative_scalars",
                reason: "degradation coefficients and slopes must be finite non-negative",
            });
        }
        if ![
            self.drho_dt_kg_m3_k,
            self.dk_dt_w_m_k2,
            self.dsigma_dt_n_m_k,
        ]
        .iter()
        .all(|v| v.is_finite())
        {
            return Err(VesselErrorV1::InvalidCard {
                context: "temperature_slopes",
                reason: "property temperature slopes must be finite",
            });
        }
        if !(self.t_smoke_fresh_k < self.t_flash_k && self.t_flash_k < self.t_autoignition_k) {
            return Err(VesselErrorV1::InvalidCard {
                context: "safety_ordering",
                reason: "require T_smoke_fresh < T_flash < T_autoignition",
            });
        }
        Ok(())
    }

    /// Evaluate oil density `rho(T)` [kg/m^3], safely bounded positive.
    #[must_use]
    pub fn density_kg_m3(&self, temperature_k: f64) -> f64 {
        (self.rho_ref_kg_m3 + self.drho_dt_kg_m3_k * (temperature_k - self.t_ref_k)).max(400.0)
    }

    /// Evaluate oil viscosity `mu(T, TPM) = mu_0 * exp(E_visc / (R T)) * exp(gamma_TPM * TPM)` [Pa s].
    #[must_use]
    pub fn viscosity_pa_s(&self, temperature_k: f64, tpm_fraction: f64) -> f64 {
        let t = temperature_k.clamp(250.0, 900.0);
        let tpm = tpm_fraction.clamp(0.0, 1.0);
        self.mu_0_pa_s
            * exp(self.e_visc_j_mol / (UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * t))
            * exp(self.gamma_tpm * tpm)
    }

    /// Evaluate degraded smoke point `T_smoke = T_smoke,0 - beta_FFA * [FFA] - beta_TPM * [TPM]` [K].
    #[must_use]
    pub fn smoke_point_k(&self, ffa_fraction: f64, tpm_fraction: f64) -> f64 {
        (self.t_smoke_fresh_k
            - self.beta_smoke_ffa_k * ffa_fraction.max(0.0)
            - self.beta_smoke_tpm_k * tpm_fraction.max(0.0))
        .max(300.0)
    }

    fn append_identity_words(&self, words: &mut Vec<u64>) {
        let scalars = [
            self.lipid_mass_fraction,
            self.water_mass_fraction,
            self.milk_solids_mass_fraction,
            self.salt_nacl_mass_fraction,
            self.t_ref_k,
            self.rho_ref_kg_m3,
            self.drho_dt_kg_m3_k,
            self.cp_ref_j_kg_k,
            self.dcp_dt_j_kg_k2,
            self.k_ref_w_m_k,
            self.dk_dt_w_m_k2,
            self.mu_0_pa_s,
            self.e_visc_j_mol,
            self.sigma_ref_n_m,
            self.dsigma_dt_n_m_k,
            self.t_smoke_fresh_k,
            self.t_flash_k,
            self.t_autoignition_k,
            self.k_hyd_s_inv,
            self.ea_hyd_j_mol,
            self.k_ox_s_inv,
            self.ea_ox_j_mol,
            self.k_poly_s_inv,
            self.ea_poly_j_mol,
            self.beta_smoke_ffa_k,
            self.beta_smoke_tpm_k,
            self.gamma_tpm,
            self.ffa_fresh,
            self.tpm_discard_limit,
            self.ffa_discard_limit,
            self.hydrolysis_contact_time_s,
            self.oxidation_surface_factor_kg_m2,
        ];
        words.extend(scalars.iter().map(|v| v.to_bits()));
    }
}

/// Complete physical and geometric configuration for [`MultiphaseVesselComponent`].
#[derive(Clone, Debug, PartialEq)]
pub struct MultiphaseVesselConfig {
    /// Heated bottom wetting area `A_base` [m^2].
    pub base_area_m2: f64,
    /// Open free-surface / oil-air interface area `A_air` [m^2].
    pub free_surface_area_m2: f64,
    /// Internal vessel brim height `H_pot` [m] (for Bikerman boil-over detection).
    pub vessel_height_m: f64,
    /// Number of immersed food batch channels on `immersed_food_heat` (`M >= 1`).
    pub immersed_batch_count: usize,
    /// Dryout wetting regularization film thickness `delta_min` [m] (default `50e-6 m`).
    pub dryout_delta_min_m: f64,
    /// Volatile mass floor below which dryout clamps cleanly instead of shortening [kg].
    pub dryout_floor_kg: f64,
    /// Smallest step size `dt_min` [s] offered by `Proposal::Shorten`.
    pub min_dt_s: f64,
    /// Optional dry vessel liner residual heat capacity [J/K] (`0.0` for pure-fluid accounting).
    pub dry_vessel_heat_capacity_j_k: f64,
    /// Pan-floor fond contact area `A_fond` [m^2].
    pub fond_area_m2: f64,
    /// Fond dissolution reference mass-transfer coefficient `k_diss,0` at `80 C` [m/s].
    pub fond_k_diss_ref_m_s: f64,
    /// Fond dissolution activation energy `Ea_diss` [J/mol].
    pub fond_ea_diss_j_mol: f64,
    /// Fond aqueous saturation concentration `C_sat,ref` at `25 C` [kg/m^3].
    pub fond_c_sat_ref_kg_m3: f64,
    /// Bikerman foam relaxation time constant `tau_foam` [s].
    pub foam_tau_s: f64,
    /// Bikerman starch foaminess ceiling `Sigma_starch` [s].
    pub foam_sigma_starch_s: f64,
    /// Bikerman half-saturation starch concentration `c_starch,1/2` [kg/m^3].
    pub foam_c_starch_half_kg_m3: f64,
    /// Bikerman protein/milk-solid foaminess ceiling `Sigma_protein` [s].
    pub foam_sigma_protein_s: f64,
    /// Bikerman half-saturation protein concentration `c_protein,1/2` [kg/m^3].
    pub foam_c_protein_half_kg_m3: f64,
    /// Carreau/Herschel-Bulkley aqueous solids thickening exponential factor `k_conc`.
    pub sauce_viscosity_conc_coeff: f64,
    /// Carreau shear-thinning flow index `n_flow in (0.1, 1.0]`.
    pub sauce_carreau_n: f64,
    /// Carreau relaxation time `lambda_carreau` [s].
    pub sauce_carreau_lambda_s: f64,
    /// Natural-convection baseline wall heat-transfer coefficient `h_wall,0` [W/(m^2 K)].
    pub natural_wall_htc_w_m2_k: f64,
    /// Natural-convection baseline immersed-food HTC `h_food,0` [W/(m^2 K)].
    pub natural_immersed_htc_w_m2_k: f64,
    /// Shear-stirring HTC gain coefficient `c_stir` (gives 3–5x curd-mixing boost at 120–180 RPM).
    pub stir_htc_gain: f64,
    /// Vapor-bubble agitation HTC gain coefficient `c_bubble` for frying/boiling immersion.
    pub bubble_htc_gain: f64,
    /// Reference superficial vapor velocity `j_g,ref` [m/s] for bubble-agitation HTC boost.
    pub bubble_jg_ref_m_s: f64,
    /// Wall-film scorching warning threshold temperature [K] (`~418.15 K` / `145 C`).
    pub scorch_wall_threshold_k: f64,
    /// Enable deep-fryer 2-zone sub-tube cold zone thermal stratification.
    pub cold_zone_enabled: bool,
    /// Stratified thermal conductance `UA_coldzone` [W/K] between active vat and sub-tube cold zone.
    pub cold_zone_ua_w_k: f64,
    /// Enable stoichiometric triglyceride oil hydrolysis (`triglyceride + 3 H2O -> 3 ffa + glycerol`).
    pub hydrolysis_enabled: bool,
    /// Optional cooking-oil / butter medium card.
    pub oil: Option<OilCardV1>,
}

impl MultiphaseVesselConfig {
    /// Standard commercial saucepan (`D = 20 cm`, `H = 12 cm`, wine deglazing / sauce reduction / eggs).
    #[must_use]
    pub const fn saucepan_default() -> Self {
        Self {
            base_area_m2: 0.031_416,
            free_surface_area_m2: 0.031_416,
            vessel_height_m: 0.12,
            immersed_batch_count: 1,
            dryout_delta_min_m: 50.0e-6,
            dryout_floor_kg: 1.0e-6,
            min_dt_s: 1.0e-4,
            dry_vessel_heat_capacity_j_k: 0.0,
            fond_area_m2: 0.031_416,
            fond_k_diss_ref_m_s: 2.5e-4,
            fond_ea_diss_j_mol: 28_000.0,
            fond_c_sat_ref_kg_m3: 180.0,
            foam_tau_s: 8.0,
            foam_sigma_starch_s: 20.0,
            foam_c_starch_half_kg_m3: 4.0,
            foam_sigma_protein_s: 16.0,
            foam_c_protein_half_kg_m3: 3.0,
            sauce_viscosity_conc_coeff: 8.5,
            sauce_carreau_n: 0.55,
            sauce_carreau_lambda_s: 0.8,
            natural_wall_htc_w_m2_k: 450.0,
            natural_immersed_htc_w_m2_k: 400.0,
            stir_htc_gain: 1.85,
            bubble_htc_gain: 2.5,
            bubble_jg_ref_m_s: 0.015,
            scorch_wall_threshold_k: 418.15,
            cold_zone_enabled: false,
            cold_zone_ua_w_k: 0.0,
            hydrolysis_enabled: false,
            oil: None,
        }
    }

    /// Standard commercial pasta / stockpot (`D = 26 cm`, `H = 22 cm`).
    #[must_use]
    pub const fn stockpot_default() -> Self {
        Self {
            base_area_m2: 0.053_093,
            free_surface_area_m2: 0.053_093,
            vessel_height_m: 0.22,
            fond_area_m2: 0.053_093,
            ..Self::saucepan_default()
        }
    }

    /// Wide sauté skillet (`D = 26 cm`, `H = 5.5 cm`) for pan deglazing, searing/fond, or scrambled eggs.
    #[must_use]
    pub const fn skillet_default() -> Self {
        Self {
            base_area_m2: 0.053_093,
            free_surface_area_m2: 0.053_093,
            vessel_height_m: 0.055,
            fond_area_m2: 0.053_093,
            ..Self::saucepan_default()
        }
    }

    /// Butter / *beurre noisette* sauté skillet configured with an oil or butter card.
    #[must_use]
    pub const fn butter_skillet_default(oil: OilCardV1) -> Self {
        Self {
            base_area_m2: 0.045,
            free_surface_area_m2: 0.045,
            vessel_height_m: 0.055,
            fond_area_m2: 0.045,
            natural_wall_htc_w_m2_k: 320.0,
            natural_immersed_htc_w_m2_k: 280.0,
            hydrolysis_enabled: false,
            oil: Some(oil),
            ..Self::saucepan_default()
        }
    }

    /// Commercial twin-basket deep fryer (`0.40 m x 0.40 m` vat, `2` immersed basket channels,
    /// stably-stratified sub-tube cold zone `UA = 0.20 W/K`, and active stoichiometric oil hydrolysis).
    #[must_use]
    pub const fn deep_fryer_default(oil: OilCardV1) -> Self {
        Self {
            base_area_m2: 0.16,
            free_surface_area_m2: 0.16,
            vessel_height_m: 0.30,
            immersed_batch_count: 2,
            fond_area_m2: 0.16,
            natural_wall_htc_w_m2_k: 260.0,
            natural_immersed_htc_w_m2_k: 250.0,
            bubble_htc_gain: 2.8,
            bubble_jg_ref_m_s: 0.012,
            cold_zone_enabled: true,
            cold_zone_ua_w_k: 0.20,
            hydrolysis_enabled: true,
            oil: Some(oil),
            ..Self::saucepan_default()
        }
    }

    /// Validate every geometric, numerical, and material configuration parameter fail-closed.
    pub fn validate(&self) -> Result<(), VesselErrorV1> {
        let strict_pos = [
            ("base_area_m2", self.base_area_m2),
            ("free_surface_area_m2", self.free_surface_area_m2),
            ("vessel_height_m", self.vessel_height_m),
            ("dryout_delta_min_m", self.dryout_delta_min_m),
            ("dryout_floor_kg", self.dryout_floor_kg),
            ("min_dt_s", self.min_dt_s),
            ("fond_c_sat_ref_kg_m3", self.fond_c_sat_ref_kg_m3),
            ("foam_tau_s", self.foam_tau_s),
            ("foam_c_starch_half_kg_m3", self.foam_c_starch_half_kg_m3),
            ("foam_c_protein_half_kg_m3", self.foam_c_protein_half_kg_m3),
            ("sauce_carreau_lambda_s", self.sauce_carreau_lambda_s),
            ("natural_wall_htc_w_m2_k", self.natural_wall_htc_w_m2_k),
            (
                "natural_immersed_htc_w_m2_k",
                self.natural_immersed_htc_w_m2_k,
            ),
            ("bubble_jg_ref_m_s", self.bubble_jg_ref_m_s),
            ("scorch_wall_threshold_k", self.scorch_wall_threshold_k),
        ];
        for (parameter, value) in strict_pos {
            if !value.is_finite() || value <= 0.0 {
                return Err(VesselErrorV1::InvalidConfig {
                    parameter,
                    reason: "must be finite and strictly positive",
                });
            }
        }
        let non_neg = [
            (
                "dry_vessel_heat_capacity_j_k",
                self.dry_vessel_heat_capacity_j_k,
            ),
            ("fond_area_m2", self.fond_area_m2),
            ("fond_k_diss_ref_m_s", self.fond_k_diss_ref_m_s),
            ("fond_ea_diss_j_mol", self.fond_ea_diss_j_mol),
            ("foam_sigma_starch_s", self.foam_sigma_starch_s),
            ("foam_sigma_protein_s", self.foam_sigma_protein_s),
            (
                "sauce_viscosity_conc_coeff",
                self.sauce_viscosity_conc_coeff,
            ),
            ("stir_htc_gain", self.stir_htc_gain),
            ("bubble_htc_gain", self.bubble_htc_gain),
            ("cold_zone_ua_w_k", self.cold_zone_ua_w_k),
        ];
        for (parameter, value) in non_neg {
            if !value.is_finite() || value < 0.0 {
                return Err(VesselErrorV1::InvalidConfig {
                    parameter,
                    reason: "must be finite and non-negative",
                });
            }
        }
        if self.immersed_batch_count == 0 {
            return Err(VesselErrorV1::InvalidConfig {
                parameter: "immersed_batch_count",
                reason: "must be >= 1",
            });
        }
        if !self.sauce_carreau_n.is_finite()
            || self.sauce_carreau_n <= 0.1
            || self.sauce_carreau_n > 1.0
        {
            return Err(VesselErrorV1::InvalidConfig {
                parameter: "sauce_carreau_n",
                reason: "Carreau flow index must lie in (0.1, 1.0]",
            });
        }
        if let Some(oil) = &self.oil {
            oil.validate()?;
        } else if self.hydrolysis_enabled {
            return Err(VesselErrorV1::InvalidConfig {
                parameter: "hydrolysis_enabled",
                reason: "oil hydrolysis requires an attached OilCardV1",
            });
        }
        Ok(())
    }

    fn lipid_cp_params(&self) -> (f64, f64, f64) {
        if let Some(oil) = &self.oil {
            (oil.cp_ref_j_kg_k, oil.dcp_dt_j_kg_k2, oil.t_ref_k)
        } else {
            (
                CP_TRIG_DEFAULT_J_KG_K,
                DCP_DT_TRIG_DEFAULT_J_KG_K2,
                JANAF_REFERENCE_TEMPERATURE_K,
            )
        }
    }
}

/// Structured initial inventory descriptor converted into a thermodynamically
/// consistent [`MultiphaseVesselState`] on `janaf-298.15K-formation`.
#[derive(Clone, Debug, PartialEq)]
pub struct VesselInitialStateV1 {
    /// Initial active liquid temperature [K].
    pub t_liq_k: f64,
    /// Initial sub-tube cold-zone temperature [K] (ignored unless `cold_zone_enabled`).
    pub t_coldzone_k: f64,
    /// Liquid water mass [kg].
    pub m_h2o_liq_kg: f64,
    /// Liquid ethanol mass [kg].
    pub m_etoh_liq_kg: f64,
    /// Active liquid triglyceride oil mass [kg].
    pub m_triglyceride_kg: f64,
    /// Active liquid free fatty acid mass [kg].
    pub m_ffa_kg: f64,
    /// Active liquid glycerol mass [kg].
    pub m_glycerol_kg: f64,
    /// Dissolved salt (`NaCl`) mass [kg].
    pub m_salt_kg: f64,
    /// Dissolved sugar (`sucrose`) mass [kg].
    pub m_sugar_kg: f64,
    /// Dissolved fond solute mass [kg].
    pub m_fond_dissolved_kg: f64,
    /// Undissolved adhered pan fond mass [kg] (included in `stored_mass_kg` and `solute_fond`).
    pub m_fond_adhered_kg: f64,
    /// Suspended butter milk-solids mass [kg].
    pub m_milk_solids_kg: f64,
    /// Suspended starch mass [kg].
    pub m_starch_kg: f64,
    /// Suspended protein mass [kg].
    pub m_protein_kg: f64,
    /// Sub-tube cold-zone triglyceride oil mass [kg].
    pub m_oil_coldzone_kg: f64,
    /// Initial headspace water vapor mass [kg].
    pub m_h2o_vap_kg: f64,
    /// Initial headspace ethanol vapor mass [kg].
    pub m_etoh_vap_kg: f64,
    /// Initial auxiliary headspace gases `[CO2, O2, N2, CH4, C3H8, CO, grease_vapor]` [kg].
    pub aux_gas_kg: [f64; AUX_GAS_COUNT],
    /// Initial oil TPM mass fraction (`0.0` for fresh oil).
    pub tpm_fraction: f64,
}

impl VesselInitialStateV1 {
    /// Pure water pot/pan initial state.
    #[must_use]
    pub const fn water(mass_kg: f64, temperature_k: f64) -> Self {
        Self {
            t_liq_k: temperature_k,
            t_coldzone_k: temperature_k,
            m_h2o_liq_kg: mass_kg,
            m_etoh_liq_kg: 0.0,
            m_triglyceride_kg: 0.0,
            m_ffa_kg: 0.0,
            m_glycerol_kg: 0.0,
            m_salt_kg: 0.0,
            m_sugar_kg: 0.0,
            m_fond_dissolved_kg: 0.0,
            m_fond_adhered_kg: 0.0,
            m_milk_solids_kg: 0.0,
            m_starch_kg: 0.0,
            m_protein_kg: 0.0,
            m_oil_coldzone_kg: 0.0,
            m_h2o_vap_kg: 0.0,
            m_etoh_vap_kg: 0.0,
            aux_gas_kg: [0.0; AUX_GAS_COUNT],
            tpm_fraction: 0.0,
        }
    }

    /// Wine deglazing pan with adhered fond and ethanol-water liquid.
    #[must_use]
    pub fn wine_deglaze(
        wine_mass_kg: f64,
        ethanol_mass_frac: f64,
        fond_adhered_kg: f64,
        temperature_k: f64,
    ) -> Self {
        let mut s = Self::water(wine_mass_kg * (1.0 - ethanol_mass_frac), temperature_k);
        s.m_etoh_liq_kg = wine_mass_kg * ethanol_mass_frac;
        s.m_fond_adhered_kg = fond_adhered_kg;
        s
    }

    /// Starchy salted pasta-boiling pot initial state.
    #[must_use]
    pub const fn pasta_pot(
        water_mass_kg: f64,
        starch_mass_kg: f64,
        salt_mass_kg: f64,
        temperature_k: f64,
    ) -> Self {
        let mut s = Self::water(water_mass_kg, temperature_k);
        s.m_starch_kg = starch_mass_kg;
        s.m_salt_kg = salt_mass_kg;
        s
    }

    /// Reducing pan sauce / stock with dissolved fond, starch, salt, and sugar.
    #[must_use]
    pub const fn sauce(
        water_mass_kg: f64,
        fond_kg: f64,
        starch_kg: f64,
        salt_kg: f64,
        sugar_kg: f64,
        temperature_k: f64,
    ) -> Self {
        let mut s = Self::water(water_mass_kg, temperature_k);
        s.m_fond_dissolved_kg = fond_kg;
        s.m_starch_kg = starch_kg;
        s.m_salt_kg = salt_kg;
        s.m_sugar_kg = sugar_kg;
        s
    }

    /// Commercial deep-fryer vat with active oil, optional trace moisture, fresh FFA,
    /// and sub-tube cold-zone oil.
    #[must_use]
    pub fn oil_fryer(
        oil: &OilCardV1,
        active_oil_mass_kg: f64,
        coldzone_oil_mass_kg: f64,
        active_temp_k: f64,
        coldzone_temp_k: f64,
    ) -> Self {
        let ffa_kg = active_oil_mass_kg * oil.lipid_mass_fraction * oil.ffa_fresh;
        let trig_kg = (active_oil_mass_kg * oil.lipid_mass_fraction - ffa_kg).max(0.0);
        Self {
            t_liq_k: active_temp_k,
            t_coldzone_k: coldzone_temp_k,
            m_h2o_liq_kg: active_oil_mass_kg * oil.water_mass_fraction,
            m_etoh_liq_kg: 0.0,
            m_triglyceride_kg: trig_kg,
            m_ffa_kg: ffa_kg,
            m_glycerol_kg: 0.0,
            m_salt_kg: active_oil_mass_kg * oil.salt_nacl_mass_fraction,
            m_sugar_kg: 0.0,
            m_fond_dissolved_kg: 0.0,
            m_fond_adhered_kg: 0.0,
            m_milk_solids_kg: active_oil_mass_kg * oil.milk_solids_mass_fraction,
            m_starch_kg: 0.0,
            m_protein_kg: 0.0,
            m_oil_coldzone_kg: coldzone_oil_mass_kg,
            m_h2o_vap_kg: 0.0,
            m_etoh_vap_kg: 0.0,
            aux_gas_kg: [0.0; AUX_GAS_COUNT],
            tpm_fraction: 0.0,
        }
    }

    /// Butter, clarified butter, or pan-fat initial state built directly from an [`OilCardV1`].
    #[must_use]
    pub fn butter_pan(butter: &OilCardV1, total_mass_kg: f64, temperature_k: f64) -> Self {
        Self::oil_fryer(butter, total_mass_kg, 0.0, temperature_k, temperature_k)
    }
}

/// Complete accepted/candidate state owned by `fs-system` for [`MultiphaseVesselComponent`].
#[derive(Clone, Debug, PartialEq)]
pub struct MultiphaseVesselState {
    /// Elapsed simulation time [s].
    pub time_s: f64,
    /// Active bulk liquid/bath temperature [K].
    pub t_liq_k: f64,
    /// Sub-tube cold-zone temperature [K] (equals `t_liq_k` when disabled).
    pub t_coldzone_k: f64,
    /// Total active-zone stored enthalpy on `janaf-298.15K-formation` [J].
    pub active_enthalpy_j: f64,
    /// Total cold-zone stored enthalpy on `janaf-298.15K-formation` [J].
    pub coldzone_enthalpy_j: f64,
    /// Active liquid water mass [kg].
    pub m_h2o_liq_kg: f64,
    /// Active liquid ethanol mass [kg].
    pub m_etoh_liq_kg: f64,
    /// Active triglyceride oil mass [kg].
    pub m_triglyceride_kg: f64,
    /// Active free fatty acid mass [kg].
    pub m_ffa_kg: f64,
    /// Active glycerol mass [kg].
    pub m_glycerol_kg: f64,
    /// Dissolved salt (`NaCl`) mass [kg].
    pub m_salt_kg: f64,
    /// Dissolved sugar (`sucrose`) mass [kg].
    pub m_sugar_kg: f64,
    /// Dissolved fond solute mass [kg].
    pub m_fond_dissolved_kg: f64,
    /// Undissolved adhered pan fond mass [kg] (included in total mass & `solute_fond`).
    pub m_fond_adhered_kg: f64,
    /// Suspended butter milk-solids mass [kg].
    pub m_milk_solids_kg: f64,
    /// Suspended starch mass [kg].
    pub m_starch_kg: f64,
    /// Suspended protein mass [kg].
    pub m_protein_kg: f64,
    /// Sub-tube cold-zone triglyceride oil mass [kg].
    pub m_oil_coldzone_kg: f64,
    /// Headspace water vapor mass [kg].
    pub m_h2o_vap_kg: f64,
    /// Headspace ethanol vapor mass [kg].
    pub m_etoh_vap_kg: f64,
    /// Headspace auxiliary gas species `[CO2, O2, N2, CH4, C3H8, CO, grease_vapor]` [kg].
    pub aux_gas_kg: [f64; AUX_GAS_COUNT],
    /// Oil total polar materials quality fraction `[0, 1]`.
    pub tpm_fraction: f64,
    /// Bikerman foam layer height [m].
    pub foam_height_m: f64,
    /// Milk-solid *beurre noisette* Maillard browning progress `[0, 1]`.
    pub milk_solid_browning: f64,
    /// Cumulative molar extent of triglyceride hydrolysis `triglyceride + 3 H2O -> 3 ffa + glycerol` [mol].
    pub cumulative_hydrolysis_mol: f64,
}

impl MultiphaseVesselState {
    /// Total active free liquid mass [kg] (excluding adhered fond, cold zone, and vapor headspace).
    #[must_use]
    pub fn active_liquid_mass_kg(&self) -> f64 {
        self.m_h2o_liq_kg
            + self.m_etoh_liq_kg
            + self.m_triglyceride_kg
            + self.m_ffa_kg
            + self.m_glycerol_kg
            + self.m_salt_kg
            + self.m_sugar_kg
            + self.m_fond_dissolved_kg
            + self.m_milk_solids_kg
            + self.m_starch_kg
            + self.m_protein_kg
    }

    /// Total free liquid mass [kg] (active + cold-zone liquid).
    #[must_use]
    pub fn total_liquid_mass_kg(&self) -> f64 {
        self.active_liquid_mass_kg() + self.m_oil_coldzone_kg
    }

    /// Total headspace vapor + gas mass [kg].
    #[must_use]
    pub fn total_vapor_mass_kg(&self) -> f64 {
        self.m_h2o_vap_kg + self.m_etoh_vap_kg + self.aux_gas_kg.iter().sum::<f64>()
    }

    /// Total volatile liquid mass (`H2O + C2H5OH`) [kg].
    #[must_use]
    pub fn volatile_liquid_mass_kg(&self) -> f64 {
        self.m_h2o_liq_kg + self.m_etoh_liq_kg
    }

    /// Total active lipid mass (`triglyceride + ffa + glycerol`) [kg].
    #[must_use]
    pub fn active_lipid_mass_kg(&self) -> f64 {
        self.m_triglyceride_kg + self.m_ffa_kg + self.m_glycerol_kg
    }

    /// Total stored vessel mass [kg] across liquid, cold zone, adhered fond, and headspace vapor.
    #[must_use]
    pub fn total_mass_kg(&self) -> f64 {
        self.active_liquid_mass_kg()
            + self.m_fond_adhered_kg
            + self.m_oil_coldzone_kg
            + self.total_vapor_mass_kg()
    }

    /// Total stored enthalpy on `janaf-298.15K-formation` [J].
    #[must_use]
    pub fn total_enthalpy_j(&self) -> f64 {
        self.active_enthalpy_j + self.coldzone_enthalpy_j
    }

    /// Liquid ethanol mole fraction `x_EtOH` in the volatile `H2O + C2H5OH` phase (`0.0` if empty).
    #[must_use]
    pub fn ethanol_mole_fraction(&self) -> f64 {
        let n_e = self.m_etoh_liq_kg / MOLAR_MASS_ETOH_KG_MOL;
        let n_w = self.m_h2o_liq_kg / MOLAR_MASS_H2O_KG_MOL;
        let n_tot = n_e + n_w;
        if n_tot > 1.0e-18 {
            (n_e / n_tot).clamp(0.0, 1.0)
        } else {
            0.0
        }
    }

    /// Non-volatile solids mass fraction in the active liquid phase.
    #[must_use]
    pub fn solids_mass_fraction(&self) -> f64 {
        let m_liq = self.active_liquid_mass_kg();
        if m_liq > 1.0e-15 {
            let m_solids = self.m_salt_kg
                + self.m_sugar_kg
                + self.m_fond_dissolved_kg
                + self.m_milk_solids_kg
                + self.m_starch_kg
                + self.m_protein_kg;
            (m_solids / m_liq).clamp(0.0, 1.0)
        } else {
            0.0
        }
    }

    /// Free fatty acid mass fraction `[FFA] = m_ffa / m_lipid_total` (`0.0` when no lipid).
    #[must_use]
    pub fn ffa_fraction(&self) -> f64 {
        let m_lip = self.active_lipid_mass_kg() + self.m_oil_coldzone_kg;
        if m_lip > 1.0e-15 {
            (self.m_ffa_kg / m_lip).clamp(0.0, 1.0)
        } else {
            0.0
        }
    }

    /// Validate all scalar invariants of this state fail-closed.
    #[allow(clippy::too_many_lines)]
    pub fn validate(&self, config: &MultiphaseVesselConfig) -> Result<(), VesselErrorV1> {
        let non_neg = [
            ("time_s", self.time_s),
            ("m_h2o_liq_kg", self.m_h2o_liq_kg),
            ("m_etoh_liq_kg", self.m_etoh_liq_kg),
            ("m_triglyceride_kg", self.m_triglyceride_kg),
            ("m_ffa_kg", self.m_ffa_kg),
            ("m_glycerol_kg", self.m_glycerol_kg),
            ("m_salt_kg", self.m_salt_kg),
            ("m_sugar_kg", self.m_sugar_kg),
            ("m_fond_dissolved_kg", self.m_fond_dissolved_kg),
            ("m_fond_adhered_kg", self.m_fond_adhered_kg),
            ("m_milk_solids_kg", self.m_milk_solids_kg),
            ("m_starch_kg", self.m_starch_kg),
            ("m_protein_kg", self.m_protein_kg),
            ("m_oil_coldzone_kg", self.m_oil_coldzone_kg),
            ("m_h2o_vap_kg", self.m_h2o_vap_kg),
            ("m_etoh_vap_kg", self.m_etoh_vap_kg),
            ("tpm_fraction", self.tpm_fraction),
            ("foam_height_m", self.foam_height_m),
            ("milk_solid_browning", self.milk_solid_browning),
            ("cumulative_hydrolysis_mol", self.cumulative_hydrolysis_mol),
        ];
        for (quantity, value) in non_neg {
            if !value.is_finite() || value < 0.0 {
                return Err(VesselErrorV1::InvalidScalar {
                    quantity,
                    value,
                    reason: "must be finite and non-negative",
                });
            }
        }
        for &m in &self.aux_gas_kg {
            if !m.is_finite() || m < 0.0 {
                return Err(VesselErrorV1::InvalidScalar {
                    quantity: "aux_gas_kg",
                    value: m,
                    reason: "must be finite and non-negative",
                });
            }
        }
        if !self.t_liq_k.is_finite()
            || !(TEMP_MIN_VALID_K..=TEMP_MAX_VALID_K).contains(&self.t_liq_k)
        {
            return Err(VesselErrorV1::InvalidScalar {
                quantity: "t_liq_k",
                value: self.t_liq_k,
                reason: "must lie inside [200, 1500] K",
            });
        }
        if !self.t_coldzone_k.is_finite()
            || !(TEMP_MIN_VALID_K..=TEMP_MAX_VALID_K).contains(&self.t_coldzone_k)
        {
            return Err(VesselErrorV1::InvalidScalar {
                quantity: "t_coldzone_k",
                value: self.t_coldzone_k,
                reason: "must lie inside [200, 1500] K",
            });
        }
        if !self.active_enthalpy_j.is_finite() || !self.coldzone_enthalpy_j.is_finite() {
            return Err(VesselErrorV1::InvalidScalar {
                quantity: "enthalpy_j",
                value: self.total_enthalpy_j(),
                reason: "must be finite",
            });
        }
        if self.tpm_fraction > 1.0 || self.milk_solid_browning > 1.0 {
            return Err(VesselErrorV1::InvalidScalar {
                quantity: "tpm_or_browning",
                value: self.tpm_fraction.max(self.milk_solid_browning),
                reason: "quality progress fraction must be <= 1.0",
            });
        }
        let expected_active_h = compute_active_enthalpy_j(self, config, self.t_liq_k);
        let active_scale = expected_active_h
            .abs()
            .max(self.active_enthalpy_j.abs())
            .max(1.0);
        if (expected_active_h - self.active_enthalpy_j).abs() > 1.0e-6 * active_scale {
            return Err(VesselErrorV1::InvalidScalar {
                quantity: "active_enthalpy_j",
                value: self.active_enthalpy_j,
                reason: "active enthalpy disagrees with state temperature and masses",
            });
        }
        let expected_cz_h = compute_coldzone_enthalpy_j(self, config, self.t_coldzone_k);
        let cz_scale = expected_cz_h
            .abs()
            .max(self.coldzone_enthalpy_j.abs())
            .max(1.0);
        if (expected_cz_h - self.coldzone_enthalpy_j).abs() > 1.0e-6 * cz_scale {
            return Err(VesselErrorV1::InvalidScalar {
                quantity: "coldzone_enthalpy_j",
                value: self.coldzone_enthalpy_j,
                reason: "cold-zone enthalpy disagrees with cold-zone temperature and mass",
            });
        }
        Ok(())
    }
}

/// Decoded safety, boil-over, oil-quality, and wall-scorching flags.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[allow(clippy::struct_excessive_bools)]
pub struct VesselSafetyFlagsV1 {
    /// Bikerman boil-over active (`H_liq + H_foam > H_pot`).
    pub boil_over: bool,
    /// Oil temperature reached or exceeded degraded smoke point (`T_liq >= T_smoke`).
    pub smoking: bool,
    /// Oil temperature reached or exceeded open-cup flashpoint (`T_liq >= T_flash`).
    pub flashpoint_exceeded: bool,
    /// Oil temperature reached or exceeded autoignition (`T_liq >= T_autoignition`).
    pub autoignition_exceeded: bool,
    /// Frying oil exceeded commercial TPM (`>= 25%`) or FFA (`>= 2%`) discard limit.
    pub oil_discard_exceeded: bool,
    /// Concentrated aqueous sauce wall-film temperature exceeded scorching threshold.
    pub scorch_risk: bool,
}

impl VesselSafetyFlagsV1 {
    /// Encode flags as a deterministic bitmask float for `probes[11]`.
    ///
    /// When only boil-over is active, this equals `1.0` identically.
    #[must_use]
    pub fn as_flag_f64(self) -> f64 {
        let mut mask = 0u32;
        if self.boil_over {
            mask |= 1;
        }
        if self.smoking {
            mask |= 2;
        }
        if self.flashpoint_exceeded {
            mask |= 4;
        }
        if self.autoignition_exceeded {
            mask |= 8;
        }
        if self.oil_discard_exceeded {
            mask |= 16;
        }
        if self.scorch_risk {
            mask |= 32;
        }
        f64::from(mask)
    }
}

/// Donor-rule headspace venting and boil-over overflow specification for open-vessel
/// or relief-valve coupled stepping via [`MultiphaseVesselComponent::open_vent_inputs`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OpenVentSpecV1 {
    /// Fraction of current headspace vapor vented through `vapor_stream` over `dt_s` in `[0, 1]`.
    pub vent_headspace_fraction: f64,
    /// Optional pressure-cooker relief threshold [Pa]: if `Some(p_rel)`, vapor vents only
    /// when `p_vessel_pa >= p_rel`.
    pub relief_pressure_pa: Option<f64>,
    /// If `true`, liquid displaced above `H_pot` during Bikerman boil-over exits through
    /// `liquid_stream` using exact donor-rule species & enthalpy accounting.
    pub allow_boilover_overflow: bool,
}

impl OpenVentSpecV1 {
    /// Fully open saucepan/stockpot/skillet (all generated headspace vapor vents each step).
    #[must_use]
    pub const fn open_pot() -> Self {
        Self {
            vent_headspace_fraction: 1.0,
            relief_pressure_pa: None,
            allow_boilover_overflow: true,
        }
    }

    /// Open commercial deep fryer (all food steam vents cleanly; no liquid overflow).
    #[must_use]
    pub const fn open_fryer() -> Self {
        Self {
            vent_headspace_fraction: 1.0,
            relief_pressure_pa: None,
            allow_boilover_overflow: false,
        }
    }

    /// Sealed vessel (zero automatic venting).
    #[must_use]
    pub const fn sealed() -> Self {
        Self {
            vent_headspace_fraction: 0.0,
            relief_pressure_pa: None,
            allow_boilover_overflow: false,
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct PortIndices {
    wall_heat: usize,
    immersed_food_heat: usize,
    liquid_stream: usize,
    vapor_stream: usize,
    stir_rpm: usize,
    p_vessel_pa: usize,
    probes: usize,
}

/// Stirred multiphase liquid/vapor/oil/emulsion vessel implementing [`Component`].
#[derive(Clone, Debug)]
pub struct MultiphaseVesselComponent {
    config: MultiphaseVesselConfig,
    ports: Vec<PortSpec>,
    layout: PortLayout,
    indices: PortIndices,
    reaction_elemental: Option<ElementalMatrix>,
    reaction_molar_masses: Vec<(SpeciesId, MolarMass)>,
}

impl MultiphaseVesselComponent {
    /// Construct a validated [`MultiphaseVesselComponent`].
    pub fn new(config: MultiphaseVesselConfig) -> Result<Self, VesselErrorV1> {
        config.validate()?;
        let (ports, indices) = build_vessel_ports(&config)?;
        let layout = PortLayout::new(&ports)?;
        let (reaction_elemental, reaction_molar_masses) = if config.hydrolysis_enabled {
            let (elem, masses) = build_hydrolysis_reaction_basis()?;
            (Some(elem), masses)
        } else {
            (None, Vec::new())
        };
        Ok(Self {
            config,
            ports,
            layout,
            indices,
            reaction_elemental,
            reaction_molar_masses,
        })
    }

    /// Borrow the vessel configuration.
    #[must_use]
    pub const fn config(&self) -> &MultiphaseVesselConfig {
        &self.config
    }

    /// Borrow the validated [`PortLayout`].
    #[must_use]
    pub const fn layout(&self) -> &PortLayout {
        &self.layout
    }

    /// Port index of `wall_heat`.
    #[must_use]
    pub const fn wall_heat_port(&self) -> usize {
        self.indices.wall_heat
    }

    /// Port index of `immersed_food_heat`.
    #[must_use]
    pub const fn immersed_food_heat_port(&self) -> usize {
        self.indices.immersed_food_heat
    }

    /// Port index of `liquid_stream`.
    #[must_use]
    pub const fn liquid_stream_port(&self) -> usize {
        self.indices.liquid_stream
    }

    /// Port index of `vapor_stream`.
    #[must_use]
    pub const fn vapor_stream_port(&self) -> usize {
        self.indices.vapor_stream
    }

    /// Port index of `stir_rpm`.
    #[must_use]
    pub const fn stir_rpm_port(&self) -> usize {
        self.indices.stir_rpm
    }

    /// Port index of `p_vessel_pa`.
    #[must_use]
    pub const fn p_vessel_pa_port(&self) -> usize {
        self.indices.p_vessel_pa
    }

    /// Port index of `probes`.
    #[must_use]
    pub const fn probes_port(&self) -> usize {
        self.indices.probes
    }

    /// Build a thermodynamically consistent initial [`MultiphaseVesselState`] on the
    /// `janaf-298.15K-formation` enthalpy datum.
    pub fn initial_state(
        &self,
        spec: VesselInitialStateV1,
    ) -> Result<MultiphaseVesselState, VesselErrorV1> {
        let t_cz = if self.config.cold_zone_enabled {
            spec.t_coldzone_k
        } else {
            spec.t_liq_k
        };
        let m_cz = if self.config.cold_zone_enabled {
            spec.m_oil_coldzone_kg
        } else {
            0.0
        };
        let mut state = MultiphaseVesselState {
            time_s: 0.0,
            t_liq_k: spec.t_liq_k,
            t_coldzone_k: t_cz,
            active_enthalpy_j: 0.0,
            coldzone_enthalpy_j: 0.0,
            m_h2o_liq_kg: spec.m_h2o_liq_kg,
            m_etoh_liq_kg: spec.m_etoh_liq_kg,
            m_triglyceride_kg: spec.m_triglyceride_kg,
            m_ffa_kg: spec.m_ffa_kg,
            m_glycerol_kg: spec.m_glycerol_kg,
            m_salt_kg: spec.m_salt_kg,
            m_sugar_kg: spec.m_sugar_kg,
            m_fond_dissolved_kg: spec.m_fond_dissolved_kg,
            m_fond_adhered_kg: spec.m_fond_adhered_kg,
            m_milk_solids_kg: spec.m_milk_solids_kg,
            m_starch_kg: spec.m_starch_kg,
            m_protein_kg: spec.m_protein_kg,
            m_oil_coldzone_kg: m_cz,
            m_h2o_vap_kg: spec.m_h2o_vap_kg,
            m_etoh_vap_kg: spec.m_etoh_vap_kg,
            aux_gas_kg: spec.aux_gas_kg,
            tpm_fraction: spec.tpm_fraction,
            foam_height_m: 0.0,
            milk_solid_browning: 0.0,
            cumulative_hydrolysis_mol: 0.0,
        };
        state.active_enthalpy_j = compute_active_enthalpy_j(&state, &self.config, state.t_liq_k);
        state.coldzone_enthalpy_j =
            compute_coldzone_enthalpy_j(&state, &self.config, state.t_coldzone_k);
        state.validate(&self.config)?;
        Ok(state)
    }

    /// Compute effective liquid density `rho_liq` [kg/m^3] of the active liquid phase.
    #[must_use]
    pub fn liquid_density_kg_m3(&self, state: &MultiphaseVesselState) -> f64 {
        let m_liq = state.active_liquid_mass_kg();
        if m_liq <= 1.0e-15 {
            return if let Some(oil) = &self.config.oil {
                oil.density_kg_m3(state.t_liq_k)
            } else {
                958.0
            };
        }
        let rho_w = (1_000.0 - 0.35 * (state.t_liq_k - 277.15).max(0.0)).clamp(700.0, 1_005.0);
        let rho_e = (789.0 - 0.85 * (state.t_liq_k - 293.15).max(0.0)).clamp(600.0, 810.0);
        let rho_oil = self.config.oil.as_ref().map_or_else(
            || (914.0 - 0.66 * (state.t_liq_k - JANAF_REFERENCE_TEMPERATURE_K)).clamp(700.0, 960.0),
            |card| card.density_kg_m3(state.t_liq_k),
        );
        let rho_sol = 1_550.0;
        let vol = state.m_h2o_liq_kg / rho_w
            + state.m_etoh_liq_kg / rho_e
            + (state.m_triglyceride_kg + state.m_ffa_kg + state.m_glycerol_kg) / rho_oil
            + (state.m_salt_kg
                + state.m_sugar_kg
                + state.m_fond_dissolved_kg
                + state.m_milk_solids_kg
                + state.m_starch_kg
                + state.m_protein_kg)
                / rho_sol;
        if vol > 1.0e-18 {
            (m_liq / vol).clamp(500.0, 2_000.0)
        } else {
            958.0
        }
    }

    /// Smooth dryout wetting fraction `f_wet = tanh(m_liq / (rho_l * A_base * delta_min))` in `[0, 1]`.
    #[must_use]
    pub fn wetting_fraction(&self, state: &MultiphaseVesselState) -> f64 {
        let m_liq = state.active_liquid_mass_kg();
        if m_liq <= 0.0 {
            return 0.0;
        }
        let rho_l = self.liquid_density_kg_m3(state);
        let m_ref = (rho_l * self.config.base_area_m2 * self.config.dryout_delta_min_m).max(1.0e-9);
        tanh(m_liq / m_ref)
    }

    /// Evaluate volatile bubble point `T_bub(P, x_etoh, solutes)` [K], exact NRTL equilibrium
    /// vapor ethanol mole fraction `y_etoh`, and solute boiling-point elevation `Delta T_BPE` [K]
    /// at caller-supplied headspace pressure `p_vessel_pa` [Pa].
    pub fn bubble_point_and_vapor_composition(
        &self,
        state: &MultiphaseVesselState,
        p_vessel_pa: f64,
    ) -> Result<(f64, f64, f64), VesselErrorV1> {
        if !p_vessel_pa.is_finite() || p_vessel_pa <= 0.0 {
            return Err(VesselErrorV1::InvalidPressure {
                pressure_pa: p_vessel_pa,
            });
        }
        let pressure = Pressure::new(p_vessel_pa);
        let x_etoh = state.ethanol_mole_fraction();
        let (t_bub_base, y_etoh) = if x_etoh > 1.0e-9 {
            let vle = NrtlModelV1::ethanol_water_vle_v1();
            let bp = vle
                .bubble_point_t(
                    pressure,
                    &[x_etoh, (1.0 - x_etoh).max(0.0)],
                    VleSolverOptionsV1::default(),
                )
                .map_err(|_| VesselErrorV1::Thermochem("NRTL ethanol-water bubble point failed"))?;
            (bp.temperature.value(), bp.vapor_y[0].clamp(0.0, 1.0))
        } else {
            let sat = If97SaturationPair::at_pressure(pressure).map_err(|_| {
                VesselErrorV1::InvalidPressure {
                    pressure_pa: p_vessel_pa,
                }
            })?;
            (sat.temperature().value(), 0.0)
        };

        // Exact Raoult / Clausius-Clapeyron boiling-point elevation from dissolved non-volatile solutes.
        let dt_bpe = if state.m_h2o_liq_kg > 1.0e-12 {
            let n_w = state.m_h2o_liq_kg / MOLAR_MASS_H2O_KG_MOL;
            let n_eff_solute = 2.0 * (state.m_salt_kg / MOLAR_MASS_NACL_KG_MOL)
                + (state.m_sugar_kg / MOLAR_MASS_SUCROSE_KG_MOL)
                + (state.m_fond_dissolved_kg / MOLAR_MASS_FOND_EFF_KG_MOL);
            if n_eff_solute > 0.0 {
                let x_w_eff = (n_w / (n_w + n_eff_solute)).clamp(0.25, 1.0);
                let h_fg_mol = (specific_enthalpy_h2o_vap(t_bub_base)
                    - specific_enthalpy_h2o_liq(t_bub_base))
                    * MOLAR_MASS_H2O_KG_MOL;
                let raw_bpe = -(UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * t_bub_base * t_bub_base
                    / h_fg_mol.max(10_000.0))
                    * ln(x_w_eff);
                raw_bpe.clamp(0.0, 45.0)
            } else {
                0.0
            }
        } else {
            0.0
        };
        Ok((t_bub_base + dt_bpe, y_etoh, dt_bpe))
    }

    /// Evaluate effective liquid dynamic viscosity [Pa s].
    ///
    /// - Lipid-dominated baths (`m_lipid >= 0.5 * m_active_liq`): oil card Arrhenius-TPM law
    ///   `mu(T, TPM) = mu_0 * exp(E_visc / (R T)) * exp(gamma_TPM * TPM)`.
    /// - Aqueous/sauce/egg baths: Vogel-Fulcher-Tammann base water viscosity times concentration
    ///   thickening with Carreau shear-thinning driven by `stir_rpm`.
    #[must_use]
    pub fn effective_viscosity_pa_s(&self, state: &MultiphaseVesselState, stir_rpm: f64) -> f64 {
        let rpm = stir_rpm.max(0.0);
        let m_liq = state.active_liquid_mass_kg();
        if m_liq <= 1.0e-15 {
            return 0.0;
        }
        let m_lip = state.active_lipid_mass_kg();
        if m_lip >= 0.5 * m_liq {
            if let Some(oil) = &self.config.oil {
                return oil.viscosity_pa_s(state.t_liq_k, state.tpm_fraction);
            }
            return OilCardV1::canola_refined().viscosity_pa_s(state.t_liq_k, state.tpm_fraction);
        }

        // Aqueous / sauce Carreau + Herschel-Bulkley shear-thinning viscosity.
        let t_clamped = state.t_liq_k.clamp(274.0, 500.0);
        let mu_w = 2.414e-5 * pow(10.0, 247.8 / (t_clamped - 140.0));
        let x_sol = state.solids_mass_fraction();
        let denom = (1.0 - x_sol).max(0.08);
        let exponent = (self.config.sauce_viscosity_conc_coeff * x_sol / denom).min(10.0);
        let mu_zero = mu_w * exp(exponent);
        let mu_inf = mu_w * (1.0 + 2.5 * x_sol);

        // Metzner-Otto effective shear rate gamma_dot = max(1.0, 11.5 * N_rps).
        let gamma_dot = (11.5 * (rpm / 60.0)).max(1.0);
        let lam_g = self.config.sauce_carreau_lambda_s * gamma_dot;
        let shear_factor = pow(
            1.0 + lam_g * lam_g,
            0.5 * (self.config.sauce_carreau_n - 1.0),
        );
        let mu_carreau = mu_inf + (mu_zero - mu_inf) * shear_factor;

        // Herschel-Bulkley yield contribution for high-solids thickened reductions.
        let tau_yield = if x_sol > 0.12 {
            18.0 * pow((x_sol - 0.12).min(0.60), 2.0)
        } else {
            0.0
        };
        mu_carreau + tau_yield / (gamma_dot + 0.25)
    }

    /// Evaluate degraded oil smoke point [K] (`0.0` if no oil/butter card is configured).
    #[must_use]
    pub fn smoke_point_k(&self, state: &MultiphaseVesselState) -> f64 {
        self.config.oil.as_ref().map_or(0.0, |oil| {
            oil.smoke_point_k(state.ffa_fraction(), state.tpm_fraction)
        })
    }

    /// Effective wall-to-fluid heat-transfer coefficient [W/(m^2 K)] incorporating
    /// stirring shear enhancement (3–5x for scrambled eggs/sauces), vapor-bubble agitation,
    /// and viscous fouling damping.
    #[must_use]
    pub fn effective_wall_htc_w_m2_k(
        &self,
        state: &MultiphaseVesselState,
        stir_rpm: f64,
        j_g_m_s: f64,
    ) -> f64 {
        let rpm = stir_rpm.max(0.0);
        let jg = j_g_m_s.max(0.0);
        let n_rps = rpm / 60.0;
        let stir_mult = 1.0 + self.config.stir_htc_gain * pow(n_rps, 0.67);
        let bubble_mult =
            1.0 + self.config.bubble_htc_gain * tanh(jg / self.config.bubble_jg_ref_m_s);

        let mu = self.effective_viscosity_pa_s(state, rpm);
        let mu_ref = if state.active_lipid_mass_kg() >= 0.5 * state.active_liquid_mass_kg() {
            0.0020 // hot frying oil reference ~2 mPa*s
        } else {
            0.0005 // hot water reference ~0.5 mPa*s
        };
        let visc_factor = if mu > 0.0 {
            pow((mu_ref / mu).clamp(0.01, 2.0), 0.30)
        } else {
            1.0
        };
        let f_wet = self.wetting_fraction(state);
        self.config.natural_wall_htc_w_m2_k
            * stir_mult
            * bubble_mult
            * visc_factor
            * f_wet.max(0.02)
    }

    /// Effective fluid-to-immersed-food heat-transfer coefficient [W/(m^2 K)]
    /// incorporating stirring and fryer vapor-bubble agitation.
    #[must_use]
    pub fn effective_immersed_htc_w_m2_k(
        &self,
        state: &MultiphaseVesselState,
        stir_rpm: f64,
        j_g_m_s: f64,
    ) -> f64 {
        let rpm = stir_rpm.max(0.0);
        let jg = j_g_m_s.max(0.0);
        let n_rps = rpm / 60.0;
        let stir_mult = 1.0 + self.config.stir_htc_gain * pow(n_rps, 0.67);
        let bubble_mult =
            1.0 + self.config.bubble_htc_gain * tanh(jg / self.config.bubble_jg_ref_m_s);
        let mu = self.effective_viscosity_pa_s(state, rpm);
        let mu_ref = if state.active_lipid_mass_kg() >= 0.5 * state.active_liquid_mass_kg() {
            0.0020
        } else {
            0.0005
        };
        let visc_factor = if mu > 0.0 {
            pow((mu_ref / mu).clamp(0.02, 2.0), 0.25)
        } else {
            1.0
        };
        self.config.natural_immersed_htc_w_m2_k * stir_mult * bubble_mult * visc_factor
    }

    /// Estimate wall-contact liquid film temperature [K] under wall heat load `q_wall_w`.
    #[must_use]
    pub fn estimated_wall_film_temp_k(
        &self,
        state: &MultiphaseVesselState,
        q_wall_w: f64,
        stir_rpm: f64,
        j_g_m_s: f64,
    ) -> f64 {
        let h_wall = self
            .effective_wall_htc_w_m2_k(state, stir_rpm, j_g_m_s)
            .max(10.0);
        let dt_film = (q_wall_w / (h_wall * self.config.base_area_m2)).clamp(-100.0, 350.0);
        (state.t_liq_k + dt_film).max(TEMP_MIN_VALID_K)
    }

    /// Liquid column height `H_liq` [m] in the vessel from active+coldzone liquid volume.
    #[must_use]
    pub fn liquid_height_m(&self, state: &MultiphaseVesselState) -> f64 {
        let rho = self.liquid_density_kg_m3(state);
        state.total_liquid_mass_kg() / (rho * self.config.free_surface_area_m2)
    }

    /// Evaluate structured boil-over, oil-safety, and scorching flags.
    #[must_use]
    pub fn safety_flags(
        &self,
        state: &MultiphaseVesselState,
        q_wall_w: f64,
        stir_rpm: f64,
    ) -> VesselSafetyFlagsV1 {
        let h_liq = self.liquid_height_m(state);
        let boil_over = state.foam_height_m > 1.0e-5
            && (h_liq + state.foam_height_m) > self.config.vessel_height_m;
        let (smoking, flashpoint_exceeded, autoignition_exceeded, oil_discard_exceeded) =
            if let Some(oil) = &self.config.oil {
                let t_smoke = oil.smoke_point_k(state.ffa_fraction(), state.tpm_fraction);
                (
                    state.t_liq_k >= t_smoke,
                    state.t_liq_k >= oil.t_flash_k,
                    state.t_liq_k >= oil.t_autoignition_k,
                    state.tpm_fraction >= oil.tpm_discard_limit
                        || state.ffa_fraction() >= oil.ffa_discard_limit,
                )
            } else {
                (false, false, false, false)
            };

        let is_aqueous_sauce = state.active_liquid_mass_kg() > 1.0e-6
            && state.active_lipid_mass_kg() < 0.5 * state.active_liquid_mass_kg()
            && state.solids_mass_fraction() >= 0.05;
        let scorch_risk = if is_aqueous_sauce && q_wall_w > 0.0 {
            self.estimated_wall_film_temp_k(state, q_wall_w, stir_rpm, 0.0)
                >= self.config.scorch_wall_threshold_k
        } else {
            false
        };

        VesselSafetyFlagsV1 {
            boil_over,
            smoking,
            flashpoint_exceeded,
            autoignition_exceeded,
            oil_discard_exceeded,
            scorch_risk,
        }
    }

    /// Specific enthalpy [J/kg] of the active liquid mixture on `janaf-298.15K-formation`.
    #[must_use]
    pub fn specific_enthalpy_janaf_liquid(&self, state: &MultiphaseVesselState) -> f64 {
        let m_liq = state.active_liquid_mass_kg();
        let t = state.t_liq_k;
        if m_liq <= 1.0e-15 {
            return specific_enthalpy_h2o_liq(t);
        }
        let h_liq = state.m_h2o_liq_kg * specific_enthalpy_h2o_liq(t)
            + state.m_etoh_liq_kg * specific_enthalpy_etoh_liq(t)
            + state.m_triglyceride_kg * specific_enthalpy_trig(&self.config, t)
            + state.m_ffa_kg * specific_enthalpy_ffa(t)
            + state.m_glycerol_kg * specific_enthalpy_glycerol(t)
            + state.m_salt_kg * specific_enthalpy_salt(t)
            + state.m_sugar_kg * specific_enthalpy_sugar(t)
            + state.m_fond_dissolved_kg * specific_enthalpy_fond(t)
            + state.m_milk_solids_kg * specific_enthalpy_milk_solids(t)
            + state.m_starch_kg * specific_enthalpy_starch(t)
            + state.m_protein_kg * specific_enthalpy_protein(t);
        h_liq / m_liq
    }

    /// Specific enthalpy [J/kg] of the headspace vapor mixture on `janaf-298.15K-formation`.
    #[must_use]
    pub fn specific_enthalpy_janaf_vapor(&self, state: &MultiphaseVesselState) -> f64 {
        let m_vap = state.total_vapor_mass_kg();
        let t = state.t_liq_k;
        if m_vap <= 1.0e-15 {
            return specific_enthalpy_h2o_vap(t);
        }
        let mut h_vap = state.m_h2o_vap_kg * specific_enthalpy_h2o_vap(t)
            + state.m_etoh_vap_kg * specific_enthalpy_etoh_vap(t);
        for k in 0..AUX_GAS_COUNT {
            h_vap += state.aux_gas_kg[k] * specific_enthalpy_aux_gas(k, t);
        }
        h_vap / m_vap
    }

    /// Construct a valid `kitchen_liquid_v1` input wire slice `[m_dot, h_dot, species[8]]`
    /// carrying a single liquid species at `temperature_k` on `janaf-298.15K-formation`.
    ///
    /// `species_index`: `0=H2O, 1=C2H5OH, 2=triglyceride, 3=ffa, 4=glycerol, 5=solute_salt, 6=solute_sugar, 7=solute_fond`.
    #[must_use]
    pub fn liquid_stream_for_pure_species(
        species_index: usize,
        mass_rate_kg_s: f64,
        temperature_k: f64,
        oil: Option<&OilCardV1>,
    ) -> [f64; LIQUID_STREAM_WIDTH] {
        let mut wire = [0.0; LIQUID_STREAM_WIDTH];
        if species_index >= 8 || mass_rate_kg_s == 0.0 {
            return wire;
        }
        let h_spec = match species_index {
            0 => specific_enthalpy_h2o_liq(temperature_k),
            1 => specific_enthalpy_etoh_liq(temperature_k),
            2 => {
                let (cp_ref, dcp_dt, t_ref) = oil.map_or(
                    (
                        CP_TRIG_DEFAULT_J_KG_K,
                        DCP_DT_TRIG_DEFAULT_J_KG_K2,
                        JANAF_REFERENCE_TEMPERATURE_K,
                    ),
                    |c| (c.cp_ref_j_kg_k, c.dcp_dt_j_kg_k2, c.t_ref_k),
                );
                let u_ref = JANAF_REFERENCE_TEMPERATURE_K - t_ref;
                let u = temperature_k - t_ref;
                H_F0_TRIGLYCERIDE_J_KG
                    + cp_ref * (u - u_ref)
                    + 0.5 * dcp_dt * (u * u - u_ref * u_ref)
            }
            3 => specific_enthalpy_ffa(temperature_k),
            4 => specific_enthalpy_glycerol(temperature_k),
            5 => specific_enthalpy_salt(temperature_k),
            6 => specific_enthalpy_sugar(temperature_k),
            7 => specific_enthalpy_fond(temperature_k),
            _ => 0.0,
        };
        wire[STREAM_MASS_RATE] = mass_rate_kg_s;
        wire[STREAM_ENTHALPY_RATE] = mass_rate_kg_s * h_spec;
        wire[STREAM_CARRIED + species_index] = mass_rate_kg_s;
        wire
    }

    /// Construct a valid `kitchen_gas_v1` input wire slice `[m_dot, h_dot, species[9]]`
    /// carrying pure water vapor (steam) at `temperature_k` on `janaf-298.15K-formation`.
    #[must_use]
    pub fn vapor_stream_for_steam(
        mass_rate_kg_s: f64,
        temperature_k: f64,
    ) -> [f64; VAPOR_STREAM_WIDTH] {
        let mut wire = [0.0; VAPOR_STREAM_WIDTH];
        wire[STREAM_MASS_RATE] = mass_rate_kg_s;
        wire[STREAM_ENTHALPY_RATE] = mass_rate_kg_s * specific_enthalpy_h2o_vap(temperature_k);
        wire[STREAM_CARRIED] = mass_rate_kg_s;
        wire
    }

    /// Assemble a complete `PortLayout`-ordered input vector for [`Component::propose`].
    pub fn build_inputs(
        &self,
        q_wall_w: f64,
        q_food_w: &[f64],
        liq_stream_in: [f64; LIQUID_STREAM_WIDTH],
        vap_stream_in: [f64; VAPOR_STREAM_WIDTH],
        stir_rpm: f64,
        p_vessel_pa: f64,
    ) -> Result<Vec<f64>, VesselErrorV1> {
        if q_food_w.len() != self.config.immersed_batch_count {
            return Err(VesselErrorV1::InvalidConfig {
                parameter: "q_food_w",
                reason: "length must match immersed_batch_count",
            });
        }
        let mut inputs = vec![0.0; self.layout.input_len()];
        inputs[self.layout.input_range(self.indices.wall_heat)][0] = q_wall_w;
        inputs[self.layout.input_range(self.indices.immersed_food_heat)].copy_from_slice(q_food_w);
        inputs[self.layout.input_range(self.indices.liquid_stream)].copy_from_slice(&liq_stream_in);
        inputs[self.layout.input_range(self.indices.vapor_stream)].copy_from_slice(&vap_stream_in);
        inputs[self.layout.input_range(self.indices.stir_rpm)][0] = stir_rpm;
        inputs[self.layout.input_range(self.indices.p_vessel_pa)][0] = p_vessel_pa;
        Ok(inputs)
    }

    /// Build donor-rule open-vent (and optional boil-over overflow) `PortLayout` inputs
    /// that remove accumulated headspace vapor through `vapor_stream` with exact
    /// mass, species, and `janaf-298.15K-formation` enthalpy conservation.
    #[allow(clippy::too_many_arguments)]
    pub fn open_vent_inputs(
        &self,
        state: &MultiphaseVesselState,
        dt_s: f64,
        q_wall_w: f64,
        q_food_w: &[f64],
        extra_liq_in: [f64; LIQUID_STREAM_WIDTH],
        extra_vap_in: [f64; VAPOR_STREAM_WIDTH],
        stir_rpm: f64,
        p_vessel_pa: f64,
        vent: OpenVentSpecV1,
    ) -> Result<Vec<f64>, VesselErrorV1> {
        if !dt_s.is_finite() || dt_s <= 0.0 {
            return Err(VesselErrorV1::InvalidScalar {
                quantity: "dt_s",
                value: dt_s,
                reason: "must be finite and strictly positive",
            });
        }
        if !vent.vent_headspace_fraction.is_finite()
            || !(0.0..=1.0).contains(&vent.vent_headspace_fraction)
        {
            return Err(VesselErrorV1::InvalidScalar {
                quantity: "vent_headspace_fraction",
                value: vent.vent_headspace_fraction,
                reason: "must lie in [0, 1]",
            });
        }

        let mut vap_in = extra_vap_in;
        let relief_open = vent
            .relief_pressure_pa
            .is_none_or(|p_rel| p_vessel_pa >= p_rel);
        if relief_open && vent.vent_headspace_fraction > 0.0 {
            let f = vent.vent_headspace_fraction;
            let t = state.t_liq_k;
            let out_h2o = -f * state.m_h2o_vap_kg / dt_s;
            let out_etoh = -f * state.m_etoh_vap_kg / dt_s;
            let mut sum_m = out_h2o + out_etoh;
            let mut sum_h =
                out_h2o * specific_enthalpy_h2o_vap(t) + out_etoh * specific_enthalpy_etoh_vap(t);
            vap_in[STREAM_CARRIED] += out_h2o;
            vap_in[STREAM_CARRIED + 1] += out_etoh;
            for k in 0..AUX_GAS_COUNT {
                let out_k = -f * state.aux_gas_kg[k] / dt_s;
                vap_in[STREAM_CARRIED + 2 + k] += out_k;
                sum_m += out_k;
                sum_h += out_k * specific_enthalpy_aux_gas(k, t);
            }
            vap_in[STREAM_MASS_RATE] += sum_m;
            vap_in[STREAM_ENTHALPY_RATE] += sum_h;
        }

        let mut liq_in = extra_liq_in;
        if vent.allow_boilover_overflow {
            let h_liq = self.liquid_height_m(state);
            let excess_h = h_liq + state.foam_height_m - self.config.vessel_height_m;
            if state.foam_height_m > 1.0e-5 && excess_h > 0.0 {
                // Only liquid species represented in kitchen_liquid_v1 leave via overflow.
                let m_streamable = state.m_h2o_liq_kg
                    + state.m_etoh_liq_kg
                    + state.m_triglyceride_kg
                    + state.m_ffa_kg
                    + state.m_glycerol_kg
                    + state.m_salt_kg
                    + state.m_sugar_kg
                    + state.m_fond_dissolved_kg;
                if m_streamable > 1.0e-9 {
                    let rho = self.liquid_density_kg_m3(state);
                    // Foam is ~5% liquid volume fraction; cap overflow at 15% of liquid per step.
                    let m_overflow = (0.05 * rho * self.config.free_surface_area_m2 * excess_h)
                        .min(0.15 * m_streamable);
                    let frac = m_overflow / m_streamable;
                    let t = state.t_liq_k;
                    let species = [
                        (state.m_h2o_liq_kg, specific_enthalpy_h2o_liq(t)),
                        (state.m_etoh_liq_kg, specific_enthalpy_etoh_liq(t)),
                        (
                            state.m_triglyceride_kg,
                            specific_enthalpy_trig(&self.config, t),
                        ),
                        (state.m_ffa_kg, specific_enthalpy_ffa(t)),
                        (state.m_glycerol_kg, specific_enthalpy_glycerol(t)),
                        (state.m_salt_kg, specific_enthalpy_salt(t)),
                        (state.m_sugar_kg, specific_enthalpy_sugar(t)),
                        (state.m_fond_dissolved_kg, specific_enthalpy_fond(t)),
                    ];
                    let mut sum_m = 0.0;
                    let mut sum_h = 0.0;
                    for (idx, (m_s, h_s)) in species.into_iter().enumerate() {
                        let out_s = -frac * m_s / dt_s;
                        liq_in[STREAM_CARRIED + idx] += out_s;
                        sum_m += out_s;
                        sum_h += out_s * h_s;
                    }
                    liq_in[STREAM_MASS_RATE] += sum_m;
                    liq_in[STREAM_ENTHALPY_RATE] += sum_h;
                }
            }
        }

        self.build_inputs(q_wall_w, q_food_w, liq_in, vap_in, stir_rpm, p_vessel_pa)
    }

    #[allow(clippy::too_many_arguments)]
    fn write_port_outputs(
        &self,
        old: &MultiphaseVesselState,
        new: &MultiphaseVesselState,
        q_wall_w: f64,
        stir_rpm: f64,
        p_vessel_pa: f64,
        outputs: &mut [f64],
    ) {
        outputs.fill(0.0);
        let t_avg = 0.5 * (old.t_liq_k + new.t_liq_k);

        // 0. wall_heat output: step-average fluid temperature at wall [K].
        outputs[self.layout.output_range(self.indices.wall_heat)][0] = t_avg;

        // 1. immersed_food_heat output: step-average bath temperature [K] per batch.
        let food_out = &mut outputs[self.layout.output_range(self.indices.immersed_food_heat)];
        food_out.fill(t_avg);

        // 2. liquid_stream output: [P_liq, h_liq_janaf, y_H2O, y_C2H5OH, y_trig, y_ffa, y_glyc, y_salt, y_sugar, y_fond].
        let liq_out = &mut outputs[self.layout.output_range(self.indices.liquid_stream)];
        let rho_l = self.liquid_density_kg_m3(new);
        let h_liq_m = self.liquid_height_m(new);
        let p_liq = p_vessel_pa + 0.5 * rho_l * 9.806_65 * h_liq_m;
        let m_stream_liq = new.m_h2o_liq_kg
            + new.m_etoh_liq_kg
            + new.m_triglyceride_kg
            + new.m_ffa_kg
            + new.m_glycerol_kg
            + new.m_salt_kg
            + new.m_sugar_kg
            + new.m_fond_dissolved_kg;
        liq_out[STREAM_PRESSURE] = p_liq;
        liq_out[STREAM_SPECIFIC_ENTHALPY] = self.specific_enthalpy_janaf_liquid(new);
        if m_stream_liq > 1.0e-18 {
            liq_out[STREAM_CARRIED] = new.m_h2o_liq_kg / m_stream_liq;
            liq_out[STREAM_CARRIED + 1] = new.m_etoh_liq_kg / m_stream_liq;
            liq_out[STREAM_CARRIED + 2] = new.m_triglyceride_kg / m_stream_liq;
            liq_out[STREAM_CARRIED + 3] = new.m_ffa_kg / m_stream_liq;
            liq_out[STREAM_CARRIED + 4] = new.m_glycerol_kg / m_stream_liq;
            liq_out[STREAM_CARRIED + 5] = new.m_salt_kg / m_stream_liq;
            liq_out[STREAM_CARRIED + 6] = new.m_sugar_kg / m_stream_liq;
            liq_out[STREAM_CARRIED + 7] = new.m_fond_dissolved_kg / m_stream_liq;
        } else {
            liq_out[STREAM_CARRIED] = 1.0;
        }

        // 3. vapor_stream output: [P_vap, h_vap_janaf, y_H2O, y_C2H5OH, y_CO2, y_O2, y_N2, y_CH4, y_C3H8, y_CO, y_grease].
        let vap_out = &mut outputs[self.layout.output_range(self.indices.vapor_stream)];
        let m_vap = new.total_vapor_mass_kg();
        vap_out[STREAM_PRESSURE] = p_vessel_pa;
        vap_out[STREAM_SPECIFIC_ENTHALPY] = self.specific_enthalpy_janaf_vapor(new);
        if m_vap > 1.0e-18 {
            vap_out[STREAM_CARRIED] = new.m_h2o_vap_kg / m_vap;
            vap_out[STREAM_CARRIED + 1] = new.m_etoh_vap_kg / m_vap;
            for k in 0..AUX_GAS_COUNT {
                vap_out[STREAM_CARRIED + 2 + k] = new.aux_gas_kg[k] / m_vap;
            }
        } else {
            vap_out[STREAM_CARRIED] = 1.0;
        }

        // 6. probes output (12-wide SignalOut).
        let flags = self.safety_flags(new, q_wall_w, stir_rpm);
        let probes = &mut outputs[self.layout.output_range(self.indices.probes)];
        probes[PROBE_T_LIQ_K] = new.t_liq_k;
        probes[PROBE_T_COLDZONE_K] = new.t_coldzone_k;
        probes[PROBE_M_LIQ_KG] = new.total_liquid_mass_kg();
        probes[PROBE_X_ETOH] = new.ethanol_mole_fraction();
        probes[PROBE_X_SOLIDS] = new.solids_mass_fraction();
        probes[PROBE_VISC_PA_S] = self.effective_viscosity_pa_s(new, stir_rpm);
        probes[PROBE_FFA_PCT] = 100.0 * new.ffa_fraction();
        probes[PROBE_TPM_PCT] = 100.0 * new.tpm_fraction;
        probes[PROBE_T_SMOKE_K] = self.smoke_point_k(new);
        probes[PROBE_FOAM_H_M] = new.foam_height_m;
        probes[PROBE_FOND_ADHERED_KG] = new.m_fond_adhered_kg;
        probes[PROBE_SAFETY_OR_BOILOVER_FLAG] = flags.as_flag_f64();
    }
}

impl Component for MultiphaseVesselComponent {
    type State = MultiphaseVesselState;
    type Error = VesselErrorV1;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Smooth
    }

    fn stored_energy_j(&self, state: &Self::State) -> Result<f64, Self::Error> {
        state.validate(&self.config)?;
        Ok(state.total_enthalpy_j())
    }

    fn stored_mass_kg(&self, state: &Self::State) -> Result<Option<f64>, Self::Error> {
        state.validate(&self.config)?;
        Ok(Some(state.total_mass_kg()))
    }

    fn stored_carried_quantity(
        &self,
        state: &Self::State,
        carried: usize,
    ) -> Result<Option<f64>, Self::Error> {
        state.validate(&self.config)?;
        let Some(StreamQuantity::SpeciesMass(species)) =
            self.layout.carried_quantities().get(carried)
        else {
            return Ok(None);
        };
        let amount = match species.as_str() {
            "H2O" => state.m_h2o_liq_kg + state.m_h2o_vap_kg,
            "C2H5OH" => state.m_etoh_liq_kg + state.m_etoh_vap_kg,
            "triglyceride" => state.m_triglyceride_kg + state.m_oil_coldzone_kg,
            "ffa" => state.m_ffa_kg,
            "glycerol" => state.m_glycerol_kg,
            "solute_salt" => state.m_salt_kg,
            "solute_sugar" => state.m_sugar_kg,
            // Adhered pan fond is owned inside MultiphaseVesselState so dissolution closes identically.
            "solute_fond" => state.m_fond_dissolved_kg + state.m_fond_adhered_kg,
            "CO2" => state.aux_gas_kg[0],
            "O2" => state.aux_gas_kg[1],
            "N2" => state.aux_gas_kg[2],
            "CH4" => state.aux_gas_kg[3],
            "C3H8" => state.aux_gas_kg[4],
            "CO" => state.aux_gas_kg[5],
            "grease_vapor" => state.aux_gas_kg[6],
            _ => return Ok(None),
        };
        Ok(Some(amount))
    }

    fn carried_production(
        &self,
        old: &Self::State,
        new: &Self::State,
        dt_s: f64,
        carried: usize,
    ) -> Result<Option<f64>, Self::Error> {
        if !self.config.hydrolysis_enabled {
            return Ok(None);
        }
        if !dt_s.is_finite() || dt_s <= 0.0 {
            return Err(VesselErrorV1::InvalidScalar {
                quantity: "dt_s",
                value: dt_s,
                reason: "must be finite and strictly positive",
            });
        }
        let d_xi = (new.cumulative_hydrolysis_mol - old.cumulative_hydrolysis_mol).max(0.0);
        if d_xi == 0.0 {
            return Ok(Some(0.0));
        }
        let Some(StreamQuantity::SpeciesMass(species)) =
            self.layout.carried_quantities().get(carried)
        else {
            return Ok(Some(0.0));
        };
        // Exact stoichiometry: C57H104O6 (triglyceride) + 3 H2O -> 3 C18H34O2 (ffa) + C3H8O3 (glycerol).
        let rate = match species.as_str() {
            "triglyceride" => -(d_xi * MOLAR_MASS_TRIGLYCERIDE_KG_MOL) / dt_s,
            "H2O" => -(3.0 * d_xi * MOLAR_MASS_H2O_KG_MOL) / dt_s,
            "ffa" => (3.0 * d_xi * MOLAR_MASS_FFA_KG_MOL) / dt_s,
            "glycerol" => (d_xi * MOLAR_MASS_GLYCEROL_KG_MOL) / dt_s,
            _ => 0.0,
        };
        Ok(Some(rate))
    }

    fn reaction_basis(&self) -> Option<(&ElementalMatrix, &[(SpeciesId, MolarMass)])> {
        self.reaction_elemental
            .as_ref()
            .map(|elem| (elem, self.reaction_molar_masses.as_slice()))
    }

    fn limit_dt_s(&self, state: &Self::State, proposed_s: f64) -> Result<f64, Self::Error> {
        state.validate(&self.config)?;
        if !proposed_s.is_finite() || proposed_s <= 0.0 {
            return Err(VesselErrorV1::InvalidScalar {
                quantity: "proposed_s",
                value: proposed_s,
                reason: "must be finite and strictly positive",
            });
        }
        Ok(proposed_s)
    }

    #[allow(clippy::too_many_lines)]
    fn propose(
        &self,
        old: &Self::State,
        dt_s: f64,
        inputs: &[f64],
        new: &mut Self::State,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, Self::Error> {
        cx.checkpoint()?;
        old.validate(&self.config)?;
        if !dt_s.is_finite() || dt_s <= 0.0 {
            return Err(VesselErrorV1::InvalidScalar {
                quantity: "dt_s",
                value: dt_s,
                reason: "step must be finite and strictly positive",
            });
        }
        if inputs.len() != self.layout.input_len() || outputs.len() != self.layout.output_len() {
            return Err(VesselErrorV1::InvalidConfig {
                parameter: "port_vectors",
                reason: "inputs or outputs length does not match PortLayout",
            });
        }
        if !inputs.iter().all(|v| v.is_finite()) {
            return Err(VesselErrorV1::InvalidConfig {
                parameter: "inputs",
                reason: "all port input coefficients must be finite",
            });
        }

        let q_wall_w = inputs[self.layout.input_range(self.indices.wall_heat)][0];
        let q_food_slice = &inputs[self.layout.input_range(self.indices.immersed_food_heat)];
        let q_food_sum_w: f64 = q_food_slice.iter().sum();
        let liq_in = &inputs[self.layout.input_range(self.indices.liquid_stream)];
        let vap_in = &inputs[self.layout.input_range(self.indices.vapor_stream)];
        let stir_rpm = inputs[self.layout.input_range(self.indices.stir_rpm)][0];
        let p_vessel_pa = inputs[self.layout.input_range(self.indices.p_vessel_pa)][0];

        if stir_rpm < 0.0 {
            return Err(VesselErrorV1::InvalidScalar {
                quantity: "stir_rpm",
                value: stir_rpm,
                reason: "stir speed must be non-negative",
            });
        }
        if !p_vessel_pa.is_finite() || p_vessel_pa <= 0.0 {
            return Err(VesselErrorV1::InvalidPressure {
                pressure_pa: p_vessel_pa,
            });
        }

        // Validate stream mass-rate consistency against carried species sums.
        validate_stream_mass_consistency("liquid_stream", liq_in)?;
        validate_stream_mass_consistency("vapor_stream", vap_in)?;

        let h_dot_liq_w = liq_in[STREAM_ENTHALPY_RATE];
        let h_dot_vap_w = vap_in[STREAM_ENTHALPY_RATE];
        let liq_c = &liq_in[STREAM_CARRIED..];
        let vap_c = &vap_in[STREAM_CARRIED..];

        // 1. Apply boundary liquid & vapor stream transfers with clean dryout/depletion shortening.
        let mut shorten_dt: Option<f64> = None;
        let apply_delta = |m_old: f64,
                           rate_kg_s: f64,
                           species: &'static str,
                           shorten: &mut Option<f64>|
         -> Result<f64, VesselErrorV1> {
            let m_cand = m_old + dt_s * rate_kg_s;
            let tol = 1.0e-12 * m_old.max(1.0);
            if m_cand < -tol {
                if m_old > 1.0e-12 && rate_kg_s < 0.0 {
                    let dt_deplete = (0.5 * m_old / (-rate_kg_s)).max(self.config.min_dt_s);
                    if dt_deplete < dt_s {
                        *shorten = Some(shorten.map_or(dt_deplete, |prev| prev.min(dt_deplete)));
                        return Ok(0.0);
                    }
                }
                return Err(VesselErrorV1::DepletedInventory { species });
            }
            Ok(m_cand.max(0.0))
        };

        let mut m_h2o_liq = apply_delta(old.m_h2o_liq_kg, liq_c[0], "H2O_liq", &mut shorten_dt)?;
        let mut m_etoh_liq =
            apply_delta(old.m_etoh_liq_kg, liq_c[1], "C2H5OH_liq", &mut shorten_dt)?;
        let mut m_trig = apply_delta(
            old.m_triglyceride_kg,
            liq_c[2],
            "triglyceride",
            &mut shorten_dt,
        )?;
        let mut m_ffa = apply_delta(old.m_ffa_kg, liq_c[3], "ffa", &mut shorten_dt)?;
        let mut m_glyc = apply_delta(old.m_glycerol_kg, liq_c[4], "glycerol", &mut shorten_dt)?;
        let m_salt = apply_delta(old.m_salt_kg, liq_c[5], "solute_salt", &mut shorten_dt)?;
        let m_sugar = apply_delta(old.m_sugar_kg, liq_c[6], "solute_sugar", &mut shorten_dt)?;
        let mut m_fond_diss = apply_delta(
            old.m_fond_dissolved_kg,
            liq_c[7],
            "solute_fond",
            &mut shorten_dt,
        )?;
        let mut m_fond_adh = old.m_fond_adhered_kg;

        let mut m_h2o_vap = apply_delta(old.m_h2o_vap_kg, vap_c[0], "H2O_vap", &mut shorten_dt)?;
        let mut m_etoh_vap =
            apply_delta(old.m_etoh_vap_kg, vap_c[1], "C2H5OH_vap", &mut shorten_dt)?;
        let mut aux_gas = old.aux_gas_kg;
        for k in 0..AUX_GAS_COUNT {
            aux_gas[k] = apply_delta(old.aux_gas_kg[k], vap_c[2 + k], "aux_gas", &mut shorten_dt)?;
        }

        if let Some(dt_short) = shorten_dt {
            return Ok(Proposal::Shorten { dt_s: dt_short });
        }

        // 2. Conservative internal heat exchange between active vat and sub-tube cold zone.
        let q_active_to_cz_w = if self.config.cold_zone_enabled
            && old.m_oil_coldzone_kg > 1.0e-12
            && self.config.cold_zone_ua_w_k > 0.0
        {
            // Implicitly damped linear exchange prevents oscillation on large dt_s.
            let c_cz =
                coldzone_heat_capacity_j_k(old.m_oil_coldzone_kg, &self.config, old.t_coldzone_k);
            let c_act = active_heat_capacity_j_k(old, &self.config, old.t_liq_k).max(1.0);
            let ua = self.config.cold_zone_ua_w_k;
            let rate = ua * (1.0 / c_act + 1.0 / c_cz.max(1.0));
            let eff_factor = if rate * dt_s > 1.0e-9 {
                -expm1(-rate * dt_s) / (rate * dt_s)
            } else {
                1.0
            };
            ua * eff_factor * (old.t_liq_k - old.t_coldzone_k)
        } else {
            0.0
        };

        let port_energy_in_j = dt_s * (q_wall_w + q_food_sum_w + h_dot_liq_w + h_dot_vap_w);
        let q_cz_j = dt_s * q_active_to_cz_w;
        let active_enthalpy_new_j = old.active_enthalpy_j + port_energy_in_j - q_cz_j;
        let coldzone_enthalpy_new_j = old.coldzone_enthalpy_j + q_cz_j;

        // 3. Internal fond dissolution (m_fond_adhered -> m_fond_dissolved; mass & solute_fond invariant).
        let mut work = MultiphaseVesselState {
            time_s: old.time_s + dt_s,
            t_liq_k: old.t_liq_k,
            t_coldzone_k: old.t_coldzone_k,
            active_enthalpy_j: active_enthalpy_new_j,
            coldzone_enthalpy_j: coldzone_enthalpy_new_j,
            m_h2o_liq_kg: m_h2o_liq,
            m_etoh_liq_kg: m_etoh_liq,
            m_triglyceride_kg: m_trig,
            m_ffa_kg: m_ffa,
            m_glycerol_kg: m_glyc,
            m_salt_kg: m_salt,
            m_sugar_kg: m_sugar,
            m_fond_dissolved_kg: m_fond_diss,
            m_fond_adhered_kg: m_fond_adh,
            m_milk_solids_kg: old.m_milk_solids_kg,
            m_starch_kg: old.m_starch_kg,
            m_protein_kg: old.m_protein_kg,
            m_oil_coldzone_kg: old.m_oil_coldzone_kg,
            m_h2o_vap_kg: m_h2o_vap,
            m_etoh_vap_kg: m_etoh_vap,
            aux_gas_kg: aux_gas,
            tpm_fraction: old.tpm_fraction,
            foam_height_m: old.foam_height_m,
            milk_solid_browning: old.milk_solid_browning,
            cumulative_hydrolysis_mol: old.cumulative_hydrolysis_mol,
        };

        let f_wet_pre = self.wetting_fraction(&work);
        let m_solvent = m_h2o_liq + m_etoh_liq;
        if m_fond_adh > 0.0 && m_solvent > 1.0e-9 && f_wet_pre > 0.0 {
            let rho_l = self.liquid_density_kg_m3(&work);
            let vol_liq = (work.active_liquid_mass_kg() / rho_l).max(1.0e-12);
            let x_etoh = work.ethanol_mole_fraction();
            let c_fond = m_fond_diss / vol_liq;
            let c_sat = self.config.fond_c_sat_ref_kg_m3
                * (1.0 + 0.85 * ((old.t_liq_k - JANAF_REFERENCE_TEMPERATURE_K).max(0.0) / 75.0))
                * (1.0 + 0.65 * x_etoh);
            if c_sat > c_fond {
                let k_diss = self.config.fond_k_diss_ref_m_s
                    * exp(
                        -(self.config.fond_ea_diss_j_mol / UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K)
                            * (1.0 / old.t_liq_k - 1.0 / 353.15),
                    )
                    * (1.0 + 0.75 * sqrt(stir_rpm / 60.0))
                    * (1.0 + 0.80 * x_etoh);
                let dm_fond =
                    (f_wet_pre * k_diss * self.config.fond_area_m2 * (c_sat - c_fond) * dt_s)
                        .min(m_fond_adh)
                        .max(0.0);
                m_fond_adh -= dm_fond;
                m_fond_diss += dm_fond;
                work.m_fond_adhered_kg = m_fond_adh;
                work.m_fond_dissolved_kg = m_fond_diss;
            }
        }

        // 4. Stoichiometric frying-oil hydrolysis + TPM aging ODEs (before thermal inversion so
        //    reaction enthalpy is captured identically in JANAF single-accounting total enthalpy).
        let m_steam_food_rate_in = liq_c[0].max(0.0) + vap_c[0].max(0.0);
        if let Some(oil) = &self.config.oil {
            let m_lip_act = m_trig + m_ffa + m_glyc;
            if m_lip_act > 1.0e-9 {
                let t_kin = old.t_liq_k.clamp(273.15, 800.0);
                let rt = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * t_kin;
                let k_hyd = oil.k_hyd_s_inv * exp(-oil.ea_hyd_j_mol / rt);
                let k_ox = oil.k_ox_s_inv * exp(-oil.ea_ox_j_mol / rt);
                let k_poly = oil.k_poly_s_inv * exp(-oil.ea_poly_j_mol / rt);

                if self.config.hydrolysis_enabled {
                    let m_h2o_avail = m_h2o_liq + m_h2o_vap;
                    if m_trig > 0.0 && m_h2o_avail > 0.0 && m_steam_food_rate_in > 0.0 {
                        let r_ffa_per_s = k_hyd
                            * oil.hydrolysis_contact_time_s
                            * (m_steam_food_rate_in / m_lip_act);
                        let dm_ffa_target = r_ffa_per_s * dt_s * m_lip_act;
                        let xi_target = dm_ffa_target / (3.0 * MOLAR_MASS_FFA_KG_MOL);
                        let xi_max_trig = m_trig / MOLAR_MASS_TRIGLYCERIDE_KG_MOL;
                        let xi_max_h2o = m_h2o_avail / (3.0 * MOLAR_MASS_H2O_KG_MOL);
                        let d_xi_mol = xi_target.min(xi_max_trig).min(xi_max_h2o).max(0.0);

                        if d_xi_mol > 0.0 {
                            let dm_h2o_rxn = 3.0 * d_xi_mol * MOLAR_MASS_H2O_KG_MOL;
                            let from_liq = dm_h2o_rxn.min(m_h2o_liq);
                            m_h2o_liq = (m_h2o_liq - from_liq).max(0.0);
                            let rem_h2o = (dm_h2o_rxn - from_liq).max(0.0);
                            m_h2o_vap = (m_h2o_vap - rem_h2o).max(0.0);

                            m_trig = (m_trig - d_xi_mol * MOLAR_MASS_TRIGLYCERIDE_KG_MOL).max(0.0);
                            m_ffa += 3.0 * d_xi_mol * MOLAR_MASS_FFA_KG_MOL;
                            m_glyc += d_xi_mol * MOLAR_MASS_GLYCEROL_KG_MOL;

                            work.m_h2o_liq_kg = m_h2o_liq;
                            work.m_h2o_vap_kg = m_h2o_vap;
                            work.m_triglyceride_kg = m_trig;
                            work.m_ffa_kg = m_ffa;
                            work.m_glycerol_kg = m_glyc;
                            work.cumulative_hydrolysis_mol =
                                old.cumulative_hydrolysis_mol + d_xi_mol;
                        }
                    }
                }

                // TPM oxidation + polymerization ODE (only significant at frying temperatures).
                let ffa_frac_now = work.ffa_fraction();
                let d_tpm = (k_ox
                    * oil.oxidation_surface_factor_kg_m2
                    * (self.config.free_surface_area_m2 / m_lip_act)
                    + k_poly * ffa_frac_now)
                    * dt_s;
                work.tpm_fraction = (old.tpm_fraction + d_tpm).clamp(0.0, 0.99);
            }
        }

        // 5. Trial sensible temperature inversion, then equilibrium VLE vaporization with
        //    smooth dryout regularization f_wet = tanh(m_vol / m_wet_ref) and Proposal::Shorten.
        let t_trial =
            invert_active_temperature(&work, &self.config, active_enthalpy_new_j, old.t_liq_k)?;
        let mut dm_evap_total = 0.0;
        let m_vol = m_h2o_liq + m_etoh_liq;
        if m_vol > 0.0 {
            let (t_bub_k, y_etoh_mol, _dt_bpe) =
                self.bubble_point_and_vapor_composition(&work, p_vessel_pa)?;
            if t_trial > t_bub_k {
                let h_at_bub = compute_active_enthalpy_j(&work, &self.config, t_bub_k);
                let e_excess_j = (active_enthalpy_new_j - h_at_bub).max(0.0);
                let w_etoh_v = if y_etoh_mol > 0.0 {
                    let num = y_etoh_mol * MOLAR_MASS_ETOH_KG_MOL;
                    let den = num + (1.0 - y_etoh_mol) * MOLAR_MASS_H2O_KG_MOL;
                    (num / den).clamp(0.0, 1.0)
                } else {
                    0.0
                };
                let w_h2o_v = 1.0 - w_etoh_v;
                let h_fg_h2o =
                    specific_enthalpy_h2o_vap(t_bub_k) - specific_enthalpy_h2o_liq(t_bub_k);
                let h_fg_etoh =
                    specific_enthalpy_etoh_vap(t_bub_k) - specific_enthalpy_etoh_liq(t_bub_k);
                let h_fg_mix = (w_etoh_v * h_fg_etoh + w_h2o_v * h_fg_h2o).max(1.0e5);

                // Smooth wetting regularization on the volatile phase.
                let rho_l = self.liquid_density_kg_m3(&work);
                let m_wet_ref =
                    (rho_l * self.config.base_area_m2 * self.config.dryout_delta_min_m).max(1.0e-9);
                let f_wet_vol = tanh(m_vol / m_wet_ref);
                let dm_evap_raw = f_wet_vol * (e_excess_j / h_fg_mix);

                // Emit Proposal::Shorten on true vessel dryout approach (volatile-dominated bath
                // where a single step would consume > 55% of remaining volatile inventory).
                // Lipid-dominated fryers/butter pans flashing trace water stay wet in oil and
                // clamp flashed water cleanly without spurious micro-steps.
                let m_active_liq = work.active_liquid_mass_kg();
                let is_volatile_dominated = work.active_lipid_mass_kg() < 0.25 * m_active_liq;
                if is_volatile_dominated
                    && m_vol > self.config.dryout_floor_kg
                    && dm_evap_raw > 0.55 * m_vol
                    && dt_s > self.config.min_dt_s
                {
                    let dt_short = (0.45 * dt_s * (m_vol / dm_evap_raw))
                        .clamp(self.config.min_dt_s, 0.5 * dt_s);
                    if dt_short < dt_s {
                        return Ok(Proposal::Shorten { dt_s: dt_short });
                    }
                }

                // Clamp strictly to available liquid inventories (zero negative mass ever).
                let dm_etoh = (dm_evap_raw * w_etoh_v).min(m_etoh_liq).max(0.0);
                let dm_h2o = (dm_evap_raw * w_h2o_v).min(m_h2o_liq).max(0.0);
                m_etoh_liq = (m_etoh_liq - dm_etoh).max(0.0);
                m_h2o_liq = (m_h2o_liq - dm_h2o).max(0.0);
                m_etoh_vap += dm_etoh;
                m_h2o_vap += dm_h2o;
                dm_evap_total = dm_etoh + dm_h2o;

                work.m_etoh_liq_kg = m_etoh_liq;
                work.m_h2o_liq_kg = m_h2o_liq;
                work.m_etoh_vap_kg = m_etoh_vap;
                work.m_h2o_vap_kg = m_h2o_vap;
            }
        }

        // 6. Final exact analytical temperature inversion on active & cold zones.
        work.t_liq_k =
            invert_active_temperature(&work, &self.config, active_enthalpy_new_j, t_trial)?;
        work.t_coldzone_k = if self.config.cold_zone_enabled && work.m_oil_coldzone_kg > 1.0e-12 {
            invert_coldzone_temperature(
                &work,
                &self.config,
                coldzone_enthalpy_new_j,
                old.t_coldzone_k,
            )?
        } else {
            work.t_liq_k
        };

        // 7. Milk-solid *beurre noisette* Maillard browning kinetics (activates smoothly above ~115 C
        //    once butter water flashes off and bath temperature climbs into 120-150 C).
        if work.m_milk_solids_kg > 1.0e-9 && work.t_liq_k > 383.15 {
            let t_mid = 0.5 * (old.t_liq_k + work.t_liq_k);
            // Arrhenius Maillard milk-protein/lactose browning (Ea ~ 105 kJ/mol, tau ~ 90 s at 135 C).
            let k_brown = 2.8e11 * exp(-105_000.0 / (UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * t_mid));
            let frac_step = -expm1(-k_brown * dt_s);
            work.milk_solid_browning = (old.milk_solid_browning
                + (1.0 - old.milk_solid_browning) * frac_step)
                .clamp(0.0, 1.0);
        }

        // 8. Bikerman boil-over foam ODE: tau_foam * dH_foam/dt = Sigma_foam(c_starch, c_prot) * j_g - H_foam.
        //    In lipid-continuous baths (e.g. melted butter or frying oil), liquid triglyceride
        //    acts as a lamella-breaking antifoam relative to aqueous pasta/milk boils, and
        //    liquid-inventory conservation caps a ~5%-liquid-fraction foam at <= 10 * H_liq.
        let m_liq_now = work.active_liquid_mass_kg();
        let rho_l_now = self.liquid_density_kg_m3(&work);
        let vol_liq_now = (m_liq_now / rho_l_now).max(1.0e-9);
        let c_starch = work.m_starch_kg / vol_liq_now;
        let c_protein = (work.m_protein_kg + work.m_milk_solids_kg) / vol_liq_now;
        let w_lipid_now = if m_liq_now > 1.0e-15 {
            (work.active_lipid_mass_kg() / m_liq_now).clamp(0.0, 1.0)
        } else {
            0.0
        };
        let lipid_antifoam = 1.0 - 0.85 * w_lipid_now;
        let sigma_foam = lipid_antifoam
            * (self.config.foam_sigma_starch_s
                * (c_starch / (c_starch + self.config.foam_c_starch_half_kg_m3))
                + self.config.foam_sigma_protein_s
                    * (c_protein / (c_protein + self.config.foam_c_protein_half_kg_m3)));

        let rho_vap_est = (p_vessel_pa * MOLAR_MASS_H2O_KG_MOL
            / (UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * work.t_liq_k))
            .max(0.05);
        let m_dot_vapor_bubbling = dm_evap_total / dt_s + m_steam_food_rate_in;
        let j_g = m_dot_vapor_bubbling / (rho_vap_est * self.config.free_surface_area_m2);
        let h_liq_now = vol_liq_now / self.config.free_surface_area_m2;
        let h_foam_eq = (sigma_foam * j_g).min(10.0 * h_liq_now);
        let alpha_foam = -expm1(-dt_s / self.config.foam_tau_s);
        work.foam_height_m =
            (old.foam_height_m + alpha_foam * (h_foam_eq - old.foam_height_m)).max(0.0);

        work.validate(&self.config)?;
        *new = work;
        self.write_port_outputs(old, new, q_wall_w, stir_rpm, p_vessel_pa, outputs);

        let e_old = old.total_enthalpy_j();
        let e_new = new.total_enthalpy_j();
        let energy_tol =
            256.0 * f64::EPSILON * (e_old.abs() + e_new.abs() + port_energy_in_j.abs() + 1.0);
        let receipt = Receipt {
            stored_change_j: e_new - e_old,
            dissipated_j: 0.0,
            numerical_loss_j: 0.0,
            entropy_generated_j_per_k: 0.0,
            energy_tolerance_j: energy_tol,
            residual: 0.0,
            tolerance: 0.0,
        };
        Ok(Proposal::Step(receipt))
    }

    fn admit(
        &self,
        _old: &Self::State,
        new: &Self::State,
        _receipt: &Receipt,
    ) -> Result<(), Self::Error> {
        new.validate(&self.config)
    }

    const STATE_CHARTER: Option<StateIdentityCharterV2> = Some(StateIdentityCharterV2 {
        owner: "fs-multiphase::MultiphaseVessel",
        state_family: "stirred-vle-v1",
        schema_grammar: "time_s, t_liq_k, t_coldzone_k, active_enthalpy_j, coldzone_enthalpy_j, m_h2o_liq_kg, m_etoh_liq_kg, m_triglyceride_kg, m_ffa_kg, m_glycerol_kg, m_salt_kg, m_sugar_kg, m_fond_dissolved_kg, m_fond_adhered_kg, m_milk_solids_kg, m_starch_kg, m_protein_kg, m_oil_coldzone_kg, m_h2o_vap_kg, m_etoh_vap_kg, aux_gas_kg[7], tpm_fraction, foam_height_m, milk_solid_browning, cumulative_hydrolysis_mol in coherent SI on janaf-298.15K-formation",
        codec_grammar: "little-endian u64 words: codec tag, field count (31), followed by 31 IEEE-754 f64 bit patterns in declaration order",
        codec_version: 1,
    });

    fn encode_state(&self, state: &Self::State, words: &mut Vec<u64>) -> bool {
        if state.validate(&self.config).is_err() {
            return false;
        }
        words.push(VESSEL_STATE_CODEC_TAG);
        words.push(VESSEL_STATE_SCALAR_COUNT as u64);
        let scalars = [
            state.time_s,
            state.t_liq_k,
            state.t_coldzone_k,
            state.active_enthalpy_j,
            state.coldzone_enthalpy_j,
            state.m_h2o_liq_kg,
            state.m_etoh_liq_kg,
            state.m_triglyceride_kg,
            state.m_ffa_kg,
            state.m_glycerol_kg,
            state.m_salt_kg,
            state.m_sugar_kg,
            state.m_fond_dissolved_kg,
            state.m_fond_adhered_kg,
            state.m_milk_solids_kg,
            state.m_starch_kg,
            state.m_protein_kg,
            state.m_oil_coldzone_kg,
            state.m_h2o_vap_kg,
            state.m_etoh_vap_kg,
            state.aux_gas_kg[0],
            state.aux_gas_kg[1],
            state.aux_gas_kg[2],
            state.aux_gas_kg[3],
            state.aux_gas_kg[4],
            state.aux_gas_kg[5],
            state.aux_gas_kg[6],
            state.tpm_fraction,
            state.foam_height_m,
            state.milk_solid_browning,
            state.cumulative_hydrolysis_mol,
        ];
        words.extend(scalars.iter().map(|v| v.to_bits()));
        true
    }

    fn decode_state(&self, words: &[u64]) -> Option<Self::State> {
        if words.len() != 2 + VESSEL_STATE_SCALAR_COUNT {
            return None;
        }
        if words[0] != VESSEL_STATE_CODEC_TAG || words[1] != VESSEL_STATE_SCALAR_COUNT as u64 {
            return None;
        }
        let s = |idx: usize| f64::from_bits(words[2 + idx]);
        let state = MultiphaseVesselState {
            time_s: s(0),
            t_liq_k: s(1),
            t_coldzone_k: s(2),
            active_enthalpy_j: s(3),
            coldzone_enthalpy_j: s(4),
            m_h2o_liq_kg: s(5),
            m_etoh_liq_kg: s(6),
            m_triglyceride_kg: s(7),
            m_ffa_kg: s(8),
            m_glycerol_kg: s(9),
            m_salt_kg: s(10),
            m_sugar_kg: s(11),
            m_fond_dissolved_kg: s(12),
            m_fond_adhered_kg: s(13),
            m_milk_solids_kg: s(14),
            m_starch_kg: s(15),
            m_protein_kg: s(16),
            m_oil_coldzone_kg: s(17),
            m_h2o_vap_kg: s(18),
            m_etoh_vap_kg: s(19),
            aux_gas_kg: [s(20), s(21), s(22), s(23), s(24), s(25), s(26)],
            tpm_fraction: s(27),
            foam_height_m: s(28),
            milk_solid_browning: s(29),
            cumulative_hydrolysis_mol: s(30),
        };
        state.validate(&self.config).ok().map(|()| state)
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        let c = &self.config;
        words.push(c.immersed_batch_count as u64);
        words.push(u64::from(c.cold_zone_enabled));
        words.push(u64::from(c.hydrolysis_enabled));
        words.push(u64::from(c.oil.is_some()));
        let scalars = [
            c.base_area_m2,
            c.free_surface_area_m2,
            c.vessel_height_m,
            c.dryout_delta_min_m,
            c.dryout_floor_kg,
            c.min_dt_s,
            c.dry_vessel_heat_capacity_j_k,
            c.fond_area_m2,
            c.fond_k_diss_ref_m_s,
            c.fond_ea_diss_j_mol,
            c.fond_c_sat_ref_kg_m3,
            c.foam_tau_s,
            c.foam_sigma_starch_s,
            c.foam_c_starch_half_kg_m3,
            c.foam_sigma_protein_s,
            c.foam_c_protein_half_kg_m3,
            c.sauce_viscosity_conc_coeff,
            c.sauce_carreau_n,
            c.sauce_carreau_lambda_s,
            c.natural_wall_htc_w_m2_k,
            c.natural_immersed_htc_w_m2_k,
            c.stir_htc_gain,
            c.bubble_htc_gain,
            c.bubble_jg_ref_m_s,
            c.scorch_wall_threshold_k,
            c.cold_zone_ua_w_k,
        ];
        words.extend(scalars.iter().map(|v| v.to_bits()));
        if let Some(oil) = &c.oil {
            oil.append_identity_words(words);
        }
    }
}

fn build_vessel_ports(
    config: &MultiphaseVesselConfig,
) -> Result<(Vec<PortSpec>, PortIndices), VesselErrorV1> {
    let mk_id = |name: &'static str| {
        StableId::new(name).map_err(|_| VesselErrorV1::InvalidConfig {
            parameter: name,
            reason: "failed to construct StableId",
        })
    };

    let mut ports = Vec::with_capacity(7);

    let wall_heat = ports.len();
    ports.push(PortSpec::heat(mk_id("wall_heat")?, PortRole::FlowIn)?);

    let immersed_food_heat = ports.len();
    let batch_width =
        NonZeroUsize::new(config.immersed_batch_count).ok_or(VesselErrorV1::InvalidConfig {
            parameter: "immersed_batch_count",
            reason: "must be >= 1",
        })?;
    ports.push(PortSpec::heat_field(
        mk_id("immersed_food_heat")?,
        PortRole::FlowIn,
        batch_width,
    )?);

    let liquid_stream = ports.len();
    ports.push(PortSpec::identified_stream(
        mk_id("liquid_stream")?,
        PortRole::FlowIn,
        kitchen_liquid_v1()?,
    )?);

    let vapor_stream = ports.len();
    ports.push(PortSpec::identified_stream(
        mk_id("vapor_stream")?,
        PortRole::FlowIn,
        kitchen_gas_v1()?,
    )?);

    let stir_rpm = ports.len();
    ports.push(PortSpec::signal_in(mk_id("stir_rpm")?, NonZeroUsize::MIN));

    let p_vessel_pa = ports.len();
    ports.push(PortSpec::signal_in(
        mk_id("p_vessel_pa")?,
        NonZeroUsize::MIN,
    ));

    let probes = ports.len();
    let probe_width =
        NonZeroUsize::new(VESSEL_PROBE_COUNT).ok_or(VesselErrorV1::InvalidConfig {
            parameter: "probes",
            reason: "probe width must be nonzero",
        })?;
    ports.push(PortSpec::signal_out(mk_id("probes")?, probe_width));

    Ok((
        ports,
        PortIndices {
            wall_heat,
            immersed_food_heat,
            liquid_stream,
            vapor_stream,
            stir_rpm,
            p_vessel_pa,
            probes,
        },
    ))
}

fn build_hydrolysis_reaction_basis()
-> Result<(ElementalMatrix, Vec<(SpeciesId, MolarMass)>), VesselErrorV1> {
    let el = |s: &str| {
        ElementId::new(s).map_err(|_| VesselErrorV1::InvalidConfig {
            parameter: "reaction_basis",
            reason: "invalid element symbol",
        })
    };
    let sp = |s: &str| {
        SpeciesId::new(s).map_err(|_| VesselErrorV1::InvalidConfig {
            parameter: "reaction_basis",
            reason: "invalid species identifier",
        })
    };
    // Exact triolein hydrolysis: C57H104O6 (triglyceride) + 3 H2O -> 3 C18H34O2 (ffa) + C3H8O3 (glycerol).
    // Canonical ASCII-sorted species order: [H2O, ffa, glycerol, triglyceride].
    let elemental = ElementalMatrix::new(
        vec![el("C")?, el("H")?, el("O")?],
        vec![sp("H2O")?, sp("ffa")?, sp("glycerol")?, sp("triglyceride")?],
        vec![
            vec![0, 18, 3, 57],  // C
            vec![2, 34, 8, 104], // H
            vec![1, 2, 3, 6],    // O
        ],
    )
    .map_err(|_| VesselErrorV1::InvalidConfig {
        parameter: "reaction_basis",
        reason: "failed to construct triglyceride hydrolysis ElementalMatrix",
    })?;
    let molar_masses = vec![
        (sp("H2O")?, MolarMass::new(MOLAR_MASS_H2O_KG_MOL)),
        (sp("ffa")?, MolarMass::new(MOLAR_MASS_FFA_KG_MOL)),
        (sp("glycerol")?, MolarMass::new(MOLAR_MASS_GLYCEROL_KG_MOL)),
        (
            sp("triglyceride")?,
            MolarMass::new(MOLAR_MASS_TRIGLYCERIDE_KG_MOL),
        ),
    ];
    Ok((elemental, molar_masses))
}

fn validate_stream_mass_consistency(port: &'static str, wire: &[f64]) -> Result<(), VesselErrorV1> {
    let mdot_tot = wire[STREAM_MASS_RATE];
    let carried = &wire[STREAM_CARRIED..];
    let species_sum: f64 = carried.iter().sum();
    let abs_scale = mdot_tot.abs() + carried.iter().map(|v| v.abs()).sum::<f64>() + 1.0e-12;
    if (mdot_tot - species_sum).abs() > 1.0e-10 * abs_scale {
        return Err(VesselErrorV1::InconsistentStreamMassFlow {
            port,
            declared_total_kg_s: mdot_tot,
            species_sum_kg_s: species_sum,
        });
    }
    Ok(())
}

// --- Exact JANAF-298.15K-formation specific enthalpies [J/kg] ---

#[inline]
fn specific_enthalpy_h2o_liq(t_k: f64) -> f64 {
    H_F0_H2O_LIQ_J_KG + CP_H2O_LIQ_J_KG_K * (t_k - JANAF_REFERENCE_TEMPERATURE_K)
}

#[inline]
fn specific_enthalpy_h2o_vap(t_k: f64) -> f64 {
    H_F0_H2O_VAP_J_KG + CP_H2O_VAP_J_KG_K * (t_k - JANAF_REFERENCE_TEMPERATURE_K)
}

#[inline]
fn specific_enthalpy_etoh_liq(t_k: f64) -> f64 {
    H_F0_ETOH_LIQ_J_KG + CP_ETOH_LIQ_J_KG_K * (t_k - JANAF_REFERENCE_TEMPERATURE_K)
}

#[inline]
fn specific_enthalpy_etoh_vap(t_k: f64) -> f64 {
    H_F0_ETOH_VAP_J_KG + CP_ETOH_VAP_J_KG_K * (t_k - JANAF_REFERENCE_TEMPERATURE_K)
}

#[inline]
fn specific_enthalpy_trig(config: &MultiphaseVesselConfig, t_k: f64) -> f64 {
    let (cp_ref, dcp_dt, t_ref) = config.lipid_cp_params();
    let u_ref = JANAF_REFERENCE_TEMPERATURE_K - t_ref;
    let u = t_k - t_ref;
    let delta_h = cp_ref * (u - u_ref) + 0.5 * dcp_dt * (u * u - u_ref * u_ref);
    H_F0_TRIGLYCERIDE_J_KG + delta_h
}

#[inline]
fn specific_enthalpy_ffa(t_k: f64) -> f64 {
    H_F0_FFA_J_KG + CP_FFA_J_KG_K * (t_k - JANAF_REFERENCE_TEMPERATURE_K)
}

#[inline]
fn specific_enthalpy_glycerol(t_k: f64) -> f64 {
    H_F0_GLYCEROL_J_KG + CP_GLYCEROL_J_KG_K * (t_k - JANAF_REFERENCE_TEMPERATURE_K)
}

#[inline]
fn specific_enthalpy_salt(t_k: f64) -> f64 {
    H_F0_SALT_J_KG + CP_SALT_J_KG_K * (t_k - JANAF_REFERENCE_TEMPERATURE_K)
}

#[inline]
fn specific_enthalpy_sugar(t_k: f64) -> f64 {
    H_F0_SUGAR_J_KG + CP_SUGAR_J_KG_K * (t_k - JANAF_REFERENCE_TEMPERATURE_K)
}

#[inline]
fn specific_enthalpy_fond(t_k: f64) -> f64 {
    H_F0_FOND_J_KG + CP_FOND_J_KG_K * (t_k - JANAF_REFERENCE_TEMPERATURE_K)
}

#[inline]
fn specific_enthalpy_milk_solids(t_k: f64) -> f64 {
    H_F0_MILK_SOLIDS_J_KG + CP_MILK_SOLIDS_J_KG_K * (t_k - JANAF_REFERENCE_TEMPERATURE_K)
}

#[inline]
fn specific_enthalpy_starch(t_k: f64) -> f64 {
    H_F0_STARCH_J_KG + CP_STARCH_J_KG_K * (t_k - JANAF_REFERENCE_TEMPERATURE_K)
}

#[inline]
fn specific_enthalpy_protein(t_k: f64) -> f64 {
    H_F0_PROTEIN_J_KG + CP_PROTEIN_J_KG_K * (t_k - JANAF_REFERENCE_TEMPERATURE_K)
}

#[inline]
fn specific_enthalpy_aux_gas(index: usize, t_k: f64) -> f64 {
    AUX_GAS_H_F0_J_KG[index] + AUX_GAS_CP_J_KG_K[index] * (t_k - JANAF_REFERENCE_TEMPERATURE_K)
}

/// Active-zone enthalpy polynomial coefficients `(H_298, C_1, C_2)` such that
/// `H(T) = H_298 + C_1 * u + 0.5 * C_2 * u^2` where `u = T - 298.15 K`.
fn active_enthalpy_polynomial(
    state: &MultiphaseVesselState,
    config: &MultiphaseVesselConfig,
) -> (f64, f64, f64) {
    let (cp_ref, dcp_dt, t_ref) = config.lipid_cp_params();
    // Shift oil cp(T) = cp_ref + dcp_dt*(T - t_ref) = cp_at_298 + dcp_dt*(T - 298.15).
    let cp_trig_298 = cp_ref + dcp_dt * (JANAF_REFERENCE_TEMPERATURE_K - t_ref);

    let mut h0 = state.m_h2o_liq_kg * H_F0_H2O_LIQ_J_KG
        + state.m_etoh_liq_kg * H_F0_ETOH_LIQ_J_KG
        + state.m_triglyceride_kg * H_F0_TRIGLYCERIDE_J_KG
        + state.m_ffa_kg * H_F0_FFA_J_KG
        + state.m_glycerol_kg * H_F0_GLYCEROL_J_KG
        + state.m_salt_kg * H_F0_SALT_J_KG
        + state.m_sugar_kg * H_F0_SUGAR_J_KG
        + (state.m_fond_dissolved_kg + state.m_fond_adhered_kg) * H_F0_FOND_J_KG
        + state.m_milk_solids_kg * H_F0_MILK_SOLIDS_J_KG
        + state.m_starch_kg * H_F0_STARCH_J_KG
        + state.m_protein_kg * H_F0_PROTEIN_J_KG
        + state.m_h2o_vap_kg * H_F0_H2O_VAP_J_KG
        + state.m_etoh_vap_kg * H_F0_ETOH_VAP_J_KG;

    let mut c1 = config.dry_vessel_heat_capacity_j_k
        + state.m_h2o_liq_kg * CP_H2O_LIQ_J_KG_K
        + state.m_etoh_liq_kg * CP_ETOH_LIQ_J_KG_K
        + state.m_triglyceride_kg * cp_trig_298
        + state.m_ffa_kg * CP_FFA_J_KG_K
        + state.m_glycerol_kg * CP_GLYCEROL_J_KG_K
        + state.m_salt_kg * CP_SALT_J_KG_K
        + state.m_sugar_kg * CP_SUGAR_J_KG_K
        + (state.m_fond_dissolved_kg + state.m_fond_adhered_kg) * CP_FOND_J_KG_K
        + state.m_milk_solids_kg * CP_MILK_SOLIDS_J_KG_K
        + state.m_starch_kg * CP_STARCH_J_KG_K
        + state.m_protein_kg * CP_PROTEIN_J_KG_K
        + state.m_h2o_vap_kg * CP_H2O_VAP_J_KG_K
        + state.m_etoh_vap_kg * CP_ETOH_VAP_J_KG_K;

    for k in 0..AUX_GAS_COUNT {
        h0 += state.aux_gas_kg[k] * AUX_GAS_H_F0_J_KG[k];
        c1 += state.aux_gas_kg[k] * AUX_GAS_CP_J_KG_K[k];
    }

    let c2 = state.m_triglyceride_kg * dcp_dt;
    (h0, c1, c2)
}

fn compute_active_enthalpy_j(
    state: &MultiphaseVesselState,
    config: &MultiphaseVesselConfig,
    t_k: f64,
) -> f64 {
    let (h0, c1, c2) = active_enthalpy_polynomial(state, config);
    let u = t_k - JANAF_REFERENCE_TEMPERATURE_K;
    h0 + c1 * u + 0.5 * c2 * u * u
}

fn active_heat_capacity_j_k(
    state: &MultiphaseVesselState,
    config: &MultiphaseVesselConfig,
    t_k: f64,
) -> f64 {
    let (_, c1, c2) = active_enthalpy_polynomial(state, config);
    (c1 + c2 * (t_k - JANAF_REFERENCE_TEMPERATURE_K)).max(0.0)
}

fn compute_coldzone_enthalpy_j(
    state: &MultiphaseVesselState,
    config: &MultiphaseVesselConfig,
    t_k: f64,
) -> f64 {
    if state.m_oil_coldzone_kg <= 0.0 {
        return 0.0;
    }
    state.m_oil_coldzone_kg * specific_enthalpy_trig(config, t_k)
}

fn coldzone_heat_capacity_j_k(m_oil_cz_kg: f64, config: &MultiphaseVesselConfig, t_k: f64) -> f64 {
    let (cp_ref, dcp_dt, t_ref) = config.lipid_cp_params();
    m_oil_cz_kg * (cp_ref + dcp_dt * (t_k - t_ref)).max(100.0)
}

fn invert_quadratic_enthalpy(
    h0: f64,
    c1: f64,
    c2: f64,
    h_target: f64,
    fallback_t_k: f64,
    stage: &'static str,
) -> Result<f64, VesselErrorV1> {
    let dh = h_target - h0;
    if c1 <= 1.0e-12 && c2 <= 1.0e-12 {
        if dh.abs() <= 1.0e-9 {
            return Ok(fallback_t_k);
        }
        return Err(VesselErrorV1::EmptyVesselThermalSingularity);
    }
    let u = if c2 <= 1.0e-12 {
        dh / c1
    } else {
        let disc = c1 * c1 + 2.0 * c2 * dh;
        if disc < 0.0 {
            return Err(VesselErrorV1::SolverConvergence {
                stage,
                residual: -disc,
                tolerance: 0.0,
            });
        }
        // Cancellation-free Citardauq root for 0.5*c2*u^2 + c1*u - dh = 0.
        2.0 * dh / (c1 + sqrt(disc))
    };
    let t_k = JANAF_REFERENCE_TEMPERATURE_K + u;
    if !t_k.is_finite() || !(TEMP_MIN_VALID_K..=TEMP_MAX_VALID_K).contains(&t_k) {
        return Err(VesselErrorV1::SolverConvergence {
            stage,
            residual: t_k,
            tolerance: TEMP_MAX_VALID_K,
        });
    }
    Ok(t_k)
}

fn invert_active_temperature(
    state: &MultiphaseVesselState,
    config: &MultiphaseVesselConfig,
    h_target: f64,
    fallback_t_k: f64,
) -> Result<f64, VesselErrorV1> {
    let (h0, c1, c2) = active_enthalpy_polynomial(state, config);
    invert_quadratic_enthalpy(h0, c1, c2, h_target, fallback_t_k, "active_temperature")
}

fn invert_coldzone_temperature(
    state: &MultiphaseVesselState,
    config: &MultiphaseVesselConfig,
    h_target: f64,
    fallback_t_k: f64,
) -> Result<f64, VesselErrorV1> {
    if state.m_oil_coldzone_kg <= 1.0e-15 {
        return Ok(fallback_t_k);
    }
    let (cp_ref, dcp_dt, t_ref) = config.lipid_cp_params();
    let cp_298 = cp_ref + dcp_dt * (JANAF_REFERENCE_TEMPERATURE_K - t_ref);
    let h0 = state.m_oil_coldzone_kg * H_F0_TRIGLYCERIDE_J_KG;
    let c1 = state.m_oil_coldzone_kg * cp_298;
    let c2 = state.m_oil_coldzone_kg * dcp_dt;
    invert_quadratic_enthalpy(h0, c1, c2, h_target, fallback_t_k, "coldzone_temperature")
}

/// Expose auxiliary gas molar mass lookup for diagnostics (compile-checked constant use).
#[must_use]
pub const fn aux_gas_molar_mass_kg_mol(index: usize) -> Option<f64> {
    if index < AUX_GAS_COUNT {
        Some(AUX_GAS_MOLAR_MASS_KG_MOL[index])
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

    fn with_cx_and_gate(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
        let gate = CancelGate::new_clock_free();
        let pool = ArenaPool::new(ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 1,
                    kernel_id: 1,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            f(&gate, &cx);
        });
    }

    #[test]
    fn hydrolysis_stoichiometry_closes_mass_and_atoms_to_machine_precision() {
        let lhs = MOLAR_MASS_TRIGLYCERIDE_KG_MOL + 3.0 * MOLAR_MASS_H2O_KG_MOL;
        let rhs = 3.0 * MOLAR_MASS_FFA_KG_MOL + MOLAR_MASS_GLYCEROL_KG_MOL;
        assert!((lhs - rhs).abs() < 1.0e-15);

        let (elem, masses) = build_hydrolysis_reaction_basis().expect("basis");
        assert_eq!(elem.elements().len(), 3);
        assert_eq!(elem.species().len(), 4);
        assert_eq!(masses.len(), 4);
    }

    #[test]
    fn quadratic_enthalpy_inversion_round_trips_water_and_oil_to_machine_precision() {
        let cfg = MultiphaseVesselConfig::deep_fryer_default(OilCardV1::canola_refined());
        let comp = MultiphaseVesselComponent::new(cfg.clone()).expect("fryer");
        let state = comp
            .initial_state(VesselInitialStateV1::oil_fryer(
                &OilCardV1::canola_refined(),
                12.0,
                3.5,
                450.15,
                355.15,
            ))
            .expect("state");
        let t_act =
            invert_active_temperature(&state, &cfg, state.active_enthalpy_j, 300.0).expect("inv");
        let t_cz = invert_coldzone_temperature(&state, &cfg, state.coldzone_enthalpy_j, 300.0)
            .expect("inv_cz");
        assert!((t_act - 450.15).abs() < 1.0e-12);
        assert!((t_cz - 355.15).abs() < 1.0e-12);
    }

    #[test]
    fn codec_round_trip_and_identity_words_are_deterministic() {
        let cfg = MultiphaseVesselConfig::saucepan_default();
        let comp = MultiphaseVesselComponent::new(cfg).expect("vessel");
        let state = comp
            .initial_state(VesselInitialStateV1::wine_deglaze(0.5, 0.11, 0.025, 355.0))
            .expect("state");

        let mut words = Vec::new();
        assert!(comp.encode_state(&state, &mut words));
        let decoded = comp.decode_state(&words).expect("decode");
        assert_eq!(decoded, state);

        // Truncated or corrupted tag is rejected cleanly.
        assert!(comp.decode_state(&words[..words.len() - 1]).is_none());
        let mut bad = words.clone();
        bad[0] ^= 1;
        assert!(comp.decode_state(&bad).is_none());

        let mut id1 = Vec::new();
        let mut id2 = Vec::new();
        comp.identity_words(&mut id1);
        comp.identity_words(&mut id2);
        assert_eq!(id1, id2);
        assert!(!id1.is_empty());
    }

    #[test]
    fn stirring_boosts_curd_mixing_htc_three_to_five_fold() {
        let comp =
            MultiphaseVesselComponent::new(MultiphaseVesselConfig::skillet_default()).expect("pan");
        let state = comp
            .initial_state(VesselInitialStateV1::sauce(
                0.25, 0.0, 0.01, 0.002, 0.0, 345.0,
            ))
            .expect("state");
        let h_still = comp.effective_wall_htc_w_m2_k(&state, 0.0, 0.0);
        let h_whisk_120 = comp.effective_wall_htc_w_m2_k(&state, 120.0, 0.0);
        let h_whisk_180 = comp.effective_wall_htc_w_m2_k(&state, 180.0, 0.0);
        let ratio_120 = h_whisk_120 / h_still;
        let ratio_180 = h_whisk_180 / h_still;
        assert!(
            (3.0..=5.5).contains(&ratio_120),
            "expected 3-5x boost at 120 RPM, got {ratio_120}"
        );
        assert!(ratio_180 > ratio_120);
    }

    #[test]
    fn propose_refuses_non_positive_pressure_and_cancelled_context() {
        let comp =
            MultiphaseVesselComponent::new(MultiphaseVesselConfig::saucepan_default()).expect("v");
        let state = comp
            .initial_state(VesselInitialStateV1::water(1.0, 350.0))
            .expect("s");
        with_cx_and_gate(|gate, cx| {
            let bad_inputs = comp
                .build_inputs(500.0, &[0.0], [0.0; 10], [0.0; 11], 0.0, 0.0)
                .expect("inputs");
            let mut next = state.clone();
            let mut out = vec![0.0; comp.layout().output_len()];
            let err = comp
                .propose(&state, 1.0, &bad_inputs, &mut next, &mut out, cx)
                .expect_err("zero pressure must refuse");
            assert!(matches!(err, VesselErrorV1::InvalidPressure { .. }));

            gate.request();
            let ok_inputs = comp
                .build_inputs(500.0, &[0.0], [0.0; 10], [0.0; 11], 0.0, 101_325.0)
                .expect("inputs");
            let err_c = comp
                .propose(&state, 1.0, &ok_inputs, &mut next, &mut out, cx)
                .expect_err("cancelled context must refuse");
            assert!(matches!(err_c, VesselErrorV1::Cancelled(_)));
        });
    }
}
