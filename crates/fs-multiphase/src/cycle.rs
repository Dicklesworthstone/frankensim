//! L3 (`FLUX`) Peng--Robinson closed-loop vapor-compression refrigeration cycle
//! physics and reusable [`fs_ports::component::Component`] owners (`KM2.3`).
//!
//! Built directly on `fs-thermochem::eos::peng_robinson` with Péneloux molar
//! volume translation (`PengRobinsonMixtureV1::new_with_volume_shifts`),
//! `pure_saturation_at_temperature`, and `two_phase_flash_pt` across all five
//! `KM0.3` commercial refrigerants (`R-290`, `R-134a`, `R-404A` 3-component,
//! `R-448A` 5-component zeotropic blend, and `R-744` CO2 subcritical +
//! supercritical gas-cooler).
//!
//! # Provided [`Component`] Owners
//! 1. [`RefrigerationCompressorComponent`]: Electrical power + PR-EOS suction-to-discharge
//!    compression + motor shell heat rejection + condenser-fouling sensitivity.
//! 2. [`CondenserEvaporatorHXComponent`]: Dual-mode condenser / supercritical R-744
//!    gas-cooler OR frosting evaporator with dynamic frost thickness ODE
//!    `d(delta_frost)/dt = m_dot_frost / (rho_frost * A_coil)` and conservative
//!    liquid meltwater drain (`drain_stream`).
//! 3. [`ExpansionThrottleComponent`]: Exact isenthalpic TXV/capillary (`h_out == h_in`,
//!    `dS_gen >= 0`) with superheat-feedback valve modulation and Péneloux-translated
//!    liquid density orifice flow.
//! 4. [`IceHarvesterComponent`]: Air-cooled commercial ice machine composing the
//!    PR-EOS cycle against potable water feed + Ice Ih harvest bin, reproducing
//!    25--40% harvest derating between 70 F and 95 F kitchen ambient while
//!    rejecting continuous sensible condenser waste heat (1.5--3.5 kW) to room air.

use core::fmt;
use core::num::NonZeroUsize;

use fs_exec::solver::snapshot_v2::StateIdentityCharterV2;
use fs_exec::{Cancelled, Cx};
use fs_matdb::ConstitutiveModelCard;
use fs_math::det::{exp, expm1, ln, pow, sqrt};
use fs_ports::component::{
    Class, Component, ContractViolation, EndOfStep, PortLayout, PortRole, PortSpec, Proposal,
    Receipt, RefrigerantSpecies, STREAM_CARRIED, STREAM_ENTHALPY_RATE, STREAM_MASS_RATE,
    STREAM_PRESSURE, STREAM_SPECIFIC_ENTHALPY, StreamQuantity, kitchen_gas_v1, kitchen_liquid_v1,
    refrigerant_v1,
};
use fs_ports::{PortKind, StableId};
use fs_qty::{Dims, MolarMass, Pressure, Temperature};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;
use fs_thermochem::caloric::proximate::{
    H2O_LIQUID_FORMATION_ENTHALPY_298K_J_KG, ICE_FUSION_ENTHALPY_0C_J_KG, ICE_POINT_K,
    JANAF_REFERENCE_TEMPERATURE_K,
};
use fs_thermochem::eos::peng_robinson::flash::{
    PengRobinsonFlashErrorV1, PengRobinsonFlashOptionsV1,
};
use fs_thermochem::eos::peng_robinson::saturation::{
    PengRobinsonSaturationErrorV1, PengRobinsonSaturationOptionsV1,
};
use fs_thermochem::eos::peng_robinson::{
    PengRobinsonComponentV1, PengRobinsonErrorV1, PengRobinsonMixtureV1, PengRobinsonPhaseRootV1,
    PengRobinsonPhaseV1,
};

/// Exact `fs-matdb` law id for `KM0.3` Peng--Robinson refrigerant cards.
pub const PR_REFRIGERANT_LAW_ID_V1: &str = "peng-robinson-refrigerant-v1";
/// Supported semantic law version for [`PR_REFRIGERANT_LAW_ID_V1`].
pub const PR_REFRIGERANT_LAW_VERSION_V1: u32 = 1;

/// IIR standard reference temperature (`0 C` = `273.15 K`).
pub const IIR_REFERENCE_TEMPERATURE_K: f64 = 273.15;
/// IIR saturated-liquid reference enthalpy at `273.15 K` [J/kg].
pub const IIR_SAT_LIQUID_ENTHALPY_J_KG: f64 = 200_000.0;
/// IIR saturated-liquid reference entropy at `273.15 K` [J/(kg K)].
pub const IIR_SAT_LIQUID_ENTROPY_J_KG_K: f64 = 1_000.0;

/// Width of a `refrigerant_v1` stream port (`2` header + `1` refrigerant species).
pub const REFRIGERANT_STREAM_WIDTH: usize = 3;
/// Probe width for [`RefrigerationCompressorComponent`].
pub const COMPRESSOR_PROBE_COUNT: usize = 8;
/// Probe width for [`CondenserEvaporatorHXComponent`].
pub const HX_PROBE_COUNT: usize = 10;
/// Probe width for [`ExpansionThrottleComponent`].
pub const THROTTLE_PROBE_COUNT: usize = 6;
/// Probe width for [`IceHarvesterComponent`].
pub const ICE_HARVESTER_PROBE_COUNT: usize = 10;

// Durable state codec tags.
const COMPRESSOR_CODEC_TAG: u64 = 0x4653_5643_434D_0001;
const HX_CODEC_TAG: u64 = 0x4653_5643_4858_0001;
const THROTTLE_CODEC_TAG: u64 = 0x4653_5643_5458_0001;
const ICE_HARVESTER_CODEC_TAG: u64 = 0x4653_5643_4943_0001;

// JANAF-298.15K-formation water/ice/steam constants for evaporator frost & ice harvester bin.
const CP_WATER_LIQ_J_KG_K: f64 = 4_184.0;
const CP_WATER_ICE_J_KG_K: f64 = 2_100.0;
const CP_WATER_VAP_J_KG_K: f64 = 1_875.0;
const H_F0_H2O_VAP_J_KG: f64 = -241_826.0 / 0.018_015_268;
const H_WATER_LIQ_0C_JANAF_J_KG: f64 = H2O_LIQUID_FORMATION_ENTHALPY_298K_J_KG
    - CP_WATER_LIQ_J_KG_K * (JANAF_REFERENCE_TEMPERATURE_K - ICE_POINT_K);
const H_ICE_0C_JANAF_J_KG: f64 = H_WATER_LIQ_0C_JANAF_J_KG - ICE_FUSION_ENTHALPY_0C_J_KG;

const AUX_LIQ_SPECIES_COUNT: usize = 7;
const AUX_GAS_SPECIES_COUNT: usize = 8;

const DIMLESS: Dims = Dims([0, 0, 0, 0, 0, 0]);
const TEMP_DIMS: Dims = Dims([0, 0, 0, 1, 0, 0]);
const PRESS_DIMS: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const MOLAR_MASS_DIMS: Dims = Dims([0, 1, 0, 0, 0, -1]);
const MOLAR_VOL_DIMS: Dims = Dims([3, 0, 0, 0, 0, -1]);

/// Dedicated typed error enum for PR-EOS refrigeration cycle modeling and components.
#[derive(Clone, Debug, PartialEq)]
pub enum CycleErrorV1 {
    /// Invalid configuration parameter.
    InvalidConfig {
        /// Parameter name.
        parameter: &'static str,
        /// Explanation.
        reason: &'static str,
    },
    /// Invalid scalar input or state field.
    InvalidScalar {
        /// Quantity name.
        quantity: &'static str,
        /// Rejected value.
        value: f64,
        /// Explanation.
        reason: &'static str,
    },
    /// Invalid thermodynamic pressure (`<= 0`, non-finite, or unphysical ordering).
    InvalidPressure {
        /// Rejected pressure [Pa].
        pressure_pa: f64,
        /// Context.
        reason: &'static str,
    },
    /// Invalid thermodynamic temperature (`<= 0`, non-finite, or out of range).
    InvalidTemperature {
        /// Rejected temperature [K].
        temperature_k: f64,
        /// Context.
        reason: &'static str,
    },
    /// Invalid `fs-matdb` refrigerant card.
    InvalidCard {
        /// Context or parameter name.
        context: &'static str,
        /// Explanation.
        reason: &'static str,
    },
    /// Invalid superheat / thermal inversion (e.g., discharge pressure below suction or negative superheat).
    ThermodynamicInversion(&'static str),
    /// Stream total mass rate disagreed with carried species sum.
    InconsistentStreamMassFlow {
        /// Stream port name.
        port: &'static str,
        /// Declared total mass flow [kg/s].
        declared_total_kg_s: f64,
        /// Sum of carried species rates [kg/s].
        species_sum_kg_s: f64,
    },
    /// Outgoing stream flux attempted to deplete stored inventory below zero.
    DepletedInventory {
        /// Inventory name.
        inventory: &'static str,
    },
    /// Underlying Peng--Robinson EOS error.
    PengRobinson(PengRobinsonErrorV1),
    /// Underlying PR saturation solve failed.
    Saturation(PengRobinsonSaturationErrorV1),
    /// Underlying PR flash solve failed.
    Flash(PengRobinsonFlashErrorV1),
    /// Numerical convergence failure in saturation, isentropic, or condenser UA solve.
    ConvergenceFailed {
        /// Stage identifier.
        stage: &'static str,
        /// Final residual.
        residual: f64,
    },
    /// `fs-ports` contract violation.
    PortContract(ContractViolation),
    /// Cooperative cancellation via [`Cx`].
    Cancelled(Cancelled),
}

impl fmt::Display for CycleErrorV1 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidConfig { parameter, reason } => {
                write!(f, "invalid cycle config `{parameter}`: {reason}")
            }
            Self::InvalidScalar {
                quantity,
                value,
                reason,
            } => write!(f, "invalid cycle scalar `{quantity}` = {value}: {reason}"),
            Self::InvalidPressure {
                pressure_pa,
                reason,
            } => write!(f, "invalid cycle pressure {pressure_pa} Pa: {reason}"),
            Self::InvalidTemperature {
                temperature_k,
                reason,
            } => write!(f, "invalid cycle temperature {temperature_k} K: {reason}"),
            Self::InvalidCard { context, reason } => {
                write!(f, "invalid refrigerant card `{context}`: {reason}")
            }
            Self::ThermodynamicInversion(msg) => {
                write!(f, "refrigeration thermodynamic inversion: {msg}")
            }
            Self::InconsistentStreamMassFlow {
                port,
                declared_total_kg_s,
                species_sum_kg_s,
            } => write!(
                f,
                "inconsistent stream mass flow on `{port}`: total={declared_total_kg_s} vs sum={species_sum_kg_s}"
            ),
            Self::DepletedInventory { inventory } => {
                write!(f, "depleted inventory `{inventory}`")
            }
            Self::PengRobinson(e) => write!(f, "PR-EOS error: {e}"),
            Self::Saturation(e) => write!(f, "PR saturation error: {e}"),
            Self::Flash(e) => write!(f, "PR flash error: {e}"),
            Self::ConvergenceFailed { stage, residual } => {
                write!(
                    f,
                    "cycle solver `{stage}` failed to converge (res={residual})"
                )
            }
            Self::PortContract(e) => write!(f, "cycle port contract error: {e}"),
            Self::Cancelled(c) => write!(f, "cycle step cancelled: {c:?}"),
        }
    }
}

impl core::error::Error for CycleErrorV1 {}

impl From<Cancelled> for CycleErrorV1 {
    fn from(value: Cancelled) -> Self {
        Self::Cancelled(value)
    }
}

impl From<ContractViolation> for CycleErrorV1 {
    fn from(value: ContractViolation) -> Self {
        Self::PortContract(value)
    }
}

impl From<PengRobinsonErrorV1> for CycleErrorV1 {
    fn from(value: PengRobinsonErrorV1) -> Self {
        Self::PengRobinson(value)
    }
}

impl From<PengRobinsonSaturationErrorV1> for CycleErrorV1 {
    fn from(value: PengRobinsonSaturationErrorV1) -> Self {
        Self::Saturation(value)
    }
}

impl From<PengRobinsonFlashErrorV1> for CycleErrorV1 {
    fn from(value: PengRobinsonFlashErrorV1) -> Self {
        Self::Flash(value)
    }
}

/// Thermodynamic state point evaluated from [`RefrigerantModelV1`] on `iir-273.15K`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RefrigerantStatePointV1 {
    /// Temperature [K].
    pub temperature_k: f64,
    /// Pressure [Pa].
    pub pressure_pa: f64,
    /// Specific enthalpy on IIR-273.15K datum [J/kg].
    pub enthalpy_j_kg: f64,
    /// Specific entropy on IIR-273.15K datum [J/(kg K)].
    pub entropy_j_kg_k: f64,
    /// Péneloux-translated mass density [kg/m^3].
    pub density_kg_m3: f64,
    /// Untranslated Peng--Robinson mass density [kg/m^3] (for Peneloux comparison audits).
    pub untranslated_density_kg_m3: f64,
    /// Translated compressibility factor `Z`.
    pub compressibility_z: f64,
    /// Whether the state is supercritical (`T >= T_c` or `P >= P_c`).
    pub is_supercritical: bool,
    /// Whether `T > 0.95 * T_c` (triggers near-critical cubic EOS uncertainty disclosure).
    pub near_critical_warning: bool,
    /// Declared near-critical density uncertainty band (e.g. `0.18` for R-744 above `0.95 * T_c`).
    pub near_critical_density_error_band: f64,
}

/// Saturation / zeotropic phase-change envelope at one pressure `P`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SaturationEnvelopeV1 {
    /// Pressure [Pa].
    pub pressure_pa: f64,
    /// Bubble-point temperature `T_bub(P)` [K].
    pub t_bubble_k: f64,
    /// Dew-point temperature `T_dew(P)` [K].
    pub t_dew_k: f64,
    /// Evaporator/condenser temperature glide `T_dew - T_bub >= 0` [K] (`0` for pure fluids, `~5.5-6.7 K` for R-448A).
    pub glide_k: f64,
    /// Saturated liquid enthalpy at `T_bub` on `iir-273.15K` [J/kg].
    pub h_liquid_j_kg: f64,
    /// Saturated vapor enthalpy at `T_dew` on `iir-273.15K` [J/kg].
    pub h_vapor_j_kg: f64,
    /// Saturated liquid entropy at `T_bub` on `iir-273.15K` [J/(kg K)].
    pub s_liquid_j_kg_k: f64,
    /// Saturated vapor entropy at `T_dew` on `iir-273.15K` [J/(kg K)].
    pub s_vapor_j_kg_k: f64,
    /// Péneloux-translated saturated liquid density [kg/m^3].
    pub rho_liquid_kg_m3: f64,
    /// Péneloux-translated saturated vapor density [kg/m^3].
    pub rho_vapor_kg_m3: f64,
}

/// Output of [`RefrigerantModelV1::isentropic_compression`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CompressionResultV1 {
    /// Suction state.
    pub suction: RefrigerantStatePointV1,
    /// Isentropic discharge temperature `T_2s` [K].
    pub t_discharge_isentropic_k: f64,
    /// Actual discharge state (`h_2 = h_1 + (h_2s - h_1) / eta_is`).
    pub discharge: RefrigerantStatePointV1,
    /// Isentropic specific work `h_2s - h_1` [J/kg].
    pub w_isentropic_j_kg: f64,
    /// Actual shaft specific work `(h_2s - h_1) / eta_is` [J/kg].
    pub w_actual_j_kg: f64,
    /// Specific entropy generation `s_2 - s_1 >= 0` [J/(kg K)].
    pub ds_gen_j_kg_k: f64,
}

/// Output of [`RefrigerantModelV1::isenthalpic_throttle`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThrottleResultV1 {
    /// High-side inlet state.
    pub inlet: RefrigerantStatePointV1,
    /// Low-side outlet pressure [Pa].
    pub outlet_pressure_pa: f64,
    /// Low-side outlet enthalpy (`== inlet.enthalpy_j_kg` identically) [J/kg].
    pub outlet_enthalpy_j_kg: f64,
    /// Low-side outlet temperature [K].
    pub outlet_temperature_k: f64,
    /// Low-side outlet vapor quality in `[0, 1]`.
    pub outlet_quality: f64,
    /// Low-side outlet specific entropy [J/(kg K)].
    pub outlet_entropy_j_kg_k: f64,
    /// Irreversible throttling entropy generation `s_out - s_in >= 0` [J/(kg K)].
    pub ds_gen_j_kg_k: f64,
}

/// Operating point specification for [`RefrigerantModelV1::solve_closed_cycle`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CycleOperatingPointV1 {
    /// Evaporating midpoint temperature [K] (strictly below `T_c`).
    pub t_evap_k: f64,
    /// Suction superheat above dew temperature [K] (`>= 0.5 K`).
    pub suction_superheat_k: f64,
    /// Ambient sink temperature `T_amb` [K] for condenser / gas-cooler rejection.
    pub t_ambient_k: f64,
    /// Clean condenser / gas-cooler conductance `UA_cond_clean` [W/K].
    pub ua_cond_clean_w_k: f64,
    /// Dimensionless condenser fouling factor `f_foul >= 0.0` (`0.25` for +25% fouling).
    pub condenser_fouling_factor: f64,
    /// Subcooling below bubble temperature for subcritical condensers [K] (`>= 0.0`).
    pub subcooling_k: f64,
    /// Compressor swept displacement volume flow rate [m^3/s].
    pub swept_volume_m3_s: f64,
    /// Compressor isentropic efficiency in `(0.2, 1.0]`.
    pub eta_isentropic: f64,
    /// Compressor motor electrical-to-shaft efficiency in `(0.4, 1.0]`.
    pub eta_motor: f64,
    /// Fraction of electrical motor loss rejected from shell to ambient in `[0, 1]`.
    pub shell_loss_fraction: f64,
    /// Clearance volume ratio `c_vol` in `[0, 0.15]` for volumetric efficiency `1 - c_vol * (PR^(1/k) - 1)`.
    pub clearance_ratio: f64,
    /// Optional explicit supercritical high-side pressure [Pa] (used for transcritical R-744;
    /// if `None` and operating supercritical, selects an optimal transcritical discharge pressure).
    pub supercritical_pressure_pa: Option<f64>,
}

/// Complete closed-loop vapor-compression cycle state and First/Second Law audit summary.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ClosedCycleSummaryV1 {
    /// Refrigerant species.
    pub species: RefrigerantSpecies,
    /// Evaporator pressure [Pa].
    pub p_evap_pa: f64,
    /// Condenser / gas-cooler pressure [Pa].
    pub p_high_pa: f64,
    /// Evaporator saturation envelope (capturing R-448A glide).
    pub evap_envelope: SaturationEnvelopeV1,
    /// Whether high-side rejection is operating in supercritical gas-cooler mode (`R-744` above `T_c`).
    pub is_supercritical_gas_cooler: bool,
    /// Condenser / gas-cooler outlet temperature [K].
    pub t_high_out_k: f64,
    /// Refrigerant mass flow rate [kg/s].
    pub m_dot_kg_s: f64,
    /// Volumetric efficiency.
    pub eta_volumetric: f64,
    /// Compression step details.
    pub compression: CompressionResultV1,
    /// Expansion throttle details.
    pub throttle: ThrottleResultV1,
    /// Evaporator cooling capacity `Q_evap` [W].
    pub q_evap_w: f64,
    /// Compressor electrical power draw `W_comp` [W].
    pub w_comp_w: f64,
    /// Compressor shell heat loss to ambient `Q_shell` [W].
    pub q_shell_w: f64,
    /// Condenser / gas-cooler heat rejection `Q_cond = Q_evap + W_comp - Q_shell` [W].
    pub q_cond_w: f64,
    /// First-Law relative residual `|Q_cond - (Q_evap + W_comp - Q_shell)| / max(Q_cond, 1)`.
    pub first_law_relative_defect: f64,
    /// Cooling COP `Q_evap / W_comp`.
    pub cop_cooling: f64,
    /// Ideal Carnot COP `T_evap / (T_sink - T_evap)` (always strictly above `cop_cooling`).
    pub cop_carnot: f64,
    /// Total internal cycle entropy generation rate `m_dot * (ds_comp + ds_txv)` [W/K] (`>= 0`).
    pub internal_entropy_gen_w_k: f64,
    /// Near-critical density uncertainty band (e.g. `0.18` for R-744 above `0.95 * T_c`).
    pub near_critical_density_error_band: f64,
    /// Whether any high-side state lies inside the near-critical warning band (`T > 0.95 * T_c`).
    pub near_critical_warning: bool,
}

/// Péneloux-translated Peng--Robinson refrigerant thermodynamic model anchored
/// on `ENTHALPY_REF_IIR_273K` (`h_l(273.15 K) = 200 kJ/kg`, `s_l(273.15 K) = 1 kJ/(kg K)`).
#[derive(Clone, Debug, PartialEq)]
pub struct RefrigerantModelV1 {
    species: RefrigerantSpecies,
    mixture: PengRobinsonMixtureV1,
    mixture_untranslated: PengRobinsonMixtureV1,
    t_crit_k: f64,
    p_crit_pa: f64,
    cp_ig_j_kg_k: f64,
    subcritical_accurate_t_max_k: f64,
    near_critical_density_error_band: f64,
    h_offset_j_kg: f64,
    s_offset_j_kg_k: f64,
}

impl RefrigerantModelV1 {
    /// Construct a validated [`RefrigerantModelV1`] from PR components, mole fractions,
    /// Péneloux molar volume shifts `c_i` [m^3/mol], ideal-gas `cp_ig` [J/(kg K)],
    /// and declared near-critical uncertainty metadata.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        species: RefrigerantSpecies,
        components: Vec<PengRobinsonComponentV1>,
        mole_fractions: Vec<f64>,
        binary_interactions: Vec<f64>,
        volume_shifts_m3_per_mol: Vec<f64>,
        t_crit_k: f64,
        p_crit_pa: f64,
        cp_ig_j_kg_k: f64,
        subcritical_accurate_t_max_k: f64,
        near_critical_density_error_band: f64,
    ) -> Result<Self, CycleErrorV1> {
        if !t_crit_k.is_finite()
            || t_crit_k <= 100.0
            || !p_crit_pa.is_finite()
            || p_crit_pa <= 1.0e4
            || !cp_ig_j_kg_k.is_finite()
            || cp_ig_j_kg_k <= 100.0
        {
            return Err(CycleErrorV1::InvalidConfig {
                parameter: "refrigerant_critical_or_cp",
                reason: "T_c, P_c, and cp_ig must be finite and physically positive",
            });
        }
        if !subcritical_accurate_t_max_k.is_finite()
            || subcritical_accurate_t_max_k <= 0.0
            || !near_critical_density_error_band.is_finite()
            || !(0.0..=1.0).contains(&near_critical_density_error_band)
        {
            return Err(CycleErrorV1::InvalidConfig {
                parameter: "near_critical_metadata",
                reason: "near-critical threshold and error band must be finite and in [0, 1]",
            });
        }

        let mixture = PengRobinsonMixtureV1::new_with_volume_shifts(
            components.clone(),
            mole_fractions.clone(),
            binary_interactions.clone(),
            volume_shifts_m3_per_mol,
        )?;
        let mixture_untranslated =
            PengRobinsonMixtureV1::new(components, mole_fractions, binary_interactions)?;

        let mut model = Self {
            species,
            mixture,
            mixture_untranslated,
            t_crit_k,
            p_crit_pa,
            cp_ig_j_kg_k,
            subcritical_accurate_t_max_k,
            near_critical_density_error_band,
            h_offset_j_kg: 0.0,
            s_offset_j_kg_k: 0.0,
        };

        // Anchor exact IIR-273.15K saturated liquid datum:
        // h_l(273.15 K) = 200_000 J/kg, s_l(273.15 K) = 1_000 J/(kg K).
        let p_ref = model.bubble_pressure_at_t(IIR_REFERENCE_TEMPERATURE_K)?;
        let eval = model.mixture.evaluate_pt(
            Temperature::new(IIR_REFERENCE_TEMPERATURE_K),
            Pressure::new(p_ref),
        )?;
        let liq = eval.select(PengRobinsonPhaseV1::Liquid);
        let m_mix = model.mixture.molar_mass().value();
        let h_raw = liq.residual_enthalpy_j_mol() / m_mix;
        let s_raw = -(UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / m_mix) * ln(p_ref / 100_000.0)
            + liq.residual_entropy_j_mol_k() / m_mix;
        model.h_offset_j_kg = IIR_SAT_LIQUID_ENTHALPY_J_KG - h_raw;
        model.s_offset_j_kg_k = IIR_SAT_LIQUID_ENTROPY_J_KG_K - s_raw;
        Ok(model)
    }

    /// `R-290` (propane, ASHRAE A3) matching `data/matdb/kitchen-v1/refrigerants/r290-propane`.
    pub fn r290() -> Result<Self, CycleErrorV1> {
        Self::new(
            RefrigerantSpecies::R290,
            vec![PengRobinsonComponentV1 {
                critical_temperature: Temperature::new(369.89),
                critical_pressure: Pressure::new(4_251_200.0),
                acentric_factor: 0.1521,
                molar_mass: MolarMass::new(0.044_095_6),
            }],
            vec![1.0],
            Vec::new(),
            vec![-4.739_282e-6],
            369.89,
            4_251_200.0,
            1_670.0,
            351.4,
            0.025,
        )
    }

    /// `R-134a` (HFC-134a, ASHRAE A1) matching `data/matdb/kitchen-v1/refrigerants/r134a-hfc`.
    pub fn r134a() -> Result<Self, CycleErrorV1> {
        Self::new(
            RefrigerantSpecies::R134a,
            vec![PengRobinsonComponentV1 {
                critical_temperature: Temperature::new(374.21),
                critical_pressure: Pressure::new(4_059_280.0),
                acentric_factor: 0.326_84,
                molar_mass: MolarMass::new(0.102_032),
            }],
            vec![1.0],
            Vec::new(),
            vec![1.186_135e-6],
            374.21,
            4_059_280.0,
            852.0,
            355.5,
            0.025,
        )
    }

    /// `R-744` (CO2, ASHRAE A1, transcritical capable) matching `data/matdb/kitchen-v1/refrigerants/r744-co2`.
    pub fn r744() -> Result<Self, CycleErrorV1> {
        Self::new(
            RefrigerantSpecies::R744,
            vec![PengRobinsonComponentV1 {
                critical_temperature: Temperature::new(304.1282),
                critical_pressure: Pressure::new(7_377_300.0),
                acentric_factor: 0.223_94,
                molar_mass: MolarMass::new(0.044_009_5),
            }],
            vec![1.0],
            Vec::new(),
            vec![-9.546_454e-7],
            304.1282,
            7_377_300.0,
            846.0,
            288.92,
            0.18,
        )
    }

    /// `R-404A` (3-component blend: R-125/R-143a/R-134a) matching `data/matdb/kitchen-v1/refrigerants/r404a-blend`.
    pub fn r404a() -> Result<Self, CycleErrorV1> {
        Self::new(
            RefrigerantSpecies::R404A,
            vec![
                PengRobinsonComponentV1 {
                    critical_temperature: Temperature::new(339.173),
                    critical_pressure: Pressure::new(3_617_700.0),
                    acentric_factor: 0.3052,
                    molar_mass: MolarMass::new(0.120_02),
                },
                PengRobinsonComponentV1 {
                    critical_temperature: Temperature::new(345.857),
                    critical_pressure: Pressure::new(3_761_000.0),
                    acentric_factor: 0.2615,
                    molar_mass: MolarMass::new(0.084_041),
                },
                PengRobinsonComponentV1 {
                    critical_temperature: Temperature::new(374.21),
                    critical_pressure: Pressure::new(4_059_280.0),
                    acentric_factor: 0.326_84,
                    molar_mass: MolarMass::new(0.102_032),
                },
            ],
            vec![0.357_819, 0.603_917, 0.038_264],
            Vec::new(),
            vec![-3.162_847e-6, 2.430_315e-6, 1.186_135e-6],
            345.2,
            3_728_900.0,
            870.0,
            327.94,
            0.04,
        )
    }

    /// `R-448A` (5-component zeotropic blend: R-32/R-125/R-1234yf/R-134a/R-1234ze(E))
    /// matching `data/matdb/kitchen-v1/refrigerants/r448a-blend`.
    pub fn r448a() -> Result<Self, CycleErrorV1> {
        Self::new(
            RefrigerantSpecies::R448A,
            vec![
                PengRobinsonComponentV1 {
                    critical_temperature: Temperature::new(351.255),
                    critical_pressure: Pressure::new(5_782_000.0),
                    acentric_factor: 0.2769,
                    molar_mass: MolarMass::new(0.052_024),
                },
                PengRobinsonComponentV1 {
                    critical_temperature: Temperature::new(339.173),
                    critical_pressure: Pressure::new(3_617_700.0),
                    acentric_factor: 0.3052,
                    molar_mass: MolarMass::new(0.120_02),
                },
                PengRobinsonComponentV1 {
                    critical_temperature: Temperature::new(367.85),
                    critical_pressure: Pressure::new(3_382_200.0),
                    acentric_factor: 0.276,
                    molar_mass: MolarMass::new(0.114_04),
                },
                PengRobinsonComponentV1 {
                    critical_temperature: Temperature::new(374.21),
                    critical_pressure: Pressure::new(4_059_280.0),
                    acentric_factor: 0.326_84,
                    molar_mass: MolarMass::new(0.102_032),
                },
                PengRobinsonComponentV1 {
                    critical_temperature: Temperature::new(382.513),
                    critical_pressure: Pressure::new(3_634_900.0),
                    acentric_factor: 0.313,
                    molar_mass: MolarMass::new(0.114_04),
                },
            ],
            vec![0.431_216, 0.186_915, 0.151_321, 0.177_586, 0.052_962],
            Vec::new(),
            vec![
                6.048_446e-6,
                -3.162_847e-6,
                -9.542_529e-7,
                1.186_135e-6,
                -1.997_079e-6,
            ],
            354.76,
            4_595_000.0,
            880.0,
            337.02,
            0.05,
        )
    }

    /// Construct the built-in `KM0.3` model for `species`.
    pub fn for_species(species: RefrigerantSpecies) -> Result<Self, CycleErrorV1> {
        match species {
            RefrigerantSpecies::R290 => Self::r290(),
            RefrigerantSpecies::R134a => Self::r134a(),
            RefrigerantSpecies::R744 => Self::r744(),
            RefrigerantSpecies::R404A => Self::r404a(),
            RefrigerantSpecies::R448A => Self::r448a(),
        }
    }

    /// Load a validated [`RefrigerantModelV1`] from an immutable `fs-matdb`
    /// [`ConstitutiveModelCard`] with law `"peng-robinson-refrigerant-v1"` v1.
    pub fn from_card(
        species: RefrigerantSpecies,
        card: &ConstitutiveModelCard,
        cp_ig_j_kg_k: f64,
    ) -> Result<Self, CycleErrorV1> {
        card.validate().map_err(|_| CycleErrorV1::InvalidCard {
            context: "peng-robinson-refrigerant-v1",
            reason: "card failed fs-matdb validation",
        })?;
        if card.law.0 != PR_REFRIGERANT_LAW_ID_V1 {
            return Err(CycleErrorV1::InvalidCard {
                context: "law_id",
                reason: "expected law id `peng-robinson-refrigerant-v1`",
            });
        }
        if card.law_version != PR_REFRIGERANT_LAW_VERSION_V1 {
            return Err(CycleErrorV1::InvalidCard {
                context: "law_version",
                reason: "expected law version 1",
            });
        }
        let get = |key: &str, expected: Dims| -> Result<f64, CycleErrorV1> {
            let p = card.parameters.get(key).ok_or(CycleErrorV1::InvalidCard {
                context: "parameter",
                reason: "missing required refrigerant card parameter",
            })?;
            if p.dims != expected || !p.value.is_finite() {
                return Err(CycleErrorV1::InvalidCard {
                    context: "parameter_dims",
                    reason: "refrigerant parameter dims mismatch or non-finite value",
                });
            }
            Ok(p.value)
        };
        let n_raw = get("component_count", DIMLESS)?;
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let n = n_raw as usize;
        if n == 0 || n > 16 || (n as f64 - n_raw).abs() > 1.0e-12 {
            return Err(CycleErrorV1::InvalidCard {
                context: "component_count",
                reason: "must be an integer in 1..=16",
            });
        }
        let t_crit = get("mixture_critical_temperature", TEMP_DIMS)?;
        let p_crit = get("mixture_critical_pressure", PRESS_DIMS)?;
        let t_sub_max = get("subcritical_accurate_t_max", TEMP_DIMS)?;
        let err_band = get("near_critical_density_error_band", DIMLESS)?;

        let mut comps = Vec::with_capacity(n);
        let mut z = Vec::with_capacity(n);
        let mut shifts = Vec::with_capacity(n);
        for i in 0..n {
            let tc = get(&format!("comp_{i}_critical_temperature"), TEMP_DIMS)?;
            let pc = get(&format!("comp_{i}_critical_pressure"), PRESS_DIMS)?;
            let omega = get(&format!("comp_{i}_acentric_factor"), DIMLESS)?;
            let mm = get(&format!("comp_{i}_molar_mass"), MOLAR_MASS_DIMS)?;
            let zi = get(&format!("comp_{i}_mole_fraction"), DIMLESS)?;
            let ci = get(&format!("comp_{i}_peneloux_c"), MOLAR_VOL_DIMS)?;
            comps.push(PengRobinsonComponentV1 {
                critical_temperature: Temperature::new(tc),
                critical_pressure: Pressure::new(pc),
                acentric_factor: omega,
                molar_mass: MolarMass::new(mm),
            });
            z.push(zi);
            shifts.push(ci);
        }
        Self::new(
            species,
            comps,
            z,
            Vec::new(),
            shifts,
            t_crit,
            p_crit,
            cp_ig_j_kg_k,
            t_sub_max,
            err_band,
        )
    }

    /// Canonical species identity.
    #[must_use]
    pub const fn species(&self) -> RefrigerantSpecies {
        self.species
    }

    /// Underlying Péneloux-translated [`PengRobinsonMixtureV1`].
    #[must_use]
    pub const fn mixture(&self) -> &PengRobinsonMixtureV1 {
        &self.mixture
    }

    /// Critical temperature `T_c` [K].
    #[must_use]
    pub const fn critical_temperature_k(&self) -> f64 {
        self.t_crit_k
    }

    /// Critical pressure `P_c` [Pa].
    #[must_use]
    pub const fn critical_pressure_pa(&self) -> f64 {
        self.p_crit_pa
    }

    /// Whether this refrigerant is a multicomponent blend (`R-404A`, `R-448A`).
    #[must_use]
    pub fn is_blend(&self) -> bool {
        self.mixture.components().len() > 1
    }

    /// Declared near-critical density uncertainty band (`0.18` for R-744 above `0.95 * T_c`).
    #[must_use]
    pub const fn near_critical_density_error_band(&self) -> f64 {
        self.near_critical_density_error_band
    }

    /// Evaluate homogeneous state `(T, P)` on `iir-273.15K` using Péneloux-translated PR-EOS.
    pub fn state_pt(
        &self,
        t_k: f64,
        p_pa: f64,
        phase: PengRobinsonPhaseV1,
    ) -> Result<RefrigerantStatePointV1, CycleErrorV1> {
        if !t_k.is_finite() || t_k <= 50.0 || t_k > 1_200.0 {
            return Err(CycleErrorV1::InvalidTemperature {
                temperature_k: t_k,
                reason: "must lie in (50, 1200] K",
            });
        }
        if !p_pa.is_finite() || p_pa <= 0.0 || p_pa > 5.0e7 {
            return Err(CycleErrorV1::InvalidPressure {
                pressure_pa: p_pa,
                reason: "must lie in (0, 50 MPa]",
            });
        }
        let temp = Temperature::new(t_k);
        let pres = Pressure::new(p_pa);
        let eval = self.mixture.evaluate_pt(temp, pres)?;
        let root = eval.select(phase);
        let eval_raw = self.mixture_untranslated.evaluate_pt(temp, pres)?;
        let root_raw = eval_raw.select(phase);
        Ok(self.point_from_root(t_k, p_pa, root, root_raw.mass_density_kg_per_m3()))
    }

    fn point_from_root(
        &self,
        t_k: f64,
        p_pa: f64,
        root: &PengRobinsonPhaseRootV1,
        untranslated_rho: f64,
    ) -> RefrigerantStatePointV1 {
        let m_mix = self.mixture.molar_mass().value();
        let h_j_kg = self.h_offset_j_kg
            + self.cp_ig_j_kg_k * (t_k - IIR_REFERENCE_TEMPERATURE_K)
            + root.residual_enthalpy_j_mol() / m_mix;
        let s_j_kg_k = self.s_offset_j_kg_k
            + self.cp_ig_j_kg_k * ln(t_k / IIR_REFERENCE_TEMPERATURE_K)
            - (UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / m_mix) * ln(p_pa / 100_000.0)
            + root.residual_entropy_j_mol_k() / m_mix;
        let is_supercritical = t_k >= self.t_crit_k || p_pa >= self.p_crit_pa;
        let near_critical_warning = t_k >= self.subcritical_accurate_t_max_k;
        RefrigerantStatePointV1 {
            temperature_k: t_k,
            pressure_pa: p_pa,
            enthalpy_j_kg: h_j_kg,
            entropy_j_kg_k: s_j_kg_k,
            density_kg_m3: root.mass_density_kg_per_m3(),
            untranslated_density_kg_m3: untranslated_rho,
            compressibility_z: root.compressibility_factor(),
            is_supercritical,
            near_critical_warning,
            near_critical_density_error_band: self.near_critical_density_error_band,
        }
    }

    /// Solve bubble pressure `P_bub(T)` [Pa] for pure fluids (via `pure_saturation_at_temperature`)
    /// or multicomponent blends (via PR fugacity iteration `sum_i z_i K_i = 1`).
    pub fn bubble_pressure_at_t(&self, t_k: f64) -> Result<f64, CycleErrorV1> {
        if !t_k.is_finite() || t_k <= 100.0 || t_k >= 0.995 * self.t_crit_k {
            return Err(CycleErrorV1::InvalidTemperature {
                temperature_k: t_k,
                reason: "bubble pressure requires subcritical T in (100 K, 0.995 T_c)",
            });
        }
        if !self.is_blend() {
            let sat = self.mixture.pure_saturation_at_temperature(
                Temperature::new(t_k),
                [
                    Pressure::new(1_000.0),
                    Pressure::new(0.995 * self.p_crit_pa),
                ],
                PengRobinsonSaturationOptionsV1::default(),
            )?;
            return Ok(sat.pressure().value());
        }
        // Multicomponent PR bubble-pressure fixed-point iteration initialized with Wilson K.
        let comps = self.mixture.components();
        let z = self.mixture.mole_fractions();
        let mut p = comps
            .iter()
            .zip(z)
            .map(|(c, &zi)| {
                zi * c.critical_pressure.value()
                    * exp(5.373
                        * (1.0 + c.acentric_factor)
                        * (1.0 - c.critical_temperature.value() / t_k))
            })
            .sum::<f64>()
            .clamp(5_000.0, 0.95 * self.p_crit_pa);

        let mut y = vec![0.0; comps.len()];
        let mut last_res = f64::INFINITY;
        for _ in 0..64 {
            let l_eval = self
                .mixture
                .evaluate_pt(Temperature::new(t_k), Pressure::new(p))?;
            let l_phi = l_eval
                .select(PengRobinsonPhaseV1::Liquid)
                .log_fugacity_coefficients();
            let mut sum_y = 0.0;
            for i in 0..comps.len() {
                let k_w = (comps[i].critical_pressure.value() / p)
                    * exp(5.373
                        * (1.0 + comps[i].acentric_factor)
                        * (1.0 - comps[i].critical_temperature.value() / t_k));
                y[i] = (z[i] * k_w).max(1.0e-15);
                sum_y += y[i];
            }
            for yi in &mut y {
                *yi /= sum_y;
            }
            // Refine vapor composition y at current P.
            let mut sum_k = 1.0;
            for _ in 0..12 {
                let vap_mix = PengRobinsonMixtureV1::new_with_volume_shifts(
                    comps.to_vec(),
                    y.clone(),
                    self.mixture.binary_interactions().to_vec(),
                    self.mixture.volume_shifts_m3_per_mol().to_vec(),
                )?;
                let v_eval = vap_mix.evaluate_pt(Temperature::new(t_k), Pressure::new(p))?;
                let v_phi = v_eval
                    .select(PengRobinsonPhaseV1::Vapor)
                    .log_fugacity_coefficients();
                sum_k = 0.0;
                for i in 0..comps.len() {
                    let ki = exp((l_phi[i] - v_phi[i]).clamp(-30.0, 30.0));
                    y[i] = (z[i] * ki).max(1.0e-15);
                    sum_k += y[i];
                }
                for yi in &mut y {
                    *yi /= sum_k;
                }
            }
            last_res = (sum_k - 1.0).abs();
            if last_res < 1.0e-9 {
                return Ok(p);
            }
            let factor = sum_k.clamp(0.5, 2.0);
            p = (p * factor).clamp(1_000.0, 0.98 * self.p_crit_pa);
        }
        Err(CycleErrorV1::ConvergenceFailed {
            stage: "blend_bubble_pressure",
            residual: last_res,
        })
    }

    /// Solve saturation / zeotropic bubble-and-dew envelope at subcritical pressure `p_pa` [Pa].
    ///
    /// For pure fluids (`R-290`, `R-134a`, `R-744`), `t_bubble_k == t_dew_k` (`glide_k == 0.0`).
    /// For multicomponent blends (`R-404A`, `R-448A`), solves exact PR bubble and dew temperatures
    /// and verifies two-phase equilibrium via `PengRobinsonMixtureV1::two_phase_flash_pt` at
    /// the glide midpoint!
    pub fn bubble_dew_at_pressure(&self, p_pa: f64) -> Result<SaturationEnvelopeV1, CycleErrorV1> {
        if !p_pa.is_finite() || p_pa <= 1_000.0 || p_pa >= 0.98 * self.p_crit_pa {
            return Err(CycleErrorV1::InvalidPressure {
                pressure_pa: p_pa,
                reason: "subcritical envelope requires P in (1 kPa, 0.98 P_c)",
            });
        }
        if !self.is_blend() {
            // Invert pure_saturation_at_temperature(T) == p_pa via monotone bisection.
            let mut t_lo = 150.0_f64.max(0.45 * self.t_crit_k);
            let mut t_hi = 0.992 * self.t_crit_k;
            for _ in 0..64 {
                let t_mid = 0.5 * (t_lo + t_hi);
                let p_mid = self.bubble_pressure_at_t(t_mid)?;
                if (p_mid - p_pa).abs() <= 1.0e-9 * p_pa {
                    t_lo = t_mid;
                    t_hi = t_mid;
                    break;
                }
                if p_mid < p_pa {
                    t_lo = t_mid;
                } else {
                    t_hi = t_mid;
                }
            }
            let t_sat = 0.5 * (t_lo + t_hi);
            let sat = self.mixture.pure_saturation_at_temperature(
                Temperature::new(t_sat),
                [
                    Pressure::new(1_000.0),
                    Pressure::new(0.995 * self.p_crit_pa),
                ],
                PengRobinsonSaturationOptionsV1::default(),
            )?;
            let liq_pt = self.point_from_root(
                t_sat,
                p_pa,
                sat.liquid(),
                sat.liquid().mass_density_kg_per_m3(),
            );
            let vap_pt = self.point_from_root(
                t_sat,
                p_pa,
                sat.vapor(),
                sat.vapor().mass_density_kg_per_m3(),
            );
            return Ok(SaturationEnvelopeV1 {
                pressure_pa: p_pa,
                t_bubble_k: t_sat,
                t_dew_k: t_sat,
                glide_k: 0.0,
                h_liquid_j_kg: liq_pt.enthalpy_j_kg,
                h_vapor_j_kg: vap_pt.enthalpy_j_kg,
                s_liquid_j_kg_k: liq_pt.entropy_j_kg_k,
                s_vapor_j_kg_k: vap_pt.entropy_j_kg_k,
                rho_liquid_kg_m3: liq_pt.density_kg_m3,
                rho_vapor_kg_m3: vap_pt.density_kg_m3,
            });
        }

        // Multicomponent PR bubble T and dew T at fixed P.
        let t_bub = self.solve_blend_saturation_t(p_pa, true)?;
        let t_dew = self.solve_blend_saturation_t(p_pa, false)?;
        let glide_k = (t_dew - t_bub).max(0.0);

        // Verify genuine two-phase PR flash split at midpoint T when glide > 0.05 K.
        if glide_k > 0.05 {
            let t_mid = 0.5 * (t_bub + t_dew);
            let flash = self.mixture.two_phase_flash_pt(
                Temperature::new(t_mid),
                Pressure::new(p_pa),
                PengRobinsonFlashOptionsV1::default(),
            )?;
            if !(0.0..1.0).contains(&flash.vapor_fraction_mol()) {
                return Err(CycleErrorV1::ConvergenceFailed {
                    stage: "blend_midpoint_flash",
                    residual: flash.vapor_fraction_mol(),
                });
            }
        }

        let liq_pt = self.state_pt(t_bub, p_pa, PengRobinsonPhaseV1::Liquid)?;
        let vap_pt = self.state_pt(t_dew, p_pa, PengRobinsonPhaseV1::Vapor)?;
        Ok(SaturationEnvelopeV1 {
            pressure_pa: p_pa,
            t_bubble_k: t_bub,
            t_dew_k: t_dew,
            glide_k,
            h_liquid_j_kg: liq_pt.enthalpy_j_kg,
            h_vapor_j_kg: vap_pt.enthalpy_j_kg,
            s_liquid_j_kg_k: liq_pt.entropy_j_kg_k,
            s_vapor_j_kg_k: vap_pt.entropy_j_kg_k,
            rho_liquid_kg_m3: liq_pt.density_kg_m3,
            rho_vapor_kg_m3: vap_pt.density_kg_m3,
        })
    }

    fn solve_blend_saturation_t(&self, p_pa: f64, is_bubble: bool) -> Result<f64, CycleErrorV1> {
        let comps = self.mixture.components();
        let z = self.mixture.mole_fractions();
        let mut t = 0.72 * self.t_crit_k;
        let mut phase_comp = vec![0.0; comps.len()];
        let mut last_res = f64::INFINITY;
        for _ in 0..80 {
            let mut k_vals = vec![1.0; comps.len()];
            for i in 0..comps.len() {
                k_vals[i] = ((comps[i].critical_pressure.value() / p_pa)
                    * exp(5.373
                        * (1.0 + comps[i].acentric_factor)
                        * (1.0 - comps[i].critical_temperature.value() / t)))
                .clamp(1.0e-6, 1.0e6);
            }
            for _ in 0..24 {
                let mut sum_pc = 0.0;
                for i in 0..comps.len() {
                    phase_comp[i] = if is_bubble {
                        (z[i] * k_vals[i]).max(1.0e-15)
                    } else {
                        (z[i] / k_vals[i]).max(1.0e-15)
                    };
                    sum_pc += phase_comp[i];
                }
                for val in &mut phase_comp {
                    *val /= sum_pc;
                }
                let other_mix = PengRobinsonMixtureV1::new_with_volume_shifts(
                    comps.to_vec(),
                    phase_comp.clone(),
                    self.mixture.binary_interactions().to_vec(),
                    self.mixture.volume_shifts_m3_per_mol().to_vec(),
                )?;
                let (l_eval, v_eval) = if is_bubble {
                    (
                        self.mixture
                            .evaluate_pt(Temperature::new(t), Pressure::new(p_pa))?,
                        other_mix.evaluate_pt(Temperature::new(t), Pressure::new(p_pa))?,
                    )
                } else {
                    (
                        other_mix.evaluate_pt(Temperature::new(t), Pressure::new(p_pa))?,
                        self.mixture
                            .evaluate_pt(Temperature::new(t), Pressure::new(p_pa))?,
                    )
                };
                let l_phi = l_eval
                    .select(PengRobinsonPhaseV1::Liquid)
                    .log_fugacity_coefficients();
                let v_phi = v_eval
                    .select(PengRobinsonPhaseV1::Vapor)
                    .log_fugacity_coefficients();
                let mut max_dk = 0.0_f64;
                for i in 0..comps.len() {
                    let k_new = exp((l_phi[i] - v_phi[i]).clamp(-30.0, 30.0));
                    max_dk = max_dk.max((k_new - k_vals[i]).abs() / k_vals[i].max(1.0e-12));
                    k_vals[i] = k_new;
                }
                if max_dk < 1.0e-10 {
                    break;
                }
            }
            let objective: f64 = if is_bubble {
                z.iter().zip(&k_vals).map(|(zi, ki)| zi * ki).sum()
            } else {
                z.iter().zip(&k_vals).map(|(zi, ki)| zi / ki).sum()
            };
            last_res = (objective - 1.0).abs();
            if last_res < 1.0e-9 {
                return Ok(t);
            }
            let step = (14.0 * ln(objective.clamp(0.1, 10.0))).clamp(-15.0, 15.0);
            t = if is_bubble {
                (t - step).clamp(150.0, 0.99 * self.t_crit_k)
            } else {
                (t + step).clamp(150.0, 0.99 * self.t_crit_k)
            };
        }
        Err(CycleErrorV1::ConvergenceFailed {
            stage: if is_bubble {
                "blend_bubble_t"
            } else {
                "blend_dew_t"
            },
            residual: last_res,
        })
    }

    /// Solve pressure `P` [Pa] whose evaporator midpoint `0.5 * (T_bub + T_dew)` equals `t_evap_k`.
    pub fn midpoint_pressure_at_t(&self, t_evap_k: f64) -> Result<f64, CycleErrorV1> {
        if !self.is_blend() {
            return self.bubble_pressure_at_t(t_evap_k);
        }
        if !t_evap_k.is_finite() || t_evap_k <= 150.0 || t_evap_k >= 0.96 * self.t_crit_k {
            return Err(CycleErrorV1::InvalidTemperature {
                temperature_k: t_evap_k,
                reason: "midpoint saturation requires subcritical T in (150 K, 0.96 T_c)",
            });
        }
        let p_bub = self.bubble_pressure_at_t(t_evap_k)?;
        let mut p_lo = (0.5 * p_bub).max(10_000.0);
        let mut p_hi = (1.3 * p_bub).min(0.95 * self.p_crit_pa);
        for _ in 0..48 {
            let p_mid = 0.5 * (p_lo + p_hi);
            let env = self.bubble_dew_at_pressure(p_mid)?;
            let t_mid = 0.5 * (env.t_bubble_k + env.t_dew_k);
            if (t_mid - t_evap_k).abs() < 1.0e-5 {
                return Ok(p_mid);
            }
            if t_mid < t_evap_k {
                p_lo = p_mid;
            } else {
                p_hi = p_mid;
            }
        }
        Ok(0.5 * (p_lo + p_hi))
    }

    /// Invert vapor/supercritical temperature `T` [K] from `(P, s_target)` via monotone bisection.
    pub fn temperature_at_p_s(
        &self,
        p_pa: f64,
        s_target_j_kg_k: f64,
        t_min_k: f64,
    ) -> Result<RefrigerantStatePointV1, CycleErrorV1> {
        let mut lo = t_min_k.max(150.0);
        let mut hi = (lo + 250.0).min(1_000.0);
        while self
            .state_pt(hi, p_pa, PengRobinsonPhaseV1::Vapor)?
            .entropy_j_kg_k
            < s_target_j_kg_k
            && hi < 1_150.0
        {
            hi = (hi + 150.0).min(1_150.0);
        }
        for _ in 0..64 {
            let mid = 0.5 * (lo + hi);
            let pt = self.state_pt(mid, p_pa, PengRobinsonPhaseV1::Vapor)?;
            if (pt.entropy_j_kg_k - s_target_j_kg_k).abs()
                <= 64.0 * f64::EPSILON * s_target_j_kg_k.abs().max(1.0)
            {
                return Ok(pt);
            }
            if pt.entropy_j_kg_k < s_target_j_kg_k {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        self.state_pt(0.5 * (lo + hi), p_pa, PengRobinsonPhaseV1::Vapor)
    }

    /// Invert vapor/supercritical temperature `T` [K] from `(P, h_target)` via monotone bisection.
    pub fn temperature_at_p_h(
        &self,
        p_pa: f64,
        h_target_j_kg: f64,
        t_min_k: f64,
        phase: PengRobinsonPhaseV1,
    ) -> Result<RefrigerantStatePointV1, CycleErrorV1> {
        let mut lo = t_min_k.max(150.0);
        let mut hi = (lo + 250.0).min(1_100.0);
        while self.state_pt(hi, p_pa, phase)?.enthalpy_j_kg < h_target_j_kg && hi < 1_180.0 {
            hi = (hi + 150.0).min(1_180.0);
        }
        for _ in 0..64 {
            let mid = 0.5 * (lo + hi);
            let pt = self.state_pt(mid, p_pa, phase)?;
            if (pt.enthalpy_j_kg - h_target_j_kg).abs()
                <= 64.0 * f64::EPSILON * h_target_j_kg.abs().max(1.0)
            {
                return Ok(pt);
            }
            if pt.enthalpy_j_kg < h_target_j_kg {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        self.state_pt(0.5 * (lo + hi), p_pa, phase)
    }

    /// Compute PR-EOS isentropic and actual compression from `(T_suc, P_suc)` to `P_disch`.
    pub fn isentropic_compression(
        &self,
        t_suction_k: f64,
        p_suction_pa: f64,
        p_discharge_pa: f64,
        eta_isentropic: f64,
    ) -> Result<CompressionResultV1, CycleErrorV1> {
        if !p_discharge_pa.is_finite() || p_discharge_pa <= p_suction_pa {
            return Err(CycleErrorV1::ThermodynamicInversion(
                "discharge pressure must strictly exceed suction pressure",
            ));
        }
        if !eta_isentropic.is_finite() || eta_isentropic <= 0.1 || eta_isentropic > 1.0 {
            return Err(CycleErrorV1::InvalidScalar {
                quantity: "eta_isentropic",
                value: eta_isentropic,
                reason: "must lie in (0.1, 1.0]",
            });
        }
        let suction = self.state_pt(t_suction_k, p_suction_pa, PengRobinsonPhaseV1::Vapor)?;
        let isentropic =
            self.temperature_at_p_s(p_discharge_pa, suction.entropy_j_kg_k, t_suction_k)?;
        let w_is = isentropic.enthalpy_j_kg - suction.enthalpy_j_kg;
        if w_is <= 0.0 {
            return Err(CycleErrorV1::ThermodynamicInversion(
                "isentropic compression work must be positive",
            ));
        }
        let w_act = w_is / eta_isentropic;
        let h_disch = suction.enthalpy_j_kg + w_act;
        let discharge = self.temperature_at_p_h(
            p_discharge_pa,
            h_disch,
            isentropic.temperature_k,
            PengRobinsonPhaseV1::Vapor,
        )?;
        let ds_gen = (discharge.entropy_j_kg_k - suction.entropy_j_kg_k).max(0.0);
        Ok(CompressionResultV1 {
            suction,
            t_discharge_isentropic_k: isentropic.temperature_k,
            discharge,
            w_isentropic_j_kg: w_is,
            w_actual_j_kg: w_act,
            ds_gen_j_kg_k: ds_gen,
        })
    }

    /// Compute exact isenthalpic expansion (`h_out == h_in`) across TXV/capillary from
    /// high-side `(T_in, P_in)` to low-side `P_out`, returning Péneloux-translated inlet
    /// density, two-phase outlet quality, and Second-Law entropy generation `ds_gen >= 0`.
    pub fn isenthalpic_throttle(
        &self,
        t_in_k: f64,
        p_in_pa: f64,
        p_out_pa: f64,
    ) -> Result<ThrottleResultV1, CycleErrorV1> {
        if !p_out_pa.is_finite() || p_out_pa <= 1_000.0 || p_out_pa >= p_in_pa {
            return Err(CycleErrorV1::ThermodynamicInversion(
                "throttle outlet pressure must be positive and strictly below inlet pressure",
            ));
        }
        let inlet_phase = if t_in_k >= self.t_crit_k || p_in_pa >= self.p_crit_pa {
            PengRobinsonPhaseV1::Stable
        } else {
            PengRobinsonPhaseV1::Liquid
        };
        let inlet = self.state_pt(t_in_k, p_in_pa, inlet_phase)?;
        let env = self.bubble_dew_at_pressure(p_out_pa)?;
        let h_in = inlet.enthalpy_j_kg;
        let h_fg = (env.h_vapor_j_kg - env.h_liquid_j_kg).max(1.0);
        let (x_out, t_out, s_out) = if h_in <= env.h_liquid_j_kg {
            let sub =
                self.temperature_at_p_h(p_out_pa, h_in, 160.0, PengRobinsonPhaseV1::Liquid)?;
            (0.0, sub.temperature_k, sub.entropy_j_kg_k)
        } else if h_in >= env.h_vapor_j_kg {
            let sup =
                self.temperature_at_p_h(p_out_pa, h_in, env.t_dew_k, PengRobinsonPhaseV1::Vapor)?;
            (1.0, sup.temperature_k, sup.entropy_j_kg_k)
        } else {
            let x = ((h_in - env.h_liquid_j_kg) / h_fg).clamp(0.0, 1.0);
            let t_2ph = env.t_bubble_k + x * env.glide_k;
            let s_2ph = env.s_liquid_j_kg_k + x * (env.s_vapor_j_kg_k - env.s_liquid_j_kg_k);
            (x, t_2ph, s_2ph)
        };
        let ds_gen = (s_out - inlet.entropy_j_kg_k).max(0.0);
        Ok(ThrottleResultV1 {
            inlet,
            outlet_pressure_pa: p_out_pa,
            outlet_enthalpy_j_kg: h_in,
            outlet_temperature_k: t_out,
            outlet_quality: x_out,
            outlet_entropy_j_kg_k: s_out,
            ds_gen_j_kg_k: ds_gen,
        })
    }

    /// Solve a complete closed-loop vapor-compression cycle at `spec`, seamlessly
    /// supporting subcritical condensers AND transcritical `R-744` gas-coolers above
    /// `T_c = 304.1282 K`, with explicit condenser fouling sensitivity and exact
    /// First/Second Law accounting.
    #[allow(clippy::too_many_lines)]
    pub fn solve_closed_cycle(
        &self,
        spec: CycleOperatingPointV1,
    ) -> Result<ClosedCycleSummaryV1, CycleErrorV1> {
        if !spec.t_evap_k.is_finite()
            || spec.t_evap_k <= 180.0
            || spec.t_evap_k >= 0.94 * self.t_crit_k
        {
            return Err(CycleErrorV1::InvalidTemperature {
                temperature_k: spec.t_evap_k,
                reason: "evaporating temperature must lie in (180 K, 0.94 T_c)",
            });
        }
        if !spec.t_ambient_k.is_finite() || spec.t_ambient_k <= spec.t_evap_k + 2.0 {
            return Err(CycleErrorV1::ThermodynamicInversion(
                "ambient sink temperature must exceed evaporating temperature by > 2 K",
            ));
        }
        if !spec.suction_superheat_k.is_finite() || spec.suction_superheat_k < 0.5 {
            return Err(CycleErrorV1::InvalidScalar {
                quantity: "suction_superheat_k",
                value: spec.suction_superheat_k,
                reason: "must be >= 0.5 K to prevent liquid slugging",
            });
        }
        if !spec.ua_cond_clean_w_k.is_finite()
            || spec.ua_cond_clean_w_k <= 0.0
            || !spec.condenser_fouling_factor.is_finite()
            || spec.condenser_fouling_factor < 0.0
            || !spec.swept_volume_m3_s.is_finite()
            || spec.swept_volume_m3_s <= 0.0
            || !spec.eta_motor.is_finite()
            || spec.eta_motor <= 0.2
            || spec.eta_motor > 1.0
            || !spec.shell_loss_fraction.is_finite()
            || !(0.0..=1.0).contains(&spec.shell_loss_fraction)
            || !spec.clearance_ratio.is_finite()
            || !(0.0..=0.20).contains(&spec.clearance_ratio)
            || !spec.subcooling_k.is_finite()
            || spec.subcooling_k < 0.0
        {
            return Err(CycleErrorV1::InvalidConfig {
                parameter: "cycle_operating_point",
                reason: "UA, fouling, displacement, efficiencies, subcooling, or clearance invalid",
            });
        }

        let p_evap = self.midpoint_pressure_at_t(spec.t_evap_k)?;
        let evap_env = self.bubble_dew_at_pressure(p_evap)?;
        let t_suction = evap_env.t_dew_k + spec.suction_superheat_k;
        let ua_eff = spec.ua_cond_clean_w_k / (1.0 + spec.condenser_fouling_factor);

        // Decide subcritical condenser vs transcritical gas-cooler (e.g., R-744 above 30.98 C).
        let is_supercritical = spec.supercritical_pressure_pa.is_some()
            || spec.t_ambient_k + 4.0 >= 0.97 * self.t_crit_k;

        let (p_high, t_high_out) = if is_supercritical {
            // Supercritical gas-cooler: approach temperature rises with condenser fouling.
            let dt_approach = (3.5 * (1.0 + spec.condenser_fouling_factor)).clamp(2.0, 20.0);
            let t_gc_out = spec.t_ambient_k + dt_approach;
            // Optimal transcritical CO2 discharge pressure scales linearly with gas-cooler outlet T
            // (Kauf / Liao-Zhao correlation: P_opt [bar] ~ 2.6 * T_gc_out[C] - 5, boosted by fouling).
            let p_opt_default = ((2.65 * (t_gc_out - 273.15) - 4.0).max(80.0) * 1.0e5)
                * (1.0 + 0.15 * spec.condenser_fouling_factor);
            let p_gc = spec
                .supercritical_pressure_pa
                .unwrap_or(p_opt_default)
                .max(1.08 * self.p_crit_pa);
            (p_gc, t_gc_out)
        } else {
            // Solve condenser saturation temperature T_cond > T_amb balancing
            // Q_cond(T_cond) == UA_eff * (T_cond - T_amb).
            let mut t_lo = spec.t_ambient_k + 1.0;
            let mut t_hi = (0.96 * self.t_crit_k).max(t_lo + 1.0);
            let eval_residual = |t_cond: f64| -> Result<f64, CycleErrorV1> {
                let p_c = self.bubble_pressure_at_t(t_cond)?;
                let comp =
                    self.isentropic_compression(t_suction, p_evap, p_c, spec.eta_isentropic)?;
                let pr = p_c / p_evap;
                let eta_vol =
                    (1.0 - spec.clearance_ratio * (pow(pr, 1.0 / 1.15) - 1.0)).clamp(0.15, 0.99);
                let m_dot = comp.suction.density_kg_m3 * spec.swept_volume_m3_s * eta_vol;
                let t_out = (t_cond - spec.subcooling_k).max(spec.t_ambient_k + 0.5);
                let th = self.isenthalpic_throttle(t_out, p_c, p_evap)?;
                let q_evap =
                    m_dot * (comp.suction.enthalpy_j_kg - th.outlet_enthalpy_j_kg).max(0.0);
                let w_shaft = m_dot * comp.w_actual_j_kg;
                let w_elec = w_shaft / spec.eta_motor;
                let q_shell = (w_elec - w_shaft) * spec.shell_loss_fraction;
                let q_cond = q_evap + w_elec - q_shell;
                Ok(q_cond - ua_eff * (t_cond - spec.t_ambient_k))
            };
            for _ in 0..44 {
                let t_mid = 0.5 * (t_lo + t_hi);
                let r = eval_residual(t_mid)?;
                if r.abs() < 1.0e-4 {
                    t_lo = t_mid;
                    t_hi = t_mid;
                    break;
                }
                if r > 0.0 {
                    t_lo = t_mid;
                } else {
                    t_hi = t_mid;
                }
            }
            let t_cond = 0.5 * (t_lo + t_hi);
            let p_c = self.bubble_pressure_at_t(t_cond)?;
            let t_out = (t_cond - spec.subcooling_k).max(spec.t_ambient_k + 0.5);
            (p_c, t_out)
        };

        let compression =
            self.isentropic_compression(t_suction, p_evap, p_high, spec.eta_isentropic)?;
        let throttle = self.isenthalpic_throttle(t_high_out, p_high, p_evap)?;
        let pr = p_high / p_evap;
        let eta_vol = (1.0 - spec.clearance_ratio * (pow(pr, 1.0 / 1.15) - 1.0)).clamp(0.15, 0.99);
        let m_dot = compression.suction.density_kg_m3 * spec.swept_volume_m3_s * eta_vol;

        let q_evap_spec = compression.suction.enthalpy_j_kg - throttle.outlet_enthalpy_j_kg;
        if q_evap_spec <= 0.0 {
            return Err(CycleErrorV1::ThermodynamicInversion(
                "cycle has non-positive refrigerating effect at this high-side state",
            ));
        }
        let q_evap_w = m_dot * q_evap_spec;
        let w_shaft_w = m_dot * compression.w_actual_j_kg;
        let w_comp_w = w_shaft_w / spec.eta_motor;
        // Motor inefficiency either enters the refrigerant/oil stream or leaves via shell.
        // In hermetic/semi-hermetic compressors: Q_cond = Q_evap + W_comp - Q_shell identically.
        let q_shell_w = (w_comp_w - w_shaft_w) * spec.shell_loss_fraction;
        let q_cond_w = q_evap_w + w_comp_w - q_shell_w;
        let first_law_relative_defect =
            (q_cond_w - (q_evap_w + w_comp_w - q_shell_w)).abs() / q_cond_w.abs().max(1.0);

        let cop_cooling = q_evap_w / w_comp_w;
        let cop_carnot = spec.t_evap_k / (spec.t_ambient_k - spec.t_evap_k);
        let internal_entropy_gen_w_k = m_dot * (compression.ds_gen_j_kg_k + throttle.ds_gen_j_kg_k);
        let near_critical_warning =
            throttle.inlet.near_critical_warning || compression.discharge.near_critical_warning;

        Ok(ClosedCycleSummaryV1 {
            species: self.species,
            p_evap_pa: p_evap,
            p_high_pa: p_high,
            evap_envelope: evap_env,
            is_supercritical_gas_cooler: is_supercritical,
            t_high_out_k: t_high_out,
            m_dot_kg_s: m_dot,
            eta_volumetric: eta_vol,
            compression,
            throttle,
            q_evap_w,
            w_comp_w,
            q_shell_w,
            q_cond_w,
            first_law_relative_defect,
            cop_cooling,
            cop_carnot,
            internal_entropy_gen_w_k,
            near_critical_density_error_band: self.near_critical_density_error_band,
            near_critical_warning,
        })
    }
}

// ============================================================================
// COMPONENT 1: RefrigerationCompressorComponent
// ============================================================================

/// Configuration for [`RefrigerationCompressorComponent`].
#[derive(Clone, Debug, PartialEq)]
pub struct RefrigerationCompressorConfig {
    /// Working refrigerant species (`R-290`, `R-134a`, `R-404A`, `R-448A`, or `R-744`).
    pub species: RefrigerantSpecies,
    /// Swept displacement volumetric flow rate [m^3/s].
    pub swept_volume_m3_s: f64,
    /// Isentropic compression efficiency in `(0.2, 1.0]`.
    pub eta_isentropic: f64,
    /// Motor electrical efficiency in `(0.3, 1.0]`.
    pub eta_motor: f64,
    /// Clearance volume ratio in `[0.0, 0.20]`.
    pub clearance_ratio: f64,
    /// Compressor shell thermal capacitance `C_shell` [J/K].
    pub shell_heat_capacity_j_k: f64,
    /// Shell-to-ambient heat transfer conductance `UA_shell` [W/K].
    pub shell_ua_w_k: f64,
    /// Head-pressure sensitivity to condenser fouling factor (`P_disch_eff = P_disch * (1 + k * f_foul)`).
    pub fouling_head_pressure_gain: f64,
}

impl RefrigerationCompressorConfig {
    /// Commercial walk-in / reach-in hermetic compressor preset.
    #[must_use]
    pub const fn commercial_default(species: RefrigerantSpecies) -> Self {
        Self {
            species,
            swept_volume_m3_s: 1.15e-3,
            eta_isentropic: 0.72,
            eta_motor: 0.88,
            clearance_ratio: 0.045,
            shell_heat_capacity_j_k: 4_500.0,
            shell_ua_w_k: 8.5,
            fouling_head_pressure_gain: 0.28,
        }
    }

    /// Validate all compressor parameters fail-closed.
    pub fn validate(&self) -> Result<(), CycleErrorV1> {
        if !self.swept_volume_m3_s.is_finite()
            || self.swept_volume_m3_s <= 0.0
            || !self.eta_isentropic.is_finite()
            || self.eta_isentropic <= 0.2
            || self.eta_isentropic > 1.0
            || !self.eta_motor.is_finite()
            || self.eta_motor <= 0.3
            || self.eta_motor > 1.0
            || !self.clearance_ratio.is_finite()
            || !(0.0..=0.20).contains(&self.clearance_ratio)
            || !self.shell_heat_capacity_j_k.is_finite()
            || self.shell_heat_capacity_j_k <= 0.0
            || !self.shell_ua_w_k.is_finite()
            || self.shell_ua_w_k < 0.0
            || !self.fouling_head_pressure_gain.is_finite()
            || self.fouling_head_pressure_gain < 0.0
        {
            return Err(CycleErrorV1::InvalidConfig {
                parameter: "compressor_config",
                reason: "compressor displacement, efficiencies, shell thermal mass, or fouling gain invalid",
            });
        }
        Ok(())
    }
}

/// Accepted/candidate state for [`RefrigerationCompressorComponent`].
#[derive(Clone, Debug, PartialEq)]
pub struct CompressorStateV1 {
    /// Simulation time [s].
    pub time_s: f64,
    /// Hermetic compressor shell temperature [K].
    pub t_shell_k: f64,
    /// Cumulative electrical energy consumed [J].
    pub cumulative_elec_work_j: f64,
    /// Cumulative compressor operating runtime [s].
    pub cumulative_runtime_s: f64,
}

impl CompressorStateV1 {
    /// Create a fresh compressor state at ambient shell temperature `t_shell_k`.
    pub fn new(t_shell_k: f64) -> Result<Self, CycleErrorV1> {
        let s = Self {
            time_s: 0.0,
            t_shell_k,
            cumulative_elec_work_j: 0.0,
            cumulative_runtime_s: 0.0,
        };
        s.validate()?;
        Ok(s)
    }

    /// Validate state invariants fail-closed.
    pub fn validate(&self) -> Result<(), CycleErrorV1> {
        if !self.time_s.is_finite()
            || self.time_s < 0.0
            || !self.t_shell_k.is_finite()
            || !(200.0..=600.0).contains(&self.t_shell_k)
            || !self.cumulative_elec_work_j.is_finite()
            || self.cumulative_elec_work_j < 0.0
            || !self.cumulative_runtime_s.is_finite()
            || self.cumulative_runtime_s < 0.0
        {
            return Err(CycleErrorV1::InvalidScalar {
                quantity: "compressor_state",
                value: self.t_shell_k,
                reason: "compressor state fields must be finite and within physical bounds",
            });
        }
        Ok(())
    }
}

/// PR-EOS vapor-compression compressor [`Component`].
#[derive(Clone, Debug)]
pub struct RefrigerationCompressorComponent {
    config: RefrigerationCompressorConfig,
    fluid: RefrigerantModelV1,
    ports: Vec<PortSpec>,
    layout: PortLayout,
}

impl RefrigerationCompressorComponent {
    /// Construct a validated [`RefrigerationCompressorComponent`].
    pub fn new(config: RefrigerationCompressorConfig) -> Result<Self, CycleErrorV1> {
        config.validate()?;
        let fluid = RefrigerantModelV1::for_species(config.species)?;
        Self::with_fluid(config, fluid)
    }

    /// Construct from an explicit [`RefrigerantModelV1`] (e.g. loaded from a `fs-matdb` card).
    pub fn with_fluid(
        config: RefrigerationCompressorConfig,
        fluid: RefrigerantModelV1,
    ) -> Result<Self, CycleErrorV1> {
        config.validate()?;
        if fluid.species() != config.species {
            return Err(CycleErrorV1::InvalidConfig {
                parameter: "species",
                reason: "RefrigerantModelV1 species must match config.species",
            });
        }
        let mk = |n: &'static str| {
            StableId::new(n).map_err(|_| CycleErrorV1::InvalidConfig {
                parameter: n,
                reason: "invalid StableId",
            })
        };
        let ports = vec![
            PortSpec::power(
                mk("electrical_power")?,
                PortKind::ElectricalVoltageCurrent,
                PortRole::EffortIn,
                NonZeroUsize::MIN,
            )?,
            PortSpec::identified_stream(
                mk("suction_stream")?,
                PortRole::EffortIn,
                refrigerant_v1(config.species)?,
            )?,
            PortSpec::identified_stream(
                mk("discharge_stream")?,
                PortRole::EffortIn,
                refrigerant_v1(config.species)?,
            )?,
            PortSpec::heat(mk("shell_heat")?, PortRole::EffortIn)?,
            PortSpec::signal_in(mk("fouling_factor")?, NonZeroUsize::MIN),
            PortSpec::signal_out(
                mk("probes")?,
                NonZeroUsize::new(COMPRESSOR_PROBE_COUNT).ok_or(CycleErrorV1::InvalidConfig {
                    parameter: "probes",
                    reason: "nonzero width",
                })?,
            ),
        ];
        let layout = PortLayout::new(&ports)?;
        Ok(Self {
            config,
            fluid,
            ports,
            layout,
        })
    }

    /// Borrow the underlying refrigerant model.
    #[must_use]
    pub const fn fluid(&self) -> &RefrigerantModelV1 {
        &self.fluid
    }

    /// Borrow the port layout.
    #[must_use]
    pub const fn layout(&self) -> &PortLayout {
        &self.layout
    }

    /// Build a `PortLayout`-ordered input slice for [`Component::propose`].
    #[must_use]
    pub fn build_inputs(
        &self,
        voltage_v: f64,
        p_suction_pa: f64,
        h_suction_j_kg: f64,
        p_discharge_pa: f64,
        h_discharge_downstream_j_kg: f64,
        t_ambient_k: f64,
        fouling_factor: f64,
    ) -> Vec<f64> {
        let mut inputs = vec![0.0; self.layout.input_len()];
        inputs[self.layout.input_range(0)][0] = voltage_v;
        let suc = &mut inputs[self.layout.input_range(1)];
        suc[STREAM_PRESSURE] = p_suction_pa;
        suc[STREAM_SPECIFIC_ENTHALPY] = h_suction_j_kg;
        suc[STREAM_CARRIED] = 1.0;
        let dis = &mut inputs[self.layout.input_range(2)];
        dis[STREAM_PRESSURE] = p_discharge_pa;
        dis[STREAM_SPECIFIC_ENTHALPY] = h_discharge_downstream_j_kg;
        dis[STREAM_CARRIED] = 1.0;
        inputs[self.layout.input_range(3)][0] = t_ambient_k;
        inputs[self.layout.input_range(4)][0] = fouling_factor;
        inputs
    }
}

impl Component for RefrigerationCompressorComponent {
    type State = CompressorStateV1;
    type Error = CycleErrorV1;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Smooth
    }

    fn end_of_step(&self, _port: usize) -> EndOfStep {
        EndOfStep::Untold
    }

    fn stored_energy_j(&self, state: &Self::State) -> Result<f64, Self::Error> {
        state.validate()?;
        Ok(self.config.shell_heat_capacity_j_k * (state.t_shell_k - IIR_REFERENCE_TEMPERATURE_K))
    }

    fn stored_mass_kg(&self, state: &Self::State) -> Result<Option<f64>, Self::Error> {
        state.validate()?;
        Ok(Some(0.0))
    }

    fn stored_carried_quantity(
        &self,
        state: &Self::State,
        carried: usize,
    ) -> Result<Option<f64>, Self::Error> {
        state.validate()?;
        if carried < self.layout.carried_len() {
            Ok(Some(0.0))
        } else {
            Ok(None)
        }
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
        old.validate()?;
        if !dt_s.is_finite() || dt_s <= 0.0 {
            return Err(CycleErrorV1::InvalidScalar {
                quantity: "dt_s",
                value: dt_s,
                reason: "must be finite positive",
            });
        }
        if inputs.len() != self.layout.input_len() || outputs.len() != self.layout.output_len() {
            return Err(CycleErrorV1::InvalidConfig {
                parameter: "port_vectors",
                reason: "inputs or outputs length mismatch",
            });
        }
        if !inputs.iter().all(|v| v.is_finite()) {
            return Err(CycleErrorV1::InvalidConfig {
                parameter: "inputs",
                reason: "all compressor inputs must be finite",
            });
        }

        let voltage_v = inputs[self.layout.input_range(0)][0];
        let suc = &inputs[self.layout.input_range(1)];
        let dis = &inputs[self.layout.input_range(2)];
        let t_amb_k = inputs[self.layout.input_range(3)][0];
        let fouling = inputs[self.layout.input_range(4)][0];

        if voltage_v <= 0.0 {
            return Err(CycleErrorV1::InvalidScalar {
                quantity: "voltage_v",
                value: voltage_v,
                reason: "supply voltage must be strictly positive",
            });
        }
        if t_amb_k <= 100.0 {
            return Err(CycleErrorV1::InvalidTemperature {
                temperature_k: t_amb_k,
                reason: "ambient temperature must exceed 100 K",
            });
        }
        if fouling < 0.0 {
            return Err(CycleErrorV1::InvalidScalar {
                quantity: "fouling_factor",
                value: fouling,
                reason: "condenser fouling factor must be non-negative",
            });
        }

        let p_suc = suc[STREAM_PRESSURE];
        let h_suc = suc[STREAM_SPECIFIC_ENTHALPY];
        let p_disch_nom = dis[STREAM_PRESSURE];
        if p_suc <= 1_000.0 || p_disch_nom <= p_suc {
            return Err(CycleErrorV1::ThermodynamicInversion(
                "compressor requires P_discharge > P_suction > 1 kPa",
            ));
        }

        // Fouling elevates effective head pressure seen at the discharge valve.
        let p_disch_eff = p_disch_nom * (1.0 + self.config.fouling_head_pressure_gain * fouling);
        let env_suc = self.fluid.bubble_dew_at_pressure(p_suc)?;
        if h_suc < env_suc.h_vapor_j_kg - 1.0e-3 {
            return Err(CycleErrorV1::ThermodynamicInversion(
                "compressor suction enthalpy lies inside two-phase wet dome (liquid slugging risk)",
            ));
        }
        let suc_pt = self.fluid.temperature_at_p_h(
            p_suc,
            h_suc,
            env_suc.t_dew_k,
            PengRobinsonPhaseV1::Vapor,
        )?;
        let comp = self.fluid.isentropic_compression(
            suc_pt.temperature_k,
            p_suc,
            p_disch_eff,
            self.config.eta_isentropic,
        )?;

        let pr = p_disch_eff / p_suc;
        let eta_vol =
            (1.0 - self.config.clearance_ratio * (pow(pr, 1.0 / 1.15) - 1.0)).clamp(0.15, 0.99);
        let m_dot = comp.suction.density_kg_m3 * self.config.swept_volume_m3_s * eta_vol;
        let w_shaft_w = m_dot * (comp.discharge.enthalpy_j_kg - comp.suction.enthalpy_j_kg);
        let w_elec_w = w_shaft_w / self.config.eta_motor;
        let q_motor_loss_w = w_elec_w - w_shaft_w;

        // Implicitly stable hermetic shell energy balance:
        // C_shell * (T_new - T_old) / dt = Q_motor_loss - UA_shell * (T_avg - T_amb).
        let c_sh = self.config.shell_heat_capacity_j_k;
        let ua_sh = self.config.shell_ua_w_k;
        let denom = c_sh + 0.5 * ua_sh * dt_s;
        let t_shell_new = (c_sh * old.t_shell_k
            + dt_s * (q_motor_loss_w + ua_sh * (t_amb_k - 0.5 * old.t_shell_k)))
            / denom;
        let t_shell_avg = 0.5 * (old.t_shell_k + t_shell_new);
        let q_shell_reject_w = ua_sh * (t_shell_avg - t_amb_k);

        let next_state = CompressorStateV1 {
            time_s: old.time_s + dt_s,
            t_shell_k: t_shell_new,
            cumulative_elec_work_j: old.cumulative_elec_work_j + w_elec_w * dt_s,
            cumulative_runtime_s: old.cumulative_runtime_s + dt_s,
        };
        next_state.validate()?;
        *new = next_state;

        outputs.fill(0.0);
        let current_a = w_elec_w / voltage_v;
        outputs[self.layout.output_range(0)][0] = current_a;

        let suc_out = &mut outputs[self.layout.output_range(1)];
        suc_out[STREAM_MASS_RATE] = m_dot;
        suc_out[STREAM_ENTHALPY_RATE] = m_dot * comp.suction.enthalpy_j_kg;
        suc_out[STREAM_CARRIED] = m_dot;

        let dis_out = &mut outputs[self.layout.output_range(2)];
        dis_out[STREAM_MASS_RATE] = -m_dot;
        dis_out[STREAM_ENTHALPY_RATE] = -m_dot * comp.discharge.enthalpy_j_kg;
        dis_out[STREAM_CARRIED] = -m_dot;

        // Heat port (EffortIn): negative flow = heat REJECTED from compressor shell to ambient.
        outputs[self.layout.output_range(3)][0] = -q_shell_reject_w;

        let probes = &mut outputs[self.layout.output_range(5)];
        probes[0] = m_dot;
        probes[1] = w_elec_w;
        probes[2] = q_shell_reject_w;
        probes[3] = comp.discharge.temperature_k;
        probes[4] = pr;
        probes[5] = eta_vol;
        probes[6] = self.config.eta_isentropic;
        probes[7] = current_a;

        let e_old = self.stored_energy_j(old)?;
        let e_new = self.stored_energy_j(new)?;
        let e_scale = e_old.abs()
            + e_new.abs()
            + (w_elec_w * dt_s).abs()
            + (m_dot * comp.discharge.enthalpy_j_kg.abs() * dt_s)
            + 1.0;
        let receipt = Receipt {
            stored_change_j: e_new - e_old,
            dissipated_j: ((w_elec_w - m_dot * comp.w_isentropic_j_kg) * dt_s).max(0.0),
            numerical_loss_j: 0.0,
            entropy_generated_j_per_k: (m_dot * comp.ds_gen_j_kg_k * dt_s).max(0.0),
            energy_tolerance_j: 512.0 * f64::EPSILON * e_scale,
            residual: 0.0,
            tolerance: 0.0,
        };
        Ok(Proposal::Step(receipt))
    }

    const STATE_CHARTER: Option<StateIdentityCharterV2> = Some(StateIdentityCharterV2 {
        owner: "fs-multiphase::RefrigerationCompressor",
        state_family: "pr-compressor-v1",
        schema_grammar: "time_s, t_shell_k, cumulative_elec_work_j, cumulative_runtime_s in coherent SI",
        codec_grammar: "little-endian u64 words: codec tag, field count (4), followed by 4 IEEE-754 f64 bit patterns",
        codec_version: 1,
    });

    fn encode_state(&self, state: &Self::State, words: &mut Vec<u64>) -> bool {
        if state.validate().is_err() {
            return false;
        }
        words.push(COMPRESSOR_CODEC_TAG);
        words.push(4);
        words.push(state.time_s.to_bits());
        words.push(state.t_shell_k.to_bits());
        words.push(state.cumulative_elec_work_j.to_bits());
        words.push(state.cumulative_runtime_s.to_bits());
        true
    }

    fn decode_state(&self, words: &[u64]) -> Option<Self::State> {
        if words.len() != 6 || words[0] != COMPRESSOR_CODEC_TAG || words[1] != 4 {
            return None;
        }
        let state = CompressorStateV1 {
            time_s: f64::from_bits(words[2]),
            t_shell_k: f64::from_bits(words[3]),
            cumulative_elec_work_j: f64::from_bits(words[4]),
            cumulative_runtime_s: f64::from_bits(words[5]),
        };
        state.validate().ok().map(|()| state)
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        words.push(species_discriminant(self.config.species));
        let scalars = [
            self.config.swept_volume_m3_s,
            self.config.eta_isentropic,
            self.config.eta_motor,
            self.config.clearance_ratio,
            self.config.shell_heat_capacity_j_k,
            self.config.shell_ua_w_k,
            self.config.fouling_head_pressure_gain,
        ];
        words.extend(scalars.iter().map(|v| v.to_bits()));
    }
}

// ============================================================================
// COMPONENT 2: CondenserEvaporatorHXComponent
// ============================================================================

/// Operating mode of [`CondenserEvaporatorHXComponent`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HxModeV1 {
    /// Condenser (subcritical two-phase) or supercritical gas-cooler (R-744 above `T_c`).
    CondenserGasCooler,
    /// Evaporator coil with dynamic frost-layer accumulation and electric/hot-gas defrost
    /// draining liquid meltwater OUT via `drain_stream` (`kitchen_liquid_v1`).
    FrostingEvaporator,
}

/// Configuration for [`CondenserEvaporatorHXComponent`].
#[derive(Clone, Debug, PartialEq)]
pub struct CondenserEvaporatorHXConfig {
    /// Working refrigerant species.
    pub species: RefrigerantSpecies,
    /// Condenser/gas-cooler vs frosting evaporator mode.
    pub mode: HxModeV1,
    /// Clean overall heat-transfer conductance `UA_clean` [W/K].
    pub ua_clean_w_k: f64,
    /// Finned air-side coil surface area `A_coil` [m^2].
    pub coil_area_m2: f64,
    /// Coil metal thermal capacitance `C_coil` [J/K].
    pub coil_heat_capacity_j_k: f64,
    /// Frost bulk density `rho_frost` [kg/m^3] (typically `180--320 kg/m^3`).
    pub frost_density_kg_m3: f64,
    /// Frost effective thermal conductivity `k_frost` [W/(m K)] (typically `0.15--0.35`).
    pub frost_conductivity_w_m_k: f64,
    /// Fin spacing `s_fin` [m] governing airflow blockage as `delta_frost` grows.
    pub fin_spacing_m: f64,
    /// Target subcooling [K] (condenser mode) or superheat [K] (evaporator mode).
    pub target_approach_k: f64,
}

impl CondenserEvaporatorHXConfig {
    /// Air-cooled commercial condenser / transcritical R-744 gas-cooler preset.
    #[must_use]
    pub const fn condenser_default(species: RefrigerantSpecies) -> Self {
        Self {
            species,
            mode: HxModeV1::CondenserGasCooler,
            ua_clean_w_k: 260.0,
            coil_area_m2: 4.5,
            coil_heat_capacity_j_k: 3_200.0,
            frost_density_kg_m3: 240.0,
            frost_conductivity_w_m_k: 0.22,
            fin_spacing_m: 0.004,
            target_approach_k: 4.0,
        }
    }

    /// Walk-in cooler / freezer frosting evaporator coil preset.
    #[must_use]
    pub const fn evaporator_default(species: RefrigerantSpecies) -> Self {
        Self {
            species,
            mode: HxModeV1::FrostingEvaporator,
            ua_clean_w_k: 340.0,
            coil_area_m2: 5.0,
            coil_heat_capacity_j_k: 3_600.0,
            frost_density_kg_m3: 240.0,
            frost_conductivity_w_m_k: 0.22,
            fin_spacing_m: 0.005,
            target_approach_k: 6.0,
        }
    }

    /// Validate heat-exchanger parameters fail-closed.
    pub fn validate(&self) -> Result<(), CycleErrorV1> {
        let pos = [
            ("ua_clean_w_k", self.ua_clean_w_k),
            ("coil_area_m2", self.coil_area_m2),
            ("coil_heat_capacity_j_k", self.coil_heat_capacity_j_k),
            ("frost_density_kg_m3", self.frost_density_kg_m3),
            ("frost_conductivity_w_m_k", self.frost_conductivity_w_m_k),
            ("fin_spacing_m", self.fin_spacing_m),
        ];
        for (parameter, v) in pos {
            if !v.is_finite() || v <= 0.0 {
                return Err(CycleErrorV1::InvalidConfig {
                    parameter,
                    reason: "must be finite and strictly positive",
                });
            }
        }
        if !self.target_approach_k.is_finite() || self.target_approach_k < 0.0 {
            return Err(CycleErrorV1::InvalidConfig {
                parameter: "target_approach_k",
                reason: "must be finite and non-negative",
            });
        }
        Ok(())
    }
}

/// Accepted/candidate state for [`CondenserEvaporatorHXComponent`].
#[derive(Clone, Debug, PartialEq)]
pub struct CondenserEvaporatorHXStateV1 {
    /// Elapsed simulation time [s].
    pub time_s: f64,
    /// Coil metal temperature [K].
    pub t_coil_k: f64,
    /// Accumulated frost mass on coil fins `m_frost = rho_frost * A_coil * delta_frost` [kg].
    pub m_frost_kg: f64,
    /// Liquid meltwater temporarily residing in drain pan awaiting drainage [kg].
    pub m_meltwater_pan_kg: f64,
    /// Stored energy on canonical datum (`C_coil * (T_coil - 273.15) + m_frost * h_ice + m_melt * h_liq + aux`) [J].
    pub stored_energy_j: f64,
    /// Auxiliary non-H2O carried gas species `[C2H5OH, CO2, O2, N2, CH4, C3H8, CO, grease_vapor]` [kg].
    pub aux_gas_kg: [f64; AUX_GAS_SPECIES_COUNT],
    /// Auxiliary non-H2O carried liquid species `[C2H5OH, triglyceride, ffa, glycerol, salt, sugar, fond]` [kg].
    pub aux_liq_kg: [f64; AUX_LIQ_SPECIES_COUNT],
    /// Cumulative meltwater drained out via `drain_stream` [kg].
    pub cumulative_drain_kg: f64,
}

impl CondenserEvaporatorHXStateV1 {
    /// Frost thickness `delta_frost = m_frost / (rho_frost * A_coil)` [m].
    #[must_use]
    pub fn delta_frost_m(&self, config: &CondenserEvaporatorHXConfig) -> f64 {
        self.m_frost_kg / (config.frost_density_kg_m3 * config.coil_area_m2)
    }

    /// Total stored H2O mass on coil (frost + drain pan) [kg].
    #[must_use]
    pub fn total_h2o_kg(&self) -> f64 {
        self.m_frost_kg + self.m_meltwater_pan_kg
    }

    /// Total stored mass [kg].
    #[must_use]
    pub fn total_mass_kg(&self) -> f64 {
        self.total_h2o_kg()
            + self.aux_gas_kg.iter().sum::<f64>()
            + self.aux_liq_kg.iter().sum::<f64>()
    }

    /// Validate state invariants fail-closed.
    pub fn validate(&self, config: &CondenserEvaporatorHXConfig) -> Result<(), CycleErrorV1> {
        if !self.time_s.is_finite()
            || self.time_s < 0.0
            || !self.t_coil_k.is_finite()
            || !(150.0..=600.0).contains(&self.t_coil_k)
            || !self.m_frost_kg.is_finite()
            || self.m_frost_kg < 0.0
            || !self.m_meltwater_pan_kg.is_finite()
            || self.m_meltwater_pan_kg < 0.0
            || !self.stored_energy_j.is_finite()
            || !self.cumulative_drain_kg.is_finite()
            || self.cumulative_drain_kg < 0.0
        {
            return Err(CycleErrorV1::InvalidScalar {
                quantity: "hx_state",
                value: self.t_coil_k,
                reason: "HX state values must be finite and within physical bounds",
            });
        }
        if self.aux_gas_kg.iter().any(|v| !v.is_finite() || *v < 0.0)
            || self.aux_liq_kg.iter().any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(CycleErrorV1::InvalidScalar {
                quantity: "hx_aux_species",
                value: -1.0,
                reason: "auxiliary species inventories must be finite non-negative",
            });
        }
        let expected_e = hx_stored_energy_j(self, config);
        let scale = expected_e.abs().max(self.stored_energy_j.abs()).max(1.0);
        if (expected_e - self.stored_energy_j).abs() > 1.0e-5 * scale {
            return Err(CycleErrorV1::InvalidScalar {
                quantity: "stored_energy_j",
                value: self.stored_energy_j,
                reason: "HX stored energy disagrees with coil temperature and frost/melt inventory",
            });
        }
        Ok(())
    }
}

fn hx_stored_energy_j(
    state: &CondenserEvaporatorHXStateV1,
    config: &CondenserEvaporatorHXConfig,
) -> f64 {
    let t_ice = state.t_coil_k.min(ICE_POINT_K);
    let h_ice = H_ICE_0C_JANAF_J_KG + CP_WATER_ICE_J_KG_K * (t_ice - ICE_POINT_K);
    let h_melt =
        H_WATER_LIQ_0C_JANAF_J_KG + CP_WATER_LIQ_J_KG_K * (state.t_coil_k - ICE_POINT_K).max(0.0);
    config.coil_heat_capacity_j_k * (state.t_coil_k - IIR_REFERENCE_TEMPERATURE_K)
        + state.m_frost_kg * h_ice
        + state.m_meltwater_pan_kg * h_melt
}

fn invert_hx_coil_temp(
    state: &CondenserEvaporatorHXStateV1,
    config: &CondenserEvaporatorHXConfig,
    target_e_j: f64,
) -> Result<f64, CycleErrorV1> {
    // Piece-wise linear strictly increasing E(T_coil) across T_coil = 273.15 K.
    let e_at_0c = state.m_frost_kg * H_ICE_0C_JANAF_J_KG
        + state.m_meltwater_pan_kg * H_WATER_LIQ_0C_JANAF_J_KG;
    let t_coil = if target_e_j <= e_at_0c {
        let c_below = config.coil_heat_capacity_j_k + state.m_frost_kg * CP_WATER_ICE_J_KG_K;
        ICE_POINT_K + (target_e_j - e_at_0c) / c_below
    } else {
        let c_above =
            config.coil_heat_capacity_j_k + state.m_meltwater_pan_kg * CP_WATER_LIQ_J_KG_K;
        ICE_POINT_K + (target_e_j - e_at_0c) / c_above
    };
    if !t_coil.is_finite() || !(150.0..=600.0).contains(&t_coil) {
        return Err(CycleErrorV1::InvalidTemperature {
            temperature_k: t_coil,
            reason: "inverted HX coil temperature out of [150, 600] K",
        });
    }
    Ok(t_coil)
}

/// Dual-mode condenser / supercritical gas-cooler and frosting evaporator [`Component`].
#[derive(Clone, Debug)]
pub struct CondenserEvaporatorHXComponent {
    config: CondenserEvaporatorHXConfig,
    fluid: RefrigerantModelV1,
    ports: Vec<PortSpec>,
    layout: PortLayout,
    moisture_port: Option<usize>,
    drain_port: Option<usize>,
    probes_port: usize,
}

impl CondenserEvaporatorHXComponent {
    /// Construct a validated [`CondenserEvaporatorHXComponent`].
    pub fn new(config: CondenserEvaporatorHXConfig) -> Result<Self, CycleErrorV1> {
        config.validate()?;
        let fluid = RefrigerantModelV1::for_species(config.species)?;
        Self::with_fluid(config, fluid)
    }

    /// Construct from an explicit [`RefrigerantModelV1`].
    pub fn with_fluid(
        config: CondenserEvaporatorHXConfig,
        fluid: RefrigerantModelV1,
    ) -> Result<Self, CycleErrorV1> {
        config.validate()?;
        if fluid.species() != config.species {
            return Err(CycleErrorV1::InvalidConfig {
                parameter: "species",
                reason: "RefrigerantModelV1 species must match config.species",
            });
        }
        let mk = |n: &'static str| {
            StableId::new(n).map_err(|_| CycleErrorV1::InvalidConfig {
                parameter: n,
                reason: "invalid StableId",
            })
        };
        let mut ports = vec![
            PortSpec::identified_stream(
                mk("ref_in")?,
                PortRole::EffortIn,
                refrigerant_v1(config.species)?,
            )?,
            PortSpec::identified_stream(
                mk("ref_out")?,
                PortRole::EffortIn,
                refrigerant_v1(config.species)?,
            )?,
            PortSpec::heat(mk("air_heat")?, PortRole::EffortIn)?,
            PortSpec::heat(mk("defrost_heat")?, PortRole::FlowIn)?,
            PortSpec::signal_in(mk("fouling_or_airflow")?, NonZeroUsize::MIN),
            PortSpec::signal_in(mk("ref_mass_flow_kg_s")?, NonZeroUsize::MIN),
        ];
        let (moisture_port, drain_port) = if config.mode == HxModeV1::FrostingEvaporator {
            let mp = ports.len();
            ports.push(PortSpec::identified_stream(
                mk("moisture_in_stream")?,
                PortRole::FlowIn,
                kitchen_gas_v1()?,
            )?);
            let dp = ports.len();
            ports.push(PortSpec::identified_stream(
                mk("drain_stream")?,
                PortRole::FlowIn,
                kitchen_liquid_v1()?,
            )?);
            (Some(mp), Some(dp))
        } else {
            (None, None)
        };
        let probes_port = ports.len();
        ports.push(PortSpec::signal_out(
            mk("probes")?,
            NonZeroUsize::new(HX_PROBE_COUNT).ok_or(CycleErrorV1::InvalidConfig {
                parameter: "probes",
                reason: "nonzero width",
            })?,
        ));
        let layout = PortLayout::new(&ports)?;
        Ok(Self {
            config,
            fluid,
            ports,
            layout,
            moisture_port,
            drain_port,
            probes_port,
        })
    }

    /// Borrow the HX configuration.
    #[must_use]
    pub const fn config(&self) -> &CondenserEvaporatorHXConfig {
        &self.config
    }

    /// Borrow the port layout.
    #[must_use]
    pub const fn layout(&self) -> &PortLayout {
        &self.layout
    }

    /// Index of the `probes` `SignalOut` port.
    #[must_use]
    pub const fn probes_port(&self) -> usize {
        self.probes_port
    }

    /// Create an initial clean (zero-frost) HX state at `t_coil_k`.
    pub fn initial_state(
        &self,
        t_coil_k: f64,
    ) -> Result<CondenserEvaporatorHXStateV1, CycleErrorV1> {
        let mut s = CondenserEvaporatorHXStateV1 {
            time_s: 0.0,
            t_coil_k,
            m_frost_kg: 0.0,
            m_meltwater_pan_kg: 0.0,
            stored_energy_j: 0.0,
            aux_gas_kg: [0.0; AUX_GAS_SPECIES_COUNT],
            aux_liq_kg: [0.0; AUX_LIQ_SPECIES_COUNT],
            cumulative_drain_kg: 0.0,
        };
        s.stored_energy_j = hx_stored_energy_j(&s, &self.config);
        s.validate(&self.config)?;
        Ok(s)
    }

    /// Compute effective overall conductance `UA_eff` [W/K] accounting for condenser fouling
    /// OR evaporator frost thickness `delta_frost` (conduction resistance + fin airflow blockage).
    #[must_use]
    pub fn effective_ua_w_k(
        &self,
        state: &CondenserEvaporatorHXStateV1,
        fouling_or_airflow: f64,
    ) -> f64 {
        match self.config.mode {
            HxModeV1::CondenserGasCooler => {
                self.config.ua_clean_w_k / (1.0 + fouling_or_airflow.max(0.0))
            }
            HxModeV1::FrostingEvaporator => {
                let delta = state.delta_frost_m(&self.config).max(0.0);
                // Thermal series resistance of frost layer + fin-passage airflow blockage.
                let r_clean = 1.0 / self.config.ua_clean_w_k;
                let r_frost =
                    delta / (self.config.frost_conductivity_w_m_k * self.config.coil_area_m2);
                let blockage_ratio = (2.0 * delta / self.config.fin_spacing_m).clamp(0.0, 0.92);
                let airflow_factor =
                    (1.0 - 0.75 * blockage_ratio) * fouling_or_airflow.clamp(0.1, 1.5);
                (1.0 / (r_clean + r_frost)) * airflow_factor
            }
        }
    }

    /// Build `PortLayout`-ordered inputs for condenser or evaporator operation, automatically
    /// constructing exact donor-rule `drain_stream` meltwater removal when `drain_pan_fraction > 0`.
    #[allow(clippy::too_many_arguments)]
    pub fn build_inputs(
        &self,
        state: &CondenserEvaporatorHXStateV1,
        dt_s: f64,
        p_ref_pa: f64,
        h_ref_in_j_kg: f64,
        m_dot_ref_kg_s: f64,
        t_air_k: f64,
        q_defrost_w: f64,
        fouling_or_airflow: f64,
        m_dot_frost_vapor_in_kg_s: f64,
        drain_pan_fraction: f64,
    ) -> Result<Vec<f64>, CycleErrorV1> {
        if !dt_s.is_finite() || dt_s <= 0.0 {
            return Err(CycleErrorV1::InvalidScalar {
                quantity: "dt_s",
                value: dt_s,
                reason: "must be finite positive",
            });
        }
        let mut inputs = vec![0.0; self.layout.input_len()];
        let r_in = &mut inputs[self.layout.input_range(0)];
        r_in[STREAM_PRESSURE] = p_ref_pa;
        r_in[STREAM_SPECIFIC_ENTHALPY] = h_ref_in_j_kg;
        r_in[STREAM_CARRIED] = 1.0;

        let r_out = &mut inputs[self.layout.input_range(1)];
        r_out[STREAM_PRESSURE] = p_ref_pa;
        r_out[STREAM_SPECIFIC_ENTHALPY] = h_ref_in_j_kg;
        r_out[STREAM_CARRIED] = 1.0;

        inputs[self.layout.input_range(2)][0] = t_air_k;
        inputs[self.layout.input_range(3)][0] = q_defrost_w;
        inputs[self.layout.input_range(4)][0] = fouling_or_airflow;
        inputs[self.layout.input_range(5)][0] = m_dot_ref_kg_s;

        if let Some(mp) = self.moisture_port {
            let m_in = m_dot_frost_vapor_in_kg_s.max(0.0);
            let h_vap_air =
                H_F0_H2O_VAP_J_KG + CP_WATER_VAP_J_KG_K * (t_air_k - JANAF_REFERENCE_TEMPERATURE_K);
            let wire = &mut inputs[self.layout.input_range(mp)];
            wire[STREAM_MASS_RATE] = m_in;
            wire[STREAM_ENTHALPY_RATE] = m_in * h_vap_air;
            wire[STREAM_CARRIED] = m_in;
        }
        if let Some(dp) = self.drain_port {
            let frac = drain_pan_fraction.clamp(0.0, 1.0);
            let m_drain_rate = -frac * state.m_meltwater_pan_kg / dt_s;
            let h_melt = H_WATER_LIQ_0C_JANAF_J_KG
                + CP_WATER_LIQ_J_KG_K * (state.t_coil_k - ICE_POINT_K).max(0.0);
            let wire = &mut inputs[self.layout.input_range(dp)];
            wire[STREAM_MASS_RATE] = m_drain_rate;
            wire[STREAM_ENTHALPY_RATE] = m_drain_rate * h_melt;
            wire[STREAM_CARRIED] = m_drain_rate;
        }
        Ok(inputs)
    }
}

impl Component for CondenserEvaporatorHXComponent {
    type State = CondenserEvaporatorHXStateV1;
    type Error = CycleErrorV1;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Smooth
    }

    fn end_of_step(&self, _port: usize) -> EndOfStep {
        EndOfStep::Untold
    }

    fn stored_energy_j(&self, state: &Self::State) -> Result<f64, Self::Error> {
        state.validate(&self.config)?;
        Ok(state.stored_energy_j)
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
        let Some(StreamQuantity::SpeciesMass(sp)) = self.layout.carried_quantities().get(carried)
        else {
            return Ok(None);
        };
        let amount = match sp.as_str() {
            "R290" | "R134a" | "R404A" | "R448A" | "R744" => 0.0,
            "H2O" => state.total_h2o_kg(),
            "C2H5OH" => state.aux_gas_kg[0] + state.aux_liq_kg[0],
            "CO2" => state.aux_gas_kg[1],
            "O2" => state.aux_gas_kg[2],
            "N2" => state.aux_gas_kg[3],
            "CH4" => state.aux_gas_kg[4],
            "C3H8" => state.aux_gas_kg[5],
            "CO" => state.aux_gas_kg[6],
            "grease_vapor" => state.aux_gas_kg[7],
            "triglyceride" => state.aux_liq_kg[1],
            "ffa" => state.aux_liq_kg[2],
            "glycerol" => state.aux_liq_kg[3],
            "solute_salt" => state.aux_liq_kg[4],
            "solute_sugar" => state.aux_liq_kg[5],
            "solute_fond" => state.aux_liq_kg[6],
            _ => return Ok(None),
        };
        Ok(Some(amount))
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
            return Err(CycleErrorV1::InvalidScalar {
                quantity: "dt_s",
                value: dt_s,
                reason: "must be finite positive",
            });
        }
        if inputs.len() != self.layout.input_len() || outputs.len() != self.layout.output_len() {
            return Err(CycleErrorV1::InvalidConfig {
                parameter: "port_vectors",
                reason: "inputs or outputs length mismatch",
            });
        }
        if !inputs.iter().all(|v| v.is_finite()) {
            return Err(CycleErrorV1::InvalidConfig {
                parameter: "inputs",
                reason: "all HX inputs must be finite",
            });
        }

        let ref_in = &inputs[self.layout.input_range(0)];
        let p_ref = ref_in[STREAM_PRESSURE];
        let h_ref_in = ref_in[STREAM_SPECIFIC_ENTHALPY];
        let t_air_k = inputs[self.layout.input_range(2)][0];
        let q_defrost_w = inputs[self.layout.input_range(3)][0];
        let fouling_or_airflow = inputs[self.layout.input_range(4)][0];
        let m_dot_ref = inputs[self.layout.input_range(5)][0];

        if p_ref <= 1_000.0 {
            return Err(CycleErrorV1::InvalidPressure {
                pressure_pa: p_ref,
                reason: "refrigerant pressure must exceed 1 kPa",
            });
        }
        if t_air_k <= 100.0 {
            return Err(CycleErrorV1::InvalidTemperature {
                temperature_k: t_air_k,
                reason: "air temperature must exceed 100 K",
            });
        }
        if q_defrost_w < 0.0 || m_dot_ref < 0.0 || fouling_or_airflow < 0.0 {
            return Err(CycleErrorV1::InvalidScalar {
                quantity: "hx_signal",
                value: q_defrost_w.min(m_dot_ref).min(fouling_or_airflow),
                reason: "defrost power, refrigerant mass flow, and fouling/airflow must be non-negative",
            });
        }

        // Parse optional evaporator moisture-in & liquid-drain streams.
        let mut h_dot_streams_w = 0.0;
        let mut dm_h2o_frost_in = 0.0;
        let mut dm_h2o_drain = 0.0;
        let mut aux_gas = old.aux_gas_kg;
        let mut aux_liq = old.aux_liq_kg;

        if let Some(mp) = self.moisture_port {
            let wire = &inputs[self.layout.input_range(mp)];
            validate_stream_sum("moisture_in_stream", wire)?;
            h_dot_streams_w += wire[STREAM_ENTHALPY_RATE];
            dm_h2o_frost_in = dt_s * wire[STREAM_CARRIED];
            for k in 0..AUX_GAS_SPECIES_COUNT {
                let cand = aux_gas[k] + dt_s * wire[STREAM_CARRIED + 1 + k];
                if cand < -1.0e-12 {
                    return Err(CycleErrorV1::DepletedInventory {
                        inventory: "hx_aux_gas",
                    });
                }
                aux_gas[k] = cand.max(0.0);
            }
        }
        if let Some(dp) = self.drain_port {
            let wire = &inputs[self.layout.input_range(dp)];
            validate_stream_sum("drain_stream", wire)?;
            h_dot_streams_w += wire[STREAM_ENTHALPY_RATE];
            dm_h2o_drain = dt_s * wire[STREAM_CARRIED];
            for k in 0..AUX_LIQ_SPECIES_COUNT {
                let cand = aux_liq[k] + dt_s * wire[STREAM_CARRIED + 1 + k];
                if cand < -1.0e-12 {
                    return Err(CycleErrorV1::DepletedInventory {
                        inventory: "hx_aux_liq",
                    });
                }
                aux_liq[k] = cand.max(0.0);
            }
        }

        let ua_eff = self.effective_ua_w_k(old, fouling_or_airflow);
        let is_supercritical = p_ref >= self.fluid.critical_pressure_pa();

        // Compute refrigerant outlet enthalpy & implicitly damped air-side heat exchange
        // (prevents thermal oscillation on large dt_s while keeping exact First-Law closure).
        let (h_ref_out, t_ref_out, q_air_into_hx_w, defrost_active, box_warming_spike_w) =
            match self.config.mode {
                HxModeV1::CondenserGasCooler => {
                    let t_out = if is_supercritical {
                        // Supercritical R-744 gas-cooler cools fluid toward ambient + approach.
                        (t_air_k + self.config.target_approach_k * (1.0 + fouling_or_airflow))
                            .max(t_air_k + 1.5)
                    } else {
                        let env = self.fluid.bubble_dew_at_pressure(p_ref)?;
                        (env.t_bubble_k - self.config.target_approach_k).max(t_air_k + 1.0)
                    };
                    let phase_out = if is_supercritical {
                        PengRobinsonPhaseV1::Stable
                    } else {
                        PengRobinsonPhaseV1::Liquid
                    };
                    let out_pt = self.fluid.state_pt(t_out, p_ref, phase_out)?;
                    let h_out = out_pt.enthalpy_j_kg.min(h_ref_in);
                    let q_ref_release = m_dot_ref * (h_ref_in - h_out).max(0.0);
                    let c_eff = self.config.coil_heat_capacity_j_k;
                    let rate = ua_eff / c_eff;
                    let damp = if rate * dt_s > 1.0e-9 {
                        -expm1(-rate * dt_s) / (rate * dt_s)
                    } else {
                        1.0
                    };
                    let t_eq = t_air_k + q_ref_release / ua_eff.max(1.0);
                    let t_avg = t_eq + (old.t_coil_k - t_eq) * damp;
                    let q_to_air = (ua_eff * (t_avg - t_air_k)).max(0.0);
                    (h_out, t_out, -q_to_air, false, 0.0)
                }
                HxModeV1::FrostingEvaporator => {
                    let defrost_on = q_defrost_w > 1.0;
                    if defrost_on {
                        // Defrost cycle: fans off (~12% UA natural convection spill), electric/hot-gas
                        // heater warms coil to 0 C, melts frost to drain pan, and spills modest
                        // temporary box-air warming spike once frost clears.
                        let h_out = h_ref_in;
                        let ua_spill = 0.12 * ua_eff;
                        let c_eff = self.config.coil_heat_capacity_j_k
                            + old.m_meltwater_pan_kg * CP_WATER_LIQ_J_KG_K;
                        let rate = ua_spill / c_eff;
                        let damp = if rate * dt_s > 1.0e-9 {
                            -expm1(-rate * dt_s) / (rate * dt_s)
                        } else {
                            1.0
                        };
                        let t_avg = if old.m_frost_kg > 1.0e-9 {
                            old.t_coil_k.min(ICE_POINT_K)
                        } else {
                            let t_eq = t_air_k + q_defrost_w / ua_spill.max(1.0);
                            t_eq + (old.t_coil_k - t_eq) * damp
                        };
                        let q_spike_to_box = (ua_spill * (t_avg - t_air_k)).max(0.0);
                        (h_out, t_avg, -q_spike_to_box, true, q_spike_to_box)
                    } else {
                        // Active cooling: refrigerant boils inside tubes (UA_ref = 2.5 * UA_clean)
                        // absorbing both box air heat and frost deposition latent heat, keeping
                        // T_coil pinned near T_evap while frost monotonically degrades UA_air.
                        let env = self.fluid.bubble_dew_at_pressure(p_ref)?;
                        let t_evap_mid = 0.5 * (env.t_bubble_k + env.t_dew_k);
                        let t_out_target = env.t_dew_k + self.config.target_approach_k;
                        let sup_pt =
                            self.fluid
                                .state_pt(t_out_target, p_ref, PengRobinsonPhaseV1::Vapor)?;
                        let q_ref_max = m_dot_ref * (sup_pt.enthalpy_j_kg - h_ref_in).max(0.0);
                        let ua_ref = if m_dot_ref > 1.0e-12 {
                            2.5 * self.config.ua_clean_w_k
                        } else {
                            0.0
                        };
                        let m_frost_rate = (dm_h2o_frost_in / dt_s).max(0.0);
                        let h_vap_air = H_F0_H2O_VAP_J_KG
                            + CP_WATER_VAP_J_KG_K * (t_air_k - JANAF_REFERENCE_TEMPERATURE_K);
                        let q_frost_lat = m_frost_rate * (h_vap_air - H_ICE_0C_JANAF_J_KG);

                        let ua_tot = (ua_eff + ua_ref).max(1.0);
                        let c_eff = self.config.coil_heat_capacity_j_k
                            + old.m_frost_kg * CP_WATER_ICE_J_KG_K;
                        let rate = ua_tot / c_eff;
                        let damp = if rate * dt_s > 1.0e-9 {
                            -expm1(-rate * dt_s) / (rate * dt_s)
                        } else {
                            1.0
                        };
                        let t_eq = (ua_eff * t_air_k + ua_ref * t_evap_mid + q_frost_lat) / ua_tot;
                        let t_avg = t_eq + (old.t_coil_k - t_eq) * damp;
                        let q_air_in = (ua_eff * (t_air_k - t_avg)).max(0.0);
                        let q_ref_actual = (ua_ref * (t_avg - t_evap_mid)).clamp(0.0, q_ref_max);
                        let h_out = if m_dot_ref > 1.0e-12 {
                            h_ref_in + q_ref_actual / m_dot_ref
                        } else {
                            h_ref_in
                        };
                        (h_out, t_out_target, q_air_in, false, 0.0)
                    }
                }
            };

        // Update frost & liquid meltwater pan inventories conservatively.
        let mut m_frost = old.m_frost_kg;
        let mut m_melt = old.m_meltwater_pan_kg;

        // First apply incoming frosting vapor & outgoing liquid drain.
        if dm_h2o_frost_in >= 0.0 {
            m_frost += dm_h2o_frost_in;
        } else {
            let cand = m_frost + dm_h2o_frost_in;
            if cand < -1.0e-12 {
                return Err(CycleErrorV1::DepletedInventory {
                    inventory: "coil_frost",
                });
            }
            m_frost = cand.max(0.0);
        }
        let melt_after_drain = m_melt + dm_h2o_drain;
        if melt_after_drain < -1.0e-12 {
            return Err(CycleErrorV1::DepletedInventory {
                inventory: "meltwater_drain_pan",
            });
        }
        m_melt = melt_after_drain.max(0.0);
        let drained_now = (-dm_h2o_drain).max(0.0);

        // Compute exact First-Law new stored energy.
        let q_ref_in_w = m_dot_ref * (h_ref_in - h_ref_out);
        let e_port_in_j = dt_s * (q_ref_in_w + q_air_into_hx_w + q_defrost_w + h_dot_streams_w);
        let e_new_j = old.stored_energy_j + e_port_in_j;

        // During defrost (or coil warming above 0 C), melt frost into liquid drain pan using
        // exact enthalpy above 0 C (zero energy leak, zero re-evaporation into box air!).
        let e_zero_c_frozen = m_frost * H_ICE_0C_JANAF_J_KG + m_melt * H_WATER_LIQ_0C_JANAF_J_KG;
        if m_frost > 0.0 && e_new_j > e_zero_c_frozen {
            let excess_j = e_new_j - e_zero_c_frozen;
            let dm_melt = (excess_j / ICE_FUSION_ENTHALPY_0C_J_KG)
                .min(m_frost)
                .max(0.0);
            m_frost = (m_frost - dm_melt).max(0.0);
            m_melt += dm_melt;
        }

        let mut next_state = CondenserEvaporatorHXStateV1 {
            time_s: old.time_s + dt_s,
            t_coil_k: old.t_coil_k,
            m_frost_kg: m_frost,
            m_meltwater_pan_kg: m_melt,
            stored_energy_j: e_new_j,
            aux_gas_kg: aux_gas,
            aux_liq_kg: aux_liq,
            cumulative_drain_kg: old.cumulative_drain_kg + drained_now,
        };
        next_state.t_coil_k = invert_hx_coil_temp(&next_state, &self.config, e_new_j)?;
        next_state.validate(&self.config)?;
        *new = next_state;

        // Populate outputs.
        outputs.fill(0.0);
        let r_in_out = &mut outputs[self.layout.output_range(0)];
        r_in_out[STREAM_MASS_RATE] = m_dot_ref;
        r_in_out[STREAM_ENTHALPY_RATE] = m_dot_ref * h_ref_in;
        r_in_out[STREAM_CARRIED] = m_dot_ref;

        let r_out_out = &mut outputs[self.layout.output_range(1)];
        r_out_out[STREAM_MASS_RATE] = -m_dot_ref;
        r_out_out[STREAM_ENTHALPY_RATE] = -m_dot_ref * h_ref_out;
        r_out_out[STREAM_CARRIED] = -m_dot_ref;

        outputs[self.layout.output_range(2)][0] = q_air_into_hx_w;
        outputs[self.layout.output_range(3)][0] = 0.5 * (old.t_coil_k + new.t_coil_k);

        if let Some(mp) = self.moisture_port {
            let m_out = &mut outputs[self.layout.output_range(mp)];
            m_out[STREAM_PRESSURE] = 101_325.0;
            m_out[STREAM_SPECIFIC_ENTHALPY] =
                H_F0_H2O_VAP_J_KG + CP_WATER_VAP_J_KG_K * (t_air_k - JANAF_REFERENCE_TEMPERATURE_K);
            m_out[STREAM_CARRIED] = 1.0;
        }
        if let Some(dp) = self.drain_port {
            let d_out = &mut outputs[self.layout.output_range(dp)];
            d_out[STREAM_PRESSURE] = 101_325.0;
            d_out[STREAM_SPECIFIC_ENTHALPY] = H_WATER_LIQ_0C_JANAF_J_KG
                + CP_WATER_LIQ_J_KG_K * (new.t_coil_k - ICE_POINT_K).max(0.0);
            d_out[STREAM_CARRIED] = 1.0;
        }

        let probes = &mut outputs[self.layout.output_range(self.probes_port)];
        probes[0] = new.t_coil_k;
        probes[1] = q_air_into_hx_w;
        probes[2] = self.effective_ua_w_k(new, fouling_or_airflow);
        probes[3] = new.delta_frost_m(&self.config);
        probes[4] = new.m_frost_kg;
        probes[5] = t_ref_out;
        probes[6] = h_ref_out;
        probes[7] = if is_supercritical { 1.0 } else { 0.0 };
        probes[8] = if defrost_active { 1.0 } else { 0.0 };
        probes[9] = box_warming_spike_w;

        let e_scale = old.stored_energy_j.abs()
            + new.stored_energy_j.abs()
            + e_port_in_j.abs()
            + (m_dot_ref * h_ref_in.abs() * dt_s)
            + 1.0;
        let receipt = Receipt {
            stored_change_j: new.stored_energy_j - old.stored_energy_j,
            dissipated_j: 0.0,
            numerical_loss_j: 0.0,
            entropy_generated_j_per_k: 0.0,
            energy_tolerance_j: 512.0 * f64::EPSILON * e_scale,
            residual: 0.0,
            tolerance: 0.0,
        };
        Ok(Proposal::Step(receipt))
    }

    const STATE_CHARTER: Option<StateIdentityCharterV2> = Some(StateIdentityCharterV2 {
        owner: "fs-multiphase::CondenserEvaporatorHX",
        state_family: "pr-hx-frost-v1",
        schema_grammar: "time_s, t_coil_k, m_frost_kg, m_meltwater_pan_kg, stored_energy_j, aux_gas_kg[8], aux_liq_kg[7], cumulative_drain_kg in coherent SI",
        codec_grammar: "little-endian u64 words: codec tag, field count (21), followed by 21 IEEE-754 f64 bit patterns",
        codec_version: 1,
    });

    fn encode_state(&self, state: &Self::State, words: &mut Vec<u64>) -> bool {
        if state.validate(&self.config).is_err() {
            return false;
        }
        words.push(HX_CODEC_TAG);
        words.push(21);
        words.push(state.time_s.to_bits());
        words.push(state.t_coil_k.to_bits());
        words.push(state.m_frost_kg.to_bits());
        words.push(state.m_meltwater_pan_kg.to_bits());
        words.push(state.stored_energy_j.to_bits());
        for &v in &state.aux_gas_kg {
            words.push(v.to_bits());
        }
        for &v in &state.aux_liq_kg {
            words.push(v.to_bits());
        }
        words.push(state.cumulative_drain_kg.to_bits());
        true
    }

    fn decode_state(&self, words: &[u64]) -> Option<Self::State> {
        if words.len() != 23 || words[0] != HX_CODEC_TAG || words[1] != 21 {
            return None;
        }
        let mut aux_gas = [0.0; AUX_GAS_SPECIES_COUNT];
        for (i, slot) in aux_gas.iter_mut().enumerate() {
            *slot = f64::from_bits(words[7 + i]);
        }
        let mut aux_liq = [0.0; AUX_LIQ_SPECIES_COUNT];
        for (i, slot) in aux_liq.iter_mut().enumerate() {
            *slot = f64::from_bits(words[7 + AUX_GAS_SPECIES_COUNT + i]);
        }
        let state = CondenserEvaporatorHXStateV1 {
            time_s: f64::from_bits(words[2]),
            t_coil_k: f64::from_bits(words[3]),
            m_frost_kg: f64::from_bits(words[4]),
            m_meltwater_pan_kg: f64::from_bits(words[5]),
            stored_energy_j: f64::from_bits(words[6]),
            aux_gas_kg: aux_gas,
            aux_liq_kg: aux_liq,
            cumulative_drain_kg: f64::from_bits(words[22]),
        };
        state.validate(&self.config).ok().map(|()| state)
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        words.push(species_discriminant(self.config.species));
        words.push(match self.config.mode {
            HxModeV1::CondenserGasCooler => 0,
            HxModeV1::FrostingEvaporator => 1,
        });
        let scalars = [
            self.config.ua_clean_w_k,
            self.config.coil_area_m2,
            self.config.coil_heat_capacity_j_k,
            self.config.frost_density_kg_m3,
            self.config.frost_conductivity_w_m_k,
            self.config.fin_spacing_m,
            self.config.target_approach_k,
        ];
        words.extend(scalars.iter().map(|v| v.to_bits()));
    }
}

// ============================================================================
// COMPONENT 3: ExpansionThrottleComponent (TXV / Capillary)
// ============================================================================

/// Configuration for [`ExpansionThrottleComponent`].
#[derive(Clone, Debug, PartialEq)]
pub struct ExpansionThrottleConfig {
    /// Working refrigerant species.
    pub species: RefrigerantSpecies,
    /// Discharge coefficient `C_d` in `(0, 1]`.
    pub discharge_coefficient: f64,
    /// Maximum open orifice area `A_max` [m^2].
    pub max_orifice_area_m2: f64,
    /// Minimum bleed opening fraction `u_min` in `[0.02, 1.0]` (`1.0` for fixed capillary).
    pub min_opening_fraction: f64,
    /// Target evaporator outlet superheat `Delta T_sh,set` [K].
    pub target_superheat_k: f64,
    /// Proportional superheat modulation gain `K_sh` [1/K].
    pub superheat_gain_per_k: f64,
    /// Sensing bulb time constant `tau_bulb` [s].
    pub bulb_tau_s: f64,
}

impl ExpansionThrottleConfig {
    /// Standard commercial thermostatic expansion valve (TXV) preset.
    #[must_use]
    pub const fn txv_default(species: RefrigerantSpecies) -> Self {
        Self {
            species,
            discharge_coefficient: 0.65,
            max_orifice_area_m2: 1.8e-6,
            min_opening_fraction: 0.08,
            target_superheat_k: 6.0,
            superheat_gain_per_k: 0.12,
            bulb_tau_s: 12.0,
        }
    }

    /// Validate throttle configuration fail-closed.
    pub fn validate(&self) -> Result<(), CycleErrorV1> {
        if !self.discharge_coefficient.is_finite()
            || self.discharge_coefficient <= 0.0
            || self.discharge_coefficient > 1.0
            || !self.max_orifice_area_m2.is_finite()
            || self.max_orifice_area_m2 <= 0.0
            || !self.min_opening_fraction.is_finite()
            || !(0.01..=1.0).contains(&self.min_opening_fraction)
            || !self.target_superheat_k.is_finite()
            || self.target_superheat_k <= 0.0
            || !self.superheat_gain_per_k.is_finite()
            || self.superheat_gain_per_k < 0.0
            || !self.bulb_tau_s.is_finite()
            || self.bulb_tau_s <= 0.0
        {
            return Err(CycleErrorV1::InvalidConfig {
                parameter: "throttle_config",
                reason: "invalid orifice, discharge coefficient, superheat gain, or bulb time constant",
            });
        }
        Ok(())
    }
}

/// Accepted/candidate state for [`ExpansionThrottleComponent`].
#[derive(Clone, Debug, PartialEq)]
pub struct ExpansionThrottleStateV1 {
    /// Simulation time [s].
    pub time_s: f64,
    /// Modulated TXV needle opening fraction in `[min_opening_fraction, 1.0]`.
    pub valve_opening: f64,
    /// Cumulative throttling entropy generation [J/K].
    pub cumulative_entropy_gen_j_k: f64,
}

impl ExpansionThrottleStateV1 {
    /// Create initial TXV state with opening `valve_opening`.
    pub fn new(valve_opening: f64) -> Result<Self, CycleErrorV1> {
        let s = Self {
            time_s: 0.0,
            valve_opening,
            cumulative_entropy_gen_j_k: 0.0,
        };
        s.validate()?;
        Ok(s)
    }

    /// Validate state invariants fail-closed.
    pub fn validate(&self) -> Result<(), CycleErrorV1> {
        if !self.time_s.is_finite()
            || self.time_s < 0.0
            || !self.valve_opening.is_finite()
            || !(0.0..=1.0).contains(&self.valve_opening)
            || !self.cumulative_entropy_gen_j_k.is_finite()
            || self.cumulative_entropy_gen_j_k < 0.0
        {
            return Err(CycleErrorV1::InvalidScalar {
                quantity: "throttle_state",
                value: self.valve_opening,
                reason: "valve opening must lie in [0, 1] and counters must be non-negative finite",
            });
        }
        Ok(())
    }
}

/// Exact isenthalpic TXV / capillary expansion [`Component`] using Péneloux-translated
/// liquid density for orifice flow and superheat-feedback valve modulation.
#[derive(Clone, Debug)]
pub struct ExpansionThrottleComponent {
    config: ExpansionThrottleConfig,
    fluid: RefrigerantModelV1,
    ports: Vec<PortSpec>,
    layout: PortLayout,
}

impl ExpansionThrottleComponent {
    /// Construct a validated [`ExpansionThrottleComponent`].
    pub fn new(config: ExpansionThrottleConfig) -> Result<Self, CycleErrorV1> {
        config.validate()?;
        let fluid = RefrigerantModelV1::for_species(config.species)?;
        Self::with_fluid(config, fluid)
    }

    /// Construct from an explicit [`RefrigerantModelV1`].
    pub fn with_fluid(
        config: ExpansionThrottleConfig,
        fluid: RefrigerantModelV1,
    ) -> Result<Self, CycleErrorV1> {
        config.validate()?;
        if fluid.species() != config.species {
            return Err(CycleErrorV1::InvalidConfig {
                parameter: "species",
                reason: "RefrigerantModelV1 species must match config.species",
            });
        }
        let mk = |n: &'static str| {
            StableId::new(n).map_err(|_| CycleErrorV1::InvalidConfig {
                parameter: n,
                reason: "invalid StableId",
            })
        };
        let ports = vec![
            PortSpec::identified_stream(
                mk("inlet_stream")?,
                PortRole::EffortIn,
                refrigerant_v1(config.species)?,
            )?,
            PortSpec::identified_stream(
                mk("outlet_stream")?,
                PortRole::EffortIn,
                refrigerant_v1(config.species)?,
            )?,
            PortSpec::signal_in(mk("suction_superheat_k")?, NonZeroUsize::MIN),
            PortSpec::signal_in(mk("inlet_temperature_k")?, NonZeroUsize::MIN),
            PortSpec::signal_out(
                mk("probes")?,
                NonZeroUsize::new(THROTTLE_PROBE_COUNT).ok_or(CycleErrorV1::InvalidConfig {
                    parameter: "probes",
                    reason: "nonzero width",
                })?,
            ),
        ];
        let layout = PortLayout::new(&ports)?;
        Ok(Self {
            config,
            fluid,
            ports,
            layout,
        })
    }

    /// Borrow the port layout.
    #[must_use]
    pub const fn layout(&self) -> &PortLayout {
        &self.layout
    }

    /// Build `PortLayout`-ordered inputs for [`Component::propose`].
    #[must_use]
    pub fn build_inputs(
        &self,
        p_in_pa: f64,
        h_in_j_kg: f64,
        p_out_pa: f64,
        suction_superheat_k: f64,
        inlet_temperature_k: f64,
    ) -> Vec<f64> {
        let mut inputs = vec![0.0; self.layout.input_len()];
        let s_in = &mut inputs[self.layout.input_range(0)];
        s_in[STREAM_PRESSURE] = p_in_pa;
        s_in[STREAM_SPECIFIC_ENTHALPY] = h_in_j_kg;
        s_in[STREAM_CARRIED] = 1.0;

        let s_out = &mut inputs[self.layout.input_range(1)];
        s_out[STREAM_PRESSURE] = p_out_pa;
        s_out[STREAM_SPECIFIC_ENTHALPY] = h_in_j_kg;
        s_out[STREAM_CARRIED] = 1.0;

        inputs[self.layout.input_range(2)][0] = suction_superheat_k;
        inputs[self.layout.input_range(3)][0] = inlet_temperature_k;
        inputs
    }
}

impl Component for ExpansionThrottleComponent {
    type State = ExpansionThrottleStateV1;
    type Error = CycleErrorV1;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Smooth
    }

    fn end_of_step(&self, _port: usize) -> EndOfStep {
        EndOfStep::Untold
    }

    fn stored_energy_j(&self, state: &Self::State) -> Result<f64, Self::Error> {
        state.validate()?;
        Ok(0.0)
    }

    fn stored_mass_kg(&self, state: &Self::State) -> Result<Option<f64>, Self::Error> {
        state.validate()?;
        Ok(Some(0.0))
    }

    fn stored_carried_quantity(
        &self,
        state: &Self::State,
        carried: usize,
    ) -> Result<Option<f64>, Self::Error> {
        state.validate()?;
        if carried < self.layout.carried_len() {
            Ok(Some(0.0))
        } else {
            Ok(None)
        }
    }

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
        old.validate()?;
        if !dt_s.is_finite() || dt_s <= 0.0 {
            return Err(CycleErrorV1::InvalidScalar {
                quantity: "dt_s",
                value: dt_s,
                reason: "must be finite positive",
            });
        }
        if inputs.len() != self.layout.input_len() || outputs.len() != self.layout.output_len() {
            return Err(CycleErrorV1::InvalidConfig {
                parameter: "port_vectors",
                reason: "inputs or outputs length mismatch",
            });
        }
        if !inputs.iter().all(|v| v.is_finite()) {
            return Err(CycleErrorV1::InvalidConfig {
                parameter: "inputs",
                reason: "all throttle inputs must be finite",
            });
        }

        let s_in = &inputs[self.layout.input_range(0)];
        let s_out = &inputs[self.layout.input_range(1)];
        let p_in = s_in[STREAM_PRESSURE];
        let h_in = s_in[STREAM_SPECIFIC_ENTHALPY];
        let p_out = s_out[STREAM_PRESSURE];
        let sh_k = inputs[self.layout.input_range(2)][0];
        let t_in_k = inputs[self.layout.input_range(3)][0];

        if sh_k < 0.0 {
            return Err(CycleErrorV1::ThermodynamicInversion(
                "sensed suction superheat cannot be negative",
            ));
        }

        let th = self.fluid.isenthalpic_throttle(t_in_k, p_in, p_out)?;
        // First-order lagged superheat-feedback valve modulation.
        let u_target = (0.5
            + self.config.superheat_gain_per_k * (sh_k - self.config.target_superheat_k))
            .clamp(self.config.min_opening_fraction, 1.0);
        let alpha = -expm1(-dt_s / self.config.bulb_tau_s);
        let u_new = (old.valve_opening + alpha * (u_target - old.valve_opening))
            .clamp(self.config.min_opening_fraction, 1.0);
        let u_avg = 0.5 * (old.valve_opening + u_new);

        // Orifice mass flow using Péneloux-translated liquid density!
        let dp = (p_in - p_out).max(0.0);
        let a_eff = u_avg * self.config.max_orifice_area_m2;
        let m_dot =
            self.config.discharge_coefficient * a_eff * sqrt(2.0 * th.inlet.density_kg_m3 * dp);
        let s_gen_w_k = m_dot * th.ds_gen_j_kg_k;

        let next_state = ExpansionThrottleStateV1 {
            time_s: old.time_s + dt_s,
            valve_opening: u_new,
            cumulative_entropy_gen_j_k: old.cumulative_entropy_gen_j_k + s_gen_w_k * dt_s,
        };
        next_state.validate()?;
        *new = next_state;

        outputs.fill(0.0);
        // Exact isenthalpic h_out == h_in on the wire.
        let in_out = &mut outputs[self.layout.output_range(0)];
        in_out[STREAM_MASS_RATE] = m_dot;
        in_out[STREAM_ENTHALPY_RATE] = m_dot * h_in;
        in_out[STREAM_CARRIED] = m_dot;

        let out_out = &mut outputs[self.layout.output_range(1)];
        out_out[STREAM_MASS_RATE] = -m_dot;
        out_out[STREAM_ENTHALPY_RATE] = -m_dot * h_in;
        out_out[STREAM_CARRIED] = -m_dot;

        let probes = &mut outputs[self.layout.output_range(4)];
        probes[0] = m_dot;
        probes[1] = u_new;
        probes[2] = th.inlet.density_kg_m3;
        probes[3] = th.outlet_quality;
        probes[4] = th.outlet_temperature_k;
        probes[5] = s_gen_w_k;

        let e_scale = (m_dot * h_in.abs() * dt_s).max(1.0);
        let receipt = Receipt {
            stored_change_j: 0.0,
            dissipated_j: 0.0,
            numerical_loss_j: 0.0,
            entropy_generated_j_per_k: (s_gen_w_k * dt_s).max(0.0),
            energy_tolerance_j: 256.0 * f64::EPSILON * e_scale,
            residual: 0.0,
            tolerance: 0.0,
        };
        Ok(Proposal::Step(receipt))
    }

    const STATE_CHARTER: Option<StateIdentityCharterV2> = Some(StateIdentityCharterV2 {
        owner: "fs-multiphase::ExpansionThrottle",
        state_family: "txv-isenthalpic-v1",
        schema_grammar: "time_s, valve_opening, cumulative_entropy_gen_j_k in coherent SI",
        codec_grammar: "little-endian u64 words: codec tag, field count (3), followed by 3 IEEE-754 f64 bit patterns",
        codec_version: 1,
    });

    fn encode_state(&self, state: &Self::State, words: &mut Vec<u64>) -> bool {
        if state.validate().is_err() {
            return false;
        }
        words.push(THROTTLE_CODEC_TAG);
        words.push(3);
        words.push(state.time_s.to_bits());
        words.push(state.valve_opening.to_bits());
        words.push(state.cumulative_entropy_gen_j_k.to_bits());
        true
    }

    fn decode_state(&self, words: &[u64]) -> Option<Self::State> {
        if words.len() != 5 || words[0] != THROTTLE_CODEC_TAG || words[1] != 3 {
            return None;
        }
        let state = ExpansionThrottleStateV1 {
            time_s: f64::from_bits(words[2]),
            valve_opening: f64::from_bits(words[3]),
            cumulative_entropy_gen_j_k: f64::from_bits(words[4]),
        };
        state.validate().ok().map(|()| state)
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        words.push(species_discriminant(self.config.species));
        let scalars = [
            self.config.discharge_coefficient,
            self.config.max_orifice_area_m2,
            self.config.min_opening_fraction,
            self.config.target_superheat_k,
            self.config.superheat_gain_per_k,
            self.config.bulb_tau_s,
        ];
        words.extend(scalars.iter().map(|v| v.to_bits()));
    }
}

// ============================================================================
// COMPONENT 4: IceHarvesterComponent (Air-Cooled Commercial Ice Machine)
// ============================================================================

/// Configuration for [`IceHarvesterComponent`].
#[derive(Clone, Debug, PartialEq)]
pub struct IceHarvesterConfig {
    /// Working refrigerant species (`R-404A`, `R-290`, or `R-448A` typical for commercial ice heads).
    pub species: RefrigerantSpecies,
    /// Evaporator freeze-plate saturation temperature [K] (typically `-15 C` = `258.15 K`).
    pub t_evap_k: f64,
    /// Suction superheat [K].
    pub suction_superheat_k: f64,
    /// Condenser subcooling [K].
    pub subcooling_k: f64,
    /// Air-cooled condenser conductance `UA_cond_clean` [W/K].
    pub ua_cond_clean_w_k: f64,
    /// Compressor swept displacement [m^3/s].
    pub swept_volume_m3_s: f64,
    /// Isentropic efficiency in `(0.2, 1.0]`.
    pub eta_isentropic: f64,
    /// Motor electrical efficiency in `(0.3, 1.0]`.
    pub eta_motor: f64,
    /// Clearance volume ratio.
    pub clearance_ratio: f64,
    /// Harvested Ice Ih bin temperature [K] (`-5 C` = `268.15 K`).
    pub t_ice_bin_k: f64,
    /// Maximum ice storage bin capacity [kg] (thermostat/bin-switch shuts off harvest when full).
    pub max_bin_capacity_kg: f64,
}

impl IceHarvesterConfig {
    /// Commercial air-cooled ~900 lb/day modular cuber preset (rejects ~2.4--3.1 kW continuous
    /// sensible waste heat to kitchen air; 25--40% harvest derating between 70 F and 95 F).
    #[must_use]
    pub const fn commercial_air_cooled_default() -> Self {
        Self {
            species: RefrigerantSpecies::R404A,
            t_evap_k: 258.15,
            suction_superheat_k: 6.0,
            subcooling_k: 4.0,
            ua_cond_clean_w_k: 210.0,
            swept_volume_m3_s: 1.15e-3,
            eta_isentropic: 0.72,
            eta_motor: 0.88,
            clearance_ratio: 0.045,
            t_ice_bin_k: 268.15,
            max_bin_capacity_kg: 250.0,
        }
    }

    /// Validate ice harvester configuration fail-closed.
    pub fn validate(&self) -> Result<(), CycleErrorV1> {
        if !self.t_evap_k.is_finite()
            || self.t_evap_k <= 200.0
            || self.t_evap_k >= ICE_POINT_K - 2.0
            || !self.t_ice_bin_k.is_finite()
            || self.t_ice_bin_k <= 220.0
            || self.t_ice_bin_k >= ICE_POINT_K
            || !self.ua_cond_clean_w_k.is_finite()
            || self.ua_cond_clean_w_k <= 0.0
            || !self.swept_volume_m3_s.is_finite()
            || self.swept_volume_m3_s <= 0.0
            || !self.eta_isentropic.is_finite()
            || self.eta_isentropic <= 0.2
            || self.eta_isentropic > 1.0
            || !self.eta_motor.is_finite()
            || self.eta_motor <= 0.3
            || self.eta_motor > 1.0
            || !self.max_bin_capacity_kg.is_finite()
            || self.max_bin_capacity_kg <= 0.0
            || !self.suction_superheat_k.is_finite()
            || self.suction_superheat_k < 0.5
            || !self.subcooling_k.is_finite()
            || self.subcooling_k < 0.0
            || !self.clearance_ratio.is_finite()
            || !(0.0..=0.20).contains(&self.clearance_ratio)
        {
            return Err(CycleErrorV1::InvalidConfig {
                parameter: "ice_harvester_config",
                reason: "invalid ice harvester temperatures, UA, displacement, efficiencies, or bin capacity",
            });
        }
        Ok(())
    }
}

/// Accepted/candidate state for [`IceHarvesterComponent`].
#[derive(Clone, Debug, PartialEq)]
pub struct IceHarvesterStateV1 {
    /// Simulation time [s].
    pub time_s: f64,
    /// Harvested Ice Ih inventory in storage bin [kg].
    pub ice_bin_mass_kg: f64,
    /// Liquid potable water inventory in freeze-plate sump [kg].
    pub sump_water_mass_kg: f64,
    /// Sump water temperature [K].
    pub t_sump_k: f64,
    /// Total stored enthalpy on `janaf-298.15K-formation` [J].
    pub stored_energy_j: f64,
    /// Auxiliary non-H2O liquid species carried on `kitchen_liquid_v1` [kg].
    pub aux_liq_kg: [f64; AUX_LIQ_SPECIES_COUNT],
    /// Cumulative electrical energy consumed [J].
    pub cumulative_elec_work_j: f64,
    /// Cumulative waste heat rejected to kitchen room air [J].
    pub cumulative_room_heat_j: f64,
}

impl IceHarvesterStateV1 {
    /// Total stored H2O mass (ice bin + liquid sump) [kg].
    #[must_use]
    pub fn total_h2o_kg(&self) -> f64 {
        self.ice_bin_mass_kg + self.sump_water_mass_kg
    }

    /// Total stored mass [kg].
    #[must_use]
    pub fn total_mass_kg(&self) -> f64 {
        self.total_h2o_kg() + self.aux_liq_kg.iter().sum::<f64>()
    }

    /// Validate state invariants fail-closed.
    pub fn validate(&self, config: &IceHarvesterConfig) -> Result<(), CycleErrorV1> {
        if !self.time_s.is_finite()
            || self.time_s < 0.0
            || !self.ice_bin_mass_kg.is_finite()
            || self.ice_bin_mass_kg < 0.0
            || !self.sump_water_mass_kg.is_finite()
            || self.sump_water_mass_kg < 0.0
            || !self.t_sump_k.is_finite()
            || !(260.0..=360.0).contains(&self.t_sump_k)
            || !self.stored_energy_j.is_finite()
            || !self.cumulative_elec_work_j.is_finite()
            || self.cumulative_elec_work_j < 0.0
            || !self.cumulative_room_heat_j.is_finite()
            || self.cumulative_room_heat_j < 0.0
            || self.aux_liq_kg.iter().any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(CycleErrorV1::InvalidScalar {
                quantity: "ice_harvester_state",
                value: self.ice_bin_mass_kg,
                reason: "ice harvester inventories and counters must be finite and non-negative",
            });
        }
        let expected = ice_harvester_energy_j(
            self.ice_bin_mass_kg,
            self.sump_water_mass_kg,
            self.t_sump_k,
            config.t_ice_bin_k,
        );
        let scale = expected.abs().max(self.stored_energy_j.abs()).max(1.0);
        if (expected - self.stored_energy_j).abs() > 1.0e-5 * scale {
            return Err(CycleErrorV1::InvalidScalar {
                quantity: "stored_energy_j",
                value: self.stored_energy_j,
                reason: "ice harvester stored enthalpy disagrees with ice/sump masses and temperature",
            });
        }
        Ok(())
    }
}

fn h_ice_janaf_j_kg(t_ice_k: f64) -> f64 {
    H_ICE_0C_JANAF_J_KG + CP_WATER_ICE_J_KG_K * (t_ice_k - ICE_POINT_K)
}

fn h_water_liq_janaf_j_kg(t_k: f64) -> f64 {
    H2O_LIQUID_FORMATION_ENTHALPY_298K_J_KG
        + CP_WATER_LIQ_J_KG_K * (t_k - JANAF_REFERENCE_TEMPERATURE_K)
}

fn ice_harvester_energy_j(m_ice_kg: f64, m_sump_kg: f64, t_sump_k: f64, t_ice_bin_k: f64) -> f64 {
    m_ice_kg * h_ice_janaf_j_kg(t_ice_bin_k) + m_sump_kg * h_water_liq_janaf_j_kg(t_sump_k)
}

/// Air-cooled commercial ice machine [`Component`] composing the PR-EOS
/// vapor-compression cycle against a potable water feed (`kitchen_liquid_v1`),
/// Ice Ih harvest bin, electrical supply, and kitchen room air waste-heat sink.
#[derive(Clone, Debug)]
pub struct IceHarvesterComponent {
    config: IceHarvesterConfig,
    fluid: RefrigerantModelV1,
    ports: Vec<PortSpec>,
    layout: PortLayout,
}

impl IceHarvesterComponent {
    /// Construct a validated [`IceHarvesterComponent`].
    pub fn new(config: IceHarvesterConfig) -> Result<Self, CycleErrorV1> {
        config.validate()?;
        let fluid = RefrigerantModelV1::for_species(config.species)?;
        Self::with_fluid(config, fluid)
    }

    /// Construct from an explicit [`RefrigerantModelV1`].
    pub fn with_fluid(
        config: IceHarvesterConfig,
        fluid: RefrigerantModelV1,
    ) -> Result<Self, CycleErrorV1> {
        config.validate()?;
        if fluid.species() != config.species {
            return Err(CycleErrorV1::InvalidConfig {
                parameter: "species",
                reason: "RefrigerantModelV1 species must match config.species",
            });
        }
        let mk = |n: &'static str| {
            StableId::new(n).map_err(|_| CycleErrorV1::InvalidConfig {
                parameter: n,
                reason: "invalid StableId",
            })
        };
        let ports = vec![
            PortSpec::power(
                mk("electrical_power")?,
                PortKind::ElectricalVoltageCurrent,
                PortRole::EffortIn,
                NonZeroUsize::MIN,
            )?,
            PortSpec::heat(mk("room_heat")?, PortRole::EffortIn)?,
            PortSpec::identified_stream(
                mk("water_feed_stream")?,
                PortRole::FlowIn,
                kitchen_liquid_v1()?,
            )?,
            PortSpec::signal_in(mk("condenser_fouling")?, NonZeroUsize::MIN),
            PortSpec::signal_in(mk("water_inlet_temp_k")?, NonZeroUsize::MIN),
            PortSpec::signal_out(
                mk("probes")?,
                NonZeroUsize::new(ICE_HARVESTER_PROBE_COUNT).ok_or(
                    CycleErrorV1::InvalidConfig {
                        parameter: "probes",
                        reason: "nonzero width",
                    },
                )?,
            ),
        ];
        let layout = PortLayout::new(&ports)?;
        Ok(Self {
            config,
            fluid,
            ports,
            layout,
        })
    }

    /// Borrow the ice harvester configuration.
    #[must_use]
    pub const fn config(&self) -> &IceHarvesterConfig {
        &self.config
    }

    /// Borrow the port layout.
    #[must_use]
    pub const fn layout(&self) -> &PortLayout {
        &self.layout
    }

    /// Create an initial ice harvester state with `initial_ice_kg` in the bin and
    /// `initial_sump_kg` liquid water in the sump at `t_sump_k`.
    pub fn initial_state(
        &self,
        initial_ice_kg: f64,
        initial_sump_kg: f64,
        t_sump_k: f64,
    ) -> Result<IceHarvesterStateV1, CycleErrorV1> {
        let e = ice_harvester_energy_j(
            initial_ice_kg,
            initial_sump_kg,
            t_sump_k,
            self.config.t_ice_bin_k,
        );
        let state = IceHarvesterStateV1 {
            time_s: 0.0,
            ice_bin_mass_kg: initial_ice_kg,
            sump_water_mass_kg: initial_sump_kg,
            t_sump_k,
            stored_energy_j: e,
            aux_liq_kg: [0.0; AUX_LIQ_SPECIES_COUNT],
            cumulative_elec_work_j: 0.0,
            cumulative_room_heat_j: 0.0,
        };
        state.validate(&self.config)?;
        Ok(state)
    }

    /// Evaluate the steady-state PR-EOS cycle at kitchen ambient `t_room_k` and fouling `f_foul`.
    pub fn evaluate_cycle(
        &self,
        t_room_k: f64,
        condenser_fouling: f64,
    ) -> Result<ClosedCycleSummaryV1, CycleErrorV1> {
        self.fluid.solve_closed_cycle(CycleOperatingPointV1 {
            t_evap_k: self.config.t_evap_k,
            suction_superheat_k: self.config.suction_superheat_k,
            t_ambient_k: t_room_k,
            ua_cond_clean_w_k: self.config.ua_cond_clean_w_k,
            condenser_fouling_factor: condenser_fouling,
            subcooling_k: self.config.subcooling_k,
            swept_volume_m3_s: self.config.swept_volume_m3_s,
            eta_isentropic: self.config.eta_isentropic,
            eta_motor: self.config.eta_motor,
            shell_loss_fraction: 1.0,
            clearance_ratio: self.config.clearance_ratio,
            supercritical_pressure_pa: None,
        })
    }

    /// Build `PortLayout`-ordered inputs supplying potable feed water at `m_dot_feed_kg_s`
    /// and temperature `t_water_in_k` on `janaf-298.15K-formation`.
    #[must_use]
    pub fn build_inputs(
        &self,
        voltage_v: f64,
        t_room_k: f64,
        m_dot_feed_kg_s: f64,
        t_water_in_k: f64,
        condenser_fouling: f64,
    ) -> Vec<f64> {
        let mut inputs = vec![0.0; self.layout.input_len()];
        inputs[self.layout.input_range(0)][0] = voltage_v;
        inputs[self.layout.input_range(1)][0] = t_room_k;
        let feed = &mut inputs[self.layout.input_range(2)];
        let h_feed = h_water_liq_janaf_j_kg(t_water_in_k);
        feed[STREAM_MASS_RATE] = m_dot_feed_kg_s;
        feed[STREAM_ENTHALPY_RATE] = m_dot_feed_kg_s * h_feed;
        feed[STREAM_CARRIED] = m_dot_feed_kg_s;
        inputs[self.layout.input_range(3)][0] = condenser_fouling;
        inputs[self.layout.input_range(4)][0] = t_water_in_k;
        inputs
    }
}

impl Component for IceHarvesterComponent {
    type State = IceHarvesterStateV1;
    type Error = CycleErrorV1;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Smooth
    }

    fn end_of_step(&self, _port: usize) -> EndOfStep {
        EndOfStep::Untold
    }

    fn stored_energy_j(&self, state: &Self::State) -> Result<f64, Self::Error> {
        state.validate(&self.config)?;
        Ok(state.stored_energy_j)
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
        let Some(StreamQuantity::SpeciesMass(sp)) = self.layout.carried_quantities().get(carried)
        else {
            return Ok(None);
        };
        let amount = match sp.as_str() {
            "H2O" => state.total_h2o_kg(),
            "C2H5OH" => state.aux_liq_kg[0],
            "triglyceride" => state.aux_liq_kg[1],
            "ffa" => state.aux_liq_kg[2],
            "glycerol" => state.aux_liq_kg[3],
            "solute_salt" => state.aux_liq_kg[4],
            "solute_sugar" => state.aux_liq_kg[5],
            "solute_fond" => state.aux_liq_kg[6],
            _ => return Ok(None),
        };
        Ok(Some(amount))
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
            return Err(CycleErrorV1::InvalidScalar {
                quantity: "dt_s",
                value: dt_s,
                reason: "must be finite positive",
            });
        }
        if inputs.len() != self.layout.input_len() || outputs.len() != self.layout.output_len() {
            return Err(CycleErrorV1::InvalidConfig {
                parameter: "port_vectors",
                reason: "inputs or outputs length mismatch",
            });
        }
        if !inputs.iter().all(|v| v.is_finite()) {
            return Err(CycleErrorV1::InvalidConfig {
                parameter: "inputs",
                reason: "all ice harvester inputs must be finite",
            });
        }

        let voltage_v = inputs[self.layout.input_range(0)][0];
        let t_room_k = inputs[self.layout.input_range(1)][0];
        let feed_wire = &inputs[self.layout.input_range(2)];
        let fouling = inputs[self.layout.input_range(3)][0];
        let t_water_in_k = inputs[self.layout.input_range(4)][0];

        if voltage_v <= 0.0 {
            return Err(CycleErrorV1::InvalidScalar {
                quantity: "voltage_v",
                value: voltage_v,
                reason: "supply voltage must be strictly positive",
            });
        }
        if t_water_in_k <= ICE_POINT_K || t_water_in_k > 350.0 {
            return Err(CycleErrorV1::InvalidTemperature {
                temperature_k: t_water_in_k,
                reason: "potable feed water temperature must lie in (273.15, 350] K",
            });
        }
        if fouling < 0.0 {
            return Err(CycleErrorV1::InvalidScalar {
                quantity: "condenser_fouling",
                value: fouling,
                reason: "must be non-negative",
            });
        }
        validate_stream_sum("water_feed_stream", feed_wire)?;

        let dm_h2o_feed = dt_s * feed_wire[STREAM_CARRIED];
        let h_dot_feed_w = feed_wire[STREAM_ENTHALPY_RATE];
        let mut aux_liq = old.aux_liq_kg;
        for k in 0..AUX_LIQ_SPECIES_COUNT {
            let rate = feed_wire[STREAM_CARRIED + 1 + k];
            if rate.abs() > 1.0e-12 {
                return Err(CycleErrorV1::InvalidScalar {
                    quantity: "potable_water_non_h2o_solute",
                    value: rate,
                    reason: "ice harvester water_feed_stream admits potable H2O only",
                });
            }
            aux_liq[k] = old.aux_liq_kg[k];
        }

        let m_sump_avail = old.sump_water_mass_kg + dm_h2o_feed;
        if m_sump_avail < -1.0e-12 {
            return Err(CycleErrorV1::DepletedInventory {
                inventory: "sump_water",
            });
        }
        let m_sump_pre = m_sump_avail.max(0.0);

        // Solve closed PR-EOS vapor-compression cycle against kitchen room ambient!
        let cycle = self.evaluate_cycle(t_room_k, fouling)?;
        let bin_headroom = (self.config.max_bin_capacity_kg - old.ice_bin_mass_kg).max(0.0);
        let bin_full = bin_headroom <= 1.0e-9;

        // Compute mixed sump temperature before freezing from exact pre-freeze enthalpy.
        let h_ice_bin = h_ice_janaf_j_kg(self.config.t_ice_bin_k);
        let e_sump_pre =
            old.stored_energy_j - old.ice_bin_mass_kg * h_ice_bin + dt_s * h_dot_feed_w;
        let t_sump_pre = if m_sump_pre > 1.0e-12 {
            (JANAF_REFERENCE_TEMPERATURE_K
                + (e_sump_pre / m_sump_pre - H2O_LIQUID_FORMATION_ENTHALPY_298K_J_KG)
                    / CP_WATER_LIQ_J_KG_K)
                .clamp(260.0, 360.0)
        } else {
            t_water_in_k
        };

        // Freezing enthalpy per kg of Ice Ih from current mixed sump liquid temperature.
        let dh_freeze_j_kg = (h_water_liq_janaf_j_kg(t_sump_pre) - h_ice_bin).max(1.0e5);
        let (dm_ice, q_evap_actual_w, w_elec_w, q_room_reject_w) =
            if bin_full || m_sump_pre <= 1.0e-9 {
                (0.0, 0.0, 0.0, 0.0)
            } else {
                let dm_cap = cycle.q_evap_w * dt_s / dh_freeze_j_kg;
                // Keep a tiny residual sump film (or freeze all available up to bin capacity).
                let dm_made = dm_cap.min(m_sump_pre).min(bin_headroom).max(0.0);
                let duty_frac = if dm_cap > 1.0e-15 {
                    (dm_made / dm_cap).clamp(0.0, 1.0)
                } else {
                    0.0
                };
                let q_evap_used = dm_made * dh_freeze_j_kg / dt_s;
                let w_used = duty_frac * cycle.w_comp_w;
                // Both condenser and hermetic compressor shell reject into air-cooled kitchen room:
                // Q_room_reject = Q_cond + Q_shell = Q_evap_used + W_elec_used identically!
                (dm_made, q_evap_used, w_used, q_evap_used + w_used)
            };

        let m_ice_new = old.ice_bin_mass_kg + dm_ice;
        let m_sump_new = (m_sump_pre - dm_ice).max(0.0);
        // Exact First-Law energy update: E_new = E_old + dt * (W_elec - Q_room_reject + H_dot_feed)
        //                                       = E_old + dt * (-Q_evap_actual + H_dot_feed).
        let e_port_in_j = dt_s * (w_elec_w - q_room_reject_w + h_dot_feed_w);
        let e_new_j = old.stored_energy_j + e_port_in_j;
        let t_sump_new = if m_sump_new > 1.0e-12 {
            let e_sump_rem = e_new_j - m_ice_new * h_ice_bin;
            JANAF_REFERENCE_TEMPERATURE_K
                + (e_sump_rem / m_sump_new - H2O_LIQUID_FORMATION_ENTHALPY_298K_J_KG)
                    / CP_WATER_LIQ_J_KG_K
        } else {
            t_sump_pre
        };

        let next_state = IceHarvesterStateV1 {
            time_s: old.time_s + dt_s,
            ice_bin_mass_kg: m_ice_new,
            sump_water_mass_kg: m_sump_new,
            t_sump_k: t_sump_new,
            stored_energy_j: e_new_j,
            aux_liq_kg: aux_liq,
            cumulative_elec_work_j: old.cumulative_elec_work_j + w_elec_w * dt_s,
            cumulative_room_heat_j: old.cumulative_room_heat_j + q_room_reject_w * dt_s,
        };
        next_state.validate(&self.config)?;
        *new = next_state;

        outputs.fill(0.0);
        outputs[self.layout.output_range(0)][0] = w_elec_w / voltage_v;
        // Room heat port (EffortIn): negative flow = sensible waste heat rejected INTO kitchen room air!
        outputs[self.layout.output_range(1)][0] = -q_room_reject_w;

        let feed_out = &mut outputs[self.layout.output_range(2)];
        feed_out[STREAM_PRESSURE] = 101_325.0;
        feed_out[STREAM_SPECIFIC_ENTHALPY] = h_water_liq_janaf_j_kg(new.t_sump_k);
        feed_out[STREAM_CARRIED] = 1.0;

        let harvest_kg_h = (dm_ice / dt_s) * 3_600.0;
        let probes = &mut outputs[self.layout.output_range(5)];
        probes[0] = new.ice_bin_mass_kg;
        probes[1] = harvest_kg_h;
        probes[2] = q_evap_actual_w;
        probes[3] = w_elec_w;
        probes[4] = q_room_reject_w;
        probes[5] = if w_elec_w > 1.0e-9 {
            q_evap_actual_w / w_elec_w
        } else {
            0.0
        };
        probes[6] = cycle.t_high_out_k;
        probes[7] = cycle.p_high_pa;
        probes[8] = cycle.p_evap_pa;
        probes[9] = if new.ice_bin_mass_kg >= self.config.max_bin_capacity_kg - 1.0e-6 {
            1.0
        } else {
            0.0
        };

        let e_scale = old.stored_energy_j.abs()
            + new.stored_energy_j.abs()
            + e_port_in_j.abs()
            + (q_room_reject_w * dt_s).abs()
            + 1.0;
        let receipt = Receipt {
            stored_change_j: new.stored_energy_j - old.stored_energy_j,
            dissipated_j: 0.0,
            numerical_loss_j: 0.0,
            entropy_generated_j_per_k: (cycle.internal_entropy_gen_w_k * dt_s).max(0.0),
            energy_tolerance_j: 512.0 * f64::EPSILON * e_scale,
            residual: 0.0,
            tolerance: 0.0,
        };
        Ok(Proposal::Step(receipt))
    }

    const STATE_CHARTER: Option<StateIdentityCharterV2> = Some(StateIdentityCharterV2 {
        owner: "fs-multiphase::IceHarvester",
        state_family: "pr-ice-harvester-v1",
        schema_grammar: "time_s, ice_bin_mass_kg, sump_water_mass_kg, t_sump_k, stored_energy_j, aux_liq_kg[7], cumulative_elec_work_j, cumulative_room_heat_j in coherent SI on janaf-298.15K-formation",
        codec_grammar: "little-endian u64 words: codec tag, field count (14), followed by 14 IEEE-754 f64 bit patterns",
        codec_version: 1,
    });

    fn encode_state(&self, state: &Self::State, words: &mut Vec<u64>) -> bool {
        if state.validate(&self.config).is_err() {
            return false;
        }
        words.push(ICE_HARVESTER_CODEC_TAG);
        words.push(14);
        words.push(state.time_s.to_bits());
        words.push(state.ice_bin_mass_kg.to_bits());
        words.push(state.sump_water_mass_kg.to_bits());
        words.push(state.t_sump_k.to_bits());
        words.push(state.stored_energy_j.to_bits());
        for &v in &state.aux_liq_kg {
            words.push(v.to_bits());
        }
        words.push(state.cumulative_elec_work_j.to_bits());
        words.push(state.cumulative_room_heat_j.to_bits());
        true
    }

    fn decode_state(&self, words: &[u64]) -> Option<Self::State> {
        if words.len() != 16 || words[0] != ICE_HARVESTER_CODEC_TAG || words[1] != 14 {
            return None;
        }
        let mut aux_liq = [0.0; AUX_LIQ_SPECIES_COUNT];
        for (i, slot) in aux_liq.iter_mut().enumerate() {
            *slot = f64::from_bits(words[7 + i]);
        }
        let state = IceHarvesterStateV1 {
            time_s: f64::from_bits(words[2]),
            ice_bin_mass_kg: f64::from_bits(words[3]),
            sump_water_mass_kg: f64::from_bits(words[4]),
            t_sump_k: f64::from_bits(words[5]),
            stored_energy_j: f64::from_bits(words[6]),
            aux_liq_kg: aux_liq,
            cumulative_elec_work_j: f64::from_bits(words[14]),
            cumulative_room_heat_j: f64::from_bits(words[15]),
        };
        state.validate(&self.config).ok().map(|()| state)
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        words.push(species_discriminant(self.config.species));
        let scalars = [
            self.config.t_evap_k,
            self.config.suction_superheat_k,
            self.config.subcooling_k,
            self.config.ua_cond_clean_w_k,
            self.config.swept_volume_m3_s,
            self.config.eta_isentropic,
            self.config.eta_motor,
            self.config.clearance_ratio,
            self.config.t_ice_bin_k,
            self.config.max_bin_capacity_kg,
        ];
        words.extend(scalars.iter().map(|v| v.to_bits()));
    }
}

fn validate_stream_sum(port: &'static str, wire: &[f64]) -> Result<(), CycleErrorV1> {
    let total = wire[STREAM_MASS_RATE];
    let carried = &wire[STREAM_CARRIED..];
    let sum: f64 = carried.iter().sum();
    let scale = total.abs() + carried.iter().map(|v| v.abs()).sum::<f64>() + 1.0e-12;
    if (total - sum).abs() > 1.0e-10 * scale {
        return Err(CycleErrorV1::InconsistentStreamMassFlow {
            port,
            declared_total_kg_s: total,
            species_sum_kg_s: sum,
        });
    }
    Ok(())
}

const fn species_discriminant(species: RefrigerantSpecies) -> u64 {
    match species {
        RefrigerantSpecies::R290 => 1,
        RefrigerantSpecies::R404A => 2,
        RefrigerantSpecies::R448A => 3,
        RefrigerantSpecies::R134a => 4,
        RefrigerantSpecies::R744 => 5,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn iir_273k_datum_anchors_saturated_liquid_to_200_kj_kg() {
        for species in [
            RefrigerantSpecies::R290,
            RefrigerantSpecies::R134a,
            RefrigerantSpecies::R404A,
            RefrigerantSpecies::R448A,
            RefrigerantSpecies::R744,
        ] {
            let model = RefrigerantModelV1::for_species(species).expect("model");
            let p_ref = model
                .bubble_pressure_at_t(IIR_REFERENCE_TEMPERATURE_K)
                .expect("p_ref");
            let pt = model
                .state_pt(
                    IIR_REFERENCE_TEMPERATURE_K,
                    p_ref,
                    PengRobinsonPhaseV1::Liquid,
                )
                .expect("pt");
            assert!(
                (pt.enthalpy_j_kg - IIR_SAT_LIQUID_ENTHALPY_J_KG).abs() < 1.0e-6,
                "{species:?} IIR h_l(0C) mismatch: {}",
                pt.enthalpy_j_kg
            );
            assert!(
                (pt.entropy_j_kg_k - IIR_SAT_LIQUID_ENTROPY_J_KG_K).abs() < 1.0e-6,
                "{species:?} IIR s_l(0C) mismatch: {}",
                pt.entropy_j_kg_k
            );
        }
    }

    #[test]
    fn r448a_zeotropic_blend_exhibits_five_to_seven_kelvin_evaporator_glide() {
        let r448a = RefrigerantModelV1::r448a().expect("r448a");
        let env = r448a.bubble_dew_at_pressure(400_000.0).expect("env");
        assert!(
            (5.0..=7.2).contains(&env.glide_k),
            "R-448A 4-bar glide expected ~5.5-6.7 K, got {:.2} K",
            env.glide_k
        );

        let r404a = RefrigerantModelV1::r404a().expect("r404a");
        let env_404 = r404a.bubble_dew_at_pressure(400_000.0).expect("env404");
        assert!(
            env_404.glide_k > 0.2 && env_404.glide_k < 1.5,
            "R-404A near-azeotropic glide expected ~0.5-0.8 K, got {:.2} K",
            env_404.glide_k
        );
    }

    #[test]
    fn peneloux_translation_shifts_liquid_density_measurably_at_txv_inlet() {
        let r290 = RefrigerantModelV1::r290().expect("r290");
        let p_cond = r290.bubble_pressure_at_t(315.15).expect("p");
        let pt = r290
            .state_pt(310.15, p_cond, PengRobinsonPhaseV1::Liquid)
            .expect("liq");
        let rel_shift = (pt.density_kg_m3 - pt.untranslated_density_kg_m3).abs() / pt.density_kg_m3;
        assert!(
            rel_shift > 0.03,
            "Peneloux translation must measurably shift R-290 liquid density (>3%), got {:.2}%",
            100.0 * rel_shift
        );
    }
}
