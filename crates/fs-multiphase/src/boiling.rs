//! Full five-regime Nukiyama pool-boiling curve with safeguarded `C^1` transition
//! boiling, Schluender/Thome-Shakir binary mixture suppression, Lienhard-Dhir CHF,
//! Berenson/Bromley + Henry Leidenfrost film boiling, and explicit triglyceride
//! oil non-boiling refusal.
//!
//! All transcendentals route strictly through `fs_math::det`.

use core::fmt;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy};
use fs_math::det::{exp, expm1, ln, pow, sqrt};
use fs_ports::component::Proposal;
use fs_qty::{Dims, Pressure, Temperature};
use fs_thermochem::caloric::proximate::water_state_enthalpy_janaf;
use fs_thermochem::solution::VleSolverOptionsV1;
use fs_thermochem::solution::nrtl::NrtlModelV1;
use fs_thermochem::water::if97::regions::If97FluidState;
use fs_thermochem::water::if97::{If97Region, If97SaturationPair, If97State};
use fs_thermochem::water::transport::If97TransportState;

/// Standard terrestrial gravitational acceleration [m/s^2].
pub const STANDARD_GRAVITY_M_S2: f64 = 9.806_65;
/// Stefan-Boltzmann constant [W/(m^2 K^4)], exact post-2019 SI definition.
pub const STEFAN_BOLTZMANN_W_M2_K4: f64 = 5.670_374_419e-8;
/// Lienhard & Dhir (1973) horizontal flat-plate CHF constant (Incropera Eq. 10.6).
pub const LIENHARD_DHIR_FLAT_PLATE_C_CHF_V1: f64 = 0.149;
/// Zuber (1959) infinite horizontal flat-plate CHF constant.
pub const ZUBER_INFINITE_PLATE_C_CHF_V1: f64 = 0.131;
/// Berenson (1961) horizontal-plate minimum film-boiling heat flux constant.
pub const BERENSON_MIN_FLUX_CONSTANT_V1: f64 = 0.09;
/// Default Thome & Shakir (1987) / Schluender (1983) liquid mass-transfer coefficient [m/s].
pub const SCHLUNDER_BETA_L_M_S_V1: f64 = 3.0e-4;
/// Exact `fs-matdb` law id for a Rohsenow surface-fluid boiling card.
pub const ROHSENOW_SURFACE_LAW_ID_V1: &str = "rohsenow-boiling-surface";
/// Supported semantic version for Rohsenow surface-fluid cards.
pub const ROHSENOW_SURFACE_LAW_VERSION_V1: u32 = 1;

// Pure liquid ethanol standard formation enthalpy at 298.15 K, 1 bar [J/kg]
// (NIST WebBook: -277,600 J/mol / 0.04606844 kg/mol).
const ETHANOL_LIQUID_FORMATION_298K_J_KG: f64 = -6_025_817.240_610_057;
const WATER_MOLAR_MASS_KG_MOL: f64 = 0.018_015_28;
const ETHANOL_MOLAR_MASS_KG_MOL: f64 = 0.046_068_44;

/// Dedicated typed refusal for Nukiyama pool-boiling curve preparation and evaluation.
#[derive(Clone, Debug, PartialEq)]
pub enum BoilingErrorV1 {
    /// Round-2 Physics #3: pure triglyceride cooking oils undergo thermal cracking,
    /// smoking, and flash/autoignition without a thermodynamic boiling curve.
    NonBoilingLipidPhase {
        /// Oil identifier (e.g., `"canola_oil"`, `"peanut_oil"`).
        lipid_name: &'static str,
        /// Declared fresh smoke point [K].
        smoke_point_k: f64,
        /// Declared flash point [K].
        flash_point_k: f64,
    },
    /// Ambient/vessel barometric pressure is non-finite, non-positive, or unsupported.
    InvalidPressure {
        /// Rejected pressure [Pa].
        pressure_pa: f64,
    },
    /// Wall or reference temperature is non-finite or non-positive.
    InvalidTemperature {
        /// Rejected temperature [K].
        temperature_k: f64,
    },
    /// Surface card, fluid property, or geometry parameter is invalid.
    InvalidParameter(&'static str),
    /// Binary mixture mole fraction is outside `[0, 1]` or non-finite.
    InvalidComposition(&'static str),
    /// Boiling regime temperatures collapsed (e.g., near-critical state where `T_chf >= T_min`).
    InvalidRegimeOrdering(&'static str),
    /// Underlying `fs-thermochem` saturation, transport, or NRTL VLE query failed.
    Thermochem(&'static str),
    /// Intermediate or output quantity overflowed, underflowed, or became non-finite.
    Unrepresentable(&'static str),
}

impl fmt::Display for BoilingErrorV1 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonBoilingLipidPhase {
                lipid_name,
                smoke_point_k,
                flash_point_k,
            } => write!(
                f,
                "Nukiyama pool boiling refused on pure triglyceride lipid {lipid_name:?} \
                 (smoke={smoke_point_k:.1} K, flash={flash_point_k:.1} K); oils crack/smoke/flash \
                 without boiling"
            ),
            Self::InvalidPressure { pressure_pa } => {
                write!(
                    f,
                    "boiling pressure {pressure_pa} Pa is invalid or unsupported"
                )
            }
            Self::InvalidTemperature { temperature_k } => {
                write!(
                    f,
                    "boiling temperature {temperature_k} K must be finite and positive"
                )
            }
            Self::InvalidParameter(msg) => write!(f, "invalid boiling parameter: {msg}"),
            Self::InvalidComposition(msg) => {
                write!(f, "invalid boiling mixture composition: {msg}")
            }
            Self::InvalidRegimeOrdering(msg) => write!(f, "boiling regime ordering invalid: {msg}"),
            Self::Thermochem(msg) => {
                write!(f, "thermochemical saturation/VLE lookup failed: {msg}")
            }
            Self::Unrepresentable(msg) => write!(f, "unrepresentable boiling quantity: {msg}"),
        }
    }
}

impl core::error::Error for BoilingErrorV1 {}

/// Fluid specification routed into [`BoilingCurveV1::prepare`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum BoilingFluidKindV1 {
    /// Ordinary pure water/steam evaluated via IAPWS-IF97 + R12-08 + R15-11 + R1-76(2014).
    PureWater,
    /// Binary ethanol-water liquid mixture evaluated via `KM0.1` Aspen Plus NRTL VLE
    /// at the caller's barometric pressure (`x_etoh in [0, 1]`), supplying exact
    /// `Delta T_bp = T_dew - T_bub` for Schluender/Thome-Shakir suppression.
    EthanolWaterNrtl {
        /// Liquid ethanol mole fraction `x_EtOH` in `[0, 1]`.
        ethanol_mole_fraction: f64,
    },
    /// Caller-supplied explicit fluid saturation and transport properties.
    CustomFluid(BoilingFluidPropertiesV1),
    /// Pure triglyceride cooking oil (explicitly refused by Nukiyama boiling per
    /// Round-2 Physics #3; use cracking/smoking/flash kinetics in `MultiphaseVesselComponent`).
    TriglycerideOil {
        /// Descriptive lipid name.
        name: &'static str,
        /// Fresh smoke point [K].
        smoke_point_k: f64,
        /// Flash point [K].
        flash_point_k: f64,
    },
}

/// Complete coherent-SI fluid saturation, transport, mixture, and JANAF-datum properties
/// at one ambient barometric pressure `P`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BoilingFluidPropertiesV1 {
    /// Ambient / vessel pressure `P` [Pa].
    pub pressure_pa: f64,
    /// Bubble / saturation temperature `T_sat(P)` [K].
    pub saturation_temperature_k: f64,
    /// Binary mixture boiling range `Delta T_bp = T_dew - T_bub >= 0` [K] (`0` for pure fluids).
    pub boiling_range_k: f64,
    /// Saturated liquid density `rho_l` [kg/m^3].
    pub liquid_density_kg_m3: f64,
    /// Saturated vapor density `rho_v` [kg/m^3].
    pub vapor_density_kg_m3: f64,
    /// Saturated liquid dynamic viscosity `mu_l` [Pa s].
    pub liquid_viscosity_pa_s: f64,
    /// Vapor film dynamic viscosity `mu_v` [Pa s].
    pub vapor_film_viscosity_pa_s: f64,
    /// Vapor film density `rho_v,f` [kg/m^3] evaluated at representative film temperature.
    pub vapor_film_density_kg_m3: f64,
    /// Saturated liquid thermal conductivity `k_l` [W/(m K)].
    pub liquid_conductivity_w_m_k: f64,
    /// Vapor film thermal conductivity `k_v,f` [W/(m K)].
    pub vapor_film_conductivity_w_m_k: f64,
    /// Saturated liquid isobaric heat capacity `c_{p,l}` [J/(kg K)].
    pub liquid_cp_j_kg_k: f64,
    /// Vapor film isobaric heat capacity `c_{p,v}` [J/(kg K)].
    pub vapor_cp_j_kg_k: f64,
    /// Liquid volumetric thermal expansion coefficient `beta_l` [1/K].
    pub liquid_thermal_expansion_1_k: f64,
    /// Enthalpy of vaporization `h_fg` [J/kg].
    pub latent_heat_j_kg: f64,
    /// Liquid-vapor surface tension `sigma` [N/m].
    pub surface_tension_n_m: f64,
    /// Recommended Rohsenow Prandtl exponent `n` (`1.0` for pure water, up to `1.7` for organics).
    pub default_prandtl_exponent_n: f64,
    /// Dimensionless Rohsenow fluid-pair scale applied to water-calibrated `C_sf`
    /// (`1.0` for pure water; `exp(x_etoh * ln(0.0027 / 0.0130))` for ethanol-water per
    /// Incropera Table 10.1 ethanol-vs-water scaling).
    pub c_sf_fluid_scale: f64,
    /// Schluender/Thome-Shakir liquid mass-transfer coefficient `beta_L` [m/s].
    pub schlunder_beta_l_m_s: f64,
    /// Saturated liquid enthalpy on `janaf-298.15K-formation` datum [J/kg].
    pub liquid_enthalpy_janaf_j_kg: f64,
    /// Saturated vapor enthalpy on `janaf-298.15K-formation` datum [J/kg].
    pub vapor_enthalpy_janaf_j_kg: f64,
}

impl BoilingFluidPropertiesV1 {
    /// Validate positivity and finiteness of every physical property.
    pub fn validate(&self) -> Result<(), BoilingErrorV1> {
        if !self.pressure_pa.is_finite() || self.pressure_pa <= 0.0 {
            return Err(BoilingErrorV1::InvalidPressure {
                pressure_pa: self.pressure_pa,
            });
        }
        if !self.saturation_temperature_k.is_finite() || self.saturation_temperature_k <= 0.0 {
            return Err(BoilingErrorV1::InvalidTemperature {
                temperature_k: self.saturation_temperature_k,
            });
        }
        let positives = [
            self.liquid_density_kg_m3,
            self.vapor_density_kg_m3,
            self.liquid_viscosity_pa_s,
            self.vapor_film_viscosity_pa_s,
            self.vapor_film_density_kg_m3,
            self.liquid_conductivity_w_m_k,
            self.vapor_film_conductivity_w_m_k,
            self.liquid_cp_j_kg_k,
            self.vapor_cp_j_kg_k,
            self.liquid_thermal_expansion_1_k,
            self.latent_heat_j_kg,
            self.surface_tension_n_m,
            self.default_prandtl_exponent_n,
            self.c_sf_fluid_scale,
            self.schlunder_beta_l_m_s,
        ];
        if !positives.into_iter().all(|v| v.is_finite() && v > 0.0) {
            return Err(BoilingErrorV1::InvalidParameter(
                "all fluid boiling properties must be finite and strictly positive",
            ));
        }
        if !self.boiling_range_k.is_finite() || self.boiling_range_k < 0.0 {
            return Err(BoilingErrorV1::InvalidParameter(
                "boiling_range_k must be finite and nonnegative",
            ));
        }
        if self.liquid_density_kg_m3 <= self.vapor_density_kg_m3
            || self.liquid_density_kg_m3 <= self.vapor_film_density_kg_m3
        {
            return Err(BoilingErrorV1::InvalidParameter(
                "liquid density must strictly exceed vapor density",
            ));
        }
        if !self.liquid_enthalpy_janaf_j_kg.is_finite()
            || !self.vapor_enthalpy_janaf_j_kg.is_finite()
        {
            return Err(BoilingErrorV1::InvalidParameter(
                "JANAF liquid and vapor enthalpies must be finite",
            ));
        }
        Ok(())
    }

    /// Evaluate exact IAPWS-IF97 + R12-08 + R15-11 + R1-76(2014) pure-water properties
    /// at caller-specified barometric pressure `pressure`.
    pub fn pure_water_at_pressure(pressure: Pressure) -> Result<Self, BoilingErrorV1> {
        let p = pressure.value();
        if !p.is_finite() || p <= 0.0 {
            return Err(BoilingErrorV1::InvalidPressure { pressure_pa: p });
        }
        let pair = If97SaturationPair::at_pressure(pressure)
            .map_err(|_| BoilingErrorV1::InvalidPressure { pressure_pa: p })?;
        let t_sat = pair.temperature().value();
        let liq = *pair.liquid();
        let vap = *pair.vapor();

        let liq_tr = If97TransportState::from_state(If97FluidState::Gibbs(liq)).map_err(|_| {
            BoilingErrorV1::Thermochem("IAPWS saturated liquid transport evaluation failed")
        })?;
        let surf = pair
            .surface_tension()
            .map_err(|_| BoilingErrorV1::Thermochem("IAPWS surface tension evaluation failed"))?;

        // Representative film vapor state at T_sat + 57 K (midpoint of saturated water
        // and ~214 C Leidenfrost wall in superheated Region 2 steam at same P) for
        // Berenson/Bromley vapor-film transport properties.
        let t_film_est = (t_sat + 57.0).min(1000.0);
        let vap_film = If97State::at_pressure_temperature(
            pressure,
            Temperature::new(t_film_est),
            If97Region::Region2,
        )
        .unwrap_or(vap);
        let vap_film_tr =
            If97TransportState::from_state(If97FluidState::Gibbs(vap_film)).map_err(|_| {
                BoilingErrorV1::Thermochem("IAPWS vapor-film transport evaluation failed")
            })?;

        let h_l_janaf = water_state_enthalpy_janaf(&liq)
            .map_err(|_| BoilingErrorV1::Thermochem("JANAF water liquid enthalpy failed"))?
            .value();
        let h_v_janaf = water_state_enthalpy_janaf(&vap)
            .map_err(|_| BoilingErrorV1::Thermochem("JANAF water vapor enthalpy failed"))?
            .value();

        // Near 277 K water thermal expansion crosses zero; guard natural-convection beta_l >= 1e-5.
        let beta_l = liq.thermal_expansion().value().abs().max(1e-5);

        let props = Self {
            pressure_pa: p,
            saturation_temperature_k: t_sat,
            boiling_range_k: 0.0,
            liquid_density_kg_m3: liq.density().value(),
            vapor_density_kg_m3: vap.density().value(),
            liquid_viscosity_pa_s: liq_tr.dynamic_viscosity().value(),
            vapor_film_viscosity_pa_s: vap_film_tr.dynamic_viscosity().value(),
            vapor_film_density_kg_m3: vap_film.density().value(),
            liquid_conductivity_w_m_k: liq_tr.thermal_conductivity().value(),
            vapor_film_conductivity_w_m_k: vap_film_tr.thermal_conductivity().value(),
            liquid_cp_j_kg_k: liq.cp().value(),
            vapor_cp_j_kg_k: vap_film.cp().value(),
            liquid_thermal_expansion_1_k: beta_l,
            latent_heat_j_kg: pair.latent_enthalpy().value(),
            surface_tension_n_m: surf.surface_tension().value(),
            default_prandtl_exponent_n: 1.0,
            c_sf_fluid_scale: 1.0,
            schlunder_beta_l_m_s: SCHLUNDER_BETA_L_M_S_V1,
            liquid_enthalpy_janaf_j_kg: h_l_janaf,
            vapor_enthalpy_janaf_j_kg: h_v_janaf,
        };
        props.validate()?;
        Ok(props)
    }

    /// Evaluate ethanol-water binary mixture boiling properties at caller-specified
    /// pressure `pressure` and liquid ethanol mole fraction `ethanol_mole_fraction`
    /// using `KM0.1`'s Aspen Plus NRTL Modified Raoult VLE model.
    #[allow(clippy::too_many_lines)]
    pub fn ethanol_water_nrtl_at_pressure(
        pressure: Pressure,
        ethanol_mole_fraction: f64,
    ) -> Result<Self, BoilingErrorV1> {
        let p = pressure.value();
        if !p.is_finite() || p <= 0.0 {
            return Err(BoilingErrorV1::InvalidPressure { pressure_pa: p });
        }
        if !ethanol_mole_fraction.is_finite() || !(0.0..=1.0).contains(&ethanol_mole_fraction) {
            return Err(BoilingErrorV1::InvalidComposition(
                "ethanol mole fraction must be finite and in [0, 1]",
            ));
        }
        if ethanol_mole_fraction == 0.0 {
            return Self::pure_water_at_pressure(pressure);
        }

        let vle = NrtlModelV1::ethanol_water_vle_v1();
        let (bubble, _dew, dt_bp_raw) = vle
            .boiling_range_k(
                pressure,
                &[ethanol_mole_fraction, 1.0 - ethanol_mole_fraction],
                VleSolverOptionsV1::default(),
            )
            .map_err(|_| {
                BoilingErrorV1::Thermochem("NRTL ethanol-water boiling-range VLE solve failed")
            })?;
        let t_bub = bubble.temperature.value();
        let dt_bp = dt_bp_raw.max(0.0);
        let y_etoh = bubble.vapor_y[0].clamp(0.0, 1.0);

        // Liquid and vapor mass fractions.
        let m_liq = ethanol_mole_fraction * ETHANOL_MOLAR_MASS_KG_MOL
            + (1.0 - ethanol_mole_fraction) * WATER_MOLAR_MASS_KG_MOL;
        let w_etoh_l = (ethanol_mole_fraction * ETHANOL_MOLAR_MASS_KG_MOL) / m_liq;
        let w_h2o_l = 1.0 - w_etoh_l;

        let m_vap = y_etoh * ETHANOL_MOLAR_MASS_KG_MOL + (1.0 - y_etoh) * WATER_MOLAR_MASS_KG_MOL;
        let w_etoh_v = (y_etoh * ETHANOL_MOLAR_MASS_KG_MOL) / m_vap;
        let w_h2o_v = 1.0 - w_etoh_v;

        // Pure water reference at the same pressure (for IAPWS baseline).
        let water = Self::pure_water_at_pressure(pressure)?;

        // Published saturated ethanol thermophysical properties around 330-375 K (Poling / DIPPR):
        // rho_l(T) ≈ 789 - 0.87*(T - 293.15), mu_l ≈ 4.35e-4 Pa*s at 351 K, k_l ≈ 0.155 W/(m K),
        // cp_l ≈ 3050 J/(kg K), h_fg ≈ 840 kJ/kg, sigma_etoh ≈ 0.0177 N/m at 351 K.
        let t_shift = t_bub - 351.45;
        let rho_l_e = (736.0 - 0.92 * t_shift).clamp(620.0, 820.0);
        let rho_v_e = (p * ETHANOL_MOLAR_MASS_KG_MOL)
            / (fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * t_bub);
        let rho_vf_e = (p * ETHANOL_MOLAR_MASS_KG_MOL)
            / (fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * (t_bub + 57.0));
        let mu_l_e = (4.35e-4 * exp(-0.016 * t_shift)).clamp(1.5e-4, 1.8e-3);
        let mu_vf_e = (1.16e-5 + 2.8e-8 * (t_bub + 57.0 - 351.45)).clamp(8.0e-6, 2.5e-5);
        let k_l_e = (0.155 - 2.0e-4 * t_shift).clamp(0.11, 0.19);
        let k_vf_e = (0.0235 + 1.1e-4 * (t_bub + 57.0 - 351.45)).clamp(0.014, 0.055);
        let cp_l_e = (3050.0 + 7.5 * t_shift).clamp(2300.0, 3800.0);
        let cp_v_e = (1780.0 + 3.2 * (t_bub + 57.0 - 351.45)).clamp(1400.0, 2600.0);
        let h_fg_e = (841_000.0 - 2_350.0 * t_shift).clamp(550_000.0, 960_000.0);
        let sigma_e = (0.0177 - 9.0e-5 * t_shift).clamp(0.008, 0.026);
        let beta_l_e = 1.15e-3;

        // Mass/volume-consistent mixture rules (plus Connors-Wright surface enrichment
        // for strong ethanol surface-tension depression in aqueous solutions).
        let rho_l = 1.0 / (w_h2o_l / water.liquid_density_kg_m3 + w_etoh_l / rho_l_e);
        let rho_v = (p * m_vap) / (fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * t_bub);
        let rho_vf =
            (p * m_vap) / (fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * (t_bub + 57.0));
        let mu_l = exp(w_h2o_l * ln(water.liquid_viscosity_pa_s) + w_etoh_l * ln(mu_l_e));
        let mu_vf = w_h2o_v * water.vapor_film_viscosity_pa_s + w_etoh_v * mu_vf_e;
        let k_l = w_h2o_l * water.liquid_conductivity_w_m_k + w_etoh_l * k_l_e;
        let k_vf = w_h2o_v * water.vapor_film_conductivity_w_m_k + w_etoh_v * k_vf_e;
        let cp_l = w_h2o_l * water.liquid_cp_j_kg_k + w_etoh_l * cp_l_e;
        let cp_v = w_h2o_v * water.vapor_cp_j_kg_k + w_etoh_v * cp_v_e;
        let beta_l = w_h2o_l * water.liquid_thermal_expansion_1_k + w_etoh_l * beta_l_e;
        let h_fg = w_h2o_v * water.latent_heat_j_kg + w_etoh_v * h_fg_e;

        // Connors-Wright (1989) surface-enrichment form for aqueous ethanol surface tension:
        // sigma_mix = sigma_w - [1 + b*(1 - x)/(1 - a*(1 - x))] * x * (sigma_w - sigma_e)
        // with published Connors-Wright constants a = 0.85, b = 0.85 for ethanol-water.
        let one_minus_x = 1.0 - ethanol_mole_fraction;
        let enrichment = 1.0 + (0.85 * one_minus_x) / (1.0 - 0.85 * one_minus_x);
        let sigma = (water.surface_tension_n_m
            - enrichment * ethanol_mole_fraction * (water.surface_tension_n_m - sigma_e))
            .clamp(sigma_e, water.surface_tension_n_m);

        // JANAF-298.15K-formation enthalpies for liquid and vapor phases.
        let h_l_e_janaf = ETHANOL_LIQUID_FORMATION_298K_J_KG + 2_600.0 * (t_bub - 298.15);
        let h_v_e_janaf = h_l_e_janaf + h_fg_e;
        let h_l_janaf = w_h2o_l * water.liquid_enthalpy_janaf_j_kg + w_etoh_l * h_l_e_janaf;
        let h_v_janaf = w_h2o_v * water.vapor_enthalpy_janaf_j_kg + w_etoh_v * h_v_e_janaf;

        // Avoid unused helper bindings while keeping rho_v_e / rho_vf_e exact for pure ethanol.
        let (rho_v_final, rho_vf_final) = if (1.0 - y_etoh) < 1e-12 {
            (rho_v_e, rho_vf_e)
        } else {
            (rho_v, rho_vf)
        };

        let props = Self {
            pressure_pa: p,
            saturation_temperature_k: t_bub,
            boiling_range_k: dt_bp,
            liquid_density_kg_m3: rho_l,
            vapor_density_kg_m3: rho_v_final,
            liquid_viscosity_pa_s: mu_l,
            vapor_film_viscosity_pa_s: mu_vf,
            vapor_film_density_kg_m3: rho_vf_final,
            liquid_conductivity_w_m_k: k_l,
            vapor_film_conductivity_w_m_k: k_vf,
            liquid_cp_j_kg_k: cp_l,
            vapor_cp_j_kg_k: cp_v,
            liquid_thermal_expansion_1_k: beta_l,
            latent_heat_j_kg: h_fg,
            surface_tension_n_m: sigma,
            default_prandtl_exponent_n: 1.0 + 0.7 * ethanol_mole_fraction,
            // Incropera Table 10.1 Rohsenow ethanol-vs-water C_sf scaling (0.0027 / 0.0130).
            c_sf_fluid_scale: exp(ethanol_mole_fraction * ln(0.0027 / 0.0130)),
            schlunder_beta_l_m_s: SCHLUNDER_BETA_L_M_S_V1,
            liquid_enthalpy_janaf_j_kg: h_l_janaf,
            vapor_enthalpy_janaf_j_kg: h_v_janaf,
        };
        props.validate()?;
        Ok(props)
    }
}

/// Leidenfrost minimum film-boiling temperature mode.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum LeidenfrostModeV1 {
    /// Saturated horizontal pool minimum (Berenson 1961 + Henry 1974 wall thermal
    /// effusivity; ~195-220 C for 1 atm water on cookware metals).
    SaturatedPool,
    /// Subcooled/impinging droplet or pan-quench dynamic Leidenfrost (Henry 1974
    /// subcooled correction + impact multiplier; ~250-310 C for cool water on hot steel).
    SubcooledImpinging {
        /// Liquid subcooling `T_sat - T_liquid >= 0` [K].
        subcooling_k: f64,
        /// Dimensionless dynamic impingement multiplier (`>= 1.0`, typically `1.0..=1.25`).
        dynamic_multiplier: f64,
    },
}

/// Heated wall surface-fluid boiling parameters (Rohsenow + Lienhard-Dhir + Henry + Bromley).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BoilingSurfaceCardV1 {
    /// Rohsenow surface-fluid coefficient `C_sf > 0` (e.g., `0.013` SS304, `0.0068` scored Cu).
    pub c_sf: f64,
    /// Optional explicit Prandtl exponent override `n > 0` (`None` uses fluid default).
    pub prandtl_exponent_override: Option<f64>,
    /// CHF geometry constant (`0.149` Lienhard-Dhir flat plate, `0.131` Zuber infinite plate).
    pub chf_geometry_constant: f64,
    /// Static wetting contact angle [rad] in `(0, pi)`.
    pub contact_angle_rad: f64,
    /// Thermal emissivity of heated surface in `[0, 1]` (bare SS304 `0.20`, seasoned cast iron `0.88`).
    pub surface_emissivity: f64,
    /// Characteristic heater dimension `L_c > 0` [m].
    pub characteristic_length_m: f64,
    /// Wall solid thermal conductivity `k_w > 0` [W/(m K)] (Henry 1974 effusivity).
    pub wall_conductivity_w_m_k: f64,
    /// Wall solid density `rho_w > 0` [kg/m^3].
    pub wall_density_kg_m3: f64,
    /// Wall solid specific heat `c_{p,w} > 0` [J/(kg K)].
    pub wall_cp_j_kg_k: f64,
    /// Leidenfrost minimum mode (saturated pool vs subcooled/impinging skillet).
    pub leidenfrost_mode: LeidenfrostModeV1,
}

impl BoilingSurfaceCardV1 {
    /// Validate all surface-fluid boiling card parameters fail-closed.
    pub fn validate(&self) -> Result<(), BoilingErrorV1> {
        let positives = [
            self.c_sf,
            self.chf_geometry_constant,
            self.characteristic_length_m,
            self.wall_conductivity_w_m_k,
            self.wall_density_kg_m3,
            self.wall_cp_j_kg_k,
        ];
        if !positives.into_iter().all(|v| v.is_finite() && v > 0.0) {
            return Err(BoilingErrorV1::InvalidParameter(
                "C_sf, C_CHF, characteristic length, and wall thermal properties must be finite positive",
            ));
        }
        if let Some(n) = self.prandtl_exponent_override
            && (!n.is_finite() || n <= 0.0)
        {
            return Err(BoilingErrorV1::InvalidParameter(
                "Rohsenow Prandtl exponent override must be finite and positive",
            ));
        }
        if !self.contact_angle_rad.is_finite()
            || self.contact_angle_rad <= 0.0
            || self.contact_angle_rad >= core::f64::consts::PI
        {
            return Err(BoilingErrorV1::InvalidParameter(
                "contact angle must lie strictly inside (0, pi) rad",
            ));
        }
        if !self.surface_emissivity.is_finite() || !(0.0..=1.0).contains(&self.surface_emissivity) {
            return Err(BoilingErrorV1::InvalidParameter(
                "surface emissivity must be finite and in [0, 1]",
            ));
        }
        if let LeidenfrostModeV1::SubcooledImpinging {
            subcooling_k,
            dynamic_multiplier,
        } = self.leidenfrost_mode
            && (!subcooling_k.is_finite()
                || subcooling_k < 0.0
                || !dynamic_multiplier.is_finite()
                || dynamic_multiplier < 1.0)
        {
            return Err(BoilingErrorV1::InvalidParameter(
                "subcooled impinging Leidenfrost requires subcooling >= 0 and multiplier >= 1",
            ));
        }
        Ok(())
    }

    /// Load a validated Rohsenow boiling surface card from an immutable `fs-matdb` card.
    ///
    /// Expected parameters: `c_sf` (dimless), `c_chf` (dimless), `contact_angle_rad` (dimless),
    /// `emissivity` (dimless), `characteristic_length_m` (`[1,0,0,0,0,0]`),
    /// `wall_k` (`[1,1,-3,-1,0,0]`), `wall_rho` (`[-3,1,0,0,0,0]`), `wall_cp` (`[2,0,-2,-1,0,0]`),
    /// and optional `prandtl_n` (dimless).
    pub fn from_card(card: &ConstitutiveModelCard) -> Result<Self, BoilingErrorV1> {
        card.validate().map_err(|_| {
            BoilingErrorV1::InvalidParameter("boiling matdb card failed validation")
        })?;
        if card.law.0 != ROHSENOW_SURFACE_LAW_ID_V1 {
            return Err(BoilingErrorV1::InvalidParameter(
                "boiling matdb card law id mismatch",
            ));
        }
        if card.law_version != ROHSENOW_SURFACE_LAW_VERSION_V1 {
            return Err(BoilingErrorV1::InvalidParameter(
                "boiling matdb card law version mismatch",
            ));
        }
        if card.state_schema_version != 0
            || card.initial_state != InitialStatePolicy::ZeroInternalState
        {
            return Err(BoilingErrorV1::InvalidParameter(
                "boiling matdb card must be stateless",
            ));
        }
        let get = |name: &'static str, expected: Dims| -> Result<f64, BoilingErrorV1> {
            let p = card
                .parameters
                .get(name)
                .ok_or(BoilingErrorV1::InvalidParameter(
                    "boiling matdb card missing required parameter",
                ))?;
            if p.dims != expected {
                return Err(BoilingErrorV1::InvalidParameter(
                    "boiling matdb card parameter dimension mismatch",
                ));
            }
            Ok(p.value)
        };
        let prandtl_exponent_override = match card.parameters.get("prandtl_n") {
            Some(p) => {
                if p.dims != Dims::NONE {
                    return Err(BoilingErrorV1::InvalidParameter(
                        "boiling matdb card prandtl_n must be dimensionless",
                    ));
                }
                Some(p.value)
            }
            None => None,
        };
        let surface = Self {
            c_sf: get("c_sf", Dims::NONE)?,
            prandtl_exponent_override,
            chf_geometry_constant: get("c_chf", Dims::NONE)?,
            contact_angle_rad: get("contact_angle_rad", Dims::NONE)?,
            surface_emissivity: get("emissivity", Dims::NONE)?,
            characteristic_length_m: get("characteristic_length_m", Dims([1, 0, 0, 0, 0, 0]))?,
            wall_conductivity_w_m_k: get("wall_k", Dims([1, 1, -3, -1, 0, 0]))?,
            wall_density_kg_m3: get("wall_rho", Dims([-3, 1, 0, 0, 0, 0]))?,
            wall_cp_j_kg_k: get("wall_cp", Dims([2, 0, -2, -1, 0, 0]))?,
            leidenfrost_mode: LeidenfrostModeV1::SaturatedPool,
        };
        surface.validate()?;
        Ok(surface)
    }

    /// Mechanically polished stainless steel SS304 cookware preset (Incropera Table 10.1:
    /// `C_sf = 0.0130`, Lienhard-Dhir `C_CHF = 0.149`, bare SS304 emissivity `0.20`).
    #[must_use]
    pub const fn stainless_steel_304_v1() -> Self {
        Self {
            c_sf: 0.0130,
            prandtl_exponent_override: None,
            chf_geometry_constant: LIENHARD_DHIR_FLAT_PLATE_C_CHF_V1,
            contact_angle_rad: 1.15,
            surface_emissivity: 0.20,
            characteristic_length_m: 0.25,
            wall_conductivity_w_m_k: 16.2,
            wall_density_kg_m3: 7_900.0,
            wall_cp_j_kg_k: 500.0,
            leidenfrost_mode: LeidenfrostModeV1::SaturatedPool,
        }
    }

    /// Scored copper boiling plate preset (Incropera Table 10.1: `C_sf = 0.0068`,
    /// Lienhard-Dhir `C_CHF = 0.149`, copper `k_w = 385 W/(m K)`).
    #[must_use]
    pub const fn scored_copper_v1() -> Self {
        Self {
            c_sf: 0.0068,
            prandtl_exponent_override: None,
            chf_geometry_constant: LIENHARD_DHIR_FLAT_PLATE_C_CHF_V1,
            contact_angle_rad: 0.95,
            surface_emissivity: 0.15,
            characteristic_length_m: 0.25,
            wall_conductivity_w_m_k: 385.0,
            wall_density_kg_m3: 8_960.0,
            wall_cp_j_kg_k: 385.0,
            leidenfrost_mode: LeidenfrostModeV1::SaturatedPool,
        }
    }

    /// Seasoned cast-iron skillet preset (`C_sf = 0.0132`, seasoned IR emissivity `0.88`,
    /// cast iron `k_w = 55 W/(m K)`).
    #[must_use]
    pub const fn seasoned_cast_iron_v1() -> Self {
        Self {
            c_sf: 0.0132,
            prandtl_exponent_override: None,
            chf_geometry_constant: LIENHARD_DHIR_FLAT_PLATE_C_CHF_V1,
            contact_angle_rad: 1.35,
            surface_emissivity: 0.88,
            characteristic_length_m: 0.25,
            wall_conductivity_w_m_k: 55.0,
            wall_density_kg_m3: 7_200.0,
            wall_cp_j_kg_k: 460.0,
            leidenfrost_mode: LeidenfrostModeV1::SaturatedPool,
        }
    }

    /// Return a copy of this surface card with a custom [`LeidenfrostModeV1`].
    #[must_use]
    pub const fn with_leidenfrost_mode(mut self, mode: LeidenfrostModeV1) -> Self {
        self.leidenfrost_mode = mode;
        self
    }

    /// Return a copy of this surface card with a custom `C_CHF` geometry constant.
    #[must_use]
    pub const fn with_chf_constant(mut self, chf_geometry_constant: f64) -> Self {
        self.chf_geometry_constant = chf_geometry_constant;
        self
    }
}

/// Active regime along the five-region Nukiyama pool-boiling curve.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BoilingRegimeV1 {
    /// Wall temperature `T_w <= T_sat` (zero boiling flux; owned by single-phase convection).
    BelowSaturation,
    /// Single-phase turbulent natural convection (`0 < Delta T_e <= Delta T_onb`).
    NaturalConvection,
    /// Nucleate boiling (`Rohsenow` x `Schluender/Thome-Shakir` mixture suppression + `C^1` CHF cap).
    NucleateBoiling,
    /// Safeguarded `C^1` monotone transition boiling (`dq''/dT_w <= 0`, `T_chf < T_w < T_min`).
    TransitionBoiling,
    /// Bromley / Berenson convective + radiative film boiling (`T_w >= T_min`).
    FilmBoiling,
}

/// Key regime-transition wall temperatures [K] and anchor heat fluxes [W/m^2].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BoilingRegimeBoundariesV1 {
    /// Saturation / bubble-point temperature `T_sat` [K].
    pub t_sat_k: f64,
    /// Onset of nucleate boiling crossover temperature `T_onb` [K].
    pub t_onb_k: f64,
    /// Critical heat flux peak temperature `T_chf` [K].
    pub t_chf_k: f64,
    /// Leidenfrost minimum film-boiling temperature `T_min` [K].
    pub t_leidenfrost_k: f64,
    /// Lienhard-Dhir / Zuber critical heat flux `q''_CHF` [W/m^2].
    pub q_chf_w_m2: f64,
    /// Leidenfrost minimum heat flux `q''_min` [W/m^2].
    pub q_leidenfrost_w_m2: f64,
}

/// Point evaluation of [`BoilingCurveV1`] at a wall temperature `T_w`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BoilingEvaluationV1 {
    /// Evaluated wall temperature `T_w` [K].
    pub wall_temperature: Temperature,
    /// Saturation / bubble-point temperature `T_sat(P)` [K].
    pub saturation_temperature: Temperature,
    /// Wall superheat `Delta T_e = T_w - T_sat` [K].
    pub superheat_k: f64,
    /// Active Nukiyama regime.
    pub regime: BoilingRegimeV1,
    /// Wall-to-fluid heat flux `q''(T_w)` [W/m^2].
    pub heat_flux_w_m2: f64,
    /// Effective heat-transfer coefficient `h = q'' / max(Delta T_e, 0)` [W/(m^2 K)].
    pub htc_w_m2_k: f64,
    /// Exact piecewise-analytic tangent `dq''/dT_w` [W/(m^2 K)].
    pub dq_dtw_w_m2_k: f64,
    /// Ideal pure-fluid Rohsenow flux at `Delta T_e` [W/m^2] (diagnostic baseline).
    pub ideal_rohsenow_flux_w_m2: f64,
    /// Schluender/Thome-Shakir binary suppression factor `q''_mix / q''_ideal` in `(0, 1]`.
    pub mixture_suppression_factor: f64,
    /// Evaporation mass flux `q'' / h_fg` [kg/(m^2 s)] (`0` below ONB/saturation).
    pub vaporization_mass_flux_kg_m2_s: f64,
    /// Saturated vapor enthalpy on `janaf-298.15K-formation` datum [J/kg].
    pub vapor_enthalpy_janaf_j_kg: f64,
    /// Saturated liquid enthalpy on `janaf-298.15K-formation` datum [J/kg].
    pub liquid_enthalpy_janaf_j_kg: f64,
}

/// Time-step capping recommendation for crossing ONB, CHF, or Leidenfrost boundaries
/// cleanly inside a coupled implicit Newton integrator (`KM6.1a` / `KM2.2`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BoilingStepAdviceV1 {
    /// Whether the trial step crossed a regime boundary or exceeded `max_delta_t_k`.
    pub shorten_recommended: bool,
    /// Recommended maximum step duration [s] (`<= trial_dt_s`).
    pub recommended_max_dt_s: f64,
    /// Optional `fs_ports::component::Proposal::Shorten` payload when shortening is advised.
    pub shorten_proposal: Option<Proposal>,
}

/// Prepared, `C^1`-smooth five-regime Nukiyama pool-boiling curve for one
/// `(fluid, surface, pressure)` operating point.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BoilingCurveV1 {
    props: BoilingFluidPropertiesV1,
    surface: BoilingSurfaceCardV1,
    c_nc: f64,
    k_rohsenow: f64,
    schlunder_alpha: f64,
    film_conv_coeff: f64,
    dt_onb: f64,
    dt_onb_lo: f64,
    dt_onb_hi: f64,
    q_onb_lo: f64,
    dq_onb_lo: f64,
    q_onb_hi: f64,
    dq_onb_hi: f64,
    dt_chf_lo: f64,
    q_chf_lo: f64,
    dq_chf_lo: f64,
    dt_chf: f64,
    q_chf: f64,
    dt_min: f64,
    q_valley: f64,
    dt_film_hi: f64,
    q_film_hi: f64,
    dq_film_hi: f64,
}

impl BoilingCurveV1 {
    /// Prepare a validated, `C^1`-smooth Nukiyama boiling curve at explicit pressure `pressure`.
    ///
    /// # Errors
    /// Refuses pure triglyceride oils with [`BoilingErrorV1::NonBoilingLipidPhase`],
    /// invalid pressures/parameters, or collapsed near-critical regimes.
    #[allow(clippy::too_many_lines)]
    pub fn prepare(
        fluid: BoilingFluidKindV1,
        surface: BoilingSurfaceCardV1,
        pressure: Pressure,
    ) -> Result<Self, BoilingErrorV1> {
        surface.validate()?;
        let props = match fluid {
            BoilingFluidKindV1::TriglycerideOil {
                name,
                smoke_point_k,
                flash_point_k,
            } => {
                return Err(BoilingErrorV1::NonBoilingLipidPhase {
                    lipid_name: name,
                    smoke_point_k,
                    flash_point_k,
                });
            }
            BoilingFluidKindV1::PureWater => {
                BoilingFluidPropertiesV1::pure_water_at_pressure(pressure)?
            }
            BoilingFluidKindV1::EthanolWaterNrtl {
                ethanol_mole_fraction,
            } => BoilingFluidPropertiesV1::ethanol_water_nrtl_at_pressure(
                pressure,
                ethanol_mole_fraction,
            )?,
            BoilingFluidKindV1::CustomFluid(mut custom) => {
                let p = pressure.value();
                if !p.is_finite() || p <= 0.0 {
                    return Err(BoilingErrorV1::InvalidPressure { pressure_pa: p });
                }
                custom.pressure_pa = p;
                custom.validate()?;
                custom
            }
        };

        let g = STANDARD_GRAVITY_M_S2;
        let rho_l = props.liquid_density_kg_m3;
        let rho_v = props.vapor_density_kg_m3;
        let rho_vf = props.vapor_film_density_kg_m3;
        let drho = rho_l - rho_v;
        let drho_f = rho_l - rho_vf;
        let mu_l = props.liquid_viscosity_pa_s;
        let mu_vf = props.vapor_film_viscosity_pa_s;
        let k_l = props.liquid_conductivity_w_m_k;
        let k_vf = props.vapor_film_conductivity_w_m_k;
        let cp_l = props.liquid_cp_j_kg_k;
        let h_fg = props.latent_heat_j_kg;
        let sigma = props.surface_tension_n_m;
        let pr_l = (mu_l * cp_l) / k_l;
        let n_pr = surface
            .prandtl_exponent_override
            .unwrap_or(props.default_prandtl_exponent_n);

        // 1. Turbulent natural convection coefficient C_nc: q''_nc = C_nc * dT^(4/3).
        let nc_group =
            (g * props.liquid_thermal_expansion_1_k * rho_l * rho_l * cp_l) / (mu_l * k_l);
        let c_nc = 0.15 * k_l * pow(nc_group, 1.0 / 3.0);

        // 2. Rohsenow nucleate boiling coefficient K_Roh: q''_id = K_Roh * dT^3.
        let pr_pow = pow(pr_l, n_pr);
        let c_sf_eff = if surface.prandtl_exponent_override.is_some() {
            surface.c_sf
        } else {
            surface.c_sf * props.c_sf_fluid_scale
        };
        let roh_ratio = cp_l / (c_sf_eff * h_fg * pr_pow);
        let k_rohsenow =
            mu_l * h_fg * sqrt((g * drho) / sigma) * (roh_ratio * roh_ratio * roh_ratio);
        let schlunder_alpha = k_rohsenow / (rho_l * h_fg * props.schlunder_beta_l_m_s);

        // 3. Lienhard-Dhir / Zuber Critical Heat Flux q''_CHF.
        let q_chf = surface.chf_geometry_constant
            * h_fg
            * rho_v
            * pow((sigma * g * drho) / (rho_v * rho_v), 0.25);

        // 4. Exact nucleate boiling crossover dT_chf_raw where q''_nuc(dT) == q''_CHF.
        let dt_chf_raw = solve_nucleate_superheat_for_flux(
            q_chf,
            k_rohsenow,
            props.boiling_range_k,
            schlunder_alpha,
        )?;

        // 5. Exact ONB crossover where q''_nc(dT) == q''_nuc(dT).
        let dt_onb = solve_onb_superheat(
            c_nc,
            k_rohsenow,
            props.boiling_range_k,
            schlunder_alpha,
            dt_chf_raw,
        )?;

        // 6. Berenson (1961) + Henry (1974) Leidenfrost minimum dT_min and q''_min.
        let lambda_c = sqrt(sigma / (g * drho_f));
        let q_min_ber = BERENSON_MIN_FLUX_CONSTANT_V1
            * rho_v
            * h_fg
            * pow(
                (sigma * g * drho) / ((rho_l + rho_v) * (rho_l + rho_v)),
                0.25,
            );

        let dt_min_ber = 0.127
            * ((rho_vf * h_fg) / k_vf)
            * pow((g * drho_f) / (rho_l + rho_vf), 2.0 / 3.0)
            * lambda_c
            * pow(mu_vf / (g * drho_f), 1.0 / 3.0);

        // Henry (1974) Eq. (8) wall thermal-effusivity correction.
        let liq_eff_sq = k_l * rho_l * cp_l;
        let wall_eff_sq =
            surface.wall_conductivity_w_m_k * surface.wall_density_kg_m3 * surface.wall_cp_j_kg_k;
        let henry_factor = 0.42
            * pow(
                (liq_eff_sq / wall_eff_sq) * (h_fg / (surface.wall_cp_j_kg_k * dt_min_ber)),
                0.6,
            );

        // Saturated pool minimum (Berenson + Henry; falls naturally in ~195-220 C for 1 atm water
        // on SS304/cast iron; higher-conductivity Cu uses Berenson saturated-vapor scale floor).
        let dt_min_ber_sat = 0.127
            * ((rho_v * h_fg) / k_vf)
            * pow((g * drho) / (rho_l + rho_v), 2.0 / 3.0)
            * sqrt(sigma / (g * drho))
            * pow(mu_vf / (g * drho), 1.0 / 3.0);
        let dt_min_sat =
            (dt_min_ber * (1.0 + henry_factor)).max(dt_min_ber_sat * (1.0 + 0.5 * henry_factor));

        let dt_min = match surface.leidenfrost_mode {
            LeidenfrostModeV1::SaturatedPool => dt_min_sat,
            LeidenfrostModeV1::SubcooledImpinging {
                subcooling_k,
                dynamic_multiplier,
            } => (dt_min_ber + subcooling_k) * (1.0 + henry_factor) * dynamic_multiplier,
        };

        // 7. Bromley / Berenson horizontal-plate film-boiling convective prefactor:
        //    q''_conv_film(dT) = film_conv_coeff * dT^(3/4) * (h_fg + 0.5*cp_v*dT)^(1/4)
        let film_conv_coeff = 0.425
            * pow(
                (k_vf * k_vf * k_vf * rho_vf * drho_f * g) / (mu_vf * lambda_c),
                0.25,
            );

        // 8. Build strictly ordered, monotone C^1 knots across all 5 regimes.
        let dt_chf_lo = 0.88 * dt_chf_raw;
        let (q_chf_lo, dq_chf_lo, _) = eval_nucleate_raw(
            dt_chf_lo,
            k_rohsenow,
            props.boiling_range_k,
            schlunder_alpha,
        )?;
        // Guarantee Fritsch-Carlson monotonicity (dq_chf_lo <= 2.5 * secant) on CHF cap.
        let min_chf_span = if dq_chf_lo > 0.0 {
            1.25 * (q_chf - q_chf_lo) / dq_chf_lo
        } else {
            0.18 * dt_chf_raw
        };
        let dt_chf = dt_chf_lo + (0.18 * dt_chf_raw).max(min_chf_span);

        if !(dt_min > dt_chf * 1.05) {
            return Err(BoilingErrorV1::InvalidRegimeOrdering(
                "Leidenfrost superheat T_min must strictly exceed CHF superheat T_chf",
            ));
        }

        // Smooth ONB C^1 blend around dt_onb, strictly inside (0, dt_chf_lo).
        // Using a 0.30 * dt_onb half-width keeps third-derivative curvature bounded
        // across low-ONB binary mixtures (e.g. ethanol-water) while staying strictly monotone.
        let onb_half = (0.30 * dt_onb).min(0.45 * (dt_chf_lo - dt_onb).max(0.0));
        if !(dt_onb > 0.0 && onb_half > 0.0 && dt_onb + onb_half < dt_chf_lo) {
            return Err(BoilingErrorV1::InvalidRegimeOrdering(
                "ONB superheat must lie strictly between saturation and CHF",
            ));
        }
        let dt_onb_lo = dt_onb - onb_half;
        let dt_onb_hi = dt_onb + onb_half;
        let q_onb_lo = c_nc * pow(dt_onb_lo, 4.0 / 3.0);
        let dq_onb_lo = (4.0 / 3.0) * c_nc * pow(dt_onb_lo, 1.0 / 3.0);
        let (q_onb_hi, dq_onb_hi, _) = eval_nucleate_raw(
            dt_onb_hi,
            k_rohsenow,
            props.boiling_range_k,
            schlunder_alpha,
        )?;

        // Leidenfrost valley cup on [dt_min, dt_film_hi] matching Bromley film boiling at dt_film_hi.
        let dt_film_hi = 1.15 * dt_min;
        let (q_film_hi, dq_film_hi) = eval_film_raw(
            dt_film_hi,
            props.saturation_temperature_k,
            film_conv_coeff,
            h_fg,
            props.vapor_cp_j_kg_k,
            surface.surface_emissivity,
        );
        // Fritsch-Carlson cap ensures zero undershoot below q_valley on [dt_min, dt_film_hi].
        let fc_max_valley = q_film_hi - 0.5 * dq_film_hi * (dt_film_hi - dt_min);
        let q_valley = q_min_ber.min(fc_max_valley).min(0.5 * q_chf);
        if !q_valley.is_finite() || q_valley <= 0.0 || q_valley >= q_chf {
            return Err(BoilingErrorV1::InvalidRegimeOrdering(
                "Leidenfrost valley heat flux must lie strictly in (0, q_CHF)",
            ));
        }

        let curve = Self {
            props,
            surface,
            c_nc,
            k_rohsenow,
            schlunder_alpha,
            film_conv_coeff,
            dt_onb,
            dt_onb_lo,
            dt_onb_hi,
            q_onb_lo,
            dq_onb_lo,
            q_onb_hi,
            dq_onb_hi,
            dt_chf_lo,
            q_chf_lo,
            dq_chf_lo,
            dt_chf,
            q_chf,
            dt_min,
            q_valley,
            dt_film_hi,
            q_film_hi,
            dq_film_hi,
        };
        Ok(curve)
    }

    /// Evaluated fluid saturation and transport properties.
    #[must_use]
    pub const fn fluid_properties(&self) -> &BoilingFluidPropertiesV1 {
        &self.props
    }

    /// Retained surface boiling parameters.
    #[must_use]
    pub const fn surface(&self) -> &BoilingSurfaceCardV1 {
        &self.surface
    }

    /// Local saturation / bubble-point temperature `T_sat(P)` [K].
    #[must_use]
    pub const fn saturation_temperature(&self) -> Temperature {
        Temperature::new(self.props.saturation_temperature_k)
    }

    /// Ambient / vessel pressure `P` [Pa].
    #[must_use]
    pub const fn pressure(&self) -> Pressure {
        Pressure::new(self.props.pressure_pa)
    }

    /// Key regime boundary temperatures [K] and peak/valley fluxes [W/m^2].
    #[must_use]
    pub const fn regime_boundaries(&self) -> BoilingRegimeBoundariesV1 {
        let t_sat = self.props.saturation_temperature_k;
        BoilingRegimeBoundariesV1 {
            t_sat_k: t_sat,
            t_onb_k: t_sat + self.dt_onb,
            t_chf_k: t_sat + self.dt_chf,
            t_leidenfrost_k: t_sat + self.dt_min,
            q_chf_w_m2: self.q_chf,
            q_leidenfrost_w_m2: self.q_valley,
        }
    }

    /// Evaluate the `C^1`-smooth Nukiyama boiling heat flux and exact analytical
    /// tangent `dq''/dT_w` at `wall_temperature`.
    #[allow(clippy::too_many_lines)]
    pub fn evaluate(
        &self,
        wall_temperature: Temperature,
    ) -> Result<BoilingEvaluationV1, BoilingErrorV1> {
        let tw = wall_temperature.value();
        if !tw.is_finite() || tw <= 0.0 {
            return Err(BoilingErrorV1::InvalidTemperature { temperature_k: tw });
        }
        let t_sat = self.props.saturation_temperature_k;
        let dt_e = tw - t_sat;

        if dt_e <= 0.0 {
            return Ok(BoilingEvaluationV1 {
                wall_temperature,
                saturation_temperature: Temperature::new(t_sat),
                superheat_k: dt_e,
                regime: BoilingRegimeV1::BelowSaturation,
                heat_flux_w_m2: 0.0,
                htc_w_m2_k: 0.0,
                dq_dtw_w_m2_k: 0.0,
                ideal_rohsenow_flux_w_m2: 0.0,
                mixture_suppression_factor: 1.0,
                vaporization_mass_flux_kg_m2_s: 0.0,
                vapor_enthalpy_janaf_j_kg: self.props.vapor_enthalpy_janaf_j_kg,
                liquid_enthalpy_janaf_j_kg: self.props.liquid_enthalpy_janaf_j_kg,
            });
        }

        let q_id_roh = self.k_rohsenow * dt_e * dt_e * dt_e;
        let (regime, q, dq_dt, suppression) = if dt_e <= self.dt_onb_lo {
            let q_nc = self.c_nc * pow(dt_e, 4.0 / 3.0);
            let dq_nc = (4.0 / 3.0) * self.c_nc * pow(dt_e, 1.0 / 3.0);
            (BoilingRegimeV1::NaturalConvection, q_nc, dq_nc, 1.0)
        } else if dt_e < self.dt_onb_hi {
            let (q_b, dq_b) = cubic_hermite(
                dt_e,
                self.dt_onb_lo,
                self.dt_onb_hi,
                self.q_onb_lo,
                self.dq_onb_lo,
                self.q_onb_hi,
                self.dq_onb_hi,
            );
            let reg = if dt_e <= self.dt_onb {
                BoilingRegimeV1::NaturalConvection
            } else {
                BoilingRegimeV1::NucleateBoiling
            };
            (reg, q_b, dq_b, (q_b / q_id_roh).min(1.0))
        } else if dt_e <= self.dt_chf_lo {
            let (q_n, dq_n, supp) = eval_nucleate_raw(
                dt_e,
                self.k_rohsenow,
                self.props.boiling_range_k,
                self.schlunder_alpha,
            )?;
            (BoilingRegimeV1::NucleateBoiling, q_n, dq_n, supp)
        } else if dt_e <= self.dt_chf {
            let (q_c, dq_c) = cubic_hermite(
                dt_e,
                self.dt_chf_lo,
                self.dt_chf,
                self.q_chf_lo,
                self.dq_chf_lo,
                self.q_chf,
                0.0,
            );
            let (_, _, supp) = eval_nucleate_raw(
                dt_e,
                self.k_rohsenow,
                self.props.boiling_range_k,
                self.schlunder_alpha,
            )?;
            (BoilingRegimeV1::NucleateBoiling, q_c, dq_c, supp)
        } else if dt_e < self.dt_min {
            // Log-flux monotone cubic Hermite spline from (dt_chf, q_chf, 0) to (dt_min, q_valley, 0).
            let span = self.dt_min - self.dt_chf;
            let s = (dt_e - self.dt_chf) / span;
            let h_s = s * s * (3.0 - 2.0 * s);
            let hp_s = 6.0 * s * (1.0 - s);
            let log_ratio = ln(self.q_valley / self.q_chf);
            let q_tr = self.q_chf * exp(h_s * log_ratio);
            let dq_tr = q_tr * (hp_s / span) * log_ratio;
            (BoilingRegimeV1::TransitionBoiling, q_tr, dq_tr, 1.0)
        } else if dt_e < self.dt_film_hi {
            let (q_v, dq_v) = cubic_hermite(
                dt_e,
                self.dt_min,
                self.dt_film_hi,
                self.q_valley,
                0.0,
                self.q_film_hi,
                self.dq_film_hi,
            );
            (BoilingRegimeV1::FilmBoiling, q_v, dq_v, 1.0)
        } else {
            let (q_f, dq_f) = eval_film_raw(
                dt_e,
                t_sat,
                self.film_conv_coeff,
                self.props.latent_heat_j_kg,
                self.props.vapor_cp_j_kg_k,
                self.surface.surface_emissivity,
            );
            (BoilingRegimeV1::FilmBoiling, q_f, dq_f, 1.0)
        };

        if !q.is_finite() || q < 0.0 || !dq_dt.is_finite() {
            return Err(BoilingErrorV1::Unrepresentable(
                "evaluated boiling flux or derivative is non-finite",
            ));
        }

        let htc = q / dt_e;
        let m_dot_vap = match regime {
            BoilingRegimeV1::BelowSaturation | BoilingRegimeV1::NaturalConvection => 0.0,
            BoilingRegimeV1::NucleateBoiling
            | BoilingRegimeV1::TransitionBoiling
            | BoilingRegimeV1::FilmBoiling => q / self.props.latent_heat_j_kg,
        };

        Ok(BoilingEvaluationV1 {
            wall_temperature,
            saturation_temperature: Temperature::new(t_sat),
            superheat_k: dt_e,
            regime,
            heat_flux_w_m2: q,
            htc_w_m2_k: htc,
            dq_dtw_w_m2_k: dq_dt,
            ideal_rohsenow_flux_w_m2: q_id_roh,
            mixture_suppression_factor: suppression,
            vaporization_mass_flux_kg_m2_s: m_dot_vap,
            vapor_enthalpy_janaf_j_kg: self.props.vapor_enthalpy_janaf_j_kg,
            liquid_enthalpy_janaf_j_kg: self.props.liquid_enthalpy_janaf_j_kg,
        })
    }

    /// Recommend a shortened step `dt_s` (and [`Proposal::Shorten`]) whenever a trial
    /// wall-temperature change `t_start -> t_trial` over `trial_dt_s` exceeds
    /// `max_delta_t_k` or jumps across ONB, CHF, or Leidenfrost regime boundaries.
    pub fn recommend_step_cap(
        &self,
        t_start: Temperature,
        t_trial: Temperature,
        trial_dt_s: f64,
        max_delta_t_k: f64,
    ) -> Result<BoilingStepAdviceV1, BoilingErrorV1> {
        let t0 = t_start.value();
        let t1 = t_trial.value();
        if !t0.is_finite() || t0 <= 0.0 || !t1.is_finite() || t1 <= 0.0 {
            return Err(BoilingErrorV1::InvalidTemperature {
                temperature_k: if !t0.is_finite() || t0 <= 0.0 { t0 } else { t1 },
            });
        }
        if !trial_dt_s.is_finite()
            || trial_dt_s <= 0.0
            || !max_delta_t_k.is_finite()
            || max_delta_t_k <= 0.0
        {
            return Err(BoilingErrorV1::InvalidParameter(
                "trial_dt_s and max_delta_t_k must be finite and strictly positive",
            ));
        }

        let delta = (t1 - t0).abs();
        let mut factor = 1.0_f64;

        if delta > max_delta_t_k {
            factor = factor.min(max_delta_t_k / delta);
        }

        // Check whether [t0, t1] straddles any Nukiyama regime boundary (ONB, CHF, Leidenfrost).
        let b = self.regime_boundaries();
        let lo = t0.min(t1);
        let hi = t0.max(t1);
        for boundary in [b.t_sat_k, b.t_onb_k, b.t_chf_k, b.t_leidenfrost_k] {
            if boundary > lo && boundary < hi && delta > 1e-9 {
                let dist = (boundary - t0).abs();
                // Land slightly past or just shy of the knot only if step is wide; cap at 0.5*dist/delta + small floor.
                if delta > 0.5 * max_delta_t_k {
                    let cross_frac = ((dist / delta) * 0.85).clamp(0.10, 0.95);
                    factor = factor.min(cross_frac);
                }
            }
        }

        if factor < 1.0 {
            let rec_dt = (trial_dt_s * factor).max(trial_dt_s * 1e-3);
            Ok(BoilingStepAdviceV1 {
                shorten_recommended: true,
                recommended_max_dt_s: rec_dt,
                shorten_proposal: Some(Proposal::Shorten { dt_s: rec_dt }),
            })
        } else {
            Ok(BoilingStepAdviceV1 {
                shorten_recommended: false,
                recommended_max_dt_s: trial_dt_s,
                shorten_proposal: None,
            })
        }
    }
}

/// Evaluate exact implicit Schluender (1983) / Thome-Shakir (1987) nucleate boiling:
/// `u + Delta T_bp * (1 - exp(-alpha * u^3)) = Delta T_e`, with `q'' = K_Roh * u^3`.
///
/// Returns `(q_flux, dq_dT, suppression_factor = q / q_ideal)`.
fn eval_nucleate_raw(
    dt_e: f64,
    k_rohsenow: f64,
    dt_bp: f64,
    alpha: f64,
) -> Result<(f64, f64, f64), BoilingErrorV1> {
    if dt_bp == 0.0 {
        let q = k_rohsenow * dt_e * dt_e * dt_e;
        let dq = 3.0 * k_rohsenow * dt_e * dt_e;
        return Ok((q, dq, 1.0));
    }

    // Solve strictly monotone F(u) = u + dt_bp*(1 - exp(-alpha*u^3)) - dt_e = 0 on (0, dt_e].
    let mut lo = 0.0_f64;
    let mut hi = dt_e;
    let mut u = dt_e / (1.0 + (dt_bp / dt_e).min(2.0));
    for _ in 0..48 {
        let u2 = u * u;
        let u3 = u2 * u;
        let arg = alpha * u3;
        let one_minus_e = -expm1(-arg);
        let e = exp(-arg);
        let f = u + dt_bp * one_minus_e - dt_e;
        let fp = 1.0 + 3.0 * alpha * dt_bp * u2 * e;
        if f.abs() <= 1e-15 * dt_e.max(1.0) {
            break;
        }
        if f < 0.0 {
            lo = u;
        } else {
            hi = u;
        }
        let next = u - f / fp;
        u = if next.is_finite() && next > lo && next < hi {
            next
        } else {
            0.5 * (lo + hi)
        };
    }

    let u2 = u * u;
    let u3 = u2 * u;
    let e = exp(-alpha * u3);
    let du_ddt = 1.0 / (1.0 + 3.0 * alpha * dt_bp * u2 * e);
    let q = k_rohsenow * u3;
    let dq = 3.0 * k_rohsenow * u2 * du_ddt;
    let q_id = k_rohsenow * dt_e * dt_e * dt_e;
    let supp = if q_id > 0.0 { q / q_id } else { 1.0 };

    if !q.is_finite() || !dq.is_finite() || !supp.is_finite() {
        return Err(BoilingErrorV1::Unrepresentable(
            "Schluender-Thome nucleate boiling evaluation became non-finite",
        ));
    }
    Ok((q, dq, supp))
}

fn solve_nucleate_superheat_for_flux(
    q_target: f64,
    k_rohsenow: f64,
    dt_bp: f64,
    alpha: f64,
) -> Result<f64, BoilingErrorV1> {
    if !q_target.is_finite() || q_target <= 0.0 || !k_rohsenow.is_finite() || k_rohsenow <= 0.0 {
        return Err(BoilingErrorV1::Unrepresentable(
            "invalid CHF or Rohsenow constant",
        ));
    }
    // Exact inverse! Since q = K_Roh * u^3, u = (q_target / K_Roh)^(1/3), and
    // Delta T_e = u + Delta T_bp * (1 - exp(-alpha * u^3))!
    let u = pow(q_target / k_rohsenow, 1.0 / 3.0);
    let dt_e = u + dt_bp * (-expm1(-alpha * u * u * u));
    if !dt_e.is_finite() || dt_e <= 0.0 {
        return Err(BoilingErrorV1::Unrepresentable(
            "CHF superheat is non-positive or non-finite",
        ));
    }
    Ok(dt_e)
}

fn solve_onb_superheat(
    c_nc: f64,
    k_rohsenow: f64,
    dt_bp: f64,
    alpha: f64,
    dt_chf_raw: f64,
) -> Result<f64, BoilingErrorV1> {
    // Find unique root of q_nuc(dT) - q_nc(dT) = 0 on (0, dt_chf_raw).
    let mut lo = 1e-6_f64.min(0.01 * dt_chf_raw);
    let mut hi = dt_chf_raw;
    let f_at = |dt: f64| -> Result<f64, BoilingErrorV1> {
        let (q_n, _, _) = eval_nucleate_raw(dt, k_rohsenow, dt_bp, alpha)?;
        let q_nc = c_nc * pow(dt, 4.0 / 3.0);
        Ok(q_n - q_nc)
    };
    if f_at(lo)? >= 0.0 {
        return Ok(lo);
    }
    if f_at(hi)? <= 0.0 {
        return Err(BoilingErrorV1::InvalidRegimeOrdering(
            "natural convection exceeds CHF across nucleate domain",
        ));
    }
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if (hi - lo) <= 1e-13 * mid.max(1.0) {
            return Ok(mid);
        }
        if f_at(mid)? < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Ok(0.5 * (lo + hi))
}

fn eval_film_raw(
    dt_e: f64,
    t_sat: f64,
    film_conv_coeff: f64,
    h_fg: f64,
    cp_v: f64,
    emissivity: f64,
) -> (f64, f64) {
    let h_fg_prime = h_fg + 0.5 * cp_v * dt_e;
    let q_conv = film_conv_coeff * pow(dt_e, 0.75) * pow(h_fg_prime, 0.25);
    let dq_conv = q_conv * (0.75 / dt_e + (0.125 * cp_v) / h_fg_prime);

    let tw = t_sat + dt_e;
    let tw2 = tw * tw;
    let ts2 = t_sat * t_sat;
    let q_rad = 0.75 * emissivity * STEFAN_BOLTZMANN_W_M2_K4 * (tw2 * tw2 - ts2 * ts2);
    let dq_rad = 3.0 * emissivity * STEFAN_BOLTZMANN_W_M2_K4 * tw2 * tw;

    (q_conv + q_rad, dq_conv + dq_rad)
}

/// Evaluate exact `C^1` cubic Hermite interpolant `(value, derivative)` on `[x0, x1]`.
fn cubic_hermite(x: f64, x0: f64, x1: f64, y0: f64, m0: f64, y1: f64, m1: f64) -> (f64, f64) {
    let dx = x1 - x0;
    let s = (x - x0) / dx;
    let s2 = s * s;
    let s3 = s2 * s;
    let h00 = 2.0 * s3 - 3.0 * s2 + 1.0;
    let h10 = s3 - 2.0 * s2 + s;
    let h01 = -2.0 * s3 + 3.0 * s2;
    let h11 = s3 - s2;
    let val = h00 * y0 + h10 * (dx * m0) + h01 * y1 + h11 * (dx * m1);

    let dh00 = 6.0 * s2 - 6.0 * s;
    let dh10 = 3.0 * s2 - 4.0 * s + 1.0;
    let dh01 = -6.0 * s2 + 6.0 * s;
    let dh11 = 3.0 * s2 - 2.0 * s;
    let deriv = (dh00 * y0 + dh01 * y1) / dx + dh10 * m0 + dh11 * m1;
    (val, deriv)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn water_1atm_chf_and_berenson_henry_leidenfrost_match_textbook() {
        let curve = BoilingCurveV1::prepare(
            BoilingFluidKindV1::PureWater,
            BoilingSurfaceCardV1::stainless_steel_304_v1(),
            Pressure::new(101_325.0),
        )
        .unwrap();
        let b = curve.regime_boundaries();
        assert!((b.q_chf_w_m2 / 1e6 - 1.261).abs() < 0.02);
        let t_min_c = b.t_leidenfrost_k - 273.15;
        assert!(
            (195.0..=225.0).contains(&t_min_c),
            "expected saturated pool Leidenfrost near 195-220 C, got {t_min_c:.2} C"
        );
    }

    #[test]
    fn analytic_tangent_matches_central_fd_across_all_five_regimes() {
        let curve = BoilingCurveV1::prepare(
            BoilingFluidKindV1::EthanolWaterNrtl {
                ethanol_mole_fraction: 0.05,
            },
            BoilingSurfaceCardV1::seasoned_cast_iron_v1(),
            Pressure::new(101_325.0),
        )
        .unwrap();
        let t_sat = curve.saturation_temperature().value();
        let b = curve.regime_boundaries();

        let test_temps = [
            t_sat + 1.5,
            b.t_onb_k,
            t_sat + 12.0,
            b.t_chf_k - 1.5,
            0.5 * (b.t_chf_k + b.t_leidenfrost_k),
            b.t_leidenfrost_k + 4.0,
            b.t_leidenfrost_k + 80.0,
        ];
        for tw in test_temps {
            let ev = curve.evaluate(Temperature::new(tw)).unwrap();
            let h = 1e-4;
            let qp = curve
                .evaluate(Temperature::new(tw + h))
                .unwrap()
                .heat_flux_w_m2;
            let qm = curve
                .evaluate(Temperature::new(tw - h))
                .unwrap()
                .heat_flux_w_m2;
            let fd = (qp - qm) / (2.0 * h);
            let rel = (ev.dq_dtw_w_m2_k - fd).abs() / fd.abs().max(1.0);
            assert!(
                rel < 1e-8,
                "dq/dT mismatch at Tw={tw:.2} K ({:?}): analytic={}, fd={fd}, rel={rel:e}",
                ev.regime,
                ev.dq_dtw_w_m2_k
            );
        }
    }

    #[test]
    fn triglyceride_oil_is_explicitly_refused() {
        let res = BoilingCurveV1::prepare(
            BoilingFluidKindV1::TriglycerideOil {
                name: "refined_canola_oil",
                smoke_point_k: 493.15,
                flash_point_k: 588.15,
            },
            BoilingSurfaceCardV1::stainless_steel_304_v1(),
            Pressure::new(101_325.0),
        );
        assert!(matches!(
            res,
            Err(BoilingErrorV1::NonBoilingLipidPhase { .. })
        ));
    }
}
