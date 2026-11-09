//! Sourced horological material cards and constitutive evaluators for the
//! FrankenSim Grand Complication Pocket Watch (`Calibre FS-1903 Tourbillon`).
//!
//! All properties trace to `data/matdb/seed-v1` packs:
//! - `gold-pure-nbs-c447` & `gold-solid-nbs1989-thermochemistry`
//! - `brass-c36000-cda-reference` & `brass-c26000-nist-cryogenic`
//! - `music-wire-spring-moduli-fuchs-1968` & `aisi-52100-cvm-hot-hardness`
//! - `beryllium-copper-c17200-ngk-reference` & `beryllium-copper-nist-thermal-reference`
//! - `invar-fe36ni-nist-thermal-reference` & `nilo-36-special-metals-2004`
//! - `titanium-ti6al4v-nist-thermal-reference`
//! - `sapphire-al2o3-crystran-reference` & `alumina-alpha-nist-shomate-caloric`
//! - `pennzane-shf-x-2000-bearing-oil` & `napc-pe-5-l-1307-1553-gear-oil`
//! - `nitrogen-transport-lemmon2004` & `nitrogen-nasa9-mcbride2002`
//! - `jinshan-n42-pristine-temperature` & `ngyc-n42-sintered-nickel-coated`

use fs_blake3::DomainHasher;
pub use fs_conduction::STEFAN_BOLTZMANN_W_M2_K4 as STEFAN_BOLTZMANN;
pub use fs_matdb::{LawId, MATDB_SCHEMA_VERSION};
pub use fs_material::creep::plastic::{NortonJ2, NortonJ2State};
pub use fs_material::creep::{CreepState, NortonCreep, NortonCreepParameters};
pub use fs_material::elastic::IsotropicElastic;
pub use fs_material::plastic::{J2Plasticity, J2State};
use fs_render::conductor::{
    CONDUCTOR_IOR_SAMPLE_COUNT, CoherentDielectricCoating, ConductorDataStatus, ConductorIorSample,
    ConductorOptics, ConductorSource, fresnel_conductor,
};
use fs_render::dielectric::CauchyIor;
pub use fs_thermal::enthalpy::{DepositionError, ThermalLaw, ThermalStorageBasis};

/// Reference ambient laboratory temperature $T_0 = 293.15\text{ K}$ ($20.0^\circ\text{C}$).
pub const T_REF_K: f64 = 293.15;

/// Liquid nitrogen normal boiling point at $101{,}325\text{ Pa}$ ($77.36\text{ K}$).
pub const T_LN2_SAT_K: f64 = 77.36;

/// Liquid nitrogen latent heat of vaporization $L_v$ ($\text{J/kg}$).
pub const LN2_LATENT_VAP_J_PER_KG: f64 = 199_180.0;

/// Leidenfrost transition temperature for $\text{LN}_2$ on metal surfaces ($\text{K}$).
pub const T_LN2_LEIDENFROST_K: f64 = 126.0;

/// Canonical material identifiers used across the pocket-watch assembly.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum WatchMaterialId {
    /// 18k/24k Horological Gold (`gold-pure-nbs-c447` / `gold-solid-nbs1989-thermochemistry`).
    Gold18kCase = 0,
    /// Watchmaker's Free-Cutting / Cartridge Brass (`brass-c36000-cda-reference`).
    BrassC36000 = 1,
    /// High-Carbon Cold-Drawn Music Wire / 52100 Steel (`music-wire-spring-moduli-fuchs-1968`).
    HighCarbonSteel = 2,
    /// Beryllium-Copper Glucydur Balance Alloy (`beryllium-copper-c17200-ngk-reference`).
    GlucydurBeCu = 3,
    /// Invar / Nilo-36 Nickel-Iron Hairspring Alloy (`invar-fe36ni-nist-thermal-reference`).
    Invar36Ni = 4,
    /// Titanium Grade 5 Tourbillon Cage (`titanium-ti6al4v-nist-thermal-reference`).
    TitaniumGrade5 = 5,
    /// Monocrystalline Corundum Ruby / Sapphire $\alpha\text{-Al}_2\text{O}_3$ (`sapphire-al2o3-crystran-reference`).
    SyntheticCorundumRuby = 6,
}

/// Complete thermophysical, mechanical, phase-change, magnetic, and optical card.
#[derive(Debug, Clone, Copy)]
pub struct MaterialCard {
    pub id: WatchMaterialId,
    pub matdb_pack_id: &'static str,
    pub display_name: &'static str,
    /// Reference mass density $\rho_0$ at $293.15\text{ K}$ ($\text{kg/m}^3$).
    pub density_kg_m3: f64,
    /// Specific heat capacity $c_{p,0}$ at $293.15\text{ K}$ ($\text{J/(kg}\cdot\text{K)}$).
    pub cp_j_kg_k: f64,
    /// Thermal conductivity $k_0$ at $293.15\text{ K}$ ($\text{W/(m}\cdot\text{K)}$).
    pub conductivity_w_m_k: f64,
    /// Linear coefficient of thermal expansion $\alpha_0$ at $293.15\text{ K}$ ($\text{K}^{-1}$).
    pub cte_per_k: f64,
    /// Solidus temperature $T_{\text{sol}}$ ($\text{K}$).
    pub t_solidus_k: f64,
    /// Liquidus / melting temperature $T_{\text{liq}}$ ($\text{K}$).
    pub t_liquidus_k: f64,
    /// Latent heat of fusion $L_f$ ($\text{J/kg}$).
    pub latent_fusion_j_kg: f64,
    /// Reference Young's modulus $E_0$ at $293.15\text{ K}$ ($\text{Pa}$).
    pub youngs_modulus_pa: f64,
    /// Thermoelastic coefficient $\beta_E = \frac{1}{E_0}\frac{dE}{dT}$ ($\text{K}^{-1}$).
    pub thermoelastic_beta_per_k: f64,
    /// Relative magnetic susceptibility $\chi_m = \mu_r - 1$ (dimensionless).
    pub magnetic_susceptibility: f64,
    /// Saturation magnetization $B_{\text{sat}}$ ($\text{T}$).
    pub saturation_flux_t: f64,
    /// Electrical conductivity $\sigma_e$ at $293.15\text{ K}$ ($\text{S/m}$) for eddy-current braking.
    pub electrical_conductivity_s_m: f64,
    /// Complex index of refraction $(n, \kappa)$ at $(630\text{ nm}, 532\text{ nm}, 465\text{ nm})$.
    pub optical_n_rgb: [f64; 3],
    pub optical_k_rgb: [f64; 3],
    /// Hemispherical thermal emissivity $\varepsilon_0$.
    pub emissivity: f64,
}

pub const GOLD_18K_CARD: MaterialCard = MaterialCard {
    id: WatchMaterialId::Gold18kCase,
    matdb_pack_id: "gold-pure-nbs-c447 / gold-solid-nbs1989-thermochemistry",
    display_name: "18k Yellow Gold (750 Au-Ag-Cu) / NBS-C447",
    density_kg_m3: 19_320.0,
    cp_j_kg_k: 129.1,
    conductivity_w_m_k: 317.0,
    cte_per_k: 14.2e-6,
    t_solidus_k: 1215.0,
    t_liquidus_k: 1337.33,
    latent_fusion_j_kg: 63_700.0,
    youngs_modulus_pa: 79.0e9,
    thermoelastic_beta_per_k: -3.1e-4,
    magnetic_susceptibility: -3.4e-5,
    saturation_flux_t: 0.0,
    electrical_conductivity_s_m: 4.10e7,
    // True complex index of Au at (630nm, 532nm, 465nm)
    optical_n_rgb: [0.183, 0.467, 1.380],
    optical_k_rgb: [3.430, 2.410, 1.910],
    emissivity: 0.07,
};

pub const BRASS_C36000_CARD: MaterialCard = MaterialCard {
    id: WatchMaterialId::BrassC36000,
    matdb_pack_id: "brass-c36000-cda-reference / brass-c26000-nist-cryogenic",
    display_name: "Horological Gilt/Rhodium Brass C36000",
    density_kg_m3: 8_500.0,
    cp_j_kg_k: 380.0,
    conductivity_w_m_k: 115.0,
    cte_per_k: 20.5e-6,
    t_solidus_k: 1158.0,
    t_liquidus_k: 1188.0,
    latent_fusion_j_kg: 168_000.0,
    youngs_modulus_pa: 110.0e9,
    thermoelastic_beta_per_k: -3.8e-4,
    magnetic_susceptibility: -1.0e-5,
    saturation_flux_t: 0.0,
    electrical_conductivity_s_m: 1.59e7,
    optical_n_rgb: [0.440, 0.560, 1.120],
    optical_k_rgb: [3.520, 2.750, 2.150],
    emissivity: 0.10,
};

pub const HIGH_CARBON_STEEL_CARD: MaterialCard = MaterialCard {
    id: WatchMaterialId::HighCarbonSteel,
    matdb_pack_id: "music-wire-spring-moduli-fuchs-1968 / aisi-52100-cvm-hot-hardness",
    display_name: "High-Carbon Tempered Spring & Arbor Steel (ASTM A228 / 52100)",
    density_kg_m3: 7_850.0,
    cp_j_kg_k: 475.0,
    conductivity_w_m_k: 46.6,
    cte_per_k: 11.5e-6,
    t_solidus_k: 1670.0,
    t_liquidus_k: 1725.0,
    latent_fusion_j_kg: 270_000.0,
    youngs_modulus_pa: 207.0e9,
    thermoelastic_beta_per_k: -2.65e-4,
    magnetic_susceptibility: 185.0,
    saturation_flux_t: 1.65,
    electrical_conductivity_s_m: 5.80e6,
    optical_n_rgb: [2.880, 2.710, 2.520],
    optical_k_rgb: [3.350, 3.180, 3.010],
    emissivity: 0.24,
};

pub const GLUCYDUR_BECU_CARD: MaterialCard = MaterialCard {
    id: WatchMaterialId::GlucydurBeCu,
    matdb_pack_id: "beryllium-copper-c17200-ngk-reference / beryllium-copper-nist-thermal-reference",
    display_name: "Glucydur Beryllium-Bronze Balance Alloy (C17200)",
    density_kg_m3: 8_250.0,
    cp_j_kg_k: 420.0,
    conductivity_w_m_k: 105.0,
    cte_per_k: 17.0e-6,
    t_solidus_k: 1140.0,
    t_liquidus_k: 1255.0,
    latent_fusion_j_kg: 195_000.0,
    youngs_modulus_pa: 131.0e9,
    thermoelastic_beta_per_k: -3.2e-4,
    magnetic_susceptibility: -1.2e-5,
    saturation_flux_t: 0.0,
    electrical_conductivity_s_m: 1.28e7,
    optical_n_rgb: [0.320, 0.680, 1.240],
    optical_k_rgb: [3.410, 2.580, 2.090],
    emissivity: 0.09,
};

pub const INVAR_36NI_CARD: MaterialCard = MaterialCard {
    id: WatchMaterialId::Invar36Ni,
    matdb_pack_id: "invar-fe36ni-nist-thermal-reference / nilo-36-special-metals-2004",
    display_name: "Compensated Invar / Elinvar Fe-36Ni Hairspring Alloy",
    density_kg_m3: 8_050.0,
    cp_j_kg_k: 515.0,
    conductivity_w_m_k: 13.4,
    cte_per_k: 1.2e-6,
    t_solidus_k: 1695.0,
    t_liquidus_k: 1720.0,
    latent_fusion_j_kg: 265_000.0,
    youngs_modulus_pa: 148.0e9,
    // Near-zero first-order thermoelastic coefficient near 293 K, with quadratic Dent middle-temperature error
    thermoelastic_beta_per_k: 1.8e-5,
    magnetic_susceptibility: 95.0,
    saturation_flux_t: 1.20,
    electrical_conductivity_s_m: 1.22e6,
    optical_n_rgb: [2.620, 2.480, 2.310],
    optical_k_rgb: [3.620, 3.410, 3.190],
    emissivity: 0.22,
};

pub const TITANIUM_G5_CARD: MaterialCard = MaterialCard {
    id: WatchMaterialId::TitaniumGrade5,
    matdb_pack_id: "titanium-ti6al4v-nist-thermal-reference",
    display_name: "Grade 5 Titanium Tourbillon Carriage (Ti-6Al-4V)",
    density_kg_m3: 4_430.0,
    cp_j_kg_k: 526.0,
    conductivity_w_m_k: 6.7,
    cte_per_k: 8.6e-6,
    t_solidus_k: 1877.0,
    t_liquidus_k: 1933.0,
    latent_fusion_j_kg: 292_000.0,
    youngs_modulus_pa: 114.0e9,
    thermoelastic_beta_per_k: -4.1e-4,
    magnetic_susceptibility: 1.8e-4,
    saturation_flux_t: 0.0,
    electrical_conductivity_s_m: 5.80e5,
    optical_n_rgb: [2.710, 2.540, 2.390],
    optical_k_rgb: [3.310, 3.150, 2.970],
    emissivity: 0.31,
};

pub const SYNTHETIC_RUBY_CARD: MaterialCard = MaterialCard {
    id: WatchMaterialId::SyntheticCorundumRuby,
    matdb_pack_id: "sapphire-al2o3-crystran-reference / alumina-alpha-nist-shomate-caloric",
    display_name: "Verneuil Monocrystalline Corundum Ruby / Sapphire (alpha-Al2O3)",
    density_kg_m3: 3_980.0,
    cp_j_kg_k: 775.0,
    conductivity_w_m_k: 35.0,
    cte_per_k: 5.4e-6,
    t_solidus_k: 2315.0,
    t_liquidus_k: 2326.0,
    latent_fusion_j_kg: 1_090_000.0,
    youngs_modulus_pa: 345.0e9,
    thermoelastic_beta_per_k: -1.1e-4,
    magnetic_susceptibility: -1.8e-5,
    saturation_flux_t: 0.0,
    electrical_conductivity_s_m: 1.0e-12,
    optical_n_rgb: [1.766, 1.772, 1.778],
    optical_k_rgb: [0.0, 0.0, 0.0],
    emissivity: 0.48,
};

impl WatchMaterialId {
    #[must_use]
    pub const fn card(self) -> &'static MaterialCard {
        match self {
            Self::Gold18kCase => &GOLD_18K_CARD,
            Self::BrassC36000 => &BRASS_C36000_CARD,
            Self::HighCarbonSteel => &HIGH_CARBON_STEEL_CARD,
            Self::GlucydurBeCu => &GLUCYDUR_BECU_CARD,
            Self::Invar36Ni => &INVAR_36NI_CARD,
            Self::TitaniumGrade5 => &TITANIUM_G5_CARD,
            Self::SyntheticCorundumRuby => &SYNTHETIC_RUBY_CARD,
        }
    }
}

impl MaterialCard {
    /// Specific enthalpy $h(T)$ relative to $T_{\text{ref}} = 293.15\text{ K}$ ($\text{J/kg}$),
    /// including sensible heat and mushy-zone latent heat of fusion (`fs-material::phase::melting`).
    #[must_use]
    pub fn enthalpy_from_temperature_j_kg(&self, temp_k: f64) -> f64 {
        let t = temp_k.clamp(20.0, 3500.0);
        let sensible = self.cp_j_kg_k * (t - T_REF_K);
        if t <= self.t_solidus_k {
            sensible
        } else if t >= self.t_liquidus_k {
            sensible + self.latent_fusion_j_kg
        } else {
            let span = (self.t_liquidus_k - self.t_solidus_k).max(1.0);
            let phi = ((t - self.t_solidus_k) / span).clamp(0.0, 1.0);
            sensible + phi * self.latent_fusion_j_kg
        }
    }

    /// Inverse enthalpy map $h \mapsto (T(h), \phi_{\text{liq}}(h))$ where $\phi_{\text{liq}} \in [0, 1]$
    /// is the liquid melt fraction across the mushy solidus-liquidus interval (`MushyPhaseLaw`).
    #[must_use]
    pub fn temperature_and_melt_from_enthalpy(&self, enthalpy_j_kg: f64) -> (f64, f64) {
        let h_sol = self.cp_j_kg_k * (self.t_solidus_k - T_REF_K);
        let span_t = (self.t_liquidus_k - self.t_solidus_k).max(1.0);
        let h_liq = self.cp_j_kg_k * (self.t_liquidus_k - T_REF_K) + self.latent_fusion_j_kg;

        if enthalpy_j_kg <= h_sol {
            let t = (T_REF_K + enthalpy_j_kg / self.cp_j_kg_k).max(20.0);
            (t, 0.0)
        } else if enthalpy_j_kg >= h_liq {
            let t = self.t_liquidus_k + (enthalpy_j_kg - h_liq) / self.cp_j_kg_k;
            (t.min(3500.0), 1.0)
        } else {
            let phi = ((enthalpy_j_kg - h_sol) / (h_liq - h_sol)).clamp(0.0, 1.0);
            let t = self.t_solidus_k + phi * span_t;
            (t, phi)
        }
    }

    /// Carman-Kozeny mushy-zone fluid mobility $\mathcal{M}(\phi_{\text{liq}}) \in [0, 1]$
    /// (`fs-material::phase::mushy`). Zero in solid state; rises rapidly past percolation.
    #[must_use]
    pub fn carman_kozeny_mobility(melt_fraction: f64) -> f64 {
        let phi = melt_fraction.clamp(0.0, 1.0);
        if phi <= 1.0e-4 {
            return 0.0;
        }
        let num = phi * phi * phi;
        let den = num + 0.15 * (1.0 - phi) * (1.0 - phi);
        (num / den).clamp(0.0, 1.0)
    }

    /// Integrated thermal expansion strain $\epsilon_{\text{th}}(T) = \int_{T_0}^T \alpha(T')\,dT'$,
    /// capturing cryogenic Debye flattening below $150\text{ K}$ and high-temperature nonlinearity.
    #[must_use]
    pub fn thermal_strain(&self, temp_k: f64) -> f64 {
        let dt = temp_k - T_REF_K;
        // Debye suppression factor at cryogenic temperatures (77 K LN2 regime)
        let cryo_scale = if temp_k < T_REF_K {
            0.68 + 0.32 * (temp_k / T_REF_K).clamp(0.1, 1.0)
        } else {
            1.0 + 1.8e-4 * (dt / 100.0).min(12.0)
        };
        // Invar loses its low-expansion anomaly above its Curie point (553 K)
        if matches!(self.id, WatchMaterialId::Invar36Ni) && temp_k > 553.0 {
            let below_curie = self.cte_per_k * (553.0 - T_REF_K);
            let above_curie = 13.2e-6 * (temp_k - 553.0);
            return below_curie + above_curie;
        }
        self.cte_per_k * dt * cryo_scale
    }

    /// Temperature-dependent Young's modulus ratio $E(T)/E_0$, including second-order
    /// Dent middle-temperature curvature and Curie softening for Invar.
    #[must_use]
    pub fn youngs_modulus_ratio(&self, temp_k: f64) -> f64 {
        let dt = temp_k - T_REF_K;
        if matches!(self.id, WatchMaterialId::Invar36Ni) {
            if temp_k <= 553.0 {
                // Classic Dent quadratic anomaly (-2.2e-7 K^-2) + small linear term
                let ratio = 1.0 + self.thermoelastic_beta_per_k * dt - 2.2e-7 * dt * dt;
                return ratio.clamp(0.25, 1.15);
            }
            // Above Curie point (553 K), Invar reverts to normal ferromagnetic-loss steel slope
            let at_curie = 1.0 + self.thermoelastic_beta_per_k * (553.0 - T_REF_K)
                - 2.2e-7 * (553.0 - T_REF_K) * (553.0 - T_REF_K);
            return (at_curie - 2.8e-4 * (temp_k - 553.0)).clamp(0.15, 1.15);
        }
        let quad = -6.5e-8 * dt * dt;
        (1.0 + self.thermoelastic_beta_per_k * dt + quad).clamp(0.10, 1.18)
    }

    /// Canonical `fs_matdb::LawId` and `MATDB_SCHEMA_VERSION` for this material card.
    #[must_use]
    pub fn matdb_law_identity(&self) -> (LawId, u32) {
        (LawId(self.matdb_pack_id.to_string()), MATDB_SCHEMA_VERSION)
    }
}

impl ThermalLaw for MaterialCard {
    fn storage_basis(&self) -> ThermalStorageBasis {
        ThermalStorageBasis::Enthalpy
    }

    fn temperature_k(&self, specific_h_j_kg: f64) -> Result<f64, DepositionError> {
        Ok(self.temperature_and_melt_from_enthalpy(specific_h_j_kg).0)
    }

    fn bulk_density_kg_m3(&self, _specific_h_j_kg: f64) -> Result<Option<f64>, DepositionError> {
        Ok(Some(self.density_kg_m3))
    }

    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / self.cp_j_kg_k.max(1.0)
    }
}

/// Dynamic viscosity $\mu(T)$ ($\text{Pa}\cdot\text{s}$) of horological bearing oil
/// (`pennzane-shf-x-2000-bearing-oil` / `napc-pe-5-l-1307-1553-gear-oil`) via
/// Vogel-Fulcher-Tammann (VFT) equation.
///
/// - At $T = 293.15\text{ K}$ ($20^\circ\text{C}$): $\mu \approx 0.082\text{ Pa}\cdot\text{s}$ ($82\text{ cP}$).
/// - At $T = 350\text{ K}$ ($77^\circ\text{C}$): $\mu \approx 0.011\text{ Pa}\cdot\text{s}$ (thinned oil).
/// - Below $T_{\text{pour}} = 215\text{ K}$ down to $77.36\text{ K}$ ($\text{LN}_2$): $\mu$ spikes by $>10^5\times$ (glassy solid wax).
#[must_use]
pub fn watch_oil_viscosity_pa_s(temp_k: f64) -> f64 {
    let t_vft = 158.0;
    let b_vft = 895.0;
    let mu_inf = 1.10e-4;
    let denom = (temp_k - t_vft).max(11.5);
    (mu_inf * (b_vft / denom).exp()).clamp(1.0e-4, 5.0e4)
}

/// Surface tension $\gamma(T)$ ($\text{N/m}$) of watch oil for capillary meniscus bridges
/// (`fs-contact::capillary_bridge`).
#[must_use]
pub fn watch_oil_surface_tension_n_m(temp_k: f64) -> f64 {
    (0.0315 - 7.2e-5 * (temp_k - T_REF_K)).clamp(0.008, 0.048)
}

/// Eötvös temperature derivative $\frac{d\sigma_{\text{Au}}}{dT}$ ($\text{N/(m}\cdot\text{K)}$)
/// of molten 18k Gold (`gold-pure-nbs-c447`), driving thermocapillary Marangoni shear flow.
pub const GOLD_MARANGONI_DSIGMA_DT_N_M_K: f64 = -5.2e-4;

/// Surface tension $\sigma_{\text{Au}}(T)$ ($\text{N/m}$) of molten 18k Gold via the Eötvös / Guggenheim law:
/// $\sigma_{\text{Au}}(T) = 1.145 - 5.2\times 10^{-4} (T - 1337.33\text{ K})\text{ N/m}$.
#[must_use]
pub fn gold_liquid_surface_tension_n_m(temp_k: f64) -> f64 {
    let dt = (temp_k - GOLD_18K_CARD.t_solidus_k).max(0.0);
    (1.145 + GOLD_MARANGONI_DSIGMA_DT_N_M_K * dt).clamp(0.45, 1.25)
}

/// Temperature-dependent equilibrium wetting cosine $\cos\theta_c(T_{\text{sub}})$ of molten 18k Gold
/// on horological brass / steel / sapphire substrates:
/// - Below $900\text{ K}$, surface oxide films make molten gold strongly non-wetting ($\theta_c \approx 118^\circ$, $\cos\theta_c \approx -0.47$),
///   causing liquid gold to ball up into high-curvature beads with capillary height $h_{\text{cap}} \approx 2.4\text{--}3.8\text{ mm}$.
/// - Above $1050\text{ K}$, thermal dissociation & intermetallic diffusion trigger reactive wetting ($\theta_c \approx 50^\circ$, $\cos\theta_c \approx +0.64$).
#[must_use]
pub fn gold_substrate_contact_angle_cos(substrate_temp_k: f64) -> f64 {
    let xi = ((substrate_temp_k - 900.0) / 280.0).clamp(0.0, 1.0);
    let smooth = xi * xi * (3.0 - 2.0 * xi);
    -0.47 + 1.11 * smooth
}

/// Dynamic viscosity $\mu_{\text{eff}}(T, f_L)$ ($\text{Pa}\cdot\text{s}$) of molten and mushy-zone 18k Gold:
/// - Above liquidus ($f_L = 1$): Andrade Arrhenius equation $\mu_{\text{liq}}(T) = \mu_0 \exp(E_\mu / (R T)) \approx 5.1\times 10^{-3}\text{ Pa}\cdot\text{s}$ at $T_m$.
/// - In the mushy solidus-liquidus interval ($0 < f_L < 1$): Krieger-Dougherty / Carman-Kozeny suspension viscosity
///   $\mu_{\text{eff}}(T, f_L) = \mu_{\text{liq}}(T) \cdot \max(f_L, 0.02)^{-2.8}$.
#[must_use]
pub fn gold_effective_viscosity_pa_s(temp_k: f64, melt_fraction: f64) -> f64 {
    let t_eff = temp_k.max(300.0);
    // Andrade parameters for liquid Au: E_mu / R = 2650 K, mu(1337 K) = 5.13e-3 Pa*s
    let mu_liq = 7.08e-4 * (2650.0 / t_eff).exp();
    let fl = melt_fraction.clamp(0.015, 1.0);
    (mu_liq * fl.powf(-2.8)).clamp(2.5e-3, 1.0e5)
}

/// Johnson-Cook + Ludwik thermo-viscoplastic yield stress $\sigma_y(T, f_L, \bar{\varepsilon}^p)$ ($\text{Pa}$)
/// of 18k Gold (`Au-750`):
/// $$\sigma_y = \left(\sigma_{y0} + K_{\text{hard}} (\bar{\varepsilon}^p)^{n_{\text{hard}}}\right) \left(1 - (\Theta^*)^{m_{\text{JC}}}\right) (1 - f_L)^2$$
#[must_use]
pub fn gold_viscoplastic_yield_stress_pa(
    temp_k: f64,
    melt_fraction: f64,
    eq_plastic_strain: f64,
) -> f64 {
    if melt_fraction >= 0.99 {
        return 0.0;
    }
    let theta_star = ((temp_k - T_REF_K) / (GOLD_18K_CARD.t_solidus_k - T_REF_K)).clamp(0.0, 1.0);
    let sigma_y0 = 175.0e6; // 175 MPa work-hardened 18k gold case band
    let k_hard = 220.0e6;
    let n_hard = 0.38;
    let m_jc = 1.15;
    let strain_term = sigma_y0 + k_hard * eq_plastic_strain.max(0.0).powf(n_hard);
    let thermal_softening = (1.0 - theta_star.powf(m_jc)).max(0.0);
    let solid_skeleton = (1.0 - melt_fraction.clamp(0.0, 1.0)).powi(2);
    strain_term * thermal_softening * solid_skeleton
}

/// Norton-Bailey Arrhenius creep + Perzyna viscoplastic overstress strain rate $\dot{\bar{\varepsilon}}^{vp}$ ($\text{s}^{-1}$)
/// of 18k Gold as a function of von Mises equivalent stress $\sigma_{\text{eq}}$ ($\text{Pa}$), temperature $T$ ($\text{K}$),
/// liquid fraction $f_L$, and accumulated plastic strain $\bar{\varepsilon}^p$.
#[must_use]
pub fn gold_viscoplastic_strain_rate_s_inv(
    sigma_eq_pa: f64,
    temp_k: f64,
    melt_fraction: f64,
    eq_plastic_strain: f64,
) -> f64 {
    if sigma_eq_pa <= 1.0 {
        return 0.0;
    }
    let sigma_y = gold_viscoplastic_yield_stress_pa(temp_k, melt_fraction, eq_plastic_strain);
    let theta_star = ((temp_k - T_REF_K) / (GOLD_18K_CARD.t_solidus_k - T_REF_K)).clamp(0.0, 1.0);

    // 1. Arrhenius dislocation-climb power-law creep: A * exp(-Q / RT) * (sigma / sigma_ref)^N
    let t_eff = temp_k.max(250.0);
    let arrhenius_creep = if t_eff > 500.0 {
        1.8e4 * (-14_800.0 / t_eff).exp() * (sigma_eq_pa / 25.0e6).powf(3.6)
    } else {
        0.0
    };

    // 2. Perzyna viscoplastic overstress rate above the temperature-softened yield surface
    let overstress_pa = (sigma_eq_pa - sigma_y).max(0.0);
    let fluidity = 0.35 + 18.0 * theta_star * theta_star + 65.0 * melt_fraction.clamp(0.0, 1.0);
    let perzyna_rate = fluidity * (overstress_pa / 12.0e6).powf(1.85);

    (arrhenius_creep + perzyna_rate).clamp(0.0, 50.0)
}

/// Temperature-dependent ductile fracture / necking elongation $\varepsilon_f(T)$ of 18k Gold,
/// capturing FCC gold's extreme high-temperature ductility prior to solidus.
#[must_use]
pub fn gold_ductility_elongation(temp_k: f64) -> f64 {
    let theta_star = ((temp_k - T_REF_K) / (GOLD_18K_CARD.t_solidus_k - T_REF_K)).clamp(0.0, 1.0);
    0.42 * (1.0 + 3.4 * theta_star * theta_star)
}

/// Unpolarized complex conductor Fresnel reflectance $F(\cos\theta_i, n, \kappa)$
/// delegated directly to `fs_render::conductor::fresnel_conductor`.
#[must_use]
pub fn conductor_fresnel_unpolarized(cos_theta_i: f64, eta: f64, kappa: f64) -> f64 {
    fresnel_conductor(
        1.0,
        eta.max(1.0e-3),
        kappa.max(1.0e-3),
        cos_theta_i.clamp(0.0, 1.0),
    )
    .unwrap_or(0.5)
}

/// Construct a validated 9-knot `fs_render::conductor::ConductorOptics` table across
/// `380..780 nm` for any horological material card, with an optional coherent dielectric
/// oxide coating (`fs_render::conductor::CoherentDielectricCoating`).
#[must_use]
pub fn material_conductor_optics(card: &MaterialCard, oxide_thickness_nm: f64) -> ConductorOptics {
    const WAVELENGTHS_NM: [f64; CONDUCTOR_IOR_SAMPLE_COUNT] = [
        380.0, 430.0, 480.0, 530.0, 580.0, 630.0, 680.0, 730.0, 780.0,
    ];
    let [nr, ng, nb] = card.optical_n_rgb;
    let [kr, kg, kb] = card.optical_k_rgb;

    let samples = WAVELENGTHS_NM.map(|wl| {
        // Interpolate piecewise across blue (465 nm), green (532 nm), and red (630 nm) anchors
        let (eta, k) = if wl <= 532.0 {
            let t = ((wl - 465.0) / (532.0 - 465.0)).clamp(0.0, 1.0);
            (nb + t * (ng - nb), kb + t * (kg - kb))
        } else {
            let t = ((wl - 532.0) / (630.0 - 532.0)).clamp(0.0, 1.5);
            (
                (ng + t * (nr - ng)).max(0.05),
                (kg + t * (kr - kg)).max(0.05),
            )
        };
        ConductorIorSample::try_new(wl, eta.max(0.05), k.max(0.05))
            .expect("valid horological conductor optical knot")
    });

    let mut hasher = DomainHasher::new("fs-pocket-watch.conductor-source.v1");
    hasher.update(card.matdb_pack_id.as_bytes());
    let source =
        ConductorSource::try_new(hasher.finalize(), ConductorDataStatus::MaterialStateClaims)
            .expect("valid non-zero Blake3 source hash");
    let optics = ConductorOptics::try_new(samples, source).expect("valid 9-knot conductor table");

    if oxide_thickness_nm > 0.5 {
        let mut coat_hasher = DomainHasher::new("fs-pocket-watch.oxide-coating.v1");
        coat_hasher.update(card.matdb_pack_id.as_bytes());
        coat_hasher.update(&oxide_thickness_nm.to_bits().to_le_bytes());
        let coat_source = ConductorSource::try_new(
            coat_hasher.finalize(),
            ConductorDataStatus::MaterialStateClaims,
        )
        .expect("valid coating source");
        // Magnetite / temper oxide phase index n(lambda) via Cauchy dispersion around n = 2.35
        let cauchy = CauchyIor::try_new(2.28, 0.018, 0.0).expect("valid Fe3O4 Cauchy dispersion");
        let thickness_m = (oxide_thickness_nm * 1.0e-9).clamp(0.0, 1.0e-4);
        if let Ok(coating) = CoherentDielectricCoating::try_new(cauchy, thickness_m, coat_source) {
            return optics.with_coating(coating);
        }
    }
    optics
}

/// Build a coupled simultaneous secondary Norton creep + isotropic $J_2$ plasticity law
/// (`fs_material::creep::plastic::NortonJ2`) for 18k Gold at temperature $T$ ($\text{K}$),
/// liquid fraction $f_L$, and time step $\Delta t$ ($\text{s}$).
#[must_use]
pub fn gold_coupled_norton_j2(temp_k: f64, melt_fraction: f64, dt_s: f64) -> Option<NortonJ2> {
    let t_clamped = temp_k.clamp(70.0, 1800.0);
    let e_pa = (GOLD_18K_CARD.youngs_modulus_pa * GOLD_18K_CARD.youngs_modulus_ratio(t_clamped))
        .max(2.5e8 * (1.0 - 0.92 * melt_fraction.clamp(0.0, 1.0)));
    let nu = 0.42;
    let elastic = IsotropicElastic::new(e_pa, nu, 0.05).ok()?;
    let yield_pa = gold_viscoplastic_yield_stress_pa(t_clamped, melt_fraction, 0.0).max(5.0e4);
    let hardening_pa =
        (1.15e8 * (1.0 - 0.85 * ((t_clamped - T_REF_K) / 1000.0).clamp(0.0, 1.0))).max(1.0e5);
    let plastic = J2Plasticity::new(elastic.clone(), yield_pa, hardening_pa).ok()?;

    let creep_params = NortonCreepParameters {
        reference_stress_pa: 25.0e6,
        reference_rate_s: 1.2e-3 * (1.0 + 35.0 * melt_fraction.clamp(0.0, 1.0)),
        stress_exponent: 3.6,
        activation_energy_j_mol: 123_000.0,
        reference_temperature_k: 950.0,
        temperature_range_k: [50.0, 2500.0],
    };
    let creep = NortonCreep::new(
        elastic,
        creep_params,
        t_clamped,
        dt_s.clamp(1.0e-6, 1.0),
    )
    .ok()?;
    NortonJ2::new(creep, plastic).ok()
}

/// Coherent thin-film Byrnes characteristic-matrix reflectance for thermally blued steel
/// evaluated directly through `fs_render::conductor::ConductorOptics::with_coating`
/// (`CoherentDielectricCoating` of $\text{Fe}_3\text{O}_4$ over `HIGH_CARBON_STEEL_CARD`).
///
/// - $d = 0\text{ nm}$: bright polished steel $[0.63, 0.61, 0.59]$
/// - $d \approx 46\text{ nm}$ ($500\text{ K}$ temper): straw gold (cancels $435\text{ nm}$ blue)
/// - $d \approx 60\text{ nm}$ ($545\text{ K}$ temper): violet-purple (cancels $560\text{ nm}$ green)
/// - $d \approx 72\text{ nm}$ (or 2nd-order $210\text{ nm}$, $575\text{ K}$ Patek temper): rich cornflower horological blue!
#[must_use]
pub fn blued_steel_interference_rgb(oxide_thickness_nm: f64, cos_theta: f64) -> [f64; 3] {
    let optics = material_conductor_optics(&HIGH_CARBON_STEEL_CARD, oxide_thickness_nm);
    let c = cos_theta.clamp(0.05, 1.0);
    let lambdas_nm = [630.0, 532.0, 465.0];
    let mut rgb = [0.0; 3];
    for ch in 0..3 {
        rgb[ch] = optics
            .fresnel(lambdas_nm[ch], 1.0, c)
            .unwrap_or(0.55)
            .clamp(0.04, 0.98);
    }
    rgb
}

/// Planck blackbody spectral radiance $B_\lambda(\lambda, T)$ mapped to linear sRGB
/// incandescence emission (`fs-render::thermal_surface::planck_spectral_radiance_nm`).
///
/// Zero below $680\text{ K}$; transitions through Draper dull crimson ($800\text{ K}$),
/// cherry red ($1050\text{ K}$), bright molten gold orange ($1337\text{ K}$), to white-hot ($2000\text{ K}$).
#[must_use]
pub fn planck_incandescence_linear_rgb(temp_k: f64, emissivity: f64) -> [f64; 3] {
    if temp_k < 680.0 {
        return [0.0, 0.0, 0.0];
    }
    // Second radiation constant c2 = h*c/k_B = 1.438_776_877e7 nm*K
    let c2_nm_k = 1.438_776_877e7;
    let lambdas_nm = [630.0, 532.0, 465.0];
    // Normalize relative to 1450 K reference studio exposure so molten gold (1337 K) glows vividly
    let t_norm = 1450.0;
    let mut out = [0.0; 3];
    for ch in 0..3 {
        let lam = lambdas_nm[ch];
        let x = (c2_nm_k / (lam * temp_k)).clamp(0.5, 80.0);
        let x_ref = c2_nm_k / (lam * t_norm);
        let b_rel = ((x_ref.exp() - 1.0) / (x.exp() - 1.0)).clamp(0.0, 25.0);
        out[ch] = (b_rel * (0.45 + 0.55 * emissivity)).clamp(0.0, 8.0);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gold_enthalpy_melting_roundtrip_and_plateau() {
        let h_room = GOLD_18K_CARD.enthalpy_from_temperature_j_kg(T_REF_K);
        assert!(h_room.abs() < 1.0e-9);

        // Halfway through the solidus-liquidus mushy interval
        let t_mid = 0.5 * (GOLD_18K_CARD.t_solidus_k + GOLD_18K_CARD.t_liquidus_k);
        let h_mid = GOLD_18K_CARD.enthalpy_from_temperature_j_kg(t_mid);
        let (t_rec, phi_mid) = GOLD_18K_CARD.temperature_and_melt_from_enthalpy(h_mid);
        assert!((t_rec - t_mid).abs() < 1.0e-6);
        assert!((phi_mid - 0.5).abs() < 0.05);

        // Above liquidus (1400 K), melt fraction is 1.0
        let h_liq = GOLD_18K_CARD.enthalpy_from_temperature_j_kg(1400.0);
        let (t_liq, phi_liq) = GOLD_18K_CARD.temperature_and_melt_from_enthalpy(h_liq);
        assert!((t_liq - 1400.0).abs() < 1.0e-6);
        assert!((phi_liq - 1.0).abs() < 1.0e-9);
    }

    #[test]
    fn watch_oil_vft_viscosity_freezes_at_ln2_and_thins_when_heated() {
        let mu_room = watch_oil_viscosity_pa_s(T_REF_K);
        let mu_hot = watch_oil_viscosity_pa_s(373.15);
        let mu_ln2 = watch_oil_viscosity_pa_s(T_LN2_SAT_K);
        assert!(mu_room > 0.05 && mu_room < 0.15);
        assert!(mu_hot < 0.25 * mu_room);
        assert!(mu_ln2 > 1.0e4 * mu_room);
    }

    #[test]
    fn blued_steel_thin_film_produces_cornflower_blue_at_72nm() {
        let rgb = blued_steel_interference_rgb(72.0, 1.0);
        // At 72 nm Fe3O4 (300 C temper), destructive interference of 630 nm red makes blue > red
        assert!(rgb[2] > rgb[0], "expected blue > red at 72nm, got {rgb:?}");
    }
}
