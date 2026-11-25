//! Sourced material library for the Siemens SST-600 steam turbine digital twin.
//!
//! Every mechanical, thermal, tribological, and electromagnetic constant is
//! loaded from FrankenSim's `data/matdb/seed-v1` packs embedded via [`crate::matdb_src`],
//! with thermal blackbody chromaticity computed through `fs_render::thermal_surface`
//! and high-temperature rotor/casing creep evaluated via `fs_material::creep`.

use std::sync::OnceLock;

use fs_blake3::DomainHasher;
use fs_material::creep::{NortonCreep, NortonCreepParameters};
use fs_material::elastic::IsotropicElastic;
use fs_render::thermal_surface::{DiffuseThermalSurface, GrayDiffuseThermalSurface};

use crate::matdb_src::{self as db, Curve};

/// Standard reference temperature [K] (20 °C).
pub const T_REF_K: f64 = 293.15;
/// Standard atmospheric pressure [Pa].
pub const P_ATM_PA: f64 = 101_325.0;
/// Stefan–Boltzmann constant [W/(m²·K⁴)].
pub const STEFAN_BOLTZMANN: f64 = 5.670_374_419e-8;
/// Standard gravity [m/s²].
pub const G0_M_S2: f64 = 9.806_65;

/// Sourced turbine solid alloy with optional temperature-dependent matdb curves.
#[derive(Debug, Clone)]
pub struct TurbineSolid {
    pub name: &'static str,
    pub matdb_sources: &'static str,
    pub rho_kg_m3: f64,
    pub cp_j_kg_k: f64,
    pub k_w_m_k: f64,
    pub young_modulus_pa: f64,
    pub poisson_ratio: f64,
    pub alpha_per_k: f64,
    pub yield_strength_pa: f64,
    pub emissivity: f64,
    e_curve: Option<Curve>,
    k_curve: Option<Curve>,
    cp_curve: Option<Curve>,
    alpha_curve: Option<Curve>,
}

impl TurbineSolid {
    /// Temperature-dependent Young's modulus [Pa].
    #[must_use]
    pub fn young_modulus_at(&self, t_k: f64) -> f64 {
        if let Some(ref c) = self.e_curve {
            // Scale curve to this alloy's room-temperature anchor if curve anchor differs.
            let ref_val = c.eval(T_REF_K).max(1.0);
            self.young_modulus_pa * (c.eval(t_k) / ref_val)
        } else {
            // Mild ferritic steel thermal softening (-3.0e-4 / K above 293 K)
            let dt = (t_k - T_REF_K).max(0.0);
            self.young_modulus_pa * (1.0 - 3.0e-4 * dt).clamp(0.55, 1.0)
        }
    }

    /// Temperature-dependent thermal conductivity [W/(m·K)].
    #[must_use]
    pub fn conductivity_at(&self, t_k: f64) -> f64 {
        if let Some(ref c) = self.k_curve {
            c.eval(t_k)
        } else {
            self.k_w_m_k
        }
    }

    /// Temperature-dependent specific heat capacity [J/(kg·K)].
    #[must_use]
    pub fn cp_at(&self, t_k: f64) -> f64 {
        if let Some(ref c) = self.cp_curve {
            c.eval(t_k)
        } else {
            self.cp_j_kg_k
        }
    }

    /// Mean linear thermal expansion coefficient between `T_REF_K` and `t_k` [1/K].
    #[must_use]
    pub fn mean_alpha_at(&self, t_k: f64) -> f64 {
        if let Some(ref c) = self.alpha_curve {
            c.eval(t_k)
        } else {
            self.alpha_per_k
        }
    }

    /// Linear thermal strain `ε_th = α(T) * (T - T_ref)`.
    #[must_use]
    pub fn thermal_strain(&self, t_k: f64) -> f64 {
        self.mean_alpha_at(t_k) * (t_k - T_REF_K)
    }

    /// Thermal diffusivity `k / (ρ c_p)` [m²/s].
    #[must_use]
    pub fn diffusivity_at(&self, t_k: f64) -> f64 {
        self.conductivity_at(t_k) / (self.rho_kg_m3 * self.cp_at(t_k))
    }

    /// Construct a FrankenSim `IsotropicElastic` law at temperature `t_k`.
    #[must_use]
    pub fn isotropic_elastic_at(&self, t_k: f64) -> IsotropicElastic {
        IsotropicElastic::new(self.young_modulus_at(t_k), self.poisson_ratio, 0.02)
            .expect("valid turbine solid elastic constants")
    }
}

/// Bearing lubricating mineral oil from `shell-heat-transfer-oil-s2-2011` (ISO VG 32 class).
///
/// Uses log-linear viscosity interpolation between the 0 / 40 / 100 / 200 °C matdb
/// knots so Reynolds film viscosity stays smooth and strictly positive across the full
/// turbine bearing temperature range (20–160 °C).
#[derive(Debug, Clone)]
pub struct TurbineLubeOil {
    rho_curve: Curve,
    cp_curve: Curve,
    k_curve: Curve,
    nu_curve: Curve,
    /// ISO VG 32 scaling ratio relative to the Shell S2 40 °C knot (`25 cSt -> 32 cSt`).
    vg32_scale: f64,
}

impl TurbineLubeOil {
    fn build() -> Self {
        let rho_curve = db::SHELL_OIL_S2.require_curve("rho-0-200c");
        let cp_curve = db::SHELL_OIL_S2.require_curve("cp-0-200c");
        let k_curve = db::SHELL_OIL_S2.require_curve("k-0-200c");
        let nu_curve = db::SHELL_OIL_S2.require_curve("nu-0-200c");
        // Shell Heat Transfer Oil S2 has nu(40 C) = 25 mm^2/s; standard Siemens SST-600
        // turbine oil is ISO VG 32 (32 mm^2/s at 40 C), scaled cleanly on the same mineral base curve.
        let vg32_scale = 32.0e-6 / nu_curve.eval(313.15);
        Self {
            rho_curve,
            cp_curve,
            k_curve,
            nu_curve,
            vg32_scale,
        }
    }

    /// Oil density [kg/m³].
    #[must_use]
    pub fn density(&self, t_k: f64) -> f64 {
        self.rho_curve.eval(t_k.clamp(273.15, 473.15))
    }

    /// Oil specific heat [J/(kg·K)].
    #[must_use]
    pub fn cp(&self, t_k: f64) -> f64 {
        self.cp_curve.eval(t_k.clamp(273.15, 473.15))
    }

    /// Oil thermal conductivity [W/(m·K)].
    #[must_use]
    pub fn conductivity(&self, t_k: f64) -> f64 {
        self.k_curve.eval(t_k.clamp(273.15, 473.15))
    }

    /// Kinematic viscosity [m²/s] via piecewise log-linear interpolation in `1/T`
    /// across the 0 / 40 / 100 / 200 °C matdb knots (exact at every knot).
    #[must_use]
    pub fn kinematic_viscosity(&self, t_k: f64) -> f64 {
        let t = t_k.clamp(273.15, 473.15);
        let knots = &self.nu_curve.knots;
        if t <= knots[0].0 {
            return knots[0].1 * self.vg32_scale;
        }
        let last = knots.len() - 1;
        if t >= knots[last].0 {
            return knots[last].1 * self.vg32_scale;
        }
        let i = knots
            .partition_point(|&(x, _)| x < t)
            .saturating_sub(1)
            .min(last - 1);
        let (t0, nu0) = knots[i];
        let (t1, nu1) = knots[i + 1];
        let f = (1.0 / t - 1.0 / t0) / (1.0 / t1 - 1.0 / t0);
        let ln_nu = nu0.ln() * (1.0 - f) + nu1.ln() * f;
        ln_nu.exp() * self.vg32_scale
    }

    /// Dynamic viscosity `μ = ρ ν` [Pa·s].
    #[must_use]
    pub fn dynamic_viscosity(&self, t_k: f64) -> f64 {
        self.density(t_k) * self.kinematic_viscosity(t_k)
    }
}

/// Condenser cooling water curves from `water-liquid-iapws-sr6-08`.
#[derive(Debug, Clone)]
pub struct CoolingWaterCurves {
    rho: Curve,
    cp: Curve,
    k: Curve,
    mu: Curve,
}

impl CoolingWaterCurves {
    fn build() -> Self {
        Self {
            rho: db::WATER_LIQ.require_curve("density-283.15-353.15k"),
            cp: db::WATER_LIQ.require_curve("specific-heat-capacity-283.15-353.15k"),
            k: db::WATER_LIQ.require_curve("thermal-conductivity-283.15-353.15k"),
            mu: db::WATER_LIQ.require_curve("dynamic-viscosity-283.15-353.15k"),
        }
    }

    #[must_use]
    pub fn density(&self, t_k: f64) -> f64 {
        self.rho.eval(t_k)
    }
    #[must_use]
    pub fn cp(&self, t_k: f64) -> f64 {
        self.cp.eval(t_k)
    }
    #[must_use]
    pub fn conductivity(&self, t_k: f64) -> f64 {
        self.k.eval(t_k)
    }
    #[must_use]
    pub fn dynamic_viscosity(&self, t_k: f64) -> f64 {
        self.mu.eval(t_k)
    }
}

/// Synchronous generator electromagnetic materials from `nasa-tn-d-8184-m19-material-deck`,
/// `torrent-2018-m19-steinmetz-inputs`, and `copper-annealed-iacs-nbs-hb100`.
#[derive(Debug, Clone)]
pub struct GeneratorMaterials {
    pub m19_bh_curve: Curve,
    pub m19_specific_core_loss_ref_w_kg: f64,
    pub m19_freq_exponent: f64,
    pub steinmetz_kh: f64,
    pub steinmetz_a: f64,
    pub steinmetz_n: f64,
    pub copper_rho_20c_ohm_m: f64,
    pub copper_rho_slope_ohm_m_per_k: f64,
    pub copper_density_kg_m3: f64,
    pub copper_melting_point_k: f64,
}

impl GeneratorMaterials {
    fn build() -> Self {
        let cu_curve = db::COPPER_IACS.require_curve("rho-10-30c");
        let rho_10 = cu_curve.eval(283.15);
        let rho_20 = cu_curve.eval(293.15);
        let rho_30 = cu_curve.eval(303.15);
        let slope = (rho_30 - rho_10) / 20.0;
        Self {
            m19_bh_curve: db::M19_STEEL.require_curve("m19-magnetization-curve"),
            m19_specific_core_loss_ref_w_kg: db::M19_STEEL
                .require("m19-specific-core-loss-reference"),
            m19_freq_exponent: db::M19_STEEL.require("m19-core-loss-frequency-exponent"),
            steinmetz_kh: db::M19_STEINMETZ.require("eq4-kh"),
            steinmetz_a: db::M19_STEINMETZ.require("eq4-frequency-exponent"),
            steinmetz_n: db::M19_STEINMETZ.require("eq4-flux-exponent"),
            copper_rho_20c_ohm_m: rho_20,
            copper_rho_slope_ohm_m_per_k: slope,
            copper_density_kg_m3: db::COPPER_PURE.require("density"),
            copper_melting_point_k: db::COPPER_PURE.require("melting-point"),
        }
    }

    /// Copper winding resistivity at temperature `t_k` [Ω·m].
    #[must_use]
    pub fn copper_resistivity(&self, t_k: f64) -> f64 {
        (self.copper_rho_20c_ohm_m + self.copper_rho_slope_ohm_m_per_k * (t_k - T_REF_K))
            .max(1.0e-9)
    }

    /// M19 stator iron core loss per kg at electrical frequency `f_hz` and peak flux `b_peak_t` [W/kg].
    #[must_use]
    pub fn m19_specific_core_loss(&self, f_hz: f64, b_peak_t: f64) -> f64 {
        if f_hz <= 0.0 || b_peak_t <= 0.0 {
            return 0.0;
        }
        // Torrent (2018) Eq. 4 hysteresis + classical eddy contribution scaled to 50 Hz M19 core.
        let p_hyst = 1.0e-2
            * self.steinmetz_kh
            * f_hz.powf(self.steinmetz_a)
            * b_peak_t.powf(self.steinmetz_n);
        let p_ref_scale = self.m19_specific_core_loss_ref_w_kg
            * (f_hz / 400.0).powf(self.m19_freq_exponent)
            * (b_peak_t / 1.0).powi(2);
        0.5 * (p_hyst + p_ref_scale)
    }
}

/// Dry ambient powerhouse air parameters from `air-dry-ussa1976`.
#[derive(Debug, Clone, Copy)]
pub struct AmbientAir {
    pub molar_mass_kg_mol: f64,
    pub gamma: f64,
    pub mu_ref_pa_s: f64,
    pub t_ref_k: f64,
    pub sutherland_s_k: f64,
}

impl AmbientAir {
    fn build() -> Self {
        Self {
            molar_mass_kg_mol: db::DRY_AIR.require("molar-mass"),
            gamma: db::DRY_AIR.require("gamma"),
            mu_ref_pa_s: db::DRY_AIR.require("mu-reference"),
            t_ref_k: db::DRY_AIR.require("temperature-reference"),
            sutherland_s_k: db::DRY_AIR.require("sutherland-s"),
        }
    }

    /// Sutherland dynamic viscosity [Pa·s].
    #[must_use]
    pub fn dynamic_viscosity(&self, t_k: f64) -> f64 {
        let t = t_k.max(150.0);
        self.mu_ref_pa_s
            * (t / self.t_ref_k).powf(1.5)
            * ((self.t_ref_k + self.sutherland_s_k) / (t + self.sutherland_s_k))
    }

    /// Ideal-gas density [kg/m³] at `(t_k, p_pa)`.
    #[must_use]
    pub fn density(&self, t_k: f64, p_pa: f64) -> f64 {
        let r_spec = fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / self.molar_mass_kg_mol;
        p_pa / (r_spec * t_k.max(150.0))
    }

    /// Acoustic speed of sound [m/s].
    #[must_use]
    pub fn sound_speed(&self, t_k: f64) -> f64 {
        let r_spec = fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / self.molar_mass_kg_mol;
        (self.gamma * r_spec * t_k.max(150.0)).sqrt()
    }
}

/// Physically computed CIE 1931 → linear sRGB blackbody incandescence LUT
/// built via `fs_render::thermal_surface::GrayDiffuseThermalSurface`.
#[derive(Debug, Clone)]
pub struct BlackbodyLut {
    pub t0_k: f64,
    pub dt_k: f64,
    pub rgb: Vec<[f64; 3]>,
    pub luminance_cd_m2: Vec<f64>,
}

impl BlackbodyLut {
    fn build() -> Self {
        let t0_k = 500.0;
        let dt_k = 20.0;
        let n = 66; // 500 K .. 1800 K
        let mut hasher = DomainHasher::new("fs-steam-turbine.blackbody-lut.v1");
        hasher.update(b"gray-diffuse-turbine-metal");
        let id = hasher.finalize();
        let mut rgb = Vec::with_capacity(n);
        let mut lum = Vec::with_capacity(n);
        for i in 0..n {
            let t = t0_k + dt_k * i as f64;
            let obs = GrayDiffuseThermalSurface::try_new(t, 1.0, id)
                .ok()
                .map(DiffuseThermalSurface::from)
                .and_then(|s| s.visible_observables().ok());
            match obs {
                Some(o) => {
                    let c = fs_render::spectral::xyz_to_linear_srgb(o.xyz());
                    let m = c[0].max(c[1]).max(c[2]).max(1.0e-30);
                    rgb.push([
                        (c[0] / m).max(0.0),
                        (c[1] / m).max(0.0),
                        (c[2] / m).max(0.0),
                    ]);
                    lum.push(o.luminance_cd_m2());
                }
                None => {
                    rgb.push([1.0, 0.25, 0.04]);
                    lum.push(0.0);
                }
            }
        }
        Self {
            t0_k,
            dt_k,
            rgb,
            luminance_cd_m2: lum,
        }
    }

    /// Sample normalized chromaticity and physical luminance [cd/m²] at `t_k` with emissivity `eps`.
    #[must_use]
    pub fn sample(&self, t_k: f64, eps: f64) -> ([f64; 3], f64) {
        if t_k < self.t0_k {
            return ([0.0; 3], 0.0);
        }
        let x = ((t_k - self.t0_k) / self.dt_k).min((self.rgb.len() - 1) as f64);
        let i = (x.floor() as usize).min(self.rgb.len() - 2);
        let f = x - i as f64;
        let mut c = [0.0; 3];
        for ch in 0..3 {
            c[ch] = self.rgb[i][ch] * (1.0 - f) + self.rgb[i + 1][ch] * f;
        }
        let l = self.luminance_cd_m2[i] * (1.0 - f) + self.luminance_cd_m2[i + 1] * f;
        (c, l * eps)
    }

    /// Display-referred blackbody glow for WebGL shaders (0 below Draper point ~770 K).
    #[must_use]
    pub fn display_rgb(&self, t_k: f64, eps: f64) -> [f64; 3] {
        let (c, l) = self.sample(t_k, eps);
        let exposure = ((l / 1.0e-3).max(1.0).ln() / (1.0e5_f64).ln()).clamp(0.0, 2.0);
        [c[0] * exposure, c[1] * exposure, c[2] * exposure]
    }
}

/// Complete material library for the Siemens SST-600 steam turbine digital twin.
pub struct TurbineMaterialLibrary {
    /// Forged CrMoV rotor shaft (AISI 4140/4130 + 1018 thermal + 1.4016 high-T E/alpha curves).
    pub rotor_steel: TurbineSolid,
    /// Cast ferritic CrMoV HP/IP split casing (`stainless-14016` + `steel-4130`).
    pub hp_casing_steel: TurbineSolid,
    /// Fabricated structural steel axial exhaust hood (`steel-a36` + `steel-1018`).
    pub lp_exhaust_steel: TurbineSolid,
    /// 12% Cr ferritic-martensitic HP/IP blading (`stainless-14016` + `17-4PH` yield).
    pub hp_blade_steel: TurbineSolid,
    /// Precipitation-hardened martensitic LP last-stage blading (`stainless-17-4ph-h1025`).
    pub lp_blade_steel: TurbineSolid,
    /// High-temperature labyrinth seal fin strips & casing bolts (`inconel-718` + `inconel-625`).
    pub seal_inconel: TurbineSolid,
    /// Bearing shell bronze backing (`bronze-c93200`) + tin-base Babbitt melting limit (`tin-pure`).
    pub bearing_bronze: TurbineSolid,
    /// Tin-base Babbitt melting point [K] (`tin-pure-nbs-c447`: 505.05 K).
    pub babbitt_melt_k: f64,
    /// Bearing oil-film / Babbitt trip temperature [K] (120 °C = 393.15 K, industrial alarm 110 °C).
    pub babbitt_wipe_onset_k: f64,
    /// Baseplate / pedestal cast iron (`gray-cast-iron-s2-s`).
    pub pedestal_cast_iron: TurbineSolid,
    /// Turbine lube oil (ISO VG 32 anchored to `shell-heat-transfer-oil-s2-2011`).
    pub lube_oil: TurbineLubeOil,
    /// Condenser cooling water (`water-liquid-iapws-sr6-08`).
    pub cooling_water: CoolingWaterCurves,
    /// Synchronous generator M19 silicon steel & IACS annealed copper.
    pub generator: GeneratorMaterials,
    /// Ambient turbine-hall air (`air-dry-ussa1976`).
    pub ambient_air: AmbientAir,
    /// Physically computed blackbody incandescence lookup table (`fs_render`).
    pub blackbody: BlackbodyLut,
    /// CrMoV high-temperature secondary creep parameters (`fs_material::creep`).
    pub rotor_creep_params: NortonCreepParameters,
}

impl TurbineMaterialLibrary {
    fn build() -> Self {
        let ss14016_e = db::SS_14016.require_curve("E");
        let ss14016_alpha = db::SS_14016.require_curve("mean-alpha");
        let inc718_k = db::INCONEL_718.require_curve("conductivity");
        let inc625_cp = db::INCONEL_625.require_curve("heat-capacity");

        let rotor_steel = TurbineSolid {
            name: "CrMoV Forged Rotor Steel (AISI 4140/4130)",
            matdb_sources: "aisi-4140-rc33, steel-4130-sheet-normalized-mil-hdbk-5j, steel-1018-onlinemetals-reference, stainless-14016-thyssenkrupp-reference",
            rho_kg_m3: db::STEEL_4130.require("density"),
            cp_j_kg_k: db::STEEL_1018.require("heat-capacity"),
            k_w_m_k: db::STEEL_1018.require("conductivity"),
            young_modulus_pa: db::STEEL_4130.require("young-modulus"),
            poisson_ratio: db::STEEL_4130.require("poisson-ratio"),
            alpha_per_k: ss14016_alpha.eval(573.15),
            yield_strength_pa: db::AISI_4140.require("yield-strength-0p2-offset-26p7C"),
            emissivity: 0.35,
            e_curve: Some(ss14016_e.clone()),
            k_curve: None,
            cp_curve: None,
            alpha_curve: Some(ss14016_alpha.clone()),
        };

        let hp_casing_steel = TurbineSolid {
            name: "Cast CrMoV HP/IP Casing (1.4016 / 4130 class)",
            matdb_sources: "stainless-14016-thyssenkrupp-reference, steel-4130-sheet-normalized-mil-hdbk-5j",
            rho_kg_m3: db::SS_14016.require("rho"),
            cp_j_kg_k: db::SS_14016.require("cp"),
            k_w_m_k: db::SS_14016.require("k"),
            young_modulus_pa: ss14016_e.eval(T_REF_K),
            poisson_ratio: db::STEEL_4130.require("poisson-ratio"),
            alpha_per_k: ss14016_alpha.eval(573.15),
            yield_strength_pa: db::STEEL_4130.require("tensile-yield"),
            emissivity: 0.65,
            e_curve: Some(ss14016_e.clone()),
            k_curve: None,
            cp_curve: None,
            alpha_curve: Some(ss14016_alpha.clone()),
        };

        let lp_exhaust_steel = TurbineSolid {
            name: "Fabricated Structural Steel Exhaust Hood (ASTM A36)",
            matdb_sources: "steel-a36-onlinemetals-reference, steel-1018-onlinemetals-reference",
            rho_kg_m3: db::STEEL_A36.require("density"),
            cp_j_kg_k: db::STEEL_A36.require("heat-capacity"),
            k_w_m_k: db::STEEL_A36.require("conductivity"),
            young_modulus_pa: db::STEEL_A36.require("young-modulus"),
            poisson_ratio: db::STEEL_1018.require("poisson-ratio"),
            alpha_per_k: 12.0e-6,
            yield_strength_pa: db::STEEL_A36.require("yield-strength"),
            emissivity: 0.78, // Industrial slate-blue coated casing
            e_curve: None,
            k_curve: None,
            cp_curve: None,
            alpha_curve: None,
        };

        let hp_blade_steel = TurbineSolid {
            name: "12% Cr Martensitic HP/IP Blade Steel (1.4016/17-4PH)",
            matdb_sources: "stainless-14016-thyssenkrupp-reference, stainless-17-4ph-h1025-bar-mil-hdbk-5j",
            rho_kg_m3: db::SS_14016.require("rho"),
            cp_j_kg_k: db::SS_14016.require("cp"),
            k_w_m_k: db::SS_14016.require("k"),
            young_modulus_pa: db::SS_17_4PH.require("young-modulus"),
            poisson_ratio: db::SS_17_4PH.require("poisson-ratio"),
            alpha_per_k: ss14016_alpha.eval(573.15),
            yield_strength_pa: db::SS_17_4PH.require("tensile-yield-l"),
            emissivity: 0.28,
            e_curve: Some(ss14016_e.clone()),
            k_curve: None,
            cp_curve: None,
            alpha_curve: Some(ss14016_alpha.clone()),
        };

        let lp_blade_steel = TurbineSolid {
            name: "17-4PH H1025 Precipitation-Hardened LP Last-Stage Blade",
            matdb_sources: "stainless-17-4ph-h1025-bar-mil-hdbk-5j, stainless-14016-thyssenkrupp-reference",
            rho_kg_m3: db::SS_14016.require("rho"),
            cp_j_kg_k: db::SS_14016.require("cp"),
            k_w_m_k: db::SS_14016.require("k"),
            young_modulus_pa: db::SS_17_4PH.require("young-modulus"),
            poisson_ratio: db::SS_17_4PH.require("poisson-ratio"),
            alpha_per_k: ss14016_alpha.eval(373.15),
            yield_strength_pa: db::SS_17_4PH.require("tensile-yield-l"),
            emissivity: 0.25,
            e_curve: Some(ss14016_e),
            k_curve: None,
            cp_curve: None,
            alpha_curve: Some(ss14016_alpha),
        };

        let seal_inconel = TurbineSolid {
            name: "Inconel 718/625 Labyrinth Seal Strips & HP Bolts",
            matdb_sources: "inconel-718-aged-special-metals-2007, inconel-625-special-metals-2013",
            rho_kg_m3: db::INCONEL_718.require("density"),
            cp_j_kg_k: db::INCONEL_718.require("heat-capacity"),
            k_w_m_k: inc718_k.eval(T_REF_K),
            young_modulus_pa: db::INCONEL_718.require("dynamic-young-modulus"),
            poisson_ratio: 0.29,
            alpha_per_k: 13.0e-6,
            yield_strength_pa: 1.03e9,
            emissivity: 0.32,
            e_curve: None,
            k_curve: Some(inc718_k),
            cp_curve: Some(inc625_cp),
            alpha_curve: None,
        };

        let bearing_bronze = TurbineSolid {
            name: "C93200 Bearing Bronze Shell with Sn-Babbitt Lining",
            matdb_sources: "bronze-c93200-cda-reference, tin-pure-nbs-c447",
            rho_kg_m3: db::BRONZE_C93200.require("density"),
            cp_j_kg_k: db::BRONZE_C93200.require("heat-capacity"),
            k_w_m_k: db::BRONZE_C93200.require("conductivity"),
            young_modulus_pa: db::BRONZE_C93200.require("young-modulus"),
            poisson_ratio: 0.34,
            alpha_per_k: 18.0e-6,
            yield_strength_pa: 124.0e6,
            emissivity: 0.40,
            e_curve: None,
            k_curve: None,
            cp_curve: None,
            alpha_curve: None,
        };

        let babbitt_melt_k = db::TIN_PURE.require("melting-point");

        let pedestal_cast_iron = TurbineSolid {
            name: "Gray Cast Iron Bearing Pedestal & Baseplate",
            matdb_sources: "gray-cast-iron-s2-s, steel-a36-onlinemetals-reference",
            rho_kg_m3: 7_200.0,
            cp_j_kg_k: db::STEEL_A36.require("heat-capacity"),
            k_w_m_k: db::CAST_IRON.require("thermal-conductivity"),
            young_modulus_pa: 115.0e9,
            poisson_ratio: 0.26,
            alpha_per_k: 11.0e-6,
            yield_strength_pa: db::CAST_IRON.require("ultimate-tensile-strength"),
            emissivity: 0.85,
            e_curve: None,
            k_curve: None,
            cp_curve: None,
            alpha_curve: None,
        };

        // Secondary Norton-Arrhenius creep parameters representative of 1Cr-1Mo-0.25V
        // steam turbine rotor forgings (Viswanathan, Damage Mechanisms & Life Assessment).
        let rotor_creep_params = NortonCreepParameters {
            reference_stress_pa: 150.0e6,
            reference_rate_s: 1.0e-10,
            stress_exponent: 6.5,
            activation_energy_j_mol: 320_000.0,
            reference_temperature_k: 813.15, // 540 C
            temperature_range_k: [273.15, 1100.0],
        };

        Self {
            rotor_steel,
            hp_casing_steel,
            lp_exhaust_steel,
            hp_blade_steel,
            lp_blade_steel,
            seal_inconel,
            bearing_bronze,
            babbitt_melt_k,
            babbitt_wipe_onset_k: 403.15, // 130 C Babbitt softening / wipe onset (well below Sn melt 505 K)
            pedestal_cast_iron,
            lube_oil: TurbineLubeOil::build(),
            cooling_water: CoolingWaterCurves::build(),
            generator: GeneratorMaterials::build(),
            ambient_air: AmbientAir::build(),
            blackbody: BlackbodyLut::build(),
            rotor_creep_params,
        }
    }

    /// Evaluate instantaneous secondary Norton creep strain rate [1/s] using
    /// `fs_material::creep::NortonCreep` at von Mises stress `sigma_eq_pa` and temperature `t_k`.
    #[must_use]
    pub fn rotor_creep_rate_per_s(&self, sigma_eq_pa: f64, t_k: f64) -> f64 {
        let t = t_k.clamp(
            self.rotor_creep_params.temperature_range_k[0],
            self.rotor_creep_params.temperature_range_k[1],
        );
        if sigma_eq_pa <= 0.0 || t < 550.0 {
            return 0.0;
        }
        let elastic = self.rotor_steel.isotropic_elastic_at(t);
        if let Ok(law) = NortonCreep::new(elastic, self.rotor_creep_params.clone(), t, 1.0) {
            law.equivalent_creep_rate_s(sigma_eq_pa).unwrap_or(0.0)
        } else {
            0.0
        }
    }
}

static LIB: OnceLock<TurbineMaterialLibrary> = OnceLock::new();

/// Shared turbine material library initialized once from compile-time embedded matdb packs.
pub fn library() -> &'static TurbineMaterialLibrary {
    LIB.get_or_init(TurbineMaterialLibrary::build)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn all_21_matdb_packs_load_and_produce_physical_values() {
        let lib = library();

        // Rotor forged steel
        assert!((lib.rotor_steel.rho_kg_m3 - 7833.4).abs() < 20.0);
        assert!(lib.rotor_steel.yield_strength_pa > 900.0e6);
        assert!(lib.rotor_steel.young_modulus_at(773.15) < lib.rotor_steel.young_modulus_at(293.15));
        assert!(lib.rotor_steel.thermal_strain(798.15) > 0.005);

        // Blade and seal alloys
        assert!(lib.lp_blade_steel.yield_strength_pa > 950.0e6);
        assert!(lib.seal_inconel.conductivity_at(800.0) > lib.seal_inconel.conductivity_at(300.0));

        // Babbitt / tin melting point from NBS C447
        assert!((lib.babbitt_melt_k - 505.05).abs() < 1.0);

        // ISO VG 32 lube oil anchored to Shell S2 knots
        let nu_40c = lib.lube_oil.kinematic_viscosity(313.15);
        assert!((nu_40c - 32.0e-6).abs() < 1.0e-7, "nu_40c={nu_40c}");
        let mu_50c = lib.lube_oil.dynamic_viscosity(323.15);
        let mu_90c = lib.lube_oil.dynamic_viscosity(363.15);
        assert!(mu_50c > mu_90c && mu_90c > 0.003);

        // Cooling water IAPWS SR6-08
        assert!((lib.cooling_water.density(293.15) - 998.2).abs() < 0.5);
        assert!((lib.cooling_water.cp(293.15) - 4184.0).abs() < 5.0);

        // Generator M19 + IACS copper
        assert!(lib.generator.copper_resistivity(373.15) > lib.generator.copper_resistivity(293.15));
        assert!(lib.generator.m19_specific_core_loss(50.0, 1.4) > 0.5);

        // Ambient air USSA-1976
        assert!((lib.ambient_air.sound_speed(293.15) - 343.2).abs() < 1.5);

        // Blackbody LUT from fs_render
        let (_rgb_cold, lum_cold) = lib.blackbody.sample(400.0, 0.8);
        let (rgb_hot, lum_hot) = lib.blackbody.sample(1100.0, 0.8);
        assert_eq!(lum_cold, 0.0);
        assert!(lum_hot > 1.0 && rgb_hot[0] > rgb_hot[2]);

        // Norton creep rate rises sharply with temperature
        let rate_450c = lib.rotor_creep_rate_per_s(180.0e6, 723.15);
        let rate_550c = lib.rotor_creep_rate_per_s(180.0e6, 823.15);
        assert!(rate_550c > 100.0 * rate_450c && rate_550c > 0.0);
    }
}
