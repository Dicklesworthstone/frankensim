//! Material models for the FS-1903 watch, sourced from the matdb seed packs
//! (see [`crate::matdb_src`]) wherever a row exists, and from cited handbook
//! constants where the seed database has no row yet.
//!
//! What is derived rather than tabulated:
//! - 18k yellow gold (Au 750 / Ag 125 / Cu 125 by mass) from its three
//!   constituents: density by volume additivity, Debye-shaped heat capacity
//!   (Kopp–Neumann mass mixing, each element scaled to its matdb room-temperature
//!   value), latent heat by mass mixing, and thermal conductivity from
//!   Wiedemann–Franz with a Matthiessen residual resistivity. The enthalpy ↔
//!   temperature ↔ liquid-fraction map is an `fs_material` equilibrium phase curve.
//! - Watch-oil viscosity: a Vogel–Fulcher–Tammann law fitted exactly through the
//!   three Pennzane SHF X-2000 kinematic-viscosity rows (−40, 40, 100 °C).
//! - Incandescent colour: CIE 1931 tristimulus of a gray Planck emitter through
//!   `fs_render::thermal_surface`, tabulated once.
//! - Liquid-nitrogen saturation temperature versus chamber pressure: Clausius–Clapeyron.

use std::sync::OnceLock;

use fs_blake3::DomainHasher;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_render::conductor::{
    CONDUCTOR_IOR_SAMPLE_COUNT, CoherentDielectricCoating, ConductorDataStatus,
    ConductorIorSample, ConductorOptics, ConductorSource,
};
use fs_render::dielectric::CauchyIor;
use fs_render::thermal_surface::{DiffuseThermalSurface, GrayDiffuseThermalSurface};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R_GAS;

use crate::matdb_src::{self as db, Curve};

pub const T_REF_K: f64 = 293.15;
pub const STEFAN_BOLTZMANN: f64 = 5.670_374_419e-8;
pub const G0_M_S2: f64 = 9.806_65;
pub const P_ATM_PA: f64 = 101_325.0;
pub const MU0: f64 = 1.256_637_062_12e-6;
/// Sommerfeld value of the Lorenz number [W Ω K⁻²].
pub const LORENZ_L0: f64 = 2.443e-8;

/// Liquid nitrogen at 1 atm (NIST Chemistry WebBook fluid tables; not yet in matdb).
pub mod ln2 {
    pub const T_SAT_1ATM_K: f64 = 77.355;
    pub const LATENT_J_KG: f64 = 199_180.0;
    pub const RHO_LIQUID: f64 = 806.1;
    pub const CP_LIQUID: f64 = 2_041.0;
    pub const K_LIQUID: f64 = 0.1396;
    pub const MU_LIQUID: f64 = 1.61e-4;
    pub const SURFACE_TENSION: f64 = 8.85e-3;
    pub const MOLAR_MASS: f64 = 0.028_013_4;
    /// Triple point; below this pressure-dependent limit the pool freezes.
    pub const T_TRIPLE_K: f64 = 63.15;
    pub const P_TRIPLE_PA: f64 = 12_520.0;
}

/// Saturation temperature of LN₂ at chamber pressure `p_pa` (Clausius–Clapeyron
/// through the 1 atm boiling point, clamped at the triple point).
#[must_use]
pub fn ln2_saturation_temperature_k(p_pa: f64) -> f64 {
    let p = p_pa.max(ln2::P_TRIPLE_PA);
    let b = ln2::LATENT_J_KG * ln2::MOLAR_MASS / R_GAS;
    let inv_t = 1.0 / ln2::T_SAT_1ATM_K - (p / P_ATM_PA).ln() / b;
    (1.0 / inv_t).max(ln2::T_TRIPLE_K)
}

/// Kinds of material a grid column or component can be made of.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum WatchMaterialId {
    Gold18k = 0,
    BrassC36000 = 1,
    CarbonSteel = 2,
    BerylliumCopper = 3,
    Sapphire = 4,
    Nivarox = 5,
}

impl WatchMaterialId {
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::Gold18k => "18k yellow gold Au750/Ag125/Cu125",
            Self::BrassC36000 => "brass C36000",
            Self::CarbonSteel => "carbon spring steel",
            Self::BerylliumCopper => "Glucydur (Be-Cu C17200)",
            Self::Sapphire => "sapphire Al2O3",
            Self::Nivarox => "Nivarox-type Fe-Ni-Cr-Be",
        }
    }
}

/// Room-temperature isotropic solid properties in SI.
#[derive(Debug, Clone, Copy)]
pub struct Solid {
    pub rho: f64,
    pub cp: f64,
    pub k: f64,
    pub e: f64,
    pub alpha: f64,
    pub nu: f64,
}

impl Solid {
    #[must_use]
    pub fn diffusivity(&self) -> f64 {
        self.k / (self.rho * self.cp)
    }
}

/// Debye heat-capacity function `D(Θ/T)` normalised to 1 at high temperature.
/// Debye heat-capacity fraction c_v/(3R) at θ_D/T.
#[must_use]
pub fn debye_fraction(theta_over_t: f64) -> f64 {
    let y = theta_over_t.max(1.0e-6);
    if y > 60.0 {
        // Low-temperature T³ limit: 3·(4π⁴/15) / y³.
        return 77.927_272 / (y * y * y);
    }
    let n = 200;
    let h = y / n as f64;
    let f = |x: f64| {
        if x < 1.0e-8 {
            return x * x; // x⁴eˣ/(eˣ−1)² → x² as x → 0
        }
        let ex = x.exp();
        x.powi(4) * ex / ((ex - 1.0) * (ex - 1.0))
    };
    let mut s = f(0.0) + f(y);
    for i in 1..n {
        let x = i as f64 * h;
        s += if i % 2 == 1 { 4.0 } else { 2.0 } * f(x);
    }
    3.0 / (y * y * y) * s * h / 3.0
}

/// One element of the 18k alloy: matdb-scaled Debye heat capacity.
#[derive(Debug, Clone, Copy)]
struct Element {
    mass_fraction: f64,
    rho: f64,
    theta_debye_k: f64,
    /// Scale that makes the Debye curve hit the matdb room-temperature cp.
    cp_high_t: f64,
    cp_liquid: f64,
    latent_j_kg: f64,
}

impl Element {
    fn new(
        mass_fraction: f64,
        rho: f64,
        cp_room: f64,
        theta: f64,
        cp_liquid: f64,
        latent: f64,
    ) -> Self {
        Self {
            mass_fraction,
            rho,
            theta_debye_k: theta,
            cp_high_t: cp_room / debye_fraction(theta / T_REF_K),
            cp_liquid,
            latent_j_kg: latent,
        }
    }

    fn cp_solid(&self, t: f64) -> f64 {
        // Anharmonic + electronic rise above room temperature (~+10 % per 1000 K
        // for the noble metals, Hultgren 1973).
        let rise = if t > T_REF_K { 1.0 + 1.0e-4 * (t - T_REF_K) } else { 1.0 };
        self.cp_high_t * debye_fraction(self.theta_debye_k / t) * rise
    }
}

/// 18k yellow gold with a full equilibrium enthalpy phase curve.
pub struct GoldAlloy {
    elements: [Element; 3],
    pub rho_solid_ref: f64,
    pub alpha: f64,
    pub e: f64,
    pub nu: f64,
    pub latent_j_kg: f64,
    pub t_solidus_k: f64,
    pub t_liquidus_k: f64,
    pub cp_liquid: f64,
    pub residual_resistivity_ohm_m: f64,
    pub phonon_resistivity_ref_ohm_m: f64,
    curve: EquilibriumEnthalpyPhaseCurve,
    /// (T, h) samples for the inverse map, increasing in both.
    t_to_h: Vec<(f64, f64)>,
    pub h_min: f64,
    pub h_max: f64,
}

/// Evaluated gold state at one specific enthalpy.
#[derive(Debug, Clone, Copy)]
pub struct GoldState {
    pub temperature_k: f64,
    pub liquid_fraction: f64,
    pub density: f64,
}

impl GoldAlloy {
    fn build() -> Self {
        // Mass fractions of 750 yellow gold (3N colour standard ISO 8654).
        let (w_au, w_ag, w_cu) = (0.75, 0.125, 0.125);
        // Debye temperatures (Kittel, Introduction to Solid State Physics, Table 5.1).
        // Liquid heat capacities: Au 33.4, Ag 30.5 J/(mol K) (Hultgren et al. 1973);
        // Cu from the matdb NASA TP-3287 liquid curve.
        let cu_liquid_cp = db::COPPER_LIQUID
            .curve("cp-liquid")
            .map(|c| c.eval(1400.0))
            .unwrap_or(516.0);
        // Latent heats of fusion: Au 12.55 kJ/mol, Ag 11.28 kJ/mol (CRC Handbook);
        // Cu from matdb NASA TP-3287.
        let au = Element::new(
            w_au,
            db::GOLD.require("density"),
            db::GOLD.require("specific-heat"),
            165.0,
            33.4 / 0.196_967,
            12_550.0 / 0.196_967,
        );
        let ag = Element::new(
            w_ag,
            db::SILVER.require("density"),
            db::SILVER.require("specific-heat"),
            225.0,
            30.5 / 0.107_868,
            11_280.0 / 0.107_868,
        );
        // Copper room-temperature cp is not in the c447 row set; 385 J/(kg K) (CRC).
        let cu = Element::new(
            w_cu,
            db::COPPER.require("density"),
            385.0,
            343.0,
            cu_liquid_cp,
            db::COPPER_FUSION.require("latent-heat-fusion"),
        );
        let elements = [au, ag, cu];
        let rho_solid_ref = 1.0 / elements.iter().map(|e| e.mass_fraction / e.rho).sum::<f64>();
        let latent_j_kg = elements.iter().map(|e| e.mass_fraction * e.latent_j_kg).sum::<f64>();
        let cp_liquid = elements.iter().map(|e| e.mass_fraction * e.cp_liquid).sum::<f64>();
        // Volume-fraction mixing of expansion and stiffness.
        let vol = |e: &Element| e.mass_fraction / e.rho * rho_solid_ref;
        let alpha = vol(&au) * db::GOLD.require("thermal-expansion")
            + vol(&ag) * db::SILVER.require("thermal-expansion")
            + vol(&cu) * 16.5e-6;
        let e_mod = vol(&au) * db::GOLD.require("young-modulus")
            + vol(&ag) * db::SILVER.require("young-modulus")
            + vol(&cu) * 117.0e9;
        // Solidus / liquidus of 3N 750 yellow gold: 885 / 895 °C (alloy maker data;
        // not in matdb).
        let t_solidus_k = 1158.15;
        let t_liquidus_k = 1168.15;
        // Pure-gold phonon resistivity from its matdb conductivity via Wiedemann–Franz;
        // 18k yellow gold conducts ~16 % IACS, which fixes the residual term.
        let phonon_resistivity_ref_ohm_m =
            LORENZ_L0 * T_REF_K / db::GOLD.require("thermal-conductivity");
        let residual_resistivity_ohm_m = 1.7241e-8 / 0.16 - phonon_resistivity_ref_ohm_m;

        let rho_at = |t: f64| rho_solid_ref / (1.0 + 3.0 * alpha * (t - T_REF_K));
        let cp_mix = |t: f64| elements.iter().map(|e| e.mass_fraction * e.cp_solid(t)).sum::<f64>();

        // Integrate h(T) = ∫ cp dT from T_REF in fine steps, emitting knots.
        let t_min = 40.0;
        let mut knots: Vec<EnthalpyPhaseKnot> = Vec::new();
        let mut t_to_h: Vec<(f64, f64)> = Vec::new();
        let integrate = |a: f64, b: f64| {
            let n = ((b - a).abs() / 0.5).ceil().max(1.0) as usize;
            let h = (b - a) / n as f64;
            (0..n)
                .map(|i| {
                    let t0 = a + i as f64 * h;
                    0.5 * (cp_mix(t0) + cp_mix(t0 + h)) * h
                })
                .sum::<f64>()
        };
        let mut solid_temps = Vec::new();
        let mut t = t_min;
        while t < t_solidus_k - 1.0e-9 {
            solid_temps.push(t);
            t += if t < 400.0 { 10.0 } else { 40.0 };
        }
        solid_temps.push(t_solidus_k);
        for &ts in &solid_temps {
            let h = -integrate(ts, T_REF_K);
            knots.push(EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: h,
                temperature_k: ts,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: rho_at(ts),
            });
            t_to_h.push((ts, h));
        }
        let h_sol = t_to_h.last().map(|p| p.1).unwrap_or(0.0);
        let cp_sol = cp_mix(t_solidus_k);
        let rho_liq_at_liquidus = rho_at(t_solidus_k) * (1.0 - 0.051); // ~5 % volume rise on melting
        for i in 1..=4 {
            let f = i as f64 / 4.0;
            let ts = t_solidus_k + f * (t_liquidus_k - t_solidus_k);
            let cp_mush = (1.0 - f) * cp_sol + f * cp_liquid;
            let h = h_sol + f * latent_j_kg + cp_mush * (ts - t_solidus_k);
            let rho = 1.0 / ((1.0 - f) / rho_at(ts) + f / rho_liq_at_liquidus);
            knots.push(EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: h,
                temperature_k: ts,
                liquid_mass_fraction: f,
                bulk_density_kg_m3: rho,
            });
            t_to_h.push((ts, h));
        }
        let h_liq = t_to_h.last().map(|p| p.1).unwrap_or(0.0);
        for &tl in &[1600.0, 2400.0, 3600.0] {
            let h = h_liq + cp_liquid * (tl - t_liquidus_k);
            // Liquid gold dρ/dT ≈ −1.2 kg/(m³ K) (Iida & Guthrie), scaled to the alloy.
            let rho = rho_liq_at_liquidus - 0.96 * (tl - t_liquidus_k);
            knots.push(EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: h,
                temperature_k: tl,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: rho,
            });
            t_to_h.push((tl, h));
        }
        let h_min = knots[0].specific_enthalpy_j_kg;
        let h_max = knots[knots.len() - 1].specific_enthalpy_j_kg;
        let mut hasher = DomainHasher::new("fs-pocket-watch.gold-18k-yellow.v2");
        hasher.update(db::GOLD.id.as_bytes());
        hasher.update(db::SILVER.id.as_bytes());
        hasher.update(db::COPPER.id.as_bytes());
        let curve = EquilibriumEnthalpyPhaseCurve::try_new(hasher.finalize(), knots)
            .expect("18k gold phase curve knots are monotone and span solid→liquid");
        Self {
            elements,
            rho_solid_ref,
            alpha,
            e: e_mod,
            nu: 0.42,
            latent_j_kg,
            t_solidus_k,
            t_liquidus_k,
            cp_liquid,
            residual_resistivity_ohm_m,
            phonon_resistivity_ref_ohm_m,
            curve,
            t_to_h,
            h_min,
            h_max,
        }
    }

    /// Temperature, liquid fraction and density at specific enthalpy `h`
    /// (J/kg relative to the solid at 293.15 K).
    #[must_use]
    pub fn state(&self, h: f64) -> GoldState {
        let hc = h.clamp(self.h_min, self.h_max);
        match self.curve.values_at_specific_enthalpy(hc) {
            Ok(k) => GoldState {
                temperature_k: k.temperature_k,
                liquid_fraction: k.liquid_mass_fraction,
                density: k.bulk_density_kg_m3,
            },
            Err(_) => GoldState {
                temperature_k: T_REF_K,
                liquid_fraction: 0.0,
                density: self.rho_solid_ref,
            },
        }
    }

    /// Inverse map: specific enthalpy of fully solid / liquid gold at `t_k`
    /// (inside the mushy interval returns the equilibrium value).
    #[must_use]
    pub fn enthalpy_at(&self, t_k: f64) -> f64 {
        let k = &self.t_to_h;
        let t = t_k.clamp(k[0].0, k[k.len() - 1].0);
        let i = k.partition_point(|p| p.0 <= t).clamp(1, k.len() - 1);
        let (t0, h0) = k[i - 1];
        let (t1, h1) = k[i];
        if t1 <= t0 {
            return h1;
        }
        h0 + (h1 - h0) * (t - t0) / (t1 - t0)
    }

    /// (T, h) knots of the inverse map, increasing in both.
    #[must_use]
    pub fn t_to_h(&self) -> &[(f64, f64)] {
        &self.t_to_h
    }

    /// Equilibrium liquid fraction at temperature `t_k` (lever rule across the
    /// mushy interval, which the enthalpy curve treats as linear in T).
    #[must_use]
    pub fn liquid_fraction_at(&self, t_k: f64) -> f64 {
        ((t_k - self.t_solidus_k) / (self.t_liquidus_k - self.t_solidus_k)).clamp(0.0, 1.0)
    }

    /// Solid specific heat at `t_k` (Debye-shaped, matdb-scaled).
    #[must_use]
    pub fn cp_solid(&self, t_k: f64) -> f64 {
        self.elements.iter().map(|e| e.mass_fraction * e.cp_solid(t_k.max(20.0))).sum()
    }

    /// Thermal conductivity from Wiedemann–Franz with Matthiessen's rule:
    /// residual alloy scattering plus phonon scattering ∝ T, and the ~1.5×
    /// resistivity jump of the liquid.
    #[must_use]
    pub fn conductivity(&self, t_k: f64, liquid_fraction: f64) -> f64 {
        let t = t_k.max(30.0);
        let rho_solid =
            self.residual_resistivity_ohm_m + self.phonon_resistivity_ref_ohm_m * t / T_REF_K;
        let rho = rho_solid * (1.0 + 0.5 * liquid_fraction.clamp(0.0, 1.0));
        LORENZ_L0 * t / rho
    }

    /// Liquid viscosity: Arrhenius with gold's 15.9 kJ/mol activation energy and
    /// 5.38 mPa·s at gold's melting point (Iida & Guthrie, *Physical Properties of
    /// Liquid Metals*), continued to the alloy's liquidus.
    #[must_use]
    pub fn liquid_viscosity(&self, t_k: f64) -> f64 {
        let t = t_k.max(self.t_solidus_k);
        5.38e-3 * (15_900.0 / R_GAS * (1.0 / t - 1.0 / 1337.33)).exp()
    }

    /// Effective viscosity of the mushy slurry: Krieger–Dougherty divergence as
    /// the solid fraction approaches dendrite coherency (φ_max ≈ 0.6); fully
    /// solid gold is returned as a creep viscosity (Arrhenius power-law creep at
    /// the low stresses of self-weight, Frost & Ashby deformation maps).
    #[must_use]
    pub fn effective_viscosity(&self, t_k: f64, liquid_fraction: f64) -> f64 {
        let fs = 1.0 - liquid_fraction.clamp(0.0, 1.0);
        let phi_max = 0.6;
        if fs < phi_max * 0.98 {
            self.liquid_viscosity(t_k) * (1.0 - fs / phi_max).powf(-2.5 * phi_max)
        } else {
            // Diffusional creep viscosity η ≈ k_B T d² / (42 Ω D_v) collapses to an
            // Arrhenius law; ~1e9 Pa·s just below the solidus, ~1e15 at 900 K.
            let t = t_k.clamp(60.0, self.t_solidus_k);
            1.0e9 * (174_000.0 / R_GAS * (1.0 / t - 1.0 / self.t_solidus_k)).exp().min(1.0e12)
        }
    }

    /// Liquid surface tension (Au 1.145 N/m at its melting point, dσ/dT = −0.52 mN/(m K);
    /// silver and copper pull the alloy to ~1.10 N/m at its liquidus).
    #[must_use]
    pub fn surface_tension(&self, t_k: f64) -> f64 {
        (1.10 - 0.52e-3 * (t_k - self.t_liquidus_k)).clamp(0.6, 1.3)
    }

    /// Total hemispherical emissivity: polished gold ~0.02, rising with copper
    /// oxide (Cu₂O/CuO) on the 12.5 % copper as `oxide_nm` grows, and ~0.10 for the
    /// clean liquid (Touloukian & DeWitt, *Thermal Radiative Properties*, Vol. 7).
    #[must_use]
    pub fn emissivity(&self, liquid_fraction: f64, oxide_nm: f64) -> f64 {
        self.emissivity_with_oxide(liquid_fraction, Self::oxide_emissivity(oxide_nm))
    }

    /// The oxide film's share of [`Self::emissivity`].
    #[must_use]
    pub fn oxide_emissivity(oxide_nm: f64) -> f64 {
        0.5 * (1.0 - (-oxide_nm / 400.0).exp())
    }

    /// [`Self::emissivity`] given the film's share from [`Self::oxide_emissivity`].
    #[must_use]
    pub fn emissivity_with_oxide(&self, liquid_fraction: f64, oxide_part: f64) -> f64 {
        let clean = 0.02 + 0.08 * liquid_fraction.clamp(0.0, 1.0);
        (clean + oxide_part).clamp(0.02, 0.65)
    }
}

/// Watch oil: Pennzane SHF X-2000 (multiply-alkylated cyclopentane) from matdb.
#[derive(Debug, Clone, Copy)]
pub struct WatchOil {
    /// VFT parameters for kinematic viscosity: ln ν = a + b / (T − t0), ν in m²/s.
    pub vft_a: f64,
    pub vft_b: f64,
    pub vft_t0: f64,
    pub density_ref: f64,
    pub pour_point_k: f64,
    pub flash_point_k: f64,
}

impl WatchOil {
    fn build() -> Self {
        let p = db::PENNZANE;
        let pts = [
            (233.15, p.require("kinematic-viscosity-minus40degC")),
            (313.15, p.require("kinematic-viscosity-40degC")),
            (373.15, p.require("kinematic-viscosity-100degC")),
        ];
        // Exact 3-point VFT: solve for t0 by bisection on the ratio of log differences.
        let l = pts.map(|(_, nu)| nu.ln());
        let ratio = (l[0] - l[1]) / (l[1] - l[2]);
        let g = |t0: f64| {
            let d01 = 1.0 / (pts[0].0 - t0) - 1.0 / (pts[1].0 - t0);
            let d12 = 1.0 / (pts[1].0 - t0) - 1.0 / (pts[2].0 - t0);
            d01 / d12 - ratio
        };
        let (mut lo, mut hi) = (0.0_f64, pts[0].0 - 1.0);
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if g(lo) * g(mid) <= 0.0 {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        let t0 = 0.5 * (lo + hi);
        let b = (l[1] - l[2]) / (1.0 / (pts[1].0 - t0) - 1.0 / (pts[2].0 - t0));
        let a = l[1] - b / (pts[1].0 - t0);
        Self {
            vft_a: a,
            vft_b: b,
            vft_t0: t0,
            density_ref: p.require("density-25degC"),
            pour_point_k: p.require("pour-point"),
            flash_point_k: p.require("flash-point"),
        }
    }

    #[must_use]
    pub fn density(&self, t_k: f64) -> f64 {
        // Typical PAO/MAC thermal expansion 7e-4 /K about 25 °C.
        self.density_ref * (1.0 - 7.0e-4 * (t_k - 298.15)).clamp(0.6, 1.25)
    }

    /// Dynamic viscosity [Pa·s]. Below the VFT divergence the oil is a glass;
    /// the value is capped at 1e12 Pa·s (glass-transition viscosity).
    #[must_use]
    pub fn viscosity(&self, t_k: f64) -> f64 {
        let dt = t_k - self.vft_t0;
        if dt <= 1.0 {
            return 1.0e12;
        }
        let nu = (self.vft_a + self.vft_b / dt).exp();
        (nu * self.density(t_k)).min(1.0e12)
    }

    /// Surface tension of a cyclopentane/PAO base stock (~30 mN/m at 25 °C,
    /// −0.07 mN/(m K)).
    #[must_use]
    pub fn surface_tension(&self, t_k: f64) -> f64 {
        (0.030 - 7.0e-5 * (t_k - 298.15)).clamp(0.01, 0.05)
    }
}

/// Sintered NdFeB N42 (Jinshan pristine endpoints from matdb).
#[derive(Debug, Clone, Copy)]
pub struct MagnetN42 {
    pub br_25c_t: f64,
    pub br_temp_coeff: f64,
    pub hcj_25c_a_m: f64,
    pub hcj_temp_coeff: f64,
}

impl MagnetN42 {
    fn build() -> Self {
        let p = db::N42;
        Self {
            br_25c_t: p.require("remanence-25degC"),
            br_temp_coeff: p.require("remanence-temperature-coefficient"),
            hcj_25c_a_m: p.require("intrinsic-coercivity-25degC"),
            hcj_temp_coeff: p.require("intrinsic-coercivity-temperature-coefficient"),
        }
    }

    /// Reversible remanence at magnet temperature `t_k` (zero above the 583 K Curie point).
    #[must_use]
    pub fn remanence_t(&self, t_k: f64) -> f64 {
        if t_k >= 583.0 {
            return 0.0;
        }
        (self.br_25c_t * (1.0 + self.br_temp_coeff * (t_k - 298.15))).max(0.0)
    }
}

/// Hairspring alloy model (balance-spring thermoelastic and magnetic behaviour).
#[derive(Debug, Clone, Copy)]
pub struct HairspringAlloy {
    pub id: WatchMaterialId,
    pub e_ref: f64,
    pub rho: f64,
    pub alpha: f64,
    /// First-order thermoelastic coefficient (1/E)(dE/dT) at T_REF [1/K].
    pub thermoelastic_1: f64,
    /// Second-order coefficient [1/K²] (Elinvar "middle temperature error").
    pub thermoelastic_2: f64,
    /// Internal friction (loss tangent) of the strip.
    pub loss_tangent: f64,
    /// Saturation magnetisation [T] and retentivity fraction for the play model.
    pub saturation_t: f64,
    pub coercive_a_m: f64,
    /// Relative stiffness change from coil–coil magnetic attraction at saturation.
    pub magnetic_stiffening_at_sat: f64,
}

impl HairspringAlloy {
    /// Young's modulus at temperature `t_k`.
    #[must_use]
    pub fn youngs_modulus(&self, t_k: f64) -> f64 {
        let dt = t_k - T_REF_K;
        self.e_ref * (1.0 + self.thermoelastic_1 * dt + self.thermoelastic_2 * dt * dt).max(0.2)
    }
}

/// Blackbody colour lookup: linear-sRGB chromaticity (max channel = 1) and
/// luminance of a unit-emissivity surface, from 500 K to 3500 K.
pub struct BlackbodyLut {
    pub t0: f64,
    pub dt: f64,
    pub rgb: Vec<[f64; 3]>,
    pub luminance_cd_m2: Vec<f64>,
}

impl BlackbodyLut {
    fn build() -> Self {
        let t0 = 500.0;
        let dt = 25.0;
        let n = 121;
        let mut hasher = DomainHasher::new("fs-pocket-watch.blackbody-lut.v1");
        hasher.update(b"gray-unit-emissivity");
        let id = hasher.finalize();
        let mut rgb = Vec::with_capacity(n);
        let mut lum = Vec::with_capacity(n);
        for i in 0..n {
            let t = t0 + dt * i as f64;
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
                    rgb.push([1.0, 0.3, 0.05]);
                    lum.push(0.0);
                }
            }
        }
        Self { t0, dt, rgb, luminance_cd_m2: lum }
    }

    /// Colour and luminance of a gray emitter of emissivity `eps` at `t_k`.
    #[must_use]
    pub fn sample(&self, t_k: f64, eps: f64) -> ([f64; 3], f64) {
        if t_k < self.t0 {
            return ([0.0; 3], 0.0);
        }
        let x = ((t_k - self.t0) / self.dt).min((self.rgb.len() - 1) as f64);
        let i = (x.floor() as usize).min(self.rgb.len() - 2);
        let f = x - i as f64;
        let mut c = [0.0; 3];
        for ch in 0..3 {
            c[ch] = self.rgb[i][ch] * (1.0 - f) + self.rgb[i + 1][ch] * f;
        }
        let l = self.luminance_cd_m2[i] * (1.0 - f) + self.luminance_cd_m2[i + 1] * f;
        (c, l * eps)
    }

    /// Display-referred incandescence: chromaticity × a log-compressed exposure
    /// of the true luminance (0 below visibility, ~1 at bright orange heat).
    #[must_use]
    pub fn display_rgb(&self, t_k: f64, eps: f64) -> [f64; 3] {
        let (c, l) = self.sample(t_k, eps);
        // Draper point (~798 K) is where glow becomes visible in the dark: ~1e-3 cd/m².
        let exposure = ((l / 1.0e-3).max(1.0).ln() / (1.0e5_f64).ln()).clamp(0.0, 2.0);
        [c[0] * exposure, c[1] * exposure, c[2] * exposure]
    }
}

/// All material models the watch uses, built once.
pub struct Library {
    pub gold: GoldAlloy,
    pub brass: Solid,
    pub steel: Solid,
    pub becu: Solid,
    pub sapphire: Solid,
    pub invar_e: Curve,
    pub invar_k: Curve,
    pub oil: WatchOil,
    pub n42: MagnetN42,
    pub carbon_steel_hairspring: HairspringAlloy,
    pub nivarox_hairspring: HairspringAlloy,
    pub blackbody: BlackbodyLut,
}

impl Library {
    fn build() -> Self {
        let gold = GoldAlloy::build();
        let brass = Solid {
            rho: db::BRASS_C36000.require("rho"),
            cp: db::BRASS_C36000.require("cp"),
            k: db::BRASS_C36000.require("k"),
            e: db::BRASS_C36000.require("e"),
            alpha: db::BRASS_C36000.require("alpha-mean"),
            nu: 0.31,
        };
        let steel = Solid {
            rho: db::MUSIC_WIRE.require("density"),
            cp: db::STEEL_1018.require("heat-capacity"),
            k: db::STEEL_1018.require("conductivity"),
            e: db::MUSIC_WIRE.require("young-modulus"),
            alpha: 11.5e-6,
            nu: db::STEEL_1018.require("poisson-ratio"),
        };
        let becu = Solid {
            rho: db::BERYLLIUM_COPPER.require("density"),
            cp: db::BERYLLIUM_COPPER.require("heat-capacity"),
            k: db::BERYLLIUM_COPPER.require("conductivity"),
            e: db::BERYLLIUM_COPPER.require("young-modulus"),
            alpha: 17.8e-6,
            nu: 0.3,
        };
        let sapphire = Solid {
            rho: db::SAPPHIRE.require("density"),
            cp: db::SAPPHIRE.require("heat-capacity"),
            k: db::SAPPHIRE.require("conductivity"),
            e: 345.0e9,
            alpha: 5.6e-6,
            nu: 0.29,
        };
        let invar_e = db::INVAR.curve("young-modulus").expect("invar E(T) curve");
        let invar_k = db::INVAR.curve("conductivity").expect("invar k(T) curve");
        // Carbon-steel balance spring: E and loss tangent from the LIGO music-wire
        // rows; thermoelastic coefficient −2.6e-4 /K (Ashby & Jones; Reed-Hill),
        // which makes an uncompensated watch lose ~11 s/day per kelvin.
        let carbon_steel_hairspring = HairspringAlloy {
            id: WatchMaterialId::CarbonSteel,
            e_ref: steel.e,
            rho: steel.rho,
            alpha: steel.alpha,
            thermoelastic_1: -2.6e-4,
            thermoelastic_2: -2.0e-7,
            loss_tangent: db::MUSIC_WIRE.require("loss-tangent-465um"),
            saturation_t: 2.0,
            coercive_a_m: 4_000.0,
            magnetic_stiffening_at_sat: 0.020,
        };
        // Nivarox-type Fe-Ni-Cr-Be: the alloy maker tunes the Invar-type anomaly so
        // the oscillator's first-order temperature coefficient vanishes against the
        // Glucydur balance, i.e. TEC = 2 α_balance − 3 α_spring. The residual
        // second-order "middle temperature error" is ~2 s/day across ±15 K
        // (COSC test data), i.e. (1/E) d²E/dT² ≈ −2e-7 /K².
        let alpha_niv = 8.0e-6;
        let tec = 2.0 * becu.alpha - 3.0 * alpha_niv;
        let nivarox_hairspring = HairspringAlloy {
            id: WatchMaterialId::Nivarox,
            e_ref: 186.0e9,
            rho: 8_150.0,
            alpha: alpha_niv,
            thermoelastic_1: tec,
            thermoelastic_2: -2.0e-7,
            loss_tangent: 2.5e-4,
            saturation_t: 0.9,
            coercive_a_m: 300.0,
            magnetic_stiffening_at_sat: 0.0015,
        };
        Self {
            gold,
            brass,
            steel,
            becu,
            sapphire,
            invar_e,
            invar_k,
            oil: WatchOil::build(),
            n42: MagnetN42::build(),
            carbon_steel_hairspring,
            nivarox_hairspring,
            blackbody: BlackbodyLut::build(),
        }
    }

    pub fn solid(&self, id: WatchMaterialId) -> Solid {
        match id {
            WatchMaterialId::Gold18k => Solid {
                rho: self.gold.rho_solid_ref,
                cp: self.gold.cp_solid(T_REF_K),
                k: self.gold.conductivity(T_REF_K, 0.0),
                e: self.gold.e,
                alpha: self.gold.alpha,
                nu: self.gold.nu,
            },
            WatchMaterialId::BrassC36000 => self.brass,
            WatchMaterialId::CarbonSteel => self.steel,
            WatchMaterialId::BerylliumCopper => self.becu,
            WatchMaterialId::Sapphire => self.sapphire,
            WatchMaterialId::Nivarox => Solid {
                rho: self.nivarox_hairspring.rho,
                cp: 460.0,
                k: 12.0,
                e: self.nivarox_hairspring.e_ref,
                alpha: self.nivarox_hairspring.alpha,
                nu: 0.3,
            },
        }
    }
}

static LIBRARY: OnceLock<Library> = OnceLock::new();

/// The shared material library (built on first use).
pub fn library() -> &'static Library {
    LIBRARY.get_or_init(Library::build)
}

/// Air (or any `GasSpec`) properties at temperature and pressure; `None` in vacuum.
#[must_use]
pub fn air_state(t_k: f64, p_pa: f64) -> Option<fs_material::gas::GasState> {
    if p_pa < 1.0 {
        return None;
    }
    fs_material::gas::GasState::try_new(
        &fs_material::gas::GasSpec::dry_air_ussa1976(),
        t_k.clamp(60.0, 1990.0),
        p_pa,
    )
    .ok()
}

/// Water-ice saturation vapour pressure [Pa] (`fs_material::sublimation`), 0 when
/// the surface is too cold for the curve (vapour pressure is then negligible).
#[must_use]
pub fn ice_saturation_pressure_pa(t_k: f64) -> f64 {
    ice_saturation_table().pressure_pa(t_k)
}

/// `ln p_sat(T)` of water ice tabulated once from `fs_material::sublimation`
/// (the frost model asks for it per cell per thermal substep). Linear
/// interpolation in `ln p` on a 0.22 K pitch is good to ~2e-4 relative even at
/// 77 K, where Clausius–Clapeyron curvature is largest.
struct IceSaturationTable {
    ln_p: Vec<f64>,
}

const ICE_TABLE_T0_K: f64 = 50.0;
const ICE_TABLE_T1_K: f64 = 273.16;
const ICE_TABLE_N: usize = 1024;

impl IceSaturationTable {
    fn build() -> Self {
        let dt = (ICE_TABLE_T1_K - ICE_TABLE_T0_K) / (ICE_TABLE_N - 1) as f64;
        let ln_p = (0..ICE_TABLE_N)
            .map(|i| {
                let t = ICE_TABLE_T0_K + dt * i as f64;
                let p = fs_material::sublimation::WaterIceSublimation::saturation_pressure_pa(t).unwrap_or(0.0);
                if p > 0.0 { p.ln() } else { f64::NEG_INFINITY }
            })
            .collect();
        Self { ln_p }
    }

    fn pressure_pa(&self, t_k: f64) -> f64 {
        if !(t_k >= ICE_TABLE_T0_K) {
            return 0.0;
        }
        if t_k >= ICE_TABLE_T1_K {
            return fs_material::sublimation::WaterIceSublimation::saturation_pressure_pa(t_k).unwrap_or(0.0);
        }
        let x = (t_k - ICE_TABLE_T0_K) / (ICE_TABLE_T1_K - ICE_TABLE_T0_K) * (ICE_TABLE_N - 1) as f64;
        let i = (x as usize).min(ICE_TABLE_N - 2);
        let f = x - i as f64;
        let (a, b) = (self.ln_p[i], self.ln_p[i + 1]);
        if !a.is_finite() || !b.is_finite() {
            return 0.0;
        }
        (a + f * (b - a)).exp()
    }
}

fn ice_saturation_table() -> &'static IceSaturationTable {
    static TABLE: OnceLock<IceSaturationTable> = OnceLock::new();
    TABLE.get_or_init(IceSaturationTable::build)
}

/// Parabolic (Wagner) growth constant of the temper oxide on carbon steel
/// [nm²/s]: k_p = A exp(−Q/RT) with Q = 110 kJ/mol, A fitted so 573 K gives
/// first-order blue (~70 nm) in five minutes, 523 K straw in the same time.
#[must_use]
pub fn steel_temper_oxide_rate_nm2_s(t_k: f64) -> f64 {
    if t_k < 380.0 {
        return 0.0;
    }
    1.7e11 * (-110_000.0 / (R_GAS * t_k)).exp()
}

/// Thin-film reflectance of temper-blued steel at the three display primaries,
/// through `fs_render::conductor::ConductorOptics` with a coherent Fe₃O₄ coating.
#[must_use]
pub fn blued_steel_interference_rgb(oxide_thickness_nm: f64, cos_theta: f64) -> [f64; 3] {
    const WAVELENGTHS_NM: [f64; CONDUCTOR_IOR_SAMPLE_COUNT] =
        [380.0, 430.0, 480.0, 530.0, 580.0, 630.0, 680.0, 730.0, 780.0];
    // Iron optical constants (Johnson & Christy 1974), linearly interpolated.
    let fe = |wl: f64| {
        let t = ((wl - 380.0) / 400.0).clamp(0.0, 1.0);
        (2.6 + 0.35 * t, 2.9 + 0.9 * t)
    };
    let samples = WAVELENGTHS_NM.map(|wl| {
        let (n, k) = fe(wl);
        ConductorIorSample::try_new(wl, n, k).expect("valid iron n,k knot")
    });
    let mut hasher = DomainHasher::new("fs-pocket-watch.iron-johnson-christy.v1");
    hasher.update(b"Fe");
    let Ok(source) =
        ConductorSource::try_new(hasher.finalize(), ConductorDataStatus::MaterialStateClaims)
    else {
        return [0.6; 3];
    };
    let Ok(mut optics) = ConductorOptics::try_new(samples, source) else {
        return [0.6; 3];
    };
    if oxide_thickness_nm > 0.5 {
        let mut coat_hasher = DomainHasher::new("fs-pocket-watch.fe3o4-coating.v1");
        coat_hasher.update(&oxide_thickness_nm.to_bits().to_le_bytes());
        if let (Ok(cauchy), Ok(coat_source)) = (
            CauchyIor::try_new(2.28, 0.018, 0.0),
            ConductorSource::try_new(coat_hasher.finalize(), ConductorDataStatus::MaterialStateClaims),
        ) && let Ok(coating) = CoherentDielectricCoating::try_new(
            cauchy,
            (oxide_thickness_nm * 1.0e-9).clamp(0.0, 1.0e-4),
            coat_source,
        ) {
            optics = optics.with_coating(coating);
        }
    }
    let c = cos_theta.clamp(0.05, 1.0);
    let lambdas_nm = [630.0, 532.0, 465.0];
    lambdas_nm.map(|wl| optics.fresnel(wl, 1.0, c).unwrap_or(0.55).clamp(0.02, 0.98))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gold_alloy_properties_come_out_of_the_mixing_rules() {
        let g = &library().gold;
        // 18k yellow gold is ~15.5 g/cm³ and conducts ~60–70 W/(m K).
        assert!((g.rho_solid_ref - 15_430.0).abs() < 60.0, "{}", g.rho_solid_ref);
        let k = g.conductivity(T_REF_K, 0.0);
        assert!(k > 55.0 && k < 75.0, "k = {k}");
        assert!(g.latent_j_kg > 80_000.0 && g.latent_j_kg < 95_000.0);
        let cp = g.cp_solid(T_REF_K);
        assert!(cp > 165.0 && cp < 185.0, "cp = {cp}");
        // Debye: much lower heat capacity at LN₂ temperature.
        assert!(g.cp_solid(77.0) < 0.8 * cp);
    }

    #[test]
    fn gold_phase_curve_round_trips_and_has_a_latent_plateau() {
        let g = &library().gold;
        assert!(g.enthalpy_at(T_REF_K).abs() < 1.0);
        for t in [80.0, 293.15, 700.0, 1100.0, 1400.0, 2000.0] {
            let s = g.state(g.enthalpy_at(t));
            assert!((s.temperature_k - t).abs() < 0.5, "{t} -> {}", s.temperature_k);
        }
        let h_sol = g.enthalpy_at(g.t_solidus_k);
        let mid = g.state(h_sol + 0.5 * g.latent_j_kg);
        assert!(mid.liquid_fraction > 0.3 && mid.liquid_fraction < 0.7);
        assert!(mid.temperature_k > g.t_solidus_k && mid.temperature_k < g.t_liquidus_k);
        assert!(g.state(g.enthalpy_at(1500.0)).liquid_fraction == 1.0);
    }

    #[test]
    fn pennzane_vft_reproduces_its_three_matdb_points_and_glasses_at_ln2() {
        let oil = &library().oil;
        let nu = |t: f64| oil.viscosity(t) / oil.density(t);
        assert!((nu(313.15) / 107.0e-6 - 1.0).abs() < 0.03);
        assert!((nu(373.15) / 14.3e-6 - 1.0).abs() < 0.03);
        assert!((nu(233.15) / 80_500.0e-6 - 1.0).abs() < 0.05);
        assert!(oil.viscosity(77.0) >= 1.0e11);
    }

    #[test]
    fn nivarox_cancels_the_first_order_temperature_error_and_steel_does_not() {
        let lib = library();
        let rate_per_k = |hs: &HairspringAlloy| {
            // ω ∝ sqrt(E b h³ / (L I)), I ∝ (1 + α_bal ΔT)²
            0.5 * (hs.thermoelastic_1 + 3.0 * hs.alpha - 2.0 * lib.becu.alpha) * 86_400.0
        };
        assert!(rate_per_k(&lib.carbon_steel_hairspring) < -9.0);
        assert!(rate_per_k(&lib.nivarox_hairspring).abs() < 0.05);
    }

    #[test]
    fn ln2_boils_colder_in_a_partial_vacuum() {
        assert!((ln2_saturation_temperature_k(P_ATM_PA) - 77.355).abs() < 1e-6);
        let t = ln2_saturation_temperature_k(0.1 * P_ATM_PA);
        assert!(t > 63.0 && t < 68.0, "{t}");
    }

    #[test]
    fn ice_saturation_table_matches_the_sublimation_law() {
        let mut worst: f64 = 0.0;
        for i in 0..400 {
            let t = 60.0 + 0.53 * f64::from(i) + 0.017;
            let exact = fs_material::sublimation::WaterIceSublimation::saturation_pressure_pa(t).unwrap_or(0.0);
            let tab = ice_saturation_pressure_pa(t);
            if exact > 0.0 {
                worst = worst.max((tab / exact - 1.0).abs());
            } else {
                assert!(tab == 0.0, "T={t}: table {tab}, law has no value");
            }
        }
        assert!(worst < 1.0e-3, "worst relative error {worst}");
        assert!((ice_saturation_pressure_pa(273.15) - 611.0).abs() < 5.0);
    }

    #[test]
    fn blackbody_lut_glows_red_then_white() {
        let bb = &library().blackbody;
        let (c800, l800) = bb.sample(800.0, 1.0);
        let (c2500, l2500) = bb.sample(2500.0, 1.0);
        assert!(c800[0] > 5.0 * c800[2]);
        assert!(c2500[2] > 0.04 && c2500[2] < 0.2 && c2500[1] > 0.25, "{c2500:?}");
        assert!(l2500 > 1.0e4 * l800);
    }

    #[test]
    fn blued_steel_thin_film_is_blue_at_first_order() {
        let rgb = blued_steel_interference_rgb(72.0, 1.0);
        assert!(rgb[2] > rgb[0], "{rgb:?}");
    }
}
