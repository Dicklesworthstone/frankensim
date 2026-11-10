//! Thermal, fluid, cryogenic and magnetic fields of the FS-1903 pocket watch.
//!
//! The watch is a disc of columns on a 64 × 64 grid in its own plane (x towards
//! 3 o'clock, y towards 12, z out of the dial). Each column carries an
//! enthalpy over its layers — the 18k case back and band, the brass movement,
//! any gold that has run into it — so melting and freezing follow the gold's
//! equilibrium phase curve, and the brass and gold heat capacities follow
//! their Debye shapes down to liquid-nitrogen temperatures.
//!
//! Heat moves by in-plane conduction (Wiedemann–Franz conductivity of the
//! gold), and leaves or enters through:
//! - flame impingement from the torch: a stagnation-point film coefficient
//!   between the adiabatic flame temperature and the surface;
//! - the compressed-air jet, natural convection (Churchill–Chu correlations
//!   from `fs-convection`, vanishing in vacuum) and grey-body radiation with an
//!   emissivity that rises as the copper in the alloy oxidises;
//! - liquid nitrogen boiling on the full pool-boiling curve: Rohsenow nucleate
//!   boiling up to the Zuber critical flux, transition boiling, and Berenson
//!   film boiling above the Leidenfrost point, with T_sat from
//!   Clausius–Clapeyron at the bell-jar pressure;
//! - frost deposited from the room's water vapour onto surfaces below the ice
//!   saturation point (mass transfer by the heat–mass analogy).
//!
//! Molten gold and poured oil flow as shallow layers by virtual pipes driven by
//! gravity (any direction), the hydrostatic head of the layer, capillary
//! pressure, and the stagnation pressure of the torch and air jets, slowed by
//! laminar and turbulent wall friction. The gold's viscosity is the
//! Krieger-mushy and creep law of [`crate::materials::GoldAlloy`], so gold
//! stops where it freezes. Liquid that reaches the edge of the case drips off.
//!
//! The magnet is an N42 dipole whose remanence follows its matdb temperature
//! curve; its field gradient pulls on the steel pallet lever (a torque that
//! depends on where the tourbillon cage has carried the lever), its gradient
//! damps the balance rim by eddy currents, and its field magnetises the
//! hairspring through a play-operator hysteresis that an AC degauss erases.

use core::f64::consts::{PI, TAU};
use std::sync::OnceLock;

use fs_convection::{CorrelationId, CorrelationInputs, ThermalConductivity, evaluate as eval_convection};
use fs_qty::Length as ConvLength;

use crate::horology::{MovementEnvironment, PartField};
use crate::materials::{
    G0_M_S2, GoldAlloy, MU0, P_ATM_PA, STEFAN_BOLTZMANN, T_REF_K, air_state, blued_steel_interference_rgb, debye_fraction,
    ice_saturation_pressure_pa, library, ln2, ln2_saturation_temperature_k, steel_temper_oxide_rate_nm2_s,
};

pub const NODE_COUNT: usize = 12;
pub const GRID_N: usize = 64;
pub const GRID_CELLS: usize = GRID_N * GRID_N;
pub const GRID_HALF_SPAN_MM: f64 = 28.8;
pub const GRID_DX_MM: f64 = 2.0 * GRID_HALF_SPAN_MM / GRID_N as f64;
const DX: f64 = GRID_DX_MM * 1.0e-3;
const CELL_AREA: f64 = DX * DX;

pub const CASE_SHELL_N_THETA: usize = 64;
pub const CASE_SHELL_N_Z: usize = 8;
pub const CASE_SHELL_NODES: usize = CASE_SHELL_N_THETA * CASE_SHELL_N_Z;

/// Layout of the movement in the watch frame [mm].
pub const TOURBILLON_CENTER_MM: [f64; 2] = [0.0, -11.0];
pub const BARREL_CENTER_MM: [f64; 2] = [0.0, 11.5];
const CASE_INNER_R_MM: f64 = 22.2;
const CASE_OUTER_R_MM: f64 = 27.6;
const WELL_R_MM: f64 = 7.5;
/// Floor of the tourbillon aperture below the dial surface (z = 0) [mm].
const WELL_FLOOR_MM: f64 = -1.0;
/// Open-face movement: liquids may pile this high above the dial [m].
const OPEN_CEILING_M: f64 = 6.0e-3;
/// Layers of an interior column: 18k case back, brass plate and bridges.
const BACK_GOLD_M: f64 = 0.8e-3;
const BRASS_M: f64 = 2.9e-3;
const BRASS_DEBYE_K: f64 = 320.0;
/// Room humidity for frost (50 % RH at 20 °C).
const ROOM_RH: f64 = 0.5;
const WATER_P_SAT_293_PA: f64 = 2_339.0;
const M_WATER: f64 = 0.018_015;
const R_GAS: f64 = 8.314_462_618;
const FROST_DENSITY: f64 = 200.0;
/// Degauss coil: peak field [T], decay constant and run time [s].
const DEGAUSS_PEAK_T: f64 = 0.040;
const DEGAUSS_TAU_S: f64 = 0.35;
const DEGAUSS_DURATION_S: f64 = 3.0;
const ICE_LATENT_SUB_J_KG: f64 = 2.834e6;
const ICE_LATENT_FUS_J_KG: f64 = 3.34e5;

/// Canonical indices into the 12 component regions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(usize)]
pub enum WatchNodeId {
    CaseTop18kGold = 0,
    CaseLeft18kGold = 1,
    CaseRight18kGold = 2,
    CaseBottom18kGold = 3,
    CrystalSapphire = 4,
    MainplateBrass = 5,
    BarrelMainspring = 6,
    GoingTrainWheels = 7,
    TourbillonCage = 8,
    EscapementForkWheel = 9,
    BalanceWheelRim = 10,
    HairspringBreguet = 11,
}

/// Aggregated state of one of the 12 regions, for telemetry and rendering.
#[derive(Debug, Clone, Copy, Default)]
pub struct ThermalPhaseNode {
    pub temperature_k: f64,
    pub enthalpy_j_kg: f64,
    pub melt_fraction: f64,
    pub plastic_slump: f64,
    pub molten_gold_g: f64,
    pub frozen_gold_g: f64,
    pub ln2_film_g: f64,
    pub frost_thickness_um: f64,
    pub oil_film_um: f64,
    pub incandescence_rgb: [f64; 3],
}

/// The bench tools, as the user sets them.
#[derive(Debug, Clone, Copy)]
pub struct PerturbationInputs {
    pub blowtorch_power: f64,
    pub blowtorch_pos_mm: [f64; 3],
    pub ln2_pour_rate: f64,
    pub ln2_pos_mm: [f64; 3],
    pub oil_pour_rate: f64,
    pub oil_pos_mm: [f64; 3],
    pub magnet_enabled: bool,
    pub magnet_pos_mm: [f64; 3],
    pub magnet_strength: f64,
    pub compensated_hairspring: bool,
    pub air_jet_power: f64,
    pub air_jet_pos_mm: [f64; 3],
    pub ambient_pressure_atm: f64,
    /// Gravity / centrifuge acceleration in g, watch frame.
    pub gravity_g_vec: [f64; 3],
}

impl Default for PerturbationInputs {
    fn default() -> Self {
        Self {
            blowtorch_power: 0.0,
            blowtorch_pos_mm: [0.0, -24.0, 2.0],
            ln2_pour_rate: 0.0,
            ln2_pos_mm: [0.0, -16.0, 4.0],
            oil_pour_rate: 0.0,
            oil_pos_mm: [0.0, -11.0, 2.0],
            magnet_enabled: false,
            magnet_pos_mm: [65.0, -11.0, 5.0],
            magnet_strength: 1.0,
            compensated_hairspring: true,
            air_jet_power: 0.0,
            air_jet_pos_mm: [0.0, -20.0, 3.0],
            ambient_pressure_atm: 1.0,
            gravity_g_vec: [0.0, -1.0, 0.0],
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CellKind {
    Outside,
    Band,
    Interior,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Region {
    Outside,
    CaseTop,
    CaseLeft,
    CaseRight,
    CaseBottom,
    Plate,
    Barrel,
    Train,
    Cage,
    Escapement,
    Balance,
}

/// One column of the watch on the plan grid.
#[derive(Debug, Clone, Copy)]
pub struct SpatialCell {
    pub x_mm: f64,
    pub y_mm: f64,
    pub r_mm: f64,
    pub kind: CellKind,
    pub region: Region,
    /// Surface liquids rest on: the dial (0) or the tourbillon well floor [mm].
    pub bed_mm: f64,
    pub in_hairspring: bool,
    pub in_crown: bool,
    /// Structural case gold still standing in the column, and as built [kg].
    pub struct_gold_kg: f64,
    pub struct_gold0_kg: f64,
    /// Heat capacity of the non-gold layers at 293 K [J/K].
    c_other: f64,
    /// Thickness of brass for conduction [m].
    brass_m: f64,
    /// Column enthalpy relative to 293.15 K [J].
    pub enthalpy_j: f64,
    pub temperature_k: f64,
    /// Liquid fraction of the gold in the column.
    pub melt_fraction: f64,
    /// Gold that has melted out of the case and run here [kg].
    pub gold_kg: f64,
    /// Poured oil on the surface [m³].
    pub oil_m3: f64,
    pub ln2_kg: f64,
    pub frost_kg_m2: f64,
    /// Copper-oxide film on the gold surface [nm].
    pub oxide_nm: f64,
}

impl SpatialCell {
    #[must_use]
    pub fn gold_depth_m(&self) -> f64 {
        self.gold_kg / (library().gold.rho_solid_ref * 0.95 * CELL_AREA)
    }

    #[must_use]
    pub fn oil_depth_m(&self) -> f64 {
        self.oil_m3 / CELL_AREA
    }

    /// Fraction of the deposited gold that is mobile (above coherency).
    #[must_use]
    pub fn mobile_fraction(&self) -> f64 {
        ((self.melt_fraction - 0.3) / 0.3).clamp(0.0, 1.0)
    }

    #[must_use]
    pub fn liquid_gold_mm(&self) -> f64 {
        self.gold_depth_m() * self.mobile_fraction() * 1.0e3
    }

    #[must_use]
    pub fn frozen_gold_mm(&self) -> f64 {
        self.gold_depth_m() * (1.0 - self.mobile_fraction()) * 1.0e3
    }

    #[must_use]
    pub fn frost_um(&self) -> f64 {
        self.frost_kg_m2 / FROST_DENSITY * 1.0e6
    }

    /// Fraction of the structural case gold that has melted out.
    #[must_use]
    pub fn case_slump_frac(&self) -> f64 {
        if self.struct_gold0_kg > 0.0 { 1.0 - self.struct_gold_kg / self.struct_gold0_kg } else { 0.0 }
    }
}

/// One node of the bassine case band mesh.
#[derive(Debug, Clone, Copy)]
pub struct CaseShellNode {
    pub theta_rad: f64,
    pub ref_pos_mm: [f64; 3],
    pub pos_mm: [f64; 3],
    pub ref_thickness_mm: f64,
    pub thickness_mm: f64,
    pub temperature_k: f64,
    pub melt_fraction: f64,
    /// Accumulated creep sag along gravity [mm].
    pub sag_mm: f64,
    cell: usize,
}

/// Merged (T, h_gold, h_other) knots for inverting column enthalpy.
struct ColumnTable {
    t: Vec<f64>,
    hg: Vec<f64>,
    /// Debye-shaped brass enthalpy per J/K of room-temperature capacity [K].
    ho: Vec<f64>,
}

fn column_table() -> &'static ColumnTable {
    static TABLE: OnceLock<ColumnTable> = OnceLock::new();
    TABLE.get_or_init(|| {
        let gold = &library().gold;
        let mut t: Vec<f64> = gold.t_to_h().iter().map(|p| p.0).collect();
        let mut x = 40.0;
        while x < 3500.0 {
            t.push(x);
            x += if x < 400.0 { 10.0 } else { 50.0 };
        }
        // The room-temperature reference must be a node, so that a column at
        // zero enthalpy inverts to exactly 293.15 K and an idle watch is still.
        t.push(T_REF_K);
        t.sort_by(f64::total_cmp);
        t.dedup_by(|a, b| (*a - *b).abs() < 1e-6);
        let hg: Vec<f64> = t.iter().map(|&tk| gold.enthalpy_at(tk)).collect();
        // ∫ D(θ/T)/D(θ/293) dT from 293 K.
        let shape = |tk: f64| debye_fraction(BRASS_DEBYE_K / tk.max(10.0)) / debye_fraction(BRASS_DEBYE_K / T_REF_K);
        let integrate = |a: f64, b: f64| {
            let n = ((b - a).abs() / 2.0).ceil().max(1.0) as usize;
            let h = (b - a) / n as f64;
            (0..n).map(|i| 0.5 * (shape(a + i as f64 * h) + shape(a + (i as f64 + 1.0) * h)) * h).sum::<f64>()
        };
        // One cumulative pass outward from the reference node (which is in
        // the table), so the cost is linear in the temperature range.
        let k_ref = t.iter().position(|&x| x == T_REF_K).unwrap_or(0);
        let mut ho = vec![0.0; t.len()];
        for k in k_ref + 1..t.len() {
            ho[k] = ho[k - 1] + integrate(t[k - 1], t[k]);
        }
        for k in (0..k_ref).rev() {
            ho[k] = ho[k + 1] + integrate(t[k + 1], t[k]);
        }
        ColumnTable { t, hg, ho }
    })
}

impl ColumnTable {
    fn enthalpy(&self, m_gold: f64, c_other: f64, t_k: f64) -> f64 {
        let i = self.t.partition_point(|&x| x <= t_k).clamp(1, self.t.len() - 1);
        let f = ((t_k - self.t[i - 1]) / (self.t[i] - self.t[i - 1])).clamp(0.0, 1.0);
        let hg = self.hg[i - 1] + f * (self.hg[i] - self.hg[i - 1]);
        let ho = self.ho[i - 1] + f * (self.ho[i] - self.ho[i - 1]);
        m_gold * hg + c_other * ho
    }

    fn temperature(&self, m_gold: f64, c_other: f64, h: f64) -> f64 {
        let at = |k: usize| m_gold * self.hg[k] + c_other * self.ho[k];
        let n = self.t.len();
        if h <= at(0) {
            return self.t[0];
        }
        if h >= at(n - 1) {
            return self.t[n - 1];
        }
        let (mut lo, mut hi) = (0usize, n - 1);
        while hi - lo > 1 {
            let mid = (lo + hi) / 2;
            if at(mid) <= h {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let (h0, h1) = (at(lo), at(hi));
        if h1 <= h0 {
            return self.t[hi];
        }
        self.t[lo] + (self.t[hi] - self.t[lo]) * (h - h0) / (h1 - h0)
    }

    /// [`Self::temperature`] by hunting from the segment `seg` the column last
    /// inverted in (updated). The column enthalpy is non-decreasing over the
    /// knots, so the bracketing segment is unique and the result is identical
    /// to the bisection's; between substeps a column rarely leaves its segment.
    fn temperature_hunt(&self, m_gold: f64, c_other: f64, h: f64, seg: &mut u32) -> f64 {
        let at = |k: usize| m_gold * self.hg[k] + c_other * self.ho[k];
        let n = self.t.len();
        if h <= at(0) {
            *seg = 0;
            return self.t[0];
        }
        if h >= at(n - 1) {
            *seg = (n - 2) as u32;
            return self.t[n - 1];
        }
        let mut lo = (*seg as usize).min(n - 2);
        while lo > 0 && at(lo) > h {
            lo -= 1;
        }
        while at(lo + 1) <= h {
            lo += 1;
        }
        *seg = lo as u32;
        let (h0, h1) = (at(lo), at(lo + 1));
        if h1 <= h0 {
            return self.t[lo + 1];
        }
        self.t[lo] + (self.t[lo + 1] - self.t[lo]) * (h - h0) / (h1 - h0)
    }

    /// Smallest dH/dT over segment `seg` and its two neighbours [J/K]: a
    /// lower bound on the column's capacity for the next small change in T,
    /// which ignores the latent heat of a mushy segment the column may be
    /// about to leave.
    fn capacity_floor(&self, m_gold: f64, c_other: f64, seg: u32) -> f64 {
        let n = self.t.len();
        let s = (seg as usize).min(n - 2);
        let (a, b) = (s.saturating_sub(1), (s + 1).min(n - 2));
        (a..=b)
            .map(|j| (m_gold * (self.hg[j + 1] - self.hg[j]) + c_other * (self.ho[j + 1] - self.ho[j])) / (self.t[j + 1] - self.t[j]))
            .fold(f64::INFINITY, f64::min)
    }

    /// Effective heat capacity dH/dT at `t_k` [J/K].
    #[cfg(test)]
    fn capacity(&self, m_gold: f64, c_other: f64, t_k: f64) -> f64 {
        let i = self.t.partition_point(|&x| x <= t_k).clamp(1, self.t.len() - 1);
        let dt = self.t[i] - self.t[i - 1];
        (m_gold * (self.hg[i] - self.hg[i - 1]) + c_other * (self.ho[i] - self.ho[i - 1])) / dt
    }
}

/// Liquid-nitrogen pool-boiling heat flux [W/m²] at wall superheat `dt_sat`
/// [K], pressure `p_pa` and gravity `g`.
#[must_use]
pub fn ln2_boiling_flux(dt_sat: f64, p_pa: f64, g: f64) -> f64 {
    if dt_sat <= 0.0 {
        return 0.0;
    }
    let g = g.max(0.1 * G0_M_S2);
    let t_sat = ln2_saturation_temperature_k(p_pa);
    let (rl, mul, cpl, kl, sigma, hfg) =
        (ln2::RHO_LIQUID, ln2::MU_LIQUID, ln2::CP_LIQUID, ln2::K_LIQUID, ln2::SURFACE_TENSION, ln2::LATENT_J_KG);
    // Nitrogen vapour: ideal-gas density, kinetic-theory μ, k ∝ T^0.8.
    let vapour = |t: f64| {
        let rho = p_pa.max(ln2::P_TRIPLE_PA) * ln2::MOLAR_MASS / (R_GAS * t);
        let mu = 5.4e-6 * (t / 77.0).powf(0.75);
        let k = 7.6e-3 * (t / 77.0).powf(0.85);
        (rho, mu, k)
    };
    let (rv_sat, _, _) = vapour(t_sat);
    let drho = rl - rv_sat;
    let pr_l = mul * cpl / kl;
    // Rohsenow nucleate boiling (C_sf = 0.013, n = 1.7).
    let rohsenow = |dt: f64| {
        mul * hfg * (g * drho / sigma).sqrt() * (cpl * dt / (0.013 * hfg * pr_l.powf(1.7))).powi(3)
    };
    // Zuber critical heat flux (Lienhard coefficient 0.149).
    let q_max = 0.149 * hfg * rv_sat.sqrt() * (sigma * g * drho).powf(0.25);
    let dt_chf = (q_max / rohsenow(1.0)).cbrt();
    // Zuber minimum (Leidenfrost) flux and Berenson's superheat for it.
    let q_min = 0.09 * rv_sat * hfg * (sigma * g * drho / (rl + rv_sat).powi(2)).powf(0.25);
    let (rvf, muvf, kvf) = vapour(t_sat + 15.0);
    let dt_min = 0.127 * rvf * hfg / kvf * (g * drho / (rl + rvf)).powf(2.0 / 3.0) * (sigma / (g * drho)).sqrt()
        * (muvf / (g * drho)).cbrt();
    let berenson = |dt: f64| {
        let tf = t_sat + 0.5 * dt;
        let (rv, muv, kv) = vapour(tf);
        let lam = (sigma / (g * drho)).sqrt();
        let hfg_eff = hfg + 0.4 * 1_040.0 * dt;
        let h = 0.425 * (kv.powi(3) * rv * g * (rl - rv) * hfg_eff / (muv * dt * lam)).powf(0.25);
        // Radiation across the vapour film (ε ≈ 0.3 metal, 1 liquid).
        let ts = t_sat + dt;
        h * dt + 0.3 * STEFAN_BOLTZMANN * (ts.powi(4) - t_sat.powi(4))
    };
    if dt_sat <= dt_chf {
        rohsenow(dt_sat)
    } else if dt_sat >= dt_min.max(dt_chf * 1.5) {
        berenson(dt_sat).max(q_min)
    } else {
        // Transition boiling: log-linear between the two extrema.
        let d_min = dt_min.max(dt_chf * 1.5);
        let f = (dt_sat.ln() - dt_chf.ln()) / (d_min.ln() - dt_chf.ln());
        (q_max.ln() + f * (berenson(d_min).max(q_min).ln() - q_max.ln())).exp()
    }
}

/// The LN₂ boiling curve for one (pressure, gravity), tabulated on log ΔT:
/// evaluating the correlations per wet cell per substep cost more than the
/// rest of the field solve.
#[derive(Debug, Clone)]
struct BoilingTable {
    p_pa: f64,
    g: f64,
    q: Vec<f64>,
    /// Steepest slope of the curve [W/(m² K)].
    max_slope: f64,
}

impl BoilingTable {
    const LN_LO: f64 = -2.302_585_092_994_046; // ln 0.1 K
    const LN_HI: f64 = 8.006_367_567_650_246; // ln 3000 K
    const N: usize = 256;

    fn build(p_pa: f64, g: f64) -> Self {
        let q: Vec<f64> = (0..Self::N)
            .map(|k| {
                let x = Self::LN_LO + (Self::LN_HI - Self::LN_LO) * k as f64 / (Self::N - 1) as f64;
                ln2_boiling_flux(x.exp(), p_pa, g)
            })
            .collect();
        // Steepest |dq/dΔT| of the curve (the 0.1 K linear ramp included):
        // the stiffest a wetted cell's boiling can be.
        let mut max_slope = q[0] / 0.1;
        for k in 0..Self::N - 1 {
            let x0 = (Self::LN_LO + (Self::LN_HI - Self::LN_LO) * k as f64 / (Self::N - 1) as f64).exp();
            let x1 = (Self::LN_LO + (Self::LN_HI - Self::LN_LO) * (k + 1) as f64 / (Self::N - 1) as f64).exp();
            max_slope = max_slope.max(((q[k + 1] - q[k]) / (x1 - x0)).abs());
        }
        Self { p_pa, g, q, max_slope }
    }

    fn flux(&self, dt_sat: f64) -> f64 {
        if dt_sat <= 0.1 {
            return if dt_sat > 0.0 { self.q[0] * dt_sat / 0.1 } else { 0.0 };
        }
        let f = (dt_sat.ln() - Self::LN_LO) / (Self::LN_HI - Self::LN_LO) * (Self::N - 1) as f64;
        let i = (f.floor() as usize).min(Self::N - 2);
        let w = (f - i as f64).clamp(0.0, 1.0);
        self.q[i] + w * (self.q[i + 1] - self.q[i])
    }
}

/// Natural-convection film coefficient [W/(m² K)] for the whole watch, from
/// `fs-convection` (Churchill–Chu vertical plate, or the horizontal-plate
/// correlations when it lies flat), tabulated against ΔT at the current
/// pressure.
#[derive(Debug, Clone)]
struct ConvectionTable {
    pressure_pa: f64,
    vertical: bool,
    dts: [f64; 9],
    h: [f64; 9],
}

impl ConvectionTable {
    fn build(pressure_pa: f64, vertical: bool, g: f64) -> Self {
        let dts = [0.5, 2.0, 8.0, 30.0, 80.0, 200.0, 400.0, 700.0, 1100.0];
        let mut h = [0.0; 9];
        let l = 2.0 * CASE_OUTER_R_MM * 1.0e-3;
        for (slot, &dt) in h.iter_mut().zip(&dts) {
            let tf = T_REF_K + 0.5 * dt;
            let Some(air) = air_state(tf, pressure_pa) else { continue };
            let nu = air.dynamic_viscosity / air.density;
            let alpha = air.thermal_conductivity / (air.density * air.specific_heat_cp);
            let (id, length) = if vertical {
                (CorrelationId::ChurchillChuVerticalPlate, l)
            } else {
                (CorrelationId::HorizontalPlateUnstableLaminar, l / 4.0)
            };
            let ra = (g * dt / tf * length.powi(3) / (nu * alpha)).max(1.0e-3);
            let k = air.thermal_conductivity;
            let correlated = eval_convection(id, CorrelationInputs::natural(ra, air.prandtl))
                .and_then(|nu_eval| {
                    nu_eval.heat_transfer_coefficient(ThermalConductivity::new(k), ConvLength::new(length))
                })
                .map(|ev| ev.value.value());
            // Outside a correlation's validity window fall back to its
            // laminar form (Nu = 0.68 + 0.67 Ra^¼ / [1 + (0.492/Pr)^(9/16)]^(4/9)).
            *slot = correlated.unwrap_or_else(|_| {
                let pr = air.prandtl;
                let nu_l = 0.68 + 0.670 * ra.powf(0.25) / (1.0 + (0.492 / pr).powf(9.0 / 16.0)).powf(4.0 / 9.0);
                nu_l * k / length
            });
        }
        Self { pressure_pa, vertical, dts, h }
    }

    fn eval(&self, dt: f64) -> f64 {
        let d = dt.abs();
        if d <= self.dts[0] {
            return self.h[0];
        }
        for i in 1..self.dts.len() {
            if d <= self.dts[i] {
                let f = (d - self.dts[i - 1]) / (self.dts[i] - self.dts[i - 1]);
                return self.h[i - 1] + f * (self.h[i] - self.h[i - 1]);
            }
        }
        self.h[self.h.len() - 1]
    }
}

/// Reused per-cell buffers of the virtual-pipe solver (one substep would
/// otherwise allocate and zero seven grid-sized arrays).
#[derive(Clone)]
struct PipeScratch {
    depth: Vec<f64>,
    head: Vec<f64>,
    nu: Vec<f64>,
    out_v: Vec<f64>,
    in_v: Vec<f64>,
    s_out: Vec<f64>,
    s_in: Vec<f64>,
    moves: Vec<(usize, usize, f64)>,
}

impl Default for PipeScratch {
    fn default() -> Self {
        let z = || vec![0.0; GRID_CELLS];
        Self { depth: z(), head: z(), nu: z(), out_v: z(), in_v: z(), s_out: z(), s_in: z(), moves: Vec::with_capacity(512) }
    }
}

/// Reused per-cell buffers of the explicit heat solve, and the parts of it
/// that hold still within one step (the nozzles and the oxide film do not
/// move between substeps).
#[derive(Clone, Default)]
struct ThermalScratch {
    /// Conductance thickness k·t of each column [W/K].
    kt: Vec<f64>,
    dq: Vec<f64>,
    /// Torch and air-jet film coefficients [W/(m² K)].
    torch_h: Vec<f64>,
    jet_h: Vec<f64>,
    /// Oxide part of the gold emissivity.
    eps_ox: Vec<f64>,
    /// Enthalpy-table segment each column last inverted in.
    seg: Vec<u32>,
    /// Columns inside the watch.
    active: Vec<u32>,
}

impl ThermalScratch {
    fn new(grid: &[SpatialCell]) -> Self {
        let table = column_table();
        let z = || vec![0.0; GRID_CELLS];
        let mut seg = vec![0u32; GRID_CELLS];
        for (c, s) in grid.iter().zip(&mut seg) {
            table.temperature_hunt(c.struct_gold_kg + c.gold_kg, c.c_other.max(1.0e-5), c.enthalpy_j, s);
        }
        let active = (0..GRID_CELLS).filter(|&i| grid[i].kind != CellKind::Outside).map(|i| i as u32).collect();
        Self { kt: z(), dq: z(), torch_h: z(), jet_h: z(), eps_ox: z(), seg, active }
    }
}

/// Longest explicit thermal substep [s]; shorter when the stability bound
/// of the current fields is.
const H_THERMAL_MAX_S: f64 = 10.0e-3;

/// Complete field state of the watch.
#[derive(Clone)]
pub struct MultiphysicsState {
    pub grid: Vec<SpatialCell>,
    pub case_shell: Vec<CaseShellNode>,
    pub nodes: [ThermalPhaseNode; NODE_COUNT],
    pub inputs: PerturbationInputs,
    /// Face fluxes of the virtual-pipe layers [m³/s]: x faces (i,i+1) and y
    /// faces (i,i+N), indexed by the lower cell.
    gold_qx: Vec<f64>,
    gold_qy: Vec<f64>,
    oil_qx: Vec<f64>,
    oil_qy: Vec<f64>,
    /// Torch + air-jet dynamic pressure on each cell [Pa], refreshed once
    /// per flow call (the nozzles do not move within it).
    gas_q: Vec<f64>,
    pipe: PipeScratch,
    th: ThermalScratch,
    part_cells: [Vec<usize>; 5],
    crown_cells: Vec<usize>,
    convection: ConvectionTable,
    boiling: BoilingTable,
    /// Hairspring magnetisation M/Mₛ (play-operator state).
    pub hairspring_magnetisation: f64,
    degauss_t: f64,
    pub shock_envelope_g: f64,
    pub hairspring_oxide_nm: f64,
    // Outputs.
    pub environment: MovementEnvironment,
    pub b_field_escapement_mt: f64,
    pub remanent_magnetization_mt: f64,
    pub ln2_boiling_power_w: f64,
    pub blowtorch_heating_power_w: f64,
    pub air_cooling_power_w: f64,
    pub dripped_gold_g: f64,
    pub blued_steel_rgb: [f64; 3],
    quiet: bool,
}

/// Parts sampled for the movement: barrel, train (cage frame), cage-well
/// escapement, balance, and the third wheel under the dial.
const PART_BARREL: usize = 0;
const PART_TRAIN: usize = 1;
const PART_ESCAPEMENT: usize = 2;
const PART_BALANCE: usize = 3;
const PART_PLATE: usize = 4;

impl MultiphysicsState {
    #[must_use]
    pub fn new() -> Self {
        let lib = library();
        let rho_gold = lib.gold.rho_solid_ref;
        let brass = lib.brass;
        // Bassine case band profile (radius, height, wall) as rendered.
        let profile_z_mm: [f64; CASE_SHELL_N_Z] = [-3.8, -2.8, -1.6, -0.5, 0.5, 1.6, 2.8, 3.8];
        let profile_r_mm: [f64; CASE_SHELL_N_Z] = [23.0, 25.0, 26.3, 26.8, 26.8, 26.3, 25.0, 23.0];
        let profile_h_mm: [f64; CASE_SHELL_N_Z] = [2.1, 2.3, 2.6, 2.8, 2.8, 2.6, 2.3, 2.1];
        // Gold volume of the band per unit plan area: meridional length × wall
        // thickness, spread over the band's radial extent.
        let mut merid_wall = 0.0;
        for iz in 0..CASE_SHELL_N_Z - 1 {
            let dr = profile_r_mm[iz + 1] - profile_r_mm[iz];
            let dz = profile_z_mm[iz + 1] - profile_z_mm[iz];
            merid_wall += (dr * dr + dz * dz).sqrt() * 0.5 * (profile_h_mm[iz] + profile_h_mm[iz + 1]);
        }
        let band_gold_thick_m = merid_wall / (CASE_OUTER_R_MM - CASE_INNER_R_MM) * 1.0e-3 * 0.55;

        let mut grid = Vec::with_capacity(GRID_CELLS);
        for idx in 0..GRID_CELLS {
            let gx = idx % GRID_N;
            let gy = idx / GRID_N;
            let x = -GRID_HALF_SPAN_MM + (gx as f64 + 0.5) * GRID_DX_MM;
            let y = -GRID_HALF_SPAN_MM + (gy as f64 + 0.5) * GRID_DX_MM;
            let r = x.hypot(y);
            let in_pendant = x.abs() <= 4.8 && y > 22.0 && y <= 28.8;
            let kind = if r <= CASE_INNER_R_MM && !in_pendant {
                CellKind::Interior
            } else if r <= CASE_OUTER_R_MM || in_pendant {
                CellKind::Band
            } else {
                CellKind::Outside
            };
            let dt_ = (x - TOURBILLON_CENTER_MM[0]).hypot(y - TOURBILLON_CENTER_MM[1]);
            let db = (x - BARREL_CENTER_MM[0]).hypot(y - BARREL_CENTER_MM[1]);
            let ang = y.atan2(x);
            let region = match kind {
                CellKind::Outside => Region::Outside,
                CellKind::Band => {
                    if (PI / 4.0..3.0 * PI / 4.0).contains(&ang) {
                        Region::CaseTop
                    } else if (-3.0 * PI / 4.0..-PI / 4.0).contains(&ang) {
                        Region::CaseBottom
                    } else if ang.abs() >= 3.0 * PI / 4.0 {
                        Region::CaseLeft
                    } else {
                        Region::CaseRight
                    }
                }
                CellKind::Interior => {
                    if dt_ <= 5.6 {
                        Region::Balance
                    } else if dt_ <= 6.6 {
                        Region::Escapement
                    } else if dt_ <= WELL_R_MM {
                        Region::Cage
                    } else if db <= 7.6 {
                        Region::Barrel
                    } else if x.abs() <= 6.5 && y > -3.5 && y < 4.0 {
                        Region::Train
                    } else {
                        Region::Plate
                    }
                }
            };
            let in_well = kind == CellKind::Interior && dt_ <= WELL_R_MM;
            let (struct_gold, c_other, brass_m) = match kind {
                CellKind::Outside => (0.0, 0.0, 0.0),
                CellKind::Band => (rho_gold * band_gold_thick_m * CELL_AREA, 0.0, 0.0),
                CellKind::Interior => {
                    let brass_t = if in_well { 0.5 * BRASS_M } else { BRASS_M };
                    (rho_gold * BACK_GOLD_M * CELL_AREA, brass.rho * brass.cp * brass_t * CELL_AREA, brass_t)
                }
            };
            grid.push(SpatialCell {
                x_mm: x,
                y_mm: y,
                r_mm: r,
                kind,
                region,
                bed_mm: if in_well { WELL_FLOOR_MM } else { 0.0 },
                in_hairspring: kind == CellKind::Interior && dt_ <= 4.0,
                in_crown: in_pendant || (kind == CellKind::Band && x.abs() < 5.0 && y > 0.0),
                struct_gold_kg: struct_gold,
                struct_gold0_kg: struct_gold,
                c_other,
                brass_m,
                enthalpy_j: 0.0,
                temperature_k: T_REF_K,
                melt_fraction: 0.0,
                gold_kg: 0.0,
                oil_m3: 0.0,
                ln2_kg: 0.0,
                frost_kg_m2: 0.0,
                oxide_nm: 2.0,
            });
        }

        // Each column starts at room temperature on its own enthalpy curve.
        for c in &mut grid {
            if c.kind != CellKind::Outside {
                c.enthalpy_j = column_enthalpy(c, T_REF_K);
            }
        }

        let d_theta = TAU / CASE_SHELL_N_THETA as f64;
        let case_shell = (0..CASE_SHELL_NODES)
            .map(|idx| {
                let it = idx % CASE_SHELL_N_THETA;
                let iz = idx / CASE_SHELL_N_THETA;
                let theta = it as f64 * d_theta;
                let pendant_lobe = if (theta - 0.5 * PI).abs() < 0.22 {
                    let u = 1.0 - (theta - 0.5 * PI).abs() / 0.22;
                    2.2 * u * u * (1.0 - 0.12 * profile_z_mm[iz].abs())
                } else {
                    0.0
                };
                let r0 = profile_r_mm[iz] + pendant_lobe * theta.sin().max(0.0);
                let pos = [r0 * theta.cos(), r0 * theta.sin(), profile_z_mm[iz]];
                // The band column this node belongs to (mid-band radius).
                let rm = 0.5 * (CASE_INNER_R_MM + CASE_OUTER_R_MM) + 0.4 * pendant_lobe;
                let cell = cell_at(rm * theta.cos(), rm * theta.sin()).unwrap_or(0);
                CaseShellNode {
                    theta_rad: theta,
                    ref_pos_mm: pos,
                    pos_mm: pos,
                    ref_thickness_mm: profile_h_mm[iz] + 0.35 * pendant_lobe,
                    thickness_mm: profile_h_mm[iz] + 0.35 * pendant_lobe,
                    temperature_k: T_REF_K,
                    melt_fraction: 0.0,
                    sag_mm: 0.0,
                    cell,
                }
            })
            .collect();

        let mut part_cells: [Vec<usize>; 5] = Default::default();
        let mut crown_cells = Vec::new();
        for (i, c) in grid.iter().enumerate() {
            match c.region {
                Region::Barrel => part_cells[PART_BARREL].push(i),
                Region::Cage => part_cells[PART_TRAIN].push(i),
                Region::Escapement => part_cells[PART_ESCAPEMENT].push(i),
                Region::Balance => part_cells[PART_BALANCE].push(i),
                Region::Train => part_cells[PART_PLATE].push(i),
                _ => {}
            }
            if c.in_crown {
                crown_cells.push(i);
            }
        }
        let th = ThermalScratch::new(&grid);
        let mut s = Self {
            grid,
            th,
            case_shell,
            nodes: [ThermalPhaseNode::default(); NODE_COUNT],
            inputs: PerturbationInputs::default(),
            gold_qx: vec![0.0; GRID_CELLS],
            gold_qy: vec![0.0; GRID_CELLS],
            oil_qx: vec![0.0; GRID_CELLS],
            oil_qy: vec![0.0; GRID_CELLS],
            gas_q: vec![0.0; GRID_CELLS],
            pipe: PipeScratch::default(),
            part_cells,
            crown_cells,
            convection: ConvectionTable::build(P_ATM_PA, true, G0_M_S2),
            boiling: BoilingTable::build(P_ATM_PA, G0_M_S2),
            hairspring_magnetisation: 0.0,
            degauss_t: -1.0,
            shock_envelope_g: 0.0,
            hairspring_oxide_nm: 2.0,
            environment: MovementEnvironment::default(),
            b_field_escapement_mt: 0.0,
            remanent_magnetization_mt: 0.0,
            ln2_boiling_power_w: 0.0,
            blowtorch_heating_power_w: 0.0,
            air_cooling_power_w: 0.0,
            dripped_gold_g: 0.0,
            blued_steel_rgb: [0.0; 3],
            quiet: false,
        };
        s.refresh_outputs(0.0, 0.0);
        s
    }

    fn gravity(&self) -> [f64; 3] {
        self.inputs.gravity_g_vec.map(|g| g * G0_M_S2)
    }

    fn pressure_pa(&self) -> f64 {
        self.inputs.ambient_pressure_atm.max(0.0) * P_ATM_PA
    }

    /// Start an AC degauss: a 50 Hz field from 40 mT decaying geometrically
    /// (τ = 0.35 s, ~3 % per half cycle) for 3 s, so the last swings through
    /// the coercive field are small and leave almost no remanence.
    pub fn trigger_degauss(&mut self) {
        self.degauss_t = 0.0;
    }

    /// Slosh the liquids and ring the case for a bench drop.
    pub fn apply_mechanical_shock(&mut self, g_peak: f64) {
        self.shock_envelope_g = self.shock_envelope_g.max(g_peak);
        // The frame jumps by Δv; the liquid layers lag behind it.
        let dv = 2.0 / PI * g_peak * G0_M_S2 * 3.0e-4;
        for i in 0..GRID_CELLS {
            let hg = self.grid[i].gold_depth_m() * self.grid[i].mobile_fraction();
            let ho = self.grid[i].oil_depth_m();
            self.gold_qy[i] -= dv * hg * DX;
            self.oil_qy[i] -= dv * ho * DX;
        }
        self.quiet = false;
    }

    /// Ultrasonic clean: wash out all poured oil and melt off the frost.
    pub fn ultrasonic_clean_and_relubricate(&mut self) {
        for c in &mut self.grid {
            c.oil_m3 = 0.0;
        }
        self.oil_qx.fill(0.0);
        self.oil_qy.fill(0.0);
    }

    /// Re-run the field solve on the next step (inputs such as gravity or
    /// pressure changed while the grid was quiescent).
    pub fn wake(&mut self) {
        self.quiet = false;
    }

    /// Re-sample the movement environment from the current fields without
    /// advancing time (a setter changed pressure, gravity or the air jet).
    pub fn refresh_environment(&mut self) {
        self.refresh_outputs(0.0, 0.0);
    }

    fn tools_active(&self) -> bool {
        let i = &self.inputs;
        i.blowtorch_power > 0.0 || i.ln2_pour_rate > 0.0 || i.oil_pour_rate > 0.0 || i.air_jet_power > 0.0
    }

    /// Advance the fields by `dt`. `cage_angle` and `lever_angle` place the
    /// steel lever for the magnetic torque.
    pub fn step(&mut self, dt: f64, cage_angle: f64, lever_angle: f64) {
        let vertical = {
            let g = self.inputs.gravity_g_vec;
            g[0].hypot(g[1]) >= g[2].abs()
        };
        let g_mag = self.gravity().iter().map(|x| x * x).sum::<f64>().sqrt();
        if (self.convection.pressure_pa - self.pressure_pa()).abs() > 1.0 || self.convection.vertical != vertical {
            self.convection = ConvectionTable::build(self.pressure_pa(), vertical, g_mag.max(0.1));
        }
        if self.inputs.ln2_pour_rate > 0.0 || self.grid.iter().any(|c| c.ln2_kg > 0.0) {
            let (p, g) = (self.pressure_pa(), g_mag.max(0.1 * G0_M_S2));
            if (self.boiling.p_pa - p).abs() > 1e-3 * p.max(1.0) || (self.boiling.g - g).abs() > 1e-3 * g {
                self.boiling = BoilingTable::build(p, g);
            }
        }
        if self.quiet && !self.tools_active() {
            self.update_magnetics(dt, cage_angle, lever_angle);
            self.shock_envelope_g *= (-dt / 0.3).exp();
            self.refresh_outputs(cage_angle, lever_angle);
            return;
        }
        // Thermal substeps within the explicit stability limit of the
        // current fields.
        let mut th = std::mem::take(&mut self.th);
        let h_stable = self.prepare_thermal(&mut th);
        let n_th = (dt / h_stable).ceil().max(1.0) as usize;
        let h_th = dt / n_th as f64;
        let (mut torch_w, mut ln2_w, mut air_w) = (0.0, 0.0, 0.0);
        for _ in 0..n_th {
            let (a, b, c) = self.thermal_substep(h_th, &mut th);
            torch_w += a * h_th;
            ln2_w += b * h_th;
            air_w += c * h_th;
            self.melt_case();
        }
        self.th = th;
        self.blowtorch_heating_power_w = torch_w / dt;
        self.ln2_boiling_power_w = ln2_w / dt;
        self.air_cooling_power_w = air_w / dt;
        // Copper oxide on hot gold (parabolic growth; slow, so once per step).
        for c in &mut self.grid {
            if c.kind != CellKind::Outside && c.temperature_k > 420.0 {
                let k_ox = 0.3 * steel_temper_oxide_rate_nm2_s(c.temperature_k);
                if k_ox > 0.0 {
                    c.oxide_nm = (c.oxide_nm * c.oxide_nm + k_ox * dt).sqrt();
                }
            }
        }
        self.flow_liquids(dt);
        self.flow_ln2(dt);
        self.update_magnetics(dt, cage_angle, lever_angle);
        self.update_case_shell(dt);
        self.shock_envelope_g *= (-dt / 0.3).exp();
        // Hairspring temper colour.
        let tb = self.environment.balance.temperature_k;
        let x = self.hairspring_oxide_nm;
        self.hairspring_oxide_nm = (x * x + steel_temper_oxide_rate_nm2_s(tb) * dt).sqrt();
        self.refresh_outputs(cage_angle, lever_angle);
        // Quiescent once nothing is happening: skip the grid until a tool acts.
        let calm = self.grid.iter().all(|c| {
            c.kind == CellKind::Outside
                || ((c.temperature_k - T_REF_K).abs() < 0.05 && c.ln2_kg <= 0.0 && c.frost_kg_m2 <= 0.0)
        });
        let still = self.gold_qx.iter().chain(&self.gold_qy).chain(&self.oil_qx).chain(&self.oil_qy).all(|q| q.abs() < 1e-14);
        self.quiet = calm && still && !self.tools_active();
    }

    /// Torch and jet footprints: film coefficient [W/(m² K)] and stagnation
    /// pressure [Pa] at a cell.
    fn torch_at(&self, c: &SpatialCell) -> (f64, f64) {
        let i = &self.inputs;
        if i.blowtorch_power <= 0.0 {
            return (0.0, 0.0);
        }
        let p = i.blowtorch_pos_mm;
        // Air-fuel torch: no flame without the room's oxygen; denser jets at
        // higher pressure carry more heat.
        let p_atm = i.ambient_pressure_atm;
        let o2 = (p_atm / 0.5).clamp(0.0, 1.0) * p_atm.max(1.0).powf(0.3);
        let standoff = (p[2] - 0.0).abs() + 4.0;
        let sigma = 1.8 + 0.25 * standoff;
        let r2 = (c.x_mm - p[0]).powi(2) + (c.y_mm - p[1]).powi(2);
        if r2 > 25.0 * sigma * sigma {
            return (0.0, 0.0);
        }
        let shape = (-r2 / (2.0 * sigma * sigma)).exp();
        let h = 4_000.0 * i.blowtorch_power * o2 * (6.0 / standoff).sqrt().min(1.3) * shape;
        let q_dyn = 220.0 * i.blowtorch_power * o2 * shape;
        (h, q_dyn)
    }

    fn jet_at(&self, c: &SpatialCell) -> (f64, f64) {
        let i = &self.inputs;
        if i.air_jet_power <= 0.0 {
            return (0.0, 0.0);
        }
        let p = i.air_jet_pos_mm;
        let d_noz = 1.5;
        let standoff = p[2].abs().max(1.0) + 2.0;
        let u0 = 150.0 * i.air_jet_power;
        let u = u0 * (6.2 * d_noz / standoff).min(1.0);
        let sigma = 1.0 + 0.2 * standoff;
        let r2 = (c.x_mm - p[0]).powi(2) + (c.y_mm - p[1]).powi(2);
        if r2 > 25.0 * sigma * sigma {
            return (0.0, 0.0);
        }
        let shape = (-r2 / (2.0 * sigma * sigma)).exp();
        let rho = 1.2;
        let nu = 1.5e-5;
        let re = u * d_noz * 1.0e-3 / nu;
        let h0 = 0.9 * re.sqrt() * 0.71_f64.powf(0.4) * 0.026 / (d_noz * 1.0e-3);
        (h0 * shape, 0.5 * rho * u * u * shape)
    }

    /// Conductance thickness k·t of each column at its temperature [W/K].
    fn fill_conductance(&self, th: &mut ThermalScratch) {
        let lib = library();
        let gold = &lib.gold;
        for &iu in &th.active {
            let i = iu as usize;
            let c = &self.grid[i];
            let tg = (c.struct_gold_kg + c.gold_kg) / (gold.rho_solid_ref * CELL_AREA);
            th.kt[i] = gold.conductivity(c.temperature_k, c.melt_fraction) * tg + lib.brass.k * c.brass_m;
        }
    }

    /// Fill the per-step parts of the heat solve (nozzle footprints, oxide
    /// emissivity, table segments) and return the explicit stability bound
    /// [s]: forward Euler on C dT/dt = −D T + … is stable for h < C/D in
    /// every column (Gershgorin), with D the column's face conductances plus
    /// the linearised film, radiation and boiling coefficients and C a floor
    /// on its capacity (sensible heat only, so a column leaving the mushy
    /// zone mid-step stays inside the bound).
    fn prepare_thermal(&self, th: &mut ThermalScratch) -> f64 {
        let table = column_table();
        let n = GRID_N;
        let a = CELL_AREA;
        let torch_on = self.inputs.blowtorch_power > 0.0;
        let jet_on = self.inputs.air_jet_power > 0.0;
        self.fill_conductance(th);
        let mut h_min = H_THERMAL_MAX_S;
        for &iu in &th.active {
            let i = iu as usize;
            let c = &self.grid[i];
            th.torch_h[i] = if torch_on { self.torch_at(c).0 } else { 0.0 };
            th.jet_h[i] = if jet_on { self.jet_at(c).0 } else { 0.0 };
            th.eps_ox[i] = GoldAlloy::oxide_emissivity(c.oxide_nm);
            // Gold may have run in or out since the last substep.
            let m_gold = c.struct_gold_kg + c.gold_kg;
            let c_other = c.c_other.max(1.0e-5);
            table.temperature_hunt(m_gold, c_other, c.enthalpy_j, &mut th.seg[i]);
            let (gx, gy) = (i % n, i / n);
            let k0 = th.kt[i];
            let face = |j: usize| {
                let kj = th.kt[j];
                if kj > 0.0 {
                    2.0 * k0 * kj / (k0 + kj)
                } else {
                    0.0
                }
            };
            let mut d = 0.0;
            if gx > 0 {
                d += face(i - 1);
            }
            if gx + 1 < n {
                d += face(i + 1);
            }
            if gy > 0 {
                d += face(i - n);
            }
            if gy + 1 < n {
                d += face(i + n);
            }
            let t = c.temperature_k;
            let h_nat = self.convection.eval(t - T_REF_K);
            d += (th.torch_h[i] + th.jet_h[i] + 2.0 * 1.5 * h_nat) * a;
            d += 4.0 * 1.3 * STEFAN_BOLTZMANN * t.powi(3) * a;
            if c.ln2_kg > 0.0 {
                d += self.boiling.max_slope * a;
            }
            let cap = table.capacity_floor(m_gold, c_other, th.seg[i]);
            if d > 0.0 {
                h_min = h_min.min(cap / d);
            }
        }
        (0.9 * h_min).max(2.0e-5)
    }

    /// One explicit thermal step. Returns (torch, LN₂, air-jet) powers [W].
    fn thermal_substep(&mut self, h: f64, th: &mut ThermalScratch) -> (f64, f64, f64) {
        let lib = library();
        let table = column_table();
        let gold = &lib.gold;
        let n = GRID_N;
        let t_amb = T_REF_K;
        let p_pa = self.pressure_pa();
        let t_sat = ln2_saturation_temperature_k(p_pa);
        let flame_t = 2_250.0;
        let rho_v_amb = ROOM_RH * WATER_P_SAT_293_PA * (p_pa / P_ATM_PA) * M_WATER / (R_GAS * t_amb);
        let air = air_state(t_amb, p_pa);
        let rho_cp_air = air.map_or(0.0, |a| a.density * a.specific_heat_cp);
        self.fill_conductance(th);
        let kt = &th.kt;
        let dq = &mut th.dq;
        dq.fill(0.0);
        for gy in 0..n {
            for gx in 0..n {
                let i = gy * n + gx;
                if kt[i] <= 0.0 {
                    continue;
                }
                if gx + 1 < n && kt[i + 1] > 0.0 {
                    let gface = 2.0 * kt[i] * kt[i + 1] / (kt[i] + kt[i + 1]);
                    let q = gface * (self.grid[i + 1].temperature_k - self.grid[i].temperature_k);
                    dq[i] += q;
                    dq[i + 1] -= q;
                }
                if gy + 1 < n && kt[i + n] > 0.0 {
                    let gface = 2.0 * kt[i] * kt[i + n] / (kt[i] + kt[i + n]);
                    let q = gface * (self.grid[i + n].temperature_k - self.grid[i].temperature_k);
                    dq[i] += q;
                    dq[i + n] -= q;
                }
            }
        }
        let (mut torch_w, mut ln2_w, mut air_w) = (0.0, 0.0, 0.0);
        for &iu in &th.active {
            let i = iu as usize;
            let (h_torch, h_jet) = (th.torch_h[i], th.jet_h[i]);
            let c = &mut self.grid[i];
            let t = c.temperature_k;
            let a = CELL_AREA;
            let mut q = dq[i];
            if h_torch > 0.0 {
                let p = h_torch * (flame_t - t) * a;
                q += p;
                torch_w += p;
            }
            if h_jet > 0.0 {
                let p = h_jet * (t_amb - t) * a;
                q += p;
                air_w += p;
            }
            // LN₂ boiling under the covered part of the cell.
            let cover = (c.ln2_kg / (ln2::RHO_LIQUID * a * 0.2e-3)).clamp(0.0, 1.0);
            if c.ln2_kg > 0.0 {
                let flux = self.boiling.flux(t - t_sat);
                let mut p = flux * a * cover;
                // Each kilogram boiled also carries off the sensible heat of
                // its vapour, which leaves near the film temperature (the
                // same modified latent heat the Berenson flux uses).
                let hfg_eff = ln2::LATENT_J_KG + 0.4 * 1_040.0 * (t - t_sat).max(0.0);
                let boil_max = c.ln2_kg * hfg_eff / h;
                p = p.min(boil_max);
                c.ln2_kg -= p * h / hfg_eff;
                if c.ln2_kg < 1e-12 {
                    c.ln2_kg = 0.0;
                }
                q -= p;
                ln2_w += p;
            }
            // Natural convection on both faces and grey-body radiation.
            let gold_frac = if c.kind == CellKind::Band { 1.0 } else { 0.0 };
            let eps_gold = gold.emissivity_with_oxide(c.melt_fraction, th.eps_ox[i]);
            let eps = gold_frac * 2.0 * eps_gold + (1.0 - gold_frac) * (0.35 + eps_gold);
            let h_nat = self.convection.eval(t - t_amb);
            q -= h_nat * (t - t_amb) * 2.0 * a * (1.0 - 0.5 * cover);
            q -= eps * STEFAN_BOLTZMANN * (t.powi(4) - t_amb.powi(4)) * a;
            // Frost: vapour deposits below the ice saturation point.
            if rho_cp_air > 0.0 {
                let p_ice = if t < 273.15 { ice_saturation_pressure_pa(t) } else { f64::INFINITY };
                let rho_v_s = if p_ice.is_finite() { p_ice * M_WATER / (R_GAS * t.max(1.0)) } else { f64::INFINITY };
                let h_m = h_nat / rho_cp_air * 1.1;
                if rho_v_s < rho_v_amb {
                    let m = h_m * (rho_v_amb - rho_v_s) * (1.0 - cover);
                    c.frost_kg_m2 += m * h;
                    q += m * ICE_LATENT_SUB_J_KG * a;
                } else if c.frost_kg_m2 > 0.0 {
                    if t >= 273.15 {
                        let melt = (c.frost_kg_m2 * 0.2).max(1e-5 * h);
                        let dm = melt.min(c.frost_kg_m2);
                        c.frost_kg_m2 -= dm;
                        q -= dm * ICE_LATENT_FUS_J_KG * a / h;
                    } else {
                        let m = (h_m * (rho_v_s - rho_v_amb)).min(c.frost_kg_m2 / h);
                        c.frost_kg_m2 -= m * h;
                        q -= m * ICE_LATENT_SUB_J_KG * a;
                    }
                }
            }
            c.enthalpy_j += q * h;
            let m_gold = c.struct_gold_kg + c.gold_kg;
            let c_other = c.c_other.max(1.0e-5);
            c.temperature_k = table.temperature_hunt(m_gold, c_other, c.enthalpy_j, &mut th.seg[i]);
            c.melt_fraction = if m_gold > 0.0 { gold.liquid_fraction_at(c.temperature_k) } else { 0.0 };
        }
        (torch_w, ln2_w, air_w)
    }

    /// Structural case gold whose column is half molten loses its shape and
    /// becomes free gold that can run.
    fn melt_case(&mut self) {
        for c in &mut self.grid {
            if c.struct_gold_kg > 0.0 && c.melt_fraction > 0.5 && c.kind == CellKind::Band {
                c.gold_kg += c.struct_gold_kg;
                c.struct_gold_kg = 0.0;
            }
        }
    }

    /// Room in a column for liquids above its bed [m].
    fn capacity_m(c: &SpatialCell) -> f64 {
        match c.kind {
            CellKind::Outside => f64::INFINITY,
            CellKind::Interior => OPEN_CEILING_M,
            CellKind::Band => {
                if c.struct_gold_kg <= 0.05 * c.struct_gold0_kg {
                    OPEN_CEILING_M
                } else {
                    (c.struct_gold0_kg - c.struct_gold_kg) / (library().gold.rho_solid_ref * CELL_AREA)
                }
            }
        }
    }

    /// Virtual-pipe shallow-layer flow of the gold and the oil.
    fn flow_liquids(&mut self, dt: f64) {
        let any = self.grid.iter().any(|c| (c.gold_kg > 0.0 && c.mobile_fraction() > 0.0) || c.oil_m3 > 0.0)
            || self.inputs.oil_pour_rate > 0.0;
        if !any {
            self.gold_qx.fill(0.0);
            self.gold_qy.fill(0.0);
            return;
        }
        // Oil pour source: a dropper stream of 0.25 mL/s at full rate.
        if self.inputs.oil_pour_rate > 0.0 {
            let p = self.inputs.oil_pos_mm;
            let rate = 0.25e-6 * self.inputs.oil_pour_rate;
            self.deposit(p, 1.2, |c, w| c.oil_m3 += rate * dt * w);
        }
        // Substeps: CFL on the gravity-wave and flow speeds.
        let g = self.gravity();
        let gmag = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt().max(1.0);
        let mut hmax: f64 = 1e-4;
        for c in &self.grid {
            hmax = hmax.max(c.gold_depth_m() * c.mobile_fraction()).max(c.oil_depth_m());
        }
        let c_wave = (gmag * hmax).sqrt() + 1.0;
        let n_sub = ((dt / (0.4 * DX / c_wave)).ceil() as usize).clamp(1, 400);
        let h = dt / n_sub as f64;
        for i in 0..GRID_CELLS {
            let c = &self.grid[i];
            self.gas_q[i] = self.torch_at(c).1 + self.jet_at(c).1;
        }
        // A layer with nothing that can move is skipped (its face fluxes are
        // zero); gold that freezes mid-call simply stops moving.
        let gold_moves = self.grid.iter().any(|c| c.gold_kg > 0.0 && c.mobile_fraction() > 0.0);
        let oil_moves = self.grid.iter().any(|c| c.oil_m3 > 0.0);
        if !gold_moves {
            self.gold_qx.fill(0.0);
            self.gold_qy.fill(0.0);
        }
        if !oil_moves {
            self.oil_qx.fill(0.0);
            self.oil_qy.fill(0.0);
        }
        let mut sc = std::mem::take(&mut self.pipe);
        for _ in 0..n_sub {
            if gold_moves {
                self.pipe_substep(h, true, &mut sc);
            }
            if oil_moves {
                self.pipe_substep(h, false, &mut sc);
            }
        }
        self.pipe = sc;
        // Liquids fall off a dial that faces down.
        if g[2] > 0.5 {
            let k = (g[2] / G0_M_S2) * 2.0;
            for c in &mut self.grid {
                if c.kind != CellKind::Outside {
                    let f = (-k * dt).exp();
                    let lost = c.gold_kg * c.mobile_fraction() * (1.0 - f);
                    c.gold_kg -= lost;
                    self.dripped_gold_g += lost * 1e3;
                    c.oil_m3 *= f;
                }
            }
        }
    }

    fn pipe_substep(&mut self, h: f64, is_gold: bool, sc: &mut PipeScratch) {
        let lib = library();
        let gold = &lib.gold;
        let table = column_table();
        let n = GRID_N;
        let g = self.gravity();
        let g_ax = g[2].abs();
        let (rho_g, rho_o) = (gold.rho_solid_ref * 0.95, lib.oil.density_ref);
        let rho = if is_gold { rho_g } else { rho_o };
        // Depths everywhere (cheap); the wet bounding box bounds the rest.
        let PipeScratch { depth, head, nu, out_v, in_v, s_out, s_in, moves } = sc;
        let mut lo = (n, n);
        let mut hi = (0usize, 0usize);
        for i in 0..GRID_CELLS {
            let c = &self.grid[i];
            let d = if is_gold { c.gold_depth_m() } else { c.oil_depth_m() };
            depth[i] = d;
            if d > 0.0 && c.kind != CellKind::Outside {
                let (gx, gy) = (i % n, i / n);
                lo = (lo.0.min(gx), lo.1.min(gy));
                hi = (hi.0.max(gx), hi.1.max(gy));
            }
        }
        if lo.0 > hi.0 {
            return;
        }
        let (x0, y0) = (lo.0.saturating_sub(2), lo.1.saturating_sub(2));
        let (x1, y1) = ((hi.0 + 2).min(n - 1), (hi.1 + 2).min(n - 1));
        // Head over the box and one cell beyond (the box's outer faces);
        // mobility only where there is liquid to move.
        for gy in y0..=(y1 + 1).min(n - 1) {
            for gx in x0..=(x1 + 1).min(n - 1) {
                let i = gy * n + gx;
                let c = &self.grid[i];
                let d = depth[i];
                let below = if is_gold { 0.0 } else { c.gold_depth_m() };
                nu[i] = if d >= 1.0e-7 {
                    let mu = if is_gold {
                        gold.effective_viscosity(c.temperature_k, c.melt_fraction)
                    } else {
                        lib.oil.viscosity(c.temperature_k)
                    };
                    mu / rho
                } else {
                    0.0
                };
                let bed = if c.kind == CellKind::Outside { -10.0e-3 } else { c.bed_mm * 1.0e-3 + below };
                head[i] = -(g[0] * c.x_mm + g[1] * c.y_mm) * 1.0e-3 + g_ax * (bed + d) + self.gas_q[i] / rho;
            }
        }
        // Capillary pressure of the free surface: −(σ/ρ)∇²h over wet cells.
        let sigma_over_rho = if is_gold { gold.surface_tension(1_200.0) / rho_g } else { 0.03 / rho_o };
        for gy in y0..=y1 {
            for gx in x0..=x1 {
                let i = gy * n + gx;
                if depth[i] <= 0.0 {
                    continue;
                }
                let mut lap = 0.0;
                for (j, ok) in [(i.wrapping_sub(1), gx > 0), (i + 1, gx + 1 < n), (i.wrapping_sub(n), gy > 0), (i + n, gy + 1 < n)] {
                    if ok && depth[j] > 0.0 {
                        lap += depth[j] - depth[i];
                    }
                }
                head[i] -= sigma_over_rho * lap / (DX * DX);
            }
        }
        let h_pin = if is_gold { 1.0e-3 } else { 0.0 };
        let (qx, qy) = if is_gold { (&mut self.gold_qx, &mut self.gold_qy) } else { (&mut self.oil_qx, &mut self.oil_qy) };
        let face = |q: &mut f64, a: usize, b: usize| {
            let dh = head[a] - head[b];
            let (src, dst) = if dh >= 0.0 { (a, b) } else { (b, a) };
            let hf = depth[src];
            let is_out = |k: usize| self.grid[k].kind == CellKind::Outside;
            if hf < 1.0e-7 || is_out(src) || (depth[dst] < 1.0e-6 && hf < h_pin && !is_out(dst)) {
                *q = 0.0;
                return;
            }
            let u = q.abs() / (DX * hf);
            let fric = 3.0 * nu[src] / (hf * hf) + 0.01 * u / hf;
            *q = (*q + h * hf * (head[a] - head[b])) / (1.0 + h * fric);
            // A layer cannot flow against its head.
            if (*q > 0.0) != (dh > 0.0) {
                *q = 0.0;
            }
        };
        for gy in y0..=y1 {
            for gx in x0..=x1 {
                let i = gy * n + gx;
                if gx + 1 < n {
                    face(&mut qx[i], i, i + 1);
                }
                if gy + 1 < n {
                    face(&mut qy[i], i, i + n);
                }
            }
        }
        // Outflow limited by what each cell holds; inflow by the room it has.
        for gy in y0..=y1 {
            let row = gy * n;
            out_v[row + x0..=row + x1].fill(0.0);
            in_v[row + x0..=row + x1].fill(0.0);
        }
        let mut tally = |q: f64, a: usize, b: usize| {
            let (s, d) = if q > 0.0 { (a, b) } else { (b, a) };
            out_v[s] += q.abs() * h;
            in_v[d] += q.abs() * h;
        };
        for gy in y0..=y1 {
            for gx in x0..=x1 {
                let i = gy * n + gx;
                if gx + 1 < n && qx[i] != 0.0 {
                    tally(qx[i], i, i + 1);
                }
                if gy + 1 < n && qy[i] != 0.0 {
                    tally(qy[i], i, i + n);
                }
            }
        }
        for gy in y0..=y1 {
            for gx in x0..=x1 {
                let i = gy * n + gx;
                s_out[i] = 1.0;
                s_in[i] = 1.0;
                let avail = depth[i] * CELL_AREA;
                if out_v[i] > avail {
                    s_out[i] = avail / out_v[i];
                }
                let c = &self.grid[i];
                let room = (Self::capacity_m(c) - c.gold_depth_m() - c.oil_depth_m()).max(0.0) * CELL_AREA;
                if in_v[i] > room {
                    s_in[i] = room / in_v[i];
                }
            }
        }
        // Move volume, and with gold its enthalpy.
        moves.clear();
        for gy in y0..=y1 {
            for gx in x0..=x1 {
                let i = gy * n + gx;
                for (q, j, ok) in [(&mut qx[i], i + 1, gx + 1 < n), (&mut qy[i], i + n, gy + 1 < n)] {
                    if !ok || *q == 0.0 {
                        continue;
                    }
                    let (s, d) = if *q > 0.0 { (i, j) } else { (j, i) };
                    let f = s_out[s].min(s_in[d]);
                    *q *= f;
                    moves.push((s, d, q.abs() * h));
                }
            }
        }
        for &(s, d, vol) in moves.iter() {
            if is_gold {
                let rho = rho_g;
                let m = (vol * rho).min(self.grid[s].gold_kg);
                let ts = self.grid[s].temperature_k;
                let hm = gold.enthalpy_at(ts) * m;
                self.grid[s].gold_kg -= m;
                self.grid[s].enthalpy_j -= hm;
                if self.grid[d].kind == CellKind::Outside {
                    self.dripped_gold_g += m * 1.0e3;
                } else {
                    self.grid[d].gold_kg += m;
                    self.grid[d].enthalpy_j += hm;
                }
            } else {
                let v = vol.min(self.grid[s].oil_m3);
                self.grid[s].oil_m3 -= v;
                if self.grid[d].kind != CellKind::Outside {
                    self.grid[d].oil_m3 += v;
                }
            }
        }
        // Columns that gained or lost gold re-equilibrate their temperature.
        if is_gold {
            for &(s, d, _) in moves.iter() {
                for k in [s, d] {
                    let c = &mut self.grid[k];
                    if c.kind == CellKind::Outside {
                        continue;
                    }
                    let m = c.struct_gold_kg + c.gold_kg;
                    c.temperature_k = table.temperature(m, c.c_other.max(1e-5), c.enthalpy_j);
                    c.melt_fraction = if m > 0.0 { gold.liquid_fraction_at(c.temperature_k) } else { 0.0 };
                }
            }
        }
    }

    /// Spread `f(cell, weight)` over a Gaussian footprint of width `sigma_mm`
    /// centred under `p` (weights sum to one over the case).
    fn deposit(&mut self, p: [f64; 3], sigma_mm: f64, mut f: impl FnMut(&mut SpatialCell, f64)) {
        let mut w = Vec::new();
        let mut sum = 0.0;
        for (i, c) in self.grid.iter().enumerate() {
            if c.kind == CellKind::Outside {
                continue;
            }
            let r2 = (c.x_mm - p[0]).powi(2) + (c.y_mm - p[1]).powi(2);
            if r2 < 16.0 * sigma_mm * sigma_mm {
                let wi = (-r2 / (2.0 * sigma_mm * sigma_mm)).exp();
                w.push((i, wi));
                sum += wi;
            }
        }
        if sum <= 0.0 {
            return;
        }
        for (i, wi) in w {
            f(&mut self.grid[i], wi / sum);
        }
        self.quiet = false;
    }

    /// Liquid nitrogen: the pour, and its run-off as Leidenfrost drops and
    /// rivulets along gravity, dripping off the edge.
    fn flow_ln2(&mut self, dt: f64) {
        if self.inputs.ln2_pour_rate > 0.0 {
            let p = self.inputs.ln2_pos_mm;
            let rate = 4.0e-3 * self.inputs.ln2_pour_rate;
            // A 5 mL/s stream landing at ~1.5 m/s splashes into a sheet
            // several jet diameters across.
            self.deposit(p, 4.0, |c, w| c.ln2_kg += rate * dt * w);
        }
        if !self.grid.iter().any(|c| c.ln2_kg > 0.0) {
            return;
        }
        let g = self.gravity();
        let g_ip = (g[0] * g[0] + g[1] * g[1]).sqrt();
        let g_mag = (g_ip * g_ip + g[2] * g[2]).sqrt().max(1e-9);
        // Leidenfrost rivulets run down a tilted face at ~0.15 m/s per g.
        let u = (0.15 * g_ip / G0_M_S2).min(0.5);
        let n = GRID_N;
        let frac = (u * dt / DX).min(0.9);
        let (fx, fy) = if g_ip > 0.0 { (g[0] / g_ip, g[1] / g_ip) } else { (0.0, 0.0) };
        // The face-normal share of gravity levels the free surface: liquid
        // runs down head differences (bed + depth), fills recesses such as
        // the tourbillon well, and spills over the bezel lip only when deeper
        // than it.
        let normal = (g[2].abs() / g_mag).clamp(0.0, 1.0);
        const LIP_M: f64 = 0.6e-3;
        let rho_a = ln2::RHO_LIQUID * CELL_AREA;
        let head: Vec<f64> = self.grid.iter().map(|c| c.bed_mm * 1e-3 + c.ln2_kg / rho_a).collect();
        let mut next: Vec<f64> = self.grid.iter().map(|c| c.ln2_kg).collect();
        for gy in 0..n {
            for gx in 0..n {
                let i = gy * n + gx;
                let m = self.grid[i].ln2_kg;
                if m <= 0.0 {
                    continue;
                }
                let send = |j: Option<usize>, amt: f64, next: &mut Vec<f64>| {
                    next[i] -= amt;
                    if let Some(j) = j {
                        if self.grid[j].kind != CellKind::Outside {
                            next[j] += amt;
                        }
                    }
                };
                let xn = if fx >= 0.0 { (gx + 1 < n).then(|| i + 1) } else { (gx > 0).then(|| i - 1) };
                let yn = if fy >= 0.0 { (gy + 1 < n).then(|| i + n) } else { (gy > 0).then(|| i - n) };
                send(xn, m * frac * fx.abs(), &mut next);
                send(yn, m * frac * fy.abs(), &mut next);
                if normal > 0.05 {
                    let depth = m / rho_a;
                    for j in [(gx > 0).then(|| i - 1), (gx + 1 < n).then(|| i + 1), (gy > 0).then(|| i - n), (gy + 1 < n).then(|| i + n)]
                        .into_iter()
                        .flatten()
                    {
                        let dh = if self.grid[j].kind == CellKind::Outside {
                            depth - LIP_M
                        } else {
                            head[i] - head[j]
                        };
                        if dh > 0.0 {
                            // Relax a fifth of the step in head per step (stable for 4 neighbours).
                            send(Some(j), (0.2 * normal * dh * rho_a).min(0.24 * m), &mut next);
                        }
                    }
                }
            }
        }
        for (c, m) in self.grid.iter_mut().zip(next) {
            c.ln2_kg = if c.kind == CellKind::Outside { 0.0 } else { m.max(0.0) };
        }
    }

    /// Dipole field of the N42 magnet [T] at `p` (watch frame, mm).
    fn magnet_field(&self, p: [f64; 3]) -> [f64; 3] {
        let i = &self.inputs;
        if !i.magnet_enabled || i.magnet_strength <= 0.0 {
            return [0.0; 3];
        }
        let lib = library();
        // Ø12 × 6 mm cylinder at room temperature; strength scales its volume.
        let vol = PI * 36.0e-6 * 6.0e-3 * i.magnet_strength;
        let moment = lib.n42.remanence_t(T_REF_K) * vol / MU0;
        let m_pos = i.magnet_pos_mm;
        // Axis points from the magnet at the tourbillon.
        let mut axis = [TOURBILLON_CENTER_MM[0] - m_pos[0], TOURBILLON_CENTER_MM[1] - m_pos[1], -m_pos[2]];
        let al = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt().max(1e-9);
        axis = axis.map(|a| a / al);
        let r = [(p[0] - m_pos[0]) * 1e-3, (p[1] - m_pos[1]) * 1e-3, (p[2] - m_pos[2]) * 1e-3];
        let rl = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]).sqrt().max(4.0e-3);
        let rh = r.map(|x| x / rl);
        let mdotr = axis[0] * rh[0] + axis[1] * rh[1] + axis[2] * rh[2];
        let k = MU0 * moment / (4.0 * PI * rl.powi(3));
        [0, 1, 2].map(|j| k * (3.0 * mdotr * rh[j] - axis[j]))
    }

    fn b_mag(&self, p: [f64; 3]) -> f64 {
        let b = self.magnet_field(p);
        (b[0] * b[0] + b[1] * b[1] + b[2] * b[2]).sqrt()
    }

    fn update_magnetics(&mut self, dt: f64, cage_angle: f64, lever_angle: f64) {
        let tc = TOURBILLON_CENTER_MM;
        let p_bal = [tc[0], tc[1], 0.6];
        let b = self.magnet_field(p_bal);
        let b_signed = b[0];
        // AC degauss: a decaying 50 Hz field swept through the play operator.
        let hs = if self.inputs.compensated_hairspring {
            &library().nivarox_hairspring
        } else {
            &library().carbon_steel_hairspring
        };
        let h_sat = 5.0 * hs.coercive_a_m;
        let play = |h_field: f64, m: &mut f64| {
            let lo = (h_field - hs.coercive_a_m) / h_sat;
            let hi = (h_field + hs.coercive_a_m) / h_sat;
            *m = m.clamp(lo, hi).clamp(-1.0, 1.0);
        };
        let mut m = self.hairspring_magnetisation;
        if self.degauss_t >= 0.0 {
            let n = (dt / 5.0e-4).ceil().max(1.0) as usize;
            let hstep = dt / n as f64;
            for _ in 0..n {
                self.degauss_t += hstep;
                let amp = DEGAUSS_PEAK_T * (-self.degauss_t / DEGAUSS_TAU_S).exp();
                let bac = amp * (TAU * 50.0 * self.degauss_t).sin();
                play((b_signed + bac) / MU0, &mut m);
            }
            if self.degauss_t >= DEGAUSS_DURATION_S {
                self.degauss_t = -1.0;
            }
        } else {
            play(b_signed / MU0, &mut m);
        }
        self.hairspring_magnetisation = m;
        let b_abs = (b[0] * b[0] + b[1] * b[1] + b[2] * b[2]).sqrt();
        self.b_field_escapement_mt = b_abs * 1.0e3;
        self.remanent_magnetization_mt = m.abs() * hs.saturation_t * 1.0e3;
        // Gradient magnitude at the balance (central differences, 0.2 mm).
        let d = 0.2;
        let grad = [
            (self.b_mag([p_bal[0] + d, p_bal[1], p_bal[2]]) - self.b_mag([p_bal[0] - d, p_bal[1], p_bal[2]])) / (2.0 * d * 1e-3),
            (self.b_mag([p_bal[0], p_bal[1] + d, p_bal[2]]) - self.b_mag([p_bal[0], p_bal[1] - d, p_bal[2]])) / (2.0 * d * 1e-3),
            (self.b_mag([p_bal[0], p_bal[1], p_bal[2] + d]) - self.b_mag([p_bal[0], p_bal[1], p_bal[2] - d])) / (2.0 * d * 1e-3),
        ];
        let grad_mag = (grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2]).sqrt();
        self.environment.b_balance_t = b_abs;
        self.environment.grad_b_balance_t_m = grad_mag;
        self.environment.hairspring_magnetisation = m;
        self.environment.lever_magnetic_torque_nm = self.lever_torque(cage_angle, lever_angle);
    }

    /// Magnetic torque on the soft-steel pallet lever: each lump of steel is
    /// magnetised by the local field (demagnetising factor 0.3, saturating at
    /// 1.6 T) and pulled up the gradient of |B|; the torque is taken about the
    /// lever arbor, in the cage that carries it.
    fn lever_torque(&self, cage_angle: f64, lever_angle: f64) -> f64 {
        if !self.inputs.magnet_enabled {
            return 0.0;
        }
        let tc = TOURBILLON_CENTER_MM;
        let steel_rho = library().steel.rho;
        // Lever lumps in the cage frame relative to its arbor [mm] and mass [kg].
        let lumps: [([f64; 2], f64); 4] = [
            ([0.0, 3.0], 2.0e-6),
            ([0.0, 1.2], 2.0e-6),
            ([1.4, -0.65], 2.0e-6),
            ([-1.4, -0.65], 2.0e-6),
        ];
        let arbor_cage = [0.0, -4.2522];
        let rot = |v: [f64; 2], a: f64| [v[0] * a.cos() - v[1] * a.sin(), v[0] * a.sin() + v[1] * a.cos()];
        let to_watch = |v: [f64; 2]| {
            let w = rot(v, -cage_angle);
            [tc[0] + w[0], tc[1] + w[1], 0.2]
        };
        let pivot = to_watch(arbor_cage);
        let mut tau = 0.0;
        for (off, mass) in lumps {
            let local = rot(off, lever_angle);
            let p = to_watch([arbor_cage[0] + local[0], arbor_cage[1] + local[1]]);
            let bm = self.b_mag(p);
            let m_mag = (bm / (MU0 * 0.3)).min(1.6 / MU0);
            let d = 0.1;
            let gx = (self.b_mag([p[0] + d, p[1], p[2]]) - self.b_mag([p[0] - d, p[1], p[2]])) / (2.0 * d * 1e-3);
            let gy = (self.b_mag([p[0], p[1] + d, p[2]]) - self.b_mag([p[0], p[1] - d, p[2]])) / (2.0 * d * 1e-3);
            let v = mass / steel_rho;
            let f = [v * m_mag * gx, v * m_mag * gy];
            let r = [(p[0] - pivot[0]) * 1e-3, (p[1] - pivot[1]) * 1e-3];
            // Watch-frame torque; the cage frame is turned by −cage_angle,
            // which leaves a z torque unchanged.
            tau += r[0] * f[1] - r[1] * f[0];
        }
        tau
    }

    fn update_case_shell(&mut self, dt: f64) {
        let gold = &library().gold;
        let g = self.gravity();
        let gm = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt().max(1e-9);
        let gh = g.map(|x| x / gm);
        let l = 8.0e-3;
        for node in &mut self.case_shell {
            let c = &self.grid[node.cell];
            node.temperature_k = c.temperature_k;
            node.melt_fraction = c.melt_fraction;
            let frac = if c.struct_gold0_kg > 0.0 { c.struct_gold_kg / c.struct_gold0_kg } else { 1.0 };
            node.thickness_mm = node.ref_thickness_mm * frac.max(0.0);
            // Viscous slab sag under self-weight: v ≈ ρ g L² / (8 μ).
            let mu = gold.effective_viscosity(c.temperature_k, c.melt_fraction);
            let v = (gold.rho_solid_ref * gm * l * l / (8.0 * mu)).min(0.05);
            node.sag_mm = (node.sag_mm + v * dt * 1.0e3).min(if frac < 0.05 { 9.0 } else { 6.0 });
            for k in 0..3 {
                node.pos_mm[k] = node.ref_pos_mm[k] + gh[k] * node.sag_mm;
            }
        }
    }

    /// Sample a part's footprint into the movement's environment.
    fn sample_part(&self, cells: &[usize], floor_mm: f64) -> PartField {
        if cells.is_empty() {
            return PartField::default();
        }
        let n = cells.len() as f64;
        let mut f = PartField { temperature_k: 0.0, gold_temperature_k: 0.0, ..PartField::default() };
        let mut fl = 0.0;
        for &i in cells {
            let c = &self.grid[i];
            f.temperature_k += c.temperature_k / n;
            let lift = (c.bed_mm - floor_mm).max(0.0) * 1.0e-3;
            let gd = c.gold_depth_m();
            let mob = c.mobile_fraction();
            if gd > 0.0 {
                f.liquid_gold_m += (lift + gd) * mob / n;
                f.solid_gold_m += (lift + gd) * (1.0 - mob) / n;
            }
            let od = c.oil_depth_m();
            if od > 0.0 {
                f.oil_m += (lift + gd + od) / n;
            }
            f.frost_m += c.frost_kg_m2 / FROST_DENSITY / n;
            fl += c.melt_fraction / n;
            f.gold_temperature_k += c.temperature_k / n;
        }
        f.gold_liquid_fraction = fl;
        f
    }

    fn refresh_outputs(&mut self, _cage_angle: f64, _lever_angle: f64) {
        let lib = library();
        let gold = &lib.gold;
        let floor = WELL_FLOOR_MM;
        let mut env = self.environment;
        env.balance = self.sample_part(&self.part_cells[PART_BALANCE], floor);
        env.escapement = self.sample_part(&self.part_cells[PART_ESCAPEMENT], floor);
        env.train = self.sample_part(&self.part_cells[PART_TRAIN], floor);
        let mut barrel = self.sample_part(&self.part_cells[PART_BARREL], 0.0);
        // Liquids on the dial do not reach the barrel; frozen gold round the
        // pendant does block the keyless works.
        let crown = self.sample_part(&self.crown_cells, 0.0);
        barrel.liquid_gold_m = 0.0;
        barrel.oil_m = 0.0;
        barrel.solid_gold_m = crown.solid_gold_m;
        env.barrel = barrel;
        env.air_pressure_pa = self.pressure_pa();
        env.air_temperature_k = env.balance.temperature_k;
        env.gravity_m_s2 = self.gravity();
        // Air-jet speed where it crosses the balance (free-jet decay).
        let i = &self.inputs;
        env.air_jet_speed_m_s = if i.air_jet_power > 0.0 {
            let p = i.air_jet_pos_mm;
            let d = ((p[0] - TOURBILLON_CENTER_MM[0]).powi(2) + (p[1] - TOURBILLON_CENTER_MM[1]).powi(2) + (p[2] - 0.6).powi(2))
                .sqrt()
                .max(1.0);
            150.0 * i.air_jet_power * (6.2 * 1.5 / d).min(1.0)
        } else {
            0.0
        };
        self.environment = env;

        // Twelve region nodes.
        let mut acc = [ThermalPhaseNode::default(); NODE_COUNT];
        let mut cnt = [0.0_f64; NODE_COUNT];
        for c in &self.grid {
            let mut slots: [Option<usize>; 2] = [None, None];
            slots[0] = match c.region {
                Region::Outside => None,
                Region::CaseTop => Some(0),
                Region::CaseLeft => Some(1),
                Region::CaseRight => Some(2),
                Region::CaseBottom => Some(3),
                Region::Plate => Some(5),
                Region::Barrel => Some(6),
                Region::Train => Some(7),
                Region::Cage => Some(8),
                Region::Escapement => Some(9),
                Region::Balance => Some(10),
            };
            if c.in_hairspring {
                slots[1] = Some(11);
            } else if c.kind == CellKind::Interior {
                slots[1] = Some(4);
            }
            for k in slots.into_iter().flatten() {
                let a = &mut acc[k];
                a.temperature_k += c.temperature_k;
                a.melt_fraction += c.melt_fraction;
                a.plastic_slump += c.case_slump_frac();
                a.molten_gold_g += c.gold_kg * c.mobile_fraction() * 1.0e3;
                a.frozen_gold_g += c.gold_kg * (1.0 - c.mobile_fraction()) * 1.0e3;
                a.ln2_film_g += c.ln2_kg * 1.0e3;
                a.frost_thickness_um += c.frost_um();
                a.oil_film_um += c.oil_depth_m() * 1.0e6;
                cnt[k] += 1.0;
            }
        }
        for k in 0..NODE_COUNT {
            let n = cnt[k].max(1.0);
            let a = &mut acc[k];
            a.temperature_k /= n;
            a.melt_fraction /= n;
            a.plastic_slump /= n;
            a.frost_thickness_um /= n;
            a.oil_film_um /= n;
            a.enthalpy_j_kg = gold.enthalpy_at(a.temperature_k);
            let eps = if k <= 3 { gold.emissivity(a.melt_fraction, 50.0) } else { 0.4 };
            a.incandescence_rgb = lib.blackbody.display_rgb(a.temperature_k, eps);
        }
        self.nodes = acc;
        self.blued_steel_rgb = blued_steel_interference_rgb(self.hairspring_oxide_nm, 0.85);
    }

    /// Loudness of the degauss coil, 0–1, following its field amplitude.
    #[must_use]
    pub fn degauss_level(&self) -> f64 {
        if self.degauss_t >= 0.0 { (-self.degauss_t / DEGAUSS_TAU_S).exp() } else { 0.0 }
    }

    /// Average temperature of the balance region [K].
    #[must_use]
    pub fn balance_temperature_k(&self) -> f64 {
        self.environment.balance.temperature_k
    }
}

impl Default for MultiphysicsState {
    fn default() -> Self {
        Self::new()
    }
}

/// Grid cell containing watch-frame point (x, y) [mm].
#[must_use]
pub fn cell_at(x_mm: f64, y_mm: f64) -> Option<usize> {
    let gx = ((x_mm + GRID_HALF_SPAN_MM) / GRID_DX_MM).floor();
    let gy = ((y_mm + GRID_HALF_SPAN_MM) / GRID_DX_MM).floor();
    if gx < 0.0 || gy < 0.0 || gx >= GRID_N as f64 || gy >= GRID_N as f64 {
        return None;
    }
    Some(gy as usize * GRID_N + gx as usize)
}

/// Column enthalpy at `t_k` for a cell's current contents (used by tests and
/// for setting initial states).
#[must_use]
pub fn column_enthalpy(cell: &SpatialCell, t_k: f64) -> f64 {
    column_table().enthalpy(cell.struct_gold_kg + cell.gold_kg, cell.c_other.max(1e-5), t_k)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn column_enthalpy_inverts_through_the_mushy_zone() {
        let t = column_table();
        for &tk in &[77.0, 200.0, 293.15, 900.0, 1158.15, 1162.0, 1168.15, 1400.0] {
            let h = t.enthalpy(1.0e-3, 2.0e-3, tk);
            let back = t.temperature(1.0e-3, 2.0e-3, h);
            assert!((back - tk).abs() < 1e-6, "{tk} → {back}");
        }
        let mushy = t.capacity(1e-3, 0.0, 1162.0);
        let solid = t.capacity(1e-3, 0.0, 900.0);
        assert!(mushy > 15.0 * solid, "latent heat shows as capacity: {mushy} vs {solid} J/K");
    }

    #[test]
    fn ln2_boiling_curve_has_its_three_regimes() {
        let p = P_ATM_PA;
        let q5 = ln2_boiling_flux(5.0, p, G0_M_S2);
        let q_chf = (1..60).map(|d| ln2_boiling_flux(d as f64, p, G0_M_S2)).fold(0.0_f64, f64::max);
        let q200 = ln2_boiling_flux(200.0, p, G0_M_S2);
        assert!(q_chf > 1.2e5 && q_chf < 2.6e5, "CHF {q_chf}");
        assert!(q5 > 1.0e4 && q5 < q_chf, "nucleate {q5}");
        assert!(q200 > 5.0e3 && q200 < 6.0e4, "film {q200}");
        // The Leidenfrost minimum lies between them.
        let q_min = (20..120).map(|d| ln2_boiling_flux(d as f64, p, G0_M_S2)).fold(f64::INFINITY, f64::min);
        assert!(q_min < q200 && q_min < q_chf / 5.0, "minimum {q_min}");
    }

    #[test]
    fn idle_watch_stays_at_room_temperature() {
        let mut m = MultiphysicsState::new();
        for _ in 0..20 {
            m.step(0.05, 0.0, 0.0);
        }
        let t = m.environment.balance.temperature_k;
        assert!((t - T_REF_K).abs() < 1e-6, "balance {t}");
        assert!(m.quiet, "idle grid goes quiet");
    }

    #[test]
    fn torch_heats_and_ln2_cools_the_movement() {
        let mut m = MultiphysicsState::new();
        m.inputs.blowtorch_power = 1.0;
        m.inputs.blowtorch_pos_mm = [0.0, -11.0, 3.0];
        for _ in 0..100 {
            m.step(0.05, 0.0, 0.0);
        }
        let hot = m.environment.balance.temperature_k;
        assert!(hot > 340.0, "torch on the tourbillon for 5 s: {hot} K");
        assert!(m.blowtorch_heating_power_w > 20.0, "torch power {}", m.blowtorch_heating_power_w);
        m.inputs.blowtorch_power = 0.0;
        m.inputs.ln2_pour_rate = 1.0;
        m.inputs.ln2_pos_mm = [0.0, -11.0, 3.0];
        for _ in 0..200 {
            m.step(0.05, 0.0, 0.0);
        }
        // A hot part under poured LN₂ sits in film boiling (Leidenfrost): a few
        // 10⁴ W/m² over the wetted area, so it cools by hundreds of kelvin in
        // ten seconds but not to cryogenic temperatures.
        let cold = m.environment.balance.temperature_k;
        assert!(cold < hot - 150.0, "LN₂ quench: {hot} K → {cold} K");
        let wetted: f64 = m.grid.iter().map(|c| (c.ln2_kg / (ln2::RHO_LIQUID * CELL_AREA * 0.2e-3)).clamp(0.0, 1.0)).sum::<f64>()
            * CELL_AREA;
        let flux = m.ln2_boiling_power_w / wetted.max(1e-12);
        assert!((5.0e3..2.0e5).contains(&flux), "film-boiling flux {flux} W/m²");
    }

    #[test]
    fn ln2_pooled_on_a_flat_watch_chills_the_tourbillon() {
        // On the bench, dial up: the nitrogen pools in the tourbillon recess.
        let mut m = MultiphysicsState::new();
        m.inputs.gravity_g_vec = [0.0, 0.0, -1.0];
        m.inputs.ln2_pour_rate = 1.0;
        m.inputs.ln2_pos_mm = [0.0, -11.0, 3.0];
        let mut trace = Vec::new();
        for i in 0..600 {
            m.step(0.05, 0.0, 0.0);
            if i % 100 == 99 {
                trace.push((m.environment.balance.temperature_k, m.ln2_boiling_power_w, m.environment.balance.frost_m * 1e6));
            }
        }
        println!("LN₂ on a cold watch, every 5 s (T K, P W, frost µm): {trace:.1?}");
        // Film boiling over one face of a ~16 J/K watch: tens of watts, a
        // couple of kelvin per second, steadily.
        let t = m.environment.balance.temperature_k;
        assert!(t < 250.0, "30 s of LN₂ on the tourbillon: {t} K");
        assert!((20.0..200.0).contains(&m.ln2_boiling_power_w), "boiling power {} W", m.ln2_boiling_power_w);
    }

    #[test]
    fn magnet_pulls_on_the_lever_and_magnetises_a_steel_spring() {
        let mut m = MultiphysicsState::new();
        m.inputs.compensated_hairspring = false;
        m.inputs.magnet_enabled = true;
        m.inputs.magnet_pos_mm = [13.0, -11.0, 2.5];
        m.step(0.05, 0.3, -0.09);
        assert!(m.b_field_escapement_mt > 20.0, "B {} mT", m.b_field_escapement_mt);
        assert!(m.environment.lever_magnetic_torque_nm.abs() > 1.0e-8);
        let mag = m.hairspring_magnetisation;
        m.inputs.magnet_enabled = false;
        m.step(0.05, 0.3, -0.09);
        assert!(m.hairspring_magnetisation.abs() > 0.0 && mag != 0.0, "remanence after removal");
        m.trigger_degauss();
        for _ in 0..60 {
            m.step(0.05, 0.3, -0.09);
        }
        assert!(m.hairspring_magnetisation.abs() < 0.05, "degaussed: {}", m.hairspring_magnetisation);
    }
}
