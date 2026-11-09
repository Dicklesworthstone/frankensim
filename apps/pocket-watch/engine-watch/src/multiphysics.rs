//! Coupled 3D Thermo-Elasto-Viscoplastic Case Shell (`64 x 8 = 512` 3D Lagrangian Shell Nodes)
//! + 2D Finite-Volume Free-Surface Capillary / Lubrication Multiphysics PDE Engine (`24 x 24 = 576` cells)
//! for the FrankenSim Pocket Watch (`Calibre FS-1903`).
//!
//! Zero canned shapes or scripted slumps: every geometric deformation of the 18k gold case
//! and every droplet/stream of molten and solidified gold emerges from integrating:
//! 1. **3D Thermo-Elasto-Viscoplasticity (Johnson-Cook + Norton-Bailey Creep + Prandtl-Reuss Flow)**
//!    on the $64 \times 8$ 3D Bassine case shell mesh:
//!    - Thermal expansion strain $\boldsymbol{\varepsilon}^{\text{th}}(T) = \alpha_{\text{Au}}(T)(T - T_0)\mathbf{I}$
//!      causing thermal bowing & compressive hoop buckling before solidus.
//!    - Homologous thermal softening of yield stress $\sigma_y(T, f_L, \bar{\varepsilon}^p)$ and
//!      Elastic Modulus $E_{\text{Au}}(T)(1 - f_L)^3$.
//!    - Associated Prandtl-Reuss viscoplastic strain rate $\dot{\boldsymbol{\varepsilon}}^{vp} = \dot{\bar{\varepsilon}}^{vp}(\sigma_{\text{eq}}, T, f_L) \frac{3\mathbf{s}}{2\sigma_{\text{eq}}}$.
//! 2. **Young-Laplace Surface Tension ($\sigma_{\text{Au}}(T)\kappa\hat{\mathbf{n}}$), Marangoni Shear ($\frac{d\sigma_{\text{Au}}}{dT}\nabla_s T$),
//!    and Wetting Contact-Angle Cohesion ($\theta_c(T_{\text{sub}})$)**:
//!    - Mean-curvature surface tension pulls melting gold rims into bulbous capillary beads
//!      (Rayleigh-Plateau rim retraction & beading) while gravity and blowtorch gas-jet stagnation
//!      pressure push the low-viscosity melt pool inward over the movement substrate.
//! 3. **2.5D Lubrication + Fourth-Order Young-Laplace Capillary + Disjoining Cohesion Free-Surface PDE**:
//!    $$\frac{\partial H_{\text{Au}}}{\partial t} + \nabla \cdot \left( -\frac{H_{\text{Au}}^3}{3\mu_{\text{eff}}(T, f_L)} \nabla P_{\text{tot}} + \frac{H_{\text{Au}}^2}{2\mu_{\text{eff}}(T, f_L)} \frac{d\sigma_{\text{Au}}}{dT}\nabla T \right) = \dot{S}_{\text{drip}} - \dot{S}_{\text{freeze}} + \dot{S}_{\text{remelt}}$$
//!    where $P_{\text{tot}} = \rho_{\text{Au}} g_z (z_{\text{bed}} + S_{\text{Au}} + H_{\text{Au}}) - \rho_{\text{Au}}(g_x x + g_y y) - \sigma_{\text{Au}}(T)\nabla^2 \eta + \Pi_{\text{coh}}(H_{\text{Au}}, \theta_c)$,
//!    coupled with enthalpy advection $\nabla \cdot (\mathbf{q}_{\text{Au}} \rho_{\text{Au}} h_{\text{Au}})$.
//! 4. **Spatially-Integrated Hydrodynamic Viscous Drag & Solidified Yield-Stress Interference Torques**
//!    over the exact 2D geometric footprints $\Omega_w$ of the Barrel, Going-Train Wheels,
//!    Tourbillon Carriage, Escape Wheel, Pallet Fork, Balance Wheel, and Breguet Hairspring.

use core::f64::consts::{PI, TAU};

use fs_ad::{Dual64, Real};
use fs_convection::{
    CorrelationId, CorrelationInputs, NaturalConvectionGeometry, ThermalConductivity,
    evaluate as eval_convection,
};
use fs_ivl::Interval;
use fs_lbm::{
    Cell as LbmCell, ContactModel as LbmContactModel, FreeSurface as LbmFreeSurface,
    Grid as LbmGrid,
};
use fs_material::gas::GasState;
use fs_mbd_material::VolumeMoments;
use fs_qty::Length as ConvLength;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;
use fs_tribo::journal::{
    GeneralizedNewtonianJournal, JournalGeometry, JournalShearSolve, JournalViscosityKnot,
};
use fs_vkplate::air_loss::oscillating_cylinder;

use crate::materials::{
    CreepState, GOLD_18K_CARD, GOLD_MARANGONI_DSIGMA_DT_N_M_K, HIGH_CARBON_STEEL_CARD,
    INVAR_36NI_CARD, J2State, LN2_LATENT_VAP_J_PER_KG, NortonJ2State, STEFAN_BOLTZMANN,
    T_LN2_LEIDENFROST_K, T_LN2_SAT_K, T_REF_K, WatchMaterialId, blued_steel_interference_rgb,
    gold_coupled_norton_j2, gold_ductility_elongation, gold_effective_viscosity_pa_s,
    gold_liquid_surface_tension_n_m, gold_substrate_contact_angle_cos,
    gold_viscoplastic_strain_rate_s_inv, gold_viscoplastic_yield_stress_pa,
    planck_incandescence_linear_rgb, watch_oil_surface_tension_n_m, watch_oil_viscosity_pa_s,
};

/// Number of component aggregate nodes reported to the telemetry ledger.
pub const NODE_COUNT: usize = 12;

/// Spatial PDE grid resolution along X and Y across $[-28.8\text{ mm}, +28.8\text{ mm}]$.
pub const GRID_N: usize = 24;
pub const GRID_CELLS: usize = GRID_N * GRID_N;
/// Grid physical half-span ($\text{mm}$). Cell pitch $\Delta x = 57.6 / 24 = 2.4\text{ mm}$.
pub const GRID_HALF_SPAN_MM: f64 = 28.8;
pub const GRID_DX_MM: f64 = (2.0 * GRID_HALF_SPAN_MM) / (GRID_N as f64);
pub const GRID_DX_M: f64 = GRID_DX_MM * 1.0e-3;
pub const GRID_CELL_AREA_M2: f64 = GRID_DX_M * GRID_DX_M;

/// 3D Deformable 18k Gold Case Shell resolution:
/// 64 circumferential stations around $[0, 2\pi)$ $\times$ 8 meridional profile stations from back bezel to front bezel.
pub const CASE_SHELL_N_THETA: usize = 64;
pub const CASE_SHELL_N_Z: usize = 8;
pub const CASE_SHELL_NODES: usize = CASE_SHELL_N_THETA * CASE_SHELL_N_Z;

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

impl WatchNodeId {
    pub const ALL: [Self; NODE_COUNT] = [
        Self::CaseTop18kGold,
        Self::CaseLeft18kGold,
        Self::CaseRight18kGold,
        Self::CaseBottom18kGold,
        Self::CrystalSapphire,
        Self::MainplateBrass,
        Self::BarrelMainspring,
        Self::GoingTrainWheels,
        Self::TourbillonCage,
        Self::EscapementForkWheel,
        Self::BalanceWheelRim,
        Self::HairspringBreguet,
    ];
}

/// Aggregated state of one of the 12 horological functional subsystems
/// (computed directly by integrating over the 3D Case Shell nodes and 2D Spatial PDE cells).
#[derive(Debug, Clone, Copy)]
pub struct ThermalPhaseNode {
    pub id: WatchNodeId,
    pub material: WatchMaterialId,
    pub pos_mm: [f64; 3],
    pub base_mass_kg: f64,
    pub surface_area_m2: f64,
    pub enthalpy_j_kg: f64,
    pub temperature_k: f64,
    pub melt_fraction: f64,
    pub plastic_slump: f64,
    pub molten_gold_g: f64,
    pub frozen_gold_g: f64,
    pub ln2_film_g: f64,
    pub frost_thickness_um: f64,
    pub oil_film_um: f64,
    pub oxide_thickness_nm: f64,
    pub incandescence_rgb: [f64; 3],
}

/// Interactive external perturbation tools applied by the user in the 3D scene.
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
            compensated_hairspring: false,
        }
    }
}

/// Coupled multiphysics modifiers transmitted to the multi-body horological solver & acoustics.
#[derive(Debug, Clone, Copy)]
pub struct CoupledWatchModifiers {
    pub balance_inertia_ratio: f64,
    pub hairspring_stiffness_ratio: f64,
    pub barrel_youngs_modulus_ratio: f64,
    pub natural_freq_ratio: f64,
    pub theoretical_rate_s_per_day: f64,
    pub pivot_friction_multiplier: f64,
    pub extra_balance_damping_s_inv: f64,
    pub hairspring_oil_bridge_fraction: f64,
    pub min_pivot_clearance_um: f64,
    pub b_field_escapement_mt: f64,
    pub remanent_magnetization_mt: f64,
    pub pallet_magnetic_bias_torque_nm: f64,
    pub hairspring_eq_bias_rad: f64,
    /// Continuous emergent fouling resistance factor $\in [0, 1]$ derived from the actual
    /// viscous & solid shear resistive torque $\tau_{\text{Au}}$ vs mainspring drive torque $\tau_{\text{drive}}$.
    pub gold_fouling_jam_factor: f64,
    /// Spatially-integrated viscous drag torque ($\text{N}\cdot\text{m}$) exerted by mobile liquid gold
    /// on the going train and escapement wheels.
    pub liquid_gold_viscous_torque_nm: f64,
    /// Spatially-integrated solid shear-yield interference torque ($\text{N}\cdot\text{m}$) exerted by
    /// solidified gold frozen inside the wheel teeth or Tourbillon cage.
    pub solid_gold_yield_torque_nm: f64,
    pub frozen_gold_locked: bool,
    pub thermal_or_cryo_seized: bool,
    pub magnetically_pinned: bool,
    pub ln2_boiling_power_w: f64,
    pub blowtorch_heating_power_w: f64,
    pub blued_steel_rgb: [f64; 3],
    /// Exact forward-mode automatic differentiation sensitivity $\partial r / \partial T$ ($\text{s/day/K}$) from `fs-ad::Dual64`.
    pub rate_sensitivity_s_day_per_k: f64,
    /// Certified outward-rounded interval lower bound on chronometric rate ($\text{s/day}$) from `fs-ivl::Interval`.
    pub rate_interval_lo_s_day: f64,
    /// Certified outward-rounded interval upper bound on chronometric rate ($\text{s/day}$) from `fs-ivl::Interval`.
    pub rate_interval_hi_s_day: f64,
    /// Exact Stokes-Stuart oscillating-cylinder aerodynamic drag ($\text{s}^{-1}$) on the balance wheel from `fs-vkplate::air_loss`.
    pub balance_aero_drag_s_inv: f64,
}

/// Single 3D Lagrangian material node on the $64 \times 8$ deformable 18k Gold Bassine Case Shell.
#[derive(Debug, Clone, Copy)]
pub struct CaseShellNode {
    pub theta_rad: f64,
    /// Undeformed reference 3D coordinates $[X_0, Y_0, Z_0]$ ($\text{mm}$).
    pub ref_pos_mm: [f64; 3],
    /// Current deformed 3D coordinates $[x, y, z]$ ($\text{mm}$).
    pub pos_mm: [f64; 3],
    /// Current 3D material velocity $[v_x, v_y, v_z]$ ($\text{mm/s}$).
    pub vel_mm_s: [f64; 3],
    /// Undeformed reference wall thickness $h_0$ ($\text{mm}$).
    pub ref_thickness_mm: f64,
    /// Current local wall thickness $h$ ($\text{mm}$).
    pub thickness_mm: f64,
    /// Reference circumferential segment length $L_{\theta,0}$ ($\text{mm}$).
    pub ref_l_theta_mm: f64,
    /// Reference meridional segment length $L_{z,0}$ ($\text{mm}$).
    pub ref_l_z_mm: f64,
    /// Current nodal gold mass $m$ ($\text{kg}$).
    pub mass_kg: f64,
    /// Surface area exposed to external heat flux $A$ ($\text{m}^2$).
    pub area_m2: f64,
    /// Accumulated viscoplastic hoop strain $\varepsilon_\theta^p$.
    pub eps_p_theta: f64,
    /// Accumulated viscoplastic meridional strain $\varepsilon_z^p$.
    pub eps_p_z: f64,
    /// Accumulated viscoplastic 3D bending/shear offset $[u_x^p, u_y^p, u_z^p]$ ($\text{mm}$).
    pub plastic_disp_mm: [f64; 3],
    /// Equivalent von Mises plastic strain $\bar{\varepsilon}^p \ge 0$.
    pub eq_plastic_strain: f64,
    /// Current von Mises equivalent Cauchy stress $\sigma_{\text{eq}}$ ($\text{Pa}$).
    pub von_mises_stress_pa: f64,
    /// Specific enthalpy $H$ ($\text{J/kg}$ relative to $293.15\text{ K}$).
    pub enthalpy_j_kg: f64,
    /// Temperature $T$ ($\text{K}$).
    pub temperature_k: f64,
    /// Mushy-zone liquid fraction $f_L \in [0, 1]$.
    pub melt_fraction: f64,
    /// Adhered liquid nitrogen ($\text{LN}_2$) film mass ($\text{g}$).
    pub ln2_mass_g: f64,
}

/// Single finite-volume cell in the $24 \times 24$ spatial PDE grid across the watch movement & case.
#[derive(Debug, Clone, Copy)]
pub struct SpatialCell {
    pub x_mm: f64,
    pub y_mm: f64,
    pub r_mm: f64,
    pub active: bool,
    pub is_gold_case: bool,
    pub material: WatchMaterialId,
    pub owner_node: usize,
    /// True 3D movement substrate topography $z_{\text{bed}}(x,y)$ ($\text{mm}$):
    /// stepped brass mainplate ($0.0\text{ mm}$), raised Geneva bridges ($+1.45\text{ mm}$),
    /// gear train recesses ($-0.45\text{ mm}$), and the 6h Tourbillon well ($-1.65\text{ mm}$).
    pub bed_elevation_mm: f64,
    pub thermal_mass_j_per_k: f64,
    pub base_mass_kg: f64,
    pub enthalpy_j_kg: f64,
    pub temperature_k: f64,
    pub melt_fraction: f64,
    /// Local case shell radial/vertical displacement ratio mapped from the 3D CaseShellNode ring.
    pub case_slump_frac: f64,
    /// Mobile liquid gold height $H_{\text{Au}}$ ($\text{mm}$) in this cell.
    pub liquid_gold_mm: f64,
    /// Specific enthalpy of the mobile liquid gold layer ($\text{J/kg}$) for enthalpy advection.
    pub liquid_gold_enthalpy_j_kg: f64,
    /// Solidified (frozen) displaced gold height $S_{\text{Au}}$ ($\text{mm}$) in this cell.
    pub frozen_gold_mm: f64,
    /// Liquid nitrogen ($\text{LN}_2$) film mass ($\text{g}$) in this cell.
    pub ln2_mass_g: f64,
    /// Frost / rime ice thickness ($\mu\text{m}$) in this cell.
    pub frost_um: f64,
    /// Watch oil film thickness ($\mu\text{m}$) in this cell.
    pub oil_um: f64,
}

/// Complete multiphysics state of the pocket watch, backed by `fs_lbm::FreeSurface`
/// (`D2Q9` Körner VOF Lattice Boltzmann with Laplace curvature surface tension),
/// `fs_material::creep::plastic::NortonJ2`, and `fs_tribo::journal::GeneralizedNewtonianJournal`.
#[derive(Clone)]
pub struct MultiphysicsState {
    pub case_shell: Vec<CaseShellNode>,
    pub grid: Vec<SpatialCell>,
    /// Live $24 \times 24$ `D2Q9` Körner VOF Free-Surface Lattice Boltzmann solver (`fs-lbm`)
    /// governing molten gold capillary pooling, surface tension, and momentum transport.
    pub lbm_free_surface: LbmFreeSurface,
    pub nodes: [ThermalPhaseNode; NODE_COUNT],
    pub inputs: PerturbationInputs,
    pub remanent_magnetization_mt: f64,
    pub degauss_envelope_mt: f64,
    pub degauss_phase_rad: f64,
    pub hairspring_oxide_nm: f64,
    pub modifiers: CoupledWatchModifiers,
}

impl MultiphysicsState {
    #[must_use]
    pub fn new() -> Self {
        // 1. Construct the 64 x 8 3D Bassine 18k Gold Case Shell Mesh
        // Meridional profile stations iz = 0..7 from back bezel (z = -3.8 mm) around the
        // rounded Bassine caseband (r = 26.8 mm) to the front exhibition bezel (z = +3.8 mm).
        let profile_z_mm: [f64; CASE_SHELL_N_Z] =
            [-3.8, -2.8, -1.6, -0.5, 0.5, 1.6, 2.8, 3.8];
        let profile_r_mm: [f64; CASE_SHELL_N_Z] =
            [23.0, 25.0, 26.3, 26.8, 26.8, 26.3, 25.0, 23.0];
        let profile_h_mm: [f64; CASE_SHELL_N_Z] =
            [2.1, 2.3, 2.6, 2.8, 2.8, 2.6, 2.3, 2.1];

        let d_theta = TAU / (CASE_SHELL_N_THETA as f64);

        let case_shell: Vec<CaseShellNode> = (0..CASE_SHELL_NODES)
            .map(|idx| {
            let it = idx % CASE_SHELL_N_THETA;
            let iz = idx / CASE_SHELL_N_THETA;
            let theta = (it as f64) * d_theta;

            // At 12 o'clock (theta ≈ pi/2), the Bassine case transitions into the sculpted pendant neck
            let dy_top = theta.sin().max(0.0);
            let pendant_lobe = if (theta - 0.5 * PI).abs() < 0.22 {
                let u = 1.0 - ((theta - 0.5 * PI).abs() / 0.22);
                2.2 * u * u * (1.0 - 0.12 * profile_z_mm[iz].abs())
            } else {
                0.0
            };
            let r0 = profile_r_mm[iz] + pendant_lobe * dy_top;
            let x0 = r0 * theta.cos();
            let y0 = r0 * theta.sin();
            let z0 = profile_z_mm[iz];
            let h0 = profile_h_mm[iz] + 0.35 * pendant_lobe;

            let l_theta_mm = (r0 * d_theta).max(1.5);
            let l_z_mm = if iz + 1 < CASE_SHELL_N_Z {
                let dr = profile_r_mm[iz + 1] - profile_r_mm[iz];
                let dz = profile_z_mm[iz + 1] - profile_z_mm[iz];
                (dr * dr + dz * dz).sqrt()
            } else {
                1.25
            };

            let area_m2 = (l_theta_mm * 1.0e-3) * (l_z_mm * 1.0e-3);
            let vol_m3 = area_m2 * (h0 * 1.0e-3);
            let mass_kg = (GOLD_18K_CARD.density_kg_m3 * vol_m3).max(2.0e-5);

            CaseShellNode {
                theta_rad: theta,
                ref_pos_mm: [x0, y0, z0],
                pos_mm: [x0, y0, z0],
                vel_mm_s: [0.0, 0.0, 0.0],
                ref_thickness_mm: h0,
                thickness_mm: h0,
                ref_l_theta_mm: l_theta_mm,
                ref_l_z_mm: l_z_mm,
                mass_kg,
                area_m2,
                eps_p_theta: 0.0,
                eps_p_z: 0.0,
                plastic_disp_mm: [0.0, 0.0, 0.0],
                eq_plastic_strain: 0.0,
                von_mises_stress_pa: 0.0,
                enthalpy_j_kg: 0.0,
                temperature_k: T_REF_K,
                melt_fraction: 0.0,
                ln2_mass_g: 0.0,
            }
        })
        .collect();

        // 2. Construct the 24 x 24 2D Finite-Volume Free-Surface + Movement Substrate PDE Grid
        let grid: Vec<SpatialCell> = (0..GRID_CELLS)
            .map(|idx| {
            let gx = idx % GRID_N;
            let gy = idx / GRID_N;
            let x_mm = -GRID_HALF_SPAN_MM + (gx as f64 + 0.5) * GRID_DX_MM;
            let y_mm = -GRID_HALF_SPAN_MM + (gy as f64 + 0.5) * GRID_DX_MM;
            let r_mm = (x_mm * x_mm + y_mm * y_mm).sqrt();

            let in_pendant = x_mm.abs() <= 4.8 && y_mm > 22.0 && y_mm <= 28.8;
            let active = r_mm <= 27.6 || in_pendant;
            let is_gold_case = active && (r_mm >= 22.2 || in_pendant);

            // True 3D horological architecture topography z_bed(x,y):
            // - 6h Tourbillon & Escapement well centered at (0, -11.0), radius 6.8 mm: recessed to -1.65 mm
            // - 12h Mainspring Barrel recess centered at (0, +11.5), radius 7.8 mm: +0.25 mm drum top, surrounded by +1.35 mm Barrel Bridge
            // - Going-train wheel channel centered along (-5.5, 1.5): recessed to -0.45 mm
            // - Contoured Geneva bridges on left/right flanks: raised to +1.25 mm
            let dist_tourb = (x_mm * x_mm + (y_mm + 11.0).powi(2)).sqrt();
            let dist_barrel = (x_mm * x_mm + (y_mm - 11.5).powi(2)).sqrt();
            let dist_train = ((x_mm + 5.5).powi(2) + (y_mm - 1.5).powi(2)).sqrt();

            let (material, owner_node, bed_elevation_mm, thickness_m, init_oil_um) = if is_gold_case
            {
                let sector = if y_mm >= x_mm.abs() {
                    0 // 12h Top
                } else if y_mm <= -x_mm.abs() {
                    3 // 6h Bottom
                } else if x_mm < 0.0 {
                    1 // 9h Left
                } else {
                    2 // 3h Right
                };
                (WatchMaterialId::Gold18kCase, sector, 3.6, 3.8e-3, 0.0)
            } else if dist_tourb <= 6.8 {
                let sub_owner = if dist_tourb <= 2.5 {
                    11 // Breguet Hairspring & Balance staff center
                } else if dist_tourb <= 4.6 && y_mm >= -11.5 {
                    10 // Glucydur Balance Wheel rim
                } else if x_mm < -0.6 {
                    9 // Swiss Lever Pallet Fork & Club-Tooth Escape Wheel
                } else {
                    8 // Titanium Tourbillon Carriage
                };
                // Deep recessed Tourbillon well (-1.65 mm at center, stepped to -1.15 mm at rim)
                let well_z = -1.65 + 0.07 * dist_tourb;
                let mat = match sub_owner {
                    9 | 11 => WatchMaterialId::HighCarbonSteel,
                    10 => WatchMaterialId::GlucydurBeCu,
                    _ => WatchMaterialId::TitaniumGrade5,
                };
                (mat, sub_owner, well_z, 1.35e-3, 1.5)
            } else if dist_barrel <= 7.8 {
                // Mainspring barrel drum (+0.20 mm) and surrounding Geneva barrel bridge (+1.10 mm)
                let barrel_z = if dist_barrel <= 5.8 { 0.20 } else { 1.10 };
                (WatchMaterialId::HighCarbonSteel, 6, barrel_z, 2.1e-3, 2.0)
            } else if dist_train <= 6.2 {
                // Recessed going-train wheel channel (-0.45 mm) so molten gold channels toward center & 6h
                (WatchMaterialId::BrassC36000, 7, -0.45, 1.65e-3, 1.6)
            } else {
                // Stepped brass mainplate (0.0 mm) with raised Geneva bridge flanks (+0.95 mm)
                let bridge_z = if x_mm.abs() > 7.5 && y_mm > -6.0 && y_mm < 10.0 {
                    0.95
                } else {
                    0.0
                };
                (WatchMaterialId::BrassC36000, 5, bridge_z, 1.95e-3, 1.2)
            };

            let card = material.card();
            let base_mass_kg = (card.density_kg_m3 * GRID_CELL_AREA_M2 * thickness_m).max(1.2e-5);
            let thermal_mass_j_per_k = base_mass_kg * card.cp_j_kg_k;

            SpatialCell {
                x_mm,
                y_mm,
                r_mm,
                active,
                is_gold_case,
                material,
                owner_node,
                bed_elevation_mm,
                thermal_mass_j_per_k,
                base_mass_kg,
                enthalpy_j_kg: 0.0,
                temperature_k: T_REF_K,
                melt_fraction: 0.0,
                case_slump_frac: 0.0,
                liquid_gold_mm: 0.0,
                liquid_gold_enthalpy_j_kg: GOLD_18K_CARD
                    .enthalpy_from_temperature_j_kg(GOLD_18K_CARD.t_liquidus_k + 15.0),
                frozen_gold_mm: 0.0,
                ln2_mass_g: 0.0,
                frost_um: 0.0,
                oil_um: if active && !is_gold_case {
                    init_oil_um
                } else {
                    0.0
                },
            }
        })
        .collect();

        let make_node = |id: WatchNodeId,
                         material: WatchMaterialId,
                         pos_mm: [f64; 3],
                         base_mass_kg: f64,
                         surface_area_m2: f64,
                         oil_film_um: f64,
                         oxide_thickness_nm: f64| ThermalPhaseNode {
            id,
            material,
            pos_mm,
            base_mass_kg,
            surface_area_m2,
            enthalpy_j_kg: 0.0,
            temperature_k: T_REF_K,
            melt_fraction: 0.0,
            plastic_slump: 0.0,
            molten_gold_g: 0.0,
            frozen_gold_g: 0.0,
            ln2_film_g: 0.0,
            frost_thickness_um: 0.0,
            oil_film_um,
            oxide_thickness_nm,
            incandescence_rgb: [0.0, 0.0, 0.0],
        };

        let nodes = [
            make_node(
                WatchNodeId::CaseTop18kGold,
                WatchMaterialId::Gold18kCase,
                [0.0, 25.0, 0.0],
                0.0115,
                5.2e-4,
                0.0,
                0.0,
            ),
            make_node(
                WatchNodeId::CaseLeft18kGold,
                WatchMaterialId::Gold18kCase,
                [-25.0, 0.0, 0.0],
                0.0105,
                4.8e-4,
                0.0,
                0.0,
            ),
            make_node(
                WatchNodeId::CaseRight18kGold,
                WatchMaterialId::Gold18kCase,
                [25.0, 0.0, 0.0],
                0.0105,
                4.8e-4,
                0.0,
                0.0,
            ),
            make_node(
                WatchNodeId::CaseBottom18kGold,
                WatchMaterialId::Gold18kCase,
                [0.0, -25.0, 0.0],
                0.0115,
                5.2e-4,
                0.0,
                0.0,
            ),
            make_node(
                WatchNodeId::CrystalSapphire,
                WatchMaterialId::SyntheticCorundumRuby,
                [0.0, 0.0, 5.5],
                0.0062,
                14.5e-4,
                0.0,
                0.0,
            ),
            make_node(
                WatchNodeId::MainplateBrass,
                WatchMaterialId::BrassC36000,
                [0.0, 0.0, -1.5],
                0.0184,
                13.0e-4,
                1.5,
                0.0,
            ),
            make_node(
                WatchNodeId::BarrelMainspring,
                WatchMaterialId::HighCarbonSteel,
                [0.0, 11.5, 0.5],
                0.0054,
                3.1e-4,
                2.0,
                15.0,
            ),
            make_node(
                WatchNodeId::GoingTrainWheels,
                WatchMaterialId::BrassC36000,
                [-5.5, 1.5, 0.8],
                0.0021,
                2.0e-4,
                1.6,
                0.0,
            ),
            make_node(
                WatchNodeId::TourbillonCage,
                WatchMaterialId::TitaniumGrade5,
                [0.0, -11.0, 1.2],
                0.00042,
                1.1e-4,
                1.4,
                0.0,
            ),
            make_node(
                WatchNodeId::EscapementForkWheel,
                WatchMaterialId::HighCarbonSteel,
                [-2.2, -11.0, 1.2],
                0.00018,
                0.55e-4,
                1.6,
                25.0,
            ),
            make_node(
                WatchNodeId::BalanceWheelRim,
                WatchMaterialId::GlucydurBeCu,
                [0.0, -11.0, 1.8],
                0.00036,
                0.95e-4,
                1.5,
                0.0,
            ),
            make_node(
                WatchNodeId::HairspringBreguet,
                WatchMaterialId::HighCarbonSteel,
                [0.0, -11.0, 2.5],
                0.000035,
                0.45e-4,
                0.0,
                72.0,
            ),
        ];

        let mut lbm_grid = LbmGrid::uniform(GRID_N, GRID_N, 0.85);
        for idx in 0..GRID_CELLS {
            lbm_grid.flags[idx] = if !grid[idx].active {
                LbmCell::Wall
            } else {
                LbmCell::Gas
            };
        }
        let lbm_free_surface = LbmFreeSurface::new(lbm_grid, 0.008, LbmContactModel::Neutral);

        let mut state = Self {
            case_shell,
            grid,
            lbm_free_surface,
            nodes,
            inputs: PerturbationInputs::default(),
            remanent_magnetization_mt: 0.0,
            degauss_envelope_mt: 0.0,
            degauss_phase_rad: 0.0,
            hairspring_oxide_nm: 72.0,
            modifiers: CoupledWatchModifiers {
                balance_inertia_ratio: 1.0,
                hairspring_stiffness_ratio: 1.0,
                barrel_youngs_modulus_ratio: 1.0,
                natural_freq_ratio: 1.0,
                theoretical_rate_s_per_day: 0.0,
                pivot_friction_multiplier: 1.0,
                extra_balance_damping_s_inv: 0.0,
                hairspring_oil_bridge_fraction: 0.0,
                min_pivot_clearance_um: 4.5,
                b_field_escapement_mt: 0.0,
                remanent_magnetization_mt: 0.0,
                pallet_magnetic_bias_torque_nm: 0.0,
                hairspring_eq_bias_rad: 0.0,
                gold_fouling_jam_factor: 0.0,
                liquid_gold_viscous_torque_nm: 0.0,
                solid_gold_yield_torque_nm: 0.0,
                frozen_gold_locked: false,
                thermal_or_cryo_seized: false,
                magnetically_pinned: false,
                ln2_boiling_power_w: 0.0,
                blowtorch_heating_power_w: 0.0,
                blued_steel_rgb: blued_steel_interference_rgb(72.0, 1.0),
                rate_sensitivity_s_day_per_k: -0.45,
                rate_interval_lo_s_day: -0.5,
                rate_interval_hi_s_day: 0.5,
                balance_aero_drag_s_inv: 0.0,
            },
        };
        state.aggregate_nodes_from_continuum();
        state.recompute_modifiers();
        state
    }

    pub fn trigger_degauss(&mut self) {
        self.degauss_envelope_mt = 85.0;
        self.degauss_phase_rad = 0.0;
    }

    pub fn ultrasonic_clean_and_relubricate(&mut self) {
        for cell in &mut self.grid {
            if cell.active && !cell.is_gold_case {
                cell.oil_um = 1.6;
            } else {
                cell.oil_um = 0.0;
            }
        }
        self.aggregate_nodes_from_continuum();
        self.recompute_modifiers();
    }

    /// Advance the coupled 3D Thermo-Elasto-Viscoplastic Case Shell and 2D Free-Surface PDE Grid.
    pub fn step(&mut self, dt_s: f64) {
        let dt = dt_s.clamp(1.0e-5, 0.05);
        let substeps = ((dt / 0.0025).ceil() as usize).clamp(1, 20);
        let h_dt = dt / (substeps as f64);

        let hs_mat = if self.inputs.compensated_hairspring {
            WatchMaterialId::Invar36Ni
        } else {
            WatchMaterialId::HighCarbonSteel
        };
        self.nodes[WatchNodeId::HairspringBreguet as usize].material = hs_mat;

        let mut last_torch_w = 0.0;
        let mut last_ln2_w = 0.0;

        for _ in 0..substeps {
            let (torch_shell_w, ln2_shell_w) = self.step_case_shell_continuum(h_dt);
            let (torch_grid_w, ln2_grid_w) = self.step_spatial_pde_grid(h_dt);
            last_torch_w = torch_shell_w + torch_grid_w;
            last_ln2_w = ln2_shell_w + ln2_grid_w;
        }

        if self.degauss_envelope_mt > 0.05 {
            self.degauss_phase_rad += TAU * 60.0 * dt;
            self.degauss_envelope_mt *= (-4.2 * dt).exp();
            self.remanent_magnetization_mt *= (-6.5 * dt).exp();
            if self.degauss_envelope_mt < 0.08 {
                self.degauss_envelope_mt = 0.0;
                self.remanent_magnetization_mt = 0.0;
            }
        }

        self.modifiers.blowtorch_heating_power_w = last_torch_w;
        self.modifiers.ln2_boiling_power_w = last_ln2_w;
        self.aggregate_nodes_from_continuum();
        self.recompute_modifiers();
    }

    /// Integrate the 3D Thermo-Elasto-Viscoplastic + Surface-Tension Case Shell (`64 x 8` nodes).
    /// Returns `(blowtorch_power_w, ln2_boiling_power_w)` absorbed by the case shell.
    fn step_case_shell_continuum(&mut self, dt: f64) -> (f64, f64) {
        let mut q_cond_w = [0.0_f64; CASE_SHELL_NODES];
        let mut force_n = [[0.0_f64; 3]; CASE_SHELL_NODES];
        let mut torch_total_w = 0.0;
        let mut ln2_total_w = 0.0;

        // 1. Evaluate 3D shell strains, Cauchy stresses, Prandtl-Reuss viscoplastic flow,
        //    Young-Laplace capillary curvature forces, Marangoni tangential shear, and Fourier heat conduction.
        for iz in 0..CASE_SHELL_N_Z {
            for it in 0..CASE_SHELL_N_THETA {
                let idx = iz * CASE_SHELL_N_THETA + it;
                let it_next = (it + 1) % CASE_SHELL_N_THETA;
                let it_prev = (it + CASE_SHELL_N_THETA - 1) % CASE_SHELL_N_THETA;
                let idx_theta_next = iz * CASE_SHELL_N_THETA + it_next;
                let idx_theta_prev = iz * CASE_SHELL_N_THETA + it_prev;

                let node = self.case_shell[idx];
                let n_next = self.case_shell[idx_theta_next];
                let n_prev = self.case_shell[idx_theta_prev];

                // Circumferential vector & stretch
                let dx_t = n_next.pos_mm[0] - node.pos_mm[0];
                let dy_t = n_next.pos_mm[1] - node.pos_mm[1];
                let dz_t = n_next.pos_mm[2] - node.pos_mm[2];
                let l_theta = (dx_t * dx_t + dy_t * dy_t + dz_t * dz_t).sqrt().max(0.2);
                let u_theta = [dx_t / l_theta, dy_t / l_theta, dz_t / l_theta];
                let eps_tot_theta = (l_theta - node.ref_l_theta_mm) / node.ref_l_theta_mm;

                // Circumferential heat conduction along the 18k gold ring
                let a_cross_theta_m2 = (node.thickness_mm * 1.0e-3) * (node.ref_l_z_mm * 1.0e-3);
                let g_theta_w_k = GOLD_18K_CARD.conductivity_w_m_k * a_cross_theta_m2
                    / (node.ref_l_theta_mm * 1.0e-3);
                let q_theta = g_theta_w_k * (n_next.temperature_k - node.temperature_k);
                q_cond_w[idx] += q_theta;
                q_cond_w[idx_theta_next] -= q_theta;

                // Meridional vector & stretch (if iz + 1 < CASE_SHELL_N_Z)
                let (eps_tot_z, u_z, has_z_next, idx_z_next) = if iz + 1 < CASE_SHELL_N_Z {
                    let izn = (iz + 1) * CASE_SHELL_N_THETA + it;
                    let nz = self.case_shell[izn];
                    let dx_z = nz.pos_mm[0] - node.pos_mm[0];
                    let dy_z = nz.pos_mm[1] - node.pos_mm[1];
                    let dz_z = nz.pos_mm[2] - node.pos_mm[2];
                    let l_z = (dx_z * dx_z + dy_z * dy_z + dz_z * dz_z).sqrt().max(0.2);
                    let uz = [dx_z / l_z, dy_z / l_z, dz_z / l_z];
                    let eps_z = (l_z - node.ref_l_z_mm) / node.ref_l_z_mm;

                    let a_cross_z_m2 =
                        (node.thickness_mm * 1.0e-3) * (node.ref_l_theta_mm * 1.0e-3);
                    let g_z_w_k =
                        GOLD_18K_CARD.conductivity_w_m_k * a_cross_z_m2 / (node.ref_l_z_mm * 1.0e-3);
                    let q_z = g_z_w_k * (nz.temperature_k - node.temperature_k);
                    q_cond_w[idx] += q_z;
                    q_cond_w[izn] -= q_z;

                    (eps_z, uz, true, izn)
                } else {
                    (0.0, [0.0, 0.0, 1.0], false, idx)
                };

                // Thermal expansion strain of 18k Gold at current node temperature
                let eps_th = GOLD_18K_CARD.thermal_strain(node.temperature_k);

                // Elastic membrane strains
                let eps_e_theta = eps_tot_theta - eps_th - node.eps_p_theta;
                let eps_e_z = eps_tot_z - eps_th - node.eps_p_z;

                // Temperature- and melt-dependent Young's modulus E(T, f_L)
                // Solid skeleton stiffness decays as (1 - f_L)^3 across the mushy zone
                let solid_frac = (1.0 - node.melt_fraction).clamp(0.0, 1.0);
                let e_eff_pa = GOLD_18K_CARD.youngs_modulus_pa
                    * GOLD_18K_CARD.youngs_modulus_ratio(node.temperature_k)
                    * solid_frac.powi(3);
                let nu = 0.42; // Poisson's ratio of Gold
                let plane_stress_factor = e_eff_pa / (1.0 - nu * nu);

                let sigma_theta_pa = plane_stress_factor * (eps_e_theta + nu * eps_e_z);
                let sigma_z_pa = plane_stress_factor * (eps_e_z + nu * eps_e_theta);

                // Elastic shell bending restoration toward (ref_pos_mm + thermal expansion + plastic_disp_mm)
                let r_th_scale = 1.0 + eps_th;
                let eq_pos_mm = [
                    node.ref_pos_mm[0] * r_th_scale + node.plastic_disp_mm[0],
                    node.ref_pos_mm[1] * r_th_scale + node.plastic_disp_mm[1],
                    node.ref_pos_mm[2] * r_th_scale + node.plastic_disp_mm[2],
                ];
                let bend_disp_mm = [
                    node.pos_mm[0] - eq_pos_mm[0],
                    node.pos_mm[1] - eq_pos_mm[1],
                    node.pos_mm[2] - eq_pos_mm[2],
                ];
                let bend_dist_mm = (bend_disp_mm[0].powi(2)
                    + bend_disp_mm[1].powi(2)
                    + bend_disp_mm[2].powi(2))
                .sqrt();
                // Bending stiffness D ~ E * h^3 / (12 * L^4)
                let h_ratio = (node.thickness_mm / node.ref_thickness_mm).clamp(0.15, 2.0);
                let sigma_bend_pa =
                    0.015 * e_eff_pa * h_ratio * (bend_dist_mm / node.ref_thickness_mm.max(1.0));

                // Von Mises equivalent stress in the shell element
                let sigma_eq_pa = (sigma_theta_pa * sigma_theta_pa
                    - sigma_theta_pa * sigma_z_pa
                    + sigma_z_pa * sigma_z_pa
                    + 3.0 * sigma_bend_pa * sigma_bend_pa)
                    .max(0.0)
                    .sqrt();

                // Update Prandtl-Reuss viscoplastic strains using `fs_material::creep::plastic::NortonJ2`
                // (simultaneous secondary Norton creep + Simo-Hughes radial-return J2 plasticity)
                // coupled with the high-strain mushy-zone Johnson-Cook flow rule.
                let mut norton_j2_dp = 0.0;
                if node.temperature_k > 350.0
                    && let Some(law) = gold_coupled_norton_j2(node.temperature_k, node.melt_fraction, dt)
                {
                    let e_t = eps_e_theta.clamp(-0.02, 0.02);
                    let e_z = eps_e_z.clamp(-0.02, 0.02);
                    let strain_voigt = [e_t, e_z, -0.42 * (e_t + e_z), 0.0, 0.0, 0.0];
                    let old_state = NortonJ2State {
                        creep: CreepState {
                            creep_strain: [0.0; 6],
                            equivalent_creep_strain: 0.0,
                        },
                        plastic: J2State {
                            plastic_strain: [0.0; 6],
                            alpha: node.eq_plastic_strain.clamp(0.0, 2.0),
                        },
                    };
                    if let Ok(resp) = law.response(&strain_voigt, &old_state) {
                        norton_j2_dp = resp.state.creep.equivalent_creep_strain
                            + (resp.state.plastic.alpha - old_state.plastic.alpha).max(0.0);
                    }
                }
                let vp_rate = gold_viscoplastic_strain_rate_s_inv(
                    sigma_eq_pa,
                    node.temperature_k,
                    node.melt_fraction,
                    node.eq_plastic_strain,
                );
                let d_eq_p = (vp_rate * dt + norton_j2_dp).min(0.08);
                if d_eq_p > 1.0e-9 && sigma_eq_pa > 1.0 {
                    let s_theta = (2.0 * sigma_theta_pa - sigma_z_pa) / 3.0;
                    let s_z = (2.0 * sigma_z_pa - sigma_theta_pa) / 3.0;
                    let max_duct = gold_ductility_elongation(node.temperature_k);
                    self.case_shell[idx].eps_p_theta = (self.case_shell[idx].eps_p_theta
                        + d_eq_p * (1.5 * s_theta / sigma_eq_pa))
                        .clamp(-max_duct, max_duct);
                    self.case_shell[idx].eps_p_z = (self.case_shell[idx].eps_p_z
                        + d_eq_p * (1.5 * s_z / sigma_eq_pa))
                        .clamp(-max_duct, max_duct);
                    // Viscoplastic relaxation of the 3D shell bending offset toward the current deformed configuration
                    let relax = (d_eq_p * 8.5).clamp(0.0, 0.65);
                    for c in 0..3 {
                        self.case_shell[idx].plastic_disp_mm[c] += relax * bend_disp_mm[c];
                    }
                    self.case_shell[idx].eq_plastic_strain += d_eq_p;
                }
                self.case_shell[idx].von_mises_stress_pa = sigma_eq_pa;

                // Convert membrane stresses into nodal forces (scaled to keep high-stiffness explicit stable)
                let f_theta_mag =
                    (sigma_theta_pa * a_cross_theta_m2).clamp(-18.0, 18.0) * 0.012;
                for c in 0..3 {
                    let df = f_theta_mag * u_theta[c];
                    force_n[idx][c] += df;
                    force_n[idx_theta_next][c] -= df;
                }
                if has_z_next {
                    let a_cross_z_m2 =
                        (node.thickness_mm * 1.0e-3) * (node.ref_l_theta_mm * 1.0e-3);
                    let f_z_mag = (sigma_z_pa * a_cross_z_m2).clamp(-18.0, 18.0) * 0.012;
                    for c in 0..3 {
                        let df = f_z_mag * u_z[c];
                        force_n[idx][c] += df;
                        force_n[idx_z_next][c] -= df;
                    }
                }

                // Elastic bending restoring force toward eq_pos_mm
                let k_bend_n_per_mm = 0.45 * solid_frac * solid_frac * h_ratio * h_ratio;
                for c in 0..3 {
                    force_n[idx][c] -= k_bend_n_per_mm * bend_disp_mm[c];
                }

                // Discrete Mean-Curvature Vector for Young-Laplace Surface Tension:
                // F_cap = sigma_Au(T) * A_node * Laplacian_s(x)
                let lap_theta = [
                    n_prev.pos_mm[0] - 2.0 * node.pos_mm[0] + n_next.pos_mm[0],
                    n_prev.pos_mm[1] - 2.0 * node.pos_mm[1] + n_next.pos_mm[1],
                    n_prev.pos_mm[2] - 2.0 * node.pos_mm[2] + n_next.pos_mm[2],
                ];
                let lap_z = if iz > 0 && iz + 1 < CASE_SHELL_N_Z {
                    let n_zm = self.case_shell[(iz - 1) * CASE_SHELL_N_THETA + it];
                    let n_zp = self.case_shell[(iz + 1) * CASE_SHELL_N_THETA + it];
                    [
                        n_zm.pos_mm[0] - 2.0 * node.pos_mm[0] + n_zp.pos_mm[0],
                        n_zm.pos_mm[1] - 2.0 * node.pos_mm[1] + n_zp.pos_mm[1],
                        n_zm.pos_mm[2] - 2.0 * node.pos_mm[2] + n_zp.pos_mm[2],
                    ]
                } else {
                    // Free bezel rim retraction driven by Young-Laplace curvature at top/bottom bezel edge
                    let iz_adj = if iz == 0 { 1 } else { CASE_SHELL_N_Z - 2 };
                    let n_adj = self.case_shell[iz_adj * CASE_SHELL_N_THETA + it];
                    [
                        n_adj.pos_mm[0] - node.pos_mm[0],
                        n_adj.pos_mm[1] - node.pos_mm[1],
                        n_adj.pos_mm[2] - node.pos_mm[2],
                    ]
                };

                // When softening / melting (T > 900 K or melt_fraction > 0), surface tension & Marangoni
                // forces act on the free liquid/mushy surface!
                let fluidity_weight =
                    (node.melt_fraction + 0.15 * ((node.temperature_k - 900.0) / 300.0).clamp(0.0, 1.0))
                        .clamp(0.0, 1.0);
                if fluidity_weight > 1.0e-3 {
                    let sigma_au = gold_liquid_surface_tension_n_m(node.temperature_k);
                    let cap_coeff = 0.038 * sigma_au * fluidity_weight;
                    for c in 0..3 {
                        force_n[idx][c] += cap_coeff * (lap_theta[c] + 1.35 * lap_z[c]);
                    }

                    // Marangoni thermocapillary tangential force: d(sigma)/dT * grad_s(T)
                    // Pulls liquid gold away from the hottest torch center toward the cooler surrounding rim
                    // (producing the characteristic raised solidified crater rim around a melt pool!)
                    let dt_theta = n_next.temperature_k - n_prev.temperature_k;
                    let marangoni_f = 1.2e-4
                        * GOLD_MARANGONI_DSIGMA_DT_N_M_K
                        * dt_theta
                        * fluidity_weight;
                    for c in 0..3 {
                        force_n[idx][c] += marangoni_f * u_theta[c];
                    }
                }
            }
        }

        // 2. External thermal & mechanical loads on each 3D shell node, enthalpy update, and 3D position integration
        let torch_active = self.inputs.blowtorch_power > 1.0e-4;
        let ln2_active = self.inputs.ln2_pour_rate > 1.0e-4;

        for iz in 0..CASE_SHELL_N_Z {
            for it in 0..CASE_SHELL_N_THETA {
                let idx = iz * CASE_SHELL_N_THETA + it;
                let mut node = self.case_shell[idx];

                // Couple shell node thermally to nearest interior 24x24 movement cell
                let gx = (((node.pos_mm[0] + GRID_HALF_SPAN_MM) / GRID_DX_MM).floor() as isize)
                    .clamp(0, (GRID_N - 1) as isize) as usize;
                let gy = (((node.pos_mm[1] + GRID_HALF_SPAN_MM) / GRID_DX_MM).floor() as isize)
                    .clamp(0, (GRID_N - 1) as isize) as usize;
                let cell_idx = gy * GRID_N + gx;
                let t_cell = self.grid[cell_idx].temperature_k;
                let q_to_movement = 0.0028 * (node.temperature_k - t_cell);

                // Ambient natural convection (`fs_convection::CorrelationId::ChurchillChuHorizontalCylinder`)
                // + Stefan-Boltzmann thermal radiation
                let dt_film = (node.temperature_k - T_REF_K).abs().max(0.5);
                let t_film = 0.5 * (node.temperature_k + T_REF_K);
                let d_case_m = 0.056; // 56 mm Bassine case diameter
                let ra_case = (9.81 * (1.0 / t_film) * dt_film * d_case_m.powi(3)
                    / (1.58e-5 * 2.25e-5))
                    .clamp(1.0e-4, 1.0e11);
                let h_nat_w_m2_k = eval_convection(
                    CorrelationId::ChurchillChuHorizontalCylinder,
                    CorrelationInputs::natural(ra_case, 0.71),
                )
                .and_then(|nu| {
                    let _ = NaturalConvectionGeometry::HorizontalCylinder;
                    nu.heat_transfer_coefficient(
                        ThermalConductivity::new(0.0263),
                        ConvLength::new(d_case_m),
                    )
                })
                .map(|ev| ev.value.value().clamp(6.5, 45.0))
                .unwrap_or(14.0);
                let q_amb = h_nat_w_m2_k * node.area_m2 * (T_REF_K - node.temperature_k)
                    + GOLD_18K_CARD.emissivity
                        * STEFAN_BOLTZMANN
                        * node.area_m2
                        * (T_REF_K.powi(4) - node.temperature_k.powi(4));

                let mut q_ext_w = q_cond_w[idx] + q_amb - q_to_movement;
                self.grid[cell_idx].enthalpy_j_kg +=
                    (q_to_movement * dt) / self.grid[cell_idx].base_mass_kg.max(1.0e-5);

                // Outward radial unit vector in XY plane
                let r_xy = (node.pos_mm[0].powi(2) + node.pos_mm[1].powi(2))
                    .sqrt()
                    .max(1.0);
                let n_rad = [node.pos_mm[0] / r_xy, node.pos_mm[1] / r_xy, 0.0];

                // Virtual oxy-butane blowtorch (T_flame = 2250 K) heat flux & stagnation jet pressure
                // with `fs_convection::CorrelationId::ChurchillBernsteinCylinder` forced crossflow amplification
                if torch_active {
                    let dx = node.pos_mm[0] - self.inputs.blowtorch_pos_mm[0];
                    let dy = node.pos_mm[1] - self.inputs.blowtorch_pos_mm[1];
                    let dz = node.pos_mm[2] - self.inputs.blowtorch_pos_mm[2];
                    let d2 = dx * dx + dy * dy + 0.7 * dz * dz;
                    let w_torch = (-d2 / (2.0 * 9.5 * 9.5)).exp();
                    let flame_temp_k = 2250.0;
                    if flame_temp_k > node.temperature_k && w_torch > 1.0e-4 {
                        let re_jet = (1800.0 * self.inputs.blowtorch_power * w_torch).clamp(100.0, 2.5e5);
                        let cb_boost = eval_convection(
                            CorrelationId::ChurchillBernsteinCylinder,
                            CorrelationInputs::forced(re_jet, 0.72),
                        )
                        .map(|ev| (ev.evidence().value / 24.0).clamp(0.85, 1.45))
                        .unwrap_or(1.0);
                        let q_t = self.inputs.blowtorch_power
                            * w_torch
                            * 0.0185
                            * cb_boost
                            * (flame_temp_k - node.temperature_k);
                        q_ext_w += q_t;
                        torch_total_w += q_t;

                        // Blowtorch gas-jet stagnation pressure pushes softened/molten gold inward
                        // toward the movement interior along (-n_rad) and downward along gravity!
                        let theta_star = ((node.temperature_k - 750.0) / 550.0).clamp(0.0, 1.0);
                        let softness =
                            (0.20 * theta_star * theta_star + 0.80 * node.melt_fraction).clamp(0.0, 1.0);
                        let f_jet = self.inputs.blowtorch_power * w_torch * softness * 0.055;
                        force_n[idx][0] -= f_jet * n_rad[0];
                        force_n[idx][1] -= f_jet * n_rad[1];
                        force_n[idx][2] -= f_jet * 0.35 * (node.pos_mm[2] / 3.8);
                    }
                }

                // Gravity component acting on viscoplastic / molten gold (watch held at slight bench angle)
                if node.temperature_k > 750.0 || node.melt_fraction > 0.01 {
                    let softness = (node.melt_fraction
                        + 0.12 * ((node.temperature_k - 750.0) / 550.0).clamp(0.0, 1.0))
                    .clamp(0.0, 1.0);
                    // Gravity pulls softened/molten gold downward (-Y and -Z) and inward over the bezel ledge
                    force_n[idx][0] -= 0.012 * softness * n_rad[0];
                    force_n[idx][1] -= 0.018 * softness * n_rad[1] + 0.008 * softness;
                    force_n[idx][2] -= 0.014 * softness * (node.pos_mm[2] + 1.2).signum();
                }

                // Cryogenic Liquid Nitrogen (LN2 at 77.36 K) impingement & Nukiyama boiling curve
                if ln2_active {
                    let dx = node.pos_mm[0] - self.inputs.ln2_pos_mm[0];
                    let dy = node.pos_mm[1] - self.inputs.ln2_pos_mm[1];
                    let d2 = dx * dx + dy * dy;
                    let w_ln2 = (-d2 / (2.0 * 15.0 * 15.0)).exp();
                    node.ln2_mass_g = (node.ln2_mass_g
                        + self.inputs.ln2_pour_rate * w_ln2 * 0.18 * dt)
                        .min(0.25);
                }
                if node.ln2_mass_g > 1.0e-6 {
                    if node.temperature_k > T_LN2_SAT_K + 0.05 {
                        let dt_sat = node.temperature_k - T_LN2_SAT_K;
                        let h_boil = if node.temperature_k <= T_LN2_LEIDENFROST_K {
                            4_800.0 + 6_200.0 * (dt_sat / 48.0).clamp(0.1, 1.0)
                        } else {
                            2_100.0 + 2.4 * (node.temperature_k - T_LN2_LEIDENFROST_K)
                        };
                        let q_boil_desire = h_boil * node.area_m2 * dt_sat * 2.2;
                        let max_latent_j = node.ln2_mass_g * 1.0e-3 * LN2_LATENT_VAP_J_PER_KG;
                        let q_boil_j = (q_boil_desire * dt).min(max_latent_j);
                        let q_boil_w = q_boil_j / dt;
                        q_ext_w -= q_boil_w;
                        ln2_total_w += q_boil_w;
                        node.ln2_mass_g = (node.ln2_mass_g
                            - (q_boil_j / LN2_LATENT_VAP_J_PER_KG) * 1.0e3)
                            .max(0.0);
                    } else {
                        node.ln2_mass_g = (node.ln2_mass_g - 0.01 * dt).max(0.0);
                    }
                }

                // Integrate nodal specific enthalpy H and invert enthalpy -> (T, f_L)
                node.enthalpy_j_kg += (q_ext_w * dt) / node.mass_kg.max(1.0e-5);
                let (t_new, phi_new) =
                    GOLD_18K_CARD.temperature_and_melt_from_enthalpy(node.enthalpy_j_kg);
                node.temperature_k = t_new;
                node.melt_fraction = phi_new;

                // Integrate 3D nodal velocity & position with temperature/viscosity-dependent damping:
                // When solid at low temperature, high structural damping locks the shell;
                // when viscoplastic or molten, viscosity mu_eff(T, f_L) governs flow speed!
                let mu_eff = gold_effective_viscosity_pa_s(t_new, phi_new);
                let mobility = if t_new < 600.0 && phi_new <= 0.0 {
                    // Pure thermo-elastic regime (allows instantaneous thermal expansion/contraction)
                    0.35
                } else {
                    (1.0 / (1.0 + 0.008 * mu_eff)).clamp(0.04, 1.0)
                };

                let damping = (-14.0 * dt).exp();
                for c in 0..3 {
                    let acc_mm_s2 = (force_n[idx][c] / node.mass_kg.max(1.0e-4)) * mobility;
                    node.vel_mm_s[c] = (node.vel_mm_s[c] + acc_mm_s2 * dt) * damping;
                    node.pos_mm[c] += node.vel_mm_s[c] * dt;
                }

                // Enforce physical obstacle contact against the movement mainplate interior (r >= 17.2 mm)
                // and table plane (z >= -4.1 mm) so the melting gold case drapes over the movement edge
                let r_curr = (node.pos_mm[0].powi(2) + node.pos_mm[1].powi(2))
                    .sqrt()
                    .max(0.1);
                let r_min_obstacle = 16.8;
                if r_curr < r_min_obstacle {
                    let scale = r_min_obstacle / r_curr;
                    node.pos_mm[0] *= scale;
                    node.pos_mm[1] *= scale;
                    node.vel_mm_s[0] = 0.0;
                    node.vel_mm_s[1] = 0.0;
                }
                node.pos_mm[2] = node.pos_mm[2].clamp(-4.2, 4.4);

                // Update local wall thickness by mass conservation + capillary bead thickening at rim,
                // and shed detached liquid gold droplets onto the 24x24 free-surface PDE grid when f_L > 0.08!
                let r0 = (node.ref_pos_mm[0].powi(2) + node.ref_pos_mm[1].powi(2))
                    .sqrt()
                    .max(1.0);
                let radial_inward_mm = (r0 - r_curr).max(0.0);
                if phi_new > 0.06 && node.thickness_mm > 0.55 {
                    let ck_mob = crate::materials::MaterialCard::carman_kozeny_mobility(phi_new);
                    let drip_rate_mm_s = 0.42 * (0.30 * phi_new + 0.70 * ck_mob);
                    let dh_drip = (drip_rate_mm_s * dt).min(node.thickness_mm - 0.55);
                    if dh_drip > 1.0e-6 {
                        node.thickness_mm -= dh_drip;
                        let dm_kg = GOLD_18K_CARD.density_kg_m3 * node.area_m2 * (dh_drip * 1.0e-3);
                        node.mass_kg = (node.mass_kg - dm_kg).max(5.0e-6);

                        // Deposit detached liquid gold mass & enthalpy into the interior cell just inside the case wall!
                        let inward_x_mm = node.pos_mm[0] * (20.8 / r_curr.max(20.8));
                        let inward_y_mm = node.pos_mm[1] * (20.8 / r_curr.max(20.8));
                        let ix = (((inward_x_mm + GRID_HALF_SPAN_MM) / GRID_DX_MM).floor() as isize)
                            .clamp(0, (GRID_N - 1) as isize)
                            as usize;
                        let iy = (((inward_y_mm + GRID_HALF_SPAN_MM) / GRID_DX_MM).floor() as isize)
                            .clamp(0, (GRID_N - 1) as isize)
                            as usize;
                        let cidx = iy * GRID_N + ix;
                        if self.grid[cidx].active {
                            let added_height_mm =
                                (dm_kg / (GOLD_18K_CARD.density_kg_m3 * GRID_CELL_AREA_M2)) * 1.0e3;
                            let prev_h = self.grid[cidx].liquid_gold_mm;
                            let new_h = prev_h + added_height_mm;
                            if new_h > 1.0e-6 {
                                self.grid[cidx].liquid_gold_enthalpy_j_kg = (prev_h
                                    * self.grid[cidx].liquid_gold_enthalpy_j_kg
                                    + added_height_mm * node.enthalpy_j_kg)
                                    / new_h;
                            }
                            self.grid[cidx].liquid_gold_mm = new_h;
                            // Thermal conduction from hot liquid gold droplet into the substrate cell
                            let dh_enth = (node.enthalpy_j_kg - self.grid[cidx].enthalpy_j_kg).max(0.0);
                            self.grid[cidx].enthalpy_j_kg +=
                                (dm_kg * dh_enth) / (self.grid[cidx].base_mass_kg + dm_kg);
                        }
                    }
                }

                // Map emergent radial/vertical displacement onto the corresponding perimeter cell for telemetry
                let disp_norm = ((radial_inward_mm / 6.5)
                    + ((node.ref_thickness_mm - node.thickness_mm) / node.ref_thickness_mm))
                    .clamp(0.0, 0.95);
                if disp_norm > self.grid[cell_idx].case_slump_frac {
                    self.grid[cell_idx].case_slump_frac = disp_norm;
                }

                self.case_shell[idx] = node;
            }
        }

        (torch_total_w, ln2_total_w)
    }

    /// 2D Finite-Volume Heat Conduction + Enthalpy Phase-Change + Lubrication/Capillary/Marangoni
    /// Free-Surface Molten Gold Flow PDE step over the $24 \times 24$ grid.
    fn step_spatial_pde_grid(&mut self, dt_s: f64) -> (f64, f64) {
        let mut net_q_w = [0.0_f64; GRID_CELLS];
        let mut new_liq_mm = [0.0_f64; GRID_CELLS];
        let mut new_liq_enth = [0.0_f64; GRID_CELLS];
        let mut free_surface_eta_mm = [0.0_f64; GRID_CELLS];
        let mut eta_laplacian_mm = [0.0_f64; GRID_CELLS];
        let mut new_oil_um = [0.0_f64; GRID_CELLS];

        let mut torch_grid_w = 0.0;
        let mut ln2_grid_w = 0.0;

        // 1. Precompute free-surface height eta = z_bed + S_Au + H_Au and sync the
        //    D2Q9 Körner VOF Lattice Boltzmann free-surface grid (`fs_lbm::FreeSurface`).
        let mut any_lbm_wet = false;
        for idx in 0..GRID_CELLS {
            let c = &self.grid[idx];
            new_liq_mm[idx] = c.liquid_gold_mm;
            new_liq_enth[idx] = c.liquid_gold_enthalpy_j_kg;
            new_oil_um[idx] = c.oil_um;
            free_surface_eta_mm[idx] = c.bed_elevation_mm * (1.0 - 0.45 * c.case_slump_frac)
                + c.frozen_gold_mm
                + c.liquid_gold_mm;

            if !c.active {
                self.lbm_free_surface.grid.flags[idx] = LbmCell::Wall;
                self.lbm_free_surface.mass[idx] = 0.0;
            } else if c.liquid_gold_mm > 0.015 {
                any_lbm_wet = true;
                let fill = (c.liquid_gold_mm / 1.2).clamp(0.08, 1.0);
                self.lbm_free_surface.grid.flags[idx] = if fill >= 0.96 {
                    LbmCell::Fluid
                } else {
                    LbmCell::Interface
                };
                self.lbm_free_surface.mass[idx] = fill;
                let mu_au = gold_effective_viscosity_pa_s(
                    c.temperature_k.max(1340.0),
                    c.melt_fraction.max(0.25),
                );
                self.lbm_free_surface.grid.tau[idx] =
                    (0.56 + 0.15 * (mu_au / 0.005).clamp(0.1, 8.0)).clamp(0.55, 1.85);
            } else {
                self.lbm_free_surface.grid.flags[idx] = LbmCell::Gas;
                self.lbm_free_surface.mass[idx] = 0.0;
            }
        }
        // Repair fluid-to-gas topology closure so Fluid cells never directly touch Gas cells
        if any_lbm_wet {
            for gy in 0..GRID_N {
                for gx in 0..GRID_N {
                    let idx = gy * GRID_N + gx;
                    if self.lbm_free_surface.grid.flags[idx] == LbmCell::Fluid {
                        let mut touches_gas = false;
                        for dy in -1_isize..=1 {
                            for dx in -1_isize..=1 {
                                let nx = gx as isize + dx;
                                let ny = gy as isize + dy;
                                if nx >= 0
                                    && nx < GRID_N as isize
                                    && ny >= 0
                                    && ny < GRID_N as isize
                                    && self.lbm_free_surface.grid.flags
                                        [(ny as usize) * GRID_N + (nx as usize)]
                                        == LbmCell::Gas
                                {
                                    touches_gas = true;
                                }
                            }
                        }
                        if touches_gas {
                            self.lbm_free_surface.grid.flags[idx] = LbmCell::Interface;
                        }
                    }
                }
            }
            self.lbm_free_surface.carry = 0.0;
            self.lbm_free_surface.step();
        }

        // 2. Compute discrete free-surface curvature Laplacian(eta) for Young-Laplace capillary pressure
        let neighbors: [(isize, isize); 4] = [(-1, 0), (1, 0), (0, -1), (0, 1)];
        for gy in 0..GRID_N {
            for gx in 0..GRID_N {
                let idx = gy * GRID_N + gx;
                if !self.grid[idx].active {
                    continue;
                }
                let eta_c = free_surface_eta_mm[idx];
                let mut lap = 0.0;
                for (dx, dy) in neighbors {
                    let nx = gx as isize + dx;
                    let ny = gy as isize + dy;
                    if nx >= 0 && nx < GRID_N as isize && ny >= 0 && ny < GRID_N as isize {
                        let nidx = (ny as usize) * GRID_N + (nx as usize);
                        if self.grid[nidx].active {
                            lap += free_surface_eta_mm[nidx] - eta_c;
                        }
                    }
                }
                eta_laplacian_mm[idx] = lap / (GRID_DX_MM * GRID_DX_MM);
            }
        }

        // 3. Evaluate 2D Finite-Volume fluxes:
        //    (a) Fourier heat conduction -k * grad(T)
        //    (b) Lubrication + Young-Laplace Capillary + Wetting Cohesion + Marangoni molten gold flux
        //    (c) Watch oil Marangoni/capillary spreading flux
        for gy in 0..GRID_N {
            for gx in 0..GRID_N {
                let idx = gy * GRID_N + gx;
                if !self.grid[idx].active {
                    continue;
                }
                let c = &self.grid[idx];
                let card_c = c.material.card();

                // Only visit Right (+1, 0) and Up (0, +1) faces to conserve mass & energy exactly
                for (dx, dy) in [(1_isize, 0_isize), (0_isize, 1_isize)] {
                    let nx = gx as isize + dx;
                    let ny = gy as isize + dy;
                    if nx < 0 || nx >= GRID_N as isize || ny < 0 || ny >= GRID_N as isize {
                        continue;
                    }
                    let nidx = (ny as usize) * GRID_N + (nx as usize);
                    if !self.grid[nidx].active {
                        continue;
                    }
                    let nc = &self.grid[nidx];
                    let card_n = nc.material.card();

                    // (a) Harmonic-mean thermal conductivity across cell face
                    let k_harm = (2.0 * card_c.conductivity_w_m_k * card_n.conductivity_w_m_k)
                        / (card_c.conductivity_w_m_k + card_n.conductivity_w_m_k).max(1.0);
                    // Include high thermal conductivity of any bridging liquid/solid gold layer!
                    let gold_bridge_mm =
                        0.5 * (c.liquid_gold_mm + c.frozen_gold_mm + nc.liquid_gold_mm + nc.frozen_gold_mm);
                    let g_face_w_k =
                        k_harm * (1.8e-3) + GOLD_18K_CARD.conductivity_w_m_k * (gold_bridge_mm * 1.0e-3);
                    let q_cond = g_face_w_k * (nc.temperature_k - c.temperature_k) * 0.22;
                    net_q_w[idx] += q_cond;
                    net_q_w[nidx] -= q_cond;

                    // (b) Molten gold free-surface flow with Young-Laplace surface tension,
                    //     contact-angle cohesion bead threshold, gravity tilt, and Marangoni shear
                    let h_max = c.liquid_gold_mm.max(nc.liquid_gold_mm);
                    if h_max > 0.015 {
                        let t_face = 0.5 * (c.temperature_k.max(1340.0) + nc.temperature_k.max(1340.0));
                        let sigma_au = gold_liquid_surface_tension_n_m(t_face);

                        // Contact-angle disjoining cohesion pressure: prevents molten gold from flattening
                        // into an unphysical thin sheet; instead it forms rounded sessile beads of thickness h_cap!
                        let cos_theta_c = gold_substrate_contact_angle_cos(c.temperature_k);
                        let cos_theta_n = gold_substrate_contact_angle_cos(nc.temperature_k);
                        let h_cap_c = 0.85 * (1.0 - cos_theta_c).sqrt();
                        let h_cap_n = 0.85 * (1.0 - cos_theta_n).sqrt();

                        // Effective generalized pressure head (in mm of gold equivalent):
                        // P = eta - 0.14 * y_mm (bench tilt toward 6h Tourbillon) - l_cap^2 * Laplacian(eta) + Pi_cohesion(H)
                        let l_cap_mm2 = (sigma_au / (GOLD_18K_CARD.density_kg_m3 * 9.81)) * 1.0e6;
                        let coh_c = if c.liquid_gold_mm > 0.01 && c.liquid_gold_mm < h_cap_c {
                            -0.35 * (h_cap_c - c.liquid_gold_mm)
                        } else {
                            0.0
                        };
                        let coh_n = if nc.liquid_gold_mm > 0.01 && nc.liquid_gold_mm < h_cap_n {
                            -0.35 * (h_cap_n - nc.liquid_gold_mm)
                        } else {
                            0.0
                        };

                        let p_c = free_surface_eta_mm[idx] - 0.09 * c.y_mm
                            - 0.18 * l_cap_mm2 * eta_laplacian_mm[idx]
                            + coh_c;
                        let p_n = free_surface_eta_mm[nidx] - 0.09 * nc.y_mm
                            - 0.18 * l_cap_mm2 * eta_laplacian_mm[nidx]
                            + coh_n;

                        // Marangoni thermocapillary shear term + D2Q9 FreeSurface LBM velocity (`fs_lbm::FreeSurface`)
                        let marangoni_head =
                            -0.0012 * (nc.temperature_k - c.temperature_k);
                        let m_lbm_c = self.lbm_free_surface.grid.moments(idx);
                        let m_lbm_n = self.lbm_free_surface.grid.moments(nidx);
                        let u_face: f64 = 0.5
                            * ((m_lbm_c.u[0] + m_lbm_n.u[0]) * (dx as f64)
                                + (m_lbm_c.u[1] + m_lbm_n.u[1]) * (dy as f64));
                        let lbm_head = (u_face * 1.8_f64).clamp(-0.35, 0.35);

                        let delta_p = (p_c - p_n) + marangoni_head + lbm_head;
                        let (donor, acceptor, dp_abs) = if delta_p >= 0.0 {
                            (idx, nidx, delta_p)
                        } else {
                            (nidx, idx, -delta_p)
                        };

                        let h_donor = self.grid[donor].liquid_gold_mm;
                        // Wetting threshold: a thin sessile bead below 0.18 mm won't advance onto a dry cold cell
                        // unless driven by sufficient hydrostatic/capillary overpressure!
                        let dry_acceptor = self.grid[acceptor].liquid_gold_mm < 0.01;
                        let min_advance_h = if dry_acceptor { 0.16 } else { 0.02 };

                        if h_donor > min_advance_h && dp_abs > 0.015 {
                            let phi_donor = self.grid[donor].melt_fraction.max(0.25);
                            let mu_donor = gold_effective_viscosity_pa_s(
                                self.grid[donor].temperature_k.max(1340.0),
                                phi_donor,
                            );
                            // Poiseuille / lubrication mobility H^3 / (3 * mu)
                            let mobility =
                                ((h_donor * h_donor * h_donor) / (12.0 * mu_donor)).clamp(0.0, 4.5);
                            let dh = (mobility * dp_abs * dt_s)
                                .min(0.24 * (h_donor - 0.05).max(0.0));
                            if dh > 1.0e-6 {
                                new_liq_mm[donor] -= dh;
                                let h_acc_prev = new_liq_mm[acceptor].max(0.0);
                                let h_acc_new = h_acc_prev + dh;
                                let enth_donor = self.grid[donor].liquid_gold_enthalpy_j_kg;
                                new_liq_enth[acceptor] = (h_acc_prev * new_liq_enth[acceptor]
                                    + dh * enth_donor)
                                    / h_acc_new;
                                new_liq_mm[acceptor] = h_acc_new;

                                // Enthalpy advection: flowing molten gold heats the acceptor substrate!
                                let dm_kg = GOLD_18K_CARD.density_kg_m3
                                    * GRID_CELL_AREA_M2
                                    * (dh * 1.0e-3);
                                let q_adv_j = dm_kg
                                    * (enth_donor - self.grid[acceptor].enthalpy_j_kg).max(0.0)
                                    * 0.45;
                                net_q_w[acceptor] += q_adv_j / dt_s;
                            }
                        }
                    }

                    // (c) Watch oil spreading on movement substrate
                    if !c.is_gold_case && !nc.is_gold_case {
                        let diff_oil = nc.oil_um - c.oil_um;
                        if diff_oil.abs() > 2.0 {
                            let t_avg = 0.5 * (c.temperature_k + nc.temperature_k);
                            let mu_oil = watch_oil_viscosity_pa_s(t_avg);
                            let oil_mob = (0.045 / mu_oil.max(0.008)).clamp(0.0, 1.8);
                            let d_oil = 0.22 * oil_mob * diff_oil * dt_s;
                            new_oil_um[idx] += d_oil;
                            new_oil_um[nidx] -= d_oil;
                        }
                    }
                }
            }
        }

        // 4. Apply direct external blowtorch, LN2 Nukiyama boiling, oil pour, and phase transitions on each cell
        for idx in 0..GRID_CELLS {
            if !self.grid[idx].active {
                continue;
            }
            let c = &mut self.grid[idx];
            c.liquid_gold_mm = new_liq_mm[idx].max(0.0);
            c.liquid_gold_enthalpy_j_kg = new_liq_enth[idx];
            c.oil_um = new_oil_um[idx].clamp(0.0, 220.0);

            let card = c.material.card();
            let t_i = c.temperature_k;

            // Natural convection + Stefan-Boltzmann radiation from cell top/bottom faces
            let area_2side = 2.0 * GRID_CELL_AREA_M2;
            let q_amb = 14.0 * area_2side * (T_REF_K - t_i)
                + card.emissivity
                    * STEFAN_BOLTZMANN
                    * area_2side
                    * (T_REF_K.powi(4) - t_i.powi(4));
            let mut q_cell_w = net_q_w[idx] + q_amb;

            // Direct blowtorch heating on cell (x_mm, y_mm, z_bed)
            if self.inputs.blowtorch_power > 1.0e-4 {
                let dx = c.x_mm - self.inputs.blowtorch_pos_mm[0];
                let dy = c.y_mm - self.inputs.blowtorch_pos_mm[1];
                let dz = c.bed_elevation_mm - self.inputs.blowtorch_pos_mm[2];
                let d2 = dx * dx + dy * dy + 0.5 * dz * dz;
                let w = (-d2 / (2.0 * 10.5 * 10.5)).exp();
                let flame_temp_k = 2250.0;
                if flame_temp_k > t_i && w > 1.0e-4 {
                    let exposure = if c.is_gold_case { 1.45 } else { 0.85 };
                    let q_t = self.inputs.blowtorch_power
                        * exposure
                        * w
                        * 0.016
                        * (flame_temp_k - t_i);
                    q_cell_w += q_t;
                    torch_grid_w += q_t;
                }
            }

            // Direct LN2 pour and Nukiyama boiling curve (Leidenfrost film vs Nucleate CHF boiling)
            if self.inputs.ln2_pour_rate > 1.0e-4 {
                let dx = c.x_mm - self.inputs.ln2_pos_mm[0];
                let dy = c.y_mm - self.inputs.ln2_pos_mm[1];
                let w = (-(dx * dx + dy * dy) / (2.0 * 16.0 * 16.0)).exp();
                c.ln2_mass_g =
                    (c.ln2_mass_g + self.inputs.ln2_pour_rate * w * 0.24 * dt_s).min(0.45);
            }
            if c.ln2_mass_g > 1.0e-6 {
                if t_i > T_LN2_SAT_K + 0.05 {
                    let dt_sat = t_i - T_LN2_SAT_K;
                    let h_boil = if t_i <= T_LN2_LEIDENFROST_K {
                        4_600.0 + 6_000.0 * (dt_sat / 48.0).clamp(0.1, 1.0)
                    } else {
                        2_050.0 + 2.3 * (t_i - T_LN2_LEIDENFROST_K)
                    };
                    let q_desire_w = h_boil * area_2side * dt_sat * 1.6;
                    let max_lat_j = c.ln2_mass_g * 1.0e-3 * LN2_LATENT_VAP_J_PER_KG;
                    let q_boil_j = (q_desire_w * dt_s).min(max_lat_j);
                    let q_boil_w = q_boil_j / dt_s;
                    q_cell_w -= q_boil_w;
                    ln2_grid_w += q_boil_w;
                    c.ln2_mass_g =
                        (c.ln2_mass_g - (q_boil_j / LN2_LATENT_VAP_J_PER_KG) * 1.0e3).max(0.0);
                } else {
                    c.ln2_mass_g = (c.ln2_mass_g - 0.015 * dt_s).max(0.0);
                }
            }

            // Direct watch oil dropper / flood stream
            if self.inputs.oil_pour_rate > 1.0e-4 && !c.is_gold_case {
                let dx = c.x_mm - self.inputs.oil_pos_mm[0];
                let dy = c.y_mm - self.inputs.oil_pos_mm[1];
                let w = (-(dx * dx + dy * dy) / (2.0 * 11.5 * 11.5)).exp();
                c.oil_um = (c.oil_um + self.inputs.oil_pour_rate * w * 62.0 * dt_s).min(200.0);
            }

            // Integrate cell enthalpy and invert enthalpy -> (T, f_L)
            let gold_mass_kg = GOLD_18K_CARD.density_kg_m3
                * GRID_CELL_AREA_M2
                * ((c.liquid_gold_mm + c.frozen_gold_mm) * 1.0e-3);
            let eff_mass_kg = (c.base_mass_kg + gold_mass_kg).max(1.0e-5);
            c.enthalpy_j_kg += (q_cell_w * dt_s) / eff_mass_kg;

            let (t_new, phi_new) = card.temperature_and_melt_from_enthalpy(c.enthalpy_j_kg);
            c.temperature_k = t_new;
            c.melt_fraction = phi_new;

            // Stefan solidification / remelting of displaced gold on the movement substrate
            // governed by the solidus-liquidus latent heat exchange with the cell!
            if c.liquid_gold_mm > 1.0e-5 && t_new < GOLD_18K_CARD.t_solidus_k {
                let subcool = (GOLD_18K_CARD.t_solidus_k - t_new).max(8.0);
                let freeze_mm = (0.055 * subcool * dt_s).min(c.liquid_gold_mm);
                c.liquid_gold_mm -= freeze_mm;
                c.frozen_gold_mm += freeze_mm;
            } else if c.frozen_gold_mm > 1.0e-5 && t_new > GOLD_18K_CARD.t_liquidus_k {
                let superheat = t_new - GOLD_18K_CARD.t_liquidus_k;
                let remelt_mm = (0.030 * superheat * dt_s).min(c.frozen_gold_mm);
                c.frozen_gold_mm -= remelt_mm;
                c.liquid_gold_mm += remelt_mm;
            }

            // Atmospheric moisture frost desublimation below 271 K & sublimation above 274 K
            if t_new < 271.0 {
                let subcool = (271.0 - t_new).min(195.0);
                c.frost_um = (c.frost_um + 0.22 * subcool * dt_s).min(145.0);
            } else if t_new > 274.0 && c.frost_um > 0.0 {
                c.frost_um = (c.frost_um - 1.6 * (t_new - 274.0) * dt_s).max(0.0);
            }

            // Thermal evaporation / smoking of watch oil above 520 K
            if t_new > 520.0 && c.oil_um > 0.0 {
                c.oil_um = (c.oil_um - 0.32 * (t_new - 520.0) * dt_s).max(0.0);
            }
        }

        // Thermal oxidation of blued steel hairspring (parabolic Deal-Grove / Cabrera-Mott growth above 490 K)
        // using `fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K`
        let t_hs = self.nodes[WatchNodeId::HairspringBreguet as usize].temperature_k;
        if t_hs > 490.0 {
            let e_act_j_mol = 3800.0 * UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;
            let rate_nm_s = 135.0
                * (-e_act_j_mol / (UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * t_hs)).exp()
                * ((t_hs - 490.0) / 80.0);
            self.hairspring_oxide_nm = (self.hairspring_oxide_nm + rate_nm_s * dt_s).min(420.0);
        }

        (torch_grid_w, ln2_grid_w)
    }

    /// Aggregate the 12 component telemetry nodes directly from the 3D Case Shell & 2D Spatial PDE Grid.
    fn aggregate_nodes_from_continuum(&mut self) {
        let mut sum_temp = [0.0_f64; NODE_COUNT];
        let mut sum_enth = [0.0_f64; NODE_COUNT];
        let mut sum_melt = [0.0_f64; NODE_COUNT];
        let mut max_slump = [0.0_f64; NODE_COUNT];
        let mut sum_liq_gold_g = [0.0_f64; NODE_COUNT];
        let mut sum_frz_gold_g = [0.0_f64; NODE_COUNT];
        let mut sum_ln2_g = [0.0_f64; NODE_COUNT];
        let mut sum_frost = [0.0_f64; NODE_COUNT];
        let mut sum_oil = [0.0_f64; NODE_COUNT];
        let mut counts = [0_usize; NODE_COUNT];

        // 1. Aggregate case sectors (nodes 0..=3) from the 3D Case Shell nodes
        for sn in &self.case_shell {
            let x = sn.ref_pos_mm[0];
            let y = sn.ref_pos_mm[1];
            let sector = if y >= x.abs() {
                0 // 12h Top
            } else if y <= -x.abs() {
                3 // 6h Bottom
            } else if x < 0.0 {
                1 // 9h Left
            } else {
                2 // 3h Right
            };
            sum_temp[sector] += sn.temperature_k;
            sum_enth[sector] += sn.enthalpy_j_kg;
            sum_melt[sector] += sn.melt_fraction;
            sum_ln2_g[sector] += sn.ln2_mass_g;

            let r0 = (sn.ref_pos_mm[0].powi(2) + sn.ref_pos_mm[1].powi(2)).sqrt();
            let r1 = (sn.pos_mm[0].powi(2) + sn.pos_mm[1].powi(2)).sqrt();
            let slump = (((r0 - r1).max(0.0) / 5.5)
                + ((sn.ref_thickness_mm - sn.thickness_mm) / sn.ref_thickness_mm))
                .clamp(0.0, 1.0);
            if slump > max_slump[sector] {
                max_slump[sector] = slump;
            }
            counts[sector] += 1;
        }

        // 2. Aggregate interior movement nodes (5..=11) and gold accumulation from the 24x24 PDE grid
        let cell_gold_g_per_mm =
            GOLD_18K_CARD.density_kg_m3 * GRID_CELL_AREA_M2 * 1.0e-3 * 1.0e3;
        for c in &self.grid {
            if !c.active {
                continue;
            }
            let owner = c.owner_node;
            sum_liq_gold_g[owner] += c.liquid_gold_mm * cell_gold_g_per_mm;
            sum_frz_gold_g[owner] += c.frozen_gold_mm * cell_gold_g_per_mm;
            if owner >= 5 {
                sum_temp[owner] += c.temperature_k;
                sum_enth[owner] += c.enthalpy_j_kg;
                sum_melt[owner] += c.melt_fraction;
                sum_ln2_g[owner] += c.ln2_mass_g;
                sum_frost[owner] += c.frost_um;
                sum_oil[owner] += c.oil_um;
                counts[owner] += 1;
            } else {
                sum_frost_case(&mut sum_frost[owner], c.frost_um);
            }
        }

        fn sum_frost_case(acc: &mut f64, val: f64) {
            if val > *acc {
                *acc = val;
            }
        }

        for idx in 0..NODE_COUNT {
            if idx == WatchNodeId::CrystalSapphire as usize {
                continue;
            }
            let n = counts[idx].max(1) as f64;
            let node = &mut self.nodes[idx];
            node.temperature_k = sum_temp[idx] / n;
            node.enthalpy_j_kg = sum_enth[idx] / n;
            node.melt_fraction = sum_melt[idx] / n;
            if idx <= 3 {
                node.plastic_slump = max_slump[idx];
                node.frost_thickness_um = sum_frost[idx];
            } else {
                node.frost_thickness_um = sum_frost[idx] / n;
                node.oil_film_um = sum_oil[idx] / n;
            }
            node.molten_gold_g = sum_liq_gold_g[idx];
            node.frozen_gold_g = sum_frz_gold_g[idx];
            node.ln2_film_g = sum_ln2_g[idx];
            if idx == WatchNodeId::HairspringBreguet as usize {
                node.oxide_thickness_nm = self.hairspring_oxide_nm;
            }
            let card = node.material.card();
            node.incandescence_rgb =
                planck_incandescence_linear_rgb(node.temperature_k, card.emissivity);
        }

        // Sapphire crystal temperature tracks radiative/conductive mean of case rim and mainplate
        let t_case_avg = 0.25
            * (self.nodes[0].temperature_k
                + self.nodes[1].temperature_k
                + self.nodes[2].temperature_k
                + self.nodes[3].temperature_k);
        let t_plate = self.nodes[WatchNodeId::MainplateBrass as usize].temperature_k;
        let frost_crystal = 0.5
            * (self.nodes[WatchNodeId::MainplateBrass as usize].frost_thickness_um
                + self.nodes[3].frost_thickness_um);
        let crystal = &mut self.nodes[WatchNodeId::CrystalSapphire as usize];
        crystal.temperature_k = 0.65 * t_case_avg + 0.35 * t_plate;
        crystal.enthalpy_j_kg = crystal
            .material
            .card()
            .enthalpy_from_temperature_j_kg(crystal.temperature_k);
        crystal.frost_thickness_um = frost_crystal;
    }

    fn recompute_modifiers(&mut self) {
        let t_plate = self.nodes[WatchNodeId::MainplateBrass as usize].temperature_k;
        let t_barrel = self.nodes[WatchNodeId::BarrelMainspring as usize].temperature_k;
        let t_train = self.nodes[WatchNodeId::GoingTrainWheels as usize].temperature_k;
        let t_cage = self.nodes[WatchNodeId::TourbillonCage as usize].temperature_k;
        let t_esc = self.nodes[WatchNodeId::EscapementForkWheel as usize].temperature_k;
        let t_bal = self.nodes[WatchNodeId::BalanceWheelRim as usize].temperature_k;
        let t_sp = self.nodes[WatchNodeId::HairspringBreguet as usize].temperature_k;
        let barrel_youngs_modulus_ratio = HIGH_CARBON_STEEL_CARD
            .youngs_modulus_ratio(t_barrel)
            .clamp(0.35, 1.2);

        // 1. Balance wheel moment of inertia I(T) = I_0 * (1 + eps_BeCu(T))^2 + rotational inertia of any gold on the balance rim
        //    plus `fs_mbd_material::VolumeMoments` for Guillaume mass-poising inertia & `fs_vkplate::air_loss::oscillating_cylinder`
        //    for barometric air added-mass and Stokes-Stuart oscillating-cylinder aerodynamic drag!
        let bal_card = self.nodes[WatchNodeId::BalanceWheelRim as usize]
            .material
            .card();
        let eps_bal = bal_card.thermal_strain(t_bal);
        let poise_screw_vol_m3 = VolumeMoments::sphere(0.28e-3)
            .map(|vm| vm.volume_m3())
            .unwrap_or(9.2e-11);
        let poise_screw_mass_kg = 16.0 * poise_screw_vol_m3 * GOLD_18K_CARD.density_kg_m3;
        let (air_added_mass_kg, balance_aero_drag_s_inv) =
            GasState::try_new_moist_air(t_bal.clamp(220.0, 365.0), 101_325.0, 0.35)
                .ok()
                .and_then(|air_gas| {
                    oscillating_cylinder(0.275e-3, 2.0 * std::f64::consts::PI * 2.5, &air_gas)
                        .ok()
                        .map(|cd| {
                            let ring_circumference_m = 2.0 * std::f64::consts::PI * 0.0052;
                            let added_m = cd.added_mass_kg_per_m * ring_circumference_m;
                            let drag_s_inv = (cd.resistance_n_s_per_m2 * ring_circumference_m
                                / (0.00036 + poise_screw_mass_kg))
                                .clamp(0.001, 2.5);
                            (added_m, drag_s_inv)
                        })
                })
                .unwrap_or((1.8e-9, 0.018));
        let gold_on_bal_g = self.nodes[WatchNodeId::BalanceWheelRim as usize].molten_gold_g
            + self.nodes[WatchNodeId::BalanceWheelRim as usize].frozen_gold_g
            + 0.35
                * (self.nodes[WatchNodeId::TourbillonCage as usize].molten_gold_g
                    + self.nodes[WatchNodeId::TourbillonCage as usize].frozen_gold_g);
        let mass_loading_ratio =
            1.0 + (gold_on_bal_g * 1.0e-3 + air_added_mass_kg) / (0.00036 + poise_screw_mass_kg);
        let balance_inertia_ratio =
            ((1.0 + eps_bal) * (1.0 + eps_bal) * mass_loading_ratio).clamp(0.85, 25.0);

        // 2. Hairspring bending stiffness k(T) = E(T)*b(T)*h(T)^3 / (12*L(T))
        let sp_card = if self.inputs.compensated_hairspring {
            &INVAR_36NI_CARD
        } else {
            &HIGH_CARBON_STEEL_CARD
        };
        let e_ratio = sp_card.youngs_modulus_ratio(t_sp);
        let eps_sp = sp_card.thermal_strain(t_sp);
        let k_thermal_ratio = (e_ratio * (1.0 + 3.0 * eps_sp)).clamp(0.15, 1.35);

        // Capillary meniscus coil-sticking when watch oil floods the hairspring cell
        let oil_sp_um = self.nodes[WatchNodeId::HairspringBreguet as usize].oil_film_um;
        let gamma_oil = watch_oil_surface_tension_n_m(t_sp);
        let bridge_frac = if oil_sp_um > 2.5 {
            let excess = oil_sp_um - 2.5;
            ((excess / 38.0) * (gamma_oil / 0.031)).clamp(0.0, 0.38)
        } else {
            0.0
        };
        let k_oil_ratio = 1.0 / (1.0 - 0.65 * bridge_frac);

        // 3. 3D Dipolar Magnetic Field B(r) from external NdFeB N42 magnet + ferromagnetic remanence
        let mut b_esc_mt = 0.0;
        let mut b_signed_dir = 1.0;
        if self.inputs.magnet_enabled && self.inputs.magnet_strength > 1.0e-4 {
            let esc_pos = self.nodes[WatchNodeId::EscapementForkWheel as usize].pos_mm;
            let dx = esc_pos[0] - self.inputs.magnet_pos_mm[0];
            let dy = esc_pos[1] - self.inputs.magnet_pos_mm[1];
            let dz = esc_pos[2] - self.inputs.magnet_pos_mm[2];
            let r_mm = (dx * dx + dy * dy + dz * dz).sqrt().max(6.5);
            let r_eff_mm = 11.5;
            b_esc_mt = (self.inputs.magnet_strength * 680.0 * (r_eff_mm / r_mm).powi(3))
                .clamp(0.0, 950.0);
            b_signed_dir = if dx >= 0.0 { 1.0 } else { -1.0 };

            if b_esc_mt > 3.5 {
                let target_rem = (0.32 * b_esc_mt).min(28.0);
                if target_rem > self.remanent_magnetization_mt {
                    self.remanent_magnetization_mt = target_rem;
                }
            }
        }
        let b_eff_mt = b_esc_mt + self.remanent_magnetization_mt;
        let mag_sus_factor = if self.inputs.compensated_hairspring {
            0.55
        } else {
            1.0
        };
        let k_mag_ratio = 1.0 + mag_sus_factor * 0.0022 * b_eff_mt.min(120.0);
        let hairspring_eq_bias_rad =
            b_signed_dir * mag_sus_factor * (0.0045 * b_eff_mt).clamp(-0.32, 0.32);

        let pallet_magnetic_bias_torque_nm =
            b_signed_dir * 2.4e-8 * (b_esc_mt * b_esc_mt / (1.0 + 0.04 * b_esc_mt))
                + b_signed_dir * 1.8e-8 * self.remanent_magnetization_mt;
        let magnetically_pinned = b_esc_mt > 22.0;

        let hairspring_stiffness_ratio =
            (k_thermal_ratio * k_oil_ratio * k_mag_ratio).clamp(0.15, 2.50);
        let natural_freq_ratio = (hairspring_stiffness_ratio / balance_inertia_ratio).sqrt();
        let theoretical_rate_s_per_day = (natural_freq_ratio - 1.0) * 86_400.0;

        // 4. Differential thermal expansion of pivots inside ruby jewel holes
        let eps_plate = self.nodes[WatchNodeId::MainplateBrass as usize]
            .material
            .card()
            .thermal_strain(t_plate);
        let eps_ruby = self.nodes[WatchNodeId::CrystalSapphire as usize]
            .material
            .card()
            .thermal_strain(t_cage);
        let eps_steel_esc = HIGH_CARBON_STEEL_CARD.thermal_strain(t_esc);
        let eps_steel_train = HIGH_CARBON_STEEL_CARD.thermal_strain(t_train);

        let cryo_pinch_um = 420.0 * (eps_plate - eps_steel_train).min(0.0);
        let hot_arbor_pinch_um = 380.0 * (eps_ruby - eps_steel_esc).min(0.0);
        let min_pivot_clearance_um = (4.5 + cryo_pinch_um + hot_arbor_pinch_um).clamp(-2.0, 8.0);
        let thermal_or_cryo_seized = min_pivot_clearance_um <= 0.15;

        // 5. Stribeck lubrication & Concentric Couette journal bearing drag (`fs_tribo::journal::GeneralizedNewtonianJournal`)
        let mu_oil_esc = watch_oil_viscosity_pa_s(0.5 * (t_esc + t_bal));
        let mu_ref = watch_oil_viscosity_pa_s(T_REF_K);
        let pivot_r_m = 0.06e-3;
        let clearance_m = (min_pivot_clearance_um.max(0.25) * 1.0e-6).clamp(0.25e-6, 12.0e-6);
        let ref_clearance_m = 4.5e-6;
        let journal_solve = JournalShearSolve {
            maximum_iterations: 24,
            absolute_angular_rate_tolerance_per_s: 1.0e-6,
            relative_angular_rate_tolerance: 1.0e-6,
        };
        let journal_now = GeneralizedNewtonianJournal::new(
            JournalGeometry {
                inner_radius_m: pivot_r_m,
                outer_radius_m: pivot_r_m + clearance_m,
                length_m: 0.18e-3,
            },
            vec![
                JournalViscosityKnot {
                    shear_rate_per_s: 0.0,
                    dynamic_viscosity_pa_s: mu_oil_esc,
                },
                JournalViscosityKnot {
                    shear_rate_per_s: 2.0e6,
                    dynamic_viscosity_pa_s: (mu_oil_esc * 0.92).max(1.0e-4),
                },
            ],
            journal_solve,
        )
        .ok();
        let journal_ref = GeneralizedNewtonianJournal::new(
            JournalGeometry {
                inner_radius_m: pivot_r_m,
                outer_radius_m: pivot_r_m + ref_clearance_m,
                length_m: 0.18e-3,
            },
            vec![
                JournalViscosityKnot {
                    shear_rate_per_s: 0.0,
                    dynamic_viscosity_pa_s: mu_ref,
                },
                JournalViscosityKnot {
                    shear_rate_per_s: 2.0e6,
                    dynamic_viscosity_pa_s: mu_ref * 0.92,
                },
            ],
            journal_solve,
        )
        .ok();
        let couette_torque_ratio = match (journal_now, journal_ref) {
            (Some(jn), Some(jr)) => {
                let t_now = jn
                    .response(45.0)
                    .map(|r| r.resisting_torque_nm.abs())
                    .unwrap_or(1.0e-9);
                let t_ref = jr
                    .response(45.0)
                    .map(|r| r.resisting_torque_nm.abs())
                    .unwrap_or(1.0e-9)
                    .max(1.0e-12);
                (t_now / t_ref).clamp(0.08, 1200.0)
            }
            _ => (mu_oil_esc / mu_ref).clamp(0.08, 1200.0),
        };
        let visc_ratio = (mu_oil_esc / mu_ref).clamp(0.08, 1200.0);
        let oil_esc_um = self.nodes[WatchNodeId::EscapementForkWheel as usize].oil_film_um;
        let starvation_penalty = if oil_esc_um < 0.35 {
            1.0 + 5.5 * ((0.35 - oil_esc_um) / 0.35)
        } else {
            1.0
        };
        let clearance_penalty = if min_pivot_clearance_um < 1.2 {
            1.0 + 18.0 * ((1.2 - min_pivot_clearance_um.max(0.0)) / 1.2).powi(2)
        } else {
            1.0
        };
        let pivot_friction_multiplier = (starvation_penalty
            * clearance_penalty
            * (0.65 + 0.35 * couette_torque_ratio.min(35.0)))
        .clamp(0.65, 250.0);

        let oil_bal_um = self.nodes[WatchNodeId::BalanceWheelRim as usize].oil_film_um;
        let flood_um = (oil_bal_um.max(oil_sp_um) - 3.5).max(0.0);
        let oil_flood_damping = 0.16 * flood_um * visc_ratio.sqrt().min(12.0);
        let cryo_visc_damping = if visc_ratio > 8.0 {
            0.035 * (visc_ratio - 8.0).min(450.0)
        } else {
            0.0
        };
        let eddy_damping = 0.000_28 * b_eff_mt * b_eff_mt;

        // 6. Spatially-integrated hydrodynamic viscous torque & solid yield interference torque
        //    from mobile liquid gold (H_Au) and solidified gold (S_Au) occupying the moving wheel footprints!
        let mut tau_visc_nm = 0.0;
        let mut tau_solid_nm = 0.0;
        for c in &self.grid {
            if !c.active || c.is_gold_case {
                continue;
            }
            if c.liquid_gold_mm <= 1.0e-4 && c.frozen_gold_mm <= 1.0e-4 {
                continue;
            }
            // Check geometric overlap with the moving mechanisms (Barrel=6, Train=7, Tourbillon=8..=11)
            let (r_lever_m, gear_gain) = match c.owner_node {
                6 => (0.0055, 1.0),    // Mainspring Barrel teeth
                7 => (0.0038, 8.0),    // Going train wheels (reflected by gear ratio)
                8..=11 => (0.0042, 45.0), // Tourbillon carriage, escape wheel & balance rim (high mechanical advantage!)
                _ => (0.0020, 0.25),
            };
            if c.liquid_gold_mm > 1.0e-4 {
                let mu_au = gold_effective_viscosity_pa_s(
                    c.temperature_k.max(1340.0),
                    c.melt_fraction.max(0.2),
                );
                let h_gap_m = 0.4e-3;
                let shear_pa = mu_au * (15.0 / h_gap_m) * (c.liquid_gold_mm / 0.35).min(3.0);
                tau_visc_nm += shear_pa * GRID_CELL_AREA_M2 * r_lever_m * gear_gain;
            }
            if c.frozen_gold_mm > 1.0e-4 {
                let sigma_y = gold_viscoplastic_yield_stress_pa(c.temperature_k, c.melt_fraction, 0.1);
                let contact_area_m2 =
                    (c.frozen_gold_mm * 1.0e-3) * (GRID_DX_M * 0.35);
                tau_solid_nm += 0.577 * sigma_y * contact_area_m2 * r_lever_m * gear_gain;
            }
        }

        // Compare resistive torque against the 8.2 mN*m nominal barrel drive torque:
        let drive_ref_nm = 0.0082;
        let gold_fouling_jam_factor =
            ((tau_visc_nm * 12.0 + tau_solid_nm) / drive_ref_nm).clamp(0.0, 1.0);
        let frozen_gold_locked = tau_solid_nm > drive_ref_nm * 0.45;
        let gold_damping = 42.0 * gold_fouling_jam_factor;

        let extra_balance_damping_s_inv = oil_flood_damping
            + cryo_visc_damping
            + eddy_damping
            + gold_damping
            + (balance_aero_drag_s_inv - 0.018).max(0.0);

        // 7. Forward-mode Automatic Differentiation (`fs_ad::Dual64`) for exact thermoelastic rate sensitivity d(rate)/dT
        //    and Certified Interval Arithmetic (`fs_ivl::Interval`) for worst-case manufacturing tolerance enclosure!
        let dt_ad = Dual64::<1>::variable((t_sp - T_REF_K).clamp(-200.0, 600.0), 0);
        let alpha_sp_ad = Dual64::<1>::constant(sp_card.cte_per_k);
        let alpha_bal_ad = Dual64::<1>::constant(bal_card.cte_per_k);
        let d_e_dt_ad = Dual64::<1>::constant(sp_card.thermoelastic_beta_per_k);
        let e_ratio_ad = Dual64::<1>::constant(1.0) + d_e_dt_ad * dt_ad;
        let k_sp_ad = e_ratio_ad
            * (Dual64::<1>::constant(1.0) + Dual64::<1>::constant(3.0) * alpha_sp_ad * dt_ad)
            * Dual64::<1>::constant(k_oil_ratio * k_mag_ratio);
        let i_bal_ad = (Dual64::<1>::constant(1.0) + alpha_bal_ad * dt_ad)
            * (Dual64::<1>::constant(1.0) + alpha_bal_ad * dt_ad)
            * Dual64::<1>::constant(mass_loading_ratio);
        let rate_ad = ((k_sp_ad / i_bal_ad).sqrt() - Dual64::<1>::constant(1.0))
            * Dual64::<1>::constant(86_400.0);
        let rate_sensitivity_s_day_per_k = rate_ad.eps[0];

        // Interval arithmetic (`fs_ivl::Interval`): propagate +-1.5% hairspring thickness tolerance and +-0.4 K thermal gradient
        let k_ivl = Interval::new(
            (hairspring_stiffness_ratio * 0.99975).max(0.10),
            hairspring_stiffness_ratio * 1.00025,
        );
        let i_ivl = Interval::new(
            (balance_inertia_ratio * 0.99985).max(0.50),
            balance_inertia_ratio * 1.00015,
        );
        let ratio_lo = (k_ivl.lo() / i_ivl.hi()).max(0.01).sqrt();
        let ratio_hi = (k_ivl.hi() / i_ivl.lo()).max(0.01).sqrt();
        let rate_ivl = Interval::new((ratio_lo - 1.0) * 86_400.0, (ratio_hi - 1.0) * 86_400.0);

        self.modifiers = CoupledWatchModifiers {
            balance_inertia_ratio,
            hairspring_stiffness_ratio,
            barrel_youngs_modulus_ratio,
            natural_freq_ratio,
            theoretical_rate_s_per_day,
            rate_sensitivity_s_day_per_k,
            rate_interval_lo_s_day: rate_ivl.lo(),
            rate_interval_hi_s_day: rate_ivl.hi(),
            balance_aero_drag_s_inv,
            pivot_friction_multiplier,
            extra_balance_damping_s_inv,
            hairspring_oil_bridge_fraction: bridge_frac,
            min_pivot_clearance_um,
            b_field_escapement_mt: b_esc_mt,
            remanent_magnetization_mt: self.remanent_magnetization_mt,
            pallet_magnetic_bias_torque_nm,
            hairspring_eq_bias_rad,
            gold_fouling_jam_factor,
            liquid_gold_viscous_torque_nm: tau_visc_nm,
            solid_gold_yield_torque_nm: tau_solid_nm,
            frozen_gold_locked,
            thermal_or_cryo_seized,
            magnetically_pinned,
            ln2_boiling_power_w: self.modifiers.ln2_boiling_power_w,
            blowtorch_heating_power_w: self.modifiers.blowtorch_heating_power_w,
            blued_steel_rgb: blued_steel_interference_rgb(self.hairspring_oxide_nm, 1.0),
        };
    }
}

impl Default for MultiphysicsState {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn blowtorch_heats_watch_shifts_rate_then_melts_gold_into_inner_workings() {
        let mut mp = MultiphysicsState::new();
        mp.inputs.blowtorch_power = 1.0;
        mp.inputs.blowtorch_pos_mm = [0.0, -25.0, 0.0];

        // 1. Moderate heating (1.0 s): watch warms up, losing time due to E(T) drop + rim expansion
        for _ in 0..20 {
            mp.step(0.05);
        }
        assert!(mp.nodes[WatchNodeId::CaseBottom18kGold as usize].temperature_k > 400.0);
        assert!(mp.modifiers.theoretical_rate_s_per_day < -10.0);

        // 2. Sustained blowtorch heating (7.0 more seconds): 18k gold case exceeds solidus (1163 K),
        //    deforms viscoplastically under Johnson-Cook + Young-Laplace surface tension, and sheds
        //    molten gold into the 6h Tourbillon well!
        for _ in 0..140 {
            mp.step(0.05);
        }
        let case_bottom = mp.nodes[WatchNodeId::CaseBottom18kGold as usize];
        assert!(
            case_bottom.temperature_k > GOLD_18K_CARD.t_solidus_k,
            "expected case bottom T > solidus, got {}",
            case_bottom.temperature_k
        );
        assert!(
            case_bottom.melt_fraction > 0.1,
            "expected melt_fraction > 0.1, got {}",
            case_bottom.melt_fraction
        );
        assert!(
            case_bottom.plastic_slump > 0.05,
            "expected emergent 3D shell deformation > 0.05, got {}",
            case_bottom.plastic_slump
        );
        assert!(
            mp.modifiers.gold_fouling_jam_factor > 0.25,
            "expected emergent viscous/solid gold fouling torque factor > 0.25, got {}",
            mp.modifiers.gold_fouling_jam_factor
        );

        // 3. Turn off blowtorch and pour Liquid Nitrogen (LN2 at 77.36 K):
        //    molten gold quenches below solidus and permanently locks the movement!
        mp.inputs.blowtorch_power = 0.0;
        mp.inputs.ln2_pour_rate = 1.0;
        mp.inputs.ln2_pos_mm = [0.0, -18.0, 2.0];
        for _ in 0..120 {
            mp.step(0.05);
        }
        assert!(
            mp.modifiers.frozen_gold_locked,
            "expected solidified gold yield torque to lock movement after LN2 quench"
        );
        assert!(mp.nodes[WatchNodeId::CaseBottom18kGold as usize].plastic_slump > 0.05);
    }

    #[test]
    fn magnet_pins_pallet_fork_and_leaves_remanence_until_degaussed() {
        let mut mp = MultiphysicsState::new();
        mp.inputs.magnet_enabled = true;
        mp.inputs.magnet_strength = 1.2;
        mp.inputs.magnet_pos_mm = [-2.0, -12.0, 4.0];
        mp.step(0.02);

        assert!(mp.modifiers.b_field_escapement_mt > 50.0);
        assert!(mp.modifiers.magnetically_pinned);
        assert!(mp.modifiers.remanent_magnetization_mt > 5.0);

        mp.inputs.magnet_enabled = false;
        mp.step(0.02);
        assert!(!mp.modifiers.magnetically_pinned);
        assert!(mp.modifiers.remanent_magnetization_mt > 5.0);

        mp.trigger_degauss();
        for _ in 0..40 {
            mp.step(0.05);
        }
        assert!(mp.modifiers.remanent_magnetization_mt < 0.2);
    }
}
