//! Going train, mainspring, balance and timing of the FS-1903 tourbillon.
//!
//! Nothing in this module sets a rate, an amplitude or a stop condition. They
//! come out of three pieces of physics:
//!
//! - the contact-dynamics lever escapement ([`crate::escapement`]), which turns
//!   the escape-wheel torque into impulses on the balance and produces the
//!   lock drops a timegrapher hears;
//! - the loads on the balance, lever and escape wheel, assembled here from the
//!   geometry of the parts and the local fields that the multiphysics model
//!   samples at each part (temperature, air, oil, frost, liquid and frozen
//!   gold, magnetic field, gravity);
//! - the mainspring and going train, which set the escape-wheel torque from
//!   the state of wind, the steel's modulus at the barrel temperature and the
//!   tooth-friction efficiency of each mesh.
//!
//! Rate, beat error and amplitude are then *measured* from the lock-drop times
//! and the balance turning points, as a timegrapher would, and a stopped watch
//! is diagnosed from whichever resistance dominates.

use core::f64::consts::{PI, TAU};
use std::collections::VecDeque;

use fs_eproc::BettingEProcess;
use fs_tribo::journal::JournalGeometry;

use crate::escapement::{ContactImpact, ContactKind, Escapement, EscapementDesign, EscapementLoads};
use crate::materials::{G0_M_S2, HairspringAlloy, P_ATM_PA, T_REF_K, air_state, library};

pub const NOMINAL_FREQ_HZ: f64 = 2.5;
pub const NOMINAL_BEAT_PERIOD_S: f64 = 1.0 / (2.0 * NOMINAL_FREQ_HZ);
pub const MAINSPRING_COILS: usize = 18;
const DEG: f64 = PI / 180.0;
const SECONDS_PER_DAY: f64 = 86_400.0;
const R_GAS: f64 = 8.314_462_618;

/// Tooth counts of the going train. The fourth wheel is fixed to the plate;
/// the escape pinion rolls round it, so the escape wheel turns ten times per
/// turn of the one-minute cage.
pub mod train {
    pub const BARREL: f64 = 78.0;
    pub const CENTER_PINION: f64 = 10.0;
    pub const CENTER: f64 = 80.0;
    pub const THIRD_PINION: f64 = 10.0;
    pub const THIRD: f64 = 75.0;
    pub const CAGE_PINION: f64 = 10.0;
    pub const FIXED_FOURTH: f64 = 70.0;
    pub const ESCAPE_PINION: f64 = 7.0;
    pub const BARREL_TO_CENTER: f64 = BARREL / CENTER_PINION;
    pub const CENTER_TO_THIRD: f64 = CENTER / THIRD_PINION;
    pub const THIRD_TO_CAGE: f64 = THIRD / CAGE_PINION;
    pub const CAGE_TO_ESCAPE: f64 = FIXED_FOURTH / ESCAPE_PINION;
    /// Escape-wheel turns (in the cage) per barrel turn.
    pub const RATIO: f64 = BARREL_TO_CENTER * CENTER_TO_THIRD * THIRD_TO_CAGE * CAGE_TO_ESCAPE;
    /// Unlubricated brass-wheel / steel-pinion sliding friction.
    pub const MESH_MU: f64 = 0.2;
}

/// Efficiency of one external mesh with sliding friction `mu` over the
/// approach and recess arcs: `1 − μ π (1/z₁ + 1/z₂)` (Merritt; Buckingham).
#[must_use]
pub fn mesh_efficiency(z1: f64, z2: f64, mu: f64) -> f64 {
    (1.0 - mu * PI * (1.0 / z1 + 1.0 / z2)).clamp(0.0, 1.0)
}

/// Product of the four mesh efficiencies from barrel to escape pinion.
#[must_use]
pub fn train_efficiency(mu: f64) -> f64 {
    use train::*;
    mesh_efficiency(BARREL, CENTER_PINION, mu)
        * mesh_efficiency(CENTER, THIRD_PINION, mu)
        * mesh_efficiency(THIRD, CAGE_PINION, mu)
        * mesh_efficiency(FIXED_FOURTH, ESCAPE_PINION, mu)
}

/// Physical fields the multiphysics model samples at one part of the movement.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PartField {
    pub temperature_k: f64,
    /// Depth of mobile (liquid or mushy) gold standing on the plate under the part [m].
    pub liquid_gold_m: f64,
    /// Thickness of frozen gold on the plate under the part [m].
    pub solid_gold_m: f64,
    pub gold_liquid_fraction: f64,
    pub gold_temperature_k: f64,
    /// Depth of poured oil standing on the plate [m].
    pub oil_m: f64,
    /// Frost or ice thickness grown on the part's surfaces [m].
    pub frost_m: f64,
}

impl Default for PartField {
    fn default() -> Self {
        Self {
            temperature_k: T_REF_K,
            liquid_gold_m: 0.0,
            solid_gold_m: 0.0,
            gold_liquid_fraction: 0.0,
            gold_temperature_k: T_REF_K,
            oil_m: 0.0,
            frost_m: 0.0,
        }
    }
}

/// Everything the movement feels from the outside world.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MovementEnvironment {
    pub balance: PartField,
    pub escapement: PartField,
    pub train: PartField,
    pub barrel: PartField,
    pub air_pressure_pa: f64,
    pub air_temperature_k: f64,
    /// Speed of the compressed-air jet where it crosses the balance [m/s].
    pub air_jet_speed_m_s: f64,
    /// Gravity in the watch frame [m/s²]: x towards 3 o'clock, y towards 12,
    /// z out of the dial.
    pub gravity_m_s2: [f64; 3],
    /// Flux density at the balance [T] and the magnitude of its gradient [T/m].
    pub b_balance_t: f64,
    pub grad_b_balance_t_m: f64,
    /// Signed hairspring magnetisation M/Mₛ (remanent plus induced).
    pub hairspring_magnetisation: f64,
    /// Magnetic torque on the steel pallet lever about its arbor [N m].
    pub lever_magnetic_torque_nm: f64,
}

impl Default for MovementEnvironment {
    fn default() -> Self {
        Self {
            balance: PartField::default(),
            escapement: PartField::default(),
            train: PartField::default(),
            barrel: PartField::default(),
            air_pressure_pa: P_ATM_PA,
            air_temperature_k: T_REF_K,
            air_jet_speed_m_s: 0.0,
            gravity_m_s2: [0.0, -G0_M_S2, 0.0],
            b_balance_t: 0.0,
            grad_b_balance_t_m: 0.0,
            hairspring_magnetisation: 0.0,
            lever_magnetic_torque_nm: 0.0,
        }
    }
}

/// Glucydur balance with gold timing screws, its jewelled pivots and the
/// hairspring above it. Heights are measured from the main plate.
#[derive(Clone, Copy, Debug)]
pub struct BalanceDesign {
    pub rim_outer_m: f64,
    pub rim_inner_m: f64,
    pub rim_height_m: f64,
    pub arm_count: f64,
    pub arm_width_m: f64,
    pub arm_height_m: f64,
    pub hub_radius_m: f64,
    pub screw_count: f64,
    pub screw_radius_m: f64,
    pub screw_diameter_m: f64,
    pub screw_length_m: f64,
    pub staff_mass_kg: f64,
    pub pivot_radius_m: f64,
    pub pivot_length_m: f64,
    pub pivot_clearance_m: f64,
    pub jewel_mu_oiled: f64,
    pub jewel_mu_dry: f64,
    pub rim_bottom_z_m: f64,
    pub gap_below_m: f64,
    pub gap_above_m: f64,
    pub gap_side_m: f64,
    /// Residual out-of-poise (mass × eccentricity) and its direction.
    pub poise_kg_m: f64,
    pub poise_angle_rad: f64,
    pub hairspring_width_m: f64,
    pub hairspring_thickness_m: f64,
    pub hairspring_coils: f64,
    pub hairspring_inner_m: f64,
    pub hairspring_outer_m: f64,
    pub hairspring_bottom_z_m: f64,
    /// Glucydur electrical conductivity [S/m] (C17200, ~22 % IACS).
    pub rim_conductivity_s_m: f64,
}

impl BalanceDesign {
    #[must_use]
    pub fn pocket_tourbillon() -> Self {
        Self {
            rim_outer_m: 5.6e-3,
            rim_inner_m: 4.9e-3,
            rim_height_m: 0.45e-3,
            arm_count: 2.0,
            arm_width_m: 0.5e-3,
            arm_height_m: 0.35e-3,
            hub_radius_m: 0.5e-3,
            screw_count: 8.0,
            screw_radius_m: 5.9e-3,
            screw_diameter_m: 0.4e-3,
            screw_length_m: 0.6e-3,
            staff_mass_kg: 5.0e-6,
            pivot_radius_m: 0.04e-3,
            pivot_length_m: 0.06e-3,
            pivot_clearance_m: 3.0e-6,
            jewel_mu_oiled: 0.15,
            jewel_mu_dry: 0.25,
            rim_bottom_z_m: 1.6e-3,
            gap_below_m: 0.3e-3,
            gap_above_m: 0.45e-3,
            gap_side_m: 0.35e-3,
            poise_kg_m: 1.0e-10,
            poise_angle_rad: 0.7,
            hairspring_width_m: 0.15e-3,
            hairspring_thickness_m: 0.035e-3,
            hairspring_coils: 13.0,
            hairspring_inner_m: 1.0e-3,
            hairspring_outer_m: 4.0e-3,
            hairspring_bottom_z_m: 2.3e-3,
            rim_conductivity_s_m: 1.28e7,
        }
    }

    fn rim_volume(&self) -> f64 {
        PI * (self.rim_outer_m.powi(2) - self.rim_inner_m.powi(2)) * self.rim_height_m
    }

    /// Mass [kg] and moment of inertia [kg m²] at the reference temperature,
    /// from the Glucydur rim and arms and the 18k gold screws.
    #[must_use]
    pub fn mass_and_inertia(&self) -> (f64, f64) {
        let lib = library();
        let rho_rim = lib.becu.rho;
        let rho_screw = lib.gold.rho_solid_ref;
        let m_rim = rho_rim * self.rim_volume();
        let i_rim = 0.5 * m_rim * (self.rim_outer_m.powi(2) + self.rim_inner_m.powi(2));
        let lin = rho_rim * self.arm_width_m * self.arm_height_m;
        let (r1, r2) = (self.hub_radius_m, self.rim_inner_m);
        let m_arms = self.arm_count * lin * (r2 - r1);
        let i_arms = self.arm_count * lin * (r2.powi(3) - r1.powi(3)) / 3.0;
        let m_screw = rho_screw * PI * (0.5 * self.screw_diameter_m).powi(2) * self.screw_length_m;
        let m_screws = self.screw_count * m_screw;
        let i_screws = m_screws * self.screw_radius_m.powi(2);
        (m_rim + m_arms + m_screws + self.staff_mass_kg, i_rim + i_arms + i_screws)
    }

    /// Σ r² · (length across the flow) over arms and screws, for linear
    /// (Stokes-regime) drag, and Σ r³ · (frontal area) for form drag.
    fn drag_moments(&self) -> (f64, f64) {
        let (r1, r2) = (self.hub_radius_m, self.rim_inner_m);
        let lin_arms = self.arm_count * (r2.powi(3) - r1.powi(3)) / 3.0;
        let quad_arms = self.arm_count * self.arm_height_m * (r2.powi(4) - r1.powi(4)) / 4.0;
        let lin_screws = self.screw_count * self.screw_length_m * self.screw_radius_m.powi(2);
        let quad_screws =
            self.screw_count * self.screw_diameter_m * self.screw_length_m * self.screw_radius_m.powi(3);
        (lin_arms + lin_screws, quad_arms + quad_screws)
    }

    fn coil_gap_m(&self) -> f64 {
        let pitch = (self.hairspring_outer_m - self.hairspring_inner_m) / self.hairspring_coils;
        (pitch - self.hairspring_thickness_m).max(1e-6)
    }
}

/// Lumped barrel model of a preformed spiral mainspring.
#[derive(Clone, Copy, Debug)]
pub struct Mainspring {
    pub width_m: f64,
    pub thickness_m: f64,
    pub length_m: f64,
    pub barrel_radius_m: f64,
    pub arbor_radius_m: f64,
    /// Free-shape (reverse-curve) preform, as equivalent turns of bending
    /// already present when the spring lies against the barrel wall.
    pub preform_turns: f64,
    /// Turn scale over which the coils peel off the wall at let-down.
    pub peel_turns: f64,
    /// Inter-coil friction fraction (output torque M(1−f), winding M(1+f)).
    pub coil_friction: f64,
    /// Torque at which the slipping bridle lets the outer end slide.
    pub slip_torque_nm: f64,
    /// Turns wound up from the let-down state.
    pub turns: f64,
    /// Permanent set from over-tempering (0 = as made).
    pub set_fraction: f64,
}

impl Mainspring {
    #[must_use]
    pub fn pocket_tourbillon() -> Self {
        Self {
            width_m: 2.8e-3,
            thickness_m: 0.13e-3,
            length_m: 0.50,
            barrel_radius_m: 7.4e-3,
            arbor_radius_m: 1.9e-3,
            preform_turns: 3.0,
            peel_turns: 0.4,
            coil_friction: 0.08,
            slip_torque_nm: 0.018,
            turns: 0.0,
            set_fraction: 0.0,
        }
    }

    /// Bending stiffness E I / L [N m/rad] at temperature `t_k`.
    #[must_use]
    pub fn stiffness(&self, t_k: f64) -> f64 {
        let steel = &library().carbon_steel_hairspring;
        let e = steel.youngs_modulus(t_k);
        e * self.width_m * self.thickness_m.powi(3) / 12.0 / self.length_m
    }

    /// Static spring torque at `turns` of wind [N m].
    #[must_use]
    pub fn torque_at(&self, turns: f64, t_k: f64) -> f64 {
        let n = turns.max(0.0);
        let peel = 1.0 - (-n / self.peel_turns).exp();
        self.stiffness(t_k) * TAU * (n + self.preform_turns) * peel * (1.0 - self.set_fraction)
    }

    /// Turns at which the bridle starts to slip.
    #[must_use]
    pub fn slip_turns(&self, t_k: f64) -> f64 {
        let (mut lo, mut hi) = (0.0, self.geometric_turns());
        for _ in 0..50 {
            let mid = 0.5 * (lo + hi);
            if self.torque_at(mid, t_k) * (1.0 + self.coil_friction) < self.slip_torque_nm {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        lo
    }

    fn packed_turns(&self, inner: bool, length: f64) -> f64 {
        let t = self.thickness_m;
        let a = length * t / PI;
        if inner {
            ((self.arbor_radius_m.powi(2) + a).sqrt() - self.arbor_radius_m) / t
        } else {
            (self.barrel_radius_m - (self.barrel_radius_m.powi(2) - a).max(0.0).sqrt()) / t
        }
    }

    /// Turns between "all against the wall" and "all on the arbor".
    #[must_use]
    pub fn geometric_turns(&self) -> f64 {
        self.packed_turns(true, self.length_m) - self.packed_turns(false, self.length_m)
    }

    /// Radii of 18 coils from the arbor outwards [mm]: the spring length is
    /// split between a packed spiral on the arbor and one against the wall so
    /// that the coil count difference equals the turns wound.
    #[must_use]
    pub fn coil_radii_mm(&self) -> [f64; MAINSPRING_COILS] {
        let l = self.length_m;
        let n_wall_all = self.packed_turns(false, l);
        let target = self.turns.clamp(0.0, self.geometric_turns());
        let (mut lo, mut hi) = (0.0, l);
        for _ in 0..50 {
            let la = 0.5 * (lo + hi);
            let n = self.packed_turns(true, la) + self.packed_turns(false, l - la) - n_wall_all;
            if n < target {
                lo = la;
            } else {
                hi = la;
            }
        }
        let la = 0.5 * (lo + hi);
        let t = self.thickness_m;
        let mut out = [0.0; MAINSPRING_COILS];
        for (k, slot) in out.iter_mut().enumerate() {
            let s = (k as f64 + 0.5) / MAINSPRING_COILS as f64 * l;
            let r = if s <= la {
                (self.arbor_radius_m.powi(2) + s * t / PI).sqrt()
            } else {
                (self.barrel_radius_m.powi(2) - (l - s) * t / PI).max(0.0).sqrt()
            };
            *slot = r * 1.0e3;
        }
        out
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum WatchOperatingStatus {
    RunningChronometer = 0,
    RunningPerturbed = 1,
    StoppedMainspringUnwound = 2,
    JammedByMoltenGold = 3,
    LockedInSolidifiedGold = 4,
    SeizedByThermalOrCryoClearance = 5,
    PinnedByMagneticField = 6,
    FloodedStalledByOil = 7,
}

impl WatchOperatingStatus {
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::RunningChronometer => "RUNNING — CHRONOMETER CERTIFIED",
            Self::RunningPerturbed => "RUNNING — MULTIPHYSICS PERTURBED",
            Self::StoppedMainspringUnwound => "STOPPED — MAINSPRING UNWOUND",
            Self::JammedByMoltenGold => "STOPPED — FOULED BY LIQUID MOLTEN GOLD",
            Self::LockedInSolidifiedGold => "LOCKED — ENCASED IN SOLIDIFIED GOLD",
            Self::SeizedByThermalOrCryoClearance => "SEIZED — THERMAL / CRYOGENIC CLEARANCE LOSS",
            Self::PinnedByMagneticField => "STALLED — PALLET FORK MAGNETICALLY PINNED",
            Self::FloodedStalledByOil => "STALLED — VISCOUS OIL FLOOD OVERDAMPING",
        }
    }
}

/// Sound sources for the acoustic model, all from real contacts.
#[derive(Debug, Clone, Copy)]
pub enum AcousticEvent {
    /// An escapement contact impact (pallet/tooth, pin/fork, banking, guard).
    Contact {
        kind: ContactKind,
        is_lock: bool,
        lab_time_s: f64,
        dissipated_j: f64,
        impulse_n_s: f64,
    },
    /// The click dropping into the next ratchet tooth while winding.
    RatchetClick { lab_time_s: f64, energy_j: f64 },
    /// The bridle sliding round the barrel wall at full wind.
    BridleSlip { lab_time_s: f64, energy_j: f64 },
    /// A bench drop or knock on the case.
    Shock { lab_time_s: f64, g_peak: f64 },
}

/// Energy lost per balance cycle at the current amplitude, by mechanism [J].
#[derive(Debug, Clone, Copy, Default)]
pub struct LossBreakdown {
    pub air: f64,
    pub pivots: f64,
    pub pivot_oil: f64,
    pub hairspring: f64,
    pub eddy: f64,
    pub oil_flood: f64,
    pub liquid_gold: f64,
    pub solid_gold: f64,
    pub ice: f64,
    /// Static torques gripping the balance [N m].
    pub balance_solid_gold_nm: f64,
    pub balance_ice_nm: f64,
    /// Train-side resistances referred to the escape wheel [N m].
    pub train_liquid_gold_nm: f64,
    pub train_solid_gold_nm: f64,
    pub train_ice_nm: f64,
}

impl LossBreakdown {
    #[must_use]
    pub fn total(&self) -> f64 {
        self.air
            + self.pivots
            + self.pivot_oil
            + self.hairspring
            + self.eddy
            + self.oil_flood
            + self.liquid_gold
            + self.solid_gold
            + self.ice
    }
}

#[derive(Debug, Clone, Copy)]
struct ShockPulse {
    start_s: f64,
    duration_s: f64,
    peak_m_s2: f64,
    dir: [f64; 2],
}

#[derive(Debug, Clone, Copy, Default)]
struct BeatRecord {
    lab_time_s: f64,
    side: u8,
}

/// Fluid surrounding a part: density, viscosity.
#[derive(Debug, Clone, Copy)]
struct Fluid {
    rho: f64,
    mu: f64,
}

/// Factory regulation, done once at build time: the hairspring's active
/// length relative to its free (elastic-only) length that brings each alloy
/// [compensated, carbon steel] to zero rate in the factory position (20 °C,
/// pendant up, delivery wind). [`Horology::regulate`] finds it by running the
/// full escapement and measuring the rate it keeps; the test
/// `factory_regulation_constants_are_current` re-runs it and fails if the
/// escapement has changed and these need updating.
pub const REGULATION_FACTOR: [f64; 2] = [0.999_859_933, 0.999_863_198];

fn regulated_length(compensated: bool) -> f64 {
    let balance = BalanceDesign::pocket_tourbillon();
    Horology::free_length(&balance, Horology::alloy_for(compensated)) * REGULATION_FACTOR[usize::from(!compensated)]
}

/// Default state of wind on delivery, as a fraction of the slip turns.
pub const DELIVERY_WIND: f64 = 0.88;

#[derive(Debug, Clone)]
pub struct Horology {
    pub esc: Escapement,
    pub mainspring: Mainspring,
    pub balance: BalanceDesign,
    pub compensated: bool,
    hairspring_length_m: [f64; 2],
    balance_mass_kg: f64,
    balance_inertia_ref: f64,
    pub k_ref: f64,
    lab_origin_s: f64,
    psi_origin: f64,
    pub lab_time_s: f64,
    pub watch_displayed_time_s: f64,
    pub barrel_arbor_angle_rad: f64,
    pub barrel_wheel_angle_rad: f64,
    pub center_wheel_angle_rad: f64,
    pub third_wheel_angle_rad: f64,
    pub tourbillon_cage_angle_rad: f64,
    pub moonphase_angle_rad: f64,
    pub click_lever_angle_rad: f64,
    pub crown_wheel_angle_rad: f64,
    pub mainspring_coil_radii_mm: [f64; MAINSPRING_COILS],
    pub mainspring_torque_nm: f64,
    pub escape_drive_nm: f64,
    /// Factory oil remaining at the jewels (0 = burnt off).
    pub pivot_oil: f64,
    // Timing.
    beats: VecDeque<BeatRecord>,
    last_lock_s: f64,
    pub total_beats: u64,
    pub measured_rate_s_per_day: f64,
    pub rate_sigma_s_per_day: f64,
    rate_valid: bool,
    pub measured_beat_error_ms: f64,
    pub theoretical_rate_s_per_day: f64,
    pub balance_amplitude_deg: f64,
    swing_peak: f64,
    swing_dir: f64,
    last_peaks: [f64; 2],
    last_peak_s: f64,
    drain: Vec<f64>,
    impacts: Vec<ContactImpact>,
    // Per-beat escapement metrics.
    beat_drive_work0: f64,
    beat_impulse_work0: f64,
    pub last_impulse_energy_j: f64,
    pub escapement_efficiency: f64,
    pin_arc: Option<(f64, f64)>,
    pub measured_lift_deg: f64,
    // Loads and diagnosis.
    pub loads: EscapementLoads,
    pub losses: LossBreakdown,
    pub balance_q: f64,
    pub balance_pivot_friction_nm: f64,
    pub inertia_ratio: f64,
    pub stiffness_ratio: f64,
    pub pivot_friction_multiplier: f64,
    pivot_friction_ref: f64,
    pub min_clearance_um: f64,
    pub gold_jam_factor: f64,
    pub balance_aero_drag_s_inv: f64,
    pub rate_sensitivity_s_day_per_k: f64,
    pub status: WatchOperatingStatus,
    cosc: BettingEProcess,
    beats_since_cosc: u32,
    pub cosc_e_value: f64,
    shock: Option<ShockPulse>,
    shock_count: u32,
    /// What stopped the watch, kept while it stays stopped.
    stop_cause: Option<(WatchOperatingStatus, f64)>,
}

impl Horology {
    /// A regulated watch, wound to [`DELIVERY_WIND`], at 10:10:30.
    #[must_use]
    pub fn new() -> Self {
        let length = regulated_length(true);
        Self::build([length, f64::NAN], true)
    }

    fn alloy(&self) -> &'static HairspringAlloy {
        Self::alloy_for(self.compensated)
    }

    fn alloy_for(compensated: bool) -> &'static HairspringAlloy {
        let lib = library();
        if compensated { &lib.nivarox_hairspring } else { &lib.carbon_steel_hairspring }
    }

    /// Hairspring length giving `NOMINAL_FREQ_HZ` for the free balance at
    /// 20 °C: the starting guess before the escapement error is regulated out.
    fn free_length(balance: &BalanceDesign, alloy: &HairspringAlloy) -> f64 {
        let (_, i) = balance.mass_and_inertia();
        let w = TAU * NOMINAL_FREQ_HZ;
        let k = i * w * w;
        alloy.youngs_modulus(T_REF_K) * balance.hairspring_width_m * balance.hairspring_thickness_m.powi(3)
            / (12.0 * k)
    }

    /// Run the movement in the factory position and adjust the hairspring's
    /// active length (the regulator index) until the measured rate is zero.
    pub fn regulate(compensated: bool) -> f64 {
        let balance = BalanceDesign::pocket_tourbillon();
        let alloy = Self::alloy_for(compensated);
        let mut length = Self::free_length(&balance, alloy);
        let env = MovementEnvironment::default();
        let mut events = Vec::new();
        for _ in 0..3 {
            let mut h = Self::build([length, length], compensated);
            for _ in 0..140 {
                h.step(0.05, &env, &mut events);
                events.clear();
            }
            let (rate, _) = h.fit_rate(true).unwrap_or((0.0, 0.0));
            // ω ∝ L^(-1/2): a watch gaining r s/day needs L·(1 + 2r/86400).
            length *= 1.0 + 2.0 * rate / SECONDS_PER_DAY;
            if rate.abs() < 0.05 {
                break;
            }
        }
        length
    }

    fn build(lengths: [f64; 2], compensated: bool) -> Self {
        let balance = BalanceDesign::pocket_tourbillon();
        let (mass, inertia) = balance.mass_and_inertia();
        let mut mainspring = Mainspring::pocket_tourbillon();
        mainspring.turns = DELIVERY_WIND * mainspring.slip_turns(T_REF_K);
        let alloy = Self::alloy_for(compensated);
        let li = if compensated { 0 } else { 1 };
        let k_ref = alloy.youngs_modulus(T_REF_K) * balance.hairspring_width_m
            * balance.hairspring_thickness_m.powi(3)
            / (12.0 * lengths[li]);
        let start = 10.0 * 3600.0 + 10.0 * 60.0 + 30.0;
        let esc = Escapement::new(EscapementDesign::swiss_lever_pocket(), 280.0 * DEG, 0.0);
        let psi0 = esc.q[2];
        let mut h = Self {
            esc,
            mainspring,
            balance,
            compensated,
            hairspring_length_m: lengths,
            balance_mass_kg: mass,
            balance_inertia_ref: inertia,
            k_ref,
            lab_origin_s: start,
            psi_origin: psi0,
            lab_time_s: start,
            watch_displayed_time_s: start,
            barrel_arbor_angle_rad: 0.0,
            barrel_wheel_angle_rad: 0.0,
            center_wheel_angle_rad: 0.0,
            third_wheel_angle_rad: 0.0,
            tourbillon_cage_angle_rad: 0.0,
            moonphase_angle_rad: 0.35 * TAU,
            click_lever_angle_rad: 0.0,
            crown_wheel_angle_rad: 0.0,
            mainspring_coil_radii_mm: [0.0; MAINSPRING_COILS],
            mainspring_torque_nm: 0.0,
            escape_drive_nm: 0.0,
            pivot_oil: 1.0,
            beats: VecDeque::with_capacity(96),
            last_lock_s: start,
            total_beats: 0,
            measured_rate_s_per_day: 0.0,
            rate_sigma_s_per_day: 0.0,
            rate_valid: false,
            measured_beat_error_ms: 0.0,
            theoretical_rate_s_per_day: 0.0,
            balance_amplitude_deg: 280.0,
            swing_peak: 0.0,
            swing_dir: 0.0,
            last_peaks: [280.0 * DEG; 2],
            last_peak_s: start,
            drain: Vec::with_capacity(256),
            impacts: Vec::with_capacity(64),
            beat_drive_work0: 0.0,
            beat_impulse_work0: 0.0,
            last_impulse_energy_j: 0.0,
            escapement_efficiency: 0.0,
            pin_arc: None,
            measured_lift_deg: 0.0,
            loads: EscapementLoads::default(),
            losses: LossBreakdown::default(),
            balance_q: 0.0,
            balance_pivot_friction_nm: 0.0,
            inertia_ratio: 1.0,
            stiffness_ratio: 1.0,
            pivot_friction_multiplier: 1.0,
            pivot_friction_ref: 0.0,
            min_clearance_um: 0.0,
            gold_jam_factor: 0.0,
            balance_aero_drag_s_inv: 0.0,
            rate_sensitivity_s_day_per_k: 0.0,
            status: WatchOperatingStatus::RunningChronometer,
            cosc: BettingEProcess::new(0.25),
            beats_since_cosc: 0,
            cosc_e_value: 1.0,
            shock: None,
            shock_count: 0,
            stop_cause: None,
        };
        let env = MovementEnvironment::default();
        h.loads = h.assemble_loads(&env);
        h.pivot_friction_ref = h.balance_pivot_friction_nm.max(1e-15);
        h.update_kinematics();
        h
    }

    /// Fit `t_k = a + T k + d s_k` (s_k = ±1 by pallet side, so a beat error
    /// does not leak into the residual) to the recorded lock drops. Returns
    /// the rate (s/day) and its standard error. With `skip_transient`, beats
    /// in the first 2 s of the record are ignored.
    fn fit_rate(&self, skip_transient: bool) -> Option<(f64, f64)> {
        let t_first = self.beats.front()?.lab_time_s;
        let pts: Vec<(f64, f64, f64)> = self
            .beats
            .iter()
            .enumerate()
            .filter(|(_, b)| !skip_transient || b.lab_time_s - t_first > 2.0)
            .map(|(k, b)| (k as f64, if b.side == 0 { 1.0 } else { -1.0 }, b.lab_time_s))
            .collect();
        if pts.len() < 8 {
            return None;
        }
        let t0 = pts[0].2;
        // Normal equations for the basis [1, k, s].
        let mut a = [[0.0_f64; 3]; 3];
        let mut rhs = [0.0_f64; 3];
        for &(k, sd, t) in &pts {
            let basis = [1.0, k, sd];
            for i in 0..3 {
                rhs[i] += basis[i] * (t - t0);
                for j in 0..3 {
                    a[i][j] += basis[i] * basis[j];
                }
            }
        }
        let det3 = |m: &[[f64; 3]; 3]| {
            m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
        };
        let d = det3(&a);
        if d.abs() < 1e-12 {
            return None;
        }
        let solve = |col: usize| {
            let mut m = a;
            for i in 0..3 {
                m[i][col] = rhs[i];
            }
            det3(&m) / d
        };
        let (c0, slope, c2) = (solve(0), solve(1), solve(2));
        if slope <= 0.0 {
            return None;
        }
        let n = pts.len() as f64;
        let res: f64 = pts
            .iter()
            .map(|&(k, sd, t)| (t - t0) - (c0 + slope * k + c2 * sd))
            .map(|r| r * r)
            .sum();
        let sigma_t = (res / (n - 3.0).max(1.0)).sqrt();
        // Var(slope) = σ² (AᵀA)⁻¹₁₁, the (1,1) cofactor over the determinant.
        let cof11 = a[0][0] * a[2][2] - a[0][2] * a[2][0];
        let sigma_slope = sigma_t * (cof11 / d).max(0.0).sqrt();
        let rate = (NOMINAL_BEAT_PERIOD_S / slope - 1.0) * SECONDS_PER_DAY;
        let sigma = NOMINAL_BEAT_PERIOD_S / (slope * slope) * sigma_slope * SECONDS_PER_DAY;
        Some((rate, sigma))
    }

    /// Wind the crown: `delta_turns` of the barrel arbor (positive winds).
    pub fn wind_crown(&mut self, delta_turns: f64, env: &MovementEnvironment, events: &mut Vec<AcousticEvent>) {
        if delta_turns <= 0.0 {
            return;
        }
        // Frozen gold round the barrel and keyless works blocks the crown.
        if env.barrel.solid_gold_m > 0.4e-3 {
            return;
        }
        const RATCHET_TEETH: f64 = 30.0;
        let t_k = env.barrel.temperature_k;
        let max_turns = self.mainspring.slip_turns(t_k);
        let before = self.mainspring.turns;
        let after = (before + delta_turns).min(max_turns.max(before));
        let slipped = before + delta_turns - after;
        self.mainspring.turns = after;
        let teeth_before = (self.barrel_arbor_angle_rad / TAU * RATCHET_TEETH).floor();
        self.barrel_arbor_angle_rad += delta_turns * TAU;
        self.crown_wheel_angle_rad -= delta_turns * TAU * (78.0 / 42.0);
        let teeth_after = (self.barrel_arbor_angle_rad / TAU * RATCHET_TEETH).floor();
        let clicks = (teeth_after - teeth_before).max(0.0) as usize;
        // Each click: the click spring's stored energy, ½ k δ², released into
        // the ratchet tooth (click spring 2 mN m/rad, 0.12 rad lift).
        let e_click = 0.5 * 2.0e-3 * 0.12 * 0.12;
        for i in 0..clicks.min(64) {
            events.push(AcousticEvent::RatchetClick {
                lab_time_s: self.lab_time_s + i as f64 * 0.004,
                energy_j: e_click,
            });
        }
        if slipped > 1e-6 {
            events.push(AcousticEvent::BridleSlip {
                lab_time_s: self.lab_time_s,
                energy_j: self.mainspring.slip_torque_nm * slipped * TAU * 0.02,
            });
        }
        self.click_lever_angle_rad = 0.12 * (self.barrel_arbor_angle_rad / TAU * RATCHET_TEETH).fract();
        self.update_kinematics();
    }

    /// Fit the alternative hairspring (compensated Nivarox or carbon steel),
    /// regulated in the factory position like the original.
    pub fn set_compensated(&mut self, compensated: bool) {
        if self.compensated == compensated {
            return;
        }
        self.compensated = compensated;
        let li = if compensated { 0 } else { 1 };
        if self.hairspring_length_m[li].is_nan() {
            self.hairspring_length_m[li] = regulated_length(compensated);
        }
        self.k_ref = self.alloy().youngs_modulus(T_REF_K) * self.balance.hairspring_width_m
            * self.balance.hairspring_thickness_m.powi(3)
            / (12.0 * self.hairspring_length_m[li]);
    }

    /// A bench drop of peak acceleration `g_peak` (in g): a 0.3 ms half-sine
    /// in the plane of the movement plus the case rotation it starts.
    pub fn apply_mechanical_shock(&mut self, g_peak: f64, events: &mut Vec<AcousticEvent>) {
        let a = g_peak * G0_M_S2;
        let duration = 3.0e-4;
        self.shock_count += 1;
        let ang = PI / 2.0 + 0.35 * (self.shock_count as f64 * 2.399_963).sin();
        self.shock = Some(ShockPulse {
            start_s: self.lab_time_s,
            duration_s: duration,
            peak_m_s2: a,
            dir: [ang.cos(), ang.sin()],
        });
        // An off-centre landing spins the case: Δω ≈ 0.25 Δv / r_case. Every
        // body pivoted in the movement keeps its absolute speed, so relative
        // to the movement it gains −Δω.
        let dv = 2.0 / PI * a * duration;
        let sign = if self.shock_count % 2 == 0 { 1.0 } else { -1.0 };
        let d_omega = sign * 0.25 * dv / 0.025;
        self.esc.v[0] -= d_omega;
        self.esc.v[1] -= d_omega;
        self.esc.asleep = false;
        self.esc.frozen = false;
        events.push(AcousticEvent::Shock { lab_time_s: self.lab_time_s, g_peak });
    }

    /// Advance the movement by `dt` in the environment `env`.
    pub fn step(&mut self, dt: f64, env: &MovementEnvironment, events: &mut Vec<AcousticEvent>) {
        self.age_materials(env, dt);
        let base = self.assemble_loads(env);
        self.loads = base;
        let mut remaining = dt;
        while remaining > 1e-12 {
            let in_shock = self.shock.is_some_and(|s| self.lab_time_s - s.start_s < s.duration_s);
            let h = (if in_shock { 2.0e-5_f64 } else { 1.0e-3 }).min(remaining);
            let mut l = base;
            if in_shock {
                self.add_shock_loads(&mut l, env);
            } else if self.shock.is_some() {
                self.shock = None;
            }
            let psi0 = self.esc.q[2];
            self.impacts.clear();
            self.esc.advance(h, &l, &mut self.impacts);
            remaining -= h;
            self.lab_time_s = self.lab_origin_s + self.esc.time_s;
            // The barrel lets down as the escape wheel advances.
            let d_psi = self.esc.q[2] - psi0;
            self.mainspring.turns = (self.mainspring.turns - d_psi / (TAU * train::RATIO)).max(0.0);
            self.track_pin_arc();
            self.process_impacts(events);
            self.track_swing();
        }
        self.update_measurements();
        self.update_kinematics();
        self.diagnose(env);
    }

    fn add_shock_loads(&self, l: &mut EscapementLoads, env: &MovementEnvironment) {
        let Some(s) = self.shock else { return };
        let t = self.lab_time_s - s.start_s;
        let a = s.peak_m_s2 * (PI * t / s.duration_s).sin().max(0.0);
        // Inertial acceleration in the watch frame, rotated into the cage.
        let g_eff = [env.gravity_m_s2[0] - a * s.dir[0], env.gravity_m_s2[1] - a * s.dir[1]];
        let (pa, ph) = self.poise_terms(g_eff);
        l.balance_poise_amp = pa;
        l.balance_poise_phase = ph;
        // Lever: 8 mg with its centre of mass 0.05 mm along the fork.
        let me_l = 8.0e-6 * 0.05e-3;
        let c = -self.tourbillon_cage_angle_rad;
        let gx = g_eff[0] * c.cos() - g_eff[1] * c.sin();
        let gy = g_eff[0] * c.sin() + g_eff[1] * c.cos();
        let phi = self.esc.q[1];
        let (rx, ry) = (-phi.sin() * me_l, phi.cos() * me_l);
        l.lever_torque += rx * gy - ry * gx;
    }

    /// Out-of-poise torque `a sin(θ + phase)` for in-plane gravity `g` (watch
    /// frame), seen in the rotating cage.
    fn poise_terms(&self, g: [f64; 2]) -> (f64, f64) {
        let g_mag = (g[0] * g[0] + g[1] * g[1]).sqrt();
        if g_mag < 1e-9 {
            return (0.0, 0.0);
        }
        // Gravity direction in the cage frame (cage turns by −angle on screen).
        let beta = g[1].atan2(g[0]) + self.tourbillon_cage_angle_rad;
        (-self.balance.poise_kg_m * g_mag, self.balance.poise_angle_rad - beta)
    }

    /// Slow material changes: factory oil burning off, mainspring tempering.
    fn age_materials(&mut self, env: &MovementEnvironment, dt: f64) {
        let t = env.balance.temperature_k.max(env.escapement.temperature_k);
        // Pennzane evaporation/decomposition (activation 130 kJ/mol).
        let k_evap = 3.0e9 * (-1.3e5 / (R_GAS * t)).exp();
        self.pivot_oil *= (-k_evap * dt).exp();
        // Stress relaxation of the hardened spring steel (activation 180 kJ/mol):
        // a blued mainspring held above ~700 K loses its set.
        let tb = env.barrel.temperature_k;
        let k_set = 1.0e10 * (-1.8e5 / (R_GAS * tb)).exp();
        let s = self.mainspring.set_fraction;
        self.mainspring.set_fraction = (s + (0.9 - s) * (1.0 - (-k_set * dt).exp())).min(0.9);
    }

    /// Recompute the derived readouts (rate sensitivity, Q, aero drag, …)
    /// for `env` without advancing the oscillator.
    pub fn refresh_derived(&mut self, env: &MovementEnvironment) {
        let _ = self.assemble_loads(env);
    }

    /// Relubricate the jewels after the excess oil is cleaned out.
    pub fn relubricate(&mut self) {
        self.pivot_oil = 1.0;
    }

    /// Hairspring stiffness [N m/rad] at temperature `t_k` (no field, no oil).
    fn hairspring_k(&self, t_k: f64) -> f64 {
        let alloy = self.alloy();
        let dt = t_k - T_REF_K;
        // E(T) and the strip's b h³ / L growing as (1 + α ΔT)³.
        self.k_ref * alloy.youngs_modulus(t_k) / alloy.youngs_modulus(T_REF_K) * (1.0 + alloy.alpha * dt).powi(3)
    }

    fn balance_inertia(&self, t_k: f64) -> f64 {
        let a = library().becu.alpha;
        self.balance_inertia_ref * (1.0 + a * (t_k - T_REF_K)).powi(2)
    }

    /// Assemble every load on the three escapement bodies from the part
    /// geometry and the environment.
    fn assemble_loads(&mut self, env: &MovementEnvironment) -> EscapementLoads {
        let lib = library();
        let b = self.balance;
        let tb = env.balance.temperature_k;
        let i_bal = self.balance_inertia(tb);
        let alloy = *self.alloy();

        // --- Hairspring stiffness: temperature, magnetism, oil and ice bridges.
        let mut k = self.hairspring_k(tb);
        let m = env.hairspring_magnetisation.clamp(-1.0, 1.0);
        k *= 1.0 + alloy.magnetic_stiffening_at_sat * m * m;
        let gap = b.coil_gap_m();
        let spring_wet = |depth: f64| ((depth - b.hairspring_bottom_z_m) / b.hairspring_width_m).clamp(0.0, 1.0);
        let oil_bridge = 0.5 * spring_wet(env.balance.oil_m);
        let ice_bridge = ((2.0 * env.balance.frost_m - 0.5 * gap) / (0.5 * gap)).clamp(0.0, 1.0);
        let gold_bridge = spring_wet(env.balance.liquid_gold_m + env.balance.solid_gold_m);
        let bridged = (oil_bridge + ice_bridge + gold_bridge).min(0.98);
        k /= 1.0 - bridged;
        let w0 = (k / i_bal).sqrt();

        // Amplitude for the amplitude-dependent terms (Reynolds number, losses).
        let amp = (self.balance_amplitude_deg * DEG).max(30.0 * DEG);

        // --- Fluids round the balance.
        let air = air_state(env.air_temperature_k, env.air_pressure_pa).map(|g| Fluid {
            rho: g.density,
            mu: g.dynamic_viscosity,
        });
        let gold = &lib.gold;
        let gold_f = Fluid {
            rho: gold.state(gold.enthalpy_at(env.balance.gold_temperature_k)).density,
            mu: gold.effective_viscosity(env.balance.gold_temperature_k, env.balance.gold_liquid_fraction),
        };
        let oil_f = Fluid {
            rho: lib.oil.density(tb),
            mu: lib.oil.viscosity(tb),
        };
        let rim_wet = |depth: f64| ((depth - b.rim_bottom_z_m) / b.rim_height_m).clamp(0.0, 1.0);
        let below_wet = |depth: f64| if depth > b.rim_bottom_z_m - b.gap_below_m { 1.0 } else { 0.0 };
        let above_wet = |depth: f64| {
            if depth > b.rim_bottom_z_m + b.rim_height_m { 1.0 } else { 0.0 }
        };
        let face = PI * (b.rim_outer_m.powi(4) - b.rim_inner_m.powi(4)) / 2.0
            + b.arm_count * b.arm_width_m * (b.rim_inner_m.powi(4) - b.hub_radius_m.powi(4)) / 4.0;
        let (lin_m, quad_m) = b.drag_moments();
        // Linear (viscous) and quadratic (form) coefficients for a fluid that
        // covers the lower face, the rim and the upper face by the given
        // fractions. Face shear uses the gap or the Stokes layer, whichever is
        // thinner; the rim is a Couette journal against the cage band.
        let balance_drag = |f: Fluid, below: f64, side: f64, above: f64| -> (f64, f64) {
            if below + side + above <= 0.0 {
                return (0.0, 0.0);
            }
            let delta = (2.0 * f.mu / (f.rho * w0)).sqrt();
            let mut c = f.mu * face * (below / b.gap_below_m.min(delta) + above / b.gap_above_m.min(delta));
            if side > 0.0 {
                let journal = JournalGeometry {
                    inner_radius_m: b.rim_outer_m,
                    outer_radius_m: b.rim_outer_m + b.gap_side_m.min(delta),
                    length_m: side * b.rim_height_m,
                };
                c += journal.resistance_factor_m3().unwrap_or(0.0) * f.mu;
            }
            // Arms and screws: Stokes-regime drag ~6 μ per unit length plus
            // form drag with C_D ≈ 1.2.
            let wet = side.max(below.min(1.0) * 0.5);
            c += 6.0 * f.mu * lin_m * wet;
            let cq = 0.6 * f.rho * quad_m * wet;
            (c, cq)
        };

        let mut losses = LossBreakdown::default();
        let per_cycle = |c: f64, cq: f64, tau: f64| -> f64 {
            PI * c * w0 * amp * amp + 8.0 / 3.0 * cq * w0 * w0 * amp.powi(3) + 4.0 * tau * amp
        };

        // Air fills whatever the liquids leave.
        let liquid_top = env.balance.oil_m.max(env.balance.liquid_gold_m).max(env.balance.solid_gold_m);
        let (mut c_bal, mut cq_bal) = (0.0, 0.0);
        if let Some(a) = air {
            let (c, cq) = balance_drag(
                a,
                1.0 - below_wet(liquid_top),
                1.0 - rim_wet(liquid_top),
                1.0 - above_wet(liquid_top),
            );
            // Crossflow from the air jet linearises its form drag:
            // d/dv ½ρC_D A (u+v)|u+v| = ρ C_D A |u|, averaged over the rim.
            let c_jet = 2.0 / PI * 1.2 * a.rho * env.air_jet_speed_m_s
                * (b.screw_count * b.screw_diameter_m * b.screw_length_m)
                * b.screw_radius_m.powi(2);
            c_bal += c + c_jet;
            cq_bal += cq;
            losses.air = per_cycle(c + c_jet, cq, 0.0);
        }
        // Liquid gold (mobile part) round the balance.
        let gold_depth = env.balance.liquid_gold_m;
        if gold_depth > 0.0 {
            let (c, cq) = balance_drag(gold_f, below_wet(gold_depth), rim_wet(gold_depth), above_wet(gold_depth));
            c_bal += c;
            cq_bal += cq;
            losses.liquid_gold = per_cycle(c, cq, 0.0);
        }
        // Poured oil round the balance.
        let oil_depth = env.balance.oil_m;
        if oil_depth > 0.0 {
            let (c, cq) = balance_drag(oil_f, below_wet(oil_depth), rim_wet(oil_depth), above_wet(oil_depth));
            c_bal += c;
            cq_bal += cq;
            losses.oil_flood = per_cycle(c, cq, 0.0);
        }

        // --- Jewelled pivots: Coulomb friction under the dead weight, Petroff
        // shear of the factory oil in the journal.
        let g = env.gravity_m_s2;
        let g_in = (g[0] * g[0] + g[1] * g[1]).sqrt();
        let mu_j = b.jewel_mu_dry + (b.jewel_mu_oiled - b.jewel_mu_dry) * self.pivot_oil;
        let w_radial = self.balance_mass_kg * g_in;
        let w_axial = self.balance_mass_kg * g[2].abs();
        // Hertz contact of the rounded pivot end on the sapphire endstone.
        let steel = lib.steel;
        let sap = lib.sapphire;
        let e_star = 1.0 / ((1.0 - steel.nu * steel.nu) / steel.e + (1.0 - sap.nu * sap.nu) / sap.e);
        let a_hertz = (3.0 * w_axial * b.pivot_radius_m / (4.0 * e_star)).cbrt();
        let tau_pivot = mu_j * w_radial * b.pivot_radius_m + 3.0 * PI / 16.0 * mu_j * w_axial * a_hertz;
        self.balance_pivot_friction_nm = tau_pivot;
        losses.pivots = per_cycle(0.0, 0.0, tau_pivot);
        let clearance = b.pivot_clearance_m
            + b.pivot_radius_m * (sap.alpha - steel.alpha) * (tb - T_REF_K);
        let c_petroff = if clearance > 0.0 {
            2.0 * JournalGeometry {
                inner_radius_m: b.pivot_radius_m,
                outer_radius_m: b.pivot_radius_m + clearance,
                length_m: b.pivot_length_m,
            }
            .resistance_factor_m3()
            .unwrap_or(0.0)
                * lib.oil.viscosity(tb)
                * self.pivot_oil
        } else {
            0.0
        };
        c_bal += c_petroff;
        losses.pivot_oil = per_cycle(c_petroff, 0.0, 0.0);

        // Internal friction of the strip: loss tangent η gives c = k η / ω.
        let c_hs = k * alloy.loss_tangent / w0;
        c_bal += c_hs;
        losses.hairspring = per_cycle(c_hs, 0.0, 0.0);

        // Eddy currents in the Glucydur rim moving through a field gradient.
        let r_rim = 0.5 * (b.rim_outer_m + b.rim_inner_m);
        let c_eddy = 0.5 * b.rim_conductivity_s_m * b.rim_volume() * r_rim.powi(2)
            * (r_rim * env.grad_b_balance_t_m).powi(2);
        c_bal += c_eddy;
        losses.eddy = per_cycle(c_eddy, 0.0, 0.0);

        // --- Solids gripping the balance: frozen gold and ice bridges.
        let tau_gold_yield = gold_shear_yield_pa(env.balance.gold_temperature_k, env.balance.gold_liquid_fraction);
        let solid_depth = env.balance.solid_gold_m;
        let grip = |wet_below: f64, wet_side: f64| -> f64 {
            let a_face = PI * (b.rim_outer_m.powi(2) - b.rim_inner_m.powi(2));
            wet_below * a_face * r_rim + wet_side * TAU * b.rim_outer_m * b.rim_height_m * b.rim_outer_m
        };
        let tau_solid_gold = if solid_depth > b.rim_bottom_z_m - b.gap_below_m {
            let solid_yield = gold_shear_yield_pa(env.balance.temperature_k, 0.0);
            solid_yield * grip(1.0, rim_wet(solid_depth))
        } else {
            0.0
        };
        // Mushy gold above its coherency point also carries a yield stress.
        let tau_mushy = if gold_depth > b.rim_bottom_z_m - b.gap_below_m {
            tau_gold_yield * grip(below_wet(gold_depth), rim_wet(gold_depth))
        } else {
            0.0
        };
        // Frost growing from both faces bridges the gap below the rim when 2s > g.
        let ice_fill = ((2.0 * env.balance.frost_m - b.gap_below_m) / (0.1 * b.gap_below_m)).clamp(0.0, 1.0);
        let tau_ice = ice_fill * ICE_SHEAR_PA * grip(1.0, 0.0) * 0.1;
        let tau_bal = tau_pivot + tau_solid_gold + tau_mushy + tau_ice;
        losses.solid_gold = per_cycle(0.0, 0.0, tau_solid_gold + tau_mushy);
        losses.ice = per_cycle(0.0, 0.0, tau_ice);
        losses.balance_solid_gold_nm = tau_solid_gold + tau_mushy;
        losses.balance_ice_nm = tau_ice;

        // --- Poise.
        let (pa, ph) = self.poise_terms([g[0], g[1]]);

        // --- Going train and escape wheel.
        let t_barrel = env.barrel.temperature_k;
        let ms = self.mainspring.torque_at(self.mainspring.turns, t_barrel);
        self.mainspring_torque_nm = ms;
        let eta = train_efficiency(train::MESH_MU);
        let drive = ms * (1.0 - self.mainspring.coil_friction) * eta / train::RATIO;
        self.escape_drive_nm = drive;
        // Pivot friction of the escape wheel and of the cage (referred down by
        // the cage ratio).
        let mu_p = 0.15;
        let esc_pivots = mu_p * 12.0e-6 * g_in * 0.05e-3 + mu_p * 0.35e-3 * g_in * 0.1e-3 / train::CAGE_TO_ESCAPE;
        let lever_pivots = mu_p * 8.0e-6 * g_in * 0.035e-3;
        // Fluids and solids in the train: escape wheel (ratio 1), cage (1/10),
        // third wheel (1/75) — torque refers to the escape wheel × ratio,
        // viscous × ratio², form drag × ratio³.
        let rotors = [
            Rotor { radius: 2.4e-3, z0: 1.0e-3, thick: 0.15e-3, teeth: 15.0, ratio: 1.0, part: env.escapement },
            Rotor {
                radius: 7.0e-3,
                z0: 0.3e-3,
                thick: 0.25e-3,
                teeth: 0.0,
                ratio: 1.0 / train::CAGE_TO_ESCAPE,
                part: env.train,
            },
            Rotor {
                radius: 5.0e-3,
                z0: 0.2e-3,
                thick: 0.3e-3,
                teeth: train::THIRD,
                ratio: 1.0 / (train::CAGE_TO_ESCAPE * train::THIRD_TO_CAGE),
                // Under the dial: only its temperature reaches it.
                part: PartField { temperature_k: env.train.temperature_k, ..PartField::default() },
            },
        ];
        let (mut esc_c, mut esc_cq, mut esc_yield_gold, mut esc_mushy, mut esc_ice) = (1.0e-14, 0.0, 0.0, 0.0, 0.0);
        for r in rotors {
            let p = r.part;
            let gold_mu = gold.effective_viscosity(p.gold_temperature_k, p.gold_liquid_fraction);
            let gold_rho = gold.state(gold.enthalpy_at(p.gold_temperature_k)).density;
            let oil_mu = lib.oil.viscosity(p.temperature_k);
            let oil_rho = lib.oil.density(p.temperature_k);
            for (depth, mu, rho) in [(p.liquid_gold_m, gold_mu, gold_rho), (p.oil_m, oil_mu, oil_rho)] {
                let wet = ((depth - r.z0) / r.thick).clamp(0.0, 1.0);
                if depth <= r.z0 - 0.2e-3 {
                    continue;
                }
                // Spoked wheel: half of the disc area shears a 0.2 mm film.
                let c = mu * 0.5 * PI * r.radius.powi(4) / (2.0 * 0.2e-3) * (1.0 + wet);
                let a_teeth = r.teeth.max(6.0) * 0.3e-3 * r.thick * wet;
                let cq = 0.6 * rho * a_teeth * r.radius.powi(3);
                esc_c += c * r.ratio * r.ratio;
                esc_cq += cq * r.ratio.powi(3);
            }
            let grip_area = TAU * r.radius * r.thick + 0.5 * PI * r.radius.powi(2);
            if p.solid_gold_m > r.z0 {
                let fill = ((p.solid_gold_m - r.z0) / r.thick).clamp(0.2, 1.0);
                esc_yield_gold += gold_shear_yield_pa(p.temperature_k, 0.0) * grip_area * fill * r.radius * r.ratio;
            }
            if p.liquid_gold_m > r.z0 {
                let fill = ((p.liquid_gold_m - r.z0) / r.thick).clamp(0.2, 1.0);
                esc_mushy += gold_shear_yield_pa(p.gold_temperature_k, p.gold_liquid_fraction)
                    * grip_area
                    * fill
                    * r.radius
                    * r.ratio;
            }
            let ice = ((2.0 * p.frost_m - 0.1e-3) / 0.02e-3).clamp(0.0, 1.0);
            esc_ice += ice * ICE_SHEAR_PA * grip_area * 0.1 * r.radius * r.ratio;
        }
        // Liquid resistance at a typical escape-wheel speed during impulse
        // and drop (12° in ~10 ms ≈ 20 rad/s).
        let w_esc = 20.0;
        losses.train_liquid_gold_nm = esc_mushy + esc_c * w_esc + esc_cq * w_esc * w_esc;
        losses.train_solid_gold_nm = esc_yield_gold;
        losses.train_ice_nm = esc_ice;
        self.losses = losses;

        // Pallet friction from the state of the oil on the stones.
        let te = env.escapement.temperature_k;
        let oil_visc = lib.oil.viscosity(te);
        let (mu_s, mu_k) = if oil_visc > 1.0e5 {
            (0.45, 0.3) // glassy oil glues the stones
        } else if self.pivot_oil < 0.2 {
            (0.28, 0.2)
        } else {
            (0.17, 0.12)
        };

        // Derived quantities for the panels.
        let e_cycle = 0.5 * k * amp * amp;
        let total = losses.total().max(1e-30);
        self.balance_q = TAU * e_cycle / total;
        // Air's share of the decay rate γ = c_eq / 2I, with c_eq the viscous
        // coefficient that would dissipate the same energy per cycle.
        self.balance_aero_drag_s_inv = losses.air / (PI * w0 * amp * amp) / (2.0 * i_bal);
        self.inertia_ratio = i_bal / self.balance_inertia_ref;
        self.stiffness_ratio = k / self.k_ref;
        self.pivot_friction_multiplier = tau_pivot / self.pivot_friction_ref.max(1e-15);
        let frost_gap = b.gap_below_m - 2.0 * env.balance.frost_m;
        let coil_gap = gap - 2.0 * env.balance.frost_m;
        self.min_clearance_um = clearance.min(frost_gap).min(coil_gap) * 1.0e6;
        let train_block = esc_yield_gold + esc_mushy + esc_ice;
        self.gold_jam_factor = ((esc_yield_gold + esc_mushy + esc_cq) / (drive + 1e-15))
            .max((tau_solid_gold + tau_mushy) / (k * 30.0 * DEG))
            .max(losses.liquid_gold / (losses.total() + 1e-30))
            .clamp(0.0, 1.0);
        let w_nom = TAU * NOMINAL_FREQ_HZ;
        self.theoretical_rate_s_per_day = (w0 / w_nom - 1.0) * SECONDS_PER_DAY;
        // d(rate)/dT = ½ (d ln E/dT + 3 α_spring − 2 α_balance) · 86400.
        let dt = tb - T_REF_K;
        let dln_e = (alloy.thermoelastic_1 + 2.0 * alloy.thermoelastic_2 * dt)
            / (1.0 + alloy.thermoelastic_1 * dt + alloy.thermoelastic_2 * dt * dt);
        let a_b = lib.becu.alpha;
        self.rate_sensitivity_s_day_per_k = 0.5
            * (dln_e + 3.0 * alloy.alpha / (1.0 + alloy.alpha * dt) - 2.0 * a_b / (1.0 + a_b * dt))
            * SECONDS_PER_DAY
            * (w0 / w_nom);

        EscapementLoads {
            i_balance: i_bal,
            k_hairspring: k,
            theta_eq: 0.0,
            balance_viscous: c_bal,
            balance_quadratic: cq_bal,
            balance_coulomb: tau_bal,
            balance_poise_amp: pa,
            balance_poise_phase: ph,
            lever_torque: env.lever_magnetic_torque_nm,
            lever_coulomb: lever_pivots,
            lever_viscous: 1.0e-14,
            escape_drive: drive,
            escape_viscous: esc_c,
            escape_coulomb: esc_pivots + train_block,
            pallet_mu_static: mu_s,
            pallet_mu_kinetic: mu_k,
            pin_mu: 0.15,
        }
    }

    fn track_pin_arc(&mut self) {
        if self.esc.active_mask & (1 << ContactKind::PinFork as u32) != 0 {
            let th = self.esc.q[0];
            self.pin_arc = Some(match self.pin_arc {
                None => (th, th),
                Some((a, _)) => (a, th),
            });
        }
    }

    fn process_impacts(&mut self, events: &mut Vec<AcousticEvent>) {
        for im in self.impacts.drain(..) {
            let lab = self.lab_origin_s + im.time_s;
            events.push(AcousticEvent::Contact {
                kind: im.kind,
                is_lock: im.is_lock,
                lab_time_s: lab,
                dissipated_j: im.dissipated_j,
                impulse_n_s: im.normal_impulse_n_s,
            });
            let is_drop = im.is_lock && matches!(im.kind, ContactKind::EntryStoneTooth | ContactKind::ExitStoneTooth);
            if !is_drop || lab - self.last_lock_s < 0.02 {
                continue;
            }
            let side = u8::from(matches!(im.kind, ContactKind::ExitStoneTooth));
            // A stall, or two drops on the same stone (a tripped escapement),
            // breaks the beat train: the rate fit restarts from here.
            let gap = lab - self.last_lock_s;
            if gap > 3.0 * NOMINAL_BEAT_PERIOD_S || self.beats.back().is_some_and(|b| b.side == side) {
                self.beats.clear();
            }
            self.last_lock_s = lab;
            self.total_beats += 1;
            self.beats.push_back(BeatRecord { lab_time_s: lab, side });
            while self.beats.len() > 60 {
                self.beats.pop_front();
            }
            if self.drain.len() < 3 * 512 {
                self.drain.extend_from_slice(&[lab, f64::from(side), im.dissipated_j * 1.0e9]);
            }
            // Per-beat energetics.
            let e_imp = self.esc.impulse_work_on_balance_j - self.beat_impulse_work0;
            let e_drv = self.esc.drive_work_j - self.beat_drive_work0;
            self.beat_impulse_work0 = self.esc.impulse_work_on_balance_j;
            self.beat_drive_work0 = self.esc.drive_work_j;
            if e_drv > 0.0 && e_imp > 0.0 {
                self.last_impulse_energy_j = e_imp;
                self.escapement_efficiency = (e_imp / e_drv).clamp(0.0, 1.0);
            }
            if let Some((a, b)) = self.pin_arc.take() {
                self.measured_lift_deg = (b - a).abs() / DEG;
            }
        }
    }

    fn track_swing(&mut self) {
        let th = self.esc.q[0] - self.loads.theta_eq;
        let w = self.esc.v[0];
        let dir = w.signum();
        self.swing_peak = self.swing_peak.max(th.abs());
        if dir != 0.0 && self.swing_dir != 0.0 && dir != self.swing_dir {
            self.last_peaks = [self.last_peaks[1], self.swing_peak];
            self.last_peak_s = self.lab_time_s;
            self.swing_peak = th.abs();
        }
        if dir != 0.0 {
            self.swing_dir = dir;
        }
    }

    fn update_measurements(&mut self) {
        // Amplitude from the last two turning points; a balance that has not
        // turned for half a second is reported by how far it has strayed.
        if self.lab_time_s - self.last_peak_s < 0.5 {
            self.balance_amplitude_deg = 0.5 * (self.last_peaks[0] + self.last_peaks[1]) / DEG;
        } else {
            let th = (self.esc.q[0] - self.loads.theta_eq).abs();
            self.balance_amplitude_deg = self.swing_peak.max(th) / DEG;
            self.swing_peak = th;
            self.last_peaks = [th, th];
        }
        if self.lab_time_s - self.last_lock_s > 1.0 {
            self.escapement_efficiency = 0.0;
        }
        if let Some((rate, sigma)) = self.fit_rate(false) {
            self.measured_rate_s_per_day = rate;
            self.rate_sigma_s_per_day = sigma;
            self.rate_valid = self.beats.len() >= 20;
        }
        // Beat error: half the difference of the mean tick and tock intervals.
        let (mut s0, mut n0, mut s1, mut n1) = (0.0, 0.0, 0.0, 0.0);
        for w in self.beats.iter().collect::<Vec<_>>().windows(2) {
            let d = w[1].lab_time_s - w[0].lab_time_s;
            if d > 0.5 * NOMINAL_BEAT_PERIOD_S && d < 2.0 * NOMINAL_BEAT_PERIOD_S {
                if w[0].side == 0 {
                    s0 += d;
                    n0 += 1.0;
                } else {
                    s1 += d;
                    n1 += 1.0;
                }
            }
        }
        if n0 > 0.0 && n1 > 0.0 {
            self.measured_beat_error_ms = 0.5 * (s0 / n0 - s1 / n1).abs() * 1.0e3;
        }
        // COSC compliance e-process, one observation per 25 beats.
        if self.rate_valid && self.total_beats as u32 >= self.beats_since_cosc + 25 {
            self.beats_since_cosc = self.total_beats as u32;
            let out = !(-4.0..=6.0).contains(&self.measured_rate_s_per_day);
            self.cosc.observe(if out { 1.0 } else { 0.0 });
            self.cosc_e_value = self.cosc.e_value();
        }
    }

    fn update_kinematics(&mut self) {
        let psi = self.esc.q[2] - self.psi_origin;
        // Escape wheel turns 10 per cage turn; one cage turn per minute.
        let elapsed_watch = psi / TAU * 60.0 / train::CAGE_TO_ESCAPE;
        let t = self.lab_origin_s + elapsed_watch;
        self.watch_displayed_time_s = t;
        self.tourbillon_cage_angle_rad = t / 60.0 * TAU;
        self.third_wheel_angle_rad = -self.tourbillon_cage_angle_rad / train::THIRD_TO_CAGE;
        self.center_wheel_angle_rad = t / 3600.0 * TAU;
        self.barrel_wheel_angle_rad = -self.center_wheel_angle_rad / train::BARREL_TO_CENTER;
        self.moonphase_angle_rad = (0.35 + t / (29.530_589 * SECONDS_PER_DAY)) * TAU;
        self.mainspring_coil_radii_mm = self.mainspring.coil_radii_mm();
        self.click_lever_angle_rad *= 0.6;
    }

    fn diagnose(&mut self, env: &MovementEnvironment) {
        let running = self.lab_time_s - self.last_lock_s < 4.0 * NOMINAL_BEAT_PERIOD_S;
        if running {
            let in_band = !self.rate_valid || (-4.0..=6.0).contains(&self.measured_rate_s_per_day);
            self.status = if in_band && self.balance_amplitude_deg >= 180.0 {
                WatchOperatingStatus::RunningChronometer
            } else {
                WatchOperatingStatus::RunningPerturbed
            };
            self.stop_cause = None;
            return;
        }
        // Stopped: name the resistance that dominates. Train-side blocks are
        // compared with the escape-wheel drive; balance-side losses with the
        // energy the escapement could deliver per cycle at 200°.
        let l = &self.losses;
        let drive = self.escape_drive_nm.max(1e-15);
        let deliverable = 2.0 * drive * (TAU / 30.0) * 0.3;
        let k = self.loads.k_hairspring.max(1e-12);
        let hold = 26.0 * DEG * k;
        let mut best = (WatchOperatingStatus::StoppedMainspringUnwound, 0.0);
        let mut consider = |s: WatchOperatingStatus, score: f64| {
            if score > best.1 {
                best = (s, score);
            }
        };
        consider(
            WatchOperatingStatus::LockedInSolidifiedGold,
            (l.train_solid_gold_nm / drive).max(l.balance_solid_gold_nm / hold) * 10.0,
        );
        consider(
            WatchOperatingStatus::JammedByMoltenGold,
            (l.train_liquid_gold_nm / drive).max(l.liquid_gold / deliverable),
        );
        let oil_glass = l.pivot_oil / deliverable;
        consider(
            WatchOperatingStatus::SeizedByThermalOrCryoClearance,
            (l.train_ice_nm / drive).max(l.balance_ice_nm / hold).max(oil_glass).max(if self.min_clearance_um <= 0.0 {
                2.0
            } else {
                0.0
            }),
        );
        let r_lock = self.esc.design.entry_stone[0].len().max(1e-4);
        let draw_torque = drive / self.esc.design.r_tip * self.esc.design.draw.sin() * r_lock;
        // A lever pulled across by a magnet (and left overbanked), or a
        // magnetised hairspring whose coils cling together.
        consider(
            WatchOperatingStatus::PinnedByMagneticField,
            (env.lever_magnetic_torque_nm.abs() / (draw_torque + 1e-15)).max(env.hairspring_magnetisation.abs() / 0.1),
        );
        consider(WatchOperatingStatus::FloodedStalledByOil, l.oil_flood / deliverable);
        let reserve = self.power_reserve_h();
        let unwound = if reserve <= 0.0 { 50.0 } else { (0.5 / reserve).min(5.0) };
        consider(WatchOperatingStatus::StoppedMainspringUnwound, unwound);
        // A clear present cause wins; otherwise the cause it stopped from
        // stands (removing the magnet does not un-jam an overbanked lever).
        let status = match self.stop_cause {
            Some((cause, _)) if best.1 < 1.0 => cause,
            _ => best.0,
        };
        if best.1 >= 1.0 || self.stop_cause.is_none() {
            self.stop_cause = Some((status, best.1));
        }
        self.status = status;
    }

    /// Lock drops since the last call: `[lab_time_s, side, energy_nJ]` triples.
    pub fn drain_timegrapher(&mut self) -> Vec<f64> {
        core::mem::take(&mut self.drain)
    }

    /// Power reserve [h]: barrel turns left above the torque at which the
    /// amplitude would fall below ~150°, at one barrel turn per 7.8 h.
    #[must_use]
    pub fn power_reserve_h(&self) -> f64 {
        let stop = self.mainspring.slip_torque_nm * 0.35;
        let ms = &self.mainspring;
        let delivered = |n: f64| ms.torque_at(n, T_REF_K) * (1.0 - ms.coil_friction);
        // Below the stopping torque at the present wind (a spring that has lost
        // its set may never reach it again) there is no reserve at all.
        if delivered(ms.turns) < stop {
            return 0.0;
        }
        // The output torque rises monotonically with wind: bisect for the
        // turns at which it falls to the stopping torque.
        let (mut lo, mut hi) = (0.0, ms.turns);
        for _ in 0..40 {
            let mid = 0.5 * (lo + hi);
            if delivered(mid) >= stop {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        (ms.turns - hi).max(0.0) * train::BARREL_TO_CENTER
    }

    /// Wind state as a fraction of the slip turns.
    #[must_use]
    pub fn wind_fraction(&self) -> f64 {
        (self.mainspring.turns / self.mainspring.slip_turns(T_REF_K).max(1e-9)).clamp(0.0, 1.0)
    }

    /// 16-value escapement state for the front end (see `escapement_state`).
    #[must_use]
    pub fn escapement_state(&self) -> [f64; 16] {
        let l = &self.losses;
        let total = l.total().max(1e-30);
        [
            self.esc.q[0],
            self.esc.v[0],
            self.esc.q[1],
            self.esc.v[1],
            self.esc.q[2],
            self.esc.v[2],
            f64::from(self.esc.active_mask),
            self.last_impulse_energy_j * 1.0e6,
            self.escapement_efficiency,
            self.measured_lift_deg,
            self.balance_q,
            self.escape_drive_nm * 1.0e6,
            self.mainspring_torque_nm * 1.0e3,
            self.power_reserve_h(),
            self.balance_pivot_friction_nm * 1.0e9,
            l.air / total,
        ]
    }

    #[must_use]
    pub fn rate_interval(&self) -> (f64, f64) {
        let r = self.measured_rate_s_per_day;
        let s = 2.0 * self.rate_sigma_s_per_day;
        (r - s, r + s)
    }
}

impl Default for Horology {
    fn default() -> Self {
        Self::new()
    }
}

/// Shear strength of frost-bonded ice [Pa] (hoar frost to clear ice spans
/// ~0.1–1 MPa; adhesion of ice on steel ~0.5 MPa).
const ICE_SHEAR_PA: f64 = 0.5e6;

/// Shear yield stress of 18k gold [Pa] at `t_k` with liquid fraction `fl`:
/// ~150 MPa annealed at room temperature (tensile ~260 MPa / √3), falling with
/// homologous temperature, and for a mush only above its coherency point
/// (solid fraction > 0.5).
#[must_use]
pub fn gold_shear_yield_pa(t_k: f64, fl: f64) -> f64 {
    let gold = &library().gold;
    let th = ((t_k - T_REF_K) / (gold.t_solidus_k - T_REF_K)).clamp(0.0, 1.0);
    let solid = 150.0e6 * (1.0 - th).powf(1.5) + 2.0e6 * (1.0 - th);
    let coherent = ((0.5 - fl) / 0.5).clamp(0.0, 1.0);
    solid * coherent.powi(2)
}

#[derive(Clone, Copy)]
struct Rotor {
    radius: f64,
    z0: f64,
    thick: f64,
    teeth: f64,
    ratio: f64,
    part: PartField,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn train_ratio_and_efficiency() {
        assert!((train::RATIO - 4680.0).abs() < 1e-9);
        let eta = train_efficiency(train::MESH_MU);
        assert!(eta > 0.68 && eta < 0.76, "train efficiency {eta}");
    }

    #[test]
    fn balance_and_mainspring_geometry() {
        let b = BalanceDesign::pocket_tourbillon();
        let (m, i) = b.mass_and_inertia();
        assert!(m > 0.08e-3 && m < 0.2e-3, "balance mass {m}");
        assert!(i > 2.0e-9 && i < 4.0e-9, "balance inertia {i}");
        let ms = Mainspring::pocket_tourbillon();
        let n = ms.geometric_turns();
        assert!(n > 8.0 && n < 16.0, "usable turns {n}");
        let slip = ms.slip_turns(T_REF_K);
        assert!(slip > 4.0 && slip <= n, "slip turns {slip}");
        let mut ms2 = ms;
        ms2.turns = 0.0;
        let r0 = ms2.coil_radii_mm();
        ms2.turns = n;
        let r1 = ms2.coil_radii_mm();
        // Let down: coils on the wall; wound: coils on the arbor.
        assert!(r0[0] > 5.0 && r1[MAINSPRING_COILS - 1] < 5.5, "{r0:?} {r1:?}");
    }

    #[test]
    fn power_reserve_falls_with_wind_and_with_lost_set() {
        let mut h = Horology::new();
        let fresh = h.power_reserve_h();
        assert!(fresh > 30.0 && fresh < 60.0, "fresh reserve {fresh} h");
        // Less wind, less reserve.
        let full = h.mainspring.turns;
        h.mainspring.turns = full * 0.6;
        let partial = h.power_reserve_h();
        assert!(partial < fresh && partial >= 0.0, "partial {partial} h vs fresh {fresh} h");
        h.mainspring.turns = full;
        // An over-tempered spring that lost its set never holds more than a fresh one, and
        // one too weak to reach the stopping torque holds none.
        let mut last = fresh;
        for set in [0.2, 0.4, 0.6, 0.8, 0.9] {
            h.mainspring.set_fraction = set;
            let r = h.power_reserve_h();
            assert!(r <= last + 1e-9, "set {set}: reserve {r} h rose above {last} h");
            last = r;
        }
        assert!(last == 0.0, "a spring at 90 % set still claims {last} h");
    }

    #[test]
    fn regulated_watch_keeps_time_and_measures_itself() {
        let mut h = Horology::new();
        let env = MovementEnvironment::default();
        let mut ev = Vec::new();
        for _ in 0..240 {
            h.step(0.05, &env, &mut ev);
        }
        assert!(h.total_beats > 50, "beats {}", h.total_beats);
        assert!(h.measured_rate_s_per_day.abs() < 3.0, "rate {}", h.measured_rate_s_per_day);
        assert!(
            h.balance_amplitude_deg > 220.0 && h.balance_amplitude_deg < 320.0,
            "amplitude {}",
            h.balance_amplitude_deg
        );
        assert_eq!(h.status, WatchOperatingStatus::RunningChronometer);
        assert!(h.balance_q > 80.0 && h.balance_q < 600.0, "Q {}", h.balance_q);
        assert!(h.escapement_efficiency > 0.1, "efficiency {}", h.escapement_efficiency);
        assert!(!h.drain_timegrapher().is_empty());
    }

    #[test]
    fn uncompensated_steel_spring_loses_about_eleven_seconds_per_kelvin() {
        let mut h = Horology::new();
        h.set_compensated(false);
        let mut env = MovementEnvironment::default();
        env.balance.temperature_k = T_REF_K + 10.0;
        let _ = h.assemble_loads(&env);
        let s = h.rate_sensitivity_s_day_per_k;
        assert!(s < -9.0 && s > -13.0, "steel sensitivity {s}");
        h.set_compensated(true);
        let _ = h.assemble_loads(&env);
        assert!(h.rate_sensitivity_s_day_per_k.abs() < 0.5, "nivarox {}", h.rate_sensitivity_s_day_per_k);
    }

    #[test]
    #[ignore = "diagnostic dump"]
    fn dump_running_watch() {
        let mut h = Horology::new();
        let env = MovementEnvironment::default();
        let mut ev = Vec::new();
        let t0 = std::time::Instant::now();
        for k in 0..400 {
            h.step(0.05, &env, &mut ev);
            if k % 40 == 0 {
                println!(
                    "t={:.1} beats {} rate {:.2}±{:.2} amp {:.1} be {:.3} eff {:.2} lift {:.1} Q {:.0} drive {:.3e} status {:?}",
                    h.lab_time_s - h.lab_origin_s,
                    h.total_beats,
                    h.measured_rate_s_per_day,
                    h.rate_sigma_s_per_day,
                    h.balance_amplitude_deg,
                    h.measured_beat_error_ms,
                    h.escapement_efficiency,
                    h.measured_lift_deg,
                    h.balance_q,
                    h.escape_drive_nm,
                    h.status
                );
            }
        }
        println!("losses {:?}", h.losses);
        println!("20 s of movement in {:?}", t0.elapsed());
    }

    #[test]
    fn factory_regulation_constants_are_current() {
        let balance = BalanceDesign::pocket_tourbillon();
        for (i, compensated) in [true, false].into_iter().enumerate() {
            let free = Horology::free_length(&balance, Horology::alloy_for(compensated));
            let factor = Horology::regulate(compensated) / free;
            // 0.5 s/day of rate is a length error of 2 × 0.5 / 86400.
            let tol = 2.0 * 0.5 / SECONDS_PER_DAY;
            println!("REGULATION_FACTOR[{i}] = {factor:.9}");
            assert!(
                (factor - REGULATION_FACTOR[i]).abs() < tol,
                "re-bake REGULATION_FACTOR[{i}] = {factor:.9} (was {})",
                REGULATION_FACTOR[i]
            );
        }
    }
}
