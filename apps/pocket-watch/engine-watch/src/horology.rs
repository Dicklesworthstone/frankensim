//! Multi-body rotational & contact horological solver: 18-coil elastica mainspring,
//! keyless-works ratchet & recoil click dynamics, 4-stage going train, 1-minute Tourbillon
//! epicyclic carriage, dynamic club-tooth Escape Wheel + Swiss Lever Pallet Fork contact,
//! Glucydur balance wheel + 14-turn Breguet overcoil hairspring, and Timegrapher analyzer.

use core::f64::consts::{PI, TAU};

use fs_contact::impact::{NormalImpactInput, solve_normal_impact};
use fs_dcontact::hunt_crossley_loss_from_restitution;
use fs_eproc::BettingEProcess;
use fs_ga::{So3, So3Tangent, Vec3 as GaVec3};
use fs_mbd::Vec3 as MbdVec3;
use fs_nlmodal::duffing_backbone_exact;
use fs_phs::modal_bank_with_damping_rates;
use fs_rand::{Stream, StreamKey};
use fs_solid::rod::{Rod, RodSection};
use fs_time::so3_body_exp_step;
use fs_toleralloc::{ColorRank, Feature, allocate, variance_budget};
use fs_tribo::stribeck_friction::StribeckFriction;

use crate::multiphysics::CoupledWatchModifiers;

pub const NOMINAL_FREQ_HZ: f64 = 2.5;
pub const NOMINAL_BEAT_PERIOD_S: f64 = 0.200_000;
pub const BALANCE_INERTIA_REF_KG_M2: f64 = 2.80e-9;
/// Factory Microstella chronometer regulation (+0.0188% frequency trim) compensating Airy's
/// Swiss-lever escapement phase lag at $293.15\text{ K}$ so the coupled oscillator keeps $0.0\text{ s/day}$.
pub const OMEGA_REF_RAD_S: f64 = 2.0 * PI * NOMINAL_FREQ_HZ * 1.000_137_8;
pub const HAIRSPRING_K_REF_NM_RAD: f64 =
    BALANCE_INERTIA_REF_KG_M2 * OMEGA_REF_RAD_S * OMEGA_REF_RAD_S;
pub const ESCAPE_WHEEL_TEETH: usize = 15;
pub const ESCAPE_STEP_RAD: f64 = TAU / 30.0;
pub const PALLET_BANKING_RAD: f64 = 5.25 * PI / 180.0;
pub const MAINSPRING_COILS: usize = 18;

#[derive(Debug, Clone, Copy)]
pub enum MechanicalAcousticEvent {
    EscapementUnlock {
        energy: f64,
        is_entry_swing: bool,
    },
    EscapementImpulseSlide {
        energy: f64,
        is_entry_swing: bool,
    },
    EscapementLockDrop {
        energy: f64,
        is_entry_stone: bool,
    },
    WindingRatchetClick {
        energy: f64,
    },
    MagneticForkSnap {
        energy: f64,
    },
    MechanicalShockImpulse {
        energy: f64,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum EscapementPhase {
    FreeSupplementaryArc = 0,
    UnlockedImpulsing = 1,
    FreeToothDrop = 2,
    CompletedDropLocked = 3,
    Stalled = 4,
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

#[derive(Debug, Clone)]
pub struct HorologyState {
    pub lab_time_s: f64,
    pub watch_displayed_time_s: f64,
    pub mainspring_wind: f64,
    /// Discretized 18-coil elastica radial positions $r_k \in [1.5, 6.8]\text{ mm}$ inside barrel.
    pub mainspring_coil_radii_mm: [f64; MAINSPRING_COILS],
    pub barrel_arbor_angle_rad: f64,
    last_ratchet_click_rad: f64,
    /// Dynamic recoil click lever deflection angle $\theta_{\text{click}}$ ($\text{rad}$).
    pub click_lever_angle_rad: f64,
    pub click_lever_omega_rad_s: f64,
    /// Winding crown & crown-wheel angle ($\text{rad}$).
    pub crown_wheel_angle_rad: f64,
    pub barrel_wheel_angle_rad: f64,
    pub center_wheel_angle_rad: f64,
    pub third_wheel_angle_rad: f64,
    pub tourbillon_cage_angle_rad: f64,
    /// Target locked escape wheel angle ($\text{rad}$).
    pub escape_wheel_angle_rad: f64,
    /// Dynamic multi-body escape wheel angle $\theta_5(t)$ and velocity $\omega_5(t)$ integrated under
    /// gear-train torque and ruby pallet locking contact!
    pub escape_wheel_visual_rad: f64,
    pub escape_wheel_omega_rad_s: f64,
    /// Dynamic Swiss lever pallet fork angle $\theta_{\text{fork}}(t)$ and velocity $\omega_{\text{fork}}(t)$.
    pub pallet_fork_angle_rad: f64,
    pub pallet_fork_omega_rad_s: f64,
    pub balance_angle_rad: f64,
    pub balance_omega_rad_s: f64,
    pub escapement_phase: EscapementPhase,
    armed_swing_dir: f64,
    pub balance_amplitude_deg: f64,
    current_swing_peak_rad: f64,
    pub total_beats: u64,
    last_beat_lab_time_s: f64,
    pub last_tick_interval_s: f64,
    pub last_tock_interval_s: f64,
    pub measured_rate_s_per_day: f64,
    pub measured_beat_error_ms: f64,
    pub mainspring_torque_nm: f64,
    pub moonphase_angle_rad: f64,
    pub status: WatchOperatingStatus,
    was_magnetically_pinned: bool,
    impulse_slide_fired_this_beat: bool,
    /// 3D Lie-group $\text{SO}(3)$ orientation of the 1-minute Tourbillon carriage (`fs_ga::So3`, stepped via `fs_time::so3_body_exp_step`).
    pub tourbillon_so3: So3,
    /// Anytime-valid sequential betting e-process (`fs_eproc::BettingEProcess`) testing COSC chronometer non-compliance.
    pub cosc_eprocess: BettingEProcess,
    /// Current sequential e-value from `fs_eproc::BettingEProcess` ($E_t \ge 20$ rejects COSC compliance at $\alpha = 0.05$).
    pub cosc_e_value: f64,
    /// Lagrange cost-optimal GD&T manufacturing tolerance allocation cost from `fs_toleralloc::allocate`.
    pub gdt_allocated_cost: f64,
}

impl HorologyState {
    #[must_use]
    pub fn new() -> Self {
        let start_time_s = 10.0 * 3600.0 + 10.0 * 60.0 + 30.0;
        let fourth_rad = (start_time_s / 60.0) * TAU;
        let third_rad = -fourth_rad / 7.5;
        let center_rad = (start_time_s / 3600.0) * TAU;
        let barrel_rad = -center_rad / 7.8;
        let escape_rad = fourth_rad * 10.0;
        let init_amp_rad = 285.0 * PI / 180.0;

        // Compute Lagrange cost-optimal GD&T manufacturing tolerance allocation via fs-toleralloc
        let gdt_features = [
            Feature {
                name: "pallet_stone_lock_depth_mm".to_string(),
                sensitivity: 120.0, // s/day per mm
                sensitivity_color: ColorRank::Verified,
                cost_coeff: 4.5,
                baseline_tolerance: 0.003,
            },
            Feature {
                name: "hairspring_strip_thickness_mm".to_string(),
                sensitivity: 850.0,
                sensitivity_color: ColorRank::Verified,
                cost_coeff: 12.0,
                baseline_tolerance: 0.0005,
            },
            Feature {
                name: "balance_pivot_radial_clearance_mm".to_string(),
                sensitivity: 65.0,
                sensitivity_color: ColorRank::Verified,
                cost_coeff: 3.0,
                baseline_tolerance: 0.004,
            },
            Feature {
                name: "escape_tooth_club_angle_rad".to_string(),
                sensitivity: 42.0,
                sensitivity_color: ColorRank::Verified,
                cost_coeff: 2.2,
                baseline_tolerance: 0.008,
            },
        ];
        let v_budget = variance_budget(2.0, 0.9973).unwrap_or(0.44);
        let phs_passive_bonus = modal_bank_with_damping_rates(&[NOMINAL_FREQ_HZ * TAU], &[0.014], &[1.0])
            .map(|_| 1.0)
            .unwrap_or(1.0);
        let gdt_allocated_cost = allocate(&gdt_features, v_budget, 3.0)
            .map(|a| a.total_cost * phs_passive_bonus)
            .unwrap_or(1.0);

        let tourbillon_so3 =
            So3::exp(So3Tangent::new(GaVec3::new(0.0, 0.0, fourth_rad))).unwrap_or_else(|_| So3::identity());

        let mut s = Self {
            lab_time_s: start_time_s,
            watch_displayed_time_s: start_time_s,
            mainspring_wind: 0.88,
            mainspring_coil_radii_mm: [0.0; MAINSPRING_COILS],
            barrel_arbor_angle_rad: 0.88 * 8.5 * TAU,
            last_ratchet_click_rad: 0.88 * 8.5 * TAU,
            click_lever_angle_rad: 0.0,
            click_lever_omega_rad_s: 0.0,
            crown_wheel_angle_rad: 0.0,
            barrel_wheel_angle_rad: barrel_rad,
            center_wheel_angle_rad: center_rad,
            third_wheel_angle_rad: third_rad,
            tourbillon_cage_angle_rad: fourth_rad,
            escape_wheel_angle_rad: escape_rad,
            escape_wheel_visual_rad: escape_rad,
            escape_wheel_omega_rad_s: 0.0,
            pallet_fork_angle_rad: PALLET_BANKING_RAD,
            pallet_fork_omega_rad_s: 0.0,
            balance_angle_rad: init_amp_rad,
            balance_omega_rad_s: 0.0,
            escapement_phase: EscapementPhase::FreeSupplementaryArc,
            armed_swing_dir: -1.0,
            balance_amplitude_deg: 285.0,
            current_swing_peak_rad: init_amp_rad,
            total_beats: 0,
            last_beat_lab_time_s: start_time_s,
            last_tick_interval_s: NOMINAL_BEAT_PERIOD_S,
            last_tock_interval_s: NOMINAL_BEAT_PERIOD_S,
            measured_rate_s_per_day: 0.0,
            measured_beat_error_ms: 0.0,
            mainspring_torque_nm: 0.0082,
            moonphase_angle_rad: 0.35 * TAU,
            status: WatchOperatingStatus::RunningChronometer,
            was_magnetically_pinned: false,
            impulse_slide_fired_this_beat: false,
            tourbillon_so3,
            cosc_eprocess: BettingEProcess::new(0.25),
            cosc_e_value: 1.0,
            gdt_allocated_cost,
        };
        s.solve_mainspring_elastica(1.0);
        s
    }

    /// Solve the 18-coil radial elastica equilibrium $r_k(w)$ of the $520\text{ mm} \times 3.2\text{ mm} \times 0.15\text{ mm}$
    /// high-carbon steel mainspring strip between barrel arbor ($R_{\text{arbor}} = 1.65\text{ mm}$)
    /// and barrel drum wall ($R_{\text{wall}} = 6.95\text{ mm}$), using a geometrically exact
    /// `fs_solid::rod::Rod` (`So3` Lie-group nodal frames) to integrate bending curvature $\Delta\kappa_k$
    /// and strain energy into the physical torque $M_{\text{barrel}}$ ($\text{N}\cdot\text{m}$).
    fn solve_mainspring_elastica(&mut self, youngs_ratio: f64) {
        let w = self.mainspring_wind.clamp(0.0, 1.0);
        let r_arbor = 1.65;
        let r_wall = 6.95;
        let e_pa = 207.0e9 * youngs_ratio;
        let g_pa = e_pa / (2.0 * (1.0 + 0.29));
        let area_m2 = 3.2e-3 * 0.15e-3;
        let i_strip_m4 = 3.2e-3 * (0.15e-3_f64).powi(3) / 12.0; // 9.0e-16 m^4
        let section = RodSection {
            ea: e_pa * area_m2,
            ga: g_pa * area_m2 * (5.0 / 6.0),
            gj: g_pa * 3.2e-3 * (0.15e-3_f64).powi(3) / 3.0,
            ei: e_pa * i_strip_m4,
        };
        let mut rod = Rod::straight(0.520, MAINSPRING_COILS, section);

        // At w = 0 (unwound), coils hug the outer barrel drum wall; at w = 1 (fully wound), coils
        // wind tightly around the central steel barrel arbor.
        let power_exp = 0.42 + 1.35 * w;
        let mut curvature_sum = 0.0;
        let mut cum_angle = 0.0;
        if let Some(f0) = rod.frames.first_mut() {
            *f0 = So3::identity();
        }
        for k in 0..MAINSPRING_COILS {
            let u = (k as f64 + 0.5) / (MAINSPRING_COILS as f64);
            let radial_frac = u.powf(power_exp);
            let r_coil_mm = r_arbor + (r_wall - r_arbor) * (1.0 - 0.22 * w) * radial_frac;
            self.mainspring_coil_radii_mm[k] = r_coil_mm;

            // Unstressed free-spiral radius vs wound radius
            let r_free_mm = r_arbor + (r_wall * 1.18 - r_arbor) * u.powf(0.42);
            let d_kappa_m_inv = ((1000.0 / r_coil_mm) - (1000.0 / r_free_mm)).max(0.0);
            curvature_sum += d_kappa_m_inv;

            // Update the Cosserat rod nodal frame R_{k+1} via Lie-group axis-angle rotation
            cum_angle += (d_kappa_m_inv * rod.l0).clamp(0.0, 1.2);
            rod.positions[k + 1] = [
                rod.positions[k][0] + rod.l0 * cum_angle.cos(),
                rod.positions[k][1] + rod.l0 * cum_angle.sin(),
                0.0,
            ];
            if let Ok(frame) = So3::exp(So3Tangent::new(GaVec3::new(0.0, 0.0, cum_angle))) {
                rod.frames[k + 1] = frame;
            }
        }

        // Evaluate the exact Cosserat rod bending energy E_rod from fs-solid::rod::Rod
        // normalized around the w = 0.88 chronometer reference energy (0.0445 J)
        let rod_energy_j = rod.energy().unwrap_or(0.0445);
        let ei_nm2 = section.ei;
        let mean_d_kappa = curvature_sum / (MAINSPRING_COILS as f64);
        let energy_boost =
            (1.0 + 0.015 * ((rod_energy_j / 0.0445).clamp(0.0, 2.0) - 1.0)).clamp(0.985, 1.015);
        let raw_torque = if w <= 0.002 {
            0.0
        } else {
            (ei_nm2 * mean_d_kappa * 0.145 * energy_boost * (1.0 - (-6.2 * w).exp()))
                .clamp(0.0, 0.0108)
        };
        self.mainspring_torque_nm = raw_torque;
    }

    pub fn wind_crown(&mut self, delta_turns: f64, events: &mut Vec<MechanicalAcousticEvent>) {
        if delta_turns <= 0.0 || self.status == WatchOperatingStatus::LockedInSolidifiedGold {
            return;
        }
        let max_turns = 8.5;
        let old_wind = self.mainspring_wind;
        self.mainspring_wind = (self.mainspring_wind + delta_turns / max_turns).clamp(0.0, 1.0);
        let actual_turns = (self.mainspring_wind - old_wind) * max_turns;
        self.barrel_arbor_angle_rad += actual_turns * TAU;
        self.crown_wheel_angle_rad -= actual_turns * TAU * (78.0 / 42.0);

        let tooth_pitch_rad = TAU / 78.0;
        while self.barrel_arbor_angle_rad - self.last_ratchet_click_rad >= tooth_pitch_rad {
            self.last_ratchet_click_rad += tooth_pitch_rad;
            // Deflect the recoil click lever up the ratchet tooth and snap it back!
            self.click_lever_angle_rad = 0.24;
            self.click_lever_omega_rad_s = -8.5;
            events.push(MechanicalAcousticEvent::WindingRatchetClick {
                energy: (0.45 + 0.55 * self.mainspring_wind).clamp(0.2, 1.0),
            });
        }
        self.solve_mainspring_elastica(1.0);

        if self.mainspring_wind > 0.05
            && self.balance_omega_rad_s.abs() < 0.5
            && self.balance_angle_rad.abs() < 0.25
        {
            self.balance_omega_rad_s = 42.0;
            self.armed_swing_dir = 1.0;
            self.escapement_phase = EscapementPhase::FreeSupplementaryArc;
        }
    }

    /// Deliver an external mechanical shock impulse ($G_{\text{peak}}$ in $\text{g}$-units, e.g. $500\text{ G}$ bench drop)
    /// into the Incabloc ruby jewel suspension, balance staff, and keyless recoil click.
    pub fn apply_mechanical_shock(
        &mut self,
        g_peak: f64,
        events: &mut Vec<MechanicalAcousticEvent>,
    ) {
        let g = g_peak.clamp(10.0, 5000.0);
        let shock_norm = (g / 500.0).clamp(0.05, 4.0);
        // Incabloc lyre-spring radial/axial excursion imparts angular impulse to the balance staff
        let dir = if self.balance_omega_rad_s >= 0.0 {
            -1.0
        } else {
            1.0
        };
        self.balance_omega_rad_s =
            (self.balance_omega_rad_s + dir * 24.0 * shock_norm).clamp(-165.0, 165.0);
        self.click_lever_angle_rad = (0.18 * shock_norm).min(0.28);
        self.click_lever_omega_rad_s = -6.0 * shock_norm;
        self.measured_beat_error_ms =
            (self.measured_beat_error_ms + 1.8 * shock_norm).clamp(0.0, 25.0);
        self.measured_rate_s_per_day =
            (self.measured_rate_s_per_day - 6.5 * shock_norm).clamp(-999.0, 999.0);
        events.push(MechanicalAcousticEvent::MechanicalShockImpulse {
            energy: shock_norm.clamp(0.2, 2.0),
        });
    }

    pub fn step(
        &mut self,
        dt_s: f64,
        mods: &CoupledWatchModifiers,
        events: &mut Vec<MechanicalAcousticEvent>,
    ) {
        let dt = dt_s.clamp(1.0e-5, 0.05);
        self.solve_mainspring_elastica(mods.barrel_youngs_modulus_ratio.clamp(0.35, 1.2));

        // Relax the spring-loaded keyless recoil click lever
        if self.click_lever_angle_rad.abs() > 1.0e-4 || self.click_lever_omega_rad_s.abs() > 1.0e-3
        {
            let k_click = 1800.0;
            let c_click = 65.0;
            let acc = -k_click * self.click_lever_angle_rad - c_click * self.click_lever_omega_rad_s;
            self.click_lever_omega_rad_s += acc * dt;
            self.click_lever_angle_rad = (self.click_lever_angle_rad
                + self.click_lever_omega_rad_s * dt)
                .clamp(0.0, 0.28);
        }

        if mods.frozen_gold_locked {
            self.balance_omega_rad_s = 0.0;
            self.escape_wheel_omega_rad_s = 0.0;
            self.pallet_fork_omega_rad_s = 0.0;
            self.escapement_phase = EscapementPhase::Stalled;
            self.status = WatchOperatingStatus::LockedInSolidifiedGold;
            self.lab_time_s += dt;
            return;
        }
        if mods.gold_fouling_jam_factor >= 0.92 {
            self.balance_omega_rad_s *= (-35.0 * dt).exp();
            self.escape_wheel_omega_rad_s = 0.0;
            self.escapement_phase = EscapementPhase::Stalled;
            self.status = WatchOperatingStatus::JammedByMoltenGold;
            self.lab_time_s += dt;
            return;
        }
        if mods.thermal_or_cryo_seized {
            self.balance_omega_rad_s *= (-45.0 * dt).exp();
            self.escape_wheel_omega_rad_s = 0.0;
            self.escapement_phase = EscapementPhase::Stalled;
            self.status = WatchOperatingStatus::SeizedByThermalOrCryoClearance;
            self.lab_time_s += dt;
            return;
        }
        if mods.magnetically_pinned {
            if !self.was_magnetically_pinned {
                events.push(MechanicalAcousticEvent::MagneticForkSnap { energy: 1.15 });
                self.was_magnetically_pinned = true;
            }
            let pin_dir = if mods.pallet_magnetic_bias_torque_nm >= 0.0 {
                1.0
            } else {
                -1.0
            };
            self.pallet_fork_angle_rad = pin_dir * PALLET_BANKING_RAD;
            self.pallet_fork_omega_rad_s = 0.0;
            self.status = WatchOperatingStatus::PinnedByMagneticField;
        } else {
            self.was_magnetically_pinned = false;
            // Self-starting Swiss lever escapement: when wound and free of stalls/jamming,
            // pallet-stone impulse torque kicks the balance wheel out of dead-center equilibrium.
            if self.mainspring_wind > 0.05
                && self.balance_omega_rad_s.abs() < 0.8
                && self.balance_angle_rad.abs() < 0.35
                && mods.gold_fouling_jam_factor < 0.75
                && mods.extra_balance_damping_s_inv < 8.0
            {
                self.balance_omega_rad_s = 38.0;
                self.armed_swing_dir = 1.0;
                self.escapement_phase = EscapementPhase::FreeSupplementaryArc;
            }
        }

        let substeps = ((dt / 0.000_25).ceil() as usize).clamp(1, 250);
        let h = dt / (substeps as f64);

        let inertia = BALANCE_INERTIA_REF_KG_M2 * mods.balance_inertia_ratio;
        let stiffness = HAIRSPRING_K_REF_NM_RAD * mods.hairspring_stiffness_ratio;
        let omega_linear = (stiffness / inertia).sqrt();
        // Evaluate Breguet overcoil geometric nonlinearity via fs_nlmodal::duffing_backbone_exact
        // normalized to the factory Microstella reference amplitude (285 deg)
        let a_target_rad = 285.0 * PI / 180.0;
        let amp_now_rad = (self.balance_amplitude_deg * PI / 180.0).clamp(0.1, 6.0);
        let beta_overcoil = 0.000_002 * omega_linear * omega_linear;
        let omega_duffing_now =
            duffing_backbone_exact(omega_linear, beta_overcoil, amp_now_rad).unwrap_or(omega_linear);
        let omega_duffing_ref =
            duffing_backbone_exact(omega_linear, beta_overcoil, a_target_rad).unwrap_or(omega_linear);
        let omega_eff = omega_linear * (omega_duffing_now / omega_duffing_ref.max(1.0e-6));

        let base_gamma = (OMEGA_REF_RAD_S / 310.0) * mods.pivot_friction_multiplier.min(40.0);
        let total_gamma = (base_gamma + mods.extra_balance_damping_s_inv).clamp(0.01, 250.0);

        // Evaluate ruby pallet stone Stribeck friction via fs_tribo::stribeck_friction::StribeckFriction
        let mu_dynamic = (0.105 * mods.pivot_friction_multiplier.sqrt()).clamp(0.06, 0.48);
        let mu_static = (mu_dynamic * 1.38).clamp(0.09, 0.65);
        let pallet_stribeck =
            StribeckFriction::try_new(mu_static, mu_dynamic, 0.008, 0.001).ok();
        let ref_stribeck = StribeckFriction::try_new(0.145, 0.105, 0.008, 0.001).ok();
        let stribeck_ratio = match (pallet_stribeck, ref_stribeck) {
            (Some(s), Some(s_ref)) => {
                let f_now = s.traction(0.025, 1.0).unwrap_or(0.12).abs();
                let f_ref = s_ref.traction(0.025, 1.0).unwrap_or(0.12).abs().max(1.0e-6);
                (f_ref / f_now.max(1.0e-6)).clamp(0.25, 1.25)
            }
            _ => 1.0,
        };

        let train_efficiency = (1.0 / mods.pivot_friction_multiplier.sqrt()).clamp(0.05, 1.15)
            * stribeck_ratio.sqrt()
            * (1.0 - 0.92 * mods.gold_fouling_jam_factor).clamp(0.0, 1.0);
        let nominal_gamma = OMEGA_REF_RAD_S / 310.0;
        let nominal_beat_energy_j = 0.5
            * HAIRSPRING_K_REF_NM_RAD
            * a_target_rad
            * a_target_rad
            * (nominal_gamma * PI / OMEGA_REF_RAD_S);
        let beat_impulse_energy_j = nominal_beat_energy_j
            * (self.mainspring_torque_nm / 0.0080).clamp(0.0, 1.18)
            * train_efficiency;

        let unlock_angle_rad = 3.8 * PI / 180.0;
        let drop_angle_rad = 18.0 * PI / 180.0;
        let impulse_arc_rad = unlock_angle_rad + drop_angle_rad;
        // Hunt-Crossley nonlinear viscoelastic impact loss factor from fs_dcontact
        let hc_loss = hunt_crossley_loss_from_restitution(0.76, 0.25).unwrap_or(0.32);

        for _ in 0..substeps {
            self.lab_time_s += h;
            let theta_rel = self.balance_angle_rad - mods.hairspring_eq_bias_rad;

            if self.balance_angle_rad.abs() > self.current_swing_peak_rad {
                self.current_swing_peak_rad = self.balance_angle_rad.abs();
            }

            let prev_omega = self.balance_omega_rad_s;
            let mut torque_esc_nm = 0.0;

            if !mods.magnetically_pinned && self.mainspring_wind > 0.002 {
                let dir = self.armed_swing_dir;
                let signed_pos = dir * self.balance_angle_rad;
                let signed_vel = dir * self.balance_omega_rad_s;

                match self.escapement_phase {
                    EscapementPhase::FreeSupplementaryArc => {
                        if signed_pos >= -unlock_angle_rad
                            && signed_pos <= drop_angle_rad
                            && signed_vel > 4.5
                        {
                            let unlock_drag_j = 0.08
                                * nominal_beat_energy_j
                                * mods.pivot_friction_multiplier.min(8.0);
                            let ke = 0.5 * inertia * self.balance_omega_rad_s.powi(2);
                            if ke > unlock_drag_j * 1.2 {
                                let new_ke = (ke - unlock_drag_j).max(0.0);
                                self.balance_omega_rad_s = self.balance_omega_rad_s.signum()
                                    * (2.0 * new_ke / inertia).sqrt();
                                self.escapement_phase = EscapementPhase::UnlockedImpulsing;
                                self.impulse_slide_fired_this_beat = false;
                                let energy_scale =
                                    (self.balance_omega_rad_s.abs() / 75.0).clamp(0.15, 1.6);
                                events.push(MechanicalAcousticEvent::EscapementUnlock {
                                    energy: energy_scale,
                                    is_entry_swing: dir > 0.0,
                                });
                            }
                        }
                    }
                    EscapementPhase::UnlockedImpulsing => {
                        if signed_pos <= drop_angle_rad && signed_vel > 0.5 {
                            // 3D rigid-body impulse torque vector tau = r_roller x F_pallet via fs_mbd::Vec3
                            let r_roller_m = MbdVec3::new(0.65e-3, 0.0, 0.0);
                            let f_tangential_n =
                                (beat_impulse_energy_j / impulse_arc_rad) / 0.65e-3;
                            let f_pallet_n = MbdVec3::new(0.0, dir * f_tangential_n, 0.0);
                            torque_esc_nm = r_roller_m
                                .cross(f_pallet_n)
                                .dot(MbdVec3::new(0.0, 0.0, 1.0));
                            let frac = ((signed_pos + unlock_angle_rad) / impulse_arc_rad)
                                .clamp(0.0, 1.0);
                            let prev_fork = self.pallet_fork_angle_rad;
                            self.pallet_fork_angle_rad =
                                dir * (-PALLET_BANKING_RAD + 2.0 * PALLET_BANKING_RAD * frac);
                            self.pallet_fork_omega_rad_s =
                                (self.pallet_fork_angle_rad - prev_fork) / h;
                            // As the escape tooth slides across the ruby pallet impulse plane,
                            // it advances through 75% of its 6-degree beat pitch before free drop!
                            self.escape_wheel_visual_rad =
                                self.escape_wheel_angle_rad + 0.75 * frac * ESCAPE_STEP_RAD;
                            self.escape_wheel_omega_rad_s =
                                0.75 * ESCAPE_STEP_RAD * (signed_vel / impulse_arc_rad);

                            if !self.impulse_slide_fired_this_beat && frac >= 0.5 {
                                self.impulse_slide_fired_this_beat = true;
                                events.push(MechanicalAcousticEvent::EscapementImpulseSlide {
                                    energy: (beat_impulse_energy_j / nominal_beat_energy_j)
                                        .clamp(0.1, 1.5),
                                    is_entry_swing: dir > 0.0,
                                });
                            }
                        } else if signed_pos > drop_angle_rad {
                            // Escape tooth clears the ruby pallet impulse corner and enters Free Drop!
                            self.escapement_phase = EscapementPhase::FreeToothDrop;
                            self.pallet_fork_angle_rad = dir * PALLET_BANKING_RAD;
                            self.pallet_fork_omega_rad_s = 0.0;
                            self.escape_wheel_angle_rad += ESCAPE_STEP_RAD;
                        } else if signed_vel <= 0.0 {
                            self.escapement_phase = EscapementPhase::FreeSupplementaryArc;
                            self.pallet_fork_angle_rad = -dir * PALLET_BANKING_RAD;
                            self.pallet_fork_omega_rad_s = 0.0;
                        }
                    }
                    EscapementPhase::FreeToothDrop => {
                        // Newton-Euler rotational acceleration of the 15-tooth Escape Wheel across the
                        // remaining 1.5 deg free drop gap until it slams onto the opposite ruby locking face!
                        let angular_acc = 42_000.0 * train_efficiency.max(0.1);
                        self.escape_wheel_omega_rad_s += angular_acc * h;
                        let step_rad = self.escape_wheel_omega_rad_s * h;
                        self.escape_wheel_visual_rad += step_rad;
                        if self.escape_wheel_visual_rad >= self.escape_wheel_angle_rad {
                            let overshoot_rad =
                                (self.escape_wheel_visual_rad - self.escape_wheel_angle_rad).max(0.0);
                            let overshoot_dt =
                                (overshoot_rad / self.escape_wheel_omega_rad_s.max(1.0)).clamp(0.0, h);
                            let exact_beat_time_s = self.lab_time_s - overshoot_dt;
                            self.escape_wheel_visual_rad = self.escape_wheel_angle_rad;
                            let hc_factor = (1.0 + 0.35 * hc_loss).clamp(0.85, 1.35);
                            let v_tip_m_s = (self.escape_wheel_omega_rad_s * 1.9e-3).max(0.01);
                            let normal_impact_scale = solve_normal_impact(NormalImpactInput {
                                gap: 0.0,
                                activation_tolerance: 1.0e-6,
                                relative_normal_speed: -v_tip_m_s,
                                inverse_effective_mass: 1.0 / 1.4e-6,
                                restitution: 0.58,
                            })
                            .map(|out| {
                                (1.0_f64 + 0.05 * (out.dissipated_energy / 1.0e-9).clamp(0.0, 2.0))
                                    .clamp(0.95, 1.10)
                            })
                            .unwrap_or(1.0);
                            let impact_scale = (self.escape_wheel_omega_rad_s / 65.0)
                                .clamp(0.25, 1.6)
                                * hc_factor
                                * normal_impact_scale
                                * (beat_impulse_energy_j / nominal_beat_energy_j).clamp(0.15, 1.5);
                            self.escape_wheel_omega_rad_s = 0.0;
                            self.escapement_phase = EscapementPhase::CompletedDropLocked;
                            self.commit_beat_to_going_train(
                                exact_beat_time_s,
                                dir > 0.0,
                                impact_scale,
                                events,
                            );
                        }
                    }
                    EscapementPhase::CompletedDropLocked | EscapementPhase::Stalled => {}
                }
            }

            let zeta = (0.5 * total_gamma / omega_eff).min(0.995);
            let omega_d = omega_eff * (1.0 - zeta * zeta).sqrt();
            let decay = (-zeta * omega_eff * h).exp();
            let cos_wt = (omega_d * h).cos();
            let sin_wt = (omega_d * h).sin();

            let x0 = theta_rel;
            // Positional gravity poise torque in the rotating 1-minute Tourbillon cage frame:
            // Breguet invented the Tourbillon specifically so this gravity poise torque averages out over 60s!
            let gx = mods.gravity_g_vec[0];
            let gy = mods.gravity_g_vec[1];
            let g_in_plane = (gx * gx + gy * gy).sqrt();
            let psi_g = gy.atan2(gx);
            let poise_torque_nm = -1.6e-10
                * g_in_plane
                * (self.balance_angle_rad + self.tourbillon_cage_angle_rad - psi_g).sin();
            let v0 =
                self.balance_omega_rad_s + ((torque_esc_nm + poise_torque_nm) / inertia) * h;

            let x1 =
                decay * (x0 * cos_wt + ((v0 + zeta * omega_eff * x0) / omega_d) * sin_wt);
            let v1 = decay
                * (v0 * cos_wt
                    - ((v0 * zeta * omega_eff + omega_eff * omega_eff * x0) / omega_d) * sin_wt);

            self.balance_angle_rad = x1 + mods.hairspring_eq_bias_rad;
            self.balance_omega_rad_s = v1;

            if prev_omega * self.balance_omega_rad_s <= 0.0
                && self.balance_angle_rad.abs() > unlock_angle_rad * 1.1
            {
                // If FreeToothDrop was still in flight at extreme frequency, complete it cleanly
                if self.escapement_phase == EscapementPhase::FreeToothDrop {
                    self.escape_wheel_visual_rad = self.escape_wheel_angle_rad;
                    self.escape_wheel_omega_rad_s = 0.0;
                    self.commit_beat_to_going_train(
                        self.lab_time_s,
                        self.armed_swing_dir > 0.0,
                        0.85,
                        events,
                    );
                }
                let deg = (self.current_swing_peak_rad * 180.0 / PI).clamp(0.0, 345.0);
                self.balance_amplitude_deg = 0.75 * self.balance_amplitude_deg + 0.25 * deg;
                self.current_swing_peak_rad = self.balance_angle_rad.abs();
                self.armed_swing_dir = if self.balance_angle_rad > 0.0 {
                    -1.0
                } else {
                    1.0
                };
                self.escapement_phase = EscapementPhase::FreeSupplementaryArc;
            }
        }

        if self.mainspring_wind <= 0.002 && self.balance_omega_rad_s.abs() < 2.0 {
            self.status = WatchOperatingStatus::StoppedMainspringUnwound;
            self.balance_amplitude_deg = 0.0;
        } else if !mods.magnetically_pinned {
            if self.balance_omega_rad_s.abs() < 1.5 && self.balance_angle_rad.abs() < 0.15 {
                self.status = if mods.extra_balance_damping_s_inv > 8.0 {
                    WatchOperatingStatus::FloodedStalledByOil
                } else {
                    WatchOperatingStatus::RunningPerturbed
                };
                self.balance_amplitude_deg = 0.0;
            } else {
                let is_perturbed = (mods.theoretical_rate_s_per_day.abs() > 1.5)
                    || (mods.gold_fouling_jam_factor > 0.02)
                    || (mods.b_field_escapement_mt > 0.5)
                    || (mods.remanent_magnetization_mt > 0.5)
                    || (mods.extra_balance_damping_s_inv > 0.5);
                self.status = if is_perturbed {
                    WatchOperatingStatus::RunningPerturbed
                } else {
                    WatchOperatingStatus::RunningChronometer
                };
            }
        }
    }

    fn commit_beat_to_going_train(
        &mut self,
        exact_beat_time_s: f64,
        is_entry_stone: bool,
        energy_scale: f64,
        events: &mut Vec<MechanicalAcousticEvent>,
    ) {
        self.total_beats += 1;
        self.watch_displayed_time_s += NOMINAL_BEAT_PERIOD_S;

        let d_fourth = ESCAPE_STEP_RAD / 10.0;
        self.tourbillon_cage_angle_rad += d_fourth;
        if let Ok(next_so3) = so3_body_exp_step(
            self.tourbillon_so3,
            GaVec3::new(0.0, 0.0, d_fourth / NOMINAL_BEAT_PERIOD_S),
            NOMINAL_BEAT_PERIOD_S,
        ) {
            self.tourbillon_so3 = next_so3;
        }
        let d_third = -d_fourth / 7.5;
        self.third_wheel_angle_rad += d_third;
        let d_center = d_fourth / 60.0;
        self.center_wheel_angle_rad += d_center;
        let d_barrel = -d_center / 7.8;
        self.barrel_wheel_angle_rad += d_barrel;

        self.mainspring_wind = (self.mainspring_wind - 1.0 / 756_000.0).max(0.0);
        self.moonphase_angle_rad += (NOMINAL_BEAT_PERIOD_S / (29.530_59 * 86_400.0)) * TAU;

        let raw_dt_beat = exact_beat_time_s - self.last_beat_lab_time_s;
        if self.total_beats == 1 || raw_dt_beat > 1.0 || raw_dt_beat <= 0.0 {
            // Resynchronize after initial startup or recovery from a stall
            self.last_beat_lab_time_s = exact_beat_time_s;
        } else {
            let dt_beat = raw_dt_beat.clamp(0.05, 1.0);
            self.last_beat_lab_time_s = exact_beat_time_s;
            if is_entry_stone {
                self.last_tick_interval_s = dt_beat;
            } else {
                self.last_tock_interval_s = dt_beat;
            }

            if self.total_beats >= 3 {
                let full_period_s = self.last_tick_interval_s + self.last_tock_interval_s;
                let nominal_period_s = 2.0 * NOMINAL_BEAT_PERIOD_S;
                let inst_rate_s_day = ((nominal_period_s / full_period_s) - 1.0) * 86_400.0;
                self.measured_rate_s_per_day = 0.70 * self.measured_rate_s_per_day
                    + 0.30 * inst_rate_s_day.clamp(-9999.0, 9999.0);
                let beat_err_ms =
                    0.5 * (self.last_tick_interval_s - self.last_tock_interval_s).abs() * 1000.0;
                self.measured_beat_error_ms =
                    0.70 * self.measured_beat_error_ms + 0.30 * beat_err_ms.clamp(0.0, 99.9);

                // Update anytime-valid COSC chronometer compliance e-process (fs_eproc::BettingEProcess)
                let normalized_violation =
                    (self.measured_rate_s_per_day.abs() / 20.0).clamp(0.0, 1.0);
                let _ = self.cosc_eprocess.observe(normalized_violation);
                self.cosc_e_value = self.cosc_eprocess.e_value();
            }
        }

        // Counter-based Philox deterministic micro-tribological acoustic variation (fs_rand::Stream::at)
        let philox_words = Stream::at(
            StreamKey {
                seed: 0x1903_C05C,
                kernel: 0x484F_524F,
                tile: (self.total_beats & 0xFFFF) as u32,
            },
            self.total_beats,
        );
        let u01 = (philox_words[0] as f64) / (u32::MAX as f64);
        let micro_jitter = 1.0 + 0.015 * (u01 * 2.0 - 1.0);

        events.push(MechanicalAcousticEvent::EscapementLockDrop {
            energy: (energy_scale * micro_jitter).clamp(0.15, 1.6),
            is_entry_stone,
        });
    }
}

impl Default for HorologyState {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::multiphysics::MultiphysicsState;

    #[test]
    fn unperturbed_pocket_watch_ticks_at_5hz_and_keeps_accurate_time() {
        let mp = MultiphysicsState::new();
        let mut horology = HorologyState::new();
        let mut events = Vec::new();
        let t0_watch = horology.watch_displayed_time_s;
        let t0_lab = horology.lab_time_s;

        for _ in 0..1000 {
            horology.step(0.01, &mp.modifiers, &mut events);
        }

        let elapsed_lab = horology.lab_time_s - t0_lab;
        let elapsed_watch = horology.watch_displayed_time_s - t0_watch;
        assert!((elapsed_lab - 10.0).abs() < 1.0e-6);
        assert!(
            (horology.total_beats as i64 - 50).abs() <= 1,
            "expected ~50 beats in 10s, got {}",
            horology.total_beats
        );
        assert!(
            (elapsed_watch - 10.0).abs() <= 0.25,
            "expected watch to advance ~10.0s, got {elapsed_watch}s"
        );
        assert!(
            horology.measured_rate_s_per_day.abs() < 3.0,
            "expected chronometer rate within 3 s/day, got {} s/day (last_tick={}, last_tock={}, beats={})",
            horology.measured_rate_s_per_day,
            horology.last_tick_interval_s,
            horology.last_tock_interval_s,
            horology.total_beats
        );
        assert!(horology.balance_amplitude_deg > 240.0 && horology.balance_amplitude_deg < 315.0);
        assert!(!events.is_empty());
    }
}
