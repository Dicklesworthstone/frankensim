//! Convergence-architecture `fs_ports::component::Component` implementations for the
//! industrial steam turbine simulation.
//!
//! Every physics domain owns its own accepted state, typed ports (`Power`, `Stream`,
//! `Heat`, `Signal`), and pure `propose` transition, and is composed exclusively via
//! `fs_system::SystemBuilder` (`bond`, `common_effort`, `common_flow`, `heat_route`,
//! `signal`, `boundary`). No pairing modules or ad-hoc `*_with_*` entry points exist.
//!
//! Components implemented here:
//! 1. [`SteamControlValveComponent`] — live-steam CV/ESV & extraction throttling valves (IF97 + Stodola)
//! 2. [`SteamPlenumComponent`] — unsteady IF97 steam control volumes (HP chest, 12-bar extraction header)
//! 3. [`TurbineStageGroupComponent`] — HP (Curtis + 5 stages) and IP/LP (8 IP + 3 LP stages):
//!    IF97 expansion, velocity-triangle efficiency, Baumann wetness penalty, windage heating,
//!    balance-piston axial thrust, and Alford cross-coupled steam whirl
//! 4. [`SurfaceCondenserComponent`] — LP exhaust shell coupled to `fs_multiphase::condensation::CondensationModelV1`
//!    (Nusselt horizontal-tube film condensation + Colburn-Hougen air ingress) and
//!    `fs_convection::evaluate(CorrelationId::Gnielinski)` cooling-water tube side
//! 5. [`RotorShaftComponent`] — exact midpoint torsional inertia `J_p` from `AxisymmetricChart::mass_properties`
//! 6. [`JournalBearingComponent`] — front radial/thrust + rear radial hydrodynamic bearings
//!    (`fs_tribo::journal::JournalGeometry` + Ocvirk Sommerfeld + Shell Oil S2)
//! 7. [`SynchronousGeneratorComponent`] — 2-pole 50 Hz turbo-generator backed by
//!    `fs_phs::circuit::synchronous::SynchronousMotor` + M19 Steinmetz core & IACS copper heat losses
//! 8. [`TurbineCasingThermalComponent`] — 8-node metal/oil thermal mass network with
//!    `fs_convection::evaluate(CorrelationId::ChurchillChuHorizontalCylinder)`, differential
//!    axial expansion, hoop thermal stress, and `fs_material::creep::NortonCreep`
//! 9. [`GovernorControlComponent`] — 4% speed droop + `fs_control::speed::SpeedController` PI +
//!    `fs_control::thermal::ThermalDeratingPolicy` trip latches + 110% overspeed ESV trip

use core::f64::consts::PI;
use core::num::NonZeroUsize;
use std::cell::{Cell, RefCell};

use fs_control::speed::{
    ShaftPowerLimits, SpeedCommand, SpeedController, SpeedControllerParameters,
    SpeedControllerState,
};
use fs_control::thermal::{ThermalDeratingBand, ThermalDeratingPolicy, ThermalDeratingState};
use fs_convection::{CorrelationId, CorrelationInputs};
use fs_exec::Cx;
use fs_multiphase::condensation::{
    CondensationGeometryV1, CondensationModelV1, NonCondensableGasV1,
};
use fs_phs::circuit::synchronous::{SynchronousMotor, SynchronousMotorParameters};
use fs_ports::component::{
    Class, Component, PortRole, PortSpec, PortTangent, Proposal, Receipt, STREAM_CARRIED,
    STREAM_ENTHALPY_RATE, STREAM_MASS_RATE, STREAM_PRESSURE, STREAM_SPECIFIC_ENTHALPY,
};
use fs_ports::{PortKind, StableId};
use fs_qty::{Length, Pressure, Temperature};
use fs_thermochem::water::if97::If97SaturationPair;

use crate::materials::{STEFAN_BOLTZMANN, T_REF_K, library};
use crate::rotordynamics::{BearingDynamicState, evaluate_journal_bearing};
use crate::steam::{SteamPoint, eval_ph, eval_ps, eval_uv, stodola_mass_flow, water_if97_v1};

/// Typed error returned by turbine `Component`s.
#[derive(Debug, Clone, PartialEq)]
pub enum TurbineComponentError {
    Cancelled,
    InvalidState(&'static str),
}

#[inline]
fn sid(name: &str) -> StableId {
    StableId::new(name).expect("valid static StableId")
}

#[inline]
fn nz(n: usize) -> NonZeroUsize {
    NonZeroUsize::new(n).expect("nonzero port width")
}

#[inline]
fn check_cx(cx: &Cx<'_>) -> Result<(), TurbineComponentError> {
    if cx.checkpoint().is_err() {
        Err(TurbineComponentError::Cancelled)
    } else {
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// 1. SteamControlValveComponent (Live Steam CV/ESV & Controlled Extraction Valve)
// ---------------------------------------------------------------------------

/// Stateless isenthalpic steam control valve with Stodola choking and exact
/// IF97 throttling entropy generation.
pub struct SteamControlValveComponent {
    ports: [PortSpec; 3],
    /// Maximum effective Stodola swallowing area at 100% valve lift [m²].
    pub k_max_m2: f64,
    /// Critical sonic pressure ratio `p_out / p_in` (~0.546 for superheated steam).
    pub r_crit: f64,
}

#[derive(Debug, Clone, Copy, Default)]
pub struct ValveState {
    pub m_dot_kg_s: f64,
    pub lift: f64,
    pub dp_pa: f64,
    pub s_gen_w_k: f64,
}

impl SteamControlValveComponent {
    #[must_use]
    pub fn new(k_max_m2: f64, r_crit: f64) -> Self {
        let basis = water_if97_v1();
        Self {
            ports: [
                PortSpec::identified_stream(sid("inlet"), PortRole::EffortIn, basis.clone())
                    .expect("valve inlet stream port"),
                PortSpec::identified_stream(sid("outlet"), PortRole::EffortIn, basis)
                    .expect("valve outlet stream port"),
                PortSpec::signal_in(sid("lift"), nz(1)),
            ],
            k_max_m2,
            r_crit,
        }
    }
}

impl Component for SteamControlValveComponent {
    type State = ValveState;
    type Error = TurbineComponentError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Smooth
    }

    fn stored_energy_j(&self, _state: &Self::State) -> Result<f64, Self::Error> {
        Ok(0.0)
    }

    fn stored_mass_kg(&self, _state: &Self::State) -> Result<Option<f64>, Self::Error> {
        Ok(Some(0.0))
    }

    fn stored_carried_quantity(
        &self,
        _state: &Self::State,
        _carried: usize,
    ) -> Result<Option<f64>, Self::Error> {
        Ok(Some(0.0))
    }

    fn propose(
        &self,
        _old: &Self::State,
        dt_s: f64,
        inputs: &[f64],
        new: &mut Self::State,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, Self::Error> {
        check_cx(cx)?;
        // Inputs layout: inlet [P, h, x] (0..3), outlet [P, h, x] (3..6), lift [u] (6)
        let p_in = inputs[STREAM_PRESSURE];
        let h_in = inputs[STREAM_SPECIFIC_ENTHALPY];
        let p_out = inputs[3 + STREAM_PRESSURE];
        let h_out = inputs[3 + STREAM_SPECIFIC_ENTHALPY];
        let lift = inputs[6].clamp(0.0, 1.0);

        // Smooth equal-percentage-like valve characteristic with tiny closed-seat leak floor
        // so the Newton Jacobian stays nonsingular even after an ESV trip!
        let cv_frac = 1.0e-4 + (1.0 - 1.0e-4) * (0.25 * lift + 0.75 * lift * lift);
        let k_eff = self.k_max_m2 * cv_frac;

        let st_in = eval_ph(p_in, h_in);
        let m_dot = stodola_mass_flow(k_eff, p_in, st_in.v_m3_kg, p_out, self.r_crit);
        let h_donor = if m_dot >= 0.0 { h_in } else { h_out };

        // Positive rates enter the component: +m_dot at inlet, -m_dot at outlet.
        outputs[STREAM_MASS_RATE] = m_dot;
        outputs[STREAM_ENTHALPY_RATE] = m_dot * h_donor;
        outputs[STREAM_CARRIED] = m_dot;

        outputs[3 + STREAM_MASS_RATE] = -m_dot;
        outputs[3 + STREAM_ENTHALPY_RATE] = -m_dot * h_donor;
        outputs[3 + STREAM_CARRIED] = -m_dot;

        // Exact Second-Law throttling entropy production: s_dot = m_dot * (s_out - s_in) >= 0
        let s_gen_w_k = if m_dot > 0.0 && p_in > p_out {
            let st_out = eval_ph(p_out, h_donor);
            m_dot * (st_out.s_j_kg_k - st_in.s_j_kg_k)
        } else {
            0.0
        };

        *new = ValveState {
            m_dot_kg_s: m_dot,
            lift,
            dp_pa: (p_in - p_out).max(0.0),
            s_gen_w_k,
        };

        let e_scale = (dt_s * m_dot.abs() * h_donor.abs()).max(1.0);
        Ok(Proposal::Step(Receipt {
            stored_change_j: 0.0,
            dissipated_j: 0.0,
            numerical_loss_j: 0.0,
            entropy_generated_j_per_k: dt_s * s_gen_w_k,
            energy_tolerance_j: 1.0e-9 * e_scale + 1.0e-6,
            residual: 0.0,
            tolerance: 0.0,
        }))
    }
}

// ---------------------------------------------------------------------------
// 2. SteamPlenumComponent (HP Steam Chest & 12-Bar Controlled Extraction Header)
// ---------------------------------------------------------------------------

/// Unsteady IAPWS-IF97 steam control volume with exact mass, H2O species, and
/// internal-energy storage.
pub struct SteamPlenumComponent {
    ports: [PortSpec; 2],
    pub volume_m3: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct PlenumState {
    pub mass_kg: f64,
    pub total_u_j: f64,
    pub p_pa: f64,
    pub t_k: f64,
    pub h_j_kg: f64,
    pub s_j_kg_k: f64,
    pub quality_x: f64,
}

impl PlenumState {
    #[must_use]
    pub fn from_steam_point(st: SteamPoint, volume_m3: f64) -> Self {
        let mass_kg = (volume_m3 / st.v_m3_kg).max(1.0e-3);
        Self {
            mass_kg,
            total_u_j: mass_kg * st.u_j_kg,
            p_pa: st.p_pa,
            t_k: st.t_k,
            h_j_kg: st.h_j_kg,
            s_j_kg_k: st.s_j_kg_k,
            quality_x: st.quality_x,
        }
    }
}

impl SteamPlenumComponent {
    #[must_use]
    pub fn new(volume_m3: f64) -> Self {
        Self {
            ports: [
                PortSpec::identified_stream(sid("steam"), PortRole::FlowIn, water_if97_v1())
                    .expect("plenum steam port"),
                PortSpec::heat(sid("wall_heat"), PortRole::FlowIn).expect("plenum heat port"),
            ],
            volume_m3,
        }
    }
}

impl Component for SteamPlenumComponent {
    type State = PlenumState;
    type Error = TurbineComponentError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Smooth
    }

    fn stored_energy_j(&self, state: &Self::State) -> Result<f64, Self::Error> {
        Ok(state.total_u_j)
    }

    fn stored_mass_kg(&self, state: &Self::State) -> Result<Option<f64>, Self::Error> {
        Ok(Some(state.mass_kg))
    }

    fn stored_carried_quantity(
        &self,
        state: &Self::State,
        _carried: usize,
    ) -> Result<Option<f64>, Self::Error> {
        Ok(Some(state.mass_kg))
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
        check_cx(cx)?;
        // Inputs: steam [m_dot_in, H_dot_in, m_h2o_in] (0..3), wall_heat [Q_in] (3)
        let m_dot_in = inputs[STREAM_MASS_RATE];
        let h_dot_in = inputs[STREAM_ENTHALPY_RATE];
        let q_in = inputs[3];

        let m_new = old.mass_kg + dt_s * m_dot_in;
        // Deep vacuum (0.8 kPa, v ≈ 160 m³/kg) in a 0.30 m³ chest holds ~1.8 g of steam,
        // so the positivity guard must sit well below vacuum equilibrium (~0.1 g).
        let m_min = (self.volume_m3 / 2_000.0).min(1.0e-4);
        if m_new <= m_min || !m_new.is_finite() {
            return Ok(Proposal::Shorten { dt_s: 0.5 * dt_s });
        }
        let u_new_total = old.total_u_j + dt_s * (h_dot_in + q_in);
        if u_new_total <= 50.0 || !u_new_total.is_finite() {
            return Ok(Proposal::Shorten { dt_s: 0.5 * dt_s });
        }

        let u_spec = u_new_total / m_new;
        let v_spec = self.volume_m3 / m_new;
        let st_new = eval_uv(u_spec, v_spec, old.p_pa);

        *new = PlenumState {
            mass_kg: m_new,
            total_u_j: u_new_total,
            p_pa: st_new.p_pa,
            t_k: st_new.t_k,
            h_j_kg: st_new.h_j_kg,
            s_j_kg_k: st_new.s_j_kg_k,
            quality_x: st_new.quality_x,
        };

        // Implicit-endpoint (or midpoint) state returned to connected transport elements
        let p_mid = 0.5 * (old.p_pa + st_new.p_pa);
        let h_mid = 0.5 * (old.h_j_kg + st_new.h_j_kg);
        let t_mid = 0.5 * (old.t_k + st_new.t_k);
        outputs[STREAM_PRESSURE] = p_mid;
        outputs[STREAM_SPECIFIC_ENTHALPY] = h_mid;
        outputs[STREAM_CARRIED] = 1.0;
        outputs[3] = t_mid;

        let du = u_new_total - old.total_u_j;
        let e_tol = 1.0e-7 * (u_new_total.abs() + old.total_u_j.abs() + 1.0);
        Ok(Proposal::Step(Receipt {
            stored_change_j: du,
            dissipated_j: 0.0,
            numerical_loss_j: 0.0,
            entropy_generated_j_per_k: 0.0,
            energy_tolerance_j: e_tol,
            residual: 0.0,
            tolerance: 0.0,
        }))
    }
}

// ---------------------------------------------------------------------------
// 3. TurbineStageGroupComponent (HP Curtis+5 Stages & IP/LP 8+3 Stages)
// ---------------------------------------------------------------------------

/// Multi-stage steam turbine expansion section implementing real IAPWS-IF97 expansion,
/// velocity-triangle efficiency, Baumann Region-4 moisture loss, labyrinth leakage,
/// low-flow windage heating, twin-balance-piston axial thrust, and Alford steam whirl.
///
/// The section's isentropic efficiency is not a parameter: it arrives on the `eta` signal
/// port, sampled each 20 ms sub-step from the 17-stage row-by-row meanline
/// ([`crate::stage_meanline::StageMeanlineReport::section_efficiencies`]) at the plant's
/// own operating point, so speed, wetness, blade geometry and tip clearance reach the
/// power through one loss model.
pub struct TurbineStageGroupComponent {
    ports: [PortSpec; 5],
    pub n_stages: u32,
    pub k_stodola_m2: f64,
    pub r_crit: f64,
    pub r_mean_m: f64,
    pub mean_blade_h_m: f64,
    pub leakage_fraction: f64,
    pub k_windage: f64,
    pub ua_wall_w_k: f64,
    pub m_dot_nom_kg_s: f64,
    pub eff_thrust_area_m2: f64,
    pub balance_piston_area_m2: f64,
    pub alford_beta: f64,
}

#[derive(Debug, Clone, Copy, Default)]
pub struct StageGroupState {
    pub m_dot_kg_s: f64,
    pub eta_is: f64,
    pub shaft_power_w: f64,
    pub shaft_torque_nm: f64,
    pub axial_thrust_n: f64,
    pub alford_kxy_n_per_m: f64,
    pub t_in_k: f64,
    pub t_out_k: f64,
    pub p_in_pa: f64,
    pub p_out_pa: f64,
    pub h_in_j_kg: f64,
    pub h_out_j_kg: f64,
    pub s_in_j_kg_k: f64,
    pub s_out_j_kg_k: f64,
    pub quality_out: f64,
    pub leakage_kg_s: f64,
    pub windage_w: f64,
    pub q_to_metal_w: f64,
    /// Entropy produced in the section over the last step [W/K] (heat leaves at the wall temperature).
    pub s_gen_w_k: f64,
}

impl TurbineStageGroupComponent {
    /// HP section: Curtis impulse control wheel + 5 HP reaction stages (`r_mean = 0.318 m`,
    /// twin balance pistons at `r = 0.258 m, 0.304 m` from `geometry::build_rotor_chart()`).
    #[must_use]
    pub fn hp_section() -> Self {
        let basis = water_if97_v1();
        Self {
            ports: [
                PortSpec::identified_stream(sid("inlet"), PortRole::EffortIn, basis.clone())
                    .expect("hp inlet"),
                PortSpec::identified_stream(sid("outlet"), PortRole::EffortIn, basis)
                    .expect("hp outlet"),
                PortSpec::power(
                    sid("shaft"),
                    PortKind::RotationalTorqueAngularVelocity,
                    PortRole::FlowIn,
                    nz(1),
                )
                .expect("hp shaft"),
                PortSpec::heat(sid("casing_heat"), PortRole::EffortIn).expect("hp heat"),
                PortSpec::signal_in(sid("eta"), nz(1)),
            ],
            n_stages: 6,
            k_stodola_m2: 0.003_45,
            r_crit: 0.14,
            r_mean_m: 0.318,
            mean_blade_h_m: 0.072,
            leakage_fraction: 0.018,
            k_windage: 2.2e-4,
            ua_wall_w_k: 4_800.0,
            m_dot_nom_kg_s: 52.0,
            eff_thrust_area_m2: 0.245,
            // Twin balance pistons (r=0.304 m vs r=0.156 m gland -> annular area ~0.214 m²)
            balance_piston_area_m2: PI * (0.304_f64.powi(2) - 0.156_f64.powi(2)),
            alford_beta: 2.2,
        }
    }

    /// IP/LP section: 8 IP reaction stages + 3 LP stages ending in the giant 0.935 m last stage.
    #[must_use]
    pub fn ip_lp_section() -> Self {
        let basis = water_if97_v1();
        Self {
            ports: [
                PortSpec::identified_stream(sid("inlet"), PortRole::EffortIn, basis.clone())
                    .expect("iplp inlet"),
                PortSpec::identified_stream(sid("outlet"), PortRole::EffortIn, basis)
                    .expect("iplp outlet"),
                PortSpec::power(
                    sid("shaft"),
                    PortKind::RotationalTorqueAngularVelocity,
                    PortRole::FlowIn,
                    nz(1),
                )
                .expect("iplp shaft"),
                PortSpec::heat(sid("casing_heat"), PortRole::EffortIn).expect("iplp heat"),
                PortSpec::signal_in(sid("eta"), nz(1)),
            ],
            n_stages: 11,
            k_stodola_m2: 0.0180,
            r_crit: 0.06,
            r_mean_m: 0.485,
            mean_blade_h_m: 0.265,
            leakage_fraction: 0.014,
            k_windage: 2.6e-3, // Large LP last-stage blades churn heavily at zero flow!
            ua_wall_w_k: 5_200.0,
            m_dot_nom_kg_s: 42.0,
            eff_thrust_area_m2: 0.360,
            balance_piston_area_m2: 0.195,
            alford_beta: 1.6,
        }
    }

    /// HP section parameterized by live engineer-editable [`crate::stage_meanline::TurbineDesignParams`].
    #[must_use]
    pub fn hp_section_with_design(d: &crate::stage_meanline::TurbineDesignParams) -> Self {
        let mut s = Self::hp_section();
        let d = d.clamped();
        s.leakage_fraction = d.hp_leakage_fraction();
        s.eff_thrust_area_m2 *= d.thrust_reaction_factor();
        s
    }

    /// IP/LP section parameterized by live engineer-editable [`crate::stage_meanline::TurbineDesignParams`].
    #[must_use]
    pub fn ip_lp_section_with_design(d: &crate::stage_meanline::TurbineDesignParams) -> Self {
        let mut s = Self::ip_lp_section();
        let d = d.clamped();
        s.leakage_fraction = d.iplp_leakage_fraction();
        s.eff_thrust_area_m2 *= d.thrust_reaction_factor();
        s
    }
}

impl Component for TurbineStageGroupComponent {
    type State = StageGroupState;
    type Error = TurbineComponentError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Smooth
    }

    fn stored_energy_j(&self, _state: &Self::State) -> Result<f64, Self::Error> {
        Ok(0.0)
    }

    fn stored_mass_kg(&self, _state: &Self::State) -> Result<Option<f64>, Self::Error> {
        Ok(Some(0.0))
    }

    fn stored_carried_quantity(
        &self,
        _state: &Self::State,
        _carried: usize,
    ) -> Result<Option<f64>, Self::Error> {
        Ok(Some(0.0))
    }

    fn propose(
        &self,
        _old: &Self::State,
        dt_s: f64,
        inputs: &[f64],
        new: &mut Self::State,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, Self::Error> {
        check_cx(cx)?;
        // Inputs: inlet (0..3), outlet (3..6), shaft omega_in = -omega_shaft (6), casing T_wall (7),
        // section isentropic efficiency from the meanline (8)
        let p_in = inputs[STREAM_PRESSURE];
        let h_in = inputs[STREAM_SPECIFIC_ENTHALPY];
        let p_out = inputs[3 + STREAM_PRESSURE];
        let h_out_down = inputs[3 + STREAM_SPECIFIC_ENTHALPY];
        let omega_in = inputs[6];
        let omega_shaft = -omega_in;
        let t_wall = inputs[7].clamp(280.0, 1_100.0);
        let eta_is = inputs[8];
        if !(eta_is > 0.0 && eta_is < 1.0) {
            return Err(TurbineComponentError::InvalidState(
                "stage section efficiency from the meanline must lie in (0, 1)",
            ));
        }

        let st_in = eval_ph(p_in, h_in);
        let m_dot = stodola_mass_flow(self.k_stodola_m2, p_in, st_in.v_m3_kg, p_out, self.r_crit);
        let st_is = eval_ps(p_out, st_in.s_j_kg_k);
        let dh_is = (h_in - st_is.h_j_kg).max(0.0);

        let c0 = (2.0 * dh_is / f64::from(self.n_stages)).sqrt().max(10.0);

        let m_pos = m_dot.max(0.0);
        let leakage = m_pos * self.leakage_fraction;
        let m_work = m_pos - leakage;
        let p_gross = m_work * eta_is * dh_is;

        // Low-flow blade windage heating (turns kinetic churning power into steam enthalpy!)
        let flow_shield = 1.0 / (1.0 + (m_pos / 5.0).powi(2));
        let p_windage = self.k_windage * st_is.rho_kg_m3 * omega_shaft.abs().powi(3) * flow_shield;

        // Compute driving and windage torques cleanly across all (m_dot, omega) with zero singularity
        let denom = omega_shaft
            .abs()
            .max(2.0)
            .copysign(if omega_shaft >= 0.0 { 1.0 } else { -1.0 });
        let tau_wind = -p_windage / denom;
        let tau_drive = if omega_shaft.abs() > 2.0 {
            p_gross / omega_shaft
        } else {
            // Starting stall torque proportional to momentum flux of steam jet
            (m_work * self.r_mean_m * c0 * 0.35).min(2.2e5)
        };
        let tau_exact = tau_drive + tau_wind;
        let p_gross_eff = tau_drive * omega_shaft;
        let p_wind_eff = -tau_wind * omega_shaft; // >= 0 everywhere; == p_windage for |omega| >= 2
        let p_shaft_exact = tau_exact * omega_shaft; // == p_gross_eff - p_wind_eff identically

        // Convective heat transfer between steam path and surrounding inner casing wall
        let t_steam_avg = 0.5 * (st_in.t_k + st_is.t_k);
        let flow_scale = 0.18 + 0.82 * (m_pos / self.m_dot_nom_kg_s).clamp(0.0, 2.0).powf(0.8);
        let q_to_metal = self.ua_wall_w_k * flow_scale * (t_steam_avg - t_wall);

        // C^1-smooth First-Law energy balance across m_dot = 0:
        // Blend weight w = m_pos^2 / (m_pos^2 + m_scale^2) smoothly transfers windage heating
        // from through-flow steam (at normal flow) directly to the casing wall (at zero flow)
        // with zero discontinuity and zero 1/m_dot singularity:
        //   q_target = w * q_to_metal + (1 - w) * p_wind_eff
        //   m_dot * (h_in - h_out) = p_shaft_exact + q_target
        //                          = p_gross_eff + w * (q_to_metal - p_wind_eff)
        let m_scale = 0.5_f64;
        let inv_m_smooth = m_pos / (m_pos * m_pos + m_scale * m_scale);
        let (h_in_port, h_out, q_exact) = if m_dot > 0.0 {
            let dh_drive = p_gross_eff / m_dot;
            let h_raw = h_in - dh_drive - inv_m_smooth * (q_to_metal - p_wind_eff);
            let h_clamped = h_raw.clamp(120_000.0, 4_050_000.0);
            // Exact First-Law closure (equals q_target to machine epsilon whenever unclamped)
            let q_fin = m_dot * (h_in - h_clamped) - p_shaft_exact;
            (h_in, h_clamped, q_fin)
        } else {
            // Zero or slight reverse flow: donor enthalpy h_out_down cancels identically across ports
            // and 100% of windage power heats the casing wall (C^1 continuous with m_dot -> 0+!)
            (h_out_down, h_out_down, p_wind_eff)
        };

        let st_out = eval_ph(p_out, h_out);

        // Port outputs (all rates positive INTO this component):
        outputs[STREAM_MASS_RATE] = m_dot;
        outputs[STREAM_ENTHALPY_RATE] = m_dot * h_in_port;
        outputs[STREAM_CARRIED] = m_dot;

        outputs[3 + STREAM_MASS_RATE] = -m_dot;
        outputs[3 + STREAM_ENTHALPY_RATE] = -m_dot * h_out;
        outputs[3 + STREAM_CARRIED] = -m_dot;

        // Shaft port (FlowIn): returns torque exerted ON the shaft (so power into stage = omega_in * tau = -p_shaft)
        outputs[6] = tau_exact;
        // Casing heat port (EffortIn): heat entering stage = -q_exact (heat leaves stage to metal)
        outputs[7] = -q_exact;

        // Axial thrust minus balance-piston counter-thrust [N]
        let axial_thrust_n = (p_in - p_out) * (self.eff_thrust_area_m2 - self.balance_piston_area_m2);
        // Alford cross-coupled aerodynamic steam-whirl stiffness [N/m]
        let alford_kxy =
            self.alford_beta * tau_exact.abs() / (2.0 * self.r_mean_m * self.mean_blade_h_m);

        // Internal aerodynamic + moisture dissipation & Second-Law entropy generation
        let aero_diss_w = (m_pos * (h_out - st_is.h_j_kg).max(0.0) + p_windage).max(0.0);
        let s_gen_w_k = m_pos * (st_out.s_j_kg_k - st_in.s_j_kg_k) + q_exact / t_wall;

        *new = StageGroupState {
            m_dot_kg_s: m_dot,
            eta_is,
            shaft_power_w: p_shaft_exact,
            shaft_torque_nm: tau_exact,
            axial_thrust_n,
            alford_kxy_n_per_m: alford_kxy,
            t_in_k: st_in.t_k,
            t_out_k: st_out.t_k,
            p_in_pa: p_in,
            p_out_pa: p_out,
            h_in_j_kg: h_in,
            h_out_j_kg: h_out,
            s_in_j_kg_k: st_in.s_j_kg_k,
            s_out_j_kg_k: st_out.s_j_kg_k,
            quality_out: st_out.quality_x,
            leakage_kg_s: leakage,
            windage_w: p_windage,
            q_to_metal_w: q_exact,
            s_gen_w_k,
        };

        let e_scale = (dt_s * (m_dot.abs() * h_in.abs() + p_shaft_exact.abs() + q_exact.abs())).max(1.0);
        Ok(Proposal::Step(Receipt {
            stored_change_j: 0.0,
            dissipated_j: dt_s * aero_diss_w,
            numerical_loss_j: 0.0,
            entropy_generated_j_per_k: dt_s * s_gen_w_k,
            energy_tolerance_j: 1.0e-8 * e_scale + 1.0e-4,
            residual: 0.0,
            tolerance: 0.0,
        }))
    }
}

// ---------------------------------------------------------------------------
// 4. SurfaceCondenserComponent (fs-multiphase CondensationModelV1 + fs-convection)
// ---------------------------------------------------------------------------

/// Surface condenser vacuum shell combining:
/// - `fs_multiphase::condensation::CondensationModelV1::steam_air` on horizontal tubes
///   (Nusselt filmwise condensation + Colburn-Hougen non-condensable air blanketing!)
/// - `fs_convection::evaluate(CorrelationId::Gnielinski)` for turbulent cooling-water flow
///   inside 3,200 condenser tubes (`D_i = 0.023 m`, `D_o = 0.0254 m`, `A_total = 1,550 m²`).
pub struct SurfaceCondenserComponent {
    ports: [PortSpec; 4],
    pub vapor_volume_m3: f64,
    pub tube_outer_diam_m: f64,
    pub tube_inner_diam_m: f64,
    pub tube_area_m2: f64,
    pub tube_count: f64,
    pub cw_nom_flow_kg_s: f64,
    /// Exact memoization of `evaluate_heat_transfer` keyed by `(p_cond, t_cw_in, cw_flow, w_air)`
    /// so repeated Newton finite-difference calls within the same step skip redundant 26-step
    /// shell-side `CondensationModelV1` solves (bit-exact identical result, ~40x speedup).
    ht_cache: Cell<((f64, f64, f64, f64), (f64, f64, f64, f64, f64, f64))>,
}

#[derive(Debug, Clone, Copy)]
pub struct CondenserState {
    pub vapor_mass_kg: f64,
    pub total_u_j: f64,
    pub p_cond_pa: f64,
    pub t_sat_k: f64,
    pub t_wall_k: f64,
    pub t_cw_out_k: f64,
    pub q_cond_w: f64,
    pub h_cond_w_m2_k: f64,
    pub h_cw_w_m2_k: f64,
    pub u_overall_w_m2_k: f64,
    pub m_condensate_kg_s: f64,
    pub air_mass_fraction: f64,
    /// Entropy produced inside the shell over the last step [W/K], from the exact
    /// control-volume balance (heat leaves at the cooling water's mean temperature).
    pub s_gen_w_k: f64,
}

impl CondenserState {
    #[must_use]
    pub fn at_vacuum_pressure(p_cond_pa: f64, vapor_volume_m3: f64) -> Self {
        let pair = If97SaturationPair::at_pressure(Pressure::new(p_cond_pa.clamp(1_500.0, 2.0e5)))
            .expect("IF97 condenser saturation pair");
        let vap = pair.vapor();
        let rho_g = 1.0 / vap.specific_volume().value();
        let mass = (vapor_volume_m3 * rho_g).max(0.2);
        let u_g = vap.specific_internal_energy().value();
        Self {
            vapor_mass_kg: mass,
            total_u_j: mass * u_g,
            p_cond_pa,
            t_sat_k: vap.temperature().value(),
            t_wall_k: vap.temperature().value() - 6.0,
            t_cw_out_k: 299.15,
            q_cond_w: 85.0e6,
            h_cond_w_m2_k: 11_000.0,
            h_cw_w_m2_k: 8_500.0,
            u_overall_w_m2_k: 3_600.0,
            m_condensate_kg_s: 38.0,
            air_mass_fraction: 0.0005,
            s_gen_w_k: 0.0,
        }
    }
}

/// Pressure range of the condenser shell model [Pa]: from just above the triple point
/// (IF97 Region 4 starts at 611.657 Pa) to a 4 bar relief limit.
const COND_P_MIN_PA: f64 = 700.0;
const COND_P_MAX_PA: f64 = 4.0e5;

/// Energy function of the saturated vapour space, `G(p) = V rho_g(p) (u_g(p) - h_f)`,
/// with the liquid enthalpy `h_f` held at its start-of-step value.
fn condenser_energy_function(
    volume_m3: f64,
    h_f_frozen: f64,
    p_pa: f64,
) -> (f64, If97SaturationPair) {
    let pair = If97SaturationPair::at_pressure(Pressure::new(p_pa))
        .expect("IF97 condenser saturation pair");
    let vap = pair.vapor();
    let rho_g = 1.0 / vap.specific_volume().value();
    let u_g = vap.specific_internal_energy().value();
    (volume_m3 * rho_g * (u_g - h_f_frozen), pair)
}

/// Solve `G(p) = g_target` for the shell pressure by safeguarded secant iteration in
/// `ln p` (G rises monotonically with p). Returns `None` when the target lies outside
/// `[COND_P_MIN_PA, COND_P_MAX_PA]`, so the caller shortens the step instead of clamping.
fn solve_condenser_pressure(
    volume_m3: f64,
    h_f_frozen: f64,
    g_target: f64,
    p_guess_pa: f64,
) -> Option<(f64, If97SaturationPair)> {
    if !(g_target.is_finite() && g_target > 0.0) {
        return None;
    }
    let ln_g_t = g_target.ln();
    let f = |ln_p: f64| {
        let (g, pair) = condenser_energy_function(volume_m3, h_f_frozen, ln_p.exp());
        (g.max(1.0e-300).ln() - ln_g_t, pair)
    };
    let mut lo = COND_P_MIN_PA.ln();
    let mut hi = COND_P_MAX_PA.ln();
    let (f_lo, _) = f(lo);
    let (f_hi, _) = f(hi);
    if f_lo > 0.0 || f_hi < 0.0 {
        return None;
    }
    // ln G grows about as ln p (rho_g ~ p^0.95), so a unit-slope Newton guess starts the secant.
    let mut x0 = p_guess_pa.clamp(COND_P_MIN_PA, COND_P_MAX_PA).ln();
    let (mut f0, mut pair0) = f(x0);
    let mut x1 = (x0 - f0).clamp(lo, hi);
    for _ in 0..60 {
        if f0.abs() < 1.0e-14 {
            return Some((x0.exp(), pair0));
        }
        if f0 < 0.0 {
            lo = lo.max(x0);
        } else {
            hi = hi.min(x0);
        }
        if (x1 - x0).abs() < 1.0e-15 || !(x1 > lo && x1 < hi) {
            x1 = 0.5 * (lo + hi);
        }
        let (f1, pair1) = f(x1);
        let slope = (f1 - f0) / (x1 - x0);
        let x2 = if slope.is_finite() && slope > 0.0 {
            x1 - f1 / slope
        } else {
            0.5 * (lo + hi)
        };
        x0 = x1;
        f0 = f1;
        pair0 = pair1;
        x1 = x2;
    }
    (f0.abs() < 1.0e-10).then(|| (x0.exp(), pair0))
}

impl SurfaceCondenserComponent {
    #[must_use]
    pub fn new() -> Self {
        let basis = water_if97_v1();
        Self {
            ports: [
                PortSpec::identified_stream(sid("steam_in"), PortRole::FlowIn, basis.clone())
                    .expect("condenser steam_in"),
                PortSpec::identified_stream(sid("condensate_out"), PortRole::EffortIn, basis)
                    .expect("condenser condensate_out"),
                PortSpec::heat(sid("cooling_heat"), PortRole::EffortIn)
                    .expect("condenser cooling_heat"),
                PortSpec::signal_in(sid("controls"), nz(3)),
            ],
            // Effective condenser + exhaust hood vapor space [m³] (gives smooth ~2-5 s vacuum dynamics)
            vapor_volume_m3: 85.0,
            tube_outer_diam_m: 0.0254,
            tube_inner_diam_m: 0.0230,
            tube_area_m2: 1_950.0,
            tube_count: 3_600.0,
            cw_nom_flow_kg_s: 2_600.0,
            ht_cache: Cell::new((
                (f64::NAN, f64::NAN, f64::NAN, f64::NAN),
                (0.0, 300.0, 300.0, 0.0, 0.0, 0.0),
            )),
        }
    }

    /// Evaluate shell-side film condensation (`fs_multiphase`) + tube-side Gnielinski (`fs_convection`)
    /// coupled across the condenser tube wall!
    fn evaluate_heat_transfer(
        &self,
        p_cond_pa: f64,
        t_sat_k: f64,
        t_cw_in_k: f64,
        cw_flow_frac: f64,
        air_mass_fraction: f64,
    ) -> (f64, f64, f64, f64, f64, f64) {
        let key = (p_cond_pa, t_cw_in_k, cw_flow_frac, air_mass_fraction);
        let cached = self.ht_cache.get();
        if cached.0 == key {
            return cached.1;
        }
        let lib = library();
        let cw_flow = (self.cw_nom_flow_kg_s * cw_flow_frac.clamp(0.15, 1.5)).max(100.0);
        let t_cw_ref = (t_cw_in_k + 4.0).clamp(283.15, 345.0);
        let rho_w = lib.cooling_water.density(t_cw_ref);
        let cp_w = lib.cooling_water.cp(t_cw_ref);
        let k_w = lib.cooling_water.conductivity(t_cw_ref);
        let mu_w = lib.cooling_water.dynamic_viscosity(t_cw_ref);

        // Cooling-water velocity inside 3,600 tubes (2-pass -> 1,800 parallel tubes per pass)
        let a_pass = (0.5 * self.tube_count) * 0.25 * PI * self.tube_inner_diam_m.powi(2);
        let vel_w = cw_flow / (rho_w * a_pass);
        let re_w = (rho_w * vel_w * self.tube_inner_diam_m / mu_w).clamp(3_500.0, 4.0e6);
        let pr_w = (cp_w * mu_w / k_w).clamp(0.6, 200.0);

        // Genuine fs_convection Gnielinski correlation!
        let h_cw = fs_convection::evaluate(
            CorrelationId::Gnielinski,
            CorrelationInputs::forced(re_w, pr_w),
        )
        .and_then(|nu| {
            nu.heat_transfer_coefficient(
                fs_convection::ThermalConductivity::new(k_w),
                Length::new(self.tube_inner_diam_m),
            )
        })
        .map(|ev| ev.value.value())
        .unwrap_or(7_500.0);

        // Wall + light fouling resistance [m²·K/W]
        let r_wall_foul = 8.5e-5;
        let u_cold = 1.0 / (1.0 / h_cw.max(200.0) + r_wall_foul);

        // Shell-side film condensation via genuine fs_multiphase::condensation::CondensationModelV1!
        let w_air = air_mass_fraction.clamp(0.0, 0.45);
        let geom = CondensationGeometryV1::HorizontalCylinder {
            diameter_m: self.tube_outer_diam_m,
        };
        let model_opt = if w_air < 1.0e-6 {
            CondensationModelV1::pure_steam(Pressure::new(p_cond_pa), geom).ok()
        } else {
            NonCondensableGasV1::steam_air(w_air, Temperature::new(t_sat_k + 0.5), 45.0)
                .ok()
                .and_then(|gas| {
                    CondensationModelV1::steam_air(Pressure::new(p_cond_pa), geom, gas).ok()
                })
        };

        let Some(cond_model) = model_opt else {
            return (0.0, t_cw_in_k, t_cw_in_k, 0.0, h_cw, 0.0);
        };
        let t_dew = cond_model.dew_temperature().value();
        if t_dew <= t_cw_in_k + 0.15 {
            return (0.0, t_cw_in_k, t_cw_in_k, 0.0, h_cw, 0.0);
        }

        // Solve coupled wall temperature T_wall in (t_cw_in, t_dew) matching shell flux q''_cond(T_wall)
        // with tube-side ε-NTU heat removal to cooling water!
        let c_cw = cw_flow * cp_w;
        let mut lo = t_cw_in_k + 0.05;
        let mut hi = t_dew - 0.05;
        let mut best_q = 0.0;
        let mut best_h_cond = 0.0;
        for _ in 0..26 {
            let tw = 0.5 * (lo + hi);
            let eval = match cond_model.evaluate(Temperature::new(tw)) {
                Ok(e) => e,
                Err(_) => break,
            };
            // Tube bundle inundation factor (~0.72 for multi-row condenser bank)
            let q_shell = 0.72 * self.tube_area_m2 * eval.heat_flux_w_m2;
            // Tube-side effectiveness for constant wall temperature tw
            let ntu = (u_cold * self.tube_area_m2 / c_cw).clamp(0.01, 8.0);
            let eff = 1.0 - (-ntu).exp();
            let q_water = c_cw * eff * (tw - t_cw_in_k);
            best_q = 0.5 * (q_shell + q_water);
            best_h_cond = 0.72 * eval.effective_htc_w_m2_k;
            if q_shell > q_water {
                lo = tw;
            } else {
                hi = tw;
            }
        }
        let t_wall = 0.5 * (lo + hi);
        let t_cw_out = t_cw_in_k + best_q / c_cw;
        let lmtd = (t_dew - 0.5 * (t_cw_in_k + t_cw_out)).max(0.5);
        let u_overall = best_q / (self.tube_area_m2 * lmtd);
        let res = (best_q, t_wall, t_cw_out, best_h_cond, h_cw, u_overall);
        self.ht_cache.set((key, res));
        res
    }
}

impl Default for SurfaceCondenserComponent {
    fn default() -> Self {
        Self::new()
    }
}

impl Component for SurfaceCondenserComponent {
    type State = CondenserState;
    type Error = TurbineComponentError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Smooth
    }

    fn stored_energy_j(&self, state: &Self::State) -> Result<f64, Self::Error> {
        Ok(state.total_u_j)
    }

    fn stored_mass_kg(&self, state: &Self::State) -> Result<Option<f64>, Self::Error> {
        Ok(Some(state.vapor_mass_kg))
    }

    fn stored_carried_quantity(
        &self,
        state: &Self::State,
        _carried: usize,
    ) -> Result<Option<f64>, Self::Error> {
        Ok(Some(state.vapor_mass_kg))
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
        check_cx(cx)?;
        // Inputs: steam_in (0..3), condensate_out [P_hotwell, h_hotwell, 1] (3..6),
        // cooling_heat [T_sink] (6), controls [t_cw_in_k, cw_flow_frac, air_mass_frac] (7..10)
        let m_in = inputs[STREAM_MASS_RATE];
        let h_dot_in = inputs[STREAM_ENTHALPY_RATE];
        let t_cw_in = inputs[7].clamp(278.0, 325.0);
        let cw_flow_frac = inputs[8].clamp(0.15, 1.5);
        let w_air = inputs[9].clamp(0.0, 0.45);

        // The shell is a saturated control volume: vapour of mass M = V rho_g(p) over its own
        // condensate. Mass:   dM/dt = m_in - m_out
        //             Energy: dU/dt = H_in - m_out h_f - q,   U = M u_g(p).
        // Eliminating the condensate flow m_out gives one equation for p,
        //   d/dt [V rho_g (u_g - h_f)] = H_in - m_in h_f - q,
        // so wet exhaust (h_in < h_g) needs no separate rule: its liquid simply leaves as
        // condensate, and the duty settles at m_in (h_in - h_f).
        let pair_old = If97SaturationPair::at_pressure(Pressure::new(
            old.p_cond_pa.clamp(COND_P_MIN_PA, COND_P_MAX_PA),
        ))
        .expect("condenser old sat pair");
        let h_f = pair_old.liquid().specific_enthalpy().value();
        let s_f = pair_old.liquid().specific_entropy().value();
        let h_g = pair_old.vapor().specific_enthalpy().value();
        let t_sat_old = pair_old.vapor().temperature().value();

        let (q_cond, t_wall, t_cw_out, h_cond, h_cw, u_overall) =
            self.evaluate_heat_transfer(old.p_cond_pa, t_sat_old, t_cw_in, cw_flow_frac, w_air);

        let g_old = old.total_u_j - h_f * old.vapor_mass_kg;
        let g_target = g_old + dt_s * (h_dot_in - m_in * h_f - q_cond);
        let Some((p_new, pair_new)) =
            solve_condenser_pressure(self.vapor_volume_m3, h_f, g_target, old.p_cond_pa)
        else {
            return Ok(Proposal::Shorten { dt_s: 0.5 * dt_s });
        };
        let vap_new = pair_new.vapor();
        let m_new = self.vapor_volume_m3 / vap_new.specific_volume().value();
        let u_new = m_new * vap_new.specific_internal_energy().value();
        let t_sat_new = vap_new.temperature().value();
        let h_g_new = vap_new.specific_enthalpy().value();
        // Condensate leaving to the hotwell [kg/s] (negative while the shell flashes liquid back).
        let m_out = m_in - (m_new - old.vapor_mass_kg) / dt_s;

        // Entropy produced in the shell, from its exact balance
        //   S_gen = dS - dt (m_in s_in - m_out s_f - q / T_cw),   S = M s_g(p),
        // with the heat leaving at the cooling water's mean temperature.
        let p_mid = 0.5 * (old.p_cond_pa + p_new);
        let s_in = if m_in.abs() > 1.0e-12 {
            eval_ph(p_mid, h_dot_in / m_in).s_j_kg_k
        } else {
            0.0
        };
        let s_old = old.vapor_mass_kg * pair_old.vapor().specific_entropy().value();
        let s_new = m_new * vap_new.specific_entropy().value();
        let t_cw_mean = 0.5 * (t_cw_in + t_cw_out);
        let s_gen_w_k = (s_new - s_old) / dt_s - m_in * s_in + m_out * s_f + q_cond / t_cw_mean;
        // (The heat route still delivers q to the cooling-water sink at its fixed route
        // temperature, input 6; production between the water and that sink is the route's.)

        *new = CondenserState {
            vapor_mass_kg: m_new,
            total_u_j: u_new,
            p_cond_pa: p_new,
            t_sat_k: t_sat_new,
            t_wall_k: t_wall,
            t_cw_out_k: t_cw_out,
            q_cond_w: q_cond,
            h_cond_w_m2_k: h_cond,
            h_cw_w_m2_k: h_cw,
            u_overall_w_m2_k: u_overall,
            m_condensate_kg_s: m_out,
            air_mass_fraction: w_air,
            s_gen_w_k,
        };

        // Outputs:
        // 0..3 steam_in (FlowIn -> returns [P_mid, h_g_mid, 1.0])
        outputs[STREAM_PRESSURE] = p_mid;
        outputs[STREAM_SPECIFIC_ENTHALPY] = 0.5 * (h_g + h_g_new);
        outputs[STREAM_CARRIED] = 1.0;

        // 3..6 condensate_out (EffortIn -> returns [-m_out, -m_out*h_f, -m_out], leaving condenser!)
        outputs[3 + STREAM_MASS_RATE] = -m_out;
        outputs[3 + STREAM_ENTHALPY_RATE] = -m_out * h_f;
        outputs[3 + STREAM_CARRIED] = -m_out;

        // 6 cooling_heat (EffortIn -> returns -q_cond leaving condenser into cooling-water sink!)
        outputs[6] = -q_cond;

        let du = u_new - old.total_u_j;
        let e_tol = 1.0e-7 * (u_new.abs() + old.total_u_j.abs() + dt_s * q_cond.abs() + 1.0);
        Ok(Proposal::Step(Receipt {
            stored_change_j: du,
            dissipated_j: 0.0,
            numerical_loss_j: 0.0,
            entropy_generated_j_per_k: dt_s * s_gen_w_k,
            energy_tolerance_j: e_tol,
            residual: 0.0,
            tolerance: 0.0,
        }))
    }
}

// ---------------------------------------------------------------------------
// 5. RotorShaftComponent (Exact Geometry-Derived Polar Inertia J_p)
// ---------------------------------------------------------------------------

/// Rotational shaft inertia component whose polar moment `J_p` comes directly from
/// `geometry::build_turbine_geometry().rotor_props.polar_inertia_kg_m2`.
///
/// Uses exact midpoint integration so discrete kinetic energy change equals
/// `dt * tau_net * omega_mid` to machine precision (zero numerical energy drift!).
pub struct RotorShaftComponent {
    ports: [PortSpec; 1],
    /// Combined polar moment of inertia (forged rotor + 17 NURBS blade rows + generator rotor) [kg·m²].
    pub polar_inertia_kg_m2: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct RotorShaftState {
    pub omega_rad_s: f64,
    pub theta_rad: f64,
    pub net_torque_nm: f64,
}

impl RotorShaftComponent {
    #[must_use]
    pub fn new(polar_inertia_kg_m2: f64) -> Self {
        Self {
            ports: [PortSpec::power(
                sid("shaft"),
                PortKind::RotationalTorqueAngularVelocity,
                PortRole::EffortIn,
                nz(1),
            )
            .expect("rotor shaft port")],
            polar_inertia_kg_m2,
        }
    }
}

impl Component for RotorShaftComponent {
    type State = RotorShaftState;
    type Error = TurbineComponentError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn stored_energy_j(&self, state: &Self::State) -> Result<f64, Self::Error> {
        Ok(0.5 * self.polar_inertia_kg_m2 * state.omega_rad_s * state.omega_rad_s)
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
        check_cx(cx)?;
        let tau_net = inputs[0];
        let d_omega = dt_s * tau_net / self.polar_inertia_kg_m2;
        let omega_new = old.omega_rad_s + d_omega;
        let omega_mid = 0.5 * (old.omega_rad_s + omega_new);
        let theta_new = (old.theta_rad + dt_s * omega_mid).rem_euclid(2.0 * PI);

        *new = RotorShaftState {
            omega_rad_s: omega_new,
            theta_rad: theta_new,
            net_torque_nm: tau_net,
        };
        outputs[0] = omega_mid;

        // Exact discrete work identity: 0.5 J (w1^2 - w0^2) = J * w_mid * (w1 - w0) = dt * tau * w_mid
        let de = dt_s * tau_net * omega_mid;
        let e_scale = (0.5 * self.polar_inertia_kg_m2 * (omega_new * omega_new + old.omega_rad_s * old.omega_rad_s))
            .max(1.0);
        Ok(Proposal::Step(Receipt::lossless(de, 1.0e-10 * e_scale + 1.0e-6)))
    }
}

impl PortTangent for RotorShaftComponent {
    fn jvp(
        &self,
        _old: &Self::State,
        dt_s: f64,
        _inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), Self::Error> {
        check_cx(cx)?;
        d_outputs[0] = 0.5 * dt_s * d_inputs[0] / self.polar_inertia_kg_m2;
        Ok(())
    }

    fn vjp(
        &self,
        _old: &Self::State,
        dt_s: f64,
        _inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), Self::Error> {
        check_cx(cx)?;
        bar_inputs[0] = 0.5 * dt_s * bar_outputs[0] / self.polar_inertia_kg_m2;
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// 6. JournalBearingComponent (Front Radial/Thrust + Rear Radial Bearings)
// ---------------------------------------------------------------------------

/// Hydrodynamic journal & Kingsbury thrust bearing component computing real
/// `fs_tribo::journal::JournalGeometry` Couette shear drag + Ocvirk Sommerfeld
/// film thickness & 8-coefficient `(K_ij, C_ij)` + active Stribeck labyrinth rub torque.
pub struct JournalBearingComponent {
    ports: [PortSpec; 3],
    pub front_static_load_n: f64,
    pub rear_static_load_n: f64,
    pub clearance_m: f64,
    pub preload: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct BearingSystemState {
    pub front: BearingDynamicState,
    pub rear: BearingDynamicState,
    pub thrust_loss_w: f64,
    pub rub_loss_w: f64,
    pub total_drag_torque_nm: f64,
    pub total_heat_w: f64,
}

impl JournalBearingComponent {
    #[must_use]
    pub fn new(front_static_load_n: f64, rear_static_load_n: f64) -> Self {
        Self::with_design(
            front_static_load_n,
            rear_static_load_n,
            &crate::stage_meanline::TurbineDesignParams::default(),
        )
    }

    #[must_use]
    pub fn with_design(
        front_static_load_n: f64,
        rear_static_load_n: f64,
        d: &crate::stage_meanline::TurbineDesignParams,
    ) -> Self {
        let d = d.clamped();
        Self {
            ports: [
                PortSpec::power(
                    sid("shaft"),
                    PortKind::RotationalTorqueAngularVelocity,
                    PortRole::FlowIn,
                    nz(1),
                )
                .expect("bearing shaft port"),
                PortSpec::heat(sid("oil_heat"), PortRole::EffortIn).expect("bearing heat port"),
                PortSpec::signal_in(sid("faults"), nz(4)),
            ],
            front_static_load_n,
            rear_static_load_n,
            clearance_m: d.bearing_clearance_um * 1.0e-6,
            preload: d.bearing_preload,
        }
    }

    #[must_use]
    pub fn initial_state(&self, omega_rad_s: f64) -> BearingSystemState {
        let eff_destab = (0.32 - self.preload).max(0.0) * 1.1;
        let front = evaluate_journal_bearing(
            0.128,
            0.24,
            self.clearance_m,
            self.front_static_load_n,
            omega_rad_s,
            333.15,
            1.0,
            eff_destab,
        );
        let rear = evaluate_journal_bearing(
            0.132,
            0.26,
            self.clearance_m,
            self.rear_static_load_n,
            omega_rad_s,
            335.15,
            1.0,
            eff_destab,
        );
        BearingSystemState {
            front,
            rear,
            thrust_loss_w: 18_000.0,
            rub_loss_w: 0.0,
            total_drag_torque_nm: front.drag_torque_nm + rear.drag_torque_nm + 58.0,
            total_heat_w: front.shear_loss_w + rear.shear_loss_w + 18_000.0,
        }
    }
}

impl Component for JournalBearingComponent {
    type State = BearingSystemState;
    type Error = TurbineComponentError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Smooth
    }

    fn stored_energy_j(&self, _state: &Self::State) -> Result<f64, Self::Error> {
        Ok(0.0)
    }

    fn propose(
        &self,
        _old: &Self::State,
        dt_s: f64,
        inputs: &[f64],
        new: &mut Self::State,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, Self::Error> {
        check_cx(cx)?;
        // Inputs: shaft omega_in = -omega_shaft (0), oil_heat T_oil (1),
        // faults [oil_pressure_frac, oil_whip_destab, rub_torque_nm, axial_thrust_n] (2..6)
        let omega_in = inputs[0];
        let omega_shaft = -omega_in;
        let t_oil = inputs[1].clamp(285.0, 520.0);
        let oil_p_frac = inputs[2].clamp(0.15, 1.5);
        let oil_whip_destab = inputs[3].clamp(0.0, 2.5);
        let rub_torque_nm = inputs[4].max(0.0);
        let axial_thrust_n = inputs[5].abs();

        let eff_destab = (oil_whip_destab + (0.32 - self.preload).max(0.0) * 1.1).clamp(0.0, 2.5);
        let front = evaluate_journal_bearing(
            0.128,
            0.24,
            self.clearance_m,
            self.front_static_load_n,
            omega_shaft,
            t_oil,
            oil_p_frac,
            eff_destab,
        );
        let rear = evaluate_journal_bearing(
            0.132,
            0.26,
            self.clearance_m,
            self.rear_static_load_n,
            omega_shaft,
            t_oil + 2.0,
            oil_p_frac,
            eff_destab,
        );

        // Kingsbury tilting-pad thrust collar shear loss (r_in = 0.128, r_out = 0.214 from rotor chart!)
        let mu = front.viscosity_pa_s;
        let h_thrust = (65.0e-6 / (1.0 + axial_thrust_n / 1.5e5).sqrt()).clamp(14.0e-6, 90.0e-6);
        let thrust_torque =
            (0.5 * PI * mu * omega_shaft.abs() * (0.214_f64.powi(4) - 0.128_f64.powi(4)) / h_thrust)
                * (1.0 / oil_p_frac.clamp(0.3, 1.2));
        let thrust_loss_w = thrust_torque * omega_shaft.abs();
        let rub_loss_w = rub_torque_nm * omega_shaft.abs();

        let preload_shear_factor = 1.0 + 0.35 * self.preload;
        let total_drag_mag = (front.drag_torque_nm + rear.drag_torque_nm) * preload_shear_factor
            + thrust_torque
            + rub_torque_nm;
        // Braking torque opposes shaft rotation: when omega_shaft > 0, tau_on_shaft < 0
        let tau_on_shaft = -total_drag_mag * omega_shaft.signum();
        // Mechanical power entering bearing = omega_in * tau_on_shaft = (-omega_shaft) * (-drag * sgn) >= 0
        let p_heat = (omega_in * tau_on_shaft).max(0.0);

        outputs[0] = tau_on_shaft;
        outputs[1] = -p_heat; // Exact dissipated mechanical power exits via heat port!

        *new = BearingSystemState {
            front,
            rear,
            thrust_loss_w,
            rub_loss_w,
            total_drag_torque_nm: total_drag_mag,
            total_heat_w: p_heat,
        };

        let s_gen = dt_s * p_heat / t_oil;
        Ok(Proposal::Step(Receipt {
            stored_change_j: 0.0,
            dissipated_j: dt_s * p_heat,
            numerical_loss_j: 0.0,
            entropy_generated_j_per_k: s_gen,
            energy_tolerance_j: 1.0e-8 * (dt_s * p_heat + 1.0) + 1.0e-6,
            residual: 0.0,
            tolerance: 0.0,
        }))
    }
}

// ---------------------------------------------------------------------------
// 7. SynchronousGeneratorComponent (fs-phs SynchronousMotor + M19 Steinmetz + IACS Cu)
// ---------------------------------------------------------------------------

/// 2-pole, 50 Hz, 11 kV / 60 MVA synchronous turbo-generator on an infinite bus, in the
/// classical stability model: a constant EMF `E'` behind the transient reactance `x'_d` plus
/// the step-up transformer's `x_t`, with stator transients neglected (quasi-steady d-q).
///
/// The machine is an [`fs_phs::circuit::synchronous::SynchronousMotor`] whose field linkage is
/// `psi = E' / omega_s` and whose (round-rotor) inductance is `(x'_d + x_t) / omega_s`. Each step
/// solves the stator voltage equations against the bus phasor at the load angle `delta`,
///   `v_d = R i_d - X i_q`,  `v_q = R i_q + X i_d + omega psi`,  `(v_d, v_q) = V (sin delta, cos delta)`,
/// builds the machine state from the resulting flux linkages and takes the electromagnetic torque,
/// currents and copper loss from `SynchronousMotor::operating_point`. Reactive power, line
/// current and pull-out follow from the same equations; nothing is a fitted torque map. The
/// amortisseur damping is the one lumped machine coefficient. A pole slip (`|delta| > pi`)
/// trips the machine off the bus (out-of-step protection) until the breaker is opened and
/// closed again.
pub struct SynchronousGeneratorComponent {
    ports: [PortSpec; 4],
    /// fs-phs machine for the last field linkage used, rebuilt when the excitation changes.
    phs_machine: RefCell<Option<(f64, SynchronousMotor)>>,
    pub sync_omega_rad_s: f64,
    pub nominal_voltage_v: f64,
    pub rated_va: f64,
    /// Transient d-axis reactance [pu] (turbo-generators: 0.2..0.3).
    pub x_d_transient_pu: f64,
    /// Step-up transformer leakage reactance [pu].
    pub x_transformer_pu: f64,
    /// Armature resistance at `T_REF_K` [pu].
    pub r_armature_pu: f64,
    pub stator_core_mass_kg: f64,
    pub amortisseur_damping_nm_s_per_rad: f64,
}

/// Electrical solution of the generator at one load angle.
#[derive(Debug, Clone, Copy)]
struct GeneratorSolution {
    /// Electromagnetic torque ON the rotor [N m] (negative while generating).
    torque_nm: f64,
    current_dq_a: [f64; 2],
    flux_dq_wb: [f64; 2],
    /// Reactive power delivered to the bus [VAr].
    q_out_var: f64,
    copper_loss_w: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct GeneratorState {
    /// fs-phs machine state at the operating point, rotor-aligned (`theta = 0`):
    /// `[lambda_d, lambda_q, 0, p_rotor, p_stator]`.
    pub phs_x: [f64; 5],
    /// Power/load angle `δ` relative to the 50 Hz infinite bus [rad].
    pub load_angle_rad: f64,
    /// Active electrical power exported to grid [W].
    pub elec_power_w: f64,
    /// Reactive power [VAr].
    pub reactive_power_var: f64,
    /// Stator line current RMS [A].
    pub line_current_a: f64,
    /// Copper I²R + M19 Steinmetz core loss rejected as heat [W].
    pub total_loss_w: f64,
    pub core_loss_w: f64,
    pub copper_loss_w: f64,
    pub breaker_closed: bool,
    pub em_torque_nm: f64,
    /// Latched by a pole slip (`|delta| > pi`); cleared when the breaker is opened.
    pub out_of_step_tripped: bool,
}

impl SynchronousGeneratorComponent {
    #[must_use]
    pub fn new() -> Self {
        Self {
            ports: [
                PortSpec::power(
                    sid("shaft"),
                    PortKind::RotationalTorqueAngularVelocity,
                    PortRole::FlowIn,
                    nz(1),
                )
                .expect("gen shaft port"),
                PortSpec::power(
                    sid("grid_elec"),
                    PortKind::ElectricalVoltageCurrent,
                    PortRole::EffortIn,
                    nz(1),
                )
                .expect("gen grid port"),
                PortSpec::heat(sid("cooler_heat"), PortRole::EffortIn).expect("gen cooler port"),
                PortSpec::signal_in(sid("controls"), nz(3)),
            ],
            phs_machine: RefCell::new(None),
            sync_omega_rad_s: 100.0 * PI, // 3,000 rpm = 50 Hz 2-pole
            nominal_voltage_v: 11_000.0,
            rated_va: 60.0e6,
            x_d_transient_pu: 0.25,
            x_transformer_pu: 0.10,
            r_armature_pu: 0.003,
            stator_core_mass_kg: 14_500.0,
            // Damper-winding torque coefficient D = 15 pu (P_pu per unit speed deviation; machines
            // with amortisseurs: about 2..20), i.e. D_pu * S / omega_s^2. It was 1.4e5 N m s/rad
            // (230 pu), which over-damped the electromechanical mode (damping ratio 1.85).
            amortisseur_damping_nm_s_per_rad: 15.0 * 60.0e6 / (100.0 * PI * 100.0 * PI),
        }
    }

    /// Excitation range the field system accepts [pu].
    const EXCITATION_RANGE_PU: (f64, f64) = (0.4, 1.5);

    /// M19 silicon-steel stator core loss [W] at shaft speed `omega_shaft` and field
    /// `excitation_pu` (Steinmetz parameters from matdb).
    fn core_loss_w(&self, omega_shaft: f64, excitation_pu: f64) -> f64 {
        let f_elec = omega_shaft.abs() / (2.0 * PI);
        let b_peak_t = (1.35 * excitation_pu).clamp(0.2, 1.75);
        self.stator_core_mass_kg * library().generator.m19_specific_core_loss(f_elec, b_peak_t)
    }

    fn base_impedance_ohm(&self) -> f64 {
        self.nominal_voltage_v * self.nominal_voltage_v / self.rated_va
    }

    /// Series reactance from the EMF to the bus at grid frequency [ohm].
    #[must_use]
    pub fn reactance_ohm(&self) -> f64 {
        (self.x_d_transient_pu + self.x_transformer_pu) * self.base_impedance_ohm()
    }

    /// Field linkage giving an EMF of `excitation_pu` times the nominal bus voltage at
    /// synchronous speed [Wb, power-invariant].
    #[must_use]
    pub fn field_linkage_wb(&self, excitation_pu: f64) -> f64 {
        excitation_pu * self.nominal_voltage_v / self.sync_omega_rad_s
    }

    fn with_machine<R>(
        &self,
        psi_wb: f64,
        f: impl FnOnce(&SynchronousMotor) -> R,
    ) -> Result<R, TurbineComponentError> {
        let mut slot = self.phs_machine.borrow_mut();
        let stale = !matches!(&*slot, Some((psi, _)) if *psi == psi_wb);
        if stale {
            let l = self.reactance_ohm() / self.sync_omega_rad_s;
            let params = SynchronousMotorParameters {
                d_inductance_h: l,
                q_inductance_h: l,
                winding_resistance_ohm: self.r_armature_pu * self.base_impedance_ohm(),
                pm_flux_linkage_wb: psi_wb,
                pole_pairs: 1,
                rotor_inertia_kg_m2: 180.0,
                stator_inertia_kg_m2: 1.0e7,
                bearing_drag_nm_s_per_rad: 0.0,
            };
            let motor = SynchronousMotor::new(params).map_err(|_| {
                TurbineComponentError::InvalidState("fs-phs generator parameters refused")
            })?;
            *slot = Some((psi_wb, motor));
        }
        Ok(f(&slot.as_ref().expect("machine just built").1))
    }

    /// Quasi-steady stator solution at load angle `delta_rad`, shaft speed `omega_shaft`, bus
    /// voltage `v_bus` and winding resistance `r_ohm`.
    fn solve_electrical(
        &self,
        delta_rad: f64,
        omega_shaft: f64,
        v_bus: f64,
        excitation_pu: f64,
        r_ohm: f64,
    ) -> Result<GeneratorSolution, TurbineComponentError> {
        let x = self.reactance_ohm();
        let l = x / self.sync_omega_rad_s;
        let psi = self.field_linkage_wb(excitation_pu);
        let e = omega_shaft * psi;
        let (sd, cd) = delta_rad.sin_cos();
        let (v_d, v_q) = (v_bus * sd, v_bus * cd);
        let det = r_ohm * r_ohm + x * x;
        let i_d = (r_ohm * v_d + x * (v_q - e)) / det;
        let i_q = (r_ohm * (v_q - e) - x * v_d) / det;
        let flux = [psi + l * i_d, l * i_q];
        // Machine state at the rotor-aligned frame (theta = 0): [lambda_alpha, lambda_beta,
        // theta, p_rotor, p_stator].
        let state = [flux[0], flux[1], 0.0, 180.0 * omega_shaft, 0.0];
        let op = self
            .with_machine(psi, |m| m.operating_point(&state))?
            .map_err(|_| {
                TurbineComponentError::InvalidState("fs-phs generator operating point refused")
            })?;
        let [id, iq] = op.current_dq_a;
        Ok(GeneratorSolution {
            torque_nm: op.electromagnetic_torque_nm,
            current_dq_a: op.current_dq_a,
            flux_dq_wb: op.flux_dq_wb,
            q_out_var: v_d * iq - v_q * id,
            copper_loss_w: r_ohm * (id * id + iq * iq),
        })
    }

    /// Steady operating state exporting `initial_power_w` at `excitation_pu` on the nominal bus.
    #[must_use]
    pub fn initial_state(
        &self,
        omega_rad_s: f64,
        initial_power_w: f64,
        excitation_pu: f64,
    ) -> GeneratorState {
        let (exc_lo, exc_hi) = Self::EXCITATION_RANGE_PU;
        let excitation_pu = excitation_pu.clamp(exc_lo, exc_hi);
        let v = self.nominal_voltage_v;
        let e = excitation_pu * v;
        let delta = (initial_power_w * self.reactance_ohm() / (e * v))
            .clamp(-0.95, 0.95)
            .asin();
        let r = self.r_armature_pu * self.base_impedance_ohm();
        let sol = self
            .solve_electrical(delta, self.sync_omega_rad_s, v, excitation_pu, r)
            .expect("generator initial operating point");
        let [id, iq] = sol.current_dq_a;
        let core_loss_w = self.core_loss_w(omega_rad_s, excitation_pu);
        GeneratorState {
            phs_x: [
                sol.flux_dq_wb[0],
                sol.flux_dq_wb[1],
                0.0,
                180.0 * omega_rad_s,
                0.0,
            ],
            load_angle_rad: delta,
            elec_power_w: initial_power_w,
            reactive_power_var: sol.q_out_var,
            line_current_a: id.hypot(iq) / 3.0_f64.sqrt(),
            total_loss_w: sol.copper_loss_w + core_loss_w,
            core_loss_w,
            copper_loss_w: sol.copper_loss_w,
            breaker_closed: true,
            em_torque_nm: sol.torque_nm,
            out_of_step_tripped: false,
        }
    }
}

impl Default for SynchronousGeneratorComponent {
    fn default() -> Self {
        Self::new()
    }
}

impl Component for SynchronousGeneratorComponent {
    type State = GeneratorState;
    type Error = TurbineComponentError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Smooth
    }

    fn stored_energy_j(&self, _state: &Self::State) -> Result<f64, Self::Error> {
        Ok(0.0)
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
        check_cx(cx)?;
        // Inputs: shaft omega_in = -omega_shaft (0), grid_elec V_bus (1), cooler T_cool (2),
        // controls [breaker_closed (0/1), excitation_pu, grid_freq_hz] (3..6)
        let omega_in = inputs[0];
        let omega_shaft = -omega_in;
        let v_bus = inputs[1].max(100.0);
        let t_cool = inputs[2].clamp(285.0, 500.0);
        let breaker_closed = inputs[3] >= 0.5;
        let (exc_lo, exc_hi) = Self::EXCITATION_RANGE_PU;
        let excitation_pu = inputs[4].clamp(exc_lo, exc_hi);
        let grid_omega = 2.0 * PI * inputs[5].clamp(45.0, 65.0);

        let lib = library();
        let core_loss_w = self.core_loss_w(omega_shaft, excitation_pu);

        // Out-of-step protection: a pole slip latches the machine off the bus until the breaker
        // is opened by the operator.
        let tripped = breaker_closed && old.out_of_step_tripped;
        let synced = breaker_closed && !tripped && omega_shaft > 20.0;
        let r_cu = self.r_armature_pu
            * self.base_impedance_ohm()
            * (lib.generator.copper_resistivity(t_cool + 35.0)
                / lib.generator.copper_resistivity(T_REF_K));
        let (delta_new, tau_on_shaft, q_var, i_rms, cu_loss_w, damper_heat_w, phs_next, slipped) =
            if synced {
                // Load angle against the infinite bus: d(delta)/dt = omega_shaft - omega_grid.
                let delta_new = old.load_angle_rad + dt_s * (omega_shaft - grid_omega);
                let delta_mid = 0.5 * (old.load_angle_rad + delta_new);
                let sol =
                    self.solve_electrical(delta_mid, omega_shaft, v_bus, excitation_pu, r_cu)?;
                // Amortisseur (damper-winding) torque on the slip; its slip power is heat.
                let slip = omega_shaft - grid_omega;
                let tau_damper = -self.amortisseur_damping_nm_s_per_rad * slip;
                let tau_core = -core_loss_w / omega_shaft;
                let [id, iq] = sol.current_dq_a;
                (
                    delta_new,
                    sol.torque_nm + tau_damper + tau_core,
                    sol.q_out_var,
                    id.hypot(iq) / 3.0_f64.sqrt(),
                    sol.copper_loss_w,
                    self.amortisseur_damping_nm_s_per_rad * slip * slip,
                    [
                        sol.flux_dq_wb[0],
                        sol.flux_dq_wb[1],
                        0.0,
                        180.0 * omega_shaft,
                        0.0,
                    ],
                    delta_new.abs() > PI,
                )
            } else {
                // Breaker open, tripped or standstill: no stator current, only open-circuit core drag.
                let tau_core = if omega_shaft.abs() > 1.0 {
                    -core_loss_w / omega_shaft
                } else {
                    0.0
                };
                let open_circuit = [
                    self.field_linkage_wb(excitation_pu),
                    0.0,
                    0.0,
                    180.0 * omega_shaft,
                    0.0,
                ];
                (0.0, tau_core, 0.0, 0.0, 0.0, 0.0, open_circuit, false)
            };

        // Exact power balance: shaft power in = grid export + heat to the cooler.
        let p_shaft_in = omega_in * tau_on_shaft;
        let p_heat_out = if synced {
            cu_loss_w + core_loss_w + damper_heat_w
        } else {
            p_shaft_in.max(0.0)
        };
        let p_grid_export = p_shaft_in - p_heat_out; // Positive when generating into the grid
        let i_grid_into_gen = -p_grid_export / v_bus;

        outputs[0] = tau_on_shaft;
        outputs[1] = i_grid_into_gen;
        outputs[2] = -p_heat_out;

        let out_of_step = breaker_closed && (tripped || slipped);
        *new = GeneratorState {
            phs_x: phs_next,
            load_angle_rad: if out_of_step { 0.0 } else { delta_new },
            elec_power_w: p_grid_export,
            reactive_power_var: q_var,
            line_current_a: i_rms,
            total_loss_w: p_heat_out,
            core_loss_w,
            copper_loss_w: cu_loss_w,
            breaker_closed: breaker_closed && !out_of_step,
            em_torque_nm: tau_on_shaft,
            out_of_step_tripped: out_of_step,
        };

        let s_gen = dt_s * p_heat_out / t_cool;
        let e_scale = (dt_s * (p_shaft_in.abs() + p_grid_export.abs() + p_heat_out.abs())).max(1.0);
        Ok(Proposal::Step(Receipt {
            stored_change_j: 0.0,
            dissipated_j: dt_s * p_heat_out,
            numerical_loss_j: 0.0,
            entropy_generated_j_per_k: s_gen,
            energy_tolerance_j: 1.0e-8 * e_scale + 1.0e-4,
            residual: 0.0,
            tolerance: 0.0,
        }))
    }
}

// ---------------------------------------------------------------------------
// 8. TurbineCasingThermalComponent (8-Node Metal/Oil Network + NortonCreep)
// ---------------------------------------------------------------------------

/// Thermal node indices inside the 8-node turbine metal & oil thermal network.
pub const TH_HP_ROTOR: usize = 0;
pub const TH_HP_CASING: usize = 1;
pub const TH_IP_ROTOR: usize = 2;
pub const TH_IP_CASING: usize = 3;
pub const TH_LP_HOOD: usize = 4;
pub const TH_FRONT_BEARING: usize = 5;
pub const TH_REAR_BEARING: usize = 6;
pub const TH_LUBE_OIL: usize = 7;

/// Multi-node metal and lube-oil thermal mass network of the model turbine.
///
/// Computes exact heat storage across 8 physical thermal capacitances (derived from
/// `geometry::build_turbine_geometry()` casing/rotor masses and matdb `cp`), ambient
/// free convection via `fs_convection::evaluate(CorrelationId::ChurchillChuHorizontalCylinder)`,
/// rotor-vs-casing differential axial expansion `ΔL_diff`, casing thermal hoop stress, and
/// live high-temperature CrMoV creep rate via `fs_material::creep::NortonCreep`.
pub struct TurbineCasingThermalComponent {
    ports: [PortSpec; 5],
    /// Thermal heat capacities `C_i = m_i * cp_i` [J/K] for the 8 nodes.
    pub caps_j_k: [f64; 8],
    /// External casing area for ambient convection + radiation [m²].
    pub casing_ext_area_m2: f64,
    /// Lube-oil cooler UA conductance [W/K] at nominal oil flow.
    pub oil_cooler_ua_w_k: f64,
    /// Share of the HP steam-path heat that enters the HP rotor; the rest enters the HP
    /// casing. From the wetted areas of the flow path (`r_hub L` against `r_tip L`).
    pub hp_rotor_share: f64,
    /// Shares of the IP/LP steam-path heat entering the IP/LP rotor, the IP casing and the
    /// LP hood, from the same wetted areas.
    pub iplp_shares: [f64; 3],
}

/// Flow-path section of the blade rows: axial length and mean hub/tip radii [m].
struct FlowSection {
    length_m: f64,
    r_hub_m: f64,
    r_tip_m: f64,
}

fn flow_section(
    rows: &[crate::geometry::BladeRowSpec],
    stages: core::ops::RangeInclusive<usize>,
) -> FlowSection {
    let sec: Vec<_> = rows
        .iter()
        .filter(|r| stages.contains(&r.stage_index))
        .collect();
    let z_min = sec
        .iter()
        .map(|r| r.z_center_m - 0.5 * r.axial_chord_m)
        .fold(f64::INFINITY, f64::min);
    let z_max = sec
        .iter()
        .map(|r| r.z_center_m + 0.5 * r.axial_chord_m)
        .fold(f64::NEG_INFINITY, f64::max);
    let n = sec.len().max(1) as f64;
    FlowSection {
        length_m: (z_max - z_min).max(0.0),
        r_hub_m: sec.iter().map(|r| r.r_hub_m).sum::<f64>() / n,
        r_tip_m: sec.iter().map(|r| r.r_tip_m).sum::<f64>() / n,
    }
}

#[derive(Debug, Clone, Copy)]
pub struct ThermalNetworkState {
    /// Temperatures of all 8 nodes [K].
    pub temps_k: [f64; 8],
    /// Rotor-minus-casing differential axial thermal expansion [mm] (+ = rotor expanding faster).
    pub diff_expansion_mm: f64,
    /// Peak HP casing thermal bending/hoop stress [MPa].
    pub hp_thermal_stress_mpa: f64,
    /// Instantaneous CrMoV Norton creep strain rate [1/s] from `fs_material::creep`.
    pub hp_creep_rate_per_s: f64,
    /// Cumulative HP rotor creep strain [%].
    pub cumulative_creep_pct: f64,
    /// Ambient + oil-cooler total heat rejected [W].
    pub heat_rejected_w: f64,
}

impl TurbineCasingThermalComponent {
    /// Lumped thermal network whose node capacities are the full `m cp` of the parts.
    ///
    /// The rotor forging is split between the HP and IP/LP rotor nodes in proportion to
    /// `r_hub^2 L` of each flow-path section (a solid forging), and the HP/IP casing between
    /// the HP and IP casing nodes in proportion to `r_tip L` (a shell around the flow path).
    /// No capacity is scaled down to make transients quicker: thermal time constants are
    /// what the masses give, and the page offers faster simulation rates instead.
    #[must_use]
    pub fn new(
        rotor_mass_kg: f64,
        hp_casing_mass_kg: f64,
        lp_hood_mass_kg: f64,
        rows: &[crate::geometry::BladeRowSpec],
    ) -> Self {
        let lib = library();
        let cp_r = lib.rotor_steel.cp_j_kg_k;
        let cp_c = lib.hp_casing_steel.cp_j_kg_k;
        let cp_lp = lib.lp_exhaust_steel.cp_j_kg_k;
        let cp_brg = lib.bearing_bronze.cp_j_kg_k;
        let cp_oil = lib.lube_oil.cp(333.15);

        let hp = flow_section(rows, 0..=5);
        let ip = flow_section(rows, 6..=13);
        let lp = flow_section(rows, 14..=16);
        let w_rotor = |s: &FlowSection| s.r_hub_m * s.r_hub_m * s.length_m;
        let w_shell = |s: &FlowSection| s.r_tip_m * s.length_m;
        let rotor_hp_frac = w_rotor(&hp) / (w_rotor(&hp) + w_rotor(&ip) + w_rotor(&lp));
        let casing_hp_frac = w_shell(&hp) / (w_shell(&hp) + w_shell(&ip));

        let caps_j_k = [
            rotor_hp_frac * rotor_mass_kg * cp_r,              // 0: HP rotor
            casing_hp_frac * hp_casing_mass_kg * cp_c,         // 1: HP casing
            (1.0 - rotor_hp_frac) * rotor_mass_kg * cp_r,      // 2: IP/LP rotor
            (1.0 - casing_hp_frac) * hp_casing_mass_kg * cp_c, // 3: IP casing
            lp_hood_mass_kg * cp_lp,                           // 4: LP exhaust hood
            65.0 * cp_brg,                                     // 5: Front Babbitt/bronze pads
            70.0 * cp_brg,                                     // 6: Rear Babbitt/bronze pads
            420.0 * cp_oil,                                    // 7: Lube-oil inventory
        ];

        // Steam-side heat splits by wetted area of each section.
        let hp_rotor_share = hp.r_hub_m / (hp.r_hub_m + hp.r_tip_m);
        let a_rotor = ip.r_hub_m * ip.length_m + lp.r_hub_m * lp.length_m;
        let a_ip_cas = ip.r_tip_m * ip.length_m;
        let a_lp_hood = lp.r_tip_m * lp.length_m;
        let a_tot = a_rotor + a_ip_cas + a_lp_hood;
        let iplp_shares = [a_rotor / a_tot, a_ip_cas / a_tot, a_lp_hood / a_tot];

        Self {
            ports: [
                PortSpec::heat(sid("hp_heat"), PortRole::FlowIn).expect("thermal hp_heat"),
                PortSpec::heat(sid("iplp_heat"), PortRole::FlowIn).expect("thermal iplp_heat"),
                PortSpec::heat(sid("bearing_heat"), PortRole::FlowIn)
                    .expect("thermal bearing_heat"),
                PortSpec::heat(sid("ambient_out"), PortRole::EffortIn)
                    .expect("thermal ambient_out"),
                PortSpec::signal_in(sid("controls"), nz(2)), // [oil_cooler_frac, ambient_temp_k]
            ],
            caps_j_k,
            casing_ext_area_m2: 24.0,
            oil_cooler_ua_w_k: 9_500.0,
            hp_rotor_share,
            iplp_shares,
        }
    }

    #[must_use]
    pub fn warm_nominal_state() -> ThermalNetworkState {
        ThermalNetworkState {
            temps_k: [
                720.0, // HP rotor (~447 °C)
                695.0, // HP casing (~422 °C)
                535.0, // IP rotor (~262 °C)
                520.0, // IP casing (~247 °C)
                332.0, // LP exhaust hood (~59 °C)
                348.0, // Front bearing (~75 °C)
                352.0, // Rear bearing (~79 °C)
                326.0, // Lube oil sump (~53 °C)
            ],
            diff_expansion_mm: 0.42,
            hp_thermal_stress_mpa: 68.0,
            hp_creep_rate_per_s: 1.0e-12,
            cumulative_creep_pct: 0.0,
            heat_rejected_w: 155_000.0,
        }
    }

    #[must_use]
    pub fn cold_start_state() -> ThermalNetworkState {
        ThermalNetworkState {
            temps_k: [320.0, 315.0, 312.0, 310.0, 302.0, 315.0, 315.0, 313.15],
            diff_expansion_mm: 0.05,
            hp_thermal_stress_mpa: 8.0,
            hp_creep_rate_per_s: 0.0,
            cumulative_creep_pct: 0.0,
            heat_rejected_w: 12_000.0,
        }
    }
}

impl Component for TurbineCasingThermalComponent {
    type State = ThermalNetworkState;
    type Error = TurbineComponentError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Smooth
    }

    fn stored_energy_j(&self, state: &Self::State) -> Result<f64, Self::Error> {
        let mut e = 0.0;
        for i in 0..8 {
            e += self.caps_j_k[i] * state.temps_k[i];
        }
        Ok(e)
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
        check_cx(cx)?;
        // Inputs: hp_heat Q_hp (0), iplp_heat Q_iplp (1), bearing_heat Q_brg (2),
        // ambient_out T_amb_sink (3), controls [oil_cooler_frac, oil_supply_t_k] (4..6)
        let q_hp = inputs[0];
        let q_iplp = inputs[1];
        let q_brg = inputs[2];
        let t_amb = inputs[3].clamp(273.15, 340.0);
        let oil_cooler_frac = inputs[4].clamp(0.05, 1.5);
        let t_oil_cw = inputs[5].clamp(283.15, 330.0);

        let lib = library();
        let t = old.temps_k;

        // External horizontal-cylinder natural convection via fs_convection::ChurchillChuHorizontalCylinder!
        let t_surf = 0.5 * (t[TH_HP_CASING] + t[TH_IP_CASING]);
        let t_film = 0.5 * (t_surf + t_amb);
        let mu_air = lib.ambient_air.dynamic_viscosity(t_film);
        let rho_air = lib.ambient_air.density(t_film, 101_325.0);
        let nu_air = mu_air / rho_air;
        let alpha_air = nu_air / 0.71;
        let d_cyl: f64 = 1.45;
        let ra = (9.80665 * (1.0 / t_film) * (t_surf - t_amb).abs().max(0.5) * d_cyl.powi(3)
            / (nu_air * alpha_air))
            .clamp(1.0e4, 1.0e11);
        let h_amb_conv = fs_convection::evaluate(
            CorrelationId::ChurchillChuHorizontalCylinder,
            CorrelationInputs::natural(ra, 0.71),
        )
        .and_then(|nu| {
            nu.heat_transfer_coefficient(
                fs_convection::ThermalConductivity::new(0.032),
                Length::new(d_cyl),
            )
        })
        .map(|ev| ev.value.value())
        .unwrap_or(8.0);

        // Radiative + convective ambient loss from HP/IP casing & LP hood
        let q_amb_hp = 0.45
            * self.casing_ext_area_m2
            * (h_amb_conv * (t[TH_HP_CASING] - t_amb)
                + lib.hp_casing_steel.emissivity
                    * STEFAN_BOLTZMANN
                    * (t[TH_HP_CASING].powi(4) - t_amb.powi(4)));
        let q_amb_ip = 0.30
            * self.casing_ext_area_m2
            * (h_amb_conv * (t[TH_IP_CASING] - t_amb)
                + lib.hp_casing_steel.emissivity
                    * STEFAN_BOLTZMANN
                    * (t[TH_IP_CASING].powi(4) - t_amb.powi(4)));
        let q_amb_lp = 0.25
            * self.casing_ext_area_m2
            * (h_amb_conv * (t[TH_LP_HOOD] - t_amb)
                + lib.lp_exhaust_steel.emissivity
                    * STEFAN_BOLTZMANN
                    * (t[TH_LP_HOOD].powi(4) - t_amb.powi(4)));
        // Lube-oil cooler heat rejection to oil cooling water
        let q_oil_cooler = self.oil_cooler_ua_w_k * oil_cooler_frac * (t[TH_LUBE_OIL] - t_oil_cw);

        let q_total_out = q_amb_hp + q_amb_ip + q_amb_lp + q_oil_cooler;

        // Internal conductive couplings between adjacent metal nodes (strictly conservative pairwise fluxes!)
        let q_hpr_ipr = 420.0 * (t[TH_HP_ROTOR] - t[TH_IP_ROTOR]);
        let q_hpc_ipc = 850.0 * (t[TH_HP_CASING] - t[TH_IP_CASING]);
        let q_ipc_lph = 600.0 * (t[TH_IP_CASING] - t[TH_LP_HOOD]);
        let q_hpr_fb = 110.0 * (t[TH_HP_ROTOR] - t[TH_FRONT_BEARING]);
        let q_ipr_rb = 95.0 * (t[TH_IP_ROTOR] - t[TH_REAR_BEARING]);
        let q_fb_oil = 4_200.0 * oil_cooler_frac.powf(0.85) * (t[TH_FRONT_BEARING] - t[TH_LUBE_OIL]);
        let q_rb_oil = 4_200.0 * oil_cooler_frac.powf(0.85) * (t[TH_REAR_BEARING] - t[TH_LUBE_OIL]);

        // Net heat rate into each of the 8 nodes (sum over all 8 == q_hp + q_iplp + q_brg - q_total_out IDENTICALLY!)
        let q_node = [
            self.hp_rotor_share * q_hp - q_hpr_ipr - q_hpr_fb,               // 0: HP rotor
            (1.0 - self.hp_rotor_share) * q_hp - q_hpc_ipc - q_amb_hp,       // 1: HP casing
            self.iplp_shares[0] * q_iplp + q_hpr_ipr - q_ipr_rb,             // 2: IP/LP rotor
            self.iplp_shares[1] * q_iplp + q_hpc_ipc - q_ipc_lph - q_amb_ip, // 3: IP casing
            self.iplp_shares[2] * q_iplp + q_ipc_lph - q_amb_lp,             // 4: LP hood
            0.48 * q_brg + q_hpr_fb - q_fb_oil,                              // 5: Front bearing
            0.52 * q_brg + q_ipr_rb - q_rb_oil,                              // 6: Rear bearing
            q_fb_oil + q_rb_oil - q_oil_cooler,                              // 7: Lube oil sump
        ];

        let mut t_new = [0.0_f64; 8];
        for i in 0..8 {
            t_new[i] = t[i] + dt_s * q_node[i] / self.caps_j_k[i];
            if !(t_new[i].is_finite() && t_new[i] > 200.0 && t_new[i] < 1_800.0) {
                return Ok(Proposal::Shorten { dt_s: 0.5 * dt_s });
            }
        }

        // Exact differential axial expansion ΔL = L_hp*(ε_rotor - ε_casing) + L_ip*(ε_ipr - ε_ipc) [mm]
        let dl_hp = 1.60
            * (lib.rotor_steel.thermal_strain(t_new[TH_HP_ROTOR])
                - lib.hp_casing_steel.thermal_strain(t_new[TH_HP_CASING]));
        let dl_ip = 2.00
            * (lib.rotor_steel.thermal_strain(t_new[TH_IP_ROTOR])
                - lib.hp_casing_steel.thermal_strain(t_new[TH_IP_CASING]));
        let diff_expansion_mm = (dl_hp + dl_ip) * 1.0e3;

        // Thermal hoop stress on thick HP casing from rotor-to-casing ΔT gradient [MPa]
        let e_hp = lib.hp_casing_steel.young_modulus_at(t_new[TH_HP_CASING]);
        let a_hp = lib.hp_casing_steel.mean_alpha_at(t_new[TH_HP_CASING]);
        let nu_hp = lib.hp_casing_steel.poisson_ratio;
        let dt_grad = (t_new[TH_HP_ROTOR] - t_new[TH_HP_CASING]).abs();
        let sigma_th_pa = e_hp * a_hp * dt_grad / (2.0 * (1.0 - nu_hp));
        // Add nominal centrifugal + pressure base von Mises stress (~140 MPa) for Norton creep evaluation
        let sigma_eq_pa = (140.0e6_f64).hypot(sigma_th_pa);
        let creep_rate = lib.rotor_creep_rate_per_s(sigma_eq_pa, t_new[TH_HP_ROTOR]);

        *new = ThermalNetworkState {
            temps_k: t_new,
            diff_expansion_mm,
            hp_thermal_stress_mpa: sigma_th_pa * 1.0e-6,
            hp_creep_rate_per_s: creep_rate,
            cumulative_creep_pct: old.cumulative_creep_pct + creep_rate * dt_s * 100.0,
            heat_rejected_w: q_total_out,
        };

        // Outputs: midpoint temperatures on FlowIn heat ports; -q_total_out on ambient_out EffortIn port!
        let t_brg_old = 0.5 * (t[TH_FRONT_BEARING] + t[TH_REAR_BEARING]);
        let t_brg_new = 0.5 * (t_new[TH_FRONT_BEARING] + t_new[TH_REAR_BEARING]);
        outputs[0] = 0.5 * (t[TH_HP_CASING] + t_new[TH_HP_CASING]);
        outputs[1] = 0.5 * (t[TH_IP_CASING] + t_new[TH_IP_CASING]);
        outputs[2] = 0.5 * (t_brg_old + t_brg_new);
        outputs[3] = -q_total_out;

        // Exact First-Law stored energy change
        let de = dt_s * (q_hp + q_iplp + q_brg - q_total_out);
        // Second-Law internal conductive & ambient entropy production >= 0
        let mut s_gen_w_k = 0.0_f64;
        for (q, ta, tb) in [
            (q_hpr_ipr, t[TH_HP_ROTOR], t[TH_IP_ROTOR]),
            (q_hpc_ipc, t[TH_HP_CASING], t[TH_IP_CASING]),
            (q_ipc_lph, t[TH_IP_CASING], t[TH_LP_HOOD]),
            (q_hpr_fb, t[TH_HP_ROTOR], t[TH_FRONT_BEARING]),
            (q_ipr_rb, t[TH_IP_ROTOR], t[TH_REAR_BEARING]),
            (q_fb_oil, t[TH_FRONT_BEARING], t[TH_LUBE_OIL]),
            (q_rb_oil, t[TH_REAR_BEARING], t[TH_LUBE_OIL]),
        ] {
            // q has the sign of (ta - tb), so each term is >= 0 without clamping.
            s_gen_w_k += q * (1.0 / tb - 1.0 / ta);
        }
        let e_total = self.stored_energy_j(new)?;
        Ok(Proposal::Step(Receipt {
            stored_change_j: de,
            dissipated_j: 0.0,
            numerical_loss_j: 0.0,
            entropy_generated_j_per_k: dt_s * s_gen_w_k,
            energy_tolerance_j: 1.0e-8 * e_total + 1.0e-3,
            residual: 0.0,
            tolerance: 0.0,
        }))
    }
}

impl PortTangent for TurbineCasingThermalComponent {
    fn jvp(
        &self,
        _old: &Self::State,
        dt_s: f64,
        _inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), Self::Error> {
        check_cx(cx)?;
        // Exact analytical linearization wrt the 3 coupled FlowIn heat ports!
        // (ambient_out T_amb and controls are fixed boundaries, not loop unknowns.)
        let k0 = 0.5 * dt_s * (1.0 - self.hp_rotor_share) / self.caps_j_k[TH_HP_CASING];
        let k1 = 0.5 * dt_s * self.iplp_shares[1] / self.caps_j_k[TH_IP_CASING];
        let k2 = 0.5 * dt_s * 0.5 * (0.48 / self.caps_j_k[TH_FRONT_BEARING] + 0.52 / self.caps_j_k[TH_REAR_BEARING]);
        d_outputs[0] = k0 * d_inputs[0];
        d_outputs[1] = k1 * d_inputs[1];
        d_outputs[2] = k2 * d_inputs[2];
        d_outputs[3] = 0.0;
        Ok(())
    }

    fn vjp(
        &self,
        _old: &Self::State,
        dt_s: f64,
        _inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), Self::Error> {
        check_cx(cx)?;
        let k0 = 0.5 * dt_s * (1.0 - self.hp_rotor_share) / self.caps_j_k[TH_HP_CASING];
        let k1 = 0.5 * dt_s * self.iplp_shares[1] / self.caps_j_k[TH_IP_CASING];
        let k2 = 0.5 * dt_s * 0.5 * (0.48 / self.caps_j_k[TH_FRONT_BEARING] + 0.52 / self.caps_j_k[TH_REAR_BEARING]);
        bar_inputs.fill(0.0);
        bar_inputs[0] = k0 * bar_outputs[0];
        bar_inputs[1] = k1 * bar_outputs[1];
        bar_inputs[2] = k2 * bar_outputs[2];
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// 9. GovernorControlComponent (fs-control SpeedController + ThermalDeratingPolicy)
// ---------------------------------------------------------------------------

/// Electro-hydraulic steam turbine governor & protection system backed by
/// [`fs_control::speed::SpeedController`] and [`fs_control::thermal::ThermalDeratingPolicy`].
///
/// Implements:
/// - 4% ISO speed droop + PI load-frequency control
/// - 110% overspeed (3,300 rpm) Emergency Stop Valve (ESV) hard trip latch
/// - Bearing Babbitt over-temperature (> 120 °C / 393 K) & LP exhaust over-temperature (> 125 °C / 398 K)
///   derating and trip latches via `ThermalDeratingPolicy`
/// - Hydraulic servo rate limits (0.25 s full ESV slam close; 1.8 s CV opening ramp)
pub struct GovernorControlComponent {
    ports: [PortSpec; 2],
    speed_ctrl: SpeedController,
    thermal_bands: [ThermalDeratingBand; 2],
    pub rated_omega_rad_s: f64,
    pub overspeed_trip_rad_s: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct GovernorState {
    pub cv_lift: f64,
    pub esv_open: f64,
    pub ext_lift: f64,
    pub overspeed_tripped: bool,
    pub vibration_tripped: bool,
    pub manual_tripped: bool,
    pub thermal_state: ThermalDeratingState,
    pub pi_state: SpeedControllerState,
    pub derating_fraction: f64,
}

impl GovernorControlComponent {
    #[must_use]
    pub fn new() -> Self {
        let speed_ctrl = SpeedController::new(SpeedControllerParameters {
            sample_period_s: 0.02,
            proportional_gain_nm_s_per_rad: 2.8e4,
            integral_gain_nm_per_rad: 1.2e4,
            anti_windup_time_s: 0.25,
        })
        .expect("valid fs-control SpeedController");
        // Channel 0: Bearing Babbitt temp (derate > 105 C, trip >= 120 C = 393.15 K matching matdb wipe onset)
        // Channel 1: LP Exhaust hood temp (derate > 100 C, trip >= 130 C = 403.15 K)
        let thermal_bands = [
            ThermalDeratingBand {
                full_output_temperature_k: 378.15,
                trip_temperature_k: 393.15,
                reset_temperature_k: 365.15,
            },
            ThermalDeratingBand {
                full_output_temperature_k: 373.15,
                trip_temperature_k: 403.15,
                reset_temperature_k: 358.15,
            },
        ];
        Self {
            ports: [
                // Inputs: [omega_rad_s, speed_ref_rad_s, load_setpoint_frac, ext_demand_frac,
                //          bearing_max_t_k, lp_hood_t_k, max_vib_um_pp, manual_trip, reset_trip]
                PortSpec::signal_in(sid("sense"), nz(9)),
                // Outputs: [effective_hp_valve_lift, ext_valve_lift]
                PortSpec::signal_out(sid("valves"), nz(2)),
            ],
            speed_ctrl,
            thermal_bands,
            rated_omega_rad_s: 100.0 * PI,        // 3,000 rpm
            overspeed_trip_rad_s: 110.0 * PI,     // 3,300 rpm (110% ISO overspeed trip!)
        }
    }

    #[must_use]
    pub fn initial_state(cv_lift: f64, ext_lift: f64) -> GovernorState {
        GovernorState {
            cv_lift,
            esv_open: 1.0,
            ext_lift,
            overspeed_tripped: false,
            vibration_tripped: false,
            manual_tripped: false,
            thermal_state: ThermalDeratingState { tripped: false },
            pi_state: SpeedControllerState::default(),
            derating_fraction: 1.0,
        }
    }
}

impl Default for GovernorControlComponent {
    fn default() -> Self {
        Self::new()
    }
}

impl Component for GovernorControlComponent {
    type State = GovernorState;
    type Error = TurbineComponentError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Smooth
    }

    fn stored_energy_j(&self, _state: &Self::State) -> Result<f64, Self::Error> {
        Ok(0.0)
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
        check_cx(cx)?;
        let omega = inputs[0];
        let omega_ref = inputs[1].clamp(0.0, 130.0 * PI);
        let load_sp = inputs[2].clamp(0.0, 1.15);
        let ext_sp = inputs[3].clamp(0.0, 1.0);
        let t_brg = inputs[4].clamp(274.0, 600.0);
        let t_lp = inputs[5].clamp(274.0, 700.0);
        let max_vib = inputs[6].max(0.0);
        let manual_trip = inputs[7] >= 0.5;
        let reset_req = inputs[8] >= 0.5;

        // 1. Evaluate fs_control::thermal::ThermalDeratingPolicy
        let policy = ThermalDeratingPolicy::new(&self.thermal_bands)
            .expect("valid thermal derating policy");
        let base_limits = ShaftPowerLimits {
            positive_torque_limit_nm: 2.5e5,
            negative_torque_limit_nm: 5.0e4,
            motoring_power_w: 80.0e6,
            regeneration_power_w: 20.0e6,
        };
        let th_prop = policy
            .propose_controlled(
                &old.thermal_state,
                &[t_brg, t_lp],
                base_limits,
                reset_req,
                || false,
            )
            .expect("valid thermal derating proposal");

        // 2. Overspeed (110% = 3,300 rpm) & ISO 7919 Zone D (> 240 µm) trip latches
        let overspeed_tripped = if reset_req && omega < 0.96 * self.rated_omega_rad_s {
            false
        } else {
            old.overspeed_tripped || omega >= self.overspeed_trip_rad_s
        };
        let vibration_tripped = if reset_req && max_vib < 120.0 {
            false
        } else {
            old.vibration_tripped || max_vib >= crate::rotordynamics::ISO_7919_TRIP_UM_PP
        };
        let any_trip =
            overspeed_tripped || vibration_tripped || manual_trip || th_prop.next_state.tripped;

        // 3. 4% ISO droop + fs_control::speed::SpeedController PI trim
        let droop_frac = 0.04;
        let speed_err_pu = (omega_ref - omega) / self.rated_omega_rad_s;
        let droop_cmd = load_sp + speed_err_pu / droop_frac;

        let demand = self.speed_ctrl.demand_controlled(
            &old.pi_state,
            [omega, 0.0],
            SpeedCommand {
                relative_speed_rad_s: omega_ref,
                feedforward_torque_nm: load_sp * 1.6e5,
            },
            th_prop.limits,
            || false,
        );
        let (pi_next, pi_trim) = if let Ok(d) = demand {
            let prop = d.finish_controlled(d.requested_torque_nm, || false);
            let next_st = prop.map(|p| p.next_state).unwrap_or(old.pi_state);
            let trim = ((d.requested_torque_nm - load_sp * 1.6e5) / 2.5e5).clamp(-0.35, 0.35);
            (next_st, trim)
        } else {
            (old.pi_state, 0.0)
        };

        let cv_target = if any_trip {
            0.0
        } else {
            ((droop_cmd + 0.15 * pi_trim) * th_prop.output_fraction).clamp(0.0, 1.0)
        };
        let esv_target = if any_trip { 0.0 } else { 1.0 };
        let ext_target = if any_trip { 0.0 } else { ext_sp };

        // Realistic IEC 60045-1 / API 612 quick-closing electro-hydraulic servo dynamics:
        // ESV slams shut in 0.125 s (8.0/s), HP CV closes in 0.167 s (6.0/s), smooth opening (~1.5 s).
        let slew = |cur: f64, tgt: f64, open_rate: f64, close_rate: f64| -> f64 {
            let diff = tgt - cur;
            let max_step = if diff >= 0.0 { open_rate * dt_s } else { close_rate * dt_s };
            (cur + diff.clamp(-max_step, max_step)).clamp(0.0, 1.0)
        };
        let cv_new = slew(old.cv_lift, cv_target, 0.65, 6.0);
        let esv_new = slew(old.esv_open, esv_target, 0.45, 8.0);
        let ext_new = slew(old.ext_lift, ext_target, 0.55, 3.0);

        *new = GovernorState {
            cv_lift: cv_new,
            esv_open: esv_new,
            ext_lift: ext_new,
            overspeed_tripped,
            vibration_tripped,
            manual_tripped: manual_trip,
            thermal_state: th_prop.next_state,
            pi_state: pi_next,
            derating_fraction: th_prop.output_fraction,
        };

        outputs[0] = cv_new * esv_new;
        outputs[1] = ext_new;

        Ok(Proposal::Step(Receipt::lossless(0.0, 1.0e-12)))
    }
}
