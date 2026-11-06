//! Adiabatic homogeneous-equilibrium response of original IF97 liquid/vapor.
//!
//! Both phases share pressure, temperature and velocity, and phase exchange is
//! fast on the acoustic timescale. This is the equilibrium limit, complementary
//! to the frozen phase masses in `FrozenHomogeneousMixture`; it is not a fitted
//! interpolation, finite-rate phase-transfer law or bubble-resonance model.
//!
//! For saturation temperature T and vapor mass fraction x, the ORIGINAL source
//! gives u(T,x) and v(T,x). A reversible adiabatic volume perturbation obeys
//! du + p dv = 0, hence x_T = -(u_T + p v_T)/(u_v-u_l+p(v_v-v_l)).
//! The denominator uses the original internal-energy and volume gaps at the
//! pair's declared pressure. It equals latent enthalpy when both phase EOS
//! pressures equal that pressure. Dense-root residuals remain explicit, and
//! each phase's caloric tangent includes its pressure-work correction. Together
//! with the original Region 4 p_T this gives K = -v p_T/(v_T + (v_v-v_l) x_T).
//! No finite differences, inverse solves, altered caloric datums or allocations
//! are used in production evaluation.
//!
//! IF97 has finite saturation matching defects (release section 10.2). The
//! tangent conserves its original u/v pressure work; its small RAW-source
//! entropy derivative and Clapeyron/Gibbs defects remain visible. They are
//! source inconsistency, not physical entropy production or added heat.
//! In exact thermodynamic coexistence this is the usual isentropic HEM speed.
//! References: IAPWS R7-97(2012), sections 5, 6, 8 and 10.2,
//! <https://iapws.org/technical-guidance/release/IF97-Rev.download>;
//! Brennen, Cavitation and Bubble Dynamics, section 6.3,
//! <https://media.library.caltech.edu/CaltechBOOK:1995.001/chap6.htm>.

use crate::MaterialError;
use crate::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_thermochem::water::if97::regions::{If97FluidSaturationPair, If97FluidState};
use fs_thermochem::water::if97::{If97SaturationPair, saturation_pressure_derivative};

/// Local adiabatic tangent and unmodified source consistency diagnostics.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquilibriumWaterAcousticTangent {
    /// Original Region 4 saturation derivative [Pa/K].
    pub saturation_pressure_derivative_pa_k: f64,
    /// Vapor MASS fraction derivative on this energy-conserving tangent [1/Pa].
    pub vapor_mass_fraction_derivative_per_pa: f64,
    /// Specific-volume derivative on this tangent [m³/(kg Pa)], negative.
    pub specific_volume_derivative_m3_kg_pa: f64,
    /// Specific-internal-energy derivative [J/(kg Pa)] = [m³/kg].
    pub specific_internal_energy_derivative_m3_kg: f64,
    /// RAW source ds/dp [J/(kg K Pa)]. Ideally zero; not an entropy source.
    pub source_entropy_derivative_j_kg_k_pa: f64,
    /// Original g_v - g_l [J/kg]; no endpoint datum correction is applied.
    pub source_gibbs_gap_j_kg: f64,
    /// (p_T (v_v-v_l) - (s_v-s_l))/(s_v-s_l), ideally zero.
    pub source_relative_clapeyron_defect: f64,
    /// Actual phase EOS pressures minus the pair's declared pressure [liquid,vapor].
    pub source_phase_pressure_residuals_pa: [f64; 2],
    /// Original Region 4 p(T) minus the pair's declared pressure [Pa].
    pub source_saturation_pressure_residual_pa: f64,
    /// Original h_v-h_l minus [u_v-u_l+p*(v_v-v_l)] [J/kg].
    /// Dense phase pressure residuals need not make this exactly zero.
    pub source_latent_enthalpy_defect_j_kg: f64,
}

/// Original saturation pair and an explicit homogeneous-equilibrium closure.
///
/// Both phase masses must be strictly positive. At a saturation endpoint the
/// phase-changing and single-phase tangents are different one-sided responses;
/// this type refuses to choose one silently. A finite perturbation must remain
/// inside coexistence and the original saturation temperature domain. The
/// caller must establish rapid phase/thermal equilibration and negligible slip,
/// interface energy and noncondensable gas for the intended frequency band.
///
/// The default source retains the original Region 1/2 pair. The explicit
/// regional constructor retains an `If97FluidSaturationPair`, including its
/// Region 3 roots and solve diagnostics. Both use the same response algebra;
/// no new source query, root solve or near-critical extrapolation is performed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HomogeneousEquilibriumWater<S = If97SaturationPair> {
    source: S,
    vapor_mass_fraction: f64,
    vapor_volume_fraction: f64,
    acoustic: FluidAcousticState,
    tangent: EquilibriumWaterAcousticTangent,
}

fn invalid(what: &str) -> MaterialError {
    MaterialError::State { what: what.into() }
}

// The second component is du/dT + p_pair dv/dT along the ORIGINAL saturation
// line. The Gibbs identity first gives du/dT + p_phase dv/dT; retain the dense
// root's pressure-work correction separately from the original entropy slope.
fn along_saturation(
    state: If97FluidState,
    pressure: f64,
    pressure_slope: f64,
) -> (f64, f64, f64) {
    let vt = state.dv_dt_at_pressure().value();
    let vp = state.dv_dp_at_temperature().value();
    let temperature = state.temperature().value();
    let volume_slope = vt + vp * pressure_slope;
    let adiabatic_caloric_slope = state.cp().value() - temperature * vt * pressure_slope;
    (
        volume_slope,
        adiabatic_caloric_slope + (pressure - state.pressure().value()) * volume_slope,
        adiabatic_caloric_slope / temperature,
    )
}

impl HomogeneousEquilibriumWater {
    /// Equilibrium response retaining the original Region 1/2 saturation pair.
    pub fn try_new(
        source: If97SaturationPair,
        vapor_mass_fraction: f64,
    ) -> Result<Self, MaterialError> {
        Self::from_source(
            source,
            If97FluidSaturationPair::Gibbs(source),
            vapor_mass_fraction,
        )
    }
}

impl HomogeneousEquilibriumWater<If97FluidSaturationPair> {
    /// Equilibrium response throughout the original admitted coexistence domain.
    ///
    /// The existing regional owner supplies both phases, including the dense
    /// roots above 623.15 K. Work uses its declared pair pressure; the actual
    /// phase pressures, Region 4 round-trip defect and caloric mismatch remain
    /// visible. Exact critical or unresolved source states must already refuse
    /// in that owner, and nonpositive/unrepresentable acoustic response refuses
    /// here. No source energy, entropy, pressure or density is replaced.
    pub fn try_from_regional(
        source: If97FluidSaturationPair,
        vapor_mass_fraction: f64,
    ) -> Result<Self, MaterialError> {
        Self::from_source(source, source, vapor_mass_fraction)
    }
}

impl<S> HomogeneousEquilibriumWater<S> {
    fn from_source(
        retained_source: S,
        source: If97FluidSaturationPair,
        vapor_mass_fraction: f64,
    ) -> Result<Self, MaterialError> {
        let x = vapor_mass_fraction;
        if !x.is_finite() || x <= 0.0 || x >= 1.0 {
            return Err(invalid(
                "equilibrium two-phase acoustics requires both phase masses",
            ));
        }
        let liquid = source.liquid();
        let vapor = source.vapor();
        let vl = liquid.specific_volume().value();
        let vv = vapor.specific_volume().value();
        let volume_gap = vv - vl;
        let enthalpy_gap = source.latent_enthalpy().value();
        let entropy_gap = source.latent_entropy().value();
        let pressure = source.pressure().value();
        let energy_gap =
            vapor.specific_internal_energy().value() - liquid.specific_internal_energy().value();
        let phase_change_energy = energy_gap + pressure * volume_gap;
        let volume = (1.0 - x) * vl + x * vv;
        let pressure_slope = saturation_pressure_derivative(source.temperature())
            .map_err(|error| invalid(&format!("equilibrium acoustic source: {error}")))?
            .value();
        if [
            volume_gap,
            enthalpy_gap,
            entropy_gap,
            phase_change_energy,
            volume,
            pressure_slope,
        ]
        .into_iter()
        .any(|v| !v.is_finite() || v <= 0.0)
        {
            return Err(invalid("unresolved original water coexistence response"));
        }
        let (vl_t, ql_t, sl_t) = along_saturation(liquid, pressure, pressure_slope);
        let (vv_t, qv_t, sv_t) = along_saturation(vapor, pressure, pressure_slope);
        let mix_v_t = (1.0 - x) * vl_t + x * vv_t;
        let mix_q_t = (1.0 - x) * ql_t + x * qv_t;
        let quality_t = -mix_q_t / phase_change_energy;
        let volume_t = mix_v_t + volume_gap * quality_t;
        let volume_p = volume_t / pressure_slope;
        // Evaluate source du separately so the displayed account does not hide
        // roundoff by assigning du=-p*dv after the fact.
        let ul_t = ql_t - pressure * vl_t;
        let uv_t = qv_t - pressure * vv_t;
        let energy_t = (1.0 - x) * ul_t + x * uv_t + energy_gap * quality_t;
        let tangent = EquilibriumWaterAcousticTangent {
            saturation_pressure_derivative_pa_k: pressure_slope,
            vapor_mass_fraction_derivative_per_pa: quality_t / pressure_slope,
            specific_volume_derivative_m3_kg_pa: volume_p,
            specific_internal_energy_derivative_m3_kg: energy_t / pressure_slope,
            source_entropy_derivative_j_kg_k_pa: ((1.0 - x) * sl_t
                + x * sv_t
                + entropy_gap * quality_t)
                / pressure_slope,
            source_gibbs_gap_j_kg: source.gibbs_energy_gap().value(),
            source_relative_clapeyron_defect: (pressure_slope * volume_gap - entropy_gap)
                / entropy_gap,
            source_phase_pressure_residuals_pa: [
                liquid.pressure().value() - pressure,
                vapor.pressure().value() - pressure,
            ],
            source_saturation_pressure_residual_pa: source.saturation_pressure_residual_pa(),
            source_latent_enthalpy_defect_j_kg: enthalpy_gap - phase_change_energy,
        };
        let alpha = x * vv / volume;
        if !volume_p.is_finite()
            || volume_p >= 0.0
            || !alpha.is_finite()
            || alpha <= 0.0
            || alpha >= 1.0
            || [
                tangent.vapor_mass_fraction_derivative_per_pa,
                tangent.specific_internal_energy_derivative_m3_kg,
                tangent.source_entropy_derivative_j_kg_k_pa,
                tangent.source_gibbs_gap_j_kg,
                tangent.source_relative_clapeyron_defect,
                tangent.source_phase_pressure_residuals_pa[0],
                tangent.source_phase_pressure_residuals_pa[1],
                tangent.source_saturation_pressure_residual_pa,
                tangent.source_latent_enthalpy_defect_j_kg,
            ]
            .into_iter()
            .any(|v| !v.is_finite())
        {
            return Err(invalid(
                "unrepresentable stable equilibrium acoustic tangent",
            ));
        }
        let acoustic =
            FluidAcousticState::from_adiabatic_bulk_modulus(1.0 / volume, -volume / volume_p)?;
        Ok(Self {
            source: retained_source,
            vapor_mass_fraction: x,
            vapor_volume_fraction: alpha,
            acoustic,
            tangent,
        })
    }

    /// Original unmodified source, including regional root diagnostics when used.
    pub const fn source(&self) -> &S {
        &self.source
    }
    pub const fn vapor_mass_fraction(&self) -> f64 {
        self.vapor_mass_fraction
    }
    pub const fn vapor_volume_fraction(&self) -> f64 {
        self.vapor_volume_fraction
    }
    pub const fn state(&self) -> FluidAcousticState {
        self.acoustic
    }
    pub const fn tangent(&self) -> EquilibriumWaterAcousticTangent {
        self.tangent
    }
}

impl<S> FluidAcoustics for HomogeneousEquilibriumWater<S> {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        Ok(self.acoustic)
    }
}

#[cfg(test)]
#[path = "acoustics/tests.rs"]
mod tests;
