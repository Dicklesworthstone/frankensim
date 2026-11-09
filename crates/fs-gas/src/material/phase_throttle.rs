//! Adiabatic throttling with original mixture enthalpy and phase equilibrium.
//!
//! A passive restriction dissipates pressure without shaft work. Its steady
//! energy balance is `mass_flow * (h_out - h_in) = 0` when kinetic and potential
//! energy changes are negligible. The original pressure/enthalpy owner resolves
//! Joule--Thomson heating or cooling and phase redistribution without replacing
//! the species calorics by a constant Cp or adding a separate latent source.
//!
//! Process convention: MIT 5.60, lecture 4, Joule--Thomson expansion,
//! <https://ocw.mit.edu/courses/5-60-thermodynamics-kinetics-spring-2008/e1f8054833e8077bf9428c27cbf0465c_5_60_lecture4.pdf>.

use fs_material::thermochemical_gas::nonideal::flash::{
    ThermochemicalFlashError,
    equilibrium::{PhaseEnthalpyControls, PhaseEnthalpySolution},
};
use fs_qty::{MassFlowRate, Pressure, Temperature};

/// Original source, inlet state and positive supplied mass flow. The same feed
/// type is accepted by the existing heat-exchanger and mixing consumers.
pub use super::phase_heat_exchanger::MixtureHeatExchangerInlet as MixtureThrottleInlet;

/// Original phase/caloric policy and physical steady-flow acceptance allowances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MixtureThrottleControl {
    /// Original temperature bracket, phase policy and whole-call PT work limit.
    /// The source enthalpy tolerance is additionally bounded by the power
    /// allowance divided by the actual mass flow before the solve begins.
    pub equilibrium: PhaseEnthalpyControls,
    /// Allowed outgoing-minus-incoming enthalpy flow [W].
    pub absolute_power_tolerance_w: f64,
    /// Allowed numerical shortfall from nonnegative entropy production [W/K].
    /// Actual negative roundoff is reported, never clipped to zero.
    pub entropy_tolerance_w_k: f64,
}
impl MixtureThrottleControl {
    /// An explicit original-source outlet temperature interval, which may cross
    /// liquid, coexistence and vapor regions under the selected phase policy.
    #[must_use]
    pub fn new(temperature_bracket: [Temperature; 2]) -> Self {
        Self {
            equilibrium: PhaseEnthalpyControls::new(temperature_bracket),
            absolute_power_tolerance_w: 1e-5,
            entropy_tolerance_w_k: 1e-7,
        }
    }

    /// Admit physical tolerances and the bounded source policy before work.
    pub fn validate(self) -> Result<(), MixtureThrottleError> {
        self.equilibrium.validate().map_err(source_error)?;
        if !positive(self.absolute_power_tolerance_w)
            || !self.entropy_tolerance_w_k.is_finite()
            || self.entropy_tolerance_w_k < 0.0
        {
            return Err(MixtureThrottleError::Invalid(
                "positive power and nonnegative entropy allowances required",
            ));
        }
        Ok(())
    }
}

/// A refused throttle publishes no outlet or replacement source model.
#[derive(Clone, Debug, PartialEq)]
pub enum MixtureThrottleError {
    /// Invalid state, uphill pressure direction or unrepresentable flow account.
    Invalid(&'static str),
    /// Original phase, caloric, source-domain or bounded-work refusal.
    Equilibrium(ThermochemicalFlashError),
    /// Actual stream enthalpy failed the declared steady first-law allowance.
    Energy { residual_w: f64, tolerance_w: f64 },
    /// Actual original entropy flow violated the adiabatic second-law allowance.
    Entropy {
        production_w_k: f64,
        tolerance_w_k: f64,
    },
    /// Cancellation observed before publication of the complete outlet.
    Cancelled,
}
impl core::fmt::Display for MixtureThrottleError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adiabatic mixture throttle: {self:?}")
    }
}
impl std::error::Error for MixtureThrottleError {}

/// Original accepted outlet and its steady mass, energy and entropy accounts.
#[derive(Clone, Debug, PartialEq)]
pub struct MixtureThrottleStep {
    /// Original phase outlet and source-integrated Delta h/Delta s. Its work
    /// counters cover the entire solve, including one inlet phase selection.
    pub equilibrium: PhaseEnthalpySolution,
    /// Identical positive inlet/outlet mass flow; species remain unchanged.
    pub mass_flow: MassFlowRate,
    /// Actual `mass_flow * Delta h` [W], retained as a numerical defect. It is
    /// neither delivered heat nor shaft work; both are zero in this model.
    pub energy_residual_w: f64,
    /// Actual `mass_flow * Delta s` [W/K], without clipping roundoff.
    pub entropy_production_w_k: f64,
}

/// Reduce a steady stream's pressure adiabatically with no shaft work.
///
/// Resolves `Delta h = 0` through the retained species laws, EOS and equilibrium
/// phase owner. Original formation and entropy reference offsets cancel before
/// rounding. Flow work is already part of stream enthalpy: no extra `p*v` debit,
/// isentropic work, fitted heat capacity or latent-heat source is applied.
/// Temperature may rise or fall, and the original phase owner may change the
/// phase amounts. There is no forced cooling or temperature average.
///
/// Outlet pressure and positive mass flow are supplied; this API does not
/// determine valve area, pressure loss, discharge coefficient or choking. It
/// assumes steady flow with negligible inlet/outlet kinetic and potential
/// energy differences, no external heat, no shaft work, no reaction and phase
/// equilibrium. It supplies no slip/phase kinetics or global-stability proof.
/// A pressure increase is refused. At equal pressure the native exact-rest path
/// uses one inlet PT query when its temperature lies in the declared bracket.
pub fn throttle_adiabatic_mixture_stream(
    inlet: MixtureThrottleInlet<'_>,
    outlet_pressure: Pressure,
    mut controls: MixtureThrottleControl,
    mut cancelled: impl FnMut() -> bool,
) -> Result<MixtureThrottleStep, MixtureThrottleError> {
    poll(&mut cancelled)?;
    controls.validate()?;
    let flow = inlet.mass_flow.value();
    if !positive(inlet.temperature.value())
        || !positive(inlet.pressure.value())
        || !positive(outlet_pressure.value())
        || !positive(flow)
        || outlet_pressure.value() > inlet.pressure.value()
    {
        return Err(MixtureThrottleError::Invalid(
            "positive inlet state/flow and a positive nonincreasing pressure required",
        ));
    }
    controls.equilibrium.enthalpy.enthalpy_tolerance_j_kg = controls
        .equilibrium
        .enthalpy
        .enthalpy_tolerance_j_kg
        .min(0.5 * controls.absolute_power_tolerance_w / flow);
    controls.validate()?;
    let equilibrium = inlet
        .model
        .phase_equilibrium_pressure_enthalpy_increment(
            inlet.temperature,
            inlet.pressure,
            outlet_pressure,
            0.0,
            controls.equilibrium,
            &mut cancelled,
        )
        .map_err(source_error)?;
    let delta_h = equilibrium
        .accepted_increment_j_kg()
        .ok_or(MixtureThrottleError::Invalid(
            "missing original throttle enthalpy increment",
        ))?;
    let delta_s = equilibrium
        .entropy_change_j_kg_k()
        .ok_or(MixtureThrottleError::Invalid(
            "missing original throttle entropy increment",
        ))?;
    let energy_residual_w = flow * delta_h;
    if !energy_residual_w.is_finite()
        || energy_residual_w.abs() > controls.absolute_power_tolerance_w
    {
        return Err(MixtureThrottleError::Energy {
            residual_w: energy_residual_w,
            tolerance_w: controls.absolute_power_tolerance_w,
        });
    }
    let entropy_production_w_k = flow * delta_s;
    if !entropy_production_w_k.is_finite()
        || entropy_production_w_k < -controls.entropy_tolerance_w_k
    {
        return Err(MixtureThrottleError::Entropy {
            production_w_k: entropy_production_w_k,
            tolerance_w_k: controls.entropy_tolerance_w_k,
        });
    }
    poll(&mut cancelled)?;
    Ok(MixtureThrottleStep {
        equilibrium,
        mass_flow: inlet.mass_flow,
        energy_residual_w,
        entropy_production_w_k,
    })
}

fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), MixtureThrottleError> {
    if cancelled() {
        Err(MixtureThrottleError::Cancelled)
    } else {
        Ok(())
    }
}
fn source_error(error: ThermochemicalFlashError) -> MixtureThrottleError {
    if error == ThermochemicalFlashError::Cancelled {
        MixtureThrottleError::Cancelled
    } else {
        MixtureThrottleError::Equilibrium(error)
    }
}
