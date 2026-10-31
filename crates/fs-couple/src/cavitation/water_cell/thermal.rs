//! Finite wall storage and equilibrium water coupled by native conduction.
//!
//! One signed heat Q leaves the original wall enthalpy curve and enters the
//! conserved internal energy of the rigid water cell. The existing conducting
//! stack evaluates its original k(T) integrals at BOTH proposed endpoint
//! temperatures. The shared bounded scalar solver finds Q = dt * P_end;
//! this is backward Euler, with first-order time accuracy.
//!
//! The wall is a declared uniform finite thermal body with negligible pV work.
//! Its original enthalpy/phase chart remains the sole caloric owner. The
//! intervening layers are separate, source-free, quasi-steady thermal
//! resistances with negligible storage. Their masses must not be counted again
//! as the uniform body's storage. Uniform water temperature and negligible
//! unresolved fluid-film resistance are explicit assumptions. This is not a
//! boiling/condensation film coefficient, a moving phase-front calculation,
//! or a justification of lumping from geometry alone.
//!
//! The water retains its original IF97 equilibrium, finite liquid/vapor
//! inventory, fixed volume and conserved U. No latent heat is added separately.
//! Source inversion, interface-law and conserved-energy residuals are distinct.
//! External heat has no specified boundary temperature, so the exposed entropy
//! remainder is not certified entropy production.
//!
//! Fourier resistance reference: MIT Unified Engineering, section 16.4,
//! https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node118.html

use super::{EquilibriumWaterCell, WaterCellError, WaterCellStep};
pub use crate::cavity_phs::piston::thermal::HeatExchangeSolve;
use crate::cavity_phs::piston::thermal::feedback;
use fs_material::phase::{
    EquilibriumPhaseState, PhaseStateError, UniformEnthalpyStep, UniformEnthalpyStepInput,
};
use fs_material::thermal::conduction::{
    LayerConductionConfig, LayerStackConduction, LayerStackResponse,
};
use fs_thermal::enthalpy::DepositionError;

/// Whole-step work, endpoint temperature and native conduction admission.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterWallControls {
    /// Q is positive from wall to water. Every queried bracket state must be
    /// supported by both original material laws and the temperature limits.
    pub exchange: HeatExchangeSolve,
    pub conduction: LayerConductionConfig,
    /// Positive interval ceiling; not an estimate of temporal truncation error.
    pub maximum_step_s: f64,
    /// Positive absolute temperature-change limits [wall, water] in kelvin.
    pub maximum_temperature_increment_k: [f64; 2],
}

/// Original material/flow refusals and bounded coupled-solve failures.
#[derive(Clone, Debug, PartialEq)]
pub enum WaterWallError {
    Invalid(&'static str),
    Water(WaterCellError),
    Wall(PhaseStateError),
    Conduction(DepositionError),
    TemperatureIncrement {
        body: &'static str,
        increment_k: f64,
        maximum_k: f64,
    },
    Budget {
        quantity: &'static str,
        residual_j: f64,
        tolerance_j: f64,
    },
    OutsideBracket {
        lower_residual_j: f64,
        upper_residual_j: f64,
    },
    NotConverged {
        iterations: usize,
        residual_j: f64,
    },
    Cancelled,
}
impl From<WaterCellError> for WaterWallError {
    fn from(error: WaterCellError) -> Self {
        match error {
            WaterCellError::Cancelled => Self::Cancelled,
            other => Self::Water(other),
        }
    }
}
impl From<PhaseStateError> for WaterWallError {
    fn from(error: PhaseStateError) -> Self {
        Self::Wall(error)
    }
}
impl From<DepositionError> for WaterWallError {
    fn from(error: DepositionError) -> Self {
        match error {
            DepositionError::Cancelled => Self::Cancelled,
            other => Self::Conduction(other),
        }
    }
}
impl core::fmt::Display for WaterWallError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite wall/water thermal exchange: {self:?}")
    }
}
impl std::error::Error for WaterWallError {}

/// Native entropy changes and the endpoint conductive entropy account.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterWallEntropy {
    /// Exact integral M integral(dh/T) on the original wall chart.
    pub wall_change_j_k: f64,
    /// Original IF97 endpoint entropy change, including its source seams.
    pub water_source_change_j_k: f64,
    /// dt * P_end * (1/T_water,end - 1/T_wall,end), nonnegative.
    pub conduction_generation_j_k: f64,
    /// Q * (1/T_water,end - 1/T_wall,end).
    pub transferred_interface_entropy_j_k: f64,
    /// Transferred interface entropy minus native conduction generation.
    /// This is the exposed heat-root error times the reciprocal-T difference.
    pub interface_source_residual_j_k: f64,
    /// Sum of original storage entropy changes minus the transferred interface
    /// entropy. Contains time-discretization/source defects and any unspecified
    /// external heat entropy; it is not labelled physical entropy production.
    pub unaccounted_change_j_k: f64,
}

/// Accepted original owner responses and separate coupled physical accounts.
#[derive(Clone, Debug, PartialEq)]
pub struct WaterWallStep {
    pub initial_wall: EquilibriumPhaseState,
    pub water: WaterCellStep,
    /// Original wall endpoint and native conducting-layer response. Boundary
    /// heat into the wall is exactly -wall_to_water_heat_j.
    pub wall: UniformEnthalpyStep<LayerStackResponse>,
    pub wall_to_water_heat_j: f64,
    /// Q - dt times the original stack's endpoint power.
    pub source_residual_j: f64,
    pub wall_energy_change_j: f64,
    /// Wall storage + Q - authored internal wall heat.
    pub wall_energy_residual_j: f64,
    /// Conserved water delta U + wall delta H - the two external heat inputs.
    pub combined_energy_residual_j: f64,
    /// Same account using the actual recovered IF97 source delta U.
    /// The difference from the conserved residual is the change in the two
    /// exposed water source-inversion defects; no error is relabelled as heat.
    pub source_combined_energy_residual_j: f64,
    pub entropy: WaterWallEntropy,
    /// Complete coupled proposals, at most exchange.maximum_iterations + 2.
    pub trials: usize,
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), WaterWallError> {
    if cancelled() {
        Err(WaterWallError::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(value: f64) -> Result<f64, WaterWallError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(WaterWallError::Invalid(
            "nonfinite thermal exchange arithmetic",
        ))
    }
}
fn energy_budget(
    quantity: &'static str,
    residual_j: f64,
    tolerance_j: f64,
) -> Result<(), WaterWallError> {
    if residual_j.is_finite() && residual_j.abs() <= tolerance_j {
        Ok(())
    } else {
        Err(WaterWallError::Budget {
            quantity,
            residual_j,
            tolerance_j,
        })
    }
}
fn temperature_budget(
    body: &'static str,
    increment_k: f64,
    maximum_k: f64,
) -> Result<(), WaterWallError> {
    if increment_k.is_finite() && increment_k.abs() <= maximum_k {
        Ok(())
    } else {
        Err(WaterWallError::TemperatureIncrement {
            body,
            increment_k,
            maximum_k,
        })
    }
}

/// Advance one finite water cell and one original uniform wall state together.
///
/// The wall's existing `UniformEnthalpyStepInput` supplies immutable calorics,
/// mass, geometry and external internal heating. Its `initial` must equal
/// `*wall`, and its duration must be the actual cell-clock increment. The
/// returned original wall state can therefore continue any existing uniform
/// thermal consumer without introducing another wall state or EOS.
///
/// `water_heat_j` is prescribed heat into the water EXCLUDING the solved
/// interface transfer. The wall's `internal_heat_j` likewise excludes that
/// transfer. Neither history changes on cancellation, a source-domain refusal,
/// a missing root, or a failed conservation/temperature budget. Publication
/// occurs after the final cancellation poll; callback side effects are not
/// rolled back. Both source laws remain immutable during the solve.
///
/// Time and temperature ceilings are admission bounds. Refine the interval to
/// establish transient accuracy; this solve does not resolve wall diffusion,
/// liquid mixing, nucleation, or boiling-film resistance.
#[allow(clippy::too_many_arguments)]
pub fn advance_water_wall(
    water: &mut EquilibriumWaterCell,
    wall: &mut EquilibriumPhaseState,
    body: UniformEnthalpyStepInput<'_>,
    layers: LayerStackConduction<'_>,
    water_heat_j: f64,
    end_time_s: f64,
    controls: WaterWallControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<WaterWallStep, WaterWallError> {
    poll(&mut cancelled)?;
    let duration = end_time_s - water.state().time_s();
    let [lo, hi] = controls.exchange.heat_bracket_j;
    if !positive(duration)
        || !end_time_s.is_finite()
        || duration != body.duration_s
        || duration > controls.maximum_step_s
        || !positive(controls.maximum_step_s)
        || !controls
            .maximum_temperature_increment_k
            .into_iter()
            .all(positive)
        || !positive(body.mass_kg)
        || !positive(body.volume_m3)
        || !positive(body.surface_area_m2)
        || !body.internal_heat_j.is_finite()
        || !water_heat_j.is_finite()
        || body.initial != *wall
        || !lo.is_finite()
        || !hi.is_finite()
        || lo > hi
        || lo > 0.0
        || hi < 0.0
        || !positive(controls.exchange.source_tolerance_j)
        || !positive(controls.exchange.energy_tolerance_j)
        || !(1..=128).contains(&controls.exchange.maximum_iterations)
        || !(1..=128).contains(&controls.conduction.max_iterations)
        || !controls.conduction.absolute_heat_tolerance_w.is_finite()
        || controls.conduction.absolute_heat_tolerance_w < 0.0
        || !controls.conduction.relative_heat_tolerance.is_finite()
        || controls.conduction.relative_heat_tolerance < 0.0
    {
        return Err(WaterWallError::Invalid(
            "matching wall/clock, supported heat bracket, finite sources and bounded controls required",
        ));
    }
    let h0 = body.initial.specific_enthalpy_j_kg();
    if body.curve.state_at_specific_enthalpy(h0)? != body.initial {
        return Err(WaterWallError::Invalid(
            "wall must retain its exact original enthalpy/phase chart",
        ));
    }
    let initial_water_temperature = water.state().source_state().temperature_k();
    let step = feedback::solve_scalar(
        controls.exchange.heat_bracket_j,
        controls.exchange.source_tolerance_j,
        controls.exchange.maximum_iterations,
        |q, trials, cancelled| {
            let wall_heat = finite(body.internal_heat_j - q)?;
            let water_heat = finite(water_heat_j + q)?;
            let h1 = finite(h0 + wall_heat / body.mass_kg)?;
            if (wall_heat != 0.0 && h1 == h0)
                || (q != 0.0 && (wall_heat == body.internal_heat_j || water_heat == water_heat_j))
            {
                return Err(WaterWallError::Invalid(
                    "interface heat disappears in a conserved energy coordinate",
                ));
            }
            let wall_state = body.curve.state_at_specific_enthalpy(h1)?;
            temperature_budget(
                "wall",
                wall_state.temperature_k() - body.initial.temperature_k(),
                controls.maximum_temperature_increment_k[0],
            )?;
            // Reuse the same native proposal as ordinary water-cell transport;
            // all trials start from the accepted conserved M,U and fixed V.
            let water_step = water.propose(&[], 0.0, water_heat, end_time_s, cancelled)?;
            let water_temperature = water_step.state.source_state().temperature_k();
            temperature_budget(
                "water",
                water_temperature - initial_water_temperature,
                controls.maximum_temperature_increment_k[1],
            )?;
            poll(cancelled)?;
            let conduction = layers.evaluate(
                [wall_state.temperature_k(), water_temperature],
                controls.conduction,
                || {
                    if cancelled() {
                        Err(DepositionError::Cancelled)
                    } else {
                        Ok(())
                    }
                },
            )?;
            let source_heat = finite(duration * conduction.power_a_to_b_w)?;
            if conduction.power_a_to_b_w != 0.0 && source_heat == 0.0 {
                return Err(WaterWallError::Invalid("conductive heat underflows"));
            }
            let source_residual = finite(q - source_heat)?;
            // Layer heat-continuity error is a separate native solve residual.
            energy_budget(
                "conducting-layer heat continuity",
                duration * conduction.max_layer_heat_residual_w,
                controls.exchange.source_tolerance_j,
            )?;
            let stored_wall = finite(body.mass_kg * (h1 - h0))?;
            let wall_residual = finite((stored_wall + q) - body.internal_heat_j)?;
            let combined_residual = finite(
                (water_step.internal_energy_change_j + stored_wall)
                    - water_heat_j
                    - body.internal_heat_j,
            )?;
            energy_budget(
                "wall storage",
                wall_residual,
                controls.exchange.energy_tolerance_j,
            )?;
            energy_budget(
                "combined conserved storage",
                combined_residual,
                controls.exchange.energy_tolerance_j,
            )?;
            let source_combined = finite(water_step.source_first_law_residual_j + wall_residual)?;
            let wall_entropy =
                finite(body.mass_kg * body.curve.specific_entropy_change_j_kg_k(h0, h1)?)?;
            let reciprocal_difference = (wall_state.temperature_k() - water_temperature)
                / wall_state.temperature_k()
                / water_temperature;
            let conduction_entropy = finite(source_heat * reciprocal_difference)?;
            if conduction_entropy < 0.0 {
                return Err(WaterWallError::Invalid(
                    "native conduction has negative entropy generation",
                ));
            }
            let transferred_entropy = finite(q * reciprocal_difference)?;
            let entropy = WaterWallEntropy {
                wall_change_j_k: wall_entropy,
                water_source_change_j_k: water_step.source_entropy_change_j_k,
                conduction_generation_j_k: conduction_entropy,
                transferred_interface_entropy_j_k: transferred_entropy,
                interface_source_residual_j_k: finite(transferred_entropy - conduction_entropy)?,
                unaccounted_change_j_k: finite(
                    (wall_entropy + water_step.source_entropy_change_j_k) - transferred_entropy,
                )?,
            };
            poll(cancelled)?;
            Ok((
                WaterWallStep {
                    initial_wall: body.initial,
                    water: water_step,
                    wall: UniformEnthalpyStep {
                        state: wall_state,
                        external_heat_j: -q,
                        energy_residual_tolerance_j: controls.exchange.energy_tolerance_j,
                        report: conduction,
                    },
                    wall_to_water_heat_j: q,
                    source_residual_j: source_residual,
                    wall_energy_change_j: stored_wall,
                    wall_energy_residual_j: wall_residual,
                    combined_energy_residual_j: combined_residual,
                    source_combined_energy_residual_j: source_combined,
                    entropy,
                    trials,
                },
                source_residual,
            ))
        },
        &mut cancelled,
    )
    .map_err(|error| match error {
        feedback::ScalarError::Evaluation(error) => error,
        feedback::ScalarError::Invalid(reason) => WaterWallError::Invalid(reason),
        feedback::ScalarError::Cancelled => WaterWallError::Cancelled,
        feedback::ScalarError::OutsideBracket { lower, upper } => WaterWallError::OutsideBracket {
            lower_residual_j: lower,
            upper_residual_j: upper,
        },
        feedback::ScalarError::NotConverged {
            iterations,
            residual,
        } => WaterWallError::NotConverged {
            iterations,
            residual_j: residual,
        },
    })?;
    poll(&mut cancelled)?;
    // All work and validation is finished; nothing fallible follows publication.
    water.state = step.water.state;
    *wall = step.wall.state;
    Ok(step)
}

#[cfg(test)]
#[path = "thermal/tests.rs"]
mod tests;
