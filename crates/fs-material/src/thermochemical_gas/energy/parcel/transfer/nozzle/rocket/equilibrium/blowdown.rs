//! Finite reacting reservoirs and physical-time rocket discharge.
//!
//! The reusable reservoir leaf owns the immutable birth isentrope and actual
//! finite mass/energy exchange. Rocket and higher-layer converging restrictions
//! use that same owner and the shared bounded midpoint mass equation.

use super::*;
use reservoir::check_cancel;

mod reservoir;
pub use reservoir::{
    EquilibriumBlowdownCandidate, EquilibriumBlowdownError, EquilibriumBlowdownSolve,
    EquilibriumBlowdownStep, EquilibriumGasReservoir,
};
/// Shared bounded midpoint mass equation for original finite-outflow consumers.
pub mod midpoint;
/// Reuse the original reservoir/density/finite-energy owner across flow laws.
pub mod prepared;

/// Physical-time equilibrium blowdown on one reusable birth isentrope.
pub mod time {
//! Implicit midpoint closes dm=dt*mdot(m0-dm/2). At an ordinary sonic maximum,
//! d(mdot)/dh_chamber=mdot/a_throat^2 and dh_chamber/dm=a_chamber^2/m.
//! These reacting derivatives supply a safeguarded analytic Newton slope.
//! Mass trials solve only the throat; full initial/midpoint/final flow admission
//! precedes atomic physical publication.

use super::*;

/// Physical-time work and outgoing-inventory limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquilibriumRocketBlowdownSolve {
    /// Original density and extensive-inventory admission.
    pub reservoir: EquilibriumBlowdownSolve,
    /// Reacting nozzle controls; thermochemistry must match tank birth.
    pub nozzle: EquilibriumRocketNozzleSolve,
    /// Maximum fraction of OLD tank mass expelled, in (0,0.5]. No clipping.
    pub maximum_outgoing_fraction: f64,
    /// At most 128 mass iterations, each with bounded density/sonic roots.
    pub maximum_mass_iterations: usize,
}

/// Failure returns neither a new inventory nor an accepted impulse.
#[derive(Clone, Debug, PartialEq)]
pub enum EquilibriumRocketBlowdownError {
    /// Original finite reservoir refusal.
    Reservoir(EquilibriumBlowdownError),
    /// Original nozzle/source/regime refusal.
    Nozzle(EquilibriumRocketNozzleError),
    /// Invalid duration, control mismatch or unrepresentable exchange.
    Invalid(&'static str),
    /// Requested interval exceeds the admitted mass fraction.
    InventoryLimit { maximum_duration_s: f64 },
    /// Actual mass-rate equation failed to close within the budget.
    NoConvergence { mass_iterations: usize, mass_residual_kg: f64 },
}
impl From<EquilibriumBlowdownError> for EquilibriumRocketBlowdownError {
    fn from(value: EquilibriumBlowdownError) -> Self { Self::Reservoir(value) }
}
impl From<EquilibriumExpansionErrorV1> for EquilibriumRocketBlowdownError {
    fn from(value: EquilibriumExpansionErrorV1) -> Self {
        Self::Reservoir(EquilibriumBlowdownError::Thermochemistry(value))
    }
}
impl From<EquilibriumRocketNozzleError> for EquilibriumRocketBlowdownError {
    fn from(value: EquilibriumRocketNozzleError) -> Self { Self::Nozzle(value) }
}
impl core::fmt::Display for EquilibriumRocketBlowdownError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "equilibrium rocket blowdown: {self:?}")
    }
}
impl std::error::Error for EquilibriumRocketBlowdownError {}
type AdvanceError = EquilibriumRocketBlowdownError;

/// Accepted physical interval with separate thermal and propulsive accounts.
#[derive(Clone, Debug, PartialEq)]
pub struct EquilibriumRocketBlowdownStep {
    /// Positive shared physical duration [s].
    pub duration_s: f64,
    /// Actual reservoir mass and integrated exported stagnation enthalpy.
    pub gas: EquilibriumBlowdownStep,
    /// Full source-bearing midpoint chamber/throat/exit and performance.
    pub midpoint_nozzle: EquilibriumRocketNozzleFlow,
    /// Midpoint total thrust times duration [N s], including pressure thrust once.
    pub axial_impulse_n_s: f64,
    /// Actual expelled mass minus dt*accepted_midpoint_mdot [kg].
    pub mass_quadrature_residual_kg: f64,
    /// Number of bounded implicit mass iterations.
    pub mass_iterations: usize,
    /// Pressure queries in all density inversions, including final reconstruction.
    pub density_pressure_queries: usize,
    /// SP queries in all initial/trial/midpoint/final nozzle solves.
    pub nozzle_pressure_queries: usize,
}
/// Reservoir and propulsive exchange are accepted together.
#[derive(Clone, Debug, PartialEq)]
pub struct EquilibriumRocketBlowdownCandidate {
    endpoint: EquilibriumGasReservoir,
    step: EquilibriumRocketBlowdownStep,
}
impl EquilibriumRocketBlowdownCandidate {
    /// Original-source remaining gas; input reservoir is unchanged.
    pub const fn endpoint(&self) -> &EquilibriumGasReservoir { &self.endpoint }
    /// Same interval's actual mass, energy and midpoint impulse.
    pub const fn step(&self) -> &EquilibriumRocketBlowdownStep { &self.step }
    /// Accept without another thermodynamic or nozzle solve.
    pub fn into_parts(self) -> (EquilibriumGasReservoir, EquilibriumRocketBlowdownStep) {
        (self.endpoint, self.step)
    }
}

impl EquilibriumGasReservoir {
    /// Advance finite equilibrium propellant with source-derived nozzle thrust.
    /// The reservoir and nozzle share ONE prepared birth isentrope and source
    /// catalog throughout the interval. Smooth time accuracy is second order;
    /// the outgoing-fraction guard is NOT a temporal error estimate.
    /// Initial, midpoint and endpoint nozzle regimes must all be admitted.
    /// Any late refusal leaves the old reservoir untouched. No valve, wall heat,
    /// injector, finite-rate chemistry or exhaust species quadrature is inferred.
    pub fn rocket_blowdown(&self, ambient_pressure_pa: f64, nozzle: RocketNozzle,
        duration_s: f64, c: EquilibriumRocketBlowdownSolve, mut cancelled: impl FnMut() -> bool)
        -> Result<EquilibriumRocketBlowdownCandidate, AdvanceError> {
        check_cancel(&mut cancelled)?;
        c.reservoir.validate(self.isentropic_sample().pressure.value())?;
        validate_inputs(self.isentropic_sample().pressure.value(), ambient_pressure_pa, nozzle, c.nozzle)?;
        if !positive(duration_s) || !positive(c.maximum_outgoing_fraction)
            || c.maximum_outgoing_fraction > 0.5 || !(1..=128).contains(&c.maximum_mass_iterations)
            || c.nozzle.thermochemistry != self.thermochemical_controls() {
            return Err(AdvanceError::Invalid("invalid time/inventory controls or changed birth thermochemistry"));
        }
        let mut path = self.path(&mut cancelled)?;
        let mut work = Work { composition: path.reference_state().equilibrium().composition_iterations(),
            ..Work::default() };
        let initial = flow_at(&mut path, self.thermodynamics().clone(), self.isentropic_sample().specific_enthalpy_change_j_kg,
            ambient_pressure_pa, nozzle, c.nozzle, &mut work, &mut cancelled)?;
        let mut density_queries = 0;
        let quadrature = midpoint::solve(self.mass_kg(), initial.mass_flow_kg_s,
            midpoint::OutflowMidpointSolve { duration_s,
                maximum_outgoing_fraction: c.maximum_outgoing_fraction,
                mass_tolerance_kg: c.reservoir.mass_tolerance_kg,
                maximum_iterations: c.maximum_mass_iterations },
            |midpoint_mass, mut cancelled| -> Result<midpoint::OutflowRate, AdvanceError> {
                let query = self.density_at_mass(&mut path, midpoint_mass, c.reservoir, &mut cancelled)?;
                density_queries += query.pressure_queries;
                work.composition += query.composition_iterations;
                work.temperature += query.temperature_iterations;
                let chamber = query.sample;
                let (_, throat) = sonic_throat(&mut path, chamber.pressure.value(),
                    chamber.specific_enthalpy_change_j_kg, c.nozzle, &mut work, &mut cancelled)?;
                let rate = throat.flux * nozzle.throat_area_m2;
                if !positive(rate) {
                    return Err(AdvanceError::Invalid("unrepresentable midpoint mass rate"));
                }
                let ratio = chamber.response.sound_speed_m_s / throat.sample.response.sound_speed_m_s;
                let derivative = (rate / midpoint_mass) * ratio * ratio;
                Ok(midpoint::OutflowRate { mass_flow_kg_s: rate,
                    derivative_per_s: derivative.is_finite().then_some(derivative) })
            }, &mut cancelled).map_err(|error| match error {
                midpoint::OutflowMidpointError::Invalid(reason) => AdvanceError::Invalid(reason),
                midpoint::OutflowMidpointError::Source(source) => source,
                midpoint::OutflowMidpointError::InventoryLimit { maximum_duration_s } =>
                    AdvanceError::InventoryLimit { maximum_duration_s },
                midpoint::OutflowMidpointError::NoConvergence { iterations, mass_residual_kg } =>
                    AdvanceError::NoConvergence { mass_iterations: iterations, mass_residual_kg },
                midpoint::OutflowMidpointError::Cancelled =>
                    AdvanceError::Reservoir(EquilibriumBlowdownError::Cancelled),
            })?;
        let remaining = quadrature.remaining_mass_kg;
        let mass_iterations = quadrature.iterations;
        let actual = self.mass_kg() - remaining;
        let midpoint_mass = quadrature.midpoint_mass_kg;
        let midpoint = self.density_at_mass(&mut path, midpoint_mass, c.reservoir, &mut cancelled)?;
        density_queries += midpoint.pressure_queries;
        work.composition += midpoint.composition_iterations;
        work.temperature += midpoint.temperature_iterations;
        let chamber = path.publish_current(&mut cancelled)?;
        let midpoint_nozzle = flow_at(&mut path, chamber, midpoint.sample.specific_enthalpy_change_j_kg,
            ambient_pressure_pa, nozzle, c.nozzle, &mut work, &mut cancelled)?;
        let residual = actual - duration_s * midpoint_nozzle.mass_flow_kg_s;
        if !residual.is_finite() || residual.abs() > c.reservoir.mass_tolerance_kg {
            return Err(AdvanceError::NoConvergence { mass_iterations, mass_residual_kg: residual });
        }
        let candidate = self.candidate_at_mass(&mut path, remaining, c.reservoir, &mut cancelled)?;
        density_queries += candidate.step().pressure_queries;
        flow_at(&mut path, candidate.endpoint().thermodynamics().clone(),
            candidate.endpoint().isentropic_sample().specific_enthalpy_change_j_kg,
            ambient_pressure_pa, nozzle, c.nozzle, &mut work, &mut cancelled)?;
        let impulse = duration_s * midpoint_nozzle.thrust_n;
        if !positive(impulse) { return Err(AdvanceError::Invalid("unrepresentable nozzle impulse")); }
        let (endpoint, gas) = candidate.into_parts();
        check_cancel(&mut cancelled)?;
        Ok(EquilibriumRocketBlowdownCandidate { endpoint, step: EquilibriumRocketBlowdownStep {
            duration_s, gas, midpoint_nozzle, axial_impulse_n_s: impulse,
            mass_quadrature_residual_kg: residual, mass_iterations,
            density_pressure_queries: density_queries, nozzle_pressure_queries: work.samples,
        } })
    }
}
}
