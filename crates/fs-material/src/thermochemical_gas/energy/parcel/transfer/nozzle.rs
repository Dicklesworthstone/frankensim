//! Pressure-driven mass flow through a quasi-steady converging restriction.
//!
//! Use the original species h(T), s(T,p) and frozen acoustic c(T). At the exit,
//! h0=h+v^2/2 and s=s0. Lowering back pressure reaches the sonic state v=c;
//! the converging throat then remains choked, instead of using an incompressible
//! Bernoulli rate or freezing gamma at the upstream temperature. No new caloric
//! coefficients, species registry or fixed-air constants are introduced.
//!
//! The isentropic nozzle is an ideal inviscid model with an explicit discharge
//! coefficient, not a viscous valve calibration. It assumes a continuous,
//! ordinary gas branch with one sonic crossing over the admitted temperature
//! bracket. It models neither supersonic diverging geometry nor shocks. The
//! receiving well-mixed tank thermalizes the jet's total enthalpy; that kinetic
//! energy must NOT be deposited again as extra heat.
//! Reference: NASA Glenn, Isentropic Flow Equation Derivations / Mass Flow Choking.
//! https://www.grc.nasa.gov/www/BGH/isndrv.html
//! https://www.grc.nasa.gov/www/BGH/mflchk.html

/// Conservative pressure-driven exchange on a graph of existing gas parcels.
pub mod network;
/// Frozen-composition supersonic expansion and physical rocket thrust.
pub mod rocket;
/// Conservative finite-reservoir exhaust and physical-time blowdown.
pub mod blowdown;

use super::{ClosedGasParcel, GasEnergyError, GasEnergySolve, GasMassTransferCandidate,
    GasMassTransferSolve, invalid, poll, positive, source_error, updated_point};
use crate::thermochemical_gas::ResolvedThermochemicalGasState;
use fs_math::det;
use fs_qty::{Pressure, Temperature};

/// Explicit effective-flow geometry. Zero area or coefficient closes the path.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasNozzle {
    /// Minimum physical flow area [m2].
    pub area_m2: f64,
    /// Caller-declared coefficient in [0,1]; one is the ideal nozzle.
    pub discharge_coefficient: f64,
}
impl GasNozzle {
    fn validate(self) -> Result<(), GasEnergyError> {
        if !self.area_m2.is_finite() || self.area_m2 < 0.0
            || !self.discharge_coefficient.is_finite()
            || !(0.0..=1.0).contains(&self.discharge_coefficient) {
            return Err(invalid("finite nonnegative nozzle area and discharge coefficient in [0,1] required"));
        }
        Ok(())
    }
}

/// A source-resolved instantaneous nozzle flow, before any inventory update.
#[derive(Clone, Debug, PartialEq)]
pub struct GasNozzleFlow {
    throat: ResolvedThermochemicalGasState,
    /// Positive mass rate from upstream reservoir toward the receiver [kg/s].
    pub mass_flow_kg_s: f64,
    /// Ideal gas speed at the converging exit/throat [m/s].
    pub velocity_m_s: f64,
    /// Velocity divided by the SAME throat material sound speed.
    pub mach_number: f64,
    /// Whether the sonic condition, rather than back pressure, sets the exit.
    pub choked: bool,
    /// Original stagnation enthalpy [J/kg], including the source reference.
    pub upstream_enthalpy_j_kg: f64,
    /// Actual exit entropy error times upstream Kelvin [J/kg].
    pub entropy_energy_residual_j_kg: f64,
    /// Interior temperature samples, at most the specified iteration cap.
    pub iterations: usize,
}
impl GasNozzleFlow {
    /// Complete flow-state medium, including original material/pack receipts.
    #[must_use]
    pub const fn throat(&self) -> &ResolvedThermochemicalGasState { &self.throat }
}

/// Work and physical-step bounds, separate from numerical root convergence.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasNozzleStepSolve {
    /// Bracket and absolute J/kg tolerance for the physical nozzle residual.
    pub nozzle: GasEnergySolve,
    /// Both endpoint energy inversions and extensive conservation tolerances.
    pub transfer: GasMassTransferSolve,
    /// Largest accepted fraction of upstream inventory in one explicit step.
    /// Must lie in (0,1); this bound is not a time-error estimate.
    pub maximum_donor_fraction: f64,
}

struct Sample {
    pressure_residual: f64,
    sonic_residual: f64,
    entropy_shift: f64,
}
impl Sample {
    fn residual(&self) -> f64 { self.pressure_residual.min(self.sonic_residual) }
}

impl ClosedGasParcel {
    /// Resolve an ideal converging-nozzle flow at the unchanged upstream state.
    /// Supply an explicit lower-pressure boundary and geometry. Vacuum is a
    /// valid back-pressure boundary; the actual sonic throat stays positive.
    /// Reverse flow
    /// requires swapping donor/receiver, not extrapolating this upwind law.
    ///
    /// The root is min(T0*(s(T,p0)-s0-R*ln(pb/p0)), h(T)-h0+c(T)^2/2).
    /// Both terms have units J/kg. Its highest-temperature zero is the first
    /// restriction encountered from the reservoir: back-pressure matching OR
    /// sonic choking. Bisection uses only the existing material evaluations;
    /// no derivative of a second cp/gamma approximation is required.
    ///
    /// A lower bracket need only reach the active exit condition. Small
    /// pressure drops do not require transport data at a hypothetical sonic
    /// temperature. Endpoint transport and entropy are checked at the actual
    /// accepted exit. Work exhaustion and source gaps are refusals, not clamps.
    pub fn nozzle_flow_to(
        &self, back_pressure_pa: f64, nozzle: GasNozzle, controls: GasEnergySolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasNozzleFlow, GasEnergyError> {
        controls.validate()?;
        nozzle.validate()?;
        poll(&mut cancelled)?;
        let gas = self.state().state();
        let (t0, p0, r) = (gas.temperature, gas.pressure, gas.specific_gas_constant);
        if !back_pressure_pa.is_finite() || back_pressure_pa < 0.0 || back_pressure_pa > p0
            || t0 < controls.minimum_temperature_k || t0 > controls.maximum_temperature_k {
            return Err(invalid("nozzle requires nonnegative finite downstream pressure no higher than upstream and an admitted upstream temperature"));
        }
        let h0 = self.state().specific_enthalpy_j_kg();
        if back_pressure_pa == p0 || nozzle.area_m2 == 0.0 || nozzle.discharge_coefficient == 0.0 {
            poll(&mut cancelled)?;
            return Ok(GasNozzleFlow { throat: self.state().clone(), mass_flow_kg_s: 0.0,
                velocity_m_s: 0.0, mach_number: 0.0, choked: false,
                upstream_enthalpy_j_kg: h0, entropy_energy_residual_j_kg: 0.0, iterations: 0 });
        }
        let thermodynamics = self.parcel.model.thermodynamics();
        // Entropy references and ideal mixing cancel BEFORE rounding. This also
        // preserves the pressure signal for very small reservoir pressure drops.
        // Vacuum is a boundary condition, never an EOS query at p=0. It can
        // only reach the sonic constraint in a converging nozzle.
        let log_back = if back_pressure_pa == 0.0 { None } else {
            Some(-thermodynamics.specific_entropy_change(
                Temperature::new(t0), Pressure::new(p0),
                Temperature::new(t0), Pressure::new(back_pressure_pa))
                .map_err(source_error)?.value()/r)
        };
        let sample = |t: f64| -> Result<Sample, GasEnergyError> {
            let properties = thermodynamics.evaluate_properties(Temperature::new(t), Pressure::new(p0))
                .map_err(source_error)?;
            let cp = properties.mass_specific_cp().value();
            let cv = properties.mass_specific_cv().value();
            let gamma = cp / cv;
            if !positive(cp) || !positive(cv) || !gamma.is_finite()
                || gamma <= 1.0 || gamma > 5.0/3.0+1e-9 {
                return Err(invalid("nozzle caloric sample has inadmissible heat capacities"));
            }
            let entropy_shift = thermodynamics.specific_entropy_change(
                Temperature::new(t0), Pressure::new(p0),
                Temperature::new(t), Pressure::new(p0))
                .map_err(source_error)?.value()/r;
            let enthalpy_change = thermodynamics.specific_enthalpy_change(
                Temperature::new(t0), Temperature::new(t))
                .map_err(source_error)?.value();
            let sonic_residual = enthalpy_change+0.5*gamma*r*t;
            let pressure_residual = log_back.map_or(f64::INFINITY, |back| r*t0*(entropy_shift-back));
            if !entropy_shift.is_finite() || !sonic_residual.is_finite()
                || (log_back.is_some() && !pressure_residual.is_finite()) {
                return Err(invalid("nozzle thermal residual is unrepresentable"));
            }
            Ok(Sample { pressure_residual, sonic_residual, entropy_shift })
        };
        // The sonic derivative depends on the original temperature-varying
        // caloric law. A bounded derivative-free bracket avoids inventing it.
        let (mut lo, mut hi) = (controls.minimum_temperature_k, t0);
        poll(&mut cancelled)?;
        let lower_value = sample(lo)?;
        let lower = lower_value.residual();
        poll(&mut cancelled)?;
        let upper = sample(hi)?.residual();
        if lower > controls.energy_tolerance_j_kg || upper < -controls.energy_tolerance_j_kg {
            return Err(GasEnergyError::OutsideBracket { lower_residual_j_kg: lower, upper_residual_j_kg: upper });
        }
        let mut temperature = lo;
        let mut value = lower_value;
        let mut iterations = 0;
        while value.residual().abs() > controls.energy_tolerance_j_kg {
            if iterations == controls.maximum_iterations {
                return Err(GasEnergyError::NotConverged { iterations, residual_j_kg: value.residual() });
            }
            let next = lo+0.5*(hi-lo);
            if next <= lo || next >= hi {
                return Err(GasEnergyError::NotConverged { iterations, residual_j_kg: value.residual() });
            }
            temperature = next;
            poll(&mut cancelled)?;
            value = sample(temperature)?;
            iterations += 1;
            if value.residual() < 0.0 { lo = temperature; } else { hi = temperature; }
        }
        let choked = value.sonic_residual <= value.pressure_residual;
        let pressure = if choked { p0*det::exp(value.entropy_shift) } else { back_pressure_pa };
        poll(&mut cancelled)?;
        let throat = self.parcel.model.evaluate(&updated_point(self.point(), temperature, pressure)?)?;
        // Formation references cancel in the shared caloric owner BEFORE
        // rounding, just as for entropy. No reference energy accelerates a jet.
        let speed_squared = -2.0*thermodynamics.specific_enthalpy_change(
            Temperature::new(t0), Temperature::new(temperature))
            .map_err(source_error)?.value();
        let c = throat.state().sound_speed;
        let sonic = 0.5*(c*c-speed_squared);
        let entropy_residual = t0*thermodynamics.specific_entropy_change(
            Temperature::new(t0), Pressure::new(p0),
            Temperature::new(temperature), Pressure::new(pressure)).map_err(source_error)?.value();
        if !positive(speed_squared) || !sonic.is_finite() || !entropy_residual.is_finite()
            || entropy_residual.abs() > 2.0*controls.energy_tolerance_j_kg
            || (choked && sonic.abs() > 2.0*controls.energy_tolerance_j_kg)
            || (!choked && sonic < -2.0*controls.energy_tolerance_j_kg) {
            return Err(invalid("resolved nozzle state failed entropy, energy or sonic admission"));
        }
        let velocity = det::sqrt(speed_squared);
        let rate = (nozzle.discharge_coefficient*nozzle.area_m2)*throat.state().density*velocity;
        let mach = velocity/c;
        if !positive(rate) || !positive(mach) {
            return Err(invalid("nonzero nozzle mass flow is unrepresentable"));
        }
        poll(&mut cancelled)?;
        Ok(GasNozzleFlow { throat, mass_flow_kg_s: rate, velocity_m_s: velocity, mach_number: mach,
            choked, upstream_enthalpy_j_kg: h0, entropy_energy_residual_j_kg: entropy_residual, iterations })
    }

    /// Transfer a pressure-driven mass over an explicit physical duration.
    /// All rate, species and energy calculations start from the SAME accepted
    /// upstream state. Return neither endpoint on source/budget/cancellation
    /// failure, excessive withdrawal or pressure-order reversal during the step.
    /// No mass is clipped to fit a bound: shorten/refine the caller's timestep.
    ///
    /// This is first-order quasi-steady tank exchange, not valve-actuator or
    /// resolved pipe/jet dynamics. Volumes stay fixed during transfer; compose
    /// subsequent heat/motion steps explicitly rather than double-depositing
    /// enthalpy, pressure work or jet kinetic energy.
    pub fn through_nozzle_to(
        &self, receiver: &Self, nozzle: GasNozzle, duration_s: f64,
        controls: GasNozzleStepSolve, mut cancelled: impl FnMut() -> bool,
    ) -> Result<(GasMassTransferCandidate, GasNozzleFlow), GasEnergyError> {
        if !positive(duration_s) || !positive(controls.maximum_donor_fraction)
            || controls.maximum_donor_fraction >= 1.0 {
            return Err(invalid("positive nozzle step duration and donor fraction in (0,1) required"));
        }
        controls.transfer.validate()?;
        let flow = self.nozzle_flow_to(receiver.state().state().pressure, nozzle, controls.nozzle, &mut cancelled)?;
        let mass = flow.mass_flow_kg_s*duration_s;
        if !mass.is_finite() || (flow.mass_flow_kg_s > 0.0 && mass == 0.0)
            || mass/self.mass_kg() > controls.maximum_donor_fraction {
            return Err(invalid("nozzle timestep exceeds the donor inventory bound or mass resolution"));
        }
        let candidate = self.transfer_mass_to(receiver, mass, controls.transfer, &mut cancelled)?;
        if candidate.donor().state().state().pressure < candidate.receiver().state().state().pressure {
            return Err(invalid("nozzle timestep crosses pressure equilibrium; reduce the physical step"));
        }
        poll(&mut cancelled)?;
        Ok((candidate, flow))
    }
}

#[cfg(test)]
mod tests;
