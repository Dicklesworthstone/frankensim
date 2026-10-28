//! Frozen-composition, thermally perfect converging-diverging rocket nozzles.
//!
//! The existing sonic closure sets the throat mass flux. The diverging branch
//! solves rho(T) v(T) A_exit = mdot using the SAME NASA-9 entropy and enthalpy
//! differences, with v²/2 = h(T0)-h(T). No upstream-gamma substitution, fitted
//! chamber temperature, or independently prescribed incompatible mass rate.
//!
//! This is an ideal, quasi-one-dimensional, adiabatic, inviscid nozzle with an
//! ordinary continuous gas branch and one sonic crossing. The present lane
//! admits pressure-matched and underexpanded exhaust, including vacuum ambient.
//! Overexpanded operation requires a shock/separation model and is refused.
//! Frozen composition does not imply combustion equilibrium or real-gas validity.
//!
//! References: NASA Glenn, Rocket Thrust Equations and Isentropic Flow Equations:
//! <https://www.grc.nasa.gov/www/BGH/rktthsum.html>
//! <https://www.grc.nasa.gov/www/BGH/isentrop.html>

use super::{
    ClosedGasParcel, GasEnergyError, GasEnergySolve, GasNozzle, ResolvedThermochemicalGasState,
    det, invalid, poll, positive, source_error, updated_point,
};
use fs_qty::{Pressure, Temperature};

/// Standard gravity used only to express specific impulse in seconds [m/s²].
/// It is not an assumed gravitational acceleration acting on a vehicle.
pub const STANDARD_GRAVITY_M_S2: f64 = 9.80665;

/// Physical flow areas of an ideal rocket nozzle, in square metres.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RocketNozzle {
    /// Sonic minimum area, strictly positive [m²].
    pub throat_area_m2: f64,
    /// Exit area, at least the throat area [m²].
    pub exit_area_m2: f64,
}
impl RocketNozzle {
    fn validate(self) -> Result<(), GasEnergyError> {
        if !positive(self.throat_area_m2)
            || !positive(self.exit_area_m2)
            || self.exit_area_m2 < self.throat_area_m2
            || !(self.exit_area_m2 / self.throat_area_m2).is_finite()
        {
            return Err(invalid(
                "rocket nozzle requires finite positive areas with exit area at least throat area",
            ));
        }
        Ok(())
    }
}

/// Explicit source-temperature, residual and root-work limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RocketNozzleSolve {
    /// Source bracket and absolute energy/entropy-energy tolerance [J/kg].
    /// The lower temperature must reach the requested supersonic exit branch.
    /// The iteration cap applies separately to the throat and exit solves.
    pub nozzle: GasEnergySolve,
    /// Maximum absolute relative mismatch of exit and throat mass rates.
    pub relative_mass_flow_tolerance: f64,
}
impl RocketNozzleSolve {
    fn validate(self) -> Result<(), GasEnergyError> {
        self.nozzle.validate()?;
        if !positive(self.relative_mass_flow_tolerance) || self.relative_mass_flow_tolerance >= 1.0
        {
            return Err(invalid(
                "rocket nozzle relative mass-flow tolerance must lie in (0,1)",
            ));
        }
        Ok(())
    }
}

/// One accepted source-backed rocket flow at the unchanged reservoir state.
#[derive(Clone, Debug, PartialEq)]
pub struct RocketNozzleFlow {
    throat: ResolvedThermochemicalGasState,
    exit: ResolvedThermochemicalGasState,
    /// Choked propellant mass rate fixed by chamber state and throat [kg/s].
    pub mass_flow_kg_s: f64,
    /// Exit speed relative to the nozzle [m/s].
    pub exit_velocity_m_s: f64,
    /// Exit velocity divided by the local frozen sound speed.
    pub exit_mach_number: f64,
    /// Axial exhaust momentum contribution, mdot*v_exit [N].
    pub momentum_thrust_n: f64,
    /// Axial pressure contribution, (p_exit-p_ambient)*A_exit [N].
    pub pressure_thrust_n: f64,
    /// Total axial thrust on the vehicle, opposite exhaust direction [N].
    pub thrust_n: f64,
    /// Total thrust divided by actual propellant mass rate [m/s].
    pub effective_exhaust_velocity_m_s: f64,
    /// Characteristic velocity p_chamber*A_throat/mdot [m/s].
    pub characteristic_velocity_m_s: f64,
    /// Effective exhaust velocity divided by standard gravity [s].
    pub specific_impulse_s: f64,
    /// Accepted exit mass rate divided by throat mass rate, minus one.
    pub relative_mass_flow_residual: f64,
    /// Actual h_exit-h_chamber + v_exit²/2 [J/kg], reference independent.
    pub energy_residual_j_kg: f64,
    /// Chamber temperature times the actual exit entropy change [J/kg].
    pub entropy_energy_residual_j_kg: f64,
    /// Total interior iterations over the throat and exit roots.
    pub iterations: usize,
}
impl RocketNozzleFlow {
    /// Source-resolved sonic throat, including the original material receipts.
    #[must_use]
    pub const fn throat(&self) -> &ResolvedThermochemicalGasState {
        &self.throat
    }
    /// Source-resolved supersonic exit (sonic when the two areas coincide).
    #[must_use]
    pub const fn exit(&self) -> &ResolvedThermochemicalGasState {
        &self.exit
    }
}

impl ClosedGasParcel {
    /// Expand through an ideal choked rocket nozzle into explicit ambient pressure.
    ///
    /// Both root solves poll cancellation and use allocation-free caloric
    /// evaluations. Only accepted endpoints resolve full transport/source data.
    /// A source gap, exhausted budget or unsupported back pressure returns no
    /// flow candidate. No inventory, momentum or thermal state is mutated.
    pub fn rocket_nozzle_flow(
        &self,
        ambient_pressure_pa: f64,
        nozzle: RocketNozzle,
        controls: RocketNozzleSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RocketNozzleFlow, GasEnergyError> {
        nozzle.validate()?;
        controls.validate()?;
        if !ambient_pressure_pa.is_finite() || ambient_pressure_pa < 0.0 {
            return Err(invalid(
                "rocket ambient pressure must be finite and nonnegative",
            ));
        }
        poll(&mut cancelled)?;
        // Vacuum requests the sonic branch without evaluating an EOS at p=0.
        let sonic = self.nozzle_flow_to(
            0.0,
            GasNozzle {
                area_m2: nozzle.throat_area_m2,
                discharge_coefficient: 1.0,
            },
            controls.nozzle,
            &mut cancelled,
        )?;
        let reservoir = self.state().state();
        let (t0, p0, r) = (
            reservoir.temperature,
            reservoir.pressure,
            reservoir.specific_gas_constant,
        );
        let thermodynamics = self.parcel.model.thermodynamics();
        let target_flux = sonic.mass_flow_kg_s / nozzle.exit_area_m2;
        if !positive(target_flux) {
            return Err(invalid("rocket exit mass flux is unrepresentable"));
        }

        // Positive residual at the sonic endpoint, negative below the requested
        // supersonic exit. This bracket explicitly excludes the subsonic root.
        let sample = |t: f64| -> Result<(f64, f64, f64), GasEnergyError> {
            let dh = thermodynamics
                .specific_enthalpy_change(Temperature::new(t0), Temperature::new(t))
                .map_err(source_error)?
                .value();
            let ds = thermodynamics
                .specific_entropy_change(
                    Temperature::new(t0),
                    Pressure::new(p0),
                    Temperature::new(t),
                    Pressure::new(p0),
                )
                .map_err(source_error)?
                .value();
            let pressure = p0 * det::exp(ds / r);
            let velocity_squared = -2.0 * dh;
            if !positive(pressure) || !positive(velocity_squared) {
                return Err(invalid(
                    "rocket expansion pressure or kinetic energy is unrepresentable",
                ));
            }
            let velocity = det::sqrt(velocity_squared);
            let residual = (pressure / (r * t)) * velocity / target_flux - 1.0;
            if !residual.is_finite() {
                return Err(invalid("rocket mass-flux residual is unrepresentable"));
            }
            Ok((residual, pressure, velocity))
        };
        let mut temperature = sonic.throat().state().temperature;
        let mut iterations = 0;
        let (pressure, velocity) = if nozzle.exit_area_m2 == nozzle.throat_area_m2 {
            (sonic.throat().state().pressure, sonic.velocity_m_s)
        } else {
            let (mut lo, mut hi) = (controls.nozzle.minimum_temperature_k, temperature);
            poll(&mut cancelled)?;
            let mut value = sample(lo)?;
            if value.0 > controls.relative_mass_flow_tolerance {
                return Err(invalid(
                    "rocket supersonic exit lies below the admitted temperature bracket",
                ));
            }
            temperature = lo;
            while value.0.abs() > controls.relative_mass_flow_tolerance {
                if iterations == controls.nozzle.maximum_iterations {
                    return Err(invalid(
                        "rocket exit mass-flow root exhausted its iteration budget",
                    ));
                }
                let next = lo + 0.5 * (hi - lo);
                if next <= lo || next >= hi {
                    return Err(invalid(
                        "rocket exit mass-flow root exhausted floating-point progress",
                    ));
                }
                poll(&mut cancelled)?;
                temperature = next;
                value = sample(temperature)?;
                iterations += 1;
                if value.0 < 0.0 {
                    lo = temperature;
                } else {
                    hi = temperature;
                }
            }
            (value.1, value.2)
        };
        // Do not use an isentropic expansion as a silent shock/separation model.
        if ambient_pressure_pa > pressure {
            return Err(invalid(
                "overexpanded rocket nozzle requires a shock/separation model",
            ));
        }
        poll(&mut cancelled)?;
        let exit = if nozzle.exit_area_m2 == nozzle.throat_area_m2 {
            sonic.throat().clone()
        } else {
            self.parcel
                .model
                .evaluate(&updated_point(self.point(), temperature, pressure)?)?
        };
        let dh = thermodynamics
            .specific_enthalpy_change(Temperature::new(t0), Temperature::new(temperature))
            .map_err(source_error)?
            .value();
        let energy_residual = dh + 0.5 * velocity * velocity;
        let entropy_residual = t0
            * thermodynamics
                .specific_entropy_change(
                    Temperature::new(t0),
                    Pressure::new(p0),
                    Temperature::new(temperature),
                    Pressure::new(pressure),
                )
                .map_err(source_error)?
                .value();
        let mass_residual = exit.state().density * velocity / target_flux - 1.0;
        let mach = velocity / exit.state().sound_speed;
        if !mach.is_finite()
            || mach <= 0.0
            || (nozzle.exit_area_m2 > nozzle.throat_area_m2 && mach < 1.0)
            || !mass_residual.is_finite()
            || mass_residual.abs() > controls.relative_mass_flow_tolerance
            || !energy_residual.is_finite()
            || energy_residual.abs() > 2.0 * controls.nozzle.energy_tolerance_j_kg
            || !entropy_residual.is_finite()
            || entropy_residual.abs() > 2.0 * controls.nozzle.energy_tolerance_j_kg
        {
            return Err(invalid(
                "rocket exit failed mass, energy, entropy or supersonic admission",
            ));
        }
        let momentum = sonic.mass_flow_kg_s * velocity;
        let pressure_thrust = (pressure - ambient_pressure_pa) * nozzle.exit_area_m2;
        let thrust = momentum + pressure_thrust;
        let effective = thrust / sonic.mass_flow_kg_s;
        let characteristic = (p0 / target_flux) * (nozzle.throat_area_m2 / nozzle.exit_area_m2);
        let specific_impulse = effective / STANDARD_GRAVITY_M_S2;
        if !positive(momentum)
            || !pressure_thrust.is_finite()
            || !positive(thrust)
            || !positive(effective)
            || !positive(characteristic)
            || !positive(specific_impulse)
        {
            return Err(invalid("rocket thrust or performance is unrepresentable"));
        }
        poll(&mut cancelled)?;
        Ok(RocketNozzleFlow {
            throat: sonic.throat().clone(),
            exit,
            mass_flow_kg_s: sonic.mass_flow_kg_s,
            exit_velocity_m_s: velocity,
            exit_mach_number: mach,
            momentum_thrust_n: momentum,
            pressure_thrust_n: pressure_thrust,
            thrust_n: thrust,
            effective_exhaust_velocity_m_s: effective,
            characteristic_velocity_m_s: characteristic,
            specific_impulse_s: specific_impulse,
            relative_mass_flow_residual: mass_residual,
            energy_residual_j_kg: energy_residual,
            entropy_energy_residual_j_kg: entropy_residual,
            iterations: sonic.iterations + iterations,
        })
    }
}

#[cfg(test)]
mod tests;
