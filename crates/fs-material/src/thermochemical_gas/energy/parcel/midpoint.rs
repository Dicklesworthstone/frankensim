//! Simultaneous closed-gas heating and volume work in conservative coordinates.
//!
//! This is the implicit midpoint material update in (specific internal energy,
//! volume), NOT an average of endpoint temperatures or an isentropic shortcut.
//! At Vm=(V0+V1)/2, solve 2*(u(Tm)-u0)+pm*(V1-V0)/m=Q/m. The endpoint then uses
//! the existing density/energy closure with m*(u1-u0)=Q-pm*(V1-V0). Both stages
//! resolve the original NASA-9, EOS and transport sources. Mechanical consumers
//! must use THIS pm for their work-conjugate load, not endpoint pressure.
//!
//! Frozen composition and spatially uniform gas only. Midpoint is a time
//! discretization, not an exact finite isentrope: entropy is reported, not
//! constrained to be conserved or nondecreasing. A caller still owns physical
//! heat transfer, motion, entropy accounting and timestep-refinement checks.
//! First-law sign convention: https://www.grc.nasa.gov/www/k-12/airplane/thermo1.html

/// Discrete mechanical force and sensitivity pullbacks for the same midpoint state.
pub mod mechanical_load;

use super::ClosedGasParcel;
use super::super::{GasEnergyError, GasEnergySolve, ResolvedThermochemicalGasState,
    map_root, poll, positive, root, updated_point};
use fs_qty::Temperature;

/// Accepted heat/work accounting and derivatives of the DISCRETE pressure law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasHeatVolumeStep {
    /// Initial volume [m3].
    pub initial_volume_m3: f64,
    /// Final volume [m3].
    pub final_volume_m3: f64,
    /// Prescribed heat INTO the gas over this whole step [J].
    pub supplied_heat_j: f64,
    /// Work OUT of the gas: midpoint pressure times volume change [J].
    pub work_by_gas_j: f64,
    /// Actual m*(u1-u0) from the accepted material states [J].
    pub stored_energy_change_j: f64,
    /// Actual stored energy plus work minus heat [J].
    pub energy_residual_j: f64,
    /// Actual m*(s1-s0) [J/K]; NOT entropy production or an isentropy certificate.
    pub entropy_change_j_k: f64,
    /// Actual 2*(um-u0)-(u1-u0) [J/kg].
    pub midpoint_energy_residual_j_kg: f64,
    /// d(pm)/d(V1) with initial state and total heat held fixed [Pa/m3].
    pub pressure_final_volume_derivative_pa_m3: f64,
    /// d(pm)/d(Q) with initial state and endpoint volume held fixed [Pa/J].
    pub pressure_heat_derivative_pa_j: f64,
    /// Interior iterations summed over midpoint and endpoint inversions.
    pub iterations: usize,
}

/// A tentative, complete gas update. Input storage is never mutated by a trial.
/// Publish the endpoint only after the mechanical/thermal coupling also accepts.
#[derive(Clone, Debug, PartialEq)]
pub struct GasHeatVolumeCandidate {
    endpoint: ClosedGasParcel,
    midpoint: ResolvedThermochemicalGasState,
    step: GasHeatVolumeStep,
}
impl GasHeatVolumeCandidate {
    /// Complete endpoint medium/history for the next accepted physical step.
    #[must_use]
    pub const fn endpoint(&self) -> &ClosedGasParcel { &self.endpoint }
    /// Common midpoint medium. Its pressure is the effort used in boundary work.
    #[must_use]
    pub const fn midpoint(&self) -> &ResolvedThermochemicalGasState { &self.midpoint }
    /// Actual balance and analytic discrete derivatives.
    #[must_use]
    pub const fn step(&self) -> &GasHeatVolumeStep { &self.step }
    /// Move the accepted endpoint and accounting without a second material query.
    #[must_use]
    pub fn into_parts(self) -> (ClosedGasParcel, GasHeatVolumeStep) { (self.endpoint, self.step) }
}

fn invalid(message: &'static str) -> GasEnergyError { GasEnergyError::Invalid(message) }

impl ClosedGasParcel {
    /// Tentatively advance heat and volume TOGETHER using the shared material law.
    /// `heat_j` is whole-step energy, not watts; positive heat enters the gas.
    /// The caller owns the time interval and its work-conjugate motion solve.
    ///
    /// Each bounded inversion uses at most maximum_iterations+2 thermal samples;
    /// each accepted stage then resolves the original transport cards. Both
    /// stages, the old temperature and all source conditions must be admitted.
    /// Absolute joule and stage-energy gates precede return. No bracket-width,
    /// nominal-gamma, clipped-temperature or partial-state fallback is present.
    ///
    /// For fixed volume this agrees with the existing nonlinear add_heat result,
    /// while retaining the midpoint needed by other physical consumers. At zero
    /// heat AND zero motion the endpoint is preserved bit-for-bit.
    pub fn midpoint_heat_volume(
        &self, volume_m3: f64, heat_j: f64, controls: GasEnergySolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasHeatVolumeCandidate, GasEnergyError> {
        controls.validate()?;
        poll(&mut cancelled)?;
        let old = &self.parcel;
        let temperature = old.state.state().temperature;
        if !positive(volume_m3) || !heat_j.is_finite()
            || temperature < controls.minimum_temperature_k || temperature > controls.maximum_temperature_k {
            return Err(invalid("positive finite volume, finite heat and bracket containing old temperature required"));
        }
        let vm = 0.5 * old.volume_m3 + 0.5 * volume_m3;
        let dv = volume_m3 - old.volume_m3;
        let density = old.mass_kg / vm;
        let q = heat_j / old.mass_kg;
        let tolerance_j = old.mass_kg * controls.energy_tolerance_j_kg;
        let inner = GasEnergySolve { energy_tolerance_j_kg: 0.25 * controls.energy_tolerance_j_kg, ..controls };
        inner.validate()?;
        if !positive(vm) || !positive(density) || !q.is_finite()
            || (heat_j != 0.0 && q == 0.0) || !positive(tolerance_j) {
            return Err(invalid("midpoint density, heat or absolute joule tolerance is unrepresentable"));
        }
        let ratio = dv / vm;
        let initial_u = old.state.specific_internal_energy_j_kg();
        let unchanged = dv == 0.0 && heat_j == 0.0;
        let (midpoint, midpoint_iterations) = if unchanged {
            (old.state.clone(), 0)
        } else {
            let thermal = |t: f64| -> Result<root::Sample, GasEnergyError> {
                let reference = old.model.thermodynamics().components()[0].convention().reference_pressure();
                let evaluated = old.model.thermodynamics().evaluate(Temperature::new(t), reference)
                    .map_err(|error| super::super::super::failure("midpoint energy", error))?;
                let properties = evaluated.properties();
                let r = properties.mass_specific_gas_constant().value();
                let u = properties.u().value() / properties.molar_mass().value();
                Ok(root::Sample {
                    residual: (u - initial_u) + 0.5 * (ratio * r * t - q),
                    tangent: properties.mass_specific_cv().value() + 0.5 * ratio * r,
                    pressure_pa: density * r * t,
                })
            };
            let solved = root::solve(inner.minimum_temperature_k, inner.maximum_temperature_k,
                inner.energy_tolerance_j_kg, inner.maximum_iterations, thermal, &mut cancelled).map_err(map_root)?;
            poll(&mut cancelled)?;
            let state = old.model.evaluate(&updated_point(&old.point, solved.temperature_k, solved.sample.pressure_pa)?)?;
            if (state.state().density / density - 1.0).abs() > 64.0 * f64::EPSILON {
                return Err(invalid("midpoint EOS did not retain the conserved mass/volume density"));
            }
            (state, solved.iterations)
        };
        let gas = midpoint.state();
        let work = gas.pressure * dv;
        let change = (heat_j - work) / old.mass_kg;
        if !work.is_finite() || !change.is_finite() || (dv != 0.0 && work == 0.0) {
            return Err(invalid("midpoint pressure work or specific energy change is unrepresentable"));
        }
        let (endpoint, endpoint_iterations) = if unchanged {
            (self.clone(), 0)
        } else {
            let solved = old.model.solve_energy_change(&old.point, old.mass_kg / volume_m3,
                initial_u, change, inner, &mut cancelled)?;
            let mut parcel = old.clone();
            parcel.point = updated_point(&old.point, solved.state().state().temperature, solved.state().state().pressure)?;
            parcel.state = solved.state().clone();
            parcel.volume_m3 = volume_m3;
            (ClosedGasParcel::from_isochoric(parcel), solved.iterations())
        };
        let du = endpoint.state().specific_internal_energy_j_kg() - initial_u;
        let stored = old.mass_kg * du;
        let residual = (stored + work) - heat_j;
        let stage_residual = 2.0 * (midpoint.specific_internal_energy_j_kg() - initial_u) - du;
        let entropy_change = old.mass_kg * old.model.specific_entropy_change_j_kg_k(
            &old.point, endpoint.point())?;
        if !stored.is_finite() || !residual.is_finite() || residual.abs() > tolerance_j
            || !stage_residual.is_finite() || stage_residual.abs() > controls.energy_tolerance_j_kg
            || !entropy_change.is_finite() {
            return Err(invalid("midpoint/endpoint material states did not close heat-work energy balance"));
        }
        // Differentiate the same residual, using the owner's actual cv(Tm).
        // These are DISCRETE derivatives; the isentropic static tangent is not
        // interchangeable with them during a heated midpoint coupling solve.
        let cv = midpoint.thermodynamics().properties().mass_specific_cv().value();
        let r = gas.specific_gas_constant;
        let denominator = cv + 0.5 * ratio * r;
        let dpdv = -(gas.pressure / vm) * (0.5 + 0.5 * r * (old.volume_m3 / vm) / denominator);
        let dpdq = 0.5 * (gas.pressure / gas.temperature / denominator) / old.mass_kg;
        if !positive(denominator) || !dpdv.is_finite() || dpdv >= 0.0 || !positive(dpdq) {
            return Err(invalid("discrete material pressure derivatives are unrepresentable"));
        }
        let step = GasHeatVolumeStep {
            initial_volume_m3: old.volume_m3, final_volume_m3: volume_m3, supplied_heat_j: heat_j,
            work_by_gas_j: work, stored_energy_change_j: stored, energy_residual_j: residual,
            entropy_change_j_k: entropy_change, midpoint_energy_residual_j_kg: stage_residual,
            pressure_final_volume_derivative_pa_m3: dpdv, pressure_heat_derivative_pa_j: dpdq,
            iterations: midpoint_iterations + endpoint_iterations,
        };
        poll(&mut cancelled)?;
        Ok(GasHeatVolumeCandidate { endpoint, midpoint, step })
    }
}

#[cfg(test)]
mod tests;
