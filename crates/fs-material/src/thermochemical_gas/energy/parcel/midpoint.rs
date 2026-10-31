//! Simultaneous closed-gas heating and volume work in conservative coordinates.
//!
//! This is implicit midpoint in (specific internal energy, volume), NOT an
//! average of endpoint temperatures or an isentropic shortcut. At Vm=(V0+V1)/2,
//! solve Delta u(T0,Tm)+pm*(V1-V0)/(2m)=Q/(2m), then reconstruct the endpoint
//! from Delta u(T0,T1)=2*Delta u(T0,Tm), independently auditing Q-pm*(V1-V0).
//! Both increments use the original NASA-9
//! difference law, so formation references cannot erase heat or boundary work.
//! Mechanical consumers must use THIS pm, not endpoint pressure, in their load.
//!
//! Frozen composition and spatially uniform gas only. Midpoint is a time
//! discretization, not an exact finite isentrope: entropy is reported, not
//! constrained to be conserved or nondecreasing. The caller owns physical heat
//! transfer, motion, entropy accounting and timestep-refinement checks.
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
    /// Original-law m*Delta u between accepted thermodynamic states [J].
    pub stored_energy_change_j: f64,
    /// Actual stored energy plus work minus heat [J].
    pub energy_residual_j: f64,
    /// Actual m*(s1-s0) [J/K]; NOT entropy production or an isentropy certificate.
    pub entropy_change_j_k: f64,
    /// Original-law 2*(u_mid-u_initial)-(u_final-u_initial) [J/kg].
    /// Each increment includes any source-owned density dependence.
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
pub struct GasHeatVolumeCandidate<P = ClosedGasParcel, S = ResolvedThermochemicalGasState> {
    endpoint: P,
    midpoint: S,
    step: GasHeatVolumeStep,
}
impl<P, S> GasHeatVolumeCandidate<P, S> {
    /// Shared construction boundary for original ideal and real material owners.
    pub(in crate::thermochemical_gas) fn from_parts(endpoint: P, midpoint: S, step: GasHeatVolumeStep) -> Self {
        Self { endpoint, midpoint, step }
    }
    /// Complete endpoint medium/history for the next accepted physical step.
    #[must_use]
    pub const fn endpoint(&self) -> &P { &self.endpoint }
    /// Common midpoint medium. Its pressure is the effort used in boundary work.
    #[must_use]
    pub const fn midpoint(&self) -> &S { &self.midpoint }
    /// Actual balance and analytic discrete derivatives.
    #[must_use]
    pub const fn step(&self) -> &GasHeatVolumeStep { &self.step }
    /// Move the accepted endpoint and accounting without a second material query.
    #[must_use]
    pub fn into_parts(self) -> (P, GasHeatVolumeStep) { (self.endpoint, self.step) }
}

fn invalid(message: &'static str) -> GasEnergyError { GasEnergyError::Invalid(message) }

impl ClosedGasParcel {
    /// Tentatively advance heat and volume TOGETHER using the shared material law.
    /// `heat_j` is whole-step energy, not watts; positive heat enters the gas.
    /// The caller owns the time interval and its work-conjugate motion solve.
    ///
    /// Prepare source continuity and positive-capacity bounds once for both
    /// stages, then evaluate original caloric increments without trial receipts.
    /// Each bounded inversion uses at most maximum_iterations+2 samples; each
    /// accepted stage resolves original transport cards. Both stages and the old
    /// temperature must lie in the continuous admitted branch. Joule and stage-
    /// energy gates precede return, without bracket-width or clipped fallbacks.
    /// For fixed volume this agrees with nonlinear add_heat. At zero heat AND
    /// zero motion the endpoint is preserved bit-for-bit.
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
        let unchanged = dv == 0.0 && heat_j == 0.0;
        let branch = if unchanged { None }
            else { Some(old.model.prepare_caloric_branch(inner, &mut cancelled)?) };
        let increment = |t| -> Result<f64, GasEnergyError> {
            old.model.thermodynamics().specific_internal_energy_change(
                Temperature::new(temperature), Temperature::new(t),
            ).map(|value| value.value())
                .map_err(|error| super::super::super::failure("midpoint energy increment", error).into())
        };
        let (midpoint, midpoint_iterations) = if unchanged {
            (old.state.clone(), 0)
        } else {
            let thermal = |t: f64| -> Result<root::Sample, GasEnergyError> {
                let reference = old.model.thermodynamics().components()[0].convention().reference_pressure();
                let properties = old.model.thermodynamics().evaluate_properties(Temperature::new(t), reference)
                    .map_err(|error| super::super::super::failure("midpoint energy", error))?;
                let r = properties.mass_specific_gas_constant().value();
                Ok(root::Sample {
                    residual: increment(t)? + 0.5 * (ratio * r * t - q),
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
        // Enforce the actual caloric midpoint relation at the endpoint. The
        // independent first-law gate below still checks the requested heat and
        // actual pressure work, not a recomputed or rounded heat input.
        let change = 2.0 * increment(gas.temperature)?;
        if !work.is_finite() || !change.is_finite() || (dv != 0.0 && work == 0.0) {
            return Err(invalid("midpoint pressure work or specific energy change is unrepresentable"));
        }
        let (endpoint, endpoint_iterations) = if let Some(branch) = branch {
            let solved = branch.reconstruct_density_internal_energy_increment(
                &old.point, old.mass_kg / volume_m3, change, &mut cancelled)?;
            let mut parcel = old.clone();
            parcel.point = updated_point(&old.point, solved.state().state().temperature, solved.state().state().pressure)?;
            parcel.state = solved.state().clone();
            parcel.volume_m3 = volume_m3;
            (ClosedGasParcel::from_isochoric(parcel), solved.iterations())
        } else { (self.clone(), 0) };
        let du = increment(endpoint.state().state().temperature)?;
        let stored = old.mass_kg * du;
        let residual = (stored + work) - heat_j;
        let stage_residual = 2.0 * increment(midpoint.state().temperature)? - du;
        let entropy_change = old.mass_kg * old.model.specific_entropy_change_j_kg_k(
            &old.point, endpoint.point())?;
        if !stored.is_finite() || !residual.is_finite() || residual.abs() > tolerance_j
            || !stage_residual.is_finite() || stage_residual.abs() > controls.energy_tolerance_j_kg
            || !entropy_change.is_finite() {
            return Err(invalid("midpoint/endpoint material states did not close heat-work energy balance"));
        }
        // Differentiate the same discrete residual, not a static isentropic law.
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
#[cfg(test)]
mod increments;
