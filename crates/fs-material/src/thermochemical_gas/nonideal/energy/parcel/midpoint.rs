//! Conservative midpoint heat/work updates of the original real-gas parcel.
//!
//! At Vm=(V0+V1)/2, solve u(Tm,m/Vm)-u(T0,rho0)+pm*(V1-V0)/(2m)=Q/(2m).
//! The endpoint has u(T1,m/V1)-u(T0,rho0)=2*(um-u0). Both differences retain
//! the original initial density in the PR departure, not a fictitious initial
//! state already compressed to the final density. NASA reference constants
//! cancel before evaluating either increment. This is midpoint in energy and
//! volume, not in temperature, and not an exact finite-step isentrope.

use super::super::CaloricSample;
use super::{IsochoricRealGasParcel, ResolvedThermochemicalRealGasState};
use crate::thermochemical_gas::energy::parcel::midpoint::{
    GasHeatVolumeCandidate, GasHeatVolumeStep,
};
use crate::thermochemical_gas::energy::{
    GasEnergyError, GasEnergySolve, map_root, poll, root, updated_point,
};
use crate::thermochemical_gas::{failure, positive};
use fs_qty::{Pressure, Temperature};

/// The original common candidate/report owner, retaining complete real states.
pub type RealGasHeatVolumeCandidate =
    GasHeatVolumeCandidate<IsochoricRealGasParcel, ResolvedThermochemicalRealGasState>;

fn invalid(message: &'static str) -> GasEnergyError {
    GasEnergyError::Invalid(message)
}

impl IsochoricRealGasParcel {
    /// Construct an immutable simultaneous heat/volume midpoint candidate.
    /// Positive `heat_j` enters the gas; mechanical work uses the returned
    /// MIDPOINT pressure times the actual volume change. The source-dependent
    /// endpoint pressure must not replace that effort in the coupled integrator.
    ///
    /// The same PR/NASA source supplies both accepted material states, complete
    /// Cv, pressure derivatives and `du/dnu = T*(dp/dT)_nu-p`, where nu=V/m.
    /// The shared bounded solver inverts both stages. Both stages must retain
    /// exact mass and pass complete material/transport admission before return.
    /// Original composition, sources and auxiliary query coordinates survive.
    ///
    /// A continuous homogeneous frozen-composition branch across the caller's
    /// temperature bracket and physical path is required. Local mechanical/Cv
    /// stability does not establish phase equilibrium or prevent coexistence;
    /// no latent heat or phase switching is inferred. Reported Delta S is not
    /// constrained to vanish or be nonnegative for a finite midpoint step. The
    /// caller owns heat-transfer physics, mechanical dynamics and refinement.
    /// Discrete dp/dQ retains its source-derived sign, including negative
    /// thermal expansion; it is not forced to an ideal-gas positive value.
    /// Cancellation and any failed trial leave this entire parcel unchanged.
    pub fn midpoint_heat_volume(
        &self,
        volume_m3: f64,
        heat_j: f64,
        controls: GasEnergySolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RealGasHeatVolumeCandidate, GasEnergyError> {
        controls.validate()?;
        poll(&mut cancelled)?;
        let temperature = self.state.temperature_k();
        if !positive(volume_m3)
            || !heat_j.is_finite()
            || temperature < controls.minimum_temperature_k
            || temperature > controls.maximum_temperature_k
        {
            return Err(invalid(
                "real-gas midpoint requires positive volume, finite heat and bracketed initial temperature",
            ));
        }
        let vm = 0.5 * self.volume_m3 + 0.5 * volume_m3;
        let dv = volume_m3 - self.volume_m3;
        let midpoint_density = if dv == 0.0 {
            self.density_kg_m3
        } else {
            self.mass_kg / vm
        };
        let final_density = if dv == 0.0 {
            self.density_kg_m3
        } else {
            self.mass_kg / volume_m3
        };
        let q = heat_j / self.mass_kg;
        let tolerance_j = self.mass_kg * controls.energy_tolerance_j_kg;
        let mut inner = GasEnergySolve {
            energy_tolerance_j_kg: 0.25 * controls.energy_tolerance_j_kg,
            ..controls
        };
        if dv == 0.0 && heat_j != 0.0 {
            inner.energy_tolerance_j_kg = inner.energy_tolerance_j_kg.min(0.25 * q.abs());
        }
        inner.validate()?;
        if !positive(vm)
            || !positive(midpoint_density)
            || !positive(final_density)
            || !q.is_finite()
            || (heat_j != 0.0 && q == 0.0)
            || !positive(tolerance_j)
            || (dv != 0.0 && final_density == self.density_kg_m3)
        {
            return Err(invalid(
                "real-gas midpoint density, heat, volume change or tolerance is unrepresentable",
            ));
        }
        let half_specific_volume_change = 0.5 * (dv / self.mass_kg);
        if !half_specific_volume_change.is_finite() {
            return Err(invalid(
                "real-gas midpoint specific volume change is unrepresentable",
            ));
        }
        let molar_mass = self.model.mechanics().molar_mass().value();
        let initial_residual = self.state.root().residual_internal_energy_j_mol() / molar_mass;
        let increment = |t, departure| {
            self.model
                .caloric_increment(temperature, t, initial_residual, departure)
        };
        let unchanged = dv == 0.0 && heat_j == 0.0;
        let (midpoint, midpoint_iterations) = if unchanged {
            (self.state.clone(), 0)
        } else {
            let solved = root::solve(
                inner.minimum_temperature_k,
                inner.maximum_temperature_k,
                inner.energy_tolerance_j_kg,
                inner.maximum_iterations,
                |t| {
                    let sample: CaloricSample =
                        self.model.caloric_at_density(midpoint_density, t)?;
                    Ok::<_, GasEnergyError>(root::Sample {
                        residual: increment(t, sample.residual_energy_j_kg)?
                            + sample.pressure_pa * half_specific_volume_change
                            - 0.5 * q,
                        tangent: sample.cv_j_kg_k
                            + sample.pressure_temperature_derivative * half_specific_volume_change,
                        pressure_pa: sample.pressure_pa,
                    })
                },
                &mut cancelled,
            )
            .map_err(map_root)?;
            poll(&mut cancelled)?;
            let point =
                updated_point(&self.point, solved.temperature_k, solved.sample.pressure_pa)?;
            let state = self
                .model
                .evaluate_density_temperature(&point, midpoint_density)?;
            if state.density_kg_m3().to_bits() != midpoint_density.to_bits() {
                return Err(invalid(
                    "real-gas midpoint reconstruction changed conserved density",
                ));
            }
            (state, solved.iterations)
        };
        let midpoint_increment = increment(
            midpoint.temperature_k(),
            midpoint.root().residual_internal_energy_j_mol() / molar_mass,
        )?;
        let change = 2.0 * midpoint_increment;
        let work = midpoint.pressure_pa() * dv;
        if !change.is_finite() || !work.is_finite() || (dv != 0.0 && work == 0.0) {
            return Err(invalid(
                "real-gas midpoint pressure work or energy change is unrepresentable",
            ));
        }
        let (endpoint, endpoint_iterations) = if unchanged {
            (self.clone(), 0)
        } else {
            // The shared increment inversion retains the ORIGINAL density as
            // its energy reference, including density-dependent PR departures.
            let solved = self
                .model
                .reconstruct_density_internal_energy_increment_between_densities(
                    &self.point,
                    self.density_kg_m3,
                    final_density,
                    change,
                    inner,
                    &mut cancelled,
                )?;
            let iterations = solved.iterations();
            let state = solved.into_state();
            let point = updated_point(&self.point, state.temperature_k(), state.pressure_pa())?;
            let mut endpoint = self.clone();
            endpoint.point = point;
            endpoint.state = state;
            endpoint.density_kg_m3 = final_density;
            endpoint.volume_m3 = volume_m3;
            (endpoint, iterations)
        };
        let final_residual = endpoint.state.root().residual_internal_energy_j_mol() / molar_mass;
        let du = increment(endpoint.state.temperature_k(), final_residual)?;
        let stored = self.mass_kg * du;
        let energy_residual = (stored + work) - heat_j;
        let stage_residual = change - du;
        let ideal_entropy = self
            .model
            .source()
            .thermodynamics()
            .specific_entropy_change(
                Temperature::new(temperature),
                Pressure::new(self.state.pressure_pa()),
                Temperature::new(endpoint.state.temperature_k()),
                Pressure::new(endpoint.state.pressure_pa()),
            )
            .map_err(|error| failure("real-gas midpoint relative entropy", error))?
            .value();
        let entropy_change = self.mass_kg
            * (ideal_entropy
                + (endpoint.state.root().residual_entropy_j_mol_k()
                    - self.state.root().residual_entropy_j_mol_k())
                    / molar_mass);
        if !stored.is_finite()
            || !energy_residual.is_finite()
            || energy_residual.abs() > tolerance_j
            || !stage_residual.is_finite()
            || stage_residual.abs() > controls.energy_tolerance_j_kg
            || !entropy_change.is_finite()
        {
            return Err(invalid(
                "real-gas midpoint/endpoint did not close source heat-work energy balance",
            ));
        }
        // Differentiate the same conservative midpoint residual. The source
        // Maxwell identity u_nu=T*p_T-p cancels the explicit pressure terms.
        // This is the discrete heated law, not static isentropic stiffness.
        let cv = midpoint.specific_heat_cv_j_kg_k();
        let (p_t, p_v_molar) = midpoint.root().pressure_derivatives();
        let p_nu = molar_mass * p_v_molar;
        let denominator = cv + half_specific_volume_change * p_t;
        let dpdv =
            0.5 * ((p_nu * cv - midpoint.temperature_k() * p_t * p_t) / denominator) / self.mass_kg;
        let dpdq = 0.5 * (p_t / denominator) / self.mass_kg;
        if !positive(denominator)
            || !p_nu.is_finite()
            || p_nu >= 0.0
            || !dpdv.is_finite()
            || dpdv >= 0.0
            || !dpdq.is_finite()
        {
            return Err(invalid(
                "real-gas discrete midpoint pressure derivatives are unrepresentable",
            ));
        }
        let step = GasHeatVolumeStep {
            initial_volume_m3: self.volume_m3,
            final_volume_m3: volume_m3,
            supplied_heat_j: heat_j,
            work_by_gas_j: work,
            stored_energy_change_j: stored,
            energy_residual_j: energy_residual,
            entropy_change_j_k: entropy_change,
            midpoint_energy_residual_j_kg: stage_residual,
            pressure_final_volume_derivative_pa_m3: dpdv,
            pressure_heat_derivative_pa_j: dpdq,
            iterations: midpoint_iterations + endpoint_iterations,
        };
        poll(&mut cancelled)?;
        Ok(GasHeatVolumeCandidate::from_parts(endpoint, midpoint, step))
    }
}
