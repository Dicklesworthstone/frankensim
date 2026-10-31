//! Reversible volume work on the existing source-bound real-gas parcel.
use super::IsochoricRealGasParcel;
use crate::thermochemical_gas::energy::entropy::{GasEntropyError, GasEntropySolve, map_context};
use crate::thermochemical_gas::energy::parcel::GasVolumeStep;
use crate::thermochemical_gas::energy::{poll, updated_point};
use crate::thermochemical_gas::positive;

impl IsochoricRealGasParcel {
    /// Construct a reversible adiabatic volume endpoint at conserved exact mass.
    /// The original parcel is immutable, including on cancellation or refusal.
    /// The candidate retains the same composition, sources and extra query axes;
    /// subsequent `add_heat` holds its newly accepted volume and density fixed.
    ///
    /// The caller admits one continuous reversible homogeneous phase path and a
    /// temperature bracket containing the initial state. Existing PR mechanical
    /// and complete Cv checks admit each evaluated state. These local checks do
    /// not prove global phase stability or exclude coexistence along the path;
    /// no latent heat, phase conversion or equilibrium flash is supplied.
    ///
    /// The shared entropy reconstruction solves Delta s = 0 with derivative
    /// Cv/T at the final density. NASA reference constants cancel before
    /// evaluating the ideal entropy/energy changes; source-owned PR departures
    /// supply the real corrections. Work delivered by the gas is exactly minus
    /// the returned source-law Delta U, without subtracting large absolute
    /// stored energies. This is thermodynamic quasistatic work, not a pressure-
    /// quadrature or mechanical validation certificate.
    /// `work_by_gas_j / (V_final - V_initial)` is the work-conjugate secant
    /// pressure; `candidate.state().pressure_pa()` is the endpoint pressure.
    /// External loads, piston inertia and friction remain the composing
    /// mechanical owner's responsibility.
    pub fn at_volume_isentropic(
        &self,
        volume_m3: f64,
        controls: GasEntropySolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(Self, GasVolumeStep), GasEntropyError> {
        controls.validate()?;
        poll(&mut cancelled).map_err(map_context)?;
        let initial_temperature = self.state.temperature_k();
        if !positive(volume_m3)
            || initial_temperature < controls.minimum_temperature_k
            || initial_temperature > controls.maximum_temperature_k
        {
            return Err(GasEntropyError::Invalid(
                "positive real-gas volume and a bracket containing the initial temperature required",
            ));
        }
        let tolerance = self.mass_kg * controls.entropy_tolerance_j_kg_k;
        if !positive(tolerance) {
            return Err(GasEntropyError::Invalid(
                "total real-gas entropy tolerance is unrepresentable",
            ));
        }
        let mut report = GasVolumeStep {
            initial_volume_m3: self.volume_m3,
            final_volume_m3: volume_m3,
            stored_energy_change_j: 0.0,
            work_by_gas_j: 0.0,
            entropy_change_j_k: 0.0,
            iterations: 0,
        };
        if volume_m3 == self.volume_m3 {
            let candidate = self.clone();
            poll(&mut cancelled).map_err(map_context)?;
            return Ok((candidate, report));
        }
        let density = self.mass_kg / volume_m3;
        if !positive(density) || density == self.density_kg_m3 {
            return Err(GasEntropyError::Invalid(
                "real-gas volume change requires a representable positive density change",
            ));
        }
        let solved = self.model.reconstruct_density_entropy_change(
            &self.point,
            self.density_kg_m3,
            density,
            0.0,
            controls,
            &mut cancelled,
        )?;
        let entropy = self.mass_kg * solved.residual_j_kg_k();
        let stored = self.mass_kg
            * self
                .model
                .specific_internal_energy_change_between_densities(
                    self.density_kg_m3,
                    initial_temperature,
                    density,
                    solved.state().temperature_k(),
                )
                .map_err(map_context)?;
        let secant_pressure = -stored / (volume_m3 - self.volume_m3);
        if !stored.is_finite()
            || !entropy.is_finite()
            || entropy.abs() > tolerance
            || !positive(secant_pressure)
        {
            return Err(GasEntropyError::Invalid(
                "real-gas isentropic work or entropy balance is unrepresentable/inconsistent",
            ));
        }
        let point = updated_point(
            &self.point,
            solved.state().temperature_k(),
            solved.state().pressure_pa(),
        )
        .map_err(map_context)?;
        report.stored_energy_change_j = stored;
        report.work_by_gas_j = -stored;
        report.entropy_change_j_k = entropy;
        report.iterations = solved.iterations();
        let mut candidate = self.clone();
        candidate.point = point;
        candidate.state = solved.into_state();
        candidate.density_kg_m3 = density;
        candidate.volume_m3 = volume_m3;
        poll(&mut cancelled).map_err(map_context)?;
        Ok((candidate, report))
    }
}
