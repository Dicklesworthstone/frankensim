//! Discrete force and heat sensitivities of the already accepted midpoint gas step.
//!
//! Generalized pressure loading is the pullback of the actual midpoint pressure:
//! F_i = (p_mid - p_ext) dV/dq_i. Its Jacobian includes BOTH the gas law's
//! discrete d(p_mid)/d(V_final) and the geometry Hessian. In particular this
//! is not the static isentropic stiffness, which is only valid for an
//! infinitesimal reversible equilibrium perturbation, not an implicit heated
//! time step. The same callback also supplies an exact heat-to-force tangent.
//! No extra EOS, fluid model, dense matrix, or geometry derivative is invented.

use super::GasHeatVolumeCandidate;
use super::super::GasEnergyError;

/// The maximum caller-owned generalized coordinates per discrete gas load.
pub const MAX_DISCRETE_GAS_COORDINATES: usize = 256;

/// Borrowed fixed-material, fixed-geometry trial. The midpoint and endpoint
/// cannot change while the mechanical solver uses this tangent.
#[derive(Debug)]
pub struct DiscreteGasLoad<'a> {
    candidate: &'a GasHeatVolumeCandidate,
    gradient: &'a [f64],
    net_pressure_pa: f64,
}

fn invalid(msg: &'static str) -> GasEnergyError { GasEnergyError::Invalid(msg) }

fn finite(value: f64) -> Result<f64, GasEnergyError> {
    if value.is_finite() { Ok(value) } else { Err(invalid("nonfinite discrete gas mechanical response")) }
}

fn vector<'a>(values: &'a [f64], len: usize) -> Result<&'a [f64], GasEnergyError> {
    if values.len() != len || values.iter().any(|v| !v.is_finite()) {
        return Err(invalid("discrete gas coordinate vector has invalid length or entries"));
    }
    Ok(values)
}

fn dot(a: &[f64], b: &[f64]) -> Result<f64, GasEnergyError> {
    // Deterministic fixed-order compensated sum independent of worker count.
    let mut sum = 0.0;
    let mut correction = 0.0;
    for (&x, &y) in a.iter().zip(b) {
        let term = finite(x * y)?;
        let y = term - correction;
        let next = finite(sum + y)?;
        correction = (next - sum) - y;
        sum = next;
    }
    Ok(sum)
}

impl GasHeatVolumeCandidate {
    /// Bind the supplied *final* geometry to the material midpoint pressure.
    ///
    /// The gradient is dV_final/dq. It must use the same ordered generalized
    /// coordinates as the displacement, velocities and Hessian action. No
    /// default atmosphere is inferred: the exterior pressure is explicit and
    /// constant during this linearization (zero means a vacuum).
    ///
    /// This is an allocation-free binding, not a state update. The caller
    /// owns geometry and its exact derivatives, as well as any non-gas forces.
    pub fn discrete_mechanical_load<'a>(
        &'a self,
        final_volume_m3: f64,
        volume_gradient: &'a [f64],
        exterior_pressure_pa: f64,
    ) -> Result<DiscreteGasLoad<'a>, GasEnergyError> {
        if final_volume_m3 != self.step().final_volume_m3
            || !exterior_pressure_pa.is_finite() || exterior_pressure_pa < 0.0
            || volume_gradient.is_empty()
            || volume_gradient.len() > MAX_DISCRETE_GAS_COORDINATES
        {
            return Err(invalid("discrete load requires exact final volume, nonnegative exterior pressure and bounded coordinate count"));
        }
        vector(volume_gradient, volume_gradient.len())?;
        let net_pressure_pa = finite(self.midpoint().state().pressure - exterior_pressure_pa)?;
        Ok(DiscreteGasLoad { candidate: self, gradient: volume_gradient, net_pressure_pa })
    }
}

impl DiscreteGasLoad<'_> {
    /// Exactly the material candidate used to evaluate midpoint gas work.
    #[must_use]
    pub const fn candidate(&self) -> &GasHeatVolumeCandidate { self.candidate }

    /// Signed pressure difference, gas minus explicitly declared exterior [Pa].
    #[must_use]
    pub const fn net_pressure_pa(&self) -> f64 { self.net_pressure_pa }

    /// Force in the caller's generalized coordinate chart: pressure * dV/dq.
    pub fn generalized_forces(&self) -> Result<Vec<f64>, GasEnergyError> {
        self.gradient.iter().map(|g| finite(self.net_pressure_pa * g)).collect()
    }

    /// Matrix-free derivative of the complete gas force with respect to the
    /// final mechanical coordinates, keeping supplied step heat constant.
    ///
    /// The Hessian action is D²V(q)[direction], not a frozen volume gradient.
    /// This accounts for geometry curvature and the original discrete
    /// nonlinear thermodynamic law in O(number_of_coordinates) work.
    pub fn tangent_action(
        &self, direction: &[f64], volume_hessian_action: &[f64],
    ) -> Result<Vec<f64>, GasEnergyError> {
        vector(direction, self.gradient.len())?;
        vector(volume_hessian_action, self.gradient.len())?;
        let pressure_change = finite(
            self.candidate.step().pressure_final_volume_derivative_pa_m3
                * dot(self.gradient, direction)?
        )?;
        self.gradient.iter().zip(volume_hessian_action).map(|(&g, &h)| {
            finite(finite(pressure_change * g)? + finite(self.net_pressure_pa * h)?)
        }).collect()
    }

    /// Derivative of generalized force with respect to the total signed heat
    /// input [N/J for a translation chart], at fixed final geometry.
    pub fn heat_tangent(&self) -> Result<Vec<f64>, GasEnergyError> {
        self.gradient.iter().map(|g| {
            finite(self.candidate.step().pressure_heat_derivative_pa_j * g)
        }).collect()
    }

    /// Mechanical power transferred from gas minus exterior pressure [W].
    /// The whole-step gas work remains the candidate's midpoint p*dV; a
    /// caller integrates this instantaneous rate using its own time model.
    pub fn mechanical_power_w(&self, coordinate_rates: &[f64]) -> Result<f64, GasEnergyError> {
        vector(coordinate_rates, self.gradient.len())?;
        finite(self.net_pressure_pa * dot(self.gradient, coordinate_rates)?)
    }
}

#[cfg(test)]
mod tests;
