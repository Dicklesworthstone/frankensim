//! Advective-diffusive thermal contacts on the original characteristic faces.
//!
//! The acoustic characteristic still owns pressure and normal velocity. For
//! inflow, the original upstream specific enthalpy and molecular half-cell
//! Fourier flux determine one face temperature. This is a Danckwerts thermal
//! inlet condition, not a fitted entropy blend or an extra volume heat source.
//! Tangential donor selection is unchanged; arbitrary reversal is unqualified.

use super::{NsfError, NsfGrid3};
use crate::lattice3::Lattice3;

/// Thermal condition at each original characteristic face, low then high.
/// A selected face and its localized patch share the same thermal policy,
/// while keeping independent pressure/temperature drive data.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NsfCharacteristicThermalPolicy {
    /// Original convected entropy trace, including its original Fourier flux.
    UpwindEntropy,
    /// Upstream enthalpy enters with inward mass flow. The original material
    /// Fourier flux and face enthalpy balance; outward/stationary conduction
    /// is zero. No finite external inventory or species source is introduced.
    EnthalpyFlux,
}

impl<L: Lattice3<27>> NsfGrid3<L> {
    /// Select physical thermal contacts after attaching the characteristic axis
    /// and before the first step. The original gas, face quadrature, pressure
    /// and population/energy owners remain in force. Active SGS is refused
    /// until its boundary gradient and conductivity can be solved together.
    /// The exact-zero SGS control remains the molecular route.
    pub fn with_characteristic_thermal_policy(
        mut self,
        policies: [NsfCharacteristicThermalPolicy; 2],
    ) -> Result<Self, NsfError> {
        let open = self.open.ok_or(NsfError::Input(
            "characteristic thermal policy requires an attached characteristic axis",
        ))?;
        if self.steps != 0 {
            return Err(NsfError::Input(
                "characteristic thermal policy is cold-only",
            ));
        }
        if self.sgs.is_some() && policies.contains(&NsfCharacteristicThermalPolicy::EnthalpyFlux) {
            return Err(NsfError::Input(
                "enthalpy-flux characteristic thermal contacts do not yet support active SGS",
            ));
        }
        if policies == self.characteristic_thermal_policy {
            return Ok(self);
        }
        self.characteristic_thermal_policy = policies;
        self.rebind_characteristics(
            open,
            self.characteristic_pressure_offset_pa,
            self.characteristic_inflow_temperature_k,
            self.characteristic_patch,
        )?;
        Ok(self)
    }

    /// Selected low/high policies. An unattached grid retains the original
    /// defaults; selection itself requires a characteristic axis.
    #[must_use]
    pub const fn characteristic_thermal_policy(&self) -> [NsfCharacteristicThermalPolicy; 2] {
        self.characteristic_thermal_policy
    }

    // Use the same caloric chart as OpenFace::internal. In particular, the
    // legacy GasSpec route freezes gamma at construction; querying cp(T) here
    // would prescribe a different energy from the one the original FV stores.
    fn contact_enthalpy(&self, theta: f64, pressure_pa: f64) -> Result<(f64, f64), NsfError> {
        let temperature = self.law.kelvin_of_theta(theta);
        let rho = pressure_pa / (self.law.gas_constant * temperature * self.law.rho_ref_kg_m3);
        let material = self.spec.gas.resolve(&self.law, rho, theta)?;
        let internal = material
            .as_ref()
            .map_or(theta / (self.law.gamma - 1.0), |material| {
                material.sensible_energy_lattice
            });
        let enthalpy = internal + theta;
        if !(rho.is_finite() && rho > 0.0 && enthalpy.is_finite()) {
            return Err(NsfError::Input("nonfinite characteristic contact enthalpy"));
        }
        Ok((rho, enthalpy))
    }

    // Lattice form of |rho*u_n| (h_f-h_up) + (2*k/dx)(T_f-T_i) = 0.
    // The same face conductivity/mean-pressure query is used by
    // open_energy_flux. Its outward heat is -2*k*(theta_f-theta_i), so the
    // actual incoming total flux equals rho*u_n*(h_up + |u_f|^2/2), apart
    // from the separately retained viscous work and floating-point rounding.
    pub(super) fn enthalpy_flux_theta(
        &self,
        adjacent_theta: f64,
        adjacent_pressure_pa: f64,
        pressure_pa: f64,
        outward_velocity_lu: f64,
        upstream_temperature_k: f64,
    ) -> Result<f64, NsfError> {
        if !(pressure_pa.is_finite()
            && pressure_pa > 0.0
            && adjacent_pressure_pa.is_finite()
            && adjacent_pressure_pa > 0.0
            && outward_velocity_lu.is_finite())
        {
            return Err(NsfError::Input(
                "invalid characteristic thermal contact state",
            ));
        }
        if outward_velocity_lu >= 0.0 {
            return Ok(adjacent_theta);
        }
        if !upstream_temperature_k.is_finite()
            || upstream_temperature_k < self.spec.temperature_window_k[0]
            || upstream_temperature_k > self.spec.temperature_window_k[1]
        {
            return Err(NsfError::Input(
                "characteristic contact upstream temperature leaves its material window",
            ));
        }
        let upstream_theta = self.law.theta_of_kelvin(upstream_temperature_k);
        if upstream_theta == adjacent_theta {
            return Ok(adjacent_theta);
        }
        let (_, upstream_enthalpy) = self.contact_enthalpy(upstream_theta, pressure_pa)?;
        let residual = |theta: f64| -> Result<f64, NsfError> {
            let (rho, enthalpy) = self.contact_enthalpy(theta, pressure_pa)?;
            let k = self.spec.gas.face_conductivity(
                &self.law,
                0.5 * (adjacent_theta + theta),
                0.5 * (adjacent_pressure_pa + pressure_pa),
            )?;
            let value = (-rho * outward_velocity_lu) * (enthalpy - upstream_enthalpy)
                + 2.0 * k * (theta - adjacent_theta);
            if !(k.is_finite() && k > 0.0 && value.is_finite()) {
                return Err(NsfError::Input(
                    "characteristic thermal contact requires positive finite conductivity and finite flux",
                ));
            }
            Ok(value)
        };
        let mut low = adjacent_theta.min(upstream_theta);
        let mut high = adjacent_theta.max(upstream_theta);
        let mut low_value = residual(low)?;
        let mut high_value = residual(high)?;
        if low_value > 0.0 || high_value < 0.0 {
            return Err(NsfError::Input(
                "characteristic thermal contact is not bracketed",
            ));
        }
        if low_value == 0.0 {
            return Ok(low);
        }
        if high_value == 0.0 {
            return Ok(high);
        }
        // The admitted temperature interval is close to the reference chart.
        // Sixty-four bisections bound work and reach adjacent f64 values;
        // no blending width or authored tolerance changes the heat.
        for _ in 0..64 {
            let midpoint = low + 0.5 * (high - low);
            if midpoint == low || midpoint == high {
                return Ok(if low_value.abs() <= high_value.abs() {
                    low
                } else {
                    high
                });
            }
            let value = residual(midpoint)?;
            if value == 0.0 {
                return Ok(midpoint);
            }
            if value < 0.0 {
                low = midpoint;
                low_value = value;
            } else {
                high = midpoint;
                high_value = value;
            }
        }
        Err(NsfError::Input(
            "characteristic thermal contact exhausted its bounded solve",
        ))
    }
}

#[cfg(test)]
mod tests;
