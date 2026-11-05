//! Original molecular reservoirs in the complete axial NSF field.
//!
//! For s=-i*omega, the original law supplies Cp(s)=Cv(s)+R. Only the energy
//! equation changes: Q'=A*s*(p-rho*Cp(s)*T'). Momentum, Fourier conduction,
//! ideal-gas continuity and the thermodynamic-pressure/traction map retain
//! their original owners. The complex caloric law and its equilibrium/frozen
//! limits are described by Petculescu & Lueptow, PRL 94, 238301 (2005),
//! doi:10.1103/PhysRevLett.94.238301, equations 2/4; their phasor convention is
//! conjugated here. No constant bulk viscosity or added absorption represents
//! these same reservoirs a second time.

use super::{
    BulkThermalCylinder, BulkThermalError, BulkThermoviscousEquilibrium, ThermalEndCondition,
    ThermalNetwork, poll,
};
use fs_material::relaxation::RelaxationState;

impl BulkThermalCylinder {
    /// Retain both acoustic and thermal spatial modes of an original relaxing
    /// ideal gas, with independent thermal end conditions and traction ports.
    ///
    /// The supplied state owns equilibrium thermodynamics, shear viscosity,
    /// conductivity and molecular heat capacities/times. Its complex Cp is
    /// evaluated once at this frequency. There is no additional bulk-viscosity,
    /// ISO or Stokes--Kirchhoff attenuation term. An empty/zero-strength law
    /// recovers the native classical cylinder with zero bulk viscosity exactly.
    ///
    /// This is homogeneous frozen-equilibrium, constant-area linear NSF with
    /// slip walls. It adds no lateral wall layer, mean flow, chemistry, changing
    /// material state or nonlinear heating. Original material validity remains
    /// the caller's responsibility; no measured-data accuracy is inferred.
    pub fn from_relaxation(
        state: &RelaxationState,
        radius_m: f64,
        length_m: f64,
        omega_rad_s: f64,
        ends: [ThermalEndCondition; 2],
        maximum_condition: f64,
    ) -> Result<Self, BulkThermalError> {
        Self::from_relaxation_controlled(
            state,
            radius_m,
            length_m,
            omega_rad_s,
            ends,
            maximum_condition,
            || false,
        )
    }

    /// Same original constitutive preparation, with cancellation around the
    /// bounded exponential and inside the native scattering doublings. Neither
    /// a rejected material state nor cancellation publishes a partial cylinder.
    pub fn from_relaxation_controlled(
        state: &RelaxationState,
        radius_m: f64,
        length_m: f64,
        omega_rad_s: f64,
        ends: [ThermalEndCondition; 2],
        maximum_condition: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, BulkThermalError> {
        poll(&mut cancelled, 0)?;
        let equilibrium = BulkThermoviscousEquilibrium::from_material(state.gas(), 0.0)?;
        let cp = if state
            .processes()
            .iter()
            .all(|p| p.heat_capacity_j_kg_k == 0.0)
        {
            None
        } else {
            Some(state.caloric_response(omega_rad_s)?.cp_j_kg_k)
        };
        let network = ThermalNetwork::with_caloric_response(
            &equilibrium,
            radius_m,
            length_m,
            omega_rad_s,
            maximum_condition,
            cp,
            &mut cancelled,
        )?;
        let prepared = Self::from_network(network, ends, maximum_condition)?;
        poll(&mut cancelled, 1)?;
        Ok(prepared)
    }
}

#[cfg(test)]
mod tests;
