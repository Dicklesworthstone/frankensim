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
    BulkThermalCylinder, BulkThermalError, BulkThermoviscousEquilibrium, ComposeOptions, DuctError,
    Segment, ThermalEndCondition, ThermalNetwork, bad, poll,
};
use fs_material::relaxation::RelaxationState;

/// One original relaxing ideal-gas region in left-to-right physical order.
/// The borrowed source owns its gas, transport and complete reservoir spectrum.
#[derive(Clone, Copy, Debug)]
pub struct BulkThermalRelaxationRegion<'a> {
    pub material: &'a RelaxationState,
    /// Positive finite physical axial length [m].
    pub length_m: f64,
}

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
        let network = ThermalNetwork::from_relaxation(
            state,
            radius_m,
            length_m,
            omega_rad_s,
            maximum_condition,
            &mut cancelled,
        )?;
        let prepared = Self::from_network(network, ends, maximum_condition)?;
        poll(&mut cancelled, 1)?;
        Ok(prepared)
    }

    /// Compose complete NSF fields through independently relaxing ideal gases.
    /// Each original region supplies its own density, shear viscosity,
    /// conductivity and complex Cp at this frequency; no coefficients or
    /// reservoirs are averaged. An empty spectrum uses the original classical
    /// zero-bulk-viscosity region exactly.
    ///
    /// Ideal internal contacts preserve total traction, volume flow, physical
    /// temperature perturbation and heat power. Only the TWO EXTERNAL ends
    /// receive the specified thermal conditions. The first region's positive
    /// real power references serve both outer ports; each endpoint retains its
    /// own thermodynamic-pressure map. Pressure itself need not be continuous
    /// across different longitudinal viscosities.
    ///
    /// All regions must share the same equilibrium temperature. These are
    /// frozen, prescribed, constant-area regions with ideal contacts and slip
    /// walls; the caller establishes mechanical and phase compatibility.
    /// Different source compositions or reservoir identities are supported,
    /// but no mixing, material diffusion, stationary temperature-gradient
    /// source, contact resistance, lateral wall layers or mean flow is inferred.
    /// `maximum_blocks` bounds the complete authored chain before preparation.
    pub fn from_relaxation_regions(
        regions: &[BulkThermalRelaxationRegion<'_>],
        radius_m: f64,
        omega_rad_s: f64,
        ends: [ThermalEndCondition; 2],
        options: ComposeOptions,
    ) -> Result<Self, BulkThermalError> {
        Self::from_relaxation_regions_controlled(
            regions,
            radius_m,
            omega_rad_s,
            ends,
            options,
            || false,
        )
    }

    /// Same complete composition, polling inside each original exponential
    /// preparation and scattering doubling and at every material contact.
    /// No partial chain or thermally terminated candidate escapes a refusal.
    pub fn from_relaxation_regions_controlled(
        regions: &[BulkThermalRelaxationRegion<'_>],
        radius_m: f64,
        omega_rad_s: f64,
        ends: [ThermalEndCondition; 2],
        options: ComposeOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, BulkThermalError> {
        options.validate()?;
        let first = regions.first().ok_or(DuctError::EmptyDuct)?;
        if regions.len() > options.maximum_blocks {
            return Err(bad("regional molecular cylinder exceeds its declared region cap").into());
        }
        let temperature = first.material.gas().temperature;
        for (index, region) in regions.iter().enumerate() {
            poll(&mut cancelled, index)?;
            Segment::Cylinder {
                radius: radius_m,
                length: region.length_m,
            }
            .validate()?;
            if region.material.gas().temperature != temperature {
                return Err(bad(
                    "regional molecular contacts require one common equilibrium temperature",
                )
                .into());
            }
        }
        let mut network = ThermalNetwork::from_relaxation(
            first.material,
            radius_m,
            first.length_m,
            omega_rad_s,
            options.maximum_condition,
            &mut cancelled,
        )?;
        for region in &regions[1..] {
            let right = ThermalNetwork::from_relaxation(
                region.material,
                radius_m,
                region.length_m,
                omega_rad_s,
                options.maximum_condition,
                &mut cancelled,
            )?;
            network.append(right, options.maximum_condition, &mut cancelled)?;
        }
        let prepared = Self::from_network(network, ends, options.maximum_condition)?;
        poll(&mut cancelled, regions.len())?;
        Ok(prepared)
    }
}

impl ThermalNetwork {
    fn from_relaxation(
        state: &RelaxationState,
        radius_m: f64,
        length_m: f64,
        omega_rad_s: f64,
        maximum_condition: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, BulkThermalError> {
        poll(cancelled, 0)?;
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
        Self::with_caloric_response(
            &equilibrium,
            radius_m,
            length_m,
            omega_rad_s,
            maximum_condition,
            cp,
            cancelled,
        )
    }
}

#[cfg(test)]
mod tests;
