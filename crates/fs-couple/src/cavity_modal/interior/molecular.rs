//! Source-owned molecular calorics in the original spatial pressure/heat FE.
//!
//! The original RelaxationState supplies cp at each frequency. Its exp(-iwt)
//! response is conjugated exactly once for the FE exp(+iwt) convention. Both
//! molecular bulk compression and thermal diffusion then follow from that one
//! heat equation; total plane-wave absorption is never a bulk-loss input.

use super::{
    AcousticFeError, BoundaryLayerGeometry, ConductiveThermalWall, HelmholtzFe, HelmholtzModes,
    ThermalTransportSample, ThermalWall, ThermalWallDepth, boundary_layer_transport,
    conductive_wall_transport,
};
use fs_duct::compose::TwoPort;
use fs_feec::acoustics::helmholtz::molecular::{
    CaloricGas, CaloricInterior, PreparedCaloricInterior,
};
use fs_feec::acoustics::helmholtz::multiport::BoundaryMode;
use fs_material::{MaterialError, relaxation::RelaxationState};
use fs_math::c64::C64;

/// The thermal wall source; conductivity direction and caloric constraint were
/// selected by the original prepared material owner before constructing this.
pub enum MolecularThermalWall<'a> {
    Isothermal,
    Adiabatic,
    Conductive {
        sample: &'a ThermalTransportSample,
        depth: ThermalWallDepth,
        minimum_tangential_scale_m: f64,
    },
}

/// Original material or spatial FE refusal, without discarding its cause.
#[derive(Debug)]
pub enum MolecularInteriorError {
    Material(MaterialError),
    FiniteElement(AcousticFeError),
}
impl core::fmt::Display for MolecularInteriorError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "molecular interior: {self:?}")
    }
}
impl std::error::Error for MolecularInteriorError {}
impl From<MaterialError> for MolecularInteriorError {
    fn from(error: MaterialError) -> Self {
        Self::Material(error)
    }
}
impl From<AcousticFeError> for MolecularInteriorError {
    fn from(error: AcousticFeError) -> Self {
        Self::FiniteElement(error)
    }
}

/// Uniform frozen-mean source gas, spatial pressure and resolved heat field.
/// Pressure, heat and wall work retain their original FE coordinates. Optional
/// viscosity remains the original admitted thin-wall approximation, not LNS.
pub struct MolecularInterior {
    source: RelaxationState,
    wall_source: Option<ThermalTransportSample>,
    model: CaloricInterior,
}

impl MolecularInterior {
    /// Bind original source calorics and transport to the SAME FE volume.
    /// Source and FE rho/c must agree bitwise; conductive walls also require
    /// one exact common equilibrium temperature. Caller-provided pressure
    /// basis and full thermal FE dimension remain explicit numerical choices.
    pub fn new(
        fe: &HelmholtzFe<'_>,
        modes: &HelmholtzModes,
        ports: &[BoundaryMode],
        wall_faces: &[[u32; 3]],
        source: RelaxationState,
        wall: MolecularThermalWall<'_>,
        geometry: BoundaryLayerGeometry,
        cancelled: &dyn Fn() -> bool,
    ) -> Result<Self, MolecularInteriorError> {
        let gas = source.gas();
        let transport = boundary_layer_transport(gas)?;
        let medium = fe.medium();
        if medium.density_kg_m3.to_bits() != gas.density.to_bits()
            || medium.sound_speed_m_s.to_bits() != gas.sound_speed.to_bits()
        {
            return Err(AcousticFeError::InvalidInput {
                what: "molecular interior source differs from original FE medium",
            }
            .into());
        }
        let (wall, wall_source) = match wall {
            MolecularThermalWall::Isothermal => (ThermalWall::Isothermal, None),
            MolecularThermalWall::Adiabatic => (ThermalWall::Adiabatic, None),
            MolecularThermalWall::Conductive {
                sample,
                depth,
                minimum_tangential_scale_m,
            } => (
                ThermalWall::Conductive(ConductiveThermalWall {
                    transport: conductive_wall_transport(gas, sample)?,
                    depth,
                    minimum_tangential_scale_m,
                }),
                Some(*sample),
            ),
        };
        let model = CaloricInterior::new(
            fe,
            modes,
            ports,
            wall_faces,
            CaloricGas {
                pressure_pa: gas.pressure,
                temperature_k: gas.temperature,
                transport,
            },
            wall,
            geometry,
            cancelled,
        )?;
        Ok(Self {
            source,
            wall_source,
            model,
        })
    }

    #[must_use]
    pub fn source(&self) -> &RelaxationState {
        &self.source
    }
    #[must_use]
    pub fn wall_source(&self) -> Option<&ThermalTransportSample> {
        self.wall_source.as_ref()
    }
    #[must_use]
    pub fn port_mode_count(&self) -> usize {
        self.model.port_mode_count()
    }

    fn cp(&self, frequency_hz: f64) -> Result<C64, MolecularInteriorError> {
        Ok(self
            .source
            .caloric_response(core::f64::consts::TAU * frequency_hz)?
            .cp_j_kg_k
            .conj())
    }

    /// One original caloric evaluation and two original factorizations prepare
    /// every simultaneous excitation; returned response includes positive loss
    /// terms and the full original temperature field.
    pub fn prepare(
        &self,
        frequency_hz: f64,
        cancelled: &dyn Fn() -> bool,
    ) -> Result<PreparedCaloricInterior<'_>, MolecularInteriorError> {
        Ok(self
            .model
            .prepare(frequency_hz, self.cp(frequency_hz)?, cancelled)?)
    }

    /// Full row-major real-reference scattering in the duct exp(-iwt) convention.
    /// Modal work coordinates and inward flow orientation do not change.
    pub fn duct_scattering(
        &self,
        frequency_hz: f64,
        reference_resistance: &[f64],
        cancelled: &dyn Fn() -> bool,
    ) -> Result<Vec<C64>, MolecularInteriorError> {
        let mut s = self.model.scattering(
            frequency_hz,
            self.cp(frequency_hz)?,
            reference_resistance,
            cancelled,
        )?;
        for z in &mut s {
            *z = z.conj();
        }
        Ok(s)
    }

    /// Exactly two retained coordinates in the existing duct composition owner.
    /// No transverse mode is silently dropped and no surrogate loss is fitted.
    pub fn duct_two_port(
        &self,
        frequency_hz: f64,
        reference: f64,
        cancelled: &dyn Fn() -> bool,
    ) -> Result<TwoPort, MolecularInteriorError> {
        if self.port_mode_count() != 2 {
            return Err(AcousticFeError::InvalidInput {
                what: "molecular duct TwoPort requires exactly two retained FE modes",
            }
            .into());
        }
        let s = self.duct_scattering(frequency_hz, &[reference; 2], cancelled)?;
        TwoPort::new([s[0], s[1], s[2], s[3]], reference).map_err(|_| {
            AcousticFeError::InvalidInput {
                what: "molecular interior duct scattering block",
            }
            .into()
        })
    }
}

#[cfg(test)]
mod tests;
