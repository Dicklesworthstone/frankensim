//! One original source fluid across pressure and resolved thermoviscous regions.
//!
//! The pressure FE uses the source's actual density and adiabatic sound speed;
//! the resolved owner retains its original fluid/solid material sources and
//! full thermal expansion/compressibility closure. Gases and liquids meet
//! through the same geometry-verified modal traction/flow interface. Neither
//! sound speed nor thermal expansion is replaced by an ideal-gas assumption.
//! The pressure region retains its adiabatic pressure approximation and any
//! explicitly declared thin-wall law; full viscous/thermal fields belong to
//! the resolved region. The source frequency band remains explicit.
//!
//! An original molecular source is evaluated once per frequency. Its complex
//! bulk modulus enters pressure-region compression and its cp enters the
//! resolved gas heat equation, with both loss channels independently reported.
//! Molecular pressure regions currently require insulated/slip acoustic walls:
//! their approximate thermal boundary-layer overlay is not combined with this
//! bulk law. Finite solids in the resolved region remain fully supported.
//!
//! Fields use exp(+i omega t); exported duct scattering is conjugated once to
//! exp(-i omega t). Retained spatial modes, real power references and thermal
//! interface/placement limitations remain those of the original hybrid owner.

use super::thermoviscous::ThermoviscousInterior;
use fs_duct::{DuctError, compose::TwoPort};
use fs_feec::acoustics::helmholtz::multiport::InteriorMultiport;
use fs_feec::acoustics::hybrid::{HybridError, HybridInterior as SpatialHybrid, PreparedHybrid};
use fs_material::fluid_acoustics::thermoviscous::ThermoviscousFluid;
use fs_material::{MaterialError, gas::GasState};
use fs_math::c64::C64;

/// Source admission and original finite-element/duct failures remain distinct.
#[derive(Debug)]
pub enum HybridSourceError {
    InvalidInput(&'static str),
    Material(MaterialError),
    Hybrid(HybridError),
    Duct(DuctError),
}
impl core::fmt::Display for HybridSourceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidInput(what) => write!(f, "source hybrid acoustics: {what}"),
            Self::Material(error) => error.fmt(f),
            Self::Hybrid(error) => error.fmt(f),
            Self::Duct(error) => error.fmt(f),
        }
    }
}
impl std::error::Error for HybridSourceError {}
impl From<MaterialError> for HybridSourceError {
    fn from(error: MaterialError) -> Self {
        Self::Material(error)
    }
}
impl From<HybridError> for HybridSourceError {
    fn from(error: HybridError) -> Self {
        Self::Hybrid(error)
    }
}
impl From<DuctError> for HybridSourceError {
    fn from(error: DuctError) -> Self {
        Self::Duct(error)
    }
}

/// Borrowed source materials and original spatial solvers. No acoustic field,
/// finite wall, retained mode or source identity is replaced by this adapter.
pub struct HybridInterior<'a, 'c, S = GasState> {
    source: &'a ThermoviscousInterior<'c, S>,
    model: SpatialHybrid<'a, 'c>,
    helmholtz_port_count: usize,
}

impl<'a, 'c, S: ThermoviscousFluid> HybridInterior<'a, 'c, S> {
    /// Bind the same original equilibrium fluid to both owners. Density,
    /// adiabatic sound speed and any pressure-wall transport must match the
    /// original material snapshot exactly. The resolved source has already
    /// admitted the full EOS derivatives through its own constructor.
    /// Selected interface modes are verified from their physical traces.
    /// External ordering is unselected pressure modes then unselected resolved
    /// modes, preserving each owner's original order.
    pub fn new(
        helmholtz: &'a InteriorMultiport,
        source: &'a ThermoviscousInterior<'c, S>,
        helmholtz_modes: &[usize],
        thermoviscous_modes: &[usize],
        cancelled: &dyn Fn() -> bool,
    ) -> Result<Self, HybridSourceError> {
        if cancelled() {
            return Err(HybridError::Cancelled.into());
        }
        if source.relaxation().is_some() && helmholtz.boundary_layer_transport().is_some() {
            return Err(HybridSourceError::InvalidInput(
                "molecular pressure region cannot also use an independent boundary-layer overlay",
            ));
        }
        let fluid = source.material_source().thermoviscous_state()?;
        let acoustic = fluid.acoustic();
        let transport = fluid.transport();
        let medium = helmholtz.medium();
        if medium.density_kg_m3.to_bits() != acoustic.density_kg_m3().to_bits()
            || medium.sound_speed_m_s.to_bits() != acoustic.sound_speed_m_s().to_bits()
        {
            return Err(HybridSourceError::InvalidInput(
                "hybrid pressure region differs from the original source fluid",
            ));
        }
        let resolved = source.inner().medium();
        for (actual, expected) in [
            (resolved.density_kg_m3, acoustic.density_kg_m3()),
            (resolved.cp_j_kg_k, transport.specific_heat_j_kg_k()),
            (
                resolved.dynamic_viscosity_pa_s,
                transport.dynamic_viscosity_pa_s(),
            ),
            (
                resolved.thermal_conductivity_w_m_k,
                transport.thermal_conductivity_w_m_k(),
            ),
        ] {
            if actual.to_bits() != expected.to_bits() {
                return Err(HybridSourceError::InvalidInput(
                    "hybrid resolved region differs from the original source fluid",
                ));
            }
        }
        if let Some(wall) = helmholtz.boundary_layer_transport() {
            for (actual, expected) in [
                (
                    wall.dynamic_viscosity_pa_s,
                    transport.dynamic_viscosity_pa_s(),
                ),
                (
                    wall.thermal_conductivity_w_m_k,
                    transport.thermal_conductivity_w_m_k(),
                ),
                (wall.cp_j_kg_k, transport.specific_heat_j_kg_k()),
                (wall.heat_capacity_ratio, fluid.heat_capacity_ratio()),
            ] {
                if actual.to_bits() != expected.to_bits() {
                    return Err(HybridSourceError::InvalidInput(
                        "hybrid pressure-wall transport differs from the original source fluid",
                    ));
                }
            }
        }
        let model = SpatialHybrid::new(
            helmholtz,
            source.inner(),
            helmholtz_modes,
            thermoviscous_modes,
            cancelled,
        )?;
        Ok(Self {
            source,
            model,
            helmholtz_port_count: helmholtz.port_mode_count(),
        })
    }
}

impl<'a, 'c, S> HybridInterior<'a, 'c, S> {
    /// Original fluid source and complete resolved material/field owner.
    pub fn source(&self) -> &ThermoviscousInterior<'c, S> {
        self.source
    }

    pub fn inner(&self) -> &SpatialHybrid<'a, 'c> {
        &self.model
    }

    /// Original loaded volume factors and bounded interface factor, reused for
    /// all simultaneous excitations and scattering columns at this frequency.
    pub fn prepare(
        &self,
        frequency_hz: f64,
        helmholtz_reference: &[f64],
        thermoviscous_reference: &[f64],
        cancelled: &dyn Fn() -> bool,
    ) -> Result<PreparedHybrid<'_, 'c>, HybridSourceError> {
        if cancelled() {
            return Err(HybridError::Cancelled.into());
        }
        let band = self.source.band_hz();
        if !frequency_hz.is_finite() || frequency_hz < band[0] || frequency_hz > band[1] {
            return Err(HybridSourceError::InvalidInput(
                "hybrid source frequency outside declared band",
            ));
        }
        match self.source.relaxation() {
            Some(source) => {
                let caloric = source.caloric_response(core::f64::consts::TAU * frequency_hz)?;
                // Original material response is exp(-iwt); conjugate both
                // coefficients once for the shared FE exp(+iwt) convention.
                Ok(self.model.prepare_with_bulk_modulus(
                    frequency_hz,
                    caloric.cp_j_kg_k.conj(),
                    caloric.bulk_modulus_pa.conj(),
                    helmholtz_reference,
                    thermoviscous_reference,
                    cancelled,
                )?)
            }
            None => Ok(self.model.prepare(
                frequency_hz,
                C64::from_re(self.source.inner().medium().cp_j_kg_k),
                helmholtz_reference,
                thermoviscous_reference,
                cancelled,
            )?),
        }
    }

    /// Complete external power scattering in the original duct convention.
    /// All source-bound spatial modes and cross-couplings remain present.
    pub fn duct_scattering(
        &self,
        frequency_hz: f64,
        helmholtz_reference: &[f64],
        thermoviscous_reference: &[f64],
        cancelled: &dyn Fn() -> bool,
    ) -> Result<Vec<C64>, HybridSourceError> {
        let prepared = self.prepare(
            frequency_hz,
            helmholtz_reference,
            thermoviscous_reference,
            cancelled,
        )?;
        let mut scattering = prepared.scattering(cancelled)?;
        for row in scattering.chunks_mut(self.model.external_port_count()) {
            if cancelled() {
                return Err(HybridSourceError::Hybrid(HybridError::Cancelled));
            }
            for entry in row {
                *entry = entry.conj();
            }
        }
        Ok(scattering)
    }

    /// Exactly two external work coordinates in the original duct TwoPort.
    /// A common real resistance normalizes every port; the original hybrid
    /// interface accounts for its actual spatial dual map. No mode is dropped.
    pub fn duct_two_port(
        &self,
        frequency_hz: f64,
        reference_impedance: f64,
        cancelled: &dyn Fn() -> bool,
    ) -> Result<TwoPort, HybridSourceError> {
        if self.model.external_port_count() != 2 {
            return Err(HybridSourceError::InvalidInput(
                "hybrid duct TwoPort requires exactly two external coordinates",
            ));
        }
        let h = vec![reference_impedance; self.helmholtz_port_count];
        let l = vec![reference_impedance; self.source.inner().port_mode_count()];
        let s = self.duct_scattering(frequency_hz, &h, &l, cancelled)?;
        Ok(TwoPort::new([s[0], s[1], s[2], s[3]], reference_impedance)?)
    }
}

#[cfg(test)]
mod tests;
