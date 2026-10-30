//! One original source gas across pressure and resolved thermoviscous regions.
//!
//! The pressure FE uses the source gas's equilibrium density and sound speed;
//! the resolved owner retains its original gas/solid material sources. Both
//! regions meet through the original geometry-verified modal traction/flow
//! interface. The source frequency band remains explicit.
//!
//! This adapter currently admits calorically perfect source gas only. A
//! molecular source would require its dispersive bulk response in the pressure
//! region as well as cp in the resolved region; applying it to only one side
//! would silently change the physical model and therefore refuses here.
//!
//! Fields use exp(+i omega t); exported duct scattering is conjugated once to
//! exp(-i omega t). Retained spatial modes, real power references and thermal
//! interface/placement limitations remain those of the original hybrid owner.

use super::thermoviscous::ThermoviscousInterior;
use fs_duct::{DuctError, compose::TwoPort};
use fs_feec::acoustics::helmholtz::multiport::InteriorMultiport;
use fs_feec::acoustics::hybrid::{HybridError, HybridInterior as SpatialHybrid, PreparedHybrid};
use fs_math::c64::C64;

/// Source admission and original finite-element/duct failures remain distinct.
#[derive(Debug)]
pub enum HybridSourceError {
    InvalidInput(&'static str),
    Hybrid(HybridError),
    Duct(DuctError),
}
impl core::fmt::Display for HybridSourceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidInput(what) => write!(f, "source hybrid acoustics: {what}"),
            Self::Hybrid(error) => error.fmt(f),
            Self::Duct(error) => error.fmt(f),
        }
    }
}
impl std::error::Error for HybridSourceError {}
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
pub struct HybridInterior<'a, 'c> {
    source: &'a ThermoviscousInterior<'c>,
    model: SpatialHybrid<'a, 'c>,
    helmholtz_port_count: usize,
}

impl<'a, 'c> HybridInterior<'a, 'c> {
    /// Bind the same original equilibrium gas to both owners. Selected
    /// interface modes are verified from their original physical traces.
    /// External ordering is unselected pressure modes then unselected resolved
    /// modes, preserving each owner's original order.
    pub fn new(
        helmholtz: &'a InteriorMultiport,
        source: &'a ThermoviscousInterior<'c>,
        helmholtz_modes: &[usize],
        thermoviscous_modes: &[usize],
        cancelled: &dyn Fn() -> bool,
    ) -> Result<Self, HybridSourceError> {
        if source.relaxation().is_some() {
            return Err(HybridSourceError::InvalidInput(
                "hybrid molecular gas requires dispersive compression in both regions",
            ));
        }
        let medium = helmholtz.medium();
        let gas = source.gas();
        if medium.density_kg_m3.to_bits() != gas.density.to_bits()
            || medium.sound_speed_m_s.to_bits() != gas.sound_speed.to_bits()
        {
            return Err(HybridSourceError::InvalidInput(
                "hybrid pressure region differs from the original source gas",
            ));
        }
        if let Some(transport) = helmholtz.boundary_layer_transport() {
            for (actual, expected) in [
                (transport.dynamic_viscosity_pa_s, gas.dynamic_viscosity),
                (
                    transport.thermal_conductivity_w_m_k,
                    gas.thermal_conductivity,
                ),
                (transport.cp_j_kg_k, gas.specific_heat_cp),
                (transport.heat_capacity_ratio, gas.gamma),
            ] {
                if actual.to_bits() != expected.to_bits() {
                    return Err(HybridSourceError::InvalidInput(
                        "hybrid pressure-wall transport differs from the original source gas",
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

    pub fn source(&self) -> &ThermoviscousInterior<'c> {
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
        let band = self.source.band_hz();
        if !frequency_hz.is_finite() || frequency_hz < band[0] || frequency_hz > band[1] {
            return Err(HybridSourceError::InvalidInput(
                "hybrid source frequency outside declared band",
            ));
        }
        Ok(self.model.prepare(
            frequency_hz,
            C64::from_re(self.source.gas().specific_heat_cp),
            helmholtz_reference,
            thermoviscous_reference,
            cancelled,
        )?)
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
