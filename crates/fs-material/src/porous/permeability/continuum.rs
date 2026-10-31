//! Reuse the original intrinsic-permeability geometry with continuum fluids.
//!
//! This view adds no transport law. The consuming hydraulic owner supplies its
//! actual fluid density and viscosity and retains its source-backed state. The
//! original rigid skeleton, layer order, tensor axes and source support remain
//! borrowed here. A single fluid occupies the connected pores: calibrated gas
//! slip and gas-only relative mobility cannot become liquid properties.

use super::{DarcyGasBarrier, MaterialError, PorousGasSupport, invalid, positive};

/// A checked, borrowed view of the original no-slip, single-fluid pore geometry.
///
/// The original gas barrier also owns the fluid-independent intrinsic tensors
/// and series geometry. This view reuses those exact data without constructing
/// a gas state or treating the gas EOS as an incompressible-liquid law.
/// Permeability keeps the original barrier's explicit transverse constraint:
/// planar normal pressure gradient or fully developed confined normal flow.
/// This view never changes that geometry to suit a different fluid.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ContinuumDarcyGeometry<'a> {
    barrier: &'a DarcyGasBarrier,
}

impl DarcyGasBarrier {
    /// Borrow the existing continuum geometry without applying a gas-only law.
    ///
    /// # Errors
    /// Refuses every nonzero Klinkenberg calibration and every prescribed wet
    /// gas state, including a dry endpoint with unit gas relative permeability.
    /// Such states do not qualify a liquid or general-fluid saturation law.
    pub fn continuum_geometry(&self) -> Result<ContinuumDarcyGeometry<'_>, MaterialError> {
        if self.maximum_slip_pa != 0.0 || self.layers.iter().any(|layer| layer.wet_state.is_some())
        {
            return Err(invalid(
                "continuum geometry requires intrinsic permeability without gas slip or gas-only wet mobility",
            ));
        }
        Ok(ContinuumDarcyGeometry { barrier: self })
    }
}

impl<'a> ContinuumDarcyGeometry<'a> {
    /// Original layers, tensor frames, source receipts and declared support.
    pub const fn barrier(&self) -> &'a DarcyGasBarrier {
        self.barrier
    }

    /// Original sum(L/k_effective), in 1/m. None means an exactly sealed layer.
    ///
    /// Divide by the actual superficial flow area for the viscous geometric
    /// coefficient in dp = mu * sum(L/k_effective) * Q/area. The existing fluid
    /// transport and hydraulic owners supply mu and Q; no coefficient is frozen
    /// at the permeability source's original query temperature.
    pub const fn geometric_resistance_per_m(&self) -> Option<f64> {
        self.barrier.geometric_resistance_per_m
    }

    /// Original calibrated sum(beta_i*L_i), dimensionless, in
    /// dp = mu*sum(L/k)*Q/area + rho*sum(beta*L)*Q*|Q|/area².
    /// The original fluid owner supplies its actual density, viscosity and Q.
    /// A consumer must retain this term; zero is the exact Darcy limit.
    pub const fn inertial_geometric_resistance(&self) -> f64 {
        self.barrier.inertial_geometric_resistance()
    }

    /// Original model limits intersected with all retained source T/p limits.
    pub const fn support(&self) -> PorousGasSupport {
        self.barrier.support
    }

    /// Admit an actual local fluid state and signed superficial velocity.
    ///
    /// Both absolute face pressures must fit the original permeability support.
    /// The consumer owns its EOS/compressibility and temperature approximation;
    /// this admission does not replace those checks or the fluid's source limits.
    ///
    /// # Errors
    /// Refuses nonfinite/out-of-domain states, excessive superficial speed, and
    /// any nonzero velocity through a sealed layer. No epsilon opens a seal.
    pub fn admit_flow(
        &self,
        temperature_k: f64,
        left_pressure_pa: f64,
        right_pressure_pa: f64,
        superficial_speed_m_s: f64,
    ) -> Result<(), MaterialError> {
        self.admit_temperature_speed(temperature_k, superficial_speed_m_s)?;
        self.barrier
            .support
            .admit_state(temperature_k, &[left_pressure_pa, right_pressure_pa])
    }

    /// Admit local transport samples before a network has reconstructed pressure.
    ///
    /// This validates only the actual temperature, signed superficial speed and
    /// exact seal. A consumer must still call `admit_flow` with both actual
    /// absolute face pressures before accepting its complete hydraulic trial.
    /// No reference, gas or surrogate pressure is substituted here.
    ///
    /// # Errors
    /// Refuses out-of-domain temperature, invalid/excessive speed or nonzero
    /// velocity through a sealed layer.
    pub fn admit_temperature_speed(
        &self,
        temperature_k: f64,
        superficial_speed_m_s: f64,
    ) -> Result<(), MaterialError> {
        self.barrier.support.admit_state(temperature_k, &[])?;
        if !superficial_speed_m_s.is_finite()
            || superficial_speed_m_s.abs() > self.barrier.support.maximum_superficial_speed_m_s
            || self.barrier.geometric_resistance_per_m.is_none() && superficial_speed_m_s != 0.0
        {
            return Err(invalid(
                "continuum Darcy flow exceeds original speed support or crosses a sealed layer",
            ));
        }
        Ok(())
    }
}

impl PorousGasSupport {
    // Both the original compressible gas owner and continuum consumers enforce
    // these same original source coordinates. No GasState is fabricated.
    pub(super) fn admit_state(
        self,
        temperature_k: f64,
        absolute_pressures_pa: &[f64],
    ) -> Result<(), MaterialError> {
        if !positive(temperature_k)
            || temperature_k < self.temperature_k.0
            || temperature_k > self.temperature_k.1
            || absolute_pressures_pa.iter().any(|&pressure| {
                !pressure.is_finite()
                    || pressure < self.pressure_pa.0
                    || pressure > self.pressure_pa.1
            })
        {
            return Err(invalid(
                "permeability query lies outside source temperature/absolute-pressure support",
            ));
        }
        Ok(())
    }
}
