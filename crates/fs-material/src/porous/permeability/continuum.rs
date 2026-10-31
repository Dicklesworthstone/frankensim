//! Reuse the original intrinsic-permeability geometry with continuum fluids.
//!
//! This view adds no transport law. The consuming hydraulic owner supplies its
//! actual fluid density and viscosity and retains its source-backed state. The
//! original rigid skeleton, layer order, tensor axes and source support remain
//! borrowed here. Full single-fluid occupation or an explicitly calibrated frozen
//! liquid saturation is supported; gas slip and gas-only relative mobility
//! cannot become liquid properties.

use super::{DarcyGasBarrier, MaterialError, PorousGasSupport, invalid, positive};
use crate::porous::wet::liquid::WetPoreLiquidState;

/// A checked, borrowed view of the original no-slip continuum pore geometry.
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
    liquid_states: Option<&'a [WetPoreLiquidState]>,
    geometric_resistance_per_m: Option<f64>,
    support: PorousGasSupport,
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
        Ok(ContinuumDarcyGeometry {
            barrier: self,
            liquid_states: None,
            geometric_resistance_per_m: self.geometric_resistance_per_m,
            support: self.support,
        })
    }
}

impl<'a> ContinuumDarcyGeometry<'a> {
    /// Bind one explicitly calibrated liquid state to each original layer.
    /// The scalar kr_l changes liquid mobility, not intrinsic k or its axes.
    /// Rebinding starts from the original geometry, so updates never compound.
    /// No capillary pressure is added to hydraulic drag; this frozen snapshot
    /// neither solves gas pressure nor advances saturation or liquid storage.
    ///
    /// # Errors
    /// Refuses an incomplete state list, disjoint support, unrepresentable
    /// positive resistance, or any positive original Forchheimer calibration:
    /// a fully filled inertial law does not qualify partial-saturation inertia.
    pub fn with_liquid_states<'b>(
        self,
        states: &'b [WetPoreLiquidState],
    ) -> Result<ContinuumDarcyGeometry<'b>, MaterialError>
    where
        'a: 'b,
    {
        if states.len() != self.barrier.layers.len() {
            return Err(invalid(
                "liquid mobility requires exactly one state per original porous layer",
            ));
        }
        if self.barrier.inertial_geometric_resistance() != 0.0 {
            return Err(invalid(
                "original Forchheimer inertia is not calibrated for prescribed liquid saturation",
            ));
        }
        let mut support = self.barrier.support;
        let mut sealed = self.barrier.geometric_resistance_per_m.is_none();
        for (layer, state) in self.barrier.layers.iter().zip(states) {
            state.law().admit_layer_source(layer)?;
            let active = state.support().flow_support();
            support.temperature_k = (
                support.temperature_k.0.max(active.temperature_k.0),
                support.temperature_k.1.min(active.temperature_k.1),
            );
            support.pressure_pa = (
                support.pressure_pa.0.max(active.pressure_pa.0),
                support.pressure_pa.1.min(active.pressure_pa.1),
            );
            support.maximum_superficial_speed_m_s = support
                .maximum_superficial_speed_m_s
                .min(active.maximum_superficial_speed_m_s);
            sealed |= state.liquid_relative_permeability() == 0.0;
        }
        support.validate()?;
        let mut resistance = 0.0;
        if !sealed {
            for (layer, state) in self.barrier.layers.iter().zip(states) {
                let intrinsic = self
                    .barrier
                    .flow_constraint
                    .permeability(&layer.permeability, self.barrier.normal)?;
                let part = (layer.thickness_m / intrinsic) / state.liquid_relative_permeability();
                if !positive(part) || !positive(resistance + part) {
                    return Err(invalid(
                        "positive unsaturated liquid series resistance is unrepresentable",
                    ));
                }
                resistance += part;
            }
        }
        Ok(ContinuumDarcyGeometry {
            barrier: self.barrier,
            liquid_states: Some(states),
            geometric_resistance_per_m: if sealed { None } else { Some(resistance) },
            support,
        })
    }

    /// Frozen liquid samples with their complete retained calibrations, if any.
    pub const fn liquid_states(&self) -> Option<&'a [WetPoreLiquidState]> {
        self.liquid_states
    }
    /// Original layers, tensor frames, source receipts and declared support.
    pub const fn barrier(&self) -> &'a DarcyGasBarrier {
        self.barrier
    }

    /// Sum(L/(k_effective*kr_l)), in 1/m; kr_l=1 without liquid states.
    /// None means an exactly sealed intrinsic or liquid-mobility layer.
    ///
    /// Divide by the actual superficial flow area for the viscous geometric
    /// coefficient in dp = mu * sum(L/k_effective) * Q/area. The existing fluid
    /// transport and hydraulic owners supply mu and Q; no coefficient is frozen
    /// at the permeability source's original query temperature.
    pub const fn geometric_resistance_per_m(&self) -> Option<f64> {
        self.geometric_resistance_per_m
    }

    /// Original calibrated sum(beta_i*L_i), dimensionless, in
    /// dp = mu*sum(L/k)*Q/area + rho*sum(beta*L)*Q*|Q|/area².
    /// The original fluid owner supplies its actual density, viscosity and Q.
    /// A consumer must retain this term; zero is the exact Darcy limit.
    pub const fn inertial_geometric_resistance(&self) -> f64 {
        self.barrier.inertial_geometric_resistance()
    }

    /// Original source limits intersected with every active liquid calibration.
    pub const fn support(&self) -> PorousGasSupport {
        self.support
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
        self.support
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
        self.support.admit_state(temperature_k, &[])?;
        if !superficial_speed_m_s.is_finite()
            || superficial_speed_m_s.abs() > self.support.maximum_superficial_speed_m_s
            || self.geometric_resistance_per_m.is_none() && superficial_speed_m_s != 0.0
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
