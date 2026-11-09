//! Shared prescribed stream mass flux from the original finite gas EOS state.
use super::{MaterialConvectionError, PistonGasState, invalid};

/// Velocity used by one convection trial, resolved from its bulk gas density.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasMassFlowState {
    pub mass_flow_kg_s: f64,
    /// Stream cross section, distinct from the exposed thermal area.
    pub flow_area_m2: f64,
    pub bulk_density_kg_m3: f64,
    pub speed_m_s: f64,
}

#[derive(Clone, Copy, Debug)]
pub(super) struct MassFlowInput {
    mass_flow_kg_s: f64,
    flow_area_m2: f64,
    mass_flux_kg_m2_s: f64,
}

impl MassFlowInput {
    pub(super) fn new(
        mass_flow_kg_s: f64,
        flow_area_m2: f64,
    ) -> Result<Self, MaterialConvectionError> {
        let flux = mass_flow_kg_s / flow_area_m2;
        if !mass_flow_kg_s.is_finite()
            || mass_flow_kg_s <= 0.0
            || !flow_area_m2.is_finite()
            || flow_area_m2 <= 0.0
            || !flux.is_finite()
            || flux <= 0.0
        {
            return Err(invalid(
                "finite-gas convection requires positive finite mass rate and flow area with representable mass flux",
            ));
        }
        Ok(Self {
            mass_flow_kg_s,
            flow_area_m2,
            mass_flux_kg_m2_s: flux,
        })
    }

    pub(super) const fn mass_flow_kg_s(self) -> f64 {
        self.mass_flow_kg_s
    }
    pub(super) const fn flow_area_m2(self) -> f64 {
        self.flow_area_m2
    }

    pub(super) fn admit_flow_area(self, expected: f64) -> Result<(), MaterialConvectionError> {
        if !expected.is_finite()
            || expected <= 0.0
            || (self.flow_area_m2 - expected).abs() > 64.0 * f64::EPSILON * expected
        {
            return Err(invalid(
                "prescribed duct mass flow requires the cross section of its actual circular or rectangular passage",
            ));
        }
        Ok(())
    }

    pub(super) fn admit_circular_area(
        self,
        diameter_m: f64,
    ) -> Result<(), MaterialConvectionError> {
        if !diameter_m.is_finite() || diameter_m <= 0.0 {
            return Err(invalid(
                "prescribed tube mass flow requires a positive inside diameter",
            ));
        }
        self.admit_flow_area(core::f64::consts::PI * diameter_m * diameter_m / 4.0)
    }

    pub(super) fn resolve(
        self,
        gas: &impl PistonGasState,
    ) -> Result<GasMassFlowState, MaterialConvectionError> {
        // Read actual EOS density directly. No film or wall transport sample
        // may replace this state, and unused bulk mu/k need not be queried.
        let density = gas.density_kg_m3();
        if !density.is_finite() || density <= 0.0 {
            return Err(invalid(
                "prescribed convection requires the original positive bulk density",
            ));
        }
        let speed = self.mass_flux_kg_m2_s / density;
        if !speed.is_finite() || speed <= 0.0 {
            return Err(invalid(
                "bulk-density convection velocity is unrepresentable",
            ));
        }
        Ok(GasMassFlowState {
            mass_flow_kg_s: self.mass_flow_kg_s,
            flow_area_m2: self.flow_area_m2,
            bulk_density_kg_m3: density,
            speed_m_s: speed,
        })
    }
}
