//! Caloric linearizations on the existing conservative energy face stencil.

use super::{
    ChtError, ConvectionScheme, EnergyConfig, FlowField, ThermalFace, ThermalSetup, VoxelDomain,
};
use fs_exec::CancelGate;

/// A frozen local caloric tangent, in SI units: `h(T) = h0 + Cp * (T - T0)`.
///
/// The material owner supplies the value and slope from ONE admitted enthalpy
/// chart, then iterates and checks the true law at the accepted temperature.
/// This is a numerical coefficient image, not a material model or source claim.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AffineSpecificEnthalpy {
    /// Linearization temperature, positive K.
    pub temperature_k: f64,
    /// Specific enthalpy at that temperature, J/kg, using the common datum.
    pub enthalpy_j_kg: f64,
    /// Positive local chart slope, J/(kg K).
    pub heat_capacity_j_kg_k: f64,
}

impl AffineSpecificEnthalpy {
    /// Evaluate the frozen tangent in centered form. This does not evaluate the
    /// nonlinear material chart or establish that `temperature_k` is admissible.
    #[must_use]
    pub fn at_temperature(self, temperature_k: f64) -> f64 {
        self.heat_capacity_j_kg_k
            .mul_add(temperature_k - self.temperature_k, self.enthalpy_j_kg)
    }

    pub(super) fn intercept(self) -> f64 {
        self.heat_capacity_j_kg_k
            .mul_add(-self.temperature_k, self.enthalpy_j_kg)
    }
}

/// One constant-density fluid's frozen caloric field for first-order upwind
/// transport. Interior faces use the SAME upstream tangent on both sides, so
/// `rho * Q * h_upstream` cancels between neighboring control volumes even when
/// their Cp values differ. Molecular conduction remains the existing stencil.
///
/// All entries and boundary values must share one enthalpy datum. The original
/// `FluidProperties` supplies density; its reference Cp is not consumed by this
/// transport path. The supplied flow must conserve volume locally. No species,
/// changing density, compressibility, pressure work or phase transfer is added.
#[derive(Debug, Clone, Copy)]
pub struct AdvectiveEnthalpyField<'a> {
    /// Exactly one entry per cell: `Some` for fluid and `None` for solid.
    pub cells: &'a [Option<AffineSpecificEnthalpy>],
    /// Exact inlet/backflow enthalpy at each open face's declared temperature;
    /// `None` on closed faces. Face order is x-, x+, y-, y+, z-, z+.
    pub boundary_j_kg: [Option<f64>; 6],
}

impl AdvectiveEnthalpyField<'_> {
    pub(super) fn validate(
        &self,
        domain: &VoxelDomain,
        flow: &FlowField,
        setup: &ThermalSetup,
        config: &EnergyConfig,
        gate: &CancelGate,
    ) -> Result<(), ChtError> {
        let invalid = |reason: &str| ChtError::InvalidInput {
            field: "thermal.advective_enthalpy",
            reason: reason.into(),
        };
        if config.scheme != ConvectionScheme::Upwind {
            return Err(invalid(
                "sourced enthalpy transport requires first-order upwind",
            ));
        }
        if self.cells.len() != domain.cell_count() || flow.dims() != domain.dims() {
            return Err(invalid(
                "enthalpy and flow fields must match the thermal domain",
            ));
        }
        for (cell, value) in self.cells.iter().enumerate() {
            if cell.is_multiple_of(4096) {
                super::super::poll(gate)?;
            }
            match (domain.is_fluid(cell), value) {
                (true, Some(h))
                    if h.temperature_k.is_finite()
                        && h.temperature_k > 0.0
                        && h.enthalpy_j_kg.is_finite()
                        && h.heat_capacity_j_kg_k.is_finite()
                        && h.heat_capacity_j_kg_k > 0.0
                        && h.intercept().is_finite() => {}
                (false, None) => {}
                _ => {
                    return Err(invalid(
                        "need finite enthalpy, positive temperature/Cp on each fluid cell, and no solid entries",
                    ));
                }
            }
            for f in 0..6 {
                if flow.outward(domain, cell, f) != 0.0
                    && (!domain.is_fluid(cell)
                        || domain
                            .neighbor(cell, f)
                            .is_some_and(|n| !domain.is_fluid(n)))
                {
                    return Err(invalid("advected fluid cannot cross a solid cell face"));
                }
            }
            let scale: f64 = (0..6).map(|f| flow.outward(domain, cell, f).abs()).sum();
            let net = flow.net_outflow(domain, cell);
            if !scale.is_finite() || !net.is_finite() || net.abs() > 128.0 * f64::EPSILON * scale {
                return Err(invalid(
                    "constant-density enthalpy transport requires locally volume-conservative flow",
                ));
            }
        }
        for (face, h) in setup.faces.iter().zip(self.boundary_j_kg) {
            match (face, h) {
                (ThermalFace::Inflow { temperature: t }, Some(h))
                | (
                    ThermalFace::Outflow {
                        backflow_temperature: t,
                    },
                    Some(h),
                ) if t.is_finite() && *t > 0.0 && h.is_finite() => {}
                (ThermalFace::Inflow { .. } | ThermalFace::Outflow { .. }, _) => {
                    return Err(invalid(
                        "every open face needs enthalpy at its positive inlet/backflow temperature",
                    ));
                }
                (_, None) => {}
                _ => return Err(invalid("closed faces cannot declare advected enthalpy")),
            }
        }
        super::super::poll(gate)
    }
}

#[cfg(test)]
mod tests;
