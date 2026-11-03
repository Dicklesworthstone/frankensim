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

/// One constant-density fluid's frozen caloric field for conservative upwind
/// or limited linear enthalpy transport. Interior faces use the SAME upstream
/// tangent and reconstruction on both sides, so
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
        if !matches!(config.scheme, ConvectionScheme::Upwind | ConvectionScheme::EnthalpyMinmod) {
            return Err(invalid(
                "sourced enthalpy transport requires upwind or limited enthalpy reconstruction",
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

    /// Half-cell changes in specific enthalpy for the positive x/y/z faces.
    /// Compute once per frozen caloric image, outside matrix assembly/Krylov
    /// work. Every face then shares its donor's same three cached offsets.
    pub(crate) fn limited_offsets(
        &self,
        domain: &VoxelDomain,
        gate: &CancelGate,
    ) -> Result<Vec<[f64; 3]>, ChtError> {
        let invalid = || ChtError::InvalidInput {
            field: "thermal.advective_enthalpy",
            reason: "limited enthalpy reconstruction is unrepresentable".into(),
        };
        if self.cells.len() != domain.cell_count() {
            return Err(invalid());
        }
        let mut offsets = vec![[0.0; 3]; domain.cell_count()];
        for (cell, changes) in offsets.iter_mut().enumerate() {
            if cell.is_multiple_of(256) {
                super::super::poll(gate)?;
            }
            let Some(center) = self.cells[cell] else { continue };
            for (axis, change) in changes.iter_mut().enumerate() {
                let (Some(left), Some(right)) = (
                    domain.neighbor(cell, 2 * axis),
                    domain.neighbor(cell, 2 * axis + 1),
                ) else { continue };
                let (Some(left_h), Some(right_h)) = (self.cells[left], self.cells[right])
                else { continue };
                let half = 0.5 * domain.widths(cell)[axis];
                let left_distance = half + 0.5 * domain.widths(left)[axis];
                let right_distance = half + 0.5 * domain.widths(right)[axis];
                if !left_distance.is_finite() || left_distance <= 0.0
                    || !right_distance.is_finite() || right_distance <= 0.0
                {
                    return Err(invalid());
                }
                let left_slope = (center.enthalpy_j_kg - left_h.enthalpy_j_kg) / left_distance;
                let right_slope = (right_h.enthalpy_j_kg - center.enthalpy_j_kg) / right_distance;
                if !left_slope.is_finite() || !right_slope.is_finite() {
                    return Err(invalid());
                }
                // Compare signs directly: their product can overflow. A local
                // extremum has zero slope; otherwise the smaller magnitude
                // bounds both face traces by the adjacent source enthalpies.
                let slope = if left_slope > 0.0 && right_slope > 0.0 {
                    left_slope.min(right_slope)
                } else if left_slope < 0.0 && right_slope < 0.0 {
                    left_slope.max(right_slope)
                } else {
                    0.0
                };
                *change = half * slope;
                if !change.is_finite() {
                    return Err(invalid());
                }
            }
        }
        super::super::poll(gate)?;
        Ok(offsets)
    }
}

#[cfg(test)]
mod tests;
