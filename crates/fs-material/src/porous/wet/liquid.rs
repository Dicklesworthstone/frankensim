//! Prescribed wetting-liquid saturation and capillary retention.
//!
//! An explicitly calibrated van Genuchten--Mualem branch supplies
//! Se=(Sl-Slr)/(1-Slr-Sgr), kr_l=Se^ell*[1-(1-Se^(1/m))^m]^2 and
//! Pc=P0*(Se^(-1/m)-1)^(1-m), with Pc=pg-pl in Pa. The connectivity
//! exponent ell is explicit and nonnegative; ell=0.5 is the usual Mualem
//! instance. The matching retention exponent is n=1/(1-m).
//! <https://documentation.pflotran.org/theory_guide/constitutive_relations.html>
//!
//! This is one calibrated instantaneous branch, not a drainage/imbibition
//! history or a Richards solver. Saturation is prescribed. No liquid inventory,
//! gas pressure, capillary work, porosity or wet Forchheimer coefficient is
//! inferred. A hydraulic consumer uses kr_l with the original intrinsic tensor;
//! it must not add Pc to its pressure loss as another drag term.

use super::{contains, fraction, invalid, log_one_minus, positive};
use crate::{MaterialError, porous::permeability::PorousGasSupport};
use fs_math::det;
use std::sync::Arc;

mod source;
pub use source::WetPoreLiquidProperties;
/// Conservative one-dimensional liquid saturation/storage evolution.
pub mod transient;

/// Inclusive applicability of this liquid calibration. Pressure means the
/// actual absolute LIQUID pressure, independently of the capillary difference.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WetPoreLiquidSupport {
    pub liquid_saturation: (f64, f64),
    pub temperature_k: (f64, f64),
    pub liquid_pressure_pa: (f64, f64),
    pub maximum_superficial_speed_m_s: f64,
}
impl WetPoreLiquidSupport {
    fn validate(self) -> Result<(), MaterialError> {
        if !fraction(self.liquid_saturation.0)
            || !fraction(self.liquid_saturation.1)
            || self.liquid_saturation.0 > self.liquid_saturation.1
            || !positive(self.temperature_k.0)
            || !positive(self.temperature_k.1)
            || self.temperature_k.0 > self.temperature_k.1
            || !positive(self.liquid_pressure_pa.0)
            || !positive(self.liquid_pressure_pa.1)
            || self.liquid_pressure_pa.0 > self.liquid_pressure_pa.1
            || !positive(self.maximum_superficial_speed_m_s)
        {
            return Err(invalid(
                "liquid retention requires ordered saturation, positive temperature/absolute-liquid-pressure and speed support",
            ));
        }
        Ok(())
    }
    pub(crate) const fn flow_support(self) -> PorousGasSupport {
        PorousGasSupport {
            temperature_k: self.temperature_k,
            pressure_pa: self.liquid_pressure_pa,
            maximum_superficial_speed_m_s: self.maximum_superficial_speed_m_s,
        }
    }
}

/// A dry endpoint has no finite equilibrium capillary pressure under this law.
/// No numerical cutoff manufactures a finite pressure or a positive mobility.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum LiquidCapillaryPressure {
    Finite(f64),
    DryLimit,
}

#[derive(Debug, PartialEq)]
struct LiquidCalibration {
    residual_liquid_saturation: f64,
    residual_gas_saturation: f64,
    maximum_liquid_saturation: f64,
    mobile_saturation_span: f64,
    van_genuchten_m: f64,
    connectivity_exponent: f64,
    pressure_scale_pa: f64,
    support: WetPoreLiquidSupport,
    provenance: String,
}

/// Immutable authored or source-bound calibration, distinct from the gas
/// wet-pore law. A literature citation alone supplies no measured authority.
#[derive(Clone, Debug, PartialEq)]
pub struct WetPoreLiquidLaw(Arc<LiquidCalibration>, Option<Arc<source::LiquidSource>>);
impl WetPoreLiquidLaw {
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        residual_liquid_saturation: f64,
        residual_gas_saturation: f64,
        van_genuchten_m: f64,
        connectivity_exponent: f64,
        pressure_scale_pa: f64,
        support: WetPoreLiquidSupport,
        provenance: impl Into<String>,
    ) -> Result<Self, MaterialError> {
        support.validate()?;
        let maximum_liquid_saturation = 1.0 - residual_gas_saturation;
        let mobile_saturation_span = maximum_liquid_saturation - residual_liquid_saturation;
        let provenance = provenance.into();
        if !fraction(residual_liquid_saturation)
            || !fraction(residual_gas_saturation)
            || !positive(mobile_saturation_span)
            || !positive(van_genuchten_m)
            || van_genuchten_m >= 1.0
            || !connectivity_exponent.is_finite()
            || connectivity_exponent < 0.0
            || !positive(pressure_scale_pa)
            || provenance.trim().is_empty()
            || provenance.len() > 4096
        {
            return Err(invalid(
                "liquid retention needs residual saturations with sum <1, m in (0,1), explicit nonnegative connectivity exponent, positive pressure scale and bounded provenance",
            ));
        }
        Ok(Self(
            Arc::new(LiquidCalibration {
                residual_liquid_saturation,
                residual_gas_saturation,
                maximum_liquid_saturation,
                mobile_saturation_span,
                van_genuchten_m,
                connectivity_exponent,
                pressure_scale_pa,
                support,
                provenance,
            }),
            None,
        ))
    }
    pub fn residual_liquid_saturation(&self) -> f64 {
        self.0.residual_liquid_saturation
    }
    pub fn residual_gas_saturation(&self) -> f64 {
        self.0.residual_gas_saturation
    }
    pub fn van_genuchten_m(&self) -> f64 {
        self.0.van_genuchten_m
    }
    pub fn connectivity_exponent(&self) -> f64 {
        self.0.connectivity_exponent
    }
    pub fn pressure_scale_pa(&self) -> f64 {
        self.0.pressure_scale_pa
    }
    pub fn support(&self) -> WetPoreLiquidSupport {
        self.0.support
    }
    pub fn provenance(&self) -> &str {
        &self.0.provenance
    }

    /// Sample a prescribed saturation without advancing physical state.
    /// Domain checks precede the explicit endpoint extensions: Sl<=Slr is
    /// immobile with no finite Pc; Sl>=1-Sgr has kr_l=1 and Pc=0.
    pub fn sample(&self, liquid_saturation: f64) -> Result<WetPoreLiquidState, MaterialError> {
        let c = &self.0;
        if !contains(c.support.liquid_saturation, liquid_saturation) {
            return Err(invalid(
                "liquid saturation is outside its retention calibration",
            ));
        }
        let (effective, relative, capillary) = if liquid_saturation <= c.residual_liquid_saturation
        {
            (0.0, 0.0, LiquidCapillaryPressure::DryLimit)
        } else if liquid_saturation >= c.maximum_liquid_saturation {
            (1.0, 1.0, LiquidCapillaryPressure::Finite(0.0))
        } else {
            let effective =
                (liquid_saturation - c.residual_liquid_saturation) / c.mobile_saturation_span;
            let remaining =
                (c.maximum_liquid_saturation - liquid_saturation) / c.mobile_saturation_span;
            let log_effective = if remaining < 0.5 {
                log_one_minus(remaining)
            } else {
                det::ln(effective)
            };
            let log_root = log_effective / c.van_genuchten_m;
            let root = det::exp(log_root);
            // Both dry cancellation 1-(1-z)^m and the wet gap 1-z retain
            // their represented positive increments rather than rounding away.
            let log_complement = if root < 0.5 {
                log_one_minus(root)
            } else {
                det::ln(-det::expm1(log_root))
            };
            let integral = -det::expm1(c.van_genuchten_m * log_complement);
            let relative =
                det::exp(c.connectivity_exponent * log_effective + 2.0 * det::ln(integral));
            let argument = -log_root;
            let log_retention = if argument < 0.5 {
                det::ln(det::expm1(argument))
            } else {
                argument + log_one_minus(det::exp(-argument))
            };
            let pressure =
                det::exp(det::ln(c.pressure_scale_pa) + (1.0 - c.van_genuchten_m) * log_retention);
            if !positive(effective) || !positive(relative) || relative > 1.0 || !positive(pressure)
            {
                return Err(invalid(
                    "positive liquid mobility or capillary pressure is unrepresentable",
                ));
            }
            (
                effective,
                relative,
                LiquidCapillaryPressure::Finite(pressure),
            )
        };
        Ok(WetPoreLiquidState {
            law: self.clone(),
            liquid_saturation,
            effective_liquid_saturation: effective,
            liquid_relative_permeability: relative,
            capillary_pressure: capillary,
        })
    }

    /// Invert the same retention branch for finite Pc=pg-pl>=0, in Pa.
    /// Negative pressures and unrepresentable positive-distance endpoints
    /// refuse instead of clipping into a saturated or numerically dry state.
    pub fn at_capillary_pressure(
        &self,
        capillary_pressure_pa: f64,
    ) -> Result<WetPoreLiquidState, MaterialError> {
        if !capillary_pressure_pa.is_finite() || capillary_pressure_pa < 0.0 {
            return Err(invalid(
                "retention inversion requires finite nonnegative Pc=pg-pl in Pa",
            ));
        }
        let c = &self.0;
        if capillary_pressure_pa == 0.0 {
            return self.sample(c.maximum_liquid_saturation);
        }
        let exponent = (det::ln(capillary_pressure_pa) - det::ln(c.pressure_scale_pa))
            / (1.0 - c.van_genuchten_m);
        // log(1+exp(exponent)) with no overflow or loss of small arguments.
        let tail = det::exp(-exponent.abs());
        let logarithm = exponent.max(0.0) - log_one_minus(tail / (1.0 + tail));
        let log_effective = -c.van_genuchten_m * logarithm;
        let effective = det::exp(log_effective);
        let saturation = if effective < 0.5 {
            c.residual_liquid_saturation + c.mobile_saturation_span * effective
        } else {
            c.maximum_liquid_saturation - c.mobile_saturation_span * (-det::expm1(log_effective))
        };
        if saturation <= c.residual_liquid_saturation || saturation >= c.maximum_liquid_saturation {
            return Err(invalid(
                "finite capillary pressure has an unrepresentable interior saturation",
            ));
        }
        self.sample(saturation)
    }
}

/// One frozen liquid constitutive state. Cloning retains the whole immutable
/// calibration; sampling performs bounded scalar work and no heap allocation.
#[derive(Clone, Debug, PartialEq)]
pub struct WetPoreLiquidState {
    law: WetPoreLiquidLaw,
    liquid_saturation: f64,
    effective_liquid_saturation: f64,
    liquid_relative_permeability: f64,
    capillary_pressure: LiquidCapillaryPressure,
}
impl WetPoreLiquidState {
    pub const fn law(&self) -> &WetPoreLiquidLaw {
        &self.law
    }
    pub const fn liquid_saturation(&self) -> f64 {
        self.liquid_saturation
    }
    pub const fn effective_liquid_saturation(&self) -> f64 {
        self.effective_liquid_saturation
    }
    pub const fn liquid_relative_permeability(&self) -> f64 {
        self.liquid_relative_permeability
    }
    pub const fn capillary_pressure(&self) -> LiquidCapillaryPressure {
        self.capillary_pressure
    }
    pub fn support(&self) -> WetPoreLiquidSupport {
        self.law.support()
    }
}

#[cfg(test)]
mod tests;
