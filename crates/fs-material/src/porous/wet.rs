//! Prescribed liquid saturation and pore geometry for gas transport.
//!
//! This owner supplies instantaneous constitutive coefficients. The caller
//! supplies the current porosity, gas-conducting pore radius and saturation;
//! liquid mass, capillary pressure, swelling, damage and hysteresis are not
//! evolved here. Repeated sampling does not advance or mutate physical state.
//!
//! The gas relative permeability is the van Genuchten--Mualem nonwetting
//! relation, with explicitly calibrated residual liquid/gas saturations:
//! Se=(Sl-Slr)/(1-Slr-Sgr), krg=sqrt(1-Se)*(1-Se^(1/m))^(2m).
//! The declared endpoint extensions are krg=1 at Sl<=Slr and krg=0 at
//! Sl>=1-Sgr, within the caller's supported saturation interval.
//! This is a single constitutive branch, not a hysteresis model.
//! Formulation: Parker, Lenhard and Kuppusamy (1987),
//! <https://doi.org/10.1029/WR023i004p00618>; explicit residual-gas
//! normalization and exponents: Ruprecht et al. (2014), Eqs. (4),(6), Table 1,
//! <https://doi.org/10.1016/j.ijggc.2014.05.003>.
//!
//! Diffusion is selected independently: gas volume fraction theta=phi*(1-Sl)
//! gives either theta/tau with explicit transport tortuosity, or the
//! Millington--Quirk 1961 factor theta^(10/3)/phi^2.
//! <https://doi.org/10.1039/TF9615701200>; gas-diffusion form independently
//! stated by Lakshani et al. (2023), Eq. (6),
//! <https://doi.org/10.1002/vzj2.20236>.
//! Neither closure imposes a finite gas-connectivity threshold. Residual gas
//! immobility in the Darcy relation therefore does not automatically turn off
//! diffusion. Applicability to disconnected/trapped gas must be established
//! separately, or the supported saturation interval must exclude that region.
//!
//! These correlations supply no universal material calibration. Constructor
//! coefficients, applicability bounds and provenance are mandatory model
//! inputs. They do not acquire measured authority from a literature citation.
//! Temperature/pressure/composition applicability remains the consuming gas
//! law's explicit responsibility. No intrinsic permeability is inferred from
//! porosity or radius; the original skeleton permeability remains independent.

pub mod liquid;

use crate::MaterialError;
use fs_math::det;
use std::sync::Arc;

fn invalid(what: &str) -> MaterialError {
    MaterialError::Parameters { what: what.into() }
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn fraction(x: f64) -> bool {
    x.is_finite() && (0.0..=1.0).contains(&x)
}
fn contains(bounds: (f64, f64), x: f64) -> bool {
    x.is_finite() && x >= bounds.0 && x <= bounds.1
}

/// Explicit gas-diffusion model choice. There is intentionally no default.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum GasDiffusivityLaw {
    /// D_eff/D_bulk=theta_g^(10/3)/phi^2, with the current total porosity phi.
    /// The implied transport tortuosity is phi^(-1/3)*(1-Sl)^(-7/3).
    MillingtonQuirk1961,
    /// D_eff/D_bulk=theta_g/factor. The transport factor is >=1; it is not
    /// a path-length ratio that the implementation silently squares.
    ConstantTransportTortuosity { factor: f64 },
}

/// Inclusive geometric/saturation applicability of one fixed calibration.
/// Bounds are checked before applying either endpoint extension.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WetPoreSupport {
    pub liquid_saturation: (f64, f64),
    pub total_porosity: (f64, f64),
    pub mean_gas_pore_radius_m: (f64, f64),
}
impl WetPoreSupport {
    fn validate(self) -> Result<(), MaterialError> {
        if !fraction(self.liquid_saturation.0)
            || !fraction(self.liquid_saturation.1)
            || self.liquid_saturation.0 > self.liquid_saturation.1
            || !positive(self.total_porosity.0)
            || !fraction(self.total_porosity.1)
            || self.total_porosity.0 > self.total_porosity.1
            || !positive(self.mean_gas_pore_radius_m.0)
            || !positive(self.mean_gas_pore_radius_m.1)
            || self.mean_gas_pore_radius_m.0 > self.mean_gas_pore_radius_m.1
        {
            return Err(invalid(
                "wet pores require ordered saturation bounds in [0,1], porosity in (0,1], and positive gas-pore radius bounds",
            ));
        }
        Ok(())
    }
}

#[derive(Debug, PartialEq)]
struct WetPoreCalibration {
    residual_liquid_saturation: f64,
    residual_gas_saturation: f64,
    van_genuchten_m: f64,
    maximum_mobile_liquid_saturation: f64,
    mobile_saturation_span: f64,
    diffusivity: GasDiffusivityLaw,
    support: WetPoreSupport,
    provenance: String,
}

/// One explicit wet-pore calibration. Clones and successful samples share the
/// immutable calibration; sampling performs bounded scalar work and allocates
/// no memory. There are no hidden drainage/imbibition history variables.
#[derive(Clone, Debug, PartialEq)]
pub struct WetPoreGasLaw(Arc<WetPoreCalibration>);

impl WetPoreGasLaw {
    pub fn try_new(
        residual_liquid_saturation: f64,
        residual_gas_saturation: f64,
        van_genuchten_m: f64,
        diffusivity: GasDiffusivityLaw,
        support: WetPoreSupport,
        provenance: impl Into<String>,
    ) -> Result<Self, MaterialError> {
        support.validate()?;
        let maximum_mobile_liquid_saturation = 1.0 - residual_gas_saturation;
        let mobile_saturation_span = maximum_mobile_liquid_saturation - residual_liquid_saturation;
        let provenance = provenance.into();
        if !fraction(residual_liquid_saturation)
            || !fraction(residual_gas_saturation)
            || !positive(mobile_saturation_span)
            || !positive(van_genuchten_m)
            || van_genuchten_m >= 1.0
            || provenance.trim().is_empty()
            || provenance.len() > 4096
        {
            return Err(invalid(
                "wet pores require explicit residual saturations with sum <1, van Genuchten m in (0,1), and bounded provenance",
            ));
        }
        if let GasDiffusivityLaw::ConstantTransportTortuosity { factor } = diffusivity {
            if !factor.is_finite() || factor < 1.0 {
                return Err(invalid(
                    "wet-pore transport tortuosity must be finite and >=1",
                ));
            }
        }
        Ok(Self(Arc::new(WetPoreCalibration {
            residual_liquid_saturation,
            residual_gas_saturation,
            van_genuchten_m,
            maximum_mobile_liquid_saturation,
            mobile_saturation_span,
            diffusivity,
            support,
            provenance,
        })))
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
    pub fn diffusivity_law(&self) -> GasDiffusivityLaw {
        self.0.diffusivity
    }
    pub fn support(&self) -> WetPoreSupport {
        self.0.support
    }
    pub fn provenance(&self) -> &str {
        &self.0.provenance
    }

    /// Sample the prescribed current geometry. The radius is an independent
    /// model input: saturation does not imply a particular pore radius,
    /// Kozeny--Carman permeability, gas storage or displaced-liquid transfer.
    /// At Sl=1 the retained radius is inactive and both gas pathways close.
    pub fn sample(
        &self,
        liquid_saturation: f64,
        total_porosity: f64,
        mean_gas_pore_radius_m: f64,
    ) -> Result<WetPoreGasState, MaterialError> {
        let calibration = &self.0;
        let support = calibration.support;
        if !contains(support.liquid_saturation, liquid_saturation)
            || !contains(support.total_porosity, total_porosity)
            || !contains(support.mean_gas_pore_radius_m, mean_gas_pore_radius_m)
        {
            return Err(invalid(
                "wet-pore state is outside its saturation/porosity/radius calibration",
            ));
        }
        let (effective_liquid_saturation, gas_relative_permeability) =
            if liquid_saturation <= calibration.residual_liquid_saturation {
                (0.0, 1.0)
            } else if liquid_saturation >= calibration.maximum_mobile_liquid_saturation {
                (1.0, 0.0)
            } else {
                let liquid = (liquid_saturation - calibration.residual_liquid_saturation)
                    / calibration.mobile_saturation_span;
                // Use the actual remaining mobile-gas gap, not 1-Se, which
                // can round to zero for an interior state next to wet closure.
                let gas = (calibration.maximum_mobile_liquid_saturation - liquid_saturation)
                    / calibration.mobile_saturation_span;
                let log_liquid = if gas < 0.5 {
                    log_one_minus(gas)
                } else {
                    det::ln(liquid)
                };
                // expm1 preserves the capillary term near Se=1.
                let capillary = -det::expm1(log_liquid / calibration.van_genuchten_m);
                let permeability = det::sqrt(gas)
                    * det::exp((2.0 * calibration.van_genuchten_m) * det::ln(capillary));
                if !positive(permeability) || permeability > 1.0 {
                    return Err(invalid(
                        "nonzero wet-pore gas permeability is unrepresentable",
                    ));
                }
                (liquid, permeability)
            };
        let gas_saturation = 1.0 - liquid_saturation;
        let gas_filled_porosity = total_porosity * gas_saturation;
        let (effective_diffusivity_factor, transport_tortuosity_factor) = if gas_saturation == 0.0 {
            (0.0, None)
        } else {
            if !positive(gas_filled_porosity) {
                return Err(invalid("nonzero gas-filled pore volume is unrepresentable"));
            }
            let tortuosity = match calibration.diffusivity {
                GasDiffusivityLaw::ConstantTransportTortuosity { factor } => factor,
                GasDiffusivityLaw::MillingtonQuirk1961 =>
                // Equivalent to theta^(10/3)/phi^2 without forming
                // small powers whose ratio could underflow prematurely.
                {
                    det::exp(
                        (-1.0 / 3.0) * det::ln(total_porosity)
                            + (-7.0 / 3.0) * det::ln(gas_saturation),
                    )
                }
            };
            let factor = gas_filled_porosity / tortuosity;
            if !positive(tortuosity)
                || tortuosity < 1.0
                || !positive(factor)
                || factor > gas_filled_porosity
            {
                return Err(invalid(
                    "nonzero wet-pore gas diffusivity is unrepresentable",
                ));
            }
            (factor, Some(tortuosity))
        };
        Ok(WetPoreGasState {
            law: self.clone(),
            liquid_saturation,
            total_porosity,
            mean_gas_pore_radius_m,
            effective_liquid_saturation,
            gas_relative_permeability,
            gas_filled_porosity,
            effective_diffusivity_factor,
            transport_tortuosity_factor,
        })
    }
}

// Cancellation-correct ln(1-x), only called for x in (0,0.5). The small
// correction accounts for rounding 1-x; no platform libm or iteration.
fn log_one_minus(x: f64) -> f64 {
    let rounded = 1.0 - x;
    if rounded == 1.0 {
        -x
    } else {
        det::ln(rounded) * (x / (1.0 - rounded))
    }
}

/// Immutable instantaneous response retaining its exact calibration. Its
/// factors are not new intrinsic-material measurements. The Darcy consumer
/// applies krg to the original tensor, preserving axes and material sources.
/// The diffusion consumer uses the absolute factor once, not krg*factor.
#[derive(Clone, Debug, PartialEq)]
pub struct WetPoreGasState {
    law: WetPoreGasLaw,
    liquid_saturation: f64,
    total_porosity: f64,
    mean_gas_pore_radius_m: f64,
    effective_liquid_saturation: f64,
    gas_relative_permeability: f64,
    gas_filled_porosity: f64,
    effective_diffusivity_factor: f64,
    transport_tortuosity_factor: Option<f64>,
}
impl WetPoreGasState {
    pub const fn law(&self) -> &WetPoreGasLaw {
        &self.law
    }
    pub const fn liquid_saturation(&self) -> f64 {
        self.liquid_saturation
    }
    pub const fn total_porosity(&self) -> f64 {
        self.total_porosity
    }
    pub const fn mean_gas_pore_radius_m(&self) -> f64 {
        self.mean_gas_pore_radius_m
    }
    pub const fn effective_liquid_saturation(&self) -> f64 {
        self.effective_liquid_saturation
    }
    pub const fn gas_relative_permeability(&self) -> f64 {
        self.gas_relative_permeability
    }
    pub const fn gas_filled_porosity(&self) -> f64 {
        self.gas_filled_porosity
    }
    pub const fn effective_diffusivity_factor(&self) -> f64 {
        self.effective_diffusivity_factor
    }
    /// None only at full liquid saturation; no finite artificial transport
    /// factor or positive leakage floor is manufactured for the closed pores.
    pub const fn transport_tortuosity_factor(&self) -> Option<f64> {
        self.transport_tortuosity_factor
    }
}

#[cfg(test)]
mod tests;
