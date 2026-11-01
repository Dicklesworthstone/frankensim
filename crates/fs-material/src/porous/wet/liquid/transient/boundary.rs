//! Hydraulic boundary conditions on the original porous-liquid column.

use super::{
    LiquidSaturationCell, LiquidSaturationStepInput, WetPoreLiquidState, admit_state, face_flux,
    hydraulic_mobility, invalid, positive,
};
use crate::MaterialError;

/// One boundary of the same one-dimensional porous column.
///
/// All superficial fluxes [m/s] are positive left-to-right, including at the
/// right boundary: positive left flux enters; positive right flux leaves.
/// Resistance relates a pressure jump to superficial flux, not total flow.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum LiquidSaturationBoundary {
    /// Exactly zero external flux; does not change the material permeability.
    Sealed,
    /// An imposed signed superficial liquid flux. A nonzero value cannot cross
    /// an exact intrinsic or dry seal. Its required material-side face pressure
    /// and actual speed must fit the original liquid source support.
    Flux { superficial_flux_m_s: f64 },
    /// Absolute external liquid pressure with a nonnegative surface hydraulic
    /// resistance [Pa s/m]. Zero resistance is the prescribed-pressure limit.
    /// The resistance has no thickness, liquid storage or body-force increment;
    /// the original porous half-cell resistance and body force remain separate.
    Pressure {
        liquid_pressure_pa: f64,
        surface_resistance_pa_s_m: f64,
    },
}
impl LiquidSaturationBoundary {
    fn validate(self) -> Result<(), MaterialError> {
        match self {
            Self::Sealed => Ok(()),
            Self::Flux {
                superficial_flux_m_s,
            } if superficial_flux_m_s.is_finite() => Ok(()),
            Self::Pressure {
                liquid_pressure_pa,
                surface_resistance_pa_s_m,
            } if positive(liquid_pressure_pa)
                && surface_resistance_pa_s_m.is_finite()
                && surface_resistance_pa_s_m >= 0.0 =>
            {
                Ok(())
            }
            _ => Err(invalid(
                "liquid boundary needs finite signed superficial flux or positive absolute pressure and nonnegative finite surface resistance",
            )),
        }
    }

    pub(super) fn flow(
        self,
        mobility: f64,
        pressure: Option<f64>,
        half_width: f64,
        body_gradient: f64,
        left: bool,
    ) -> Result<BoundaryFlow, MaterialError> {
        match self {
            Self::Sealed => Ok(BoundaryFlow::ZERO),
            Self::Flux {
                superficial_flux_m_s,
            } => {
                if mobility == 0.0 && superficial_flux_m_s != 0.0 {
                    return Err(invalid(
                        "nonzero prescribed liquid flux cannot cross a dry or intrinsic seal",
                    ));
                }
                Ok(BoundaryFlow {
                    flux: superficial_flux_m_s,
                    ..BoundaryFlow::ZERO
                })
            }
            Self::Pressure {
                liquid_pressure_pa,
                surface_resistance_pa_s_m,
            } => {
                if mobility == 0.0 {
                    return Ok(BoundaryFlow::ZERO);
                }
                if !positive(half_width) {
                    return Err(invalid("liquid boundary half-width is unrepresentable"));
                }
                let pressure = pressure
                    .ok_or_else(|| invalid("mobile liquid boundary pressure is missing"))?;
                let bare_conductance = mobility / half_width;
                let (conductance, fraction, flux) = if surface_resistance_pa_s_m == 0.0 {
                    // Retain the original prescribed-pressure arithmetic.
                    let (p_left, p_right) = if left {
                        (liquid_pressure_pa, pressure)
                    } else {
                        (pressure, liquid_pressure_pa)
                    };
                    (
                        bare_conductance,
                        1.0,
                        face_flux(mobility, p_left, p_right, half_width, body_gradient)?,
                    )
                } else {
                    // C = 1/(d/M + R); f = (d/M)/(d/M + R). Avoid overflowing
                    // the series sum or the product R*C0 for a dominant surface.
                    let product = surface_resistance_pa_s_m * bare_conductance;
                    let (conductance, fraction) = if bare_conductance.is_finite() && product <= 1.0
                    {
                        let fraction = 1.0 / (1.0 + product);
                        (bare_conductance * fraction, fraction)
                    } else {
                        let inverse = (half_width / mobility) / surface_resistance_pa_s_m;
                        let denominator = 1.0 + inverse;
                        (
                            (1.0 / surface_resistance_pa_s_m) / denominator,
                            inverse / denominator,
                        )
                    };
                    if !positive(conductance)
                        || !fraction.is_finite()
                        || !(0.0..=1.0).contains(&fraction)
                    {
                        return Err(invalid(
                            "liquid boundary series conductance is unrepresentable",
                        ));
                    }
                    let difference = if left {
                        liquid_pressure_pa - pressure
                    } else {
                        pressure - liquid_pressure_pa
                    };
                    let flux = conductance * (difference + body_gradient * half_width);
                    (conductance, fraction, flux)
                };
                if !flux.is_finite() {
                    return Err(invalid("liquid boundary flux is unrepresentable"));
                }
                Ok(BoundaryFlow {
                    flux,
                    pressure_conductance: conductance,
                    mobility_fraction: fraction,
                })
            }
        }
    }

    /// Actual porous-side pressure for this already evaluated flux. External
    /// pressure across a contact is not the porous material's source coordinate.
    fn material_pressure(
        self,
        mobility: f64,
        pressure: Option<f64>,
        half_width: f64,
        body_gradient: f64,
        left: bool,
        flux: f64,
    ) -> Result<Option<f64>, MaterialError> {
        if mobility == 0.0 {
            if flux != 0.0 {
                return Err(invalid(
                    "nonzero liquid boundary flux crosses an exact seal",
                ));
            }
            return Ok(None);
        }
        if !positive(half_width) {
            return Err(invalid("liquid boundary half-width is unrepresentable"));
        }
        let sign = if left { 1.0 } else { -1.0 };
        let material = match self {
            Self::Pressure {
                liquid_pressure_pa,
                surface_resistance_pa_s_m,
            } => liquid_pressure_pa - sign * flux * surface_resistance_pa_s_m,
            Self::Sealed | Self::Flux { .. } => {
                let pressure = pressure
                    .ok_or_else(|| invalid("mobile liquid boundary pressure is missing"))?;
                pressure + sign * half_width * (flux / mobility - body_gradient)
            }
        };
        if !positive(material) {
            return Err(invalid(
                "required material-side liquid boundary pressure is unrepresentable or nonpositive",
            ));
        }
        Ok(Some(material))
    }
}

/// Prescribed thermohydraulic coordinates with independently typed left/right
/// boundary conditions, held over one physical interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LiquidSaturationBoundaryStepInput {
    pub duration_s: f64,
    pub temperature_k: f64,
    pub gas_pressure_pa: f64,
    pub dynamic_viscosity_pa_s: f64,
    /// Signed body-force contribution rho*g·x_hat [Pa/m], inside porous cells.
    pub body_force_pressure_gradient_pa_m: f64,
    pub boundaries: [LiquidSaturationBoundary; 2],
}
impl From<LiquidSaturationStepInput> for LiquidSaturationBoundaryStepInput {
    fn from(input: LiquidSaturationStepInput) -> Self {
        Self {
            duration_s: input.duration_s,
            temperature_k: input.temperature_k,
            gas_pressure_pa: input.gas_pressure_pa,
            dynamic_viscosity_pa_s: input.dynamic_viscosity_pa_s,
            body_force_pressure_gradient_pa_m: input.body_force_pressure_gradient_pa_m,
            boundaries: input.boundary_liquid_pressure_pa.map(|liquid_pressure_pa| {
                LiquidSaturationBoundary::Pressure {
                    liquid_pressure_pa,
                    surface_resistance_pa_s_m: 0.0,
                }
            }),
        }
    }
}
impl LiquidSaturationBoundaryStepInput {
    pub(super) fn validate(self) -> Result<(), MaterialError> {
        if !positive(self.duration_s)
            || !positive(self.temperature_k)
            || !positive(self.gas_pressure_pa)
            || !positive(self.dynamic_viscosity_pa_s)
            || !self.body_force_pressure_gradient_pa_m.is_finite()
        {
            return Err(invalid(
                "transient liquid step requires positive dt/T/gas pressure/viscosity and finite body-force gradient",
            ));
        }
        for boundary in self.boundaries {
            boundary.validate()?;
        }
        Ok(())
    }

    pub(super) fn admit_faces(
        self,
        cells: &[LiquidSaturationCell],
        states: &[WetPoreLiquidState],
        pressures: &[Option<f64>],
        faces: &[f64],
    ) -> Result<(), MaterialError> {
        let last = cells.len() - 1;
        for (side, cell, face) in [(0, 0, 0), (1, last, cells.len())] {
            let mobility =
                hydraulic_mobility(&cells[cell], &states[cell], self.dynamic_viscosity_pa_s)?;
            let material_pressure = self.boundaries[side].material_pressure(
                mobility,
                pressures[cell],
                0.5 * cells[cell].thickness_m,
                self.body_force_pressure_gradient_pa_m,
                side == 0,
                faces[face],
            )?;
            if let Some(pressure) = material_pressure {
                admit_state(
                    &states[cell],
                    self.temperature_k,
                    pressure,
                    faces[face].abs(),
                )?;
            }
        }
        Ok(())
    }
}

pub(super) struct BoundaryFlow {
    pub(super) flux: f64,
    pub(super) pressure_conductance: f64,
    /// Sensitivity of ln(conductance) to ln(cell mobility), in [0,1].
    pub(super) mobility_fraction: f64,
}
impl BoundaryFlow {
    const ZERO: Self = Self {
        flux: 0.0,
        pressure_conductance: 0.0,
        mobility_fraction: 0.0,
    };
}
