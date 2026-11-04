//! Thermal laws on the existing stationary, impermeable halfway slab walls.
//!
//! The scalar may be absolute temperature or a declared affine temperature
//! coordinate. Heat flux is an inward scalar increment per wall link and time
//! step. For `theta = (T - T0) / delta_T`, dimensional inputs scale as
//! `q = q'' dt / (rho Cp dx delta_T)` and `transfer = h dt / (rho Cp dx)`.
//! The original collision owns diffusivity; no relaxation factor is applied
//! again to the incoming-minus-outgoing boundary exchange.

/// One spatially uniform thermal law on a planar slab wall.
///
/// All temperatures use the same coordinates as the thermal populations.
/// Positive heat flux enters the fluid at either wall. Geometry, material
/// properties, source scaling and validity remain with the existing slab owner.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ThermalWall {
    /// Prescribed halfway wall temperature; original anti-bounce-back law.
    Temperature(f64),
    /// Exactly zero normal scalar exchange, including nonequilibrium states.
    Adiabatic,
    /// Prescribed inward scalar increment per link over one lattice step.
    HeatFlux(f64),
    /// A film to a prescribed ambient, `q = transfer * (ambient - wall)`.
    /// The actual halfway wall temperature is solved with the reflected pair.
    Convective {
        /// Finite nonnegative `h dt / (rho Cp dx)`; zero is insulation.
        transfer: f64,
        /// Ambient temperature in the populations' scalar coordinates.
        ambient: f64,
    },
}

/// Invalid wall data; admission precedes mutation of the original slab.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ThermalWallError {
    /// A coefficient, temperature or prescribed flux cannot be admitted.
    Invalid(&'static str),
}

impl core::fmt::Display for ThermalWallError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(reason) => write!(f, "thermal wall: {reason}"),
        }
    }
}
impl std::error::Error for ThermalWallError {}

/// Measured exchange and surface trace for the accepted reflected populations.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ThermalWallExchange {
    /// Sum of actual incoming-minus-outgoing populations; positive into fluid.
    /// Multiply by the original physical cell heat-capacity/temperature scale
    /// to obtain energy. This is not an additional source to deposit.
    pub heat_in: f64,
    /// Arithmetic mean of halfway wall temperatures over its normal links.
    /// This is the surface trace, not the neighboring fluid-cell temperature.
    /// A prescribed temperature retains its exact authored value; the other
    /// wall laws reconstruct their trace from the accepted normal pair.
    pub mean_temperature: f64,
    /// Smallest individual wall trace; the checked owner admits local extrema.
    pub min_temperature: f64,
    /// Largest individual wall trace; a mean alone can conceal unsupported heat.
    pub max_temperature: f64,
}

impl ThermalWall {
    /// Check the authored scalar data. Negative temperature coordinates and
    /// cooling fluxes are allowed; a physical material owner admits absolute
    /// temperature and the resolved wall trace in its own validity interval.
    pub fn validate(self) -> Result<(), ThermalWallError> {
        let invalid = ThermalWallError::Invalid;
        match self {
            Self::Temperature(t) if !t.is_finite() => {
                Err(invalid("wall temperature must be finite"))
            }
            Self::HeatFlux(q) if !q.is_finite() => {
                Err(invalid("inward wall heat flux must be finite"))
            }
            Self::Convective { transfer, ambient }
                if !transfer.is_finite() || transfer < 0.0 || !ambient.is_finite() =>
            {
                Err(invalid(
                    "film transfer must be finite nonnegative and ambient finite",
                ))
            }
            _ => Ok(()),
        }
    }

    /// Reflect one post-collision outgoing population without allocating.
    ///
    /// The original ABB pair gives `T_w = (g_in + g_out)/(2w)`, while its
    /// conservative link exchange is `q = g_in - g_out`. Eliminating `T_w`
    /// from `q = H (T_ambient - T_w)` gives
    /// `g_in = g_out + 2H/(2w + H) * (w T_ambient - g_out)`.
    /// Two scaled branches retain both the insulating and prescribed-T limits
    /// without forming an overflowing `H/(2w)` or `2w + H`. This is a local
    /// planar boundary closure, not a curved-wall or unsteady-order theorem.
    /// Inputs have been admitted by the calling owner; nonrepresentable output
    /// remains visible to its existing numerical-state/temperature checks.
    pub(crate) fn reflect(self, outgoing: f64, weight: f64) -> (f64, f64) {
        let incoming = match self {
            // Retain the original fixed-temperature arithmetic exactly.
            Self::Temperature(temperature) => {
                return (2.0f64.mul_add(weight * temperature, -outgoing), temperature);
            }
            Self::Adiabatic => outgoing,
            Self::HeatFlux(flux) => outgoing + flux,
            Self::Convective { transfer, .. } if transfer == 0.0 => outgoing,
            Self::Convective { transfer, ambient } => {
                let target = weight * ambient;
                let pair_weight = 2.0 * weight;
                if transfer <= pair_weight {
                    let ratio = transfer / pair_weight;
                    let fraction = (2.0 * ratio) / (1.0 + ratio);
                    let difference = target - outgoing;
                    if difference.is_finite() {
                        fraction.mul_add(difference, outgoing)
                    } else {
                        // A convex reflection can remain finite even when
                        // subtracting opposite extreme populations overflows.
                        (1.0 - fraction).mul_add(outgoing, fraction * target)
                    }
                } else {
                    let ratio = pair_weight / transfer;
                    let fraction = (2.0 * ratio) / (1.0 + ratio);
                    let prescribed = 2.0f64.mul_add(target, -outgoing);
                    let difference = outgoing - target;
                    if prescribed.is_finite() && difference.is_finite() {
                        fraction.mul_add(difference, prescribed)
                    } else {
                        ((ratio - 1.0) / (1.0 + ratio))
                            .mul_add(outgoing, (2.0 / (1.0 + ratio)) * target)
                    }
                }
            }
        };
        (incoming, (0.5 * incoming + 0.5 * outgoing) / weight)
    }
}
